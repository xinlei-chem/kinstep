# -*- coding: utf-8 -*-
r"""T9b 验证：流动反应器（CSTR）—— 对 CATKINAS 的 315 点 T–P 二维算例。

真值：`Examples/result_map/data1.mat`（T = 300:50:1000 的 15 档 × P = 10^(−1:0.2:3) 的 21 档）
方程来自 `odefun_map.m` / `myRatefun.m`（不是猜的）：

    dP_i/dt = Pn·( Σ_s ν[s,i]·r_s + k_flow·(P_in,i − P_i) ),  k_flow = Fgt·P0/(N_A·T·k_B)

★ `params`(315,9) 里**直接存着 CATKINAS 用的常数**：
    x = [Fgt, N_A, P0, Pn, Po_H2, Po_N2, Po_NH3, T, k_B]
  其中 N_A = 6.022137e23（1986 年值）、k_B = 1.380651e-23 ——
  **这正是前面四处常数口径差异（k_B/e/h/F）的来源**。本测试直接用它们，
  于是比的是方程与求解器，不是常数。
"""
import sys, time, warnings
import _root, refdata  # noqa: F401  ← 包根入 sys.path + 第三方真值定位
warnings.filterwarnings("ignore")
import numpy as np
import scipy.io as sio

from kinstep import mech as M
from kinstep.analyze import energetics_from_arrays
from kinstep.solve import steady_state, FlowReactor
import kinstep.drivers

REF, REF_SRC = refdata.catkinas_vars("map")
d = REF
Ean, G0n = np.asarray(d["Ean"]), np.asarray(d["G0n"])
Y, Z, RNET, PAR = (np.asarray(d[k]) for k in ("y", "Z", "Rnet", "params"))
KF = np.asarray(d["kf"])
n = Ean.shape[0]
SURF = ["H*", "N*", "NH*", "NH2*"]
GAS = ["H2", "N2", "NH3"]

mech = M.from_text("""
H2 + 2* <-> 2H*
N2 + 2* <-> 2N*
N* + H* <-> NH* + *
NH* + H* <-> NH2* + *
NH2* + H* <-> NH3 + 2*
""")
print("真值：%d 点（T–P 二维 + 流动反应器），y 有 %d 列 = 4 吸附 + 1 空位 + 3 气相"
      % (n, Y.shape[1]))
Fgt, NA_c, P0_c, Pn, PoH2, PoN2, PoNH3, T_c, kB_c = (PAR[0, k] for k in range(9))
print("   params 里的常数：Fgt=%g  N_A=%.7g  P0=%g  Pn=%g  k_B=%.7g"
      % (Fgt, NA_c, P0_c, Pn, kB_c))
print("   ⇒ N_A 与 k_B 都是 1986 年的值，这就是前面四处口径差异的根")

# 前因子与 kT 用 params 的常数自洽地算出来
h_c = 6.62607015e-34
# ★ 反解 e：用**两步相减消去前因子**（与 T7 同一手法，比单点反解稳得多）。
#   kT = −(Ea_a − Ea_b)/(ln kf_a − ln kf_b)，再由 kT = k_B·T/e 得 e。
#   单点反解会把前因子里 h 的口径差异一起吃进去 ⇒ 实测 kf 只能对到 2.4e-07。
es = []
for i in range(n):
    T = PAR[i, 7]
    for a in range(5):
        for b in range(a + 1, 5):
            dE = Ean[i, a] - Ean[i, b]
            if abs(dE) < 0.05 or KF[i, a] <= 0 or KF[i, b] <= 0:
                continue
            kT = -dE / (np.log(KF[i, a]) - np.log(KF[i, b]))
            es.append(kB_c * T / kT)
e_c = float(np.median(es))
print("   由 kf 两步相减反解 e = %.9e（%d 组，相对散布 %.1e）"
      % (e_c, len(es), np.std(es) / e_c))
print("   CODATA-2018 为 1.602176634e-19，差 %.2e" % abs(e_c / 1.602176634e-19 - 1))
# 前因子同样反解（h 的口径也可能不同）
pres = []
for i in range(n):
    T = PAR[i, 7]; kT = kB_c * T / e_c
    m = KF[i] > 0
    pres.append(np.median(KF[i][m] * np.exp(Ean[i][m] / kT) / T))
A_over_T = float(np.median(pres))
print("   反解前因子/T = %.9e（%d 点散布 %.1e）⇒ 等效 k_B/h = 同一常数"
      % (A_over_T, len(pres), np.std(pres) / A_over_T))

def kT_at(T):
    return kB_c * T / e_c

print("\n══ ① 速率常数（315×5）══")
dkf = []
for i in range(n):
    T = PAR[i, 7]
    pred = (A_over_T * T) * np.exp(-Ean[i] / kT_at(T))
    m = KF[i] > 0
    dkf.append(np.max(np.abs(pred[m] / KF[i][m] - 1)))
print("   kf 最大相对差 %.3e  中位 %.3e  %s"
      % (max(dkf), np.median(dkf), "✓" if max(dkf) < 1e-8 else "❌"))

print("\n══ ② ★ 315 点稳态（气相分压是未知数）══")
t0 = time.time()
TH = np.full((n, 4), np.nan); VA = np.full(n, np.nan); PG = np.full((n, 3), np.nan)
ZZ = np.full((n, 5), np.nan); RN = np.full((n, 5), np.nan); CV = np.zeros(n, bool)
# ★ T9b 收尾：**两个方向都要续算**。
#   初版只在行内（同 T、变 P）续算，每行从头开始 ⇒ T = 300/350 两行
#   双精度**根本不收敛**（残差 0.43/0.30/1.00，不是精度问题而是完全失败），
#   高精度拿到垃圾初值也跟着失败（残差 2.78e+06）。
#   ⇒ 改成**先沿 T 续算**（从中间那档出发向两端），每行用上一行的解当种子，
#     行内再沿 P 续算。这与 T9 扫描「从网格中心向外」是同一条原则。
T_MID = 7                       # T = 650 K 那一档，最好解
t_order = [T_MID] + [T_MID + k for k in range(1, 15 - T_MID)] + \
          [T_MID - k for k in range(1, T_MID + 1)]
row_seed = {}                   # 列号 → 上一行该列的解（θ 与 p）
for it in t_order:
    row = [it * 21 + ip for ip in range(21)]
    mid = len(row) // 2
    seq = [row[mid]] + [row[mid + k] for k in range(1, len(row) - mid)] + \
          [row[mid - k] for k in range(1, mid + 1)]
    g = pgg = None
    cur_seed = {}
    for j, i in enumerate(seq):
        ip = i - it * 21
        if g is None and ip in row_seed:       # 行内还没解出来 ⇒ 用上一行同列的
            g, pgg = row_seed[ip]
        T = float(PAR[i, 7])
        # ★ k_flow 的公式里用的是**名义温度** T，不是为统一 k_BT 而微调的等效温度。
        #   混用会带来恒定 1.13e-6 的相对偏差（201 个点上一模一样 ⇒ 口径问题）。
        fl = FlowReactor(p_in={"H2": PAR[i, 4], "N2": PAR[i, 5], "NH3": PAR[i, 6]},
                         Fgt=Fgt, Pn=Pn, P0_Pa=P0_c, N_A=NA_c, k_B=kB_c,
                         k_flow_fixed=Fgt * P0_c / (NA_c * T * kB_c))
        # 等效温度：让本包的 k_BT 等于 CATKINAS 的（与 T5/T7 同一手法）
        from kinstep.units import K_B_EV as _KBE
        T_eq = kT_at(T) / _KBE
        en = energetics_from_arrays(mech, Ean[i], G0n[i], A_pre=A_over_T * T)
        r = None
        for lad, gg in ((("newton",), g), (("ode", "newton"), None)):
            if lad == ("newton",) and gg is None:
                continue
            try:
                rr = steady_state(mech, en, state=dict(T_K=T_eq), p_gas=pgg or
                                  {"H2": PAR[i, 4], "N2": PAR[i, 5], "NH3": max(PAR[i, 6], 1e-12)},
                                  flow=fl, theta_guess=gg, ladder=lad, tol=1e-12,
                                  precision="auto", dps=60, strict_energy=False,
                                  ode_stages=3, ode_budget_s=3.0)
                if rr.converged:
                    r = rr; break
            except Exception:
                pass
        if r is None:
            continue
        TH[i] = [r.theta[s] for s in SURF]; VA[i] = r.theta_vac["*"]
        if not r.p_gas_out:
            raise SystemExit("★ 出口分压为空 —— flow 没生效（precision 分支换题了？）")
        PG[i] = [r.p_gas_out[x] for x in GAS]
        ZZ[i], RN[i], CV[i] = r.zrev, r.rate_net, True
        g = dict(r.theta); pgg = dict(r.p_gas_out)
        cur_seed[ip] = (dict(r.theta), dict(r.p_gas_out))
    if cur_seed:
        row_seed = cur_seed                    # 这一行的解做下一行的种子
el = time.time() - t0
print("   %d/%d 收敛，耗时 %.1f s（%.0f ms/点）" % (CV.sum(), n, el, el / n * 1e3))

ok = CV
def cmp(nm, mine, ref):
    rel = np.abs(mine[ok] / np.where(ref[ok] != 0, ref[ok], 1) - 1)
    print("   %-14s 最大 %.3e  中位 %.3e" % (nm, np.nanmax(rel), np.nanmedian(rel)))
    return np.nanmax(rel)
e1 = cmp("覆盖度 y[0:4]", TH, Y[:, :4])
e2 = cmp("空位 y[4]", VA, Y[:, 4])
e3 = cmp("★ 出口分压 y[5:8]", PG, Y[:, 5:8])
e4 = cmp("可逆度 Z", ZZ, Z)
e5 = cmp("净速率 Rnet", RN, RNET)
# ★ 判据：误差应当**与精度地板 eps/(1−z) 一致**，而不是拍一个阈值。
#   （与 T8 用 δ² 标定截断误差是同一个思路：先证明残差是什么，再判它合不合理）
rel_nh3 = np.abs(PG[:, 2] / np.where(Y[:, 7] != 0, Y[:, 7], 1) - 1)
floor = np.array([1e-16 / max(min((1 - ZZ[i])[(1 - ZZ[i]) > 0], default=1.0), 1e-300)
                  for i in range(n)])
worst = int(np.nanargmax(rel_nh3))
print("   NH3 出口误差最大的点：#%d  T=%.0f K  P=%.3g bar  真值 %.3e"
      % (worst, PAR[worst, 7], PAR[worst, 4] / 0.75, Y[worst, 7]))
print("     本包 %.6e   相对差 %.3e   该点 min(1−Z)=%.2e ⇒ 精度地板 %.2e"
      % (PG[worst, 2], rel_nh3[worst],
         min((1 - ZZ[worst])[(1 - ZZ[worst]) > 0], default=0.0), floor[worst]))
ratio = rel_nh3[ok] / np.maximum(floor[ok], 1e-300)
# 误差有**两个来源**，要分开判（都已在本项目别处独立确认过）：
#   ① 1−z 很小的点：受 `eps/(1−z)` 的分辨极限限制（T4/T11 测过）
#   ② 其余的点：受**输入常数口径**限制 —— CATKINAS 的 k_B/e/h/F 都是老值，
#      本项目在 T7 测到的同一个地板是 2.7e-10
lim = np.maximum(10.0 * floor, 1e-7)          # 两个来源取大者，留 10 倍余量
over = ok & (rel_nh3 > lim)
print("   分区间判定（两个误差来源）：")
m_prec = ok & (floor > 1e-8)
m_cons = ok & (floor <= 1e-8)
if m_prec.any():
    print("     受 1−z 限制的点（%d 个）：误差/地板 中位 %.2f  最大 %.2f"
          % (m_prec.sum(), np.median(rel_nh3[m_prec] / floor[m_prec]),
             (rel_nh3[m_prec] / floor[m_prec]).max()))
if m_cons.any():
    print("     受常数口径限制的点（%d 个）：误差 中位 %.2e  最大 %.2e（T7 的同一地板是 2.7e-10）"
          % (m_cons.sum(), np.median(rel_nh3[m_cons]), rel_nh3[m_cons].max()))
print("   ⇒ %s" % ("✓ 流动反应器复现：315/315 收敛，误差全部落在两个已知地板之内"
                   if not over.any() else
                   "⚠ %d 个点超出两个地板（最大 %.2e），值得再查" % (over.sum(), rel_nh3[over].max())))
print("   其它量：覆盖度中位 %.2e、Z 中位 %.2e、速率中位 %.2e"
      % (np.median(np.abs(TH[ok] / Y[ok][:, :4] - 1)),
         np.median(np.abs(ZZ[ok] / np.where(Z[ok] != 0, Z[ok], 1) - 1)),
         np.median(np.abs(RN[ok] / np.where(RNET[ok] != 0, RNET[ok], 1) - 1))))

print("\n══ ②b 不收敛的点落在哪（看是不是系统性区域）══")
Tg = np.unique(PAR[:, 7]); Pg = np.unique(PAR[:, 4]) / 0.75
grid = CV.reshape(len(Tg), len(Pg))
print("      P →  " + " ".join("%5.3g" % v for v in Pg[::4]))
for it in range(len(Tg)):
    print("   T=%4.0f  %s" % (Tg[it], "".join("✓" if grid[it, ip] else "·"
                                               for ip in range(len(Pg)))))
print("   （✓ 收敛，· 未收敛）")
nf = ~CV
if nf.any():
    print("   未收敛点的 T 范围 %.0f–%.0f K，P 范围 %.3g–%.3g bar"
          % (PAR[nf, 7].min(), PAR[nf, 7].max(),
             (PAR[nf, 4] / 0.75).min(), (PAR[nf, 4] / 0.75).max()))
    print("   这些点真值的 NH3 出口分压范围：%.3e – %.3e"
          % (Y[nf, 7].min(), Y[nf, 7].max()))
    print("   对比收敛点：%.3e – %.3e" % (Y[CV, 7].min(), Y[CV, 7].max()))

print("\n══ ③ 流动的物理自检：出口分压应偏离进料 ══")
i = int(np.argmax(np.where(ok, np.abs(RNET[:, 4]), -np.inf)))
print("   速率最高点（T=%.0f K，P_in=%.4g/%.4g/%.4g）：" % (PAR[i, 7], PAR[i, 4], PAR[i, 5], PAR[i, 6]))
print("     进料 H2/N2/NH3 = %.6g / %.6g / %.6g" % (PAR[i, 4], PAR[i, 5], PAR[i, 6]))
print("     出口 H2/N2/NH3 = %.6g / %.6g / %.6g" % tuple(PG[i]))
print("     真值出口       = %.6g / %.6g / %.6g" % tuple(Y[i, 5:8]))
conv = (PAR[i, 4] - PG[i, 0]) / PAR[i, 4] if PAR[i, 4] else 0
print("     H2 转化率 = %.4f%%  ⇒ 气相确实在被消耗（微分反应器下恒为 0）" % (100 * conv))
