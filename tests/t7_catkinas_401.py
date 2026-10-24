# -*- coding: utf-8 -*-
r"""T7：★ CATKINAS 401 点逐数字回归（本项目的主回归）

真值
----
`工具/小程序1/Examples/result_curve/data1.mat` —— CATKINAS 自带算例的输出，
合成氨 5 步，描述符 E_N 从 −2 到 +2 步长 0.01 共 **401 个采样点**。
输入见同目录 `INPUT_curve.m`：T=673 K，P(H2,N2,NH3)=(75,25,1)。

★ 口径是**反推+核对**出来的，不是猜的（见 notes/DEVLOG.md T7 条）
-------------------------------------------------------------------
  kf   = (k_BT/h)·exp(−Ean/k_BT)        ← 对全部 5 步成立（含吸附步），pred/kf = 1.0000
  kf/kr= exp(−G0n/k_BT)                 ← ratio = 1.0
  压力 = 无量纲 bar 数（P0 = 1 bar）    ← 用准平衡的步 1 反推得 75.000000
  y    = [H*, N*, NH*, NH2*, *]（含空位，和为 1）
  Z    = 可逆度（本包的 zrev）

⇒ `Ean/G0n` 已把 ThermoMode=4 的熵校正与 BarrierMode=1 的碰撞理论**全吸收进去**，
  所以本回归**绕过**还没做的热力学库（T10），**干净地只测内核**：
  速率常数 → 稳态求解 → 可逆度 → 通量。T10 做完后再回归 Eao/G0o → Ean/G0n 那一层。

★ 这个算例正是通量投影的主场
----------------------------
步1 有 **267/401** 个点 1−Z < 1e-6（步3: 174，步4: 162，步5: 180）。
CATKINAS 为此默认开 **1000 位算术**；本包用 double + 通量投影，看能走到哪。
"""
import sys, time
import _root, refdata  # noqa: F401  ← 包根入 sys.path + 第三方真值定位
import numpy as np
import scipy.io as sio

from kinstep import mech as M
from kinstep.units import K_B_EV
from kinstep.rates import Energetics, compile_rates, apply_state
from kinstep.solve import steady_state
import kinstep.drivers

REF, REF_SRC = refdata.catkinas_vars("curve")

# ★ 等效温度：让两边的 k_BT 完全相同（与 T5 里缩放 U 同一手法）
#
# CATKINAS 硬编码的是 **CODATA-2006 的** k_B ≈ 8.617343e-5 eV/K
# （从 401×5 个 kf 用"两步相减消去前因子"反解得 kT = 0.057994719627 eV，
#   103 组数据的相对散布只有 6.4e-16 ⇒ 确实是单一常数）。
# 本包用 CODATA-2018，kT(673K) = 0.057994652854 eV，相对差 1.151e-06。
# 这点差异被 Ea/kT ≈ 25 放大成 kf 的 ~4e-5 —— **不处理就测的是常数差异，不是公式差异**。
#
# 办法：本包按 T' = T·(kT_cat/kT_mine) 计算，则 k_BT 逐位相同，
# 而 **前因子 k_BT/h 与指数 Ea/kT 同时对上**（实测前因子 1.402306e+13 逐位一致）。
T_NOMINAL = 673.0
KT_CATKINAS = 0.057994719627
from kinstep.units import K_B_EV as _KBE
T_K = T_NOMINAL * KT_CATKINAS / (_KBE * T_NOMINAL)
print("等效温度 T' = %.9f K（名义 %.1f K，偏 %.3e）—— 统一两边的 k_BT"
      % (T_K, T_NOMINAL, T_K / T_NOMINAL - 1))
P = {"H2": 75.0, "N2": 25.0, "NH3": 1.0}
SP = ["H*", "N*", "NH*", "NH2*"]          # y 的前 4 列；第 5 列是空位

# ★ 前因子也从真值反解（CATKINAS 的 k_B、h 都是老常数）：
#   pre = kf·exp(Ean/kT)，2005 个样本的相对散布 ~1e-16 ⇒ 确实是单一常数。
#   本包 k_BT/h = 1.402304467e13，CATKINAS = 1.402306113158e13，相对差 1.174e-06。
#   显式传 A_pre ⇒ 两边常数完全相同，**测的就是公式与求解器本身**。
A_PRE_CATKINAS = 1.402306113158e13

TXT = """
H2 + 2* <-> 2H*
N2 + 2* <-> 2N*
N* + H* <-> NH* + *
NH* + H* <-> NH2* + *
NH2* + H* <-> NH3 + 2*
"""

d = REF
print("   参考值来源 = %s" % REF_SRC)

# ★ 存档副本会不会跑偏？两边都在场时逐位比一遍 —— 否则"自足"只是看起来自足：
#   发行物升级、或重新生成 npz 时漏了某个数组，回归照样全绿，而比的已经不是同一组数。
if REF_SRC == "mat":
    import os as _os
    _npz = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)),
                         "catkinas_refdata.npz")
    if _os.path.exists(_npz):
        _z = np.load(_npz)
        _bad = []
        for _k in ("Ean", "G0n", "kf", "kr", "Z", "y", "Rnet"):
            _a, _b = np.asarray(d[_k], float), np.asarray(_z["curve__" + _k], float)
            if _a.shape != _b.shape or not np.array_equal(_a, _b):
                _bad.append(_k)
        print("   存档副本 vs 发行物原件：%s"
              % ("逐位相同 ✓" if not _bad else "❌ 不一致：%s" % ", ".join(_bad)))
        if _bad:
            # 这是完整性问题，不是精度问题：继续往下比没有意义，直接判失败。
            sys.stderr.write("catkinas_refdata.npz 与发行物原件不一致：%s\n"
                             "  重新生成：python tests/make_catkinas_refdata.py\n"
                             % ", ".join(_bad))
            sys.exit(1)

Ean, G0n = np.asarray(d["Ean"]), np.asarray(d["G0n"])
KF, KR, Z, Y, RNET = (np.asarray(d[k]) for k in ("kf", "kr", "Z", "y", "Rnet"))
E1 = np.arange(-2.0, 2.0001, 0.01)
n = Ean.shape[0]
assert n == 401 == len(E1), (n, len(E1))

mech = M.from_text(TXT)
assert mech.surface_species == SP, mech.surface_species
print("机理解析：%d 步，表面物种 %s，空位 %s"
      % (mech.n_step, mech.surface_species, [mech.vacancy_of(s) for s in mech.sites]))
print("真值：%d 个采样点（E_N = %.2f … %.2f eV），T = %.0f K，P = %s"
      % (n, E1[0], E1[-1], T_K, P))

# ══════════ ① 速率常数逐点对数（401×5 = 2005 个 kf + 2005 个 kr）══════════
print("\n══ ① 速率常数 kf / kr：401×5 逐点 ══")
dkf = np.zeros_like(KF); dkr = np.zeros_like(KR)
def mk_en(i):
    return Energetics(mech, {s + 1: dict(Ga_fwd=float(Ean[i, s]), dG_rxn=float(G0n[i, s]),
                                         A_pre=A_PRE_CATKINAS) for s in range(5)})

for i in range(n):
    rc = apply_state(compile_rates(mk_en(i)), T_K)
    dkf[i] = np.abs(rc.kf / KF[i] - 1.0)
    dkr[i] = np.abs(rc.kr / KR[i] - 1.0)
print("   kf 相对差：最大 %.3e   中位 %.3e" % (dkf.max(), np.median(dkf)))
print("   kr 相对差：最大 %.3e   中位 %.3e" % (dkr.max(), np.median(dkr)))
print("   ⇒ %s" % ("✅ 逐点一致" if max(dkf.max(), dkr.max()) < 1e-9
                   else "❌ 有偏差，先查口径"))
print("   判据为何是 1e-9 而不是机器精度：kT 与前因子是**从真值反解**的"
      "（kT 用两步相减消前因子，103 组散布 6.4e-16；前因子 2005 样本），"
      "\n     反解本身就有 ~1e-11 的精度，再叠加 .mat 以 float64 存 1000 位算术结果的截断 "
      "⇒ ~1e-10 就是这个比法的地板。")

# ══════════ ② 稳态：覆盖度 / 可逆度 / 净速率 ══════════
print("\n══ ② 稳态求解：覆盖度 y、可逆度 Z、净速率 Rnet ══")
t0 = time.time()
TH = np.zeros((n, 4)); VAC = np.zeros(n); ZZ = np.zeros((n, 5)); RN = np.zeros((n, 5))
FIX = np.zeros(n); RES = np.zeros(n); NIT = np.zeros(n, int)
fail = []

def solve_at(i, guess):
    """先续算（用邻点的解当初值），不成再退回完整阶梯。"""
    for attempt in (("newton",), ("ode", "newton")):
        if attempt == ("newton",) and guess is None:
            continue
        try:
            r = steady_state(mech, mk_en(i), state=dict(T_K=T_K), p_gas=P,
                             theta_guess=guess if attempt == ("newton",) else None,
                             ladder=attempt, tol=1e-12,
                             ode_stages=4, ode_budget_s=6.0)
            if r.converged:
                return r
        except Exception:
            pass
    return None

def store(i, r):
    TH[i] = [r.theta[s] for s in SP]
    VAC[i] = r.theta_vac["*"]
    ZZ[i], RN[i] = r.zrev, r.rate_net
    FIX[i], RES[i], NIT[i] = r.flux_fix, np.max(np.abs(r.residual)), r.n_iter

# ★ 从**中间**出发向两侧续算：中间点最好解，两端（强/弱吸附极限）最难。
#   初版从 i=0（E_N=−2，表面中毒、θ_vac~1e-30）起步，**第一个点就卡住**。
MID = n // 2
r0 = solve_at(MID, None)
if r0 is None:
    raise SystemExit("中间点就解不出来，先查口径")
store(MID, r0)
for direction in (+1, -1):
    guess = dict(r0.theta)
    i = MID + direction
    while 0 <= i < n:
        r = solve_at(i, guess)
        if r is None:
            fail.append(i)
            guess = None
        else:
            store(i, r)
            guess = dict(r.theta)
        i += direction
fail.sort()
el = time.time() - t0
ok = np.ones(n, bool); ok[fail] = False
print("   401 点求解完毕：%.1f s（%.0f ms/点），未收敛 %d 个 %s"
      % (el, el / n * 1e3, len(fail), fail[:10] if fail else ""))
print("   牛顿迭代中位 %d 次；缩放残差中位 %.2e 最大 %.2e"
      % (int(np.median(NIT[ok])), np.median(RES[ok]), RES[ok].max()))

def cmp(name, mine, ref, mask):
    rel = np.abs(np.where(np.abs(ref) > 0, mine / np.where(ref != 0, ref, 1) - 1,
                          np.abs(mine)))
    r = rel[mask]
    print("   %-12s 最大 %.3e   中位 %.3e   95%% 分位 %.3e"
          % (name, r.max(), np.median(r), np.percentile(r, 95)))
    return r

print("\n   逐点相对差（%d 个收敛点）：" % ok.sum())
r_th = cmp("y[H*,N*,NH*,NH2*]", TH, Y[:, :4], ok)
r_vac = cmp("y[空位]", VAC, Y[:, 4], ok)
r_z = cmp("Z (可逆度)", ZZ, Z, ok)
r_rn = cmp("Rnet (净速率)", RN, RNET, ok)

# ══════════ ②b Rnet 要换度量：它跨 300 个数量级 ══════════
print("\n══ ②b 净速率：跨 300 个数量级，相对差在 1e-300 附近没有意义 ══")
mag = np.abs(RNET[:, 4])
print("   真值 NH3 生成速率的量级分布：")
for lo, hi, nm in ((1e-30, 1e30, ">1e-30"), (1e-100, 1e-30, "1e-100~1e-30"),
                   (1e-200, 1e-100, "1e-200~1e-100"), (0.0, 1e-200, "<1e-200")):
    print("      %-16s %3d 个点" % (nm, int(((mag >= lo) & (mag < hi)).sum())))
meaningful = ok & (mag > 1e-250)
print("   ⇒ 只比速率在物理上有意义的点（>1e-250，%d 个）：" % meaningful.sum())
rm = np.abs(RN[meaningful] / np.where(RNET[meaningful] != 0, RNET[meaningful], 1) - 1)
print("      最大 %.3e   中位 %.3e   95%% 分位 %.3e"
      % (rm.max(), np.median(rm), np.percentile(rm, 95)))
both = ok[:, None] & (np.abs(RNET) > 0) & (np.abs(RN) > 0)
dlog = np.abs(np.log10(np.abs(RN[both])) - np.log10(np.abs(RNET[both])))
print("   ⇒ log10(速率) 的偏差（跨 300 数量级时这才是有意义的度量）：")
print("      最大 %.3e 个数量级   中位 %.3e" % (dlog.max(), np.median(dlog)))

print("\n══ ③ 通量投影起了多大作用 ══")
print("   flux_fix（原始解的通量不自洽程度）：中位 %.2e  最大 %.2e"
      % (np.median(FIX[ok]), FIX[ok].max()))
print("   flux_fix > 1e-8 的点：%d/%d —— 这些点若不投影，净速率就带着这个量级的错"
      % (int((FIX[ok] > 1e-8).sum()), int(ok.sum())))
fc = np.abs(2 * RNET[:, 1] / np.where(RNET[:, 4] != 0, RNET[:, 4], 1) - 1)
print("   真值自身的 N 通量守恒：最大偏差 %.2e（CATKINAS 用 1000 位算术）" % fc.max())
mk = ok & (RN[:, 4] != 0)
fcm = np.abs(2 * RN[mk][:, 1] / RN[mk][:, 4] - 1)
print("   本包投影后的 N 通量守恒：最大偏差 %.2e" % fcm.max())

print("\n══ ④ 偏差最大的几个点长什么样 ══")
relall = np.abs(RN / np.where(RNET != 0, RNET, 1) - 1).max(axis=1)
relall[~ok] = 0.0
worst = np.argsort(-relall)[:6]
print("   点号  E_N     Rnet相对差   真值Rnet[5]   1−Z最小    flux_fix")
for i in worst:
    omz = 1.0 - ZZ[i]
    print("   %3d   %+.2f   %.3e    %.3e   %.2e   %.2e"
          % (i, E1[i], relall[i], RNET[i, 4],
             omz[omz > 0].min() if (omz > 0).any() else 0.0, FIX[i]))

print("\n══ 判定 ══")
hard = max(r_th.max(), r_vac.max(), r_z.max())
print("   覆盖度 / 可逆度 最大相对差 = %.3e  %s" % (hard, "✅" if hard < 1e-8 else "⚠"))
print("   净速率（>1e-250 的点）    = %.3e  %s" % (rm.max(), "✅" if rm.max() < 1e-8 else "⚠"))
print("   log10(速率) 最大偏差      = %.3e 数量级  %s"
      % (dlog.max(), "✅" if dlog.max() < 1e-8 else "⚠"))
print("   收敛率 %d/401，耗时 %.1f s（%.0f ms/点）" % (int(ok.sum()), el, el / n * 1e3))
np.savez(_root.ROOT + "/notes/T7_catkinas401_result.npz",
         E1=E1, theta=TH, vac=VAC, zrev=ZZ, rate_net=RN, flux_fix=FIX,
         resid=RES, n_iter=NIT, ok=ok,
         ref_y=Y, ref_Z=Z, ref_Rnet=RNET, ref_kf=KF, ref_kr=KR)
print("   逐点结果已存 notes/T7_catkinas401_result.npz（供后续作图与 DRC 用）")
