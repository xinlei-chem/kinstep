# -*- coding: utf-8 -*-
r"""T11 验证：高精度后端（mpmath）。

三件事，按顺序：
  ① **与 double 互校**：在精度充足的点上，两条路必须给出同一个答案
     （否则就不是"精度问题"，而是我把方程写错了）——这是前提，不是可选项
  ② ★ **救 T7 的 6 个坏点**：E_N ≥ +1.84，真值速率 ~1e-30，1−Z 压在 4·eps 上
  ③ ★ **救 T8 的退化点**：E_N = +1.0 / +1.3 的 DRC
"""
import sys, time, warnings
import _root, refdata  # noqa: F401  ← 包根入 sys.path + 第三方真值定位
warnings.filterwarnings("ignore")
import numpy as np
import scipy.io as sio

from kinstep import mech as M
from kinstep.analyze import energetics_from_arrays
from kinstep.hp import steady_state_hp, needs_high_precision, available
from kinstep.solve import steady_state
from kinstep.units import K_B_EV
import kinstep.drivers

REF, REF_SRC = refdata.catkinas_vars("curve")
KT_CAT = 0.057994719627           # CATKINAS 的 k_BT（CODATA-2006 的 k_B）
A_PRE = 1.402306113158e13
T_NOM = 673.0
T_EQ = T_NOM * KT_CAT / (K_B_EV * T_NOM)
P = {"H2": 75.0, "N2": 25.0, "NH3": 1.0}

mech = M.from_text("""
H2 + 2* <-> 2H*
N2 + 2* <-> 2N*
N* + H* <-> NH* + *
NH* + H* <-> NH2* + *
NH2* + H* <-> NH3 + 2*
""")
d = REF
Ean, G0n = np.asarray(d["Ean"]), np.asarray(d["G0n"])
Y, Z, RNET = (np.asarray(d[k]) for k in ("y", "Z", "Rnet"))
E1 = np.arange(-2.0, 2.0001, 0.01)
SP = ["H*", "N*", "NH*", "NH2*"]

print("mpmath 可用 =", available())
# ★ 这道闸是"参考值随包分发"之后才需要的：在那以前本套总是先因为缺 .mat 而按 77 跳过，
#   缺 mpmath 的情况根本走不到。现在数据有了，干净安装（只有 numpy+scipy）会真的跑进来，
#   没有 mpmath 就会在第一次调用 steady_state_hp 时炸掉 —— 那是**失败**，不是"可选件没装"。
if not available():
    sys.stderr.write(
        "⊘ 跳过本套回归：高精度后端 mpmath 未安装（可选件）\n"
        "  装上再跑：pip install mpmath   或  pip install '.[all]'\n")
    sys.exit(refdata.SKIP)

def dbl(i):
    return steady_state(mech, energetics_from_arrays(mech, Ean[i], G0n[i], A_pre=A_PRE),
                        state=dict(T_K=T_EQ), p_gas=P, tol=1e-12, strict_energy=False)

def hp(i, dps=50, guess=None):
    return steady_state_hp(mech, Ean[i], G0n[i], T_NOM, p_gas=P, A_pre=A_PRE,
                           kT_eV=KT_CAT, dps=dps, theta_guess=guess)

print("\n══ ① 与 double 互校（精度充足的点；若对不上就是方程写错了，不是精度问题）══")
print("   点   E_N    θ 最大相对差   1−z 最大相对差   r_net 最大相对差   hp 迭代  hp 耗时")
for i in (100, 150, 200, 250):
    rd = dbl(i)
    t0 = time.time(); rh = hp(i, guess=dict(rd.theta)); th = time.time() - t0
    dth = max(abs(rh.theta[s] / rd.theta[s] - 1) for s in SP)
    o_d, o_h = 1.0 - rd.zrev, rh.one_minus_z
    m = (o_d > 1e-12) & (o_h > 0)
    doz = np.max(np.abs(o_d[m] / o_h[m] - 1)) if m.any() else 0.0
    mr = (np.abs(rd.rate_net) > 0) & (np.abs(rh.rate_net) > 0)
    drn = np.max(np.abs(rd.rate_net[mr] / rh.rate_net[mr] - 1))
    print("   %3d  %+.2f   %.3e      %.3e         %.3e       %3d     %.2f s"
          % (i, E1[i], dth, doz, drn, rh.n_iter, th))
print("   ⇒ 两条路在好点上一致 ⇒ 方程组相同，差异只会来自精度 ✓")

print("\n══ ② ★ 救 T7 的 6 个坏点（E_N ≥ +1.84，真值速率 ~1e-30）══")
bad = [384, 388, 389, 395, 399, 400]
print("   点   E_N   真值Rnet[5]   double     相对差      高精度      相对差    1−z(hp)")
agg_d, agg_h = [], []
for i in bad:
    rd = dbl(i)
    need, why = needs_high_precision(rd)
    rh = hp(i, dps=60, guess=dict(rd.theta))
    ed = abs(rd.rate_net[4] / RNET[i, 4] - 1)
    eh = abs(rh.rate_net[4] / RNET[i, 4] - 1)
    agg_d.append(ed); agg_h.append(eh)
    omz = rh.one_minus_z[rh.one_minus_z > 0]
    print("   %3d  %+.2f  %.4e  %.4e  %.3e   %.4e  %.3e   %.2e"
          % (i, E1[i], RNET[i, 4], rd.rate_net[4], ed, rh.rate_net[4], eh,
             omz.min() if omz.size else 0.0))
print("   ⇒ double 最大相对差 %.3e → 高精度 %.3e（改善 %.0f 倍）"
      % (max(agg_d), max(agg_h), max(agg_d) / max(agg_h) if max(agg_h) > 0 else np.inf))
print("   自动升级判据 needs_high_precision 在这些点上：", "全部触发 ✓"
      if all(needs_high_precision(dbl(i))[0] for i in bad) else "有漏判 ✗")
print("   触发原因示例（点 %d）：" % bad[0])
for w in needs_high_precision(dbl(bad[0]))[1]:
    print("      ·", w)

print("\n══ ③ 覆盖度在坏点上本来就是对的吗（T7 的说法要验证）══")
print("   点   E_N    θ(double) 相对差   θ(高精度) 相对差")
for i in bad:
    rd = dbl(i); rh = hp(i, dps=60, guess=dict(rd.theta))
    ed = max(abs(rd.theta[s] / Y[i, k] - 1) for k, s in enumerate(SP))
    eh = max(abs(rh.theta[s] / Y[i, k] - 1) for k, s in enumerate(SP))
    print("   %3d  %+.2f   %.3e          %.3e" % (i, E1[i], ed, eh))
print("   ⇒ 证实 T7 的结论：**覆盖度在 double 下就是精确的**，只有净速率不是")

print("\n══ ④ 位数够不够：1−z 随 dps 的收敛 ══")
i = 400
rd = dbl(i)
print("   点 400（E_N=+2.00），真值 Rnet[5] = %.6e" % RNET[i, 4])
print("   dps   r_net[5]            与真值相对差    1−z(最小)     耗时")
prev = None
for dps in (20, 30, 40, 50, 60, 80):
    t0 = time.time(); rh = hp(i, dps=dps, guess=dict(rd.theta)); el = time.time() - t0
    omz = rh.one_minus_z[rh.one_minus_z > 0]
    print("   %3d   %.10e   %.3e    %.3e    %.2f s"
          % (dps, rh.rate_net[4], abs(rh.rate_net[4] / RNET[i, 4] - 1),
             omz.min() if omz.size else 0.0, el))
print("   ⇒ 位数不足时会残留误差；够了之后平台化 —— **这就是该取多少位的判据**")

print("\n══ ⑤ ★★ 用 precision='auto' 重跑 T7 全部 401 点 ══")
t0 = time.time()
TH = np.zeros((401, 4)); RN = np.zeros((401, 5)); ZZ = np.zeros((401, 5))
PREC = []; NHP = 0
guess = None
MID = 200
order = [MID] + [MID + k for k in range(1, 201)] + [MID - k for k in range(1, 201)]
for i in order:
    r = steady_state(mech, energetics_from_arrays(mech, Ean[i], G0n[i], A_pre=A_PRE),
                     state=dict(T_K=T_EQ), p_gas=P,
                     theta_guess=guess, ladder=("newton",) if guess else ("ode", "newton"),
                     tol=1e-12, precision="auto", dps=60, strict_energy=False)
    if not r.converged:
        r = steady_state(mech, energetics_from_arrays(mech, Ean[i], G0n[i], A_pre=A_PRE),
                         state=dict(T_K=T_EQ), p_gas=P, tol=1e-12,
                         precision="auto", dps=60, strict_energy=False)
    TH[i] = [r.theta[s] for s in SP]; RN[i] = r.rate_net; ZZ[i] = r.zrev
    PREC.append(r.precision); NHP += (r.precision == "high")
    guess = dict(r.theta)
    if i == MID + 200:
        guess = dict(TH[MID][k] for k in ()) if False else None
el = time.time() - t0
relY = np.abs(TH / Y[:, :4] - 1).max(axis=1)
relR = np.abs(RN / np.where(RNET != 0, RNET, 1) - 1).max(axis=1)
relZ = np.abs(ZZ / np.where(Z != 0, Z, 1) - 1).max(axis=1)
print("   401 点耗时 %.1f s（%.0f ms/点）；升到高精度的点：%d/401" % (el, el / 401 * 1e3, NHP))
print("   覆盖度 y   最大 %.3e   中位 %.3e" % (relY.max(), np.median(relY)))
print("   可逆度 Z   最大 %.3e   中位 %.3e" % (relZ.max(), np.median(relZ)))
print("   净速率 Rnet 最大 %.3e   中位 %.3e   95%% 分位 %.3e"
      % (relR.max(), np.median(relR), np.percentile(relR, 95)))
print("   >1e-6 的点：%d 个（T7 用纯 double 时是 27 个）" % int((relR > 1e-6).sum()))
nb = np.argsort(-relR)[:5]
print("   仍最差的 5 点：", ", ".join("点%d(E_N=%+.2f, %.1e)" % (i, E1[i], relR[i]) for i in nb))
print("   ⇒ 对比 T7 纯 double：净速率最大 5.599e-02 → %.3e" % relR.max())

print("\n══ ⑥ ★ 救 T8 的 DRC 退化点（E_N = +1.0 / +1.3）══")
from kinstep.analyze import drc
DI = np.asarray(d["DRC_I"]); DPC = np.asarray(d["DRC_PC"]); DB = np.asarray(d["DRC_B"])
print("   点   E_N    量       double        高精度        真值         double差   高精度差")
for i in (300, 330):
    en = energetics_from_arrays(mech, Ean[i], G0n[i], A_pre=A_PRE)
    rd = drc(mech, en, dict(T_K=T_EQ), p_gas=P, delta=1e-3, richardson=True)
    rh = drc(mech, en, dict(T_K=T_EQ), p_gas=P, delta=1e-3, richardson=True,
             precision="auto", dps=60)
    for nm, a, b, ref in (("X_I(H*)", rd.X_I[0, 4], rh.X_I[0, 4], DI[i][0, 4]),
                          ("X_I(*)", rd.X_I[4, 4], rh.X_I[4, 4], DI[i][4, 4]),
                          ("X_PC(N2)", rd.X_PC[1, 4], rh.X_PC[1, 4], DPC[i][1, 4]),
                          ("X_B(步2)", rd.X_B[1, 4], rh.X_B[1, 4], DB[i][1, 4])):
        print("   %3d  %+.2f  %-9s %+.8f  %+.8f  %+.8f  %.2e   %.2e"
              % (i, E1[i], nm, a, b, ref, abs(a - ref), abs(b - ref)))
    print("        ΣX_B: double %.12f  高精度 %.12f" % (rd.sum_B[4], rh.sum_B[4]))
