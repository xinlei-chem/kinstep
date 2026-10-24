# -*- coding: utf-8 -*-
r"""T9 验证：scaling 层 + 1D/2D 扫描 + 火山图。

真值全部来自 CATKINAS（不需要别的程序，理由见台账 2026-10-22"决定：不装 CatMAP"）：
  ① scaling `ΔH⁰ = a·E_N + b` 与 BEP `Ea = α·ΔH⁰ + β` 的系数 ← `INPUT_curve.m` 的能量表
  ② 401 点的 `Eao/G0o`（原始层）与 `Rnet/y/Z`（求解结果）← `result_curve/data1.mat`
  ③ 熵校正是逐步常数偏移（实测跨度 < 1.4e-15）

★ 最硬的检验：**从 INPUT 的 4 个系数出发，不碰 .mat 里的能量，
  要能原地复现 T7 那 401 点的速率曲线** —— 这把 scaling + BEP + 偏移 + 求解器全链路串起来了。
"""
import sys, time, warnings
import _root, refdata  # noqa: F401  ← 包根入 sys.path + 第三方真值定位
warnings.filterwarnings("ignore")
import numpy as np
import scipy.io as sio

from kinstep import mech as M
from kinstep.scaling import from_table, fit_scaling
from kinstep.scan import scan
from kinstep.units import K_B_EV
from kinstep import viz
import kinstep.drivers

REF, REF_SRC = refdata.catkinas_vars("curve")
OUT = _root.ROOT + "/notes/T9_figures"
import os; os.makedirs(OUT, exist_ok=True)

KT_CAT = 0.057994719627
A_PRE = 1.402306113158e13   # 从真值反解；CATKINAS 的 e/h 也是老常数（见台账 T9 条）
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
Eao, G0o, Ean, G0n = (np.asarray(d[k]) for k in ("Eao", "G0o", "Ean", "G0n"))
RNET, Y, Z = (np.asarray(d[k]) for k in ("Rnet", "y", "Z"))
E1 = np.arange(-2.0, 2.0001, 0.01)
SP = ["H*", "N*", "NH*", "NH2*"]

# INPUT_curve.m 的能量表：[ Ea的BEP: α β | H0的scaling: a b ]
TAB = {1: dict(bep=(0.54, 0.60), scale=(0.30, -0.68)),
       2: dict(bep=(0.78, 1.54), scale=(2.00,  0.00)),
       3: dict(bep=(0.26, 1.45), scale=(-0.15, -0.18)),
       4: dict(bep=(0.43, 1.66), scale=(-0.45, -0.45)),
       5: dict(bep=(0.31, 1.45), scale=(-0.85,  1.17))}

print("══ ① scaling + BEP 两层对 401 点（原始层 Eao/G0o）══")
sc0 = from_table(mech, TAB, desc_names=["E_N"])
eg = [np.max(np.abs(sc0.energies(x)[1] - G0o[i])) for i, x in enumerate(E1)]
ea = [np.max(np.abs(sc0.energies(x)[0] - Eao[i])) for i, x in enumerate(E1)]
print("   ΔG⁰ 最大绝对差 %.3e ；Ga 最大绝对差 %.3e   %s"
      % (max(eg), max(ea), "✅ 逐点精确" if max(max(eg), max(ea)) < 1e-12 else "✗"))

print("\n══ ② 热化学偏移：从 .mat 测出来（T10 之前手工给）══")
off_Ga = (Ean - Eao)[0]
off_dG = (G0n - G0o)[0]
print("   ΔGa = %s" % np.array2string(off_Ga, precision=6))
print("   ΔG⁰ = %s" % np.array2string(off_dG, precision=6))
print("   401 点上的跨度：%.2e / %.2e ⇒ 确是常数（与描述符无关）"
      % (np.ptp(Ean - Eao, axis=0).max(), np.ptp(G0n - G0o, axis=0).max()))

print("\n══ ③ ★★ 全链路：只用 INPUT 的 4 个系数 + 偏移，复现 T7 的 401 点 ══")
sc = from_table(mech, TAB, desc_names=["E_N"], offset_Ga=off_Ga, offset_dG=off_dG)
# 核对加完偏移后与 Ean/G0n 一致
eg2 = [np.max(np.abs(sc.energies(x)[1] - G0n[i])) for i, x in enumerate(E1)]
print("   加偏移后 ΔG⁰ 对 G0n：最大 %.3e" % max(eg2))
t0 = time.time()
res = scan(mech, sc, dict(T_K=T_EQ), [E1], p_gas=P, A_pre=A_PRE,
           precision="auto", dps=60)
print("   " + res.summary().replace("\n", "\n   "))
relR = np.abs(res.rate_of(-1) / RNET[:, 4] - 1)
relY = np.array([max(abs(res.theta[s][i] / Y[i, k] - 1) for k, s in enumerate(SP))
                 for i in range(401)])
relZ = np.abs(res.zrev / np.where(Z != 0, Z, 1) - 1).max(axis=1)
print("   逐点对 CATKINAS：")
print("     净速率 Rnet  最大 %.3e  中位 %.3e" % (np.nanmax(relR), np.nanmedian(relR)))
print("     覆盖度 y     最大 %.3e  中位 %.3e" % (np.nanmax(relY), np.nanmedian(relY)))
print("     可逆度 Z     最大 %.3e  中位 %.3e" % (np.nanmax(relZ), np.nanmedian(relZ)))
print("   ⇒ %s（对比 T7 直接喂 Ean/G0n：2.706e-10）"
      % ("✅ 全链路复现" if np.nanmax(relR) < 1e-8 else "⚠ 见上"))

print("\n══ ④ 2D 扫描（造第二个描述符：让步 3/4/5 也依赖 E_H）══")
TAB2 = {1: dict(bep=(0.54, 0.60), scale=((0.30, 1.00), -0.68)),
        2: dict(bep=(0.78, 1.54), scale=((2.00, 0.00),  0.00)),
        3: dict(bep=(0.26, 1.45), scale=((-0.15, -0.50), -0.18)),
        4: dict(bep=(0.43, 1.66), scale=((-0.45, -0.50), -0.45)),
        5: dict(bep=(0.31, 1.45), scale=((-0.85, -0.50),  1.17))}
sc2 = from_table(mech, TAB2, desc_names=["E_N", "E_H"],
                 offset_Ga=off_Ga, offset_dG=off_dG)
ax_N = np.linspace(-1.5, 1.5, 31)
ax_H = np.linspace(-0.4, 0.4, 21)
res2 = scan(mech, sc2, dict(T_K=T_EQ), [ax_N, ax_H], p_gas=P, A_pre=A_PRE,
            precision="auto", dps=60)
print("   " + res2.summary().replace("\n", "\n   "))
# ★ 退化检验：E_H = 0 的那一列应当与 1D 结果一致（scaling 在 E_H=0 处退化为 1D）
j0 = int(np.argmin(np.abs(ax_H)))
print("   ★ 退化检验：E_H = %.3f 这一列 vs 1D 扫描同 E_N 处" % ax_H[j0])
r1d = np.array([res.rate_of(-1)[np.argmin(np.abs(E1 - v))] for v in ax_N])
r2d = res2.rate_of(-1)[:, j0]
m = np.isfinite(r1d) & np.isfinite(r2d) & (r1d != 0)
# ★ 1e-6 这个判据**建立在高精度升级之上**：强吸附支上 1−z 低于 double 地板，
#   没装 mpmath 时两条路径只能一致到 ~1e-4（纯粹是双精度相消，不是扫描逻辑的问题）。
#   mpmath 是可选件 ⇒ 缺它时按双精度能达到的水平判，并把原因印出来。
#   （参考值随包分发之前，本套总先因缺 .mat 而按 77 跳过，这条路走不到。）
from kinstep.hp import available as _hp_ok          # noqa: E402
_tol9, _why9 = (1e-6, "") if _hp_ok() else (
    1e-3, "（mpmath 未装 ⇒ 仅双精度，判据放宽到 1e-3）")
_dev9 = float(np.max(np.abs(r2d[m] / r1d[m] - 1)))
print("     最大相对差 %.3e  %s%s"
      % (_dev9, "✅ 2D 在 E_H=0 处退化为 1D" if _dev9 < _tol9 else "✗", _why9))

print("\n══ ⑤ scaling 拟合（给自己有 DFT 数据的用户）══")
rng = np.random.default_rng(7)
dsamp = rng.uniform(-1.5, 1.5, 12)
gsamp = 2.0 * dsamp + 0.0 + rng.normal(0, 0.03, 12)
a, b, r2 = fit_scaling(dsamp, gsamp)
print("   真参数 a=2.00 b=0.00（加 0.03 eV 噪声）⇒ 拟合 a=%.4f b=%+.4f  R²=%.6f"
      % (a[0], b, r2))

print("\n══ ⑥ 作图 + 导出原始数据 ══")
print("   matplotlib 可用 =", viz.available())
# ★ matplotlib 是**可选件**：没装就只导数据、不出图，不能让整套回归炸掉。
#   这条路以前走不到 —— 本套总是先因为缺第三方参考数据而按 77 跳过；
#   参考值随包分发之后它才真的会在"只有 numpy+scipy"的干净安装里跑起来。
if viz.available():
    viz.volcano_1d(res, path=OUT + "/volcano_1d.png")
    viz.coverage_1d(res, path=OUT + "/coverage_1d.png")
    viz.volcano_2d(res2, path=OUT + "/volcano_2d.png")
else:
    print("   ⊘ 跳过出图（matplotlib 未装，属可选件）；下面的 .dat 仍然照导")
res.save_dat(OUT + "/scan_1d.dat")
res2.save_dat(OUT + "/scan_2d.dat")
import os
for f in sorted(os.listdir(OUT)):
    print("   %-18s %8.1f KB" % (f, os.path.getsize(os.path.join(OUT, f)) / 1024))
