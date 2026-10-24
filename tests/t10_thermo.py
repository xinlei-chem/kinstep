# -*- coding: utf-8 -*-
r"""T10 验证：气相热化学（Shomate / NIST）→ 熵校正。

目标是**钉死的**（台账 2026-10-22"决定：不装 CatMAP"条）：
  复现 CATKINAS 合成氨算例在 673 K 下的 5 个常数偏移
      ΔGa = [+1.077301, +1.504256, 0, 0, 0] eV
      ΔG⁰ = [+1.077301, +1.504256, 0, 0, −1.575974] eV
不是"整理一个热力学库"这种没有判据的说法。

数据来源：NIST Chemistry WebBook（Shomate），**不从二手手册抄表**。
"""
import sys, warnings
import _root, refdata  # noqa: F401  ← 包根入 sys.path + 第三方真值定位
warnings.filterwarnings("ignore")
import numpy as np
import scipy.io as sio

from kinstep import mech as M
from kinstep.thermo import self_check, entropy_offsets, DB, ThermoError, CITATION
from kinstep.scaling import from_table
from kinstep.scan import scan
from kinstep.units import K_B_EV
import kinstep.drivers

print("══ ① 数据自检：系数能否现算出各自的参考 S°/Cp°(298.15 K) ══")
ok, rows = self_check(verbose=True)
print("   ⇒ %s" % ("✓ 全过" if ok else "❌ 有对不上的 —— 多半是系数抄错"))

print("\n══ ② 温区：不许静默外推 ══")
try:
    DB["N2"]._range_at(50.0)
    print("   ❌ 没拦住")
except ThermoError as e:
    print("   ✓ 拦住了：%s" % str(e).splitlines()[0])
print("   N2 在 298 K 用对区间 = %.2f（文献 191.61）" % DB["N2"].S(298.15))
import math
bad = DB["N2"].ranges[1]          # 500–2000 K 段
print("   若错用 500–2000 K 段外推到 298 K = %.2f  ← 差 %.2f，看着挺像但它是错的"
      % (bad.S(298.15), bad.S(298.15) - 191.61))

print("\n══ ③ ★★ 复现 CATKINAS 的 5 个常数偏移（T10 的唯一判据）══")
mech = M.from_text("""
H2 + 2* <-> 2H*
N2 + 2* <-> 2N*
N* + H* <-> NH* + *
NH* + H* <-> NH2* + *
NH2* + H* <-> NH3 + 2*
""")
T = 673.0
# ★ 与 CATKINAS 对数时统一 J/mol→eV 的换算因子（第四处常数口径差异；
#   三个物种反解出的值散布只有 3.8e-10，确是它代码里的单一常数）
F_CATKINAS = 96485.2608
off_Ga, off_dG, notes = entropy_offsets(mech, T, J_per_mol_per_eV=F_CATKINAS)
off_Ga_own, off_dG_own, _ = entropy_offsets(mech, T)        # 本包自己的常数
d, REF_SRC = refdata.catkinas_vars("curve")
ref_Ga = (np.asarray(d["Ean"]) - np.asarray(d["Eao"]))[0]
ref_dG = (np.asarray(d["G0n"]) - np.asarray(d["G0o"]))[0]
print("   步  反应                      本包 ΔGa      CATKINAS     差")
for s in range(5):
    print("   %d   %-26s %+.9f  %+.9f  %.2e"
          % (s + 1, mech.steps[s].raw[:26], off_Ga[s], ref_Ga[s], abs(off_Ga[s] - ref_Ga[s])))
print("   步  反应                      本包 ΔG⁰      CATKINAS     差")
for s in range(5):
    print("   %d   %-26s %+.9f  %+.9f  %.2e"
          % (s + 1, mech.steps[s].raw[:26], off_dG[s], ref_dG[s], abs(off_dG[s] - ref_dG[s])))
eGa, edG = np.max(np.abs(off_Ga - ref_Ga)), np.max(np.abs(off_dG - ref_dG))
print("   ⇒ ΔGa 最大差 %.3e eV ；ΔG⁰ 最大差 %.3e eV   %s"
      % (eGa, edG, "✓ 复现" if max(eGa, edG) < 1e-8 else "❌ 对不上"))
e_own = max(np.max(np.abs(off_Ga_own - ref_Ga)), np.max(np.abs(off_dG_own - ref_dG)))
print("   （若用本包自己的 CODATA-2018 换算因子：最大差 %.3e eV —— 纯常数口径差异，"
      "换算因子相对差 7.39e-07）" % e_own)
for n in notes[1:]:
    print("  " + n)

print("\n══ ④ 缺数据时必须报错，不能按 0 处理 ══")
m2 = M.from_text("C2H6 + 2* <-> 2CH3*")
try:
    entropy_offsets(m2, 500.0)
    print("   ❌ 没报错")
except ThermoError as e:
    print("   ✓ 报错：%s" % str(e).splitlines()[0])
    print("     " + str(e).splitlines()[1].strip())
oG, odG, nt = entropy_offsets(m2, 500.0, strict=False)
print("   strict=False 时按 0 处理，但**说明里标出来**：%s" % nt[0])

print("\n══ ⑤ ★ 全链路：熵校正由热力学库自动算，复现 401 点 ══")
KT_CAT = 0.057994719627
A_PRE = 1.402306113158e13
T_EQ = T * KT_CAT / (K_B_EV * T)
TAB = {1: dict(bep=(0.54, 0.60), scale=(0.30, -0.68)),
       2: dict(bep=(0.78, 1.54), scale=(2.00,  0.00)),
       3: dict(bep=(0.26, 1.45), scale=(-0.15, -0.18)),
       4: dict(bep=(0.43, 1.66), scale=(-0.45, -0.45)),
       5: dict(bep=(0.31, 1.45), scale=(-0.85,  1.17))}
sc = from_table(mech, TAB, desc_names=["E_N"], offset_Ga=off_Ga, offset_dG=off_dG)
E1 = np.arange(-2.0, 2.0001, 0.01)
res = scan(mech, sc, dict(T_K=T_EQ), [E1], p_gas={"H2": 75.0, "N2": 25.0, "NH3": 1.0},
           A_pre=A_PRE, precision="auto", dps=60)
RNET, Y = np.asarray(d["Rnet"]), np.asarray(d["y"])
rel = np.abs(res.rate_of(-1) / RNET[:, 4] - 1)
relY = np.array([max(abs(res.theta[s][i] / Y[i, k] - 1)
                     for k, s in enumerate(["H*", "N*", "NH*", "NH2*"])) for i in range(401)])
print("   %s" % res.summary().splitlines()[1].strip())
print("   净速率最大相对差 %.3e ；覆盖度 %.3e" % (np.nanmax(rel), np.nanmax(relY)))
# ★ 这条判据的 1e-8 是**建立在高精度升级之上**的：强吸附支上 1−z 低于 double 地板，
#   没有 mpmath 就只能停在 ~1e-2（见正文 Fig.2 的 double-only 曲线）。
#   mpmath 是可选件，所以缺它时改用"双精度本来能达到的水平"做判据，并说清楚，
#   而不是报一个其实是"没装可选件"的失败。
#   （这条路以前走不到：本套总先因缺第三方参考数据而按 77 跳过。）
from kinstep.hp import available as _hp_ok          # noqa: E402
_tol, _why = (1e-8, "") if _hp_ok() else (
    1e-1, "（mpmath 未装 ⇒ 仅双精度，判据放宽到 1e-1）")
print("   ⇒ %s（T9 手工给偏移时是 2.705e-10）%s"
      % ("✓ 热力学库接上后全链路不变" if np.nanmax(rel) < _tol else "❌", _why))
print("\n   数据来源：%s" % CITATION)
