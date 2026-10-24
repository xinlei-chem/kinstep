# -*- coding: utf-8 -*-
"""案例 1：合成氨的描述符扫描与火山图（热催化）

展示：机理解析 → scaling/BEP → 熵校正 → 1D 扫描 → 火山图 + 原始数据
真值：CATKINAS 401 点（本案例的数值已逐点对过，见 tests/t9_scan.py）
"""
import numpy as np
from kinstep import mech as M
from kinstep.scaling import from_table
from kinstep.thermo import entropy_offsets
from kinstep.scan import scan
from kinstep import viz

# ① 机理：照着纸上写。四道闸自动跑（元素/电荷配平、位点守恒、撞名、箭头）
mech = M.from_text("""
N2 + 2* <-> 2N*
H2 + 2* <-> 2H*
N* + H* <-> NH* + *
NH* + H* <-> NH2* + *
NH2* + H* <-> NH3 + 2*
""")
print(mech.summary())

# ② 气相熵校正：由 NIST Shomate 自动算（673 K 下 H2 这一项就有 1.08 eV，不能省）
T = 673.0
off_Ga, off_dG, notes = entropy_offsets(mech, T)
print("\n熵校正（eV）：ΔGa =", np.round(off_Ga, 4), " ΔG⁰ =", np.round(off_dG, 4))

# ③ scaling + BEP：描述符 E_N → 逐步能量
#    ΔG⁰ = a·E_N + b ；Ga = α·ΔG⁰ + β
sc = from_table(mech, {
    1: dict(bep=(0.78, 1.54), scale=(2.00,  0.00)),   # N2 解离
    2: dict(bep=(0.54, 0.60), scale=(0.30, -0.68)),   # H2 解离
    3: dict(bep=(0.26, 1.45), scale=(-0.15, -0.18)),
    4: dict(bep=(0.43, 1.66), scale=(-0.45, -0.45)),
    5: dict(bep=(0.31, 1.45), scale=(-0.85,  1.17)),  # NH3 脱附
}, desc_names=["E_N"], offset_Ga=off_Ga, offset_dG=off_dG)
print("\n" + sc.describe())

# ④ 扫描：从网格中心向两侧续算；precision="auto" 在需要时自动升到高精度
res = scan(mech, sc, dict(T_K=T), [np.linspace(-2.0, 2.0, 201)],
           p_gas={"N2": 25.0, "H2": 75.0, "NH3": 1.0}, precision="auto")
print("\n" + res.summary())

# ⑤ 出图 + 导出可作图的原始数据（图给人看，数据给后续分析用）
if viz.available():
    viz.volcano_1d(res, path="ex1_volcano.png")
    viz.coverage_1d(res, path="ex1_coverage.png")
    print("\n图 → ex1_volcano.png / ex1_coverage.png")
res.save_dat("ex1_data.dat")
print("数据 → ex1_data.dat（Origin / gnuplot 直接读）")
