# -*- coding: utf-8 -*-
"""显式物种前缀（Sp_ / E_）的回归测试 —— 老板 2026-10-22 提议。"""
import sys
import _root, refdata  # noqa: F401  ← 包根入 sys.path + 第三方真值定位
from kinstep import mech as M
from kinstep.naming import strip_species_prefix

def show(t, **kw):
    try:
        m = M.from_text(t, **kw); print(m.summary()); return m
    except Exception as e:
        print("  ⛔", str(e).replace("\n", "\n    ")); return None

print("══ ① 前缀剥离 ══")
for n in ("Sp_CO", "E_C", "E_CO2", "sp_Cat", "CO", "E_", "E_T"):
    print("   %-8s → %s" % (n, strip_species_prefix(n)))

print("\n══ ② ★ 混写不能产生两个物种（E_CO 与 CO 必须是同一个）══")
m = show("E_CO + * <-> CO*\nCO* + E_O* -> CO2 + 2*\nE_O2 + 2* <-> 2O*")
print("   物种表：", m.species_order)
print("   ⇒ 'CO' 只出现一次 =", m.species_order.count("CO") == 1,
      " 'O*' 只出现一次 =", m.species_order.count("O*") == 1)

print("\n══ ③ 前缀豁免保留名 —— 一个真的叫 T 的物种 ══")
print("  不带前缀（应被拦）：")
show("T + * <-> T*")
print("  带前缀（应放行）：")
m3 = show("Sp_T + * <-> Sp_T*", species_formula={"T": "H3"})   # 氚示意
if m3: print("   物种表：", m3.species_order)

print("\n══ ④ 非元素式名字：Cat / Dye ══")
print("  (a) 不带前缀 —— 解析不了（原来就这样）：")
show("Cat + * <-> Cat*")
print("  (b) 带前缀但没给组成 —— 应报错并给出两种修法：")
show("Sp_Cat + * <-> Sp_Cat*")
print("  (c) 带前缀 + species_formula —— 应通过，且元素配平仍然生效：")
m4 = show("Sp_Cat + * <-> Sp_Cat*", species_formula={"Cat": "C6H4O2"})
print("  (d) 组成给错导致不配平 —— 配平闸应抓住：")
show("Sp_Cat + * <-> Sp_Dye*", species_formula={"Cat": "C6H4O2", "Dye": "C6H4O3"})
print("  (e) exempt_balance 豁免 —— 应通过，但豁免必须在 summary 里可见：")
m5 = show("Sp_Nafion + * <-> Sp_Nafion*", exempt_balance=["Nafion"])

print("\n══ ⑤ 带电物种默认液相（Koper 坑的修复）══")
m6 = M.from_text("H2O(l) + * <-> OH* + H+ + e-   {beta_ct: 0.5}")
print("   H+ 的相态 =", m6.species["H+"].phase, " 电荷 =", m6.species["H+"].charge)
print("   ⇒ 液相物种：", m6.liq_species, "（以前 H+ 落在气相，p_gas 里找不到就静默取 0）")

print("\n══ ⑥ 原有机理不受影响（回归）══")
m7 = M.from_text("""
N2 + 2* <-> 2N*
H2 + 2* <-> 2H*
N* + H* <-> NH* + *
NH* + H* <-> NH2* + *
NH2* + H* <-> NH3 + 2*
""")
print("   合成氨 5 步仍然解析正常：%d 步 %d 物种" % (m7.n_step, len(m7.species_order)))
print("   ν 列 =", m7.species_order)
