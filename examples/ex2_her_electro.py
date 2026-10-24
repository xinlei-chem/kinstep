# -*- coding: utf-8 -*-
"""案例 2：HER 的电位扫描与 Tafel 斜率（电催化）

展示：electro driver（CHE/BV）→ 电位扫描 → 电流密度 → Tafel 斜率
★ 机理里显式写了 H+ ⇒ 电位必须用 SHE 标度，否则 pH 会被重复计算
"""
import numpy as np
from kinstep import mech as M
from kinstep.analyze import energetics_from_arrays, tafel_slope
from kinstep.solve import steady_state
from kinstep.driver import DriverStack, make_driver
import kinstep.drivers

mech = M.from_text("""
H+ + e- + * <-> H*            {beta_ct: 0.5}
H* + H+ + e- -> H2 + *        {beta_ct: 0.5}
2H* -> H2 + 2*
""")
print(mech.summary())
print("\n注意第 1、2 步自动数出 nu_elec —— 这个符号最容易写错，所以不让用户填")

# 能量：这是**只有你能提供的部分**（换成你自己体系的 DFT 值）
Ga = np.array([0.45, 0.60, 0.75])
dG = np.array([-0.10, -0.30, -0.20])
en = energetics_from_arrays(mech, Ga, dG)

el = DriverStack([make_driver("electro", U_ref="SHE", gamma_site=2.2e-5)])
US = np.linspace(-0.30, 0.0, 31)
j, guess = [], None
for U in US:
    r = steady_state(mech, en, state=dict(T_K=298.15, U_elec=U), driver=el,
                     p_gas={"H2": 1.0}, c_liq={"H+": 1.0},
                     theta_guess=guess, ladder=("newton",) if guess else ("ode", "newton"),
                     tol=1e-12, precision="auto", strict_energy=False)
    j.append(r.obs["j_curr"]); guess = dict(r.theta)
j = np.array(j)
ts = tafel_slope(US, j)
m = np.isfinite(ts) & (np.abs(j) > 0)
print("\n   U(V vs SHE)   j (A/m²)        Tafel (mV/dec)")
for i in range(0, len(US), 5):
    print("   %+.3f        %+.4e     %s" % (US[i], j[i],
          "%.1f" % (1000 * ts[i]) if m[i] else "—"))
print("\n   强极化区的 Tafel 斜率 ≈ %.1f mV/dec" % (1000 * np.median(ts[m][:8])))
print("   参考：Volmer 控速 ~120、Heyrovsky ~40、Tafel ~30 mV/dec")
