# -*- coding: utf-8 -*-
"""案例 3：光催化的光强幂律（光催化）

展示：准费米能级分裂 + 载流子收支自洽 → 光强幂律自己从 1 走到 0.5
★ 幂指数是**算出来的结果**，不是输入的假设 —— 这正是 MKM 相对经验幂律的价值
"""
import numpy as np
from kinstep import mech as M
from kinstep.analyze import energetics_from_arrays
from kinstep.drivers.photo import photo_steady
import kinstep.drivers

# 控速步只转移**一个**载流子 —— 经典半级定律隐含这个条件
mech = M.from_text("""
H2O(l) + 2* <-> OH* + H*
OH* -> OH(l) + *            {nu_hole_vb: 1}
2H* -> H2 + 2*              {nu_elec_cb: 2}
""")
print(mech.summary())
print("\n每轮 1 空穴 + 2×(1/2) 电子 ⇒ 电荷自动守恒（本包会自检）")

en = energetics_from_arrays(mech, np.array([0.20, 1.35, 0.50]),
                            np.array([-0.30, 0.25, -0.35]))
PG, CL = {"H2": 1e-6}, {"H2O": 1.0, "OH": 1e-9}

def st(flux):
    return dict(T_K=300.0, photon_flux=flux, absorptance=0.5,
                gamma_site_photo=1.5e19, photon_energy_eV=3.2, qfl_frac_n=0.5,
                recomb_rate_const=1.0e-17, photo_carrier_mode="steady_carrier")

flux = np.logspace(13, 26, 14)
rate, qfl, aqy, imb = [], [], [], []
for F in flux:
    r, x, info = photo_steady(mech, en, st(F), p_gas=PG, c_liq=CL, tol=1e-9)
    rate.append(abs(r.rate_net[1])); qfl.append(info["qfl_eV"])
    aqy.append(info["AQY"]); imb.append(info["charge_imbalance"])
rate = np.array(rate)
slope = np.gradient(np.log10(rate), np.log10(flux))
print("\n   光子流(m⁻²s⁻¹)   ΔE_qfl(eV)   速率(s⁻¹)      AQY        局域幂指数")
for i in range(0, len(flux), 2):
    print("   %.2e        %.4f      %.4e   %.3e   %.4f"
          % (flux[i], qfl[i], rate[i], aqy[i], slope[i]))
print("\n   低光强端幂指数 %.4f（文献：一级）" % slope[:3].mean())
print("   高光强端幂指数 %.4f（文献：半级，Turchi & Ollis 1990）" % slope[-3:].mean())
print("   AQY 从 %.3f 降到 %.2e —— 高光强下光子被复合浪费，正是半级的物理来源"
      % (aqy[0], aqy[-1]))
print("   电荷守恒自检：最大失衡 %.2e" % max(imb))
