# -*- coding: utf-8 -*-
"""案例 4：超声空化热点的多工况加权（超声催化）

★ 核心：r(⟨T⟩) ≠ ⟨r(T)⟩ —— 速率对温度指数敏感，**不能用平均工况代替逐工况求解**
"""
import numpy as np
from kinstep import mech as M
from kinstep.analyze import energetics_from_arrays
from kinstep.solve import steady_state, steady_state_multi
from kinstep.driver import DriverStack, make_driver
from kinstep.drivers.sono import calibrate_fraction_from_yield
import kinstep.drivers

mech = M.from_text("""
C2H6 + 2* <-> C2H5* + H*
C2H5* + * -> C2H4* + H*
C2H4* -> C2H4 + *
C2H5* -> C2H5 + *
2H* -> H2 + 2*
""")
en = energetics_from_arrays(mech, np.array([0.55, 1.95, 0.60, 0.85, 0.70]),
                            np.array([-0.20, 0.35, 0.25, 0.40, -0.30]))
PG = {"C2H6": 1.0, "C2H4": 1e-8, "C2H5": 1e-8, "H2": 1e-6}
T_BULK, T_HOT, F_HOT = 350.0, 5000.0, 1e-3
DS = DriverStack([make_driver("sono")])

print("体系：乙烷脱氢。裂解步 Ga=1.95 eV（高垒）vs C2H5 脱附 Ga=0.85 eV（低垒）")
res = steady_state_multi(mech, en, dict(T_K=T_BULK, hotspot_T_K=T_HOT,
                                        hotspot_fraction=F_HOT), DS,
                         p_gas=PG, tol=1e-12, strict_energy=False)
print("\n" + res.report(step=2))

T_avg = F_HOT * T_HOT + (1 - F_HOT) * T_BULK
r_avg = steady_state(mech, en, state=dict(T_K=T_avg), p_gas=PG,
                     tol=1e-12, strict_energy=False)
print("\n★ 对比两种做法（热点份额 %.0e，热点 %.0f K，本体 %.0f K）：" % (F_HOT, T_HOT, T_BULK))
print("   ⟨r(T)⟩ 逐工况加权 = %.6e   ← 对的做法" % res.rate_net[2])
print("   r(⟨T⟩) 用平均温度 %.2f K = %.6e   ← 常见的错做法" % (T_avg, r_avg.rate_net[2]))
print("   ⇒ 相差 %.3e 倍（%.1f 个数量级）" % (res.rate_net[2] / r_avg.rate_net[2],
                                             np.log10(res.rate_net[2] / r_avg.rate_net[2])))

r_hot = steady_state(mech, en, state=dict(T_K=T_HOT), p_gas=PG,
                     tol=1e-12, strict_energy=False)
r_bulk = steady_state(mech, en, state=dict(T_K=T_BULK), p_gas=PG,
                      tol=1e-12, strict_energy=False)
print("\n★ 交叉份额（热点何时开始主导）：f* = r_bulk/r_hot = %.3e"
      % (abs(r_bulk.rate_net[2]) / abs(r_hot.rate_net[2])))
print("   ⇒ 对 1.95 eV 这种高垒反应，空化份额小到 1e-27 热点就已主导。")
print("     这正是声化学能把常温下根本不走的反应推动起来的原因。")

print("\n★ 由实测产率反推热点份额（文献的剂量计路线，而不是猜）：")
r_meas = 0.004 * r_hot.rate_net[2] + 0.996 * r_bulk.rate_net[2]
w, msg = calibrate_fraction_from_yield(r_meas, r_hot.rate_net[2], r_bulk.rate_net[2])
print("   " + msg)
