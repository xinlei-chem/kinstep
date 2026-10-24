# -*- coding: utf-8 -*-
"""案例 6：连续流动反应器（CSTR）—— 气相组成也是未知数

微分反应器假定气相活度恒定（转化率极低时成立）；
要算**转化率随温度/压力怎么变**，就必须让气相随反应演化。
真值：CATKINAS 315 点 T–P 二维（见 tests/t9b_flow.py，315/315 收敛）
"""
import numpy as np
from kinstep import mech as M
from kinstep.analyze import energetics_from_arrays
from kinstep.solve import steady_state, FlowReactor
import kinstep.drivers

mech = M.from_text("""
H2 + 2* <-> 2H*
N2 + 2* <-> 2N*
N* + H* <-> NH* + *
NH* + H* <-> NH2* + *
NH2* + H* <-> NH3 + 2*
""")
Ga = np.array([1.310, 3.044, 1.403, 1.467, 1.813])
dG = np.array([0.397, 1.504, -0.18, -0.45, -0.406])
en = energetics_from_arrays(mech, Ga, dG)
P_TOT, T = 100.0, 900.0
p_in = {"H2": 0.75 * P_TOT, "N2": 0.25 * P_TOT, "NH3": 0.0}

print("进料 H2/N2/NH3 = %.1f / %.1f / %.1f bar，T = %.0f K"
      % (p_in["H2"], p_in["N2"], p_in["NH3"], T))

# ① 微分反应器：气相钉在进料 ⇒ 转化率恒为 0
r0 = steady_state(mech, en, state=dict(T_K=T),
                  p_gas={k: max(v, 1e-12) for k, v in p_in.items()},
                  tol=1e-12, strict_energy=False, precision="auto")
print("\n微分反应器（气相固定）：NH3 生成速率 %.4e /site/s，转化率恒为 0"
      % r0.rate_net[4])

# ② 流动反应器：气相是未知数
flow = FlowReactor(p_in=p_in, Fgt=0.02, Pn=1e-3)
r = steady_state(mech, en, state=dict(T_K=T), p_gas={k: max(v, 1e-12)
                                                     for k, v in p_in.items()},
                 flow=flow, tol=1e-12, strict_energy=False, precision="auto")
print("\n流动反应器（CSTR，气相是未知数）：%s" % ("收敛" if r.converged else "★未收敛"))
print("   出口 H2/N2/NH3 = %.4f / %.4f / %.6f bar"
      % (r.p_gas_out["H2"], r.p_gas_out["N2"], r.p_gas_out["NH3"]))
conv = (p_in["H2"] - r.p_gas_out["H2"]) / p_in["H2"]
print("   H2 转化率 = %.4f%%" % (100 * conv))
print("   求解路径：%s" % r.method)

# ③ 扫温度看转化率曲线（沿温度续算 —— 这是必需的，不是优化）
print("\n温度扫描（沿 T 续算）：")
print("   T(K)    出口 NH3 (bar)   H2 转化率")
guess = dict(r.theta); pg = dict(r.p_gas_out)
for Tx in (900, 800, 700, 600, 500):
    rr = steady_state(mech, en, state=dict(T_K=float(Tx)), p_gas=pg, flow=flow,
                      theta_guess=guess, ladder=("newton",), tol=1e-12,
                      strict_energy=False, precision="auto")
    if rr.converged:
        guess, pg = dict(rr.theta), dict(rr.p_gas_out)
        print("   %4d    %.6e     %.4f%%"
              % (Tx, rr.p_gas_out["NH3"],
                 100 * (p_in["H2"] - rr.p_gas_out["H2"]) / p_in["H2"]))
    else:
        print("   %4d    未收敛" % Tx)
