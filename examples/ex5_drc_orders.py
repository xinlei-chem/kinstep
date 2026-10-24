# -*- coding: utf-8 -*-
"""案例 5：速率控制度 DRC、反应级数、表观活化能

展示：四种 DRC 统一一个定义 X = ∂ln r/∂(−G_x/k_BT)
★ Campbell 求和规则 ΣX_B = 1 既是正确性检查，也用来**自动选差分步长**
"""
import numpy as np
from kinstep import mech as M
from kinstep.analyze import energetics_from_arrays, drc, calibrate_delta
import kinstep.drivers

mech = M.from_text("""
N2 + 2* <-> 2N*
H2 + 2* <-> 2H*
N* + H* <-> NH* + *
NH* + H* <-> NH2* + *
NH2* + H* <-> NH3 + 2*
""")
en = energetics_from_arrays(mech, np.array([3.044, 1.310, 1.403, 1.467, 1.813]),
                            np.array([1.504, 0.397, -0.18, -0.45, -0.406]))
state = dict(T_K=673.0)
PG = {"N2": 25.0, "H2": 75.0, "NH3": 1.0}

# ① 差分步长不靠经验值 —— 用 Campbell 求和规则自己选
best, rows = calibrate_delta(mech, en, state, p_gas=PG, verbose=True)
print("   ⇒ 自动选中 δ = %.1e eV" % best)

# ② 全套 DRC（Richardson 外推消去二阶截断误差）
d = drc(mech, en, state, p_gas=PG, delta=best, richardson=True)
print("\n" + d.report(target=4))
print("\n★ 怎么读这张表：")
top = int(np.argmax(d.X_B[:, 4])) + 1
print("   · X_B 最大的是第 %d 步 ⇒ 它是控速步（要提速就改它的势垒）" % top)
print("   · ΣX_B = %.12f —— 偏差 %.1e，这条规则不依赖机理也不依赖能量，"
      % (d.sum_B[4], abs(d.sum_B[4] - 1)))
print("     所以它**不对就别信结论**")
print("   · X_I 是中间体的热力学 DRC：负值表示稳定该中间体会**降低**速率")
print("   · 反应级数与表观活化能都是**实验上直接可测**的量，适合做 validation")
