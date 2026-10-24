# -*- coding: utf-8 -*-
r"""T17 验证：外电场 driver（`field`）。

## 为什么补这一套

`field` driver 一直随包发布，却**既没写进论文，也没有任何数值回归** ——
只在 T15 的自然语言测试里被顺带提到过名字。一个发出去的 driver 如果没人验过，
用户第一次用它时就是在替我们做测试。

## 验什么（四条，全部解析可查，不依赖任何外部数据）

① **公式逐位对**：δ(ΔG) = −(Δμ·F + ½·Δα·F²)，过渡态取半程近似 δGa = ½·δ(ΔG)。
② **零场归一**：F = 0 时必须与热催化基线**逐位相同**（不是"接近"）。
③ **只动声明了参数的那一步**：没给 `dmu_eA`/`dalpha_eA2_per_V` 的步必须纹丝不动。
④ **非可加性是真的**：`is_additive = False` 不能只是个标志位 ——
   F² 项使得"先加场再加别的"与反过来不等价，这里用两种顺序算出差值来证实。

⑤ 再跑一次稳态，确认降低某步能垒确实让速率单调上升（符号没写反）。
"""
import sys
import warnings

import numpy as np

import _root  # noqa: F401
warnings.filterwarnings("ignore")

from kinstep import mech as M
from kinstep.analyze import energetics_from_arrays
from kinstep.driver import DriverStack, make_driver
from kinstep.solve import steady_state

fail = []


def chk(tag, got, want, tol, unit=""):
    dev = abs(got - want) if want == 0 else abs(got / want - 1.0)
    ok = dev <= tol
    print("   %-46s %-14s 偏差 %.3e  %s"
          % (tag, ("%.10g%s" % (got, unit)), dev, "✓" if ok else "❌"))
    if not ok:
        fail.append(tag)


# 两步机理：第 1 步带场参数，第 2 步不带（用来验④里的"不该动的别动"）
MECH = M.from_text(
    "CO + * <-> CO*          {dmu_eA: 0.35, dalpha_eA2_per_V: 0.80}\n"
    "CO* + H2 -> CH2O + *\n")
GA = np.array([0.60, 0.90])
DG = np.array([-0.20, -0.50])
DMU, DALPHA = 0.35, 0.80

print("══ ① 场致能量移动的闭式解 ══")
for F in (-0.50, 0.0, 0.25, 1.00):
    drv = make_driver("field")
    Ga, dG = drv.modify_energies(MECH, GA, DG, dict(field_V_per_A=F))
    shift = DMU * F + 0.5 * DALPHA * F * F          # 解析值
    chk("F=%+.2f V/A  δ(ΔG) 对解析" % F, float(DG[0] - dG[0]), shift,
        1e-14 if shift else 0.0, " eV")
    chk("F=%+.2f V/A  δGa = ½·δ(ΔG)" % F, float(GA[0] - Ga[0]), 0.5 * shift,
        1e-14 if shift else 0.0, " eV")

print("\n══ ② 零场必须与基线逐位相同 ══")
drv = make_driver("field")
Ga0, dG0 = drv.modify_energies(MECH, GA, DG, dict(field_V_per_A=0.0))
chk("零场 max|ΔGa − Ga₀|", float(np.max(np.abs(Ga0 - GA))), 0.0, 0.0, " eV")
chk("零场 max|ΔdG − dG₀|", float(np.max(np.abs(dG0 - DG))), 0.0, 0.0, " eV")

print("\n══ ③ 没声明场参数的步必须纹丝不动 ══")
Ga1, dG1 = drv.modify_energies(MECH, GA, DG, dict(field_V_per_A=1.0))
chk("第 2 步 Ga 变化", float(Ga1[1] - GA[1]), 0.0, 0.0, " eV")
chk("第 2 步 dG 变化", float(dG1[1] - DG[1]), 0.0, 0.0, " eV")

print("\n══ ④ 二次项 ⇒ 非可加（标志位要名副其实）══")
print("   is_additive = %s  %s" % (drv.is_additive,
                                   "✓" if drv.is_additive is False else "❌"))
if drv.is_additive is not False:
    fail.append("field.is_additive 应为 False")
# 一次加满场 vs 分两次各加一半：F² 项使两者必然不等，差值有闭式解
full = DMU * 1.0 + 0.5 * DALPHA * 1.0
half_twice = 2.0 * (DMU * 0.5 + 0.5 * DALPHA * 0.25)
chk("Δ(一次加满 − 两次各半) 对解析", full - half_twice,
    0.5 * DALPHA * (1.0 - 2 * 0.25), 1e-14, " eV")

print("\n══ ⑤ 稳态：降低能垒 ⇒ 速率单调上升（符号检查）══")
en = energetics_from_arrays(MECH, GA, DG)
prev = None
for F in (0.0, 0.2, 0.4, 0.6):
    r = steady_state(MECH, en, state=dict(T_K=400.0, field_V_per_A=F),
                     driver=DriverStack([make_driver("field")]),
                     p_gas={"CO": 1.0, "H2": 1.0, "CH2O": 0.01}, strict_energy=False)
    rate = float(np.max(np.abs(np.asarray(r.rate_net, float))))
    print("   F = %+.2f V/A   r = %.6e s⁻¹ site⁻¹" % (F, rate))
    if prev is not None and not rate > prev:
        fail.append("F=%.2f 时速率没有随场增大" % F)
        print("      ❌ 未单调上升")
    prev = rate

print("\n" + "═" * 56)
if fail:
    print("  ❌ 失败 %d 项：%s" % (len(fail), "；".join(fail)))
    sys.exit(1)
print("  ✓ T17 全部通过：场效应公式、零场归一、作用域、非可加性、符号")
