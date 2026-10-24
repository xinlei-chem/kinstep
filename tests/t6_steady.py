# -*- coding: utf-8 -*-
r"""T6 验证：稳态求解器。

四套验证。**残差小只说明"自洽"，不说明"对"** —— 所以前两套是独立真值：
  ① **解析闭式解**：甲酸分解（单表面物种），θ 与 r 都有手写闭式解
  ② **长时瞬态极限**：完全独立的数值路径（ODE 积到稳）
  ③ 合成氨 5 步：位点守恒 + **元素通量守恒**（线性机理的免费自检）
  ④ ★ 准平衡极限：把第 1 步推到 1−z ~ 1e-13，**与闭式解逐档对比**，
     看求解器在哪一档开始失真、诊断有没有诚实报出来
"""
import sys, math
import _root, refdata  # noqa: F401  ← 包根入 sys.path + 第三方真值定位
import numpy as np
from kinstep import mech as M
from kinstep.rates import Energetics, compile_rates, apply_state
from kinstep.solve import steady_state, transient, Constant
import kinstep.drivers

# ── 甲酸分解：HCOOH + * <-> HCOOH* ;  HCOOH* -> CO2 + H2 + * ──
#    稳态： k1·p·θv − k₋₁·θ − k2·θ = 0 ,  θ + θv = 1
#      ⇒ θ = k1·p / (k1·p + k₋₁ + k2) ,  r = k2·θ
TXT1 = """
HCOOH + * <-> HCOOH*
HCOOH* -> CO2 + H2 + *
"""
m1 = M.from_text(TXT1)
T, p_HA = 500.0, 0.25

def closed_form(en, T, p):
    cr = compile_rates(en); rc = apply_state(cr, T)
    k1, km1, k2 = rc.kf[0], rc.kr[0], rc.kf[1]
    th = k1 * p / (k1 * p + km1 + k2)
    return th, k2 * th, (k1, km1, k2)

print("══ ① 解析闭式解：甲酸分解（单表面物种）══")
en1 = Energetics(m1, {1: dict(Ga_fwd=0.60, dG_rxn=-0.30),
                      2: dict(Ga_fwd=1.10, dG_rxn=-2.0)})
th_x, r_x, (k1, km1, k2) = closed_form(en1, T, p_HA)
print("   kf1=%.6e  kr1=%.6e  kf2=%.6e" % (k1, km1, k2))
r1 = steady_state(m1, en1, state=dict(T_K=T),
                  p_gas={"HCOOH": p_HA, "CO2": 0.0, "H2": 0.0})
print("   θ(HCOOH*) 本包 %.16e" % r1.theta["HCOOH*"])
print("             闭式 %.16e   相对差 %.2e" % (th_x, abs(r1.theta["HCOOH*"]/th_x - 1)))
print("   r_net(步2) 本包 %.16e" % r1.rate_net[1])
print("             闭式 %.16e   相对差 %.2e" % (r_x, abs(r1.rate_net[1]/r_x - 1)))
print("   位点守恒误差 %.2e ；缩放残差 %.2e" % (max(r1.site_error.values()),
                                               np.max(np.abs(r1.residual))))
print("   路径：%s（%d 次残差求值）" % (r1.method, r1.n_iter))

print("\n══ ② 长时瞬态极限（独立数值路径）══")
tr = transient(m1, en1, {}, (0.0, 1.0), protocol=Constant(T_K=T),
               p_gas={"HCOOH": p_HA, "CO2": 0.0, "H2": 0.0},
               n_eval=5, method="BDF", rtol=1e-12, atol=1e-18)
print("   ODE 末值 θ = %.16e" % tr.theta["HCOOH*"][-1])
print("   稳态解   θ = %.16e   相对差 %.2e"
      % (r1.theta["HCOOH*"], abs(tr.theta["HCOOH*"][-1]/r1.theta["HCOOH*"] - 1)))

print("\n══ ③ 合成氨 5 步（CATKINAS 算例的机理，T7 要用）══")
m3 = M.from_text("""
N2 + 2* <-> 2N*
H2 + 2* <-> 2H*
N* + H* <-> NH* + *
NH* + H* <-> NH2* + *
NH2* + H* <-> NH3 + 2*
""")
en3 = Energetics(m3, {1: dict(Ga_fwd=1.50, dG_rxn=-0.60),
                      2: dict(Ga_fwd=0.30, dG_rxn=-0.80),
                      3: dict(Ga_fwd=1.00, dG_rxn=0.20),
                      4: dict(Ga_fwd=1.05, dG_rxn=0.15),
                      5: dict(Ga_fwd=0.95, dG_rxn=-0.10)})
st3 = dict(T_K=700.0); pg3 = {"N2": 0.25, "H2": 0.75, "NH3": 1e-4}
r3 = steady_state(m3, en3, state=st3, p_gas=pg3)
print(r3.report())
print("\n   ★ 免费自检 1：元素通量守恒（线性机理，稳态下各元素进出相等）")
rr = r3.rate_net
print("      N：2·r1 = %.9e    r5 = %.9e    相对差 %.2e"
      % (2*rr[0], rr[4], abs(2*rr[0]/rr[4] - 1)))
print("      H：2·r2 = %.9e  3·r5 = %.9e    相对差 %.2e"
      % (2*rr[1], 3*rr[4], abs(2*rr[1]/(3*rr[4]) - 1)))
ok, rel = r3.flux_check(m3)
print("      内建 flux_check：最大相对不自洽 %.2e  %s" % (rel, "✓" if ok else "✗"))
print("      投影前的原始速率：", " ".join("%.9e" % v for v in r3.rate_net_raw))
print("      H 通量（投影前）：2·r2 = %.9e  3·r5 = %.9e  相对差 %.2e"
      % (2*r3.rate_net_raw[1], 3*r3.rate_net_raw[4],
         abs(2*r3.rate_net_raw[1]/(3*r3.rate_net_raw[4]) - 1)))
for w in r3.precision_warnings():
    print("      ⚠ " + w)
print("   ★ 免费自检 2：与长时瞬态交叉核对")
tr3 = transient(m3, en3, {}, (0.0, 1e4), protocol=Constant(**st3), p_gas=pg3,
                n_eval=5, method="BDF", rtol=1e-11, atol=1e-20)
for n in m3.surface_species:
    a, b = r3.theta[n], tr3.theta[n][-1]
    print("      %-6s 稳态 %.9e   ODE末 %.9e   相对差 %.2e"
          % (n, a, b, abs(b/a - 1) if a else 0.0))

print("\n══ ④ ★ 准平衡极限：与闭式解逐档对比（看在哪一档失真）══")
print("   做法：把步 2 的势垒逐步抬高 ⇒ 步 1 越来越接近平衡 ⇒ 1−z → 0")
print("   Ga2    1−z(步1)    θ 本包/闭式 相对差   r_net 本包/闭式 相对差   诊断")
for Ga2 in (1.10, 1.60, 2.10, 2.60, 3.10, 3.60):
    en4 = Energetics(m1, {1: dict(Ga_fwd=0.60, dG_rxn=-0.30),
                          2: dict(Ga_fwd=Ga2, dG_rxn=-2.0)})
    th_x, r_x, _ = closed_form(en4, T, p_HA)
    r4 = steady_state(m1, en4, state=dict(T_K=T),
                      p_gas={"HCOOH": p_HA, "CO2": 0.0, "H2": 0.0}, tol=1e-13)
    omz = 1.0 - r4.zrev[0]
    dth = abs(r4.theta["HCOOH*"]/th_x - 1)
    dr = abs(r4.rate_net[1]/r_x - 1)
    tag = "⚠已标出" if r4.precision_warnings() else "—"
    print("   %.2f   %.3e   %.3e              %.3e              %s"
          % (Ga2, omz, dth, dr, tag))
print("   ⇒ θ 与 r 一直对得上闭式解（它们不靠 1−z）；")
print("     1−z 本身是被 ln z 的绝对精度限制的那个量，所以诊断盯的是它。")
