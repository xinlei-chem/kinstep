# -*- coding: utf-8 -*-
"""T4 验证：三条取常数的路 + Koper 逐点对数 + 因式形式 + Wegscheider 闸。"""
import sys, math
import _root, refdata  # noqa: F401  ← 包根入 sys.path + 第三方真值定位
import numpy as np
from kinstep import mech as M
from kinstep.rates import (Energetics, rate_constants, activities, rates_factored,
                           cycle_consistency, k_tst, k_collision)
from kinstep.driver import DriverStack, make_driver
from kinstep.units import K_B, K_B_EV, H_PLANCK, Q_E
import kinstep.drivers  # 触发注册

print("══ ① TST / 碰撞理论 基准值 ══")
T = 500.0
print("  k_BT/h (500K)        = %.6e s⁻¹" % (K_B*T/H_PLANCK))
print("  k_tst(Ga=1.00eV,500K)= %.6e s⁻¹" % k_tst(1.00, T))
print("  手算核对              = %.6e" % ((K_B*T/H_PLANCK)*math.exp(-1.0/(K_B_EV*T))))
kc = k_collision(28.0, T)   # CO
print("  k_collision(CO,500K) = %.6e s⁻¹  (S0=1, A=9Å², P0=1bar)" % kc)

print("\n══ ② ★ Koper/Ager 逐点对数（T4 的硬指标）══")
# Koper 的 CO stripping：两步，都含一个电子
txt = """
H2O(l) + * <-> OH* + H+ + e-        {beta_ct: 0.5}
CO* + OH* -> CO2 + H+ + e- + 2*     {beta_ct: 0.5}
"""
mk = M.from_text(txt)
print(mk.summary())

T298 = 298.15
kT = K_B_EV * T298
alpha = 0.5
k1zero, kminus1zero, k2zero = 0.02, 1.0e4, 0.8234
en = Energetics(mk, {
    1: dict(mode="given", kf0=k1zero, kr0=kminus1zero),
    2: dict(mode="given", kf0=k2zero, dG_rxn=-1.0),   # 不可逆步，dG 只用于 Keq（不参与）
})
el = DriverStack([make_driver("electro", U_ref="SHE", gamma_site=2.2e-5)])
print("  由 kf0/kr0 反推的隐含热力学：")
rc0 = rate_constants(en, T298, el, state=dict(U_elec=0.0, T_K=T298))
for nt in rc0.notes:
    if "反推" in nt: print("   ", nt)

print("\n   U(V)    kf1(本包)     k1(notebook)   kr1(本包)     k-1(notebook)  kf2(本包)     k2(notebook)")
worst = 0.0
for U in (0.1, 0.2, 0.3, 0.4, 0.5):
    rc = rate_constants(en, T298, el, state=dict(U_elec=U, T_K=T298))
    # notebook 原式（α_1 用于全部三个常数，含那处 k2 误用 α_1 的地方）
    nb_k1  = k1zero      * math.exp( (alpha*Q_E*U)/(K_B*T298))
    nb_km1 = kminus1zero * math.exp(-((1-alpha)*Q_E*U)/(K_B*T298))
    nb_k2  = k2zero      * math.exp( (alpha*Q_E*U)/(K_B*T298))
    rel = max(abs(rc.kf[0]/nb_k1-1), abs(rc.kr[0]/nb_km1-1), abs(rc.kf[1]/nb_k2-1))
    worst = max(worst, rel)
    print("   %.1f   %.6e  %.6e   %.6e  %.6e   %.6e  %.6e"
          % (U, rc.kf[0], nb_k1, rc.kr[0], nb_km1, rc.kf[1], nb_k2))
print("   ⇒ 最大相对偏差 = %.3e（残差全部来自 notebook 的 q 是截断值）" % worst)
# 用 notebook 自己的 q 复算一遍 —— 公式若一致，应到机器精度
q_nb = 1.60217662e-19
worst2 = 0.0
for U in (0.1, 0.2, 0.3, 0.4, 0.5):
    rc = rate_constants(en, T298, el, state=dict(U_elec=U, T_K=T298))
    scale = q_nb / Q_E
    nb_k1  = k1zero      * math.exp( (alpha*q_nb*U)/(K_B*T298))
    nb_km1 = kminus1zero * math.exp(-((1-alpha)*q_nb*U)/(K_B*T298))
    nb_k2  = k2zero      * math.exp( (alpha*q_nb*U)/(K_B*T298))
    # 本包用 U'=U·(q_nb/q_codata) 等价于把 q 换成 notebook 的值
    rc2 = rate_constants(en, T298, el, state=dict(U_elec=U*scale, T_K=T298))
    worst2 = max(worst2, abs(rc2.kf[0]/nb_k1-1), abs(rc2.kr[0]/nb_km1-1), abs(rc2.kf[1]/nb_k2-1))
print("   ⇒ 换成 notebook 的 q 后最大相对偏差 = %.3e  %s"
      % (worst2, "✅ 公式逐字一致（到机器精度）" if worst2 < 1e-13 else "❌ 公式有差异"))

print("\n══ ③ 因式形式速率 + 可逆度 ══")
act = activities(mk, theta={"OH*": 0.10, "CO*": 0.75},
                 p_gas={"CO2": 1e-6}, c_liq={"H2O": 1.0, "H+": 1.0})
print("  活度向量：", dict(zip(mk.species_order, np.round(act, 8))))
rc = rate_constants(en, T298, el, state=dict(U_elec=0.3, T_K=T298))
rf, z, rn = rates_factored(mk, rc, act)
for i in range(mk.n_step):
    print("   步%d  r_fwd=%.6e  z=%.6e  r_net=%.6e" % (i+1, rf[i], z[i], rn[i]))
# 与 notebook 的 J 对照：v1 - v_minus1 应 = r_net[0]
nb_v1  = 0.02*math.exp((alpha*Q_E*0.3)/(K_B*T298))*(1-0.10-0.75)
nb_vm1 = 1e4*math.exp(-((1-alpha)*Q_E*0.3)/(K_B*T298))*0.10
print("   notebook: v1−v₋₁ = %.6e   本包 r_net[0] = %.6e   相对差 %.2e"
      % (nb_v1-nb_vm1, rn[0], abs(rn[0]/(nb_v1-nb_vm1)-1)))

print("\n══ ④ ★★ 准平衡步的灾难性相消（因式形式的存在理由）══")
# 构造一个 z→1 的步：kf 与 kr 都很大，净速率极小
txt2 = "A2 + 2* <-> 2A*"
m2 = M.from_text(txt2.replace("A2","H2").replace("A*","H*"))
en2 = Energetics(m2, {1: dict(mode="given", kf0=1e10, kr0=1e10)})
rc2 = rate_constants(en2, T298)
# θ_H 选在接近平衡处
th = 0.5
a2 = activities(m2, theta={"H*": th}, p_gas={"H2": 0.25})
rf2, z2, rn2 = rates_factored(m2, rc2, a2)
naive = rc2.kf[0]*0.25*(1-th)**2 - rc2.kr[0]*th**2
print("  kf=kr=1e10, p(H2)=0.25, θ_H=0.5")
print("   因式形式 r_fwd=%.10e  z=%.16f  r_net=%.10e" % (rf2[0], z2[0], rn2[0]))
print("   朴素相减 kf·A − kr·B = %.10e" % naive)
print("   ⇒ 两者此处相符（z 离 1 还不够近）；下面把 z 推到 1−1e-14：")
# 让 kr 稍变使 z = 1 - 1e-14
kr_t = rc2.kf[0]*0.25*(1-th)**2/(th**2)*(1-1e-14)
en3 = Energetics(m2, {1: dict(mode="given", kf0=1e10, kr0=kr_t)})
rc3 = rate_constants(en3, T298)
rf3, z3, rn3 = rates_factored(m2, rc3, a2)
naive3 = rc3.kf[0]*0.25*(1-th)**2 - rc3.kr[0]*th**2
exact3 = rc3.kf[0]*0.25*(1-th)**2*1e-14
print("   因式形式 r_net = %.10e   (1−z = %.3e)" % (rn3[0], 1-z3[0]))
print("   朴素相减 r_net = %.10e" % naive3)
print("   解析真值 r_net = %.10e" % exact3)
print("   ⇒ 因式形式误差 %.2e ；朴素相减误差 %.2e"
      % (abs(rn3[0]/exact3-1), abs(naive3/exact3-1)))

print("\n══ ⑤ Wegscheider 热力学循环闸 ══")
# 平行路径：A 可以经 B 到 C，也可以直达 C ⇒ 存在循环
txt4 = """
N2 + 2* <-> 2N*
N* + H* <-> NH* + *
NH* + H* <-> NH2* + *
H2 + 2* <-> 2H*
N* + 2H* <-> NH2* + 2*
"""
m4 = M.from_text(txt4)
# ★ 标记语义要以"测试通过与否"为准，不能以"闸返回 True/False"为准。
#   初版这里对"闸成功拦住错误"打了个 ✗，于是全量回归脚本（它按 ❌/✗ 判失败）
#   把一次**成功的拦截**误判成失败。—— 标记是给机器读的，语义必须统一。
ok, rep, _ = cycle_consistency(m4, [-0.5, -0.3, -0.2, -0.8, -0.5])
print("  自洽情形(−0.3−0.2 = −0.5)：", "✓ 判为自洽" if ok else "❌ 误报",
      rep.splitlines()[0])
ok2, rep2, _ = cycle_consistency(m4, [-0.5, -0.3, -0.2, -0.8, -0.4])
print("  填错一个数(最后一步 −0.4)：", "✓ 被闸抓住（符合预期）" if not ok2 else "❌ 漏判")
print("   ", rep2.replace("\n", "\n    "))

print("\n══ ⑥ 物理合理性提示 ══")
m5 = M.from_text("CO + * <-> CO*")
en5 = Energetics(m5, {1: dict(Ga_fwd=0.20, dG_rxn=0.80)})   # Ga < dG ⇒ 逆向负势垒
rc5 = rate_constants(en5, 500.0)
for nt in rc5.notes: print("  ", nt.replace("\n","\n   "))
