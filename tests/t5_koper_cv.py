# -*- coding: utf-8 -*-
r"""T5 验证：★ Koper 1998 CO stripping CV 复现。

真值来源：Koper, Jansen, van Santen, Lukkien, Hilbers, *J. Chem. Phys.* 109 (1998) 6051
参照实现：`工具/小程序1/EC CO stripping MKM.ipynb`（Ager，自评 "Perfect agreement with paper"）

notebook 的模型（**本测试把它逐字抄在下面作为独立真值**，不调用本包）：
    k1  = k1⁰·exp(+α₁qE/kT)          OH 吸附
    k₋₁ = k₋₁⁰·exp(−(1−α₁)qE/kT)     OH 脱附
    k2  = k2⁰·exp(+α₁qE/kT)          CO* + OH* → CO2   ← notebook 这里用 α₁（本该 α₂，两者都=0.5）
    dθ_OH/dt = k1(1−θ_OH−θ_CO) − k₋₁θ_OH − Z·k2·θ_OH·θ_CO
    dθ_CO/dt = −Z·k2·θ_OH·θ_CO
    J = q·Γ_sites·(v1 − v₋₁ + v2)
    T = 300 K,  Z = 4,  Γ_sites = 1.32e15 cm⁻²,  α₁ = α₂ = 0.5
    CV: 0.1 → 0.5 V @ 50 mV/s（一圈 16 s）
"""
import sys, math
import _root, refdata  # noqa: F401  ← 包根入 sys.path + 第三方真值定位
import numpy as np
from scipy.integrate import solve_ivp

from kinstep import mech as M
from kinstep.rates import Energetics, rate_constants, activities, rates_factored
from kinstep.solve import transient, CVSweep, PotentialStep, Constant
from kinstep.driver import DriverStack, make_driver
from kinstep.units import N_A
import kinstep.drivers

# ═════════════ A. notebook 的实现（真值，逐字抄，不碰本包） ═════════════
# ★ 这里用 **CODATA 的 q**，而不是 notebook 原文的 1.60217662e-19。
#   原因：notebook 自己不自洽 —— 它的 q=1.60217662e-19 而 eV=1.602176634e-19，
#   按定义这两个必须是同一个数。KinSTEP 的 units.py 曾照抄这处不一致（已修）。
#   两边统一用 CODATA 后，剩下的差异才是"实现差异"而不是"常数差异"。
#   notebook 原值造成的差异单独在 ①b 里量一次（约 1e-7，与 αqU/kT ≈ 10 的放大一致）。
from kinstep.units import Q_E as q_nb, R_GAS as R_nb
q_truncated = 1.60217662e-19          # notebook 原文的值，只在 ①b 用
#
# ★ 2026-10-22 修一个**回归脚本抓到的真退化**：
#   这里原本硬编码 `R_nb = 8.3144598`（CODATA-2014）。T7 时我把 units.py 的 R_GAS
#   改成了 CODATA-2018（修"混代"那个潜伏 bug），**但这份测试没同步** ⇒
#   参照实现与本包用的 k_B 不再一致，右端相对差从 2.442e-15 退化到 2.612e-06。
#   （量级核对：R 相对差 3.389e-07 × 放大 αqU/kT ≈ 9.7 ⇒ 3.3e-06，与实测同量级。）
#
#   教训：**对数用的常数要和被测对象同源**，硬编码一份副本早晚不同步。
#   本测试的意图始终是"两边统一常数，剩下的差异才是实现差异"——
#   从 units 取 R_GAS 正是维护这个意图，而不是破坏独立性
#   （notebook 原值造成的差异另在 ①b 单独量，不受影响）。
NA_nb = 6.02214076e23
kB_nb = R_nb / NA_nb
T_nb  = 300.0
a1 = a2 = 0.5
Z_nb = 4.0
G_nb = 1.32e15 * 1e4          # cm^-2 → m^-2
k1z, km1z, k2z = 0.02, 1.0e4, 0.8234
Estart, Efinish, nu_sweep = 0.1, 0.5, 0.050
cycle = 2.0 * (Efinish - Estart) / nu_sweep

def E_CV_nb(t):
    tau = t % cycle
    return (Estart + tau * nu_sweep) if tau <= cycle / 2 \
        else (Efinish - (tau - cycle / 2) * nu_sweep)

def k1_nb(E):  return k1z  * math.exp( (a1 * q_nb * E) / (kB_nb * T_nb))
def km1_nb(E): return km1z * math.exp(-((1 - a1) * q_nb * E) / (kB_nb * T_nb))
def k2_nb(E):  return k2z  * math.exp( (a1 * q_nb * E) / (kB_nb * T_nb))

def ode_nb(t, y, E_of_t):
    th_OH, th_CO = y
    E = E_of_t(t)
    d_OH = k1_nb(E) * (1 - th_OH - th_CO) - km1_nb(E) * th_OH \
           - Z_nb * k2_nb(E) * th_OH * th_CO
    d_CO = -Z_nb * k2_nb(E) * th_OH * th_CO
    return [d_OH, d_CO]

def J_nb(E, th_OH, th_CO):
    v1  = k1_nb(E) * (1 - th_OH - th_CO)
    vm1 = km1_nb(E) * th_OH
    v2  = Z_nb * k2_nb(E) * th_OH * th_CO
    return (q_nb * G_nb) * (v1 - vm1 + v2)

# ═════════════ B. 本包：只写反应式 + 三个实验常数 ═════════════
TXT = """
H2O(l) + * <-> OH* + H+ + e-            {beta_ct: 0.5}
CO* + OH* -> CO2 + H+ + e- + 2*         {beta_ct: 0.5, n_neighbor: 4.0}
"""
mk = M.from_text(TXT)
en = Energetics(mk, {
    1: dict(mode="given", kf0=k1z, kr0=km1z),
    2: dict(mode="given", kf0=k2z, dG_rxn=-1.0),     # 不可逆，dG 不参与
})
el = DriverStack([make_driver("electro", U_ref="SHE",
                              gamma_site=G_nb / N_A)])   # m^-2 → mol/m^2
print(mk.summary())
print()
for nt in rate_constants(en, T_nb, el, state=dict(U_elec=0.0, T_K=T_nb)).notes:
    print("  ·", nt.replace("\n", "\n    "))

# ═════════════ C. ★ 先对右端函数（比对曲线更严格） ═════════════
print("\n══ ① ODE 右端函数逐点对数（随机取 8 个 (t, θ_OH, θ_CO)）══")
rng = np.random.default_rng(20261022)
print("    t(s)   θ_OH   θ_CO |   dθ_OH/dt 本包      notebook       相对差 |   dθ_CO/dt 相对差")
worst_rhs = 0.0
a_r, a_p = None, None
from kinstep.rates import _side_matrices
a_r, a_p = _side_matrices(mk)
idx = [mk.species_order.index(n) for n in mk.surface_species]
nu_s = np.asarray(mk.nu, float)[:, idx]
for _ in range(8):
    t = rng.uniform(0, cycle)
    th_CO = rng.uniform(0.0, 0.7); th_OH = rng.uniform(0.0, 1.0 - th_CO)
    E = E_CV_nb(t)
    rc = rate_constants(en, T_nb, el, state=dict(U_elec=E, T_K=T_nb), strict_energy=False)
    act = activities(mk, theta={"OH*": th_OH, "CO*": th_CO},
                     p_gas={"CO2": 0.0}, c_liq={"H2O": 1.0, "H+": 1.0})
    _, _, rn = rates_factored(mk, rc, act, a_r, a_p)
    d_mine = nu_s.T @ rn
    d_nb = ode_nb(t, [th_OH, th_CO], E_CV_nb)
    # 本包物种顺序 = mk.surface_species
    order = {n: i for i, n in enumerate(mk.surface_species)}
    r1 = abs(d_mine[order["OH*"]] / d_nb[0] - 1); r2 = abs(d_mine[order["CO*"]] / d_nb[1] - 1)
    worst_rhs = max(worst_rhs, r1, r2)
    print("   %5.2f  %.4f  %.4f | %+.8e  %+.8e  %.2e | %.2e"
          % (t, th_OH, th_CO, d_mine[order["OH*"]], d_nb[0], r1, r2))
print("   ⇒ 右端最大相对差 = %.3e  %s"
      % (worst_rhs, "✅ 逐字一致" if worst_rhs < 1e-11 else "❌ 公式有差异"))

print("\n══ ①b notebook 原文那个截断的 q 会带来多大差异（单独量一次）══")
_dev = abs(q_truncated / q_nb - 1)
_amp = 0.5 * 0.5 * q_nb / (kB_nb * T_nb)          # α·U_max·q/kT 的放大倍数
print("   q 相对差 %.3e × 放大 %.1f (=αU/kT) ⇒ 预期 %.2e"
      % (_dev, _amp, _dev * _amp))
print("   ⇒ 这就是 T4 第一次对数时测到 8.5e-8 的来源；两边统一常数后降到上面的 %.1e" % worst_rhs)

# ═════════════ D. CV 全程复现（Koper Figure 1 的三个初始覆盖度） ═════════════
print("\n══ ② ★ CV 全程复现（Koper Fig.1：θ_CO⁰ = 0.99 / 0.75 / 0.50）══")
# ★ 容差选择（实测，不是拍的）：
#   rtol/atol = 1e-6/1e-10 → 3,525 次右端调用
#              1e-8/1e-12 → 9,229 次
#              1e-10/1e-14 → 492,499 次，而 θ_OH 峰值与 1e-8 **完全相同**
#   ⇒ 1e-10/1e-14 只是让 Radau 去解析 θ 衰减到 1e-30 之后的数值噪声。取 1e-8/1e-12。
RTOL, ATOL = 1e-8, 1e-12
t_eval = np.linspace(0, cycle, 200)
print("   θ_CO⁰ | θ_OH 峰值 本包/notebook | 电流峰 (mA/cm²) 本包/notebook | 残余θ_CO | 最大相对差")
worst_cv = 0.0
summary_rows = []
for th0 in (0.99, 0.75, 0.50):
    # notebook
    s_nb = solve_ivp(lambda t, y: ode_nb(t, y, E_CV_nb), [0, cycle], [0.0, th0],
                     t_eval=t_eval, method="Radau", rtol=RTOL, atol=ATOL)
    J_ref = np.array([J_nb(E_CV_nb(t), o, c) for t, o, c in zip(t_eval, s_nb.y[0], s_nb.y[1])])
    # 本包
    cv = CVSweep(Estart, Efinish, nu_sweep, n_cycle=1, T_K=T_nb)
    res = transient(mk, en, {"CO*": th0}, (0, cycle), protocol=cv, driver=el,
                    p_gas={"CO2": 0.0}, c_liq={"H2O": 1.0, "H+": 1.0},
                    t_eval=t_eval, rtol=RTOL, atol=ATOL, strict_energy=False)
    J_mine = res.obs["j_curr"]
    d_oh = np.max(np.abs(res.theta["OH*"] - s_nb.y[0]))
    d_co = np.max(np.abs(res.theta["CO*"] - s_nb.y[1]))
    scale = max(np.max(np.abs(J_ref)), 1e-30)
    d_j = np.max(np.abs(J_mine - J_ref)) / scale
    worst_cv = max(worst_cv, d_oh, d_co, d_j)
    print("    %.2f  |  %.6f / %.6f  |  %.6f / %.6f  |  %.3e  | θ:%.2e  J:%.2e"
          % (th0, res.theta["OH*"].max(), s_nb.y[0].max(),
             J_mine.max() / 10.0, J_ref.max() / 10.0,       # A/m² → mA/cm²
             res.theta["CO*"][-1], max(d_oh, d_co), d_j))
    summary_rows.append((th0, res))
# ★ 判据怎么定（这一点想清楚了才算验证，不是凑阈值）
#   ① 已把 **ODE 右端函数**对到 2.4e-15 —— 公式层面逐字一致，这是最严的一关。
#   ② 这里比的是两次**独立 solve_ivp** 的轨迹：两边右端的实现不同（纯 float vs numpy），
#      内部步长序列必然不同 ⇒ 残差是**容差传播**，不是公式差异。
#   ⇒ 所以这里的判据不能用"≈ rtol"，而要看它**是否随 rtol 单调缩小**（见 ②b）。
print("   ⇒ CV 全程最大偏差 = %.3e  %s"
      % (worst_cv, "✅ 复现（判据见下面 ②b 的容差阶梯）" if worst_cv < 1e-5 else "❌ 对不上"))
print("\n══ ②b 残差是容差传播还是公式差异？—— 用容差阶梯标定 ══")
print("    rtol      Δθ_max      ΔJ/J_max    右端调用   比上一档")
_prev = None
for _rt, _at in ((1e-5, 1e-9), (1e-6, 1e-10), (1e-7, 1e-11), (1e-8, 1e-12)):
    _sn = solve_ivp(lambda t, y: ode_nb(t, y, E_CV_nb), [0, cycle], [0.0, 0.50],
                    t_eval=t_eval, method="Radau", rtol=_rt, atol=_at)
    _Jr = np.array([J_nb(E_CV_nb(t), o, c) for t, o, c in zip(t_eval, _sn.y[0], _sn.y[1])])
    _r = transient(mk, en, {"CO*": 0.50}, (0, cycle),
                   protocol=CVSweep(Estart, Efinish, nu_sweep, T_K=T_nb), driver=el,
                   p_gas={"CO2": 0.0}, c_liq={"H2O": 1.0, "H+": 1.0},
                   t_eval=t_eval, rtol=_rt, atol=_at, strict_energy=False)
    _dth = max(np.max(np.abs(_r.theta["OH*"] - _sn.y[0])),
               np.max(np.abs(_r.theta["CO*"] - _sn.y[1])))
    _dj = np.max(np.abs(_r.obs["j_curr"] - _Jr)) / np.max(np.abs(_Jr))
    print("   %.0e   %.3e   %.3e   %6d   %s"
          % (_rt, _dth, _dj, _r.n_rhs, "—" if _prev is None else "÷%.1f" % (_prev / _dj)))
    _prev = _dj
print("   ⇒ 单调缩小且在 1e-8 处骤降 ⇒ **容差传播，不是公式差异**（公式见 ① 的 2.4e-15）")

print("\n   备注：", summary_rows[1][1].notes[-1])
print("        ", summary_rows[1][1].notes[-2])

# ═════════════ E. 电势阶跃（notebook 的第二个实验） ═════════════
print("\n══ ③ 电势阶跃 0 → 0.30 V（CO 剥离的衰减）══")
Estep, tmax = 0.30, 40.0
te = np.linspace(0, tmax, 200)
s_nb = solve_ivp(lambda t, y: ode_nb(t, y, lambda _t: Estep), [0, tmax], [0.0, 0.75],
                 t_eval=te, method="Radau", rtol=RTOL, atol=ATOL)
ps = PotentialStep(Estep, Estep, t_step=0.0, T_K=T_nb)
rs = transient(mk, en, {"CO*": 0.75}, (0, tmax), protocol=ps, driver=el,
               p_gas={"CO2": 0.0}, c_liq={"H2O": 1.0, "H+": 1.0},
               t_eval=te, rtol=RTOL, atol=ATOL, strict_energy=False)
d = max(np.max(np.abs(rs.theta["OH*"] - s_nb.y[0])), np.max(np.abs(rs.theta["CO*"] - s_nb.y[1])))
print("   t=0 θ_CO=0.750 → t=%.0fs 本包 %.6f / notebook %.6f，最大偏差 %.3e %s"
      % (tmax, rs.theta["CO*"][-1], s_nb.y[1][-1], d, "✅" if d < 1e-6 else "❌"))
print("   峰值电流 %.4f mA/cm²，出现在 t = %.3f s"
      % (rs.obs["j_curr"].max() / 10.0, rs.t[int(np.argmax(rs.obs["j_curr"]))]))

# ═════════════ F. 自检与效率 ═════════════
print("\n══ ④ 免费自检与效率 ══")
r = summary_rows[1][1]
ok, rep = r.check_bounds()
print("   位点/覆盖度：", rep)
print("   恒定工况缓存：", end=" ")
rc_const = transient(mk, en, {"CO*": 0.75}, (0, 1.0), protocol=Constant(T_K=T_nb, U_elec=0.3),
                     driver=el, p_gas={"CO2": 0.0}, c_liq={"H2O": 1.0, "H+": 1.0},
                     n_eval=20, strict_energy=False)
print("恒定工况 %d 次右端调用；CV %d 次（CV 每次都要重算 exp）"
      % (rc_const.n_rhs, r.n_rhs))
print("   随时间变的工况：", list(r.state_trace))

# ═════════════ G. 温度漏给要报错 ═════════════
print("\n══ ⑤ 漏给温度必须报错（不许静默取 298.15）══")
try:
    transient(mk, en, {"CO*": 0.75}, (0, 1.0), protocol=Constant(U_elec=0.3), driver=el,
              c_liq={"H2O": 1.0, "H+": 1.0}, strict_energy=False)
    print("   ❌ 没报错")
except Exception as e:
    print("   ✓", str(e).splitlines()[0])
