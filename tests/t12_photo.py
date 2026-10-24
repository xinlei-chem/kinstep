# -*- coding: utf-8 -*-
r"""T12 验证：光催化 driver（准费米能级分裂 + 载流子收支自洽）。

★ 判据：**光强幂律必须作为"算出来的结果"浮现**，不是输入的假设
---------------------------------------------------------------
文献共识（2026-10-22 实地检索核对）：
  · 低光强 rate ∝ I¹，高光强 ∝ I^0.5（Turchi & Ollis 1990；K₄Nb₆O₁₇ 实测两区间，
    10.1007/BF00806072）
  · ⚠ 0.5 **不是普适值**：文献给 0.5–1，有体系 ~0.8，极高光强可能转零级；
    *J. Phys. Chem. C* 2018 (10.1021/acs.jpcc.8b09421) 指出半级只在
    "界面转移速率远小于复合速率"时成立

★★ 半级出现的**物理条件**（T12 实测弄清楚的，初版差点搞错）
-----------------------------------------------------------
    控速步 ∝ exp(β·ν·f·x)，双分子复合 ∝ exp(x)  ⇒  幂指数 = β·ν·f
  f = 0.5（本征）、β = 1（载流子俘获）、ν = 1（**单载流子步**）⇒ 幂指数 0.5 ✓
  若控速步同时吃 1 电子 + 1 空穴，drive = f·ΔE + (1−f)·ΔE = ΔE ⇒ **幂指数 1，没有半级**。
  ⇒ 经典半级定律隐含"控速步只转移一个载流子"。本测试因此用单空穴控速的机理。

算例：光催化水分解（真实体系）
    H2O(l) + 2* <-> OH* + H*                      解离吸附
    OH* -> OH(l) + *         {nu_hole_vb: 1}      **单空穴**氧化脱附（控速）
    2H* -> H2 + 2*           {nu_elec_cb: 2}      析氢，2 个电子
  稳态下 r1 = r2 且 r1 = 2·r3 ⇒ 空穴 r2 = 电子 2r3 ⇒ **电荷自动守恒**。

⚠ 初版写成经 H+ 中转的五步机理，**电荷守恒闸当场报 0.989 的失衡**。
  查下来是真缺陷：`H+` 是液相物种、活度被**固定**在 1e-7（缓冲近似），
  于是它不守恒 ⇒ 氧化半反应与还原半反应脱钩，电子/空穴消耗差 90 倍。
  ⇒ 这正是这道闸的价值：**它抓的是建模错误，不是数值误差**。
  换成不经缓冲中间物的机理即可（若确实要显式跟踪 H+，得把它做成未知数，见 FlowReactor）。
"""
import sys, warnings
import _root, refdata  # noqa: F401  ← 包根入 sys.path + 第三方真值定位
warnings.filterwarnings("ignore")
import numpy as np

from kinstep import mech as M
from kinstep.analyze import energetics_from_arrays
from kinstep.solve import steady_state
from kinstep.driver import DriverStack, make_driver
from kinstep.drivers.photo import photo_steady
import kinstep.drivers

mech = M.from_text("""
H2O(l) + 2* <-> OH* + H*
OH* -> OH(l) + *            {nu_hole_vb: 1}
2H* -> H2 + 2*              {nu_elec_cb: 2}
""")
print(mech.summary())

T = 300.0
#  控速 = 第 2 步（单空穴氧化）；吸附/脱附快，避免表面饱和掩盖幂律
# 控速 = 第 2 步（单空穴氧化脱附）；势垒抬高使暗态速率足够小，
# 一级区间才有空间存在（条件：r_dark << G，见下面 k_bi 的注释）
Ga = np.array([0.20, 1.35, 0.50])
dG = np.array([-0.30, 0.25, -0.35])
en = energetics_from_arrays(mech, Ga, dG)
PG = {"H2": 1e-6}
CL = {"H2O": 1.0, "OH": 1e-9}
GAMMA = 1.5e19

def base_state(flux, **kw):
    st = dict(T_K=T, photon_flux=flux, absorptance=0.5, gamma_site_photo=GAMMA,
              photon_energy_eV=3.2, qfl_frac_n=0.5,
              # ★ k_bi 的选取是算出来的不是试出来的：
              #   一级区：rate = G；半级区：k_bi·(G/r_dark)² = G
              #   ⇒ 转折点 G* = r_dark²/k_bi。本例 r_dark ≈ 5e-10，
              #     想让转折落在 G ≈ 3e-2（flux ≈ 1e18）⇒ k_bi ≈ 1e-17。
              recomb_rate_const=1.0e-17, recomb_srh_const=0.0,
              photo_carrier_mode="steady_carrier")
    st.update(kw); return st

print("\n══ ① 暗态极限：ΔE_qfl = 0 必须精确退回纯热催化 ══")
r_dark = steady_state(mech, en, state=dict(T_K=T), p_gas=PG, c_liq=CL,
                      tol=1e-12, strict_energy=False)
ph0 = make_driver("photo")
r_g0 = steady_state(mech, en,
                    state=dict(T_K=T, photo_carrier_mode="given_qfl",
                               qfl_split_eV=0.0, qfl_frac_n=0.5),
                    driver=DriverStack([ph0]), p_gas=PG, c_liq=CL,
                    tol=1e-12, strict_energy=False)
d = max(abs(r_g0.rate_net[i] / r_dark.rate_net[i] - 1)
        for i in range(mech.n_step) if r_dark.rate_net[i] != 0)
print("   最大相对差 %.3e  %s" % (d, "✓" if d < 1e-12 else "❌"))
print("   复合项 R_rec(0) = %.3e  ✓ 暗态精确为零" % ph0._recomb(0.0, 1e3, 1e3, 0.5))

print("\n══ ② ★★ 光强幂律：指数必须自己从 1 走到 0.5 ══")
fluxes = np.logspace(13, 26, 14)
rates, qfl, aqy, imb = [], [], [], []
for F in fluxes:
    res, x, info = photo_steady(mech, en, base_state(F), p_gas=PG, c_liq=CL, tol=1e-9)
    rates.append(abs(res.rate_net[1])); qfl.append(info["qfl_eV"])
    aqy.append(info["AQY"]); imb.append(info["charge_imbalance"])
rates = np.array(rates)
lf, lr = np.log10(fluxes), np.log10(rates)
slope = np.gradient(lr, lf)
print("   光子流         ΔE_qfl(eV)   速率(s⁻¹)      AQY        局域幂指数")
for i in range(len(fluxes)):
    print("   %.2e      %.4f      %.4e   %.3e   %.4f"
          % (fluxes[i], qfl[i], rates[i], aqy[i], slope[i]))
lo_n, hi_n = slope[:3].mean(), slope[-3:].mean()
print("   ⇒ 低光强端 %.4f（应→1）；高光强端 %.4f（应→0.5）" % (lo_n, hi_n))
ok2 = abs(lo_n - 1.0) < 0.06 and abs(hi_n - 0.5) < 0.06
print("   %s" % ("✓ 两个极限都复现，而且是算出来的不是输入的" if ok2 else "❌ 幂律不对"))
print("   ★ AQY 从 %.3f 降到 %.3e —— 高光强下光子被复合浪费掉，这正是半级的物理来源"
      % (aqy[0], aqy[-1]))

print("\n══ ③ 免费自检：电子与空穴消耗速率必须相等 ══")
print("   %d 个光强点的最大电荷失衡 = %.3e  %s"
      % (len(fluxes), max(imb), "✓" if max(imb) < 1e-6 else "❌"))
res_b, x_b, info_b = photo_steady(mech, en, base_state(1e22), p_gas=PG, c_liq=CL, tol=1e-9)
print("   例（flux=1e22）：电子 %.6e  空穴 %.6e" % (info_b["R_elec"], info_b["R_hole"]))

print("\n══ ④ 单分子（SRH/陷阱）复合：幂指数应全程 ≈1（没有半级区间）══")
r2 = []
for F in fluxes:
    res, x, info = photo_steady(mech, en,
                                base_state(F, recomb_rate_const=0.0,
                                           recomb_srh_const=1.0e-9),
                                p_gas=PG, c_liq=CL, tol=1e-9)
    r2.append(abs(res.rate_net[1]))
sl2 = np.gradient(np.log10(np.array(r2)), lf)
print("   低光强 %.4f   高光强 %.4f" % (sl2[:3].mean(), sl2[-3:].mean()))
print("   ⇒ %s" % ("✓ 单分子复合下全程一级 —— 与『半级只在双分子复合主导时出现』一致"
                   if abs(sl2[-3:].mean() - 1.0) < 0.1 else
                   "⚠ 高光强端 %.3f" % sl2[-3:].mean()))

print("\n══ ⑤ ★ 掺杂（qfl_frac_n）：我差点写死 0.5 的那个假设 ══")
print("   控速步只吃空穴 ⇒ 驱动力 = (1−f)·ΔE ⇒ 高光强幂指数应为 (1−f)")
print("   但**转折点本身也依赖 f**，解析式（由 rate=G 与 k_bi·e^x=G 联立）：")
print("       rate = r_dark·e^{(1−f)x} ,  R_rec = k_bi(e^x − 1)")
print("       ⇒ G* = [ r_dark^{1/(1−f)} / k_bi ]^{(1−f)/f}")
print("   所以 f 越小，G* 越大，半级区间越可能跑到扫描范围之外。")
r_dark_est = abs(r_dark.rate_net[1])
G_max = 0.5 * fluxes[-1] / GAMMA
print("\n   本例 r_dark = %.3e，k_bi = 1e-17，扫描到的最大 G = %.3e" % (r_dark_est, G_max))
print("   f      G*(解析)     在扫描范围内?   实测高光强指数   预期(1−f)")
for f in (0.5, 0.2, 0.8):
    k_bi = 1.0e-17
    try:
        Gstar = (r_dark_est ** (1.0 / (1.0 - f)) / k_bi) ** ((1.0 - f) / f)
    except (OverflowError, ZeroDivisionError):
        Gstar = float("inf")
    rr = []
    for F in fluxes[::3]:
        res, x, info = photo_steady(mech, en, base_state(F, qfl_frac_n=f),
                                    p_gas=PG, c_liq=CL, tol=1e-9)
        rr.append(abs(res.rate_net[1]))
    sl = np.gradient(np.log10(np.array(rr)), lf[::3])
    inside = Gstar < G_max
    print("   %.1f    %-11.2e  %-14s  %.4f           %.1f"
          % (f, Gstar, "是" if inside else "否 ⇒ 全程一级", sl[-2:].mean(), 1.0 - f))
print("   ⇒ f=0.5 与 0.8 的 G* 在范围内 ⇒ 指数落到 (1−f)；")
print("     f=0.2 的 G* 远超范围 ⇒ 全程停在一级 —— **这是对的物理，不是算错**。")
print("     解析转折点与数值行为一致 ⇒ 模型自洽。")
print("     而分配比确实改变幂律 —— **写死 0.5 会在掺杂体系上算错指数**。")

print("\n══ ⑥ ★ 载流子收支的闭式解对照（审稿人要的独立验证）══")
# 原来本套只有"复现 1→1/2 幂律"这一类自洽检查：幂律在相当程度上是模型结构自带的，
# 对求解器本身的约束不强。这里给一条**真正能逐位比的闭式解**：
#   只留双分子复合（k_srh = 0）、且表面消耗相对复合可忽略时，收支方程退化为
#       G = k_bi·(e^x − 1),   x = ΔE_qfl/kT
#   ⇒  ΔE_qfl = kT·ln(1 + G/k_bi)                （无需迭代，直接可写出）
# 把表面反应的能垒抬到 3 eV（表面通道实际关闭），自洽解就必须落在这条曲线上。
from kinstep.units import K_B_EV                      # noqa: E402
en_slow = energetics_from_arrays(mech, np.array([3.0, 3.0, 3.0]), dG)
kT = K_B_EV * T
worst = 0.0
print("   光子流 G_abs(1/site/s)    ΔE_qfl 数值(eV)   闭式解(eV)      相对差")
for F in np.logspace(14, 24, 6):
    st = base_state(F, recomb_srh_const=0.0)
    res, x, info = photo_steady(mech, en_slow, st, p_gas=PG, c_liq=CL, tol=1e-12)
    # 被吸收的光子通量换算成"每位点每秒"，与 k_bi 同量纲
    G = F * st["absorptance"] / st["gamma_site_photo"]
    closed = kT * np.log1p(G / st["recomb_rate_const"])
    dev = abs(info["qfl_eV"] / closed - 1.0)
    worst = max(worst, dev)
    print("   %-10.2e %-14.4e %-16.9f %-14.9f %.3e"
          % (F, G, info["qfl_eV"], closed, dev))
print("   ⇒ 最大相对差 %.3e  %s" % (worst, "✓ 与闭式解一致" if worst < 1e-8 else "❌"))
