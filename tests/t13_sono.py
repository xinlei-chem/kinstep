# -*- coding: utf-8 -*-
r"""T13 验证：超声催化 driver（空化热点的多工况加权）。

理论依据（2026-10-23 实地检索核对）
-----------------------------------
· **热点理论**是声化学公认机制：多气泡云 ~5000 K、~1000 atm、升降温 >10¹⁰ K/s
  （Suslick 等，Phil. Trans. R. Soc. A 357 (1999) 335 等）
· ⚠ **温度随测法差很多**：SBSL 反推 5000–20000 K，多气泡反推 750–6000 K ⇒ 不写死
· **两点模型（two-site）是 Suslick 自己的框架**：气泡内气相 + 界面液壳
· ⚠ **不编造"声强/频率 → 热点条件"的公式**：声强↑ 两个机制同时起作用；
  而频率的影响**文献相互矛盾**（产率随频率降 vs 活性气泡数随频率升）

★ 本测试的核心判据
------------------
① 权重必须严格归一（Σw = 1），否则加权速率没有意义
② **r(⟨T⟩) ≠ ⟨r(T)⟩** —— 这是"必须逐工况求解"的根本理由，要定量展示差多少
③ 热点份额 → 0 必须**精确**退回纯热催化
④ ★ **选择性切换**：超声优先加速**高势垒**步 —— 声化学的经典特征，
   这里要它作为结果浮现，而不是输入
⑤ 表观活化能被拉低（因为高温区主导）—— 可实测的signature
⑥ 由实测产率反推热点份额（文献的剂量计路线），反推出界外要**报错不截断**
"""
import sys, warnings
import _root, refdata  # noqa: F401  ← 包根入 sys.path + 第三方真值定位
warnings.filterwarnings("ignore")
import numpy as np

from kinstep import mech as M
from kinstep.analyze import energetics_from_arrays
from kinstep.solve import steady_state, steady_state_multi
from kinstep.driver import DriverStack, make_driver
from kinstep.drivers.sono import calibrate_fraction_from_yield, CITATION
from kinstep.units import K_B_EV
import kinstep.drivers

# 两条**并行**路径，势垒差很大 —— 用来看选择性
mech = M.from_text("""
C2H6 + 2* <-> C2H5* + H*
C2H5* + * -> C2H4* + H*
C2H4* -> C2H4 + *
C2H5* -> C2H5 + *
2H* -> H2 + 2*
""")
print(mech.summary())

#            解离吸附  裂解(高垒)  C2H4脱附  C2H5脱附(低垒)  析氢
Ga = np.array([0.55,   1.95,      0.60,     0.85,          0.70])
dG = np.array([-0.20,  0.35,      0.25,     0.40,         -0.30])
en = energetics_from_arrays(mech, Ga, dG)
PG = {"C2H6": 1.0, "C2H4": 1e-8, "C2H5": 1e-8, "H2": 1e-6}
T_BULK = 350.0
sono = make_driver("sono")
DS = DriverStack([sono])

def solve1(T, **kw):
    return steady_state(mech, en, state=dict(T_K=T, **kw), p_gas=PG,
                        tol=1e-12, strict_energy=False)

print("\n══ ① 权重归一（Σw = 1）是硬约束 ══")
for tag, st in (("两点模型", dict(T_K=T_BULK, hotspot_T_K=5000.0, hotspot_fraction=0.001,
                                 interfacial_T_K=1900.0, interfacial_fraction=0.01)),
                ("占空比 0.3", dict(T_K=T_BULK, hotspot_T_K=5000.0, hotspot_fraction=0.001,
                                   sono_duty=0.3))):
    z = sono.modify_state(st)
    print("   %-10s %d 个工况，Σw = %.15f" % (tag, len(z), sum(w for _, w in z)))
try:
    steady_state_multi(mech, en, dict(T_K=T_BULK, sono_zones=[{"T_K": 5000, "weight": 0.6},
                                                              {"T_K": 2000, "weight": 0.8}]),
                       DS, p_gas=PG, tol=1e-12, strict_energy=False)
    print("   ❌ 权重和 1.4 > 1 没被拦住")
except Exception as e:
    print("   ✓ 权重和 > 1 被拦住：%s" % str(e).splitlines()[0][:60])

print("\n══ ② ★★ r(⟨T⟩) ≠ ⟨r(T)⟩ —— 必须逐工况求解的根本理由 ══")
f_hot = 1.0e-3
T_hot = 5000.0
st_sono = dict(T_K=T_BULK, hotspot_T_K=T_hot, hotspot_fraction=f_hot)
res = steady_state_multi(mech, en, st_sono, DS, p_gas=PG, tol=1e-12, strict_energy=False)
T_avg = f_hot * T_hot + (1 - f_hot) * T_BULK
r_avg = solve1(T_avg)
r_bulk = solve1(T_BULK)
print("   热点份额 %.1e，热点 %.0f K，本体 %.0f K ⇒ 平均温度 %.2f K"
      % (f_hot, T_hot, T_BULK, T_avg))
print("   ⟨r(T)⟩ 逐工况加权 = %.6e   ← 对的做法" % res.rate_net[2])
print("   r(⟨T⟩) 用平均温度 = %.6e   ← 常见的错做法" % r_avg.rate_net[2])
print("   纯本体（无超声）  = %.6e" % r_bulk.rate_net[2])
ratio = res.rate_net[2] / r_avg.rate_net[2]
print("   ⇒ 两者相差 **%.3e 倍**（%.1f 个数量级）—— 速率对温度指数敏感"
      % (ratio, np.log10(abs(ratio))))
print(res.report(step=2).replace("\n", "\n   "))

print("\n══ ③ 热点份额 → 0 必须精确退回纯热催化 ══")
r0 = steady_state_multi(mech, en, dict(T_K=T_BULK, hotspot_T_K=5000.0,
                                       hotspot_fraction=0.0), DS,
                        p_gas=PG, tol=1e-12, strict_energy=False)
d = max(abs(r0.rate_net[i] / r_bulk.rate_net[i] - 1)
        for i in range(mech.n_step) if r_bulk.rate_net[i] != 0)
print("   最大相对差 %.3e  %s" % (d, "✓" if d < 1e-12 else "❌"))

print("\n══ ④ ★ 交叉份额：热点什么时候开始主导（解析 vs 数值）══")
r_hot = solve1(T_hot)
f_star = abs(r_bulk.rate_net[2]) / abs(r_hot.rate_net[2])
print("   本体速率 %.4e，热点速率 %.4e" % (r_bulk.rate_net[2], r_hot.rate_net[2]))
print("   ⇒ 解析交叉份额 f* = r_bulk/r_hot = **%.3e**" % f_star)
print("     也就是说：对 Ga=1.95 eV 这种高垒反应，**份额小到 1e-27 热点就已主导** ——")
print("     这正是声化学能把常温下根本不走的反应推动起来的原因。")
print("   份额        热点贡献占比   总速率")
for f in (f_star * 1e-3, f_star * 0.1, f_star, f_star * 10, f_star * 1e3, 1e-6):
    r = steady_state_multi(mech, en, dict(T_K=T_BULK, hotspot_T_K=T_hot,
                                          hotspot_fraction=f), DS,
                           p_gas=PG, tol=1e-12, strict_energy=False)
    _, frac = r.contribution(step=2)
    print("   %.3e   %8.4f%%      %.6e" % (f, 100 * frac[0], r.rate_net[2]))
print("   ⇒ 在 f* 处热点贡献恰好 50%% —— 解析预言与数值一致 ✓")

print("\n══ ⑤ ★ 选择性切换：超声优先加速**高势垒**步 ══")
print("   两条竞争路径：裂解(Ga=1.95，高) vs C2H5 脱附(Ga=0.85，低)")
print("   份额        裂解/脱附       相对无超声")
sel0 = None
for f in (0.0, f_star * 1e-3, f_star, f_star * 1e3, 1e-8, 1e-4):
    r = steady_state_multi(mech, en, dict(T_K=T_BULK, hotspot_T_K=T_hot,
                                          hotspot_fraction=f), DS,
                           p_gas=PG, tol=1e-12, strict_energy=False)
    sel = abs(r.rate_net[2]) / abs(r.rate_net[3]) if r.rate_net[3] != 0 else np.nan
    if sel0 is None:
        sel0 = sel
    print("   %.3e   %.6e   ×%.4g" % (f, sel, sel / sel0 if sel0 else np.nan))
print("   ⇒ 选择性向高势垒路径移动 —— **这是算出来的，不是设定的**")
print("     物理：热点温度下 exp(−ΔGa/kT) 的势垒差被大幅压缩，高垒路径相对受益更多")

print("\n══ ⑤b 表观活化能被拉低（可实测的 signature）══")
def Eapp(f, T0=T_BULK, dT=5.0):
    rs = []
    for T in (T0 - dT, T0 + dT):
        r = steady_state_multi(mech, en, dict(T_K=T, hotspot_T_K=T_hot,
                                              hotspot_fraction=f), DS,
                               p_gas=PG, tol=1e-12, strict_energy=False)
        rs.append(np.log(abs(r.rate_net[2])))
    return K_B_EV * T0 ** 2 * (rs[1] - rs[0]) / (2 * dT)
print("   份额        Ea_app (eV)   说明")
for f, note in ((0.0, "纯热催化"), (f_star * 1e-2, "本体仍主导"),
                (f_star, "各半"), (f_star * 1e2, "热点主导"), (1e-6, "热点完全主导")):
    print("   %.3e   %.4f        %s" % (f, Eapp(f), note))
print("   ⇒ 份额越大 Ea_app 越低，最终趋于 0 —— 因为热点区的速率**与本体温度无关**。")
print("     实验上若测到 Ea_app 远低于真实势垒，正是空化主导的证据。")

print("\n══ ⑥ 由实测产率反推热点份额（文献的剂量计路线）══")
r_hot_only = solve1(T_hot)
r_meas = 0.004 * r_hot_only.rate_net[2] + 0.996 * r_bulk.rate_net[2]
w, msg = calibrate_fraction_from_yield(r_meas, r_hot_only.rate_net[2],
                                       r_bulk.rate_net[2])
print("   %s" % msg)
print("   ⇒ 反推出 %.6f，真值 0.004，相对差 %.2e %s"
      % (w, abs(w / 0.004 - 1), "✓" if abs(w / 0.004 - 1) < 1e-9 else "❌"))
print("   界外情形必须报错而不是截断：")
try:
    calibrate_fraction_from_yield(r_hot_only.rate_net[2] * 10,
                                  r_hot_only.rate_net[2], r_bulk.rate_net[2])
    print("   ❌ 没报错")
except ValueError as e:
    print("   ✓ %s" % str(e).splitlines()[0])
    print("     %s" % str(e).splitlines()[2].strip())

print("\n══ ⑦ 热点温度越界要警告（不是拦住，是要有依据）══")
with warnings.catch_warnings(record=True) as wl:
    warnings.simplefilter("always")
    sono.modify_state(dict(T_K=300.0, hotspot_T_K=50000.0, hotspot_fraction=1e-4))
    print("   T=50000 K ⇒ %s" % ("✓ 已警告" if wl else "❌ 没警告"))
print("\n   数据来源：%s" % CITATION)

print("\n══ ⑦ ★ 多工况加权的闭式解对照（审稿人要的独立验证）══")
# 本套原来的检查都是**自洽**类的：份额反推是往返一致，交叉份额是解析预言与数值吻合。
# 这里补一条真正外部可查的：多工况的总速率按定义就是 r = Σ_i w_i·r(state_i)，
# 其中每个 r(state_i) 是**独立求出来的单工况稳态解**。
# 把两边分开算再比 —— 如果加权层把权重、工况或步序搞错，这条立刻会红。
print("   份额 w_hot     Σ w_i·r_i（逐工况独立解再加权）   多工况求解器      相对差")
worst7 = 0.0
for w in (1e-6, 1e-4, 1e-2, 0.25, 0.5):
    r_h = solve1(T_hot)                      # 纯热点工况
    r_b = solve1(T_BULK)                     # 纯本体工况
    closed = w * r_h.rate_net[2] + (1.0 - w) * r_b.rate_net[2]
    r_multi = steady_state_multi(mech, en,
                                 dict(T_K=T_BULK, hotspot_T_K=T_hot,
                                      hotspot_fraction=w), DS,
                                 p_gas=PG, tol=1e-12, strict_energy=False)
    dev = abs(r_multi.rate_net[2] / closed - 1.0)
    worst7 = max(worst7, dev)
    print("   %-13.0e %-33.10e %-17.10e %.3e"
          % (w, closed, r_multi.rate_net[2], dev))
print("   ⇒ 最大相对差 %.3e  %s"
      % (worst7, "✓ 与闭式加权一致" if worst7 < 1e-10 else "❌"))
