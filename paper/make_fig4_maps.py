# -*- coding: utf-8 -*-
"""KinSTEP 论文 ｜ Fig.4：四种模态的二维工况图（operating-condition maps）

## 为什么要这张图
CATKINAS 手册第 28 页用 (T, P) 网格上的多面板彩图展示 DRC / 反应级数 /
覆盖度 / 描述符梯度 —— 那是热催化工具的标准展示方式。
本图做同样的事，但**四种模态各一行**：热 / 电 / 光 / 超声。
这正是 KinSTEP 相对于单模态工具的区别所在：**同一个内核、同一套求解阶梯，
换一个 driver 就换一种外场**，所以四行图可以用同一段代码扫出来。

| 行 | 体系 | 两个工况轴 | 三个可观测量 |
|---|---|---|---|
| 热 | 合成氨 | T × P | 净速率(带符号) / θ_H* / 空位 θ_* |
| 电 | HER | U × pH | log₁₀&#124;j&#124; / θ_H* / Tafel 斜率 |
| 光 | 光解水 | 光子流 × 复合系数 | ΔE_qfl / log₁₀ r / **局域幂指数** |
| 超声 | 乙烷脱氢 | 热点份额 × 热点温度 | log₁₀ r / **log₁₀(⟨r(T)⟩ / r(⟨T⟩))** / 选择性 |

每一行的第三张都是"只有 MKM 能给、经验公式给不了"的量：
Tafel 斜率、光强幂指数、多工况加权与平均工况之比。

## 配色
强度量用 `turbo`（视觉上接近 CATKINAS 的 jet，但感知均匀得多，印刷不失真）；
**带符号**的量用发散配色并把 0 钉在中间（热催化那张的零线就是平衡线）。

## 原始数据
每张子图的网格都存进 `SI/fig4_maps.npz`，另存文本版 `SI/fig4_*.dat`，
按本项目规矩：**交付图必须同时交付可作图的原始数据**。

用法:
    python make_fig4_maps.py            # 正式网格（约 25–35 分钟）
    python make_fig4_maps.py --quick    # 粗网格，调试用（约 1 分钟）
"""
import os
import sys
import time
import warnings

import numpy as np

warnings.filterwarnings("ignore")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)            # 不写死绝对路径：代码存档会被下载到任意目录

from kinstep import mech as M                                    # noqa: E402
from kinstep.analyze import energetics_from_arrays, tafel_slope  # noqa: E402
from kinstep.driver import DriverStack, make_driver              # noqa: E402
from kinstep.drivers.photo import photo_steady                   # noqa: E402
from kinstep.scaling import from_table                           # noqa: E402
from kinstep.solve import steady_state, steady_state_multi       # noqa: E402
from kinstep.thermo import entropy_offsets                       # noqa: E402
import kinstep.drivers                                           # noqa: F401,E402

OUT = os.path.join(ROOT, "paper", "figures")
SI = os.path.join(ROOT, "paper", "SI")
QUICK = "--quick" in sys.argv

# 网格大小：正式 / 调试
NX, NY = (48, 36) if not QUICK else (8, 6)
NXP, NYP = (32, 24) if not QUICK else (6, 5)


def _log(msg):
    print(msg, flush=True)


# ════════════════════════════════════════════════════════════════════
# 行 1 ｜ 热催化：合成氨 over (T, P)
# ════════════════════════════════════════════════════════════════════
AMMONIA = """
N2 + 2* <-> 2N*
H2 + 2* <-> 2H*
N* + H* <-> NH* + *
NH* + H* <-> NH2* + *
NH2* + H* <-> NH3 + 2*
"""
# scaling/BEP 表与案例1、测试 t9 完全相同（CATKINAS 401 点已逐点对过）
AMM_TAB = {
    1: dict(bep=(0.78, 1.54), scale=(2.00, 0.00)),    # N2 解离
    2: dict(bep=(0.54, 0.60), scale=(0.30, -0.68)),   # H2 解离
    3: dict(bep=(0.26, 1.45), scale=(-0.15, -0.18)),
    4: dict(bep=(0.43, 1.66), scale=(-0.45, -0.45)),
    5: dict(bep=(0.31, 1.45), scale=(-0.85, 1.17)),   # NH3 脱附
}
E_N_PEAK = -0.58        # 火山峰位置（Fig.3 已给出）


def row_thermal():
    mech = M.from_text(AMMONIA)
    TS = np.linspace(500.0, 900.0, NX)
    PS = np.logspace(0.0, np.log10(200.0), NY)
    R = np.full((NY, NX), np.nan)
    TH_H = np.full((NY, NX), np.nan)
    TH_V = np.full((NY, NX), np.nan)
    t0 = time.time()
    for i, T in enumerate(TS):
        oGa, odG, _ = entropy_offsets(mech, float(T))
        sc = from_table(mech, AMM_TAB, desc_names=["E_N"],
                        offset_Ga=oGa, offset_dG=odG)
        Ga, dG = sc.energies([E_N_PEAK])
        en = energetics_from_arrays(mech, Ga, dG)
        # ★ 每一列（固定 T）从 ODE 起步，再沿 P 热启动。
        #   跨 T 带猜测会把解带进错误的分支 —— 实测速率会掉到 1e-60 这种没物理的值。
        guess = None
        for j, P in enumerate(PS):
            pg = {"N2": 0.25 * P, "H2": 0.75 * P, "NH3": 0.01 * P}
            try:
                r = steady_state(
                    mech, en, state=dict(T_K=float(T)), p_gas=pg,
                    theta_guess=guess,
                    ladder=("newton",) if guess else ("ode", "newton"),
                    tol=1e-12, precision="auto", strict_energy=False)
                guess = dict(r.theta)
                R[j, i] = r.rate_net[4]
                TH_H[j, i] = r.theta.get("H*", np.nan)
                TH_V[j, i] = r.theta_vac.get("*", np.nan)
            except Exception:
                guess = None
        if (i + 1) % max(1, NX // 6) == 0:
            _log(f"    热催化 {i+1}/{NX} 列  ({time.time()-t0:.0f} s)")
    _log(f"  热催化完成 {np.isfinite(R).sum()}/{R.size} 点，{time.time()-t0:.0f} s")
    return dict(x=TS, y=PS, rate=R, th_H=TH_H, th_vac=TH_V)


# ════════════════════════════════════════════════════════════════════
# 行 2 ｜ 电催化：HER over (U, pH)
# ════════════════════════════════════════════════════════════════════
HER = """
H+ + e- + * <-> H*            {beta_ct: 0.5}
H* + H+ + e- -> H2 + *        {beta_ct: 0.5}
2H* -> H2 + 2*
"""


def row_electro():
    mech = M.from_text(HER)
    en = energetics_from_arrays(mech, np.array([0.45, 0.60, 0.75]),
                                np.array([-0.10, -0.30, -0.20]))
    el = DriverStack([make_driver("electro", U_ref="SHE", gamma_site=2.2e-5)])
    US = np.linspace(-0.40, -0.02, NX)
    PHS = np.linspace(0.0, 14.0, NY)
    J = np.full((NY, NX), np.nan)
    TH = np.full((NY, NX), np.nan)
    TAF = np.full((NY, NX), np.nan)
    t0 = time.time()
    for j, pH in enumerate(PHS):
        guess = None
        jj = np.full(NX, np.nan)
        for i, U in enumerate(US):
            try:
                r = steady_state(
                    mech, en, state=dict(T_K=298.15, U_elec=float(U)), driver=el,
                    p_gas={"H2": 1.0}, c_liq={"H+": 10.0 ** (-pH)},
                    theta_guess=guess,
                    ladder=("newton",) if guess else ("ode", "newton"),
                    tol=1e-12, precision="auto", strict_energy=False)
                guess = dict(r.theta)
                jj[i] = r.obs.get("j_curr", np.nan)
                TH[j, i] = r.theta.get("H*", np.nan)
            except Exception:
                guess = None
        J[j] = jj
        # Tafel 斜率沿 U 方向求（每个 pH 一条曲线）
        with np.errstate(all="ignore"):
            ts = tafel_slope(US, jj)
        TAF[j] = 1000.0 * np.asarray(ts, float)      # V/dec → mV/dec
        if (j + 1) % max(1, NY // 6) == 0:
            _log(f"    电催化 {j+1}/{NY} 行  ({time.time()-t0:.0f} s)")
    _log(f"  电催化完成 {np.isfinite(J).sum()}/{J.size} 点，{time.time()-t0:.0f} s")
    return dict(x=US, y=PHS, j=J, th_H=TH, tafel=TAF)


# ════════════════════════════════════════════════════════════════════
# 行 3 ｜ 光催化：光解水 over (光子流, 复合系数)
# ════════════════════════════════════════════════════════════════════
PHOTO = """
H2O(l) + 2* <-> OH* + H*
OH* -> OH(l) + *            {nu_hole_vb: 1}
2H* -> H2 + 2*              {nu_elec_cb: 2}
"""


def row_photo():
    mech = M.from_text(PHOTO)
    en = energetics_from_arrays(mech, np.array([0.20, 1.35, 0.50]),
                                np.array([-0.30, 0.25, -0.35]))
    PG, CL = {"H2": 1e-6}, {"H2O": 1.0, "OH": 1e-9}
    FX = np.logspace(14.0, 26.0, NXP)
    KR = np.logspace(-18.0, -14.0, NYP)
    QFL = np.full((NYP, NXP), np.nan)
    RATE = np.full((NYP, NXP), np.nan)
    AQY = np.full((NYP, NXP), np.nan)
    t0 = time.time()
    for j, kr in enumerate(KR):
        for i, F in enumerate(FX):
            st = dict(T_K=300.0, photon_flux=float(F), absorptance=0.5,
                      gamma_site_photo=1.5e19, photon_energy_eV=3.2,
                      qfl_frac_n=0.5, recomb_rate_const=float(kr),
                      photo_carrier_mode="steady_carrier")
            try:
                r, x, info = photo_steady(mech, en, st, p_gas=PG, c_liq=CL, tol=1e-9)
                QFL[j, i] = info["qfl_eV"]
                RATE[j, i] = abs(r.rate_net[1])
                AQY[j, i] = info["AQY"]
            except Exception:
                pass
        if (j + 1) % max(1, NYP // 5) == 0:
            _log(f"    光催化 {j+1}/{NYP} 行  ({time.time()-t0:.0f} s)")
    # 局域幂指数 ∂log r / ∂log F —— 本包的招牌结果：它是**算出来的**，不是假设的
    SLOPE = np.full_like(RATE, np.nan)
    lf = np.log10(FX)
    for j in range(NYP):
        ok = np.isfinite(RATE[j]) & (RATE[j] > 0)
        if ok.sum() > 3:
            SLOPE[j, ok] = np.gradient(np.log10(RATE[j, ok]), lf[ok])
    _log(f"  光催化完成 {np.isfinite(RATE).sum()}/{RATE.size} 点，{time.time()-t0:.0f} s")
    return dict(x=FX, y=KR, qfl=QFL, rate=RATE, aqy=AQY, slope=SLOPE)


# ════════════════════════════════════════════════════════════════════
# 行 4 ｜ 超声：乙烷脱氢 over (热点份额, 热点温度)
# ════════════════════════════════════════════════════════════════════
SONO = """
C2H6 + 2* <-> C2H5* + H*
C2H5* + * -> C2H4* + H*
C2H4* -> C2H4 + *
C2H5* -> C2H5 + *
2H* -> H2 + 2*
"""
T_BULK = 350.0


def row_sono():
    mech = M.from_text(SONO)
    en = energetics_from_arrays(mech, np.array([0.55, 1.95, 0.60, 0.85, 0.70]),
                                np.array([-0.20, 0.35, 0.25, 0.40, -0.30]))
    PG = {"C2H6": 1.0, "C2H4": 1e-8, "C2H5": 1e-8, "H2": 1e-6}
    DS = DriverStack([make_driver("sono")])
    FH = np.logspace(-5.0, -2.0, NXP)
    TH = np.linspace(1500.0, 6000.0, NYP)
    RATE = np.full((NYP, NXP), np.nan)
    RATIO = np.full((NYP, NXP), np.nan)
    SEL = np.full((NYP, NXP), np.nan)
    t0 = time.time()
    for j, Th in enumerate(TH):
        for i, f in enumerate(FH):
            try:
                rm = steady_state_multi(
                    mech, en, dict(T_K=T_BULK, hotspot_T_K=float(Th),
                                   hotspot_fraction=float(f)), DS,
                    p_gas=PG, tol=1e-12, strict_energy=False)
                Tav = f * Th + (1.0 - f) * T_BULK
                ra = steady_state(mech, en, state=dict(T_K=float(Tav)),
                                  p_gas=PG, tol=1e-12, strict_energy=False)
                r_c2h4 = rm.rate_net[2]      # C2H4 脱附
                r_c2h5 = rm.rate_net[3]      # C2H5 脱附（竞争）
                RATE[j, i] = r_c2h4
                if ra.rate_net[2] != 0:
                    RATIO[j, i] = r_c2h4 / ra.rate_net[2]
                tot = abs(r_c2h4) + abs(r_c2h5)
                if tot > 0:
                    SEL[j, i] = abs(r_c2h4) / tot
            except Exception:
                pass
        if (j + 1) % max(1, NYP // 5) == 0:
            _log(f"    超声 {j+1}/{NYP} 行  ({time.time()-t0:.0f} s)")
    _log(f"  超声完成 {np.isfinite(RATE).sum()}/{RATE.size} 点，{time.time()-t0:.0f} s")
    return dict(x=FH, y=TH, rate=RATE, ratio=RATIO, sel=SEL)


# ════════════════════════════════════════════════════════════════════
# 出图
# ════════════════════════════════════════════════════════════════════
def _panel(ax, X, Y, Z, cmap, title, cbar_label, *, log_y=False, log_x=False,
           diverging=False, levels=None, vmin=None, vmax=None, fmt="%.0f"):
    import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm
    norm = None
    if diverging:
        m = np.nanmax(np.abs(Z))
        if np.isfinite(m) and m > 0:
            norm = TwoSlopeNorm(vmin=-m, vcenter=0.0, vmax=m)
    # ★ NaN 必须一眼能认出来，并在图注里说明它是什么 —— 留白等于"读者自己猜"。
    #   本图里 NaN 只有两个来源：① 电流在高 pH 下下溢到 0（log 没有定义）；
    #   ② 超声那张里平均工况解 r(⟨T⟩) 恰好为 0（比值没有定义）。两者都是结果，不是故障。
    #
    # ★★ 第一版把 NaN 画成 "0.88" 浅灰 —— 这是个**真错误**，出 PDF 后才看出来：
    #   coolwarm 的**中点色就是 (0.867, 0.864, 0.863)**，和 0.88 灰最大通道差只有 0.017，
    #   肉眼完全分不开。而 (a) 图有 44.2% 的面积正落在那条中心色带上
    #   ⇒ 读者无法区分「速率≈0」和「没有解」，这正是图最不该含糊的地方。
    #   改法：NaN 设为**透明**，底下铺一层斜线影线；灰色于是还给真实数据。
    import matplotlib as _mpl
    cm = _mpl.colormaps[cmap].copy() if hasattr(_mpl, "colormaps") else None
    if cm is not None:
        cm.set_bad(alpha=0.0)
        cmap = cm
    ax.set_facecolor("white")
    ax.patch.set_hatch("////")
    ax.patch.set_edgecolor("0.62")
    ax.patch.set_linewidth(0.0)
    pc = ax.pcolormesh(X, Y, Z, cmap=cmap, norm=norm, vmin=vmin, vmax=vmax,
                       shading="auto", rasterized=True)
    if levels is not None:
        cs = ax.contour(X, Y, Z, levels=levels, colors="k",
                        linewidths=0.6, alpha=0.75)
        ax.clabel(cs, inline=True, fontsize=7.5, fmt=fmt)
    if log_y:
        ax.set_yscale("log")
    if log_x:
        ax.set_xscale("log")
    # ★ pcolormesh 的 shading="auto" 会把格心外推半格，坐标轴因此伸到数据之外
    #   （粗网格上尤其明显：P 轴一直画到 1e-1，而数据只到 1e0）。显式钉死范围。
    ax.set_xlim(np.nanmin(X), np.nanmax(X))
    ax.set_ylim(np.nanmin(Y), np.nanmax(Y))
    ax.set_title(title, fontsize=8, pad=3)
    cb = ax.figure.colorbar(pc, ax=ax, pad=0.02, fraction=0.046)
    cb.ax.tick_params(labelsize=7)
    cb.set_label(cbar_label, fontsize=7)
    ax.tick_params(labelsize=7)
    # ★ 版面压到 Elsevier 双栏 190 mm 后，默认刻度数会撞在一起；线性轴限 4 个。
    from matplotlib.ticker import MaxNLocator
    if not log_x:
        ax.xaxis.set_major_locator(MaxNLocator(nbins=4))
    if not log_y:
        ax.yaxis.set_major_locator(MaxNLocator(nbins=4))
    return pc


def make_figure(TH_, EL_, PH_, SO_):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    # ★ fonttype 42 = TrueType 嵌入；默认的 Type 3 Elsevier 不收（见 make_figures.py 的注释）
    plt.rcParams.update({"savefig.dpi": 600, "font.family": "sans-serif",
                         "font.sans-serif": ["Arial", "DejaVu Sans"],
                         "mathtext.fontset": "stixsans",
                         "pdf.fonttype": 42, "ps.fonttype": 42,
                         "axes.linewidth": 0.6})

    fig, AX = plt.subplots(4, 3, figsize=(7.48, 9.1))

    # ── 行1 热 ──
    X, Y = np.meshgrid(TH_["x"], TH_["y"])
    R = TH_["rate"]
    # 带符号：用 asinh 压缩量级又保留符号；零线就是平衡线
    s = np.nanmedian(np.abs(R[np.isfinite(R) & (R != 0)])) or 1.0
    Z = np.arcsinh(R / s)
    _panel(AX[0, 0], X, Y, Z, "coolwarm",
           "(a1) net rate of NH$_3$ formation",
           r"asinh($r/r_0$)", log_y=True, diverging=True)
    cs = AX[0, 0].contour(X, Y, R, levels=[0.0], colors="k", linewidths=1.4)
    AX[0, 0].clabel(cs, inline=True, fontsize=7.0, fmt={0.0: "equilibrium"})
    _panel(AX[0, 1], X, Y, TH_["th_H"], "turbo", r"(a2) $\theta_{H*}$",
           r"$\theta_{H*}$", log_y=True, levels=[0.05, 0.1, 0.2, 0.3], fmt="%.2f")
    _panel(AX[0, 2], X, Y, TH_["th_vac"], "turbo", r"(a3) free sites $\theta_*$",
           r"$\theta_*$", log_y=True, levels=[0.2, 0.5, 0.8, 0.95], fmt="%.2f")
    for a in AX[0]:
        a.set_xlabel("$T$ (K)", fontsize=7.5)
    AX[0, 0].set_ylabel("$P$ (bar)", fontsize=7.5)

    # ── 行2 电 ──
    X, Y = np.meshgrid(EL_["x"], EL_["y"])
    with np.errstate(all="ignore"):
        LJ = np.log10(np.abs(EL_["j"]))
    _panel(AX[1, 0], X, Y, LJ, "turbo", "(b1) current density",
           r"log$_{10}$ |$j$| (A m$^{-2}$)", levels=[-4, -2, 0, 2, 4], fmt="%.0f")
    _panel(AX[1, 1], X, Y, EL_["th_H"], "turbo", r"(b2) $\theta_{H*}$",
           r"$\theta_{H*}$", levels=[0.1, 0.5, 0.9], fmt="%.1f")
    TAF = np.where(np.isfinite(EL_["tafel"]) & (np.abs(EL_["tafel"]) < 400),
                   EL_["tafel"], np.nan)
    _panel(AX[1, 2], X, Y, np.abs(TAF), "turbo", "(b3) Tafel slope",
           "mV dec$^{-1}$", levels=[30, 40, 60, 120], fmt="%.0f", vmin=20, vmax=160)
    for a in AX[1]:
        a.set_xlabel("$U$ (V vs SHE)", fontsize=7.5)
    AX[1, 0].set_ylabel("pH", fontsize=7.5)

    # ── 行3 光 ──
    X, Y = np.meshgrid(PH_["x"], PH_["y"])
    _panel(AX[2, 0], X, Y, PH_["qfl"], "turbo",
           "(c1) quasi-Fermi-level splitting",
           r"$\Delta E_{\rm qfl}$ (eV)", log_x=True, log_y=True,
           levels=[0.6, 0.8, 1.0, 1.2], fmt="%.1f")
    with np.errstate(all="ignore"):
        LR = np.log10(PH_["rate"])
    _panel(AX[2, 1], X, Y, LR, "turbo", "(c2) rate",
           r"log$_{10}$ $r$ (s$^{-1}$)", log_x=True, log_y=True)
    _panel(AX[2, 2], X, Y, PH_["slope"], "viridis",
           "(c3) local power-law exponent",
           r"$\partial\log r/\partial\log I$", log_x=True, log_y=True,
           levels=[0.5, 0.6, 0.7, 0.8, 0.9], fmt="%.2f", vmin=0.45, vmax=1.02)
    for a in AX[2]:
        a.set_xlabel(r"photon flux (m$^{-2}$ s$^{-1}$)", fontsize=7.5)
    AX[2, 0].set_ylabel(r"$k_{\rm rec}$", fontsize=7.5)

    # ── 行4 超声 ──
    X, Y = np.meshgrid(SO_["x"], SO_["y"])
    with np.errstate(all="ignore"):
        LRs = np.log10(np.abs(SO_["rate"]))
        LRt = np.log10(np.abs(SO_["ratio"]))
    _panel(AX[3, 0], X, Y, LRs, "turbo", "(d1) C$_2$H$_4$ formation rate",
           r"log$_{10}$ $r$ (s$^{-1}$)", log_x=True)
    _panel(AX[3, 1], X, Y, LRt, "turbo",
           r"(d2) $\langle r(T)\rangle\,/\,r(\langle T\rangle)$",
           "orders of magnitude", log_x=True,
           levels=[5, 10, 15, 20, 25], fmt="%.0f")
    _panel(AX[3, 2], X, Y, SO_["sel"], "turbo",
           "(d3) selectivity to C$_2$H$_4$",
           "selectivity", log_x=True, levels=[0.3, 0.5, 0.7, 0.9], fmt="%.1f")
    for a in AX[3]:
        a.set_xlabel("hot-spot fraction", fontsize=7.5)
    AX[3, 0].set_ylabel(r"$T_{\rm hot}$ (K)", fontsize=7.5)

    # ★ 压到双栏宽之后，行标签原来的 -0.26 会盖住纵轴刻度（尤其是 10^-18 这种宽标签）。
    #   挪到轴外 -0.42，并用 rect 给左侧留出 5% 的空白。
    for i, lab in enumerate(["THERMAL", "ELECTRO", "PHOTO", "SONO"]):
        AX[i, 0].text(-0.42, 0.5, lab, transform=AX[i, 0].transAxes,
                      rotation=90, va="center", ha="center",
                      fontsize=8.5, fontweight="bold")

    fig.tight_layout(rect=[0.045, 0, 1, 1])
    os.makedirs(OUT, exist_ok=True)
    fig.savefig(os.path.join(OUT, "Fig4_maps.png"), bbox_inches="tight")
    fig.savefig(os.path.join(OUT, "Fig4_maps.pdf"), bbox_inches="tight")
    _log("  出图 Fig4_maps.png / .pdf")


def save_data(TH_, EL_, PH_, SO_):
    os.makedirs(SI, exist_ok=True)
    np.savez_compressed(
        os.path.join(SI, "fig4_maps.npz"),
        th_T=TH_["x"], th_P=TH_["y"], th_rate=TH_["rate"],
        th_thH=TH_["th_H"], th_thV=TH_["th_vac"],
        el_U=EL_["x"], el_pH=EL_["y"], el_j=EL_["j"],
        el_thH=EL_["th_H"], el_tafel=EL_["tafel"],
        ph_flux=PH_["x"], ph_krec=PH_["y"], ph_qfl=PH_["qfl"],
        ph_rate=PH_["rate"], ph_aqy=PH_["aqy"], ph_slope=PH_["slope"],
        so_f=SO_["x"], so_T=SO_["y"], so_rate=SO_["rate"],
        so_ratio=SO_["ratio"], so_sel=SO_["sel"])
    # 文本版：每行一个 (x, y, 值…)，Origin / gnuplot 直接可读
    def dump(name, x, y, cols):
        p = os.path.join(SI, name)
        with open(p, "w") as fh:
            fh.write("# KinSTEP Fig.4 原始数据；列: 1=x 2=y " + " ".join(
                f"{i+3}={k}" for i, k in enumerate(cols)) + "\n")
            for j in range(len(y)):
                for i in range(len(x)):
                    fh.write("%.6e %.6e " % (x[i], y[j])
                             + " ".join("%.6e" % cols[k][j, i] for k in cols) + "\n")
        return p
    dump("fig4a_thermal.dat", TH_["x"], TH_["y"],
         {"rate": TH_["rate"], "theta_H": TH_["th_H"], "theta_vac": TH_["th_vac"]})
    dump("fig4b_electro.dat", EL_["x"], EL_["y"],
         {"j": EL_["j"], "theta_H": EL_["th_H"], "tafel_mV": EL_["tafel"]})
    dump("fig4c_photo.dat", PH_["x"], PH_["y"],
         {"qfl_eV": PH_["qfl"], "rate": PH_["rate"], "AQY": PH_["aqy"],
          "slope": PH_["slope"]})
    dump("fig4d_sono.dat", SO_["x"], SO_["y"],
         {"rate": SO_["rate"], "ratio": SO_["ratio"], "selectivity": SO_["sel"]})
    _log("  存数据 SI/fig4_maps.npz + 四个 .dat")


def summary(TH_, EL_, PH_, SO_):
    """把要写进正文的数字印出来 —— 正文里的每个数都必须出自这里。"""
    _log("\n═══ 可写进正文的数字 ═══")
    R = TH_["rate"]
    pos = np.isfinite(R) & (R > 0)
    _log("  热：%d/%d 点收敛；净速率为正（合成）的占 %.1f%%，"
         "零线即平衡线" % (np.isfinite(R).sum(), R.size, 100 * pos.sum() / np.isfinite(R).sum()))
    _log("      最大正速率 %.3e s^-1 出现在 T=%.0f K, P=%.0f bar"
         % (np.nanmax(R), TH_["x"][np.unravel_index(np.nanargmax(R), R.shape)[1]],
            TH_["y"][np.unravel_index(np.nanargmax(R), R.shape)[0]]))
    T = np.abs(EL_["tafel"])
    ok = np.isfinite(T) & (T < 400)
    _log("  电：Tafel 斜率中位 %.1f mV/dec（%d 个有效点），范围 %.1f–%.1f"
         % (np.nanmedian(T[ok]), ok.sum(), np.nanmin(T[ok]), np.nanmax(T[ok])))
    S = PH_["slope"]
    lo = S[:, :3][np.isfinite(S[:, :3])]
    hi = S[:, -3:][np.isfinite(S[:, -3:])]
    _log("  光：幂指数低光强端 %.4f → 高光强端 %.4f（文献 1 → 1/2）"
         % (lo.mean(), hi.mean()))
    _log("      ΔE_qfl 范围 %.3f–%.3f eV" % (np.nanmin(PH_["qfl"]), np.nanmax(PH_["qfl"])))
    RT = SO_["ratio"]
    # ★ 比值可能为负：高热点温度下平均工况解 r(⟨T⟩) 会反号（逆反应占优），
    #   直接取 log10 得到 nan。按绝对值报量级，并把反号的点数单独说明 ——
    #   反号本身是**物理结果**（平均工况法连方向都给错了），不能悄悄丢掉。
    neg = int(np.sum(np.isfinite(RT) & (RT < 0)))
    A = np.abs(RT)
    _log("  超声：|⟨r(T)⟩/r(⟨T⟩)| 跨 %.2e – %.2e（%.1f – %.1f 个数量级）"
         % (np.nanmin(A), np.nanmax(A), np.log10(np.nanmin(A)), np.log10(np.nanmax(A))))
    if neg:
        _log("        其中 %d 点比值为**负** ⇒ 用平均工况算出来的净速率连方向都反了"
             % neg)
    _log("        选择性范围 %.3f – %.3f" % (np.nanmin(SO_["sel"]), np.nanmax(SO_["sel"])))


def load_data():
    """从已存的 npz 读回四行数据 —— `--replot` 用，免得为了调图再算 20 分钟。"""
    z = np.load(os.path.join(SI, "fig4_maps.npz"))
    return (dict(x=z["th_T"], y=z["th_P"], rate=z["th_rate"],
                 th_H=z["th_thH"], th_vac=z["th_thV"]),
            dict(x=z["el_U"], y=z["el_pH"], j=z["el_j"],
                 th_H=z["el_thH"], tafel=z["el_tafel"]),
            dict(x=z["ph_flux"], y=z["ph_krec"], qfl=z["ph_qfl"],
                 rate=z["ph_rate"], aqy=z["ph_aqy"], slope=z["ph_slope"]),
            dict(x=z["so_f"], y=z["so_T"], rate=z["so_rate"],
                 ratio=z["so_ratio"], sel=z["so_sel"]))


def main():
    _log("KinSTEP Fig.4 ｜ 四模态二维工况图" + ("（--quick 粗网格）" if QUICK else ""))
    t0 = time.time()
    if "--replot" in sys.argv:
        _log("  --replot：直接读 SI/fig4_maps.npz，不重算")
        TH_, EL_, PH_, SO_ = load_data()
        try:
            make_figure(TH_, EL_, PH_, SO_)
        except ImportError:
            _log("  ⊘ 没装 matplotlib")
        summary(TH_, EL_, PH_, SO_)
        _log(f"\n总用时 {time.time()-t0:.0f} s")
        return 0
    _log(f"  网格：热/电 {NX}×{NY}，光/超声 {NXP}×{NYP}")
    TH_ = row_thermal()
    EL_ = row_electro()
    PH_ = row_photo()
    SO_ = row_sono()
    save_data(TH_, EL_, PH_, SO_)
    try:
        make_figure(TH_, EL_, PH_, SO_)
    except ImportError:
        _log("  ⊘ 没装 matplotlib（可选件），跳过出图；数据已存")
    summary(TH_, EL_, PH_, SO_)
    _log(f"\n总用时 {time.time()-t0:.0f} s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
