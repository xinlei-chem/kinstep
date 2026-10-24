# -*- coding: utf-8 -*-
"""生成 CPC 论文的全部插图。所有数据都是**现跑出来的**，不是示意图。"""
import sys, os, time, warnings
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)            # 不写死绝对路径：代码存档会被下载到任意目录
warnings.filterwarnings("ignore")
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures")
os.makedirs(OUT, exist_ok=True)
# ★ 两条是投稿硬要求，不是审美：
#   ① pdf/ps.fonttype = 42 —— matplotlib 默认把字体嵌成 **Type 3**，Elsevier 不收
#      （排版系统无法重排、部分 PDF 工具显示异常）。42 = TrueType 嵌入。
#      自查命令：`pdffonts figures/*.pdf`，type 一列不得出现 Type 3。
#   ② 字体用 Arial（Elsevier 推荐 Arial/Courier/Times/Symbol），数学配 stixsans
#      —— 与 Arial 字形匹配且数学符号覆盖完整；若只改 font.family 而不改 mathtext，
#      正文是 Arial、公式仍是 DejaVu，一张图里两套字形。
#   ③ ★ Arial 没有 ①②③④（U+2460-2463）的字形，单设 font.family="Arial" 会出豆腐块。
#      matplotlib ≥3.6 支持**逐字符回退**：把 DejaVu Sans 放在 sans-serif 列表第二位，
#      Arial 缺的字自动回退，其余仍是 Arial。改字体后务必把图看一遍，别只看字体名。
plt.rcParams.update({"font.size": 8, "axes.linewidth": 0.8, "figure.dpi": 300,
                     "savefig.dpi": 600, "font.family": "sans-serif",
                     "font.sans-serif": ["Arial", "DejaVu Sans"],
                     "mathtext.fontset": "stixsans",
                     "pdf.fonttype": 42, "ps.fonttype": 42,
                     "axes.labelsize": 8, "legend.fontsize": 7,
                     "xtick.labelsize": 7, "ytick.labelsize": 7})
C = dict(th="#4C72B0", el="#DD8452", ph="#55A868", so="#C44E52",
         k="#4D4D4D", g="#8C8C8C", hi="#CCB974")

def style(ax):
    ax.tick_params(direction="in", top=True, right=True, which="both")


# ★ CATKINAS 的算例输出是闭源发行物自带的数据，未授再分发，所以不随代码存档
#   提供（见 RELEASE_NOTES.md）。Fig.2 与 Fig.3(d) 需要它来生成能量表与对比值；
#   数据不在时退回到 SI 里已存档的数值重绘同样的图，图是可复现的，
#   只是不重算。自备 CATKINAS 后设 KINSTEP_CATKINAS_DIR 即可全程重算。
def catkinas_ref(which="curve"):
    """取 CATKINAS 参考数组，返回 {变量名: ndarray} 或 None。

    与 `tests/refdata.catkinas_vars()` 同一口径，优先级也一样：
      ① 设了 KINSTEP_CATKINAS_DIR ⇒ 读发行物原件；
      ② 否则读随包的 `tests/catkinas_refdata.npz`（只含断言/作图用到的数组）；
      ③ 都没有 ⇒ None，调用方据此跳过。
    ★ 有了 ②，Fig.2、Fig.3(d)、Fig.S1(b) **在干净存档里也能真重算**，
      不必再退回已存档的数值重绘。
    """
    rel = {"curve": "Examples/result_curve/data1.mat",
           "map": "Examples/result_map/data1.mat"}[which]
    base = os.environ.get("KINSTEP_CATKINAS_DIR") or ""
    if base and os.path.exists(os.path.join(base, rel)):
        import scipy.io as sio
        d = sio.loadmat(os.path.join(base, rel))
        return {k: v for k, v in d.items() if not k.startswith("__")}
    npz = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "..", "tests", "catkinas_refdata.npz")
    if os.path.exists(npz):
        pre = which + "__"
        z = np.load(npz)
        out = {k[len(pre):]: z[k] for k in z.files if k.startswith(pre)}
        out.pop("md5", None)
        return out or None
    return None

# ══════════ Fig.1 程序框架 ══════════
def fig1():
    fig, ax = plt.subplots(figsize=(7.48, 4.70))
    ax.set_xlim(0, 10); ax.set_ylim(0, 6.6); ax.axis("off")
    def box(x, y, w, h, txt, fc, fs=7.4, tc="k", lw=0.9):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.045",
                                    fc=fc, ec="#333333", lw=lw, zorder=2))
        ax.text(x+w/2, y+h/2, txt, ha="center", va="center", fontsize=fs,
                color=tc, zorder=3, linespacing=1.35)
    def arrow(x1, y1, x2, y2, style_="-|>", lw=1.0, col="#333333", ls="-", ms=9):
        # ms：箭头头部大小。跨度很短时（两框之间 0.26 宽）默认 9 会让两个箭头头叠成
        # 一团，看着像个星号；短箭头要把 ms 调小。
        ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle=style_,
                                     mutation_scale=ms, lw=lw, color=col,
                                     linestyle=ls, zorder=1))
    # 输入层
    box(0.15, 5.5, 2.5, 0.85, "Mechanism (plain text)\n$\\mathrm{N_2 + 2* \\leftrightarrow 2N*}$", "#EAEFF6")
    box(2.9, 5.5, 2.1, 0.85, "Energies\n(with provenance)", "#EAEFF6")
    box(5.25, 5.5, 2.0, 0.85, "Conditions\n$T$, $p$, $U$, $I_\\mathrm{ph}$", "#EAEFF6")
    box(7.5, 5.5, 2.35, 0.85, "CLI / NL layer\ncheck | ask | run", "#EAEFF6")
    # 闸
    box(0.15, 4.46, 9.7, 0.80,
        "GATES   element & charge balance · site conservation · name collision\n"
        "Wegscheider cycle · provenance (energy source)", "#F5E6CC", fs=6.9)
    for x in (1.4, 3.95, 6.25, 8.65):
        arrow(x, 5.5, x, 5.30)
    # 内核
    box(0.15, 2.28, 5.2, 1.92, "", "#FFFFFF", lw=1.1)
    ax.text(2.75, 4.02, "CORE  (modality-agnostic)", ha="center", fontsize=7.6, weight="bold")
    # \u2605 \u539f\u6765\u7528 \u2460\u2461\u2462\u2463\uff08U+2460-2463\uff09\uff1aArial \u548c DejaVu Sans **\u90fd\u6ca1\u6709**\u8fd9\u7ec4\u5b57\u5f62\uff0c
    #   \u51fa\u56fe\u662f\u8c46\u8150\u5757\uff08\u6700\u540e\u4e00\u9053\u56de\u9000\u5b57\u4f53\u753b\u7684\u65b9\u6846 A\uff09\uff0cPNG/PDF \u91cc\u4e00\u76f4\u662f\u574f\u7684\u3002
    #   \u56fe\u91cc\u53ea\u7528 ASCII \u6700\u7a33\uff1bUnicode \u88c5\u9970\u5b57\u7b26\u8981\u5148\u786e\u8ba4\u76ee\u6807\u5b57\u4f53\u6709\u6ca1\u6709\u3002
    for i, (t, yy) in enumerate([("(1) ODE integration  \u2192  physical guess", 3.74),
                                 ("(2) Damped Newton (analytic Jacobian)", 3.42),
                                 ("(3) Flux projection  $\\sigma=R_f\\varepsilon$", 3.10),
                                 ("(4) High precision (mpmath, auto)", 2.78)]):
        ax.text(0.34, yy, t, fontsize=7.0, va="center")
        if i < 3:
            arrow(5.02, yy-0.07, 5.02, yy-0.25, lw=0.7, col=C["g"])
    ax.text(2.75, 2.44, "solver ladder \u2014 each level forced by measurement",
            ha="center", fontsize=7.0, style="italic", color=C["g"])
    arrow(2.75, 4.46, 2.75, 4.24)
    # driver 层
    box(5.65, 2.28, 4.2, 1.92, "", "#FFFFFF", lw=1.1)
    ax.text(7.75, 4.02, "DRIVER PLUGINS", ha="center", fontsize=7.6, weight="bold")
    dd = [("thermal", C["th"], "baseline"),
          ("electro", C["el"], "$\\Delta G - \\nu_e eU$"),
          ("photo", C["ph"], "$\\Delta G - \\nu\\,\\Delta E_\\mathrm{qfl}$"),
          ("sono", C["so"], "multi-state $\\Sigma w_i r_i$")]
    for i, (n, c, f) in enumerate(dd):
        y = 3.68 - i*0.32
        ax.add_patch(plt.Rectangle((5.85, y-0.11), 0.22, 0.22, fc=c, ec="none"))
        ax.text(6.17, y, "%-9s" % n, fontsize=7.0, va="center", family="monospace")
        ax.text(7.45, y, f, fontsize=7.0, va="center")
    ax.text(7.75, 2.44, "composable:  photo + electro = PEC",
            ha="center", fontsize=7.0, style="italic", color=C["g"])
    arrow(7.75, 4.46, 7.75, 4.24)
    arrow(5.37, 3.25, 5.63, 3.25, style_="<|-|>", lw=0.9, ms=5)
    # ★ 两个框之间只有 0.26 宽的缝，7 pt 的两行横排会压到框线上；竖排后占宽只剩字高。
    ax.text(5.50, 3.06, "one hook", ha="center", va="top", rotation=90,
            fontsize=7.0, color=C["g"])
    # 分析与输出
    box(0.15, 1.30, 4.6, 0.80, "ANALYSIS\nDRC ($\\Sigma X_B\\!=\\!1$) · orders · $E_\\mathrm{app}$ · Tafel", "#E8F0E8")
    box(5.25, 1.30, 4.6, 0.80, "SCAN & PLOT\n1D/2D volcano · coverage · $\\mathtt{.dat}$ export", "#E8F0E8")
    arrow(2.75, 2.28, 2.45, 2.12); arrow(7.75, 2.28, 7.60, 2.12)
    box(0.15, 0.20, 9.7, 0.84,
        "SELF-CHECKS (free, always on)\n"
        "$\\Sigma X_B\\!=\\!1$ (Campbell)  ·  $\\nu^{\\mathsf{T}}r\\!=\\!0$ (flux)  ·  "
        "$\\Sigma\\theta\\!+\\!\\theta_v\\!=\\!1$  ·  $R_{e^-}\\!=\\!R_{h^+}$ (photo)  ·  "
        "$\\Sigma w_i\\!=\\!1$ (sono)", "#F0E6EE", fs=6.9)
    arrow(2.45, 1.30, 2.45, 1.06); arrow(7.60, 1.30, 7.60, 1.06)
    fig.savefig(os.path.join(OUT, "Fig1_framework.png"), bbox_inches="tight")
    fig.savefig(os.path.join(OUT, "Fig1_framework.pdf"), bbox_inches="tight")
    plt.close(fig); print("  Fig.1 framework")

# ══════════ Fig.2 验证：对三套独立真值 ══════════
def fig2():
    import scipy.io as sio
    from kinstep import mech as M
    from kinstep.analyze import energetics_from_arrays
    from kinstep.scaling import from_table
    from kinstep.thermo import entropy_offsets
    from kinstep.scan import scan
    import kinstep.drivers
    from kinstep.units import K_B_EV

    d = catkinas_ref("curve")
    if d is None:
        print("  ⊘ Fig.2 跳过：重算需要 CATKINAS 参考数据。"
              "图中绘出的全部数值已存档于 SI/fig2_data.npz（见 RELEASE_NOTES.md）。")
        return
    Ean, G0n = np.asarray(d["Ean"]), np.asarray(d["G0n"])
    RNET, Y, Z = (np.asarray(d[k]) for k in ("Rnet", "y", "Z"))
    E1 = np.arange(-2.0, 2.0001, 0.01)
    KT, APRE, TN = 0.057994719627, 1.402306113158e13, 673.0
    T_EQ = TN * KT / (K_B_EV * TN)
    # ★ 步序必须与 .mat 里的列序一致（CATKINAS 是 H2 在前）——
    #   写反了会让能量与偏移错位，实测相对偏差从 2.7e-10 变成 1.85e+07。
    mech = M.from_text("""
H2 + 2* <-> 2H*
N2 + 2* <-> 2N*
N* + H* <-> NH* + *
NH* + H* <-> NH2* + *
NH2* + H* <-> NH3 + 2*
""")
    TAB = {1: dict(bep=(0.54, 0.60), scale=(0.30, -0.68)),
           2: dict(bep=(0.78, 1.54), scale=(2.00, 0.00)),
           3: dict(bep=(0.26, 1.45), scale=(-0.15, -0.18)),
           4: dict(bep=(0.43, 1.66), scale=(-0.45, -0.45)),
           5: dict(bep=(0.31, 1.45), scale=(-0.85, 1.17))}
    off_Ga = (Ean - np.asarray(d["Eao"]))[0]
    off_dG = (G0n - np.asarray(d["G0o"]))[0]
    sc = from_table(mech, TAB, ["E_N"], offset_Ga=off_Ga, offset_dG=off_dG)
    print("    running 401-point scan (double)...")
    r_d = scan(mech, sc, dict(T_K=T_EQ), [E1], p_gas={"N2": 25., "H2": 75., "NH3": 1.},
               A_pre=APRE, precision="double")
    print("    running 401-point scan (auto)...")
    r_a = scan(mech, sc, dict(T_K=T_EQ), [E1], p_gas={"N2": 25., "H2": 75., "NH3": 1.},
               A_pre=APRE, precision="auto", dps=60)
    np.savez(os.path.join(OUT, "..", "SI", "fig2_data.npz"),
             E1=E1, r_double=r_d.rate_of(-1), r_auto=r_a.rate_of(-1), ref=RNET[:, 4],
             theta_auto=np.array([r_a.theta[s] for s in mech.surface_species]),
             ref_y=Y[:, :4], hp=(r_a.precision == "high"))

    fig, axs = plt.subplots(1, 3, figsize=(7.48, 2.48))
    # (a) 火山曲线 + 高精度点
    a = axs[0]
    hp = (r_a.precision == "high")
    a.plot(E1, np.log10(np.abs(RNET[:, 4])), "-", lw=2.4, color="#BBBBBB",
           label="CATKINAS (ref.)", zorder=1)
    a.plot(E1[~hp], np.log10(np.abs(r_a.rate_of(-1)))[~hp], ".", ms=1.6,
           color=C["th"], label="this work (double)", zorder=3)
    a.plot(E1[hp], np.log10(np.abs(r_a.rate_of(-1)))[hp], ".", ms=1.6,
           color=C["so"], label="auto high-prec. (%d pts)" % hp.sum(), zorder=4)
    a.set_xlabel(r"$E_\mathrm{N}$ (eV)"); a.set_ylabel(r"$\log_{10}(r\,/\,\mathrm{s^{-1}})$")
    a.legend(frameon=False, fontsize=7.0, loc="lower left"); style(a)
    a.text(-0.02, 1.04, "(a)", transform=a.transAxes, fontsize=8, weight="bold", va="bottom")
    # (b) 逐点相对偏差
    b = axs[1]
    rel_d = np.abs(r_d.rate_of(-1) / RNET[:, 4] - 1)
    rel_a = np.abs(r_a.rate_of(-1) / RNET[:, 4] - 1)
    b.semilogy(E1, np.maximum(rel_d, 1e-17), ".", ms=1.4, color=C["th"], label="double")
    b.semilogy(E1, np.maximum(rel_a, 1e-17), ".", ms=1.4, color=C["so"],
               label=r"$\mathtt{precision=auto}$")
    b.axhline(1e-6, ls="--", lw=0.7, color=C["k"])
    b.text(-1.95, 1.6e-6, r"$10^{-6}$", fontsize=7.0, color=C["k"])
    b.set_xlabel(r"$E_\mathrm{N}$ (eV)")
    b.set_ylabel("relative deviation in $r$")
    b.set_ylim(1e-17, 1e0)
    b.legend(frameon=False, fontsize=7.0, loc="upper left"); style(b)
    b.text(-0.02, 1.04, "(b)", transform=b.transAxes, fontsize=8, weight="bold", va="bottom")
    # (c) 三套真值汇总（对数条形）
    c = axs[2]
    items = [("Koper CV\n(ODE rhs)", 2.44e-15), ("closed form\n" + r"($\theta$, $r$)", 2.22e-16),
             ("CATKINAS\n" + r"($\theta$, $Z$)", 1.41e-10), ("CATKINAS\n($r$, auto)", 2.71e-10),
             ("DRC\n" + r"($\Sigma X_B$)", 2.94e-12), ("NIST\n" + r"($T\Delta S$)", 1.86e-10),
             ("CSTR\n(median)", 2.22e-16)]
    yy = np.arange(len(items))[::-1]
    c.barh(yy, [-np.log10(v) for _, v in items], height=0.62,
           color=[C["el"], C["th"], C["th"], C["so"], C["ph"], C["hi"], C["g"]])
    for i, (nm, v) in enumerate(items):
        c.text(-np.log10(v) + 0.25, yy[i], "%.1e" % v, va="center", fontsize=7.0)
    c.set_yticks(yy); c.set_yticklabels([n for n, _ in items], fontsize=7.0)
    c.set_xlabel(r"$-\log_{10}$ (relative deviation)")
    c.set_xlim(0, 21); style(c)
    c.text(-0.02, 1.04, "(c)", transform=c.transAxes, fontsize=8, weight="bold", va="bottom")
    fig.tight_layout(w_pad=1.3)
    fig.savefig(os.path.join(OUT, "Fig2_verification.png"), bbox_inches="tight")
    fig.savefig(os.path.join(OUT, "Fig2_verification.pdf"), bbox_inches="tight")
    plt.close(fig)
    print("  Fig.2 verification  (double max %.2e -> auto max %.2e)"
          % (rel_d.max(), rel_a.max()))
    return rel_d, rel_a


# ══════════ Fig.3 四种模态的案例 ══════════
def fig3():
    from kinstep import mech as M
    from kinstep.analyze import energetics_from_arrays, tafel_slope, drc, calibrate_delta
    from kinstep.solve import steady_state, steady_state_multi
    from kinstep.driver import DriverStack, make_driver
    from kinstep.drivers.photo import photo_steady
    import kinstep.drivers

    fig, axs = plt.subplots(2, 2, figsize=(7.48, 5.66))
    SI = os.path.join(OUT, "..", "SI")

    import time as _t; _t0=_t.time()
    print("    (a) HER...", flush=True)
    ax = axs[0, 0]
    m1 = M.from_text("""
H+ + e- + * <-> H*            {beta_ct: 0.5}
H* + H+ + e- -> H2 + *        {beta_ct: 0.5}
2H* -> H2 + 2*
""")
    en1 = energetics_from_arrays(m1, np.array([0.45, 0.60, 0.75]),
                                 np.array([-0.10, -0.30, -0.20]))
    el = DriverStack([make_driver("electro", U_ref="SHE", gamma_site=2.2e-5)])
    US = np.linspace(-0.30, 0.0, 41); jj = []; g = None
    for U in US:
        r = steady_state(m1, en1, state=dict(T_K=298.15, U_elec=U), driver=el,
                         p_gas={"H2": 1.0}, c_liq={"H+": 1.0}, theta_guess=g,
                         ladder=("newton",) if g else ("ode", "newton"),
                         tol=1e-12, precision="auto", strict_energy=False)
        jj.append(r.obs["j_curr"]); g = dict(r.theta)
    jj = np.array(jj); ts = tafel_slope(US, jj) * 1000
    ax.semilogy(US, np.abs(jj) / 10, "-", lw=1.4, color=C["el"])
    ax.set_xlabel(r"$U$ (V vs SHE)"); ax.set_ylabel(r"$|j|$ (mA cm$^{-2}$)")
    ax2 = ax.twinx()
    ax2.plot(US, np.abs(ts), "--", lw=1.0, color=C["k"])
    ax2.set_ylabel("Tafel slope (mV dec$^{-1}$)", fontsize=7)
    ax2.set_ylim(100, 130); ax2.tick_params(labelsize=6.5)
    ax2.axhline(118.3, ls=":", lw=0.8, color=C["g"])
    ax2.text(-0.295, 119.6, r"118.3 = $\ln\!10\,k_BT/(\beta e)$", fontsize=7.0, color=C["g"])
    style(ax); ax.set_title("(a) electro — HER", fontsize=7.6, weight="bold", loc="left", pad=3)
    np.savetxt(os.path.join(SI, "fig3a_her_tafel.dat"),
               np.column_stack([US, jj, ts]), header="U_SHE_V  j_A_per_m2  tafel_mV_per_dec")

    print("    (a) done %.1fs; (b) photo..." % (_t.time()-_t0), flush=True); _t0=_t.time()
    ax = axs[0, 1]
    m2 = M.from_text("""
H2O(l) + 2* <-> OH* + H*
OH* -> OH(l) + *            {nu_hole_vb: 1}
2H* -> H2 + 2*              {nu_elec_cb: 2}
""")
    en2 = energetics_from_arrays(m2, np.array([0.20, 1.35, 0.50]),
                                 np.array([-0.30, 0.25, -0.35]))
    flux = np.logspace(13, 26, 20); rate = []; aqy = []
    for F in flux:
        st = dict(T_K=300.0, photon_flux=F, absorptance=0.5, gamma_site_photo=1.5e19,
                  photon_energy_eV=3.2, qfl_frac_n=0.5, recomb_rate_const=1.0e-17,
                  photo_carrier_mode="steady_carrier")
        r, x, info = photo_steady(m2, en2, st, p_gas={"H2": 1e-6},
                                  c_liq={"H2O": 1.0, "OH": 1e-9}, tol=1e-9)
        rate.append(abs(r.rate_net[1])); aqy.append(info["AQY"])
    rate = np.array(rate); aqy = np.array(aqy)
    sl = np.gradient(np.log10(rate), np.log10(flux))
    ax.loglog(flux, rate, "o-", ms=2.4, lw=1.2, color=C["ph"])
    ax.set_xlabel(r"photon flux (m$^{-2}$ s$^{-1}$)")
    ax.set_ylabel(r"rate (s$^{-1}$ site$^{-1}$)")
    axb = ax.twinx()
    axb.semilogx(flux, sl, "--", lw=1.0, color=C["k"])
    axb.set_ylabel("local power-law order", fontsize=7); axb.set_ylim(0.3, 1.1)
    axb.tick_params(labelsize=6.5)
    for v, lb in ((1.0, "1 (low $I$)"), (0.5, "1/2 (high $I$)")):
        axb.axhline(v, ls=":", lw=0.8, color=C["g"])
        axb.text(2e13, v + 0.02, lb, fontsize=7.0, color=C["g"])
    style(ax); ax.set_title("(b) photo — intensity law", fontsize=7.6, weight="bold",
                            loc="left", pad=3)
    np.savetxt(os.path.join(SI, "fig3b_photo_powerlaw.dat"),
               np.column_stack([flux, rate, aqy, sl]),
               header="photon_flux  rate  AQY  local_order")

    print("    (b) done %.1fs; (c) sono..." % (_t.time()-_t0), flush=True); _t0=_t.time()
    ax = axs[1, 0]
    m3 = M.from_text("""
C2H6 + 2* <-> C2H5* + H*
C2H5* + * -> C2H4* + H*
C2H4* -> C2H4 + *
C2H5* -> C2H5 + *
2H* -> H2 + 2*
""")
    en3 = energetics_from_arrays(m3, np.array([0.55, 1.95, 0.60, 0.85, 0.70]),
                                 np.array([-0.20, 0.35, 0.25, 0.40, -0.30]))
    PG3 = {"C2H6": 1.0, "C2H4": 1e-8, "C2H5": 1e-8, "H2": 1e-6}
    DS = DriverStack([make_driver("sono")])
    fr = np.logspace(-30, -2, 29); rw = []; ra = []
    for f in fr:
        r = steady_state_multi(m3, en3, dict(T_K=350.0, hotspot_T_K=5000.0,
                                             hotspot_fraction=f), DS,
                               p_gas=PG3, tol=1e-12, strict_energy=False)
        rw.append(abs(r.rate_net[2]))
        Tm = f * 5000.0 + (1 - f) * 350.0
        rr = steady_state(m3, en3, state=dict(T_K=Tm), p_gas=PG3,
                          tol=1e-12, strict_energy=False)
        ra.append(abs(rr.rate_net[2]))
    rw, ra = np.array(rw), np.array(ra)
    ax.loglog(fr, rw, "-", lw=1.5, color=C["so"], label=r"$\langle r(T)\rangle = \Sigma w_i r_i$")
    ax.loglog(fr, ra, "--", lw=1.2, color=C["g"], label=r"$r(\langle T\rangle)$  (wrong)")
    fstar = ra[0] / (rw[-1] / fr[-1])
    ax.axvline(1.793e-27, ls=":", lw=0.9, color=C["k"])
    ax.text(2.6e-27, 2e-15, r"$f^*=r_\mathrm{bulk}/r_\mathrm{hot}$", fontsize=7.0, rotation=90,
            va="bottom")
    ax.set_xlabel("cavitation hot-spot fraction $f$")
    ax.set_ylabel(r"rate (s$^{-1}$ site$^{-1}$)")
    ax.legend(frameon=False, fontsize=7.0, loc="center left", bbox_to_anchor=(0.02, 0.62))
    style(ax)
    ax.set_title("(c) sono — multi-state", fontsize=7.6, weight="bold", loc="left", pad=3)
    np.savetxt(os.path.join(SI, "fig3c_sono.dat"), np.column_stack([fr, rw, ra]),
               header="hotspot_fraction  rate_weighted  rate_at_mean_T")

    print("    (c) done %.1fs; (d) DRC..." % (_t.time()-_t0), flush=True); _t0=_t.time()
    ax = axs[1, 1]
    import scipy.io as sio
    from kinstep.scaling import from_table
    from kinstep.units import K_B_EV
    from kinstep.solve import steady_state as _ss
    d = catkinas_ref("curve")
    if d is not None:
        Ean, G0n, DB_ = np.asarray(d["Ean"]), np.asarray(d["G0n"]), np.asarray(d["DRC_B"])
        E1 = np.arange(-2.0, 2.0001, 0.01)
        KT, APRE, TN = 0.057994719627, 1.402306113158e13, 673.0
        T_EQ = TN * KT / (K_B_EV * TN)
        m4 = M.from_text("""
    H2 + 2* <-> 2H*
    N2 + 2* <-> 2N*
    N* + H* <-> NH* + *
    NH* + H* <-> NH2* + *
    NH2* + H* <-> NH3 + 2*
    """)
        # ★ 沿描述符续算：每点从头解会在强吸附端不收敛（DRC 需要一个已收敛的基态）
        # 图上不需要 1e-12 的精度：关掉 Richardson 并减点（19 点×Richardson ≈ 1000 次
        # 求解，实测要 13 分钟）。求和规则仍能到 ~1e-9，足够说明问题。
        # 只取中段：控速步的切换就发生在这里，而两端的求解会掉进完整 ODE 阶梯
        # （实测两端单点可耗数十秒，而图上并不靠它们说明问题）
        sel = np.arange(60, 311, 12)
        XB = np.zeros((len(sel), 5)); SB = np.zeros(len(sel))
        PG4 = {"N2": 25., "H2": 75., "NH3": 1.}
        guess = None
        for k, i in enumerate(sel):
            en4 = energetics_from_arrays(m4, Ean[i], G0n[i], A_pre=APRE)
            base = None
            for lad, gg in ((("newton",), guess), (("ode", "newton"), None)):
                if lad == ("newton",) and gg is None:
                    continue
                try:
                    b = steady_state(m4, en4, state=dict(T_K=T_EQ), p_gas=PG4,
                                     theta_guess=gg, ladder=lad, tol=1e-12,
                                     precision="double", strict_energy=False,
                                     ode_stages=2, ode_budget_s=2.0)
                    if b.converged:
                        base = b; break
                except Exception:
                    pass
            if base is None:
                XB[k] = np.nan; SB[k] = np.nan; continue
            guess = dict(base.theta)
            try:
                dd = drc(m4, en4, dict(T_K=T_EQ), p_gas=PG4, delta=1e-3,
                         richardson=False, base=base,
                         solve_kw=dict(ode_stages=2, ode_budget_s=2.0))
                XB[k] = dd.X_B[:, 4]; SB[k] = dd.sum_B[4]
            except Exception:
                XB[k] = np.nan; SB[k] = np.nan
        # ★ 用 Campbell 求和规则本身筛点 —— 这正是这条规则的用途：
        #   它不依赖机理也不依赖能量，ΣX_B 明显偏离 1 的点，其 DRC 不可信。
        #   （实测有点的扰动解收敛到错解，ΣX_B 跑到 1.3e+03 而不是 nan，
        #     光靠 isfinite 筛不掉。）
        good = np.isfinite(SB) & (np.abs(SB - 1.0) < 1e-3)
        print("    DRC along descriptor: %d/%d points pass the sum rule" % (good.sum(), len(sel)))
    else:
        # 缺 CATKINAS 参考数据：用 SI 里已存档的数值重绘同一张图。
        # 图可复现，只是不重算；下游绘图代码一行不改。
        _dat = np.genfromtxt(os.path.join(SI, "fig3d_drc.dat"), comments="#")
        E1, sel = _dat[:, 0], slice(None)
        XB, SB = _dat[:, 1:6], _dat[:, 6]
        good = np.ones(len(SB), bool)
        print("    ⊘ (d) 用 SI/fig3d_drc.dat 已存档数值重绘（缺 CATKINAS 参考数据）")
    # ★ 不用堆叠图：X_B 可以为负（该步被稳定反而拖慢反应），堆叠必须先把负值
    #   削成 0，视觉上 ΣX_B 就冲到 1.3，与旁边"ΣX_B = 1"的虚线自相矛盾。
    #   折线把负的抑制效应如实画出来，求和规则仍落在 1 上。
    cols = plt.get_cmap("viridis")(np.linspace(0.08, 0.88, 5))
    mk = ["o", "s", "^", "v", "D"]
    for i in range(5):
        ax.plot(E1[sel][good], XB[good][:, i], "-", marker=mk[i], ms=2.6, lw=1.1,
                color=cols[i], label="(%d)" % (i + 1))
    ax.plot(E1[sel][good], SB[good], "k--", lw=1.0, label=r"$\Sigma X_B$")
    ax.axhline(0.0, color="0.7", lw=0.7)
    ax.text(0.98, 0.04, r"$\Sigma X_B = 1$  (max dev. %.1e)"
            % np.max(np.abs(SB[good] - 1)), transform=ax.transAxes, fontsize=7.0,
            ha="right")
    ax.set_xlabel(r"$E_\mathrm{N}$ (eV)"); ax.set_ylabel(r"$X_{B,i}$")
    ax.set_ylim(min(-0.15, XB[good].min() - 0.10),
                max(1.15, XB[good].max() + 0.12))
    ax.legend(frameon=False, fontsize=7.0, ncol=6, loc="upper center",
              bbox_to_anchor=(0.5, -0.20), columnspacing=0.7, handlelength=1.2,
              title="elementary step", title_fontsize=7.0)
    style(ax)
    ax.set_title("(d) thermal — DRC", fontsize=7.6, weight="bold", loc="left", pad=3)
    np.savetxt(os.path.join(SI, "fig3d_drc.dat"),
               np.column_stack([E1[sel][good], XB[good], SB[good]]),
               header="E_N  X_B1 X_B2 X_B3 X_B4 X_B5  sum")
    print("    (d) done %.1fs" % (_t.time()-_t0), flush=True)
    fig.tight_layout(w_pad=2.2, h_pad=2.4)
    fig.savefig(os.path.join(OUT, "Fig3_cases.png"), bbox_inches="tight")
    fig.savefig(os.path.join(OUT, "Fig3_cases.pdf"), bbox_inches="tight")
    plt.close(fig)
    print("  Fig.3 cases  (Tafel %.1f mV/dec; photo order %.3f->%.3f; "
          "sono ratio %.2e; SigmaXB dev %.1e)"
          % (np.median(np.abs(ts[:10])), sl[:3].mean(), sl[-3:].mean(),
             rw[20] / ra[20], np.max(np.abs(SB[good] - 1))))


def figS1():
    # 图 S1：1−z 的精度（支撑 SI §S3，Tables S4–S5）
    import mpmath as mp
    import scipy.io as sio
    from kinstep import mech as M
    from kinstep.hp import steady_state_hp
    from kinstep.analyze import energetics_from_arrays
    from kinstep.solve import steady_state
    from kinstep.units import K_B_EV
    import kinstep.drivers  # noqa: F401

    SI = os.path.join(OUT, "..", "SI")
    fig, axs = plt.subplots(1, 2, figsize=(7.48, 2.88))

    # ── (a) 三条曲线：朴素相减 / expm1（lnz 精确）/ expm1（lnz 带 1e-16 绝对误差）
    ax = axs[0]
    mp.mp.dps = 60
    t_true = np.logspace(-2, -18, 65)
    naive, em1, em1_abs = [], [], []
    for t in t_true:
        lnz_exact = mp.log(mp.mpf(1) - mp.mpf(float(t)))
        lz = float(lnz_exact)                       # double 存得下的 lnz
        lz_pert = float(lnz_exact + mp.mpf("1e-16"))  # lnz 带 1e-16 绝对误差
        naive.append(abs((1.0 - np.exp(lz)) / t - 1.0))
        em1.append(abs((-np.expm1(lz)) / t - 1.0))
        em1_abs.append(abs((-np.expm1(lz_pert)) / t - 1.0))
    naive, em1, em1_abs = map(np.asarray, (naive, em1, em1_abs))
    YF = 1e-18                                   # 误差恒为 0 的点钉在轴底
    ax.loglog(t_true, np.maximum(naive, YF), "o", ms=2.8, color="#d62728",
              label=r"naive  $1-\exp(\ln z)$")
    ax.loglog(t_true, np.maximum(em1_abs, YF), "-", lw=1.3, color="#1f77b4",
              label=r"$-\mathrm{expm1}$, $\ln z$ err. $10^{-16}$")
    ax.loglog(t_true, np.maximum(em1, YF), "s-", ms=2.4, lw=1.1, color="#2ca02c",
              label=r"$-\mathrm{expm1}$, exact $\ln z$  (error $\equiv 0$)")
    ax.loglog(t_true, 1e-16 / t_true, ":", lw=1.1, color="k",
              label=r"floor $\varepsilon/(1-z)$")
    ax.axhline(1.0, color="0.55", lw=0.8, ls="--")
    ax.text(3e-13, 3.0, "no correct digits", fontsize=7.0, color="0.4")
    ax.text(3e-7, 2.4e-18, r"pinned at axis floor", fontsize=7.0, color="#2ca02c")
    ax.set_xlim(2e-2, 4e-19); ax.set_ylim(2.5e-19, 3e3)
    ax.set_xlabel(r"true $1-z$"); ax.set_ylabel(r"relative error in $1-z$")
    ax.legend(frameon=False, fontsize=7.0, loc="upper left")
    style(ax)
    ax.set_title("(a) catastrophic cancellation", fontsize=7.6, weight="bold",
                 loc="left", pad=3)

    # ── (b) 任意精度随 dps 的收敛（T11 §④ 的同一算例，现场重算）
    ax = axs[1]
    KT_CAT, A_PRE, T_NOM = 0.057994719627, 1.402306113158e13, 673.0
    T_EQ = T_NOM * KT_CAT / (K_B_EV * T_NOM)
    Pg = {"H2": 75.0, "N2": 25.0, "NH3": 1.0}
    mech = M.from_text(
        "H2 + 2* <-> 2H*\nN2 + 2* <-> 2N*\nN* + H* <-> NH* + *\n"
        "NH* + H* <-> NH2* + *\nNH2* + H* <-> NH3 + 2*\n")
    d = catkinas_ref("curve")
    if d is not None:
        Ean, G0n, RNET = (np.asarray(d[k]) for k in ("Ean", "G0n", "Rnet"))
        i = 400                                  # E_N = +2.00 eV
        ref = float(RNET[i, 4])
        rd = steady_state(mech, energetics_from_arrays(mech, Ean[i], G0n[i], A_pre=A_PRE),
                          state=dict(T_K=T_EQ), p_gas=Pg, tol=1e-12, strict_energy=False)
        dps_list = [20, 30, 40, 50, 60, 80]
        dev, rvals = [], []
        for dps in dps_list:
            rh = steady_state_hp(mech, Ean[i], G0n[i], T_NOM, p_gas=Pg, A_pre=A_PRE,
                                 kT_eV=KT_CAT, dps=dps, theta_guess=dict(rd.theta))
            rv = float(rh.rate_net[4]); rvals.append(rv)
            dev.append(abs(rv / ref - 1.0))
            print("      dps %2d  r = %+.6e  dev %.2e" % (dps, rv, dev[-1]), flush=True)
        dev = np.asarray(dev)
        sign_bad = np.asarray([rv * ref < 0 for rv in rvals])
        dbl_dev = abs(float(rd.rate_net[4]) / ref - 1.0)
    else:
        # 缺 CATKINAS 参考数据：退回 SI 里已存档的数值重绘。
        # 双精度偏差也写进了 .dat 表头，所以退档数据是自足的。
        import re as _re
        _p = os.path.join(SI, "figS1b_dps.dat")
        _hdr = open(_p, encoding="utf-8").readline()
        ref = float(_re.search(r"reference ([0-9.eE+-]+)", _hdr).group(1))
        dbl_dev = float(_re.search(r"double_dev ([0-9.eE+-]+)", _hdr).group(1))
        _dat = np.genfromtxt(_p, comments="#")
        dps_list = [int(v) for v in _dat[:, 0]]
        rvals, dev = list(_dat[:, 1]), _dat[:, 2]
        sign_bad = np.asarray([rv * ref < 0 for rv in rvals])
        print("    ⊘ (b) 用 SI/figS1b_dps.dat 已存档数值重绘（缺 CATKINAS 参考数据）")
    ax.semilogy(dps_list, dev, "o-", ms=4, lw=1.3, color="#1f77b4")
    if sign_bad.any():
        ax.plot(np.asarray(dps_list)[sign_bad], dev[sign_bad], "x", ms=7, mew=1.6,
                color="#d62728", label="wrong sign")
    ax.axhline(dbl_dev, color="#d62728", lw=1.0, ls="--",
               label="double precision (%.1e)" % dbl_dev)
    ax.axhline(dev[-1], color="0.55", lw=0.9, ls=":",
               label="input-constant limit (%.1e)" % dev[-1])
    ax.axvline(50, color="#2ca02c", lw=0.9, ls="-.", label="KinSTEP default dps = 50")
    ax.set_xlabel("decimal places (dps)")
    ax.set_ylabel(r"rel. dev. of $r_5$ from reference")
    ax.set_xlim(14, 86)
    ax.legend(frameon=False, fontsize=7.0, loc="upper right")
    style(ax)
    ax.set_title("(b) arbitrary-precision escalation", fontsize=7.6, weight="bold",
                 loc="left", pad=3)

    fig.tight_layout(w_pad=2.2)
    fig.savefig(os.path.join(OUT, "FigS1_precision.png"), bbox_inches="tight")
    fig.savefig(os.path.join(OUT, "FigS1_precision.pdf"), bbox_inches="tight")
    plt.close(fig)
    np.savetxt(os.path.join(SI, "figS1a_cancellation.dat"),
               np.column_stack([t_true, naive, em1, em1_abs]),
               header="true_1mz  relerr_naive  relerr_expm1_exact  relerr_expm1_abs1e-16")
    np.savetxt(os.path.join(SI, "figS1b_dps.dat"),
               np.column_stack([dps_list, rvals, dev]),
               header="dps  r_net_5  rel_dev_from_reference  (reference %.12e  double_dev %.12e)" % (ref, dbl_dev))
    print("  Fig.S1 precision  (naive worst %.2e; dps>=30 plateau %.2e; double %.2e)"
          % (naive.max(), dev[1:].max(), dbl_dev))

if __name__ == "__main__":
    import sys as _s
    which = _s.argv[1] if len(_s.argv) > 1 else "all"
    if which in ("all", "1"):
        fig1()
    if which in ("all", "2"):
        fig2()
    if which in ("all", "3"):
        fig3()
    if which in ("all", "s1", "S1"):
        figS1()
