# -*- coding: utf-8 -*-
r"""Fig.5 ｜ 命令行界面：四条命令走完一次完整建模

## 为什么这张图值得占一个图位

审稿人读 CPC 的程序论文，第一个问题是「这东西好不好用」。
正文可以写「输入是纸上那样的反应式」，但**不如直接把屏幕给他看**。
这张图要说明的是四件事，每格一件：

  (a) 一条命令就拿到**机理已填好**的骨架，而且每个空槽都写着要算什么
  (b) `check` 先过四道结构闸，再被**来源闸拦住**并逐条说缺什么（退出码 3）
  (c) 一句自然语言 → 结构 + 问题，**零 LLM**；能量一律变成问题而不是假设
  (d) `run --drc` 一次给出覆盖度、精度边界、DRC 与 Campbell 求和自检

## ★ 设计决定：文本**现场真跑**，不写死

第一版想把输出粘成字符串常量 —— 否决了。理由：
程序一改，图就和程序脱节，而且**没有任何机制会发现**。
现在每次出图都用 `subprocess` 真跑一遍 `KINSTEP_LANG=en` 的命令再排版，
所以这张图永远是当下代码的真实输出。代价是出图要几十秒，值得。

## ★ 两个踩过的坑

1. **工作目录必须是纯 ASCII 路径**。CLI 会把文件路径回显出来，而本仓库路径里有汉字
   （`…/工具/微动力学/…`），图上就会出现中文 —— 和 T16 第一版踩的是同一个坑。
   所以用 `tempfile.mkdtemp()`，并且命令里只用**裸文件名**。
2. **行宽**。终端输出有些行上百字符，直接排进双栏图会糊成一团。
   每格都显式给 `keep`（保留哪些行）和 `width`（截断宽度），**截断在图题里说明**，
   不装作那就是全部输出。

用法:
    python3 make_fig5_cli.py            # 出图
    python3 make_fig5_cli.py --dump     # 只打印抓到的原始输出，排版前先肉眼看一遍
"""
import os
import re
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
PKG_ROOT = os.path.dirname(HERE)                 # …/工具/微动力学
FIGDIR = os.path.join(HERE, "figures")
SIDIR = os.path.join(HERE, "SI")

# ── 配色：深色终端。不用纯黑，印刷时纯黑底会把细字吃掉 ──
C_BG = "#1b1f27"       # 面板底
C_BAR = "#2b313c"      # 标题栏
C_FG = "#d7dce5"       # 正文
C_DIM = "#8b93a1"      # 分隔线与次要信息
C_OK = "#6fcf7f"       # ✓
C_BAD = "#f2777a"      # ✗
C_WARN = "#f0c674"     # ⚠
C_STAR = "#8ab8f0"     # ★ 与被强调的数字
C_CMD = "#9fe0c4"      # 命令行
C_KEY = "#c9a7e8"      # 小节标题


def _runner():
    """建一个纯 ASCII 工作目录，返回 (run, workdir)。"""
    wd = tempfile.mkdtemp(prefix="kinstep_fig5_")
    env = dict(os.environ)
    env["KINSTEP_LANG"] = "en"
    env["PYTHONPATH"] = PKG_ROOT + os.pathsep + env.get("PYTHONPATH", "")
    env["COLUMNS"] = "110"

    def run(args):
        p = subprocess.run([sys.executable, "-m", "kinstep.cli"] + args,
                           cwd=wd, env=env, capture_output=True, text=True,
                           timeout=900)
        return p.returncode, (p.stdout + p.stderr)
    return run, wd


def capture():
    """真跑四条命令，返回 {键: (命令字符串, 退出码, 输出)}。"""
    run, wd = _runner()
    out = {}

    rc, s = run(["new", "--template", "her", "--force", "her.kinstep.toml"])
    # ★ (a) 格还要把**生成出来的文件**贴上去 —— 这正是坑 29 修的那份。
    #   只截屏回显是不够的：读者真正要动手改的是这个文件。
    skel = open(os.path.join(wd, "her.kinstep.toml"), encoding="utf-8").read()
    out["a"] = ("kinstep new --template her her.kinstep.toml", rc,
                s + "\n@@@CAT her.kinstep.toml\n" + skel)

    rc, s = run(["check", "her.kinstep.toml"])
    out["b"] = ("kinstep check her.kinstep.toml", rc, s)

    rc, s = run(["ask", "hydrogen evolution at pH 0, -0.1 V vs SHE, 25 C"])
    out["c"] = ('kinstep ask "hydrogen evolution at pH 0, -0.1 V vs SHE, 25 C"', rc, s)

    run(["new", "ammonia.kinstep.toml", "--force"])
    rc, s = run(["run", "ammonia.kinstep.toml", "--drc", "--no-plot"])
    out["d"] = ("kinstep run ammonia.kinstep.toml --drc", rc, s)

    # ★ 这张图的硬判据：英文模式下一个汉字都不许有。
    #   图是要印在论文上的，这里再查一遍，不依赖 T16 跑过。
    cjk = re.compile(r"[一-鿿　-〿＀-￯]")
    for k, (cmd, rc_, s_) in out.items():
        hit = cjk.findall(s_)
        if hit:
            raise SystemExit("✗ (%s) 输出里有 %d 个 CJK 字符：%s —— 先修翻译再出图"
                             % (k, len(hit), "".join(sorted(set(hit))[:12])))
    return out, wd


# ── 每格保留哪些行。用**子串匹配**而不是行号：行号会随输出微调而错位 ──
PANELS = {
    "a": dict(
        title="(a)  one command, then every blank says what to compute",
        keep=["wrote her.kinstep.toml", "energy slots are left",
              "Tafel slope identifies", "@@@CAT",
              "# her - hydrogen evolution", "# ! This template carries NO energy",
              "[mechanism]", "H+ + e- + * <-> H*", "H* + H+ + e- -> H2 + *",
              "2H* -> H2 + 2*", "[conditions]", "T_K = 298.15",
              "# (1) H+ + e- + * <-> H*", "#     what to compute: Volmer step",
              "# Ga_fwd = ?", "# dG_rxn = ?"],
        width=82),
    "b": dict(
        title="(b)  structural gates pass, the provenance gate refuses  (exit code 3)",
        keep=["mechanism: 3 steps", "(1) H+", "(2) H*", "(3) 2H*",
              "all four gates passed", "temperature T_K",
              "provenance gate rejected", "value(s) missing",
              "step 1 / Ga_fwd", "step 1 / dG_rxn"],
        width=82, tail="  (steps 2 and 3 list the same two slots; 6 in total)"),
    "c": dict(
        title="(c)  one sentence in, structure and questions out  (rule-based, no LLM)",
        keep=["matched system:", "values extracted from your sentence",
              "U_elec =", "U_ref =", "pH_val =", "still needed from you",
              "the energies of all 3 steps", "(1) H+ + e-", "temperature (K or degC)",
              "active-site areal density"],
        width=82),
    "d": dict(
        title="(d)  solve, then report coverages, precision limits, DRC and self-checks",
        keep=["steady state (ode", "coverages:", "vacancies:",
              "scaled residual", "step 1  1-z=", "flux-conservation self-check",
              "step 1  H2 + 2*", "step 2  N2 + 2*", "step 5  NH2*",
              "sum rule  sum X_B = 1.000000000 (", "reaction orders:",
              "apparent activation energy"],
        width=82),
}


def select(text, keep, width, tail=None):
    """按子串挑行、去 ANSI、截断。缺行就报错 —— 不要默默少画一行。"""
    lines = [re.sub(r"\x1b\[[0-9;]*m", "", ln).rstrip() for ln in text.splitlines()]
    picked, missing = [], []
    for needle in keep:
        hit = next((ln for ln in lines if needle in ln), None)
        if hit is None:
            missing.append(needle)
        else:
            picked.append(hit)
    if missing:
        raise SystemExit("✗ 这些行没抓到（输出变了？）：%s" % missing)
    out = []
    for ln in picked:
        # `@@@CAT x` 是 capture() 插的分隔标记，排版时还原成一条命令行
        if ln.startswith("@@@CAT "):
            ln = "$ cat " + ln[len("@@@CAT "):]
        out.append(ln if len(ln) <= width else ln[:width - 1] + "…")
    if tail:
        out.append(tail)
    return out


def colour_of(line):
    """按行首标记上色。只认程序真会打印的那几个标记。"""
    s = line.strip()
    if s.startswith("$ "):
        return C_CMD
    if s.startswith("✓"):
        return C_OK
    if s.startswith("✗"):
        return C_BAD
    if s.startswith("⚠"):
        return C_WARN
    if s.startswith("★") or s.startswith("·"):
        return C_STAR
    if s.startswith("["):
        return C_KEY
    if s.startswith("#"):
        return C_DIM
    return C_FG


def draw(cap, out_png, out_pdf):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyBboxPatch, Rectangle

    plt.rcParams.update({"font.family": "monospace", "pdf.fonttype": 42,
                         "ps.fonttype": 42})

    rows = {k: select(cap[k][2], **{kk: vv for kk, vv in PANELS[k].items()
                                    if kk != "title"})
            for k in "abcd"}
    # ★ 版面为什么是单列四格，而不是原来的 2×2：
    #   Elsevier 要求图内文字终稿不小于 7 pt。终端行宽 82 字符，等宽字体在 7 pt 下
    #   需要 82×7×0.6023 ≈ 346 pt ≈ 4.8 in 的文字宽度；2×2 时每格只有 3.7 in，
    #   **几何上就放不下**，原图只好用 4.5 pt —— 那是不可读的。
    #   改成单列后每格可用宽度 ≈ 5.4 in，7 pt 放得下，图宽取 Elsevier 的 1.5 栏（140 mm）。
    FS = 7.0                      # 终端正文字号（= 下限，不再小）
    LINE_IN = FS * 1.35 / 72.0    # 行距（英寸）
    PAD_LINES = 3.4               # 每格的标题栏 + 上下留白，按"行"折算
    nline = {k: len(rows[k]) for k in "abcd"}
    H = (sum(nline.values()) + 4 * PAD_LINES) * LINE_IN + 0.10

    fig = plt.figure(figsize=(5.51, H), dpi=400)
    gs = fig.add_gridspec(4, 1, left=0.012, right=0.988, top=0.985, bottom=0.015,
                          hspace=0.09,
                          height_ratios=[nline[k] + PAD_LINES for k in "abcd"])

    for idx, k in enumerate("abcd"):
        ax = fig.add_subplot(gs[idx, 0])
        nmax = nline[k]
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis("off")
        # 面板底 + 标题栏
        ax.add_patch(FancyBboxPatch((0.004, 0.004), 0.992, 0.992,
                                    boxstyle="round,pad=0,rounding_size=0.012",
                                    transform=ax.transAxes, facecolor=C_BG,
                                    edgecolor="#3a4250", linewidth=0.6, zorder=0))
        ax.add_patch(Rectangle((0.004, 0.905), 0.992, 0.091, transform=ax.transAxes,
                               facecolor=C_BAR, edgecolor="none", zorder=1))
        # 三个窗口圆点，纯装饰，告诉读者"这是一个终端窗口"
        for j, col in enumerate(("#e06b6b", "#e0c06b", "#7fc47f")):
            ax.add_patch(plt.Circle((0.022 + j * 0.019, 0.950), 0.0065,
                                    transform=ax.transAxes, facecolor=col,
                                    edgecolor="none", zorder=2))
        ax.text(0.082, 0.949, PANELS[k]["title"], transform=ax.transAxes,
                va="center", ha="left", fontsize=FS + 0.5, color="#e7ecf4",
                family="sans-serif", fontweight="bold", zorder=3)

        # ★ parse_math=False 不能省：matplotlib 默认把 `$…$` 当数学公式，
        #   于是 "$ echo $?" 被吃成斜体的 "echo?" —— 第一版图上就是这么错的。
        y, dy = 0.872, 0.868 / (nmax + 2.6)
        ax.text(0.022, y, "$ " + cap[k][0], transform=ax.transAxes, va="top",
                ha="left", fontsize=FS, color=C_CMD, zorder=3, parse_math=False)
        y -= dy * 1.35
        for ln in rows[k]:
            ax.text(0.022, y, ln, transform=ax.transAxes, va="top", ha="left",
                    fontsize=FS, color=colour_of(ln), zorder=3, parse_math=False)
            y -= dy
        ax.text(0.022, y - dy * 0.1, "exit code %d" % cap[k][1],
                transform=ax.transAxes, va="top", ha="left", fontsize=FS,
                color=C_DIM, zorder=3, parse_math=False)

    os.makedirs(FIGDIR, exist_ok=True)
    fig.savefig(out_png, facecolor="white")
    fig.savefig(out_pdf, facecolor="white")
    plt.close(fig)


def main():
    cap, wd = capture()
    if "--dump" in sys.argv:
        for k in "abcd":
            cmd, rc, s = cap[k]
            print("\n" + "=" * 78 + "\n$ %s   [exit %d]\n" % (cmd, rc) + "=" * 78)
            print(s)
        return 0
    png = os.path.join(FIGDIR, "Fig5_cli.png")
    pdf = os.path.join(FIGDIR, "Fig5_cli.pdf")
    draw(cap, png, pdf)
    # 原始输出一并存进 SI：图上是截过的，审稿人要看全的就在这里
    os.makedirs(SIDIR, exist_ok=True)
    dump = os.path.join(SIDIR, "fig5_cli_session.txt")
    with open(dump, "w", encoding="utf-8") as f:
        f.write("# Fig.5 raw terminal session, captured by make_fig5_cli.py\n"
                "# KINSTEP_LANG=en, run in a clean temporary directory.\n"
                "# The figure shows a selection of these lines; this file is complete.\n")
        for k in "abcd":
            cmd, rc, s = cap[k]
            f.write("\n" + "=" * 78 + "\n$ %s\n" % cmd + "=" * 78 + "\n")
            f.write(s)
            f.write("exit code: %d\n" % rc)
    for p in (png, pdf, dump):
        print("  %-28s %7.1f KB" % (os.path.basename(p), os.path.getsize(p) / 1024))
    print("  workdir: %s" % wd)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
