# -*- coding: utf-8 -*-
"""把 manuscript.md / SI.md 编译成 Word（CPC 投稿格式）。

CPC 的要求与一般期刊不同的地方：
  · **PROGRAM SUMMARY 是强制项**，而且是一张双列表 —— 这里做成真 Word 表格
  · 正文要连续行号（审稿要求）
  · 图随文插入并带编号题注
Markdown 里的 `**粗体**`、`$公式$`、表格、列表都会被正确转换。
"""
import os, re, sys
from PIL import Image          # 只用来读图的像素宽高，决定按宽还是按高插
from docx import Document
from docx.shared import Pt, Inches, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_LINE_SPACING
from docx.enum.section import WD_SECTION
from docx.oxml.ns import qn
from docx.oxml import OxmlElement


SIFIGCAP = {
 "FigS1_precision.png":
   ("Fig. S1.", "Precision of the reversibility factor $1-z$. "
    "(a) Relative error in $1-z$ at double precision. Evaluating $1-\\exp(\\ln z)$ directly "
    "(red) loses every correct digit once $1-z$ falls below $\\sim10^{-16}$: the subtraction "
    "cancels. With $-\\mathrm{expm1}(\\ln z)$ and an exactly representable $\\ln z$ (green) "
    "the error is identically zero at every point shown, so the curve is pinned at the axis "
    "floor. The practically attainable accuracy (blue) is set instead by the absolute accuracy "
    "of $\\ln z$ itself: because $\\ln z = \\Delta G/k_BT + \\sum_j \\nu_j \\ln a_j$ is a "
    "near-cancelling sum of terms of order unity, that accuracy is absolute rather than "
    "relative, and the relative error in $1-z$ therefore follows the floor "
    "$\\varepsilon/(1-z)$ (dotted). This floor, not the subtraction, is what forces the "
    "escalation to arbitrary precision. "
    "(b) Convergence of the arbitrary-precision solution with working precision, at the "
    "strong-binding point $E_\\mathrm{N}=+2.00$ eV of the CATKINAS scan where the reference "
    "net rate is $4.801095\\times10^{-32}$ s$^{-1}$ site$^{-1}$. Double precision is wrong by "
    "$6.7\\times10^{-3}$ and 20 digits give the wrong sign; the result plateaus from 30 digits "
    "onward at $\\approx3\\times10^{-10}$, which is the precision of the reverse-engineered "
    "input "
    "constants and not of the arithmetic. The KinSTEP default of 50 digits leaves margin."),
}

HERE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.join(HERE, "figures")

FIGCAP = {
 "Fig1_framework.png":
   ("Fig. 1.", "Architecture of KinSTEP. Input passes through a fixed set of gates before "
    "reaching a modality-agnostic core; external stimuli enter only through the driver plug-in "
    "layer. The solver is a four-level ladder, each level of which was added in response to a "
    "measured failure. Five consistency relations are evaluated at no extra cost."),
 "Fig2_verification.png":
   ("Fig. 2.", "Verification against the CATKINAS 401-point descriptor scan. "
    "(a) Volcano curve; points at which the solver escalated automatically to arbitrary precision "
    "are shown in red. (b) Point-by-point relative deviation in the net rate: in double precision "
    "the deviation rises above $10^{-6}$ for 27 points on the strong-binding branch, whereas with "
    "automatic escalation it stays below $10^{-9}$ throughout. (c) Summary of all verification "
    "targets; longer bars denote closer agreement."),
 "Fig3_cases.png":
   ("Fig. 3.", "One example per modality, all computed with the same core. "
    "(a) Hydrogen evolution: the Tafel slope converges to 118.3 mV dec$^{-1}$, exactly "
    "$\\ln 10\\,k_BT/(\\beta e)$, although this value is nowhere supplied as input. "
    "(b) Photocatalytic water splitting: the local power-law order in photon flux moves from 1 to "
    "1/2 as recombination takes over, the apparent quantum yield falling by five orders of "
    "magnitude across the same range. (c) Cavitation: the correctly weighted rate "
    "$\\sum_i w_i r_i$ and the rate evaluated at the mean temperature differ by up to 23 orders of "
    "magnitude; the vertical line marks the analytically predicted crossover fraction. "
    "(d) Step degrees of rate control along the ammonia-synthesis descriptor. Rate control "
    "passes from the final hydrogenation step (5), which governs the strong-binding branch, to "
    "dissociative H$_2$ adsorption (2) on the weak-binding branch, the crossover falling between "
    "$E_\\mathrm{N}=-0.68$ and $-0.56$ eV. On the strong-binding side $X_{B,2}$ is negative (down "
    "to $-0.36$): lowering that barrier slows the overall reaction. The Campbell sum rule "
    "$\\sum_i X_{B,i}=1$ is reproduced to $7.7\\times10^{-7}$ at every point shown, and is used "
    "here as the acceptance criterion for the perturbed solutions rather than as a "
    "post-hoc check."),
 "Fig4_maps.png":
   ("Fig. 4.", "Two-dimensional operating maps, one row per modality, 4992 converged "
    "solutions obtained with the same core. Rows are labelled (a)–(d) and the three panels of "
    "each row (a1)–(a3), (b1)–(b3) and so on; a reference such as Fig. 4a in the text means the "
    "whole modality row. "
    "(a) Ammonia synthesis over $T=500$–900 K and $P=1$–200 bar (1728 points): signed net rate "
    "on an $\\mathrm{asinh}$ scale, with the $r=0$ equilibrium contour drawn explicitly rather "
    "than left to the colour scale; the net rate is positive at 60.2\\% of the points and peaks "
    "at $0.206$ s$^{-1}$ site$^{-1}$ at 857 K and 200 bar (a1). Coverage panels (a2, a3) show "
    "$\\theta_{\\mathrm{H}^*}$ and the vacancy fraction trading against each other. "
    "(b) Hydrogen evolution over $U=-0.40$ to $-0.02$ V vs SHE and pH 0–14: the current density "
    "spans thirty orders of magnitude (b1), and the Tafel-slope panel (b3) is bimodal — of the 1394 points "
    "at which the slope is defined and below 400 mV dec$^{-1}$, 49.2\\% lie in the "
    "25–45 mV dec$^{-1}$ band and 34.1\\% in the 100–140 mV dec$^{-1}$ band, so the "
    "rate-determining step changes across the operating plane. "
    "(c) Photocatalysis over twelve decades of photon flux and four decades of recombination "
    "coefficient: the local power-law order in flux (c3), computed as a local logarithmic derivative "
    "rather than assumed, spans 0.5000 to 0.9999, averaging 0.932 over the three lowest-flux "
    "columns and 0.500 over the three highest, so the whole first-order-to-half-order transition "
    "is recovered. "
    "(d) Cavitation over hot-spot fraction $10^{-5}$–$10^{-2}$ and hot-spot temperature "
    "1500–6000 K: the correctly weighted rate and the rate evaluated at the mean condition (d2) "
    "differ by 16.1 to 23.9 orders of magnitude, and at 3 of the 721 points where the comparison "
    "is defined the ratio is negative, i.e. the mean-condition estimate reverses the direction of "
    "reaction. Hatched cells mark quantities that are undefined rather than unconverged: the "
    "current density where it underflows to zero so that its logarithm does not exist, the Tafel "
    "slope where it is not defined or exceeds the plotted bound, and the rate ratio where the "
    "mean-condition rate underflows to zero. Every one of the 4992 steady states converged. Note "
    "that the pale grey at the centre of the diverging scale in the first panel is a rate close "
    "to zero, not missing data."),
 "Fig5_cli.png":
   ("Fig. 5.", "The command-line interface, reproduced from real terminal output. The figure is "
    "not a transcript pasted into the manuscript: the plotting script re-runs the four commands "
    "with $\\mathtt{KINSTEP\\_LANG{=}en}$ in a clean temporary directory and typesets whatever "
    "they print, so it cannot drift away from the code. "
    "(a) One command writes a skeleton with the mechanism already filled in; every blank energy "
    "slot carries a sentence stating what to compute, and the generated file is shown beneath the "
    "command. (b) The structural gates pass, after which the provenance gate refuses the model "
    "because energies are missing, listing each slot and exiting with code 3. (c) One sentence of "
    "plain English becomes a structure plus a list of questions, by rule and with no language "
    "model: stated quantities are extracted and tagged, whereas every energy becomes a question "
    "rather than an assumption. (d) The solver reports coverages, the attainable precision of each "
    "step, degrees of rate control and the Campbell sum rule in a single pass. Long lines are "
    "truncated to fit the column; the complete session is deposited as "
    "$\\mathtt{fig5\\_cli\\_session.txt}$."),
}

PROGRAM_SUMMARY = [
 ("Program Title", "KinSTEP"),
 ("CPC Library link to program files", "(to be added by Technical Editor)"),
 ("Developer's repository link",
  "https://github.com/xinlei-chem/kinstep (version v0.1.0 archived at "
  "https://doi.org/10.5281/zenodo.23613004)"),
 ("Code Ocean capsule", "(to be added by Technical Editor)"),
 ("Licensing provisions", "GPL-3.0-or-later"),
 ("Programming language", "Python 3.9+ (tested on CPython 3.12)"),
 ("Supplementary material",
  "Supplementary Information: verification tables, constant conventions, benchmark data, "
  "full self-check output"),
 ("External routines/libraries",
  "NumPy and SciPy (required); mpmath and matplotlib (optional extras, used for "
  "arbitrary-precision escalation and for plotting respectively). No compiled extension "
  "and no build step."),
 ("Typical running time",
  "A single steady state solves in well under a second. On one core of a desktop CPU the "
  "401-point descriptor scan of Section 4 takes 4.2 s (1.7 s with arbitrary precision "
  "disabled), the 315-point flow-reactor map 7.8 s, and the 4992-point four-modality map "
  "of Fig. 4 about 20 min. Peak resident memory stays below 100 MB in all of these, "
  "dominated by the NumPy/SciPy import footprint rather than by the model."),
 ("Nature of problem",
  "Steady-state and transient microkinetic modelling of heterogeneous catalysis under different "
  "external stimuli (temperature, electrode potential, illumination, ultrasound). The governing "
  "equations are stiff, the unknowns span tens of orders of magnitude, and quasi-equilibrated "
  "elementary steps cause catastrophic cancellation in the net rate. Existing packages address "
  "one modality at a time and require modality-specific input scripts."),
 ("Solution method",
  "Elementary steps are parsed from plain text and compiled to a stoichiometric matrix; electron "
  "stoichiometry is counted automatically. Coverages are solved in ln(theta) space with an "
  "analytic Jacobian using a four-level ladder: stiff ODE continuation for a physical initial "
  "guess, damped Newton with step limiting and backtracking, weighted projection of the net rates "
  "onto the flux-conservation manifold, and an automatically triggered arbitrary-precision stage. "
  "External stimuli enter through a driver plug-in returning perturbed energies, or through a "
  "second hook returning weighted sub-states. Degrees of rate control, reaction orders, apparent "
  "activation energies and Tafel slopes are obtained by central differences with Richardson "
  "extrapolation, the step size being selected from the Campbell sum rule."),
 ("Additional comments",
  "Core dependencies are NumPy and SciPy only; mpmath and matplotlib are optional extras. A "
  "command-line interface provides input validation, plain-language explanation of a model, and a "
  "rule-based natural-language front end that produces structure and questions but never "
  "numerical energies. The deposit carries a regression suite of 82 assertions, all of which run "
  "on a clean download: 77 pass with NumPy and SciPy alone and five are skipped, with their "
  "reasons, because they need an optional extra (four mpmath, one matplotlib); with both extras "
  "installed all 82 pass. The reference values compared against in Section 4 are archived with "
  "the deposit, so no third-party software is needed to reproduce any assertion or figure. The "
  "interface is bilingual: all user-facing output is available in English via KINSTEP_LANG=en, "
  "and one regression asserts that no CJK character appears in English mode."),
]


def set_font(run, name="Times New Roman", size=11, bold=False, italic=False, mono=False):
    run.font.name = "Consolas" if mono else name
    run.font.size = Pt(size)
    run.bold = bold
    run.italic = italic
    r = run._element.rPr.rFonts
    r.set(qn("w:eastAsia"), "Times New Roman")


def add_line_numbers(section):
    sectPr = section._sectPr
    ln = OxmlElement("w:lnNumType")
    ln.set(qn("w:countBy"), "1"); ln.set(qn("w:restart"), "continuous")
    ln.set(qn("w:distance"), "360")
    sectPr.append(ln)


_BLOCK = re.compile(r"^(#|```|\* |- |\[\d+\]|\||---$|\$\$)")
INLINE = re.compile(r"(\*\*[^*]+\*\*|\*[^*]+\*|`[^`]+`|\$[^$]+\$)")


# ── LaTeX → Unicode：直接删 \cmd 会把 \nu、\Delta 整个吞掉，
#    正文里就会出现「-, E_qfl」这种读不懂的东西（已在 PDF 里亲眼看到）
MATHMAP = {
    r"\alpha": "α", r"\beta": "β", r"\gamma": "γ", r"\Gamma": "Γ",
    r"\delta": "δ", r"\Delta": "Δ", r"\epsilon": "ε", r"\varepsilon": "ε",
    r"\zeta": "ζ", r"\eta": "η", r"\theta": "θ", r"\Theta": "Θ",
    r"\kappa": "κ", r"\lambda": "λ", r"\Lambda": "Λ", r"\mu": "μ",
    r"\nu": "ν", r"\xi": "ξ", r"\pi": "π", r"\Pi": "Π", r"\rho": "ρ",
    r"\sigma": "σ", r"\Sigma": "Σ", r"\tau": "τ", r"\phi": "φ",
    r"\Phi": "Φ", r"\varphi": "φ", r"\chi": "χ", r"\psi": "ψ",
    r"\omega": "ω", r"\Omega": "Ω",
    r"\times": "×", r"\cdot": "·", r"\approx": "≈", r"\simeq": "≃",
    r"\sim": "~", r"\pm": "±", r"\mp": "∓", r"\leq": "≤", r"\le": "≤",
    r"\geq": "≥", r"\ge": "≥", r"\ll": "≪", r"\gg": "≫", r"\neq": "≠",
    r"\equiv": "≡", r"\propto": "∝", r"\infty": "∞", r"\partial": "∂",
    r"\nabla": "∇", r"\int": "∫", r"\sum": "Σ", r"\prod": "Π",
    r"\sqrt": "√", r"\rightarrow": "→", r"\leftarrow": "←",
    r"\leftrightarrow": "⇌", r"\Rightarrow": "⇒", r"\to": "→",
    r"\langle": "⟨", r"\rangle": "⟩", r"\lvert": "|", r"\rvert": "|",
    r"\mathrm": "", r"\mathbf": "", r"\mathit": "", r"\text": "",
    # 函数名要留住后面的空格（\ln z → "ln z"），符号则要吃掉（\Delta E → "ΔE"）
    r"\ln": "ln ", r"\log": "log ", r"\exp": "exp ",
    r"\min": "min ", r"\max": "max ", r"\lim": "lim ",
    r"\det": "det ", r"\sin": "sin ", r"\cos": "cos ",
    r"\left": "", r"\right": "", r"\!": "", r"\,": " ", r"\;": " ",
    r"\ ": " ", r"\quad": "  ", r"\qquad": "    ",
    r"\ddagger": "‡", r"\circ": "°", r"\dots": "…",
    r"\ldots": "…", r"\cdots": "⋯", r"\ne": "≠",
    r"\longrightarrow": "⟶", r"\longleftarrow": "⟵",
    r"\mathcal": "", r"\mathsf": "", r"\big": "", r"\Big": "",
    r"\bigg": "", r"\Bigg": "", r"\displaystyle": "",
}
_UP_A, _UP_B = "\x01", "\x02"     # \mathrm{...} 的内容用哨兵包起来，排正体
_SUPSUB = re.compile(r"([_^])\{([^{}]*)\}"
                     r"|([_^])(\x01[^\x02]*\x02|\\?[A-Za-z0-9+\-−]|\S)")
_LATIN = re.compile(r"([A-Za-z]+)")
# 函数名按期刊惯例排正体，只有变量才斜体
_UPRIGHT = {"ln", "log", "exp", "min", "max", "lim", "det", "sin", "cos", "tan",
            "d", "e"}


def _grab_group(s, i):
    """s[i] 必须是 '{'；返回 (组内文本, 右花括号之后的位置)。"""
    d = 0
    for j in range(i, len(s)):
        if s[j] == "{":
            d += 1
        elif s[j] == "}":
            d -= 1
            if d == 0:
                return s[i + 1:j], j + 1
    return s[i + 1:], len(s)


def _expand_frac(tex):
    r"""\frac{a}{b} → (a)/(b)。

    不能用正则：分母里有嵌套花括号（\frac{e^{x}-1}{e^{fx}+e^{(1-f)x}}），
    `[^{}]*` 接不住。而直接把 \frac 删掉会把分子分母并成一串，意思就变了。
    """
    out, i = [], 0
    while True:
        k = tex.find("\\frac", i)
        if k < 0:
            out.append(tex[i:])
            return "".join(out)
        out.append(tex[i:k])
        j = k + 5
        while j < len(tex) and tex[j] == " ":
            j += 1
        if j < len(tex) and tex[j] == "{":
            num, j = _grab_group(tex, j)
            while j < len(tex) and tex[j] == " ":
                j += 1
            if j < len(tex) and tex[j] == "{":
                den, j = _grab_group(tex, j)
                out.append("(%s)/(%s)" % (_expand_frac(num), _expand_frac(den)))
                i = j
                continue
        out.append(tex[k:k + 5])        # 不是标准写法，原样留着
        i = k + 5


def _plain(par, seg, size, sup=False, sub=False):
    """拉丁字母用斜体（变量），数字与符号用正体 —— 期刊排版的通行约定。"""
    seg = seg.replace("{", "").replace("}", "")
    seg = re.sub(r"\\[a-zA-Z]+", "", seg).replace("\\", "")
    for chunk in re.split(r"(\x01[^\x02]*\x02)", seg):
        if not chunk:
            continue
        upright = chunk.startswith(_UP_A)
        parts = [chunk[1:-1]] if upright else [q for q in _LATIN.split(chunk) if q]
        for part in parts:
            if not part:
                continue
            r = par.add_run(part)
            set_font(r, size=size, italic=(not upright)
                     and bool(_LATIN.fullmatch(part)) and part not in _UPRIGHT)
            if sup:
                r.font.superscript = True
            if sub:
                r.font.subscript = True


def math_runs(par, tex, size=11):
    """把 $...$ 渲染成真正的上/下标 run，LaTeX 命令转 Unicode。"""
    # ★ LaTeX 里 \cmd 后面的那个空格是命令终结符，不是真空格：
    #   直接 replace 会得到「Δ E_qfl」「⟨ r(T)⟩」。按长度倒序替换并吃掉终结空格；
    #   负向前瞻防止 \sim 先吃掉 \simeq、\le 先吃掉 \leftrightarrow。
    tex = _expand_frac(tex)
    tex = re.sub(r"\\(?:mathrm|mathbf|mathit|text|operatorname)\{([^{}]*)\}",
                 lambda m: _UP_A + m.group(1) + _UP_B, tex)
    for k in sorted(MATHMAP, key=len, reverse=True):
        v = MATHMAP[k]
        if k[1:].isalpha():
            tex = re.sub(re.escape(k) + r"(?![A-Za-z]) ?", v.replace("\\", "\\\\"), tex)
        else:
            tex = tex.replace(k, v)
    tex = tex.replace("\u2009(", "(")          # exp(x) 不要写成 exp (x)
    pos = 0
    for m in _SUPSUB.finditer(tex):
        if m.start() > pos:
            _plain(par, tex[pos:m.start()], size)
        mark = m.group(1) or m.group(3)
        body = m.group(2) if m.group(2) is not None else m.group(4)
        _plain(par, body, size, sup=(mark == "^"), sub=(mark == "_"))
        pos = m.end()
    if pos < len(tex):
        _plain(par, tex[pos:], size)


def emit(par, text, size=11, base_italic=False):
    """处理行内 **粗体** / *斜体* / `等宽` / $公式$。"""
    for tok in INLINE.split(text):
        if not tok:
            continue
        if tok.startswith("**") and tok.endswith("**"):
            set_font(par.add_run(tok[2:-2]), size=size, bold=True)
        elif tok.startswith("`") and tok.endswith("`"):
            set_font(par.add_run(tok[1:-1]), size=size - 0.5, mono=True)
        elif tok.startswith("$") and tok.endswith("$"):
            math_runs(par, tok[1:-1], size=size)
        elif tok.startswith("*") and tok.endswith("*") and len(tok) > 2:
            set_font(par.add_run(tok[1:-1]), size=size, italic=True)
        else:
            set_font(par.add_run(tok), size=size, italic=base_italic)


def md_table(doc, rows, size=9):
    ncol = max(len(r) for r in rows)
    t = doc.add_table(rows=0, cols=ncol)
    t.style = "Table Grid"
    for i, row in enumerate(rows):
        cells = t.add_row().cells
        for j in range(ncol):
            cells[j].text = ""
            p = cells[j].paragraphs[0]
            p.paragraph_format.space_before = Pt(1)
            p.paragraph_format.space_after = Pt(1)
            if j < len(row):
                emit(p, row[j], size=size)
                if i == 0:
                    for r in p.runs:
                        r.bold = True
    return t


def build(md_path, out_path, title_prefix="", with_summary=False, with_figs=False,
          figmap=None):
    doc = Document()
    sec = doc.sections[0]
    sec.top_margin = sec.bottom_margin = Inches(1.0)
    sec.left_margin = sec.right_margin = Inches(1.0)
    add_line_numbers(sec)
    st = doc.styles["Normal"]
    st.font.name = "Times New Roman"; st.font.size = Pt(11)
    st.paragraph_format.line_spacing_rule = WD_LINE_SPACING.DOUBLE
    st.paragraph_format.space_after = Pt(0)

    lines = open(md_path, encoding="utf-8").read().splitlines()

    def gather(k):
        """把 markdown 的软换行并回一段。

        不并的话，跨行的 **粗体** 会在每一行里各剩一半、两个星号都原样印出来
        （摘要下面那句「**photo and electro drivers compose … **」就是这么漏的）。
        """
        buf, j = [lines[k].strip()], k + 1
        while j < len(lines):
            t = lines[j].strip()
            if not t or _BLOCK.match(t):
                break
            buf.append(t)
            j += 1
        return " ".join(buf), j

    i, nfig, table_buf = 0, 0, []
    while i < len(lines):
        ln = lines[i]
        # 表格累积
        if ln.strip().startswith("|"):
            table_buf.append([c.strip() for c in ln.strip().strip("|").split("|")])
            i += 1
            continue
        if table_buf:
            rows = [r for r in table_buf
                    if not all(set(c) <= set("-: ") for c in r if c)]
            if rows:
                md_table(doc, rows)
                doc.add_paragraph()
            table_buf = []
        s = ln.strip()
        if not s:
            i += 1; continue
        if s == "---":
            i += 1; continue
        if s.startswith("#"):
            lvl = len(s) - len(s.lstrip("#"))
            txt = s.lstrip("# ").strip()
            p = doc.add_paragraph()
            p.paragraph_format.space_before = Pt(10 if lvl > 1 else 14)
            p.paragraph_format.space_after = Pt(4)
            emit(p, txt, size=14 if lvl == 1 else (12 if lvl == 2 else 11))
            for r in p.runs:
                r.bold = True
            if lvl == 1:
                p.alignment = WD_ALIGN_PARAGRAPH.CENTER
            # PROGRAM SUMMARY 用真表格
            if with_summary and txt.upper().startswith("PROGRAM SUMMARY"):
                md_table(doc, [["Item", "Content"]] + [list(x) for x in PROGRAM_SUMMARY], size=9)
                doc.add_paragraph()
                while i + 1 < len(lines) and not lines[i + 1].startswith("##"):
                    i += 1
            i += 1; continue
        if s.startswith("```"):
            i += 1
            code = []
            while i < len(lines) and not lines[i].strip().startswith("```"):
                code.append(lines[i]); i += 1
            i += 1
            for c in code:
                p = doc.add_paragraph()
                p.paragraph_format.line_spacing_rule = WD_LINE_SPACING.SINGLE
                p.paragraph_format.left_indent = Inches(0.3)
                set_font(p.add_run(c), size=9.5, mono=True)
            doc.add_paragraph()
            continue
        if s.startswith("$$"):            # 行间公式，可跨多行
            buf = []
            while True:
                buf.append(lines[i].strip())
                done = len(buf) > 1 or len(buf[0]) > 2
                if done and buf[-1].endswith("$$"):
                    break
                i += 1
                if i >= len(lines):
                    break
            i += 1
            tex = " ".join(buf).strip()
            tex = tex[2:] if tex.startswith("$$") else tex
            tex = tex[:-2] if tex.endswith("$$") else tex
            pm = doc.add_paragraph()
            pm.alignment = WD_ALIGN_PARAGRAPH.CENTER
            pm.paragraph_format.space_before = Pt(6)
            pm.paragraph_format.space_after = Pt(6)
            math_runs(pm, tex.strip(), size=11)
            continue
        if s.startswith("* ") or s.startswith("- "):
            txt, i = gather(i)
            p = doc.add_paragraph(style="List Bullet")
            p.paragraph_format.line_spacing_rule = WD_LINE_SPACING.DOUBLE
            emit(p, txt[2:])
            continue
        if re.match(r"^\[\d+\]", s):          # 参考文献
            txt, i = gather(i)
            p = doc.add_paragraph()
            p.paragraph_format.left_indent = Inches(0.35)
            p.paragraph_format.first_line_indent = Inches(-0.35)
            emit(p, txt, size=10)
            continue
        txt, i = gather(i)
        p = doc.add_paragraph()
        p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
        emit(p, txt)

    # 图：统一放在正文之后（CPC 允许随文或集中）
    if with_figs:
        doc.add_page_break()
        h = doc.add_paragraph(); set_font(h.add_run("Figures"), size=13, bold=True)
        for fn, (num, cap) in (figmap or FIGCAP).items():
            path = os.path.join(FIGDIR, fn)
            if not os.path.exists(path):
                print("   ⚠ 缺图 %s" % fn); continue
            # ★ 不能一律按宽度插：Fig.5 改成单列四格后是 140×211 mm 的竖图，
            #   按 6.4 in 宽放进来高度会到 9.7 in，**超出版心**（信纸 1 in 页边距下约 9 in）。
            #   所以按宽高比取"宽 6.4 in"与"高 8.2 in"里更小的那个。
            with Image.open(path) as _im:
                _w, _h = _im.size
            if _h / _w * 6.4 > 8.2:
                doc.add_picture(path, height=Inches(8.2))
            else:
                doc.add_picture(path, width=Inches(6.4))
            doc.paragraphs[-1].alignment = WD_ALIGN_PARAGRAPH.CENTER
            cp = doc.add_paragraph()
            cp.paragraph_format.line_spacing_rule = WD_LINE_SPACING.SINGLE
            cp.paragraph_format.space_after = Pt(16)
            set_font(cp.add_run(num + " "), size=10, bold=True)
            emit(cp, cap, size=10)
    doc.save(out_path)
    return out_path


if __name__ == "__main__":
    a = build(os.path.join(HERE, "manuscript.md"),
              os.path.join(HERE, "KinSTEP_CPC_manuscript.docx"),
              with_summary=True, with_figs=True)
    b = build(os.path.join(HERE, "SI.md"),
              os.path.join(HERE, "KinSTEP_CPC_SI.docx"),
              with_figs=True, figmap=SIFIGCAP)
    for f in (a, b):
        print("  %s  (%.1f KB)" % (os.path.basename(f), os.path.getsize(f) / 1024))
