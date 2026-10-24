# KinSTEP deposit notes (CPC)

> This file answers the three questions a referee or reader will ask: **what is in the deposit**,
> **how third-party material is handled**, and **how far the archive runs once downloaded**.
> Created 2026-10-23. A Chinese version is kept as `RELEASE_NOTES.zh-CN.md`; when the distribution
> policy changes, this file and `ATTRIBUTION.md` are changed together.

---

## 1. Licence

| object | licence |
|---|---|
| KinSTEP itself (`kinstep/`, `tests/`, `examples/`, `run_all_tests.sh`, documentation) | **GPL-3.0-or-later**; the full text is in `LICENSE`, the unmodified GNU original |

`pyproject.toml` declares `license = "GPL-3.0-or-later"` and `license-files = ["LICENSE"]`, and
every source file under `kinstep/` carries the GPL notice in its header.

---

## 2. Third-party material: what the deposit does **not** contain, and why

This was the only judgement call that had to be made before release.

### 2.1 CATKINAS example data — not redistributed

CATKINAS (Chen, Wang and Hu, East China University of Science and Technology) is a **closed-source
distribution** (MATLAB pcode, `.p`), obtained from `catkinas.com/downloads/`. Section 2.2 of its
manual grants **use subject to citation in publications**, and **grants no redistribution right**.

Its example outputs `Examples/result_curve/data1.mat` and `Examples/result_map/data1.mat` are
therefore **not shipped with this deposit**.

What *is* shipped is `tests/catkinas_refdata.npz` (0.2 MB): the reference **arrays that this
suite compares against**, extracted from those two files — fourteen of the twenty-three variables
of the first and seven of the twenty-three of the second — by `tests/make_catkinas_refdata.py`.
The distinction is the one drawn above and in the paper: publishing the values one has compared
against is normal scientific practice, and is a different act from redistributing the software or
its example files. The practical consequence is that **every assertion and every figure in this
deposit runs on a clean download**, which was not true before: six suites and Fig. 2 used to be
skipped. When the original distribution is present the tests read it instead, and one assertion
compares the two element by element, so the archived copy cannot drift away from its source
unnoticed.

> Verifiability is unaffected: every comparison value reported digit by digit in the paper is
> archived with the Supplementary Information (`paper/SI/fig2_data.npz` holds the 401-point CATKINAS
> reference rates and coverages; `paper/SI/fig3d_drc.dat` and `figS1b_dps.dat` likewise). Publishing
> the reference values one has compared against is normal and necessary in a scientific paper, and
> is a different matter from redistributing the software.

To re-run those suites, obtain CATKINAS and then

```bash
export KINSTEP_CATKINAS_DIR=/path/to/CATKINAS     # the level that contains Examples/
./run_all_tests.sh
```

### 2.2 The upstream CatMAP clone — not shipped, but not for licence reasons

A full shallow clone of CatMAP (**GPLv3**, SUNCAT), 9.3 MB, exists in the development tree. GPLv3
permits redistribution, so it could legally be included; it is left out for two other reasons.

1. **It is useless for the regression suite.** Not one byte of the reference values comes from
   CatMAP: its repository contains reference output for two coverage points only, with no TOF, rate
   or DRC data. CatMAP is a **functional comparison**, not a source of reference values.
2. **Including it would invite the misreading** that KinSTEP is derived from CatMAP. In fact no code
   was copied and `catmap` is never imported — and without an import there is no derived work in the
   sense of the GPL.

Section 3 of `ATTRIBUTION.md` states the functional overlap with CatMAP honestly, along with the
differences that justify the project.

### 2.3 The Koper CO-stripping reference implementation — not shipped

The reference for suite T5 is the model of Koper et al. 1998 (reference [29],
doi:10.1063/1.477230) and the public Jupyter notebook implementation by J. W. Ager. T5 **reads no
external file**: the model is re-implemented from the paper and each reference value is written in
the script with its source. T5 is therefore self-contained.

### 2.4 NIST thermochemical data

Shomate coefficients are taken from the **primary entries of the NIST Chemistry WebBook (SRD 69)**
and annotated one by one; nothing is transcribed from a secondary handbook. Section 2 of
`ATTRIBUTION.md` records this self-correction.

---

## 3. How far the archive runs

```bash
pip install -e .                 # needs only numpy + scipy
pip install -e ".[all]"          # plus mpmath (arbitrary precision) and matplotlib (plotting)
./run_all_tests.sh
```

The suite contains **82 assertions**, and **all of them run on a clean download** — the tier that
used to require third-party software is gone, because the reference *values* it compares against are
now archived with the deposit (section 2.1). What remains is a smaller distinction: five assertions
need an optional extra.

| assertions | what they need | behaviour without it |
|---|---|---|
| 77 | NumPy and SciPy only | — |
| 4 | the optional `mpmath` extra (values below the double-precision floor) | reported as *skipped*, with the reason, never as failed |
| 1 | the optional `matplotlib` extra (a figure export) | as above |

> Skipping is not silent: the summary line prints the number skipped, and each skipped assertion
> prints why. Setting `KINSTEP_CATKINAS_DIR` makes the comparisons read the original distribution
> instead of the archived arrays, and adds one assertion's worth of element-by-element checking
> that the two agree.

**The three configurations, all measured** (the last two by installing a freshly built source
distribution into an empty virtual environment):

| configuration | result |
|---|---|
| optional extras installed, reference distribution present | **82 passed / 0 failed / 0 skipped** |
| optional extras installed, reference distribution absent | **82 passed / 0 failed / 0 skipped** |
| only numpy and scipy | **77 passed / 0 failed / 5 skipped** |

The five skips in the third configuration are the designed behaviour of optional extras, not
defects: four assertions compare against values that cannot be resolved in double precision at all
(so without `mpmath` there is nothing to compare), and one exports a figure (so without
`matplotlib` there is no figure). Each prints its own reason.

**Figures.** Three scripts produce the five figures: `paper/make_figures.py` (Figs. 1, 2, 3, S1),
`paper/make_fig4_maps.py` (Fig. 4) and `paper/make_fig5_cli.py` (Fig. 5).

- Fig. 1, Fig. 3(a–c) and Fig. S1(a) are fully self-contained.
- Fig. 2 needs the CATKINAS reference data to be **recomputed** (without it the whole figure is
  skipped; the values are in `SI/fig2_data.npz`).
- Fig. 3(d) and Fig. S1(b) **fall back to the archived values in SI and redraw the same panel** when
  the data is absent: the figure is reproducible, it is simply not recomputed.
- Fig. 4 is fully self-contained (about 20 minutes for 4992 points); `--replot` redraws it from
  `SI/fig4_maps.npz` so that adjusting the layout costs nothing.
- Fig. 5 is fully self-contained and **stores no text**: every run re-executes the four commands
  with `KINSTEP_LANG=en` and typesets whatever they print, so the figure cannot drift away from the
  code. It carries its own gate — if any Chinese character appears in the output it refuses to draw.
  The complete session is archived as `SI/fig5_cli_session.txt`.

---

## 4. Contents of the deposit

```
LICENSE                  full GPLv3 text
README.md                installation and usage (English; README.zh-CN.md is the Chinese version)
EXAMPLES.md              the six examples, one by one, with numbers from runnable scripts
ATTRIBUTION.md           what is borrowed and whom to cite
RELEASE_NOTES.md         this file
CITATION.cff             how to cite the software
pyproject.toml           GPL-3.0-or-later; core dependencies are numpy and scipy only
MANIFEST.in              prune refdata_local — keeps the material of section 2 out of the distribution
run_all_tests.sh         82 assertions
kinstep/                 the package, 18 modules, about 6900 lines
tests/                   14 numerical suites, about 2200 lines; _root.py and refdata.py make paths portable
examples/                six runnable examples (ex1 … ex6)
paper/
  manuscript.md          manuscript text
  SI.md                  Supplementary Information
  make_figures.py        Figs. 1, 2, 3, S1
  make_fig4_maps.py      Fig. 4 (two-dimensional maps for four modalities; --replot redraws from npz)
  make_fig5_cli.py       Fig. 5 (command-line interface; re-runs the commands, stores no fixed text)
  build_docx.py          markdown to Word
  highlights.txt         highlights for submission
  figures/               all figures (*.pdf is shipped, *.png is not)
  SI/                    plottable raw data *.dat / *.npz and the full fig. 5 session
```

> **Every path is ASCII** (changed 2026-10-23). The directories and files previously had Chinese
> names, which is fine for local development, but a CPC deposit is downloaded and unpacked
> worldwide, and **a ZIP without the UTF-8 flag mangles such names**. File *contents* still contain
> Chinese comments, which is unaffected — only filename encoding was at issue.

`paper/` carries the **figure scripts, the plottable raw data of the SI (`*.dat`, `*.npz`) and the
vector figures (`*.pdf`)**. The `.docx` manuscript and the `.png` bitmaps are not included: the
manuscript goes to the journal, and the bitmaps can be regenerated by the scripts.

**Not in the deposit:** the third-party reference data of section 2, `notes/` (an internal
development log in Chinese), `notes/FAILED_*.log`, `__pycache__`, `paper/*.docx` and
`paper/figures/*.png`.

> `kinstep.egg-info/PKG-INFO` and `SOURCES.txt` do appear in the source distribution: setuptools
> regenerates and includes that metadata itself, `prune` in `MANIFEST.in` does not reach it, and it
> is harmless.

Packaging and self-check:

```bash
python -m build --sdist
A=dist/kinstep-*.tar.gz
tar tzf $A | grep -cE "refdata_local|catmap|data1\.mat|\.docx"   # forbidden items; must be 0
tar tzf $A | grep -c 'tests/catkinas_refdata.npz'       # the archived reference arrays; must be 1
tar tzf $A | grep -cP '[^\x00-\x7F]'                    # non-ASCII filenames; must be 0
tar tzf $A | grep -c 'paper/SI/fig5_cli_session.txt'    # the SI says it is deposited; must be 1
```

> All three checks are read as **counts**. `tar` also prints 0 when it cannot find the archive at
> all, which is a failure rather than a pass — confirm that `$A` exists first.

---

## 5. Portability (fixed 2026-10-23)

Three ways in which the code would only have run on the author's machine were found and removed
before release.

1. All 14 test scripts contained a hard-coded `sys.path.insert(0, "<developer's absolute path>")`.
   They now `import _root`, and `tests/_root.py` derives the package root from `__file__`.
2. Six test and figure scripts contained the absolute path of the CATKINAS installation. They now
   go through `refdata.catkinas(...)` / `catkinas_mat(...)`, which read `KINSTEP_CATKINAS_DIR` and
   either skip or fall back to the archived values when it is unset.
3. Those two helpers originally kept the developer's own path as a fallback. That fallback has been
   removed: the deposit is public, a path containing someone's user name is of no use to anyone
   else, and it would turn "the data is missing" into "the path does not match".

A related defect was fixed in the regression driver: in `if ! cmd; then … exit=$?`, the `!` resets
`$?` to 0, so the failure message always printed `exit=0`.
