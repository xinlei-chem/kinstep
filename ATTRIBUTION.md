# Attribution

> This file records what KinSTEP borrows and whom it must cite. **Any new borrowing is registered
> here before the code is written.** Created 2026-10-22. A Chinese version is kept as
> `ATTRIBUTION.zh-CN.md`.

## 0. Ground rules

1. **No copied code.** As of 2026-10-22 everything under `kinstep/` was written for this project;
   no code was taken from CatMAP, PyCatKin, OpenMKM, CATKINAS or Zacros. Before any third-party
   code is reused, its licence is checked and the reuse is registered here.
2. **Methods must be cited.** A published method taken from someone else is cited in three places
   at once: the README, the methods section of the paper, and the header comment of the source file
   that implements it.
3. **Data comes from the primary source.** Data tables are never transcribed from a secondary
   handbook.

## 1. Methods that must be cited

### 1.1 Reversibility iteration and the factorised rate form (solver core)

```
r_i = kf_i · ∏c · (1 − z_i),    z_i ≡ r_i(reverse) / r_i(forward)
```

**Source:** Chen, J.; Mao, Y.; Wang, H.; Hu, P. "Reversibility Iteration Method for Understanding
Reaction Networks and for Solving Microkinetics in Heterogeneous Catalysis." *ACS Catal.* **2016**,
*6*(10), 7078–7087. **doi:10.1021/acscatal.6b02405** (the solver core of CATKINAS).

> Correction, 2026-10-23: this entry previously read "Chen, J.; Wang, H.; Hu, P. et al.", which
> **omitted the second author, Mao, Y.**, and paraphrased the title from memory. It has been checked
> field by field against Crossref and corrected, as has reference [16] of `paper/manuscript.md`.
> Lesson: verify every citation against the Crossref API rather than writing it from memory.

**How we learned of it:** reading `Examples/result_single/log1` of CATKINAS, where the
automatically derived rate expression appears as `r(1) = −P_H2·Q_v²·kf(1)·(z1−1)` rather than
`kf·A − kr·B`, which is what made the role of the factorised form in avoiding catastrophic
cancellation clear.

**Handling:** cited in the header of `kinstep/solve.py`, in the README, and in the methods section
of any paper. **It must never be presented as an invention of this project.**

### 1.2 The solver ladder (ODE for an initial guess, then Newton)

Standard numerical practice, but the **specific order that works** was read off the CATKINAS log:

```
direct high-precision Newton  (fails)
pure reversibility Newton     (fails)
ODE15S for an initial guess + reversibility Newton   (converges in 9 iterations)
```

**Handling:** the source file states "solver ladder inspired by CATKINAS (Chen & Wang, ECUST)".

### 1.3 Degree of rate control

Campbell's degree of rate control; a published method, cited normally. The **self-check** that the
degrees of rate control sum to 1 (the Campbell sum rule) also comes from that work.

### 1.4 Computational hydrogen electrode

Nørskov, J. K. et al., *J. Phys. Chem. B* **2004**, *108*, 17886. The standard method in
electrocatalysis.

### 1.5 Verification benchmark for the electrochemical core

Koper, M. T. M.; Jansen, A. P. J.; van Santen, R. A.; Lukkien, J. J.; Hilbers, P. A. J.
*J. Chem. Phys.* **1998**, *109*, 6051–6062, doi:10.1063/1.477230 (mean-field model plus lattice
Monte Carlo for CO oxidation on Pt), together with its public Python implementation by Joel W. Ager
(LBL), `EC CO stripping MKM.ipynb`. **This project uses it as the reference values of a regression
test** and does not copy its code.

## 2. Data sources

| data | not to be taken from | taken from |
|---|---|---|
| Shomate parameters (A–H), heat capacity and entropy | **not transcribed from CATKINAS manual appendix 7.2** | **NIST Chemistry WebBook**, primary entries, annotated one by one |
| vibrational frequencies (for ZPE) | as above | NIST / the original paper |
| relative atomic masses | — | IUPAC |

> Self-correction, 2026-10-22: an earlier version of this file stated that appendix 7.2 of the
> manual "can be turned into JSON directly". **That statement is withdrawn.** The primary source of
> those data is NIST; transcribing them from a secondary handbook carries both a copyright risk and
> the risk of importing transcription errors.

## 3. Functional overlap with comparable tools (the honest version)

### 3.1 What overlaps — **deliberately**, so that their examples can serve as regression references

| capability | already available in |
|---|---|
| descriptor scan, scaling relations, 1D/2D volcano plots | **CatMAP** (Medford et al. 2015, SUNCAT) |
| DRC, reaction orders, apparent activation energy | CatMAP / MKMCXX / CATKINAS |
| mean-field steady state and CSTR | PyCatKin / OpenMKM / CATKINAS |
| electrochemistry (CHE and potential dependence) | CatMAP's `hbond_electrochemical` mode and ORR tutorial |

**Conclusion:** if KinSTEP did no more than "core + scan + DRC", it would be essentially a Python
rewrite of CatMAP and would not be worth publishing. The overlap is a **means** — it is what makes
reference values available — and not the point.

### 3.2 What we could not find an existing implementation of (the reason the project exists)

| feature | state of the art |
|---|---|
| **a multi-modality driver layer** in which thermal, electro, photo, sono and field effects **compose** | CatMAP treats electrochemistry as a correction mode, not a pluggable layer; the other tools have nothing comparable |
| **photocatalytic and sonocatalytic MKM** | **searched for specifically; no tool found** |
| **provenance tracking, with unverified energies refusing to solve** | every tool assumes the numbers the user typed are right |
| **no input script** (interactive wizard plus `need()`) | CatMAP needs an `.mkm` file plus a Python script; CATKINAS needs a table of magic numbers |
| name-collision checking at registration; multi-state weighting | not available |

## 4. Conventions adopted to avoid confusion

| item | convention | reason |
|---|---|---|
| input file extension | **`.kinstep.toml`**, never `.mkm` | `.mkm` is CatMAP's input suffix; reusing it would confuse tool chains |
| package name | `kinstep` | `mkm` and `pymkm` are already taken on PyPI |
| command name | `kinstep` | no conflict |
| wording in the documentation | never "first" or "novel"; instead "we are not aware of" | a literature search is not a proof of absence |

## 5. Distribution policy (fixed 2026-10-23)

See `RELEASE_NOTES.md` for the full statement. In summary:

| third-party material | licence | shipped in the deposit? | reason |
|---|---|---|---|
| CATKINAS example output `.mat` (the files themselves) | closed source; manual §2.2 grants use subject to citation only | **no** | **no redistribution right granted** |
| the reference arrays this suite compares against, extracted from those files (0.2 MB) | as above | **yes** | publishing the values one compares against is ordinary scientific practice and is a different act from redistributing the software; it is what makes every assertion and figure reproducible from a clean download |
| full upstream clone of CatMAP | GPLv3 (redistribution permitted) | **no** | 9.3 MB of source that KinSTEP never imports; including it would invite the misreading that KinSTEP is derived from it |
| the solver output of CatMAP's own CO-oxidation tutorial (225 points: rate constants, coverages, step rates), used for the comparison of Section 4 | GPLv3 | **yes**, as `paper/SI/catmap_co_oxidation.json` | GPLv3 permits redistribution; the file carries its provenance and the required citation (Medford et al., *Catal. Lett.* **145** (2015) 794), and without it the comparison could not be reproduced from the deposit |
| NIST Shomate coefficients | public data, annotated entry by entry | yes (inside the code) | taken from the primary WebBook entries, not transcribed from a secondary source |
| Koper 1998 model | a published method | not applicable | re-implemented from the paper in T5; no external file is read |

KinSTEP itself is released under **GPL-3.0-or-later** (`LICENSE` is the unmodified GNU text).
`prune refdata_local` in `MANIFEST.in` keeps the non-redistributable material out of the source
distribution.

The deposit ships the reference **arrays** the suite compares against
(`tests/catkinas_refdata.npz`, 0.2 MB) but not the software or its example files, so **all 82
assertions run on a clean download**. Measured by installing a freshly built source distribution
into an empty virtual environment: **77 passed / 0 failed / 5 skipped** with only numpy and scipy,
the five skips being assertions that need an optional extra (`mpmath` or `matplotlib`); with both
extras installed, **82 passed / 0 failed**. All three configurations are tabulated in section 3 of
`RELEASE_NOTES.md`.

## 6. Open items

- [x] check and register the licence of any third-party code before reusing it — done, see section 5
- [x] take thermochemical data from NIST and annotate each entry — done in T11
- [ ] before publication: ask the CATKINAS authors whether the wording of entries 1.1 and 1.2 is
      acceptable to them
