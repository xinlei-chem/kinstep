# KinSTEP

**Kin**etics + **S**ono / **T**hermal / **E**lectro / **P**hoto — one core, four catalytic modalities.

Microkinetic modelling in pure Python: start from elementary reaction equations written the way you
would write them on paper, and obtain coverages, rates, degrees of rate control (DRC), reaction
orders, apparent activation energies, Tafel slopes and volcano plots.

```bash
git clone https://github.com/xinlei-chem/kinstep.git
cd kinstep

pip install -e .            # core needs only numpy + scipy
pip install -e ".[all]"     # plus mpmath (arbitrary precision) and matplotlib (plotting)

kinstep templates                # list the built-in mechanisms
kinstep new --template her       # write a skeleton with the mechanism already filled in
kinstep check her.kinstep.toml   # it tells you, slot by slot, what is still missing
kinstep run   her.kinstep.toml --drc
```

A Chinese version of this file is kept as [`README.zh-CN.md`](README.zh-CN.md). The command-line
interface is bilingual: set `KINSTEP_LANG=en` for English output (the default is Chinese).

---

## 1. What it is for

Existing microkinetic tools share two difficulties.

**1. They require an elaborate input script.** In CATKINAS the input is a piece of MATLAB code in
which a species name *is* a variable name, so a species called `T` silently overwrites the
temperature. KinSTEP takes the equations themselves:

```
N2 + 2* <-> 2N*
NH2* + H* <-> NH3 + 2*
H2O(l) + * <-> OH* + H+ + e-    {beta_ct: 0.5}
```

`*` is a site, `<->` reversible, `->` irreversible, `(l)` liquid phase, `e-` an electron. The
**electron stoichiometry `nu_elec` is counted automatically from the equation** — it is the single
parameter whose sign is easiest to get wrong by hand.

**2. They cover one modality.** A thermal-catalysis code does not know about electrochemistry, and an
electrochemical one does not know about photocatalysis. KinSTEP isolates *how an external stimulus
perturbs an elementary step* behind a **driver plug-in**, so the core never needs to know which kind
of catalysis it is solving:

| driver | what it changes | typical observables |
|---|---|---|
| `thermal` | nothing (baseline) | — |
| `electro` | `ΔG(U) = ΔG(0) − ν_e·e·U` (CHE / Butler–Volmer) | current density, Tafel slope |
| `photo` | `ΔG = ΔG(dark) − ν·ΔE_qfl` (quasi-Fermi-level splitting) | apparent quantum yield, carrier densities |
| `sono` | **not the energies but the operating point itself** (cavitation hot spots as weighted sub-states) | — |
| `field` | dipole + polarizability field effect | — |

Because `−ν·ΔE_qfl` and `−ν_e·e·U` are **mathematically isomorphic**, `photo + electro` compose
directly into photoelectrocatalysis (PEC) without a third body of code.

---

## 2. Three design principles that run through the whole package

### 2.1 A wrong energy does not raise an error — it returns a plausible-looking curve

Every number therefore carries a **provenance tag** (`provenance.py`):

```python
Val(1.21, "dft_file", detail="./OUTCAR")      # trusted
Val(0.80, "llm_suggested")                    # ★ for an energy ⇒ the solver refuses to run
Val(0.50, "default")                          # allowed, but a default must be visible
```

There are nine provenance levels. For quantities of energy type, `llm_suggested`, `guess` and
`REQUIRED` **cause the solver to refuse**, with a slot-by-slot explanation of how to fix it.

### 2.2 An acceptance criterion must be about whether the result can actually be used

Each gate in this package targets a class of error that **does not raise, and only bends the
conclusion**:

| gate | what it catches | what it actually caught |
|---|---|---|
| element / charge balance | typos | `H2O* -> OH* + H*` in a template: sites 1 ≠ 2 |
| site conservation | desorption that forgets to release the site | `NH2* + H* <-> NH3 + *` |
| reserved-name collision | a species named like a physical quantity | `T`, `theta`, `E` |
| **Wegscheider cycle** | one wrong energy in one step | a 0.1 eV change is caught immediately |
| **Campbell sum rule ΣX_B = 1** | a wrong DRC | also used to **choose the finite-difference step automatically** |
| **flux conservation νᵀr = 0** | inconsistent rates of quasi-equilibrated steps | 1.55e-07 → 2.22e-16 |
| **charge balance (photo)** | mismatched `nu_elec_cb` / `nu_hole_vb` | electron and hole consumption differing by 90× |
| **weight normalisation Σw = 1 (sono)** | a meaningless multi-state average | a weight sum of 1.4 was stopped |

### 2.3 The precision limit has to be measured, not assumed

In the factorised rate law the reversibility factor `1−z` suffers catastrophic cancellation on
quasi-equilibrated steps. We measured where the limit actually is:

- the bottleneck is **not** `expm1` (which is already full precision) but the **absolute accuracy of
  `ln z = ΔG/k_BT + Σν·ln a`**;
- in double precision, `1−z ≈ 1e-30` leaves not one correct significant digit;
- hence the solver is a **four-level ladder**, every level of which was forced by a measured failure:

```
(1) ODE integration (adaptive continuation, wall-clock budget)  -> physical initial guess
(2) damped Newton (analytic Jacobian, step limiting, backtracking) -> double-precision solution
(3) flux projection (weighted, sigma = R_f * eps)  -> rates of quasi-equilibrated steps from conservation
(4) arbitrary-precision Newton (mpmath, triggered automatically) -> 1-z actually resolved
```

`precision="auto"` decides for itself whether to escalate and records the reason in `notes`.

---

## 3. How far it is verified

| reference | quantity | result |
|---|---|---|
| **Koper 1998** (*J. Chem. Phys.* **109**, 6051), CO stripping | ODE right-hand side, point by point | **2.4e-15** |
| same | peak current (mA/cm²) | **all six digits** identical to the reference implementation |
| analytic closed form (formic-acid decomposition) | θ and r | **2.2e-16** |
| **CATKINAS, 401 points** (ammonia-synthesis descriptor scan) | coverages / reversibility | **1.4e-10** |
| same | net rate (with automatic escalation) | **2.7e-10**; points in error by more than 1e-6: **27 → 0** |
| same | DRC / Campbell sum | **3.8e-12 / 2.9e-12** |
| **CATKINAS, 315 points** (T–P map, flow reactor) | median outlet partial pressure | **2.2e-16**, 315/315 converged |
| NIST Shomate | five entropy-correction constants | **1.9e-10 eV** |
| photocatalysis literature | light-intensity power law (**computed, not assumed**) | 0.996 → **0.500** |

These comparisons are not "close enough by eye": they are digit-by-digit down to the accuracy floor
set by the input constants. Along the way, **five differences in constant conventions** relative to
CATKINAS were found (k_B, e, h, F, and the T inside k_flow); each was traced to the specific value
and handled explicitly, because leaving them in would mean comparing constants rather than formulae.

Full regression suite: `./run_all_tests.sh` (14 numerical suites, key-number regression checks, CLI
behaviour, and error paths).

---

## 4. Five commands

```bash
kinstep templates [name]      # 11 built-in systems (thermal/electro/photo/sono): mechanism only, no energies
kinstep new --template her    # write a skeleton
kinstep ask "oxygen evolution on IrO2 at 1.5 V vs SHE, 25 degrees C" --write
                              # one sentence -> structure + questions (no LLM involved)
kinstep check   xxx.toml      # four structural gates + provenance gate + driver; every error says how to fix it
kinstep explain xxx.toml      # plain language: what is being solved, which values are defaults
kinstep run     xxx.toml --drc  # solve -> report -> plot -> export .dat
kinstep drivers               # modalities, parameters, and which parameter is required by whom
```

Exit codes: `0` success, `2` bad input, `3` stopped by the provenance gate, `4` driver
misconfiguration.

### The red line in the natural-language layer

`kinstep ask` turns one sentence into **a structure plus a list of questions**. Numbers fall into two
classes and are treated completely differently:

- **numbers the user states** ("673 K", "pH = 0", "−100 mV") are extracted, tagged `user_given`, and
  **every unit conversion is recorded explicitly** ("400 degrees C" → `673.15 K (400 °C + 273.15)`);
- **numbers the system would have to invent** (any energy) **always become questions**.

"100 mW/cm²" needs a photon energy before it becomes a photon flux; without one it is **a question,
not an assumption**.

The LLM back end is **off by default**, requires an explicit `confirm=True`, and whatever it returns
must pass three gates: whitelist → mechanism parsing → provenance gate. Measured: for a hostile back
end that returned `Ga_fwd=0.42`, all that survived was `template_key`, `questions` and `modality`.

---

## 5. Documentation

- [`EXAMPLES.md`](EXAMPLES.md) — six runnable examples (thermal / electro / photo / sono / scan /
  flow reactor).
- [`ATTRIBUTION.md`](ATTRIBUTION.md) — which published methods are used, and the avoidance clauses.
- [`RELEASE_NOTES.md`](RELEASE_NOTES.md) — what is and is not redistributed, and why.
- `paper/` — the manuscript sources, the figure scripts and the raw data behind every figure.

## 6. Citing

Please cite the original literature for the methods, not this package:

- reversibility iteration and the factorised rate form: Chen, Wang, Hu, *ACS Catal.* **2016**, 6, 7078
- computational hydrogen electrode: Nørskov et al., *J. Phys. Chem. B* **2004**, 108, 17886
- degree of rate control: Campbell; Stegelmann, Andreasen, Campbell, *J. Am. Chem. Soc.* **2009**
- gas-phase thermochemistry: NIST Chemistry WebBook (SRD 69); Chase, *NIST-JANAF*, 4th ed., 1998
- electrochemical regression reference: Koper et al., *J. Chem. Phys.* **1998**, 109, 6051
- thermal-catalysis regression reference: CATKINAS, *J. Comput. Chem.* **2021**, 42, 379

If you need to cite the software itself, see [`CITATION.cff`](CITATION.cff).

## 7. Licence and distribution

**GPL-3.0-or-later**; the full text is in [`LICENSE`](LICENSE).

The policy on third-party material is in [`RELEASE_NOTES.md`](RELEASE_NOTES.md):

- The CATKINAS example output (`.mat`) is data shipped with a closed-source distribution and is **not
  licensed for redistribution**, so neither it nor the software is included here. What *is* included
  is `tests/catkinas_refdata.npz` (0.2 MB): the reference **arrays this suite compares against**,
  extracted from those files. Publishing the values one has compared against is ordinary practice
  and is a different act from redistributing the software — and it is what makes **every assertion
  and every figure reproducible from a clean download**. With the original distribution present
  (`KINSTEP_CATKINAS_DIR`), the tests read it instead and one assertion compares the two element by
  element.
- The upstream CatMAP clone is GPLv3 and could legally be included, but **not one byte** of this
  suite's reference values comes from it, so it is left out to avoid any suggestion that KinSTEP is
  derived from it.
- The suite contains 82 assertions, all of which run on a clean download. Measured: with the
  optional extras installed, **82 passed / 0 failed**, with or without the reference distribution;
  on a clean install with only numpy and scipy, **77 passed / 5 skipped**, the five skips being
  four that need `mpmath` to resolve values below the double-precision floor and one figure export
  that needs `matplotlib`.
- **The interface is bilingual.** `KINSTEP_LANG=en` switches all output to English (the default is
  Chinese), and one regression scans every command in English mode and asserts that **no Chinese
  character appears anywhere in the output, including in the files the program writes**.
