# Supplementary Information

**KinSTEP: a verified microkinetic solver spanning thermal, electro-, photo- and sonocatalysis**

---

## S1. Verification summary

All comparisons are digit-by-digit against the stated reference. Where the reference implementation
uses different physical constants, those constants were matched explicitly before comparing;
otherwise the comparison measures constants rather than formulations (see S2).

**Table S1.** Verification against four independent references.

| # | Reference | Quantity compared | Points | Max. relative deviation |
|---|---|---|---|---|
| 1 | Koper et al. 1998 [29] | ODE right-hand side | 8 random $(t,\theta)$ | $2.44\times10^{-15}$ |
| 2 | Koper et al. 1998 [29] | peak current, 3 initial coverages | 3 | six digits identical |
| 3 | Koper et al. 1998 [29] | potential-step transient | 200 | $1.41\times10^{-14}$ |
| 4 | analytic closed form | $\theta$, $r$ (formic-acid decomposition) | 6 regimes | $2.22\times10^{-16}$ |
| 5 | analytic closed form | long-time transient limit | 1 | $0.00$ (bitwise) |
| 6 | CATKINAS [11] | $k_f$, $k_r$ | $401\times5$ | $2.71\times10^{-10}$ |
| 7 | CATKINAS [11] | coverages $y$ | $401\times5$ | $1.41\times10^{-10}$ |
| 8 | CATKINAS [11] | reversibility $Z$ | $401\times5$ | $3.59\times10^{-11}$ |
| 9 | CATKINAS [11] | net rate (double only) | 401 | $5.60\times10^{-2}$ |
| 10 | CATKINAS [11] | net rate (`precision=auto`) | 401 | $2.71\times10^{-10}$ |
| 11 | CATKINAS [11] | step DRC $X_B$ | 7 descriptor points | $3.81\times10^{-12}$ |
| 12 | CATKINAS [11] | intermediate DRC $X_I$ | 7 descriptor points | $4.23\times10^{-12}$ |
| 13 | CATKINAS [11] | reaction orders | 7 descriptor points | $2.44\times10^{-10}$ |
| 14 | CATKINAS [11] | apparent activation energy | 1 | $3.21\times10^{-8}$ eV |
| 15 | CATKINAS [11] | CSTR outlet pressures (median) | 315 | $2.22\times10^{-16}$ |
| 16 | NIST Shomate [21] | entropy corrections, 5 steps | 5 | $1.86\times10^{-10}$ eV |
| 17 | internal | double vs. arbitrary precision, $\theta$ | 4 points | $1.55\times10^{-15}$ |
| 18 | internal | Campbell sum rule $\sum X_B - 1$ | 7 points | $2.94\times10^{-12}$ |
| 19 | CatMAP [5] | coverages, same solution branch | 193 | $8.35\times10^{-10}$ (median $3.08\times10^{-31}$) |
| 20 | CatMAP [5] | step rates, where CatMAP's own mass balance holds to $10^{-10}$ | 91 | median $2.2$–$2.4\times10^{-15}$ (3 steps) |
| 21 | internal | mass-balance residual $\max(\lvert r_1-r_3\rvert,\lvert 2r_2-r_3\rvert)$ on the CatMAP map | 193 | $0$ (exact) |

**Table S2.** Full-chain reproduction of the CATKINAS 401-point scan starting **only** from the
four published scaling/BEP coefficients per step (no energies taken from the reference output).

| Quantity | Max. | Median |
|---|---|---|
| net rate $r$ | $2.705\times10^{-10}$ | $1.129\times10^{-10}$ |
| coverages $y$ | $1.412\times10^{-10}$ | $4.593\times10^{-11}$ |
| reversibility $Z$ | $3.590\times10^{-11}$ | $3.587\times10^{-11}$ |

The maximum agrees with that obtained when the reference energies are supplied directly
($2.706\times10^{-10}$), showing that the scaling, BEP and thermochemistry layers introduce no
additional error.

---

## S2. Constant conventions

Five differences in physical constants between this work (CODATA-2018) and CATKINAS were located
during verification. Each was identified by the signature of a *constant* relative offset across
many points rather than scatter.

**Table S3.** Constant-convention differences and their effect.

| Constant | CATKINAS | This work | Relative diff. | Effect on comparison |
|---|---|---|---|---|
| $k_B$ | $8.617343\times10^{-5}$ eV K$^{-1}$ (CODATA-2006) | CODATA-2018 | $1.15\times10^{-6}$ | $\times E_a/k_BT\approx25$ ⇒ $4\times10^{-5}$ in $k_f$ |
| $e$ | $1.60217653\times10^{-19}$ C | $1.602176634\times10^{-19}$ C | $6.5\times10^{-8}$ | enters $k_BT$ in eV |
| prefactor $k_BT/h$ | $1.402306113158\times10^{13}$ s$^{-1}$ (673 K) | $1.402304467\times10^{13}$ | $1.17\times10^{-6}$ | directly proportional in $r$ |
| $F$ (J mol$^{-1}$ eV$^{-1}$) | $96485.2608$ | $96485.3321$ | $7.39\times10^{-7}$ | $\times T\Delta S/k_BT$ ⇒ $1.9\times10^{-5}$ |
| $T$ in $k_\mathrm{flow}$ | nominal $T$ | equivalent $T$ | $1.6\times10^{-6}$ | constant $1.13\times10^{-6}$ offset in outlet $p$ |

The `params` array stored in the CATKINAS output files contains $N_A = 6.022137\times10^{23}$ and
$k_B = 1.380651\times10^{-23}$, i.e. 1986-era values, which is the origin of the first four entries.

Two errors in this package's own constants were found and corrected in the same process: `Q_E` and
`EV` had been given different numerical values although they are the same constant by definition
($1\ \mathrm{eV}\equiv1\ e\times1\ \mathrm{V}$), and `R_GAS` was a CODATA-2014 value while $N_A$ and
$h$ were CODATA-2018, so that $k_B=R/N_A$ missed the exact 2018 value by $3.4\times10^{-7}$.

---

## S3. Precision of $1-z$

**Table S4.** Relative error in $1-z$ computed two ways (double precision).

| true $1-z$ | $1-\exp(\ln z)$ | rel. error | $-\mathrm{expm1}(\ln z)$ | rel. error |
|---|---|---|---|---|
| $10^{-4}$ | $1.000000\times10^{-4}$ | $1.10\times10^{-13}$ | $1.000000\times10^{-4}$ | $0$ |
| $10^{-8}$ | $1.000000\times10^{-8}$ | $5.02\times10^{-9}$ | $1.000000\times10^{-8}$ | $0$ |
| $10^{-12}$ | $9.999779\times10^{-13}$ | $2.21\times10^{-5}$ | $1.000000\times10^{-12}$ | $0$ |
| $10^{-14}$ | $9.992007\times10^{-15}$ | $7.99\times10^{-4}$ | $1.000000\times10^{-14}$ | $0$ |
| $10^{-15}$ | $9.992007\times10^{-16}$ | $7.99\times10^{-4}$ | $1.000000\times10^{-15}$ | $0$ |

`expm1` removes the subtraction error entirely; the remaining limit is the absolute accuracy of
$\ln z$ ($\sim10^{-16}$), giving a relative error in $1-z$ of $\varepsilon/(1-z)$.

Fig. S1(a) shows the three regimes over the full range. Note what it also shows: once $\ln z$
carries only absolute accuracy — which is the generic case here, since
$\ln z = \Delta G/k_BT + \sum_j \nu_j \ln a_j$ is a near-cancelling sum of terms of order unity —
`expm1` and the naive subtraction become equally accurate, both tracking
$\varepsilon/(1-z)$. `expm1` is therefore necessary but not sufficient: it guarantees that the
accuracy of $1-z$ is never worse than that of $\ln z$, and it is what makes the escalation
criterion of §S3 meaningful, but it cannot by itself recover a rate whose $1-z$ has fallen below
$\varepsilon$. That is what the arbitrary-precision level is for, and Fig. S1(b) shows it
working.

**Table S5.** Convergence of the arbitrary-precision result with decimal places, at
$E_\mathrm{N}=+2.00$ eV where the reference net rate is $4.801095\times10^{-32}$.

| dps | $r$ | rel. dev. from reference |
|---|---|---|
| 20 | $-1.394\times10^{-27}$ | $2.90\times10^{+4}$ (sign wrong) |
| 30 | $4.8010950942\times10^{-32}$ | $2.94\times10^{-10}$ |
| 40 | $4.8010950942\times10^{-32}$ | $3.07\times10^{-10}$ |
| 50 | $4.8010950940\times10^{-32}$ | $2.71\times10^{-10}$ |
| 60 | $4.8010950940\times10^{-32}$ | $2.71\times10^{-10}$ |
| 80 | $4.8010950940\times10^{-32}$ | $2.71\times10^{-10}$ |

The result plateaus from 30 digits; the default of 50 provides margin. The plateau value equals the
precision of the reverse-engineered input constants, not of the arithmetic. These data are plotted
in Fig. S1(b).

---

## S4. Flux projection

**Table S6.** Effect of the error-bar choice in the weighted projection onto
$\{r:\nu^{\mathsf{T}}r=0\}$, ammonia synthesis at 700 K.

| Quantity | No projection | $\sigma=|r|\varepsilon/(1-z)$ | $\sigma=R_f\varepsilon$ |
|---|---|---|---|
| N-element flux residual | $1.4\times10^{-15}$ | — | $4.4\times10^{-16}$ |
| H-element flux residual | $1.55\times10^{-7}$ | — | $2.22\times10^{-16}$ |
| rates at a point where $1-z$ underflows | — | all low by factor 14 (ratios exact) | correct |

With $\sigma=|r|\varepsilon/(1-z)$, a step whose $1-z$ has underflowed has $r=0$ and therefore a
vanishing error bar, receiving maximum weight although it carries no information. The absolute
uncertainty of a net rate is $R_f\,\delta(1-z)\approx R_f\varepsilon$, which is well defined in that
limit.

---

## S5. Degree-of-rate-control truncation

**Table S7.** Deviation of the intermediate DRC from the reference as a function of the
central-difference step $\delta$, demonstrating pure second-order truncation.

| $\delta$ (eV) | deviation | deviation$/\delta^2$ |
|---|---|---|
| $4.00\times10^{-3}$ | $1.524\times10^{-4}$ | $9.523$ |
| $2.00\times10^{-3}$ | $3.810\times10^{-5}$ | $9.525$ |
| $1.00\times10^{-3}$ | $9.526\times10^{-6}$ | $9.526$ |
| $5.00\times10^{-4}$ | $2.381\times10^{-6}$ | $9.526$ |
| $2.50\times10^{-4}$ | $5.954\times10^{-7}$ | $9.526$ |
| $1.25\times10^{-4}$ | $1.488\times10^{-7}$ | $9.526$ |

Constancy of the third column over a 32-fold range establishes that the residual is truncation, not
a definitional difference. Richardson extrapolation $X\approx[4X(\delta/2)-X(\delta)]/3$ then
reduces the deviation by a factor of $2.3\times10^{5}$ (Table S8).

**Table S8.** Effect of Richardson extrapolation at $E_\mathrm{N}=0$.

| Quantity | $\delta=10^{-3}$ | extrapolated | improvement |
|---|---|---|---|
| $X_B$ (step 2) | $9.94\times10^{-11}$ | $5.63\times10^{-13}$ | $\times177$ |
| $X_I$ (H*) | $9.53\times10^{-6}$ | $4.11\times10^{-11}$ | $\times2.3\times10^{5}$ |
| $X_I$ (*) | $9.52\times10^{-6}$ | $4.21\times10^{-11}$ | $\times2.3\times10^{5}$ |
| $X_{PC}$ (N$_2$) | $2.26\times10^{-7}$ | $1.21\times10^{-12}$ | $\times1.9\times10^{5}$ |

---

## S6. Gas-phase thermochemistry self-check

**Table S9.** Shomate coefficients recomputed against each species' own reference values at
298.15 K. Tolerance $0.05$ J mol$^{-1}$ K$^{-1}$.

| Species | $S^\circ$ computed | $S^\circ$ ref. | $\Delta S$ | $C_p^\circ$ computed | $C_p^\circ$ ref. | $\Delta C_p$ |
|---|---|---|---|---|---|---|
| H$_2$ | 130.68 | 130.68 | $+0.000$ | 28.84 | 28.84 | $-0.003$ |
| N$_2$ | 191.61 | 191.61 | $-0.001$ | 29.12 | 29.12 | $+0.004$ |
| O$_2$ | 205.15 | 205.15 | $-0.003$ | 29.38 | 29.38 | $+0.003$ |
| CO | 197.66 | 197.66 | $+0.003$ | 29.15 | 29.14 | $+0.009$ |
| CO$_2$ | 213.79 | 213.79 | $-0.002$ | 37.13 | 37.13 | $-0.000$ |
| NH$_3$ | 192.77 | 192.77 | $+0.000$ | 35.65 | 35.06 | $+0.000$ * |
| CH$_4$ | 186.25 | 186.25 | $+0.005$ | 35.65 | 35.69 | $-0.042$ |
| NO | 210.76 | 210.76 | $+0.003$ | 29.86 | 29.86 | $-0.003$ |
| N$_2$O | 219.96 | 219.96 | $-0.004$ | 38.61 | 38.62 | $-0.009$ |
| H$_2$O (g) | — | — | skipped | — | — | range starts at 500 K |

\* NH$_3$ carries a documented end-point deviation of $+0.587$ J mol$^{-1}$ K$^{-1}$ in $C_p$,
which is subtracted before comparison. This is an intrinsic property of the NIST fit at the lower
limit of its range (298–1400 K): the deviation falls to $\le0.3\%$ between 400 and 1000 K, while
$S^\circ(298.15)$ is reproduced exactly and $S^\circ(673)$ agrees with the value back-derived from
the CATKINAS entropy correction. Recording the known deviation in the data allowed the overall
tolerance to be tightened sevenfold (0.35 → 0.05) rather than loosened.

**Table S10.** Entropy corrections at 673 K obtained from one rule (Section 2.3).

| Step | Reaction | $\delta G^\ddagger$ (eV) | $\delta\Delta G$ (eV) | reference | deviation |
|---|---|---|---|---|---|
| 1 | H$_2$ + 2* $\leftrightarrow$ 2H* | $+1.077301$ | $+1.077301$ | $+1.077301$ | $1.9\times10^{-10}$ |
| 2 | N$_2$ + 2* $\leftrightarrow$ 2N* | $+1.504256$ | $+1.504256$ | $+1.504256$ | $1.9\times10^{-10}$ |
| 3 | N* + H* $\leftrightarrow$ NH* + * | $0$ | $0$ | $0$ | $0$ |
| 4 | NH* + H* $\leftrightarrow$ NH$_2$* + * | $0$ | $0$ | $0$ | $0$ |
| 5 | NH$_2$* + H* $\leftrightarrow$ NH$_3$ + 2* | $0$ | $-1.575974$ | $-1.575974$ | $1.9\times10^{-10}$ |

---

## S7. Photocatalysis: conditions for the half-order law

The high-intensity power-law exponent is $\beta\nu f$ (or $\beta\nu(1-f)$ for a hole-driven step),
where $\beta$ is the transfer coefficient, $\nu$ the number of carriers transferred in the
rate-determining step and $f$ the fraction of the quasi-Fermi-level splitting carried by the
electron level.

**Table S11.** Measured exponents (20-point flux scans over 13 decades).

| Case | low-$I$ exponent | high-$I$ exponent | predicted | note |
|---|---|---|---|---|
| $\nu=1$ (hole), $f=0.5$, bimolecular recombination | $0.9963$ | $0.5000$ | $0.5$ | classical case |
| $\nu=1$ (hole), $f=0.8$ | — | $0.2000$ | $0.2$ | $p$-type |
| $\nu=1$ (hole), $f=0.2$ | $1.000$ | $1.0000$ | $0.8$ | crossover flux $G^*=3.8\times10^{18}$ lies outside the scan |
| $\nu=1$, monomolecular (SRH) recombination | $1.0000$ | $0.9982$ | $1$ | no half-order regime |
| $\nu=2$ (1 e$^-$ + 1 h$^+$) | $1.000$ | $1.000$ | $1$ | full $\Delta E_\mathrm{qfl}$ acts |

The crossover generation rate follows from equating $r=G$ and $k_\mathrm{bi}e^{x}=G$:

$$ G^{*} = \Big[r_\mathrm{dark}^{1/(1-f)}\big/k_\mathrm{bi}\Big]^{(1-f)/f}, $$

which for $f=0.2$ gives $3.8\times10^{18}$, outside the scanned range — the model therefore
correctly remains first order there.

Carrier balance is checked at every point; the maximum imbalance
$|R_{e^-}-R_{h^+}|/\max(R_{e^-},R_{h^+})$ over the 14-point scan is $7.5\times10^{-16}$.

---

## S8. Sonocatalysis: multi-state weighting

**Table S12.** Weighted rate versus rate at the mean temperature, hot spot at 5000 K,
bulk at 350 K, barrier 1.95 eV.

| $f$ | hot-spot contribution | $\langle r(T)\rangle$ | $r(\langle T\rangle)$ | ratio |
|---|---|---|---|---|
| $1.793\times10^{-30}$ | $0.10\%$ | $2.2227\times10^{-16}$ | — | — |
| $1.793\times10^{-27}$ | $50.00\%$ | $4.4409\times10^{-16}$ | — | — |
| $1.793\times10^{-24}$ | $99.90\%$ | $2.2227\times10^{-13}$ | — | — |
| $1.0\times10^{-3}$ | $100.00\%$ | $1.2382\times10^{+8}$ | $8.8818\times10^{-16}$ | $1.39\times10^{23}$ |

The crossover fraction $f^{*}=r_\mathrm{bulk}/r_\mathrm{hot}=1.793\times10^{-27}$ predicted
analytically coincides with the numerically observed 50 % contribution point.

**Table S13.** Apparent activation energy versus hot-spot fraction.

| $f$ | $E_\mathrm{app}$ (eV) | regime |
|---|---|---|
| $0$ | $1.4634$ | purely thermal |
| $1.79\times10^{-29}$ | $1.4555$ | bulk dominates |
| $1.79\times10^{-27}$ | $0.9673$ | equal contributions |
| $1.79\times10^{-25}$ | $0.0309$ | hot spot dominates |
| $10^{-6}$ | $0.0000$ | hot spot only |

Literature hot-spot temperatures vary substantially with the measurement method: 5000–20000 K from
single-bubble sonoluminescence spectra [40], 750–6000 K from multi-bubble sonochemistry [38,39].
KinSTEP
treats these as inputs and warns outside both ranges.

---

## S9. Self-check output and regression

The distributed test suite comprises 15 numerical regression scripts plus a repository-wide
syntax sweep, command-line checks, interface-language checks and error-path checks (82 assertions in total), and verifies not
only exit status but the key numerical values quoted in this paper, so that a silent regression in
accuracy — or a silent drift between the code and the manuscript — is detected.

All 82 assertions run on a clean download. The comparisons of Section 4 read the reference arrays
archived with the deposit (`tests/catkinas_refdata.npz`, 0.2 MB, the fourteen variables of the
descriptor-scan output and the seven of the flow-reactor output that the assertions actually use);
the software that produced them is closed-source and is not redistributed. When a user has that
distribution and sets `KINSTEP_CATKINAS_DIR`, the tests read it directly instead and one assertion
compares the archive against it element by element, so the archived copy cannot drift from its
source unnoticed.

Five assertions depend on an optional extra and are reported as *skipped*, not failed, when it is
absent, each printing its reason: four compare against values below the double-precision floor and
therefore need `mpmath`, and one exports a figure and needs `matplotlib`.

Three configurations were measured before submission, the latter two from a freshly built source
archive installed into an empty virtual environment:

| Configuration | Result |
|---|---|
| Optional extras installed, reference distribution present | 82 passed / 0 failed / 0 skipped |
| Optional extras installed, reference distribution absent | 82 passed / 0 failed / 0 skipped |
| NumPy and SciPy only | 77 passed / 0 failed / 5 skipped |

The following were caught by the
self-checks during development and are listed as evidence that they are not decorative:

| Self-check | What it caught |
|---|---|
| element/charge balance | an unbalanced step in one of the distributed templates |
| site conservation | a dissociative adsorption step missing a vacancy |
| flux conservation | H-element flux error of $1.55\times10^{-7}$ from unprojected rates |
| carrier balance | electron/hole consumption differing by $90\times$ because a buffered intermediate decoupled the half-reactions |
| Campbell sum rule | also used to select the DRC finite-difference step |
| full residual vector | high-$T$/low-$p$ CSTR points reported as unconverged although the surface residual was $10^{-15}$; the gas-phase block was the culprit and had not been stored |

---

## S10. Data files

| File | Content |
|---|---|
| `fig2_data.npz` | 401-point scan: descriptor, double and auto rates, reference, coverages, high-precision flags |
| `fig3a_her_tafel.dat` | HER: potential, current density, Tafel slope |
| `fig3b_photo_powerlaw.dat` | photon flux, rate, AQY, local power-law order |
| `fig3c_sono.dat` | hot-spot fraction, weighted rate, rate at mean temperature |
| `fig3d_drc.dat` | descriptor, five step DRCs, sum |
| `figS1a_cancellation.dat` | true $1-z$; relative error of the naive form, of `expm1` with exact $\ln z$, and of `expm1` with $\ln z$ perturbed by $10^{-16}$ |
| `figS1b_dps.dat` | working precision, net rate of step 5, relative deviation from the reference |
| `fig4_maps.npz` | all four operating maps of Fig. 4: axes, net rate, coverages, current density, Tafel slope, quasi-Fermi splitting, local power-law order, weighted-to-mean rate ratio and selectivity (4992 solutions) |
| `fig4a_thermal.dat` | ammonia synthesis: $T$, $P$, net rate, $\theta_{\mathrm{H}^*}$, $\theta_*$ |
| `fig4b_electro.dat` | HER: potential, pH, current density, $\theta_{\mathrm{H}^*}$, Tafel slope |
| `fig4c_photo.dat` | photon flux, recombination coefficient, quasi-Fermi splitting, rate, local power-law order |
| `fig4d_sono.dat` | hot-spot fraction, hot-spot temperature, rate, weighted-to-mean ratio, selectivity |
| `catmap_co_oxidation.json` | the CatMAP CO-oxidation comparison of Section 4: for each of the 225 descriptor points, the rate constants, coverages and step rates exported from CatMAP (GPLv3; cite Medford et al., *Catal. Lett.* **145** (2015) 794), kept to 34 significant digits. CatMAP itself is not redistributed; only the solver output compared against is. |
| `fig5_cli_session.txt` | the complete terminal session of Fig. 5: the four commands, their full output and their exit codes, captured with `KINSTEP_LANG=en` in a clean temporary directory |

Figs. 1–3 and S1 are regenerated with `paper/make_figures.py`, Fig. 4 with
`paper/make_fig4_maps.py` and Fig. 5 with `paper/make_fig5_cli.py`. The last of these does not
store any text: it re-runs the four commands and typesets their real output, so the figure cannot
drift away from the code. All numerical claims can be reproduced with `./run_all_tests.sh`.

Two further files are part of the deposit rather than of a figure. `tests/catkinas_refdata.npz`
holds the reference arrays of Section 4 — the fourteen variables of the descriptor-scan output and
the seven of the flow-reactor output that the assertions use, 0.2 MB — so that every comparison can
be repeated without the reference software; `tests/make_catkinas_refdata.py` regenerates it from
that software for anyone who has it, and records the md5 of each source file.
`paper/threshold_sensitivity.py` produces the table of Section S11 in about twenty seconds, and `paper/compare_catmap.py` reproduces the CatMAP comparison of Section 4.

---

## S11. Sensitivity of the escalation threshold

Section 3.3 states the criteria by which the solver decides, point by point, whether to escalate to
arbitrary precision: any reversible step with $1-z\le0$ or $1-z<10^{-10}$, a flux projection that
moves the rates by more than $10^{-6}$ in relative terms, or a full residual above the convergence
tolerance. The first of these thresholds is the one that fires in practice, so it is the one whose
choice has to be defended.

The table below repeats the 401-point descriptor comparison of Section 4 — the same route, the same
continuation order, the same `dps = 60` — with that threshold moved over two orders of magnitude.
The script is `paper/threshold_sensitivity.py` and runs in about twenty seconds.

| escalation threshold on $1-z$ | points escalated (of 401) | max. relative deviation from reference | wall clock |
|---|---|---|---|
| none (double precision throughout) | 0 | $5.599\times10^{-2}$ | 1.7 s |
| $10^{-8}$ | 264 | $2.706\times10^{-10}$ | 4.6 s |
| $10^{-10}$ (default) | 236 | $2.706\times10^{-10}$ | 4.2 s |
| $10^{-12}$ | 235 | $2.706\times10^{-10}$ | 4.1 s |

Two things follow. First, the agreement with the reference is **identical to four significant
figures** across the whole range: $2.706\times10^{-10}$, which is the storage precision of the
reference values themselves. The result therefore does not depend on where the threshold is put,
only on whether escalation happens at all — without it the deviation is $5.6\times10^{-2}$, eight
orders of magnitude worse. Second, what the threshold does change is how much work is done: a
threshold two orders of magnitude looser escalates 264 points instead of 235, for no gain. The
default sits where the deviation has already reached the reference's own precision and further
escalation buys nothing.

Note that running the same study through the full scaling-and-thermochemistry chain of
Section 5.1, rather than reading the reference energies point by point, gives 234 escalated points
and a deviation of $2.992\times10^{-9}$. That is the stricter test — it reconstructs the energies
instead of reading them — and it leads to the same conclusion; the two routes are quoted separately
because mixing them would make the numbers in this paper impossible to reconcile.
