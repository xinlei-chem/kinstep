# KinSTEP: a verified microkinetic solver spanning thermal, electro-, photo- and sonocatalysis

**Authors** — Xin Lei$^{a,*}$

$^a$ School of Chemical Engineering and Pharmacy, Wuhan Institute of Technology, Wuhan 430205, China

$^*$ Corresponding author. E-mail: xin.lei@wit.edu.cn

---

## Abstract

Microkinetic modelling turns first-principles energetics into measurable rates, coverages and
selectivities, but existing packages are built around a single catalytic modality and require users
to write elaborate input scripts. We present **KinSTEP**, a pure-Python package whose numerical core
is *modality-agnostic*: how an external stimulus perturbs an elementary step is isolated behind a
small driver interface, so thermal, electrochemical, photochemical and sonochemical catalysis share
one verified solver. Because the photochemical driving force enters as $-\nu\,\Delta
E_\mathrm{qfl}$ in exact analogy with the electrochemical $-\nu_e eU$, photo- and electro-drivers
compose into photoelectrocatalysis without additional code; a second hook splits a nominal operating
point into weighted sub-states, which sonocatalysis requires because the rate is exponentially
sensitive to temperature — for a representative cavitation hot spot $\langle r(T)\rangle$ and
$r(\langle T\rangle)$ differ by $1.4\times10^{23}$. The solver is a four-level ladder — stiff ODE
integration for a physical initial guess, a damped Newton iteration with an analytic Jacobian in
$\ln\theta$ space, a weighted projection of the net rates onto the flux-conservation manifold, and
an automatically triggered arbitrary-precision stage — each level of which was added in response to
a measured failure rather than by design. KinSTEP is verified digit-by-digit against four
independent references: a published CO-stripping model ($2.4\times10^{-15}$), an analytic
closed-form steady state ($2.2\times10^{-16}$), the 401-point descriptor scan and 315-point
$T$–$p$ maps of CATKINAS (coverages to $1.4\times10^{-10}$, net rates to $2.7\times10^{-10}$ with
automatic escalation, which reduces the number of points whose error exceeds $10^{-6}$ from 27 to
zero), and a CatMAP CO-oxidation map (coverages to $3.1\times10^{-31}$). Package-wide self-checks are free consequences of the formulation and are always evaluated.

**Keywords:** microkinetic modelling; heterogeneous catalysis; electrocatalysis; photocatalysis;
sonocatalysis; degree of rate control; volcano plot; arbitrary precision

---

## PROGRAM SUMMARY

| | |
|---|---|
| *Program Title* | KinSTEP |
| *CPC Library link to program files* | (to be added by Technical Editor) |
| *Developer's repository link* | https://github.com/xinlei-chem/kinstep (version v0.1.0 archived at https://doi.org/10.5281/zenodo.23613004) |
| *Code Ocean capsule* | (to be added by Technical Editor) |
| *Licensing provisions* | GPL-3.0-or-later |
| *Programming language* | Python 3.9+ (tested on CPython 3.12) |
| *Supplementary material* | Supplementary Information (verification tables, constant conventions, benchmark data, full self-check output) |
| *External routines/libraries* | NumPy and SciPy (required); `mpmath` and `matplotlib` (optional extras, used for arbitrary-precision escalation and for plotting respectively). No compiled extension, no build step. |
| *Typical running time* | A single steady state solves in well under a second. On one core of a desktop CPU the 401-point descriptor scan of Section 4 takes 4.2 s (1.7 s if arbitrary precision is disabled), the 315-point flow-reactor map 7.8 s, and the 4992-point four-modality map of Fig. 4 about 20 min. Peak resident memory stays below 100 MB in all of these, dominated by the NumPy/SciPy import footprint rather than by the model. |
| *Nature of problem* | Steady-state and transient microkinetic modelling of heterogeneous catalysis under different external stimuli (temperature, electrode potential, illumination, ultrasound). The governing equations are stiff, the unknowns span tens of orders of magnitude, and quasi-equilibrated elementary steps cause catastrophic cancellation in the net rate. Existing packages address one modality at a time and require modality-specific input scripts. |
| *Solution method* | Elementary steps are parsed from plain text and compiled to a stoichiometric matrix; electron stoichiometry is counted automatically. Coverages are solved in $\ln\theta$ space with an analytic Jacobian, using a four-level ladder (stiff ODE continuation $\to$ damped Newton $\to$ weighted flux projection $\to$ arbitrary-precision Newton, triggered by measured criteria). External stimuli enter through a driver plug-in that returns perturbed energies, or through a second hook returning weighted sub-states. Degrees of rate control, reaction orders, apparent activation energies and Tafel slopes are obtained by central differences with Richardson extrapolation, the step size being selected automatically from the Campbell sum rule. |
| *Additional comments* | Core dependencies are NumPy and SciPy only; `mpmath` (arbitrary precision) and `matplotlib` (plotting) are optional extras. A command-line interface provides input validation, plain-language explanation of a model, and a rule-based natural-language front end that produces *structure and questions* but never numerical energies. The deposit carries a regression suite of 82 assertions, all of which run on a clean download: 77 pass with NumPy and SciPy alone and five are skipped, with their reasons, because they need an optional extra (four need `mpmath`, one needs `matplotlib`); with both extras installed all 82 pass. The reference values compared against in Section 4 are archived with the deposit, so no third-party software is needed to reproduce any assertion or figure. The command-line interface is bilingual: all user-facing output is available in English via `KINSTEP_LANG=en`, and one regression asserts that no CJK character appears in English mode. |

---

## 1. Introduction

Microkinetic modelling connects computed reaction energetics to quantities that experiments actually
measure — turnover frequencies, coverages, reaction orders, apparent activation energies and
selectivities [1–4]. Mature and widely used packages exist: CatMAP introduced descriptor-based
screening built on adsorption-energy scaling relations [5–7], later gaining a robust
numbers-based solver [8]; MKMCXX has been applied to Fischer–Tropsch and syngas-conversion
networks [9]; OpenMKM
builds on Cantera for reactor-scale simulation [10]; CATKINAS contributed a multi-level
root-finding strategy together with automated sensitivity analysis [11]. Where lateral interactions
or adlayer ordering make the mean-field approximation untenable, graph-theoretical kinetic Monte
Carlo is the established alternative [12], and the trade-off between the two levels of description
has been reviewed in detail [13].

Two limitations recur across these tools.

**(i) One modality at a time.** Thermal, electrochemical, photochemical and sonochemical catalysis
are described by the same surface mass-action kinetics, differing only in how an external stimulus
perturbs the elementary steps. In practice, however, codes are organised around one modality.
Electrochemical extensions of thermal codes exist, but photo- and sonocatalysis remain largely
outside the scope of general MKM packages; a recent perspective notes that stimulus-dependent
microkinetic models are only beginning to appear [14].

**(ii) Modality-specific input scripts.** Mechanisms are frequently written in a domain-specific
format or, in the case of CATKINAS, directly as code, so that a species named `T` can overwrite the
temperature variable. The most error-prone parameter — the signed electron stoichiometry of each
electrochemical step — is typically supplied by hand.

A third, less often discussed, issue is numerical. The net rate of an elementary step is
$r = R_f(1-z)$, with $z$ the reversibility. Quasi-equilibrated steps have $1-z \to 0$, and in double
precision the net rate of such a step eventually becomes noise. CATKINAS addresses this by running
at 1000-digit precision by default [11]. That is effective, but it is paid for everywhere rather
than where it is needed — the run log distributed with its own ammonia example records 657 s on
eight parallel workers for the same 401 descriptor points that take 4.2 s on one thread here
(Section 4) — and running at high precision throughout also obscures *where* the double-precision
result is actually trustworthy.

KinSTEP addresses all three. Section 2 sets out the formulation and the driver abstraction;
Section 3 describes the solver ladder and the measurements that forced each level; Section 4 reports
verification against four independent references; Section 5 presents one example per modality.

## 2. Formulation

### 2.1 Mechanism and rate law

Elementary steps are written as they appear on paper,

```
N2 + 2* <-> 2N*
H2O(l) + * <-> OH* + H+ + e-        {beta_ct: 0.5}
CO*1 + O*2 -> CO2 + *1 + *2
```

where `*` denotes a site, `<->`/`->` reversible/irreversible steps, `(l)` the liquid phase and `e-`
an electron. The parser builds the stoichiometric matrix $\nu$ and counts the electron
stoichiometry $\nu_e$ **from the equation itself** (products minus reactants, positive for oxidation),
removing the single most sign-error-prone input.

Rates are evaluated in the factorised form

$$r_s = \mathcal{P}_s R_{f,s}\,(1-z_s), \qquad
R_{f,s} = k_{f,s}\prod_j a_j^{\,a^{r}_{sj}}, \qquad
\ln z_s = \frac{\Delta G_s}{k_BT} + \sum_j (a^{p}_{sj}-a^{r}_{sj})\ln a_j ,$$

with $a_j$ the activity (coverage, dimensionless pressure or concentration), $a^{r}$ and $a^{p}$ the
reactant- and product-side exponent matrices, and $\mathcal{P}_s$ an optional ensemble prefactor
giving the number of equivalent site pairs available to step $s$. $\mathcal{P}_s$ is a constant per
step and is not a coverage-dependent lateral-interaction term; adsorbate–adsorbate interactions are
not treated here (Section 6). Forward rate constants follow from transition-state theory,
$k_f=(k_BT/h)\exp(-G^\ddagger/k_BT)$ [15], or, for non-activated adsorption, from collision theory
with an explicit sticking coefficient; the factorised form itself follows the reversibility
formulation of Chen et al. [16]. Writing $\ln z$ directly from $\Delta G$ and the activities, rather than as
$\ln k_r - \ln k_f + \dots$, both halves the number of terms and avoids the overflow of $k_r$ for
strongly exergonic adsorption steps. The reverse rate constant is never supplied independently:
$k_r = k_f/K_\mathrm{eq}$ with $K_\mathrm{eq}=\exp(-\Delta G/k_BT)$, which is the only way to
guarantee thermodynamic consistency. Models that specify $k_f$ and $k_r$ separately — as several
published phenomenological models do — are accepted, but the implied $\Delta G$ is then back-derived
and reported for the user to check.

### 2.2 Drivers: one hook for the modality

The core never learns which modality is being simulated. A driver supplies perturbed energies,

$$ (G^\ddagger, \Delta G) \;\longrightarrow\; (G^\ddagger + \delta G^\ddagger,\; \Delta G + \delta\Delta G), $$

and the rate constants follow from the *same* correction regardless of how $k_f$ was obtained:
$k_f = k_{f,0}\exp(-\delta G^\ddagger/k_BT)$, $K_\mathrm{eq}=K_{\mathrm{eq},0}\exp(-\delta\Delta G/k_BT)$.

* **electro** — computational hydrogen electrode and Butler–Volmer [17,18]:
  $\delta\Delta G = -\nu_e eU$, $\delta G^\ddagger = -\beta\nu_e eU$.
* **photo** — the quasi-Fermi-level splitting $\Delta E_\mathrm{qfl}=E_{F\!n}-E_{F\!p}$ is the
  thermodynamic driving force for interfacial transfer [14,19,20]:
  $\delta\Delta G = -(\nu_{cb}\,\delta E_{F\!n} + \nu_{vb}\,\delta E_{F\!p})$. Since this is
  mathematically isomorphic to the electrochemical term, **photo and electro drivers compose
  directly into photoelectrocatalysis** with no additional code.
* **sono** — cavitation does not shift energies; it makes the operating point non-uniform. A second
  hook returns weighted sub-states, and the total rate is $r=\sum_i w_i r(\mathrm{state}_i)$ with
  $\sum_i w_i = 1$ enforced. The hook is not specific to cavitation: any distribution over operating
  conditions — a temperature profile along a bed, a distribution of local potentials — enters the
  same way.
* **field** — a static field shifts each step by its own dipole and polarizability change,
  $\delta\Delta G=-(\Delta\mu F+\tfrac12\Delta\alpha F^2)$, with $\Delta\mu$ and $\Delta\alpha$
  supplied per step from electronic-structure calculations. The quadratic term makes this driver
  non-additive, which the stack detects and reports when it is combined with another driver. It is
  included in the package and covered by the regression suite, but no application of it is reported
  in this paper.

For the photochemical case the splitting is obtained self-consistently from the carrier balance
$G = R_\mathrm{rec}(x) + R_\mathrm{surf}(x,\theta)$ with $x=\Delta E_\mathrm{qfl}/k_BT$ and

$$R_\mathrm{rec}(x) = k_\mathrm{bi}\,(e^{x}-1)
 + k_\mathrm{srh}\,\frac{e^{x}-1}{e^{fx}+e^{(1-f)x}} ,$$

both terms vanishing exactly in the dark. The splitting fraction $f$ between electron and hole
quasi-Fermi levels is a parameter: $f=1/2$ holds for intrinsic or high-injection conditions, whereas
in a doped semiconductor the majority-carrier level barely moves and the splitting is carried almost
entirely by the minority carrier. We show in Section 5.2 that the high-intensity power-law exponent
equals $\beta\nu f$, so fixing $f=1/2$ would give the wrong exponent for doped systems.

### 2.3 Gas-phase thermochemistry

Density-functional energies are electronic; gas-phase entropy is not included. At 673 K the entropy
term for H$_2$ alone is 1.08 eV, larger than most chemical effects. KinSTEP computes it from NIST
Shomate coefficients [21] using one rule,

$$ \delta\Delta G_s = -T\sum_j \nu_{sj} S_j(T), \qquad
   \delta G^\ddagger_s = +T\sum_j a^{r}_{sj} S_j(T), $$

the barrier correction counting only reactant-side gases because the transition state is adsorbed:
a gas molecule has already surrendered its translational and rotational freedom by the time it
reaches the transition state, whereas a gaseous product appears only in the final state. Each
species carries its own reference $S^\circ$ and $C_p^\circ$ at 298.15 K, and a self-check recomputes
them from the coefficients; a transcription error therefore cannot be released silently.
Extrapolation outside a Shomate temperature range raises an error rather than returning a plausible
but wrong number.

## 3. Numerical method

### 3.1 Unknowns and Jacobian

Surface coverages are solved in $\ln\theta$ space, with vacancies treated as additional unknowns and
site conservation imposed as an equation. Three reasons: coverages span twenty orders of magnitude;
$\theta=e^{x}>0$ holds automatically, so no bound constraints are needed; and mass action is
*linear* in log space, so the Jacobian is analytic,

$$ \frac{\partial r_s}{\partial x_j} = R_{f,s}\Big[a^{r}_{sj}(1-z_s) - z_s\big(a^{p}_{sj}-a^{r}_{sj}\big)\Big], \qquad
\frac{\partial F_i}{\partial x_j} = \sum_s \nu_{si}\frac{\partial r_s}{\partial x_j}. $$

An analytic Jacobian matters here beyond speed: a finite-difference Jacobian would inject noise of
order $10^{-8}$ into precisely the quantity the formulation is designed to preserve. Residuals and
Jacobians are assembled as NumPy arrays [22] and passed to the stiff integrators and least-squares
routines of SciPy [23]; figures are produced with Matplotlib [24].

For a continuous-flow (CSTR) calculation the gas partial pressures become unknowns as well, with
$\mathrm{d}p_i/\mathrm{d}t = \mathcal{P}_n[\sum_s\nu_{si}r_s + k_\mathrm{flow}(p_{\mathrm{in},i}-p_i)]=0$,
in which $k_\mathrm{flow}$ (s$^{-1}$) is the reciprocal residence time of the reactor and
$\mathcal{P}_n$ converts a per-site rate into a rate of change of partial pressure through the
site density and the gas volume.
Because $\ln p$ enters the rate law exactly as $\ln\theta$ does, the Jacobian expression above
applies unchanged to the gas columns.

### 3.2 The precision limit of $1-z$

Evaluating $1-z$ as $1-\exp(\ln z)$ loses all significance as $z\to1$ — the textbook case of
catastrophic cancellation [25,26]: at a true $1-z$ of
$10^{-12}$ the relative error is already $2.2\times10^{-5}$, and $8.0\times10^{-4}$ at $10^{-14}$.
Using $-\mathrm{expm1}(\ln z)$ removes this error entirely (measured relative error zero at both
values). The bottleneck is therefore **not** the evaluation of $1-z$ but the absolute accuracy of
$\ln z$ itself, which is a sum of terms of order ten and is consequently accurate to about
$10^{-16}$ in absolute terms. No reformulation in double precision can recover $1-z\sim10^{-30}$.

This observation has a second, more useful consequence. A quasi-equilibrated step pins the coverages
on an equilibrium relation — exact in log space — while the flux through it is set by the slow
steps. Its net rate should therefore be obtained from flux conservation, not from $R_f(1-z)$. We
implement this as a weighted least-squares projection onto $\{r:\nu^{\mathsf{T}}_\mathrm{surf}r=0\}$
with error bars $\sigma_s = R_{f,s}\,\varepsilon$, where $\varepsilon$ is the double-precision unit
roundoff and $R_{f,s}$ is the forward rate of the step, so that each step is weighted by the
absolute error of its own net rate. Using $|r_s|\varepsilon/(1-z_s)$ instead is
algebraically equivalent where $1-z>0$ but assigns maximum weight to a step carrying no information
when $1-z$ underflows; in a 401-point scan this produced rates uniformly low by a factor of 14 while
the flux *ratios* remained exact. With $\sigma=R_f\varepsilon$ the element-flux residual of the
ammonia test case falls from $1.55\times10^{-7}$ to $2.22\times10^{-16}$.

### 3.3 The ladder

| Level | Role | Why it is needed (measured) |
|---|---|---|
| (1) stiff ODE, adaptive continuation | physical initial guess | Newton from an empty surface diverges; the residual criterion replaces guessing an integration horizon. Estimating the horizon from $\min(k)$ fails by 18 orders of magnitude. |
| (2) damped Newton, analytic Jacobian | double-precision solution | A general-purpose solver does not know $\theta\le1$ and overflows $\exp$; step limiting plus backtracking is required. |
| (3) flux projection | net rates of quasi-equilibrated steps | Section 3.2. |
| (4) arbitrary precision (`mpmath`) | $1-z$ below the double-precision floor | Triggered automatically from $1-z$, the projection magnitude and the full residual (thresholds below). |

Level (4) is entered when any reversible step has $1-z\le0$ or $1-z<10^{-10}$ — at which point the
attainable relative precision of that step's net rate, $\varepsilon/(1-z)$, has already fallen to
$10^{-6}$ — or when the flux projection has to move the rates by more than $10^{-6}$ in relative
terms, or when the full residual, gas-phase block included, exceeds the convergence tolerance. The
conclusions do not depend on the particular value of the first threshold: moving it between
$10^{-8}$ and $10^{-12}$ changes how many points are escalated but leaves the final agreement with
the reference at the same level, as tabulated in SI Section S11.

The high-precision stage solves the *same* equation system with the *same* analytic Jacobian, so the
two paths can be cross-checked digit-by-digit; agreement is $1.6\times10^{-15}$ in $\theta$ and
$2.2\times10^{-15}$ in the net rates at points where double precision is sufficient. Level (4) raised
the 401-point maximum deviation from $5.599\times10^{-2}$ to $2.706\times10^{-10}$ — the storage
precision of the reference values themselves — while increasing the total runtime only from 1.7 s
to 4.2 s, because it is invoked at the 236 points where the criteria fire and starts each of them
from the double-precision solution.

### 3.4 Self-checks

Six consistency relations follow from the formulation and cost nothing to evaluate: the Campbell
sum rule $\sum_s X_{B,s}=1$ [27,28]; flux conservation $\nu^{\mathsf{T}}r=0$; site conservation;
thermodynamic consistency of every independent cycle (the Wegscheider conditions
$\sigma^{\mathsf{T}}\Delta G=0$ for each $\sigma$ in the left null space of $\nu$, which catch a
single mistyped reaction energy — a 0.1 eV change in one step is detected immediately and reported
together with the cycle that fails, although a mechanism whose steps are linearly independent has
no non-trivial cycle and the check then reports that it does not apply rather than that it passed);
equality of electron and hole consumption rates in photocatalysis; and
normalisation of multi-state weights. They are not decoration. The carrier-balance check, for example, immediately exposed a
mechanism in which a buffered intermediate decoupled the oxidation and reduction half-reactions,
giving electron and hole consumption rates differing by a factor of 90 — a *modelling* error that
produces a perfectly smooth and entirely wrong curve. The Campbell sum rule is additionally used to
select the finite-difference step size for the degree of rate control, since the attainable accuracy
of $\ln r$ varies by orders of magnitude between systems and a fixed step size must eventually fail.

## 4. Verification

Four independent references were used; all comparisons are digit-by-digit, and in each case the
constants of the reference implementation were matched explicitly before comparing, since otherwise
one compares constants rather than formulations. Five separate constant-convention differences from
CATKINAS were located and handled this way ($k_B$, $e$, $h$, $F$, and the temperature entering
$k_\mathrm{flow}$); the package's own `units` module was corrected twice in the process, having
mixed CODATA generations.

**Koper CO-stripping model [29].** The ODE right-hand side agrees to $2.44\times10^{-15}$ at random
$(t,\theta_\mathrm{OH},\theta_\mathrm{CO})$; peak currents for the three initial coverages of the
original figure agree to six digits. The residual difference in the integrated CV scales with the
integrator tolerance (a factor of 667 reduction from rtol $10^{-7}$ to $10^{-8}$), confirming it is
tolerance propagation rather than a formulation difference.

**Analytic closed form.** For formic-acid decomposition, $\theta$ and $r$ agree with
$\theta=k_1p/(k_1p+k_{-1}+k_2)$ to $2.22\times10^{-16}$, including in the quasi-equilibrium limit
where $1-z$ underflows.

**CATKINAS 401-point descriptor scan and 315-point flow-reactor map [11].** Starting only from the
four scaling/BEP coefficients of the published input and the entropy corrections computed from NIST
data, the full chain reproduces coverages to $1.41\times10^{-10}$, reversibilities to
$3.59\times10^{-11}$ and net rates to $2.71\times10^{-10}$; the latter equals the reference value's
own storage precision. All 401 points converge in 4.2 s. The flow-reactor map converges at all 315
points in 7.8 s with a median outlet-pressure deviation of $2.2\times10^{-16}$. Fig. 2 summarises.

Still against CATKINAS, the sensitivity measures agree as well: degrees of rate control to
$3.8\times10^{-12}$ (step), with the Campbell sum rule satisfied to $2.9\times10^{-12}$, and the
apparent activation energy to $3.2\times10^{-8}$ eV. The residual difference in the intermediate
DRC was shown to be pure second-order truncation by verifying that
$\text{deviation}/\delta^2$ is constant to four significant figures over a 32-fold range of $\delta$,
after which Richardson extrapolation reduced it by a factor of $2\times10^{5}$.

The reference values of these comparisons are archived with the deposit, so that the comparison can
be repeated without the reference software. What is archived is not that software and not its
example files, but the arrays this suite actually compares against — fourteen of the twenty-three
variables of the descriptor-scan output and seven of the twenty-three of the flow-reactor output,
0.2 MB in total. Where the original distribution is present the tests read it directly instead, and
one assertion then compares the two element by element, so the archived copy cannot silently drift
away from the distribution it was taken from.

The cost of the two precision strategies can be compared on exactly these examples, because the
reference distribution ships the run log of each. At its default setting of 1000 significant digits
the reference code reports 657 s for the solve stage of the 401-point scan, executed on eight
parallel workers (`npar = 8`), and 556 s for the 315-point flow-reactor map; the corresponding
KinSTEP times are 4.2 s and 7.8 s on a single thread. The machines are not the same — the reference
log is from the original 2020 run — so these are order-of-magnitude statements and not a benchmark.
What they bound is the cost of applying arbitrary precision everywhere instead of at the points
where a measured criterion requires it: here that criterion fires at 236 of the 401 points and
raises the total from 1.7 s to 4.2 s.

**CatMAP CO-oxidation map [5].** The three references above share a limitation: none of them is an
open-source package of comparable scope. We therefore added a fourth, on a different chemistry and
against an independent implementation. CatMAP was run unmodified on its own CO-oxidation tutorial
input (225 descriptor points, 500 K, its default 1000-digit working precision), and its rate
constants, coverages and step rates were exported. KinSTEP was then given the *same* rate constants
through its `mode="given"` interface — no energies, no scaling, so that only the rate law and the
steady-state solution are being compared. Conventions were checked first: CatMAP's own rates are
reproduced to the last digit from its own constants and coverages by $r=k_f\prod a - k_r\prod a$ at
all points tested, so what follows is not a convention difference.

Coverages agree to a median of $3.1\times10^{-31}$ over the 193 points at which both codes select
the same steady state. For the step rates the comparison has to be conditioned, and the reason is
itself the subject of Section 3.2. At 91 of those 193 points CatMAP's three step rates satisfy the
mass balance implied by its own mechanism ($r_1=r_3$, $2r_2=r_3$) to better than $10^{-10}$; there
the two codes agree to a median of $2.2\times10^{-15}$, $2.4\times10^{-15}$ and
$2.4\times10^{-15}$ for the three steps — machine precision. At the remaining points the
discrepancy is confined to the quasi-equilibrated CO adsorption step, whose net rate is a
difference of two large numbers: there the error comes from the residual convergence of the
coverages and is therefore not removed by carrying more digits in the arithmetic. The flux-projected
value is mass-balance consistent by construction — KinSTEP's residual
$\max(|r_1-r_3|,|2r_2-r_3|)$ is exactly zero at every one of the 193 points — which is the
practical form of the argument of Section 3.2. The comparison is reproduced by
`paper/compare_catmap.py` from the exported values deposited with the Supplementary Information.

## 5. Example applications

### 5.1 Ammonia synthesis: descriptor scan (thermal)

Starting from the five-step mechanism and scaling/BEP coefficients, a 401-point scan over the
nitrogen binding energy produces the volcano of Fig. 2(a), with the peak at
$E_\mathrm{N}=-0.58$ eV and rates spanning 28 orders of magnitude. Fig. 3(d) shows the step degree
of rate control [27,28] along the same descriptor. Control changes hands between
$E_\mathrm{N}=-0.68$ and $-0.56$ eV, that is, essentially at the volcano maximum: on the
strong-binding branch the final hydrogenation step NH$_2^*+$H$^*\to$NH$_3$ (step 5) is
rate-controlling, whereas on the weak-binding branch control transfers to dissociative H$_2$
adsorption (step 2). The components are drawn as lines rather than stacked, because $X_{B,2}$ is
*negative* on the strong-binding branch, reaching $-0.36$ at $E_\mathrm{N}=-0.68$ eV where
$X_{B,5}=1.36$: lowering the H$_2$ adsorption barrier there slows ammonia production, the sites
gained by H$^*$ being taken from the nitrogen-removal path. A stacked area would have to clip that
component, and would then display a total of 1.36 beside an annotation asserting
$\sum_s X_{B,s}=1$. The sum rule is in fact satisfied to $7.7\times10^{-7}$ at every point of the
scan, negative components included, and is used as the acceptance criterion for the perturbed
solutions rather than as a post-hoc check.

### 5.2 Hydrogen evolution: Tafel analysis (electro)

For the Volmer–Heyrovsky–Tafel mechanism [30,31] the computed Tafel slope in the strongly polarised
region is $118.3$ mV dec$^{-1}$ (Fig. 3a), which is exactly $\ln\!10\,k_BT/(\beta e)$ at $\beta=0.5$
and 298.15 K. This value is not an input: it emerges from the mechanism, the energies and the CHE/BV
relations. We describe it as an end-to-end consistency check rather than an independent validation,
since the analytic limit follows from the same Butler–Volmer form that the driver implements; what
it verifies is that the whole chain — electron counting, potential referencing, the barrier
correction and the steady-state solution — reproduces that limit without any value having been
tuned to it. The symmetry
factor is a driver parameter rather than a fitted constant; first-principles treatments of the
potential dependence of electrochemical barriers give values near $0.5$ for Volmer-type steps [18].
The same CHE construction underlies published oxygen-reduction analyses [32].

### 5.3 Photocatalytic water splitting: light-intensity law (photo)

Semiconductor photocatalysis [33,34] depends characteristically on photon flux. Scanning the flux
over thirteen decades (Fig. 3b), the local power-law exponent moves from
$0.996$ at low intensity to $0.5000$ at high intensity, in agreement with the classical
transition from first- to half-order [35,36]. The exponent is a *result*, not an input, and the
apparent quantum yield falls from 1.000 to $4.5\times10^{-5}$ across the same range, identifying
recombination as the origin of the half-order regime. The quantity is reported as an *apparent*
yield throughout, since converting it to a true quantum efficiency requires the absorbed rather than
the incident photon flux [37]. Two conditions emerge from the model that are
not usually stated explicitly: the half-order law requires the rate-determining step to transfer
**one** carrier (a step consuming one electron and one hole experiences the full
$\Delta E_\mathrm{qfl}$ and remains first order), and the exponent is $\beta\nu f$, so a doped
semiconductor with $f\ne 1/2$ follows a different law. The analytically predicted crossover flux
agrees with the numerical result.

### 5.4 Cavitation hot spots (sono)

With a hot spot at 5000 K — the value established spectroscopically for multibubble cavitation
[38,39], single-bubble collapse reaching higher still [40] — occupying a fraction $10^{-3}$ of the
system and a bulk at 350 K, the
weighted rate $\sum_i w_i r_i$ exceeds the rate evaluated at the mean temperature of 354.65 K by
$1.39\times10^{23}$ (Fig. 3c). The crossover fraction at which the hot spot begins to dominate is
$f^{*}=r_\mathrm{bulk}/r_\mathrm{hot}=1.79\times10^{-27}$ for a 1.95 eV barrier, which is why an
extremely small cavitating volume produces a measurable chemical effect, and the apparent activation
energy falls from 1.46 eV to zero as $f$ increases — an experimentally accessible signature of
cavitation control. KinSTEP does not fabricate a relation between acoustic intensity or frequency
and hot-spot conditions; the literature is not consistent on the frequency dependence, and
single-bubble yields require full bubble dynamics. Hot-spot conditions are inputs, and a helper
converts a measured sonochemical yield into the hot-spot fraction following established dosimetry
practice [38].

### 5.5 Two-dimensional operating maps

Each example above fixes one condition and varies another. Because the same core serves all four
modalities, a two-dimensional map over the operating plane costs no more than a line scan, and it is
in such maps that the qualitative behaviour becomes visible. Fig. 4 shows one map per modality,
4992 solutions in total, every one of which converged.

For ammonia synthesis over the plane $T=500$–900 K, $P=1$–200 bar (1728 points, Fig. 4a) the net
rate is positive at 60.2% of the points; the explicit $r=0$ contour, not a colour-scale boundary,
separates synthesis from decomposition, and the highest net rate, $0.206$ s$^{-1}$ site$^{-1}$,
occurs at the hot, high-pressure corner (857 K, 200 bar). The accompanying coverage panels (a2, a3) show
why: $\theta_{\mathrm{H}^*}$ and the vacancy fraction trade against each other across the plane,
with $\theta_*$ ranging from $4.1\times10^{-5}$ to $0.987$.

For hydrogen evolution over $U=-0.40$ to $-0.02$ V vs SHE and pH 0–14 (Fig. 4b) the current density
spans thirty orders of magnitude, from $1.1\times10^{-23}$ to $4.6\times10^{6}$ A m$^{-2}$ at the
assumed active-site areal density of $2.2\times10^{-5}$ mol m$^{-2}$, which converts the
per-site rate to a current density; the upper end is far beyond what mass transport would sustain
and is reported as the kinetic prediction of the model, not as an attainable current. The Tafel panel (b3) is the more informative one. Over the 1394 points at which the slope is
defined and below 400 mV dec$^{-1}$ (two further points exceed that bound and are clipped in the
display) its magnitude has median 46.5 mV dec$^{-1}$ with interquartile range 31.0–118.3, and the
distribution is bimodal: 49.2% of those points fall in the 25–45 mV dec$^{-1}$ band characteristic
of Heyrovsky or Tafel control and 34.1% in the 100–140 mV dec$^{-1}$ band characteristic of Volmer
control. The rate-determining step therefore *changes across the operating plane*, which a Tafel
slope measured at one pH cannot reveal and which the map makes immediate.

For photocatalysis over twelve decades of photon flux and four decades of recombination coefficient
(Fig. 4c) the local power-law order in flux — computed as a local logarithmic derivative, not
assumed — spans 0.5000 to 0.9999 across the plane, averaging 0.932 over the three lowest-flux
columns and 0.500 over the three highest, so the entire first-order-to-half-order transition is
recovered rather than the two limits alone; the quasi-Fermi-level splitting varies from 0.473 to
1.460 eV over the same plane.

For cavitation over hot-spot fraction $10^{-5}$–$10^{-2}$ and hot-spot temperature 1500–6000 K
(Fig. 4d) the ratio between the correctly weighted rate and the rate evaluated at the mean condition
spans 16.1 to 23.9 orders of magnitude. More seriously, at 3 of the 721 points where the comparison
is defined that ratio is *negative*: evaluating the model at the mean condition returns not merely
the wrong magnitude but the wrong direction of reaction. The remaining 47 cells are hatched because
the rate at the mean condition underflows to zero, so the ratio is undefined there; they are not
solver failures. Hatching marks an undefined quantity throughout Fig. 4, and is used in preference
to a flat colour because the centre of a diverging scale is itself a pale grey — every one of the
4992 steady states converged.

## 6. Conclusions

KinSTEP demonstrates that the numerical core of microkinetic modelling is genuinely independent of
catalytic modality, and that isolating the stimulus behind a small plug-in interface is enough to
cover thermal, electrochemical, photochemical and sonochemical catalysis with a single verified
solver. The arrangement also composes: photoelectrocatalysis requires no code beyond the two
existing drivers.

Two numerical results are of broader interest. First, the precision limit of quasi-equilibrated
steps resides in the absolute accuracy of $\ln z$ rather than in the evaluation of $1-z$, which
means no double-precision reformulation can recover the deep quasi-equilibrium regime; an automatic,
measured escalation to arbitrary precision is both sufficient and cheap. Second, the net rate of a
quasi-equilibrated step is better recovered from flux conservation than computed directly, with the
correct error weighting being $\sigma=R_f\varepsilon$.

One limitation of the ladder is worth stating because the CatMAP comparison exposed it. At the
CO-saturated corner of that map the coverage of the free site falls to $10^{-31}$, and the
arbitrary-precision stage fails to converge there at 50, 80 or 150 digits, so the solver falls back
to the double-precision solution and says so in its notes. The correct net rates at those points
come from level (3), the flux projection, not from level (4). The ladder therefore degrades
gracefully rather than silently, but its top level is not universal.

Current limitations: lateral interactions are not yet implemented (the interface accommodates them
through a coverage-dependent driver); the arbitrary-precision backend is a scalar implementation and
is therefore slow for very large mechanisms; the gas-phase thermochemistry database currently
contains ten common species, to which further species are added by supplying their Shomate
coefficients with a literature reference; and the model is purely kinetic — mass transport to and
from the surface, the diffusion layer of an electrode and light attenuation through a suspension
are all outside it, so the rates reported here are the kinetically attainable ones and must be
compared with experiment only where transport is not limiting.

### 6.1 Using the program

A microkinetic study in KinSTEP is four commands, shown in Fig. 5 exactly as they appear on the
terminal. `kinstep new --template her` writes a skeleton in which the mechanism is already filled in
and every blank energy slot carries a sentence saying what to compute (Fig. 5a); the number of
blanks it reports is the number the validator will ask for. `kinstep check` then runs the structural
gates and, when they pass, hands the file to the provenance gate, which refuses the model while any
energy is missing or untrustworthy and exits with code 3 after listing each slot and what is needed
(Fig. 5b). `kinstep ask` converts one sentence of plain English into a structure plus a list of
questions, by rule and with no language model involved: quantities the user stated are extracted and
tagged `user_given`, whereas every energy becomes a question rather than an assumption (Fig. 5c).
Finally `kinstep run --drc` solves the model and reports coverages, the attainable precision of each
step, degrees of rate control and the Campbell sum rule in one pass (Fig. 5d). Exit codes are
distinct — 0, 2, 3 and 4 for success, malformed input, a provenance refusal and a driver
misconfiguration — so the commands compose in scripts.

All user-facing output is available in English by setting `KINSTEP_LANG=en`, which is how Fig. 5
was produced; one regression asserts that no CJK character appears anywhere in English-mode output,
including in the files the program writes. The complete terminal session behind Fig. 5, of which
the figure shows a selection of lines, is deposited with the Supplementary Information.

### 6.2 Code availability and reproducibility

KinSTEP is released under GPL-3.0-or-later. The deposit contains the package, the regression suite,
the six worked examples, and the scripts and raw data behind every figure in this paper and its
Supplementary Information.

Two points about reproducibility are stated explicitly rather than left for the reader to discover.
First, the regression suite consists of 82 assertions, **all of which run on a clean download**.
They cover the analytic closed-form comparisons, the published CO-stripping model [29], the
closed-form carrier balance of the photocatalytic driver, the closed-form multi-state average of the
sonochemical driver, the field driver, all command-line behaviour, the six examples and the error
paths, together with the comparisons against the reference values of Section 4. Nothing here
requires third-party software: the reference *arrays* those comparisons need are archived with the
deposit, while the software that produced them is not redistributed. Installed from a freshly built
source archive into an empty virtual environment with only NumPy and SciPy, 77 of the 82 pass and
five are reported as *skipped*, each with its reason: four compare against values that only the
optional `mpmath` extra can resolve below the double-precision floor, and one exports a figure and
therefore needs the optional `matplotlib` extra. With both extras installed, all 82 pass. All three
configurations were measured before submission rather than inferred.

Second, every figure in this paper regenerates from the deposit alone, Fig. 2 and panels 3(d) and
S1(b) included: they are recomputed from the archived reference arrays rather than redrawn from
stored coordinates. No figure contains a number that is not either computed by the deposited code or
deposited alongside it.

## CRediT authorship contribution statement

**Xin Lei:** Conceptualization, Methodology, Software, Validation, Formal analysis, Investigation,
Data curation, Writing – original draft, Writing – review & editing, Visualization.

## Declaration of competing interest

The author declares no known competing financial interests or personal relationships that could have
appeared to influence the work reported in this paper.

## Data availability

All data are available with the program. The deposited archive contains the source code, the
regression suite, the six worked examples, and the scripts and raw data behind every figure in this
paper and in the Supplementary Information, so that every figure and every assertion can be
reproduced from the deposit alone. This includes the reference arrays of Section 4: CATKINAS [11]
is third-party closed-source software whose licence does not permit redistribution, so the software
and its example files are not included, but the numerical values this work compares against — 0.2 MB
of arrays, listed in the Supplementary Information — are archived with the deposit, as the values
compared against in a paper normally are. Where a user has the original distribution, the tests read
it directly and one assertion compares the two element by element.

The version described in this paper, v0.1.0, is archived at
https://doi.org/10.5281/zenodo.23613004; the concept DOI
https://doi.org/10.5281/zenodo.23613003 always resolves to the most recent release.

## Funding

This research did not receive any specific grant from funding agencies in the public, commercial, or
not-for-profit sectors.

## References

[1] J.A. Dumesic, D.F. Rudd, L.M. Aparicio, J.E. Rekoske, A.A. Treviño, *The Microkinetics of
Heterogeneous Catalysis*, American Chemical Society, Washington, DC, 1993.

[2] P. Stoltze, Microkinetic simulation of catalytic reactions, *Prog. Surf. Sci.* **65** (2000)
65–150. doi:10.1016/S0079-6816(00)00019-8

[3] A.H. Motagamwala, J.A. Dumesic, Microkinetic modeling: a tool for rational catalyst design,
*Chem. Rev.* **121** (2021) 1049–1076. doi:10.1021/acs.chemrev.0c00394

[4] L.C. Grabow, Computational catalyst screening, in: *Computational Catalysis*, Royal Society of
Chemistry, Cambridge, 2013, pp. 1–58. doi:10.1039/9781849734905-00001

[5] A.J. Medford, C. Shi, M.J. Hoffmann, A.C. Lausche, S.R. Fitzgibbon, T. Bligaard, J.K. Nørskov,
CatMAP: a software package for descriptor-based microkinetic mapping of catalytic trends,
*Catal. Lett.* **145** (2015) 794–807. doi:10.1007/s10562-015-1495-6

[6] F. Abild-Pedersen, J. Greeley, F. Studt, J. Rossmeisl, T.R. Munter, P.G. Moses, E. Skúlason,
T. Bligaard, J.K. Nørskov, Scaling properties of adsorption energies for hydrogen-containing
molecules on transition-metal surfaces, *Phys. Rev. Lett.* **99** (2007) 016105.
doi:10.1103/PhysRevLett.99.016105

[7] T. Bligaard, J.K. Nørskov, S. Dahl, J. Matthiesen, C.H. Christensen, J. Sehested, The
Brønsted–Evans–Polanyi relation and the volcano curve in heterogeneous catalysis,
*J. Catal.* **224** (2004) 206–217. doi:10.1016/j.jcat.2004.02.034

[8] S. Vijay, H.H. Heenen, A.R. Singh, K. Chan, J. Voss, Number of sites-based solver for
determining coverages from steady-state mean-field micro-kinetic models,
*J. Comput. Chem.* **45** (2024) 546–551. doi:10.1002/jcc.27263

[9] I.A.W. Filot, R.A. van Santen, E.J.M. Hensen, The optimally performing Fischer–Tropsch catalyst,
*Angew. Chem. Int. Ed.* **53** (2014) 12746–12750. doi:10.1002/anie.201406521

[10] B. Medasani, S. Kasiraju, D.G. Vlachos, OpenMKM: an open-source C++ multiscale modeling
simulator for homogeneous and heterogeneous catalytic reactions, *J. Chem. Inf. Model.*
**63** (2023) 3377–3391. doi:10.1021/acs.jcim.3c00088

[11] J. Chen, H. Wang, Y. Hu, CATKINAS: a large-scale catalytic microkinetic analysis software for
mechanism auto-analysis and catalyst screening, *J. Comput. Chem.* **42** (2021) 379–391.
doi:10.1002/jcc.26464

[12] M. Stamatakis, D.G. Vlachos, A graph-theoretical kinetic Monte Carlo framework for on-lattice
chemical kinetics, *J. Chem. Phys.* **134** (2011) 214115. doi:10.1063/1.3596751

[13] M. Salciccioli, M. Stamatakis, S. Caratzoulas, D.G. Vlachos, A review of multiscale modeling of
metal-catalyzed reactions: mechanism development for complexity and emergent behavior,
*Chem. Eng. Sci.* **66** (2011) 4319–4355. doi:10.1016/j.ces.2011.05.050

[14] M. Monai, W. Albrecht, A. Alkemper, N. Artrith, A. Baldi, et al., Grand challenges and
opportunities in stimulated dynamic and resonant catalysis, *ACS Catal.* **16** (2026) 4077–4112.
doi:10.1021/acscatal.5c07014

[15] H. Eyring, The activated complex in chemical reactions, *J. Chem. Phys.* **3** (1935) 107–115.
doi:10.1063/1.1749604

[16] J. Chen, Y. Mao, H. Wang, P. Hu, Reversibility iteration method for understanding reaction
networks and for solving microkinetics in heterogeneous catalysis, *ACS Catal.* **6** (2016)
7078–7087. doi:10.1021/acscatal.6b02405

[17] J.K. Nørskov, J. Rossmeisl, A. Logadottir, L. Lindqvist, J.R. Kitchin, T. Bligaard, H. Jónsson,
Origin of the overpotential for oxygen reduction at a fuel-cell cathode,
*J. Phys. Chem. B* **108** (2004) 17886–17892. doi:10.1021/jp047349j

[18] K. Chan, J.K. Nørskov, Electrochemical barriers made simple,
*J. Phys. Chem. Lett.* **6** (2015) 2663–2668. doi:10.1021/acs.jpclett.5b01043

[19] S. Suo, C. Sheehan, F. Zhao, L. Xiao, Z. Xu, J. Meng, T.E. Mallouk, T. Lian, Direct vibrational
Stark shift probe of quasi-Fermi level alignment in metal nanoparticle catalyst-based
metal–insulator–semiconductor junction photoelectrodes, *J. Am. Chem. Soc.* **145** (2023)
14260–14266. doi:10.1021/jacs.3c02333

[20] Y. Shioiri, K. Obata, Y. Kawase, T. Higashi, M. Katayama, M. Schleuning, R. van de Krol,
D. Friedrich, F.F. Abdi, K. Takanabe, Estimating the quasi-Fermi level of holes at the surface of
semiconductor photoanodes using outer-sphere redox couples, *Nat. Commun.* **16** (2025) 3688.
doi:10.1038/s41467-025-58837-9

[21] M.W. Chase, Jr., *NIST-JANAF Thermochemical Tables*, 4th ed., *J. Phys. Chem. Ref. Data*
Monograph 9 (1998) 1–1951; NIST Chemistry WebBook, SRD 69.

[22] C.R. Harris, K.J. Millman, S.J. van der Walt, R. Gommers, et al., Array programming with NumPy,
*Nature* **585** (2020) 357–362. doi:10.1038/s41586-020-2649-2

[23] P. Virtanen, R. Gommers, T.E. Oliphant, M. Haberland, et al., SciPy 1.0: fundamental algorithms
for scientific computing in Python, *Nat. Methods* **17** (2020) 261–272.
doi:10.1038/s41592-019-0686-2

[24] J.D. Hunter, Matplotlib: a 2D graphics environment, *Comput. Sci. Eng.* **9** (2007) 90–95.
doi:10.1109/MCSE.2007.55

[25] D. Goldberg, What every computer scientist should know about floating-point arithmetic,
*ACM Comput. Surv.* **23** (1991) 5–48. doi:10.1145/103162.103163

[26] N.J. Higham, *Accuracy and Stability of Numerical Algorithms*, 2nd ed., SIAM, Philadelphia,
2002. doi:10.1137/1.9780898718027

[27] C. Stegelmann, A. Andreasen, C.T. Campbell, Degree of rate control: how much the energies of
intermediates and transition states control rates, *J. Am. Chem. Soc.* **131** (2009) 8077–8082.
doi:10.1021/ja9000097

[28] C.T. Campbell, The degree of rate control: a powerful tool for catalysis research,
*ACS Catal.* **7** (2017) 2770–2779. doi:10.1021/acscatal.7b00115

[29] M.T.M. Koper, A.P.J. Jansen, R.A. van Santen, J.J. Lukkien, P.A.J. Hilbers, Monte Carlo
simulations of a simple model for the electrocatalytic CO oxidation on platinum,
*J. Chem. Phys.* **109** (1998) 6051–6062. doi:10.1063/1.477230

[30] J.K. Nørskov, T. Bligaard, A. Logadottir, J.R. Kitchin, J.G. Chen, S. Pandelov, U. Stimming,
Trends in the exchange current for hydrogen evolution,
*J. Electrochem. Soc.* **152** (2005) J23–J26. doi:10.1149/1.1856988

[31] E. Skúlason, V. Tripkovic, M.E. Björketun, S. Gudmundsdóttir, G. Karlberg, J. Rossmeisl,
T. Bligaard, H. Jónsson, J.K. Nørskov, Modeling the electrochemical hydrogen oxidation and
evolution reactions on the basis of density functional theory calculations,
*J. Phys. Chem. C* **114** (2010) 18182–18197. doi:10.1021/jp1048887

[32] H.A. Hansen, V. Viswanathan, J.K. Nørskov, Unifying kinetic and thermodynamic analysis of
2 e$^-$ and 4 e$^-$ reduction of oxygen on metal surfaces,
*J. Phys. Chem. C* **118** (2014) 6706–6718. doi:10.1021/jp4100608

[33] A. Fujishima, K. Honda, Electrochemical photolysis of water at a semiconductor electrode,
*Nature* **238** (1972) 37–38. doi:10.1038/238037a0

[34] M.R. Hoffmann, S.T. Martin, W. Choi, D.W. Bahnemann, Environmental applications of
semiconductor photocatalysis, *Chem. Rev.* **95** (1995) 69–96. doi:10.1021/cr00033a004

[35] C.S. Turchi, D.F. Ollis, Photocatalytic degradation of organic water contaminants: mechanisms
involving hydroxyl radical attack, *J. Catal.* **122** (1990) 178–192.
doi:10.1016/0021-9517(90)90269-P

[36] S. Tabata, H. Ohnishi, E. Yagasaki, M. Ippommatsu, K. Domen, Light-intensity dependence in
photocatalytic decomposition of water over K$_4$Nb$_6$O$_{17}$ catalyst,
*Catal. Lett.* **28** (1994) 417–422. doi:10.1007/BF00806072

[37] B. Ohtani, Photocatalysis A to Z — what we know and what we do not know in a scientific sense,
*J. Photochem. Photobiol. C* **11** (2010) 157–178. doi:10.1016/j.jphotochemrev.2011.02.001

[38] K.S. Suslick, Y. Didenko, M.M. Fang, T. Hyeon, K.J. Kolbeck, W.B. McNamara III, M.M. Mdleleni,
M. Wong, Acoustic cavitation and its chemical consequences,
*Phil. Trans. R. Soc. A* **357** (1999) 335–353. doi:10.1098/rsta.1999.0330

[39] E.B. Flint, K.S. Suslick, The temperature of cavitation,
*Science* **253** (1991) 1397–1399. doi:10.1126/science.253.5026.1397

[40] Y.T. Didenko, K.S. Suslick, The energy efficiency of formation of photons, radicals and ions
during single-bubble cavitation, *Nature* **418** (2002) 394–397. doi:10.1038/nature00895
