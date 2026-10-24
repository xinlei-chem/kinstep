# KinSTEP examples

Six **runnable** examples covering the four catalytic modalities and every core capability. They
live in `examples/` and are run directly, e.g. `python examples/ex1_ammonia_volcano.py`, after
`pip install -e ".[all]"`.

Every number below is **output from an actual run**, not an illustration. A Chinese version of this
file is kept as `EXAMPLES.zh-CN.md`.

---

## Example 1: the ammonia-synthesis volcano (thermal)

`examples/ex1_ammonia_volcano.py` — the complete descriptor-scan chain

```
mechanism parsing -> NIST entropy correction -> scaling + BEP -> 1D scan -> volcano plot + raw data
```

**What it shows**

- the mechanism is written as it is written on paper, and the four gates run by themselves;
- **the entropy correction is computed automatically from the NIST Shomate data**: at 673 K the H₂
  term alone is **1.08 eV**, larger than most chemical effects, so it cannot be left out;
- the scaling relation `ΔG⁰ = a·E_N + b` and the BEP relation `Ga = α·ΔG⁰ + β` are written as two
  separate layers, so it is visible which part is chemistry and which part is thermochemistry.

**Actual output**

```
entropy correction (eV): dGa = [1.5043 1.0773 0. 0. 0.]  dG0 = [1.5043 1.0773 0. 0. -1.576]
scan: E_N[201], 201 points, 2.34 s (11.7 ms/point)
  order: 1D continuation outward from the centre (100); converged 201/201; escalated to high precision at 118 points
  rate maximum: E_N = -0.580, r = 1.041326e-03
  the rate spans 28.3 orders of magnitude
```

**Worth noting:** **118 of the 201 points escalated automatically to arbitrary precision** — on the
strong-binding branch the net rate cannot be obtained in double precision at all. Those points are
drawn in red, which **makes the region where double precision fails visible on the figure itself**.

Products: `ex1_volcano.png`, `ex1_coverage.png`, `ex1_data.dat` (read directly by Origin or
gnuplot).

---

## Example 2: the Tafel slope of the HER (electro)

`examples/ex2_her_electro.py`

**What it shows**

- the `electro` driver: `ΔG(U) = ΔG(0) − ν_e·e·U` (CHE) and `Ga(U) = Ga(0) − β·ν_e·e·U`
  (Butler–Volmer);
- **`nu_elec` is counted from the equation** (electrons on the product side minus electrons on the
  reactant side), which is the parameter whose sign is easiest to get wrong;
- `H+` appears explicitly in the mechanism, so the potential **must be on the SHE scale** — using
  RHE would count the pH twice.

**Actual output**

```
U(V vs SHE)   j (A/m^2)       Tafel (mV/dec)
-0.300        -6.5064e+05     -118.3
-0.200        -9.2936e+04     -118.3
-0.100        -1.3274e+04     -118.3

Tafel slope in the strongly polarised region ~ -118.3 mV/dec
```

**Why that number matters:** 118.3 mV/dec is exactly `ln10·k_BT/(β·n·e)` for β = 0.5, n = 1 and
298 K. It is **not an input**; it comes out of the mechanism, the energies and the CHE/BV relations,
which makes it an independent check of the whole electrochemical chain.

Experimentally the Tafel slope distinguishes the rate-determining step: Volmer about 120,
Heyrovsky about 40, Tafel about 30 mV/dec.

---

## Example 3: the light-intensity power law in photocatalysis (photo)

`examples/ex3_photo_intensity.py`

**What it shows**

The quasi-Fermi-level splitting `ΔE_qfl` acts as the driving force and is solved self-consistently
from the carrier balance:

```
G = R_rec(ΔE_qfl) + R_surf(ΔE_qfl, θ)
```

**Actual output**

```
photon flux (m^-2 s^-1)   dE_qfl(eV)   rate(s^-1)     AQY         local exponent
1.00e+13                  0.3696       1.6666e-07     1.000e+00   0.9998
1.00e+17                  0.8300       1.2272e-03     7.363e-01   0.7893
1.00e+21                  1.1022       2.3728e-01     1.424e-02   0.5044
1.00e+25                  1.3407       2.3895e+01     1.434e-04   0.5000

low-flux exponent  0.9963  (literature: first order)
high-flux exponent 0.5000  (literature: half order, Turchi & Ollis 1990)
AQY falls from 1.000 to 4.53e-05 - at high flux photons are wasted to recombination,
which is the physical origin of the half-order regime
charge-balance self-check: largest imbalance 7.51e-16
```

**Why this example makes the case for MKM:** the empirical practice is to fit a power law
`rate ∝ I^β`; here **β is a computed result**, both limits emerge on their own, and the transition
between them is continuous.

The half-order law carries an implicit condition: the rate-determining step must transfer **one**
carrier. If it consumes one electron and one hole at the same time, the driving force is the full
ΔE_qfl and the exponent is 1, with no half-order regime at all. (This is not stated explicitly in
the literature; it was established by measurement in this project.)

A second caveat: the default `qfl_frac_n = 0.5` holds only for intrinsic material or high
injection. In a doped semiconductor almost all of the splitting falls on the minority-carrier side,
and the measured high-flux exponent is (1−f): f = 0.5 gives 0.5027 and f = 0.8 gives 0.2000, digit
for digit with the analytic prediction.

---

## Example 4: cavitation hot spots in sonocatalysis (sono)

`examples/ex4_sono_cavitation.py`

**In one line: `r(⟨T⟩) ≠ ⟨r(T)⟩`.**

**Actual output**

```
state    T(K)      weight      rate of that state   contribution   share
 1        5000.0   0.001000    1.238233e+11         1.238233e+08   100.000%
 2         350.0   0.999000    2.220446e-16         2.218226e-16     0.000%

<r(T)>, weighted state by state = 1.238233e+08   <- the correct procedure
r(<T>), at the mean temperature 354.65 K = 8.881784e-16   <- the common mistake
=> a factor of 1.394e+23 (23.1 orders of magnitude)

crossover fraction f* = r_bulk/r_hot = 1.793e-27
```

**Two results worth remembering**

1. **Using the mean temperature is wrong by 23 orders of magnitude.** The rate is exponentially
   sensitive to temperature, so a non-uniform system must be solved state by state and then
   averaged.
2. **f\* = 1.79e-27** — for a barrier as high as 1.95 eV, a cavitation fraction of only 10⁻²⁷ is
   already enough for the hot spot to dominate. This is why sonochemistry can drive reactions that
   do not proceed at all at ambient temperature.

**Where the hot-spot conditions come from:** this package **does not invent** a formula mapping
acoustic intensity and frequency onto hot-spot conditions. (The literature disagrees about the
effect of frequency, and the single-bubble yield requires full bubble dynamics.) The conditions are
supplied by the user from bubble dynamics or from **a measured dosimeter**, and
`calibrate_fraction_from_yield()` implements the literature route as a function; the measured
round-trip error of that inversion is **0.00e+00**.

Hot-spot temperatures differ greatly between methods: single-bubble sonoluminescence gives
5000–20000 K, multi-bubble systems 750–6000 K. Values outside that range raise a warning.

---

## Example 5: degrees of rate control, reaction orders, apparent activation energy

`examples/ex5_drc_orders.py`

Four kinds of DRC share one definition: `X = ∂ln r / ∂(−G_x/k_BT)`.

**Actual output**

```
step 1  N2 + 2* <-> 2N*       X_B = +0.999998     <- rate-determining
step 2  H2 + 2* <-> 2H*       X_B = +0.000000
sum rule: Sigma X_B = 1.000000000  (deviation 2.0e-13)

intermediate DRC: H* -0.439470   * -1.555504   sum = -2.000000
reaction orders:  H2: -0.2040   N2: +1.0044   NH3: -0.0138
apparent activation energy Ea_app = 3.016460 eV
```

**How to read it**

- the step with the largest `X_B` is rate-determining: to speed the reaction up, change that
  barrier;
- `Σ X_B = 1` (Campbell) **holds independently of the mechanism and of the energies**, so if it
  fails the result should not be trusted. This package also uses it to **choose the
  finite-difference step automatically**: a fixed δ will fail on some system sooner or later,
  whereas this rule is a free yardstick;
- `X_I` is the thermodynamic DRC of an intermediate: **a negative value means that stabilising that
  intermediate slows the reaction down**;
- **reaction orders and apparent activation energies are directly measurable**, which makes them
  the cheapest route to validation.

Also: `richardson=True` removes the second-order truncation error by Richardson extrapolation,
which improved the measured accuracy of `X_I` by a factor of **2.3 × 10⁵**.

---

## Example 6: a continuous-flow reactor (CSTR)

`examples/ex6_flow_reactor.py`

A differential reactor assumes the gas-phase activities are constant, which holds at very low
conversion. To obtain **how conversion varies with temperature and pressure**, the gas composition
must itself become an unknown.

**Actual output**

```
feed H2/N2/NH3 = 75.0 / 25.0 / 0.0 bar, T = 900 K

differential reactor (gas phase fixed): NH3 formation rate 3.0271e-03 /site/s, conversion identically 0

flow reactor (CSTR): converged
   outlet H2/N2/NH3 = 74.9830 / 24.9943 / 0.011323 bar
   H2 conversion = 0.0226%
   solver path: ode (3-stage continuation) -> flow initial guess (differential reactor
                + 2 gas-phase self-consistency steps) -> damped Newton (5 steps)
                -> high precision (50 digits)

temperature scan (continuation along T):
   900 K  NH3 = 1.132320e-02 bar   conversion 0.0226%
   700 K  NH3 = 1.473116e-07 bar
   500 K  NH3 = 1.907026e-16 bar
```

**Two points**

1. **The solver-path line is worth reading**: all four levels of the ladder plus the flow-specific
   initial guess were used. Starting directly from "gas phase = feed" falls into a **spurious
   steady state** (the product comes out as 1e-37 where the true value is 1.93, **and it reports
   convergence**), so the differential reactor is solved first to obtain a physical surface state,
   after which the gas phase is fixed algebraically by the mass balance.
2. **Continuation along temperature is necessary, not an optimisation**: starting each point from
   scratch does not converge at all at 300–350 K, whereas with continuation **315/315 points
   converge and the time falls from 50 s to 7.8 s**.

---

## Where to start

```bash
kinstep templates              # 11 built-in systems
kinstep new --template her     # a skeleton with the mechanism filled in and each energy slot saying what to compute
kinstep check her.kinstep.toml # tells you, slot by slot, what is still missing
# (fill in the DFT energies of your own system)
kinstep run her.kinstep.toml --drc
```

Or start from one sentence:

```bash
kinstep ask "oxygen evolution on IrO2 at 1.5 V vs SHE, 25 degrees C" --write
```

This produces **a structure plus a list of questions**, writes the conditions you stated into the
file with their provenance, and **leaves every energy slot empty** — that is the part only you can
supply.
