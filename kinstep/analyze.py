# -*- coding: utf-8 -*-
# KinSTEP - a modality-agnostic microkinetic solver.
# Copyright (C) 2026 The KinSTEP authors.
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.  It is distributed WITHOUT ANY WARRANTY; see the GNU General Public
# License (LICENSE) for details.
r"""[English] Analysis layer: degree of rate control, reaction orders, apparent activation
energy and Tafel slope, by central differences with Richardson extrapolation and
a step size chosen from the Campbell sum rule.

(The commentary below is in Chinese: it records why each choice was made
and which measurement forced it.)

分析层（T8）：速率控制度 DRC、反应级数、表观活化能、Tafel 斜率。

★ 全包**一套统一定义**（与 CATKINAS 对齐，口径反推自它的 401 点输出）
--------------------------------------------------------------------
    X_x  =  ∂ln r / ∂(−G_x / k_BT)  =  −k_BT · ∂ln r / ∂G_x

对不同的 x 就得到不同的 DRC：

| x | 名字 | 怎么扰动（这是定义的全部内容）|
|---|---|---|
| 过渡态 s | **台阶 DRC** `X_B` | `Ga_fwd[s] += δ`，**`dG_rxn` 不动** ⇒ 只抬过渡态 |
| 中间体 k | **中间体 DRC** `X_I` | 抬 G_k：`dG_rxn[s] += δ·ν[s,k]`，`Ga_fwd[s] −= δ·a_r[s,k]`（过渡态不动）|
| 气相 g | **−反应级数** `X_PC` | 等价于扰动 `ln p_g`（见下方证明）|
| 描述符 | `X_G` | 沿扫描曲线差分（见 `drc_descriptor`）|

★ 为什么"扰动 ln p"与"抬气相 μ°"等价（省掉一套代码）
    ln z   = ΔG°/k_BT + Σν·ln a        ⇒ 抬 ln p_g 得 ν[s,g]·ε；抬 μ°_g 得 ν[s,g]·δ/k_BT
    ln r_f = ln kf   + Σa_r·ln a       ⇒ 抬 ln p_g 得 a_r[s,g]·ε；
                                          抬 μ°_g（过渡态不动）⇒ Ga 降 a_r[s,g]·δ ⇒ 同样
    ⇒ 取 ε = δ/k_BT 两者逐项相同。所以 `X_PC = −∂ln r/∂ln p_g = −(反应级数)`。

★★ 免费自检：**Campbell 求和规则 Σ_s X_B = 1**
---------------------------------------------
这条规则不依赖机理、不依赖能量，所以它既是**正确性检查**，
也是**选差分步长 δ 的标尺** —— 本模块的 `calibrate_delta()` 就是用它来定 δ 的，
不靠经验值（见该函数的说明）。

★ 数值上的要点
--------------
中心差分的总误差 ≈ 截断 O(δ²) + 舍入 O(ε_r/δ)，其中 ε_r 是 ln r 的可达相对精度。
T7 实测 ε_r ≈ 1e-10（好点）⇒ 最优 δ ≈ (ε_r)^(1/3) ≈ 5e-4 eV 量级。
但在**准平衡区** ε_r 会退化到 1e-2，那里 DRC 也跟着退化 —— `drc_steps` 会把
每个点的求和规则残差一起返回，**让精度可见，不假装算准了**。
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from .driver import DriverStack
from .mech import Mechanism
from .provenance import Val
from .rates import Energetics, StepEnergy, compile_rates, apply_state
from .solve import SteadyResult, steady_state
from .units import K_B_EV


class AnalyzeError(RuntimeError):
    """DRC / 分析层的输入或前提有问题。"""


# ══════════════════════ 从数组造 Energetics（扰动要用）══════════════════════
def energetics_from_arrays(mech: Mechanism, Ga_fwd, dG_rxn, A_pre=None,
                           source="user_given") -> Energetics:
    """把 (Ga_fwd, dG_rxn) 两个数组包成 `Energetics`。

    DRC 要反复**微扰能量再重解**，走 `Energetics(mech, {1: dict(...)})` 那条
    面向用户的路会不断重建 Val 并过闸；这里直接给数值层用。
    """
    Ga = np.asarray(Ga_fwd, float)
    dG = np.asarray(dG_rxn, float)
    tab = {}
    for s in range(mech.n_step):
        spec = dict(mode="tst", Ga_fwd=Val(float(Ga[s]), source),
                    dG_rxn=Val(float(dG[s]), source))
        if A_pre is not None:
            spec["A_pre"] = Val(float(A_pre), "user_given")
        tab[s + 1] = StepEnergy(**spec)
    return Energetics(mech, tab)


def _base_arrays(en: Energetics) -> Tuple[np.ndarray, np.ndarray, Optional[float]]:
    """取出 Energetics 的 (Ga_fwd, dG_rxn, A_pre)。只支持全 tst 的输入。"""
    cr = compile_rates(en, strict_energy=False)
    if np.any(cr.mode != 0):
        raise AnalyzeError(
            "DRC 目前只支持全部步都用 mode='tst' 的输入。\n"
            "碰撞理论/直接给 k⁰ 的步没有"
            "可微扰的 Ga_fwd —— 请先把它们折算成等效的 Ga_fwd（CATKINAS 的 Ean 就是这么做的）。")
    A = None
    if np.any(np.isfinite(cr.A_pre)):
        a = cr.A_pre[np.isfinite(cr.A_pre)]
        if np.ptp(a) > 0:
            raise AnalyzeError("各步的 A_pre 不同，DRC 暂不支持（请统一前因子）")
        A = float(a[0])
    return cr.Ga_in.copy(), cr.dG_in.copy(), A


# ══════════════════════════════ 结果容器 ══════════════════════════════
@dataclass
class DRCResult:
    X_B: np.ndarray                 # (n_step, n_target) 台阶 DRC
    X_I: np.ndarray                 # (n_surf+n_vac, n_target) 中间体 DRC
    X_PC: np.ndarray                # (n_gas, n_target) = −反应级数
    order: np.ndarray               # (n_gas, n_target) 反应级数本身
    Ea_app: np.ndarray              # (n_target,) 表观活化能 eV
    sum_B: np.ndarray               # (n_target,) Σ_s X_B，应为 1
    sum_I: np.ndarray               # (n_target,) Σ_k X_I
    delta: float
    step_names: List[str]
    inter_names: List[str]
    gas_names: List[str]
    notes: List[str] = field(default_factory=list)

    def sum_rule_ok(self, tol=1e-3):
        return bool(np.all(np.abs(self.sum_B - 1.0) <= tol))

    def report(self, target=0) -> str:
        from .i18n import t
        tgt = self.step_names[target] if target < len(self.step_names) else target
        L = [t("DRC（δ = %.2e eV，目标 = %s）" % (self.delta, tgt),
               "DRC (delta = %.2e eV, target = %s)" % (self.delta, tgt))]
        L.append("  " + t("台阶 DRC（Campbell）：",
                          "step-wise DRC (Campbell):"))
        for i, nm in enumerate(self.step_names):
            L.append("     " + t("步%-2d" % (i + 1), "step %-2d" % (i + 1))
                     + " %-28s X_B = %+.6f" % (nm[:28], self.X_B[i, target]))
        L.append("     ★ " + t("求和规则 Σ X_B = %.9f（应为 1，偏差 %.1e）"
                               % (self.sum_B[target], abs(self.sum_B[target] - 1.0)),
                               "sum rule  sum X_B = %.9f (should be 1; deviation %.1e)"
                               % (self.sum_B[target], abs(self.sum_B[target] - 1.0))))
        L.append("  " + t("中间体 DRC：", "intermediate DRC:"))
        for i, nm in enumerate(self.inter_names):
            L.append("     %-8s X_I = %+.6f" % (nm, self.X_I[i, target]))
        L.append("     " + t("Σ X_I = %+.6f", "sum X_I = %+.6f") % self.sum_I[target])
        L.append("  " + t("反应级数：", "reaction orders: ")
                 + "  ".join("%s: %+.4f" % (g, self.order[i, target])
                             for i, g in enumerate(self.gas_names)))
        L.append("  " + t("表观活化能 Ea_app = %.6f eV",
                          "apparent activation energy Ea_app = %.6f eV")
                 % self.Ea_app[target])
        return "\n".join(L)


# ══════════════════════════════ 主计算 ══════════════════════════════
def _solve(mech, en, state, driver, p_gas, c_liq, guess, tol,
           precision="double", dps=50, solve_kw=None):
    kw = dict(state=state, driver=driver, p_gas=p_gas, c_liq=c_liq, tol=tol,
              precision=precision, dps=dps, strict_energy=False)
    # ★ 扰动解若 Newton 失手就会掉进完整 ODE 阶梯；不给出口的话，描述符两端的
    #   单点 DRC 能跑掉数十秒（20 次扰动解 × 每次一次完整 ODE）。出图/扫描这类
    #   不需要极限精度的场合，用 solve_kw 把 ode_stages/ode_budget_s 压住即可。
    if solve_kw:
        kw.update(solve_kw)
    r = steady_state(mech, en, theta_guess=guess,
                     ladder=("newton",) if guess else ("ode", "newton"), **kw)
    if not r.converged:
        r2 = steady_state(mech, en, ladder=("ode", "newton"), **kw)
        if r2.converged:
            return r2
    return r


def drc(mech: Mechanism, en: Energetics, state: Dict, driver: DriverStack = None,
        p_gas: Dict = None, c_liq: Dict = None, delta: float = 1e-3,
        base: SteadyResult = None, tol: float = 1e-12,
        dT: float = None, richardson: bool = False,
        precision: str = "double", dps: int = 50,
        solve_kw: Dict = None) -> DRCResult:
    """算全套 DRC：台阶 / 中间体 / 气相（反应级数）/ 表观活化能。

    每个参数用**中心差分**（两次重解），初值一律取未扰动解 ⇒ 快且稳。
    目标量是**各步的净速率**（线性机理下它们成比例，DRC 相同；并行路径下不同）。

    `richardson=True`：在 δ 与 δ/2 各算一次再做 **Richardson 外推**
    `X ≈ (4·X(δ/2) − X(δ))/3`，把二阶截断误差消掉（代价是 2 倍重解次数）。

    ★ 为什么确定能这么做（2026-10-22 实测，不是套公式）
    ----------------------------------------------------
    拿 CATKINAS 的 `DRC_I` 当真值，扫 δ 从 4e-3 到 1.25e-4（32 倍范围）：

        δ=4.0e-3 → 差 1.524e-04 ;  差/δ² = 9.523
        δ=2.0e-3 → 差 3.810e-05 ;  差/δ² = 9.525
        δ=1.0e-3 → 差 9.526e-06 ;  差/δ² = 9.526
        δ=5.0e-4 → 差 2.381e-06 ;  差/δ² = 9.526
        δ=2.5e-4 → 差 5.954e-07 ;  差/δ² = 9.526
        δ=1.25e-4→ 差 1.488e-07 ;  差/δ² = 9.526

    `差/δ²` 恒定到 4 位有效数字 ⇒ 残差**纯粹是中心差分的二阶截断**，
    不是定义差异（若是定义差异，差不会随 δ 变）。
    ⇒ 外推后实测：X_I(H*) = −0.438578034465 vs CATKINAS −0.438578034，差 ~5e-10。
    """
    T_K = float(state["T_K"])
    kT = K_B_EV * T_K
    Ga0, dG0, A_pre = _base_arrays(en)
    surf = mech.surface_species
    vacs = [mech.vacancy_of(s) for s in mech.sites if mech.vacancy_of(s) in mech.species]
    inter = surf + vacs
    gases = mech.gas_species
    idx = {n: i for i, n in enumerate(mech.species_order)}
    nu = np.asarray(mech.nu, float)
    from .rates import _side_matrices
    a_r, _ = _side_matrices(mech)

    if base is None:
        base = _solve(mech, en, state, driver, p_gas, c_liq, None, tol, precision, dps,
                      solve_kw)
    if not base.converged:
        raise AnalyzeError("未扰动的稳态就没收敛，DRC 无从谈起")
    guess = dict(base.theta)
    n_t = mech.n_step

    def lnr(Ga, dG, st, pg):
        r = _solve(mech, energetics_from_arrays(mech, Ga, dG, A_pre),
                   st, driver, pg, c_liq, guess, tol, precision, dps, solve_kw)
        if not r.converged:
            return np.full(n_t, np.nan)
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.log(np.abs(r.rate_net))

    notes = []
    # ── ① 台阶 DRC：只抬过渡态（dG 不动）──
    X_B = np.zeros((mech.n_step, n_t))
    for s in range(mech.n_step):
        gp, gm = Ga0.copy(), Ga0.copy()
        gp[s] += delta; gm[s] -= delta
        X_B[s] = -kT * (lnr(gp, dG0, state, p_gas) - lnr(gm, dG0, state, p_gas)) / (2 * delta)

    # ── ② 中间体 DRC：抬 G_k，过渡态不动 ──
    #   dG_rxn[s] += δ·ν[s,k]   （产物侧为正）
    #   Ga_fwd[s] −= δ·a_r[s,k] （从初态量的势垒变矮）
    X_I = np.zeros((len(inter), n_t))
    for j, nm in enumerate(inter):
        c = idx[nm]
        dnu, dar = nu[:, c], a_r[:, c]
        gp = Ga0 - delta * dar; dp = dG0 + delta * dnu
        gm = Ga0 + delta * dar; dm = dG0 - delta * dnu
        X_I[j] = -kT * (lnr(gp, dp, state, p_gas) - lnr(gm, dm, state, p_gas)) / (2 * delta)

    # ── ③ 气相：X_PC = −∂ln r/∂ln p（= −反应级数）──
    eps_p = delta / kT          # 与能量扰动等价的 ln p 步长（见模块 docstring）
    order = np.zeros((len(gases), n_t))
    pg0 = dict(p_gas or {})
    for j, g in enumerate(gases):
        key = g if g in pg0 else mech.species[g].formula
        p0 = float(pg0.get(key, 0.0))
        if p0 <= 0:
            notes.append("气相 %s 的分压为 0，反应级数置 0" % g)
            continue
        pp = dict(pg0); pp[key] = p0 * np.exp(eps_p)
        pm = dict(pg0); pm[key] = p0 * np.exp(-eps_p)
        order[j] = (lnr(Ga0, dG0, state, pp) - lnr(Ga0, dG0, state, pm)) / (2 * eps_p)
    X_PC = -order

    # ── ④ 表观活化能：Ea_app = k_BT²·∂ln r/∂T ──
    if dT is None:
        dT = max(T_K * 1e-4, 1e-3)
    sp = dict(state); sp["T_K"] = T_K + dT
    sm = dict(state); sm["T_K"] = T_K - dT
    Ea_app = K_B_EV * T_K ** 2 * (lnr(Ga0, dG0, sp, p_gas)
                                  - lnr(Ga0, dG0, sm, p_gas)) / (2 * dT)

    res = DRCResult(X_B=X_B, X_I=X_I, X_PC=X_PC, order=order, Ea_app=Ea_app,
                    sum_B=X_B.sum(axis=0), sum_I=X_I.sum(axis=0), delta=delta,
                    step_names=[st.raw for st in mech.steps],
                    inter_names=inter, gas_names=gases, notes=notes)
    if richardson:
        half = drc(mech, en, state, driver, p_gas, c_liq, delta=delta / 2.0,
                   base=base, tol=tol, dT=dT, richardson=False,
                   precision=precision, dps=dps)
        rich = lambda a, b: (4.0 * b - a) / 3.0        # X ≈ (4X(δ/2) − X(δ))/3
        res = DRCResult(
            X_B=rich(res.X_B, half.X_B), X_I=rich(res.X_I, half.X_I),
            X_PC=rich(res.X_PC, half.X_PC), order=rich(res.order, half.order),
            Ea_app=rich(res.Ea_app, half.Ea_app),
            sum_B=rich(res.X_B, half.X_B).sum(axis=0),
            sum_I=rich(res.X_I, half.X_I).sum(axis=0),
            delta=delta, step_names=res.step_names, inter_names=inter,
            gas_names=gases,
            notes=notes + ["Richardson 外推（δ=%.2e 与 %.2e）⇒ 消去二阶截断误差"
                           % (delta, delta / 2.0)])
    return res


def calibrate_delta(mech, en, state, driver=None, p_gas=None, c_liq=None,
                    candidates=(5e-3, 2e-3, 1e-3, 5e-4, 2e-4, 1e-4, 5e-5),
                    tol=1e-12, verbose=False):
    """★ 用 **Campbell 求和规则 Σ X_B = 1** 来选差分步长 δ。

    为什么不用经验值
    ----------------
    中心差分的误差 = 截断 O(δ²) + 舍入 O(ε_r/δ)，而 ε_r（ln r 的可达精度）
    **随体系和工况变化好几个数量级**（T7 实测：好点 1e-10，准平衡区退化到 1e-2）。
    固定一个 δ 必然在某些体系上失手。

    而求和规则 Σ_s X_B = 1 **不依赖机理、不依赖能量、不依赖是否准平衡** ——
    它是一个免费的、永远成立的标尺。拿它挑 δ，等于让数据自己说话。

    返回 `(最佳 δ, [(δ, |ΣX_B − 1|), ...])`。
    """
    base = _solve(mech, en, state, driver, p_gas, c_liq, None, tol)
    rows = []  # (δ, |ΣX_B − 1|)
    for dl in candidates:
        try:
            r = drc(mech, en, state, driver, p_gas, c_liq, delta=dl, base=base, tol=tol)
            err = float(np.nanmax(np.abs(r.sum_B - 1.0)))
        except Exception as e:
            err = np.inf
        rows.append((dl, err))
        if verbose:
            print("   δ = %.1e ⇒ |ΣX_B − 1| = %.3e" % (dl, err))
    best = min(rows, key=lambda t: t[1])[0]
    return best, rows


def drc_descriptor(desc_values, rate_curve, T_K):
    """描述符 DRC：`X_G = −k_BT·d ln r / d(描述符)`，沿扫描曲线中心差分。

    与 CATKINAS 的 `DRC_G` 同一口径（已核对：中心差分 1.490668 vs 它的 1.490678）。
    端点用单侧差分。
    """
    x = np.asarray(desc_values, float)
    r = np.asarray(rate_curve, float)
    with np.errstate(divide="ignore", invalid="ignore"):
        lr = np.log(np.abs(r))
    g = np.gradient(lr, x, edge_order=2)
    return -K_B_EV * T_K * g


def tafel_slope(U_values, j_values):
    """Tafel 斜率 ∂U/∂log₁₀|j|（V/decade）。返回逐点值（中心差分）。"""
    U = np.asarray(U_values, float)
    j = np.abs(np.asarray(j_values, float))
    with np.errstate(divide="ignore", invalid="ignore"):
        lg = np.log10(np.where(j > 0, j, np.nan))
    return np.gradient(U, lg, edge_order=2)
