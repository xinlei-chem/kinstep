# -*- coding: utf-8 -*-
# KinSTEP - a modality-agnostic microkinetic solver.
# Copyright (C) 2026 The KinSTEP authors.
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.  It is distributed WITHOUT ANY WARRANTY; see the GNU General Public
# License (LICENSE) for details.
r"""[English] Arbitrary-precision backend (mpmath): resolves 1-z below the double-precision
floor. Entered automatically from measured criteria, solving the same equations
with the same analytic Jacobian as the double path.

(The commentary below is in Chinese: it records why each choice was made
and which measurement forced it.)

高精度后端（T11）—— 用 mpmath 把 double 救不回来的那段算出来。

为什么需要它：**两轮测量把它从"锦上添花"变成了"必需"**
--------------------------------------------------------
原计划排在最后（台账第二节第 5 条："先用 double，遇到不收敛再上 mpmath"）。
但 T7 与 T8 在**同一个边界**上各撞了一次：

  · T7（CATKINAS 401 点）：`E_N ≥ +1.84` 的 6 个点上，**所有步的最小 1−Z 中位数
    = 8.88e-16 = 4·eps** —— `1−z` 压在 double 的最后一位，净速率不可得；
  · T8（DRC）：点 300/330（E_N=+1.0/+1.3）X_I 退化到 1.6e-06 / 5.2e-05，
    因为 `ln r` 本身精度退化，**对它求导只会更差**。

★ 精度到底丢在哪一步（这决定了高精度要用在哪）
----------------------------------------------
    ln z  = ΔG/k_BT + Σ ν·ln a        ← **O(10) 量级项之和**
    1 − z = −expm1(ln z)

`expm1` 已经把"从 ln z 到 1−z"做到了全精度（T6 实测：1−z 真值 1e-12/1e-14 时
相对误差均为 0，而 `1 − exp(ln z)` 分别是 2.2e-05 / 8.0e-04）。
所以**瓶颈不在 expm1，而在 ln z 自身的绝对精度**：
double 下它是 ~1e-16，于是 1−z ≈ 1e-30 时有效数字一位不剩。

⇒ 高精度**只需要用在求和这一步**。50 位十进制下 ln z 的绝对精度 ~1e-49，
  足以分辨 1−z ~ 1e-40。本模块因此把整条"ln 空间牛顿"搬到 mpmath 上。

★ 变量选择与 double 版**完全一致**（solve.steady_state 的那套）
    未知数 x = ln θ（表面物种 + 每个位点的空位）
    方程 = n_surf 个稳态 + n_site 个位点守恒
    解析雅可比：∂r_net[s]/∂x_j = R_f[s]·[a_r[s,j]·(1−z_s) − z_s·(a_p−a_r)[s,j]]
  —— 这样 double 与高精度两条路**算的是同一个方程组**，可以逐数字互校（本模块有测试）。

★ mpmath 是**可选依赖**
本包核心只要 numpy + scipy。没装 mpmath 时 `available()` 返回 False，
`steady_state(..., precision="auto")` 会退回 double 并**明确告知哪些步精度不足**，
而不是静默给出坏数。
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from .mech import Mechanism
from .rates import _side_matrices

_MP = None


def available() -> bool:
    """mpmath 在不在。**不在也能跑**，只是没有高精度这一级。"""
    global _MP
    if _MP is None:
        try:
            import mpmath as mp
            _MP = mp
        except ImportError:
            _MP = False
    return _MP is not False


def _mp():
    if not available():
        raise RuntimeError(
            "高精度后端需要 mpmath，但没装。\n"
            "  装： pip install mpmath   （纯 Python，无编译依赖）\n"
            "  或： 用 precision='double'，并读 SteadyResult.precision_table() "
            "了解哪些步的净速率精度不足。")
    return _MP


@dataclass
class HPResult:
    """高精度稳态解。数值同时给 float（好用）与 mpf（不丢精度）两份。"""
    theta: Dict[str, float]
    theta_vac: Dict[str, float]
    rate_fwd: np.ndarray
    zrev: np.ndarray
    one_minus_z: np.ndarray            # ★ 这才是 double 算不出来的那个量
    rate_net: np.ndarray
    ln_z: np.ndarray
    residual: np.ndarray
    n_iter: int
    dps: int
    converged: bool
    theta_mp: Dict[str, object] = field(default_factory=dict)
    #: 流动反应器时的出口气相分压
    p_gas_out: Dict[str, float] = field(default_factory=dict)
    rate_net_mp: List[object] = field(default_factory=list)
    one_minus_z_mp: List[object] = field(default_factory=list)
    notes: List[str] = field(default_factory=list)

    def report(self) -> str:
        L = ["高精度稳态解（%d 位十进制，%d 次牛顿，%s）"
             % (self.dps, self.n_iter, "收敛" if self.converged else "★未收敛")]
        L.append("  覆盖度：" + "  ".join("%s=%.10g" % (k, v) for k, v in self.theta.items()))
        L.append("  空位　：" + "  ".join("%s=%.10g" % (k, v) for k, v in self.theta_vac.items()))
        L.append("  逐步 1−z（double 下多半是 0 或噪声）：")
        for i, v in enumerate(self.one_minus_z_mp):
            L.append("     步%d  1−z = %s" % (i + 1, _MP.nstr(v, 12)))
        L.append("  逐步净速率：")
        for i, v in enumerate(self.rate_net_mp):
            L.append("     步%d  r = %s" % (i + 1, _MP.nstr(v, 12)))
        L.append("  残差 ‖F‖∞ = %.3e" % np.max(np.abs(self.residual)))
        return "\n".join(L)


def steady_state_hp(mech: Mechanism, Ga_fwd, dG_rxn, T_K: float,
                    p_gas: Dict = None, c_liq: Dict = None,
                    A_pre: float = None, prefactor=None, kT_eV: float = None,
                    dps: int = 50, theta_guess: Dict[str, float] = None,
                    max_iter: int = 200, tol_exp: int = 12,
                    reversible=None, flow=None, gas_guess: Dict = None,
                    verbose: bool = False) -> HPResult:
    """高精度稳态求解。

    参数
    ----
    Ga_fwd, dG_rxn : 逐步的正向活化自由能与反应自由能（eV）
    kT_eV : 直接给 k_BT（eV）。**与别的程序对数时用这个**，免得常数口径不同
            （T7 就是因为 CATKINAS 用 CODATA-2006 的 k_B 才要显式指定）。
    A_pre : TST 前因子；不给就用 k_BT/h（以 mpmath 的常数算）
    dps : 十进制位数。50 位可分辨 1−z ~ 1e-40；要更深就加大
    tol_exp : 收敛判据取 10^(−dps+tol_exp)

    返回 `HPResult`。`one_minus_z_mp` 是本模块存在的理由。
    """
    mp = _mp()
    old_dps = mp.mp.dps
    mp.mp.dps = int(dps)
    try:
        return _solve_hp(mp, mech, Ga_fwd, dG_rxn, T_K, p_gas, c_liq, A_pre,
                         prefactor, kT_eV, dps, theta_guess, max_iter, tol_exp,
                         reversible, flow, gas_guess, verbose)
    finally:
        mp.mp.dps = old_dps


def _solve_hp(mp, mech, Ga_fwd, dG_rxn, T_K, p_gas, c_liq, A_pre, prefactor,
              kT_eV, dps, theta_guess, max_iter, tol_exp, reversible, flow,
              gas_guess, verbose):
    surf = mech.surface_species
    vacs = [mech.vacancy_of(s) for s in mech.sites if mech.vacancy_of(s) in mech.species]
    # ★ T9b 收尾加的：流动反应器时气相分压也是未知数。
    #   与 double 版**完全同构**（同样的未知数排列、同样的方程、同样的解析雅可比），
    #   所以两条路仍可逐数字互校。
    gas = list(mech.gas_species) if flow is not None else []
    names = surf + vacs + gas
    n_s, n_v, n_g = len(surf), len(vacs), len(gas)
    n_x = n_s + n_v + n_g
    if n_x == 0:
        raise RuntimeError("机理里没有表面物种")
    site_of = {mech.vacancy_of(s): s for s in mech.sites}
    idx = {n: i for i, n in enumerate(mech.species_order)}
    a_r_np, a_p_np = _side_matrices(mech)
    nu_np = np.asarray(mech.nu, float)
    n_step = mech.n_step

    # ── 常数（高精度）──
    if kT_eV is not None:
        kT = mp.mpf(repr(float(kT_eV)))
    else:
        # CODATA-2018：k_B = R/N_A，e = 1.602176634e-19
        K_B = mp.mpf("8.31446261815324") / mp.mpf("6.02214076e23")
        kT = K_B * mp.mpf(repr(float(T_K))) / mp.mpf("1.602176634e-19")
    if A_pre is not None:
        A = mp.mpf(repr(float(A_pre)))
    else:
        K_B = mp.mpf("8.31446261815324") / mp.mpf("6.02214076e23")
        A = K_B * mp.mpf(repr(float(T_K))) / mp.mpf("6.62607015e-34")

    Ga = [mp.mpf(repr(float(v))) for v in np.asarray(Ga_fwd, float)]
    dG = [mp.mpf(repr(float(v))) for v in np.asarray(dG_rxn, float)]
    pref = ([mp.mpf(repr(float(v))) for v in np.asarray(prefactor, float)]
            if prefactor is not None else [mp.mpf(1) for _ in range(n_step)])
    rev = (np.asarray(reversible, bool) if reversible is not None
           else np.array([st.reversible for st in mech.steps], bool))

    ln_kf = [mp.log(A) - Ga[s] / kT for s in range(n_step)]

    # ── 固定活度（气/液）的 ln ──
    ln_a_fixed = [None] * len(mech.species_order)
    p_gas = p_gas or {}
    c_liq = c_liq or {}
    for nm, sp in mech.species.items():
        j = idx[nm]
        if sp.phase == "gas":
            v = p_gas.get(nm, p_gas.get(sp.formula, 0.0))
        elif sp.phase == "liq":
            v = c_liq.get(nm, c_liq.get(sp.formula, 0.0))
        elif sp.phase == "electron":
            v = 1.0
        else:
            continue                                  # 表面/空位是未知数
        if nm in gas:
            continue                                  # 流动时气相也是未知数
        ln_a_fixed[j] = mp.log(mp.mpf(repr(float(v)))) if v > 0 else mp.mpf("-inf")

    # 预存每步非零指数的 (物种列, 指数)，避免在内层乘零
    rows_r, rows_p = [], []
    for s in range(n_step):
        rows_r.append([(j, mp.mpf(repr(float(a_r_np[s, j]))))
                       for j in range(len(mech.species_order)) if a_r_np[s, j] != 0.0])
        rows_p.append([(j, mp.mpf(repr(float(a_p_np[s, j]))))
                       for j in range(len(mech.species_order)) if a_p_np[s, j] != 0.0])

    col_of = [idx[nm] for nm in names]

    def residual(x):
        """x = [ln θ_surf..., ln θ_vac...]（mpf 列表）→ (F, J, aux)。"""
        ln_a = list(ln_a_fixed)
        for k, nm in enumerate(names):
            ln_a[idx[nm]] = x[k]
        th = [mp.e ** xi for xi in x]
        R_f, z, omz, r_net, ln_z = [], [], [], [], []
        for s in range(n_step):
            lf = ln_kf[s]
            bad = False
            for j, a in rows_r[s]:
                if ln_a[j] == mp.mpf("-inf"):
                    bad = True; break
                lf += a * ln_a[j]
            if bad:
                R_f.append(mp.mpf(0)); z.append(mp.mpf(0)); omz.append(mp.mpf(1))
                r_net.append(mp.mpf(0)); ln_z.append(mp.mpf("-inf")); continue
            lz = dG[s] / kT
            zbad = False
            for j, a in rows_p[s]:
                if ln_a[j] == mp.mpf("-inf"):
                    zbad = True; break
                lz += a * ln_a[j]
            for j, a in rows_r[s]:
                lz -= a * ln_a[j]
            if (not rev[s]) or zbad:
                zz, om = mp.mpf(0), mp.mpf(1)
                lz = mp.mpf("-inf")
            else:
                zz = mp.e ** lz
                om = -mp.expm1(lz)            # ★ 关键：高精度的 expm1
            rf = pref[s] * mp.e ** lf
            R_f.append(rf); z.append(zz); omz.append(om)
            r_net.append(rf * om); ln_z.append(lz)

        if n_g:
            kfl = flow.k_flow(T_K)
            kfl = mp.mpf(repr(float(kfl)))
            p_in_v = [mp.mpf(repr(float(flow.p_in.get(nm, flow.p_in.get(
                mech.species[nm].formula, 0.0))))) for nm in gas]

        F = [mp.mpf(0)] * n_x
        for i in range(n_s):
            acc = mp.mpf(0)
            ci = col_of[i]
            for s in range(n_step):
                v = nu_np[s, ci]
                if v != 0.0:
                    acc += mp.mpf(repr(float(v))) * r_net[s]
            F[i] = acc
        for v_i, vn in enumerate(vacs):
            site = site_of[vn]
            acc = th[n_s + v_i]
            for k, nm in enumerate(surf):
                if mech.species[nm].site == site:
                    acc += th[k]
            F[n_s + v_i] = acc - 1

        J = mp.matrix(n_x, n_x)
        for i in range(n_s):
            ci = col_of[i]
            for k in range(n_x):
                ck = col_of[k]
                acc = mp.mpf(0)
                for s in range(n_step):
                    nv = nu_np[s, ci]
                    if nv == 0.0:
                        continue
                    ar = mp.mpf(repr(float(a_r_np[s, ck])))
                    ap = mp.mpf(repr(float(a_p_np[s, ck])))
                    acc += mp.mpf(repr(float(nv))) * R_f[s] * (ar * omz[s] - z[s] * (ap - ar))
                J[i, k] = acc
        for v_i, vn in enumerate(vacs):
            site = site_of[vn]
            J[n_s + v_i, n_s + v_i] = th[n_s + v_i]
            for k, nm in enumerate(surf):
                if mech.species[nm].site == site:
                    J[n_s + v_i, k] = th[k]
        for g_i, gn in enumerate(gas):
            cg = idx[gn]
            for k in range(n_x):
                ck = col_of[k]
                acc = mp.mpf(0)
                for s in range(n_step):
                    nv = nu_np[s, cg]
                    if nv == 0.0:
                        continue
                    ar = mp.mpf(repr(float(a_r_np[s, ck])))
                    ap = mp.mpf(repr(float(a_p_np[s, ck])))
                    acc += mp.mpf(repr(float(nv))) * R_f[s] * (ar * omz[s] - z[s] * (ap - ar))
                J[n_s + n_v + g_i, k] = acc
            J[n_s + n_v + g_i, n_s + n_v + g_i] -= kfl * th[n_s + n_v + g_i]
        return F, J, dict(R_f=R_f, z=z, omz=omz, r_net=r_net, ln_z=ln_z, th=th)

    # ── 初值 ──
    g = dict(theta_guess or {})
    x = []
    for nm in surf:
        x.append(mp.log(mp.mpf(repr(max(float(g.get(nm, 1e-8)), 1e-300)))))
    for vn in vacs:
        site = site_of[vn]
        used = sum(float(g.get(nm, 0.0)) for nm in surf
                   if mech.species[nm].site == site)
        x.append(mp.log(mp.mpf(repr(max(1.0 - used, 1e-300)))))
    gg = dict(gas_guess or {})
    for gn in gas:
        v = gg.get(gn, gg.get(mech.species[gn].formula,
                              flow.p_in.get(gn, flow.p_in.get(
                                  mech.species[gn].formula, 1e-12))))
        x.append(mp.log(mp.mpf(repr(max(float(v), 1e-300)))))

    # 行缩放（固定，用初值处的 Σ|ν|·R_f）
    _, _, aux0 = residual(x)
    rs = []
    for i in range(n_s):
        ci = col_of[i]
        acc = mp.mpf(0)
        for s in range(n_step):
            acc += abs(mp.mpf(repr(float(nu_np[s, ci])))) * aux0["R_f"][s]
        rs.append(acc if acc > 0 else mp.mpf(1))
    rs += [mp.mpf(1)] * n_v
    if n_g:
        kfl0 = mp.mpf(repr(float(flow.k_flow(T_K))))
        for g_i, gn in enumerate(gas):
            cg = idx[gn]
            acc = mp.mpf(0)
            for s in range(n_step):
                acc += abs(mp.mpf(repr(float(nu_np[s, cg])))) * aux0["R_f"][s]
            pin = mp.mpf(repr(float(flow.p_in.get(gn, flow.p_in.get(
                mech.species[gn].formula, 0.0)))))
            sc = acc + kfl0 * max(pin, aux0["th"][n_s + n_v + g_i])
            rs.append(sc if sc > 0 else mp.mpf(1))

    tol = mp.mpf(10) ** (-(int(dps) - int(tol_exp)))
    step_max = mp.mpf(4)
    F, J, aux = residual(x)
    nF = max(abs(F[i] / rs[i]) for i in range(n_x))
    it = 0
    for it in range(1, max_iter + 1):
        if nF <= tol:
            break
        Fs = mp.matrix([[-F[i] / rs[i]] for i in range(n_x)])
        Js = mp.matrix(n_x, n_x)
        for i in range(n_x):
            for k in range(n_x):
                Js[i, k] = J[i, k] / rs[i]
        try:
            dx = mp.lu_solve(Js, Fs)
        except Exception:
            try:
                dx = mp.qr_solve(Js, Fs)[0]
            except Exception:
                break
        mx = max(abs(dx[i]) for i in range(n_x))
        if mx == 0:
            break
        scl = min(mp.mpf(1), step_max / mx)
        lam, ok_ls = mp.mpf(1), False
        for _ in range(60):
            # 上界：θ ≤ 1 ⇒ ln θ ≤ 0；但**气相分压可以 >1**（算例到 750 bar）
            xn = [min(max(x[i] + lam * scl * dx[i], mp.mpf(-1500)),
                      mp.mpf(50) if i >= n_s + n_v else mp.mpf(0))
                  for i in range(n_x)]
            Fn, Jn, auxn = residual(xn)
            nn = max(abs(Fn[i] / rs[i]) for i in range(n_x))
            if nn < nF:
                x, F, J, aux, nF, ok_ls = xn, Fn, Jn, auxn, nn, True
                break
            lam /= 2
        if verbose:
            print("   hp 牛顿第 %d 步：‖F‖ = %s" % (it, mp.nstr(nF, 6)))
        if not ok_ls:
            break

    th = aux["th"]
    theta = {nm: float(th[i]) for i, nm in enumerate(surf)}
    tvac = {vn: float(th[n_s + i]) for i, vn in enumerate(vacs)}
    f2 = lambda seq: np.array([float(v) for v in seq])
    return HPResult(
        theta=theta, theta_vac=tvac,
        rate_fwd=f2(aux["R_f"]), zrev=f2(aux["z"]), one_minus_z=f2(aux["omz"]),
        rate_net=f2(aux["r_net"]),
        ln_z=np.array([float(v) if v != mp.mpf("-inf") else -np.inf
                       for v in aux["ln_z"]]),
        residual=np.array([float(F[i] / rs[i]) for i in range(n_x)]),
        n_iter=it, dps=int(dps), converged=bool(nF <= tol),
        theta_mp={nm: th[i] for i, nm in enumerate(surf)},
        p_gas_out={nm: float(th[n_s + n_v + i]) for i, nm in enumerate(gas)},
        rate_net_mp=list(aux["r_net"]), one_minus_z_mp=list(aux["omz"]),
        notes=["未知数与方程组与 double 版**完全相同**（ln θ + 位点守恒 + 解析雅可比），"
               "所以两条路可以逐数字互校",
               "收敛判据 ‖F‖∞ ≤ 1e-%d（缩放残差）" % (int(dps) - int(tol_exp))])


def needs_high_precision(res, omz_floor: float = 1e-10,
                         flux_floor: float = 1e-6) -> Tuple[bool, List[str]]:
    """判断一个 double 解是否该升到高精度。返回 (要不要, 原因列表)。

    两条触发条件，都是**测出来的**而不是猜的：
      ① 某个可逆步的 `1−z` 下溢为 0，或小于 `omz_floor`
         ⇒ 该步净速率的相对精度 ≈ eps/(1−z) 已不可用；
      ② `flux_fix > flux_floor` —— 通量投影改动过大，说明原始解的通量严重不自洽。
    """
    why = []
    for i, z in enumerate(np.asarray(res.zrev, float)):
        omz = 1.0 - z
        if omz <= 0.0:
            why.append("第%d步 1−z 下溢为 0" % (i + 1))
        elif omz < omz_floor:
            why.append("第%d步 1−z = %.2e < %.0e（净速率相对精度仅 %.1e）"
                       % (i + 1, omz, omz_floor, 1e-16 / omz))
    fx = float(getattr(res, "flux_fix", 0.0))
    if fx > flux_floor:
        why.append("通量投影改动 %.2e > %.0e（原始解通量严重不自洽）" % (fx, flux_floor))
    # ★ T9b 收尾加的第三条：**气相物料平衡的相消**也会把残差卡住。
    #   实测高温低压角：表面残差 1e-15、位点守恒 0，而**气相 2e-11 ~ 3e-11** ——
    #   `prod + k_flow·(p_in − p) = 0` 两项相消，double 下就是降不到 1e-11。
    #   所以那批点「解其实是对的，只是判据判它没收敛」。
    #
    #   ⚠ 诊断时还被另一个盲点坑过：`SteadyResult.residual` **只存表面那一段**，
    #     打印出来是 1e-15 却判未收敛 —— 真正大的气相项根本没被存下来。
    #     现在 `residual_full` 把完整残差留住了，这条判据才写得出来。
    ff = getattr(res, "residual_full", None)
    if ff is not None and len(ff) > 0:
        fmax = float(np.max(np.abs(ff)))
        if fmax > 1e-11:
            why.append("完整残差 %.2e > 1e-11（多半是气相物料平衡的相消）" % fmax)
    return bool(why), why
