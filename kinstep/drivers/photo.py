# -*- coding: utf-8 -*-
# KinSTEP - a modality-agnostic microkinetic solver.
# Copyright (C) 2026 The KinSTEP authors.
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.  It is distributed WITHOUT ANY WARRANTY; see the GNU General Public
# License (LICENSE) for details.
r"""[English] Photocatalysis driver: quasi-Fermi-level splitting solved self-consistently from
the carrier balance.

(The commentary below is in Chinese: it records why each choice was made
and which measurement forced it.)

光催化 driver（T12）：准费米能级分裂 + 载流子收支自洽。

理论依据（2026-10-22 实地检索核对过，不是凭推导）
--------------------------------------------------
① **准费米能级分裂 ΔE_qfl = E_Fn − E_Fp 就是表面反应的热力学驱动力。**
   这一条在多个独立来源上一致：
   · Beh/…, *JACS* 2023, 10.1021/jacs.3c02333 —— "the splitting between the quasi-Fermi
     levels of the minority and majority carriers provides the thermodynamic driving
     force for catalytic reactions on the metal surface"
   · *Nat. Commun.* 2025, 10.1038/s41467-025-58837-9（用外球氧化还原对测光阳极表面的空穴准费米能级）
   · *ACS Catal.* perspective 10.1021/acscatal.5c07014 —— 驱动力 = 载流子准费米能级与
     吸附质氧化还原电位之差；并指出 "stimulus-dependent microkinetic models … are
     beginning to offer insights"（领域自认这块刚起步）

   ⇒ 能量修正与电催化的 `−ν_e·e·U` **数学同构**，所以 photo 与 electro 可直接叠加成 PEC。

② **光强幂律：低光强一级 → 高光强半级。** 经典出处 Turchi & Ollis 1990；
   K₄Nb₆O₁₇ 分解水实测两个区间都观察到（10.1007/BF00806072）。
   ⚠ **但 0.5 不是普适值**：文献给的是 **0.5–1 之间**，有体系报 ~0.8，
   极高光强还可能转零级；*J. Phys. Chem. C* 2018 (10.1021/acs.jpcc.8b09421) 指出
   半级**只在"界面转移速率远小于复合速率"时**才成立。
   ⇒ 本模块**不把 0.5 写死**，而是让它从载流子收支自洽地解出来 ——
     指数是**算出来的结果**，不是输入的假设。这正是 MKM 相对经验幂律的价值。

★★ 我差点写进去的一个洞（检索后修掉）
--------------------------------------
初版骨架里写死了"电子与空穴各占分裂的一半"。**那只对本征/高注入成立**：
掺杂半导体里多数载流子浓度几乎不随光照变，分裂**几乎全部落在少子一侧**，
幂律也随之不同（n 型里少子驱动的反应全程一级，没有半级区间）。
⇒ 分配比做成参数 `qfl_frac_n`（默认 0.5 = 本征/高注入），掺杂情形自己给。

模型（都取对暗态的相对量，所以暗态下一切归零）
----------------------------------------------
    x ≡ ΔE_qfl / k_BT
    n/n_dark = exp(f·x) ,   p/p_dark = exp((1−f)·x) ,   f = qfl_frac_n
    ⇒ np/(n_dark·p_dark) = exp(x)                      （与 n·p = n_i²·e^x 一致）

    产生（每位点每秒）  G = absorptance · photon_flux / gamma_site_photo
    复合  R_rec(x) = k_bi·(e^x − 1)  +  k_srh·(e^x − 1)/(e^{f x} + e^{(1−f)x})
          —— 两项**在暗态都精确为零**（这是正确的暗极限，不是近似）；
             高注入下第二项 → e^{0.5x}（即 ∝ 载流子浓度），正是 SRH 的高注入行为。
    界面消耗  R_surf(x, θ) = Σ_s (ν_cb[s] + ν_vb[s])·r_net[s]

    稳态：**G = R_rec(x) + R_surf(x, θ)**  → 解 x

两种模式
--------
`given_qfl`     直接给 ΔE_qfl（最简；与 electro 完全同构，可用来互校）
`steady_carrier` 由上面的收支方程**自洽解出** ΔE_qfl —— 这时 `needs_theta=True`，
                 内核会在牛顿内层反复调用（T2.5 留的那个口子正是为此）。
"""
import numpy as np

from ..driver import Driver, register_driver
from ..units import K_B_EV


@register_driver
class PhotoDriver(Driver):
    driver_name = "photo"
    state_params = (
        "photon_flux",         # 入射光子流 (m^-2 s^-1)
        "photon_energy_eV",    # 光子能量（单色）或等效平均值
        "absorptance",         # 吸收率 0~1
        "gamma_site_photo",    # 活性位面密度 (m^-2)，把光子流折算成"每位点"
        "qfl_split_eV",        # given_qfl 模式：直接给 ΔE_qfl
        "qfl_frac_n",          # ★ 分裂中落在电子一侧的比例 f（0.5 = 本征/高注入）
        "recomb_rate_const",   # 双分子（辐射/带带）复合系数，每位点每秒
        "recomb_srh_const",    # 单分子（陷阱/SRH）复合系数
        "beta_photo",          # 光生驱动力进入势垒的传递系数
        "photo_carrier_mode",  # 'given_qfl' | 'steady_carrier'
        "qfl_max_eV",          # 分裂上限（物理上不超过带隙），默认 photon_energy_eV
    )
    step_params = ("nu_elec_cb", "nu_hole_vb")
    #: ★ 设计修正（T12 实现时发现）：**自洽放在求解器外层，不放进 modify_energies**。
    #   若在 driver 内部重算界面消耗速率，就得自己拿前因子、气相活度……
    #   那是**隐式耦合**，极易和内核算的不一致，而且不一致时不会报错。
    #   ⇒ driver 只管"给定分裂 → 改能量"（与 electro 完全同构），
    #     载流子收支由 `photo_steady()` 在外层做不动点迭代，**复用已验证的内核**。
    needs_theta = False
    is_additive = True
    provides = ("AQY", "qfl_eV", "carrier_n_rel", "carrier_p_rel", "G_photo", "R_recomb")

    # ── 内部：把 x = ΔE/kT 换成能量修正 ──
    @staticmethod
    def _shift(mech, Ga_fwd, dG_rxn, dE_n, dE_p, beta):
        """δE_Fn 抬高电子的化学势；δE_Fp 降低空穴侧。与 CHE 同构。"""
        Ga = np.asarray(Ga_fwd, float).copy()
        dG = np.asarray(dG_rxn, float).copy()
        for sid, step in enumerate(mech.steps):
            n_cb = float(step.params.get("nu_elec_cb", 0.0))
            n_vb = float(step.params.get("nu_hole_vb", 0.0))
            if n_cb == 0.0 and n_vb == 0.0:
                continue
            drive = n_cb * dE_n + n_vb * dE_p
            dG[sid] -= drive
            Ga[sid] -= beta * drive
        return Ga, dG

    def _params(self, state):
        f = float(self.get_param(state, "qfl_frac_n", 0.5))
        if not 0.0 <= f <= 1.0:
            raise ValueError("qfl_frac_n 要在 [0,1]：0 = 全部落在空穴侧（重 n 掺杂），"
                             "0.5 = 本征/高注入，1 = 全部落在电子侧（重 p 掺杂）。收到 %g" % f)
        # ★★ 默认 1.0，不是电催化那个 0.5 —— 这是 T12 实测暴露出来的一个真物理错。
        #
        # 电催化取 β=0.5 的理由是"电子从电极转移、过渡态在中途"。
        # 但**界面载流子俘获的标准动力学形式是 `rate = k·n·θ`（对载流子浓度一级）**，
        # 而载流子浓度 n/n_dark = exp(f·ΔE/k_BT)。要让速率律复现这个一级依赖，
        # 驱动力必须**全额**进入势垒，即 β = 1：
        #     rate ∝ exp(β·ν_cb·f·ΔE/k_BT)  --β=1,ν=1,f=0.5-->  exp(0.5·x) = n/n_dark ✓
        # 取 0.5 会把载流子的质量作用依赖**少算一半**，于是高光强幂指数算成 0.25 而非 0.5。
        #
        # β < 1 留给"载流子俘获之外还有额外活化势垒"的情形，但那要用户明确给。
        beta = float(self.get_param(state, "beta_photo", 1.0))
        T = float(self.get_param(state, "T_K", required=True))
        return f, beta, K_B_EV * T

    def _generation(self, state):
        flux = self.get_param(state, "photon_flux")
        if flux is None:
            raise ValueError(
                "steady_carrier 模式需要 photon_flux（光子流, m⁻² s⁻¹）。\n"
                "  只想给定分裂就用 photo_carrier_mode='given_qfl' + qfl_split_eV。")
        absb = float(self.get_param(state, "absorptance", 1.0))
        gs = self.get_param(state, "gamma_site_photo")
        if not gs:
            raise ValueError(
                "steady_carrier 模式需要 gamma_site_photo（活性位面密度, m⁻²）——\n"
                "  光子流是**每面积**的，表面速率是**每位点**的，没有它两者没法比。")
        return float(absb) * float(flux) / float(gs)

    @staticmethod
    def _recomb(x, k_bi, k_srh, f):
        """R_rec(x)，**暗态（x=0）精确为零**。"""
        if x <= 0.0:
            return 0.0
        em1 = np.expm1(min(x, 700.0))                 # e^x − 1，小 x 也准
        out = k_bi * em1
        if k_srh:
            denom = np.exp(min(f * x, 700.0)) + np.exp(min((1.0 - f) * x, 700.0))
            out += k_srh * em1 / denom
        return float(out)

    def modify_energies(self, mech, Ga_fwd, dG_rxn, state, theta=None):
        mode = str(self.get_param(state, "photo_carrier_mode", "given_qfl"))
        f, beta, kT = self._params(state)

        if mode == "given_qfl":
            dq = self.get_param(state, "qfl_split_eV")
            if dq is None:
                raise ValueError("'given_qfl' 模式需要 qfl_split_eV（准费米能级分裂, eV）。")
            dq = float(dq)
            return self._shift(mech, Ga_fwd, dG_rxn, f * dq, (1.0 - f) * dq, beta)

        if mode != "steady_carrier":
            raise ValueError("photo_carrier_mode 只能是 'given_qfl' 或 'steady_carrier'，"
                             "收到 %r" % mode)
        raise RuntimeError(
            "'steady_carrier' 不能直接走 steady_state —— 载流子收支要在**外层**自洽。\n"
            "  用： from kinstep.drivers.photo import photo_steady\n"
            "       res = photo_steady(mech, en, state, ...)\n"
            "  （这样复用已验证的内核算速率，而不是在 driver 里重算一遍 —— "
            "重算的版本一旦与内核不一致，是不会报错的。）")

    def observables(self, mech, solution, state):
        out = {}
        last = getattr(self, "_last", None)
        rate = np.asarray(solution["rate_net"], float)
        n_e = sum(float(s.params.get("nu_elec_cb", 0.0)) * rate[i]
                  for i, s in enumerate(mech.steps))
        if last:
            f, x, kT = last["f"], last["x"], last["kT"]
            out["qfl_eV"] = x * kT
            out["carrier_n_rel"] = float(np.exp(min(f * x, 700.0)))
            out["carrier_p_rel"] = float(np.exp(min((1.0 - f) * x, 700.0)))
            out["G_photo"] = last["G"]
            out["R_recomb"] = last["R_rec"]
            out["AQY"] = abs(n_e) / last["G"] if last["G"] else float("nan")
        else:
            dq = self.get_param(state, "qfl_split_eV")
            if dq is not None:
                out["qfl_eV"] = float(dq)
            flux = self.get_param(state, "photon_flux")
            gs = self.get_param(state, "gamma_site_photo")
            if flux and gs:
                absb = float(self.get_param(state, "absorptance", 1.0))
                G = absb * float(flux) / float(gs)
                out["G_photo"] = G
                out["AQY"] = abs(n_e) / G if G else float("nan")
        return out

    def validate(self, mech):
        if not any(s.params.get("nu_elec_cb") or s.params.get("nu_hole_vb")
                   for s in mech.steps):
            raise ValueError(
                "启用了 photo driver，但没有任何一步标了 nu_elec_cb / nu_hole_vb。\n"
                "  在机理里这样写： A* + e-(cb) -> A-*   {nu_elec_cb: 1}")


# ══════════════════════ 外层：载流子收支自洽 ══════════════════════
def photo_steady(mech, en, state, p_gas=None, c_liq=None, drivers=None,
                 max_iter=60, tol=1e-10, damping=0.5, verbose=False, **kw):
    r"""解 **G = R_rec(x) + R_surf(x, θ)**，返回稳态解（x = ΔE_qfl/k_BT）。

    为什么用不动点而不是把它塞进 `modify_energies`
    ----------------------------------------------
    界面消耗 `R_surf` 要用到**内核算出来的速率**。在 driver 里重算一遍，
    就得自己拿前因子、气相/液相活度、通量投影…… 任何一处与内核不一致都**不会报错**，
    只会悄悄给出一个看着合理的 ΔE_qfl。⇒ 外层迭代，**内核算速率，driver 只改能量**。

    收敛性：`R_rec(x) + R_surf(x)` 对 x 单调增（复合随分裂增大；界面消耗也是），
    所以根唯一；用二分法定位、再带阻尼更新，不依赖初值。
    """
    import numpy as _np
    from ..driver import DriverStack, make_driver
    from ..solve import steady_state
    from ..units import K_B_EV

    ph = None
    stack = list(drivers or [])
    for d in stack:
        if getattr(d, "driver_name", "") == "photo":
            ph = d
    if ph is None:
        ph = make_driver("photo")
        stack.append(ph)
    f, beta, kT = ph._params(state)
    G = ph._generation(state)
    k_bi = float(ph.get_param(state, "recomb_rate_const", 0.0))
    k_srh = float(ph.get_param(state, "recomb_srh_const", 0.0))
    if k_bi <= 0 and k_srh <= 0:
        raise ValueError(
            "载流子收支自洽至少要给一个复合系数（recomb_rate_const 双分子 / "
            "recomb_srh_const 单分子）。\n"
            "  全不给 = 载流子永不复合 ⇒ 分裂会一直涨到带隙上限，**那不是物理**。")
    x_max = float(ph.get_param(state, "qfl_max_eV",
                               ph.get_param(state, "photon_energy_eV", 3.0))) / kT
    nu_cb = _np.array([float(s.params.get("nu_elec_cb", 0.0)) for s in mech.steps])
    nu_vb = _np.array([float(s.params.get("nu_hole_vb", 0.0)) for s in mech.steps])
    nu_tot = nu_cb + nu_vb

    def solve_at(x, guess=None):
        st = dict(state)
        st["photo_carrier_mode"] = "given_qfl"
        st["qfl_split_eV"] = x * kT
        return steady_state(mech, en, state=st,
                            driver=DriverStack(stack, warn_order=False),
                            p_gas=p_gas, c_liq=c_liq, theta_guess=guess,
                            ladder=("newton",) if guess else ("ode", "newton"),
                            strict_energy=False, **kw)

    def drain_at(x, guess=None):
        r = solve_at(x, guess)
        if not r.converged:
            r = solve_at(x, None)
        return float(_np.sum(nu_tot * r.rate_net)), r

    d0, r0 = drain_at(0.0)
    if ph._recomb(0.0, k_bi, k_srh, f) + d0 >= G:
        if verbose:
            print("   暗态界面消耗已 ≥ G ⇒ ΔE_qfl = 0")
        return r0, 0.0, dict(G=G, R_rec=0.0, R_surf=d0, n_iter=0,
                             qfl_eV=0.0, AQY=abs(d0) / G if G else float("nan"),
                             R_elec=0.0, R_hole=0.0, charge_imbalance=0.0)
    guess = dict(r0.theta)

    # ★ 先用"只有复合"的方程算上界 —— 纯标量，**一次内核都不调**。
    #   因为界面消耗 R_surf ≥ 0，所以真根一定 ≤ 这个 x_rec。
    #   不做这一步就得从 [0, x_max=带隙/kT≈120] 开始二分，而根可能在 0.003 ——
    #   实测 44 层二分里有 30 多层纯粹在缩区间上浪费（每层都要解一次稳态）。
    lo_s, hi_s = 0.0, x_max
    for _ in range(200):
        m = 0.5 * (lo_s + hi_s)
        if ph._recomb(m, k_bi, k_srh, f) >= G:
            hi_s = m
        else:
            lo_s = m
        if hi_s - lo_s < 1e-14 * max(1.0, hi_s):
            break
    x_rec = hi_s
    if verbose:
        print("   仅复合时的上界 x_rec = %.6g（ΔE = %.6f eV）" % (x_rec, x_rec * kT))

    # Brent 求根：同样的括号下比纯二分少用 4~5 倍的函数求值
    from scipy.optimize import brentq
    calls = {"n": 0, "res": r0, "guess": guess}

    def bal(x):
        calls["n"] += 1
        d, r = drain_at(x, calls["guess"])
        if r.converged:
            calls["guess"] = dict(r.theta)
            calls["res"] = r
        v = ph._recomb(x, k_bi, k_srh, f) + d - G
        if verbose:
            print("   #%2d  x=%.8g (ΔE=%.6f eV)  R_rec=%.4e  R_surf=%.4e  失衡 %+.3e"
                  % (calls["n"], x, x * kT, ph._recomb(x, k_bi, k_srh, f), d, v))
        return v

    f_hi = bal(x_rec)
    if f_hi < 0:                      # 上界处仍不足 ⇒ 根在带隙上限
        x = x_max
    else:
        x = brentq(bal, 0.0, x_rec, xtol=tol * max(1.0, x_rec), rtol=max(tol, 1e-15),
                   maxiter=max_iter)
    it = calls["n"]
    d, res = drain_at(x, calls["guess"])
    # ★ 免费自检：**电子与空穴的界面消耗速率必须相等**（光催化剂不能净充电）。
    #   与 Campbell 求和规则同性质：不依赖能量、不依赖机理细节，永远成立。
    #   不相等 ⇒ 机理里的 nu_elec_cb / nu_hole_vb 配错了（这是最常见的一类错，
    #   而且**不会报错**，只会给出一条看着合理的曲线）。
    R_e = float(_np.sum(nu_cb * res.rate_net))
    R_h = float(_np.sum(nu_vb * res.rate_net))
    scale = max(abs(R_e), abs(R_h), 1e-300)
    charge_imbalance = abs(R_e - R_h) / scale
    info = dict(G=G, R_rec=ph._recomb(x, k_bi, k_srh, f), R_surf=d, n_iter=it,
                qfl_eV=x * kT, AQY=abs(d) / G if G else float("nan"),
                R_elec=R_e, R_hole=R_h, charge_imbalance=charge_imbalance)
    if charge_imbalance > 1e-6:
        info["warning"] = (
            "⚠ 电子与空穴的消耗速率不相等（%.3e vs %.3e，相对失衡 %.2e）——\n"
            "   光催化剂会净充电，这在稳态下不可能。多半是某一步的 "
            "nu_elec_cb / nu_hole_vb 配错了。" % (R_e, R_h, charge_imbalance))
    ph._last = dict(x=x, G=G, kT=kT, f=f, R_rec=info["R_rec"])
    return res, x, info
