# -*- coding: utf-8 -*-
# KinSTEP - a modality-agnostic microkinetic solver.
# Copyright (C) 2026 The KinSTEP authors.
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.  It is distributed WITHOUT ANY WARRANTY; see the GNU General Public
# License (LICENSE) for details.
r"""[English] Electrocatalysis driver: computational hydrogen electrode and Butler-Volmer.

(The commentary below is in Chinese: it records why each choice was made
and which measurement forced it.)

电催化 driver —— CATKINAS 完全没有的那一半。

符号约定（**全包唯一**，见 naming.py）
--------------------------------------
基元步里**显式写出 e-**；令 `nu_elec` = e- 的计量系数，**产物为正 = 氧化**。

    ΔG(U)   = ΔG(0)   − nu_elec · e · U_elec          （CHE）
    Ga_f(U) = Ga_f(0) − beta_ct · nu_elec · e · U_elec （Butler–Volmer）
    Ga_r(U) = Ga_r(0) + (1−beta_ct) · nu_elec · e · U_elec
              —— 由 Ga_r = Ga_f − ΔG 自动保证，不单独给，**免得两者不自洽**

    j_curr  = Σ_steps nu_elec · F · gamma_site · rate_net        (A/m²)
    U_RHE   = U_SHE + (k_BT·ln10/e)·pH_val

与 Koper(1998)/Ager notebook 的对齐校验
---------------------------------------
    (1) H2O + * ⇌ OH* + H+ + e-   (nu_elec = +1)
    ⇒ Ga_f(U) = Ga_f(0) − beta_ct·e·U  ⇒ k1 ∝ exp(+beta_ct·e·U/k_BT)   ✅
    ⇒ Ga_r(U) = Ga_r(0) + (1−beta_ct)·e·U ⇒ k-1 ∝ exp(−(1−beta_ct)e·U/k_BT) ✅
  与 notebook 的 k1 / kminus1 逐字一致。

⚠ `beta_ct` **按步绑定**，不设全局默认覆盖 —— notebook 里 `k2` 误用 `α_1`
  正是因为用了全局名而没人发现（两者恰好相等）。
"""
import numpy as np

from ..driver import Driver, register_driver
from ..units import K_B_EV, LN10, F_FARADAY


@register_driver
class ElectroDriver(Driver):
    driver_name = "electro"
    #: 工况参数：电极电势、参比、pH、位点面密度（算电流用）
    state_params = ("U_elec", "U_ref", "pH_val", "gamma_site")
    #: 逐步参数：电子计量系数、电荷传递系数
    step_params = ("nu_elec", "beta_ct")
    #: 这几个保留名本来就是电化学的，driver 是它们的主人
    owns_reserved = ("U_elec", "U_ref", "pH_val", "gamma_site", "nu_elec", "beta_ct")
    #: CHE/BV 的修正对 Ga/dG 是**线性可加**的 ⇒ 与别的可加 driver 换序结果不变
    is_additive = True
    #: 能量不依赖覆盖度（若将来要加覆盖度依赖的双电层效应，改成 True 即可，内核自动改走内层）
    needs_theta = False
    provides = ("j_curr", "U_she")

    def _U_she(self, state):
        """统一折算到 SHE —— 内部一律用 SHE，避免 RHE/SHE 混用算错。"""
        U = float(self.get_param(state, "U_elec", required=True))
        ref = str(self.get_param(state, "U_ref", "SHE")).upper()
        if ref == "SHE":
            return U
        if ref == "RHE":
            return U - LN10 * K_B_EV * float(self.get_param(state, "T_K", 298.15)) \
                       * float(self.get_param(state, "pH_val", 0.0))
        raise ValueError("U_ref 只能是 'SHE' 或 'RHE'，收到 %r" % ref)

    def modify_energies(self, mech, Ga_fwd, dG_rxn, state, theta=None):
        Ga = np.asarray(Ga_fwd, float).copy()
        dG = np.asarray(dG_rxn, float).copy()
        U = self._U_she(state)
        for sid, step in enumerate(mech.steps):
            nu = float(step.params.get("nu_elec", 0.0))
            if nu == 0.0:
                continue                      # 化学步，不受电势影响
            beta = float(step.params["beta_ct"])
            dG[sid] -= nu * U                 # eV：e·U 的数值就是 U(V)
            Ga[sid] -= beta * nu * U
        return Ga, dG

    def observables(self, mech, solution, state):
        """电流密度与（若给了扫描点）Tafel 斜率。"""
        gs = float(self.get_param(state, "gamma_site", 0.0))
        rate = np.asarray(solution["rate_net"], float)
        j = 0.0
        for sid, step in enumerate(mech.steps):
            nu = float(step.params.get("nu_elec", 0.0))
            if nu:
                j += nu * F_FARADAY * gs * rate[sid]
        return {"j_curr": j, "U_she": self._U_she(state)}

    def validate(self, mech):
        n_ec = 0
        for sid, step in enumerate(mech.steps):
            nu = step.params.get("nu_elec", 0.0)
            if nu:
                n_ec += 1
                if "beta_ct" not in step.params:
                    raise ValueError(
                        "第 %d 步是电化学步（nu_elec=%s）却没给 beta_ct。\n"
                        "⚠ 本包**不设全局默认**：Ager 的 notebook 正因全局 α 而把 "
                        "k2 的 α_2 写成 α_1 且无人察觉。" % (sid + 1, nu))
                b = float(step.params["beta_ct"])
                if not 0.0 < b < 1.0:
                    raise ValueError("第 %d 步 beta_ct=%g 不在 (0,1)" % (sid + 1, b))
        if n_ec == 0:
            raise ValueError("启用了 electro driver，但没有任何一步含 e- —— 机理里忘了写电子？")
