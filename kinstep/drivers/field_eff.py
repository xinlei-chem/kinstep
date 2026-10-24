# -*- coding: utf-8 -*-
# KinSTEP - a modality-agnostic microkinetic solver.
# Copyright (C) 2026 The KinSTEP authors.
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.  It is distributed WITHOUT ANY WARRANTY; see the GNU General Public
# License (LICENSE) for details.
r"""[English] Static-field driver: dipole and polarizability shift of each step; quadratic in
the field and therefore non-additive when stacked with another driver.

(The commentary below is in Chinese: it records why each choice was made
and which measurement forced it.)

外电场 / 界面电场 driver（可与 electro 叠加）。

    ΔG(F) = ΔG(0) − Δmu · F − 0.5 · Δalpha · F²

`Δmu`（偶极矩变化，e·Å）与 `Δalpha`（极化率变化，e·Å²/V）逐步给定。
这是研究20"曲率改变界面电场→改变吸附"那条链在微动力学里的落点。

★ 公式已实现（一次+二次项），但**Δmu/Δalpha 要用户自己从 DFT 算**，
  本 driver 不代为估计。
"""
import numpy as np

from ..driver import Driver, register_driver


@register_driver
class FieldDriver(Driver):
    driver_name = "field"
    state_params = ("field_V_per_A",)
    step_params = ("dmu_eA", "dalpha_eA2_per_V")
    #: ★ 含 F² 二次项 ⇒ **非可加**，与别的 driver 叠加时顺序影响结果，Stack 会警告
    is_additive = False
    needs_theta = False
    provides = ()
    owns_reserved = ("field_V_per_A",)

    def modify_energies(self, mech, Ga_fwd, dG_rxn, state, theta=None):
        Ga = np.asarray(Ga_fwd, float).copy()
        dG = np.asarray(dG_rxn, float).copy()
        Fv = float(self.get_param(state, "field_V_per_A", required=True))
        for sid, step in enumerate(mech.steps):
            dmu = float(step.params.get("dmu_eA", 0.0))
            dal = float(step.params.get("dalpha_eA2_per_V", 0.0))
            if dmu == 0.0 and dal == 0.0:
                continue
            shift = dmu * Fv + 0.5 * dal * Fv * Fv
            dG[sid] -= shift
            Ga[sid] -= 0.5 * shift      # 过渡态按半程近似；需要更准时用户自己给 TS 的 μ
        return Ga, dG

    def validate(self, mech):
        return None
