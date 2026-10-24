# -*- coding: utf-8 -*-
# KinSTEP - a modality-agnostic microkinetic solver.
# Copyright (C) 2026 The KinSTEP authors.
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.  It is distributed WITHOUT ANY WARRANTY; see the GNU General Public
# License (LICENSE) for details.
"""[English] Thermal baseline driver: changes nothing, making 'no stimulus' explicit.

(The commentary below is in Chinese: it records why each choice was made
and which measurement forced it.)

热催化：基线 driver —— 不改能量，只是把"什么都不加"显式化。

存在的意义：让"纯热催化"也走同一条 DriverStack 代码路径，
避免内核里出现 `if driver is None:` 这种分支（分支越少，错越少）。
"""
import numpy as np
from ..driver import Driver, register_driver


@register_driver
class ThermalDriver(Driver):
    driver_name = "thermal"
    state_params = ()
    step_params = ()
    needs_theta = False
    is_additive = True
    provides = ()

    def modify_energies(self, mech, Ga_fwd, dG_rxn, state, theta=None):
        return np.asarray(Ga_fwd, float).copy(), np.asarray(dG_rxn, float).copy()
