# -*- coding: utf-8 -*-
# KinSTEP - a modality-agnostic microkinetic solver.
# Copyright (C) 2026 The KinSTEP authors.
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.  It is distributed WITHOUT ANY WARRANTY; see the GNU General Public
# License (LICENSE) for details.
"""[English] KinSTEP: microkinetic modelling across thermal, electro-, photo- and
sonocatalysis in pure Python.

(The commentary below is in Chinese: it records why each choice was made
and which measurement forced it.)

KinSTEP：微动力学模拟（热 / 电 / 光 / 超声 …），纯 Python。"""
from . import units, naming, driver        # noqa: F401
from . import drivers                      # noqa: F401  导入即注册所有 driver
from .driver import (available_drivers, make_driver, param_owner_table,   # noqa: F401
                     Driver, DriverStack, register_driver)

__version__ = "0.1.0.dev"
