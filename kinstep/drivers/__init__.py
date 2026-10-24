# -*- coding: utf-8 -*-
# KinSTEP - a modality-agnostic microkinetic solver.
# Copyright (C) 2026 The KinSTEP authors.
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.  It is distributed WITHOUT ANY WARRANTY; see the GNU General Public
# License (LICENSE) for details.
"""[English] Drivers for the catalytic modalities; importing this package registers them.

(The commentary below is in Chinese: it records why each choice was made
and which measurement forced it.)

各催化模态的 driver。import 本包即完成注册（并触发撞名检查）。"""
from . import thermal, electro, photo, sono, field_eff   # noqa: F401
