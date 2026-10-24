# -*- coding: utf-8 -*-
# KinSTEP - a modality-agnostic microkinetic solver.
# Copyright (C) 2026 The KinSTEP authors.
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.  It is distributed WITHOUT ANY WARRANTY; see the GNU General Public
# License (LICENSE) for details.
"""[English] Physical constants and units, with the CODATA generation stated explicitly.

(The commentary below is in Chinese: it records why each choice was made
and which measurement forced it.)

物理常数与单位。

口径与 Ager 的 `EC CO stripping MKM.ipynb` 对齐（T2 已逐个对数），
这样 Koper 的复现不会因为常数差异而对不上。

★ 全部用 **SI**，但能量对外一律以 **eV** 表达（催化界惯例），
  内部换算只在一处（`eV`），避免到处乘除出错。
"""
import math

# ── 基本常数（CODATA-2018 精确值）──
#
# ★ 2026-10-22 修正一个潜伏 bug：原先 `Q_E = 1.60217662e-19`（抄 notebook 的截断值）
#   而 `EV = 1.602176634e-19`，**两者按定义必须是同一个数**（1 eV ≡ 1 e × 1 V）。
#   两个值相差 8.7e-9，看着无害，但电化学里要乘 U/k_BT ≈ 19.5 的放大 ⇒ 8.5e-8，
#   而且**只在电势一大就显出来**，单点测试永远看不见。
#   ⇒ 现在只定义一次，`EV = Q_E`，物理上不可能再不自洽。
# ★ 2026-10-22 第二次修同类错（T7 的 401 点回归逼出来的）：
#   原先 `R_GAS = 8.3144598` 是 **CODATA-2014** 的值，而 N_A / H_PLANCK / Q_E
#   用的是 **CODATA-2018 精确值** —— 又混了两代。
#   后果：K_B = R_GAS/N_A = 1.3806485e-23，而 2018 的精确值是 1.380649e-23，
#   相对差 3.4e-07。本身很小，但 `Ea/k_BT ≈ 25` 会把它放大到 ~1e-5 的 kf 偏差。
#   ⇒ 用 2018 的 R，使 K_B = R_GAS/N_A **恰好等于** 1.380649e-23。
#   教训与 Q_E/EV 那次完全同类：**常数要整代取，不许各抄一处**。
R_GAS = 8.31446261815324   # J mol^-1 K^-1（CODATA-2018：= N_A·K_B 精确）
Q_E = 1.602176634e-19      # C，元电荷（CODATA-2018 定义值）
EV = Q_E                   # J，1 eV —— **同一个数，不许各写一遍**
N_A = 6.02214076e23        # mol^-1（定义值）
F_FARADAY = N_A * Q_E      # C mol^-1
K_B = R_GAS / N_A          # J K^-1
H_PLANCK = 6.62607015e-34  # J s（定义值）

# ── 常用导出量 ──
K_B_EV = K_B / EV          # eV K^-1，8.617e-5
LN10 = math.log(10.0)

# ── 标准态（与 CATKINAS 一致：压力除以 P0，浓度除以 C0）──
P0_PA = 1.0e5              # 0.1 MPa
C0_MOL_L = 1.0

# ── 单位换算便捷量 ──
CM = 1e-2                  # m
MA = 1e-3                  # A
MV = 1e-3                  # V


def kT_eV(T):
    """k_B·T，单位 eV。"""
    return K_B_EV * T


def thermal_voltage(T):
    """k_B·T/e，单位 V（= kT_eV 的数值，但语义不同，分开写以免误用）。"""
    return K_B_EV * T


def tafel_ideal(T, beta=0.5, n=1):
    """理想 Tafel 斜率 (V/decade) = ln10·k_BT/(β·n·e)。298 K、β=0.5、n=1 时 ≈ 0.118 V/dec。"""
    return LN10 * K_B_EV * T / (beta * n)


def U_she_to_rhe(U_she, pH, T=298.15):
    """SHE → RHE：U_RHE = U_SHE + (k_BT·ln10/e)·pH。"""
    return U_she + LN10 * K_B_EV * T * pH


def U_rhe_to_she(U_rhe, pH, T=298.15):
    return U_rhe - LN10 * K_B_EV * T * pH
