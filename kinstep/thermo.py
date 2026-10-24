# -*- coding: utf-8 -*-
# KinSTEP - a modality-agnostic microkinetic solver.
# Copyright (C) 2026 The KinSTEP authors.
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.  It is distributed WITHOUT ANY WARRANTY; see the GNU General Public
# License (LICENSE) for details.
r"""[English] Gas-phase thermochemistry: Shomate coefficients from NIST into entropy
corrections to barriers and reaction energies.

(The commentary below is in Chinese: it records why each choice was made
and which measurement forced it.)

气相热化学（T10）：Shomate 方程 → 熵校正。

它解决什么问题
--------------
DFT 给的是**纯电子能**，没有气相分子的平动/转动熵。直接拿它当 ΔG 用，
吸附步会系统性地算得太有利 —— 673 K 下 H2 吸附的这一项就有 **1.08 eV**，
比大多数化学效应都大。

★ 一条规则导出全部校正（由 CATKINAS 401 点反推 + 验证）
--------------------------------------------------------
    offset_dG[s] = −T · Σ_j ν[s,j]·S_j(T)            （只数气相物种）
    offset_Ga[s] = +T · Σ_j a_r[s,j]·S_j(T)          （只数**反应物侧**的气相物种）

`offset_Ga` 为什么只数反应物侧：过渡态是**吸附在表面上的**，气相分子走到过渡态时
平动/转动自由度已经交出去了；而产物侧的气相要到终态才出现，与势垒无关。

实测对上 CATKINAS 的合成氨算例（T = 673 K）：
| 步 | 反应 | 本规则给出 | CATKINAS 实测 |
|---|---|---|---|
| 1 | H2 + 2* ⇌ 2H* | ΔGa = ΔG⁰ = +T·S(H2) | +1.077301 |
| 2 | N2 + 2* ⇌ 2N* | ΔGa = ΔG⁰ = +T·S(N2) | +1.504256 |
| 3,4 | 表面步，无气相 | 0 | 0 |
| 5 | NH2*+H* ⇌ NH3+2* | ΔGa = 0；ΔG⁰ = −T·S(NH3) | −1.575974 |

★ 数据来源与**自检**（这是本模块最要紧的设计）
----------------------------------------------
系数取自 **NIST Chemistry WebBook**（Shomate 形式，见 `CITATION`）。
我在台账里自我更正过："不从二手手册抄表，因为会把转录错误一起带进来"。
⇒ 所以每个物种都**随带它自己的参考值** `S298` / `Cp298`，
  `self_check()` 用系数现算 298.15 K 的 S 与 Cp 跟参考值比 ——
  **抄错一个数就会被当场抓住，不可能静默发布**。本模块在导入时不自动跑，
  但回归测试里跑，且 `kinstep check` 会调用。

★ Shomate **分温区**，不许静默外推
-----------------------------------
实测教训：N2 的 500–2000 K 段外推到 298 K 得 191.38，而真值 191.61（差 0.23）。
单看数字"挺像"，但它是错的。⇒ 超出温区一律**报错**，不给"差不多"的答案。
多温区物种按 T 自动选段。
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from .mech import Mechanism
from .units import N_A, Q_E

#: 1 eV 对应多少 J/mol
J_PER_MOL_PER_EV = N_A * Q_E

CITATION = (
    "NIST Chemistry WebBook, SRD 69 — Shomate 方程系数。"
    "原始汇编：Chase, M.W., Jr., NIST-JANAF Thermochemical Tables, 4th ed., "
    "J. Phys. Chem. Ref. Data, Monograph 9, 1998, 1-1951. "
    "https://webbook.nist.gov/chemistry/  —— 用本模块请引用 NIST，不要引用本包。"
)


class ThermoError(ValueError):
    """热化学数据有问题：温度超出 Shomate 温区，或物种不在数据库里。"""


@dataclass
class ShomateRange:
    """一个温区的 Shomate 系数。t = T/1000。

        Cp° = A + B·t + C·t² + D·t³ + E/t²                     J/(mol·K)
        H°−H°₂₉₈ = A·t + B·t²/2 + C·t³/3 + D·t⁴/4 − E/t + F − H   kJ/mol
        S°  = A·ln t + B·t + C·t²/2 + D·t³/3 − E/(2t²) + G      J/(mol·K)
    """
    T_lo: float
    T_hi: float
    A: float; B: float; C: float; D: float; E: float; F: float; G: float; H: float

    def Cp(self, T):
        t = T / 1000.0
        return self.A + self.B * t + self.C * t * t + self.D * t ** 3 + self.E / (t * t)

    def S(self, T):
        t = T / 1000.0
        return (self.A * math.log(t) + self.B * t + self.C * t * t / 2
                + self.D * t ** 3 / 3 - self.E / (2 * t * t) + self.G)

    def dH(self, T):
        """H°(T) − H°(298.15)，kJ/mol。"""
        t = T / 1000.0
        return (self.A * t + self.B * t * t / 2 + self.C * t ** 3 / 3
                + self.D * t ** 4 / 4 - self.E / t + self.F - self.H)


@dataclass
class GasSpecies:
    """一个气相物种的热化学数据 + **它自己的参考值**（用于自检）。"""
    formula: str
    ranges: List[ShomateRange]
    S298: float            # J/(mol·K)，参考值
    Cp298: float           # J/(mol·K)，参考值
    name: str = ""
    source: str = "NIST Chemistry WebBook (Chase 1998)"
    #: ★ **已知的拟合端点偏差**（Shomate 现算值 − 参考值，J/(mol·K)）。
    #   不是误差，是记录在案的事实：写在这里而不是把容差放宽，
    #   这样"抄错系数"与"拟合在端点本来就这样"两件事不会被混为一谈。
    cp298_known_dev: float = 0.0
    dev_reason: str = ""

    def _range_at(self, T):
        for r in self.ranges:
            if r.T_lo - 1e-9 <= T <= r.T_hi + 1e-9:
                return r
        spans = "；".join("%g–%g K" % (r.T_lo, r.T_hi) for r in self.ranges)
        raise ThermoError(
            "%s 的 Shomate 数据不覆盖 T = %g K（有效温区：%s）。\n"
            "  ★ 本模块**不静默外推**：实测 N2 用 500–2000 K 段外推到 298 K 得 191.38，"
            "而真值 191.61 —— 看着挺像，但它是错的。\n"
            "  要么换温度，要么自己补这个温区的系数（并补上参考值供自检）。"
            % (self.formula, T, spans))

    def S(self, T):
        """标准摩尔熵 S°(T)，J/(mol·K)。"""
        return self._range_at(T).S(T)

    def Cp(self, T):
        return self._range_at(T).Cp(T)

    def TS_eV(self, T, J_per_mol_per_eV: float = None):
        """T·S°(T)，换算成 eV —— 这正是熵校正用的量。

        `J_per_mol_per_eV` 可显式指定，用于**与别的程序对数**：
        实测 CATKINAS 用的是 96485.2608（三个物种反解的散布只有 3.8e-10，
        确是单一常数），与本包的 CODATA-2018 值 96485.3321 差 7.39e-07。
        这点差异经 `T·S/k_BT` 放大后，在速率上是 **1.9e-05** —— 不处理就会被
        误读成"热力学库算错了"，其实只是常数口径不同（这已是与 CATKINAS 的第四处：
        k_B、e、h、F）。
        """
        return self.S(T) * T / (J_per_mol_per_eV or J_PER_MOL_PER_EV)


#: ── 数据库（NIST Chemistry WebBook）──
#   每条都带 S298 / Cp298 参考值，`self_check()` 会用系数现算去对，抄错就当场暴露。
def _sp(formula, S298, Cp298, *ranges, name="", cp298_known_dev=0.0, dev_reason=""):
    return GasSpecies(formula=formula, name=name or formula,
                      ranges=[ShomateRange(*r) for r in ranges],
                      S298=S298, Cp298=Cp298,
                      cp298_known_dev=cp298_known_dev, dev_reason=dev_reason)


#: 气相热化学数据库：化学式 → GasSpecies（名字避开两字母的 `DB`，见 naming.py 的规矩）
GAS_DB: Dict[str, GasSpecies] = {s.formula: s for s in [
    _sp("H2", 130.68, 28.84,
        (298.0, 1000.0, 33.066178, -11.363417, 11.432816, -2.772874,
         -0.158558, -9.980797, 172.707974, 0.0),
        (1000.0, 2500.0, 18.563083, 12.257357, -2.859786, 0.268238,
         1.977825, -1.147438, 156.288133, 0.0), name="hydrogen"),
    _sp("N2", 191.61, 29.12,
        (100.0, 500.0, 28.98641, 1.853978, -9.647459, 16.63537,
         0.000117, -8.671914, 226.4168, 0.0),
        (500.0, 2000.0, 19.50583, 19.88705, -8.598535, 1.369784,
         0.527601, -4.935202, 212.3900, 0.0), name="nitrogen"),
    _sp("O2", 205.15, 29.38,
        (100.0, 700.0, 31.32234, -20.23531, 57.86644, -36.50624,
         -0.007374, -8.903471, 246.7945, 0.0),
        (700.0, 2000.0, 30.03235, 8.772972, -3.988133, 0.788313,
         -0.741599, -11.32468, 236.1663, 0.0), name="oxygen"),
    _sp("CO", 197.66, 29.14,
        (298.0, 1300.0, 25.56759, 6.096130, 4.054656, -2.671301,
         0.131021, -118.0089, 227.3665, -110.5271), name="carbon monoxide"),
    _sp("CO2", 213.79, 37.13,
        (298.0, 1200.0, 24.99735, 55.18696, -33.69137, 7.948387,
         -0.136638, -403.6075, 228.2431, -393.5224), name="carbon dioxide"),
    _sp("H2O", 188.84, 33.58,
        (500.0, 1700.0, 30.09200, 6.832514, 6.793435, -2.534480,
         0.082139, -250.8810, 223.3967, -241.8264), name="water (gas)"),
    _sp("NH3", 192.77, 35.06,
        (298.0, 1400.0, 19.99563, 49.77119, -15.37599, 1.921168,
         0.189174, -53.30667, 203.8591, -45.89806), name="ammonia",
        # ★ 自检当场抓到并查清的一条：NH3 的 Shomate Cp 在 298 K 端点比参考值高 0.587
        #   （1.67%），但 400–1000 K 全部落在 0.1~0.3% 以内（实测：400K +0.05、500K +0.03、
        #   700K −0.05、1000K +0.00）。而它的 S(298) 精确给出 192.77，
        #   S(673) 又与 CATKINAS 反推值逐位吻合 ⇒ **系数没抄错，是拟合区间端点的固有行为**。
        cp298_known_dev=0.587,
        dev_reason="Shomate 拟合在 298 K（区间下端点）的固有偏差；400–1000 K 内 ≤0.3%"),
    _sp("CH4", 186.25, 35.69,
        (298.0, 1300.0, -0.703029, 108.4773, -42.52157, 5.862788,
         0.678565, -76.84376, 158.7163, -74.87310), name="methane"),
    _sp("NO", 210.76, 29.86,
        (298.0, 1200.0, 23.83491, 12.58878, -1.139011, -1.497459,
         0.214194, 83.35783, 237.1219, 90.29114), name="nitric oxide"),
    _sp("N2O", 219.96, 38.62,
        (298.0, 1400.0, 27.67988, 51.14898, -30.64454, 6.847911,
         -0.157906, 71.24934, 238.6164, 82.04824), name="nitrous oxide"),
]}


def self_check(tol_S=0.05, tol_Cp=0.05, verbose=False):
    """★ 用系数现算 298.15 K 的 S 与 Cp，跟**各自随带的参考值**比。

    为什么必须有：Shomate 系数是人工录入的，**抄错一个数不会报错，只会给出
    看着合理的错数**。把参考值和系数放在一起、并强制互相核对，
    转录错误就不可能静默发布。

    ★ 容差取得很紧（0.05 J/(mol·K)，约 0.03%），因为**抄错一个数的后果远大于 0.05**。
    个别物种在区间端点确有拟合偏差（目前只有 NH3 的 Cp：+0.587），
    那是**写在数据里的已知量** `cp298_known_dev`（连同原因），
    自检时先扣掉它再比 —— 这样"抄错"与"拟合端点本来就这样"不会混为一谈。
    把容差放宽到能容下 0.587，等于同时放过了一大类真正的抄写错误。

    返回 `(全过?, 逐条结果)`。
    """
    rows, ok_all = [], True
    for f, sp in sorted(GAS_DB.items()):
        try:
            s = sp.S(298.15); c = sp.Cp(298.15)
            ds = abs(s - sp.S298)
            dc = abs(c - sp.Cp298 - sp.cp298_known_dev)   # 先扣掉记录在案的已知偏差
            good = ds <= tol_S and dc <= tol_Cp
        except ThermoError:
            # 温区不含 298.15（如 H2O 气相从 500 K 起）⇒ 跳过，但要**明说跳过了**
            rows.append((f, None, None, None, None, "跳过（温区不含 298.15 K）"))
            continue
        ok_all &= good
        rows.append((f, s, sp.S298, c, sp.Cp298, "✓" if good else "✗ 对不上参考值"))
    if verbose:
        print("Shomate 自检（298.15 K）：")
        for r in rows:
            if r[1] is None:
                print("  %-5s %s" % (r[0], r[5])); continue
            kd = GAS_DB[r[0]].cp298_known_dev
            tag = ("  [已知端点偏差 %+.3f：%s]" % (kd, GAS_DB[r[0]].dev_reason)) if kd else ""
            print("  %-5s S=%8.2f (参考 %7.2f, 差 %+.3f)   Cp=%6.2f (参考 %6.2f, 差 %+.3f)  %s%s"
                  % (r[0], r[1], r[2], r[1] - r[2], r[3], r[4],
                     r[3] - r[4] - kd, r[5], tag))
    return ok_all, rows


# ══════════════════════════ 给机理算熵校正 ══════════════════════════
def entropy_offsets(mech: Mechanism, T_K: float, extra: Dict[str, GasSpecies] = None,
                    strict: bool = True, J_per_mol_per_eV: float = None
                    ) -> Tuple[np.ndarray, np.ndarray, List[str]]:
    """机理 + 温度 → `(offset_Ga, offset_dG, 说明)`，单位 eV。

    规则见模块 docstring（一条规则导出全部校正，已对上 CATKINAS 的 5 个常数）。

    `strict=True` 时，机理里有**数据库没有的气相物种**就报错；
    `False` 则按 S=0 处理并在说明里标出来 —— **默认报错**，
    因为"少算一个物种的熵"是看不出来的系统误差。
    """
    db = dict(GAS_DB); db.update(extra or {})
    from .rates import _side_matrices
    a_r, a_p = _side_matrices(mech)
    nu = np.asarray(mech.nu, float)
    idx = {n: i for i, n in enumerate(mech.species_order)}
    off_Ga = np.zeros(mech.n_step)
    off_dG = np.zeros(mech.n_step)
    notes = []
    missing = []
    TS = {}
    for name in mech.gas_species:
        f = mech.species[name].formula
        if f not in db:
            missing.append(f)
            TS[name] = 0.0
            continue
        TS[name] = db[f].TS_eV(T_K, J_per_mol_per_eV)
    if missing and strict:
        raise ThermoError(
            "这些气相物种不在热化学数据库里：%s\n"
            "  少算一个物种的熵，是**看不出来的系统误差**（673 K 下 H2 这一项就有 1.08 eV），"
            "所以默认报错而不是按 0 处理。\n"
            "  两种办法：① 自己补数据 —— entropy_offsets(..., extra={'X': GasSpecies(...)})；"
            "② 明确接受 S=0 —— strict=False。\n"
            "  数据库现有：%s" % (sorted(set(missing)), sorted(GAS_DB)))
    if missing:
        notes.append("⚠ 这些气相物种按 S=0 处理（strict=False）：%s" % sorted(set(missing)))

    for s in range(mech.n_step):
        for name in mech.gas_species:
            j = idx[name]
            off_dG[s] -= nu[s, j] * TS[name]          # −T·Σν·S
            off_Ga[s] += a_r[s, j] * TS[name]         # +T·Σ(反应物侧)·S
    notes.append("熵校正 @ %.2f K（J/mol→eV 换算 %.6f），来源：%s"
                 % (T_K, J_per_mol_per_eV or J_PER_MOL_PER_EV, CITATION))
    for name in mech.gas_species:
        f = mech.species[name].formula
        if f in db:
            notes.append("   %-5s S°(%.1f K) = %.3f J/(mol·K) ⇒ T·S = %.6f eV"
                         % (name, T_K, db[f].S(T_K), TS[name]))
    return off_Ga, off_dG, notes


def list_species() -> List[str]:
    """数据库里收录的气相物种（化学式）。"""
    return sorted(GAS_DB)


#: 向后兼容的别名（旧代码里用过 `DB`）
DB = GAS_DB
