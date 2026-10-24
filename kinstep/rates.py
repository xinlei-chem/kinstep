# -*- coding: utf-8 -*-
# KinSTEP - a modality-agnostic microkinetic solver.
# Copyright (C) 2026 The KinSTEP authors.
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.  It is distributed WITHOUT ANY WARRANTY; see the GNU General Public
# License (LICENSE) for details.
r"""[English] Rate constants and the factorised rate law, with the free consistency checks
(Wegscheider cycles, Campbell sum rule).

(The commentary below is in Chinese: it records why each choice was made
and which measurement forced it.)

速率常数与速率律 —— 内核的第二层。

本模块的职责边界
----------------
输入：机理（mech.py）+ 每步的能量/常数（带来源标签）+ 工况 + driver 栈
输出：`kf_step` / `kr_step` / `Keq_step`，以及**因式形式**的速率

    ┌ 输入层：Val（带来源，慢但只过一次）──→ compile() ──→ ndarray（快，进内层循环）

三条取常数的路（**按步选**，不是全局选）
----------------------------------------
真实机理里三条路是**混用**的，所以 `mode` 绑在步上：

| mode        | 公式 | 用在哪 | 真值来源 |
|-------------|------|--------|----------|
| `"tst"`     | kf = A·exp(−Ga_fwd/k_BT)，A 默认 k_BT/h | 表面基元步（DFT 给 Ga） | 过渡态理论 |
| `"collision"`| kf = S0·A_site·P0/√(2π m k_BT) | **吸附步**（无过渡态） | 碰撞理论，CATKINAS 同 |
| `"given"`   | kf = kf0（实验/拟合值） | Koper 这类现象学模型 | Ager notebook |

★★★ 逆向常数**永不独立给**：`kr = kf / Keq`
------------------------------------------------
    Keq = exp(−dG_rxn / k_BT)

这是**热力学自洽**的唯一保证。独立给 kf 和 kr 的模型（很多文献里的表格就是）
会让"稳态解"偷偷违反细致平衡 —— 症状不是报错，而是**平衡位置算错**而曲线照样漂亮。

唯一例外：`mode="given"` 允许同时给 `kf0` 与 `kr0`（Koper 就是这么给的），
此时**不信任它们各自的值，而是反推** Keq0 = kf0/kr0 → dG_rxn = −k_BT·ln Keq0，
并把这个隐含 dG **打印出来让用户核对**。反推之后仍然只有一个自由度。

★ 外场修正怎么进来（与 driver 层的契约）
----------------------------------------
driver 只会改能量，**三条路共用同一个修正量**：

    ΔGa = Ga_corr − Ga_0        （活化自由能被外场改了多少）
    ΔΔG = dG_corr − dG_0        （反应自由能被外场改了多少）

    kf = kf_0 · exp(−ΔGa / k_BT)                （对 tst 路，kf_0 = A·exp(−Ga_0/kT)）
    Keq = Keq_0 · exp(−ΔΔG / k_BT)
    kr = kf / Keq

逐字验证（Koper/Ager，ν_e=+1）：ΔGa = −β·U、ΔΔG = −U ⇒

    kf = kf0·exp(+β·U/k_BT)                       ← notebook 的 k1  ✅
    kr = kr0·exp((ΔΔG−ΔGa)/kT) = kr0·exp(−(1−β)U/k_BT)  ← notebook 的 kminus1  ✅

**这样电化学修正只写在 electro driver 一处**，三条取常数的路自动都对。

★★ 因式形式速率（这是 CATKINAS 最值得学的一招）
------------------------------------------------
    r_fwd = kf·∏(反应物活度^|ν|)
    zrev  = r_rev / r_fwd                  ← 可逆度，∈[0,1]
    r_net = r_fwd·(1 − zrev)

**绝不写 `kf·A − kr·B`**：准平衡步上两项各 10¹⁰ 而差 10⁻⁵，相减把有效数字吃光。
本模块用 `ln r` 做差再取 exp 得 zrev，连 r_fwd 本身溢出都绕开了。

★ Wegscheider 闸（免费自检，别处没有）
--------------------------------------
若机理里存在**净变化为零的循环**（σ·ν = 0 的步组合），热力学要求 σ·dG_rxn = 0。
能量填错一个数，这条就会炸 —— 而其它任何检查都看不出来。
`cycle_consistency()` 用 ν 的左零空间把这件事做掉。
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from .driver import DriverStack
from .i18n import t as _t
from .mech import Mechanism
from .provenance import (ENERGY_KEYS, Val, UnverifiedDataError, as_val, gate,
                         report_defaults)
from .units import EV, H_PLANCK, K_B, K_B_EV, N_A, P0_PA

#: 原子质量单位 → kg
AMU_KG = 1.66053906660e-27

#: 取常数的三条路
RATE_MODES = ("tst", "collision", "given")


class RateModelError(ValueError):
    """速率常数的输入有问题。"""


# ══════════════════════════════════ 输入层 ══════════════════════════════════
@dataclass
class StepEnergy:
    """**一个基元步**的能量/常数输入。用户只填自己这条路需要的几项。

    `mode="tst"`        需要 Ga_fwd + dG_rxn
    `mode="collision"`  需要 dG_rxn + mass_amu（+ 可选 sticking / A_site_m2）
    `mode="given"`      需要 kf0（+ 可选 kr0；给了 kr0 就反推 dG_rxn）
    """
    mode: str = "tst"
    Ga_fwd: Val = field(default_factory=Val.required)
    dG_rxn: Val = field(default_factory=Val.required)
    #: tst 路的前因子；不给就用 k_BT/h
    A_pre: Optional[Val] = None
    #: given 路
    kf0: Optional[Val] = None
    kr0: Optional[Val] = None
    #: collision 路
    mass_amu: Optional[Val] = None
    sticking: Optional[Val] = None
    A_site_m2: Optional[Val] = None

    def __post_init__(self):
        if self.mode not in RATE_MODES:
            raise RateModelError("mode 只能是 %s，收到 %r" % (list(RATE_MODES), self.mode))
        for k in ("Ga_fwd", "dG_rxn", "A_pre", "kf0", "kr0",
                  "mass_amu", "sticking", "A_site_m2"):
            v = getattr(self, k)
            if v is not None and not isinstance(v, Val):
                setattr(self, k, as_val(v))


def _default_sticking():
    return Val(1.0, "default", detail="初始黏附系数，无数据时取 1")


def _default_site_area():
    # 1 个位点约 (3 Å)² ——  金属 (111) 面典型值；只是默认，输出里会标出来
    return Val(9.0e-20, "default", detail="位点面积 9 Å²（(111) 面典型值）")


class Energetics:
    """整套机理的能量输入 + **求解前的来源闸**。

    用法（键可以是 1 起的步号，也可以是机理原文那一行）::

        en = Energetics(mech, {
            1: dict(mode="collision", dG_rxn=-0.47, mass_amu=28.0),
            2: dict(Ga_fwd=1.21, dG_rxn=0.35),          # 默认 tst
        })
    """

    def __init__(self, mech: Mechanism, table: Dict, default_mode="tst"):
        self.mech = mech
        self.steps: List[StepEnergy] = []
        raw2idx = {st.raw: i for i, st in enumerate(mech.steps)}
        norm: Dict[int, Dict] = {}
        for key, val in table.items():
            if isinstance(key, (int, np.integer)):
                idx = int(key) - 1
                if not 0 <= idx < mech.n_step:
                    raise RateModelError(
                        "步号 %s 超出范围：机理只有 %d 步（步号从 1 数）" % (key, mech.n_step))
            elif key in raw2idx:
                idx = raw2idx[key]
            else:
                raise RateModelError(
                    "认不出键 %r。键要么是 1 起的步号，要么是机理里那一行的原文。\n"
                    "机理各行：\n  %s"
                    % (key, "\n  ".join("(%d) %s" % (i + 1, s.raw)
                                        for i, s in enumerate(mech.steps))))
            if idx in norm:
                raise RateModelError("第 %d 步被指定了两次" % (idx + 1))
            norm[idx] = val
        missing_steps = [i + 1 for i in range(mech.n_step) if i not in norm]
        if missing_steps:
            raise RateModelError(
                "这些步没给能量/常数：%s\n"
                "每一步都必须给 —— **不替用户猜**（猜一个数就能让火山图看着很合理）。"
                % missing_steps)
        for i in range(mech.n_step):
            spec = norm[i]
            if isinstance(spec, StepEnergy):
                se = spec
            else:
                spec = dict(spec)
                spec.setdefault("mode", default_mode)
                se = StepEnergy(**spec)
            self._fill_defaults(se, i)
            self.steps.append(se)

    def _fill_defaults(self, se: StepEnergy, idx: int):
        """补默认值 —— 但**每个默认值都带 'default' 来源，输出里必须可见**。"""
        if se.mode == "collision":
            if se.sticking is None:
                se.sticking = _default_sticking()
            if se.A_site_m2 is None:
                se.A_site_m2 = _default_site_area()
            if se.mass_amu is None or se.mass_amu.is_missing:
                se.mass_amu = Val.required(
                    _t("第 %d 步走碰撞理论，需要吸附分子的分子量 mass_amu（如 CO 给 28.0）"
                       % (idx + 1),
                       "step %d uses collision theory and needs the molar mass "
                       "mass_amu of the adsorbing molecule (e.g. 28.0 for CO)"
                       % (idx + 1)))
        if se.mode == "given":
            if se.kf0 is None:
                se.kf0 = Val.required(
                    _t("第 %d 步 mode='given'，需要 kf0（正向速率常数，1/s）" % (idx + 1),
                       "step %d has mode='given' and needs kf0 "
                       "(forward rate constant, 1/s)" % (idx + 1)))
            # 给了 kr0 ⇒ 反推 dG_rxn（见模块 docstring）；没给就必须给 dG_rxn
            if se.kr0 is not None and not se.kr0.is_missing:
                se.dG_rxn = Val(None, "REQUIRED")   # 占位，compile 时由 kf0/kr0 反推
        if se.mode == "tst" and (se.Ga_fwd is None or se.Ga_fwd.is_missing):
            se.Ga_fwd = Val.required(
                _t("第 %d 步走过渡态理论，需要正向活化自由能 Ga_fwd (eV)" % (idx + 1),
                   "step %d uses transition-state theory and needs the forward "
                   "activation free energy Ga_fwd (eV)" % (idx + 1)))
        # ★ dG_rxn 原本只是 `Val.required()` 的**空提示**，于是来源闸会打印「没有提示」。
        #   而 `kinstep new` 的输出承诺「每个能量槽位都标了要算什么」—— 自相矛盾，
        #   而且恰好落在审稿人跑的第一条命令上（HER 模板 6 个缺口里有 3 个是它）。
        #   `given` 模式不在此列：那边给了 kr0 就反推 dG_rxn，已有自己的占位逻辑。
        if se.mode in ("tst", "collision") and (se.dG_rxn is None
                                                or se.dG_rxn.is_missing):
            se.dG_rxn = Val.required(
                _t("第 %d 步需要反应自由能 dG_rxn (eV) = G(产物) − G(反应物)，"
                   "它定逆向速率常数 kr = kf/Keq" % (idx + 1),
                   "step %d needs the reaction free energy dG_rxn (eV) = "
                   "G(products) - G(reactants); it fixes the reverse rate "
                   "constant kr = kf/Keq" % (idx + 1)))

    # ── 求解前的闸 ──
    def audit_items(self):
        """给 provenance.gate() 的清单。"""
        items = []
        for i, se in enumerate(self.steps):
            label = _t("第%d步" % (i + 1), "step %d" % (i + 1))
            keys = {"tst": ("Ga_fwd", "dG_rxn", "A_pre"),
                    "collision": ("dG_rxn", "mass_amu", "sticking", "A_site_m2"),
                    "given": ("kf0", "kr0")}[se.mode]
            for k in keys:
                v = getattr(se, k)
                if v is None:
                    continue
                if se.mode == "given" and k == "kr0":
                    pass
                elif se.mode == "given" and k == "kf0" and v.is_missing:
                    pass
                if v.is_missing and se.mode == "given" and k == "kr0":
                    continue      # kr0 本来就是可选的
                items.append((label, k, v))
        return items

    def check(self, strict_energy=True):
        """闸：缺值或能量来源不可信 ⇒ 抛异常。返回"用了哪些默认值"的清单。"""
        return gate(self.audit_items(), strict_energy=strict_energy)


# ══════════════════════════════ 基准常数（U=0 / 无外场）══════════════════════
@dataclass
class RateConstants:
    """编译好的速率常数。**全是 ndarray，可以直接进内层循环。**"""
    kf: np.ndarray          # 正向，1/s（吸附步已把无量纲压力的因子算进去）
    kr: np.ndarray          # 逆向，1/s
    Keq: np.ndarray         # 平衡常数
    Ga_used: np.ndarray     # 实际用的 Ga（含外场修正），eV；given 路为 nan
    dG_used: np.ndarray     # 实际用的 dG（含外场修正），eV
    T_K: float
    #: 逐步的**平均场近邻因子** `n_neighbor`（默认全 1）。见 `rate_constants` 的说明
    prefactor: np.ndarray = None
    notes: List[str] = field(default_factory=list)

    def __post_init__(self):
        if self.prefactor is None:
            self.prefactor = np.ones_like(self.kf)

    def table(self) -> str:
        L = ["  步  %-12s %-12s %-12s %-10s %-10s" % ("kf", "kr", "Keq", "Ga(eV)", "dG(eV)")]
        for i in range(len(self.kf)):
            ga = "  —  " if not np.isfinite(self.Ga_used[i]) else "%+.4f" % self.Ga_used[i]
            L.append("  %2d  %-12.4e %-12.4e %-12.4e %-10s %+.4f"
                     % (i + 1, self.kf[i], self.kr[i], self.Keq[i], ga, self.dG_used[i]))
        return "\n".join(L)


def k_tst(Ga_eV, T_K, A_pre=None):
    """过渡态理论：kf = A·exp(−Ga/k_BT)，A 默认 k_BT/h。"""
    A = (K_B * T_K / H_PLANCK) if A_pre is None else A_pre
    return A * np.exp(-np.asarray(Ga_eV, float) / (K_B_EV * T_K))


def k_collision(mass_amu, T_K, sticking=1.0, A_site_m2=9.0e-20, P_ref_Pa=P0_PA):
    """碰撞理论吸附常数（配**无量纲压力** p = P/P0 使用），单位 1/s。

    k_ads = S0 · A_site · P0 / √(2π m k_BT)

    量纲核对：Pa·m²/√(kg·J) = N/(kg·m/s) = (kg·m/s²)/(kg·m/s) = 1/s ✓
    """
    m = float(mass_amu) * AMU_KG
    return (float(sticking) * float(A_site_m2) * P_ref_Pa
            / math.sqrt(2.0 * math.pi * m * K_B * T_K))


def _base_constants(en: Energetics, T_K: float):
    """算 U=0 / 无外场下的 (kf_0, Keq_0, Ga_0, dG_0)，以及提示信息。"""
    n = en.mech.n_step
    kf0 = np.zeros(n)
    dG0 = np.zeros(n)
    Ga0 = np.full(n, np.nan)
    notes: List[str] = []
    kT = K_B_EV * T_K
    for i, se in enumerate(en.steps):
        if se.mode == "given":
            kf0[i] = float(se.kf0)
            if se.kr0 is not None and not se.kr0.is_missing:
                # ★ 反推而不是直接用两个独立值（见模块 docstring）
                Keq_implied = float(se.kf0) / float(se.kr0)
                dG0[i] = -kT * math.log(Keq_implied)
                notes.append(
                    "第%d步 由 kf0=%.4g / kr0=%.4g 反推：Keq=%.4g ⇒ dG_rxn=%+.4f eV"
                    "（请核对这个隐含的热力学量是否合理）"
                    % (i + 1, float(se.kf0), float(se.kr0), Keq_implied, dG0[i]))
            else:
                dG0[i] = float(se.dG_rxn)
        elif se.mode == "collision":
            kf0[i] = k_collision(float(se.mass_amu), T_K,
                                 float(se.sticking), float(se.A_site_m2))
            dG0[i] = float(se.dG_rxn)
            notes.append("第%d步 碰撞理论：kf=%.4g s⁻¹（S0=%g, A_site=%.3g m², M=%g）"
                         % (i + 1, kf0[i], float(se.sticking),
                            float(se.A_site_m2), float(se.mass_amu)))
        else:                                    # tst
            Ga0[i] = float(se.Ga_fwd)
            A = None if se.A_pre is None or se.A_pre.is_missing else float(se.A_pre)
            kf0[i] = k_tst(Ga0[i], T_K, A)
            dG0[i] = float(se.dG_rxn)
            # ★ 物理合理性：活化能不能低于反应能（否则逆向无势垒，Ga_rev<0）
            if Ga0[i] < dG0[i] - 1e-9:
                notes.append("⚠ 第%d步 Ga_fwd=%.3f < dG_rxn=%.3f ⇒ 逆向势垒为负 %.3f eV，"
                             "这在物理上说不通，请核对能量" % (i + 1, Ga0[i], dG0[i],
                                                        Ga0[i] - dG0[i]))
    return kf0, dG0, Ga0, notes


#: mode 的数值编码（编译后内层只认整数，不认字符串）
_MODE_CODE = {"tst": 0, "collision": 1, "given": 2}


@dataclass
class CompiledRates:
    """★ 编译后的**数值层** —— 只有 ndarray，`Val` 与来源闸全留在编译之前。

    为什么必须有这一层（2026-10-22 实测教训）
    ------------------------------------------
    初版 `solve.py` 在 ODE 右端**每次都调 `rate_constants()`**，而它里面跑的是
    完整的来源闸（`gate()` → `audit()` → 逐个 `Val.__float__`）。
    症状：Koper 的 16 秒 CV **跑了 2 分 19 秒还没出**，CPU 满载。
    台账 T2.5 条目早就写明"Val 不许进内层循环"（实测慢 5 倍）—— 我在 solve.py 里
    违反了自己的设计。⇒ 现在分成两层：

        compile_rates()  ← 过闸、取值、转 ndarray，**只做一次**
        apply_state()    ← 只做 driver 修正 + exp，内层循环走这条

    温度也存成数组而非在编译时写死 ⇒ 变温协议（TPD/升温）照样走快路。
    """
    mech: Mechanism
    mode: np.ndarray            # 0=tst 1=collision 2=given
    Ga_in: np.ndarray           # tst 的 Ga_fwd（其余为 nan）
    dG_in: np.ndarray           # 直接给的 dG_rxn（给了 kr0 的步为 nan）
    Keq_given: np.ndarray       # 由 kf0/kr0 反推的 Keq（其余为 nan）
    A_pre: np.ndarray           # tst 前因子（nan = 用 k_BT/h）
    kf_given: np.ndarray        # given 的 kf0
    mass_amu: np.ndarray
    sticking: np.ndarray
    A_site: np.ndarray
    prefactor: np.ndarray       # n_neighbor
    reversible: np.ndarray
    a_r: np.ndarray             # 反应物侧指数矩阵（预算，别进内层）
    a_p: np.ndarray
    notes: List[str] = field(default_factory=list)

    def base_at_T(self, T_K):
        """算 U=0 / 无外场下的 (kf_0, dG_0, Ga_0)。**纯 numpy，无 Python 对象。**"""
        kT = K_B_EV * T_K
        kf = np.where(self.mode == 2, self.kf_given, 0.0)
        # tst
        m0 = self.mode == 0
        if m0.any():
            A = np.where(np.isnan(self.A_pre), K_B * T_K / H_PLANCK, self.A_pre)
            kf = np.where(m0, A * np.exp(-np.nan_to_num(self.Ga_in, nan=0.0) / kT), kf)
        # collision
        m1 = self.mode == 1
        if m1.any():
            m = self.mass_amu * AMU_KG
            kc = (self.sticking * self.A_site * P0_PA
                  / np.sqrt(2.0 * np.pi * np.where(m > 0, m, 1.0) * K_B * T_K))
            kf = np.where(m1, kc, kf)
        # dG：直接给的用给的；由 kf0/kr0 反推的随温度走
        dG = np.where(np.isnan(self.Keq_given), np.nan_to_num(self.dG_in, nan=0.0),
                      -kT * np.log(np.where(np.isnan(self.Keq_given), 1.0, self.Keq_given)))
        return kf, dG, self.Ga_in.copy()


def compile_rates(en: Energetics, strict_energy: bool = True) -> CompiledRates:
    """过来源闸、把 `Val` 取成 ndarray。**整个求解过程只调一次。**"""
    defaulted = en.check(strict_energy=strict_energy)
    n = en.mech.n_step
    z = lambda: np.full(n, np.nan)
    mode = np.zeros(n, int)
    Ga_in, dG_in, Keq_g, A_pre, kf_g = z(), z(), z(), z(), z()
    mass, stick, area = z(), z(), z()
    notes: List[str] = []
    for i, se in enumerate(en.steps):
        mode[i] = _MODE_CODE[se.mode]
        if se.mode == "given":
            kf_g[i] = float(se.kf0)
            if se.kr0 is not None and not se.kr0.is_missing:
                Keq_g[i] = float(se.kf0) / float(se.kr0)
                notes.append(
                    "第%d步 由 kf0=%.4g / kr0=%.4g 反推 Keq=%.4g —— ΔG 随温度走，"
                    "各温度下的值见 RateConstants.dG_used（请核对是否合理）"
                    % (i + 1, float(se.kf0), float(se.kr0), Keq_g[i]))
            else:
                dG_in[i] = float(se.dG_rxn)
        elif se.mode == "collision":
            mass[i], stick[i], area[i] = (float(se.mass_amu), float(se.sticking),
                                          float(se.A_site_m2))
            dG_in[i] = float(se.dG_rxn)
            notes.append("第%d步 碰撞理论：S0=%g, A_site=%.3g m², M=%g"
                         % (i + 1, stick[i], area[i], mass[i]))
        else:
            Ga_in[i] = float(se.Ga_fwd)
            dG_in[i] = float(se.dG_rxn)
            if se.A_pre is not None and not se.A_pre.is_missing:
                A_pre[i] = float(se.A_pre)
            if Ga_in[i] < dG_in[i] - 1e-9:
                notes.append("⚠ 第%d步 Ga_fwd=%.3f < dG_rxn=%.3f ⇒ 逆向势垒为负 %.3f eV，"
                             "这在物理上说不通，请核对能量"
                             % (i + 1, Ga_in[i], dG_in[i], Ga_in[i] - dG_in[i]))
    pref = np.array([float(st.params.get("n_neighbor", 1.0)) for st in en.mech.steps])
    if np.any(pref != 1.0):
        notes.append("近邻因子 n_neighbor = %s（只乘净速率，不改平衡位置）"
                     % np.array2string(pref, precision=3))
    if defaulted:
        notes.append(report_defaults(defaulted))
    a_r, a_p = _side_matrices(en.mech)
    return CompiledRates(mech=en.mech, mode=mode, Ga_in=Ga_in, dG_in=dG_in,
                         Keq_given=Keq_g, A_pre=A_pre, kf_given=kf_g,
                         mass_amu=mass, sticking=stick, A_site=area,
                         prefactor=pref, reversible=en_reversible(en.mech),
                         a_r=a_r, a_p=a_p, notes=notes)


def apply_state(cr: CompiledRates, T_K: float, driver: Optional[DriverStack] = None,
                state: Optional[Dict] = None, theta=None) -> RateConstants:
    """★ 内层循环走这条：编译结果 + 工况 → 速率常数。**不碰 Val、不过闸。**"""
    kf_0, dG_0, Ga_0 = cr.base_at_T(T_K)
    kT = K_B_EV * T_K
    if driver is not None:
        Ga_in = np.where(np.isfinite(Ga_0), Ga_0, 0.0)
        Ga_c_in, dG_c = driver.modify_energies(cr.mech, Ga_in, dG_0, state, theta)
        dGa = np.asarray(Ga_c_in, float) - Ga_in
        Ga_c = Ga_0 + dGa
    else:
        dGa, dG_c, Ga_c = np.zeros_like(dG_0), dG_0, Ga_0
    kf = kf_0 * np.exp(-dGa / kT)
    Keq = np.exp(-dG_c / kT)
    with np.errstate(divide="ignore", over="ignore", invalid="ignore"):
        kr = np.where(cr.reversible, kf / Keq, 0.0)
    return RateConstants(kf=kf, kr=kr, Keq=Keq, Ga_used=Ga_c, dG_used=dG_c,
                         T_K=float(T_K), prefactor=cr.prefactor, notes=[])


def rate_constants(en: Energetics, T_K: float, driver: Optional[DriverStack] = None,
                   state: Optional[Dict] = None, theta=None,
                   strict_energy: bool = True) -> RateConstants:
    """**便捷入口**（编译 + 施加工况各一次）。

    ⚠ **不要在循环里调它** —— 它每次都会跑一遍来源闸。
    循环里请用 `compile_rates()` 一次 + `apply_state()` 多次（见 CompiledRates 的说明）。
    """
    defaulted = en.check(strict_energy=strict_energy)
    kf_0, dG_0, Ga_0, notes = _base_constants(en, T_K)
    kT = K_B_EV * T_K

    Ga_c, dG_c = Ga_0.copy(), dG_0.copy()
    if driver is not None:
        # ★ Ga_0 对 given/collision 路是 nan；driver 做的是**加减**，nan 会穿透，
        #   所以先用 0 占位送进去，再把修正量取回来 —— 修正量本身是有限的。
        Ga_in = np.where(np.isfinite(Ga_0), Ga_0, 0.0)
        Ga_c_in, dG_c = driver.modify_energies(en.mech, Ga_in, dG_0, state, theta)
        dGa = np.asarray(Ga_c_in, float) - Ga_in          # ΔGa
        Ga_c = Ga_0 + dGa                                  # nan 的仍是 nan，符合语义
    else:
        dGa = np.zeros_like(dG_0)

    ddG = dG_c - dG_0                                      # ΔΔG
    kf = kf_0 * np.exp(-dGa / kT)
    Keq = np.exp(-dG_c / kT)
    with np.errstate(divide="ignore", over="ignore"):
        kr = np.where(en_reversible(en.mech), kf / Keq, 0.0)

    # ★ 平均场近邻因子 n_neighbor（可选，逐步给：`{n_neighbor: 4.0}`）
    #
    # 双吸附态反应 A* + B* → … 要求两者**相邻**；平均场下每个位点的 A–B 近邻对数
    # = Z·θ_A·θ_B，所以速率要乘一个配位数 Z。Koper(1998) 的 CO stripping 就是
    # `v2 = Z·k2·θ_OH·θ_CO`，Z = 4。
    #
    # ⚠ 它**只乘在 r_fwd 上，不改 zrev** —— 因为本包的 zrev 由 ΔG 与活度直接算
    # （见 rates_factored），与 kf/kr 无关。等价于把 kf 与 kr 同乘 Z，
    # 所以**平衡位置不变** —— 这是唯一不破坏细致平衡的放法：
    # 组合计数因子是"位点配对统计"，不是热力学量，不许动 Keq。
    prefactor = np.array([float(st.params.get("n_neighbor", 1.0))
                          for st in en.mech.steps])
    if np.any(prefactor != 1.0):
        notes.append("近邻因子 n_neighbor = %s（只乘净速率，不改平衡位置）"
                     % np.array2string(prefactor, precision=3))
    if defaulted:
        notes.append(report_defaults(defaulted))
    return RateConstants(kf=kf, kr=kr, Keq=Keq, Ga_used=Ga_c, dG_used=dG_c,
                         T_K=float(T_K), prefactor=prefactor, notes=notes)


def en_reversible(mech: Mechanism) -> np.ndarray:
    """哪些步是可逆的（`->` 写的步，逆向常数置零）。"""
    return np.array([st.reversible for st in mech.steps], bool)


# ═══════════════════════════ 活度与因式形式速率 ═══════════════════════════
def activities(mech: Mechanism, theta: Dict[str, float] = None,
               p_gas: Dict[str, float] = None, c_liq: Dict[str, float] = None,
               theta_vac: Dict[str, float] = None) -> np.ndarray:
    """按 `mech.species_order` 组装活度向量。

    表面物种用覆盖度 θ、气相用无量纲分压 p/P0、液相用无量纲浓度 c/C0、
    空位用 θ_vac（不给就按"1 − 该位点上所有吸附物种"自动算）。
    """
    theta, p_gas, c_liq = theta or {}, p_gas or {}, c_liq or {}
    theta_vac = dict(theta_vac or {})
    # 空位自动补：每种位点 θ_vac = 1 − Σθ
    for site in mech.sites:
        vac = mech.vacancy_of(site)
        if vac in theta_vac or vac not in mech.species:
            continue
        used = sum(theta.get(n, 0.0) for n in mech.surface_species
                   if mech.species[n].site == site)
        theta_vac[vac] = 1.0 - used
    act = np.zeros(len(mech.species_order))
    for i, name in enumerate(mech.species_order):
        sp = mech.species[name]
        if sp.phase == "surf":
            act[i] = theta.get(name, 0.0)
        elif sp.phase == "site":
            act[i] = theta_vac.get(name, 0.0)
        elif sp.phase == "gas":
            act[i] = p_gas.get(name, p_gas.get(sp.formula, 0.0))
        elif sp.phase == "liq":
            act[i] = c_liq.get(name, c_liq.get(sp.formula, 0.0))
        else:                                   # electron：不进质量作用定律
            act[i] = 1.0
    return act


def _side_matrices(mech: Mechanism):
    """把 ν 拆成反应物侧指数 `a_r` 与产物侧指数 `a_p`（都取正）。

    ★ 为什么不直接用 ν 的正负：同一物种可能**两边都出现**（如 `A + * <-> A* + *` 这种
      写法里的 `*`），净计量会把它抵消掉，但质量作用定律里两边各自都要算。
    """
    idx = {n: i for i, n in enumerate(mech.species_order)}
    ns, nsp = mech.n_step, len(mech.species_order)
    a_r = np.zeros((ns, nsp))
    a_p = np.zeros((ns, nsp))
    e_idx = idx.get("e-")
    for s, st in enumerate(mech.steps):
        for n, c in st.reactants.items():
            a_r[s, idx[n]] += c
        for n, c in st.products.items():
            a_p[s, idx[n]] += c
    if e_idx is not None:                 # 电子不进质量作用定律（它的作用在 U 里）
        a_r[:, e_idx] = 0.0
        a_p[:, e_idx] = 0.0
    return a_r, a_p


def _weighted_ln(expo: np.ndarray, ln_a: np.ndarray) -> np.ndarray:
    """Σ_j expo[s,j]·ln a[j]，但**指数为 0 的项一律不参与**。

    ★ 这一行是一个真 bug 的修法（2026-10-22 测试当场抓到）：
      原先写 `expo @ ln_a`。某物种活度为 0 ⇒ ln_a = −inf；若它在这一步的指数是 0，
      numpy 算出 `0 × (−inf) = nan`，**整个 ln r 被污染成 nan**，于是速率被判成 0。
      实测症状：Koper 机理里 `H+` 只出现在产物侧，却把步 1 的 `r_fwd` 清零了 ——
      而且不报错不警告，看上去就是"这一步不反应"。
    """
    masked = np.where(expo > 0.0, ln_a[None, :], 0.0)
    return (expo * masked).sum(axis=1)


def rates_factored(mech: Mechanism, rc: RateConstants, act: np.ndarray,
                   a_r=None, a_p=None):
    """★ 因式形式速率：返回 `(rate_fwd, zrev, rate_net)`。

    zrev = r_rev/r_fwd **不用两个速率相除**，而是直接从热力学量算：

        ln z = ΔG_rxn/k_BT + Σ(a_p − a_r)·ln a          （因为 kr/kf = 1/Keq）

    好处有两个：① 项数少一半、条件数更好；② **完全不碰 kr** ——
    准平衡吸附步的 kr 常常是 10³⁰ 这种会溢出的量级。

    ⚠ 精度的真实边界（必须写清楚，别自欺）
    --------------------------------------
    `zrev ∈ [0, ∞)`，**不是 [0,1]**；z>1 只是说这一步净向逆方向走。

    而当 z → 1（准平衡步）时，**本函数给出的 (1−z) 相对误差会很大**：
    z 本身只能算到 ~1e-16 的相对精度，所以 1−z ≈ 1e-14 时有效数字只剩两位。
    实测：解析真值 6.2500e-6，本函数给 6.6613e-6（差 6.6%）。
    换成 ln 做差、换成高精度常数都**救不了** —— 信息本身就不在 (kf, kr, a) 里。

    ⇒ **这不是实现问题，是变量选择问题。** CATKINAS 的"以可逆度为未知数"
      不是一个求和技巧，而是**换变量**：稳态求解时 z 是牛顿法的自变量，
      (1−z) 是一个精确的 double，从来不由相除得到。
      所以：本函数用于**瞬态/ODE（给定 θ 求速率）**；
      稳态路径必须走 `rates_from_z()`，由求解器持有 z。
    """
    if a_r is None or a_p is None:
        a_r, a_p = _side_matrices(mech)
    act = np.asarray(act, float)
    kT = K_B_EV * rc.T_K
    with np.errstate(divide="ignore", invalid="ignore"):
        ln_a = np.where(act > 0.0, np.log(np.where(act > 0.0, act, 1.0)), -np.inf)
        ln_f = np.log(np.where(rc.kf > 0.0, rc.kf, np.nan)) + _weighted_ln(a_r, ln_a)
        # ln z = ΔG/kT + Σ(a_p−a_r)·ln a —— 产物侧与反应物侧各自掩蔽后再相减，
        #        这样"指数为 0 的物种活度为 0"不会污染任何一侧
        # ★ 两侧都含零活度物种时，`(-inf) − (-inf) = nan` 会抛 RuntimeWarning。
        #   结果本来就会被 isfinite 掩掉，但警告会淹没真正要看的输出 ⇒ 干净地处理：
        #   任一侧为 -inf 即判定 z = 0（该方向走不动）。
        _lp, _lr = _weighted_ln(a_p, ln_a), _weighted_ln(a_r, ln_a)
        _both = np.isfinite(_lp) & np.isfinite(_lr)
        ln_z = np.where(_both, (rc.dG_used / kT) + np.where(_both, _lp, 0.0)
                        - np.where(_both, _lr, 0.0), -np.inf)
    r_fwd = rc.prefactor * np.where(np.isfinite(ln_f), np.exp(ln_f), 0.0)
    ln_z = np.where(rc.kr > 0.0, ln_z, -np.inf)    # 不可逆步（`->` 写的）z ≡ 0
    zrev = np.where(np.isfinite(ln_z), np.exp(np.clip(ln_z, -700.0, 700.0)), 0.0)
    # ★★ (1−z) 必须用 **expm1(ln z)**，不能写 1.0 − z（2026-10-22 实测）：
    #   1−z 真值 1e-12 时，`1 − exp(ln z)` 的相对误差已经是 2.2e-05；
    #   1e-14 时是 8.0e-04。而 `−expm1(ln z)` 在这几档**相对误差全是 0**。
    #   expm1 把瓶颈从"相减"搬到了"ln z 本身的绝对精度"（~1e-15）——
    #   所以准平衡步要真正算准，ln z 还得是**未知数**（见 solve.steady_state）。
    one_minus_z = -np.expm1(np.clip(ln_z, -700.0, 700.0))
    one_minus_z = np.where(np.isfinite(ln_z), one_minus_z, 1.0)
    r_net = r_fwd * one_minus_z
    return r_fwd, zrev, r_net


def rates_from_z(mech: Mechanism, rc: RateConstants, act: np.ndarray, zrev,
                 a_r=None, a_p=None):
    """给定可逆度 z 直接算净速率 —— T6 的稳态求解器以 z 为未知数，要用这个。"""
    if a_r is None or a_p is None:
        a_r, a_p = _side_matrices(mech)
    r_fwd, _, _ = rates_factored(mech, rc, act, a_r, a_p)
    return r_fwd * (1.0 - np.asarray(zrev, float))


# ════════════════════════════ Wegscheider 自检闸 ════════════════════════════
def cycle_consistency(mech: Mechanism, dG_rxn, tol=1e-6):
    """★ 热力学循环自检：净变化为零的步组合，其 ΔG 之和必须为零。

    返回 `(ok, 报告字符串, 违例列表)`。

    原理
    ----
    若存在整数向量 σ 使 σᵀ·ν = 0（走完这组步，所有物种数目回到原处），
    则这是一个**热力学循环**，Wegscheider 条件要求 σᵀ·dG = 0。
    这条约束**与速率无关**，纯热力学，所以是免费的 —— 而且它能抓到
    "某一步能量填错了一个数" 这种其它任何检查都看不出来的错。

    注意：只有存在平行路径的机理才有非平凡循环；线性序列机理没有，
    此时返回 ok=True 并说明"没有可检的循环"（**不是检查通过的假象**）。
    """
    nu = np.asarray(mech.nu, float)
    dG = np.asarray(dG_rxn, float)
    if nu.shape[0] < 2:
        return True, "只有 1 步，无循环可检。", []
    # ν 的左零空间 = null(νᵀ)
    u, s, vt = np.linalg.svd(nu.T, full_matrices=True)
    rank = int((s > max(nu.shape) * np.finfo(float).eps * (s[0] if s.size else 1.0)).sum())
    null_dim = nu.shape[0] - rank
    if null_dim <= 0:
        return True, ("机理的 %d 步线性无关（秩 %d），**没有热力学循环可检** —— "
                      "这不等于能量没错，只是这道闸在这个机理上用不上。"
                      % (nu.shape[0], rank)), []
    null_basis = vt[rank:].T if vt.shape[0] > rank else np.zeros((nu.shape[0], 0))
    bad = []
    for k in range(null_basis.shape[1]):
        sigma = null_basis[:, k]
        sigma = sigma / np.max(np.abs(sigma))
        resid = float(sigma @ dG)
        if abs(resid) > tol:
            bad.append((sigma.copy(), resid))
    if not bad:
        return True, "存在 %d 个热力学循环，ΔG 之和均为 0（残差 < %g）✓" % (null_dim, tol), []
    lines = ["⛔ 热力学循环不自洽（%d 个）—— 一定是某一步的 dG_rxn 填错了：" % len(bad)]
    for sigma, resid in bad:
        terms = " ".join("%+.3g·步%d" % (sigma[i], i + 1)
                         for i in range(len(sigma)) if abs(sigma[i]) > 1e-8)
        lines.append("  循环 [%s] 的 ΣΔG = %+.4g eV（应为 0）" % (terms, resid))
    lines.append("  ⇒ 这条约束与速率无关，纯热力学；**算不出来对的结果，先别往下跑**。")
    return False, "\n".join(lines), bad


def drc_sum_check(drc, tol=1e-3):
    """Campbell 求和规则：Σ DRC = 1（T8 会用；放这里是因为它也属于"免费自检"）。"""
    s = float(np.sum(drc))
    return abs(s - 1.0) <= tol, s
