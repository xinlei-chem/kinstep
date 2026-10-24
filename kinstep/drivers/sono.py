# -*- coding: utf-8 -*-
# KinSTEP - a modality-agnostic microkinetic solver.
# Copyright (C) 2026 The KinSTEP authors.
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.  It is distributed WITHOUT ANY WARRANTY; see the GNU General Public
# License (LICENSE) for details.
r"""[English] Sonocatalysis driver: weighted sub-states for cavitation hot spots, or for any
non-uniform distribution of operating conditions.

(The commentary below is in Chinese: it records why each choice was made
and which measurement forced it.)

超声催化 driver（T13）：空化热点的多工况加权。

理论依据（2026-10-23 实地检索核对，不是凭印象）
------------------------------------------------
**热点理论（hot-spot theory）** 是声化学公认的机制（Suslick 及其后续工作）：
气泡压缩快于热传导 ⇒ 形成微秒级的局部高温高压区。

    多气泡云（multibubble）：~5000 K，~1000 atm，升降温速率 > 10¹⁰ K/s
    —— Suslick, *Phil. Trans. R. Soc. A* 1999；*IEEE Trans.* 1986 等

⚠⚠ **但温度随测法差很多，绝不能写死**：
    单气泡声致发光（SBSL）光谱反推：**5000 – 20000 K**
    多气泡声化学/声致发光反推：     **750 – 6000 K**
⇒ 本模块把热点条件**当输入**，并在超出文献区间时警告，而不是给一个"典型值"。

★ **两点模型（two-site）是 Suslick 自己的框架**，不是我原以为的"三区"
----------------------------------------------------------------
反应位点有两个：**气泡内部的气相**（最热）与**气泡周围的界面液壳**（次热）；
本体液相基本处在环境温度。他们用化学剂量计验证了这个分区，
并发现**溶剂蒸气压与热点温度相关**（蒸气进入气泡会吸热、降低崩塌温度）。
⇒ 本模块支持**任意多个区**（`zones`），两点/三区只是它的特例。

★★ 本模块**不编造**"声强/频率 → 热点条件"的公式
------------------------------------------------
检索结果明确：
  · 声强↑ ⇒ 活性气泡数↑ **且** 单泡产率↑（两个机制叠加）
  · 频率的影响**文献相互矛盾**：有报氧化剂产率随频率**下降**的，
    也有报活性气泡数随频率**上升**的 —— 净效应取决于哪项主导
  · 单泡产率要靠完整气泡动力学（Rayleigh–Plesset + 泡内数十步反应）才算得出来，
    那是**另一个尺度的模型**，不属于表面 MKM

⇒ 编一个拟合式进来，只会给出"看着合理但没有依据"的趋势。
  本模块提供的是 **把热点条件与权重正确地加权进 MKM** 这一件事，
  条件本身由用户从气泡动力学或**实测剂量计**（如 H₂O₂/Weissler 法）给定。
  `calibrate_fraction_from_yield()` 就是把文献那条路（实测产率 + 单泡摩尔产率
  ⇒ 活性气泡数密度）直接做成函数。

数学形式
--------
体系在时间/空间上非均匀 ⇒ **总速率是各工况的加权和**，不是单一工况的速率：

    r_total = Σ_i w_i · r(state_i)  ,   Σ_i w_i = 1

这正是 `Driver.modify_state()` 这个口子的存在理由（T2.5 加的）。
**注意它与 `modify_energies` 的区别**：后者假定"外场只改能量"，
对超声**根本不成立** —— 热点改的是温度本身，而速率对温度是指数敏感的。
"""
import warnings

import numpy as np

from ..driver import Driver, register_driver

#: 文献报道的热点温度区间（用于越界警告，不是默认值）
T_HOTSPOT_MULTIBUBBLE = (750.0, 6000.0)     # 多气泡：声化学/声致发光反推
T_HOTSPOT_SINGLE = (5000.0, 20000.0)        # 单气泡 SBSL 光谱反推

CITATION = (
    "热点理论与条件区间：Suslick 等，Phil. Trans. R. Soc. A 357 (1999) 335；"
    "IEEE Trans. UFFC 33 (1986) 143；Flint & Suslick, Science 253 (1991) 1397。"
    "多气泡 ~5000 K / ~1000 atm / >10^10 K/s；SBSL 反推 5000–20000 K，"
    "多气泡反推 750–6000 K —— **数值随测法差异很大**。"
)


@register_driver
class SonoDriver(Driver):
    driver_name = "sono"
    state_params = (
        "sono_intensity_W_cm2",    # 声强（仅notes/供标定用，**不参与公式**）
        "sono_freq_kHz",           # 频率（同上）
        "sono_duty",               # 占空比 0~1：超声开启的时间份额
        "hotspot_T_K",             # 热点温度
        "hotspot_P_bar",           # 热点压力
        "hotspot_fraction",        # 热点的时空份额
        "interfacial_T_K",         # 界面液壳温度（两点模型的第二个位点）
        "interfacial_fraction",    # 界面液壳的份额
        "sono_zones",              # 任意多区：[{T_K, P_bar, weight}, ...]
        "sono_mass_transfer_gain", # 声流传质增益（改的是外扩散，不是本模块的主效应）
    )
    step_params = ()
    splits_state = True
    is_additive = True
    needs_theta = False
    provides = ()

    def modify_energies(self, mech, Ga_fwd, dG_rxn, state, theta=None):
        """超声**不改能量** —— 它改的是工况本身，见 `modify_state()`。

        这一点必须写明：否则后人会以为"调用了 sono driver 就算上了超声"。
        """
        return np.asarray(Ga_fwd, float).copy(), np.asarray(dG_rxn, float).copy()

    def modify_state(self, state):
        """把一个名义工况拆成**多个加权工况**，返回 `[(state_i, w_i), …]`，`Σw_i = 1`。

        三种给法（从一般到特殊）：
          · `sono_zones = [{T_K:…, P_bar:…, weight:…}, …]`  任意多区
          · `hotspot_* ` + `interfacial_*`                   Suslick 两点模型
          · 只给 `hotspot_*`                                  最简两态近似
        `sono_duty < 1` 时整体再按占空比与"超声关闭"的本体态混合。
        """
        zones = self.get_param(state, "sono_zones")
        duty = float(self.get_param(state, "sono_duty", 1.0))
        if not 0.0 <= duty <= 1.0:
            raise ValueError("sono_duty=%g 不在 [0,1]" % duty)

        if zones:
            out, tot = [], 0.0
            for z in zones:
                w = float(z.get("weight", 0.0))
                if w < 0:
                    raise ValueError("sono_zones 里出现负权重 %g" % w)
                tot += w
                st = dict(state)
                if "T_K" in z:
                    st["T_K"] = float(z["T_K"])
                    self._warn_T(float(z["T_K"]))
                if "P_bar" in z:
                    st["P_bar_sono"] = float(z["P_bar"])
                out.append((st, w))
            if tot > 1.0 + 1e-12:
                raise ValueError("sono_zones 的权重之和 = %g > 1" % tot)
            if tot < 1.0 - 1e-12:            # 余下的份额是本体态
                out.append((dict(state), 1.0 - tot))
            return self._apply_duty(out, state, duty)

        f_hot = self.get_param(state, "hotspot_fraction")
        f_int = self.get_param(state, "interfacial_fraction")
        if f_hot is None and f_int is None:
            warnings.warn(
                "sono driver 已启用，但既没给 hotspot_fraction 也没给 sono_zones\n"
                "  ⇒ 退化为单一工况，**等于完全没考虑超声**（结果会与纯热催化相同）。\n"
                "  超声的效应全在「工况非均匀」这一层，不在能量修正这一层。",
                RuntimeWarning, stacklevel=2)
            return [(dict(state), 1.0)]

        out, tot = [], 0.0
        for frac, tkey in ((f_hot, "hotspot_T_K"), (f_int, "interfacial_T_K")):
            if frac is None:
                continue
            w = float(frac)
            if not 0.0 <= w <= 1.0:
                raise ValueError("份额 %g 不在 [0,1]" % w)
            tot += w
            st = dict(state)
            T = self.get_param(state, tkey)
            if T is not None:
                st["T_K"] = float(T)
                self._warn_T(float(T))
            if tkey == "hotspot_T_K":
                P = self.get_param(state, "hotspot_P_bar")
                if P is not None:
                    st["P_bar_sono"] = float(P)
            out.append((st, w))
        if tot > 1.0 + 1e-12:
            raise ValueError(
                "热点份额 + 界面份额 = %g > 1 —— 它们是**同一个体系的不同部分**，"
                "加起来不能超过 1。" % tot)
        out.append((dict(state), 1.0 - tot))        # 本体
        return self._apply_duty(out, state, duty)

    @staticmethod
    def _apply_duty(zones, state, duty):
        """占空比：`duty` 的时间开超声（上面那套分区），其余时间全是本体态。"""
        if duty >= 1.0 - 1e-15:
            return zones
        out = [(st, w * duty) for st, w in zones if w * duty > 0.0]
        out.append((dict(state), 1.0 - duty))
        return out

    @staticmethod
    def _warn_T(T):
        lo1, hi1 = T_HOTSPOT_MULTIBUBBLE
        lo2, hi2 = T_HOTSPOT_SINGLE
        if not (lo1 <= T <= hi2):
            warnings.warn(
                "热点温度 %g K 落在文献区间之外（多气泡 %g–%g K，单气泡 SBSL %g–%g K）。\n"
                "  不是不能用，但请说明依据 —— 速率对温度是**指数敏感**的，"
                "这个数字基本决定结论。\n  %s" % (T, lo1, hi1, lo2, hi2, CITATION),
                RuntimeWarning, stacklevel=3)

    def observables(self, mech, solution, state):
        return {}

    def validate(self, mech):
        return None


# ══════════════════ 从实测声化学产率标定热点份额 ══════════════════
def calibrate_fraction_from_yield(measured_rate, single_zone_rate,
                                  bulk_rate=0.0, name="热点份额"):
    r"""由**实测总速率**反推热点份额 —— 文献那条路的直接实现。

    文献做法：用化学剂量计（H₂O₂ / Weissler 法等）测出声化学产率，
    再配合单泡摩尔产率，反推活性气泡的数密度。本函数是它在 MKM 里的对应物：

        r_measured = w·r_hot + (1 − w)·r_bulk   ⇒   w = (r_meas − r_bulk)/(r_hot − r_bulk)

    ★ **这比"猜一个份额"诚实得多**：份额是从实验标出来的，不是拟合出来的自由参数。
    返回 `(w, 说明)`；若反推出的 w 不在 [0,1]，**报错而不是截断** ——
    那说明热点温度假设有问题，截断只会把矛盾藏起来。
    """
    r_m, r_h, r_b = float(measured_rate), float(single_zone_rate), float(bulk_rate)
    if abs(r_h - r_b) < 1e-300:
        raise ValueError("热点速率与本体速率几乎相等，无法由总速率反推份额")
    w = (r_m - r_b) / (r_h - r_b)
    if not 0.0 <= w <= 1.0:
        raise ValueError(
            "反推出的%s = %.4g 不在 [0,1] —— **不截断，报错**。\n"
            "  这说明假设的热点条件与实测对不上：\n"
            "    w > 1 ⇒ 假设的热点温度偏低（单靠热点解释不了这么高的速率）\n"
            "    w < 0 ⇒ 实测比本体还慢（有抑制效应？或本体速率给错了）\n"
            "  实测 %.4g，热点态 %.4g，本体态 %.4g" % (name, w, r_m, r_h, r_b))
    return w, ("%s = %.4g（由实测 %.4g 反推：热点态 %.4g，本体态 %.4g）"
               % (name, w, r_m, r_h, r_b))
