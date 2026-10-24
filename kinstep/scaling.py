# -*- coding: utf-8 -*-
# KinSTEP - a modality-agnostic microkinetic solver.
# Copyright (C) 2026 The KinSTEP authors.
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.  It is distributed WITHOUT ANY WARRANTY; see the GNU General Public
# License (LICENSE) for details.
r"""[English] Descriptor to per-step energies: scaling relations and BEP relations kept as two
separate layers.

(The commentary below is in Chinese: it records why each choice was made
and which measurement forced it.)

描述符 → 逐步能量（T9 的第一层）。

这一层在做什么
--------------
火山图的横轴是**描述符**（如 N 的吸附能 E_N），而求解器要的是**逐步的 Ga 与 ΔG**。
中间这一步靠两组公认的线性关系：

    ① **scaling 关系**：ΔG⁰_s = Σ_d a[s,d]·desc_d + b[s]
       —— 吸附能之间的线性相关（Abild-Pedersen 2007 一类）
    ② **BEP 关系**：Ga_s = α[s]·ΔG⁰_s + β[s]
       —— Brønsted–Evans–Polanyi，势垒与反应能线性相关

★ 两层都有**精确真值**（这是 T9 不必依赖别的程序的原因）
--------------------------------------------------------
CATKINAS 的 `INPUT_curve.m` 把系数直接写在能量表里，与它 401 点的输出逐点核对：

    scaling  ΔH⁰ = a·E_N + b        401 点最大绝对差 **7.550e-15**
    BEP      Ea  = α·ΔH⁰ + β        401 点最大绝对差 **8.882e-16**

⇒ 本模块的两层公式与 CATKINAS 完全同口径（本包 `tests/t9_scan.py` 有逐点回归）。

★ 常数偏移 `offset_Ga` / `offset_dG`
------------------------------------
热化学校正（气相熵、碰撞理论的平动项等）在**固定温度下是逐步的常数**，
与描述符无关 —— 这一点也是实测的：CATKINAS 的 `Ean−Eao` 与 `G0n−G0o`
在 401 个点上跨度 < 1.4e-15：

    ΔGa = [+1.077301, +1.504256, 0, 0, 0] eV
    ΔG0 = [+1.077301, +1.504256, 0, 0, −1.575974] eV

所以把它们做成**可加的常数向量**是对的，不必混进 scaling 系数里
（混进去就看不出哪部分是化学、哪部分是热化学了）。这两个常数由 T10 的热力学库产出；
在 T10 之前可以手工给（回归时就是这么做的）。

★ 为什么 `Ga` 允许绕过 BEP 直接给
----------------------------------
真实输入里常有混搭：多数步用 BEP 估势垒，而某一步有 DFT 算好的真势垒
（CATKINAS 的 `Mr_NH3 = -2` 就是"NH3 脱附的势垒直接给"）。
⇒ `Ga_direct` 里非 nan 的步**跳过 BEP**，直接用给的值（仍加 offset）。
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from .mech import Mechanism


class ScalingError(ValueError):
    """scaling / BEP 的系数给得不对（维数、缺项、键名）。"""


@dataclass
class ScalingModel:
    """描述符 → (Ga_fwd, dG_rxn)。`n_step` 步、`n_desc` 个描述符。"""
    dG_coeff: np.ndarray            # (n_step, n_desc)  scaling 斜率 a
    dG_const: np.ndarray            # (n_step,)         scaling 截距 b
    bep_alpha: np.ndarray           # (n_step,)         BEP 斜率 α
    bep_beta: np.ndarray            # (n_step,)         BEP 截距 β
    desc_names: List[str]
    Ga_direct: np.ndarray = None    # (n_step,) 非 nan 的步跳过 BEP
    offset_Ga: np.ndarray = None    # (n_step,) 热化学常数偏移
    offset_dG: np.ndarray = None
    notes: List[str] = field(default_factory=list)

    def __post_init__(self):
        n = len(self.dG_const)
        for nm in ("dG_coeff", "dG_const", "bep_alpha", "bep_beta"):
            setattr(self, nm, np.asarray(getattr(self, nm), float))
        if self.dG_coeff.ndim == 1:
            self.dG_coeff = self.dG_coeff.reshape(n, -1)
        if self.dG_coeff.shape != (n, len(self.desc_names)):
            raise ScalingError(
                "dG_coeff 形状应为 (%d 步, %d 个描述符)，收到 %s"
                % (n, len(self.desc_names), self.dG_coeff.shape))
        for nm in ("Ga_direct", "offset_Ga", "offset_dG"):
            v = getattr(self, nm)
            setattr(self, nm, np.full(n, np.nan) if v is None and nm == "Ga_direct"
                    else (np.zeros(n) if v is None else np.asarray(v, float)))

    @property
    def n_step(self):
        return len(self.dG_const)

    @property
    def n_desc(self):
        return len(self.desc_names)

    def energies(self, desc) -> Tuple[np.ndarray, np.ndarray]:
        """描述符值 → `(Ga_fwd, dG_rxn)`，单位 eV。

        `desc` 是长度 `n_desc` 的序列（1D 扫描时就是一个数）。
        """
        d = np.atleast_1d(np.asarray(desc, float))
        if d.size != self.n_desc:
            raise ScalingError("需要 %d 个描述符 %s，收到 %d 个"
                               % (self.n_desc, self.desc_names, d.size))
        dG = self.dG_coeff @ d + self.dG_const           # ① scaling
        Ga = self.bep_alpha * dG + self.bep_beta         # ② BEP
        Ga = np.where(np.isnan(self.Ga_direct), Ga, self.Ga_direct)  # 直接给的优先
        return Ga + self.offset_Ga, dG + self.offset_dG

    def describe(self) -> str:
        L = ["scaling 模型：%d 步，描述符 %s" % (self.n_step, self.desc_names)]
        for s in range(self.n_step):
            terms = " ".join("%+.4g·%s" % (self.dG_coeff[s, k], self.desc_names[k])
                             for k in range(self.n_desc))
            direct = "" if np.isnan(self.Ga_direct[s]) else \
                "  [Ga 直接给 = %.4g，跳过 BEP]" % self.Ga_direct[s]
            L.append("  步%d  ΔG⁰ = %s %+.4g   |   Ga = %+.4g·ΔG⁰ %+.4g%s"
                     % (s + 1, terms, self.dG_const[s],
                        self.bep_alpha[s], self.bep_beta[s], direct))
        if np.any(self.offset_Ga != 0) or np.any(self.offset_dG != 0):
            L.append("  热化学常数偏移  ΔGa = %s" % np.array2string(self.offset_Ga, precision=6))
            L.append("                  ΔG⁰ = %s" % np.array2string(self.offset_dG, precision=6))
        return "\n".join(L)


def from_table(mech: Mechanism, table: Dict, desc_names: Sequence[str],
               offset_Ga=None, offset_dG=None) -> ScalingModel:
    """从**按步写的表**造 scaling 模型 —— 和机理一样，用户照着纸上写。

    用法（1D，就是 CATKINAS 那张表的写法）::

        from_table(mech, {
            1: dict(bep=(0.54, 0.60), scale=(0.30, -0.68)),
            2: dict(bep=(0.78, 1.54), scale=(2.00,  0.00)),
            ...
        }, desc_names=["E_N"])

    2D 只是把 `scale` 的斜率写成元组::

        1: dict(bep=(0.5, 0.8), scale=((0.3, 0.1), -0.68))

    某步的势垒直接给（跳过 BEP）::

        5: dict(Ga=1.21, scale=(-0.85, 1.17))
    """
    n, nd = mech.n_step, len(desc_names)
    a = np.zeros((n, nd)); b = np.zeros(n)
    al = np.zeros(n); be = np.zeros(n)
    gd = np.full(n, np.nan)
    raw2idx = {st.raw: i for i, st in enumerate(mech.steps)}
    seen = set()
    for key, spec in table.items():
        if isinstance(key, (int, np.integer)):
            i = int(key) - 1
            if not 0 <= i < n:
                raise ScalingError("步号 %s 超出范围（机理 %d 步，步号从 1 数）" % (key, n))
        elif key in raw2idx:
            i = raw2idx[key]
        else:
            raise ScalingError("认不出键 %r；用 1 起的步号或机理原文那一行" % (key,))
        if i in seen:
            raise ScalingError("第 %d 步给了两次" % (i + 1))
        seen.add(i)
        if "scale" not in spec:
            raise ScalingError("第 %d 步缺 scale=(斜率, 截距)" % (i + 1))
        sl, ic = spec["scale"]
        sl = np.atleast_1d(np.asarray(sl, float))
        if sl.size != nd:
            raise ScalingError("第 %d 步的 scaling 斜率要给 %d 个（描述符 %s），收到 %d 个"
                               % (i + 1, nd, list(desc_names), sl.size))
        a[i], b[i] = sl, float(ic)
        if "Ga" in spec:
            gd[i] = float(spec["Ga"])
        elif "bep" in spec:
            al[i], be[i] = (float(v) for v in spec["bep"])
        else:
            raise ScalingError("第 %d 步既没给 bep=(α, β) 也没给 Ga=势垒" % (i + 1))
    missing = [i + 1 for i in range(n) if i not in seen]
    if missing:
        raise ScalingError("这些步没给 scaling/BEP：%s —— **不替你猜**" % missing)
    return ScalingModel(dG_coeff=a, dG_const=b, bep_alpha=al, bep_beta=be,
                        desc_names=list(desc_names), Ga_direct=gd,
                        offset_Ga=offset_Ga, offset_dG=offset_dG)


def fit_scaling(desc_samples, dG_samples):
    """从 (描述符, ΔG⁰) 样本最小二乘拟合 scaling 系数，返回 `(a, b, R²)`。

    给"手上有一组金属的 DFT 数据，想自己定 scaling"的用户。
    **R² 一并返回并要求用户看** —— scaling 关系好不好，决定整张火山图可不可信。
    """
    X = np.atleast_2d(np.asarray(desc_samples, float))
    if X.shape[0] < X.shape[1]:
        X = X.T
    y = np.asarray(dG_samples, float)
    A = np.hstack([X, np.ones((X.shape[0], 1))])
    sol, *_ = np.linalg.lstsq(A, y, rcond=None)
    pred = A @ sol
    ss_res = float(np.sum((y - pred) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
    return sol[:-1], float(sol[-1]), r2
