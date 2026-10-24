# -*- coding: utf-8 -*-
# KinSTEP - a modality-agnostic microkinetic solver.
# Copyright (C) 2026 The KinSTEP authors.
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.  It is distributed WITHOUT ANY WARRANTY; see the GNU General Public
# License (LICENSE) for details.
r"""[English] Provenance tags and the unverified-data gate: every number carries where it came
from, and an energy of untrusted origin refuses to solve.

(The commentary below is in Chinese: it records why each choice was made
and which measurement forced it.)

参数**来源标签**与"未核实数据"闸。

为什么要这个
------------
MKM 的输入里，**机理是公开知识，能量是私有数据**。两者的可信度天差地别：
机理写错了肉眼能看出来，能量填错了**只会得到一条看着很合理的火山曲线**。

所以每个数值都带来源，并且：

    ★ 只要有任何一个【能量】参数的来源是 "llm_suggested" 或 "guess"，
      求解器**拒绝运行**，并逐条列出是哪一步、要用户怎么处理。

这条闸同时服务三件事：
  ① 自然语言插件（它只许产出结构与问题，**不许产出能量数字**）；
  ② 用户自己随手填的占位值（`0.5 # 先随便填一个` 这种）；
  ③ 模板库给的"典型值"（可以跑，但输出里必须标明不是你的体系的值）。
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Optional

#: 来源等级：越靠前越可信
SOURCE_LEVELS = (
    "user_given",      # 用户明确给的
    "dft_file",        # 从 DFT 输出读的（带文件路径）
    "literature",      # 文献（**必须带 DOI**）
    "experiment",      # 实验测的
    "default",         # 本包的物理默认值（如 beta_ct=0.5）
    "template",        # 模板库的典型值 —— 能跑，但不是你的体系
    "llm_suggested",   # ★ 自然语言层建议 —— **能量类一律拦截**
    "guess",           # ★ 明确标注的猜测 —— 一律拦截
    "REQUIRED",        # 还没有值
)

#: 这些来源的**能量**参数会被闸拦住
UNTRUSTED_FOR_ENERGY = {"llm_suggested", "guess", "REQUIRED"}

#: 哪些参数算"能量类"（拦截只针对它们；beta_ct=0.5 这种有公认默认值的放行）
ENERGY_KEYS = {"Ga_fwd", "Ga_rev", "dG_rxn", "Ea", "G0", "H0", "desc_energy"}


class UnverifiedDataError(RuntimeError):
    """存在未核实的能量数据，拒绝求解。"""


@dataclass
class Val:
    """一个带来源的数值。

    用法::

        Val(1.21, "dft_file", detail="./OUTCAR")
        Val.required("第 3 步的正向活化自由能")
        Val(0.5, "default")
    """
    value: Optional[float] = None
    source: str = "REQUIRED"
    detail: str = ""           # 文件路径 / DOI / 说明
    verified: bool = False     # 用户显式确认过
    ask: str = ""              # 还缺时，给用户看的提示

    def __post_init__(self):
        if self.source not in SOURCE_LEVELS:
            raise ValueError("未知来源 '%s'；可选：%s" % (self.source, list(SOURCE_LEVELS)))
        if self.source in ("user_given", "dft_file", "experiment"):
            self.verified = True

    @classmethod
    def required(cls, ask=""):
        return cls(None, "REQUIRED", ask=ask)

    @property
    def is_missing(self):
        return self.value is None

    def trusted_for_energy(self):
        return self.verified or self.source not in UNTRUSTED_FOR_ENERGY

    def __float__(self):
        if self.value is None:
            raise UnverifiedDataError("该参数还没有值：%s" % (self.ask or "(无提示)"))
        return float(self.value)

    def __repr__(self):
        v = "缺失" if self.value is None else ("%g" % self.value)
        tag = self.source + ("✓" if self.verified else "")
        return "Val(%s, %s%s)" % (v, tag, (", " + self.detail) if self.detail else "")


def as_val(x, default_source="user_given"):
    """把裸数字 / dict / Val 统一成 Val。让用户可以直接写 `1.21`。"""
    if isinstance(x, Val):
        return x
    if x is None:
        return Val.required()
    if isinstance(x, dict):
        return Val(**x)
    return Val(float(x), default_source)


def audit(items, strict_energy=True):
    """检查一批 (标签, 键, Val)，返回 (缺失列表, 不可信列表, 默认值列表)。"""
    missing, untrusted, defaulted = [], [], []
    for label, key, v in items:
        if v.is_missing:
            missing.append((label, key, v))
        elif key in ENERGY_KEYS and strict_energy and not v.trusted_for_energy():
            untrusted.append((label, key, v))
        elif v.source in ("default", "template"):
            defaulted.append((label, key, v))
    return missing, untrusted, defaulted


def gate(items, strict_energy=True):
    """求解前的闸：缺值或能量不可信 ⇒ 抛异常，并**把怎么修说清楚**。"""
    missing, untrusted, defaulted = audit(items, strict_energy)
    if not missing and not untrusted:
        return defaulted
    from .i18n import t
    msg = []
    if missing:
        msg.append(t("【缺 %d 项】" % len(missing),
                     "[%d value(s) missing]" % len(missing)))
        for label, key, v in missing:
            msg.append("  " + t("%s 的 %s —— %s" % (label, key,
                                                    v.ask or "没有提示"),
                                "%s / %s -- %s" % (label, key,
                                                   v.ask or "no hint available")))
    if untrusted:
        msg.append(t("【%d 项能量数据未核实，拒绝求解】" % len(untrusted),
                     "[%d energy value(s) unverified; refusing to solve]"
                     % len(untrusted)))
        for label, key, v in untrusted:
            det = ("（" + v.detail + "）") if v.detail else ""
            msg.append("  " + t("%s 的 %s = %g  ← 来源 '%s'%s"
                                % (label, key, v.value, v.source, det),
                                "%s / %s = %g  <- provenance '%s'%s"
                                % (label, key, v.value, v.source,
                                   (" (" + v.detail + ")") if v.detail else "")))
        msg.append("  " + t("⇒ 能量决定一切结论。请用 DFT/文献值替换；",
                            "=> the energies decide every conclusion. Replace them "
                            "with DFT or literature values;"))
        msg.append("    " + t("若确实想用这个值，显式标注：Val(值, 'user_given') 或 verified=True。",
                              "if you really mean to use this number, say so "
                              "explicitly: Val(x, 'user_given') or verified=True."))
    raise UnverifiedDataError("\n".join(msg))


def report_defaults(defaulted):
    """把"用了默认值"的清单格式化 —— **默认值必须可见**，这是与 CATKINAS 的关键差别。"""
    from .i18n import t
    if not defaulted:
        return t("（没有使用任何默认值）", "(no default values were used)")
    out = [t("以下 %d 项使用了默认/模板值（不改也能跑，但请确认合理）：" % len(defaulted),
             "%d item(s) below use a default/template value - the run will work, "
             "but please confirm they are sensible:" % len(defaulted))]
    for label, key, v in defaulted:
        out.append("  %-24s %-12s = %-10g [%s]" % (label, key, v.value, v.source))
    return "\n".join(out)
