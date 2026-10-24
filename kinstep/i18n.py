# -*- coding: utf-8 -*-
# KinSTEP - a modality-agnostic microkinetic solver.
# Copyright (C) 2026 The KinSTEP authors.
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.  It is distributed WITHOUT ANY WARRANTY; see the GNU General Public
# License (LICENSE) for details.
r"""[English] Interface language switch (Chinese / English). KINSTEP_LANG=en gives English
output; one regression asserts that no CJK character survives in that mode.

(The commentary below is in Chinese: it records why each choice was made
and which measurement forced it.)

界面语言开关（中文 / English）。

## 为什么要有这个模块
本包的命令行输出原本只有中文。投 CPC 这类国际期刊时，**审稿人要真的把程序跑起来**，
一个只说中文的界面是真实的使用门槛，不只是插图好不好看的问题。

## 设计取舍（三条，都是为了不破坏既有的 75 项回归）
1. **默认仍是中文**。`KINSTEP_LANG` 不设或设成别的值 ⇒ 中文。
   既有回归逐字比对中文输出，默认不变它们就不受影响。
2. **行内双语 `t("中文", "English")`**，不另起一张 key 表。
   key 表要求改写所有调用点、而且 key 和内容会慢慢对不上；
   行内双语把两种说法摆在同一行，改了一个忘了另一个**肉眼就能看见**。
3. **惰性读环境变量**，并提供 `set_lang()`，这样测试里可以来回切换，
   不必在导入之前就把环境变量设好。

用法::

    from .i18n import t
    print(t("求解完成", "solved"))

    KINSTEP_LANG=en kinstep check her.kinstep.toml
"""
from __future__ import annotations

import os

_LANG = None


def lang() -> str:
    """当前界面语言：``"zh"`` 或 ``"en"``。首次调用时读环境变量。"""
    global _LANG
    if _LANG is None:
        v = (os.environ.get("KINSTEP_LANG") or "").strip().lower()
        _LANG = "en" if v.startswith("en") else "zh"
    return _LANG


def set_lang(value) -> str:
    """显式设定语言（测试与库调用用）。返回设定后的值。"""
    global _LANG
    v = str(value).strip().lower()
    _LANG = "en" if v.startswith("en") else "zh"
    return _LANG


def reset_lang() -> None:
    """清掉缓存，下次 `lang()` 重新读环境变量。"""
    global _LANG
    _LANG = None


def t(zh: str, en: str) -> str:
    """按当前语言二选一。两个参数都必须给 —— 漏了一个就会在这里报错，
    而不是在用户的屏幕上留下半句中文。"""
    return en if lang() == "en" else zh
