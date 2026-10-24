# -*- coding: utf-8 -*-
# KinSTEP - a modality-agnostic microkinetic solver.
# Copyright (C) 2026 The KinSTEP authors.
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.  It is distributed WITHOUT ANY WARRANTY; see the GNU General Public
# License (LICENSE) for details.
r"""[English] The input file (.kinstep.toml): a plain-text model description, so that nobody
has to write an input script.

(The commentary below is in Chinese: it records why each choice was made
and which measurement forced it.)

输入文件（`.kinstep.toml`）—— 把"别再填复杂输入脚本"落到实处。

格式为什么选 TOML
-----------------
① **`tomllib` 在 Python 3.11+ 的标准库里** ⇒ 读文件零新依赖，与"核心只要 numpy+scipy"一致；
② 它是通行格式，编辑器有高亮、有校验，不必让用户学一套本包自创的语法；
③ **机理本身仍然是自由文本**（多行字符串），照着纸上写 —— 这是本包的首要承诺，
   不能因为用了配置格式就把反应式拆成键值对。

扩展名用 **`.kinstep.toml`**，不用 `.mkm`（那是 CatMAP 的后缀，见 ATTRIBUTION.md 的避让条款）。

最小可跑的文件长这样
--------------------
```toml
[mechanism]
text = '''
H2 + 2* <-> 2H*
N2 + 2* <-> 2N*
N* + H* <-> NH* + *
NH* + H* <-> NH2* + *
NH2* + H* <-> NH3 + 2*
'''

[conditions]
T_K = 673.0
p_gas = { H2 = 75.0, N2 = 25.0, NH3 = 1.0 }

# 逐步能量：顺序与机理各行一一对应
[[steps]]
Ga_fwd = 1.310101
dG_rxn = 0.397301
# …（每步一个 [[steps]]）
```

要扫描就再加（两种给法二选一：逐步能量 **或** scaling+BEP）
```toml
[scan]
descriptors = ["E_N"]
axes = [[-2.0, 2.0, 401]]        # [起, 止, 点数]

[[scaling]]                       # 每步一个，顺序同机理
bep   = [0.54, 0.60]              # Ga = α·ΔG⁰ + β
scale = [[0.30], -0.68]           # ΔG⁰ = Σ a_d·desc_d + b
```
"""
from __future__ import annotations

import os
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from . import mech as _mech
from .provenance import Val
from .rates import Energetics, StepEnergy
from .scaling import ScalingModel, from_table as _scaling_from_table

SUFFIX = ".kinstep.toml"


class InputError(ValueError):
    """输入文件有问题。**报错要说清楚怎么改**，不只是说哪里错。"""


@dataclass
class Spec:
    """一份输入文件解析后的全部内容。"""
    path: str
    mech: _mech.Mechanism
    conditions: Dict[str, Any]
    p_gas: Dict[str, float]
    c_liq: Dict[str, float]
    energetics: Optional[Energetics] = None
    scaling: Optional[ScalingModel] = None
    scan_axes: Optional[List[np.ndarray]] = None
    drivers: List[Tuple[str, Dict]] = field(default_factory=list)
    solver: Dict[str, Any] = field(default_factory=dict)
    raw: Dict[str, Any] = field(default_factory=dict)
    notes: List[str] = field(default_factory=list)

    @property
    def is_scan(self):
        return self.scan_axes is not None


def _need(d, key, where, hint=""):
    if key not in d:
        raise InputError("%s 里缺 `%s`。%s" % (where, key, hint))
    return d[key]


def load(path: str) -> Spec:
    """读一份 `.kinstep.toml`。每一处报错都带**怎么改**。"""
    # ★ 3.11 起 `tomllib` 在标准库里；3.9/3.10 上回退到 `tomli`（同一套 API，
    #   由 pyproject 的条件依赖自动装上）。
    #   为什么要回退：很多机群默认还是 3.9/3.10，为了一个读配置的函数把整个包的
    #   最低版本钉在 3.11，对使用者是实打实的门槛，而代价只是一行 try/except。
    try:
        import tomllib
    except ModuleNotFoundError:                      # Python < 3.11
        try:
            import tomli as tomllib                  # type: ignore
        except ModuleNotFoundError:
            raise InputError(
                "读 .kinstep.toml 需要 TOML 解析器：Python 3.11+ 自带 tomllib；"
                "更早的版本请装 tomli（pip install tomli）") from None
    if not os.path.exists(path):
        from .i18n import t as _t
        raise InputError(_t("找不到文件：%s", "file not found: %s") % path)
    with open(path, "rb") as f:
        try:
            raw = tomllib.load(f)
        except Exception as e:
            raise InputError("TOML 解析失败：%s\n  常见原因：多行字符串要用三引号 '''…'''；"
                             "行内表用 { a = 1, b = 2 }" % e) from None

    # ── 机理 ──
    mb = _need(raw, "mechanism", "文件", "写成：\n  [mechanism]\n  text = '''\n  A + * <-> A*\n  '''")
    txt = _need(mb, "text", "[mechanism]", "它是多行字符串，每行一个基元反应。")
    kw = {}
    if "species_formula" in mb:
        kw["species_formula"] = mb["species_formula"]
    if "exempt_balance" in mb:
        kw["exempt_balance"] = mb["exempt_balance"]
    try:
        m = _mech.from_text(txt, **kw)
    except Exception as e:
        # ★ 机理层的异常（MechanismError / NameCollisionError）统一包成 InputError，
        #   这样 CLI 的 check 能走同一条"输入有问题 + 怎么改"的呈现路径，
        #   而不是掉进通用异常分支（首版实测：不配平的机理被报成"解析失败"，
        #   信息其实是全的，但呈现不一致、也没带上定位信息）。
        raise InputError("[mechanism] 里的机理有问题：\n%s\n"
                         "  （机理文本来自 %s 的 [mechanism].text）" % (e, path)) from None

    # ── 工况 ──
    cond = dict(raw.get("conditions", {}))
    p_gas = dict(cond.pop("p_gas", {}) or {})
    c_liq = dict(cond.pop("c_liq", {}) or {})
    if "T_K" not in cond:
        raise InputError("[conditions] 里缺 `T_K`（温度，K）。\n"
                         "  速率常数对温度是**指数敏感**的，本包不设默认值。")

    # ── driver ──
    drivers = []
    for dname, dcfg in (raw.get("drivers", {}) or {}).items():
        drivers.append((dname, dict(dcfg) if isinstance(dcfg, dict) else {}))

    # ── 能量：逐步给 或 scaling+BEP ──
    steps = raw.get("steps", [])
    scal = raw.get("scaling", [])
    if steps and scal:
        raise InputError("[[steps]] 与 [[scaling]] 只能给一种：\n"
                         "  · [[steps]]   —— 直接给每步的 Ga_fwd / dG_rxn（单点计算）\n"
                         "  · [[scaling]] —— 给 scaling+BEP 系数（描述符扫描）")
    en = sc = None
    if steps:
        if len(steps) != m.n_step:
            raise InputError("机理有 %d 步，但只给了 %d 个 [[steps]]。**每步都要给，不替你猜。**\n"
                             "  机理各行：\n    %s"
                             % (m.n_step, len(steps),
                                "\n    ".join("(%d) %s" % (i + 1, s.raw)
                                              for i, s in enumerate(m.steps))))
        tab = {}
        for i, sp in enumerate(steps):
            spec = dict(sp)
            src = spec.pop("source", "user_given")
            mode = spec.pop("mode", "tst")
            d = dict(mode=mode)
            for k in ("Ga_fwd", "dG_rxn", "A_pre", "kf0", "kr0",
                      "mass_amu", "sticking", "A_site_m2"):
                if k in spec:
                    d[k] = Val(float(spec.pop(k)), src)
            if spec:
                raise InputError("第 %d 个 [[steps]] 里有认不得的键：%s" % (i + 1, list(spec)))
            tab[i + 1] = StepEnergy(**d)
        en = Energetics(m, tab)
    if scal:
        sb = raw.get("scan", {})
        names = sb.get("descriptors") or ["desc"]
        if len(scal) != m.n_step:
            raise InputError("机理有 %d 步，但只给了 %d 个 [[scaling]]" % (m.n_step, len(scal)))
        tab = {}
        for i, sp in enumerate(scal):
            d = {}
            if "bep" in sp:
                d["bep"] = tuple(float(v) for v in sp["bep"])
            if "Ga" in sp:
                d["Ga"] = float(sp["Ga"])
            if "scale" not in sp:
                raise InputError("第 %d 个 [[scaling]] 缺 `scale = [[斜率…], 截距]`" % (i + 1))
            sl, ic = sp["scale"]
            d["scale"] = (np.atleast_1d(np.asarray(sl, float)), float(ic))
            tab[i + 1] = d
        off_Ga = raw.get("thermo", {}).get("offset_Ga")
        off_dG = raw.get("thermo", {}).get("offset_dG")
        sc = _scaling_from_table(m, tab, names, offset_Ga=off_Ga, offset_dG=off_dG)

    # ── 扫描轴 ──
    axes = None
    if "scan" in raw:
        ax = _need(raw["scan"], "axes", "[scan]", "写成 axes = [[起, 止, 点数]]")
        axes = []
        for a in ax:
            if len(a) != 3:
                raise InputError("[scan].axes 的每一项要写成 [起, 止, 点数]，收到 %s" % (a,))
            axes.append(np.linspace(float(a[0]), float(a[1]), int(a[2])))
        if sc is None:
            raise InputError("给了 [scan] 却没给 [[scaling]] —— 扫描需要"
                             "描述符到能量的关系。")
        if len(axes) != sc.n_desc:
            raise InputError("[scan].descriptors 有 %d 个，但 axes 给了 %d 条"
                             % (sc.n_desc, len(axes)))
    if en is None and sc is None:
        raise InputError("既没有 [[steps]] 也没有 [[scaling]] —— 不知道每步的能量。\n"
                         "  单点算：每步一个 [[steps]]，给 Ga_fwd 与 dG_rxn\n"
                         "  扫描　：每步一个 [[scaling]]，给 bep 与 scale，并加 [scan]")

    return Spec(path=path, mech=m, conditions=cond, p_gas=p_gas, c_liq=c_liq,
                energetics=en, scaling=sc, scan_axes=axes, drivers=drivers,
                solver=dict(raw.get("solver", {})), raw=raw)


def template(kind: str = "single") -> str:
    """生成一份可直接跑的模板（`kinstep new` 用它当骨架）。

    ★ 这份骨架的注释也必须跟着 `KINSTEP_LANG` 走 —— 坑 29 的同族。
      `new --template xxx` 走 `templates.Template.to_toml()`，
      而**不带 --template 时走这里**；两条路径要各自本地化，
      只修了一条就会出现「回显英文、文件中文」。T16 ④ 现在两条都扫。
    """
    from .i18n import t
    _kind = (t("单点稳态", "single steady state") if kind == "single"
             else t("描述符扫描", "descriptor scan"))
    _c1 = t("# 机理照着纸上写：`*` 是位点，`<->` 可逆，`->` 不可逆，`(l)` 液相，`e-` 电子。",
            "# Write the mechanism just as you would on paper: `*` is a site, "
            "`<->` reversible, `->` irreversible, `(l)` liquid, `e-` electron.")
    _c2 = t("# nu_elec（电子计量系数）**自动从方程里数**，不用填。",
            "# nu_elec (the electron stoichiometry) is **counted from the equation "
            "automatically** - you do not fill it in.")
    head = '''%s
%s
%s

[mechanism]
text = \'\'\'
H2 + 2* <-> 2H*
N2 + 2* <-> 2N*
N* + H* <-> NH* + *
NH* + H* <-> NH2* + *
NH2* + H* <-> NH3 + 2*
\'\'\'

[conditions]
T_K = 673.0
p_gas = { H2 = 75.0, N2 = 25.0, NH3 = 1.0 }
''' % (t("# KinSTEP 输入文件（%s）" % _kind,
         "# KinSTEP input file (%s)" % _kind), _c1, _c2)
    if kind == "single":
        body = "\n" + t("# 逐步能量（eV），顺序与机理各行一一对应",
                        "# Per-step energies (eV), in the same order as the "
                        "mechanism lines") + "\n"
        for i, (g, d) in enumerate([(1.310101, 0.397301), (3.044256, 1.504256),
                                    (1.4032, -0.18), (1.4665, -0.45), (1.8127, -0.405974)]):
            body += "\n[[steps]]\nGa_fwd = %.6f\ndG_rxn = %+.6f\n" % (g, d)
        body += '\n[solver]\nprecision = "auto"      # %s\n' % t(
            "double / auto / high —— auto 会在需要时自动升高精度",
            "double / auto / high - auto escalates the precision when it is needed")
    else:
        body = '''
[scan]
descriptors = ["E_N"]
axes = [[-2.0, 2.0, 401]]        # %s

[thermo]
# %s
offset_Ga = [1.077301, 1.504256, 0.0, 0.0, 0.0]
offset_dG = [1.077301, 1.504256, 0.0, 0.0, -1.575974]

''' % (t("[起, 止, 点数]", "[from, to, number of points]"),
       t("热化学常数偏移（逐步，eV）。T10 的热力学库做好后可以自动算；在那之前手工给。",
         "Thermochemical offsets (per step, eV). The thermo module can compute "
         "these; until then supply them by hand."))
        for bep, sl, ic in [((0.54, 0.60), 0.30, -0.68), ((0.78, 1.54), 2.00, 0.00),
                            ((0.26, 1.45), -0.15, -0.18), ((0.43, 1.66), -0.45, -0.45),
                            ((0.31, 1.45), -0.85, 1.17)]:
            # α·ΔG⁰·β 是通用数学符号、不是汉字，中英两版都保留原样
            body += ("[[scaling]]\nbep   = [%.2f, %.2f]      # Ga = α·ΔG⁰ + β\n"
                     "scale = [[%.2f], %+.2f]   # ΔG⁰ = a·E_N + b\n\n" % (bep[0], bep[1], sl, ic))
        body += '[solver]\nprecision = "auto"\n'
    return head + body
