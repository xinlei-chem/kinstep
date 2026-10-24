# -*- coding: utf-8 -*-
# KinSTEP - a modality-agnostic microkinetic solver.
# Copyright (C) 2026 The KinSTEP authors.
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.  It is distributed WITHOUT ANY WARRANTY; see the GNU General Public
# License (LICENSE) for details.
r"""[English] Mechanism parsing: elementary steps written as on paper into a stoichiometric
matrix, with the electron stoichiometry counted from the equation itself.

(The commentary below is in Chinese: it records why each choice was made
and which measurement forced it.)

机理解析：**照着纸上写的反应式** → 化学计量矩阵。

设计目标（承"别再写复杂输入脚本"）
------------------------------------
用户只写反应式本身，**不写任何 Mode、不写电子数、不写计量矩阵**：

    H2 + 2* <-> 2H*
    H2O + * <-> OH* + H+ + e-          {beta_ct: 0.5}
    CO*1 + O*2 -> CO2 + *1 + *2

约定（都是催化文献里的通行写法，不是本包发明的）
    *          默认位点；  *1 / *2 / *M / *O   多种位点
    CO*        吸附在默认位点的 CO；  CO*1  吸附在 1 号位点
    (g) (l)    显式指定相态；不写则气相
    e-         电子；  H+  质子
    <->        可逆；  ->  不可逆（逆向速率常数置零）
    { ... }    该步的 driver 参数（如 beta_ct）

★ `nu_elec` **自动从方程里数**，不让用户填 —— 它是最容易写错符号的参数
  （产物侧 e- 计数 − 反应物侧 e- 计数；正 = 氧化）。

四道闸（都是"静默错误"，不设闸不会报错）
------------------------------------------
① **撞名**：物种/位点名不得撞保留名（见 naming.py）
② **元素与电荷守恒**：每一步左右两边逐元素 + 总电荷必须配平 —— 抓笔误
③ **位点守恒**：每种位点在每一步左右两边数目必须相等
④ **pH 双重计数**：机理里写了 H+ 却又用 RHE 标度 ⇒ 报错，二选一（在 solve 期检查）
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from .naming import (ELECTRON_TOKENS, HYDROXIDE_TOKENS, PROTON_TOKENS,
                     NameCollisionError, check_user_names, strip_species_prefix)

#: 相态标记
PHASE_SUFFIX = {"(g)": "gas", "(gas)": "gas", "(l)": "liq", "(liq)": "liq",
                "(aq)": "liq", "(c)": "liq", "(p)": "gas"}

#: 元素符号（按长度降序，保证先匹配两字母元素）
_ELEMENTS = sorted(
    ("H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn "
     "Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La "
     "Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po "
     "At Rn Fr Ra Ac Th Pa U Np Pu").split(), key=len, reverse=True)
_ELEM_RE = re.compile("(" + "|".join(_ELEMENTS) + r")(\d*)")


class MechanismError(ValueError):
    """机理写法有问题。**一律报错，不猜用户的意思。**"""


@dataclass
class Species:
    name: str                  # 规范名，如 'CO*1'、'H2O(g)'、'e-'
    formula: str               # 去掉位点与相态，如 'CO'
    phase: str                 # 'gas' | 'liq' | 'surf' | 'site' | 'electron'
    site: Optional[str] = None # 吸附/空位所在位点名（'' = 默认位点）
    charge: int = 0
    #: 是否用 `Sp_`/`E_` 前缀显式声明过（豁免保留名检查、允许非元素式名字）
    declared: bool = False
    #: 元素组成；None = 由 formula 现场解析。非元素式名字必须由 species_formula 给
    composition: Optional[Dict[str, int]] = None

    @property
    def is_surface(self):
        return self.phase in ("surf", "site")


@dataclass
class Step:
    raw: str
    reactants: Dict[str, float]
    products: Dict[str, float]
    reversible: bool = True
    params: Dict[str, float] = field(default_factory=dict)

    def net(self) -> Dict[str, float]:
        d = {k: -v for k, v in self.reactants.items()}
        for k, v in self.products.items():
            d[k] = d.get(k, 0.0) + v
        return {k: v for k, v in d.items() if v != 0.0}


def _parse_formula(formula: str) -> Tuple[Dict[str, int], int]:
    """化学式 → (元素计数, 电荷)。`OH-` → ({'O':1,'H':1}, -1)。"""
    charge = 0
    f = formula
    m = re.search(r"([+-]\d*)$", f)
    if m:
        sign = 1 if m.group(1)[0] == "+" else -1
        mag = m.group(1)[1:]
        charge = sign * (int(mag) if mag else 1)
        f = f[: m.start()]
    counts: Dict[str, int] = {}
    pos = 0
    while pos < len(f):
        mm = _ELEM_RE.match(f, pos)
        if not mm:
            raise MechanismError("化学式 '%s' 在第 %d 位解析不了（'%s'）。\n"
                                 "请用标准元素符号；不规则名字请改写，例如 'Cat' → 'C1a1t1' 不可取，"
                                 "应给它一个纯符号名并在 species_formula 里显式提供组成。"
                                 % (formula, pos, f[pos:pos + 4]))
        el, n = mm.group(1), mm.group(2)
        counts[el] = counts.get(el, 0) + (int(n) if n else 1)
        pos = mm.end()
    return counts, charge


def _parse_token(tok: str, species_formula: Dict[str, str] = None) -> Tuple[float, Species]:
    """'2H*1' → (2.0, Species('H*1', 'H', 'surf', site='1'))

    可选的 `Sp_`/`E_` 前缀在**最前面**先剥掉（`2Sp_Cat*` 也认），
    剥掉后的裸名就是规范名 —— 所以 `E_CO` 与 `CO` 是同一个物种。
    """
    species_formula = species_formula or {}
    tok = tok.strip()
    if not tok:
        raise MechanismError("空物种项")
    m = re.match(r"^(\d+(?:\.\d+)?)\s*(.+)$", tok)
    coef, body = (float(m.group(1)), m.group(2).strip()) if m else (1.0, tok)

    if body in ELECTRON_TOKENS:
        return coef, Species("e-", "", "electron", charge=-1)

    # ★ 显式物种前缀先剥（在相态/位点之前，因为它在最前面）
    body, declared = strip_species_prefix(body)

    phase = None
    for suf, ph in PHASE_SUFFIX.items():
        if body.lower().endswith(suf):
            phase = ph
            body = body[: -len(suf)].strip()
            break

    site = None
    if "*" in body:
        head, _, tail = body.partition("*")
        if "*" in tail:
            raise MechanismError("物种 '%s' 里有多个 '*'，一个物种只能占一个位点标记" % tok)
        site = tail.strip()          # '' = 默认位点
        body = head.strip()
        if body == "":               # 裸 '*' = 空位
            return coef, Species("*" + site, "", "site", site=site)
        comp = _composition_of(body, declared, species_formula, tok)
        return coef, Species(body + "*" + site, body, "surf", site=site,
                             declared=declared, composition=comp)

    comp = _composition_of(body, declared, species_formula, tok)
    charge = 0
    if comp is not None:
        # 组成已知时顺带把电荷解析出来（`OH-` 这种）
        try:
            _, charge = _parse_formula(body)
        except MechanismError:
            charge = 0
    # ★ 带电物种默认**液相**：催化机理里的 H+/OH- 都是水合离子，
    #   "气相裸质子"在物理上说不通，默认成气相会让用户在 p_gas 里找不到它而静默取 0
    #   （这正是 Koper 算例上踩到的坑）
    if phase is None:
        phase = "liq" if charge != 0 else "gas"
    return coef, Species(body + ("(%s)" % phase if phase == "liq" and charge == 0 else ""),
                         body, phase, charge=charge,
                         declared=declared, composition=comp)


def _composition_of(body, declared, species_formula, tok):
    """确定元素组成：`species_formula` 优先；否则按化学式解析。

    非元素式名字（`Cat`/`Dye`/`Nafion`）只有**带前缀声明过**才允许，
    且必须在 `species_formula` 里给组成，或列入 `exempt_balance` —— 二者都没有
    就报错，**不静默跳过元素配平**（配平是本包最有效的查错手段，不能偷偷放弃）。
    """
    if body in species_formula:
        comp, _ = _parse_formula(str(species_formula[body]))
        return comp
    try:
        comp, _ = _parse_formula(body)
        return comp
    except MechanismError:
        if declared:
            return None          # 组成未知，留给 from_text 统一报错/豁免
        raise


def _parse_side(side: str, species_formula: Dict[str, str] = None
                ) -> Tuple[Dict[str, float], Dict[str, Species]]:
    out: Dict[str, float] = {}
    reg: Dict[str, Species] = {}
    for tok in re.split(r"\s\+\s|^\+\s|\s\+$", side.strip()):
        tok = tok.strip()
        if not tok:
            continue
        coef, sp = _parse_token(tok, species_formula)
        out[sp.name] = out.get(sp.name, 0.0) + coef
        reg[sp.name] = sp
    if not out:
        raise MechanismError("反应式一侧是空的：'%s'" % side)
    return out, reg


class Mechanism:
    """解析好的机理。`steps` / `species` / `sites` / `nu`（化学计量矩阵）。"""

    def __init__(self, steps: List[Step], species: Dict[str, Species]):
        self.steps = steps
        self.species = species
        self.sites = sorted({s.site for s in species.values()
                             if s.site is not None})
        self.species_order = sorted(species)
        #: 明确豁免元素配平的物种（由 from_text 填）—— 必须在输出里可见
        self.exempt_balance: List[str] = []
        self.nu = self._build_nu()

    # ── 化学计量矩阵 ν[step, species] ──
    def _build_nu(self):
        idx = {n: i for i, n in enumerate(self.species_order)}
        nu = np.zeros((len(self.steps), len(self.species_order)))
        for si, st in enumerate(self.steps):
            for n, v in st.net().items():
                nu[si, idx[n]] = v
        return nu

    @property
    def n_step(self):
        return len(self.steps)

    @property
    def surface_species(self):
        return [n for n in self.species_order if self.species[n].phase == "surf"]

    @property
    def gas_species(self):
        return [n for n in self.species_order if self.species[n].phase == "gas"]

    @property
    def liq_species(self):
        return [n for n in self.species_order if self.species[n].phase == "liq"]

    def vacancy_of(self, site):
        return "*" + site

    def has_proton(self):
        return any(n in PROTON_TOKENS or self.species[n].formula == "H"
                   and self.species[n].charge == 1 for n in self.species_order)

    def has_electron(self):
        return "e-" in self.species

    def summary(self) -> str:
        from .i18n import t
        none = t("无", "none")
        sites = [(t("默认", "default") if s == "" else s) for s in self.sites] or none
        L = [t("机理：%d 步，%d 个物种，位点 %s"
               % (self.n_step, len(self.species_order), sites),
               "mechanism: %d steps, %d species, site type(s) %s"
               % (self.n_step, len(self.species_order), sites))]
        for i, st in enumerate(self.steps):
            nu_e = st.params.get("nu_elec", 0.0)
            tag = (t("  [电化学 nu_elec=%+g]", "  [electrochemical, nu_elec=%+g]")
                   % nu_e) if nu_e else ""
            L.append("  (%d) %s%s" % (i + 1, st.raw, tag))
        L.append("%s %s" % (t("气相:", "gas:"), self.gas_species or none))
        L.append("%s %s" % (t("液相:", "liquid:"), self.liq_species or none))
        L.append("%s %s" % (t("表面:", "surface:"), self.surface_species or none))
        decl = [n for n, s in self.species.items() if s.declared]
        if decl:
            L.append("%s %s" % (t("前缀声明的物种:", "prefix-declared species:"),
                                sorted(decl)))
        if self.exempt_balance:
            L.append(t("⚠ 豁免元素配平: %s —— 这些步的元素配平闸已失效"
                       % self.exempt_balance,
                       "WARNING element balance waived for: %s - the element-balance "
                       "gate is disabled for those steps" % self.exempt_balance))
        return "\n".join(L)


# ══════════════════════════════════ 四道闸 ══════════════════════════════════
def _check_balance(step: Step, species: Dict[str, Species], exempt=()):
    """闸②：逐元素 + 总电荷配平。**抓笔误最有效的一道**。"""
    exempt = set(exempt)

    def side_counts(d):
        el: Dict[str, float] = {}
        ch = 0.0
        for n, c in d.items():
            sp = species[n]
            if sp.phase in ("site",):
                continue                        # 空位不参与元素配平
            if sp.phase == "electron":
                ch += c * sp.charge
                continue
            # 豁免按**化学式**匹配，不是按规范名 —— 否则写 exempt_balance=['Nafion']
            # 只豁免了气相的 'Nafion'，吸附态 'Nafion*' 照样被 _parse_formula 炸掉
            if n in exempt or sp.formula in exempt:
                ch += c * sp.charge      # 豁免元素配平，但电荷仍然要配
                continue
            if sp.composition is not None:
                cnt, q = dict(sp.composition), sp.charge
            else:
                cnt, q = _parse_formula(sp.formula) if sp.formula else ({}, 0)
            for e, k in cnt.items():
                el[e] = el.get(e, 0.0) + c * k
            ch += c * (sp.charge or q)
        return el, ch

    lel, lch = side_counts(step.reactants)
    rel, rch = side_counts(step.products)
    bad = []
    for e in set(lel) | set(rel):
        if abs(lel.get(e, 0.0) - rel.get(e, 0.0)) > 1e-9:
            bad.append("元素 %s：左 %g ≠ 右 %g" % (e, lel.get(e, 0.0), rel.get(e, 0.0)))
    if abs(lch - rch) > 1e-9:
        bad.append("电荷：左 %+g ≠ 右 %+g" % (lch, rch))
    if bad:
        raise MechanismError("这一步没配平：%s\n  %s" % (step.raw, "\n  ".join(bad)))


def _check_site_balance(step: Step, species: Dict[str, Species], sites):
    """闸③：每种位点在左右两侧的总数必须相等（吸附态也占一个位点）。"""
    def side_sites(d):
        c: Dict[str, float] = {}
        for n, k in d.items():
            sp = species[n]
            if sp.site is None:
                continue
            c[sp.site] = c.get(sp.site, 0.0) + k
        return c

    l, r = side_sites(step.reactants), side_sites(step.products)
    bad = [("位点 %s：左 %g ≠ 右 %g" % (("默认" if s == "" else s), l.get(s, 0.0), r.get(s, 0.0)))
           for s in sites if abs(l.get(s, 0.0) - r.get(s, 0.0)) > 1e-9]
    if bad:
        raise MechanismError("这一步位点不守恒：%s\n  %s\n"
                             "  提示：吸附态 CO* 也占一个位点；脱附要把 * 写回产物侧。"
                             % (step.raw, "\n  ".join(bad)))


def from_text(text: str, strict: bool = True,
              species_formula: Dict[str, str] = None,
              exempt_balance: Sequence[str] = ()) -> Mechanism:
    """解析机理文本。每行一个基元反应；`%` 或 `!` 或 `#` 之后为注释。

    参数
    ----
    species_formula : {物种名: 化学式}
        给**非元素式名字**提供组成，例如 ``{"Cat": "C6H4O2", "Dye": "C20H12N2O4"}``。
        这些名字必须在机理里用 ``Sp_``/``E_`` 前缀声明过（见 naming.py）。
    exempt_balance : 物种名列表
        明确豁免**元素**配平的物种（电荷仍然配）—— 给"组成确实说不清"的情形，
        例如聚合物链段。**一旦豁免，这一步的元素配平闸就失效**，所以
        `summary()` 会把豁免清单打出来，让它在输出里可见。
    """
    steps: List[Step] = []
    species: Dict[str, Species] = {}
    for lineno, raw in enumerate(text.splitlines(), 1):
        line = re.split(r"[%!]|(?<!\S)#(?!\d|\w)", raw)[0].strip()
        if not line:
            continue
        # 可选的 {参数}
        params: Dict[str, float] = {}
        mp = re.search(r"\{([^}]*)\}\s*$", line)
        if mp:
            for kv in mp.group(1).split(","):
                if not kv.strip():
                    continue
                k, _, v = kv.partition(":")
                if not _:
                    raise MechanismError("第 %d 行参数要写成 {键: 值}：'%s'" % (lineno, kv))
                params[k.strip()] = float(v)
            line = line[: mp.start()].strip()
        # 可逆 / 不可逆
        if "<->" in line:
            lhs, rhs, rev = (*line.split("<->", 1), True)
        elif "->" in line:
            lhs, rhs, rev = (*line.split("->", 1), False)
        else:
            raise MechanismError("第 %d 行没有 '<->' 或 '->'：%s" % (lineno, raw))
        try:
            rct, r1 = _parse_side(lhs, species_formula)
            prd, r2 = _parse_side(rhs, species_formula)
        except MechanismError as e:
            raise MechanismError("第 %d 行：%s\n  原文：%s" % (lineno, e, raw)) from None
        species.update(r1)
        species.update(r2)
        st = Step(raw=line.strip(), reactants=rct, products=prd,
                  reversible=rev, params=params)
        # ★ nu_elec 自动数出来，不让用户填（最易写错符号的参数）
        ne = prd.get("e-", 0.0) - rct.get("e-", 0.0)
        if ne:
            st.params.setdefault("nu_elec", ne)
        steps.append(st)

    if not steps:
        raise MechanismError("没解析出任何基元反应")

    # 闸①：撞名（物种名与位点名）—— **带前缀声明过的名字豁免**
    user_names = [s.formula for s in species.values() if s.formula and not s.declared] + \
                 [s for s in {sp.site for sp in species.values() if sp.site}]
    check_user_names(sorted(set(user_names)), kind="物种/位点")

    # 闸①b：声明了却组成不明、又没豁免 ⇒ 报错（不静默放弃元素配平）
    exempt = set(exempt_balance)
    # ★ 按**化学式**去重报：`Cat` 与 `Cat*` 是同一个式子，用户只需给一次组成
    #   （初版按规范名报，于是让用户给 'Cat*' 的组成 —— 多余且会把人带偏）
    unknown = sorted({s.formula for s in species.values()
                      if s.declared and s.composition is None
                      and s.formula and s.formula not in exempt})
    if unknown:
        raise MechanismError(
            "这些物种的元素组成不明：%s\n"
            "它们用前缀声明过，所以名字本身是合法的，但**元素配平这道闸需要组成**。\n"
            "两种给法（选一种）：\n"
            "  ① 给组成： from_text(txt, species_formula={%s})\n"
            "  ② 明确豁免元素配平（电荷仍配）： from_text(txt, exempt_balance=%s)\n"
            "  —— 不替你猜，因为配平是本包最有效的查错手段。"
            % (unknown, ", ".join("%r: '…'" % n for n in unknown), unknown))

    sites = sorted({s.site for s in species.values() if s.site is not None})
    mech = Mechanism(steps, species)
    mech.exempt_balance = sorted(exempt)
    if strict:
        for st in steps:
            _check_balance(st, species, exempt)      # 闸②
            _check_site_balance(st, species, sites)   # 闸③
    return mech
