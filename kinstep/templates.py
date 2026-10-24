# -*- coding: utf-8 -*-
# KinSTEP - a modality-agnostic microkinetic solver.
# Copyright (C) 2026 The KinSTEP authors.
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.  It is distributed WITHOUT ANY WARRANTY; see the GNU General Public
# License (LICENSE) for details.
r"""[English] Offline template library: ready-made mechanisms for common systems. A template
carries structure only, never energy values.

(The commentary below is in Chinese: it records why each choice was made
and which measurement forced it.)

离线模板库（T14）：给常见体系一份**现成的机理骨架**。

红线：**只给「结构」和「问题」，不给「数字」**
----------------------------------------------
这是台账里定过的（`provenance.py` 的设计依据）：

  · **机理是公开知识** —— 它写错了肉眼能看出来，而且本来就发表在文献里；
  · **能量是私有数据** —— 它填错了**只会得到一条看着很合理的火山曲线**。

所以模板给的是：机理文本、哪些步要哪些参数、典型工况、以及**出处**；
**绝不给任何能量数值**。每个能量槽位都是 `Val.required()` 并带一句
「你需要算/查什么」的提示。

⚠ 这条红线由 `tests/t14_templates.py` **机器检查**：
   任何模板里出现能量数值，测试直接失败。不靠自觉。

与"别再填复杂输入脚本"的关系
----------------------------
    kinstep templates                 # 看有哪些体系
    kinstep new --template her        # 生成骨架，机理已填好
    kinstep check xxx.kinstep.toml    # 它会明确告诉你还缺哪几个能量

—— 用户要做的只剩「把自己体系的 DFT 能量填进去」，这正是**只有他能做**的那部分。

⚠ 出处说明
----------
下列文献是经典机理的常见出处，**引用前请自行核实**（本包不替你背书）。
给了 DOI 的是本项目在开发中实地核过的；没给的只给作者/期刊/年份。
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence

from .mech import from_text as _parse


@dataclass
class Template:
    """一个体系的机理骨架。**不含任何能量数值。**"""
    key: str
    title: str
    modality: str                 # thermal / electro / photo / sono
    mechanism: str                # 机理文本（自由格式，照纸上写）
    conditions: Dict[str, object] # 典型工况（温度/压力等，**不是能量**）
    drivers: Dict[str, Dict] = field(default_factory=dict)
    reference: str = ""
    doi: str = ""                 # 只在本项目实地核过时才填
    notes: List[str] = field(default_factory=list)
    energy_hints: Dict[int, str] = field(default_factory=dict)  # 步号 → 要算什么

    def parse(self):
        """解析机理（四道闸都会跑）——模板自己必须是合法的。"""
        return _parse(self.mechanism)

    def describe(self) -> str:
        from .i18n import t
        m = self.parse()
        L = [t("%s —— %s（%s）", "%s - %s (%s)")
             % (self.key, _en(TITLE_EN, self.key, self.title), self.modality)]
        L.append("  " + t("机理 %d 步，%d 个物种" % (m.n_step, len(m.species_order)),
                          "mechanism: %d steps, %d species"
                          % (m.n_step, len(m.species_order))))
        for i, st in enumerate(m.steps):
            hint = hint_localized(self.key, i + 1,
                                  self.energy_hints.get(i + 1, ""))
            L.append("    (%d) %-38s %s" % (i + 1, st.raw, hint))
        L.append("  " + t("典型工况：%s" % self.conditions,
                          "typical conditions: %s" % self.conditions))
        if self.drivers:
            L.append("  " + t("需要的 driver：%s" % list(self.drivers),
                              "driver(s) required: %s" % list(self.drivers)))
        L.append("  " + t("出处：", "source: ")
                 + _en(REF_EN, self.key, self.reference)
                 + (("  DOI %s" % self.doi) if self.doi else ""))
        for n in self.notes_localized():
            L.append("  · " + n)
        return "\n".join(L)

    def notes_localized(self) -> List[str]:
        """按当前语言返回 notes；英文表里缺这一项就退回中文。"""
        from .i18n import lang
        if lang() == "en":
            return NOTES_EN.get(self.key) or self.notes
        return self.notes

    def title_localized(self) -> str:
        return _en(TITLE_EN, self.key, self.title)

    def reference_localized(self) -> str:
        return _en(REF_EN, self.key, self.reference)

    @staticmethod
    def _k(key) -> str:
        """TOML 键：裸键只允许 A-Za-z0-9_-，其余必须加引号。

        ★ 这是一个真 bug 的修法（T14 测试抓到）：物种名里常有 `+` `-` `(` `)`
          （`H+`、`OH-`、`H2O(liq)`），直接写成裸键会让生成的文件**根本读不了** ——
          OER / ORR / CO-stripping 三个模板都中招。
        """
        import re as _re
        return key if _re.fullmatch(r"[A-Za-z0-9_-]+", str(key)) else '"%s"' % key

    def to_toml(self, overrides: Dict[str, object] = None,
                override_note: str = None) -> str:
        """生成 `.kinstep.toml` 骨架。**能量槽位留空并写明要填什么。**

        `overrides` **覆盖** `[conditions]` 里的同名项，而不是在前面**追加** ——
        追加会产生重复键，TOML 直接报 "Cannot overwrite a value"。
        这是 T15 的 `ask --write` 当场踩到的真 bug。
        被覆盖的项会标注来源，让用户一眼看出哪些是他自己说的。
        """
        # ★ 这份骨架文件**本身也要跟着语言走**。
        #   T16 第一版只扫 stdout，于是 `KINSTEP_LANG=en` 下
        #   `kinstep new` 回显是英文、而写出来的 .toml 整个是中文 ——
        #   审稿人拿到的恰恰是这份文件。现在 ④ 会连文件一起扫。
        from .i18n import t
        if override_note is None:
            override_note = t("自你的描述", "from your description")
        m = self.parse()
        ov = dict(overrides or {})
        _title = self.title_localized()
        _ref = self.reference_localized()
        _doi = ("  DOI %s" % self.doi) if self.doi else ""
        out = [t("# %s —— %s" % (self.key, _title),
                 "# %s - %s" % (self.key, _title)),
               t("# 出处：%s%s" % (_ref, _doi),
                 "# source: %s%s" % (_ref, _doi)),
               t("# ⚠ 本模板**不含能量数值** —— 机理是公开知识，能量是你的体系的私有数据。",
                 "# ! This template carries NO energy values - the mechanism is public "
                 "knowledge, the energies are private data for your own system."),
               t("#    下面每个 [[steps]] 都要你填 Ga_fwd 与 dG_rxn（eV）。",
                 "#    every [[steps]] block below needs Ga_fwd and dG_rxn (eV) "
                 "from you."),
               "", "[mechanism]", "text = '''"]
        out += [ln for ln in self.mechanism.strip().splitlines()]
        out += ["'''", "", "[conditions]"]
        merged = dict(self.conditions)
        merged.update(ov)                      # ★ 覆盖，不是追加
        for k, v in merged.items():
            if k in ov:
                out.append(t("# ← %s" % override_note,
                             "# <- %s" % override_note))
            if isinstance(v, dict):
                out.append("%s = { %s }"
                           % (self._k(k), ", ".join("%s = %r" % (self._k(kk), vv)
                                                    for kk, vv in v.items())))
            else:
                out.append("%s = %r" % (self._k(k), v))
        for dname, dcfg in self.drivers.items():
            out += ["", "[drivers.%s]" % dname]
            for k, v in dcfg.items():
                out.append("%s = %r" % (self._k(k), v))
        out.append("")
        for i, st in enumerate(m.steps):
            hint = hint_localized(self.key, i + 1,
                                  self.energy_hints.get(i + 1, ""))
            out += ["[[steps]]",
                    "# (%d) %s" % (i + 1, st.raw)]
            if hint:
                out.append(t("#     要算什么：%s" % hint,
                             "#     what to compute: %s" % hint))
            out += [t("# Ga_fwd = ?   # 正向活化自由能 (eV) —— 填完去掉注释",
                      "# Ga_fwd = ?   # forward activation free energy (eV) - "
                      "uncomment once filled in"),
                    t("# dG_rxn = ?   # 反应自由能 (eV)",
                      "# dG_rxn = ?   # reaction free energy (eV)"),
                    ""]
        out += ["[solver]", 'precision = "auto"', ""]
        return "\n".join(out)


def _T(**kw):
    return Template(**kw)


# ════════════════════════ 英文版文案（KINSTEP_LANG=en 时用）════════════════════════
# ★ 为什么集中放一张表，而不是给每个模板加 title_en / notes_en 字段：
#   11 个模板的英文分散在 160 行里，改了中文忘了英文就看不出来；
#   集中成表之后，中英两版可以并排核对，缺哪个 key 一眼就能发现。
#   `_en()` 在缺条目时**退回中文**（宁可显示中文，也不显示空字符串）。
TITLE_EN = {
    "ammonia": "ammonia synthesis (Haber-Bosch)",
    "co_oxidation": "CO oxidation (Langmuir-Hinshelwood)",
    "wgs": "water-gas shift (carboxyl route)",
    "methanation": "CO methanation",
    "her": "hydrogen evolution, HER (Volmer-Heyrovsky-Tafel)",
    "oer": "oxygen evolution, OER (four-electron route)",
    "orr": "oxygen reduction, ORR (associative route)",
    "co_stripping": "CO stripping (CV / potential step)",
    "photo_water": "photocatalytic water splitting",
    "photo_methanol": "photocatalytic methanol reforming to H2",
    "sono_radical": "sonochemical water splitting to OH radicals (cavitation hot spot)",
}

NOTES_EN = {
    "ammonia": [
        "the usual descriptor is E_N (binding energy of N); the scaling and BEP "
        "coefficients must be fitted for your own system",
        "this system was used for the 401-point digit-by-digit regression of this "
        "package (test T7)"],
    "co_oxidation": [
        "two-dimensional volcano plots usually take (E_CO, E_O) as descriptors"],
    "wgs": [
        "there is also a redox route, CO* + O* -> CO2 + 2*, which can be added as a "
        "parallel path"],
    "methanation": [
        "if your system follows H-assisted dissociation, replace the first step by "
        "CO* + H* <-> CHO* + *"],
    "her": [
        "H+ appears explicitly in the mechanism, so the potential MUST be on the SHE "
        "scale - using RHE would count pH twice",
        "the Tafel slope identifies the rate-determining step: Volmer ~120, "
        "Heyrovsky ~40, Tafel ~30 mV/dec"],
    "oer": [
        "H+ is explicit, so use the SHE scale; to use RHE, remove H+ from the mechanism",
        "the theoretical overpotential is usually taken as max(dG_i) - 1.23 eV"],
    "orr": [
        "the usual descriptor is dG_OH; the volcano peak lies near "
        "dG_OH = 0.1-0.2 eV relative to Pt",
        "the CatMAP ORR tutorial reproduces Hansen et al. J. Phys. Chem. C 2014, "
        "118, 6706",
        "H+ is explicit, so use the SHE scale"],
    "co_stripping": [
        "this system was used to regression-test the transient core of this package: "
        "the ODE right-hand side matches to 2.4e-15 (test T5)",
        "suitable for both the CVSweep and PotentialStep protocols"],
    "photo_water": [
        "the electron and hole consumption rates must balance (this package checks it) "
        "- a wrong nu_* does not raise an error, it just gives a wrong curve",
        "the default qfl_frac_n = 0.5 only holds for intrinsic or high-injection "
        "conditions; for a doped system you must supply your own value",
        "the high-intensity power-law exponent is beta*nu*f; the classical value of "
        "0.5 presupposes that the rate-determining step transfers one carrier"],
    "photo_methanol": [
        "if the rate-determining step transfers an electron AND a hole, the driving "
        "force is the full qfl splitting, so there is NO half-order regime"],
    "sono_radical": [
        "reported hot-spot temperatures differ widely with the method: 5000-20000 K "
        "from single-bubble sonoluminescence, 750-6000 K for multi-bubble fields",
        "calibrate the hot-spot fraction against a measured dosimeter yield "
        "(calibrate_fraction_from_yield) rather than guessing it",
        "two-site model: gas phase inside the bubble plus the interfacial liquid "
        "shell; the bulk is essentially unreactive"],
}

REF_EN = {
    "ammonia": "Ertl, G. Catal. Rev. Sci. Eng. 1980, 21, 201 "
               "(N2 dissociation as the rate-determining step)",
    "co_oxidation": "Langmuir-Hinshelwood mechanism; for descriptor selection see "
                    "Falsig et al. Angew. Chem. Int. Ed. 2008, 47, 4835",
    "wgs": "carboxyl route, see Gokhale, Dumesic, Mavrikakis, "
           "J. Am. Chem. Soc. 2008, 130, 1402",
    "methanation": "direct CO dissociation route; descriptor analysis in "
                   "Andersson et al. J. Catal. 2008, 255, 6",
    "her": "Volmer-Heyrovsky-Tafel mechanism; for the dG_H* volcano see "
           "Norskov et al. J. Electrochem. Soc. 2005, 152, J23",
    "oer": "Rossmeisl, Qu, Zhu, Kroes, Norskov, J. Electroanal. Chem. 2007, 607, 83",
    "orr": "Norskov et al. J. Phys. Chem. B 2004, 108, 17886 "
           "(the original computational hydrogen electrode paper)",
    "co_stripping": "Koper, Jansen, van Santen, Lukkien, Hilbers, "
                    "J. Chem. Phys. 1998, 109, 6051",
    "photo_water": "quasi-Fermi-level splitting as the driving force: see the "
                   "ACS Catal. 2025 perspective 10.1021/acscatal.5c07014",
    "photo_methanol": "photocatalytic H2 evolution with methanol as hole scavenger; "
                      "please verify the mechanism partitioning before citing",
    "sono_radical": "hot-spot theory: Suslick et al. Phil. Trans. R. Soc. A 1999, "
                    "357, 335; Flint & Suslick, Science 1991, 253, 1397",
}


HINTS_EN = {
    "ammonia": {1: "N2 dissociative adsorption - usually rate-determining, "
                   "the barrier matters most",
                2: "H2 dissociative adsorption - quasi-equilibrated on most metals",
                5: "NH3 desorption - remember the gas-phase entropy "
                   "(see the thermo module)"},
    "co_oxidation": {2: "O2 dissociative adsorption - quasi-equilibrated on noble "
                        "metals, possibly rate-determining on transition metals",
                     3: "surface coupling - the two adsorbates must be adjacent; "
                        "consider n_neighbor"},
    "wgs": {3: "COOH formation - this is what separates the carboxyl route from "
               "the redox route",
            4: "COOH dissociation - usually rate-determining"},
    "methanation": {2: "CO dissociation - competes with the H-assisted route and is "
                       "usually rate-determining"},
    "her": {1: "Volmer step - the key quantity is dG_H* (hydrogen adsorption free "
               "energy)",
            2: "Heyrovsky step",
            3: "Tafel step (chemical recombination, no electron transfer)"},
    "oer": {1: "*OH formation",
            2: "*O formation",
            3: "*OOH formation - the dG(OOH)-dG(OH) = 3.2 eV scaling relation is the "
               "origin of the OER overpotential",
            4: "O2 release"},
    "orr": {1: "*OOH formation - limiting step on the weak-binding side",
            4: "*OH reductive desorption - limiting step on the strong-binding side"},
    "co_stripping": {1: "OH adsorption - usually supplied as a rate constant "
                        "(mode='given')",
                     2: "CO + OH coupling - n_neighbor is the coordination number "
                        "(Koper uses 4)"},
    "photo_water": {2: "hole oxidation and desorption - the classical half-order "
                       "light dependence requires that only ONE carrier is transferred",
                    3: "hydrogen evolution - two electrons"},
    "photo_methanol": {2: "carrier-mediated dehydrogenation - if it consumes one "
                          "electron AND one hole, the light dependence is first order"},
    "sono_radical": {2: "water dissociation - the barrier is far too high at ambient "
                        "temperature; this is exactly what the hot spot is for"},
}


def _en(table, key, zh):
    """英文模式下取表里的英文；缺条目就退回中文（**不返回空串**）。"""
    from .i18n import t
    return t(zh, table.get(key) or zh)


def hint_localized(key, step, zh):
    """逐步能量提示的本地化；英文表里缺这一条就退回中文。"""
    from .i18n import lang
    if lang() == "en":
        return (HINTS_EN.get(key) or {}).get(step) or zh
    return zh


LIBRARY: Dict[str, Template] = {t.key: t for t in [
    # ───────────────────────── 热催化 ─────────────────────────
    _T(key="ammonia", title="合成氨（Haber–Bosch）", modality="thermal",
       mechanism="""
N2 + 2* <-> 2N*
H2 + 2* <-> 2H*
N* + H* <-> NH* + *
NH* + H* <-> NH2* + *
NH2* + H* <-> NH3 + 2*
""",
       conditions=dict(T_K=673.0, p_gas={"N2": 25.0, "H2": 75.0, "NH3": 1.0}),
       reference="Ertl, G. Catal. Rev. Sci. Eng. 1980, 21, 201（N2 解离为控速步）",
       energy_hints={1: "N2 解离吸附 —— 通常是控速步，势垒最要紧",
                     2: "H2 解离吸附 —— 多数金属上准平衡",
                     5: "NH3 脱附 —— 注意算上气相熵（见 thermo 模块）"},
       notes=["描述符常取 E_N（N 的吸附能）；scaling 与 BEP 系数要自己拟合",
              "本项目用它做了 401 点逐数字回归（见台账 T7）"]),

    _T(key="co_oxidation", title="CO 氧化（Langmuir–Hinshelwood）", modality="thermal",
       mechanism="""
CO + * <-> CO*
O2 + 2* <-> 2O*
CO* + O* -> CO2 + 2*
""",
       conditions=dict(T_K=500.0, p_gas={"CO": 0.1, "O2": 0.1, "CO2": 1e-6}),
       reference="Langmuir–Hinshelwood 机理；描述符筛选见 Falsig et al. "
                 "Angew. Chem. Int. Ed. 2008, 47, 4835",
       doi="10.1002/anie.200801479",
       energy_hints={2: "O2 解离吸附 —— 贵金属上常为准平衡，过渡金属上可能控速",
                     3: "表面偶联 —— 两个吸附态要相邻，可考虑 n_neighbor"},
       notes=["二维火山图常取 (E_CO, E_O) 作描述符"]),

    _T(key="wgs", title="水煤气变换（羧基路径）", modality="thermal",
       mechanism="""
CO + * <-> CO*
H2O + 2* <-> OH* + H*
CO* + OH* -> COOH* + *
COOH* + * -> CO2 + H* + *
2H* -> H2 + 2*
""",
       conditions=dict(T_K=573.0, p_gas={"CO": 0.1, "H2O": 0.1, "CO2": 1e-3, "H2": 1e-3}),
       reference="羧基（carboxyl）路径，见 Gokhale, Dumesic, Mavrikakis, "
                 "J. Am. Chem. Soc. 2008, 130, 1402",
       energy_hints={3: "COOH 生成 —— 羧基路径与氧化还原路径的分水岭",
                     4: "COOH 解离 —— 常为控速步"},
       notes=["另有氧化还原（redox）路径：CO* + O* -> CO2 + 2*，可并入做并行路径"]),

    _T(key="methanation", title="CO 甲烷化", modality="thermal",
       mechanism="""
CO + * <-> CO*
CO* + * -> C* + O*
H2 + 2* <-> 2H*
C* + H* <-> CH* + *
CH* + H* <-> CH2* + *
CH2* + H* <-> CH3* + *
CH3* + H* <-> CH4 + 2*
O* + H* <-> OH* + *
OH* + H* <-> H2O + 2*
""",
       conditions=dict(T_K=573.0, p_gas={"CO": 0.1, "H2": 0.3, "CH4": 1e-4, "H2O": 1e-4}),
       reference="CO 直接解离路径；描述符分析见 Andersson et al. J. Catal. 2008, 255, 6",
       energy_hints={2: "CO 解离 —— 与 H 辅助解离路径竞争，常为控速步"},
       notes=["若体系走 H 辅助解离，要改成 CO* + H* <-> CHO* + * 的路径"]),

    # ───────────────────────── 电催化 ─────────────────────────
    _T(key="her", title="析氢反应 HER（Volmer–Heyrovsky–Tafel）", modality="electro",
       mechanism="""
H+ + e- + * <-> H*            {beta_ct: 0.5}
H* + H+ + e- -> H2 + *        {beta_ct: 0.5}
2H* -> H2 + 2*
""",
       conditions=dict(T_K=298.15, c_liq={"H+": 1.0}, p_gas={"H2": 1.0}),
       drivers={"electro": dict(U_elec=0.0, U_ref="SHE", pH_val=0.0, gamma_site=2.2e-5)},
       reference="Volmer–Heyrovsky–Tafel 机理；描述符 ΔG_H* 的火山关系见 "
                 "Nørskov et al. J. Electrochem. Soc. 2005, 152, J23",
       energy_hints={1: "Volmer 步 —— 关键量是 ΔG_H*（氢吸附自由能）",
                     2: "Heyrovsky 步", 3: "Tafel 步（化学复合，无电子转移）"},
       notes=["★ 机理里显式写了 H+，所以电位标度必须用 SHE —— 用 RHE 会**重复计 pH**",
              "Tafel 斜率可区分控速步：Volmer ~120、Heyrovsky ~40、Tafel ~30 mV/dec"]),

    _T(key="oer", title="析氧反应 OER（四电子路径）", modality="electro",
       mechanism="""
H2O(l) + * <-> OH* + H+ + e-      {beta_ct: 0.5}
OH* <-> O* + H+ + e-              {beta_ct: 0.5}
O* + H2O(l) <-> OOH* + H+ + e-    {beta_ct: 0.5}
OOH* -> O2 + * + H+ + e-          {beta_ct: 0.5}
""",
       conditions=dict(T_K=298.15, c_liq={"H2O": 1.0, "H+": 1.0}, p_gas={"O2": 1e-3}),
       drivers={"electro": dict(U_elec=1.23, U_ref="SHE", gamma_site=2.2e-5)},
       reference="Rossmeisl, Qu, Zhu, Kroes, Nørskov, J. Electroanal. Chem. 2007, 607, 83",
       energy_hints={1: "*OH 生成", 2: "*O 生成", 3: "*OOH 生成 —— ΔG(OOH)−ΔG(OH)≈3.2 eV 的 scaling 是 OER 过电位的根源",
                     4: "O2 释放"},
       notes=["★ 显式写了 H+ ⇒ 用 SHE 标度；要用 RHE 请把 H+ 从机理里去掉",
              "理论过电位常由 max(ΔG_i) − 1.23 eV 给出"]),

    _T(key="orr", title="氧还原反应 ORR（缔合路径）", modality="electro",
       mechanism="""
O2 + * + H+ + e- <-> OOH*          {beta_ct: 0.5}
OOH* + H+ + e- -> O* + H2O(l)      {beta_ct: 0.5}
O* + H+ + e- <-> OH*               {beta_ct: 0.5}
OH* + H+ + e- -> H2O(l) + *        {beta_ct: 0.5}
""",
       conditions=dict(T_K=298.15, p_gas={"O2": 0.21}, c_liq={"H2O": 1.0, "H+": 1.0}),
       drivers={"electro": dict(U_elec=0.9, U_ref="SHE", gamma_site=2.2e-5)},
       reference="Nørskov et al. J. Phys. Chem. B 2004, 108, 17886（CHE 方法的原始文献）",
       energy_hints={1: "*OOH 生成 —— 弱吸附端的限制步",
                     4: "*OH 还原脱附 —— 强吸附端的限制步"},
       notes=["描述符常取 ΔG_OH；火山顶在 ΔG_OH ≈ 0.1–0.2 eV（相对 Pt）",
              "CatMAP 的 ORR 教程复现的是 Hansen et al. J. Phys. Chem. C 2014, 118, 6706",
              "★ 显式写了 H+ ⇒ 用 SHE 标度"]),

    _T(key="co_stripping", title="CO 剥离（CV / 电势阶跃）", modality="electro",
       mechanism="""
H2O(l) + * <-> OH* + H+ + e-            {beta_ct: 0.5}
CO* + OH* -> CO2 + H+ + e- + 2*         {beta_ct: 0.5, n_neighbor: 4.0}
""",
       conditions=dict(T_K=300.0, c_liq={"H2O": 1.0, "H+": 1.0}, p_gas={"CO2": 0.0}),
       drivers={"electro": dict(U_elec=0.3, U_ref="SHE", gamma_site=2.1919e-5)},
       reference="Koper, Jansen, van Santen, Lukkien, Hilbers, J. Chem. Phys. 1998, 109, 6051",
       doi="10.1063/1.477228",
       energy_hints={1: "OH 吸附 —— 常以 k⁰ 形式给（mode='given'）",
                     2: "CO + OH 偶联 —— n_neighbor = 近邻配位数（Koper 取 4）"},
       notes=["★ 本项目用它做了瞬态内核的回归：ODE 右端对到 2.4e-15（见台账 T5）",
              "适合用 CVSweep / PotentialStep 两种协议"]),

    # ───────────────────────── 光催化 ─────────────────────────
    _T(key="photo_water", title="光催化水分解", modality="photo",
       mechanism="""
H2O(l) + 2* <-> OH* + H*
OH* -> OH(l) + *            {nu_hole_vb: 1}
2H* -> H2 + 2*              {nu_elec_cb: 2}
""",
       conditions=dict(T_K=300.0, c_liq={"H2O": 1.0, "OH": 1e-9}, p_gas={"H2": 1e-6}),
       drivers={"photo": dict(photo_carrier_mode="steady_carrier", photon_flux=1e21,
                              absorptance=0.5, gamma_site_photo=1.5e19,
                              photon_energy_eV=3.2, qfl_frac_n=0.5,
                              recomb_rate_const=1e-10)},
       reference="准费米能级分裂作为驱动力：见 ACS Catal. 2025 perspective "
                 "10.1021/acscatal.5c07014；以及 Nat. Commun. 2025, 16 (表面空穴 QFL 测定)",
       doi="10.1021/acscatal.5c07014",
       energy_hints={2: "空穴氧化脱附 —— **只转移一个载流子**时才有经典半级光强依赖",
                     3: "析氢 —— 两个电子"},
       notes=["★ 电子与空穴的消耗速率必须相等（本包有自检）—— 配错 nu_* 不会报错只会给错曲线",
              "★ qfl_frac_n 默认 0.5 只对本征/高注入成立；掺杂体系要自己给",
              "高光强幂指数 = β·ν·f，经典 0.5 隐含「控速步只转移一个载流子」"]),

    _T(key="photo_methanol", title="光催化甲醇重整产氢", modality="photo",
       mechanism="""
CH3OH + * <-> CH3OH*
CH3OH* -> CH2O* + H2       {nu_hole_vb: 1, nu_elec_cb: 1}
CH2O* -> CH2O + *
""",
       conditions=dict(T_K=300.0, p_gas={"CH3OH": 0.1, "CH2O": 1e-6, "H2": 1e-6}),
       drivers={"photo": dict(photo_carrier_mode="steady_carrier", photon_flux=1e21,
                              absorptance=0.5, gamma_site_photo=1.5e19,
                              photon_energy_eV=3.2, recomb_rate_const=1e-10)},
       reference="甲醇作空穴牺牲剂的光催化产氢；机理分区见相关综述，引用前请核实",
       energy_hints={2: "载流子介导的脱氢 —— 同时吃 1 电子 + 1 空穴时光强依赖为一级"},
       notes=["⚠ 若控速步同时转移电子与空穴，驱动力 = 完整的 ΔE_qfl ⇒ **没有半级区间**"]),

    # ───────────────────────── 超声催化 ─────────────────────────
    _T(key="sono_radical", title="超声水解离产 •OH（空化热点）", modality="sono",
       mechanism="""
H2O(l) + * <-> H2O*
H2O* + * -> OH* + H*
OH* -> OH(l) + *
2H* -> H2 + 2*
""",
       conditions=dict(T_K=300.0, c_liq={"H2O": 1.0, "OH": 1e-9}, p_gas={"H2": 1e-6},
                       hotspot_T_K=5000.0, hotspot_fraction=1e-6,
                       interfacial_T_K=1900.0, interfacial_fraction=1e-4),
       drivers={"sono": {}},
       reference="热点理论：Suslick 等, Phil. Trans. R. Soc. A 1999, 357, 335；"
                 "Flint & Suslick, Science 1991, 253, 1397",
       energy_hints={2: "水的解离 —— 常温下势垒很高，正是靠热点才走得动"},
       notes=["★ 热点温度随测法差很多：SBSL 反推 5000–20000 K，多气泡 750–6000 K",
              "★ 份额建议用实测剂量计反推（calibrate_fraction_from_yield），不要猜",
              "两点模型：气泡内气相 + 界面液壳；本体基本不反应"]),
]}


def list_templates(modality: str = None) -> List[str]:
    """列出模板名；给了 `modality` 就只列该模态的。"""
    return sorted(k for k, t in LIBRARY.items()
                  if modality is None or t.modality == modality)


def get(key: str) -> Template:
    """按名字取模板；取不到时把现有的名字一并列出来。"""
    if key not in LIBRARY:
        raise KeyError("没有模板 %r。现有：%s" % (key, list_templates()))
    return LIBRARY[key]


def summary() -> str:
    """按模态分组的一览表（CLI 的 `kinstep templates` 用它）。"""
    from .i18n import t
    L = [t("离线模板库：%d 个体系（**只给机理与问题，不给能量数值**）" % len(LIBRARY),
           "offline template library: %d systems "
           "(mechanism and questions only, never energy values)" % len(LIBRARY))]
    for mod in ("thermal", "electro", "photo", "sono"):
        ks = list_templates(mod)
        if not ks:
            continue
        L.append("  " + t("【%s】" % mod, "[%s]" % mod))
        for k in ks:
            L.append("    %-16s %s" % (k, LIBRARY[k].title_localized()))
    L.append("\n  " + t("用法： kinstep new --template <名字>   然后 kinstep check 会告诉你还缺什么",
                        "usage: kinstep new --template <name>, then kinstep check "
                        "tells you what is still missing"))
    return "\n".join(L)
