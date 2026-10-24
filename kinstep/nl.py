# -*- coding: utf-8 -*-
# KinSTEP - a modality-agnostic microkinetic solver.
# Copyright (C) 2026 The KinSTEP authors.
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.  It is distributed WITHOUT ANY WARRANTY; see the GNU General Public
# License (LICENSE) for details.
r"""[English] Natural-language layer: one sentence becomes a structure plus a list of
questions. Stated quantities are extracted and tagged; energies always become
questions, never assumptions.

(The commentary below is in Chinese: it records why each choice was made
and which measurement forced it.)

自然语言层（T15）：把一句话变成**结构 + 问题**，绝不变成数字。

分三层，风险递增，默认只开前两层
--------------------------------
    ① 离线模板库（T14）      零 LLM，机理来自文献，**无编数字风险**
    ② **规则翻译（本模块）** 零 LLM，关键词匹配 + 带单位的数字抽取
    ③ LLM 后端               **默认关**，且产出一律打 `llm_suggested` 标签

★★ 红线：**只许产出「结构」和「问题」，不许产出「数字」**
---------------------------------------------------------
更准确地说，数字分两类，处理方式完全不同：

  · **用户自己说出来的数字** —— 「673 K」「pH=0」「−0.1 V」
    ⇒ 可以抽取，标 `user_given`。它本来就是用户的数据，抽取只是省去手填。
  · **系统要去猜的数字** —— 任何能量、任何没说的参数
    ⇒ **一律变成问题**，标 `REQUIRED`，由 `provenance.gate()` 拦住求解。

为什么这条界线必须划得这么死：机理写错了肉眼能看出来，
**能量填错了只会得到一条看着很合理的火山曲线**。
—— 这正是 `provenance.py` 当初的设计依据，本模块只是它的上游入口。

★ 单位换算一律**显式记录**
--------------------------
「400 °C」→ 673.15 K 这种换算必须写进 `conversions`，让用户能复核。
而「100 mW/cm²」→ 光子流**需要光子能量**，没给就是**一个问题**，不是一个假设。
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

from .provenance import Val
from .templates import LIBRARY, Template, get as _get_template


# ══════════════════════════ 关键词 → 模板 ══════════════════════════
#: 中英文都收；**一个模板多组关键词**，匹配到几组就算几分
KEYWORDS: Dict[str, Sequence[str]] = {
    "ammonia":        ("合成氨", "氨合成", "固氮", "ammonia", "haber", "n2 + 3h2", "nh3 synth"),
    "co_oxidation":   ("co 氧化", "co氧化", "一氧化碳氧化", "co oxidation", "co+o2"),
    "wgs":            ("水煤气", "变换反应", "water gas shift", "wgs", "co + h2o"),
    "methanation":    ("甲烷化", "methanation", "sabatier", "co + 3h2"),
    "her":            ("析氢", "产氢", "her", "hydrogen evolution", "volmer", "heyrovsky"),
    "oer":            ("析氧", "产氧", "oer", "oxygen evolution", "水氧化", "oxygen evol"),
    "orr":            ("氧还原", "orr", "oxygen reduction", "燃料电池阴极"),
    "co_stripping":   ("co 剥离", "co剥离", "co stripping", "循环伏安", "cv 扫描"),
    "photo_water":    ("光催化水分解", "光解水", "光催化分解水", "photocatalytic water splitting",
                       "水分解", "water splitting"),
    "photo_methanol": ("光催化甲醇", "甲醇重整", "牺牲剂", "methanol reforming",
                       "photocatalytic hydrogen"),
    "sono_radical":   ("超声", "声化学", "空化", "sonocatal", "sonochem", "cavitation",
                       "ultrasound", "ultrasonic"),
}

#: 模态关键词（用于在没匹配到具体体系时，至少说清楚是哪一类）
MODALITY_HINTS = {
    "electro": ("电催化", "电化学", "电极", "电位", "电势", "overpotential", "electro",
                "vs rhe", "vs she", "tafel"),
    "photo":   ("光催化", "光照", "可见光", "紫外", "光子", "photo", "irradiat", "light"),
    "sono":    ("超声", "声化学", "空化", "sono", "ultrasound", "cavitation"),
    "thermal": ("热催化", "固定床", "thermal", "fixed bed"),
}


@dataclass
class Extracted:
    """一个抽取出来的量：值、单位、原文、以及**换算过程**。"""
    key: str
    value: float
    unit: str
    raw: str
    conversion: str = ""

    def __str__(self):
        from .i18n import t as _t
        c = (_t("（%s）", " (%s)") % self.conversion) if self.conversion else ""
        v = self.unit if self.value != self.value else "%g" % self.value   # nan ⇒ 显示单位
        return _t("%s = %s  ← 原文「%s」%s", "%s = %s  <- from \"%s\"%s") % (
            self.key, v, self.raw, c)


@dataclass
class Request:
    """一句话理解的结果：**结构 + 问题**，没有猜出来的数字。"""
    text: str
    template_key: Optional[str]
    candidates: List[Tuple[str, int]]          # (模板, 命中关键词数)
    modality: Optional[str]
    extracted: List[Extracted] = field(default_factory=list)
    questions: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)

    @property
    def confident(self) -> bool:
        """只有唯一最高分且分数 ≥1 才算有把握。**平手一律不替用户选。**"""
        if not self.candidates:
            return False
        top = self.candidates[0][1]
        return top >= 1 and sum(1 for _, s in self.candidates if s == top) == 1

    def conditions(self) -> Dict[str, object]:
        """抽取到的工况，**全部标 `user_given`** —— 它们本来就是用户说的。

        参比标度这类**非数值**的量（U_ref）直接给字符串，不包成 Val。
        """
        out = {}
        for e in self.extracted:
            if e.value != e.value:              # nan ⇒ 是个标度/标签
                out[e.key] = e.unit
            else:
                out[e.key] = Val(e.value, "user_given", detail="自原文「%s」" % e.raw)
        return out

    def report(self) -> str:
        from .i18n import t as _t
        L = [_t("理解结果（规则翻译，零 LLM）",
                "what was understood (rule-based translation, no LLM)")]
        L.append("  " + _t("原文：", "input: ") + self.text.strip())
        if self.template_key:
            tp = _get_template(self.template_key)
            L.append("  " + _t("匹配体系：", "matched system: ")
                     + "%s —— %s [%s]" % (tp.key, tp.title_localized(), tp.modality))
        else:
            L.append("  ✗ " + _t("没有匹配到具体体系",
                                  "no specific system could be matched"))
        if len(self.candidates) > 1:
            L.append("  " + _t("其它候选：", "other candidates: ")
                     + _t("，", ", ").join("%s(%d)" % c for c in self.candidates[1:4]))
        if self.modality and not self.template_key:
            L.append("  " + _t("但能判断模态：", "but the modality is clear: ")
                     + str(self.modality))
        if self.extracted:
            L.append("  " + _t("从原文抽到的量（全部标记为 user_given）：",
                               "values extracted from your sentence "
                               "(all tagged user_given):"))
            for e in self.extracted:
                L.append("     " + str(e))
        if self.warnings:
            L.append("  " + _t("提醒：", "warnings:"))
            for w in self.warnings:
                L.append("     ⚠ " + w)
        L.append("  ★ " + _t("还需要你提供（**这些本模块不会替你猜**）：",
                             "still needed from you "
                             "(this module will never guess these):"))
        for q in self.questions:
            L.append("     · " + q)
        return "\n".join(L)


# ══════════════════════════ 数字 + 单位抽取 ══════════════════════════
_NUM = r"([-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)"

def _extract(text: str) -> Tuple[List[Extracted], List[str]]:
    from .i18n import t as _t
    out: List[Extracted] = []
    warn: List[str] = []
    t = text.lower().replace("，", ",")

    # 温度：K / ℃ / 摄氏度（★ 换算显式记录）
    for m in re.finditer(_NUM + r"\s*(k\b|开尔文)", t):
        out.append(Extracted("T_K", float(m.group(1)), "K", m.group(0)))
    for m in re.finditer(_NUM + r"\s*(°c|℃|摄氏度|deg\s*c)", t):
        c = float(m.group(1))
        out.append(Extracted("T_K", c + 273.15, "K", m.group(0),
                             "%.6g °C + 273.15" % c))
    # 电位：V / mV，并记下参比
    ref = "RHE" if re.search(r"vs\.?\s*rhe|相对\s*rhe|rhe", t) else (
          "SHE" if re.search(r"vs\.?\s*she|相对\s*she|she|nhe", t) else None)
    for m in re.finditer(_NUM + r"\s*(mv)\b", t):
        v = float(m.group(1))
        out.append(Extracted("U_elec", v * 1e-3, "V", m.group(0), "%g mV ÷ 1000" % v))
    for m in re.finditer(_NUM + r"\s*v\b(?!\s*/)", t):
        out.append(Extracted("U_elec", float(m.group(1)), "V", m.group(0)))
    if ref:
        e = Extracted("U_ref", float("nan"), ref, "vs %s" % ref)
        e.conversion = _t("参比标度（不是数值）", "reference scale (not a number)")
        out.append(e)
    # pH
    for m in re.finditer(r"ph\s*=?\s*" + _NUM, t):
        out.append(Extracted("pH_val", float(m.group(1)), "-", m.group(0)))
    # 压力：bar / atm / MPa / kPa
    for pat, fac, note in ((r"\s*bar\b", 1.0, ""), (r"\s*atm\b", 1.01325, "×1.01325"),
                           (r"\s*mpa\b", 10.0, "×10"), (r"\s*kpa\b", 0.01, "×0.01")):
        for m in re.finditer(_NUM + pat, t):
            v = float(m.group(1))
            out.append(Extracted("p_total_bar", v * fac, "bar", m.group(0),
                                 ("%g %s" % (v, note)) if note else ""))
    # 光强：mW/cm² 或 W/m² —— ★ 换成光子流**需要光子能量**，没给就是一个问题
    for m in re.finditer(_NUM + r"\s*(mw\s*/?\s*cm\^?2|mw·cm-2|mw/cm2)", t):
        out.append(Extracted("irradiance_W_m2", float(m.group(1)) * 10.0, "W/m²",
                             m.group(0), "%g mW/cm² × 10" % float(m.group(1))))
        warn.append("给的是**辐照度**（W/m²），而求解要的是**光子流**（m⁻² s⁻¹）——"
                    "换算需要光子能量（或光谱），这一项本模块不替你假设。")
    # 频率（超声）
    for m in re.finditer(_NUM + r"\s*(khz)", t):
        out.append(Extracted("sono_freq_kHz", float(m.group(1)), "kHz", m.group(0)))
    # 去重：同一个 key 多次出现时保留第一个并提醒
    seen, uniq = {}, []
    for e in out:
        if e.key in seen:
            warn.append("原文里 %s 出现了不止一次（「%s」与「%s」），只取第一个，请核对。"
                        % (e.key, seen[e.key], e.raw))
            continue
        seen[e.key] = e.raw
        uniq.append(e)
    return uniq, warn


def _score(text: str) -> List[Tuple[str, int]]:
    t = text.lower()
    scores = []
    for key, words in KEYWORDS.items():
        n = sum(1 for w in words if w in t)
        if n:
            scores.append((key, n))
    return sorted(scores, key=lambda kv: -kv[1])


def _modality(text: str) -> Optional[str]:
    t = text.lower()
    best, n_best = None, 0
    for mod, words in MODALITY_HINTS.items():
        n = sum(1 for w in words if w in t)
        if n > n_best:
            best, n_best = mod, n
    return best


def understand(text: str) -> Request:
    """★ 主入口：一句话 → `Request`（结构 + 问题）。**零 LLM。**"""
    cands = _score(text)
    extracted, warn = _extract(text)
    mod = _modality(text)
    key = cands[0][0] if cands else None
    req = Request(text=text, template_key=key, candidates=cands, modality=mod,
                  extracted=extracted, warnings=list(warn))
    if cands and len(cands) > 1 and cands[0][1] == cands[1][1]:
        req.template_key = None
        req.warnings.append(
            "有多个体系同分（%s）——**不替你选**。请明确说是哪一个。"
            % "、".join(c[0] for c in cands if c[1] == cands[0][1]))
    req.warnings.extend(_consistency_warnings(req))
    req.questions = _questions(req)
    return req


def _consistency_warnings(req: "Request") -> List[str]:
    """★ 跨「用户说的」与「模板机理」的一致性检查 —— 规则翻译最有价值的一部分。

    这里抓的是**不会报错、只会把结论带偏**的那类矛盾。
    """
    out: List[str] = []
    if not req.template_key:
        return out
    t = _get_template(req.template_key)
    got = {e.key: e for e in req.extracted}

    # ① pH 双重计数：机理里显式写了 H+，而用户说的是 RHE 标度
    has_hp = "H+" in t.mechanism
    uref = got["U_ref"].unit if "U_ref" in got else None
    if has_hp and uref == "RHE":
        out.append("⛔ **pH 可能被重复计算**：%s 的机理里显式写了 H+，"
                   "而你说的是 vs RHE。两者只能选一个 ——\n"
                   "       写 H+ ⇒ 用 SHE 标度；用 RHE ⇒ 把 H+ 从机理里去掉。\n"
                   "       （U_RHE = U_SHE + 0.0592·pH，写了 H+ 又用 RHE 等于把 pH 算两遍）"
                   % t.key)
    if has_hp and uref is None and any(e.key == "U_elec" for e in req.extracted):
        out.append("给了电位却没说参比标度。%s 的机理含 H+ ⇒ **应当用 SHE**；"
                   "若你本意是 RHE，机理要相应改写。" % t.key)

    # ② 模态对不上：说的是光/超声，却匹配到热/电催化的体系（或反之）
    if req.modality and req.modality != t.modality:
        out.append("你的描述像是 **%s**，但匹配到的体系 %s 是 **%s** 的 —— 请确认没选错。"
                   % (req.modality, t.key, t.modality))

    # ③ 温度离该体系典型值太远
    if "T_K" in got:
        T = got["T_K"].value
        T0 = float(t.conditions.get("T_K", 298.15))
        if T > 0 and (T / T0 > 2.5 or T0 / T > 2.5):
            out.append("温度 %g K 与 %s 的典型值 %g K 差得较多 —— 不是不行，"
                       "但请确认机理在这个温度下仍然适用。" % (T, t.key, T0))
    return out


def _questions(req: Request) -> List[str]:
    """★ 本模块真正的产出：**问题清单**。"""
    from .i18n import t as _t
    qs: List[str] = []
    if not req.template_key:
        if req.modality:
            qs.append(_t("具体是哪个反应体系？（模态看起来是 %s；可用 `kinstep templates` "
                         "看现有体系，或直接把机理写出来）",
                         "which reaction system exactly? (the modality looks like %s; "
                         "run `kinstep templates` to see what is available, or just "
                         "write the mechanism out)") % req.modality)
        else:
            qs.append(_t("具体是哪个反应体系？可用 `kinstep templates` 查看，"
                         "或直接把基元反应式写出来（本包就认纸上那种写法）。",
                         "which reaction system exactly? Run `kinstep templates`, or "
                         "just write the elementary steps out - this package reads "
                         "them the way you write them on paper."))
        return qs
    t = _get_template(req.template_key)
    m = t.parse()
    got = {e.key for e in req.extracted}

    # ★ 能量：永远是问题，永远不猜
    qs.append(_t("**%d 步的能量**（每步 Ga_fwd 与 dG_rxn，单位 eV）—— "
                 "这是只有你能提供的部分：",
                 "the energies of all %d steps (Ga_fwd and dG_rxn per step, in eV) "
                 "- this is the part only you can supply:") % m.n_step)
    from .templates import hint_localized as _hl
    for i, st in enumerate(m.steps):
        hint = _hl(t.key, i + 1, t.energy_hints.get(i + 1, ""))
        qs.append("    (%d) %s%s" % (i + 1, st.raw,
                                     (_t("   ← ", "   <- ") + hint) if hint else ""))

    if "T_K" not in got:
        qs.append(_t("温度（K 或 °C）；该体系的典型值是 %g K，但请按你的实验给",
                     "temperature (K or degC); the typical value for this system is "
                     "%g K, but use your own experiment")
                  % t.conditions.get("T_K", 298.15))
    if t.modality == "electro":
        if "U_elec" not in got:
            qs.append(_t("电极电位（V）及参比标度（SHE / RHE）",
                         "electrode potential (V) and its reference scale (SHE / RHE)"))
        elif "U_ref" not in got:
            qs.append(_t("电位的参比标度是 SHE 还是 RHE？—— "
                         "**机理里写了 H+ 就必须用 SHE**，否则 pH 会被重复计算",
                         "is the potential on the SHE or the RHE scale? With H+ "
                         "written explicitly in the mechanism it MUST be SHE, "
                         "otherwise pH is counted twice"))
        if "gamma_site" not in got:
            qs.append(_t("活性位面密度 gamma_site（mol/m²），用来把表面速率换成电流密度",
                         "active-site areal density gamma_site (mol/m2), used to "
                         "convert the surface rate into a current density"))
    if t.modality == "photo":
        qs.append(_t("光子流（m⁻² s⁻¹）或辐照度 + 光子能量；吸收率；活性位面密度",
                     "photon flux (m-2 s-1), or irradiance plus photon energy; "
                     "absorptance; active-site areal density"))
        qs.append(_t("复合系数（双分子 recomb_rate_const 或单分子 recomb_srh_const）——"
                     "两个都不给的话载流子永不复合，那不是物理",
                     "a recombination constant (bimolecular recomb_rate_const or "
                     "monomolecular recomb_srh_const) - with neither, carriers never "
                     "recombine, which is not physical"))
        qs.append(_t("准费米能级分裂的分配比 qfl_frac_n（本征/高注入取 0.5；掺杂体系要另给）",
                     "the partition ratio qfl_frac_n of the quasi-Fermi-level "
                     "splitting (0.5 for intrinsic/high injection; a doped system "
                     "needs its own value)"))
    if t.modality == "sono":
        qs.append(_t("热点温度与份额 —— **份额建议用实测剂量计反推**"
                     "（见 calibrate_fraction_from_yield），不要猜",
                     "the hot-spot temperature and fraction - calibrate the fraction "
                     "against a measured dosimeter yield "
                     "(see calibrate_fraction_from_yield) rather than guessing"))
        qs.append(_t("⚠ 热点温度随测法差异很大：SBSL 反推 5000–20000 K，多气泡 750–6000 K",
                     "WARNING reported hot-spot temperatures differ widely with the "
                     "method: 5000-20000 K from single-bubble sonoluminescence, "
                     "750-6000 K for multi-bubble fields"))
    return qs


# ══════════════════════════ LLM 后端（默认关）══════════════════════════
#: 允许后端返回的字段 —— **白名单**，其余一律丢弃
ALLOWED_BACKEND_FIELDS = {"template_key", "mechanism", "questions", "notes", "modality"}

_BACKEND = None


def set_backend(fn, confirm: bool = False):
    """装一个 LLM 后端。**必须显式 `confirm=True`**。

    ★ 为什么默认关、而且要二次确认
    ------------------------------
    LLM 产出"看着很像"的能量值几乎是零成本的，而**能量填错不会报错，
    只会给出一条看着很合理的火山曲线**。所以：
      ① 默认不装；
      ② 装的时候要显式确认；
      ③ 装上之后，后端的产出**仍然**要过 `sanitize_backend_output()`：
         白名单之外的字段一律丢弃；任何数值一律打 `llm_suggested`，
         而能量类的 `llm_suggested` 会被 `provenance.gate()` 直接拒绝求解。
    三道都过不去，才算安全。
    """
    global _BACKEND
    if fn is not None and not confirm:
        raise PermissionError(
            "装 LLM 后端必须显式确认：set_backend(fn, confirm=True)。\n"
            "  理由：LLM 产出『看着很像』的能量值几乎零成本，而能量填错**不会报错**，"
            "只会给出一条看着很合理的火山曲线。\n"
            "  装上之后它的产出仍然只能是**结构与问题**；任何数值都会被标成 "
            "llm_suggested 并被求解器拒绝。")
    _BACKEND = fn
    return _BACKEND


def backend_enabled() -> bool:
    """LLM 后端装没装。**默认是没装**（见 `set_backend` 的说明）。"""
    return _BACKEND is not None


def sanitize_backend_output(raw: Dict) -> Tuple[Dict, List[str]]:
    """★ 把后端返回的东西过一遍白名单与标签。返回 `(净化后的, 丢弃说明)`。

    规则（都不可绕过）：
      · 白名单之外的字段**丢弃**；
      · `mechanism` 必须能被 `mech.from_text` 解析（四道闸照跑），否则丢弃；
      · **任何数值一律不收**（除非是白名单里本来就不含数值的字段）——
        数值要用户自己给，或走规则翻译从原文抽。
    """
    import numbers
    from . import mech as _mech
    kept, dropped = {}, []
    for k, v in (raw or {}).items():
        if k not in ALLOWED_BACKEND_FIELDS:
            dropped.append("丢弃字段 %r（不在白名单 %s 里）" % (k, sorted(ALLOWED_BACKEND_FIELDS)))
            continue
        if isinstance(v, numbers.Number) and not isinstance(v, bool):
            dropped.append("丢弃 %r = %r：**后端不许产出数值**" % (k, v))
            continue
        if k == "mechanism":
            try:
                _mech.from_text(str(v))
            except Exception as e:
                dropped.append("丢弃 mechanism：解析不通过（%s）"
                               % str(e).splitlines()[0][:60])
                continue
        if k == "template_key" and str(v) not in LIBRARY:
            dropped.append("丢弃 template_key=%r：模板库里没有这个" % v)
            continue
        kept[k] = v
    return kept, dropped


def understand_llm(text: str) -> Request:
    """先走规则翻译；装了后端就再让它补充**结构与问题**（数值一律丢弃）。"""
    req = understand(text)
    if _BACKEND is None:
        req.warnings.append("LLM 后端未启用（默认关）——"
                            "以上全部来自规则翻译，零 LLM，无编数字风险。")
        return req
    try:
        raw = _BACKEND(text)
    except Exception as e:
        req.warnings.append("LLM 后端调用失败（%s），已退回规则翻译的结果。" % e)
        return req
    kept, dropped = sanitize_backend_output(raw if isinstance(raw, dict) else {})
    for d in dropped:
        req.warnings.append("后端产出被过滤：" + d)
    if not req.template_key and kept.get("template_key"):
        req.template_key = kept["template_key"]
        req.warnings.append("体系由 LLM 后端建议为 %r —— **请自行确认**。"
                            % req.template_key)
        req.questions = _questions(req)
    for q in kept.get("questions", []) or []:
        if str(q) not in req.questions:
            req.questions.append("（LLM 补充）" + str(q))
    return req
