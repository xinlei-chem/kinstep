# -*- coding: utf-8 -*-
# KinSTEP - a modality-agnostic microkinetic solver.
# Copyright (C) 2026 The KinSTEP authors.
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.  It is distributed WITHOUT ANY WARRANTY; see the GNU General Public
# License (LICENSE) for details.
r"""[English] Naming conventions and collision checking, the first gate: a species may not be
named like a physical quantity the solver uses.

(The commentary below is in Chinese: it records why each choice was made
and which measurement forced it.)

命名约定与**撞名检查** —— 本包的第一道安全闸。

为什么单独做一个模块
--------------------
催化/电化学这两套符号体系天生互相踩：

    k  : 速率常数 / 玻尔兹曼常数 / 循环下标
    R  : 气体常数 / 反应速率 / reverse（逆向）
    F  : 法拉第常数 / 电场 / 力
    E  : 能量 / **电极电势** / 电场          ← 最危险的一个
    Z  : 近邻配位数（Koper）/ **可逆度 z**（CATKINAS）/ 配分函数
    Q  : **覆盖度**（CATKINAS 用 Q_）/ 电荷 / 热
    n  : 电子数 / 下标 / 反应级数
    j  : **电流密度** / 循环下标
    α β: 传递系数 / BEP 斜率

实证：Ager 的 notebook 里

    def k2(E, p): return p.k2zero*np.exp((α_1*q*E)/(kB*T))
                                          ^^^ 本该是 α_2

因为算例中 α_1 == α_2 == 0.5，**错了也看不出来**；谁改了 α_2 就会被静默忽略。

⇒ 本包的三条硬规矩
------------------
1. **不用单字母做有物理含义的变量名**（下标 i/j/k 除外，且只在局部循环里）。
2. **同一概念全包只用一个名字**，见下面 `GLOSSARY`。
3. **用户输入的物种名/位点名若撞上保留名，解析阶段直接报错**，不做"自动改名"。
   —— CATKINAS 的 INPUT 是 MATLAB 代码，物种若叫 `T` 会直接把温度覆盖掉，
     这种坑在本包里用 `check_user_names()` 堵住。
"""

#: 本包的统一词汇表：概念 → 唯一变量名（**改名必须改这里，不许就地另起**）
GLOSSARY = {
    # —— 物理常数（全大写，带后缀，不与任何物理量重名）——
    "玻尔兹曼常数": "K_B / K_B_EV",
    "元电荷": "Q_E",
    "法拉第常数": "F_FARADAY",
    "气体常数": "R_GAS",
    "阿伏伽德罗常数": "N_A",
    "普朗克常数": "H_PLANCK",

    # —— 状态变量 ——
    "温度(K)": "T_K",
    "覆盖度": "theta",            # 不用 Q（撞电荷/热）
    "空位覆盖度": "theta_vac",
    "气相分压(无量纲)": "p_gas",
    "液相浓度(无量纲)": "c_liq",

    # —— 动力学 ——
    "正向速率常数": "kf_step",     # 不用裸 k
    "逆向速率常数": "kr_step",
    "平衡常数": "Keq_step",
    "净速率": "rate_net",
    "正向速率": "rate_fwd",
    "逆向速率": "rate_rev",
    "可逆度 r_rev/r_fwd": "zrev",  # 不用 Z（撞配位数）

    # —— 能量（一律 eV）——
    "正向活化自由能": "Ga_fwd",
    "逆向活化自由能": "Ga_rev",
    "反应自由能": "dG_rxn",
    "表观活化能": "Ea_app",

    # —— 电化学（CATKINAS 完全没有的部分）——
    "电极电势(V)": "U_elec",       # **不用 E**（E 永远指能量）
    "参比标度": "U_ref",           # 'SHE' / 'RHE'
    "电子计量系数(产物为正)": "nu_elec",
    "电荷传递系数": "beta_ct",     # 不用裸 beta（撞 BEP 斜率）
    "电流密度": "j_curr",          # 不用裸 j（撞下标）
    "位点面密度": "gamma_site",
    "pH": "pH_val",
    "电场(V/Å)": "field_V_per_A",
    "Tafel 斜率(V/dec)": "tafel_slope",

    # —— 结构/统计 ——
    "近邻配位数": "n_neighbor",    # 不用 Z
    "位点编号": "site_id",
    "基元步编号": "step_id",

    # —— BEP / scaling ——
    "BEP 斜率": "alpha_bep",
    "BEP 截距": "beta_bep",
    "能量描述符": "desc_energy",   # 不用 E1/E2
}

#: 保留名：用户的物种名/位点名**不得**与之相同（大小写不敏感比较）
RESERVED_NAMES = {
    # 常数
    "k_b", "k_b_ev", "q_e", "f_faraday", "r_gas", "n_a", "h_planck", "ev",
    # 状态与动力学
    "t_k", "theta", "theta_vac", "p_gas", "c_liq",
    "kf_step", "kr_step", "keq_step", "rate_net", "rate_fwd", "rate_rev", "zrev",
    "ga_fwd", "ga_rev", "dg_rxn", "ea_app",
    # 电化学
    "u_elec", "u_ref", "nu_elec", "beta_ct", "j_curr", "gamma_site",
    "ph_val", "field_v_per_a", "tafel_slope",
    # 结构
    "n_neighbor", "site_id", "step_id", "alpha_bep", "beta_bep", "desc_energy",
    # ★ 2026-10-22 修正：**不能一刀切禁所有单字母** ——
    #   `H*`（吸附氢）、`N*`（吸附氮）、`O*`、`C*` 是催化里最常见的物种名，
    #   初版把它们全拦了，合成氨机理一行都解析不了（测试当场抓到）。
    #
    #   根本差别：在 CATKINAS 里物种名**就是 MATLAB 变量**，`T` 会把温度覆盖掉；
    #   在 KinSTEP 里物种名**永远只是字典的键，不会被当代码求值** ⇒ 风险小得多。
    #
    #   新规则：**单字母只有在它是合法元素符号时才放行**。
    #   下面这些单字母不是元素符号，却是物理量的通用记号 ⇒ 仍然禁止。
    "t", "r", "e", "z", "q", "j", "g", "a", "d", "l", "m", "x",
    # python 关键字里最容易被当成分子名的
    "lambda", "none", "true", "false", "in", "is", "as", "or", "and", "not",
}

#: 单字母的合法元素符号 —— 这些可以当物种名
SINGLE_LETTER_ELEMENTS = {"h", "b", "c", "n", "o", "f", "p", "s", "k", "v", "y", "i", "w", "u"}


def is_element_symbol(name) -> bool:
    """是不是单字母元素符号（大小写不敏感）。"""
    return str(name).strip().lower() in SINGLE_LETTER_ELEMENTS


# ══════════════════════ 显式物种前缀（2026-10-22 老板提议）══════════════════════
#
# 老板的建议：给元素/化学式加一个首字母前缀（如 `E_C` 表示碳、`E_CO` 表示一氧化碳），
# 从而彻底避免物种名与物理量撞车。
#
# ★ 评估结论：**作为"可选逃生门"非常有用，作为"强制写法"不合适。** 三条理由：
#
#   ① 本包的撞名风险本来就比 CATKINAS 小一个量级 —— 在 CATKINAS 里物种名**就是
#      MATLAB 变量**（叫 `T` 会把温度覆盖掉）；在 KinSTEP 里物种名**永远只是字典的键**，
#      不会被当代码求值。所以"强制加前缀"付的代价换不回等量的安全。
#   ② 本包的首要设计目标是"照着纸上写的反应式写"（老板原话：别再像 CATKINAS 那样
#      填复杂输入脚本）。强制 `E_CO*1 + E_O*2 -> E_CO2 + *1 + *2` 既难写也难读，
#      和这个目标直接冲突。
#   ③ ⚠ **`E_` 这个前缀本身有个真问题**：催化文献里 `E_CO` 的标准含义是
#      "CO 的吸附能"（本包 naming.py 开头就把 `E` 列为"最危险的一个"）。
#      所以这里**同时接受 `Sp_`**（species）作为无歧义的首选写法，`E_` 作为别名保留。
#
# ★ 但这个建议解决了一个**我原本留着的真洞**：`_parse_formula` 要求名字能分解成
#   元素符号，所以 `Cat`（催化剂）、`Dye`（染料）、`Nafion`、`Enz` 这类名字**根本解析不了**；
#   我在报错信息里承诺了一个 `species_formula` 参数，而它**当时并不存在**。
#   现在前缀 + `species_formula` 一起把这个洞补上了。
#
#: 可选的"这是一个化学物种"前缀。首选 `Sp_`；`E_` 为别名（见上面第③条）
SPECIES_PREFIXES = ("Sp_", "sp_", "E_")


def strip_species_prefix(name):
    """去掉显式物种前缀。返回 `(裸名, 是否带过前缀)`。

    带过前缀的名字有两项特权：
      · **豁免保留名检查** —— `Sp_T` 可以正当地表示一个叫 T 的物种；
      · **允许非元素式名字** —— `Sp_Cat`，组成另给（`species_formula`）。

    ★ 裸名就是规范名：`E_CO` 与 `CO` 是**同一个物种**，可以在一份机理里混写
      而不会变成两个物种（这一点有回归测试）。
    """
    s = str(name).strip()
    for p in SPECIES_PREFIXES:
        if s.startswith(p) and len(s) > len(p):
            return s[len(p):], True
    return s, False


#: 电子在基元反应里的写法（只认这几种，避免和物种名混淆）
ELECTRON_TOKENS = {"e-", "e−", "e^-", "electron"}

#: 质子/氢氧根的标准写法（便于自动套用 pH）
PROTON_TOKENS = {"H+", "H^+"}
HYDROXIDE_TOKENS = {"OH-", "OH^-", "OH−"}


class NameCollisionError(ValueError):
    """用户给的名字撞上了保留名。**不自动改名，直接报错。**"""


def check_user_names(names, kind="物种"):
    """检查用户输入的物种名/位点名是否撞保留名；撞了就抛异常。

    参数
    ----
    names : 可迭代的字符串
    kind  : 用于报错信息的中文名（"物种" / "位点"）

    为什么不自动改名
    ----------------
    自动改名会让输入文件与输出结果对不上号，而且用户不会注意到。
    **报错逼用户改**，是这里唯一安全的做法。
    """
    bad = []
    seen = {}
    for nm in names:
        key = str(nm).strip().lower()
        # ★ 单字母元素符号（H/N/O/C/S/…）一律放行 —— 它们是正当的化学物种名
        if len(key) == 1 and is_element_symbol(key):
            seen[key] = nm
            continue
        if key in RESERVED_NAMES:
            bad.append((nm, "撞保留名"))
        if key in seen and seen[key] != nm:
            bad.append((nm, "与 '%s' 仅大小写不同，易混" % seen[key]))
        seen[key] = nm
    if bad:
        msg = ["%s命名冲突（%d 处）：" % (kind, len(bad))]
        for nm, why in bad:
            msg.append("  '%s' —— %s" % (nm, why))
        msg.append("保留名清单见 kinstep/naming.py 的 RESERVED_NAMES；请改用更明确的名字，"
                   "例如用 'Tcat' 代替 'T'、'Hplus' 代替 'H'。")
        raise NameCollisionError("\n".join(msg))
    return True


def audit_module_namespace(mod, allow=()):
    """审计一个模块的顶层命名：挑出"有物理含义却只有 1–2 个字符"的名字。

    这是给**开发者自己**用的回归测试 —— 防止以后手一滑又写出 `E = ...`。
    """
    offenders = []
    for nm in dir(mod):
        if nm.startswith("_") or nm in allow:
            continue
        obj = getattr(mod, nm)
        if callable(obj) or isinstance(obj, type):
            continue
        if len(nm) <= 2 and nm not in ("CM", "MA", "MV", "EV"):
            offenders.append(nm)
    return offenders
