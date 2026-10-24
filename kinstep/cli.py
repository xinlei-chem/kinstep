# -*- coding: utf-8 -*-
# KinSTEP - a modality-agnostic microkinetic solver.
# Copyright (C) 2026 The KinSTEP authors.
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.  It is distributed WITHOUT ANY WARRANTY; see the GNU General Public
# License (LICENSE) for details.
r"""[English] Command-line entry point: templates, new, check, explain, ask, run, drivers.
Exit codes: 0 success, 2 malformed input, 3 provenance refusal, 4 driver
misconfiguration.

(The commentary below is in Chinese: it records why each choice was made
and which measurement forced it.)

命令行入口（T8.5–8.7）。

本包的首要承诺是"**别再像 CATKINAS 那样填复杂输入脚本**"。
光有 Python API 兑现不了这句话 —— 使用者里很多人不写 Python。
所以 CLI 的设计目标不是"包一层参数解析"，而是：

  kinstep check   <file>   把能查的都查了，**每条报错都给出怎么改**
  kinstep explain <file>   用人话说清楚"这份输入在算什么、哪些值是默认的"
  kinstep run     <file>   求解 → 报告 → 出图 → 导出可作图的原始数据
  kinstep new     [--scan] 生成一份**能直接跑**的模板（不是空壳）
  kinstep drivers          列出可用模态及其参数（和它们各自的"主人"）

★ `explain` 为什么单独做一个命令
--------------------------------
MKM 最危险的不是报错，是**不报错但结果不对**：某个默认值被静默采用、
某步其实处在准平衡而用户不知道、pH 被重复计了一次……
`explain` 把这些**在算之前**摊开说，是 `provenance.py` 那套来源标签的出口。

★ 界面语言（2026-10-23 加）
---------------------------
所有面向用户的文字都包在 `i18n.t("中文", "English")` 里。
**默认仍是中文**（既有回归逐字比对中文输出，不受影响）；
`KINSTEP_LANG=en` 切英文。理由见 `i18n.py` 的模块说明。
"""
from __future__ import annotations

import argparse
import os
import sys
from typing import Dict, List, Optional

import numpy as np

from .i18n import t

#（版本号在 __init__.py 里，CLI 不直接依赖）

OK, WARN, BAD = "✓", "⚠", "✗"


def _bar(title):
    return "\n" + title + "\n" + "─" * max(len(title), 30)


# ════════════════════════════ check ════════════════════════════
def cmd_check(args) -> int:
    from .inputs import load, InputError
    from .provenance import audit, report_defaults
    try:
        sp = load(args.file)
    except InputError as e:
        print("%s %s\n%s" % (BAD, t("输入有问题：", "invalid input:"), e))
        return 2
    except Exception as e:
        print("%s %s%s" % (BAD, t("解析失败：", "could not parse: "), e))
        return 2

    print(_bar(t("① 机理", "[1/5] mechanism")))
    print(sp.mech.summary())
    print("  %s %s" % (OK, t("四道闸全过（元素/电荷配平、位点守恒、撞保留名、箭头）",
                             "all four gates passed "
                             "(element/charge balance, site balance, "
                             "reserved-name clash, arrow syntax)")))

    print(_bar(t("② 工况", "[2/5] operating conditions")))
    print("  %s T_K = %g K" % (t("温度", "temperature"), sp.conditions["T_K"]))
    for nm, d in ((t("气相分压", "gas partial pressures"), sp.p_gas),
                  (t("液相浓度", "liquid concentrations"), sp.c_liq)):
        if d:
            print("  %s%s%s" % (nm, t("：", ": "),
                                "  ".join("%s=%g" % kv for kv in d.items())))
    # ★ 液相也要查（T14 抓到：HER 模板漏了 c_liq ⇒ H+ 活度为 0 ⇒ 速率恒为 0，
    #   而 check 只查气相，一声不响地放过了）
    for kind, species, given in (
            (t("气相分压", "a gas partial pressure"), sp.mech.gas_species, sp.p_gas),
            (t("液相浓度", "a liquid concentration"), sp.mech.liq_species, sp.c_liq)):
        miss = [s for s in species
                if s not in given and sp.mech.species[s].formula not in given]
        if not miss:
            continue
        react = [s for s in miss
                 if any(s in st.reactants for st in sp.mech.steps)]
        print("  %s %s" % (WARN, t("这些物种没给%s，将按 0 处理：%s" % (kind, miss),
                                   "no %s given for these species; "
                                   "they will be taken as 0: %s" % (kind, miss))))
        if react:
            print("     ⛔ %s" % t("其中 %s 出现在**反应物侧** ⇒ 相应步的速率会恒为 0。" % react,
                                   "%s appear on the REACTANT side, so the rate of "
                                   "those steps will be identically zero." % react))
        else:
            print("     %s" % t("它们只在产物侧，按 0 处理通常没问题。",
                                "they only appear as products, so 0 is usually fine."))
    for nm, d in ((t("气相", "gas phase"), sp.p_gas),
                  (t("液相", "liquid phase"), sp.c_liq)):
        known = {s for s in sp.mech.species_order} | {sp.mech.species[s].formula
                                                      for s in sp.mech.species_order}
        extra = [k for k in d if k not in known]
        if extra:
            print("  %s %s" % (WARN, t("%s里有机理中不存在的物种：%s —— 多半是拼写错了"
                                       % (nm, extra),
                                       "species in the %s that do not occur in the "
                                       "mechanism: %s - most likely a typo"
                                       % (nm, extra))))

    print(_bar(t("③ 能量与来源", "[3/5] energetics and provenance")))
    if sp.energetics is not None:
        try:
            defaulted = sp.energetics.check(strict_energy=True)
            print("  %s %s" % (OK, t("来源闸通过（没有未核实的能量数据）",
                                     "provenance gate passed "
                                     "(no unverified energy data)")))
            print("  " + report_defaults(defaulted).replace("\n", "\n  "))
        except Exception as e:
            print("%s %s\n%s" % (BAD, t("来源闸拦住了：",
                                        "provenance gate rejected the input:"), e))
            return 3
    else:
        print("  %s" % t("走 scaling + BEP：", "using scaling + BEP:"))
        print("  " + sp.scaling.describe().replace("\n", "\n  "))

    print(_bar(t("④ driver（模态）", "[4/5] driver (modality)")))
    if not sp.drivers:
        print("  %s" % t("没启用 driver ⇒ 纯热催化",
                         "no driver enabled => purely thermal catalysis"))
    else:
        from .driver import make_driver, DriverStack
        try:
            st = DriverStack([make_driver(n, **c) for n, c in sp.drivers])
            st.validate(sp.mech)
            print("  %s %s" % (OK, st.describe()))
        except Exception as e:
            print("%s %s%s" % (BAD, t("driver 配置有问题：",
                                      "driver configuration is invalid: "), e))
            return 4

    print(_bar(t("⑤ 要做什么", "[5/5] what will be computed")))
    if sp.is_scan:
        tot = int(np.prod([len(a) for a in sp.scan_axes]))
        print("  %s" % t("%d 维描述符扫描，共 %d 点：" % (len(sp.scan_axes), tot),
                         "%dD descriptor scan, %d points in total:"
                         % (len(sp.scan_axes), tot)))
        for nm, a in zip(sp.scaling.desc_names, sp.scan_axes):
            print("     %s: %.3f … %.3f（%d %s）"
                  % (nm, a[0], a[-1], len(a), t("点", "points")))
    else:
        print("  %s" % t("单点稳态求解", "single-point steady-state solve"))
    print("\n%s %s" % (OK, t("这份输入可以跑。", "this input is ready to run.")))
    return 0


# ════════════════════════════ explain ════════════════════════════
def cmd_explain(args) -> int:
    from .inputs import load
    from .provenance import report_defaults
    sp = load(args.file)
    m = sp.mech
    print(_bar(t("这份输入在算什么", "what this input computes")))
    print("  %s" % t("体系：%d 个基元步，%d 个物种，%d 种位点"
                     % (m.n_step, len(m.species_order), len(m.sites)),
                     "system: %d elementary steps, %d species, %d site type(s)"
                     % (m.n_step, len(m.species_order), len(m.sites))))
    print("  %s" % t("表面物种 %s；空位 %s"
                     % (m.surface_species, [m.vacancy_of(s) for s in m.sites]),
                     "surface species %s; vacancies %s"
                     % (m.surface_species, [m.vacancy_of(s) for s in m.sites])))
    ec = [i + 1 for i, s in enumerate(m.steps) if s.params.get("nu_elec")]
    if ec:
        print("  %s" % t("电化学步：%s（nu_elec 由方程自动数出，没让你填）" % ec,
                         "electrochemical steps: %s (nu_elec is counted from the "
                         "equation, you never type it)" % ec))
    irr = [i + 1 for i, s in enumerate(m.steps) if not s.reversible]
    if irr:
        print("  %s" % t("不可逆步（用 `->` 写的，逆向常数置 0）：%s" % irr,
                         "irreversible steps (written with `->`, reverse constant "
                         "set to 0): %s" % irr))
    goal_zh = ("扫描 %s，画火山图" % sp.scaling.desc_names) if sp.is_scan \
        else "求稳态覆盖度与各步速率"
    goal_en = ("scan %s and draw the volcano plot" % sp.scaling.desc_names) if sp.is_scan \
        else "solve for steady-state coverages and per-step rates"
    print("  %s%s" % (t("目标：", "goal: "), t(goal_zh, goal_en)))

    print(_bar(t("哪些值是你给的，哪些是默认的",
                 "which values you supplied, and which are defaults")))
    if sp.energetics is not None:
        d = sp.energetics.check(strict_energy=False)
        print("  " + report_defaults(d).replace("\n", "\n  "))
    else:
        print("  %s" % t("能量由 scaling + BEP 生成，系数全部来自输入文件（无默认值）",
                         "energies come from scaling + BEP; every coefficient is "
                         "taken from the input file (no defaults)"))
    off = sp.raw.get("thermo", {})
    if off:
        print("  %s" % t("热化学常数偏移由输入文件给定：",
                         "thermochemical offsets given in the input file:"))
        for k, v in off.items():
            print("     %s = %s" % (k, v))
    else:
        print("  %s %s" % (WARN, t("没给热化学偏移 ⇒ 按 0 处理。",
                                   "no thermochemical offsets given => taken as 0.")))
        print("     %s" % t("若你的 Ga/ΔG 是**纯电子能**（没含气相熵、零点能），结果会系统性偏。",
                            "if your Ga/dG are bare electronic energies (no gas-phase "
                            "entropy, no zero-point energy), the result will be "
                            "systematically biased."))

    print(_bar(t("算之前要知道的三件事", "three things to know before you run")))
    print("  1. %s" % t("速率跨几十个数量级是常态 —— 判断好坏看 log10(速率)，不看绝对值。",
                        "rates spanning tens of orders of magnitude are normal - "
                        "judge by log10(rate), not by the absolute value."))
    print("  2. %s" % t("准平衡步（1−z → 0）的净速率在双精度下可能算不准；",
                        "net rates of quasi-equilibrated steps (1-z -> 0) can lose "
                        "all significant digits in double precision;"))
    print("     %s" % t("solver.precision = \"auto\" 会**自动升到高精度**并在报告里标出来。",
                        "solver.precision = \"auto\" escalates to extended precision "
                        "automatically and says so in the report."))
    print("  3. %s" % t("结果里的 `ΣX_B = 1`（Campbell 求和规则）是免费自检 —— 它不对，别信结论。",
                        "the `sum X_B = 1` line (Campbell's sum rule) is a free "
                        "self-check - if it fails, do not trust the conclusion."))
    ph = any(n in ("H+", "OH-") for n in m.species_order)
    ur = str(sp.conditions.get("U_ref", "")).upper()
    if ph and ur == "RHE":
        print("\n%s %s" % (WARN, t("**pH 可能被重复计算**：机理里显式写了 H+，而电位又用 RHE 标度。",
                                   "pH may be counted TWICE: H+ appears explicitly in "
                                   "the mechanism while the potential is on the RHE "
                                   "scale.")))
        print("   %s" % t("二者只能选一个：写 H+ 就用 SHE；用 RHE 就别在机理里写 H+。",
                          "pick one: with explicit H+ use SHE; with RHE do not write "
                          "H+ in the mechanism."))
    return 0


# ════════════════════════════ run ════════════════════════════
def cmd_run(args) -> int:
    from .inputs import load
    from .solve import steady_state
    from .driver import make_driver, DriverStack
    sp = load(args.file)
    drv = DriverStack([make_driver(n, **c) for n, c in sp.drivers]) if sp.drivers else None
    state = dict(sp.conditions)
    prec = sp.solver.get("precision", args.precision)
    dps = int(sp.solver.get("dps", args.dps))
    outdir = args.out or (os.path.splitext(os.path.basename(args.file))[0] + "_out")
    os.makedirs(outdir, exist_ok=True)

    if sp.is_scan:
        from .scan import scan
        from . import viz
        res = scan(sp.mech, sp.scaling, state, sp.scan_axes, driver=drv,
                   p_gas=sp.p_gas, c_liq=sp.c_liq,
                   A_pre=sp.solver.get("A_pre"), precision=prec, dps=dps)
        print(_bar(t("扫描结果", "scan results")))
        print(res.summary())
        dat = os.path.join(outdir, "scan.dat")
        res.save_dat(dat)
        print("\n  %s%s" % (t("原始数据 → ", "raw data -> "), dat))
        if viz.available() and not args.no_plot:
            if len(res.shape) == 1:
                viz.volcano_1d(res, path=os.path.join(outdir, "volcano.png"))
                viz.coverage_1d(res, path=os.path.join(outdir, "coverage.png"))
                print("  %s%s/volcano.png, coverage.png"
                      % (t("图 → ", "figures -> "), outdir))
            else:
                viz.volcano_2d(res, path=os.path.join(outdir, "volcano_2d.png"))
                print("  %s%s/volcano_2d.png" % (t("图 → ", "figure -> "), outdir))
        elif not viz.available():
            print("  %s" % t("（没装 matplotlib，只出数据不出图：pip install matplotlib）",
                             "(matplotlib not installed: data only, no figures - "
                             "pip install matplotlib)"))
        if not res.converged.all():
            print("\n%s %s" % (WARN, t("有 %d 个点没收敛 —— 结论要小心"
                                       % int((~res.converged).sum()),
                                       "%d point(s) did not converge - be careful "
                                       "with the conclusion"
                                       % int((~res.converged).sum()))))
            return 1
        return 0

    r = steady_state(sp.mech, sp.energetics, state=state, driver=drv,
                     p_gas=sp.p_gas, c_liq=sp.c_liq,
                     tol=float(sp.solver.get("tol", 1e-12)),
                     precision=prec, dps=dps)
    print(_bar(t("稳态结果", "steady-state result")))
    print(r.report())
    if r.obs:
        print("  %s" % t("模态可观测量：", "modality observables: ")
              + "  ".join("%s=%.6g" % kv for kv in r.obs.items()))
    ok, rel = r.flux_check(sp.mech)
    print("  %s %s%.2e" % (OK if ok else WARN,
                           t("通量守恒自检：最大相对不自洽 ",
                             "flux-conservation self-check: max relative residual "),
                           rel))
    for w in r.precision_warnings():
        print("  %s %s" % (WARN, w))
    if args.drc:
        from .analyze import drc, calibrate_delta
        print(_bar(t("速率控制度（DRC）", "degree of rate control (DRC)")))
        best, _ = calibrate_delta(sp.mech, sp.energetics, state, drv, sp.p_gas, sp.c_liq)
        d = drc(sp.mech, sp.energetics, state, drv, sp.p_gas, sp.c_liq,
                delta=best, richardson=True, precision=prec, dps=dps)
        print("  %s" % t("（差分步长 δ=%.1e 由 Campbell 求和规则自动选出）" % best,
                         "(finite-difference step delta=%.1e chosen automatically "
                         "by Campbell's sum rule)" % best))
        print(d.report(target=sp.mech.n_step - 1))
        s = d.sum_B[sp.mech.n_step - 1]
        print("  %s %s%.12f%s%s%.1e%s"
              % (OK if abs(s - 1) < 1e-3 else BAD,
                 t("求和规则 ΣX_B = ", "sum rule  sum X_B = "), s,
                 t("（", " ("), t("偏差 ", "deviation "), abs(s - 1),
                 t("）", ")")))
    with open(os.path.join(outdir, "result.txt"), "w", encoding="utf-8") as f:
        f.write(r.report() + "\n")
    print("\n  %s%s/result.txt" % (t("文字报告 → ", "text report -> "), outdir))
    return 0 if r.converged else 1


# ════════════════════════════ ask ════════════════════════════
def cmd_ask(args) -> int:
    """用一句话描述需求 → 结构 + 问题（**零 LLM**，默认不编任何数字）。"""
    from .nl import understand, understand_llm, backend_enabled
    text = " ".join(args.text)
    req = (understand_llm if args.llm else understand)(text)
    print(_bar(t("理解结果", "what was understood")))
    print(req.report())
    if args.llm and not backend_enabled():
        print("\n  %s" % t("（--llm 给了，但没装后端 ⇒ 仍是纯规则翻译。"
                           "装法见 kinstep.nl.set_backend，默认关是有意的）",
                           "(--llm was given but no backend is installed, so this is "
                           "still the pure rule-based translation. See "
                           "kinstep.nl.set_backend; off by default on purpose.)"))
    if req.template_key and args.write:
        from .templates import get
        from .inputs import SUFFIX
        path = args.write if args.write != "-" else (req.template_key + SUFFIX)
        tpl = get(req.template_key)
        # ★ 把**用户自己说出来的**工况覆盖进去（标 user_given）；能量一律不写。
        #   用 to_toml(overrides=…) 而不是在外面插行 —— 插行会造成重复键，
        #   TOML 直接报 "Cannot overwrite a value"（T15 当场踩到的真 bug）。
        got = {k: (v if isinstance(v, str) else v.value)
               for k, v in req.conditions().items()}
        open(path, "w", encoding="utf-8").write(tpl.to_toml(overrides=got))
        print("\n%s %s" % (OK, t("已写出 %s（机理 + 你说过的工况；**能量槽位仍留空**）" % path,
                                 "wrote %s (mechanism + the conditions you stated; "
                                 "the energy slots are still empty)" % path)))
        print("   %s kinstep check %s" % (t("下一步：", "next:"), path))
    return 0


def cmd_templates(args) -> int:
    from .templates import summary, get, list_templates
    if args.name:
        try:
            print(get(args.name).describe())
        except KeyError as e:
            print("%s %s" % (BAD, e)); return 2
        return 0
    print(_bar(t("离线模板库", "offline template library")))
    print(summary())
    print("\n★ %s" % t("红线：模板**只给机理与问题，不给能量数值**。",
                       "hard rule: a template gives you the mechanism and the "
                       "questions, never energy values."))
    print("  %s" % t("机理是公开知识（写错肉眼能看出来）；能量是你体系的私有数据",
                     "a mechanism is public knowledge (a mistake is visible by eye); "
                     "energies are private data for your own system"))
    print("  %s" % t("（填错**只会得到一条看着很合理的火山曲线**）。",
                     "(a wrong energy only gives you a perfectly plausible-looking "
                     "volcano curve)."))
    return 0


def cmd_new(args) -> int:
    from .inputs import template, SUFFIX
    if args.template:
        from .templates import get
        try:
            tpl = get(args.template)
        except KeyError as e:
            print("%s %s" % (BAD, e)); return 2
        path = args.file or (args.template + SUFFIX)
        if os.path.exists(path) and not args.force:
            print("%s %s" % (BAD, t("%s 已存在。加 --force 覆盖。" % path,
                                    "%s already exists. Use --force to overwrite."
                                    % path)))
            return 2
        with open(path, "w", encoding="utf-8") as f:
            f.write(tpl.to_toml())
        m = tpl.parse()
        title = tpl.title_localized()
        print("%s %s" % (OK, t("已生成 %s —— %s" % (path, title),
                               "wrote %s - %s" % (path, title))))
        # ★ 槽位数是**每步两个**（Ga_fwd + dG_rxn，`to_toml` 就是这么写的）。
        #   原来这里填的是 m.n_step，于是 HER 上 `new` 说「3 个槽位」、
        #   `check` 说「缺 6 项」—— 同一份输入上两条命令自己打自己。
        n_slot = 2 * m.n_step
        print("   %s" % t("机理 %d 步已填好；**%d 个能量槽位要你填**（每个都标了要算什么）"
                          % (m.n_step, n_slot),
                          "all %d mechanism steps are filled in; %d energy slots are "
                          "left for you (each one says what to compute)"
                          % (m.n_step, n_slot)))
        print("   %s%s" % (t("出处：", "source: "), tpl.reference_localized()))
        for n in tpl.notes_localized():
            print("   · %s" % n)
        print("\n   %s kinstep check %s   # %s"
              % (t("下一步：", "next:"), path,
                 t("它会逐条告诉你还缺什么",
                   "it tells you, item by item, what is still missing")))
        return 0
    kind = "scan" if args.scan else "single"
    path = args.file or ("ammonia" + SUFFIX)
    if os.path.exists(path) and not args.force:
        print("%s %s" % (BAD, t("%s 已存在。加 --force 覆盖，或换个文件名。" % path,
                                "%s already exists. Use --force, or pick another "
                                "file name." % path)))
        return 2
    with open(path, "w", encoding="utf-8") as f:
        f.write(template(kind))
    print("%s %s" % (OK, t("已生成 %s（%s模板，**可以直接跑**）"
                           % (path, "扫描" if args.scan else "单点"),
                           "wrote %s (%s template, ready to run as-is)"
                           % (path, "scan" if args.scan else "single-point"))))
    print("  %s" % t("下一步：", "next:"))
    print("    kinstep check   %s     # %s" % (path, t("先查一遍", "check it first")))
    print("    kinstep explain %s     # %s" % (path, t("看它要算什么",
                                                       "see what it will compute")))
    print("    kinstep run     %s     # %s" % (path, t("开算", "run it")))
    return 0


# ════════════════════════════ drivers ════════════════════════════
def cmd_drivers(args) -> int:
    from . import drivers as _d                      # 触发注册
    from .driver import available_drivers, make_driver, param_owner_table
    print(_bar(t("已注册的模态（driver）", "registered modalities (drivers)")))
    for n in available_drivers():
        try:
            d = make_driver(n)
        except Exception:
            from .driver import _REGISTRY
            d = _REGISTRY[n]
        cap = d.capabilities() if hasattr(d, "capabilities") else {}
        print("  %-9s %s%-44s %s%s"
              % (n, t("工况参数 ", "state params "),
                 str(list(getattr(d, "state_params", ()))),
                 t("逐步参数 ", "per-step params "),
                 str(list(getattr(d, "step_params", ())))))
        flags = [k for k in ("needs_theta", "is_additive", "splits_state")
                 if cap.get(k) is not None]
        if flags:
            print("            %s" % t("能力：", "capabilities: ") + "  ".join(
                "%s=%s" % (k, cap.get(k)) for k in flags)
                + ((t("  产出：", "  provides: ") + str(list(cap.get("provides", ()))))
                   if cap.get("provides") else ""))
    print(_bar(t("参数名归属（查“这个参数是谁要的”）",
                 "parameter ownership (which driver asks for what)")))
    for p, owner in sorted(param_owner_table().items()):
        print("  %-16s ← %s" % (p, owner))
    print("\n%s" % t("在输入文件里这样启用：", "enable it in the input file like this:"))
    print("  [drivers.electro]")
    print("  U_elec = 0.3")
    print("  U_ref  = \"SHE\"")
    print("  gamma_site = 2.19e-5")
    return 0


# ════════════════════════════ 入口 ════════════════════════════
def build_parser():
    p = argparse.ArgumentParser(
        prog="kinstep",
        description=t("KinSTEP —— 微动力学：热 / 电 / 光 / 超声催化"
                      "（Sono-Thermal-Electro-Photo）",
                      "KinSTEP - microkinetic modelling for thermal / electro / "
                      "photo / sono catalysis (Sono-Thermal-Electro-Photo)"),
        epilog=t("先 `kinstep new` 生成模板，再 check / explain / run。",
                 "start with `kinstep new`, then check / explain / run."))
    sub = p.add_subparsers(dest="cmd", required=True)

    c = sub.add_parser("check", help=t("检查输入文件：机理、工况、能量来源、driver",
                                       "check an input file: mechanism, conditions, "
                                       "energy provenance, driver"))
    c.add_argument("file")
    c.set_defaults(func=cmd_check)

    e = sub.add_parser("explain", help=t("用人话说清楚这份输入在算什么、哪些值是默认的",
                                         "explain in plain words what the input "
                                         "computes and which values are defaults"))
    e.add_argument("file")
    e.set_defaults(func=cmd_explain)

    r = sub.add_parser("run", help=t("求解并出图、导出原始数据",
                                     "solve, plot, and export the raw data"))
    r.add_argument("file")
    r.add_argument("--out", help=t("输出目录（默认 <文件名>_out）",
                                   "output directory (default <file>_out)"))
    r.add_argument("--drc", action="store_true",
                   help=t("同时算速率控制度", "also compute the degree of rate control"))
    r.add_argument("--precision", default="auto", choices=("double", "auto", "high"))
    r.add_argument("--dps", type=int, default=50,
                   help=t("高精度位数（默认 50）",
                          "digits for extended precision (default 50)"))
    r.add_argument("--no-plot", action="store_true",
                   help=t("只出数据不出图", "export data only, no figures"))
    r.set_defaults(func=cmd_run)

    n = sub.add_parser("new", help=t("生成一份可直接跑的输入模板",
                                     "generate a ready-to-run input template"))
    n.add_argument("file", nargs="?")
    n.add_argument("--scan", action="store_true",
                   help=t("生成描述符扫描的模板",
                          "generate a descriptor-scan template"))
    n.add_argument("--template", help=t("从离线模板库生成（见 kinstep templates）",
                                        "generate from the offline library "
                                        "(see kinstep templates)"))
    n.add_argument("--force", action="store_true",
                   help=t("覆盖已存在的文件", "overwrite an existing file"))
    n.set_defaults(func=cmd_new)

    d = sub.add_parser("drivers", help=t("列出可用模态与它们的参数",
                                         "list the available modalities and their "
                                         "parameters"))
    d.set_defaults(func=cmd_drivers)

    tp = sub.add_parser("templates", help=t("列出离线模板库（常见体系的机理骨架）",
                                            "list the offline template library "
                                            "(mechanism skeletons for common systems)"))
    tp.add_argument("name", nargs="?", help=t("只看某一个模板的详情",
                                              "show the details of one template only"))
    tp.set_defaults(func=cmd_templates)

    a = sub.add_parser("ask", help=t("用一句话描述需求，得到结构 + 还缺什么（零 LLM）",
                                     "describe what you want in one sentence; get the "
                                     "structure plus what is still missing (no LLM)"))
    a.add_argument("text", nargs="+",
                   help=t("例：析氢反应 pH=0 电位 -0.1 V vs SHE",
                          "e.g. hydrogen evolution at pH 0, -0.1 V vs SHE"))
    a.add_argument("--write", nargs="?", const="-",
                   help=t("把结果写成 .kinstep.toml（能量槽位仍留空）",
                          "write the result to .kinstep.toml (energy slots stay empty)"))
    a.add_argument("--llm", action="store_true",
                   help=t("允许 LLM 后端补充结构（默认关；它的产出仍只能是结构与问题）",
                          "allow an LLM backend to fill in structure (off by default; "
                          "its output is still restricted to structure and questions)"))
    a.set_defaults(func=cmd_ask)
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return args.func(args)
    except KeyboardInterrupt:
        print("\n%s" % t("已中断", "interrupted"))
        return 130
    except Exception as e:
        print("%s %s: %s" % (BAD, type(e).__name__, e))
        if os.environ.get("KINSTEP_TRACEBACK"):
            raise
        print("  %s" % t("（设 KINSTEP_TRACEBACK=1 可看完整调用栈）",
                         "(set KINSTEP_TRACEBACK=1 for the full traceback)"))
        return 1


if __name__ == "__main__":
    sys.exit(main())
