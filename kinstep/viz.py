# -*- coding: utf-8 -*-
# KinSTEP - a modality-agnostic microkinetic solver.
# Copyright (C) 2026 The KinSTEP authors.
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.  It is distributed WITHOUT ANY WARRANTY; see the GNU General Public
# License (LICENSE) for details.
r"""[English] Plotting (optional, needs matplotlib): volcano curves and maps, coverages,
degrees of rate control.

(The commentary below is in Chinese: it records why each choice was made
and which measurement forced it.)

作图（T9）：1D 火山曲线、2D 火山图、覆盖度图、DRC 图。

三条规矩
--------
⓪ **图注默认英文。** 本包的图是要进论文的，而期刊图注一律英文；
   中文图注还要依赖系统装了中日韩字体，没装就渲染成一排豆腐块（首版实测就是这样，
   曲线完全正确却看不懂图例）。需要中文时调 `use_cjk_font()` 一次即可（见该函数）。
① **配色用 viridis 一类的 perceptually-uniform colormap，不用 jet。**
   jet 的亮度不单调，会在数据里**造出不存在的条纹和边界**，
   而且对色觉障碍读者不可读 —— 这是作图错误，不是审美偏好。
② **图之外一定给数据。** 每个作图函数都提示对应的 `ScanResult.save_dat()`；
   图是给人看的，数据是给后续分析和重画用的。

matplotlib 是**可选依赖**（核心只要 numpy+scipy）。没装时给出清楚的提示。
"""
from __future__ import annotations

from typing import Dict, List, Optional, Sequence

import numpy as np

_PLT = None


def available() -> bool:
    """matplotlib 在不在。**不在也能跑**，只是不出图（数据照常导出）。"""
    global _PLT
    if _PLT is None:
        try:
            import matplotlib
            matplotlib.use("Agg")          # 无显示环境也能出图
            import matplotlib.pyplot as plt
            _PLT = plt
        except ImportError:
            _PLT = False
    return _PLT is not False


def _plt():
    if not available():
        raise RuntimeError("作图需要 matplotlib：pip install matplotlib\n"
                           "（不想装也行：用 ScanResult.save_dat() 导出列数据，"
                           "拿 Origin/gnuplot 画）")
    return _PLT


def use_cjk_font(name=None, verbose=True):
    """切到中文图注。挑一个系统里装了的中日韩字体设进 matplotlib。

    为什么**不默认开**：① 论文图注要英文；② 字体是系统依赖，换台机器就可能没有，
    那时中文会静默渲染成豆腐块 —— **静默的坏输出比报错更糟**。
    所以默认英文（到哪都对），要中文时显式调用，并且**查不到字体就明确报错**。
    """
    plt = _plt()
    import matplotlib.font_manager as fm
    have = {f.name for f in fm.fontManager.ttflist}
    prefer = [name] if name else ["Noto Sans CJK SC", "Source Han Sans SC",
                                  "WenQuanYi Zen Hei", "Microsoft YaHei", "SimHei",
                                  "Droid Sans Fallback", "AR PL UMing CN"]
    for f in prefer:
        if f and f in have:
            plt.rcParams["font.sans-serif"] = [f] + list(
                plt.rcParams.get("font.sans-serif", []))
            plt.rcParams["axes.unicode_minus"] = False      # 负号也要能显示
            if verbose:
                print("图注字体已切到：%s" % f)
            return f
    raise RuntimeError(
        "系统里找不到中日韩字体（试过 %s）。\n"
        "  装一个： sudo apt install fonts-noto-cjk\n"
        "  或者就用默认的英文图注（到哪台机器都对）。" % prefer)


def _style(ax):
    ax.tick_params(direction="in", top=True, right=True, which="both")
    for sp in ax.spines.values():
        sp.set_linewidth(0.8)


def volcano_1d(res, step=-1, path=None, ylabel=r"$\log_{10}$ (rate / s$^{-1}$)",
               mark_peak=True, show_precision=True, figsize=(5.2, 3.8), dpi=200):
    """1D 火山曲线：log₁₀(速率) vs 描述符。

    `show_precision=True` 时把**自动升到高精度的点**标出来 ——
    让"哪一段是 double 算不出来的"在图上可见，而不是埋在日志里。
    """
    plt = _plt()
    if len(res.shape) != 1:
        raise ValueError("volcano_1d 只画 1D 扫描；2D 用 volcano_2d")
    x = res.desc_axes[0]
    r = np.abs(res.rate_of(step))
    ok = res.converged & np.isfinite(r) & (r > 0)
    with np.errstate(divide="ignore"):
        y = np.where(ok, np.log10(np.where(r > 0, r, 1.0)), np.nan)
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    ax.plot(x[ok], y[ok], "-", lw=1.6, color="#2c5f8a", zorder=3)
    if show_precision:
        hp = ok & (res.precision == "high")
        if hp.any():
            ax.plot(x[hp], y[hp], "o", ms=2.4, mfc="#d1495b", mec="none", zorder=4,
                    label="high-precision (%d pts)" % int(hp.sum()))
    if mark_peak and ok.any():
        i = int(np.nanargmax(np.where(ok, y, -np.inf)))
        ax.plot([x[i]], [y[i]], "*", ms=11, color="#edae49", mec="#8a6d1f",
                mew=0.6, zorder=5, label="peak: %s = %.3f" % (res.desc_names[0], x[i]))
    ax.set_xlabel("%s (eV)" % res.desc_names[0])
    ax.set_ylabel(ylabel)
    if ax.get_legend_handles_labels()[0]:
        ax.legend(frameon=False, fontsize=7.5, loc="best")
    _style(ax)
    fig.tight_layout()
    if path:
        fig.savefig(path, bbox_inches="tight")
    return fig, ax


def volcano_2d(res, step=-1, path=None, cmap="viridis", levels=24,
               clabel=r"$\log_{10}$ (rate / s$^{-1}$)", mark_peak=True,
               figsize=(5.4, 4.2), dpi=200):
    """2D 火山图（等值填充 + 等值线）。**配色默认 viridis，不用 jet**（理由见模块 docstring）。"""
    plt = _plt()
    if len(res.shape) != 2:
        raise ValueError("volcano_2d 只画 2D 扫描")
    X, Y = np.meshgrid(res.desc_axes[0], res.desc_axes[1], indexing="ij")
    r = np.abs(res.rate_of(step))
    ok = res.converged & np.isfinite(r) & (r > 0)
    with np.errstate(divide="ignore"):
        Z = np.where(ok, np.log10(np.where(r > 0, r, 1.0)), np.nan)
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    cf = ax.contourf(X, Y, Z, levels=levels, cmap=cmap)
    cs = ax.contour(X, Y, Z, levels=levels, colors="k", linewidths=0.25, alpha=0.5)
    cb = fig.colorbar(cf, ax=ax, pad=0.02)
    cb.set_label(clabel)
    if mark_peak and ok.any():
        ij = np.unravel_index(np.nanargmax(np.where(ok, Z, -np.inf)), Z.shape)
        ax.plot([X[ij]], [Y[ij]], "*", ms=13, color="#ffffff", mec="k", mew=0.7,
                zorder=6)
        ax.annotate("peak (%.2f, %.2f)" % (X[ij], Y[ij]), (X[ij], Y[ij]),
                    textcoords="offset points", xytext=(8, 6), fontsize=7.5,
                    color="w")
    if not ok.all():
        ax.contourf(X, Y, np.where(ok, np.nan, 1.0), levels=[0.5, 1.5],
                    colors="none", hatches=["xx"])
    ax.set_xlabel("%s (eV)" % res.desc_names[0])
    ax.set_ylabel("%s (eV)" % res.desc_names[1])
    _style(ax)
    fig.tight_layout()
    if path:
        fig.savefig(path, bbox_inches="tight")
    return fig, ax


def coverage_1d(res, path=None, figsize=(5.2, 3.8), dpi=200, logy=True):
    """覆盖度随描述符变化。覆盖度跨很多数量级，默认用对数纵轴。"""
    plt = _plt()
    x = res.desc_axes[0]
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    cols = plt.get_cmap("viridis")(np.linspace(0.08, 0.88, len(res.theta) + len(res.theta_vac)))
    for k, (nm, arr) in enumerate(list(res.theta.items()) + list(res.theta_vac.items())):
        ax.plot(x, arr, "-", lw=1.4, color=cols[k], label=nm)
    if logy:
        ax.set_yscale("log")
    ax.set_xlabel("%s (eV)" % res.desc_names[0])
    ax.set_ylabel(r"coverage $\theta$")
    ax.legend(frameon=False, fontsize=7.5, ncol=2)
    _style(ax)
    fig.tight_layout()
    if path:
        fig.savefig(path, bbox_inches="tight")
    return fig, ax


def drc_1d(desc, X_B, step_labels, path=None, figsize=(5.2, 3.8), dpi=200):
    """台阶 DRC 随描述符变化（堆叠面积图：Σ = 1，正好能用堆叠表现）。"""
    plt = _plt()
    X = np.asarray(X_B, float)
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    cols = plt.get_cmap("viridis")(np.linspace(0.05, 0.9, X.shape[1]))
    ax.stackplot(desc, np.clip(X.T, 0, None), colors=cols,
                 labels=[s[:24] for s in step_labels], ec="w", lw=0.3)
    ax.set_xlabel("descriptor (eV)")
    ax.set_ylabel(r"degree of rate control $X_{\mathrm{B}}$")
    ax.set_ylim(0, 1.02)
    ax.legend(frameon=False, fontsize=6.8, loc="upper center", ncol=2)
    ax.text(0.015, 0.97, r"$\sum X_{\mathrm{B}} = 1$ (Campbell sum rule)", transform=ax.transAxes,
            fontsize=7, va="top")
    _style(ax)
    fig.tight_layout()
    if path:
        fig.savefig(path, bbox_inches="tight")
    return fig, ax
