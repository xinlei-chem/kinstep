# -*- coding: utf-8 -*-
# KinSTEP - a modality-agnostic microkinetic solver.
# Copyright (C) 2026 The KinSTEP authors.
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.  It is distributed WITHOUT ANY WARRANTY; see the GNU General Public
# License (LICENSE) for details.
r"""[English] Descriptor scans: 1D curves and 2D maps, with continuation between points.

(The commentary below is in Chinese: it records why each choice was made
and which measurement forced it.)

描述符扫描（T9）：1D 曲线与 2D 火山图。

扫描不是"循环调一遍求解器"
--------------------------
★ **续算（continuation）是必需的，不是优化。** T7 实测：
  · 从最难的一端（E_N = −2，表面中毒、θ_vac ~ 1e-30）起步 ⇒ **第一个点就卡住**；
  · 从中间出发向两侧走、每点用邻点的解当初值 ⇒ 401 点 1.0 s，**牛顿迭代中位 0 次**。
  差别不是快慢，是**能不能算出来**。

所以本模块的扫描顺序是**从网格中心向外**：
  1D：从中点向左右各走一遍；
  2D：先在中心行按 1D 走通，再以该行为种子向上下逐行推进，行内再左右续算。
  —— 每个点的初值都来自**物理上相邻**的已解点。

★ 并行怎么处理（为什么默认不并行）
----------------------------------
续算天生是**串行依赖**的：下一个点要用上一个点的解。
可以按"行"并行（2D 时每行独立续算），但行与行之间就失去了种子 ⇒ 难点行容易失败。
实测 401 点只要 1.0 s（double）/ 4.0 s（含自动高精度），**并行收益远小于稳健性损失**。
⇒ 默认串行续算；`n_jobs>1` 时按行并行，并在结果里记下用了哪种（可复现）。

★ 输出一定带**可作图的原始数据**
`ScanResult.save_dat()` 直接写出 Origin/gnuplot 能读的列数据 ——
图是给人看的，数据是给后续分析和重画用的，两者都要给。
"""
from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np

from .analyze import energetics_from_arrays
from .driver import DriverStack
from .mech import Mechanism
from .scaling import ScalingModel
from .solve import SteadyResult, steady_state


class ScanError(RuntimeError):
    """扫描本身没法进行（坐标轴与描述符对不上、中心点解不出来等）。"""


@dataclass
class ScanResult:
    desc_axes: List[np.ndarray]        # 每个描述符的坐标轴
    desc_names: List[str]
    shape: Tuple[int, ...]
    theta: Dict[str, np.ndarray]       # 物种名 → 与 shape 同形
    theta_vac: Dict[str, np.ndarray]
    rate_net: np.ndarray               # shape + (n_step,)
    zrev: np.ndarray
    converged: np.ndarray              # bool，shape
    precision: np.ndarray              # 'double'/'high'，shape（object）
    n_iter: np.ndarray
    flux_fix: np.ndarray
    wall_s: float
    n_high: int
    order: str
    notes: List[str] = field(default_factory=list)

    @property
    def n_point(self):
        return int(np.prod(self.shape))

    def rate_of(self, step=-1):
        """取某一步的净速率（默认最后一步，通常是产物生成步）。"""
        return self.rate_net[..., step]

    def summary(self) -> str:
        L = ["扫描：%s，共 %d 点，耗时 %.2f s（%.1f ms/点）"
             % (" × ".join("%s[%d]" % (n, len(a))
                           for n, a in zip(self.desc_names, self.desc_axes)),
                self.n_point, self.wall_s, self.wall_s / max(self.n_point, 1) * 1e3)]
        L.append("  顺序：%s；收敛 %d/%d；升到高精度 %d 点"
                 % (self.order, int(self.converged.sum()), self.n_point, self.n_high))
        r = self.rate_of()
        ok = self.converged & np.isfinite(r) & (np.abs(r) > 0)
        if ok.any():
            i = np.unravel_index(np.nanargmax(np.where(ok, np.abs(r), -np.inf)), self.shape)
            peak = [float(self.desc_axes[k][i[k]]) for k in range(len(self.shape))]
            L.append("  速率峰：%s = %s，r = %.6e"
                     % (self.desc_names, ["%.3f" % v for v in peak], float(r[i])))
            L.append("  速率跨 %.1f 个数量级"
                     % (np.log10(np.abs(r[ok]).max() / np.abs(r[ok]).min())))
        if not self.converged.all():
            bad = np.argwhere(~self.converged)
            L.append("  ⚠ 未收敛 %d 点，首个在索引 %s" % (len(bad), tuple(bad[0])))
        return "\n".join(L)

    def save_dat(self, path, step=-1, extra_species=True):
        """写出列数据（Origin / gnuplot 直接读）。**图之外必须给数据。**"""
        axes = np.meshgrid(*self.desc_axes, indexing="ij")
        cols, head = [], []
        for k, nm in enumerate(self.desc_names):
            cols.append(axes[k].ravel()); head.append(nm)
        r = self.rate_of(step).ravel()
        cols.append(r); head.append("rate_net")
        with np.errstate(divide="ignore", invalid="ignore"):
            cols.append(np.log10(np.abs(r))); head.append("log10_abs_rate")
        if extra_species:
            for nm, arr in list(self.theta.items()) + list(self.theta_vac.items()):
                cols.append(arr.ravel()); head.append("theta_" + nm.replace("*", "s"))
        cols.append(self.converged.ravel().astype(int)); head.append("converged")
        dat = np.column_stack(cols)
        np.savetxt(path, dat, header="  ".join(head), fmt="%.10e", comments="# ")
        return path


def _solve_one(mech, sc, desc, state, driver, p_gas, c_liq, guess, tol,
               precision, dps, A_pre=None, omz_floor=None):
    # ★ 2026-10-22（T9 抓到）：初版这里写死 `A_pre=None`，**让 scan() 的 A_pre 参数形同虚设**。
    #   症状极其隐蔽：能量层对到 7.1e-15，可速率整体差 2.249e-08 ——
    #   因为速率 ∝ 前因子，而默认前因子 k_BT/h 与 CATKINAS 的老常数差 2e-8。
    #   "签名有参数但实现忽略它" 是最难查的一类错：不报错、量级又小得像数值噪声。
    Ga, dG = sc.energies(desc)
    en = energetics_from_arrays(mech, Ga, dG, A_pre=A_pre)
    kw = dict(state=state, driver=driver, p_gas=p_gas, c_liq=c_liq, tol=tol,
              precision=precision, dps=dps, omz_floor=omz_floor,
              strict_energy=False)
    for lad, g in ((("newton",), guess), (("ode", "newton"), None)):
        if lad == ("newton",) and g is None:
            continue
        try:
            r = steady_state(mech, en, theta_guess=g, ladder=lad,
                             ode_stages=4, ode_budget_s=6.0, **kw)
            if r.converged:
                return r
        except Exception:
            pass
    return None


def scan(mech: Mechanism, sc: ScalingModel, state: Dict,
         desc_axes: Sequence[Sequence[float]], driver: DriverStack = None,
         p_gas: Dict = None, c_liq: Dict = None, A_pre: float = None,
         tol: float = 1e-12, precision: str = "auto", dps: int = 50,
         omz_floor: float = None,
         verbose: bool = False) -> ScanResult:
    """1D 或 2D 描述符扫描。`desc_axes` 给每个描述符的取值序列。

    返回 `ScanResult`。顺序一律**从网格中心向外续算**（理由见模块 docstring）。
    """
    axes = [np.asarray(a, float) for a in desc_axes]
    if len(axes) != sc.n_desc:
        raise ScanError("scaling 模型有 %d 个描述符 %s，但给了 %d 个坐标轴"
                        % (sc.n_desc, sc.desc_names, len(axes)))
    if len(axes) > 2:
        raise ScanError("目前支持 1D 与 2D 扫描（%d 维请自己外层循环）" % len(axes))
    shape = tuple(len(a) for a in axes)
    surf = mech.surface_species
    vacs = [mech.vacancy_of(s) for s in mech.sites if mech.vacancy_of(s) in mech.species]

    TH = {n: np.full(shape, np.nan) for n in surf}
    VA = {v: np.full(shape, np.nan) for v in vacs}
    RN = np.full(shape + (mech.n_step,), np.nan)
    ZZ = np.full(shape + (mech.n_step,), np.nan)
    CV = np.zeros(shape, bool)
    PR = np.empty(shape, dtype=object); PR[...] = "—"
    NI = np.zeros(shape, int)
    FX = np.full(shape, np.nan)
    n_high = 0

    def store(ix, r):
        nonlocal n_high
        for n in surf:
            TH[n][ix] = r.theta[n]
        for v in vacs:
            VA[v][ix] = r.theta_vac[v]
        RN[ix] = r.rate_net; ZZ[ix] = r.zrev
        CV[ix] = True; PR[ix] = r.precision; NI[ix] = r.n_iter
        FX[ix] = r.flux_fix
        n_high += (r.precision == "high")

    t0 = time.time()
    if len(axes) == 1:
        mid = shape[0] // 2
        r0 = _solve_one(mech, sc, axes[0][mid], state, driver, p_gas, c_liq,
                        None, tol, precision, dps, A_pre, omz_floor)
        if r0 is None:
            raise ScanError("网格中心点就解不出来 —— 先查 scaling 系数与工况")
        store((mid,), r0)
        for d in (+1, -1):
            g = dict(r0.theta); i = mid + d
            while 0 <= i < shape[0]:
                r = _solve_one(mech, sc, axes[0][i], state, driver, p_gas, c_liq,
                               g, tol, precision, dps, A_pre, omz_floor)
                if r is None:
                    g = None
                else:
                    store((i,), r); g = dict(r.theta)
                i += d
        order = "1D 由中心 (%d) 向两侧续算" % mid
    else:
        mi, mj = shape[0] // 2, shape[1] // 2
        r0 = _solve_one(mech, sc, [axes[0][mi], axes[1][mj]], state, driver,
                        p_gas, c_liq, None, tol, precision, dps, A_pre, omz_floor)
        if r0 is None:
            raise ScanError("网格中心点就解不出来 —— 先查 scaling 系数与工况")
        store((mi, mj), r0)
        # 中心行先走通（它是后面所有行的种子）
        row_seed = {mj: dict(r0.theta)}
        for d in (+1, -1):
            g = dict(r0.theta); j = mj + d
            while 0 <= j < shape[1]:
                r = _solve_one(mech, sc, [axes[0][mi], axes[1][j]], state, driver,
                               p_gas, c_liq, g, tol, precision, dps, A_pre, omz_floor)
                if r is None:
                    g = None
                else:
                    store((mi, j), r); g = dict(r.theta); row_seed[j] = dict(r.theta)
                j += d
        # 以中心行为种子，向上下逐行推进；行内再左右续算
        for di in (+1, -1):
            prev = dict(row_seed)
            i = mi + di
            while 0 <= i < shape[0]:
                cur = {}
                for d in (0, +1, -1):
                    js = [mj] if d == 0 else range(mj + d, shape[1] if d > 0 else -1, d)
                    g = prev.get(mj)
                    for j in js:
                        seed = cur.get(j - d) if d else None
                        r = _solve_one(mech, sc, [axes[0][i], axes[1][j]], state,
                                       driver, p_gas, c_liq,
                                       seed or prev.get(j) or g, tol, precision, dps,
                                       A_pre)
                        if r is not None:
                            store((i, j), r); cur[j] = dict(r.theta)
                if cur:
                    prev = cur
                i += di
        order = "2D 由中心 (%d,%d) 先走中心行，再逐行向外；行内左右续算" % (mi, mj)
    wall = time.time() - t0

    res = ScanResult(desc_axes=axes, desc_names=list(sc.desc_names), shape=shape,
                     theta=TH, theta_vac=VA, rate_net=RN, zrev=ZZ, converged=CV,
                     precision=PR, n_iter=NI, flux_fix=FX, wall_s=wall,
                     n_high=n_high, order=order,
                     notes=["续算是必需的不是优化：T7 实测从最难一端起步第一个点就卡住",
                            "precision='%s'，高精度位数 %d" % (precision, dps)])
    if verbose:
        print(res.summary())
    return res
