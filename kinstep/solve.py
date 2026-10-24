# -*- coding: utf-8 -*-
# KinSTEP - a modality-agnostic microkinetic solver.
# Copyright (C) 2026 The KinSTEP authors.
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.  It is distributed WITHOUT ANY WARRANTY; see the GNU General Public
# License (LICENSE) for details.
r"""[English] The solver: transient integration and the steady-state ladder (ODE guess, damped
Newton in ln-theta with an analytic Jacobian, flux projection, arbitrary
precision).

(The commentary below is in Chinese: it records why each choice was made
and which measurement forced it.)

求解器 —— 瞬态（T5）。稳态在 T6。

瞬态要解的就是一组常微分方程
-----------------------------
    dθ_i/dt = Σ_s ν[s, i] · r_net[s](θ, 工况(t))          i ∈ 表面物种

空位**不进方程组**，而是由位点守恒现算：θ_vac = 1 − Σθ（按位点分别算）。
这样做的两个理由：
  ① 少一个未知数，而且**守恒恒成立**（不会因积分误差漂移）；
  ② 机理解析时的位点守恒闸（闸③）已保证 Σν 对每种位点为 0，
     所以 d(Σθ + θ_vac)/dt ≡ 0 是**恒等式而不是需要监控的量**。

气相/液相活度默认**保持恒定**（微分反应器 / 开放体系）。
批次反应器（气相随反应消耗）留到 T6 之后，接口已经留好（`consume_gas`）。

★ 工况随时间变化：Protocol
--------------------------
电化学的两个经典实验都是"工况随时间走"：

    CVSweep        三角波扫电势（循环伏安）
    PotentialStep  电势阶跃
    Constant       恒定工况（热催化、恒电势）

`Protocol.state_at(t)` 返回那一刻的 state 字典，求解器据此**每一步重算速率常数**。
`is_constant=True` 的协议只算一次 —— 否则纯热催化会白算几千次 exp。

★ 为什么用隐式方法（Radau/BDF）而不是 RK45
-------------------------------------------
表面覆盖度方程天生**刚性**：准平衡的吸附步 k ~ 10¹⁰ s⁻¹，而反应步 ~ 10⁰ s⁻¹，
时间尺度差 10 个数量级。显式方法会被最快的那一步逼到极小步长而几乎不前进。
（Koper 的 notebook 也用 Radau，这一点上照着来。）
"""
from __future__ import annotations

import abc
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np

from .driver import DriverStack
from .i18n import t as _t
from .mech import Mechanism
from .rates import (CompiledRates, Energetics, RateConstants, activities,
                    apply_state, compile_rates, rates_factored)
from .units import K_B_EV


class SolverError(RuntimeError):
    """求解失败。**不返回半成品**，宁可报错。"""


# ══════════════════════════════ 工况协议 ══════════════════════════════
class Protocol(abc.ABC):
    """工况随时间的变化。`state_at(t)` → state 字典。"""

    #: True ⇒ 工况不随时间变 ⇒ 速率常数只算一次
    is_constant = False

    def __init__(self, **base):
        #: 不随时间变的那部分工况（T_K、pH_val、gamma_site …）
        self.base = dict(base)

    @abc.abstractmethod
    def state_at(self, t: float) -> Dict:
        raise NotImplementedError

    def describe(self) -> str:
        return self.__class__.__name__


class Constant(Protocol):
    """恒定工况。"""
    is_constant = True

    def state_at(self, t):
        return dict(self.base)

    def describe(self):
        return "恒定工况 %s" % self.base


class CVSweep(Protocol):
    """三角波扫电势（循环伏安）。

    参数
    ----
    U_start, U_end : V，扫描的两端
    sweep_rate     : V/s，**正数**（方向由 start→end 决定）
    n_cycle        : 扫几圈（默认 1 圈 = 正扫 + 回扫）

    一圈的时长 = 2·|U_end − U_start| / sweep_rate。
    与 Koper notebook 的 `E_CV(t)` 等价，但用 `np.interp` 写，
    **避免 notebook 里 `np.piecewise` + `and` 那种只能传标量的写法**。
    """

    def __init__(self, U_start, U_end, sweep_rate, n_cycle=1, **base):
        super().__init__(**base)
        self.U_start, self.U_end = float(U_start), float(U_end)
        self.sweep_rate = abs(float(sweep_rate))
        self.n_cycle = int(n_cycle)
        if self.sweep_rate <= 0:
            raise ValueError("sweep_rate 必须 > 0")
        self.half = abs(self.U_end - self.U_start) / self.sweep_rate
        self.cycle_time = 2.0 * self.half
        self.t_total = self.cycle_time * self.n_cycle

    def U_at(self, t):
        """标量或数组都接受。"""
        t = np.asarray(t, float)
        tau = np.mod(t, self.cycle_time)
        up = tau <= self.half
        U = np.where(up,
                     self.U_start + np.sign(self.U_end - self.U_start) * tau * self.sweep_rate,
                     self.U_end - np.sign(self.U_end - self.U_start)
                     * (tau - self.half) * self.sweep_rate)
        return U if U.ndim else float(U)

    def state_at(self, t):
        st = dict(self.base)
        st["U_elec"] = self.U_at(t)
        return st

    def describe(self):
        return ("CV 扫描 %.3f→%.3f V @ %.1f mV/s，%d 圈（一圈 %.2f s，总 %.2f s）"
                % (self.U_start, self.U_end, self.sweep_rate * 1e3,
                   self.n_cycle, self.cycle_time, self.t_total))


class PotentialStep(Protocol):
    """电势阶跃：t < t_step 时 U = U_before，之后 U = U_after。"""

    def __init__(self, U_before, U_after, t_step=0.0, **base):
        super().__init__(**base)
        self.U_before, self.U_after = float(U_before), float(U_after)
        self.t_step = float(t_step)

    def state_at(self, t):
        st = dict(self.base)
        st["U_elec"] = self.U_before if t < self.t_step else self.U_after
        return st

    def describe(self):
        return "电势阶跃 %.3f → %.3f V（t = %.3g s 处）" % (
            self.U_before, self.U_after, self.t_step)


class Ramp(Protocol):
    """任意单调升温/升压 —— 给 TPD / TPR 留的口子（物理同 CV，只是换了变量）。"""

    def __init__(self, key, start, rate, stop=None, **base):
        super().__init__(**base)
        self.key, self.start, self.rate, self.stop = key, float(start), float(rate), stop

    def state_at(self, t):
        st = dict(self.base)
        v = self.start + self.rate * t
        st[self.key] = v if self.stop is None else min(v, self.stop)
        return st

    def describe(self):
        return "%s 从 %g 以 %g /s 线性变化" % (self.key, self.start, self.rate)


class Custom(Protocol):
    """自己给一个 t → state 的函数（别的都不够用时）。"""

    def __init__(self, fn: Callable[[float], Dict], **base):
        super().__init__(**base)
        self.fn = fn

    def state_at(self, t):
        st = dict(self.base)
        st.update(self.fn(t))
        return st


# ══════════════════════════════ 流动反应器 ══════════════════════════════
@dataclass
class FlowReactor:
    r"""连续流动（CSTR）：气相组成**也是未知数**，不再是固定的工况。

    微分反应器（本包默认）假定气相活度恒定 —— 转化率极低时成立。
    一旦要算**转化率随温度/压力怎么变**（工业条件、积分反应器），
    就必须让气相随反应演化：

        dP_i/dt = Pn·( Σ_s ν[s,i]·r_s + k_flow·(P_in,i − P_i) )
        k_flow  = Fgt·P0 / (N_A·T·k_B)

    `Pn` 把"每位点每秒"的表面速率换算成压力变化率；`k_flow` 是归一化的进出料速率。
    （形式与 CATKINAS 的 `odefun_map.m` 逐项一致，已逐数字回归。）

    ★ 为什么能直接接进现有求解器
    ----------------------------
    **ln p 进速率方程的方式与 ln θ 完全相同**（都是 `Σ a·ln(活度)` 里的一项），
    所以解析雅可比的那条公式
        ∂r_net[s]/∂x_j = R_f[s]·[a_r[s,j]·(1−z) − z·(a_p−a_r)[s,j]]
    对气相列**原样成立** —— 只需把气相加进未知数向量，再补上各自的物料平衡方程。
    不必另写一套求解器。

    参数
    ----
    p_in     : {气相物种: 进料分压}（无量纲，P/P0）
    Fgt      : 归一化体积流量
    P0_Pa    : 参比压力（默认 1e5 Pa）
    Pn       : 表面速率 → 压力变化率的换算系数
    N_A, k_B : 可显式指定，用于**与别的程序对数**（CATKINAS 用的是 1986 年的值）
    """
    p_in: Dict[str, float]
    Fgt: float
    Pn: float = 1.0
    P0_Pa: float = 1.0e5
    N_A: float = None
    k_B: float = None
    #: ★ 直接给 k_flow（1/s），绕过公式。**与别的程序对数时用这个。**
    #   实测教训：公式里的 T 该用**名义温度**，而求解器传进来的是「等效温度」
    #   （为统一 k_BT 而微调过的那个，差 1.6e-6）。两者混用会让出口分压
    #   带上一个**恒定的 1.13e-6 相对偏差** —— 201 个点上一模一样，
    #   是常数而非噪声，典型的口径问题。
    k_flow_fixed: float = None

    def k_flow(self, T_K):
        """归一化进出料速率（1/s）。给了 `k_flow_fixed` 就直接用它。"""
        if self.k_flow_fixed is not None:
            return float(self.k_flow_fixed)
        from .units import N_A as _NA, K_B as _KB
        na = self.N_A if self.N_A else _NA
        kb = self.k_B if self.k_B else _KB
        return self.Fgt * self.P0_Pa / (na * float(T_K) * kb)

    def describe(self, T_K=None):
        t = "" if T_K is None else "，k_flow(%.0f K) = %.6g" % (T_K, self.k_flow(T_K))
        return ("CSTR：进料 %s，Fgt = %g，Pn = %g%s"
                % ({k: "%g" % v for k, v in self.p_in.items()}, self.Fgt, self.Pn, t))


# ══════════════════════════════ 结果容器 ══════════════════════════════
@dataclass
class TransientResult:
    t: np.ndarray                       # (n_t,)
    theta: Dict[str, np.ndarray]        # 表面物种 → (n_t,)
    theta_vac: Dict[str, np.ndarray]    # 空位名 → (n_t,)
    rate_net: np.ndarray                # (n_t, n_step)
    zrev: np.ndarray                    # (n_t, n_step)
    obs: Dict[str, np.ndarray]          # driver 的可观测量 → (n_t,)
    state_trace: Dict[str, np.ndarray]  # 随时间变的工况 → (n_t,)
    n_rhs: int = 0                      # 右端函数调用次数（看效率）
    notes: List[str] = field(default_factory=list)

    def check_bounds(self, tol=1e-6):
        """★ 免费自检：所有覆盖度与空位都该落在 [0,1]。

        越界说明**机理或能量有问题**（或步长太大），不是"数值小毛病"。
        返回 (ok, 报告)。
        """
        bad = []
        for name, arr in list(self.theta.items()) + list(self.theta_vac.items()):
            lo, hi = float(np.min(arr)), float(np.max(arr))
            if lo < -tol or hi > 1.0 + tol:
                bad.append("%s ∈ [%.3e, %.6f]" % (name, lo, hi))
        if not bad:
            return True, "覆盖度与空位全程落在 [0,1]（容差 %g）✓" % tol
        return False, ("⚠ 覆盖度越界：%s\n  ⇒ 多半是位点守恒没写对，或某一步的能量离谱。"
                       % "；".join(bad))


# ══════════════════════════════ 瞬态求解 ══════════════════════════════
def transient(mech: Mechanism, en: Energetics, theta0: Dict[str, float],
              t_span, protocol: Protocol = None, driver: DriverStack = None,
              p_gas: Dict[str, float] = None, c_liq: Dict[str, float] = None,
              n_eval: int = 200, t_eval=None, method: str = "Radau",
              rtol: float = 1e-8, atol: float = 1e-12,
              strict_energy: bool = True) -> TransientResult:
    """积分 dθ/dt = νᵀ·r_net。

    参数
    ----
    theta0 : {表面物种: 初始覆盖度}；没写的按 0
    t_span : (t0, t1) 或 一个数（= (0, 该数)）
    protocol : 工况协议；不给就用 `Constant(T_K=298.15)`
    method : 'Radau'(默认) / 'BDF' / 'LSODA' —— **都是隐式的**，原因见模块 docstring

    ★ 温度从哪来
    ------------
    一律从 `protocol.base["T_K"]` 或 `state_at(t)["T_K"]` 取；**本函数不设温度默认**，
    漏给就报错 —— 速率常数对温度是指数敏感的，静默取 298.15 K 太危险。
    """
    from scipy.integrate import solve_ivp

    protocol = protocol or Constant()
    surf = mech.surface_species
    if not surf:
        raise SolverError("机理里没有表面物种，没有覆盖度可积分。")
    unknown = set(theta0 or {}) - set(surf)
    if unknown:
        raise SolverError("theta0 里这些名字不是表面物种：%s\n  表面物种只有：%s"
                          % (sorted(unknown), surf))
    y0 = np.array([float((theta0 or {}).get(n, 0.0)) for n in surf])

    t0, t1 = (0.0, float(t_span)) if np.isscalar(t_span) else (float(t_span[0]), float(t_span[1]))
    if t_eval is None:
        t_eval = np.linspace(t0, t1, int(n_eval))

    # 位点容量自检：初值就超了，后面不可能对
    for site in mech.sites:
        s = sum(y0[i] for i, n in enumerate(surf) if mech.species[n].site == site)
        if s > 1.0 + 1e-12:
            raise SolverError("位点 %s 的初始覆盖度之和 = %.6f > 1 —— 位点装不下"
                              % (("默认" if site == "" else site), s))

    idx_surf = [mech.species_order.index(n) for n in surf]
    nu_surf = np.asarray(mech.nu, float)[:, idx_surf]        # (n_step, n_surf)

    def _T_of(st):
        if "T_K" not in st or st["T_K"] is None:
            raise SolverError(
                "工况里没有温度 T_K。速率常数对温度是**指数敏感**的，本包不设默认值。\n"
                "给法： Constant(T_K=300.0)  或  CVSweep(..., T_K=300.0)")
        return float(st["T_K"])

    # ★★ 两层分离：来源闸 + Val 取值**只做这一次**（compile_rates）。
    #    初版在 ODE 右端每次都调 rate_constants() ⇒ 每次都跑一遍来源闸，
    #    Koper 的 16 s CV **跑了 2 分 19 秒还没出**。详见 CompiledRates 的说明。
    cr = compile_rates(en, strict_energy=strict_energy)
    st0 = protocol.state_at(t0)
    T0 = _T_of(st0)

    # 工况恒定 ⇒ 速率常数连 exp 都只算一次
    cache = {"rc": apply_state(cr, T0, driver, st0)} if protocol.is_constant else {}
    counter = {"n": 0}

    def _rc_at(t):
        st = protocol.state_at(t)
        if "rc" in cache:
            return cache["rc"], st
        return apply_state(cr, _T_of(st), driver, st), st

    # 内层只用到这几个局部名，避免每次都走属性查找
    _af, _ar, _ap = activities, cr.a_r, cr.a_p
    _nuT = nu_surf.T

    def rhs(t, y):
        counter["n"] += 1
        rc, _ = _rc_at(t)
        act = _af(mech, dict(zip(surf, y)), p_gas, c_liq)
        _, _, rn = rates_factored(mech, rc, act, _ar, _ap)
        return _nuT @ rn

    sol = solve_ivp(rhs, (t0, t1), y0, t_eval=t_eval, method=method,
                    rtol=rtol, atol=atol)
    if not sol.success:
        raise SolverError(
            "积分失败：%s\n"
            "常见原因：① 某步能量离谱导致刚性极强（先用 rate_constants().table() 看 kf/kr 跨几个数量级）；"
            "② 初始覆盖度之和 > 1；③ t_span 比最快的时间尺度还长太多。" % sol.message)

    # ── 逐时间点回算速率、空位与可观测量 ──
    theta = {n: sol.y[i] for i, n in enumerate(surf)}
    n_t = len(sol.t)
    rate_net = np.zeros((n_t, mech.n_step))
    zrev = np.zeros((n_t, mech.n_step))
    vac = {mech.vacancy_of(s): np.zeros(n_t) for s in mech.sites
           if mech.vacancy_of(s) in mech.species}
    obs: Dict[str, List[float]] = {}
    trace: Dict[str, List[float]] = {}
    for k, tk in enumerate(sol.t):
        th = {n: theta[n][k] for n in surf}
        rc, st = _rc_at(tk)
        act = activities(mech, th, p_gas, c_liq)
        _, zk, rk = rates_factored(mech, rc, act, cr.a_r, cr.a_p)
        rate_net[k], zrev[k] = rk, zk
        for site in mech.sites:
            vn = mech.vacancy_of(site)
            if vn in vac:
                vac[vn][k] = act[mech.species_order.index(vn)]
        if driver is not None:
            for kk, vv in driver.observables(mech, {"rate_net": rk}, st).items():
                obs.setdefault(kk, []).append(float(vv))
        for kk, vv in st.items():
            if isinstance(vv, (int, float, np.floating)):
                trace.setdefault(kk, []).append(float(vv))

    res = TransientResult(
        t=sol.t, theta=theta, theta_vac=vac, rate_net=rate_net, zrev=zrev,
        obs={k: np.asarray(v) for k, v in obs.items()},
        state_trace={k: np.asarray(v) for k, v in trace.items()
                     if np.ptp(v) > 0},            # 只留**真的在变**的工况
        n_rhs=counter["n"],
        notes=cr.notes + ["协议：" + protocol.describe(),
               "方法：%s（隐式，表面覆盖度方程天生刚性）" % method,
               "右端函数调用 %d 次，输出 %d 个时间点" % (counter["n"], n_t)])
    ok, rep = res.check_bounds()
    res.notes.append(rep)
    return res


# ══════════════════════════════ 稳态求解（T6）══════════════════════════════
r"""稳态部分的设计说明（T4 的测量直接决定了这里的变量选择）

未知数为什么取 ln θ 而不是 θ
----------------------------
① θ 跨 20 个数量级（10⁻²⁰ 的空位覆盖度是常事），线性变量下牛顿步长没法同时照顾两端；
② **取 exp 保证 θ > 0 恒成立**，不必在迭代里加边界约束（约束一加，雅可比就不干净了）；
③ 质量作用定律在 log 空间是**线性**的：ln r_fwd = ln kf + a_r·ln a，
   所以雅可比有干净的解析式（见下），不用数值差分 ——
   这一点很要紧：1e-8 的差分噪声本身就会摧毁我们想保住的那点精度。

空位：当方程，不当代入
----------------------
θ_vac 也作为未知数（ln θ_vac），位点守恒 Σθ + θ_vac = 1 作为**一条方程**。
这样未知数 = n_surf + n_site，方程 = n_surf 个稳态 + n_site 个守恒，恰好配平。

空位的稳态方程**不需要**单独写：机理闸③保证每一步对每种位点 Σν = 0，
所以 Σ_{含空位} dθ/dt ≡ 0 —— 吸附态的稳态方程都满足时，空位自动稳态。
（这是恒等式而非近似，所以"用守恒替掉空位方程"是精确的。）

★ (1−z) 怎么算（T4/T5 两轮测量的结论）
--------------------------------------
    ln z = ΔG/k_BT + (a_p − a_r)·ln a      ← log 空间里是线性的
    1 − z = −expm1(ln z)                   ← **必须 expm1**

实测：1−z 真值 1e-12 时 `1 − exp(ln z)` 相对误差 2.2e-05，而 `−expm1` 是 0。
剩下的瓶颈是 ln z 自身的绝对精度（~1e-15）—— 它是 O(10) 量级项之和。
⇒ **本求解器能可靠分辨到 1−z ≈ 1e-10；更接近平衡的步会被诊断标出来**，
  不假装算得准（真要更准需要把 ln z 本身当未知数或上 mpmath，见 notes）。

解析雅可比
----------
    r_net[s] = R_f[s]·(1−z_s),   R_f = pref·exp(ln kf + a_r·ln a)
    ∂R_f/∂x_j     = R_f·a_r[s,j]
    ∂(1−z)/∂x_j   = −z_s·(a_p − a_r)[s,j]
    ⇒ ∂r_net[s]/∂x_j = R_f[s]·[ a_r[s,j]·(1−z_s) − z_s·(a_p−a_r)[s,j] ]
    ⇒ ∂F_i/∂x_j = Σ_s ν[s,i]·∂r_net[s]/∂x_j
守恒方程：∂g_site/∂x_i = θ_i（该位点上的物种与空位），其余为 0。
"""


@dataclass
class SteadyResult:
    theta: Dict[str, float]
    theta_vac: Dict[str, float]
    rate_fwd: np.ndarray
    zrev: np.ndarray
    rate_net: np.ndarray          # ★ 通量投影后的（该用这个）
    rate_net_raw: np.ndarray      # R_f·(1−z) 直接算的（准平衡步上不可靠）
    flux_fix: float               # 投影改动了多少（相对），= 原始解的通量不自洽程度
    ln_z: np.ndarray
    residual: np.ndarray          # 缩放后的稳态残差（**仅表面物种那几项**）
    resid_abs: np.ndarray         # 未缩放（物理量纲）
    site_error: Dict[str, float]  # 位点守恒误差
    n_iter: int
    method: str
    converged: bool
    precision: str = "double"     # 'double' | 'high'
    hp: object = None             # 升到高精度时，这里挂 hp.HPResult
    #: 流动反应器时的**出口气相分压**（气相是未知数时才有）
    p_gas_out: Dict[str, float] = field(default_factory=dict)
    #: ★ 完整残差向量（表面 + 位点守恒 + 气相物料平衡）。
    #  收敛判据看的是**这个**，而 `residual` 只有表面那一段 ——
    #  T9b 诊断时被这个盲点坑过：打印出 1e-15 的残差却判未收敛，
    #  真正大的是气相那几项，而它们**根本没被存下来**。
    residual_full: np.ndarray = None
    obs: Dict[str, float] = field(default_factory=dict)
    notes: List[str] = field(default_factory=list)

    def report(self) -> str:
        from .i18n import t
        conv = t("收敛", "converged") if self.converged \
            else t("★未收敛", "*** NOT converged")
        L = [t("稳态解（%s，%d 次迭代，%s）" % (self.method, self.n_iter, conv),
               "steady state (%s, %d iterations, %s)"
               % (self.method, self.n_iter, conv))]
        L.append("  " + t("覆盖度：", "coverages:  ")
                 + "  ".join("%s=%.6g" % (k, v) for k, v in self.theta.items()))
        L.append("  " + t("空位　：", "vacancies:  ")
                 + "  ".join("%s=%.6g" % (k, v) for k, v in self.theta_vac.items()))
        L.append("  " + t("位点守恒误差：", "site-balance error: ")
                 + "  ".join("%s:%.2e" % (k, v) for k, v in self.site_error.items()))
        L.append("  " + t("缩放残差 ‖F‖∞ = %.3e" % np.max(np.abs(self.residual)),
                          "scaled residual ||F||_inf = %.3e"
                          % np.max(np.abs(self.residual))))
        L.append("  " + t("逐步净速率(通量投影后)： ",
                          "net rate per step (after flux projection): ")
                 + "  ".join("r%d=%.6g" % (i + 1, r) for i, r in enumerate(self.rate_net)))
        L.append("  " + t("逐步 1−z 与可达精度：",
                          "per-step 1-z and attainable precision:"))
        for sid, omz, rel, dig in self.precision_table():
            note = t("  ← 由通量守恒推定", "  <- inferred from flux conservation") \
                if rel > 1e-6 else ""
            L.append("     " + t("步%d  1−z=%-11.4g  自身可达相对精度 %-9.1e (~%d 位)%s"
                                 % (sid, omz, rel, dig, note),
                                 "step %d  1-z=%-11.4g  attainable rel. precision "
                                 "%-9.1e (~%d digits)%s" % (sid, omz, rel, dig, note)))
        L.append("  " + t("通量投影改动量 = %.2e（即原始解的通量不自洽程度）" % self.flux_fix,
                          "flux-projection correction = %.2e (i.e. how far the raw "
                          "solution was from flux consistency)" % self.flux_fix))
        return "\n".join(L)

    def precision_table(self, eps=1e-16):
        """每一步净速率的**可达相对精度** ≈ eps/(1−z)。永远算出来，不设开关。

        为什么要永远报
        --------------
        初版只在 `0 < 1−z < 1e-10` 时才告警，结果两头都漏：
          · 合成氨算例里第 2 步 1−z = 2.3e-09，刚好在阈值之上 ⇒ **一声不响**，
            而实测 H 元素通量就差了 1.55e-07（预期 4.4e-08，同量级）；
          · 1−z **恰好下溢成 0** 时（最严重的情形）`0 < 1−z` 为假 ⇒ 也不告警。
        ⇒ 阈值式告警在这件事上靠不住；把数字直接摆出来才诚实。
        """
        out = []
        for i, z in enumerate(self.zrev):
            omz = 1.0 - z
            if omz >= 1.0:                      # 不可逆步，没有相消问题
                out.append((i + 1, omz, eps, 16))
                continue
            rel = (eps / omz) if omz > 0 else np.inf
            digits = 0 if not np.isfinite(rel) else max(0, int(-np.log10(max(rel, 1e-17))))
            out.append((i + 1, omz, rel, digits))
        return out

    def precision_warnings(self, want_rel=1e-6):
        """净速率精度达不到 `want_rel` 的步。**通量投影之后这些步已被修正**，
        这里保留是为了让用户知道"哪一步的速率是靠守恒推出来的，不是自己算出来的"。"""
        from .i18n import t
        out = []
        for sid, omz, rel, dig in self.precision_table():
            if rel > want_rel:
                if omz <= 0:
                    how = t("1−z 下溢为 0，净速率完全无法由 (1−z) 得到",
                            "1-z underflowed to 0, so the net rate cannot be obtained "
                            "from (1-z) at all")
                else:
                    how = t("1−z = %.3e ⇒ 自身只能算到 %.1e 的相对精度（约 %d 位有效数字）"
                            % (omz, rel, dig),
                            "1-z = %.3e, so on its own it reaches only %.1e relative "
                            "precision (~%d significant digits)" % (omz, rel, dig))
                out.append(t("第%d步 已准平衡：%s；该步净速率已由**通量守恒**从其它步推定。"
                             % (sid, how),
                             "step %d is quasi-equilibrated: %s; its net rate was "
                             "inferred from the other steps by flux conservation."
                             % (sid, how)))
        return out

    def flux_check(self, mech):
        """★ 免费自检：稳态下 νᵀ·r 必须为 0（逐表面物种）。返回 (ok, 最大相对不自洽)。"""
        idx = [mech.species_order.index(n) for n in mech.surface_species]
        A = np.asarray(mech.nu, float)[:, idx].T
        scale = np.abs(A) @ np.abs(self.rate_net)
        resid = A @ self.rate_net
        rel = np.max(np.abs(resid) / np.where(scale > 0, scale, 1.0))
        return rel < 1e-12, float(rel)


def _steady_residual(cr, rc, mech, surf, vac_names, site_of, x, row_scale,
                     flow=None, gas_names=(), T_K=None):
    """返回 (F, J, aux)。

    未知数排列：[ln θ_surf …, ln θ_vac …, (ln p_gas …)]
    方程排列　：[n_surf 个稳态, n_site 个位点守恒, (n_gas 个物料平衡)]
    气相那一段只在给了 `flow`（流动反应器）时才有。
    """
    n_s, n_v = len(surf), len(vac_names)
    n_g = len(gas_names) if flow is not None else 0
    names = surf + vac_names + (list(gas_names) if flow is not None else [])
    th = np.exp(x)
    # 组装活度（气/液部分由 cr 外部已固定，放在 _fixed_ln_a 里）
    ln_a = cr._ln_a_fixed.copy()
    for k, nm in enumerate(names):
        ln_a[cr._idx_of[nm]] = x[k]
    kT = K_B_EV * rc.T_K
    a_r, a_p = cr.a_r, cr.a_p
    mask_r = a_r > 0
    mask_p = a_p > 0
    sum_r = (a_r * np.where(mask_r, ln_a[None, :], 0.0)).sum(axis=1)
    sum_p = (a_p * np.where(mask_p, ln_a[None, :], 0.0)).sum(axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        ln_f = np.log(np.where(rc.kf > 0.0, rc.kf, np.nan)) + sum_r
    ln_z = rc.dG_used / kT + sum_p - sum_r
    ln_z = np.where(rc.kr > 0.0, ln_z, -np.inf)
    R_f = rc.prefactor * np.where(np.isfinite(ln_f), np.exp(ln_f), 0.0)
    lz = np.clip(ln_z, -700.0, 700.0)
    z = np.where(np.isfinite(ln_z), np.exp(lz), 0.0)
    omz = np.where(np.isfinite(ln_z), -np.expm1(lz), 1.0)      # ★ expm1
    r_net = R_f * omz

    nu_cols = np.array([cr._idx_of[nm] for nm in names])
    nu_sub = np.asarray(mech.nu, float)[:, nu_cols]            # (n_step, n_unknown)
    F_ss = nu_sub[:, :n_s].T @ r_net                            # 只对吸附态写稳态方程
    # 位点守恒
    F_cons = np.zeros(n_v)
    for v, vn in enumerate(vac_names):
        site = site_of[vn]
        ssum = th[n_s + v] + sum(th[k] for k, nm in enumerate(surf)
                                 if mech.species[nm].site == site)
        F_cons[v] = ssum - 1.0
    parts = [F_ss / row_scale, F_cons]
    if n_g:
        kfl = flow.k_flow(T_K)
        prod = nu_sub[:, n_s + n_v:].T @ r_net                  # Σ ν·r（逐气相）
        p_now = th[n_s + n_v:]
        p_in = np.array([float(flow.p_in.get(nm, flow.p_in.get(
            mech.species[nm].formula, 0.0))) for nm in gas_names])
        F_gas = prod + kfl * (p_in - p_now)
        gs = np.maximum(np.abs(prod) + kfl * np.maximum(p_in, p_now), 1e-300)
        parts.append(F_gas / gs)
    F = np.concatenate(parts)

    # ── 解析雅可比 ──
    # ★ 同一条公式同时覆盖表面与气相列（ln p 与 ln θ 在速率律里地位相同）
    dnet = R_f[:, None] * (a_r[:, nu_cols] * omz[:, None]
                           - z[:, None] * (a_p[:, nu_cols] - a_r[:, nu_cols]))
    n_x = n_s + n_v + n_g
    J_ss = nu_sub[:, :n_s].T @ dnet                             # (n_s, n_x)
    J_cons = np.zeros((n_v, n_x))
    for v, vn in enumerate(vac_names):
        site = site_of[vn]
        J_cons[v, n_s + v] = th[n_s + v]
        for k, nm in enumerate(surf):
            if mech.species[nm].site == site:
                J_cons[v, k] = th[k]
    blocks = [J_ss / row_scale[:, None], J_cons]
    if n_g:
        J_gas = nu_sub[:, n_s + n_v:].T @ dnet                  # (n_g, n_x)
        for g in range(n_g):                                    # −k_flow·P_i（只打对角）
            J_gas[g, n_s + n_v + g] -= kfl * p_now[g]
        blocks.append(J_gas / gs[:, None])
    J = np.vstack(blocks)
    return F, J, dict(R_f=R_f, z=z, omz=omz, ln_z=ln_z, r_net=r_net, theta=th)


def _hp_inputs(cr):
    """从 CompiledRates 取出高精度后端要的 (Ga, dG, A_pre)。

    `given`/`collision` 的步没有 Ga ⇒ 折算成**等效 Ga**（与 CATKINAS 的 Ean 同一思路）：
    `kf = A·exp(−Ga/kT)` 反解 `Ga = −kT·ln(kf/A)`，这样高精度后端只需处理一种形式。
    """
    import numpy as _np
    from .units import K_B, H_PLANCK, K_B_EV
    kT = K_B_EV * 298.15                 # 只用于形状；真正的 kT 在 hp 内按 T_K 重算
    A = None
    if _np.any(_np.isfinite(cr.A_pre)):
        a = cr.A_pre[_np.isfinite(cr.A_pre)]
        if _np.ptp(a) == 0:
            A = float(a[0])
    return cr.Ga_in.copy(), cr.dG_in.copy(), A


def _flux_project(mech, surf, r0, R_f, eps=1e-16):
    r"""★ 把净速率投影到**严格满足通量守恒**的子空间，按各步可信度加权。

    为什么必须有这一步（2026-10-22 的关键认识）
    --------------------------------------------
    准平衡步的净速率**本来就不该由它自己的 (1−z) 决定**。
    物理上：一个准平衡步的作用是把覆盖度**钉在平衡关系上**（这在 log 空间是精确的），
    而流过它的净通量是由**慢步**决定的 —— 这正是"准平衡"的定义。

    数值上：1−z 的相对精度只有 eps/(1−z)，所以 R_f·(1−z) 在准平衡步上几乎是噪声；
    而同一个通量从守恒关系 νᵀr = 0 推过来却是**精确**的。
    实测合成氨算例：H2 吸附步 1−z = 2.3e-09，直接算的 r₂ 让 H 元素通量差 1.55e-07；
    投影之后降到机器精度。

    做法：在 {r : νᵀ_surf·r = 0} 里找加权最小二乘最近点，权重 w_s = 1/σ_s²。

    ★ 误差棒必须取 **σ_s = R_f,s · eps**（正向速率 × 机器精度），不是 |r⁰_s|·eps/(1−z_s)
    -------------------------------------------------------------------------------
    两者在 1−z > 0 时代数相同（R_f·eps = |r⁰|·eps/(1−z)），但在 **1−z 下溢为 0** 时
    天差地别：那时 r⁰ = 0，`|r⁰|·eps/(1−z)` 兜底成极小值 ⇒ **给了一个毫无信息的步
    最大权重**，把整条解拽向零。

    T7 实测（点 344，E_N=+1.44）：步 1 与步 5 的 1−z 恰好下溢为 0 ⇒
    本包算出的全部 5 个速率**整体偏小 1/14**（0.071429，而通量比例 1.5:0.5:1:1:1
    依然完全精确）—— 典型的"方向对、标度错"。

    换成 σ = R_f·eps 后：净速率 r = R_f·(1−z) 的绝对不确定度就是 R_f·δ(1−z) ≈ R_f·eps
    （因为 (1−z) 由 ln z 得来，其**绝对**精度是 eps）。1−z 越小，R_f 相对 r 就越大，
    权重自动降下去 —— 物理与数值在这里是同一件事。
    """
    idx = [mech.species_order.index(n) for n in surf]
    A = np.asarray(mech.nu, float)[:, idx].T          # (n_surf, n_step)
    r0 = np.asarray(r0, float)
    # ★ 2026-10-22（T7 抓到）：原来写 `np.allclose(r0, 0.0)` —— **错得很隐蔽**。
    #   np.allclose 的默认 atol = 1e-8，于是**任何速率小于 1e-8 的体系都被判成"全零"**，
    #   投影被静默跳过。合成氨 673 K 算例的速率是 1e-8 ~ 1e-31，
    #   401 个点里一大半的 flux_fix 因此恒为 0 —— 看上去"不需要投影"，其实是没执行。
    #   ⇒ 只在**真正全零**时跳过（绝对零，不带容差）。
    if r0.size == 0 or not np.any(np.isfinite(r0)) or np.max(np.abs(r0)) == 0.0:
        return r0.copy(), 0.0
    u, sv, vt = np.linalg.svd(A, full_matrices=True)
    tol = max(A.shape) * np.finfo(float).eps * (sv[0] if sv.size else 1.0)
    rank = int((sv > tol).sum())
    if vt.shape[0] <= rank:                            # 零空间是 {0} ⇒ 无可投影
        return r0.copy(), 0.0
    N = vt[rank:].T                                    # (n_step, d) 通量锥的基
    ref = np.max(np.abs(r0))
    sig = np.abs(np.asarray(R_f, float)) * eps         # ★ σ = R_f·eps（见上）
    sig = np.maximum(sig, ref * eps * 1e-3)            # 下限只为避免除零
    W = 1.0 / sig ** 2
    M_ = N.T @ (W[:, None] * N)
    b = N.T @ (W * r0)
    try:
        c = np.linalg.solve(M_, b)
    except np.linalg.LinAlgError:
        c = np.linalg.lstsq(M_, b, rcond=None)[0]
    r = N @ c
    fix = float(np.max(np.abs(r - r0)) / ref) if ref > 0 else 0.0
    return r, fix


@dataclass
class MultiStateResult:
    """多工况加权的结果（超声热点等）。

    ★ 为什么不能只给一个"平均工况"的解
    ----------------------------------
    速率对温度是**指数敏感**的，所以 `r(⟨T⟩) ≠ ⟨r(T)⟩` —— 差距可以是好几个数量级。
    体系在时间/空间上非均匀时，**必须逐工况求解再加权**：
        r_total = Σ_i w_i·r(state_i) ,  Σ_i w_i = 1
    本容器同时保留各工况的解，好让用户看清"总速率里有多少来自哪个区"。
    """
    rate_net: np.ndarray                  # 加权总净速率
    theta: Dict[str, float]               # 加权平均覆盖度（**仅供参考**，见 note）
    zones: List[Tuple[Dict, float, "SteadyResult"]]   # (工况, 权重, 该工况的解)
    weights: np.ndarray
    converged: bool
    notes: List[str] = field(default_factory=list)

    @property
    def n_zone(self):
        return len(self.zones)

    def contribution(self, step=-1):
        """各工况对某一步总速率的贡献（已乘权重）与占比。"""
        c = np.array([w * r.rate_net[step] for _, w, r in self.zones])
        tot = c.sum()
        return c, (c / tot if tot != 0 else np.full_like(c, np.nan))

    def report(self, step=-1) -> str:
        L = ["多工况加权解（%d 个工况，%s）"
             % (self.n_zone, "全部收敛" if self.converged else "★有工况未收敛")]
        L.append("  权重之和 = %.15f" % float(self.weights.sum()))
        c, frac = self.contribution(step)
        L.append("  工况     T(K)      权重        该工况速率      贡献        占比")
        for i, (st, w, r) in enumerate(self.zones):
            L.append("   %-6d %8.1f  %.6f  %.6e  %.6e  %7.3f%%"
                     % (i + 1, st.get("T_K", float("nan")), w,
                        r.rate_net[step], c[i], 100 * frac[i]))
        L.append("  加权总速率 = %.6e" % self.rate_net[step])
        L.append("  ⚠ 覆盖度是按权重平均的，**没有物理意义上的「总覆盖度」** ——"
                 "各工况的表面状态本来就不同，要看请逐工况看。")
        return "\n".join(L)


def steady_state_multi(mech: Mechanism, en: Energetics, state: Dict,
                       driver: DriverStack, **kw) -> MultiStateResult:
    r"""**多工况加权**求解（T13）：逐工况解稳态，再按权重加和。

    由 `driver.modify_state()` 决定怎么拆（目前只有 sono driver 会拆）。
    `Σw_i = 1` 在这里**强制校验**：权重不归一，总速率就没有意义。
    """
    combos = driver.modify_state(dict(state))
    w = np.array([c[1] for c in combos], float)
    if abs(w.sum() - 1.0) > 1e-9:
        raise SolverError(
            "多工况权重之和 = %.12g ≠ 1 —— 总速率 Σw·r 没有意义。\n"
            "  各工况是**同一个体系的不同部分**，份额必须归一。" % w.sum())
    if np.any(w < 0):
        raise SolverError("出现负权重：%s" % w)

    zones, rate, guess = [], None, None
    all_ok = True
    for st_i, w_i in combos:
        if w_i <= 0.0:
            continue
        r = steady_state(mech, en, state=st_i, driver=driver, **kw)
        if not r.converged:
            # 换个初值再试一次（不同工况的解可能差很多）
            kw2 = dict(kw); kw2.pop("theta_guess", None)
            r2 = steady_state(mech, en, state=st_i, driver=driver, **kw2)
            if r2.converged:
                r = r2
        all_ok &= r.converged
        zones.append((st_i, float(w_i), r))
        rate = (w_i * r.rate_net) if rate is None else rate + w_i * r.rate_net
        guess = dict(r.theta)

    th = {}
    for nm in mech.surface_species:
        th[nm] = float(sum(w_i * r.theta[nm] for _, w_i, r in zones))
    return MultiStateResult(
        rate_net=rate, theta=th, zones=zones,
        weights=np.array([z[1] for z in zones]), converged=all_ok,
        notes=["逐工况求解再加权：r_total = Σ w_i·r(state_i)",
               "★ 不能用 r(⟨T⟩) 代替 ⟨r(T)⟩ —— 速率对温度指数敏感，两者可差几个数量级"])


def steady_state(mech: Mechanism, en: Energetics, state: Dict = None,
                 driver: DriverStack = None, theta_guess: Dict[str, float] = None,
                 p_gas: Dict[str, float] = None, c_liq: Dict[str, float] = None,
                 ladder: Sequence[str] = ("ode", "newton"),
                 t_ode: float = None, ode_stages: int = 4, ode_budget_s: float = 10.0,
                 tol: float = 1e-10, max_iter: int = 200,
                 precision: str = "double", dps: int = 50,
                 flow: "FlowReactor" = None,
                 omz_floor: float = None,
                 strict_energy: bool = True, verbose: bool = False) -> SteadyResult:
    """求稳态覆盖度与各步速率。

    求解是**阶梯**不是单一算法（CATKINAS 的 log 也是这么走的）：

        "ode"    先把瞬态积到足够长，拿一个**物理上可行**的初值
        "newton" 从该初值做解析雅可比的牛顿迭代，精确收敛

    直接从空表面做牛顿常常发散（残差跨 20 个数量级，雅可比在 θ→0 处病态），
    而刚性 ODE 永远不会跑到非物理区 —— 两者**分工**，不是二选一。

    参数
    ----
    state : 工况字典，必须含 `T_K`
    ladder : 走哪几级；给 `("newton",)` 可跳过 ODE（已有好初值时）
    t_ode : ODE 阶段积到多久；不给就按最慢的速率常数自动估
    tol : 缩放残差的 ‖·‖∞ 收敛判据
    flow : 给了就是**流动反应器（CSTR）**：气相分压也成为未知数，
        多出 n_gas 个物料平衡方程（见 `FlowReactor`）。不给 = 微分反应器，
        气相活度由 `p_gas` 固定。
    precision : 求解阶梯的**第三级**（T11 加）
        'double' —— 只用双精度（默认，最快）
        'auto'   —— 先 double，再用 `hp.needs_high_precision()` 判断；
                    判定不足且装了 mpmath ⇒ **自动升到高精度重解**；
                    没装 mpmath ⇒ 退回 double 并在 notes 里**明确说哪些步不可信**
        'high'   —— 直接高精度
    dps : 高精度的十进制位数。实测（T11 ④）：dps=20 完全不够（点 400 上误差 2.9e+04、
          连符号都错），dps=30 起结果平台化；默认 50 留足余量。

    ★ 为什么要有这一级（两轮测量的结论，不是预设）
        T7：E_N ≥ +1.84 的 6 个点，1−Z 压在 4·eps 上，净速率在 double 下错 2.8e-02；
        T8：同一边界让 DRC 的 X_I 退化到 5.2e-05。
        升到高精度后：净速率相对差 2.76e-02 → **2.71e-10**（恰好等于输入常数的精度地板）。
    """
    if precision not in ("double", "auto", "high"):
        raise SolverError("precision 只能是 'double' / 'auto' / 'high'，收到 %r" % precision)
    from scipy.optimize import root

    state = dict(state or {})
    if "T_K" not in state or state["T_K"] is None:
        raise SolverError("state 里没有温度 T_K —— 本包不设温度默认值。")
    T_K = float(state["T_K"])
    surf = mech.surface_species
    vac_names = [mech.vacancy_of(s) for s in mech.sites
                 if mech.vacancy_of(s) in mech.species]
    if not surf:
        raise SolverError("机理里没有表面物种。")
    site_of = {mech.vacancy_of(s): s for s in mech.sites}
    gas_names = list(mech.gas_species) if flow is not None else []

    cr = compile_rates(en, strict_energy=strict_energy)
    rc = apply_state(cr, T_K, driver, state)
    notes = list(cr.notes)

    # 固定（气/液/电子）活度的 ln，预算一次
    cr._idx_of = {n: i for i, n in enumerate(mech.species_order)}
    act_fixed = activities(mech, theta={n: 1.0 for n in surf},
                           p_gas=p_gas, c_liq=c_liq,
                           theta_vac={v: 1.0 for v in vac_names})
    with np.errstate(divide="ignore"):
        cr._ln_a_fixed = np.where(act_fixed > 0, np.log(np.where(act_fixed > 0, act_fixed, 1.0)),
                                  -np.inf)

    def _x_from(g, pg=None):
        xs = ([np.log(max(g.get(n, 1e-8), 1e-300)) for n in surf]
              + [np.log(max(1.0 - sum(g.get(n, 0.0) for n in surf
                                      if mech.species[n].site == site_of[v]), 1e-300))
                 for v in vac_names])
        if gas_names:
            src = pg if pg is not None else (p_gas or {})
            for nm in gas_names:
                v = src.get(nm, src.get(mech.species[nm].formula,
                                        flow.p_in.get(nm, flow.p_in.get(
                                            mech.species[nm].formula, 1e-12))))
                xs.append(np.log(max(float(v), 1e-300)))
        return np.array(xs)

    # ── 第一级：ODE 拿物理初值 —— **自适应时域，不猜** ──
    #
    # ★ 2026-10-22 踩的坑：初版用 `t_ode = 100/min(k)` 估时域。
    #   `min(k)` 是**完全错的尺度**：甲酸算例里慢步 k₂ = 5e-24 s⁻¹，于是估出 10¹² s，
    #   刚性 ODE 从 0 积到 10¹² 直接卡住（测试跑了 2 分钟没出）。
    #   而覆盖度的松弛其实由**快过程**定：τ ≈ 1/(k₁p + k₋₁) ≈ 4e-7 s —— 差 18 个数量级。
    #   物理上：k₂ 极小只说明"产物几乎不生成"，不说明"覆盖度要很久才稳"。
    #
    # ⇒ 改成：从最快时间尺度起步，每级 ×1000 并**用稳态残差本身当判据**
    #   （接着上一级的末态继续积 = 续算，比每级从头来便宜得多）。
    guess = dict(theta_guess or {})
    used = []
    row_scale = np.ones(len(surf))          # 先用 1，拿到初值后再定
    if "ode" in ladder and not theta_guess:
        k_all = np.concatenate([rc.kf[np.isfinite(rc.kf) & (rc.kf > 0)],
                                rc.kr[np.isfinite(rc.kr) & (rc.kr > 0)]])
        k_max = float(np.max(k_all)) if k_all.size else 1.0
        t_stage = (t_ode if t_ode else 10.0 / k_max)
        cur = dict(guess)
        # ★ 2026-10-22（T7）：阶梯必须**有上限**。
        #   ODE 这一级的职责只是给牛顿一个**物理上可行**的初值，不是自己收敛。
        #   T7 的 401 点扫描里，E_N = −1.95 处 N 强吸附、表面中毒（θ_vac ~ 1e-30），
        #   相对残差永远降不到 1e-3，于是 8 级全走完、每级 ×1000 ⇒ 积到天文时间，卡住。
        #   ⇒ 级数上限 + 墙钟预算；超了就把当前最好的点交给牛顿（牛顿能从粗初值走）。
        import time as _time
        _t_begin = _time.time()
        for stage in range(max(1, int(ode_stages))):
            try:
                tr = transient(mech, en, cur, (0.0, t_stage),
                               protocol=Constant(**state), driver=driver,
                               p_gas=p_gas, c_liq=c_liq, n_eval=3,
                               method="BDF", rtol=1e-9, atol=1e-16,
                               strict_energy=False)
            except Exception as _exc:
                used.append("ode(第%d级积分失败:%s ⇒ 用上一级初值)"
                            % (stage + 1, str(_exc).splitlines()[0][:36]))
                break
            cur = {n: float(min(max(tr.theta[n][-1], 1e-300), 1.0)) for n in surf}
            Fp, _, aux_p = _steady_residual(cr, rc, mech, surf, vac_names, site_of,
                                            _x_from(cur), np.ones(len(surf)),
                                            flow, gas_names, T_K)
            sc = np.abs(np.asarray(mech.nu, float)[
                :, [cr._idx_of[n] for n in surf]]).T @ aux_p["R_f"]
            rel = float(np.max(np.abs(Fp[:len(surf)]) / np.where(sc > 0, sc, 1.0)))
            if verbose:
                print("  ODE 第%d级 t=%.3g s ⇒ 相对残差 %.3e  θ=%s"
                      % (stage + 1, t_stage, rel,
                         {k: "%.4g" % v for k, v in cur.items()}))
            used.append("ode(%.2g s, 残差 %.1e)" % (t_stage, rel))
            if rel < 1e-3:
                break
            if _time.time() - _t_begin > ode_budget_s:
                used.append("ode(墙钟预算 %.0fs 用尽 ⇒ 交给牛顿)" % ode_budget_s)
                break
            t_stage *= 1000.0
        guess = cur
        if t_ode is None and len(used) > 1:
            used = [used[0].split("(")[0] + _t("(%d 级续算至 %.2g s)", "(%d continuation stages up to %.2g s)") % (len(used), t_stage)]

    # ★ 流动反应器专用的一级阶梯（T9b 实测逼出来的）
    #
    # 直接从"气相 = 进料"起步会掉进**伪稳态**：产物进料分压是 0，而真解可能是 1.93，
    # 在 ln 空间要跨 28，而牛顿步长限幅是 4 —— 走不过去，回溯又让它停在
    # "什么都不反应"那个同样满足方程的点上。实测 315 点里这样坏掉的近一半，
    # 而且**是"收敛"的**（残差确实小），比不收敛更危险。
    #
    # 修法：先解**微分反应器**（气相钉在进料，这是已经验证过的问题），
    # 拿到物理上可行的表面态 θ；再由物料平衡**直接算出**气相的稳态估计
    #     Σν·r + k_flow·(p_in − p) = 0  ⇒  p ≈ p_in + (Σν·r)/k_flow
    # 这一步是代数的、没有迭代，却一步就把气相送到正确的数量级。
    pg_guess = None
    if flow is not None and gas_names:
        p_in_vec = {nm: float(flow.p_in.get(nm, flow.p_in.get(
            mech.species[nm].formula, 0.0))) for nm in gas_names}
        try:
            pre = steady_state(mech, en, state=state, driver=driver,
                               p_gas={**(p_gas or {}), **p_in_vec}, c_liq=c_liq,
                               theta_guess=theta_guess, ladder=ladder,
                               t_ode=t_ode, ode_stages=ode_stages,
                               ode_budget_s=ode_budget_s, tol=tol,
                               max_iter=max_iter, precision="double",
                               flow=None, strict_energy=False)
        except Exception:
            pre = None
        if pre is not None and pre.converged:
            guess = dict(pre.theta)
            kfl0 = flow.k_flow(T_K)
            idx_all = {n: i for i, n in enumerate(mech.species_order)}
            nu_all = np.asarray(mech.nu, float)
            # ★ 一次代数估计在**高转化率**下不够：`p_in + Σνr/k_flow` 可能直接算成负数
            #   （实测高温低压角就是这样坏的）。⇒ 做几步**自洽迭代**：
            #   在新气相下重解表面 → 再更新气相，并对更新量做阻尼，避免来回振荡。
            pg_guess = dict(p_in_vec)
            last = pre
            for _it in range(4):
                newp, moved = {}, 0.0
                for nm in gas_names:
                    prod = float(nu_all[:, idx_all[nm]] @ last.rate_net)
                    est = p_in_vec[nm] + prod / kfl0 if kfl0 > 0 else p_in_vec[nm]
                    # 负压没有物理意义：说明这一步跨太大 ⇒ 退到当前值与 0 的几何中点
                    if est <= 0.0:
                        est = max(pg_guess[nm] * 0.1, 1e-300)
                    newp[nm] = est
                    ref = max(abs(est), abs(pg_guess[nm]), 1e-300)
                    moved = max(moved, abs(est - pg_guess[nm]) / ref)
                pg_guess = newp
                if moved < 1e-3:
                    break
                try:
                    nxt = steady_state(mech, en, state=state, driver=driver,
                                       p_gas={**(p_gas or {}), **pg_guess}, c_liq=c_liq,
                                       theta_guess=dict(last.theta), ladder=("newton",),
                                       tol=tol, max_iter=max_iter, precision="double",
                                       flow=None, strict_energy=False)
                except Exception:
                    break
                if not nxt.converged:
                    break
                last = nxt
                guess = dict(nxt.theta)
            used.append("flow 初值（微分反应器 + %d 步气相自洽）" % (_it + 1))

    x0 = _x_from(guess, pg_guess)

    # 行缩放：用初值处的正向速率定一次，**之后不再变**（否则雅可比就不是解析的了）
    _, _, aux0 = _steady_residual(cr, rc, mech, surf, vac_names, site_of, x0,
                                  np.ones(len(surf)), flow, gas_names, T_K)
    nu_cols = np.array([cr._idx_of[n] for n in surf])
    scale = np.abs(np.asarray(mech.nu, float)[:, nu_cols]).T @ aux0["R_f"]
    row_scale = np.where(scale > 0, scale, 1.0)

    # ── 第二级：解析雅可比牛顿 ──
    counter = {"n": 0}

    def fun(x):
        counter["n"] += 1
        F, _, _ = _steady_residual(cr, rc, mech, surf, vac_names, site_of, x, row_scale,
                                   flow, gas_names, T_K)
        return F

    def jac(x):
        _, J, _ = _steady_residual(cr, rc, mech, surf, vac_names, site_of, x, row_scale,
                                   flow, gas_names, T_K)
        return J

    # ★ 自写阻尼牛顿（2026-10-22，T7 的 401 点回归逼出来的）
    #
    # 为什么不直接用 scipy 的 hybr：**未知数是 ln θ，而 θ ≤ 1 是硬物理约束**。
    # hybr 不知道这件事，一步就能把 ln θ 推到 +700 以上 ⇒ `np.exp(x)` 溢出，
    # 残差与雅可比整片变 nan，然后它在 nan 里瞎走（T7 首跑 401 点时就是这样炸的）。
    #
    # 三件事缺一不可：
    #   ① **硬截上界 x ≤ 0**（θ ≤ 1）与下界 x ≥ −700（θ 再小也没有意义）；
    #   ② **步长限幅** ‖Δx‖∞ ≤ step_max：一次迭代最多让 θ 变 e^4 ≈ 55 倍，
    #      不然在 θ 跨 20 个数量级的体系里牛顿一步就跳出物理区；
    #   ③ **回溯线搜索**：‖F‖ 不下降就把步长减半，保证单调。
    # 用 lstsq 解牛顿方程（J 在准平衡区会接近奇异，lstsq 不会直接崩）。
    step_max = 4.0
    _hi0 = np.zeros_like(x0)
    if gas_names:
        _hi0[len(surf) + len(vac_names):] = 50.0
    x = np.clip(x0, -700.0, _hi0)
    F, J, _ = _steady_residual(cr, rc, mech, surf, vac_names, site_of, x, row_scale,
                               flow, gas_names, T_K)
    nF = float(np.max(np.abs(F)))
    it = 0
    for it in range(1, max_iter + 1):
        if nF <= tol:
            break
        if not np.all(np.isfinite(J)) or not np.all(np.isfinite(F)):
            break
        try:
            dx = np.linalg.lstsq(J, -F, rcond=None)[0]
        except np.linalg.LinAlgError:
            break
        if not np.all(np.isfinite(dx)) or np.max(np.abs(dx)) == 0.0:
            break
        dx *= min(1.0, step_max / np.max(np.abs(dx)))        # ② 限幅
        lam, ok_ls = 1.0, False
        for _ in range(40):                                   # ③ 回溯
            # ① 硬截：θ ≤ 1 ⇒ ln θ ≤ 0；但**气相分压可以 > 1**（本算例到 750 bar），
            #    所以上界只打在表面/空位那几维上。
            hi = np.zeros_like(x)
            if gas_names:
                hi[len(surf) + len(vac_names):] = 50.0
            xn = np.clip(x + lam * dx, -700.0, hi)
            Fn, Jn, _ = _steady_residual(cr, rc, mech, surf, vac_names, site_of,
                                         xn, row_scale, flow, gas_names, T_K)
            if np.all(np.isfinite(Fn)) and float(np.max(np.abs(Fn))) < nF:
                x, F, J, nF, ok_ls = xn, Fn, Jn, float(np.max(np.abs(Fn))), True
                break
            lam *= 0.5
        if not ok_ls:
            break
    used.append(_t("damped-newton(%d 步, 解析雅可比, 限幅+回溯)",
                   "damped-newton(%d steps, analytic Jacobian, clamp+backtracking)") % it)
    # 兜底：自写牛顿停在高残差上时，再让 hybr 从当前点试一把
    if nF > tol * 10:
        try:
            sol = root(fun, x, jac=jac, method="hybr", tol=tol,
                       options=dict(maxfev=max_iter * (len(x0) + 1)))
            xf = np.clip(sol.x, -700.0, 0.0)
            Ff, _, _ = _steady_residual(cr, rc, mech, surf, vac_names, site_of,
                                        xf, row_scale, flow, gas_names, T_K)
            if np.all(np.isfinite(Ff)) and np.max(np.abs(Ff)) < nF:
                x, nF = xf, float(np.max(np.abs(Ff)))
                used.append("hybr 兜底")
        except Exception:
            pass
    F, _, aux = _steady_residual(cr, rc, mech, surf, vac_names, site_of, x, row_scale,
                                 flow, gas_names, T_K)
    conv = bool(np.max(np.abs(F)) <= max(tol, 1e-12) * 10)

    th = aux["theta"]
    theta = {n: float(th[i]) for i, n in enumerate(surf)}
    tvac = {v: float(th[len(surf) + i]) for i, v in enumerate(vac_names)}
    p_out = ({nm: float(th[len(surf) + len(vac_names) + i])
              for i, nm in enumerate(gas_names)} if gas_names else {})
    site_err = {}
    for v in vac_names:
        s = tvac[v] + sum(theta[n] for n in surf if mech.species[n].site == site_of[v])
        site_err[v] = abs(s - 1.0)

    # ★ 通量投影：准平衡步的净速率由守恒关系推定（见 _flux_project 的说明）
    r_raw = aux["r_net"].copy()
    r_proj, fix = _flux_project(mech, surf, r_raw, aux["R_f"])

    obs = {}
    if driver is not None:
        obs = {k: float(v) for k, v in
               driver.observables(mech, {"rate_net": r_proj}, state).items()}

    res = SteadyResult(
        theta=theta, theta_vac=tvac, rate_fwd=aux["R_f"], zrev=aux["z"],
        rate_net=r_proj, rate_net_raw=r_raw, flux_fix=fix, ln_z=aux["ln_z"],
        residual=F[:len(surf)], resid_abs=F[:len(surf)] * row_scale,
        site_error=site_err, n_iter=counter["n"], method=" → ".join(used),
        converged=conv, obs=obs, notes=notes, p_gas_out=p_out,
        residual_full=F.copy())
    res.notes.append("残差行缩放用的是初值处的 Σ|ν|·r_fwd（固定不变，保证雅可比仍是解析的）")

    # ── 第三级：高精度（T11）────────────────────────────────────────────
    if precision in ("auto", "high"):
        from . import hp as _hp
        need, why = (True, ["precision='high'，直接走高精度"]) if precision == "high" \
            else _hp.needs_high_precision(res, **(
                {} if omz_floor is None else dict(omz_floor=omz_floor)))
        # ★ 2026-10-23（T9b 收尾）：高精度后端**现在支持流动反应器**了。
        #   在此之前它的方程组里没有气相未知数，若放它跑会**静默地解一个不含流动的问题**
        #   —— 表面上「升级成功、精度更高」，实则换了题。当时的做法是明确拒绝；
        #   现在 hp 的未知数排列与本函数完全同构，可以放心升级。
        if need and not _hp.available():
            res.notes.append(
                "★ 判定 double 精度不足，但没装 mpmath ⇒ **仍返回 double 结果**。"
                "原因：" + "；".join(why) + "  装 mpmath 可自动升级。")
        elif need:
            Ga_in, dG_in, A_pre = _hp_inputs(cr)
            hr = _hp.steady_state_hp(
                mech, Ga_in, dG_in, T_K, p_gas=p_gas, c_liq=c_liq, A_pre=A_pre,
                prefactor=cr.prefactor, dps=dps, theta_guess=theta,
                reversible=cr.reversible, flow=flow,
                gas_guess=p_out if p_out else None)
            if hr.converged:
                res = SteadyResult(
                    theta=hr.theta, theta_vac=hr.theta_vac, rate_fwd=hr.rate_fwd,
                    zrev=hr.zrev, rate_net=hr.rate_net, rate_net_raw=aux["r_net"],
                    flux_fix=0.0, ln_z=hr.ln_z, residual=hr.residual[:len(surf)],
                    resid_abs=hr.residual[:len(surf)] * row_scale,
                    site_error={v: abs(hr.theta_vac[v]
                                       + sum(hr.theta[n] for n in surf
                                             if mech.species[n].site == site_of[v]) - 1.0)
                                for v in vac_names},
                    n_iter=hr.n_iter, method=res.method + _t(" → high-precision(%d 位)", " -> high-precision(%d digits)") % dps,
                    converged=True, precision="high", hp=hr,
                    p_gas_out=dict(hr.p_gas_out or p_out),
                    obs=res.obs, notes=res.notes + [
                        "★ 已自动升到高精度（%d 位）。触发原因：%s" % (dps, "；".join(why)),
                        "高精度下 1−z 的最小值 = %.3e（double 给不出这个量级）"
                        % (hr.one_minus_z[hr.one_minus_z > 0].min()
                           if (hr.one_minus_z > 0).any() else 0.0)] + hr.notes)
                if driver is not None:
                    res.obs = {k: float(v) for k, v in driver.observables(
                        mech, {"rate_net": hr.rate_net}, state).items()}
            else:
                res.notes.append("★ 高精度也没收敛（%d 位），返回 double 结果" % dps)
    for w in res.precision_warnings():
        res.notes.append("⚠ " + w)
    if not conv:
        res.notes.append("★ 未达收敛判据：‖F‖∞ = %.3e > %.1e。scipy 消息：%s"
                         % (np.max(np.abs(F)), tol * 10, sol.message))
    return res
