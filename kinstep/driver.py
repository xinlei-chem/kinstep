# -*- coding: utf-8 -*-
# KinSTEP - a modality-agnostic microkinetic solver.
# Copyright (C) 2026 The KinSTEP authors.
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.  It is distributed WITHOUT ANY WARRANTY; see the GNU General Public
# License (LICENSE) for details.
r"""[English] The driver plug-in layer, the extension point of the package. A driver returns
perturbed energies (and optionally weighted sub-states); that is the only thing
the solver core knows about an external stimulus.

(The commentary below is in Chinese: it records why each choice was made
and which measurement forced it.)

**驱动（Driver）插件层** —— 本包的扩展骨架。

设计立意
--------
微动力学的**内核**（机理解析 → 速率律 → 稳态/瞬态求解 → DRC）对所有催化模态
**完全相同**；不同的只有两件事：

  ① 外场怎么**改变基元步的能量**（Ga_fwd / dG_rxn）；
  ② 该模态**独有的可观测量**（电流密度、量子效率、声化学产额 …）。

所以把这两件事抽成 `Driver` 接口，内核永远不必知道"这是电催化还是光催化"。

    ┌─────────────────────────── core（永不随模态改动）──────────────────────────┐
    │  mech.py 机理→化学计量矩阵    rates.py 速率律    solve.py 求解    analyze.py │
    └───────────────────────────────────┬──────────────────────────────────────┘
                                        │ 只通过这一个口子
    ┌───────────────────────────────────▼──────────────────────────────────────┐
    │  Driver:  modify_energies() / extra_species_steps() / observables()      │
    │  thermal │ electro │ photo │ sono │ field │ plasma │ …（新增不动内核）      │
    └──────────────────────────────────────────────────────────────────────────┘

★ 可叠加（这是必须的，不是锦上添花）
    光电催化 PEC = photo + electro 同时作用；
    场效应电催化 = electro + field。
  因此用 `DriverStack` 串联，按注册顺序依次修正能量。

★ 参数撞名在**注册时**就查（承 naming.py 的规矩）
  每个 driver 必须声明自己要用的"工况参数名"与"逐步参数名"，
  两个 driver 若声明了同名参数 ⇒ 注册即报错，不等到运行时才发现。
"""
from __future__ import annotations

import abc
import warnings
from typing import Any, Dict, Iterable, List, Sequence, Tuple

from .naming import RESERVED_NAMES, NameCollisionError

#: 全局 driver 注册表：名字 → 类
_REGISTRY: Dict[str, type] = {}
#: 参数名 → 声明它的 driver 名（用于撞名检查）
_PARAM_OWNER: Dict[str, str] = {}
#: 可观测量名 → 声明它的 driver 名
_OBS_OWNER: Dict[str, str] = {}

#: 本包支持的 driver API 版本。老插件声明低版本时，内核按低版本能力降级运行。
DRIVER_API_VERSION = 1


class DriverParamError(ValueError):
    """driver 需要的参数缺失或类型不对。"""


def register_driver(cls):
    """类装饰器：把 driver 登记进注册表，并在**注册时**做撞名检查。

    为什么在注册时查而不是运行时
    ----------------------------
    运行时才发现"两个 driver 都用了 `intensity`"，往往已经算了一半，
    而且症状是"结果不对"而不是"报错"。注册时查 = 导入包就暴露。
    """
    name = getattr(cls, "driver_name", None)
    if not name:
        raise ValueError("%s 必须定义类属性 driver_name" % cls.__name__)
    if name in _REGISTRY:
        raise NameCollisionError("driver 名 '%s' 已被 %s 占用" % (name, _REGISTRY[name].__name__))

    declared = set(cls.state_params) | set(cls.step_params)
    # ① 不得撞包内保留名
    hit_reserved = {p for p in declared if p.lower() in RESERVED_NAMES}
    # —— 例外：driver 本来就是这些量的**主人**（如 electro 拥有 U_elec）——
    hit_reserved -= set(getattr(cls, "owns_reserved", ()))
    if hit_reserved:
        raise NameCollisionError(
            "driver '%s' 的参数撞保留名：%s\n"
            "（若该 driver 确实是这些量的主人，请在类里写 owns_reserved = (...)）"
            % (name, sorted(hit_reserved)))
    # ② 不得与别的 driver 撞名
    for p in declared:
        owner = _PARAM_OWNER.get(p)
        if owner and owner != name:
            raise NameCollisionError(
                "参数名 '%s' 冲突：driver '%s' 与 '%s' 都声明了它。\n"
                "请给其中一个加前缀，例如 '%s_%s'。" % (p, name, owner, name, p))
    # ③ 可观测量名查重（T2.5 加：原先只在运行时查，等于算完才报错）
    for ob in getattr(cls, "provides", ()):
        owner = _OBS_OWNER.get(ob)
        if owner and owner != name:
            raise NameCollisionError(
                "可观测量 '%s' 冲突：driver '%s' 与 '%s' 都声明要产出它。" % (ob, name, owner))
    for ob in getattr(cls, "provides", ()):
        _OBS_OWNER[ob] = name
    for p in declared:
        _PARAM_OWNER[p] = name
    _REGISTRY[name] = cls
    return cls


def available_drivers() -> List[str]:
    """已注册的 driver 名字（导入 `kinstep.drivers` 即全部注册）。"""
    return sorted(_REGISTRY)


def make_driver(name: str, **kwargs):
    """按名字造一个 driver；`kwargs` 作为**默认工况**存进 `cfg`
    （运行时 `state` 里的同名项优先，见 `Driver.get_param`）。"""
    if name not in _REGISTRY:
        raise KeyError("未注册的 driver '%s'；已注册：%s" % (name, available_drivers()))
    return _REGISTRY[name](**kwargs)


def param_owner_table() -> Dict[str, str]:
    """参数名 → 所属 driver。给用户排查"这个参数是谁要的"。"""
    return dict(_PARAM_OWNER)


class Driver(abc.ABC):
    """所有催化模态驱动的基类。

    子类必须提供
    ------------
    driver_name  : str              唯一名字
    state_params : tuple[str]       **工况级**参数（整个体系一个值，如 U_elec、photon_flux）
    step_params  : tuple[str]       **逐基元步**参数（每步一个值，如 nu_elec、beta_ct）
    owns_reserved: tuple[str]       （可选）本 driver 合法拥有的保留名

    子类可覆盖
    ----------
    modify_energies()      ★ 主接口：改 Ga_fwd / dG_rxn
    extra_species_steps()  额外的物种与基元步（如光生载流子的产生/复合）
    observables()          本模态独有的可观测量（j_curr / AQY / …）
    validate()             投产前的自检
    """

    # ── 身份 ──
    driver_name: str = ""
    state_params: Tuple[str, ...] = ()
    step_params: Tuple[str, ...] = ()
    owns_reserved: Tuple[str, ...] = ()

    # ── 能力声明（T2.5 加）：内核据此决定怎么调用，老插件也能在新内核上跑 ──
    #: 本 driver 写作时针对的 API 版本
    api_version: int = 1
    #: True ⇒ 能量依赖覆盖度（侧向相互作用等）⇒ 内核必须在**牛顿内层**每次重算能量；
    #: False ⇒ 能量与 theta 无关 ⇒ 求解前算一次即可（快得多）
    needs_theta: bool = False
    #: True ⇒ 本 driver 的修正是**可加的**，与其它 additive driver 换序结果不变；
    #: False ⇒ 含非线性项（如 field 的二次项），顺序影响结果，Stack 会警告
    is_additive: bool = True
    #: 本 driver 产出的可观测量名（注册时查重）
    provides: Tuple[str, ...] = ()
    #: True ⇒ 本 driver 会改"工况"（如超声热点）⇒ 内核需按多工况加权求解
    splits_state: bool = False

    def __init__(self, **kwargs):
        unknown = set(kwargs) - set(self.state_params)
        if unknown:
            raise DriverParamError(
                "driver '%s' 不认识参数 %s；它只接受 %s"
                % (self.driver_name, sorted(unknown), list(self.state_params)))
        #: 构造时给的值 = **默认值**；运行时 state 里的同名项优先级更高
        self.cfg = dict(kwargs)

    # ── 参数取值的唯一入口（T2.5 定清楚）────────────────────────────────
    def get_param(self, state, key, default=None, required=False):
        """按 **state ▸ cfg ▸ default** 的优先级取工况参数。

        为什么必须定这个优先级
        ----------------------
        扫描时 `U_elec` 是**逐点变化**的，只能来自 `state`；而构造时给的值
        应当是"没特别说明时用这个"的默认。初版两处来源没分清楚，
        `make_driver(U_elec=0.3)` 放进 cfg，driver 却去 `state["U_elec"]` 读 —— 
        单点能跑、一扫描就 KeyError（本测试当场抓到）。
        """
        if state is not None and key in state and state[key] is not None:
            return state[key]
        if key in self.cfg and self.cfg[key] is not None:
            return self.cfg[key]
        if required:
            raise DriverParamError(
                "driver '%s' 需要工况参数 '%s'，但 state 和构造参数里都没有。\n"
                "给法：make_driver('%s', %s=值)  或  run(..., %s=值)"
                % (self.driver_name, key, self.driver_name, key, key))
        return default

    # ── 主接口 ────────────────────────────────────────────────────────────
    @abc.abstractmethod
    def modify_energies(self, mech, Ga_fwd, dG_rxn, state, theta=None):
        """返回修正后的 (Ga_fwd, dG_rxn)。

        参数
        ----
        Ga_fwd, dG_rxn : ndarray —— **纯 numpy 数组，不是 Val**
            来源标签在"输入/校验层"已被 `compile()` 剥掉。
            实测 Val 进内层循环比 ndarray 慢 5 倍（20 步×20 万次：0.923 s vs 0.172 s），
            而真实内层还有 exp/连乘/雅可比，差距更大。
        theta : ndarray | None
            覆盖度。**只有 `needs_theta=True` 的 driver 才会拿到非 None**；
            此时内核会在牛顿迭代内层反复调用本方法（自洽）。
            `needs_theta=False` 时内核只在求解前调一次，theta 恒为 None。

        约定
        ----
        **只返回新数组，不就地改** —— 就地改会让多 driver 叠加时出现
        "谁先谁后结果不同、还查不出来"。
        """
        raise NotImplementedError

    # ── 可选接口（默认什么都不做）────────────────────────────────────────
    def modify_state(self, state):
        """把一个工况拆成**多个加权工况** —— 返回 [(state_i, weight_i), ...]。

        为什么需要这个口子（T2.5 加）
        ----------------------------
        `modify_energies` 假定"外场只改能量"，但这对超声**根本不成立**：
        空化热点制造的是微秒级的局部 T/P 脉冲，体系实际处在**时间/空间上非均匀**
        的工况里 ⇒ 总速率 = Σ wᵢ·r(stateᵢ)，而不是 r(单一修正后的 state)。

        这个口子同时服务：
          · 超声热点两态近似
          · 非均匀表面（台阶/平台位点按比例混合）
          · 温度/压力分布（而非单一值）

        默认返回 [(state, 1.0)]，即"单一工况"，对绝大多数 driver 正确。
        """
        return [(state, 1.0)]

    def extra_species_steps(self, mech, state):
        """返回 (额外物种列表, 额外基元步列表)。

        例：光催化里的 `e-(cb)` / `h+(vb)` 与"光生载流子产生""体相复合"两步，
        它们不是化学基元步，但在动力学上必须进方程组。
        """
        return [], []

    def observables(self, mech, solution, state) -> Dict[str, float]:
        """本模态独有的可观测量。内核把各 driver 的结果 merge 进总输出。"""
        return {}

    def validate(self, mech) -> None:
        """投产前自检：逐步参数是否齐全、取值是否合理。"""
        for sid, step in enumerate(mech.steps):
            for p in self.step_params:
                if p not in step.params:
                    raise DriverParamError(
                        "driver '%s' 需要第 %d 步提供参数 '%s'，但没找到。\n"
                        "在机理里这样写： {%s: 值}" % (self.driver_name, sid + 1, p, p))

    def capabilities(self) -> Dict[str, Any]:
        """能力摘要 —— 内核与 CLI 都读它，也写进输出以便复现。"""
        return dict(api_version=self.api_version, needs_theta=self.needs_theta,
                    is_additive=self.is_additive, splits_state=self.splits_state,
                    provides=tuple(self.provides))

    def describe(self) -> str:
        flags = []
        if self.needs_theta:
            flags.append("需θ自洽")
        if not self.is_additive:
            flags.append("非可加")
        if self.splits_state:
            flags.append("多工况")
        tag = ("[" + "/".join(flags) + "]") if flags else ""
        return "%s%s(工况=%s, 逐步=%s)" % (
            self.driver_name, tag, list(self.state_params), list(self.step_params))


class DriverStack:
    """把多个 driver 串起来 —— 光电催化 = photo + electro 就靠它。

    能量修正**按注册顺序依次施加**；每个 driver 拿到的是上一个的输出。
    ⚠ 顺序会影响结果（非线性修正时），所以 `describe()` 会把顺序打出来，
      并写进输出，**让顺序可追溯**。
    """

    def __init__(self, drivers: Sequence[Driver], warn_order: bool = True):
        self.drivers = list(drivers)
        names = [d.driver_name for d in self.drivers]
        if len(set(names)) != len(names):
            raise NameCollisionError("同一个 driver 被重复加入：%s" % names)
        # 顺序敏感性检查（T2.5 加）：只要有一个非可加 driver 且总数 >1，顺序就可能影响结果
        self.order_sensitive = (len(self.drivers) > 1
                                and any(not d.is_additive for d in self.drivers))
        if self.order_sensitive and warn_order:
            bad = [d.driver_name for d in self.drivers if not d.is_additive]
            warnings.warn(
                "driver 叠加顺序会影响结果：%s 含非线性修正（is_additive=False）。\n"
                "当前顺序：%s —— 该顺序已写入输出以便复现。"
                % (bad, names), RuntimeWarning, stacklevel=2)
        # 老插件兼容：按最低 api_version 降级
        self.api_version = min([d.api_version for d in self.drivers] + [DRIVER_API_VERSION])

    #: 本 stack 里是否有 driver 要 theta ⇒ 内核据此决定能量在内层还是外层算
    @property
    def needs_theta(self) -> bool:
        return any(d.needs_theta for d in self.drivers)

    @property
    def splits_state(self) -> bool:
        return any(d.splits_state for d in self.drivers)

    def modify_energies(self, mech, Ga_fwd, dG_rxn, state, theta=None):
        for d in self.drivers:
            # 不需要 theta 的 driver 一律拿不到它 —— 免得有人"顺手"用了却没声明，
            # 导致内核以为可以在外层算一次，结果其实是自洽的
            th = theta if d.needs_theta else None
            Ga_fwd, dG_rxn = d.modify_energies(mech, Ga_fwd, dG_rxn, state, th)
        return Ga_fwd, dG_rxn

    def modify_state(self, state):
        """多个 driver 各自拆工况 ⇒ **笛卡尔积**，权重相乘。"""
        combos = [(state, 1.0)]
        for d in self.drivers:
            if not d.splits_state:
                continue
            new = []
            for st, w in combos:
                for st2, w2 in d.modify_state(st):
                    new.append((st2, w * w2))
            combos = new
        return combos

    def extra_species_steps(self, mech, state):
        sp_all, st_all = [], []
        for d in self.drivers:
            sp, st = d.extra_species_steps(mech, state)
            sp_all += list(sp)
            st_all += list(st)
        return sp_all, st_all

    def observables(self, mech, solution, state) -> Dict[str, float]:
        out: Dict[str, float] = {}
        for d in self.drivers:
            for k, v in d.observables(mech, solution, state).items():
                if k in out:
                    raise NameCollisionError(
                        "可观测量 '%s' 被多个 driver 同时输出（%s）" % (k, self.describe()))
                out[k] = v
        return out

    def validate(self, mech) -> None:
        for d in self.drivers:
            d.validate(mech)

    def required_state_params(self) -> List[str]:
        return sorted({p for d in self.drivers for p in d.state_params})

    def describe(self) -> str:
        return " ▸ ".join(d.describe() for d in self.drivers) or "(无 driver：纯热催化)"
