# -*- coding: utf-8 -*-
r"""升精度判据的阈值敏感性 —— 回答"阈值换一个，结论会不会变"。

论文 §3.3 写明自动升精度由 `1−z`、通量投影幅度和完整残差三条判据触发（默认
`1−z < 1e-10`、`flux_fix > 1e-6`），但审稿人一定会追问：**阈值挪一个量级呢**？
这个脚本给出答案：在 401 点描述符扫描上把 `omz_floor` 取 1e-8 / 1e-10 / 1e-12，
报告每档

    ① 有多少点被判为需要升精度；
    ② 混合后对 CATKINAS 参考值的**最大相对偏差**。

★ 算法上要省：三档的差别**只在"判不判它升"**，并不改变任何一个解本身。
  所以 double 解算一遍、高精度解算一遍（两次全场扫描，各几秒），
  三档只是从这两套结果里按判据挑选，不必反复求解。
  第一版逐点独立求解（不续算）跑了 13 分钟还没出第一行 —— 教训是
  **扫描要用 `scan()` 的续算，别自己写 for 循环从零起算**。

    python paper/threshold_sensitivity.py
"""
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "tests"))

import refdata                                        # noqa: E402
from kinstep import mech as M                         # noqa: E402
from kinstep.scaling import from_table                # noqa: E402
from kinstep.scan import scan                         # noqa: E402
from kinstep.units import K_B_EV                      # noqa: E402
import kinstep.drivers                                # noqa: E402,F401

MECH = M.from_text(
    "H2 + 2* <-> 2H*\nN2 + 2* <-> 2N*\nN* + H* <-> NH* + *\n"
    "NH* + H* <-> NH2* + *\nNH2* + H* <-> NH3 + 2*\n")
KT, A_PRE, T_NOM = 0.057994719627, 1.402306113158e13, 673.0
T_EQ = T_NOM * KT / (K_B_EV * T_NOM)
P = {"H2": 75.0, "N2": 25.0, "NH3": 1.0}
# 与 T7/T9 同一套 scaling / BEP 系数（CATKINAS 发表输入里的四个数）
TAB = {1: dict(bep=(0.54, 0.60), scale=(0.30, -0.68)),
       2: dict(bep=(0.78, 1.54), scale=(2.00,  0.00)),
       3: dict(bep=(0.26, 1.45), scale=(-0.15, -0.18)),
       4: dict(bep=(0.43, 1.66), scale=(-0.45, -0.45)),
       5: dict(bep=(0.31, 1.45), scale=(-0.85,  1.17))}

FLOORS = (1e-8, 1e-10, 1e-12)
FLUX_FLOOR = 1e-6          # 通量投影判据，三档都用默认值（本研究只动 omz_floor）
F_CATKINAS = 96485.2608    # 对齐发行物的 J/mol→eV 换算常数（见 T10 §②）


def main():
    d, src = refdata.catkinas_vars("curve")
    Ean, G0n, RNET = (np.asarray(d[k], float) for k in ("Ean", "G0n", "Rnet"))
    print("参考值来源：%s；%d 个描述符点" % (src, Ean.shape[0]))

    # 熵校正用 T9/T10 同一条路：直接取发行物逐点能量与标度关系的差
    from kinstep.thermo import entropy_offsets
    # ★ 要用 CATKINAS 的 F 常数，否则比的是常数不是公式（见正文 §4）
    off_Ga, off_dG, _ = entropy_offsets(MECH, T_NOM,
                                        J_per_mol_per_eV=F_CATKINAS)
    sc = from_table(MECH, TAB, desc_names=["E_N"],
                    offset_Ga=off_Ga, offset_dG=off_dG)
    E1 = np.arange(-2.0, 2.0001, 0.01)

    # ★ 走**与正文 §4 完全相同的那条路**：逐点用发行物自己的能量 (Ean, G0n)，
    #   由中心 (200) 向两侧续算，precision="auto"。
    #   走 scaling→热化学重建那条路（T10 §⑤）会得到 234 点 / 2.99e-9 —— 也对，
    #   但那是另一个口径，与正文引用的数字不是一回事，混用会让读者对不上。
    from kinstep.analyze import energetics_from_arrays       # noqa: E402
    from kinstep.solve import steady_state                   # noqa: E402

    ref = RNET[:, 4]
    m = np.isfinite(ref) & (np.abs(ref) > 0)
    MID = 200
    order = [MID] + [MID + k for k in range(1, 201)] + [MID - k for k in range(1, 201)]

    def run(floor):
        rn = np.zeros(401)
        n_hi = 0
        guess = None
        t0 = time.time()
        for i in order:
            en = energetics_from_arrays(MECH, Ean[i], G0n[i], A_pre=A_PRE)
            kw = dict(state=dict(T_K=T_EQ), p_gas=P, tol=1e-12, dps=60,
                      strict_energy=False)
            if floor is None:                       # 对照：全程 double
                r = steady_state(MECH, en, theta_guess=guess,
                                 ladder=("newton",) if guess else ("ode", "newton"),
                                 precision="double", **kw)
            else:
                r = steady_state(MECH, en, theta_guess=guess,
                                 ladder=("newton",) if guess else ("ode", "newton"),
                                 precision="auto", omz_floor=floor, **kw)
            if not r.converged:
                r = steady_state(MECH, en,
                                 precision="double" if floor is None else "auto",
                                 omz_floor=floor, **kw)
            rn[i] = r.rate_net[4]
            n_hi += (r.precision == "high")
            guess = dict(r.theta)
        dev = float(np.nanmax(np.abs(rn[m] / ref[m] - 1.0)))
        return n_hi, dev, time.time() - t0

    rows = []
    n_hi, dev, dt = run(None)
    print("  （对照）全程 double        ：升精度 %3d / 401，最大相对偏差 %.3e，%.1f s"
          % (n_hi, dev, dt))
    for floor in FLOORS:
        n_hi, dev, dt = run(floor)
        rows.append((floor, n_hi, dev, dt))
        print("  omz_floor = %-8.0e      ：升精度 %3d / 401，最大相对偏差 %.3e，%.1f s"
              % (floor, n_hi, dev, dt))

    print("\n  %-12s %12s %18s %10s" % ("omz_floor", "升精度点/401", "最大相对偏差", "耗时"))
    for floor, nh, dev, dt in rows:
        print("  %-12.0e %12d %18.3e %9.1f s" % (floor, nh, dev, dt))

    print("\n判读：三档的最大偏差若落在同一量级，说明结论不依赖阈值的具体取值；")
    print("      默认 1e-10 取在『偏差已到参考值的存储精度』与『别多算』之间。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
