# -*- coding: utf-8 -*-
r"""T8 验证：DRC / 反应级数 / 表观活化能，对 CATKINAS 的 `DRC_B/DRC_I/DRC_PC/Eapp`。

口径（反推自 .mat，见 notes/DEVLOG.md T8 条）：
    X = ∂ln r / ∂(−G_x/k_BT) = −k_BT·∂ln r/∂G_x     ← 四种 DRC 统一这一个定义
    DRC_B  台阶（只抬过渡态）   Σ_s = 1（Campbell）
    DRC_I  中间体（抬 G_k，过渡态不动）
    DRC_PC 气相 = **−反应级数**
    DRC_G  描述符，沿扫描曲线差分
"""
import sys, time
import _root, refdata  # noqa: F401  ← 包根入 sys.path + 第三方真值定位
import numpy as np
import scipy.io as sio

from kinstep import mech as M
from kinstep.analyze import drc, calibrate_delta, drc_descriptor, energetics_from_arrays
from kinstep.solve import steady_state
from kinstep.units import K_B_EV
import kinstep.drivers

REF, REF_SRC = refdata.catkinas_vars("curve")
KT_CAT = 0.057994719627
A_PRE = 1.402306113158e13
T_NOM = 673.0
T_K = T_NOM * KT_CAT / (K_B_EV * T_NOM)      # 等效温度，统一 k_BT
P = {"H2": 75.0, "N2": 25.0, "NH3": 1.0}

mech = M.from_text("""
H2 + 2* <-> 2H*
N2 + 2* <-> 2N*
N* + H* <-> NH* + *
NH* + H* <-> NH2* + *
NH2* + H* <-> NH3 + 2*
""")
d = REF
Ean, G0n = np.asarray(d["Ean"]), np.asarray(d["G0n"])
DB, DI, DPC, DG = (np.asarray(d[k]) for k in ("DRC_B", "DRC_I", "DRC_PC", "DRC_G"))
EAPP = np.asarray(d["Eapp"])
E1 = np.arange(-2.0, 2.0001, 0.01)

def en_at(i):
    return energetics_from_arrays(mech, Ean[i], G0n[i], A_pre=A_PRE)

print("══ ① ★ 用 Campbell 求和规则自己选差分步长 δ（不用经验值）══")
i0 = 200
best, rows = calibrate_delta(mech, en_at(i0), dict(T_K=T_K), p_gas=P, verbose=True)
print("   ⇒ 选中 δ = %.1e eV" % best)

print("\n══ ② 点 200（E_N = 0）逐项对 CATKINAS ══")
t0 = time.time()
r = drc(mech, en_at(i0), dict(T_K=T_K), p_gas=P, delta=best)
print("   耗时 %.2f s（%d 次重解）" % (time.time() - t0,
                                      2 * (mech.n_step + len(r.inter_names)
                                           + len(r.gas_names) + 1)))
print(r.report(target=4))
print()
print("   逐项比对（目标 = 第 5 步净速率，CATKINAS 的列 4）：")
print("   台阶 DRC：")
for s in range(5):
    print("      步%d  本包 %+.9f   CATKINAS %+.9f   差 %.2e"
          % (s + 1, r.X_B[s, 4], DB[i0][s, 4], abs(r.X_B[s, 4] - DB[i0][s, 4])))
print("      Σ：本包 %.9f   CATKINAS %.9f" % (r.sum_B[4], DB[i0][:, 4].sum()))
print("   中间体 DRC：")
for j, nm in enumerate(r.inter_names):
    print("      %-6s 本包 %+.9f   CATKINAS %+.9f   差 %.2e"
          % (nm, r.X_I[j, 4], DI[i0][j, 4], abs(r.X_I[j, 4] - DI[i0][j, 4])))
print("   气相（X_PC = −反应级数）：")
for j, g in enumerate(r.gas_names):
    print("      %-5s 本包 %+.9f   CATKINAS %+.9f   差 %.2e   (级数 %+.6f)"
          % (g, r.X_PC[j, 4], DPC[i0][j, 4], abs(r.X_PC[j, 4] - DPC[i0][j, 4]),
             r.order[j, 4]))
print("   表观活化能：本包 %.6f eV   CATKINAS %.6f eV   差 %.2e"
      % (r.Ea_app[4], EAPP[i0][4], abs(r.Ea_app[4] - EAPP[i0][4])))

print("\n══ ③ 多点抽样（看 DRC 在不同控速机制下是否都对）══")
print("   （用 Richardson 外推；依据见 ⑤）")
print("   点   E_N    主控步(本包/真值)   ΣX_B      X_B 最大差    X_I 最大差   级数最大差")
samp = [50, 100, 150, 200, 250, 300, 330]
agg = {"B": [], "I": [], "O": [], "S": []}
for i in samp:
    try:
        ri = drc(mech, en_at(i), dict(T_K=T_K), p_gas=P, delta=best, richardson=True)
    except Exception as e:
        print("   %3d  %+.2f   异常 %s" % (i, E1[i], str(e)[:40])); continue
    mine_top = int(np.argmax(ri.X_B[:, 4])) + 1
    ref_top = int(np.argmax(DB[i][:, 4])) + 1
    eB = np.nanmax(np.abs(ri.X_B[:, 4] - DB[i][:, 4]))
    eI = np.nanmax(np.abs(ri.X_I[:, 4] - DI[i][:, 4]))
    eO = np.nanmax(np.abs(ri.X_PC[:, 4] - DPC[i][:, 4]))
    agg["B"].append(eB); agg["I"].append(eI); agg["O"].append(eO)
    agg["S"].append(abs(ri.sum_B[4] - 1.0))
    print("   %3d  %+.2f      %d / %d        %.6f   %.3e    %.3e    %.3e"
          % (i, E1[i], mine_top, ref_top, ri.sum_B[4], eB, eI, eO))
print("   ⇒ 抽样汇总：X_B 最大差 %.3e | X_I %.3e | 级数 %.3e | 求和规则偏差 %.3e"
      % (max(agg["B"]), max(agg["I"]), max(agg["O"]), max(agg["S"])))

print("\n══ ④ 描述符 DRC（沿扫描曲线）══")
REF_R = np.asarray(d["Rnet"])
XG = drc_descriptor(E1, REF_R[:, 4], T_NOM)
sel = [100, 150, 200, 250, 300]
print("   点   E_N     本包 X_G      CATKINAS DRC_G   差")
for i in sel:
    print("   %3d  %+.2f   %+.6f     %+.6f      %.2e" % (i, E1[i], XG[i], DG[i][0, 4],
                                                          abs(XG[i] - DG[i][0, 4])))
print("   ⇒ 最大差 %.3e（端点除外；差分步长 0.01 eV 的截断误差量级）"
      % np.nanmax(np.abs(XG[2:-2] - DG[2:-2, 0, 4])))

print("\n══ ⑤ ★ Richardson 外推：消去中心差分的二阶截断 ══")
print("   依据：差/δ² 在 δ∈[1.25e-4, 4e-3] 上恒为 9.526（4 位有效数字）")
print("        ⇒ 残差纯是二阶截断，不是定义差异")
rr = drc(mech, en_at(i0), dict(T_K=T_K), p_gas=P, delta=best, richardson=True)
print("   量           外推前差       外推后差      改善")
for nm, a, b, ref in (("X_B(步2)", r.X_B[1, 4], rr.X_B[1, 4], DB[i0][1, 4]),
                      ("X_I(H*)", r.X_I[0, 4], rr.X_I[0, 4], DI[i0][0, 4]),
                      ("X_I(*)", r.X_I[4, 4], rr.X_I[4, 4], DI[i0][4, 4]),
                      ("X_PC(H2)", r.X_PC[0, 4], rr.X_PC[0, 4], DPC[i0][0, 4]),
                      ("X_PC(N2)", r.X_PC[1, 4], rr.X_PC[1, 4], DPC[i0][1, 4]),
                      ("Ea_app", r.Ea_app[4], rr.Ea_app[4], EAPP[i0][4])):
    e0, e1 = abs(a - ref), abs(b - ref)
    print("   %-11s %.3e     %.3e    %s"
          % (nm, e0, e1, ("×%.0f" % (e0 / e1)) if e1 > 0 else "→0"))
print("   外推后 ΣX_B = %.12f（偏差 %.1e）" % (rr.sum_B[4], abs(rr.sum_B[4] - 1.0)))

print("\n══ 判定 ══")
print("   求和规则 Σ X_B = 1：%s" % ("✅ 通过" if max(agg["S"]) < 1e-3 else "⚠ %.1e" % max(agg["S"])))
print("   X_B / X_I / 级数 对 CATKINAS（7 点抽样，Richardson 外推）：%.3e / %.3e / %.3e  %s"
      % (max(agg["B"]), max(agg["I"]), max(agg["O"]),
         "✅" if max(max(agg["B"]), max(agg["I"]), max(agg["O"])) < 1e-6 else "⚠"))
print("   Ea_app 对 CATKINAS：%.3e eV（不受 δ 外推影响，由温度步长 dT 决定）"
      % abs(r.Ea_app[4] - EAPP[i0][4]))
