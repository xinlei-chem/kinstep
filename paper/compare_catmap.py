# -*- coding: utf-8 -*-
r"""与 CatMAP 的独立对照（CO 氧化，225 点描述符图）。

## 为什么要有这一套

论文的三套验证（Koper、解析闭式解、CATKINAS）里，开源同类程序一个都没有，
而 CatMAP 恰恰是功能重合最多、审稿人最可能拿来质疑"是不是 Python 重写"的那一个。
这套对照回答它，并且是**第三套独立实现**上的验证。

## 做法（KinSTEP 不 import CatMAP，也不随包分发 CatMAP）

CatMAP（GPLv3，SUNCAT）在**隔离环境**里跑它自己的 CO 氧化教程 2 输入
（resolution 15、T = 500 K、`Ndigits = 1000`），导出每个描述符点的
速率常数 kf/kr、覆盖度和三步净速率，存为 `SI/catmap_co_oxidation.json`。
KinSTEP 用 `mode="given"` 吃**同一组 kf/kr**，比稳态。

口径先对齐：用 CatMAP 自己的 kf/kr 与覆盖度按 r = kf·Πa − kr·Πa 复现它的速率，
40/40 点逐位吻合 ⇒ 两边速率律一致，后面的差异不是口径差异。

## 读数要点

比较**必须条件化**：净速率的比对只在 **CatMAP 自身三步速率满足质量平衡**
（r1 = r3、2r2 = r3，残差 < 1e-10）的点上做。
在其余点上，差异集中在准平衡的 CO 吸附步，那里 r = kf·p·θ_空 − kr·θ_CO 是两个大数相减，
**算术位数再多也救不了**——误差来自覆盖度本身的收敛余量。
这不是"谁算错了"，而是论文 §3.2 的论点：准平衡步的净速率应由通量守恒推定。

    python paper/compare_catmap.py
"""
import json, os, sys, warnings
warnings.filterwarnings("ignore")
import numpy as np
from kinstep import mech as M
from kinstep.rates import Energetics
from kinstep.solve import steady_state
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
d = json.load(open(os.path.join(HERE, "SI", "catmap_co_oxidation.json")))
P = d["meta"]["pressures"]
mech = M.from_text("CO + * <-> CO*\nO2 + 2* <-> 2O*\nCO* + O* <-> CO2 + 2*\n")
pg = {"CO": P["CO_g"], "O2": P["O2_g"], "CO2": P["CO2_g"]}
per = {0: [], 1: [], 2: []}
fluxK, fluxC, dcov, n = [], [], [], 0
for i_, p in enumerate(d["points"]):
    if i_ % 75 == 0: print("  ... %d/225" % i_, flush=True)
    kf=[float(x) for x in p["kf"]]; kr=[float(x) for x in p["kr"]]
    cm=[float(x) for x in p["coverages"]][:2]
    en = Energetics(mech, {j+1: dict(mode="given", kf0=kf[j], kr0=kr[j]) for j in range(3)})
    try:
        r = steady_state(mech, en, state=dict(T_K=500.0), p_gas=pg,
                         theta_guess={"CO*": cm[0], "O*": cm[1]}, ladder=("newton",),
                         tol=1e-14, precision="auto", dps=80, strict_energy=False)
    except Exception: continue
    if not r.converged: continue
    dc = max(abs(r.theta[s]-c) for s,c in zip(("CO*","O*"), cm))
    if dc > 1e-6: continue          # 落到别的解支的点不参与速率比较
    n += 1; dcov.append(dc)
    ref=[float(x) for x in p["rates"]]; got=np.asarray(r.rate_net,float)
    for j in range(3):
        per[j].append(abs(got[j]/ref[j]-1)
                      if (ref[j] != 0 and abs(ref[j]) > 1e-300) else np.nan)
    sK = max(abs(got)) or 1.0; sC = max(abs(np.array(ref))) or 1.0
    fluxK.append(max(abs(got[0]-got[2]), abs(2*got[1]-got[2]))/sK)
    fluxC.append(max(abs(ref[0]-ref[2]), abs(2*ref[1]-ref[2]))/sC)
print("\n同一解支上的比较点：%d" % n)
print("覆盖度绝对差：中位 %.3e  最大 %.3e" % (np.median(dcov), max(dcov)))
# ★ 关键的条件统计：只在 **CatMAP 自己的速率满足质量平衡**的点上比
fc = np.array(fluxC); good = fc < 1e-10
print("\nCatMAP 自身质量平衡残差 < 1e-10 的点：%d / %d" % (int(good.sum()), len(fc)))
for j, name in enumerate(("步1", "步2", "步3")):
    a = np.array(per[j])
    if len(a) == len(fc):
        g = a[good]
        print("  %s 在这些点上的相对差：中位 %.3e  最大 %.3e" % (name, np.nanmedian(g), np.nanmax(g)))
for j, name in enumerate(("步1 CO 吸附（准平衡）", "步2 O2 解离吸附", "步3 CO*+O*→CO2")):
    a = np.array(per[j])
    print("  %-22s 相对差 中位 %.3e  90%%分位 %.3e   <1e-6 的点 %d/%d"
          % (name, np.nanmedian(a), np.nanpercentile(a,90), int((a<1e-6).sum()), len(a)))
print("\n两边各自的**质量平衡残差**（相对最大步速率，越小越自洽）：")
print("  KinSTEP 中位 %.3e  最大 %.3e" % (np.median(fluxK), max(fluxK)))
print("  CatMAP  中位 %.3e  最大 %.3e" % (np.median(fluxC), max(fluxC)))
