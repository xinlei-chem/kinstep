# -*- coding: utf-8 -*-
r"""T14 验证：离线模板库。

★ 红线（台账里定过的）：**只给「结构」和「问题」，不给「数字」**
----------------------------------------------------------------
  · 机理是公开知识 —— 写错了肉眼能看出来，而且本来就发表在文献里
  · 能量是私有数据 —— 填错了**只会得到一条看着很合理的火山曲线**

本测试把这条红线做成**机器检查**：模板里出现任何能量数值，直接失败。不靠自觉。
"""
import sys, os, re, tempfile, warnings
import _root, refdata  # noqa: F401  ← 包根入 sys.path + 第三方真值定位
warnings.filterwarnings("ignore")

from kinstep import templates as TP
from kinstep.inputs import load, InputError
from kinstep.provenance import ENERGY_KEYS
from kinstep.driver import make_driver, available_drivers

fail = []

print("══ ① 每个模板的机理必须通过四道闸 ══")
for k in TP.list_templates():
    t = TP.get(k)
    try:
        m = t.parse()
        print("   %-16s ✓ %d 步 %2d 物种  [%s]" % (k, m.n_step, len(m.species_order), t.modality))
    except Exception as e:
        fail.append((k, "机理不合法: %s" % str(e).splitlines()[0]))
        print("   %-16s ❌ %s" % (k, str(e).splitlines()[0]))

print("\n══ ② ★★ 红线机器检查：模板里不许出现任何能量数值 ══")
# 能量槽位的键名 + 形如 "Ga_fwd = 1.23" 的赋值，一律不许在模板里出现
pat = re.compile(r"(%s)\s*[=:]\s*[-+0-9]" % "|".join(sorted(ENERGY_KEYS)))
for k in TP.list_templates():
    t = TP.get(k)
    blobs = [t.mechanism, t.to_toml(), repr(t.conditions), repr(t.drivers),
             repr(t.energy_hints), " ".join(t.notes)]
    hits = [m.group(0) for b in blobs for m in pat.finditer(b)]
    # Template 的数据字段里也不许直接带能量键
    for d in (t.conditions, *(t.drivers.values())):
        for key in d:
            if key in ENERGY_KEYS:
                hits.append("字段 %s" % key)
    if hits:
        fail.append((k, "出现能量数值: %s" % hits))
        print("   %-16s ❌ %s" % (k, hits))
    else:
        print("   %-16s ✓ 无能量数值" % k)
print("   判据说明：匹配 %s 这些键后面跟数字的写法" % sorted(ENERGY_KEYS))

print("\n══ ②b 反向验证：故意塞一个能量数值进去，检查必须报警 ══")
probe = TP.Template(key="_probe", title="探针", modality="thermal",
                    mechanism="CO + * <-> CO*",
                    conditions=dict(T_K=500.0, Ga_fwd=1.23))   # ← 违规
hits = [key for key in probe.conditions if key in ENERGY_KEYS]
print("   塞入 Ga_fwd=1.23 ⇒ %s" % ("✓ 被抓到：%s" % hits if hits else "❌ 没抓到"))
if not hits:
    fail.append(("_probe", "红线检查失效"))

print("\n══ ③ 每个模板都要有出处 ══")
for k in TP.list_templates():
    t = TP.get(k)
    if not t.reference.strip():
        fail.append((k, "没有出处"))
        print("   %-16s ❌ 没有出处" % k)
    else:
        print("   %-16s ✓ %s%s" % (k, t.reference[:52],
                                   "  DOI %s" % t.doi if t.doi else ""))

print("\n══ ④ 声明了 driver 的模板，机理里必须真有对应的逐步参数 ══")
for k in TP.list_templates():
    t = TP.get(k)
    if not t.drivers:
        continue
    m = t.parse()
    for dn in t.drivers:
        try:
            d = make_driver(dn)
            d.validate(m)
            print("   %-16s ✓ driver=%s 的逐步参数齐全" % (k, dn))
        except Exception as e:
            fail.append((k, "driver %s: %s" % (dn, str(e).splitlines()[0])))
            print("   %-16s ❌ driver=%s: %s" % (k, dn, str(e).splitlines()[0][:60]))

print("\n══ ⑤ ★ 生成的 TOML：load() 应当成功，缺能量由 check() 报 ══")
print("   （职责划分：load 负责「文件读得懂吗」，check 负责「还缺什么」——")
print("    这正是 `kinstep check` 存在的意义）")
from kinstep.provenance import UnverifiedDataError
for k in TP.list_templates():
    t = TP.get(k)
    with tempfile.TemporaryDirectory() as td:
        fp = os.path.join(td, k + ".kinstep.toml")
        open(fp, "w", encoding="utf-8").write(t.to_toml())
        try:
            sp = load(fp)
        except Exception as e:
            fail.append((k, "load 失败: %s" % str(e).splitlines()[0][:50]))
            print("   %-16s ❌ load 失败：%s" % (k, str(e).splitlines()[0][:52]))
            continue
        try:
            sp.energetics.check()
            fail.append((k, "check 没报缺能量"))
            print("   %-16s ❌ check 居然通过了（能量明明没填）" % k)
        except UnverifiedDataError as e:
            n = str(e).count("—— ")
            print("   %-16s ✓ load 通过；check 报出缺 %d 项能量" % (k, n))
        except Exception as e:
            fail.append((k, "check 抛了 %s" % type(e).__name__))
            print("   %-16s ❌ check 抛 %s" % (k, type(e).__name__))

print("\n══ ⑥ pH 双重计数的提醒：显式写 H+ 的模板必须注明用 SHE ══")
for k in TP.list_templates("electro"):
    t = TP.get(k)
    has_hp = "H+" in t.mechanism
    uref = t.drivers.get("electro", {}).get("U_ref", "")
    noted = any("SHE" in n for n in t.notes)
    if has_hp:
        ok = (uref == "SHE")
        print("   %-16s 机理含 H+，U_ref=%s %s%s"
              % (k, uref or "(未给)", "✓" if ok else "❌",
                 "  且 notes 里已提醒" if noted else ""))
        if not ok:
            fail.append((k, "含 H+ 却没用 SHE 标度 ⇒ pH 会被重复计算"))

print("\n══ ⑦ ★★ 最强的检查：每一步的正向速率都必须非零 ══")
print("   要抓的是「活度漏给」这类缺陷（HER 模板漏了 c_liq ⇒ 全零，前六项检查全放过了）。")
print("   做法：给所有表面物种一个相同的占位覆盖度，直接看逐步 r_fwd ——")
print("   **不做完整稳态求解**：那既慢（9 步的 methanation 卡了 4 分钟），")
print("   又会被「能量取值是否合理」干扰，而这里要测的根本不是能量。")
import numpy as np
from kinstep.analyze import energetics_from_arrays
from kinstep.rates import compile_rates, apply_state, activities, rates_factored
from kinstep.driver import DriverStack, make_driver
for k in TP.list_templates():
    t = TP.get(k)
    m = t.parse()
    en = energetics_from_arrays(m, np.full(m.n_step, 0.8),
                                np.linspace(-0.3, 0.3, m.n_step))
    st = dict(t.conditions)
    pg = st.pop("p_gas", {}) or {}
    cl = st.pop("c_liq", {}) or {}
    for dn, cfg in t.drivers.items():
        st.update(cfg)
    st["photo_carrier_mode"] = "given_qfl"
    st.setdefault("qfl_split_eV", 0.1)
    drv = None
    if t.drivers:
        try:
            drv = DriverStack([make_driver(dn) for dn in t.drivers], warn_order=False)
        except Exception:
            drv = None
    try:
        rc = apply_state(compile_rates(en, strict_energy=False),
                         float(st["T_K"]), drv, st)
        n_surf = max(len(m.surface_species), 1)
        th = {nm: 0.9 / n_surf for nm in m.surface_species}
        act = activities(m, th, pg, cl)
        rf, _, _ = rates_factored(m, rc, act)
        zero = [i + 1 for i, v in enumerate(rf) if not (v > 0)]
        if zero:
            bad = [nm for i, nm in enumerate(m.species_order)
                   if act[i] <= 0 and m.species[nm].phase in ("gas", "liq")]
            fail.append((k, "第 %s 步 r_fwd = 0，活度为 0 的物种：%s" % (zero, bad)))
            print("   %-16s X 第 %s 步 r_fwd=0；活度为 0 的物种 %s" % (k, zero, bad))
        else:
            print("   %-16s OK %d 步全部 r_fwd > 0（最小 %.3e）"
                  % (k, m.n_step, float(np.min(rf))))
    except Exception as exc:
        fail.append((k, "算不出速率: %s" % str(exc).splitlines()[0][:50]))
        print("   %-16s X %s" % (k, str(exc).splitlines()[0][:58]))

print("\n══ 判定 ══")
if fail:
    print("   ❌ %d 项不合格：" % len(fail))
    for k, why in fail:
        print("      %-16s %s" % (k, why))
    sys.exit(1)
print("   ✓ %d 个模板全部合格（机理合法 / 无能量数值 / 有出处 / driver 参数齐 / "
      "读写正确 / pH 标度一致 / **占位能量下能跑出非零速率**）" % len(TP.list_templates()))
