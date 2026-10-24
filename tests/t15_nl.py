# -*- coding: utf-8 -*-
r"""T15 验证：自然语言层（规则翻译 + LLM 后端闸）。

★ 红线：**只许产出「结构」和「问题」，不许产出「数字」**
---------------------------------------------------------
更准确地说，数字分两类：
  · **用户自己说出来的** —— 可以抽取，标 user_given（本来就是他的数据）
  · **系统要去猜的**（任何能量、任何没说的参数）—— 一律变成问题，标 REQUIRED

本测试把这条红线做成机器检查，并且**反向验证**：
喂进含能量数值的输入、以及一个恶意的 LLM 后端，看闸是否拦得住。
"""
import sys, warnings
import _root, refdata  # noqa: F401  ← 包根入 sys.path + 第三方真值定位
warnings.filterwarnings("ignore")

from kinstep import nl
from kinstep.provenance import Val, ENERGY_KEYS, gate, UnverifiedDataError
from kinstep.templates import LIBRARY

fail = []

print("══ ① 规则翻译：各种说法都要落到对的体系上 ══")
CASES = [
    ("我想算 Pt 上的析氢反应", "her"),
    ("hydrogen evolution reaction on MoS2", "her"),
    ("Fe 催化剂 673 K 下的合成氨", "ammonia"),
    ("ammonia synthesis volcano plot", "ammonia"),
    ("IrO2 的析氧过电位", "oer"),
    ("燃料电池阴极的氧还原", "orr"),
    ("TiO2 光催化水分解", "photo_water"),
    ("超声空化促进的水解离", "sono_radical"),
    ("Pd 上的 CO 氧化", "co_oxidation"),
    ("水煤气变换反应 Cu/ZnO", "wgs"),
]
for text, want in CASES:
    r = nl.understand(text)
    ok = r.template_key == want
    if not ok:
        fail.append((text, "匹配到 %r，期望 %r" % (r.template_key, want)))
    print("   %-32s → %-14s %s" % (text[:32], r.template_key, "✓" if ok else "❌ 期望 " + want))

print("\n══ ② 数字抽取 + 单位换算必须显式记录 ══")
r = nl.understand("400 摄氏度，50 bar，电位 -0.15 V vs SHE，pH=1，超声 20 kHz")
for e in r.extracted:
    print("   " + str(e))
chk = {e.key: e.value for e in r.extracted}
exp = {"T_K": 673.15, "p_total_bar": 50.0, "U_elec": -0.15, "pH_val": 1.0,
       "sono_freq_kHz": 20.0}
for k, v in exp.items():
    got = chk.get(k)
    ok = got is not None and abs(got - v) < 1e-9
    if not ok:
        fail.append((k, "抽到 %r，期望 %r" % (got, v)))
    print("   %-16s %s" % (k, "✓" if ok else "❌ 抽到 %r 期望 %r" % (got, v)))

print("\n══ ③ ★★ 红线：产出里不许有任何「猜出来的」能量 ══")
for text, _ in CASES:
    r = nl.understand(text)
    conds = r.conditions()
    bad = [k for k in conds if k in ENERGY_KEYS]
    if bad:
        fail.append((text, "产出里出现能量键 %s" % bad))
        print("   %-32s ❌ %s" % (text[:32], bad))
# 反向验证：用户原文里就带能量数字，也不许被当成能量抽走
r2 = nl.understand("合成氨，Ga_fwd = 1.31 eV，dG_rxn = 0.4 eV，673 K")
c2 = r2.conditions()
bad2 = [k for k in c2 if k in ENERGY_KEYS]
print("   反向验证：原文里写了 Ga_fwd=1.31 ⇒ 抽到的键 = %s" % sorted(c2))
print("   %s" % ("✓ 能量没被抽走（仍是问题）" if not bad2 else "❌ 能量被当成已知了：%s" % bad2))
if bad2:
    fail.append(("反向验证", "能量被抽走 %s" % bad2))
has_q = any("能量" in q for q in r2.questions)
print("   %s" % ("✓ 且仍把能量列为问题" if has_q else "❌ 没把能量列为问题"))
if not has_q:
    fail.append(("反向验证", "没把能量列为问题"))

print("\n══ ④ ★ 一致性检查：pH 双重计数必须拦住 ══")
r3 = nl.understand("析氢反应 pH=0 电位 -0.1 V vs RHE")
hit = any("重复计算" in w for w in r3.warnings)
print("   机理含 H+ 又说 RHE ⇒ %s" % ("✓ 已警告" if hit else "❌ 没警告"))
if not hit:
    fail.append(("pH", "双重计数没警告"))
r4 = nl.understand("析氢反应 pH=0 电位 -0.1 V vs SHE")
print("   改成 SHE ⇒ %s" % ("✓ 不再警告" if not any("重复计算" in w for w in r4.warnings)
                            else "❌ 误报"))
r5 = nl.understand("光催化的合成氨")
print("   模态不符（说光催化却匹配到热催化体系）⇒ %s"
      % ("✓ 已提醒" if any("请确认没选错" in w for w in r5.warnings) else "❌ 没提醒"))

print("\n══ ⑤ LLM 后端：默认关，且要显式确认才能装 ══")
print("   默认状态 backend_enabled() = %s %s"
      % (nl.backend_enabled(), "✓" if not nl.backend_enabled() else "❌ 居然默认开着"))
if nl.backend_enabled():
    fail.append(("backend", "默认开着"))
try:
    nl.set_backend(lambda t: {}, )        # 不加 confirm
    print("   ❌ 没有二次确认也能装")
    fail.append(("backend", "没二次确认也能装"))
except PermissionError as e:
    print("   ✓ 不加 confirm=True 装不上：%s" % str(e).splitlines()[0][:50])

print("\n══ ⑥ ★★ 恶意后端：返回能量、返回垃圾字段，都必须被过滤 ══")
def evil_backend(text):
    return {"template_key": "her",
            "Ga_fwd": 0.42,                    # ← 想塞能量
            "dG_rxn": -0.1,                    # ← 想塞能量
            "T_K": 9999,                       # ← 想塞工况数值
            "mechanism": "A + * <-> B*",       # ← 元素不配平，应被机理闸刷掉
            "questions": ["你的催化剂是什么？"],
            "secret_field": "rm -rf /",        # ← 白名单之外
            "modality": "electro"}
kept, dropped = nl.sanitize_backend_output(evil_backend("x"))
print("   保留的字段：%s" % sorted(kept))
for d in dropped:
    print("     丢弃：%s" % d)
bad6 = [k for k in kept if k in ENERGY_KEYS] + \
       [k for k in kept if k not in nl.ALLOWED_BACKEND_FIELDS]
num6 = [k for k, v in kept.items() if isinstance(v, (int, float)) and not isinstance(v, bool)]
ok6 = not bad6 and not num6 and "mechanism" not in kept
print("   %s" % ("✓ 能量、数值、垃圾字段、不配平的机理 全部被挡" if ok6
                 else "❌ 漏了：%s %s" % (bad6, num6)))
if not ok6:
    fail.append(("sanitize", "漏了 %s %s" % (bad6, num6)))

print("\n══ ⑦ 最后一道：就算有人硬塞 llm_suggested 的能量，求解器也必须拒绝 ══")
items = [("第1步", "Ga_fwd", Val(0.42, "llm_suggested")),
         ("第1步", "dG_rxn", Val(-0.1, "llm_suggested"))]
try:
    gate(items)
    print("   ❌ 居然放行了")
    fail.append(("gate", "llm_suggested 的能量被放行"))
except UnverifiedDataError as e:
    print("   ✓ 被拒绝：%s" % str(e).splitlines()[0])
    print("     %s" % str(e).splitlines()[1].strip())

print("\n══ 判定 ══")
if fail:
    print("   ❌ %d 项不合格：" % len(fail))
    for k, why in fail:
        print("      %-24s %s" % (str(k)[:24], why))
    sys.exit(1)
print("   ✓ 全部合格：体系匹配 / 单位换算 / **能量永不被猜** / pH 一致性 / "
      "后端默认关 / 恶意产出被过滤 / 求解器最后一道闸")
