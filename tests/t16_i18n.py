# -*- coding: utf-8 -*-
r"""T16 ｜ 界面语言：英文模式下**不得出现一个汉字**，中文模式必须原样不动

## 为什么要这道回归
`KINSTEP_LANG=en` 是为 CPC 投稿加的：审稿人要真的把程序跑起来。
但双语是**最容易悄悄烂掉**的东西 —— 以后随便谁加一行 `print("……")`，
英文模式就会冒出半句中文，而且**不会报错**。
所以这里不查"翻得好不好"（那是人读的事），只查两条机器能判的：

  ① 英文模式下，所有命令的输出里**一个 CJK 字符都不能有**；
  ② 中文模式下，若干关键串必须**逐字不变**（证明默认行为没被动过）。

②  这一条同样要紧：75 项既有回归都逐字比对中文输出，
   如果有人"顺手"把中文也改了，这里会先炸，而不是等到那 75 项里去猜是谁改的。

用法: python t16_i18n.py
"""
import os
import re
import subprocess
import sys

import _root, refdata  # noqa: F401  ← 包根入 sys.path + 第三方真值定位

CJK = re.compile(r"[一-鿿　-〿＀-￯]")
PY = sys.executable
ROOT = _root.ROOT
# ★ 临时目录必须落在**纯 ASCII 路径**下。
#   第一版放在 ROOT/notes/_t16_tmp，而 ROOT 本身就是「…/工具/微动力学/…」——
#   CLI 把文件路径回显出来，于是这道"不得有汉字"的闸被**自己的工作目录**打响了。
#   这是误报，不是翻译漏了。用 tempfile 建一个 ASCII 路径。
import tempfile
TMP = tempfile.mkdtemp(prefix="kinstep_t16_")

# 英文模式要扫的命令。覆盖全部 7 个子命令 + 一条会**报错退出**的路径
#（报错文案最容易被漏掉 —— 它平时不出现）
CMDS_EN = [
    ["templates"],
    ["templates", "her"],
    ["drivers"],
    ["new", "--template", "her", "--force", "%(her)s"],
    ["check", "%(her)s"],                 # 会以退出码 3 结束（能量槽位空着）
    ["explain", "%(her)s"],
    ["ask", "hydrogen evolution at pH 0, -0.1 V vs SHE, 25 C"],
    ["new", "%(amm)s", "--force"],
    ["run", "%(amm)s", "--drc", "--no-plot", "--out", "%(out)s"],
    ["check", "%(nope)s"],                # 文件不存在 ⇒ 退出码 2
]


def run(args, lang):
    env = dict(os.environ)
    if lang:
        env["KINSTEP_LANG"] = lang
    else:
        env.pop("KINSTEP_LANG", None)
    env["PYTHONPATH"] = ROOT + os.pathsep + env.get("PYTHONPATH", "")
    p = subprocess.run([PY, "-m", "kinstep.cli"] + args, env=env,
                       capture_output=True, text=True, cwd=TMP, timeout=600)
    return p.returncode, p.stdout + p.stderr


def main():
    os.makedirs(TMP, exist_ok=True)
    subst = dict(her=os.path.join(TMP, "her.kinstep.toml"),
                 amm=os.path.join(TMP, "ammonia.kinstep.toml"),
                 out=os.path.join(TMP, "out"),
                 ask=os.path.join(TMP, "asked.kinstep.toml"),
                 demo=os.path.join(TMP, "demo.kinstep.toml"),
                 scan=os.path.join(TMP, "scan.kinstep.toml"),
                 nope=os.path.join(TMP, "does_not_exist.toml"))
    bad = 0

    print("═══ ① 英文模式：输出里不得有汉字 ═══")
    for raw in CMDS_EN:
        args = [a % subst if "%(" in a else a for a in raw]
        rc, out = run(args, "en")
        hits = CJK.findall(out)
        shown = " ".join(args[:3])
        if hits:
            bad += 1
            uniq = sorted(set(hits))[:12]
            print("  ✗ %-42s 出现 %d 个 CJK 字符：%s"
                  % (shown, len(hits), "".join(uniq)))
            for ln in out.splitlines():
                if CJK.search(ln):
                    print("      | " + ln[:110])
                    break
        else:
            print("  ✓ %-42s 无汉字（退出码 %d，%d 行）"
                  % (shown, rc, len(out.splitlines())))

    print("\n═══ ② 中文模式：关键串必须逐字不变 ═══")
    # 这些串就是既有 75 项回归所依赖的那几条；它们一变，那边也会跟着变
    MUST_ZH = [
        (["templates"], "离线模板库"),
        (["templates"], "★ 红线：模板**只给机理与问题，不给能量数值**。"),
        (["drivers"], "已注册的模态（driver）"),
        (["check", subst["her"]], "① 机理"),
        (["check", subst["her"]], "四道闸全过"),
    ]
    for args, needle in MUST_ZH:
        rc, out = run(args, None)
        if needle in out:
            print("  ✓ %-28s 含「%s」" % (args[0], needle[:26]))
        else:
            bad += 1
            print("  ✗ %-28s **缺**「%s」—— 中文默认被改动了" % (args[0], needle[:26]))

    print("\n═══ ③ set_lang / 环境变量的优先级 ═══")
    sys.path.insert(0, ROOT)
    from kinstep import i18n
    i18n.set_lang("en")
    a = i18n.t("中", "en")
    i18n.set_lang("zh")
    b = i18n.t("中", "en")
    i18n.reset_lang()
    ok = (a == "en" and b == "中")
    print("  %s set_lang('en')→%r  set_lang('zh')→%r" % ("✓" if ok else "✗", a, b))
    bad += 0 if ok else 1
    # 没设 / 设成奇怪的值 ⇒ 一律退回中文
    for v in ("", "fr", "EN_US", "english"):
        i18n.set_lang(v)
        got = i18n.lang()
        want = "en" if v.lower().startswith("en") else "zh"
        if got != want:
            bad += 1
            print("  ✗ set_lang(%r) 得到 %r，应为 %r" % (v, got, want))
    i18n.reset_lang()
    print("  ✓ 未知语言值一律退回中文")


    print("\n\u2550\u2550\u2550 \u2463 \u82f1\u6587\u6a21\u5f0f\uff1a**\u751f\u6210\u51fa\u6765\u7684\u6587\u4ef6**\u91cc\u4e5f\u4e0d\u5f97\u6709\u6c49\u5b57 \u2550\u2550\u2550")
    # \u2605 \u8fd9\u4e00\u8282\u662f\u88ab\u771f bug \u903c\u51fa\u6765\u7684\uff1a\u7b2c\u4e00\u7248\u53ea\u626f stdout\uff0c\u4e8e\u662f
    #   `new` \u7684\u56de\u663e\u662f\u82f1\u6587\u3001\u800c\u5199\u51fa\u6765\u7684 .toml \u6574\u4e2a\u662f\u4e2d\u6587\u3002
    #   \u800c\u5ba1\u7a3f\u4eba\u8981\u6539\u7684\u6b63\u662f\u90a3\u4efd\u6587\u4ef6\uff0c\u4e0d\u662f\u5c4f\u5e55\u56de\u663e\u3002
    #   ★★ 第二轮又漏了一条：最先只测了带 `--template` 的路径，
    #   而 `new <path>`（**不带** --template）走的是 `inputs.template()`，
    #   那是另一份硬编码骸架，当时还有 84 个汉字。
    #   闸的覆盖面不全，和没有闸是一样的。四条路径现在全测。
    FILE_CASES = [
        (["new", "--template", "her", "--force", "%(her)s"], subst["her"]),
        (["new", "--template", "ammonia", "--force", "%(amm)s"], subst["amm"]),
        (["new", "--force", "%(demo)s"], subst["demo"]),
        (["new", "--scan", "--force", "%(scan)s"], subst["scan"]),
        (["ask", "hydrogen evolution at pH 0, -0.1 V vs SHE, 25 C",
          "--write", "%(ask)s"], subst["ask"]),
    ]
    for raw, path in FILE_CASES:
        args = [a % subst if "%(" in a else a for a in raw]
        run(args, "en")
        try:
            text = open(path, encoding="utf-8").read()
        except OSError as e:
            bad += 1
            print("  \u2717 %-42s \u6ca1\u751f\u6210\u6587\u4ef6\uff1a%s" % (" ".join(args[:3]), e))
            continue
        hits = CJK.findall(text)
        if hits:
            bad += 1
            print("  \u2717 %-42s \u6587\u4ef6\u91cc\u6709 %d \u4e2a CJK\uff1a%s"
                  % (" ".join(args[:3]), len(hits), "".join(sorted(set(hits))[:12])))
            for ln in text.splitlines():
                if CJK.search(ln):
                    print("      | " + ln[:110])
                    break
        else:
            print("  \u2713 %-42s \u6587\u4ef6\u65e0\u6c49\u5b57\uff08%d \u884c\uff09"
                  % (" ".join(args[:3]), len(text.splitlines())))

    print("\n\u2550\u2550\u2550 \u2464 `new` \u8bf4\u7684\u69fd\u4f4d\u6570 \u5fc5\u987b\u7b49\u4e8e `check` \u62a5\u7684\u7f3a\u53e3\u6570 \u2550\u2550\u2550")
    # \u2605 HER \u4e0a\u66fe\u7ecf `new` \u8bf4 3 \u4e2a\u3001`check` \u8bf4\u7f3a 6 \u9879 \u2014\u2014 \u4e24\u6761\u547d\u4ee4\u81ea\u5df1\u6253\u81ea\u5df1\u3002
    for tpl_key, fp in (("her", subst["her"]), ("ammonia", subst["amm"])):
        _, o_new = run(["new", "--template", tpl_key, "--force", fp], "en")
        _, o_chk = run(["check", fp], "en")
        m_new = re.search(r"(\d+) energy slots", o_new)
        m_chk = re.search(r"\[(\d+) value\(s\) missing\]", o_chk)
        if not (m_new and m_chk):
            bad += 1
            print("  \u2717 %-10s \u6ca1\u89e3\u6790\u51fa\u6570\u5b57 new=%r check=%r"
                  % (tpl_key, bool(m_new), bool(m_chk)))
            continue
        a, b = int(m_new.group(1)), int(m_chk.group(1))
        if a == b:
            print("  \u2713 %-10s new \u8bf4 %d \u4e2a\u69fd\u4f4d\uff0ccheck \u7f3a %d \u9879\uff0c\u4e00\u81f4" % (tpl_key, a, b))
        else:
            bad += 1
            print("  \u2717 %-10s new \u8bf4 %d \u4e2a\u69fd\u4f4d\uff0c\u4f46 check \u7f3a %d \u9879" % (tpl_key, a, b))

    print("\n\u2550\u2550\u2550 \u2465 \u6bcf\u4e2a\u7f3a\u5931\u503c\u90fd\u5f97\u6709\u63d0\u793a\uff08\u4e0d\u5f97\u6709 no hint\uff09 \u2550\u2550\u2550")
    # \u2605 `new` \u627f\u8bfa\u300c\u6bcf\u4e2a\u69fd\u4f4d\u90fd\u6807\u4e86\u8981\u7b97\u4ec0\u4e48\u300d\uff0c\u90a3\u5c31\u4e0d\u80fd\u6709\u4e00\u6761\u5199\u7740 no hint\u3002
    for tpl_key, fp in (("her", subst["her"]), ("ammonia", subst["amm"])):
        run(["new", "--template", tpl_key, "--force", fp], "en")
        _, o = run(["check", fp], "en")
        n_nohint = o.count("no hint available")
        if n_nohint:
            bad += 1
            print("  \u2717 %-10s \u6709 %d \u6761\u300cno hint available\u300d" % (tpl_key, n_nohint))
        else:
            print("  \u2713 %-10s \u6bcf\u4e2a\u7f3a\u5931\u503c\u90fd\u5e26\u63d0\u793a" % tpl_key)

    print("\n" + ("✓ T16 全过" if bad == 0 else "✗ T16 有 %d 项不过" % bad))
    return 0 if bad == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
