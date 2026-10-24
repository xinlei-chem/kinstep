#!/usr/bin/env python3
"""重 build 后必跑：检查行内公式有没有原样漏进 Word/PDF。

════ 为什么需要这个脚本 ════
2026-10-24 投稿前发现渲染出来的 PDF 里残留了原始 LaTeX，共 3 处，两类原因：

  ① **粗体跨换行且内含 $...$** —— 粗体标记本身被正确去掉了（PDF 里看不到 `**`），
     但里面的 `$\\nu_e$` / `$1-z$` / `$\\ln z$` 原样漏出。
     同一行里 `$\\nu$` 在粗体之外却渲染正常 ⇒ 是"粗体 run 内不再做公式替换"。
  ② **表格单元格内的公式含 `|`** —— `$\\sigma=|r|\\varepsilon/(1-z)$` 里的竖线
     被当成列分隔符，该单元格在 PDF 里裂成三列。

两处都不会报错、不会 build 失败，只在**渲染结果**里看得见 ⇒ 必须查成品，不能只查源文件。
所以本脚本两头都查：源文件查"会出问题的写法"，PDF 查"是否真的漏了"。

用法:
    python3 check_math_render.py           # 查源文件 + （若有 PDF）查渲染结果
退出码 0 = 干净；1 = 有问题。
"""
from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SOURCES = ["manuscript.md", "SI.md"]
PDFS = ["KinSTEP_CPC_manuscript.pdf", "KinSTEP_CPC_SI.pdf"]

# PDF 里出现这些就说明公式没转：常见的 LaTeX 命令，且必须带反斜杠。
# 不要用裸词（"times" 会命中 sometimes、"frac" 会命中 fraction）—— 这个误报踩过。
TEX_CMD = re.compile(r"\\(?:ln|log|exp|theta|nu|mu|sigma|varepsilon|epsilon|Delta|delta|"
                     r"alpha|beta|gamma|lambda|times|cdot|frac|mathrm|mathsf|text|sim|"
                     r"approx|lvert|rvert|langle|rangle|partial|sum|prod|to|infty)\b")


def scan_sources() -> int:
    bad = 0
    for fn in SOURCES:
        p = HERE / fn
        if not p.exists():
            print("  ! 找不到 %s" % fn)
            continue
        src = p.read_text(encoding="utf-8")

        # ① 含公式的粗体（跨行的尤其危险）
        for m in re.finditer(r"\*\*(.+?)\*\*", src, flags=re.S):
            inner = m.group(1)
            if "$" in inner:
                ln = src[: m.start()].count("\n") + 1
                print("  ✗ %s:%d 粗体内含公式%s：%s"
                      % (fn, ln, "（跨行）" if "\n" in inner else "",
                         inner.replace("\n", " ")[:70]))
                bad += 1

        # ② 表格单元格内公式含 |
        for i, line in enumerate(src.splitlines(), 1):
            if line.lstrip().startswith("|"):
                for m in re.finditer(r"\$[^$]*\$", line):
                    if "|" in m.group(0):
                        print("  ✗ %s:%d 表格内公式含竖线：%s" % (fn, i, m.group(0)))
                        bad += 1
    return bad


def scan_pdfs() -> int:
    bad = 0
    for fn in PDFS:
        p = HERE / fn
        if not p.exists():
            print("  ⊘ %s 不存在（先 build 再查渲染结果）" % fn)
            continue
        try:
            txt = subprocess.run(["pdftotext", "-layout", str(p), "-"],
                                 capture_output=True, text=True, check=True).stdout
        except (FileNotFoundError, subprocess.CalledProcessError) as exc:
            print("  ! 无法解析 %s：%s" % (fn, exc))
            bad += 1
            continue

        for i, line in enumerate(txt.splitlines(), 1):
            cmds = TEX_CMD.findall(line)
            dollars = line.count("$")
            if cmds or dollars:
                print("  ✗ %s 第%d行残留 LaTeX（$×%d, 命令 %s）：%s"
                      % (fn, i, dollars, ",".join(sorted(set(cmds))) or "无", line.strip()[:76]))
                bad += 1
    return bad


def main() -> int:
    print("① 源文件：会导致公式漏出的写法")
    s = scan_sources()
    print("   干净" if not s else "   共 %d 处" % s)

    print("\n② 渲染结果：PDF 里是否真的残留原始 LaTeX")
    r = scan_pdfs()
    print("   干净" if not r else "   共 %d 处" % r)

    total = s + r
    print("\n%s" % ("✅ 公式渲染检查通过" if not total else "⛔ 共 %d 处问题，改完重 build 再查" % total))
    return 0 if not total else 1


if __name__ == "__main__":
    sys.exit(main())
