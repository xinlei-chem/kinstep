# -*- coding: utf-8 -*-
r"""静态扫描：源码里有没有 Python 3.10/3.11 才有的构造。

## 为什么要有它

`pyproject.toml` 把最低版本从 3.11 放宽到 3.9（`tomllib` 用 `tomli` 回退），
但**开发机上没有 3.9 解释器**，没法真跑。声明一个自己没验过的版本下限是不负责任的，
所以至少把"语言层面有没有用到 3.10+ 的东西"用静态扫描查实，并在文档里如实写明
验证到哪一档：**3.11/3.12 是实跑的，3.9/3.10 只过了这一关**。

查什么：
  · `match` / `except*` 语句（3.10 / 3.11）
  · `typing` 里 3.10 之后才有的名字（TypeAlias、ParamSpec、Self …）
  · `zip(strict=)`、`dataclass(slots=/kw_only=)`、`itertools.pairwise`、`int.bit_count`
  · 运行期求值的 `X | Y` 注解 —— 只在**没有** `from __future__ import annotations`
    的模块里才算问题

查不到的（已知局限，不假装覆盖）：
  · 标准库函数的新参数（如 3.9 没有的某个 kwarg）
  · 带括号的多项 with（3.9 的解析器不支持，但 AST 层面看不出来）
  ⇒ 真要背书 3.9，正确做法是在 CI 里装一个 3.9 跑全套；本脚本只是下限。

    python tests/py_compat_scan.py        # 0 = 没发现问题
"""
import ast
import glob
import os
import sys

POST39_TYPING = {"TypeAlias", "ParamSpec", "Concatenate", "TypeGuard", "Self",
                 "Never", "LiteralString", "Required", "NotRequired",
                 "assert_type", "dataclass_transform"}
POST39_ATTRS = {"pairwise", "bit_count"}

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def scan(path):
    out = []
    tree = ast.parse(open(path, encoding="utf-8").read(), path)
    has_future = any(
        isinstance(n, ast.ImportFrom) and n.module == "__future__"
        and any(a.name == "annotations" for a in n.names) for n in tree.body)
    for n in ast.walk(tree):
        cname = n.__class__.__name__
        if cname in ("Match", "TryStar"):
            out.append((n.lineno, "%s 语句（3.10/3.11）" % cname))
        if isinstance(n, ast.ImportFrom) and n.module == "typing":
            for a in n.names:
                if a.name in POST39_TYPING:
                    out.append((n.lineno, "typing.%s（3.10+）" % a.name))
        if isinstance(n, ast.Attribute) and n.attr in POST39_ATTRS:
            out.append((n.lineno, ".%s（3.10+）" % n.attr))
        if isinstance(n, ast.Call):
            for kw in n.keywords:
                if kw.arg == "strict" and getattr(n.func, "id", "") == "zip":
                    out.append((n.lineno, "zip(strict=)（3.10+）"))
                if kw.arg in ("slots", "kw_only"):
                    out.append((n.lineno, "dataclass(%s=)（3.10+）" % kw.arg))
        if not has_future:
            ann = getattr(n, "annotation", None) or getattr(n, "returns", None)
            if isinstance(ann, ast.BinOp) and isinstance(ann.op, ast.BitOr):
                out.append((n.lineno, "注解里的 X|Y 且本模块无 from __future__"
                                      " import annotations（3.10+）"))
    return out


def main():
    files = []
    for pat in ("kinstep/*.py", "kinstep/drivers/*.py", "tests/*.py",
                "examples/*.py", "paper/*.py"):
        files += sorted(glob.glob(os.path.join(ROOT, pat)))
    bad = 0
    for f in files:
        for lineno, why in scan(f):
            print("  ⚠ %s:%d  %s" % (os.path.relpath(f, ROOT), lineno, why))
            bad += 1
    print("扫描 %d 个文件，命中 %d 处%s"
          % (len(files), bad, "" if bad else "（未发现 3.10+ 专属构造）"))
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
