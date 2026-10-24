# -*- coding: utf-8 -*-
r"""把回归真正比对的那几个 CATKINAS 参考数组抽出来，存成一个 npz。

## 为什么要有这一步

原先 T7/T8/T9/T9b/T10/T11 直接 `loadmat` CATKINAS 发行物自带的 `data1.mat`。
后果是：**存档里这 14 项断言跑不了（按 77 跳过），Fig.2 也出不来**——
而 CPC 的审稿人是会把存档下载下来跑的，跑不全就是减分项。

论文的立场一向是「**登载我们比对的参考值**，与**再分发那个软件**是两回事」
（`RELEASE_NOTES.md` 第二节）：SI 里已经存了 401 点的参考速率与覆盖度。
这个脚本把同一条口径推到其余几套回归上——**只抽断言真正用到的数组**，
不是把 `.mat` 整个搬过来：

    result_curve/data1.mat   23 个变量 → 只取 14 个
    result_map/data1.mat     23 个变量 → 只取 7 个

## 用法（只有作者需要跑；存档用户直接用打好的 npz）

    export KINSTEP_CATKINAS_DIR=/path/to/CATKINAS
    python tests/make_catkinas_refdata.py

生成 `tests/catkinas_refdata.npz`（压缩后 ~0.2 MB），随包分发。

## 口径

npz 的键名是 `curve__Ean` 这种 `<来源>__<变量>` 形式，免得两个 `.mat` 的同名变量撞车。
每次重新生成都会把来源文件的 md5 一并写进去，便于核对是不是同一份发行物。
"""
import hashlib
import os
import sys

import numpy as np

import _root  # noqa: F401
import refdata

# 每个来源 → 断言真正读到的变量
NEED = {
    "curve": ("Examples/result_curve/data1.mat",
              ["Ean", "G0n", "Eao", "G0o", "kf", "kr", "Z", "y", "Rnet",
               "DRC_B", "DRC_I", "DRC_PC", "DRC_G", "Eapp"]),
    "map":   ("Examples/result_map/data1.mat",
              ["Ean", "G0n", "y", "Z", "Rnet", "params", "kf"]),
}

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   "catkinas_refdata.npz")


def main():
    import scipy.io as sio
    if not os.environ.get(refdata.ENV):
        sys.stderr.write(
            "本脚本只有作者需要跑：它要读 CATKINAS 发行物自带的算例输出。\n"
            "  export %s=/path/to/CATKINAS\n" % refdata.ENV)
        return refdata.SKIP

    blob = {}
    for tag, (rel, keys) in NEED.items():
        path = refdata.catkinas(rel)          # 缺文件时它自己按 77 退出
        with open(path, "rb") as fh:
            blob["%s__md5" % tag] = np.array(
                hashlib.md5(fh.read()).hexdigest())
        d = sio.loadmat(path)
        for k in keys:
            blob["%s__%s" % (tag, k)] = np.asarray(d[k], dtype=float)
        print("  %-6s %-34s 取 %2d / %2d 个变量"
              % (tag, rel, len(keys),
                 len([x for x in d if not x.startswith("__")])))

    np.savez_compressed(OUT, **blob)
    print("\n  -> %s  (%.2f MB)" % (os.path.basename(OUT),
                                    os.path.getsize(OUT) / 1e6))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
