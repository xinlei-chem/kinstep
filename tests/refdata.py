# -*- coding: utf-8 -*-
"""第三方参考数据（"真值"）的定位，以及缺失时的优雅跳过。

★ 为什么不把这些数据打包进来
----------------------------
CATKINAS 是闭源发行物（MATLAB pcode `.p`），其手册 §2.2 只授予
"使用时必须引用"，**没有授予再分发许可**。所以它自带算例的输出
`Examples/result_curve/data1.mat` 等文件不随 KinSTEP 分发。

CatMAP 是 GPLv3，再分发是允许的，但它的上游完整克隆（9.3 MB）对回归
没有用处——本套回归的真值一个字节都不来自 CatMAP（见 ATTRIBUTION.md
第三节），所以同样不进发行物。

依赖这些数据的回归（T7/T8/T9/T9b/T10/T11）在数据缺失时**以退出码 77
跳过**，不算失败；其余 8 套与 6 个案例、CLI、错误路径检查完全自足。
论文里逐数字报告的对比值已随 SI 存档（paper/SI/*.npz、*.dat），
所以结论可核，只是重跑那几套需要自备 CATKINAS。
"""
import os
import sys

import numpy as np

SKIP = 77          # 和 automake 的约定一致：77 = 跳过，不是失败
ENV = "KINSTEP_CATKINAS_DIR"

# ★ 这里**不留开发机路径作为兜底**：存档要公开发布，兜底路径里带着作者的用户名
#   与中文目录，对别人毫无用处，还会让"数据不在"退化成"路径对不上"的误判。
#   位置只认环境变量；没设就按 77 跳过，并把设法打印出来。


def catkinas(rel):
    """返回 CATKINAS 算例文件的绝对路径；没有就打印办法并以 77 退出。"""
    base = os.environ.get(ENV) or ""
    path = os.path.join(base, rel)
    if os.path.exists(path):
        return path
    sys.stderr.write(
        "⊘ 跳过本套回归：缺第三方参考数据\n"
        "    需要 : <CATKINAS>/%s\n"
        "    找过 : %s\n"
        "  CATKINAS 为闭源发行物，未授再分发许可，故不随本包提供。\n"
        "  自备一份后设环境变量再跑：\n"
        "    export %s=/path/to/CATKINAS\n" % (rel, path, ENV))
    sys.exit(SKIP)


# ── 参考值的两条取法 ────────────────────────────────────────────────
#
# ★ 为什么不再直接 loadmat：存档里没有 `.mat`，那 14 项断言与 Fig.2 就全跳过，
#   而审稿人是会把存档跑一遍的。论文一贯的口径是「登载比对值 ≠ 再分发软件」
#   （RELEASE_NOTES 第二节，SI 里本来就存了 401 点的参考速率与覆盖度），
#   所以把**断言真正用到的那些数组**抽出来随包分发，存档因此自足。
#
# 优先级：
#   ① 设了 KINSTEP_CATKINAS_DIR ⇒ 直接读发行物原件（权威，可完全重算）；
#   ② 否则读随包的 catkinas_refdata.npz（由 make_catkinas_refdata.py 生成）；
#   ③ 两者都没有 ⇒ 按 77 跳过。
#
# 两条路返回同一种东西：一个 {变量名: ndarray} 的 dict，调用方不必知道来源。

_NPZ = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "catkinas_refdata.npz")

_REL = {"curve": "Examples/result_curve/data1.mat",
        "map":   "Examples/result_map/data1.mat"}


def catkinas_vars(which):
    """取一组 CATKINAS 参考数组。which ∈ {'curve', 'map'}。

    返回 (dict, source)，source 是 'mat' 或 'npz'，供回归在日志里标明口径。
    """
    if which not in _REL:
        raise ValueError("which 只能是 %s，收到 %r" % (sorted(_REL), which))

    if os.environ.get(ENV):
        import scipy.io as sio
        d = sio.loadmat(catkinas(_REL[which]))
        return {k: v for k, v in d.items() if not k.startswith("__")}, "mat"

    if os.path.exists(_NPZ):
        pre = which + "__"
        z = np.load(_NPZ)
        out = {k[len(pre):]: z[k] for k in z.files if k.startswith(pre)}
        out.pop("md5", None)
        if out:
            return out, "npz"

    sys.stderr.write(
        "⊘ 跳过本套回归：既没有随包的参考值存档，也没有 CATKINAS 原件\n"
        "    找过 : %s\n"
        "  自备 CATKINAS 后：export %s=/path/to/CATKINAS\n" % (_NPZ, ENV))
    sys.exit(SKIP)
