# -*- coding: utf-8 -*-
"""把包根塞进 sys.path。

★ 为什么要有这个文件：测试脚本原先把**开发机的绝对路径**写死在 sys.path 里
  —— 这台机器之外一律跑不起来。CPC 的代码存档会被审稿人下载到任意路径下，
  写死绝对路径等于交付一个只有作者能跑的回归套件。

用法：测试脚本开头 `import _root` 即可（Python 会自动把脚本所在目录加入
sys.path，所以同目录的本模块总能被 import 到，不需要再写任何路径）。
"""
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
