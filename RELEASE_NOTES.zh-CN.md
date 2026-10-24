# KinSTEP 代码存档说明（CPC deposit）

> 这份文件回答审稿人/读者会问的三个问题：**存档里有什么**、**第三方材料怎么处理的**、
> **拿到压缩包能跑到什么程度**。
> 建于 2026-10-23。改动发布口径时必须同时改这里和 `ATTRIBUTION.md`。

---

## 一、许可证

| 对象 | 许可证 |
|---|---|
| KinSTEP 本体（`kinstep/`、`tests/`、`examples/`、`run_all_tests.sh`、文档） | **GPL-3.0-or-later**，全文见 `LICENSE`（md5 `d32239bc…`，GNU 官方原文未改动） |

`pyproject.toml` 里 `license = "GPL-3.0-or-later"` 与 `license-files = ["LICENSE"]` 已对齐。

---

## 二、★ 第三方材料：存档里**不含**哪些，为什么

这是本项目发布前唯一需要判断的问题，结论如下。

### 2.1 CATKINAS 的算例数据 —— 不分发

CATKINAS（Chen、Wang、Hu，华东理工）是 **闭源发行物**（MATLAB pcode `.p`），
从 `catkinas.com/downloads/` 取得。其手册 **§2.2「版权应用」只授予"使用时必须在
出版物中引用"**，**没有授予再分发许可**。

因此其自带算例的输出
`Examples/result_curve/data1.mat`、`Examples/result_map/data1.mat`
**不随本存档分发**。

> 可核性不受影响：论文里逐数字报告的对比值已随 SI 存档
> （`paper/SI/fig2_data.npz` 含 401 点的 CATKINAS 参考速率与覆盖度，
> `paper/SI/fig3d_drc.dat`、`figS1b_dps.dat` 同理）。
> 在科学论文中登载"我们对照的参考值"是正常且必要的做法，与再分发该软件是两回事。

想自己重跑这几套回归：自备 CATKINAS，然后

```bash
export KINSTEP_CATKINAS_DIR=/path/to/CATKINAS     # 解包后含 Examples/ 的那一层
./run_all_tests.sh
```

### 2.2 CatMAP 上游克隆 —— 不分发（但不是许可证原因）

`真值/catmap/` 下有 CatMAP（**GPLv3**，SUNCAT）的完整浅克隆，9.3 MB。
GPLv3 允许再分发，所以**法律上可以带**，但本存档仍然不带，理由有两条：

1. **对回归没有用处**。本套回归的真值**一个字节都不来自 CatMAP**
   （CatMAP 仓库里只有 2 个点的覆盖度参考输出，没有 TOF / 速率 / DRC；详见
   `真值/catmap/PROVENANCE.txt` 第二节）。CatMAP 是**功能对照**对象，不是真值来源。
2. **带着反而容易被误读**成 KinSTEP 衍生自 CatMAP。实际情况是零代码复制
   （`ATTRIBUTION.md` 第〇节），且从未 `import catmap`——没有 import，就不构成
   GPL 意义上的衍生作品。

`ATTRIBUTION.md` 第三节如实列出了与 CatMAP 的功能重合，以及立项依据所在的差异项。

### 2.3 Koper CO stripping 的参照实现 —— 不分发

T5 的参照是 Koper 等 1998 年论文（`[29]`，doi:10.1063/1.477230）的模型，
以及 J. W. Ager 公开的 Jupyter notebook 实现。T5 **不读任何外部文件**：
模型按论文重新实现，真值在脚本里逐数字写明来源。所以 T5 自足。

### 2.4 NIST 热化学数据

Shomate 系数取自 **NIST Chemistry WebBook（SRD 69）原始条目**，逐条标注来源，
不从任何二手手册转抄（`ATTRIBUTION.md` 第二节记录了这条自我更正）。

---

## 三、拿到存档能跑到什么程度

```bash
pip install -e .                 # 只需要 numpy + scipy
pip install -e ".[all]"          # 另加 mpmath（高精度）与 matplotlib（作图）
./run_all_tests.sh
```

回归共 **82 项断言，全部能在干净下载后跑**（原先那个「需第三方软件」的档已经没有了，
因为它要比的**参考值**现在随包分发，见第二节）。剩下的区别只是：5 项依赖可选件。

| 档 | 项数 | 内容 | 干净下载后能否跑 |
|---|---|---|---|
| **A 自足** | **66** | 环境与全包语法扫描；T3b/T4/T5/T6/T12/T13/T14/T15 八套数值回归（含 Koper 逐点右端、稳态 vs 解析闭式解、光催化电荷守恒、超声份额反推）；CLI 20 项；6 个案例；21 项文档数字一致性；5 条错误路径 | **能**。装了可选件 66 项全过；只有 numpy+scipy 时其中"出图 volcano.png"一项按跳过计（matplotlib 属可选件），为 65 过 |
| **B 需第三方参考数据** | **14** | T7/T8/T9/T9b/T10/T11 六套（CATKINAS 401 点描述符扫描、315 点流动反应器图、DRC、热力学链路、高精度阶梯），以及 ③ 段里依赖它们的 8 项关键数字检查 | 自备 CATKINAS 并设 `KINSTEP_CATKINAS_DIR` 后能跑；否则**以退出码 77 跳过**，计为"跳过"而非失败 |

> 跳过不是静默的：汇总行会印「跳过 N 项」，每套跳过的回归也会打印需要哪个文件、
> 去哪里找、怎么设环境变量。

**实测的三种配置**（后两种都是把新打的 sdist 装进空 venv 裸跑的）：

| 配置 | 结果 |
|---|---|
| 装了可选件，且备有参考原件 | **82 过 / 0 失败 / 0 跳过** |
| 装了可选件，没有参考原件 | **82 过 / 0 失败 / 0 跳过** |
| 只有 numpy + scipy | **77 过 / 0 失败 / 5 跳过** |

第三种多跳过的那一项是"出图 volcano.png"——matplotlib 是**可选**件，
没装就不出图本来就是设计如此，所以按跳过计，不报失败。

**插图**：三个脚本分管五张图 ——
`paper/make_figures.py`（Fig.1/2/3/S1）、`paper/make_fig4_maps.py`（Fig.4）、
`paper/make_fig5_cli.py`（Fig.5）。

- Fig.1、Fig.3(a–c)、Fig.S1(a) 完全自足；
- Fig.2 需要 CATKINAS 参考数据才能**重算**（缺数据时整图跳过，数值见 `SI/fig2_data.npz`）；
- Fig.3(d) 与 Fig.S1(b) 在缺数据时**退回 SI 里已存档的数值重绘同一张图**——图可复现，只是不重算；
- Fig.4 完全自足（约 20 分钟算 4992 个点）；加 `--replot` 直接读 `SI/fig4_maps.npz` 重绘，
  调图不必重算；
- Fig.5 完全自足，且**不存任何文本**：每次出图都真跑四条 `KINSTEP_LANG=en` 命令再排版，
  所以图不会和代码脱节。它自带一道闸——输出里出现任何汉字就拒绝出图。
  完整会话另存 `SI/fig5_cli_session.txt`。

---

## 四、存档清单

```
LICENSE                  GPLv3 全文
README.md                安装与使用说明
EXAMPLES.md              六个案例逐个讲解（数字都来自可跑的脚本）
ATTRIBUTION.md           借鉴了什么、必须引用谁
RELEASE_NOTES.md         本文件
pyproject.toml           GPL-3.0-or-later；核心依赖只有 numpy + scipy
MANIFEST.in              ★ prune 真值 —— 把第二节说的东西排除在分发之外
run_all_tests.sh         82 项断言
kinstep/                 包本体 18 个模块 ≈ 6900 行
tests/                   14 套数值回归 ≈ 2200 行；_root.py / refdata.py 负责路径可移植
examples/                6 个可直接跑的案例（ex1…ex6）
notes/DEVLOG.md          逐条的坑：现象 → 根因 → 修复 → 防复发
paper/
  manuscript.md          手稿正文
  SI.md                  Supplementary Information
  make_figures.py        Fig.1 / 2 / 3 / S1
  make_fig4_maps.py      Fig.4（四模态二维彩图；--replot 从 npz 重绘）
  make_fig5_cli.py       Fig.5（命令行界面；每次真跑命令，不存死文本）
  build_docx.py          markdown → Word
  figures/               全部插图（*.pdf 随分发，*.png 不随）
  SI/                    可作图原始数据 *.dat / *.npz + fig5 完整会话
```

> ★ **全部路径都是 ASCII**（2026-10-23 改）。
> 之前目录与文件名用中文（`论文/ 测试/ 案例/ 记录/`、`正文.md`、`脚本_出图.py` 等），
> 对本地开发没问题，但 CPC 的存档要被全球下载解压，**ZIP 不带 UTF-8 标志时文件名会乱码**。
> 文件**内容**仍是中文注释，那不受影响——出问题的只有文件名编码。

`paper/` 里收的是**出图脚本、SI 的可作图原始数据（`*.dat`/`*.npz`）和矢量图（`*.pdf`）**；
`.docx` 手稿与 `.png` 位图不收——手稿走期刊，位图由脚本可重出。

**不在存档里**：`真值/`（第二节）、`notes/FAILED_*.log`、`__pycache__`、`paper/*.docx`、`paper/figures/*.png`。

> `kinstep.egg-info/PKG-INFO` 与 `SOURCES.txt` 会出现在 sdist 里——那是 setuptools
> 自己每次重新生成并强制打包的元数据，`MANIFEST.in` 的 `prune` 管不到，也无害。

打包与自检：

```bash
python -m build --sdist
A=dist/kinstep-*.tar.gz
tar tzf $A | grep -cE "真值|catmap|data1\.mat|\.docx"   # 禁入项，必须是 0
tar tzf $A | grep -cP '[^\x00-\x7F]'                    # 非 ASCII 文件名，必须是 0
tar tzf $A | grep -c 'paper/SI/fig5_cli_session.txt'    # SI 说收了它，必须是 1
```

> ⚠ 三条都要**看计数**。`tar` 找不到归档时也会打印 0，
> 那是失败不是通过 —— 先确认 `$A` 真的存在。

---

## 五、可移植性（2026-10-23 修）

发布前发现并修掉的两处"只有作者机器能跑"：

1. 14 个测试脚本都写死了 `sys.path.insert(0, "/home/<dev>/...")`
   → 改为 `import _root`，由 `tests/_root.py` 从 `__file__` 推包根。
2. 6 个测试与出图脚本写死了 CATKINAS 的绝对路径
   → 改为 `refdata.catkinas(...)` / `catkinas_mat(...)`，走 `KINSTEP_CATKINAS_DIR`，
   缺失时跳过或退回存档数值。

同时修掉回归脚本里 `if ! cmd; then … exit=$?` 的老毛病——`!` 会把 `$?` 翻成 0，
失败信息里永远印 `exit=0`。
