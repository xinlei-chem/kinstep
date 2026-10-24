# KinSTEP（中文版）

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23613003.svg)](https://doi.org/10.5281/zenodo.23613003)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)

> 英文版 [`README.md`](README.md) 是发布正本（CPC 的读者看英文）；本文与它内容一一对应。

**Kin**etics + **S**ono / **T**hermal / **E**lectro / **P**hoto —— 一个内核，四种催化模态。

纯 Python 微动力学（microkinetic modelling）：从基元反应式出发，算覆盖度、速率、
速率控制度（DRC）、反应级数、表观活化能、Tafel 斜率、火山图。

```bash
git clone https://github.com/xinlei-chem/kinstep.git
cd kinstep

pip install -e .            # 核心只要 numpy + scipy
pip install -e ".[all]"     # 外加 mpmath（高精度）与 matplotlib（作图）

kinstep templates           # 看有哪些现成体系
kinstep new --template her  # 生成骨架，机理已填好
kinstep check her.kinstep.toml   # 它会逐条告诉你还缺什么
kinstep run   her.kinstep.toml --drc
```

---

## 一、它想解决什么

现有的微动力学工具有两个共同的麻烦：

1. **要填复杂的输入脚本。** CATKINAS 的输入是一段 MATLAB 代码，物种名就是变量名
   （物种叫 `T` 会把温度覆盖掉）。KinSTEP 让你**照着纸上写反应式**：

   ```
   N2 + 2* <-> 2N*
   NH2* + H* <-> NH3 + 2*
   H2O(l) + * <-> OH* + H+ + e-    {beta_ct: 0.5}
   ```
   `*` 是位点，`<->` 可逆，`->` 不可逆，`(l)` 液相，`e-` 电子。
   **电子计量系数 `nu_elec` 自动从方程里数出来** —— 它是最容易写错符号的参数。

2. **只覆盖一种模态。** 热催化的工具不懂电化学，电化学的不懂光催化。
   KinSTEP 把"外场如何改变基元步"抽成 **driver 插件**，内核不必知道这是哪种催化：

   | driver | 改什么 | 典型可观测量 |
   |---|---|---|
   | `thermal` | 什么都不改（基线） | — |
   | `electro` | `ΔG(U) = ΔG(0) − ν_e·e·U`（CHE/BV） | 电流密度、Tafel 斜率 |
   | `photo` | `ΔG = ΔG(暗) − ν·ΔE_qfl`（准费米能级分裂） | AQY、载流子浓度 |
   | `sono` | **不改能量，改工况本身**（空化热点多工况加权） | — |
   | `field` | 偶极 + 极化率的场效应 | — |

   因为 `−ν·ΔE_qfl` 与 `−ν·e·U` **数学同构**，`photo + electro` 直接叠加就是光电催化（PEC），
   不需要第三套代码。

---

## 二、三条贯穿始终的设计原则

### ① 能量填错不会报错，只会给出一条看着很合理的曲线

所以每个数值都带**来源标签**（`provenance.py`）：

```python
Val(1.21, "dft_file", detail="./OUTCAR")      # 可信
Val(0.80, "llm_suggested")                     # ★ 能量类 ⇒ 求解器直接拒绝
Val(0.50, "default")                           # 能跑，但**默认值必须可见**
```

九级来源。**能量类的 `llm_suggested` / `guess` / `REQUIRED` ⇒ 拒绝求解**，并逐条说明怎么修。

### ② 判据要落在"能不能真的用"上

本项目的每一道闸都对着一类**不会报错、只会把结论带偏**的错误：

| 闸 | 抓什么 | 实际抓到过 |
|---|---|---|
| 元素/电荷配平 | 笔误 | 模板里 `H2O* -> OH* + H*` 位点 1≠2 |
| 位点守恒 | 脱附忘了写回空位 | `NH2* + H* <-> NH3 + *` |
| 撞保留名 | 物种名与物理量撞车 | `T`、`theta`、`E` |
| **Wegscheider 循环** | 某一步能量填错一个数 | 改 0.1 eV 当场被抓 |
| **Campbell 求和 ΣX_B=1** | DRC 算错 | 也用来**自动选差分步长** |
| **通量守恒 νᵀr=0** | 准平衡步的速率不自洽 | 1.55e-07 → 2.22e-16 |
| **电荷守恒（光催化）** | `nu_elec_cb`/`nu_hole_vb` 配错 | 电子/空穴消耗差 90 倍 |
| **权重归一 Σw=1（超声）** | 多工况加权没意义 | 权重和 1.4 被拦住 |

### ③ 精度边界要测出来，不能假装没有

速率律里 `1−z`（z 是可逆度）在准平衡步上会被相消吃掉。本项目实测清楚了：

- 瓶颈**不在** `expm1`（它已经全精度），而在 `ln z = ΔG/k_BT + Σν·ln a` 的**绝对精度**
- double 下 `1−z ≈ 1e-30` 时一位有效数字都不剩
- ⇒ 求解器是**四级阶梯**，每一级都是被实测逼出来的：

```
① ODE 积分（自适应续算 + 墙钟预算）    → 物理初值
② 阻尼牛顿（解析雅可比 + 限幅 + 回溯）  → double 解
③ 通量投影（σ = R_f·eps 加权）        → 准平衡步的速率由守恒推定
④ 高精度牛顿（mpmath，自动触发）       → 把 1−z 真正算出来
```

`precision="auto"` 会自己判断要不要升级，并把原因写进 `notes`。

---

## 三、验证到什么程度

| 对照 | 量 | 结果 |
|---|---|---|
| **Koper 1998** (JCP 109, 6051) CO stripping | ODE 右端逐点 | **2.4e-15** |
| 同上 | 电流峰值 (mA/cm²) | 与参照实现**六位全同** |
| 解析闭式解（甲酸分解） | θ 与 r | **2.2e-16** |
| **CATKINAS 401 点**（合成氨描述符扫描） | 覆盖度 / 可逆度 | **1.4e-10** |
| 同上 | 净速率（高精度） | **2.7e-10**，误差 >1e-6 的点 **27 → 0** |
| 同上 | DRC / Campbell 求和 | **3.8e-12 / 2.9e-12** |
| **CATKINAS 315 点**（T–P 二维 + 流动反应器） | 出口分压中位 | **2.2e-16**，315/315 收敛 |
| NIST Shomate | 熵校正 5 个常数 | **1.9e-10 eV** |
| 光催化文献 | 光强幂律（**算出来的**） | 0.996 → **0.500** |

> 这些对照**不是"看着差不多"**，是逐数字比到输入常数的精度地板。
> 过程中还发现与 CATKINAS 有**五处常数口径差异**（k_B、e、h、F、k_flow 里的 T），
> 每一处都定位到具体数值并显式处理 —— 不消掉它们，比的就是常数而不是公式。

全量回归：`./run_all_tests.sh`（14 套数值 + 关键数字退化检查 + CLI + 错误路径）。

---

## 四、五个命令

```bash
kinstep templates [名字]      # 11 个现成体系（热/电/光/超声），只给机理不给能量
kinstep new --template her    # 生成骨架
kinstep ask "IrO2 的析氧，1.5 V vs SHE，25 摄氏度" --write
                              # 一句话 → 结构 + 问题（零 LLM）
kinstep check  xxx.toml       # 四道闸 + 来源闸 + driver，每条报错都给"怎么改"
kinstep explain xxx.toml      # 用人话说清在算什么、哪些值是默认的
kinstep run    xxx.toml --drc # 求解 → 报告 → 出图 → 导出 .dat
kinstep drivers               # 模态、参数、"这个参数是谁要的"
```

退出码：`0` 正常 / `2` 输入有问题 / `3` 来源闸拦住 / `4` driver 配置错。

### 自然语言层的红线

`kinstep ask` 把一句话变成**结构 + 问题**。数字分两类，处理完全不同：

- **你自己说出来的**（「673 K」「pH=0」「−100 mV」）⇒ 抽取，标 `user_given`，
  **单位换算显式记录**（「400 摄氏度」→ `673.15 K（400 °C + 273.15）`）
- **系统要去猜的**（任何能量）⇒ **一律变成问题**

而「100 mW/cm²」→ 光子流**需要光子能量**，没给就是**一个问题，不是一个假设**。

LLM 后端**默认关**，装它要显式 `confirm=True`，且产出要过三道闸：
白名单 → 机理解析 → 来源闸。实测：一个返回 `Ga_fwd=0.42` 的恶意后端，
保留下来的只剩 `template_key` / `questions` / `modality`。

---

## 五、文档

- `EXAMPLES.md` —— 六个可直接跑的案例（热 / 电 / 光 / 超声 / 扫描 / 流动反应器）
- `notes/DEVLOG.md` —— 2000+ 行开发台账：每个设计决定的依据、每个 bug 的现象→根因→修复。
  **仅在开发树里，不随存档发布**（内部记录、中文、体量大），理由见 `RELEASE_NOTES.md` 第四节
- `ATTRIBUTION.md` —— 借鉴了哪些已发表方法，以及避让条款

## 六、引用

用到的方法请引原始文献，不要引本包：

- 可逆度迭代 + 因式形式速率：Chen, Wang, Hu, *ACS Catal.* **2016**, 6, 7078
- CHE（计算氢电极）：Nørskov et al., *J. Phys. Chem. B* **2004**, 108, 17886
- 速率控制度 DRC：Campbell；Stegelmann, Andreasen, Campbell, *JACS* **2009**
- 气相热化学：NIST Chemistry WebBook (SRD 69)；Chase, *NIST-JANAF*, 4th ed., 1998
- 电化学回归真值：Koper et al., *J. Chem. Phys.* **1998**, 109, 6051
- 热催化回归真值：CATKINAS, *J. Comput. Chem.* **2021**, 42, 379

## 七、许可与发布

**GPL-3.0-or-later**，全文见 `LICENSE`。

第三方材料的发布口径见 **`RELEASE_NOTES.md`**，要点：

- CATKINAS 的算例输出 `.mat` 是闭源发行物自带数据、**未授再分发**，所以不随本包提供；
  但**本套比对用到的那些参考数组**（`tests/catkinas_refdata.npz`，0.2 MB）随包分发 ——
  登载比对值与再分发软件是两回事，这一步让**全部断言与全部插图在干净下载后都能跑**。
  备有原件（设 `KINSTEP_CATKINAS_DIR`）时测试改读原件，并有一项断言逐位核对两者一致。
  旧说明保留供对照：数据缺失时按退出码 77 跳过，并打印需要哪个文件、
  怎么设 `KINSTEP_CATKINAS_DIR`。
- CatMAP 上游克隆是 GPLv3、法律上可带，但本套回归的真值**一个字节都不来自它**，
  所以也不带（避免被误读成衍生作品）。
- 回归共 **82 项，全部能在干净下载后跑**。实测：装了可选件 **82 过 / 0 失败**（有无原件都一样）；
  只有 numpy+scipy 时 **77 过 / 5 跳过** —— 5 项跳过里 4 项需要 mpmath（要比的数在双精度地板以下），
  1 项是出图需要 matplotlib。
- **界面双语**：`KINSTEP_LANG=en` 切英文输出（默认中文）；T16 回归逐命令扫描，
  英文模式下输出里**不得出现一个汉字**。
