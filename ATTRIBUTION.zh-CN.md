# 来源与引用（ATTRIBUTION）

> 本文件记录 KinSTEP 借鉴了什么、必须引用谁。**任何新增借鉴都要先更新本文件再写代码。**
> 建于 2026-10-22。

## 〇、总原则

1. **零代码复制。** 截至 2026-10-22，`kinstep/` 下全部代码为自行编写，
   未从 CatMAP / PyCatKin / OpenMKM / CATKINAS / Zacros 任何一方复制代码。
   将来若要复用任何代码，**先核其许可证并在此登记**。
2. **方法必须引用。** 用了别人发表的方法，就在 README、论文方法段、
   以及**实现该方法的源文件顶部注释**三处同时引用。
3. **数据取原始出处。** 不从二手手册抄数据表。

## 一、★ 必须引用的方法

### 1.1 可逆度迭代 + 因式形式速率（求解内核）

```
r_i = kf_i · ∏c · (1 − z_i),    z_i ≡ r_i(逆)/r_i(正)
```

**出处**：Chen, J.; Mao, Y.; Wang, H.; Hu, P. "Reversibility Iteration Method for Understanding
Reaction Networks and for Solving Microkinetics in Heterogeneous Catalysis."
*ACS Catal.* **2016**, *6*(10), 7078–7087. **doi:10.1021/acscatal.6b02405**
（CATKINAS 的求解核心）。

> 2026-10-23 更正：先前本条写作「Chen, J.; Wang, H.; Hu, P. 等」，**漏了第二作者 Mao, Y.**，
> 且标题是凭印象转述的。已按 Crossref 逐字段核对改正。论文 `paper/manuscript.md` 的 [16] 同步改正。
> 教训：引文一律过 Crossref API 核，不靠记忆写。

**我是怎么知道的**：读 `工具/小程序1/Examples/result_single/log1` 时，
看到 CATKINAS 自动推出的速率式写成 `r(1) = −P_H2·Q_v²·kf(1)·(z1−1)` 而非
`kf·A − kr·B`，由此理解到因式形式是避免灾难性相消的关键。

**处理**：在 `kinstep/solve.py` 顶部注释、README、以及将来任何论文的方法段引用上述文献。
**绝不可把它写成本项目的发明。**

### 1.2 求解阶梯（ODE 拿初值 → 牛顿精收敛）

通用数值实践，但**具体可行的顺序**是从 CATKINAS 的 log 里读出来的：

```
直接高精度牛顿 ❌ → 纯可逆度牛顿 ❌ → ODE15S 拿初值 + 可逆度牛顿 ✅（9 次迭代）
```

**处理**：源文件注明 "solver ladder inspired by CATKINAS (Chen & Wang, ECUST)"。

### 1.3 DRC（速率控制度）

Campbell, C. T. 的 degree of rate control。公开方法，正常引用即可。
**自检**：DRC 求和 ≈ 1（Campbell 求和规则）——这条也来自原文。

### 1.4 CHE（计算氢电极）

Nørskov, J. K. 等, *J. Phys. Chem. B* **2004**, *108*, 17886。电催化标准方法。

### 1.5 电化学内核的校验基准

Koper, M. T. M.; Jansen, A. P. J.; van Santen, R. A.; Lukkien, J. J.; Hilbers, P. A. J.
*J. Chem. Phys.* **1998**, *109*, 6051–6062. DOI: 10.1063/1.477230
（Pt 上 CO 氧化的平均场模型 + 格子 MC）

以及其 Python 实现：Joel W. Ager (LBL), `EC CO stripping MKM.ipynb`。
**本项目用它做回归测试的真值**，不复制其代码。

## 二、★ 数据来源

| 数据 | ❌ 不要从哪拿 | ✅ 从哪拿 |
|---|---|---|
| Shomate 7 参数（A–H）、热容/熵 | **不从 CATKINAS 手册附录 7.2 抄表** | **NIST Chemistry WebBook** 原始条目，并逐条标注 |
| 振动频率（算 ZPE） | 同上 | NIST / 原始文献 |
| 元素相对分子质量 | — | IUPAC |

> 2026-10-22 自我更正：先前曾说"手册附录 7.2 可直接整理成 JSON 用"——**该说法作废**。
> 那些数据的原始出处是 NIST，从二手手册抄既有版权风险又容易带进转录错误。

## 三、与同类工具的功能重合（诚实版）

### 3.1 重合的部分（**刻意为之**，为了能用它们的算例做回归验证）

| 功能 | 已有于 |
|---|---|
| 描述符扫描 + scaling 关系 + 1D/2D 火山图 | **CatMAP**（Medford et al. 2015，SUNCAT） |
| DRC、反应级数、表观活化能 | CatMAP / MKMCXX / CATKINAS |
| 平均场稳态 + CSTR | PyCatKin / OpenMKM / CATKINAS |
| 电化学（CHE + 势依赖） | CatMAP 的 `hbond_electrochemical` 模式 + ORR 教程 |

**⚠ 结论**：若 KinSTEP 只做到"内核 + 扫描 + DRC"，
**基本等于 CatMAP 的 Python 重写**，没有发表价值。
重合是**手段**（为了有真值可验证），不是目的。

### 3.2 查不到任何已有实现的部分（立项依据）

| 特性 | 现状 |
|---|---|
| **多模态 driver 插件层**，热/电/光/超声/场**可叠加** | CatMAP 把电化学做成"校正模式"，不是可插拔层；其余工具无 |
| **光催化 / 超声催化 MKM** | **专门检索，未发现任何工具** |
| **数据来源追溯 + 未核实能量拒绝求解** | 所有工具都默认用户填的数是对的 |
| **免输入脚本**（交互向导 + `need()`） | CatMAP 要 `.mkm` + py 脚本；CATKINAS 要查魔数表 |
| 撞名在注册时拦截 / 多工况加权 | 无 |

## 四、避免混淆的具体约定

| 事项 | 约定 | 原因 |
|---|---|---|
| 输入文件扩展名 | **`.kinstep.yaml`**，**不用 `.mkm`** | `.mkm` 是 CatMAP 的输入后缀，撞了会造成工具链混淆 |
| 包名 | `kinstep` | `mkm` / `pymkm` 在 PyPI 已被占 |
| 命令名 | `kinstep` | 无冲突 |
| 文档措辞 | 不写"first"/"novel"，写"据我们检索未见" | 检索不等于穷尽 |

## 五、发布口径（2026-10-23 定）

详见 **`RELEASE_NOTES.md`**。要点：

| 第三方材料 | 许可证 | 随代码存档分发？ | 理由 |
|---|---|---|---|
| CATKINAS 算例输出 `.mat` | 闭源，手册 §2.2 只授「使用须引用」 | **否** | **未授再分发许可** |
| CatMAP 上游完整克隆 | GPLv3（允许再分发） | **否** | 对回归无用（真值零字节来自它），带着反而易被误读为衍生 |
| NIST Shomate 系数 | 公开数据，逐条标注来源 | 是（在代码里） | 取自 WebBook 原始条目，非二手转抄 |
| Koper 1998 模型 | 文献方法 | 不适用 | T5 按论文重新实现，不读外部文件 |

KinSTEP 本体以 **GPL-3.0-or-later** 发布（`LICENSE` = GNU 官方原文，md5 `d32239bc…`）。
`MANIFEST.in` 里 `prune 真值` 把上表中不分发的部分排除在 sdist 之外。

存档分发的是**本套比对用到的参考数组**（`tests/catkinas_refdata.npz`，0.2 MB），
不含软件本体与其算例文件 ⇒ **82 项断言全部能在干净下载后跑**。
把新打的 sdist 装进空 venv 实测：只有 numpy+scipy 时 **77 过 / 0 失败 / 5 跳过**
（5 项跳过都是缺可选件 mpmath / matplotlib）；装齐可选件 **82 过 / 0 失败**。
三种配置的实测结果见 `RELEASE_NOTES.md` 第三节。

## 六、待办

- [x] 真要复用任何第三方代码前：逐条核许可证并登记于此 —— 已逐条核完，见第五节
- [x] T11 做热力学库时：从 NIST 取数并逐条标注来源
- [ ] 发表前：请 CATKINAS 作者确认 1.1/1.2 的引用措辞是否得体
