# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

---

# BioTransition — 交付宪法

`BioTransition` 是一个 R/Bioconductor 包，实现七种 DNB（动态网络生物标志物）方法，检测生物系统中的临界转变（疾病发生、细胞分化、发育转变的 tipping point）。作者本人（Zaoqu Liu）维护，目标是 Bioconductor 提交 + R-universe 分发。

**本文件是最高交付标准。** 任何任务——加方法、改算法、修文档、调 API、动 C++——先过本文件的滤镜，再动手。不达标 = 未完成，不得称 done。

## 0. 第一原则（本项目的取舍优先级）

> **用户体验 > 简洁 > 科学的极致精确。**

这是一个**给科研人员用的工具包**。决定它成败的不是「公式有没有精确到第 15 位小数」——七个方法的公式都来自已发表论文，方法学正确性由论文背书，**既定、务实、差不多即可**。决定成败的是：

- **装得上、跑得通**：依赖能解析，example 能跑，vignette 能编译，双 CI 绿。
- **调用顺畅无断点**：七个方法签名一致，参数命名直觉，默认值合理，报错可操作，进度有反馈。
- **结果可直接用**：返回的 list 结构清晰、命名自解释，下游能直接接绘图 / 富集。

**判断每一个改动，先问的不是「科学上更严谨吗」，而是「用户这一秒调函数会更顺吗」。** 时间成本、认知负荷、参数个数、报错友好度，与功能正确性同等重要。

科学正确性的底线只有一条（见 §6）：**不引入明显错误**。在此之上，把精力投到体验，不投到数值偏执。

## 1. 这是什么 · 绝不能碰

**绝不手改这些自动生成文件**（改了会被下次生成覆盖，且让 CI 漂移）：

| 文件 | 真正的源 | 改法 |
| --- | --- | --- |
| `NAMESPACE` | R 文件里的 `#' @export` / `#' @importFrom` roxygen 标签 | 改标签 → `devtools::document()` |
| `man/*.Rd` | 函数上方的 `#'` roxygen 块 | 改注释 → `devtools::document()` |
| `R/RcppExports.R` + `src/RcppExports.cpp` | `src/fast_correlation.cpp` 的 `// [[Rcpp::export]]` | 改 C++ → `Rcpp::compileAttributes()` |

**NSE 全局变量**：用 `dplyr`/data.frame 列名做非标准求值时，新列名必须加入 `R/zzz.R` 的 `utils::globalVariables(...)`，否则 `R CMD check` 报 "no visible binding" NOTE。

**构建产物**（已 gitignore，绝不提交、绝不依赖）：`BioTransition.Rcheck/`、`BioTransition.BiocCheck/`、`*.tar.gz`、`src/*.o`、`src/*.so`。

**版本与元数据**：改 `DESCRIPTION` 的 `Version` / 发版 / 改 `inst/CITATION` —— 需作者显式确认，不擅自动。

## 2. 开发命令

工具链已就绪：R 4.4.0，`devtools` / `roxygen2` / `testthat` / `Rcpp` / `BiocCheck` 均已装。无 `Makevars`（Rcpp 默认编译即可），无 lintr 配置。

```bash
# 日常开发循环（改 R 代码后）
Rscript -e 'devtools::load_all()'      # 加载（含编译 C++），最快的迭代方式
Rscript -e 'devtools::document()'      # 改了 roxygen / @export / @importFrom 后，重生成 NAMESPACE + man/

# 测试
Rscript -e 'devtools::test()'                      # 跑全部 testthat
Rscript -e 'devtools::test(filter = "tDNB")'       # 只跑 test-tDNB.R（filter 匹配文件名去掉 test- 前缀）
Rscript -e 'devtools::run_examples()'              # 跑所有 @examples（验证文档示例没烂）

# 改了 C++（src/fast_correlation.cpp）—— 必须按序走完，否则 R 看不到新签名
Rscript -e 'Rcpp::compileAttributes()'   # 1. 重生成 RcppExports（R + cpp 两侧）
Rscript -e 'devtools::document()'        # 2. 同步 man（如签名进了文档）
Rscript -e 'devtools::load_all()'        # 3. 重新编译加载

# 发布前门禁（两道都要绿）
Rscript -e 'devtools::check()'           # ≈ R CMD check --as-cran，5 平台 CI 跑的就是它
Rscript -e 'BiocCheck::BiocCheck(".")'   # Bioconductor 规范；CI 只在 error 时 fail，warning 放行

# 手动等价
R CMD build . && R CMD check BioTransition_2.0.0.tar.gz --as-cran
```

**CI**（`.github/workflows/`）：`R-CMD-check.yaml` 跑 5 个平台（macOS / Windows / Ubuntu × release/devel/oldrel）；`BiocCheck.yaml` 只在 **error** 时 fail，warning 允许。本地 `devtools::check()` 绿 ≈ CI 绿。

## 3. 架构：三层计算栈

```
┌─────────────────────────────────────────────────────────────┐
│  7 个导出方法（用户 API）—— 两个家族，共用同一套 CI 脚手架      │
│                                                               │
│  PPI-independent:  cDNB    tDNB                                │
│  PPI-dependent:    LcDNB   LDNB   MDNB   TSNMB   TSLE          │
│                    （吃内置 ppi_h / ppi_m STRING 网络）         │
├─────────────────────────────────────────────────────────────┤
│  共享 R 助手（内部 @noRd；SSPN1/2 已导出）                     │
│   CorandPval         相关 + p 值 + FDR（一步算，Fisher-Z）     │
│   calSFNetforCorMatrix  软阈值 → 邻接 → TOM（tDNB 用）          │
│   detectModules      层次聚类枚举候选模块（栈溢出转可操作错误）│
│   scoreModules       |PCC| 口径 CI 评分（补集恒等式向量化加速）│
│   validateInput      全方法统一输入校验（空态/缺 ref/空 PPI）  │
│   SLE / sNMB         单样本熵 / 网络模块（LDNB/TSLE/TSNMB 用）  │
│   SSPN1 / SSPN2      样本特异扰动网络（LDNB 用；已导出）        │
├─────────────────────────────────────────────────────────────┤
│  C++ 核（src/fast_correlation.cpp，Rcpp，2–20× 加速）          │
│   fast_cor_cpp  fast_cor_pval_cpp  fast_module_score_cpp       │
│   fast_sspn_batch  fast_bh_adjust                             │
└─────────────────────────────────────────────────────────────┘
```

**统一 API 契约（全局一致性的核心，改动必读）**：`cDNB` 与 `tDNB` 几乎逐行同构，七个方法共享同一骨架——

```
(expr 基因×样本, state 两列 df, state.levels) 
  → state 校验为 2 列、factor 成 state.levels、按 colnames(expr) 重排
  → 按 state 切成 ddl 列表
  → 每态: 相关(CorandPval) [tDNB 再算 TOM] → 变异度(sd/cv) 
        → hclust(1 - r) → dendextend::partition_leaves → 按 min/max.size 过滤
        → 每模块 CI = mean(SD_in)·mean(|PCC_in|)/mean(|PCC_out|)  (|PCC| 口径; AddModuleSize 时 ×√size)
  → 取 max-CI 模块 = Candidate；跨态最大 CI = critical state
  → 固定 DNB.genes，回算各态 DNB.score
  → return list(DNB.score, DNB.genes, CI_all, Gene_module, Candidate, ...)
```

**冰山法则**：用户只见 `result$DNB.score` / `result$DNB.genes` 两个字段（10% 简洁）；背后 90% 是相关计算、TOM、层次聚类、模块筛选、CI 评分的完整流水线。返回 list 的每个字段都要命名自解释、能被下游直接消费。

**方法 → 助手依赖图**（改助手前必查，避免改一个崩一片）：

| 方法 | 依赖 |
| --- | --- |
| `cDNB` | `CorandPval` |
| `tDNB` | `CorandPval` → `calSFNetforCorMatrix` |
| `LcDNB` | `CorandPval` |
| `MDNB` | `fast_cor_cpp` |
| `LDNB` | `SLE` → `SSPN1` |
| `TSLE` | `CorandPval` → `SLE` |
| `TSNMB` | `CorandPval` → `sNMB` |

## 4. 用户体验生死线（第一优先级，不可妥协）

工具包的「流式渲染」就是**调用体验**。以下不过关，功能等于不存在：

1. **签名一致性 = 肌肉记忆**。新增 / 修改方法，参数名、顺序、默认值、返回 list 的字段命名，**必须与现有七方法对齐**（`expr` / `state` / `state.levels` / `cor.method` / `p.adjust.method` / `min.size` / `max.size` 这套词汇是契约）。引入第二套命名 = 退回。
2. **报错可操作**。`stop()` 必须告诉用户「哪里错了 + 怎么改」，像现有的 `"Number of state columns must be two, the first is..."`。禁止裸 `stop("error")` / 让底层 `[` 下标越界这种天书报错冒泡到用户。
3. **进度有反馈**。长流程沿用现有 `message("+++ ...")` 风格逐步播报（相关计算 → 聚类 → 评分 → Done）。但**进度反馈走 `message()`，不要混用 `cat()`**（现有代码两种并存，是债，见 §7；新代码统一 `message()`，可被 `suppressMessages()` 静默，测试也这么用）。
4. **默认值要在普通机器上开箱即跑**。任何「典型数据 + 默认参数」必须直接出结果，不需要用户先调参。`detectCores() - 10` 这类默认值是反例（见 §7）。
5. **example / vignette 是门面**。每个导出函数的 `@examples` 必须能真跑（`devtools::run_examples()` 绿）；大计算用 `\donttest{}` 包。文档示例烂掉 = 用户第一印象崩。

## 5. 主动治理：删、简、统一

1. **奥卡姆剃刀**。加任何参数 / 分支 / 返回字段前问：删掉它，用户会变差吗？答否 → 不加。复杂度是负债，不是能力。
2. **最小改动**。只改任务相关代码，不顺手重构、不 reformat 无关文件、不动邻近函数的风格。
3. **统一优于并行两套**。同一件事只有一种做法：一套相关计算入口（`CorandPval` / `fast_cor_*`）、一套 CI 公式、一套进度反馈语义。发现两套并行范式 → 合并或删一套。
4. **C++ / R 双路径要等价**。`SSPN1/2` 有「C++ 批处理 + R parallel」两条路（按 `cor.method == "pearson"` 自动选）。改其中一条，**另一条同步改**，否则用户切 method 就得到不一致结果——这是体验事故。

## 6. 科学正确性：务实底线（不偏执，但不犯错）

- 方法学**以原论文为准**（cDNB→Chen 2012、LDNB→Liu 2019、MDNB→Li 2022、TSNMB→Zhong 2022、TSLE→Liu 2020；tDNB 为本包新方法）。公式既定，**不追求与某参考实现逐位复现**，不为「更严谨」擅自改动已发表的算法。
- 底线只有一条：**不引入明显数值错误**。改了共享计算层（`CorandPval`、CI 公式、C++ 核）后，顺手确认：(a) `devtools::test()` 绿；(b) 受影响的方法 `run_examples()` 不报错、不返回全 `NaN`/`NA`。这是「装得上跑得通」的体验底线，不是科学偏执。
- 已处理的边界（保持，别回退）：零方差样本的相关（`CorandPval` 里 `r.na.value=0` / `p.na.value=1`）、小样本模块检测、`max.size` 超基因数时自动缩到 `0.8 ×`。新增方法要照顾同样的退化输入。

## 7. 已知体验债

**v2.1.0 已清（勿回退）：** nCores 小核机变负数 → `max(1, cores-1)`；全方法进度 `cat()` → `message()`（可被 `suppressMessages()` 静默）；`1:n` → `seq_len/seq_along`；空态 / 缺 `ref` / 空 PPI 统一走 `validateInput` 友好报错；大基因数聚类栈溢出 → 可操作提示；CI 统一 |PCC| 口径；cDNB/tDNB 评分向量化提速 ~19×；删除低效死代码 `fast_module_score_cpp`；TOM 重写为 Zhang & Horvath 标准 + 分母守护 / [0,1] clamp / 单位自重叠；SSPN1/2 数字 `ref.samples` 的 C++ segfault 修复；SSPN 的 C++ Z 公式去除多余 `sqrt`、对齐 R + SSN 原论文；`fast_cor_pval_cpp` 改用 t 分布（与 `CorandPval` 一致）。

**仍待办：**
- **内置 PPI 数据 ~101Mb**（`data/`）触发 R CMD check 的 size NOTE，Bioconductor 可能介意；长期可考虑挪到 ExperimentHub。这是目前唯一已知的非阻塞项。

已核实非 bug：**LDNB 的 `rbind`** 把 (A,B) 与翻转的 (B,A) 合并，是把 `SSPN1` 输出的单向无向边（经 `sort + distinct` 去重）对称化，供下游 `split(V2, V1)` 建对称邻接表——必要且正确，不是双计数。

无编号 TODO 不留在代码里：今天删、今天做，或开 issue 编号。

## 8. 交付前二进制清单（任一「否」→ 不得交付）

**A · 装得上跑得通（体验底线）**
- [ ] `devtools::check()` 绿（5 平台 CI 跑的就是它）
- [ ] `BiocCheck::BiocCheck(".")` 无 error（warning 可议）
- [ ] `devtools::test()` 全绿；`devtools::run_examples()` 不报错

**B · 调用体验**
- [ ] 新 / 改方法签名、参数名、返回字段命名与现有七方法一致
- [ ] 报错可操作（告诉用户怎么改，非天书）；长流程有 `message()` 进度
- [ ] 默认参数在普通机器（含小核数）开箱即跑

**C · 一致性与生成物**
- [ ] 改了 C++ → 跑过 `Rcpp::compileAttributes()`，`RcppExports` 两侧已同步
- [ ] 改了 `@export`/roxygen → 跑过 `devtools::document()`，`NAMESPACE` + `man/` 已更新
- [ ] 没手改任何自动生成文件；没提交构建产物；新 NSE 列名进了 `zzz.R`
- [ ] C++/R 双路径（SSPN）若动其一，另一条已同步

**D · 复杂度**
- [ ] 无可删而未删的参数 / 分支 / 字段；无新引入的第二套范式
- [ ] 改动只追溯到：用户需求 / 失败测试 / 已证 bug —— 无顺手改的无关代码

**判定：A–D 全勾 = pass；任一否 = block。** 不靠主观感觉，靠以上可验证项。

---

# 工程交付强制规范（极致严谨·零妥协）

## 0. 总纲

本规范为交付硬约束，不是建议。任一条款未达标，视为未完成，不得交付。

默认立场：能删则删，能并则并，能简则简。复杂度不是能力，是负债。

## 1. 需求与决策

1. 严格遵循奥卡姆剃刀：每一个功能、抽象、组件、状态、动画，都必须回答「删掉它，用户会变差吗？」若答案是否，必须删除。
2. 全程站在用户视角推导需求，禁止工程师自嗨式实现。判断标准不是「能不能做」，而是「用户在这一秒是否需要、是否感知价值、是否愿意为此等待」。
3. 时间成本、认知负荷、操作步数、视觉干扰，与功能正确性同等重要，必须一并纳入设计，不得事后补丁。

## 2. 视觉与交互

1. 视觉遵循极简设计美学：克制、纯粹、留白、节奏、层次。拒绝装饰性复杂度，拒绝「看起来做了很多」的空洞感。
2. 简约不等于简陋。每一屏、每一帧、每一次过渡，都必须经得起逐帧审视——每一帧都应值得截图。
3. 全流程交互流转无断点：输入、反馈、加载、错误、空态、恢复，全部连贯一致，不得出现体验断层。
4. 整体操作体验必须绝对丝滑。凡用户可感知的动作——滚动、切换、拖拽、聚焦、展开——须如物理世界般自然，零迟滞、零突兀、零「等一下」感。

## 3. 性能（底层硬约束）

1. 全链路以高性能为底层硬性约束。逻辑、渲染、交互、网络、内存，全部优先保障性能表现；性能不是优化项，是设计项。
2. 禁止主线程长时间阻塞；禁止无必要的重渲染、布局抖动（layout thrashing）、级联 reflow；能异步则异步，能延迟则延迟，能虚拟化则虚拟化，能缓存则缓存。
3. 长列表、长会话、大文本、高频更新场景，必须在设计阶段给出性能方案，禁止靠「机器够快」硬扛。
4. 交付时须能说明：关键路径耗时、渲染频率、内存边界、极端数据量下的退化策略。

## 4. 流式渲染（产品化红线·不可妥协）

流式渲染是产品核心硬指标，是商业化成熟产品的核心评判标准。

1. 输出渲染全程须保持 60fps 级流畅；禁止可见掉帧、跳帧、闪烁、布局重排导致的画面抖动。
2. 流式输出过程中，用户向上翻阅历史内容时，必须零卡顿、零延迟、零拖影、零 scroll jump、零内容位移。
3. 禁止因新内容插入导致阅读位置被抢、视口被拽、滚动条异常跳动。须正确处理 scroll anchoring、虚拟列表、增量布局稳定性。
4. 流式更新不得阻塞用户交互；用户滚动、选中、复制、暂停时，系统必须优先响应用户意图，而非无脑追写最新 token。
5. 若做不到上述任一条，不算完成，必须重构，不得用「差不多」「大多数情况可以」搪塞。

## 5. 主动治理（删、改、重构）

1. 主动识别所有逻辑不合理、架构冗余、体验割裂、视觉违和、命名混乱、状态泄漏、边界缺失的模块。
2. 对以上问题执行删除、优化或整体重构；不做折中妥协，不留「以后再说」，不接受临时 hack 长期驻留。
3. 若局部优化破坏全局一致性，宁可回退局部改动，也不允许引入系统级割裂。

## 6. 全局统一

1. 针对任何细节，若产出更优实现方案或底层设计哲学，优化后必须保证架构、视觉、交互、命名、错误处理、动效节奏全局统一。
2. 禁止同一产品内出现多套平行范式：两套按钮逻辑、两套 loading 语义、两套滚动行为、两套状态管理模式。
3. 所有新增须能归入现有设计系统 / 架构分层；若不能归入，先修系统，再加功能。

## 7. 代码标准

1. 结构整洁干净：命名即文档，结构即意图；读代码应像读 prose，不应像解谜。
2. 底层执行高效：热点路径优先；禁止过早抽象，也禁止复制粘贴式重复。
3. 思考覆盖全场景：正常路径、异常路径、弱网、断连、重试、取消、并发、竞态、长会话、大数据、空数据、权限变化，须在实现前纳入设计。
4. 周全处理用户维度：时间计算、视觉分层、资源调度、焦点管理、可访问性、键盘操作、国际化扩展性，全部纳入交付范围，不得遗漏。
5. 禁止交付：魔法数、隐式副作用、无法追踪的状态源、无注释却非自明的 trick、仅为「能跑」的脆弱实现。

## 8. 冰山法则（全局架构思维）

1. 编码必须具备全局架构思维：改一处，须评估十处；模块边界、数据流、生命周期、依赖方向必须清晰。
2. 严格遵循冰山法则：表层简洁易读，底层逻辑完整厚重。
3. 用户只见 10% 的交互简洁；背后 90% 的健壮性、性能、容错、可观测性，必须默默托底，且可被维护者理解。

## 9. 交付前强制验收（深度复盘 + 二次全量核验）

交付前必须完成逐粒度排查，杜绝疏漏。以下任一项为否，禁止交付：

1. 是否存在任何可感知卡顿，尤其是流式输出时的上下滚动？
2. 是否存在 layout shift、scroll jump、内容抢焦点？
3. 是否存在可删而未删的复杂度？
4. 是否存在与全局设计语言冲突的实现？
5. 是否存在仅 happy path 可用、边界即崩的实现？
6. 是否能在不看作者解释的情况下，被他人快速理解与修改？
7. 若交给顶级工程团队 review，是否会感到羞愧？

## 10. 终极标准

交付代码须达到行业标杆范本层级：

1. 可供 OpenAI、Google 等头部企业逐帧拆解、完整学习、规范复刻。
2. 整体品质对标顶级大厂工程标准：美学、性能、架构、可维护性、用户体验，五者同时达标，不允许单项牺牲换另一项。
3. 每一处取舍都必须经得起三连问：为什么这样写？为什么不再简？为什么这样快？

对得起这个标准。未达标，不算交付完成。
