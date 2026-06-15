# DLLifting：`isleq = 0`（`>=` 背包）测试算例说明

## 算法与数学背景（简要）

### 问题设定

- **`isleq = 1`**：背包约束为 \(\sum_i w_i x_i \le b\)，升维子问题在剩余容量上求 **最大** 利润（DP：`dp[j] = max`）。
- **`isleq = 0`**：背包约束为 \(\sum_i w_i x_i \ge b\)，与上式对偶/对称；子问题在“达到至少某重量”的意义下求 **最小** 利润（DP：`dp[0]=0`，`dp[j>0]=+\infty` 初始化，转移取 `min`）。

代码入口：`lifting(..., isLeq, ...)`，`DLLifting.h` 中 `isleq` 注释。

### 升维（Lifting）在做什么

在覆盖不等式（种子不等式）\(\sum_i p_i x_i \le \text{rhs}\) 上，按 `liftingorder` 依次对变量做 **向上升维**（`isuseub=0`，`Lifting_Up`）或 **向下升维**（`isuseub=1`，`Lifting_Down`），得到仍对背包可行域有效的更强不等式。

对 **`>=` 背包**，`Lifting_Up` / `Lifting_Down` 与 **`<=`** 对称：

| 操作 | `isleq=1` (`<=`) | `isleq=0` (`>=`) |
|------|------------------|------------------|
| Up   | \(\alpha \leftarrow \min\) | \(\alpha \leftarrow \max\) |
| Down | \(\alpha \leftarrow \max\) | \(\alpha \leftarrow \min\) |

见 `Lifting_Up` / `Lifting_Down` 中 `if(lift->isleq)` 分支。

### 数据结构

1. **DP 表** `dplist[j]`：`isleq=0` 时为达到重量 \(j\) 的最小 \(\sum p\)（`Lifting_DPiter` 中 `min` 转移）。
2. **DL 表** `(wsum[k], psum[k])`：非支配点列；`Lifting_Findind(..., isleq=0)` 返回 **最小的** \(k\) 使 `wsum[k] >= cap`（与 `<=` 时 “最大 `<= cap`” 相反）。

`threshold` 控制 DL / DP 切换：`<100` 偏 DL，`>100` 偏 DP。

### 正确性检验（测试用）

对 `>=` 实例，枚举 \(x_i \in \{0,\ldots,u_i\}\)，对所有满足 \(\sum w_i x_i \ge b\) 的 \(x\) 检查：

\[
\sum_i p_i x_i \le \text{rhs}.
\]

---

## 如何运行测试

```bash
cd DLLifting
make run          # 编译并运行 test_isgeq
./test_isgeq      # 仅运行
```

---

## 手算/小规模算例（推荐优先调试）

### 算例 G1：单物品，仅种子（DP 路径）

| 字段 | 值 |
|------|-----|
| 约束 | \(5 x_0 \ge 5\) |
| `n` | 1 |
| `w` | `[5]` |
| `u` | `[1]` |
| `p` 初值 | `[1]`（种子） |
| `cap` | `5` |
| `seed` | `[0]` |
| `liftingorder` | `[]` |
| `isleq` | `0` |
| `threshold` | `200`（走 DP） |

期望：升维后仍满足可行性；`rhs` 由种子表在 `subcap=5` 处查询得到。

### 算例 G2：三物品，先 Up 再 Up（`test_isgeq` 中 `geq_tiny_3var`）

| 变量 | `w` | `u` | 初值 `p` | 说明 |
|------|-----|-----|----------|------|
| 0 | 5 | 1 | 1 | 种子 |
| 1 | 3 | 1 | 0 | Up 升维 |
| 2 | 2 | 1 | 0 | Up 升维 |

- 背包：\(5x_0 + 3x_1 + 2x_2 \ge 5\)
- `seed = [0]`，`liftingorder = [1, 2]`，`isuseub = [0,0,0]`

### 算例 G3：含向下升维（固定于上界）

| 变量 | `w` | `u` | `isuseub` |
|------|-----|-----|-----------|
| 0 | 4 | 1 | 0 |
| 1 | 3 | 2 | **1**（Down） |
| 2 | 5 | 1 | 0 |

- 约束：\(4x_0 + 3x_1 + 5x_2 \ge 4\)
- `seed = [0]`，`liftingorder = [1, 2]`

### 算例 G4：大上界（`Mergesortinf` / 多副本）

| 字段 | 值 |
|------|-----|
| `w` | `[3, 7]` |
| `u` | `[1, 8]` |
| 约束 | \(3x_0 + 7x_1 \ge 10\) |
| `seed` | `[0]`，`liftingorder` | `[1]` |

用于触发 `!isleq` 下 `w*u >= maxcap` 时的 `Lifting_Mergesortinf` 分支。

### 算例 G5：`Findind` 回归（不调用完整 lifting）

手工 DL 表 `wsum = {0,3,7,12,20}`，`cap=7`，`isleq=0` → 应返回索引 `2`（重量 7）。

---

## 自动化测试项（`test_isgeq.cpp`）

| 测试 | 目的 |
|------|------|
| `test_findind_geq` | `Lifting_Findind` 在 `isleq=0` 的二分语义 |
| `test_dp_geq_*` | `Lifting_DPiter` 与参考最小化 DP 一致 |
| `test_up_lifting_step_geq` | 单步 `Lifting_Up` 与暴力公式一致 |
| `test_dp_mode_switch` | `threshold>100` 时 `isDL=false` |
| `test_full_lifting_geq_*` | 完整 `lifting(..., isleq=0, threshold=200)` + 全枚举有效性 |
| `test_dl_dp_paths_agree` | 同一实例 DL(`thr=10`) 与 DP(`thr=200`) 结果一致 |
| `test_random_geq_validity` | 随机小微实例有效性 |

---

## 已知问题（测试中发现）

1. **`Lifting_Mergesort` 的 `!isleq` 分支**曾被人为改写成简化版，导致 DL 表 `wsum` 乱序；已按 `HiGHS/highs/mip/DLLifting.cpp` 恢复旧版合并逻辑。若 DL/DP 仍不一致，请对比 HiGHS 版本继续排查。
2. 仓库中 **`RandomTest` / `Integer/main.cpp` 默认 `isleq=1`**，不会覆盖 `>=` 路径；测 `>=` 必须显式传 `isleq=0`。
3. **`MinimalCover`** 针对 **`<=` 覆盖** 设计；`>=` 实例应 **手工指定** `seed` / `liftingorder`（本文算例即如此）。

---

## 在业务代码中调用示例

```cpp
int isleq = 0;   // >= 背包
double rhs;
lifting(&lift, p, w, u, isuseub, cap, 0,
        seed, n_seed, liftingorder, n_liftingorder,
        &rhs, isleq, x, n, threshold, duration);
```

`cap` 为背包右端 \(b\)；`issubcap=0` 表示传入的是原约束容量。
