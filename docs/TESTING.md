# DLLifting 测试说明

## 运行

```bash
cd DLLifting
make run          # 主测试套件 test_dllifting
make run-all      # test_dllifting + test_isgeq
```

编译带 `-DNDEBUG`，避免 `assert` 中断，便于一次跑完所有用例。

## 测试内容（`test_dllifting.cpp`）

| 类别 | 内容 |
|------|------|
| **DP 表** | `Lifting_DPiter` / `threshold=200` 与参考 DP（含二进制拆分 `u`）逐位比较 |
| **DL 表** | `threshold=10` 建表 → 单调性 → `Lifting_Expand` 与参考 DP 比较 |
| **Findsol** | DL 展开后 `Lifting_Findsol` 与参考 DP 查询一致 |
| **Compress/Expand** | DP 表压缩再展开（`isleq=0/1`） |
| **完整升维** | `lifting()`：`thr=10`（DL）vs `thr=200`（DP）系数/rhs 一致 + 枚举可行性 |

参数网格：

- `isleq ∈ {0,1}`（`>=` / `<=` 背包）
- `threshold ∈ {10, 200}`（DL / DP 路径）
- `isdl_mode ∈ {AUTO, DL, DP}`（v1.1.0：`test_mixed_vars` 用 forced DL/DP 对比）
- 物品数、上界 `u`、种子/升维顺序多种手算例 + 随机小例

## 当前结论（请以最新 `./test_dllifting` 输出为准）

### 通过（`isleq=1`，`<=` 背包）

- **DP 与 DL（thr=10）** 在全部手算表用例上 `Expand(DL) == DP == 参考实现`
- **thr=200** 时 DL 路径会切到 DP，与纯 DP 路径一致
- **完整升维** DL/DP 一致且枚举可行域有效

### 失败（`isleq=0`，`>=` 背包）— 代码问题

1. **`Lifting_Mergesort`（`!isleq`）**  
   DL 表 `wsum/psum` 非单调，常出现 `w=0` 重复或 `w=INF` 等异常点。  
   已修 `k=0` 时访问 `newwsum[k-1]` 的 UB，但 **合并逻辑仍有误**。

2. **完整升维 `thr=10` vs `thr=200`**  
   `geq` 实例上 rhs/系数不一致（DL 路径错误）。

3. **个别 `>=` 手算例**（如 `lift_geq_down`）  
   若仍报 `DP lift cut invalid`，需核对种子/容量是否构成合法 `>=` 覆盖初值。

### 建议修复方向

- 重点审查 `DLLifting.cpp` 中 `Lifting_Mergesort` 的 `else` 分支（`!isleq`），与 `Lifting_Expand` 在 `geq` 初值 `dp[j]=INF` 上的交互。
- 修复前：**`>=` 升维请使用 `threshold>100`（走 DP）**，不要用 DL 路径。

## 其它文件

- `test_mixed_vars.cpp`：混合有界/无界变量；分别用 `DLLIFTING_MODE_DL` 与 `_DP` 验证 forced mode
- `test_isgeq.cpp`：`>=` 专项较早用例
- `TEST_CASES_GEQ.md`：手算例参数表
