# 精简构建说明（node11 / zqw）

本文件记录本次「精简编译」的全部改动与复现方法。
日期：2026-10-11

## 1. 构建结果

| 项目 | 路径 / 值 |
|---|---|
| 源码 | `$HOME/crest-modify` |
| 安装（二进制） | `$HOME/crest-modify-bin/bin/crest` |
| 静态 OpenBLAS | `$HOME/opt/openblas-static`（89 MB） |
| 构建工具（cmake/ninja） | `$HOME/opt/buildtools`（pip 安装） |
| 版本 | CREST 3.0.2（locally modified source） |
| 清洁构建耗时 | **约 8 分钟**（48 核，515 个编译步骤） |

特性开关（已写入 `_build/crest_metadata.fh` 验证）：

```
WITH_TOMLF    = ON      # TOML 输入解析，必需
WITH_GFN0     = ON      # 构象搜索开头的灵活性/拓扑测量用内部 GFN0，必需
WITH_TBLITE   = OFF     # 内部 GFN2-xTB（~6200 文件，编译变慢的主因）
WITH_GFNFF    = OFF     # 内部 GFN-FF
WITH_LIBPVOL  = OFF     # 外压计算
WITH_LWONIOM  = OFF     # ONIOM
```

`ldd crest` 结果：**无缺失库，也没有动态 BLAS/LAPACK**（OpenBLAS 已静态链入）。

### 为什么必须静态链 OpenBLAS

本机（Rocky 9）**没有安装任何 BLAS/LAPACK**，也没有 root 权限。所用方案：
从 Rocky CRB 仓库取 `openblas-static` RPM（0.3.29），`rpm2cpio` 解包后把
`libopenblas.a`（DYNAMIC_ARCH 构建，含 3622 个内核对象，可在任意 x86-64 CPU 运行）
放到 `$HOME/opt/openblas-static/lib/`。

> `cmake`/`ninja` 也不在系统里，用 `python3 -m venv` + `pip install cmake ninja`
> 装在 `~/opt/buildtools`，不需要 root。

## 2. 如何重新编译

```bash
export PATH=$HOME/opt/buildtools/bin:$PATH
cd $HOME/crest-modify
NO_SUBMODULE=1 JOBS=48 STATIC_BLAS=1 \
  OPENBLAS_DIR=$HOME/opt/openblas-static \
  bash build_crest.sh $HOME/crest-modify $HOME/crest-modify-bin
```

只想快速迭代（不重新 configure，只编改动过的文件）：

```bash
export PATH=$HOME/opt/buildtools/bin:$PATH
ninja -C $HOME/crest-modify/_build -j48
ninja -C $HOME/crest-modify/_build install
```

### 如果以后又需要「进程内 tblite」（即 `-gfn2`）

源码已在服务器上，不用重新上传，只要加一个开关重编：

```bash
export PATH=$HOME/opt/buildtools/bin:$PATH
cd $HOME/crest-modify
NO_SUBMODULE=1 JOBS=48 STATIC_BLAS=1 WITH_TBLITE=ON \
  OPENBLAS_DIR=$HOME/opt/openblas-static \
  bash build_crest.sh $HOME/crest-modify $HOME/crest-modify-bin
```

（`FULL=1` 可一次性恢复上游全部特性：tblite + gfnff + libpvol + lwoniom。）

## 3. 计算级别怎么配（重要变化）

| 想算什么 | 配置 |
|---|---|
| **ORCA** | `method = "orca"` —— 你现有的 `orca.toml` / `your-system-orca.toml` **完全不用改** |
| **g-xTB** | `method = "gxtb_dev"`（走外部 xtb 子进程，不受本次改动影响） |
| **gfn2** | ⚠️ `crest x.xyz -gfn2` 现在**不可用**（它走进程内 tblite，已关闭） |

### gfn2 改走外部 xtb 的写法（已验证）

```toml
input = "water.xyz"

[calculation.level]
method = "xtb"                                  # -> 外部 xtb 子进程
binary = "/opt/software/xtb-6.7.1-0/bin/xtb"    # xtb 可执行文件
flags  = "--gfn2"                               # 可选：xtb 默认就是 GFN2
```

实测（水分子单点）：

```
 : xTB calculation via external binary
 :   selected binary : /opt/software/xtb-6.7.1-0/bin/xtb
 TOTAL ENERGY          -5.0703251286 Eh
 CREST terminated normally.
```

关闭 tblite 后 `-gfn2` 的报错是明确的，不会静默算错：

```
Error: Compiled without tblite support!
Use -DWITH_TBLITE=true in the setup to enable this function
ERROR STOP
```

## 4. 为编译通过而修复的 3 处源码问题

这三处都是**原来就存在**的问题，与开关无关；只是关闭 tblite/gfnff 后才暴露出来。
修复已同时应用到本地源码树和服务器。

1. **`src/algos/parallel.f90`（2 处）**
   注释里写了 `MDFILES/*`。CREST 用 `gfortran -cpp` 编译，`/*` 会开启一个
   永不结束的 C 风格注释：
   `错误: 未结束的注释`。
   → 改成 `anything under MDFILES/`。

2. **`src/calculator/gfnff_api.F90`**
   `gfnff_dump_sasa()` 内用了 `ff_dat%solvation`，但该分量只存在于真正的
   GFN-FF 数据类型里；`WITH_GFNFF=OFF` 时用的是占位类型：
   `'solvation' at (1) is not a member of the 'gfnff_data' structure`。
   该函数唯一的调用点（`api_helpers.F90` 的 `gfnff_properties`）本身就在
   `#ifdef WITH_GFNFF` 内，所以关掉时函数不可达。
   → 用 `#ifdef WITH_GFNFF` 包住函数体（含两个局部变量）。

3. **`src/dynamics/dynamics_module.f90`**
   `md_stop_seen()` 与 `md_stop_cleanup()` 被写在 `subroutine dynamics`
   （179–859 行）**内部**，却在模块声明区写了
   `public :: md_stop_seen,md_stop_cleanup`。`public` 只能指向模块实体：
   `Symbol 'md_stop_seen' at (1) has no IMPLICIT type`
   → 把这三个（含辅助的 `md_stopfile_delete`）移到 `end subroutine dynamics`
   之后的模块层。它们只访问模块变量 `stopfile_hit/stopfile_name/stopfile_launchdir`，
   行为不变。

## 5. 运行提示

```bash
export PATH=$HOME/crest-modify-bin/bin:$PATH
# 每个作业独占的临时目录（强烈建议，ORCA/xtb 的 scratch 写入很吃 NFS）
export CREST_SCRATCH_DIR=/tmp/crest-$SLURM_JOB_ID
# 一块 BLAS 线程池就够，CREST 自己用 OpenMP 并行
export OPENBLAS_NUM_THREADS=1
```

验证当前二进制实际带有哪些特性：

```bash
ldd  $HOME/crest-modify-bin/bin/crest | grep -iE 'blas|lapack|not found'
grep -E 'tblitevar|gfn0var|gfnffvar|tomlfvar|libpvolvar|lwoniomvar' \
     $HOME/crest-modify/_build/crest_metadata.fh
```

## 6. 保留的可复用文件（别误删）

| 路径 | 用途 |
|---|---|
| `$HOME/opt/openblas-static/lib/libopenblas.a` | 重新编译必需 |
| `$HOME/opt/buildtools/bin/{cmake,ninja}` | 重新编译必需 |
| `$HOME/rpms2/*.rpm` | OpenBLAS 静态包来源（60 MB，可选） |
| `$HOME/crest-slim.tar.gz` | 上传用的源码包（14 MB，可选） |

已清理的临时目录：`~/rpmstage`（810 MB 解包中间产物）、`~/openblas-build`（截断的源码包）。
