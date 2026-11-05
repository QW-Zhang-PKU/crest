# CREST 3.0.2 本地修改说明

## 一、快照信息

本目录是 CREST 3.0.2 的完整源码树，包含我们为 6011 集群所做的本地修改、离线构建脚本和补丁存档。

- 上游基线：CREST 3.0.2 原始源码包 `crest-latest-backup.zip`
- 基线 SHA-256：`E1D466149F0192761BD07B802B2A966F7AC22746921965B904B6794DB93EE19F`
- 本地源码来源：`/home/workspace/zqw/crest-modify`
- 最新源码文件时间：2026-10-28
- 许可证：LGPL-3.0-or-later，详见 `COPYING`、`COPYING.LESSER`、`LICENSE`、`NOTICE`
- 本快照**不包含**：`_build/`、编译产物、二进制、`.bak-*` 工作副本、日志、测试运行输出、本机路径缓存

## 二、改动总览

以原始 CREST 3.0.2 的 281 个文件为基线，对每个共同文件做 SHA-256 比对，并做路径级集合比对：

| 项目 | 数量 |
|---|---:|
| 修改的上游文件 | **56** |
| 新增源码文件 | **3** |
| 新增本地构建、文档、许可证与补丁文件 | **23** |
| 内置 subproject 依赖文件 | **7,048** |
| 删除的上游文件 | **0** |

`subprojects/` 下的 tblite、toml-f、gfn0、gfnff、lwoniom、libpvol、test-drive 是为**离线编译**而内置的完整依赖源码，不属于我们对 CREST 本体的修改。每个依赖保留自己的许可证文件。

## 三、主要改动分类

### 1. 构建、离线依赖与可移植性

- 内置全部 subproject 源码，`build-now.sh` 可在无网络环境完成构建。
- `build_crest.sh` 提供精简构建开关，可关闭 tblite/gfnff/libpvol/lwoniom，缩短编译时间。
- 支持静态 OpenBLAS 链接，避免“编译机有 OpenBLAS、计算节点没有”的运行时错误。
- 修正关闭 GFN-FF 时 `gfnff_api.F90` 的占位类型编译错误。
- 更新 CMake/Meson 源文件列表，加入新增源码。
- 修正若干现代 gfortran 下的生命周期、长度和类型问题。

### 2. 外部计算器、子进程、信号与线程

- `src/signal.c` 增加 `fork/exec/waitpid` 路径，使外部 ORCA/xtb 子进程处于独立进程组，可被整体终止。
- 信号处理器改为异步信号安全实现，只使用 `write(2)`、`kill()`、`_exit()`；捕获 SIGINT、SIGQUIT、SIGTERM、SIGHUP、SIGXCPU。
- 保留 `nohup` 继承的 SIG_IGN，不会强行覆盖。
- 修正 POSIX `setenv()` 的 3 参数绑定，确保子进程真正拿到 `OMP_NUM_THREADS`、`MKL_NUM_THREADS`、`OPENBLAS_NUM_THREADS`。
- 子进程额外设置 `OMP_MAX_ACTIVE_LEVELS=1`，避免“线程里再开线程”。
- 外部计算器梯度/能量文件路径改为每次按当前 `calcspace` 重新推导，避免并行 worker 读同一个陈旧文件。
- ORCA 模板改为定长字符数组，规避 gcc PR89200/PR113793 导致的堆损坏；ORCA stderr 写入 `ORCA.err`。
- ORCA 的线程数、MPI 参数、`%pal nprocs` 与 `--oversubscribe` 处理更一致。

### 3. 输入、几何与约束安全

- TOML 在坐标读取前解析时增加 `env%ref` 安全闸，避免空引用导致 SIGSEGV。
- `--cinp` 的用户键约束默认登记为刚性 SHAKE 约束；`-cbonds/-cmetal/-cheavy/-clight` 自动生成的键仍为谐振约束。
- 支持 `atoms:`/`elements:` 的 pairwise 约束路径，并增加冗余约束提示。
- SHAKE 修正扫描由 Jacobi 改为 Gauss-Seidel；扫掠预算从 250 提高到 1000。
- SHAKE 失败报告改为 `es12.4` 格式，显示迭代数、偏差、容差、违反约束的键，并区分“目标不可行、迭代爆炸、预算不足、冗余停滞”。
- 拓扑检查前增加预筛；默认排除金属参与的键，减少金属体系误删构象。

### 4. 采样、优化与文件管理

- `--onebatch`：只跑一轮 MTD，适合快速探索。
- `--loadconformer`：从 `readconformers.xyz` 读取已知构象，从第一轮 MTD 开始施加静态 RMSD 偏置。
- `-iru`：把已有构象作为第二个静态偏置叠加到动态 MTD。
- 动态核数再分配：先结束的并行作业把核让给仍在运行的作业；支持 `CREST_DYNAMIC_SHARE`、`CREST_MAX_CORES_PER_JOB`。
- trial MTD 默认缩短为 ORCA 0.125 ps、其他计算器 0.25 ps；支持 `CREST_TRIAL_MD_PS` 覆盖。
- trial MTD 的 Vbias dump 默认固定 50 fs，并保证至少两次更新。
- `stop` 文件支持：在动力学阶段检测到运行期间新建的 `stop` 文件后，干净结束当前动力学并把已采结构交给后续步骤。
- 高频中间文件可迁移到 `CREST_SCRATCH_DIR`；低频中间文件进入 `CREST_WORK_DIR`（默认 `crest_work/`）。
- 轨迹后缀改为 `.xyz`，优化日志改为 `crestopt.xyz`。
- MD/MTD 进度只在 25/50/75% 打印，去掉 99/100% 行；进度行缩短，并显示 elapsed/remaining。
- 修复超过 1 小时时计时打印“模 1 小时”的错误；修复 `60 sec` 进位问题。
- 初始几何优化结束打印圈数和耗时；ORCA 初始优化开始打印 `optimization starting`。

## 四、修改的上游文件（56 个）

```text
docs/man/crest.adoc
src/CMakeLists.txt
src/algos/dynamics.f90
src/algos/numhess.f90
src/algos/parallel.f90
src/algos/refine.f90
src/algos/search_1.f90
src/algos/search_conformers.f90
src/algos/search_entropy.f90
src/algos/search_mecp.f90
src/algos/search_newnci.f90
src/algos/setuptest.f90
src/algos/singlepoint.f90
src/bondconstraint.f90
src/calculator/calc_type.f90
src/calculator/calculator.F90
src/calculator/gfnff_api.F90
src/calculator/orca_sc.f90
src/calculator/subprocess_engrad.f90
src/calculator/subprocess_types.f90
src/calculator/turbom_sc.f90
src/calculator/xtb_sc.f90
src/choose_settings.f90
src/classes.f90
src/cleanup.f90
src/confparse.f90
src/cregen.f90
src/crest_main.f90
src/dynamics/dynamics_module.f90
src/dynamics/metadynamics_module.f90
src/dynamics/shake_module.f90
src/flexi.F90
src/iomod.F90
src/legacy_algos/confscript2_misc.f90
src/legacy_wrappers.f90
src/ls_rmsd.f90
src/meson.build
src/minitools.f90
src/ompmklset.F90
src/optimize/CMakeLists.txt
src/optimize/ancopt.f90
src/optimize/gd.f90
src/optimize/meson.build
src/optimize/optimize_module.f90
src/optimize/optutils.f90
src/optimize/rfo.f90
src/parsing/confparse2.f90
src/parsing/parse_calcdata.f90
src/parsing/parse_maindata.f90
src/parsing/parse_xtbinput.f90
src/printouts.f90
src/restartlog.f90
src/signal.c
src/sigterm.f90
src/strucreader.f90
src/timer.f90
```

## 五、新增源码文件（3 个）

```text
src/exclfile.f90
src/minregistry.f90
src/optimize/lbfgs.f90
```

## 六、新增本地构建、文档、许可证与补丁文件（23 个）

```text
BUILD-NOTES-zqw.md
CHANGED-FILES.txt
LICENSE
LOCAL-CHANGES.md
NOTICE
ORCA-HOWTO.md
README-BUILD.md
README-YUGROUP.md
bench-nmtd.sh
build-now-intel.sh
build-now.sh
build-openblas-static.sh
build-single-file.sh
build_crest.sh
check-orca-threads.sh
check-xtb-threads.sh
crest_env.sh
patches/crest-omp-env-fix.patch
patches/crest-signal-fix.patch
patches/crest-threadfix.patch
patches/crest-toml-segv.patch
patches/crest-xtbsys-mtd-fix.patch
verify-fixes.sh
```

`README-BUILD.md` 和 `BUILD-NOTES-zqw.md` 是历史工作笔记，里面出现较早的“25 个改动文件”等计数属于旧快照；**本文件和 `CHANGED-FILES.txt` 才是当前快照的权威清单**。

## 七、构建与验证

推荐离线构建：

```bash
bash build-now.sh
# 或
PREFIX=$HOME/opt/crest JOBS=16 bash build-now.sh
```

静态 BLAS：

```bash
bash build-openblas-static.sh
OPENBLAS_DIR=$HOME/opt/openblas-static STATIC_BLAS=1 bash build-now.sh
```

常用开关：

```bash
FULL=1 bash build-now.sh                 # 恢复上游可选依赖
WITH_TBLITE=ON bash build-now.sh         # 进程内 GFN2-xTB
STATICBUILD=ON bash build-now.sh         # 全静态，不推荐
```

编译前自检：

```bash
bash verify-fixes.sh
```

## 八、发布与开源合规说明

- 根目录保留上游 `COPYING`（GPL-3.0）和 `COPYING.LESSER`（LGPL-3.0），并新增 `LICENSE`（LGPL-3.0 全文）和 `NOTICE`。
- 所有修改文件和新增源码均保留或补充许可证信息；`src/signal.c` 是少数原本没有文件头的源文件，本次补充了 SPDX 和 LGPL 声明。
- `subprojects/` 下的每个第三方组件都保留自己的 `COPYING`、`LICENSE-*` 等许可证文件；分发时必须一并保留。
- 本快照不包含编译产物或二进制，因此不存在“只发二进制、不提供对应源码”的问题。
- `src/signal.c` 原文件头写有“taken from the internet”，其上游来源没有单独的文件级许可证声明。本 fork 的修改按 LGPL-3.0-or-later 发布；若要做正式法律审查，建议向上游维护者确认该片段的最初来源与授权。
- `README-BUILD.md` 和 `BUILD-NOTES-zqw.md` 中的历史本机路径已统一改为占位符或相对路径；这些文件属于开发笔记，不包含密码或密钥。

---

# English

## 1. Snapshot

This directory is a complete CREST 3.0.2 source tree with the YuGroup 6011 cluster's local fixes, offline build scripts, and patch archive.

- Upstream base: CREST 3.0.2 pristine source archive `crest-latest-backup.zip`
- Base SHA-256: `E1D466149F0192761BD07B802B2A966F7AC22746921965B904B6794DB93EE19F`
- Local source origin: `/home/workspace/zqw/crest-modify`
- Latest source file timestamp: 2026-10-28
- License: LGPL-3.0-or-later; see `COPYING`, `COPYING.LESSER`, `LICENSE`, and `NOTICE`
- Not included: `_build/`, binaries, `.bak-*` working copies, logs, test-run outputs, and machine-specific caches

## 2. Change inventory

Against the 281 files in the pristine archive, every common file was compared by SHA-256 and the path sets were compared:

| Item | Count |
|---|---:|
| Modified upstream files | **56** |
| Added source files | **3** |
| Added local build/documentation/license/patch files | **23** |
| Vendored subproject files | **7,048** |
| Removed upstream files | **0** |

The complete `subprojects/` tree (tblite, toml-f, gfn0, gfnff, lwoniom, libpvol, test-drive) is vendored for **offline builds**. It is not part of the modifications to CREST itself, and each dependency retains its own license files.

## 3. Main areas of change

### 3.1 Build, offline dependencies, and portability

- All subproject sources are vendored; `build-now.sh` can build without network access.
- `build_crest.sh` provides slim-build switches to disable tblite/gfnff/libpvol/lwoniom and shorten the build.
- Static OpenBLAS linking is supported to avoid runtime failures on compute nodes that lack the build host's BLAS.
- A compile error in `gfnff_api.F90` with GFN-FF disabled was fixed.
- CMake/Meson source lists were updated for the new files.
- Several modern-gfortran lifetime, length, and type issues were fixed.

### 3.2 External calculators, child processes, signals, and threads

- `src/signal.c` now uses `fork/exec/waitpid`, puts external ORCA/xtb children in their own process groups, and terminates the whole group.
- Signal handling is async-signal-safe and uses only `write(2)`, `kill()`, and `_exit()`; it handles SIGINT, SIGQUIT, SIGTERM, SIGHUP, and SIGXCPU.
- Inherited `SIG_IGN` is preserved, so `nohup` workflows still work.
- The POSIX `setenv()` binding now has the correct three-argument signature, so child processes receive `OMP_NUM_THREADS`, `MKL_NUM_THREADS`, and `OPENBLAS_NUM_THREADS`.
- Children also receive `OMP_MAX_ACTIVE_LEVELS=1` to avoid nested thread teams.
- External-calculator gradient/energy paths are re-derived from the current `calcspace` on every call, preventing parallel workers from reading a stale file.
- ORCA templates use fixed-length character arrays to avoid gcc PR89200/PR113793 heap corruption; ORCA stderr is captured in `ORCA.err`.
- ORCA thread counts, MPI arguments, `%pal nprocs`, and `--oversubscribe` handling are more consistent.

### 3.3 Input, geometry, and constraint safety

- A safety guard prevents the `env%ref` null-reference SIGSEGV when TOML is parsed before coordinates.
- User bonds from `--cinp` are registered as rigid SHAKE constraints by default; bonds generated by `-cbonds/-cmetal/-cheavy/-clight` remain harmonic restraints.
- The pairwise `atoms:`/`elements:` constraint path is supported, with redundant-constraint diagnostics.
- SHAKE sweeps changed from Jacobi to Gauss-Seidel; the sweep budget increased from 250 to 1000.
- SHAKE failure reports use `es12.4` formatting and show iteration count, deviation, tolerance, violated bonds, and distinguish infeasible targets, divergence, insufficient budget, and redundant stagnation.
- A topology prescreen runs before refinement; metal-involving bonds are excluded from the topology comparison by default to reduce false conformer removal.

### 3.4 Sampling, optimization, and file management

- `--onebatch`: run one MTD round only, for quick exploration.
- `--loadconformer`: read known conformers from `readconformers.xyz` and apply a static RMSD bias from the first MTD round.
- `-iru`: add existing conformers as a second static bias on top of the dynamic MTD.
- Dynamic core redistribution: a finished parallel job hands its cores to jobs still running; controlled by `CREST_DYNAMIC_SHARE` and `CREST_MAX_CORES_PER_JOB`.
- Trial MTD defaults shortened to 0.125 ps for ORCA and 0.25 ps for other calculators; `CREST_TRIAL_MD_PS` overrides.
- Trial-MTD Vbias dump defaults to a fixed 50 fs with at least two updates.
- `stop` file support: during dynamics, a newly created `stop` file makes the next checkpoint end the current dynamics cleanly and pass collected structures to subsequent steps.
- High-frequency intermediates can move to `CREST_SCRATCH_DIR`; low-frequency intermediates go to `CREST_WORK_DIR` (default `crest_work/`).
- Trajectory suffixes changed to `.xyz`; optimization logs changed to `crestopt.xyz`.
- MD/MTD progress is reported only at 25/50/75%; the 99/100% lines are gone. Progress lines are shorter and include elapsed/remaining time.
- The timing bug that printed runtimes modulo one hour was fixed, as was the `60 sec` carry bug.
- Initial geometry optimization now reports cycles and wall time; ORCA initial optimization prints `optimization starting`.

## 4. Modified upstream files (56)

```text
docs/man/crest.adoc
src/CMakeLists.txt
src/algos/dynamics.f90
src/algos/numhess.f90
src/algos/parallel.f90
src/algos/refine.f90
src/algos/search_1.f90
src/algos/search_conformers.f90
src/algos/search_entropy.f90
src/algos/search_mecp.f90
src/algos/search_newnci.f90
src/algos/setuptest.f90
src/algos/singlepoint.f90
src/bondconstraint.f90
src/calculator/calc_type.f90
src/calculator/calculator.F90
src/calculator/gfnff_api.F90
src/calculator/orca_sc.f90
src/calculator/subprocess_engrad.f90
src/calculator/subprocess_types.f90
src/calculator/turbom_sc.f90
src/calculator/xtb_sc.f90
src/choose_settings.f90
src/classes.f90
src/cleanup.f90
src/confparse.f90
src/cregen.f90
src/crest_main.f90
src/dynamics/dynamics_module.f90
src/dynamics/metadynamics_module.f90
src/dynamics/shake_module.f90
src/flexi.F90
src/iomod.F90
src/legacy_algos/confscript2_misc.f90
src/legacy_wrappers.f90
src/ls_rmsd.f90
src/meson.build
src/minitools.f90
src/ompmklset.F90
src/optimize/CMakeLists.txt
src/optimize/ancopt.f90
src/optimize/gd.f90
src/optimize/meson.build
src/optimize/optimize_module.f90
src/optimize/optutils.f90
src/optimize/rfo.f90
src/parsing/confparse2.f90
src/parsing/parse_calcdata.f90
src/parsing/parse_maindata.f90
src/parsing/parse_xtbinput.f90
src/printouts.f90
src/restartlog.f90
src/signal.c
src/sigterm.f90
src/strucreader.f90
src/timer.f90
```

## 5. Added source files (3)

```text
src/exclfile.f90
src/minregistry.f90
src/optimize/lbfgs.f90
```

## 6. Added local build, documentation, license, and patch files (23)

```text
BUILD-NOTES-zqw.md
CHANGED-FILES.txt
LICENSE
LOCAL-CHANGES.md
NOTICE
ORCA-HOWTO.md
README-BUILD.md
README-YUGROUP.md
bench-nmtd.sh
build-now-intel.sh
build-now.sh
build-openblas-static.sh
build-single-file.sh
build_crest.sh
check-orca-threads.sh
check-xtb-threads.sh
crest_env.sh
patches/crest-omp-env-fix.patch
patches/crest-signal-fix.patch
patches/crest-threadfix.patch
patches/crest-toml-segv.patch
patches/crest-xtbsys-mtd-fix.patch
verify-fixes.sh
```

`README-BUILD.md` and `BUILD-NOTES-zqw.md` are historical working notes. Older counts such as “25 files” refer to an earlier snapshot; this file and `CHANGED-FILES.txt` are the authoritative inventory for the current snapshot.

## 7. Build and verification

Recommended offline build:

```bash
bash build-now.sh
# or
PREFIX=$HOME/opt/crest JOBS=16 bash build-now.sh
```

Static BLAS:

```bash
bash build-openblas-static.sh
OPENBLAS_DIR=$HOME/opt/openblas-static STATIC_BLAS=1 bash build-now.sh
```

Common switches:

```bash
FULL=1 bash build-now.sh                 # restore upstream optional dependencies
WITH_TBLITE=ON bash build-now.sh         # in-process GFN2-xTB
STATICBUILD=ON bash build-now.sh         # fully static, not recommended
```

Self-check before building:

```bash
bash verify-fixes.sh
```

## 8. Release and open-source compliance notes

- The root retains upstream `COPYING` (GPL-3.0) and `COPYING.LESSER` (LGPL-3.0), and adds `LICENSE` (full LGPL-3.0 text) and `NOTICE`.
- Modified files and new source files retain or add license information; `src/signal.c` was one of the few files without a file-level header, and now carries an SPDX/LGPL notice.
- Every vendored subproject keeps its own `COPYING`, `LICENSE-*`, or equivalent license files; those must remain with any redistribution.
- This snapshot contains no build products or binaries, so there is no binary-only distribution without corresponding source.
- The original `src/signal.c` header says “taken from the internet” and did not have a separate file-level license declaration. This fork's modifications are released under LGPL-3.0-or-later. For a formal legal review, confirm the provenance and licensing of that original snippet with the upstream maintainers.
- Historical machine-specific paths in `README-BUILD.md` and `BUILD-NOTES-zqw.md` have been replaced with placeholders or relative paths; these files are development notes and contain no passwords or keys.
