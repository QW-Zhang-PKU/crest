# CREST 3.0.2 · QW-Zhang-PKU modified (at Yu group, PKU)

English REDAME is given below the Chinese content.

本目录是 CREST 3.0.2 的修改版，由QW-Zhang完成。本项目旨在解决CREST的长期存在的一些痛点，让其更适配 [有机（尤其是金属有机）] 的机理研究。修改代码时使用了Hy4和Deepseek，修改后的功能已经过测试，具体内容如下。本文件绝大部分内容为作者本人撰写，请放心阅读。

代码修改内容细节在 `LOCAL-CHANGES.md`，许可证与第三方声明在 `NOTICE`、`LICENSE`、`COPYING` 和 `COPYING.LESSER`。

## 一、使用方法

## 二、新增功能（对比CREST 3.0.2）

### 1. 增加 g-xtb、ORCA 接口

**原因** ：CREST内置的方法包括gfn-xtb、gfn2-xtb等。在这个项目开始前，我通常用gfn2-xtb比较多。但由于gfn2-xtb对于过渡金属配合物描述较差，有时会跑出完全错误的结果（比如金属上配位的一氧化碳被插羰，等），所以迫切的需要更高级的方法整合入CREST。
**g-xtb：** g-xtb时Grimme课题组发展的新一代半经验方法，其精度显著好于gfn2-xtb（尤其是针对过渡金属配合物），且计算速度比gfn2-xtb增加不多，所以是我本人目前在CREST中最常用的方法。注意g-xtb目前仅发表在预印本，后续可能仍有修改。且对溶剂化支持不好，所以g-xtb模式下不允许设置溶剂。
**ORCA：** 有时某些体系可能更加棘手，必须DFT作为计算级别（即使很耗时），可以使用ORCA。我们完善并开启了CREST的ORCA接口。经过测试，ORCA下使用纯泛函+极小基（如：! PBE D3 MINIX GCP(DFT/MINIX) LooseSCF）是可以算的动的，但是耗时非常长，仅推荐在疑难杂症时使用。目前只测试了ORCA 6.1.1，老版本ORCA应该兼容。注意ORCA运行时建议设定临时文件夹在SSD上，否则在读写较慢的硬盘（机械硬盘或超算分布式存储）上运行会被拖累。
**更改** ：

- 用 -gxtb_dev 开启g-xtb。需要指定外部的含g-xtb的xtb可执行文件
- ORCA 模板的读写修复了 gfortran 下的生命周期/长度问题，减少生成损坏 `ORCA.in` 的概率。
- ORCA 的标准错误现在会写入 `ORCA.err`，不再被直接丢掉；启动失败时更容易看到真正原因。
- 初始几何优化开始时会打印 `optimization starting`，不再长时间完全静默。

### 1. 长时间 MD/MTD 更不容易白跑

- 修复了并行 worker 读到错误/冻结梯度文件的问题；此前典型表现是 `*MTD n terminated EARLY`，看起来像构象搜索很快就失败。
- 修复了外部 xtb/ORCA 子进程拿不到正确线程数的问题，减少 CPU 空转和线程过度订阅。
- Ctrl+C、`kill`、SLURM 超时信号现在能真正终止 CREST 及其外部子进程，不再出现终端“假死”、只能 `kill -9` 的情况。
- 增加 `stop` 文件：任务运行中在工作目录执行 `touch stop`，最近的检查点会结束当前动力学，把已经采到的结构交给后续优化/筛选，而不是直接丢弃。

**对你的好处**：长任务中途想收手时，不必再等它跑完或强杀；已经得到的构象可以保留下来。因信号处理、线程预算或梯度文件导致的莫名失败会明显减少。

### 3. 输出更适合盯任务

- MD/MTD 进度固定显示 **25% / 50% / 75%**，不再出现没有意义的 99%/100% 行。
- 进度行更短，适合窄终端；并行任务的多行输出不容易折行错位。
- 显示已用时间和预计剩余时间（`elapsed ... remaining ...`）。
- 轨迹和优化日志改用通用格式后缀：
  - `crest_dynamics.trj` → `crest_dynamics.xyz`
  - `crest_<n>.trj` → `crest_<n>.xyz`
  - `crestopt.log` → `crestopt.xyz`
- 初始几何优化结束后会报告优化圈数和耗时，例如：

```text
Initial geometry optimization: 37 cycles, 12 sec   (crestopt.xyz)
```

**对你的好处**：VMD、Avogadro、Mercury、ASE 等工具可以直接打开轨迹；不用等任务结束才知道优化走了多少圈；判断“卡住还是慢”更容易。

### 4. 工作目录更干净，NFS 压力更小

- 高频读写的中间文件可以放到节点本地临时目录：

```bash
export CREST_SCRATCH_DIR=/tmp/crest-$SLURM_JOB_ID
```

- 低频读写的中间文件默认收进工作目录下的 `crest_work/`：

```bash
export CREST_WORK_DIR=crest_work    # 默认
export CREST_WORK_DIR=none          # 回到上游布局
```

- `crest_conformers.xyz`、`crest_rotamers.xyz`、`crest_best.xyz` 等结果文件仍然留在当前工作目录。

**对你的好处**：工作目录不会被大量 `coord`、`wbo`、`MDFILES`、`OPTIM` 等中间文件淹没；把 scratch 放到 `/tmp` 还能减少 NFS 往返，外部计算器较慢时尤其明显。

### 5. 约束更符合化学直觉

- `--cinp` 里的键约束默认按**刚性 SHAKE 约束**处理，而不是很硬的谐振弹簧。固定反应键时不需要再猜力常数，也不容易因为步长过大导致 `SHAKE did not converge`。
- `-cbonds`、`-cmetal`、`-cheavy`、`-clight` 自动生成的键约束仍然是**谐振约束**，不会把整个分子冻成刚体。
- 默认在拓扑检查中排除金属原子参与的键，减少金属体系因为配位数取整跳变而被误删构象。
- SHAKE 的迭代格式、收敛预算和失败报告都做了改进；失败时会打印偏差、容差、违反约束的键等信息。

**对你的好处**：过渡态反应键、金属配合物和含自动键约束的体系更不容易出现“约束太硬导致轨迹炸掉”或“构象被误判拓扑变化而删掉”。

### 6. 采样流程更省时间

- trial MTD 更短：ORCA 默认 0.125 ps，其他计算器默认 0.25 ps；可用 `CREST_TRIAL_MD_PS` 覆盖。
- trial MTD 的 Vbias dump 默认固定 50 fs，并保证至少两次更新。
- 并行 MTD 中，先完成的作业可以把核让给仍在运行的作业（`CREST_DYNAMIC_SHARE=0` 可关闭；`CREST_MAX_CORES_PER_JOB` 可设上限）。
- 新增/接通：
  - `--onebatch`：只跑一轮 MTD，适合想快速探索、不想等完整多轮迭代的任务。
  - `--loadconformer`：从 `readconformers.xyz` 读入已知构象，从第一轮 MTD 就加静态 RMSD 偏置。
  - `-iru`：把已找到的构象作为第二个静态偏置叠加到动态 MTD 上。

**对你的好处**：在外部 ORCA 这种“每次梯度都贵”的场景下，可以减少 trial 和空置核心的浪费；有已知构象时，可以让新一轮搜索从一开始就往未探索区域走。

## 二、升级后需要知道的差异

1. **轨迹后缀变成 `.xyz`**  
   旧脚本如果硬编码 `.trj` 或 `crestopt.log`，需要改成 `.xyz`。

2. **trial MTD 更短**  
   ORCA 默认 0.125 ps，其余默认 0.25 ps。优点是省时间；代价是时间预估会更抖，而且只在大约 50 步之后才暴露的不稳定性可能漏检。需要保守时用 `CREST_TRIAL_MD_PS=0.5` 或更长。

3. **`--cinp` 默认是刚性约束**  
   如果你确实需要旧的谐振软约束，用 `-softcons`，或在 TOML 中设 `hardconstraints = false`。

4. **Vbias dump 默认值变了**  
   生产 MTD 默认 500 fs；trial 默认 50 fs。旧任务如果依赖原来的 1000 fs 间隔，请显式设置 `-vbdump`/`CREST_TRIAL_VBDUMP_FS`。

5. **动态核数再分配**  
   对 ORCA/MPI 不一定总是核越多越快；48 个 MPI rank 在小体系上可能比 16–24 个更慢。可以用：

```bash
export CREST_MAX_CORES_PER_JOB=24
```

6. **`stop` 文件只对动力学阶段生效**  
   多层优化和 cregen 阶段不会检查它；已经进入这些阶段后，`touch stop` 不会立刻结束整个任务。

## 三、快速开始

编译或取得可执行文件后：

```bash
source ./crest_env.sh 40
crest molecule.xyz -T 40 --quick
```

使用外部 xtb 的 TOML：

```toml
input = "molecule.xyz"

[calculation.level]
method = "xtb"
binary = "/path/to/xtb"
flags  = "--gfn2"
```

使用 ORCA 时，按 `ORCA-HOWTO.md` 配置 `method = "orca"` 和模板。长任务建议在每个作业内使用独立的 scratch：

```bash
export CREST_SCRATCH_DIR=/tmp/crest-$SLURM_JOB_ID
export OPENBLAS_NUM_THREADS=1
```

运行前可以自检：

```bash
bash verify-fixes.sh
```

## 四、适用范围与引用

- 适用于 Linux 集群上的 CREST 3.0.2 工作流，尤其是外部 xtb / ORCA 计算器。
- 这些改动主要是稳定性、性能、诊断和文件管理改进；在算法未改变的部分，科学结果应与上游一致。SHAKE 的迭代格式和预算变化属于数值收敛改进，涉及强约束体系时建议用自己的测试体系复核。
- 如果用于发表，请按上游 `CITATION.cff` 和 CREST 论文引用 CREST；本地修改不改变上游 authorship。
- 许可证：LGPL-3.0-or-later。详见 `COPYING`、`COPYING.LESSER`、`LICENSE`、`NOTICE`。

---

# CREST 3.0.2 · YuGroup Edition

This directory is a locally modified CREST 3.0.2 tree for chemists who use external xtb or ORCA for conformer searches. It does not introduce a new conformer-search theory. The changes focus on **job robustness, external-program integration, readable output, file management, and parts of the sampling workflow**.

If you only want to know what changes for you as a user, this document is enough. File-by-file technical details are in `LOCAL-CHANGES.md`; licensing and third-party notices are in `NOTICE`, `LICENSE`, `COPYING`, and `COPYING.LESSER`.

## 1. What changes in practice

### 1.1 Long MD/MTD runs fail less often

- Parallel workers no longer read a stale or frozen gradient file. The old symptom was an early `*MTD n terminated EARLY`, which looked like a quick failure of the conformer search.
- External xtb/ORCA child processes now receive the intended thread budget, reducing idle cores and oversubscription.
- Ctrl+C, `kill`, and SLURM timeout signals now actually stop CREST and its external child processes. The terminal no longer appears frozen until a `kill -9`.
- A `stop` file is supported: while the job is running, `touch stop` in the working directory makes the next checkpoint end the current dynamics and pass the structures already collected to the following optimization/filtering steps.

**Benefit:** you can stop a long run early and keep useful structures instead of losing everything. Failures caused by signal handling, thread budgets, or gradient-file reuse should be much less frequent.

### 1.2 More reliable ORCA use

- ORCA thread counts, MPI arguments, and `%pal nprocs` handling are more consistent.
- ORCA template read/write handling was fixed for a gfortran lifetime/length problem, reducing the chance of a corrupted `ORCA.in`.
- ORCA standard error is now written to `ORCA.err` instead of being discarded.
- `optimization starting` is printed before the initial geometry optimization, so the log no longer goes completely silent for a long time.

**Benefit:** ORCA failures are easier to diagnose; you no longer have to guess from a generic “Initial geometry optimization failed” message.

### 1.3 Output that is easier to monitor

- MD/MTD progress is reported at **25% / 50% / 75%**; the old meaningless 99%/100% lines are gone.
- Progress lines are shorter and less likely to wrap on a narrow terminal or to interleave badly across parallel jobs.
- Elapsed time and remaining time are shown (`elapsed ... remaining ...`).
- Trajectory and optimization-log suffixes now use common formats:
  - `crest_dynamics.trj` → `crest_dynamics.xyz`
  - `crest_<n>.trj` → `crest_<n>.xyz`
  - `crestopt.log` → `crestopt.xyz`
- The initial geometry optimization now reports cycles and wall time, for example:

```text
Initial geometry optimization: 37 cycles, 12 sec   (crestopt.xyz)
```

**Benefit:** VMD, Avogadro, Mercury, ASE, and similar tools can open the trajectories directly. You can see how far an optimization has progressed and distinguish “slow” from “stuck” much more easily.

### 1.4 Cleaner working directories and less NFS traffic

- High-frequency intermediate files can be placed on node-local scratch:

```bash
export CREST_SCRATCH_DIR=/tmp/crest-$SLURM_JOB_ID
```

- Low-frequency intermediate files go into `crest_work/` by default:

```bash
export CREST_WORK_DIR=crest_work    # default
export CREST_WORK_DIR=none          # restore the upstream layout
```

- Result files such as `crest_conformers.xyz`, `crest_rotamers.xyz`, and `crest_best.xyz` stay in the working directory.

**Benefit:** the working directory is no longer flooded with `coord`, `wbo`, `MDFILES`, `OPTIM`, and other intermediate files. Putting scratch on `/tmp` also reduces NFS round trips, which is especially noticeable when the external calculator is slow.

### 1.5 Constraints that match chemical intuition better

- Bond constraints read from `--cinp` are now **rigid SHAKE constraints** by default, instead of very stiff harmonic springs. Fixing a reaction coordinate no longer requires guessing a force constant, and it is less likely to cause `SHAKE did not converge` when the time step is large.
- Bonds generated automatically by `-cbonds`, `-cmetal`, `-cheavy`, and `-clight` remain **harmonic restraints**, so they do not freeze the whole molecule into a rigid body.
- Metal-involving bonds are excluded from the topology check by default, reducing the loss of conformers when a metal coordination number changes by rounding.
- SHAKE iteration, convergence budget, and failure reporting were improved; failure reports now include deviations, tolerances, and the affected bonds.

**Benefit:** transition-state reaction coordinates, metal complexes, and systems with automatic bond constraints are less likely to fail from an over-stiff constraint or to lose conformers because of a topology false positive.

### 1.6 Less wasted sampling time

- Trial MTD is shorter: 0.125 ps for ORCA and 0.25 ps for other calculators by default. `CREST_TRIAL_MD_PS` overrides it.
- The trial-MTD Vbias dump is fixed at 50 fs by default, with at least two updates.
- In parallel MTD, a finished job can hand its cores to jobs still running (`CREST_DYNAMIC_SHARE=0` disables this; `CREST_MAX_CORES_PER_JOB` caps the per-job share).
- New or newly wired options:
  - `--onebatch`: run only one MTD round, for a quick exploration without the full multi-round workflow.
  - `--loadconformer`: read known conformers from `readconformers.xyz` and apply a static RMSD bias from the first MTD round.
  - `-iru`: add the conformers found so far as a second static bias on top of the dynamic MTD.

**Benefit:** when each ORCA gradient is expensive, trial-MTD time and idle-core time are reduced. If you already have known conformers, a new search can start by exploring regions not covered by them.

## 2. Compatibility notes

1. **Trajectory suffixes changed to `.xyz`**  
   Scripts that hard-code `.trj` or `crestopt.log` must be updated.

2. **Trial MTD is shorter**  
   Defaults are 0.125 ps (ORCA) and 0.25 ps (other calculators). This saves time, but runtime estimates are noisier, and instabilities that appear only after roughly 50 steps may be missed. Use `CREST_TRIAL_MD_PS=0.5` or longer when conservative behavior is required.

3. **`--cinp` uses rigid constraints by default**  
   Use `-softcons` or TOML `hardconstraints = false` if you need the old harmonic behavior.

4. **Vbias dump defaults changed**  
   Production MTD defaults to 500 fs; trial MTD defaults to 50 fs. If a workflow depends on the old 1000 fs interval, set `-vbdump` or `CREST_TRIAL_VBDUMP_FS` explicitly.

5. **Dynamic core redistribution**  
   More MPI ranks are not always faster for ORCA, especially for small systems. A cap can be set:

```bash
export CREST_MAX_CORES_PER_JOB=24
```

6. **`stop` applies to dynamics only**  
   Multilevel optimization and cregen do not check it. Once those stages have started, `touch stop` will not immediately end the whole job.

## 3. Quick start

After building or obtaining the executable:

```bash
source ./crest_env.sh 40
crest molecule.xyz -T 40 --quick
```

External xtb TOML example:

```toml
input = "molecule.xyz"

[calculation.level]
method = "xtb"
binary = "/path/to/xtb"
flags  = "--gfn2"
```

For ORCA, configure `method = "orca"` and the template as described in `ORCA-HOWTO.md`. For long jobs, use a per-job scratch directory:

```bash
export CREST_SCRATCH_DIR=/tmp/crest-$SLURM_JOB_ID
export OPENBLAS_NUM_THREADS=1
```

A self-check can be run before production:

```bash
bash verify-fixes.sh
```

## 4. Scope, citation, and license

- Intended for Linux clusters running CREST 3.0.2, particularly with external xtb or ORCA calculators.
- The changes are mainly robustness, performance, diagnostics, and file-management improvements. Where the algorithm was not changed, scientific results should agree with upstream CREST. The SHAKE iteration/budget changes are numerical-convergence improvements; validate them on your own strongly constrained systems.
- For publications, cite CREST according to the upstream `CITATION.cff` and the CREST papers. Local modifications do not change upstream authorship.
- License: LGPL-3.0-or-later. See `COPYING`, `COPYING.LESSER`, `LICENSE`, and `NOTICE`.
