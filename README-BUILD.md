# crest-now — 已应用全部本地修复的 CREST 3.0.2 源码

这个目录是**可以直接拿去 Linux 编译**的完整 CREST 源码树：
所有本地修复都已经写进 `src/`，并且依赖子项目（tblite / toml-f / gfn0 / gfnff /
lwoniom / pvol / test-drive）的真实源码已经就位，**编译全程不需要联网**。

- 文件数：7329，体积约 74 MB
- 校验状态：与参考树逐字节一致（25 个改动文件全部比对通过）

---

## 一、快速开始

### 方式 A（推荐）：先把整个目录拷到 Linux 家目录再编译

```bash
# 在 Linux 上
cp -a /PATH/TO/crest-code/crest-now ~/crest-now
cd ~/crest-now
bash build-now.sh                 # 默认安装到 $HOME/crest-dist
```

### 方式 B：直接在 /PATH/TO/mount 下跑

```bash
cd /PATH/TO/crest-code/crest-now
bash build-now.sh
```

`build-now.sh` 会自动识别「当前目录在 `/mnt/` 之类的慢速挂载上」，
先把源码 stage 到 `$HOME/crest-now` 再在那里编译 —— 万把个文件走 9p 桥接编译会慢到没法忍。

### 可选参数

```bash
PREFIX=~/opt/crest bash build-now.sh    # 自定义安装前缀
JOBS=16 bash build-now.sh               # 限制编译并行度
```

### 静态链接（编译机 ≠ 计算节点时必须看）

CREST 默认**动态**链接 OpenBLAS。如果编译机和跑任务的节点软件栈不同，
计算节点上会直接起不来：

```
crest: error while loading shared libraries: libopenblas.so.0: cannot open shared object file
```

两种解法，按推荐度排序：

**A. 只把 BLAS/LAPACK 静态链进去（推荐）** —— 先编一份静态 OpenBLAS，
再重编 CREST。`libc` / `libgomp` / `libgfortran` 保持动态，
它们属于操作系统，同集群各节点都有，几乎不会编挂。

```bash
bash build-openblas-static.sh                       # 装到 ~/opt/openblas-static
OPENBLAS_DIR=$HOME/opt/openblas-static STATIC_BLAS=1 bash build-now.sh
```

若链接时报 `undefined reference to 'pthread_create'`（pthread 版 `libopenblas.a`
偶发），加上线程库再重来：

```bash
OPENBLAS_DIR=$HOME/opt/openblas-static STATIC_BLAS=1 \
  STATIC_BLAS_EXTRA_LIBS="-lpthread -lm" bash build-now.sh
```

`build-openblas-static.sh` 里两个开关不可省：

| 开关 | 作用 |
|---|---|
| `NO_SHARED=1` | 只产 `libopenblas.a`，避免 CMake 又链回 `.so` |
| `DYNAMIC_ARCH=1` | **运行时**按 CPU 派发内核。编译机是 AVX512 而节点只有 AVX2 时，没有它会 `illegal instruction` |

**B. 全静态（CREST 自带开关，不推荐）** —— `STATICBUILD=ON` 会连 glibc
一起静态（`CMakeLists.txt:48`），需要 `glibc-static` 包，且静态二进制配
OpenMP 会有告警。只在 A 走不通时试：

```bash
STATICBUILD=ON bash build-now.sh
```

编完务必确认（期望：没有 `blas`/`lapack`，也没有 `not found`）：

```bash
ldd $PREFIX/bin/crest | grep -i -E 'blas|lapack|not found'
```

`build_crest.sh` 在装完后会自动跑一遍这个检查。

**C. 不想重编的临时办法** —— 把编译机上的 `libopenblas.so.0`（含符号链接）
拷到 crest 安装目录的 `../lib/` 下，SLURM 脚本已经把 `$CRESTROOT/lib`
加进了 `LD_LIBRARY_PATH`。两分钟就能跑起来，但换节点可能再犯。

**D. 一步到位：产出单个自包含可执行文件** —— `bash build-single-file.sh`。
它不重新实现 cmake，只是把链接模式钉死（`STATIC=intel` +
`LA=mkl-static`）后调 `build-now-intel.sh`，然后**验收并打包**：

```bash
bash build-single-file.sh                 # 结果：~/crest-single/crest
DIST=~/dist JOBS=16 bash build-single-file.sh
TOOLCHAIN=gnu bash build-single-file.sh   # gfortran + 静态 OpenBLAS
BIN=$PREFIX/bin/crest bash build-single-file.sh   # 不重编，只验收+打包
```

三道关卡，任一道不过就非零退出：

| 关卡 | 判据 | 典型失败 |
|---|---|---|
| 依赖门禁 | `ldd` 里只剩 glibc（`libc/libm/libdl/libpthread`） | 出现 `libmkl_*` / `libiomp5` / `libifcore` / `libopenblas` / `libgomp` → 说明 `STATIC`/`LA` 没生效 |
| OpenMP 检查 | `flags.make` 里有 `-qopenmp`（`-fopenmp`） | 静默变成单线程，这是最坏的结果 |
| 干净环境启动 | `env -i ./crest --version` 能跑 | 还需 `LD_LIBRARY_PATH` 就不算自包含 |

**"单文件"的边界（别误会）**：`crest` 本来就只有一个 ELF
（`CMakeLists.txt:281-286`，`crest-exe` 链 `lib-crest-static`，安装项里
只有 `bin/crest` 是运行必需的）。差别在**链接方式**，不在文件个数。
产物仍依赖 glibc（⇒ 在最老的节点上编译），并且 **xtb / orca 永远是外部
程序**，由 `system()` 调起（`env%ProgName`，`src/confparse.f90:137`），
不可能链进去。

不要用 `-DSTATICBUILD=ON`（`config/CMakeLists.txt:41`）来"省事"：它等价于
全静态 `-static`，而且会关掉 object 库与共享库、让 CREST 自己的目标不再
直接链 `OpenMP::OpenMP_Fortran`（`CMakeLists.txt:190-200`）——`-qopenmp`
只剩从 `tblite::tblite` 传递过来这一条路
（`subprojects/tblite/CMakeLists.txt:63-95` 的 `PUBLIC ${lib-deps}`）。
本树带 tblite 所以没事，但一旦去掉 tblite 就会静默变成单线程。

### 编译完成后

```bash
source ./crest_env.sh 40                # 把 crest 加进 PATH 并设置线程预算
cd ~/work/your-system                      # 工作目录务必放原生 ext4，不要放 /PATH/TO/mount
crest your-system-xtb.toml -T 40 --quick --noreftopo
```

运行时另开一个终端确认子进程真的拿到了线程数：

```bash
watch -n 2 ./check-xtb-threads.sh
```

期望看到：`OMP_NUM_THREADS=6  MKL_NUM_THREADS=6  OPENBLAS_NUM_THREADS=6`
`OMP_MAX_ACTIVE_LEVELS=1`，每个 xtb 进程约 6~8 个线程。

---

## 二、目录内容

```
crest-now/
├── src/  include/  config/  docs/  examples/  test/  assets/   # CREST 本体（已打补丁）
├── subprojects/                     # 7 个依赖的真实源码（原为空占位目录）
│   ├── tblite/  toml-f/  gfn0/  gfnff/  lwoniom/  pvol/  test-drive/
├── patches/                         # 本次全部改动的补丁存档
│   ├── crest-xtbsys-mtd-fix.patch   # 外部计算器梯度文件被并行 worker 复用
│   ├── crest-omp-env-fix.patch      # 子进程线程预算没送出去（setenv 绑定少一个参数）
│   └── crest-threadfix.patch        # 25 文件聚合补丁（= 上面两个 + 早期线程修复）
├── build-now.sh                     # ← 一键编译入口
├── build_crest.sh                   # 通用 CMake 构建脚本（被 build-now.sh 调用）
├── build-openblas-static.sh         # 编静态 OpenBLAS（编译机 ≠ 计算节点时用）
├── verify-fixes.sh                  # ← 编译前自检：确认补丁都在、依赖齐全
├── crest_env.sh                     # 运行时环境（PATH / OMP / BLAS / XTBPATH）
├── check-xtb-threads.sh             # ← 运行时诊断：看每个 xtb 实际有几个线程
├── bench-nmtd.sh                    # ← 扫描 -nmtd，找出真正吃满核的 MTD 数量
├── your-system-xtb.toml                # 外部 xtb 输入示例
├── your-system-orca.toml               # ← 外部 ORCA 输入示例
├── orca.inp                         # ← ORCA 输入模板（计算级别 / 基组）
└── ORCA-HOWTO.md                    # ← 用 ORCA 当计算器的完整说明
```

**建议先跑一次自检**（`build-now.sh` 也会自动调用，所以并非必须）：

```bash
bash verify-fixes.sh
# RESULT: PASS - 说明补丁齐全、依赖齐全，可以离线编译
```

---

## 三、已应用的修复

| # | 文件 | 问题 |
|---|---|---|
| 1 | `src/parsing/parse_calcdata.f90` | TOML 在读取坐标之前被解析，`env%ref` 为空时调用 `axis()` 导致 SIGSEGV。加了安全闸，并要求 TOML 里写 `input = "xxx.xyz"` |
| 2 | `src/calculator/xtb_sc.f90`<br>`orca_sc.f90` / `turbom_sc.f90` | 梯度/能量文件路径被缓存，`parallel.f90` 换 `calcspace` 时没清掉，导致所有并行 worker 读同一个陈旧文件 → MD 用"冻结梯度"积分 → `*MTD n terminated EARLY`。改为每次按当前 `calcspace` 重新拼路径 |
| 3 | `src/algos/parallel.f90` | 新增 `crest_calc_is_subprocess()`：外部计算器模式下 CREST 线程本来就阻塞在 `waitpid`，再开嵌套 OMP team 纯属重复占线程 |
| 4 | `src/dynamics/metadynamics_module.f90` | MTD 的 RMSD 计算每个 MD 步都要跑，参考结构少时开 OMP team 反而更慢（线程空转吃 sys 时间）。只有参考结构 ≥ 32 个才开 |
| 5 | `src/sigterm.f90` | 信号处理在每个 energy/gradient 调用里都重装一次，多线程下反复抢进程级 sighand 锁 → 纯内核时间。改为只装一次 |
| 6 | `src/iomod.F90` | **`setenv` 的 C 绑定只声明了两个参数**，POSIX 实际是三个；`overwrite` 从脏寄存器读，一旦为 0 就拒绝覆盖已存在的变量 → `ompenvset()` 改不动 `OMP/MKL/OPENBLAS_NUM_THREADS` |
| 7 | `src/ompmklset.F90` | 配合 #6，并额外给子进程导出 `OMP_MAX_ACTIVE_LEVELS=1`，防止 xtb「线程套线程」把外层 N 变成 N×N |
| 8 | `src/calculator/orca_sc.f90` | **ORCA 的 `rd_ORCA_wbo()` 原来只是个报错桩**（`**ERROR** ORCA bond order reader not yet implemented`）。`-shake`、NCI 模式以及 `wbo`/`flexi` 工具一旦向 ORCA 要 WBO 就会直接崩。现在改为解析 `ORCA.out` 里的 Mayer 键级块（`B( 0-C , 1-S ) : 1.0550`），把 ORCA 的 **0 基**原子编号转成 CREST 的 1 基约定、对称化矩阵，找不到键级时只告警不中止 |
| 9 | `src/calculator/orca_sc.f90` | **ORCA 的 stderr 被 `2>/dev/null` 丢掉**（沿用了其他子进程计算器的 `dev0`）。ORCA 启动失败时（`LD_LIBRARY_PATH` 没设 → `libmpi.so.40: cannot open shared object file`、二进制路径写错、授权问题）消息全部走 stderr，屏幕上什么都看不到，最终只报一句毫无信息的 `Initial geometry optimization failed!`。改为 `2> ORCA.err`，错误落在 `calculation.level.<n>/ORCA.err` |
| 10 | `src/calculator/subprocess_types.f90` | **CREST 把 `ORCA.in` 写成二进制垃圾**。`orca_input` 里 `character(len=:),allocatable :: input(:)` 是「**延迟长度字符数组**」作为派生类型分量；派生类型赋值时 GCC 会把它静默拷错（**gcc PR fortran/89200，GCC 9 的回归**，8.2.0 正常；ifort/ifx 也有同类问题）。CREST 会拷贝 `calculation_settings`（`setuptest.f90:306 tmpcalc = env%calc`、解析 TOML 时的 `calc%add`），拷贝后那些元素不再指向模板而是指向陈旧的堆内存 → `ORCA.in` 里出现 NUL 字节、CREST 自己的 banner、`jobdescription()` 表、systemcall 字符串。**全代码库只有这一处「延迟长度字符数组」分量，所以只有 ORCA 这一条路中招**（xtb/turbomole/generic 不受影响，`coord`、`crest_input_copy.xyz` 也都正常，可作对照）。改法：元素长度定长 `character(len=orca_linelen),allocatable :: input(:)`，纯 memcpy 拷贝，不可能拷错 |
| 11 | `src/calculator/subprocess_types.f90` | **修完 #10 后立刻 `double free or corruption (out)` / SIGABRT**。把元素改成定长（1056）后，原 `allocate(self%input(nlines),source=repeat(" ",width))` 的 SOURCE 长度（来自 `getlines()`，约 85）与声明的元素长度不再相等。**gcc PR fortran/113793**：gfortran ≤ 13 对此不报错，反而按 SOURCE 长度分配元素，随后 `self%input(k) = trim(atmp)` 每次写 1056 字节进 ~85 字节的元素 → 堆越界 → 下一次 `free()` 直接 abort（GCC 14 修）。改法：SOURCE 长度必须等于元素长度，即 `source=repeat(" ",orca_linelen)`。**#10 和 #11 必须成对修**，只修一半就会看到 `double free` 落在 `subprocess_types.f90` |
| 12 | `src/calculator/subprocess_types.f90`<br>`src/calculator/orca_sc.f90` | **ORCA 始终只用 1 个核**。CREST 是通过 `OMP_NUM_THREADS` 把「每任务核数预算」交给外部程序的（`new_ompautoset()`：`max`→T，`auto`,n→T/n），但 **ORCA 走 MPI、压根不看这个变量**，不写 `%pal nprocs N` 就永远是串行。改法：写 `ORCA.in` 时把预算读回来并自动生成 `%pal` 块，于是初始优化拿到 T 个核、MTD/MD 每个 worker 拿到 T/nsim 个核、批量优化拿到 `CREST_CORES_PER_JOB` 个核，全程无需配置；可用 `CREST_ORCA_NPROCS` 手动封顶。同时把子进程的 OMP/MKL 层钉成 1，避免 N 个 rank × N 个线程 = N² |
| 13 | `src/ompmklset.F90`<br>`src/algos/parallel.f90` | **批量优化退化成「每个任务 1 核」**。`crest_oloop` / `crest_sploop` 用 `'auto',nall` 分配，而 `'auto'` 只在 `T > maxjobs` 时才重新分配——构象数 `nall` 几乎总是远大于 T，条件不成立 → `cores_per_job = 1`。新增 `'fixed'` 模式：每任务 `CREST_CORES_PER_JOB` 核（默认 4），并发数 `T/核数`，其余任务由 OpenMP task 队列自然排队。结构数少于 `T/4` 时余核分给已有任务，不浪费 |
| 14 | `src/ompmklset.F90`<br>`src/algos/setuptest.f90` | **trial MTD 强制单核 + 时间预估失真**。原代码 `new_ompautoset(env,'min',...)`（注释写着 "to simulate max parallelization"）→ trial MTD 全程串行；而随后的 `timeestimate` 又 **假设** `floor(T/nmetadyn)` 倍完美线性加速，两头都不准。改法：新增 `'share'` 模式，让 trial 拿生产阶段同样的 `T/nmetadyn` 核；估算器改为「直接用测得值」+ 按 `parallel`（= 实际并发数）算波次，不再凭空除一个数。对 ORCA 这类几乎无收益的标定步骤，这是纯粹的省时间 |
| 15 | `src/algos/setuptest.f90` | **trial MTD 缩短到 1/4（1.0 ps → 0.25 ps）**。这一步只做两件事：校验时间步/SHAKE 设置，以及给时间预估提供实测值，本身没有科学产出，所以不必跑满。配套三处改动：① `MDSTART%length_ps` 用变量 `trial_md_ps`（默认 0.25）；② `MTD%cvdump_fs` 从写死的 550 fs 改成 `550*trial_md_ps` —— **否则 0.25 ps 内一次 Vbias 更新都不会发生，trial 就退化成纯 MD 了**；③ `timeestimate()` 新增 `trialps` 参数，按 `mdtime/trialps` 缩放，而不是沿用「测得值就是 1 ps」的隐含假设。可用 `CREST_TRIAL_MD_PS=0.5` 覆盖，不用重编 |
| 16 | `src/confparse.f90` | **正式 MTD 的 Vbias dump 默认值 1000 fs → 500 fs（1 ps → 0.5 ps）**。注意这是**生产**设置，和 #15 里 trial 的那个不是同一个变量：生产走 `env%mddump`（默认值在 `confparse.f90:222`）→ `parallel.f90` 读成 `mtd%cvdump_fs`；trial 走 `setuptest.f90` 里的局部值（上游就是 550 fs，与生产的 1000 fs 本来就不一致）。命令行 `-vbdump <ps>` 依旧可覆盖。副作用：`cvdumpstep` 减半，`metadynamics_module.f90` 里按 `log(3)/(0.5*dum2)` 自动算的 ramp 参数随之加倍 —— 这是设计如此，保证阻尼相对 dump 间隔保持一致 |
| 17 | `src/calculator/orca_sc.f90`<br>`src/calculator/subprocess_types.f90` | **ORCA 子进程的线程数 / MPI 参数可外部调节**（原来写死 `OMP_NUM_THREADS=1`，调用固定带 `--oversubscribe`）。用于排查「CREST 里 ORCA 比手工跑慢 2~3 倍」：① `CREST_ORCA_OMP_THREADS=N` 控制子进程线程数（默认仍 1）；② `CREST_ORCA_NO_OVERSUBSCRIBE=1` 去掉 `--oversubscribe`（该标志会让 OpenMPI 关掉进程绑定）；③ `CREST_ORCA_MPI_ARGS="--bind-to core"` 追加任意参数。配合已有的 `CREST_ORCA_NPROCS=1` 即可测「1 进程 + 8 线程」对小基组是否比「8 个 MPI rank」更快。详见 `ORCA-HOWTO.md` 第 10 节 |
| 18 | `src/iomod.F90`<br>`src/calculator/calc_type.f90`<br>`src/cleanup.f90` | **支持把 `calculation.level.N` 挪到节点本地盘**（新增 `CREST_SCRATCH_DIR`）。CREST 默认在**当前工作目录**下给每个单点建 `calculation.level.N`，工作目录在 NFS 机械盘上时，ORCA 每次写 scratch 都是网络往返，对 MINIX 这种便宜级别往往是主要开销。设 `CREST_SCRATCH_DIR=/tmp/crest-$SLURM_JOB_ID` 后目录变成 `<root>/calculation.level.N`。实现要点：① `iomod::crest_scratch_dir()` 用 `mkdir -p` 建根（**`iomod%makedir()` 是裸 `mkdir()`，不会建父目录**）；② `calc_type.f90` 的 `calculation_settings_autocomplete` 里拼前缀，TOML 里显式写 `calcspace` 仍然优先；③ `cleanup.f90` 只删 `<root>/calculation.level.*`，**绝不删 root 本身**（否则 `CREST_SCRATCH_DIR=/tmp` 会出事）。**root 必须每个作业唯一**，否则两个 CREST 会抢同一个 `.../calculation.level.1`。**先 `df -hT /tmp` 确认 /tmp 不是 tmpfs**，tmpfs 占的是内存配额，`%maxcore` 大的 ORCA 会把作业撑爆 |
| 19 | `src/algos/setuptest.f90` | **初始优化（Initial Geometry Optimization）结束后报告耗时与优化圈数**。上游在 `smallhead('Initial Geometry Optimization')` 和 `Geometry successfully optimized.` 之间**一个字都不打**（`pr=.false.`），用慢的外部程序时看起来就像卡死了，事后也没有任何数字可比。改法：用 `timer` 计时，并数 `crestopt.log` 里的帧数——三种优化引擎（ANCOPT / RFO / 梯度下降）都是**每做一次梯度计算就 append 一帧**，所以帧数 == 优化圈数。输出形如 `Initial geometry optimization: 37 cycles, 12 sec   (crestopt.xyz)`。顺带把 `time2string` 从 `trialMD_calculator` 的 `contains` 提到文件作用域并**改名为 `crest_time2string`** —— 原名 `time2string` 在 `legacy_algos/trialmd_legacy.f90` 里已经有一个同名的文件作用域版本，直接提上去会**链接期符号重复** |
| 20 | `src/algos/setuptest.f90`<br>`src/dynamics/dynamics_module.f90` | **trial MTD 打印所用参数 + 25/50/75/100% 进度**。① 参数：trial 是以 `pr=.false.` 跑的（`pr` 控制 MD 的全部输出），所以日志里从来看不到时间步、时长、Vbias 参数——而这几个正是判断后面时间预估是否可信的依据。现在在 `Starting trial MTD` 标题后打印 length / time step / steps / temperature / SHAKE / Vbias k,alpha,dump / cores per job。② 进度：`mddata` 新增 `progress` / `progress_steps` / `progress_tag` 三个分量，`dynamics()` 里当 `progress_steps <= 0` 时自动取 `length_steps/4`，在 MD 主循环里**独立于 `pr`** 打印，并带 `flush` —— 不加 flush 的话在重定向到 SLURM 日志时整行都压在缓冲区里，直到任务结束才出现，进度提示就失去意义了。输出形如 `Trial MTD  25%  (   0.06 /   0.25 ps,  elapsed 12 sec, ETA 36 sec)`。**只给 `MDSTART` 打了标记，`env%mddat`（生产 MTD）不受影响**；`progress_steps` 在 `!$omp critical` 内推导，末步单独判断（因为 `length_steps` 一般不是 4 的整数倍） |
| 21 | `src/optimize/ancopt.f90`<br>`rfo.f90` / `gd.f90`<br>`src/algos/dynamics.f90`<br>`parallel.f90` / `search_*.f90`<br>`src/dynamics/dynamics_module.f90` | **轨迹文件后缀 `.trj` → `.xyz`，优化日志 `crestopt.log` → `crestopt.xyz`**。内容格式一点没变（本来就是标准多帧 xyz），只是扩展名改成 `.xyz`，这样 VMD / Avogadro / Mercury / ASE 等工具双击即可打开，不必再手动改后缀。改动清单：① 三个优化引擎里 `open(newunit=ilog,file='crestopt.log')` 改成 `'crestopt.xyz'`；② `MDFILES/crest_<i>.trj` → `crest_<i>.xyz`（`parallel.f90` 里 `trajectoryfile` 的拼装 + `dynamics_module.f90` 里写轨迹头时的 `commentline`）；③ 汇总轨迹 `crest_dynamics.trj` → `crest_dynamics.xyz`（`algos/dynamics.f90` 的 `trjf` 常量与 6 个 `search_*.f90` 里的 `ensnam`）；④ `search_entropy.f90` 的 `crest_<i>.trj` → `.xyz`。**只动现代路径**：`legacy_algos/` 下的老算法仍写 `.trj`，不受影响。注意 `setuptest.f90` 里数优化圈数时读的就是这个文件，已同步改为 `crestopt.xyz` |
| 22 | `src/iomod.F90`<br>`src/crest_main.f90`<br>`src/algos/parallel.f90`<br>`src/algos/setuptest.f90` | **中间文件分级落盘：高频读写进 `CREST_SCRATCH_DIR`，低频读写进 `CREST_WORK_DIR`，用户关心的输出留在当前目录**。上游把所有中间产物（`coord`、`wbo`、`crestopt.xyz`、`gfnff_adjacency`、`gfnff_topo`、`.CHRG`、`.UHF`、`MDFILES/`、`METADYN/`、`NORMMD/`、`MRMSD/`、`OPTIM/`…）一律平铺在工作目录里，既慢又乱。现在按「写入频率」分两级：**热（hot）** = 每几个 MD 步就写一次的目录 `MDFILES` / `METADYN` / `NORMMD` / `MRMSD` / `OPTIM` → 落到 `CREST_SCRATCH_DIR`（节点本地盘）；**冷（cold）** = 一个任务只写几次的文件 `coord` / `wbo` / `crestopt.xyz` / `gfnff_adjacency` / `gfnff_topo` / `.CHRG` / `.UHF` → 落到 `CREST_WORK_DIR`（默认 `crest_work/`，工作目录下的子目录）。`crest_rotamers.xyz`、`crest_conformers.xyz`、`crest_best.xyz`、`crest_ensemble.xyz` 等**结果文件一律不动**，仍在工作目录。实现要点：① `crest_stage_files()` 在 `crest_main.f90` 启动时调用一次，用**符号链接**把上述名字指过去（`ln -s`），而不是改写 CREST 里几十处 `coord`/`wbo` 的调用点——漏掉一个最坏的结果只是那个文件回到工作目录，**运行依旧正确**，属于良性退化；② `inquire()` 判定已存在则**一律不碰**，绝不覆盖已有文件；③ `MDFILES` 与 `TRIALMD` 是**按路径搬迁**的（`crest_hot_dir()`），不能用软链——这两个目录 CREST 会自己 `rmrf()` + `makedir()` 重建，软链会被删掉；且 `iomod%makedir()` 是裸 `mkdir()`，所以父目录必须由 `crest_scratch_dir()` 的 `mkdir -p` 先建好；④ 关闭开关：`CREST_STAGE_FILES=0` 全关，`CREST_WORK_DIR=none` 只关冷文件（恢复上游原布局） |
| 23 | `src/algos/parallel.f90`<br>`src/confparse.f90` | **把上游的死开关 `-iru` 真正接上线：把已找到的构象作为「第二个（静态）偏置」叠到 MTD 上**。上游声明了 `env%iru`（`classes.f90:560`）、给了默认值（`confparse.f90:247`）、也解析了命令行（`confparse.f90:917-918`），但**全代码库再没有任何地方读它**——打开和不开完全一样。现在的语义：`-iru` 打开后，从**第 2 轮 MTD 迭代**起（第 1 轮还没有构象集合），每个 MTD 作业的 `npot` 从 1 变成 2：`mtd(1)` = 原来的动态 RMSD MTD（行为不变，照旧沉积高斯丘）；`mtd(2)` = `cv_rmsd_static`，参考结构取自 `crest_conformers.xyz`（即到目前为止找到的构象）。两者能量/梯度各自计算后相加（`md_calc_mtd` 本来就遍历 `npot`），所以这是**叠加一个记忆项，而不是把动态 MTD 换掉**。实现要点：① 参考集合只在 `crest_search_multimd_init2` 里读**一次**，直接填进每个作业的 `mtd(2)%cvxyz`（**必须换算成 Bohr**，`rdensemble` 给的是 Å）；这样 `mtd_ini()` 走的是「`cvxyz` 已分配」的分支，nsim 个作业不会并发读同一个文件，也不会和本轮后面 cregen 重写 `crest_conformers.xyz` 抢文件；② 静态偏置用**正常的 k/α**（就是该作业自己的 `env%metadfac(i)` / `env%metadexp(i)`，即开头打印的那张 `$metadyn k alpha` 表），不打折；③ `CREST_IRU_KSCALE` 可整体缩放静态偏置的 k（默认 1.0，不用重编即可调小）；④ 因为静态偏置是**对所有已知构象求和**，其量级随集合大小增长，构象多到上百个时可能把结构推散——先小规模试，`CREST_IRU_KSCALE=0.1` 起调；⑤ `clean_V2i` 不删 `crest_conformers.xyz`，所以跨迭代可用；⑥ 每个作业的 banner 会多打两行 `static bias (-iru)` / `static bias k`，日志里 `*sMTD` 前缀只在**只有**静态偏置时出现（本实现里 `cvtype(1)` 仍是 `cv_rmsd`，所以进度行前缀仍是 `*MTD`） |
| 24 | `src/timer.f90` | **修掉计时打印「模 1 小时」的 bug：超过 1 小时的耗时会被印成去掉小时的部分**。`write_timing()` 先把秒数换算成 `walldays` / `wallhours` / `wallmins`，并**依次从 `walltime` 里减掉**天和小时，剩下的 `walltime` 就是「不足 1 分钟的秒数」。但非 verbose 分支（默认，即 `crest` 不加 `-verbose` 时的所有 `... xxx min, yyy sec` 行）只印 `wallmins, walltime`，把已经剥出来的 `walldays` / `wallhours` **直接丢掉**。于是 `1 h 24 min 49.396 s` 印成 `24 min, 49.396 sec`（看起来慢了 3.4 倍），`2 h 14 min 8.978 s` 印成 `14 min, 8.978 sec`（慢了近 10 倍）。外部计算器跑得慢时几乎每轮 MTD 都超过 1 小时，所以这个 bug 会让日志上的耗时全部失真。改法：`walldays > 0 .or. wallhours > 0` 时换成带 `d / h / min / sec` 四段的新格式，否则保持原来的 `min, sec` 格式（短作业的输出排版不变）。verbose 分支本来就印全了，未改动 |
| 25 | `src/calculator/subprocess_types.f90`<br>`src/confparse.f90` | **输入文件回显的版式统一为「`****` 横幅 + 行首 `*`」**。CREST 现在会回显三种输入文件：TOML（`parse_datastruct.f90` 的 `root_print2`）、每个 ORCA 模板（`subprocess_types.f90` 的 `read_orca_input`）、以及 `-cinp` 指定的约束文件（`confparse.f90` 的 `parseRC2`）。上游/先前它们各用一种前缀：TOML 用 `*`，ORCA 模板用 `| `，约束文件用 `> `，一份日志里三种风格并排，很难读。现在统一成 `repeat('*',80)` 夹标题 `* INPUT FILE <名字> content (...)`，正文每行 `* `，再以一行 `****` 收尾——与 TOML 完全同形。注意两处细节：① ORCA 模板的回显必须在 `atmp = adjustl(lowercase(atmp))` **之前**打，否则拿到的是被改过大小写/缩进的版本，不是磁盘上的原样；② `parseRC2` 原来用 `write(*,*)`，现在统一走 `stdout`（`iomod` 里 `stdout => output_unit`，是同一个单元，只是不再依赖 `*` 单元号） |
| 26 | `src/calculator/orca_sc.f90` | **初始几何优化开始前打印 `optimization starting`**。这一步用 `pr=.false.` 跑（`optimize_geometry` 全程不输出），所以用慢的外部计算器时是长达十几分钟到几十分钟的完全静默，看起来像卡死；已有的 `Initial geometry optimization: N cycles, T` 只在**结束**时出现。写在哪：CREST 每次运行的**第一个** ORCA 调用一定来自这一步的第一次梯度计算，所以把标记打在 `ORCA_setup()` 里命令行回显之后（用 `logical,save :: optstartprinted` 保证整个进程只打一次），日志顺序就是 `ORCA call: ...` → 空行 → `optimization starting` →（优化结束）`Initial geometry optimization: 44 cycles, 17 min 32 sec`。即使不设 `CREST_ORCA_PRINT_CALL`，这一行也会出现在「Initial Geometry Optimization」标题和 ORCA 自己的输出之间，正好补上原来那段空白。注意不要挪到 `setuptest.f90` 的 `optimize_geometry` **之前**：那样会排在 `ORCA call:` 前面 |
| 27 | `src/dynamics/dynamics_module.f90`<br>`src/algos/parallel.f90`<br>`src/algos/setuptest.f90` | **MD/MTD 进度只保留 25 / 50 / 75 %，彻底去掉 99 % 和 100 %**。原来的触发条件是 `mod(t,progress_steps) == 0 .and. t < length_steps`，而 `progress_steps = floor(length_steps/4)`：当 `length_steps` 不是 4 的整数倍时会**多触发一次**第 4 个标记，位置在 `t = 4*floor(length_steps/4)`，这个值仍然 `< length_steps`，于是 `pct = 100*4k/(4k+1)` 四舍五入后就是 `99%` / `100%`。日志里于是出现 `*MTD 2 100% ( 3.14 / 3.15 ps, elapsed 19 h 32 min, expected remaining 2 min 14 sec)` 这种既没用（MD 已经结束）又误导（马上又来一行 `*MTD 2 completed successfully ... 34 min`）的行。改法不再依赖 `mod`，而是直接比对三个目标步：`t == progress_steps .or. t == 2*progress_steps .or. t == 3*progress_steps`（外加 `t < length_steps` 兜底），100 % 标记**在构造上就不可能出现**。之前那次「只在 `t < length_steps` 处拦一刀」的修法是不够的，因为被算出 100 % 的那一步本身就 `< length_steps` |
| 28 | `src/classes.f90`<br>`src/parsing/parse_maindata.f90`<br>`src/cregen.f90`<br>`src/algos/refine.f90` | **在 refinement（单点重排序）之前加一道「只做拓扑检查」的预筛，有多级（有 `refine`）时默认打开**。背景：`sort_and_check` → `newcregen` → `cregen_topocheck` 会在 refinement **之后**删掉那些邻接表与 `coord` 不一致的结构，而这些结构的单点能量已经算完了——实测 245 个结构里 136 个（55%）被删，约 54 分钟机时白烧。判据本身与能量无关，所以完全可以提前。改动：① 新增 `env%refine_presort_topo`（`classes.f90`，默认 **`.true.`**），只有 `allocated(env%refine_queue)` 为真时才真正执行，所以「没有 refine 时不会打开」；② TOML 根键 `refine_presort_topo = false` 可关；③ 环境变量 `CREST_REFINE_PRESORT_TOPO=0` / `=1` 覆盖 TOML 与默认值；④ `env%checktopo` 仍然生效，即 `-notopo` / TOML `notopo = true` 会一并关掉它；⑤ 把 `cregen_topocheck` 通过 `cregen_interface` 暴露出来（它原本是 `cregen.f90` 里的外部子程序，签名里的 `env` **没有 intent**，接口必须保持一致）；⑥ 调用点必须在 `xyz = xyz/bohr` **之前**——`cregen_topocheck` 内部自己除以 `bohr`，要的是 `rdensemble` 刚读出来的 Å 坐标。实现要点：`cregen_topocheck` 会把幸存结构搬移到数组前部并对 `xyz` 和 `comments` 施加**同一个置换**，但它不认识实数能量数组，所以把能量写进 `comments`（`f18.8`，和 CREGEN 自己的写法一致），筛完再用 `grepenergy()` 读回来；全军覆没（`nallnew < 1`）时**保留全部结构并告警**，绝不产出空集合。**与现有的 `refine_presort` 的区别**：后者跑的是**完整** cregen（含 EWIN 能量窗 + RMSD 去重），而这两步依赖 refinement **之后**才有的能量，默认开会有误砍风险；本开关只做拓扑，科学结果与不筛时**完全等价**，纯粹省机时 |
| 29 | `src/dynamics/shake_module.f90`<br>`src/classes.f90`<br>`src/cregen.f90`<br>`src/confparse.f90`<br>`src/parsing/parse_maindata.f90` | **拓扑判据优化：默认把「含金属原子的键」排除在拓扑比较之外，判据与 SHAKE 完全一致**。背景：`bondtotopo`（`ztopology.f90:270`）用**整数化配位数**建邻接表（`icn = nint(cn2(i))`，再只保留最强的 icn 条键）。金属中心的 CN 大、软、随几何漂移，取整时一跳就多/少一条键 → `cregen_topocheck` 判为"拓扑变了" → 整片构象被误删（实测 245 个删 136，55%）。SHAKE 早就因为同一个理由拒绝约束金属键（`init_shake`：`metalbond = metal(at(i))==1 .or. metal(at(j))==1`，**唯一例外是 Li/Be**，它们的键照样约束）。现在把这张表（原本是 `shake_module` 的 private）以 `crest_is_metal(z)` 导出，CREGEN 复用**同一个定义**：① 新开关 `env%excludeTOPOmetals`（`classes.f90`，默认 **`.true.`**）；② `cregen_topocheck` 把 `-notopo <list>` 和金属原子 OR 成一个 `topoexcl(nat)` 掩码，然后**统一走 `bondtotopo_excl`**——该函数只在 `bondtotopo` 基础上多一行 `if (excl(i).or.excl(j)) neighbourmat(i,j) = .false.`，且 `excl(i).or.excl(j)` 使**双向**都置零，所以金属自己的键、以及后面那段"碳特例"都会被一起屏蔽；掩码全 false 时它与 `bondtotopo` **逐字等价**，因此无金属体系的行为完全不变；③ 关掉：`-nometaltopo` / TOML `metaltopo = false`（`nometaltopo = true` 亦可）。**关键坑**：`topoexcl` 必须在构建**参考**拓扑之前就填好——参考结构是后面所有比较的基准，先把 `bondtotopo_excl` 调用写在掩码分配之前会导致未分配就使用。日志会多两行 `# atoms excluded from topo` / `=> of which metals`（写进 `cregen_<n>.out`，并在屏幕上额外打一行 `CREGEN> N metal atom(s) excluded from the topology check`） |

| 30 | `src/classes.f90`<br>`src/confparse.f90`<br>`src/parsing/parse_maindata.f90`<br>`src/printouts.f90` | **新增 `--onebatch` 模式：只跑一批 MTD，然后走常规 MD，不再有第二轮及以后的迭代 MTD，能量窗口放宽到 8 kcal/mol**。定位与 `--quick` / `--squick` / `-mquick` / `--extensive` 同类（都在 `confparse.f90` 的 V2 分支里），但**不改 `runver`**，因此沿用默认的 iMTD-GC MTD 参数网格（12+2 条），可以与 `--quick` 叠加。三个效果各由一个**已有**开关承担，没有新增控制流：① `env%Maxrestart = 1` —— `mtdloop: do i = 1, env%Maxrestart` 只进 i=1，而 `i==1 .and. start` 分支直接 `cycle` 出去，**只执行一轮 MTD**；② `env%allowrestart = .false.` —— `elowcheck`（`confscript2_misc.f90:1127-1130`）会把 `lower` 强制置回 `.false.`，于是 `mtdloop` 的 `exit` 与 `MAINLOOP` 里 rotamermds / GC 之后的 `cycle MAINLOOP` **都不会再触发**，MD（默认开）照旧执行，之后直接进入最终优化；③ `env%ewin = 8.0d0`（与 `--extensive` 同值）。**绝对不要顺手把 `env%iterativeV2` 也关掉**：它在 `choose_settings.f90:293` 被读，关掉会把 12+2 的 k/α 网格换成 24 条并改掉两个增量，是一次与意图无关的采样变更；它剩下的唯一作用只是 GC 块里 `cycle MAINLOOP` 的一道闸（`search_conformers.f90:226`），而那道闸已经被 `allowrestart=.false.` 覆盖了。命令行 `--onebatch` / `-onebatch` 均可（`--` 前缀在 arg 循环里被剥成 `-`）；TOML 根键 `onebatch = true` 同效。**顺序陷阱**：与 `-quick`(→5.0) / `-squick`(→5.0) / `-mquick`(→2.5) / `-extensive`(→8.0) 一样，`ewin` 是在解析到该参数时立即写入的，所以 `--onebatch` 要写在其它模式开关**之后**才会生效；反过来想让 `-ewin 12` 压过它，就要把 `-ewin` 写在更后面。该模式不会动 `rotamermds`（默认 `.true.`），所以 `-norotmd` 仍可单独关掉 MD。日志的 Job settings 段会多一行 `one-batch mode : "-onebatch" (1 MTD round, no iteration)` |
| 31 | `src/dynamics/dynamics_module.f90`<br>`src/sigterm.f90` | **MD 轨迹与 stdout 的「写了但没写完」问题：在每一帧/每一次打印后显式 `flush`**。症状：运行时去看 `crest.out` 或 `TRIALMD/crest_trial_md_N.xyz`，内容**停在某一帧的中间**（例如 97 个原子只写出前 59 个），看起来像"输出被截断"。根因是 Fortran 运行时 I/O 缓冲，不是 CREST 写错了：轨迹文件在 `dynamics()` 开头 `open(newunit=trj,file=trajectory)`（:308）、直到 MD 结束才 `close(trj)`（:514），中间**一次 `flush` 都没有**；写进去的字节全部堆在运行时缓冲区里，只有缓冲区凑满或 `close` 时才真正落盘，而缓冲区边界几乎必然落在某一帧的中间。同理，`stdout`（`output_unit`，被重定向到文件时是**全缓冲**）也会一块一块地出现在 `crest.out` 里。改动：① `wrxyz(trj,...)` 之后立刻 `flush (trj)` —— 每帧一次 `write(2)`，相对产生这一帧的 engrad 调用可以忽略；② 同一处 `flush (stdout)` —— 因为 `mddata%progress` 默认 `.false.`，并行 MTD（`pr=.false.`）跑起来**完全没有心跳**，unit 6 会几个小时不动，而缓冲是 FIFO 的，所以这一下会把之前所有积压的日志一起推出去；③ MD 表格打印后 `flush (stdout)`；④ `creststop()` 在 `call exit(io)` 之前 `flush (stdout)` —— `exit(3)` 不保证冲刷 Fortran 单元缓冲，否则"CREST terminated ..." 这行自己就可能丢失。**关键事实**：Intel 的 `FORT_BUFFERED` 环境变量**明确排除 stdout（`*` 与 unit 6）**，所以 `crest.out` 的问题只能靠显式 `flush` 解决，设 `FORT_BUFFERED=false` 只能影响文件单元；文件单元在 Intel 下默认是**不缓冲**的，因此如果看到 xyz 也被截断，先 `env \| grep FORT_` 确认运行环境/登录 shell 是否把它设成了 `true`（很多超算默认设），或二进制是否为 gfortran 构建。副作用：被 Ctrl+C / walltime 杀掉时最多丢一帧，而不是丢整个缓冲区 |
| 32 | `src/classes.f90`<br>`src/parsing/parse_xtbinput.f90`<br>`src/dynamics/shake_module.f90`<br>`src/dynamics/dynamics_module.f90`<br>`src/choose_settings.f90`<br>`src/algos/parallel.f90`<br>`src/algos/setuptest.f90`<br>`src/confparse.f90`<br>`src/parsing/parse_maindata.f90`<br>`src/printouts.f90` | **`--cinp` 读进来的键约束默认改为刚性（SHAKE）约束，不再是谐振子软约束**。背景：软约束是 `V(r)=½k(r−r_ref)²`，k 以 **Eh/Bohr²** 计（`calculator.F90:257-289` → `calc_constraint` → `constraints.f90:545 bond_constraint`，注释明写"r and k are in atomic units"）；`--cinp` 文件里的 `force constant=` 默认 0.05、CREST 内部 `fcdefault=0.01`。k=3.5 是默认值的 70–350 倍、比真实 C–C 键（≈4.5 mdyn/Å）硬 12 倍 → 周期 `T=2π√(μ/k)` 对 C–C 只有 **8.5 fs**、C–H **3.3 fs**，而 MD 默认步长 **5 fs**（`confparse.f90:218`），稳定积分要求 `dt ≲ T/10` ⇒ 轨迹发散，唯一可见的症状就是 `SHAKE did not converge` / `Trial MTD did not converge`。刚性约束没有力常数、没有频率、不限制步长，这才是"固定过渡态反应键"真正要的东西。实现链路（**只在解析约束文件时登记，绝不用 `calc%cons` 反推**——`-cbonds` 经 `autoconstraint_internal` 也会为每根键生成 type=1 约束，全变刚性会把整个分子冻住）：① `parse_xtbinput.f90` 里每读到一个 `distance`/`bond` 就 `call add_hardcons(env,i1,i2,dist)`（`dist` 已由 `aatoau` 转成 Bohr），存进 `env%nhardcons/hcatms(2,n)/hcdist(n)`，按 16 个一档增长；② `choose_settings.f90` 的 `env_to_mddat()` 在原有 SHAKE 段之前调 `crest_hard_bonds_to_shake(env%mddat%shk,env%nhardcons,env%hcatms,env%hcdist,.true.)` 填 `shk%nusr/conslistu/distconsu`；`env%mddat` 是模板，`setuptest.f90:121 MDSTART=env%mddat`、`parallel.f90 mddat=env%mddat`、`dynamics.f90 mddat=env%mddat` 都会带走它；③ `shake_module.f90` 里 `distcons` 默认取**传入几何的当前距离平方**（:287），现在若 `distconsu` 存在就改用它的平方 —— 否则用户写的约束值会被静默忽略、键被钉在起始结构碰巧的那个长度上；同时合并用户键与自动探测键时**去重**（`shk%ncons = shk%nusr+nadd`，默认 `shake_mode=2` 时反应键本来就在自动列表里）；④ 刚性约束**必须**让 SHAKE 开着：`choose_settings.f90` 里 `if (env%mddat%shk%nusr > 0) env%mddat%shake = .true.`（`-shake 0` 也压不过），`setuptest.f90` 里 `MD%shake = (shakemode > 0).or.hardcons` 且**不允许**自动重试链把 SHAKE 关掉（`if (tstep <= 1.0d0.and.hardcons) … status_safety`，否则原始故障会让它降级成"关 SHAKE"从而丢掉约束）；`shake_mode==0` 时 `init_shake` 不分配 `cons2` ⇒ 结果恰好是"只有用户那几根键，别的都不约束"，这正是 `-shake 0 --cinp` 想要的语义；⑤ `parallel.f90` 里 `mddat%shk = shk` 会**覆盖**掉 `mddat=env%mddat` 带来的那份 shk，所以必须显式把 `nusr/conslistu/distconsu` 搬过去；⑥ `env%mddat%shk%shake_mode` 那行**保持原样**（不能因为 nusr>0 就跳过），否则单 MD 路径会退化成"只有用户键"。**软约束刻意保留**：SHAKE 把 r 钉在 r_ref 上，`(r−r_ref)` 在 tolshake 内为 0 ⇒ MD 期间弹簧贡献为零；而几何优化与 CREGEN **不跑 SHAKE**，那时它是唯一拉住这根键的东西。开关：默认硬约束，`-softcons` / `-soft` 或 TOML `hardconstraints = false` 回到旧的软约束行为，`-hardcons` / `-hard` 显式打开。**已知边界**：`legacy = true`（TOML）时 `confparse.f90:2336` 会跳过 `parse_xtbinputfile`，因此登记不到，此时该功能静默不生效；TOML `[calculation]` 的约束块（`parse_calcdata.f90:766 rdbondconstraint`）同理未接入 |
| 33 | `src/dynamics/dynamics_module.f90` | **编译修复（Round H 的收尾）：`allocate` 里整型数组与实型数组共用同一个 `source=`**。症状：`[ 83%] Building Fortran object … src/dynamics/dynamics_module.f90.o` 之后报「`Error: Type of entity at (1) is type incompatible with source-expr at (2)`」，停在`allocate (shk%conslistu(2,n),shk%distconsu(n),source=0)`（:1377）。根因：`allocate` 的`source=` / `mold=` 是**对每一个分配对象都要成立**的一致性约束，不是只管第一个；`shk%conslistu` 是 `integer`、`shk%distconsu` 是 `real(wp)`，而 `0` 是整型字面量，于是第二个对象的类型不匹配。**fparser 只做语法分析、不做类型检查，所以本地那道`f90_syntax_check.py` 关卡完全看不到它**，只有真正调用 gfortran 才会炸。修法是每种类型各一条 `allocate`：`allocate (shk%conslistu(2,n),source=0)` + `allocate (shk%distconsu(n),source=0.0_wp)`（注意实型那条的 `0.0_wp`，写成 `0` 又会是整型）。防复发：`/PATH/TO/crest-code/_orig/f90_semantic_check.py` 专门查这一类，运行时请传**目录**而不是文件清单——派生类型表必须在审计任何文件之前先建全，否则跨文件的 `type(T) :: x` 组件查不到。它已被并进 `verify-fixes.sh` 的 C 节（`D. dependencies` 顺延为 D 节），全树 164 个文件当前`0 ERROR, 0 UNRESOLVED`；`verify-fixes.sh` 的 B 节也加了两条断言专门盯这两行，防止再被手贱合并回一条 |
| 34 | `src/parsing/parse_xtbinput.f90`<br>`src/dynamics/dynamics_module.f90` | **补上 `atoms:` / `elements:` 这条路径：pairwise 距离约束现在也登记为刚性（SHAKE）**。发现经过：用户的约束文件写的是 `$constrain` + `atoms: 1,2,6,7,97` + `force constant=3.5`，而不是 `distance:`；而 `add_hardcons` 原先**只在 `case ('distance','bond')` 分支里被调用过一次**（`:193`）。`atoms:`/`elements:` 走的是**完全独立的 pairwise 分支**——先在 `:270-291` 标记 `pairwise(:)` 掩码，真正的约束要到 `:300-329` 才生成（`conslist(k)%bondconstraint(j,i,dist,force_constant)` 再 `env%calc%add(k,conslist)`），全程没有登记 ⇒ 这类文件 100% 保持谐振软约束，第 32 项的"默认硬约束"对它**完全没生效**，块里的 `force constant=3.5` 反而仍是唯一决定刚度的参数。改动就是在 `:320` 之后补一句 `call add_hardcons(env,j,i,dist)`。**两个必须知道的事实**：① pairwise 生成的是**全部两两组合** C(n,2)，写 5 个原子 = **10 根**约束而不是 5 根；② 距离约束的振动频率是 `sqrt(k/mu)`，**与两点相距多远无关**——所以即使是很远的非成键 pair，在 k=3.5 下也和真实化学键一样硬（μ=6 amu → T≈8.5 fs，要求 dt ≲ 0.85 fs；含氢 μ≈0.9 → T≈3 fs）。因此 `--tstep 3.0` 配 k=3.5 之前必然发散，这很可能正是 trial MTD 不收敛的直接原因。冗余风险：5 个原子全互连给出 10 个距离，而该刚体只有 3n−6 = 9 个内部自由度，约束集**相容但线性相关**，SHAKE 必须做投影、收敛变慢 ⇒ `crest_hard_bonds_to_shake` 现在会在开头打印一段 `note: these bonds tie 5 atoms together with 10 distances ... REDUNDANT` 的诊断，并建议收窄原子列表或退回 `-softcons`。另注：`$metadyn` 块里的 `atoms:` 是另一回事（选 MTDRMSD 的集体变量原子，`get_xtb_metadyn_block` 的 `:569`），**不是约束** |

| 35 | `src/dynamics/shake_module.f90`<br>`src/dynamics/dynamics_module.f90` | **MD 失败现在会说清是"哪里"失败 + 是哪根键（Round P 的 P1+P2，纯诊断，不改任何数值行为）**。背景：`dynamics()` 只有两个提前退出口，却对外表现为**同一个整数**——`engrad()` 失败（`:340`，外部 xtb/ORCA 子进程出错）与 `do_shake()` 返回 `io/=0`（`:504`，约束集解不开）；调用方 `setuptest.f90:229` 只判 `if (io == 0)` ⇒ 两种情况在日志里都只剩一句 `Trial MTD n did not converge!`。更糟的是 SHAKE 那句诊断原文`else if (pr) then write (*,*) 'SHAKE did not converge! maxdev=',maxdev`（`shake_module.f90:479`），而 trial MTD 是**用 `pr=.false.` 跑的**（`setuptest.f90:164`）⇒ 最关键的信息被静默吞掉。**判断依据**：4 / 3 / 2 / 1 fs 四档全挂 ⇒ 与步长无关，不可能是刚度问题（刚度问题一定会随步长减小而好转），只可能是"约束集不可行"或"外部计算器报错"。改动：① 去掉 `pr` 门控，无条件调用新增的 `shake_failreport`，打印**迭代数/上限**（默认 `maxcyc=250`）、**最大偏差/容差**（默认 `tolshake=1e-7`，注意是**距离平方**的相对偏差，即约 1e-7 Å）、**违反约束的条数**，以及**最差的 8 根键的原目标与实际距离（Å）**——这列数字是区分"不可行"（r 根本到不了 r₀，通常是 `reference=coord.ref` 与 coord 不一致）与"冗余"（能到，但求解器在 250 次 sweep 内到不了）的唯一凭据；② `dynamics()` 两个出口各调一次新增的 internal procedure `md_failreport`，打印**失败来源**、**发生在第几步 / 共几步**、**term 状态**，并且 `term==1`（第一步就挂）时直接提示"去看 coord 与 coord.ref 是否一致"。**为什么可以无条件打印**：两条路径都会立刻 `exit MD`，所以一次 MD 最多输出一份报告，不会把日志刷爆；写操作套了 `!$omp critical` 并 `flush`，避免并行 MTD 的输出交错串行、以及重定向日志滞留缓冲区。副作用为零——没有改任何阈值、步长或约束列表 |

| 36 | `src/dynamics/shake_module.f90` | **SHAKE 失败报告现在可读、且能自我判别（Round Q，仍是纯诊断）**。第 35 项的报告第一次真机触发就暴露了两个自己身上的 bug：① 格式用 `e12.4` 打印时**把指数符 E 吃掉了**，于是 `0.4404E+291` 变成 `0.4404+291`，看着像 0.44；`f9.5` 更直接溢出成 `*********`（这正是日志里 `r=*********` 的来源）。**全部数值格式改为 `es12.4`**（ scientific with explicit E，`es` 保证指数位永远带符号和 E）。② 报告只给了**最后一次**迭代的偏差 ⇒ 无法区分"目标不可行"（残差是**有限**的，量级≈目标错配）与"迭代发散"（从 1e-9 放大到 1e290）。修法：新增 `real(wp) :: maxdev0`，在 `icyc == 0` 时捕获——注意 `icyc` 是在循环体**末尾**才自增的，所以第一次进入体时它仍是 0，捕获到的正是**未经任何修正的原始偏差**；报告里同时打印 `deviation BEFORE correcting` 与 `AFTER correcting`。再加一段自动判别：`maxdev0 <= 1e-3 且 maxdev > 1e3` ⇒ 打印 `ITERATION BLEW IT UP`（求解器不稳，看第 37 项）；`maxdev0 > 1e-3` ⇒ 目标不可行/步长过大（去看 `reference=` 与 coord 是否一致）；否则 ⇒ 冗余约束导致的停滞。放大倍数 `amp` 计算加了 Inf/NaN 保护，避免溢出后再除以 0 打出 `NaN`。**真机结论**：一次 97 原子体系的报告是 `constraints violated: 25 of 110`、`BEFORE=1e-9`、`AFTER=0.4404E+291` ⇒ **不是冗余、是发散**；而且用日志里 `Input structure:` 的坐标反算，发现 `bond i-j fixed at` 打印的目标与 coord 相差最大 5.2e-3 Å ⇒ `reference=coord.ref` 生效且 coord.ref ≠ coord，110 根约束里 100 根自动键按 coord 钉、10 根用户键按 coord.ref 钉，5e-3 Å 的缝对 1e-7 的容差是无法弥合的 |

| 37 | `src/dynamics/shake_module.f90` | **SHAKE 修正扫描由 Jacobi 改为 Gauss-Seidel（Round R 的 P4，第一个真正改数值的修复）**。根因：`do_shake` 的迭代体在循环**顶部**一次性算好 `shk%dr(1:3,i)` / `shk%dr(4,i)`（:413-414），下面的修正 sweep 却要用这同一批**陈旧**的键向量去处理全部 ncons 根约束——而 sweep 走到第 i 根时，前面 i−1 根已经把原子挪过了。也就是说整个 sweep 在**同一个、已经过期的几何**上做线性化，这是 Jacobi 迭代；强耦合约束集（110 根共享原子，外加 `--cinp` 的 10 根）上它的迭代矩阵谱半径远大于串行版本，后果**不只是慢，而是放大**：真机上一步之内偏差从 1e-9 涨到 1e290。修法：在**两条**修正分支（非冻结分支 ~:470、冻结分支 ~:489）里，紧贴 `denom` 计算之前各补两行，用**当前**的 `shk%xyzt` 重算 `shk%dr(1:3,i)` 与 `shk%dr(4,i)`，这就是 GROMACS / AMBER / LAMMPS 用的经典串行 SHAKE。**刻意不改的三件事**：① 修正**方向**仍是 `dro`（步进前的键向量），经典 SHAKE 本来就固定这个参考，只刷新线性化用的 `dr`；② 质量权重 `rmi/rmj` 不动，所以投影仍是最小质量加权位移、拉格朗日乘子不变；③ 后面的速度/加速度修正用的是**实际施加的位移** `shk%xyzt - xyz`，自动跟随。因此迭代的**不动点**（"所有约束都满足"）与原来**完全相同的方程组**——两种算法都收敛时结果在 tolshake 内一致，Gauss-Seidel 只是收敛得更快、更稳，**不改变采样系综**；只有在需要 bit-exact 复现旧轨迹时才算"变化"。`verify-fixes.sh` 除了两条逐分支断言，还加了一条计数检查：`grep -cF "shk%dr(1:3,i) = shk%xyzt(1:3,iat)-shk%xyzt(1:3,jat)"` 必须 ≥ 3（测试循环 1 处 + 两个分支各 1 处），防止以后只改了一支就以为改完了。**重要边界：P4 治不了"目标不可行"**——第 36 项那个 5.2e-3 Å 的目标错配属于约束集本身无解，任何迭代格式都救不回来。正确的处理顺序是**先去掉 `reference=coord.ref`（零改码）让 coord 与 coord.ref 一致**，若仍然出现 `BEFORE` 很小而 `AFTER` 爆炸，才是 P4 该上场的时候 |

| 38 | `src/dynamics/shake_module.f90` | **SHAKE 扫掠预算 `maxcyc` 由 250 提到 1000（Round S 的 P3-a，第二个改数值的修复）**。修完 P4 后真机日志变成了：同一体系、同样的 10 根键、`constraints violated : 3 of 110`、`deviation BEFORE correcting : 2.5213E-01`、`deviation AFTER correcting : 1.2887E-07`、`tolerance : 1.0000E-07`。**1e290 消失了** —— P4 是对的。剩下的是纯预算问题：把端点拟合一下，`rho = (1.2887E-07/0.2521)^(1/250) = 0.9437`，即 `dev(n) = 0.2521 * 0.9437**n`，解 `dev = 1e-07` 得 **n ≈ 254**；校验 `0.2521*0.9437**250 = 1.288e-07` 与日志**逐位吻合**。`maxcyc` 是 `shake_module.f90:88` 的类型默认值（`integer :: maxcyc = 250`，与 `tolshake = 1.d-7` 并列，**全树没有任何 CLI/TOML 旋钮，也没有任何赋值点**——`%maxcyc` 只在 `:507` 被读、`:641` 被打印）。⇒ 它既不是不可行（那会停在有限残差）、不是发散（那是 1e290）、也不是冗余停滞（残差一直在动），**就是差了 4 次扫掠**。改为 1000 的理由是余量：收敛是几何级数，多 17.3 次扫掠多换一个数量级，1000 次能把同样起点推到 ~1e-19，再也不会成为约束。**语义上零变更**：这个改动只可能把"放弃"变成"收敛"，不可能改变任何在 250 次内已收敛的步；额外代价只在**本来会直接中止整个 run** 的那些步上发生，且是 ncons 条约束上的纯算术，相对一次 engrad 调用可忽略。**同时要纠正 Round Q 的一个过度结论**：当时说 `reference=coord.ref` 的 5.2e-3 Å 错配"永远不可满足"，过强了 —— SHAKE 本来就能吸收大的初始违背（`maxdev0 = 0.25` 正是它的量级），而扫掠需求对 `maxdev0` 只呈**对数**依赖；本次日志里该错配**依然存在**（1-2：target `1.46187` vs coord 算出 `1.46637`，差 4.5e-3 Å），SHAKE 照样收敛。**顺带解释"为什么降步长的阶梯只是勉强过"**：`n ≈ ln(maxdev0/tol)/(-ln rho)`，ρ=0.9437 ⇒ 每十倍一个数量级只需 17.3 次扫掠；`maxdev0 ∝ dt²` 时 4→3 fs 只省 `ln(1.78)/0.05795 ≈ 9` 次（254→245，**2% 余量**），4→1 fs 约省 48 次。所以阶梯看起来"能过但总在临界的另一侧"，真正杠杆是 `maxcyc`。**一个实测警告**：不要图省事"少写一条约束"——实测删掉最长的那条 `1-7`（3.364 Å）会让约束集**失稳**，即使已是 Gauss-Seidel，dev 仍冲到 `9.3e+88`；删 `2-6` 才稳。复现脚本（纯 stdlib，`/PATH/TO/crest-code/_orig/`）：`shake_conv_test.py`（冗余 vs 去冗余 ρ 对比）、`shake_sweeps_needed.py`（达标扫掠数标定）、`shake_floor_test.py`（残差地板来源判别：10 条距离对 5 个原子是**超定**，任何舍入都留下永久地板；9 条则是**方阵**，泛型可精确实现） |

| 39 | `src/dynamics/shake_module.f90` | **修正 SHAKE 失败报告的 DIAGNOSIS 误判，并增印 ρ 与"还差几次扫掠"（Round S 的 P3-c，仍不改任何数值行为）**。问题：原来的三分支是按 `maxdev0` 判的，于是真机那次（`maxdev0=2.5213E-01`、`maxdev=1.2887E-07`、`tol=1.0000E-07`）落进了 `maxdev0 > 1e-3` 那一支，打出 **"the step is too large - here a shorter time step DOES help"** —— 完全反了：残差已经降了**六个数量级**、离容差只差 1.29 倍，它失败是因为**扫掠次数用完了**，而缩短时间步对扫掠需求只有**对数**效应（4→3 fs 只省约 9 次）。改法：① 新增拟合 —— `rho = (maxdev/maxdev0)**(1/icyc)`（实测 0.9437），`n_more = ln(tolshake/maxdev)/ln(rho)`（实测 **4.4**），并在报告里以 `contraction per sweep (rho)` / `more sweeps needed to reach` 两行打印，日志因此**自证**是"几何收敛但预算不够"，而不是发散或停滞；② 判据改为按 `maxdev` 与 `tolshake` 的关系分四支：**爆炸**（`maxdev > 1e3` 且 **`maxdev >= max(10*maxdev0,1e-3)`**，即真的**长大了** —— 旧写法要求 `maxdev0 <= 1e-3`，会把"起点就很差 + 爆炸"误判成停滞）、**预算不足**（`maxdev <= 1e2*tol`，明说 NOT a feasibility problem 并给出还差多少次扫掠）、**停滞**（`rho >= 1`，残差不再下降 ⇒ 不可行/冗余，提示去读 `r -> r0` 那一列）、**仍在下降但离容差很远**（优先提 `maxcyc`，其次才考虑步长）。③ 全部对 `maxdev0 <= 0`、非有限值、`rho >= 1` 做了保护。**验证**：`/PATH/TO/crest-code/_orig/p3c_branch_test.py` 把 6 个真实/构造场景回放一遍，全部命中预期分支（真机那次 → 分支 2，ρ=0.9437、n_more=4.4，与日志反演逐位一致）；门禁新增 6 条存在断言 + 1 条**否定断言**（`time step DOES help.` 必须消失），`verify-fixes.sh` 178 OK / 0 FAIL。**踩坑记录**：改这个断言数组时，一次编辑把 7 行的**结尾双引号全部吞掉**，`bash -n` 报的是**下面几行**的 `syntax error near unexpected token |` —— 数组字面量里，报错行往上看才是真凶，改完一定要先 `bash -n` |

| 40 | `src/algos/setuptest.f90`<br>`src/algos/parallel.f90` | **trial MTD 的默认长度按计算器区分：调用 ORCA 时降到 0.125 ps，其余保持 0.25 ps**。背景：trial MTD 只干两件事 —— (a) 验证所选时间步/SHAKE 设置下 MD 能否跑完，(b) 提供一个 wall time 样本供外推到 production MTD；对内部 GFN 方法它只花几秒，但对 ORCA 这种**完整 SCF**（一次梯度几十秒到几分钟）它纯粹是开销。位置：`setuptest.f90:111` 附近的 `trial_md_ps = 0.25_wp`（上游是 1.0 ps，本树此前已砍到 1/4，并留有 `CREST_TRIAL_MD_PS` 环境变量覆盖）。改法：① 新增谓词 **`logical function crest_calc_is_orca(calc)`**（放在 `parallel.f90` 里 `crest_calc_is_subprocess` 的正下方，风格一致；`setuptest.f90` 用 `logical,external` 声明链接，与它的兄弟完全同款做法）。判据是遍历 `calc%calcs(1:size(calc%calcs))%id == jobtype%orca` 的 **any** 语义而不是只看第一个 —— ONIOM 这类复合计算器可以把 ORCA 放在非首位，而成本仍由 ORCA 子进程主导；开头有 `if (.not.allocated(calc%calcs)) return` 防止计算器未初始化时崩掉；用 `size()` 而非 `ncalculations`，谓词不该依赖计数器与数组同步。② `trial_md_ps` 变成 `if (crest_calc_is_orca(env%calc)) then 0.125_wp else 0.25_wp`，**环境变量覆盖仍然压过两者**。③ 新增 `character(len=64) :: trial_src` 并在 settings 横幅多打一行 `length set by      :`，值为 `ORCA default (0.125 ps)` / `internal default (0.25 ps)` / `CREST_TRIAL_MD_PS env override`，免得日志里出现 0.125 被当成默认值没生效。**随之变化的数字**：默认 4 fs 下 trial 步数 63 → **31**（3 fs: 84 → 42，2 fs: 126 → 63，1 fs: 252 → 125）；`MTD%cvdump_fs` 当时还是 `550*trial_md_ps`（137.5 → 68.75 fs）—— **这一项已在第 41 行被改为固定的 50 fs，不再随长度缩放**。运行时估计是 `rtime*(mdtime/trialps)`，按每 ps 归一，所以估计**依然成立**，只是样本变短、估计更抖。**取舍（必须知道）**：trial 更短 ⇒ wall time 样本更少 ⇒ ETA 更抖；且只在约 50 步之后才显现的不稳定性会漏检。0.125 ps/4 fs = 31 步对"第 1 步就挂"这类问题绰绰有余（本次真机故障就死在 step 1），但它是**标定而非证明** —— 若哪天遇到 trial 通过而 production MTD 炸掉，应该把它调回去，而不是继续缩短。**不受影响**：`legacy = true` 走 `trialMD_legacy`（`legacy_algos/trialmd_legacy.f90:57` 硬编码 `length = 1.0d0`），那是老的 xtb-only 子进程路径，会写 xtb 的 MD 输入、**无法驱动 ORCA**，所以 ORCA 流程必然走 `trialMD_calculator`，本次改动覆盖完整 |

| 41 | `src/algos/setuptest.f90` | **trial MTD 的 Vbias dump 默认改为固定 50 fs（不再随 trial 长度缩放）**。原状：`MTD%cvdump_fs = 550.0_wp*trial_md_ps` —— 550 fs 是上游为它当年那个 1 ps trial 定的，本树把 trial 缩短后为"避免 0 次 bias 更新"而跟着缩放。**缩放是错的杠杆**：这个间隔回答的是"偏置多久自适应一次"，答案不该取决于 trial 恰好有多长。改成固定 **50 fs**（`trial_vbdump_fs = 50.0_wp`，可用 `CREST_TRIAL_VBDUMP_FS` 覆盖，与 `CREST_TRIAL_MD_PS` 对称）。**实测对比（trial 内 bias 更新次数）**：ORCA 默认 0.125 ps 时旧 68.75 fs → **1 次**，新 50 fs → **2 次**；内部默认 0.25 ps 时旧 137.5 fs → **1 次**，新 50 fs → **5 次** —— 也就是说旧写法下 trial 其实**几乎没在测试偏置**，这本身就是个隐性缺陷。**保留原有的安全性、但不用缩放来实现**：`MTD%cvdump_fs = min(trial_vbdump_fs,trial_md_ps*500.0_wp)` —— 间隔不得超 trial 的一半，于是**任何** trial 都至少有 2 次 bias 更新。这条守卫是必需的，因为 `CREST_TRIAL_MD_PS` 可以把 trial 缩到任意短；没有它，`CREST_TRIAL_MD_PS=0.02` 会让 trial 退化成纯 MD 却仍被当成 metadynamics 测试。已验证 `0.5 / 0.25 / 0.125 / 0.05 / 0.02` ps 五档全部 ≥2 次。**注意单位换算的下游行为**：fs 值会转成整数步（`metadynamics_module.f90:126`，`dum2 = max(1, cvdump_fs/tstep)`），所以 4 fs 下实际间隔是 13 步 = 52 fs，**横幅打印的是请求值 50 fs**。**production 侧不受影响**：那走 `env%mddump = 500` fs（`confparse.f90:232`，可由 `-vbdump <ps>` 覆盖），与本行无关。**门禁教训**：本次改动让一条旧断言（钉 `550.0_wp*trial_md_ps` 这个缩放公式）变红 —— 校验脚本按预期工作，已把它改写成断言新的 `min(...)` 守卫 |

| 42 | `src/choose_settings.f90`<br>`src/classes.f90`<br>`src/printouts.f90` | **`--onebatch` 改为固定跑 2×2 = 4 条并行 MTD（比 `-quick` 的 6 条再少一级），k = 0.0020 / 0.0010、α = 1.00 / 0.50**。背景：`-onebatch` 此前只靠两个既有杠杆实现（`Maxrestart = 1` + `allowrestart = .false.` + `ewin = 8`），**完全没有碰 k/α 栅格**，所以它实际跑的是默认 V2 的 12+2 = 14 条 —— 单轮 + 14 条，对 ORCA 这种按次收费的外部计算器恰好是最贵的组合。位置与做法：`defaultGF()` 里 `select case (env%runver)` **结束之后**、栅格生成循环（`m = 0` 之前）插入一段覆盖块，**不新增 controlflow、不新增 runver**（沿用本扩展一贯的做法）：`if (env%onebatch.and.env%runver.ne.787878) then na=2; nk=2; nrem=0; nmtdyn=na*nk; alp=1.0d0; kstart=0.002d0; alpinc=2.0d0; kinc=2.0d0`。**`nrem = 0` 不可省**：`-nmtd N` 会经 `gfdistribute` 留下非零 `nrem`，而下面 `nmtdynmax = nmtdyn-nrem`，不清零会让 4−nrem 掉到 1 以下。**选点的四条理由**：① **拓扑用 2×2 而不是 4×1** —— 4×1（一个 k、四个 α）会把四条轨迹压在同一推力上，这个 k 一旦选错四条全废；2×2 在两个维度上同时对冲。② **k 取 `-quick` 那一对** 0.0020/0.0010（`rednat ≈ 92` ⇒ `metadfac = 0.184 / 0.092 Eh`），是 CREST 在"故意少跑 MTD"时唯一用过的 k 区间；它比默认 V2 的上限 0.003 温和 1.5 倍、比 `-nmtd` 的上限 0.005 温和 **2.5 倍** —— 这一点对本体系要紧：推力越大、每步位移越大、SHAKE 的**初始违背 `maxdev0` 越大**（`maxdev0 ∝ 力 × dt²`），而本体系有 10 根 `--cinp` 刚性键，第 38 项实测 `maxdev0 = 0.2521`、需要 ~254 次扫掠才到 1e-7。③ **α 取 1.00 / 0.50** —— 偏置高斯半宽 `sqrt(ln2/α)` = **0.83 / 1.18 Å**。α ≤ 0.3（半宽 1.5 Å）对 90 原子体系已经分不清构象；α ≥ 1.3（默认 V2 的上限）比"单轮不重复"需要的更激进。两者几何中心 0.707，落在默认 V2（`0.281…1.3`）与 `-quick`（`0.3…1.2`）**共同认可区间**的偏上位置。④ **刻意不动 `rfac`** —— `runver` 保持 1，单条长度仍是 10 ps（基数 `anint(9.9)`），Σ = 4×10 = **40 ps**，比 `-quick` 的 30 ps 还多；单轮没有第二轮来补偿，靠**更长的单条**补，而不是靠更宽的栅格。**`-nmtd <N>` 仍然赢**（`runver 787878` 被排除在覆盖之外），因为那是显式用户请求；注意 `-nmtd 4` 走的 `gfdistribute` 给的是 k=0.005/0.00333、α=1.3/0.78 且 `rfac = 14/4 = 3.5`（35 ps/条），既更激进又长得多，**不是**本项的等价物。**不改 `-nmtd` 之外的任何行为**：`env%nmetadyn` 经 `env%allocate(4)` 置 4 ⇒ `nsim = 4`（`parallel.f90:945`），日志里 `list of applied metadynamics Vbias parameters` 会如实列出这 4 组；`search_conformers.f90:157` 的 `nmetadyn-2` 只在第 2 轮之后才生效，而 `Maxrestart=1` 下没有第 2 轮。**调参不用重编译**：`-mtdkscal <f>`（CLI，`confparse.f90:1228`）或 TOML `kscal`（`parse_xtbinput.f90:603`）可整体缩放 k，若 0.002 真的顶到 SHAKE，先试 `-mtdkscal 0.75` 而不是改源码。**门禁**：新增 8 条存在断言（含"覆盖块存在""被 `-nmtd` 让位""`nrem` 清零"）+ 1 条否定断言（旧的 `(1 MTD round, no iteration)` 文案必须消失） |
| 43 | `src/classes.f90`<br>`src/confparse.f90`<br>`src/parsing/parse_maindata.f90`<br>`src/algos/parallel.f90`<br>`src/printouts.f90` | **新增 `--loadconformer`：从工作目录读取 `readconformers.xyz`，把里面的构象作为第二个「静态 RMSD 偏置」加进每一条 MTD**。**机制与 `-iru` 完全相同**（同一段代码，不是复制一份）：`crest_search_multimd_init2()` 里 `npot = 2`，`mtd(1)` 还是原来的动态 RMSD 偏置（照旧堆高斯），`mtd(2)` 是 `cv_rmsd_static`，`cvxyz(3,nat,nall)` 在**这里一次性读好并以 Bohr 为单位**发给每个 job（`iruxyz = iruxyz/bohr`，与 `mol%xyz` 同单位），于是 `mtd_ini()` 走它的 `if (allocated(pot%cvxyz))` 分支，**没有任何 job 会再去开这个文件**——既避免 `nsim` 个读者并发，也避免与 cregen 重写该文件的竞争。**和 `-iru` 的唯一实质差别是「文件来源 + 生效时机」**：`-iru` 读 `conformerfile` = `crest_conformers.xyz`，那个文件是**多层优化（multilevel optimization）在跑完第一轮 MTD 之后才写出来的**，所以 `-iru` 最早只能从第 2 轮 MTD 开始生效；`--loadconformer` 读新增的公开参数 `loadconformerfile` = `readconformers.xyz`（`classes.f90`，紧挨着 `conformerfile`），由**用户**在运行前放到 CWD，因此**第一条 MTD 就带着偏置**——这正是加这个选项的目的：手上已经有同一分子的一批已知构象（上一次 crest 的结果、手动搭的构象、MD 快照）时，让这次采样从一开始就往没走过的区域推。**解析位置刻意放在 V2 分支（`confparse.f90`，紧挨 `-onebatch`）而不是 `-iru` 所在的 V2i 子块**，因为 V2i 子块外面还套着 `if (env%iterativeV2)`；放 V2 分支才能让 `--loadconformer` 与 `--quick` / `-squick` / `-mquick` / `--onebatch` 任意叠加。**两个开关同时给时 `--loadconformer` 赢**（`irufile` 先被赋成 `conformerfile`，再在 `env%loadconformer` 为真时改赋 `loadconformerfile`），因为用户显式提供的文件是更具体的请求、且它必然是完整的。**偏置强度**：`mtd(2)%kpush = mtds(i)%kpush * irukscal`，即沿用该 MTD 自己在 `$metadyn k alpha` 列表里的 k；`irukscal` 依次探测 `CREST_LOADCONF_KSCALE` → `CREST_IRU_KSCALE` → 默认 `1.0`，所以用户系综可以单独缩放，不必和 `-iru` 共用同一个值。**注意静态偏置是对系综里所有构象求和的，幅值随构象数线性增长**：几十~上百个构象时建议先设 `CREST_LOADCONF_KSCALE=0.1`（日志里 `k (static, max)` 那一行会打印实际值），否则结构容易被直接推开。**两道兜底**：① 解析时 `inquire (file=loadconformerfile,exist=ex)` 检查，文件不存在就打印 `**WARNING**` 且**不置位** `env%loadconformer`（不会静默变成无偏置的空跑）；② 读到文件但 `irunat /= env%ref%nat`（原子数不符）时打印 `**WARNING**` 并禁用偏置，不会用错配的坐标去算 RMSD。**原子掩码同样施加**：下面 `if (sum(env%includeRMSD) /= env%ref%nat)` 那段已经改成 `do k = 1, mddats(i)%npot`，所以第二个势和第一个势用同一套 `atinclude`；另外 `mtd(2)%cvdump_fs = 0.0` / `cvdumpstep = 0`，静态势从不清空/写出新的 CV。**入口三种**：CLI `--loadconformer` / `-loadconformer` / `-loadconf`；TOML 根键 `loadconformer = true`（`parse_main_bool`，同样带 `inquire` 兜底）；`-h` 帮助里有 6 行说明。**日志**：Job settings 段新增 `static RMSD bias     : "-loadconformer", ensemble read from "readconformers.xyz" (active from the 1st MTD)`（`env%iru` 为真而 `--loadconformer` 未给时，同一行会如实写成 `-iru` / `crest_conformers.xyz` / `from the 2nd MTD iteration on`）；MTD 启动前另有一组 `bias ensemble` / `bias structures` / `k (static, max)` / `k-scaling`。**门禁**：新增 16 条存在断言 + 2 条否定断言（`rdensembleparam(conformerfile` 与把 `conformerfile` 硬编码进 banner 的旧写法必须消失，否则 `--loadconformer` 会被悄悄忽略） |
| 44 | `src/algos/parallel.f90` | **新增「动态核数再分配」：某条 MTD/MD 提前跑完后，它那份核立刻交给还在跑的作业，而不是整轮空置**。背景：`crest_search_multimd*()` 给每条作业的是**一次性算好的固定切片** —— `new_ompautoset(env,'auto',nsim,T,Tn)` 算出 `Tn = T/nsim`，由 `ompenvset()` 导出成 `OMP_NUM_THREADS`，ORCA 再把它翻译成 `%pal nprocs Tn`。这个切片**此后从不再计算**。于是先结束的作业把它那份核一直闲置到整轮结束：4 条 MTD / 48 核 = 每条 12 核，早结束 1 条就空 12 核、2 条就空 24 核；而同一轮 MTD 步数相同，快出来的正是**尾部**那段时间。**做法**：新增 `module crest_dynshare`（放在 `parallel.f90` 里 `end module parallel_interface` 之后，不新增文件、不动构建清单），三个例程 + 7 个模块级共享状态（`ds_on / ds_njobs / ds_nslots / ds_total / ds_cap / ds_share / ds_share0 / ds_nfin`），全部只在 `!$omp critical (crest_dynshare_lock)` 内触碰。接线三处：① `new_ompautoset` 之后、OMP 区之前 `if (subproc) call crest_dynshare_init(nsim,T,env%Threads,Tn)`；② 每个作业 `dynamics()` 返回后（无论成功失败）`call crest_dynshare_release()`；③ `!$omp end parallel` 之后 `call crest_dynshare_reset()` 还原初始切片。**生效通道只有一条** —— `OMP_NUM_THREADS` 环境变量：ORCA 是 `orca_sc.f90` **每一步**都重写 `ORCA.in`（`write_orca_input` → `orca_nprocs()` → 读 `OMP_NUM_THREADS` → 决定这一步的 `%pal nprocs N`），xtb / turbomole / generic 则是子进程直接继承该变量。所以新份额发布后，**还在跑的作业在下一个单点步就吃到**，而**正在外部调用里的作业保持它已有的 ranks 不变** —— 这是交接，不是给活的计算改尺寸。**波形（nsim > slots）正确处理**：`nrunning = max(1,min(ds_nslots,ds_njobs-ds_nfin))`，作业排着队时一个跑完马上被补上，只有队列排空后份额才上涨。**只对子进程计算器生效**：进程内计算器（tblite / gfnff ...）把预算用在 `dynamics()` **内部**创建的嵌套 OMP 组上，那个尺寸从外面改不了，所以 `subproc` 为假时整个机制不启动、行为完全不变。**线程安全的一个细节**：`setenv()` 被工作线程调用时，其它线程可能正在 fork/exec。glibc 的 `setenv()` 在新值能装进旧槽位时是**原地改写**（只有值变长才会替换指针，而那正是会和并发 `getenv()` 抢的情况）。`crest_dynshare_init()` 因此**先把本轮可能出现的最大值 `ds_total` 写一次**，再写实际值 —— 之后所有值都 ≤ `ds_total`，于是每一次 `setenv()` 都是原地写，永不换指针。**三个旋钮**：`CREST_DYNAMIC_SHARE=0` 关掉（默认开）；`CREST_MAX_CORES_PER_JOB=<n>` 给单作业份额封顶（默认 0 = 不封顶）；`CREST_ORCA_NPROCS=<n>` 是既有旋钮，它会**钉死** ORCA 的 rank 数从而**覆盖**本机制 —— 两者同时生效时会打印告警，不会静默失效。**ORCA/MPI 务必注意**：48 个 MPI rank 跑 MINIX/PBE 这种小基组**很可能比 12 个还慢**（每轮 SCF 都要 Allreduce 一次 Fock 矩阵），建议先试 `CREST_MAX_CORES_PER_JOB=16`~`24`；不封顶就是用户要求的「平分全部核数」。**日志**会打印 `[cores] dynamic sharing: 4 jobs, 12 cores each - ...` 以及每次变化时 `[cores] 1 of 4 jobs done -> 3 still running, 16 cores/job (was 12)`。**`crest_search_multimd` 与 `crest_search_multimd2`（`-nci` / `-entropy`）两个驱动都接了线**，门禁里有 >=2 次计数断言防止只接一个。**门禁**：新增 9 条存在断言 + 3 条「恰好 2 次」计数断言（init / release / reset 各要在两个驱动里出现一次）。**未做（同类但非本次需求）**：`crest_search_multimd` 开头那段**串行**的 worker 初始化（`do i=1,T; call engrad(...)`，每个槽一次完整外部调用）同样只用 12 核，理论上可以先临时提到 48 核；它和提前结束无关，没有动 |
| 45 | `src/dynamics/dynamics_module.f90` | **进度行压到 73 列（原 108 列），顺便修掉 `60 sec` 这个进位 bug**。背景：4 条 MTD 在 `!$omp critical` 里**交错**打印进度，原行 108 列在窄终端上会折行，一折行就完全读不出是哪个 job 的哪一段。改法三处，都在 `md_progress()` / `md_sec2str()` 里：① **删多余空格** —— `"%  ("`→`"% ("`、`" ps,  elapsed "`→`" ps, elapsed "`；② **`f7.2`→`f6.2`**（两个 ps 字段各省 1 列，仍可表示到 999.99 ps，生产 MTD 最长几十 ps，不会溢出成 `******`）；③ **时间串紧凑化** —— `3 h 32 min 21 sec`→`3h32m21s`，分钟和秒用 `i2.2` **补零**，于是长度固定（`1h08m21s`、`10h37m03s` 都是同一宽度），列不会跳。**刻意没动的空格**：`1x,a`（tag）与 `i3`（百分比）的 padding **必须保留** —— 它们正是 1..N 号 job 之间列对齐的唯一来源；`*MTD   2` 里那 3 个空格是 `write(atmp,'(a,1x,i3)')` 的固定宽度，删掉会让 `*MTD 12` 与 `*MTD  2` 错位一列。**标签 `expected remaining time` → `remaining`**：在 `elapsed X, remaining Y` 这个成对语境里语义无歧义，且**没有**退回最早那个被否掉的 `ETA`（断言里仍钉着"必须保留显式标签"这件事）。**顺带修的 bug**：`md_sec2str` 用 `anint()` 取秒，59.6 s 会进到 60，于是出现过 `3 h 31 min 60 sec`；现在在**选格式之前**先做进位规范化（`seconds>=60`→进分，`minutes>=60`→进时），旧格式下这个缺陷只是不明显，不是不存在。**实测列数（同一批样本）**：`*MTD   2  25%` 那行 108→**73**；`Trial MTD  25%` 87→**63**；`*MTD  12  75%` 103→**70**；最长样本 108→**73**（省 32%）。**门禁**：新增 3 条存在断言（新格式串、`i2.2` 紧凑时间串、60 秒进位规范化）+ 3 条否定断言（`" ps,  elapsed "` 双空格、`i0," h ",i0," min ",i0," sec"` 旧时间串、`", expected remaining time "` 长标签必须消失）；原来那条钉 `expected remaining time` 的存在断言已改写。**一个坑**：源码注释里保留了旧格式作为 before 对照，里面自然含 `ps,  elapsed` 字面量，所以否定断言的 pattern **必须带上前后双引号**写成 `" ps,  elapsed "`，否则会匹配到注释而假红（本次就踩了） |
| 46 | `src/calculator/orca_sc.f90`<br>`src/algos/parallel.f90` | **修「环境变量开关只看存在、不看值」的 bug：`CREST_ORCA_PRINT_CALL=0` 会把本来要关掉的回显打印出来**。现象：用户 `export CREST_ORCA_PRINT_CALL=0`，日志里照样出现 `ORCA call: cd /tmp/... && OMP_NUM_THREADS=1 ... orca ORCA.in "--bind-to none" > ORCA.out 2> ORCA.err`。**根因**：`ORCA_setup()` 里的判据是 `if (len_trim(envstr) > 0.and..not.callprinted)` —— 只判断变量**是否非空**，**从不解析它的值**，所以 `0` / `false` / `no` 全都算"开"。同一个函数里 `CREST_ORCA_NO_OVERSUBSCRIBE`（:217）是同款 bug 且方向相反：`if (len_trim(envstr) == 0)` 判据下，设 `=0` 反而**删掉** `--oversubscribe`，与变量名含义完全相反。**改法**：统一成树里既有的整数解析约定（与 `parallel.f90` 的 `CREST_DYNAMIC_SHARE` 完全一致）——`ival = 0` → `read (envstr,*,iostat=io) ival` → `if (io == 0.and.ival /= 0.and..not.callprinted)`（PRINT_CALL，默认关）/`if (io /= 0.or.ival == 0)`（NO_OVERSUBSCRIBE，默认开）。只认整数，未设置或不可解析时回落默认值。**顺带修的一处未定义行为**：Fortran 的内部 `read(...,iostat=)` **失败时不会给变量赋值**，而 `.and.` / `.or.` **不保证短路求值**，所以 `if (io == 0.and.ival /= 0)` 在 `io /= 0` 时仍可能读到未初始化的 `ival` —— 必须先 `ival = 0`。`parallel.f90:183` 的 `CREST_DYNAMIC_SHARE` 有同样的隐患（未设置时 `r /= 0`、垃圾 `ival` 若为 0 会误判成"关闭"），一并加上 `ival = 0`；这一行**不改变**逻辑（read 失败时 `r == 0` 已是假），只消除 UB。**不重编译的临时绕过**：`unset CREST_ORCA_PRINT_CALL` 或 `export CREST_ORCA_PRINT_CALL=`（赋空串即等价于未设置）。**门禁**：新增 3 条存在断言（整数解析、`=0` 真的关闭、`=0` 保留 `--oversubscribe`）+ 2 条否定断言（两个 `len_trim(envstr)` 旧判据必须消失）。ORCA-HOWTO.md 的开关表补上 `CREST_ORCA_PRINT_CALL` 并加了一句「按值解析，`=0` 就是关」的说明 |
| 47 | `src/parsing/parse_xtbinput.f90`<br>`src/parsing/confparse2.f90`<br>`src/confparse.f90`<br>`src/printouts.f90` | **`--cinp` 的键保持 RIGID（SHAKE），而 `-cbonds` / `-cmetal` / `-cheavy` / `-clight` 恢复成谐振约束**。背景：`get_xtb_constraint_block()` 是**两个读者共用**的——`confparse.f90:2361` 读 `--cinp`（或 TOML `constraints`）指名的文件，`confparse2.f90:157`（`internal_constraint_repair`，**上游代码、本树未改**）读 `-cbonds` 自动生成的 `bondlengths`。本树加的 `add_hardcons()` 挂在它里面，于是**两个来源的键都被注册成刚性**。`--cinp` 那样是对的（用户手写的那几根键正是要钉死），但 `bondlengths` 是**每根键一条 `distance:`**，全变 SHAKE 会把整个分子（`-cmetal` 时是整个第一配位球）冻成刚体——这恰好违反 `classes.f90:396-399` 与 `parse_xtbinput.f90` 里我们自己写的注释。**改法**：给 `parse_xtb_inputfile()` 加一个 **opt-in** 的 `logical,optional,intent(in) :: hard`，**默认 `.false.`（= 上游的谐振行为，fail-safe）**，向下透传给 `get_xtb_constraint_block(env,blk,hard)`，两处 `add_hardcons()` 调用改成 `if (hard) call add_hardcons(...)`。两个调用点都**显式**写清意图：`confparse.f90` 传 `hard=.true.`（用户手写文件），`confparse2.f90` 传 `hard=.false.`（自动生成文件），这样读代码的人不必去查默认值。**为什么不用 `hard=env%hardcons`**：`env%hardcons`（`-softcons` / TOML `hardcons`）本来就在 **apply 时**把关（`choose_settings.f90:546 if (env%hardcons.and.env%nhardcons > 0)`），提前到注册时判断会引入对 TOML/CLI 解析顺序的依赖；这里只改"哪些文件有资格变刚性"，`-softcons` 的语义完全不变。**实际效果**：`-cmetal 0.05` 现在在 **MD 里**也真的按 `V=½k(r−r₀)²` 起作用了（此前 SHAKE 把键钉死、弹簧恒为 0，力常数只在优化里有效）；`--cinp` 不受任何影响。**代价**：`src/parsing/confparse2.f90` 此前是未改动的，本次为写显式 `hard=.false.` 而进入补丁集，**本地改动文件数 36 → 37**。**`-h` 帮助**在 `-cbonds` 下补了两行说明「harmonic restraints - never rigid, unlike -cinp / (the same holds for -cmetal/-cheavy/-clight)」。**门禁**：新增 9 条存在断言（`optional :: hard` 签名、`hardcons = .false.` 默认值、`present(hard)` 覆盖、向下透传、`get_xtb_constraint_block` 收到、两处 `if (hard)` 守卫、两个调用点的显式 `hard=`）+ 3 条否定断言（两处无守卫的 `call add_hardcons(...)` 与未带参数的 `parse_xtbinputfile(env,'bondlengths')` 必须消失，防止有人把守卫删掉退回旧行为）。**校检**：fparser 3 文件 OK；semantic 0 ERROR / 0 UNRESOLVED（164）；mkpatch **37 文件**；verify_pal **37 IDENTICAL / PASS**；verify-fixes **253 OK / 0 FAIL / PASS**。未做（同一族、仍需确认）：`confparse.f90:1416` 的 `-fc` 守卫写成 `if (i+1 >= nra)` 而文件里其余 ~15 处同型守卫都是 `i+1 <= nra`，即 **`-fc` 只有出现在命令行最后两位时才真正赋值**（上游 bug，在 `crest-threadfix.patch` 里是未改动的上下文行；另外 `-fc` 作为最后一个 token 时 `arg(i+1)` 越界）——所以设力常数请用 `-cmetal 0.05` 这种紧跟写法 |
| 48 | `src/dynamics/dynamics_module.f90`<br>`src/algos/setuptest.f90`<br>`src/algos/dynamics.f90`<br>`src/algos/parallel.f90`<br>`src/printouts.f90` | **新增「stop 文件」：动力学跑着跑着可以在工作目录里 `touch stop`，采样立刻收尾并把已采到的结构交给下一步优化**。背景：ORCA 这种外部计算器下 6 ps 的 MTD 要跑几天；中途 Ctrl+C / 杀作业什么都拿不到——因为**还没有任何结构被优化或写出**。现在只要在工作目录出现一个名为 `stop` 的文件，最近的那个检查点就会结束该条 MD，轨迹照常写盘、照常进 `collect` → 多层优化 → cregen，用户几分钟到几十分钟后就能拿到一套（较小的）构象系综。**轮询位置刻意放在 MD 步的开头**（`MD: do t = 1,dat%length_steps` 里 `call initsignal()` 之后），而不是结尾：目的正是**跳过下一步的 `engrad()`**，对外部计算器那是整个循环里最贵的一步；检查本身每步只多一次 `system_clock()`（几纳秒），每间隔才一次 `stat()`，相对它省掉的那次梯度调用完全不可见。**按墙钟而不是按步数定时**：GFN2 一步几毫秒、ORCA 一步几十秒到几分钟，按步数设间隔在这两者之间会差好几个数量级。**两个环境变量**（沿用树里既有的 `CREST_*` 约定，且**按值解析**而非只看存在）：`CREST_STOP_CHECK_SEC=<秒>` 间隔，默认 **300**（5 分钟），`0` 或负数 = 关闭；`CREST_STOP_FILE=<名字>` 改文件名，默认 `stop`。单位写成 `_SEC` 是有意的——`5` 到底是 5 秒还是 5 分钟必须从名字就能看出来。**触发条件是「文件是在 crest 运行期间出现的」**：`md_stopcheck_setup()` 在**本进程第一次** `dynamics()` 时记录一次 `stopfile_seen0`，之后判据是 `hit = (ex .and. (.not. pre))`。这一条是必需的——否则上周停过的那次留下的 `stop` 会让今天在同一目录新起的那次跑 5 分钟就死，而日志里唯一的线索只是「MD stopped on user request」，用户根本无从猜起。**不删文件**：触发时只打印提示「Remove the file before the next run」，不去做删除/改名这种用户没要求的副作用。**新的终止状态 `term = 3`**：不能复用 `term = 2`（"terminated during MD"，那同时是 SHAKE 不收敛 / engrad 失败的状态），因为调用方分不开；`crest_moleculardynamics()` 原来只有 `if (io == 0) ... else ...` 两支，`term=3` 会落进 else 被报成 `WARNING: MD run terminated ABNORMALLY` 并置 `status_failed`——**要的轨迹明明写出来了却说失败**。所以改成 `select case`，`case (3)` 单独一条；`parallel_md_finish_printout()` 同样加 `case (3)`，把完成行写成 `*MTD   2 stopped ("stop" file)`，让「为什么这条 5 分钟就结束了」在日志里自解释。**`mddata` 新增 `logical :: stopcheck = .true.` 只给调用方一个 opt-out**，唯一用到它的是 **`setuptest.f90`（trial MTD）必须关掉**：那个循环把 `io /= 0` 读成「没跑通 → 时间步减半重试」，六次之后 `status_safety` **中止整个 crest 运行**。一个 stop 文件若能让 trial 触发，用户会得到六行 `Trial MTD n did not converge!` 加一个彻底失败的 run——与 stop 文件的目的恰好相反。**报告无条件打印**（`md_stopreport()`，自带 `!$omp critical`）：并行 MTD 是 `pr=.false.` 的静默运行，没有它日志里就什么都不剩。**补做（同轮追加）**：`stop` 会在**工作目录**和**启动目录**两处都找，所以放在 toml/xyz 那个目录里在任何模式下都生效。判据是 `md_stopfile_exists()`：`inquire` 当前目录 → 若未命中且 `stopfile_launchdir` 非空且名字是相对路径，再 `inquire` `trim(launchdir)//'/'//trim(fname)`；绝对路径不做拼接。启动目录由 `parseflags()` 在 `getcwd(env%homedir)`（`confparse.f90:133`）之后立刻用 `md_stopcheck_setdir()` 交给 `dynamics_module` —— **必须在任何 `chdir` 之前**。`confparse.f90` 本来就 `use dynamics_module`，所以没有新增模块依赖。**澄清一个易混点**：`CREST_SCRATCH_DIR`（环境变量，`iomod.F90:1205`）只是把 `calculation.level.N` 用符号链接挪走，**从不 `chdir`**，因此它本来就不影响 stop 文件；真正会 `chdir` 的是 CLI 的 **`-scratch`**（`env%scratch` → `scratch.f90:54`）。现在两种模式走同一条代码路径。**已知边界（必须知道）**：① 一次运行里若有多段动力学（V2 的 MTD 迭代、`crest_rotamermds` 那段附加常规 MD），每段都会在**各自跑满一个间隔之后**停下，不会在第 0 步就空手结束——因为这个检查只在间隔到达后才做。③ 只作用于动力学；多层优化 / cregen 阶段不检查。`-h` 帮助在 `-vbdump` 下面加了 9 行说明。**本地改动文件数不变（37）**，五个文件都在既有补丁集里。**门禁**：新增 22 条存在断言（opt-out 字段、`stopfile_initialized` 进程级参考、两个 `get_environment_variable`、`if (rd == 0)` 按值解析、`< 1.0_wp` 下限、`stopfile_seen0 = ex`、轮询调用点、`hit = (ex.and.(.not.pre))`、`term = 3`、`md_stopreport`、trial MTD 的 `.false.`、`case (3)` 两处、`-h` 文档，以及启动目录那 6 条：`stopfile_launchdir` 字段、`public :: md_stopcheck_setdir`、`confparse.f90` 里的 `md_stopcheck_setdir(env%homedir)`、`md_stopfile_exists` 函数、两处调用点改用它、绝对路径不拼接）+ 3 条否定断言（`MD%stopcheck = .true.`、`if (io == 0) then` 旧两分支、以及只看 CWD 的 `inquire (file=fname,exist=ex)` 必须消失）。**校检**（无本地 Fortran 编译器，以下仅为语法/语义/补丁级）：fparser 5 文件 OK；semantic 0 ERROR / 0 UNRESOLVED（164）；mkpatch **37 文件**；verify_pal **37 IDENTICAL / PASS**；verify-fixes **278 OK / 0 FAIL / PASS**（上轮 253，+25） |
| 49 | `src/choose_settings.f90`<br>`src/printouts.f90` | **`--onebatch` 的 Vbias 栅格第二个 α 由 0.50 改成 0.67**（k 不变，仍是 0.0020/0.0010）。背景：α 决定**单个高斯包的宽度**，半高半宽 = `sqrt(ln2/α)` —— 0.50 → **1.177 Å**，0.67 → **1.017 Å**（收窄 13.6%）。α=0.50 时一个包太宽，沉积几次就把一整片 RMSD 区域填满，结构还没把所在势阱采样完就被推出去，落到不理想的构象上；收窄到 0.67 可以在**同样的峰高 k** 下换取每个势阱内更多的步数。**改法**：`alpinc` 是**除数**（循环里是 `alp = alp/alpinc`），所以落进栅格的值是 `alp/alpinc`；写成 `alpinc = 1.0d0/alp2`（`alp2 = 0.67d0`）而不是直接写 1.4925…，因为 0.67 才是有意义的那个数、倒数对读者毫无信息。**新增 `CREST_ONEBATCH_ALPHA2=<real>`**（按值解析，读失败或非正回落 0.67）：调 α 不用重编译 —— 这个参数本来就是要反复试的。**两个独立的 if 而不是 `.or.`**：`read` 失败时 `alp2` 变成未定义，而 `.or.` 不保证短路求值（与本树第 46 行修 `CREST_ORCA_PRINT_CALL` 时踩的是同一个坑）。**实际栅格（已按循环逐点验算）**：MTD 1/2 → k=0.0020/0.0010、α=**1.00**；MTD 3/4 → k=0.0020/0.0010、α=**0.67** —— 改动正好落在「两个 α=0.5」那两条上，另外两条不变。**k 未动**，所以 `metadfac` 仍是 0.184 / 0.092 Eh（`rednat ≈ 92`），SHAKE 的初始违背量级不变。**栅格几何中心** `sqrt(1.00×α₂)` 由 0.707 升到 **0.819**，仍落在默认 V2（`0.281…1.3`）与 `-quick`（`0.3…1.2`）共同认可区间内、且离「过窄」那条边界（α ≤ 0.3 时 ~90 原子体系已分不清构象）更远。`-h` 帮助补了 3 行，直接写出 `k = 0.0020/0.0010, alpha = 1.00/0.67 (hill half width 0.83/1.02 A)`，不用再让人去源码里数。**门禁**：新增 4 条存在断言（`alp2 = 0.67d0`、`CREST_ONEBATCH_ALPHA2`、`alpinc = 1.0d0/alp2`、`-h` 里的 `alpha = 1.00/0.67`）+ 1 条否定断言（`alpinc = 2.0d0` 必须消失）；原来钉 `0.83 / 1.18 Ang` 注释的那条断言已改写为 `0.83 / 1.02 Ang`。**校检**（无本地 Fortran 编译器，仅为语法/语义/补丁级）：fparser 2 文件 OK；semantic 0 ERROR / 0 UNRESOLVED（164）；mkpatch **37 文件**；verify_pal **37 IDENTICAL / PASS**；verify-fixes **283 OK / 0 FAIL / PASS**（278 → 283） |

> 附带说明：早期的 `crest-toml-segv.patch` 内容与 `crest-xtbsys-mtd-fix.patch`> 里 `parse_calcdata.f90` 的两个 hunk **完全相同**，属于冗余，因此没有放进 `patches/`。

> **流程提醒（本轮踩过一次）**：改完 `src/` 之后，务必按上面「重新生成补丁」那一段的
> **三步顺序**走完——`mkpatch.py` 会先把新补丁写到 `/PATH/TO/crest-code/`（不是 `patches/`！），
> 再由 `verify_pal.py` 同步进 `crest-now/patches/` 并做「原始 zip + 补丁 == 工作树」往返比对，
> 最后才是 `verify-fixes.sh`。如果跳过前两步直接跑 `verify-fixes.sh`，A 段会报
> `crest-threadfix.patch` / `crest-xtbsys-mtd-fix.patch` **NOT-APPLIED**——看着像补丁丢了，
> 其实只是 `patches/` 里的副本还没跟上 `src/`。补一次 `verify_pal.py` 立刻转绿（本轮 104 OK / 0 FAIL、29/29 IDENTICAL）。

---

## 四、用 ORCA 代替 xtb / tblite

**结论先说：CREST 3.0.2 本来就内置 ORCA 子进程接口，不需要新写计算器代码**，
只需要配置。完整说明见 `ORCA-HOWTO.md`，这里给最短路径。

需要三个文件：

```
your-molecule.xyz      坐标
orca.inp               计算级别 / 基组   ← 模板（本目录已给示例）
crest-orca.toml        把两者接起来      ← 本目录 your-system-orca.toml
```

**推荐用法：ORCA 只做最后精修，MTD 用便宜的方法。**
一次构象搜索要在每个 MD 步求一次能量+梯度，也就是成千上万次单点；每次都跑一遍
完整 ORCA 通常跑不动。CREST 为此提供了 refinement 层级 —— 写两个 level：

```toml
input = "your-molecule.xyz"

# level 1：便宜，驱动成千上万步 MTD / 预优化
[[calculation.level]]
method = "gff"                  # 或 "gfn2" / "gfn"（外部 xtb）

# level 2：ORCA，只在精修阶段用
[[calculation.level]]
method        = "orca"
orca_cmd      = "/PATH/TO/orca"
orca_template = "orca.inp"
chrg          = 1
multiplicity  = 1               # 2S+1，注意不是 nα-nβ
refine        = "opt"           # "opt" = 对整个集合重新做几何优化
                                # "sp"  = 只做单点重新排序（最省）
```

`refine` 的值：`"sp"`(1) / `"add"`(2) / `"opt"`(3)；不写则这个 level 不参与精修。
日志里出现 `> Geometry optimization of N structures` 就说明 ORCA 级别真的生效了。

**如果你就想让 ORCA 全程驱动**（小分子 + 便宜方法），把 ORCA 写成唯一的
`[calculation.level]` 即可，见 `ORCA-HOWTO.md` 第 3 节。先跑 trial MTD
（现在默认 0.25 ps，见下节）估时间，别直接上完整搜索。

### 两个必须知道的坑

1. **自旋用 `multiplicity`，不要用 `uhf`。**
   CREST 把 `uhf`/`multiplicity` 存在同一个字段里，写出 `*xyz chrg mult` 之前只做
   "不足 1 就取 1" 的钳位。所以 `multiplicity = 2` 给出正确的二重态；而
   `uhf = 1` 同样会被钳成 1，即**变成单重态** —— 错的。开壳层务必用
   `multiplicity = 2S+1`。

2. **模板里一旦出现 `%pal`，CREST 会往命令行追加 `--oversubscribe`。**
   它拼出来的是 `<orca_cmd> ORCA.in --oversubscribe > ORCA.out`；普通的 `orca`
   二进制可能不认这个参数。两个安全做法：模板里**不写 `%pal`**（用 `-T` 做进程级
   并行），或者把启动器连同 `--oversubscribe` 一起写进 `orca_cmd`，例如
   `orca_cmd = "mpirun --oversubscribe -np 4 orca"`。

### 电荷 / 多重度 / 坐标从哪来

CREST 会把模板里所有坐标相关的内容剥掉，每次调用自己追加

```
*xyz <chrg> <multiplicity>
...
*
```

所以**坐标和电荷/自旋都由 CREST 输入文件决定，`orca.inp` 里写的 `*xyz` 会被丢弃**。
模板里也不要出现 `%coords`、`$new_job`、`! Opt`/`! Freq`/`! MD`，否则会被直接拒绝。
`EnGrad` 若缺失会自动补到第一条 `!` 行上。

### 输出与排错

每次调用的工作目录是 `calculation.level.<n>`（并行 worker 是
`calculation.level.<n>_<i>`），最后一次调用留下三个文件：`ORCA.out`（ORCA 正常输出）、
`ORCA.err`（ORCA 的 stderr —— **启动类错误只能在这里看到**）、`ORCA.engrad`（梯度）。
`crestopt.log` 只有 CREST 自己的优化记录，不含 ORCA 的报错。

**遇到 `Initial geometry optimization failed!` 先看 `calculation.level.1/ORCA.err`。**
最典型的原因是 `*_shared_openmpi418*` 版的 ORCA 需要 `LD_LIBRARY_PATH`
（否则报 `libmpi.so.40: cannot open shared object file`）。另外，用
`bash run.sh` 跑脚本时**不会**读 `~/.bashrc`，所以只在 `.bashrc` 里 export 过的变量
到不了 CREST 子进程里。完整排错表见 `ORCA-HOWTO.md` 第 8、9 节。

**如果 `ORCA.err` 是空的、而 `ORCA.out` 里有**
`ERROR: non-ASCII character(s) found in words ... [qcscan1.cpp, line 264]`，
那不是环境问题，而是修复 #10 那个 bug：`ORCA.in` 被写成了二进制垃圾。
用 `file ORCA.in` / `grep -c $'\x00' ORCA.in` 一看就知道。这种情况**必须重新编译**
（或者先用 `ORCA-HOWTO.md` 第 10 节的「零行模板 + 包装器」临时绕过，不用重编）。

---

## 五、性能：45 秒那个数字到底是什么

**它不是预测，是"理想线性加速"的下界。** `src/algos/setuptest.f90` 的
`timeestimate()` 就是拿单线程实测时间除以 `floor(线程数 / MTD 数)`：

```fortran
w = float(threads)/float(nmetadyn)      ! 40/6 = 6.67
t = float(floor(w))                     ! = 6
rtime = rtime/t                         ! 268 s / 6 = 44.7 s
```

它**假设 xtb 在 6 线程上能加速 6 倍**。实测只有约 **1.56 倍**，所以 2 分 51 秒
才是这套配置的真实成绩。

### 症结：6 个 MTD × 6 线程，把 40 核用成了等效 9 个核

子进程模式下 CREST 这样分线程（`parallel.f90` → `new_ompautoset`）：

```
parallel_jobs = nmetadyn                每个 MTD 一个 CREST 线程
cores_per_job = floor(-T / nmetadyn)    → OMP_NUM_THREADS 给 xtb
```

用你日志里的数字算（500 步/MTD，dt=5 fs）：

| 配置 | 单次 xtb 调用 | 每核每秒调用数 |
|---|---|---|
| 1 线程（trial MTD 实测） | **536 ms** | 1.87 |
| 6 线程（批量实测） | **343 ms** | 0.49 |

**xtb 线程级并行每核效率只有 1 线程时的 1/3.8**，而**进程级并行是完美的**
（每个 xtb 独立进程）。现在的配置是 6 个作业 × 6 线程 = 36 线程，
等效吞吐只有 `3000 × 0.536 / 171.7 ≈ 9.4` 个核 —— 40 核里浪费了 3/4。

### 解法：让 `cores_per_job = 1`，把线程换成人手 —— `-nmtd`

CREST 支持 `-nmtd <N>`（`confparse.f90:1283`，`choose_settings.f90: gfdistribute`），
它把 N 分配到 k×α 的 Vbias 网格上（`a=ceiling(sqrt(N))`，`k=nint(sqrt(N))`）。
总 MD 时间 `env%mdtime`（默认 15 ps）会**平分**给 N 个 MTD：

| `-nmtd` | 线程/MTD | 每 MTD 步数 | 单步耗时 | 批次墙钟（同 15 ps 总量） |
|---|---|---|---|---|
| 6（现状） | 6 | 500 | 343 ms | **172 s** |
| 12 | 3 | 250 | ~390 ms | ~98 s |
| 20 | 2 | 150 | ~410 ms | ~62 s |
| 40 | 1 | 75 | 536 ms | **~40 s** |

> 直接跑 `bash bench-nmtd.sh 40 6 12 20 40` 实测，别信表里的插值。

```bash
crest your-system-xtb.toml -T 40 --quick --noreftopo -nmtd 40
```

⚠️ **`-nmtd` 必须放在 `--quick` 后面**：两者都写 `env%runver`，
`--quick` 会把它设回 2，从而覆盖 `-nmtd` 的 MTD 网格选择（`nmetadyn` 仍会被设成 40，
但网格只有 6 个 → 数组不一致）。

⚠️ MTD 长度会随 N 变短（`t(MTD) = mdtime/N`）。40 × 0.375 ps 可能太短而削弱 MTD 的
驱动效果 —— 若要保留每个 MTD 的实际时长，用 `-mdtime` 把总量一起放大
（例如 `-mdtime 60 -nmtd 40` → 每个 1.5 ps，总工作量是现在的 4 倍）。

### 关于 `ratio c/w < 1`：这不是变慢，是子进程模式的必然

`src/timer.f90: timing()` 用的是 Fortran `CPU_TIME`：

```fortran
call cpu_time(time_cpu)      ! glibc: clock() = 只统计【本进程及其线程】
```

它**不包含 fork 出去的 xtb 子进程**。子进程模式下 40 个 CREST 线程全堵在
`waitpid` 上，所以 cpu-time 必然远小于墙钟 —— `0.034` 恰恰说明
「CREST 自己在等 xtb」，而不是程序变慢了。在 tblite 进程内模式下才会看到 c/w > 1。
**判断有没有变慢要看墙钟（wall），不要看 c/w。**

（同理，`htop` 里的红色是内核态时间，来自 `fork /bin/sh` + `fork/exec xtb` +
每步删 3 个文件 + 写/读 xtb 的各个文件。把工作目录放在原生 ext4 而不是 `/mnt/...`
可以显著削掉这一块。）

---

## 六、排错

| 现象 | 处理 |
|---|---|
| `verify-fixes.sh` 报 `NOT-APPLIED` | 说明这不是打完补丁的树 —— 别编译，重新拿一份 `crest-now` |
| 找不到 `gfortran` / `cmake` | `build-now.sh` 会先做工具链检查；gfortran 需要 ≥ 10 |
| 编译很慢 / 卡住 | 确认不是在 `/PATH/TO/mount` 下编译；装个 `ninja` 会明显更快 |
| 运行时报 `command not found: crest` | 先 `source ./crest_env.sh 40`（或把 `$PREFIX/bin` 加进 PATH） |
| `*MTD n terminated EARLY` 又出现 | 跑 `bash verify-fixes.sh` 确认补丁在树里；再看 `grep -i thread calculation.level.*/_*/xtb.out` |
| 还是比预期慢、CPU 不满 | 用 `check-xtb-threads.sh` 看子进程线程数与 env 是否符合预期（应全为 6） |
| 编译停在 `src/dynamics/dynamics_module.f90.o`，报 `Type of entity at (1) is type incompatible with source-expr at (2)` | 同一个 `allocate` 里混用了不同类型的分配对象（`source=` / `mold=` 要对**每一个**对象都成立，不是只管第一个）。每种类型拆成一条 `allocate`，实型那条用 `0.0_wp` 而不是 `0`。详见第三节第 33 行 |
| 想在不编译的情况下先排查同类错误 | `py /PATH/TO/crest-code/_orig/f90_semantic_check.py src`（传**目录**，不要传文件清单）。这是纯标准库脚本，不需要 fparser，164 个文件约 0.8 秒 |

---

## 七、备份与来源

- 未打补丁的原始源码：`/PATH/TO/crest-code/crest-latest-backup.zip`
- 参考工作树：`/PATH/TO/crest-code/crest-now`
  （2026-09-14：原开发树 `/PATH/TO/crest-code/crest-latest/crest-latest` 已被删除；
   补丁改由 `crest-now/src` 与原始 zip 直接 diff 生成。两者对 `src/` 是等价的，
   重生成后仍是同样的 25 个改动文件，可据此确认没有丢失修复。）
- 重新生成补丁（**顺序不能变**，否则 `verify-fixes.sh` 会假红）：
  ```bash
  python /PATH/TO/crest-code/_orig/mkpatch.py     # 重生成 /PATH/TO/crest-code/*.patch
  python /PATH/TO/crest-code/_orig/verify_pal.py  # 同步到 patches/ + 原始 zip 往返比对
  bash ./verify-fixes.sh                    # 期望 RESULT: PASS
  ```
  **顺序敏感**：只要动过 `src/` 就必须重跑 `mkpatch.py`。`verify-fixes.sh` 的 A 节做的是
  「补丁反向试运行」，源码和补丁一旦不同步就会立刻报 `NOT-APPLIED`（第 33 项那次就是这么
  发现的）——这是**预期行为**，不是树坏了，重生成补丁即可。
  另外还有一道不依赖编译器的关卡：`/PATH/TO/crest-code/_orig/f90_semantic_check.py`，查 fparser
  看不见的类型语义错误（`allocate` 的 `source=` / `mold=` 一致性），已并入 `verify-fixes.sh`
  的 C 节。
  `verify_pal.py` 会顺手把三个补丁复制进 `crest-now/patches/`，并做
  「原始 zip + 补丁 == 工作树」的往返校验。它校验的文件清单是**从补丁里
  解析出来的**，不再手工维护（以前手工清单过期导致过两次假 FAIL）。
  注：`mk_crest_now.py` 依赖已删除的开发树，目前**不可用**；
  `crest-now` 本身就是可直接编译的完整源码树。
