# CREST 3.0.2 · QW-Zhang-PKU modified (@ Yu group, PKU)

English README is given below the Chinese part (translated by AI).

本目录是 CREST 3.0.2 的修改版，由QW-Zhang完成。本项目旨在解决CREST的长期存在的一些痛点，让其更适配 [有机（尤其是金属有机）] 的机理研究。修改代码时使用了Hy4和Deepseek，修改后的功能已经过测试，具体内容如下。本文件绝大部分内容为作者本人撰写，请放心阅读。

联系邮箱：<2201110483@stu.pku.edu.cn>

代码修改内容细节见 `LOCAL-CHANGES.md`，许可证与第三方声明： `NOTICE`、`LICENSE`、`COPYING` 和 `COPYING.LESSER`。

## 一、编译与运行

### 编译

源码里已内置全部 subprojects，GNU 工具链下可以完全离线编译：

```bash
bash build-now.sh                               # 装到 $HOME/crest-dist/bin/crest
```

`build-now.sh` 自检后委托 `build_crest.sh` 做 CMake + Ninja 编译。默认精简构建会关掉 tblite / GFN-FF / libpvol / lwoniom，原因见第二节第 5 条；需要恢复上游全部特性用 `FULL=1`。

编译机与计算节点的 BLAS 环境不一致时，可先编一份静态 OpenBLAS：

```bash
bash build-openblas-static.sh
OPENBLAS_DIR=$HOME/opt/openblas-static STATIC_BLAS=1 bash build-now.sh
```

完整构建说明见 `README-BUILD.md`。

### 快速开始

```bash
source ./crest_env.sh 40
crest molecule.xyz -T 40 --quick
```

计算级别主要介绍 gfn2-xtb、g-xTB、ORCA。

**GFN2**（本修改版走外部 xtb 子进程）可以直接用命令行指定，并确保环境变量 PATH 里有 xtb 可执行文件，且`XTBPATH=`指向 xtb 可执行文件：

```bash
crest molecule.xyz -T 40 -gxtb_dev --quick
```

或者用toml格式作为输入文件：

```toml
input = "molecule.xyz"

[calculation.level]
method = "gfn2"
binary = "/PATH/TO/xtb/bin/xtb"   # 可省略，用环境变量 PATH 的 xtb
flags  = "--gfn2"
```

**g-xTB** 可以直接用命令行指定，并确保环境变量 PATH 里有 g-xtb 对应的 xtb 可执行文件：

```bash
crest molecule.xyz -T 40 -gxtb_dev --quick
```

**ORCA** TOML 指定

```toml
input = "molecule.xyz"

[[calculation.level]]
method = "orca"                          # 采样级别

[[calculation.level]]
method        = "orca"
orca_cmd      = "/PATH/TO/orca"   # 改为实际的ORCA位置
orca_template = "orca.inp"               # ORCA的输入文件模板
chrg          = 0                        # 改为实际的电荷
multiplicity  = 1                        # 改为实际的多重度
refine        = "sp"
```

`orca_template` 是自己写的 ORCA 输入模板；CREST 每次调用自动补上 `EnGrad`，不需要在模板里指定 `%pal` ，CREST会根据实际情况填入。  
更多配置（多级别、MPI 参数、约束文件写法）见 `ORCA-HOWTO.md`。  
如果ORCA跑在 机械硬盘/nfs挂载/超算分布式储存 会变慢，建议给每个作业单独的 scratch (尽量是SSD)：

```bash
export CREST_SCRATCH_DIR=/tmp/crest-$SLURM_JOB_ID
```

以下是一个orca.inp示例：

```inp
! PBE D3 MINIX GCP(DFT/MINIS) CPCM(THF) looseSCF nopop miniprint

%maxcore 2000

```

## 二、新增功能（对比CREST 3.0.2，共 11+ 项）

### 1. 增加 g-xtb、ORCA 接口

**原因：** CREST内置的方法包括gfn-xtb、gfn2-xtb等。在这个项目开始前，我通常用gfn2-xtb比较多。但由于gfn2-xtb对于过渡金属配合物描述较差，有时会跑出完全错误的结果（比如金属上配位的一氧化碳被插羰，等），所以迫切的需要更高级的方法整合入CREST。

**g-xtb：** g-xtb 是 Grimme 课题组发展的新一代半经验方法，其精度显著好于 gfn2-xtb（尤其是针对过渡金属配合物），且计算速度比 gfn2-xtb 增加不多，所以是我本人目前在 CREST 中最常用的方法。注意 g-xtb 目前仅发表在预印本，后续可能仍有修改。且对溶剂化支持不好，所以 g-xtb 模式下不建议设置溶剂，默认保持真空即可。

**ORCA：** 有时某些体系可能更加棘手，必须DFT作为计算级别（即使很耗时），可以使用ORCA。我们完善并开启了CREST的ORCA接口。经过测试，ORCA下使用纯泛函+极小基（如：! PBE D3 MINIX GCP(DFT/MINIS) LooseSCF）是可以算的动的，但是耗时非常长，仅推荐在疑难杂症时使用。目前只测试了ORCA 6.1.1，老版本ORCA应该兼容。注意ORCA运行时建议设定**临时文件夹**在**SSD**上，否则在读写较慢的硬盘（机械硬盘或超算分布式存储）上运行会被拖累。

**更改** ：

- 用 `-gxtb_dev` 开启 g-xtb。
- 需要下载 g-xtb 的可执行文件（<https://github.com/grimme-lab/g-xtb>），需要让 `XTBPATH=` 指向其`share/xtb`
- 需要使用 g-xtb 溶剂化时，通过环境变量透传，例如 `export CREST_GXTB_SOLV='--gbe thf'` 仅做调试使用，一般打开后任务就会失败。
- ORCA通过toml文件进行对接，见上方使用方法。
- `CREST_ORCA_PRINT_CALL=1` 可以调出CREST是如何调用ORCA的，供调试使用。



### 2. 约束（-cinp）改为刚性约束

**原因：** 搜索过渡态构象时，需要锁定部分原子（比如正在形成C-C键的两个原子），因此需要用`-cinp`锁定它们。原版默认`-cinp`用谐振模型约束，且强度默认是0.5。这会导致过渡态在经过构象搜索后非常容易变回产物或反应物，（回到反应物/产物的放热很大，克服了约束的弹簧强度）从而无法得到有效结果。如果把谐振强度改为5，虽然可以有效地锁住反应中的原子，但是过大的谐振强度会让 MTD 不易收敛而报错，不得不减小步长为1~2 fs。这是可以说是填一个坑同时挖另一个坑，因为这样任务的总耗时会急剧增加。

**SHAKE：** cinp改用用刚性的SHAKE约束，这样确保谐振强度改为5或6这样的大值时，更不容易报错。

**更改：**

- SHAKE 的修正扫描从 Jacobi 改成 **Gauss–Seidel**（`src/dynamics/shake_module.f90:511-558`）。
- 最大扫掠次数 `maxcyc = 6000`（上游 250，实测某个真实体系要 >3000 次扫掠），收敛容差 `tolshake = 1e-5`（上游 1e-7）。
- 失败报告用 `es12.4` 打印迭代数/上限、违反约束的键数、原始步偏差、修正后偏差、容差等信息供调试；
- 保留退回谐振模型的开关：`-softcons`。
- `atoms:` / `elements:` 生成的 pairwise 键也补上了刚性登记。这类写法会把原子列表里**每一对**都展开成键（n 个原子是 C(n,2) 根），超过 3n−6 个自由度时约束集线性相关，代码会打印提示。
- 预优化和批量优化的环节并**未改为刚性约束**（后续可以改进，但需要修改优化器，且不是那么迫切，所以暂时没改）通过设一个很大的力常数（>5）可以做到几乎锁住（被锁定的原子仍会有非常微小的位移，但是几乎可以忽略）。
- `-cbonds`、`-cmetal`、`-cheavy`、`-clight` 自动生成的键约束仍然是**谐振约束**，不会把整个分子完全锁住。`-cmetal`是一个很有用的功能，可以防止配体从金属中心解离。



### 3. 改进多轮MTD的实用性（-iru等）

**原因：** 以下仅针对我涉及的有机化学机理研究：CREST在大部分任务下，至少会跑2轮MTD动力学采样。其本意应该是通过第一轮找到能量更低的构象，用这个构象进行第二轮，如果又发现能量更低的，则进行第三轮，不断迭代。但是我发现绝大部分的任务（单过渡金属中心+常见配体），一轮MTD已经可以把构象基本找全。第二次MTD通常是无用功（90%以上情况），只是起到“保险”的作用。如果是这种情况，实际上就导致：**跑两次参数完全相同的 N ps 的 MTD**，意义不大，甚至还不如跑一次 2N ps 的 MTD。

**-iru：**&#x4E00;方面，新增了只跑一轮MTD的模式（`-onebatch `, 下述）。另一方面，增加了-iru，使第二轮MTD可以读取第一轮已探索过的极小点区域，通过对其增加bias，可以更快走出第一轮的“已知空间”，将更多时间留给“未探索空间”，从而将第二次MTD充分利用。此外，还增加了readconformers的功能可以跨任务读取上一次任务已探索过的构象，对于想给某体系跑第二次CREST的情况提供。  
此外我们做到了，每轮MTD，偏置势强度k都逐轮放大，让每轮的探索强度递增，防止在已知区域不断回旋。

**更改：**

- iru 的静态偏置的 k 通过环境变量 `CREST_IRU_KSCALE=`，默认 **1.0**。
- 每轮 MTD 的 k 放大 **1.25 倍**（`CREST_MTD_BOOST`），累计上限 **4.0**（`CREST_MTD_BOOST_MAX`），从第 2 轮起乘到之前的值。 `CREST_MTD_BOOST=1` 可关掉。这个放大也会作用到 `-iru` 的静态势。
- 跨任务读构象确认是 `--loadconformer`（别名 `-loadconf`，或 TOML 中写 `loadconformer`）：从**当前工作目录**的 `readconformers.xyz` 读入已知构象，机制和 `-iru` 完全一样。



### 4. 改进分子柔性的计算方式

**原因：**&#x5206;子柔性值 flex 决定了MTD的时长： `t(MTD) = 3.0·exp(0.10·flex·max(1,rednat−8))` ps，下限 5 ps、上限 500 ps。MTD时长又间接影响了整个任务的长度。

在对过渡态找构象时，设置了`-cinp `时容易出现柔性虚高的问题。  
具体来说，如果发生分子内成环反应，而我用-cinp锁定即将成环的原子A和原子B。CREST会因为A-B键长达不到单键，识别不出这是一个环，从而判断这是链状分子。

这导致CRSET认为柔性很大，需要长时间的MTD。但实际上这个是环状分子，构象少，不需要那么久的MTD。  
此外，原版CREST中分子柔性的值容易随着键长出现跳变（两个看着几乎一样的结构，判断出的柔性差别极大），因此我们改进了其计算公式。

新版做了三个改进：

- **环判定连续化**：用 sigmoid `w = 1/(1 + exp(−16·(rco/r − 1)))` 给每根键一个 0~1 的权重，再用 Dijkstra 在 `−log w` 上找绕开 i–j 键的最强闭环路径，得到环尺寸 L 和连续隶属度 ρ（`RHOFULL = 0.98`），`ringf = 1 − (1 − 0.5·(1 − exp(−0.06·L)))·ρ`。原来非黑即白的判据变成插值，键长稍变不会让柔性突变。`CREST_FLEXI_DEBUG=1` 会逐键打印 `r/rcov`、是否被强制、ρ 等。
- **让 `-cinp` 声明的键参与成环**：这些键即使长于共价半径和也被强制插入 bond graph（用 `ffac = 1.3` 门控滤掉 `atoms:` 展开出来的伪键），并在环权重里直接给 `w = 1`。这样"即将成环"的长键能闭合环，柔性不再虚高。注意这里不是给被锁定原子打折，而是让环被正确识别后由 `ringf` 自动降柔性。
- **单键截断 + 换聚合方式**：`val = min(val_raw, 0.9)`，再把均方根换成 power mean `flex = (Σ val^p / m)^(1/p)`，`p = 1.5`。cap 削掉个别异常键的主导，p 小于 2 让整体不再被大值拉着走。设 `CREST_FLEXI_VALCAP=1.0` 加 `CREST_FLEXI_PMEAN=2.0` 可还原上游公式做 A/B 对比。
- 如果算出的MTD时长还是较长（较短），可以用 `-mdlen` 可以设置一个系数乘到其上。我实际工作中用 `-mdlen=0.8`



### 5. 改进CREST 3.0.2中线程挤占的问题（这是一个非常严重的问题）

**原因：**&#x539F;版CREST 3.0.2通过内置的tblite进行xtb计算。我们发现运行中，即使设置`-T 40`，在MTD阶段也会产生>80个进程，且互相抢占资源，造成运行速度不可接受地慢！  
这一点从Trail MTD预估的MTD耗时和实际耗时的区别就可以看出。实际MTD耗时往往是预估的10倍以上，大部分的时间是因为bug被浪费的，没有起到任何作用（纯浪费）。  
根因是 tblite 在CREST 进程内部另开自己的 OpenMP team；`-T` 虽然管得住外部 xtb/ORCA，但管不住它。

**更改：**

- 我们弃用内部的tblite，改为调用外部的xtb程序（就像CREST 2.x 老版本那样），避开了这个问题。修复后，MTD的速度是正常的，不会被浪费了。
- `gfn2` / `--gfn2` 从 `jobtype%tblite` 重定向到 `jobtype%xtbsys`，即调用外部 xtb 子进程
- gfn1因为用的较少，没有做重定向
- 构建脚本默认 `WITH_TBLITE=OFF`，编译时就不带进程内 GFN2，从源头去掉这个嵌套线程源；需要时 `WITH_TBLITE=ON` 可以开回来。

**以下其他bug的修复内容是AI写的：**  
线程预算的下发也一并修了：`ompenvset()` 把 `OMP_NUM_THREADS` / `MKL_NUM_THREADS` / `OPENBLAS_NUM_THREADS` 导出为 `cores_per_job`，并固定 `OMP_MAX_ACTIVE_LEVELS=1` 防止"线程里再开线程"；原先 POSIX `setenv()` 的 Fortran 绑定少声明了第三个 `overwrite` 参数，读到未定义寄存器，变量已存在时 `setenv` 会静默拒绝更新，子进程于是沿用调用 shell 的值。现在绑定改成三参数并显式传 1（补丁 `crest-omp-env-fix.patch` 记录的那次实测：40 线程 6 个 MTD，从 2m48s 回到接近预估的 44s）。

另外还修了一个表现很像的 bug：`xtb_sc` / `orca_sc` / `turbom_sc` 把梯度文件的**完整路径**缓存在 `calc%gradfile` 里，而并行 MTD 会给每个 worker 换 `calcspace` 后缀，结果所有 worker 都去读 trial MTD 那份陈旧梯度，积分几步就 SHAKE 失败、报 `*MTD n terminated EARLY`。现在每次调用按当前 `calcspace` 重新推导路径（见 `patches/crest-xtbsys-mtd-fix.patch`）。



### 6. 动态核数分配

**原因：**&#x4EE5;MTD举例，当我用36核同时跑6个MTD，每个MTD分到6核。如果1个MTD提早结束，那么就会有6核完全闲置。

现在我们将其改为：剩下5个任务平分这36核，每个任务分到7核，仅1核闲置。这可以提高整体CPU的利用率，并减少任务耗时。

对于xtb半经验任务，这么做对效率的提升存在，但是并不算特别大。但对于大体系+ORCA计算的任务来说，可以显著减少耗时（有时可以减少几小时）

**更改：**&#x76EE;前默认开启动态核数，输出的out文件中会显示：`[cores] dynamic sharing: ...`只在首次触发动态核数时打印一次。后续不打印，因为没什么有效信息。



### 7. 修复CREST中剔除结构的Bug

- 默认在拓扑检查中排除金属原子参与的键，减少金属体系因为配位数取整跳变而被误删构象。
- `cregen_topocheck()` 先构造 `topoexcl(:)` 掩码（`-notopo` 指定的原子 + 所有金属），再用 `bondtotopo_excl()` 做拓扑比较（`src/cregen.f90:896-930`，实现在 `src/ztopology.f90:339`）。金属判定用显式的 86 元素表，另外排除 Li、Be（`src/dynamics/shake_module.f90:27-32,188-204`），不是简单的 Z > 20。`excludeTOPOmetals` 默认开，`-nometaltopo` 或 TOML `metaltopo = false` 可关。
- `CREST_SKIP_TOPO` 默认开启：拓扑对不上的结构改成**保留**（并记一个 `SKIP` 标签），而不是直接丢掉；`=0` / `false` / `off` 回到上游行为。这是 108 原子 Rh 体系实测后改的默认值。
- EWIN 的能量零点改用"最低的拓扑**有效**结构"，不再是全局最低（`cregen.f90:1290`），避免一个拓扑坏掉但能量很低的结构把整批构象卡掉。
- refine 阶段新增 `crest_refine_topopresort()`（`src/algos/refine.f90:176-229`），在跑昂贵的单点之前先做一遍拓扑预筛，只筛拓扑，不做能量排序和 RMSD 去重。
- 被剔除的结构可以导出审计：`CREST_EXCL_XYZ=1` 会写 `crest_excluded_structures.xyz`，每帧注释里带 `CLASH` / `TOPO` / `EWIN` / `DUP` / `SKIP` 标签和原因（`src/exclfile.f90`，默认关）。查"这个构象去哪了"比翻日志快。
- `cregen_CRE` 的去重循环做了几处性能改动：缓存 masked 坐标避免 O(nall²) 内层重复计算、利用能量升序缩小 j 的扫描窗口、重原子子集每结构只抽一次。RMSD / 能量阈值本身没动。
- 新增的 `src/minregistry.f90` 是极小点登记表：同一优化阶段内，如果某个还在跑的优化被判定正落向已登记的极小点，就提前停掉并直接采用该极小点的几何与能量。每个优化阶段开始会清空（防止下一阶段更紧的阈值被自己上一阶段的结果误剪），`CREST_PRUNE_ANCHORS=0` 可关。



### 8. 优化了输出信息

- MD/MTD 进度固定显示 **25% / 50% / 75%**。
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

- CREST原版中 MTD 和 MD 启动时，输出参数块乱序，已被修复。
- ORCA 初始优化开始前打印 `optimization starting`，不至于长时间一点输出都没有 （有时预优化要好几十分钟）。



### 9. 信号与 stop 文件

**原因：**（尤其是在用ORCA作为计算引擎时，）有时动力学未跑到程序预设的 MTD 时长，其结构已经采样够了（比如：BINAP螯合配位的一个磷不配位了，甩出去了）继续跑 MTD 已无意义。

增加一个方式可以随时停止 MTD ，进入下一阶段（批量优化）。

- `stop` 文件：动力学阶段每个 MD step 顶部按墙钟间隔检查（默认每 300 秒一次 `inquire`），只认**运行期间新建**的 `stop` 文件，遗留的旧文件不会误触发。命中后以 `term = 3` 干净结束当前动力学，轨迹逐帧 `flush`，已采到的结构交给后续优化/筛选，整轮并行结束后自动删掉 stop 文件。
- 变量名 `CREST_STOP_FILE`，间隔 `CREST_STOP_CHECK_SEC`（`0` 关闭）。
- Trial MTD 不适用这个机制；多层优化和 CREGEN 阶段也不检查，进入这些阶段后 `touch stop` 不会结束。

此外一些别的修改：

- Ctrl+C、`kill`、SLURM 超时信号现在能真正终止 CREST 及其外部子进程，不再出现终端“假死”、只能 `kill -9` 的情况。



### 10. 工作目录优化

我在集群上使用crest，原版CREST默认的工作目录无法更改。比如我的文件在nfs挂载的局域网的机械硬盘上，这导致文件访问会拖慢运行速度。  
现在改为：高频读写的中间文件（`MDFILES`、`OPTIM`、各 `calculation.level.*` 计算目录）可以放到节点本地临时目录；低频的（`coord`、`wbo`、`crestopt.xyz` 等）默认收进工作目录下的 `crest_work/`：

```bash
export CREST_SCRATCH_DIR=/tmp/crest-$SLURM_JOB_ID   # 热文件，建议 SSD
export CREST_WORK_DIR=crest_work                    # 冷文件，默认值
export CREST_WORK_DIR=none                          # 回到上游布局
```

实现方式是符号链接，没有去改各个 `open()` 调用点。`crest_conformers.xyz`、`crest_rotamers.xyz`、`crest_best.xyz` 等结果文件仍留在当前工作目录。`CREST_STAGE_FILES=0` 可以整体关掉。



### 11. 采样流程优化

- trial MTD 更短：ORCA 默认 0.125 ps，其他计算器默认 0.25 ps；可用 `CREST_TRIAL_MD_PS` 覆盖。
- trial MTD 的 Vbias dump 默认固定 50 fs，并保证至少两次更新。
- 新增/接通：
  - `--onebatch`：只跑一轮 MTD。仍然包含 trial MTD + MTD + MD + crosscof。
  - `--loadconformer`：从 `readconformers.xyz` 读入已知构象，从第一轮 MTD 就加静态 RMSD 偏置。（见上面）
  - `-iru`：把已找到的构象作为第二个静态偏置叠加到动态 MTD 上。（见上面）



## 三、环境变量

| 变量                                          | 默认值                  | 作用                             |
| ------------------------------------------- | -------------------- | ------------------------------ |
| `CREST_DYNAMIC_SHARE`                       | 开                    | `=0` 关闭动态核数再分配                 |
| `CREST_MAX_CORES_PER_JOB`                   | 0（不限）                | 单作业核数上限                        |
| `CREST_CORES_PER_JOB`                       | 4                    | fixed 模式下每作业核数                 |
| `CREST_MTD_BOOST`                           | 1.25                 | 每轮 MTD 的 k 放大倍数，`=1` 禁用        |
| `CREST_MTD_BOOST_MAX`                       | 4.0                  | 累计放大上限                         |
| `CREST_IRU_KSCALE`                          | 1.0                  | `-iru` 静态偏置的 k 缩放              |
| `CREST_LOADCONF_KSCALE`                     | 回退上一项                | `--loadconformer` 专用缩放         |
| `CREST_IRU_ROTAMERS`                        | 开                    | `=0` 改读 `crest_conformers.xyz` |
| `CREST_TRIAL_MD_PS`                         | ORCA 0.125 / 其他 0.25 | trial MTD 时长（ps）               |
| `CREST_TRIAL_VBDUMP_FS`                     | 50                   | trial Vbias dump 间隔（fs）        |
| `CREST_FLEXI_VALCAP` / `CREST_FLEXI_PMEAN`  | 0.9 / 1.5            | 柔性单键上限 / power mean 指数         |
| `CREST_FLEXI_FORCEFAC`                      | 1.3                  | 声明键的距离门控                       |
| `CREST_SCRATCH_DIR`                         | 未设                   | 高频中间文件目录                       |
| `CREST_WORK_DIR`                            | `crest_work`         | 低频中间文件目录，`none` 回上游            |
| `CREST_STOP_CHECK_SEC` / `CREST_STOP_FILE`  | 300 / `stop`         | stop 文件轮询间隔、文件名                |
| `CREST_EXCL_XYZ`                            | 关                    | 导出被 CREGEN 剔除的结构               |
| `CREST_SKIP_TOPO`                           | 开                    | 拓扑不匹配时保留结构                     |
| `CREST_GXTB_SOLV`                           | 空                    | g-xTB 的 `--gbe` / `--cosmo` 透传 |
| `CREST_ORCA_MPI_ARGS` / `CREST_ORCA_NPROCS` | 空 / 未设               | ORCA 额外 MPI 参数 / 固定 rank 数     |

一个 48 核节点的常见配置：

```bash
source ./crest_env.sh 48
export CREST_SCRATCH_DIR=/tmp/crest-$SLURM_JOB_ID
export OPENBLAS_NUM_THREADS=1
export CREST_MAX_CORES_PER_JOB=24
export CREST_GXTB_SOLV='--gbe thf'
crest molecule.xyz -T 48 -gxtb_dev
```

## 四、升级后需要知道的差异

1. **轨迹后缀变成 `.xyz`**：旧脚本如果硬编码 `.trj` 或 `crestopt.log`，需要改成 `.xyz`。
2. **GFN2 走外部 xtb**：PATH 上必须有 xtb 可执行文件；`-gfn1` 仍走进程内 tblite，精简构建没编 tblite 时不可用。
3. **trial MTD 更短**：ORCA 默认 0.125 ps，其余 0.25 ps。省时间，但时间预估会更抖，而且只在大约 50 步之后才暴露的不稳定性可能漏检。需要保守时用 `CREST_TRIAL_MD_PS=0.5` 或更长。
4. **`--cinp` 默认是刚性约束**：确实需要旧的谐振软约束，用 `-softcons`，或 TOML `hardconstraints = false`。
5. **Vbias dump 默认值变了**：生产 MTD 默认 500 fs，trial 默认 50 fs。旧任务依赖原来的 1000 fs 间隔，请显式设 `-vbdump` / `CREST_TRIAL_VBDUMP_FS`。
6. **动态核数再分配默认开启**：`CREST_DYNAMIC_SHARE=0` 可关。ORCA/MPI 不一定核越多越快，小体系建议 `CREST_MAX_CORES_PER_JOB=24` 一类的上限。
7. **`stop` 文件只对动力学阶段生效**：多层优化和 CREGEN 不检查；已经进入这些阶段后，`touch stop` 不会立刻结束整个任务。
8. **柔性与拓扑判定改过**：柔性公式、金属键的拓扑排除、`CREST_SKIP_TOPO` 都会影响构象集合和 MTD 时长。要复现上游结果，用 `CREST_FLEXI_VALCAP=1.0` + `CREST_FLEXI_PMEAN=2.0`、`-nometaltopo`、`CREST_SKIP_TOPO=0`。

## 五、适用范围与引用

- 适用于 Linux 集群上的 CREST 3.0.2 工作流，尤其是外部 xtb / g-xTB / ORCA 计算器。
- 这些改动主要是稳定性、性能、诊断和文件管理改进；算法未动的部分，科学结果应与上游一致。SHAKE 的迭代格式与预算、柔性公式、柔性和拓扑的默认行为属于数值/流程改动，涉及强约束或金属体系时建议用自己的测试体系复核。
- 如果用于发表，请按上游 `CITATION.cff` 和 CREST 论文引用 CREST；本地修改不改变上游 authorship。
- 许可证：LGPL-3.0-or-later。详见 `COPYING`、`COPYING.LESSER`、`LICENSE`、`NOTICE`。

---

# CREST 3.0.2 · QW-Zhang-PKU modified (@ Yu group, PKU) — English

*This is the English translation of the Chinese part above (translated by AI).*

This directory contains a modified version of CREST 3.0.2, made by QW-Zhang. The project aims to fix a number of long-standing pain points in CREST and make it better suited to mechanistic studies in organic chemistry, especially organometallic chemistry. Hy4 and Deepseek were used while modifying the code; the resulting features have been tested and are described below. Most of this document was written by the author himself.

Contact: <2201110483@stu.pku.edu.cn>

File-by-file technical details are in `LOCAL-CHANGES.md`; licensing and third-party notices are in `NOTICE`, `LICENSE`, `COPYING` and `COPYING.LESSER`.

## 1. Build and Run

### Build

All subprojects are vendored in the source tree, so a completely offline build with the GNU toolchain is possible:

```bash
bash build-now.sh                               # installs to $HOME/crest-dist/bin/crest
```

`build-now.sh` runs a self-check and then delegates to `build_crest.sh` for the CMake + Ninja build. The default slim build turns off tblite / GFN-FF / libpvol / lwoniom; see item 5 in section 2 for the reason. Use `FULL=1` to restore the full upstream feature set.

If the build machine and the compute nodes have different BLAS environments, you can build a static OpenBLAS first:

```bash
bash build-openblas-static.sh
OPENBLAS_DIR=$HOME/opt/openblas-static STATIC_BLAS=1 bash build-now.sh
```

Full build notes are in `README-BUILD.md`.

### Quick start

```bash
source ./crest_env.sh 40
crest molecule.xyz -T 40 --quick
```

The calculation levels introduced here are gfn2-xtb, g-xTB and ORCA.

**GFN2** (which runs through an external xtb subprocess in this version) can be selected directly on the command line; make sure the xtb executable is in `PATH` and that `XTBPATH=` points at the xtb executable:

```bash
crest molecule.xyz -T 40 -gxtb_dev --quick
```

Or use a TOML file as the input:

```toml
input = "molecule.xyz"

[calculation.level]
method = "gfn2"
binary = "/PATH/TO/xtb/bin/xtb"   # optional; use xtb from PATH
flags  = "--gfn2"
```

**g-xTB** can likewise be selected directly on the command line; make sure the xtb executable belonging to g-xTB is in `PATH`:

```bash
crest molecule.xyz -T 40 -gxtb_dev --quick
```

**ORCA** is specified in TOML:

```toml
input = "molecule.xyz"

[[calculation.level]]
method = "orca"                          # sampling level

[[calculation.level]]
method        = "orca"
orca_cmd      = "/PATH/TO/orca"   # change to the real ORCA location
orca_template = "orca.inp"               # ORCA input template
chrg          = 0                        # change to the real charge
multiplicity  = 1                        # change to the real multiplicity
refine        = "sp"
```

`orca_template` is an ORCA input template you write yourself. CREST appends `EnGrad` on every call, so you do not need to specify `%pal` in the template; CREST fills it in to match the actual run.  
More configuration (multiple levels, MPI arguments, constraint-file syntax) is in `ORCA-HOWTO.md`.  
ORCA slows down on a mechanical disk, an NFS mount or a cluster's distributed storage, so give each job its own scratch directory (preferably on an SSD):

```bash
export CREST_SCRATCH_DIR=/tmp/crest-$SLURM_JOB_ID
```

An example `orca.inp`:

```inp
! PBE D3 MINIX GCP(DFT/MINIS) CPCM(THF) looseSCF nopop miniprint

%maxcore 2000

```

## 2. New Features (vs. CREST 3.0.2, 11+ items)

### 2.1 g-xTB and ORCA interfaces

**Reason:** CREST ships with methods such as gfn-xtb and gfn2-xtb. Before this project I mostly used gfn2-xtb. But gfn2-xtb describes transition-metal complexes poorly and sometimes produces completely wrong results (for example, a CO coordinated to the metal undergoing migratory insertion), so higher-level methods urgently needed to be integrated into CREST.

**g-xTB:** g-xTB is a new-generation semiempirical method developed by the Grimme group. It is considerably more accurate than gfn2-xtb, especially for transition-metal complexes, and costs only modestly more compute time, so it is the method I now use most often in CREST. Note that g-xTB has so far only appeared as a preprint and may still change. Its solvation support is not good either, so setting a solvent in g-xTB mode is not recommended; leave it in vacuum by default.

**ORCA:** Sometimes a system is stubborn enough that DFT is required as the calculation level, however expensive; ORCA can be used then. We completed and enabled CREST's ORCA interface. Tests show that a pure functional with a minimal basis (e.g. `! PBE D3 MINIX GCP(DFT/MINIS) LooseSCF`) is affordable, but it takes a very long time and is only recommended for difficult cases. Only ORCA 6.1.1 has been tested so far; older versions should be compatible. Note that ORCA's **temporary directory** should be placed on an **SSD**, otherwise slow storage (a mechanical disk or a cluster's distributed file system) will drag the run down.

**Changes:**

- Use `-gxtb_dev` to enable g-xTB.
- You need to download the g-xTB executable (<https://github.com/grimme-lab/g-xtb>) and point `XTBPATH=` at its `share/xtb` directory.
- If you need g-xTB solvation, pass it through an environment variable, e.g. `export CREST_GXTB_SOLV='--gbe thf'`. This is for debugging only; turning it on usually makes the job fail.
- ORCA is connected through a TOML file; see the usage section above.
- `CREST_ORCA_PRINT_CALL=1` prints how CREST invokes ORCA, for debugging.

### 2.2 `-cinp` constraints become rigid

**Reason:** When searching transition-state conformers, some atoms have to be locked (for example the two atoms forming a C–C bond), so `-cinp` is used to lock them. Upstream treats `-cinp` constraints with a harmonic model whose strength defaults to 0.5. As a result the transition state very easily relaxes back to product or reactant during the conformer search (the relaxation back to reactant/product is strongly exothermic and overcomes the spring), and no useful result is obtained. Raising the harmonic strength to 5 does lock the reacting atoms, but such a stiff spring makes MTD fail to converge and forces the time step down to 1–2 fs. This is a case of filling one hole while digging another, because the total runtime then increases sharply.

**SHAKE:** `-cinp` now uses rigid SHAKE constraints, so setting the harmonic strength to a large value such as 5 or 6 is much less likely to cause an error.

**Changes:**

- The SHAKE correction sweep was changed from Jacobi to **Gauss–Seidel** (`src/dynamics/shake_module.f90:511-558`).
- Maximum sweeps `maxcyc = 6000` (upstream 250; one real system needed more than 3000 sweeps), convergence tolerance `tolshake = 1e-5` (upstream 1e-7).
- The failure report prints the iteration count/limit, the number of violated bonds, the deviation of the raw step, the deviation after correction, the tolerance, etc. in `es12.4` format for debugging.
- The switch back to the harmonic model is kept: `-softcons`.
- Bonds generated pairwise by `atoms:` / `elements:` are registered as rigid too. This syntax expands **every pair** in the atom list into a bond (n atoms give C(n,2) bonds); once the set exceeds 3n−6 degrees of freedom it becomes linearly dependent and the code prints a warning.
- Pre-optimization and batch optimization are **not** switched to rigid constraints (this could be improved later, but it requires changing the optimizer and is not urgent, so it was left alone). Setting a very large force constant (>5) can lock them almost completely (the locked atoms still move very slightly, but the displacement is negligible).
- Bonds generated automatically by `-cbonds`, `-cmetal`, `-cheavy` and `-clight` are still **harmonic constraints** and do not lock the whole molecule. `-cmetal` is a very useful option: it prevents ligands from dissociating from the metal center.

### 2.3 Making multiple MTD rounds useful (`-iru` etc.)

**Reason:** The following applies only to the organic-mechanism work I do. In most jobs CREST runs at least 2 rounds of MTD sampling. The idea seems to be: the first round finds a lower-energy conformer, the second round starts from it, and if an even lower one is found a third round follows, and so on. But I found that for the vast majority of my jobs (a single transition-metal center plus common ligands) one MTD round already finds essentially all conformers. The second MTD is usually wasted work (in more than 90% of cases) and only acts as insurance. In that situation you are effectively **running the same N ps MTD twice with identical parameters**, which is not very meaningful and is even worse than running one 2N ps MTD.

**`-iru`:** On one hand, a mode that runs only one MTD round was added (`-onebatch `, see below). On the other hand, `-iru` was added so that the second MTD can read the minimum regions already explored in the first round and apply a bias to them, leaving the "known space" of the first round faster and giving more time to "unexplored space", so the second MTD is put to good use. In addition, the readconformers feature was added so that conformers explored in a previous job can be read across jobs, for cases where you want to run CREST a second time on the same system.  
We also made the bias strength k grow round by round, so the exploration strength increases every round and the run does not keep circling in known regions.

**Changes:**

- The k of the `-iru` static bias is set through the environment variable `CREST_IRU_KSCALE=`, default **1.0**.
- k is multiplied by **1.25** each MTD round (`CREST_MTD_BOOST`), with a cumulative cap of **4.0** (`CREST_MTD_BOOST_MAX`), applied from the 2nd round on to the previous value. `CREST_MTD_BOOST=1` turns it off. This boost also applies to the `-iru` static potential.
- The cross-run conformer reading is `--loadconformer` (alias `-loadconf`, or `loadconformer` in TOML): it reads known conformers from `readconformers.xyz` in the **current working directory**, using exactly the same mechanism as `-iru`.

### 2.4 A better flexibility measure

**Reason:** The flexibility value flex determines the MTD length: `t(MTD) = 3.0·exp(0.10·flex·max(1,rednat−8))` ps, with a floor of 5 ps and a ceiling of 500 ps. The MTD length in turn affects the length of the whole job.

When searching conformers for a transition state, setting `-cinp ` tends to inflate the flexibility.  
Specifically, if the reaction is an intramolecular cyclization and I use -cinp to lock the two atoms A and B that are about to form the ring, CREST sees that the A–B distance is not yet a single bond, fails to recognize a ring, and treats the molecule as a chain.

This makes CRSET think the molecule is very flexible and needs a long MTD, when in fact it is a cyclic molecule with few conformers and does not need such a long MTD.  
In addition, the upstream flexibility value jumps as bond lengths change (two structures that look almost identical can get very different flexibilities), so the formula was improved.

The new version makes three improvements:

- **Continuous ring detection**: a sigmoid `w = 1/(1 + exp(−16·(rco/r − 1)))` gives every bond a weight between 0 and 1, then Dijkstra on `−log w` finds the strongest closed path that avoids the i–j bond, giving the ring size L and a continuous membership ρ (`RHOFULL = 0.98`), with `ringf = 1 − (1 − 0.5·(1 − exp(−0.06·L)))·ρ`. The former black-or-white criterion becomes an interpolation, so a small change in bond length no longer makes the flexibility jump. `CREST_FLEXI_DEBUG=1` prints `r/rcov`, whether the bond was forced, ρ, etc. for every bond.
- **Let the bonds declared by `-cinp` take part in ring closure**: these bonds are inserted into the bond graph even when longer than the sum of covalent radii (an `ffac = 1.3` gate filters out bogus bonds produced by `atoms:` expansion) and are given full weight `w = 1` in the ring weight. A long "about to close" bond can then close the ring and the flexibility is no longer inflated. Note that this does not discount the locked atoms; the ring is recognized correctly and `ringf` lowers the flexibility automatically.
- **Cap each bond and change the aggregation**: `val = min(val_raw, 0.9)`, and the root mean square is replaced by a power mean `flex = (Σ val^p / m)^(1/p)` with `p = 1.5`. The cap removes the dominance of a few anomalous bonds, and p < 2 stops the result from being pulled around by large values. Set `CREST_FLEXI_VALCAP=1.0` and `CREST_FLEXI_PMEAN=2.0` to restore the upstream formula for an A/B comparison.
- If the resulting MTD length is still too long (or too short), `-mdlen` multiplies it by a coefficient. In my own work I use `-mdlen=0.8`.

### 2.5 Fixing thread oversubscription in CREST 3.0.2 (a very serious problem)

**Reason:** Upstream CREST 3.0.2 runs xtb calculations through the built-in tblite. We found that even with `-T 40`, the MTD stage spawns more than 80 processes that compete with each other, making the run unacceptably slow.  
You can see it in the difference between the MTD time estimated by the Trail MTD and the actual time: the actual MTD time is often more than 10× the estimate, and most of that time is wasted by the bug, doing nothing at all (pure waste).  
The root cause is that tblite opens its own OpenMP team **inside** the CREST process; `-T` controls external xtb/ORCA, but not that team.

**Changes:**

- We dropped the in-process tblite and call an external xtb program instead (as in the old CREST 2.x versions), which avoids the problem. After the fix the MTD runs at normal speed and time is no longer wasted.
- `gfn2` / `--gfn2` was redirected from `jobtype%tblite` to `jobtype%xtbsys`, i.e. it calls an external xtb subprocess.
- gfn1 was not redirected because it is used less often.
- The build scripts default to `WITH_TBLITE=OFF`, so in-process GFN2 is not compiled at all and the nested thread source is removed at its origin; `WITH_TBLITE=ON` turns it back on when needed.

**The descriptions of the other bug fixes below were written by AI:**  
Thread-budget propagation was fixed at the same time: `ompenvset()` exports `OMP_NUM_THREADS` / `MKL_NUM_THREADS` / `OPENBLAS_NUM_THREADS` as `cores_per_job` and fixes `OMP_MAX_ACTIVE_LEVELS=1` to prevent "threads inside threads". The original Fortran binding for POSIX `setenv()` declared one argument too few (the third `overwrite` argument), read an undefined register, and when the variable already existed `setenv` silently refused to update it, so children inherited the calling shell's value instead. The binding is now three-argument with an explicit 1 (the measurement recorded in patch `crest-omp-env-fix.patch`: 40 threads, 6 MTDs, from 2m48s back to close to the estimated 44s).

Another bug with a very similar symptom was fixed too: `xtb_sc` / `orca_sc` / `turbom_sc` cached the **full path** of the gradient file in `calc%gradfile`, while parallel MTD gives every worker a different `calcspace` suffix, so all workers read the stale trial-MTD gradient, SHAKE failed after a few integration steps, and the run reported `*MTD n terminated EARLY`. The path is now re-derived from the current `calcspace` on every call (see `patches/crest-xtbsys-mtd-fix.patch`).

### 2.6 Dynamic core allocation

**Reason:** Take MTD as an example: I use 36 cores to run 6 MTDs at once, 6 cores each. If one MTD finishes early, those 6 cores sit completely idle.

Now the remaining 5 jobs split all 36 cores, 7 cores each, with only 1 core idle. This raises overall CPU utilization and shortens the job.

For semiempirical xtb jobs the gain exists but is not especially large. For large systems with ORCA it can cut the runtime noticeably (sometimes by several hours).

**Changes:** Dynamic core sharing is now on by default. The output file shows `[cores] dynamic sharing: ...` only once, when dynamic sharing is triggered for the first time; it is not printed afterwards because there is no useful information in repeating it.

### 2.7 Fixing the bug that discarded structures

- By default, bonds involving metal atoms are excluded from the topology check, so metal systems lose fewer conformers because of jumps in the rounded coordination number.
- `cregen_topocheck()` first builds a `topoexcl(:)` mask (atoms given by `-notopo` plus all metals) and then compares topologies with `bondtotopo_excl()` (`src/cregen.f90:896-930`, implemented in `src/ztopology.f90:339`). Metals are identified with an explicit 86-element table, with Li and Be additionally excluded (`src/dynamics/shake_module.f90:27-32,188-204`), not simply by Z > 20. `excludeTOPOmetals` is on by default; `-nometaltopo` or TOML `metaltopo = false` turns it off.
- `CREST_SKIP_TOPO` is on by default: a structure whose topology does not match is now **kept** (and tagged `SKIP`) instead of being discarded; `=0` / `false` / `off` restores the upstream behavior. This default was changed after testing on a 108-atom Rh system.
- The EWIN energy zero is now the "lowest topologically **valid** structure" instead of the global minimum (`cregen.f90:1290`), so one structure with a broken topology but a very low energy cannot discard the whole batch.
- `crest_refine_topopresort()` was added to the refine stage (`src/algos/refine.f90:176-229`): a cheap topology prescreen runs before the expensive single points. It screens topology only, with no energy sorting or RMSD deduplication.
- Rejected structures can be exported for auditing: `CREST_EXCL_XYZ=1` writes `crest_excluded_structures.xyz`, with a `CLASH` / `TOPO` / `EWIN` / `DUP` / `SKIP` tag and reason in each frame's comment line (`src/exclfile.f90`, off by default). This is faster than digging through the log to find out "where did this conformer go".
- Several performance changes were made to the deduplication loop in `cregen_CRE`: masked coordinates are cached to avoid repeated O(nall²) inner computations, the energy-sorted j window is narrowed, and the heavy-atom subset is extracted only once per structure. The RMSD / energy thresholds themselves were not changed.
- The new `src/minregistry.f90` is a registry of located minima: within one optimization stage, if a running optimization is judged to be heading toward an already registered minimum, it is stopped early and adopts that minimum's geometry and energy. It is cleared at the start of every optimization stage (so that a later stage with tighter thresholds is not pruned by the previous stage's results); `CREST_PRUNE_ANCHORS=0` turns it off.

### 2.8 Better output

- MD/MTD progress is reported at fixed **25% / 50% / 75%** marks.
- Progress lines are shorter and fit narrow terminals; the multiple output lines of parallel jobs are less likely to wrap out of alignment.
- Elapsed time and estimated remaining time are shown (`elapsed ... remaining ...`).
- Trajectories and optimization logs use common filename suffixes:
  - `crest_dynamics.trj` → `crest_dynamics.xyz`
  - `crest_<n>.trj` → `crest_<n>.xyz`
  - `crestopt.log` → `crestopt.xyz`
- The initial geometry optimization reports the number of cycles and the wall time, for example:

```text
Initial geometry optimization: 37 cycles, 12 sec   (crestopt.xyz)
```

- In upstream CREST the parameter blocks printed at MTD and MD startup were out of order; this has been fixed.
- `optimization starting` is printed before the ORCA initial optimization, so the log is not completely silent for a long time (the pre-optimization can sometimes take tens of minutes).

### 2.9 Signals and the `stop` file

**Reason:** (Especially when ORCA is the calculation engine,) sometimes the dynamics has not yet reached the preset MTD length but enough structures have already been sampled (for example, one phosphorus of a BINAP chelating ligand has decoordinated and swung away), so continuing the MTD is pointless.

A way was added to stop the MTD at any time and move on to the next stage (batch optimization).

- `stop` file: at the top of every MD step in the dynamics stage the file is checked on a wall-clock interval (one `inquire` every 300 s by default), and only a `stop` file **created during the run** counts, so a stale file left over from before will not trigger it. On a hit the current dynamics ends cleanly with `term = 3`, the trajectory is flushed frame by frame, the structures collected so far are handed to the following optimization/filtering, and the stop file is deleted automatically after the whole parallel round finishes.
- The variables are `CREST_STOP_FILE` for the name and `CREST_STOP_CHECK_SEC` for the interval (`0` disables it).
- Trial MTD does not use this mechanism, and the multilevel optimization and CREGEN stages do not check it either; once those stages have started, `touch stop` will not end the run.

Some other changes:

- Ctrl+C, `kill` and SLURM timeout signals now really terminate CREST and its external child processes; the terminal no longer appears frozen until you resort to `kill -9`.

### 2.10 Cleaner working directories

I run crest on a cluster, and upstream CREST gives you no way to change the working directory. My files sit on a mechanical disk on an NFS-mounted LAN, which slows file access down.  
Now high-frequency intermediate files (`MDFILES`, `OPTIM`, the `calculation.level.*` calculation directories) can be placed on a node-local temporary directory, while low-frequency ones (`coord`, `wbo`, `crestopt.xyz`, etc.) go into `crest_work/` under the working directory by default:

```bash
export CREST_SCRATCH_DIR=/tmp/crest-$SLURM_JOB_ID   # hot files, SSD recommended
export CREST_WORK_DIR=crest_work                    # cold files, default
export CREST_WORK_DIR=none                          # back to the upstream layout
```

This is implemented with symlinks; none of the `open()` call sites were changed. Result files such as `crest_conformers.xyz`, `crest_rotamers.xyz` and `crest_best.xyz` still stay in the working directory. `CREST_STAGE_FILES=0` turns the whole mechanism off.

### 2.11 A cheaper sampling workflow

- Trial MTD is shorter: 0.125 ps by default for ORCA, 0.25 ps for other calculators; `CREST_TRIAL_MD_PS` overrides it.
- The trial-MTD Vbias dump is fixed at 50 fs by default, with at least two updates guaranteed.
- Newly added / wired up:
  - `--onebatch`: run only one MTD round. It still includes trial MTD + MTD + MD + crosscof.
  - `--loadconformer`: read known conformers from `readconformers.xyz` and apply a static RMSD bias from the first MTD round (see above).
  - `-iru`: add the conformers found so far as a second static bias on top of the dynamic MTD (see above).

## 3. Environment Variables

| Variable | Default | Purpose |
| --- | --- | --- |
| `CREST_DYNAMIC_SHARE` | on | `=0` disables dynamic core redistribution |
| `CREST_MAX_CORES_PER_JOB` | 0 (unlimited) | cap on cores per job |
| `CREST_CORES_PER_JOB` | 4 | cores per job in fixed mode |
| `CREST_MTD_BOOST` | 1.25 | per-round multiplier for MTD k; `=1` disables it |
| `CREST_MTD_BOOST_MAX` | 4.0 | cumulative boost cap |
| `CREST_IRU_KSCALE` | 1.0 | k scaling of the `-iru` static bias |
| `CREST_LOADCONF_KSCALE` | falls back to the previous item | scaling specific to `--loadconformer` |
| `CREST_IRU_ROTAMERS` | on | `=0` reads `crest_conformers.xyz` instead |
| `CREST_TRIAL_MD_PS` | ORCA 0.125 / others 0.25 | trial MTD length (ps) |
| `CREST_TRIAL_VBDUMP_FS` | 50 | trial Vbias dump interval (fs) |
| `CREST_FLEXI_VALCAP` / `CREST_FLEXI_PMEAN` | 0.9 / 1.5 | per-bond flexibility cap / power-mean exponent |
| `CREST_FLEXI_FORCEFAC` | 1.3 | distance gate for declared bonds |
| `CREST_SCRATCH_DIR` | unset | directory for high-frequency intermediate files |
| `CREST_WORK_DIR` | `crest_work` | directory for low-frequency intermediate files; `none` restores upstream |
| `CREST_STOP_CHECK_SEC` / `CREST_STOP_FILE` | 300 / `stop` | stop-file polling interval and filename |
| `CREST_EXCL_XYZ` | off | export structures rejected by CREGEN |
| `CREST_SKIP_TOPO` | on | keep structures whose topology does not match |
| `CREST_GXTB_SOLV` | empty | pass through g-xTB's `--gbe` / `--cosmo` |
| `CREST_ORCA_MPI_ARGS` / `CREST_ORCA_NPROCS` | empty / unset | extra ORCA MPI arguments / fixed rank count |

A typical configuration on a 48-core node:

```bash
source ./crest_env.sh 48
export CREST_SCRATCH_DIR=/tmp/crest-$SLURM_JOB_ID
export OPENBLAS_NUM_THREADS=1
export CREST_MAX_CORES_PER_JOB=24
export CREST_GXTB_SOLV='--gbe thf'
crest molecule.xyz -T 48 -gxtb_dev
```

## 4. Compatibility Notes

1. **Trajectory suffixes changed to `.xyz`**: scripts that hard-code `.trj` or `crestopt.log` need to be updated.
2. **GFN2 goes through an external xtb**: an xtb executable must be on `PATH`; `-gfn1` still uses in-process tblite and is unavailable in a slim build compiled without tblite.
3. **Trial MTD is shorter**: 0.125 ps for ORCA, 0.25 ps otherwise. This saves time, but runtime estimates are noisier and instabilities that only show up after about 50 steps may be missed. Use `CREST_TRIAL_MD_PS=0.5` or longer for conservative behavior.
4. **`--cinp` is rigid by default**: if you really need the old harmonic soft constraints, use `-softcons`, or set `hardconstraints = false` in TOML.
5. **Vbias dump defaults changed**: 500 fs for production MTD, 50 fs for trial. If an old job depended on the original 1000 fs interval, set `-vbdump` / `CREST_TRIAL_VBDUMP_FS` explicitly.
6. **Dynamic core redistribution is on by default**: `CREST_DYNAMIC_SHARE=0` turns it off. More ORCA/MPI ranks are not necessarily faster; for small systems an upper bound such as `CREST_MAX_CORES_PER_JOB=24` is recommended.
7. **The `stop` file only affects the dynamics stage**: multilevel optimization and CREGEN do not check it; once those stages have started, `touch stop` will not end the job immediately.
8. **Flexibility and topology decisions have changed**: the flexibility formula, the exclusion of metal bonds from topology checks and `CREST_SKIP_TOPO` all affect the conformer set and the MTD length. To reproduce upstream results, use `CREST_FLEXI_VALCAP=1.0` + `CREST_FLEXI_PMEAN=2.0`, `-nometaltopo` and `CREST_SKIP_TOPO=0`.

## 5. Scope, Citation, and License

- Intended for CREST 3.0.2 workflows on Linux clusters, especially with external xtb / g-xTB / ORCA calculators.
- These changes are mainly stability, performance, diagnostics and file-management improvements; where the algorithm was not touched, scientific results should agree with upstream. The SHAKE iteration scheme and budget, the flexibility formula, and the new flexibility/topology defaults are numerical/workflow changes; if you work with strong constraints or metal-containing systems, validate them on your own test systems.
- For publications, cite CREST according to the upstream `CITATION.cff` and the CREST papers; local modifications do not change upstream authorship.
- License: LGPL-3.0-or-later. See `COPYING`, `COPYING.LESSER`, `LICENSE` and `NOTICE`.
