# 🧬 SeqForge — 零依赖生物序列智能工具箱

**简体中文** | [繁體中文](README.zh-TW.md) | [English](README.en.md)

[![CI](https://github.com/gitstq/seqforge/actions/workflows/ci.yml/badge.svg)](https://github.com/gitstq/seqforge/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python](https://img.shields.io/badge/Python-3.8%2B-blue.svg)](https://www.python.org/downloads/)
[![Dependencies](https://img.shields.io/badge/dependencies-0-brightgreen.svg)](https://github.com/gitstq/seqforge)

> 🔥 **一个文件、一条命令，搞定 FASTA/FASTQ 的日常分析** —— 纯 Python 标准库实现，无需安装任何第三方包，Windows / macOS / Linux 全平台开箱即用。

---

## 🎉 项目介绍

**SeqForge** 是一个轻量级的生物信息学序列分析工具箱。做序列分析的你是否遇到过这些烦恼：

- 😫 `biopython` 功能强大，但依赖安装繁琐，版本冲突频发
- 😤 `seqkit` 性能出色，却需要下载二进制，内网 / 老旧 Windows 环境处处受限
- 🤯 只想算个 GC 含量、找个 ORF、翻译一下序列，却要打开沉重的 IDE 或写几十行脚本

SeqForge 就是为终结这些烦恼而生的 **100% 自研** 项目：整个工具仅依赖 Python 标准库，一个 `pip install`（甚至不安装）即可运行，把最常用的序列分析能力浓缩进 16 个清晰的子命令。

**自研差异化亮点**：

- 🪶 **零依赖**：不引用任何第三方包，克隆即可运行，彻底告别"环境地狱"
- 🌏 **中文编码友好**：自动识别 UTF-8 / GB18030 / Latin-1 编码，国内实验室的老数据不再乱码
- 📦 **`.gz` 透明支持**：压缩格式的 FASTA/FASTQ 直接读取，无需手动解压
- 🤖 **CI 友好**：所有命令支持 `--json` 机器可读输出，`validate` 提供 0/2 退出码，可直接作为流水线质量门禁
- 🧪 **云端质量门禁**：GitHub Actions 在 Windows / macOS / Linux × Python 3.8–3.12 矩阵上执行语法检查与全量单元测试

**灵感来源**：GitHub Trending 持续火热的"本地优先 / 零依赖数据处理"工具方向（浏览器端大数据表格、本地索引引擎等高热项目），叠加生物信息社区对"轻量、免安装、跨平台"序列工具的长期诉求。

---

## ✨ 核心特性

| 能力 | 说明 |
| --- | --- |
| 📊 **统计引擎** | 序列长度分布、**N50 / L50**、GC 含量、**GC 偏斜**、滑窗 GC 图谱 |
| ✅ **格式校验** | 结构校验、非法字符检测、重复 ID 检查、空序列检查，编码自动嗅探 |
| 🧬 **遗传密码** | 内置 **5 张 NCBI 遗传密码表**（标准 / 脊椎动物线粒体 / 酵母线粒体 / 无脊椎动物线粒体 / 纤毛虫），支持 **六框翻译** |
| 🔎 **ORF 查找** | 双链六框扫描，IUPAC 兼容，自动坐标回填到原始链 |
| 🎯 **基序搜索** | IUPAC 歧义码匹配（如 `GASTY`），支持双链检索 |
| ✂️ **酶切位点** | 内置 **20 种常用限制性内切酶**（EcoRI、BamHI、NotI……） |
| 🌡️ **引物 Tm** | Wallace 法则 + 盐浓度校正公式，一键评估退火温度 |
| 🧹 **FASTQ 质控** | Phred+33/+64 偏移自动检测、**BWA 式滑窗修剪**、长度 / 质量 / N 含量过滤 |
| 🔄 **格式互转** | FASTA / FASTQ / CSV / TSV / NDJSON 自由转换，FASTA→FASTQ 一键补全 |
| 🗜️ **压缩支持** | `.gz` 文件透明读写 |
| 🎲 **可复现工具** | 带随机种子的抽样、大小写不敏感的去重 |
| 🖥️ **跨平台** | 纯 Python 3.8+，零编译、零运行时依赖 |

---

## 🚀 快速开始

### 环境要求

- 🐍 **Python 3.8 或更高版本**（仅此一项，无任何第三方依赖）

### 方式一：pip 安装（推荐）

```bash
pip install seqforge
seqforge --version
```

### 方式二：免安装，直接运行

```bash
git clone https://github.com/gitstq/seqforge.git
cd seqforge
python -m seqforge --version
```

### 方式三：作为库引入

```python
from seqforge.seqio import read_fasta
from seqforge.genetic import translate

records = read_fasta("demo.fasta")
for record in records:
    print(record.id, len(record.seq), translate(record.seq))
```

### 30 秒上手

```bash
# 1. 准备一份示例数据
cat > demo.fasta <<'EOF'
>gene_demo 示例序列
ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG
>gene_short
GGCCGG
EOF

# 2. 查看统计摘要（N50、GC、长度分布一目了然）
seqforge stats demo.fasta

# 3. 六框翻译，寻找潜在编码序列
seqforge translate demo.fasta --frames 6
```

---

## 📖 详细使用指南

### 命令总览

| 命令 | 功能 |
| --- | --- |
| `stats` | 长度 / GC / N50 / L50 摘要统计 |
| `validate` | FASTA/FASTQ 结构校验（CI 门禁） |
| `gc` | 逐序列 GC 含量、GC 偏斜、滑窗 GC |
| `revcomp` | 反向互补（DNA/RNA） |
| `transcribe` | DNA → RNA 转录 |
| `translate` | 翻译为蛋白（1/3/6 框、5 种密码表） |
| `orf` | 开放阅读框查找（双链） |
| `find` | IUPAC 基序检索 |
| `enzymes` | 限制性酶切位点分析 |
| `tm` | 引物退火温度计算 |
| `fqstat` | FASTQ 质量统计 |
| `fqtrim` | FASTQ 修剪与过滤 |
| `convert` | 格式互转（fasta/fastq/csv/tsv/ndjson） |
| `sample` | 可复现随机抽样 |
| `dedup` | 重复序列去除 |
| `info` | 格式探测与记录预览 |

> 💡 所有命令均支持 `-o/--output` 写入文件、`--encoding` 强制指定编码；统计 / 分析类命令支持 `--json` 机器可读输出。

### 📊 示例一：数据质控（`stats` + `validate`）

```bash
$ seqforge stats demo.fasta
+---------------+---------+
|    metric     |  value  |
+---------------+---------+
| format        | fasta   |
| records       | 2       |
| total bases   | 45      |
| min length    | 6       |
| max length    | 39      |
| mean length   | 22.5    |
| median length | 22.5    |
| length stdev  | 16.5    |
| N50           | 39      |
| L50           | 1       |
| GC content(%) | 64.44   |
+---------------+---------+

$ seqforge validate demo.fasta
# 校验通过时退出码为 0；发现错误时退出码为 2，可直接用于 CI 门禁
```

### 🧬 示例二：ORF 查找与翻译

```bash
$ seqforge orf gene.fasta --min-len 100
+-----------+--------+-------+-------+-----+--------+----------+
|    id     | strand | frame | start | end | length |  protein |
+-----------+--------+-------+-------+-----+--------+----------+
| gene_demo |   +    |   1   |   1   | 39  |   39   | MAIVNGL...|
+-----------+--------+-------+-------+-----+--------+----------+

$ seqforge translate gene.fasta --table 2 --frames 6
# --table 2 使用脊椎动物线粒体密码表；--frames 6 同时翻译两条链
```

### 🎯 示例三：TATA 框基序检索

```bash
$ seqforge find promoter.fasta TATAAA --both-strands
+----------+-------+-----+--------+
|    id    | start | end | strand |
+----------+-------+-----+--------+
| promoter |   53  | 59  |   +    |
| promoter |  120  | 126 |   -    |
+----------+-------+-----+--------+
```

### ✂️ 示例四：FASTQ 质控流水线

```bash
# 查看质量概况（自动识别 Phred+33/+64）
$ seqforge fqstat reads.fastq

# 滑窗修剪 + 过滤，一步到位
$ seqforge fqtrim reads.fastq \
    --quality 20 --window 4 \
    --min-len 50 --max-n 0.1 \
    -o reads.clean.fastq
```

### 🌡️ 示例五：引物 Tm 计算

```bash
$ seqforge tm ATGCGTATGCGTATGCGT --na 0.05
+---------------------+--------+--------+----------------+-----------------+
|       primer        | length | gc(%)  | Tm Wallace(C)  | Tm adjusted(C)  |
+---------------------+--------+--------+----------------+-----------------+
| ATGCGTATGCGTATGCGT  |   18   | 50.0   |       84       |      42.9       |
+---------------------+--------+--------+----------------+-----------------+
```

### 🤖 示例六：CI 中的 JSON 门禁

```bash
# 输出机器可读 JSON，配合 jq 或 CI 脚本断言
seqforge stats data.fasta --json | jq '.summary.n50'
seqforge validate data.fasta --json   # 退出码 2 即校验失败
```

### 📦 详细参数说明

以 `fqtrim` 为例：

| 参数 | 默认值 | 说明 |
| --- | --- | --- |
| `--quality` | 20 | 滑窗内的最低平均质量 |
| `--window` | 4 | 滑窗大小（碱基） |
| `--min-len` | 0 | 修剪后短于该长度的读段被丢弃 |
| `--max-len` | 无 | 长于该长度的读段被丢弃 |
| `--min-mean-q` | 0 | 平均质量低于该值的读段被丢弃 |
| `--max-n` | 1.0 | N 碱基占比上限 |

> 📌 **运行截图 / 演示动图占位**：终端输出为纯 ASCII 表格，可直接复制上方示例复现。

---

## 💡 设计思路与迭代规划

### 设计理念

1. **零依赖不是口号，是架构约束** —— 所有能力（含 Phred 偏移检测、N50、IUPAC 匹配、遗传密码表）均由标准库实现，代码全量可审计
2. **确定性优先** —— 库代码不访问网络、不依赖随机源（抽样必须显式传种子），保证科学结果可复现
3. **失败要响亮** —— 非法字符、截断的 FASTQ、长度不匹配的质量串都会立即报错并指出行号，绝不静默产出错误生物学结论
4. **CI 优先设计** —— `--json` 输出、稳定退出码、纯 ASCII 表格，让工具天然适合嵌入自动化流水线

### 技术选型原因

- **Python 标准库**：生物信息社区事实标准语言，`csv` / `gzip` / `statistics` / `argparse` 覆盖全部需求，零依赖分发成本最低
- **argparse 子命令**：无需引入 CLI 框架即可获得完整的 `--help` 生态与 Shell 补全潜力
- **生成器流式解析**：大文件逐记录产出，内存占用与记录数解耦

### 迭代规划

- **v1.1**：多行折叠 FASTQ 支持、FASTA 注释行（`;`）完整支持、更多 NCBI 密码表
- **v1.2**：SAM/BAM flag 解码、限制酶切割位点片段长度预测
- **v1.3**：交互式 TUI 仪表盘（纯标准库 curses）
- **v2.0**：插件机制，允许社区以零依赖方式扩展新命令

---

## 📦 打包与部署指南

### 本地打包（由 GitHub Actions 云端执行）

```bash
python -m pip install --upgrade build
python -m build        # 产出 dist/seqforge-1.0.0.tar.gz 与 .whl
```

> ☁️ 本项目**不在本地执行构建**：`build` 作业由云端 Runner 在 Ubuntu 上完成，产物自动上传为 Workflow Artifact。

### 安装部署

```bash
# 从源码目录安装（含 console script）
pip install .

# 或直接免安装使用
python -m seqforge <command>
```

### 兼容环境

| 平台 | 状态 |
| --- | --- |
| Windows 10/11 | ✅ 云端 CI 验证 |
| macOS 12+ | ✅ 云端 CI 验证 |
| Ubuntu 20.04+ | ✅ 云端 CI 验证 |
| Python 3.8 / 3.9 / 3.10 / 3.11 / 3.12 | ✅ 全版本 CI 矩阵 |

### 作为库集成

```python
from seqforge.seqio import read_fastq, write_fasta
from seqforge.quality import filter_fastq

records = read_fastq("raw.fastq")
clean = filter_fastq(records, min_len=50, min_mean_q=20)
write_fasta(clean, "clean.fasta")
```

---

## 🤝 贡献指南

欢迎各类贡献！本项目有两条铁律：**零第三方依赖** 与 **测试先行**。

```bash
# 本地验证（与 CI 完全一致）
python -m compileall -q seqforge tests
python -m unittest discover -s tests -v
```

- 🐛 发现 Bug？请提交 Issue 并附上最小复现输入
- 💡 有新想法？先开 Issue 讨论方案，避免零依赖红线被突破
- 🔧 提交代码请遵循 Angular 规范：`feat:` / `fix:` / `docs:` / `test:` / `refactor:` / `chore:`

详见 [CONTRIBUTING.md](CONTRIBUTING.md)。

---

## 📄 开源协议说明

本项目基于 [MIT License](LICENSE) 开源 —— 你可以自由地使用、修改、分发和商用，唯一要求是在副本中保留版权声明与许可声明。

```
Copyright (c) 2026 gitstq
```

---

<div align="center">

**🧬 SeqForge** — 让序列分析回归简单

[简体中文](README.md) | [繁體中文](README.zh-TW.md) | [English](README.en.md)

</div>
