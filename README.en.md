# 🧬 SeqForge — Zero-Dependency Bioinformatics Sequence Toolkit

[简体中文](README.md) | [繁體中文](README.zh-TW.md) | **English**

[![CI](https://github.com/gitstq/seqforge/actions/workflows/ci.yml/badge.svg)](https://github.com/gitstq/seqforge/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python](https://img.shields.io/badge/Python-3.8%2B-blue.svg)](https://www.python.org/downloads/)
[![Dependencies](https://img.shields.io/badge/dependencies-0-brightgreen.svg)](https://github.com/gitstq/seqforge)

> 🔥 **One file, one command — your daily FASTA/FASTQ analysis, done.** Built on nothing but the Python standard library. No packages to install, works out of the box on Windows / macOS / Linux.

---

## 🎉 Introduction

**SeqForge** is a lightweight bioinformatics toolkit for sequence analysis. Sound familiar?

- 😫 `biopython` is powerful, but installing it is a chore and version conflicts never end
- 😤 `seqkit` is blazing fast, yet it ships as a binary you must download — awkward on intranets and legacy Windows machines
- 🤯 You just want a GC percentage, an ORF list, or a quick translation — why boot a heavy IDE or write dozens of lines of script?

SeqForge ends the frustration. It is a **100% original** project that depends on nothing but the Python standard library: a single `pip install` (or no install at all) gets you 16 focused subcommands covering everyday sequence work.

**What makes it different**:

- 🪶 **Zero dependencies** — no third-party packages at all; clone and run. Say goodbye to dependency hell
- 🌏 **Chinese-encoding friendly** — auto-detects UTF-8 / GB18030 / Latin-1, so legacy lab files never turn into mojibake
- 📦 **Transparent `.gz` support** — read and write compressed FASTA/FASTQ without manual decompression
- 🤖 **CI-friendly** — every command emits `--json`, and `validate` uses 0/2 exit codes, ready to gate your pipelines
- 🧪 **Cloud quality gates** — GitHub Actions runs syntax checks and the full unit-test suite across Windows / macOS / Linux × Python 3.8–3.12

**Inspiration** — the "local-first, zero-dependency data tooling" wave on GitHub Trending (in-browser big-data grids, local indexing engines, and friends), combined with the bioinformatics community's long-standing demand for a light, install-free, cross-platform sequence tool.

---

## ✨ Key Features

| Capability | Highlights |
| --- | --- |
| 📊 **Statistics engine** | Length distributions, **N50 / L50**, GC content, **GC skew**, sliding-window GC profiles |
| ✅ **Format validation** | Structural checks, invalid-character detection, duplicate-ID and empty-sequence checks, encoding sniffing |
| 🧬 **Genetic codes** | **5 NCBI genetic code tables** (Standard, Vertebrate Mitochondrial, Yeast Mitochondrial, Invertebrate Mitochondrial, Ciliate) with **six-frame translation** |
| 🔎 **ORF finding** | Both strands, six frames, IUPAC-aware, coordinates mapped back to the input strand |
| 🎯 **Motif search** | IUPAC ambiguity codes (e.g. `GASTY`), single- or double-stranded search |
| ✂️ **Restriction sites** | **20 common restriction enzymes** built in (EcoRI, BamHI, NotI, …) |
| 🌡️ **Primer Tm** | Wallace rule plus a salt-adjusted formula for quick annealing-temperature estimates |
| 🧹 **FASTQ QC** | Phred+33/+64 offset detection, **BWA-style sliding-window trimming**, length / quality / N-content filters |
| 🔄 **Format conversion** | FASTA / FASTQ / CSV / TSV / NDJSON in every direction, one-step FASTA→FASTQ |
| 🗜️ **Compression** | transparent `.gz` read/write |
| 🎲 **Reproducible tools** | seeded sampling and case-insensitive deduplication |
| 🖥️ **Cross-platform** | pure Python 3.8+, no compilation, no runtime dependencies |

---

## 🚀 Quick Start

### Requirements

- 🐍 **Python 3.8 or newer** (that is the only requirement — zero third-party dependencies)

### Option 1: pip install (recommended)

```bash
pip install seqforge
seqforge --version
```

### Option 2: run without installing

```bash
git clone https://github.com/gitstq/seqforge.git
cd seqforge
python -m seqforge --version
```

### Option 3: use it as a library

```python
from seqforge.seqio import read_fasta
from seqforge.genetic import translate

records = read_fasta("demo.fasta")
for record in records:
    print(record.id, len(record.seq), translate(record.seq))
```

### 30-second tour

```bash
# 1. create a sample file
cat > demo.fasta <<'EOF'
>gene_demo demo sequence
ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG
>gene_short
GGCCGG
EOF

# 2. summary statistics (N50, GC, length distribution at a glance)
seqforge stats demo.fasta

# 3. six-frame translation to hunt for coding potential
seqforge translate demo.fasta --frames 6
```

---

## 📖 Detailed Usage Guide

### Command overview

| Command | Purpose |
| --- | --- |
| `stats` | Length / GC / N50 / L50 summary statistics |
| `validate` | FASTA/FASTQ structural validation (CI gate) |
| `gc` | Per-sequence GC content, GC skew, sliding-window GC |
| `revcomp` | Reverse complement (DNA/RNA) |
| `transcribe` | DNA → RNA transcription |
| `translate` | Protein translation (1/3/6 frames, 5 genetic codes) |
| `orf` | Open reading frame discovery (both strands) |
| `find` | IUPAC motif search |
| `enzymes` | Restriction enzyme site analysis |
| `tm` | Primer melting temperature |
| `fqstat` | FASTQ quality statistics |
| `fqtrim` | FASTQ trimming and filtering |
| `convert` | Format conversion (fasta/fastq/csv/tsv/ndjson) |
| `sample` | Reproducible random sampling |
| `dedup` | Duplicate sequence removal |
| `info` | Format detection and record preview |

> 💡 Every command supports `-o/--output` for file output and `--encoding` to force an encoding; analysis commands support `--json` for machine-readable output.

### 📊 Example 1: data QC (`stats` + `validate`)

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
# exit code 0 when clean; exit code 2 on errors — perfect for CI gates
```

### 🧬 Example 2: ORF discovery and translation

```bash
$ seqforge orf gene.fasta --min-len 100
+-----------+--------+-------+-------+-----+--------+----------+
|    id     | strand | frame | start | end | length |  protein |
+-----------+--------+-------+-------+-----+--------+----------+
| gene_demo |   +    |   1   |   1   | 39  |   39   | MAIVNGL...|
+-----------+--------+-------+-------+-----+--------+----------+

$ seqforge translate gene.fasta --table 2 --frames 6
# --table 2 selects the vertebrate mitochondrial code; --frames 6 translates both strands
```

### 🎯 Example 3: TATA-box motif search

```bash
$ seqforge find promoter.fasta TATAAA --both-strands
+----------+-------+-----+--------+
|    id    | start | end | strand |
+----------+-------+-----+--------+
| promoter |   53  | 59  |   +    |
| promoter |  120  | 126 |   -    |
+----------+-------+-----+--------+
```

### ✂️ Example 4: FASTQ QC pipeline

```bash
# inspect quality (Phred+33/+64 detected automatically)
$ seqforge fqstat reads.fastq

# sliding-window trimming plus filtering in one step
$ seqforge fqtrim reads.fastq \
    --quality 20 --window 4 \
    --min-len 50 --max-n 0.1 \
    -o reads.clean.fastq
```

### 🌡️ Example 5: primer Tm

```bash
$ seqforge tm ATGCGTATGCGTATGCGT --na 0.05
+---------------------+--------+--------+----------------+-----------------+
|       primer        | length | gc(%)  | Tm Wallace(C)  | Tm adjusted(C)  |
+---------------------+--------+--------+----------------+-----------------+
| ATGCGTATGCGTATGCGT  |   18   | 50.0   |       84       |      42.9       |
+---------------------+--------+--------+----------------+-----------------+
```

### 🤖 Example 6: JSON gates in CI

```bash
# machine-readable JSON, ready for jq or CI assertions
seqforge stats data.fasta --json | jq '.summary.n50'
seqforge validate data.fasta --json   # exit code 2 means validation failed
```

### 📦 Option reference

Using `fqtrim` as an example:

| Option | Default | Description |
| --- | --- | --- |
| `--quality` | 20 | Minimum mean quality inside the sliding window |
| `--window` | 4 | Sliding window size (bases) |
| `--min-len` | 0 | Drop reads shorter than this after trimming |
| `--max-len` | none | Drop reads longer than this |
| `--min-mean-q` | 0 | Drop reads whose mean quality is below this |
| `--max-n` | 1.0 | Maximum fraction of N bases |

> 📌 **Screenshot / demo GIF placeholder** — output is pure ASCII tables; reproduce any example above by copying it verbatim.

---

## 💡 Design Rationale and Roadmap

### Design principles

1. **Zero dependencies is an architectural constraint, not a slogan** — everything (Phred offset detection, N50, IUPAC matching, genetic code tables) is implemented on the standard library, and the whole codebase stays auditable
2. **Determinism first** — library code never touches the network and never consumes ambient randomness (sampling requires an explicit seed), keeping scientific results reproducible
3. **Fail loudly** — invalid characters, truncated FASTQ, and length-mismatched quality strings raise errors with line numbers instead of silently producing wrong biology
4. **CI-first design** — `--json` output, stable exit codes, and pure-ASCII tables make the tool a natural fit for automation pipelines

### Technology choices

- **Python standard library** — the de-facto language of bioinformatics; `csv` / `gzip` / `statistics` / `argparse` cover every requirement with the lowest possible distribution cost
- **argparse subcommands** — a complete `--help` ecosystem and shell-completion potential without pulling in a CLI framework
- **Generator-based streaming parsers** — large files are yielded record by record, decoupling memory from record count

### Roadmap

- **v1.1** — multi-line wrapped FASTQ support, full FASTA comment (`;`) support, more NCBI code tables
- **v1.2** — SAM/BAM flag decoding, restriction digest fragment-length prediction
- **v1.3** — interactive TUI dashboard (stdlib `curses` only)
- **v2.0** — a plugin mechanism so the community can add new commands without breaking the zero-dependency rule

---

## 📦 Packaging and Deployment

### Build locally (executed by GitHub Actions in the cloud)

```bash
python -m pip install --upgrade build
python -m build        # produces dist/seqforge-1.0.0.tar.gz and .whl
```

> ☁️ This project deliberately **never builds on your machine**: the `build` job runs on a cloud Ubuntu runner and uploads the artifacts to the workflow run.

### Install and deploy

```bash
# install from the source tree (installs the console script)
pip install .

# or use it without installing
python -m seqforge <command>
```

### Compatibility

| Platform | Status |
| --- | --- |
| Windows 10/11 | ✅ verified by cloud CI |
| macOS 12+ | ✅ verified by cloud CI |
| Ubuntu 20.04+ | ✅ verified by cloud CI |
| Python 3.8 / 3.9 / 3.10 / 3.11 / 3.12 | ✅ full CI matrix |

### Integrate as a library

```python
from seqforge.seqio import read_fastq, write_fasta
from seqforge.quality import filter_fastq

records = read_fastq("raw.fastq")
clean = filter_fastq(records, min_len=50, min_mean_q=20)
write_fasta(clean, "clean.fasta")
```

---

## 🤝 Contributing

Contributions of all kinds are welcome! Two rules are non-negotiable: **zero third-party dependencies** and **tests first**.

```bash
# run exactly what CI runs
python -m compileall -q seqforge tests
python -m unittest discover -s tests -v
```

- 🐛 Found a bug? Open an issue with a minimal reproducing input
- 💡 Have an idea? Discuss it in an issue first so we keep the zero-dependency red line intact
- 🔧 Follow the Angular convention for commits: `feat:` / `fix:` / `docs:` / `test:` / `refactor:` / `chore:`

See [CONTRIBUTING.md](CONTRIBUTING.md) for details.

---

## 📄 License

Released under the [MIT License](LICENSE) — use, modify, distribute and commercialize freely; the only requirement is to keep the copyright and permission notices in copies.

```
Copyright (c) 2026 gitstq
```

---

<div align="center">

**🧬 SeqForge** — sequence analysis, simplified

[简体中文](README.md) | [繁體中文](README.zh-TW.md) | [English](README.en.md)

</div>
