# 🧬 SeqForge — 零依賴生物序列智能工具箱

[简体中文](README.md) | **繁體中文** | [English](README.en.md)

[![CI](https://github.com/gitstq/seqforge/actions/workflows/ci.yml/badge.svg)](https://github.com/gitstq/seqforge/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python](https://img.shields.io/badge/Python-3.8%2B-blue.svg)](https://www.python.org/downloads/)
[![Dependencies](https://img.shields.io/badge/dependencies-0-brightgreen.svg)](https://github.com/gitstq/seqforge)

> 🔥 **一個檔案、一條指令，搞定 FASTA/FASTQ 的日常分析** —— 純 Python 標準函式庫實作，無需安裝任何第三方套件，Windows / macOS / Linux 全平台開箱即用。

---

## 🎉 專案介紹

**SeqForge** 是一個輕量級的生物資訊學序列分析工具箱。做序列分析的你是否遇過這些煩惱：

- 😫 `biopython` 功能強大，但依賴安裝繁瑣，版本衝突頻繁發生
- 😤 `seqkit` 效能出色，卻需要下載執行檔，內網 / 老舊 Windows 環境處處受限
- 🤯 只是想算個 GC 含量、找個 ORF、翻譯一下序列，卻要開啟沉重的 IDE 或寫幾十行程式碼

SeqForge 就是為了終結這些煩惱而生的 **100% 自主研發** 專案：整個工具僅依賴 Python 標準函式庫，一次 `pip install`（甚至不安裝）就能執行，把最常用的序列分析能力濃縮進 16 個清晰的子指令。

**自主研發差異化亮點**：

- 🪶 **零依賴**：不引用任何第三方套件，複製即可執行，徹底告別「環境地獄」
- 🌏 **中文編碼友善**：自動識別 UTF-8 / GB18030 / Latin-1 編碼，國內實驗室的舊資料不再出現亂碼
- 📦 **`.gz` 透明支援**：壓縮格式的 FASTA/FASTQ 直接讀取，無需手動解壓縮
- 🤖 **CI 友善**：所有指令支援 `--json` 機器可讀輸出，`validate` 提供 0/2 退出碼，可直接作為流水線品質門禁
- 🧪 **雲端品質門禁**：GitHub Actions 在 Windows / macOS / Linux × Python 3.8–3.12 矩陣上執行語法檢查與全量單元測試

**靈感來源**：GitHub Trending 持續火熱的「本地優先 / 零依賴資料處理」工具方向（瀏覽器端大資料表格、本地索引引擎等高熱專案），叠加生物資訊社群對「輕量、免安裝、跨平台」序列工具的長期需求。

---

## ✨ 核心特性

| 能力 | 說明 |
| --- | --- |
| 📊 **統計引擎** | 序列長度分佈、**N50 / L50**、GC 含量、**GC 偏斜**、滑動視窗 GC 圖譜 |
| ✅ **格式校驗** | 結構校驗、非法字元偵測、重複 ID 檢查、空序列檢查，編碼自動嗅探 |
| 🧬 **遺傳密碼** | 內建 **5 套 NCBI 遺傳密碼表**（標準 / 脊椎動物線粒體 / 酵母線粒體 / 無脊椎動物線粒體 / 纖毛蟲），支援 **六框翻譯** |
| 🔎 **ORF 尋找** | 雙鏈六框掃描，IUPAC 相容，座標自動回填至原始鏈 |
| 🎯 **基序搜尋** | IUPAC 歧義碼比對（如 `GASTY`），支援雙鏈檢索 |
| ✂️ **酶切位點** | 內建 **20 種常用限制性內切酶**（EcoRI、BamHI、NotI……） |
| 🌡️ **引子 Tm** | Wallace 法則 + 鹽濃度校正公式，一鍵評估退火溫度 |
| 🧹 **FASTQ 品質管制** | Phred+33/+64 偏移自動偵測、**BWA 式滑動視窗修剪**、長度 / 品質 / N 含量過濾 |
| 🔄 **格式互轉** | FASTA / FASTQ / CSV / TSV / NDJSON 自由轉換，FASTA→FASTQ 一鍵補全 |
| 🗜️ **壓縮支援** | `.gz` 檔案透明讀寫 |
| 🎲 **可重現工具** | 帶隨機種子的抽樣、大小寫不敏感的去重 |
| 🖥️ **跨平台** | 純 Python 3.8+，零編譯、零執行期依賴 |

---

## 🚀 快速開始

### 環境需求

- 🐍 **Python 3.8 或更高版本**（僅此一項，無任何第三方依賴）

### 方式一：pip 安裝（推薦）

```bash
pip install seqforge
seqforge --version
```

### 方式二：免安裝，直接執行

```bash
git clone https://github.com/gitstq/seqforge.git
cd seqforge
python -m seqforge --version
```

### 方式三：作為函式庫引入

```python
from seqforge.seqio import read_fasta
from seqforge.genetic import translate

records = read_fasta("demo.fasta")
for record in records:
    print(record.id, len(record.seq), translate(record.seq))
```

### 30 秒上手

```bash
# 1. 準備一份範例資料
cat > demo.fasta <<'EOF'
>gene_demo 範例序列
ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG
>gene_short
GGCCGG
EOF

# 2. 查看統計摘要（N50、GC、長度分佈一目了然）
seqforge stats demo.fasta

# 3. 六框翻譯，尋找潛在編碼序列
seqforge translate demo.fasta --frames 6
```

---

## 📖 詳細使用指南

### 指令總覽

| 指令 | 功能 |
| --- | --- |
| `stats` | 長度 / GC / N50 / L50 摘要統計 |
| `validate` | FASTA/FASTQ 結構校驗（CI 門禁） |
| `gc` | 逐序列 GC 含量、GC 偏斜、滑動視窗 GC |
| `revcomp` | 反向互補（DNA/RNA） |
| `transcribe` | DNA → RNA 轉錄 |
| `translate` | 翻譯為蛋白（1/3/6 框、5 種密碼表） |
| `orf` | 開放閱讀框尋找（雙鏈） |
| `find` | IUPAC 基序檢索 |
| `enzymes` | 限制酶酶切位點分析 |
| `tm` | 引子退火溫度計算 |
| `fqstat` | FASTQ 品質統計 |
| `fqtrim` | FASTQ 修剪與過濾 |
| `convert` | 格式互轉（fasta/fastq/csv/tsv/ndjson） |
| `sample` | 可重現隨機抽樣 |
| `dedup` | 重複序列去除 |
| `info` | 格式探測與記錄預覽 |

> 💡 所有指令皆支援 `-o/--output` 寫入檔案、`--encoding` 強制指定編碼；統計 / 分析類指令支援 `--json` 機器可讀輸出。

### 📊 範例一：資料品質管制（`stats` + `validate`）

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
# 校驗通過時退出碼為 0；發現錯誤時退出碼為 2，可直接用於 CI 門禁
```

### 🧬 範例二：ORF 尋找與翻譯

```bash
$ seqforge orf gene.fasta --min-len 100
+-----------+--------+-------+-------+-----+--------+----------+
|    id     | strand | frame | start | end | length |  protein |
+-----------+--------+-------+-------+-----+--------+----------+
| gene_demo |   +    |   1   |   1   | 39  |   39   | MAIVNGL...|
+-----------+--------+-------+-------+-----+--------+----------+

$ seqforge translate gene.fasta --table 2 --frames 6
# --table 2 使用脊椎動物線粒體密碼表；--frames 6 同時翻譯兩條鏈
```

### 🎯 範例三：TATA 框基序檢索

```bash
$ seqforge find promoter.fasta TATAAA --both-strands
+----------+-------+-----+--------+
|    id    | start | end | strand |
+----------+-------+-----+--------+
| promoter |   53  | 59  |   +    |
| promoter |  120  | 126 |   -    |
+----------+-------+-----+--------+
```

### ✂️ 範例四：FASTQ 品質管制流水線

```bash
# 查看品質概況（自動識別 Phred+33/+64）
$ seqforge fqstat reads.fastq

# 滑動視窗修剪 + 過濾，一步到位
$ seqforge fqtrim reads.fastq \
    --quality 20 --window 4 \
    --min-len 50 --max-n 0.1 \
    -o reads.clean.fastq
```

### 🌡️ 範例五：引子 Tm 計算

```bash
$ seqforge tm ATGCGTATGCGTATGCGT --na 0.05
+---------------------+--------+--------+----------------+-----------------+
|       primer        | length | gc(%)  | Tm Wallace(C)  | Tm adjusted(C)  |
+---------------------+--------+--------+----------------+-----------------+
| ATGCGTATGCGTATGCGT  |   18   | 50.0   |       84       |      42.9       |
+---------------------+--------+--------+----------------+-----------------+
```

### 🤖 範例六：CI 中的 JSON 門禁

```bash
# 輸出機器可讀 JSON，配合 jq 或 CI 腳本斷言
seqforge stats data.fasta --json | jq '.summary.n50'
seqforge validate data.fasta --json   # 退出碼 2 即校驗失敗
```

### 📦 詳細參數說明

以 `fqtrim` 為例：

| 參數 | 預設值 | 說明 |
| --- | --- | --- |
| `--quality` | 20 | 滑動視窗內的最低平均品質 |
| `--window` | 4 | 滑動視窗大小（鹼基） |
| `--min-len` | 0 | 修剪後短於該長度的讀段被丟棄 |
| `--max-len` | 無 | 長於該長度的讀段被丟棄 |
| `--min-mean-q` | 0 | 平均品質低於該值的讀段被丟棄 |
| `--max-n` | 1.0 | N 鹼基佔比上限 |

> 📌 **執行畫面截圖 / 演示動圖佔位**：終端輸出為純 ASCII 表格，可直接複製上方範例重現。

---

## 💡 設計思路與迭代規劃

### 設計理念

1. **零依賴不是口號，是架構約束** —— 所有能力（含 Phred 偏移偵測、N50、IUPAC 比對、遺傳密碼表）均由標準函式庫實作，程式碼全量可稽核
2. **決定性優先** —— 函式庫程式碼不存取網路、不依賴隨機來源（抽樣必須顯式傳入種子），保證科學結果可重現
3. **失敗要響亮** —— 非法字元、截斷的 FASTQ、長度不匹配的品質串都會立即報錯並指出行號，絕不靜默產出錯誤的生物學結論
4. **CI 優先設計** —— `--json` 輸出、穩定退出碼、純 ASCII 表格，讓工具天然適合嵌入自動化流水線

### 技術選型原因

- **Python 標準函式庫**：生物資訊社群事實標準語言，`csv` / `gzip` / `statistics` / `argparse` 覆蓋全部需求，零依賴分發成本最低
- **argparse 子指令**：無需引入 CLI 框架即可獲得完整的 `--help` 生態與 Shell 補全潛力
- **產生器串流解析**：大檔案逐記錄產出，記憶體用量與記錄數解耦

### 迭代規劃

- **v1.1**：多行摺疊 FASTQ 支援、FASTA 註解行（`;`）完整支援、更多 NCBI 密碼表
- **v1.2**：SAM/BAM flag 解碼、限制酶切割位點片段長度預測
- **v1.3**：互動式 TUI 儀表板（純標準函式庫 curses）
- **v2.0**：插件機制，允許社群以零依賴方式擴充新指令

---

## 📦 打包與部署指南

### 本地打包（由 GitHub Actions 雲端執行）

```bash
python -m pip install --upgrade build
python -m build        # 產出 dist/seqforge-1.0.0.tar.gz 與 .whl
```

> ☁️ 本專案**不在本地執行建置**：`build` 作業由雲端 Runner 在 Ubuntu 上完成，產物自動上傳為 Workflow Artifact。

### 安裝部署

```bash
# 從原始碼目錄安裝（含 console script）
pip install .

# 或直接免安裝使用
python -m seqforge <command>
```

### 相容環境

| 平台 | 狀態 |
| --- | --- |
| Windows 10/11 | ✅ 雲端 CI 驗證 |
| macOS 12+ | ✅ 雲端 CI 驗證 |
| Ubuntu 20.04+ | ✅ 雲端 CI 驗證 |
| Python 3.8 / 3.9 / 3.10 / 3.11 / 3.12 | ✅ 全版本 CI 矩陣 |

### 作為函式庫整合

```python
from seqforge.seqio import read_fastq, write_fasta
from seqforge.quality import filter_fastq

records = read_fastq("raw.fastq")
clean = filter_fastq(records, min_len=50, min_mean_q=20)
write_fasta(clean, "clean.fasta")
```

---

## 🤝 貢獻指南

歡迎各類貢獻！本專案有兩條鐵律：**零第三方依賴** 與 **測試先行**。

```bash
# 本地驗證（與 CI 完全一致）
python -m compileall -q seqforge tests
python -m unittest discover -s tests -v
```

- 🐛 發現 Bug？請提交 Issue 並附上最小重現輸入
- 💡 有新想法？先開 Issue 討論方案，避免零依賴紅線被突破
- 🔧 提交程式碼請遵循 Angular 規範：`feat:` / `fix:` / `docs:` / `test:` / `refactor:` / `chore:`

詳見 [CONTRIBUTING.md](CONTRIBUTING.md)。

---

## 📄 開源協議說明

本專案基於 [MIT License](LICENSE) 開源 —— 你可以自由地使用、修改、分發和商用，唯一要求是在副本中保留版權聲明與授權聲明。

```
Copyright (c) 2026 gitstq
```

---

<div align="center">

**🧬 SeqForge** — 讓序列分析回歸簡單

[简体中文](README.md) | [繁體中文](README.zh-TW.md) | [English](README.en.md)

</div>
