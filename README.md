# WGS SNP Calling Pipeline

基于 Snakemake 的全基因组重测序 SNP/INDEL 检测流程，覆盖从原始 FASTQ 到过滤后 VCF 的完整分析步骤，并附有质控统计工具。

---

## 目录结构

```
.
├── Snakefile                        # 主流程（Snakemake）
├── config/
│   ├── config.yaml                  # 主配置文件（路径、软件、参数）
│   └── chroms.txt                   # 染色体列表（1-22, X, Y, MT）
├── bin/
│   ├── ngs_qc_summary.py            # BAM/fastq QC 汇总 → Excel
│   ├── calc_qc_stats.sh             # VCF 样本质控统计（TiTv、Het/Hom、性别预测）
│   └── calc_qc_stats_readme.md      # calc_qc_stats.sh 详细说明
└── README.md
```

---

## 流程概览

```
fastp（质控过滤）
  └── bwa mem + samtools sort（比对与排序）
        └── sambamba markdup（标记 PCR 重复）
              ├── samtools idxstats / flagstat / stats（BAM 统计）
              ├── mosdepth（覆盖度统计）
              └── GATK HaplotypeCaller（-ERC GVCF，按染色体并行）
                    └── GATK CombineGVCFs（单样本各染色体 GVCF 合并）
                          └── GATK GenomicsDBImport（按染色体，多样本联合导入）
                                └── GATK GenotypeGVCFs（按染色体联合基因分型）
                                      ├── GATK SelectVariants + VariantFiltration（SNP Hard Filter）
                                      └── GATK SelectVariants + VariantFiltration（INDEL Hard Filter）
                                            └── GATK MergeVcfs（全基因组合并）→ tabix 建索引
                                                  └── MultiQC（汇总 QC 报告）
```

**辅助工具（流程结束后手动运行）**

- `bin/ngs_qc_summary.py`：汇总各样本 fastp / BAM / mosdepth QC 指标 → Excel
- `bin/calc_qc_stats.sh`：对最终 VCF 逐样本计算 TiTv、Het/Hom、性别预测等质控指标

---

## 软件版本

| 软件 | 版本 | 用途 |
|------|------|------|
| Snakemake | 9.16.3 | 工作流引擎 |
| fastp | 0.20.1 | FASTQ 质控与接头过滤 |
| BWA | 0.7.18 | 短读段比对 |
| samtools | 1.21 (htslib 1.21) | BAM 处理与统计 |
| sambamba | 1.0.1 | 标记 PCR 重复、BAM 索引 |
| mosdepth | 0.3.10 | 测序深度/覆盖度统计 |
| GATK | 4.1.8.1 | 变异检测（HaplotypeCaller / GenomicsDB / GenotypeGVCFs / 过滤） |
| tabix | 1.21 (htslib) | VCF/GVCF 索引 |
| MultiQC | 1.27.1 | 多样本 QC 汇总报告 |
| bcftools | 1.17 | VCF 统计（calc_qc_stats.sh 使用） |
| Python | ≥ 3.8 | ngs_qc_summary.py（需 pandas、openpyxl） |

---

## 各步骤详细参数说明

### 1. fastp — FASTQ 质控

**线程：** `threads.fastp`（默认 8）  
**内存：** 8 GB  
**主要参数：**

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--detect_adapter_for_pe` | — | 自动检测 PE 接头序列 |
| `--correction` | — | 开启碱基纠错（overlap 校正） |
| `--overrepresentation_analysis` | — | 过表达序列分析 |
| `--length_required` | 50 | 过滤短于 50 bp 的 reads |
| `--qualified_quality_phred` | 20 | 碱基质量低于 Q20 视为低质量 |
| `--unqualified_percent_limit` | 40 | 低质量碱基超过 40% 则过滤该 read |
| `--n_base_limit` | 5 | N 碱基超过 5 个则过滤 |
| `--thread` | 8 | 线程数 |

**输出：** 过滤后 FASTQ、HTML/JSON 报告（供 MultiQC 使用）

---

### 2. BWA MEM — 短读段比对

**线程：** `threads.bwa`（默认 16）  
**内存：** 32 GB  
**主要参数：**

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `-t` | 16 | 比对线程数 |
| `-K 100000000` | — | 每批次处理的碱基数（提高确定性，便于并行复现） |
| `-Y` | — | 次级比对使用软裁切（推荐用于 GATK 流程） |
| `-R` | `@RG\tID:{sample}\tSM:{sample}\tPL:ILLUMINA\tLB:{sample}` | Read Group 信息（GATK 必须） |

**管道接 samtools sort：**

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `-@` | `threads.sort`（默认 4）| 排序附加线程数 |

**输出：** 排序后 BAM（临时文件，流程完成后自动清理）

---

### 3. sambamba markdup — 标记 PCR 重复

**线程：** `threads.markdup`（默认 8）  
**内存：** 32 GB  
**主要参数：**

| 参数 | 值 | 说明 |
|------|-----|------|
| `-t` | 8 | 线程数 |
| `--hash-table-size` | 1000000 | 哈希表大小（大样本建议调大） |
| `--overflow-list-size` | 1000000 | 溢出链表大小 |

**输出：** markdup BAM + BAM 索引（.bai）、重复统计文件（.markdup_metrics.txt）

---

### 4. samtools — BAM 质控统计

**线程：** 4  
**内存：** 8 GB  
**运行三个子命令：**

| 子命令 | 输出 | 说明 |
|--------|------|------|
| `idxstats` | `{sample}.idxstats.txt` | 各染色体比对 reads 数 |
| `flagstat` | `{sample}.flagstat.txt` | 比对率、配对率、重复率等汇总 |
| `stats --reference {ref}` | `{sample}.stats.txt` | 详细统计，含插入片段长度分布、错误率等 |

---

### 5. mosdepth — 测序深度/覆盖度

**线程：** `threads.mosdepth`（默认 4）  
**内存：** 16 GB  
**主要参数：**

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--threads` | 4 | 线程数 |
| `--fast-mode` | — | 跳过 per-base 深度文件，仅统计 window 和分布（速度更快） |
| `--by` | 500 | 统计窗口大小（bp） |
| `--thresholds` | `1,5,10,15,20,25,30` | 输出各深度阈值下的覆盖率 |

**输出：** `{sample}.mosdepth.summary.txt`（全基因组平均深度）、`{sample}.mosdepth.global.dist.txt`（累积覆盖率分布）、`{sample}.regions.bed.gz`（各窗口深度）

---

### 6. GATK HaplotypeCaller — 单样本变异检测（GVCF 模式）

**线程：** 2（`--native-pair-hmm-threads`）  
**内存：** 16 GB  
**JVM 参数：** `-Xmx14g -XX:+UseParallelGC -XX:ParallelGCThreads=4`  
**主要参数：**

| 参数 | 值 | 说明 |
|------|-----|------|
| `-R` | 参考基因组 | 参考序列 |
| `-I` | markdup BAM | 输入 BAM |
| `-L` | 染色体名 | 区间限定（按染色体并行） |
| `-ERC GVCF` | — | 输出 GVCF 格式（含参考基因型，供联合基因分型使用） |
| `--native-pair-hmm-threads` | 2 | PairHMM 加速线程数 |

**并行策略：** 每个样本 × 每条染色体独立运行，充分利用集群资源；per-chrom GVCF 标记为临时文件，合并后自动清理。

---

### 7. GATK CombineGVCFs — 单样本 GVCF 合并

**线程：** 4  
**内存：** 32 GB  
**JVM 参数：** `-Xmx28g -XX:+UseParallelGC -XX:ParallelGCThreads=4`  
**说明：** 将同一样本各染色体的 GVCF 合并为全基因组单文件，供 GenomicsDBImport 使用。

---

### 8. GATK GenomicsDBImport — 多样本 GVCF 联合导入

**线程：** 4（`--reader-threads`）  
**内存：** 32 GB  
**JVM 参数：** `-Xmx28g -XX:+UseParallelGC -XX:ParallelGCThreads=4`  
**主要参数：**

| 参数 | 值 | 说明 |
|------|-----|------|
| `-V` | 各样本 GVCF | 输入 GVCF（每个样本一个 `-V`） |
| `--genomicsdb-workspace-path` | `genomicsdb/{chrom}` | GenomicsDB 输出路径 |
| `-L` | 染色体名 | 区间限定 |
| `--reader-threads` | 4 | 并行读取线程 |
| `--tmp-dir` | `results/tmp` | 临时文件目录 |

---

### 9. GATK GenotypeGVCFs — 联合基因分型

**线程：** 4  
**内存：** 32 GB  
**JVM 参数：** `-Xmx28g -XX:+UseParallelGC -XX:ParallelGCThreads=4`  
**主要参数：**

| 参数 | 值 | 说明 |
|------|-----|------|
| `-R` | 参考基因组 | 参考序列 |
| `-V` | `gendb://{chrom}` | GenomicsDB 输入 |
| `-L` | 染色体名 | 区间限定 |
| `-O` | 输出 VCF | 联合基因分型结果 |

---

### 10. GATK Hard Filter — SNP/INDEL 硬过滤

**线程：** 2  
**内存：** 16 GB  
**JVM 参数：** `-Xmx14g -XX:+UseParallelGC -XX:ParallelGCThreads=4`  

每个染色体的过滤分三步：
1. `SelectVariants` 提取 SNP 或 INDEL
2. `VariantFiltration` 打标记（不符合条件的写入 FILTER 列）
3. `SelectVariants --exclude-filtered` 仅保留 PASS 位点

**SNP 过滤条件（`snp_filter_expr`）：**

| 条件 | 阈值 | 说明 |
|------|------|------|
| `QD < 2.0` | — | 质量深度比过低 |
| `FS > 60.0` | — | 链偏差（Fisher Strand）过高 |
| `MQ < 40.0` | — | 比对质量过低 |
| `MQRankSum < -12.5` | — | Mapping Quality Rank Sum Test |
| `ReadPosRankSum < -8.0` | — | 位点在 read 末端聚集 |

**INDEL 过滤条件（`indel_filter_expr`）：**

| 条件 | 阈值 | 说明 |
|------|------|------|
| `QD < 2.0` | — | 质量深度比过低 |
| `FS > 200.0` | — | 链偏差过高（INDEL 阈值更宽松） |
| `ReadPosRankSum < -20.0` | — | 位点在 read 末端聚集 |

---

### 11. GATK MergeVcfs — 全基因组 VCF 合并

**线程：** 4  
**内存：** 32 GB  
**JVM 参数：** `-Xmx28g -XX:+UseParallelGC -XX:ParallelGCThreads=4`  
**说明：** 将所有染色体的过滤后 SNP VCF 和 INDEL VCF 合并为全基因组单文件，并由 tabix 建立 .tbi 索引。

---

### 12. MultiQC — QC 报告汇总

**线程：** 1  
**内存：** 8 GB  
**主要参数：**

| 参数 | 值 | 说明 |
|------|-----|------|
| `--outdir` | `results/qc/multiqc` | 输出目录 |
| `--filename` | `multiqc_report` | 输出文件名 |
| `--interactive` | — | 生成交互式图表 |
| `--force` | — | 覆盖已有结果 |

**汇总内容：** fastp 报告、samtools flagstat/idxstats/stats、sambamba markdup 统计、mosdepth 深度统计。

---

## 辅助工具

### ngs_qc_summary.py — BAM QC 汇总

从 `results/qc/` 目录自动提取各样本质控指标，输出带条件格式的 Excel 文件（含 QC汇总、fastp质控、比对与深度三个 Sheet）。

**依赖：** Python ≥ 3.8，pandas，openpyxl

**用法：**
```bash
python3 bin/ngs_qc_summary.py -i results/ -o qc_summary.xlsx

# 自定义深度覆盖率阈值
python3 bin/ngs_qc_summary.py -i results/ -o qc_summary.xlsx --depth-thresholds 10 20 30 50
```

**参数：**

| 参数 | 必填 | 默认值 | 说明 |
|------|------|--------|------|
| `-i` | 是 | — | results 目录（含 qc/ 子目录） |
| `-o` | 否 | `qc_summary.xlsx` | 输出 Excel 路径 |
| `--depth-thresholds` | 否 | `1 5 10 20 30 50` | mosdepth 覆盖率统计阈值（X） |
| `--qc-subdir` | 否 | `qc` | QC 子目录名 |

---

### calc_qc_stats.sh — VCF 样本质控统计

使用 bcftools v1.17 对最终 VCF 逐样本统计 SNP/INDEL 数量、Ti/Tv 比值、Het/Hom 比值，并根据 X 染色体杂合率预测样本性别。

**用法：**
```bash
bash bin/calc_qc_stats.sh -i results/vcf/final/all_samples.filtered.vcf.gz

# 自定义输出文件名
bash bin/calc_qc_stats.sh -i all_samples.filtered.vcf.gz -o cohort_qc.txt

# 指定 bcftools 路径
bash bin/calc_qc_stats.sh -i all_samples.filtered.vcf.gz -b /path/to/bcftools
```

**参数：**

| 参数 | 必填 | 默认值 | 说明 |
|------|------|--------|------|
| `-i` | 是 | — | 输入 VCF（需 bgzip + .tbi 索引） |
| `-o` | 否 | `sample_vcfqc_stats.txt` | 输出文件路径 |
| `-b` | 否 | PATH 中的 bcftools | 指定 bcftools 路径 |

**输出字段：** Sample、SNPs、InDels、Ts、Tv、TiTv_ratio、Hets_SNP、HomAlt_SNP、Het_Hom_ratio、X_Het、X_HomAlt、X_Het_rate、Predicted_Sex

**Ti/Tv 参考范围：** WGS 1.9～2.1，WES 2.8～3.0  
**性别判断阈值：** X 杂合率 ≥ 0.25 → Female；≤ 0.15 → Male；其余 → Unknown

---

## 使用方法

### 1. 修改配置文件

编辑 `config/config.yaml`，至少需要设置以下路径：

```yaml
sample_dir:  /path/to/your/samples    # 样本根目录
reference:   /path/to/genome.fa       # 参考基因组（需有 .fai 索引）
chrom_list:  config/chroms.txt        # 染色体列表
outdir:      results                  # 输出根目录
```

样本目录格式（自动识别，无需样本表）：
```
sample_dir/
  sample1/
    xxx_R1.fq.gz  xxx_R2.fq.gz
  sample2/
    xxx_1.fq.gz   xxx_2.fq.gz
```

### 2. 修改染色体列表

编辑 `config/chroms.txt`，每行一个染色体名，需与参考基因组 FASTA header 完全一致（默认为 GRCh38 无 `chr` 前缀，1-22, X, Y, MT）。

### 3. 运行流程

```bash
# 试运行（检查任务依赖，不实际执行）
snakemake -s Snakefile --configfile config/config.yaml -n

# 本地多核运行
snakemake -s Snakefile --configfile config/config.yaml --cores 32

# SGE/PBS 集群投递（示例）
snakemake -s Snakefile --configfile config/config.yaml \
    --cluster "qsub -l h_vmem={resources.mem_gb}G -pe smp {threads}" \
    --jobs 50
```

---

## 输出结构

```
results/
├── clean_fastq/{sample}/                        # fastp 过滤后 FASTQ
├── bam/
│   ├── sorted/{sample}.sorted.bam               # 排序 BAM（临时，自动清理）
│   └── markdup/{sample}.markdup.bam(.bai)       # 标记重复 BAM + 索引
├── gvcf/{sample}/{chrom}.g.vcf.gz(.tbi)         # per-sample per-chrom GVCF（临时）
├── gvcf_merged/{sample}.g.vcf.gz(.tbi)          # per-sample 全基因组 GVCF
├── genomicsdb/{chrom}/                          # GenomicsDB 数据库
├── vcf/
│   ├── genotyped/{chrom}.vcf.gz(.tbi)           # 联合基因分型结果
│   ├── filtered/
│   │   ├── {chrom}.snp.filtered.vcf.gz(.tbi)    # SNP 过滤后 VCF
│   │   └── {chrom}.indel.filtered.vcf.gz(.tbi)  # INDEL 过滤后 VCF
│   └── final/
│       └── all_samples.filtered.vcf.gz(.tbi)    ★ 最终结果（SNP + INDEL 全基因组合并）
├── qc/
│   ├── fastp/{sample}/{sample}.html(.json)      # fastp 质控报告
│   ├── bam/{sample}.{flagstat,idxstats,stats}.txt
│   ├── markdup/{sample}.markdup_metrics.txt
│   ├── mosdepth/{sample}.*                      # 深度统计文件
│   └── multiqc/multiqc_report.html              ★ QC 总报告
└── logs/                                        # 所有步骤日志
```

---

## 注意事项

- **BWA 索引**：若参考基因组旁已有 `.bwt` 等索引文件，`bwa_index` rule 会自动跳过构建
- **并行粒度**：HaplotypeCaller 按样本 × 染色体并行，GenomicsDB/GenotypeGVCFs 按染色体并行，可充分利用集群资源
- **临时文件**：排序后 BAM 和 per-chrom GVCF 标记为 `temp()`，流程完成后自动清理
- **过滤策略**：默认使用 GATK Hard Filter；若有已知变异集（dbSNP、Mills 等），可在 `config.yaml` 中将 `snp_filter_expr` / `indel_filter_expr` 替换为 VQSR 流程
- **GATK 索引**：GenotypeGVCFs、过滤、合并等步骤均通过独立的 `tabix` rule 强制构建 `.tbi` 索引，避免下游步骤找不到索引文件
- **参考基因组**：默认使用 GRCh38（染色体命名无 `chr` 前缀），若使用 hg38（有 `chr` 前缀），需同步修改 `config/chroms.txt`
