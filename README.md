# WGS SNP Calling Pipeline

## 目录结构

```
.
├── Snakefile                  # 主流程
├── run.sh                     # 启动脚本
├── config/
│   ├── config.yaml            # 主配置（路径、软件、参数）
│   └── chroms.txt             # 染色体列表
└── profiles/
    └── slurm/
        └── config.yaml        # SLURM 集群配置
```

## 流程概览

```
fastp（质控）
  └── bwa mem + samtools sort（比对）
        └── sambamba markdup（标记重复）
              ├── bam_stats（idxstats / flagstat / stats）
              ├── mosdepth（覆盖度）
              └── GATK HaplotypeCaller（按染色体，-ERC GVCF）
                    └── combine_sample_gvcf（单样本 GVCF 合并）
                          └── GenomicsDBImport（按染色体，多样本）
                                └── GenotypeGVCFs（按染色体）
                                      ├── filter_snp（SelectVariants + VariantFiltration）
                                      └── filter_indel（SelectVariants + VariantFiltration）
                                            └── merge_filtered_vcf（全基因组合并）
                                                  └── MultiQC（汇总报告）
```

## 使用方法

### 1. 修改配置

编辑 `config/config.yaml`，至少需要设置：

```yaml
sample_dir:  /path/to/your/samples   # 样本根目录
reference:   /path/to/genome.fa      # 参考基因组
chrom_list:  config/chroms.txt       # 染色体列表
outdir:      results                 # 输出目录
```

样本目录格式（自动识别，无需样本表）：
```
sample_dir/
  sample1/
    xxx_1.fq.gz  xxx_2.fq.gz
  sample2/
    xxx_R1.fq.gz  xxx_R2.fq.gz
```

### 2. 修改染色体列表

编辑 `config/chroms.txt`，确保与参考基因组 fasta header 一致。

### 3. 运行

```bash
# 先做 dry run 检查
bash run.sh dry

# 本地运行
bash run.sh run

# SLURM 集群（先修改 profiles/slurm/config.yaml）
bash run.sh cluster

# 生成流程 DAG 图
bash run.sh dag
```

## 输出结构

```
results/
├── clean_fastq/{sample}/           # fastp 过滤后的 fastq
├── bam/
│   ├── sorted/{sample}.sorted.bam  # 排序后（临时，自动删除）
│   └── markdup/{sample}.markdup.bam
├── gvcf/{sample}/{chrom}.g.vcf.gz  # per-sample per-chrom GVCF（临时）
├── gvcf_merged/{sample}.g.vcf.gz   # per-sample 合并 GVCF
├── genomicsdb/{chrom}/             # GenomicsDB
├── vcf/
│   ├── genotyped/{chrom}.vcf.gz    # 联合 genotype 结果
│   ├── filtered/
│   │   ├── {chrom}.snp.filtered.vcf.gz
│   │   └── {chrom}.indel.filtered.vcf.gz
│   └── final/
│       └── all_samples.filtered.vcf.gz   ★ 最终结果
├── qc/
│   ├── fastp/{sample}/
│   ├── bam/{sample}.{flagstat,idxstats,stats}.txt
│   ├── markdup/{sample}.markdup_metrics.txt
│   ├── mosdepth/{sample}.*
│   └── multiqc/multiqc_report.html       ★ QC 总报告
└── logs/                           # 所有日志
```

## 注意事项

- **BWA 索引**：若参考基因组旁已有 `.bwt` 等索引文件，`bwa_index` rule 会自动跳过构建
- **并行粒度**：HaplotypeCaller 按样本 × 染色体并行，GenomicsDB 按染色体并行，充分利用集群资源
- **临时文件**：排序后的 BAM 和 per-chrom GVCF 标记为 `temp()`，流程完成后自动清理
- **过滤策略**：默认使用 GATK Hard Filter，如有已知变异集可在 `config.yaml` 中替换为 VQSR 参数
