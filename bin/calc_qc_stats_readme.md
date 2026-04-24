# 样本 VCF 质控统计说明

## 概述

本文档说明对多样本 VCF 文件（`output_annotated.vcf.gz`）逐样本计算质控指标的方法、来源字段和判断依据。统计脚本为 `calc_qc_stats.sh`，工具版本为 **bcftools v1.17**，输出结果保存于 `sample_qc_stats.txt`（Tab 分隔）。

---

## 运行方式

### 参数说明

| 参数 | 是否必填 | 默认值 | 说明 |
|------|---------|--------|------|
| `-i` | **必填** | — | 输入 VCF 文件（需 bgzip 压缩并建立 `.tbi` 索引） |
| `-o` | 可选 | `sample_qc_stats.txt` | 输出文件路径 |
| `-b` | 可选 | `PATH` 中的 `bcftools` | 指定 bcftools 可执行文件路径 |
| `-h` | — | — | 显示帮助信息并退出 |

### 示例

```bash
# 最简用法（仅指定输入 VCF）
bash calc_qc_stats.sh -i output_annotated.vcf.gz

# 自定义输出文件名
bash calc_qc_stats.sh -i output_annotated.vcf.gz -o cohort_qc.txt

# 指定 bcftools 路径（当环境 PATH 中无 bcftools 时使用）
bash calc_qc_stats.sh -i output_annotated.vcf.gz -o cohort_qc.txt -b /path/to/bcftools

# 后台运行并记录日志
bash calc_qc_stats.sh -i output_annotated.vcf.gz > calc_qc_stats.log 2>&1 &
```

### 运行流程

脚本对 VCF 中的每个样本依次执行两次 `bcftools stats`：
1. **全基因组**统计（获取 SNP、InDel、Ti/Tv、Het/Hom 等指标）
2. **仅限 X 染色体**统计（获取 X 染色体杂合率，用于性别预测）

运行开始时会打印当前参数配置，每处理一个样本会输出带时间戳的进度信息，同时实时追加写入输出文件。

---

## 输出字段说明

| 列名 | 说明 |
|------|------|
| `Sample` | 样本名 |
| `SNPs` | 该样本携带的 SNP 总数（杂合 + 纯合变异） |
| `InDels` | 该样本携带的 InDel 总数 |
| `Ts` | Transitions（转换突变）数量 |
| `Tv` | Transversions（颠换突变）数量 |
| `TiTv_ratio` | Ti/Tv 比值 |
| `Hets_SNP` | SNP 中杂合位点数（0/1） |
| `HomAlt_SNP` | SNP 中纯合变异位点数（1/1） |
| `Het_Hom_ratio` | Het/Hom 比值 |
| `X_Het` | X 染色体上杂合 SNP 数 |
| `X_HomAlt` | X 染色体上纯合变异 SNP 数 |
| `X_Het_rate` | X 染色体杂合率 |
| `Predicted_Sex` | 预测性别（Male / Female / Unknown） |

---

## 各指标计算方法与依据

### 1. SNP 数量（#SNPs）

**计算方法：**

```
#SNPs = nHets + nNonRefHom
```

**来源：** `bcftools stats` 输出中的 `PSC`（Per-Sample Counts）行，第 6 列（`nHets`）和第 5 列（`nNonRefHom`）。

**说明：** PSC 行中的 `nHets` 和 `nNonRefHom` 仅统计 SNP 类型的变异（InDel 另计），因此两者之和即为该样本真实携带的 SNP 总数。注意此处统计的是该样本实际具有变异基因型的位点，而非 VCF 文件中 SNP 记录的总行数。

---

### 2. InDel 数量（#InDels）

**计算方法：**

```
#InDels = nIndels（PSC 第 9 列）
```

**来源：** `PSC` 行第 9 列（`nIndels`）。

**说明：** 该值由 `bcftools stats` 在 per-sample 模式下单独计算，统计该样本中基因型为杂合或纯合变异的 InDel 位点总数（包括插入和缺失）。

---

### 3. Ti/Tv 比值（TiTv_ratio）

**计算方法：**

```
Ti/Tv = Ts / Tv
```

**来源：** `PSC` 行第 7 列（`nTransitions`，Ts）和第 8 列（`nTransversions`，Tv）。

**说明：** Transitions（Ts）指嘌呤之间（A↔G）或嘧啶之间（C↔T）的替换；Transversions（Tv）指嘌呤与嘧啶之间的替换（A/G↔C/T）。由于生物学上 Ts 发生概率约为 Tv 的两倍，真实变异的 Ti/Tv 比值应显著高于随机突变（随机期望值约为 0.5）。

> **注意：** 本脚本使用 `PSC` 行的 Ts/Tv，而非 `TSTV` 行。`TSTV` 行是全局位点统计（不区分样本），对多样本 VCF 每个样本输出值完全相同；`PSC` 行才是每个样本自身基因型对应的转换/颠换计数，反映该样本的真实情况。

**参考范围：**

| 数据类型 | 期望 Ti/Tv |
|---------|-----------|
| WGS 全基因组 | 1.9 ~ 2.1 |
| WES 外显子组 | 2.8 ~ 3.0 |
| 随机突变（理论值） | ~0.5 |

---

### 4. Het/Hom 比值（Het_Hom_ratio）

**计算方法：**

```
Het/Hom = nHets / nNonRefHom
```

**来源：** `PSC` 行第 6 列（`nHets`）和第 5 列（`nNonRefHom`）。

**说明：** 在二倍体基因组中，杂合位点（0/1）数量通常多于纯合变异位点（1/1）。该比值反映了样本基因组的杂合程度，可用于检测样本污染、亲缘关系异常或测序质量问题。

**参考范围（WGS）：**

| 情况 | 期望 Het/Hom |
|------|------------|
| 正常二倍体样本 | 1.5 ~ 2.0 |
| 偏低（< 1.0） | 可能存在样本纯化、近亲繁殖或数据问题 |
| 偏高（> 2.5） | 可能存在样本污染或 CNV 导致的假杂合 |

---

### 5. 性别预测（Predicted_Sex）

**计算方法：**

```
X_Het_rate = X_nHets / (X_nHets + X_nHomAlt)

若 X_Het_rate >= 0.25  → Female
若 X_Het_rate <= 0.15  → Male
否则                   → Unknown
```

**来源：** 仅限 X 染色体（`-r X`）的 `bcftools stats` 输出中 `PSC` 行第 5、6 列。

**原理：** 男性（46,XY）的 X 染色体为半合子（hemizygous），正常情况下 X 染色体上几乎不存在杂合变异（真实杂合位点应为测序错误或伪常染色体区域 PAR 的信号），杂合率极低（通常 < 5%）。女性（46,XX）的 X 染色体为二倍体，杂合率与常染色体相近（通常 > 30%）。

**阈值设置依据：**

| 预测性别 | X 杂合率范围 | 生物学解释 |
|---------|------------|----------|
| Male | ≤ 0.15 | X 染色体为半合子，杂合率极低 |
| Female | ≥ 0.25 | X 染色体为二倍体，杂合率正常 |
| Unknown | 0.15 ~ 0.25 | 介于两者之间，需人工核查 |

> **注意：** X 染色体伪常染色体区域（PAR1/PAR2）会表现出正常二倍体杂合率，可能略微抬高男性样本的 X 杂合率。若需更精确的性别判断，可额外排除 PAR 区域，或结合 Y 染色体覆盖度综合判断。

---

## 工具与参考

- **工具：** [bcftools v1.17](https://samtools.github.io/bcftools/)
- **PSC 字段定义：** `bcftools stats` 输出的 `# PSC` 注释行
  ```
  PSC [id] [sample] [nRefHom] [nNonRefHom] [nHets] [nTs] [nTv] [nIndels] [avgDP] [nSingletons] [nHapRef] [nHapAlt] [nMissing]
  ```
- **参考基因组：** GRCh38（染色体命名无 `chr` 前缀）
