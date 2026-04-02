# =============================================================================
# WGS SNP Calling Pipeline
# =============================================================================
# 流程：fastp → bwa mem → markdup → GATK HaplotypeCaller → GenomicsDBImport
#       → GenotypeGVCFs → VQSR/Hard Filter → MultiQC
# =============================================================================

configfile: "config/config.yaml"

import os
import glob
from pathlib import Path

# ---------------------------------------------------------------------------
# 加载样本信息
# ---------------------------------------------------------------------------
def get_samples(sample_dir):
    """
    自动扫描样本目录，支持如下结构：
        sample_dir/
            sampleA/
                *_1.fq.gz  *_R1.fq.gz  *_R1_001.fq.gz
                *_2.fq.gz  *_R2.fq.gz  *_R2_001.fq.gz
    返回 dict: {sample_name: {"r1": path, "r2": path}}
    """
    samples = {}
    sample_dir = Path(sample_dir)
    for subdir in sorted(sample_dir.iterdir()):
        if not subdir.is_dir():
            continue
        name = subdir.name
        # 匹配 R1 / 1 结尾的 fastq
        r1_list = (sorted(subdir.glob("*_R1_*.fq.gz")) or
                   sorted(subdir.glob("*_R1.fq.gz"))   or
                   sorted(subdir.glob("*_1.fq.gz"))    or
                   sorted(subdir.glob("*1.fq.gz")))
        r2_list = (sorted(subdir.glob("*_R2_*.fq.gz")) or
                   sorted(subdir.glob("*_R2.fq.gz"))   or
                   sorted(subdir.glob("*_2.fq.gz"))    or
                   sorted(subdir.glob("*2.fq.gz")))
        if not r1_list or not r2_list:
            print(f"[WARN] 跳过 {subdir}：未找到配对的 fastq 文件")
            continue
        samples[name] = {"r1": str(r1_list[0]), "r2": str(r2_list[0])}
    if not samples:
        raise ValueError(f"在 {sample_dir} 下未找到任何样本，请检查目录结构")
    return samples

SAMPLES = get_samples(config["sample_dir"])
SAMPLE_NAMES = list(SAMPLES.keys())

# ---------------------------------------------------------------------------
# 加载染色体/interval 列表
# ---------------------------------------------------------------------------
def load_chroms(chrom_file):
    """每行一个染色体名，支持 # 注释行"""
    chroms = []
    with open(chrom_file) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                chroms.append(line.split()[0])
    return chroms

CHROMS = load_chroms(config["chrom_list"])

# ---------------------------------------------------------------------------
# 路径辅助函数
# ---------------------------------------------------------------------------
REF   = config["reference"]
OUTDIR = config["outdir"]

def ref_index_files(ref):
    """bwa 索引文件列表"""
    return [f"{ref}.{ext}" for ext in ["amb", "ann", "bwt", "pac", "sa"]]

# =============================================================================
# rule all：声明最终目标文件
# =============================================================================
rule all:
    input:
        # QC 报告
        expand("{outdir}/qc/fastp/{sample}/{sample}.html", outdir=OUTDIR, sample=SAMPLE_NAMES),
        # BAM 统计
        expand("{outdir}/qc/bam/{sample}.flagstat.txt",  outdir=OUTDIR, sample=SAMPLE_NAMES),
        expand("{outdir}/qc/bam/{sample}.idxstats.txt",  outdir=OUTDIR, sample=SAMPLE_NAMES),
        expand("{outdir}/qc/bam/{sample}.stats.txt",     outdir=OUTDIR, sample=SAMPLE_NAMES),
        # 覆盖度
        expand("{outdir}/qc/mosdepth/{sample}.mosdepth.summary.txt", outdir=OUTDIR, sample=SAMPLE_NAMES),
        # 过滤后最终 VCF（SNP + INDEL 分开）
        expand("{outdir}/vcf/filtered/{chrom}.snp.filtered.vcf.gz",   outdir=OUTDIR, chrom=CHROMS),
        expand("{outdir}/vcf/filtered/{chrom}.indel.filtered.vcf.gz", outdir=OUTDIR, chrom=CHROMS),
        # 合并后的全基因组 VCF
        f"{OUTDIR}/vcf/final/all_samples.filtered.vcf.gz",
        # MultiQC
        f"{OUTDIR}/qc/multiqc/multiqc_report.html"

# =============================================================================
# 1. 构建 BWA 参考基因组索引（checkpointing：已存在则跳过）
# =============================================================================
rule bwa_index:
    input:
        ref = REF
    output:
        ref_index_files(REF)
    log:
        f"{OUTDIR}/logs/bwa_index/bwa_index.log"
    threads: 1
    resources:
        mem_gb = 16,
        walltime = "2:00:00"
    shell:
        """
        if [ -f "{input.ref}.bwt" ]; then
            echo "[INFO] BWA index already exists, skipping." | tee {log}
        else
            {config[software][bwa]} index {input.ref} > {log} 2>&1
        fi
        """

# =============================================================================
# 2. fastp 质控
# =============================================================================
rule fastp:
    input:
        r1 = lambda wc: SAMPLES[wc.sample]["r1"],
        r2 = lambda wc: SAMPLES[wc.sample]["r2"]
    output:
        r1   = f"{OUTDIR}/clean_fastq/{{sample}}/{{sample}}_R1.fq.gz",
        r2   = f"{OUTDIR}/clean_fastq/{{sample}}/{{sample}}_R2.fq.gz",
        html = f"{OUTDIR}/qc/fastp/{{sample}}/{{sample}}.html",
        json = f"{OUTDIR}/qc/fastp/{{sample}}/{{sample}}.json"
    log:
        f"{OUTDIR}/logs/fastp/{{sample}}.log"
    threads: config["threads"]["fastp"]
    resources:
        mem_gb   = 8,
        walltime = "2:00:00"
    params:
        extra = config.get("fastp_extra", "--detect_adapter_for_pe --correction --overrepresentation_analysis")
    shell:
        """
        {config[software][fastp]} \
            --in1 {input.r1} --in2 {input.r2} \
            --out1 {output.r1} --out2 {output.r2} \
            --html {output.html} --json {output.json} \
            --thread {threads} \
            {params.extra} \
            > {log} 2>&1
        """

# =============================================================================
# 3. BWA-MEM 比对 + samtools sort
# =============================================================================
rule bwa_mem:
    input:
        r1  = rules.fastp.output.r1,
        r2  = rules.fastp.output.r2,
        ref = REF,
        idx = ref_index_files(REF)   # 确保索引存在
    output:
        bam = temp(f"{OUTDIR}/bam/sorted/{{sample}}.sorted.bam")
    log:
        f"{OUTDIR}/logs/bwa_mem/{{sample}}.log"
    threads: config["threads"]["bwa"]
    resources:
        mem_gb   = 32,
        walltime = "6:00:00"
    params:
        rg     = "@RG\\tID:{sample}\\tSM:{sample}\\tPL:ILLUMINA\\tLB:{sample}",
        extra  = config.get("bwa_extra", "-K 100000000 -Y"),
        sort_t = config["threads"].get("sort", 4)
    shell:
        """
        {config[software][bwa]} mem \
            -t {threads} \
            {params.extra} \
            -R '{params.rg}' \
            {input.ref} {input.r1} {input.r2} \
            2> {log} | \
        {config[software][samtools]} sort \
            -@ {params.sort_t} \
            -o {output.bam} \
            - >> {log} 2>&1
        """

# =============================================================================
# 4. 标记重复（sambamba markdup）
# =============================================================================
rule markdup:
    input:
        bam = rules.bwa_mem.output.bam
    output:
        bam     = f"{OUTDIR}/bam/markdup/{{sample}}.markdup.bam",
        metrics = f"{OUTDIR}/qc/markdup/{{sample}}.markdup_metrics.txt"
    log:
        f"{OUTDIR}/logs/markdup/{{sample}}.log"
    threads: config["threads"]["markdup"]
    resources:
        mem_gb   = 32,
        walltime = "4:00:00"
    shell:
        """
        {config[software][sambamba]} markdup \
            -t {threads} \
            --hash-table-size 1000000 \
            --overflow-list-size 1000000 \
            {input.bam} {output.bam} \
            2> {output.metrics}
        mv {output.bam}.log {log} 2>/dev/null || true
        """

# =============================================================================
# 5. BAM 索引
# =============================================================================
rule index_bam:
    input:
        bam = rules.markdup.output.bam
    output:
        bai = f"{OUTDIR}/bam/markdup/{{sample}}.markdup.bam.bai"
    threads: 4
    resources:
        mem_gb   = 8,
        walltime = "1:00:00"
    shell:
        "{config[software][sambamba]} index -t {threads} {input.bam}"

# =============================================================================
# 6. BAM QC：idxstats / flagstat / stats
# =============================================================================
rule bam_stats:
    input:
        bam = rules.markdup.output.bam,
        bai = rules.index_bam.output.bai
    output:
        idxstats = f"{OUTDIR}/qc/bam/{{sample}}.idxstats.txt",
        flagstat = f"{OUTDIR}/qc/bam/{{sample}}.flagstat.txt",
        stats    = f"{OUTDIR}/qc/bam/{{sample}}.stats.txt"
    log:
        f"{OUTDIR}/logs/bam_stats/{{sample}}.log"
    threads: 4
    resources:
        mem_gb   = 8,
        walltime = "2:00:00"
    params:
        ref = REF
    shell:
        """
        {config[software][samtools]} idxstats  {input.bam} > {output.idxstats} 2>> {log}
        {config[software][samtools]} flagstat   {input.bam} > {output.flagstat} 2>> {log}
        {config[software][samtools]} stats \
            --reference {params.ref} {input.bam} > {output.stats} 2>> {log}
        """

# =============================================================================
# 7. 覆盖度统计（mosdepth）
# =============================================================================
rule mosdepth:
    input:
        bam = rules.markdup.output.bam,
        bai = rules.index_bam.output.bai
    output:
        global_dist = f"{OUTDIR}/qc/mosdepth/{{sample}}.mosdepth.global.dist.txt",
        summary     = f"{OUTDIR}/qc/mosdepth/{{sample}}.mosdepth.summary.txt",
        regions     = f"{OUTDIR}/qc/mosdepth/{{sample}}.regions.bed.gz"
    log:
        f"{OUTDIR}/logs/mosdepth/{{sample}}.log"
    threads: config["threads"].get("mosdepth", 4)
    resources:
        mem_gb   = 16,
        walltime = "3:00:00"
    params:
        prefix     = f"{OUTDIR}/qc/mosdepth/{{sample}}",
        thresholds = config.get("mosdepth_thresholds", "1,5,10,20,30"),
        window     = config.get("mosdepth_window", 500)
    shell:
        """
        {config[software][mosdepth]} \
            --threads {threads} \
            --fast-mode \
            --by {params.window} \
            --thresholds {params.thresholds} \
            --no-abbrev \
            {params.prefix} \
            {input.bam} \
            > {log} 2>&1
        """

# =============================================================================
# 8. GATK HaplotypeCaller（按染色体并行，输出 GVCF）
# =============================================================================
rule haplotype_caller:
    input:
        bam = rules.markdup.output.bam,
        bai = rules.index_bam.output.bai,
        ref = REF
    output:
        gvcf = temp(f"{OUTDIR}/gvcf/{{sample}}/{{chrom}}.g.vcf.gz"),
        tbi  = temp(f"{OUTDIR}/gvcf/{{sample}}/{{chrom}}.g.vcf.gz.tbi")
    log:
        f"{OUTDIR}/logs/haplotype_caller/{{sample}}/{{chrom}}.log"
    threads: 2
    resources:
        mem_gb   = 16,
        walltime = "12:00:00"
    params:
        java_opts = config.get("gatk_java_opts", "-Xmx14g -XX:+UseParallelGC")
    shell:
        """
        {config[software][gatk]} --java-options "{params.java_opts}" \
            HaplotypeCaller \
            -R {input.ref} \
            -I {input.bam} \
            -L {wildcards.chrom} \
            -O {output.gvcf} \
            -ERC GVCF \
            --native-pair-hmm-threads {threads} \
            > {log} 2>&1
        """

# =============================================================================
# 9. 合并单样本各染色体 GVCF（CombineGVCFs per sample）
# =============================================================================
rule combine_sample_gvcf:
    input:
        gvcfs = expand(
            "{outdir}/gvcf/{sample}/{chrom}.g.vcf.gz",
            outdir=OUTDIR, chrom=CHROMS, allow_missing=True
        ),
        ref   = REF
    output:
        gvcf = f"{OUTDIR}/gvcf_merged/{{sample}}.g.vcf.gz",
        tbi  = f"{OUTDIR}/gvcf_merged/{{sample}}.g.vcf.gz.tbi"
    log:
        f"{OUTDIR}/logs/combine_sample_gvcf/{{sample}}.log"
    threads: 4
    resources:
        mem_gb   = 32,
        walltime = "6:00:00"
    params:
        java_opts = config.get("gatk_java_opts", "-Xmx28g"),
        variants  = lambda wc, input: " ".join(f"--variant {v}" for v in input.gvcfs)
    shell:
        """
        {config[software][gatk]} --java-options "{params.java_opts}" \
            CombineGVCFs \
            -R {input.ref} \
            {params.variants} \
            -O {output.gvcf} \
            > {log} 2>&1
        """

# =============================================================================
# 10. GenomicsDBImport（按染色体，多样本联合）
# =============================================================================
rule genomics_db_import:
    input:
        gvcfs = expand(
            "{outdir}/gvcf_merged/{sample}.g.vcf.gz",
            outdir=OUTDIR, sample=SAMPLE_NAMES
        )
    output:
        db = directory(f"{OUTDIR}/genomicsdb/{{chrom}}")
    log:
        f"{OUTDIR}/logs/genomics_db_import/{{chrom}}.log"
    threads: 4
    resources:
        mem_gb   = 32,
        walltime = "12:00:00"
    params:
        java_opts = config.get("gatk_java_opts", "-Xmx28g"),
        variants  = lambda wc, input: " ".join(f"-V {v}" for v in input.gvcfs),
        tmp_dir   = f"{OUTDIR}/tmp"
    shell:
        """
        mkdir -p {params.tmp_dir}
        # 若已存在则先删除（snakemake rerun 时避免冲突）
        rm -rf {output.db}

        {config[software][gatk]} --java-options "{params.java_opts}" \
            GenomicsDBImport \
            {params.variants} \
            --genomicsdb-workspace-path {output.db} \
            -L {wildcards.chrom} \
            --tmp-dir {params.tmp_dir} \
            --reader-threads {threads} \
            > {log} 2>&1
        """

# =============================================================================
# 11. GenotypeGVCFs（按染色体）
# =============================================================================
rule genotype_gvcfs:
    input:
        db  = rules.genomics_db_import.output.db,
        ref = REF
    output:
        vcf = f"{OUTDIR}/vcf/genotyped/{{chrom}}.vcf.gz",
        tbi = f"{OUTDIR}/vcf/genotyped/{{chrom}}.vcf.gz.tbi"
    log:
        f"{OUTDIR}/logs/genotype_gvcfs/{{chrom}}.log"
    threads: 4
    resources:
        mem_gb   = 32,
        walltime = "8:00:00"
    params:
        java_opts = config.get("gatk_java_opts", "-Xmx28g")
    shell:
        """
        {config[software][gatk]} --java-options "{params.java_opts}" \
            GenotypeGVCFs \
            -R {input.ref} \
            -V gendb://{input.db} \
            -L {wildcards.chrom} \
            -O {output.vcf} \
            > {log} 2>&1
        """

# =============================================================================
# 12. Hard Filter：按染色体拆分 SNP / INDEL 并质控
# =============================================================================
rule filter_snp:
    input:
        vcf = rules.genotype_gvcfs.output.vcf,
        ref = REF
    output:
        vcf = f"{OUTDIR}/vcf/filtered/{{chrom}}.snp.filtered.vcf.gz",
        tbi = f"{OUTDIR}/vcf/filtered/{{chrom}}.snp.filtered.vcf.gz.tbi"
    log:
        f"{OUTDIR}/logs/filter_snp/{{chrom}}.log"
    threads: 2
    resources:
        mem_gb   = 16,
        walltime = "4:00:00"
    params:
        java_opts  = config.get("gatk_java_opts", "-Xmx14g"),
        snp_filter = config.get(
            "snp_filter_expr",
            "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"
        )
    shell:
        """
        TMP_RAW={OUTDIR}/vcf/filtered/{wildcards.chrom}.snp.raw.vcf.gz
        TMP_MARKED={OUTDIR}/vcf/filtered/{wildcards.chrom}.snp.marked.vcf.gz

        # Step 1：提取 SNP
        {config[software][gatk]} --java-options "{params.java_opts}" \
            SelectVariants \
            -R {input.ref} \
            -V {input.vcf} \
            --select-type-to-include SNP \
            -O $TMP_RAW \
            >> {log} 2>&1

        # Step 2：打标记（FILTER 列写入过滤原因，PASS 保留原样）
        {config[software][gatk]} --java-options "{params.java_opts}" \
            VariantFiltration \
            -R {input.ref} \
            -V $TMP_RAW \
            --filter-expression "{params.snp_filter}" \
            --filter-name "SNP_HARD_FILTER" \
            -O $TMP_MARKED \
            >> {log} 2>&1

        # Step 3：仅保留 FILTER == PASS 的位点
        {config[software][gatk]} --java-options "{params.java_opts}" \
            SelectVariants \
            -R {input.ref} \
            -V $TMP_MARKED \
            --exclude-filtered \
            -O {output.vcf} \
            >> {log} 2>&1

        rm -f $TMP_RAW $TMP_RAW.tbi $TMP_MARKED $TMP_MARKED.tbi
        """

rule filter_indel:
    input:
        vcf = rules.genotype_gvcfs.output.vcf,
        ref = REF
    output:
        vcf = f"{OUTDIR}/vcf/filtered/{{chrom}}.indel.filtered.vcf.gz",
        tbi = f"{OUTDIR}/vcf/filtered/{{chrom}}.indel.filtered.vcf.gz.tbi"
    log:
        f"{OUTDIR}/logs/filter_indel/{{chrom}}.log"
    threads: 2
    resources:
        mem_gb   = 16,
        walltime = "4:00:00"
    params:
        java_opts    = config.get("gatk_java_opts", "-Xmx14g"),
        indel_filter = config.get(
            "indel_filter_expr",
            "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0"
        )
    shell:
        """
        TMP_RAW={OUTDIR}/vcf/filtered/{wildcards.chrom}.indel.raw.vcf.gz
        TMP_MARKED={OUTDIR}/vcf/filtered/{wildcards.chrom}.indel.marked.vcf.gz

        # Step 1：提取 INDEL
        {config[software][gatk]} --java-options "{params.java_opts}" \
            SelectVariants \
            -R {input.ref} \
            -V {input.vcf} \
            --select-type-to-include INDEL \
            -O $TMP_RAW \
            >> {log} 2>&1

        # Step 2：打标记
        {config[software][gatk]} --java-options "{params.java_opts}" \
            VariantFiltration \
            -R {input.ref} \
            -V $TMP_RAW \
            --filter-expression "{params.indel_filter}" \
            --filter-name "INDEL_HARD_FILTER" \
            -O $TMP_MARKED \
            >> {log} 2>&1

        # Step 3：仅保留 FILTER == PASS 的位点
        {config[software][gatk]} --java-options "{params.java_opts}" \
            SelectVariants \
            -R {input.ref} \
            -V $TMP_MARKED \
            --exclude-filtered \
            -O {output.vcf} \
            >> {log} 2>&1

        rm -f $TMP_RAW $TMP_RAW.tbi $TMP_MARKED $TMP_MARKED.tbi
        """

# =============================================================================
# 13. 合并所有染色体的过滤后 VCF（SNP + INDEL 合并）
# =============================================================================
rule merge_filtered_vcf:
    input:
        snps   = expand("{outdir}/vcf/filtered/{chrom}.snp.filtered.vcf.gz",   outdir=OUTDIR, chrom=CHROMS),
        indels = expand("{outdir}/vcf/filtered/{chrom}.indel.filtered.vcf.gz", outdir=OUTDIR, chrom=CHROMS)
    output:
        vcf = f"{OUTDIR}/vcf/final/all_samples.filtered.vcf.gz",
        tbi = f"{OUTDIR}/vcf/final/all_samples.filtered.vcf.gz.tbi"
    log:
        f"{OUTDIR}/logs/merge_filtered_vcf/merge.log"
    threads: 4
    resources:
        mem_gb   = 32,
        walltime = "4:00:00"
    params:
        # 按染色体顺序排列所有 vcf（先 snp 后 indel，再用 MergeVcfs 合并）
        java_opts = config.get("gatk_java_opts", "-Xmx28g"),
        all_vcfs  = lambda wc, input: " ".join(
            f"-I {v}" for v in sorted(input.snps) + sorted(input.indels)
        )
    shell:
        """
        {config[software][gatk]} --java-options "{params.java_opts}" \
            MergeVcfs \
            {params.all_vcfs} \
            -O {output.vcf} \
            > {log} 2>&1
        """

# =============================================================================
# 14. MultiQC 汇总报告
# =============================================================================
rule multiqc:
    input:
        fastp_json  = expand("{outdir}/qc/fastp/{sample}/{sample}.json",          outdir=OUTDIR, sample=SAMPLE_NAMES),
        flagstat    = expand("{outdir}/qc/bam/{sample}.flagstat.txt",              outdir=OUTDIR, sample=SAMPLE_NAMES),
        idxstats    = expand("{outdir}/qc/bam/{sample}.idxstats.txt",              outdir=OUTDIR, sample=SAMPLE_NAMES),
        stats       = expand("{outdir}/qc/bam/{sample}.stats.txt",                 outdir=OUTDIR, sample=SAMPLE_NAMES),
        markdup     = expand("{outdir}/qc/markdup/{sample}.markdup_metrics.txt",   outdir=OUTDIR, sample=SAMPLE_NAMES),
        mosdepth    = expand("{outdir}/qc/mosdepth/{sample}.mosdepth.summary.txt", outdir=OUTDIR, sample=SAMPLE_NAMES)
    output:
        html = f"{OUTDIR}/qc/multiqc/multiqc_report.html"
    log:
        f"{OUTDIR}/logs/multiqc/multiqc.log"
    threads: 1
    resources:
        mem_gb   = 8,
        walltime = "1:00:00"
    params:
        outdir   = f"{OUTDIR}/qc/multiqc",
        data_dir = OUTDIR
    shell:
        """
        {config[software][multiqc]} \
            {params.data_dir} \
            --outdir {params.outdir} \
            --filename multiqc_report \
            --interactive \
            --force \
            > {log} 2>&1
        """
