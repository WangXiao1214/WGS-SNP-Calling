#!/bin/bash

usage() {
    echo "Usage: $0 -i <input.vcf.gz> [-o <output.txt>] [-b <bcftools_path>]"
    echo ""
    echo "  -i  输入 VCF 文件（必填，需 bgzip 压缩并建立 .tbi 索引）"
    echo "  -o  输出文件路径（可选，默认：sample_qc_stats.txt）"
    echo "  -b  bcftools 可执行文件路径（可选，默认使用 PATH 中的 bcftools）"
    echo "  -h  显示帮助信息"
    echo ""
    echo "Example:"
    echo "  $0 -i output_annotated.vcf.gz"
    echo "  $0 -i output_annotated.vcf.gz -o my_qc.txt -b /path/to/bcftools"
    exit 1
}

# 默认值
OUTPUT="sample_vcfqc_stats.txt"
BCFTOOLS=$(which bcftools 2>/dev/null || echo "/mnt/beegfs/Research/Users/wangxiao01/tools/bcftools_v1.17/bin/bcftools")
VCF=""

while getopts "i:o:b:h" opt; do
    case ${opt} in
        i) VCF="${OPTARG}" ;;
        o) OUTPUT="${OPTARG}" ;;
        b) BCFTOOLS="${OPTARG}" ;;
        h) usage ;;
        *) usage ;;
    esac
done

# 参数校验
if [[ -z "${VCF}" ]]; then
    echo "Error: 输入 VCF 文件（-i）为必填参数。"
    echo ""
    usage
fi

if [[ ! -f "${VCF}" ]]; then
    echo "Error: 找不到文件 '${VCF}'"
    exit 1
fi

if [[ ! -x "${BCFTOOLS}" ]]; then
    echo "Error: bcftools 不可执行，请通过 -b 参数指定正确路径（当前：${BCFTOOLS}）"
    exit 1
fi

echo "Input VCF  : ${VCF}"
echo "Output     : ${OUTPUT}"
echo "bcftools   : ${BCFTOOLS} ($("${BCFTOOLS}" --version | head -1))"
echo ""

echo -e "Sample\tSNPs\tInDels\tTs\tTv\tTiTv_ratio\tHets_SNP\tHomAlt_SNP\tHet_Hom_ratio\tX_Het\tX_HomAlt\tX_Het_rate\tPredicted_Sex" > "${OUTPUT}"

SAMPLES=$("${BCFTOOLS}" query -l "${VCF}")

for SAMPLE in ${SAMPLES}; do
    echo "[$(date '+%H:%M:%S')] Processing: ${SAMPLE}"

    # ---- 全基因组统计 ----
    STATS=$("${BCFTOOLS}" stats -s "${SAMPLE}" "${VCF}")

    # PSC: id sample nRefHom nNonRefHom nHets nTs nTv nIndels avgDP nSingletons nHapRef nHapAlt nMissing
    PSC_LINE=$(echo "${STATS}" | grep "^PSC")
    HOMALT=$(echo "${PSC_LINE}"  | awk '{print $5}')
    HETS=$(echo "${PSC_LINE}"    | awk '{print $6}')
    TS=$(echo "${PSC_LINE}"      | awk '{print $7}')
    TV=$(echo "${PSC_LINE}"      | awk '{print $8}')
    INDELS=$(echo "${PSC_LINE}"  | awk '{print $9}')

    # SNPs = Hets + HomAlt (PSC 中 nHets/nNonRefHom 均只计 SNPs)
    SNPS=$(( HETS + HOMALT ))

    # Ti/Tv ratio
    TSTV_RATIO=$(awk -v ts="${TS}" -v tv="${TV}" \
        'BEGIN { if (tv+0 > 0) printf "%.4f", ts/tv; else print "NA" }')

    # Het/Hom ratio
    HET_HOM=$(awk -v het="${HETS}" -v hom="${HOMALT}" \
        'BEGIN { if (hom+0 > 0) printf "%.4f", het/hom; else print "NA" }')

    # ---- X 染色体性别预测 ----
    STATS_X=$("${BCFTOOLS}" stats -s "${SAMPLE}" -r X "${VCF}")
    PSC_X=$(echo "${STATS_X}" | grep "^PSC")
    X_HOMALT=$(echo "${PSC_X}" | awk '{print $5}')
    X_HETS=$(echo "${PSC_X}"   | awk '{print $6}')

    # X 杂合率 = nHets / (nHets + nHomAlt)
    X_HET_RATE=$(awk -v het="${X_HETS}" -v hom="${X_HOMALT}" \
        'BEGIN { tot=het+hom; if (tot>0) printf "%.4f", het/tot; else print "NA" }')

    # 判断性别：X 杂合率 >= 0.25 → Female，<= 0.15 → Male，其余 Unknown
    SEX=$(awk -v r="${X_HET_RATE}" \
        'BEGIN {
            if (r=="NA") { print "Unknown" }
            else if (r+0 >= 0.25) { print "Female" }
            else if (r+0 <= 0.15) { print "Male" }
            else { print "Unknown" }
        }')

    echo -e "${SAMPLE}\t${SNPS}\t${INDELS}\t${TS}\t${TV}\t${TSTV_RATIO}\t${HETS}\t${HOMALT}\t${HET_HOM}\t${X_HETS}\t${X_HOMALT}\t${X_HET_RATE}\t${SEX}"
    echo -e "${SAMPLE}\t${SNPS}\t${INDELS}\t${TS}\t${TV}\t${TSTV_RATIO}\t${HETS}\t${HOMALT}\t${HET_HOM}\t${X_HETS}\t${X_HOMALT}\t${X_HET_RATE}\t${SEX}" >> "${OUTPUT}"
done

echo ""
echo "Done! Results saved to: ${OUTPUT}"
