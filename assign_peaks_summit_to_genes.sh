#!/bin/bash

# 用法: ./assign_peaks_to_genes.sh <输入peaks文件完整路径> <基因注释文件完整路径> <输出目录>
# 示例: ./assign_peaks_to_genes.sh /public/home/xinwang/atac/hc/tobias/SAM.peaks.bed /public/home/chaohe/db/geneR1.bed /public/home/chaohe/output

if [ $# -ne 3 ]; then
    echo "用法: $0 <输入peaks文件完整路径> <基因注释文件完整路径> <输出目录>"
    echo "示例: $0 /public/home/xinwang/atac/hc/tobias/SAM.peaks.bed /public/home/chaohe/db/geneR1.bed /public/home/chaohe/output"
    exit 1
fi

INPUT_PEAKS=$1
GENE_BED=$2
OUTPUT_DIR=$3
PREFIX=$(basename "$INPUT_PEAKS" .peaks.bed)

echo "输入peaks文件: $INPUT_PEAKS"
echo "基因注释文件: $GENE_BED"
echo "输出目录: $OUTPUT_DIR"
echo "文件前缀: $PREFIX"

# 创建输出目录
mkdir -p "$OUTPUT_DIR"

# 创建临时工作目录
TEMP_DIR=$(mktemp -d)
echo "临时工作目录: $TEMP_DIR"

# 复制输入文件到临时目录（不修改原始文件）
cp "$INPUT_PEAKS" "$TEMP_DIR/${PREFIX}.peaks.bed"
cp "$GENE_BED" "$TEMP_DIR/gene_annotation.bed"

cd "$TEMP_DIR"

# 步骤1: 检查并添加chr前缀到peaks文件（如果需要）
echo "步骤1: 检查并添加chr前缀到peaks文件"

# 统计包含chr的行数
chr_count=$(grep -c "^chr" "${PREFIX}.peaks.bed" || true)
total_count=$(wc -l < "${PREFIX}.peaks.bed" || true)

if [ "$chr_count" -eq "$total_count" ]; then
    echo "所有行都已包含chr前缀，跳过添加"
    cp "${PREFIX}.peaks.bed" "${PREFIX}.peaks_chr.bed"
elif [ "$chr_count" -eq 0 ]; then
    echo "没有行包含chr前缀，正在添加..."
    sed 's/^/chr&/g' "${PREFIX}.peaks.bed" > "${PREFIX}.peaks_chr.bed"
else
    echo "部分行包含chr前缀，统一添加..."
    sed 's/^/chr&/g' "${PREFIX}.peaks.bed" > "${PREFIX}.peaks_chr.bed"
fi

# 步骤2: 使用bedtools找到最近的基因
echo "步骤2: 寻找最近的基因"
bedtools closest -D ref -t all -mdb all -a "${PREFIX}.peaks_chr.bed" -b gene_annotation.bed > "$OUTPUT_DIR/${PREFIX}.txt"

# 步骤3: 处理正链和负链基因
echo "步骤3: 处理正负链基因"
awk '{if($16=="+") print $1"\t"$2"\t"$3"\t"$4"\t"$10+$2"\t"$14"\t"$12"\t"$13"\t"$16"\t"$17"\t"$2+$10-$12}'  "$OUTPUT_DIR/${PREFIX}.txt" | uniq > "$OUTPUT_DIR/${PREFIX}_minus.txt"
awk '{if($16=="-") print $1"\t"$2"\t"$3"\t"$4"\t"$10+$2"\t"$14"\t"$12"\t"$13"\t"$16"\t"$17"\t"$13-$2-$10}'  "$OUTPUT_DIR/${PREFIX}.txt" | uniq > "$OUTPUT_DIR/${PREFIX}_plus.txt"

# 步骤4: 合并结果
echo "步骤4: 合并正负链结果"
cat "$OUTPUT_DIR/${PREFIX}_minus.txt" "$OUTPUT_DIR/${PREFIX}_plus.txt" > "$OUTPUT_DIR/${PREFIX}_gene.peaks.bed"

# 步骤5: 使用R进行注释
echo "步骤5: 使用R进行peak注释"

RSCRIPT=$(cat << 'EOF'
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("请提供文件前缀和输出目录作为参数")
}
prefix <- args[1]
output_dir <- args[2]

cat("R工作目录:", getwd(), "\n")
cat("文件前缀:", prefix, "\n")
cat("输出目录:", output_dir, "\n")

# 构建文件路径
input_file <- file.path(output_dir, paste0(prefix, "_gene.peaks.bed"))
output_file <- file.path(output_dir, paste0(prefix, "_gene_peak.bed"))

cat("输入文件路径:", input_file, "\n")
cat("输出文件路径:", output_file, "\n")

# 检查文件是否存在
cat("检查输入文件是否存在...\n")
if (!file.exists(input_file)) {
  cat("当前目录文件列表:\n")
  print(list.files(output_dir, pattern = paste0(prefix, ".*")))
  stop("输入文件不存在: ", input_file)
}

cat("输入文件存在，文件大小:", file.info(input_file)$size, "bytes\n")

# 读取数据
# 在读取数据后添加列数检测
cat("读取数据...\n")
x <- read.table(input_file, sep="\t", header=FALSE, stringsAsFactors=FALSE)
cat("数据维度:", dim(x), "\n")

# 获取原始列数
original_cols <- ncol(x)
cat("原始列数:", original_cols, "\n")

# 显示所有列的前几个值来调试
cat("各列示例值:\n")
for(i in 1:original_cols) {
  cat("列", i, ":", head(x[,i], 3), "\n")
}

# 修正列索引计算
col_bedtools_distance <- original_cols - 1    # 距离列（倒数第二列）
col_tss_distance <- original_cols          # 链方向列（最后一列）
# 计算位置列可能不存在，使用距离列代替或者重新计算

cat("修正后的列索引:\n")
cat("bedtools距离列:", col_bedtools_distance, "\n")
cat("TSS距离列:", col_tss_distance, "\n")

if (nrow(x) == 0) {
  stop("读取的数据为空")
}

# 使用TSS距离的绝对值进行排序
cat("计算距离...\n")
x$fz <- abs(x[, col_tss_distance])

# 排序并去重
cat("排序和去重...\n")
x <- x[order(x[,1], x[,4], x$fz),]  # 使用fz列排序
x <- x[!duplicated(x[,4]),]
x <- x[, -ncol(x)]  # 删除最后一列（fz列）
x <- x[order(x[,1], x[,2], x[, col_tss_distance]),]  # 使用距离列排序

# 分类peak - 使用修正的列索引
cat("分类peak...\n")
genebody <- unique(x[which((x[, col_bedtools_distance]==0) & (abs(x[, col_tss_distance])>1500)),])
gene_downstream <- unique(x[which((x[, col_bedtools_distance] != 0) & (x[, col_tss_distance]>0)),])
promoter <- unique(x[which(((x[, col_bedtools_distance]!=0) & (x[, col_tss_distance]>= -3500)  & (x[, col_tss_distance] <= 0)) | ((x[, col_bedtools_distance] == 0)  & (abs(x[, col_tss_distance]) <= 1500))),])
intergeric <- unique(x[which((x[, col_bedtools_distance]!=0) & (x[, col_tss_distance]< -3500)),])
intergeric <- rbind(gene_downstream, intergeric)
genebody$type<-c(rep("genebody",nrow(genebody)))
promoter$type<-c(rep("promoter",nrow(promoter)))  
intergeric$type<-c(rep("intergeric",nrow(intergeric)))  
T<-rbind(genebody, promoter, intergeric)
# 合并结果
T <- rbind(genebody, promoter, intergeric)

cat("注释完成!\n")
cat("各类型peak数量:\n")
cat("genebody:", nrow(genebody), "\n")
cat("promoter:", nrow(promoter), "\n")
cat("intergeric:", nrow(intergeric), "\n")
cat("总计:", nrow(T), "\n")

# 输出结果
write.table(T, output_file, row.names=FALSE, sep="\t", quote=FALSE, col.names=FALSE)
cat("结果已保存到:", output_file, "\n")

# 验证输出文件
if (file.exists(output_file)) {
  cat("输出文件创建成功，大小:", file.info(output_file)$size, "bytes\n")
} else {
  cat("警告: 输出文件未创建\n")
}
EOF
)

# 运行R脚本，传递输出目录路径
echo "运行R分析..."
cd "$OUTPUT_DIR"
echo "当前工作目录: $(pwd)"
echo "目录内容:"
ls -la

echo "$RSCRIPT" | R --vanilla --args "$PREFIX" "$OUTPUT_DIR"

# 清理临时目录
echo "清理临时目录..."
rm -rf "$TEMP_DIR"

# 清理中间文件（可选）
echo "清理中间文件..."
rm -f "$OUTPUT_DIR/${PREFIX}_gene.peaks.bed"

# 检查最终输出
if [ -f "$OUTPUT_DIR/${PREFIX}_gene_peak.bed" ]; then
    echo "最终结果文件: $OUTPUT_DIR/${PREFIX}_gene_peak.bed"
    echo "结果文件行数: $(wc -l < "$OUTPUT_DIR/${PREFIX}_gene_peak.bed")"
else
    echo "警告: 最终输出文件未生成"
    echo "当前输出目录内容:"
    ls -la "$OUTPUT_DIR"
fi

echo "处理完成!"
