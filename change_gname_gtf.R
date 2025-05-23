library(rtracklayer)
library(tidyverse)
library(Biostrings)

gtf <- import("/data/input/Files/husasa/Ref/arahy.Tifrunner.gnm2.ann2.PVFB.gene_models_main.gtf", format = "gtf") %>%
  as.data.frame()

head(gtf)
# 找到 gene_id 不是以 'arahy.Tifrunner.gnm2.ann2.' 开头的行
filtered_gtf <- gtf %>%
  filter(!grepl("^arahy.Tifrunner.gnm2.ann2.", gene_id))

head(filtered_gtf)
# 提取 transcript_id 和 gene_id，保留唯一的 transcript_id
unique_transcripts <- gtf %>%
  select(transcript_id, gene_id) %>%
  distinct(transcript_id, .keep_all = TRUE)

head(unique_transcripts)
# 查看结果
head(unique_transcripts)
seu <- readRDS("/data/work/0.peanut/convert/peanut_dataget_Anno_concat.cg.rds")
# 查看基因名
gene_names <- rownames(seu)
head(gene_names)

# 创建基因名到转录本 ID 的映射
gene_to_transcript <- setNames(unique_transcripts$transcript_id, unique_transcripts$gene_id)

# 重命名基因名
new_gene_names <- gene_to_transcript[gene_names]
head(new_gene_names)
# 将命名向量转换为 data.frame
new_gene_names_df <- data.frame(
  original_gene_id = names(new_gene_names),
  new_transcript_id = new_gene_names,
  stringsAsFactors = FALSE
)
# 查看转换后的 data.frame
head(new_gene_names_df)

# 检查是否有未匹配的基因名
unmatched_genes <- is.na(new_gene_names)
if (any(unmatched_genes)) {
  warning("Some genes were not matched to a transcript ID.")
  new_gene_names[unmatched_genes] <- gene_names[unmatched_genes]  # 保留未匹配的基因名
}

current_assay <- DefaultAssay(seu)
rownames(seu[[current_assay]]@counts) <- new_gene_names_df$new_transcript_id
rownames(seu[[current_assay]]@data) <- new_gene_names_df$new_transcript_id
saveRDS(seu,"/data/work/0.peanut/convert/peanut_dataget_Anno_concat.cg_cgn.rds")

gtf$gene_id <- gtf$transcript_id
head(gtf)
# 加载必要的包
library(dplyr)

# 创建 attributes 列
gtf$attributes <- paste0('transcript_id "', gtf$transcript_id, '"; gene_id "', gtf$gene_id, '";')

# 选择标准 GTF 文件所需的列
standard_gtf <- gtf %>% 
  select(seqnames, source, type, start, end, score, strand, phase, attributes)

# 保存为标准 GTF 文件
output_path <- "/data/work/0.peanut/GRN/output/updated_gtf_file_standard.gtf"
write.table(standard_gtf, file = output_path, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
gtf2 <- import("/data/work/0.peanut/GRN/output/updated_gtf_file_standard.gtf", format = "gtf") %>%
  as.data.frame()

head(gtf2)