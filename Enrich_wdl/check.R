library(optparse)
option_list <- list(
    make_option(c("-g", "--genus"), type = "character", default = "Arachis",
                            help = "Genus name [default %default]"),
    make_option(c("-s", "--species"), type = "character", default = "hypogaea",
                            help = "Species name [default %default]"),
    make_option(c("-c", "--gene_csv"), type = "character", default = "/data/work/findmarker/allmarkers_test.csv",
                            help = "Path to gene CSV file [default %default]"),
    make_option(c("-d", "--DB"), type = "character", default = "/data/work/0.peanut/orgdb/org.Ahypogaea.eg.db",
                            help = "Path to orgdb database [default %default]")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
genus <- opt$genus
species <- opt$species
gene_csv <- opt$gene_csv
DB <- opt$DB

# library
db_name <- paste0("org.", substr(genus, 1, 1), species, ".eg.db")
print(db_name)
#install.packages(paste0(db_name,"/"), repos = NULL, type = "sources")
install.packages(DB, repos = NULL, type = "sources")
do.call(library, list(db_name))
db <- get(db_name)
db_gid <- keys(db, keytype = "GID")
# gene_set
gene_set <- read.csv(gene_csv, header = TRUE, sep = ",")
head(gene_set)
common_genes <- gene_set$gene[gene_set$gene %in% db_gid]
num_common_genes <- length(common_genes)
total_genes <- length(gene_set$gene)
percentage <- (num_common_genes / total_genes) * 100
# print result of alignment
sink("log.txt")
cat("数据库中的前10个 GID:\n")
print(head(db_gid, 10))
cat("数据库中 GID 的总数:", length(db_gid), "\n")
cat("\nhead(gene_set$gene):\n")
print(head(gene_set$gene, 10))
cat("gene_csv的总基因数量:", total_genes, "\n")
cat("\n存在于db_gid中的基因数量:", num_common_genes, "\n")
cat("gene_csv基因比对到db_gid的百分占比:", round(percentage, 2), "%\n")
sink()
cat("日志已保存到log.txt文件中。\n")
############################ make gene map to gid ##############################
# 使用正则表达式匹配以 "arahy.Tifrunner.gnm2.ann2." 开头的字符串，并在其最后追加 ".1"
gene_set$gene <- ifelse(grepl("^arahy\\.Tifrunner\\.gnm2\\.ann2\\.", gene_set$gene), 
                        paste0(gene_set$gene, ".1"), 
                        gene_set$gene)
head(gene_set$gene)
gene_set$gene_id <- gene_set$gene
head(gene_set)
################################################################################
write.csv(gene_set, "preprocess.csv", row.names = FALSE)
common_genes <- gene_set$gene[gene_set$gene %in% db_gid]
num_common_genes <- length(common_genes)
total_genes <- length(gene_set$gene)
percentage <- (num_common_genes / total_genes) * 100
# print result of alignment
sink("log2.txt")
cat("数据库中的前10个 GID:\n")
print(head(db_gid, 10))
cat("数据库中 GID 的总数:", length(db_gid), "\n")
cat("\nhead(gene_set$gene):\n")
print(head(gene_set$gene, 10))
cat("gene_csv的总基因数量:", total_genes, "\n")
cat("\n存在于db_gid中的基因数量:", num_common_genes, "\n")
cat("gene_csv基因比对到db_gid的百分占比:", round(percentage, 2), "%\n")
sink()
cat("日志已保存到log2.txt文件中。\n")