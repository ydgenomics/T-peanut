# Title: volcano_visual_T-peanut.R
# Date: 2025-04-10

library(ggplot2)
library(dplyr)
library(ggrepel)

markers <- read.csv("/data/work/peanut/DEA/markers_peanut.csv", header = TRUE, stringsAsFactors = FALSE)
head(markers)
setwd("/data/work/peanut/DEA/volcanos")
# table: avg_log2FC, p_val_adj, gene_id, cluster

plot_volcano <- function(data, coef_col, pval_col, coef_threshold, pval_threshold, title, filename) {
  # 添加 'Significance' 列
  data <- data %>%
    mutate(Significance = case_when(
      (!!sym(coef_col) > coef_threshold & !!sym(pval_col) < pval_threshold) ~ "Significant Up",
      (!!sym(coef_col) < -coef_threshold & !!sym(pval_col) < pval_threshold) ~ "Significant Down",
      TRUE ~ "Not Significant"
    ))
  
  # 处理零 p-value
  zero_pval_genes <- filter(data, !!sym(pval_col) == 0)
  if (nrow(zero_pval_genes) > 0) {
    min_nonzero_pval <- min(filter(data, !!sym(pval_col) > 0)[[pval_col]], na.rm = TRUE)
    if (is.na(min_nonzero_pval)) {
      min_nonzero_pval <- 1e-10
    }
    data <- mutate(data, !!sym(pval_col) := ifelse(!!sym(pval_col) == 0, min_nonzero_pval, !!sym(pval_col)))
  }
  
  # 计算 -log10(p-value)
  data <- mutate(data, neg_log_pval = -log(!!sym(pval_col), base = 10))
  
  # 绘制火山图
  ggplot(data, aes(x = !!sym(coef_col), y = neg_log_pval, color = Significance)) +
    geom_point(alpha = 0.6) +
    geom_vline(xintercept = c(-coef_threshold, coef_threshold), linetype = "dashed", color = "black", alpha = 0.5) +
    geom_hline(yintercept = -log(pval_threshold, base = 10), linetype = "dashed", color = "black") +
    labs(title = title, x = coef_col, y = "-Log10(P-value)", color = "Significance") +
    theme_minimal() +
    theme(legend.position = "right") +
    ggrepel::geom_text_repel(
      data = data %>%
        filter(Significance != "Not Significant") %>%
        top_n(10, -!!sym(pval_col)),  # 选择 p_val_adj 最低的前十个基因
      aes(label = gene_id),
      size = 3,
      box.padding = 0.35,
      point.padding = 0.5,
      color = "black",
      segment.size = 0.2,
      segment.alpha = 0.4
    )
  
  # 保存到文件
  ggsave(filename, width = 8, height = 6, dpi = 1200)
  
  # 打印零 p-value 的基因
  print("zero pval genes:")
  print(zero_pval_genes)
}

for(i in unique(markers$cluster)){
    marker_subset <- filter(markers, cluster == i)
    #marker_subset <- marker_subset %>% filter(p_val_adj < 0.05)
    plot_volcano(data = marker_subset, coef_col = "avg_log2FC", pval_col = "p_val_adj", coef_threshold = 1, pval_threshold = 0.05, title = paste0("Volcano Plot in ", i), filename = paste0("volcano_plot_", i, ".pdf"))
}