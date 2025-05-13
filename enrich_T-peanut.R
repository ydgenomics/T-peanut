# Title: enrich.R
# Date: 2025-04-10
library(ggplot2)
library(tidyverse)
library(dplyr)
library(clusterProfiler)
#library(org.Ahypogaea.eg.db)
#length(keys(org.Ahypogaea.eg.db))
library(tidyverse)
library(optparse)

option_list <- list(
  make_option(c("--gene_csv"), type = "character", default = "/data/work/peanut/DEA/markers_peanut.csv", help = "input the csv of leiden_0.5"),
  make_option(c("--minp"), type = "numeric", default = 0.05, help = "filter marker gene limited by min pvalue_adj"),
  make_option(c("--species"),type = "character", default = "peanut",help = "differ different species use different database"),
  make_option(c("--db"),type = "character", default = "org.Ahypogaea.eg.db",help = "Name of built db for enrich")
)
opt <- parse_args(OptionParser(option_list = option_list))
do.call(library, list(opt$db)) # 动态加载包

filepath <- paste0(opt$species, "_enrich")
dir.create(filepath)
setwd(filepath)

markers <- read.csv(opt$gene_csv, header = TRUE, stringsAsFactors = FALSE)
head(markers) # gene_id, cluster, p_val_adj
#if (opt$species == "Cer") {markers$gene_id <- sub(".v2.1$", ".1.p", markers$names)}
#if (opt$species == "Pog") {markers$gene_id <- gsub("_", "-", markers$names)}
pathway2gene <- AnnotationDbi::select(org.Ahypogaea.eg.db,keys = keys(org.Ahypogaea.eg.db),columns = c("Pathway","Ko")) %>%
  na.omit() %>%
  dplyr::select(Pathway, GID)
load("/script/build_orgdb/kegg_info.RData")

for(i in unique(markers$cluster)){
    marker_subset <- filter(markers, cluster == i)
    length(marker_subset$gene_id)
    gene_list <- marker_subset %>% filter(p_val_adj < opt$minp)
    gene_list <- gene_list$gene_id
    gene_list <- paste0(gene_list, ".1") # gene_id == GID
    length(gene_list)
    # run enrichGO
    #if (opt$species == "Cer") {
    #data <- enrichGO(gene = gene_list,OrgDb = org.Cthalictroides.eg.db,keyType = 'GID',ont = 'ALL',qvalueCutoff = 0.05,pvalueCutoff = 0.05)}
    #if (opt$species == "Pog") {
    #data <- enrichGO(gene = gene_list,OrgDb = org.Pcirratum.eg.db,keyType = 'GID',ont = 'ALL',qvalueCutoff = 0.05,pvalueCutoff = 0.05)}
    if (opt$species == "peanut") {
    go_data <- enrichGO(gene = gene_list,OrgDb = org.Ahypogaea.eg.db,keyType = 'GID',ont = 'ALL',qvalueCutoff = 0.05,pvalueCutoff = 0.05)
    go_data <- as.data.frame(go_data)
    kegg_result <- enricher(gene_list,TERM2GENE = pathway2gene, TERM2NAME = pathway2name,pvalueCutoff = 0.05,qvalueCutoff = 0.05)
    kegg_data <- as.data.frame(kegg_result)
    dim(kegg_data)
    kegg_data$ONTOLOGY <- "KEGG"
    col_names <- names(kegg_data)
    kegg_data <- kegg_data[, c("ONTOLOGY", col_names[!col_names %in% "ONTOLOGY"])]}
    data <- rbind(go_data, kegg_data)
    if (nrow(data) > 0) {
        print("Data is not empty. Proceeding with the code.")
        length(data$ID)
        data$Name <- paste0(data$ID,"_",data$Description)
        #write.table(data, file = paste0(i,"_enrichGO_data.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
        data_subset <- data %>% arrange(desc(Count)) %>% head(60)
        length(data_subset$ID)
        data_subset <- data_subset %>% mutate(GeneRatio = as.numeric(gsub("/.*", "", GeneRatio)) / as.numeric(gsub(".*/", "", GeneRatio)))
        data_subset$Name <- ifelse(nchar(data_subset$Name) > 100, substr(data_subset$Name, 1, 100), data_subset$Name)
        # viusal1
        if (length(data_subset$ID) > 0) {
        pdf(paste0(i,"_plot1.pdf"))
        plot1 <- ggplot(data_subset, aes(y = GeneRatio, x = reorder(Name, GeneRatio))) + 
            geom_bar(stat = "identity", aes(fill = p.adjust), width = 0.8) +  
            scale_fill_gradient(low = "red", high = "blue") +  
            facet_grid(ONTOLOGY ~ ., scales = "free", space = "free") +  
            coord_flip() + xlab("Name") + ylab("GeneRatio") + labs(title = paste0("Group ", i, " GO and KEGG Enrich")) + 
            theme(
                axis.text.x = element_text(size = 10), 
                axis.text.y = element_text(size = 5), 
                axis.title.x = element_text(size = 12),  
                axis.title.y = element_text(size = 12)) +
            geom_text(aes(label = Count), vjust = 0, size = 1.5) +
            scale_size_continuous(range = c(0.1, 3)) 
        print(plot1)
        dev.off()}
    } else {
        print("Data is empty. Skipping the code.")
    }
}