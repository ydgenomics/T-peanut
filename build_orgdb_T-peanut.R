# Title: build_orgdb.R
# Date: 2025-04-10
# Description: This script is used to build orgdb database for peanut
# Output: org.Ahypogaea.eg.db

library(clusterProfiler)
library(tidyverse)
library(stringr)
library(KEGGREST)
library(AnnotationForge)
library(clusterProfiler)
library(dplyr)
library(jsonlite)
library(purrr)
library(RCurl)
library(data.table)
library(readxl)
library(optparse)

option_list <- list(
    make_option(c("-e", "--emapper_xlsx"), type = "character", default = "/data/work/peanut/orgdb/out.emapper.annotations.xlsx", help = "Path to emapper annotations xlsx file", metavar = "character"),
    #make_option(c("-k", "--kojson"), type = "character", default = "ko.json", help = "Path to KEGG ko JSON file", metavar = "character"),
    make_option(c("-t", "--taxid"), type = "character", default = "3818", help = "Taxonomy ID", metavar = "character"),
    make_option(c("-g", "--genus"), type = "character", default = "Arachis", help = "Genus name", metavar = "character"),
    make_option(c("-s", "--species"), type = "character", default = "hypogaea", help = "Species name", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
emapper_annotations_xlsx <- opt$emapper_xlsx
#ko_json <- opt$kojson # Download the ko.json file from KEGG website[https://www.kegg.jp/kegg-bin/download_htext?htext=ko00001&format=json&filedir=]
tax_id <- opt$taxid
genus <- opt$genus
species <- opt$species

emapper <- read_excel(emapper_annotations_xlsx, skip = 2) # Front 2 rows are annotations
head(emapper)
options(stringsAsFactors = F)
emapper[emapper==""] <- NA

# gene_info
gene_info <- emapper %>% 
  dplyr::select(GID = query,GENENAME = Preferred_name) %>% 
  na.omit()
head(gene_info)

# gene2go
gos <- emapper %>% 
  dplyr::select(query, GOs) %>% 
  na.omit()
head(gos)
gene2go = data.frame(GID =character(),
                     GO = character(),
                     EVIDENCE = character())
setDT(gos)
gene2go <- gos[, {
  the_gid <- query
  the_gos <- unlist(str_split(GOs, ","))
  data.table(
    GID = rep(the_gid, length(the_gos)),
    GO = the_gos,
    EVIDENCE = rep("IEA", length(the_gos))
  )
}, by = 1:nrow(gos)]
gene2go <- gene2go[, c("GID", "GO", "EVIDENCE"), drop = FALSE]
gene2go$GO[gene2go$GO=="-"] <- NA
gene2go <- na.omit(gene2go)
head(gene2go)

# gene2ko
gene2ko <- emapper %>% 
  dplyr::select(GID = query, Ko = KEGG_ko) %>% 
  na.omit()
gene2ko$Ko <- gsub("ko:","",gene2ko$Ko)
head(gene2ko)

# gene2pathway
#update_kegg <- function(json = ko_json) {
#  pathway2name <- tibble(Pathway = character(), Name = character())
# #  ko2pathway <- tibble(Ko = character(), Pathway = character())
#   kegg <- fromJSON(json)
#   for (a in seq_along(kegg[["children"]][["children"]])) {
#     A <- kegg[["children"]][["name"]][[a]]
#     for (b in seq_along(kegg[["children"]][["children"]][[a]][["children"]])) {
#       B <- kegg[["children"]][["children"]][[a]][["name"]][[b]] 
#       for (c in seq_along(kegg[["children"]][["children"]][[a]][["children"]][[b]][["children"]])) {
#         pathway_info <- kegg[["children"]][["children"]][[a]][["children"]][[b]][["name"]][[c]]
#         pathway_id <- str_match(pathway_info, "ko[0-9]{5}")[1]
#         pathway_name <- str_replace(pathway_info, " \\[PATH:ko[0-9]{5}\\]", "") %>% str_replace("[0-9]{5} ", "")
#         pathway2name <- rbind(pathway2name, tibble(Pathway = pathway_id, Name = pathway_name))
#         kos_info <- kegg[["children"]][["children"]][[a]][["children"]][[b]][["children"]][[c]][["name"]]
#         kos <- str_match(kos_info, "K[0-9]*")[,1]
#         ko2pathway <- rbind(ko2pathway, tibble(Ko = kos, Pathway = rep(pathway_id, length(kos))))
#       }
#     }
#   }
#   save(pathway2name, ko2pathway, file = "kegg_info.RData")
# }
# update_kegg(ko_json)
load(file = "/script/build_orgdb/kegg_info.RData") # stored in the image GO-R--02
gene2pathway <- gene2ko  %>% 
  left_join(ko2pathway,by = "Ko") %>% 
  dplyr::select(GID, Pathway) %>% 
  na.omit()

# delete duplication
gene2go <- unique(gene2go)
gene2go <- gene2go[!duplicated(gene2go),]
gene2ko <- gene2ko[!duplicated(gene2ko),]
gene2pathway <- gene2pathway[!duplicated(gene2pathway),]


# Check the information of species [https://www.ncbi.nlm.nih.gov/taxonomy]
makeOrgPackage(gene_info=gene_info,
               go=gene2go,
               ko=gene2ko,
               pathway=gene2pathway,
               version="4.0",  #版本
               maintainer = "yd<2144752653@qq.com>",  #修改为你的名字和邮箱
               author = "yd<2144752653@qq.com>",  #修改为你的名字和邮箱
               outputDir = ".",  #输出文件位置
               tax_id=tax_id,
               genus=genus,
               species=species,
               goTable="go")
#install.packages("org.Ahypogaea.eg.db/", repos = NULL, type = "sources")
#columns(org.Ahypogaea.eg.db)