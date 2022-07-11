## Function: enrichment analysis of genes 
# -----------------------------------------
# Author: Cao Luolong, caoluolong@outlook.com
# Date created: 09-04-2021
# Date updated: 2022.07.05
# Warning: 
# @HEU.@FUDAN.
# -----------------------------------------

rm(list=ls())
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.15")
# BiocManager::install("clusterProfiler")  #for enrichment
# BiocManager::install("topGO")  #for plot
# BiocManager::install("Rgraphviz")
# BiocManager::install("pathview") #for KEGG pathway
# BiocManager::install("org.Hs.eg.db") #for gene annotation.

library(BiocManager)
library(clusterProfiler)
library(topGO)
library(Rgraphviz)
library(pathview)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)

dir.root <- dirname(dirname(rstudioapi::getActiveDocumentContext()$path))
dir.work <- paste0(dir.root,'/file_meta/F_HC_cort_/')
out_probe_join11 <- read.csv(file.path(paste0(dir.work,'PLS_out_genes_F_HC_cort_PLS_variance_max.csv')))

egenes.symble <- as.character(out_probe_join11$geneSymbol[which(out_probe_join11$PLS3Z>=1.96)])#set a threshhood of gene-Z score
file_name <- paste0(dir.work,'F_HC_cort_neg_')
# egenes.symble <- as.character(out_probe_join11$geneSymbol[which(out_probe_join11$PLS3Z<=-1.96)])#set a threshhood of gene-Z score

# gene_num = round(dim(out_probe_join11)[1]/100)
# egenes.symble <- as.character(out_probe_join11$geneSymbol[1:gene_num])#set a threshhood of proporation
# egenes.symble <- as.character(out_probe_join11$geneSymbol[(15408-gene_num):15408])#set a threshhood of proporation

egenes.entrez <- mapIds(x = org.Hs.eg.db, keys = egenes.symble, keytype = "SYMBOL", column="ENTREZID")
egenes.entrez <- na.omit(egenes.entrez)

#######
# columns(org.Hs.eg.db)
erich.go.BP = enrichGO(gene = egenes.entrez,
                       OrgDb = org.Hs.eg.db,
                       keyType = "ENTREZID",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       # pAdjustMethod = 'BH',
                       qvalueCutoff = 0.2,
                       readable = TRUE)

erich.go.CC = enrichGO(gene = egenes.entrez,
                       OrgDb = org.Hs.eg.db,
                       keyType = "ENTREZID",
                       ont = "CC",
                       pvalueCutoff = 0.05,
                       # pAdjustMethod = 'BH',
                       qvalueCutoff = 0.2,
                       readable = TRUE)

erich.go.MF = enrichGO(gene = egenes.entrez,
                       OrgDb = org.Hs.eg.db,
                       keyType = "ENTREZID",
                       ont = "MF",
                       pvalueCutoff = 0.05,
                       # pAdjustMethod = 'BH',
                       qvalueCutoff = 0.2,
                       readable = TRUE)

enrich.kegg <- enrichKEGG(
  egenes.entrez,# a vector of entrez gene id.
  organism = "hsa",
  keyType = "kegg", # one of "kegg", 'ncbi-geneid', 'ncib-proteinid' and 'uniprot'
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  pAdjustMethod = "BH",
  # universe,
  minGSSize = 10,
  maxGSSize = 500,
  use_internal_data = FALSE
)

##########
#BP
write.table(erich.go.BP@result, file = file.path(paste0(file_name,'BP.csv')),row.names = FALSE, quote = FALSE,sep = ',')
dotplot(erich.go.BP, showCategory=20, title = 'Biological process')
# barplot(erich.go.BP, title = 'Biological process')
#CC
write.table(erich.go.CC@result, file = file.path(paste0(file_name,'CC.csv')),row.names = FALSE, quote = FALSE,sep = ',')
dotplot(erich.go.CC, showCategory=20, title = 'Cellular component')
#MF
write.table(erich.go.MF@result, file = file.path(paste0(file_name,'MF.csv')),row.names = FALSE, quote = FALSE,sep = ',')
dotplot(erich.go.MF, showCategory=20, title = 'Molecular function')
#kegg
write.table(enrich.kegg@result, file = file.path(paste0(file_name,'kegg_z1.96.csv')),row.names = FALSE, quote = FALSE,sep = ',')
pic_dot <- dotplot(enrich.kegg, showCategory=20, title = 'KEGG Enrichment')
pic_dot
pdf(file = file.path(paste0(file_name,'kegg_Z1.96.pdf')),width=7,height=7,fonts = NULL,family="Times")
pic_dot
dev.off()

######
go <- enrichGO(egenes.entrez,
               OrgDb = org.Hs.eg.db, 
               ont='ALL',
               pAdjustMethod = 'BH',
               pvalueCutoff = 0.05, 
               qvalueCutoff = 0.2,
               keyType = 'ENTREZID')

pic_bar <- barplot(go, split="ONTOLOGY")+ facet_grid(ONTOLOGY~.,scale="free")
pic_dot <- dotplot(go, split="ONTOLOGY")+ facet_grid(ONTOLOGY~.,scale="free")
# print to pdf
pic_dot
pdf(file = file.path(paste0(file_name,'goz1.96.pdf')),width=10,height=7,family="Times")
pic_dot
dev.off()
write.table(go@result, file = file.path(paste0(file_name,'goZ1.96.csv')),
            row.names = FALSE, quote = FALSE,sep = ',')


# compare enrich result under different Z-threshold.
kegg1 <- read.csv(file = 'F:\\work_dir\\AHBAenrich\\file_meta\\F_HC_cort_\\F_HC_cort_pos_kegg_z1.96_new.csv')
kegg2 <- read.csv(file = 'F:\\work_dir\\AHBAenrich\\file_meta\\F_MDD_cort_\\F_MDD_cort_pos_kegg_z1.96_new.csv')
inr2 <- left_join(select(kegg1,ID,Description,p.adjust),select(kegg2,ID,Description,p.adjust),by = 'ID') %>%filter(as.numeric(p.adjust.y)>=0)
write.table(inr2,file = "F:\\work_dir\\AHBAenrich\\file_meta\\F_HC_cort_\\pos_keggZ1.96_integ_HC+MDD.csv",row.names = FALSE,sep = ",")
