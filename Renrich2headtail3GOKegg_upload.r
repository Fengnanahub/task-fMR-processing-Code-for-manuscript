######--enrich analysis -compute and save (GO/KEGG).csv quickly;

library(clusterProfiler)
library(topGO)
library(Rgraphviz)
library(pathview)
library(org.Hs.eg.db)

out_probe_join11 <- read.csv(file.path('D:\\DATA_Fudan\\SCZ2010data2nd\\SCZanalysisCorrectTrial220413\\argenbrain\\zscorePLSXin\\M2subcortex\\PLS_out_genes_M2_subcort_PLS_variance1.csv'))
egenes.symble <- as.character(out_probe_join11$geneSymbol[which(out_probe_join11$PLS1Z<=-3)])# set threshold 3/-3.

egenes.symble <- as.character((unique(egenes.symble))) # ENSEMBL id
columns(org.Hs.eg.db)

egenes.entrez <- mapIds(x = org.Hs.eg.db, keys = egenes.symble, keytype = "SYMBOL", column="ENTREZID") #trans Entrez gene ID 
egenes.entrez <- na.omit(egenes.entrez)

# KEGG
enrich.kegg <- enrichKEGG(
  egenes.entrez,# a vector of entrez gene id.
  organism = "hsa",
  keyType = "kegg", # one of "kegg", 'ncbi-geneid', 'ncib-proteinid' and 'uniprot'
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  # universe,
  minGSSize = 10,
  maxGSSize = 500,
  qvalueCutoff = 0.05,
  use_internal_data = FALSE
)
# seeres.kegg <- as.data.frame(enrich.kegg@result)
write.table(enrich.kegg@result, file = file.path('D:\\DATA_Fudan\\SCZ2010data2nd\\SCZanalysisCorrectTrial220413\\argenbrain\\zscorePLSXin\\M2subcortex\\negative3\\PLS_out_genes_M2_subcort_PLS_variance1_neg_Kegg.csv'),
            row.names = FALSE, quote = FALSE,sep = ',')

# save all
go <- enrichGO(gene = egenes.entrez,
               OrgDb = org.Hs.eg.db,
               keyType = "ENTREZID",
               ont = "ALL",
               pvalueCutoff = 0.05,
               pAdjustMethod = 'BH',
               qvalueCutoff = 0.05,
               readable = TRUE)
write.table(go@result, file = file.path('D:\\DATA_Fudan\\SCZ2010data2nd\\SCZanalysisCorrectTrial220413\\argenbrain\\zscorePLSXin\\M2subcortex\\negative3\\PLS_out_genes_M2_subcort_PLS_variance1_neg_GO.csv'),
            row.names = FALSE, quote = FALSE,sep = ',')

