#ORA analysis
#install packages
install.packages("data.table")
install.packages("ggupset")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biomaRt")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler")

if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("BiocUpgrade") ## you may need this
BiocManager::install("enrichplot")

# load libraries
library(biomaRt)
library(data.table)
library(tidyverse)
library(clusterProfiler)
library(enrichplot)
library(ggupset)
library(devtools)

genes <- fread("/path/to/ref_functional_annot.csv")
go_terms <- fread("/path/to/ref_df_goterms.txt", header = FALSE)

#transfrom to data.frame
genes <- as.data.frame(genes)
go_terms <- as.data.frame(go_terms) %>%
  setNames(c("GOs", "Go_desc"))

genes_godes <- merge(genes, go_terms, by = "GOs", all.x = TRUE)
head(genes_godes)

# remove empty GO terms
genes_godes <- genes_godes %>% 
  filter(!is.na(Go_desc)) # there are 37k genes with no GO terms


#List of genes to analyse
#----
#Candidate genes top 5% fst
top5 <- fread("/path/to/pixy_output/gene_list_fst_top5.txt", header = FALSE)

top5 <- as.data.frame(top5) %>%
  setNames(c("query"))

top5 <- top5 %>% 
  filter(query %in% genes_godes$query)

#ORA analysis
ora_top5_fst <- enricher(gene = top5$query, # your candidate genes
                         TERM2GENE = genes_godes[,c("Go_desc", "query")], # All genes that you use in the analysis.
                         pAdjustMethod = "fdr",
                         pvalueCutoff = 0.1)

# Different plots
barplot(ora_top5_fst)
dotplot(ora_top5_fst)
upsetplot(ora_top5_fst)
cnetplot(ora_top5_fst)

# Save results as dataframe
ora_top5_fst_table <- as.data.frame(ora_top5_fst)
head(ora_top5_fst_table)
write.table(ora_top5_fst_table, "/path/to/pixy_output/ora_fst_top5.txt", sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)

#----
#Candidate genes positive selection
pos_sel <- fread("/path/to/egglib_stats/candidate_genes_pos_sel_Ids.txt", header = FALSE)
head(pos_sel)

pos_sel <- as.data.frame(pos_sel) %>%
  setNames(c("query"))

pos_sel <- pos_sel %>% 
  filter(query %in% genes_godes$query)

#ORA analysis
ora_pos_sel <- enricher(gene = pos_sel$query, # your candidate genes
                         TERM2GENE = genes_godes[,c("Go_desc", "query")], # All genes that you use in the analysis.
                         pAdjustMethod = "none",
                         pvalueCutoff = 0.1)

# Different plots
barplot(ora_pos_sel)
dotplot(ora_pos_sel)
upsetplot(ora_pos_sel)
cnetplot(ora_pos_sel)

# Save results as dataframe
ora_pos_sel_table <- as.data.frame(ora_pos_sel)
head(ora_pos_sel_table)

#----
#Candidate genes selcan
genes_sel_AR <- fread("/path/to/selscan_files/list_genes_sel_AR.txt", header = FALSE)
head(genes_sel_AR)

genes_sel_AR <- as.data.frame(genes_sel_AR) %>%
  setNames(c("query"))

genes_sel_AR <- genes_sel_AR %>% 
  filter(query %in% genes_godes$query)

#ORA analysis
ora_sel_AR <- enricher(gene = genes_sel_AR$query, # your candidate genes
                        TERM2GENE = genes_godes[,c("Go_desc", "query")], # All genes that you use in the analysis.
                        pAdjustMethod = "fdr",
                        pvalueCutoff = 0.5)

# Different plots
barplot(ora_sel_AR)
dotplot(ora_sel_AR)
upsetplot(ora_sel_AR)
cnetplot(ora_sel_AR)

# Save results as dataframe
ora_sel_AR_table <- as.data.frame(ora_sel_AR)
head(ora_pos_sel_table)

#LA pop
genes_sel_LA <- fread("/path/to/selscan_files/list_genes_sel_LA.txt", header = FALSE)
head(genes_sel_LA)

genes_sel_LA <- as.data.frame(genes_sel_LA) %>%
  setNames(c("query"))

genes_sel_LA <- genes_sel_LA %>% 
  filter(query %in% genes_godes$query)

#ORA analysis
ora_sel_LA <- enricher(gene = genes_sel_LA$query, # your candidate genes
                       TERM2GENE = genes_godes[,c("Go_desc", "query")], # All genes that you use in the analysis.
                       pAdjustMethod = "fdr",
                       pvalueCutoff = 0.1)

# Different plots
barplot(ora_sel_LA)
dotplot(ora_sel_LA)
upsetplot(ora_sel_LA)
cnetplot(ora_sel_LA)

# Save results as dataframe
ora_sel_LA_table <- as.data.frame(ora_sel_LA)
head(ora_sel_LA_table)
write.table(ora_sel_LA_table, "/path/to/outdir/ora_sel_LA.txt", sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
