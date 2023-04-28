library(tidyverse)
library(org.Hs.eg.db)
library(hrbrthemes)
library(ggthemes)
library(extrafont)
library(clusterProfiler)
library(clusterProfiler.dplyr)
library(RColorBrewer)

### Analysis of annotated genes (apoptosis, immune response and inflammation)

# immune response annotated genes
length(intersect(rownames(mat), row_anno1$row))
immune_genes <- intersect(rownames(mat_cut_tree), row_anno1$row)
# [1] 153

N_EA_v_C_uniprot %>% dplyr::filter(row %in% immune_genes, padj <= cutoff, abs(log2FoldChange) > 1.5) # 67
N_LA_v_C_uniprot %>% dplyr::filter(row %in% immune_genes, padj <= cutoff, abs(log2FoldChange) > 1.5) # 109
N_EA_v_LA_uniprot %>% dplyr::filter(row %in% immune_genes, padj <= cutoff, abs(log2FoldChange) > 1.5) # 21

# inflammation annotated
length(intersect(rownames(mat_cut_tree), row_anno2$row))
# [1] 56

inflam_genes <- intersect(rownames(mat_cut_tree), row_anno2$row)
N_EA_v_C_uniprot %>% dplyr::filter(row %in% inflam_genes, padj <= cutoff, abs(log2FoldChange) > 1.5) # 24
N_LA_v_C_uniprot %>% dplyr::filter(row %in% inflam_genes, padj <= cutoff, abs(log2FoldChange) > 1.5) # 37
N_EA_v_LA_uniprot %>% dplyr::filter(row %in% inflam_genes, padj <= cutoff, abs(log2FoldChange) > 1.5) # 13

# apoptosis annotates
length(intersect(rownames(mat_cut_tree), row_anno3$row))
# [1] 71

apop_genes <- intersect(rownames(mat_cut_tree), row_anno3$row)
N_EA_v_C_uniprot %>% dplyr::filter(row %in% apop_genes, padj <= cutoff, abs(log2FoldChange) > 1.5) # 36
N_LA_v_C_uniprot %>% dplyr::filter(row %in% apop_genes, padj <= cutoff, abs(log2FoldChange) > 1.5) # 52
N_EA_v_LA_uniprot %>% dplyr::filter(row %in% apop_genes, padj <= cutoff, abs(log2FoldChange) > 1.5) # 7

## Exract clusters form heatmap and annotate with symbols, change 'k' to set number of clusters
clusters <- cutree(heatmap$tree_row, k = 4) %>% tibble::enframe(name = 'ENS', value = 'cluster') %>%
  dplyr::mutate(SYMBOL = mapIds(org.Hs.eg.db,
                                keys = ENS, 
                                keytype ="ENSEMBL",
                                column = "SYMBOL")) %>% 
  dplyr::mutate(ENTREZ = mapIds(org.Hs.eg.db,
                                keys = ENS, 
                                keytype ="ENSEMBL",
                                column = "ENTREZID"))

## take gene counts for each group from dds object
normalized_counts <- counts(dds, normalized=TRUE) %>% 
  data.frame(check.names = F) %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  dplyr::select(contains('N')) %>% 
  tidyr::gather('sample', 'counts', -gene) %>% 
  dplyr::mutate(
    group = case_when(
      grepl('EA', sample) ~ 'EA',
      grepl('K', sample) ~ 'LA',
      grepl('^\\d\\d\\d_\\d', sample) ~ 'C')) %>% 
  dplyr::filter(gene %in% clusters$ENS) %>% full_join(clusters, by = c('gene' = 'ENS')) %>% 
  dplyr::mutate(counts = round(counts, 2))



## compare counts of annotated genes in clusters
annotated_genes_entrez <- mapIds(org.Hs.eg.db,
                                 keys = rownames(row_anno), 
                                 keytype ="ENSEMBL",
                                 column = "ENTREZID") %>% unname()
## plot
normalized_counts %>% group_by(group) %>% 
  dplyr::mutate(group = factor(group, levels = c('EA', 'LA', 'C')),
                counts = log10(counts+0.1)) %>%
  dplyr::filter(counts>0, ENTREZ %in% annotated_genes_entrez) %>% 
  ggplot(aes(x= group, y = counts, fill = group)) +
  # geom_jitter(position=position_jitter(0.2)) +
  geom_violin(trim = F, adjust = 1) +
  geom_boxplot(outlier.shape = NA, width = 0.1, fill = 'white') +
  labs(y = 'Log Normalized Counts', x= '') +
  ggthemes::scale_fill_tableau(name=NULL) +
  facet_wrap(~factor(cluster, levels = c(4,2,3,1)), nrow = 1) +
  theme_ipsum() +
  theme(legend.position = 'none') +
  geom_signif(comparisons = list(c("C", "EA"),
                                 c("C", "LA"),
                                 c("EA", "LA")), test = 't.test', color = 'black', 
              map_signif_level=TRUE, step_increase = 0.2)

########## Doplot profiles ##########
# extract annotated genes from culusters
c_list <- list(X1 = clusters[clusters$cluster == 4,4] %>% dplyr::filter(ENTREZ %in% annotated_genes_entrez) %>% pull(ENTREZ) %>% unname(),
               X2 = clusters[clusters$cluster == 2,4] %>% dplyr::filter(ENTREZ %in% annotated_genes_entrez) %>% pull(ENTREZ) %>% unname(),
               X3 = clusters[clusters$cluster == 3,4] %>% dplyr::filter(ENTREZ %in% annotated_genes_entrez) %>% pull(ENTREZ) %>% unname(),
               X4 = clusters[clusters$cluster == 1,4] %>% dplyr::filter(ENTREZ %in% annotated_genes_entrez) %>% pull(ENTREZ) %>% unname()
               )


#### GO enrichment analysis
ck <- compareCluster(geneCluster = c_list, fun = "enrichGO", 
                     OrgDb = org.Hs.eg.db, 
                     keyType = 'ENTREZID', 
                     readable = T,
                     ont = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05)

# generate dotplot
dotplot(clusterProfiler.dplyr::filter(ck, parse_ratio(GeneRatio) > 0.1))
# Save 7x7
