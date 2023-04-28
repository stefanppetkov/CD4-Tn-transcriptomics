library(pheatmap)
library(RColorBrewer)
library(tidyverse)
library(org.Hs.eg.db)
library(hrbrthemes)
library(ggthemes)
library(extrafont)
library(clusterProfiler)
library(clusterProfiler.dplyr)
library(RColorBrewer)

### Heatmap of DEGs (adj. p < 0.01 and log2FC > 1.5

merged_data <- merge(N_EA_v_C_uniprot, N_LA_v_C_uniprot, by = 'row')
merged_data <- merge(merged_data, N_EA_v_LA_uniprot, by = 'row')

gseaRes_ <- rbind(gseaRes_EAvC, gseaRes_LAvC, gseaRes_EAvLA)
cells <- 'NTC'
cutoff <- 0.01


leading_edge <- gseaRes_ %>% 
  dplyr::filter(padj<0.01) %>%
  unnest(leadingEdge) %>%
  pull(leadingEdge) %>%  
  unique()


my_genes <- merged_data %>% 
  dplyr::filter(padj.x <= cutoff | padj.y <= cutoff | padj <= cutoff,
                abs(log2FoldChange.x)>1.5 | abs(log2FoldChange.y)>1.5 | abs(log2FoldChange)>1.5) %>% # stricter threshold
  pull(row) %>% 
  unique()

my_genes_anno <- merged_data %>% 
  dplyr::filter(padj.x <= cutoff | padj.y <= cutoff | padj <= cutoff,
                abs(log2FoldChange.x)>1.5 | abs(log2FoldChange.y)>1.5 | abs(log2FoldChange)>1.5) %>% # stricter threshold
  dplyr::select(row, SYMBOL.x)

# remove duplicated rows
my_genes_anno <- my_genes_anno[!duplicated(my_genes_anno$row), ]
my_genes_anno <- my_genes_anno %>% as_tibble %>% column_to_rownames(var = 'row')

vsd <- vst(dds)
ids <- metaData %>% dplyr::filter(cellType %in% cells) %>% 
  rownames_to_column('sample') %>% 
  mutate(condition = factor(condition, levels = c('Late.Treated', 'Early.Treated', 'Control'))) %>%  
  arrange(condition) 
mat <- assay(vsd)[my_genes, ids$sample]
mat <- mat - rowMeans(mat)

mat_cut_tree = mat
dim(mat)
anno <- metaData %>% dplyr::filter(cellType %in% cells) %>% dplyr::rename(condition, group = condition) %>%  dplyr::select(group) 
anno_colors <- list(group = c(Late.Treated = '#ef9472', Early.Treated = '#9277f1', Control = '#7bec96'))


#row annotation
row_pathways <- gseaRes_ %>% 
  unnest(leadingEdge)


# general search for genes involved in immune response processes
row_anno1 <- row_pathways %>% dplyr::select(pathway, leadingEdge) %>% dplyr::rename(SYMBOL = leadingEdge) %>% dplyr::filter(grepl('immune', pathway)) %>% dplyr::mutate(pathway = 'immune system')

# general search for genes involved in inflammation
row_anno2 <- row_pathways %>% dplyr::select(pathway, leadingEdge) %>% dplyr::rename(SYMBOL = leadingEdge) %>% dplyr::filter(grepl('inflam', pathway)) %>% dplyr::mutate(pathway = 'Inflammation')

# general search for apoptosis genes
row_anno3 <- row_pathways %>% dplyr::select(pathway, leadingEdge) %>% dplyr::rename(SYMBOL = leadingEdge) %>% dplyr::filter(grepl('apop', pathway)) %>% dplyr::mutate(pathway = 'apoptotic signaling pathway')


#convert symbols to ensembl
row_anno1['row'] <- AnnotationDbi::mapIds(org.Hs.eg.db,
                                    keys = row_anno1$SYMBOL, 
                                    keytype ="SYMBOL",
                                    column = "ENSEMBL")
row_anno2['row'] <- AnnotationDbi::mapIds(org.Hs.eg.db,
                                         keys = row_anno2$SYMBOL, 
                                         keytype ="SYMBOL",
                                         column = "ENSEMBL")
row_anno3['row'] <- AnnotationDbi::mapIds(org.Hs.eg.db,
                                         keys = row_anno3$SYMBOL, 
                                         keytype ="SYMBOL",
                                         column = "ENSEMBL")


row_anno1 <- row_anno1[!duplicated(row_anno1$row),]
row_anno2 <- row_anno2[!duplicated(row_anno2$row),]
row_anno3 <- row_anno3[!duplicated(row_anno3$row),]
row_anno1 <- row_anno1 %>% dplyr::select(row, pathway)
row_anno2 <- row_anno2 %>% dplyr::select(row, pathway)
row_anno3 <- row_anno3 %>% dplyr::select(row, pathway)


row_anno <- merge(row_anno1, row_anno2, by = 'row', all = T) %>% merge(row_anno3, by ='row', all = T) %>% column_to_rownames('row')

row_anno_symbs <- row_anno %>% dplyr::mutate(symbols = AnnotationDbi::mapIds(org.Hs.eg.db,
                                                                             keys = rownames(row_anno), 
                                                                             keytype ="ENSEMBL",
                                                                             column = "SYMBOL"))


#### OLO method clustering
cl_cb_new <- function(hcl, mat){
  # Recalculate manhattan distances for reorder method
  dists <- dist(mat, method = "manhattan")

  # Perform reordering according to OLO method
  hclust_olo <- reorder(hcl, dists)
  return(hclust_olo)
}

breaksList = seq(-2, 2, by = 0.01)

### Rearrange mat colums with kmeans
  # cluster each annotation segment
  kmeans_la <- kmeans(t(scale(mat[,grep('K', colnames(mat))])), centers = 3, nstart = 20)$cluster %>% enframe() %>% arrange(value, name) %>% pull(name)
  kmeans_ea <- kmeans(t(scale(mat[,grep('EA', colnames(mat))])), centers = 3, nstart = 20)$cluster %>% enframe() %>% arrange(value, name) %>% pull(name)
  kmeans_c <- kmeans(t(scale(mat[,grep('K|EA', colnames(mat), invert = T)])), centers = 3, nstart = 20)$cluster %>% enframe() %>% arrange(value, name) %>% pull(name)
  
  # order of mat columns
  kmeans_mat <- c(kmeans_la, kmeans_ea, kmeans_c)
  

# number of genes
  dim(mat)
# [1] 1134  39


heatmap <- pheatmap(mat[, kmeans_mat],
         clustering_distance_cols = "manhattan",
         clustering_callback = cl_cb_new,
         cutree_rows = 4,
         cluster_cols = F,
         breaks = breaksList,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),
         annotation_col = anno,
         annotation_row = row_anno,
         show_colnames = F,
         show_rownames = T,
         annotation_colors = anno_colors[1],
         clustering_method = "ward.D2",
         cluster_rows = T,
         labels_row = my_genes_anno$SYMBOL,
         fontsize_row = 0.5,
         cellheight = 0.8,
         treeheight_row = 100,
         treeheight_col = 0,
         width = 6,
         border_color = NA
)
dev.off()
# annotated genes in the matrix
length(intersect(rownames(mat_cut_tree), rownames(row_anno)))
# [1] 211
