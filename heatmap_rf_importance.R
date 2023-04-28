library(cowplot)
library(tidyverse)
library(org.Hs.eg.db)
library(hrbrthemes)
library(ggthemes)
library(extrafont)
library(seriation)
library(Boruta)
library(RColorBrewer)

### Heatmap of apoptosis-related genes

## helpter functions
# Cluster purity function
# To compute purity , each cluster is assigned to the class which is most frequent in the cluster, and then the accuracy of this assignment is measured by counting the number of correctly assigned cases and dividing by number of total cases.

cluster.purity <- function(clusters, classes) {

  sum(apply(table(classes, clusters), 2, max)) / length(clusters)
}

hc = cutree(p1$tree_col, 2)
cls = anno %>% rownames_to_column('sample') %>% dplyr::filter(sample %in% p1$tree_col$labels) %>% dplyr::mutate(group = as.numeric(group)) %>% deframe()


cluster.purity(hc, cls)

#### OLO clustering
cl_cb_new <- function(hcl, new){
  # Recalculate manhattan distances for reorder method
  dists <- dist(mat, method = "manhattan")
  
  # Perform reordering according to OLO method
  hclust_olo <- reorder(hcl, dists)
  return(hclust_olo)
}


# EA v C 
# Apoptotic pathway heatmap
data <- N_EA_v_C_uniprot
gseaRes <- gseaRes_EAvC
cells <- 'NTC'
comparison <- c('Early.Treated', 'Control')
selection <- gseaRes %>% unnest(leadingEdge) %>% dplyr::select(pathway, leadingEdge) %>% dplyr::rename(SYMBOL = leadingEdge) %>% dplyr::filter(grepl('apop', pathway)) %>% pull(pathway) %>% unique()
cutoff <- 0.01


leading_edge <- gseaRes %>% 
  dplyr::filter(pathway %in% selection) %>% 
  unnest(leadingEdge) %>%
  pull(leadingEdge) %>% 
  unique()

apop_genes <- intersect(rownames(mat_cut_tree), row_anno3$row)

my_genes <- data %>% dplyr::filter(row %in% apop_genes) %>% 
  dplyr::filter(padj <= cutoff, abs(log2FoldChange) > 1.5) %>% 
  pull(row) %>% 
  unique()
my_genes_anno <- data %>% dplyr::filter(row %in% apop_genes) %>% 
  dplyr::filter(padj <= cutoff, abs(log2FoldChange) > 1.5) %>% 
  dplyr::select(row, SYMBOL)

# remove duplicated rows
my_genes_anno <- my_genes_anno[!duplicated(my_genes_anno$row), ]
my_genes_anno <- my_genes_anno %>% as_tibble %>% column_to_rownames(var = 'row')

vsd <- vst(dds)
ids <- metaData %>% dplyr::filter(cellType %in% cells, condition %in% comparison)

mat <- assay(vsd)[my_genes, rownames(ids)]
mat <- mat - rowMeans(mat)
dim(mat)

anno <- metaData %>% dplyr::filter(cellType %in% cells) %>% dplyr::rename(condition, group = condition) %>%  dplyr::select(group) 
anno_colors <- list(group = c(Early.Treated = '#9277f1', Control = '#7bec96'))

#row annotation
row_pathways <- gseaRes %>% 
  dplyr::filter(pathway %in% selection) %>% 
  dplyr::filter(padj < 0.01) %>% 
  unnest(leadingEdge)

row_pathways <- row_pathways %>% dplyr::arrange(leadingEdge, desc(size)) # desc, keeps the largest pathways

row_pathways <- row_pathways[!duplicated(row_pathways$leadingEdge),]
row_anno <- row_pathways %>% dplyr::select(pathway, leadingEdge) %>% dplyr::rename(SYMBOL = leadingEdge)
row_anno <- merge(row_anno, data, by = 'SYMBOL')
row_anno <- row_anno[!duplicated(row_anno$row),]
row_anno <- row_anno %>% dplyr::filter(row %in% my_genes) %>%  dplyr::select(row, pathway) %>% column_to_rownames('row')


breaksList = seq(-2, 2, by = 0.01)


p1 <- pheatmap(mat,
               breaks = breaksList,
               color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),
               clustering_distance_cols = "manhattan",
               clustering_callback = cl_cb_new,
               annotation_col = anno,
               show_colnames = F,
               annotation_colors = anno_colors[1],
               clustering_method = "ward.D",
               cluster_rows = T,
               labels_row = my_genes_anno$SYMBOL,
               fontsize_row = fsize(),
               fontsize_col = 3,
               cellwidth = fsize(),
               cellheight = fsize(),
               main = paste0('EA vs. C,', ' n = ', dim(mat)[1], ', cp = ',  round(cluster.purity(hc, cls), 2)),
               treeheight_row = 0, 
               treeheight_col = 50,
               width = 6.6,
               height = 8,
               border_color = NA
)
hc = cutree(p1$tree_col, 2)
cls = anno %>% rownames_to_column('sample') %>% dplyr::filter(sample %in% p1$tree_col$labels) %>% dplyr::mutate(group = as.numeric(group)) %>% deframe()
p2 = pheatmap(mat,
              breaks = breaksList,
              color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),
              clustering_distance_cols = "manhattan",
              clustering_callback = cl_cb_new,
              annotation_col = anno,
              show_colnames = F,
              annotation_colors = anno_colors[1],
              clustering_method = "ward.D",
              cluster_rows = T,
              labels_row = my_genes_anno$SYMBOL,
              fontsize_col = 3,
              main = paste0('EA vs. C,', ' n = ', dim(mat)[1], ', cp = ',  round(cluster.purity(hc, cls), 2)),
              treeheight_row = 0, 
              treeheight_col = 50,
              width = 6,
              height = 6,
              border_color = NA
)
dev.off()
################## Importance
set.seed(123)
data_ <- data_my[, c(my_genes_anno$SYMBOL, 'group')] %>% 
  dplyr::filter(group %in% c('EA', 'C')) #!!!! change appropriately
boruta_output <- Boruta(x = data_[, my_genes_anno$SYMBOL != "group"], # protect() fix
                        y = data_$group, doTrace=0)
names(boruta_output)

# 3. Get significant variables including tentatives
boruta_signif <- getSelectedAttributes(boruta_output, withTentative = TRUE)
print(boruta_signif)  

# If you are not sure about the tentative variables being selected for granted, you can choose a TentativeRoughFix on boruta_output.
# Do a tentative rough fix
roughFixMod <- TentativeRoughFix(boruta_output)
boruta_signif <- getSelectedAttributes(roughFixMod)
print(boruta_signif)

# Variable Importance Scores
imps <- attStats(roughFixMod)
imps2 <-  imps[imps$decision != 'Rejected', c('meanImp', 'decision')]
imps2 <-  imps2[!row.names(imps2) %in% 'group', ]
head(imps2[order(-imps2$meanImp), ], 100)  # descending sort

# EA v C important genes
apop_ImpList_EAvC <- imps2 %>% dplyr::filter(meanImp >= 3) %>%  rownames_to_column() %>% pull(rowname)


# theme_set(theme_bw(base_size = 12))
imps2 %>% rownames_to_column() %>% dplyr::filter(meanImp >= 3) %>% 
  arrange(meanImp) %>% 
  mutate(name= factor(rowname, levels= rowname)) %>%
  ggplot(aes(y=meanImp, x= name)) + #use d_top30 for top 30 proteins; size is 4x6
  geom_segment(aes(y = 0, x = name, yend = meanImp, xend = name)) + 
  geom_point(stat='identity', color="orange", size=2) +
  labs(x = '', y = 'Importance') +
  ggtitle('Feature importance') +
  theme_ipsum() +
  ylim(0,10) +
  theme(panel.grid = element_blank()) +
  coord_flip()
# ggsave(plot = last_plot(), device = 'pdf', 'APOP_importance_EA.pdf', width = 2.5, height = 3, scale = 1.2)

N_EA_v_C_uniprot %>% dplyr::filter(SYMBOL %in% apop_ImpList_EAvC) %>% view()


##################################################################################### EA v LA #####################################################################################
data <- N_EA_v_LA_uniprot
gseaRes <- gseaRes_EAvLA
cells <- 'NTC'
comparison <- c('Early.Treated', 'Late.Treated')
selection <- gseaRes %>% unnest(leadingEdge) %>% dplyr::select(pathway, leadingEdge) %>% dplyr::rename(SYMBOL = leadingEdge) %>% dplyr::filter(grepl('apop', pathway)) %>% pull(pathway) %>% unique()
cutoff <- 0.01


leading_edge <- gseaRes %>% 
  dplyr::filter(pathway %in% selection) %>% 
  unnest(leadingEdge) %>%
  pull(leadingEdge) %>% 
  unique()


apop_genes <- intersect(rownames(mat_cut_tree), row_anno3$row)

my_genes <- data %>% dplyr::filter(row %in% apop_genes) %>% 
  dplyr::filter(padj <= cutoff, abs(log2FoldChange) > 1.5) %>% 
  pull(row) %>% 
  unique()
my_genes_anno <- data %>% dplyr::filter(row %in% apop_genes) %>% 
  dplyr::filter(padj <= cutoff, abs(log2FoldChange) > 1.5) %>% 
  dplyr::select(row, SYMBOL)

# remove duplicated rows
my_genes_anno <- my_genes_anno[!duplicated(my_genes_anno$row), ]
my_genes_anno <- my_genes_anno %>% as_tibble %>% column_to_rownames(var = 'row')

vsd <- vst(dds)
ids <- metaData %>% dplyr::filter(cellType %in% cells, condition %in% comparison)
mat <- assay(vsd)[my_genes, rownames(ids)]
mat <- mat - rowMeans(mat)
dim(mat)

anno <- metaData %>% dplyr::filter(cellType %in% cells) %>% dplyr::rename(condition, group = condition) %>%  dplyr::select(group) 
anno_colors <- list(group = c(Early.Treated = '#9277f1', Late.Treated = '#ef9472'))

#row annotation
row_pathways <- gseaRes %>% 
  dplyr::filter(pathway %in% selection) %>% 
  dplyr::filter(padj < 0.01) %>% 
  unnest(leadingEdge)

row_pathways <- row_pathways %>% dplyr::arrange(leadingEdge, desc(size)) # desc, keeps the largest pathways

row_pathways <- row_pathways[!duplicated(row_pathways$leadingEdge),]
row_anno <- row_pathways %>% dplyr::select(pathway, leadingEdge) %>% dplyr::rename(SYMBOL = leadingEdge)
row_anno <- merge(row_anno, data, by = 'SYMBOL')
row_anno <- row_anno[!duplicated(row_anno$row),]
row_anno <- row_anno %>% dplyr::filter(row %in% my_genes) %>%  dplyr::select(row, pathway) %>% column_to_rownames('row')



breaksList = seq(-2, 2, by = 0.01)


p1 <- pheatmap(mat,
               breaks = breaksList,
               color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),
               clustering_distance_cols = "manhattan",
               clustering_callback = cl_cb_new,
               annotation_col = anno,
               show_colnames = F,
               annotation_colors = anno_colors[1],
               clustering_method = "ward.D",
               cluster_rows = T,
               labels_row = my_genes_anno$SYMBOL,
               fontsize_col = 3,
               cellwidth = fsize(),
               cellheight = fsize(),
               main = paste0('EA vs. LA,', ' n = ', dim(mat)[1], ', cp = ',  round(cluster.purity(hc, cls), 2)),
               treeheight_row = 0, 
               treeheight_col = 50,
               width = 6,
               height = 8,
               border_color = NA
)
dev.off()
hc = cutree(p1$tree_col, 2)
cls = anno %>% rownames_to_column('sample') %>% dplyr::filter(sample %in% p1$tree_col$labels) %>% dplyr::mutate(group = as.numeric(group)) %>% deframe()
p2 = pheatmap(mat,
              breaks = breaksList,
              color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),
              clustering_distance_cols = "manhattan",
              clustering_callback = cl_cb_new,
              annotation_col = anno,
              show_colnames = F,
              annotation_colors = anno_colors[1],
              clustering_method = "ward.D",
              cluster_rows = T,
              labels_row = my_genes_anno$SYMBOL,
              fontsize_col = 3,
              main = paste0('EA vs. LA,', ' n = ', dim(mat)[1], ', cp = ',  round(cluster.purity(hc, cls), 2)),
              treeheight_row = 0, 
              treeheight_col = 50,
              width = 6,
              height = 6,
              border_color = NA
)
dev.off()

################## Importance
set.seed(123)
data_ <- data_my[, c(my_genes_anno$SYMBOL, 'group')] %>% 
  dplyr::filter(group %in% c('EA', 'LA')) 
boruta_output <- Boruta(x = data_[, my_genes_anno$SYMBOL != "group"], # protect() fix
                        y = data_$group, doTrace=0)
names(boruta_output)

# 3. Get significant variables including tentatives
boruta_signif <- getSelectedAttributes(boruta_output, withTentative = TRUE)
print(boruta_signif)  

# If you are not sure about the tentative variables being selected for granted, you can choose a TentativeRoughFix on boruta_output.
# Do a tentative rough fix
roughFixMod <- TentativeRoughFix(boruta_output)
boruta_signif <- getSelectedAttributes(roughFixMod)
print(boruta_signif)

# Variable Importance Scores
imps <- attStats(roughFixMod)
imps2 <-  imps[imps$decision != 'Rejected', c('meanImp', 'decision')]
imps2 <-  imps2[!row.names(imps2) %in% 'group', ]
head(imps2[order(-imps2$meanImp), ], 100)  # descending sort

# EA v LA important genes
apop_ImpList_EAvLA <- imps2 %>% dplyr::filter(meanImp >= 3) %>%  rownames_to_column() %>% pull(rowname)


# theme_set(theme_bw(base_size = 12))
imps2 %>% rownames_to_column() %>% dplyr::filter(meanImp >= 3) %>% 
  arrange(meanImp) %>%
  mutate(name= factor(rowname, levels= rowname)) %>%
  ggplot(aes(y=meanImp, x= name)) + #use d_top30 for top 30 proteins; size is 4x6
  geom_segment(aes(y = 0, x = name, yend = meanImp, xend = name)) + 
  geom_point(stat='identity', color="orange", size=2) +
  labs(x = '', y = 'Importance') +
  ggtitle('Feature importance') +
  theme_ipsum() +
  ylim(0,10) +
  theme(panel.grid = element_blank()) +
  coord_flip()
# ggsave(plot = last_plot(), device = 'pdf', 'APOP_importance_EAvLA.pdf', width = 2.5, height = 3, scale = 1.2)

N_EA_v_LA_uniprot %>% dplyr::filter(SYMBOL %in% apop_ImpList_EAvLA) %>% view()



########################################################################################### LA v C ###########################################################################################

data <- N_LA_v_C_uniprot
gseaRes <- gseaRes_LAvC 
cells <- 'NTC'
comparison <- c('Late.Treated', 'Control')
selection <- gseaRes %>% unnest(leadingEdge) %>% dplyr::select(pathway, leadingEdge) %>% dplyr::rename(SYMBOL = leadingEdge) %>% dplyr::filter(grepl('apop', pathway)) %>% pull(pathway) %>% unique()
cutoff <- 0.01


leading_edge <- gseaRes %>% 
  dplyr::filter(pathway %in% selection) %>% 
  unnest(leadingEdge) %>%
  pull(leadingEdge) %>% 
  unique()

# have to generate cut_tree heatmap for t e matrix
apop_genes <- intersect(rownames(mat_cut_tree), row_anno3$row)

my_genes <- data %>% dplyr::filter(row %in% apop_genes) %>% 
  dplyr::filter(padj <= cutoff, abs(log2FoldChange) > 1.5) %>% 
  pull(row) %>% 
  unique()
my_genes_anno <- data %>% dplyr::filter(row %in% apop_genes) %>% 
  dplyr::filter(padj <= cutoff, abs(log2FoldChange) > 1.5) %>% 
  dplyr::select(row, SYMBOL)

# remove duplicated rows
my_genes_anno <- my_genes_anno[!duplicated(my_genes_anno$row), ]
my_genes_anno <- my_genes_anno %>% as_tibble %>% column_to_rownames(var = 'row')

vsd <- vst(dds)
ids <- metaData %>% dplyr::filter(cellType %in% cells, condition %in% comparison)
mat <- assay(vsd)[my_genes, rownames(ids)]
mat <- mat - rowMeans(mat)
dim(mat)

anno <- metaData %>% dplyr::filter(cellType %in% cells) %>% dplyr::rename(condition, group = condition) %>%  dplyr::select(group) 
anno_colors <- list(group = c(Late.Treated = '#ef9472', Control = '#7bec96'))

#row annotation
row_pathways <- gseaRes %>% 
  dplyr::filter(pathway %in% selection) %>% 
  dplyr::filter(padj < 0.01) %>% 
  unnest(leadingEdge)
row_pathways <- row_pathways %>% dplyr::arrange(leadingEdge, desc(size)) # desc, keeps the largest pathways

row_pathways <- row_pathways[!duplicated(row_pathways$leadingEdge),]
row_anno <- row_pathways %>% dplyr::select(pathway, leadingEdge) %>% dplyr::rename(SYMBOL = leadingEdge)
row_anno <- merge(row_anno, data, by = 'SYMBOL')
row_anno <- row_anno[!duplicated(row_anno$row),]
row_anno <- row_anno %>% dplyr::filter(row %in% my_genes) %>%  dplyr::select(row, pathway) %>% column_to_rownames('row')



# dev.off()

breaksList = seq(-2, 2, by = 0.01)


p1 <- pheatmap(mat,
               breaks = breaksList,
               color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),
               clustering_distance_cols = "manhattan",
               clustering_callback = cl_cb_new,
               annotation_col = anno,
               show_colnames = F,
               annotation_colors = anno_colors[1],
               clustering_method = "ward.D",
               cluster_rows = T,
               labels_row = my_genes_anno$SYMBOL,
               fontsize_row = fsize(),
               fontsize_col = 3,
               cellwidth = fsize(),
               cellheight = fsize(),
               main = paste0('LA vs. C,', ' n = ', dim(mat)[1], ', cp = ',  round(cluster.purity(hc, cls), 2)),
               treeheight_row = 0, 
               treeheight_col = 50,
               border_color = NA,
               width = 6,
               height = 8
)
dev.off()

hc = cutree(p1$tree_col, 2)
cls = anno %>% rownames_to_column('sample') %>% dplyr::filter(sample %in% p1$tree_col$labels) %>% dplyr::mutate(group = as.numeric(group)) %>% deframe()
p2 = pheatmap(mat,
              breaks = breaksList,
              color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),
              clustering_distance_cols = "manhattan",
              clustering_callback = cl_cb_new,
              annotation_col = anno,
              show_colnames = F,
              annotation_colors = anno_colors[1],
              clustering_method = "ward.D",
              cluster_rows = T,
              labels_row = my_genes_anno$SYMBOL,
              fontsize_col = 3,
              main = paste0('LA vs. C,', ' n = ', dim(mat)[1], ', cp = ',  round(cluster.purity(hc, cls), 2)),
              treeheight_row = 0, 
              treeheight_col = 50,
              width = 6,
              height = 6,
              border_color = NA
)
dev.off()
################## Importance
set.seed(123)
data_ <- data_my[, c(my_genes_anno$SYMBOL, 'group')] %>% 
  dplyr::filter(group %in% c('LA', 'C'))
boruta_output <- Boruta(x = data_[, my_genes_anno$SYMBOL != "group"], # protect() fix
                        y = data_$group, doTrace=0)
names(boruta_output)

# 3. Get significant variables including tentatives
boruta_signif <- getSelectedAttributes(boruta_output, withTentative = TRUE)
print(boruta_signif)  

# If you are not sure about the tentative variables being selected for granted, you can choose a TentativeRoughFix on boruta_output.
# Do a tentative rough fix
roughFixMod <- TentativeRoughFix(boruta_output)
boruta_signif <- getSelectedAttributes(roughFixMod)
print(boruta_signif)

# Variable Importance Scores
imps <- attStats(roughFixMod)
imps2 <-  imps[imps$decision != 'Rejected', c('meanImp', 'decision')]
imps2 <-  imps2[!row.names(imps2) %in% 'group', ]
head(imps2[order(-imps2$meanImp), ], 100)  # descending sort

# LA v C important genes
apop_ImpList_LAvC <- imps2 %>% dplyr::filter(meanImp >= 3) %>%  rownames_to_column() %>% pull(rowname)


# theme_set(theme_bw(base_size = 12))
imps2 %>% rownames_to_column() %>% dplyr::filter(meanImp >= 3) %>% 
  arrange(meanImp) %>% 
  mutate(name= factor(rowname, levels= rowname)) %>%
  ggplot(aes(y=meanImp, x= name)) + #use d_top30 for top 30 proteins; size is 4x6
  geom_segment(aes(y = 0, x = name, yend = meanImp, xend = name)) + 
  geom_point(stat='identity', color="orange", size=2) +
  labs(x = '', y = 'Importance') +
  ggtitle('Feature importance') +
  theme_ipsum() +
  ylim(0,10) +
  theme(panel.grid = element_blank()) +
  coord_flip()
# ggsave(plot = last_plot(), device = 'pdf', 'APOP_importance_LA.pdf', width = 2.5, height = 3, scale = 1.2)



### Heatmap of inflammation-related genes

# Inflammatory response heatmap
### helper functions
#font size
fsize <- function(x=length(my_genes)){
  n <- sqrt(x^-2)*350
  if (n>10) {
    n=10
  }
  return(n)
}

# Cluster purity function
cluster.purity <- function(clusters, classes) {
  
  sum(apply(table(classes, clusters), 2, max)) / length(clusters)
}

hc = cutree(p1$tree_col, 2)
cls = anno %>% rownames_to_column('sample') %>% dplyr::filter(sample %in% p1$tree_col$labels) %>% dplyr::mutate(group = as.numeric(group)) %>% deframe()


cluster.purity(hc, cls)

#### for better arrangement of cols
library(seriation)
cl_cb_new <- function(hcl, mat){
  # Recalculate manhattan distances for reorder method
  dists <- dist(mat, method = "manhattan")
  
  # Perform reordering according to OLO method
  hclust_olo <- reorder(hcl, dists)
  return(hclust_olo)
}

##################################################################################### EA v C #####################################################################################
data <- N_EA_v_C_uniprot
gseaRes <- gseaRes_EAvC #use GO_pathays of interest to generate
cells <- 'NTC'
comparison <- c('Early.Treated', 'Control')
selection <- gseaRes %>% unnest(leadingEdge) %>% dplyr::select(pathway, leadingEdge) %>% dplyr::rename(SYMBOL = leadingEdge) %>% dplyr::filter(grepl('inflam', pathway)) %>% pull(pathway) %>% unique()
cutoff <- 0.01


# based on pathways containing "immune"
leading_edge <- gseaRes %>% 
  dplyr::filter(pathway %in% selection) %>% 
  unnest(leadingEdge) %>%
  pull(leadingEdge) %>% 
  unique()

inflam_genes <- intersect(rownames(mat_cut_tree), row_anno2$row)

my_genes <- data %>% dplyr::filter(row %in% inflam_genes) %>% 
  dplyr::filter(padj <= cutoff, abs(log2FoldChange) > 1.5) %>% 
  pull(row) %>% 
  unique()

my_genes_anno <- data %>% dplyr::filter(row %in% inflam_genes) %>% 
  dplyr::filter(padj <= cutoff, abs(log2FoldChange) > 1.5) %>% 
  dplyr::select(row, SYMBOL)

# remove duplicated rows
my_genes_anno <- my_genes_anno[!duplicated(my_genes_anno$row), ]
my_genes_anno <- my_genes_anno %>% as_tibble %>% column_to_rownames(var = 'row')

vsd <- vst(dds)
ids <- metaData %>% dplyr::filter(cellType %in% cells, condition %in% comparison)
mat <- assay(vsd)[my_genes, rownames(ids)]
mat <- mat - rowMeans(mat)
dim(mat)

anno <- metaData %>% dplyr::filter(cellType %in% cells) %>% dplyr::rename(condition, group = condition) %>%  dplyr::select(group) 
anno_colors <- list(group = c(Early.Treated = '#9277f1', Control = '#7bec96'))

#row annotation
row_pathways <- gseaRes %>% 
  dplyr::filter(pathway %in% selection) %>% # was GO_id
  dplyr::filter(padj < 0.01) %>% 
  unnest(leadingEdge)

row_pathways <- row_pathways %>% dplyr::arrange(leadingEdge, desc(size)) # desc, keeps the largest pathways

row_pathways <- row_pathways[!duplicated(row_pathways$leadingEdge),]
row_anno <- row_pathways %>% dplyr::select(pathway, leadingEdge) %>% dplyr::rename(SYMBOL = leadingEdge)
row_anno <- merge(row_anno, data, by = 'SYMBOL')
row_anno <- row_anno[!duplicated(row_anno$row),]
row_anno <- row_anno %>% dplyr::filter(row %in% my_genes) %>%  dplyr::select(row, pathway) %>% column_to_rownames('row')


breaksList = seq(-2, 2, by = 0.01)

p1 <- pheatmap(mat,
               breaks = breaksList,
               color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),
               clustering_distance_cols = "manhattan",
               clustering_callback = cl_cb_new,
               annotation_col = anno,
               show_colnames = F,
               annotation_colors = anno_colors[1],
               clustering_method = "ward.D",
               cluster_rows = T,
               labels_row = my_genes_anno$SYMBOL,
               fontsize_col = 3,
               main = paste0('EA vs. C,', ' n = ', dim(mat)[1], ', cp = ',  round(cluster.purity(hc, cls), 2)),
               treeheight_row = 0, 
               treeheight_col = 50,
               width = 6,
               height = 6,
               border_color = NA
)
dev.off()
hc = cutree(p1$tree_col, 2)
cls = anno %>% rownames_to_column('sample') %>% dplyr::filter(sample %in% p1$tree_col$labels) %>% dplyr::mutate(group = as.numeric(group)) %>% deframe()
p2 = pheatmap(mat,
              breaks = breaksList,
              color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),
              clustering_distance_cols = "manhattan",
              clustering_callback = cl_cb_new,
              annotation_col = anno,
              show_colnames = F,
              annotation_colors = anno_colors[1],
              clustering_method = "ward.D",
              cluster_rows = T,
              labels_row = my_genes_anno$SYMBOL,
              fontsize_col = 3,
              main = paste0('EA vs. C,', ' n = ', dim(mat)[1], ', cp = ',  round(cluster.purity(hc, cls), 2)),
              treeheight_row = 0, 
              treeheight_col = 50,
              width = 6,
              height = 6,
              border_color = NA
)
dev.off()
################## Importance

set.seed(123)
data_ <- data_my[, c(my_genes_anno$SYMBOL, 'group')] %>% 
  dplyr::filter(group %in% c('EA', 'C')) #!!!! change appropriately
boruta_output <- Boruta(x = data_[, my_genes_anno$SYMBOL != "group"], # protect() solution
                        y = data_$group, doTrace=0)
names(boruta_output)

# 3. Get significant variables including tentatives
boruta_signif <- getSelectedAttributes(boruta_output, withTentative = TRUE)
print(boruta_signif)  

# If you are not sure about the tentative variables being selected for granted, you can choose a TentativeRoughFix on boruta_output.
# Do a tentative rough fix
roughFixMod <- TentativeRoughFix(boruta_output)
boruta_signif <- getSelectedAttributes(roughFixMod)
print(boruta_signif)

# Variable Importance Scores
imps <- attStats(roughFixMod)
imps2 <-  imps[imps$decision != 'Rejected', c('meanImp', 'decision')]
imps2 <-  imps2[!row.names(imps2) %in% 'group', ]
head(imps2[order(-imps2$meanImp), ], 100)  # descending sort

# EA v C important genes
inflam_ImpList_EAvC <- imps2 %>% dplyr::filter(meanImp >= 3) %>%  rownames_to_column() %>% pull(rowname)

library(hrbrthemes)
imps2 %>% rownames_to_column() %>% dplyr::filter(meanImp >= 3) %>% 
  arrange(meanImp) %>% 
  mutate(name= factor(rowname, levels= rowname)) %>%
  ggplot(aes(y=meanImp, x= name)) + #use d_top30 for top 30 proteins; size is 4x6
  geom_segment(aes(y = 0, x = name, yend = meanImp, xend = name)) + 
  geom_point(stat='identity', color="orange", size=2) +
  labs(x = '', y = 'Importance') +
  ggtitle('Feature importance') +
  theme_ipsum() +
  ylim(0, 10) +
  theme(panel.grid = element_blank()) +
  coord_flip()

  # ggsave(plot = last_plot(), device = 'pdf', 'INFLAM_importance_EA.pdf', width = 2.5, height = 3, scale = 1.2)
  

##################################################################################### EA v LA #####################################################################################
# Inflammatory response heatmap
data <- N_EA_v_LA_uniprot
gseaRes <- gseaRes_EAvLA 
cells <- 'NTC'
comparison <- c('Early.Treated', 'Late.Treated')
selection <- 'inflammatory response'
cutoff <- 0.01

leading_edge <- gseaRes %>% 
  dplyr::filter(pathway %in% selection) %>% 
  unnest(leadingEdge) %>%
  pull(leadingEdge) %>%  #pathways within search matching farancesca's list
  unique()

inflam_genes <- intersect(rownames(mat_cut_tree), row_anno2$row)

my_genes <- data %>% dplyr::filter(row %in% inflam_genes) %>% 
  dplyr::filter(padj <= cutoff, abs(log2FoldChange) > 1.5) %>% 
  pull(row) %>% 
  unique()

my_genes_anno <- data %>% dplyr::filter(row %in% inflam_genes) %>% 
  dplyr::filter(padj <= cutoff, abs(log2FoldChange) > 1.5) %>% 
  dplyr::select(row, SYMBOL)


# remove duplicated rows
my_genes_anno <- my_genes_anno[!duplicated(my_genes_anno$row), ]
my_genes_anno <- my_genes_anno %>% as_tibble %>% column_to_rownames(var = 'row')

vsd <- vst(dds)
ids <- metaData %>% dplyr::filter(cellType %in% cells, condition %in% comparison)
mat <- assay(vsd)[my_genes, rownames(ids)]
mat <- mat - rowMeans(mat)
dim(mat)

anno <- metaData %>% dplyr::filter(cellType %in% cells) %>% dplyr::rename(condition, group = condition) %>%  dplyr::select(group) 
anno_colors <- list(group = c(Early.Treated = '#9277f1', Late.Treated = '#ef9472'))

#row annotation
row_pathways <- gseaRes %>% 
  dplyr::filter(pathway %in% selection) %>% 
  unnest(leadingEdge)

row_pathways <- row_pathways %>% dplyr::arrange(leadingEdge, desc(size)) # desc, keeps the largest pathways

row_pathways <- row_pathways[!duplicated(row_pathways$leadingEdge),]
row_anno <- row_pathways %>% dplyr::select(pathway, leadingEdge) %>% dplyr::rename(SYMBOL = leadingEdge)
row_anno <- merge(row_anno, data, by = 'SYMBOL')
row_anno <- row_anno[!duplicated(row_anno$row),]
row_anno <- row_anno %>% dplyr::filter(row %in% my_genes) %>%  dplyr::select(row, pathway) %>% column_to_rownames('row')



breaksList = seq(-2, 2, by = 0.01)


p1 <- pheatmap(mat,
               breaks = breaksList,
               color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),
               clustering_distance_cols = "manhattan",
               clustering_callback = cl_cb_new,
               annotation_col = anno,
               show_colnames = F,
               annotation_colors = anno_colors[1],
               clustering_method = "ward.D",
               cluster_rows = T,
               labels_row = my_genes_anno$SYMBOL,
               fontsize_col = 3,
               cellheight = fsize(),
               main = paste0('EA vs. LA,', ' n = ', dim(mat)[1], ', cp = ',  round(cluster.purity(hc, cls), 2)),
               treeheight_row = 0, 
               treeheight_col = 50,
               width = 6,
               height = 8,
               border_color = NA
)
dev.off()
hc = cutree(p1$tree_col, 2)
cls = anno %>% rownames_to_column('sample') %>% dplyr::filter(sample %in% p1$tree_col$labels) %>% dplyr::mutate(group = as.numeric(group)) %>% deframe()
p2 = pheatmap(mat,
              breaks = breaksList,
              color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),
              clustering_distance_cols = "manhattan",
              clustering_callback = cl_cb_new,
              annotation_col = anno,
              show_colnames = F,
              annotation_colors = anno_colors[1],
              clustering_method = "ward.D",
              cluster_rows = T,
              labels_row = my_genes_anno$SYMBOL,
              fontsize_col = 3,
              main = paste0('EA vs. LA,', ' n = ', dim(mat)[1], ', cp = ',  round(cluster.purity(hc, cls), 2)),
              treeheight_row = 0, 
              treeheight_col = 50,
              width = 6,
              height = 6,
              border_color = NA
)
dev.off()
################## Importance
set.seed(123)

data_ <- data_my[, c(my_genes_anno$SYMBOL, 'group')] %>% 
  dplyr::filter(group %in% c('EA', 'LA')) #!!!! change appropriately
boruta_output <- Boruta(x = data_[, my_genes_anno$SYMBOL != "group"], # protect() solution
                        y = data_$group, doTrace=0)
names(boruta_output)

# 3. Get significant variables including tentatives
boruta_signif <- getSelectedAttributes(boruta_output, withTentative = TRUE)
print(boruta_signif)  

# If you are not sure about the tentative variables being selected for granted, you can choose a TentativeRoughFix on boruta_output.
# Do a tentative rough fix
roughFixMod <- TentativeRoughFix(boruta_output)
boruta_signif <- getSelectedAttributes(roughFixMod)
print(boruta_signif)

# Variable Importance Scores
imps <- attStats(roughFixMod)
imps2 <-  imps[imps$decision != 'Rejected', c('meanImp', 'decision')]
imps2 <-  imps2[!row.names(imps2) %in% 'group', ]
head(imps2[order(-imps2$meanImp), ], 100)  # descending sort

# EA v LA important genes
inflam_ImpList_EAvLA <- imps2 %>% dplyr::filter(meanImp >= 3) %>%  rownames_to_column() %>% pull(rowname)


imps2 %>% rownames_to_column() %>% dplyr::filter(meanImp >= 3) %>% 
  arrange(meanImp) %>% 
  mutate(name= factor(rowname, levels= rowname)) %>%
  ggplot(aes(y=meanImp, x= name)) + #use d_top30 for top 30 proteins; size is 4x6
  geom_segment(aes(y = 0, x = name, yend = meanImp, xend = name)) + 
  geom_point(stat='identity', color="orange", size=2) +
  labs(x = '', y = 'Importance') +
  ggtitle('Feature importance') +
  theme_ipsum() +
  ylim(0,10) +
  theme(panel.grid = element_blank()) +
  coord_flip()
# ggsave(plot = last_plot(), device = 'pdf', 'INFLAM_importance_EAvLA.pdf', width = 2.5, height = 3, scale = 1.2)


########################################################################################### LA v C ###########################################################################################
# Inflammatory response heatmap
data <- N_LA_v_C_uniprot
gseaRes <- gseaRes_LAvC
cells <- 'NTC'
comparison <- c('Late.Treated', 'Control')
selection <- gseaRes %>% unnest(leadingEdge) %>% dplyr::select(pathway, leadingEdge) %>% dplyr::rename(SYMBOL = leadingEdge) %>% dplyr::filter(grepl('inflam', pathway)) %>% pull(pathway) %>% unique()
cutoff <- 0.01


leading_edge <- gseaRes %>% 
  dplyr::filter(pathway %in% selection) %>% 
  unnest(leadingEdge) %>%
  pull(leadingEdge) %>% 
  unique()


inflam_genes <- intersect(rownames(mat_cut_tree), row_anno2$row)


my_genes <- data %>% dplyr::filter(row %in% inflam_genes) %>% 
  dplyr::filter(padj <= cutoff, abs(log2FoldChange) > 1.5) %>% 
  pull(row) %>% 
  unique()

my_genes_anno <- data %>% dplyr::filter(row %in% inflam_genes) %>% 
  dplyr::filter(padj <= cutoff, abs(log2FoldChange) > 1.5) %>% 
  dplyr::select(row, SYMBOL)

# remove duplicated rows
my_genes_anno <- my_genes_anno[!duplicated(my_genes_anno$row), ]
my_genes_anno <- my_genes_anno %>% as_tibble %>% column_to_rownames(var = 'row')

vsd <- vst(dds)
ids <- metaData %>% dplyr::filter(cellType %in% cells, condition %in% comparison)
mat <- assay(vsd)[my_genes, rownames(ids)]
mat <- mat - rowMeans(mat)
dim(mat)

anno <- metaData %>% dplyr::filter(cellType %in% cells) %>% dplyr::rename(condition, group = condition) %>%  dplyr::select(group) 
anno_colors <- list(group = c(Late.Treated = '#ef9472', Control = '#7bec96'))

#row annotation
row_pathways <- gseaRes %>% 
  dplyr::filter(pathway %in% selection) %>% 
  dplyr::filter(padj < 0.01) %>% 
  unnest(leadingEdge)

row_pathways <- row_pathways %>% dplyr::arrange(leadingEdge, desc(size)) # desc, keeps the largest pathways

row_pathways <- row_pathways[!duplicated(row_pathways$leadingEdge),]
row_anno <- row_pathways %>% dplyr::select(pathway, leadingEdge) %>% dplyr::rename(SYMBOL = leadingEdge)
row_anno <- merge(row_anno, data, by = 'SYMBOL')
row_anno <- row_anno[!duplicated(row_anno$row),]
row_anno <- row_anno %>% dplyr::filter(row %in% my_genes) %>%  dplyr::select(row, pathway) %>% column_to_rownames('row')



breaksList = seq(-2, 2, by = 0.01)


p1 <- pheatmap(mat,
         breaks = breaksList,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),
         clustering_distance_cols = "manhattan",
         clustering_callback = cl_cb_new,
         annotation_col = anno,
         show_colnames = F,
         annotation_colors = anno_colors[1],
         clustering_method = "ward.D",
         cluster_rows = T,
         labels_row = my_genes_anno$SYMBOL,
         fontsize_row = fsize(),
         fontsize_col = 3,
         cellwidth = fsize(),
         cellheight = fsize(),
         main = paste0('LA vs. C,', ' n = ', dim(mat)[1], ', cp = ',  round(cluster.purity(hc, cls), 2)),
         treeheight_row = 0, 
         treeheight_col = 50,
         border_color = NA,
         width = 6,
         height = 8
)
dev.off()
hc = cutree(p1$tree_col, 2)
cls = anno %>% rownames_to_column('sample') %>% dplyr::filter(sample %in% p1$tree_col$labels) %>% dplyr::mutate(group = as.numeric(group)) %>% deframe()
p2 <-  pheatmap(mat,
              breaks = breaksList,
              color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),
              clustering_distance_cols = "manhattan",
              clustering_callback = cl_cb_new,
              annotation_col = anno,
              show_colnames = F,
              annotation_colors = anno_colors[1],
              clustering_method = "ward.D",
              cluster_rows = T,
              labels_row = my_genes_anno$SYMBOL,
              fontsize_col = 3,
              main = paste0('LA vs. C,', ' n = ', dim(mat)[1], ', cp = ',  round(cluster.purity(hc, cls), 2)),
              treeheight_row = 0, 
              treeheight_col = 50,
              width = 6,
              height = 6,
              border_color = NA
)
dev.off()
# This creates the annotation of importance ranked genes only
# Must run the boruta first and then come back 
imp_anno <- imps2 %>% rownames_to_column() %>% pull(rowname)
rows <- my_genes_anno %>% dplyr::mutate(SYMBOL = ifelse(SYMBOL %in% imp_anno, SYMBOL, NA))

dev.off()
################## Importance
set.seed(123)
# data_ <- data_my[, c(intersect(my_genes_anno$SYMBOL, names(data_my)), 'group')] # to exclude undefined column

data_ <- data_my[, c(my_genes_anno$SYMBOL, 'group')] %>% 
  dplyr::filter(group %in% c('LA', 'C')) #!!!! change appropriately
boruta_output <- Boruta(x = data_[, my_genes_anno$SYMBOL != "group"], # protect() solution
                        y = data_$group, doTrace=0)
names(boruta_output)

# 3. Get significant variables including tentatives
boruta_signif <- getSelectedAttributes(boruta_output, withTentative = TRUE)
print(boruta_signif)  

# If you are not sure about the tentative variables being selected for granted, you can choose a TentativeRoughFix on boruta_output.
# Do a tentative rough fix
roughFixMod <- TentativeRoughFix(boruta_output)
boruta_signif <- getSelectedAttributes(roughFixMod)
print(boruta_signif)

# Variable Importance Scores
imps <- attStats(roughFixMod)
imps2 <-  imps[imps$decision != 'Rejected', c('meanImp', 'decision')]
imps2 <-  imps2[!row.names(imps2) %in% 'group', ]
head(imps2[order(-imps2$meanImp), ], 100)  # descending sort

# LA v C important genes
inflam_ImpList_LAvC <- imps2 %>% dplyr::filter(meanImp >= 3) %>%  rownames_to_column() %>% pull(rowname)


# theme_set(theme_bw(base_size = 12))
imps2 %>% rownames_to_column() %>% dplyr::filter(meanImp >= 3) %>% 
  arrange(meanImp) %>%
  mutate(name= factor(rowname, levels= rowname)) %>%
  ggplot(aes(y=meanImp, x= name)) + #use d_top30 for top 30 proteins; size is 4x6
  geom_segment(aes(y = 0, x = name, yend = meanImp, xend = name)) + 
  geom_point(stat='identity', color="orange", size=2) +
  labs(x = '', y = 'Importance') +
  ggtitle('Feature importance') +
  theme_ipsum() +
  ylim(0,10) +
  theme(panel.grid = element_blank()) +
  coord_flip()
# ggsave(plot = last_plot(), device = 'pdf', 'INFLAM_importance_LA.pdf', width = 2.5, height = 3, scale = 1.2)




### Heatmap of genes involved in immune response

# Cluster purity function
cluster.purity <- function(clusters, classes) {
  
  sum(apply(table(classes, clusters), 2, max)) / length(clusters)
}

hc = cutree(p1$tree_col, 2)
cls = anno %>% rownames_to_column('sample') %>% dplyr::filter(sample %in% p1$tree_col$labels) %>% dplyr::mutate(group = as.numeric(group)) %>% deframe()


cluster.purity(hc, cls)

#### for better arrangement of cols
cl_cb_new <- function(hcl, mat){
  # Recalculate manhattan distances for reorder method
  dists <- dist(mat, method = "manhattan")
  
  # Perform reordering according to OLO method
  hclust_olo <- reorder(hcl, dists)
  return(hclust_olo)
}

##################################################################################### EA v C #####################################################################################
# Immune response heatmap
data <- N_EA_v_C_uniprot
gseaRes <- gseaRes_EAvC
cells <- 'NTC'
comparison <- c('Early.Treated', 'Control')
selection <- gseaRes %>% unnest(leadingEdge) %>% dplyr::select(pathway, leadingEdge) %>% dplyr::rename(SYMBOL = leadingEdge) %>% dplyr::filter(grepl('immune', pathway)) %>% pull(pathway) %>% unique()
cutoff <- 0.01


leading_edge <- gseaRes %>% 
  dplyr::filter(pathway %in% selection) %>% 
  unnest(leadingEdge) %>%
  pull(leadingEdge) %>% 
  unique()

immune_genes <- intersect(rownames(mat_cut_tree), row_anno1$row)

my_genes <- data %>% dplyr::filter(row %in% immune_genes) %>% 
  dplyr::filter(padj <= cutoff, abs(log2FoldChange) > 1.5) %>% 
  pull(row) %>% 
  unique()

my_genes_anno <- data %>% dplyr::filter(row %in% immune_genes) %>% 
  dplyr::filter(padj <= cutoff, abs(log2FoldChange) > 1.5) %>% 
  dplyr::select(row, SYMBOL)


# remove duplicated rows
my_genes_anno <- my_genes_anno[!duplicated(my_genes_anno$row), ]
my_genes_anno <- my_genes_anno %>% as_tibble %>% column_to_rownames(var = 'row')

vsd <- vst(dds)
ids <- metaData %>% dplyr::filter(cellType %in% cells, condition %in% comparison)
mat <- assay(vsd)[my_genes, rownames(ids)]
mat <- mat - rowMeans(mat)
dim(mat)

anno <- metaData %>% dplyr::filter(cellType %in% cells) %>% dplyr::rename(condition, group = condition) %>%  dplyr::select(group) 
anno_colors <- list(group = c(Early.Treated = '#9277f1', Control = '#7bec96'))

#row annotation
row_pathways <- gseaRes %>% 
  dplyr::filter(pathway %in% selection) %>% 
  dplyr::filter(padj < 0.01) %>% 
  unnest(leadingEdge)


row_pathways <- row_pathways %>% dplyr::arrange(leadingEdge, desc(size)) # desc, keeps the largest pathways

row_pathways <- row_pathways[!duplicated(row_pathways$leadingEdge),]
row_anno <- row_pathways %>% dplyr::select(pathway, leadingEdge) %>% dplyr::rename(SYMBOL = leadingEdge)
row_anno <- merge(row_anno, data, by = 'SYMBOL')
row_anno <- row_anno[!duplicated(row_anno$row),]
row_anno <- row_anno %>% dplyr::filter(row %in% my_genes) %>%  dplyr::select(row, pathway) %>% column_to_rownames('row')


breaksList = seq(-2, 2, by = 0.01)


p1 <- pheatmap(mat,
               breaks = breaksList,
               color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),
               clustering_distance_cols = "manhattan",
               clustering_callback = cl_cb_new,
               annotation_col = anno,
               show_colnames = F,
               annotation_colors = anno_colors[1],
               cluster_rows = T,
               labels_row = my_genes_anno$SYMBOL,
               fontsize_row = fsize(),
               fontsize_col = 3,
               cellwidth = fsize(),
               cellheight = fsize(),
               main = paste0('EA vs. C,', ' n = ', dim(mat)[1], ', cp = ',  round(cluster.purity(hc, cls), 2)),
               treeheight_row = 0, 
               treeheight_col = 50,
               width = 6.6,
               height = 8,
               border_color = NA
)
dev.off()
hc = cutree(p1$tree_col, 2)
cls = anno %>% rownames_to_column('sample') %>% dplyr::filter(sample %in% p1$tree_col$labels) %>% dplyr::mutate(group = as.numeric(group)) %>% deframe()
p2 <- pheatmap(mat,
               filename = 'IMMUNE_EA.pdf',
               breaks = breaksList,
               color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),
               clustering_distance_cols = "manhattan",
               clustering_callback = cl_cb_new,
               annotation_col = anno,
               show_colnames = F,
               annotation_colors = anno_colors[1],
               cluster_rows = T,
               labels_row = my_genes_anno$SYMBOL,
               fontsize_row = fsize(),
               fontsize_col = 3,
               cellwidth = fsize(),
               cellheight = fsize(),
               main = paste0('EA vs. C,', ' n = ', dim(mat)[1], ', cp = ',  round(cluster.purity(hc, cls), 2)),
               treeheight_row = 0, 
               treeheight_col = 50,
               width = 6.6,
               height = 8,
               border_color = NA
)
dev.off()
################## Importance
set.seed(123)
data_ <- data_my[, c(intersect(my_genes_anno$SYMBOL, names(data_my)), 'group')] %>% 
  dplyr::filter(group %in% c('EA', 'C')) #!!!! change appropriately
boruta_output <- Boruta(x = data_[, intersect(my_genes_anno$SYMBOL, names(data_my)) != "group"], # protect() solution
                        y = data_$group, doTrace=0)
names(boruta_output)

# 3. Get significant variables including tentatives
boruta_signif <- getSelectedAttributes(boruta_output, withTentative = TRUE)
print(boruta_signif)  

# If you are not sure about the tentative variables being selected for granted, you can choose a TentativeRoughFix on boruta_output.
# Do a tentative rough fix
roughFixMod <- TentativeRoughFix(boruta_output)
boruta_signif <- getSelectedAttributes(roughFixMod)
print(boruta_signif)

# Variable Importance Scores
imps <- attStats(roughFixMod)
imps2 <-  imps[imps$decision != 'Rejected', c('meanImp', 'decision')]
imps2 <-  imps2[!row.names(imps2) %in% 'group', ]
head(imps2[order(-imps2$meanImp), ], 100)  # descending sort

# EA v C important genes
immune_ImpList_EAvC <- imps2 %>% dplyr::filter(meanImp >= 3) %>%  rownames_to_column() %>% pull(rowname)


# theme_set(theme_bw(base_size = 12))
imps2 %>% rownames_to_column() %>% dplyr::filter(meanImp >= 3) %>% 
  arrange(meanImp) %>% 
  mutate(name= factor(rowname, levels= rowname)) %>%
  ggplot(aes(y=meanImp, x= name)) + 
  geom_segment(aes(y = 0, x = name, yend = meanImp, xend = name)) + 
  geom_point(stat='identity', color="orange", size=2) +
  labs(x = '', y = 'Importance') +
  ggtitle('Feature importance') +
  theme_ipsum() +
  ylim(0,10) +
  theme(panel.grid = element_blank()) +
  coord_flip()
# ggsave(plot = last_plot(), device = 'pdf', 'IMMUNE_importance_EA.pdf', width = 2.5, height = 3, scale = 1.2)

##################################################################################### EA v LA #####################################################################################
data <- N_EA_v_LA_uniprot
gseaRes <- gseaRes_EAvLA 
cells <- 'NTC'
comparison <- c('Early.Treated', 'Late.Treated')
selection <- gseaRes %>% unnest(leadingEdge) %>% dplyr::select(pathway, leadingEdge) %>% dplyr::rename(SYMBOL = leadingEdge) %>% dplyr::filter(grepl('immune', pathway)) %>% pull(pathway) %>% unique()
cutoff <- 0.01


leading_edge <- gseaRes %>% 
  dplyr::filter(pathway %in% selection) %>% 
  unnest(leadingEdge) %>%
  pull(leadingEdge) %>% 
  unique()

immune_genes <- intersect(rownames(mat_cut_tree), row_anno1$row)

my_genes <- data %>% dplyr::filter(row %in% immune_genes) %>% 
  dplyr::filter(padj <= cutoff, abs(log2FoldChange) > 1.5) %>% 
  pull(row) %>% 
  unique()

my_genes_anno <- data %>% dplyr::filter(row %in% immune_genes) %>% 
  dplyr::filter(padj <= cutoff, abs(log2FoldChange) > 1.5) %>% 
  dplyr::select(row, SYMBOL)


# remove duplicated rows
my_genes_anno <- my_genes_anno[!duplicated(my_genes_anno$row), ]
my_genes_anno <- my_genes_anno %>% as_tibble %>% column_to_rownames(var = 'row')

vsd <- vst(dds)
ids <- metaData %>% dplyr::filter(cellType %in% cells, condition %in% comparison)
mat <- assay(vsd)[my_genes, rownames(ids)]
mat <- mat - rowMeans(mat)
dim(mat)

anno <- metaData %>% dplyr::filter(cellType %in% cells) %>% dplyr::rename(condition, group = condition) %>%  dplyr::select(group) 
anno_colors <- list(group = c(Early.Treated = '#9277f1', Late.Treated = '#ef9472'))

#row annotation
row_pathways <- gseaRes %>% 
  dplyr::filter(pathway %in% selection) %>% 
  dplyr::filter(padj < 0.01) %>% 
  unnest(leadingEdge)

row_pathways <- row_pathways %>% dplyr::arrange(leadingEdge, desc(size)) # desc, keeps the largest pathways

row_pathways <- row_pathways[!duplicated(row_pathways$leadingEdge),]
row_anno <- row_pathways %>% dplyr::select(pathway, leadingEdge) %>% dplyr::rename(SYMBOL = leadingEdge)
row_anno <- merge(row_anno, data, by = 'SYMBOL')
row_anno <- row_anno[!duplicated(row_anno$row),]
row_anno <- row_anno %>% dplyr::filter(row %in% my_genes) %>%  dplyr::select(row, pathway) %>% column_to_rownames('row')


breaksList = seq(-2, 2, by = 0.01)


p1 <- pheatmap(mat,
               breaks = breaksList,
               color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),
               clustering_distance_cols = "manhattan",
               clustering_callback = cl_cb_new,
               annotation_col = anno,
               show_colnames = F,
               annotation_colors = anno_colors[1],
               clustering_method = "ward.D",
               cluster_rows = T,
               labels_row = my_genes_anno$SYMBOL,
               fontsize_row = fsize(),
               fontsize_col = 3,
               cellwidth = fsize(),
               cellheight = fsize(),
               main = paste0('EA vs. LA,', ' n = ', dim(mat)[1], ', cp = ',  round(cluster.purity(hc, cls), 2)),
               treeheight_row = 0, 
               treeheight_col = 50,
               width = 6,
               height = 8,
               border_color = NA
)
dev.off()
hc = cutree(p1$tree_col, 2)
cls = anno %>% rownames_to_column('sample') %>% dplyr::filter(sample %in% p1$tree_col$labels) %>% dplyr::mutate(group = as.numeric(group)) %>% deframe()
p2 <- pheatmap(mat,
               breaks = breaksList,
               color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),
               clustering_distance_cols = "manhattan",
               clustering_callback = cl_cb_new,
               annotation_col = anno,
               show_colnames = F,
               annotation_colors = anno_colors[1],
               cluster_rows = T,
               labels_row = my_genes_anno$SYMBOL,
               fontsize_row = fsize(),
               fontsize_col = 3,
               cellwidth = fsize(),
               cellheight = fsize(),
               main = paste0('EA vs. LA,', ' n = ', dim(mat)[1], ', cp = ',  round(cluster.purity(hc, cls), 2)),
               treeheight_row = 0, 
               treeheight_col = 50,
               width = 6.6,
               height = 8,
               border_color = NA
)
dev.off()

################## Importance
set.seed(123)
# my_genes_anno <- my_genes_anno
# data_ <- data_my[, c(intersect(my_genes_anno$SYMBOL, names(data_my)), 'group')]

data_ <- data_my[, c(my_genes_anno$SYMBOL, 'group')] %>% 
  dplyr::filter(group %in% c('EA', 'LA'))
boruta_output <- Boruta(x = data_[, my_genes_anno$SYMBOL != "group"], # protect() solution
                        y = data_$group, doTrace=0)
names(boruta_output)

# 3. Get significant variables including tentatives
boruta_signif <- getSelectedAttributes(boruta_output, withTentative = TRUE)
print(boruta_signif)  

# If you are not sure about the tentative variables being selected for granted, you can choose a TentativeRoughFix on boruta_output.
# Do a tentative rough fix
roughFixMod <- TentativeRoughFix(boruta_output)
boruta_signif <- getSelectedAttributes(roughFixMod)
print(boruta_signif)

# Variable Importance Scores
imps <- attStats(roughFixMod)
imps2 <-  imps[imps$decision != 'Rejected', c('meanImp', 'decision')]
imps2 <-  imps2[!row.names(imps2) %in% 'group', ]
head(imps2[order(-imps2$meanImp), ], 100)  # descending sort

# EA v LA important genes
immune_ImpList_EAvLA <- imps2 %>% dplyr::filter(meanImp >= 3) %>%  rownames_to_column() %>% pull(rowname)


# theme_set(theme_bw(base_size = 12))
imps2 %>% rownames_to_column() %>% dplyr::filter(meanImp >= 3) %>% 
  arrange(meanImp) %>% 
  mutate(name= factor(rowname, levels= rowname)) %>%
  ggplot(aes(y=meanImp, x= name)) + #use d_top30 for top 30 proteins; size is 4x6
  geom_segment(aes(y = 0, x = name, yend = meanImp, xend = name)) + 
  geom_point(stat='identity', color="orange", size=2) +
  labs(x = '', y = 'Importance') +
  ggtitle('Feature importance') +
  theme_ipsum() +
  ylim(0,10) +
  theme(panel.grid = element_blank()) +
  coord_flip()
# ggsave(plot = last_plot(), device = 'pdf', 'IMMUNE_importance_EAvLA.pdf', width = 2.5, height = 3, scale = 1.2)


########################################################################################### LA v C ###########################################################################################
data <- N_LA_v_C_uniprot
gseaRes <- gseaRes_LAvC 
cells <- 'NTC'
comparison <- c('Late.Treated', 'Control')
selection <- gseaRes %>% unnest(leadingEdge) %>% dplyr::select(pathway, leadingEdge) %>% dplyr::rename(SYMBOL = leadingEdge) %>% dplyr::filter(grepl('immune', pathway)) %>% pull(pathway) %>% unique()
cutoff <- 0.01


leading_edge <- gseaRes %>% 
  dplyr::filter(pathway %in% selection) %>% 
  unnest(leadingEdge) %>%
  pull(leadingEdge) %>% 
  unique()



immune_genes <- intersect(rownames(mat_cut_tree), row_anno1$row)

my_genes <- data %>% dplyr::filter(row %in% immune_genes) %>% 
  dplyr::filter(padj <= cutoff, abs(log2FoldChange) > 1.5) %>% 
  pull(row) %>% 
  unique()

my_genes_anno <- data %>% dplyr::filter(row %in% immune_genes) %>% 
  dplyr::filter(padj <= cutoff, abs(log2FoldChange) > 1.5) %>% 
  dplyr::select(row, SYMBOL)


# remove duplicated rows
my_genes_anno <- my_genes_anno[!duplicated(my_genes_anno$row), ]
my_genes_anno <- my_genes_anno %>% as_tibble %>% column_to_rownames(var = 'row')

vsd <- vst(dds)
ids <- metaData %>% dplyr::filter(cellType %in% cells, condition %in% comparison)
mat <- assay(vsd)[my_genes, rownames(ids)]
mat <- mat - rowMeans(mat)
dim(mat)

anno <- metaData %>% dplyr::filter(cellType %in% cells) %>% dplyr::rename(condition, group = condition) %>%  dplyr::select(group) 
anno_colors <- list(group = c(Late.Treated = '#ef9472', Control = '#7bec96'))

#row annotation
row_pathways <- gseaRes %>% 
  dplyr::filter(pathway %in% selection) %>% 
  dplyr::filter(padj < 0.01) %>% 
  unnest(leadingEdge)

row_pathways <- row_pathways %>% dplyr::arrange(leadingEdge, desc(size)) # desc, keeps the largest pathways

row_pathways <- row_pathways[!duplicated(row_pathways$leadingEdge),]
row_anno <- row_pathways %>% dplyr::select(pathway, leadingEdge) %>% dplyr::rename(SYMBOL = leadingEdge)
row_anno <- merge(row_anno, data, by = 'SYMBOL')
row_anno <- row_anno[!duplicated(row_anno$row),]
row_anno <- row_anno %>% dplyr::filter(row %in% my_genes) %>%  dplyr::select(row, pathway) %>% column_to_rownames('row')


breaksList = seq(-2, 2, by = 0.01)


p1 <- pheatmap(mat,
               breaks = breaksList,
               color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),
               clustering_distance_cols = "manhattan",
               clustering_callback = cl_cb_new,
               annotation_col = anno,
               show_colnames = F,
               annotation_colors = anno_colors[1],
               cluster_rows = T,
               labels_row = my_genes_anno$SYMBOL,
               fontsize_row = fsize(),
               fontsize_col = 3,
               cellwidth = fsize(),
               cellheight = fsize(),
               main = paste0('LA vs. C,', ' n = ', dim(mat)[1], ', cp = ',  round(cluster.purity(hc, cls), 2)),
               treeheight_row = 0, 
               treeheight_col = 50,
               border_color = NA,
               width = 6,
               height = 8
)
dev.off()
hc = cutree(p1$tree_col, 2)
cls = anno %>% rownames_to_column('sample') %>% dplyr::filter(sample %in% p1$tree_col$labels) %>% dplyr::mutate(group = as.numeric(group)) %>% deframe()
p2 <- pheatmap(mat,
               breaks = breaksList,
               color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),
               clustering_distance_cols = "manhattan",
               clustering_callback = cl_cb_new,
               annotation_col = anno,
               show_colnames = F,
               annotation_colors = anno_colors[1],
               cluster_rows = T,
               labels_row = my_genes_anno$SYMBOL,
               fontsize_row = fsize(),
               fontsize_col = 3,
               cellwidth = fsize(),
               cellheight = fsize(),
               main = paste0('LA vs. C,', ' n = ', dim(mat)[1], ', cp = ',  round(cluster.purity(hc, cls), 2)),
               treeheight_row = 0, 
               treeheight_col = 50,
               border_color = NA,
               width = 6,
               height = 8
)
dev.off()

################## Importance
set.seed(123)
data_ <- data_my[, c(my_genes_anno$SYMBOL, 'group')] %>% 
  dplyr::filter(group %in% c('LA', 'C'))
boruta_output <- Boruta(x = data_[, my_genes_anno$SYMBOL != "group"], # protect() solution
                        y = data_$group, doTrace=0)
names(boruta_output)

# 3. Get significant variables including tentatives
boruta_signif <- getSelectedAttributes(boruta_output, withTentative = TRUE)
print(boruta_signif)  

# If you are not sure about the tentative variables being selected for granted, you can choose a TentativeRoughFix on boruta_output.
# Do a tentative rough fix
roughFixMod <- TentativeRoughFix(boruta_output)
boruta_signif <- getSelectedAttributes(roughFixMod)
print(boruta_signif)

# Variable Importance Scores
imps <- attStats(roughFixMod)
imps2 <-  imps[imps$decision != 'Rejected', c('meanImp', 'decision')]
imps2 <-  imps2[!row.names(imps2) %in% 'group', ]
head(imps2[order(-imps2$meanImp), ], 100)  # descending sort

# LA v C important genes
immune_ImpList_LAvC <- imps2 %>% dplyr::filter(meanImp >= 3) %>%  rownames_to_column() %>% pull(rowname)


imps2 %>% rownames_to_column() %>% dplyr::filter(meanImp >= 3) %>% 
  arrange(meanImp) %>% 
  mutate(name= factor(rowname, levels= rowname)) %>%
  ggplot(aes(y=meanImp, x= name)) + 
  geom_segment(aes(y = 0, x = name, yend = meanImp, xend = name)) + 
  geom_point(stat='identity', color="orange", size=2) +
  labs(x = '', y = 'Importance') +
  ggtitle('Feature importance') +
  theme_ipsum() +
  ylim(0,10) +
  theme(panel.grid = element_blank()) +
  coord_flip()

# ggsave(plot = last_plot(), device = 'pdf', 'IMMUNE_importance_LA.pdf', width = 2.5, height = 3, scale = 1.2)
