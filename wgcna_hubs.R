library(tidyverse)
library(org.Hs.eg.db)
library(hrbrthemes)
library(ggthemes)
library(extrafont)
library('Homo.sapiens')
library(WGCNA)
library(corrr)
library(gplots)
library(RColorBrewer)

### Create network data for analysis with cytoscape (to identify hub genes)

# helper function to export network data
export_network_to_graphml <- function (adj_mat, filename=NULL, weighted=TRUE,
                                       threshold=0.5, max_edge_ratio=3,
                                       nodeAttr=NULL, nodeAttrDataFrame=NULL,
                                       edgeAttributes=NULL, verbose=FALSE) {
  library('igraph')
  
  # Determine filename to use
  if (is.null(filename)) {
    filename='network.graphml'
  }
  
  
  # Adjust threshold if needed to limit remaining edges
  max_edges <- max_edge_ratio * nrow(adj_mat)
  
  edge_to_total_ratio <- max_edges / length(adj_mat)
  edge_limit_cutoff <- as.numeric(quantile(abs(adj_mat), 1 - edge_to_total_ratio))
  
  # Also choose a minimum threshold to make sure that at least some edges
  # are left
  min_threshold <- as.numeric(quantile(abs(adj_mat), 0.9999))
  
  threshold <- min(min_threshold, max(threshold, edge_limit_cutoff))
  
  # Remove edges with weights lower than the cutoff
  adj_mat[abs(adj_mat) < threshold] <- 0
  
  # Drop any genes with no edges
  orphaned <- (colSums(adj_mat) == 0)
  adj_mat <- adj_mat[!orphaned, !orphaned]
  
  # Also remove annotation entries
  if (!is.null(nodeAttr)) {
    nodeAttr <- nodeAttr[!orphaned]
  }
  if (!is.null(nodeAttrDataFrame)) {
    nodeAttrDataFrame <- nodeAttrDataFrame[!orphaned,]
  }
  
  # Keep track of non-positive edges and rescale to range 0,1
  is_zero     <- adj_mat == 0
  is_negative <- adj_mat < 0
  
  adj_mat <- (abs(adj_mat) - threshold) / (max(adj_mat) - threshold)
  adj_mat[is_zero] <- 0
  adj_mat[is_negative] <- -adj_mat[is_negative]
  
  if (verbose) {
    message(sprintf("Outputting matrix with %d nodes and %d edges", 
                    nrow(adj_mat), sum(adj_mat > 0)))
  }
  
  # Create a new graph and add vertices
  # Weighted graph
  if (weighted) {
    g <- graph.adjacency(adj_mat, mode='undirected', weighted=TRUE, diag=FALSE)
  } else {
    adj_mat[adj_mat != 0] <- 1
    g <- graph.adjacency(adj_mat, mode='undirected', diag=FALSE)
  }
  
  # Add single node annotation from vector
  if (!is.null(nodeAttr)) {
    g <- set.vertex.attribute(g, "attr", value=nodeAttr)
  }
  
  # Add node one or more node annotations from a data frame
  if (!is.null(nodeAttrDataFrame)) {
    for (colname in colnames(nodeAttrDataFrame)) {
      g <- set.vertex.attribute(g, colname, value=nodeAttrDataFrame[,colname])
    }
  }
  
  edge_correlation_negative <- c()
  
  # neg_correlations[edge_list]
  edge_list <- get.edgelist(g)
  
  for (i in 1:nrow(edge_list)) {
    from <- edge_list[i, 1]    
    to   <- edge_list[i, 2]    
  }
  
  # Save graph to a file
  write.graph(g, filename, format='graphml')
  
  # return igraph
  return(g)
}

##############                      ##############
##############      EA Analysis     ############## 
##############                      ##############

EAvC_DEG <- N_EA_v_C_uniprot %>% dplyr::filter(padj<= 0.01, abs(log2FoldChange) > 1.5) %>% pull(row)
EAvLA_DEG <- N_EA_v_LA_uniprot %>% dplyr::filter(padj<= 0.01, abs(log2FoldChange) > 1.5) %>% pull(row)
EA_DEG <- c(EAvC_DEG, EAvLA_DEG) %>% unique()
anno_genes_ensembl <- intersect(rownames(mat_cut_tree), rownames(row_anno))
EA_anno_genes <- intersect(EA_DEG, anno_genes_ensembl)
length(EA_anno_genes)
# >>112 genes

# annotated genes matrix
all_anno_EA <- counts(dds, normalized=TRUE) %>%
  data.frame(check.names = F) %>%
  dplyr::select(contains('N')) %>% 
  dplyr::select(matches('EA')) %>% 
  rownames_to_column(var="gene") %>% 
  dplyr::filter(gene %in% EA_anno_genes) %>% # EA annotated genes
  mutate(symbol= mapIds(org.Hs.eg.db,
                        keys = gene, 
                        keytype ="ENSEMBL",
                        column = "SYMBOL")) %>% 
  dplyr::select(!gene) %>% 
  distinct(symbol, .keep_all = T) %>% 
  column_to_rownames('symbol') %>% t()

# log transform
all_anno_EA <- log2(all_anno_EA+1)

matrixALL_anno_EA <- all_anno_EA # annotated genes


# Construct adjacency matrix for annonated genes
sq_EA <- matrixALL_anno_EA %>% correlate(use = 'complete.obs', method = 'spearman') %>% rearrange() %>% dplyr::select(!term)
adj_matrix_EA <- adjacency.fromSimilarity(as.matrix(sq_EA), power=12, type='signed')
adj_matrix_EA <- matrix(adj_matrix_EA, nrow=nrow(adj_matrix_EA))

rownames(adj_matrix_EA) <- colnames(matrixALL_anno_EA)
colnames(adj_matrix_EA) <- colnames(matrixALL_anno_EA)

adj_matrix_EA[is.na(adj_matrix_EA)] <- 1

g_EA_noMods <- export_network_to_graphml(adj_matrix_EA, filename='network_ЕА_annotated_0.5_noMods.graphml',
                               threshold=0.5)

##############                      ##############
##############      LA Analysis     ############## 
##############                      ##############

# annotated genes LA

LAvC_DEG <- N_LA_v_C_uniprot %>% dplyr::filter(padj<= 0.01, abs(log2FoldChange) > 1.5) %>% pull(row)
LA_DEG <- c(LAvC_DEG, EAvLA_DEG) %>% unique()
anno_genes_ensembl <- intersect(rownames(row_anno), rownames(mat))
LA_anno_genes <- intersect(LA_DEG, anno_genes_ensembl)
length(LA_anno_genes)
# >>176 genes

all_anno_LA <- counts(dds, normalized=TRUE) %>%
  data.frame(check.names = F) %>%
  dplyr::select(contains('N')) %>% 
  dplyr::select(matches('K')) %>% 
  rownames_to_column(var="gene") %>% 
  dplyr::filter(gene %in% LA_anno_genes) %>% # LA annotaed genes
  mutate(symbol= mapIds(org.Hs.eg.db,
                        keys = gene, 
                        keytype ="ENSEMBL",
                        column = "SYMBOL")) %>% 
  dplyr::select(!gene) %>% 
  distinct(symbol, .keep_all = T) %>% 
  column_to_rownames('symbol') %>% t()

# log transform
all_anno_LA <- log2(all_anno_LA+1)

matrixALL_anno_LA <- all_anno_LA # annotated genes


# Construct adjacency matrix for annonated genes

sq_LA <- matrixALL_anno_LA %>% correlate(use = 'complete.obs', method = 'spearman') %>% rearrange() %>% dplyr::select(!term)
adj_matrix_LA <- adjacency.fromSimilarity(as.matrix(sq_LA), power=12, type='signed')
adj_matrix_LA <- matrix(adj_matrix_LA, nrow=nrow(adj_matrix_LA))

rownames(adj_matrix_LA) <- colnames(matrixALL_anno_LA)
colnames(adj_matrix_LA) <- colnames(matrixALL_anno_LA)


adj_matrix_LA[is.na(adj_matrix_LA)] <- 1


g_LA_noMods <- export_network_to_graphml(adj_matrix_LA, filename='network_LА_annotated_0.5_noMods.graphml',
                                       threshold= 0.5)


##############                      ##############
##############      EA-LA Analysis  ############## 
##############                      ##############


EAvLA_DEG <- N_EA_v_LA_uniprot %>% dplyr::filter(padj<= 0.01, abs(log2FoldChange) > 1.5) %>% pull(row)
EAvLA_DEG <- EAvLA_DEG %>% unique()
anno_genes_ensembl <- intersect(rownames(row_anno), rownames(mat))
EAvLA_anno_genes <- intersect(EAvLA_DEG, anno_genes_ensembl)
length(EAvLA_anno_genes)
# >>30 genes

all_anno_EALA <- counts(dds, normalized=TRUE) %>%
  data.frame(check.names = F) %>%
  dplyr::select(contains('N')) %>% 
  dplyr::select(matches('EA')) %>% 
  rownames_to_column(var="gene") %>% 
  dplyr::filter(gene %in% EAvLA_anno_genes) %>% # LA annotaed genes
  mutate(symbol= mapIds(org.Hs.eg.db,
                        keys = gene, 
                        keytype ="ENSEMBL",
                        column = "SYMBOL")) %>% 
  dplyr::select(!gene) %>% 
  distinct(symbol, .keep_all = T) %>% 
  column_to_rownames('symbol') %>% t()

# log transform
all_anno_EALA <- log2(all_anno_EALA+1)


matrixALL_anno_EALA <- all_anno_EALA # annotated genes


# Construct adjacency matrix for annonated genes

sq_EALA <- matrixALL_anno_EALA %>% correlate(use = 'complete.obs', method = 'spearman') %>% rearrange() %>% dplyr::select(!term)
adj_matrix_EALA <- adjacency.fromSimilarity(as.matrix(sq_EALA), power=12, type='signed')
adj_matrix_EALA <- matrix(adj_matrix_EALA, nrow=nrow(adj_matrix_EALA))

rownames(adj_matrix_EALA) <- colnames(matrixALL_anno_EALA)
colnames(adj_matrix_EALA) <- colnames(matrixALL_anno_EALA)

adj_matrix_EALA[is.na(adj_matrix_EALA)] <- 1



g_EALA_noMods <- export_network_to_graphml(adj_matrix_EALA, filename='network_EALА_annotated_0.5_noMods.graphml',
                                         threshold= 0.5)
#############################################################################################################
# hub gene table
############################################################################################################

########## hub genes enrichment analysis EA
analyzed_net <- read.csv('~/EA_annotated_0.5.csv')
hubs <- analyzed_net %>% 
  dplyr::filter(Degree >= 2) %>%
  dplyr::select(name, Degree) %>% 
  arrange(desc(Degree)) 

hubs_EA <- analyzed_net %>% dplyr::filter(Degree >= 2) %>% dplyr::select(name, Degree) %>% arrange(desc(Degree)) # Degree was > 1 in original
sif_EA <- read.delim('~/EA_annotated_0.5.sif', sep = '\t', header = F)
network <-sif_EA %>% dplyr::select(V1, V3)

# network with hubs
sif_EA %>% dplyr::filter(V1 %in% hubs_EA$name | V3 %in% hubs_EA$name) %>% write_tsv(., '~/EA_FC1.5_hubs.sif')


dbs <- c('GO_Biological_Process_2018')
hubGO <- data.frame(gene = NA,
                    GO = NA)

for (i in hubs$name) {
  
  networkSelect <- network %>% dplyr::filter(V1 %in% i | V3 %in% i) 
  networkSelect <- c(networkSelect$V1, networkSelect$V3) %>% unique()
  
  a <- as.data.frame(enrichR::enrichr(networkSelect, dbs))
  a <- a %>% dplyr::filter(GO_Biological_Process_2018.Adjusted.P.value <= 0.05) %>% 
    arrange(GO_Biological_Process_2018.Adjusted.P.value)
  # arrange(GO_Biological_Process_2018.Combined.Score)
  b <- sif_LA %>% dplyr::filter(V1 %in% i | V3  %in% i)
  b <- c(b$V1, b$V3) %>% unique()
  b <- b[b != i]
  b <- paste(b, collapse = ", ")
  temp = data.frame(gene = i,
                    GO = a[1,1],
                    nodes = b)
  temp = data.frame(gene = i,
                    GO = a[1,1])
  hubGO = rbind(temp, hubGO)
  
}

hubGO$GO <- gsub(' \\(GO:[0-9]+\\)', '', hubGO$GO)
hubGO[1:length(row.names(hubGO)) -1,] %>% view()

hub_table_ea <- hubs %>% left_join(hubGO, by = c('name' = 'gene')) %>% 
  dplyr::rename(gene = 'name') %>% 
  dplyr::rename(degree = 'Degree') %>% 
  relocate(gene) %>% 
  na.omit()



############ EA v LA
## hub genes enrichment analysis
analyzed_net <- read.csv('~/EAvLA_annotated_0.5.csv')
hubs <- analyzed_net %>% 
  dplyr::filter(Degree >= 2) %>%
  dplyr::select(name, Degree) %>% 
  arrange(desc(Degree))  # Degree was > 1 in original
  # head(5)

hubs_EALA <- analyzed_net %>% dplyr::filter(Degree >= 2) %>% dplyr::select(name, Degree) %>% arrange(desc(Degree)) # Degree was > 1 in original
sif_EALA <- read.delim('~/EAvLA_annotated_0.5.sif', sep = '\t', header = F)
network <-sif_EALA  %>% dplyr::select(V1, V3)

# network with hubs
sif_EALA %>% dplyr::filter(V1 %in% hubs_EALA$name | V3 %in% hubs_EALA$name) %>% write_tsv(., '~/EAvLA_FC1.5_hubs.sif')

dbs <- c('GO_Biological_Process_2018')
hubGO <- data.frame(gene = NA,
                    GO = NA)

# i = 'E2F1'
for (i in hubs$name) {
  
  networkSelect <- network %>% dplyr::filter(V1 %in% i | V3 %in% i) 
  networkSelect <- c(networkSelect$V1, networkSelect$V3) %>% unique()
  
  a <- as.data.frame(enrichR::enrichr(networkSelect, dbs))
  a <- a %>% dplyr::filter(GO_Biological_Process_2018.Adjusted.P.value <= 0.05) %>%
    arrange(GO_Biological_Process_2018.Adjusted.P.value)
  # arrange(GO_Biological_Process_2018.Combined.Score)
  b <- sif_LA %>% dplyr::filter(V1 %in% i | V3  %in% i)
  b <- c(b$V1, b$V3) %>% unique()
  b <- b[b != i]
  b <- paste(b, collapse = ", ")
  temp = data.frame(gene = i,
                    GO = a[1,1],
                    nodes = b)
  temp = data.frame(gene = i,
                    GO = a[1,1])
  hubGO = rbind(temp, hubGO)
  
}

hubGO$GO <- gsub(' \\(GO:[0-9]+\\)', '', hubGO$GO)
hubGO[1:length(row.names(hubGO)) -1,] %>% view()

hub_table_eala <- hubs %>% left_join(hubGO, by = c('name' = 'gene')) %>% 
  dplyr::rename(gene = 'name') %>% 
  dplyr::rename(degree = 'Degree') %>% 
  relocate(gene) %>% 
  na.omit()

# top hubs with WCGNA commands
# kin_EALA <- softConnectivity(adj_matrix_EALA)
# kin_hubs_EALA = (rank (-kin) <= 10) # top 10 hubs
# rownames(adj_matrix_EALA)[kin_hubs_EALA]


############ LA
## hub genes enrichment analysis
analyzed_net <- read.csv('~/LA_annoated_0.5.csv')
hubs <- analyzed_net %>% 
  dplyr::filter(Degree >= 2) %>%
  dplyr::select(name, Degree) %>% 
  arrange(desc(Degree))  # Degree was > 1 in original
# head(5)

hubs_LA <- analyzed_net %>% dplyr::filter(Degree >= 2) %>% dplyr::select(name, Degree) %>% arrange(desc(Degree)) # Degree was > 1 in original
sif_LA <- read.delim('~/LA_annoated_0.5.sif', sep = '\t', header = F)
network <-sif_LA  %>% dplyr::select(V1, V3)

# network with hubs
sif_LA %>% dplyr::filter(V1 %in% hubs_LA$name | V3 %in% hubs_LA$name) %>% write_tsv(., '~/LA_FC1.5_hubs.sif')


dbs <- c('GO_Biological_Process_2018')
hubGO <- data.frame(gene = NA,
                    GO = NA,
                    nodes = NA)

for (i in hubs$name) {
  
  networkSelect <- network %>% dplyr::filter(V1 %in% i | V3 %in% i) 
  networkSelect <- c(networkSelect$V1, networkSelect$V3) %>% unique()
  
  a <- as.data.frame(enrichR::enrichr(networkSelect, dbs))
  a <- a %>% dplyr::filter(GO_Biological_Process_2018.Adjusted.P.value <= 0.05) %>%
    arrange(GO_Biological_Process_2018.Adjusted.P.value)
  b <- sif_LA %>% dplyr::filter(V1 %in% i | V3  %in% i)
  b <- c(b$V1, b$V3) %>% unique()
  b <- b[b != i]
  b <- paste(b, collapse = ", ")
  temp = data.frame(gene = i,
                    GO = a[1,1],
                    nodes = b)
  hubGO = rbind(temp, hubGO)
  
}

hubGO$GO <- gsub(' \\(GO:[0-9]+\\)', '', hubGO$GO)
hubGO[1:length(row.names(hubGO)) -1,] %>% view()

hub_table_la <- hubs %>% left_join(hubGO, by = c('name' = 'gene')) %>% 
  dplyr::rename(gene = 'name') %>% 
  dplyr::rename(degree = 'Degree') %>% 
  relocate(gene) %>% 
  na.omit()



hubs_ea_10 <- hub_table_ea %>% dplyr::mutate(group = 'EAvC') %>% relocate(nodes, .after = group) %>% head(4)
hubs_la_10 <- hub_table_la %>%  dplyr::mutate(group = 'LAvC') %>% relocate(nodes, .after = group) %>%  head(4)
hubs_eala_10 <- hub_table_eala %>% dplyr::mutate(group = 'EAvLA') %>% relocate(nodes, .after = group) %>% head(1)

hubs_both <- rbind(hubs_ea_10, hubs_la_10, hubs_eala_10)
library(knitr)
library(rmarkdown)

table<-kable(hubs_both, format="markdown")
cat(table, sep="\n", file="table.Rmd")
render("table.Rmd",output_format = "pdf_document")
