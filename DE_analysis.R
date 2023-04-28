### Load required packages:

library(DESeq2)
library(readxl)
library(cowplot)
library(tidyverse)
library(org.Hs.eg.db)
library(EnhancedVolcano)
library(eulerr)
library(hrbrthemes)
library(ggthemes)
library(extrafont)
library('Homo.sapiens')
library(RColorBrewer)


### Load raw counts data/metadata and clean

countData <- as.matrix(read.csv('degData_naive.csv', sep = ';', row.names = 'gene_id'))
metaData <- read.csv('colData_naive.csv', sep = ';', row.names = 1)
metaData <- metaData[,c("cellType","condition")]
metaData$condition <- factor(metaData$condition)
metaData$cellType <- factor(metaData$cellType)
colnames(countData) <- sub("X", "", colnames(countData))
all(rownames(metaData) == colnames(countData))



### Construct DESEQDataSet Object

dds <- DESeqDataSetFromMatrix(countData = countData, 
                              colData = metaData, 
                              design =  ~ cellType + cellType:condition)

dds$group <- factor(paste0(dds$cellType, dds$condition))
design(dds) <- ~ group

dds <- DESeq(dds)
resultsNames(dds)

# Results for Naive CD4 T cells
N_LA_v_C <- results(dds, contrast=c("group", "NTCLate.Treated", "NTCControl"), tidy = T)
N_EA_v_C<- results(dds, contrast=c("group", "NTCEarly.Treated", "NTCControl"), tidy = T)
N_EA_v_LA <- results(dds, contrast=c("group","NTCEarly.Treated","NTCLate.Treated"), tidy = T)


### Volcano plot of DEGs between naive CD4+ T cells of Early treated (EA) and controls (C)

# ENSEMBL annotation
ens2symbol <- AnnotationDbi::select(org.Hs.eg.db,
                                    key=N_EA_v_C$row, 
                                    columns="SYMBOL",
                                    keytype="ENSEMBL")
# UNIPROT annotation
ens2uniprot <- AnnotationDbi::select(org.Hs.eg.db,
                                    key=N_EA_v_C$row,
                                    columns="UNIPROT",
                                    keytype="ENSEMBL")

# ENTREZ annotation
ens2entrez <- AnnotationDbi::select(org.Hs.eg.db,
                                    key=N_EA_v_C$row,
                                    columns="ENTREZID",
                                    keytype="ENSEMBL")

ens2symbol <- as_tibble(ens2symbol)
ens2uniprot <- as_tibble(ens2uniprot)
ens2entrez <- as_tibble(ens2entrez)

N_EA_v_C <- inner_join(N_EA_v_C, ens2symbol, by=c("row"="ENSEMBL"))
N_EA_v_C <- N_EA_v_C %>% na.omit()
N_EA_v_C <- N_EA_v_C[!duplicated(N_EA_v_C$SYMBOL), ]

# PROTEIN CODING GENES
N_EA_v_C_uniprot <- inner_join(N_EA_v_C, ens2uniprot, by=c("row"="ENSEMBL"))
N_EA_v_C_uniprot <- N_EA_v_C_uniprot %>% na.omit()
N_EA_v_C_uniprot <- N_EA_v_C_uniprot[!duplicated(N_EA_v_C_uniprot$SYMBOL), ]

# ENTREZ
N_EA_v_C_entrez <- inner_join(N_EA_v_C, ens2entrez, by=c("row"="ENSEMBL"))
N_EA_v_C_entrez <- N_EA_v_C_entrez %>% na.omit()
N_EA_v_C_entrez <- N_EA_v_C_entrez[!duplicated(N_EA_v_C_entrez$SYMBOL), ]

# summary parameters
pCutoff = 0.01
FCcutoff = 1.5

vol_eavc <- EnhancedVolcano(N_EA_v_C_uniprot,
                lab = NA, 
                x = 'log2FoldChange', 
                y = 'padj', 
                pCutoff = 0.01, #pCutoff = 0.005,
                FCcutoff = 1.5, #FCcutoff = 1
                col = c('grey', 'pink','#00509d', '#bf0603'),
                title = 'EA vs. Control',
                subtitle = '',
                caption = '',
                xlim = c(-5, 5),
                ylim = c(0, 10),
                pointSize = 1,
                gridlines.major = F,
                gridlines.minor = F,
                legendPosition = 'none',
                legendLabSize = 10
  )
  
#statistics to add to plot quadrants
N_EA_v_C_uniprot %>% 
  summarise(
    dn_red = sum(log2FoldChange < -(FCcutoff) & padj <= pCutoff),
    dn_pink = sum(log2FoldChange < -(FCcutoff) & padj >= pCutoff),
    dn_blue = sum(log2FoldChange> -(FCcutoff) & log2FoldChange<0 & padj <= pCutoff),
    dn_grey = sum(log2FoldChange> -(FCcutoff) & log2FoldChange<0 & padj >= pCutoff),
    up_red = sum(log2FoldChange > FCcutoff & padj <= pCutoff),
    up_pink = sum(log2FoldChange > FCcutoff & padj >= pCutoff),
    up_blue = sum(log2FoldChange < FCcutoff & log2FoldChange > 0 & padj <= pCutoff),
    up_grey = sum(log2FoldChange < FCcutoff & log2FoldChange > 0 & padj >= pCutoff)
  )



### Volcano plot of DEGs between naive CD4+ T cells of Early treated (EA) and Late treated (LA)

# ENSEMBL annotation
ens2symbol <- AnnotationDbi::select(org.Hs.eg.db,
                                    key=N_EA_v_LA$row, 
                                    columns="SYMBOL",
                                    keytype="ENSEMBL")
# UNIPROT anno
ens2uniprot <- AnnotationDbi::select(org.Hs.eg.db,
                                    key=N_EA_v_LA$row,
                                    columns="UNIPROT",
                                    keytype="ENSEMBL")

# ENTREZ anno
ens2entrez <- AnnotationDbi::select(org.Hs.eg.db,
                                    key=N_EA_v_LA$row,
                                    columns="ENTREZID",
                                    keytype="ENSEMBL")

ens2symbol <- as_tibble(ens2symbol)
ens2uniprot <- as_tibble(ens2uniprot)
ens2entrez <- as_tibble(ens2entrez)


N_EA_v_LA <- inner_join(N_EA_v_LA, ens2symbol, by=c("row"="ENSEMBL"))
N_EA_v_LA <- N_EA_v_LA %>% na.omit()
N_EA_v_LA <- N_EA_v_LA[!duplicated(N_EA_v_LA$SYMBOL), ]

# PROTEIN CODING GENES
N_EA_v_LA_uniprot <- inner_join(N_EA_v_LA, ens2uniprot, by=c("row"="ENSEMBL"))
N_EA_v_LA_uniprot <- N_EA_v_LA_uniprot %>% na.omit()
N_EA_v_LA_uniprot <- N_EA_v_LA_uniprot[!duplicated(N_EA_v_LA_uniprot$SYMBOL), ]

# ENTREZ
N_EA_v_LA_entrez <- inner_join(N_EA_v_LA, ens2entrez, by=c("row"="ENSEMBL"))
N_EA_v_LA_entrez <- N_EA_v_LA_entrez %>% na.omit()
N_EA_v_LA_entrez <- N_EA_v_LA_entrez[!duplicated(N_EA_v_LA_entrez$SYMBOL), ]

vol_eavla <- EnhancedVolcano(N_EA_v_LA_uniprot,
                lab = NA, 
                x = 'log2FoldChange', 
                y = 'padj', 
                col = c('grey', 'pink','#00509d', '#bf0603'),
                pCutoff = 0.01, #pCutoff = 0.005,
                FCcutoff = 1.5, #FCcutoff = 1
                title = 'EA vs. LA',
                subtitle = '',
                caption = '',
                xlim = c(-5, 5),
                ylim = c(0, 10),
                pointSize = 1,
                gridlines.major = F,
                gridlines.minor = F,
                legendPosition = 'none',
                legendLabSize = 10
  )
  
#statistics to add to plot quadrants
N_EA_v_LA_uniprot %>% 
  summarise(
    dn_red = sum(log2FoldChange < -(FCcutoff) & padj <= pCutoff),
    dn_pink = sum(log2FoldChange < -(FCcutoff) & padj >= pCutoff),
    dn_blue = sum(log2FoldChange> -(FCcutoff) & log2FoldChange<0 & padj <= pCutoff),
    dn_grey = sum(log2FoldChange> -(FCcutoff) & log2FoldChange<0 & padj >= pCutoff),
    up_red = sum(log2FoldChange > FCcutoff & padj <= pCutoff),
    up_pink = sum(log2FoldChange > FCcutoff & padj >= pCutoff),
    up_blue = sum(log2FoldChange < FCcutoff & log2FoldChange > 0 & padj <= pCutoff),
    up_grey = sum(log2FoldChange < FCcutoff & log2FoldChange > 0 & padj >= pCutoff)
  )


### Volcano plot of DEGs between naive CD4+ T cells of Late treated (LA) and controls (C)

# ENSEMBL annotation
ens2symbol <- AnnotationDbi::select(org.Hs.eg.db,
                                    key=N_LA_v_C$row, 
                                    columns="SYMBOL",
                                    keytype="ENSEMBL")
# UNIPROT anno
ens2uniprot <- AnnotationDbi::select(org.Hs.eg.db,
                                    key=N_LA_v_C$row,
                                    columns="UNIPROT",
                                    keytype="ENSEMBL")

# ENTREZ anno
ens2entrez <- AnnotationDbi::select(org.Hs.eg.db,
                                    key=N_LA_v_C$row,
                                    columns="ENTREZID",
                                    keytype="ENSEMBL")

ens2symbol <- as_tibble(ens2symbol)
ens2uniprot <- as_tibble(ens2uniprot)
ens2entrez <- as_tibble(ens2entrez)


N_LA_v_C <- inner_join(N_LA_v_C, ens2symbol, by=c("row"="ENSEMBL"))
N_LA_v_C <- N_LA_v_C %>% na.omit()
N_LA_v_C <- N_LA_v_C[!duplicated(N_LA_v_C$SYMBOL), ]

# PROTEIN CODING GENES
N_LA_v_C_uniprot <- inner_join(N_LA_v_C, ens2uniprot, by=c("row"="ENSEMBL"))
N_LA_v_C_uniprot <- N_LA_v_C_uniprot %>% na.omit()
N_LA_v_C_uniprot <- N_LA_v_C_uniprot[!duplicated(N_LA_v_C_uniprot$SYMBOL), ]

# ENTREZ
N_LA_v_C_entrez <- inner_join(N_LA_v_C, ens2entrez, by=c("row"="ENSEMBL"))
N_LA_v_C_entrez <- N_LA_v_C_entrez %>% na.omit()
N_LA_v_C_entrez <- N_LA_v_C_entrez[!duplicated(N_LA_v_C_entrez$SYMBOL), ]

vol_lavc <- EnhancedVolcano(N_LA_v_C_uniprot,
                lab = NA, 
                x = 'log2FoldChange', 
                y = 'padj', 
                pCutoff = 0.01, #pCutoff = 0.005,
                FCcutoff = 1.5, #FCcutoff = 1
                col = c('grey', 'pink','#00509d', '#bf0603'),
                title = 'LA vs. Control',
                subtitle = '',
                caption = '',
                xlim = c(-5, 5),
                ylim = c(0, 10),
                pointSize = 1,
                gridlines.major = F,
                gridlines.minor = F,
                legendPosition = 'none',
                legendLabSize = 10
  )
  
#statistics to add to plot quadrants
N_LA_v_C_uniprot %>% 
  summarise(
    dn_red = sum(log2FoldChange < -(FCcutoff) & padj <= pCutoff),
    dn_pink = sum(log2FoldChange < -(FCcutoff) & padj >= pCutoff),
    dn_blue = sum(log2FoldChange> -(FCcutoff) & log2FoldChange<0 & padj <= pCutoff),
    dn_grey = sum(log2FoldChange> -(FCcutoff) & log2FoldChange<0 & padj >= pCutoff),
    up_red = sum(log2FoldChange > FCcutoff & padj <= pCutoff),
    up_pink = sum(log2FoldChange > FCcutoff & padj >= pCutoff),
    up_blue = sum(log2FoldChange < FCcutoff & log2FoldChange > 0 & padj <= pCutoff),
    up_grey = sum(log2FoldChange < FCcutoff & log2FoldChange > 0 & padj >= pCutoff)
  )
