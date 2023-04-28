# Transcriptomic profile of naïve CD4+ T cells in HIV-1 infected patients initiating antiretroviral therapy at acute or chronic phase of infection

We analyzed the whole transcriptome characteristics of blood CD4+ T naïve (T<sub>N</sub>) cells isolated from HIV-1 infected patients starting ART at acute (early ART = EA; n = 13) or chronic (late ART = LA; n = 11) phase of infection and controls (C; n = 15). RNA sequencing revealed 389 differentially expressed genes (DEGs) in EA and 810 in LA group in relation to controls. Comparison of the two groups of patients showed 183 DEGs. We focused on DEGs involved in apoptosis, inflammation and immune response. Clustering showed a poor separation of EA from C suggesting that these two groups present a similar transcriptomic profile of CD4+ T<sub>N</sub> cells. The comparison of EA and LA patients resulted in a high cluster purity revealing that different biological dysfunctions characterize EA and LA patients. The upregulated expression of several inflammatory chemokine genes distinguished the patient groups from C; CCL2 and CCL7, however, were downregulated in EA compared to LA patients. BCL2, an anti-apoptotic factor pivotal for naïve T cell homeostasis, distinguished both EA and LA from C. The expression of several DEGs involved in different inflammatory processes (TLR4, PTGS2, RAG1, IFNA16) was lower in EA compared LA. We conclude that although the transcriptome of CD4+ TN cells isolated from patients initiating ART at acute infection reveals a more quiescent phenotype, the survival profile of these cells still appears to be affected. Our results show that the detrimental process of inflammation is under more efficient control in EA patients.

## Differential gene expression analysis
![DE analysis](https://ars.els-cdn.com/content/image/1-s2.0-S0888754321003219-gr1.jpg "DE analysis")

A. Volcano plots showing the expression magnitude (log<sub>2</sub>FC) of differentially expressed genes (DEGs) and statistical significance. Genes that are statistically significant (adj. p < 0.01) are shown in blue (absolute log2FC < 1.5) and red (absolute log<sub>2</sub>FC > 1.5); nonsignificant genes are shown in grey. Highly modulated (absolute log<sub>2</sub>FC > 1.5) but not statistically significant (adj. p < 0.01) are colored in pink.

B. Euler plot of the overlap in significantly DEGs between EA and LA patients, EA and controls and LA and controls.

C. Violin plot showing the median overall difference in gene expression (absolute log<sub>2</sub>FC) in EA and LA patients compared to healthy controls. **** = p < 0.0001.

## Clustering and ontology of differntially expressed genes
![Clustering](https://ars.els-cdn.com/content/image/1-s2.0-S0888754321003219-gr2.jpg "a title")

Heatmap and pie chart visualizing 1134 highly significant (adj.p < 0.01, log2FC > 1.5) DEGs in HIV-1-infected patients compared to controls.

A. The dendrogram labels denote the 4 different clusters defined by hierarchical clustering. Colored lanes to the left indicate the functional annotation of genes involved in apoptosis, inflammation and immune response (n = 211).

B. Significantly enriched GO terms in the list of DEGs not involved in apoptosis, inflammation and immune response. The list of 923 genes was analyzed for GO term enrichment. Similarity between significantly enriched terms (p < 0.01) was calculated using rrvigo and visualized as a pie chart. The percentage indicates the slice's proportion of the whole pie. The number of genes associated with each GO term is shown in parentheses.

## Functional characterization of DEGs within clusters
![functional analysis](https://ars.els-cdn.com/content/image/1-s2.0-S0888754321003219-gr3.jpg "a title")

Expression levels and functional characterization of DEGs within the clusters associated with apoptosis, inflammation and immune response.

A. Violin plots representing the overall gene expression in each of the 4 defined clusters containing annotated genes.

B. Dot-plot visualization of the GO analysis of genes within the 4 clusters. The size of the dot corresponds to the genes associated with the respective GO term divided by the total number of genes in the cluster (gene ratio) and the color represents the significance level.

## Heatmap visualization and random forest feature selection
We performed random forest feature selection implemented by the Boruta R package to single out essential genes controlling the processes of apoptosis, inflammation and immune response during early and late treated HIV-1 infection. The DEGs between groups are shown in heat-maps and the functions of the proteins coded by DEGs according to Uniprot.

## Weighted network correlation analysis
In order to find possible DEGs potentially involved in pathological processes in each group of patients we analyzed co-regulation and co-expression of the gene-sets specific for the comparisons EA vs C, LA vs C and EA v LA patients. This method identified hub genes, which we defined as genes whose expression significantly correlated (r > 0.5) with at least 2 additional genes. To gain insight into the function of these gene hubs we performed enrichment analyses using the co-expressed genes.

## Full article

Petkov S, Chiodi F. Distinct transcriptomic profiles of naïve CD4+ T cells distinguish HIV-1 infected patients initiating antiretroviral therapy at acute or chronic phase of infection. Genomics. 2021 Nov;113(6):3487-3500. doi: 10.1016/j.ygeno.2021.08.014. Epub 2021 Aug 20. PMID: 34425224.
