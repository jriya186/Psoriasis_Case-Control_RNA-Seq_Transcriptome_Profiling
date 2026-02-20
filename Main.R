source("/Users/jriya/Desktop/RNA-seq_psoriasis/functions.R")

library(tidyverse)

# Loading counts data
counts <- read.delim("data/GSE54456_raw_counts_GRCh38.p13_NCBI.tsv", row.names = 1)

# Loading metadata file from SRA for sample information
metadata <- read.csv("data/SraRunTable.csv")

# checks
dim(counts) # 39376 X 171
head(counts[1:5, 1:5])
head(metadata)

# Setting GSM ID as row names in metadata
# rownames(metadata) <- metadata$Sample.Name

# sum(duplicated(metadata$Sample.Name)) # 5 duplicated samples
# ‘GSM1315618’, ‘GSM1315621’, ‘GSM1315637’, ‘GSM1315754’, ‘GSM1315768’

# removing duplicate samples
metadata <- metadata[!duplicated(metadata$Sample.Name), ]

# setting row names
rownames(metadata) <- metadata$Sample.Name

# re-verifying
dim(metadata)
all(colnames(counts) %in% rownames(metadata))
all(rownames(metadata) %in% colnames(counts))

# 3 extra samples in metadata- 3 samples were eliminated in raw count generation by NCBI

# subsetting metadata to only keep samples present in counts file

metadata <- metadata[colnames(counts), ]

# re-verifying 
dim(metadata)
all(colnames(counts) %in% rownames(metadata))
all(rownames(metadata) %in% colnames(counts))

table(metadata$tissue_type)
# lesional psoriatic skin - 90; normal skin - 81

# DESeq2 analysis
library(DESeq2)

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata,
                              design = ~ tissue_type)

# checking baseline
levels(dds$tissue_type)
# lesional psoriatic skin set as baseline
# changing the reference
dds$tissue_type <- relevel(dds$tissue_type, ref="normal skin")

# re-checking baseline
levels(dds$tissue_type)

# running DESeq 
dds <- DESeq(dds)

# extracting results
res <- results(dds)
summary(res)

# adjusting threshold to padj < 0.05
res_sig <- results(dds, alpha = 0.05)
summary(res_sig)

#adding a log2fc filter
res_sig<- results(dds, alpha = 0.05, lfcThreshold = 1)
summary(res_sig)

# creating a dataframe
res_df <- as.data.frame(res_sig)
res_df <- res_df[order(res_df$padj), ]
head(res_df) 

# converting gene IDs to gene symbols

library(org.Hs.eg.db)
library(AnnotationDbi)

res_df$symbol <- mapIds(org.Hs.eg.db,
                        keys = rownames(res_df),
                        column = "SYMBOL",
                        keytype = "ENTREZID",
                        multiVals = "first")
head(res_df[, c("symbol", "log2FoldChange", "padj")])

# PCA Visulization
vsd <- vst(dds, blind = TRUE)
PCA_plot <- plotPCA(vsd, intgroup = "tissue_type")
ggsave(filename="PCA.png", plot=PCA_plot, path="/Users/jriya/Desktop/RNA-seq_psoriasis/plots", width = 6, height = 4, units = "in")

# Volcano plot
library(EnhancedVolcano)
Volcano_plot <- EnhancedVolcano(res_df,
                                lab = res_df$symbol,
                                x = 'log2FoldChange',
                                y = 'padj',
                                title = 'Psoriatic legions vs Normal Skin gene expression',
                                pCutoff = 0.05,
                                FCcutoff = 1,
                                pointSize = 2,
                                labSize = 3)
ggsave(filename="Volcano.png", plot=Volcano_plot, path="/Users/jriya/Desktop/RNA-seq_psoriasis/plots", width = 6, height = 4, units = "in")

Volcano_plot

# Heatmap
library(pheatmap)
# top 50 significant genes 
top50sig <- head(res_df[!is.na(res_df$padj) & res_df$padj <0.05 & abs(res_df$log2FoldChange) > 1, ], 50)
plot_gene_heatmap(vsd, res_df, 
                  gene_ids = rownames(top50sig),
                  title = "Top 50 DEGs",
                  metadata = metadata,
                  cluster_rows = TRUE)

# dispersion plot .. model QC
dispersion_plot <- plotDispEsts(dds)

# Sample distance heatmap
sampleDistance_heatmap <- dist(t(assay(vsd)))
pheatmap(as.matrix(sampleDistance_heatmap),
         annotation_col = annotation_col,
         show_rownames = FALSE,
         show_colnames = FALSE,
         main = "Sample Distances")

# Oxidative stress analysis
oxidative_stress_genes <- c("SOD1","SOD2","SOD3","CAT","GPX1","GPX2","GPX4",
                            "TXN", "TXNRD1", "TXNRD2", "HMOX1", "HMOX2",
                            "NOX1", "NOX4", "CYBB","NFE2L2", "KEAP1",
                            "PRDX1", "PRDX2", "PRDX3", "NOS1", "NOS2", "NOS3",
                            "MPO", "CYBA", "NCF1", "NCF2", "NCF4",
                            "SIRT1", "SESN1", "SESN2", "MTOR", "CYCS", "PON1")

ox_table <- get_oxidative_stress_results(res_df, oxidative_stress_genes)
ox_table

write.csv(ox_table, "/Users/jriya/Desktop/RNA-seq_psoriasis/oxidative_stress_results.csv", 
          row.names=FALSE)

# gene-specfic heatmap of significant oxidative stress related genes found

sig_ox_genes <- rownames(res_df[res_df$symbol %in% ox_table$Gene, ])
# Oxidative stress heatmap
plot_gene_heatmap(vsd, res_df,
                  gene_ids = sig_ox_genes,
                  title = "Oxidative Stress Genes",
                  metadata = metadata,
                  cluster_rows = FALSE)

dot_plot <- plot_oxidative_dotplot(ox_table, title = "Oxidative Stress Genes")
dot_plot

intersect(ox_table$Gene, top50sig$symbol)
# none in the top50

library(clusterProfiler)

# significant DE genes
sig_genes <- rownames(res_df[!is.na(res_df$padj) &
                               res_df$padj < 0.05 &
                               abs(res_df$log2FoldChange) > 1, ])
# GO enrichment with ORA (pathway analysis for significant genes)
# Biological Processes
ora_BP <- enrichGO(gene=sig_genes,
                   OrgDb=org.Hs.eg.db,
                   keyType="ENTREZID",
                   ont="BP",
                   pAdjustMethod="BH",
                   pvalueCutoff=0.05,
                   qvalueCutoff = 0.05)
summary(ora_BP)

# visualizing ora_BP
barplot(ora_BP,
        showCategory = 20,
        title = "Top 20 Enriched GO Biological Processes",
        x = "GeneRatio") +
  theme(axis.text.y = element_text(size = 8),
        plot.title = element_text(hjust = 0.5))

# molecular functions
ora_MF <- enrichGO(gene = sig_genes,
                   OrgDb = org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05)

# cellular component
ora_CC <- enrichGO(gene = sig_genes,
                   OrgDb = org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05)

barplot(ora_MF,
        showCategory = 20,
        title = "Top 20 Enriched GO Molecular Functions",
        x = "GeneRatio") +
  theme(axis.text.y = element_text(size = 8),
        plot.title = element_text(hjust = 0.5))

barplot(ora_CC,
        showCategory = 20,
        title = "Top 20 Enriched GO Cellular Components",
        x = "GeneRatio") +
  theme(axis.text.y = element_text(size = 8),
        plot.title = element_text(hjust = 0.5))

library(fgsea)
library(msigdbr)

# creating a ranked gene list using Wald Statistic from DESeq2 results

ranked_genes <- res_df$stat
names(ranked_genes) <- rownames(res_df)
ranked_genes <- ranked_genes[!is.na(ranked_genes)]
ranked_genes <- sort(ranked_genes, decreasing = TRUE)
head(ranked_genes)
tail(ranked_genes)

# loading hallmarks from MsigDB
hallmarks <- msigdbr(species = "Homo sapiens", category = "H")
head(hallmarks[, c("gs_name","entrez_gene")])

# converting to match gsea format
hallmark_list <- split(hallmarks$entrez_gene, hallmarks$gs_name)
names(ranked_genes) <- as.character(names(ranked_genes))

# Running fgsea
set.seed(42)
fgsea_results <- fgsea(pathways = hallmark_list,
                       stats = ranked_genes,
                       minSize = 15,
                       maxSize = 500)
head(fgsea_results[order(fgsea_results$padj), c("pathway","NES", "padj")])

#significant fgsea pathways
sum(fgsea_results$padj < 0.05) #30

# Upregulated in psoriasis i.e. positive NES
sig_up <- fgsea_results[fgsea_results$padj <0.05 & fgsea_results$NES > 0, ]
sig_up <- sig_up[order(sig_up$NES, decreasing = TRUE), ]
nrow(sig_up) #26

sig_down <- fgsea_results[fgsea_results$padj < 0.05 & fgsea_results$NES < 0, ]
sig_down <- sig_down[order(sig_down$NES), ]
nrow(sig_down) #4

print(sig_down[, c("pathway", "NES","padj")])































