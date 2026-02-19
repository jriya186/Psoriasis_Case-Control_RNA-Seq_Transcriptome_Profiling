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
# normalized top 50 counts
top50_counts <- assay(vsd)[rownames(top50sig), ]
rownames(top50_counts) <- top50sig$symbol
annotation_col <- data.frame(Condition = metadata$tissue_type)
rownames(annotation_col) <- colnames(top50_counts)

# plotting heatmap
heatmap_plot <- pheatmap(top50_counts,
                         annotation_col = annotation_col,
                         show_colnames = FALSE,
                         scale = "row",
                         main = "Top 50 DEGs",
                         fontsize_row = 7)

# dispersion plot .. model QC
dispersion_plot <- plotDispEsts(dds)

# Sample distance heatmap
sampleDistance_heatmap <- dist(t(assay(vsd)))
pheatmap(as.matrix(sampleDistance_heatmap),
         annotation_col = annotation_col,
         show_rownames = FALSE,
         show_colnames = FALSE,
         main = "Sample Distances")



















