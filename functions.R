# Functions

# Heatmap function
plot_gene_heatmap <- function(vsd, res_df, gene_ids, title = "Heatmap", metadata, cluster_rows = TRUE,
                              fontsize = 10) {
  counts <- assay(vsd)[gene_ids, ]
  rownames(counts) <- res_df[gene_ids, "symbol"]
  
  annotation_col <- data.frame(Condition = metadata$tissue_type)
  rownames(annotation_col) <- colnames(counts)
  
  pheatmap(counts,
           annotation_col = annotation_col,
           show_colnames = FALSE,
           scale = "row",
           main = title,
           fontsize_row = 7,
           cluster_rows = cluster_rows)
}

get_oxidative_stress_results <- function(res_df, gene_list) {
  ox_results <- res_df[res_df$symbol %in% gene_list, ]
  ox_results <- ox_results[order(ox_results$log2FoldChange, decreasing = TRUE), ]
  ox_results <- ox_results[!is.na(ox_results$padj) & ox_results$padj < 0.05, ] # filtering for significance
  
  ox_table <- ox_results[, c("symbol", "baseMean", "log2FoldChange", "lfcSE","padj")]
  ox_table$baseMean <- round(ox_table$baseMean, 2)
  ox_table$log2FoldChange <- round(ox_table$log2FoldChange, 3)
  ox_table$lfcSE <- round(ox_table$lfcSE, 3)
  ox_table$padj <- formatC(ox_table$padj, format = "e", digits = 3)
  
  colnames(ox_table) <- c("Gene","Base Mean", "Log2FC", "SE", "Adjusted P-value")
  return(ox_table)
}

plot_oxidative_dotplot <- function(ox_table, title = "Oxidative Stress Genes") {
  # converting back to numeric
  plot_data <- data.frame(
    Gene = ox_table$Gene,
    Log2FC = ox_table$Log2FC,
    padj = as.numeric(ox_table$`Adjusted P-value`)
  )
  
  ggplot(plot_data, aes(x= Log2FC, y = reorder(Gene, Log2FC))) +
    geom_point(aes(size = -log10(padj), color = Log2FC)) +
    scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    scale_size_continuous(name = "-log10(padj)") +
    labs(title = title,
         x = "Log2 Fold Change",
         y = "Gene",
         color = "Log2FC") +
    theme_bw() +
    theme(axis.text.y = element_text(size = 11),
          plot.title = element_text(hjust = 0.5))
}