# -------------------------
# Install required packages
# -------------------------
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("limma", "org.Hs.eg.db", "AnnotationDbi", "ggplot2"))

# -------------------------
# Load libraries
# -------------------------
library(limma)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ggplot2)

# -------------------------
# Step 1: Read datasets from local TSV files
# -------------------------
file_paths <- c(
  "C:/Users/2024/Downloads/code/GSE119290.tsv",
  "C:/Users/2024/Downloads/code/GSE53239.tsv",
  "C:/Users/2024/Downloads/code/GSE83601.tsv"
)

dataset_names <- c("GSE119290","GSE53239","GSE83601")
expr_list <- list()

for (i in seq_along(file_paths)) {
  expr <- read.table(file_paths[i], header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
  expr_list[[dataset_names[i]]] <- expr
}

# -------------------------
# Step 2: Process datasets individually
# -------------------------
deg_results <- list()

for (dataset in names(expr_list)) {
  
  expr <- expr_list[[dataset]]
  
  # -------------------------
  # Assign groups (Modify according to your dataset)
  # -------------------------
  group <- factor(c(rep("Tumor", ncol(expr)/2), rep("Normal", ncol(expr)/2)))
  design <- model.matrix(~ group)
  
  # -------------------------
  # Normalization (log2 if needed)
  # -------------------------
  if (max(expr) > 100) { # assume raw counts
    expr <- log2(expr + 1)
  }
  
  # -------------------------
  # Limma DEG analysis
  # -------------------------
  fit <- lmFit(expr, design)
  fit <- eBayes(fit)
  
  deg <- topTable(fit, adjust.method = "BH", number = Inf)
  
  # Filter DEGs
  deg_sig <- deg[deg$adj.P.Val < 0.05 & abs(deg$logFC) > 1, ]
  
  # -------------------------
  # Map gene symbols
  # -------------------------
  gene_ids <- rownames(deg_sig)
  
  # Check if rownames are Entrez IDs
  if (all(grepl("^[0-9]+$", gene_ids))) {
    symbols <- mapIds(org.Hs.eg.db,
                      keys = gene_ids,
                      column = "SYMBOL",
                      keytype = "ENTREZID",
                      multiVals = "first")
  } else {
    symbols <- gene_ids
  }
  
  # -------------------------
  # Keep only required columns
  # -------------------------
  deg_final <- data.frame(
    padj = deg_sig$adj.P.Val,
    pvalue = deg_sig$P.Value,
    log2FoldChange = deg_sig$logFC,
    Symbol = symbols
  )
  
  deg_results[[dataset]] <- deg_final
  
  # -------------------------
  # Save simplified CSV
  # -------------------------
  write.csv(deg_final, paste0(dataset, "_DEGs_simplified.csv"), row.names = TRUE)
  
  # -------------------------
  # Volcano plot
  # -------------------------
  deg_plot <- deg
  deg_plot$Significant <- ifelse(deg$adj.P.Val < 0.05 & abs(deg$logFC) > 1, "Yes", "No")
  
  volcano <- ggplot(deg_plot, aes(x = logFC, y = -log10(adj.P.Val), color = Significant)) +
    geom_point(alpha = 0.6) +
    scale_color_manual(values = c("No" = "grey", "Yes" = "red")) +
    theme_minimal() +
    labs(title = paste0("Volcano Plot: ", dataset),
         x = "log2 Fold Change",
         y = "-log10(adjusted P-value)") +
    geom_text(aes(label=ifelse(Significant=="Yes", rownames(deg_plot), "")),
              size=2, vjust=1.2, check_overlap = TRUE)
  
  ggsave(paste0(dataset, "_VolcanoPlot.png"), volcano, width = 8, height = 6)
}

# -------------------------
# Step 3: Identify common DEGs across all datasets
# -------------------------
common_genes <- Reduce(intersect, lapply(deg_results, function(x) x$Symbol))
write.csv(common_genes, "Common_DEGs_across_datasets.csv", row.names = FALSE)

cat("All datasets processed, DEGs saved, volcano plots generated, and common DEGs identified.\n")
