# ==============================================================
# Final Integrated Single-cell pipeline (Seurat v5 compatible)
# Includes: QC filtering, normalization, SingleR annotation,
# Z-score computation, log2(Obs/Exp), UMAP & result saving
# ==============================================================

# ---- Libraries ----
suppressPackageStartupMessages({
  library(Seurat)
  library(SingleR)
  library(celldex)
  library(ggplot2)
  library(dplyr)
  library(scales)
})

# ---- Settings ----
setwd("/path/to/your/folder")   # <<< change this to your working directory
dir.create("Results", showWarnings = FALSE)

# ---- Read gene lists ----
up_genes <- read.delim("up_gene", header = FALSE, stringsAsFactors = FALSE)[,1]
down_genes <- read.delim("down_gene", header = FALSE, stringsAsFactors = FALSE)[,1]

# ---- List 10X files ----
matrix_files <- list.files(pattern = "matrix.mtx", full.names = TRUE)
feature_files <- list.files(pattern = "features.tsv", full.names = TRUE)
barcode_files <- list.files(pattern = "barcodes.tsv", full.names = TRUE)

# ---- Reference for SingleR ----
hpca_ref <- celldex::HumanPrimaryCellAtlasData()

# ---- Function to process each sample ----
process_sample <- function(matrix_file, feature_file, barcode_file, sample_name) {
  message("---- Processing: ", sample_name, " ----")

  # --- Read data ---
  data <- ReadMtx(mtx = matrix_file,
                  features = feature_file,
                  cells = barcode_file,
                  feature.column = 2)

  seu <- CreateSeuratObject(counts = data, project = sample_name, min.cells = 3, min.features = 200)

  # --- Calculate mitochondrial percentage ---
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

  # --- QC filtering ---
  seu <- subset(seu, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 10)

  # --- Normalization & scaling ---
  if (!"data" %in% Layers(seu[["RNA"]])) {
    seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 1e4)
  }
  seu <- FindVariableFeatures(seu)
  seu <- ScaleData(seu)

  # --- Dimensionality reduction & clustering ---
  seu <- RunPCA(seu, npcs = 20)
  seu <- RunUMAP(seu, dims = 1:10)
  seu <- FindNeighbors(seu, dims = 1:10)
  seu <- FindClusters(seu, resolution = 0.5)

  # --- SingleR annotation ---
  singler_res <- SingleR(test = GetAssayData(seu, layer = "data"),
                         ref = hpca_ref, labels = hpca_ref$label.main)
  seu$celltype <- singler_res$labels

  # --- Compute Z-scores ---
  expr <- as.matrix(GetAssayData(seu, layer = "data"))
  gene_means <- rowMeans(expr)
  gene_sds <- apply(expr, 1, sd)
  z_expr <- sweep(expr, 1, gene_means, "-")
  z_expr <- sweep(z_expr, 1, gene_sds, "/")

  write.csv(z_expr, file = paste0("Results/", sample_name, "_Zscores.csv"))

  # --- Obs/Exp computation for up/down genes ---
  celltypes <- unique(seu$celltype)
  obs_up <- exp_up <- obs_down <- exp_down <- setNames(rep(NA, length(celltypes)), celltypes)

  for (ct in celltypes) {
    cells <- WhichCells(seu, expression = celltype == ct)
    avg_z <- rowMeans(z_expr[, cells, drop = FALSE])
    expressed <- names(avg_z[avg_z > 0])

    obs_up[ct] <- sum(expressed %in% up_genes) / length(up_genes)
    obs_down[ct] <- sum(expressed %in% down_genes) / length(down_genes)
    exp_up[ct] <- exp_down[ct] <- length(expressed) / nrow(z_expr)
  }

  # --- Compute & normalize log2(Obs/Exp) ---
  log_ratio_up <- log2(obs_up / exp_up)
  log_ratio_down <- log2(obs_down / exp_down)
  log_ratio_up[!is.finite(log_ratio_up)] <- 0
  log_ratio_down[!is.finite(log_ratio_down)] <- 0

  log_ratio_up_norm <- rescale(log_ratio_up, to = c(-1, 1))
  log_ratio_down_norm <- rescale(log_ratio_down, to = c(-1, 1))

  df_up <- data.frame(CellType = names(log_ratio_up_norm),
                      Value = log_ratio_up_norm,
                      Direction = "Upregulated")
  df_down <- data.frame(CellType = names(log_ratio_down_norm),
                        Value = log_ratio_down_norm,
                        Direction = "Downregulated")
  df_combined <- rbind(df_up, df_down)

  write.csv(df_combined,
            paste0("Results/", sample_name, "_Normalized_LogObsExp.csv"),
            row.names = FALSE)

  # --- Plot Obs/Exp normalized ---
  p <- ggplot(df_combined, aes(x = CellType, y = Value, fill = Value)) +
    geom_col() +
    facet_wrap(~Direction, scales = "free_y") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    theme_bw(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
          axis.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold")) +
    labs(y = "Normalized log2(Obs/Exp)", x = "Cell Type",
         title = paste("Normalized log2(Obs/Exp) -", sample_name))

  ggsave(paste0("Results/", sample_name, "_ObsExpPlot.png"), p, width = 9, height = 5)

  # --- UMAP with annotated cell types ---
  p_umap <- DimPlot(seu, group.by = "celltype", label = TRUE, repel = TRUE) +
    ggtitle(paste(sample_name, " - UMAP (SingleR annotated)")) +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))

  ggsave(paste0("Results/", sample_name, "_UMAP_Annotated.png"), p_umap, width = 6, height = 5)

  message("✅ ", sample_name, " completed successfully.")
}

# ---- Run for all samples ----
for (i in seq_along(matrix_files)) {
  sample_name <- gsub("_matrix.mtx", "", basename(matrix_files[i]))
  process_sample(matrix_files[i], feature_files[i], barcode_files[i], sample_name)
}

message("✅ All samples processed. Results saved in 'Results/' directory.")

