
# This script is used to create correlation heatmaps for the top 500 genes associated with PC2 in the hMSC dataset.

# --- START OF USER CONFIGURATION ---
# Root directory for all analysis outputs.
dir_output <- "/path/to/your/project/output"
# --- END OF USER CONFIGURATION ---

rm(list = ls(all = TRUE))
ls()

library(Hmisc)
library(ComplexHeatmap)
library(dendextend)
library(colorRamp2)

setwd(dir_output)
cor_dir <- file.path(dir_output, "03_CorHeat")
if (!dir.exists(cor_dir)) dir.create(cor_dir)

##### Setup the Seurat Object ####

hMSC <- readRDS("01_SeuratQC/hMSC.rds")
load(file = "01_SeuratQC/hMSC.Seurat.Rdata")

gene_pc2.500 <- TopFeatures(object = hMSC[["pca"]], dim = 2, nfeatures = 500, balanced = FALSE)
gene_pick <- c("TNFRSF12A", "PVR", "ATP1B3", "TNFRSF1A", "IFNGR1")
match(gene_pick, gene_pc2.500)
# [1] 485 379 115  58  87

gene_path <- c("CFL1", "RHOA", "SSH2")
match(gene_path, gene_pc2.500)
# [1] NA NA NA

##### Create correlation table and cluster ####

cor_mat <- rcorr(t(as.matrix(hMSC.data.num[gene_pc2.500, ])), type = c("pearson"))
hc <- hclust(dist(cor_mat$r, method = "euclidean"), method = "complete")

cor_mat_path <- rcorr(t(as.matrix(hMSC.data.num[union(gene_path, gene_pc2.500), ])), type = c("pearson"))
hc_path <- hclust(dist(cor_mat_path$r, method = "euclidean"), method = "complete")

##### Define heatmap function ####

drawCorrHeatmap <- function(cor.mat, cluster_distance, cluster_method, gene.label, title = "Correlation Heatmap") {
  
  row_anno <- rowAnnotation(link = anno_mark(at = seq(1:dim(cor.mat)[1])[rownames(cor.mat) %in% gene.label],
                                             labels = rownames(cor.mat)[rownames(cor.mat) %in% gene.label],
                                             labels_gp = gpar(fontsize = 10)))
  
  range(cor_mat$r)
  # [1] -0.6372605  1.0000000
  
  Heatmap(cor.mat,
          col = colorRamp2(c(-0.65,-0.1, 0.45, 1),
                           c("#264385", "#DFDDE4", "#E23E35", "#B22A2A")),
          name = "Correlation",
          
          column_title = title,
          column_title_gp = gpar(fontsize = 12),
          
          clustering_distance_rows = cluster_distance,
          clustering_method_rows = cluster_method,
          clustering_distance_columns = cluster_distance,
          clustering_method_columns = cluster_method,
          
          show_row_names = FALSE,
          show_column_names = FALSE,
          
          right_annotation = row_anno,
          
          width = unit(6, "cm"),
          height = unit(6, "cm"),
          
          show_heatmap_legend = TRUE,
          heatmap_legend_param = list(title_gp = gpar(fontsize = 12),
                                      labels_gp = gpar(fontsize = 10)),
          
          use_raster = TRUE,
          raster_quality = 5)
}

##### Draw heatmap ####

png(filename = "03_CorHeat/Heatmap.png", width = 11, height = 8, units = "cm", res = 600, family = "Times")
drawCorrHeatmap(cor_mat$r,
                cluster_distance = "euclidean",
                cluster_method = "complete",
                gene.label = NULL,
                "PC2 Top 500 Genes Correlation Heatmap")
dev.off()

png(filename = "03_CorHeat/Heatmap_Pick.png", width = 13, height = 8, units = "cm", res = 600, family = "Times")
drawCorrHeatmap(cor_mat$r,
                cluster_distance = "euclidean",
                cluster_method = "complete",
                gene.label = gene_pick,
                "PC2 Top 500 Genes Correlation Heatmap")
dev.off()

png(filename = "03_CorHeat/Heatmap_Path.png", width = 12, height = 8, units = "cm", res = 600, family = "Times")
drawCorrHeatmap(cor_mat_path$r,
                cluster_distance = "euclidean",
                cluster_method = "complete",
                gene.label = gene_path,
                "PC2 Top 500 Genes Correlation Heatmap")
dev.off()

png(filename = "03_CorHeat/Heatmap_Path_Pick.png", width = 13, height = 8, units = "cm", res = 600, family = "Times")
drawCorrHeatmap(cor_mat_path$r,
                cluster_distance = "euclidean",
                cluster_method = "complete",
                gene.label = union(gene_path, gene_pick),
                "PC2 Top 500 Genes Correlation Heatmap")
dev.off()

##### Save Data ####

save(cor_mat, cor_mat_path, hc, hc_path, gene_path, gene_pick, gene_pc2.500, file = "03_CorHeat/hMSC.CorHeat.Rdata")
