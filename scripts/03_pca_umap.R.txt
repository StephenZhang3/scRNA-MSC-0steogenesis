
# This script is used to perform PCA and UMAP analysis on the hMSC dataset using Seurat package.

# --- START OF USER CONFIGURATION ---
# Root directory for all analysis outputs.
dir_output <- "/path/to/your/project/output"
# --- END OF USER CONFIGURATION ---

rm(list = ls(all = TRUE))
ls()

library(Seurat)
library(ggpubr)
library(ggplot2)

setwd(dir_output)
pca_dir <- file.path(dir_output, "02_PCA")
if (!dir.exists(pca_dir)) dir.create(pca_dir)

# Setup the Seurat Object
hMSC <- readRDS("01_SeuratQC/hMSC.rds")

my_col <- c("#FF0000", "#55A0FB", "#50AD9F")
my_comparisons <- list(c("1-d AM", "1-d OM"))

##### Draw UMAP and PCA ####

d_umap <- DimPlot(hMSC,
                  cols = my_col,
                  label = FALSE, 
                  group.by = "medium",
                  pt.size = 1) +
  ggtitle("Culture Medium") +
  labs(color = "Cell Type")

d_pca <- DimPlot(hMSC,
                 dims = c(1, 2),
                 cols = my_col,
                 pt.size = 1,
                 reduction = "pca",
                 group.by = "medium") +
  ggtitle("Culture Medium") +
  labs(color = "Cell Type")


png(filename = "02_PCA/UMAP.png", width = 6.5, height = 6.5, units = "cm", res = 600, family = "Times")
d_umap + theme(plot.title = element_blank(),
               axis.title = element_text(size = 12),
               axis.text = element_text(size = 10),
               legend.position = "none",
               axis.line = element_line(linewidth = 0.5),
               axis.ticks = element_line(linewidth = 0.5))
dev.off()

png(filename = "02_PCA/PCA.png", width = 6.5, height = 6.5, units = "cm", res = 600, family = "Times")
d_pca + theme(plot.title = element_blank(),
              axis.title = element_text(size = 12),
              axis.text = element_text(size = 10),
              legend.position = "none",
              axis.line = element_line(linewidth = 0.5),
              axis.ticks = element_line(linewidth = 0.5))
dev.off()

png(filename = "02_PCA/PCA_Legend.png", width = 8.5, height = 6.5, units = "cm", res = 600, family = "Times")
d_pca + theme(plot.title = element_blank(),
              axis.title = element_text(size = 12),
              axis.text = element_text(size = 10),
              legend.title = element_text(size = 12),
              legend.text = element_text(size = 10),
              axis.line = element_line(linewidth = 0.5),
              axis.ticks = element_line(linewidth = 0.5))
dev.off()

##### Draw Violin PC Plot ####

v_pc1 <- VlnPlot(hMSC, cols = my_col, features = "PC_1", group.by = "medium", y.max = 50) + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") + xlab("Cell Type") + labs(fill = "Cell Type")
v_pc2 <- VlnPlot(hMSC, cols = my_col, features = "PC_2", group.by = "medium", y.max = 35) + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") + xlab("Cell Type") + labs(fill = "Cell Type")
v_pc3 <- VlnPlot(hMSC, cols = my_col, features = "PC_3", group.by = "medium", y.max = 30) + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") + xlab("Cell Type") + labs(fill = "Cell Type")
v_pc4 <- VlnPlot(hMSC, cols = my_col, features = "PC_4", group.by = "medium", y.max = 25) + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") + xlab("Cell Type") + labs(fill = "Cell Type")

png(filename = "02_PCA/Vln_PC1.png", width = 6.5, height = 6.5, units = "cm", res = 600, family = "Times")
v_pc1 + theme(plot.title = element_text(size = 12, face = "plain"),
              axis.title.x = element_text(size = 12),
              axis.text.x = element_blank(),
              axis.line = element_line(linewidth = 0.5),
              axis.ticks = element_line(linewidth = 0.5),
              legend.position = "none")
dev.off()

png(filename = "02_PCA/Vln_PC2.png", width = 6.5, height = 6, units = "cm", res = 600, family = "Times")
v_pc2 + theme(plot.title = element_text(size = 12, face = "plain"),
              axis.title.x = element_text(size = 12),
              axis.text.x = element_blank(),
              axis.line = element_line(linewidth = 0.5),
              axis.ticks = element_line(linewidth = 0.5),
              legend.position = "none")
dev.off()

png(filename = "02_PCA/Vln_PC3.png", width = 6.5, height = 6, units = "cm", res = 600, family = "Times")
v_pc3 + theme(plot.title = element_text(size = 12, face = "plain"),
              axis.title.x = element_text(size = 12),
              axis.text.x = element_blank(),
              axis.line = element_line(linewidth = 0.5),
              axis.ticks = element_line(linewidth = 0.5),
              legend.position = "none")
dev.off()

png(filename = "02_PCA/Vln_PC4.png", width = 6.5, height = 6, units = "cm", res = 600, family = "Times")
v_pc4 + theme(plot.title = element_text(size = 12, face = "plain"),
              axis.title.x = element_text(size = 12),
              axis.text.x = element_blank(),
              axis.line = element_line(linewidth = 0.5),
              axis.ticks = element_line(linewidth = 0.5),
              legend.position = "none")
dev.off()

png(filename = "02_PCA/Vln_PC4_Legend.png", width = 8.5, height = 6, units = "cm", res = 600, family = "Times")
v_pc4 + theme(plot.title = element_text(size = 12, face = "plain"),
              axis.title.x = element_text(size = 12),
              axis.text.x = element_blank(),
              axis.line = element_line(linewidth = 0.5),
              axis.ticks = element_line(linewidth = 0.5),
              legend.title = element_text(size = 12),
              legend.text = element_text(size = 10))
dev.off()

##### Draw Violin Path Plot ####

v_CFL1 <- VlnPlot(hMSC, cols = my_col, features = "CFL1", group.by = "medium", y.max = 8) + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") + xlab("Cell Type") + labs(fill = "Cell Type")
v_RHOA <- VlnPlot(hMSC, cols = my_col, features = "RHOA", group.by = "medium", y.max = 8) + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") + xlab("Cell Type") + labs(fill = "Cell Type")
v_SSH2 <- VlnPlot(hMSC, cols = my_col, features = "SSH2", group.by = "medium", y.max = 5.5) + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") + xlab("Cell Type") + labs(fill = "Cell Type")

png(filename = "02_PCA/Vln_Path_CFL1.png", width = 6, height = 5.5, units = "cm", res = 600, family = "Times")
v_CFL1 + theme(plot.title = element_text(size = 12, face = "plain"),
               axis.title.x = element_text(size = 12),
               axis.text.x = element_blank(),
               axis.line = element_line(linewidth = 0.5),
               axis.ticks = element_line(linewidth = 0.5),
               legend.position = "none")
dev.off()

png(filename = "02_PCA/Vln_Path_RHOA.png", width = 6, height = 5.5, units = "cm", res = 600, family = "Times")
v_RHOA + theme(plot.title = element_text(size = 12, face = "plain"),
               axis.title.x = element_text(size = 12),
               axis.text.x = element_blank(),
               axis.line = element_line(linewidth = 0.5),
               axis.ticks = element_line(linewidth = 0.5),
               legend.position = "none")
dev.off()

png(filename = "02_PCA/Vln_Path_SSH2.png", width = 6, height = 5.5, units = "cm", res = 600, family = "Times")
v_SSH2 + theme(plot.title = element_text(size = 12, face = "plain"),
               axis.title.x = element_text(size = 12),
               axis.text.x = element_blank(),
               axis.line = element_line(linewidth = 0.5),
               axis.ticks = element_line(linewidth = 0.5),
               legend.position = "none")
dev.off()

##### Draw Violin 1d_OM Plot ####

v_TNFRSF12A <- VlnPlot(hMSC, cols = my_col, features = "TNFRSF12A", group.by = "medium", y.max = 8) + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") + xlab("Cell Type") + labs(fill = "Cell Type")
v_PVR <- VlnPlot(hMSC, cols = my_col, features = "PVR", group.by = "medium", y.max = 7) + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") + xlab("Cell Type") + labs(fill = "Cell Type")

png(filename = "02_PCA/Vln_OM_TNFRSF12A.png", width = 6, height = 5.5, units = "cm", res = 600, family = "Times")
v_TNFRSF12A + theme(plot.title = element_text(size = 12, face = "plain"),
                    axis.title.x = element_text(size = 12),
                    axis.text.x = element_blank(),
                    axis.line = element_line(linewidth = 0.5),
                    axis.ticks = element_line(linewidth = 0.5),
                    legend.position = "none")
dev.off()
png(filename = "02_PCA/Vln_OM_PVR.png", width = 6, height = 5.5, units = "cm", res = 600, family = "Times")
v_PVR + theme(plot.title = element_text(size = 12, face = "plain"),
              axis.title.x = element_text(size = 12),
              axis.text.x = element_blank(),
              axis.line = element_line(linewidth = 0.5),
              axis.ticks = element_line(linewidth = 0.5),
              legend.position = "none")
dev.off()

png(filename = "02_PCA/Vln_OM_PVR_Legend.png", width = 8, height = 5.5, units = "cm", res = 600, family = "Times")
v_PVR + theme(plot.title = element_text(size = 12, face = "plain"),
              axis.title.x = element_text(size = 12),
              axis.text.x = element_blank(),
              axis.line = element_line(linewidth = 0.5),
              axis.ticks = element_line(linewidth = 0.5),
              legend.title = element_text(size = 12),
              legend.text = element_text(size = 10))
dev.off()

##### Draw Violin 1d_AM Plot ####

v_ATP1B3 <- VlnPlot(hMSC, cols = my_col, features = "ATP1B3", group.by = "medium", y.max = 7.5) + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") + xlab("Cell Type") + labs(fill = "Cell Type")
v_TNFRSF1A <- VlnPlot(hMSC, cols = my_col, features = "TNFRSF1A", group.by = "medium", y.max = 7.5) + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") + xlab("Cell Type") + labs(fill = "Cell Type")
v_IFNGR1 <- VlnPlot(hMSC, cols = my_col, features = "IFNGR1", group.by = "medium", y.max = 7.5) + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") + xlab("Cell Type") + labs(fill = "Cell Type")

png(filename = "02_PCA/Vln_AM_ATP1B3.png", width = 6, height = 5.5, units = "cm", res = 600, family = "Times")
v_ATP1B3 + theme(plot.title = element_text(size = 12, face = "plain"),
              axis.title.x = element_text(size = 12),
              axis.text.x = element_blank(),
              axis.line = element_line(linewidth = 0.5),
              axis.ticks = element_line(linewidth = 0.5),
              legend.position = "none")
dev.off()

png(filename = "02_PCA/Vln_AM_TNFRSF1A.png", width = 6, height = 5.5, units = "cm", res = 600, family = "Times")
v_TNFRSF1A + theme(plot.title = element_text(size = 12, face = "plain"),
              axis.title.x = element_text(size = 12),
              axis.text.x = element_blank(),
              axis.line = element_line(linewidth = 0.5),
              axis.ticks = element_line(linewidth = 0.5),
              legend.position = "none")
dev.off()

png(filename = "02_PCA/Vln_AM_IFNGR1.png", width = 6, height = 5.5, units = "cm", res = 600, family = "Times")
v_IFNGR1 + theme(plot.title = element_text(size = 12, face = "plain"),
              axis.title.x = element_text(size = 12),
              axis.text.x = element_blank(),
              axis.line = element_line(linewidth = 0.5),
              axis.ticks = element_line(linewidth = 0.5),
              legend.position = "none")
dev.off()

##### Draw Violin MSC Positive Plot ####

v_ITGB1 <- VlnPlot(hMSC, cols = my_col, features = "ITGB1", group.by = "medium") + xlab("Cell Type") + labs(fill = "Cell Type")
v_CD44 <- VlnPlot(hMSC, cols = my_col, features = "CD44", group.by = "medium") + xlab("Cell Type") + labs(fill = "Cell Type")
v_NT5E <- VlnPlot(hMSC, cols = my_col, features = "NT5E", group.by = "medium") + xlab("Cell Type") + labs(fill = "Cell Type")
v_THY1 <- VlnPlot(hMSC, cols = my_col, features = "THY1", group.by = "medium") + xlab("Cell Type") + labs(fill = "Cell Type")
v_ENG <- VlnPlot(hMSC, cols = my_col, features = "ENG", group.by = "medium") + xlab("Cell Type") + labs(fill = "Cell Type")
v_ALCAM <- VlnPlot(hMSC, cols = my_col, features = "ALCAM", group.by = "medium") + xlab("Cell Type") + labs(fill = "Cell Type")

png(filename = "02_PCA/Vln_MSCposi_ITGB1.png", width = 6.5, height = 6, units = "cm", res = 600, family = "Times")
v_ITGB1 + theme(plot.title = element_text(size = 12, face = "plain"),
                axis.title.x = element_text(size = 12),
                axis.text.x = element_blank(),
                axis.line = element_line(linewidth = 0.5),
                axis.ticks = element_line(linewidth = 0.5),
                legend.position = "none")
dev.off()

png(filename = "02_PCA/Vln_MSCposi_CD44.png", width = 6.5, height = 6, units = "cm", res = 600, family = "Times")
v_CD44 + theme(plot.title = element_text(size = 12, face = "plain"),
               axis.title.x = element_text(size = 12),
               axis.text.x = element_blank(),
               axis.line = element_line(linewidth = 0.5),
               axis.ticks = element_line(linewidth = 0.5),
               legend.position = "none")
dev.off()

png(filename = "02_PCA/Vln_MSCposi_NT5E.png", width = 6.5, height = 6, units = "cm", res = 600, family = "Times")
v_NT5E + theme(plot.title = element_text(size = 12, face = "plain"),
               axis.title.x = element_text(size = 12),
               axis.text.x = element_blank(),
               axis.line = element_line(linewidth = 0.5),
               axis.ticks = element_line(linewidth = 0.5),
               legend.position = "none")
dev.off()

png(filename = "02_PCA/Vln_MSCposi_THY1.png", width = 6.5, height = 6, units = "cm", res = 600, family = "Times")
v_THY1 + theme(plot.title = element_text(size = 12, face = "plain"),
               axis.title.x = element_text(size = 12),
               axis.text.x = element_blank(),
               axis.line = element_line(linewidth = 0.5),
               axis.ticks = element_line(linewidth = 0.5),
               legend.position = "none")
dev.off()

png(filename = "02_PCA/Vln_MSCposi_ENG.png", width = 6.5, height = 6, units = "cm", res = 600, family = "Times")
v_ENG + theme(plot.title = element_text(size = 12, face = "plain"),
              axis.title.x = element_text(size = 12),
              axis.text.x = element_blank(),
              axis.line = element_line(linewidth = 0.5),
              axis.ticks = element_line(linewidth = 0.5),
              legend.position = "none")
dev.off()

png(filename = "02_PCA/Vln_MSCposi_ALCAM.png", width = 6.5, height = 6, units = "cm", res = 600, family = "Times")
v_ALCAM + theme(plot.title = element_text(size = 12, face = "plain"),
                axis.title.x = element_text(size = 12),
                axis.text.x = element_blank(),
                axis.line = element_line(linewidth = 0.5),
                axis.ticks = element_line(linewidth = 0.5),
                legend.position = "none")
dev.off()

##### Draw Violin MSC Negative Plot ####

v_CD14 <- VlnPlot(hMSC, cols = my_col, features = "CD14", group.by = "medium") + xlab("Cell Type") + labs(fill = "Cell Type")
v_CD19 <- VlnPlot(hMSC, cols = my_col, features = "CD19", group.by = "medium") + xlab("Cell Type") + labs(fill = "Cell Type")
v_CD34 <- VlnPlot(hMSC, cols = my_col, features = "CD34", group.by = "medium") + xlab("Cell Type") + labs(fill = "Cell Type")
v_PTPRC <- VlnPlot(hMSC, cols = my_col, features = "PTPRC", group.by = "medium") + xlab("Cell Type") + labs(fill = "Cell Type")
v_rna_CD79A <- VlnPlot(hMSC, cols = my_col, features = "rna_CD79A", group.by = "medium") + xlab("Cell Type") + labs(fill = "Cell Type")

png(filename = "02_PCA/Vln_MSCnega_CD14.png", width = 6.5, height = 6, units = "cm", res = 600, family = "Times")
v_CD14 + theme(plot.title = element_text(size = 12, face = "plain"),
               axis.title.x = element_text(size = 12),
               axis.text.x = element_blank(),
               axis.line = element_line(linewidth = 0.5),
               axis.ticks = element_line(linewidth = 0.5),
               legend.position = "none")
dev.off()

png(filename = "02_PCA/Vln_MSCnega_CD19.png", width = 6.5, height = 6, units = "cm", res = 600, family = "Times")
v_CD19 + theme(plot.title = element_text(size = 12, face = "plain"),
               axis.title.x = element_text(size = 12),
               axis.text.x = element_blank(),
               axis.line = element_line(linewidth = 0.5),
               axis.ticks = element_line(linewidth = 0.5),
               legend.position = "none")
dev.off()

png(filename = "02_PCA/Vln_MSCnega_CD34.png", width = 6.5, height = 6, units = "cm", res = 600, family = "Times")
v_CD34 + theme(plot.title = element_text(size = 12, face = "plain"),
               axis.title.x = element_text(size = 12),
               axis.text.x = element_blank(),
               axis.line = element_line(linewidth = 0.5),
               axis.ticks = element_line(linewidth = 0.5),
               legend.position = "none")
dev.off()

png(filename = "02_PCA/Vln_MSCnega_PTPRC.png", width = 6.5, height = 6, units = "cm", res = 600, family = "Times")
v_PTPRC + theme(plot.title = element_text(size = 12, face = "plain"),
                axis.title.x = element_text(size = 12),
                axis.text.x = element_blank(),
                axis.line = element_line(linewidth = 0.5),
                axis.ticks = element_line(linewidth = 0.5),
                legend.position = "none")
dev.off()

png(filename = "02_PCA/Vln_MSCnega_rna_CD79A.png", width = 6.5, height = 6, units = "cm", res = 600, family = "Times")
v_rna_CD79A + theme(plot.title = element_text(size = 12, face = "plain"),
                    axis.title.x = element_text(size = 12),
                    axis.text.x = element_blank(),
                    axis.line = element_line(linewidth = 0.5),
                    axis.ticks = element_line(linewidth = 0.5),
                    legend.position = "none")
dev.off()

##### Draw Gene Expression Dot Plot ####

png(filename = "02_PCA/Dot_Gene.png", width = 15, height = 6, units = "cm", res = 600, family = "Times")
DotPlot(hMSC,
        features = c("CFL1", "RHOA", "SSH2", "TNFRSF12A", "PVR", "ATP1B3", "TNFRSF1A", "IFNGR1"),
        # cols = c("#264385", "#6A8CB2", "#DFDDE4", "#E23E35", "#B22A2A"),
        group.by = "medium",
        dot.scale = 8) +
  scale_colour_gradient2(low = "#264385", mid = "#DFDDE4", high = "#B22A2A") +
  theme(plot.title = element_blank(),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 30, hjust = 1),
        axis.line = element_line(linewidth = 0.5),
        axis.ticks = element_line(linewidth = 0.5),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))
dev.off()

png(filename = "02_PCA/Dot_MSC.png", width = 15, height = 6, units = "cm", res = 600, family = "Times")
DotPlot(hMSC,
        features = c("ITGB1", "CD44", "NT5E", "THY1", "ENG", "ALCAM", "CD14", "CD19", "CD34", "PTPRC", "rna_CD79A"),
        # cols = c("#264385", "#6A8CB2", "#DFDDE4", "#E23E35", "#B22A2A"),
        group.by = "medium",
        dot.scale = 8) +
  scale_colour_gradient2(low = "#264385", mid = "#DFDDE4", high = "#B22A2A") +
  theme(plot.title = element_blank(),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 30, hjust = 1),
        axis.line = element_line(linewidth = 0.5),
        axis.ticks = element_line(linewidth = 0.5),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))
dev.off()
