
# This script is used to generate a heatmap of gene expression for selected genes in the hMSC dataset.

# --- START OF USER CONFIGURATION ---
# Root directory for all analysis outputs.
dir_output <- "/path/to/your/project/output"
# --- END OF USER CONFIGURATION ---

rm(list = ls(all = TRUE))
ls()

library(ComplexHeatmap)
library(dendextend)
library(colorRamp2)

setwd(dir_output)
cellheat_dir <- file.path(dir_output, "05_CellHeat")
if (!dir.exists(cellheat_dir)) dir.create(cellheat_dir)

##### Setup the Seurat Object ####
load(file = "01_SeuratQC/hMSC.Seurat.Rdata")
load(file = "04_GO/hMSC.GO.Rdata")

##### Setup the heatmap data metrix ####

cell_medium <- data.frame(cellID = colnames(hMSC.data.num)[order(colnames(hMSC.data.num))])
cell_medium$medium[grepl("1-d AM", cell_medium$cellID)] <- "1-d AM"
cell_medium$medium[grepl("1-d OM", cell_medium$cellID)] <- "1-d OM"
cell_medium$medium[grepl("GM", cell_medium$cellID)] <- "GM"

hc_gene.c1 <- hclust(dist(hMSC.data.num[gene.c1,], method = "euclidean"), method = "complete")
hc_gene.c3 <- hclust(dist(hMSC.data.num[gene.c3,], method = "euclidean"), method = "complete")

gene_cluster <- data.frame(geneName = union(gene.c1[hc_gene.c1$order], gene.c3[hc_gene.c3$order]))
gene_cluster$cluster[gene_cluster$geneName %in% gene.c3] <- "Adhesion Block"
gene_cluster$cluster[gene_cluster$geneName %in% gene.c1] <- "ECM Block"

hMSC.data.heat <- t(scale(t(hMSC.data.num[gene_cluster$geneName, cell_medium$cellID]), center = TRUE, scale = TRUE))
range(hMSC.data.heat)
# [1] -1.904077 14.128269

##### Define heatmap function ####

drawCellHeatmap <- function(hdata, anno_type, gene.label) {
  top_anno <- HeatmapAnnotation(Type = anno_type$medium,
                                col = list(Type = c("1-d AM" = "#E02B35",
                                                    "1-d OM" = "#2066A8",
                                                    "GM" = "#50AD9F")),
                                
                                annotation_legend_param = list(Type = list(title_gp = gpar(fontsize = 12),
                                                                           labels_gp = gpar(fontsize = 10))),
                                
                                show_annotation_name = FALSE,
                                simple_anno_size = unit(2, "mm"))
  
  row_anno <- rowAnnotation(link = anno_mark(at = seq(1:dim(hdata)[1])[rownames(hdata) %in% gene.label],
                                             labels = rownames(hdata)[rownames(hdata) %in% gene.label],
                                             side = "left",
                                             labels_gp = gpar(fontsize = 10)))

  hp_list <- Heatmap(hdata,
                     col = colorRamp2(c(-2, -1, 0, 1, 2),
                                      c("#264385", "#6A8CB2", "#DFDDE4", "#E23E35", "#B22A2A")),
                     name = "Expression",
                     column_title = "Expression Heatmap",
                     
                     column_title_gp = gpar(fontsize = 12),
                     
                     cluster_rows = TRUE,
                     show_row_dend = FALSE,
                     row_dend_reorder = TRUE,
                     cluster_columns = FALSE,
                     row_title = NULL,
                     show_row_names = FALSE,
                     show_column_names = FALSE,
                     top_annotation = top_anno,
                     left_annotation = row_anno,
                     width = unit(5, "cm"),
                     height = unit(5, "cm"),
                     
                     heatmap_legend_param = list(title_gp = gpar(fontsize = 12),
                                                 labels_gp = gpar(fontsize = 10))) +
    Heatmap(gene_cluster$cluster,
            name = "Block",
            show_column_names = FALSE,
            col = c("ECM Block" = "#F9AE78",
                    "Adhesion Block" = "#E47159"),
            width = unit(2, "mm"),
            
            heatmap_legend_param = list(title_gp = gpar(fontsize = 12),
                                        labels_gp = gpar(fontsize = 10)))
  
  draw(hp_list, row_km = 1, row_split = gene_cluster$cluster)
}

##### Define labled genes ####

gene_path <- c("CFL1", "RHOA", "SSH2")
gene_pick <- c("TNFRSF12A", "PVR", "ATP1B3", "TNFRSF1A", "IFNGR1")

##### Draw heatmap ####

png(filename = "05_CellHeat/CellHeatmap.png", width = 12, height = 7, units = "cm", res = 600, family = "Times")
drawCellHeatmap(hMSC.data.heat, cell_medium, gene.label = NULL)
dev.off()

png(filename = "05_CellHeat/CellHeatmap_Pick.png", width = 13, height = 7, units = "cm", res = 600, family = "Times")
drawCellHeatmap(hMSC.data.heat, cell_medium, gene_pick)
dev.off()

##### Save data ####

save(cell_medium, gene_cluster, hc_gene.c1, hc_gene.c3, hMSC.data.heat, file = "05_CellHeat/July2025.CellHeat.Rdata")
