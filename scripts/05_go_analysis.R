
# This script is used to perform Gene Ontology (GO) analysis on the top 500 genes associated with PC2 in the hMSC dataset.

# --- START OF USER CONFIGURATION ---
# Root directory for all analysis outputs.
dir_output <- "/path/to/your/project/output"
# --- END OF USER CONFIGURATION ---

rm(list = ls(all = TRUE))
ls()

library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)

setwd(dir_output)
go_dir <- file.path(dir_output, "04_GO")
if (!dir.exists(go_dir)) dir.create(go_dir)

##### Setup the Seurat Object ####

load(file = "03_CorHeat/hMSC.CorHeat.Rdata")

##### Define gene select function and get genes ####

getGeneFromCluster <- function(hc, cluster_k, cluster_gene) {
  clusters <- dendextend::cutree(hc, k = cluster_k)
  gene <- hc$labels[clusters == cluster_gene]
  return(gene)
}

gene.c1 <- getGeneFromCluster(hc, cluster_k = 4, cluster_gene = 1)
gene.c2 <- getGeneFromCluster(hc, cluster_k = 4, cluster_gene = 2)
gene.c3 <- getGeneFromCluster(hc, cluster_k = 4, cluster_gene = 3)
gene.c4 <- getGeneFromCluster(hc, cluster_k = 4, cluster_gene = 4)

gene_pick
# [1] "TNFRSF12A" "PVR"       "ATP1B3"    "TNFRSF1A"  "IFNGR1"   
match(gene_pick, gene.c1)
# [1] NA NA 93 52 73
match(gene_pick, gene.c3)
# [1] 53 19 NA NA NA

##### Generate GO information ####

gene.c1.GO <- enrichGO(gene.c1, OrgDb = org.Hs.eg.db, keyType = "SYMBOL") # ECM Block (3)
gene.c2.GO <- enrichGO(gene.c2, OrgDb = org.Hs.eg.db, keyType = "SYMBOL") # Microtubule Block (1)
gene.c3.GO <- enrichGO(gene.c3, OrgDb = org.Hs.eg.db, keyType = "SYMBOL") # Adhesion Block (2)
gene.c4.GO <- enrichGO(gene.c4, OrgDb = org.Hs.eg.db, keyType = "SYMBOL") # Actin Block (4)

write.table(gene.c1.GO@result[which(gene.c1.GO@result$qvalue <= 0.05), c(1, 2, 3, 4, 6, 7)], file = "04_GO/gene.c1.GO.txt", row.names = FALSE, col.names = TRUE, sep = "\t")
write.table(gene.c2.GO@result[which(gene.c2.GO@result$qvalue <= 0.05), c(1, 2, 3, 4, 6, 7)], file = "04_GO/gene.c2.GO.txt", row.names = FALSE, col.names = TRUE, sep = "\t")
write.table(gene.c3.GO@result[which(gene.c3.GO@result$qvalue <= 0.05), c(1, 2, 3, 4, 6, 7)], file = "04_GO/gene.c3.GO.txt", row.names = FALSE, col.names = TRUE, sep = "\t")
write.table(gene.c4.GO@result[which(gene.c4.GO@result$qvalue <= 0.05), c(1, 2, 3, 4, 6, 7)], file = "04_GO/gene.c4.GO.txt", row.names = FALSE, col.names = TRUE, sep = "\t")

##### Draw GO barplot ####

png(filename = "04_GO/Block1.GO.png", width = 15, height = 8, units = "cm", res = 600)
barplot(filter(gene.c2.GO, qvalue <= 0.05), showCategory = 10) +
  ggtitle("Gene Ontology Enrichment: Microtubule Block Genes (PC2 Top 500)") +
  theme_bw(base_family = "tnr") +
  theme(text = element_text(family = "tnr", size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10, lineheight = 0.7),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        plot.title = element_text(hjust = 1, size = 12),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1))
dev.off()

png(filename = "04_GO/Block2.GO.png", width = 15, height = 3, units = "cm", res = 600)
barplot(filter(gene.c3.GO, qvalue <= 0.05), showCategory = 1) +
  ggtitle("Gene Ontology Enrichment: Adhesion Block Genes (PC2 Top 500)") +
  theme_bw(base_family = "Times New Roman") +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y =element_text(size = 10),
        legend.key.height = unit(0.2, "cm"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        plot.title = element_text(hjust = 1, size = 12),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1))
dev.off()

png(filename = "04_GO/Block3.GO.png", width = 15, height = 8, units = "cm", res = 600, family = "Times")
barplot(filter(gene.c1.GO, qvalue <= 0.05), showCategory = 10) +
  ggtitle("Gene Ontology Enrichment: ECM Block Genes (PC2 Top 500)") +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10, lineheight = 0.7),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        plot.title = element_text(hjust = 1, size = 12),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1))
dev.off()

png(filename = "04_GO/Block3.GO_short.png", width = 15, height = 8, units = "cm", res = 600, family = "Times")
gene.c1.GO_short <- gene.c1.GO
gene.c1.GO_short@result <- gene.c1.GO@result[c(-7, -8, -11),]
barplot(gene.c1.GO_short, showCategory = 10) +
  ggtitle("Gene Ontology Enrichment: ECM Block Genes (PC2 Top 500)") +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10, lineheight = 0.7),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        plot.title = element_text(hjust = 1, size = 12),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1))
dev.off()

png(filename = "04_GO/Block4.GO.png", width = 15, height = 4, units = "cm", res = 600, family = "Times")
barplot(filter(gene.c4.GO, qvalue <= 0.05), showCategory = 2) +
  ggtitle("Gene Ontology Enrichment: Actin Block Genes (PC2 Top 500)") +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10, lineheight = 0.7),
        legend.key.height = unit(0.2, "cm"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        plot.title = element_text(hjust = 1, size = 12),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1))
dev.off()

##### Save genes ####

save(gene.c1, gene.c2, gene.c3, gene.c4, file = "04_GO/hMSC.GO.Rdata")
