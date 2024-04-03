library(tidyverse)
library(ggplot2)
library(Seurat)


# inspect the downloaded .RDS file
seu <- readRDS(file = "data/data_9_9_annotated_seurat_all_ut.rds")
View(seu[[]])
dim(seu)
DimPlot(seu) + NoAxes()
DimPlot(seu, group.by = "major_clusters") + NoAxes()
DimPlot(seu, group.by = "subcluster") + NoAxes()
DimPlot(seu, group.by = "subcluster.v2") + NoAxes()
FeaturePlot(seu,features = "PRKCI", cols = c("lightgrey","darkred"), label = F, repel = F) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "PRKCZ", cols = c("lightgrey","darkred"), label = F, repel = F) + NoAxes() + NoLegend()

#########################################################






DimPlot(seu, group.by = "cell_type") + NoAxes()
DimPlot(seu, group.by = "Columnar_clusters") + NoAxes()
DimPlot(seu, group.by = "tissue") + NoAxes()
DimPlot(seu, group.by = "Tissue_in_paper") + NoAxes()
DimPlot(seu, group.by = "disease") + NoAxes() # normal includes various tissues
DimPlot(seu, group.by = "donor_id") + NoAxes()
table(seu$orig.ident) # body of stomach includes 41003 cells
table(seu$Tissue_in_paper) # CAG, GIM, NAG and NGB consist body of stomach
Idents(seu) <- "Celltypes_global"

FeaturePlot(seu,features = "EPCAM", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()











seu[["RNA"]]$data <- NULL
seu[["RNA"]] <- split(seu[["RNA"]], f = seu$sample_ID)
seu <- NormalizeData(seu, verbose = FALSE)
seu <- FindVariableFeatures(seu, verbose = FALSE)
seu <- ScaleData(seu, verbose = FALSE)
seu <- RunPCA(seu, npcs = 30, verbose = FALSE)

# Dimplot with integration
seu <- IntegrateLayers(
  object = seu, method = scVIIntegration,
  new.reduction = "integrated.scvi",
  # reference = c(6, 7, 9, 28),
  conda_env = "/usr/local/Caskroom/miniforge/base/envs/scvi-env2",
  verbose = TRUE
)
seu <- FindNeighbors(seu, reduction = "integrated.scvi", dims = 1:30, verbose = FALSE)
seu <- FindClusters(seu, resolution = 1.2, cluster.name = "scvi_clusters", verbose = FALSE)
seu <- RunUMAP(seu, reduction = "integrated.scvi", dims = 1:30, reduction.name = "umap.scvi", verbose = FALSE)
# DimPlot(seu, reduction = "umap.scvi", group.by = c("sample_ID", "scvi_clusters"), combine = FALSE)
DimPlot(seu, label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
DimPlot(seu, group.by = "sample_ID") + NoAxes()

saveRDS(seu, file = "RDSfiles/seu_DF.RDS")

# Check which cluster is epithelial, immune or stroma ----

FeaturePlot(seu,features = "EPCAM", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "PTPRC", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "FN1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()

FeaturePlot(seu,features = "MUC5AC", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "TFF1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "MUC6", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "LIPF", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "PGA3", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "REG4", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "TFF3", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "CHGA", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()

FeaturePlot(seu,features = "RGS5", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "NOTCH3", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "LUM", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "DCN", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "PLVAP", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "ACKR1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()

FeaturePlot(seu,features = "CD68", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "CD3D", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "CD79A", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()

FeaturePlot(seu,features = "IL2RA", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "KLRD1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "MS4A1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "TNFRSF17", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "KIT", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()

FeaturePlot(seu,features = "LRG1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "CD38", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "ENG", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()

FeaturePlot(seu,features = "TRPM5", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()

# subset epithelial, stromal and immune clusters
seu_epi <- subset(seu, idents = c(11,15,16,24,28,29,31,32,34,38,45,50)) # c18 and c37 positive fro both immune and epithelial markers, mostly from sample 29 > omit
# seu_str <- subset(seu, idents = c(6,20,23,26,33,36,41,48,49))
seu_imm <- subset(seu, idents = c(0:5,7:9,12:14,17,21,22,25,27,30,35,39,40,42:44,46,47,51)) # c45 and c52 negative for all the markers > omit

saveRDS(seu_epi, file = "RDSfiles/seu_epi_DF.RDS")
# saveRDS(seu_str, file = "RDSfiles/seu_str_DF.RDS")
# saveRDS(seu_imm, file = "RDSfiles/seu_imm_DF.RDS")