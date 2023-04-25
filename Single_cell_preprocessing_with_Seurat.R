

library(dplyr)
library(Seurat)
library(patchwork)

nameplot <- "PDNew.png"
expression_matrix <- ReadMtx(
  mtx = "./data/PDNew/GSE178265_Homo_matrix.mtx.gz", features = "./data/PDNew/GSE178265_Homo_features.tsv.gz",
  cells = "./data/PDNew/GSE178265_Homo_bcd.tsv.gz"
)
pbmc <- CreateSeuratObject(counts = expression_matrix, project = "PDNew", min.cells = 3, min.features = 200)



# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")


# Show QC metrics for the first 5 cells
head(pbmc@meta.data, 5)

# Visualize QC metrics as a violin plot
png(file = paste("~/Research_Practical/Seurat/Results/Vln_",nameplot))
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()


# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
png(file = paste("~/Research_Practical/Seurat/Results/Feature_",nameplot))
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
dev.off()

pbmc <- subset(pbmc, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 5)

#Normalize data 
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
#The standard parameters are in this function this is not necessary -> pbmc <- NormalizeData(pbmc) This works too

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)


#SAVE
#plot variable features with and without labels
png(file = paste("~/Research_Practical/Seurat/Results/FeatureVar_",nameplot))
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2
dev.off()

#linear transformation to scale the data. Mean expression set to 0 and mean variance to 1
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

#PCA  input can be defined using features argument
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

#Seurat provides several useful ways of visualizing both cells and features that define the PCA, including VizDimReduction(), DimPlot(), and DimHeatmap()
# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

png(file = paste("~/Research_Practical/Seurat/Results/PCA_",nameplot),width = 1920, height = 1080, res = 120)
DimPlot(pbmc, reduction = "pca")
dev.off()

DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
png(file = paste("~/Research_Practical/Seurat/Results/HM_",nameplot),width = 1920, height = 1080, res = 120)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()

# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
# pbmc <- JackStraw(pbmc, num.replicate = 100)
# pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
# 
# #"significant PC will be stongly enriched with low p value like in here PC 1-6
# png(file = paste("C:/Users/user/Documents/RWorkspace/Research_Practical/Results/JackStraw_",nameplot))
# JackStrawPlot(pbmc, dims = 1:15)
# dev.off()
# 
# ElbowPlot(pbmc)
#The bend(elbow) indicates the cutoff of PC where most true signals are
#Here around 9

pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.6)

save.image(file = "./PDData.RData")

# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
#reticulate::py_install(packages = 'umap-learn')

pbmc <- RunUMAP(pbmc, dims = 1:10)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters


#SAVE
png(file = paste("~/Research_Practical/Seurat/Results/umap_",nameplot),width = 1920, height = 1080, res = 120)
DimPlot(pbmc, reduction = "umap")
dev.off()

#saveRDS(pbmc, file = "../output/pbmc_tutorial.rds")

# find all markers of cluster 2
# cluster2.markers <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25)
# head(cluster2.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
#cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
#head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

#Marker depending on sample region

#Substantia Nigra
png(file = paste("~/Research_Practical/Seurat/Results/VlnCluster_SN_",nameplot),width = 1920, height = 1080, res = 120)
VlnPlot(pbmc, features = c("GFAP", "OLR1", "GINS3", "TH", "SLC6A3", "RGS5", "GAD1", "GAD2", "CSF1R", "MOG", "MOBP", "PALM2", "LGALS1", "PPM1G", "VCAN"))
dev.off()



png(file = paste("~/Research_Practical/Seurat/Results/FeaturePlot_",nameplot),width = 1920, height = 1080, res = 120)
FeaturePlot(pbmc, features = c("GFAP", "OLR1", "GINS3", "TH", "SLC6A3", "RGS5", "GAD1", "GAD2", "CSF1R", "MOG", "MOBP", "PALM2", "LGALS1", "PPM1G", "VCAN"))
dev.off()


#DoHeatmap() generates an expression heatmap for given cells and features. 
#In this case, we are plotting the top 20 markers (or all markers if less than 20) for each cluster.

pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
png(file = paste("~/Research_Practical/Seurat/Results/HMMarker_",nameplot),width = 1920, height = 1080, res = 120)
DoHeatmap(pbmc, features = top10$gene) + NoLegend()
dev.off()

#Map markers to known cell type. Need to be adapted to fit the data

new.cluster.ids <- c("ODC", "ODC", "ODC", "AC", "Microglia", "ExNeuron",
                     "ODC", "OPC", "ODC", "ODC", "InhiNeuron", "AC", "Microglia", "?", "Endothelial", "DA", "ExNeuron", "ODC", "AC", "AC",
                     "AC", "OPC", "Microglia", "Microglia", "Endo/AC" , "Endo/AC", "InhiNeuron", "OPC", "AC", "AC")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)

png(file = paste("~/Research_Practical/Seurat/Results/umaplabel_",nameplot),width = 1920, height = 1080, res = 120)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()
save.image(file = "./PDData.RData")
# saveRDS(pbmc, file = "C:/Users/user/Documents/RWorkspace/Research_Practical/Agarwal_dataSeurat.rds")

# Separating the clusters into different variables

for (i in 0:length(pbmc$seurat_clusters)) {
  match <- !is.na(match(pbmc$seurat_clusters,as.character(i)))
  match <- grep(TRUE, match)
  Cluster <- pbmc[,match]
  RNAData <- Cluster@assays$RNA[,]
  TheMatrix2 <- as.matrix(RNAData)
  write.table(TheMatrix2, file = paste("./NewModel/Cluster_",as.character(i),".txt",sep = ""))
}

