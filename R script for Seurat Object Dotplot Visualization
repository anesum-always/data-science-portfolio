#Libraries
library(Seurat)
library(ggplot2) 
library(biomaRt)
library(dplyr)
# Function to load and preprocess a cluster 
load_preprocess_cluster <- function(path, cluster_name) { 
obj <- readRDS(path)
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj) 
obj <- RunPCA(obj)
obj <- FindNeighbors(obj, dims = 1:10) 
obj <- FindClusters(obj, resolution = 0.5) 
obj <- RunUMAP(obj, dims = 1:10) 
obj$ParentCluster <- cluster_name 
return(obj) 
}
# Load and preprocess clusters
cluster3_cells <- load_preprocess_cluster("C:/Users/Documents/Cluster3/seurat_obj_mcherry-gfp-SCT_snn_res.0.4-cluster3.rds", "Cluster3") 
cluster4_cells <- load_preprocess_cluster("C:/Users/Documents/Cluster 4/seurat_obj_mcherry-gfp-SCT_snn_res.0.4-cluster4.rds", "Cluster4") 
cluster5_cells <- load_preprocess_cluster("C:/Users/Documents/Cluster 5/seurat_obj_mcherry-gfp-SCT_snn_res.0.4-cluster5.rds", "Cluster5") 
cluster8_cells <- load_preprocess_cluster("C:/Users/Documents/Cluster 8/seurat_obj_mcherry-gfp-SCT_snn_res.0.4-cluster8.rds", "Cluster8") 
cluster9_cells <- load_preprocess_cluster("C:/Users/Documents/Cluster 9/seurat_obj_mcherry-gfp-SCT_snn_res.0.4-cluster9.rds", "Cluster9")
# Merge all clusters
merged_clusters <- merge(cluster3_cells, y = list(cluster4_cells, cluster5_cells, cluster8_cells, cluster9_cells))
merged_clusters$Group <- paste0(merged_clusters$ParentCluster, "_", merged_clusters$sample) # use 'sample' column to separate GFP/mCherry
Idents(merged_clusters) <- "Group" 
#Define marker gene list by cell type 
marker_genes <- list( 
"Retinal Progenitors" = c("NOTCH1", "NOTCH3", "GLI2", "PAX6"),
"Bipolar" = c("VSX1", "GNB3", "CALB2B"),
# Flatten gene list and get unique gene symbols
all_genes <- unique(unlist(marker_genes))
# Map gene symbols to Ensembl IDs (zebrafish)
ensembl <- useEnsembl(biomart = "genes", dataset = "drerio_gene_ensembl") 
gene_mapping <- getBM(
attributes = c("ensembl_gene_id", "external_gene_name"),
filters = "external_gene_name", 
values = all_genes, 
mart = ensembl
)
#Build gene symbol → Ensembl ID map
symbol_to_ensembl <- setNames(gene_mapping$ensembl_gene_id, gene_mapping$external_gene_name)
#Filter and rename genes to Ensembl IDs 
valid_symbols <- intersect(all_genes, names(symbol_to_ensembl)) 
ensembl_genes <- symbol_to_ensembl[valid_symbols] 
# Only keep genes present in merged data
genes_present <- ensembl_genes[ensembl_genes %in% rownames(merged_clusters)] 
# Build Ensembl ID → Gene Symbol and → CellType maps 
ensembl_to_symbol <- setNames(names(ensembl_genes), ensembl_genes)
gene_to_type <- unlist(lapply(names(marker_genes), function(celltype) {
  genes <- marker_genes[[celltype]] 
  found <- genes[genes %in% names(symbol_to_ensembl)] 
  setNames(rep(celltype, length(found)), found) 
}))
ensembl_to_type <- gene_to_type[names(ensembl_genes)]
names(ensembl_to_type) <- ensembl_genes
# Generate DotPlot
dp <- DotPlot(merged_clusters, features = genes_present, group.by = "sample") 
plot_data <- dp$data 
#Annotate with gene symbol and cell type 
plot_data$Gene <- ensembl_to_symbol[as.character(plot_data$features.plot)] 
plot_data$CellType <- ensembl_to_type[as.character(plot_data$features.plot)] 
#Drop rows with missing info 
plot_data <- plot_data[!is.na(plot_data$Gene) & !is.na(plot_data$CellType), ] 
#Order axes 
gene_order <- rev(unique(unlist(marker_genes)))
plot_data$Gene <- factor(plot_data$Gene, levels = gene_order) 
plot_data$CellType <- factor(plot_data$CellType, levels = names(marker_genes))
#Use 'id' as sample identifier from DotPlot data
plot_data$sample <- plot_data$id 
plot_data$sample <- factor(plot_data$sample, levels = c("GFP_pos", "mCherry_pos")) 
#Plot with sample on x-axis and cell types on y-axis 
ggplot(plot_data, aes(x = sample, y = CellType)) + 
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) + 
  scale_color_gradient(low = "lightgrey", high = "blue") + 
  scale_size(range = c(1, 6)) + 
labs( 
  title = "Expression of Marker Genes by Sample (GFP vs mCherry)",
  x = "Sample",
  y = "Cell Type",
  size = "% Expressed",
  color = "Avg Scaled Expression" 
) + 
theme_minimal(base_size = 14) +
theme( 
  axis.text.x = element_text(angle = 45, hjust = 1, size = 12), 
  axis.text.y = element_text(size = 11),
  plot.title = element_text(hjust = 0.5) 
)
