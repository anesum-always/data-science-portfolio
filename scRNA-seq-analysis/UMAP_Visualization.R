#Load Libraries library(Seurat) 
library(ggplot2) 
library(cowplot)
library(dplyr)
# Load Seurat object 
cluster4_cells <- readRDS("C:/Users/Documents/Cluster 4/seurat_obj_mcherry-gfp-SCT_snn_res.0.4-cluster4.rds")
# Basic preprocessing 
cluster4_cells <- NormalizeData(cluster4_cells) 
cluster4_cells<- FindVariableFeatures(cluster4_cells) 
cluster4_cells <- ScaleData(cluster4_cells)
cluster4_cells <- RunPCA(cluster4_cells)
cluster4_cells <- FindNeighbors(cluster4_cells, dims = 1:10) 
cluster4_cells <- FindClusters(cluster4_cells, resolution = 0.5)
cluster4_cells <- RunUMAP(cluster4_cells, dims = 1:10)
# Set identity 
Idents(cluster4_cells) <- "seurat_clusters"
# Tag each cell with its group
gfp_ids <- colnames(cluster4_cells)[cluster4_cells$orig.ident == "GFP"]
mcherry_ids <- colnames(cluster4_cells)[cluster4_cells$orig.ident == "mCherry"]
cluster4_cells$group_highlight <- "Other"
cluster4_cells$group_highlight[gfp_ids] <- paste0("GFP_", FetchData(cluster4_cells, "seurat_clusters")[gfp_ids, 1]) 
cluster4_cells$group_highlight[mcherry_ids] <- paste0("mCherry_", FetchData(cluster4_cells, "seurat_clusters")[mcherry_ids, 1])
# Biological labels for clusters 
cluster_labels <- c(
"mCherry_1" = "Amacrine",
"GFP_1" = "Retinal Ganglion", 
) 
# Color Scheme for unique cell types 
cell_type_colors <- c(
"Amacrine" = "#d95f02", "Retinal Ganglion" = "#DC143C", 
)
# Extract UMAP coordinates 
umap_coords <- as.data.frame(Embeddings(cluster4_cells, "umap")) 
umap_coords$cell <- rownames(umap_coords)
colnames(umap_coords)[1:2] <- c("UMAP_1", "UMAP_2") 
# Get metadata
meta <- cluster4_cells@meta.data
meta$cell <- rownames(meta)
# GFP data gfp_meta <- meta[meta$orig.ident == "GFP", ] 
gfp_meta$group <- paste0("GFP_", gfp_meta$seurat_clusters) 
gfp_meta$label <- cluster_labels[gfp_meta$group] 
gfp_plot_data <- left_join(gfp_meta, umap_coords, by = "cell")
# mCherry data mcherry_meta <- meta[meta$orig.ident == "mCherry", ]
mcherry_meta$group <- paste0("mCherry_", mcherry_meta$seurat_clusters)
mcherry_meta$label <- cluster_labels[mcherry_meta$group]
mcherry_plot_data <- left_join(mcherry_meta, umap_coords, by = "cell") 
# Define light grey for background 
cells bg_cells <- meta
bg_cells$cell <- rownames(bg_cells) 
bg_plot_data <- left_join(bg_cells, umap_coords, by = "cell") 
# GFP plot 
p_gfp <- ggplot() +
geom_point(data = bg_plot_data, aes(x = UMAP_1, y = UMAP_2), color = "#505050", size = 1.5) +
geom_point(data = gfp_plot_data, aes(x = UMAP_1, y = UMAP_2, color = label), size = 1.5) + 
scale_color_manual(values = cell_type_colors) +
ggtitle("GFP+ Cells") +
DarkTheme() + 
theme( 
legend.position = "none",
plot.title = element_text(hjust = 0.5, color = "white"),
panel.grid = element_blank(),  # remove all grid lines
panel.grid.major = element_blank(),
panel.grid.minor = element_blank()
)
# mCherry plot 
p_mcherry <- ggplot() + 
geom_point(data = bg_plot_data, aes(x = UMAP_1, y = UMAP_2), color = "#505050", size = 1.5) + 
geom_point(data = mcherry_plot_data, aes(x = UMAP_1, y = UMAP_2, color = label), size = 1.5) +
scale_color_manual(values = cell_type_colors) +
ggtitle("mCherry+ Cells") + 
DarkTheme() + 
theme(
legend.position = "none",
plot.title = element_text(hjust = 0.5, color = "white"),
panel.grid = element_blank(), # remove all grid lines 
panel.grid.major = element_blank(), 
panel.grid.minor = element_blank()
) 
# Combine plots side by side
plot_grid(p_gfp, p_mcherry, ncol = 2)
