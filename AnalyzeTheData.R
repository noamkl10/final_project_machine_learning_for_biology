############################################################
# Cell-cell communication analysis pipeline
# Using CellCall + Seurat + visualization (Bubble, Circos, Heatmap, Sankey, UMAP)
# Author: <Your Name>
# Date: <Insert Date>
############################################################

##############################
# 1. Load Required Libraries #
##############################

# List of packages (general R, Bioconductor, and visualization packages)
packages <- c(
  "AnnotationDbi", "assertthat", "backports", "Biobase", "BiocGenerics", "stringr", "BiocManager",
  "BiocParallel", "bit", "bit64", "blob", "callr", "circlize", "cli", "clue", "cluster",
  "clusterProfiler", "colorspace", "ComplexHeatmap", "cowplot", "crayon", "data.table", "DBI", "desc", "devtools",
  "digest", "DO.db", "DOSE", "dplyr", "ellipsis", "enrichplot", "europepmc", "fansi", "farver",
  "fastmatch", "fgsea", "fs", "generics", "GetoptLong", "ggalluvial", "ggforce", "ggplot2",
  "ggplotify", "ggraph", "ggrepel", "ggridges", "GlobalOptions", "glue", "GO.db", "GOSemSim", "graphlayouts",
  "gridBase", "gridExtra", "gridGraphics", "gtable", "hms", "htmltools", "htmlwidgets", "httr", "igraph",
  "IRanges", "jsonlite", "knitr", "labeling", "lattice", "lifecycle", "magrittr", "MASS", "Matrix",
  "memoise", "mnormt", "munsell", "networkD3", "nlme", "pheatmap", "pillar", "pkgbuild", "pkgconfig",
  "pkgload", "plyr", "png", "polyclip", "prettyunits", "processx", "progress", "ps", "psych",
  "purrr", "qvalue", "R6", "RColorBrewer", "Rcpp", "remotes", "reshape2", "rjson", "rlang",
  "roxygen2", "rprojroot", "RSQLite", "rstudioapi", "rvcheck", "S4Vectors", "scales", "sessioninfo", "shape",
  "stringi", "testthat", "tibble", "tidygraph", "tidyr", "tidyselect", "triebeard", "tweenr", "urltools",
  "usethis", "vctrs", "viridis", "viridisLite", "withr", "xfun", "xml2", "readr"
)

# Load CellCall library (required for main analysis)
library(cellcall)

# Load all other packages
lapply(packages, library, character.only = TRUE)


###################################
# 2. Load and Prepare Input Matrix #
###################################

# Load input expression matrix (genes x cells/samples)
in.content <- read_tsv("final_matrix.tsv")

# Extract gene names (first column)
gene_names <- as.character(in.content[[1]])

# Remove gene name column to keep only numeric values
in.content <- in.content[,-1]

# Convert to data.frame and set row names to gene names
in.content <- as.data.frame(in.content)
rownames(in.content) <- gene_names

# Quick check of dimensions and preview
dim(in.content)
in.content[1:4, 1:4]


#######################################
# 3. Create CellCall NichCon Object   #
#######################################

mt <- CreateNichConObject(
  data = in.content,
  min.feature = 3,
  names.field = 1,          # Column split index for cell names
  names.delim = "_",        # Delimiter in column names
  source = "TPM",           # Set "counts" if raw counts
  scale.factor = 1e6,
  Org = "Homo sapiens",
  project = "YourProject"
)

# Run cell-cell communication inference
mt <- TransCommuProfile(
  object = mt,
  pValueCor = 0.05,
  CorValue = 0.1,
  topTargetCor = 1,
  p.adjust = 0.05,
  use.type = "median",
  probs = 0.9,
  method = "weighted",
  IS_core = TRUE,
  Org = 'Homo sapiens'
)

# Save intermediate results
save(mt, file = "mt_after_transcommu.RData")

#######################
# 8. PCA & UMAP       #
#######################

# Load expression counts
counts_matrix <- read.table("final_matrix.tsv", header = TRUE, row.names = 1, sep = "\t")
counts_matrix <- as.matrix(counts_matrix)

# Create Seurat object
library(Seurat)
seurat_obj <- CreateSeuratObject(
  counts = counts_matrix,
  project = "MyProject",
  min.cells = 3,
  min.features = 200
)

# Normalize, select features, scale, PCA
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))

# Run UMAP on first PCs
seurat_obj <- RunUMAP(seurat_obj, dims = 1:9)

# Define colors for cell types
my.colors <- c(
  "D"="#5e2a84", "Diff.ed.SPG"="pink","Diff.ing.SPG"="#1f78b4","L1"="#33a02c",
  "L2"="#bdbdbd","L3"="#fb9a99","P"="#cab2d6","S1"="#a6cee3","S2"="#40E0D0",
  "S3"="#b15928","S4"="gray","SPC7"="#b2df8a","SSC"="#e31a1c","ST"="#525252","Z"="#ff7f00"
)

# Assign cell type metadata
cell_types <- sub("_.*", "", colnames(seurat_obj))
seurat_obj <- AddMetaData(seurat_obj, metadata = cell_types, col.name = "cell_type")

# Compute cluster centers for label placement
umap_data <- Embeddings(seurat_obj, "umap") %>% as.data.frame()
umap_data$cell_type <- seurat_obj$cell_type
colnames(umap_data)[1:2] <- c("UMAP_1", "UMAP_2")
centers <- do.call(rbind, lapply(unique(umap_data$cell_type), function(ct){
  subset_data <- umap_data[umap_data$cell_type == ct, ]
  data.frame(cell_type=ct, UMAP_1=mean(subset_data$UMAP_1)-1.5, UMAP_2=mean(subset_data$UMAP_2))
}))

# Plot UMAP with labels
p <- DimPlot(seurat_obj, reduction = "umap", group.by = "cell_type", cols = my.colors) +
  geom_text(data = centers, aes(x = UMAP_1, y = UMAP_2, label = cell_type), inherit.aes = FALSE, size = 5)

# Save UMAP plot
ggsave("umap_with_labels.png", plot = p, width = 8, height = 6, dpi = 300)


#######################
# 4. Bubble Plot      #
#######################

# Reload object if needed
load("mt_after_transcommu.RData")

# Extract ligand-receptor matrix
n <- mt@data$expr_l_r_log2_scale

# Select sender-target interactions starting with "ST"
st_interactions <- grep("^ST-", colnames(n), value = TRUE)
st_interactions <- st_interactions[!grepl("-ST$", st_interactions)]

# Calculate enriched pathways per interaction
pathway.hyper.list <- lapply(st_interactions, function(i){
  print(i)
  getHyperPathway(data = n, object = mt, cella_cellb = i, Org="Homo sapiens")
})
names(pathway.hyper.list) <- st_interactions

# Prepare data for bubble plot
myPub.df <- getForBubble(pathway.hyper.list, cella_cellb = st_interactions)

# Generate bubble plot
p <- plotBubble(myPub.df)

# Beautify plot appearance
if("ggplot" %in% class(p)) {
  p <- p + labs(x = "ST") +
    scale_x_discrete(labels = function(x) gsub("^ST-", "", x)) +
    theme(
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      plot.title = element_text(size = 18),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12),
      strip.text = element_text(size = 12)
    )
}

# Save bubble plot
ggsave("bubble_plot.png", plot = p, width = 10, height = 10)


#######################
# 5. Circle Plot      #
#######################

# Define cell type colors
cell_color <- data.frame(
  color=c('#f35d58','#e31a1c','#1f78b4','#e78ac3','#33a02c','#808080','#fdbf6f','#ff7f00',
          '#cab2d6','#6a3d9a','#b2df8a','#a6cee3','#1f78b4','#ff7f00','#cab2d6'),
  stringsAsFactors = FALSE
)
rownames(cell_color) <- c("ST","SSC","Diff.ing.SPG","Diff.ed.SPG","L1","L2","L3",
                          "Z","P","D","SPC7","S1","S2","S3","S4")

# Subset ST interactions
n_sub <- n[, st_interactions, drop = FALSE]

# Save circos plot as PNG
png("circos_plot.png", width = 14, height = 10, units = "in", res = 300)
ViewInterCircos(
  object = n_sub,
  font = 1,
  cellColor = cell_color,
  lrColor = c("#F16B6F", "#84B1ED"),
  arr.type = "big.arrow",
  arr.length = 0.04,
  trackhight1 = 0.05,
  slot = "expr_l_r_log2_scale",
  linkcolor.from.sender = TRUE,
  gap.degree = 1,
  order.vector = rownames(cell_color),
  trackhight2 = 0.032,
  track.margin2 = c(0.01,0.12),
  DIY = TRUE
)
dev.off()


#######################
# 6. Heatmap          #
#######################

# Filter for stronger signals
threshold <- 0.5
max_values <- apply(n_sub, 1, max, na.rm = TRUE)
filtered_rows <- names(max_values[max_values > threshold])
st_data_filtered <- n_sub[filtered_rows, , drop = FALSE]

# Transpose for heatmap
st_data_filtered_transposed <- t(st_data_filtered)

# Update object for heatmap
mt_filtered_transposed <- mt
mt_filtered_transposed@data$expr_l_r_log2_scale <- st_data_filtered_transposed

# Plot heatmap
png("heatmap_output.png", width = 12, height = 8, units = "in", res = 300)
viewPheatmap(
  object = mt_filtered_transposed, 
  slot = "expr_l_r_log2_scale", 
  show_rownames = TRUE,
  show_colnames = TRUE,
  treeheight_row = 10,
  treeheight_col = 0,
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  fontsize = 12,
  angle_col = "90",
  main = "ST ➡ Germ cells"
)
dev.off()


#######################
# 7. Sankey Plot      #
#######################

# Infer ligand → receptor → TF signaling
mt <- LR2TF(
  object = mt,
  sender_cell="ST",
  recevier_cell="SSC",
  slot="expr_l_r_log2_scale",
  org="Homo sapiens"
)

# Extract relations
tmp <- mt@reductions$sankey
tmp1 <- dplyr::filter(tmp, weight1 > 0)
tmp.df <- trans2tripleScore(tmp1)

# Filter by weight threshold
weight_threshold <- 0.4
tmp.df.filtered <- tmp.df %>%
  filter(value > weight_threshold) %>%
  mutate(
    Ligand   = sub("^sender:",   "", Ligand),
    Receptor = sub("^receiver:", "", Receptor),
    TF       = sub("^TF:",       "", TF)
  )

# Color palette for nodes
mycol.vector <- c('#d62728','#1F78B4','#33A02C','#ff8c00','#ffb347','#A6CEE3','#B2DF8A','#FB9A99',
                  '#90ee90','#32cd32','#228b22','#B15928','#FF1493','#00CED1','#1f77b4','#6A3D9A',
                  '#DC143C','#4169E1','#228B22','#FF4500')

elments.num <- length(unique(tmp.df.filtered$Ligand))
mycol.vector.list <- rep(mycol.vector, times=ceiling(elments.num/length(mycol.vector)))

# Save Sankey plot
png("sankey_plot.png", width = 21, height = 20, units = "in", res = 300)
sankey_graph(
  df = tmp.df.filtered, 
  axes = 1:3, 
  mycol = mycol.vector.list[1:elments.num], 
  font.size = 6, 
  boder.col = "white", 
  isGrandSon = TRUE,
  set_alpha = 0.8
)
dev.off()


