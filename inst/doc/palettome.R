## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5,
  dpi = 96
)

## ----setup--------------------------------------------------------------------
library(palettome)

## ----example-data-------------------------------------------------------------
set.seed(42)

compartments <- list(
  T     = list(center = c(-5,  4), subclusters = c("T_naive", "T_mem", "T_reg", "T_cyto")),
  Mono  = list(center = c( 6,  2), subclusters = c("Mono_cl", "Mono_nc", "Mac")),
  B     = list(center = c( 0, -6), subclusters = c("B_naive", "B_mem", "Plasma"))
)

cells <- do.call(rbind, lapply(names(compartments), function(comp) {
  info <- compartments[[comp]]
  do.call(rbind, lapply(info$subclusters, function(cl) {
    offset <- stats::rnorm(2, sd = 1.6)
    n <- sample(150:400, 1)
    data.frame(
      UMAP_1  = stats::rnorm(n, info$center[1] + offset[1], 0.7),
      UMAP_2  = stats::rnorm(n, info$center[2] + offset[2], 0.7),
      cluster = cl,
      compartment = comp,
      stringsAsFactors = FALSE
    )
  }))
}))

str(cells)

## ----as-cluster-data----------------------------------------------------------
pdata <- as_cluster_data(
  cells,
  coord_cols  = c("UMAP_1", "UMAP_2"),
  cluster_col = "cluster"
)

pdata$coord_cols
pdata$has_family

## ----centroids----------------------------------------------------------------
compute_centroids(pdata)

## ----seurat, eval = FALSE-----------------------------------------------------
# # --- replaces the as_cluster_data() call above ---------------------------
# # pdata <- as_cluster_data(
# #   seurat_obj,
# #   reduction   = "umap",            # any reduction name (Seurat::Embeddings())
# #   dims        = 1:2,               # or 1:3 for a 3D embedding
# #   cluster_col = "seurat_clusters", # any column of seurat_obj[[]]
# #   family_col  = NULL               # optional: metadata column with an
# #                                    # explicit compartment/family label
# # )
# #
# # --- from here on, identical to the data-frame path ----------------------
# # fam     <- detect_families(pdata, k = 3)
# # session <- generate_palette(fam$assignment, mode = "harmonious",
# #                             neighbors = fam$neighbors, seed = 1)
# #
# # --- and hand the colors straight back to Seurat's own plotting ----------
# # Seurat::DimPlot(seurat_obj, reduction = "umap",
# #                 group.by = "seurat_clusters",
# #                 cols = cluster_colors(session))

## ----detect-families----------------------------------------------------------
fam <- detect_families(pdata, k = 3)
fam$method
fam$assignment

## ----dendrogram---------------------------------------------------------------
plot_dendrogram(fam$dendro, k = 3)

## ----recut--------------------------------------------------------------------
cut_families(fam$dendro, k = 5)

## ----explicit-families--------------------------------------------------------
pdata_explicit <- as_cluster_data(
  cells,
  coord_cols  = c("UMAP_1", "UMAP_2"),
  cluster_col = "cluster",
  family_col  = "compartment"
)
fam_explicit <- detect_families(pdata_explicit)
fam_explicit$method
head(fam_explicit$assignment)

## ----harmonious---------------------------------------------------------------
session <- generate_palette(
  fam$assignment,
  mode       = "harmonious",
  neighbors  = fam$neighbors,
  seed       = 1
)

cluster_colors(session)

## ----harmonious-plot----------------------------------------------------------
plot_palette_static(pdata, session)

## ----harmonious-swatches, fig.height = 4--------------------------------------
plot_swatches(session)

## ----contrast-----------------------------------------------------------------
session_contrast <- generate_palette(
  fam$assignment,
  mode      = "contrast",
  neighbors = fam$neighbors,
  seed      = 1
)

family_colors(session_contrast)

## ----contrast-swatches, fig.height = 4----------------------------------------
plot_swatches(session_contrast)

## ----seeds, fig.height = 4----------------------------------------------------
for (s in c(3, 11)) {
  plot_swatches(generate_palette(fam$assignment, mode = "harmonious",
                                 neighbors = fam$neighbors, seed = s))
}

## ----sweep-anchors, fig.height = 4--------------------------------------------
session_anchored <- generate_palette(
  fam$assignment,
  mode          = "harmonious",
  sweep_anchors = c("#0B3D91", "#FF5733", "#FFD166"),
  neighbors     = fam$neighbors
)
plot_swatches(session_anchored)

## ----per-family, fig.height = 4-----------------------------------------------
plot_swatches(generate_palette(fam$assignment, mode = "harmonious",
                               harmonious_style = "per_family",
                               neighbors = fam$neighbors, seed = 1))

## ----manual-------------------------------------------------------------------
session <- set_manual_color(session, cluster = "Plasma", color = "#2E5A87")
session$clusters[session$clusters$manual_color, c("cluster", "color", "manual_color")]

## ----manual-multi-------------------------------------------------------------
multi <- set_manual_color(session, cluster = c("B_naive", "B_mem", "Not_here"),
                          color = c("#6E6E6E", "#9A9A9A", "#000000"))
cluster_colors(multi)[c("B_naive", "B_mem")]

## ----regenerate---------------------------------------------------------------
session2 <- generate_palette(
  fam$assignment,
  session   = session,
  mode      = "contrast",
  neighbors = fam$neighbors
)
cluster_colors(session2)[["Plasma"]]

## ----manual-family------------------------------------------------------------
session2 <- set_manual_color(session2, family = "F1", color = "#7A4E9E")
session2 <- generate_palette(fam$assignment, session = session2,
                             mode = "harmonious", neighbors = fam$neighbors)

family_colors(session2)
cluster_colors(session2, format = "data.frame")

## ----optimize-----------------------------------------------------------------
optimize_color("#7B00FF", session = session2)

## ----reset--------------------------------------------------------------------
session3 <- reset_overrides(session2, cluster = "Plasma", family = FALSE)
any(session3$clusters$manual_color)

## ----distance-----------------------------------------------------------------
d <- color_distance(cluster_colors(session))
round(d[1:4, 1:4], 1)

## ----cvd, fig.height = 4------------------------------------------------------
simulate_cvd(cluster_colors(session)[1:4], type = "deutan")
plot_swatches(session, cvd = "deutan")

## ----use-colors, eval = FALSE-------------------------------------------------
# Seurat::DimPlot(seurat_obj, cols = cluster_colors(session))
# 
# ggplot2::ggplot(cells, ggplot2::aes(UMAP_1, UMAP_2, color = cluster)) +
#   ggplot2::geom_point() +
#   ggplot2::scale_color_manual(values = cluster_colors(session))

## ----color-hard-code, eval=FALSE----------------------------------------------
# cluster_colors(session) |> dput()

## ----colors-df----------------------------------------------------------------
head(cluster_colors(session, format = "data.frame"))

## ----json---------------------------------------------------------------------
path <- tempfile(fileext = ".json")
export_palette_json(session, path)

session_back <- import_palette_json(path)
identical(cluster_colors(session_back), cluster_colors(session))

## ----csv----------------------------------------------------------------------
export_palette_csv(session, clusters_file = tempfile(fileext = ".csv"))

## ----downsample---------------------------------------------------------------
preview <- downsample_stratified(pdata, max_n = 500, min_per_cluster = 20)
nrow(preview$cells)
table(preview$cells$cluster)

## ----ui, eval = FALSE---------------------------------------------------------
# session <- launch_palettome_ui(pdata, session = session)

