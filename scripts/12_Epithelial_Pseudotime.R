############################################################
## GSE294300
## Figure 8A-B REDO
##
## CLEAN epithelial-cell Monocle3 pseudotime analysis
##
## Remove previous immune-like epithelial substates:
## Ep6 -> T/NK/ILC
## Ep7 -> Plasma
## Ep8 -> B
##
## Then rebuild Monocle3 FROM SCRATCH
##
## A: Pseudotime trajectory
## B: Normal / Tumor distribution on SAME trajectory
##
## Root label "1" is intentionally retained
############################################################

rm(list = ls())
options(stringsAsFactors = FALSE)
Sys.setenv(LANGUAGE = "en")
set.seed(12345)

############################################################
## 1. Required packages
############################################################

required_packages <- c(
  "Seurat",
  "SeuratObject",
  "monocle3",
  "Matrix",
  "dplyr",
  "SummarizedExperiment",
  "SingleCellExperiment",
  "igraph",
  "ggplot2",
  "patchwork"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop("缺少以下R包：", paste(missing_packages, collapse = ", "))
}

############################################################
## 2. Paths
############################################################

data_dir <- "C:/Users/92066/Desktop/单细胞孟德尔/Figure7/GSE294300"

input_file <- file.path(
  data_dir,
  "GSE294300_scRNA_final_annotated_with_tSNE.rds"
)

## OLD Monocle object:
## ONLY used to recover the already-validated Ep6/Ep7/Ep8 cell barcodes
old_monocle_dir <- "C:/Users/92066/Desktop/单细胞孟德尔/Figure8/GSE294300_Epithelial_Monocle3_pseudotime"

## New clean analysis directory
output_dir <- "C:/Users/92066/Desktop/单细胞孟德尔/Figure8/Figure8_AB_CLEAN"

new_monocle_dir <- file.path(
  output_dir,
  "GSE294300_CLEAN_Epithelial_Monocle3_pseudotime"
)

if (!file.exists(input_file)) {
  stop("找不到Seurat对象：", input_file)
}

if (!dir.exists(old_monocle_dir)) {
  stop("找不到旧Monocle3对象：", old_monocle_dir)
}

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

############################################################
## 3. Load corrected final Seurat object
############################################################

scRNA <- readRDS(input_file)

cat("\n====================================================\n")
cat("FINAL SEURAT OBJECT\n")
cat("====================================================\n")
cat("Cells:", ncol(scRNA), "\n")
cat("Genes:", nrow(scRNA), "\n")

required_meta <- c(
  "final_cluster",
  "celltype",
  "patient",
  "group"
)

missing_meta <- setdiff(
  required_meta,
  colnames(scRNA@meta.data)
)

if (length(missing_meta) > 0) {
  stop(
    "Seurat metadata缺少：",
    paste(missing_meta, collapse = ", ")
  )
}

if ("seurat_clusters" %in% colnames(scRNA@meta.data)) {
  stop(
    "当前对象仍包含seurat_clusters；请确认读取的是CORRECTED对象。"
  )
}

############################################################
## 4. Standardize metadata
############################################################

scRNA$final_cluster <- factor(
  as.character(scRNA$final_cluster),
  levels = as.character(0:21)
)

scRNA$celltype <- as.character(
  scRNA$celltype
)

scRNA$patient <- as.character(
  scRNA$patient
)

scRNA$group <- factor(
  as.character(scRNA$group),
  levels = c("Normal", "Tumor")
)

if (anyNA(scRNA$group)) {
  stop("group中出现Normal/Tumor以外的值。")
}

Seurat::Idents(scRNA) <- "final_cluster"

############################################################
## 5. Load OLD Monocle object
##
## IMPORTANT:
## Only retrieve Ep6/Ep7/Ep8 BARCODE identities.
## No old pseudotime / UMAP / trajectory will be reused.
############################################################

old_cds <- monocle3::load_monocle_objects(
  directory_path = old_monocle_dir
)

old_clusters <- monocle3::clusters(
  old_cds,
  reduction_method = "UMAP"
)

if (is.null(names(old_clusters))) {
  stop(
    "旧Monocle3 clusters()没有cell barcode names，停止分析。"
  )
}

if (!setequal(
  names(old_clusters),
  colnames(old_cds)
)) {
  stop(
    "旧Monocle cluster assignments与cell barcodes不完全对应。"
  )
}

old_clusters <- old_clusters[
  colnames(old_cds)
]

if (!identical(
  names(old_clusters),
  colnames(old_cds)
)) {
  stop("旧Monocle cluster重新排序失败。")
}

old_ep_state <- paste0(
  "Ep",
  as.character(old_clusters)
)

names(old_ep_state) <- names(old_clusters)

cat("\n====================================================\n")
cat("OLD Ep-state counts\n")
cat("====================================================\n")

print(
  table(old_ep_state)
)

############################################################
## 6. Identify Ep6 / Ep7 / Ep8 cells
############################################################

immune_like_states <- c(
  "Ep6",
  "Ep7",
  "Ep8"
)

remove_cells <- names(old_ep_state)[
  old_ep_state %in% immune_like_states
]

cat("\n====================================================\n")
cat("IMMUNE-LIKE CELLS TO REMOVE\n")
cat("====================================================\n")

print(
  table(
    old_ep_state[remove_cells]
  )
)

cat(
  "\nTotal cells to remove:",
  length(remove_cells),
  "\n"
)

############################################################
## 7. Strict barcode validation
############################################################

if (!all(
  remove_cells %in% colnames(scRNA)
)) {
  stop(
    "部分Ep6/Ep7/Ep8 cells无法映射回final Seurat对象。"
  )
}

original_epi_cells <- rownames(scRNA@meta.data)[
  scRNA$celltype == "Epithelial"
]

if (!all(
  remove_cells %in% original_epi_cells
)) {
  stop(
    "Ep6/Ep7/Ep8中存在不属于原始Epithelial compartment的细胞。"
  )
}

cat(
  "Original epithelial cells:",
  length(original_epi_cells),
  "\n"
)

############################################################
## 8. Expected cell-number validation
############################################################

expected_removed <- 1692L

if (length(remove_cells) != expected_removed) {
  stop(
    "预期删除1692个Ep6/Ep7/Ep8 cells，但实际为：",
    length(remove_cells)
  )
}

############################################################
## 9. Create refined major-cell annotation
############################################################

scRNA$celltype_refined <- scRNA$celltype

ep6_cells <- names(old_ep_state)[
  old_ep_state == "Ep6"
]

ep7_cells <- names(old_ep_state)[
  old_ep_state == "Ep7"
]

ep8_cells <- names(old_ep_state)[
  old_ep_state == "Ep8"
]

scRNA$celltype_refined[
  ep6_cells
] <- "T/NK/ILC"

scRNA$celltype_refined[
  ep7_cells
] <- "Plasma"

scRNA$celltype_refined[
  ep8_cells
] <- "B"

celltype_order <- c(
  "Epithelial",
  "T/NK/ILC",
  "B",
  "Plasma",
  "Mast",
  "Myeloid",
  "Endothelial",
  "Glial",
  "Mural cells",
  "Fibroblast"
)

scRNA$celltype_refined <- factor(
  scRNA$celltype_refined,
  levels = celltype_order
)

if (anyNA(scRNA$celltype_refined)) {
  stop("celltype_refined中出现NA。")
}

############################################################
## 10. Check before / after cell-type counts
############################################################

cat("\n====================================================\n")
cat("CELLTYPE BEFORE REFINEMENT\n")
cat("====================================================\n")

print(
  table(scRNA$celltype)
)

cat("\n====================================================\n")
cat("CELLTYPE AFTER REFINEMENT\n")
cat("====================================================\n")

print(
  table(scRNA$celltype_refined)
)

############################################################
## 11. Extract CLEAN epithelial cells
############################################################

clean_epi_cells <- rownames(scRNA@meta.data)[
  scRNA$celltype_refined == "Epithelial"
]

expected_clean <- 12559L

cat("\n====================================================\n")
cat("CLEAN EPITHELIAL SUBSET\n")
cat("====================================================\n")

cat(
  "Original epithelial:",
  length(original_epi_cells),
  "\n"
)

cat(
  "Removed Ep6/Ep7/Ep8:",
  length(remove_cells),
  "\n"
)

cat(
  "Clean epithelial:",
  length(clean_epi_cells),
  "\n"
)

if (length(clean_epi_cells) != expected_clean) {
  stop(
    "预期clean epithelial = 12559 cells，但实际为：",
    length(clean_epi_cells)
  )
}

############################################################
## 12. Save refinement information
############################################################

refinement_df <- data.frame(
  cell = names(old_ep_state),
  old_Ep_state = as.character(old_ep_state),
  remove_from_epithelial = old_ep_state %in% immune_like_states,
  final_cluster = as.character(
    scRNA$final_cluster[
      match(
        names(old_ep_state),
        colnames(scRNA)
      )
    ]
  ),
  old_celltype = as.character(
    scRNA$celltype[
      match(
        names(old_ep_state),
        colnames(scRNA)
      )
    ]
  ),
  refined_celltype = as.character(
    scRNA$celltype_refined[
      match(
        names(old_ep_state),
        colnames(scRNA)
      )
    ]
  ),
  stringsAsFactors = FALSE
)

write.csv(
  refinement_df,
  file.path(
    output_dir,
    "Epithelial_refinement_Ep6_Ep7_Ep8.csv"
  ),
  row.names = FALSE
)

############################################################
## 13. Save refined Seurat object
############################################################

saveRDS(
  scRNA,
  file.path(
    data_dir,
    "GSE294300_scRNA_final_REFINED_celltype.rds"
  ),
  compress = FALSE
)

############################################################
## 14. Subset clean epithelial cells
############################################################

epi <- subset(
  scRNA,
  cells = clean_epi_cells
)

cat(
  "\nClean epithelial object:",
  ncol(epi),
  "cells\n"
)

############################################################
## 15. Donor representation
############################################################

epi_meta <- as.data.frame(
  epi@meta.data
)

epi_donor_summary <- epi_meta |>
  dplyr::mutate(
    patient = as.character(patient),
    group = as.character(group)
  ) |>
  dplyr::group_by(
    patient,
    group
  ) |>
  dplyr::summarise(
    n_cells = dplyr::n(),
    .groups = "drop"
  ) |>
  dplyr::arrange(
    patient,
    group
  )

cat("\n====================================================\n")
cat("CLEAN epithelial cells per patient × group\n")
cat("====================================================\n")

print(
  epi_donor_summary,
  n = Inf
)

write.csv(
  epi_donor_summary,
  file.path(
    output_dir,
    "CLEAN_Epithelial_cells_per_patient_group.csv"
  ),
  row.names = FALSE
)

############################################################
## 16. RNA layer check
############################################################

rna_layers_before <- SeuratObject::Layers(
  epi[["RNA"]]
)

count_layers <- grep(
  "^counts",
  rna_layers_before,
  value = TRUE
)

if (length(count_layers) == 0) {
  stop("RNA assay中不存在counts layer。")
}

if (
  length(count_layers) > 1 ||
  !("counts" %in% rna_layers_before)
) {
  
  epi <- SeuratObject::JoinLayers(
    object = epi,
    assay = "RNA"
  )
}

rna_layers_after <- SeuratObject::Layers(
  epi[["RNA"]]
)

if (!("counts" %in% rna_layers_after)) {
  stop(
    "JoinLayers后不存在统一counts layer。"
  )
}

############################################################
## 17. Raw RNA counts
############################################################

counts_mat <- SeuratObject::LayerData(
  object = epi,
  assay = "RNA",
  layer = "counts"
)

if (!inherits(counts_mat, "Matrix")) {
  counts_mat <- Matrix::Matrix(
    counts_mat,
    sparse = TRUE
  )
}

if (!identical(
  colnames(counts_mat),
  colnames(epi)
)) {
  stop(
    "RNA counts细胞顺序与Seurat对象不一致。"
  )
}

if (anyDuplicated(
  colnames(counts_mat)
) > 0) {
  stop("存在重复cell barcode。")
}

if (anyDuplicated(
  rownames(counts_mat)
) > 0) {
  stop("存在重复gene name。")
}

cat(
  "\nCounts matrix:",
  nrow(counts_mat),
  "genes ×",
  ncol(counts_mat),
  "cells\n"
)

############################################################
## 18. Prepare cell metadata
############################################################

cell_metadata <- as.data.frame(
  epi@meta.data[
    colnames(counts_mat),
    ,
    drop = FALSE
  ]
)

cell_metadata$patient <- as.character(
  cell_metadata$patient
)

cell_metadata$group <- factor(
  as.character(cell_metadata$group),
  levels = c(
    "Normal",
    "Tumor"
  )
)

if (anyNA(cell_metadata$group)) {
  stop("clean epithelial metadata中group存在NA。")
}

############################################################
## 19. Gene metadata
############################################################

gene_metadata <- data.frame(
  gene_short_name = rownames(counts_mat),
  row.names = rownames(counts_mat),
  stringsAsFactors = FALSE
)

############################################################
## 20. Construct NEW Monocle3 CDS
############################################################

cds_epi <- monocle3::new_cell_data_set(
  expression_data = counts_mat,
  cell_metadata = cell_metadata,
  gene_metadata = gene_metadata
)

############################################################
## 21. Detect expressed genes
############################################################

cds_epi <- monocle3::detect_genes(
  cds_epi,
  min_expr = 0.1
)

############################################################
## 22. Preprocess: PCA 30
############################################################

set.seed(12345)

cds_epi <- monocle3::preprocess_cds(
  cds_epi,
  method = "PCA",
  num_dim = 30,
  norm_method = "log",
  scaling = TRUE,
  verbose = FALSE
)

############################################################
## 23. Align by patient
##
## IMPORTANT:
## patient is donor/batch variable
## group is NOT regressed out
############################################################

set.seed(12345)

cds_epi <- monocle3::align_cds(
  cds_epi,
  alignment_group = "patient"
)

############################################################
## 24. NEW UMAP
############################################################

set.seed(12345)

cds_epi <- monocle3::reduce_dimension(
  cds_epi,
  reduction_method = "UMAP",
  preprocess_method = "Aligned",
  umap.metric = "cosine",
  umap.min_dist = 0.3,
  umap.n_neighbors = 40,
  cores = 1,
  verbose = FALSE
)

############################################################
## 25. NEW Monocle clustering
############################################################

set.seed(12345)

cds_epi <- monocle3::cluster_cells(
  cds_epi,
  reduction_method = "UMAP",
  cluster_method = "leiden",
  k = 20,
  num_iter = 2,
  partition_qval = 0.05,
  random_seed = 12345,
  verbose = FALSE
)

############################################################
## 26. NEW principal trajectory graph
############################################################

set.seed(12345)

cds_epi <- monocle3::learn_graph(
  cds_epi,
  use_partition = FALSE,
  close_loop = TRUE,
  learn_graph_control = list(
    minimal_branch_len = 10
  ),
  verbose = FALSE
)

pg_list <- monocle3::principal_graph(
  cds_epi
)

if (is.null(
  pg_list[["UMAP"]]
)) {
  stop(
    "Monocle3未成功建立UMAP principal graph。"
  )
}

############################################################
## 27. Obtain NEW Monocle clusters
############################################################

mono_clusters <- monocle3::clusters(
  cds_epi,
  reduction_method = "UMAP"
)

if (is.null(names(mono_clusters))) {
  stop(
    "新的Monocle cluster vector没有cell barcode names。"
  )
}

if (!setequal(
  names(mono_clusters),
  colnames(cds_epi)
)) {
  stop(
    "新的Monocle cluster与cell barcodes不完全对应。"
  )
}

mono_clusters <- mono_clusters[
  colnames(cds_epi)
]

if (!identical(
  names(mono_clusters),
  colnames(cds_epi)
)) {
  stop(
    "新的Monocle cluster排序失败。"
  )
}

SummarizedExperiment::colData(
  cds_epi
)$mono_cluster <- as.character(
  mono_clusters
)

cat("\n====================================================\n")
cat("NEW CLEAN EPITHELIAL MONOCLE STATES\n")
cat("====================================================\n")

print(
  table(
    SummarizedExperiment::colData(
      cds_epi
    )$mono_cluster
  )
)

############################################################
## 28. Root-state metadata
############################################################

root_meta <- as.data.frame(
  SummarizedExperiment::colData(
    cds_epi
  )
)

root_meta$cell <- rownames(
  root_meta
)

root_meta$patient <- as.character(
  root_meta$patient
)

root_meta$group <- as.character(
  root_meta$group
)

root_meta$mono_cluster <- as.character(
  root_meta$mono_cluster
)

############################################################
## 29. Root candidate summary
############################################################

root_stats <- root_meta |>
  dplyr::group_by(
    mono_cluster
  ) |>
  dplyr::summarise(
    n_cells = dplyr::n(),
    n_normal = sum(
      group == "Normal",
      na.rm = TRUE
    ),
    n_tumor = sum(
      group == "Tumor",
      na.rm = TRUE
    ),
    normal_fraction =
      n_normal / n_cells,
    n_donors =
      dplyr::n_distinct(patient),
    n_normal_donors =
      dplyr::n_distinct(
        patient[group == "Normal"]
      ),
    n_tumor_donors =
      dplyr::n_distinct(
        patient[group == "Tumor"]
      ),
    .groups = "drop"
  ) |>
  dplyr::arrange(
    dplyr::desc(normal_fraction),
    dplyr::desc(n_normal_donors),
    dplyr::desc(n_cells)
  )

cat("\n====================================================\n")
cat("ROOT-STATE CANDIDATES\n")
cat("====================================================\n")

print(
  root_stats,
  n = Inf
)

write.csv(
  root_stats,
  file.path(
    output_dir,
    "CLEAN_Epithelial_root_state_candidate_summary.csv"
  ),
  row.names = FALSE
)

############################################################
## 30. Root selection
##
## Prefer:
## >= 5 Normal donors
## >= 100 cells
############################################################

eligible_root <- root_stats |>
  dplyr::filter(
    n_normal_donors >= 5,
    n_cells >= 100
  )

if (nrow(eligible_root) == 0) {
  
  warning(
    paste0(
      "没有state同时满足>=5 Normal donors和>=100 cells；",
      "改用>=3 Normal donors。"
    )
  )
  
  eligible_root <- root_stats |>
    dplyr::filter(
      n_normal_donors >= 3
    )
}

if (nrow(eligible_root) == 0) {
  stop(
    "没有足够Normal donors支持的root state。"
  )
}

eligible_root <- eligible_root |>
  dplyr::arrange(
    dplyr::desc(normal_fraction),
    dplyr::desc(n_normal_donors),
    dplyr::desc(n_cells)
  )

root_cluster <- as.character(
  eligible_root$mono_cluster[1]
)

cat("\n====================================================\n")
cat("SELECTED ROOT STATE\n")
cat("====================================================\n")

print(
  eligible_root[1, ],
  row.names = FALSE
)

############################################################
## 31. Root cells
############################################################

root_cells <- root_meta$cell[
  root_meta$mono_cluster ==
    root_cluster
]

if (length(root_cells) == 0) {
  stop(
    "选择的root state中没有细胞。"
  )
}

############################################################
## 32. Closest principal graph vertex
############################################################

pg_aux <- cds_epi@principal_graph_aux

if (
  is.null(pg_aux) ||
  !("UMAP" %in% names(pg_aux))
) {
  stop(
    "principal_graph_aux中不存在UMAP。"
  )
}

umap_aux <- pg_aux[["UMAP"]]

closest_vertex <-
  umap_aux$pr_graph_cell_proj_closest_vertex

if (is.null(
  closest_vertex
)) {
  stop(
    "无法取得closest principal graph vertex。"
  )
}

closest_vertex <- as.matrix(
  closest_vertex
)

if (is.null(
  rownames(closest_vertex)
)) {
  stop(
    "closest_vertex没有cell names。"
  )
}

############################################################
## 33. Principal graph node names
############################################################

pg <- monocle3::principal_graph(
  cds_epi
)[["UMAP"]]

graph_vertex_names <- igraph::V(
  pg
)$name

if (length(
  graph_vertex_names
) == 0) {
  stop(
    "principal graph不存在vertex names。"
  )
}

############################################################
## 34. Match root cells to graph
############################################################

root_cells_graph <- intersect(
  root_cells,
  rownames(closest_vertex)
)

if (length(
  root_cells_graph
) == 0) {
  stop(
    "Root cells无法匹配principal graph。"
  )
}

raw_vertex <- as.character(
  closest_vertex[
    root_cells_graph,
    1
  ]
)

if (
  all(
    raw_vertex %in%
    graph_vertex_names
  )
) {
  
  mapped_vertex <- raw_vertex
  
} else {
  
  vertex_index <- suppressWarnings(
    as.integer(
      raw_vertex
    )
  )
  
  valid_index <- (
    !is.na(vertex_index) &
      vertex_index >= 1 &
      vertex_index <=
      length(graph_vertex_names)
  )
  
  if (!all(
    valid_index
  )) {
    stop(
      "无法将closest vertex映射到principal graph node。"
    )
  }
  
  mapped_vertex <-
    graph_vertex_names[
      vertex_index
    ]
}

############################################################
## 35. Select representative root graph node
############################################################

root_node_table <- sort(
  table(mapped_vertex),
  decreasing = TRUE
)

if (length(
  root_node_table
) == 0) {
  stop(
    "无法确定root principal node。"
  )
}

root_pr_node <- names(
  root_node_table
)[1]

cat("\n====================================================\n")
cat("SELECTED ROOT PRINCIPAL NODE\n")
cat("====================================================\n")

cat(
  "Root state:",
  root_cluster,
  "\n"
)

cat(
  "Root principal node:",
  root_pr_node,
  "\n"
)

############################################################
## 36. Order cells
############################################################

cds_epi <- monocle3::order_cells(
  cds_epi,
  reduction_method = "UMAP",
  root_pr_nodes = root_pr_node
)

############################################################
## 37. Pseudotime
############################################################

pt <- monocle3::pseudotime(
  cds_epi
)

SummarizedExperiment::colData(
  cds_epi
)$pseudotime <- as.numeric(
  pt
)

finite_pt <- is.finite(
  pt
)

finite_n <- sum(
  finite_pt
)

finite_pct <- mean(
  finite_pt
) * 100

cat("\n====================================================\n")
cat("PSEUDOTIME QC\n")
cat("====================================================\n")

cat(
  "Finite pseudotime:",
  finite_n,
  "/",
  length(pt),
  " (",
  round(
    finite_pct,
    2
  ),
  "%)\n",
  sep = ""
)

############################################################
## 38. Save cell-level pseudotime metadata
############################################################

pt_df <- data.frame(
  cell = colnames(cds_epi),
  patient = as.character(
    SummarizedExperiment::colData(
      cds_epi
    )$patient
  ),
  group = as.character(
    SummarizedExperiment::colData(
      cds_epi
    )$group
  ),
  mono_cluster = as.character(
    SummarizedExperiment::colData(
      cds_epi
    )$mono_cluster
  ),
  pseudotime = as.numeric(
    SummarizedExperiment::colData(
      cds_epi
    )$pseudotime
  ),
  stringsAsFactors = FALSE
)

write.csv(
  pt_df,
  file.path(
    output_dir,
    "CLEAN_Epithelial_cell_pseudotime_metadata.csv"
  ),
  row.names = FALSE
)

############################################################
## 39. Save root information
############################################################

root_info <- data.frame(
  root_cluster =
    root_cluster,
  root_pr_node =
    root_pr_node,
  n_cells =
    eligible_root$n_cells[1],
  normal_fraction =
    eligible_root$normal_fraction[1],
  n_normal_donors =
    eligible_root$n_normal_donors[1],
  n_tumor_donors =
    eligible_root$n_tumor_donors[1],
  stringsAsFactors = FALSE
)

write.csv(
  root_info,
  file.path(
    output_dir,
    "CLEAN_Epithelial_pseudotime_root_information.csv"
  ),
  row.names = FALSE
)

############################################################
## 40. Figure 8A
############################################################

p8A <- monocle3::plot_cells(
  cds_epi,
  reduction_method = "UMAP",
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE,
  trajectory_graph_color = "black",
  trajectory_graph_segment_size = 0.8,
  cell_size = 0.50,
  alpha = 0.88,
  label_cell_groups = FALSE,
  label_groups_by_cluster = FALSE,
  label_roots = TRUE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  graph_label_size = 4
) +
  ggplot2::labs(
    x = "UMAP 1",
    y = "UMAP 2",
    color = "Pseudotime"
  ) +
  ggplot2::ggtitle(NULL) +
  ggplot2::theme_classic(
    base_size = 20
  ) +
  ggplot2::theme(
    axis.title = ggplot2::element_text(
      size = 25,
      face = "bold",
      color = "black"
    ),
    axis.text = ggplot2::element_text(
      size = 18,
      color = "black"
    ),
    axis.line = ggplot2::element_line(
      linewidth = 1.1,
      color = "black"
    ),
    axis.ticks = ggplot2::element_line(
      linewidth = 1,
      color = "black"
    ),
    legend.title = ggplot2::element_text(
      size = 20,
      face = "bold",
      color = "black"
    ),
    legend.text = ggplot2::element_text(
      size = 17,
      color = "black"
    ),
    panel.grid =
      ggplot2::element_blank(),
    plot.margin =
      ggplot2::margin(
        15,
        15,
        15,
        15
      )
  )

############################################################
## 41. Figure 8B
############################################################

p8B <- monocle3::plot_cells(
  cds_epi,
  reduction_method = "UMAP",
  color_cells_by = "group",
  show_trajectory_graph = TRUE,
  trajectory_graph_color = "black",
  trajectory_graph_segment_size = 0.8,
  cell_size = 0.50,
  alpha = 0.88,
  label_cell_groups = FALSE,
  label_groups_by_cluster = FALSE,
  label_roots = TRUE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  graph_label_size = 4
) +
  ggplot2::scale_color_manual(
    values = c(
      "Normal" = "#3C8DBC",
      "Tumor" = "#D95F59"
    ),
    drop = FALSE
  ) +
  ggplot2::labs(
    x = "UMAP 1",
    y = "UMAP 2",
    color = "Group"
  ) +
  ggplot2::ggtitle(NULL) +
  ggplot2::theme_classic(
    base_size = 20
  ) +
  ggplot2::theme(
    axis.title = ggplot2::element_text(
      size = 25,
      face = "bold",
      color = "black"
    ),
    axis.text = ggplot2::element_text(
      size = 18,
      color = "black"
    ),
    axis.line = ggplot2::element_line(
      linewidth = 1.1,
      color = "black"
    ),
    axis.ticks = ggplot2::element_line(
      linewidth = 1,
      color = "black"
    ),
    legend.title = ggplot2::element_text(
      size = 20,
      face = "bold",
      color = "black"
    ),
    legend.text = ggplot2::element_text(
      size = 17,
      color = "black"
    ),
    panel.grid =
      ggplot2::element_blank(),
    plot.margin =
      ggplot2::margin(
        15,
        15,
        15,
        15
      )
  )

############################################################
## 42. Save A
############################################################

ggplot2::ggsave(
  file.path(
    output_dir,
    "Figure8A_CLEAN_Epithelial_pseudotime_trajectory.pdf"
  ),
  p8A,
  width = 7.5,
  height = 6.2,
  device = cairo_pdf,
  bg = "white"
)

ggplot2::ggsave(
  file.path(
    output_dir,
    "Figure8A_CLEAN_Epithelial_pseudotime_trajectory.png"
  ),
  p8A,
  width = 7.5,
  height = 6.2,
  dpi = 600,
  bg = "white"
)

############################################################
## 43. Save B
############################################################

ggplot2::ggsave(
  file.path(
    output_dir,
    "Figure8B_CLEAN_Epithelial_Normal_Tumor_trajectory.pdf"
  ),
  p8B,
  width = 7.5,
  height = 6.2,
  device = cairo_pdf,
  bg = "white"
)

ggplot2::ggsave(
  file.path(
    output_dir,
    "Figure8B_CLEAN_Epithelial_Normal_Tumor_trajectory.png"
  ),
  p8B,
  width = 7.5,
  height = 6.2,
  dpi = 600,
  bg = "white"
)

############################################################
## 44. Combined A-B
############################################################

p_AB <- patchwork::wrap_plots(
  p8A,
  p8B,
  ncol = 2
) +
  patchwork::plot_annotation(
    tag_levels = "A"
  ) &
  ggplot2::theme(
    plot.tag =
      ggplot2::element_text(
        size = 28,
        face = "bold",
        color = "black"
      )
  )

ggplot2::ggsave(
  file.path(
    output_dir,
    "Figure8AB_CLEAN_Epithelial_pseudotime_combined.pdf"
  ),
  p_AB,
  width = 15,
  height = 6.2,
  device = cairo_pdf,
  bg = "white"
)

ggplot2::ggsave(
  file.path(
    output_dir,
    "Figure8AB_CLEAN_Epithelial_pseudotime_combined.png"
  ),
  p_AB,
  width = 15,
  height = 6.2,
  dpi = 600,
  bg = "white"
)

############################################################
## 45. Save NEW complete Monocle3 object
############################################################

if (dir.exists(
  new_monocle_dir
)) {
  unlink(
    new_monocle_dir,
    recursive = TRUE,
    force = TRUE
  )
}

monocle3::save_monocle_objects(
  cds = cds_epi,
  directory_path = new_monocle_dir,
  comment = paste0(
    "GSE294300 CLEAN epithelial Monocle3 pseudotime; ",
    "old Ep6/Ep7/Ep8 removed before complete reprocessing"
  ),
  archive_control = list(
    archive_type = "none"
  )
)

############################################################
## 46. Session info
############################################################

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    output_dir,
    "sessionInfo.txt"
  )
)

############################################################
## 47. Final report
############################################################

cat("\n############################################################\n")
cat("## CLEAN FIGURE 8A-B COMPLETED\n")
cat("############################################################\n")

cat(
  "\nOriginal epithelial cells:",
  length(original_epi_cells),
  "\n"
)

cat(
  "Removed immune-like Ep6/Ep7/Ep8:",
  length(remove_cells),
  "\n"
)

cat(
  "Final clean epithelial cells:",
  ncol(cds_epi),
  "\n"
)

cat(
  "New Monocle states:",
  length(
    unique(
      SummarizedExperiment::colData(
        cds_epi
      )$mono_cluster
    )
  ),
  "\n"
)

cat(
  "Selected root state:",
  root_cluster,
  "\n"
)

cat(
  "Selected root principal node:",
  root_pr_node,
  "\n"
)

cat(
  "Finite pseudotime:",
  finite_n,
  "/",
  length(pt),
  " (",
  round(
    finite_pct,
    2
  ),
  "%)\n",
  sep = ""
)

cat(
  "\nOutput directory:\n",
  output_dir,
  "\n"
)

cat(
  "\nNEW Monocle directory:\n",
  new_monocle_dir,
  "\n"
)

############################################################
## GSE294300
## CLEAN epithelial Ep1-Ep8
## Map NEW states onto EXACT SAME UMAP as Figure 8A/B
############################################################

rm(list = ls())
options(stringsAsFactors = FALSE)
Sys.setenv(LANGUAGE = "en")
set.seed(12345)

############################################################
## 1. Required packages
############################################################

required_packages <- c(
  "monocle3",
  "SummarizedExperiment",
  "SingleCellExperiment",
  "dplyr",
  "ggplot2"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  stop(
    "缺少以下R包：",
    paste(missing_packages, collapse = ", ")
  )
}

############################################################
## 2. Paths
############################################################

base_dir <- "C:/Users/92066/Desktop/单细胞孟德尔/Figure8/Figure8_AB_CLEAN"

monocle_dir <- file.path(
  base_dir,
  "GSE294300_CLEAN_Epithelial_Monocle3_pseudotime"
)

output_dir <- file.path(
  base_dir,
  "Figure8_NEW_Ep1_Ep8_on_UMAP"
)

if (!dir.exists(monocle_dir)) {
  stop(
    "找不到新的clean Monocle3对象：",
    monocle_dir
  )
}

if (!dir.exists(output_dir)) {
  dir.create(
    output_dir,
    recursive = TRUE
  )
}

############################################################
## 3. Load NEW clean Monocle3 object
############################################################

cds_epi <- monocle3::load_monocle_objects(
  directory_path = monocle_dir
)

cat("\n========================================\n")
cat("CLEAN MONOCLE3 OBJECT\n")
cat("========================================\n")

cat(
  "Cells:",
  ncol(cds_epi),
  "\n"
)

cat(
  "Genes:",
  nrow(cds_epi),
  "\n"
)

############################################################
## 4. Obtain CURRENT Monocle clusters directly
##
## IMPORTANT:
## Do NOT use any old Ep1-Ep12 metadata.
############################################################

mono_clusters <- monocle3::clusters(
  cds_epi,
  reduction_method = "UMAP"
)

if (is.null(names(mono_clusters))) {
  stop(
    "monocle3::clusters()没有返回cell barcode names。"
  )
}

if (!setequal(
  names(mono_clusters),
  colnames(cds_epi)
)) {
  stop(
    "Monocle cluster assignments与cell barcodes不完全一致。"
  )
}

############################################################
## 5. Strictly reorder by cell barcode
############################################################

mono_clusters <- mono_clusters[
  colnames(cds_epi)
]

if (!identical(
  names(mono_clusters),
  colnames(cds_epi)
)) {
  stop(
    "Monocle cluster无法按照cell barcode正确排序。"
  )
}

############################################################
## 6. Define NEW Ep1-Ep8
############################################################

cluster_values <- as.character(
  mono_clusters
)

cluster_numeric <- suppressWarnings(
  as.numeric(
    unique(cluster_values)
  )
)

if (all(!is.na(cluster_numeric))) {
  
  state_levels <- as.character(
    sort(
      unique(
        as.numeric(
          cluster_values
        )
      )
    )
  )
  
} else {
  
  state_levels <- sort(
    unique(
      cluster_values
    )
  )
}

ep_levels <- paste0(
  "Ep",
  state_levels
)

if (length(ep_levels) != 8) {
  stop(
    "当前clean Monocle对象不是8个states，实际为：",
    length(ep_levels)
  )
}

Ep_state <- factor(
  paste0(
    "Ep",
    cluster_values
  ),
  levels = ep_levels
)

names(Ep_state) <- colnames(cds_epi)

SummarizedExperiment::colData(
  cds_epi
)$Ep_state <- Ep_state

############################################################
## 7. Check state counts
############################################################

state_counts <- table(
  SummarizedExperiment::colData(
    cds_epi
  )$Ep_state
)

cat("\n========================================\n")
cat("NEW Ep1-Ep8 COUNTS\n")
cat("========================================\n")

print(
  state_counts
)

cat(
  "\nTotal cells:",
  sum(state_counts),
  "\n"
)

if (sum(state_counts) != 12559) {
  stop(
    "clean epithelial cell总数不是12559，请检查对象。"
  )
}

############################################################
## 8. Extract EXACT SAME UMAP coordinates as A/B
##
## NO new UMAP is calculated here.
############################################################

umap_mat <- SingleCellExperiment::reducedDim(
  cds_epi,
  "UMAP"
)

if (is.null(umap_mat)) {
  stop(
    "clean Monocle对象中不存在UMAP。"
  )
}

if (ncol(umap_mat) != 2) {
  stop(
    "UMAP不是二维坐标。"
  )
}

if (is.null(rownames(umap_mat))) {
  stop(
    "UMAP坐标没有cell barcode rownames。"
  )
}

if (!setequal(
  rownames(umap_mat),
  colnames(cds_epi)
)) {
  stop(
    "UMAP cell barcodes与Monocle对象不一致。"
  )
}

umap_mat <- umap_mat[
  colnames(cds_epi),
  ,
  drop = FALSE
]

############################################################
## 9. Build cell-level UMAP table
############################################################

umap_df <- data.frame(
  cell = colnames(cds_epi),
  UMAP_1 = umap_mat[, 1],
  UMAP_2 = umap_mat[, 2],
  Ep_state = as.character(
    SummarizedExperiment::colData(
      cds_epi
    )$Ep_state
  ),
  stringsAsFactors = FALSE
)

############################################################
## 10. Calculate robust label positions
##
## Step 1:
## calculate median center of each state
##
## Step 2:
## find the REAL cell nearest that center
##
## Therefore each Ep label sits on a real cell region.
############################################################

state_center <- umap_df |>
  dplyr::group_by(
    Ep_state
  ) |>
  dplyr::summarise(
    center_x = stats::median(
      UMAP_1
    ),
    center_y = stats::median(
      UMAP_2
    ),
    n_cells = dplyr::n(),
    .groups = "drop"
  )

label_df <- umap_df |>
  dplyr::left_join(
    state_center,
    by = "Ep_state"
  ) |>
  dplyr::mutate(
    distance_to_center =
      (UMAP_1 - center_x)^2 +
      (UMAP_2 - center_y)^2
  ) |>
  dplyr::group_by(
    Ep_state
  ) |>
  dplyr::slice_min(
    order_by = distance_to_center,
    n = 1,
    with_ties = FALSE
  ) |>
  dplyr::ungroup() |>
  dplyr::select(
    Ep_state,
    UMAP_1,
    UMAP_2,
    n_cells
  )

label_df$Ep_state <- factor(
  label_df$Ep_state,
  levels = ep_levels
)

label_df <- label_df |>
  dplyr::arrange(
    Ep_state
  )

cat("\n========================================\n")
cat("LABEL POSITIONS\n")
cat("========================================\n")

print(
  label_df,
  n = Inf
)

############################################################
## 11. Save NEW Ep-state cell mapping
##
## Important for later donor-level composition analysis.
############################################################

state_metadata <- data.frame(
  cell = colnames(cds_epi),
  Ep_state = as.character(
    SummarizedExperiment::colData(
      cds_epi
    )$Ep_state
  ),
  mono_cluster = as.character(
    mono_clusters
  ),
  patient = as.character(
    SummarizedExperiment::colData(
      cds_epi
    )$patient
  ),
  group = as.character(
    SummarizedExperiment::colData(
      cds_epi
    )$group
  ),
  pseudotime = as.numeric(
    SummarizedExperiment::colData(
      cds_epi
    )$pseudotime
  ),
  stringsAsFactors = FALSE
)

write.csv(
  state_metadata,
  file.path(
    output_dir,
    "NEW_Ep1_Ep8_cell_metadata.csv"
  ),
  row.names = FALSE
)

############################################################
## 12. Plot NEW Ep1-Ep8 on EXACT SAME UMAP
##
## label_cell_groups = FALSE:
## avoid Monocle automatic label / viewport problem.
##
## label_roots = TRUE:
## retain root label "1", same as A/B.
############################################################

p_ep <- monocle3::plot_cells(
  cds_epi,
  reduction_method = "UMAP",
  color_cells_by = "Ep_state",
  show_trajectory_graph = TRUE,
  trajectory_graph_color = "black",
  trajectory_graph_segment_size = 0.8,
  cell_size = 0.50,
  alpha = 0.88,
  label_cell_groups = FALSE,
  label_groups_by_cluster = FALSE,
  label_roots = TRUE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  graph_label_size = 4
) +
  
  ggplot2::geom_text(
    data = label_df,
    ggplot2::aes(
      x = UMAP_1,
      y = UMAP_2,
      label = Ep_state
    ),
    inherit.aes = FALSE,
    size = 8.5,
    fontface = "bold",
    color = "black"
  ) +
  
  ggplot2::labs(
    x = "UMAP 1",
    y = "UMAP 2",
    color = "Epithelial state"
  ) +
  
  ggplot2::ggtitle(NULL) +
  
  ggplot2::theme_classic(
    base_size = 20
  ) +
  
  ggplot2::theme(
    axis.title = ggplot2::element_text(
      size = 25,
      face = "bold",
      color = "black"
    ),
    
    axis.text = ggplot2::element_text(
      size = 18,
      color = "black"
    ),
    
    axis.line = ggplot2::element_line(
      linewidth = 1.1,
      color = "black"
    ),
    
    axis.ticks = ggplot2::element_line(
      linewidth = 1,
      color = "black"
    ),
    
    legend.title = ggplot2::element_text(
      size = 20,
      face = "bold",
      color = "black"
    ),
    
    legend.text = ggplot2::element_text(
      size = 17,
      color = "black"
    ),
    
    panel.grid = ggplot2::element_blank(),
    
    plot.margin = ggplot2::margin(
      15,
      20,
      15,
      15
    )
  )

############################################################
## 13. Save plot
############################################################

ggplot2::ggsave(
  filename = file.path(
    output_dir,
    "Figure8_NEW_Ep1_Ep8_same_UMAP.pdf"
  ),
  plot = p_ep,
  width = 8.5,
  height = 7,
  device = cairo_pdf,
  bg = "white"
)

ggplot2::ggsave(
  filename = file.path(
    output_dir,
    "Figure8_NEW_Ep1_Ep8_same_UMAP.png"
  ),
  plot = p_ep,
  width = 8.5,
  height = 7,
  dpi = 600,
  bg = "white"
)

############################################################
## 14. Save updated complete Monocle3 object
##
## Ep_state is now explicitly stored in colData.
############################################################

updated_monocle_dir <- file.path(
  output_dir,
  "GSE294300_CLEAN_Epithelial_Monocle3_with_Ep_states"
)

if (dir.exists(
  updated_monocle_dir
)) {
  unlink(
    updated_monocle_dir,
    recursive = TRUE,
    force = TRUE
  )
}

monocle3::save_monocle_objects(
  cds = cds_epi,
  directory_path = updated_monocle_dir,
  comment = paste0(
    "GSE294300 clean epithelial Monocle3 object; ",
    "new Ep1-Ep8 states added to colData"
  ),
  archive_control = list(
    archive_type = "none"
  )
)

############################################################
## 15. Finished
############################################################

cat("\n############################################################\n")
cat("## NEW Ep1-Ep8 UMAP COMPLETED\n")
cat("############################################################\n")

cat(
  "\nStates:",
  paste(
    ep_levels,
    collapse = ", "
  ),
  "\n"
)

cat(
  "Total cells:",
  ncol(cds_epi),
  "\n"
)

cat(
  "\nOutput directory:\n",
  output_dir,
  "\n"
)

cat(
  "\nUpdated Monocle directory:\n",
  updated_monocle_dir,
  "\n"
)
############################################################
## GSE294300
## Figure 8C REDO
##
## CLEAN epithelial Ep1-Ep8
## Donor-level state composition analysis
##
## Unit of analysis: patient
## Comparison: paired Tumor vs Normal
## Proportion:
## cells in one Ep state / all clean epithelial cells
##
## Statistics:
## paired Wilcoxon signed-rank test
## BH-FDR across ALL 8 states
############################################################

rm(list = ls())
options(stringsAsFactors = FALSE)
Sys.setenv(LANGUAGE = "en")
set.seed(12345)

############################################################
## 1. Required packages
############################################################

required_packages <- c(
  "monocle3",
  "SummarizedExperiment",
  "dplyr",
  "tidyr",
  "ggplot2"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  stop(
    "缺少以下R包：",
    paste(missing_packages, collapse = ", ")
  )
}

############################################################
## 2. Paths
############################################################

base_dir <- "C:/Users/92066/Desktop/单细胞孟德尔/Figure8/Figure8_AB_CLEAN"

monocle_dir <- file.path(
  base_dir,
  "GSE294300_CLEAN_Epithelial_Monocle3_pseudotime"
)

output_dir <- file.path(
  base_dir,
  "Figure8C_NEW_Ep1_Ep8_composition"
)

if (!dir.exists(monocle_dir)) {
  stop(
    "找不到clean Monocle3对象：",
    monocle_dir
  )
}

if (!dir.exists(output_dir)) {
  dir.create(
    output_dir,
    recursive = TRUE
  )
}

############################################################
## 3. Load clean Monocle3 object
############################################################

cds_epi <- monocle3::load_monocle_objects(
  directory_path = monocle_dir
)

cat("\n====================================================\n")
cat("CLEAN EPITHELIAL MONOCLE3 OBJECT\n")
cat("====================================================\n")

cat(
  "Cells:",
  ncol(cds_epi),
  "\n"
)

if (ncol(cds_epi) != 12559) {
  stop(
    "clean epithelial cell数不是预期的12559。实际：",
    ncol(cds_epi)
  )
}

############################################################
## 4. Obtain CURRENT Monocle clusters directly
############################################################

mono_cluster <- monocle3::clusters(
  cds_epi,
  reduction_method = "UMAP"
)

if (is.null(names(mono_cluster))) {
  stop(
    "monocle3::clusters()没有cell barcode names。"
  )
}

if (!setequal(
  names(mono_cluster),
  colnames(cds_epi)
)) {
  stop(
    "Monocle cluster与cell barcodes不完全一致。"
  )
}

mono_cluster <- mono_cluster[
  colnames(cds_epi)
]

if (!identical(
  names(mono_cluster),
  colnames(cds_epi)
)) {
  stop(
    "Monocle cluster排序失败。"
  )
}

############################################################
## 5. Define NEW Ep1-Ep8
############################################################

cluster_values <- as.character(
  mono_cluster
)

cluster_numeric <- suppressWarnings(
  as.numeric(
    unique(cluster_values)
  )
)

if (all(!is.na(cluster_numeric))) {
  
  state_levels <- as.character(
    sort(
      unique(
        as.numeric(cluster_values)
      )
    )
  )
  
} else {
  
  state_levels <- sort(
    unique(cluster_values)
  )
}

if (length(state_levels) != 8) {
  stop(
    "当前Monocle对象不是8个states，实际：",
    length(state_levels)
  )
}

ep_levels <- paste0(
  "Ep",
  state_levels
)

Ep_state <- factor(
  paste0(
    "Ep",
    cluster_values
  ),
  levels = ep_levels
)

names(Ep_state) <- colnames(cds_epi)

cat("\n====================================================\n")
cat("NEW Ep1-Ep8 COUNTS\n")
cat("====================================================\n")

print(
  table(Ep_state)
)

############################################################
## 6. Extract metadata
############################################################

cell_meta <- as.data.frame(
  SummarizedExperiment::colData(
    cds_epi
  )
)

required_meta <- c(
  "patient",
  "group"
)

missing_meta <- setdiff(
  required_meta,
  colnames(cell_meta)
)

if (length(missing_meta) > 0) {
  stop(
    "Monocle metadata缺少：",
    paste(missing_meta, collapse = ", ")
  )
}

############################################################
## 7. Build cell-level table
############################################################

cell_df <- data.frame(
  cell = colnames(cds_epi),
  patient = as.character(
    cell_meta[
      colnames(cds_epi),
      "patient"
    ]
  ),
  group = as.character(
    cell_meta[
      colnames(cds_epi),
      "group"
    ]
  ),
  mono_cluster = cluster_values,
  Ep_state = as.character(Ep_state),
  stringsAsFactors = FALSE
)

if (anyNA(cell_df$patient)) {
  stop("patient中存在NA。")
}

if (anyNA(cell_df$Ep_state)) {
  stop("Ep_state中存在NA。")
}

if (any(
  !cell_df$group %in%
  c("Normal", "Tumor")
)) {
  stop(
    "group中出现Normal/Tumor以外的值。"
  )
}

############################################################
## 8. Total clean epithelial cells per patient × group
############################################################

sample_total <- cell_df |>
  dplyr::group_by(
    patient,
    group
  ) |>
  dplyr::summarise(
    total_epi_cells = dplyr::n(),
    .groups = "drop"
  ) |>
  dplyr::arrange(
    patient,
    group
  )

cat("\n====================================================\n")
cat("CLEAN EPITHELIAL CELLS PER PATIENT × GROUP\n")
cat("====================================================\n")

print(
  sample_total,
  n = Inf
)

write.csv(
  sample_total,
  file.path(
    output_dir,
    "CLEAN_Epithelial_total_cells_per_patient_group.csv"
  ),
  row.names = FALSE
)

############################################################
## 9. Minimum epithelial-cell requirement
############################################################

min_cells <- 20L

valid_samples <- sample_total |>
  dplyr::filter(
    total_epi_cells >= min_cells
  )

############################################################
## 10. Identify strictly paired patients
############################################################

pair_check <- valid_samples |>
  dplyr::group_by(
    patient
  ) |>
  dplyr::summarise(
    has_normal = any(
      group == "Normal"
    ),
    has_tumor = any(
      group == "Tumor"
    ),
    n_groups = dplyr::n_distinct(
      group
    ),
    .groups = "drop"
  ) |>
  dplyr::filter(
    has_normal,
    has_tumor,
    n_groups == 2
  )

paired_patients <- pair_check$patient

cat("\n====================================================\n")
cat("PAIRED PATIENTS\n")
cat("====================================================\n")

cat(
  "Number of paired patients:",
  length(paired_patients),
  "\n"
)

print(
  paired_patients
)

if (length(paired_patients) < 5) {
  stop(
    "满足>=20个clean epithelial cells且具有Normal/Tumor配对的患者少于5例。"
  )
}

############################################################
## 11. Keep valid paired samples only
############################################################

cell_paired <- cell_df |>
  dplyr::filter(
    patient %in% paired_patients
  ) |>
  dplyr::inner_join(
    valid_samples,
    by = c(
      "patient",
      "group"
    )
  )

############################################################
## 12. Count each Ep state per patient × group
############################################################

state_counts <- cell_paired |>
  dplyr::group_by(
    patient,
    group,
    Ep_state
  ) |>
  dplyr::summarise(
    n_state_cells = dplyr::n(),
    .groups = "drop"
  )

############################################################
## 13. Complete missing states with ZERO
##
## Essential:
## if one state is absent in one sample,
## its count must be explicitly 0.
############################################################

state_grid <- tidyr::expand_grid(
  patient = paired_patients,
  group = c(
    "Normal",
    "Tumor"
  ),
  Ep_state = ep_levels
)

state_counts_complete <- state_grid |>
  dplyr::left_join(
    state_counts,
    by = c(
      "patient",
      "group",
      "Ep_state"
    )
  ) |>
  dplyr::mutate(
    n_state_cells = tidyr::replace_na(
      n_state_cells,
      0L
    )
  )

############################################################
## 14. Add total epithelial cell numbers
############################################################

composition_df <- state_counts_complete |>
  dplyr::left_join(
    sample_total,
    by = c(
      "patient",
      "group"
    )
  )

if (anyNA(
  composition_df$total_epi_cells
)) {
  stop(
    "部分patient×group无法匹配total epithelial cell数。"
  )
}

############################################################
## 15. Calculate proportions
############################################################

composition_df <- composition_df |>
  dplyr::mutate(
    proportion =
      n_state_cells /
      total_epi_cells,
    
    percentage =
      proportion * 100
  )

composition_df$group <- factor(
  composition_df$group,
  levels = c(
    "Normal",
    "Tumor"
  )
)

composition_df$Ep_state <- factor(
  composition_df$Ep_state,
  levels = ep_levels
)

write.csv(
  composition_df,
  file.path(
    output_dir,
    "NEW_Ep1_Ep8_donor_level_composition.csv"
  ),
  row.names = FALSE
)

############################################################
## 16. Verify every patient has exactly 16 rows
## 8 states × 2 groups
############################################################

row_check <- composition_df |>
  dplyr::group_by(
    patient
  ) |>
  dplyr::summarise(
    n_rows = dplyr::n(),
    .groups = "drop"
  )

if (any(
  row_check$n_rows != 16
)) {
  stop(
    "部分患者并非8 states × 2 groups = 16 rows。"
  )
}

############################################################
## 17. Paired Wilcoxon for each state
############################################################

statistics_list <- lapply(
  ep_levels,
  function(state_now) {
    
    tmp <- composition_df |>
      dplyr::filter(
        Ep_state == state_now
      ) |>
      dplyr::select(
        patient,
        group,
        percentage
      ) |>
      tidyr::pivot_wider(
        names_from = group,
        values_from = percentage
      ) |>
      dplyr::filter(
        !is.na(Normal),
        !is.na(Tumor)
      ) |>
      dplyr::arrange(
        patient
      )
    
    if (nrow(tmp) < 5) {
      
      return(
        data.frame(
          Ep_state = state_now,
          n_pairs = nrow(tmp),
          median_Normal = NA_real_,
          median_Tumor = NA_real_,
          median_delta = NA_real_,
          p_value = NA_real_,
          stringsAsFactors = FALSE
        )
      )
    }
    
    delta <- tmp$Tumor -
      tmp$Normal
    
    wt <- stats::wilcox.test(
      x = tmp$Tumor,
      y = tmp$Normal,
      paired = TRUE,
      exact = FALSE
    )
    
    data.frame(
      Ep_state = state_now,
      n_pairs = nrow(tmp),
      
      median_Normal =
        stats::median(
          tmp$Normal
        ),
      
      median_Tumor =
        stats::median(
          tmp$Tumor
        ),
      
      median_delta =
        stats::median(
          delta
        ),
      
      p_value =
        wt$p.value,
      
      stringsAsFactors = FALSE
    )
  }
)

state_statistics <- do.call(
  rbind,
  statistics_list
)

############################################################
## 18. BH correction across ALL 8 states
############################################################

state_statistics$FDR <- stats::p.adjust(
  state_statistics$p_value,
  method = "BH"
)

state_statistics$significance <- dplyr::case_when(
  is.na(state_statistics$FDR) ~ "NA",
  state_statistics$FDR < 0.001 ~ "***",
  state_statistics$FDR < 0.01 ~ "**",
  state_statistics$FDR < 0.05 ~ "*",
  TRUE ~ "ns"
)

state_statistics$direction <- dplyr::case_when(
  is.na(state_statistics$median_delta) ~ NA_character_,
  state_statistics$median_delta > 0 ~ "Tumor-enriched",
  state_statistics$median_delta < 0 ~ "Tumor-depleted",
  TRUE ~ "No median change"
)

state_statistics$state_number <- suppressWarnings(
  as.numeric(
    sub(
      "^Ep",
      "",
      state_statistics$Ep_state
    )
  )
)

state_statistics <- state_statistics |>
  dplyr::arrange(
    state_number
  )

cat("\n====================================================\n")
cat("DONOR-LEVEL Ep1-Ep8 STATISTICS\n")
cat("====================================================\n")

print(
  state_statistics,
  row.names = FALSE
)

write.csv(
  state_statistics,
  file.path(
    output_dir,
    "NEW_Ep1_Ep8_donor_level_statistics.csv"
  ),
  row.names = FALSE
)

############################################################
## 19. Significant states
############################################################

significant_states <- state_statistics |>
  dplyr::filter(
    !is.na(FDR),
    FDR < 0.05
  )

cat("\n====================================================\n")
cat("SIGNIFICANT STATES: FDR < 0.05\n")
cat("====================================================\n")

print(
  significant_states,
  row.names = FALSE
)

write.csv(
  significant_states,
  file.path(
    output_dir,
    "NEW_Ep1_Ep8_significant_states_FDR005.csv"
  ),
  row.names = FALSE
)

############################################################
## 20. Prepare paired change table
############################################################

paired_change_df <- composition_df |>
  dplyr::select(
    patient,
    group,
    Ep_state,
    percentage
  ) |>
  tidyr::pivot_wider(
    names_from = group,
    values_from = percentage
  ) |>
  dplyr::mutate(
    delta = Tumor - Normal
  )

paired_change_df$Ep_state <- factor(
  paired_change_df$Ep_state,
  levels = ep_levels
)

write.csv(
  paired_change_df,
  file.path(
    output_dir,
    "NEW_Ep1_Ep8_paired_changes.csv"
  ),
  row.names = FALSE
)

############################################################
## 21. Plot annotation data
############################################################

plot_annotation <- state_statistics |>
  dplyr::mutate(
    label = paste0(
      "FDR = ",
      formatC(
        FDR,
        format = "g",
        digits = 2
      )
    )
  )

############################################################
## 22. ALL STATES paired proportion plot
## 4 columns × 2 rows
############################################################

plot_all <- composition_df

plot_all$Ep_state <- factor(
  plot_all$Ep_state,
  levels = ep_levels
)

y_positions <- plot_all |>
  dplyr::group_by(
    Ep_state
  ) |>
  dplyr::summarise(
    ymax = max(
      percentage,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) |>
  dplyr::left_join(
    plot_annotation |>
      dplyr::select(
        Ep_state,
        label
      ),
    by = "Ep_state"
  ) |>
  dplyr::mutate(
    y_position =
      ifelse(
        ymax > 0,
        ymax * 1.12,
        1
      )
  )

p_all <- ggplot2::ggplot(
  plot_all,
  ggplot2::aes(
    x = group,
    y = percentage,
    group = patient
  )
) +
  
  ggplot2::geom_line(
    color = "grey70",
    linewidth = 0.7,
    alpha = 0.8
  ) +
  
  ggplot2::geom_point(
    ggplot2::aes(
      color = group
    ),
    size = 3.3,
    alpha = 0.90
  ) +
  
  ggplot2::stat_summary(
    ggplot2::aes(
      group = group,
      color = group
    ),
    fun = median,
    geom = "crossbar",
    width = 0.45,
    linewidth = 1.1
  ) +
  
  ggplot2::geom_text(
    data = y_positions,
    ggplot2::aes(
      x = 1.5,
      y = y_position,
      label = label
    ),
    inherit.aes = FALSE,
    size = 5.8,
    color = "black"
  ) +
  
  ggplot2::facet_wrap(
    ~ Ep_state,
    ncol = 4,
    scales = "free_y"
  ) +
  
  ggplot2::scale_color_manual(
    values = c(
      "Normal" = "#3C8DBC",
      "Tumor" = "#D95F59"
    ),
    drop = FALSE
  ) +
  
  ggplot2::labs(
    x = NULL,
    y = "Epithelial-state proportion (%)",
    color = NULL
  ) +
  
  ggplot2::theme_classic(
    base_size = 20
  ) +
  
  ggplot2::theme(
    strip.background =
      ggplot2::element_blank(),
    
    strip.text =
      ggplot2::element_text(
        size = 22,
        face = "bold",
        color = "black"
      ),
    
    axis.title.y =
      ggplot2::element_text(
        size = 24,
        face = "bold",
        color = "black"
      ),
    
    axis.text.x =
      ggplot2::element_text(
        size = 18,
        face = "bold",
        color = "black"
      ),
    
    axis.text.y =
      ggplot2::element_text(
        size = 17,
        color = "black"
      ),
    
    legend.position = "none",
    
    axis.line =
      ggplot2::element_line(
        linewidth = 1,
        color = "black"
      ),
    
    panel.grid =
      ggplot2::element_blank(),
    
    plot.margin =
      ggplot2::margin(
        15,
        15,
        15,
        15
      )
  )

ggplot2::ggsave(
  file.path(
    output_dir,
    "Figure8C_NEW_ALL_Ep1_Ep8_paired.pdf"
  ),
  p_all,
  width = 15,
  height = 8.5,
  device = cairo_pdf,
  bg = "white"
)

ggplot2::ggsave(
  file.path(
    output_dir,
    "Figure8C_NEW_ALL_Ep1_Ep8_paired.png"
  ),
  p_all,
  width = 15,
  height = 8.5,
  dpi = 600,
  bg = "white"
)

############################################################
## 23. DELTA plot
## Tumor - Normal percentage points
############################################################

delta_order_df <- paired_change_df |>
  dplyr::group_by(
    Ep_state
  ) |>
  dplyr::summarise(
    median_delta = stats::median(
      delta,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) |>
  dplyr::arrange(
    dplyr::desc(
      median_delta
    )
  )

delta_levels <- as.character(
  delta_order_df$Ep_state
)

paired_change_df$Ep_state_delta <- factor(
  as.character(
    paired_change_df$Ep_state
  ),
  levels = rev(
    delta_levels
  )
)

delta_summary <- state_statistics |>
  dplyr::mutate(
    Ep_state_delta = factor(
      Ep_state,
      levels = rev(
        delta_levels
      )
    )
  )

delta_range <- range(
  paired_change_df$delta,
  na.rm = TRUE
)

delta_span <- diff(
  delta_range
)

if (!is.finite(delta_span) ||
    delta_span == 0) {
  delta_span <- 1
}

delta_summary <- delta_summary |>
  dplyr::mutate(
    label_x =
      median_delta +
      ifelse(
        median_delta >= 0,
        0.025 * delta_span,
        -0.025 * delta_span
      )
  )

p_delta <- ggplot2::ggplot(
  paired_change_df,
  ggplot2::aes(
    x = delta,
    y = Ep_state_delta
  )
) +
  
  ggplot2::geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.9,
    color = "grey45"
  ) +
  
  ggplot2::geom_point(
    color = "grey65",
    size = 3.1,
    alpha = 0.85,
    position = ggplot2::position_jitter(
      height = 0.08,
      width = 0
    )
  ) +
  
  ggplot2::geom_point(
    data = delta_summary,
    ggplot2::aes(
      x = median_delta,
      y = Ep_state_delta,
      fill = direction
    ),
    inherit.aes = FALSE,
    shape = 23,
    size = 5.8,
    stroke = 1,
    color = "black"
  ) +
  
  ggplot2::geom_text(
    data = delta_summary,
    ggplot2::aes(
      x = label_x,
      y = Ep_state_delta,
      label = significance
    ),
    inherit.aes = FALSE,
    size = 6,
    fontface = "bold",
    hjust = ifelse(
      delta_summary$median_delta >= 0,
      0,
      1
    )
  ) +
  
  ggplot2::scale_fill_manual(
    values = c(
      "Tumor-enriched" = "#D95F59",
      "Tumor-depleted" = "#3C8DBC",
      "No median change" = "white"
    ),
    drop = FALSE
  ) +
  
  ggplot2::labs(
    x = paste0(
      "Paired change in epithelial-state proportion\n",
      "(Tumor - Normal, percentage points)"
    ),
    y = NULL,
    fill = NULL
  ) +
  
  ggplot2::theme_classic(
    base_size = 20
  ) +
  
  ggplot2::theme(
    axis.title.x =
      ggplot2::element_text(
        size = 23,
        face = "bold",
        color = "black"
      ),
    
    axis.text.x =
      ggplot2::element_text(
        size = 17,
        color = "black"
      ),
    
    axis.text.y =
      ggplot2::element_text(
        size = 20,
        face = "bold",
        color = "black"
      ),
    
    legend.position = "none",
    
    axis.line =
      ggplot2::element_line(
        linewidth = 1,
        color = "black"
      ),
    
    panel.grid =
      ggplot2::element_blank(),
    
    plot.margin =
      ggplot2::margin(
        15,
        30,
        15,
        15
      )
  )

ggplot2::ggsave(
  file.path(
    output_dir,
    "Figure8C_NEW_DELTA_ALL_Ep1_Ep8.pdf"
  ),
  p_delta,
  width = 9,
  height = 7,
  device = cairo_pdf,
  bg = "white"
)

ggplot2::ggsave(
  file.path(
    output_dir,
    "Figure8C_NEW_DELTA_ALL_Ep1_Ep8.png"
  ),
  p_delta,
  width = 9,
  height = 7,
  dpi = 600,
  bg = "white"
)

############################################################
## 24. Significant-state-only paired plot
############################################################

if (nrow(significant_states) > 0) {
  
  sig_order <- significant_states$Ep_state[
    order(
      significant_states$state_number
    )
  ]
  
  sig_plot_df <- composition_df |>
    dplyr::filter(
      Ep_state %in% sig_order
    )
  
  sig_plot_df$Ep_state <- factor(
    as.character(
      sig_plot_df$Ep_state
    ),
    levels = sig_order
  )
  
  sig_annotation <- state_statistics |>
    dplyr::filter(
      Ep_state %in% sig_order
    ) |>
    dplyr::mutate(
      Ep_state = factor(
        Ep_state,
        levels = sig_order
      ),
      label = paste0(
        "FDR = ",
        formatC(
          FDR,
          format = "g",
          digits = 2
        )
      )
    )
  
  sig_y <- sig_plot_df |>
    dplyr::group_by(
      Ep_state
    ) |>
    dplyr::summarise(
      ymax = max(
        percentage,
        na.rm = TRUE
      ),
      .groups = "drop"
    ) |>
    dplyr::left_join(
      sig_annotation |>
        dplyr::select(
          Ep_state,
          label
        ),
      by = "Ep_state"
    ) |>
    dplyr::mutate(
      y_position =
        ifelse(
          ymax > 0,
          ymax * 1.12,
          1
        )
    )
  
  p_sig <- ggplot2::ggplot(
    sig_plot_df,
    ggplot2::aes(
      x = group,
      y = percentage,
      group = patient
    )
  ) +
    
    ggplot2::geom_line(
      color = "grey70",
      linewidth = 0.8,
      alpha = 0.85
    ) +
    
    ggplot2::geom_point(
      ggplot2::aes(
        color = group
      ),
      size = 4,
      alpha = 0.95
    ) +
    
    ggplot2::stat_summary(
      ggplot2::aes(
        group = group,
        color = group
      ),
      fun = median,
      geom = "crossbar",
      width = 0.45,
      linewidth = 1.3
    ) +
    
    ggplot2::geom_text(
      data = sig_y,
      ggplot2::aes(
        x = 1.5,
        y = y_position,
        label = label
      ),
      inherit.aes = FALSE,
      size = 6.3,
      color = "black"
    ) +
    
    ggplot2::facet_wrap(
      ~ Ep_state,
      scales = "free_y",
      ncol = min(
        4,
        length(sig_order)
      )
    ) +
    
    ggplot2::scale_color_manual(
      values = c(
        "Normal" = "#3C8DBC",
        "Tumor" = "#D95F59"
      ),
      drop = FALSE
    ) +
    
    ggplot2::labs(
      x = NULL,
      y = "Epithelial-state proportion (%)",
      color = NULL
    ) +
    
    ggplot2::theme_classic(
      base_size = 21
    ) +
    
    ggplot2::theme(
      strip.background =
        ggplot2::element_blank(),
      
      strip.text =
        ggplot2::element_text(
          size = 24,
          face = "bold",
          color = "black"
        ),
      
      axis.title.y =
        ggplot2::element_text(
          size = 25,
          face = "bold",
          color = "black"
        ),
      
      axis.text.x =
        ggplot2::element_text(
          size = 19,
          face = "bold",
          color = "black"
        ),
      
      axis.text.y =
        ggplot2::element_text(
          size = 18,
          color = "black"
        ),
      
      legend.position = "none",
      
      axis.line =
        ggplot2::element_line(
          linewidth = 1.1,
          color = "black"
        ),
      
      panel.grid =
        ggplot2::element_blank()
    )
  
  sig_width <- max(
    7,
    min(
      15,
      length(sig_order) * 3.7
    )
  )
  
  ggplot2::ggsave(
    file.path(
      output_dir,
      "Figure8C_NEW_SIGNIFICANT_Ep_states.pdf"
    ),
    p_sig,
    width = sig_width,
    height = 6,
    device = cairo_pdf,
    bg = "white"
  )
  
  ggplot2::ggsave(
    file.path(
      output_dir,
      "Figure8C_NEW_SIGNIFICANT_Ep_states.png"
    ),
    p_sig,
    width = sig_width,
    height = 6,
    dpi = 600,
    bg = "white"
  )
  
} else {
  
  cat(
    "\nNo Ep state reached BH-FDR < 0.05; ",
    "significant-state-only plot was not generated.\n"
  )
}

############################################################
## 25. Session info
############################################################

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    output_dir,
    "sessionInfo.txt"
  )
)

############################################################
## 26. Finished
############################################################

cat("\n############################################################\n")
cat("## NEW Ep1-Ep8 COMPOSITION ANALYSIS COMPLETED\n")
cat("############################################################\n")

cat(
  "\nClean epithelial cells:",
  ncol(cds_epi),
  "\n"
)

cat(
  "Paired patients:",
  length(paired_patients),
  "\n"
)

cat(
  "Minimum epithelial cells per patient×group:",
  min_cells,
  "\n"
)

cat(
  "States tested:",
  length(ep_levels),
  "\n"
)

cat(
  "States with FDR < 0.05:",
  nrow(significant_states),
  "\n"
)

cat(
  "\nOutput directory:\n",
  output_dir,
  "\n"
)
############################################################
## GSE294300
## CLEAN Ep1-Ep8 marker genes + state annotation
##
## Focus:
## Ep2 = Tumor-enriched state
## Ep6 = Normal-enriched / root-associated state
##
## IMPORTANT:
## Ep states are obtained directly from CURRENT
## monocle3::clusters() and mapped by exact cell barcode.
############################################################

rm(list = ls())
options(stringsAsFactors = FALSE)
Sys.setenv(LANGUAGE = "en")
set.seed(12345)

############################################################
## 1. Required packages
############################################################

required_packages <- c(
  "Seurat",
  "SeuratObject",
  "monocle3",
  "SummarizedExperiment",
  "dplyr",
  "ggplot2"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  stop(
    "缺少以下R包：",
    paste(missing_packages, collapse = ", ")
  )
}

############################################################
## 2. Paths
############################################################

seurat_file <- "C:/Users/92066/Desktop/单细胞孟德尔/Figure7/GSE294300/GSE294300_scRNA_final_annotated_with_tSNE.rds"

monocle_dir <- "C:/Users/92066/Desktop/单细胞孟德尔/Figure8/Figure8_AB_CLEAN/GSE294300_CLEAN_Epithelial_Monocle3_pseudotime"

output_dir <- "C:/Users/92066/Desktop/单细胞孟德尔/Figure8/Figure8_AB_CLEAN/Figure8_Ep1_Ep8_marker_annotation"

if (!file.exists(seurat_file)) {
  stop("找不到Seurat对象：", seurat_file)
}

if (!dir.exists(monocle_dir)) {
  stop("找不到clean Monocle3对象：", monocle_dir)
}

if (!dir.exists(output_dir)) {
  dir.create(
    output_dir,
    recursive = TRUE
  )
}

############################################################
## 3. Load CLEAN Monocle3 object
############################################################

cds_epi <- monocle3::load_monocle_objects(
  directory_path = monocle_dir
)

cat("\n========================================\n")
cat("CLEAN MONOCLE3 OBJECT\n")
cat("========================================\n")

cat(
  "Cells:",
  ncol(cds_epi),
  "\n"
)

if (ncol(cds_epi) != 12559) {
  stop(
    "clean epithelial cell数不是12559，实际：",
    ncol(cds_epi)
  )
}

############################################################
## 4. Obtain CURRENT Ep clusters directly
############################################################

mono_cluster <- monocle3::clusters(
  cds_epi,
  reduction_method = "UMAP"
)

if (is.null(names(mono_cluster))) {
  stop(
    "monocle3::clusters()没有cell barcode names。"
  )
}

if (!setequal(
  names(mono_cluster),
  colnames(cds_epi)
)) {
  stop(
    "Monocle clusters与cell barcodes不完全一致。"
  )
}

mono_cluster <- mono_cluster[
  colnames(cds_epi)
]

if (!identical(
  names(mono_cluster),
  colnames(cds_epi)
)) {
  stop(
    "Monocle cluster barcode排序失败。"
  )
}

############################################################
## 5. Define NEW Ep1-Ep8
############################################################

cluster_values <- as.character(
  mono_cluster
)

state_levels <- as.character(
  sort(
    unique(
      as.numeric(
        cluster_values
      )
    )
  )
)

ep_levels <- paste0(
  "Ep",
  state_levels
)

if (length(ep_levels) != 8) {
  stop(
    "当前不是8个Ep states，实际：",
    length(ep_levels)
  )
}

ep_state_vector <- paste0(
  "Ep",
  cluster_values
)

names(ep_state_vector) <- names(
  mono_cluster
)

cat("\n========================================\n")
cat("CURRENT CLEAN Ep1-Ep8\n")
cat("========================================\n")

print(
  table(
    factor(
      ep_state_vector,
      levels = ep_levels
    )
  )
)

############################################################
## 6. Strict expected-count validation
############################################################

expected_counts <- c(
  Ep1 = 4796,
  Ep2 = 2674,
  Ep3 = 2358,
  Ep4 = 1681,
  Ep5 = 758,
  Ep6 = 204,
  Ep7 = 45,
  Ep8 = 43
)

observed_counts <- table(
  factor(
    ep_state_vector,
    levels = names(
      expected_counts
    )
  )
)

if (!all(
  as.integer(observed_counts) ==
  as.integer(expected_counts)
)) {
  
  cat("\nExpected:\n")
  print(expected_counts)
  
  cat("\nObserved:\n")
  print(observed_counts)
  
  stop(
    "Ep1-Ep8 cell counts与当前正式结果不一致。"
  )
}

############################################################
## 7. Load final Seurat object
############################################################

scRNA <- readRDS(
  seurat_file
)

if (!("RNA" %in% names(scRNA@assays))) {
  stop(
    "Seurat对象不存在RNA assay。"
  )
}

############################################################
## 8. Strict barcode matching
############################################################

clean_cells <- names(
  ep_state_vector
)

missing_cells <- setdiff(
  clean_cells,
  colnames(scRNA)
)

if (length(missing_cells) > 0) {
  stop(
    length(missing_cells),
    "个clean epithelial cells无法映射到Seurat对象。"
  )
}

if (length(clean_cells) != 12559) {
  stop(
    "clean epithelial barcode数不是12559。"
  )
}

############################################################
## 9. Extract exact CLEAN epithelial cells
############################################################

epi <- subset(
  scRNA,
  cells = clean_cells
)

if (ncol(epi) != 12559) {
  stop(
    "Seurat clean epithelial subset不是12559 cells。"
  )
}

############################################################
## 10. Assign Ep states by exact barcode
############################################################

ep_state_vector <- ep_state_vector[
  colnames(epi)
]

if (!identical(
  names(ep_state_vector),
  colnames(epi)
)) {
  stop(
    "Ep state与Seurat barcode顺序不一致。"
  )
}

epi$Ep_state <- factor(
  ep_state_vector,
  levels = ep_levels
)

cat("\n========================================\n")
cat("Ep states mapped to Seurat\n")
cat("========================================\n")

print(
  table(
    epi$Ep_state
  )
)

############################################################
## 11. Join Seurat v5 RNA layers
############################################################

rna_layers <- SeuratObject::Layers(
  epi[["RNA"]]
)

cat("\nRNA layers before JoinLayers:\n")
print(rna_layers)

count_layers <- grep(
  "^counts",
  rna_layers,
  value = TRUE
)

data_layers <- grep(
  "^data",
  rna_layers,
  value = TRUE
)

if (
  length(count_layers) > 1 ||
  length(data_layers) > 1
) {
  
  epi <- SeuratObject::JoinLayers(
    epi,
    assay = "RNA"
  )
}

############################################################
## 12. Ensure normalized data exists
############################################################

rna_layers <- SeuratObject::Layers(
  epi[["RNA"]]
)

if (!("data" %in% rna_layers)) {
  
  cat("\nNo unified data layer; running LogNormalize...\n")
  
  epi <- Seurat::NormalizeData(
    object = epi,
    assay = "RNA",
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    verbose = FALSE
  )
}

############################################################
## 13. Set identities
############################################################

Seurat::DefaultAssay(
  epi
) <- "RNA"

Seurat::Idents(
  epi
) <- "Ep_state"

cat("\nFinal identities:\n")
print(
  table(
    Seurat::Idents(
      epi
    )
  )
)

############################################################
## 14. Find positive markers for all Ep states
##
## Marker ranking for STATE ANNOTATION.
##
## Cell-level Wilcoxon P values are NOT interpreted
## as donor-level biological inference.
############################################################

set.seed(12345)

markers_all <- Seurat::FindAllMarkers(
  object = epi,
  assay = "RNA",
  test.use = "wilcox",
  only.pos = TRUE,
  min.pct = 0.10,
  logfc.threshold = 0.25,
  min.diff.pct = 0.05,
  return.thresh = 0.05,
  verbose = TRUE
)

if (nrow(markers_all) == 0) {
  stop(
    "FindAllMarkers没有返回marker。"
  )
}

############################################################
## 15. Identify FC column
############################################################

fc_candidates <- c(
  "avg_log2FC",
  "avg_logFC"
)

fc_col <- fc_candidates[
  fc_candidates %in%
    colnames(markers_all)
][1]

if (is.na(fc_col)) {
  stop(
    "没有找到avg_log2FC / avg_logFC。"
  )
}

############################################################
## 16. Standardize state order
############################################################

markers_all$cluster <- factor(
  as.character(
    markers_all$cluster
  ),
  levels = ep_levels
)

markers_all <- markers_all |>
  dplyr::arrange(
    cluster,
    dplyr::desc(
      .data[[fc_col]]
    )
  )

write.csv(
  markers_all,
  file.path(
    output_dir,
    "CLEAN_Ep1_Ep8_ALL_markers.csv"
  ),
  row.names = FALSE
)

############################################################
## 17. Build ROBUST annotation-marker set
##
## Requirements:
## pct.1 >= 0.20
## pct.1 - pct.2 >= 0.10
## positive FC
## adjusted P < 0.05
############################################################

robust_markers <- markers_all |>
  dplyr::mutate(
    pct_difference =
      pct.1 - pct.2
  ) |>
  dplyr::filter(
    p_val_adj < 0.05,
    pct.1 >= 0.20,
    pct_difference >= 0.10,
    .data[[fc_col]] > 0
  )

write.csv(
  robust_markers,
  file.path(
    output_dir,
    "CLEAN_Ep1_Ep8_ROBUST_markers.csv"
  ),
  row.names = FALSE
)

############################################################
## 18. Top20 robust markers
############################################################

top20 <- robust_markers |>
  dplyr::group_by(
    cluster
  ) |>
  dplyr::slice_max(
    order_by =
      .data[[fc_col]],
    n = 20,
    with_ties = FALSE
  ) |>
  dplyr::ungroup()

write.csv(
  top20,
  file.path(
    output_dir,
    "CLEAN_Ep1_Ep8_TOP20_markers.csv"
  ),
  row.names = FALSE
)

############################################################
## 19. Top10 robust markers
############################################################

top10 <- robust_markers |>
  dplyr::group_by(
    cluster
  ) |>
  dplyr::slice_max(
    order_by =
      .data[[fc_col]],
    n = 10,
    with_ties = FALSE
  ) |>
  dplyr::ungroup()

write.csv(
  top10,
  file.path(
    output_dir,
    "CLEAN_Ep1_Ep8_TOP10_markers.csv"
  ),
  row.names = FALSE
)

############################################################
## 20. Top10 summary
############################################################

top10_summary <- top10 |>
  dplyr::group_by(
    cluster
  ) |>
  dplyr::summarise(
    Top10_markers = paste(
      gene,
      collapse = ", "
    ),
    .groups = "drop"
  )

cat("\n====================================================\n")
cat("TOP10 ROBUST MARKERS PER Ep STATE\n")
cat("====================================================\n")

print(
  top10_summary,
  n = Inf
)

write.csv(
  top10_summary,
  file.path(
    output_dir,
    "CLEAN_Ep1_Ep8_TOP10_marker_summary.csv"
  ),
  row.names = FALSE
)

############################################################
## 21. Focus Ep2 / Ep6 Top30 markers
############################################################

focus_markers <- robust_markers |>
  dplyr::filter(
    cluster %in%
      c(
        "Ep2",
        "Ep6"
      )
  ) |>
  dplyr::group_by(
    cluster
  ) |>
  dplyr::slice_max(
    order_by =
      .data[[fc_col]],
    n = 30,
    with_ties = FALSE
  ) |>
  dplyr::ungroup()

cat("\n====================================================\n")
cat("Ep2 / Ep6 TOP30 MARKERS\n")
cat("====================================================\n")

print(
  focus_markers |>
    dplyr::select(
      cluster,
      gene,
      dplyr::all_of(
        fc_col
      ),
      pct.1,
      pct.2,
      p_val_adj
    ),
  n = Inf
)

write.csv(
  focus_markers,
  file.path(
    output_dir,
    "CLEAN_Ep2_Ep6_TOP30_markers.csv"
  ),
  row.names = FALSE
)

############################################################
## 22. Select Top3 markers/state for DotPlot
############################################################

top3 <- robust_markers |>
  dplyr::group_by(
    cluster
  ) |>
  dplyr::slice_max(
    order_by =
      .data[[fc_col]],
    n = 3,
    with_ties = FALSE
  ) |>
  dplyr::ungroup()

top3_genes <- unique(
  top3$gene
)

top3_genes <- top3_genes[
  top3_genes %in%
    rownames(epi)
]

############################################################
## 23. Top3-marker DotPlot
############################################################

p_top3 <- Seurat::DotPlot(
  object = epi,
  features = top3_genes,
  group.by = "Ep_state",
  assay = "RNA",
  dot.scale = 10,
  cols = c(
    "grey90",
    "#B2182B"
  )
) +
  
  ggplot2::labs(
    x = NULL,
    y = NULL,
    color = "Average expression",
    size = "Percent expressed"
  ) +
  
  ggplot2::theme_classic(
    base_size = 20
  ) +
  
  ggplot2::theme(
    axis.text.x =
      ggplot2::element_text(
        angle = 45,
        hjust = 1,
        vjust = 1,
        size = 16,
        face = "bold",
        color = "black"
      ),
    
    axis.text.y =
      ggplot2::element_text(
        size = 20,
        face = "bold",
        color = "black"
      ),
    
    legend.title =
      ggplot2::element_text(
        size = 17,
        face = "bold"
      ),
    
    legend.text =
      ggplot2::element_text(
        size = 15
      ),
    
    panel.grid =
      ggplot2::element_blank()
  )

ggplot2::ggsave(
  file.path(
    output_dir,
    "CLEAN_Ep1_Ep8_TOP3_marker_DotPlot.pdf"
  ),
  p_top3,
  width = 16,
  height = 7,
  device = cairo_pdf,
  bg = "white"
)

ggplot2::ggsave(
  file.path(
    output_dir,
    "CLEAN_Ep1_Ep8_TOP3_marker_DotPlot.png"
  ),
  p_top3,
  width = 16,
  height = 7,
  dpi = 600,
  bg = "white"
)

############################################################
## 24. CRC epithelial canonical marker panel
############################################################

canonical_markers <- c(
  
  ## General epithelial
  "EPCAM",
  "CDH1",
  "KRT8",
  "KRT18",
  "KRT19",
  "KRT20",
  
  ## Stem / progenitor
  "LGR5",
  "OLFM4",
  "SMOC2",
  "ASCL2",
  "SOX9",
  
  ## Cycling / proliferation
  "MKI67",
  "TOP2A",
  "UBE2C",
  "BIRC5",
  "CENPF",
  "CDK1",
  "CCNB1",
  "CCNA2",
  
  ## Mature absorptive / colonocyte
  "AQP8",
  "SLC26A3",
  "CA1",
  "CA2",
  "FABP1",
  "FABP2",
  "APOA1",
  "APOA4",
  "GUCA2A",
  "GUCA2B",
  
  ## BEST4 / OTOP2 absorptive
  "BEST4",
  "OTOP2",
  "CA7",
  "GUCA2A",
  "GUCA2B",
  
  ## Goblet / secretory
  "MUC2",
  "CLCA1",
  "SPINK4",
  "FCGBP",
  "TFF3",
  "REG4",
  "BEST2",
  "ITLN1",
  
  ## Enteroendocrine
  "CHGA",
  "CHGB",
  "PYY",
  "GCG",
  "INSL5",
  "PCSK1",
  "NEUROD1",
  "INSM1",
  
  ## Tuft-like
  "POU2F3",
  "TRPM5",
  "GNG13",
  "GNAT3",
  "IL17RB",
  
  ## Paneth-like
  "DEFA5",
  "DEFA6",
  "LYZ",
  
  ## Regenerative / tumor-associated epithelial
  "KRT17",
  "KRT80",
  "LAMC2",
  "TACSTD2",
  "CEACAM5",
  "CEACAM6",
  "CLDN4",
  "MUC17",
  "TFF2",
  "DUOX2",
  
  ## EMT / stress
  "VIM",
  "TGFBI",
  "KRT14",
  
  ## Inflammatory / IFN response
  "ISG15",
  "IFIT1",
  "IFIT3",
  "IFI6",
  "CXCL1",
  "CXCL8"
)

canonical_markers <- unique(
  canonical_markers
)

canonical_present <- canonical_markers[
  canonical_markers %in%
    rownames(epi)
]

canonical_missing <- setdiff(
  canonical_markers,
  canonical_present
)

cat(
  "\nCanonical markers present:",
  length(canonical_present),
  "\n"
)

if (length(canonical_missing) > 0) {
  
  cat(
    "Canonical markers absent:",
    paste(
      canonical_missing,
      collapse = ", "
    ),
    "\n"
  )
}

############################################################
## 25. Canonical-marker DotPlot
############################################################

p_canonical <- Seurat::DotPlot(
  object = epi,
  features = canonical_present,
  group.by = "Ep_state",
  assay = "RNA",
  dot.scale = 10,
  cols = c(
    "grey90",
    "#B2182B"
  )
) +
  
  ggplot2::labs(
    x = NULL,
    y = NULL,
    color = "Average expression",
    size = "Percent expressed"
  ) +
  
  ggplot2::theme_classic(
    base_size = 20
  ) +
  
  ggplot2::theme(
    axis.text.x =
      ggplot2::element_text(
        angle = 45,
        hjust = 1,
        vjust = 1,
        size = 14,
        face = "bold",
        color = "black"
      ),
    
    axis.text.y =
      ggplot2::element_text(
        size = 20,
        face = "bold",
        color = "black"
      ),
    
    legend.title =
      ggplot2::element_text(
        size = 17,
        face = "bold"
      ),
    
    legend.text =
      ggplot2::element_text(
        size = 15
      ),
    
    panel.grid =
      ggplot2::element_blank()
  )

ggplot2::ggsave(
  file.path(
    output_dir,
    "CLEAN_Ep1_Ep8_CANONICAL_marker_DotPlot.pdf"
  ),
  p_canonical,
  width = 29,
  height = 7,
  device = cairo_pdf,
  bg = "white"
)

ggplot2::ggsave(
  file.path(
    output_dir,
    "CLEAN_Ep1_Ep8_CANONICAL_marker_DotPlot.png"
  ),
  p_canonical,
  width = 29,
  height = 7,
  dpi = 600,
  bg = "white"
)

############################################################
## 26. Focused Ep2 / Ep6 canonical DotPlot
############################################################

epi_focus <- subset(
  epi,
  subset = Ep_state %in%
    c(
      "Ep2",
      "Ep6"
    )
)

epi_focus$Ep_state <- droplevels(
  epi_focus$Ep_state
)

p_focus <- Seurat::DotPlot(
  object = epi_focus,
  features = canonical_present,
  group.by = "Ep_state",
  assay = "RNA",
  dot.scale = 11,
  cols = c(
    "grey90",
    "#B2182B"
  )
) +
  
  ggplot2::labs(
    x = NULL,
    y = NULL,
    color = "Average expression",
    size = "Percent expressed"
  ) +
  
  ggplot2::theme_classic(
    base_size = 20
  ) +
  
  ggplot2::theme(
    axis.text.x =
      ggplot2::element_text(
        angle = 45,
        hjust = 1,
        vjust = 1,
        size = 14,
        face = "bold",
        color = "black"
      ),
    
    axis.text.y =
      ggplot2::element_text(
        size = 24,
        face = "bold",
        color = "black"
      ),
    
    legend.title =
      ggplot2::element_text(
        size = 17,
        face = "bold"
      ),
    
    legend.text =
      ggplot2::element_text(
        size = 15
      ),
    
    panel.grid =
      ggplot2::element_blank()
  )

ggplot2::ggsave(
  file.path(
    output_dir,
    "CLEAN_Ep2_Ep6_CANONICAL_marker_DotPlot.pdf"
  ),
  p_focus,
  width = 29,
  height = 4.5,
  device = cairo_pdf,
  bg = "white"
)

ggplot2::ggsave(
  file.path(
    output_dir,
    "CLEAN_Ep2_Ep6_CANONICAL_marker_DotPlot.png"
  ),
  p_focus,
  width = 29,
  height = 4.5,
  dpi = 600,
  bg = "white"
)

############################################################
## 27. Save cell-level Ep-state metadata
############################################################

cell_metadata_out <- data.frame(
  cell = colnames(epi),
  Ep_state = as.character(
    epi$Ep_state
  ),
  patient = as.character(
    epi$patient
  ),
  group = as.character(
    epi$group
  ),
  final_cluster = as.character(
    epi$final_cluster
  ),
  stringsAsFactors = FALSE
)

write.csv(
  cell_metadata_out,
  file.path(
    output_dir,
    "CLEAN_Ep1_Ep8_cell_metadata.csv"
  ),
  row.names = FALSE
)

############################################################
## 28. Save epithelial Seurat object
############################################################

saveRDS(
  epi,
  file.path(
    output_dir,
    "GSE294300_CLEAN_Epithelial_with_Ep1_Ep8.rds"
  ),
  compress = FALSE
)

############################################################
## 29. Session info
############################################################

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    output_dir,
    "sessionInfo.txt"
  )
)

############################################################
## 30. Finished
############################################################

cat("\n############################################################\n")
cat("## CLEAN Ep1-Ep8 MARKER ANALYSIS COMPLETED\n")
cat("############################################################\n")

cat(
  "\nClean epithelial cells:",
  ncol(epi),
  "\n"
)

cat(
  "States:",
  paste(
    ep_levels,
    collapse = ", "
  ),
  "\n"
)

cat(
  "\nMost important files:\n"
)

cat(
  "1. CLEAN_Ep1_Ep8_TOP10_marker_summary.csv\n"
)

cat(
  "2. CLEAN_Ep1_Ep8_TOP20_markers.csv\n"
)

cat(
  "3. CLEAN_Ep2_Ep6_TOP30_markers.csv\n"
)

cat(
  "4. CLEAN_Ep1_Ep8_TOP3_marker_DotPlot.png\n"
)

cat(
  "5. CLEAN_Ep1_Ep8_CANONICAL_marker_DotPlot.png\n"
)

cat(
  "6. CLEAN_Ep2_Ep6_CANONICAL_marker_DotPlot.png\n"
)

cat(
  "\nOutput directory:\n",
  output_dir,
  "\n"
)