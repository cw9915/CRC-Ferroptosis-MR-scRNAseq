
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list = ls())

############################################################
## 0. Paths
############################################################

fig1_dir   <- "C:/Users/92066/Desktop/单细胞孟德尔/Figure1"
input_dir  <- "C:/Users/92066/Desktop/单细胞孟德尔/Figure6"
output_dir <- file.path(input_dir, "28")

ssgsea_csv <- file.path(
  output_dir,
  "03_ssGSEA_28cell_normalized_scores.csv"
)

meta_csv <- file.path(
  fig1_dir,
  "GSE44861_group_final.csv"
)

if (!dir.exists(output_dir)) {
  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

############################################################
## 1. Required packages
############################################################

required_pkgs <- c(
  "ggplot2",
  "dplyr",
  "tidyr",
  "RColorBrewer"
)

missing_pkgs <- required_pkgs[
  !vapply(
    required_pkgs,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_pkgs) > 0) {
  stop(
    "Missing R packages: ",
    paste(missing_pkgs, collapse = ", ")
  )
}

############################################################
## 2. Check input files
############################################################

if (!file.exists(ssgsea_csv)) {
  stop(
    "Normalized ssGSEA score file not found:\n",
    ssgsea_csv,
    "\n\nPlease run the main 28-cell ssGSEA analysis first."
  )
}

if (!file.exists(meta_csv)) {
  stop(
    "Metadata file not found:\n",
    meta_csv
  )
}

############################################################
## 3. Read normalized ssGSEA matrix
##
## Expected structure:
## rows    = 28 immune-cell types
## columns = samples
############################################################

ssgsea_scores <- read.csv(
  ssgsea_csv,
  row.names = 1,
  check.names = FALSE
)

ssgsea_scores <- as.matrix(ssgsea_scores)
mode(ssgsea_scores) <- "numeric"

cat("\n========================================\n")
cat("ssGSEA matrix\n")
cat("========================================\n")
cat(
  nrow(ssgsea_scores),
  "immune-cell types x",
  ncol(ssgsea_scores),
  "samples\n"
)

if (nrow(ssgsea_scores) != 28) {
  stop(
    "The ssGSEA matrix contains ",
    nrow(ssgsea_scores),
    " cell types rather than 28.\n",
    "Please check the upstream analysis."
  )
}

############################################################
## 4. Use EXACT cell-type names from the matrix
############################################################

cell_types_28 <- rownames(ssgsea_scores)

cat("\nActual 28 cell-type names in the ssGSEA matrix:\n")
print(cell_types_28)

############################################################
## 5. Clean names ONLY for legend display
##
## Underlying ssGSEA matrix names and results are NOT changed.
############################################################

cell_display_labels <- cell_types_28

## Fix known typo from Table S6-derived naming
cell_display_labels <- gsub(
  "memeory",
  "memory",
  cell_display_labels,
  ignore.case = TRUE
)

## Standardize plural forms only for display
display_replacements <- c(
  "^Immature B cells$" = "Immature B cell",
  "^Activated B cells$" = "Activated B cell",
  "^Activated CD4 T cells$" = "Activated CD4 T cell",
  "^Activated CD8 T cells$" = "Activated CD8 T cell",
  "^Activated dendritic cells$" = "Activated dendritic cell",
  "^Memory B cells$" = "Memory B cell",
  "^Macrophages$" = "Macrophage",
  "^Monocytes$" = "Monocyte",
  "^Neutrophils$" = "Neutrophil",
  "^Natural killer cells$" = "Natural killer cell",
  "^Natural killer T cells$" = "Natural killer T cell",
  "^Regulatory T cells$" = "Regulatory T cell",
  "^T follicular helper cells$" = "T follicular helper cell",
  "^Type 1 T helper cells$" = "Type 1 T helper cell",
  "^Type 17 T helper cells$" = "Type 17 T helper cell",
  "^Type 2 T helper cells$" = "Type 2 T helper cell"
)

for (pat in names(display_replacements)) {
  cell_display_labels <- gsub(
    pat,
    display_replacements[[pat]],
    cell_display_labels,
    ignore.case = TRUE
  )
}

cell_display_labels <- gsub(
  "\\s+",
  " ",
  cell_display_labels
)

cell_display_labels <- trimws(cell_display_labels)

names(cell_display_labels) <- cell_types_28

cat("\nDisplay labels used in the legend:\n")
print(cell_display_labels)

############################################################
## 6. Read metadata
############################################################

meta <- read.csv(
  meta_csv,
  row.names = 1,
  check.names = FALSE
)

group_candidates <- c(
  "diagnosis",
  "group", "Group", "GROUP",
  "type", "Type",
  "sample_type", "Sample_type",
  "tissue", "Tissue",
  "condition", "Condition"
)

group_hit <- group_candidates[
  group_candidates %in% colnames(meta)
]

if (length(group_hit) == 0) {
  stop(
    "No group column was found in metadata."
  )
}

group_col <- group_hit[1]

meta$diagnosis <- as.character(
  meta[[group_col]]
)

############################################################
## 7. Standardize group labels
############################################################

meta$diagnosis[
  meta$diagnosis %in% c(
    "control", "Control",
    "normal", "Normal",
    "N", "n",
    "Adjacent", "adjacent",
    "Non-tumor", "Non_tumor"
  )
] <- "Normal"

meta$diagnosis[
  meta$diagnosis %in% c(
    "case", "Case",
    "tumor", "Tumor",
    "T", "t",
    "CRC",
    "Colorectal cancer",
    "Cancer"
  )
] <- "Tumor"

valid_groups <- na.omit(
  unique(meta$diagnosis)
)

if (!all(valid_groups %in% c("Normal", "Tumor"))) {
  cat("\nUnrecognized group labels:\n")
  print(
    table(
      meta$diagnosis,
      useNA = "ifany"
    )
  )
  stop(
    "Please correct the group labels in metadata."
  )
}

############################################################
## 8. Match ssGSEA samples with metadata
##
## Keep original ssGSEA sample order.
## DO NOT sort samples by score.
############################################################

original_samples <- colnames(
  ssgsea_scores
)

common_samples <- original_samples[
  original_samples %in% rownames(meta)
]

if (length(common_samples) == 0) {
  stop(
    "No common samples were found between ssGSEA scores and metadata."
  )
}

common_samples <- common_samples[
  meta[
    common_samples,
    "diagnosis"
  ] %in% c("Normal", "Tumor")
]

############################################################
## 9. Normal first, Tumor second;
##    preserve original order WITHIN each group
############################################################

normal_samples <- common_samples[
  meta[
    common_samples,
    "diagnosis"
  ] == "Normal"
]

tumor_samples <- common_samples[
  meta[
    common_samples,
    "diagnosis"
  ] == "Tumor"
]

sample_order <- c(
  normal_samples,
  tumor_samples
)

ssgsea_scores <- ssgsea_scores[
  ,
  sample_order,
  drop = FALSE
]

cat("\nSamples retained:\n")
cat(
  "Normal:",
  length(normal_samples),
  "\n"
)
cat(
  "Tumor :",
  length(tumor_samples),
  "\n"
)
cat(
  "Total :",
  length(sample_order),
  "\n"
)

############################################################
## 10. Convert to long format
############################################################

A_df <- as.data.frame(
  t(ssgsea_scores),
  check.names = FALSE
)

A_df$Sample <- rownames(A_df)

group_df <- data.frame(
  Sample = sample_order,
  Group = factor(
    c(
      rep(
        "Normal",
        length(normal_samples)
      ),
      rep(
        "Tumor",
        length(tumor_samples)
      )
    ),
    levels = c(
      "Normal",
      "Tumor"
    )
  ),
  stringsAsFactors = FALSE
)

A_long <- A_df |>
  dplyr::left_join(
    group_df,
    by = "Sample"
  ) |>
  tidyr::pivot_longer(
    cols = dplyr::all_of(
      cell_types_28
    ),
    names_to = "CellType",
    values_to = "Score"
  ) |>
  dplyr::mutate(
    Sample = factor(
      Sample,
      levels = sample_order
    ),
    Group = factor(
      Group,
      levels = c(
        "Normal",
        "Tumor"
      )
    ),
    CellType = factor(
      CellType,
      levels = cell_types_28
    )
  )

############################################################
## 11. Save plot data
############################################################

write.csv(
  A_long,
  file.path(
    output_dir,
    "Figure6A_plot_data_compact_maintext.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)

############################################################
## 12. Diagnostics
############################################################

stack_diag <- A_long |>
  dplyr::group_by(
    Sample,
    Group
  ) |>
  dplyr::summarise(
    Positive_stack = sum(
      Score[
        Score > 0
      ],
      na.rm = TRUE
    ),
    Negative_stack = sum(
      Score[
        Score < 0
      ],
      na.rm = TRUE
    ),
    .groups = "drop"
  )

cat("\nIndividual ssGSEA score range:\n")
print(
  range(
    A_long$Score,
    na.rm = TRUE
  )
)

cat("\nStacked range:\n")
cat(
  "Negative minimum:",
  min(
    stack_diag$Negative_stack,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Positive maximum:",
  max(
    stack_diag$Positive_stack,
    na.rm = TRUE
  ),
  "\n"
)

############################################################
## 13. Colors
############################################################

cell_colors <- grDevices::colorRampPalette(
  RColorBrewer::brewer.pal(
    12,
    "Set3"
  )
)(
  length(cell_types_28)
)

names(cell_colors) <- cell_types_28

############################################################
## 14. Figure 6A
##
## LEFT SIDE:
## Same geometry/style as previous robust version.
##
## RIGHT SIDE:
## Larger one-column legend.
############################################################

pA <- ggplot2::ggplot(
  A_long,
  ggplot2::aes(
    x = Sample,
    y = Score,
    fill = CellType
  )
) +
  
  ##########################################################
## LEFT plot: unchanged
##########################################################
ggplot2::geom_col(
  width = 0.78,
  position = "stack",
  linewidth = 0
) +
  
  ggplot2::facet_grid(
    . ~ Group,
    scales = "free_x",
    space = "free_x",
    switch = "x"
  ) +
  
  ggplot2::scale_fill_manual(
    values = cell_colors,
    breaks = cell_types_28,
    labels = cell_display_labels[
      cell_types_28
    ],
    drop = FALSE
  ) +
  
  ggplot2::scale_x_discrete(
    expand = c(
      0,
      0
    )
  ) +
  
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(
      mult = c(
        0.02,
        0.035
      )
    )
  ) +
  
  ggplot2::labs(
    x = NULL,
    y = "ssGSEA score",
    fill = "Cell Type"
  ) +
  
  ##########################################################
## RIGHT legend: one vertical column
##########################################################
ggplot2::guides(
  fill = ggplot2::guide_legend(
    ncol = 1,
    byrow = FALSE,
    title.position = "top",
    keyheight = grid::unit(
      0.27,
      "cm"
    ),
    keywidth = grid::unit(
      0.38,
      "cm"
    )
  )
) +
  
  ggplot2::theme_classic(
    base_size = 16
  ) +
  
  ggplot2::theme(
    
    ########################################################
    ## LEFT plot settings — unchanged
    ########################################################
    axis.text.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank(),
    
    axis.title.y = ggplot2::element_text(
      size = 19,
      face = "bold",
      color = "black",
      margin = ggplot2::margin(
        r = 8
      )
    ),
    
    axis.text.y = ggplot2::element_text(
      size = 14,
      color = "black"
    ),
    
    axis.line = ggplot2::element_line(
      color = "black",
      linewidth = 0.75
    ),
    
    strip.placement = "outside",
    
    strip.background = ggplot2::element_rect(
      fill = "grey94",
      color = "grey50",
      linewidth = 0.55
    ),
    
    strip.text.x.bottom = ggplot2::element_text(
      size = 17,
      face = "bold",
      color = "black",
      margin = ggplot2::margin(
        t = 3,
        b = 3
      )
    ),
    
    panel.spacing.x = grid::unit(
      0.08,
      "cm"
    ),
    
    ########################################################
    ## RIGHT legend — enlarged
    ########################################################
    legend.position = "right",
    legend.direction = "vertical",
    legend.box = "vertical",
    
    legend.title = ggplot2::element_text(
      size = 15,
      face = "bold",
      color = "black",
      margin = ggplot2::margin(
        b = 3
      )
    ),
    
    legend.text = ggplot2::element_text(
      size = 10.5,
      color = "black",
      lineheight = 0.82,
      margin = ggplot2::margin(
        t = 0,
        b = 0
      )
    ),
    
    legend.key.height = grid::unit(
      0.27,
      "cm"
    ),
    
    legend.key.width = grid::unit(
      0.38,
      "cm"
    ),
    
    ## Small gap between adjacent legend color squares
    legend.key.spacing.y = grid::unit(
      0.035,
      "cm"
    ),
    
    legend.spacing.y = grid::unit(
      0.00,
      "cm"
    ),
    
    legend.box.spacing = grid::unit(
      0.35,
      "cm"
    ),
    
    ########################################################
    ## Margin
    ########################################################
    plot.margin = ggplot2::margin(
      t = 8,
      r = 8,
      b = 7,
      l = 8
    )
  )

############################################################
## 15. Display
############################################################

print(
  pA
)

############################################################
## 16. Save figure
##
## IMPORTANT:
## Compact main-text layout:
## height is reduced to 7.2 inches and width set to 19 inches.
## The left panel data/order/geometry are unchanged; only the
## overall canvas and legend are compacted for multi-panel figures.
############################################################

ggplot2::ggsave(
  filename = file.path(
    output_dir,
    "Figure6A_28immune_ssGSEA_compact_maintext.pdf"
  ),
  plot = pA,
  width = 19,
  height = 7.2,
  limitsize = FALSE
)

ggplot2::ggsave(
  filename = file.path(
    output_dir,
    "Figure6A_28immune_ssGSEA_compact_maintext.png"
  ),
  plot = pA,
  width = 19,
  height = 7.2,
  dpi = 600,
  limitsize = FALSE
)

############################################################
## 17. Save cell-name mapping
############################################################

cell_name_mapping <- data.frame(
  Original_name = cell_types_28,
  Display_name = unname(
    cell_display_labels[
      cell_types_28
    ]
  ),
  stringsAsFactors = FALSE
)

write.csv(
  cell_name_mapping,
  file.path(
    output_dir,
    "Figure6A_cell_name_mapping_compact_maintext.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)

cat("\n============================================================\n")
cat("Figure 6A completed successfully.\n")
cat("LEFT plot geometry/order kept unchanged.\n")
cat("RIGHT legend kept in one compact vertical column.\n")
cat("Output height reduced to 7.2 inches for main-text composition.\n")
cat("Output width set to 19 inches.\n")
cat("============================================================\n")
############################################################
## Figure 6B-C ONLY
## 28 immune-cell ssGSEA analysis for GSE44861
##
## Figure 6B:
##   - Complete matched tumor-normal pairs
##   - Paired Wilcoxon signed-rank test
##   - BH correction across all 28 immune-cell signatures
##   - All fonts bold
##   - Enlarged x/y-axis fonts
##
## Figure 6C:
##   - AKR1C2 / MAPK3 / MAP3K11 vs significant immune cells
##   - Spearman correlation
##   - BH correction across displayed gene-cell tests
##   - All fonts bold
##   - Enlarged x/y-axis fonts
##
## IMPORTANT:
##   - Uses normalized ssGSEA scores directly
##   - NO extra row-wise Z-score
##   - Uses actual row names from the ssGSEA matrix
##   - Only cleans cell-type names for display
############################################################

Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list = ls())

############################################################
## 0. Paths
############################################################

fig1_dir   <- "C:/Users/92066/Desktop/单细胞孟德尔/Figure1"
input_dir  <- "C:/Users/92066/Desktop/单细胞孟德尔/Figure6"
output_dir <- file.path(input_dir, "28")

ssgsea_csv <- file.path(
  output_dir,
  "03_ssGSEA_28cell_normalized_scores.csv"
)

meta_csv <- file.path(
  fig1_dir,
  "GSE44861_group_final.csv"
)

expr_csv <- file.path(
  fig1_dir,
  "GSE44861_gene_symbol_expression.csv"
)

pair_map_csv <- file.path(
  output_dir,
  "05_GSE44861_GEO_pair_mapping.csv"
)

target_genes <- c(
  "AKR1C2",
  "MAPK3",
  "MAP3K11"
)

if (!dir.exists(output_dir)) {
  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

############################################################
## 1. Required packages
############################################################

required_pkgs <- c(
  "ggplot2",
  "dplyr",
  "tidyr",
  "limma"
)

missing_pkgs <- required_pkgs[
  !vapply(
    required_pkgs,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_pkgs) > 0) {
  stop(
    "Missing R packages: ",
    paste(missing_pkgs, collapse = ", ")
  )
}

############################################################
## 2. Helper functions
############################################################

p_to_star <- function(p) {
  ifelse(
    is.na(p), "",
    ifelse(
      p < 0.001, "***",
      ifelse(
        p < 0.01, "**",
        ifelse(
          p < 0.05, "*", "ns"
        )
      )
    )
  )
}

derive_pair_id <- function(title) {
  
  x <- toupper(
    trimws(
      as.character(title)
    )
  )
  
  x <- gsub(
    "\\s+",
    "",
    x
  )
  
  x <- sub(
    "[-_]?[NT]$",
    "",
    x
  )
  
  x
}

derive_group_from_title <- function(title) {
  
  x <- toupper(
    trimws(
      as.character(title)
    )
  )
  
  x <- gsub(
    "\\s+",
    "",
    x
  )
  
  out <- rep(
    NA_character_,
    length(x)
  )
  
  out[
    grepl(
      "N$",
      x
    )
  ] <- "Normal"
  
  out[
    grepl(
      "T$",
      x
    )
  ] <- "Tumor"
  
  out
}

############################################################
## 3. Check input files
############################################################

if (!file.exists(ssgsea_csv)) {
  stop(
    "Normalized ssGSEA score file not found:\n",
    ssgsea_csv
  )
}

if (!file.exists(meta_csv)) {
  stop(
    "Metadata file not found:\n",
    meta_csv
  )
}

if (!file.exists(expr_csv)) {
  stop(
    "Expression file not found:\n",
    expr_csv
  )
}

############################################################
## 4. Read normalized ssGSEA matrix
############################################################

ssgsea_scores <- read.csv(
  ssgsea_csv,
  row.names = 1,
  check.names = FALSE
)

ssgsea_scores <- as.matrix(
  ssgsea_scores
)

mode(
  ssgsea_scores
) <- "numeric"

if (nrow(ssgsea_scores) != 28) {
  stop(
    "The ssGSEA matrix contains ",
    nrow(ssgsea_scores),
    " cell types rather than 28."
  )
}

cell_types_28 <- rownames(
  ssgsea_scores
)

cat("\nActual immune-cell names:\n")
print(
  cell_types_28
)

############################################################
## 5. Clean cell-type names ONLY for display
############################################################

cell_display_labels <- cell_types_28

## Known typo
cell_display_labels <- gsub(
  "memeory",
  "memory",
  cell_display_labels,
  ignore.case = TRUE
)

## Common plural forms
display_replacements <- c(
  "^Immature B cells$" = "Immature B cell",
  "^Activated B cells$" = "Activated B cell",
  "^Activated CD4 T cells$" = "Activated CD4 T cell",
  "^Activated CD8 T cells$" = "Activated CD8 T cell",
  "^Activated dendritic cells$" = "Activated dendritic cell",
  "^Memory B cells$" = "Memory B cell",
  "^Macrophages$" = "Macrophage",
  "^Monocytes$" = "Monocyte",
  "^Neutrophils$" = "Neutrophil",
  "^Natural killer cells$" = "Natural killer cell",
  "^Natural killer T cells$" = "Natural killer T cell",
  "^Regulatory T cells$" = "Regulatory T cell",
  "^T follicular helper cells$" = "T follicular helper cell",
  "^Type 1 T helper cells$" = "Type 1 T helper cell",
  "^Type 17 T helper cells$" = "Type 17 T helper cell",
  "^Type 2 T helper cells$" = "Type 2 T helper cell"
)

for (pat in names(display_replacements)) {
  
  cell_display_labels <- gsub(
    pat,
    display_replacements[[pat]],
    cell_display_labels,
    ignore.case = TRUE
  )
}

cell_display_labels <- gsub(
  "\\s+",
  " ",
  cell_display_labels
)

cell_display_labels <- trimws(
  cell_display_labels
)

names(
  cell_display_labels
) <- cell_types_28

############################################################
## 6. Read metadata
############################################################

meta <- read.csv(
  meta_csv,
  row.names = 1,
  check.names = FALSE
)

group_candidates <- c(
  "diagnosis",
  "group", "Group", "GROUP",
  "type", "Type",
  "sample_type", "Sample_type",
  "tissue", "Tissue",
  "condition", "Condition"
)

group_hit <- group_candidates[
  group_candidates %in% colnames(meta)
]

if (length(group_hit) == 0) {
  stop(
    "No Normal/Tumor group column was found in metadata."
  )
}

group_col <- group_hit[1]

meta$diagnosis <- as.character(
  meta[[group_col]]
)

meta$diagnosis[
  meta$diagnosis %in% c(
    "control", "Control",
    "normal", "Normal",
    "N", "n",
    "Adjacent", "adjacent",
    "Non-tumor", "Non_tumor"
  )
] <- "Normal"

meta$diagnosis[
  meta$diagnosis %in% c(
    "case", "Case",
    "tumor", "Tumor",
    "T", "t",
    "CRC",
    "Colorectal cancer",
    "Cancer"
  )
] <- "Tumor"

if (!all(
  na.omit(
    unique(
      meta$diagnosis
    )
  ) %in% c(
    "Normal",
    "Tumor"
  )
)) {
  stop(
    "Unrecognized Normal/Tumor labels remain in metadata."
  )
}

############################################################
## 7. Determine PairID — FIXED VERSION
##
## IMPORTANT BUG FIX:
## Do NOT create meta$PairID before checking whether a real PairID
## column exists. In the previous script, creating an empty PairID
## column first caused the code to mistakenly detect that empty column
## and therefore return 0 matched pairs.
##
## Priority:
## 1) existing REAL pair/patient column in local metadata
## 2) previously generated GEO pair mapping
## 3) local title/sample-title column
## 4) GEOquery fallback
############################################################

## Save the original metadata column names BEFORE adding PairID
original_meta_columns <- colnames(meta)

pair_candidates <- c(
  "pair_id",
  "PairID",
  "pairID",
  "Pair_ID",
  "pair",
  "patient_id",
  "PatientID",
  "patientID",
  "Patient_ID",
  "patient",
  "subject_id",
  "SubjectID",
  "subjectID",
  "Subject_ID",
  "subject",
  "donor_id",
  "DonorID",
  "donor"
)

## Only search columns that truly existed in the input metadata
pair_hit <- pair_candidates[
  pair_candidates %in% original_meta_columns
]

title_candidates <- c(
  "title",
  "Title",
  "sample_title",
  "Sample_title",
  "SampleTitle",
  "geo_title",
  "GEO_title"
)

title_hit <- title_candidates[
  title_candidates %in% original_meta_columns
]

## NOW create PairID
meta$PairID <- NA_character_

pair_source <- NULL

############################################################
## 7A. Existing real pair/patient column
############################################################

if (length(pair_hit) > 0) {
  
  pair_col <- pair_hit[1]
  
  meta$PairID <- trimws(
    as.character(
      meta[[pair_col]]
    )
  )
  
  meta$PairID[
    meta$PairID %in% c(
      "",
      "NA",
      "NaN",
      "NULL"
    )
  ] <- NA_character_
  
  pair_source <- paste0(
    "local metadata column: ",
    pair_col
  )
}

############################################################
## 7B. Previously generated GEO pair mapping
############################################################

if (
  all(is.na(meta$PairID)) &&
  file.exists(pair_map_csv)
) {
  
  cat(
    "\nNo usable local PairID column found.\n",
    "Trying saved pair mapping:\n",
    pair_map_csv,
    "\n"
  )
  
  pair_map <- read.csv(
    pair_map_csv,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  if (
    all(
      c(
        "Sample",
        "PairID"
      ) %in% colnames(pair_map)
    )
  ) {
    
    pair_map$Sample <- trimws(
      as.character(
        pair_map$Sample
      )
    )
    
    pair_map$PairID <- trimws(
      as.character(
        pair_map$PairID
      )
    )
    
    pair_map$PairID[
      pair_map$PairID %in% c(
        "",
        "NA",
        "NaN",
        "NULL"
      )
    ] <- NA_character_
    
    pair_lookup <- pair_map$PairID
    names(pair_lookup) <- pair_map$Sample
    
    common_pair_samples <- intersect(
      rownames(meta),
      names(pair_lookup)
    )
    
    meta[
      common_pair_samples,
      "PairID"
    ] <- pair_lookup[
      common_pair_samples
    ]
    
    if (sum(!is.na(meta$PairID)) > 0) {
      pair_source <- "saved GEO pair mapping"
    }
  }
}

############################################################
## 7C. Derive PairID from a local sample-title column
############################################################

if (
  all(is.na(meta$PairID)) &&
  length(title_hit) > 0
) {
  
  title_col <- title_hit[1]
  
  meta$PairID <- derive_pair_id(
    meta[[title_col]]
  )
  
  meta$PairID[
    is.na(meta[[title_col]]) |
      trimws(
        as.character(
          meta[[title_col]]
        )
      ) == ""
  ] <- NA_character_
  
  pair_source <- paste0(
    "local sample-title column: ",
    title_col
  )
}

############################################################
## 7D. GEOquery fallback
############################################################

if (all(is.na(meta$PairID))) {
  
  cat(
    "\nNo usable PairID was found locally.\n",
    "Retrieving official GSE44861 sample titles with GEOquery...\n"
  )
  
  if (!requireNamespace(
    "GEOquery",
    quietly = TRUE
  )) {
    
    stop(
      "No PairID could be identified.\n\n",
      "Install GEOquery and rerun:\n",
      "if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager')\n",
      "BiocManager::install('GEOquery')"
    )
  }
  
  gse_obj <- GEOquery::getGEO(
    "GSE44861",
    GSEMatrix = TRUE,
    getGPL = FALSE
  )
  
  if (length(gse_obj) == 0) {
    stop(
      "GEOquery failed to retrieve GSE44861."
    )
  }
  
  geo_pd <- Biobase::pData(
    gse_obj[[1]]
  )
  
  geo_map <- data.frame(
    Sample = rownames(
      geo_pd
    ),
    GEO_title = as.character(
      geo_pd$title
    ),
    stringsAsFactors = FALSE
  )
  
  geo_map$PairID <- derive_pair_id(
    geo_map$GEO_title
  )
  
  geo_map$GEO_group <- derive_group_from_title(
    geo_map$GEO_title
  )
  
  pair_lookup <- geo_map$PairID
  names(pair_lookup) <- geo_map$Sample
  
  common_pair_samples <- intersect(
    rownames(meta),
    names(pair_lookup)
  )
  
  meta[
    common_pair_samples,
    "PairID"
  ] <- pair_lookup[
    common_pair_samples
  ]
  
  pair_source <- "GEOquery / official GSE44861 titles"
  
  write.csv(
    data.frame(
      Sample = rownames(meta),
      PairID = meta$PairID,
      Local_group = meta$diagnosis,
      stringsAsFactors = FALSE
    ),
    pair_map_csv,
    row.names = FALSE,
    quote = FALSE
  )
}

############################################################
## 7E. PairID diagnostics
############################################################

cat(
  "\nPairID source:",
  pair_source,
  "\n"
)

cat(
  "Samples with non-missing PairID:",
  sum(
    !is.na(
      meta$PairID
    )
  ),
  "/",
  nrow(meta),
  "\n"
)

cat(
  "Unique PairIDs:",
  length(
    unique(
      na.omit(
        meta$PairID
      )
    )
  ),
  "\n"
)

pair_debug <- data.frame(
  Sample = rownames(meta),
  Group = meta$diagnosis,
  PairID = meta$PairID,
  stringsAsFactors = FALSE
)

write.csv(
  pair_debug,
  file.path(
    output_dir,
    "Figure6BC_PairID_debug.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)

############################################################
## 8. Restrict to samples available in the ssGSEA matrix
############################################################

available_samples <- intersect(
  colnames(ssgsea_scores),
  rownames(meta)
)

meta_use <- meta[
  available_samples,
  ,
  drop = FALSE
]

meta_use$Sample <- rownames(
  meta_use
)

############################################################
## 9. Identify complete 1:1 Normal-Tumor pairs
############################################################

pair_table <- meta_use |>
  dplyr::mutate(
    PairID = as.character(
      PairID
    )
  ) |>
  dplyr::filter(
    !is.na(PairID),
    PairID != "",
    diagnosis %in% c(
      "Normal",
      "Tumor"
    )
  ) |>
  dplyr::count(
    PairID,
    diagnosis,
    name = "n"
  ) |>
  tidyr::pivot_wider(
    names_from = diagnosis,
    values_from = n,
    values_fill = 0
  )

if (!"Normal" %in% colnames(pair_table)) {
  pair_table$Normal <- 0
}

if (!"Tumor" %in% colnames(pair_table)) {
  pair_table$Tumor <- 0
}

valid_pairs <- pair_table |>
  dplyr::filter(
    Normal == 1,
    Tumor == 1
  ) |>
  dplyr::pull(
    PairID
  )

cat(
  "\nComplete matched pairs:",
  length(valid_pairs),
  "\n"
)

if (length(valid_pairs) < 10) {
  
  cat("\nPair table preview:\n")
  print(
    utils::head(
      pair_table,
      20
    )
  )
  
  stop(
    "Too few complete matched pairs were identified (",
    length(valid_pairs),
    ").\n",
    "Please inspect:\n",
    file.path(
      output_dir,
      "Figure6BC_PairID_debug.csv"
    ),
    "\n\nExpected approximately 47 complete matched tumor-normal pairs for GSE44861."
  )
}

if (length(valid_pairs) != 47) {
  warning(
    "The code identified ",
    length(valid_pairs),
    " matched pairs rather than 47. ",
    "Please verify the pair mapping."
  )
}

############################################################
## 10. Build sample-level ssGSEA table
############################################################

score_df <- as.data.frame(
  t(
    ssgsea_scores
  ),
  check.names = FALSE
)

score_df$Sample <- rownames(
  score_df
)

score_df$Group <- meta[
  score_df$Sample,
  "diagnosis"
]

score_df$PairID <- meta[
  score_df$Sample,
  "PairID"
]

############################################################
## 11. Figure 6B statistics
##     Paired Wilcoxon signed-rank test
############################################################

B_stats_list <- lapply(
  cell_types_28,
  function(ct) {
    
    tmp <- score_df |>
      dplyr::filter(
        PairID %in% valid_pairs
      ) |>
      dplyr::transmute(
        PairID = PairID,
        Group = Group,
        Score = .data[[ct]]
      ) |>
      tidyr::pivot_wider(
        names_from = Group,
        values_from = Score
      ) |>
      dplyr::filter(
        !is.na(Normal),
        !is.na(Tumor)
      )
    
    if (nrow(tmp) < 3) {
      
      return(
        data.frame(
          CellType = ct,
          N_pairs = nrow(tmp),
          Normal_median = NA_real_,
          Tumor_median = NA_real_,
          Median_paired_difference = NA_real_,
          P_value = NA_real_,
          stringsAsFactors = FALSE
        )
      )
    }
    
    wt <- stats::wilcox.test(
      x = tmp$Tumor,
      y = tmp$Normal,
      paired = TRUE,
      exact = FALSE
    )
    
    data.frame(
      CellType = ct,
      N_pairs = nrow(tmp),
      Normal_median = median(
        tmp$Normal,
        na.rm = TRUE
      ),
      Tumor_median = median(
        tmp$Tumor,
        na.rm = TRUE
      ),
      Median_paired_difference = median(
        tmp$Tumor - tmp$Normal,
        na.rm = TRUE
      ),
      P_value = wt$p.value,
      stringsAsFactors = FALSE
    )
  }
)

B_stats <- dplyr::bind_rows(
  B_stats_list
) |>
  dplyr::mutate(
    FDR = stats::p.adjust(
      P_value,
      method = "BH"
    ),
    Direction = dplyr::case_when(
      Median_paired_difference > 0 ~ "Higher in Tumor",
      Median_paired_difference < 0 ~ "Lower in Tumor",
      TRUE ~ "No change"
    ),
    Signif = p_to_star(
      FDR
    )
  ) |>
  dplyr::arrange(
    FDR,
    P_value
  )

write.csv(
  B_stats,
  file.path(
    output_dir,
    "Figure6B_paired_Wilcoxon_28cell_results_BOLD.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)

############################################################
## 12. Figure 6B plot data
############################################################

B_long <- score_df |>
  dplyr::filter(
    PairID %in% valid_pairs
  ) |>
  dplyr::select(
    Sample,
    PairID,
    Group,
    dplyr::all_of(
      cell_types_28
    )
  ) |>
  tidyr::pivot_longer(
    cols = dplyr::all_of(
      cell_types_28
    ),
    names_to = "CellType",
    values_to = "Score"
  ) |>
  dplyr::mutate(
    Group = factor(
      Group,
      levels = c(
        "Normal",
        "Tumor"
      )
    ),
    CellType = factor(
      CellType,
      levels = cell_types_28
    )
  )

star_y <- B_long |>
  dplyr::group_by(
    CellType
  ) |>
  dplyr::summarise(
    ymax = max(
      Score,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

global_range <- diff(
  range(
    B_long$Score,
    na.rm = TRUE
  )
)

star_y$y <- star_y$ymax +
  0.06 * global_range

star_y <- star_y |>
  dplyr::left_join(
    B_stats |>
      dplyr::select(
        CellType,
        FDR,
        Signif
      ),
    by = "CellType"
  )

write.csv(
  B_long,
  file.path(
    output_dir,
    "Figure6B_plot_data_BOLD.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)

############################################################
## 13. Figure 6B
##     ALL fonts bold + enlarged axes
############################################################

group_colors <- c(
  "Normal" = "#A8D5C2",
  "Tumor" = "#E98B43"
)

pB <- ggplot2::ggplot(
  B_long,
  ggplot2::aes(
    x = CellType,
    y = Score,
    fill = Group
  )
) +
  
  ggplot2::geom_boxplot(
    width = 0.66,
    outlier.shape = NA,
    position = ggplot2::position_dodge(
      width = 0.76
    ),
    linewidth = 0.65
  ) +
  
  ggplot2::geom_point(
    ggplot2::aes(
      color = Group
    ),
    position = ggplot2::position_jitterdodge(
      jitter.width = 0.10,
      dodge.width = 0.76
    ),
    size = 1.0,
    alpha = 0.35,
    show.legend = FALSE
  ) +
  
  ggplot2::geom_text(
    data = star_y,
    ggplot2::aes(
      x = CellType,
      y = y,
      label = Signif
    ),
    inherit.aes = FALSE,
    size = 5.2,
    fontface = "bold"
  ) +
  
  ggplot2::scale_x_discrete(
    labels = cell_display_labels[
      cell_types_28
    ]
  ) +
  
  ggplot2::scale_fill_manual(
    values = group_colors
  ) +
  
  ggplot2::scale_color_manual(
    values = group_colors
  ) +
  
  ggplot2::labs(
    x = NULL,
    y = "ssGSEA score",
    fill = "Group"
  ) +
  
  ggplot2::theme_classic(
    base_size = 16
  ) +
  
  ggplot2::theme(
    
    ########################################################
    ## Make ALL text bold
    ########################################################
    text = ggplot2::element_text(
      face = "bold",
      color = "black"
    ),
    
    ########################################################
    ## Enlarged x-axis fonts
    ########################################################
    axis.text.x = ggplot2::element_text(
      angle = 58,
      hjust = 1,
      vjust = 1,
      size = 21,
      face = "bold",
      color = "black",
      margin = ggplot2::margin(
        t = 8
      )
    ),
    
    ########################################################
    ## Enlarged y-axis fonts
    ########################################################
    axis.text.y = ggplot2::element_text(
      size = 15,
      face = "bold",
      color = "black"
    ),
    
    axis.title.y = ggplot2::element_text(
      size = 20,
      face = "bold",
      color = "black",
      margin = ggplot2::margin(
        r = 10
      )
    ),
    
    axis.line = ggplot2::element_line(
      linewidth = 0.8,
      color = "black"
    ),
    
    axis.ticks = ggplot2::element_line(
      linewidth = 0.7,
      color = "black"
    ),
    
    ########################################################
    ## Legend
    ########################################################
    legend.position = "top",
    
    legend.title = ggplot2::element_text(
      size = 16,
      face = "bold",
      color = "black"
    ),
    
    legend.text = ggplot2::element_text(
      size = 15,
      face = "bold",
      color = "black"
    ),
    
    ########################################################
    ## Margin
    ########################################################
    plot.margin = ggplot2::margin(
      t = 12,
      r = 10,
      b = 34,
      l = 10
    )
  )

print(
  pB
)

ggplot2::ggsave(
  filename = file.path(
    output_dir,
    "Figure6B_28immune_paired_boxplot_BOLD_XXlargeX.pdf"
  ),
  plot = pB,
  width = 25,
  height = 11.5,
  limitsize = FALSE
)

ggplot2::ggsave(
  filename = file.path(
    output_dir,
    "Figure6B_28immune_paired_boxplot_BOLD_XXlargeX.png"
  ),
  plot = pB,
  width = 25,
  height = 11.5,
  dpi = 600,
  limitsize = FALSE
)

############################################################
## 14. Figure 6C preparation
##     Significant immune cells from Figure 6B
############################################################

sig_cells <- B_stats |>
  dplyr::filter(
    FDR < 0.05
  ) |>
  dplyr::arrange(
    FDR
  ) |>
  dplyr::pull(
    CellType
  )

cat(
  "\nImmune cells significant at BH FDR < 0.05:",
  length(sig_cells),
  "\n"
)

print(
  sig_cells
)

if (length(sig_cells) == 0) {
  
  writeLines(
    c(
      "Figure 6C was not generated.",
      "No immune-cell signature remained significant after BH correction in Figure 6B."
    ),
    file.path(
      output_dir,
      "Figure6C_NOT_GENERATED_BOLD.txt"
    )
  )
  
} else {
  
  ##########################################################
  ## 15. Read expression matrix
  ##########################################################
  
  expr_raw <- read.csv(
    expr_csv,
    row.names = 1,
    check.names = FALSE
  )
  
  expr_mat <- as.matrix(
    expr_raw
  )
  
  mode(
    expr_mat
  ) <- "numeric"
  
  expr_mat <- limma::avereps(
    expr_mat
  )
  
  missing_target_genes <- setdiff(
    target_genes,
    rownames(expr_mat)
  )
  
  if (length(missing_target_genes) > 0) {
    stop(
      "Target genes missing from expression matrix: ",
      paste(
        missing_target_genes,
        collapse = ", "
      )
    )
  }
  
  ##########################################################
  ## 16. Match expression and ssGSEA samples
  ##########################################################
  
  common_corr_samples <- Reduce(
    intersect,
    list(
      colnames(expr_mat),
      colnames(ssgsea_scores),
      rownames(meta)
    )
  )
  
  if (length(common_corr_samples) < 3) {
    stop(
      "Too few common samples for Figure 6C correlation."
    )
  }
  
  expr_target_df <- as.data.frame(
    t(
      expr_mat[
        target_genes,
        common_corr_samples,
        drop = FALSE
      ]
    ),
    check.names = FALSE
  )
  
  immune_corr_df <- as.data.frame(
    t(
      ssgsea_scores[
        sig_cells,
        common_corr_samples,
        drop = FALSE
      ]
    ),
    check.names = FALSE
  )
  
  ##########################################################
  ## 17. Spearman correlations
  ##########################################################
  
  corr_list <- list()
  k <- 1L
  
  for (g in target_genes) {
    
    for (ct in sig_cells) {
      
      x <- as.numeric(
        expr_target_df[[g]]
      )
      
      y <- as.numeric(
        immune_corr_df[[ct]]
      )
      
      keep <- complete.cases(
        x,
        y
      )
      
      if (
        sum(keep) < 3 ||
        stats::sd(
          x[keep]
        ) == 0 ||
        stats::sd(
          y[keep]
        ) == 0
      ) {
        next
      }
      
      tst <- stats::cor.test(
        x[
          keep
        ],
        y[
          keep
        ],
        method = "spearman",
        exact = FALSE
      )
      
      corr_list[[k]] <- data.frame(
        Gene = g,
        CellType = ct,
        Spearman_rho = as.numeric(
          tst$estimate
        ),
        P_value = tst$p.value,
        N = sum(
          keep
        ),
        stringsAsFactors = FALSE
      )
      
      k <- k + 1L
    }
  }
  
  corr_df <- dplyr::bind_rows(
    corr_list
  ) |>
    dplyr::mutate(
      FDR = stats::p.adjust(
        P_value,
        method = "BH"
      ),
      Signif = p_to_star(
        FDR
      ),
      Abs_rho = abs(
        Spearman_rho
      )
    )
  
  write.csv(
    corr_df,
    file.path(
      output_dir,
      "Figure6C_gene_immune_Spearman_results_BOLD.csv"
    ),
    row.names = FALSE,
    quote = FALSE
  )
  
  ##########################################################
  ## 18. Figure 6C factor order and display labels
  ##########################################################
  
  corr_df$Gene <- factor(
    corr_df$Gene,
    levels = rev(
      target_genes
    )
  )
  
  corr_df$CellType <- factor(
    corr_df$CellType,
    levels = sig_cells
  )
  
  sig_display_labels <- cell_display_labels[
    sig_cells
  ]
  
  ##########################################################
  ## 19. Figure 6C
  ##     ALL fonts bold + enlarged axes
  ##########################################################
  
  pC <- ggplot2::ggplot(
    corr_df,
    ggplot2::aes(
      x = CellType,
      y = Gene
    )
  ) +
    
    ggplot2::geom_tile(
      fill = "white",
      color = "grey78",
      linewidth = 0.8
    ) +
    
    ggplot2::geom_point(
      ggplot2::aes(
        size = Abs_rho,
        fill = Spearman_rho
      ),
      shape = 22,
      color = "grey35",
      stroke = 0.8
    ) +
    
    ggplot2::geom_text(
      ggplot2::aes(
        label = Signif
      ),
      size = 5.4,
      fontface = "bold",
      color = "black"
    ) +
    
    ggplot2::scale_x_discrete(
      labels = sig_display_labels
    ) +
    
    ggplot2::scale_size(
      range = c(
        6,
        17
      ),
      guide = "none"
    ) +
    
    ggplot2::scale_fill_gradient2(
      low = "#3B4CC0",
      mid = "white",
      high = "#E6994A",
      midpoint = 0,
      limits = c(
        -1,
        1
      ),
      name = "Spearman rho"
    ) +
    
    ######################################################
  ## Force each tile to be a true square
  ######################################################
  ggplot2::coord_fixed(
    ratio = 1
  ) +
    
    ggplot2::labs(
      x = NULL,
      y = NULL
    ) +
    
    ggplot2::theme_minimal(
      base_size = 17
    ) +
    
    ggplot2::theme(
      
      ######################################################
      ## Make ALL text bold
      ######################################################
      text = ggplot2::element_text(
        face = "bold",
        color = "black"
      ),
      
      panel.grid = ggplot2::element_blank(),
      
      ######################################################
      ## Enlarged x-axis fonts
      ######################################################
      axis.text.x = ggplot2::element_text(
        angle = 52,
        hjust = 1,
        vjust = 1,
        size = 21,
        face = "bold",
        color = "black",
        margin = ggplot2::margin(
          t = 10
        )
      ),
      
      ######################################################
      ## Enlarged y-axis fonts
      ######################################################
      axis.text.y = ggplot2::element_text(
        size = 18,
        face = "bold",
        color = "black"
      ),
      
      ######################################################
      ## Legend fonts
      ######################################################
      legend.title = ggplot2::element_text(
        size = 15,
        face = "bold",
        color = "black"
      ),
      
      legend.text = ggplot2::element_text(
        size = 13.5,
        face = "bold",
        color = "black"
      ),
      
      plot.margin = ggplot2::margin(
        t = 10,
        r = 12,
        b = 32,
        l = 10
      )
    )
  
  print(
    pC
  )
  
  ##########################################################
  ## Save Figure 6C
  ##
  ## Because coord_fixed(ratio = 1) is used, the panel itself
  ## preserves true square tiles. We also choose width/height
  ## based on the number of columns and rows so the overall
  ## saved figure remains proportionate and readable.
  ##########################################################
  
  n_cols_c <- length(
    sig_cells
  )
  
  n_rows_c <- length(
    target_genes
  )
  
  fig6C_width <- max(
    14,
    1.05 * n_cols_c + 6.5
  )
  
  fig6C_height <- max(
    7.5,
    fig6C_width * n_rows_c / n_cols_c + 3.4
  )
  
  ggplot2::ggsave(
    filename = file.path(
      output_dir,
      "Figure6C_gene_immune_correlation_BOLD_square_XXlargeX.pdf"
    ),
    plot = pC,
    width = fig6C_width,
    height = fig6C_height,
    limitsize = FALSE
  )
  
  ggplot2::ggsave(
    filename = file.path(
      output_dir,
      "Figure6C_gene_immune_correlation_BOLD_square_XXlargeX.png"
    ),
    plot = pC,
    width = fig6C_width,
    height = fig6C_height,
    dpi = 600,
    limitsize = FALSE
  )
}

############################################################
## 20. Save summary
############################################################

summary_lines <- c(
  "Figure 6B-C analysis summary",
  "============================",
  "",
  paste0(
    "Matched pairs used in Figure 6B: ",
    length(valid_pairs)
  ),
  paste0(
    "Immune-cell signatures tested in Figure 6B: ",
    length(cell_types_28)
  ),
  paste0(
    "Immune cells significant at BH FDR < 0.05: ",
    sum(
      B_stats$FDR < 0.05,
      na.rm = TRUE
    )
  ),
  "",
  "Figure 6B:",
  "- paired Wilcoxon signed-rank test",
  "- Benjamini-Hochberg correction across 28 immune-cell signatures",
  "- normalized ssGSEA scores displayed directly",
  "- all fonts bold; enlarged x/y-axis labels; Figure 6C uses square tiles",
  "",
  "Figure 6C:",
  "- Spearman correlation",
  "- Benjamini-Hochberg correction across displayed gene-cell tests",
  "- all fonts bold; enlarged x/y-axis labels; Figure 6C uses square tiles"
)

writeLines(
  summary_lines,
  file.path(
    output_dir,
    "Figure6BC_analysis_summary_BOLD.txt"
  )
)

cat("\n============================================================\n")
cat("Figure 6B-C analysis completed.\n")
cat("Output directory:\n")
cat(output_dir, "\n")
cat("============================================================\n")
