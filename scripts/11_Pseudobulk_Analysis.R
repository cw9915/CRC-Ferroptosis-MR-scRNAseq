############################################################
## GSE294300
## Three genes × all 10 major cell types
## Donor-level paired pseudobulk expression plots
##
## Genes:
## AKR1C2 / MAPK3 / MAP3K11
##
## IMPORTANT:
## This script DOES NOT rerun edgeR.
## It directly uses the finalized pseudobulk logCPM and
## edgeR statistics generated previously.
############################################################

rm(list = ls())
options(stringsAsFactors = FALSE)
Sys.setenv(LANGUAGE = "en")

############################################################
## 1. Load packages
############################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(grid)

############################################################
## 2. Working directory
############################################################

work_dir <- paste0(
  "C:/Users/92066/Desktop/单细胞孟德尔/Figure7/",
  "GSE294300/007pseudobulk_3gene"
)

if (!dir.exists(work_dir)) {
  stop("目录不存在：", work_dir)
}

setwd(work_dir)

############################################################
## 3. Input files
############################################################

expression_file <- file.path(
  work_dir,
  "3genes_pseudobulk_logCPM.csv"
)

statistics_file <- file.path(
  work_dir,
  "3genes_paired_pseudobulk_edgeR_statistics.csv"
)

if (!file.exists(expression_file)) {
  stop(
    "没有找到：",
    expression_file
  )
}

if (!file.exists(statistics_file)) {
  stop(
    "没有找到：",
    statistics_file
  )
}

############################################################
## 4. Output directory
############################################################

output_dir <- file.path(
  work_dir,
  "full_celltype_paired_plots"
)

if (!dir.exists(output_dir)) {
  dir.create(
    output_dir,
    recursive = TRUE
  )
}

############################################################
## 5. Read finalized results
############################################################

pb_expression <- read.csv(
  expression_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

pb_stats <- read.csv(
  statistics_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

############################################################
## 6. Check columns
############################################################

required_expr_cols <- c(
  "patient",
  "group",
  "celltype",
  "gene",
  "logCPM"
)

missing_expr <- setdiff(
  required_expr_cols,
  colnames(pb_expression)
)

if (length(missing_expr) > 0) {
  stop(
    "logCPM文件缺少以下列：",
    paste(
      missing_expr,
      collapse = ", "
    )
  )
}

required_stat_cols <- c(
  "celltype",
  "gene",
  "n_pairs",
  "logFC",
  "PValue",
  "FDR"
)

missing_stat <- setdiff(
  required_stat_cols,
  colnames(pb_stats)
)

if (length(missing_stat) > 0) {
  stop(
    "statistics文件缺少以下列：",
    paste(
      missing_stat,
      collapse = ", "
    )
  )
}

############################################################
## 7. Define genes and cell types
############################################################

target_genes <- c(
  "AKR1C2",
  "MAPK3",
  "MAP3K11"
)

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

############################################################
## 8. Keep required data
############################################################

pb_expression <- pb_expression %>%
  filter(
    gene %in% target_genes,
    celltype %in% celltype_order,
    group %in% c(
      "Normal",
      "Tumor"
    )
  )

pb_stats <- pb_stats %>%
  filter(
    gene %in% target_genes,
    celltype %in% celltype_order
  )

############################################################
## 9. Factor order
############################################################

pb_expression$gene <- factor(
  pb_expression$gene,
  levels = target_genes
)

pb_expression$celltype <- factor(
  pb_expression$celltype,
  levels = celltype_order
)

pb_expression$group <- factor(
  pb_expression$group,
  levels = c(
    "Normal",
    "Tumor"
  )
)

pb_stats$gene <- factor(
  pb_stats$gene,
  levels = target_genes
)

pb_stats$celltype <- factor(
  pb_stats$celltype,
  levels = celltype_order
)

############################################################
## 10. Check data
############################################################

cat("\n====================================================\n")
cat("Expression data\n")
cat("====================================================\n")

print(
  table(
    pb_expression$gene,
    pb_expression$celltype
  )
)

cat("\n====================================================\n")
cat("Statistics\n")
cat("====================================================\n")

print(
  pb_stats %>%
    select(
      gene,
      celltype,
      n_pairs,
      logFC,
      PValue,
      FDR
    )
)

############################################################
## 11. FDR label function
############################################################

format_fdr <- function(x) {
  
  if (is.na(x)) {
    return("FDR = NA")
  }
  
  if (x < 0.001) {
    return("FDR < 0.001")
  }
  
  paste0(
    "FDR = ",
    formatC(
      x,
      format = "fg",
      digits = 3
    )
  )
}

############################################################
## 12. Add annotation labels
############################################################

pb_stats <- pb_stats %>%
  rowwise() %>%
  mutate(
    FDR_label = format_fdr(FDR),
    
    pair_label = ifelse(
      is.na(n_pairs),
      "n = NA",
      paste0(
        "n = ",
        n_pairs,
        " paired donors"
      )
    ),
    
    significance = case_when(
      is.na(FDR) ~ "NA",
      FDR < 0.001 ~ "***",
      FDR < 0.01 ~ "**",
      FDR < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  ) %>%
  ungroup()

############################################################
## 13. Colors
############################################################

group_colors <- c(
  "Normal" = "#3C8DBC",
  "Tumor" = "#D95F59"
)

############################################################
## 14. Function for one gene
############################################################

make_gene_full_plot <- function(gene_name) {
  
  ##########################################################
  ## Expression data for one gene
  ##########################################################
  
  dat <- pb_expression %>%
    filter(
      as.character(gene) == gene_name
    )
  
  dat$celltype <- factor(
    dat$celltype,
    levels = celltype_order
  )
  
  dat$group <- factor(
    dat$group,
    levels = c(
      "Normal",
      "Tumor"
    )
  )
  
  ##########################################################
  ## Statistics for one gene
  ##########################################################
  
  stat_gene <- pb_stats %>%
    filter(
      as.character(gene) == gene_name
    )
  
  stat_gene$celltype <- factor(
    stat_gene$celltype,
    levels = celltype_order
  )
  
  ##########################################################
  ## Get plotting range for every cell type
  ##########################################################
  
  range_df <- dat %>%
    group_by(celltype) %>%
    summarise(
      ymin = min(
        logCPM,
        na.rm = TRUE
      ),
      
      ymax = max(
        logCPM,
        na.rm = TRUE
      ),
      
      .groups = "drop"
    )
  
  ##########################################################
  ## Ensure all 10 cell types are represented
  ##########################################################
  
  base_celltypes <- data.frame(
    celltype = factor(
      celltype_order,
      levels = celltype_order
    )
  )
  
  annotation_df <- base_celltypes %>%
    left_join(
      range_df,
      by = "celltype"
    ) %>%
    left_join(
      stat_gene %>%
        select(
          celltype,
          FDR,
          FDR_label,
          n_pairs,
          pair_label,
          significance
        ),
      by = "celltype"
    )
  
  ##########################################################
  ## Annotation positions
  ##########################################################
  
  annotation_df <- annotation_df %>%
    mutate(
      plot_range = ymax - ymin,
      
      plot_range = ifelse(
        is.na(plot_range) |
          !is.finite(plot_range) |
          plot_range <= 0,
        1,
        plot_range
      ),
      
      y_FDR =
        ymax +
        plot_range * 0.22,
      
      y_pair =
        ymax +
        plot_range * 0.10
    )
  
  ##########################################################
  ## For completely missing panels
  ##########################################################
  
  annotation_df$y_FDR[
    !is.finite(annotation_df$y_FDR)
  ] <- 1
  
  annotation_df$y_pair[
    !is.finite(annotation_df$y_pair)
  ] <- 0.7
  
  ##########################################################
  ## Plot
  ##########################################################
  
  p <- ggplot(
    dat,
    aes(
      x = group,
      y = logCPM
    )
  ) +
    
    ########################################################
  ## Pairing lines
  ########################################################
  
  geom_line(
    aes(
      group = patient
    ),
    color = "grey65",
    linewidth = 0.65,
    alpha = 0.82
  ) +
    
    ########################################################
  ## Donor-level points
  ########################################################
  
  geom_point(
    aes(
      color = group
    ),
    size = 3.2,
    alpha = 0.95
  ) +
    
    ########################################################
  ## Median crossbar
  ########################################################
  
  stat_summary(
    aes(
      color = group
    ),
    fun = median,
    fun.min = median,
    fun.max = median,
    geom = "crossbar",
    width = 0.30,
    linewidth = 1.05
  ) +
    
    ########################################################
  ## FDR
  ########################################################
  
  geom_text(
    data = annotation_df,
    aes(
      x = 1.5,
      y = y_FDR,
      label = FDR_label
    ),
    inherit.aes = FALSE,
    size = 4.6,
    fontface = "bold",
    color = "black"
  ) +
    
    ########################################################
  ## Paired donor number
  ########################################################
  
  geom_text(
    data = annotation_df,
    aes(
      x = 1.5,
      y = y_pair,
      label = pair_label
    ),
    inherit.aes = FALSE,
    size = 4.1,
    color = "black"
  ) +
    
    ########################################################
  ## 10 cell types
  ########################################################
  
  facet_wrap(
    ~ celltype,
    scales = "free_y",
    ncol = 5,
    drop = FALSE
  ) +
    
    ########################################################
  ## Colors
  ########################################################
  
  scale_color_manual(
    values = group_colors,
    drop = FALSE
  ) +
    
    ########################################################
  ## Y-axis expansion
  ########################################################
  
  scale_y_continuous(
    expand = expansion(
      mult = c(
        0.06,
        0.32
      )
    )
  ) +
    
    ########################################################
  ## Labels
  ########################################################
  
  labs(
    title = gene_name,
    x = NULL,
    y = "Pseudobulk logCPM"
  ) +
    
    ########################################################
  ## Theme
  ########################################################
  
  theme_classic(
    base_size = 15
  ) +
    
    theme(
      
      ######################################################
      ## Gene title
      ######################################################
      
      plot.title = element_text(
        size = 25,
        face = "bold",
        hjust = 0.5,
        color = "black",
        margin = margin(
          b = 12
        )
      ),
      
      ######################################################
      ## Cell-type titles
      ######################################################
      
      strip.background = element_blank(),
      
      strip.text = element_text(
        size = 15,
        face = "bold",
        color = "black",
        margin = margin(
          b = 5
        )
      ),
      
      ######################################################
      ## X-axis
      ######################################################
      
      axis.text.x = element_text(
        size = 14,
        face = "bold",
        color = "black"
      ),
      
      ######################################################
      ## Y-axis
      ######################################################
      
      axis.text.y = element_text(
        size = 12,
        color = "black"
      ),
      
      axis.title.y = element_text(
        size = 17,
        face = "bold",
        color = "black",
        margin = margin(
          r = 8
        )
      ),
      
      ######################################################
      ## Lines
      ######################################################
      
      axis.line = element_line(
        linewidth = 0.8,
        color = "black"
      ),
      
      axis.ticks = element_line(
        linewidth = 0.7,
        color = "black"
      ),
      
      ######################################################
      ## Other settings
      ######################################################
      
      legend.position = "none",
      
      panel.grid = element_blank(),
      
      panel.spacing = grid::unit(
        1.25,
        "lines"
      ),
      
      plot.margin = margin(
        15,
        15,
        15,
        15
      )
    ) +
    
    coord_cartesian(
      clip = "off"
    )
  
  return(p)
}

############################################################
## 15. AKR1C2
############################################################

p_AKR1C2 <- make_gene_full_plot(
  "AKR1C2"
)

print(p_AKR1C2)

ggsave(
  filename = file.path(
    output_dir,
    "AKR1C2_all_10_celltypes_paired.pdf"
  ),
  plot = p_AKR1C2,
  width = 16,
  height = 8.5,
  device = cairo_pdf,
  bg = "white"
)

ggsave(
  filename = file.path(
    output_dir,
    "AKR1C2_all_10_celltypes_paired.png"
  ),
  plot = p_AKR1C2,
  width = 16,
  height = 8.5,
  dpi = 600,
  bg = "white"
)

############################################################
## 16. MAPK3
############################################################

p_MAPK3 <- make_gene_full_plot(
  "MAPK3"
)

print(p_MAPK3)

ggsave(
  filename = file.path(
    output_dir,
    "MAPK3_all_10_celltypes_paired.pdf"
  ),
  plot = p_MAPK3,
  width = 16,
  height = 8.5,
  device = cairo_pdf,
  bg = "white"
)

ggsave(
  filename = file.path(
    output_dir,
    "MAPK3_all_10_celltypes_paired.png"
  ),
  plot = p_MAPK3,
  width = 16,
  height = 8.5,
  dpi = 600,
  bg = "white"
)

############################################################
## 17. MAP3K11
############################################################

p_MAP3K11 <- make_gene_full_plot(
  "MAP3K11"
)

print(p_MAP3K11)

ggsave(
  filename = file.path(
    output_dir,
    "MAP3K11_all_10_celltypes_paired.pdf"
  ),
  plot = p_MAP3K11,
  width = 16,
  height = 8.5,
  device = cairo_pdf,
  bg = "white"
)

ggsave(
  filename = file.path(
    output_dir,
    "MAP3K11_all_10_celltypes_paired.png"
  ),
  plot = p_MAP3K11,
  width = 16,
  height = 8.5,
  dpi = 600,
  bg = "white"
)

############################################################
## 18. Combined supplementary figure
##
## A = AKR1C2
## B = MAPK3
## C = MAP3K11
############################################################

p_combined <- (
  p_AKR1C2 +
    plot_annotation(
      tag_levels = "A"
    )
) /
  (
    p_MAPK3 +
      plot_annotation(
        tag_levels = "B"
      )
  ) /
  (
    p_MAP3K11 +
      plot_annotation(
        tag_levels = "C"
      )
  )

############################################################
## 19. Save combined PDF
############################################################

ggsave(
  filename = file.path(
    output_dir,
    "Three_genes_all_10_celltypes_paired_COMBINED.pdf"
  ),
  plot = p_combined,
  width = 16,
  height = 25,
  device = cairo_pdf,
  limitsize = FALSE,
  bg = "white"
)

############################################################
## 20. Save combined PNG
############################################################

ggsave(
  filename = file.path(
    output_dir,
    "Three_genes_all_10_celltypes_paired_COMBINED.png"
  ),
  plot = p_combined,
  width = 16,
  height = 25,
  dpi = 600,
  limitsize = FALSE,
  bg = "white"
)

############################################################
## 21. Save plotting statistics
############################################################

write.csv(
  pb_stats,
  file.path(
    output_dir,
    "Three_genes_all_celltypes_statistics_used_for_plot.csv"
  ),
  row.names = FALSE
)

############################################################
## 22. Print important results
############################################################

cat("\n====================================================\n")
cat("FDR < 0.05 results\n")
cat("====================================================\n")

print(
  pb_stats %>%
    filter(
      !is.na(FDR),
      FDR < 0.05
    ) %>%
    select(
      gene,
      celltype,
      n_pairs,
      logFC,
      PValue,
      FDR
    ),
  row.names = FALSE
)

############################################################
## 23. Check pair numbers
############################################################

cat("\n====================================================\n")
cat("Pair numbers for all gene × cell-type combinations\n")
cat("====================================================\n")

print(
  pb_stats %>%
    select(
      gene,
      celltype,
      n_pairs
    ) %>%
    arrange(
      gene,
      celltype
    ),
  row.names = FALSE
)

############################################################
## 24. Save session information
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
## 25. Finished
############################################################

cat("\n############################################################\n")
cat("## THREE-GENE FULL CELL-TYPE PLOTS COMPLETED\n")
cat("############################################################\n")

cat("\nInput directory:\n")
cat(work_dir, "\n")

cat("\nOutput directory:\n")
cat(output_dir, "\n")

cat("\nMain files:\n")
cat("AKR1C2_all_10_celltypes_paired.pdf/png\n")
cat("MAPK3_all_10_celltypes_paired.pdf/png\n")
cat("MAP3K11_all_10_celltypes_paired.pdf/png\n")
cat("Three_genes_all_10_celltypes_paired_COMBINED.pdf/png\n")