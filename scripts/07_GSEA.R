############################################################
## Figure 5：3个候选基因 Hallmark GSEA 分析
## 严格版本：先运行全部 Hallmark，再筛选三基因共同显著通路
##
## 输入表达矩阵：GSE44861_gene_symbol_expression.csv
## Hallmark文件：h.all.v2025.1.Hs.symbols.gmt
## 输出目录：Figure5
## 候选基因：AKR1C2, MAPK3, MAP3K11
##
## 核心逻辑：
## 1. 对 AKR1C2、MAPK3、MAP3K11 分别运行全部 Hallmark GSEA
## 2. 每个基因严格以 BH-adjusted P < 0.05 定义显著通路
## 3. 对三个基因的显著通路取严格交集（Reduce(intersect, ...)）
## 4. 若共同显著通路 > 5 条，则按三基因 Mean |NES| 从高到低选前5条用于 Figure 5
## 5. 若共同显著通路 <= 5 条，则全部展示，不人为补足、不降低显著性阈值
## 6. 手动绘制 running enrichment score，使上方曲线与下方 hit ticks 严格对齐
############################################################

Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list = ls())
set.seed(12345)

############################################################
## 1. 加载 R 包
############################################################

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(ggplot2)
  library(dplyr)
  library(grid)
})

if (!requireNamespace("patchwork", quietly = TRUE)) {
  install.packages("patchwork")
}
library(patchwork)

############################################################
## 2. 设置路径
############################################################

expr_file <- "C:/Users/92066/Desktop/单细胞孟德尔/Figure1/GSE44861_gene_symbol_expression.csv"

gmt_file <- "C:/Users/92066/Desktop/单细胞孟德尔/Figure5/h.all.v2025.1.Hs.symbols.gmt"

outdir <- "C:/Users/92066/Desktop/单细胞孟德尔/Figure5"

if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

setwd(outdir)

############################################################
## 3. 读取 Hallmark GMT 文件
############################################################

if (!file.exists(gmt_file)) {
  stop("没有找到 Hallmark GMT 文件，请检查路径：", gmt_file)
}

term2gene <- clusterProfiler::read.gmt(gmt_file)

colnames(term2gene) <- c("gs_name", "gene_symbol")

term2gene <- term2gene[
  !is.na(term2gene$gs_name) &
    !is.na(term2gene$gene_symbol),
]

term2gene <- term2gene[
  term2gene$gs_name != "" &
    term2gene$gene_symbol != "",
]

term2gene <- unique(term2gene)

term2name <- data.frame(
  gs_name = unique(term2gene$gs_name),
  name = unique(term2gene$gs_name),
  stringsAsFactors = FALSE
)

hallmark_gene_sets <- split(
  term2gene$gene_symbol,
  term2gene$gs_name
)

cat("Hallmark通路数量：", length(unique(term2gene$gs_name)), "\n")
cat("Hallmark基因-通路关系数量：", nrow(term2gene), "\n")

if (length(unique(term2gene$gs_name)) < 40) {
  stop("Hallmark通路数量异常，请确认 GMT 文件是否为 h.all.v2025.1.Hs.symbols.gmt")
}

############################################################
## 4. 读取表达矩阵
############################################################

if (!file.exists(expr_file)) {
  stop("没有找到表达矩阵文件，请检查路径：", expr_file)
}

expr_raw <- read.csv(
  expr_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("表达矩阵原始维度：", nrow(expr_raw), "rows ×", ncol(expr_raw), "columns\n")
cat("表达矩阵前几列名：\n")
print(head(colnames(expr_raw)))

############################################################
## 5. 整理表达矩阵
############################################################

gene_col <- colnames(expr_raw)[1]

cat("默认使用第一列作为基因名：", gene_col, "\n")

expr_raw[[gene_col]] <- as.character(expr_raw[[gene_col]])

# 去除基因名缺失或为空的行
expr_raw <- expr_raw[
  !is.na(expr_raw[[gene_col]]) &
    trimws(expr_raw[[gene_col]]) != "",
  ,
  drop = FALSE
]

# 保存基因名；表达列全部转换为 numeric
Gene <- expr_raw[[gene_col]]
expr_data <- expr_raw[, -1, drop = FALSE]
expr_data <- as.data.frame(
  lapply(expr_data, function(x) suppressWarnings(as.numeric(as.character(x)))),
  check.names = FALSE
)

# 去除所有样本均为 NA 的基因
keep_nonempty <- rowSums(!is.na(expr_data)) > 0
expr_data <- expr_data[keep_nonempty, , drop = FALSE]
Gene <- Gene[keep_nonempty]

# 先把 Gene 作为普通列处理，避免重复 gene symbol 直接设置 rownames 时报错
expr_data$Gene <- Gene
sample_cols <- setdiff(colnames(expr_data), "Gene")

# 重复 gene symbol：保留跨样本平均表达量最高的一行
expr_data$mean_exp <- rowMeans(
  expr_data[, sample_cols, drop = FALSE],
  na.rm = TRUE
)

expr_data <- expr_data %>%
  dplyr::arrange(desc(mean_exp)) %>%
  dplyr::distinct(Gene, .keep_all = TRUE)

rownames(expr_data) <- expr_data$Gene
expr_data$Gene <- NULL
expr_data$mean_exp <- NULL

expr_mat <- as.matrix(expr_data)
mode(expr_mat) <- "numeric"

if (ncol(expr_mat) < 3) {
  stop("表达矩阵样本数过少，无法可靠计算 Pearson 相关性。")
}

cat("清理后表达矩阵维度：", nrow(expr_mat), "genes ×", ncol(expr_mat), "samples\n")

############################################################
## 6. 设置候选基因
############################################################

biomarkers <- c("AKR1C2", "MAPK3", "MAP3K11")

missing_genes <- setdiff(biomarkers, rownames(expr_mat))

if (length(missing_genes) > 0) {
  stop(
    "以下候选基因在表达矩阵中没有找到：",
    paste(missing_genes, collapse = ", "),
    "\n请检查表达矩阵第一列是否为标准基因Symbol。"
  )
}

cat("候选基因均存在于表达矩阵中：\n")
print(biomarkers)

############################################################
## 7. 颜色设置
############################################################

plot_color_pool <- c(
  "#C00000",  # red
  "#2F62B3",  # blue
  "#57A639",  # green
  "#7B3FB5",  # purple
  "#D97A1B",  # orange
  "#009E73",
  "#CC79A7",
  "#0072B2"
)

############################################################
## 8. 手动计算 running enrichment score
############################################################

calc_running_score <- function(geneList, geneSet, pathway_id, exponent = 1) {
  
  geneList <- sort(geneList, decreasing = TRUE)
  all_genes <- names(geneList)
  N <- length(geneList)
  
  geneSet <- intersect(geneSet, all_genes)
  Nh <- length(geneSet)
  Nm <- N - Nh
  
  if (Nh == 0) {
    stop("通路 ", pathway_id, " 与 geneList 没有交集。")
  }
  
  if (Nm == 0) {
    stop("通路 ", pathway_id, " 覆盖了全部 geneList，无法计算 running score。")
  }
  
  hit <- all_genes %in% geneSet
  
  hit_weight <- abs(geneList[hit])^exponent
  hit_weight_sum <- sum(hit_weight, na.rm = TRUE)
  
  if (!is.finite(hit_weight_sum) || hit_weight_sum == 0) {
    hit_weight <- rep(1, sum(hit))
    hit_weight_sum <- sum(hit_weight)
  }
  
  running_increment <- rep(NA_real_, N)
  running_increment[hit] <- abs(geneList[hit])^exponent / hit_weight_sum
  running_increment[!hit] <- -1 / Nm
  
  runningScore <- cumsum(running_increment)
  
  df <- data.frame(
    x = seq_len(N),
    runningScore = as.numeric(runningScore),
    position = as.integer(hit),
    Description = pathway_id,
    stringsAsFactors = FALSE
  )
  
  return(df)
}

############################################################
## 9. 自定义 GSEA 图：上下严格对齐
############################################################

plot_gsea_aligned <- function(
    geneList,
    gene_sets,
    plot_ids,
    plot_colors,
    gene,
    png_file,
    pdf_file,
    width = 10.5,
    height = 6.2
) {
  
  if (length(geneList) == 0) {
    stop("geneList 为空，无法作图。")
  }
  
  missing_plot_ids <- setdiff(plot_ids, names(gene_sets))
  
  if (length(missing_plot_ids) > 0) {
    stop(
      "gene_sets 中缺少以下通路：",
      paste(missing_plot_ids, collapse = ", ")
    )
  }
  
  gs_list <- lapply(
    plot_ids,
    function(pid) {
      calc_running_score(
        geneList = geneList,
        geneSet = gene_sets[[pid]],
        pathway_id = pid,
        exponent = 1
      )
    }
  )
  
  gs_df <- do.call(rbind, gs_list)
  
  pathway_levels <- plot_ids
  
  gs_df$Description <- factor(
    gs_df$Description,
    levels = pathway_levels
  )
  
  x_min <- 1
  x_max <- length(geneList)
  
  plot_colors_use <- plot_colors[pathway_levels]
  
  if (any(is.na(plot_colors_use))) {
    stop("plot_colors 中存在 NA，请检查颜色名称是否与通路ID一致。")
  }
  
  ##########################################################
  ## 上方 Running Enrichment Score 曲线
  ##########################################################
  
  p_es <- ggplot(
    gs_df,
    aes(
      x = x,
      y = runningScore,
      color = Description,
      group = Description
    )
  ) +
    geom_line(
      linewidth = 0.9,
      lineend = "round"
    ) +
    scale_color_manual(
      values = plot_colors_use,
      name = NULL
    ) +
    scale_x_continuous(
      limits = c(x_min, x_max),
      expand = c(0, 0)
    ) +
    coord_cartesian(
      xlim = c(x_min, x_max),
      clip = "on"
    ) +
    labs(
      title = paste0(gene, "_hallmark"),
      x = NULL,
      y = "Running Enrichment Score"
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(
        hjust = 0,
        size = 16,
        face = "plain"
      ),
      axis.title.y = element_text(
        size = 14,
        face = "plain"
      ),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      
      legend.position = c(0.74, 1.1),
      legend.justification = c(0.5, 1),
      legend.background = element_rect(
        fill = scales::alpha("white", 0.75),
        color = NA
      ),
      legend.key = element_blank(),
      legend.text = element_text(size = 10.5),
      legend.spacing.y = unit(0.10, "cm"),
      
      panel.grid.major.x = element_line(
        color = "grey90",
        linewidth = 0.4
      ),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.border = element_blank(),
      axis.line.x = element_line(
        color = "black",
        linewidth = 0.4
      ),
      axis.line.y = element_line(
        color = "black",
        linewidth = 0.4
      ),
      plot.margin = margin(5, 10, 0, 5)
    )
  
  ##########################################################
  ## 下方 hit ticks
  ##########################################################
  
  hit_df <- gs_df %>%
    dplyr::filter(position == 1)
  
  if (nrow(hit_df) == 0) {
    stop("没有提取到 hit ticks 数据。")
  }
  
  hit_df$Description <- factor(
    as.character(hit_df$Description),
    levels = rev(pathway_levels)
  )
  
  hit_df$y_id <- as.numeric(hit_df$Description)
  
  p_hits <- ggplot(
    hit_df,
    aes(
      x = x,
      color = Description
    )
  ) +
    geom_segment(
      aes(
        xend = x,
        y = y_id - 0.35,
        yend = y_id + 0.35
      ),
      linewidth = 0.35,
      lineend = "butt"
    ) +
    scale_color_manual(
      values = plot_colors_use,
      guide = "none"
    ) +
    scale_x_continuous(
      limits = c(x_min, x_max),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      limits = c(0.5, length(pathway_levels) + 0.5),
      expand = c(0, 0)
    ) +
    coord_cartesian(
      xlim = c(x_min, x_max),
      clip = "on"
    ) +
    labs(
      x = NULL,
      y = NULL
    ) +
    theme_bw(base_size = 13) +
    theme(
      axis.title = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line.x = element_line(
        color = "black",
        linewidth = 0.4
      ),
      axis.text.x = element_text(size = 11),
      axis.ticks.x = element_line(
        color = "black",
        linewidth = 0.4
      ),
      plot.margin = margin(0, 10, 5, 5)
    )
  
  p_final <- p_es / p_hits +
    patchwork::plot_layout(
      heights = c(4.2, 0.65)
    )
  
  ggsave(
    filename = png_file,
    plot = p_final,
    width = width,
    height = height,
    dpi = 300
  )
  
  ggsave(
    filename = pdf_file,
    plot = p_final,
    width = width,
    height = height
  )
  
  return(p_final)
}

############################################################
## 10. 第一轮：对3个基因全部运行 Hallmark GSEA
############################################################

gsea_result_list <- list()
gsea_table_list <- list()
rank_score_list <- list()

for (gene in biomarkers) {
  
  cat("\n==============================\n")
  cat("第一轮GSEA分析基因：", gene, "\n")
  cat("==============================\n")
  
  tryCatch({
    
    target_exp <- as.numeric(expr_mat[gene, ])
    
    rank_score <- apply(
      expr_mat,
      1,
      function(x) {
        suppressWarnings(
          cor(
            as.numeric(x),
            target_exp,
            method = "pearson",
            use = "pairwise.complete.obs"
          )
        )
      }
    )
    
    rank_score <- rank_score[!is.na(rank_score)]
    rank_score <- rank_score[names(rank_score) != gene]
    rank_score <- rank_score[!duplicated(names(rank_score))]
    rank_score <- sort(rank_score, decreasing = TRUE)
    
    overlap_genes <- intersect(names(rank_score), term2gene$gene_symbol)
    cat("与Hallmark基因集交集数量：", length(overlap_genes), "\n")
    
    if (length(overlap_genes) < 100) {
      stop("表达矩阵基因名与 Hallmark 基因集交集过少，请检查基因名格式。")
    }
    
    ranking_df <- data.frame(
      Gene = names(rank_score),
      Correlation = as.numeric(rank_score),
      stringsAsFactors = FALSE
    )
    
    write.csv(
      ranking_df,
      file = file.path(outdir, paste0("Ranking_", gene, ".csv")),
      row.names = FALSE
    )
    
    gsea_res <- clusterProfiler::GSEA(
      geneList = rank_score,
      TERM2GENE = term2gene,
      TERM2NAME = term2name,
      minGSSize = 15,
      maxGSSize = 500,
      pvalueCutoff = 1,          # 先保留完整 Hallmark 结果，显著性在后续统一筛选
      pAdjustMethod = "BH",
      eps = 0,
      seed = TRUE,
      verbose = FALSE
    )
    
    gsea_table <- as.data.frame(gsea_res@result)
    gsea_table$Biomarker <- gene
    
    cat("GSEA结果通路数量：", nrow(gsea_table), "\n")
    
    if (nrow(gsea_table) == 0) {
      stop("GSEA 没有得到任何结果。")
    }
    
    write.csv(
      gsea_table,
      file = file.path(outdir, paste0("GSEA_", gene, "_all_results.csv")),
      row.names = FALSE
    )
    
    gsea_result_list[[gene]] <- gsea_res
    gsea_table_list[[gene]] <- gsea_table
    rank_score_list[[gene]] <- rank_score
    
    cat("完成第一轮GSEA：", gene, "\n")
    
  }, error = function(e) {
    
    cat("基因 ", gene, " 第一轮GSEA失败。\n")
    cat("错误信息：", conditionMessage(e), "\n")
    
  })
}

if (length(gsea_table_list) == 0) {
  stop("没有任何基因成功完成 GSEA，请检查前面的错误信息。")
}

# 三个候选基因必须全部成功，否则不能进行“三基因共同通路”分析
failed_genes <- setdiff(biomarkers, names(gsea_table_list))
if (length(failed_genes) > 0) {
  stop(
    "以下候选基因 GSEA 未成功完成，不能计算三基因共同通路：",
    paste(failed_genes, collapse = ", ")
  )
}

############################################################
## 11. 严格筛选三基因共同显著 Hallmark 通路
############################################################

all_gsea_df <- do.call(rbind, gsea_table_list)
rownames(all_gsea_df) <- NULL

write.csv(
  all_gsea_df,
  file = file.path(outdir, "GSEA_3_biomarkers_all_hallmark_results.csv"),
  row.names = FALSE
)

# 与 Methods 保持一致：BH-adjusted P < 0.05 定义显著富集
padj_cutoff <- 0.05

sig_table_list <- lapply(
  biomarkers,
  function(gene) {
    gsea_table_list[[gene]] %>%
      dplyr::filter(!is.na(p.adjust), p.adjust < padj_cutoff)
  }
)
names(sig_table_list) <- biomarkers

# 输出每个基因的显著通路结果及数量
sig_count_df <- data.frame(
  Biomarker = biomarkers,
  Significant_Hallmark_count = vapply(sig_table_list, nrow, integer(1)),
  stringsAsFactors = FALSE
)

write.csv(
  sig_count_df,
  file = file.path(outdir, "GSEA_3_biomarkers_significant_pathway_counts.csv"),
  row.names = FALSE
)

cat("\n各基因 BH-adjusted P < 0.05 的 Hallmark 通路数量：\n")
print(sig_count_df)

for (gene in biomarkers) {
  write.csv(
    sig_table_list[[gene]],
    file = file.path(
      outdir,
      paste0("GSEA_", gene, "_significant_padj_lt_0.05.csv")
    ),
    row.names = FALSE
  )
}

# 严格三基因交集：每条通路必须在三个基因中均满足 p.adjust < 0.05
sig_id_list <- lapply(sig_table_list, function(x) unique(x$ID))
common_sig_ids <- Reduce(intersect, sig_id_list)

cat("\n三基因共同显著 Hallmark 通路数量：", length(common_sig_ids), "\n")

if (length(common_sig_ids) == 0) {
  stop(
    "AKR1C2、MAPK3、MAP3K11 在 BH-adjusted P < 0.05 标准下没有共同显著 Hallmark 通路。",
    "请如实报告该结果，不应降低阈值或人为补足通路。"
  )
}

# 提取共同显著通路的三基因完整结果
common_long_df <- all_gsea_df %>%
  dplyr::filter(ID %in% common_sig_ids) %>%
  dplyr::arrange(ID, match(Biomarker, biomarkers))

write.csv(
  common_long_df,
  file = file.path(outdir, "GSEA_3_biomarkers_common_significant_pathways_long.csv"),
  row.names = FALSE
)

# 可复现的代表性通路排序：
# 主排序 = 三基因 Mean |NES|；次排序 = Max adjusted P、Mean adjusted P
pathway_summary <- common_long_df %>%
  dplyr::group_by(ID, Description) %>%
  dplyr::summarise(
    Gene_count = dplyr::n_distinct(Biomarker),
    Mean_abs_NES = mean(abs(NES), na.rm = TRUE),
    Mean_NES = mean(NES, na.rm = TRUE),
    Min_NES = min(NES, na.rm = TRUE),
    Max_NES = max(NES, na.rm = TRUE),
    Mean_padjust = mean(p.adjust, na.rm = TRUE),
    Max_padjust = max(p.adjust, na.rm = TRUE),
    Direction = dplyr::case_when(
      all(NES > 0, na.rm = TRUE) ~ "Consistently positive",
      all(NES < 0, na.rm = TRUE) ~ "Consistently negative",
      TRUE ~ "Mixed"
    ),
    Biomarkers = paste(biomarkers, collapse = ";"),
    .groups = "drop"
  ) %>%
  dplyr::filter(Gene_count == length(biomarkers)) %>%
  dplyr::arrange(
    desc(Mean_abs_NES),
    Max_padjust,
    Mean_padjust,
    ID
  )

if (nrow(pathway_summary) != length(common_sig_ids)) {
  stop("共同通路汇总数量与严格交集数量不一致，请检查 GSEA 结果。")
}

write.csv(
  pathway_summary,
  file = file.path(outdir, "GSEA_3_biomarkers_common_significant_pathways_ranked.csv"),
  row.names = FALSE
)

cat("\n三基因共同显著 Hallmark 通路（按 Mean |NES| 排序）：\n")
print(pathway_summary[, c(
  "ID", "Mean_abs_NES", "Max_padjust", "Direction"
)])

# Figure 5 最多展示5条；若共同显著通路少于5条，则全部展示，不人为补足
n_show <- min(5L, nrow(pathway_summary))
selected_common_ids <- pathway_summary$ID[seq_len(n_show)]

if (n_show < 5L) {
  warning(
    "三基因共同显著 Hallmark 通路仅有 ", n_show,
    " 条。Figure 5 将全部展示，不降低显著性阈值，也不补足至5条。"
  )
}

selected_pathway_summary <- pathway_summary %>%
  dplyr::filter(ID %in% selected_common_ids) %>%
  dplyr::mutate(Figure5_rank = match(ID, selected_common_ids)) %>%
  dplyr::arrange(Figure5_rank)

cat("\nFigure 5 最终展示通路（严格来自三基因共同显著交集）：\n")
print(selected_pathway_summary[, c(
  "Figure5_rank", "ID", "Mean_abs_NES", "Max_padjust", "Direction"
)])

write.csv(
  selected_pathway_summary,
  file = file.path(outdir, "GSEA_3_biomarkers_Figure5_selected_pathways.csv"),
  row.names = FALSE
)

############################################################
## 12. 第二轮：用严格共同显著通路为3个基因统一绘图
############################################################

plot_colors <- setNames(
  plot_color_pool[seq_len(length(selected_common_ids))],
  selected_common_ids
)

summary_list <- list()

for (gene in biomarkers) {
  
  cat("\n==============================\n")
  cat("第二轮绘图基因：", gene, "\n")
  cat("==============================\n")
  
  tryCatch({
    
    gsea_table <- gsea_table_list[[gene]]
    rank_score <- rank_score_list[[gene]]
    
    # selected_common_ids 来自严格三基因交集，因此三个基因使用完全相同的 plot_ids
    plot_ids <- selected_common_ids
    
    if (!all(plot_ids %in% gsea_table$ID)) {
      stop("当前基因的 GSEA 结果缺失部分共同通路，请检查结果一致性。")
    }
    
    selected_table <- gsea_table[
      gsea_table$ID %in% plot_ids,
      ,
      drop = FALSE
    ]
    
    selected_table <- selected_table[
      match(plot_ids, selected_table$ID),
      ,
      drop = FALSE
    ]
    
    if (nrow(selected_table) != length(plot_ids) ||
        any(is.na(selected_table$p.adjust)) ||
        any(selected_table$p.adjust >= padj_cutoff)) {
      stop("Figure 5 中存在非共同显著通路，与严格筛选逻辑不一致。")
    }
    
    summary_list[[gene]] <- selected_table[, c(
      "Biomarker",
      "ID",
      "Description",
      "setSize",
      "enrichmentScore",
      "NES",
      "pvalue",
      "p.adjust",
      "qvalue"
    )]
    
    write.csv(
      summary_list[[gene]],
      file = file.path(outdir, paste0("GSEA_", gene, "_final_selected_common_pathways.csv")),
      row.names = FALSE
    )
    
    png_file <- file.path(outdir, paste0("Fig5_", gene, ".png"))
    pdf_file <- file.path(outdir, paste0("Fig5_", gene, ".pdf"))
    
    p <- plot_gsea_aligned(
      geneList = rank_score,
      gene_sets = hallmark_gene_sets,
      plot_ids = plot_ids,
      plot_colors = plot_colors,
      gene = gene,
      png_file = png_file,
      pdf_file = pdf_file,
      width = 10.5,
      height = 6.2
    )
    
    cat("PNG已保存：", png_file, "\n")
    cat("PDF已保存：", pdf_file, "\n")
    cat("完成绘图：", gene, "\n")
    
  }, error = function(e) {
    
    cat("基因 ", gene, " 绘图失败。\n")
    cat("错误信息：", conditionMessage(e), "\n")
    
  })
}

############################################################
## 13. 汇总最终展示通路结果
############################################################

if (length(summary_list) > 0) {
  
  summary_df <- do.call(rbind, summary_list)
  
  write.csv(
    summary_df,
    file = file.path(outdir, "GSEA_3_biomarkers_Figure5_selected_pathways_gene_specific.csv"),
    row.names = FALSE
  )
  
  cat("\n已输出最终汇总文件：GSEA_3_biomarkers_Figure5_selected_pathways_gene_specific.csv\n")
  
} else {
  
  warning("没有任何基因成功完成最终绘图。")
  
}

############################################################
## 14. 检查输出文件
############################################################

cat("\nFigure5目录下生成的主要图片文件：\n")
print(list.files(outdir, pattern = "^Fig5_.*\\.(png|pdf)$"))

# 保存 R 会话信息，便于论文复现
capture.output(
  sessionInfo(),
  file = file.path(outdir, "Figure5_GSEA_sessionInfo.txt")
)

cat("\n全部 Figure 5 分析结束。\n")
cat("结果输出目录：", outdir, "\n")
cat("显著性标准：BH-adjusted P < ", padj_cutoff, "\n", sep = "")
cat("三基因共同显著通路数量：", length(common_sig_ids), "\n")
cat("Figure 5 实际展示通路数量：", length(selected_common_ids), "\n")