############################################################
## Figure 2A-B
## GSE44861 colorectal cancer Tumor vs Normal
## A: Volcano plot
## B: Top DEGs heatmap without density annotation
############################################################

Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list = ls())

############################################################
## 1. 加载 R 包
############################################################

library(limma)
library(dplyr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(ComplexHeatmap)
library(circlize)
library(grid)

############################################################
## 2. 设置路径
############################################################

setwd("C:/Users/92066/Desktop/单细胞孟德尔复现/Figure2")

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

expr_file  <- "C:/Users/92066/Desktop/单细胞孟德尔复现/Figure1/GSE44861_gene_symbol_expression.csv"
group_file <- "C:/Users/92066/Desktop/单细胞孟德尔复现/Figure1/02_group.csv"

if (!file.exists(expr_file)) {
  stop("找不到表达矩阵文件：GSE44861_gene_symbol_expression.csv")
}

if (!file.exists(group_file)) {
  stop("找不到分组文件：02_group.csv")
}

############################################################
## 3. 读取表达矩阵和分组文件
############################################################

expr <- read.csv(expr_file, check.names = FALSE)
group <- read.csv(group_file, check.names = FALSE)

colnames(expr)[1] <- "Gene"
colnames(group)[1:2] <- c("Sample", "Group")

expr <- expr %>%
  distinct(Gene, .keep_all = TRUE)

expr_mat <- expr %>%
  column_to_rownames("Gene") %>%
  as.matrix()

mode(expr_mat) <- "numeric"

group <- group %>%
  filter(Sample %in% colnames(expr_mat))

expr_mat <- expr_mat[, group$Sample]

group$Group <- factor(group$Group, levels = c("Normal", "Tumor"))

cat("表达矩阵维度：\n")
print(dim(expr_mat))

cat("分组情况：\n")
print(table(group$Group))

if (any(is.na(group$Group))) {
  stop("分组中存在 NA，请检查 02_group.csv。")
}

if (!all(c("Normal", "Tumor") %in% group$Group)) {
  stop("分组中没有同时包含 Normal 和 Tumor。")
}

############################################################
## 4. 检查表达矩阵是否需要 log2 转换
############################################################

cat("表达矩阵最大值：", max(expr_mat, na.rm = TRUE), "\n")

cat("表达矩阵分布：\n")
print(
  summary(
    as.vector(
      expr_mat[
        1:min(1000, nrow(expr_mat)),
        1:min(50, ncol(expr_mat))
      ]
    )
  )
)

if (max(expr_mat, na.rm = TRUE) > 50) {
  expr_mat <- log2(expr_mat + 1)
  cat("表达矩阵已进行 log2(x+1) 转换。\n")
} else {
  cat("表达矩阵看起来已经是 log2 格式，不再转换。\n")
}

############################################################
## 5. 芯片数据标准化
############################################################

expr_mat <- normalizeBetweenArrays(expr_mat)

############################################################
## 6. limma 差异分析：Tumor vs Normal
############################################################

design <- model.matrix(~ 0 + group$Group)
colnames(design) <- levels(group$Group)

fit <- lmFit(expr_mat, design)

contrast.matrix <- makeContrasts(
  Tumor_vs_Normal = Tumor - Normal,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

deg <- topTable(
  fit2,
  coef = "Tumor_vs_Normal",
  number = Inf,
  adjust.method = "BH"
)

deg <- deg %>%
  rownames_to_column("Gene") %>%
  mutate(
    change = case_when(
      adj.P.Val < 0.05 & logFC > 0.5  ~ "Up",
      adj.P.Val < 0.05 & logFC < -0.5 ~ "Down",
      TRUE ~ "Not"
    ),
    FDR_plot = pmax(adj.P.Val, .Machine$double.xmin)
  )

write.csv(
  deg,
  file = "results/Fig2AB_All_DEGs_Tumor_vs_Normal.csv",
  row.names = FALSE
)

cat("上调基因数：", sum(deg$change == "Up"), "\n")
cat("下调基因数：", sum(deg$change == "Down"), "\n")
cat("总差异基因数：", sum(deg$change != "Not"), "\n")

############################################################
## 7. Figure 2A 火山图
############################################################

deg$label <- ""

top_up <- deg %>%
  filter(change == "Up") %>%
  arrange(adj.P.Val) %>%
  head(10) %>%
  pull(Gene)

top_down <- deg %>%
  filter(change == "Down") %>%
  arrange(adj.P.Val) %>%
  head(10) %>%
  pull(Gene)

deg$label[deg$Gene %in% c(top_up, top_down)] <- deg$Gene

p_volcano <- ggplot(
  deg,
  aes(x = logFC, y = -log10(FDR_plot))
) +
  geom_point(
    aes(color = change),
    size = 1.8,
    alpha = 0.78
  ) +
  scale_color_manual(
    values = c(
      "Up" = "#E64B35",
      "Down" = "#0B4F9C",
      "Not" = "grey75"
    ),
    breaks = c("Up", "Down", "Not")
  ) +
  geom_vline(
    xintercept = c(-0.5, 0.5),
    linetype = "dashed",
    linewidth = 0.7,
    color = "grey35"
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    linewidth = 0.7,
    color = "grey35"
  ) +
  geom_text_repel(
    aes(label = label),
    size = 4.8,
    fontface = "bold",
    max.overlaps = 40,
    box.padding = 0.45,
    point.padding = 0.35,
    segment.color = "grey45",
    show.legend = FALSE
  ) +
  labs(
    x = expression(Log[2]~"(FC)"),
    y = expression(-Log[10]~"(FDR)"),
    color = ""
  ) +
  theme_bw(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(linewidth = 1.25, color = "black"),
    axis.ticks = element_line(linewidth = 1.1, color = "black"),
    axis.text = element_text(size = 15, face = "bold", color = "black"),
    axis.title = element_text(size = 19, face = "bold", color = "black"),
    legend.text = element_text(size = 15, face = "bold", color = "black"),
    legend.position = "right"
  )

ggsave(
  filename = "figures/Fig2A_volcano.pdf",
  plot = p_volcano,
  width = 7.0,
  height = 6.2
)

ggsave(
  filename = "figures/Fig2A_volcano.png",
  plot = p_volcano,
  width = 7.0,
  height = 6.2,
  dpi = 600
)

############################################################
## 8. Figure 2B Top DEGs heatmap
## 不做 density，只保留 Sample 分组色条
############################################################

############################################################
## 8.1 选择 top DEGs
## 上调 25 个 + 下调 25 个
############################################################

top_up_genes <- deg %>%
  filter(change == "Up") %>%
  arrange(adj.P.Val) %>%
  head(25) %>%
  pull(Gene)

top_down_genes <- deg %>%
  filter(change == "Down") %>%
  arrange(adj.P.Val) %>%
  head(25) %>%
  pull(Gene)

top_genes <- unique(c(top_up_genes, top_down_genes))

if (length(top_genes) < 2) {
  stop("显著差异基因少于 2 个，无法绘制热图。")
}

write.csv(
  data.frame(Gene = top_genes),
  file = "results/Fig2B_top_DEGs.csv",
  row.names = FALSE
)

############################################################
## 8.2 按分组排序样本
############################################################

group_order <- order(group$Group)

group_ordered <- group[group_order, ]
expr_mat_ordered <- expr_mat[, group_ordered$Sample]

############################################################
## 8.3 提取表达矩阵并做 Z-score
############################################################

heat_mat <- expr_mat_ordered[top_genes, ]

heat_mat_z <- t(scale(t(heat_mat)))
heat_mat_z[is.na(heat_mat_z)] <- 0

## 截断，增强颜色对比
heat_mat_z[heat_mat_z > 4] <- 4
heat_mat_z[heat_mat_z < -4] <- -4

############################################################
## 8.4 设置颜色
############################################################

heat_col <- colorRamp2(
  c(-4, 0, 4),
  c("#0B4F9C", "white", "#E64B35")
)

group_col <- c(
  "Normal" = "#55C1D9",
  "Tumor"  = "#E64B35"
)

############################################################
## 8.5 构建顶部 Sample 分组色条
## 右侧标注 Sample
############################################################

top_ha <- HeatmapAnnotation(
  Sample = group_ordered$Group,
  
  col = list(
    Sample = group_col
  ),
  
  annotation_name_side = "right",
  annotation_name_rot = 0,
  
  annotation_name_gp = gpar(
    fontsize = 16,
    fontface = "bold",
    col = "black"
  ),
  
  simple_anno_size = unit(6, "mm"),
  
  annotation_legend_param = list(
    Sample = list(
      title = "Sample",
      title_gp = gpar(fontsize = 15, fontface = "bold"),
      labels_gp = gpar(fontsize = 14, fontface = "bold")
    )
  )
)

############################################################
## 8.6 构建热图
############################################################

ht <- Heatmap(
  heat_mat_z,
  
  name = "Expr",
  col = heat_col,
  
  top_annotation = top_ha,
  
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  
  show_column_names = FALSE,
  show_row_names = TRUE,
  row_names_side = "right",
  
  row_names_gp = gpar(
    fontsize = 13,
    fontface = "bold",
    col = "black"
  ),
  
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  
  border = FALSE,
  
  column_title = NULL,
  row_title = NULL,
  
  heatmap_legend_param = list(
    title = "Expr",
    title_gp = gpar(fontsize = 15, fontface = "bold"),
    labels_gp = gpar(fontsize = 14, fontface = "bold"),
    legend_height = unit(45, "mm")
  )
)

############################################################
## 8.7 输出热图
############################################################

pdf(
  file = "figures/Fig2B_top_DEGs_heatmap.pdf",
  width = 7.8,
  height = 9.5
)

draw(
  ht,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  padding = unit(c(5, 5, 5, 5), "mm")
)

dev.off()

png(
  filename = "figures/Fig2B_top_DEGs_heatmap.png",
  width = 2600,
  height = 3200,
  res = 300
)

draw(
  ht,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  padding = unit(c(5, 5, 5, 5), "mm")
)

dev.off()

############################################################
## 9. 保存显著差异基因
############################################################

write.csv(
  deg %>% filter(change != "Not"),
  file = "results/Fig2AB_Significant_DEGs.csv",
  row.names = FALSE
)

############################################################
## 10. 输出完成提示
############################################################

cat("\n==============================\n")
cat("Figure 2A-B 完成。\n")
cat("分组情况：\n")
print(table(group$Group))
cat("\n差异分析结果：\n")
cat("上调基因数：", sum(deg$change == "Up"), "\n")
cat("下调基因数：", sum(deg$change == "Down"), "\n")
cat("总差异基因数：", sum(deg$change != "Not"), "\n")
cat("\n输出文件：\n")
cat("figures/Fig2A_volcano.pdf/png\n")
cat("figures/Fig2B_top_DEGs_heatmap.pdf/png\n")
cat("results/Fig2AB_All_DEGs_Tumor_vs_Normal.csv\n")
cat("results/Fig2AB_Significant_DEGs.csv\n")
cat("results/Fig2B_top_DEGs.csv\n")
cat("==============================\n")

############################################################
## Figure 2C-D
## C: GO enrichment analysis of DEGs
## D: KEGG enrichment analysis of DEGs
############################################################
options(timeout = 600)
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list = ls())

############################################################
## 1. 加载 R 包
############################################################

library(dplyr)
library(tibble)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(stringr)
library(forcats)

############################################################
## 2. 设置路径
############################################################

setwd("C:/Users/92066/Desktop/单细胞孟德尔复现/Figure2")

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

deg_file <- "results/Fig2AB_Significant_DEGs.csv"

if (!file.exists(deg_file)) {
  stop("找不到 results/Fig2AB_Significant_DEGs.csv，请先完成 Figure 2A-B。")
}

############################################################
## 3. 读取差异基因
############################################################

deg <- read.csv(deg_file, check.names = FALSE)

if (!"Gene" %in% colnames(deg)) {
  stop("差异基因文件中没有 Gene 列，请检查 Fig2AB_Significant_DEGs.csv。")
}

DEGs <- unique(deg$Gene)
DEGs <- DEGs[!is.na(DEGs) & DEGs != ""]

cat("用于富集分析的显著 DEGs 数量：", length(DEGs), "\n")

if (length(DEGs) < 10) {
  stop("显著差异基因少于 10 个，不建议进行 GO/KEGG 富集分析。")
}

############################################################
## 4. Symbol 转 Entrez ID
############################################################

gene_df <- bitr(
  DEGs,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

gene_df <- gene_df %>%
  distinct(SYMBOL, ENTREZID)

cat("成功转换 Entrez ID 的基因数：", nrow(gene_df), "\n")

write.csv(
  gene_df,
  file = "results/Fig2CD_DEGs_EntrezID.csv",
  row.names = FALSE
)

if (nrow(gene_df) < 10) {
  stop("成功转换 Entrez ID 的基因少于 10 个，请检查基因名格式。")
}

############################################################
## 5. Figure 2C: GO enrichment analysis
############################################################

go_all <- enrichGO(
  gene          = gene_df$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff  = 1,
  qvalueCutoff  = 1,
  readable      = TRUE
)

go_res <- as.data.frame(go_all)

write.csv(
  go_res,
  file = "results/Fig2C_GO_all_results.csv",
  row.names = FALSE
)

############################################################
## 5.1 筛选 GO 结果
############################################################

go_plot_df <- go_res %>%
  filter(pvalue < 0.05) %>%
  filter(ONTOLOGY %in% c("BP", "CC", "MF"))

if (nrow(go_plot_df) == 0) {
  cat("没有 pvalue < 0.05 的 GO 结果，改用 pvalue 最小的 GO 条目作图。\n")
  
  go_plot_df <- go_res %>%
    filter(ONTOLOGY %in% c("BP", "CC", "MF")) %>%
    arrange(pvalue)
}

## 每组显示 8 个 GO term
go_plot_df <- go_plot_df %>%
  group_by(ONTOLOGY) %>%
  arrange(pvalue, .by_group = TRUE) %>%
  slice_head(n = 8) %>%
  ungroup()

write.csv(
  go_plot_df,
  file = "results/Fig2C_GO_top_terms_for_plot.csv",
  row.names = FALSE
)

############################################################
## 5.2 整理 GO 作图数据
############################################################

go_plot_df <- go_plot_df %>%
  mutate(
    ONTOLOGY = factor(
      ONTOLOGY,
      levels = c("BP", "CC", "MF")
    )
  )

## BP / CC / MF 配色
ontology_cols <- c(
  "BP" = "#69B9DD",
  "CC" = "#7B7FB8",
  "MF" = "#E49A77"
)

## 只缩写最长的那个 term，其他 term 不改、不换行
go_plot_df <- go_plot_df %>%
  mutate(
    Description_show = as.character(Description),
    Description_show = gsub(
      "oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor",
      "oxidoreductase activity (CH-OH, NAD/NADP)",
      Description_show,
      fixed = TRUE
    )
  )

## 关键：先按 BP、CC、MF 分组，再按 Count 从高到低排列
go_plot_df <- go_plot_df %>%
  arrange(ONTOLOGY, desc(Count))

## 关键修正：TermID 必须设为 factor，并固定 levels，防止 Term_10 排到 Term_2 前面
go_plot_df$TermID <- paste0("Term_", seq_len(nrow(go_plot_df)))
go_plot_df$TermID <- factor(
  go_plot_df$TermID,
  levels = go_plot_df$TermID
)

## X 轴字体颜色与 BP/CC/MF 对应
x_label_cols <- ontology_cols[as.character(go_plot_df$ONTOLOGY)]

############################################################
## 5.3 绘制 GO 点线图
############################################################

p_go <- ggplot(
  go_plot_df,
  aes(
    x = TermID,
    y = Count,
    color = ONTOLOGY
  )
) +
  geom_segment(
    aes(
      x = TermID,
      xend = TermID,
      y = 0,
      yend = Count
    ),
    linewidth = 1.25,
    alpha = 0.95
  ) +
  geom_point(
    size = 4.8,
    alpha = 0.95
  ) +
  geom_text(
    aes(label = Count),
    vjust = -0.65,
    size = 4.6,
    fontface = "bold",
    color = "black"
  ) +
  scale_x_discrete(
    labels = go_plot_df$Description_show
  ) +
  scale_color_manual(
    values = ontology_cols,
    name = ""
  ) +
  scale_y_continuous(
    limits = c(0, max(go_plot_df$Count, na.rm = TRUE) * 1.25),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = "",
    y = "Count"
  ) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 16) +
  theme(
    panel.border = element_blank(),
    
    axis.line = element_line(
      color = "black",
      linewidth = 0.8
    ),
    axis.ticks = element_line(
      color = "black",
      linewidth = 0.8
    ),
    
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 12.5,
      face = "bold",
      color = x_label_cols,
      lineheight = 0.9,
      margin = margin(t = 3)
    ),
    axis.text.y = element_text(
      size = 13,
      face = "bold",
      color = "black"
    ),
    axis.title.y = element_text(
      size = 15,
      face = "bold",
      color = "black",
      margin = margin(r = 6)
    ),
    
    legend.position = "right",
    legend.text = element_text(
      size = 13,
      face = "bold",
      color = "black"
    ),
    legend.key.size = unit(0.45, "cm"),
    
    plot.margin = margin(
      t = 10,
      r = 12,
      b = 28,
      l = 12,
      unit = "pt"
    )
  )

ggsave(
  filename = "figures/Fig2C_GO_lollipop.pdf",
  plot = p_go,
  width = 10.5,
  height = 7.8
)

ggsave(
  filename = "figures/Fig2C_GO_lollipop.png",
  plot = p_go,
  width = 10.5,
  height = 7.8,
  dpi = 600
)
############################################################
## 6. Figure 2D: KEGG enrichment analysis
############################################################

kk <- enrichKEGG(
  gene          = gene_df$ENTREZID,
  organism      = "hsa",
  pvalueCutoff  = 1,
  qvalueCutoff  = 1,
  pAdjustMethod = "BH"
)

kk_readable <- setReadable(
  kk,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID"
)

kegg_res <- as.data.frame(kk_readable)

write.csv(
  kegg_res,
  file = "results/Fig2D_KEGG_all_results.csv",
  row.names = FALSE
)

############################################################
## 6.1 筛选 KEGG 结果
############################################################

kegg_plot_df <- kegg_res %>%
  filter(pvalue < 0.05)

if (nrow(kegg_plot_df) == 0) {
  cat("没有 pvalue < 0.05 的 KEGG 结果，改用 pvalue 最小的 KEGG 条目作图。\n")
  
  kegg_plot_df <- kegg_res %>%
    arrange(pvalue)
}

kegg_plot_df <- kegg_plot_df %>%
  arrange(pvalue) %>%
  slice_head(n = 15)

write.csv(
  kegg_plot_df,
  file = "results/Fig2D_KEGG_top_terms_for_plot.csv",
  row.names = FALSE
)

############################################################
## 6.2 处理 GeneRatio
############################################################

kegg_plot_df <- kegg_plot_df %>%
  mutate(
    GeneRatio_num = sapply(
      GeneRatio,
      function(x) {
        tmp <- unlist(strsplit(x, "/"))
        as.numeric(tmp[1]) / as.numeric(tmp[2])
      }
    ),
    Description_wrapped = str_wrap(Description, width = 45)
  )

kegg_plot_df$Description_wrapped <- factor(
  kegg_plot_df$Description_wrapped,
  levels = rev(kegg_plot_df$Description_wrapped)
)

############################################################
## 6.3 绘制 KEGG 气泡图
############################################################

p_kegg <- ggplot(
  kegg_plot_df,
  aes(
    x = GeneRatio_num,
    y = Description_wrapped,
    size = Count,
    color = pvalue
  )
) +
  geom_point(alpha = 0.95) +
  scale_color_gradientn(
    colours = c(
      "#E64B35",
      "#F39B7F",
      "#91D1C2",
      "#3C5488"
    ),
    name = "P value",
    trans = "reverse"
  ) +
  scale_size_continuous(
    range = c(5, 12),
    name = "Count"
  ) +
  labs(
    x = "Gene ratio",
    y = ""
  ) +
  theme_bw(base_size = 18) +
  theme(
    panel.grid.major = element_line(color = "grey88", linewidth = 0.45),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 1.1),
    
    axis.text.x = element_text(
      size = 16,
      face = "bold",
      color = "black"
    ),
    axis.text.y = element_text(
      size = 15,
      face = "bold",
      color = "black"
    ),
    axis.title.x = element_text(
      size = 21,
      face = "bold",
      color = "black"
    ),
    
    legend.title = element_text(
      size = 16,
      face = "bold",
      color = "black"
    ),
    legend.text = element_text(
      size = 14,
      face = "bold",
      color = "black"
    ),
    
    plot.margin = margin(
      t = 20,
      r = 30,
      b = 20,
      l = 20,
      unit = "pt"
    )
  )

ggsave(
  filename = "figures/Fig2D_KEGG_bubble.pdf",
  plot = p_kegg,
  width = 10.5,
  height = 7.8
)

ggsave(
  filename = "figures/Fig2D_KEGG_bubble.png",
  plot = p_kegg,
  width = 10.5,
  height = 7.8,
  dpi = 600
)

############################################################
## 7. 输出完成提示
############################################################

cat("\n==============================\n")
cat("Figure 2C-D 富集分析完成。\n")
cat("用于富集分析的 DEGs 数量：", length(DEGs), "\n")
cat("成功转换 Entrez ID 的基因数：", nrow(gene_df), "\n")
cat("GO 作图条目数：", nrow(go_plot_df), "\n")
cat("KEGG 作图条目数：", nrow(kegg_plot_df), "\n")
cat("\n输出文件：\n")
cat("figures/Fig2C_GO_lollipop.pdf/png\n")
cat("figures/Fig2D_KEGG_bubble.pdf/png\n")
cat("results/Fig2C_GO_all_results.csv\n")
cat("results/Fig2C_GO_top_terms_for_plot.csv\n")
cat("results/Fig2D_KEGG_all_results.csv\n")
cat("results/Fig2D_KEGG_top_terms_for_plot.csv\n")
cat("==============================\n")