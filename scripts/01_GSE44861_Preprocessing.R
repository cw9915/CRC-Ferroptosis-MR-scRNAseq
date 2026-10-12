

Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list = ls())

############################################################
## 1. 加载 R 包
############################################################

library(GEOquery)
library(dplyr)
library(tibble)
library(data.table)
library(stringr)
library(limma)

############################################################
## 2. 设置工作目录
############################################################

setwd("C:/Users/92066/Desktop/单细胞孟德尔复现/Figure1")
getwd()

############################################################
## 3. 设置输入文件
############################################################

gse_file <- "GSE44861_series_matrix.txt.gz"
gpl_file <- "GPL3921-25447.txt"

if (!file.exists(gse_file)) {
  stop("找不到 GSE44861_series_matrix.txt.gz，请检查文件名或路径。")
}

if (!file.exists(gpl_file)) {
  stop("找不到 GPL3921-25447.txt，请检查文件名或路径。")
}

############################################################
## 4. 读取 GSE44861 series matrix
############################################################

gset <- getGEO(
  filename = gse_file,
  destdir = ".",
  getGPL = FALSE
)

if (is.list(gset)) {
  e1 <- gset[[1]]
} else {
  e1 <- gset
}

############################################################
## 5. 提取表达矩阵和临床信息
############################################################

exp <- exprs(e1)
pdata <- pData(e1)

cat("原始表达矩阵维度：\n")
print(dim(exp))

cat("临床信息维度：\n")
print(dim(pdata))

############################################################
## 6. 保存原始对象和临床信息
############################################################

save(gset, file = "GSE44861_raw_gset.rdata")

write.csv(
  pdata,
  file = "GSE44861_pdata_raw.csv",
  row.names = TRUE
)

############################################################
## 7. 读取 GPL3921 平台注释文件
## 适用于 Affymetrix 老平台注释格式
############################################################

gpl_lines <- readLines(gpl_file, warn = FALSE)

cat("GPL 文件前 20 行：\n")
print(gpl_lines[1:min(20, length(gpl_lines))])

############################################################
## 7.1 找到真正表头
## 从你的截图看，真实表头以 ID GB_ACC SPOT_ID 开头
############################################################

header_line <- grep(
  pattern = "^ID\\s+GB_ACC|^ID\\tGB_ACC|^ID,GB_ACC",
  x = gpl_lines
)

if (length(header_line) == 0) {
  header_line <- grep(
    pattern = "^ID\\s+.*Gene Symbol|^ID\\t.*Gene Symbol",
    x = gpl_lines
  )
}

if (length(header_line) == 0) {
  cat("未找到表头。GPL 前 60 行如下：\n")
  print(gpl_lines[1:min(60, length(gpl_lines))])
  stop("没有找到 GPL3921 的真实表头，请检查文件格式。")
}

header_line <- header_line[1]

cat("检测到 GPL 表头所在行：", header_line, "\n")
cat("GPL 表头内容：\n")
print(gpl_lines[header_line])

############################################################
## 7.2 从表头开始读取 GPL 表格
############################################################

gpl_text <- paste(
  gpl_lines[header_line:length(gpl_lines)],
  collapse = "\n"
)

gpl <- data.table::fread(
  input = gpl_text,
  sep = "\t",
  header = TRUE,
  data.table = FALSE,
  fill = TRUE,
  quote = ""
)

############################################################
## 7.3 如果 tab 读取列数过少，尝试自动分隔
############################################################

if (ncol(gpl) <= 2) {
  cat("tab 分隔读取列数过少，尝试自动分隔读取。\n")
  
  gpl <- data.table::fread(
    input = gpl_text,
    header = TRUE,
    data.table = FALSE,
    fill = TRUE,
    quote = ""
  )
}

cat("GPL3921 注释文件维度：\n")
print(dim(gpl))

cat("GPL3921 注释文件列名：\n")
print(colnames(gpl))

write.csv(
  gpl,
  file = "GPL3921_annotation_raw.csv",
  row.names = FALSE
)

if (ncol(gpl) <= 2) {
  stop("GPL 文件仍然没有正确分列，请检查是否为制表符分隔文件。")
}

############################################################
## 8. 构建探针 ID - Gene Symbol 注释表
############################################################

colnames(gpl) <- trimws(colnames(gpl))

if (!"ID" %in% colnames(gpl)) {
  colnames(gpl)[1] <- "ID"
}

cat("清理后的 GPL 列名：\n")
print(colnames(gpl))

############################################################
## 8.1 自动识别 Gene Symbol 列
############################################################

symbol_col <- NULL

if ("Gene Symbol" %in% colnames(gpl)) {
  symbol_col <- "Gene Symbol"
} else if ("Gene.Symbol" %in% colnames(gpl)) {
  symbol_col <- "Gene.Symbol"
} else if ("GENE_SYMBOL" %in% colnames(gpl)) {
  symbol_col <- "GENE_SYMBOL"
} else if ("Symbol" %in% colnames(gpl)) {
  symbol_col <- "Symbol"
} else if ("gene_symbol" %in% colnames(gpl)) {
  symbol_col <- "gene_symbol"
}

if (is.null(symbol_col)) {
  cat("当前 GPL 文件列名如下：\n")
  print(colnames(gpl))
  stop("没有找到 Gene Symbol 列。请确认列名是否被 R 改写。")
}

cat("使用的 Gene Symbol 列为：", symbol_col, "\n")

gene <- gpl %>%
  dplyr::select(
    ID,
    Symbol = dplyr::all_of(symbol_col)
  )

############################################################
## 8.2 清理 Gene Symbol
############################################################

gene <- gene %>%
  mutate(
    ID = as.character(ID),
    Symbol = as.character(Symbol)
  ) %>%
  filter(!is.na(Symbol)) %>%
  filter(Symbol != "") %>%
  filter(Symbol != "---") %>%
  filter(Symbol != "NA")

## 多个 Symbol 时保留第一个
gene$Symbol <- gsub("\"", "", gene$Symbol)
gene$Symbol <- sapply(
  strsplit(gene$Symbol, "///|//|;|,"),
  function(x) trimws(x[1])
)

gene <- gene %>%
  filter(!is.na(Symbol)) %>%
  filter(Symbol != "") %>%
  distinct(ID, Symbol)

cat("有效探针-基因注释数量：\n")
print(nrow(gene))

cat("注释到的唯一基因数量：\n")
print(length(unique(gene$Symbol)))

write.csv(
  gene,
  file = "GPL3921_probe_gene_annotation_used.csv",
  row.names = FALSE
)

############################################################
## 9. 表达矩阵探针 ID 转为 Gene Symbol
############################################################

df_exp <- as.data.frame(exp) %>%
  rownames_to_column("ID") %>%
  mutate(ID = as.character(ID))

exp_anno <- merge(
  x = gene,
  y = df_exp,
  by = "ID"
)

cat("合并注释后的表达矩阵维度：\n")
print(dim(exp_anno))

if (nrow(exp_anno) == 0) {
  cat("表达矩阵前 10 个探针 ID：\n")
  print(head(rownames(exp), 10))
  
  cat("GPL 注释前 10 个探针 ID：\n")
  print(head(gene$ID, 10))
  
  stop("探针注释与表达矩阵没有匹配成功，请检查 GPL3921 的 ID 是否与表达矩阵行名一致。")
}

############################################################
## 10. 多探针对应同一基因时取平均表达
############################################################

expr_gene <- exp_anno %>%
  select(-ID) %>%
  group_by(Symbol) %>%
  summarise(
    across(
      .cols = everything(),
      .fns = function(x) mean(as.numeric(x), na.rm = TRUE)
    ),
    .groups = "drop"
  )

expr_gene_mat <- expr_gene %>%
  column_to_rownames("Symbol") %>%
  as.matrix()

mode(expr_gene_mat) <- "numeric"

cat("转换为 Gene Symbol 后表达矩阵维度：\n")
print(dim(expr_gene_mat))

############################################################
## 11. 检查是否需要 log2 转换
############################################################

cat("表达矩阵最大值：", max(expr_gene_mat, na.rm = TRUE), "\n")
cat("表达矩阵分布：\n")

print(
  summary(
    as.vector(
      expr_gene_mat[
        1:min(1000, nrow(expr_gene_mat)),
        1:min(50, ncol(expr_gene_mat))
      ]
    )
  )
)

if (max(expr_gene_mat, na.rm = TRUE) > 50) {
  expr_gene_mat <- log2(expr_gene_mat + 1)
  cat("表达矩阵已进行 log2(x+1) 转换。\n")
} else {
  cat("表达矩阵看起来已经是 log2 格式，不再转换。\n")
}

############################################################
## 12. 保存初步 Gene Symbol 表达矩阵
############################################################

expr_gene_df <- as.data.frame(expr_gene_mat) %>%
  rownames_to_column("Gene")

write.csv(
  expr_gene_df,
  file = "GSE44861_gene_symbol_expression.csv",
  row.names = FALSE
)

saveRDS(
  expr_gene_mat,
  file = "GSE44861_gene_symbol_expression.rds"
)

############################################################
## 13. 整理并保存样本信息
############################################################

sample_info <- pdata %>%
  rownames_to_column("Sample")

write.csv(
  sample_info,
  file = "GSE44861_sample_information.csv",
  row.names = FALSE
)

cat("sample_info 列名：\n")
print(colnames(sample_info))

############################################################
## 14. 根据 Sample_title 精确提取 GSE44861 分组
## N = Normal
## T = Tumor
## 兼容格式：
## 10167 N
## 10167 T
## 12385-T
## 13687N
## 13687T
############################################################

if (!"title" %in% colnames(sample_info)) {
  stop("sample_info 中没有 title 列，请检查 pData(e1) 的列名。")
}

cat("前 20 个样本 title：\n")
print(head(sample_info$title, 20))

group_final <- sample_info %>%
  mutate(
    title_clean = trimws(as.character(title)),
    title_clean = gsub("\\s+$", "", title_clean),
    
    Group_raw = case_when(
      grepl("[[:space:]-]N$", title_clean, ignore.case = TRUE) ~ "N",
      grepl("[[:space:]-]T$", title_clean, ignore.case = TRUE) ~ "T",
      grepl("N$", title_clean, ignore.case = TRUE) ~ "N",
      grepl("T$", title_clean, ignore.case = TRUE) ~ "T",
      TRUE ~ NA_character_
    ),
    
    Group = case_when(
      Group_raw == "N" ~ "Normal",
      Group_raw == "T" ~ "Tumor",
      TRUE ~ NA_character_
    )
  ) %>%
  select(Sample, title, Group)

############################################################
## 14.1 检查分组结果
############################################################

cat("GSE44861 分组情况：\n")
print(table(group_final$Group, useNA = "ifany"))

if (any(is.na(group_final$Group))) {
  cat("以下样本未能识别分组：\n")
  print(group_final[is.na(group_final$Group), ])
  stop("存在未识别分组的样本，请检查 title 格式。")
}

if (!all(c("Normal", "Tumor") %in% group_final$Group)) {
  stop("分组中没有同时包含 Normal 和 Tumor，请检查 title 分组。")
}

############################################################
## 14.2 保存最终分组文件
############################################################

write.csv(
  group_final,
  file = "GSE44861_group_final.csv",
  row.names = FALSE
)

cat("最终分组文件已保存：GSE44861_group_final.csv\n")

############################################################
## 14.3 保存后续差异分析用分组文件
############################################################

group_for_deg <- group_final %>%
  select(Sample, Group)

write.csv(
  group_for_deg,
  file = "02_group.csv",
  row.names = FALSE
)

cat("差异分析分组文件已保存：02_group.csv\n")

############################################################
## 14.4 检查表达矩阵与分组文件是否完全匹配
############################################################

common_samples <- intersect(colnames(expr_gene_mat), group_for_deg$Sample)

cat("表达矩阵样本数：", ncol(expr_gene_mat), "\n")
cat("分组文件样本数：", nrow(group_for_deg), "\n")
cat("共同样本数：", length(common_samples), "\n")

if (length(common_samples) != ncol(expr_gene_mat)) {
  cat("表达矩阵中不在分组文件内的样本：\n")
  print(setdiff(colnames(expr_gene_mat), group_for_deg$Sample))
  
  cat("分组文件中不在表达矩阵内的样本：\n")
  print(setdiff(group_for_deg$Sample, colnames(expr_gene_mat)))
  
  stop("表达矩阵与分组文件样本不完全一致，请检查。")
}

############################################################
## 14.5 按分组文件顺序重排表达矩阵并重新保存
############################################################

expr_gene_mat <- expr_gene_mat[, group_for_deg$Sample]

expr_gene_df <- as.data.frame(expr_gene_mat) %>%
  rownames_to_column("Gene")

write.csv(
  expr_gene_df,
  file = "GSE44861_gene_symbol_expression.csv",
  row.names = FALSE
)

saveRDS(
  expr_gene_mat,
  file = "GSE44861_gene_symbol_expression.rds"
)

############################################################
## 15. 输出完成提示
############################################################

cat("\n==============================\n")
cat("GSE44861 表达矩阵整理完成。\n")
cat("分组情况：\n")
print(table(group_for_deg$Group))
cat("\n主要输出文件：\n")
cat("1. GSE44861_gene_symbol_expression.csv\n")
cat("2. GSE44861_gene_symbol_expression.rds\n")
cat("3. GSE44861_pdata_raw.csv\n")
cat("4. GSE44861_sample_information.csv\n")
cat("5. GSE44861_group_final.csv\n")
cat("6. 02_group.csv\n")
cat("7. GPL3921_probe_gene_annotation_used.csv\n")
cat("8. GPL3921_annotation_raw.csv\n")
cat("==============================\n")