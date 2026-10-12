#!/usr/bin/env Rscript

#========================================================
# 0. 基础设置
#========================================================

Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list = ls())

#========================================================
# 1. 检查并安装所需 R 包
#========================================================

required_packages <- c("data.table")

missing <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing)) {
  install.packages(missing, repos = "https://cloud.r-project.org")
}

library(data.table)

#========================================================
# 2. 设置工作目录
#========================================================

workdir <- "C:/Users/92066/Desktop/单细胞孟德尔复现/MR"
setwd(workdir)

cat("当前工作目录：\n")
cat(getwd(), "\n\n")

#========================================================
# 3. 输入文件
#========================================================

target_gene_file <- file.path(workdir, "intersect_genes.csv")

eqtl_file <- file.path(
  workdir,
  "eQTLGen",
  "2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz"
)

af_file <- file.path(
  workdir,
  "eQTLGen",
  "2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz"
)

#========================================================
# 4. 检查文件是否存在
#========================================================

if (!file.exists(target_gene_file)) {
  stop("目标基因文件不存在：", target_gene_file)
}

if (!file.exists(eqtl_file)) {
  stop("eQTLGen 文件不存在：", eqtl_file)
}

if (!file.exists(af_file)) {
  stop("Allele frequency 文件不存在：", af_file)
}

#========================================================
# 5. 定义函数：清洗基因名
#========================================================

clean_gene <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("\ufeff", "", x, fixed = TRUE)
  x <- gsub("\u00A0", "", x, fixed = TRUE)
  x <- gsub("\"", "", x, fixed = TRUE)
  x <- gsub("'", "", x, fixed = TRUE)
  x <- gsub("\\s+", "", x)
  x <- toupper(x)
  return(x)
}

#========================================================
# 6. 定义函数：自动识别基因列
#========================================================

resolve_gene_col <- function(dt) {
  
  nms <- names(dt)
  
  hits <- c(
    "gene", "Gene", "GENE",
    "symbol", "Symbol", "SYMBOL",
    "GeneSymbol", "genesymbol", "gene_symbol",
    "HGNC", "hgnc",
    "exposure"
  )
  
  hit <- hits[hits %in% nms]
  
  if (length(hit) > 0) {
    return(hit[1])
  }
  
  idx <- grep("gene|symbol|hgnc|exposure", nms, ignore.case = TRUE)
  
  if (length(idx) == 0) {
    stop("无法识别基因列。当前列名为：", paste(nms, collapse = ", "))
  }
  
  return(nms[idx[1]])
}

#========================================================
# 7. 读取目标基因文件
#========================================================
# 说明：
# 你的 intersect_genes.csv 是标准单列表格：
# Gene
# AQP8
# PRDX6
# ...
# 因此这里强制用 fread() 读取，不再用 readLines()
#========================================================

read_target_genes <- function(path) {
  
  dt <- fread(path)
  
  cat("目标基因文件原始列名：\n")
  print(names(dt))
  
  cat("\n目标基因文件前几行：\n")
  print(head(dt))
  
  if (ncol(dt) == 1) {
    setnames(dt, names(dt)[1], "gene")
  } else {
    gene_col <- resolve_gene_col(dt)
    setnames(dt, gene_col, "gene")
  }
  
  dt[, gene := as.character(gene)]
  dt[, gene := trimws(gene)]
  dt[, gene := gsub("\ufeff", "", gene, fixed = TRUE)]
  dt[, gene := gsub("\u00A0", "", gene, fixed = TRUE)]
  dt[, gene := gsub("\"", "", gene, fixed = TRUE)]
  dt[, gene := gsub("'", "", gene, fixed = TRUE)]
  
  # 去掉可能残留的表头
  dt <- dt[
    !is.na(gene) &
      nzchar(gene) &
      !tolower(gene) %in% c("gene", "symbol", "genesymbol", "gene_symbol", "hgnc")
  ]
  
  dt[, source := "intersect_genes"]
  dt[, selection_stage := "start_intersect_genes"]
  dt[, logFC := NA_real_]
  dt[, adj.P.Val := NA_real_]
  dt[, expression_trend := NA_character_]
  dt[, gene_clean := clean_gene(gene)]
  
  dt <- unique(dt, by = "gene_clean")
  
  return(dt)
}

target_genes <- read_target_genes(target_gene_file)

fwrite(
  target_genes,
  file.path(workdir, "standardized_target_genes.tsv"),
  sep = "\t"
)

cat("\n清洗后目标基因数量：", nrow(target_genes), "\n")
cat("清洗后目标基因前 30 个：\n")
print(head(target_genes$gene_clean, 30))

#========================================================
# 8. Step 1：查看 eQTLGen 文件表头
#========================================================

message("\nStep 1: Inspect the eQTLGen header")

eqtl_header <- fread(
  eqtl_file,
  nrows = 5
)

cat("\neQTLGen 文件列名：\n")
print(names(eqtl_header))

cat("\neQTLGen 前 5 行：\n")
print(eqtl_header)

#========================================================
# 9. 自动识别 eQTLGen 中的 GeneSymbol 列
#========================================================

preferred_gene_cols <- grep(
  "genesymbol|gene_symbol|hgnc|symbol",
  names(eqtl_header),
  ignore.case = TRUE,
  value = TRUE
)

fallback_gene_cols <- grep(
  "^gene$|gene",
  names(eqtl_header),
  ignore.case = TRUE,
  value = TRUE
)

eqtl_gene_col <- if (length(preferred_gene_cols) > 0) {
  preferred_gene_cols[1]
} else {
  fallback_gene_cols[1]
}

if (is.na(eqtl_gene_col) || length(eqtl_gene_col) == 0) {
  stop("eQTLGen 文件中无法识别基因列。")
}

cat("\neQTLGen 使用的基因列为：", eqtl_gene_col, "\n")

#========================================================
# 10. Step 2：读取完整 eQTLGen 文件
#========================================================

message("\nStep 2: Read full eQTLGen file")

eqtl_all <- fread(eqtl_file)

cat("eQTLGen 总记录数：", nrow(eqtl_all), "\n")
cat("eQTLGen 总列数：", ncol(eqtl_all), "\n")

eqtl_all[, GeneSymbol_raw := as.character(get(eqtl_gene_col))]
eqtl_all[, GeneSymbol_clean := clean_gene(GeneSymbol_raw)]

eqtl_genes_unique <- unique(
  eqtl_all[!is.na(GeneSymbol_clean) & nzchar(GeneSymbol_clean), GeneSymbol_clean]
)

cat("eQTLGen 清洗后唯一基因数：", length(eqtl_genes_unique), "\n")

#========================================================
# 11. 检查目标基因与 eQTLGen 的交集
#========================================================

message("\nStep 3: Match target genes with eQTLGen GeneSymbol")

target_gene_vec <- unique(target_genes$gene_clean)

matched_genes <- intersect(target_gene_vec, eqtl_genes_unique)
unmatched_genes <- setdiff(target_gene_vec, eqtl_genes_unique)

cat("目标基因总数：", length(target_gene_vec), "\n")
cat("成功匹配到的基因数：", length(matched_genes), "\n")
cat("未匹配到的基因数：", length(unmatched_genes), "\n\n")

cat("成功匹配到的基因：\n")
print(sort(matched_genes))

cat("\n未匹配到的基因：\n")
print(sort(unmatched_genes))

fwrite(
  data.table(matched_gene = sort(matched_genes)),
  file.path(workdir, "eqtlgen_matched_genes.tsv"),
  sep = "\t"
)

fwrite(
  data.table(unmatched_gene = sort(unmatched_genes)),
  file.path(workdir, "eqtlgen_unmatched_genes.tsv"),
  sep = "\t"
)

#========================================================
# 12. 提取目标基因对应的 eQTL 记录
#========================================================

message("\nStep 4: Extract eQTL records for matched DE-FRGs")

eqtl_subset <- copy(
  eqtl_all[GeneSymbol_clean %in% target_gene_vec]
)

eqtl_subset[, gene := GeneSymbol_clean]

fwrite(
  eqtl_subset,
  file.path(workdir, "eqtlgen_target_gene_subset.tsv"),
  sep = "\t"
)

cat("提取到的 eQTL 记录数：", nrow(eqtl_subset), "\n")
cat("有 eQTL 记录的目标基因数：", length(unique(eqtl_subset$gene)), "\n")

if (nrow(eqtl_subset) == 0) {
  stop("没有从 eQTLGen 中匹配到任何目标基因，请检查 gene symbol 是否一致。")
}

#========================================================
# 13. 统计每个目标基因的 eQTL 信息
#========================================================

message("\nStep 5: Count eQTL records for each gene")

if (!"Pvalue" %in% names(eqtl_subset)) {
  stop("eQTLGen 文件中没有 Pvalue 列，请检查文件格式。")
}

eqtl_counts <- eqtl_subset[, .(
  n_eqtl_rows = .N,
  n_genomewide_sig = sum(Pvalue < 5e-8, na.rm = TRUE),
  best_eqtl_p = suppressWarnings(min(Pvalue, na.rm = TRUE))
), by = gene]

eqtl_counts[!is.finite(best_eqtl_p), best_eqtl_p := NA_real_]

setorder(eqtl_counts, -n_genomewide_sig, best_eqtl_p, gene)

fwrite(
  eqtl_counts,
  file.path(workdir, "eqtlgen_target_gene_counts.tsv"),
  sep = "\t"
)

cat("\n每个目标基因的 eQTL 统计：\n")
print(eqtl_counts)

#========================================================
# 14. 提取目标基因对应的 SNP 列表
#========================================================

message("\nStep 6: Extract SNP list for target eQTL records")

if ("SNPName" %in% names(eqtl_subset)) {
  snp_col <- "SNPName"
} else if ("SNP" %in% names(eqtl_subset)) {
  snp_col <- "SNP"
} else if ("rsid" %in% names(eqtl_subset)) {
  snp_col <- "rsid"
} else {
  snp_col <- grep("snp|rsid", names(eqtl_subset), ignore.case = TRUE, value = TRUE)[1]
}

if (is.na(snp_col) || length(snp_col) == 0) {
  stop("eQTLGen 文件中无法识别 SNP 列。")
}

cat("eQTLGen 使用的 SNP 列为：", snp_col, "\n")

target_snps <- sort(unique(eqtl_subset[[snp_col]]))

fwrite(
  data.table(snp = target_snps),
  file.path(workdir, "eqtlgen_target_snps.tsv"),
  sep = "\t"
)

cat("目标 SNP 数量：", length(target_snps), "\n")
cat("前 20 个 SNP：\n")
print(head(target_snps, 20))

#========================================================
# 15. Step 7：提取这些 SNP 的等位基因频率信息
#========================================================

message("\nStep 7: Extract allele-frequency rows for the target SNPs")

af_header <- fread(af_file, nrows = 5)

cat("\nAllele frequency 文件列名：\n")
print(names(af_header))

af_all <- fread(af_file)

af_snp_col <- grep("^snp$|snp|rsid", names(af_all), ignore.case = TRUE, value = TRUE)[1]

if (is.na(af_snp_col) || length(af_snp_col) == 0) {
  stop("Allele frequency 文件中无法识别 SNP 列。")
}

cat("Allele frequency 文件使用的 SNP 列为：", af_snp_col, "\n")

af_subset <- copy(
  af_all[get(af_snp_col) %in% target_snps]
)

fwrite(
  af_subset,
  file.path(workdir, "eqtlgen_target_snps_allele_frequency.tsv"),
  sep = "\t"
)

cat("提取到 allele frequency 记录数：", nrow(af_subset), "\n")

cat("\nallele frequency 前几行：\n")
print(head(af_subset))

#========================================================
# 16. 生成目标基因追踪表
#========================================================

message("\nStep 8: Build gene trace table")

trace_dt <- unique(target_genes[, .(
  gene,
  gene_clean,
  source,
  selection_stage,
  logFC,
  adj.P.Val,
  expression_trend
)])

trace_dt[, in_start_set := TRUE]

trace_dt <- merge(
  trace_dt,
  eqtl_counts,
  by.x = "gene_clean",
  by.y = "gene",
  all.x = TRUE
)

trace_dt[, has_eqtl_record := !is.na(n_eqtl_rows) & n_eqtl_rows > 0]
trace_dt[, has_genomewide_eqtl := !is.na(n_genomewide_sig) & n_genomewide_sig > 0]

trace_dt[is.na(n_eqtl_rows), n_eqtl_rows := 0]
trace_dt[is.na(n_genomewide_sig), n_genomewide_sig := 0]
trace_dt[is.na(best_eqtl_p), best_eqtl_p := NA_real_]

setorder(trace_dt, -has_genomewide_eqtl, -n_genomewide_sig, -has_eqtl_record, gene_clean)

fwrite(
  trace_dt,
  file.path(workdir, "eqtlgen_gene_trace.tsv"),
  sep = "\t"
)

cat("\n基因追踪表：\n")
print(trace_dt)

#========================================================
# 17. 生成总报告
#========================================================

summary_report <- data.table(
  item = c(
    "target_gene_file",
    "eqtl_file",
    "allele_frequency_file",
    "n_target_genes",
    "n_eqtlgen_unique_genes",
    "n_matched_genes",
    "n_unmatched_genes",
    "n_eqtl_subset_rows",
    "n_eqtl_subset_genes",
    "n_target_snps",
    "n_allele_frequency_rows"
  ),
  value = c(
    target_gene_file,
    eqtl_file,
    af_file,
    as.character(length(target_gene_vec)),
    as.character(length(eqtl_genes_unique)),
    as.character(length(matched_genes)),
    as.character(length(unmatched_genes)),
    as.character(nrow(eqtl_subset)),
    as.character(length(unique(eqtl_subset$gene))),
    as.character(length(target_snps)),
    as.character(nrow(af_subset))
  )
)

fwrite(
  summary_report,
  file.path(workdir, "eqtlgen_summary_report.tsv"),
  sep = "\t"
)

cat("\n总报告：\n")
print(summary_report)

#========================================================
# 18. 完成
#========================================================

message("\nFinished.")
message("输出文件包括：")
message("1. standardized_target_genes.tsv")
message("2. eqtlgen_matched_genes.tsv")
message("3. eqtlgen_unmatched_genes.tsv")
message("4. eqtlgen_target_gene_subset.tsv")
message("5. eqtlgen_target_gene_counts.tsv")
message("6. eqtlgen_target_snps.tsv")
message("7. eqtlgen_target_snps_allele_frequency.tsv")
message("8. eqtlgen_gene_trace.tsv")
message("9. eqtlgen_summary_report.tsv")

#!/usr/bin/env Rscript
## 4.prepare_mr_subsets_local.R
# 作用：按文章的 IV 筛选逻辑，在本地完成 P < 5e-8、LD clumping、F > 10、GWAS Catalog 混杂剔除、
#CRC outcome 提取，得到可进入 MR 的工具变量集合。
#========================================================
# 0. 基础设置
#========================================================

##系统报错改为英文
Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = FALSE)

##清空环境
rm(list=ls())
# 设置超时时间，避免下载或读取大文件时超时
options(timeout = max(3600, getOption("timeout")))


#========================================================
# 1. 检查并安装所需 R 包
#========================================================

# 本脚本依赖的包：
# data.table      -> 高效读取和处理表格
# jsonlite        -> JSON 相关处理
# Rsamtools       -> 读取 Tabix 索引的 VCF 文件
# GenomicRanges   -> 构建基因组区间对象
required_packages <- c("data.table", "jsonlite", "Rsamtools", "GenomicRanges")
# 找出还没安装的包
missing <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
# CRAN 包直接 install.packages
# Bioconductor 包（Rsamtools、GenomicRanges）这里默认你已经装好了
if (length(missing)) {
  install.packages(
    setdiff(missing, c("Rsamtools", "GenomicRanges")),
    repos = "https://cloud.r-project.org"
  )
}
# 加载所需包
library(data.table)
library(jsonlite)
library(Rsamtools)
library(GenomicRanges)

#========================================================
# 2. 定义项目路径和输入文件路径
#========================================================
# 设置当前工作目录
workdir <- "C:/Users/92066/Desktop/单细胞孟德尔复现/MR"
setwd(workdir)

###  上节课生成的数据
# 把 intersect_genes.csv 读进来以后，经过脚本标准化整理后生成的结果表。
standardized_target_genes_path <- file.path(workdir, "standardized_target_genes.tsv")

# 最开始拿来做后续筛选的目标基因列表
deg_gene_path <- file.path(workdir, "intersect_genes.csv")

# 从整个 eQTLGen 大数据库里，筛出来的与你目标基因有关的所有 eQTL 记录子集
eqtlgen_target_gene_subset <- file.path(workdir, "eqtlgen_target_gene_subset.tsv")

# eQTLGen 的等位基因频率文件里提取出来的目标 SNP 频率信息表
eqtlgen_target_snps_allele_frequency <- file.path(workdir, "eqtlgen_target_snps_allele_frequency.tsv")

# 这节课需要下载的数据
# LD 参考面板前缀（EUR 欧洲人群）
ld_ref_prefix <- file.path(workdir, "ldref", "EUR", "1000G.EUR.QC")

# 结局 VCF 文件路径（这里是 CRC，结直肠癌）
outcome_vcf <- file.path(workdir, "openGWAS", "ebi-a-GCST90018808.vcf.gz")

# 本地 GWAS Catalog 关联结果表
gwas_catalog_assoc_path <- file.path(workdir, "GWASCatalog", "gwas-catalog-associations_ontology-annotated.tsv")

# PLINK 可执行文件路径
# 如果系统环境变量里设置了 PLINK_BIN，就优先使用环境变量
# 否则使用项目默认路径
plink_bin <- file.path(workdir, "plink", "plink.exe")

#========================================================
# 3. 定义辅助函数：自动识别“基因列”
#========================================================
# 不同文件里基因列名称可能不一样：
# 可能叫 gene / symbol / GeneSymbol / exposure 等
# 这个函数会自动找到最可能的基因列
resolve_gene_col <- function(dt) {
  nms <- names(dt)
  # 优先匹配常见列名
  hits <- c("gene", "gene_symbol", "symbol", "Gene", "GeneSymbol", "exposure")
  hit <- hits[hits %in% nms]
  # 如果找到了，返回第一个
  if (length(hit)) {
    return(hit[1])
  }
  # 如果没找到精确列名，则模糊搜索
  idx <- grep("gene|symbol|exposure", nms, ignore.case = TRUE)
  # 如果还是找不到，就报错
  if (!length(idx)) {
    stop("Could not identify a gene column.")
  }
  nms[idx[1]]
}

#========================================================
# 3.1 定义辅助函数：清洗基因名
#========================================================
clean_gene <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("\ufeff", "", x, fixed = TRUE)
  x <- gsub("\u00A0", "", x, fixed = TRUE)
  x <- gsub("\"", "", x, fixed = TRUE)
  x <- gsub("'", "", x, fixed = TRUE)
  x <- gsub("\\s+", "", x)
  x <- toupper(x)
  x
}

#========================================================
# 4.1 定义辅助函数：读取起始基因集合
#========================================================
read_target_genes <- function(standardized_path, deg_path) {
  if (file.exists(standardized_path)) {
    dt <- fread(standardized_path)
    gene_col <- resolve_gene_col(dt)
    setnames(dt, gene_col, "gene")
  } else if (file.exists(deg_path)) {
    dt <- fread(deg_path)
    if (ncol(dt) == 1) {
      setnames(dt, names(dt)[1], "gene")
    } else {
      gene_col <- resolve_gene_col(dt)
      setnames(dt, gene_col, "gene")
    }
    dt[, source := "intersect_genes"]
    dt[, selection_stage := "start_intersect_genes"]
  } else {
    stop(sprintf(
      "Neither standardized target genes nor intersect_genes.csv were found: %s ; %s",
      normalizePath(standardized_path, winslash = "/", mustWork = FALSE),
      normalizePath(deg_path, winslash = "/", mustWork = FALSE)
    ))
  }
  dt[, gene := as.character(gene)]
  dt[, gene := clean_gene(gene)]
  dt <- dt[
    !is.na(gene) &
      nzchar(gene) &
      !tolower(gene) %in% c("gene", "symbol", "genesymbol", "gene_symbol", "hgnc")
  ]
  if (!"source" %in% names(dt)) {
    dt[, source := "intersect_genes"]
  }
  if (!"selection_stage" %in% names(dt)) {
    dt[, selection_stage := "start_intersect_genes"]
  }
  unique(dt[!is.na(gene) & nzchar(gene)], by = "gene")
}

#========================================================
# 4. 定义辅助函数：解析 PLINK clump 结果
#========================================================
# PLINK clump 后会生成 .clumped 文件
# 这个函数的作用是把 .clumped 文件整理成 data.table
parse_clumped_file <- function(path) {
  # 读取全部行
  lines <- readLines(path, warn = FALSE)
  # 去掉空行
  lines <- lines[nzchar(trimws(lines))]
  if (!length(lines)) {
    return(data.table())
  }
  # 去掉 PLINK 开头的说明文字，只保留数据部分
  body <- lines[!grepl("^\\s*PLINK", lines) & !grepl("^\\s*$", lines)]
  if (!length(body)) {
    return(data.table())
  }
  # 拼成文本，再用 fread 读入
  text <- paste(body, collapse = "\n")
  dt <- fread(text = text, fill = TRUE)
  # 去掉列名前后的空格
  setnames(dt, trimws(names(dt)))
  dt
}

#========================================================
# 5. 定义辅助函数：计算 eaf（效应等位基因频率）
#========================================================
# assessed  = 效应等位基因
# other     = 非效应等位基因
# allele_a / allele_b = 频率文件中的两个等位基因
# af_b      = allele_b 的频率
#
# 目的：把 allele frequency 转换成“效应等位基因频率”
calc_eaf <- function(assessed, other, allele_a, allele_b, af_b) {
  out <- rep(NA_real_, length(assessed))
  # 如果效应等位基因正好是 allele_b
  idx_b <- assessed == allele_b & other == allele_a
  # 如果效应等位基因正好是 allele_a
  idx_a <- assessed == allele_a & other == allele_b
  out[idx_b] <- af_b[idx_b]
  out[idx_a] <- 1 - af_b[idx_a]
  out
}

#========================================================
# 6. 定义辅助函数：由 Z 值换算 beta 和 se
#========================================================
# 利用 Z 分数、等位基因频率和样本量
# 反推出 beta 和标准误 se
z_to_beta <- function(z, eaf, n) {
  z / sqrt(2 * eaf * (1 - eaf) * (n + z^2))
}
z_to_se <- function(z, eaf, n) {
  1 / sqrt(2 * eaf * (1 - eaf) * (n + z^2))
}

#========================================================
# 7. 定义辅助函数：解析 VCF 中单条记录
#========================================================
# 这里主要是从 outcome 的 VCF 文件中提取：
# chr、pos、rsid、效应等位基因、beta、se、p 值等信息
parse_vcf_record <- function(line) {
  # 按制表符拆分 VCF 一行
  fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
  # FORMAT 列定义的字段名
  fmt_keys <- strsplit(fields[9], ":", fixed = TRUE)[[1]]
  # 样本列对应的字段值
  fmt_vals <- strsplit(fields[10], ":", fixed = TRUE)[[1]]
  # 转成命名列表，方便后面提取 ES、SE、LP
  fmt <- as.list(fmt_vals)
  names(fmt) <- fmt_keys
  # LP 一般表示 -log10(P)
  lp <- suppressWarnings(as.numeric(fmt$LP))
  # INFO 列中提取 AF
  af_match <- regmatches(fields[8], regexpr("AF=[^;]+", fields[8]))
  af_value <- NA_real_
  if (length(af_match) && nzchar(af_match[1])) {
    af_value <- suppressWarnings(as.numeric(sub("AF=", "", af_match[1])))
  }
  # 由 LP 反算 p 值
  pval <- ifelse(is.na(lp), NA_real_, 10^(-lp))
  # 返回整理后的结果
  data.table(
    id = "ebi-a-GCST90018808.vcf",
    trait = "Colorectal cancer",
    chr = suppressWarnings(as.integer(fields[1])),
    position = suppressWarnings(as.integer(fields[2])),
    rsid = fields[3],
    ea = toupper(fields[5]),   # ALT 作为 effect allele
    nea = toupper(fields[4]),  # REF 作为 non-effect allele
    eaf = af_value,
    beta = suppressWarnings(as.numeric(fmt$ES)),
    se = suppressWarnings(as.numeric(fmt$SE)),
    p = pval,
    n = NA_real_
  )
}

#========================================================
# 8. 定义辅助函数：从本地 VCF 中批量提取 outcome 信息
#========================================================
# 输入：一批 SNP（chr、pos、rsid）
# 输出：这些 SNP 在 outcome VCF 中对应的效应信息
read_outcome_from_vcf <- function(snp_dt, vcf_path) {
  # 打开 Tabix 索引文件
  tbx <- TabixFile(vcf_path)
  on.exit(close(tbx))
  open(tbx)
  # 整理唯一 SNP 列表
  regions <- unique(snp_dt[, .(
    chr = as.character(chr.exposure),
    pos = as.integer(pos.exposure),
    rsid = SNP
  )])
  # 构建基因组区间对象
  gr <- GRanges(
    seqnames = regions$chr,
    ranges = IRanges(start = regions$pos, end = regions$pos)
  )
  # 读取这些位置的 VCF 记录
  raw_list <- scanTabix(tbx, param = gr)
  out <- vector("list", length(raw_list))
  for (i in seq_along(raw_list)) {
    rsid_target <- regions$rsid[i]
    lines <- raw_list[[i]]
    if (!length(lines)) {
      next
    }
    # 解析每一条 VCF 行
    parsed <- rbindlist(lapply(lines, parse_vcf_record), fill = TRUE)
    # 只保留目标 rsid
    parsed <- parsed[rsid == rsid_target]
    if (nrow(parsed)) {
      out[[i]] <- parsed[1]
    }
  }
  rbindlist(out, fill = TRUE)
}

#========================================================
# 9. 定义辅助函数：从 GWAS Catalog 中提取 trait 文本
#========================================================
# 不同版本的 GWAS Catalog 里，trait 相关列名可能不同
# 这个函数统一提取最合适的 trait 描述文字
extract_trait_text <- function(dt) {
  if (!nrow(dt)) {
    return(rep("", 0))
  }
  # 查找所有可能与 trait 有关的列
  candidate_cols <- setdiff(
    grep("trait|mappedTrait|efoTrait", names(dt), ignore.case = TRUE, value = TRUE),
    grep("_links\\.|href$", names(dt), ignore.case = TRUE, value = TRUE)
  )
  # 去掉更新时间、映射说明等无关列
  candidate_cols <- candidate_cols[
    !grepl("last|update|mapping", candidate_cols, ignore.case = TRUE)
  ]
  if (!length(candidate_cols)) {
    return(rep("", nrow(dt)))
  }
  out <- rep("", nrow(dt))
  for (col in candidate_cols) {
    vals <- dt[[col]]
    # 如果是 list，就展开
    if (is.list(vals)) {
      vals <- vapply(vals, function(x) paste(as.character(unlist(x)), collapse = " | "), character(1))
    } else {
      vals <- as.character(vals)
    }
    vals[is.na(vals)] <- ""
    # 优先填充第一个非空值
    better <- nzchar(trimws(vals)) & !nzchar(trimws(out))
    out[better] <- vals[better]
  }
  out
}

#========================================================
# 10. 定义辅助函数：拼接多个列的文本
#========================================================
# 用于把多个 trait/SNP 相关列拼起来，方便后面统一搜索关键词
paste_nonmissing_cols <- function(dt, cols) {
  if (!length(cols) || !nrow(dt)) {
    return(rep("", nrow(dt)))
  }
  vals <- lapply(cols, function(col) {
    x <- as.character(dt[[col]])
    x[is.na(x)] <- ""
    x
  })
  trimws(do.call(paste, c(vals, sep = " | ")))
}

#========================================================
# 11. 定义辅助函数：从本地 GWAS Catalog 中查找混杂 SNP
#========================================================
# 作用：
# 检查工具变量 SNP 是否同时和教育程度、认知能力、脑容量等混杂因素相关
# 如果显著相关，则后面需要剔除
load_local_gwas_catalog_matches <- function(assoc_path, query_snps, keyword_pattern) {
  # 文件不存在就报错
  if (!file.exists(assoc_path)) {
    stop(sprintf(
      "Local GWAS Catalog associations file is missing: %s. Run 2.download_mr_data.R first.",
      normalizePath(assoc_path, winslash = "/", mustWork = FALSE)
    ))
  }
  # 先只读取表头，看看有哪些列
  header <- names(fread(assoc_path, nrows = 0, showProgress = FALSE))
  if (!length(header)) {
    stop(sprintf("Could not read the GWAS Catalog header from %s", assoc_path))
  }
  # 只保留我们关心的列
  wanted_orig <- unique(c(
    "SNPS",
    "STRONGEST SNP-RISK ALLELE",
    "DISEASE/TRAIT",
    "MAPPED_TRAIT",
    "MAPPED_TRAIT_URI",
    "P-VALUE",
    "STUDY ACCESSION",
    "PUBMEDID"
  ))
  select_orig <- intersect(wanted_orig, header)
  # 至少要有 SNP 相关列
  if (!length(intersect(c("SNPS", "STRONGEST SNP-RISK ALLELE"), select_orig))) {
    stop("The local GWAS Catalog file does not contain a usable SNP column.")
  }
  message("  Reading the local GWAS Catalog bulk associations file")
  # 读取所需列
  dt <- fread(assoc_path, select = select_orig, showProgress = FALSE)
  # 为了防止列名中有特殊字符，这里统一合法化
  orig_names <- names(dt)
  clean_names <- make.names(orig_names, unique = TRUE)
  setnames(dt, clean_names)
  name_map <- setNames(clean_names, orig_names)
  # 找到 SNP、trait、p 值、研究编号等列
  snp_cols <- unname(name_map[intersect(c("SNPS", "STRONGEST SNP-RISK ALLELE"), names(name_map))])
  trait_cols <- unname(name_map[intersect(c("DISEASE/TRAIT", "MAPPED_TRAIT", "MAPPED_TRAIT_URI"), names(name_map))])
  p_col <- unname(name_map[intersect("P-VALUE", names(name_map))])[1]
  study_col <- unname(name_map[intersect("STUDY ACCESSION", names(name_map))])[1]
  pubmed_col <- unname(name_map[intersect("PUBMEDID", names(name_map))])[1]
  # 查询 SNP 建立映射表
  query_map <- unique(data.table(
    query_snp = tolower(as.character(query_snps)),
    SNP = as.character(query_snps)
  ))
  query_map <- query_map[!is.na(query_snp) & nzchar(query_snp)]
  if (!nrow(query_map)) {
    return(list(matches = data.table(), audit = data.table()))
  }
  # 把 SNP 列拼接后统一小写
  snp_text <- tolower(paste_nonmissing_cols(dt, snp_cols))
  # 从文本中提取 rs 编号
  rs_matches <- regmatches(snp_text, gregexpr("rs[0-9]+", snp_text, perl = TRUE))
  hit_idx <- which(lengths(rs_matches) > 0L)
  # 转成长表：一行一个 rs
  match_long <- if (length(hit_idx)) {
    data.table(
      row_id = rep.int(hit_idx, lengths(rs_matches[hit_idx])),
      query_snp = unlist(rs_matches[hit_idx], use.names = FALSE)
    )
  } else {
    data.table(row_id = integer(), query_snp = character())
  }
  # 只保留我们要查的 SNP
  match_long <- unique(match_long[query_snp %in% query_map$query_snp], by = c("row_id", "query_snp"))
  # 没匹配上则返回审计信息
  if (!nrow(match_long)) {
    audit <- copy(query_map)
    audit[, `:=`(
      source_mode = "local_bulk_catalog",
      match_status = "no_local_match",
      n_associations = 0L,
      n_remove = 0L
    )]
    return(list(matches = data.table(), audit = audit))
  }
  # 加入行号并合并回原始表
  dt[, row_id := .I]
  matched_rows <- merge(match_long, dt, by = "row_id", all.x = TRUE, allow.cartesian = TRUE)
  # 拼接 trait 文本并转小写
  matched_rows[, trait_text := tolower(paste_nonmissing_cols(.SD, trait_cols))]
  # 处理 p 值
  if (!is.na(p_col) && nzchar(p_col)) {
    matched_rows[, pvalue_num := suppressWarnings(as.numeric(gsub(",", "", get(p_col))))]
  } else {
    matched_rows[, pvalue_num := NA_real_]
  }
  # 研究编号
  if (!is.na(study_col) && nzchar(study_col)) {
    matched_rows[, study_accession := as.character(get(study_col))]
  } else {
    matched_rows[, study_accession := NA_character_]
  }
  # PubMed ID
  if (!is.na(pubmed_col) && nzchar(pubmed_col)) {
    matched_rows[, pubmedid := as.character(get(pubmed_col))]
  } else {
    matched_rows[, pubmedid := NA_character_]
  }
  # 有没有 trait 文本
  matched_rows[, trait_text_available := nzchar(trimws(trait_text))]
  # 是否命中混杂关键词
  matched_rows[, confounder_keyword_match := trait_text_available & grepl(keyword_pattern, trait_text)]
  # 判定是否剔除
  # 条件：和混杂因素相关，并且 p < 1e-5
  matched_rows[, confounder_filter_status := fifelse(
    confounder_keyword_match & !is.na(pvalue_num) & pvalue_num < 1e-5,
    "remove",
    fifelse(trait_text_available, "keep_checked", "keep_unresolved_no_trait_text")
  )]
  matched_rows[, source_mode := "local_bulk_catalog"]
  # 生成审计表：每个 SNP 匹配了多少关联，其中多少需要剔除
  audit <- matched_rows[, .(
    n_associations = .N,
    n_remove = sum(confounder_filter_status == "remove", na.rm = TRUE)
  ), by = query_snp]
  audit <- merge(query_map, audit, by = "query_snp", all.x = TRUE)
  audit[is.na(n_associations), n_associations := 0L]
  audit[is.na(n_remove), n_remove := 0L]
  audit[, source_mode := "local_bulk_catalog"]
  audit[, match_status := fifelse(n_associations > 0, "matched_in_local_bulk", "no_local_match")]
  list(matches = matched_rows, audit = audit)
}

#========================================================
# 12. 定义辅助函数：统计每一步保留下来的 SNP 数量
#========================================================
# 统计整体数量
overall_count <- function(dt, step_label) {
  if (!nrow(dt)) {
    return(data.table(
      step = step_label,
      n_gene_snp = 0L,
      n_unique_snps = 0L,
      n_genes = 0L
    ))
  }
  data.table(
    step = step_label,
    n_gene_snp = as.integer(uniqueN(dt[, .(gene, SNP)])),
    n_unique_snps = as.integer(uniqueN(dt$SNP)),
    n_genes = as.integer(uniqueN(dt$gene))
  )
}
# 按基因统计 SNP 数量
by_gene_count <- function(dt) {
  if (!nrow(dt)) {
    return(data.table(gene = character(), n_snps = integer()))
  }
  dt[, .(n_snps = as.integer(uniqueN(SNP))), by = gene]
}
# 把每一步按基因统计的结果合并到总表里
merge_trace_col <- function(base, counts_dt, col_name) {
  tmp <- copy(counts_dt)
  if (!nrow(tmp)) {
    base[[col_name]] <- 0L
    return(base)
  }
  setnames(tmp, "n_snps", col_name)
  out <- merge(base, tmp, by = "gene", all.x = TRUE)
  out[is.na(get(col_name)), (col_name) := 0L]
  out
}

#========================================================
# 13. 检查关键外部文件是否存在
#========================================================
# 检查 PLINK 是否存在
if (!file.exists(plink_bin)) {
  stop(sprintf(
    "PLINK binary not found: %s",
    normalizePath(plink_bin, winslash = "/", mustWork = FALSE)
  ))
}
# 检查 LD 参考面板是否完整（.bed .bim .fam）
if (!all(file.exists(paste0(ld_ref_prefix, c(".bed", ".bim", ".fam"))))) {
  stop(sprintf(
    "LD reference panel prefix is incomplete: %s",
    normalizePath(ld_ref_prefix, winslash = "/", mustWork = FALSE)
  ))
}

#========================================================
# 14. 读取目标基因和 eQTL 数据
#========================================================
setwd(workdir)
# 读取目标基因列表
target_genes <- read_target_genes(standardized_target_genes_path, deg_gene_path)
# 读取 eQTL 数据
eqtl <- fread( eqtlgen_target_gene_subset )
# 读取等位基因频率文件
af <- fread( eqtlgen_target_snps_allele_frequency)

#========================================================
# 15. 检查 eQTL 数据是否包含关键列
#========================================================
required_eqtl_cols <- c(
  "gene", "SNP", "SNPChr", "SNPPos", "AssessedAllele",
  "OtherAllele", "Zscore", "NrSamples", "Pvalue"
)
missing_eqtl_cols <- setdiff(required_eqtl_cols, names(eqtl))
if (length(missing_eqtl_cols)) {
  stop(sprintf(
    "Missing expected eQTLGen columns: %s",
    paste(missing_eqtl_cols, collapse = ", ")
  ))
}

#========================================================
# 16. 标准化 eQTL 数据格式
#========================================================
# 这里只保留后续分析需要的列，并统一类型
eqtl <- unique(eqtl[, .(
  gene = clean_gene(gene),
  SNP = as.character(SNP),
  SNPChr = suppressWarnings(as.integer(SNPChr)),
  SNPPos = suppressWarnings(as.integer(SNPPos)),
  AssessedAllele = toupper(as.character(AssessedAllele)),
  OtherAllele = toupper(as.character(OtherAllele)),
  Zscore = suppressWarnings(as.numeric(Zscore)),
  NrSamples = suppressWarnings(as.numeric(NrSamples)),
  Pvalue = suppressWarnings(as.numeric(Pvalue))
)], by = c("gene", "SNP"))
# 只保留目标基因对应的 eQTL
eqtl <- eqtl[gene %in% target_genes$gene]

cat("目标基因数：", nrow(target_genes), "\n")
cat("保留目标基因 eQTL 记录数：", nrow(eqtl), "\n")
cat("保留目标基因数：", uniqueN(eqtl$gene), "\n")

if (nrow(eqtl) == 0) {
  stop("目标基因与 eQTL 子集没有匹配，请检查 gene symbol 或 standardized_target_genes.tsv。")
}

#========================================================
# 17. 初始化“筛选轨迹”记录对象
#========================================================
# trace_by_gene   -> 记录每个基因在每一步剩余多少 SNP
# trace_overall   -> 记录整体在每一步还剩多少 SNP
trace_by_gene <- list()
trace_overall <- list()
#========================================================
# 18. Step 1：保留全基因组显著 eQTL（P < 5e-8）
#========================================================
message("Step 1: Filter eQTLGen rows to genome-wide significant candidate IVs")
# 记录初始 eQTL 数量
trace_by_gene$all_eqtl_rows <- by_gene_count(eqtl)
trace_overall[[length(trace_overall) + 1L]] <- overall_count(eqtl, "all_eqtl_rows")
# 只保留显著 eQTL
sig_eqtl <- eqtl[!is.na(Pvalue) & Pvalue < 5e-8]
# 记录显著 SNP 数量
trace_by_gene$genomewide_sig_snps <- by_gene_count(sig_eqtl)
trace_overall[[length(trace_overall) + 1L]] <- overall_count(sig_eqtl, "genomewide_sig_snps")
#========================================================
# 19. Step 2：按基因做 LD clump，筛选独立工具变量
#========================================================
message("Step 2: LD clump within each gene using the local EUR reference")
# 建立临时目录，保存 clump 输入输出文件
tmp_dir <- file.path(tempdir(), "mr_local_clump")
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
clumped_list <- list()
clump_audit <- list()
# 对每个基因分别做 clump
for (gene_name in unique(sig_eqtl$gene)) {
  gene_eqtl <- sig_eqtl[gene == gene_name]
  # 去重后保留 SNP、P 值、染色体、位置
  dat <- unique(gene_eqtl[, .(SNP, Pvalue, SNPChr, SNPPos)])
  if (!nrow(dat)) {
    next
  }
  # 输入文件和输出前缀
  in_file <- file.path(tmp_dir, paste0(gene_name, "_clump_input.tsv"))
  out_prefix <- file.path(tmp_dir, paste0(gene_name, "_clump"))
  # 写出给 PLINK 的输入表
  fwrite(dat[, .(SNP, P = Pvalue)], in_file, sep = "\t")
  # PLINK clump 参数：
  # --clump-p1 5e-8     显著性阈值
  # --clump-r2 0.001    LD 很严格
  # --clump-kb 100      100 kb 窗口
  args <- c(
    "--bfile", ld_ref_prefix,
    "--clump", in_file,
    "--clump-snp-field", "SNP",
    "--clump-field", "P",
    "--clump-p1", "5e-8",
    "--clump-r2", "0.001",
    "--clump-kb", "100",
    "--out", out_prefix
  )
  # 运行 PLINK
  status <- system2(plink_bin, args = args, stdout = TRUE, stderr = TRUE)
  exit_code <- attr(status, "status")
  if (is.null(exit_code)) {
    exit_code <- 0L
  }
  # 保存运行日志
  writeLines(status, paste0(out_prefix, ".stdout.txt"))
  # 运行失败则报错
  if (exit_code != 0L) {
    stop(sprintf(
      "PLINK clumping failed for %s. See %s",
      gene_name, paste0(out_prefix, ".stdout.txt")
    ))
  }
  clumped_path <- paste0(out_prefix, ".clumped")
  # 如果没生成 clumped 文件，说明该基因没有保留 SNP
  if (!file.exists(clumped_path)) {
    clump_audit[[gene_name]] <- data.table(
      gene = gene_name,
      total_input = nrow(dat),
      retained = 0L,
      note = "no_clumped_file"
    )
    next
  }
  # 读取 clump 结果
  lead_dt <- parse_clumped_file(clumped_path)
  # 如果没有 lead SNP，则记录
  if (!nrow(lead_dt) || !"SNP" %in% names(lead_dt)) {
    clump_audit[[gene_name]] <- data.table(
      gene = gene_name,
      total_input = nrow(dat),
      retained = 0L,
      note = "no_lead_snps"
    )
    next
  }
  # 保留 lead SNP 对应信息
  keep_dt <- merge(
    data.table(SNP = trimws(as.character(lead_dt$SNP))),
    dat,
    by = "SNP",
    all.x = TRUE
  )
  keep_dt[, gene := gene_name]
  clumped_list[[gene_name]] <- keep_dt[]
  # 记录 clump 审计结果
  clump_audit[[gene_name]] <- data.table(
    gene = gene_name,
    total_input = nrow(dat),
    retained = nrow(keep_dt),
    note = "ok"
  )
}
# 合并所有基因的 clump 结果
clumped_dt <- rbindlist(clumped_list, fill = TRUE)
clump_audit_dt <- rbindlist(clump_audit, fill = TRUE)
# 保存 clump 审计表
fwrite(clump_audit_dt,  "local_plink_clump_audit.tsv" , sep = "\t")
# 如果一个 SNP 都没留下，就终止
if (!nrow(clumped_dt)) {
  stop("No SNPs remained after LD clumping.")
}
# 记录 clump 后的 SNP 数量
trace_by_gene$ld_clumped_snps <- by_gene_count(clumped_dt[, .(gene, SNP)])
trace_overall[[length(trace_overall) + 1L]] <- overall_count(clumped_dt[, .(gene, SNP)], "ld_clumped_snps")
#========================================================
# 20. 根据 Z 值和等位基因频率，补充 exposure 效应值
#========================================================
# 合并显著 eQTL 的等位基因、Z 值、样本量信息
clumped_exp <- merge(
  clumped_dt,
  sig_eqtl[, .(
    gene,
    SNP,
    AssessedAllele,
    OtherAllele,
    Zscore,
    NrSamples
  )],
  by = c("gene", "SNP"),
  all.x = TRUE
)
# 再合并等位基因频率信息
clumped_exp <- merge(
  clumped_exp,
  af[, .(
    SNP = as.character(SNP),
    AlleleA = toupper(as.character(AlleleA)),
    AlleleB = toupper(as.character(AlleleB)),
    eaf_b = suppressWarnings(as.numeric(AlleleB_all))
  )],
  by = "SNP",
  all.x = TRUE
)
# 计算效应等位基因频率
clumped_exp[, eaf.exposure := calc_eaf(
  AssessedAllele, OtherAllele,
  AlleleA, AlleleB, eaf_b
)]
# 用 Z 值反推 beta 和 se
clumped_exp[, beta.exposure := z_to_beta(Zscore, eaf.exposure, NrSamples)]
clumped_exp[, se.exposure := z_to_se(Zscore, eaf.exposure, NrSamples)]
# 计算 F 统计量，用于判断工具变量强弱
clumped_exp[, F_stat := (beta.exposure / se.exposure)^2]
# 先预留 exposure ID
clumped_exp[, id.exposure := NA_character_]
# 标记工具变量来源
clumped_exp[, instrument_source := "local_eqtlgen_clumped"]
#========================================================
# 21. 生成 exposure 工具变量表
#========================================================
clumped_instruments <- unique(clumped_exp[, .(
  gene,
  id.exposure,
  SNP,
  chr.exposure = SNPChr,
  pos.exposure = SNPPos,
  beta.exposure,
  se.exposure,
  pval.exposure = Pvalue,
  samplesize.exposure = NrSamples,
  effect_allele.exposure = AssessedAllele,
  other_allele.exposure = OtherAllele,
  eaf.exposure,
  F_stat,
  instrument_source
)], by = c("gene", "SNP"))
# 保存 clump 后工具变量
fwrite(
  clumped_instruments,
  "eqtlgen_target_gene_clumped_instruments.tsv" ,
  sep = "\t"
)
#========================================================
# 22. Step 3：去掉弱工具变量（F <= 10）
#========================================================
message("Step 3: Remove weak instruments with F <= 10")
# 经典 MR 中常用 F > 10 作为强工具变量标准
strong_instruments <- clumped_instruments[!is.na(F_stat) & F_stat > 10]
# 记录该步骤保留数量
trace_by_gene$strong_f_snps <- by_gene_count(strong_instruments)
trace_overall[[length(trace_overall) + 1L]] <- overall_count(strong_instruments, "strong_f_snps")

#========================================================
# 23. Step 4：去掉可能和结直肠癌混杂因素相关的 SNP
#========================================================

message("Step 4: Remove SNPs linked to colorectal cancer-related confounders using the local bulk GWAS Catalog release")

# 定义需要关注的 CRC 相关混杂关键词
# 注意：
# 1. 不建议加入 colorectal cancer / cancer 等 outcome 本身关键词
# 2. 不建议加入 inflammation / immune 等过宽关键词，容易误删真实生物学通路
keyword_pattern <- paste(
  c(
    "body mass index",
    "bmi",
    "obesity",
    "overweight",
    "adiposity",
    "waist circumference",
    "waist-to-hip ratio",
    
    "smoking",
    "tobacco",
    "cigarette",
    
    "alcohol",
    "drinking",
    
    "diet",
    "dietary",
    "fiber",
    "fibre",
    "red meat",
    "processed meat",
    "meat intake",
    "vegetable",
    "fruit",
    
    "physical activity",
    "sedentary",
    "exercise",
    
    "diabetes",
    "type 2 diabetes",
    "insulin resistance",
    "glucose",
    
    "inflammatory bowel disease",
    "crohn",
    "crohn's disease",
    "ulcerative colitis",
    
    "aspirin",
    "nonsteroidal anti-inflammatory",
    "nsaid",
    
    "height",
    "socioeconomic",
    "education"
  ),
  collapse = "|"
)

# 在本地 GWAS Catalog 中检索这些 SNP
catalog_lookup <- load_local_gwas_catalog_matches(
  gwas_catalog_assoc_path,
  unique(strong_instruments$SNP),
  keyword_pattern
)

gc_dt <- catalog_lookup$matches
cache_audit_dt <- catalog_lookup$audit

# 需要剔除的 SNP
confounder_hits <- if (nrow(gc_dt)) {
  unique(gc_dt[confounder_filter_status == "remove", query_snp])
} else {
  character()
}

# 没匹配上的 SNP
cache_miss_dt <- if (nrow(cache_audit_dt)) {
  cache_audit_dt[match_status == "no_local_match"]
} else {
  data.table()
}

# 保存 GWAS Catalog 检索结果
fwrite(
  gc_dt,
  "gwas_catalog_clumped_snp_associations.tsv",
  sep = "\t"
)

fwrite(
  cache_audit_dt,
  "gwas_catalog_cache_audit.tsv",
  sep = "\t"
)

fwrite(
  cache_miss_dt,
  "gwas_catalog_cache_miss.tsv",
  sep = "\t"
)

# 保存命中混杂因素的 SNP 列表
fwrite(
  unique(cache_audit_dt[tolower(SNP) %in% confounder_hits, .(SNP)]),
  "gwas_catalog_confounder_hit_snps.tsv",
  sep = "\t"
)

# 去掉混杂 SNP 后的工具变量
after_confounder <- strong_instruments[
  !tolower(SNP) %in% confounder_hits
]

# 保存结果
fwrite(
  after_confounder,
  "exposure_instruments_after_confounder_filter.tsv",
  sep = "\t"
)

# 记录保留数量
trace_by_gene$after_confounder_filter_snps <- by_gene_count(after_confounder)

trace_overall[[length(trace_overall) + 1L]] <- overall_count(
  after_confounder,
  "after_confounder_filter_snps"
)

cat("GWAS Catalog 命中混杂因素的 SNP 数：", length(confounder_hits), "\n")
cat("混杂因素剔除后 SNP 记录数：", nrow(after_confounder), "\n")
cat("混杂因素剔除后基因数：", uniqueN(after_confounder$gene), "\n")
cat("混杂因素剔除后唯一 SNP 数：", uniqueN(after_confounder$SNP), "\n")

if (!nrow(after_confounder)) {
  stop("混杂因素剔除后没有保留工具变量，无法继续后续分析。")
}
#========================================================
# 24. Step 5：在 outcome VCF 中查询这些 SNP 是否已与 CRC 显著相关
#========================================================

message("Step 5: Query the local colorectal cancer VCF and remove SNPs already significantly associated with CRC")

# 检查 VCF 及索引文件是否存在
if (!file.exists(outcome_vcf) || !file.exists(paste0(outcome_vcf, ".tbi"))) {
  stop("Outcome VCF or its .tbi index is missing.")
}

# 读取这些 SNP 在 CRC 结局中的效应信息
outcome_assoc <- read_outcome_from_vcf(after_confounder, outcome_vcf)

# 保存 outcome 查询结果
fwrite(
  outcome_assoc,
  "ebi-a-GCST90018808_colorectal_cancer_clumped_snp_associations.tsv",
  sep = "\t"
)

cat("在 CRC outcome VCF 中匹配到的 SNP 数：", nrow(outcome_assoc), "\n")

#========================================================
# 25. 合并 exposure 和 outcome 信息
#========================================================

instrument_with_outcome <- merge(
  after_confounder,
  outcome_assoc[, .(
    SNP = rsid,
    beta.outcome = beta,
    se.outcome = se,
    pval.outcome = p,
    effect_allele.outcome = ea,
    other_allele.outcome = nea,
    eaf.outcome = eaf
  )],
  by = "SNP",
  all.x = TRUE
)

# 保存合并结果
fwrite(
  instrument_with_outcome,
  "exposure_instruments_with_crc_lookup.tsv",
  sep = "\t"
)

# 只保留在 outcome 中能查到的 SNP
with_outcome <- instrument_with_outcome[!is.na(beta.outcome)]

# 记录数量
trace_by_gene$with_crc_outcome_snps <- by_gene_count(with_outcome)

trace_overall[[length(trace_overall) + 1L]] <- overall_count(
  with_outcome,
  "with_crc_outcome_snps"
)

cat("合并 CRC outcome 后 SNP 记录数：", nrow(with_outcome), "\n")
cat("合并 CRC outcome 后基因数：", uniqueN(with_outcome$gene), "\n")
cat("合并 CRC outcome 后唯一 SNP 数：", uniqueN(with_outcome$SNP), "\n")

if (!nrow(with_outcome)) {
  stop("没有 SNP 在 CRC outcome VCF 中匹配到，无法继续。")
}

#========================================================
# 26. 去掉已经和 CRC 显著相关的 SNP
#========================================================

# 如果某个 SNP 本身就和 CRC 显著相关，p < 5e-8，
# 它可能不适合作为“仅通过暴露影响结局”的工具变量，因此剔除。
final_iv <- with_outcome[is.na(pval.outcome) | pval.outcome >= 5e-8]

# 保存最终可用于 MR 的工具变量
fwrite(
  final_iv,
  "exposure_instruments_for_mr.tsv",
  sep = "\t"
)

# 记录数量
trace_by_gene$after_crc_sig_filter_snps <- by_gene_count(final_iv)

trace_overall[[length(trace_overall) + 1L]] <- overall_count(
  final_iv,
  "after_crc_sig_filter_snps"
)

cat("去除与 CRC 显著相关 SNP 后记录数：", nrow(final_iv), "\n")
cat("去除与 CRC 显著相关 SNP 后基因数：", uniqueN(final_iv$gene), "\n")
cat("去除与 CRC 显著相关 SNP 后唯一 SNP 数：", uniqueN(final_iv$SNP), "\n")

if (!nrow(final_iv)) {
  stop("去除与 CRC 显著相关 SNP 后没有剩余工具变量。")
}

#========================================================
# 27. 按基因统计最终工具变量数量
#========================================================

exposure_gene_counts <- final_iv[, .(
  n_snps = as.integer(.N),
  mean_F = mean(F_stat, na.rm = TRUE),
  median_F = median(F_stat, na.rm = TRUE),
  instrument_source = instrument_source[1]
), by = gene]

# 是否满足多 SNP MR 的要求，这里定义为 >=3
exposure_gene_counts[, eligible_for_multi_snp_mr := n_snps >= 3]

# 是否至少能做单 SNP MR
exposure_gene_counts[, eligible_for_any_mr := n_snps >= 1]

# 保存基因级统计结果
fwrite(
  exposure_gene_counts,
  "exposure_gene_counts.tsv",
  sep = "\t"
)

# 保存可做多 SNP MR 的候选基因
fwrite(
  exposure_gene_counts[eligible_for_multi_snp_mr == TRUE],
  "pre_harmonization_iv_gene_candidates.tsv",
  sep = "\t"
)

cat("最终至少有 1 个 IV 的基因数：", sum(exposure_gene_counts$eligible_for_any_mr), "\n")
cat("最终有 >=3 个 IV 的基因数：", sum(exposure_gene_counts$eligible_for_multi_snp_mr), "\n")

#========================================================
# 28. Step 6：保存 SNP 和基因筛选轨迹
#========================================================

message("Step 6: Save SNP and gene filter traces")

# 基因轨迹表：每个基因从一开始到最后每一步剩多少 SNP
gene_trace <- unique(target_genes[, .(gene)])

trace_order <- c(
  "all_eqtl_rows",
  "genomewide_sig_snps",
  "ld_clumped_snps",
  "strong_f_snps",
  "after_confounder_filter_snps",
  "with_crc_outcome_snps",
  "after_crc_sig_filter_snps"
)

# 依次把每一步统计结果合并到 gene_trace
for (nm in trace_order) {
  gene_trace <- merge_trace_col(gene_trace, trace_by_gene[[nm]], nm)
}

# 是否还有最终 IV
gene_trace[, has_any_final_iv := after_crc_sig_filter_snps > 0]

# 是否满足多 SNP MR
gene_trace[, eligible_for_multi_snp_mr := after_crc_sig_filter_snps >= 3]

# 给每个基因打一个最终状态标签
gene_trace[, local_selection_stage := fifelse(
  eligible_for_multi_snp_mr,
  "eligible_for_multi_snp_mr",
  fifelse(has_any_final_iv, "eligible_for_any_mr", "dropped_before_mr")
)]

# 保存基因筛选轨迹
fwrite(
  gene_trace,
  "local_gene_filter_trace.tsv",
  sep = "\t"
)

# SNP 总体筛选轨迹表
overall_trace_dt <- rbindlist(trace_overall, fill = TRUE)
overall_trace_dt[, step_order := seq_len(.N)]
setcolorder(overall_trace_dt, c("step_order", "step", "n_gene_snp", "n_unique_snps", "n_genes"))

# 保存总体轨迹
fwrite(
  overall_trace_dt,
  "local_snp_filter_trace.tsv",
  sep = "\t"
)

#========================================================
# 29. 结束提示
#========================================================

message("Finished local MR subset preparation.")

cat("\n重点输出文件：\n")
cat("1. exposure_instruments_for_mr.tsv\n")
cat("2. exposure_gene_counts.tsv\n")
cat("3. pre_harmonization_iv_gene_candidates.tsv\n")
cat("4. local_gene_filter_trace.tsv\n")
cat("5. local_snp_filter_trace.tsv\n")
cat("6. exposure_instruments_with_crc_lookup.tsv\n")

shared_snp <- final_iv[, .(
  n_genes = uniqueN(gene),
  genes = paste(sort(unique(gene)), collapse = "; "),
  n_records = .N
), by = SNP][n_genes > 1]

# 查看结果
shared_snp

#!/usr/bin/env Rscript

#========================================================
# MR分析教学版脚本
# 作用：
# 1）读取暴露数据和结局数据
# 2）按基因逐个进行 harmonization（等位基因一致化）
# 3）执行多种 MR 方法
# 4）补充异质性、多效性、MR-PRESSO、leave-one-out、single SNP 敏感性分析
# 5）输出各类中间结果表和最终汇总表
#========================================================

#========================================================
# 0. 基础设置
#========================================================

## 将系统报错改成英文，方便复制报错去搜索
Sys.setenv(LANGUAGE = "en")

## 禁止字符串自动转成因子，避免后续数据处理出问题
options(stringsAsFactors = FALSE)

## 清空当前环境中的所有对象
rm(list = ls())

#========================================================
# 1. 检查并安装所需 R 包
#========================================================

## 本脚本依赖的核心R包
required_packages <- c("data.table", "TwoSampleMR", "MRPRESSO")

## 找出当前还没有安装的包
missing <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

## 其中可以直接从CRAN安装的包
installable_from_cran <- setdiff(missing, c("TwoSampleMR", "MRPRESSO"))

## 如果有缺失的CRAN包，则自动安装
if (length(installable_from_cran)) {
  install.packages(installable_from_cran, repos = "https://cloud.r-project.org")
}

## 再次检查是否还有缺失包
still_missing <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

## 如果还有包缺失，则直接停止，并提示用户先安装
if (length(still_missing)) {
  stop(sprintf(
    "Missing required packages: %s. Install them before running this script.",
    paste(still_missing, collapse = ", ")
  ))
}

#========================================================
# 2. 加载所需 R 包
#========================================================

## data.table：高效读写和处理大表
library(data.table)

## TwoSampleMR：进行MR分析、harmonization、异质性、多效性等
library(TwoSampleMR)

## MRPRESSO：用于检测离群工具变量和潜在水平多效性
library(MRPRESSO)

#========================================================
# 3. 定义项目路径和输入/输出文件路径
#========================================================

## 设置当前工作目录
workdir <- "C:/Users/92066/Desktop/单细胞孟德尔/MR/m7"
setwd(workdir)

## 定义结果输出目录
mr_output_dir <- file.path(workdir, "m7")

## 若输出目录不存在，则自动创建
dir.create(mr_output_dir, recursive = TRUE, showWarnings = FALSE)

## m5阶段生成的起始基因流程表
start_trace_path <- file.path(workdir, "gene_selection_trace.tsv")

## 如果没有 gene_selection_trace.tsv，则使用当前实际生成的 eqtlgen_gene_trace.tsv
if (!file.exists(start_trace_path)) {
  start_trace_path <- file.path(workdir, "eqtlgen_gene_trace.tsv")
}

## m6阶段生成的暴露工具变量文件
exposure_path <- file.path(workdir, "exposure_instruments_for_mr.tsv")

## m6阶段生成的结直肠癌结局SNP关联文件
outcome_path <- file.path(
  workdir,
  "ebi-a-GCST90018808_colorectal_cancer_clumped_snp_associations.tsv"
)

## m6阶段本地筛选后的基因追踪表
local_gene_trace_path <- file.path(workdir, "local_gene_filter_trace.tsv")

## m6阶段harmonization前每个基因的SNP数量/F值统计表
preharm_counts_path <- file.path(workdir, "exposure_gene_counts.tsv")

#========================================================
# 4. 检查输入文件是否存在
#========================================================

## 如果暴露文件不存在，则停止运行
if (!file.exists(exposure_path)) {
  stop("Exposure instrument file is missing. Run prepare_mr_subsets_local.R first.")
}

## 如果结局文件不存在，则停止运行
if (!file.exists(outcome_path)) {
  stop("Outcome association file is missing. Run prepare_mr_subsets_local.R first.")
}

#========================================================
# 5. 读取输入数据
#========================================================

## 读取暴露工具变量数据
exposure <- fread(exposure_path)

## 读取结局GWAS关联数据
outcome <- fread(outcome_path)

## 如果起始流程表存在，则读取；否则给一个空表
start_trace <- if (file.exists(start_trace_path)) fread(start_trace_path) else data.table()

## 如果本地基因筛选追踪表存在，则读取；否则给一个空表
local_gene_trace <- if (file.exists(local_gene_trace_path)) fread(local_gene_trace_path) else data.table()

## 如果harmonization前统计表存在，则读取；否则给一个空表
preharm_counts <- if (file.exists(preharm_counts_path)) fread(preharm_counts_path) else data.table()

## 兼容 gene_clean 列名
if (nrow(start_trace) && !"gene" %in% names(start_trace) && "gene_clean" %in% names(start_trace)) {
  setnames(start_trace, "gene_clean", "gene")
}

if (nrow(local_gene_trace) && !"gene" %in% names(local_gene_trace) && "gene_clean" %in% names(local_gene_trace)) {
  setnames(local_gene_trace, "gene_clean", "gene")
}

if (nrow(preharm_counts) && !"gene" %in% names(preharm_counts) && "gene_clean" %in% names(preharm_counts)) {
  setnames(preharm_counts, "gene_clean", "gene")
}

#========================================================
# 6. 定义辅助函数
#========================================================

## safe_rbind：
## 作用：把列表中的多个data.table安全合并
## 如果列表为空，则返回一个空data.table，避免报错
safe_rbind <- function(x) {
  if (!length(x)) {
    return(data.table())
  }
  rbindlist(x, fill = TRUE)
}

## build_exposure_dat：
## 作用：将原始暴露数据整理成TwoSampleMR可识别的标准格式
build_exposure_dat <- function(dt) {
  
  out <- copy(dt[, .(
    gene,
    id.exposure,
    SNP,
    beta.exposure,
    se.exposure,
    pval.exposure,
    samplesize.exposure,
    effect_allele.exposure,
    other_allele.exposure,
    eaf.exposure,
    F_stat,
    instrument_source
  )])
  
  ## 把 id.exposure 强制转为字符型
  out[, id.exposure := as.character(id.exposure)]
  
  ## 如果 id.exposure 缺失，则自动按“local-eqtlgen-基因名”补齐
  out[, id.exposure := fifelse(
    !is.na(id.exposure) & nzchar(id.exposure),
    id.exposure,
    paste0("local-eqtlgen-", gene)
  )]
  
  ## TwoSampleMR中 exposure 列通常用于展示暴露名称，这里直接使用基因名
  out[, exposure := gene]
  
  ## 等位基因统一转成大写，避免后续匹配时出现大小写不一致
  out[, effect_allele.exposure := toupper(effect_allele.exposure)]
  out[, other_allele.exposure := toupper(other_allele.exposure)]
  
  ## 把关键统计量都强制转成数值型，避免字符型导致后续函数报错
  out[, samplesize.exposure := suppressWarnings(as.numeric(samplesize.exposure))]
  out[, pval.exposure := suppressWarnings(as.numeric(pval.exposure))]
  out[, beta.exposure := suppressWarnings(as.numeric(beta.exposure))]
  out[, se.exposure := suppressWarnings(as.numeric(se.exposure))]
  out[, eaf.exposure := suppressWarnings(as.numeric(eaf.exposure))]
  out[]
}

## build_outcome_lookup：
## 作用：把原始结局GWAS表整理成标准格式，并建立SNP查找表
build_outcome_lookup <- function(dt) {
  
  out <- unique(dt[, .(
    SNP = as.character(rsid),
    beta.outcome = suppressWarnings(as.numeric(beta)),
    se.outcome = suppressWarnings(as.numeric(se)),
    pval.outcome = suppressWarnings(as.numeric(p)),
    effect_allele.outcome = toupper(as.character(ea)),
    other_allele.outcome = toupper(as.character(nea)),
    eaf.outcome = suppressWarnings(as.numeric(eaf)),
    chr.outcome = suppressWarnings(as.integer(chr)),
    pos.outcome = suppressWarnings(as.integer(position)),
    samplesize.outcome = suppressWarnings(as.numeric(n))
  )], by = "SNP")
  
  ## 补充结局ID
  out[, id.outcome := "local-ebi-a-GCST90018808"]
  
  ## 补充结局名称
  out[, outcome := "Colorectal cancer"]
  
  out[]
}

## run_methods_for_gene：
## 作用：根据某个基因最终保留下来的SNP数，选择合适的MR方法进行分析
run_methods_for_gene <- function(dat_gene, primary_ivw_method) {
  
  nsnp <- nrow(dat_gene)
  
  ## 根据SNP数决定可以使用哪些MR方法
  ## 1个SNP：只能做Wald ratio
  ## 2个SNP：通常只跑IVW
  ## >=3个SNP：可以跑Egger、加权中位数、IVW、mode方法
  methods <- if (nsnp == 1) {
    c("mr_wald_ratio")
  } else if (nsnp == 2) {
    c(primary_ivw_method)
  } else {
    c(
      "mr_egger_regression",
      "mr_weighted_median",
      primary_ivw_method,
      "mr_simple_mode",
      "mr_weighted_mode"
    )
  }
  
  out_list <- list()
  
  ## 逐个方法执行MR分析
  for (method_id in methods) {
    
    res <- tryCatch(
      as.data.table(mr(dat_gene, method_list = method_id)),
      error = function(e) data.table()
    )
    
    ## 如果该方法没有成功返回结果，则跳过
    if (!nrow(res)) {
      next
    }
    
    ## 记录方法ID
    res[, method_id := method_id]
    out_list[[method_id]] <- res
  }
  
  ## 合并所有方法结果
  rbindlist(out_list, fill = TRUE)
}

## pick_primary_result：
## 作用：从所有MR方法结果中选出“主结果”
## 逻辑：
## 1个SNP时，主结果为Wald ratio
## 多个SNP时，优先选择IVW结果
pick_primary_result <- function(method_dt, nsnp) {
  
  if (!nrow(method_dt)) {
    return(data.table())
  }
  
  ## 单SNP情形：优先取Wald ratio
  if (nsnp == 1) {
    out <- method_dt[grepl("^Wald ratio", method)][1]
    if (!nrow(out)) {
      out <- method_dt[1]
    }
    return(out)
  }
  
  ## 多SNP情形：优先取IVW结果
  out <- method_dt[grepl("^Inverse variance weighted", method)][1]
  
  if (!nrow(out)) {
    out <- method_dt[1]
  }
  
  out
}

## choose_ivw_method：
## 作用：根据异质性检验结果，自动决定主分析使用固定效应IVW还是随机效应IVW
## 逻辑：
## 若异质性P >= 0.05，用固定效应IVW（mr_ivw_fe）
## 若异质性P < 0.05，用随机效应IVW（mr_ivw_mre）
choose_ivw_method <- function(dat_gene) {
  
  if (nrow(dat_gene) < 2) {
    return(list(
      method = NA_character_,
      heterogeneity = data.table(),
      heterogeneity_p = NA_real_
    ))
  }
  
  ## 计算异质性
  het <- tryCatch(
    as.data.table(mr_heterogeneity(dat_gene, method_list = c("mr_ivw", "mr_egger_regression"))),
    error = function(e) data.table()
  )
  
  ## 提取IVW异质性检验的Q_p值
  ivw_het_p <- if (nrow(het) && "method" %in% names(het) && "Q_pval" %in% names(het)) {
    het[grepl("^Inverse variance weighted", method), Q_pval][1]
  } else {
    NA_real_
  }
  
  ## 根据异质性结果选择固定效应或随机效应IVW
  if (is.na(ivw_het_p) || ivw_het_p >= 0.05) {
    return(list(
      method = "mr_ivw_fe",
      heterogeneity = het,
      heterogeneity_p = ivw_het_p
    ))
  }
  
  list(
    method = "mr_ivw_mre",
    heterogeneity = het,
    heterogeneity_p = ivw_het_p
  )
}

#========================================================
# 7. 构建标准化后的暴露和结局数据
#========================================================

## 标准化暴露数据
exposure_dat <- build_exposure_dat(exposure)

## 标准化结局数据
outcome_lookup <- build_outcome_lookup(outcome)

#========================================================
# 8. 初始化结果保存列表
#========================================================

## 保存每个基因 harmonization 后的数据
harm_all <- list()

## 保存每个基因所有MR方法的结果
method_results_all <- list()

## 保存每个基因的主结果
primary_results_all <- list()

## 保存异质性结果
heterogeneity_all <- list()

## 保存水平多效性结果
pleiotropy_all <- list()

## 保存MR-PRESSO结果
presso_all <- list()

## 保存leave-one-out结果
loo_all <- list()

## 保存single SNP分析结果
single_all <- list()

## 保存每个基因工具变量保留情况统计
instrument_summary_all <- list()

#========================================================
# 9. 核心分析：按基因逐个做 harmonization 和 MR
#========================================================

message("Step 1: Harmonize local exposure and outcome data gene by gene")

## 遍历每一个基因
for (gene_name in unique(exposure_dat$gene)) {
  
  ## 提取该基因对应的暴露工具变量
  exp_gene <- exposure_dat[gene == gene_name]
  
  ## 取出该基因涉及的SNP，并与结局查找表合并
  ## 目的是保留同时在暴露和结局中都存在的SNP
  out_gene <- merge(
    unique(exp_gene[, .(SNP)]),
    outcome_lookup,
    by = "SNP",
    all.x = FALSE,
    all.y = FALSE
  )
  
  ## 如果没有任何SNP能在结局数据中找到，则记录并跳过
  if (!nrow(out_gene)) {
    
    instrument_summary_all[[gene_name]] <- data.table(
      gene = gene_name,
      instrument_source = exp_gene$instrument_source[1],
      total_input = nrow(exp_gene),
      with_outcome = 0L,
      outcome_sig_removed = 0L,
      palindromic_removed = 0L,
      allele_mismatch_removed = 0L,
      weak_removed = sum(is.na(exp_gene$F_stat) | exp_gene$F_stat <= 10, na.rm = TRUE),
      kept_for_mr = 0L,
      eligible_for_any_mr_after_harmonization = FALSE,
      eligible_for_multi_snp_mr_after_harmonization = FALSE
    )
    
    next
  }
  
  ## 进行等位基因一致化（harmonization）
  ## action = 3 通常表示严格处理不确定方向的回文SNP
  harm <- tryCatch(
    as.data.table(
      harmonise_data(as.data.frame(exp_gene), as.data.frame(out_gene), action = 3)
    ),
    error = function(e) data.table()
  )
  
  ## 如果harmonization失败或结果为空，则记录并跳过
  if (!nrow(harm)) {
    
    instrument_summary_all[[gene_name]] <- data.table(
      gene = gene_name,
      instrument_source = exp_gene$instrument_source[1],
      total_input = nrow(exp_gene),
      with_outcome = nrow(out_gene),
      outcome_sig_removed = 0L,
      palindromic_removed = 0L,
      allele_mismatch_removed = 0L,
      weak_removed = sum(is.na(exp_gene$F_stat) | exp_gene$F_stat <= 10, na.rm = TRUE),
      kept_for_mr = 0L,
      eligible_for_any_mr_after_harmonization = FALSE,
      eligible_for_multi_snp_mr_after_harmonization = FALSE
    )
    
    next
  }
  
  ## 给harm结果补充gene列，便于后续合并
  harm[, gene := gene_name]
  
  ## 标记：如果结局侧该SNP本身已经达到全基因组显著，则后续考虑去掉
  harm[, outcome_sig_removed := !is.na(pval.outcome) & pval.outcome < 5e-8]
  
  ## 定义最终保留规则：
  ## 1）不能被harmonization标记为remove
  ## 2）不能是回文SNP
  ## 3）不能是结局显著SNP
  ## 4）F值必须>10，避免弱工具变量
  harm[, keep := remove == FALSE &
         palindromic == FALSE &
         !outcome_sig_removed &
         !is.na(F_stat) &
         F_stat > 10]
  
  ## 保存harmonization后的原始结果
  harm_all[[gene_name]] <- harm
  
  ## 只提取最终用于MR分析的SNP
  kept <- harm[keep == TRUE]
  
  ## 汇总每个基因在各个过滤步骤中的去留情况
  instrument_summary_all[[gene_name]] <- data.table(
    gene = gene_name,
    instrument_source = exp_gene$instrument_source[1],
    total_input = nrow(exp_gene),
    with_outcome = sum(!is.na(harm$beta.outcome)),
    outcome_sig_removed = sum(harm$outcome_sig_removed, na.rm = TRUE),
    palindromic_removed = sum(harm$palindromic, na.rm = TRUE),
    allele_mismatch_removed = sum(harm$remove & !harm$palindromic, na.rm = TRUE),
    weak_removed = sum(is.na(harm$F_stat) | harm$F_stat <= 10, na.rm = TRUE),
    kept_for_mr = nrow(kept),
    eligible_for_any_mr_after_harmonization = nrow(kept) >= 1,
    eligible_for_multi_snp_mr_after_harmonization = nrow(kept) >= 3
  )
  
  ## 如果没有SNP保留下来，则不做MR
  if (!nrow(kept)) {
    next
  }
  
  ## 根据异质性自动选择主分析IVW类型
  ivw_choice <- choose_ivw_method(kept)
  
  ## 记录最终进入MR的SNP数
  nsnp_final <- nrow(kept)
  
  ## 执行MR主分析和其他方法
  method_dt_gene <- run_methods_for_gene(kept, ivw_choice$method)
  
  ## 若有MR结果，则进一步整理
  if (nrow(method_dt_gene)) {
    
    ## 把beta转成OR及95%CI
    method_dt_gene <- as.data.table(generate_odds_ratios(as.data.frame(method_dt_gene)))
    
    ## 补充说明列
    method_dt_gene[, gene := gene_name]
    method_dt_gene[, nsnp_final := nsnp_final]
    method_dt_gene[, selected_primary_ivw_method := ivw_choice$method]
    method_dt_gene[, ivw_heterogeneity_p := ivw_choice$heterogeneity_p]
    
    ## 保存所有方法结果
    method_results_all[[gene_name]] <- method_dt_gene
    
    ## 选出主结果
    primary_dt_gene <- pick_primary_result(method_dt_gene, nsnp_final)
    
    if (nrow(primary_dt_gene)) {
      primary_dt_gene[, gene := gene_name]
      primary_dt_gene[, nsnp_final := nsnp_final]
      primary_dt_gene[, selected_primary_ivw_method := ivw_choice$method]
      primary_dt_gene[, ivw_heterogeneity_p := ivw_choice$heterogeneity_p]
      primary_results_all[[gene_name]] <- primary_dt_gene
    }
  }
  
  ## 保存异质性结果
  if (nrow(ivw_choice$heterogeneity)) {
    het_dt <- copy(ivw_choice$heterogeneity)
    het_dt[, gene := gene_name]
    heterogeneity_all[[gene_name]] <- het_dt
  }
  
  ## SNP数>=3时，可以做Egger截距检验以评估水平多效性
  pleio_dt <- if (nsnp_final >= 3) {
    tryCatch(
      as.data.table(mr_pleiotropy_test(kept)),
      error = function(e) data.table()
    )
  } else {
    data.table()
  }
  
  ## 如果成功得到多效性结果，则保存
  if (nrow(pleio_dt)) {
    
    pleio_dt[, gene := gene_name]
    pleiotropy_all[[gene_name]] <- pleio_dt
    
  } else {
    
    ## 否则给出占位记录，便于后续汇总
    pleiotropy_all[[gene_name]] <- data.table(
      gene = gene_name,
      egger_intercept = NA_real_,
      se = NA_real_,
      pval = NA_real_,
      pleiotropy_status = if (nsnp_final >= 3) "no_result" else "skipped_lt3_snps"
    )
  }
  
  ## SNP数>=4时，做MR-PRESSO
  if (nsnp_final >= 4) {
    
    presso_obj <- tryCatch(
      mr_presso(
        BetaOutcome = "beta.outcome",
        BetaExposure = "beta.exposure",
        SdOutcome = "se.outcome",
        SdExposure = "se.exposure",
        data = as.data.frame(kept),
        OUTLIERtest = TRUE,
        DISTORTIONtest = TRUE,
        SignifThreshold = 0.05,
        NbDistribution = 1000,
        seed = 123
      ),
      error = function(e) e
    )
    
    ## 如果MR-PRESSO报错，则记录报错信息
    if (inherits(presso_obj, "error")) {
      
      presso_all[[gene_name]] <- data.table(
        gene = gene_name,
        presso_status = conditionMessage(presso_obj),
        presso_main_p = NA_real_,
        presso_global_p = NA_real_
      )
      
    } else {
      
      ## 提取MR-PRESSO主结果P值和全局检验P值
      main_p <- tryCatch(
        presso_obj$`Main MR results`$`P-value`[1],
        error = function(e) NA_real_
      )
      
      global_p <- tryCatch(
        presso_obj$`MR-PRESSO results`$`Global Test`$Pvalue,
        error = function(e) NA_real_
      )
      
      presso_all[[gene_name]] <- data.table(
        gene = gene_name,
        presso_status = "ok",
        presso_main_p = main_p,
        presso_global_p = global_p
      )
    }
    
  } else {
    
    ## SNP不足4个则跳过MR-PRESSO
    presso_all[[gene_name]] <- data.table(
      gene = gene_name,
      presso_status = "skipped_lt4_snps",
      presso_main_p = NA_real_,
      presso_global_p = NA_real_
    )
  }
  
  ## SNP数>=2时可以做leave-one-out分析
  if (nsnp_final >= 2) {
    
    loo_dt_gene <- tryCatch(
      as.data.table(mr_leaveoneout(kept)),
      error = function(e) data.table()
    )
    
    if (nrow(loo_dt_gene)) {
      loo_dt_gene[, gene := gene_name]
      loo_all[[gene_name]] <- loo_dt_gene
    }
  }
  
  ## single SNP分析：逐个SNP单独估计效应
  single_dt_gene <- tryCatch(
    as.data.table(mr_singlesnp(kept)),
    error = function(e) data.table()
  )
  
  if (nrow(single_dt_gene)) {
    single_dt_gene[, gene := gene_name]
    single_all[[gene_name]] <- single_dt_gene
  }
}

#========================================================
# 10. 合并所有基因的分析结果
#========================================================

## 把每个基因的结果列表合并成总表
harm_dt <- safe_rbind(harm_all)
method_dt <- safe_rbind(method_results_all)
primary_dt <- safe_rbind(primary_results_all)

## 对所有基因的主MR P值进行Benjamini-Hochberg多重检验校正
if (nrow(primary_dt)) {
  primary_dt[, MR_FDR := p.adjust(pval, method = "BH")]
}

heterogeneity_dt <- safe_rbind(heterogeneity_all)
pleiotropy_dt <- safe_rbind(pleiotropy_all)
presso_dt <- safe_rbind(presso_all)
loo_dt <- safe_rbind(loo_all)
single_dt <- safe_rbind(single_all)
instrument_summary <- safe_rbind(instrument_summary_all)

#========================================================
# 11. 导出各类中间结果和最终结果
#========================================================

## 定义本次MR分析结果输出目录
mr_local_dir <- mr_output_dir

## 导出harmonization后全部SNP表
fwrite(harm_dt, file.path(mr_local_dir, "mr_harmonized_instruments.tsv"), sep = "\t")

## 导出真正用于MR分析的SNP表
fwrite(harm_dt[keep == TRUE], file.path(mr_local_dir, "mr_analysis_ready.tsv"), sep = "\t")

## 导出每个基因工具变量筛选过程汇总表
fwrite(instrument_summary, file.path(mr_local_dir, "mr_instrument_summary.tsv"), sep = "\t")

## 导出全部MR方法结果
fwrite(method_dt, file.path(mr_local_dir, "mr_results_all_methods.tsv"), sep = "\t")

## 导出主结果表（IVW或Wald）
fwrite(primary_dt, file.path(mr_local_dir, "mr_results_ivw_wald.tsv"), sep = "\t")

## 导出异质性检验结果
fwrite(heterogeneity_dt, file.path(mr_local_dir, "mr_heterogeneity.tsv"), sep = "\t")

## 导出水平多效性检验结果
fwrite(pleiotropy_dt, file.path(mr_local_dir, "mr_pleiotropy.tsv"), sep = "\t")

## 导出MR-PRESSO结果
fwrite(presso_dt, file.path(mr_local_dir, "mr_presso_summary.tsv"), sep = "\t")

## 导出leave-one-out结果
fwrite(loo_dt, file.path(mr_local_dir, "mr_leaveoneout.tsv"), sep = "\t")

## 导出single SNP结果
fwrite(single_dt, file.path(mr_local_dir, "mr_singlesnp.tsv"), sep = "\t")

#========================================================
# 12. 输出“有资格进入MR”的基因列表
#========================================================

## 至少保留1个SNP、可以进入MR分析的基因
fwrite(
  instrument_summary[eligible_for_any_mr_after_harmonization == TRUE],
  file.path(mr_local_dir, "iv_eligible_genes.tsv"),
  sep = "\t"
)

## 至少保留3个SNP、可以进行较完整敏感性分析的基因
fwrite(
  instrument_summary[eligible_for_multi_snp_mr_after_harmonization == TRUE],
  file.path(mr_local_dir, "full_sensitivity_analysis_genes.tsv"),
  sep = "\t"
)

#========================================================
# 13. 构建基因工作流总表
#========================================================

## 如果有m5阶段起始追踪表，则以它为基础
## 否则就直接从暴露表中提取gene列表
workflow_dt <- if (nrow(start_trace)) {
  unique(start_trace, by = "gene")
} else {
  unique(exposure[, .(gene)])
}

## 合并本地基因筛选轨迹表
if (nrow(local_gene_trace)) {
  workflow_dt <- merge(workflow_dt, local_gene_trace, by = "gene", all.x = TRUE)
}

## 合并harmonization前每个基因的SNP数量和F值统计
if (nrow(preharm_counts)) {
  workflow_dt <- merge(
    workflow_dt,
    preharm_counts[, .(
      gene,
      preharm_n_snps = n_snps,
      preharm_mean_F = mean_F,
      preharm_median_F = median_F
    )],
    by = "gene",
    all.x = TRUE
  )
}

## 合并工具变量筛选结果统计表
workflow_dt <- merge(workflow_dt, instrument_summary, by = "gene", all.x = TRUE)

#========================================================
# 14. 合并主MR结果
#========================================================

if (nrow(primary_dt)) {
  
  primary_keep <- unique(primary_dt[, .(
    gene,
    primary_method = method,
    primary_method_id = method_id,
    nsnp = nsnp_final,
    beta = b,
    se = se,
    p = pval,
    MR_FDR,
    OR = or,
    lci95 = or_lci95,
    uci95 = or_uci95,
    selected_primary_ivw_method,
    ivw_heterogeneity_p
  )], by = "gene")
  
  workflow_dt <- merge(workflow_dt, primary_keep, by = "gene", all.x = TRUE)
}

#========================================================
# 15. 合并多效性和MR-PRESSO结果
#========================================================

## 合并Egger截距检验结果
if (nrow(pleiotropy_dt)) {
  
  pleio_keep <- unique(pleiotropy_dt[, .(
    gene,
    egger_intercept,
    egger_intercept_se = se,
    egger_intercept_pval = pval,
    pleiotropy_status = if ("pleiotropy_status" %in% names(pleiotropy_dt)) pleiotropy_status else NA_character_
  )], by = "gene")
  
  workflow_dt <- merge(workflow_dt, pleio_keep, by = "gene", all.x = TRUE)
}

## 合并MR-PRESSO结果
if (nrow(presso_dt)) {
  workflow_dt <- merge(workflow_dt, presso_dt, by = "gene", all.x = TRUE)
}

#========================================================
# 16. 补齐缺失列，避免后续判断时报错
#========================================================

## 这些字符型列如果不存在，则补成NA字符
for (col in c(
  "primary_method",
  "primary_method_id",
  "selected_primary_ivw_method",
  "pleiotropy_status",
  "presso_status"
)) {
  if (!col %in% names(workflow_dt)) {
    workflow_dt[, (col) := NA_character_]
  }
}

## 这些数值型列如果不存在，则补成NA数值
for (col in c(
  "beta", "se", "p", "MR_FDR", "OR", "lci95", "uci95",
  "ivw_heterogeneity_p", "egger_intercept",
  "egger_intercept_se", "egger_intercept_pval",
  "presso_main_p", "presso_global_p"
)) {
  if (!col %in% names(workflow_dt)) {
    workflow_dt[, (col) := NA_real_]
  }
}

#========================================================
# 17. 把计数列中的NA填成0
#========================================================

## 这些列本质上都是计数列，缺失就按0处理
count_cols <- intersect(
  names(workflow_dt),
  c(
    "all_eqtl_rows",
    "genomewide_sig_snps",
    "ld_clumped_snps",
    "strong_f_snps",
    "after_confounder_filter_snps",
    "with_crc_outcome_snps",
    "after_crc_sig_filter_snps",
    "preharm_n_snps",
    "total_input",
    "with_outcome",
    "outcome_sig_removed",
    "palindromic_removed",
    "allele_mismatch_removed",
    "weak_removed",
    "kept_for_mr",
    "nsnp"
  )
)

for (col in count_cols) {
  workflow_dt[is.na(get(col)), (col) := 0L]
}

#========================================================
# 18. 把逻辑型列中的NA填成FALSE
#========================================================

## 这些列本质上表示“是/否”，缺失就按FALSE处理
logical_cols <- intersect(
  names(workflow_dt),
  c(
    "in_start_set",
    "has_genomewide_eqtl",
    "has_any_final_iv",
    "eligible_for_multi_snp_mr",
    "eligible_for_any_mr_after_harmonization",
    "eligible_for_multi_snp_mr_after_harmonization"
  )
)

for (col in logical_cols) {
  workflow_dt[is.na(get(col)), (col) := FALSE]
}

#========================================================
# 19. 定义最终筛选标签
#========================================================

## 主结果名义显著性：P < 0.05
workflow_dt[, primary_nominal_significant := !is.na(p) & p < 0.05]

## Benjamini-Hochberg校正后显著：FDR < 0.05
workflow_dt[, primary_FDR_significant := !is.na(MR_FDR) & MR_FDR < 0.05]

## 主结果是否显著：返修后以BH校正后的FDR < 0.05为准
workflow_dt[, primary_significant := primary_FDR_significant]

## 候选生物标志物：这里只按“主结果显著”来定义
workflow_dt[, candidate_biomarker := primary_significant]

## 记录候选标志物的筛选依据
workflow_dt[, biomarker_selection_note := fifelse(
  primary_FDR_significant,
  "selected_by_primary_FDR_significance",
  fifelse(
    primary_nominal_significant,
    "nominally_significant_exploratory",
    "primary_not_significant"
  )
)]

## 是否提示存在水平多效性
## Egger截距检验P < 0.05，说明可能存在水平多效性
workflow_dt[, horizontal_pleiotropy_flag :=
              !is.na(egger_intercept_pval) & egger_intercept_pval < 0.05
]

## 给每个基因打上一个工作流阶段标签
## 方便后续快速查看每个基因走到了哪一步
workflow_dt[, workflow_stage := fifelse(
  candidate_biomarker,
  "candidate_biomarker",
  fifelse(
    eligible_for_multi_snp_mr_after_harmonization,
    "full_sensitivity_analysis_gene",
    fifelse(
      eligible_for_any_mr_after_harmonization,
      "mr_entry_gene",
      "start_gene_only"
    )
  )
)]

#========================================================
# 20. 排序并导出最终汇总表
#========================================================

## 按P值和基因名排序，显著结果优先排在前面
setorder(workflow_dt, p, gene)

## 导出完整的基因工作流总表
fwrite(
  workflow_dt,
  file.path(mr_local_dir, "mr_gene_selection_workflow.tsv"),
  sep = "\t"
)

## 导出候选生物标志物表
fwrite(
  workflow_dt[candidate_biomarker == TRUE],
  file.path(mr_local_dir, "candidate_biomarker_panel.tsv"),
  sep = "\t"
)

## 导出显著MR基因表
fwrite(
  workflow_dt[primary_significant == TRUE],
  file.path(mr_local_dir, "significant_mr_gene_panel.tsv"),
  sep = "\t"
)

#========================================================
# 21. 结束提示
#========================================================

## 提示整个本地MR流程运行结束
message("Finished local MR analysis with multi-method MR and sensitivity analyses.")

############################################################
## M7.5 FINAL — scientifically prioritized MR + DEG integration
##
## 目的：
## 1) 主MR统计显著性：BH-FDR < 0.05
## 2) 工具变量强度：真正进入MR的每个SNP均要求 F > 10
## 3) 不再把 nsnp >= 3 作为最终候选的硬性筛选条件
## 4) IVW heterogeneity、MR-Egger intercept、MR-PRESSO、
##    多方法方向一致性、leave-one-out 均作为敏感性/稳健性证据，
##    不机械地全部AND作为淘汰条件
## 5) 异质性存在时不自动删除基因；M7主程序已根据Q检验选择
##    fixed-effect或multiplicative random-effects IVW
## 6) MR-Egger与MR-PRESSO按照SNP数和实际可运行性解释：
##    unavailable/not applicable != failed
## 7) Figure2 DEG结果在本步骤仅用于一致性/QC注释：
##    adjusted P/FDR < 0.05 且 |logFC| > 0.5
## 8) MR与表达方向一致性也仅作为整合证据记录：
##    OR > 1 + Tumor_up
##    OR < 1 + Tumor_down
##    不在共定位前作为硬性排除条件
## 9) 共定位分析前候选仅要求：
##    主MR BH-FDR < 0.05 + min F > 10
## 10) Colocalization 单独放在本步骤之后进行。
##     共定位完成后，再结合MR-DEG方向一致性进行最终优先级判断。
############################################################

rm(list = ls())
options(stringsAsFactors = FALSE)

library(dplyr)
library(readr)
library(stringr)

############################################################
## 0. 参数与路径
############################################################

## 主MR统计学阈值
MR_FDR_CUTOFF <- 0.05

## 工具变量强度阈值
MIN_F <- 10

## 以下P值阈值仅用于敏感性分析的“证据标注”
## 不作为机械的最终AND筛选条件
HETERO_P_CUTOFF <- 0.05
PLEIO_P_CUTOFF <- 0.05
PRESSO_P_CUTOFF <- 0.05

## DEG标准：与Figure2初始DEG定义一致
DEG_FDR_CUTOFF <- 0.05
DEG_ABS_LOGFC_CUTOFF <- 0.5

## 输入目录
mr_dir <- "C:/Users/92066/Desktop/单细胞孟德尔/MR/m7"
deg_dir <- "C:/Users/92066/Desktop/单细胞孟德尔/Figure2/results"

## 输出目录
out_dir <- "C:/Users/92066/Desktop/单细胞孟德尔/MR/m7.5孟德尔+差异基因"

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

############################################################
## 1. 工具函数
############################################################

## 统一基因列名
standardize_gene_col <- function(df) {
  
  cn <- colnames(df)
  
  if ("Gene" %in% cn) {
    # 不处理
  } else if ("gene" %in% cn) {
    df <- df %>% rename(Gene = gene)
  } else if ("exposure" %in% cn) {
    df <- df %>% rename(Gene = exposure)
  } else if ("symbol" %in% cn) {
    df <- df %>% rename(Gene = symbol)
  } else if ("SYMBOL" %in% cn) {
    df <- df %>% rename(Gene = SYMBOL)
  } else {
    stop("没有识别到 Gene/gene/exposure/symbol/SYMBOL 列，请检查输入文件。")
  }
  
  df %>%
    mutate(
      Gene = trimws(as.character(Gene)),
      Gene = gsub("^local-eqtlgen-", "", Gene)
    )
}


## 检查文件是否存在
check_file <- function(path) {
  
  if (!file.exists(path)) {
    stop(paste0("没有找到文件：\n", path))
  }
}


## 安全最小值
safe_min <- function(x) {
  
  x <- suppressWarnings(as.numeric(x))
  
  if (length(x) == 0 || all(is.na(x))) {
    return(NA_real_)
  }
  
  min(x, na.rm = TRUE)
}


## 安全均值
safe_mean <- function(x) {
  
  x <- suppressWarnings(as.numeric(x))
  
  if (length(x) == 0 || all(is.na(x))) {
    return(NA_real_)
  }
  
  mean(x, na.rm = TRUE)
}


## MR-PRESSO Global P 有时可能是 "<0.001" 形式
parse_p_numeric <- function(x) {
  
  y <- trimws(as.character(x))
  
  y[y %in% c("", "NA", "NaN", "NULL")] <- NA_character_
  
  suppressWarnings(
    as.numeric(
      gsub("^[<>]=?\\s*", "", y)
    )
  )
}


## 判断P值是否 > cutoff；兼容 "<0.001" / ">0.05"
p_greater_than <- function(x, cutoff = 0.05) {
  
  y <- trimws(as.character(x))
  y[y %in% c("", "NA", "NaN", "NULL")] <- NA_character_
  
  num <- parse_p_numeric(y)
  
  out <- rep(NA, length(y))
  
  plain <- !is.na(y) & !grepl("^[<>]", y)
  gt    <- !is.na(y) & grepl("^>", y)
  lt    <- !is.na(y) & grepl("^<", y)
  
  idx_plain <- which(plain & !is.na(num))
  if (length(idx_plain) > 0) {
    out[idx_plain] <- num[idx_plain] > cutoff
  }
  
  ## >0.05 / >0.10 均可判为通过
  idx_gt <- which(gt & !is.na(num))
  if (length(idx_gt) > 0) {
    out[idx_gt] <- num[idx_gt] >= cutoff
  }
  
  ## <0.001 可明确判为不通过；
  ## 若是 <0.10，则无法确定是否大于0.05，返回NA
  idx_lt <- which(lt & !is.na(num))
  if (length(idx_lt) > 0) {
    out[idx_lt] <- ifelse(
      num[idx_lt] <= cutoff,
      FALSE,
      NA
    )
  }
  
  out
}

############################################################
## 2. 定义输入文件
############################################################

ivw_file <- file.path(
  mr_dir,
  "mr_results_ivw_wald.tsv"
)

all_mr_file <- file.path(
  mr_dir,
  "mr_results_all_methods.tsv"
)

pleio_file <- file.path(
  mr_dir,
  "mr_pleiotropy.tsv"
)

presso_file <- file.path(
  mr_dir,
  "mr_presso_summary.tsv"
)

analysis_ready_file <- file.path(
  mr_dir,
  "mr_analysis_ready.tsv"
)

loo_file <- file.path(
  mr_dir,
  "mr_leaveoneout.tsv"
)

deg_file <- file.path(
  deg_dir,
  "Fig2AB_All_DEGs_Tumor_vs_Normal.csv"
)

## 必需文件
invisible(
  lapply(
    c(
      ivw_file,
      all_mr_file,
      pleio_file,
      presso_file,
      analysis_ready_file,
      deg_file
    ),
    check_file
  )
)

############################################################
## 3. 读取文件
############################################################

ivw <- read_tsv(
  ivw_file,
  show_col_types = FALSE
)

all_mr <- read_tsv(
  all_mr_file,
  show_col_types = FALSE
)

pleio <- read_tsv(
  pleio_file,
  show_col_types = FALSE
)

presso <- read_tsv(
  presso_file,
  show_col_types = FALSE
)

analysis_ready <- read_tsv(
  analysis_ready_file,
  show_col_types = FALSE
)

deg <- read_csv(
  deg_file,
  show_col_types = FALSE
)

## leave-one-out是敏感性分析文件：
## 若M7尚未输出该文件，不中断主流程，只标记为不可用
if (file.exists(loo_file)) {
  loo <- read_tsv(
    loo_file,
    show_col_types = FALSE
  )
} else {
  loo <- tibble()
}

## 统一Gene列
ivw <- standardize_gene_col(ivw)
all_mr <- standardize_gene_col(all_mr)
pleio <- standardize_gene_col(pleio)
presso <- standardize_gene_col(presso)
analysis_ready <- standardize_gene_col(analysis_ready)
deg <- standardize_gene_col(deg)

if (nrow(loo) > 0) {
  loo <- standardize_gene_col(loo)
}

############################################################
## 4. 整理主MR结果
############################################################

required_ivw_cols <- c(
  "Gene",
  "pval",
  "b",
  "se",
  "nsnp_final",
  "ivw_heterogeneity_p"
)

missing_ivw <- setdiff(
  required_ivw_cols,
  colnames(ivw)
)

if (length(missing_ivw) > 0) {
  stop(
    paste0(
      "mr_results_ivw_wald.tsv 缺少以下必要列：",
      paste(missing_ivw, collapse = ", ")
    )
  )
}

## 每个基因必须只有一个主MR结果
dup_ivw <- ivw %>%
  count(Gene, name = "n") %>%
  filter(n > 1)

if (nrow(dup_ivw) > 0) {
  
  print(dup_ivw)
  
  stop(
    "mr_results_ivw_wald.tsv 中同一基因存在多条主MR结果，请先检查。"
  )
}

ivw2 <- ivw %>%
  mutate(
    IVW_P = as.numeric(pval),
    beta = as.numeric(b),
    se = as.numeric(se),
    nsnp = as.numeric(nsnp_final),
    
    OR = if ("or" %in% colnames(ivw)) {
      as.numeric(or)
    } else {
      exp(beta)
    },
    
    OR_LCI95 = if ("or_lci95" %in% colnames(ivw)) {
      as.numeric(or_lci95)
    } else {
      exp(beta - 1.96 * se)
    },
    
    OR_UCI95 = if ("or_uci95" %in% colnames(ivw)) {
      as.numeric(or_uci95)
    } else {
      exp(beta + 1.96 * se)
    },
    
    ## 使用当前全部主MR基因重新进行BH校正
    MR_FDR = p.adjust(
      IVW_P,
      method = "BH"
    ),
    
    MR_direction = case_when(
      OR > 1 ~ "Risk",
      OR < 1 ~ "Protective",
      TRUE ~ "Neutral"
    ),
    
    pass_IVW =
      !is.na(IVW_P) &
      IVW_P < 0.05,
    
    pass_MR_FDR =
      !is.na(MR_FDR) &
      MR_FDR < MR_FDR_CUTOFF
  )

cat("\n============================================================\n")
cat("主MR基因数（用于BH-FDR校正）：", nrow(ivw2), "\n")
cat("============================================================\n")

############################################################
## 5. IVW异质性：只作为敏感性证据，不作为硬淘汰条件
##
## M7主程序已经依据 heterogeneity P：
## P >= 0.05 -> fixed-effect IVW
## P < 0.05  -> multiplicative random-effects IVW
############################################################

hetero2 <- ivw2 %>%
  transmute(
    Gene,
    heterogeneity_P =
      suppressWarnings(
        as.numeric(
          ivw_heterogeneity_p
        )
      ),
    heterogeneity_applicable =
      !is.na(nsnp) &
      nsnp >= 2,
    heterogeneity_flag =
      heterogeneity_applicable &
      !is.na(heterogeneity_P) &
      heterogeneity_P <= HETERO_P_CUTOFF,
    heterogeneity_status = case_when(
      heterogeneity_applicable == FALSE ~
        "not_applicable_lt2_snps",
      is.na(heterogeneity_P) ~
        "result_unavailable",
      heterogeneity_P <= HETERO_P_CUTOFF ~
        "heterogeneity_detected_random_effects_IVW",
      TRUE ~
        "no_significant_heterogeneity_fixed_effects_IVW"
    )
  )

############################################################
## 6. MR-Egger intercept：水平多效性敏感性分析
##
## 不能运行/无结果 -> not applicable / unavailable
## 不是“自动失败”
############################################################

required_pleio_cols <- c(
  "Gene",
  "pval",
  "pleiotropy_status"
)

missing_pleio <- setdiff(
  required_pleio_cols,
  colnames(pleio)
)

if (length(missing_pleio) > 0) {
  stop(
    paste0(
      "mr_pleiotropy.tsv 缺少以下必要列：",
      paste(missing_pleio, collapse = ", ")
    )
  )
}

if (!("egger_intercept" %in% colnames(pleio))) {
  pleio$egger_intercept <- NA_real_
}

dup_pleio <- pleio %>%
  count(Gene, name = "n") %>%
  filter(n > 1)

if (nrow(dup_pleio) > 0) {
  print(dup_pleio)
  stop(
    "mr_pleiotropy.tsv 中同一基因存在多条结果，请先检查。"
  )
}

pleio2 <- pleio %>%
  transmute(
    Gene,
    egger_intercept =
      suppressWarnings(
        as.numeric(
          egger_intercept
        )
      ),
    pleiotropy_P =
      suppressWarnings(
        as.numeric(
          pval
        )
      ),
    pleiotropy_status =
      as.character(
        pleiotropy_status
      ),
    pleiotropy_applicable =
      !str_detect(
        tolower(
          coalesce(
            pleiotropy_status,
            ""
          )
        ),
        "skipped_lt3_snps|skip|insufficient|not.?run"
      ),
    egger_pleiotropy_flag =
      pleiotropy_applicable &
      !is.na(pleiotropy_P) &
      pleiotropy_P <= PLEIO_P_CUTOFF,
    egger_status = case_when(
      pleiotropy_applicable == FALSE ~
        "not_applicable",
      is.na(pleiotropy_P) ~
        "result_unavailable",
      pleiotropy_P <= PLEIO_P_CUTOFF ~
        "directional_pleiotropy_signal",
      TRUE ~
        "no_evidence_directional_pleiotropy"
    )
  )

############################################################
## 7. MR-PRESSO：敏感性分析
##
## M7中 nsnp < 4 时记录 skipped_lt4_snps。
## 这里将其视为 not applicable，而不是失败。
## Global P <= 0.05 记为“存在global pleiotropy信号”，
## 但本步骤不机械删除该基因。
############################################################

required_presso_cols <- c(
  "Gene",
  "presso_status",
  "presso_global_p"
)

missing_presso <- setdiff(
  required_presso_cols,
  colnames(presso)
)

if (length(missing_presso) > 0) {
  stop(
    paste0(
      "mr_presso_summary.tsv 缺少以下必要列：",
      paste(missing_presso, collapse = ", ")
    )
  )
}

if (!("presso_main_p" %in% colnames(presso))) {
  presso$presso_main_p <- NA_character_
}

dup_presso <- presso %>%
  count(Gene, name = "n") %>%
  filter(n > 1)

if (nrow(dup_presso) > 0) {
  print(dup_presso)
  stop(
    "mr_presso_summary.tsv 中同一基因存在多条结果，请先检查。"
  )
}

presso2 <- presso %>%
  mutate(
    MR_PRESSO_status =
      as.character(
        presso_status
      ),
    MR_PRESSO_main_P_raw =
      as.character(
        presso_main_p
      ),
    MR_PRESSO_global_P_raw =
      as.character(
        presso_global_p
      ),
    MR_PRESSO_main_P =
      parse_p_numeric(
        presso_main_p
      ),
    MR_PRESSO_global_P =
      parse_p_numeric(
        presso_global_p
      ),
    presso_applicable =
      !str_detect(
        tolower(
          coalesce(
            MR_PRESSO_status,
            ""
          )
        ),
        "skipped_lt4_snps|skip|insufficient|not.?run"
      ),
    presso_global_pleiotropy_flag =
      case_when(
        presso_applicable == FALSE ~ FALSE,
        is.na(MR_PRESSO_global_P) ~ FALSE,
        MR_PRESSO_global_P <= PRESSO_P_CUTOFF ~ TRUE,
        TRUE ~ FALSE
      ),
    presso_status_interpretation = case_when(
      presso_applicable == FALSE ~
        "not_applicable",
      is.na(MR_PRESSO_global_P) ~
        "result_unavailable",
      MR_PRESSO_global_P <= PRESSO_P_CUTOFF ~
        "global_pleiotropy_signal",
      TRUE ~
        "no_evidence_global_pleiotropy"
    )
  ) %>%
  select(
    Gene,
    MR_PRESSO_status,
    MR_PRESSO_main_P_raw,
    MR_PRESSO_main_P,
    MR_PRESSO_global_P_raw,
    MR_PRESSO_global_P,
    presso_applicable,
    presso_global_pleiotropy_flag,
    presso_status_interpretation
  )

############################################################
## 8. F-statistic
##
## 关键修正：
## 不再从 mr_instrument_summary.tsv 找F值。
## 直接读取真正进入MR分析的 mr_analysis_ready.tsv，
## 按基因重新计算 mean/min/max F。
############################################################

required_ready_cols <- c(
  "Gene",
  "SNP",
  "F_stat"
)

missing_ready <- setdiff(
  required_ready_cols,
  colnames(analysis_ready)
)

if (length(missing_ready) > 0) {
  stop(
    paste0(
      "mr_analysis_ready.tsv 缺少以下必要列：",
      paste(missing_ready, collapse = ", ")
    )
  )
}

## 如果文件保留了keep列，再确认只使用最终keep==TRUE的SNP
if ("keep" %in% colnames(analysis_ready)) {
  
  analysis_ready <- analysis_ready %>%
    filter(
      tolower(
        as.character(
          keep
        )
      ) %in% c("true", "t", "1")
    )
}

## 同一个Gene-SNP不允许重复
dup_ready <- analysis_ready %>%
  count(
    Gene,
    SNP,
    name = "n"
  ) %>%
  filter(
    n > 1
  )

if (nrow(dup_ready) > 0) {
  
  print(dup_ready)
  
  stop(
    "mr_analysis_ready.tsv 中存在重复的 Gene-SNP 组合，请检查。"
  )
}

analysis_ready <- analysis_ready %>%
  mutate(
    F_stat_num =
      suppressWarnings(
        as.numeric(
          F_stat
        )
      )
  )

bad_F_rows <- analysis_ready %>%
  filter(
    is.na(
      F_stat_num
    )
  )

if (nrow(bad_F_rows) > 0) {
  
  print(
    bad_F_rows %>%
      select(
        Gene,
        SNP,
        F_stat
      )
  )
  
  stop(
    "mr_analysis_ready.tsv 中存在F_stat缺失或无法转为数值的SNP。"
  )
}

inst2 <- analysis_ready %>%
  group_by(
    Gene
  ) %>%
  summarise(
    n_F_SNP =
      n_distinct(
        SNP
      ),
    
    mean_F =
      mean(
        F_stat_num
      ),
    
    min_F =
      min(
        F_stat_num
      ),
    
    max_F =
      max(
        F_stat_num
      ),
    
    .groups = "drop"
  ) %>%
  mutate(
    pass_F =
      !is.na(
        min_F
      ) &
      min_F >
      MIN_F
  )

############################################################
## 9. 核对真实MR SNP数
##
## 主MR结果中的nsnp必须和mr_analysis_ready.tsv中
## 实际SNP数完全一致。
############################################################

snp_count_check <- ivw2 %>%
  select(
    Gene,
    nsnp
  ) %>%
  left_join(
    inst2 %>%
      select(
        Gene,
        n_F_SNP
      ),
    by = "Gene"
  ) %>%
  mutate(
    snp_count_match =
      !is.na(
        nsnp
      ) &
      !is.na(
        n_F_SNP
      ) &
      nsnp ==
      n_F_SNP
  )

bad_snp_count <- snp_count_check %>%
  filter(
    is.na(
      snp_count_match
    ) |
      snp_count_match ==
      FALSE
  )

if (nrow(bad_snp_count) > 0) {
  
  print(
    bad_snp_count
  )
  
  stop(
    paste0(
      "主MR结果中的nsnp与mr_analysis_ready.tsv中的实际SNP数不一致。\n",
      "为避免错误筛选，程序停止。"
    )
  )
}

############################################################
## 10. 多种MR方法方向：只作描述性稳健性证据
##
## 不再使用“>=80%才通过”的人为硬阈值。
############################################################

required_allmr_cols <- c(
  "Gene",
  "method_id",
  "nsnp",
  "b",
  "se",
  "pval"
)

missing_allmr <- setdiff(
  required_allmr_cols,
  colnames(all_mr)
)

if (length(missing_allmr) > 0) {
  stop(
    paste0(
      "mr_results_all_methods.tsv 缺少以下必要列：",
      paste(missing_allmr, collapse = ", ")
    )
  )
}

primary_direction <- ivw2 %>%
  transmute(
    Gene,
    primary_sign =
      sign(
        beta
      )
  )

method_direction <- all_mr %>%
  mutate(
    method_id =
      as.character(
        method_id
      ),
    beta_method =
      suppressWarnings(
        as.numeric(
          b
        )
      )
  ) %>%
  filter(
    method_id %in% c(
      "mr_egger_regression",
      "mr_weighted_median",
      "mr_ivw",
      "mr_ivw_fe",
      "mr_ivw_mre",
      "mr_simple_mode",
      "mr_weighted_mode",
      "mr_wald_ratio"
    ),
    !is.na(beta_method),
    beta_method != 0
  ) %>%
  group_by(
    Gene,
    method_id
  ) %>%
  summarise(
    method_sign = case_when(
      all(beta_method > 0) ~ 1,
      all(beta_method < 0) ~ -1,
      TRUE ~ 0
    ),
    .groups = "drop"
  )

direction2 <- method_direction %>%
  left_join(
    primary_direction,
    by = "Gene"
  ) %>%
  mutate(
    same_as_primary =
      !is.na(primary_sign) &
      primary_sign != 0 &
      method_sign == primary_sign
  ) %>%
  group_by(
    Gene
  ) %>%
  summarise(
    n_methods =
      n(),
    n_same_direction =
      sum(
        same_as_primary,
        na.rm = TRUE
      ),
    direction_consistency_ratio =
      ifelse(
        n_methods > 0,
        n_same_direction / n_methods,
        NA_real_
      ),
    .groups = "drop"
  ) %>%
  mutate(
    method_direction_note = case_when(
      is.na(direction_consistency_ratio) ~
        "unavailable",
      n_methods == 1 ~
        "single_method_only",
      direction_consistency_ratio == 1 ~
        "all_available_methods_concordant",
      direction_consistency_ratio > 0.5 ~
        "majority_methods_concordant",
      direction_consistency_ratio == 0.5 ~
        "mixed_directions",
      TRUE ~
        "majority_methods_discordant"
    )
  )

############################################################
## 11. Leave-one-out：只作敏感性证据，不作为硬筛选
############################################################

if (
  nrow(loo) > 0 &&
  all(
    c("Gene", "b") %in% colnames(loo)
  )
) {
  
  loo_work <- loo %>%
    mutate(
      loo_beta =
        suppressWarnings(
          as.numeric(
            b
          )
        )
    )
  
  if ("SNP" %in% colnames(loo_work)) {
    loo_work <- loo_work %>%
      filter(
        tolower(
          as.character(
            SNP
          )
        ) != "all"
      )
  }
  
  loo2 <- loo_work %>%
    filter(
      !is.na(loo_beta),
      loo_beta != 0
    ) %>%
    left_join(
      primary_direction,
      by = "Gene"
    ) %>%
    mutate(
      same_as_primary =
        !is.na(primary_sign) &
        primary_sign != 0 &
        sign(loo_beta) == primary_sign
    ) %>%
    group_by(
      Gene
    ) %>%
    summarise(
      n_loo =
        n(),
      n_loo_same_direction =
        sum(
          same_as_primary,
          na.rm = TRUE
        ),
      loo_direction_consistency_ratio =
        ifelse(
          n_loo > 0,
          n_loo_same_direction / n_loo,
          NA_real_
        ),
      loo_all_same_direction =
        ifelse(
          n_loo > 0,
          n_loo_same_direction == n_loo,
          NA
        ),
      .groups = "drop"
    )
  
} else {
  
  loo2 <- tibble(
    Gene = character(),
    n_loo = integer(),
    n_loo_same_direction = integer(),
    loo_direction_consistency_ratio = numeric(),
    loo_all_same_direction = logical()
  )
}

############################################################
## 12. 合并MR证据
##
## 重要：
## 共定位前“硬条件”只有：
##   1) 主MR BH-FDR < 0.05
##   2) 工具变量强度 min F > 10
##
## heterogeneity / Egger / PRESSO / 多方法方向 / LOO
## 均作为 sensitivity evidence 注释，不机械AND淘汰。
############################################################

mr_evidence_table <- ivw2 %>%
  left_join(
    hetero2,
    by = "Gene"
  ) %>%
  left_join(
    pleio2,
    by = "Gene"
  ) %>%
  left_join(
    presso2,
    by = "Gene"
  ) %>%
  left_join(
    inst2,
    by = "Gene"
  ) %>%
  left_join(
    direction2,
    by = "Gene"
  ) %>%
  left_join(
    loo2,
    by = "Gene"
  ) %>%
  mutate(
    pass_F =
      coalesce(
        pass_F,
        FALSE
      ),
    
    pass_primary_MR =
      pass_MR_FDR,
    
    pass_instrument_strength =
      pass_F,
    
    MR_significant_strong_instrument =
      pass_primary_MR &
      pass_instrument_strength,
    
    ## 下面均为“警告/注释”，不是自动淘汰条件
    heterogeneity_flag =
      coalesce(
        heterogeneity_flag,
        FALSE
      ),
    
    egger_pleiotropy_flag =
      coalesce(
        egger_pleiotropy_flag,
        FALSE
      ),
    
    presso_global_pleiotropy_flag =
      coalesce(
        presso_global_pleiotropy_flag,
        FALSE
      ),
    
    ## 按SNP数判断：本应可做却没有有效结果时，标记为信息不完整
    sensitivity_scope = case_when(
      is.na(nsnp) ~
        "unknown_instrument_count",
      nsnp == 1 ~
        "single_instrument_limited_sensitivity",
      nsnp == 2 ~
        "two_instruments_limited_sensitivity",
      nsnp == 3 ~
        "Egger_available_PRESSO_not_applicable",
      nsnp >= 4 ~
        "Egger_and_PRESSO_applicable",
      TRUE ~
        "unknown"
    ),
    
    sensitivity_incomplete =
      (
        !is.na(nsnp) &
          nsnp >= 2 &
          is.na(heterogeneity_P)
      ) |
      (
        !is.na(nsnp) &
          nsnp >= 2 &
          is.na(loo_direction_consistency_ratio)
      ) |
      (
        !is.na(nsnp) &
          nsnp >= 3 &
          (
            is.na(pleiotropy_applicable) |
              pleiotropy_applicable == FALSE |
              is.na(pleiotropy_P)
          )
      ) |
      (
        !is.na(nsnp) &
          nsnp >= 4 &
          (
            is.na(presso_applicable) |
              presso_applicable == FALSE |
              is.na(MR_PRESSO_global_P)
          )
      ),
    
    sensitivity_warning_count =
      as.integer(
        heterogeneity_flag
      ) +
      as.integer(
        egger_pleiotropy_flag
      ) +
      as.integer(
        presso_global_pleiotropy_flag
      ),
    
    sensitivity_summary = case_when(
      egger_pleiotropy_flag |
        presso_global_pleiotropy_flag ~
        "Caution: pleiotropy signal detected",
      
      sensitivity_incomplete ~
        "Sensitivity evidence incomplete; manual review required",
      
      heterogeneity_flag ~
        "Heterogeneity detected; M7 used random-effects IVW",
      
      TRUE ~
        "No major sensitivity-analysis warning"
    )
  ) %>%
  arrange(
    MR_FDR,
    IVW_P
  )

############################################################
## 13. BH-FDR显著的MR基因
############################################################

mr_significant_genes <- mr_evidence_table %>%
  filter(
    pass_primary_MR ==
      TRUE
  ) %>%
  arrange(
    MR_FDR,
    IVW_P
  )

############################################################
## 14. 整理Figure2 DEG
############################################################

## logFC列
if ("logFC" %in% colnames(deg)) {
  
  deg_fc_col <- "logFC"
  
} else if ("log2FC" %in% colnames(deg)) {
  
  deg_fc_col <- "log2FC"
  
} else if ("avg_log2FC" %in% colnames(deg)) {
  
  deg_fc_col <- "avg_log2FC"
  
} else if ("avg_logFC" %in% colnames(deg)) {
  
  deg_fc_col <- "avg_logFC"
  
} else {
  
  stop(
    "Figure2 DEG结果中没有找到 logFC/log2FC/avg_log2FC/avg_logFC 列。"
  )
}


## FDR列
if ("adj.P.Val" %in% colnames(deg)) {
  
  deg_fdr_col <- "adj.P.Val"
  
} else if ("adj_p_val" %in% colnames(deg)) {
  
  deg_fdr_col <- "adj_p_val"
  
} else if ("p_val_adj" %in% colnames(deg)) {
  
  deg_fdr_col <- "p_val_adj"
  
} else if ("FDR" %in% colnames(deg)) {
  
  deg_fdr_col <- "FDR"
  
} else if ("padj" %in% colnames(deg)) {
  
  deg_fdr_col <- "padj"
  
} else {
  
  stop(
    "Figure2 DEG结果中没有找到 adjusted P/FDR 列。"
  )
}


## nominal P列，可有可无
if ("P.Value" %in% colnames(deg)) {
  
  deg_p_col <- "P.Value"
  
} else if ("p_val" %in% colnames(deg)) {
  
  deg_p_col <- "p_val"
  
} else if ("pvalue" %in% colnames(deg)) {
  
  deg_p_col <- "pvalue"
  
} else if ("PValue" %in% colnames(deg)) {
  
  deg_p_col <- "PValue"
  
} else {
  
  deg_p_col <- NA_character_
}


deg2 <- deg %>%
  mutate(
    DEG_logFC =
      suppressWarnings(
        as.numeric(
          .data[[deg_fc_col]]
        )
      ),
    
    DEG_FDR =
      suppressWarnings(
        as.numeric(
          .data[[deg_fdr_col]]
        )
      )
  )

if (!is.na(deg_p_col)) {
  
  deg2 <- deg2 %>%
    mutate(
      DEG_P =
        suppressWarnings(
          as.numeric(
            .data[[deg_p_col]]
          )
        )
    )
  
} else {
  
  deg2$DEG_P <- NA_real_
}


## 不再像旧代码那样 distinct(Gene, .keep_all=TRUE) 默默取第一条
dup_deg <- deg2 %>%
  count(
    Gene,
    name = "n"
  ) %>%
  filter(
    n > 1
  )

if (nrow(dup_deg) > 0) {
  
  print(
    dup_deg
  )
  
  stop(
    paste0(
      "Figure2 DEG结果中存在重复Gene。\n",
      "请先检查probe-to-gene整理过程，不能自动保留第一条。"
    )
  )
}


deg_simple <- deg2 %>%
  select(
    Gene,
    DEG_logFC,
    DEG_P,
    DEG_FDR
  )

############################################################
## 15. 所有MR基因合并Figure2 DEG结果
##
## DEG阈值不是新的独立证据来源，因为MR暴露本身来自DEG集合；
## 此处保留作为一致性/QC核查。
## MR 与表达方向一致性在本阶段仅记录；
## 共定位完成后再用于最终整合优先级判断。
############################################################

mr_deg_all <- mr_evidence_table %>%
  left_join(
    deg_simple,
    by = "Gene"
  ) %>%
  mutate(
    DEG_direction = case_when(
      is.na(DEG_logFC) ~
        NA_character_,
      DEG_logFC > 0 ~
        "Tumor_up",
      DEG_logFC < 0 ~
        "Tumor_down",
      TRUE ~
        "No_change"
    ),
    
    pass_DEG =
      !is.na(DEG_FDR) &
      !is.na(DEG_logFC) &
      DEG_FDR < DEG_FDR_CUTOFF &
      abs(DEG_logFC) > DEG_ABS_LOGFC_CUTOFF,
    
    ## MR与表达方向一致：
    ## OR > 1 + tumor up
    ## OR < 1 + tumor down
    MR_DEG_consistent = case_when(
      is.na(OR) |
        is.na(DEG_logFC) ~
        FALSE,
      
      OR > 1 &
        DEG_logFC > 0 ~
        TRUE,
      
      OR < 1 &
        DEG_logFC < 0 ~
        TRUE,
      
      TRUE ~
        FALSE
    ),
    
    ## 共定位分析前候选：
    ## 只依据主MR多重校正显著性 + 工具变量强度
    ## DEG及MR-DEG方向一致性在共定位前仅作注释，不作为硬排除条件
    pass_pre_coloc =
      pass_primary_MR &
      pass_instrument_strength,
    
    transcriptomic_concordance_status = case_when(
      is.na(DEG_logFC) |
        is.na(DEG_FDR) ~
        "DEG_result_unavailable",
      
      pass_DEG == FALSE ~
        "does_not_meet_predefined_DEG_definition",
      
      MR_DEG_consistent == TRUE ~
        "MR_DEG_direction_concordant",
      
      TRUE ~
        "MR_DEG_direction_discordant"
    ),
    
    pre_coloc_priority_note = case_when(
      pass_pre_coloc == FALSE ~
        "Not eligible for colocalization based on primary MR/F-statistic",
      
      egger_pleiotropy_flag |
        presso_global_pleiotropy_flag ~
        "Eligible for colocalization, but pleiotropy warning requires cautious interpretation",
      
      sensitivity_incomplete ~
        "Eligible for colocalization; sensitivity evidence is incomplete",
      
      heterogeneity_flag ~
        "Eligible for colocalization; heterogeneity present and random-effects IVW was used",
      
      TRUE ~
        "Eligible for colocalization with no major sensitivity-analysis warning"
    )
  ) %>%
  arrange(
    desc(
      pass_pre_coloc
    ),
    MR_FDR,
    IVW_P
  )

############################################################
## 16. 共定位分析前候选基因
##
## 入选规则：
##   主MR BH-FDR < 0.05
##   + min F > 10
##
## DEG及MR-DEG方向一致性不在此处剔除。
## 下一步对这些基因全部进行colocalization，
## 再结合PP.H3/PP.H4及MR-DEG方向一致性进行最终优先级判断。
############################################################

pre_coloc_candidates <- mr_deg_all %>%
  filter(
    pass_pre_coloc ==
      TRUE
  ) %>%
  arrange(
    MR_FDR,
    IVW_P
  )

############################################################
## 17. 共定位前候选基因：全部MR方法详细统计结果
############################################################

precoloc_gene_list <- unique(
  pre_coloc_candidates$Gene
)

mr_methods_detail <- all_mr %>%
  filter(
    Gene %in% precoloc_gene_list
  ) %>%
  mutate(
    method_id =
      as.character(
        method_id
      ),
    
    method_name =
      if ("method" %in% colnames(all_mr)) {
        as.character(
          method
        )
      } else {
        case_when(
          method_id == "mr_egger_regression" ~
            "MR Egger",
          method_id == "mr_weighted_median" ~
            "Weighted median",
          method_id == "mr_ivw_fe" ~
            "Inverse variance weighted (fixed effects)",
          method_id == "mr_ivw_mre" ~
            "Inverse variance weighted (multiplicative random effects)",
          method_id == "mr_ivw" ~
            "Inverse variance weighted",
          method_id == "mr_simple_mode" ~
            "Simple mode",
          method_id == "mr_weighted_mode" ~
            "Weighted mode",
          method_id == "mr_wald_ratio" ~
            "Wald ratio",
          TRUE ~
            method_id
        )
      },
    
    nsnp_method =
      suppressWarnings(
        as.numeric(
          nsnp
        )
      ),
    
    beta_method =
      suppressWarnings(
        as.numeric(
          b
        )
      ),
    
    se_method =
      suppressWarnings(
        as.numeric(
          se
        )
      ),
    
    p_method =
      suppressWarnings(
        as.numeric(
          pval
        )
      ),
    
    OR_method =
      exp(
        beta_method
      ),
    
    OR_LCI95_method =
      exp(
        beta_method -
          1.96 * se_method
      ),
    
    OR_UCI95_method =
      exp(
        beta_method +
          1.96 * se_method
      ),
    
    method_order = case_when(
      method_id == "mr_egger_regression" ~ 1,
      method_id == "mr_weighted_median" ~ 2,
      method_id %in% c(
        "mr_ivw",
        "mr_ivw_fe",
        "mr_ivw_mre"
      ) ~ 3,
      method_id == "mr_simple_mode" ~ 4,
      method_id == "mr_weighted_mode" ~ 5,
      method_id == "mr_wald_ratio" ~ 6,
      TRUE ~ 99
    )
  ) %>%
  left_join(
    ivw2 %>%
      transmute(
        Gene,
        primary_beta =
          beta,
        primary_OR =
          OR,
        primary_P =
          IVW_P,
        primary_MR_FDR =
          MR_FDR,
        primary_nsnp =
          nsnp,
        primary_method_id =
          if (
            "primary_method_id" %in% colnames(ivw2)
          ) {
            as.character(
              primary_method_id
            )
          } else {
            NA_character_
          }
      ),
    by = "Gene"
  ) %>%
  mutate(
    is_primary_method = case_when(
      !is.na(primary_method_id) ~
        method_id == primary_method_id,
      
      primary_nsnp == 1 ~
        method_id == "mr_wald_ratio",
      
      primary_nsnp >= 2 ~
        method_id %in% c(
          "mr_ivw",
          "mr_ivw_fe",
          "mr_ivw_mre"
        ),
      
      TRUE ~
        FALSE
    )
  ) %>%
  arrange(
    match(
      Gene,
      precoloc_gene_list
    ),
    method_order
  )

mr_methods_count <- mr_methods_detail %>%
  group_by(
    Gene
  ) %>%
  summarise(
    nsnp =
      first(
        primary_nsnp
      ),
    
    n_MR_methods_returned =
      n_distinct(
        method_id
      ),
    
    expected_n_methods = case_when(
      nsnp == 1 ~ 1L,
      nsnp == 2 ~ 1L,
      nsnp >= 3 ~ 5L,
      TRUE ~ NA_integer_
    ),
    
    all_expected_methods_returned = case_when(
      is.na(expected_n_methods) ~
        NA,
      n_MR_methods_returned == expected_n_methods ~
        TRUE,
      TRUE ~
        FALSE
    ),
    
    .groups = "drop"
  )

missing_method_check <- mr_methods_count %>%
  filter(
    nsnp >= 3,
    all_expected_methods_returned == FALSE
  )

if (nrow(missing_method_check) > 0) {
  
  cat("\n")
  cat("注意：以下候选基因nsnp >= 3，但M7实际成功返回的MR方法少于5种。\n")
  cat("程序不会人为补齐缺失方法，请检查原始mr_results_all_methods.tsv。\n")
  
  print(
    missing_method_check,
    n = Inf,
    width = Inf
  )
}

mr_methods_manuscript <- mr_methods_detail %>%
  transmute(
    Gene,
    Method =
      method_name,
    Method_ID =
      method_id,
    SNPs =
      nsnp_method,
    Beta =
      beta_method,
    SE =
      se_method,
    OR =
      OR_method,
    OR_LCI95 =
      OR_LCI95_method,
    OR_UCI95 =
      OR_UCI95_method,
    P_value =
      p_method,
    Primary_method =
      is_primary_method,
    Primary_MR_FDR =
      primary_MR_FDR
  ) %>%
  arrange(
    match(
      Gene,
      precoloc_gene_list
    ),
    match(
      Method_ID,
      c(
        "mr_egger_regression",
        "mr_weighted_median",
        "mr_ivw",
        "mr_ivw_fe",
        "mr_ivw_mre",
        "mr_simple_mode",
        "mr_weighted_mode",
        "mr_wald_ratio"
      )
    )
  )

############################################################
## 18. 生成“硬筛选原因 + 敏感性注释”
############################################################

check_table <- mr_deg_all %>%
  mutate(
    Pre_coloc_keep =
      pass_pre_coloc,
    
    reason_MR_FDR = case_when(
      is.na(MR_FDR) ~
        "MR BH-FDR缺失",
      
      MR_FDR >= MR_FDR_CUTOFF ~
        paste0(
          "主MR经BH校正后未显著（FDR=",
          signif(MR_FDR, 4),
          "，要求<",
          MR_FDR_CUTOFF,
          "）"
        ),
      
      TRUE ~
        NA_character_
    ),
    
    reason_F = case_when(
      is.na(min_F) ~
        "F-statistic结果缺失",
      
      min_F <= MIN_F ~
        paste0(
          "存在弱工具变量（min F=",
          signif(min_F, 4),
          "，要求>",
          MIN_F,
          "）"
        ),
      
      TRUE ~
        NA_character_
    ),
    
    reason_DEG_missing = case_when(
      is.na(DEG_logFC) |
        is.na(DEG_FDR) ~
        "未在Figure2原始DEG结果中匹配到完整结果",
      
      TRUE ~
        NA_character_
    ),
    
    reason_DEG_FDR = case_when(
      !is.na(DEG_FDR) &
        DEG_FDR >= DEG_FDR_CUTOFF ~
        paste0(
          "Figure2 DEG adjusted P/FDR未通过（FDR=",
          signif(DEG_FDR, 4),
          "，要求<",
          DEG_FDR_CUTOFF,
          "）"
        ),
      
      TRUE ~
        NA_character_
    ),
    
    reason_DEG_FC = case_when(
      !is.na(DEG_logFC) &
        abs(DEG_logFC) <= DEG_ABS_LOGFC_CUTOFF ~
        paste0(
          "Figure2差异倍数未通过（logFC=",
          signif(DEG_logFC, 4),
          "，要求|logFC|>",
          DEG_ABS_LOGFC_CUTOFF,
          "）"
        ),
      
      TRUE ~
        NA_character_
    ),
    
    reason_MR_DEG_direction = case_when(
      !is.na(DEG_logFC) &
        !MR_DEG_consistent ~
        paste0(
          "MR与表达方向不一致（OR=",
          signif(OR, 4),
          "，logFC=",
          signif(DEG_logFC, 4),
          "；要求OR>1+上调，或OR<1+下调）"
        ),
      
      TRUE ~
        NA_character_
    ),
    
    sensitivity_heterogeneity_note =
      heterogeneity_status,
    
    sensitivity_egger_note =
      egger_status,
    
    sensitivity_presso_note =
      presso_status_interpretation,
    
    sensitivity_method_direction_note =
      method_direction_note,
    
    sensitivity_loo_note = case_when(
      is.na(loo_direction_consistency_ratio) ~
        "leave-one-out unavailable/not applicable",
      loo_all_same_direction == TRUE ~
        "all leave-one-out estimates retained the primary direction",
      TRUE ~
        paste0(
          "leave-one-out direction consistency=",
          round(
            loo_direction_consistency_ratio * 100,
            1
          ),
          "%"
        )
    )
  )

hard_reason_columns <- c(
  "reason_MR_FDR",
  "reason_F"
)

check_table <- check_table %>%
  rowwise() %>%
  mutate(
    hard_exclusion_reason = {
      
      rr <- c_across(
        all_of(
          hard_reason_columns
        )
      )
      
      rr <- rr[
        !is.na(rr) &
          rr != ""
      ]
      
      if (Pre_coloc_keep) {
        
        "进入共定位分析前候选集；DEG及方向一致性仅作注释"
        
      } else if (length(rr) == 0) {
        
        "未进入候选集，但未识别到硬筛选原因，请人工检查"
        
      } else {
        
        paste(
          rr,
          collapse = "; "
        )
      }
    }
  ) %>%
  ungroup()

############################################################
## 19. 最终审核表
############################################################

check_out <- check_table %>%
  select(
    Gene,
    Pre_coloc_keep,
    hard_exclusion_reason,
    pre_coloc_priority_note,
    
    ## 主MR
    nsnp,
    beta,
    se,
    OR,
    OR_LCI95,
    OR_UCI95,
    IVW_P,
    MR_FDR,
    MR_direction,
    
    ## 工具变量强度
    n_F_SNP,
    mean_F,
    min_F,
    max_F,
    pass_instrument_strength,
    
    ## 异质性
    heterogeneity_P,
    heterogeneity_applicable,
    heterogeneity_flag,
    heterogeneity_status,
    
    ## MR-Egger
    egger_intercept,
    pleiotropy_P,
    pleiotropy_status,
    pleiotropy_applicable,
    egger_pleiotropy_flag,
    egger_status,
    
    ## MR-PRESSO
    MR_PRESSO_status,
    MR_PRESSO_main_P_raw,
    MR_PRESSO_main_P,
    MR_PRESSO_global_P_raw,
    MR_PRESSO_global_P,
    presso_applicable,
    presso_global_pleiotropy_flag,
    presso_status_interpretation,
    
    ## 多方法方向
    n_methods,
    n_same_direction,
    direction_consistency_ratio,
    method_direction_note,
    
    ## leave-one-out
    n_loo,
    n_loo_same_direction,
    loo_direction_consistency_ratio,
    loo_all_same_direction,
    
    ## 敏感性整体注释
    sensitivity_scope,
    sensitivity_warning_count,
    sensitivity_incomplete,
    sensitivity_summary,
    
    ## DEG
    DEG_logFC,
    DEG_P,
    DEG_FDR,
    DEG_direction,
    pass_DEG,
    MR_DEG_consistent,
    transcriptomic_concordance_status,
    
    ## 主判定
    pass_primary_MR,
    pass_pre_coloc,
    
    ## 具体敏感性说明
    sensitivity_heterogeneity_note,
    sensitivity_egger_note,
    sensitivity_presso_note,
    sensitivity_method_direction_note,
    sensitivity_loo_note
  ) %>%
  arrange(
    desc(
      Pre_coloc_keep
    ),
    MR_FDR,
    IVW_P
  )

############################################################
## 20. 输出结果
############################################################

write_tsv(
  mr_evidence_table,
  file.path(
    out_dir,
    "M7_01_all_MR_genes_evidence_table.tsv"
  )
)

write_tsv(
  mr_significant_genes,
  file.path(
    out_dir,
    "M7_02_BH_FDR_significant_MR_genes.tsv"
  )
)

write_tsv(
  mr_deg_all,
  file.path(
    out_dir,
    "M7_03_ALL_MR_genes_with_Fig2_DEG.tsv"
  )
)

write_tsv(
  pre_coloc_candidates,
  file.path(
    out_dir,
    "M7_04_pre_colocalization_candidates.tsv"
  )
)

write_tsv(
  check_out,
  file.path(
    out_dir,
    "M7_05_prioritization_reason_and_sensitivity.tsv"
  )
)

write_tsv(
  inst2 %>%
    arrange(
      min_F
    ),
  file.path(
    out_dir,
    "M7_06_F_statistic_check.tsv"
  )
)

write_tsv(
  loo2,
  file.path(
    out_dir,
    "M7_07_leave_one_out_summary.tsv"
  )
)

write_tsv(
  mr_methods_detail,
  file.path(
    out_dir,
    "M7_08_precoloc_candidates_all_MR_methods.tsv"
  )
)

write_tsv(
  mr_methods_manuscript,
  file.path(
    out_dir,
    "M7_09_precoloc_candidates_MR_methods_for_manuscript.tsv"
  )
)

write_tsv(
  mr_methods_count,
  file.path(
    out_dir,
    "M7_10_precoloc_candidates_MR_method_count_check.tsv"
  )
)

############################################################
## 21. 打印汇总
############################################################

cat("\n")
cat("============================================================\n")
cat("M7.5 SCIENTIFIC PRIORITIZATION SUMMARY\n")
cat("============================================================\n")

cat(
  "进入主MR并进行BH-FDR校正的基因数：",
  nrow(ivw2),
  "\n"
)

cat(
  "主MR nominal P < 0.05：",
  sum(
    ivw2$pass_IVW,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "主MR BH-FDR < 0.05：",
  sum(
    mr_evidence_table$pass_primary_MR,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "min F > 10：",
  sum(
    mr_evidence_table$pass_instrument_strength,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "共定位分析前候选基因数：",
  nrow(pre_coloc_candidates),
  "\n"
)

cat(
  "其中存在heterogeneity提示：",
  sum(
    pre_coloc_candidates$heterogeneity_flag,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "其中存在Egger水平多效性提示：",
  sum(
    pre_coloc_candidates$egger_pleiotropy_flag,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "其中存在MR-PRESSO global pleiotropy提示：",
  sum(
    pre_coloc_candidates$presso_global_pleiotropy_flag,
    na.rm = TRUE
  ),
  "\n"
)

cat("============================================================\n\n")

############################################################
## 22. 打印BH-FDR显著MR基因
############################################################

cat("\n===== BH-FDR显著MR基因 =====\n")

print(
  mr_significant_genes %>%
    select(
      Gene,
      nsnp,
      OR,
      OR_LCI95,
      OR_UCI95,
      IVW_P,
      MR_FDR,
      min_F,
      heterogeneity_P,
      egger_intercept,
      pleiotropy_P,
      MR_PRESSO_global_P,
      direction_consistency_ratio,
      loo_direction_consistency_ratio,
      sensitivity_summary
    ),
  n = Inf,
  width = Inf
)

############################################################
## 23. 打印共定位分析前候选
############################################################

cat("\n===== 共定位分析前候选基因 =====\n")

print(
  pre_coloc_candidates %>%
    select(
      Gene,
      nsnp,
      OR,
      OR_LCI95,
      OR_UCI95,
      IVW_P,
      MR_FDR,
      min_F,
      DEG_logFC,
      DEG_FDR,
      MR_direction,
      DEG_direction,
      MR_DEG_consistent,
      transcriptomic_concordance_status,
      heterogeneity_P,
      pleiotropy_P,
      MR_PRESSO_global_P,
      direction_consistency_ratio,
      loo_direction_consistency_ratio,
      sensitivity_summary,
      pre_coloc_priority_note
    ),
  n = Inf,
  width = Inf
)

############################################################
## 24. 打印未进入共定位前候选的基因
############################################################

cat("\n===== 未进入共定位分析前候选的基因及硬性原因 =====\n")

print(
  check_out %>%
    filter(
      Pre_coloc_keep == FALSE
    ) %>%
    select(
      Gene,
      hard_exclusion_reason,
      sensitivity_summary
    ),
  n = Inf,
  width = Inf
)

############################################################
## 25. 打印共定位前候选基因的全部MR方法
############################################################

cat("\n===== 共定位前候选基因：全部MR方法统计结果 =====\n")

print(
  mr_methods_manuscript,
  n = Inf,
  width = Inf
)

cat("\n===== 每个候选基因实际返回的MR方法数量 =====\n")

print(
  mr_methods_count,
  n = Inf,
  width = Inf
)

############################################################
## 26. 结束提示
############################################################

cat("\n")
cat("============================================================\n")
cat("结果输出目录：\n")
cat(out_dir, "\n")
cat("============================================================\n")

cat("\n输出文件：\n")
cat("1. M7_01_all_MR_genes_evidence_table.tsv\n")
cat("2. M7_02_BH_FDR_significant_MR_genes.tsv\n")
cat("3. M7_03_ALL_MR_genes_with_Fig2_DEG.tsv\n")
cat("4. M7_04_pre_colocalization_candidates.tsv\n")
cat("5. M7_05_prioritization_reason_and_sensitivity.tsv\n")
cat("6. M7_06_F_statistic_check.tsv\n")
cat("7. M7_07_leave_one_out_summary.tsv\n")
cat("8. M7_08_precoloc_candidates_all_MR_methods.tsv\n")
cat("9. M7_09_precoloc_candidates_MR_methods_for_manuscript.tsv\n")
cat("10. M7_10_precoloc_candidates_MR_method_count_check.tsv\n")

cat("\n下一步：对 M7_04_pre_colocalization_candidates.tsv 中的全部候选基因\n")
cat("单独进行 locus-level colocalization 分析。\n")
cat("共定位完成后，再结合PP.H3/PP.H4与MR-DEG方向一致性决定最终优先候选。\n")


## ============================================================
## m8_mr_plot_focus5.R
## Purpose: MR plotting for the final five genes
## Focus genes: AKR1C2, MAPK3, MAP3K11, FADS1, BID
## ============================================================

## 系统报错改为英文
Sys.setenv(LANGUAGE = "en")

## 禁止转化为因子
options(stringsAsFactors = FALSE)

## 清空环境
rm(list = ls())

## ------------------------------------------------------------
## 1. R packages
## ------------------------------------------------------------
required_packages <- c("data.table", "ggplot2", "patchwork", "scales")

missing <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing)) {
  install.packages(missing, repos = "https://cloud.r-project.org")
}

library(data.table)
library(ggplot2)
library(patchwork)
library(scales)

## ------------------------------------------------------------
## 2. Paths
##    只需要根据你的实际文件夹修改 workdir
## ------------------------------------------------------------
workdir <- "C:/Users/92066/Desktop/单细胞孟德尔/MR"

## MR 分析结果所在目录
mr_local_dir <- file.path(workdir, "m7")

## m8 作图输出目录
results_dir <- file.path(workdir, "m8")

## 单个基因图输出目录
individual_dir <- file.path(results_dir, "individual")

## 合并图输出目录
combined_dir <- file.path(results_dir, "combined")

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(individual_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(combined_dir, recursive = TRUE, showWarnings = FALSE)

## ------------------------------------------------------------
## 3. Required input files
## ------------------------------------------------------------
required_files <- c(
  "mr_analysis_ready.tsv",
  "mr_results_all_methods.tsv",
  "mr_results_ivw_wald.tsv",
  "mr_singlesnp.tsv",
  "mr_leaveoneout.tsv"
)

missing_files <- required_files[
  !file.exists(file.path(mr_local_dir, required_files))
]

if (length(missing_files)) {
  stop(sprintf(
    "Missing MR output files in %s: %s. Please check m7 output folder.",
    mr_local_dir,
    paste(missing_files, collapse = ", ")
  ))
}

## ------------------------------------------------------------
## 4. Read data
## ------------------------------------------------------------
analysis_dat <- fread(file.path(mr_local_dir, "mr_analysis_ready.tsv"))
method_dt    <- fread(file.path(mr_local_dir, "mr_results_all_methods.tsv"))
primary_dt   <- fread(file.path(mr_local_dir, "mr_results_ivw_wald.tsv"))
single_dt    <- fread(file.path(mr_local_dir, "mr_singlesnp.tsv"))
loo_dt       <- fread(file.path(mr_local_dir, "mr_leaveoneout.tsv"))

workflow_path <- file.path(mr_local_dir, "mr_gene_selection_workflow.tsv")
workflow_dt <- if (file.exists(workflow_path)) fread(workflow_path) else data.table()

## ------------------------------------------------------------
## 5. Basic column checks
## ------------------------------------------------------------
check_cols <- function(dt, need_cols, dt_name) {
  miss <- setdiff(need_cols, names(dt))
  if (length(miss)) {
    stop(sprintf(
      "%s is missing required columns: %s",
      dt_name,
      paste(miss, collapse = ", ")
    ))
  }
}

check_cols(
  analysis_dat,
  c("gene", "SNP", "beta.exposure", "se.exposure", "beta.outcome", "se.outcome"),
  "mr_analysis_ready.tsv"
)

check_cols(
  method_dt,
  c("gene", "method", "b"),
  "mr_results_all_methods.tsv"
)

check_cols(
  primary_dt,
  c("gene", "method", "b", "or", "or_lci95", "or_uci95", "pval"),
  "mr_results_ivw_wald.tsv"
)

check_cols(
  single_dt,
  c("gene", "SNP", "b", "se"),
  "mr_singlesnp.tsv"
)

check_cols(
  loo_dt,
  c("gene", "SNP", "b", "se"),
  "mr_leaveoneout.tsv"
)

## ------------------------------------------------------------
## 6. Focus genes
## ------------------------------------------------------------
plot_genes <- c("AKR1C2", "MAPK3", "MAP3K11", "FADS1", "BID")

missing_genes <- setdiff(plot_genes, unique(analysis_dat$gene))
if (length(missing_genes)) {
  stop(sprintf(
    "The following focus genes are not present in mr_analysis_ready.tsv: %s",
    paste(missing_genes, collapse = ", ")
  ))
}

message("Plotting genes: ", paste(plot_genes, collapse = ", "))

## ------------------------------------------------------------
## 7. Helper functions
## ------------------------------------------------------------
save_plot_pair <- function(plot_obj, png_path, pdf_path, width, height, dpi = 300) {
  ggsave(
    filename = png_path,
    plot = plot_obj,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )
  ggsave(
    filename = pdf_path,
    plot = plot_obj,
    width = width,
    height = height,
    bg = "white"
  )
}

format_p <- function(x) {
  ifelse(
    is.na(x),
    "NA",
    ifelse(x < 1e-4, formatC(x, format = "e", digits = 2), sprintf("%.4f", x))
  )
}

method_label <- function(x) {
  out <- x
  out <- sub("^Inverse variance weighted \\(fixed effects\\)$", "IVW (FE)", out)
  out <- sub("^Inverse variance weighted \\(multiplicative random effects\\)$", "IVW (MRE)", out)
  out <- sub("^Inverse variance weighted$", "IVW", out)
  out <- sub("^Weighted median$", "Weighted median", out)
  out <- sub("^MR Egger$", "MR-Egger", out)
  out <- sub("^Simple mode$", "Simple mode", out)
  out <- sub("^Weighted mode$", "Weighted mode", out)
  out <- sub("^Wald ratio$", "Wald ratio", out)
  out
}

get_primary_row <- function(gene_name) {
  tmp <- primary_dt[gene == gene_name]
  if (!nrow(tmp)) return(data.table())
  tmp[1]
}

make_empty_plot <- function(gene_name, plot_type, subtitle, detail) {
  ggplot() +
    annotate("text", x = 0, y = 0, label = detail, size = 4.2) +
    labs(
      title = sprintf("%s: %s", gene_name, plot_type),
      subtitle = subtitle
    ) +
    theme_void(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )
}

gene_subtitle <- function(gene_name) {
  primary_row <- get_primary_row(gene_name)
  
  if (!nrow(primary_row)) {
    return("No primary MR result available")
  }
  
  stage_text <- NA_character_
  
  if (nrow(workflow_dt) && all(c("gene", "workflow_stage") %in% names(workflow_dt))) {
    tmp_stage <- workflow_dt[gene == gene_name, workflow_stage]
    if (length(tmp_stage)) stage_text <- tmp_stage[1]
  }
  
  if (is.na(stage_text) || !nzchar(stage_text)) {
    stage_text <- "stage unavailable"
  }
  
  nsnp_text <- if ("nsnp_final" %in% names(primary_row)) {
    as.character(primary_row$nsnp_final)
  } else if ("nsnp" %in% names(primary_row)) {
    as.character(primary_row$nsnp)
  } else {
    as.character(nrow(analysis_dat[gene == gene_name]))
  }
  
  sprintf(
    "%s | nsnp=%s | OR=%.3f | P=%s",
    method_label(primary_row$method),
    nsnp_text,
    primary_row$or,
    format_p(primary_row$pval)
  )
}

scatter_palette <- c(
  "MR-Egger" = "#d95f02",
  "Weighted median" = "#1b9e77",
  "IVW" = "#2c7fb8",
  "IVW (FE)" = "#2c7fb8",
  "IVW (MRE)" = "#2c7fb8",
  "Simple mode" = "#7570b3",
  "Weighted mode" = "#e7298a",
  "Wald ratio" = "#636363"
)

## 只用于 scatter 图图注排序
scatter_method_breaks <- c(
  "MR-Egger",
  "Weighted median",
  "IVW (FE)",
  "IVW (MRE)",
  "IVW",
  "Simple mode",
  "Weighted mode",
  "Wald ratio"
)

## ------------------------------------------------------------
## 8. Plot functions
## ------------------------------------------------------------
build_scatter_plot <- function(gene_name) {
  dat_gene <- analysis_dat[gene == gene_name]
  res_gene <- method_dt[gene == gene_name]
  subtitle <- gene_subtitle(gene_name)
  
  if (!nrow(dat_gene) || !nrow(res_gene)) {
    return(make_empty_plot(gene_name, "Scatter", subtitle, "No harmonized MR data"))
  }
  
  method_order <- c(
    "MR Egger",
    "Weighted median",
    "Inverse variance weighted (fixed effects)",
    "Inverse variance weighted (multiplicative random effects)",
    "Inverse variance weighted",
    "Simple mode",
    "Weighted mode",
    "Wald ratio"
  )
  
  res_gene <- copy(res_gene[
    order(match(method, method_order), na.last = TRUE)
  ])
  
  res_gene[, plot_method := method_label(method)]
  
  x_limits <- range(
    c(
      dat_gene$beta.exposure - dat_gene$se.exposure,
      dat_gene$beta.exposure + dat_gene$se.exposure
    ),
    na.rm = TRUE
  )
  
  if (!all(is.finite(x_limits)) || diff(x_limits) == 0) {
    x_limits <- range(dat_gene$beta.exposure, na.rm = TRUE)
    if (!all(is.finite(x_limits)) || diff(x_limits) == 0) {
      x_limits <- c(-0.1, 0.1)
    }
  }
  
  x_pad <- max(diff(x_limits) * 0.08, 0.02)
  x_vals <- seq(x_limits[1] - x_pad, x_limits[2] + x_pad, length.out = 200)
  
  line_dt <- rbindlist(lapply(seq_len(nrow(res_gene)), function(i) {
    data.table(
      plot_method = res_gene$plot_method[i],
      x = x_vals,
      y = res_gene$b[i] * x_vals
    )
  }))
  
  legend_breaks_use <- intersect(scatter_method_breaks, unique(line_dt$plot_method))
  
  ggplot(dat_gene, aes(x = beta.exposure, y = beta.outcome)) +
    geom_hline(yintercept = 0, color = "grey85", linewidth = 0.4) +
    geom_vline(xintercept = 0, color = "grey85", linewidth = 0.4) +
    geom_segment(
      aes(
        x = beta.exposure - se.exposure,
        xend = beta.exposure + se.exposure,
        y = beta.outcome,
        yend = beta.outcome
      ),
      color = "grey72",
      linewidth = 0.35
    ) +
    geom_segment(
      aes(
        x = beta.exposure,
        xend = beta.exposure,
        y = beta.outcome - se.outcome,
        yend = beta.outcome + se.outcome
      ),
      color = "grey72",
      linewidth = 0.35
    ) +
    geom_line(
      data = line_dt,
      aes(x = x, y = y, color = plot_method),
      linewidth = 0.8,
      alpha = 0.9,
      inherit.aes = FALSE
    ) +
    geom_point(
      shape = 21,
      size = 2.5,
      fill = "#1f3d5b",
      color = "white",
      stroke = 0.4
    ) +
    scale_color_manual(
      values = scatter_palette,
      breaks = legend_breaks_use,
      limits = scatter_method_breaks,
      drop = TRUE
    ) +
    guides(
      color = guide_legend(
        title = "MR method",
        nrow = 1,
        byrow = TRUE,
        title.position = "left"
      )
    ) +
    labs(
      title = gene_name,
      subtitle = subtitle,
      x = "SNP effect on gene expression (beta)",
      y = "SNP effect on outcome risk (beta)",
      color = "MR method"
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8),
      legend.key.width = grid::unit(1.0, "lines"),
      legend.spacing.x = grid::unit(0.15, "cm"),
      plot.margin = margin(6, 12, 10, 6),
      panel.grid.minor = element_blank()
    )
}

build_forest_plot <- function(gene_name) {
  dat_gene <- single_dt[gene == gene_name]
  primary_row <- get_primary_row(gene_name)
  subtitle <- gene_subtitle(gene_name)
  
  if (!nrow(dat_gene) || !nrow(primary_row)) {
    return(make_empty_plot(gene_name, "Forest", subtitle, "No single-SNP MR results"))
  }
  
  plot_dt <- copy(dat_gene[!is.na(b) & !is.na(se) & is.finite(b) & is.finite(se)])
  
  if (!nrow(plot_dt)) {
    return(make_empty_plot(gene_name, "Forest", subtitle, "No valid single-SNP MR results"))
  }
  
  plot_dt[, `:=`(
    estimate = exp(b),
    lci = exp(b - 1.96 * se),
    uci = exp(b + 1.96 * se),
    label = SNP,
    is_summary = FALSE
  )]
  
  plot_dt <- plot_dt[, .(label, estimate, lci, uci, is_summary)]
  
  summary_row <- data.table(
    label = sprintf("Overall (%s)", method_label(primary_row$method)),
    estimate = primary_row$or,
    lci = primary_row$or_lci95,
    uci = primary_row$or_uci95,
    is_summary = TRUE
  )
  
  plot_dt <- rbind(plot_dt, summary_row, fill = TRUE)
  plot_dt[, label := factor(label, levels = rev(unique(label)))]
  
  ggplot(plot_dt, aes(x = estimate, y = label)) +
    geom_vline(xintercept = 1, color = "grey60", linetype = "dashed", linewidth = 0.5) +
    geom_segment(
      aes(x = lci, xend = uci, y = label, yend = label, color = is_summary),
      linewidth = 0.75
    ) +
    geom_point(aes(color = is_summary, shape = is_summary), size = 2.3) +
    scale_color_manual(values = c("FALSE" = "#4d4d4d", "TRUE" = "#b22222"), guide = "none") +
    scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 18), guide = "none") +
    labs(
      title = gene_name,
      subtitle = subtitle,
      x = "Odds ratio per SNP / overall MR estimate",
      y = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      panel.grid.minor = element_blank()
    )
}

build_funnel_plot <- function(gene_name) {
  dat_gene <- single_dt[gene == gene_name]
  primary_row <- get_primary_row(gene_name)
  subtitle <- gene_subtitle(gene_name)
  
  dat_gene <- copy(dat_gene[!is.na(b) & !is.na(se) & is.finite(b) & is.finite(se) & se > 0])
  
  if (nrow(dat_gene) < 2 || !nrow(primary_row)) {
    return(make_empty_plot(gene_name, "Funnel", subtitle, "Need at least 2 valid SNPs"))
  }
  
  dat_gene[, precision := 1 / se]
  
  ggplot(dat_gene, aes(x = b, y = precision)) +
    geom_vline(xintercept = primary_row$b, color = "#b22222", linetype = "dashed", linewidth = 0.6) +
    geom_point(size = 2.4, color = "#1f78b4", alpha = 0.9) +
    labs(
      title = gene_name,
      subtitle = subtitle,
      x = "Single-SNP causal estimate (beta)",
      y = "Precision (1 / SE)"
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      panel.grid.minor = element_blank()
    )
}

build_loo_plot <- function(gene_name) {
  dat_gene <- loo_dt[gene == gene_name]
  primary_row <- get_primary_row(gene_name)
  subtitle <- gene_subtitle(gene_name)
  
  if (!nrow(dat_gene) || !nrow(primary_row)) {
    return(make_empty_plot(gene_name, "Leave-one-out", subtitle, "No leave-one-out results"))
  }
  
  plot_dt <- copy(dat_gene[!is.na(b) & !is.na(se) & is.finite(b) & is.finite(se)])
  
  if (!nrow(plot_dt)) {
    return(make_empty_plot(gene_name, "Leave-one-out", subtitle, "No valid leave-one-out results"))
  }
  
  plot_dt[, `:=`(
    lci = b - 1.96 * se,
    uci = b + 1.96 * se,
    label = sprintf("Leave out %s", SNP)
  )]
  
  plot_dt[, label := factor(label, levels = rev(unique(label)))]
  
  ggplot(plot_dt, aes(x = b, y = label)) +
    geom_vline(xintercept = primary_row$b, color = "#b22222", linetype = "dashed", linewidth = 0.6) +
    geom_segment(
      aes(x = lci, xend = uci, y = label, yend = label),
      color = "#365c8d",
      linewidth = 0.75
    ) +
    geom_point(color = "#365c8d", size = 2.2) +
    labs(
      title = gene_name,
      subtitle = subtitle,
      x = "Leave-one-out MR estimate (beta)",
      y = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      panel.grid.minor = element_blank()
    )
}

## ------------------------------------------------------------
## 9. Build and save plots
## ------------------------------------------------------------
scatter_plots <- list()
forest_plots  <- list()
funnel_plots  <- list()
loo_plots     <- list()
manifest      <- list()

for (gene_name in plot_genes) {
  p_scatter <- build_scatter_plot(gene_name)
  p_forest  <- build_forest_plot(gene_name)
  p_funnel  <- build_funnel_plot(gene_name)
  p_loo     <- build_loo_plot(gene_name)
  
  scatter_plots[[gene_name]] <- p_scatter
  forest_plots[[gene_name]]  <- p_forest
  funnel_plots[[gene_name]]  <- p_funnel
  loo_plots[[gene_name]]     <- p_loo
  
  gene_panel <- (p_scatter | p_forest) / (p_funnel | p_loo) +
    plot_annotation(title = sprintf("%s MR panels", gene_name))
  
  png_path <- file.path(individual_dir, sprintf("%s_mr_panels.png", gene_name))
  pdf_path <- file.path(individual_dir, sprintf("%s_mr_panels.pdf", gene_name))
  
  save_plot_pair(gene_panel, png_path, pdf_path, width = 12, height = 10)
  
  manifest[[length(manifest) + 1L]] <- data.table(
    gene = gene_name,
    plot_type = "per_gene_panel",
    png_path = png_path,
    pdf_path = pdf_path
  )
}

## ------------------------------------------------------------
## 10. Combined plots
##     这里只改 scatter_combined，其他图不动
## ------------------------------------------------------------

scatter_combined <- (
  wrap_plots(scatter_plots, ncol = 2, guides = "collect") +
    plot_annotation(title = "MR scatter plots for five focus genes")
) &
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    legend.key.width = grid::unit(1.0, "lines"),
    legend.spacing.x = grid::unit(0.15, "cm"),
    plot.margin = margin(8, 18, 12, 8)
  )

forest_combined <- wrap_plots(forest_plots, ncol = 2) +
  plot_annotation(title = "MR forest plots for five focus genes")

funnel_combined <- wrap_plots(funnel_plots, ncol = 2) +
  plot_annotation(title = "MR funnel plots for five focus genes")

loo_combined <- wrap_plots(loo_plots, ncol = 2) +
  plot_annotation(title = "MR leave-one-out plots for five focus genes")

## scatter 合并图加宽，避免图注被截断
save_plot_pair(
  scatter_combined,
  file.path(combined_dir, "focus5_genes_scatter.png"),
  file.path(combined_dir, "focus5_genes_scatter.pdf"),
  width = 15,
  height = 10.5,
  dpi = 600
)

save_plot_pair(
  forest_combined,
  file.path(combined_dir, "focus5_genes_forest.png"),
  file.path(combined_dir, "focus5_genes_forest.pdf"),
  width = 12,
  height = 10
)

save_plot_pair(
  funnel_combined,
  file.path(combined_dir, "focus5_genes_funnel.png"),
  file.path(combined_dir, "focus5_genes_funnel.pdf"),
  width = 12,
  height = 10
)

save_plot_pair(
  loo_combined,
  file.path(combined_dir, "focus5_genes_leaveoneout.png"),
  file.path(combined_dir, "focus5_genes_leaveoneout.pdf"),
  width = 12,
  height = 10
)

manifest[[length(manifest) + 1L]] <- data.table(
  gene = "ALL",
  plot_type = c(
    "scatter_combined",
    "forest_combined",
    "funnel_combined",
    "leaveoneout_combined"
  ),
  png_path = c(
    file.path(combined_dir, "focus5_genes_scatter.png"),
    file.path(combined_dir, "focus5_genes_forest.png"),
    file.path(combined_dir, "focus5_genes_funnel.png"),
    file.path(combined_dir, "focus5_genes_leaveoneout.png")
  ),
  pdf_path = c(
    file.path(combined_dir, "focus5_genes_scatter.pdf"),
    file.path(combined_dir, "focus5_genes_forest.pdf"),
    file.path(combined_dir, "focus5_genes_funnel.pdf"),
    file.path(combined_dir, "focus5_genes_leaveoneout.pdf")
  )
)

manifest_dt <- rbindlist(manifest, fill = TRUE)
fwrite(manifest_dt, file.path(results_dir, "plot_manifest.tsv"), sep = "\t")

message("Finished MR plotting. Outputs saved in: ", normalizePath(results_dir, winslash = "/", mustWork = FALSE))