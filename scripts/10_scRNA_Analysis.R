############################################################
## GSE294300 CRC scRNA-seq
## Step 1
##
## 36 samples
## 18 patients
## Paired Normal / Tumor
##
## 目标：
## 1. 读取36个raw 10X矩阵
## 2. 根据cell_batch仅保留214854个目标细胞
## 3. 建立patient-level metadata
## 4. 创建Seurat对象
## 5. 计算QC指标
## 6. 暂不进行QC过滤
##
## Candidate genes:
## AKR1C2 / MAPK3 / MAP3K11
############################################################


rm(list = ls())

options(stringsAsFactors = FALSE)
Sys.setenv(LANGUAGE = "en")

set.seed(12345)


############################################################
## 1. 加载R包
############################################################

library(Seurat)
library(Matrix)
library(data.table)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(patchwork)


############################################################
## 2. 路径
############################################################

data_dir <- "C:/Users/92066/Desktop/单细胞孟德尔/Figure7/GSE294300"

extract_dir <- file.path(
  data_dir,
  "GSE294300_RAW"
)

batch_file <- file.path(
  data_dir,
  "GSE294300_cell_batch.tsv.gz"
)


if (!dir.exists(extract_dir)) {
  stop("没有找到已经解压的GSE294300_RAW目录。")
}

if (!file.exists(batch_file)) {
  stop("没有找到GSE294300_cell_batch.tsv.gz。")
}


setwd(data_dir)


############################################################
## 3. 固定三个候选基因
############################################################

target_genes <- c(
  "AKR1C2",
  "MAPK3",
  "MAP3K11"
)


############################################################
## 4. batch名称标准化函数
############################################################

normalize_batch <- function(x) {
  
  x <- tolower(
    as.character(x)
  )
  
  type <- ifelse(
    grepl("nor", x),
    "nor",
    ifelse(
      grepl("tum", x),
      "tum",
      NA_character_
    )
  )
  
  num <- suppressWarnings(
    as.integer(
      sub(
        ".*?([0-9]+)$",
        "\\1",
        x
      )
    )
  )
  
  ifelse(
    is.na(type) | is.na(num),
    NA_character_,
    sprintf(
      "%s_%02d",
      type,
      num
    )
  )
}


############################################################
## 5. 读取cell_batch
############################################################

cell_batch <- data.table::fread(
  batch_file,
  data.table = FALSE,
  check.names = FALSE
)


if (!all(
  c("Cell", "batch") %in%
  colnames(cell_batch)
)) {
  
  stop(
    "cell_batch必须包含Cell和batch两列。"
  )
}


############################################################
## 6. 建立正式metadata基础信息
############################################################

cell_batch$batch_norm <- normalize_batch(
  cell_batch$batch
)


cell_batch$group <- ifelse(
  grepl(
    "^nor_",
    cell_batch$batch_norm
  ),
  "Normal",
  "Tumor"
)


cell_batch$patient_number <- as.integer(
  stringr::str_extract(
    cell_batch$batch_norm,
    "[0-9]+$"
  )
)


cell_batch$patient <- sprintf(
  "P%02d",
  cell_batch$patient_number
)


## patient就是配对单位
cell_batch$pair_id <- cell_batch$patient


## 例如：
## P01_Normal
## P01_Tumor

cell_batch$sample_id <- paste(
  cell_batch$patient,
  cell_batch$group,
  sep = "_"
)


############################################################
## 7. 检查基本结构
############################################################

cat("\n===============================\n")
cat("cell_batch基本信息\n")
cat("===============================\n")

cat(
  "Cells:",
  nrow(cell_batch),
  "\n"
)

cat(
  "Samples:",
  length(unique(cell_batch$batch_norm)),
  "\n"
)

cat(
  "Patients:",
  length(unique(cell_batch$patient)),
  "\n"
)

cat("\nGroups:\n")

print(
  table(
    unique(
      cell_batch[
        c(
          "batch_norm",
          "group"
        )
      ]
    )$group
  )
)


############################################################
## 8. 找到36套10X文件
############################################################

all_files <- list.files(
  extract_dir,
  recursive = TRUE,
  full.names = TRUE
)


barcode_files <- all_files[
  grepl(
    "_barcodes\\.tsv\\.gz$",
    basename(all_files)
  )
]


feature_files <- all_files[
  grepl(
    "_features\\.tsv\\.gz$",
    basename(all_files)
  )
]


matrix_files <- all_files[
  grepl(
    "_matrix\\.mtx\\.gz$",
    basename(all_files)
  )
]


if (
  length(barcode_files) != 36 ||
  length(feature_files) != 36 ||
  length(matrix_files) != 36
) {
  
  stop(
    "10X文件数量不是36套。"
  )
}


############################################################
## 9. 从文件名识别batch
############################################################

get_batch_from_file <- function(x) {
  
  x <- basename(x)
  
  x <- sub(
    "_barcodes\\.tsv\\.gz$",
    "",
    x
  )
  
  x <- sub(
    "_features\\.tsv\\.gz$",
    "",
    x
  )
  
  x <- sub(
    "_matrix\\.mtx\\.gz$",
    "",
    x
  )
  
  normalize_batch(x)
}


barcode_manifest <- data.frame(
  batch = sapply(
    barcode_files,
    get_batch_from_file
  ),
  barcode_file = barcode_files,
  stringsAsFactors = FALSE
)


feature_manifest <- data.frame(
  batch = sapply(
    feature_files,
    get_batch_from_file
  ),
  feature_file = feature_files,
  stringsAsFactors = FALSE
)


matrix_manifest <- data.frame(
  batch = sapply(
    matrix_files,
    get_batch_from_file
  ),
  matrix_file = matrix_files,
  stringsAsFactors = FALSE
)


############################################################
## 10. 合并成正式文件manifest
############################################################

file_manifest <- barcode_manifest %>%
  
  dplyr::left_join(
    feature_manifest,
    by = "batch"
  ) %>%
  
  dplyr::left_join(
    matrix_manifest,
    by = "batch"
  )


## 加入patient / group
sample_info <- cell_batch %>%
  
  dplyr::select(
    batch_norm,
    patient_number,
    patient,
    pair_id,
    group,
    sample_id
  ) %>%
  
  dplyr::distinct()


file_manifest <- file_manifest %>%
  
  dplyr::left_join(
    sample_info,
    by = c(
      "batch" = "batch_norm"
    )
  ) %>%
  
  dplyr::arrange(
    patient_number,
    group
  )


if (nrow(file_manifest) != 36) {
  stop("file_manifest不是36个样本。")
}


cat("\n===============================\n")
cat("36个样本\n")
cat("===============================\n")

print(
  file_manifest[
    c(
      "batch",
      "patient",
      "group",
      "sample_id"
    )
  ],
  row.names = FALSE
)


############################################################
## 11. 保存feature ID / gene symbol映射
############################################################

first_feature_file <- file_manifest$feature_file[1]


feature_reference <- data.table::fread(
  first_feature_file,
  header = FALSE,
  data.table = FALSE,
  showProgress = FALSE
)


if (ncol(feature_reference) < 2) {
  
  stop(
    "features.tsv.gz没有gene symbol列。"
  )
}


colnames(feature_reference)[1:2] <- c(
  "ensembl_id",
  "gene_symbol"
)


if (ncol(feature_reference) >= 3) {
  
  colnames(feature_reference)[3] <-
    "feature_type"
}


feature_reference$gene_symbol <-
  as.character(
    feature_reference$gene_symbol
  )


## Seurat ReadMtx(unique.features=TRUE)
## 会将重复symbol转换为唯一名称

feature_reference$seurat_feature <-
  make.unique(
    feature_reference$gene_symbol
  )


write.csv(
  feature_reference,
  file = "GSE294300_feature_reference.csv",
  row.names = FALSE
)


############################################################
## 12. 检查三个目标基因
############################################################

target_gene_feature_check <-
  data.frame(
    gene = target_genes,
    n_features = sapply(
      target_genes,
      function(g) {
        sum(
          feature_reference$gene_symbol == g
        )
      }
    )
  )


cat("\n===============================\n")
cat("三个候选基因feature检查\n")
cat("===============================\n")

print(
  target_gene_feature_check,
  row.names = FALSE
)


write.csv(
  target_gene_feature_check,
  file = "GSE294300_target_gene_feature_check.csv",
  row.names = FALSE
)


if (
  any(
    target_gene_feature_check$n_features == 0
  )
) {
  
  warning(
    "至少一个候选基因没有在features文件中找到。"
  )
}


############################################################
## 13. 逐样本读取10X原始矩阵
##
## 关键：
## 每次只读取一个样本
## 然后立即只保留cell_batch指定细胞
##
## 不会把5700万个barcode同时加载
############################################################

scRNA <- NULL


sample_processing_log <- list()


for (
  i in seq_len(
    nrow(file_manifest)
  )
) {
  
  ##########################################################
  ## 当前样本
  ##########################################################
  
  bn <- file_manifest$batch[i]
  
  pat <- file_manifest$patient[i]
  
  grp <- file_manifest$group[i]
  
  sid <- file_manifest$sample_id[i]
  
  
  cat("\n\n")
  cat("====================================================\n")
  
  cat(
    "Processing sample ",
    i,
    "/",
    nrow(file_manifest),
    ": ",
    bn,
    "  ",
    pat,
    "  ",
    grp,
    "\n",
    sep = ""
  )
  
  cat("====================================================\n")
  
  
  ##########################################################
  ## 13.1 读取当前样本10X矩阵
  ##########################################################
  
  mat_raw <- Seurat::ReadMtx(
    
    mtx =
      file_manifest$matrix_file[i],
    
    cells =
      file_manifest$barcode_file[i],
    
    features =
      file_manifest$feature_file[i],
    
    cell.column = 1,
    
    feature.column = 2,
    
    unique.features = TRUE
  )
  
  
  cat(
    "Raw matrix:",
    nrow(mat_raw),
    "genes x",
    ncol(mat_raw),
    "barcodes\n"
  )
  
  
  ##########################################################
  ## 13.2 找到该样本目标细胞
  ##########################################################
  
  meta_sample <- cell_batch[
    cell_batch$batch_norm == bn,
    ,
    drop = FALSE
  ]
  
  
  target_barcode_original <-
    as.character(
      meta_sample$Cell
    )
  
  
  ## 去掉cell_batch中的 -1 / -2 / ...
  target_core <- sub(
    "-[0-9]+$",
    "",
    target_barcode_original
  )
  
  
  ## 原始独立10X matrix barcode
  raw_barcode <- colnames(
    mat_raw
  )
  
  
  raw_core <- sub(
    "-[0-9]+$",
    "",
    raw_barcode
  )
  
  
  ##########################################################
  ## 13.3 barcode匹配
  ##########################################################
  
  match_index <- match(
    target_core,
    raw_core
  )
  
  
  if (
    any(
      is.na(match_index)
    )
  ) {
    
    stop(
      "Sample ",
      bn,
      " 存在无法匹配的target cells。"
    )
  }
  
  
  ##########################################################
  ## 13.4 仅保留目标细胞
  ##########################################################
  
  mat_keep <- mat_raw[
    ,
    match_index,
    drop = FALSE
  ]
  
  
  retained_raw_barcode <-
    raw_barcode[
      match_index
    ]
  
  
  cat(
    "Retained target cells:",
    ncol(mat_keep),
    "\n"
  )
  
  
  ##########################################################
  ## 13.5 建立全局唯一细胞名
  ##
  ## 例如：
  ## nor_01_AAACCCAAGACCTGGA
  ##########################################################
  
  global_cell_name <- paste(
    bn,
    target_core,
    sep = "_"
  )
  
  
  if (
    any(
      duplicated(
        global_cell_name
      )
    )
  ) {
    
    stop(
      "当前样本产生重复global cell name：",
      bn
    )
  }
  
  
  colnames(mat_keep) <-
    global_cell_name
  
  
  ##########################################################
  ## 13.6 建立每个细胞的metadata
  ##########################################################
  
  cell_meta <- data.frame(
    
    orig.ident = bn,
    
    sample = bn,
    
    sample_id = sid,
    
    patient = pat,
    
    pair_id = pat,
    
    group = grp,
    
    barcode_core =
      target_core,
    
    cell_batch_barcode =
      target_barcode_original,
    
    raw_barcode =
      retained_raw_barcode,
    
    stringsAsFactors = FALSE
  )
  
  
  rownames(cell_meta) <-
    global_cell_name
  
  
  ##########################################################
  ## 13.7 创建当前样本Seurat对象
  ##########################################################
  
  obj_sample <- CreateSeuratObject(
    
    counts = mat_keep,
    
    meta.data = cell_meta,
    
    project = bn,
    
    min.cells = 0,
    
    min.features = 0
  )
  
  
  ##########################################################
  ## 13.8 计算线粒体比例
  ##########################################################
  
  obj_sample[["percent.mt"]] <-
    PercentageFeatureSet(
      obj_sample,
      pattern = "^MT-"
    )
  
  
  ##########################################################
  ## 13.9 记录样本处理信息
  ##########################################################
  
  sample_processing_log[[i]] <-
    data.frame(
      
      sample = bn,
      
      patient = pat,
      
      group = grp,
      
      raw_barcodes =
        ncol(mat_raw),
      
      retained_cells =
        ncol(obj_sample),
      
      genes =
        nrow(obj_sample),
      
      stringsAsFactors = FALSE
    )
  
  
  ##########################################################
  ## 13.10 合并
  ##
  ## 不JoinLayers
  ## 先保留Seurat v5各sample counts layer
  ##########################################################
  
  if (is.null(scRNA)) {
    
    scRNA <- obj_sample
    
  } else {
    
    scRNA <- merge(
      
      x = scRNA,
      
      y = obj_sample,
      
      merge.data = FALSE
    )
  }
  
  
  ##########################################################
  ## 13.11 释放当前raw大矩阵
  ##########################################################
  
  rm(
    mat_raw,
    mat_keep,
    obj_sample,
    raw_barcode,
    raw_core,
    match_index
  )
  
  gc()
  
  
  cat(
    "Current merged cells:",
    ncol(scRNA),
    "\n"
  )
}


############################################################
## 14. 汇总读取日志
############################################################

sample_processing_log <-
  dplyr::bind_rows(
    sample_processing_log
  )


write.csv(
  sample_processing_log,
  file =
    "GSE294300_sample_processing_log.csv",
  row.names = FALSE
)


############################################################
## 15. 整理最终metadata
############################################################

scRNA$group <- factor(
  scRNA$group,
  levels = c(
    "Normal",
    "Tumor"
  )
)


## patient顺序
patient_levels <- sprintf(
  "P%02d",
  sort(
    unique(
      cell_batch$patient_number
    )
  )
)


scRNA$patient <- factor(
  scRNA$patient,
  levels = patient_levels
)


scRNA$pair_id <- factor(
  scRNA$pair_id,
  levels = patient_levels
)


############################################################
## 16. 最终细胞数检查
############################################################

cat("\n\n")
cat("############################################################\n")
cat("## Seurat object created\n")
cat("############################################################\n")

cat(
  "\nGenes:",
  nrow(scRNA)
)

cat(
  "\nCells:",
  ncol(scRNA)
)

cat(
  "\nPatients:",
  length(
    unique(
      scRNA$patient
    )
  )
)

cat(
  "\nSamples:",
  length(
    unique(
      scRNA$sample
    )
  )
)

cat("\n\nGroup cell numbers:\n")

print(
  table(
    scRNA$group
  )
)


if (
  ncol(scRNA) != nrow(cell_batch)
) {
  
  stop(
    "最终Seurat细胞数与cell_batch不一致。"
  )
}


############################################################
## 17. 检查三基因是否存在
############################################################

cat("\n===============================\n")
cat("Candidate genes in Seurat object\n")
cat("===============================\n")


for (g in target_genes) {
  
  gene_hits <- grep(
    paste0(
      "^",
      g,
      "(\\.[0-9]+)?$"
    ),
    rownames(scRNA),
    value = TRUE
  )
  
  cat(
    "\n",
    g,
    ": ",
    paste(
      gene_hits,
      collapse = ", "
    ),
    "\n",
    sep = ""
  )
}


############################################################
## 18. QC metadata
############################################################

qc_meta <- scRNA@meta.data


############################################################
## 19. QC总体quantile
############################################################

qc_quantiles <- data.frame(
  
  metric = c(
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mt"
  ),
  
  Min = c(
    min(qc_meta$nFeature_RNA),
    min(qc_meta$nCount_RNA),
    min(qc_meta$percent.mt)
  ),
  
  Q01 = c(
    quantile(
      qc_meta$nFeature_RNA,
      0.01
    ),
    quantile(
      qc_meta$nCount_RNA,
      0.01
    ),
    quantile(
      qc_meta$percent.mt,
      0.01
    )
  ),
  
  Q05 = c(
    quantile(
      qc_meta$nFeature_RNA,
      0.05
    ),
    quantile(
      qc_meta$nCount_RNA,
      0.05
    ),
    quantile(
      qc_meta$percent.mt,
      0.05
    )
  ),
  
  Median = c(
    median(
      qc_meta$nFeature_RNA
    ),
    median(
      qc_meta$nCount_RNA
    ),
    median(
      qc_meta$percent.mt
    )
  ),
  
  Q95 = c(
    quantile(
      qc_meta$nFeature_RNA,
      0.95
    ),
    quantile(
      qc_meta$nCount_RNA,
      0.95
    ),
    quantile(
      qc_meta$percent.mt,
      0.95
    )
  ),
  
  Q99 = c(
    quantile(
      qc_meta$nFeature_RNA,
      0.99
    ),
    quantile(
      qc_meta$nCount_RNA,
      0.99
    ),
    quantile(
      qc_meta$percent.mt,
      0.99
    )
  ),
  
  Max = c(
    max(qc_meta$nFeature_RNA),
    max(qc_meta$nCount_RNA),
    max(qc_meta$percent.mt)
  )
)


cat("\n\n")
cat("############################################################\n")
cat("## QC quantiles\n")
cat("############################################################\n")

print(
  qc_quantiles,
  row.names = FALSE
)


write.csv(
  qc_quantiles,
  file =
    "GSE294300_QC_quantiles_before_filter.csv",
  row.names = FALSE
)


############################################################
## 20. 每个sample QC统计
############################################################

qc_by_sample <- qc_meta %>%
  
  dplyr::group_by(
    sample,
    patient,
    group
  ) %>%
  
  dplyr::summarise(
    
    n_cells = dplyr::n(),
    
    median_nFeature =
      median(nFeature_RNA),
    
    median_nCount =
      median(nCount_RNA),
    
    median_percent_mt =
      median(percent.mt),
    
    Q95_nFeature =
      quantile(
        nFeature_RNA,
        0.95
      ),
    
    Q95_nCount =
      quantile(
        nCount_RNA,
        0.95
      ),
    
    Q95_percent_mt =
      quantile(
        percent.mt,
        0.95
      ),
    
    .groups = "drop"
  )


write.csv(
  qc_by_sample,
  file =
    "GSE294300_QC_by_sample_before_filter.csv",
  row.names = FALSE
)


############################################################
## 21. QC绘图数据
############################################################

qc_plot_df <- qc_meta %>%
  
  dplyr::select(
    sample,
    patient,
    group,
    nFeature_RNA,
    nCount_RNA,
    percent.mt
  ) %>%
  
  tidyr::pivot_longer(
    
    cols = c(
      nFeature_RNA,
      nCount_RNA,
      percent.mt
    ),
    
    names_to = "metric",
    
    values_to = "value"
  )


qc_plot_df$metric <- factor(
  qc_plot_df$metric,
  levels = c(
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mt"
  )
)


############################################################
## 22. Figure QC-1
## Normal vs Tumor总体QC
##
## 注意：
## 这里只是QC可视化
## 不做Normal/Tumor统计推断
############################################################

p_qc_group <- ggplot(
  
  qc_plot_df,
  
  aes(
    x = group,
    y = value,
    fill = group
  )
  
) +
  
  geom_violin(
    scale = "width",
    trim = TRUE,
    linewidth = 0.25
  ) +
  
  geom_boxplot(
    width = 0.10,
    outlier.shape = NA,
    fill = "white",
    linewidth = 0.3
  ) +
  
  facet_wrap(
    ~ metric,
    scales = "free_y",
    nrow = 1
  ) +
  
  theme_classic(
    base_size = 13
  ) +
  
  labs(
    x = NULL,
    y = NULL,
    fill = NULL
  ) +
  
  theme(
    strip.background =
      element_blank(),
    
    strip.text =
      element_text(
        size = 14,
        face = "bold"
      ),
    
    axis.text =
      element_text(
        color = "black"
      ),
    
    legend.position =
      "none",
    
    panel.grid =
      element_blank()
  )


ggsave(
  filename =
    "GSE294300_QC_before_filter_by_group.pdf",
  plot = p_qc_group,
  width = 9,
  height = 4.8
)


ggsave(
  filename =
    "GSE294300_QC_before_filter_by_group.png",
  plot = p_qc_group,
  width = 9,
  height = 4.8,
  dpi = 300
)


############################################################
## 23. Figure QC-2
## 36个样本QC
############################################################

sample_order <- file_manifest$batch


qc_plot_df$sample <- factor(
  qc_plot_df$sample,
  levels = sample_order
)


p_qc_sample <- ggplot(
  
  qc_plot_df,
  
  aes(
    x = sample,
    y = value,
    fill = group
  )
  
) +
  
  geom_boxplot(
    outlier.shape = NA,
    linewidth = 0.2,
    width = 0.75
  ) +
  
  facet_wrap(
    ~ metric,
    scales = "free_y",
    nrow = 1
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  labs(
    x = NULL,
    y = NULL,
    fill = NULL
  ) +
  
  theme(
    
    strip.background =
      element_blank(),
    
    strip.text =
      element_text(
        size = 14,
        face = "bold"
      ),
    
    axis.text.x =
      element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5,
        size = 6,
        color = "black"
      ),
    
    axis.text.y =
      element_text(
        color = "black"
      ),
    
    legend.position =
      "right",
    
    panel.grid =
      element_blank()
  )


ggsave(
  filename =
    "GSE294300_QC_before_filter_by_sample.pdf",
  plot = p_qc_sample,
  width = 15,
  height = 5.5,
  limitsize = FALSE
)


ggsave(
  filename =
    "GSE294300_QC_before_filter_by_sample.png",
  plot = p_qc_sample,
  width = 15,
  height = 5.5,
  dpi = 300,
  limitsize = FALSE
)


############################################################
## 24. Figure QC-3
## nCount vs nFeature
##
## 只随机抽30000细胞绘图
## 防止21万点导致PDF巨大
############################################################

set.seed(12345)


plot_cells <- sample(
  rownames(qc_meta),
  size = min(
    30000,
    nrow(qc_meta)
  )
)


qc_scatter_df <-
  qc_meta[
    plot_cells,
    ,
    drop = FALSE
  ]


p_qc_scatter <- ggplot(
  
  qc_scatter_df,
  
  aes(
    x = nCount_RNA,
    y = nFeature_RNA,
    color = percent.mt
  )
  
) +
  
  geom_point(
    size = 0.35,
    alpha = 0.5
  ) +
  
  theme_classic(
    base_size = 13
  ) +
  
  labs(
    x = "nCount_RNA",
    y = "nFeature_RNA",
    color = "Mitochondrial\npercentage"
  ) +
  
  theme(
    axis.text =
      element_text(
        color = "black"
      ),
    
    panel.grid =
      element_blank()
  )


ggsave(
  filename =
    "GSE294300_QC_nCount_nFeature_scatter.pdf",
  plot = p_qc_scatter,
  width = 6.5,
  height = 5.5
)


ggsave(
  filename =
    "GSE294300_QC_nCount_nFeature_scatter.png",
  plot = p_qc_scatter,
  width = 6.5,
  height = 5.5,
  dpi = 300
)


############################################################
## 25. 保存QC前metadata
############################################################

write.csv(
  scRNA@meta.data,
  file =
    "GSE294300_metadata_before_QC.csv",
  row.names = TRUE
)


############################################################
## 26. 保存QC前Seurat对象
##
## 注意：
## 此时尚未：
## - QC filtering
## - doublet removal
## - Normalize
## - Harmony
## - clustering
############################################################

saveRDS(
  
  scRNA,
  
  file =
    "GSE294300_scRNA_214854cells_before_QC.rds",
  
  compress = FALSE
)


############################################################
## 27. 最终一次性输出
############################################################

cat("\n\n")
cat("############################################################\n")
cat("## GSE294300 Step 1 completed\n")
cat("############################################################\n")


cat(
  "\nGenes:",
  nrow(scRNA)
)


cat(
  "\nCells:",
  ncol(scRNA)
)


cat(
  "\nSamples:",
  length(
    unique(
      scRNA$sample
    )
  )
)


cat(
  "\nPatients:",
  length(
    unique(
      scRNA$patient
    )
  )
)


cat("\n\nCells by group:\n")

print(
  table(
    scRNA$group
  )
)


cat("\n\nQC quantiles:\n")

print(
  qc_quantiles,
  row.names = FALSE
)


cat("\n\nCandidate genes:\n")

for (g in target_genes) {
  
  hits <- grep(
    paste0(
      "^",
      g,
      "(\\.[0-9]+)?$"
    ),
    rownames(scRNA),
    value = TRUE
  )
  
  cat(
    g,
    " -> ",
    paste(
      hits,
      collapse = ", "
    ),
    "\n",
    sep = ""
  )
}


cat("\n")
cat(
  "Saved: GSE294300_scRNA_214854cells_before_QC.rds\n"
)
############################################################
## GSE294300 CRC scRNA-seq
## Step 2
##
## QC filtering + doublet removal
##
## 18 patients
## 18 paired Normal/Tumor
##
## Candidate genes:
## AKR1C2 / MAPK3 / MAP3K11
############################################################

rm(list = ls())
options(stringsAsFactors = FALSE)
Sys.setenv(LANGUAGE = "en")

set.seed(12345)

############################################################
## 1. 加载包
############################################################

library(Seurat)
library(Matrix)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

############################################################
## scDblFinder需要Bioconductor包
############################################################

if (!requireNamespace("scDblFinder", quietly = TRUE)) {
  stop(
    "未安装scDblFinder。\n",
    "请先运行：\n",
    "if (!requireNamespace('BiocManager')) install.packages('BiocManager')\n",
    "BiocManager::install('scDblFinder')"
  )
}

if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
  stop("未安装SingleCellExperiment。")
}

if (!requireNamespace("BiocParallel", quietly = TRUE)) {
  stop("未安装BiocParallel。")
}


############################################################
## 2. 路径
############################################################

data_dir <- "C:/Users/92066/Desktop/单细胞孟德尔/Figure7/GSE294300"

setwd(data_dir)

input_file <- file.path(
  data_dir,
  "GSE294300_scRNA_214854cells_before_QC.rds"
)

if (!file.exists(input_file)) {
  stop("没有找到QC前Seurat对象。")
}


############################################################
## 3. 读取Seurat对象
############################################################

scRNA <- readRDS(
  input_file
)

cat("\n===============================\n")
cat("Initial Seurat object\n")
cat("===============================\n")

cat("Genes    :", nrow(scRNA), "\n")
cat("Cells    :", ncol(scRNA), "\n")
cat("Samples  :", length(unique(scRNA$sample)), "\n")
cat("Patients :", length(unique(scRNA$patient)), "\n")

print(
  table(scRNA$group)
)


############################################################
## 4. 固定三个候选基因
############################################################

target_genes <- c(
  "AKR1C2",
  "MAPK3",
  "MAP3K11"
)


############################################################
## 5. QC阈值
##
## GSE294300原研究：
## genes: 200–6000
## mitochondrial proportion <=20%
##
## 当前cell_batch中最低nFeature已经621，
## 因此前面的200基本不会额外删除细胞
############################################################

min_features <- 200
max_features <- 6000
max_percent_mt <- 20


############################################################
## 6. QC前metadata
############################################################

meta_before <- scRNA@meta.data


############################################################
## 7. 给每个细胞标记QC状态
############################################################

scRNA$QC_low_feature <-
  scRNA$nFeature_RNA < min_features

scRNA$QC_high_feature <-
  scRNA$nFeature_RNA > max_features

scRNA$QC_high_mt <-
  scRNA$percent.mt > max_percent_mt


scRNA$QC_keep <-
  !scRNA$QC_low_feature &
  !scRNA$QC_high_feature &
  !scRNA$QC_high_mt


############################################################
## 8. 查看每种QC原因
############################################################

cat("\n\n")
cat("############################################################\n")
cat("## QC filtering summary\n")
cat("############################################################\n")

cat(
  "\nInitial cells:",
  ncol(scRNA)
)

cat(
  "\nLow nFeature (<200):",
  sum(scRNA$QC_low_feature)
)

cat(
  "\nHigh nFeature (>6000):",
  sum(scRNA$QC_high_feature)
)

cat(
  "\nHigh mitochondrial (>20%):",
  sum(scRNA$QC_high_mt)
)

cat(
  "\nCells passing all QC:",
  sum(scRNA$QC_keep)
)

cat(
  "\nCells removed by QC:",
  sum(!scRNA$QC_keep)
)

cat("\n")


############################################################
## 9. QC前后每个样本细胞数
############################################################

qc_sample_summary <- scRNA@meta.data %>%
  
  dplyr::group_by(
    sample,
    patient,
    group
  ) %>%
  
  dplyr::summarise(
    
    cells_before = dplyr::n(),
    
    cells_pass_QC =
      sum(QC_keep),
    
    cells_removed =
      sum(!QC_keep),
    
    retention_rate =
      cells_pass_QC / cells_before,
    
    .groups = "drop"
  )


write.csv(
  qc_sample_summary,
  file =
    "GSE294300_QC_retention_by_sample.csv",
  row.names = FALSE
)


############################################################
## 10. QC过滤
############################################################

scRNA <- subset(
  scRNA,
  subset = QC_keep
)


cat("\n===============================\n")
cat("After QC filtering\n")
cat("===============================\n")

cat(
  "Cells:",
  ncol(scRNA),
  "\n"
)

print(
  table(scRNA$group)
)


############################################################
## 11. QC过滤后重新统计分位数
############################################################

meta_qc <- scRNA@meta.data


qc_after_quantiles <- data.frame(
  
  metric = c(
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mt"
  ),
  
  Min = c(
    min(meta_qc$nFeature_RNA),
    min(meta_qc$nCount_RNA),
    min(meta_qc$percent.mt)
  ),
  
  Q01 = c(
    quantile(meta_qc$nFeature_RNA, 0.01),
    quantile(meta_qc$nCount_RNA, 0.01),
    quantile(meta_qc$percent.mt, 0.01)
  ),
  
  Q05 = c(
    quantile(meta_qc$nFeature_RNA, 0.05),
    quantile(meta_qc$nCount_RNA, 0.05),
    quantile(meta_qc$percent.mt, 0.05)
  ),
  
  Median = c(
    median(meta_qc$nFeature_RNA),
    median(meta_qc$nCount_RNA),
    median(meta_qc$percent.mt)
  ),
  
  Q95 = c(
    quantile(meta_qc$nFeature_RNA, 0.95),
    quantile(meta_qc$nCount_RNA, 0.95),
    quantile(meta_qc$percent.mt, 0.95)
  ),
  
  Q99 = c(
    quantile(meta_qc$nFeature_RNA, 0.99),
    quantile(meta_qc$nCount_RNA, 0.99),
    quantile(meta_qc$percent.mt, 0.99)
  ),
  
  Max = c(
    max(meta_qc$nFeature_RNA),
    max(meta_qc$nCount_RNA),
    max(meta_qc$percent.mt)
  )
)


write.csv(
  qc_after_quantiles,
  file =
    "GSE294300_QC_quantiles_after_filter.csv",
  row.names = FALSE
)


############################################################
## 12. 保存QC后、doublet前对象
############################################################

saveRDS(
  scRNA,
  file =
    "GSE294300_scRNA_after_QC_before_doublet.rds",
  compress = FALSE
)


############################################################
## 13. JoinLayers
##
## 后面转换SingleCellExperiment并运行scDblFinder
############################################################

cat("\nJoining Seurat layers...\n")

scRNA <- JoinLayers(
  scRNA,
  assay = "RNA"
)

gc()


############################################################
## 14. 转换SingleCellExperiment
############################################################

cat("\nConverting to SingleCellExperiment...\n")

sce <- as.SingleCellExperiment(
  scRNA,
  assay = "RNA"
)


############################################################
## 15. 检查sample信息
############################################################

cat("\nSamples in SCE:\n")

print(
  table(
    sce$sample
  )
)


############################################################
## 16. scDblFinder
##
## 关键：
## samples = "sample"
##
## 意味着36个10X library分别评估doublet，
## 而不是把214k细胞视为一个capture。
############################################################

cat("\n")
cat("############################################################\n")
cat("## Running scDblFinder\n")
cat("############################################################\n")


set.seed(12345)


sce <- scDblFinder::scDblFinder(
  
  sce,
  
  samples = "sample",
  
  BPPARAM =
    BiocParallel::SerialParam()
)


############################################################
## 17. 检查scDblFinder输出
############################################################

cat("\nOverall doublet classification:\n")

print(
  table(
    sce$scDblFinder.class
  )
)


############################################################
## 18. 写回Seurat
############################################################

match_index <- match(
  Cells(scRNA),
  colnames(sce)
)


if (any(is.na(match_index))) {
  stop(
    "Seurat和SingleCellExperiment细胞名无法完全匹配。"
  )
}


scRNA$scDblFinder.score <-
  sce$scDblFinder.score[
    match_index
  ]


scRNA$scDblFinder.class <-
  as.character(
    sce$scDblFinder.class[
      match_index
    ]
  )


############################################################
## 19. 每个样本doublet情况
############################################################

doublet_by_sample <- scRNA@meta.data %>%
  
  dplyr::group_by(
    sample,
    patient,
    group
  ) %>%
  
  dplyr::summarise(
    
    total_cells = dplyr::n(),
    
    singlets =
      sum(
        scDblFinder.class == "singlet"
      ),
    
    doublets =
      sum(
        scDblFinder.class == "doublet"
      ),
    
    doublet_rate =
      doublets / total_cells,
    
    .groups = "drop"
  )


write.csv(
  doublet_by_sample,
  file =
    "GSE294300_scDblFinder_by_sample.csv",
  row.names = FALSE
)


############################################################
## 20. 打印36个样本doublet情况
############################################################

cat("\n===============================\n")
cat("Doublet rate by sample\n")
cat("===============================\n")

print(
  as.data.frame(doublet_by_sample),
  row.names = FALSE
)


############################################################
## 21. doublet rate作图
############################################################

p_doublet <- ggplot(
  
  doublet_by_sample,
  
  aes(
    x = sample,
    y = doublet_rate,
    fill = group
  )
  
) +
  
  geom_col(
    width = 0.75
  ) +
  
  scale_y_continuous(
    labels = scales::percent_format(
      accuracy = 1
    )
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  labs(
    x = NULL,
    y = "Doublet rate",
    fill = NULL
  ) +
  
  theme(
    
    axis.text.x =
      element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5,
        size = 7,
        color = "black"
      ),
    
    axis.text.y =
      element_text(
        color = "black"
      ),
    
    legend.position =
      "right",
    
    panel.grid =
      element_blank()
  )


ggsave(
  filename =
    "GSE294300_scDblFinder_doublet_rate.pdf",
  plot = p_doublet,
  width = 12,
  height = 5
)


ggsave(
  filename =
    "GSE294300_scDblFinder_doublet_rate.png",
  plot = p_doublet,
  width = 12,
  height = 5,
  dpi = 300
)


############################################################
## 22. 去除doublets
############################################################

cells_before_doublet <-
  ncol(scRNA)


scRNA <- subset(
  scRNA,
  subset =
    scDblFinder.class == "singlet"
)


cells_after_doublet <-
  ncol(scRNA)


############################################################
## 23. 最终统计
############################################################

cat("\n\n")
cat("############################################################\n")
cat("## Final QC + doublet summary\n")
cat("############################################################\n")


cat(
  "\nOriginal cells:",
  nrow(meta_before)
)


cat(
  "\nAfter QC:",
  cells_before_doublet
)


cat(
  "\nAfter doublet removal:",
  cells_after_doublet
)


cat(
  "\nQC removed:",
  nrow(meta_before) -
    cells_before_doublet
)


cat(
  "\nDoublets removed:",
  cells_before_doublet -
    cells_after_doublet
)


cat(
  "\nFinal retained fraction:",
  cells_after_doublet /
    nrow(meta_before)
)


cat("\n\nFinal cells by group:\n")

print(
  table(
    scRNA$group
  )
)


cat("\nFinal cells by patient/group:\n")

print(
  table(
    scRNA$patient,
    scRNA$group
  )
)


############################################################
## 24. 检查是否36个样本仍全部存在
############################################################

cat(
  "\nSamples remaining:",
  length(
    unique(
      scRNA$sample
    )
  ),
  "\n"
)

cat(
  "Patients remaining:",
  length(
    unique(
      scRNA$patient
    )
  ),
  "\n"
)


############################################################
## 25. 三基因检查
############################################################

cat("\n===============================\n")
cat("Three MR genes\n")
cat("===============================\n")


for (g in target_genes) {
  
  hits <- grep(
    paste0(
      "^",
      g,
      "(\\.[0-9]+)?$"
    ),
    rownames(scRNA),
    value = TRUE
  )
  
  cat(
    g,
    " -> ",
    paste(
      hits,
      collapse = ", "
    ),
    "\n",
    sep = ""
  )
}


############################################################
## 26. 删除临时SCE释放内存
############################################################

rm(sce)
gc()


############################################################
## 27. 保存最终QC + singlet对象
############################################################

saveRDS(
  
  scRNA,
  
  file =
    "GSE294300_scRNA_after_QC_singlets.rds",
  
  compress = FALSE
)


############################################################
## 28. 保存最终metadata
############################################################

write.csv(
  
  scRNA@meta.data,
  
  file =
    "GSE294300_metadata_after_QC_singlets.csv",
  
  row.names = TRUE
)


############################################################
## 29. 最终输出
############################################################

cat("\n\n")
cat("############################################################\n")
cat("## GSE294300 Step 2 completed\n")
cat("############################################################\n")

cat(
  "\nFinal cells:",
  ncol(scRNA)
)

cat(
  "\nSamples:",
  length(unique(scRNA$sample))
)

cat(
  "\nPatients:",
  length(unique(scRNA$patient))
)

cat("\n\nGroup:\n")

print(
  table(scRNA$group)
)

cat("\nQC after filtering:\n")

print(
  qc_after_quantiles,
  row.names = FALSE
)

cat(
  "\nSaved: GSE294300_scRNA_after_QC_singlets.rds\n"
)
############################################################
## GSE294300 CRC scRNA-seq
## Step 3：标准化、降维、Harmony、UMAP、多分辨率聚类
##
## Input:
## GSE294300_scRNA_after_QC_singlets.rds
##
## 当前数据：
## 194604 singlets
## 36 samples
## 18 patients
## paired Normal / Tumor
##
## Candidate genes:
## AKR1C2 / MAPK3 / MAP3K11
##
## IMPORTANT:
## Harmony仅用于降维、聚类和可视化。
## 后续Tumor vs Normal统计不使用Harmony表达矩阵，
## 而使用patient-level paired / pseudobulk分析。
############################################################


############################################################
## 0. 清空环境
############################################################

rm(list = ls())

options(stringsAsFactors = FALSE)
Sys.setenv(LANGUAGE = "en")

## Seurat在大对象操作后主动gc，降低内存压力
options(Seurat.memsafe = TRUE)

set.seed(12345)


############################################################
## 1. 检查并加载R包
############################################################

required_packages <- c(
  "Seurat",
  "Matrix",
  "dplyr",
  "ggplot2",
  "patchwork",
  "harmony",
  "clustree"
)


missing_packages <- required_packages[
  !sapply(
    required_packages,
    requireNamespace,
    quietly = TRUE
  )
]


if (length(missing_packages) > 0) {
  
  stop(
    "以下R包尚未安装：",
    paste(
      missing_packages,
      collapse = ", "
    )
  )
}


library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
library(patchwork)
library(harmony)
library(clustree)


############################################################
## 2. 工作目录
############################################################

data_dir <- "C:/Users/92066/Desktop/单细胞孟德尔/Figure7/GSE294300"

setwd(data_dir)


input_file <- file.path(
  data_dir,
  "GSE294300_scRNA_after_QC_singlets.rds"
)


if (!file.exists(input_file)) {
  
  stop(
    "没有找到文件：\n",
    input_file
  )
}


############################################################
## 3. 从QC + doublet后的干净对象重新开始
##
## 不使用当前Global Environment里可能已经部分处理过的对象
############################################################

cat("\n====================================================\n")
cat("Reading input Seurat object\n")
cat("====================================================\n")


scRNA <- readRDS(
  input_file
)


DefaultAssay(scRNA) <- "RNA"


############################################################
## 4. 基本数据结构检查
############################################################

required_meta <- c(
  "sample",
  "patient",
  "group"
)


missing_meta <- setdiff(
  required_meta,
  colnames(scRNA@meta.data)
)


if (length(missing_meta) > 0) {
  
  stop(
    "meta.data缺少：",
    paste(
      missing_meta,
      collapse = ", "
    )
  )
}


if (anyNA(scRNA$sample)) {
  stop("sample中存在NA。")
}


if (anyNA(scRNA$patient)) {
  stop("patient中存在NA。")
}


if (anyNA(scRNA$group)) {
  stop("group中存在NA。")
}


## 固定Normal / Tumor顺序
scRNA$group <- factor(
  as.character(scRNA$group),
  levels = c(
    "Normal",
    "Tumor"
  )
)


cat("\n====================================================\n")
cat("Input object summary\n")
cat("====================================================\n")

cat(
  "Genes    :",
  nrow(scRNA),
  "\n"
)

cat(
  "Cells    :",
  ncol(scRNA),
  "\n"
)

cat(
  "Samples  :",
  length(
    unique(
      scRNA$sample
    )
  ),
  "\n"
)

cat(
  "Patients :",
  length(
    unique(
      scRNA$patient
    )
  ),
  "\n"
)


cat("\nCells by group:\n")

print(
  table(
    scRNA$group
  )
)


############################################################
## 5. 检查18位患者paired结构
############################################################

pair_table <- table(
  scRNA$patient,
  scRNA$group
)


cat("\n====================================================\n")
cat("Patient × group cell numbers\n")
cat("====================================================\n")

print(pair_table)


if (
  nrow(pair_table) != 18 ||
  ncol(pair_table) != 2
) {
  
  stop(
    "patient/group数量异常，",
    "不是18位患者 × 2组。"
  )
}


if (any(pair_table == 0)) {
  
  stop(
    "至少一位患者缺少Normal或Tumor细胞。"
  )
}


cat(
  "\nPaired patient structure: PASS\n"
)


############################################################
## 6. 检查RNA layers
##
## Step 2理论上已经JoinLayers。
## 如果仍有多个counts layer，则在这里统一Join。
############################################################

rna_layers <- Layers(
  scRNA[["RNA"]]
)


cat("\nRNA layers before analysis:\n")

print(rna_layers)


count_layers <- grep(
  "^counts",
  rna_layers,
  value = TRUE
)


if (length(count_layers) > 1) {
  
  cat(
    "\nMultiple count layers detected. ",
    "Joining RNA layers...\n",
    sep = ""
  )
  
  scRNA <- JoinLayers(
    scRNA,
    assay = "RNA"
  )
  
  gc()
}


cat("\nRNA layers used for analysis:\n")

print(
  Layers(
    scRNA[["RNA"]]
  )
)


############################################################
## 7. 固定三个MR候选基因
############################################################

target_genes <- c(
  "AKR1C2",
  "MAPK3",
  "MAP3K11"
)


missing_target_genes <- setdiff(
  target_genes,
  rownames(scRNA)
)


if (length(missing_target_genes) > 0) {
  
  stop(
    "以下目标基因不存在于Seurat对象中：",
    paste(
      missing_target_genes,
      collapse = ", "
    )
  )
}


cat("\nThree MR genes:\n")

for (g in target_genes) {
  
  cat(
    g,
    " -> FOUND\n",
    sep = ""
  )
}


############################################################
## 8. 分析参数统一放在这里
############################################################

n_hvg <- 4000

n_pca <- 50

dims_use <- 1:30

umap_neighbors <- 40

umap_min_dist <- 0.3


resolution_labels <- c(
  "0.1",
  "0.2",
  "0.3",
  "0.4",
  "0.5",
  "0.6",
  "0.8",
  "1.0"
)


resolution_values <- as.numeric(
  resolution_labels
)


############################################################
## 9. NormalizeData
############################################################

cat("\n====================================================\n")
cat("NormalizeData\n")
cat("====================================================\n")


scRNA <- NormalizeData(
  object = scRNA,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = TRUE
)


############################################################
## 10. FindVariableFeatures
##
## 使用4000 HVGs
############################################################

cat("\n====================================================\n")
cat("FindVariableFeatures\n")
cat("====================================================\n")


scRNA <- FindVariableFeatures(
  object = scRNA,
  assay = "RNA",
  selection.method = "vst",
  nfeatures = n_hvg,
  verbose = TRUE
)


n_variable <- length(
  VariableFeatures(scRNA)
)


cat(
  "\nNumber of variable features:",
  n_variable,
  "\n"
)


if (n_variable != n_hvg) {
  
  warning(
    "实际获得的HVG数量不是4000：",
    n_variable
  )
}


############################################################
## 11. 保存HVG列表
############################################################

hvg_df <- data.frame(
  rank = seq_along(
    VariableFeatures(scRNA)
  ),
  gene = VariableFeatures(scRNA),
  stringsAsFactors = FALSE
)


write.csv(
  hvg_df,
  file = "GSE294300_top4000_HVG.csv",
  row.names = FALSE
)


############################################################
## 12. HVG图
############################################################

top10_hvg <- head(
  VariableFeatures(scRNA),
  10
)


p_hvg <- VariableFeaturePlot(
  scRNA
)


p_hvg <- LabelPoints(
  plot = p_hvg,
  points = top10_hvg,
  repel = TRUE
) +
  theme_classic(
    base_size = 12
  )


ggsave(
  filename = "GSE294300_HVG_top4000.pdf",
  plot = p_hvg,
  width = 7,
  height = 5.5
)


ggsave(
  filename = "GSE294300_HVG_top4000.png",
  plot = p_hvg,
  width = 7,
  height = 5.5,
  dpi = 300
)


############################################################
## 13. ScaleData
##
## 只scale 4000个HVG。
##
## 不：
## features = rownames(scRNA)
##
## 也暂时不回归percent.mt，
## 因为已经进行了线粒体比例QC，
## 避免不必要地移除真实生物学信号。
############################################################

cat("\n====================================================\n")
cat("ScaleData: 4000 HVGs\n")
cat("====================================================\n")


scRNA <- ScaleData(
  object = scRNA,
  assay = "RNA",
  features = VariableFeatures(scRNA),
  verbose = TRUE
)


############################################################
## 14. PCA：计算50 PCs
############################################################

cat("\n====================================================\n")
cat("RunPCA: 50 PCs\n")
cat("====================================================\n")


scRNA <- RunPCA(
  object = scRNA,
  assay = "RNA",
  features = VariableFeatures(scRNA),
  npcs = n_pca,
  approx = TRUE,
  seed.use = 12345,
  verbose = TRUE
)


############################################################
## 15. PCA完整性检查
############################################################

if (!"pca" %in% Reductions(scRNA)) {
  
  stop(
    "RunPCA结束后没有生成pca reduction。"
  )
}


pca_embeddings <- Embeddings(
  scRNA,
  reduction = "pca"
)


pca_dim <- dim(
  pca_embeddings
)


cat("\nPCA dimensions:\n")
print(pca_dim)


if (pca_dim[1] != ncol(scRNA)) {
  
  stop(
    "PCA细胞数与Seurat对象细胞数不一致。"
  )
}


if (pca_dim[2] < max(dims_use)) {
  
  stop(
    "PCA维度不足，无法使用PC1-30。"
  )
}


if (!all(is.finite(pca_embeddings))) {
  
  stop(
    "PCA embedding存在NA/Inf。"
  )
}


rm(pca_embeddings)
gc()


############################################################
## 16. ElbowPlot
############################################################

p_elbow <- ElbowPlot(
  object = scRNA,
  ndims = 50
) +
  theme_classic(
    base_size = 13
  )


ggsave(
  filename = "GSE294300_ElbowPlot_PC1_50.pdf",
  plot = p_elbow,
  width = 6,
  height = 4.5
)


ggsave(
  filename = "GSE294300_ElbowPlot_PC1_50.png",
  plot = p_elbow,
  width = 6,
  height = 4.5,
  dpi = 300
)


############################################################
## 17. PCA按group展示
##
## 仅描述性展示，不进行cell-level显著性检验
############################################################

set.seed(12345)


p_pca_group <- DimPlot(
  object = scRNA,
  reduction = "pca",
  group.by = "group",
  pt.size = 0.08,
  shuffle = TRUE,
  raster = TRUE
) +
  theme_classic(
    base_size = 13
  ) +
  ggtitle("")


ggsave(
  filename = "GSE294300_PCA_by_group.pdf",
  plot = p_pca_group,
  width = 7,
  height = 6
)


ggsave(
  filename = "GSE294300_PCA_by_group.png",
  plot = p_pca_group,
  width = 7,
  height = 6,
  dpi = 300
)


############################################################
## 18. PCA按patient展示
############################################################

set.seed(12345)


p_pca_patient <- DimPlot(
  object = scRNA,
  reduction = "pca",
  group.by = "patient",
  pt.size = 0.08,
  shuffle = TRUE,
  raster = TRUE
) +
  theme_classic(
    base_size = 13
  ) +
  ggtitle("")


ggsave(
  filename = "GSE294300_PCA_by_patient.pdf",
  plot = p_pca_patient,
  width = 8.5,
  height = 6
)


ggsave(
  filename = "GSE294300_PCA_by_patient.png",
  plot = p_pca_patient,
  width = 8.5,
  height = 6,
  dpi = 300
)


############################################################
## 19. 保存PCA检查点
##
## 后面即使Harmony/UMAP发生问题，
## 也无需重新Normalize/Scale/PCA。
############################################################

saveRDS(
  scRNA,
  file = "GSE294300_scRNA_after_PCA.rds",
  compress = FALSE
)


############################################################
## 20. Harmony
##
## 关键：
## group.by.vars = "patient"
##
## 原因：
## 18位患者均同时具有Normal和Tumor，
## patient是donor/batch来源；
## 不把group作为Harmony校正变量，
## 避免主动消除Normal/Tumor生物学差异。
##
## Harmony仅用于：
## clustering / visualization
##
## 后续DE/pseudobulk使用raw RNA counts，
## 不使用Harmony embedding进行表达统计。
############################################################

cat("\n====================================================\n")
cat("RunHarmony by patient\n")
cat("====================================================\n")


## 如果对象中意外存在旧harmony，
## 先删除再重新计算
if ("harmony" %in% Reductions(scRNA)) {
  
  scRNA[["harmony"]] <- NULL
}


scRNA <- harmony::RunHarmony(
  object = scRNA,
  group.by.vars = "patient",
  reduction.use = "pca",
  dims.use = dims_use,
  reduction.save = "harmony",
  project.dim = FALSE,
  verbose = TRUE
)


############################################################
## 21. Harmony严格检查
############################################################

cat("\n====================================================\n")
cat("Harmony check\n")
cat("====================================================\n")


cat("\nAvailable reductions:\n")

print(
  Reductions(scRNA)
)


if (!"harmony" %in% Reductions(scRNA)) {
  
  stop(
    "Harmony运行结束，但Seurat对象中不存在'harmony' reduction。"
  )
}


harmony_embeddings <- Embeddings(
  scRNA,
  reduction = "harmony"
)


harmony_dim <- dim(
  harmony_embeddings
)


cat("\nHarmony dimensions:\n")

print(
  harmony_dim
)


if (harmony_dim[1] != ncol(scRNA)) {
  
  stop(
    "Harmony embedding的细胞数与Seurat对象不一致。"
  )
}


if (harmony_dim[2] != length(dims_use)) {
  
  stop(
    "Harmony embedding不是30维。"
  )
}


if (!all(is.finite(harmony_embeddings))) {
  
  stop(
    "Harmony embedding存在NA或Inf。"
  )
}


cat(
  "\nHarmony reduction successfully created.\n"
)


rm(harmony_embeddings)
gc()


############################################################
## 22. 保存Harmony检查点
############################################################

saveRDS(
  scRNA,
  file = "GSE294300_scRNA_after_Harmony.rds",
  compress = FALSE
)


############################################################
## 23. UMAP
##
## 输入：
## harmony dimensions 1:30
##
## Seurat/uwot默认metric本身就是cosine；
## 这里显式写出，保证参数可复现。
############################################################

cat("\n====================================================\n")
cat("RunUMAP using Harmony dimensions 1:30\n")
cat("====================================================\n")


## 如果意外已有旧UMAP，则删除
if ("umap" %in% Reductions(scRNA)) {
  
  scRNA[["umap"]] <- NULL
}


scRNA <- RunUMAP(
  object = scRNA,
  reduction = "harmony",
  dims = dims_use,
  n.neighbors = umap_neighbors,
  min.dist = umap_min_dist,
  metric = "cosine",
  seed.use = 12345,
  reduction.name = "umap",
  verbose = TRUE
)


############################################################
## 24. UMAP检查
############################################################

if (!"umap" %in% Reductions(scRNA)) {
  
  stop(
    "RunUMAP完成后没有生成umap reduction。"
  )
}


umap_dim <- dim(
  Embeddings(
    scRNA,
    reduction = "umap"
  )
)


cat("\nUMAP dimensions:\n")

print(
  umap_dim
)


if (
  umap_dim[1] != ncol(scRNA) ||
  umap_dim[2] != 2
) {
  
  stop(
    "UMAP维度异常。"
  )
}


############################################################
## 25. UMAP by group
############################################################

set.seed(12345)


p_umap_group <- DimPlot(
  object = scRNA,
  reduction = "umap",
  group.by = "group",
  pt.size = 0.08,
  shuffle = TRUE,
  raster = TRUE
) +
  theme_classic(
    base_size = 13
  ) +
  ggtitle("")


ggsave(
  filename = "GSE294300_UMAP_by_group.pdf",
  plot = p_umap_group,
  width = 7,
  height = 6
)


ggsave(
  filename = "GSE294300_UMAP_by_group.png",
  plot = p_umap_group,
  width = 7,
  height = 6,
  dpi = 300
)


############################################################
## 26. UMAP by patient
############################################################

set.seed(12345)


p_umap_patient <- DimPlot(
  object = scRNA,
  reduction = "umap",
  group.by = "patient",
  pt.size = 0.08,
  shuffle = TRUE,
  raster = TRUE
) +
  theme_classic(
    base_size = 13
  ) +
  ggtitle("")


ggsave(
  filename = "GSE294300_UMAP_by_patient.pdf",
  plot = p_umap_patient,
  width = 8.5,
  height = 6
)


ggsave(
  filename = "GSE294300_UMAP_by_patient.png",
  plot = p_umap_patient,
  width = 8.5,
  height = 6,
  dpi = 300
)


############################################################
## 27. Normal / Tumor split UMAP
############################################################

p_umap_split <- DimPlot(
  object = scRNA,
  reduction = "umap",
  group.by = "group",
  split.by = "group",
  pt.size = 0.08,
  raster = TRUE
) +
  theme_classic(
    base_size = 13
  ) +
  ggtitle("")


ggsave(
  filename = "GSE294300_UMAP_Normal_Tumor_split.pdf",
  plot = p_umap_split,
  width = 12,
  height = 5.5
)


ggsave(
  filename = "GSE294300_UMAP_Normal_Tumor_split.png",
  plot = p_umap_split,
  width = 12,
  height = 5.5,
  dpi = 300
)


############################################################
## 28. FindNeighbors
##
## 明确指定graph名称：
##
## harmony_nn
## harmony_snn
##
## 后面FindClusters明确使用harmony_snn，
## 避免依赖Seurat自动命名。
############################################################

cat("\n====================================================\n")
cat("FindNeighbors using Harmony\n")
cat("====================================================\n")


## 删除可能存在的同名旧graph
if ("harmony_nn" %in% Graphs(scRNA)) {
  
  scRNA[["harmony_nn"]] <- NULL
}


if ("harmony_snn" %in% Graphs(scRNA)) {
  
  scRNA[["harmony_snn"]] <- NULL
}


scRNA <- FindNeighbors(
  object = scRNA,
  reduction = "harmony",
  dims = dims_use,
  k.param = 20,
  compute.SNN = TRUE,
  graph.name = c(
    "harmony_nn",
    "harmony_snn"
  ),
  verbose = TRUE
)


############################################################
## 29. Graph检查
############################################################

cat("\nGraphs in Seurat object:\n")

print(
  Graphs(scRNA)
)


if (!all(
  c(
    "harmony_nn",
    "harmony_snn"
  ) %in% Graphs(scRNA)
)) {
  
  stop(
    "FindNeighbors未正确生成harmony_nn和harmony_snn。"
  )
}


############################################################
## 30. 多分辨率聚类
##
## 为提高不同Seurat版本兼容性，
## 不一次向resolution传入向量。
##
## 逐个resolution运行，并明确给每一列命名：
##
## harmony_snn_res.0.1
## harmony_snn_res.0.2
## ...
## harmony_snn_res.1.0
############################################################

cat("\n====================================================\n")
cat("Multi-resolution clustering\n")
cat("====================================================\n")


for (
  i in seq_along(
    resolution_values
  )
) {
  
  current_res <- resolution_values[i]
  
  current_label <- resolution_labels[i]
  
  current_cluster_name <- paste0(
    "harmony_snn_res.",
    current_label
  )
  
  
  cat(
    "\nClustering resolution = ",
    current_label,
    "\n",
    sep = ""
  )
  
  
  scRNA <- FindClusters(
    object = scRNA,
    graph.name = "harmony_snn",
    cluster.name = current_cluster_name,
    resolution = current_res,
    algorithm = 1,
    n.start = 10,
    n.iter = 10,
    random.seed = 12345,
    verbose = TRUE
  )
}


############################################################
## 31. 明确构造resolution列名
############################################################

res_cols <- paste0(
  "harmony_snn_res.",
  resolution_labels
)


missing_res_cols <- setdiff(
  res_cols,
  colnames(scRNA@meta.data)
)


if (length(missing_res_cols) > 0) {
  
  stop(
    "以下resolution列未成功生成：",
    paste(
      missing_res_cols,
      collapse = ", "
    )
  )
}


############################################################
## 32. 统计不同resolution的cluster数量
############################################################

cluster_number <- data.frame(
  
  resolution = resolution_labels,
  
  resolution_column = res_cols,
  
  n_clusters = sapply(
    res_cols,
    function(x) {
      
      length(
        unique(
          scRNA@meta.data[[x]]
        )
      )
    }
  ),
  
  stringsAsFactors = FALSE
)


cat("\n====================================================\n")
cat("Cluster number by resolution\n")
cat("====================================================\n")


print(
  cluster_number,
  row.names = FALSE
)


write.csv(
  cluster_number,
  file = "GSE294300_cluster_number_by_resolution.csv",
  row.names = FALSE
)

############################################################
## 33. 每个resolution的cluster细胞数
############################################################

cluster_size_list <- lapply(
  seq_along(res_cols),
  function(i) {
    
    cluster_vector <- scRNA@meta.data[[res_cols[i]]]
    
    temp <- as.data.frame(
      table(cluster_vector),
      stringsAsFactors = FALSE
    )
    
    colnames(temp) <- c(
      "cluster",
      "n_cells"
    )
    
    temp$resolution <- resolution_labels[i]
    
    temp$resolution_column <- res_cols[i]
    
    temp
  }
)


cluster_size_df <- dplyr::bind_rows(
  cluster_size_list
)


## 调整列顺序
cluster_size_df <- cluster_size_df %>%
  dplyr::select(
    resolution,
    resolution_column,
    cluster,
    n_cells
  )


cat("\n====================================================\n")
cat("Cluster size by resolution\n")
cat("====================================================\n")

print(
  as.data.frame(cluster_size_df),
  row.names = FALSE
)


write.csv(
  cluster_size_df,
  file = "GSE294300_cluster_size_by_resolution.csv",
  row.names = FALSE
)


############################################################
## 34. 检查是否存在极小cluster
##
## 这里只用于辅助判断resolution
## 暂时不删除任何cluster
############################################################

small_cluster_summary <- cluster_size_df %>%
  dplyr::group_by(
    resolution
  ) %>%
  dplyr::summarise(
    n_clusters = dplyr::n(),
    smallest_cluster = min(n_cells),
    median_cluster_size = median(n_cells),
    n_cluster_lt_100 = sum(n_cells < 100),
    n_cluster_lt_500 = sum(n_cells < 500),
    .groups = "drop"
  )


cat("\n====================================================\n")
cat("Small-cluster summary\n")
cat("====================================================\n")

print(
  as.data.frame(small_cluster_summary),
  row.names = FALSE
)


write.csv(
  small_cluster_summary,
  file = "GSE294300_small_cluster_summary.csv",
  row.names = FALSE
)


############################################################
## 35. Clustree
##
## 前面已经明确生成：
## harmony_snn_res.0.1
## harmony_snn_res.0.2
## ...
##
## 所以prefix固定，不再自动猜
############################################################

cat("\n====================================================\n")
cat("Clustree\n")
cat("====================================================\n")


## 再次确认所有resolution列存在
missing_res_cols <- setdiff(
  res_cols,
  colnames(scRNA@meta.data)
)


if (length(missing_res_cols) > 0) {
  
  stop(
    "以下聚类列不存在：",
    paste(
      missing_res_cols,
      collapse = ", "
    )
  )
}


p_clustree <- clustree::clustree(
  scRNA,
  prefix = "harmony_snn_res."
) +
  theme(
    plot.margin = margin(
      20,
      40,
      20,
      40
    )
  )


ggsave(
  filename = "GSE294300_clustree_resolution.pdf",
  plot = p_clustree,
  width = 18,
  height = 10,
  limitsize = FALSE
)


ggsave(
  filename = "GSE294300_clustree_resolution.png",
  plot = p_clustree,
  width = 18,
  height = 10,
  dpi = 300,
  limitsize = FALSE
)


############################################################
## 36. 各resolution输出UMAP预览
############################################################

for (i in seq_along(res_cols)) {
  
  current_col <- res_cols[i]
  
  current_label <- resolution_labels[i]
  
  
  cat(
    "Plotting resolution ",
    current_label,
    "\n",
    sep = ""
  )
  
  
  set.seed(12345)
  
  
  p_temp <- DimPlot(
    object = scRNA,
    reduction = "umap",
    group.by = current_col,
    label = TRUE,
    repel = TRUE,
    pt.size = 0.06,
    raster = TRUE
  ) +
    theme_classic(
      base_size = 12
    ) +
    ggtitle(
      paste0(
        "Resolution ",
        current_label
      )
    ) +
    theme(
      legend.title = element_blank(),
      panel.grid = element_blank()
    )
  
  
  ggsave(
    filename = paste0(
      "GSE294300_UMAP_resolution_",
      current_label,
      ".pdf"
    ),
    plot = p_temp,
    width = 7.5,
    height = 6.5,
    limitsize = FALSE
  )
  
  
  ggsave(
    filename = paste0(
      "GSE294300_UMAP_resolution_",
      current_label,
      ".png"
    ),
    plot = p_temp,
    width = 7.5,
    height = 6.5,
    dpi = 300,
    limitsize = FALSE
  )
}


############################################################
## 37. 再次确认目前不选择final resolution
############################################################

cat("\n====================================================\n")
cat("IMPORTANT\n")
cat("====================================================\n")

cat(
  "当前仅完成multi-resolution clustering。\n",
  "暂时不要把seurat_clusters作为最终聚类结果。\n",
  "最终resolution将在查看cluster number、clustree、UMAP和markers后决定。\n",
  sep = ""
)


############################################################
## 38. 保存完整multi-resolution对象
############################################################

output_rds <- file.path(
  data_dir,
  "GSE294300_scRNA_Harmony_multi_resolution.rds"
)


saveRDS(
  scRNA,
  file = output_rds,
  compress = FALSE
)


############################################################
## 39. 保存分析参数
############################################################

parameter_df <- data.frame(
  
  parameter = c(
    "genes",
    "cells",
    "patients",
    "samples",
    "HVG_number",
    "PCA_number",
    "Harmony_variable",
    "Harmony_dims",
    "UMAP_n_neighbors",
    "UMAP_min_dist",
    "UMAP_metric",
    "FindNeighbors_k",
    "Clustering_algorithm",
    "Resolutions"
  ),
  
  value = c(
    nrow(scRNA),
    ncol(scRNA),
    length(unique(scRNA$patient)),
    length(unique(scRNA$sample)),
    n_hvg,
    n_pca,
    "patient",
    paste0(
      min(dims_use),
      ":",
      max(dims_use)
    ),
    umap_neighbors,
    umap_min_dist,
    "cosine",
    20,
    "Louvain algorithm 1",
    paste(
      resolution_labels,
      collapse = ", "
    )
  ),
  
  stringsAsFactors = FALSE
)


write.csv(
  parameter_df,
  file = "GSE294300_Step3_analysis_parameters.csv",
  row.names = FALSE
)


############################################################
## 40. 保存sessionInfo
## 用于方法学可重复性和审稿回复
############################################################

session_info_text <- capture.output(
  sessionInfo()
)


writeLines(
  session_info_text,
  con = "GSE294300_Step3_sessionInfo.txt"
)


############################################################
## 41. 最终完整性检查
############################################################

harmony_ok <- "harmony" %in% Reductions(scRNA)

umap_ok <- "umap" %in% Reductions(scRNA)

graph_ok <- all(
  c(
    "harmony_nn",
    "harmony_snn"
  ) %in% Graphs(scRNA)
)

resolution_ok <- all(
  res_cols %in%
    colnames(scRNA@meta.data)
)


if (!harmony_ok) {
  stop("最终检查失败：harmony reduction不存在。")
}


if (!umap_ok) {
  stop("最终检查失败：UMAP不存在。")
}


if (!graph_ok) {
  stop("最终检查失败：Harmony graph不存在。")
}


if (!resolution_ok) {
  stop("最终检查失败：部分resolution结果不存在。")
}


############################################################
## 42. 最终一次性输出
############################################################

cat("\n\n")
cat("############################################################\n")
cat("## GSE294300 Step 3 completed successfully\n")
cat("############################################################\n")


cat(
  "\nGenes:",
  nrow(scRNA)
)


cat(
  "\nCells:",
  ncol(scRNA)
)


cat(
  "\nSamples:",
  length(
    unique(
      scRNA$sample
    )
  )
)


cat(
  "\nPatients:",
  length(
    unique(
      scRNA$patient
    )
  )
)


cat(
  "\nHarmony reduction:",
  harmony_ok
)


cat(
  "\nUMAP reduction:",
  umap_ok
)


cat(
  "\nHarmony NN/SNN graphs:",
  graph_ok
)


cat(
  "\nAll resolution columns:",
  resolution_ok
)


cat("\n\nAvailable reductions:\n")

print(
  Reductions(scRNA)
)


cat("\nAvailable graphs:\n")

print(
  Graphs(scRNA)
)


cat("\n\nCluster numbers:\n")

print(
  cluster_number,
  row.names = FALSE
)


cat("\n\nSmall-cluster summary:\n")

print(
  as.data.frame(small_cluster_summary),
  row.names = FALSE
)


cat(
  "\n\nSaved object:\n",
  output_rds,
  "\n",
  sep = ""
)


cat(
  "\nStep 3完成。暂不选择final resolution。\n"
)
library(Seurat)

data_dir <- "C:/Users/92066/Desktop/单细胞孟德尔/Figure7/GSE294300"
setwd(data_dir)

scRNA <- readRDS(
  "GSE294300_scRNA_Harmony_multi_resolution.rds"
)

############################################################
## 正确提取 resolution = 0.6
############################################################

raw_cluster <- as.character(
  scRNA@meta.data[["harmony_snn_res.0.6"]]
)

## 关键：先character，再integer
cluster_ids <- sort(
  unique(
    as.integer(raw_cluster)
  )
)

print(cluster_ids)

## 应该是0:21
stopifnot(
  identical(
    cluster_ids,
    0:21
  )
)

############################################################
## 正确建立final_cluster
############################################################

scRNA$final_cluster <- factor(
  raw_cluster,
  levels = as.character(cluster_ids)
)

Idents(scRNA) <- "final_cluster"

############################################################
## 检查
############################################################

table(
  scRNA$final_cluster,
  useNA = "ifany"
)

cat(
  "Cluster number:",
  nlevels(scRNA$final_cluster),
  "\n"
)

cat(
  "NA cells:",
  sum(is.na(scRNA$final_cluster)),
  "\n"
)

cat(
  "Total assigned cells:",
  sum(table(scRNA$final_cluster)),
  "\n"
)
saveRDS(
  scRNA,
  file = "GSE294300_scRNA_final_res0.6.rds",
  compress = FALSE
)
############################################################
## GSE294300 CRC scRNA-seq
## Marker identification and major cell-type annotation
## Resolution = 0.6
## Clusters = 0-21
## Cells = 194604
############################################################

rm(list = ls())
options(stringsAsFactors = FALSE)
Sys.setenv(LANGUAGE = "en")
set.seed(12345)

library(Seurat)
library(dplyr)
library(ggplot2)

############################################################
## 1. Load data
############################################################

data_dir <- "C:/Users/92066/Desktop/单细胞孟德尔/Figure7/GSE294300"
setwd(data_dir)

input_file <- file.path(data_dir, "GSE294300_scRNA_final_res0.6.rds")
if (!file.exists(input_file)) stop("没有找到：", input_file)

scRNA <- readRDS(input_file)
DefaultAssay(scRNA) <- "RNA"

required_meta <- c("sample", "patient", "group", "final_cluster")
missing_meta <- setdiff(required_meta, colnames(scRNA@meta.data))
if (length(missing_meta) > 0) stop("缺少metadata：", paste(missing_meta, collapse = ", "))

scRNA$final_cluster <- factor(as.character(scRNA$final_cluster), levels = as.character(0:21))
Idents(scRNA) <- "final_cluster"

############################################################
## 2. Basic validation
############################################################

cat("\n====================================================\n")
cat("Basic validation\n")
cat("====================================================\n")
print(table(scRNA$final_cluster))
cat("Cells:", ncol(scRNA), "\n")
cat("Clusters:", nlevels(scRNA$final_cluster), "\n")
cat("Patients:", length(unique(scRNA$patient)), "\n")
cat("Samples:", length(unique(scRNA$sample)), "\n")
cat("NA clusters:", sum(is.na(scRNA$final_cluster)), "\n")

stopifnot(ncol(scRNA) == 194604)
stopifnot(nlevels(scRNA$final_cluster) == 22)
stopifnot(length(unique(scRNA$patient)) == 18)
stopifnot(length(unique(scRNA$sample)) == 36)
stopifnot(sum(is.na(scRNA$final_cluster)) == 0)

if (!"umap" %in% Reductions(scRNA)) stop("Seurat对象中没有UMAP reduction。")
if (!"harmony" %in% Reductions(scRNA)) stop("Seurat对象中没有Harmony reduction。")

rna_layers <- Layers(scRNA[["RNA"]])
if (!any(grepl("^data", rna_layers))) stop("RNA assay中没有normalized data layer。")

############################################################
## 3. Cluster size and donor representation
############################################################

cluster_size_df <- data.frame(
  cluster = names(table(scRNA$final_cluster)),
  n_cells = as.integer(table(scRNA$final_cluster)),
  stringsAsFactors = FALSE
)

write.csv(cluster_size_df, "GSE294300_cluster_size.csv", row.names = FALSE)

cluster_patient_counts <- scRNA@meta.data %>%
  dplyr::count(final_cluster, patient, name = "n_cells_patient")

cluster_donor_summary <- cluster_patient_counts %>%
  dplyr::group_by(final_cluster) %>%
  dplyr::summarise(
    n_cells = sum(n_cells_patient),
    n_patients = dplyr::n_distinct(patient),
    max_patient_fraction = max(n_cells_patient) / sum(n_cells_patient),
    dominant_patient = as.character(patient[which.max(n_cells_patient)]),
    .groups = "drop"
  )

cluster_sample_summary <- scRNA@meta.data %>%
  dplyr::group_by(final_cluster) %>%
  dplyr::summarise(n_samples = dplyr::n_distinct(sample), .groups = "drop")

cluster_donor_summary <- dplyr::left_join(
  cluster_donor_summary,
  cluster_sample_summary,
  by = "final_cluster"
)

cat("\n====================================================\n")
cat("Donor representation by cluster\n")
cat("====================================================\n")
print(as.data.frame(cluster_donor_summary), row.names = FALSE)

write.csv(
  cluster_donor_summary,
  "GSE294300_cluster_donor_summary.csv",
  row.names = FALSE
)

cluster20_source <- scRNA@meta.data %>%
  dplyr::filter(as.character(final_cluster) == "20") %>%
  dplyr::count(patient, sample, group, name = "n_cells") %>%
  dplyr::arrange(desc(n_cells))

write.csv(
  cluster20_source,
  "GSE294300_cluster20_source.csv",
  row.names = FALSE
)

############################################################
## 4. Find cluster markers
############################################################

marker_file <- "GSE294300_FindAllMarkers_res0.6.csv"
rerun_markers <- FALSE

if (file.exists(marker_file) && !rerun_markers) {
  cat("\nReading existing marker file...\n")
  markers_all <- read.csv(marker_file, stringsAsFactors = FALSE, check.names = FALSE)
} else {
  Idents(scRNA) <- "final_cluster"
  markers_all <- FindAllMarkers(
    object = scRNA,
    assay = "RNA",
    slot = "data",
    test.use = "wilcox",
    only.pos = TRUE,
    min.pct = 0.10,
    logfc.threshold = 0.25,
    max.cells.per.ident = 5000,
    random.seed = 12345,
    return.thresh = 0.05,
    densify = FALSE,
    verbose = TRUE
  )
  if (nrow(markers_all) == 0) stop("FindAllMarkers没有找到marker。")
  write.csv(markers_all, marker_file, row.names = FALSE)
}

required_marker_cols <- c("gene", "cluster", "avg_log2FC", "p_val_adj")
missing_marker_cols <- setdiff(required_marker_cols, colnames(markers_all))
if (length(missing_marker_cols) > 0) {
  stop("Marker文件缺少：", paste(missing_marker_cols, collapse = ", "))
}

markers_all <- markers_all %>%
  dplyr::arrange(cluster, p_val_adj, desc(avg_log2FC))

get_top_markers <- function(n) {
  markers_all %>%
    dplyr::filter(p_val_adj < 0.05) %>%
    dplyr::group_by(cluster) %>%
    dplyr::slice_max(order_by = avg_log2FC, n = n, with_ties = FALSE) %>%
    dplyr::ungroup()
}

top10_markers <- get_top_markers(10)
top20_markers <- get_top_markers(20)
top50_markers <- get_top_markers(50)

write.csv(top10_markers, "GSE294300_Top10_markers_res0.6.csv", row.names = FALSE)
write.csv(top20_markers, "GSE294300_Top20_markers_res0.6.csv", row.names = FALSE)
write.csv(top50_markers, "GSE294300_Top50_markers_res0.6.csv", row.names = FALSE)

############################################################
## 5. Canonical markers for 10 major cell types
############################################################

marker_list <- list(
  "Epithelial" = c("EPCAM", "CDH1", "KRT8", "KRT18", "KRT19", "KRT20"),
  "T/NK/ILC" = c("CD3D", "CD3E", "TRAC", "NKG7", "GNLY", "KLRD1", "IL7R", "KLRB1"),
  "B" = c("MS4A1", "CD79A", "CD79B", "CD19", "CD37"),
  "Plasma" = c("MZB1", "JCHAIN", "SDC1", "XBP1", "DERL3"),
  "Mast" = c("TPSAB1", "TPSB2", "KIT", "MS4A2", "CPA3"),
  "Myeloid" = c("LYZ", "LST1", "TYROBP", "FCER1G", "AIF1", "CTSS"),
  "Endothelial" = c("PECAM1", "CDH5", "VWF", "KDR", "EMCN", "RAMP2"),
  "Glial" = c("S100B", "SOX10", "PLP1", "SLC1A3", "PMP2", "CDH19"),
  "Mural cells" = c("RGS5", "PDGFRB", "CSPG4", "KCNJ8", "ABCC9", "COX4I2",
                    "MYH11", "CNN1", "ACTA2", "TAGLN", "LMOD1", "MYL9"),
  "Fibroblast" = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRA")
)

marker_check <- data.frame(
  celltype = rep(names(marker_list), lengths(marker_list)),
  gene = unlist(marker_list, use.names = FALSE),
  stringsAsFactors = FALSE
)

marker_check$present <- marker_check$gene %in% rownames(scRNA)
write.csv(
  marker_check,
  "GSE294300_canonical_marker_availability.csv",
  row.names = FALSE
)

marker_list_present <- lapply(
  marker_list,
  function(x) x[x %in% rownames(scRNA)]
)

if (any(lengths(marker_list_present) == 0)) {
  stop(
    "以下细胞类型没有可用canonical marker：",
    paste(names(marker_list_present)[lengths(marker_list_present) == 0], collapse = ", ")
  )
}

############################################################
## 6. Canonical marker DotPlot by cluster
############################################################

p_marker_cluster <- DotPlot(
  object = scRNA,
  features = marker_list_present,
  group.by = "final_cluster",
  assay = "RNA",
  dot.scale = 6
) +
  RotatedAxis() +
  theme_classic(base_size = 11) +
  labs(
    x = NULL,
    y = "Cluster",
    color = "Average\nexpression",
    size = "Percent\nexpressed"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black", size = 8),
    axis.text.y = element_text(color = "black", size = 10),
    panel.grid = element_blank()
  )

ggsave(
  "GSE294300_canonical_marker_DotPlot_clusters.pdf",
  p_marker_cluster,
  width = 22,
  height = 8.5,
  limitsize = FALSE
)

ggsave(
  "GSE294300_canonical_marker_DotPlot_clusters.png",
  p_marker_cluster,
  width = 22,
  height = 8.5,
  dpi = 300,
  limitsize = FALSE
)

############################################################
## 7. Manual major cell-type annotation
## ONLY edit this table when changing annotation
############################################################

annotation_table <- data.frame(
  cluster = as.character(0:21),
  celltype = c(
    "T/NK/ILC",     # 0
    "T/NK/ILC",     # 1
    "B",            # 2
    "Plasma",       # 3
    "T/NK/ILC",     # 4
    "Plasma",       # 5
    "T/NK/ILC",     # 6
    "Epithelial",   # 7
    "Myeloid",      # 8
    "Plasma",       # 9
    "Fibroblast",   # 10
    "Mural cells",  # 11
    "T/NK/ILC",     # 12
    "Mast",         # 13
    "Fibroblast",   # 14
    "B",            # 15
    "Endothelial",  # 16
    "Epithelial",   # 17
    "Epithelial",   # 18
    "T/NK/ILC",     # 19
    "Epithelial",   # 20
    "Glial"         # 21
  ),
  stringsAsFactors = FALSE
)

celltype_order <- c(
  "Epithelial", "T/NK/ILC", "B", "Plasma", "Mast",
  "Myeloid", "Endothelial", "Glial", "Mural cells", "Fibroblast"
)

if (nrow(annotation_table) != 22) stop("annotation_table必须包含22个cluster。")
if (!identical(annotation_table$cluster, as.character(0:21))) stop("cluster必须严格按照0-21排列。")
if (anyDuplicated(annotation_table$cluster) > 0) stop("annotation_table存在重复cluster。")

unexpected_celltypes <- setdiff(unique(annotation_table$celltype), celltype_order)
if (length(unexpected_celltypes) > 0) {
  stop("出现未定义cell type：", paste(unexpected_celltypes, collapse = ", "))
}

missing_celltypes <- setdiff(celltype_order, unique(annotation_table$celltype))
if (length(missing_celltypes) > 0) {
  stop("以下cell type没有对应cluster：", paste(missing_celltypes, collapse = ", "))
}

write.csv(
  annotation_table,
  "GSE294300_cluster_annotation_table_10types.csv",
  row.names = FALSE
)

############################################################
## 8. Write annotation into Seurat
############################################################

celltype_map <- setNames(annotation_table$celltype, annotation_table$cluster)
scRNA$celltype <- unname(celltype_map[as.character(scRNA$final_cluster)])

if (anyNA(scRNA$celltype)) stop("存在cluster未成功注释。")

scRNA$celltype <- factor(
  scRNA$celltype,
  levels = celltype_order
)

if (anyNA(scRNA$celltype)) stop("设置celltype factor后出现NA。")

cluster_celltype_table <- scRNA@meta.data %>%
  dplyr::count(final_cluster, celltype, name = "n_cells") %>%
  dplyr::arrange(as.integer(as.character(final_cluster)))

celltype_summary <- scRNA@meta.data %>%
  dplyr::count(celltype, name = "n_cells", .drop = FALSE) %>%
  dplyr::mutate(proportion = n_cells / sum(n_cells))

write.csv(
  cluster_celltype_table,
  "GSE294300_cluster_to_celltype_10types.csv",
  row.names = FALSE
)

write.csv(
  celltype_summary,
  "GSE294300_celltype_summary_10types.csv",
  row.names = FALSE
)

stopifnot(nrow(celltype_summary) == 10)
stopifnot(sum(celltype_summary$n_cells) == ncol(scRNA))

############################################################
## 9. Generate tSNE once
## Same tSNE coordinates are used before and after annotation
############################################################

if (!"tsne" %in% Reductions(scRNA)) {
  set.seed(12345)
  scRNA <- RunTSNE(
    object = scRNA,
    reduction = "harmony",
    dims = 1:30,
    reduction.name = "tsne",
    reduction.key = "tSNE_",
    seed.use = 12345,
    check_duplicates = FALSE,
    verbose = TRUE
  )
}

############################################################
## 10. Plot colors
############################################################

cluster_colors <- setNames(
  scales::hue_pal(h = c(15, 375), c = 100, l = 55)(22),
  as.character(0:21)
)

celltype_colors <- c(
  "Epithelial" = "#E64B35",
  "T/NK/ILC" = "#D89000",
  "B" = "#A6A800",
  "Plasma" = "#35A600",
  "Mast" = "#00A67D",
  "Myeloid" = "#00A6B8",
  "Endothelial" = "#0072B2",
  "Glial" = "#6655CC",
  "Mural cells" = "#C43ACF",
  "Fibroblast" = "#E83E8C"
)

############################################################
## 11. UMAP before annotation: 22 clusters
############################################################

p_umap_cluster <- DimPlot(
  object = scRNA,
  reduction = "umap",
  group.by = "final_cluster",
  cols = cluster_colors,
  label = TRUE,
  label.size = 4,
  repel = TRUE,
  pt.size = 1.5,
  raster = TRUE,
  shuffle = TRUE,
  seed = 12345
) +
  theme_classic(base_size = 13) +
  labs(x = "UMAP 1", y = "UMAP 2", title = NULL) +
  theme(
    legend.title = element_blank(),
    axis.text = element_text(color = "black"),
    panel.grid = element_blank()
  )

ggsave(
  "GSE294300_cluster_UMAP.pdf",
  p_umap_cluster,
  width = 7.5,
  height = 6.5
)

ggsave(
  "GSE294300_cluster_UMAP.png",
  p_umap_cluster,
  width = 7.5,
  height = 6.5,
  dpi = 300
)

############################################################
## 12. UMAP after annotation: 10 cell types
############################################################

p_umap_celltype <- DimPlot(
  object = scRNA,
  reduction = "umap",
  group.by = "celltype",
  cols = celltype_colors,
  label = TRUE,
  label.size = 4,
  repel = TRUE,
  pt.size = 1.5,
  raster = TRUE,
  raster.dpi = c(1200, 1200),
  shuffle = TRUE,
  seed = 12345
) +
  theme_classic(base_size = 13) +
  labs(x = "UMAP 1", y = "UMAP 2", title = NULL) +
  theme(
    legend.title = element_blank(),
    axis.text = element_text(color = "black"),
    panel.grid = element_blank()
  )

ggsave(
  "GSE294300_celltype_UMAP.pdf",
  p_umap_celltype,
  width = 8.5,
  height = 6.5
)

ggsave(
  "GSE294300_celltype_UMAP.png",
  p_umap_celltype,
  width = 8.5,
  height = 6.5,
  dpi = 300
)

############################################################
## 13. tSNE before annotation: 22 clusters
############################################################

p_tsne_cluster <- DimPlot(
  object = scRNA,
  reduction = "tsne",
  group.by = "final_cluster",
  cols = cluster_colors,
  label = TRUE,
  label.size = 4,
  repel = TRUE,
  pt.size = 1.5,
  raster = TRUE,
  shuffle = TRUE,
  seed = 12345
) +
  theme_classic(base_size = 13) +
  labs(x = "tSNE 1", y = "tSNE 2", title = NULL) +
  theme(
    legend.title = element_blank(),
    axis.text = element_text(color = "black"),
    panel.grid = element_blank()
  )

ggsave(
  "GSE294300_cluster_tSNE.pdf",
  p_tsne_cluster,
  width = 7.5,
  height = 6.5
)

ggsave(
  "GSE294300_cluster_tSNE.png",
  p_tsne_cluster,
  width = 7.5,
  height = 6.5,
  dpi = 300
)

############################################################
## 14. tSNE after annotation: 10 cell types
############################################################

p_tsne_celltype <- DimPlot(
  object = scRNA,
  reduction = "tsne",
  group.by = "celltype",
  cols = celltype_colors,
  label = TRUE,
  label.size = 4,
  repel = TRUE,
  pt.size = 1.5,
  raster = TRUE,
  raster.dpi = c(1200, 1200),
  shuffle = TRUE,
  seed = 12345
) +
  theme_classic(base_size = 13) +
  labs(x = "tSNE 1", y = "tSNE 2", title = NULL) +
  theme(
    legend.title = element_blank(),
    axis.text = element_text(color = "black"),
    panel.grid = element_blank()
  )

ggsave(
  "GSE294300_celltype_tSNE.pdf",
  p_tsne_celltype,
  width = 8.5,
  height = 6.5
)
ggsave(
  "GSE294300_celltype_tSNE.png",
  p_tsne_celltype,
  width = 8.5,
  height = 6.5,
  dpi = 600
)

############################################################
## 15. Canonical marker DotPlot after annotation
############################################################

p_marker_celltype <- DotPlot(
  object = scRNA,
  features = marker_list_present,
  group.by = "celltype",
  assay = "RNA",
  dot.scale = 6
) +
  RotatedAxis() +
  theme_classic(base_size = 11) +
  labs(
    x = NULL,
    y = "Cell type",
    color = "Average\nexpression",
    size = "Percent\nexpressed"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black", size = 8),
    axis.text.y = element_text(color = "black", size = 10),
    panel.grid = element_blank()
  )

ggsave(
  "GSE294300_canonical_marker_DotPlot_celltypes.pdf",
  p_marker_celltype,
  width = 22,
  height = 6.5,
  limitsize = FALSE
)

ggsave(
  "GSE294300_canonical_marker_DotPlot_celltypes.png",
  p_marker_celltype,
  width = 22,
  height = 6.5,
  dpi = 300,
  limitsize = FALSE
)

############################################################
## 16. Check MR genes
############################################################

target_genes <- c("AKR1C2", "MAPK3", "MAP3K11")
missing_target_genes <- setdiff(target_genes, rownames(scRNA))

cat("\n====================================================\n")
cat("MR genes\n")
cat("====================================================\n")

for (g in target_genes) {
  cat(g, ":", g %in% rownames(scRNA), "\n")
}

if (length(missing_target_genes) > 0) {
  stop("以下MR基因不存在：", paste(missing_target_genes, collapse = ", "))
}

############################################################
## 17. Final validation
############################################################

cat("\n====================================================\n")
cat("Final validation\n")
cat("====================================================\n")

cat("Cells:", ncol(scRNA), "\n")
cat("Clusters:", nlevels(scRNA$final_cluster), "\n")
cat("Cell types:", nlevels(scRNA$celltype), "\n")
cat("Patients:", length(unique(scRNA$patient)), "\n")
cat("Samples:", length(unique(scRNA$sample)), "\n")
cat("NA celltype:", sum(is.na(scRNA$celltype)), "\n")
cat("tSNE cells:", nrow(Embeddings(scRNA, "tsne")), "\n")

stopifnot(ncol(scRNA) == 194604)
stopifnot(nlevels(scRNA$final_cluster) == 22)
stopifnot(nlevels(scRNA$celltype) == 10)
stopifnot(length(unique(scRNA$patient)) == 18)
stopifnot(length(unique(scRNA$sample)) == 36)
stopifnot(sum(is.na(scRNA$celltype)) == 0)
stopifnot(nrow(Embeddings(scRNA, "tsne")) == ncol(scRNA))

############################################################
## 18. Save final object
############################################################

saveRDS(
  scRNA,
  "GSE294300_scRNA_final_annotated_with_tSNE.rds",
  compress = FALSE
)

writeLines(
  capture.output(sessionInfo()),
  "GSE294300_annotation_sessionInfo.txt"
)

cat("\n############################################################\n")
cat("## ANNOTATION COMPLETED\n")
cat("############################################################\n")
cat("\nCluster-to-cell-type mapping:\n")
print(as.data.frame(cluster_celltype_table), row.names = FALSE)
cat("\nCell-type summary:\n")
print(as.data.frame(celltype_summary), row.names = FALSE)
############################################################
## GSE294300
## Paired donor-level cell-type proportion analysis
## 18 paired patients: Normal vs Tumor
##
## Significant cell types (FDR < 0.05):
##   plotted separately
##
## Non-significant cell types:
##   displayed together using the original facet layout
############################################################

rm(list = ls())
options(stringsAsFactors = FALSE)
Sys.setenv(LANGUAGE = "en")
set.seed(12345)

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

############################################################
## 1. Paths
############################################################

data_dir <- "C:/Users/92066/Desktop/单细胞孟德尔/Figure7/GSE294300"
setwd(data_dir)

input_file <- file.path(
  data_dir,
  "GSE294300_scRNA_final_annotated_with_tSNE.rds"
)

if (!file.exists(input_file)) {
  stop("没有找到：", input_file)
}

output_dir <- file.path(
  data_dir,
  "GSE294300_celltype_proportion_paired_results"
)

if (!dir.exists(output_dir)) {
  dir.create(
    output_dir,
    recursive = TRUE
  )
}

############################################################
## 2. Load annotated object
############################################################

scRNA <- readRDS(input_file)

required_meta <- c(
  "sample",
  "patient",
  "group",
  "celltype"
)

missing_meta <- setdiff(
  required_meta,
  colnames(scRNA@meta.data)
)

if (length(missing_meta) > 0) {
  stop(
    "缺少metadata：",
    paste(missing_meta, collapse = ", ")
  )
}

############################################################
## 3. Cell-type and group order
############################################################

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

scRNA$group <- factor(
  as.character(scRNA$group),
  levels = c("Normal", "Tumor")
)

scRNA$celltype <- factor(
  as.character(scRNA$celltype),
  levels = celltype_order
)

if (anyNA(scRNA$group)) {
  stop("group中存在非Normal/Tumor值。")
}

if (anyNA(scRNA$celltype)) {
  stop("celltype中存在未定义细胞类型。")
}

############################################################
## 4. Check paired study design
############################################################

sample_info <- scRNA@meta.data %>%
  dplyr::select(
    sample,
    patient,
    group
  ) %>%
  dplyr::distinct()

pair_check <- sample_info %>%
  dplyr::count(
    patient,
    group,
    name = "n_sample"
  ) %>%
  tidyr::pivot_wider(
    names_from = group,
    values_from = n_sample,
    values_fill = 0
  )

if (nrow(pair_check) != 18) {
  stop("患者数不是18。")
}

if (!all(
  pair_check$Normal == 1 &
  pair_check$Tumor == 1
)) {
  stop("不是所有患者都具有1个Normal和1个Tumor样本。")
}

############################################################
## 5. Calculate cell-type proportion for each sample
############################################################

sample_total <- scRNA@meta.data %>%
  dplyr::count(
    sample,
    name = "total_cells"
  )

celltype_count <- scRNA@meta.data %>%
  dplyr::count(
    sample,
    celltype,
    name = "n_cells"
  )

celltype_prop <- tidyr::expand_grid(
  sample = sample_info$sample,
  celltype = celltype_order
) %>%
  dplyr::left_join(
    celltype_count,
    by = c(
      "sample",
      "celltype"
    )
  ) %>%
  dplyr::mutate(
    n_cells = tidyr::replace_na(
      n_cells,
      0L
    )
  ) %>%
  dplyr::left_join(
    sample_total,
    by = "sample"
  ) %>%
  dplyr::left_join(
    sample_info,
    by = "sample"
  ) %>%
  dplyr::mutate(
    proportion = n_cells / total_cells,
    group = factor(
      group,
      levels = c(
        "Normal",
        "Tumor"
      )
    ),
    celltype = factor(
      celltype,
      levels = celltype_order
    )
  )

############################################################
## 6. Validate cell proportions
############################################################

prop_check <- celltype_prop %>%
  dplyr::group_by(sample) %>%
  dplyr::summarise(
    sum_prop = sum(proportion),
    .groups = "drop"
  )

if (
  max(
    abs(
      prop_check$sum_prop - 1
    )
  ) > 1e-10
) {
  stop(
    "至少一个样本的cell-type proportions之和不等于1。"
  )
}

write.csv(
  celltype_prop,
  file.path(
    output_dir,
    "GSE294300_celltype_proportion_by_sample.csv"
  ),
  row.names = FALSE
)

############################################################
## 7. Build paired Normal-Tumor table
############################################################

paired_prop <- celltype_prop %>%
  dplyr::select(
    patient,
    group,
    celltype,
    proportion
  ) %>%
  tidyr::pivot_wider(
    names_from = group,
    values_from = proportion
  )

if (
  !all(
    c(
      "Normal",
      "Tumor"
    ) %in% colnames(paired_prop)
  )
) {
  stop(
    "paired_prop缺少Normal或Tumor。"
  )
}

############################################################
## 8. Paired Wilcoxon test
############################################################

paired_wilcox <- function(
    normal,
    tumor
) {
  
  keep <- is.finite(normal) &
    is.finite(tumor)
  
  normal <- normal[keep]
  tumor <- tumor[keep]
  
  if (length(normal) < 2) {
    return(NA_real_)
  }
  
  if (
    all(
      (tumor - normal) == 0
    )
  ) {
    return(1)
  }
  
  result <- wilcox.test(
    tumor,
    normal,
    paired = TRUE,
    exact = FALSE,
    alternative = "two.sided"
  )
  
  result$p.value
}

############################################################
## 9. Statistics
############################################################

prop_stats <- paired_prop %>%
  dplyr::group_by(celltype) %>%
  dplyr::summarise(
    n_pairs = sum(
      is.finite(Normal) &
        is.finite(Tumor)
    ),
    median_Normal = median(
      Normal,
      na.rm = TRUE
    ),
    median_Tumor = median(
      Tumor,
      na.rm = TRUE
    ),
    median_delta = median(
      Tumor - Normal,
      na.rm = TRUE
    ),
    p_value = paired_wilcox(
      Normal,
      Tumor
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    FDR = p.adjust(
      p_value,
      method = "BH"
    ),
    significance = dplyr::case_when(
      is.na(FDR) ~ "NA",
      FDR < 0.001 ~ "***",
      FDR < 0.01 ~ "**",
      FDR < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    status = dplyr::case_when(
      is.na(FDR) ~ "NA",
      FDR < 0.05 ~ "Significant",
      TRUE ~ "ns"
    ),
    celltype = factor(
      celltype,
      levels = celltype_order
    )
  )

write.csv(
  prop_stats,
  file.path(
    output_dir,
    "GSE294300_celltype_proportion_paired_statistics.csv"
  ),
  row.names = FALSE
)

############################################################
## 10. Separate significant and non-significant results
############################################################

significant_stats <- prop_stats %>%
  dplyr::filter(
    !is.na(FDR),
    FDR < 0.05
  ) %>%
  dplyr::arrange(FDR)

nonsignificant_stats <- prop_stats %>%
  dplyr::filter(
    is.na(FDR) |
      FDR >= 0.05
  )

write.csv(
  significant_stats,
  file.path(
    output_dir,
    "Significant_celltypes_FDR_lt_0.05.csv"
  ),
  row.names = FALSE
)

write.csv(
  nonsignificant_stats,
  file.path(
    output_dir,
    "Non_significant_celltypes.csv"
  ),
  row.names = FALSE
)

cat("\n====================================================\n")
cat("Paired donor-level statistics\n")
cat("====================================================\n")

print(
  as.data.frame(prop_stats),
  row.names = FALSE
)

cat("\n====================================================\n")
cat("Significant cell types\n")
cat("====================================================\n")

print(
  as.data.frame(significant_stats),
  row.names = FALSE
)

############################################################
## 11. Plot data
############################################################

plot_data <- celltype_prop %>%
  dplyr::select(
    patient,
    group,
    celltype,
    proportion
  )

############################################################
## 12. Colors
############################################################

group_colors <- c(
  "Normal" = "#3C8DBC",
  "Tumor" = "#D95F59"
)

############################################################
## 13. Non-significant overview
## Keep original facet layout
############################################################

ns_celltypes <- as.character(
  nonsignificant_stats$celltype
)

plot_data_ns <- plot_data %>%
  dplyr::filter(
    as.character(celltype) %in%
      ns_celltypes
  ) %>%
  dplyr::mutate(
    celltype = factor(
      as.character(celltype),
      levels = celltype_order[
        celltype_order %in%
          ns_celltypes
      ]
    )
  )

############################################################
## 14. Annotation for non-significant overview
############################################################

annotation_ns <- plot_data_ns %>%
  dplyr::group_by(celltype) %>%
  dplyr::summarise(
    ymax = max(
      proportion,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  dplyr::left_join(
    nonsignificant_stats %>%
      dplyr::mutate(
        celltype = as.character(
          celltype
        )
      ) %>%
      dplyr::select(
        celltype,
        FDR,
        significance
      ),
    by = "celltype"
  ) %>%
  dplyr::mutate(
    y_position = dplyr::if_else(
      ymax == 0,
      0.01,
      ymax * 1.17
    ),
    label = dplyr::case_when(
      is.na(FDR) ~ "FDR = NA",
      TRUE ~ "ns"
    )
  )

############################################################
## 15. Plot non-significant cell types
############################################################

if (nrow(plot_data_ns) > 0) {
  
  p_ns <- ggplot(
    plot_data_ns,
    aes(
      x = group,
      y = proportion
    )
  ) +
    
    geom_line(
      aes(
        group = patient
      ),
      color = "grey70",
      linewidth = 0.65,
      alpha = 0.85
    ) +
    
    geom_point(
      aes(
        color = group
      ),
      size = 3.4,
      alpha = 0.95
    ) +
    
    stat_summary(
      aes(
        color = group
      ),
      fun = median,
      geom = "crossbar",
      width = 0.28,
      linewidth = 0.9
    ) +
    
    geom_text(
      data = annotation_ns,
      aes(
        x = 1.5,
        y = y_position,
        label = label
      ),
      inherit.aes = FALSE,
      size = 5.4,
      fontface = "bold",
      color = "black"
    ) +
    
    facet_wrap(
      ~ celltype,
      scales = "free_y",
      ncol = 5
    ) +
    
    scale_color_manual(
      values = group_colors,
      drop = FALSE
    ) +
    
    scale_y_continuous(
      labels = scales::percent_format(
        accuracy = 0.1
      ),
      expand = expansion(
        mult = c(
          0.04,
          0.22
        )
      )
    ) +
    
    labs(
      x = NULL,
      y = "Cell proportion"
    ) +
    
    theme_classic(
      base_size = 16
    ) +
    
    theme(
      strip.background = element_blank(),
      
      strip.text = element_text(
        size = 16,
        face = "bold",
        color = "black",
        margin = margin(
          b = 7
        )
      ),
      
      axis.title.y = element_text(
        size = 17,
        color = "black",
        margin = margin(
          r = 8
        )
      ),
      
      axis.text.x = element_text(
        size = 14,
        color = "black"
      ),
      
      axis.text.y = element_text(
        size = 13,
        color = "black"
      ),
      
      axis.line = element_line(
        linewidth = 0.75,
        color = "black"
      ),
      
      axis.ticks = element_line(
        linewidth = 0.65,
        color = "black"
      ),
      
      axis.ticks.length = unit(
        0.15,
        "cm"
      ),
      
      legend.position = "none",
      
      panel.grid = element_blank(),
      
      panel.spacing.x = unit(
        1.0,
        "cm"
      ),
      
      panel.spacing.y = unit(
        1.1,
        "cm"
      ),
      
      plot.margin = margin(
        15,
        20,
        15,
        20
      )
    ) +
    
    coord_cartesian(
      clip = "off"
    )
  
  ggsave(
    file.path(
      output_dir,
      "GSE294300_NS_celltype_proportions_overview.pdf"
    ),
    p_ns,
    width = 16,
    height = 8.5,
    device = cairo_pdf
  )
  
  ggsave(
    file.path(
      output_dir,
      "GSE294300_NS_celltype_proportions_overview.png"
    ),
    p_ns,
    width = 16,
    height = 8.5,
    dpi = 600,
    bg = "white"
  )
}

############################################################
## 16. Function for significant cell types
## Each significant cell type is plotted separately
############################################################

plot_significant_celltype <- function(
    celltype_name
) {
  
  plot_df <- plot_data %>%
    dplyr::filter(
      as.character(celltype) ==
        celltype_name
    )
  
  stat_df <- prop_stats %>%
    dplyr::filter(
      as.character(celltype) ==
        celltype_name
    )
  
  if (nrow(stat_df) != 1) {
    stop(
      "无法唯一匹配统计结果：",
      celltype_name
    )
  }
  
  ymax <- max(
    plot_df$proportion,
    na.rm = TRUE
  )
  
  y_position <- ifelse(
    ymax == 0,
    0.01,
    ymax * 1.18
  )
  
  fdr_label <- dplyr::case_when(
    stat_df$FDR < 0.001 ~
      "FDR < 0.001",
    
    TRUE ~
      paste0(
        "FDR = ",
        formatC(
          stat_df$FDR,
          format = "fg",
          digits = 3
        )
      )
  )
  
  ggplot(
    plot_df,
    aes(
      x = group,
      y = proportion
    )
  ) +
    
    geom_line(
      aes(
        group = patient
      ),
      color = "grey68",
      linewidth = 0.75,
      alpha = 0.85
    ) +
    
    geom_point(
      aes(
        color = group
      ),
      size = 4.5,
      alpha = 0.95
    ) +
    
    stat_summary(
      aes(
        color = group
      ),
      fun = median,
      geom = "crossbar",
      width = 0.30,
      linewidth = 1.1
    ) +
    
    annotate(
      "text",
      x = 1.5,
      y = y_position,
      label = paste0(
        fdr_label,
        "\n",
        "n = ",
        stat_df$n_pairs,
        " paired patients"
      ),
      size = 5.6,
      lineheight = 1.05
    ) +
    
    scale_color_manual(
      values = group_colors,
      drop = FALSE
    ) +
    
    scale_y_continuous(
      labels = scales::percent_format(
        accuracy = 0.1
      ),
      expand = expansion(
        mult = c(
          0.05,
          0.28
        )
      )
    ) +
    
    labs(
      x = NULL,
      y = "Cell proportion",
      title = celltype_name
    ) +
    
    theme_classic(
      base_size = 17
    ) +
    
    theme(
      plot.title = element_text(
        size = 21,
        face = "bold",
        hjust = 0.5,
        margin = margin(
          b = 10
        )
      ),
      
      axis.title.y = element_text(
        size = 18,
        color = "black",
        margin = margin(
          r = 8
        )
      ),
      
      axis.text.x = element_text(
        size = 17,
        color = "black"
      ),
      
      axis.text.y = element_text(
        size = 14,
        color = "black"
      ),
      
      axis.line = element_line(
        linewidth = 0.8,
        color = "black"
      ),
      
      axis.ticks = element_line(
        linewidth = 0.7,
        color = "black"
      ),
      
      axis.ticks.length = unit(
        0.16,
        "cm"
      ),
      
      legend.position = "none",
      
      panel.grid = element_blank(),
      
      plot.margin = margin(
        20,
        30,
        20,
        25
      )
    ) +
    
    coord_cartesian(
      clip = "off"
    )
}

############################################################
## 17. Automatically plot significant cell types
############################################################

if (nrow(significant_stats) > 0) {
  
  for (
    i in seq_len(
      nrow(significant_stats)
    )
  ) {
    
    celltype_i <- as.character(
      significant_stats$celltype[i]
    )
    
    p_sig <- plot_significant_celltype(
      celltype_i
    )
    
    safe_name <- gsub(
      "[^A-Za-z0-9]+",
      "_",
      celltype_i
    )
    
    ggsave(
      file.path(
        output_dir,
        paste0(
          "SIGNIFICANT_",
          safe_name,
          "_paired_proportion.pdf"
        )
      ),
      p_sig,
      width = 6.5,
      height = 7,
      device = cairo_pdf
    )
    
    ggsave(
      file.path(
        output_dir,
        paste0(
          "SIGNIFICANT_",
          safe_name,
          "_paired_proportion.png"
        )
      ),
      p_sig,
      width = 6.5,
      height = 7,
      dpi = 600,
      bg = "white"
    )
  }
}

############################################################
## 18. Optional complete 10-cell-type plot
## Keep one complete plot for record/checking
############################################################

annotation_all <- plot_data %>%
  dplyr::group_by(celltype) %>%
  dplyr::summarise(
    ymax = max(
      proportion,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  dplyr::left_join(
    prop_stats %>%
      dplyr::select(
        celltype,
        FDR,
        significance
      ),
    by = "celltype"
  ) %>%
  dplyr::mutate(
    y_position = dplyr::if_else(
      ymax == 0,
      0.01,
      ymax * 1.17
    ),
    
    label = dplyr::case_when(
      is.na(FDR) ~ "NA",
      FDR < 0.001 ~ "***",
      FDR < 0.01 ~ "**",
      FDR < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

p_all <- ggplot(
  plot_data,
  aes(
    x = group,
    y = proportion
  )
) +
  
  geom_line(
    aes(
      group = patient
    ),
    color = "grey70",
    linewidth = 0.65,
    alpha = 0.85
  ) +
  
  geom_point(
    aes(
      color = group
    ),
    size = 3.4,
    alpha = 0.95
  ) +
  
  stat_summary(
    aes(
      color = group
    ),
    fun = median,
    geom = "crossbar",
    width = 0.28,
    linewidth = 0.9
  ) +
  
  geom_text(
    data = annotation_all,
    aes(
      x = 1.5,
      y = y_position,
      label = label
    ),
    inherit.aes = FALSE,
    size = 5.4,
    fontface = "bold"
  ) +
  
  facet_wrap(
    ~ celltype,
    scales = "free_y",
    ncol = 5
  ) +
  
  scale_color_manual(
    values = group_colors,
    drop = FALSE
  ) +
  
  scale_y_continuous(
    labels = scales::percent_format(
      accuracy = 0.1
    ),
    expand = expansion(
      mult = c(
        0.04,
        0.22
      )
    )
  ) +
  
  labs(
    x = NULL,
    y = "Cell proportion"
  ) +
  
  theme_classic(
    base_size = 16
  ) +
  
  theme(
    strip.background = element_blank(),
    
    strip.text = element_text(
      size = 16,
      face = "bold",
      color = "black"
    ),
    
    axis.title.y = element_text(
      size = 17
    ),
    
    axis.text.x = element_text(
      size = 14,
      color = "black"
    ),
    
    axis.text.y = element_text(
      size = 13,
      color = "black"
    ),
    
    axis.line = element_line(
      linewidth = 0.75,
      color = "black"
    ),
    
    axis.ticks = element_line(
      linewidth = 0.65,
      color = "black"
    ),
    
    legend.position = "none",
    
    panel.grid = element_blank(),
    
    panel.spacing.x = unit(
      1.0,
      "cm"
    ),
    
    panel.spacing.y = unit(
      1.1,
      "cm"
    )
  ) +
  
  coord_cartesian(
    clip = "off"
  )

ggsave(
  file.path(
    output_dir,
    "GSE294300_ALL_10_celltypes_paired_proportions.pdf"
  ),
  p_all,
  width = 16,
  height = 8.5,
  device = cairo_pdf
)

ggsave(
  file.path(
    output_dir,
    "GSE294300_ALL_10_celltypes_paired_proportions.png"
  ),
  p_all,
  width = 16,
  height = 8.5,
  dpi = 600,
  bg = "white"
)

############################################################
## 19. Session information
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
## 20. Final output
############################################################

cat("\n############################################################\n")
cat("## CELL-TYPE PROPORTION ANALYSIS COMPLETED\n")
cat("############################################################\n")

cat("\nOutput folder:\n")
cat(output_dir, "\n")

cat("\nNumber of significant cell types:\n")
cat(nrow(significant_stats), "\n")

if (nrow(significant_stats) > 0) {
  
  cat("\nSignificant cell types:\n")
  
  print(
    as.data.frame(
      significant_stats
    ),
    row.names = FALSE
  )
}

cat("\nMain outputs:\n")
cat("GSE294300_NS_celltype_proportions_overview.png\n")
cat("GSE294300_ALL_10_celltypes_paired_proportions.png\n")
cat("Significant cell types are saved separately with prefix: SIGNIFICANT_\n")