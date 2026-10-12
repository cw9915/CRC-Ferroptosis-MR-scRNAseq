############################################################
## Figure 6D only
## Core lncRNA-miRNA-mRNA regulatory network
## Genes: AKR1C2, MAPK3, MAP3K11
############################################################

Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list = ls())

############################################################
## 0. 加载 R 包
############################################################

pkg_need <- c(
  "dplyr",
  "data.table",
  "igraph",
  "ggraph",
  "ggplot2",
  "ggrepel"
)

for (pkg in pkg_need) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("缺少 R 包：", pkg, "。请先安装后再运行。")
  }
}

library(dplyr)
library(data.table)
library(igraph)
library(ggraph)
library(ggplot2)
library(ggrepel)

## 避免常用函数被其他包覆盖
select     <- dplyr::select
filter     <- dplyr::filter
mutate     <- dplyr::mutate
arrange    <- dplyr::arrange
summarise  <- dplyr::summarise
group_by   <- dplyr::group_by
distinct   <- dplyr::distinct
left_join  <- dplyr::left_join
slice_head <- dplyr::slice_head

############################################################
## 1. 设置路径和目标基因
############################################################

fig6_dir <- "C:/Users/92066/Desktop/单细胞孟德尔/Figure6"

miRNA_mRNA_file <- file.path(fig6_dir, "hsa_MTI.csv")
lncRNA_miRNA_file <- file.path(fig6_dir, "spongeScan.tsv")

gene <- c(
  "AKR1C2",
  "MAPK3",
  "MAP3K11"
)

## 核心网络参数
top_miRNA_per_gene <- 6
top_lnc_per_miRNA <- 4

if (!dir.exists(fig6_dir)) {
  stop("Figure6 文件夹不存在：", fig6_dir)
}

if (!file.exists(miRNA_mRNA_file)) {
  stop("找不到 hsa_MTI.csv：", miRNA_mRNA_file)
}

if (!file.exists(lncRNA_miRNA_file)) {
  stop("找不到 spongeScan.tsv：", lncRNA_miRNA_file)
}

setwd(fig6_dir)

############################################################
## 2. 读取 miRNA-mRNA 文件
############################################################

miTB <- data.table::fread(
  miRNA_mRNA_file,
  data.table = FALSE
)

cat("hsa_MTI.csv 列名如下：\n")
print(colnames(miTB))

miRNA_col_candidates <- c(
  "miRNA",
  "miRNA Name",
  "miRNA_name",
  "Mirna",
  "mirna",
  "miRNA_ID",
  "miRNA ID",
  "miRNA.id"
)

target_col_candidates <- c(
  "Target Gene",
  "Target_Gene",
  "TargetGene",
  "Target Symbol",
  "target",
  "Target",
  "Gene",
  "gene",
  "mRNA",
  "Symbol"
)

miRNA_hit <- intersect(
  miRNA_col_candidates,
  colnames(miTB)
)

target_hit <- intersect(
  target_col_candidates,
  colnames(miTB)
)

if (length(miRNA_hit) == 0) {
  stop(
    "无法识别 hsa_MTI.csv 中的 miRNA 列。\n当前列名：",
    paste(colnames(miTB), collapse = ", ")
  )
}

if (length(target_hit) == 0) {
  stop(
    "无法识别 hsa_MTI.csv 中的 Target Gene 列。\n当前列名：",
    paste(colnames(miTB), collapse = ", ")
  )
}

miRNA_col <- miRNA_hit[1]
target_col <- target_hit[1]

cat("识别到 miRNA 列：", miRNA_col, "\n")
cat("识别到 Target Gene 列：", target_col, "\n")

## 使用基础 R 改名，避免 all_of 等函数兼容性问题
miTB2 <- miTB

colnames(miTB2)[
  match(miRNA_col, colnames(miTB2))
] <- "miRNA"

colnames(miTB2)[
  match(target_col, colnames(miTB2))
] <- "Gene"

miTB2 <- miTB2 %>%
  dplyr::mutate(
    miRNA = trimws(as.character(miRNA)),
    Gene = trimws(as.character(Gene))
  ) %>%
  dplyr::filter(
    !is.na(miRNA),
    miRNA != "",
    !is.na(Gene),
    Gene != ""
  )

############################################################
## 3. 如果存在证据列，优先保留强证据
############################################################

support_col_candidates <- c(
  "Support Type",
  "Support_Type",
  "support_type",
  "Experiments",
  "Experiment",
  "Validation",
  "validation"
)

support_hit <- intersect(
  support_col_candidates,
  colnames(miTB2)
)

if (length(support_hit) > 0) {
  
  support_col <- support_hit[1]
  
  cat(
    "检测到证据列：",
    support_col,
    "\n"
  )
  
  cat(
    "该列前几个取值：\n"
  )
  
  print(
    head(
      unique(miTB2[[support_col]]),
      20
    )
  )
  
  miTB_strong <- miTB2 %>%
    dplyr::filter(
      grepl(
        "Functional",
        .data[[support_col]],
        ignore.case = TRUE
      ) |
        grepl(
          "Strong",
          .data[[support_col]],
          ignore.case = TRUE
        ) |
        grepl(
          "Reporter|Western|qPCR|Luciferase|CLIP|Sequencing",
          .data[[support_col]],
          ignore.case = TRUE
        )
    )
  
  if (nrow(miTB_strong) > 0) {
    
    miTB2 <- miTB_strong
    
    cat(
      "已优先保留强证据 miRNA-mRNA 关系，剩余行数：",
      nrow(miTB2),
      "\n"
    )
    
  } else {
    
    cat(
      "未筛到强证据关系，继续保留全部 miRNA-mRNA 关系。\n"
    )
  }
  
} else {
  
  cat(
    "未检测到证据列，直接使用全部 miRNA-mRNA 关系。\n"
  )
}

############################################################
## 4. 筛选 3 个目标基因相关 miRNA
############################################################

miTB_sub <- miTB2 %>%
  dplyr::filter(
    Gene %in% gene
  ) %>%
  dplyr::select(
    miRNA,
    Gene
  ) %>%
  dplyr::distinct()

if (nrow(miTB_sub) == 0) {
  stop(
    "hsa_MTI.csv 中没有筛到目标基因对应的 miRNA。"
  )
}

cat(
  "筛选到 miRNA-mRNA 关系数：",
  nrow(miTB_sub),
  "\n"
)

cat(
  "各目标基因对应关系数：\n"
)

print(
  table(miTB_sub$Gene)
)

missing_gene <- setdiff(
  gene,
  unique(miTB_sub$Gene)
)

if (length(missing_gene) > 0) {
  warning(
    "以下基因没有筛到 miRNA：",
    paste(
      missing_gene,
      collapse = ", "
    )
  )
}

write.csv(
  miTB_sub,
  file.path(
    fig6_dir,
    "Figure6D_miRNA_mRNA_all_edges.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)

############################################################
## 5. 读取 lncRNA-miRNA 文件
############################################################

sponge_df <- data.table::fread(
  lncRNA_miRNA_file,
  data.table = FALSE
)

cat(
  "spongeScan.tsv 列名如下：\n"
)

print(
  colnames(sponge_df)
)

lnc_col_candidates <- c(
  "lncRNA",
  "LncRNA",
  "lncRNA_symbol",
  "lncRNA_name",
  "lncRNA Name",
  "lncrna",
  "Gene",
  "gene"
)

sponge_miRNA_col_candidates <- c(
  "miRNA",
  "miRNAs",
  "miRNA_name",
  "miRNA Name",
  "mirna",
  "Mirna",
  "miRNA_ID",
  "miRNA ID"
)

lnc_hit <- intersect(
  lnc_col_candidates,
  colnames(sponge_df)
)

sponge_miRNA_hit <- intersect(
  sponge_miRNA_col_candidates,
  colnames(sponge_df)
)

if (length(lnc_hit) == 0) {
  stop(
    "无法识别 spongeScan.tsv 中的 lncRNA 列。\n当前列名：",
    paste(
      colnames(sponge_df),
      collapse = ", "
    )
  )
}

if (length(sponge_miRNA_hit) == 0) {
  stop(
    "无法识别 spongeScan.tsv 中的 miRNA 列。\n当前列名：",
    paste(
      colnames(sponge_df),
      collapse = ", "
    )
  )
}

lnc_col <- lnc_hit[1]
sponge_miRNA_col <- sponge_miRNA_hit[1]

cat(
  "识别到 lncRNA 列：",
  lnc_col,
  "\n"
)

cat(
  "识别到 sponge miRNA 列：",
  sponge_miRNA_col,
  "\n"
)

sponge2 <- sponge_df

colnames(sponge2)[
  match(
    lnc_col,
    colnames(sponge2)
  )
] <- "lncRNA"

colnames(sponge2)[
  match(
    sponge_miRNA_col,
    colnames(sponge2)
  )
] <- "miRNA"

sponge2 <- sponge2 %>%
  dplyr::mutate(
    lncRNA = trimws(
      as.character(lncRNA)
    ),
    miRNA = trimws(
      as.character(miRNA)
    )
  ) %>%
  dplyr::filter(
    !is.na(lncRNA),
    lncRNA != "",
    !is.na(miRNA),
    miRNA != ""
  )

############################################################
## 6. 筛选目标 miRNA 对应的 lncRNA
############################################################

miRNA_list <- unique(
  miTB_sub$miRNA
)

lnc_mi_sub <- sponge2 %>%
  dplyr::filter(
    miRNA %in% miRNA_list
  ) %>%
  dplyr::select(
    lncRNA,
    miRNA
  ) %>%
  dplyr::distinct()

if (nrow(lnc_mi_sub) == 0) {
  stop(
    "spongeScan.tsv 中没有筛到目标 miRNA 对应的 lncRNA。"
  )
}

cat(
  "筛选到 lncRNA-miRNA 关系数：",
  nrow(lnc_mi_sub),
  "\n"
)

cat(
  "涉及 miRNA 数：",
  length(
    unique(
      lnc_mi_sub$miRNA
    )
  ),
  "\n"
)

cat(
  "涉及 lncRNA 数：",
  length(
    unique(
      lnc_mi_sub$lncRNA
    )
  ),
  "\n"
)

write.csv(
  lnc_mi_sub,
  file.path(
    fig6_dir,
    "Figure6D_lncRNA_miRNA_all_edges.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)

############################################################
## 7. 构建完整网络边表
############################################################

edges_lnc_mi <- lnc_mi_sub %>%
  dplyr::transmute(
    from = lncRNA,
    to = miRNA,
    edge_type = "lncRNA-miRNA"
  ) %>%
  dplyr::distinct()

edges_mi_gene <- miTB_sub %>%
  dplyr::transmute(
    from = miRNA,
    to = Gene,
    edge_type = "miRNA-mRNA"
  ) %>%
  dplyr::distinct()

edges_all <- dplyr::bind_rows(
  edges_lnc_mi,
  edges_mi_gene
) %>%
  dplyr::distinct()

nodes_all <- data.frame(
  name = unique(
    c(
      edges_all$from,
      edges_all$to
    )
  ),
  stringsAsFactors = FALSE
)

nodes_all$layer <- dplyr::case_when(
  nodes_all$name %in% gene ~ "Gene",
  grepl(
    "^hsa-|^miR|^let",
    nodes_all$name,
    ignore.case = TRUE
  ) ~ "miRNA",
  TRUE ~ "lncRNA"
)

write.csv(
  edges_all,
  file.path(
    fig6_dir,
    "Figure6D_full_network_edges.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  nodes_all,
  file.path(
    fig6_dir,
    "Figure6D_full_network_nodes.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)

cat(
  "完整网络节点数：",
  nrow(nodes_all),
  "\n"
)

cat(
  "完整网络边数：",
  nrow(edges_all),
  "\n"
)

############################################################
## 8. 提取核心网络
############################################################

## 每个 Gene 最多保留 6 个 miRNA
mi_gene_core <- edges_mi_gene %>%
  dplyr::group_by(to) %>%
  dplyr::slice_head(
    n = top_miRNA_per_gene
  ) %>%
  dplyr::ungroup()

core_miRNAs <- unique(
  mi_gene_core$from
)

## 每个 miRNA 最多保留 4 个 lncRNA
lnc_mi_core <- edges_lnc_mi %>%
  dplyr::filter(
    to %in% core_miRNAs
  ) %>%
  dplyr::group_by(to) %>%
  dplyr::slice_head(
    n = top_lnc_per_miRNA
  ) %>%
  dplyr::ungroup()

edges_core <- dplyr::bind_rows(
  lnc_mi_core,
  mi_gene_core
) %>%
  dplyr::distinct()

nodes_core <- data.frame(
  name = unique(
    c(
      edges_core$from,
      edges_core$to
    )
  ),
  stringsAsFactors = FALSE
)

nodes_core$layer <- dplyr::case_when(
  nodes_core$name %in% gene ~ "Gene",
  grepl(
    "^hsa-|^miR|^let",
    nodes_core$name,
    ignore.case = TRUE
  ) ~ "miRNA",
  TRUE ~ "lncRNA"
)

############################################################
## 9. 如果核心网络太小，自动增加 lncRNA
############################################################

if (nrow(nodes_core) < 15) {
  
  cat(
    "核心网络节点少于 15，自动将每个 miRNA 的 lncRNA 数放宽到 6。\n"
  )
  
  top_lnc_per_miRNA <- 6
  
  lnc_mi_core <- edges_lnc_mi %>%
    dplyr::filter(
      to %in% core_miRNAs
    ) %>%
    dplyr::group_by(to) %>%
    dplyr::slice_head(
      n = top_lnc_per_miRNA
    ) %>%
    dplyr::ungroup()
  
  edges_core <- dplyr::bind_rows(
    lnc_mi_core,
    mi_gene_core
  ) %>%
    dplyr::distinct()
  
  nodes_core <- data.frame(
    name = unique(
      c(
        edges_core$from,
        edges_core$to
      )
    ),
    stringsAsFactors = FALSE
  )
  
  nodes_core$layer <- dplyr::case_when(
    nodes_core$name %in% gene ~ "Gene",
    grepl(
      "^hsa-|^miR|^let",
      nodes_core$name,
      ignore.case = TRUE
    ) ~ "miRNA",
    TRUE ~ "lncRNA"
  )
}

nodes_core$node_size <- dplyr::case_when(
  nodes_core$layer == "Gene" ~ 8,
  nodes_core$layer == "miRNA" ~ 6,
  TRUE ~ 4
)

## 在最终核心网络确定后再保存
write.csv(
  edges_core,
  file.path(
    fig6_dir,
    "Figure6D_core_network_edges.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  nodes_core,
  file.path(
    fig6_dir,
    "Figure6D_core_network_nodes.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)

cat(
  "最终核心网络节点数：",
  nrow(nodes_core),
  "\n"
)

cat(
  "最终核心网络边数：",
  nrow(edges_core),
  "\n"
)

cat(
  "核心网络节点分布：\n"
)

print(
  table(nodes_core$layer)
)

############################################################
## 10. 构建 igraph 对象
############################################################

if (
  !all(
    c(
      "from",
      "to",
      "edge_type"
    ) %in% colnames(edges_core)
  )
) {
  stop(
    "edges_core 缺少 from、to 或 edge_type 列。"
  )
}

g_core <- igraph::graph_from_data_frame(
  d = edges_core,
  vertices = nodes_core,
  directed = TRUE
)

## 删除孤立节点
g_core <- igraph::delete_vertices(
  g_core,
  igraph::V(g_core)[
    igraph::degree(
      g_core,
      mode = "all"
    ) == 0
  ]
)

## 检查 edge_type 是否正确进入 igraph
if (
  !"edge_type" %in%
  igraph::edge_attr_names(g_core)
) {
  stop(
    "g_core 中没有 edge_type 属性，请检查 edges_core。"
  )
}

## 保持 nodes_core 与图顶点顺序完全一致
nodes_core <- nodes_core[
  match(
    igraph::V(g_core)$name,
    nodes_core$name
  ),
  ,
  drop = FALSE
]

############################################################
## 11. 生成分层环形布局坐标
## Gene 内圈 / miRNA 中圈 / lncRNA 外圈
############################################################

## 在原始代码基础上压缩半径
r_gene <- 1.05
r_miRNA <- 2.15
r_lnc <- 3.55

## 三个基因：
## AKR1C2 顶部
## MAPK3 右下
## MAP3K11 左下
gene_order <- c(
  "AKR1C2",
  "MAPK3",
  "MAP3K11"
)

gene_present <- gene_order[
  gene_order %in%
    nodes_core$name
]

n_gene <- length(
  gene_present
)

if (n_gene == 0) {
  stop(
    "核心网络中没有检测到目标 Gene 节点。"
  )
}

## 从顶部开始顺时针均匀排列
gene_angles <- seq(
  from = pi / 2,
  to = pi / 2 -
    2 * pi +
    2 * pi / n_gene,
  length.out = n_gene
)

gene_pos <- data.frame(
  name = gene_present,
  layer = "Gene",
  angle = gene_angles,
  radius = r_gene,
  stringsAsFactors = FALSE
)

############################################################
## 12. 为 miRNA 分配位置
############################################################

mi_gene_map <- edges_core %>%
  dplyr::filter(
    edge_type == "miRNA-mRNA"
  ) %>%
  dplyr::transmute(
    miRNA = from,
    Gene = to
  ) %>%
  dplyr::distinct()

## 保留原代码逻辑：
## miRNA 同时连接多个 Gene 时，
## 取第一个 Gene 作为绘图归属
mi_gene_assign <- mi_gene_map %>%
  dplyr::group_by(miRNA) %>%
  dplyr::summarise(
    Gene = Gene[1],
    .groups = "drop"
  )

miRNA_pos_list <- list()

## 根据实际 Gene 数自动确定扇区
## 三个 Gene 时约 98°
sector_width <- 0.82 * (
  2 * pi / n_gene
)

for (g in gene_present) {
  
  mis <- mi_gene_assign$miRNA[
    mi_gene_assign$Gene == g
  ]
  
  mis <- unique(mis)
  
  if (length(mis) == 0) {
    next
  }
  
  center_angle <- gene_pos$angle[
    gene_pos$name == g
  ]
  
  if (length(mis) == 1) {
    
    angs <- center_angle
    
  } else {
    
    angs <- seq(
      from =
        center_angle -
        sector_width / 2,
      to =
        center_angle +
        sector_width / 2,
      length.out =
        length(mis)
    )
  }
  
  miRNA_pos_list[[g]] <- data.frame(
    name = mis,
    layer = "miRNA",
    angle = angs,
    radius = r_miRNA,
    stringsAsFactors = FALSE
  )
}

if (length(miRNA_pos_list) == 0) {
  stop(
    "没有成功生成 miRNA 坐标。"
  )
}

miRNA_pos <- dplyr::bind_rows(
  miRNA_pos_list
)

############################################################
## 13. 为 lncRNA 分配位置
############################################################

lnc_mi_map <- edges_core %>%
  dplyr::filter(
    edge_type == "lncRNA-miRNA"
  ) %>%
  dplyr::transmute(
    lncRNA = from,
    miRNA = to
  ) %>%
  dplyr::distinct()

lnc_pos_list <- list()

## 每个 miRNA 周围的小扇区
lnc_width <- pi / 11

for (mi in miRNA_pos$name) {
  
  lncs <- lnc_mi_map$lncRNA[
    lnc_mi_map$miRNA == mi
  ]
  
  lncs <- unique(lncs)
  
  if (length(lncs) == 0) {
    next
  }
  
  center_angle <- miRNA_pos$angle[
    miRNA_pos$name == mi
  ]
  
  if (length(lncs) == 1) {
    
    angs <- center_angle
    
  } else {
    
    angs <- seq(
      from =
        center_angle -
        lnc_width / 2,
      to =
        center_angle +
        lnc_width / 2,
      length.out =
        length(lncs)
    )
  }
  
  lnc_pos_list[[mi]] <- data.frame(
    name = lncs,
    layer = "lncRNA",
    angle = angs,
    radius = r_lnc,
    stringsAsFactors = FALSE
  )
}

if (length(lnc_pos_list) > 0) {
  
  lnc_pos <- dplyr::bind_rows(
    lnc_pos_list
  )
  
} else {
  
  lnc_pos <- data.frame(
    name = character(0),
    layer = character(0),
    angle = numeric(0),
    radius = numeric(0),
    stringsAsFactors = FALSE
  )
}

############################################################
## 14. 合并所有节点坐标
############################################################

pos_df <- dplyr::bind_rows(
  gene_pos,
  miRNA_pos,
  lnc_pos
) %>%
  dplyr::distinct(
    name,
    .keep_all = TRUE
  )

## 防止少量节点没有被自动分配坐标
missing_nodes <- setdiff(
  nodes_core$name,
  pos_df$name
)

if (length(missing_nodes) > 0) {
  
  miss_layer <- nodes_core$layer[
    match(
      missing_nodes,
      nodes_core$name
    )
  ]
  
  miss_angle <- seq(
    0,
    2 * pi,
    length.out =
      length(missing_nodes) + 1
  )[
    seq_along(missing_nodes)
  ]
  
  miss_radius <- ifelse(
    miss_layer == "Gene",
    r_gene,
    ifelse(
      miss_layer == "miRNA",
      r_miRNA,
      r_lnc
    )
  )
  
  miss_df <- data.frame(
    name = missing_nodes,
    layer = miss_layer,
    angle = miss_angle,
    radius = miss_radius,
    stringsAsFactors = FALSE
  )
  
  pos_df <- dplyr::bind_rows(
    pos_df,
    miss_df
  )
}

## 极坐标转换成 XY
pos_df <- pos_df %>%
  dplyr::mutate(
    x = radius * cos(angle),
    y = radius * sin(angle)
  )

## 合并回 nodes_core
nodes_core <- nodes_core %>%
  dplyr::select(
    -dplyr::any_of(
      c(
        "angle",
        "radius",
        "x",
        "y"
      )
    )
  ) %>%
  dplyr::left_join(
    pos_df %>%
      dplyr::select(
        name,
        angle,
        radius,
        x,
        y
      ),
    by = "name"
  )

nodes_core <- nodes_core[
  match(
    igraph::V(g_core)$name,
    nodes_core$name
  ),
  ,
  drop = FALSE
]

if (
  anyNA(nodes_core$x) ||
  anyNA(nodes_core$y)
) {
  stop(
    "部分网络节点没有成功获得坐标。"
  )
}

############################################################
## 15. 创建手动布局
############################################################

layout_manual <- ggraph::create_layout(
  g_core,
  layout = "manual",
  x = nodes_core$x,
  y = nodes_core$y
)

############################################################
## 16. 绘制 Figure 6D
############################################################

set.seed(2026)

p_network <- ggraph::ggraph(
  layout_manual
) +
  
  ## 网络边
  ## lncRNA-miRNA = 实线
  ## miRNA-mRNA   = 虚线
  ggraph::geom_edge_link(
    ggplot2::aes(
      linetype = edge_type
    ),
    color = "grey75",
    linewidth = 0.60,
    alpha = 0.80,
    arrow = grid::arrow(
      length = grid::unit(
        2.1,
        "mm"
      )
    ),
    end_cap = circle(
      2.8,
      "mm"
    ),
    show.legend = FALSE
  ) +
  
  ## 节点
  ggraph::geom_node_point(
    ggplot2::aes(
      shape = layer,
      fill = layer,
      size = layer
    ),
    color = "white",
    stroke = 0.9,
    alpha = 0.99
  ) +
  
  ## Gene 标签
  ggraph::geom_node_text(
    data = subset(
      layout_manual,
      layer == "Gene"
    ),
    ggplot2::aes(
      x = x,
      y = y,
      label = name
    ),
    color = "#E64B35",
    size = 5.8,
    fontface = "bold",
    repel = TRUE,
    max.overlaps = Inf,
    box.padding = 0.20,
    point.padding = 0.15
  ) +
  
  ## miRNA 标签
  ggraph::geom_node_text(
    data = subset(
      layout_manual,
      layer == "miRNA"
    ),
    ggplot2::aes(
      x = x,
      y = y,
      label = name
    ),
    color = "#4DBBD5",
    size = 4.5,
    fontface = "bold",
    repel = TRUE,
    max.overlaps = Inf,
    box.padding = 0.12,
    point.padding = 0.08
  ) +
  
  ## lncRNA 标签
  ggraph::geom_node_text(
    data = subset(
      layout_manual,
      layer == "lncRNA"
    ),
    ggplot2::aes(
      x = x,
      y = y,
      label = name
    ),
    color = "#00A087",
    size = 4.0,
    fontface = "bold",
    repel = TRUE,
    max.overlaps = Inf,
    box.padding = 0.10,
    point.padding = 0.06
  ) +
  
  ## 节点形状
  ggplot2::scale_shape_manual(
    values = c(
      "Gene" = 21,
      "miRNA" = 21,
      "lncRNA" = 22
    )
  ) +
  
  ## 节点颜色
  ggplot2::scale_fill_manual(
    values = c(
      "Gene" = "#E64B35",
      "miRNA" = "#4DBBD5",
      "lncRNA" = "#00A087"
    )
  ) +
  
  ## 节点大小
  ggplot2::scale_size_manual(
    values = c(
      "Gene" = 15,
      "miRNA" = 8.0,
      "lncRNA" = 5.3
    )
  ) +
  
  ## 边线型
  ggplot2::scale_linetype_manual(
    values = c(
      "lncRNA-miRNA" = "solid",
      "miRNA-mRNA" = "22"
    )
  ) +
  
  ## 固定横纵比例，避免网络被拉成长椭圆
  ggplot2::coord_fixed(
    ratio = 1,
    clip = "off"
  ) +
  
  ## 减少外围空白
  ggplot2::scale_x_continuous(
    expand = ggplot2::expansion(
      mult = 0.08
    )
  ) +
  
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(
      mult = 0.08
    )
  ) +
  
  ggplot2::theme_void() +
  
  ggplot2::theme(
    legend.position = "right",
    
    legend.title = ggplot2::element_text(
      size = 15,
      face = "bold"
    ),
    
    legend.text = ggplot2::element_text(
      size = 14
    ),
    
    legend.key.size = grid::unit(
      0.65,
      "cm"
    ),
    
    plot.margin = ggplot2::margin(
      t = 3,
      r = 3,
      b = 3,
      l = 3
    )
  ) +
  
  ggplot2::labs(
    title = NULL,
    fill = "Layer",
    shape = "Layer",
    size = "Layer"
  )

print(p_network)

############################################################
## 17. 保存 Figure 6D
############################################################

ggplot2::ggsave(
  filename = file.path(
    fig6_dir,
    "Figure6D_core_lncRNA_miRNA_mRNA_network_circular_compact.pdf"
  ),
  plot = p_network,
  width = 10,
  height = 9.2,
  limitsize = FALSE
)

ggplot2::ggsave(
  filename = file.path(
    fig6_dir,
    "Figure6D_core_lncRNA_miRNA_mRNA_network_circular_compact.png"
  ),
  plot = p_network,
  width = 10,
  height = 9.2,
  dpi = 600,
  limitsize = FALSE
)

cat(
  "\nFigure 6D 分析及作图全部完成。\n"
)