############################################################
## Figure 3A-B only
## Nomogram + Calibration curve
## Disease: Colorectal cancer
## Genes: AKR1C2, MAPK3, MAP3K11
############################################################

rm(list = ls())
options(stringsAsFactors = FALSE)

############################################################
## 1. 加载 R 包
############################################################

library(rms)
library(dplyr)

############################################################
## 2. 设置路径
############################################################

workdir <- "C:/Users/92066/Desktop/单细胞孟德尔/Figure3"
setwd(workdir)

input_file <- "Figure3_nomogram_input.csv"

############################################################
## 3. 读取数据
############################################################

dat <- read.csv(
  input_file,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("数据维度：\n")
print(dim(dat))

cat("数据列名：\n")
print(colnames(dat))

############################################################
## 4. 整理数据
############################################################

genes_use <- c("AKR1C2", "MAPK3", "MAP3K11")

need_cols <- c("sample_id", "diagnosis", "tumor_status", genes_use)

missing_cols <- setdiff(need_cols, colnames(dat))

if (length(missing_cols) > 0) {
  stop(
    paste0(
      "以下列在输入文件中不存在：",
      paste(missing_cols, collapse = ", ")
    )
  )
}

dat <- dat %>%
  dplyr::select(
    sample_id,
    diagnosis,
    tumor_status,
    dplyr::all_of(genes_use)
  )

dat$tumor_status <- as.numeric(dat$tumor_status)

dat[, genes_use] <- lapply(
  dat[, genes_use, drop = FALSE],
  function(x) as.numeric(as.character(x))
)

dat <- dat %>%
  dplyr::filter(
    !is.na(tumor_status),
    !is.na(AKR1C2),
    !is.na(MAPK3),
    !is.na(MAP3K11),
  )

cat("清洗后数据维度：\n")
print(dim(dat))

cat("分组统计：\n")
print(table(dat$diagnosis, dat$tumor_status))

if (length(unique(dat$tumor_status)) != 2) {
  stop("tumor_status 必须包含 0 和 1 两组，请检查输入文件。")
}

############################################################
## 5. 构建 rms logistic 模型
############################################################

dd <- datadist(dat)
options(datadist = "dd")

fit <- lrm(
  tumor_status ~ AKR1C2 + MAPK3 + MAP3K11,
  data = dat,
  x = TRUE,
  y = TRUE
)

cat("模型结果：\n")
print(fit)

## 保存预测概率，B 图 rug 会用到
dat$predicted_prob <- predict(fit, type = "fitted")

############################################################
## 6. Figure 3A: Nomogram
## 要求：
## 1. 删除 Linear Predictor
## 2. 最后一行改为 Predicted probability of Tumor
## 3. 使用 rms 稳定版，保证 PDF/PNG 可打开
############################################################

nom <- nomogram(
  fit,
  lp = FALSE,
  fun = plogis,
  fun.at = c(0.1, 0.3, 0.5, 0.7, 0.9),
  funlabel = "Predicted probability of Tumor"
)

plot_nomogram_A <- function() {
  
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  
  plot(
    nom,
    xfrac = 0.35,
    cex.axis = 1.05,
    cex.var = 1.15,
    lmgp = 0.25
  )
}

pdf(
  file = "Figure3A_nomogram_CRC.pdf",
  width = 10,
  height = 7
)
plot_nomogram_A()
dev.off()

png(
  filename = "Figure3A_nomogram_CRC.png",
  width = 3000,
  height = 2100,
  res = 300
)
plot_nomogram_A()
dev.off()

cat("Figure3A nomogram finished.\n")

############################################################
## 7. Figure 3B: Calibration curve
## 要求：
## 1. 彻底删除 C.L. 曲线
## 2. 彻底删除 C.L. 图注
## 3. 只保留 Apparent、Bias-corrected、Ideal
## 4. 保留 B、Mean absolute error、n
############################################################

set.seed(123)

boot_B <- 40

cal <- calibrate(
  fit,
  method = "boot",
  B = boot_B
)

############################################################
## 7.1 稳定提取 calibrate 结果中的指定列
############################################################

get_cal_col <- function(cal_obj, col_name) {
  
  v <- tryCatch(
    {
      cal_obj[, col_name]
    },
    error = function(e) {
      NULL
    }
  )
  
  if (!is.null(v)) {
    return(as.numeric(v))
  }
  
  cal_unclass <- tryCatch(
    {
      unclass(cal_obj)
    },
    error = function(e) {
      NULL
    }
  )
  
  if (is.matrix(cal_unclass) && col_name %in% colnames(cal_unclass)) {
    return(as.numeric(cal_unclass[, col_name]))
  }
  
  if (is.data.frame(cal_unclass) && col_name %in% colnames(cal_unclass)) {
    return(as.numeric(cal_unclass[[col_name]]))
  }
  
  stop(
    paste0(
      "无法从 calibrate 对象中提取列：",
      col_name,
      "。请检查 calibrate 输出结构。"
    )
  )
}

cal_pred <- get_cal_col(cal, "predy")
cal_app  <- get_cal_col(cal, "calibrated.orig")
cal_bias <- get_cal_col(cal, "calibrated.corrected")

idx <- is.finite(cal_pred) & is.finite(cal_app) & is.finite(cal_bias)

cal_pred <- cal_pred[idx]
cal_app  <- cal_app[idx]
cal_bias <- cal_bias[idx]

if (length(cal_pred) == 0) {
  stop("calibration 曲线没有有效点，请检查模型或 calibrate 结果。")
}

mae_value <- mean(
  abs(cal_pred - cal_bias),
  na.rm = TRUE
)

n_value <- nrow(dat)

cal_subtitle <- paste0(
  "B= ", boot_B, " repetitions, boot     ",
  "Mean absolute error=", sprintf("%.3f", mae_value),
  " n=", n_value
)

cat("Calibration subtitle:\n")
cat(cal_subtitle, "\n")

############################################################
## 7.2 手动画校准曲线
############################################################

plot_calibration_B <- function() {
  
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  
  par(
    mar = c(7.2, 5.2, 2.0, 2.0),
    xpd = FALSE
  )
  
  plot(
    NA,
    NA,
    xlim = c(0, 1),
    ylim = c(0, 1),
    xlab = "Nomogram-predicted probability",
    ylab = "Observed probability",
    cex.lab = 1.35,
    cex.axis = 1.15,
    bty = "o"
  )
  
  ## Ideal line
  lines(
    c(0, 1),
    c(0, 1),
    lty = 2,
    lwd = 1.8,
    col = "black"
  )
  
  ## Apparent curve
  lines(
    cal_pred,
    cal_app,
    lty = 3,
    lwd = 2.0,
    col = "black"
  )
  
  ## Bias-corrected curve
  lines(
    cal_pred,
    cal_bias,
    lty = 1,
    lwd = 2.2,
    col = "black"
  )
  
  ## 顶部 rug，显示预测概率分布
  rug(
    dat$predicted_prob,
    side = 3,
    ticksize = 0.025,
    col = "gray40"
  )
  
  ## 图注：只保留三项，不含 C.L.
  legend(
    "bottomright",
    legend = c("Apparent", "Bias-corrected", "Ideal"),
    lty = c(3, 1, 2),
    lwd = c(2.0, 2.2, 1.8),
    col = c("black", "black", "black"),
    bty = "n",
    cex = 1.15
  )
  
  ## 底部文字
  par(xpd = NA)
  mtext(
    cal_subtitle,
    side = 1,
    line = 4.8,
    cex = 1.00,
    adj = 0.5
  )
}

pdf(
  file = "Figure3B_calibration_CRC.pdf",
  width = 6.5,
  height = 6.5
)
plot_calibration_B()
dev.off()

png(
  filename = "Figure3B_calibration_CRC.png",
  width = 2100,
  height = 2100,
  res = 300
)
plot_calibration_B()
dev.off()

cat("Figure3B calibration curve finished.\n")

############################################################
## 8. 完成
############################################################

cat("\nFigure 3A-B 作图完成！输出文件包括：\n")
cat("Figure3A_nomogram_CRC.pdf / png\n")
cat("Figure3B_calibration_CRC.pdf / png\n")
cat("\nB图底部信息：\n")
cat(cal_subtitle, "\n")

############################################################
## Figure 3
## Construction of a nomogram in GSE44861
## Disease: Colorectal cancer
## Genes: AKR1C2, MAPK3, MAP3K11
############################################################

rm(list = ls())
options(stringsAsFactors = FALSE)

############################################################
## 1. 加载 R 包
############################################################

library(rms)
library(pROC)
library(ggplot2)
library(dplyr)
library(rmda)

############################################################
## 2. 设置路径
############################################################

workdir <- "C:/Users/92066/Desktop/单细胞孟德尔/Figure3"
setwd(workdir)

input_file <- "Figure3_nomogram_input.csv"

############################################################
## 3. 读取数据
############################################################

dat <- read.csv(
  input_file,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("数据维度：\n")
print(dim(dat))

cat("数据列名：\n")
print(colnames(dat))

cat("数据预览：\n")
print(head(dat))

############################################################
## 4. 检查基因和结局变量
############################################################

genes_use <- c("AKR1C2", "MAPK3", "MAP3K11")

need_cols <- c("sample_id", "diagnosis", "tumor_status", genes_use)

missing_cols <- setdiff(need_cols, colnames(dat))

if (length(missing_cols) > 0) {
  stop(
    paste0(
      "以下列在输入文件中不存在：",
      paste(missing_cols, collapse = ", ")
    )
  )
}

dat <- dat %>%
  select(sample_id, diagnosis, tumor_status, all_of(genes_use))

dat$tumor_status <- as.numeric(dat$tumor_status)

dat <- dat %>%
  filter(!is.na(tumor_status))

dat[, genes_use] <- lapply(dat[, genes_use, drop = FALSE], as.numeric)

cat("分组统计：\n")
print(table(dat$diagnosis, dat$tumor_status))

############################################################
## 5. 设置 rms 环境
############################################################

dd <- datadist(dat)
options(datadist = "dd")

############################################################
## 6. 构建 Logistic 回归模型
## 原逻辑：疾病状态 ~ 4个biomarker
## 你的研究：tumor_status ~ AKR1C2 + MAPK3 + MAP3K11
############################################################

fit <- lrm(
  tumor_status ~ AKR1C2 + MAPK3 + MAP3K11,
  data = dat,
  x = TRUE,
  y = TRUE
)

cat("模型结果：\n")
print(fit)

## 提前保存模型预测概率：Calibration 的 rug、ROC 和 DCA 都会使用
dat$predicted_prob <- as.numeric(predict(fit, type = "fitted"))

############################################################
## 7. Figure 3A: Nomogram
############################################################

nom <- nomogram(
  fit,
  fun = plogis,
  fun.at = c(0.1, 0.3, 0.5, 0.7, 0.9),
  funlabel = "Risk of Tumor"
)

pdf("Figure3A_nomogram_CRC.pdf", width = 10, height = 7)
plot(
  nom,
  xfrac = 0.35,
  cex.axis = 1.0,
  cex.var = 1.1,
  lmgp = 0.25
)
dev.off()

png("Figure3A_nomogram_CRC.png", width = 3000, height = 2100, res = 300)
plot(
  nom,
  xfrac = 0.35,
  cex.axis = 1.0,
  cex.var = 1.1,
  lmgp = 0.25
)
dev.off()

############################################################
## 8. Figure 3B: Calibration curve
## 目的：保留 rms 原始 calibration 作图逻辑，同时显示：
## B= 40 repetitions, boot    Mean absolute error=xxx n=xxx
############################################################

set.seed(123)

## 原文 Figure 3B 显示 B=40 repetitions, boot
## 因此这里将 B 改为 40；如需更稳定结果，可改为 1000，但图中会显示 B=1000
boot_B <- 40

cal <- calibrate(
  fit,
  method = "boot",
  B = boot_B
)

## 稳定提取 calibrate 结果（兼容不同 rms 版本）
## 对 lrm，常见列名为 predy / calibrated.orig / calibrated.corrected；
## 某些版本或对象可能使用 pred / calibrated，因此这里同时兼容。
get_cal_components <- function(cal_obj) {
  cal_mat <- tryCatch(unclass(cal_obj), error = function(e) NULL)
  
  if (is.data.frame(cal_mat)) {
    cal_mat <- as.matrix(cal_mat)
  }
  if (!is.matrix(cal_mat)) {
    cal_mat <- tryCatch(as.matrix(cal_obj), error = function(e) NULL)
  }
  if (is.null(cal_mat) || is.null(colnames(cal_mat))) {
    stop("无法读取 calibrate 对象的列名。请运行 colnames(cal) 和 str(cal) 检查对象结构。")
  }
  
  cn <- colnames(cal_mat)
  
  pred_col <- if ("predy" %in% cn) {
    "predy"
  } else if ("pred" %in% cn) {
    "pred"
  } else if ("Predicted" %in% cn) {
    "Predicted"
  } else {
    stop(paste0(
      "calibrate 对象中未找到预测概率列（predy/pred/Predicted）。现有列：",
      paste(cn, collapse = ", ")
    ))
  }
  
  app_col <- if ("calibrated.orig" %in% cn) {
    "calibrated.orig"
  } else if ("calibrated" %in% cn) {
    "calibrated"
  } else {
    stop(paste0(
      "calibrate 对象中未找到 apparent calibration 列。现有列：",
      paste(cn, collapse = ", ")
    ))
  }
  
  corr_col <- if ("calibrated.corrected" %in% cn) {
    "calibrated.corrected"
  } else if ("bias.corrected" %in% cn) {
    "bias.corrected"
  } else {
    stop(paste0(
      "calibrate 对象中未找到 bias-corrected calibration 列。现有列：",
      paste(cn, collapse = ", ")
    ))
  }
  
  list(
    pred = as.numeric(cal_mat[, pred_col]),
    apparent = as.numeric(cal_mat[, app_col]),
    corrected = as.numeric(cal_mat[, corr_col])
  )
}

cal_comp <- get_cal_components(cal)

## 按 rms::plot.calibrate.default 的思路计算 MAE：
## 在每个原始样本的预测概率处，对 bias-corrected calibration curve 做插值，
## 再计算 |predicted - calibrated| 的平均值。
predicted_for_mae <- attr(cal, "predicted")
if (is.null(predicted_for_mae)) {
  predicted_for_mae <- dat$predicted_prob
}
predicted_for_mae <- as.numeric(predicted_for_mae)

ok_cal <- is.finite(cal_comp$pred) & is.finite(cal_comp$corrected)
if (sum(ok_cal) < 2) {
  stop("有效 calibration 点少于 2 个，无法计算 Mean absolute error。")
}

interp_corrected <- approx(
  x = cal_comp$pred[ok_cal],
  y = cal_comp$corrected[ok_cal],
  xout = predicted_for_mae,
  ties = mean
)$y

mae_value <- mean(abs(predicted_for_mae - interp_corrected), na.rm = TRUE)
n_value <- nrow(dat)

cal_subtitle <- paste0(
  "B= ", boot_B, " repetitions, boot     ",
  "Mean absolute error=", sprintf("%.3f", mae_value),
  " n=", n_value
)

cat("Calibration columns:\n")
print(colnames(unclass(cal)))
cat("Calibration subtitle:\n")
cat(cal_subtitle, "\n")

## 统一作图函数：仅保留 Apparent、Bias-corrected、Ideal，不显示 C.L.
plot_calibration_crc <- function(cal_obj, cal_subtitle) {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  
  cc <- get_cal_components(cal_obj)
  keep <- is.finite(cc$pred) & is.finite(cc$apparent) & is.finite(cc$corrected)
  
  cal_pred <- cc$pred[keep]
  cal_app  <- cc$apparent[keep]
  cal_bias <- cc$corrected[keep]
  
  par(mar = c(6.4, 5.2, 2.0, 2.0), xpd = NA)
  
  plot(
    NA, NA,
    xlim = c(0, 1),
    ylim = c(0, 1),
    xlab = "Nomogram-predicted probability",
    ylab = "Observed probability",
    cex.lab = 1.25,
    cex.axis = 1.10,
    bty = "o"
  )
  
  ## Ideal line
  lines(c(0, 1), c(0, 1), lty = 2, lwd = 1.8, col = "black")
  
  ## Apparent curve
  lines(cal_pred, cal_app, lty = 3, lwd = 2.0, col = "black")
  
  ## Bias-corrected curve
  lines(cal_pred, cal_bias, lty = 1, lwd = 2.2, col = "black")
  
  ## 顶部 rug：模型预测概率分布
  rug(dat$predicted_prob, side = 3, ticksize = 0.025, col = "gray40")
  
  ## 图注：只保留三项，不含 C.L.
  legend(
    "bottomright",
    legend = c("Apparent", "Bias-corrected", "Ideal"),
    lty = c(3, 1, 2),
    lwd = c(2.0, 2.2, 1.8),
    col = c("black", "black", "black"),
    bty = "n",
    cex = 1.10
  )
  
  mtext(
    cal_subtitle,
    side = 1,
    line = 4.6,
    cex = 0.95,
    adj = 0.5
  )
}

pdf("Figure3B_calibration_CRC.pdf", width = 6.5, height = 6.5)
plot_calibration_crc(cal, cal_subtitle)
dev.off()

png("Figure3B_calibration_CRC.png", width = 2100, height = 2100, res = 300)
plot_calibration_crc(cal, cal_subtitle)
dev.off()

############################################################
## 9. Figure 3C: ROC curve
############################################################

roc_obj <- roc(
  response = dat$tumor_status,
  predictor = dat$predicted_prob,
  levels = c(0, 1),
  direction = "<"
)

auc_value <- auc(roc_obj)
ci_auc <- ci.auc(roc_obj)

cat("AUC:\n")
print(auc_value)

cat("95% CI:\n")
print(ci_auc)

pdf("Figure3C_ROC_CRC.pdf", width = 6, height = 6)
plot(
  roc_obj,
  col = "#D73027",
  lwd = 2.5,
  legacy.axes = TRUE,
  xlab = "1 - Specificity",
  ylab = "Sensitivity",
  main = "",
  identity = TRUE,
  identity.col = "gray",
  identity.lty = 2,
  identity.lwd = 1.5
)
text(
  0.60,
  0.20,
  labels = paste0(
    "AUC: ",
    round(as.numeric(auc_value), 3),
    " (",
    round(ci_auc[1], 3),
    "-",
    round(ci_auc[3], 3),
    ")"
  ),
  cex = 1.1,
  col = "#D73027"
)
dev.off()

png("Figure3C_ROC_CRC.png", width = 1800, height = 1800, res = 300)
plot(
  roc_obj,
  col = "#D73027",
  lwd = 2.5,
  legacy.axes = TRUE,
  xlab = "1 - Specificity",
  ylab = "Sensitivity",
  main = "",
  identity = TRUE,
  identity.col = "gray",
  identity.lty = 2,
  identity.lwd = 1.5
)
text(
  0.60,
  0.20,
  labels = paste0(
    "AUC: ",
    round(as.numeric(auc_value), 3),
    " (",
    round(ci_auc[1], 3),
    "-",
    round(ci_auc[3], 3),
    ")"
  ),
  cex = 1.1,
  col = "#D73027"
)
dev.off()

############################################################
## 10. Figure 3D: Decision curve analysis
############################################################

## 完整 Nomogram 的 predicted_prob 已经是模型预测风险，因此使用 fitted.risk = TRUE
dca_nomogram <- decision_curve(
  tumor_status ~ predicted_prob,
  data = dat,
  thresholds = seq(0, 1, by = 0.01),
  confidence.intervals = FALSE,
  study.design = "cohort",
  fitted.risk = TRUE
)

## 单基因模型
dca_AKR1C2 <- decision_curve(
  tumor_status ~ AKR1C2,
  data = dat,
  family = binomial(link = "logit"),
  thresholds = seq(0, 1, by = 0.01),
  confidence.intervals = FALSE,
  study.design = "cohort"
)

dca_MAPK3 <- decision_curve(
  tumor_status ~ MAPK3,
  data = dat,
  family = binomial(link = "logit"),
  thresholds = seq(0, 1, by = 0.01),
  confidence.intervals = FALSE,
  study.design = "cohort"
)

dca_MAP3K11 <- decision_curve(
  tumor_status ~ MAP3K11,
  data = dat,
  family = binomial(link = "logit"),
  thresholds = seq(0, 1, by = 0.01),
  confidence.intervals = FALSE,
  study.design = "cohort"
)

dca_list <- list(
  dca_nomogram,
  dca_AKR1C2,
  dca_MAPK3,
  dca_MAP3K11
)

pdf("Figure3D_DCA_CRC.pdf", width = 7, height = 6)
plot_decision_curve(
  dca_list,
  curve.names = c("Nomogram", "AKR1C2", "MAPK3", "MAP3K11"),
  xlab = "High risk threshold",
  ylab = "Net benefit",
  standardize = FALSE,
  lwd = 2.2,
  legend.position = "topright"
)
dev.off()

png("Figure3D_DCA_CRC.png", width = 2100, height = 1800, res = 300)
plot_decision_curve(
  dca_list,
  curve.names = c("Nomogram", "AKR1C2", "MAPK3", "MAP3K11"),
  xlab = "High risk threshold",
  ylab = "Net benefit",
  standardize = FALSE,
  lwd = 2.2,
  legend.position = "topright"
)
dev.off()

############################################################
## 11. 保存模型预测结果
############################################################

out_dat <- dat %>%
  select(
    sample_id,
    diagnosis,
    tumor_status,
    all_of(genes_use),
    predicted_prob
  )

write.csv(
  out_dat,
  "Figure3_nomogram_prediction_CRC.csv",
  row.names = FALSE,
  quote = FALSE
)

############################################################
## 12. 完成
############################################################

cat("\nFigure 3 作图完成！输出文件包括：\n")
cat("Figure3A_nomogram_CRC.pdf / png\n")
cat("Figure3B_calibration_CRC.pdf / png\n")
cat("Figure3C_ROC_CRC.pdf / png\n")
cat("Figure3D_DCA_CRC.pdf / png\n")
cat("Figure3_nomogram_prediction_CRC.csv\n")



############################################################
## 输出 Figure 3 logistic model 的完整系数表
## 包括 Coefficient、OR、95%CI、P value
############################################################

glm_fit <- glm(
  tumor_status ~ AKR1C2 + MAPK3 + MAP3K11,
  data = dat,
  family = binomial(link = "logit")
)

coef_table <- summary(glm_fit)$coefficients

coef_out <- as.data.frame(coef_table)
coef_out$Variable <- rownames(coef_out)

colnames(coef_out) <- c(
  "Coefficient",
  "Std_Error",
  "Z_value",
  "P_value",
  "Variable"
)

coef_out <- coef_out %>%
  select(Variable, Coefficient, Std_Error, Z_value, P_value) %>%
  mutate(
    OR = exp(Coefficient),
    OR_95CI_low = exp(Coefficient - 1.96 * Std_Error),
    OR_95CI_high = exp(Coefficient + 1.96 * Std_Error),
    Direction = case_when(
      Variable == "(Intercept)" ~ "Intercept (not interpreted biologically)",
      Coefficient > 0 ~ "Higher expression is associated with higher Tumor risk",
      Coefficient < 0 ~ "Higher expression is associated with lower Tumor risk",
      TRUE ~ "No clear direction"
    )
  )

## 不再对 Coefficient 做 round()：
## 保留 glm 模型中的原始双精度系数，供后续 fixed-model external validation 使用。
## 其他数值列也保持原始计算精度，避免 CSV 中人为截断。

## 严格检查：导出的 Coefficient 必须与 glm_fit 的系数完全一致
stopifnot(
  isTRUE(
    all.equal(
      coef_out$Coefficient,
      as.numeric(coef(glm_fit)),
      tolerance = 0,
      check.attributes = FALSE
    )
  )
)

cat("\nFull-precision glm coefficients:\n")
print(coef_out, digits = 17)

## write.csv 写出前临时提高数字显示精度；
## quote=TRUE 可保证 Direction 中即使包含标点也不会破坏 CSV 列结构。
old_digits <- getOption("digits")
options(digits = 17)

write.csv(
  coef_out,
  file = "Figure3_model_coefficients_OR_95CI_Pvalue_CRC.csv",
  row.names = FALSE,
  quote = TRUE
)

options(digits = old_digits)

## 同时保存完整 glm 模型对象，后续如需最严格地复用固定模型，
## 可直接 readRDS() 后使用 coef()，不会发生任何 CSV 显示精度损失。
saveRDS(
  glm_fit,
  file = "Figure3_glm_model_CRC.rds"
)