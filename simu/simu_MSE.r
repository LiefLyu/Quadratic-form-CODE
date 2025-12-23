## ===== 计算各图的 MSE =====

## 参数（按你的设定改）
simu_total <- 2000L
nn         <- c(1000L, 2000L, 4000L)
kappas     <- c(0, 0.1, 0.2)
sections   <- c("norm", "k", "S", "Sc")
settings   <- c("dense", "sparse")

## 方法键（S[[1]], S[[2]], S[[3]]）
method_keys <- c("MLE_naive", "MLE_corrected", "MLE_Lasso")

## 生成文件路径
make_fpath <- function(sect, sett, kappa, n, simu_total) {
  file.path("simu_est",
            sprintf("%s_%s_k=%.2f_n=%d_simu=%d.RData",
                    sect, sett, kappa, n, simu_total))
}

## 从一个文件（一个 n）读取 S 并计算三种方法的 MSE
## 返回 data.frame: (method_key, mse, n_rep_used, truth)
mse_from_file <- function(fpath) {
  if (!file.exists(fpath)) {
    warning("File not found: ", fpath)
    return(NULL)
  }
  load(fpath)  # expects object S: list of length >= 4
  if (!exists("S") || !is.list(S) || length(S) < 4) {
    warning("Malformed object S in: ", fpath)
    return(NULL)
  }
  truth <- as.numeric(S[[4]])
  out <- lapply(1:3, function(m) {
    vec <- as.numeric(S[[m]])
    ok  <- is.finite(vec)
    n_rep <- sum(ok)
    if (n_rep == 0L) {
      data.frame(method_key = method_keys[m], mse = NA_real_, n_rep_used = 0L, truth = truth)
    } else {
      mse_val <- mean((vec[ok] - truth)^2, na.rm = TRUE)
      data.frame(method_key = method_keys[m], mse = mse_val, n_rep_used = n_rep, truth = truth)
    }
  })
  do.call(rbind, out)
}

## 主循环：遍历所有 sect/sett/kappa/n 组合
rows <- list()
idx  <- 1L
for (sect in sections) {
  for (sett in settings) {
    ## 如果你确定 S/Sc 只做某个 setting，可以在这里跳过：
    ## 例如只做 sparse: if (sect %in% c("S","Sc") && sett == "dense") next
    
    for (kappa in kappas) {
      for (n in nn) {
        fpath <- make_fpath(sect, sett, kappa, n, simu_total)
        res   <- mse_from_file(fpath)
        if (is.null(res)) next
        res$section <- sect
        res$design  <- sett
        res$kappa   <- kappa
        res$n       <- n
        rows[[idx]] <- res
        idx <- idx + 1L
      }
    }
  }
}

mse_all <- do.call(rbind, rows)

## 三位小数输出
mse_all$mse_3 <- round(mse_all$mse, 3)

## 排序显示（可选）
mse_all <- mse_all[order(mse_all$section, mse_all$design, mse_all$kappa, mse_all$n,
                         match(mse_all$method_key, method_keys)), ]

## 查看汇总
print(head(mse_all, 12))
cat("Total rows:", nrow(mse_all), "\n")

## 保存 CSV（可选）
dir.create("plot", showWarnings = FALSE, recursive = TRUE)
write.csv(mse_all, file = file.path("plot", "mse_all.csv"), row.names = FALSE)

## —— 可选：生成某些小节/设计的宽表，便于直接放到 LaTeX —— ##
## 例：norm + dense，按 n 行、方法列，分别给每个 kappa 一张表
subset_and_wide <- function(df, sect, sett, kappa) {
  sub <- subset(df, section == sect & design == sett & kappa == kappa)
  if (nrow(sub) == 0L) return(NULL)
  ## 宽表：n 行，三方法列（保留三位小数）
  n_vals <- sort(unique(sub$n))
  wide <- data.frame(n = n_vals)
  for (mk in method_keys) {
    tmp <- sub[sub$method_key == mk, c("n","mse_3")]
    tmp <- tmp[match(n_vals, tmp$n), "mse_3"]
    wide[[mk]] <- tmp
  }
  wide
}

## 示例：提取 norm + dense, kappa=0.10 的宽表
wide_norm_dense_010 <- subset_and_wide(mse_all, "norm", "dense", 0.10)
print(wide_norm_dense_010)

subset_and_wide(mse_all, "norm", "dense", 0.010)
subset_and_wide(mse_all, "norm", "dense", 0.10)

subset(mse_all, section == "norm" & design == "dense" & kappa == 0.010)
subset(mse_all, section == "k")
