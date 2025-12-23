library(pbmcapply)

# -------- 单次估计（基于你的版本微改） --------
estimate <- function(i, results, sect, A_diag, kappa){
  res <- results[[i]]
  # 万一有失败样本，容错返回 NA
  if (is.null(res)) return(c(NA_real_, NA_real_, NA_real_))
  
  beta_hat   <- res$beta_hat
  beta_lasso <- res$beta_lasso
  param      <- res$param
  Sigma_diag <- res$Sigma_diag
  
  p <- length(A_diag)  # <- 你原来用到 p，但没定义
  
  # 对角 A 的二次型：x' A x = sum(A_diag * x^2)
  qform <- function(x, w) sum(w * (x^2))
  
  if (sect == "norm" || sect == "S" || sect == "Sc") {
    est_naive <- qform(beta_hat, A_diag)
    est_corrected <- max(
      ( est_naive - (1 - kappa) * (param[3]^2) * sum(A_diag * Sigma_diag) / p ) / (param[1]^2),
      0
    )
    est_lasso <- qform(beta_lasso, A_diag)
    
  } else if (sect == "k") {
    est_naive     <- beta_hat[1]^2
    est_corrected <- (beta_hat[1] / param[1])^2
    est_lasso     <- beta_lasso[1]^2
    
  } else {
    return(c(NA_real_, NA_real_, NA_real_))
  }
  
  return(c(est_naive, est_corrected, est_lasso))
}

# -------- 主循环（尽量不改你的结构） --------
setting <- c("dense", "sparse")
sections <- c("norm", "k", "S", "Sc")

nn      <- c(1000L, 2000L, 4000L)
kappas  <- c(0, 0.1, 0.2)
ga      <- sqrt(5)

simu_total <- 2000L

outdir <- "simu_data"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
if (!dir.exists("simu_est")) dir.create("simu_est", recursive = TRUE)

cores <- max(1L, parallel::detectCores() - 2L)  # 别用满核

for (sect in sections){
  for (sett in setting){
    for (kappa in kappas) {
      for (j in seq_along(nn)) {
        n <- nn[j]
        if (kappa == 0) {
          p <- 20L
        } else {
          p <- as.integer(n * kappa)
        }
        if (p < 1L) next
        
        message(sprintf("[START] sect=%s, setting=%s, kappa=%.2f, n=%d, p=%d",
                        sect, sett, kappa, n, p))
        
        # A 的对角（按你原逻辑）
        A_diag <- list()
        A_diag[["norm"]] <- rep(1, p)
        A_diag[["k"]]    <- c(1, rep(0, max(0, p-1)))
        A_diag[["S"]]    <- c(1, 1, rep(0, max(0, p-2)))
        A_diag[["Sc"]]   <- c(0, 0, rep(1, max(0, p-2)))
        
        # 载入 bundle
        fin <- file.path(outdir, sprintf("%s_k=%.2f_n=%d_simu=%d.RData",
                                         sett, kappa, n, simu_total))
        if (!file.exists(fin)) {
          message("[SKIP] Not found: ", fin)
          next
        }
        load(fin)  # 得到 bundle
        beta <- bundle$beta
        
        # 并行计算（传入 kappa）
        simu_ma_l <- pbmclapply(
          X = 1:simu_total, FUN = estimate,
          results = bundle$results,
          sect = sect, A_diag = A_diag[[sect]], kappa = kappa,
          mc.cores = cores, mc.preschedule = FALSE, mc.cleanup = TRUE
        )
        
        simu_ma <- do.call(cbind, simu_ma_l)
        
        # 删除 NA 样本
        keep_cols <- colSums(is.na(simu_ma)) == 0
        simu_ma <- simu_ma[, keep_cols, drop = FALSE]
        
        # 真值
        truth <- sum(A_diag[[sect]] * (beta^2))
        
        # 各估计量向量
        naive     <- as.numeric(simu_ma[1, ])
        corrected <- as.numeric(simu_ma[2, ])
        lasso     <- as.numeric(simu_ma[3, ])
        
        # 结果列表
        S <- list()
        S$naive     <- naive
        S$corrected <- corrected
        S$lasso     <- lasso
        S$truth     <- truth
        
        # 统计量
        stats <- list(
          naive     = c(mean = mean(naive),     var = var(naive),     mse = mean((naive - truth)^2)),
          corrected = c(mean = mean(corrected), var = var(corrected), mse = mean((corrected - truth)^2)),
          lasso     = c(mean = mean(lasso),     var = var(lasso),     mse = mean((lasso - truth)^2))
        )
        S$stats <- stats
        
        # 打印结果
        message(sprintf("  [RESULT] naive    : mean=%.4f var=%.4f mse=%.4f", 
                        stats$naive["mean"], stats$naive["var"], stats$naive["mse"]))
        message(sprintf("           corrected: mean=%.4f var=%.4f mse=%.4f", 
                        stats$corrected["mean"], stats$corrected["var"], stats$corrected["mse"]))
        message(sprintf("           lasso    : mean=%.4f var=%.4f mse=%.4f", 
                        stats$lasso["mean"], stats$lasso["var"], stats$lasso["mse"]))
        message(sprintf("           truth    : %.4f", truth))
        
        fout <- file.path("simu_est",
                          sprintf("%s_%s_k=%.2f_n=%d_simu=%d.RData",
                                  sect, sett, kappa, n, simu_total))
        save(S, file = fout)
        message("[SAVED] ", fout)
        
        rm(simu_ma_l, simu_ma, S, bundle); gc()
      }
    }
  }
}
