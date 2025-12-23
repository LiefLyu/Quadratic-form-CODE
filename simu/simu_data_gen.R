# =========================
# Sim code
# =========================
library(pbmcapply)
library(glmnet)
# 可选：限制 BLAS/OMP 线程，避免多进程×多线程争用
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(1L)
  RhpcBLASctl::omp_set_num_threads(1L)
}

# ---- 你的 SLOE 函数 ----
source("SLOE.R")

estimate <- function(i, beta, R, n, p, kappa, sett, timeout_sec = 60){
  # 小心：进程内大量打印会拖慢并行，建议只偶尔输出
  # if (i %% 200 == 0) message("i = ", i)
  
  runner <- function(){
    set.seed(i)
    
    # ---- 你原有的主体开始 ----
    X <- matrix(rnorm(n * p), n, p) %*% R
    eta <- as.vector(X %*% beta)
    Y <- rbinom(n, 1, 1 / (1 + exp(-eta)))
    
    fit <- try(glm(Y ~ X + 0, family = binomial(), x = TRUE, y = TRUE), silent = TRUE)
    if (inherits(fit, "try-error")) return(list(ok=FALSE, i=i, reason="glm_error"))
    beta_hat <- coef(fit)
    if (anyNA(beta_hat)) return(list(ok=FALSE, i=i, reason="glm_na"))
    
    etahat <- fun_SLOE_fast(X, Y, beta_hat)
    param  <- find_param_sloe(kappa = p/n, eta = sqrt(etahat), intercept = FALSE)
    
    cv.fit <- cv.glmnet(X, Y, alpha = 1, family = "binomial", intercept = FALSE)
    beta_lasso <- as.numeric(coef(cv.fit, s = cv.fit$lambda.min))[-1]
    
    # 样本协方差（除以 n-1）及其逆的对角元
    Sigma_diag <- diag(chol2inv(chol(cov(X))))
    
    list(
      ok          = TRUE,
      i           = i,
      beta_hat    = as.numeric(beta_hat),
      beta_lasso  = beta_lasso,
      param       = param,
      Sigma_diag  = Sigma_diag
    )
    # ---- 你原有的主体结束 ----
  }
  
  # 超时保护：优先使用 R.utils；若未安装则不设超时
  if (requireNamespace("R.utils", quietly = TRUE)) {
    out <- try(
      R.utils::withTimeout(expr = runner(),
                           timeout = timeout_sec, cpu = timeout_sec,
                           onTimeout = "error"),
      silent = TRUE
    )
    if (inherits(out, "try-error")) {
      return(list(ok=FALSE, i=i, reason="timeout"))
    } else {
      return(out)
    }
  } else {
    # 没装 R.utils 就直接跑（无超时）
    return(runner())
  }
}

# ---- 实验设置 ----
setting <- c("dense", "sparse")
nn      <- c(1000L, 2000L, 4000L)
# nn <- c(4000L)
# kappas  <- c(0, 0.1, 0.2)
kappas  <- 0
ga      <- sqrt(5)

# 并行核数（建议留 1 个核）
cores <- max(1L, parallel::detectCores() - 1L)

# 模拟次数
simu_total <- 2000L

# 输出目录
outdir <- "simu_data"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# ---- 主循环（每个组合只保存一次）----
for (sett in setting){
  for (kappa in kappas) {
    for (j in 1:length(nn)) {
      n <- nn[j]
      if (kappa == 0) {
        p <- 20L
      } else {
        p <- as.integer(n * kappa)
      }
      stopifnot(p >= 2L)
      
      # 构造 beta
      set.seed(123)  # 使不同组合下 beta 可复现
      if (sett == "dense"){
        b1 <- c(1,1)
        w  <- runif(p - 2L)
        # 向量范数用 sqrt(sum(w^2))，替代 norm(w,"2")
        w  <- w / sqrt(sum(w^2)) * sqrt(ga^2 - 2)
        beta <- c(b1, w)
      } else { # sparse
        beta <- c(1, 2, rep(0, p - 2L))
      }
      
      # 相关性矩阵的 Cholesky（Toeplitz 结构）
      R <- chol(toeplitz(0.2^(0:(p - 1L))))
      
      # 为每次并行任务生成独立、可复现的种子
      seeds_all <- 1:simu_total
      
      message(sprintf("[RUN] setting=%s | kappa=%.2f | n=%d | p=%d | simu=%d",
                      sett, kappa, n, p, simu_total))
      
      # 并行模拟，返回一个长度为 simu_total 的 list
      res_list <- pbmclapply(
        X = seeds_all, FUN = estimate,
        mc.cores = cores,
        mc.preschedule = FALSE,  # 动态调度，避免尾部几个慢任务拖住
        mc.cleanup = TRUE,
        beta = beta, R = R, n = n, p = p, kappa = kappa, sett = sett,
        timeout_sec = 180          # <- 每样本超时秒数
      )
      
      # 统计失败条目（一般很少/无）
      ok_flags  <- vapply(res_list, function(z) is.list(z) && isTRUE(z$ok), logical(1))
      timeouts  <- sum(vapply(res_list, function(z) is.list(z) && identical(z$reason, "timeout"), logical(1)))
      message(sprintf("OK: %d / %d, TIMEOUT: %d", sum(ok_flags), length(res_list), timeouts))
      
      # ---- 打包 & 保存 ----
      bundle <- list(
        meta = list(
          setting = sett, kappa = kappa, n = n, p = p, simu = simu_total,
          seeds = seeds_all
        ),
        beta = beta,
        R_info = list(type = "toeplitz_chol", rho = 0.2, p = p),
        results = res_list
      )
      
      # 文件名
      fpath <- file.path(outdir, sprintf("%s_k=%.2f_n=%d_simu=%d.RData",
                                         sett, kappa, n, simu_total))
      
      # 用 save() 保存为 .RData
      save(bundle, file = fpath, compress = "xz")
      message(sprintf("[SAVED] %s", fpath))
      
      # 释放内存
      rm(res_list, bundle); gc()
    }
  }
}
