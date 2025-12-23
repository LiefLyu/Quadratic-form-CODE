library(R.utils)
library(dplyr)
library(pROC)
library(caret)
library(glmnet) # 新增：加载 glmnet 包
library(tidyr)

source("SLOE.R") # 假设 SLOE.R 文件在您的工作目录中

data_gen <- function(data0_, male_ = 1, n_ = 70, Ten_div = 0.5) {
  
  sample1 <- data0_ %>%
    filter(male == male_ & TenYearCHD == 1) %>%
    slice_sample(n = round(n_*Ten_div), replace = FALSE)
  
  sample2 <- data0_ %>%
    filter(male == male_ & TenYearCHD == 0) %>%
    slice_sample(n = round(n_*(1-Ten_div)), replace = FALSE)
  
  data <- bind_rows(sample1, sample2)
  data <- data[, -1]
  
  cols_to_scale <- setdiff(names(data), c("TenYearCHD"))
  data[cols_to_scale] <- scale(data[cols_to_scale])
  
  return(data)
}

################## no cv lasso ########################

simu <- function(i, ass_ = c(1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0), 
                 data0_, male_ = 1, Sigmahat = NULL, n_ = 70, Ten_div = 0.5,
                 lambda_lasso = 0.1) { 
  
  est_NA <- c(est_null = NA, est_1 = NA, est_all = NA, 
              naive_null = NA, naive_1 = NA, naive_all = NA,
              lasso_null = NA, lasso_1 = NA, lasso_all = NA)
  
  tryCatch({
    withTimeout({
      
      set.seed(i)
      
      while (TRUE) {
        
        data <- data_gen(data0_ = data0_, male_ = male_, n_ = n_, Ten_div = Ten_div)
        
        formula <- as.formula("TenYearCHD ~ . + 0")
        X_m <- model.matrix(formula, data = data)
        y_m <- data$TenYearCHD
        p = NCOL(X_m)
        n <- NROW(X_m)
        
        if (n == 0) next
        
        # if (sum(data$BPMeds) > 3 | sum(data$diabetes) > 3) next
        
        glm_success <- TRUE
        logistic_model <- NULL
        
        withCallingHandlers(
          {
            logistic_model <- glm(y_m ~ X_m + 0, family = binomial)
          },
          warning = function(w) {
            # cat(sprintf("glm warning in iteration %d, retrying... (Warning: %s)\n", i, conditionMessage(w)))
            glm_success <<- FALSE
            invokeRestart("muffleWarning")
          }
        )
        
        if (!glm_success) {
          next
        }
        
        break
      }
      
      lasso_model <- glmnet(X_m, y_m, alpha = 1, family = "binomial", intercept = FALSE)

      beta_lasso <- as.vector(coef(lasso_model, s = lambda_lasso))[-1]
      
      beta_hat <- coef(logistic_model)
      
      etahat <- fun_SLOE_fast(X_m, y_m, beta_hat)
      
      param <- find_param_sloe(kappa = p/n,
                               eta = sqrt(etahat),
                               intercept = FALSE,
                               x_init_ = c(1, 0.5, 0.5))
      
      kappa = p/n
      if(is.null(Sigmahat)) Sigmahat <- cov(X_m)
      
      if (inherits(try(solve(Sigmahat), silent = TRUE), "try-error")) {
        return(est_NA)
      }
      Sigma.inv = solve(Sigmahat)
      
      if (length(ass_) != p) {
        ass_ <- c(1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0)
        if (length(ass_) != p) stop("Length of ass_ does not match number of predictors.")
      }
      
      A_0 <- diag(1 - ass_)
      est_null <- (as.numeric(t(beta_hat) %*% A_0 %*% beta_hat) - (1-kappa) * param[3]^2 * sum(diag(A_0 %*% Sigma.inv)) / p ) / param[1]^2
      naive_null <- as.numeric(t(beta_hat) %*% A_0 %*% beta_hat) 
      lasso_null <- as.numeric(t(beta_lasso) %*% A_0 %*% beta_lasso) 
      
      A_1 <- diag(ass_)
      est_1 <- (as.numeric(t(beta_hat) %*% A_1 %*% beta_hat) - (1-kappa) * param[3]^2 * sum(diag(A_1 %*% Sigma.inv)) / p ) / param[1]^2
      naive_1 <- as.numeric(t(beta_hat) %*% A_1 %*% beta_hat) 
      lasso_1 <- as.numeric(t(beta_lasso) %*% A_1 %*% beta_lasso) 
      
      A_all <- diag(p)
      est_all <- (as.numeric(t(beta_hat) %*% A_all %*% beta_hat) - (1-kappa) * param[3]^2 * sum(diag(A_all %*% Sigma.inv)) / p ) / param[1]^2
      naive_all <- as.numeric(t(beta_hat) %*% A_all %*% beta_hat) 
      lasso_all <- as.numeric(t(beta_lasso) %*% A_all %*% beta_lasso) 
      
      return(c(est_null = est_null, est_1 = est_1, est_all = est_all,
               naive_null = naive_null, naive_1 = naive_1, naive_all = naive_all,
               lasso_null = lasso_null, lasso_1 = lasso_1, lasso_all = lasso_all))
      
    }, timeout = 300, onTimeout = "error")
    
  }, warning = function(w) {
    cat(sprintf("Outer warning captured in iteration %d: %s\n", i, conditionMessage(w)))
    return(est_NA)
    
  }, error = function(e) {
    if (inherits(e, "TimeoutException")) {
      cat(sprintf("Timeout occurred in iteration %d\n", i))
    } else {
      cat(sprintf("Error captured in iteration %d: %s\n", i, conditionMessage(e)))
    }
    return(est_NA)
  })
}


###################### simu cv lasso ############################

simu_cv <- function(i, ass_ = c(1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0), 
                 data0_, male_ = 1, Sigmahat = NULL, n_ = 70, Ten_div = 0.4) {
  
  est_NA <- c(est_null = NA, est_1 = NA, est_all = NA, 
              naive_null = NA, naive_1 = NA, naive_all = NA,
              lasso_null = NA, lasso_1 = NA, lasso_all = NA)
  
  tryCatch({
    withTimeout({
      
      set.seed(i)
      
      while (TRUE) {
        
        data <- data_gen(data0_ = data0_, male_ = male_, n_ = n_, Ten_div = Ten_div)
        
        formula <- as.formula("TenYearCHD ~ . + 0")
        X_m <- model.matrix(formula, data = data)
        y_m <- data$TenYearCHD
        p = NCOL(X_m)
        n <- NROW(X_m)
        
        if (n == 0) {
          next
        }
        
        glm_success <- TRUE
        logistic_model <- NULL
        
        # --- 修改点 ---
        # 现在，只要 glm() 产生任何警告，都会触发重试
        withCallingHandlers(
          {
            logistic_model <- glm(y_m ~ X_m + 0, family = binomial)
          },
          warning = function(w) {
            # 不再检查特定的警告内容。任何警告都会将成功标志设为FALSE。
            cat(sprintf("glm warning in iteration %d, retrying... (Warning: %s)\n", i, conditionMessage(w))) # 打印警告信息以供调试
            glm_success <<- FALSE
            invokeRestart("muffleWarning")
          }
        )
        
        if (!glm_success) {
          next
        }
        
        break
      }
      
      cv_lasso_model <- cv.glmnet(X_m, y_m, alpha = 1, family = "binomial", intercept = FALSE)
      
      lambda_best <- cv_lasso_model$lambda.min
      beta_lasso <- as.vector(coef(cv_lasso_model, s = lambda_best))[-1]
      
      beta_hat <- coef(logistic_model)
      
      etahat <- fun_SLOE_fast(X_m, y_m, beta_hat)
      
      param <- find_param_sloe(kappa = p/n,
                               eta = sqrt(etahat),
                               intercept = FALSE,
                               x_init_ = c(1, 0.5, 0.5))
      
      kappa = p/n
      if(is.null(Sigmahat)) Sigmahat <- cov(X_m)
      
      if (inherits(try(solve(Sigmahat), silent = TRUE), "try-error")) {
        return(est_NA)
      }
      
      if (length(ass_) != p) {
        ass_ <- c(1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0)
        if (length(ass_) != p) stop("Length of ass_ does not match number of predictors.")
      }
      
      A_0 <- diag(1 - ass_)
      est_null <- (as.numeric(t(beta_hat) %*% A_0 %*% beta_hat) - (1-kappa) * param[3]^2 * sum(diag(A_0 %*% solve(Sigmahat))) / p ) / param[1]^2
      naive_null <- as.numeric(t(beta_hat) %*% A_0 %*% beta_hat) 
      lasso_null <- as.numeric(t(beta_lasso) %*% A_0 %*% beta_lasso) 
      
      A_1 <- diag(ass_)
      est_1 <- (as.numeric(t(beta_hat) %*% A_1 %*% beta_hat) - (1-kappa) * param[3]^2 * sum(diag(A_1 %*% solve(Sigmahat))) / p ) / param[1]^2
      naive_1 <- as.numeric(t(beta_hat) %*% A_1 %*% beta_hat) 
      lasso_1 <- as.numeric(t(beta_lasso) %*% A_1 %*% beta_lasso) 
      
      A_all <- diag(p)
      est_all <- (as.numeric(t(beta_hat) %*% A_all %*% beta_hat) - (1-kappa) * param[3]^2 * sum(diag(A_all %*% solve(Sigmahat))) / p ) / param[1]^2
      naive_all <- as.numeric(t(beta_hat) %*% A_all %*% beta_hat) 
      lasso_all <- as.numeric(t(beta_lasso) %*% A_all %*% beta_lasso) 
      
      return(c(est_null = est_null, est_1 = est_1, est_all = est_all,
               naive_null = naive_null, naive_1 = naive_1, naive_all = naive_all,
               lasso_null = lasso_null, lasso_1 = lasso_1, lasso_all = lasso_all))
      
    }, timeout = 60, onTimeout = "error")
    
  }, warning = function(w) {
    # 这个外层 warning handler 现在不太可能被 glm 的警告触发，
    # 因为它们已经被内层的 withCallingHandlers 处理掉了。
    # 但保留它以处理其他部分的警告是好的实践。
    cat(sprintf("Outer warning captured in iteration %d: %s\n", i, conditionMessage(w)))
    return(est_NA)
    
  }, error = function(e) {
    if (inherits(e, "TimeoutException")) {
      cat(sprintf("Timeout occurred in iteration %d\n", i))
    } else {
      cat(sprintf("Error captured in iteration %d: %s\n", i, conditionMessage(e)))
    }
    return(est_NA)
  })
}
