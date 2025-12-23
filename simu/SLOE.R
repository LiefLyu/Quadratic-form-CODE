library(pracma)
library(nleqslv)
library(cubature)

rho_prime_logistic <- function(x) 1 / (1 + exp(-x))

f_prime1_logistic <- function(x) -1 / (1 + exp(x))

f_prime0_logistic <- function(x) 1 / (1 + exp(-x))

prox_op <- function(f_prime, lambda, x){
  f <- function(z)  z + lambda * f_prime(z) - x
  uniroot(f, interval = c(-20,20), extendInt = "yes")$root
}

hinge <- function(t) max(t, 0)

integrate2_normal <- function(f, ...){
  integrand <- function(x) f(x) * dnorm(x[1]) * dnorm(x[2])
  
  hcubature(integrand, lowerLimit = c(-8,-8), upperLimit = c(8,8),
            fDim = 1, maxEval = 0, tol = 1e-5,  ...)$integral
}

solve_beta <- function(rho_prime, kappa, gamma0, verbose = FALSE){
  if(solve_kappa(rho_prime, 0, gamma0) < kappa){
    return(-1)
  }else{
    f <- function(beta){
      val <- solve_kappa(rho_prime, beta, gamma0) - kappa
      if(verbose) cat("beta = ", beta, "; diff = ", val, "\n")
      val
    }
    beta_hat <- uniroot(f, interval = c(0, 10), extendInt = "yes", tol = 0.001)$root
    if(verbose) cat("beta_hat = ", beta_hat, "\n")
    return(beta_hat)
  }
}

solve_kappa <- function(rho_prime, beta0, gamma0){
  h <- function(t){
    f1 <- function(x) hinge(t[1] + t[2] * x[1] - x[2])^2 * rho_prime(beta0 + gamma0 * x[1])
    f2 <- function(x) hinge(-t[1] - t[2] * x[1] - x[2])^2 * (1 - rho_prime(beta0 + gamma0 * x[1]))
    integrate2_normal(f1) + integrate2_normal(f2)
  }
  optim(par=c(0, 0), h, method='L-BFGS-B')$val
}

solve_gamma <- function(rho_prime, kappa, beta0, verbose = FALSE){
  if(solve_kappa(rho_prime, beta0, 0) < kappa){
    return(0)
  }else{
    f <- function(gamma) {
      val <- solve_kappa(rho_prime, beta0, gamma) - kappa
      if(verbose) cat("gamma = ", gamma, "; diff = ", val, "\n")
      val
    }
    gamma_hat <- uniroot(f, interval = c(0, 10), extendInt = "yes", tol = 0.001)$root
    if(verbose) cat("gamma_hat = ", gamma_hat, "\n")
    return(gamma_hat)
  }
}

equation_binary_sloe <- function(rho_prime, f_prime1, f_prime0, kappa, eta, beta0, intercept = TRUE){
  function(param){
    if(eta == 0) {
      alpha <-  0; lambda <-  param[1]; sigma <-  param[2];
      if(intercept == 0){
        b <-  0
        
      }else{
        b <- param[3]
        
      }
    }else{
      alpha <- param[1]; lambda <-  param[2]; sigma <-  param[3];
      if(intercept == FALSE){
        b <-  0
        
      }else{
        b <- param[4]
        
      }
    }
    
    # gamma <- sqrt(max(eta^2 - sigma^2, 0.0001)) / alpha
    
    s1_fun <- function(x) beta0 + sqrt(max(eta^2 - sigma^2, 0.0001)) / alpha * x[1]
    s2_fun <- function(x) b + alpha * sqrt(max(eta^2 - sigma^2, 0.0001)) / alpha * x[1] + sigma * x[2]
    
    h1 <- function(x){
      s1 <- s1_fun(x); s2 <- s2_fun(x)
      rho_prime(s1) * (s2 - prox_op(f_prime1, lambda, s2))^2 + (1 - rho_prime(s1)) * (s2 - prox_op(f_prime0, lambda, s2))^2
    }
    
    h2 <- function(x){
      s1 <- s1_fun(x); s2 <- s2_fun(x)
      x[2] * prox_op(f_prime1, lambda, s2) * rho_prime(s1) + x[2] * prox_op(f_prime0, lambda, s2) * (1 - rho_prime(s1))
    }
    
    h3 <- function(x){
      s1 <- s1_fun(x); s2 <- s2_fun(x)
      x[1] * prox_op(f_prime1, lambda, s2) * rho_prime(s1) + x[1] * prox_op(f_prime0, lambda, s2) * (1 - rho_prime(s1))
    }
    
    h4_1 <- function(x){
      s1 <- s1_fun(x); s2 <- s2_fun(x)
      f_prime1(prox_op(f_prime1, lambda, s2)) * rho_prime(s1)
    }
    
    h4_2 <- function(x){
      s1 <- s1_fun(x); s2 <- s2_fun(x)
      f_prime0(prox_op(f_prime0, lambda, s2)) * (1 - rho_prime(s1))
    }
    
    if(eta == 0){
      if(intercept == FALSE){
        c(
          sigma^2 * kappa - integrate2_normal(h1),
          sigma * (1 - kappa) - integrate2_normal(h2)
        )
      }else{
        c(
          sigma^2 * kappa - integrate2_normal(h1),
          sigma * (1 - kappa) - integrate2_normal(h2),
          integrate2_normal(h4_1) + integrate2_normal(h4_2)
        )
      }
    }else if(intercept == FALSE){
      c(
        sigma^2 * kappa - integrate2_normal(h1),
        sigma * (1 - kappa) - integrate2_normal(h2),
        sqrt(max(eta^2 - sigma^2, 0.0001)) / alpha * alpha - integrate2_normal(h3)
      )
    }else{
      c(
        sigma^2 * kappa - integrate2_normal(h1),
        sigma * (1 - kappa) - integrate2_normal(h2),
        sqrt(max(eta^2 - sigma^2, 0.0001)) / alpha * alpha - integrate2_normal(h3),
        integrate2_normal(h4_1) + integrate2_normal(h4_2)
      )
    }
  }
}

find_param_sloe <- function(rho_prime = rho_prime_logistic,
                            f_prime1 = f_prime1_logistic,
                            f_prime0 = f_prime0_logistic,
                            kappa,
                            eta,
                            beta0 = 0,
                            intercept = FALSE,
                            x_init_ = NULL){
  if(eta == 0){ # case of no signal
    if(intercept == FALSE){
      x_init <- c(2, 2)
    }else{
      x_init <- c(2, 2, beta0)
    }
  }else if(intercept == FALSE){
    x_init <- c(2, 2, eta / sqrt(2))
    
    if (kappa < 0.1) {
      x_init <- c(1.0, 7*kappa, 2.9*sqrt(kappa))
    }
    if (kappa == 0.1) {
      x_init <- c(1, 1, 1)
    }
    if (kappa == 0.2) {
      x_init <- c(1.5, 3, 2)
    }
    if(!is.null(x_init_)) x_init <- x_init_
  }else{
    x_init <- c(2, 2, eta / sqrt(2), beta0)

  }
  # Setup system of equations
  f_eq <- equation_binary_sloe(rho_prime, f_prime1, f_prime0,
                               kappa, eta, beta0, intercept)

  
  # sol <- fsolve(f_eq, x_init, J = NULL, maxiter = 10, tol = 1e-4, verbose)
  # sol <- nleqslv(x_init, f_eq, control=list(btol=1e-4))
  sol <- nleqslv(
    x_init, f_eq,
    method  = "Broyden",
    control = list(
      xtol  = 1e-5,   # 解的收敛容差
      ftol  = 1e-5,   # 方程残差容差
      maxit = 1000      # 增加迭代次数
    )
  )
  # print(sol$message)
  
  if(eta == 0){c(1,sol$x)}else{sol$x}
  # if(eta == 0){c(1,sol)}else{sol}
}


g <- function(t) {
  1 / (1 + exp(-t))
}

d_g <- function(t) {
  exp(-t) / (1 + exp(-t))^2
}

# fun_SLOE <- function(X, Y, betahat) {
#   # fit <- glm(Y ~ X + 0, family = binomial)
#   # betahat <- fit$coefficients
#   
#   H <- - tcrossprod(t(X) * as.vector(d_g(crossprod(betahat, t(X)))), t(X))
#   
#   U <- colSums(t(X %*% solve(H)) * t(X))
#   
#   S <- crossprod(betahat, t(X)) + U * (Y - g(crossprod(betahat, t(X)))) / (1 + d_g(crossprod(betahat, t(X))) * U)
#   
#   eta_SLOE <- mean(S^2) - (mean(S))^2
#   
#   return(eta_SLOE)
# }

fun_SLOE <- function(X, Y, beta_hat) {
  n <- nrow(X)
  p <- ncol(X)
  H <- matrix(0, nrow = p, ncol = p)
  
  for (i in 1:n) {
    xi <- X[i, ]
    g_prime_xi <- d_g(t(beta_hat) %*% xi)
    H <- H -  xi %*% t(xi) * c(g_prime_xi)
  }
  
  # 计算 Ui
  U <- rep(NA, n)
  for (i in 1:n) {
    xi <- X[i, ]
    U[i] <- t(xi) %*% solve(H) %*% xi
  }
  
  # 计算 Si
  S <- rep(NA, n)
  for (i in 1:n) {
    xi <- X[i, ]
    g_prime_xi <- d_g(t(beta_hat) %*% xi)
    S[i] <- t(beta_hat) %*% xi + U[i] / (1 + g_prime_xi * U[i]) * (Y[i] - g(t(beta_hat) %*% xi))
  }
  
  # 计算 eta_SLOE^2
  eta_SLOE2 <- (1/n) * sum(S^2) - (1/n)^2 * sum(S)^2
  
  return(eta_SLOE2)
}



fun_SLOE_fast <- function(X, Y, beta_hat) {
  n <- nrow(X)
  
  lin_pred <- X %*% beta_hat
  g_prime_vals <- d_g(lin_pred)
  
  H <- -crossprod(X, X * as.numeric(g_prime_vals))

  H_inv <- solve(H)
  
  U <- rowSums((X %*% H_inv) * X)
  
  g_vals <- g(lin_pred)
  S <- lin_pred + U / (1 + g_prime_vals * U) * (Y - g_vals)
  
  eta_SLOE2 <- var(S) * (n - 1) / n
  
  return(as.numeric(eta_SLOE2))
}
################ test #####################

# n <- 1000
# 
# kappa = 0.1
# ga = sqrt(5)
# 
# p = n * kappa
# 
# # beta = rep(0, p)
# 
# #setup 1
# beta = runif(p)
# beta = beta / norm(beta, "2") * ga
# 
# # setup 2
# # beta = c(runif(p/2), rep(0, p/2))
# # beta = beta / norm(beta, "2") * ga
# 
# R <- diag(rep(1, p))
# # R <- chol(toeplitz(0.5^(0:(p - 1))))
# ga2 <- as.numeric(beta %*% t(R) %*% R %*% beta)
# 
# 
# est <- function(i) {
#   set.seed(i)
# 
#   X <- matrix(rnorm(n * p, 0, 1), n, p) %*% R
# 
#   Y <- rbinom(n, 1, 1 / (1 + exp(- X %*% beta)))
# 
#   eta <- fun_SLOE(X, Y)
# 
#   return(eta)
# }
# 
# eta <- pbsapply(1:1000, est)
# 
# 
# boxplot(eta)

# param <- find_param_sloe(kappa = p/n,
#                     eta = sqrt(eta_real),
#                     intercept = FALSE)
#
# param1 <- find_param(kappa = p/n,
#                     gamma = sqrt(5),
#                     intercept = FALSE)
#
# eta_real <- param1[1]^2 * ga2 + param1[3]^2
#
# mean(eta)
