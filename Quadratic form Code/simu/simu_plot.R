## ====== Common setup ======
library(ggplot2)
library(latex2exp)
library(patchwork)

## 参数（按你的新设定）
simu_total <- 2000L
nn         <- c(1000L, 2000L, 4000L)
kappas     <- c(0, 0.1, 0.2)

sections   <- c("norm", "k", "S", "Sc")  # 三个子节对应的节名
setting    <- c("dense", "sparse")
out_dir    <- "plot"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## 方法键与图例标签/颜色
method_keys   <- c("MLE_naive", "MLE_corrected", "MLE_Lasso")
labels_trunc  <- c(TeX("MLE plug-in"), TeX("Corr$^{+}$"), TeX("Lasso plug-in"))
labels_plain  <- c(TeX("MLE plug-in"), TeX("Corr"),       TeX("Lasso plug-in"))
method_colors <- c("#1B9E77", "#D95F02", "#7570B3")  # 绿 / 橙 / 紫
names(method_colors) <- method_keys

## -------- 工具函数：读取一个 (sect, sett, kappa) 下三个 n 的结果，拼成长表 --------
## 路径： simu_est/<sect>_<sett>_k=%.2f_n=%d_simu=%d.RData
collect_results_one_kappa <- function(sect, sett, kappa, nn, simu_total) {
  simu_result <- vector("list", length(nn))
  for (i in seq_along(nn)) {
    n   <- nn[i]
    f   <- file.path("simu_est",
                     sprintf("%s_%s_k=%.2f_n=%d_simu=%d.RData",
                             sect, sett, kappa, n, simu_total))
    if (!file.exists(f)) stop("File not found: ", f)
    load(f)  # expects object S
    # 简单健壮性检查
    if (!is.list(S) || length(S) < 4) stop("Object S malformed in: ", f)
    simu_result[[i]] <- S
  }
  
  # 逐 n、逐方法构建 block，使用“实际长度”
  blocks <- list()
  idx <- 1L
  for (i in seq_along(nn)) {
    n_lab <- paste0("n = ", nn[i])
    
    for (mkey in c("MLE_naive","MLE_corrected","MLE_Lasso")) {
      m_idx <- switch(mkey, MLE_naive=1L, MLE_corrected=2L, MLE_Lasso=3L)
      vec   <- simu_result[[i]][[m_idx]]
      if (is.null(vec)) stop("Missing S[[", m_idx, "]] for n=", nn[i])
      
      # 如果长度和 simu_total 不同，给个提醒（不会中断）
      if (length(vec) != simu_total) {
        warning(sprintf("Length of S[[%d]] at n=%d (k=%.2f, %s-%s) is %d (!= %d).",
                        m_idx, nn[i], kappa, sect, sett, length(vec), simu_total))
      }
      
      blocks[[idx]] <- data.frame(
        n      = factor(rep(n_lab, length(vec)), levels = paste0("n = ", nn)),
        results = as.numeric(vec),
        method  = factor(mkey, levels = c("MLE_naive","MLE_corrected","MLE_Lasso"))
      )
      idx <- idx + 1L
    }
  }
  
  com_df <- do.call(rbind, blocks)
  
  # 真值：沿用第一个文件中的 S[[4]]
  truth_val <- as.numeric(simu_result[[1]][[4]])
  list(df = com_df, truth = truth_val)
}
## -------- 工具函数：单图（一个 κ）箱线图，带对数坐标、裁剪与标题 --------
make_boxplot <- function(com_df, truth_val, ylab_tex, # 不在函数里放 kappa 文本，title 直接传入
                         title_tex,
                         trans = c("log10","log1p"),
                         use_trunc_label = TRUE,
                         cap_high_q = 1,
                         legend_text_size = 12,
                         show_ylab = TRUE) {
  trans <- match.arg(trans)
  df_plot <- com_df
  if (trans == "log10") df_plot$results[df_plot$results <= 0] <- NA
  if (trans == "log1p") df_plot$results[df_plot$results <  0] <- NA
  
  # 计算上截点（全局 99.5% 分位）
  vals <- df_plot$results[is.finite(df_plot$results)]
  if (trans == "log10") vals <- vals[vals > 0]
  if (trans == "log1p") vals <- vals[vals >= 0]
  hi <- as.numeric(stats::quantile(vals, probs = cap_high_q, na.rm = TRUE))
  if (is.finite(truth_val)) hi <- max(hi, truth_val)
  
  p <- ggplot(df_plot, aes(x = n, y = results, fill = method)) +
    geom_boxplot() +
    ggtitle(latex2exp::TeX(title_tex)) +
    geom_hline(yintercept = truth_val, color = "#e20612") +
    theme_bw() +
    xlab(NULL) +
    { if (show_ylab) ylab(latex2exp::TeX(ylab_tex)) else ylab(NULL) } +
    {
      if (use_trunc_label)
        scale_fill_manual(breaks = method_keys, labels = labels_trunc, values = method_colors)
      else
        scale_fill_manual(breaks = method_keys, labels = labels_plain, values = method_colors)
    } +
    theme(
      legend.title = element_blank(),
      legend.text  = element_text(size = legend_text_size),
      legend.key.width = unit(2, "cm"),   # 控制每个图例 key 的宽度
      plot.title   = element_text(hjust = 0)
    )
  
  if (trans == "log10") {
    p + scale_y_log10(limits = c(NA, hi), oob = scales::oob_censor)
  } else {
    p + scale_y_continuous(trans = "log1p", limits = c(NA, hi), oob = scales::oob_censor)
  }
}

## ---------------- Subsection 1：估计 ||β||^2（sect = "norm"） ----------------
## 两行：上=dense，下=sparse；每行三图：κ=0.01, 0.1, 0.2
## 采用 Corr^{+}（trunc）图例标签
make_row_norm <- function(sett, kappas, nn, simu_total, cap_high_q = 1) {
  plots <- vector("list", length(kappas))
  for (i in seq_along(kappas)) {
    kappa <- kappas[i]
    res <- collect_results_one_kappa(sect = "norm", sett = sett,
                                     kappa = kappa, nn = nn, simu_total = simu_total)
    # 标题：dense 用 (a-i)/(b-i)/(c-i)，sparse 用 (a-ii)/(b-ii)/(c-ii)
    tag_letter <- c("a","b","c")[i]
    row_tag    <- if (sett == "dense") "i" else "ii"
    # title_here <- sprintf("(%s$-$%s) $\\kappa=%.2f$, %s",
    #                       tag_letter, row_tag, kappa, sett)
    title_here <- sprintf("(%s$-$%s) $\\kappa=%.2f$, %s $\\beta$",
                          tag_letter, row_tag, kappa, sett)
    plots[[i]] <- make_boxplot(
      res$df, res$truth,
      ylab_tex = "$\\widehat{\\|\\beta\\|^2}$",
      title_tex = title_here,
      trans = "log1p",
      cap_high_q = cap_high_q,
      use_trunc_label = TRUE,
      show_ylab = (i == 1)
    )
  }
  plots
}

norm_dense_row  <- make_row_norm("dense",  kappas, nn, simu_total, cap_high_q = 1)
norm_sparse_row <- make_row_norm("sparse", kappas, nn, simu_total, cap_high_q = 1)
fig_bnorm <- (norm_dense_row[[1]] + norm_dense_row[[2]] + norm_dense_row[[3]]) /
  (norm_sparse_row[[1]] + norm_sparse_row[[2]] + norm_sparse_row[[3]]) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(filename = file.path(out_dir, "fig_bnorm.pdf"),
       plot = fig_bnorm, width = 12, height = 6, device = cairo_pdf)

## ---------------- Subsection 2：估计 β_k^2（sect = "k"） ----------------
## 两行：上=dense，下=sparse；每行三图；Corr（不加 +）
make_row_bk <- function(sett, kappas, nn, simu_total) {
  plots <- vector("list", length(kappas))
  for (i in seq_along(kappas)) {
    kappa <- kappas[i]
    res <- collect_results_one_kappa(sect = "k", sett = sett,
                                     kappa = kappa, nn = nn, simu_total = simu_total)
    tag_letter <- c("a","b","c")[i]
    row_tag    <- if (sett == "dense") "i" else "ii"
    # title_here <- sprintf("(%s$-$%s) $\\kappa=%.2f$, %s",
    #                       tag_letter, row_tag, kappa, sett)
    title_here <- sprintf("(%s$-$%s) $\\kappa=%.2f$, %s $\\beta$",
                          tag_letter, row_tag, kappa, sett)
    plots[[i]] <- make_boxplot(
      res$df, res$truth,
      ylab_tex = "$\\widehat{\\beta_1^2}$",
      title_tex = title_here,
      trans = "log1p",
      cap_high_q = 1,
      use_trunc_label = FALSE,  # Corr
      show_ylab = (i == 1)
    )
  }
  plots
}

bk_dense_row  <- make_row_bk("dense",  kappas, nn, simu_total)
bk_sparse_row <- make_row_bk("sparse", kappas, nn, simu_total)
fig_bk <- (bk_dense_row[[1]] + bk_dense_row[[2]] + bk_dense_row[[3]]) /
  (bk_sparse_row[[1]] + bk_sparse_row[[2]] + bk_sparse_row[[3]]) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(filename = file.path(out_dir, "fig_bk.pdf"),
       plot = fig_bk, width = 12, height = 6, device = cairo_pdf)

## -------- Subsection 3（拆分成两个文件）：||β_S||^2 与 ||β_{S^c}||^2 --------
## 仍然使用 Corr^{+} 图例标签，三列对应 κ=0.01, 0.1, 0.2
## 每个文件做两行：上=dense，下=sparse；左图显示 y 轴标签

make_row_group <- function(sect, sett, kappas, nn, simu_total, ylab_tex) {
  plots <- vector("list", length(kappas))
  for (i in seq_along(kappas)) {
    kappa <- kappas[i]
    res <- collect_results_one_kappa(
      sect = sect, sett = sett, kappa = kappa, nn = nn, simu_total = simu_total
    )
    tag_letter <- c("a","b","c")[i]
    row_tag <- if (sett == "dense") "i" else "ii"
    sect_tex <- if (sect == "S") "$S$" else "$S^{c}$"
    # title_here <- sprintf("(%s$-$%s) $\\kappa=%.2f$, %s, %s",
    #                       tag_letter, row_tag, kappa, sett, sect_tex)
    title_here <- sprintf("(%s$-$%s) $\\kappa=%.2f$, %s $\\beta$, %s",
                          tag_letter, row_tag, kappa, sett, sect_tex)
    
    plots[[i]] <- make_boxplot(
      res$df, res$truth,
      ylab_tex   = ylab_tex,
      title_tex  = title_here,
      trans      = "log1p",        # S^c 真值可能为 0，用 log1p 稳妥
      use_trunc_label = TRUE,      # Corr^{+}
      show_ylab  = (i == 1)
    )
  }
  plots
}

## ===== 文件 1：只画 ||β_S||^2（dense + sparse）=====
row_dense_S  <- make_row_group("S",  "dense",  kappas, nn, simu_total,
                               ylab_tex = "$\\widehat{\\|\\beta_S\\|^2}$")
row_sparse_S <- make_row_group("S",  "sparse", kappas, nn, simu_total,
                               ylab_tex = "$\\widehat{\\|\\beta_S\\|^2}$")

fig_S_only <- (row_dense_S[[1]] + row_dense_S[[2]] + row_dense_S[[3]]) /
  (row_sparse_S[[1]] + row_sparse_S[[2]] + row_sparse_S[[3]]) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(filename = file.path(out_dir, "fig_bS_dense_sparse.pdf"),
       plot = fig_S_only, width = 12, height = 6, device = cairo_pdf)

## ===== 文件 2：只画 ||β_{S^c}||^2（dense + sparse）=====
row_dense_Sc  <- make_row_group("Sc", "dense",  kappas, nn, simu_total,
                                ylab_tex = "$\\widehat{\\|\\beta_{S^c}\\|^2}$")
row_sparse_Sc <- make_row_group("Sc", "sparse", kappas, nn, simu_total,
                                ylab_tex = "$\\widehat{\\|\\beta_{S^c}\\|^2}$")

fig_Sc_only <- (row_dense_Sc[[1]] + row_dense_Sc[[2]] + row_dense_Sc[[3]]) /
  (row_sparse_Sc[[1]] + row_sparse_Sc[[2]] + row_sparse_Sc[[3]]) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(filename = file.path(out_dir, "fig_bSc_dense_sparse.pdf"),
       plot = fig_Sc_only, width = 12, height = 6, device = cairo_pdf)

