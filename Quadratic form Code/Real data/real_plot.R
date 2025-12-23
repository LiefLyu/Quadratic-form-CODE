library(readr)
library(pbapply)
library(pbmcapply)
library(ggplot2)
library(scales)
library(latex2exp)
library(ggtext)
# install.packages("R.utils")

source("simu_func.R")

framingham <- read.csv("framingham.csv")
data0 <- framingham %>% drop_na()

n_values <- c(100, 140, 200, 280, 400, 600) 
replications <- 500 

# Your other fixed parameters
ass <- c(1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0)
# Make sure the 'data0' object and the 'simu' function are loaded in your environment


# --- Step 3: Main Simulation Loop ---
# This list will store the results data.frame for each n_
all_results_list <- lapply(n_values, function(current_n) {
  
  cat(sprintf("\n--- Running simulations for n = %d ---\n", current_n))
  
  # THIS IS THE KEY CHANGE FOR MAC:
  # Use pbmclapply directly. It will work and show a progress bar.
  results_list <- pbmclapply(
    1:replications, 
    simu,                  
    ass_ = ass, data0_ = data0, male_ = 1, 
    Sigmahat = NULL, Ten_div = 0.5, lambda_lasso = 0.05,
    n_ = current_n,
    mc.cores = parallel::detectCores() # Use all available cores
  )
  
  # Combine results from the 100 reps
  results_df <- do.call(rbind, results_list)
  results_df <- as.data.frame(results_df)
  
  # Add a column to identify which 'n' these results belong to
  results_df$n <- current_n
  
  return(results_df)
})

# --- Step 4: Combine all data frames into one ---
all_results_df <- bind_rows(all_results_list)

saveRDS(all_results_df, "framingham/simulation_results.rds")

# --- Step 5: Data Preparation and Plotting (This part is identical) ---

all_results_df <- readRDS("simulation_results.rds")
# Reshape the data from wide to long
results_long <- all_results_df %>%
  pivot_longer(
    cols = -n, # Pivot all columns except for 'n'
    names_to = c("method", "quantity"), # Create two new columns
    names_sep = "_", # Split the original column names at the underscore "_"
    values_to = "value"
  )

# --- 5. Clean Up Labels for a Publication-Quality Plot ---
# Here we recode the internal names to be more descriptive
# Make sure your original 'results_long' data frame is available
results_long_final <- results_long %>%
  mutate(
    n = factor(n),
    method = case_when(
      method == "est"   ~ "Corr<sup>+</sup>",
      method == "naive" ~ "MLE plug-in",
      method == "lasso" ~ "Lasso plug-in",
      TRUE ~ method
    ),
    # 设定图例（填充）的顺序：MLE plug-in, Corr^+, Lasso plug-in
    method = factor(method, levels = c("MLE plug-in", "Corr<sup>+</sup>", "Lasso plug-in")),
    quantity = case_when(
      quantity == "null" ~ "Null Coefficients",
      quantity == "1"    ~ "Signal Coefficients",
      quantity == "all"  ~ "All Coefficients",
      TRUE ~ quantity
    ),
    quantity = factor(quantity, levels = c(
      "Null Coefficients", 
      "Signal Coefficients", 
      "All Coefficients"
    ))
  )

results_long_final_no100 <- results_long_final %>% 
  filter(n != '100') %>% 
  drop_na()

results_long_final_no100_positive <- results_long_final_no100 %>% 
  mutate(value = pmax(value, 0))

new_facet_labels <- c(
  "Null Coefficients"   = TeX("$|| \\beta_{null} ||_2^2$"),
  "Signal Coefficients" = TeX("$|| \\beta_{sig} ||_2^2$"),
  "All Coefficients"    = TeX("$|| \\beta ||_2^2$")
)

new_facet_labels_ggtext <- c(
  "Null Coefficients"   = "(a)  ||β<sub>null</sub>||<sup>2</sup>",
  "Signal Coefficients" = "(b)  ||β<sub>sig</sub>||<sup>2</sup>",
  "All Coefficients"    = "(c)  ||β||<sup>2</sup>"
)

# Now create the plot with the new data frame
pp <- ggplot(results_long_final_no100_positive, aes(x = n, y = value, fill = method)) +
  
  # Create side-by-side boxplots. position_dodge() is the key here.
  geom_boxplot(position = position_dodge(width = 0.9), outlier.size = 0.5) +
  
  # Use facet_wrap to create the 3 main plot panels.
  # scales = "free_y" allows each panel to have its own y-axis range.
  facet_wrap(~ quantity, scales = "free_y",
             labeller = labeller(quantity = new_facet_labels_ggtext)) +
  scale_y_continuous(trans = scales::log1p_trans()) +
  coord_cartesian(ylim = c(0, 3)) + 
  # Improve the plot's labels
  labs(
    # title = "Comparison of Estimators by Sample Size",
    # subtitle = "Estimates are faceted by coefficient set (Null, Signal, All)",
    x = "Sample Size (n)",
    y = "Estimate Value",
    fill = "method" # This sets the title for the legend
  ) +
  
  # Use a minimal theme and adjust the legend
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    strip.text = element_markdown(size = 14, face = "bold"), # Make facet titles bold
    # axis.text.x = element_text(angle = 45, hjust = 1) # Angle x-axis text if needed
    legend.title = element_text(size = 14), 
    legend.key.width = unit(2, "cm"),
    # 调整图例项目的字体大小
    legend.text = element_markdown(size = 12) 
  )

pp

ggsave("Real_data.pdf", 
       plot = pp, 
       width = 12, 
       height = 4, 
       units = "in",
       device = cairo_pdf)
#######################################################
# ass <- c(1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0)
# 
# # simu(i = 11, ass_ = ass, male_ = 1, data0_ = data0, Sigmahat = NULL, n_ = 280, Ten_div = 0.4, lambda_lasso = 0.1)
# # results <- pbsapply(1:100, simu, ass_ = ass, male_ = 1, data0_ = data0, Sigmahat = NULL, n_ = 280, Ten_div = 0.4, lambda_lasso = 0.1)
# 
# results_list <- pbmclapply(1:100, simu, ass_ = ass, data0_ = data0, male_ = 1, 
#                            Sigmahat = NULL, n_ = 100, Ten_div = 0.5, lambda_lasso = 0.05,
#                            mc.cores = parallel::detectCores())
# results_df <- do.call(rbind, results_list)
# results_df <- as.data.frame(results_df)
# summary(results_df)
# 
# results_df <- as.data.frame(t(results))
# summary(results_df)
# 
# rowMeans(results, na.rm = TRUE)
# boxplot(results[3,])
# View(results)
# median(results[3,], na.rm = TRUE)
