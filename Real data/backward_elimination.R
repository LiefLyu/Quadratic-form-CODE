library(readr)
library(dplyr)
library(pROC)
library(caret)
library(glmnet) # 新增：加载 glmnet 包
library(tidyr)
source("SLOE.R") # 假设 SLOE.R 文件在您的工作目录中

framingham <- read.csv("framingham.csv")
data0 <- framingham %>% drop_na()

data1 <- data0 %>% filter(male == 1)
data1 <- data1[, -1]

full_model <- glm(TenYearCHD ~ . + 0, data = data1, family = binomial)
step_model <- step(full_model, direction = "backward")

summary(step_model)
selected_vars <- names(coef(step_model))

predictor_names <- setdiff(names(data1), "TenYearCHD")

# 找 selected_vars 在 predictor_names 中的位置
selected_index <- match(selected_vars, predictor_names)
ass <- rep(0, length(predictor_names))

ass[selected_index] <- 1

# 查看结果
ass
#  1 1 0 1 0 0 1 1 0 1 1 1 1 0