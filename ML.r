rm(list = ls())
library(caret)
library(randomForest)
library(pROC)

data <- read.csv("~/trp_L_strain.csv")
data <- data[, -2]  
data$Group <- factor(data$Group)  

# one metabolite
metabolite_feature <- "L.Tryptophan" 

set.seed(42)

normal_data <- data[1:92, ]
disease_data <- data[93:219, ]

normal_train_idx <- createDataPartition(normal_data$Group, p = 0.6, list = FALSE)
normal_train <- normal_data[normal_train_idx, ]
normal_test <- normal_data[-normal_train_idx, ]

disease_train_idx <- createDataPartition(disease_data$Group, p = 0.6, list = FALSE)
disease_train <- disease_data[disease_train_idx, ]
disease_test <- disease_data[-disease_train_idx, ]

train_data <- rbind(normal_train, disease_train)
test_data <- rbind(normal_test, disease_test)

# RFE
train_control <- trainControl(
  method = "repeatedcv",
  number = 5,
  repeats = 3,
  verboseIter = TRUE
)

all_features <- setdiff(names(train_data), c("Group", metabolite_feature))

# no more than 14 strains
tune_grid <- expand.grid(.mtry = seq(1, min(10, length(all_features)), by = 1))
rfe_control <- rfeControl(functions = rfFuncs, method = "cv", number = 5)
sizes <- seq(1, min(14, length(all_features)), by = 1)

# RFE 
rfe_result <- rfe(
  x = train_data[, all_features],
  y = train_data$Group,
  sizes = sizes,
  rfeControl = rfe_control,
  trControl = train_control,
  tuneGrid = tune_grid
)

top_rfe_features <- predictors(rfe_result)[1:min(14, length(predictors(rfe_result)))]
best_features <- union(top_rfe_features, metabolite_feature)

cat("final number：", length(best_features), "\n")
print(best_features)

# train
final_rf_model <- randomForest(
  Group ~ ., 
  data = train_data[, c(best_features, "Group")],
  mtry = rfe_result$bestSubset,
  ntree = 100,
  nodesize = 15
)

importance_scores <- importance(final_rf_model)
importance_df <- data.frame(
  Feature = rownames(importance_scores),
  Importance = importance_scores[, "MeanDecreaseGini"]
)
importance_df <- importance_df[order(-importance_df$Importance), ]

cat("\nrank：\n")
print(importance_df)

write.csv(importance_df, "L_Tryptophan_feature_importance.csv", row.names = FALSE)

# ROC 
y_pred_prob_train <- predict(final_rf_model, newdata = train_data[, best_features], type = "prob")[, 2]
roc_result_train <- roc(train_data$Group, y_pred_prob_train)

y_pred_prob_test <- predict(final_rf_model, newdata = test_data[, best_features], type = "prob")[, 2]
roc_result_test <- roc(test_data$Group, y_pred_prob_test)

fixed_marker <- data[, metabolite_feature]
roc_result_fixed_marker <- roc(data$Group, fixed_marker)

auc_train <- auc(roc_result_train)
auc_test <- auc(roc_result_test)
auc_fixed_marker <- auc(roc_result_fixed_marker)

cat("train AUC: ", auc_train, "\n")
cat("test AUC: ", auc_test, "\n")
cat("fix marker AUC: ", auc_fixed_marker, "\n")


# tiff("L_Tryptophan_ROC.tiff", width = 2000, height = 1000, res = 300)
pdf("L_Tryptophan_ROC.pdf", width = 10, height = 5)
# tiff("roc_plot.tiff", width = 2000, height = 1000, res = 300)
pdf("roc_plot.pdf", width = 8, height = 8)

plot(roc_result_train, main = "L.Tryptophan ROC Curves", 
     col = "red", lwd = 2, legacy.axes = TRUE,
     xlim = rev(c(0, 1)), ylim = c(0, 1))

plot(roc_result_test, col = "green", lwd = 2, add = TRUE)
plot(roc_result_fixed_marker, col = "purple", lwd = 2, add = TRUE)

legend("bottomright", legend = c(
  paste("Training (AUC =", round(auc_train, 2), ")"),
  paste("Test (AUC =", round(auc_test, 2), ")"),
  paste("Fixed Marker (AUC =", round(auc_fixed_marker, 2), ")")
), col = c("red", "green", "purple"), lwd = 2, bty = "n", cex = 1.0)

dev.off()