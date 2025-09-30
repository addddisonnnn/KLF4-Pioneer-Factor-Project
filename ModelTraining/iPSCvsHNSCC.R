library(tidyverse)
library(caret)
library(randomForest)
library(pROC)
library(knitr)

# --- 2. Load the Data ---
# Make sure your CSV file 'Combined_Human_Histone_Data.csv' is in your R working directory,
# or provide the full path to the file.
base_path <- "/Users/addisonyam/Documents/KLF4/"
df <- read.csv(file.path(base_path, "Combined_Human_Histone_Data.csv"))

# --- 3. Data Preprocessing (Common Steps) ---

# Define all original cell line columns
original_cell_line_cols <- c("BJ_epi", "BJ_iPSC", "H9_ESC", "iPSC", "LIS49_hESC",
                             "MCF7_endo", "MCF7_induc", "U87", "HN_SCC")

# Handle 'Strand' column: Convert to factor
df$Strand <- as.factor(df$Strand)

# Handle potential empty strings in numeric columns (e.g., motif scores)
numeric_score_cols <- c("Klf4.motif.score", "Total.motif.score", "FOSL2.motif.score",
                        "FOXC2.motif.score", "OTX2.motif.score", "POU6F2.motif.score",
                        "TP53.motif.score", "ZBED4.motif.score", "ZBTB6.motif.score",
                        "ZNF213.motif.score", "ZNF384.motif.score", "ZNF460.motif.score")

for (col in numeric_score_cols) {
  if (col %in% names(df)) {
    df[[col]] <- as.numeric(ifelse(df[[col]] == "", NA, df[[col]]))
    df[[col]][is.na(df[[col]])] <- 0
  }
}

# --- 4. Create Exclusive Binding Subsets for BJ_iPSC and HN_SCC ---

message("\n--- Creating exclusive binding subsets for BJ_iPSC and HN_SCC ---")

# BJ_iPSC exclusive peaks: KLF4 binds in BJ_iPSC AND NOT in HN_SCC
bj_ipsc_exclusive_peaks <- df %>%
  filter(BJ_iPSC == 1 & HN_SCC == 0) %>%
  mutate(Binding_Specificity = "BJ_iPSC_exclusive")

# HN_SCC exclusive peaks: KLF4 binds in HN_SCC AND NOT in BJ_iPSC
hn_scc_exclusive_peaks <- df %>%
  filter(HN_SCC == 1 & BJ_iPSC == 0) %>%
  mutate(Binding_Specificity = "HN_SCC_exclusive")

message(paste0("Number of BJ_iPSC exclusive peaks: ", nrow(bj_ipsc_exclusive_peaks)))
message(paste0("Number of HN_SCC exclusive peaks: ", nrow(hn_scc_exclusive_peaks)))

# --- 5. Combine Subsets and Prepare for Modeling ---

# Check if either subset is empty
if (nrow(bj_ipsc_exclusive_peaks) == 0 || nrow(hn_scc_exclusive_peaks) == 0) {
  stop("Error: One or both exclusive peak sets are empty.
        Cannot train a classification model without both classes.
        This means there are no (or very few) KLF4 peaks exclusively binding
        in one cell line (BJ_iPSC or HN_SCC) but not the other based on your current data.")
}

# Combine the exclusive peak sets
combined_exclusive_data <- bind_rows(bj_ipsc_exclusive_peaks, hn_scc_exclusive_peaks)

# Convert the new target variable to a factor
combined_exclusive_data$Binding_Specificity <- as.factor(combined_exclusive_data$Binding_Specificity)

# Define features for the model
# Exclude identifiers (chrom, start, end)
# Exclude all original cell line binding columns, as 'BJ_iPSC' and 'HN_SCC' were used to define the target.
features_to_exclude_from_model <- c("chrom", "start", "end", original_cell_line_cols)

model_data <- combined_exclusive_data %>%
  select(-all_of(features_to_exclude_from_model))

# Ensure all remaining columns are suitable for the model
if (any(is.na(model_data))) {
  message("Warning: NAs still present in the combined exclusive data. Removing rows with NAs.")
  model_data <- na.omit(model_data)
}

# --- 6. Train-Test Split ---
# Set seed for reproducibility
set.seed(123)

# Create a stratified split to ensure similar proportions of Binding_Specificity in both sets
trainIndex <- createDataPartition(model_data$Binding_Specificity, p = 0.7, list = FALSE, times = 1)
train_data <- model_data[trainIndex, ]
test_data <- model_data[-trainIndex, ]

message(paste0("Training data rows for differential binding: ", nrow(train_data)))
message(paste0("Test data rows for differential binding: ", nrow(test_data)))
message(paste0("Proportion of BJ_iPSC_exclusive in training data: ", round(mean(train_data$Binding_Specificity == "BJ_iPSC_exclusive"), 3)))
message(paste0("Proportion of HN_SCC_exclusive in training data: ", round(mean(train_data$Binding_Specificity == "HN_SCC_exclusive"), 3)))


# --- 7. Model Training (Random Forest) ---
message("\n--- Training Random Forest model for Differential KLF4 Binding (BJ_iPSC vs HN_SCC) ---")
message("This model will identify features distinguishing BJ_iPSC exclusive from HN_SCC exclusive peaks.")

rf_model_diff <- randomForest(Binding_Specificity ~ ., data = train_data,
                              ntree = 500,
                              importance = TRUE,
                              do.trace = 100) # Print progress every 100 trees

message("Differential binding model training complete.")
print(rf_model_diff)

# --- 8. Model Evaluation ---

# Make predictions on the test set
predictions_diff <- predict(rf_model_diff, newdata = test_data)

# Get prediction probabilities for ROC curve
# Let's use 'HN_SCC_exclusive' as the positive class for AUC calculation, for example.
probabilities_diff <- predict(rf_model_diff, newdata = test_data, type = "prob")[, "HN_SCC_exclusive"]

# Confusion Matrix
conf_matrix_diff <- confusionMatrix(predictions_diff, test_data$Binding_Specificity, positive = "HN_SCC_exclusive")
message("\n--- Confusion Matrix (BJ_iPSC vs HN_SCC Exclusive Binding) ---")
print(conf_matrix_diff)

# ROC Curve and AUC
roc_obj_diff <- roc(test_data$Binding_Specificity, probabilities_diff, levels = c("BJ_iPSC_exclusive", "HN_SCC_exclusive"))
message("\n--- ROC Curve and AUC (BJ_iPSC vs HN_SCC Exclusive Binding) ---")
print(roc_obj_diff)
plot(roc_obj_diff, main = "ROC Curve for KLF4 Differential Binding (BJ_iPSC vs HN_SCC)", col = "#1c61b6")
abline(a=0, b=1, lty=2, col="gray")

# Feature Importance
importance_df_diff <- as.data.frame(importance(rf_model_diff, type = 2)) # type=2 for MeanDecreaseGini
importance_df_diff$Feature <- rownames(importance_df_diff)
importance_df_diff <- importance_df_diff[order(-importance_df_diff$MeanDecreaseGini), ]

message("\n--- Top 15 Feature Importance (MeanDecreaseGini) for Differential Binding ---")
message("These features distinguish KLF4 peaks exclusive to BJ_iPSC vs. HN_SCC.")
print(head(importance_df_diff, 15))

# Visualize feature importance
varImpPlot(rf_model_diff, main = "Feature Importance Plot (Differential KLF4 Binding)", n.var = min(20, nrow(importance_df_diff)))

# --- 9. Interpretation and Next Steps ---
message("\n--- Interpretation Guidance ---")
message("The 'Top Feature Importance' list is key here. It shows which specific KLF4 motif properties,")
message("co-motifs, or histone modifications are most predictive of whether a KLF4 peak is exclusive")
message("to BJ_iPSC or HN_SCC.")
message("\nThis analysis directly helps you understand why KLF4 binds differently in these two cell lines.")

# --- 9.1 Explicitly Stating Top Variables for Each Cell Line ---
message("\n--- Explicit Association of Top Features with Cell Lines ---")

# Get the names of the top N important features (e.g., top 15)
top_features_for_analysis <- head(importance_df_diff$Feature, 15)

# Select only the relevant features and the Binding_Specificity column from the combined_exclusive_data
# Ensure 'combined_exclusive_data' is not empty and contains the necessary columns
if (nrow(combined_exclusive_data) > 0 && all(c(top_features_for_analysis, "Binding_Specificity") %in% names(combined_exclusive_data))) {
  
  analysis_data_for_explicit_check <- combined_exclusive_data %>%
    select(all_of(c(top_features_for_analysis, "Binding_Specificity")))
  
  # Calculate the mean/proportion for each feature by Binding_Specificity
  feature_means_by_specificity <- analysis_data_for_explicit_check %>%
    group_by(Binding_Specificity) %>%
    summarise(across(all_of(top_features_for_analysis), mean, na.rm = TRUE)) %>%
    pivot_longer(cols = -Binding_Specificity, names_to = "Feature", values_to = "Mean_Value") %>%
    pivot_wider(names_from = Binding_Specificity, values_from = Mean_Value)
  
  # Join with importance scores for better context
  feature_comparison_df <- feature_means_by_specificity %>%
    left_join(importance_df_diff %>% select(Feature, MeanDecreaseGini), by = "Feature") %>%
    arrange(desc(MeanDecreaseGini))
  
  message("\nComparison of Top Features by Cell Line Specificity:")
  message("Higher value indicates stronger association with that cell line's exclusive peaks.")
  
  # Print a nicely formatted table
  print(kable(feature_comparison_df, format = "simple", digits = 4))
  
  message("\nInterpretation:")
  for (i in 1:nrow(feature_comparison_df)) {
    feature_name <- feature_comparison_df$Feature[i]
    bj_ipsc_val <- feature_comparison_df$BJ_iPSC_exclusive[i]
    hn_scc_val <- feature_comparison_df$HN_SCC_exclusive[i]
    gini_val <- feature_comparison_df$MeanDecreaseGini[i]
    
    if (bj_ipsc_val > hn_scc_val) {
      message(sprintf("- '%s' (Importance: %.2f): Higher in BJ_iPSC_exclusive (%.4f) vs HN_SCC_exclusive (%.4f). Suggests it's more characteristic of BJ_iPSC-specific KLF4 binding.", feature_name, gini_val, bj_ipsc_val, hn_scc_val))
    } else if (hn_scc_val > bj_ipsc_val) {
      message(sprintf("- '%s' (Importance: %.2f): Higher in HN_SCC_exclusive (%.4f) vs BJ_iPSC_exclusive (%.4f). Suggests it's more characteristic of HN_SCC-specific KLF4 binding.", feature_name, gini_val, hn_scc_val, bj_ipsc_val))
    } else {
      message(sprintf("- '%s' (Importance: %.2f): Similar values in both (BJ_iPSC: %.4f, HN_SCC: %.4f). Its importance might come from interactions or thresholds.", feature_name, gini_val, bj_ipsc_val, hn_scc_val))
    }
  }
} else {
  message("Cannot perform explicit feature association: 'combined_exclusive_data' is empty or missing required columns.")
}