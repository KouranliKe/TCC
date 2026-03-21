# Clear environment
rm(list = ls())

# Load required libraries
library(tidyverse)

# 1. Reconstruct the target variable (y_out) to evaluate against
# Make sure these paths match your working directory
load("Brazil/data/df_new_pipeline.rda")
data <- df_final %>% select(-date) %>% as.matrix()
nwindows <- 180
y_out <- tail(data[, "PRECOS12_IPCA12"], nwindows)

# Helper function for RMSE
f_rmse <- function(x, y) {
  sqrt(mean((x - y)^2, na.rm = TRUE))
}

# 2. Get list of all forecast files in the directory
forecast_dir <- "Brazil/forecasts"
model_files <- list.files(forecast_dir, pattern = "\\.rda$", full.names = TRUE)

# 3. Extract forecasts and calculate RMSE for all models
results_list <- list()

for (file_path in model_files) {
  # Extract the model name from the file name (e.g., "Lasso" from "Lasso.rda")
  model_name <- tools::file_path_sans_ext(basename(file_path))
  
  # Load the .rda file into a temporary environment 
  # (This prevents overwriting the 'forecasts' object in the global environment in every loop)
  env <- new.env()
  load(file_path, envir = env)
  
  if (!exists("forecasts", envir = env)) {
    warning(paste("No 'forecasts' matrix found in", file_path, "- skipping."))
    next
  }
  
  fc_matrix <- env$forecasts
  
  # Calculate RMSE for each horizon.
  # Based on your pipeline: col 1 is h=1, col 2 is h=3, col 3 is h=6
  rmse_h1 <- f_rmse(fc_matrix[, 1], y_out)
  rmse_h3 <- f_rmse(fc_matrix[, 2], y_out)
  rmse_h6 <- f_rmse(fc_matrix[, 3], y_out)
  
  results_list[[model_name]] <- data.frame(
    Model = model_name,
    RMSE_h1 = rmse_h1,
    RMSE_h3 = rmse_h3,
    RMSE_h6 = rmse_h6
  )
}

# Combine all results into a single dataframe
results_df <- bind_rows(results_list)

# 4. Compare against ARIMA benchmark (Relative RMSE)
if (!"ARIMA" %in% results_df$Model) {
  warning("ARIMA.rda not found in the forecasts directory! Skipping Relative RMSE calculation.")
} else {
  # Extract ARIMA baseline values
  arima_baseline <- results_df %>% filter(Model == "SARIMA")
  
  # Calculate Relative RMSE: values < 1 mean the model outperformed ARIMA
  results_df <- results_df %>%
    mutate(
      Rel_RMSE_h1 = RMSE_h1 / arima_baseline$RMSE_h1,
      Rel_RMSE_h3 = RMSE_h3 / arima_baseline$RMSE_h3,
      Rel_RMSE_h6 = RMSE_h6 / arima_baseline$RMSE_h6
    ) %>%
    # Reorder columns for readability
    select(Model, RMSE_h1, Rel_RMSE_h1, RMSE_h3, Rel_RMSE_h3, RMSE_h6, Rel_RMSE_h6)
}

# 5. Print out the final table sorted by h=1 performance
final_results <- results_df %>% arrange(RMSE_h1)

print("Forecast Evaluation Results (Sorted by h=1 RMSE):")
print(final_results)

# Optional: Save the results table
# write.csv(final_results, "Brazil/forecasts/model_comparisons.csv", row.names = FALSE)