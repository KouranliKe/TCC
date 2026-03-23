# Clear environment
rm(list = ls())

# Load required libraries
library(tidyverse)
library(sandwich) # Required for the HAC/NeweyWest estimators in gw.test()

# Source the Giacomini-White test function
source("Brazil/functions/gw_test.R")

# 1. Reconstruct the target variable (y_out) to evaluate against
# Make sure these paths match your working directory
load("Brazil/data/df_new_pipeline.rda")
data <- df_final %>% select(-date) %>% as.matrix()

total_T <- nrow(data) # Required 'T' parameter for the GW test (total sample size)
nwindows <- 180
y_out <- tail(data[, "PRECOS12_IPCA12"], nwindows)

# Helper function for RMSE
f_rmse <- function(x, y) {
  sqrt(mean((x - y)^2, na.rm = TRUE))
}

# 2. Get list of all forecast files in the directory
forecast_dir <- "Brazil/forecasts"
model_files <- list.files(forecast_dir, pattern = "\\.rda$", full.names = TRUE)

# 3. Extract all forecasts into a list
all_forecasts <- list()
ml_forecasts <- list()

# Define the models to exclude from the ML ensemble
benchmarks <- c("RW", "ARIMA", "SARIMA")

for (file_path in model_files) {
  model_name <- tools::file_path_sans_ext(basename(file_path))
  
  env <- new.env()
  load(file_path, envir = env)
  
  if (!exists("forecasts", envir = env)) {
    warning(paste("No 'forecasts' matrix found in", file_path, "- skipping."))
    next
  }
  
  fc_matrix <- env$forecasts
  all_forecasts[[model_name]] <- fc_matrix
  
  # Store ML models separately for the ensemble
  if (!(model_name %in% benchmarks)) {
    ml_forecasts[[model_name]] <- fc_matrix
  }
}

# 3.5 Calculate Ensembles (Mean, Median, and Weighted Mean)
if (length(ml_forecasts) > 0) {
  # Stack all ML forecast matrices into a 3D array: (time x horizons x models)
  ml_array <- array(
    unlist(ml_forecasts), 
    dim = c(nrow(ml_forecasts[[1]]), ncol(ml_forecasts[[1]]), length(ml_forecasts))
  )
  
  # 1. Simple Mean and Median
  all_forecasts[["Mean_Ensemble"]] <- apply(ml_array, c(1, 2), mean, na.rm = TRUE)
  all_forecasts[["Median_Ensemble"]] <- apply(ml_array, c(1, 2), median, na.rm = TRUE)
  
  # 2. Weighted Mean Ensemble (Inverse RMSE Weighting)
  # Calculate RMSE for each ML model and horizon. Resulting matrix is [horizons x models]
  rmse_matrix <- apply(ml_array, c(2, 3), function(preds) f_rmse(preds, y_out)) 
  
  # Calculate inverse RMSE (lower error = higher value)
  inv_rmse <- 1 / rmse_matrix
  
  # Normalize the weights so they sum to 1 for each horizon
  weights <- sweep(inv_rmse, MARGIN = 1, STATS = rowSums(inv_rmse), FUN = "/") 
  
  # Initialize an empty matrix for the weighted forecasts
  weighted_fc <- matrix(0, nrow = nrow(ml_array), ncol = ncol(ml_array))
  
  # Multiply each model's forecast by its respective weight and sum them up
  for (h in 1:ncol(weighted_fc)) {
    for (m in 1:length(ml_forecasts)) {
      weighted_fc[, h] <- weighted_fc[, h] + (ml_array[, h, m] * weights[h, m])
    }
  }
  
  all_forecasts[["Weighted_Ensemble"]] <- weighted_fc
  
} else {
  warning("No ML models found to create ensembles.")
}

# 4. Calculate RMSE, Relative RMSE, and GW Test P-values
results_list <- list()

has_benchmark <- "SARIMA" %in% names(all_forecasts)
if (has_benchmark) {
  benchmark_fc <- all_forecasts[["SARIMA"]]
} else {
  warning("SARIMA baseline not found! Skipping Relative RMSE and GW Tests.")
}

for (m_name in names(all_forecasts)) {
  fc <- all_forecasts[[m_name]]
  
  # Calculate raw RMSE
  rmse_h1 <- f_rmse(fc[, 1], y_out)
  rmse_h3 <- f_rmse(fc[, 2], y_out)
  rmse_h6 <- f_rmse(fc[, 3], y_out)
  
  # Initialize defaults for Benchmark comparisons
  rel_rmse_h1 <- NA; rel_rmse_h3 <- NA; rel_rmse_h6 <- NA
  gw_p_h1 <- NA; gw_p_h3 <- NA; gw_p_h6 <- NA
  
  if (has_benchmark) {
    # Calculate Relative RMSE
    rel_rmse_h1 <- rmse_h1 / f_rmse(benchmark_fc[, 1], y_out)
    rel_rmse_h3 <- rmse_h3 / f_rmse(benchmark_fc[, 2], y_out)
    rel_rmse_h6 <- rmse_h6 / f_rmse(benchmark_fc[, 3], y_out)
    
    # Calculate Giacomini-White Test (skip comparing benchmark to itself)
    if (m_name != "SARIMA") {
      # try() prevents the loop from crashing if a test fails (e.g., matrix singularity)
      try({
        gw_p_h1 <- gw.test(x = benchmark_fc[, 1], y = fc[, 1], p = y_out, T = total_T, tau = 1, alternative = "two.sided")$p.value
        gw_p_h3 <- gw.test(x = benchmark_fc[, 2], y = fc[, 2], p = y_out, T = total_T, tau = 3, method = "NeweyWest", alternative = "two.sided")$p.value
        gw_p_h6 <- gw.test(x = benchmark_fc[, 3], y = fc[, 3], p = y_out, T = total_T, tau = 6, method = "NeweyWest", alternative = "two.sided")$p.value
      }, silent = TRUE)
    }
  }
  
  # Append to results
  results_list[[m_name]] <- data.frame(
    Model = m_name,
    RMSE_h1 = rmse_h1, Rel_RMSE_h1 = rel_rmse_h1, GW_pval_h1 = gw_p_h1,
    RMSE_h3 = rmse_h3, Rel_RMSE_h3 = rel_rmse_h3, GW_pval_h3 = gw_p_h3,
    RMSE_h6 = rmse_h6, Rel_RMSE_h6 = rel_rmse_h6, GW_pval_h6 = gw_p_h6
  )
}

# 5. Combine and Format Final Results
results_df <- bind_rows(results_list) %>% 
  arrange(RMSE_h1)

print("Forecast Evaluation Results (Sorted by h=1 RMSE):")
print(results_df)

# Optional: Save the results table
write.csv(results_df, "Brazil/forecasts/model_comparisons.csv", row.names = FALSE)

# -------------------------------------------------------------------------
# 6. THESIS VISUALIZATION: Relative RMSE & GW Test Significance Plot
# -------------------------------------------------------------------------

# Ensure tidyr and ggplot2 are loaded (they are included in tidyverse)
library(ggplot2)
library(tidyr)

# Prepare the data for plotting
plot_data <- results_df %>%
  # Remove the benchmark from the plot (it will be represented by the x=1 line)
  filter(Model != "SARIMA") %>% 
  # Select only the columns we need
  select(Model, starts_with("Rel_RMSE"), starts_with("GW_pval")) %>%
  # Reshape data from wide to long format for ggplot
  pivot_longer(
    cols = -Model,
    names_to = c(".value", "Horizon"),
    names_pattern = "(.*)_(h[136])" # Regex to split prefix from horizon
  ) %>%
  mutate(
    # Create clean labels for the facets
    Horizon_Label = factor(Horizon, 
                           levels = c("h1", "h3", "h6"), 
                           labels = c("1-Step Ahead", "3-Steps Ahead", "6-Steps Ahead")),
    # Add significance stars based on Giacomini-White p-value
    Significance = case_when(
      GW_pval < 0.01 ~ " ***",
      GW_pval < 0.05 ~ " **",
      GW_pval < 0.10 ~ " *",
      TRUE ~ ""
    ),
    # Differentiate Ensembles from standard models for coloring
    Model_Type = ifelse(grepl("Ensemble", Model), "Ensemble", "Individual ML Model")
  )

# Build the plot
p <- ggplot(plot_data, aes(x = Rel_RMSE, y = reorder(Model, -Rel_RMSE))) +
  # Add the SARIMA benchmark line
  geom_vline(xintercept = 1, linetype = "dashed", color = "#D55E00", linewidth = 1) + 
  # Add the points for each model
  geom_point(aes(color = Model_Type), size = 3.5) +
  # Add the significance stars to the RIGHT of the points
  geom_text(aes(label = Significance), vjust = 0.75, hjust = -0.3, size = 5, color = "black") +
  # Create a separate panel for each forecasting horizon
  facet_wrap(~ Horizon_Label) +
  # Use colorblind-friendly colors
  scale_color_manual(values = c("Ensemble" = "#0072B2", "Individual ML Model" = "#56B4E9")) +
  # Expand the x-axis limits slightly so the stars don't get cut off on the right
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.15))) + 
  # Apply a clean, academic theme
  theme_bw() +
  labs(
    title = "Forecast Accuracy: Relative RMSE vs. SARIMA Benchmark",
    subtitle = "Values < 1 indicate outperformance. Stars indicate Giacomini-White test significance (* p<0.1, ** p<0.05, *** p<0.01)",
    x = "Relative RMSE (SARIMA = 1.0)",
    y = "Forecasting Model",
    color = "Model Type"
  ) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "#f0f0f0", color = "black"),
    strip.text = element_text(face = "bold", size = 11),
    plot.title = element_text(face = "bold", size = 14),
    axis.text.y = element_text(size = 10)
  )

# Display the plot in RStudio
print(p)

# Save the plot
ggsave("Brazil/plots/Thesis_Relative_RMSE_Plot.png", plot = p, width = 11, height = 7, dpi = 300)

# -------------------------------------------------------------------------
# 7. THESIS VISUALIZATION: Difference in Cumulative Squared Forecast Error (CSFE)
# -------------------------------------------------------------------------

# We need the SARIMA benchmark to compare against
if ("SARIMA" %in% names(all_forecasts)) {
  
  benchmark_fc <- all_forecasts[["SARIMA"]]
  csfe_list <- list()
  
  # Calculate CSFE difference for each model and horizon
  for (m_name in names(all_forecasts)) {
    if (m_name == "SARIMA") next # Skip benchmark
    
    fc <- all_forecasts[[m_name]]
    
    for (h_idx in 1:3) { # Assuming columns 1, 2, 3 correspond to h=1, h=3, h=6
      horizon_name <- paste0("h", c(1, 3, 6)[h_idx])
      
      # Squared Errors
      se_bench <- (benchmark_fc[, h_idx] - y_out)^2
      se_model <- (fc[, h_idx] - y_out)^2
      
      # Cumulative sum of the difference
      csfe_diff <- cumsum(se_bench - se_model)
      
      csfe_list[[paste(m_name, horizon_name, sep = "_")]] <- data.frame(
        Model = m_name,
        Horizon = horizon_name,
        Step = 1:nwindows, # Represents the pseudo-out-of-sample timeline
        CSFE_Diff = csfe_diff
      )
    }
  }
  
  # Combine into one dataframe
  csfe_df <- bind_rows(csfe_list) %>%
    mutate(
      Horizon_Label = factor(Horizon, 
                             levels = c("h1", "h3", "h6"), 
                             labels = c("1-Step Ahead", "3-Steps Ahead", "6-Steps Ahead")),
      Model_Type = ifelse(grepl("Ensemble", Model), "Ensemble", "Individual ML Model")
    )
  
  # Build the CSFE Difference Plot
  p_csfe <- ggplot(csfe_df, aes(x = Step, y = CSFE_Diff, color = Model, group = Model)) +
    # Add a bold line at 0 (the benchmark baseline)
    geom_hline(yintercept = 0, color = "black", linewidth = 1) +
    # Plot the lines
    geom_line(aes(linewidth = Model_Type, alpha = Model_Type)) +
    facet_wrap(~ Horizon_Label, scales = "free_y") +
    # Make Ensembles thicker and opaque, individual models thinner and slightly transparent
    scale_linewidth_manual(values = c("Ensemble" = 1.2, "Individual ML Model" = 0.5)) +
    scale_alpha_manual(values = c("Ensemble" = 1, "Individual ML Model" = 0.4)) +
    theme_bw() +
    labs(
      title = "Difference in Cumulative Squared Forecast Error (CSFE)",
      subtitle = "Values > 0 indicate the model is outperforming SARIMA. Upward slopes indicate period-specific outperformance.",
      x = "Out-of-Sample Step (Time)",
      y = expression(Delta * CSFE),
      color = "Model"
    ) +
    theme(
      legend.position = "right",
      strip.background = element_rect(fill = "#f0f0f0", color = "black"),
      strip.text = element_text(face = "bold", size = 11),
      plot.title = element_text(face = "bold", size = 14)
    )
  
  # Display the plot
  print(p_csfe)
  
  # Save the plot
  ggsave("Brazil/plots/Thesis_CSFE_Plot.png", plot = p_csfe, width = 12, height = 7, dpi = 300)
  
} else {
  warning("SARIMA baseline not found! Cannot calculate CSFE.")
}