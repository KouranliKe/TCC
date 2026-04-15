remove(list = ls())
library(tidyverse)
library(sandwich)
library(lmtest)

source("Brazil/functions/gw_test.R")

# Add the CSFE function
f_csfe <- function(x, y_bench, y_real) {
  error_bench <- (y_bench - y_real)^2
  error_x <- (x - y_real)^2
  result <- cumsum(error_bench - error_x)
  return(result)
}

load("Brazil/data/df_new_pipeline.rda")
data <- df_final
dates <- as.Date(data$date)
nwindows <- 180
T_total <- nrow(data) # Total sample size for gw.test
y_out <- tail(data$PRECOS12_IPCA12, nwindows)
out_dates <- tail(dates, nwindows)

forecast_dir <- "Brazil/forecasts/"
files <- list.files(forecast_dir, pattern = "\\.rda$", full.names = TRUE)
all_forecasts <- list()

for (file in files) {
  model_name <- tools::file_path_sans_ext(basename(file))
  load(file)
  all_forecasts[[model_name]] <- forecasts
}

if (!"ARIMA" %in% names(all_forecasts)) {
  stop("ARIMA.rda not found. It is required as the benchmark.")
}

horizons <- c(1, 3, 6)
results_list <- list()

for (h_idx in seq_along(horizons)) {
  tau <- horizons[h_idx]
  cat("\nEvaluating Horizon:", tau, "\n")
  
  # Extract predictions for the current horizon
  h_matrix <- sapply(all_forecasts, function(x) x[, h_idx])
  
  # --- Restore Mean and Median Ensembles ---
  # Mean of all models
  h_matrix <- cbind(h_matrix, Mean_Ensemble = rowMeans(h_matrix))
  
  # Median of all models
  h_matrix <- cbind(h_matrix, Median_Ensemble = apply(h_matrix, 1, median))
  
  sarima_preds <- h_matrix[, "ARIMA"]
  sarima_rmse <- sqrt(mean((sarima_preds - y_out)^2))
  
  horizon_results <- data.frame(
    Model = colnames(h_matrix),
    Horizon = tau,
    RMSE = NA,
    Relative_RMSE = NA,
    GW_p_value = NA,
    Significance = ""
  )
  
  for (m_idx in 1:ncol(h_matrix)) {
    m_name <- colnames(h_matrix)[m_idx]
    m_preds <- h_matrix[, m_idx]
    
    # RMSE
    m_rmse <- sqrt(mean((m_preds - y_out)^2))
    horizon_results$RMSE[m_idx] <- m_rmse
    
    # Relative RMSE
    horizon_results$Relative_RMSE[m_idx] <- m_rmse / sarima_rmse
    
    # GW Test
    if (m_name == "ARIMA") {
      horizon_results$GW_p_value[m_idx] <- NA
    } else {
      tryCatch({
        gw_res <- gw.test(
          x = m_preds,      # Model 1
          y = sarima_preds, # Model 2 (Benchmark)
          p = y_out, 
          T = T_total, 
          tau = tau, 
          method = "NeweyWest", 
          alternative = "two.sided"
        )
        horizon_results$GW_p_value[m_idx] <- gw_res$p.value
        
        pval <- gw_res$p.value
        if (!is.na(pval)) {
          if (pval < 0.01) horizon_results$Significance[m_idx] <- "***"
          else if (pval < 0.05) horizon_results$Significance[m_idx] <- "**"
          else if (pval < 0.10) horizon_results$Significance[m_idx] <- "*"
        }
      }, error = function(e) {
        horizon_results$GW_p_value[m_idx] <- NA 
      })
    }
  }
  
  horizon_results <- horizon_results %>% arrange(RMSE)
  results_list[[as.character(tau)]] <- horizon_results
  
  print(horizon_results)
  
  # --- CSFE Calculation and Plotting ---
  
  # Get the top 5 models based on RMSE, excluding ARIMA and the ensembles
  top_5_models <- horizon_results %>%
    filter(!Model %in% c("ARIMA", "Mean_Ensemble", "Median_Ensemble")) %>%
    head(5) %>%
    pull(Model)
  
  csfe_df <- data.frame(Date = out_dates)
  
  # Calculate CSFE for each of the top 5 models
  for (model in top_5_models) {
    csfe_df[[model]] <- f_csfe(
      x = h_matrix[, model], 
      y_bench = sarima_preds, 
      y_real = y_out
    )
  }
  
  # Reshape data for ggplot2
  csfe_long <- csfe_df %>%
    pivot_longer(cols = -Date, names_to = "Model", values_to = "CSFE")
  
  # Plot CSFE (with thinner lines!)
  csfe_plot <- ggplot(csfe_long, aes(x = Date, y = CSFE, color = Model)) +
    geom_line(linewidth = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    labs(
      title = paste("CSFE of Top 5 Models vs ARIMA - Horizon", tau),
      subtitle = "Values > 0 indicate the model outperforms the benchmark",
      x = "Date",
      y = "Cumulative Squared Forecast Error Difference"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold")
    )
  
  print(csfe_plot)
}

# save(results_list, file = "Brazil/forecasts/evaluation_results.rda")