remove(list = ls())
library(tidyverse)
library(sandwich)
library(lmtest)

source("Brazil/functions/gw_test.R")

# Add the CSFE function
f_csfe <- function(x, y_bench, y_real) {
  error_bench <- (y_bench*100 - y_real*100)^2
  error_x <- (x*100 - y_real*100)^2
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
  Mean <- rowMeans(h_matrix)
  # Median of all models
  Median <- apply(h_matrix, 1, median)
  
  h_matrix <- cbind(h_matrix, Mean_Ensemble = Mean)
  h_matrix <- cbind(h_matrix, Median_Ensemble = Median)
  
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
  
  # Force RF into the plot for horizons 1 and 3 if it's not already there
  models_to_plot <- top_5_models
  if (tau %in% c(1, 3) && !"RF" %in% models_to_plot) {
    models_to_plot <- c(models_to_plot, "RF")
  }
  
  csfe_df <- data.frame(Date = out_dates)
  
  # Calculate CSFE for each of the selected models
  for (model in models_to_plot) {
    # Check if the model exists in h_matrix to avoid errors
    if (model %in% colnames(h_matrix)) {
      csfe_df[[model]] <- f_csfe(
        x = h_matrix[, model], 
        y_bench = sarima_preds, 
        y_real = y_out
      )
    } else {
      warning(paste("Model", model, "not found in h_matrix. Skipping CSFE calculation."))
    }
  }
  
  # Reshape data for ggplot2
  csfe_long <- csfe_df %>%
    pivot_longer(cols = -Date, names_to = "Model", values_to = "CSFE")
  
  # Plot CSFE
  csfe_plot <- ggplot(csfe_long, aes(x = Date, y = CSFE, color = Model)) +
    # Add shaded area for COVID-19 pandemic
    annotate("rect", 
             xmin = as.Date("2020-03-01"), 
             xmax = as.Date("2022-04-22"), 
             ymin = -Inf, 
             ymax = Inf, 
             alpha = 0.2, 
             fill = "gray50",
             color = NA) +
    geom_line(linewidth = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    labs(
      #title = paste("CSFE of Top Models vs ARIMA - Horizon", tau),
      x = "Time",
      y = "CSFE"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold")
    )
  
  print(csfe_plot)

  ggsave(
    filename = paste0("Brazil/plots/csfe_h", tau, ".svg"),
    plot = csfe_plot,
    device = "svg",
    width = 9,
    height = 3.5
  )
  
  # --- Forecast vs Actual Plotting ---
  forecast_df <- data.frame(Date = out_dates, Actual = y_out)
  
  for (model in models_to_plot) {
    if (model %in% colnames(h_matrix)) {
      forecast_df[[model]] <- h_matrix[, model]
    }
  }
  
  forecast_long_models <- forecast_df %>%
    select(-Actual) %>%
    pivot_longer(cols = -Date, names_to = "Model", values_to = "Value")
  
  forecast_plot <- ggplot() +
    # Add shaded area for COVID-19 pandemic
    annotate("rect", 
             xmin = as.Date("2020-03-01"), 
             xmax = as.Date("2022-04-22"), 
             ymin = -Inf, 
             ymax = Inf, 
             alpha = 0.2, 
             fill = "gray50",
             color = NA) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray") +
    geom_line(data = forecast_long_models, aes(x = Date, y = Value * 100, color = Model), linetype = "solid", linewidth = 0.3) +
    geom_line(data = forecast_df, aes(x = Date, y = Actual * 100, linetype = "Actual"), color = "black") +
    scale_linetype_manual(name = "", values = c("Actual" = "solid")) +
    labs(
      #title = paste("Forecast vs Actuals - Horizon", tau),
      x = "Time",
      y = "PRECOS12_IPCA12",
      color = "Model"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold")
    )
  
  print(forecast_plot)
  ggsave(
    filename = paste0("Brazil/plots/forecast_h", tau, ".svg"),
    plot = forecast_plot,
    device = "svg",
    width = 9,
    height = 3.5
  )
}

 save(results_list, file = "Brazil/results/evaluation_results.rda")