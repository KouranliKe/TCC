rm(list = ls())

library(dplyr)
library(ipeadatar)
library(readxl)
library(lubridate)
library(urca)
library(tidyr)
library(tseries)

apply_transform <- function(x, t_code) {
  if (t_code == 1) return(x)                                          # Level
  if (t_code == 2) return(c(NA, diff(x)))                             # 1st Diff
  if (t_code == 3) return(c(NA, NA, diff(x, differences = 2)))        # 2nd Diff
  if (t_code == 4) return(log(x))                                     # Log Level
  if (t_code == 5) return(c(NA, diff(log(x))))                        # 1st Log-Diff
  if (t_code == 6) return(c(NA, NA, diff(log(x), differences = 2)))   # 2nd Log-Diff
  if (t_code == 7) return(c(NA, x[-1] / x[-length(x)] - 1))           # Discrete Growth
  return(x)
}

is_stationary <- function(x) {
  x_clean <- na.omit(x)
  if (length(x_clean) < 15) return(FALSE) # Safety check for sufficient data
  
  # ADF Test (Dynamic Lag Selection using AIC)
  adf_res <- suppressWarnings(ur.df(x_clean, type = "drift", selectlags = "AIC"))
  adf_pass <- adf_res@teststat[1] <= adf_res@cval[1, "5pct"]
  
  # KPSS Test 
  kpss_res <- suppressWarnings(kpss.test(x_clean, null = "Level"))
  kpss_pass <- kpss_res$p.value >= 0.05
  
  return(adf_pass && kpss_pass)
}

# Function to determine and apply the best transformation
make_stationary <- function(x, col_name) {
  
  # Assess mathematical constraints
  has_neg <- any(x < 0, na.rm = TRUE)
  has_zero <- any(x == 0, na.rm = TRUE)
  all_neg <- all(x < 0, na.rm = TRUE)
  
  # Define allowed transformation sequences based on constraints
  if (!has_neg && !has_zero) {
    # Strictly positive: Prefer Logs -> Log-Diffs -> Normal Diffs
    candidates <- c(4, 5, 6, 1, 2, 3) 
  } else if (all_neg) {
    # Strictly negative: Level -> Discrete Growth -> Diffs
    candidates <- c(1, 7, 2, 3)
  } else {
    # Contains zeros or cross-zero negatives: Cannot use logs
    candidates <- c(1, 2, 3) 
  }
  
  # Test candidates sequentially
  for (t_code in candidates) {
    x_trans <- apply_transform(x, t_code)
    if (is_stationary(x_trans)) {
      cat(sprintf("Success: '%s' is stationary. Used code: %d\n", col_name, t_code))
      return(x_trans)
    }
  }
  
  # Fail-safe: If nothing worked (e.g., hyper-seasonal data), cap at max differencing
  fallback <- tail(candidates, 1)
  cat(sprintf("WARNING: '%s' failed strict stationarity. Forced fallback code: %d\n", col_name, fallback))
  return(apply_transform(x, fallback))
}


# ==============================================================================
# MAIN PIPELINE
# ==============================================================================

# 2. Get metadata and data
dataset <- read_excel("Brazil/data/dataset.xlsx", 
                      col_types = c("text", "text", "text", "date", "date", "numeric", "text"))

metadados <- metadata(dataset$codigo)
metadados_spread <- metadata("JPM366_EMBI366")

data <- ipeadata(metadados$code)
data_spread <- ipeadata(metadados_spread$code)

# 3. Format Standard Data
df <- data %>%
  pivot_wider(names_from = "code", values_from = "value") %>%
  select(-c(uname, tcode)) %>%
  arrange(date)

# 4. Format and Aggregate Spread (Daily to Monthly)
df_spread <- data_spread %>%
  pivot_wider(names_from = "code", values_from = "value") %>%
  select(-c(uname, tcode)) %>%
  mutate(year_month = format(date, "%Y-%m")) %>%  
  group_by(year_month) %>%                        
  summarise(JPM366_EMBI366 = mean(JPM366_EMBI366, na.rm = TRUE)) %>% 
  mutate(date = as.Date(paste0(year_month, "-01"))) %>%  
  select(date, JPM366_EMBI366)

# 5. Merge and Filter Dates
df_merged <- left_join(df, df_spread, by = "date") %>%
  filter(date >= as.Date("1996-01-01") & date < as.Date("2025-01-01")) %>%
  select_if(~ !any(is.na(.))) # Drops columns with missing data in this timeframe

# 6. Apply Stationarity Transformations
cat("\n--- Starting Stationarity Tests ---\n")

# Separate dates for safekeeping
df_dates <- df_merged$date
df_features <- df_merged %>% select(-date)

# Apply the make_stationary function to every column
df_transformed <- as.data.frame(lapply(names(df_features), function(col) {
  make_stationary(df_features[[col]], col)
}))
colnames(df_transformed) <- names(df_features)

# 7. Final Cleanup
# Re-attach dates
df_final <- cbind(date = df_dates, df_transformed)

# Differencing introduces NAs at the top of the dataset. 
# We remove rows with NAs to ensure a clean block for ML.
df_final <- na.omit(df_final)

# Move the Target Variable (IPCA) to the first column (after date)
if("PRECOS12_IPCA12" %in% colnames(df_final)) {
  df_final <- df_final %>%
    relocate(PRECOS12_IPCA12, .after = date)
}

# 8. Save output
save(df_final, file = "Brazil/data/df_new_pipeline.rda")
cat("\nPipeline complete. Transformed data saved to df_new_pipeline.rda\n")