remove(list = ls())
set.seed(123)
### must add package for specific models ###
# library(devtools)
# install_github("gabrielrvsc/HDeconometrics")
library(HDeconometrics)
library(glmnet)
library(tidyverse)
library(parallel)
library(forecast)

source("Brazil/functions/rolling_window.R")
source("Brazil/functions/functions.R")

#####
## The file with the forecasts will be saved with model_name
model_name = "lightGBM-ERT"
## The function called to run models is model_function, which is a function from functions.R
model_function = runlightgbm
#####



load("Brazil/data/df_new_pipeline.rda")
data <- df_final
dates = data$date
data = data%>%select(-date)%>%as.matrix()
rownames(data) = as.character(dates)

####### run rolling window ##########
nwindows = 120
y_out = tail(data[,"PRECOS12_IPCA12"],nwindows)
out_dates <- as.Date(tail(dates, nwindows))
for_ind <- c(1, 3, 6)


# Single-core for loop, used for debugging purposes since multi-core debugging is harder
model_list = list()
for(i in for_ind){
  model = rolling_window(
    fn=model_function,
    df=data,
    nwindow=nwindows+i-1,
    horizon=i,
    variable="PRECOS12_IPCA12"
    ,n_lags = 12 # comment for ARIMA
    #,adaptive = TRUE # uncomment for adaLASSO
    #,seasonal=TRUE # uncomment for arima
    #,post = TRUE # uncomment for post-lasso
    ,extra_trees = TRUE # uncomment for lightgbm-ert
  )
  model_list[[i]] = model
  cat(i,"\n")
}

# Multi-core forecasting for better performance
# num_cores <- min(length(for_ind), detectCores() - 1)
# 
# models_parallel <- mclapply(for_ind, function(i) {
#   model = rolling_window(
#     fn=model_function,
#     df=data,
#     nwindow=nwindows+i-1,
#     horizon=i,
#     variable="PRECOS12_IPCA12"
#     ,n_lags = 12 # comment for arima
#     #,seasonal=TRUE # uncomment for arima
#     #,adaptive=TRUE # uncomment for adaLASSO
#     #,post=TRUE # uncomment for post-LASSO
#     #,extra_trees = TRUE # uncomment for lightgbm-ert
#   )
#   return(model)
# }, mc.cores = num_cores)
# 
# model_list = list()
# for(k in seq_along(for_ind)) {
#   model_list[[ for_ind[k] ]] <- models_parallel[[k]]
# }

forecasts = Reduce(
  f = cbind,
  x = lapply(model_list, function(x)head(x$forecast,nwindows))
  ) %>% as.matrix()

#forecasts = accumulate_model(forecasts)

# plotting forecasted vs actuals
plot(x = out_dates, 
     y = y_out, 
     type = "l", 
     col = "black",
     xlab = "Date", 
     ylab = "IPCA", 
     main = paste(model_name, "Forecast vs Actual (1-Step Ahead)"))

lines(x = out_dates, 
      y = forecasts[,1], 
      col = "blue", # Swapped to blue for better contrast
      lwd = 1.5)

abline(h = 0, lty = 2, col = "darkgray")

# RMSE
f_rmse <- function(x, y) {
  sqrt(mean((x - y)^2))
}
rmse <- apply(forecasts, 2, f_rmse, y = y_out) %>% print()

# Save results
save(forecasts,file = paste("Brazil/forecasts/",model_name,".rda",sep = ""))
