remove(list = ls())
RNGkind("L'Ecuyer-CMRG")
set.seed(123)
### must add package for specific models ###
# library(devtools)
# install_github("gabrielrvsc/HDeconometrics")
#library(HDeconometrics)
library(tidyverse)
library(parallel)
library(RhpcBLASctl)
library(forecast)

source("Brazil/functions/rolling_window.R")
source("Brazil/functions/functions.R")

#####
## The file with the forecasts will be saved with model_name
model_name <- "ERT"
## The function called to run models is model_function, which is a function from functions.R
model_function <- runert
#####


load("Brazil/data/df_new_pipeline.rda")
data <- df_final
dates <- data$date
data <- data%>%select(-date)%>%as.matrix()
rownames(data) <- as.character(dates)

####### run rolling window ##########
nwindows <- 180
y_out <- tail(data[,"PRECOS12_IPCA12"],nwindows)
out_dates <- as.Date(tail(dates, nwindows))
for_ind <- c(1, 3, 6)
model_list <- list()

# If running NNs, use sequential loop. For everything else use parallel
if (identical(model_function, runnn3l) || identical(model_function, runnn5l)) {#|| identical(model_function, runnn8l)) {
  library(h2o)
  
  for(i in for_ind){
    h2o.init(nthreads = -1, max_mem_size = "16G") 
    h2o.no_progress()
    model = rolling_window(
      fn=model_function,
      df=data,
      nwindow=nwindows+i-1,
      horizon=i,
      variable="PRECOS12_IPCA12"
      ,n_lags = 12
    )
    model_list[[i]] = model
    cat(i,"\n")
    h2o.shutdown(prompt = FALSE)
    gc()
    Sys.sleep(2) # to give the virtual machine time to properly shut down
  }
} else {
  num_cores <- min(length(for_ind), detectCores() - 1)
  
  models_parallel <- mclapply(for_ind, function(i) {
    model = rolling_window(
      fn=model_function,
      df=data,
      nwindow=nwindows+i-1,
      horizon=i,
      variable="PRECOS12_IPCA12"
      ,n_lags = 12 # comment for (S)ARIMA
      #,seasonal=TRUE # uncomment for ARIMA
      #,adaptive=TRUE # uncomment for adaLASSO or adaElasticNet
      #,alpha=0.5 # uncomment for elastic net and adaelasticnet, set to 0 for ridge regression
      #,alpha2=0.5 # uncomment for adaelasticnet
      #,post=TRUE # uncomment for post-LASSO
      #,extra_trees = TRUE # uncomment for lightgbm-ert
    )
    return(model)
  }, mc.cores = num_cores)
  
  model_list = list()
  for(k in seq_along(for_ind)) {
    model_list[[ for_ind[k] ]] <- models_parallel[[k]]
  }
}

forecasts = Reduce(
  f = cbind,
  x = lapply(model_list, function(x)head(x$forecast,nwindows))
  ) %>% as.matrix()

outputs = lapply(model_list, function(x) head(x$outputs, nwindows))

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
      col = "blue",
      lwd = 1.5)

abline(h = 0, lty = 2, col = "darkgray")

# RMSE
f_rmse <- function(x, y) {
  sqrt(mean((x - y)^2))
}
rmse <- apply(forecasts, 2, f_rmse, y = y_out) %>% print()

# Save results
save(forecasts,
     outputs,
     file = paste("Brazil/forecasts/",model_name,".rda",sep = ""))
