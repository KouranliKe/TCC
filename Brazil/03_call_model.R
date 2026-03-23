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
model_name <- "SVR"
## The function called to run models is model_function, which is a function from functions.R
model_function <- runsvr
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

# If running NNs or LLF, use sequential loop. For everything else use parallel
if (identical(model_function, runnn) || identical(model_function, runllf)) {
  for(i in for_ind){
    
    if (identical(model_function, runnn)) {
      library(h2o)
      h2o.init(nthreads = -1, max_mem_size = "16G")
      h2o.no_progress()
    }
    
    model = rolling_window(
      fn=model_function
      ,df=data
      ,nwindow=nwindows+i-1
      ,horizon=i
      ,variable="PRECOS12_IPCA12"
      ,n_lags = 12
      #,n_layers = 8 # Uncomment for NNs, change to 3, 5 or 8 for NN-3_layers NN-5_layers or NN-8_layers respectively
    )
    model_list[[i]] = model
    cat(i,"\n")
    
    if (identical(model_function, runnn)) {
    h2o.shutdown(prompt = FALSE)
    gc()
    Sys.sleep(5) # to give the virtual machine time to properly shut down
    }
  
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
      #,seasonal=TRUE # uncomment for SARIMA
      #,adaptive=TRUE # uncomment for adaLASSO or adaElasticNet
      #,alpha=0.5 # uncomment for elasticnet and adaelasticnet, set to 0 for ridge regression
      #,alpha2=0.5 # uncomment for adaelasticnet
      #,post=TRUE # uncomment for post-LASSO
      #,extra_trees = TRUE # uncomment for lightgbm-ert
      #,method='svm' # Uncomment and change for desired tuning method in runeztune
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
