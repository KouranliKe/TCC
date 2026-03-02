### must add package for specific models ###
# library(devtools)
# install_github("gabrielrvsc/HDeconometrics")
library(HDeconometrics)
library(glmnet)
library(tidyverse)
library(randomForest)
library(forecast)

source("ORIGINAL_MEDEIROS/functions/rolling_window.R")
source("ORIGINAL_MEDEIROS/functions/functions.R")

#####
## The file with the forecasts will be saved with model_name
model_name = "AR(4)"
## The function called to run models is model_function, which is a function from functions.R
model_function = runar
#####


load("ORIGINAL_MEDEIROS/data/data.rda")
dates = data$date
data = data%>%select(-date)%>%as.matrix()
rownames(data) = as.character(dates)

####### run rolling window ##########
nwindows = 120
y_out = tail(data[,"CPIAUCSL"],nwindows)
model_list = list()
for_ind <- c(1, 3, 6)

for(i in for_ind){
  model = rolling_window(
    fn=model_function,
    df=data,
    nwindow=nwindows+i-1,
    horizon=i,
    variable="CPIAUCSL"
    #,n_lags = 12
    #,adaptive = TRUE
    )
  model_list[[i]] = model
  cat(i,"\n")
}

forecasts = Reduce(
  f = cbind,
  x = lapply(model_list, function(x)head(x$forecast,nwindows))
  ) %>% as.matrix()

#forecasts = accumulate_model(forecasts)

plot(tail(data[,"CPIAUCSL"],nwindows),type = "l")
lines(forecasts[,1],col = 3)
abline(h=0)
#forecasts_seasonal <- forecasts

# RMSE
f_rmse <- function(x, y) {
  sqrt(mean((x - y)^2))
}
rmse <- apply(forecasts, 2, f_rmse, y = y_out) %>% print()

# Salvar resultados
save(forecasts,file = paste("ORIGINAL_MEDEIROS/forecasts/",model_name,".rda",sep = ""))
