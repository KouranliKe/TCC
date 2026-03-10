dataprep = function(ind, df, variable, horizon, n_lags = 4, factonly = FALSE, nofact = FALSE) {
  df = df[ind,]
  y = df[,variable]
  
  if(nofact == TRUE){
    x = df
    base_names = colnames(df)
  } else {
    factors = princomp(scale(df))$scores[,1:4]
    colnames(factors) = paste0("PC", 1:4)
    
    if(factonly == TRUE){
      x = cbind(df[,variable], factors)
      base_names = c(variable, colnames(factors))
      colnames(x) = base_names
    } else {
      x = cbind(df, factors)
      base_names = c(colnames(df), colnames(factors))
      colnames(x) = base_names
    }
  }
  
  row_dates = as.Date(rownames(df))
  pandemic = ifelse(row_dates >= as.Date("2020-02-01") & row_dates <= as.Date("2022-05-31"), 1, 0)
  x = cbind(x, pandemic)
  base_names = c(base_names, 'pandemic')
  colnames(x) = base_names
  
  X = embed(as.matrix(x), n_lags)
  
  lags = rep(0:(n_lags - 1), each = length(base_names))
  suffixes = ifelse(lags == 0, "", paste0("_lag", lags))
  
  col_names = paste0(rep(base_names, times = n_lags), suffixes)
  colnames(X) = col_names
  
  Xin = X[-c((nrow(X)-horizon+1):nrow(X)),]
  
  Xout = X[nrow(X), , drop = FALSE] 
  
  yin = tail(y, nrow(Xin))
  
  return(list(Xin = Xin, Xout = Xout, yin = yin))
}

runlasso=function(ind,df,variable,horizon, n_lags = 4, alpha = 1, alpha2 = 1, adaptive = FALSE){
  library(glmnet)
  prep_data = dataprep(ind,df,variable,horizon,n_lags)
  Xin = prep_data$Xin
  yin = prep_data$yin
  Xout = prep_data$Xout
  
  modelest = ic.glmnet(Xin,yin,alpha = alpha)
  if(adaptive==TRUE){
    classo = coef(modelest)
    penalty = (abs(classo[-1])+1/sqrt(length(yin)))^(-1)
    modelest = ic.glmnet(Xin,yin, penalty.factor = penalty, alpha = alpha2)
  }
  
  forecast=predict(modelest,Xout)
  
  ### outputs ###
  coeflvl=coef(modelest)[-1]
  coefpar=coeflvl*apply(Xin,2,sd)
  lambda=modelest$lambda
  outputs=list(coeflvl=coeflvl,coefpar=coefpar,lambda=lambda)
  
  return(list(forecast=forecast, outputs=outputs))
}

runarima = function(ind, df, variable, horizon, seasonal = FALSE) {
  
  # Extract the historical target data for the current rolling window
  y_train = ts(df[ind, variable], frequency=12)
  
  # Fit the ARIMA model automatically, using BIC for model selection
  modelest = forecast::auto.arima(y_train, ic = "bic", seasonal=seasonal)
  
  # Forecast forward 'horizon' steps
  fcst = forecast::forecast(modelest, h = horizon)
  forecast_val = as.numeric(fcst$mean[horizon])
  
  # Optional outputs to track what AR/MA order was selected in each window
  outputs = list(order = forecast::arimaorder(modelest))
  
  return(list(forecast = forecast_val, outputs = outputs))
}

runrf=function(ind,df,variable,horizon,n_lags = 4){
  library(randomForest)
  prep_data = dataprep(ind,df,variable,horizon,n_lags)
  Xin = prep_data$Xin
  yin = prep_data$yin
  Xout = prep_data$Xout
  
  modelest=randomForest::randomForest(Xin,yin, importance = TRUE)
  forecast=predict(modelest,Xout)
  
  ## outputs
  importance = randomForest::importance(modelest)
  outputs = list(importance = importance)
  
  return(list(forecast=forecast, outputs = outputs))
}

runrfols=function(ind,df,variable,horizon,n_lags = 4){
  library(randomForest)
  prep_data = dataprep(ind,df,variable,horizon,n_lags)
  Xin = prep_data$Xin
  yin = prep_data$yin
  Xout = prep_data$Xout
  
  modelest=randomForest::randomForest(Xin,yin,keep.inbag = TRUE,maxnodes = 25,ntree=500)
  samples=modelest$inbag
  
  predaux=rep(0,ncol(samples))
  for(k in 1:ncol(samples)){
    saux=samples[,k]
    sboot=c()
    for(i in 1:length(saux)){
      sboot=c(sboot,rep(i,saux[i]))
    }
    xaux=Xin[sboot,]
    yaux=yin[sboot]
    tr=randomForest::getTree(modelest,k)
    selected=unique(tr[,3])
    selected=sort(selected[selected>0])
    modelols=lm(yaux~xaux[,selected])
    cols=coef(modelols)
    cols[is.na(cols)]=0
    predaux[k]=c(1,Xout[selected])%*%cols
  }
  
  forecast=mean(predaux)
  
  return(list(forecast=forecast))
}

runadalassorf=function(ind,df,variable,horizon,n_lags = 4){
  prep_data = dataprep(ind,df,variable,horizon,n_lags)
  Xin = prep_data$Xin
  yin = prep_data$yin
  Xout = prep_data$Xout
  
  lasso = HDeconometrics::ic.glmnet(Xin,yin)
  classo = coef(lasso)
  penalty = (abs(classo[-1])+1/sqrt(length(yin)))^(-1)
  adalasso = HDeconometrics::ic.glmnet(Xin,yin, penalty.factor = penalty)
  
  selected=which(adalasso$coef[-1]!=0)
  modelest=randomForest::randomForest(Xin[,selected],yin)
  forecast=predict(modelest,Xout[selected])
  
  return(list(forecast=forecast))
}

runbagging=function(ind,df,variable,horizon,n_lags = 4){
  prep_data = dataprep(ind,df,variable,horizon,n_lags)
  Xin = prep_data$Xin
  yin = prep_data$yin
  Xout = prep_data$Xout
  
  
  modelest=bagging(Xin,yin,R=100,l=5,pre.testing = "individual")
  forecast = predict(modelest,Xout)
  
  ## outputs
  
  nselect=modelest$coefficients
  nselect[nselect!=0]=1
  nselect[is.na(nselect)]=0
  nselect=colSums(nselect)
  
  outputs = list(nselect = nselect)
  
  return(list(forecast=forecast, outputs = outputs))
}

runcsr=function(ind,df,variable,horizon,n_lags = 4){
  prep_data = dataprep(ind,df,variable,horizon,n_lags)
  Xin = prep_data$Xin
  yin = prep_data$yin
  Xout = prep_data$Xout
  
  indice = which(colnames(df)==variable)
  f.seq=seq(indice,ncol(Xin)-1,ncol(df)+4)
  modelest=csr(Xin,yin,fixed.controls =c(f.seq,ncol(Xin)))
  forecast = predict(modelest,Xout)
  
  ## outputs
  
  nselect=modelest$coefficients
  nselect[nselect!=0]=1
  nselect[is.na(nselect)]=0
  nselect=colSums(nselect)
  
  outputs = list(nselect = nselect)
  
  return(list(forecast=forecast, outputs = outputs))
}

runfact=function(ind,df,variable,horizon,n_lags = 4){
  prep_data = dataprep(ind,df,variable,horizon,n_lags, factonly = TRUE)
  Xin = prep_data$Xin
  yin = prep_data$yin
  Xout = prep_data$Xout
  
  bb=Inf
  for(i in seq(5,20,5)){
    m=lm(yin~Xin[,1:i])
    crit=BIC(m)
    if(crit<bb){
      bb=crit
      modelest=m
      f.coef=coef(modelest)
    }
  }
  coef=rep(0,ncol(Xin)+1)
  coef[1:length(f.coef)]=f.coef
  coef[is.na(coef)]=0
  
  forecast=(cbind(1,Xout)%*%coef)[1]
  
  ## outputs
  outputs = list(coef = coef)
  
  return(list(forecast=forecast, outputs = outputs))
}

runtfact=function(ind,df,variable,horizon,n_lags = 4){
  
  dfaux = df[ind,]
  index = which(colnames(dfaux)==variable)
  
  mat = cbind(embed(dfaux[,variable],n_lags+1),tail(dfaux[,-index],nrow(dfaux)-n_lags))
  pretest=tfaux(mat,pre.testing="individual",fixed.controls = 1:n_lags)[-c(1:(n_lags+1))]
  
  pretest[pretest!=0]=1
  aux = rep(0,ncol(dfaux))
  aux[index] = 1
  aux[-index] = pretest
  selected=which(aux==1)
  dfreduced = df[,selected]
  
  prep_data = dataprep(ind,dfreduced,variable,horizon,n_lags, factonly = TRUE)
  Xin = prep_data$Xin
  yin = prep_data$yin
  Xout = prep_data$Xout
  
  bb=Inf
  for(i in seq(5,20,5)){
    m=lm(yin~Xin[,1:i])
    crit=BIC(m)
    if(crit<bb){
      bb=crit
      modelest=m
      f.coef=coef(modelest)
    }
  }
  coef=rep(0,ncol(Xin)+1)
  coef[1:length(f.coef)]=f.coef
  coef[is.na(coef)]=0
  
  forecast=(cbind(1,Xout)%*%coef)[1]
  
  ## outputs
  outputs = list(coef = coef)
  
  return(list(forecast=forecast, outputs = outputs))
}

runrlasso = function(ind, df, variable, horizon, n_lags = 4, post = FALSE) {
  library(hdm)
  prep_data = dataprep(ind, df, variable, horizon, n_lags)
  Xin = prep_data$Xin
  yin = prep_data$yin
  Xout = prep_data$Xout
  
  # Fit model
  modelest = hdm::rlasso(x = Xin, y = yin, post = post)
  
  # Generate forecast
  forecast = predict(modelest, newdata = unname(Xout))
  
  # Store coefficients as outputs
  outputs = list(coef = modelest$coefficients)
  
  return(list(forecast = as.numeric(forecast), outputs = outputs))
}

accumulate_model = function(forecasts){
  
  acc3 = c(rep(NA,2),sapply(1:(nrow(forecasts)-2), function(x){
    prod(1+diag(forecasts[x:(x+2),1:3]))-1
  })) 
  acc6 = c(rep(NA,5),sapply(1:(nrow(forecasts)-5), function(x){
    prod(1+diag(forecasts[x:(x+5),1:6]))-1
  }))
  acc12 = c(rep(NA,11),sapply(1:(nrow(forecasts)-11), function(x){
    prod(1+diag(forecasts[x:(x+11),1:12]))-1
  }))
  
  forecasts = cbind(forecasts,acc3,acc6,acc12)
  colnames(forecasts) = c(paste("t+",1:12,sep = ""),"acc3","acc6","acc12")
  
  return(forecasts)
  
}

ic.glmnet = function (x, y, crit = c("bic", "aic", "aicc", 
                                     "hqc"), alpha = 1, ...) 
{
  if (is.matrix(x) == FALSE) {
    x = as.matrix(x)
  }
  if (is.vector(y) == FALSE) {
    y = as.vector(y)
  }
  crit = match.arg(crit)
  n = length(y)
  model = glmnet(x = x, y = y, alpha = alpha,...)
  coef = coef(model)
  lambda = model$lambda
  df = model$df
  yhat = cbind(1, x) %*% coef
  residuals = (y - yhat)
  mse = colMeans(residuals^2)
  sse = colSums(residuals^2)
  nvar = df + 1
  bic = n * log(mse) + nvar * log(n)
  aic = n * log(mse) + 2 * nvar
  aicc = aic + (2 * nvar * (nvar + 1))/(n - nvar - 1)
  hqc = n * log(mse) + 2 * nvar * log(log(n))
  sst = (n - 1) * var(y)
  r2 = 1 - (sse/sst)
  adjr2 = (1 - (1 - r2) * (n - 1)/(nrow(x) - nvar - 1))
  crit = switch(crit, bic = bic, aic = aic, aicc = aicc, hqc = hqc)
  selected = best.model = which(crit == min(crit))
  ic = c(bic = bic[selected], aic = aic[selected], aicc = aicc[selected], 
         hqc = hqc[selected])
  result = list(coefficients = coef[, selected], ic = ic, lambda = lambda[selected], 
                nvar = nvar[selected], glmnet = model, residuals = residuals[, 
                                                                             selected], fitted.values = yhat[, selected], ic.range = crit, 
                df = df, call = match.call())
  class(result) = "ic.glmnet"
  return(result)
}

tfaux=function (mat, pre.testing = c("group-joint","joint","individual"), fixed.controls = NULL,
                t.stat = 1.96,ngroups=10)
{
  pre.testing=match.arg(pre.testing)
  y = mat[, 1]
  X = mat[, -1]
  if (pre.testing == "joint") {
    if (nrow(X) < ncol(X)) {
      stop("Error: Type = joint is only for data with more observations than variables")
    }
    m1 = lm(y ~ X)
    t1 = summary(m1)$coefficients[-1, 3]
    s1 = which(abs(t1) > t.stat)
    if (length(s1) == 0) {
      stop("Error: The pre-testing excluded all variables",
           "/n")
    }
  }
  if (pre.testing == "group-joint") {
    
    N=ncol(X)
    n=ceiling(N/ngroups)
    varind=1:N
    t1=rep(NA,N)
    for(i in 1:ngroups){
      selected=sample(order(varind),min(n,length(varind)),replace = FALSE)
      X0=X[,varind[selected]]
      m1=lm(y~X0)
      
      t0=rep(0,length(selected))
      aux=which(is.na(coef(m1)[-1]))
      t = summary(m1)$coefficients[-1, 3]
      if(length(aux)==0){
        t0=t
      }else{
        t0[-aux]=t
      }
      
      t1[varind[selected]]=t0
      varind=varind[-selected]
    }
    
    s1 = which(abs(t1) > t.stat)
    if (length(s1) == 0) {
      stop("Error: The pre-testing excluded all variables",
           "/n")
    }
  }
  if (pre.testing == "individual") {
    if (length(fixed.controls) > 0) {
      w = X[, fixed.controls]
      nonw = setdiff(1:ncol(X), fixed.controls)
    }
    else {
      w = rep(0, nrow(X))
      nonw = 1:ncol(X)
    }
    store.t = rep(NA, ncol(X))
    store.t[fixed.controls] = Inf
    for (i in nonw) {
      m1 = lm(y ~ X[, i] + w)
      t1 = summary(m1)$coefficients[2, 3]
      store.t[i] = t1
    }
    s1 = which(abs(store.t) > t.stat)
  }
  if (length(s1) > nrow(X)) {
    stop("Error: The pre-testing was not able to reduce the dimension to N<T")
  }
  m2 = lm(y ~ X[, s1])
  final.coef = rep(0, ncol(X))
  final.coef[s1] = coef(m2)[-1]
  names(final.coef) = colnames(X)
  final.coef = c(coef(m2)[1], final.coef)
  return(final.coef)
}

runlightgbm = function(ind, df, variable, horizon, n_lags = 4, extra_trees = FALSE) {
  prep_data = dataprep(ind, df, variable, horizon, n_lags)
  Xin = prep_data$Xin
  yin = prep_data$yin
  Xout = prep_data$Xout
  
  dtrain = lightgbm::lgb.Dataset(data = Xin, label = yin)
  
  params_point = list(objective = "regression", extra_trees = extra_trees)
  model_point = lightgbm::lightgbm(data = dtrain, params = params_point, verbose = -1)
  forecast = predict(model_point, Xout)
  
  # Extract chosen variables from the point model
  importance = lightgbm::lgb.importance(model_point)
  chosen_variables = importance$Feature
  
  # # 3. Lower Bound
  # params_lower = list(objective = "quantile", metric = "quantile", alpha = 0.025, extra_trees = extra_trees)
  # model_lower = lightgbm::lightgbm(data = dtrain, params = params_lower, verbose = -1)
  # lower = predict(model_lower, Xout)
  # 
  # # 4. Upper Bound 
  # params_upper = list(objective = "quantile", metric = "quantile", alpha = 0.975, extra_trees = extra_trees)
  # model_upper = lightgbm::lightgbm(data = dtrain, params = params_upper, verbose = -1)
  # upper = predict(model_upper, Xout)
  
  outputs = list(
    # lower = as.numeric(lower), 
    # upper = as.numeric(upper), 
    chosen_variables = chosen_variables
  )
  return(list(forecast = as.numeric(forecast), outputs = outputs))
}

runqlightgbm = function(ind, df, variable, horizon, n_lags = 4) {
  prep_data = dataprep(ind, df, variable, horizon, n_lags)
  Xin = prep_data$Xin
  yin = prep_data$yin
  Xout = prep_data$Xout
  
  dtrain = lightgbm::lgb.Dataset(data = Xin, label = yin)
  
  params = list(objective = "quantile", metric = "quantile", alpha = 0.5)
  model_point = lightgbm::lightgbm(data = dtrain, params = params, verbose = -1)
  forecast = predict(model_point, Xout)
  
  importance = lightgbm::lgb.importance(model_point)
  chosen_variables = importance$Feature
  
  # params_lower = list(objective = "quantile", metric = "quantile", alpha = 0.025)
  # model_lower = lightgbm::lightgbm(data = dtrain, params = params_lower, verbose = -1)
  # lower = predict(model_lower, Xout)
  # 
  # params_upper = list(objective = "quantile", metric = "quantile", alpha = 0.975)
  # model_upper = lightgbm::lightgbm(data = dtrain, params = params_upper, verbose = -1)
  # upper = predict(model_upper, Xout)
  
  outputs = list(
    # lower = as.numeric(lower), 
    # upper = as.numeric(upper), 
    chosen_variables = chosen_variables
  )
  
  return(list(forecast = as.numeric(forecast), outputs = outputs))
}

rungbm = function(ind, df, variable, horizon, n_lags = 4) {
  prep_data = dataprep(ind, df, variable, horizon, n_lags)
  train_data = as.data.frame(prep_data$Xin)
  train_data$y_target = prep_data$yin
  test_data = as.data.frame(prep_data$Xout)
  
  model_mean = gbm::gbm(y_target ~ ., data = train_data, distribution = "gaussian", n.trees = 100)
  forecast = predict(model_mean, newdata = test_data, n.trees = 100)
  
  importance = summary(model_mean, n.trees = 100, plotit = FALSE)
  chosen_variables = as.character(importance$var[importance$rel.inf > 0])
  
  # # 2. Lower Bound (alpha = 0.025)
  # model_lower = gbm::gbm(y_target ~ ., data = train_data, distribution = list(name = "quantile", alpha = 0.025), n.trees = 100)
  # lower = predict(model_lower, newdata = test_data, n.trees = 100)
  # 
  # # 3. Upper Bound (alpha = 0.975)
  # model_upper = gbm::gbm(y_target ~ ., data = train_data, distribution = list(name = "quantile", alpha = 0.975), n.trees = 100)
  # upper = predict(model_upper, newdata = test_data, n.trees = 100)
  
  outputs = list(
    # lower = as.numeric(lower),
    # upper = as.numeric(upper),
    chosen_variables = chosen_variables
  )
  
  return(list(forecast = as.numeric(forecast), outputs = outputs))
}

rungbmtuned = function(ind, df, variable, horizon, n_lags = 4) {
  prep_data = dataprep(ind, df, variable, horizon, n_lags)
  Xin = prep_data$Xin
  yin = prep_data$yin
  Xout = prep_data$Xout
  
  model_tuned = EZtune::eztune(x = Xin, y = yin, method = "gbm")
  forecast = predict(model_tuned, Xout)
  best_model = model_tuned$model
  
  importance = summary(best_model, n.trees = model_tuned$n.trees, plotit = FALSE)
  chosen_variables = as.character(importance$var[importance$rel.inf > 0])
  
  outputs = list(
    best_trees = model_tuned$n.trees,
    best_shrinkage = model_tuned$shrinkage,
    best_depth = model_tuned$n.minobsinnode,
    chosen_variables = chosen_variables
  )
  
  return(list(forecast = as.numeric(forecast), outputs = outputs))
}