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

runlasso = function(ind, df, variable, horizon, n_lags = 4, alpha = 1, alpha2 = 1, adaptive = FALSE){
  library(glmnet)
  prep_data = dataprep(ind, df, variable, horizon, n_lags)
  Xin = prep_data$Xin
  yin = prep_data$yin
  Xout = prep_data$Xout
  
  modelest = cv.glmnet(Xin, yin, alpha = alpha)
  
  if(adaptive == TRUE){
    classo = coef(modelest, s = "lambda.min")
    classo_vals = as.numeric(classo)[-1] 
    penalty = (abs(classo_vals) + 1/sqrt(length(yin)))^(-1)
    modelest = cv.glmnet(Xin, yin, penalty.factor = penalty, alpha = alpha2)
  }
  
  forecast = predict(modelest, newx = Xout, s = "lambda.min")
  
  coeflvl = as.numeric(coef(modelest, s = "lambda.min"))[-1]
  coefpar = coeflvl * apply(Xin, 2, sd)
  lambda = modelest$lambda.min
  
  outputs = list(coeflvl = coeflvl, coefpar = coefpar, lambda = lambda)
  
  return(list(forecast = as.numeric(forecast), outputs = outputs))
}

runarima = function(ind, df, variable, horizon, seasonal = FALSE) {
  y_train = ts(df[ind, variable], frequency=12)
  
  modelest = forecast::auto.arima(y_train, ic = "bic", seasonal=seasonal)
  
  fcst = forecast::forecast(modelest, h = horizon)
  forecast_val = as.numeric(fcst$mean[horizon])
  
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

runadalassorf = function(ind, df, variable, horizon, n_lags = 4){
  library(glmnet)
  library(randomForest)
  
  prep_data = dataprep(ind, df, variable, horizon, n_lags)
  Xin = prep_data$Xin
  yin = prep_data$yin
  Xout = prep_data$Xout
  
  cv_lasso = cv.glmnet(Xin, yin, alpha = 1)
  classo = coef(cv_lasso, s = "lambda.min")

  classo_vals = as.numeric(classo)[-1] 
  penalty = (abs(classo_vals) + 1/sqrt(length(yin)))^(-1)
  cv_adalasso = cv.glmnet(Xin, yin, penalty.factor = penalty, alpha = 1)
  
  cadalasso = coef(cv_adalasso, s = "lambda.min")
  cadalasso_vals = as.numeric(cadalasso)[-1]
  selected = which(cadalasso_vals != 0)
  
  # if(length(selected) == 0) {
  #   return(list(forecast = mean(yin)))
  # }
  
  modelest = randomForest::randomForest(x = Xin[, selected, drop = FALSE], y = yin)
  
  forecast = predict(modelest, newdata = Xout[, selected, drop = FALSE])
  
  return(list(forecast = as.numeric(forecast)))
}

runcsr=function(ind,df,variable,horizon,n_lags = 4){
  library(HDeconometrics)
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

# accumulate_model = function(forecasts){
#   
#   acc3 = c(rep(NA,2),sapply(1:(nrow(forecasts)-2), function(x){
#     prod(1+diag(forecasts[x:(x+2),1:3]))-1
#   })) 
#   acc6 = c(rep(NA,5),sapply(1:(nrow(forecasts)-5), function(x){
#     prod(1+diag(forecasts[x:(x+5),1:6]))-1
#   }))
#   acc12 = c(rep(NA,11),sapply(1:(nrow(forecasts)-11), function(x){
#     prod(1+diag(forecasts[x:(x+11),1:12]))-1
#   }))
#   
#   forecasts = cbind(forecasts,acc3,acc6,acc12)
#   colnames(forecasts) = c(paste("t+",1:12,sep = ""),"acc3","acc6","acc12")
#   
#   return(forecasts)
#   
# }

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
  
  importance = lightgbm::lgb.importance(model_point)
  chosen_variables = importance$Feature
  
  outputs = list(
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
  
  outputs = list(
    chosen_variables = chosen_variables
  )
  
  return(list(forecast = as.numeric(forecast), outputs = outputs))
}

rungbm = function(ind, df, variable, horizon, n_lags = 4) {
  library(gbm)
  prep_data = dataprep(ind, df, variable, horizon, n_lags)
  train_data = as.data.frame(prep_data$Xin)
  train_data$y_target = prep_data$yin
  test_data = as.data.frame(prep_data$Xout)
  
  model_mean = gbm::gbm(y_target ~ ., data = train_data, distribution = "gaussian", n.trees = 100)
  forecast = predict(model_mean, newdata = test_data, n.trees = 100)
  
  importance = summary(model_mean, n.trees = 100, plotit = FALSE)
  chosen_variables = as.character(importance$var[importance$rel.inf > 0])
  
  outputs = list(
    chosen_variables = chosen_variables
  )
  
  return(list(forecast = as.numeric(forecast), outputs = outputs))
}

runrw = function(ind, df, variable, horizon, n_lags = 4) {
  prep_data = dataprep(ind, df, variable, horizon, n_lags)
  forecast = prep_data$Xout[, variable]
  
  return(list(forecast = as.numeric(forecast)))
}

runert = function(ind, df, variable, horizon, n_lags = 4) {
  library(ranger)
  prep_data = dataprep(ind, df, variable, horizon, n_lags)
  train_df = as.data.frame(prep_data$Xin)
  train_df$y_target = prep_data$yin
  test_df = as.data.frame(prep_data$Xout)
  
  modelest = ranger::ranger(
    formula = y_target ~ ., 
    data = train_df, 
    num.trees = 10000, 
    importance = 'impurity', 
    splitrule = "extratrees", 
    replace = FALSE, 
    sample.fraction = 1
  )
  
  forecast = predict(modelest, data = test_df)$predictions
  
  outputs = list(importance = modelest$variable.importance)
  
  return(list(forecast = as.numeric(forecast), outputs = outputs))
}

runqert = function(ind, df, variable, horizon, n_lags = 4) {
  library(ranger)
  prep_data = dataprep(ind, df, variable, horizon, n_lags)
  train_df = as.data.frame(prep_data$Xin)
  train_df$y_target = prep_data$yin
  test_df = as.data.frame(prep_data$Xout)
  
  modelest = ranger::ranger(
    formula = y_target ~ ., 
    data = train_df, 
    num.trees = 10000, 
    importance = 'impurity', 
    splitrule = "extratrees", 
    replace = FALSE, 
    sample.fraction = 1,
    quantreg = TRUE
  )
  
  pred = predict(
    modelest, 
    data = test_df, 
    type = "quantiles", 
    quantiles = c(0.025, 0.5, 0.975)
  )
  
  lower_bound = pred$predictions[1, 1] # 2.5% quantile
  forecast    = pred$predictions[1, 2] # 50% quantile (Median)
  upper_bound = pred$predictions[1, 3] # 97.5% quantile
  
  outputs = list(
    importance = modelest$variable.importance,
    lower = lower_bound,
    upper = upper_bound
  )
  
  return(list(forecast = as.numeric(forecast), outputs = outputs))
}

runqrf = function(ind, df, variable, horizon, n_lags = 4) {
  library(ranger)
  prep_data = dataprep(ind, df, variable, horizon, n_lags)
  train_df = as.data.frame(prep_data$Xin)
  train_df$y_target = prep_data$yin
  test_df = as.data.frame(prep_data$Xout)
  
  modelest = ranger::ranger(
    formula = y_target ~ ., 
    data = train_df, 
    num.trees = 10000, 
    importance = 'impurity', 
    quantreg = TRUE
  )
  
  pred = predict(
    modelest, 
    data = test_df, 
    type = "quantiles", 
    quantiles = c(0.025, 0.5, 0.975)
  )
  
  lower_bound = pred$predictions[1, 1] # 2.5% quantile
  forecast    = pred$predictions[1, 2] # 50% quantile (Median point forecast)
  upper_bound = pred$predictions[1, 3] # 97.5% quantile
  
  outputs = list(
    importance = modelest$variable.importance,
    lower = lower_bound,
    upper = upper_bound
  )
  
  return(list(forecast = as.numeric(forecast), outputs = outputs))
}

runcatboost = function(ind, df, variable, horizon, n_lags = 4) {
  library(catboost)
  prep_data = dataprep(ind, df, variable, horizon, n_lags)
  Xin = prep_data$Xin
  yin = prep_data$yin
  Xout = prep_data$Xout
  
  train_pool = catboost.load_pool(data = Xin, label = yin)
  test_pool = catboost.load_pool(data = Xout)
  
  params_point = list(loss_function = 'RMSE', logging_level = 'Silent', thread_count = 4) 
  model_point = catboost.train(learn_pool = train_pool, params = params_point)
  forecast = catboost.predict(model_point, test_pool)
  
  # params_lower = list(loss_function = 'Quantile:alpha=0.025', logging_level = 'Silent')
  # model_lower = catboost.train(learn_pool = train_pool, params = params_lower)
  # lower = catboost.predict(model_lower, test_pool)
  # 
  # params_upper = list(loss_function = 'Quantile:alpha=0.975', logging_level = 'Silent')
  # model_upper = catboost.train(learn_pool = train_pool, params = params_upper)
  # upper = catboost.predict(model_upper, test_pool)
  
  importance = catboost.get_feature_importance(model_point)

  outputs = list(
    importance = importance
    # ,lower = as.numeric(lower)
    # ,upper = as.numeric(upper)
  )
  
  return(list(forecast = as.numeric(forecast), outputs = outputs))
}

runnn = function(ind, df, variable, horizon, n_lags = 4, n_layers = 3) {
  prep_data = dataprep(ind, df, variable, horizon, n_lags)
  train_df = as.data.frame(prep_data$Xin)
  train_df$y_target = prep_data$yin
  test_df = as.data.frame(prep_data$Xout)
  
  train_h2o = as.h2o(train_df)
  test_h2o = as.h2o(test_df)
  y_col = "y_target"
  x_cols = setdiff(names(train_h2o), y_col)
  
 if (n_layers == 3) {
   hidden = c(32, 16, 8)
   nfolds = 0
   epochs = 100
 } else if (n_layers == 5) {
   hidden = c(32, 16, 16, 16, 8)
   nfolds = 5
   epochs = 400
 } else if (n_layers == 8) {
   hidden = c(32, 16, 16, 16, 16, 16, 16, 8)
   nfolds = 5
   epochs = 400
 }
  
  modelest = h2o.deeplearning(
    x = x_cols,
    y = y_col,
    training_frame = train_h2o,
    activation = "Rectifier",
    hidden = hidden,
    nfolds = nfolds,
    epochs = epochs,
    train_samples_per_iteration = -2,
    seed = 1
  )
  
  pred_h2o = h2o.predict(modelest, test_h2o)
  forecast = as.numeric(as.vector(pred_h2o))
  importance = h2o.varimp(modelest)
  
  h2o.rm(train_h2o)
  h2o.rm(test_h2o)
  h2o.rm(modelest)
  
  outputs = list(
    importance = importance
  )
  
  return(list(forecast = forecast, outputs = outputs))
}

runllf = function(ind, df, variable, horizon, n_lags = 4) {
  library(grf)

  prep_data = dataprep(ind, df, variable, horizon, n_lags)
  Xin = prep_data$Xin
  yin = prep_data$yin
  Xout = prep_data$Xout

  modelest = grf::ll_regression_forest(
    X = Xin,
    Y = yin,
    num.trees = 10000,
    enable.ll.split = FALSE
  )

  pred = predict(modelest, newdata = Xout)
  forecast = pred$predictions[1]
  importance = grf::variable_importance(modelest)

  outputs = list(
    importance = importance
  )

  return(list(forecast = as.numeric(forecast), outputs = outputs))
}

runsvr = function(ind, df, variable, horizon, n_lags = 4) {
  library(e1071)

  prep_data = dataprep(ind, df, variable, horizon, n_lags)
  Xin = prep_data$Xin
  yin = prep_data$yin
  Xout = prep_data$Xout

  cost_grid = c(1, 10, 50, 100, 500, 1000, 2000)
  gamma_grid = c(1e-6, 1e-5, 1e-4, 5e-4, 1e-3, 1e-2, 0.1)
  epsilon_grid = c(0.001, 0.01, 0.05, 0.1, 0.2)

  train_size = floor(0.8 * nrow(Xin))
  val_idx = (train_size + 1):nrow(Xin)

  Xin_train = Xin[1:train_size, , drop = FALSE]
  yin_train = yin[1:train_size]

  Xin_val = Xin[val_idx, , drop = FALSE]
  yin_val = yin[val_idx]

  best_rmse = Inf
  best_params = list(cost = 1000, gamma = 1e-6, epsilon = 0.1)

  for (c in cost_grid) {
    for (g in gamma_grid) {
      for (e in epsilon_grid) {

        temp_model = e1071::svm(
          x = Xin_train,
          y = yin_train,
          type = "eps-regression",
          kernel = "radial",
          cost = c,
          gamma = g,
          epsilon = e
        )

        val_preds = predict(temp_model, newdata = Xin_val)
        rmse = sqrt(mean((yin_val - val_preds)^2))

        if (!is.na(rmse) && rmse < best_rmse) {
          best_rmse = rmse
          best_params = list(cost = c, gamma = g, epsilon = e)
        }
      }
    }
  }

  final_model = e1071::svm(
    x = Xin,
    y = yin,
    type = "eps-regression",
    kernel = "radial",
    cost = best_params$cost,
    gamma = best_params$gamma,
    epsilon = best_params$epsilon
  )

  forecast = predict(final_model, newdata = Xout)

  outputs = list(
    best_cost = best_params$cost,
    best_gamma = best_params$gamma,
    best_epsilon = best_params$epsilon,
    tot_SV = final_model$tot.SV
  )

  return(list(forecast = as.numeric(forecast), outputs = outputs))
}

runrvm = function(ind, df, variable, horizon, n_lags = 4) {
  library(kernlab)

  prep_data = dataprep(ind, df, variable, horizon, n_lags)
  Xin = as.matrix(prep_data$Xin)
  yin = as.numeric(prep_data$yin)
  Xout = as.matrix(prep_data$Xout)
  
  modelest = kernlab::rvm(
    x = Xin, 
    y = yin, 
    type = "regression"
  )
  
  forecast = predict(modelest, newdata = Xout)
  
  outputs = list(
    nRV = modelest@nRV
  )
  
  return(list(forecast = as.numeric(forecast), outputs = outputs))
}