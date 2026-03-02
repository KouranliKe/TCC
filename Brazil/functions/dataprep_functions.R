### DATA PREPARATION FUNCTIONS ###
preparacao = function(X, i) {
  if (tipo[i] == 0) {
    return(log(X))
  }
  if (tipo[i] %in% c(1, 2, 3)) {
    return(X)
  }
}
  
cresc_discreto = function(X) {
  Y = c()
  for (i in 2:length(X)) {
    y = X[i]/X[i-1]-1
    Y = append(Y, y)
  }
  return(Y)
}

transform_singlestep <- function (yt, transform_id) {
  
  ytlen=length(yt)
  transformed=mat.or.vec(ytlen,1)
  transformed[]=NA
  
  if(transform_id==1) {
    transformed=yt
  }
  if(transform_id==2) {
    transformed[2:ytlen]=yt[2:ytlen]-yt[1:(ytlen-1)]
  }
  if(transform_id==3) {
    transformed[3:ytlen]=(yt[3:ytlen]-yt[2:(ytlen-1)])-(yt[2:(ytlen-1)]-yt[1:(ytlen-2)])
  }
  if(transform_id==4) {
    transformed=log(yt)
  }
  if(transform_id==5) {
    transformed[2:ytlen]=log(yt[2:ytlen])-log(yt[1:(ytlen-1)])
  }
  if(transform_id==6) {
    transformed[3:ytlen]=(log(yt[3:ytlen])-log(yt[2:(ytlen-1)]))-(log(yt[2:(ytlen-1)])-log(yt[1:(ytlen-2)]))
  }
  if(transform_id==7) {
    transformed[2:ytlen]=yt[2:ytlen]/yt[1:(ytlen-1)]-1
  }
  return(transformed)
}