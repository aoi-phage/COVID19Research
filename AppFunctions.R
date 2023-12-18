
# 平均値と中央値からweibull分布のshapeとscaleを求める
find_me_from_median <-function(xmean,xmedian){
  
  g <- function(z){
    m=z[1] # shape parameter 
    e=z[2] # scale parameter
    
    f1=e*gamma(1+(1/m))-xmean # mean=4 (default)
    f2= e*(log(2))^(1/m)-xmedian # median=3 (default)
    c(f1,f2)
  } 
  
  init <- c(1,1) #初期値
  
  gsol <- nleqslv(init, g) 
  
  return(c(gsol$x[1], gsol$x[2])) #c(m,e)
  
}# find_me_from_median



# 平均値と分散からweibull分布のshapeとscaleを求める
find_me_from_variance <-function(xmean,xvar){
  
  g <- function(z){
    m=z[1] # shape parameter 
    e=z[2] # scale parameter
    
    f1=e*gamma(1+(1/m))-xmean # mean=4 (default)
    f2=(e^2)*( gamma(1+(2/m)) - gamma(1+(1/m))^2 ) -xvar  # variance(var=13) default値に近い値
    
    c(f1,f2)
  } 
  
  init <- c(1,1) #初期値
  
  gsol <- nleqslv(init, g) 
  
  return(c(gsol$x[1], gsol$x[2])) # c(m:shape,e:scale)
  
}# find_me_from_variance



#発症者推定
#generate_onset <- function(reporteddata,m,e){
generate_onset <- function(reporteddata,xmu,xval){ #xmu:平均値、xval:分散(または中央値）
  
  #引数:reporteddata=data.frame("date","reported"), xmu:平均値, xval:分散(または中央値）
  #返り値:df_rep_ons=data.frame("date","reported","onset")
  
  # 平均値、中央値を使用する場合
  # m <- find_me_from_median(xmu,xval)[1] # shape
  # e <- find_me_from_median(xmu,xval)[2] # scale
  
  # 平均値、分散を使用する場合
   m <- find_me_from_variance(xmu,xval)[1] # shape
   e <- find_me_from_variance(xmu,xval)[2] # scale
  
  colnames(reporteddata) <- c("date","reported")
  # Shape and scale parameter of the Weibull distribution with mean 4, median 3 days.
  rep_delay = list(shape = m, scale = e) # m=1.1,e=42.(default)
  generated_onset_day_list <- list()
  
  for (d in 1:length(reporteddata$date)) {
    rep.day <-  reporteddata$date[d]
    n <-  reporteddata$reported[d]
    if (n > 0) {
      r <- rweibull(n, shape = rep_delay$shape, scale = rep_delay$scale)
      generated_onset_day_list <- c(generated_onset_day_list, rep.day-r)
    }
  }
  generated_onset_day_list <- as.Date(unlist(generated_onset_day_list))
  
  df_gen_onset <- as.data.frame(table(generated_onset_day_list))
  colnames(df_gen_onset) <- c("date", "onset")
  df_gen_onset %<>% mutate(df_gen_onset, date=as.Date(date))
  
  df_rep_ons <- merge(reporteddata, df_gen_onset, by="date", all.x=T)
  df_rep_ons [is.na(df_rep_ons)] = 0 # NA in column onset is replaced by 0
  df_rep_ons$onset <- replace(df_rep_ons$onset,nrow(df_rep_ons),df_rep_ons$onset[nrow(df_rep_ons)-1]) #最終日のonsetの値を一つ前の値で近似
  
  #return(df_rep_ons)
  return(df_rep_ons[,-2])
}

# 回帰モデルによるパラメータ推定
param_estim_by_lm <- function(xdata,ydata){
  
  #引数:xdata=dataframe(日付,説明変数)、ydata=dataframe(日付,目的変数)
  #返り値:パラメータa(切片)、パラメータb(傾き), パラメータp(残差標準偏差)
  
  xy.data <- merge(xdata, ydata, by=1) #"xdata"と"ydata"を共通日付でマージ
  xy.data <- na.omit(xy.data) #NAを消去
  
  colnames(xy.data) <- c("date","x","y")
  
  edat.log <- data.frame(lx=log10(xy.data$x), ly=log10(xy.data$y)) #対数変換
  tr.lm <- lm(edat.log$ly ~ edat.log$lx, data=edat.log)
  
  tr.a <- as.numeric(coef(tr.lm)[1]) # a:coefficient
  tr.b <- as.numeric(coef(tr.lm)[2]) # b:slope
  
  ### 信頼区間・予測区間推定用の残差標準偏差
  
  est_value <- function(x){ #予測関数
    a <- tr.a # intercept
    b <- tr.b # slope
    y <- (10^a) * (x^b)
    return(y)
  }
  
  pred.y <- unlist( lapply(xy.data$x,est_value) )
  
  log.pre.y <- log10(pred.y) 
  log.obs.y <- log10(xy.data$y)
  
  #log.pre.y
  #log.obs.y 
  
  num.d <- length(log.pre.y)
  
  #残差平方和（引数は観測発症者数）
  rss <- sum( (log.obs.y - log.pre.y)^2 )
  
  #残差分散
  #myv <-  rss/num.d
  myv <-  rss/(num.d-2)
  
  #残差標準偏差 
  rsd <- (myv)^(1/2)
  
  #return(c(tr.a,tr.b,rsd))
  return(c(tr.a,tr.b,rsd,num.d))
  
  #return(c(tr.a,tr.b))
  #return(summary(tr.lm)) #回帰分析結果のsummaryを出力
}

# 回帰モデルによる予測
epi.prediction_by_lm <- function(xdata,pa,pb,pr){
  
  #引数:xdata=dataframe(日付,予測したい期間の説明データ)、pa=切片, pb=傾き, pr=残差標準偏差
  #返り値:dataframe(日付,予測したい期間の説明データ,予測値（信頼区間・予測区間）)
  
  est_value <- function(x){ #予測関数
    a <- pa # intercept
    b <- pb # slope 
    y <- (10^a) * (x^b) 
    return(y)
  }
  
  #予測値
  prediction.y <- unlist( lapply(xdata[,2],est_value) )
  prediction.y <- log10(prediction.y) #対数変換
  
  
  #信頼区間 
  ci.up <- prediction.y + 1.96*pr*( 1/length(prediction.y) )^(1/2) #t(_0.025,無限大)=1.96　
  ci.lw <- prediction.y - 1.96*pr*( 1/length(prediction.y) )^(1/2)
  
  
  #予測区間
  pi.up <- prediction.y + 1.96*pr*( 1+( 1/length(prediction.y)) )^(1/2) #t(_0.025,無限大)=1.96
  pi.lw <- prediction.y - 1.96*pr*( 1+( 1/length(prediction.y)) )^(1/2)
  
  intl <- data.frame(ci.up, ci.lw, pi.up, pi.lw) 
  
  
  prd <- data.frame(xdata[1],log10(xdata[2]),prediction.y, ci.up, ci.lw, pi.up, pi.lw)
  
  return(prd)
  
}

# 疫学モデルによるパラメータ推定
param_estim_by_em <- function(sewagedata,onsetdata){
  
  #引数:sewagedata=dataframe(日付,下水ウイルス濃度)、onsetdata=dataframe(日付,発症者数)
  #返り値:パラメータ(v, omega, gamma, rsd)
  
  xy.data <- merge(sewagedata, onsetdata, by=1) #共通日付でmarge
  #xy.data # sewageRNAのNAを消去しない！
  colnames(xy.data) <- c("date","sewage","onset") #下水ウイルス濃度と発症者の時系列データ
  
  sval <- xy.data$sewage
  oval <- xy.data$onset
  
  T <- nrow(xy.data)
  
  # 下水ウイルス濃度から発症者を計算する関数
  # Prediction by an SeeePIAR model with vI=VA=v1, vP=v0, gammaI=gammaA=gamma
  Ot_pred <- function (TT, v, omega, gamma) { 
    
    p <- 2/3 # Fraction of symptomatic infections
    
    xlist <- c()
    t <- 1
    
    while(t<TT){
      xlist <- c( xlist, exp(-gamma*((t-1):0) ) ) # 発症後ウイルス排出量
      xlist <- c( xlist, exp(-omega*(1:(TT-t)) )) # 発症前ウイルス排出量
      
      t <- t+1
    }
    xlist <- c( xlist, exp(-gamma*((TT-1):0) ) )  # TT=tのとき
    
    shedd.mat <- matrix(xlist,TT,TT,byrow = T) #  排出量行列
    
    observed.v <- sval[1:TT]
    
    #shedd.matからobserved.vの中のNaNの位置にあたる行と列を削除する(ここで NA処理)
    if( anyNA(observed.v) ){ 
      na.val <- which(is.na(observed.v)) 
      shedd.mat <- shedd.mat[-c(na.val),-c(na.val)]
    }
    
    vna <-observed.v[!is.na(observed.v)]
    
    solve.ons <- solve(shedd.mat,(p/v)*vna) # 逆行列を解く
    
    return (solve.ons)
    
  }
  
  #発症者数実測値と下水ウイルス量から計算した発症者数の残差を求める
  eval_F_log_resid_Ot <- function(z) {
    
    p <- Ot_pred(T,z[1],z[2],z[3]) #predicted onset 
    if( anyNA(sval) ){ 
      q <- oval[-which(is.na(sval))]
    }else{q <- oval} #observed onset 
    
    p <- replace(p,which(p<=0),1) 
    plog <- log10(p)
    qlog <- log10(q)
    
    return ( sum( (plog - qlog)^2 ) ) #残差 
  }
  
  # 最適化 (minimizing with log transformed variables)
  z0 <- c(10, 1, 1) #(v, omega, gamma)
  res_log_resid_Ot <- optim(z0, eval_F_log_resid_Ot)
  
  # 最適化されたパラメータ (Optimized parameters and predictor)
  # reg_log <- list(
  #   v     = res_log_resid_Ot$par[1], 
  #   omega = res_log_resid_Ot$par[2],
  #   gamma = res_log_resid_Ot$par[3]
  # )
  reg_log <- c(res_log_resid_Ot$par[1],res_log_resid_Ot$par[2],res_log_resid_Ot$par[3])
  
  ### 信頼区間・予測区間推定用の残差標準偏差を求める
  
  #残差平方和
  rss <- eval_F_log_resid_Ot(reg_log)
  
  #残差分散
  myv <-  rss/(T-2)
  
  #残差標準偏差
  rsd <- (myv)^(1/2)
  
  reg_log2 <- list(
    v     = res_log_resid_Ot$par[1],
    omega = res_log_resid_Ot$par[2],
    gamma = res_log_resid_Ot$par[3],
    rsd = rsd
  )
  
  return(reg_log2)
  
}

# 疫学モデルによる予測（発症者数）
epi.prediction_by_em <- function(sewagedata,v,omega,gamma,pr){
  
  #引数:sewagedata=dataframe(日付,予測したい期間の下水ウイルス濃度)、パラメータ(v, omega, gamma,pr)
  #返り値:pred.data=dataframe(日付,予測したい期間の下水ウイルス濃度,予測発症者数(信頼区間・予測区間)) 
  
  colnames(sewagedata) <- c("date","sewage")
  tmax <- length(sewagedata$date)
  sval <- sewagedata$sewage
  
  Ot_pred <- function (TT, v, omega, gamma) {
    
    pI <- 2/3 
    
    xlist <- c()
    t <- 1
    
    while(t<TT){
      xlist <- c( xlist, exp(-gamma*((t-1):0) ) ) # 発症後ウイルス排出量
      xlist <- c( xlist, exp(-omega*(1:(TT-t)) )) # 発症前ウイルス排出量
      
      t <- t+1
    }
    xlist <- c( xlist, exp(-gamma*((TT-1):0) ) )  # TT=tのとき
    
    shedd.mat <- matrix(xlist,TT,TT,byrow = T) 
    
    observed.v <- sval[1:TT]
    
    if( anyNA(observed.v) ){ 
      na.val <- which(is.na(observed.v)) 
      shedd.mat <- shedd.mat[-c(na.val),-c(na.val)]
    }
    
    vna <-observed.v[!is.na(observed.v)]
    
    solve.ons <- solve(shedd.mat,(pI/v)*vna) 
    
    return (solve.ons)
    
  }#Ot_pred
  
  prediction.y <- log10( Ot_pred(tmax, v, omega, gamma) ) # 対数変換
  
  #pred.data <- cbind( na.omit(sewagedata), Ot_pred(tmax, v, omega, gamma) )
  #colnames(pred.data) <- c("date","sewage","predicted_onset") 
  
  #信頼区間 
  ci.up <- prediction.y + 1.96*pr*( 1/length(prediction.y) )^(1/2) #t(_0.025,無限大)=1.96　
  ci.lw <- prediction.y - 1.96*pr*( 1/length(prediction.y) )^(1/2)
  
  #予測区間
  pi.up <- prediction.y + 1.96*pr*( 1+( 1/length(prediction.y)) )^(1/2) #t(_0.025,無限大)=1.96
  pi.lw <- prediction.y - 1.96*pr*( 1+( 1/length(prediction.y)) )^(1/2)
  
  pred.data <- data.frame("date"=na.omit(sewagedata)$date,"sewage"=log10(na.omit(sewagedata)$sewage),prediction.y, ci.up, ci.lw, pi.up, pi.lw)
  
  
  return(pred.data)
  
}

# 統計量（決定係数）
epi.statistics <- function(xdata,odata){ # logarithmic numbers are required
  
  # xdata=data.frame(date, predicted), odata=data.frame(date, observed)
  
  edat <- merge(xdata,odata,by="date")
  
  rglog <- edat[,2] # log estimated y
  yylog <- edat[,3] # log observed y
  
  ybar<-mean(yylog)
  
  sst <- sum((yylog - ybar)^2) 　#全変動平方和
  sse <- sum( (yylog - rglog)^2)　　#残差変動平方和
  mss <- sum( (rglog - ybar)^2) # 回帰変動平方和
  
  rsq <- 1-(sse/sst)　# 決定係数 
  rsq2 <- mss/(mss+sse) # 決定係数 2 
  rsq3 <- mss/sst # 決定係数 3
  
  mse <- mean((yylog - rglog)^2) # 平均二乗誤差
  rmse <- sqrt(mse) # 二乗平均平方根誤差
  mae <- mean(abs(rglog - yylog)) # 平均絶対値誤差
  
  #stat.df <- data.frame("mean"=ybar,"SSE"=sse,"SST"=sst,"MSS"=mss,"R2.1"=rsq,"R2.2"=rsq2,"R2.3"=rsq3,"MSE"=mse,"RMSE"=rmse,"MAE"=mae)
  stat.df <- data.frame("R2.1"=rsq,"R2.2"=rsq2,"R2.3"=rsq3)
  
  return(stat.df)
}
