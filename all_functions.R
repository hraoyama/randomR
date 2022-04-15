make_basis <- function(k, p) replace(numeric(p), k, 1)

nz_calc_detail_mean <- function(x, n_folds, n_stocks){
  res = c()
  for (i in (1:(n_folds-1))){
    res = c(res, nz_calc(x[(1+(i-1)*n_stocks):(i*n_stocks)]))
  }
  mean(res)
}

nz_calc_detail_sd <- function(x, n_folds, n_stocks){
  res = c()
  for (i in (1:(n_folds-1))){
    res = c(res, nz_calc(x[(1+(i-1)*n_stocks):(i*n_stocks)]))
  }
  sd(res)
}

nz_calc <- function(x){
  sum(x!=0)/length(x)
}


algo1 <- function(a, Cu, Cv, lambda_u, lambda_v, num_iter=200, eps=1e-10){
  
  M = a %*% t(a)
  m = nrow(M)
  n = ncol(M)
  
  u_init <- rep(1/m, m)
  v_init <- rep(1/n, n)
  iterate <- TRUE
  res = list()
  count <- 0
  idx = which.max(rowSums(M))
  
  lamd_u_scaler = 1 / as.numeric(sqrt(t(u_init) %*% Cu %*% u_init ))
  lamd_v_scaler = 1 / as.numeric(sqrt(t(v_init) %*% Cv %*% v_init ))
  
  while(iterate && (count < num_iter)) {
    
    u <- M %*% v_init
    u <- u / as.numeric(sqrt(t(u) %*% Cu %*% u ))
    u <- pmax(abs(u) - 0.5 * lambda_u * lamd_u_scaler, 0) * sign(u)
    if(sum(u^2) == 0){
      u_init = make_basis(idx, m)
      v_init = make_basis(idx, n)
      break
    }
    u <- u / as.numeric(sqrt(t(u) %*% Cu %*% u ))
    
    v <- t(M) %*% u
    v <- v / as.numeric(sqrt(t(v) %*% Cv %*% v ))
    v <- pmax(abs(v) - 0.5 * lambda_v * lamd_v_scaler, 0) * sign(v)
    if(sum(v^2) == 0){
      u_init = make_basis(idx, m)
      v_init = make_basis(idx, n)
      break
    } 
    v <- v / as.numeric(sqrt(t(v) %*% Cv %*% v ))
    
    
    Du = sum(abs(u - u_init))
    Dv = sum(abs(v - v_init))
    
    if ((Du < eps) && (Dv < eps)){
      iterate <- FALSE
    } else{
      u_init <- u
      v_init <- v
    }
    
    count <- count + 1
  }
  
  res[["v"]] <- v_init / sum(abs(v_init)) * sign(sum(as.numeric(v_init) * a))
  res[["u"]] <- u_init / sum(abs(u_init)) * sign(sum(as.numeric(u_init) * a))
  
  res
  
}


FBG_plusplus <- function(B, a, lambda, num_iter=100, eps=1e-10){
  
  n = ncol(B)
  v_init <- u_init <- rep(1/sqrt(n), n)

  iterate <- TRUE
  res = list()
  count <- 0
  
  sqrt_obj <- sqrtm(B)
  Bsqrt_inv <- sqrt_obj$Binv
  M <- Bsqrt_inv %*% a %*% t(a)
  

  lamd_scaler = 2*max(t(u_init) %*% M)
  idx = which.max(rowSums(M))
  
  while(iterate && (count < num_iter)) {
    
    v = pmax(abs(t(u_init)%*% M ) - 0.5*lambda*lamd_scaler, 0)*sign(t(u_init) %*% M)
    v = t(v)
    if(t(v) %*% v == 0){ 
      v_init = make_basis(idx, n)
      break 
    }
    
    K = M %*% v
    u = K / as.numeric(sqrt(t(K) %*% K))
    
    Du = sum(abs(u - u_init))
    Dv = sum(abs(v - v_init))
    
    if ((Du < eps) && (Dv < eps)){
      iterate <- FALSE
    } else{
      u_init <- u
      v_init <- v
    }
    
    count <- count + 1
  }
  
  res[["v"]] <- as.numeric(v_init) / sum(abs(v_init)) * sign(sum(as.numeric(v_init) * a))
  res[["u"]] <- u_init
  res
  
}


sparsifier <- function(B, a, lambda){
  
  B_inv = solve(B)
  n = nrow(B_inv)
  v_init <- B_inv %*% a
  res=list()
  idx = which.max(rowSums(a %*% t(a)))
  #lamd_scaler = 1 / as.numeric(sqrt(t(v_init) %*% B_inv %*% v_init ))
    
  v <- v_init / sum(abs(v_init))
  v <- pmax(abs(v) - lambda, 0) * sign(v)
  if(sum(v^2) == 0){
      v = make_basis(idx, n)
  }

  res[["v"]] = v / sum(abs(v)) * sign(sum(as.numeric(v) * as.numeric(a)))
  res
}



modified_PMD <- function(a, Cu, Cv, lambda_u, lambda_v, num_iter=100, eps=1e-10){
  
  M = a %*% t(a)
  m = nrow(M)
  n = ncol(M)
  
  u_init <- rep(1/m, m)
  v_init <- rep(1/n, n)
  iterate <- TRUE
  res = list()
  count <- 0
  idx = which.max(rowSums(M))
  
  while(iterate && (count < num_iter)) {
    
    u <- M %*% v_init
    u <- pmax(abs(u) - 0.5 * lambda_u*sum(u), 0) * sign(u)
    if(sum(u*u) == 0){
      u_init = make_basis(idx, m)
      v_init = make_basis(idx, n)
      break
    }
    u <- u / as.numeric(sqrt(t(u) %*% Cu %*% u ))
    
    v <- t(M) %*% u
    v <- pmax(abs(v) - 0.5 * lambda_v*sum(v), 0) * sign(v)
    if(sum(v*v) == 0){
      u_init = make_basis(idx, m)
      v_init = make_basis(idx, n)
      break
    }
    v <- v / as.numeric(sqrt(t(v) %*% Cv %*% v ))
    
    Du = sum(abs(u - u_init))
    Dv = sum(abs(v - v_init))
    
    if ((Du < eps) && (Dv < eps)){
      iterate <- FALSE
    } else{
      u_init <- u
      v_init <- v
    }
    
    count <- count + 1
    
  }
  
  res[["v"]] <- v_init / sum(abs(v_init)) * sign(sum(as.numeric(v_init) * as.numeric(a)))
  res[["u"]] <- u_init / sum(abs(u_init)) * sign(sum(as.numeric(u_init) * as.numeric(a)))
  res
}


ts_crossValidation <- function(tbl, n_folds){
  
  m = nrow(tbl)
  n = m %/% n_folds
  
  folds = list()
  for (i in 1:n_folds){
    
    if (i != n_folds) {
      folds[[i]] = tbl[(1+(i-1)*n):(i*n),]
    } else{
      folds[[i]] = tbl[(1+(i-1)*n):m,]
    }
  }
  folds["length_fold"] = n
  folds
}


generic_portfolio_constructor <- function(type, bigSigma_local, 
                                          smallSigma_local, mu_local){
  
  if (type == 'MD'){
    a = smallSigma_local
    B = bigSigma_local
  } else if  (type == 'MC'){
    a = smallSigma_local
    B = diag(smallSigma_local*smallSigma_local)
  } else if (type == 'SR') {
    a = mu_local
    B = bigSigma_local
  } else {
    n = nrow(bigSigma_local)
    return(rep(1/n, n))
  } 
  
  return(Rayleigh_closed_form(a, B))
}


Rayleigh_closed_form <- function(a, B){
  
  B_inv = solve(B)
  num <- B_inv %*% a 
  den <- (t(a) %*% B_inv %*% a)
  w <- num / as.numeric(den)
  w / sum(abs(w)) * sign(sum(as.numeric(w) * as.numeric(a)))
}


ts_evaluate <- function(tbl, n_folds, pmd_constr, algo1_constr, fbg_constr, sparse_constr,
                        memoize=FALSE, use_train=FALSE){
  
  folds = ts_crossValidation(tbl, n_folds)
  pmd_MD_returns = c()
  pmd_MD_weights = c()
  pmd_MC_returns = c()
  pmd_MC_weights = c()
  pmd_SR_returns = c()
  pmd_SR_weights = c()
  
  algo1_MD_returns = c()
  algo1_MD_weights = c()
  algo1_MC_returns = c()
  algo1_MC_weights = c()
  algo1_SR_returns = c()
  algo1_SR_weights = c()
  
  modpmd_MD_returns = c()
  modpmd_MD_weights = c()
  modpmd_MC_returns = c()
  modpmd_MC_weights = c()
  modpmd_SR_returns = c()
  modpmd_SR_weights = c()
  
  fbg_MD_returns = c()
  fbg_MD_weights = c()
  fbg_SR_returns = c()
  fbg_SR_weights = c()
  fbg_MC_returns = c()
  fbg_MC_weights = c()
  
  sparse_MD_returns = c()
  sparse_MD_weights = c()
  sparse_SR_returns = c()
  sparse_SR_weights = c()
  sparse_MC_returns = c()
  sparse_MC_weights = c()
  
  gen_md_returns = c()
  gen_md_weights = c()
  gen_mc_returns = c()
  gen_mc_weights = c()
  gen_sr_returns = c()
  gen_sr_weights = c()
  gen_unif_returns = c()
  gen_unif_weights = c()
  
  tr_tbl = c()
  
  if(use_train) {
    starting_fold = 1
  } else{
    starting_fold = 2
  }
  
  for (i in starting_fold:n_folds){
    
    val_tbl = folds[[i]]
    if (memoize){
      tr_tbl = rbind(tr_tbl, folds[[i-starting_fold+1]])
    } else{
      tr_tbl = folds[[i-starting_fold+1]]
    }
    
    avg_ret = apply(tr_tbl, 2, mean)
    sigma_tr = apply(tr_tbl, 2, sd)
    
    sigma2_mtrx = sigma_tr %*% t(sigma_tr)
    return2_mtrx = avg_ret %*% t(avg_ret)
    cov_mtrx = cov(tr_tbl, use="pairwise.complete.obs")
    diag_sigma2_mtrx = diag(sigma_tr*sigma_tr)
    
    
    if(FALSE){
    ws_MD_pmd = foreach(sumabs = pmd_constr,.combine = cbind) %do% {
      pmd_obj <- PMD(sigma2_mtrx, sumabsv = sumabs, K=1, center = F, trace=F)
      w <- pmd_obj$v
      #normalization
      w <- w / sum(abs(w)) * sign(sum(as.numeric(w) * as.numeric(sigma_tr)))
      return(w)
    }
    
    ws_MD_algo1 <- foreach(sumabs = algo1_constr,.combine = cbind) %do% {
      algo1_obj <- algo1(sigma_tr, cov_mtrx, cov_mtrx, lambda_u=sumabs, lambda_v=sumabs)
      w <- algo1_obj$v
      return(w)
    }
    
    ws_MD_modpmd <- foreach(sumabs = algo1_constr,.combine = cbind) %do% {
      modpmd_obj <- modified_PMD(sigma_tr, cov_mtrx, cov_mtrx, 
                                 lambda_u=sumabs, lambda_v=sumabs)
      w <- modpmd_obj$v
      return(w)
    }
    
    ws_MD_fbg <- foreach(sumabs = fbg_constr,.combine = cbind) %do% {
      fbg_obj <- FBG_plusplus(cov_mtrx, sigma_tr, lambda=sumabs)
      w <- fbg_obj$v
      
      return(w)
    }
    }
    ws_MD_sparse <- foreach(sumabs = sparse_constr,.combine = cbind) %do% {
      sparse_obj <- sparsifier(cov_mtrx, sigma_tr, lambda=sumabs)
      w <- sparse_obj$v
      
      return(w)
    }
    if(FALSE){
    ws_SR_pmd = foreach(sumabs = pmd_constr,.combine = cbind) %do% {
      pmd_obj <- PMD(return2_mtrx, sumabsv = sumabs, K=1, center = F, trace=F)
      w <- pmd_obj$v
      #normalization
      w <- w / sum(abs(w)) * sign(sum(as.numeric(w) * as.numeric(avg_ret)))
      return(w)
    }
    
    ws_SR_algo1 <- foreach(sumabs = algo1_constr,.combine = cbind) %do% {
      algo1_obj <- algo1(avg_ret, cov_mtrx, cov_mtrx, lambda_u=sumabs, lambda_v=sumabs)
      w <- algo1_obj$v
      return(w)
    }
    
    ws_SR_modpmd <- foreach(sumabs = algo1_constr,.combine = cbind) %do% {
      modpmd_obj <- modified_PMD(avg_ret, cov_mtrx, cov_mtrx,
                                 lambda_u=sumabs, lambda_v=sumabs)
      w <- modpmd_obj$v
      return(w)
    }
    
    ws_SR_fbg <- foreach(sumabs = fbg_constr,.combine = cbind) %do% {
      fbg_obj <- FBG_plusplus(cov_mtrx, avg_ret, lambda=sumabs)
      w <- fbg_obj$v
      
      return(w)
    }
    }
    ws_SR_sparse <- foreach(sumabs = sparse_constr,.combine = cbind) %do% {
      sparse_obj <- sparsifier(cov_mtrx, avg_ret, lambda=sumabs)
      w <- sparse_obj$v
      
      return(w)
    }
    if(FALSE){
    ws_MC_algo1 <- foreach(sumabs = algo1_constr,.combine = cbind) %do% {
      algo1_obj <- algo1(sigma_tr, diag_sigma2_mtrx, diag_sigma2_mtrx,
                         lambda_u=sumabs, lambda_v=sumabs)
      w <- algo1_obj$v
      return(w)
    }
    
    ws_MC_modpmd <- foreach(sumabs = algo1_constr,.combine = cbind) %do% {
      modpmd_obj <- modified_PMD(sigma_tr, diag_sigma2_mtrx, diag_sigma2_mtrx,
                                 lambda_u=sumabs, lambda_v=sumabs)
      w <- modpmd_obj$v
      return(w)
    }
    
    ws_MC_fbg <- foreach(sumabs = fbg_constr,.combine = cbind) %do% {
      fbg_obj <- FBG_plusplus(diag_sigma2_mtrx, sigma_tr, lambda=sumabs)
      w <- fbg_obj$v
      
      return(w)
    }
    }
    ws_MC_sparse <- foreach(sumabs = sparse_constr,.combine = cbind) %do% {
      sparse_obj <- sparsifier(diag_sigma2_mtrx, sigma_tr, lambda=sumabs)
      w <- sparse_obj$v
      
      return(w)
    }
    

    if(use_train){
      val_tbl <- tr_tbl
    }
    
    # additional portfolios
    ws_md = generic_portfolio_constructor("MD", cov_mtrx, sigma_tr, avg_ret)
    ws_mc = generic_portfolio_constructor("MC", cov_mtrx, sigma_tr, avg_ret)
    ws_sr = generic_portfolio_constructor("SR", cov_mtrx, sigma_tr, avg_ret)
    ws_unif = generic_portfolio_constructor("UN", cov_mtrx, sigma_tr, avg_ret)
    
    temp <- val_tbl %*% ws_md
    gen_md_returns <- rbind(gen_md_returns, temp)
    gen_md_weights <- rbind(gen_md_weights, ws_md)
    
    temp <- val_tbl %*% ws_mc
    gen_mc_returns <- rbind(gen_mc_returns, temp)
    gen_mc_weights <- rbind(gen_mc_weights, ws_mc)
    
    temp <- val_tbl %*% ws_sr
    gen_sr_returns <- rbind(gen_sr_returns, temp)
    gen_sr_weights <- rbind(gen_sr_weights, ws_sr)
    
    temp <- val_tbl %*% ws_unif
    gen_unif_returns <- rbind(gen_unif_returns, temp)
    gen_unif_weights <- rbind(gen_unif_weights, ws_unif)
    
    if(FALSE){
    temp <- val_tbl %*% ws_MD_pmd
    pmd_MD_returns <- rbind(pmd_MD_returns, temp)
    pmd_MD_weights <- rbind(pmd_MD_weights, ws_MD_pmd)
    
    temp <- val_tbl %*% ws_MD_algo1
    algo1_MD_returns <- rbind(algo1_MD_returns, temp)
    algo1_MD_weights <- rbind(algo1_MD_weights, ws_MD_algo1)
    
    temp <- val_tbl %*% ws_MD_modpmd
    modpmd_MD_returns <- rbind(modpmd_MD_returns, temp)
    modpmd_MD_weights <- rbind(modpmd_MD_weights, ws_MD_modpmd)
    
    temp <- val_tbl %*% ws_MD_fbg
    fbg_MD_returns <- rbind(fbg_MD_returns, temp)
    fbg_MD_weights <- rbind(fbg_MD_weights, ws_MD_fbg)
    }
    temp <- val_tbl %*% ws_MD_sparse
    sparse_MD_returns <- rbind(sparse_MD_returns, temp)
    sparse_MD_weights <- rbind(sparse_MD_weights, ws_MD_sparse)
    
    if(FALSE){
    temp <- val_tbl %*% ws_SR_pmd
    pmd_SR_returns <- rbind(pmd_SR_returns, temp)
    pmd_SR_weights <- rbind(pmd_SR_weights, ws_SR_pmd)
    
    temp <- val_tbl %*% ws_SR_algo1
    algo1_SR_returns <- rbind(algo1_SR_returns, temp)
    algo1_SR_weights <- rbind(algo1_SR_weights, ws_SR_algo1)
    
    temp <- val_tbl %*% ws_SR_modpmd
    modpmd_SR_returns <- rbind(modpmd_SR_returns, temp)
    modpmd_SR_weights <- rbind(modpmd_SR_weights, ws_SR_modpmd)
    
    temp <- val_tbl %*% ws_SR_fbg
    fbg_SR_returns <- rbind(fbg_SR_returns, temp)
    fbg_SR_weights <- rbind(fbg_SR_weights, ws_SR_fbg)
    }
    temp <- val_tbl %*% ws_SR_sparse
    sparse_SR_returns <- rbind(sparse_SR_returns, temp)
    sparse_SR_weights <- rbind(sparse_SR_weights, ws_SR_sparse)
    if(FALSE){
    temp <- val_tbl %*% ws_MC_algo1
    algo1_MC_returns <- rbind(algo1_MC_returns, temp)
    algo1_MC_weights <- rbind(algo1_MC_weights, ws_MC_algo1)
    
    temp <- val_tbl %*% ws_MC_modpmd
    modpmd_MC_returns <- rbind(modpmd_MC_returns, temp)
    modpmd_MC_weights <- rbind(modpmd_MC_weights, ws_MC_modpmd)
    
    temp <- val_tbl %*% ws_MC_fbg
    fbg_MC_returns <- rbind(fbg_MC_returns, temp)
    fbg_MC_weights <- rbind(fbg_MC_weights, ws_MC_fbg)
    }
    temp <- val_tbl %*% ws_MC_sparse
    sparse_MC_returns <- rbind(sparse_MC_returns, temp)
    sparse_MC_weights <- rbind(sparse_MC_weights, ws_MC_sparse)
  }
  
  
  list(#"pmd_MD"=pmd_MD_returns, "algo1_MD"=algo1_MD_returns,
       #"modpmd_MD"=modpmd_MD_returns, "fbg_MD"=fbg_MD_returns,
       "sparse_MD"=sparse_MD_returns,
       
       #"ws_MD_pmd"=pmd_MD_weights, "ws_MD_algo1"=algo1_MD_weights,
       #"ws_MD_modpmd"=modpmd_MD_weights, "ws_MD_fbg"=fbg_MD_weights,
       "ws_MD_sparse"=sparse_MD_weights,
       
       #"pmd_SR"=pmd_SR_returns, "algo1_SR"=algo1_SR_returns,
       #"modpmd_SR"=modpmd_SR_returns, "fbg_SR"=fbg_SR_returns,
       "sparse_SR"=sparse_SR_returns,
       
       #"ws_SR_pmd"=pmd_SR_weights, "ws_SR_algo1" = algo1_SR_weights,
       #"ws_SR_modpmd"=modpmd_SR_weights, "ws_SR_fbg"=fbg_SR_weights,
       "ws_SR_sparse"=sparse_SR_weights,
       
       #"pmd_MC"=pmd_MD_returns, "algo1_MC"=algo1_MC_returns,
       #"modpmd_MC"=modpmd_MC_returns, "fbg_MC"=fbg_MC_returns,
       "sparse_MC"=sparse_MC_returns,
       
       #"ws_MC_pmd"=pmd_MD_weights, "ws_MC_algo1" = algo1_MC_weights,
       #"ws_MC_modpmd"=modpmd_MC_weights, "ws_MC_fbg"=fbg_MC_weights,
       "ws_MC_sparse"=sparse_MC_weights,
      
       "gen_MD"=gen_md_returns, "ws_gen_MD"=gen_md_weights,
       "gen_MC"=gen_mc_returns, "ws_gen_MC"=gen_mc_weights,
       "gen_SR"=gen_sr_returns, "ws_gen_SR"=gen_sr_weights,
       "gen_UNIF"=gen_unif_returns, "ws_gen_UNIF"=gen_unif_weights,
       "fold_length"=folds$length_fold)
  
}




data_transformer <- function(data_daily, is_frac){
  
  daily_data = daily_data[, c("Index",colnames(daily_data)[sapply(colnames(daily_data), function(x) { str_detect(x,"Close") } )]),with=F]
  daily_data[,date:=ymd_hms(Index)]
  daily_data[,Index:=NULL]
  if("LHN.VX.Close" %in% colnames(daily_data)){ daily_data[,LHN.VX.Close:=NULL] }
  if("EL.PA.Close" %in% colnames(daily_data)){ daily_data[,EL.PA.Close:=NULL] }
  if("X1COV.DE.Close" %in% colnames(daily_data)){ daily_data[,X1COV.DE.Close:=NULL] }
  if("VNA.DE.Close" %in% colnames(daily_data)){ daily_data[,VNA.DE.Close:=NULL] }
  if("AUTO.L.Close" %in% colnames(daily_data)){ daily_data[,AUTO.L.Close:=NULL] }
  if("BHP.L.Close" %in% colnames(daily_data)){ daily_data[,BHP.L.Close:=NULL] }
  if("DLG.L.Close" %in% colnames(daily_data)){ daily_data[,DLG.L.Close:=NULL] }
  if("CCH.L.Close" %in% colnames(daily_data)){ daily_data[,CCH.L.Close:=NULL] }
  if("NMC.L.Close" %in% colnames(daily_data)){ daily_data[,NMC.L.Close:=NULL] }
  if("OCDO.L.Close" %in% colnames(daily_data)){ daily_data[,OCDO.L.Close:=NULL] }
  if("GLEN.L.Close" %in% colnames(daily_data)){ daily_data[,GLEN.L.Close:=NULL] }
  if("EVR.L.Close" %in% colnames(daily_data)){ daily_data[,EVR.L.Close:=NULL] }
  if("TUI.L.Close" %in% colnames(daily_data)){ daily_data[,TUI.L.Close:=NULL] }
  
  daily_data_carried_forward = na.locf(daily_data)
  daily_data_carried_forward = daily_data_carried_forward[,c("date",colnames(daily_data_carried_forward)[colnames(daily_data_carried_forward) != "date"]),with=F]
  
  alldata = as.matrix(daily_data_carried_forward[,2:ncol(daily_data_carried_forward),with=F],
                      ncol=(ncol(daily_data_carried_forward)-1),
                      nrow=nrow(daily_data_carried_forward))
  alldata = apply(alldata,2,as.numeric)
  
  logWorkDataAll = matrix(NA, nrow=(nrow(alldata)-1), ncol=ncol(alldata))  # will house the log returns
  colnames(logWorkDataAll) = colnames(alldata)
  for(colc in 1:ncol(alldata)) {
    logWorkDataAll[,colc] = diff(log(alldata[,colc]))
  }
  
  allDataWithZeros = logWorkDataAll
  allDataWithZeros[is.na(allDataWithZeros)] = 0.0
  
  # IS = in-sample  ; OOS = out-of-sample
  allDataWithZeros_OOS = allDataWithZeros[round(is_frac*nrow(allDataWithZeros)):nrow(allDataWithZeros),]
  allDataWithZeros_IS = allDataWithZeros[1:(round(is_frac*nrow(allDataWithZeros))-1),]
  
  list("IS"=allDataWithZeros_IS, "OOS"=allDataWithZeros_OOS, "ALL"=allDataWithZeros)
  
}


results_by_fold <- function(data, n_folds, func){
  
  N = dim(data)
  p = N %/% n_folds
  r = N %% n_folds
  
  par(mfrow=c(3,3))
  for (i in (1:n_folds)){
    if (i < n_folds){

      plot(sumabsvs, apply(data[(1+(i-1)*p):(i*p),], 2, func), cex=.3, xlab='L1 norm', ylab='avg log return', main=paste("fold", i))
    }
    else{
      plot(sumabsvs, apply(data[((1+(i-1)*p):N),], 2, func), cex=.3, xlab='L1 norm', ylab='avg log return', main=paste("fold", i))
    }
  }
  par(mfrow=c(1,1))
}


sharpe_by_fold <- function(data, fold_length, func){
  
  N = nrow(data)
  n_folds = N %/% fold_length
  r = N %% fold_length

  folds = list()
  for (i in 1:n_folds){
    folds[[i]] = data[((1+(i-1)*fold_length):(i*fold_length)),]
  }
  if (r > 0){
    folds[[i+1]] = data[((1+i*fold_length):N),]
    n_folds = n_folds + 1
  }
  
  all_folds = c()
  for (i in 1:n_folds){
    all_folds = rbind(all_folds, apply(folds[[i]], 2, func))
  }
  
  all_folds
}

optimal_constraints_searcher <- function(data, insample_fraction, n_folds, use_train=FALSE){
  
  data_splits = data_transformer(data, insample_fraction)
  allDataWithZeros_IS = data_splits$IS
  allDataWithZeros_OOS = data_splits$OOS
  n_stocks = ncol(allDataWithZeros_OOS)
  
  
  sumabsvs = linspace(1.0, sqrt(ncol(allDataWithZeros_IS)), n=100)
  algo1_lamdas = c(linspace(0.0, 0.1, n=100), linspace(0.101, 0.15, n=50))
  fbg_lamdas = linspace(0.0, 0.999, n=100)
  sparse_lamdas = linspace(0.0, 0.15, n=100)
  
  lamda_grid = sparse_lamdas
  method = "sparse"
  
  resulting_returns <- ts_evaluate(allDataWithZeros_IS, n_folds, sumabsvs, algo1_lamdas, fbg_lamdas, sparse_lamdas, use_train=use_train)
  
  plot(lamda_grid, sqrt(252)*apply(resulting_returns[[paste(method, "MD", sep="_")]], 2, mean)/apply(resulting_returns[[paste(method, "MD", sep="_")]], 2, sd), xaxp=c(0., .15, 15), cex=.3, xlab='lambda', ylab='Sharpe ratio')
  grid (NULL,NULL, lty = 6, col = "cornsilk2")
  axis(3, lamda_grid, round(apply(resulting_returns[[paste("ws_MD", method, sep="_")]], 2, nz_calc), 2))
  title("Maximum Diversification", line = 2)

  plot(lamda_grid, sqrt(252)*apply(resulting_returns[[paste(method, "SR", sep="_")]], 2, mean)/apply(resulting_returns[[paste(method, "SR", sep="_")]], 2, sd), xaxp=c(0., .15, 15), cex=.3, xlab='lambda', ylab='Sharpe ratio')
  grid (NULL,NULL, lty = 6, col = "cornsilk2")
  axis(3, lamda_grid, round(apply(resulting_returns[[paste("ws_SR", method, sep="_")]], 2, nz_calc), 2))
  title("Sharpe Ratio", line = 2)
  
  plot(lamda_grid, sqrt(252)*apply(resulting_returns[[paste(method, "MC", sep="_")]], 2, mean)/apply(resulting_returns[[paste(method, "MC", sep="_")]], 2, sd),xaxp=c(0., .15, 15),  cex=.3, xlab='lambda', ylab='Sharpe ratio')
  grid (NULL,NULL, lty = 6, col = "cornsilk2")
  axis(3, lamda_grid, round(apply(resulting_returns[[paste("ws_MC", method, sep="_")]], 2, nz_calc), 2))
  title("Minimum Concentration", line = 2)
  
  return(resulting_returns)
}


results_saver <- function(data, insample_fraction, n_folds, constraints, use_train=FALSE){
  
  data_splits = data_transformer(data, insample_fraction)
  allDataWithZeros_IS = data_splits$IS
  allDataWithZeros_OOS = data_splits$OOS
  n_stocks = ncol(allDataWithZeros_OOS)
  
  test_pred_returns <- ts_predict(allDataWithZeros_IS, allDataWithZeros_OOS, n_folds,
                                  constraints)
  method = "sparse"
  
  mean_pred_returns <- c(mean(test_pred_returns$unif), mean(test_pred_returns[[paste(method, "SR", sep="_")]]), mean(test_pred_returns[[paste(method, "MD", sep="_")]]), mean(test_pred_returns[[paste(method, "MC", sep="_")]]))
  std_pred_returns <- c(sd(test_pred_returns$unif), sd(test_pred_returns[[paste(method, "SR", sep="_")]]), sd(test_pred_returns[[paste(method, "MD", sep="_")]]), sd(test_pred_returns[[paste(method, "MC", sep="_")]]))
  sharpe_ratios <- mean_pred_returns/std_pred_returns * sqrt(252)
  
  resulting_df <- rbind(mean_pred_returns, std_pred_returns, sharpe_ratios)
  resulting_df = as.data.frame(resulting_df)
  names(resulting_df) <- c("UNIF", "SR", "MD", "MC")
  rownames(resulting_df) <- c("MR", "STD", "SR")
  resulting_df <- round(resulting_df, 5)
  resulting_df
}


results_in_sample_saver <- function(results_obj, constraints){
  
  lamdas = linspace(0.0, 0.15, n=100)
  
  method = "sparse"
  md_lam = constraints[[paste("md", method, sep="_")]] 
  sr_lam = constraints[[paste("sr", method, sep="_")]] 
  mc_lam = constraints[[paste("mc", method, sep="_")]] 
  
  idx_md = which.min(abs(lamdas - md_lam))
  idx_sr = which.min(abs(lamdas - sr_lam))
  idx_mc = which.min(abs(lamdas - mc_lam))
  
  
  mean_pred_returns_star <- c(mean(results_obj$gen_UNIF), mean(results_obj[[paste(method, "SR", sep="_")]][, idx_sr]), mean(results_obj[[paste(method, "MD", sep="_")]][, idx_md]), mean(results_obj[[paste(method, "MC", sep="_")]][, idx_mc]))
  std_pred_returns_star <- c(sd(results_obj$gen_UNIF), sd(results_obj[[paste(method, "SR", sep="_")]][, idx_sr]), sd(results_obj[[paste(method, "MD", sep="_")]][, idx_md]), sd(results_obj[[paste(method, "MC", sep="_")]][, idx_mc]))
  sharpe_ratios_star <- mean_pred_returns_star/std_pred_returns_star * sqrt(252)
  non_zero_elements <- c(nz_calc(results_obj$ws_gen_UNIF), nz_calc(results_obj[[paste("ws_SR", method, sep="_")]][, idx_sr]), nz_calc(results_obj[[paste("ws_MD", method, sep="_")]][, idx_md]), nz_calc(results_obj[[paste("ws_MC", method, sep="_")]][, idx_mc]))
  resulting_df_star <- rbind(mean_pred_returns_star, std_pred_returns_star, sharpe_ratios_star, non_zero_elements)
  
  mean_pred_returns <- c(mean(results_obj$gen_UNIF), mean(results_obj$gen_SR), mean(results_obj$gen_MD), mean(results_obj$gen_MC))
  std_pred_returns <- c(sd(results_obj$gen_UNIF), sd(results_obj$gen_SR), sd(results_obj$gen_MD), sd(results_obj$gen_MC))
  sharpe_ratios <- mean_pred_returns/std_pred_returns * sqrt(252)

  resulting_df <- rbind(mean_pred_returns, std_pred_returns, sharpe_ratios)
  total_df <- rbind(resulting_df, resulting_df_star)
  
  total_df = as.data.frame(total_df)
  names(total_df) <- c("UNIF", "SR", "MD", "MC")
  rownames(total_df) <- c("MR", "STD", "SR", "MR*", "STD*", "SR*", "NonZero")
  total_df <- round(total_df, 4)
  total_df
}
  
  



ts_multi_period_predict <- function(tbl, test_tbl, n_folds, optimal_constraints){
  
  folds = ts_crossValidation(tbl, n_folds)
  last_tr_fold = folds[[n_folds]]
  fold_length = folds[["length_fold"]]
  
  m = nrow(test_tbl)
  n = m %/% fold_length
  r = m %% fold_length
  if (r > 0){
    steps = n+1
  } else{
    steps = n
  }
  
  all_pmd_MD_returns = c()
  all_algo1_MD_returns = c()
  all_modpmd_MD_returns = c()
  all_fbg_MD_returns = c()
  
  all_pmd_SR_returns = c()
  all_algo1_SR_returns = c()
  all_modpmd_SR_returns = c()
  
  all_algo1_MC_returns = c()
  all_modpmd_MC_returns = c()
  
  all_md_returns = c()
  all_mc_returns = c()
  all_sr_returns = c()
  all_unif_returns = c()
  
  
  for (i in (1:steps)){
    
    if (i==1){
      tr_tbl = last_tr_fold
      val_tbl = test_tbl[(1+(i-1)*fold_length):(i*fold_length),]
      
    } else if (i < steps) {
      tr_tbl = val_tbl
      val_tbl = test_tbl[(1+(i-1)*fold_length):(i*fold_length),]
    } else {
      tr_tbl = val_tbl
      val_tbl = test_tbl[(1+(i-1)*fold_length):m,]
    }
    
    avg_ret = apply(tr_tbl, 2, mean) * 100
    sigma_tr = as.vector(apply(tr_tbl,2,sd))
    
    sigma2_mtrx = sigma_tr%*%t(sigma_tr)
    return2_mtrx = avg_ret%*%t(avg_ret)
    cov_mtrx <- cov(tr_tbl, use="pairwise.complete.obs")
    diag_sigma2_mtrx = diag(sigma_tr*sigma_tr)
    
    sqrt_obj <- sqrtm(cov_mtrx)
    sqrt_inv <- sqrt_obj$Binv
    sqrt_dir <- sqrt_obj$B
    sigma_over_tr <- sqrt_inv %*% sigma2_mtrx %*% sqrt_inv
    
    
    # Maximum Diversification portfolios
    pmd_obj <- PMD(sigma2_mtrx, sumabsv=optimal_constraints[["md_pmd"]], niter=100, K=1, center=F, trace=F)
    ws_MD_pmd <- pmd_obj$v
    ws_MD_pmd <- ws_MD_pmd / sum(abs(ws_MD_pmd))
    ws_MD_pmd <- ws_MD_pmd * sign(sum(ws_MD_pmd))
    
    algo1_obj <- algo1(sigma2_mtrx, cov_mtrx, cov_mtrx, lambda_u=optimal_constraints[["md_algo1"]], 
                       lambda_v=optimal_constraints[["md_algo1"]])
    ws_MD_algo1 <- algo1_obj$v
    
    modpmd_obj <- modified_PMD(sigma2_mtrx, cov_mtrx, cov_mtrx, 
                               lambda_u=optimal_constraints[["md_modpmd"]], lambda_v=optimal_constraints[["md_modpmd"]])
    ws_MD_modpmd <- modpmd_obj$v
    
    fbg_obj <- FBG_plusplus(sigma_over_tr, sqrt_dir, lambda=optimal_constraints[["md_fbg"]])
    ws_MD_fbg <- fbg_obj$v
    
    
    # Sharpe Ratio portfolios
    pmd_obj <- PMD(return2_mtrx, sumabsv=optimal_constraints[["sr_pmd"]],niter=100, K=1, center=F, trace=F)
    ws_SR_pmd <- pmd_obj$v
    ws_SR_pmd <- ws_SR_pmd / sum(abs(ws_SR_pmd))
    ws_SR_pmd <- ws_SR_pmd * sign(sum(ws_SR_pmd))
    
    algo1_obj <- algo1(return2_mtrx, cov_mtrx, cov_mtrx, lambda_u=optimal_constraints[["sr_algo1"]], 
                       lambda_v=optimal_constraints[["sr_algo1"]])
    ws_SR_algo1 <- algo1_obj$v
    
    modpmd_obj <- modified_PMD(return2_mtrx, cov_mtrx, cov_mtrx, lambda_u=optimal_constraints[["sr_modpmd"]], 
                               lambda_v=optimal_constraints[["sr_modpmd"]])
    ws_SR_modpmd <- modpmd_obj$v
    
    # Minimum Concentration portfolios
    algo1_obj <- algo1(sigma2_mtrx, diag_sigma2_mtrx, diag_sigma2_mtrx, 
                       lambda_u=optimal_constraints[["mc_algo1"]], lambda_v=optimal_constraints[["mc_algo1"]])
    ws_MC_algo1 <- algo1_obj$v
    
    modpmd_obj <- modified_PMD(sigma2_mtrx, diag_sigma2_mtrx, diag_sigma2_mtrx, 
                               lambda_u=optimal_constraints[["mc_modpmd"]], lambda_v=optimal_constraints[["mc_modpmd"]])
    ws_MC_modpmd <- modpmd_obj$v
    
    # additional portfolios
    ws_md = generic_portfolio_constructor("MD", cov_mtrx, sigma_tr, avg_ret)
    ws_mc = generic_portfolio_constructor("MC", cov_mtrx, sigma_tr, avg_ret)
    ws_sr = generic_portfolio_constructor("SR", cov_mtrx, sigma_tr, avg_ret)
    ws_unif = generic_portfolio_constructor("UN", cov_mtrx, sigma_tr, avg_ret)
    
    
    pmd_MD_returns <- val_tbl %*% ws_MD_pmd
    algo1_MD_returns <- val_tbl %*% ws_MD_algo1
    modpmd_MD_returns <- val_tbl %*% ws_MD_modpmd
    fbg_MD_returns <- val_tbl %*% ws_MD_fbg
    
    pmd_SR_returns <- val_tbl %*% ws_SR_pmd
    algo1_SR_returns <- val_tbl %*% ws_SR_algo1
    modpmd_SR_returns <- val_tbl %*% ws_SR_modpmd
    
    algo1_MC_returns <- val_tbl %*% ws_MC_algo1
    modpmd_MC_returns <- val_tbl %*% ws_MC_modpmd
    
    md_returns <- val_tbl %*% ws_md
    mc_returns <- val_tbl %*% ws_mc
    sr_returns <- val_tbl %*% ws_sr
    unif_returns <- val_tbl %*% ws_unif
    
    
    all_pmd_MD_returns = c(all_pmd_MD_returns, pmd_MD_returns)
    all_algo1_MD_returns = c(all_algo1_MD_returns, algo1_MD_returns)
    all_modpmd_MD_returns = c(all_modpmd_MD_returns, modpmd_MD_returns)
    all_fbg_MD_returns = c(all_fbg_MD_returns, fbg_MD_returns)
    
    all_pmd_SR_returns = c(all_pmd_SR_returns, pmd_SR_returns)
    all_algo1_SR_returns = c(all_algo1_SR_returns, algo1_SR_returns)
    all_modpmd_SR_returns = c(all_modpmd_SR_returns, modpmd_SR_returns)
    
    all_algo1_MC_returns = c(all_algo1_MC_returns, algo1_MC_returns)
    all_modpmd_MC_returns = c(all_modpmd_MC_returns, modpmd_MC_returns)
    
    all_md_returns = c(all_md_returns, md_returns)
    all_mc_returns = c(all_mc_returns, mc_returns)
    all_sr_returns = c(all_sr_returns, sr_returns)
    all_unif_returns = c(all_unif_returns, unif_returns)
    
    
  }
  
  list("pmd_MD"=all_pmd_MD_returns, "algo1_MD"=all_algo1_MD_returns,
       "modpmd_MD"=all_modpmd_MD_returns, "fbg_MD"=all_fbg_MD_returns,
       #"ws_MD_pmd"=pmd_MD_weights, "ws_MD_algo1" = algo1_MD_weights,
       #"ws_MD_cca"=mscca_MD_weights,
       
       "pmd_SR"=all_pmd_SR_returns, "algo1_SR"=all_algo1_SR_returns,
       "modpmd_SR"=all_modpmd_SR_returns,
       #"ws_SR_pmd"=pmd_SR_weights, "ws_SR_algo1" = algo1_SR_weights, 
       
       "pmd_MC"=all_pmd_MD_returns, "algo1_MC"=all_algo1_MC_returns,
       "modpmd_MC"=all_modpmd_MC_returns,
       #"ws_MC_pmd"=pmd_MD_weights, "ws_MC_algo1" = algo1_MC_weights)
       
       "md"=all_md_returns, "mc"=all_mc_returns,
       "sr"=all_sr_returns, "unif"=all_unif_returns)
  #"ws_sr"=ws_sr, "ws_sr_pmd"=ws_SR_pmd, "ws_sr_algo1"=ws_SR_algo1,
  #"ws_md"=ws_md, "ws_md_pmd"=ws_MD_pmd, "ws_md_algo1"=ws_MD_pmd,
  #"ws_mc_algo1"=ws_MC_algo1, "ws_mc"=ws_mc, "obj"=pmd_obj, "algo1_obj"=algo1_obj,
  #"mtrx"=return2_mtrx)
  
}



ts_predict <- function(tbl, test_tbl, n_folds, optimal_constraints){
  
  
  folds = ts_crossValidation(tbl, n_folds)
  tr_tbl = folds[[n_folds]]
  
  avg_ret = apply(tr_tbl, 2, mean) * 100
  sigma_tr = as.vector(apply(tr_tbl,2,sd))
  sigma2_mtrx = sigma_tr%*%t(sigma_tr)
  return2_mtrx = avg_ret%*%t(avg_ret)
  cov_mtrx <- cov(tr_tbl, use="pairwise.complete.obs")
  diag_sigma2_mtrx = diag(sigma_tr*sigma_tr)
  
  if(FALSE){
  # Maximum Diversification portfolios
  pmd_obj <- PMD(sigma2_mtrx, sumabsv=optimal_constraints[["md_pmd"]], niter=100, K=1, center=F, trace=F)
  ws_MD_pmd <- pmd_obj$v
  ws_MD_pmd <- ws_MD_pmd / sum(abs(ws_MD_pmd))
  ws_MD_pmd <- ws_MD_pmd * sign(sum(as.numeric(ws_MD_pmd) * avg_ret))
  #ws_MD_pmd <- ws_MD_pmd * sign(sum(ws_MD_pmd))
  
  algo1_obj <- algo1(sigma2_mtrx, cov_mtrx, cov_mtrx, avg_ret, lambda_u=optimal_constraints[["md_algo1"]], 
                     lambda_v=optimal_constraints[["md_algo1"]])
  ws_MD_algo1 <- algo1_obj$v
  
  modpmd_obj <- modified_PMD(sigma2_mtrx, cov_mtrx, cov_mtrx, avg_ret, lambda_u=optimal_constraints[["md_modpmd"]], 
                             lambda_v=optimal_constraints[["md_modpmd"]])
  ws_MD_modpmd <- modpmd_obj$v
  
  fbg_obj <- FBG_plusplus(cov_mtrx, sigma_tr, lambda=optimal_constraints[["md_fbg"]])
  ws_MD_fbg <- fbg_obj$v
  }
  sparse_obj <- sparsifier(cov_mtrx, sigma_tr, lambda=optimal_constraints[["md_sparse"]])
  ws_MD_sparse <- sparse_obj$v
  if(FALSE){
  # Sharpe Ratio portfolios
  pmd_obj <- PMD(return2_mtrx, sumabsv=optimal_constraints[["sr_pmd"]],niter=100, K=1, center=F, trace=F)
  ws_SR_pmd <- pmd_obj$v
  ws_SR_pmd <- ws_SR_pmd / sum(abs(ws_SR_pmd))
  ws_SR_pmd <- ws_SR_pmd * sign(sum(as.numeric(ws_SR_pmd) * avg_ret))
  #ws_SR_pmd <- ws_SR_pmd * sign(sum(ws_SR_pmd))
  
  algo1_obj <- algo1(return2_mtrx, cov_mtrx, cov_mtrx, avg_ret, lambda_u=optimal_constraints[["sr_algo1"]], 
                     lambda_v=optimal_constraints[["sr_algo1"]])
  ws_SR_algo1 <- algo1_obj$v
  
  modpmd_obj <- modified_PMD(return2_mtrx, cov_mtrx, cov_mtrx, avg_ret, lambda_u=optimal_constraints[["sr_modpmd"]], 
                             lambda_v=optimal_constraints[["sr_modpmd"]])
  ws_SR_modpmd <- modpmd_obj$v
  
  fbg_obj <- FBG_plusplus(cov_mtrx, avg_ret, lambda=optimal_constraints[["sr_fbg"]])
  ws_SR_fbg <- fbg_obj$v
  }
  sparse_obj <- sparsifier(cov_mtrx, avg_ret, lambda=optimal_constraints[["sr_sparse"]])
  ws_SR_sparse <- sparse_obj$v
  if(FALSE){
  # Minimum Concentration portfolios
  algo1_obj_mc <- algo1(sigma2_mtrx, diag_sigma2_mtrx, diag_sigma2_mtrx, avg_ret,
                        lambda_u=optimal_constraints[["mc_algo1"]], lambda_v=optimal_constraints[["mc_algo1"]])
  ws_MC_algo1 <- algo1_obj_mc$v
  
  modpmd_obj_mc <- modified_PMD(sigma2_mtrx, diag_sigma2_mtrx, diag_sigma2_mtrx, 
                                avg_ret, lambda_u=optimal_constraints[["mc_modpmd"]], lambda_v=optimal_constraints[["mc_modpmd"]])
  ws_MC_modpmd <- modpmd_obj_mc$v
  
  fbg_obj <- FBG_plusplus(diag_sigma2_mtrx, sigma_tr, lambda=optimal_constraints[["mc_fbg"]])
  ws_MC_fbg <- fbg_obj$v
  }
  sparse_obj <- sparsifier(diag_sigma2_mtrx, sigma_tr, lambda=optimal_constraints[["mc_sparse"]])
  ws_MC_sparse <- sparse_obj$v
  
  # additional portfolios
  ws_md = generic_portfolio_constructor("MD", cov_mtrx, sigma_tr, avg_ret)
  ws_mc = generic_portfolio_constructor("MC", cov_mtrx, sigma_tr, avg_ret)
  ws_sr = generic_portfolio_constructor("SR", cov_mtrx, sigma_tr, avg_ret)
  ws_unif = generic_portfolio_constructor("UN", cov_mtrx, sigma_tr, avg_ret)
  
  if(FALSE){
  pmd_MD_returns <- test_tbl %*% ws_MD_pmd
  algo1_MD_returns <- test_tbl %*% ws_MD_algo1
  modpmd_MD_returns <- test_tbl %*% ws_MD_modpmd
  fbg_MD_returns <- test_tbl %*% ws_MD_fbg
  }
  sparse_MD_returns <- test_tbl %*% ws_MD_sparse
  
  if(FALSE){
  pmd_SR_returns <- test_tbl %*% ws_SR_pmd
  algo1_SR_returns <- test_tbl %*% ws_SR_algo1
  modpmd_SR_returns <- test_tbl %*% ws_SR_modpmd
  fbg_SR_returns <- test_tbl %*% ws_SR_fbg
  }
  sparse_SR_returns <- test_tbl %*% ws_SR_sparse
  if(FALSE){
  algo1_MC_returns <- test_tbl %*% ws_MC_algo1
  modpmd_MC_returns <- test_tbl %*% ws_MC_modpmd
  fbg_MC_returns <- test_tbl %*% ws_MC_fbg
  }
  sparse_MC_returns <- test_tbl %*% ws_MC_sparse
  
  md_returns <- test_tbl %*% ws_md
  mc_returns <- test_tbl %*% ws_mc
  sr_returns <- test_tbl %*% ws_sr
  unif_returns <- test_tbl %*% ws_unif
  
  list(#"pmd_MD"=pmd_MD_returns, "algo1_MD"=algo1_MD_returns,
       #"modpmd_MD"=modpmd_MD_returns, "fbg_MD"=fbg_MD_returns,
       "sparse_MD"=sparse_MD_returns,
       
       #"pmd_SR"=pmd_SR_returns, "algo1_SR"=algo1_SR_returns,
       #"modpmd_SR"=modpmd_SR_returns, "fbg_SR"=fbg_SR_returns,
       "sparse_SR"=sparse_SR_returns,
       
       #"algo1_MC"=algo1_MC_returns, "modpmd_MC"=modpmd_MC_returns,
       #"fbg_MC"=fbg_MC_returns, 
       "sparse_MC"=sparse_MC_returns,
       
       "md"=md_returns, "mc"=mc_returns,
       "sr"=sr_returns, "unif"=unif_returns,
       
       #"ws_sr_pmd"=ws_SR_pmd, "ws_sr_algo1"=ws_SR_algo1,
       #"ws_sr_modpmd"=ws_SR_modpmd, "ws_sr_fbg"=ws_SR_fbg,
       "ws_sr_sparse"=ws_SR_sparse,
       
       #"ws_md_pmd"=ws_MD_pmd, "ws_md_algo1"=ws_MD_algo1,
       #"ws_md_modpmd"=ws_MD_modpmd, "ws_md_fbg"=ws_MD_fbg,
       "ws_md_sparse"=ws_MD_sparse,
       
       #"ws_mc_algo1"=ws_MC_algo1, "ws_mc_modpmd"=ws_MC_modpmd,
       #"ws_mc_fbg"=ws_MC_fbg, 
       "ws_mc_sparse"=ws_MC_sparse,
       
       "ws_mc"=ws_mc, "ws_sr"=ws_sr,"ws_md"=ws_md)
  
}
