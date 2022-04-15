
require(ggplot2)
require(ggthemes)
require(foreach)
require(data.table)
require(plyr)
require(zoo)
require(xts)
require(lubridate)
require(stringr)
require(microbenchmark)

# Get DAX data
daily_data=fread("Downloads/DAX30_20070102_20190211.csv")
# daily_data=fread("Downloads/CAC40_20070102_20190211.csv")

daily_data = daily_data[, c("Index",colnames(daily_data)[sapply(colnames(daily_data), function(x) { str_detect(x,"Close") } )]),with=F]
daily_data[,date:=ymd_hms(Index)]
daily_data[,Index:=NULL]
if("LHN.VX.Close" %in% colnames(daily_data))
{ daily_data[,LHN.VX.Close:=NULL] }


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
allDataWithZeros_OOS = allDataWithZeros[round(0.8*nrow(allDataWithZeros)):nrow(allDataWithZeros),]
allDataWithZeros_IS = allDataWithZeros[1:(round(0.8*nrow(allDataWithZeros))-1),]


bigSigma_IS = cov(allDataWithZeros_IS, use="pairwise.complete.obs")
pcaDecomp_IS = eigen(bigSigma_IS, ncol(bigSigma_IS))
smallSigmaVec_IS = as.vector(apply(allDataWithZeros_IS,2,sd))
svd_bigSigma_IS = svd(bigSigma_IS)
sqrtbigSigma_IS = svd_bigSigma_IS[["u"]]%*%diag(sqrt(svd_bigSigma_IS[["d"]]))%*%t(svd_bigSigma_IS[["u"]])
bigC_IS = diag(1.0/smallSigmaVec_IS)%*%bigSigma_IS%*%diag(1.0/smallSigmaVec_IS)
sqrtbigSigmaInv_IS = svd_bigSigma_IS[["u"]]%*%diag(1.0/sqrt(svd_bigSigma_IS[["d"]]))%*%t(svd_bigSigma_IS[["u"]])

bigSigma_OOS = cov(allDataWithZeros_OOS, use="pairwise.complete.obs")
smallSigmaVec_OOS = as.vector(apply(allDataWithZeros_OOS,2,sd))


# generalized Rayleigh ratio III.A (8)
numerator_matrix = sqrtbigSigmaInv_IS%*%smallSigmaVec_IS%*%t(smallSigmaVec_IS)%*%sqrtbigSigmaInv_IS
numerator_eig = eigen(numerator_matrix,ncol(numerator_matrix))
wstar = sqrtbigSigmaInv_IS%*%numerator_eig[["vectors"]][,1]
wstar = wstar/sqrt(sum(wstar**2))

bigSigmaInv_IS = solve(bigSigma_IS)
wstar_closed =  ( bigSigmaInv_IS%*%smallSigmaVec_IS ) / as.numeric(sqrt(t(smallSigmaVec_IS)%*%bigSigmaInv_IS%*%smallSigmaVec_IS))
wstar_closed = wstar_closed/sqrt(sum(wstar_closed**2))


make_basis <- function(k, p) replace(numeric(p), k, 1)

algo1 <- function(M, Cu, Cv, lambda_u, lambda_v, num_iter=100, eps=1e-10){
  
  m = nrow(M)
  n = ncol(M)
  
  u_init <- rep(1/sqrt(m), m)
  v_init <- rep(1/sqrt(n), n)
  iterate <- TRUE
  res = list()
  count <- 0
  idx = which.max(rowSums(M))
  
  lamd_u_scaler = 1 / as.numeric(sqrt(t(u_init) %*% Cu %*% u_init ))
  lamd_v_scaler = 1 / as.numeric(sqrt(t(v_init) %*% Cv %*% v_init ))
  
  while(iterate && (count < num_iter)) {
    
    u <- M %*% v_init
    u <- u / as.numeric(sqrt(t(u) %*% Cu %*% u ))
    u <- pmax(u - 0.5 * lambda_u * lamd_u_scaler, 0) * sign(u)
    if(sum(u^2) == 0){
      u_init = make_basis(idx, m)
      v_init = make_basis(idx, n)
      break
    }
    u <- u / as.numeric(sqrt(t(u) %*% Cu %*% u ))
    
    v <- t(M) %*% u
    v <- v / as.numeric(sqrt(t(v) %*% Cv %*% v ))
    v <- pmax(v - 0.5 * lambda_v * lamd_v_scaler, 0) * sign(v)
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

  res[["v"]] <- v_init / sum(v_init)
  res[["u"]] <- u_init / sum(u_init)
  
  res

}


FBG_plusplus <- function(sigma_over, sigma_sqrt, lambda, num_iter=100, eps=1e-10){
  
  n = ncol(sigma_over)
  v_init <- u_init <- rep(1/sqrt(n), n)
  iterate <- TRUE
  res = list()
  count <- 0
  
  M <- sigma_over %*% sigma_sqrt
  lamd_scaler = 2*max(t(u_init)%*% M)
  idx = which.max(rowSums(M))
  
  while(iterate && (count < num_iter)) {
    
    v = pmax(abs(t(M) %*% u_init) - 0.5*lambda*lamd_scaler, 0)*sign(t(M)%*%u_init)
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
  
  res[["v"]] <- as.numeric(v_init)/sum(v_init)
  res[["u"]] <- u_init
  res
  
}



modified_PMD <- function(M, Cu, Cv, lambda_u, lambda_v, num_iter=200, eps=1e-10){
  
  m = nrow(M)
  n = ncol(M)
  
  u_init <- rep(1/sqrt(m), m)
  v_init <- rep(1/sqrt(n), n)
  iterate <- TRUE
  res = list()
  count <- 0
  idx = which.max(rowSums(M))

  while(iterate && (count < num_iter)) {

    u <- M %*% v_init
    u <- pmax(u - 0.5 * lambda_u*sum(u), 0) * sign(u)
    if(sum(u*u) == 0){
      u_init = make_basis(idx, m)
      v_init = make_basis(idx, n)
      break
    }
    u <- u / as.numeric(sqrt(t(u) %*% Cu %*% u ))
    
    v <- t(M) %*% u
    v <- pmax(v - 0.5 * lambda_v*sum(v), 0) * sign(v)
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

  res[["v"]] <- v_init / sum(v_init)
  res[["u"]] <- u_init / sum(u_init)
  res
}



algo1_ans <- modified_PMD(A, B, B, 0.01, 0.01, num_iter=5)
algo1_ans$v



A <- smallSigmaVec_IS%*%t(smallSigmaVec_IS)
B <- eye(dim(bigSigma_IS)[1])
lamd_u <- 0.4# 3.11/sqrt(30)
lamd_v <- 3./sqrt(30)
#lamd_u_seq <- linspace(0, 3/sqrt(30), n=50)
# lamd_v_seq <- linspace(0, 3/sqrt(30), n=50)



algo1_ans <- algo1(A, B, B, lamd_u, lamd_u)
v = algo1_ans$v

cor(v, v2)


pmd_res <- PMD(A, sumabsv=1, sumabsu=1, K=1, center = F, trace=F)
pmd_res$v
sum(abs(pmd_res$v))

abc = eigen(A)
abc$vectors[c(11,24),1]



v3 <- pmd_res$v
v3
cor(v, v3)



A <- allDataWithZeros_IS
Bu <- eye(dim(allDataWithZeros_IS)[1])
Bv <- eye(dim(allDataWithZeros_IS)[2])
lamd_u <- 2.8/sqrt(851)
lamd_v <- 2*0.2553830#2.7975/sqrt(30)

algo1_ans <- algo1(A, Bu, Bv, lamd_u, lamd_v)
v = algo1_ans$v
v
sum(algo1_ans$u == 0)

cor(algo1_ans$u, pmd_res$u)


sumabsvs = linspace(1, sqrt(30), n=50)
pmd_res <- PMD(allDataWithZeros_IS, sumabsv=1.1, K=1, center = F, trace=F)
sum(pmd_res$u == 0)
pmd_res$v

test1 <- function(lamd_seq){
  
  res = c()
  for (lamd in lamd_seq){
    ans = algo1(A, B, B, lamd, lamd)
    res = c(res, 1)
  }
  res
}

test2 <- function(penalty_seq){
  
  res = c()
  for (penalty in penalty_seq){
    ans = PMD(A, sumabsv=penalty, K=1, center=F, trace=F)
    res = c(res, 1)
  }
  res
}


microbenchmark(test1(lamd_v_seq), test2(sumabsvs))



lamd_u <- 1.05365/sqrt(30)/100
lamd_v <- 1.05365/sqrt(30)/100
mpmd_ans <- modified_PMD(A, B, B, lamd_u, lamd_v)


system.time(mpmd_ans <- modified_PMD(A, B, B, lamd_u, lamd_v))
v = mpmd_ans$v
v


v <- v / as.numeric(sqrt(t(v) %*% v))
(t(v) %*% B %*% v)


my_res = eigen(A)

cor(my_res$vectors[, 1], v)
cor(wstar_closed, v)





FrankeWolfe <- function(M, Cu, Cv, mu_u, mu_v, tau_u, tau_v, C, num_iter=100, eps=1e-10){
  
  m = nrow(M)
  n = ncol(M)
  
  u_init <- rep(tau_u/m, m)
  v_init <- rep(tau_v/n, n)
  
  u_init <- runif(m) 
  v_init <- runif(n)
  u_init <- u_init/sum(u_init)*tau_u
  v_init <- v_init/sum(v_init)*tau_v

  iterate <- TRUE
  res = list()
  count <- 1
  
  while(iterate && (count <= num_iter)) {

    d_u <- -(M %*% v_init) + 2 * mu_u * Cu %*% u_init
    j_u <- which.max(abs(d_u))
    s_u <- -sign(d_u[j_u]) * tau_u * make_basis(j_u, m)

    d_v <- -(t(M) %*% u_init) + 2 * mu_v * Cv %*% v_init
    j_v <- which.max(abs(d_v))
    s_v <- -sign(d_v[j_v]) * tau_v * make_basis(j_v, n)
    
    g <- as.numeric((u_init - s_u) %*% d_u + (v_init - s_v) %*% d_v)

    if (g < eps) {
      iterate=FALSE
    } else {
      gamma <- min(g/C, 1)
      u <- (1-gamma) * u_init + gamma * s_u
      v <- (1-gamma) * v_init + gamma * s_v
      
      u_init <- u
      v_init <- v
      
    }
    
    count <- count + 1
    
  }
  
  res[["v"]] <- round(v_init, 3)
  res[["u"]] <- round(u_init, 3)
  print(count)
  res
  
}

A <- smallSigmaVec_IS%*%t(smallSigmaVec_IS)
B <- eye(dim(bigSigma_IS)[1])
tau_u <- 10
tau_v <- 2
mu_u <- 0.0001
mu_v <- 0.0001
C = 50

fw_ans <- FrankeWolfe(A, B, B, mu_u, mu_v, tau_u, tau_v, C, num_iter=1000)

fw_ans$u
sum(fw_ans$u)

cor(pmd_res$v, fw_ans$u, method="spearman")
cbind(pmd_res$v, fw_ans$u)


which.max(rowSums(A))
which.max(colSums(A))

rowSums(A)[24]
rowSums(A)[11]

pmd_res <- PMD(A, sumabsv=5, K=1, center = F, trace=F)
sum(abs(pmd_res$v))
sum(abs(pmd_res$v*pmd_res$v))

