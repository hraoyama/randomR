require(foreach)
require(data.table)
require(plyr)
require(zoo)
require(xts)
require(lubridate)
require(stringr)
require(elasticnet)
require(PMA)
require(genlasso)
require(pracma)
source("all_functions.R")

# Get DAX data
daily_data=fread("Downloads/DAX30_20070102_20190211.csv")
daily_data=fread("Downloads/CAC40_20070102_20190211.csv")
daily_data=fread("Downloads/FTSE100_20070102_20190211.csv")

insample_fraction = 0.8
n_folds = 5

data_splits = data_transformer(daily_data, insample_fraction)
allDataWithZeros = data_splits$ALL
allDataWithZeros_IS = data_splits$IS
allDataWithZeros_OOS = data_splits$OOS
n_stocks = ncol(allDataWithZeros_OOS)

#bigSigma_IS = cov(allDataWithZeros_IS, use="pairwise.complete.obs")
#smallSigmaVec_IS = as.vector(apply(allDataWithZeros_IS,2,sd))
#smallmuVec_IS = as.vector(apply(allDataWithZeros_IS,2,mean))



sumabsvs = linspace(1.0, sqrt(ncol(allDataWithZeros_IS)), n=100)
algo1_lamdas = c(linspace(0.0, 0.1, n=100), linspace(0.101, 0.15, n=50))
#modpmd_lamdas = c(linspace(0.0, 0.3, n=100))
fbg_lamdas = linspace(0.0, 1, n=100)




#resulting_returns_dax_all <- ts_evaluate(allDataWithZeros, n_folds+1, sumabsvs, algo1_lamdas, fbg_lamdas)
#results_by_fold(resulting_returns_dax_all$pmd_MD, n_folds, function(x) {mean(x)/sd(x)})
#results_by_fold(resulting_returns_dax_all$pmd_MD, n_folds, mean)
#xx = sharpe_by_fold(resulting_returns_dax_all$pmd_MD, n_folds, mean)
#plot(sumabsvs, apply(xx, 2, sd))
#plot(sumabsvs, apply(xx, 2, mean))



resulting_returns_dax <- ts_evaluate(allDataWithZeros_IS, n_folds, sumabsvs, algo1_lamdas, fbg_lamdas)

#resulting_returns_cac <- ts_evaluate(allDataWithZeros_IS, n_folds, sumabsvs, algo1_lamdas, fbg_lamdas)

#resulting_returns_ftse <- ts_evaluate(allDataWithZeros_IS, n_folds, sumabsvs, algo1_lamdas, fbg_lamdas)



plot(apply(resulting_returns_dax$ws_MD_pmd, 2, nz_calc), apply(resulting_returns_dax$pmd_MD, 2, sum), cex=.3, xlab='amount of non-zero components', ylab='cumulative return', main='PMD')
plot(apply(resulting_returns_dax$ws_MD_algo1, 2, nz_calc), apply(resulting_returns_dax$algo1_MD, 2, sum), cex=.3, xlab='amount of non-zero components', ylab='cumulative return', main='Modified SCCA')
plot(apply(resulting_returns_dax$ws_MD_modpmd, 2, nz_calc), apply(resulting_returns_dax$modpmd_MD, 2, sum), cex=.3, xlab='amount of non-zero components', ylab='cumulative return', main='Modified PMD')
plot(apply(resulting_returns_dax$ws_MD_fbg, 2, nz_calc), apply(resulting_returns_dax$fbg_MD, 2, sum), cex=.3, xlab='amount of non-zero components', ylab='cumulative return', main='FBG++')


plot(apply(resulting_returns_dax$ws_MC_algo1, 2, nz_calc), apply(resulting_returns_dax$algo1_MC, 2, sum), cex=.3)
plot(apply(resulting_returns_dax$ws_MC_modpmd, 2, nz_calc), apply(resulting_returns_dax$modpmd_MC, 2, sum), cex=.3)


plot(apply(resulting_returns_dax$ws_SR_pmd, 2, nz_calc), apply(resulting_returns_dax$pmd_SR, 2, sum), cex=.3)
plot(apply(resulting_returns_dax$ws_SR_algo1, 2, nz_calc), apply(resulting_returns_dax$algo1_SR, 2, sum), cex=.3)
plot(apply(resulting_returns_dax$ws_SR_modpmd, 2, nz_calc), apply(resulting_returns_dax$modpmd_SR, 2, sum), cex=.3)







#plot(sumabsvs, apply(resulting_returns_dax$pmd_MD, 2, sum), cex=.3)
plot(algo1_lamdas, apply(resulting_returns_dax$algo1_MD, 2, sum), xaxp=c(0., 0.15, 15), cex=.3)
axis(3, algo1_lamdas, round(apply(resulting_returns_dax$ws_MD_algo1, 2, nz_calc), 2))

#plot(algo1_lamdas, apply(resulting_returns_dax$modpmd_MD, 2, sum), cex=.3)
#plot(fbg_lamdas, apply(resulting_returns_dax$fbg_MD, 2, sum), cex=.3)


#plot(sumabsvs, apply(resulting_returns_dax$pmd_SR, 2, sum), cex=.3)
plot(algo1_lamdas, apply(resulting_returns_dax$algo1_SR, 2, sum), xaxp=c(0., 0.15, 15), cex=.3)
axis(3, algo1_lamdas, round(apply(resulting_returns_dax$ws_SR_algo1, 2, nz_calc), 2))


#plot(algo1_lamdas, apply(resulting_returns_dax$modpmd_SR, 2, sum), cex=.3)

plot(algo1_lamdas, apply(resulting_returns_dax$algo1_MC, 2, sum),xaxp=c(0., 0.15, 15),  cex=.3)
axis(3, algo1_lamdas, round(apply(resulting_returns_dax$ws_MC_algo1, 2, nz_calc), 2))

#plot(algo1_lamdas, apply(resulting_returns_dax$modpmd_MC, 2, sum), cex=.3)

constraints_star = list()
constraints_star[["md_pmd"]] = sumabsvs[which.max(apply(resulting_returns_dax$pmd_MD, 2, sum))]
constraints_star[["md_algo1"]]  = 0.07#algo1_lamdas[which.max(apply(resulting_returns_dax$algo1_MD, 2, sum))]
constraints_star[["md_modpmd"]]  = algo1_lamdas[which.max(apply(resulting_returns_dax$modpmd_MD, 2, sum))]
constraints_star[["md_fbg"]]  = fbg_lamdas[which.max(apply(resulting_returns_dax$fbg_MD[,-100], 2, sum))]


constraints_star[["sr_pmd"]] = sumabsvs[which.max(apply(resulting_returns_dax$pmd_SR, 2, sum))]
constraints_star[["sr_algo1"]] = 0.13#algo1_lamdas[which.max(apply(resulting_returns_dax$algo1_SR, 2, sum))]
constraints_star[["sr_modpmd"]] = algo1_lamdas[which.max(apply(resulting_returns_dax$modpmd_SR, 2, sum))]

constraints_star[["mc_algo1"]] = 0.07#algo1_lamdas[which.max(apply(resulting_returns_dax$algo1_MC, 2, sum))]
constraints_star[["mc_modpmd"]] = algo1_lamdas[which.max(apply(resulting_returns_dax$modpmd_MC, 2, sum))]



test_pred_returns <- ts_predict(allDataWithZeros_IS, allDataWithZeros_OOS, n_folds,
                                constraints_star)


#sum(test_pred_returns$pmd_MD)
sum(test_pred_returns$algo1_MD)
#sum(test_pred_returns$modpmd_MD)
#sum(test_pred_returns$fbg_MD)
sum(test_pred_returns$md)

#sum(test_pred_returns$pmd_SR)
sum(test_pred_returns$algo1_SR)
#sum(test_pred_returns$modpmd_SR)
sum(test_pred_returns$sr)


sum(test_pred_returns$algo1_MC)
#sum(test_pred_returns$modpmd_MC)
sum(test_pred_returns$mc)
sum(test_pred_returns$unif)


test_pred_returns_multi <- ts_multi_period_predict(allDataWithZeros_IS, allDataWithZeros_OOS, n_folds,
                                                   constraints_star)


sum(test_pred_returns_multi$pmd_MD)
sum(test_pred_returns_multi$algo1_MD)
sum(test_pred_returns_multi$modpmd_MD)
sum(test_pred_returns_multi$fbg_MD)
sum(test_pred_returns_multi$md)

sum(test_pred_returns_multi$pmd_SR)
sum(test_pred_returns_multi$algo1_SR)
sum(test_pred_returns_multi$modpmd_SR)
sum(test_pred_returns_multi$sr)


sum(test_pred_returns_multi$algo1_MC)
sum(test_pred_returns_multi$modpmd_MC)
sum(test_pred_returns_multi$mc)
sum(test_pred_returns_multi$unif)



length(which(colSums(allDataWithZeros_OOS) < 0))



### CAC

plot(apply(resulting_returns_cac$ws_pmd, 2, nz_calc), apply(resulting_returns_cac$pmd, 2, sum), cex=.3)
plot(apply(resulting_returns_cac$ws_algo1, 2, nz_calc), apply(resulting_returns_cac$algo1, 2, sum), cex=.3)


plot(sumabsvs, apply(resulting_returns_cac$pmd, 2, sum), cex=.3)
plot(algo1_lamdas, apply(resulting_returns_cac$algo1, 2, sum), cex=.3)

### FTSE

plot(apply(resulting_returns_ftse$ws_pmd, 2, nz_calc), apply(resulting_returns_ftse$pmd, 2, sum), cex=.3)
plot(apply(resulting_returns_ftse$ws_algo1, 2, nz_calc), apply(resulting_returns_ftse$algo1, 2, sum), cex=.3)

plot(sumabsvs, apply(resulting_returns_ftse$pmd, 2, sum), cex=.3)
plot(algo1_lamdas, apply(resulting_returns_ftse$algo1, 2, sum), cex=.3)





algo1_lamd_star = algo1_lamdas[which.max(apply(resulting_returns_cac$algo1, 2, sum))]
sumabsv_star = sumabsvs[which.max(apply(resulting_returns_cac$pmd, 2, sum))]
test_pred_returns_cac <- ts_predict(allDataWithZeros_IS, allDataWithZeros_OOS, n_folds,
                                sumabsv_star, algo1_lamd_star)

sum(test_pred_returns_cac$pmd)
sum(test_pred_returns_cac$algo1)
sum(test_pred_returns_cac$md)/sum(abs(test_pred_returns_cac$ws_md))
sum(test_pred_returns_cac$mc)
sum(test_pred_returns_cac$sr)
sum(test_pred_returns_cac$unif)


## FTSE

algo1_lamd_star = algo1_lamdas[which.max(apply(resulting_returns_ftse$algo1, 2, sum))]
sumabsv_star = sumabsvs[which.max(apply(resulting_returns_ftse$pmd, 2, sum))]
test_pred_returns_ftse <- ts_predict(allDataWithZeros_IS, allDataWithZeros_OOS, n_folds,
                                    sumabsv_star, algo1_lamd_star)

sum(test_pred_returns_ftse$pmd)
sum(test_pred_returns_ftse$algo1)
sum(test_pred_returns_ftse$md)/sum(abs(test_pred_returns_ftse$ws_md))
sum(test_pred_returns_ftse$mc)
sum(test_pred_returns_ftse$sr)
sum(test_pred_returns_ftse$unif)


### Debugging ###
              
       
folds = ts_crossValidation(allDataWithZeros_IS, 5)
dim(folds[[5]])

val_tbl = folds[[i]]
tr_tbl = folds[[i-1]]
avg_ret = apply(tr_tbl, 2, mean)
sigma_tr = as.vector(apply(tr_tbl,2,sd))

M = sigma_tr%*%t(sigma_tr)
cov_mtrx <- cov(tr_tbl, use="pairwise.complete.obs")
MN = avg_ret%*%t(avg_ret)

sqrt_obj <- sqrtm(Cu)
sqrt_inv <- sqrt_obj$Binv
sqrt_dir <- sqrt_obj$B
sigma_over_tr <- sqrt_inv %*% M %*% sqrt_inv

avg_ret
which.max(rowSums(M))
avg_ret[which.max(rowSums(M))]



ws_all_pmd = foreach(sumabs = sumabsvs,.combine = cbind) %do% {
  pmd_obj <- PMD(M, sumabsv = sumabs, K=1, center = F, trace=F)
  w <- pmd_obj$v
  #w <- w * sign(sum(as.numeric(w) * avg_ret))
  return(w)
}





resulting_returns_memoized <- ts_evaluate(allDataWithZeros_IS, allDataWithZeros_OOS, sumabsvs, 
                                          algo1_lamdas, fbg_lamdas, memoize=TRUE)
plot(sumabsvs, apply(resulting_returns_memoized$pmd, 2, sum), cex=.3)
plot(algo1_lamdas, apply(resulting_returns_memoized$algo1, 2, sum), cex=.3)


lamd_u = lamd_v = 0.5
A <- smallSigmaVec_IS%*%t(smallSigmaVec_IS)
B <- bigSigma_IS
algo1_ans <- algo1(A, B, B, lamd_u, lamd_u)
algo1_ans$v



sqrt_obj <- sqrtm(bigSigma_IS)
sqrtbigSigmaInv_IS <- sqrt_obj$Binv
sqrtbigSigma_IS <- sqrt_obj$B
sigma_overline = sqrtbigSigmaInv_IS%*%smallSigmaVec_IS%*%t(smallSigmaVec_IS)%*%sqrtbigSigmaInv_IS
M <- sigma_overline %*% sqrtbigSigma_IS
lamd = 0.3

algo1_fbg_ans <- algo1_FBG(sigma_overline, sqrtbigSigma_IS, lamd)
algo1_fbg_ans$v

u_init <- rep(1/sqrt(30), 30)
u_init*1/as.numeric(sqrt(t(u_init) %*% B %*% u_init ))
