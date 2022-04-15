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

insample_fraction = 0.9
n_folds = 20
n_stocks = 28


sumabsvs = linspace(1.0, sqrt(n_stocks), n=100)
algo1_lamdas = c(linspace(0.0, 0.1, n=100), linspace(0.101, 0.15, n=50))
fbg_lamdas = linspace(0.0, 0.999, n=100)
sparse_lamdas = linspace(0.0, 0.15, n=100)

resulting_returns_dax = optimal_constraints_searcher(daily_data, insample_fraction, n_folds)
temp = resulting_returns_dax
#xx = sharpe_by_fold(temp[["sparse_SR"]], 1966, function(x) {mean(x)/sd(x)})
#plot(sparse_lamdas, apply(xx, 2, mean))



constraints_star = list()
#constraints_star[["md_pmd"]] = sumabsvs[which.max(apply(temp$pmd_MD, 2, sum))]
#constraints_star[["md_algo1"]]  = algo1_lamdas[which.max(apply(temp$algo1_MD, 2, sum))]
#constraints_star[["md_modpmd"]]  = algo1_lamdas[which.max(apply(temp$modpmd_MD, 2, sum))]
#constraints_star[["md_fbg"]]  = fbg_lamdas[which.max(apply(temp$fbg_MD[,-100], 2, sum))]
constraints_star[["md_sparse"]]  = 0.02#sparse_lamdas[which.max(apply(temp$sparse_MD[,-100], 2, sum))]


#constraints_star[["sr_pmd"]] = sumabsvs[which.max(apply(temp$pmd_SR, 2, sum))]
#constraints_star[["sr_algo1"]] = algo1_lamdas[which.max(apply(temp$algo1_SR, 2, sum))]
#constraints_star[["sr_modpmd"]] = algo1_lamdas[which.max(apply(temp$modpmd_SR, 2, sum))]
#constraints_star[["sr_fbg"]] = algo1_lamdas[which.max(apply(temp$fbg_SR, 2, sum))]
constraints_star[["sr_sparse"]]  = 0.03 #sparse_lamdas[which.max(apply(temp$sparse_SR[,-100], 2, sum))]


#constraints_star[["mc_algo1"]] = algo1_lamdas[which.max(apply(temp$algo1_MC, 2, sum))]
#constraints_star[["mc_modpmd"]] = algo1_lamdas[which.max(apply(temp$modpmd_MC, 2, sum))]
#constraints_star[["mc_fbg"]] = algo1_lamdas[which.max(apply(temp$modpmd_MC, 2, sum))]
#constraints_star[["md_fbg"]]  = fbg_lamdas[which.max(apply(temp$fbg_MC[,-100], 2, sum))]
constraints_star[["mc_sparse"]]  = 0.035 #sparse_lamdas[which.max(apply(temp$sparse_MC[,-100], 2, sum))]


dax_oos_stats = results_saver(daily_data, insample_fraction, n_folds, constraints_star)
dax_oos_stats

dax_val_stats = results_in_sample_saver(resulting_returns_dax, constraints_star)
dax_val_stats



#resulting_returns_dax_IS = optimal_constraints_searcher(daily_data, insample_fraction, n_folds, use_train=TRUE)
#temp=resulting_returns_dax_IS

#constraints_star = list()
#constraints_star[["md_pmd"]] = sumabsvs[which.max(apply(temp$pmd_MD, 2, sum))]
#constraints_star[["md_algo1"]]  = 0.07#algo1_lamdas[which.max(apply(temp$algo1_MD, 2, sum))]
#constraints_star[["md_modpmd"]]  = algo1_lamdas[which.max(apply(temp$modpmd_MD, 2, sum))]
#constraints_star[["md_fbg"]]  = fbg_lamdas[which.max(apply(temp$fbg_MD, 2, sum))]
#constraints_star[["md_sparse"]]  = sparse_lamdas[which.max(apply(temp$sparse_MD, 2, mean)/apply(temp$sparse_MD, 2, sd))]


#constraints_star[["sr_pmd"]] = sumabsvs[which.max(apply(temp$pmd_SR, 2, sum))]
#constraints_star[["sr_algo1"]] = 0.0 #algo1_lamdas[which.max(apply(temp$algo1_SR, 2, sum))]
#constraints_star[["sr_modpmd"]] = algo1_lamdas[which.max(apply(temp$modpmd_SR, 2, sum))]
#constraints_star[["sr_fbg"]]  = fbg_lamdas[which.max(apply(temp$fbg_SR, 2, sum))]
#constraints_star[["sr_sparse"]]  = sparse_lamdas[which.max(apply(temp$sparse_SR, 2, mean)/apply(temp$sparse_SR, 2, sd))]


#constraints_star[["mc_algo1"]] = 0.06#algo1_lamdas[which.max(apply(temp$algo1_MC, 2, sum))]
#constraints_star[["mc_modpmd"]] = algo1_lamdas[which.max(apply(temp$modpmd_MC, 2, sum))]
#constraints_star[["md_fbg"]]  = fbg_lamdas[which.max(apply(temp$fbg_MC, 2, sum))]
#constraints_star[["mc_sparse"]]  = sparse_lamdas[which.max(apply(temp$sparse_MC, 2, mean)/apply(temp$sparse_MC, 2, sd))]

#dax_is_stats = results_in_sample_saver(resulting_returns_dax_IS, constraints_star)



# Get CAC data
daily_data=fread("Downloads/CAC40_20070102_20190211.csv")

resulting_returns_cac = optimal_constraints_searcher(daily_data, insample_fraction, n_folds)
temp = resulting_returns_cac

#xx = sharpe_by_fold(resulting_returns_dax[["sparse_SR"]], 1966, function(x) {mean(x)/sd(x)})
#plot(sparse_lamdas, apply(xx, 2, mean))



constraints_star = list()
#constraints_star[["md_pmd"]] = sumabsvs[which.max(apply(temp$pmd_MD, 2, sum))]
#constraints_star[["md_algo1"]]  = algo1_lamdas[which.max(apply(temp$algo1_MD, 2, sum))]
#constraints_star[["md_modpmd"]]  = algo1_lamdas[which.max(apply(temp$modpmd_MD, 2, sum))]
#constraints_star[["md_fbg"]]  = fbg_lamdas[which.max(apply(temp$fbg_MD[,-100], 2, sum))]
constraints_star[["md_sparse"]]  = 0.05#sparse_lamdas[which.max(apply(temp$sparse_MD[,-100], 2, sum))]


#constraints_star[["sr_pmd"]] = sumabsvs[which.max(apply(temp$pmd_SR, 2, sum))]
#constraints_star[["sr_algo1"]] = algo1_lamdas[which.max(apply(temp$algo1_SR, 2, sum))]
#constraints_star[["sr_modpmd"]] = algo1_lamdas[which.max(apply(temp$modpmd_SR, 2, sum))]
#constraints_star[["sr_fbg"]] = algo1_lamdas[which.max(apply(temp$fbg_SR, 2, sum))]
constraints_star[["sr_sparse"]]  = 0.02 #sparse_lamdas[which.max(apply(temp$sparse_SR[,-100], 2, sum))]


#constraints_star[["mc_algo1"]] = algo1_lamdas[which.max(apply(temp$algo1_MC, 2, sum))]
#constraints_star[["mc_modpmd"]] = algo1_lamdas[which.max(apply(temp$modpmd_MC, 2, sum))]
#constraints_star[["mc_fbg"]] = algo1_lamdas[which.max(apply(temp$modpmd_MC, 2, sum))]
#constraints_star[["md_fbg"]]  = fbg_lamdas[which.max(apply(temp$fbg_MC[,-100], 2, sum))]
constraints_star[["mc_sparse"]]  = 0.03 #sparse_lamdas[which.max(apply(temp$sparse_MC[,-100], 2, sum))]


cac_oos_stats = results_saver(daily_data, insample_fraction, n_folds, constraints_star)
cac_oos_stats

cac_val_stats = results_in_sample_saver(resulting_returns_cac, constraints_star)
cac_val_stats


#resulting_returns_cac_IS = optimal_constraints_searcher(daily_data, insample_fraction, n_folds, use_train=TRUE)
#temp=resulting_returns_cac_IS

#constraints_star = list()
#constraints_star[["md_pmd"]] = sumabsvs[which.max(apply(temp$pmd_MD, 2, sum))]
#constraints_star[["md_algo1"]]  = 0.07#algo1_lamdas[which.max(apply(temp$algo1_MD, 2, sum))]
#constraints_star[["md_modpmd"]]  = algo1_lamdas[which.max(apply(temp$modpmd_MD, 2, sum))]
#constraints_star[["md_fbg"]]  = fbg_lamdas[which.max(apply(temp$fbg_MD, 2, sum))]
#constraints_star[["md_sparse"]]  = sparse_lamdas[which.max(apply(temp$sparse_MD, 2, mean)/apply(temp$sparse_MD, 2, sd))]


#constraints_star[["sr_pmd"]] = sumabsvs[which.max(apply(temp$pmd_SR, 2, sum))]
#constraints_star[["sr_algo1"]] = 0.0 #algo1_lamdas[which.max(apply(temp$algo1_SR, 2, sum))]
#constraints_star[["sr_modpmd"]] = algo1_lamdas[which.max(apply(temp$modpmd_SR, 2, sum))]
#constraints_star[["sr_fbg"]]  = fbg_lamdas[which.max(apply(temp$fbg_SR, 2, sum))]
#constraints_star[["sr_sparse"]]  = sparse_lamdas[which.max(apply(temp$sparse_SR, 2, mean)/apply(temp$sparse_SR, 2, sd))]


#constraints_star[["mc_algo1"]] = 0.06#algo1_lamdas[which.max(apply(temp$algo1_MC, 2, sum))]
#constraints_star[["mc_modpmd"]] = algo1_lamdas[which.max(apply(temp$modpmd_MC, 2, sum))]
#constraints_star[["md_fbg"]]  = fbg_lamdas[which.max(apply(temp$fbg_MC, 2, sum))]
#constraints_star[["mc_sparse"]]  = sparse_lamdas[which.max(apply(temp$sparse_MC, 2, mean)/apply(temp$sparse_MC, 2, sd))]


#cac_is_stats = results_in_sample_saver(resulting_returns_cac_IS, constraints_star)
#cac_is_stats



# Get FTSE data
daily_data=fread("Downloads/FTSE100_20070102_20190211.csv")
n_folds = 10
resulting_returns_ftse = optimal_constraints_searcher(daily_data, insample_fraction, n_folds)

#xx = sharpe_by_fold(resulting_returns_dax[["sparse_SR"]], 1966, function(x) {mean(x)/sd(x)})
#plot(sparse_lamdas, apply(xx, 2, mean))



constraints_star = list()
#constraints_star[["md_pmd"]] = sumabsvs[which.max(apply(resulting_returns_dax$pmd_MD, 2, sum))]
#constraints_star[["md_algo1"]]  = algo1_lamdas[which.max(apply(resulting_returns_dax$algo1_MD, 2, sum))]
#constraints_star[["md_modpmd"]]  = algo1_lamdas[which.max(apply(resulting_returns_dax$modpmd_MD, 2, sum))]
#constraints_star[["md_fbg"]]  = fbg_lamdas[which.max(apply(resulting_returns_dax$fbg_MD[,-100], 2, sum))]
constraints_star[["md_sparse"]]  = 0.02#sparse_lamdas[which.max(apply(resulting_returns_dax$sparse_MD[,-100], 2, sum))]


#constraints_star[["sr_pmd"]] = sumabsvs[which.max(apply(resulting_returns_dax$pmd_SR, 2, sum))]
#constraints_star[["sr_algo1"]] = algo1_lamdas[which.max(apply(resulting_returns_dax$algo1_SR, 2, sum))]
#constraints_star[["sr_modpmd"]] = algo1_lamdas[which.max(apply(resulting_returns_dax$modpmd_SR, 2, sum))]
#constraints_star[["sr_fbg"]] = algo1_lamdas[which.max(apply(resulting_returns_dax$fbg_SR, 2, sum))]
constraints_star[["sr_sparse"]]  = 0.02 #sparse_lamdas[which.max(apply(resulting_returns_dax$sparse_SR[,-100], 2, sum))]


#constraints_star[["mc_algo1"]] = algo1_lamdas[which.max(apply(resulting_returns_dax$algo1_MC, 2, sum))]
#constraints_star[["mc_modpmd"]] = algo1_lamdas[which.max(apply(resulting_returns_dax$modpmd_MC, 2, sum))]
#constraints_star[["mc_fbg"]] = algo1_lamdas[which.max(apply(resulting_returns_dax$modpmd_MC, 2, sum))]
#constraints_star[["md_fbg"]]  = fbg_lamdas[which.max(apply(resulting_returns_dax$fbg_MC[,-100], 2, sum))]
constraints_star[["mc_sparse"]]  = 0.017 #sparse_lamdas[which.max(apply(resulting_returns_dax$sparse_MC[,-100], 2, sum))]


ftse_oos_stats = results_saver(daily_data, insample_fraction, n_folds, constraints_star)
ftse_oos_stats

ftse_val_stats = results_in_sample_saver(resulting_returns_ftse, constraints_star)
ftse_val_stats

#resulting_returns_ftse_IS = optimal_constraints_searcher(daily_data, insample_fraction, n_folds, use_train=TRUE)
#temp=resulting_returns_ftse_IS

#constraints_star = list()
#constraints_star[["md_pmd"]] = sumabsvs[which.max(apply(temp$pmd_MD, 2, sum))]
#constraints_star[["md_algo1"]]  = 0.07#algo1_lamdas[which.max(apply(temp$algo1_MD, 2, sum))]
#constraints_star[["md_modpmd"]]  = algo1_lamdas[which.max(apply(temp$modpmd_MD, 2, sum))]
#constraints_star[["md_fbg"]]  = fbg_lamdas[which.max(apply(temp$fbg_MD, 2, sum))]
#constraints_star[["md_sparse"]]  = sparse_lamdas[which.max(apply(temp$sparse_MD, 2, mean)/apply(temp$sparse_MD, 2, sd))]


#constraints_star[["sr_pmd"]] = sumabsvs[which.max(apply(temp$pmd_SR, 2, sum))]
#constraints_star[["sr_algo1"]] = 0.0 #algo1_lamdas[which.max(apply(temp$algo1_SR, 2, sum))]
#constraints_star[["sr_modpmd"]] = algo1_lamdas[which.max(apply(temp$modpmd_SR, 2, sum))]
#constraints_star[["sr_fbg"]]  = fbg_lamdas[which.max(apply(temp$fbg_SR, 2, sum))]
#constraints_star[["sr_sparse"]]  = sparse_lamdas[which.max(apply(temp$sparse_SR, 2, mean)/apply(temp$sparse_SR, 2, sd))]


#constraints_star[["mc_algo1"]] = 0.06#algo1_lamdas[which.max(apply(temp$algo1_MC, 2, sum))]
#constraints_star[["mc_modpmd"]] = algo1_lamdas[which.max(apply(temp$modpmd_MC, 2, sum))]
#constraints_star[["md_fbg"]]  = fbg_lamdas[which.max(apply(temp$fbg_MC, 2, sum))]
#constraints_star[["mc_sparse"]]  = sparse_lamdas[which.max(apply(temp$sparse_MC, 2, mean)/apply(temp$sparse_MC, 2, sd))]


#ftse_is_stats = results_in_sample_saver(resulting_returns_cac_IS, constraints_star)
#ftse_is_stats









data_splits = data_transformer(daily_data, insample_fraction)
allDataWithZeros_IS = data_splits$IS
allDataWithZeros_OOS = data_splits$OOS
n_stocks = ncol(allDataWithZeros_OOS)
folds = ts_crossValidation(allDataWithZeros_IS, n_folds)

tr_tbl = folds[[1]]
avg_ret = apply(tr_tbl, 2, mean) * 100
sigma_tr = as.vector(apply(tr_tbl, 2, sd))
cov_mtrx <- cov(tr_tbl, use="pairwise.complete.obs")

sum(abs(folds[[2]] %*% resulting_returns_dax$ws_SR_sparse[1:28,] - resulting_returns_dax$sparse_SR[1:491,])) 

resulting_returns_dax$ws_SR_sparse[1:28,1] 

resulting_returns_dax$ws_gen_SR[1:28,1] 


