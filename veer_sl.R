# install.packages(c("Depends","Imports","tidyverse","foreach","data.table"), quietly=T)
# install.packages(c("glmnet","randomForest","class","gam","gbm","nnet","polspline","MASS","e1071","stepPlr","arm","party","spls","LogicReg"), quietly=T)
# install.packages(c("nnls","multicore","SIS","BayesTree","quadprog","ipred","mlbench","rpart","caret","mda","earth"), quietly=T)
# install.packages("dsa", quietly=T)
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager", quietly=T)
# BiocManager::install(version = "3.14", quietly=T)
# BiocManager::install(c("Biobase"), quietly=T) # "GenomicFeatures", "AnnotationDbi"
# install.packages("SuperLearner")
# install.packages(c("xgboost","KernelKnn","bartMachine"))
# install.packages("rJava")
# install.packages("randomForestSRC")
# BiocManager::install(c("breastCancerNKI"), quietly=T)
# BiocManager::install(c("Biobase"), quietly=T)
# install.packages("biglasso")
# install.packages("mltools")

# sl_lib = c("SL.xgboost", "SL.randomForest", "SL.glmnet", "SL.ksvm",
#            "SL.bayesglm", "SL.caret", "SL.cforest", "SL.gam",
#            "SL.logreg", "SL.biglasso",
#            "SL.bartMachine", "SL.kernelKnn", "SL.rpartPrune", "SL.lm", "SL.mean")

# head(vdv[vdv$Censoring!=0,1:10],4)
# head(vdv[,1:10],4)

require(SuperLearner)
require(biglasso)
require(mltools)

# data(Boston, package = "MASS")
data(vdv, package = "randomForestSRC")


set.seed(1)
train_idxs = sort(sample(1:nrow(vdv),60))
vdv_train = vdv[train_idxs,]
vdv_test = vdv[-train_idxs,]


# "SL.nnet"
listWrappers()

# Fit elastic net with 5 different alphas: 0, 0.2, 0.4, 0.6, 0.8, 1.0.
# 0 corresponds to ridge and 1 to lasso.
enet = create.Learner("SL.glmnet", detailed_names = T,
                      tune = list(alpha = seq(0, 1, length.out = 5))) # enet$names,


# sl_lib = c("SL.xgboost", "SL.randomForest", enet$names, 
#             "SL.gam",  "SL.glm", "SL.bartMachine", "SL.lm", "SL.mean")
# result_all = SuperLearner(Y = vdv$Censoring, X = vdv[, 3:ncol(vdv)], SL.library = sl_lib)
# sl_xgboost = SuperLearner(Y = vdv$Censoring, X = vdv[, 3:ncol(vdv)], SL.library = "SL.xgboost")
# sl_randomForest = SuperLearner(Y = vdv$Censoring, X = vdv[, 3:ncol(vdv)],SL.library = "SL.randomForest")

sl_lib_binomial_a = c("SL.xgboost", "SL.randomForest", enet$names, "SL.lm", "SL.mean") # works
sl_lib_binomial_b = c("SL.xgboost", "SL.randomForest", enet$names, "SL.ksvm", "SL.biglasso", "SL.lm", "SL.mean") # works

# sl_lib_binomial = c("SL.xgboost", "SL.randomForest", enet$names, "SL.ksvm", "SL.biglasso", "SL.lm", "SL.mean") 
# result_all_binomial = SuperLearner(Y = vdv$Censoring, X = vdv[, 3:ncol(vdv)], 
#                                    family=binomial(), SL.library = sl_lib_binomial)
# sl_xgboost_binomial = SuperLearner(Y = vdv$Censoring, X = vdv[, 3:ncol(vdv)], family=binomial(), SL.library = "SL.xgboost")
# sl_randomForest_binomial = SuperLearner(Y = vdv$Censoring, X = vdv[, 3:ncol(vdv)], family=binomial(), SL.library = "SL.randomForest")
# sl_glmnet0_binomial = SuperLearner(Y = vdv$Censoring, X = vdv[, 3:ncol(vdv)], family=binomial(), SL.library = enet$names[1])
# sl_glmnet0p25_binomial = SuperLearner(Y = vdv$Censoring, X = vdv[, 3:ncol(vdv)], family=binomial(), SL.library = enet$names[2])
# sl_glmnet0p5_binomial = SuperLearner(Y = vdv$Censoring, X = vdv[, 3:ncol(vdv)], family=binomial(), SL.library = enet$names[3])
# sl_glmnet0p75_binomial = SuperLearner(Y = vdv$Censoring, X = vdv[, 3:ncol(vdv)], family=binomial(), SL.library = enet$names[4])
# sl_glmnet1_binomial = SuperLearner(Y = vdv$Censoring, X = vdv[, 3:ncol(vdv)], family=binomial(), SL.library = enet$names[5])
# sl_ksvm_binomial = SuperLearner(Y = vdv$Censoring, X = vdv[, 3:ncol(vdv)], family=binomial(), SL.library = "SL.ksvm")
# sl_bilasso_binomial = SuperLearner(Y = vdv$Censoring, X = vdv[, 3:ncol(vdv)], family=binomial(), SL.library = "SL.biglasso")
# sl_lm_binomial = SuperLearner(Y = vdv$Censoring, X = vdv[, 3:ncol(vdv)], family=binomial(), SL.library = "SL.lm")

require(mltools)
rmse_test <- function(pred_res, act_res) {
  return( sqrt(mean((as.vector(pred_res$pred)-as.vector(act_res$Censoring))^2))  )
}
rmse_CV_test <- function(pred_res, act_res) {
  return( sqrt(mean((as.vector(pred_res$SL.predict)-as.vector(act_res$Censoring))^2))  )
}

# SL FIT with CV
sl_lib_binomial = c("SL.xgboost", "SL.randomForest", enet$names, "SL.ksvm", "SL.biglasso", "SL.lm", "SL.mean") 
result_all_binomial_train = SuperLearner(Y = vdv_train$Censoring, 
                                         X = vdv_train[, 3:ncol(vdv_train)], 
                                         cvControl = list(V = 4, shuffle=T, stratifyCV=T),
                                         family=binomial(), 
                                         SL.library = sl_lib_binomial)

names(result_all_binomial_train$fitLibrary)

# individual model objects
xgboostObject <- result_all_binomial_train$fitLibrary$SL.xgboost$object  # xgb.Booster
#  predict(xgboostObject, as.matrix(vdv_train[, 3:ncol(vdv_train)]))
randomForestObject <- result_all_binomial_train$fitLibrary$SL.randomForest$object  
glmnet0Object <- result_all_binomial_train$fitLibrary$SL.glmnet_0_All$object  
glmnet025Object <- result_all_binomial_train$fitLibrary$SL.glmnet_0.25_All$object  
glmnet05Object <- result_all_binomial_train$fitLibrary$SL.glmnet_0.5_All$object  
glmnet075Object <- result_all_binomial_train$fitLibrary$SL.glmnet_0.75_All$object  
glmnet1Object <- result_all_binomial_train$fitLibrary$SL.glmnet_1_All$object  
ksvmObject <- result_all_binomial_train$fitLibrary$SL.ksvm$object  
biglassoObject <- result_all_binomial_train$fitLibrary$SL.biglasso$object  
lmObject <- result_all_binomial_train$fitLibrary$SL.lm$object  
meanObject <- result_all_binomial_train$fitLibrary$SL.mean$object  

# IS predictions of subsets... + GP TRAIN Data
fwrite(data.frame(cbind(result_all_binomial_train$library.predict, "Y"=vdv_train$Censoring)), 
       "/Users/hansroggeman/Code/veer_gp_train_data.csv")
# SL in sample prediction
pred_SL_train_IS = predict(result_all_binomial_train, vdv_train[, 3:ncol(vdv_train)])
as.vector(pred_SL_train_IS$pred)

# IS individual model predictions
IS_preds = apply(result_all_binomial_train$library.predict, 2, 
      function(x) {
        return(c("auc"=auc_roc(as.vector(x), vdv_train$Censoring),
                    "rmse"=rmse(as.vector(x), vdv_train$Censoring),
                    "acc"= sum(as.numeric(as.vector(x>=0.5))==vdv_train$Censoring)/nrow(vdv_train)
                    ))
      } )
IS_preds = data.frame(cbind(IS_preds, c("auc"=auc_roc(as.vector(pred_SL_train_IS$pred), vdv_train$Censoring),
                  "rmse"=rmse(as.vector(pred_SL_train_IS$pred), vdv_train$Censoring),
                  "acc"= sum(as.numeric(as.vector(pred_SL_train_IS$pred>=0.5))==vdv_train$Censoring)/nrow(vdv_train))))
colnames(IS_preds)[ncol(IS_preds)] = "SL.All"
fwrite(IS_preds, "/Users/hansroggeman/Code/veer_sl_IS_summary.csv")

# SL OOS prediction
pred_SL_train_OOS = predict(result_all_binomial_train, vdv_test[, 3:ncol(vdv_test)])
as.vector(pred_SL_train_OOS$pred)

# OOS predictions of subsets... + GP TEST Data
fwrite(data.frame(cbind(pred_SL_train_OOS$library.predict, "Y"=vdv_test$Censoring)), 
       "/Users/hansroggeman/Code/veer_gp_test_data.csv")

# OOS individual model predictions
OOS_preds = apply(pred_SL_train_OOS$library.predict, 2, 
                 function(x) {
                   return(c("auc"=auc_roc(as.vector(x), vdv_test$Censoring),
                            "rmse"=rmse(as.vector(x), vdv_test$Censoring),
                            "acc"= sum(as.numeric(as.vector(x>=0.5))==vdv_test$Censoring)/nrow(vdv_test)
                   ))
                 } )
OOS_preds = data.frame(cbind(OOS_preds, 
                             c("auc"=auc_roc(as.vector(pred_SL_train_OOS$pred), vdv_test$Censoring),
                                "rmse"=rmse(as.vector(pred_SL_train_OOS$pred), vdv_test$Censoring),
                                "acc"= sum(as.numeric(as.vector(pred_SL_train_OOS$pred>=0.5))==vdv_test$Censoring)/nrow(vdv_test))))
colnames(OOS_preds)[ncol(OOS_preds)] = "SL.All"
fwrite(OOS_preds, "/Users/hansroggeman/Code/veer_sl_OOS_summary.csv")

# 
# # # CV SL  summary run example
# # fitSL.CV <- CV.SuperLearner(Y = vdv_train$Censoring, X = vdv_train[, 3:ncol(vdv_train)], 
# #                             SL.library = sl_lib_binomial ,
# #                             V = 5, family = binomial(),
# #                             method = "method.NNLS",
# #                             cvControl = list(stratifyCV = TRUE))
# # summary(fitSL.CV)
# # plot(fitSL.CV)
# # fitSL.CV$coef
# 
# 



# result = SuperLearner(Y = Boston$medv, X = Boston[, -14], SL.library = sl_lib)
# result_binom = SuperLearner(Y = vdv$Censoring, X = vdv[, 3:ncol(vdv)], 
#                             family = binomial(),
#                             SL.library = sl_lib)




# Review performance of each algorithm and ensemble weights.
# result

# # Use external (aka nested) cross-validation to estimate ensemble accuracy.
# # This will take a while to run.
# result2 = CV.SuperLearner(Y = Boston$medv, X = Boston[, -14], SL.library = sl_lib)
# 
# # Plot performance of individual algorithms and compare to the ensemble.
# require(ggthemes)
# plot(result2) # + theme_minimal()
# 
# # Hyperparameter optimization --
# 
# sl_lib2 = c("SL.mean", "SL.lm", enet$names)
# 
# enet_sl = SuperLearner(Y = Boston$medv, X = Boston[, -14], SL.library = sl_lib2)
# 
# # Identify the best-performing alpha value or use the automatic ensemble.
# enet_sl
# 
# 
# enet_sl2 = CV.SuperLearner(Y = Boston$medv, X = Boston[, -14], SL.library = sl_lib2)
# 
# 
# plot(enet_sl2)



