install.packages(c("foreach", "data.table", "zoo", "xts", "PMA", "genlasso", 
                   "elasticnet", "lubridate", "stringr", "pracma"))

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



# Get DAX data
daily_data=fread("Downloads/DAX30_20070102_20190211.csv")
daily_data = daily_data[, c("Index",colnames(daily_data)[sapply(colnames(daily_data), function(x) { str_detect(x,"Close") } )]),with=F]
daily_data[,date:=ymd_hms(Index)]
daily_data[,Index:=NULL]

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

# center the matrix perhaps? 
# allDataWithZeros_IS  = matrix(scale(allDataWithZeros_IS, center=T,scale=F), nrow=nrow(allDataWithZeros_IS), ncol=ncol(allDataWithZeros_IS))
bigSigma_IS = cov(allDataWithZeros_IS, use="pairwise.complete.obs")
# write.table(bigSigma_IS,file="D:/papers/sparse_rayleigh/bidSigma_IS.csv", row.names=F, col.names=F, sep=",")


pcaDecomp_IS = eigen(bigSigma_IS, ncol(bigSigma_IS))
smallSigmaVec_IS = as.vector(apply(allDataWithZeros_IS,2,sd))
svd_bigSigma_IS = svd(bigSigma_IS)
sqrtbigSigma_IS = svd_bigSigma_IS[["u"]]%*%diag(sqrt(svd_bigSigma_IS[["d"]]))%*%t(svd_bigSigma_IS[["u"]])
# all(abs(elasticnet::rootmatrix(bigSigma_IS)-sqrtbigSigma_IS)<1e-12)

# all(abs(sqrtbigSigma_IS%*%sqrtbigSigma_IS - bigSigma_IS)<1e-12)
bigC_IS = diag(1.0/smallSigmaVec_IS)%*%bigSigma_IS%*%diag(1.0/smallSigmaVec_IS)
# all( (bigC_IS - cor(allDataWithZeros_IS, use="pairwise.complete.obs")) < 1e-12)

# generalized Rayleigh ratio III.A (8)
sqrtbigSigmaInv_IS = svd_bigSigma_IS[["u"]]%*%diag(1.0/sqrt(svd_bigSigma_IS[["d"]]))%*%t(svd_bigSigma_IS[["u"]])


# the power of svd...
# bigSigmaInv_IS = solve(bigSigma_IS)
# svd_bigSigmaInv_IS = svd(bigSigmaInv_IS)
# sqrtbigSigmaInv_IS_2 = svd_bigSigmaInv_IS[["u"]]%*%diag(sqrt(svd_bigSigmaInv_IS[["d"]]))%*%t(svd_bigSigmaInv_IS[["u"]])
# all(abs(sqrtbigSigmaInv_IS - sqrtbigSigmaInv_IS_2)<1e-12)
numerator_matrix = sqrtbigSigmaInv_IS%*%smallSigmaVec_IS%*%t(smallSigmaVec_IS)%*%sqrtbigSigmaInv_IS
numerator_eig = eigen(numerator_matrix,ncol(numerator_matrix))
# numerator_eig[["values"]][1]
wstar = sqrtbigSigmaInv_IS%*%numerator_eig[["vectors"]][,1]
wstar = wstar/sqrt(sum(wstar**2))
#  abs(sqrt(sum(wstar**2))-1.0)<1e-12
# cbind(wstar,wstar_closed)

bigSigmaInv_IS = solve(bigSigma_IS)
# closed form for w_star III.B (9)
# bigCInv_IS = solve(bigC_IS)
# svd_bigCInv_IS = svd(bigCInv_IS)
# sqrtbigCInv_IS = svd_bigCInv_IS[["u"]]%*%diag(sqrt(svd_bigCInv_IS[["d"]]))%*%t(svd_bigCInv_IS[["v"]])
# eVec = matrix(rep(1,nrow(sqrtbigCInv_IS)),nrow=ncol(sqrtbigCInv_IS))
# wstar_closed =  ( sqrtbigCInv_IS%*%eVec ) / as.numeric(sqrt(t(eVec)%*%bigCInv_IS%*%eVec))
#  abs(sqrt(sum(wstar_closed**2))-1.0)<1e-12

# closed form for w_star III.B (9)
wstar_closed =  ( bigSigmaInv_IS%*%smallSigmaVec_IS ) / as.numeric(sqrt(t(smallSigmaVec_IS)%*%bigSigmaInv_IS%*%smallSigmaVec_IS))
wstar_closed = wstar_closed/sqrt(sum(wstar_closed**2))

abs(wstar_closed) - abs(wstar)



rm(list=ls())
require(ggplot2)
require(ggthemes)
require(foreach)
require(data.table)
require(plyr)
require(zoo)
require(xts)
require(lubridate)
require(stringr)


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
spcaDecomp_IS = eigen(bigSigma_IS, ncol(bigSigma_IS))
smallSigmaVec_IS = as.vector(apply(allDataWithZeros_IS,2,sd))
svd_bigSigma_IS = svd(bigSigma_IS)
sqrtbigSigma_IS = svd_bigSigma_IS[["u"]]%*%diag(sqrt(svd_bigSigma_IS[["d"]]))%*%t(svd_bigSigma_IS[["u"]])
bigC_IS = diag(1.0/smallSigmaVec_IS)%*%bigSigma_IS%*%diag(1.0/smallSigmaVec_IS)
sqrtbigSigmaInv_IS = svd_bigSigma_IS[["u"]]%*%diag(1.0/sqrt(svd_bigSigma_IS[["d"]]))%*%t(svd_bigSigma_IS[["u"]])

bigSigma_OOS = cov(allDataWithZeros_OOS, use="pairwise.complete.obs")
smallSigmaVec_OOS = as.vector(apply(allDataWithZeros_OOS,2,sd))

# write.table(bigSigma_IS,file="D:/papers/sparse_rayleigh/DAX_covmatrix_IS.csv", sep=",", col.names=FALSE)
# write.table(apply(allDataWithZeros_IS,2,mean),file="D:/papers/sparse_rayleigh/DAX_meanVec_IS.csv", sep=",", col.names=T)
# write.table(bigSigma_OOS,file="D:/papers/sparse_rayleigh/DAX_covmatrix_OOS.csv", sep=",", col.names=FALSE)
# write.table(apply(allDataWithZeros_OOS,2,mean),file="D:/papers/sparse_rayleigh/DAX_meanVec_OOS.csv", sep=",", col.names=T)

# write.table(bigSigma_IS,file="D:/papers/sparse_rayleigh/CAC_covmatrix_IS.csv", sep=",", col.names=FALSE)
# write.table(apply(allDataWithZeros_IS,2,mean),file="D:/papers/sparse_rayleigh/CAC_meanVec_IS.csv", sep=",", col.names=T)
# write.table(bigSigma_OOS,file="D:/papers/sparse_rayleigh/CAC_covmatrix_OOS.csv", sep=",", col.names=FALSE)
# write.table(apply(allDataWithZeros_OOS,2,mean),file="D:/papers/sparse_rayleigh/CAC_meanVec_OOS.csv", sep=",", col.names=T)


# generalized Rayleigh ratio III.A (8)
numerator_matrix = sqrtbigSigmaInv_IS%*%smallSigmaVec_IS%*%t(smallSigmaVec_IS)%*%sqrtbigSigmaInv_IS
numerator_eig = eigen(numerator_matrix,ncol(numerator_matrix))
wstar = sqrtbigSigmaInv_IS%*%numerator_eig[["vectors"]][,1]
wstar = wstar/sqrt(sum(wstar**2))

bigSigmaInv_IS = solve(bigSigma_IS)
# closed form for w_star III.B (9)
# bigCInv_IS = solve(bigC_IS)
# svd_bigCInv_IS = svd(bigCInv_IS)
  
# sqrtbigCInv_IS = svd_bigCInv_IS[["u"]]%*%diag(sqrt(svd_bigCInv_IS[["d"]]))%*%t(svd_bigCInv_IS[["u"]])
# eVec = matrix(rep(1,nrow(sqrtbigCInv_IS)),nrow=ncol(sqrtbigCInv_IS))
# wstar_closed =  ( sqrtbigCInv_IS%*%eVec ) / as.numeric(sqrt(t(eVec)%*%bigCInv_IS%*%eVec))

wstar_closed =  ( bigSigmaInv_IS%*%smallSigmaVec_IS ) / as.numeric(sqrt(t(smallSigmaVec_IS)%*%bigSigmaInv_IS%*%smallSigmaVec_IS))
wstar_closed = wstar_closed/sqrt(sum(wstar_closed**2))


library(PMA)
require(pracma)
require(foreach)
# MDAP
sumabsvs = linspace(1.0, sqrt(ncol(allDataWithZeros_IS)), n=50)

# wstars_PMD with different restraints

xx <- PMD(allDataWithZeros_IS, sumabsv = sqrt(30), K=1, center = F)
avret <- apply(allDataWithZeros_IS, 2, mean)
dim(as.numeric(xx$v))
dim(avret)

wstars_PMD = foreach( sumabs = sumabsvs,.combine = cbind) %do% {
  out_IS <- PMD(allDataWithZeros_IS, sumabsv = sumabs, K=1, center = F)
  wstar_PMD <- sqrtbigSigmaInv_IS%*%out_IS$v
  wstar_PMD <- wstar_PMD/sqrt(sum(wstar_PMD**2))
  return(out_IS$v)
}

# OOS performance of wstars_PMD
MDPRatio <- function(wts, bigSigma, smallSigmaVec)
{
  # wts = wstar_PMD[,1,drop=F]; bigSigma= bigSigma_OOS; smallSigmaVec = as.matrix(smallSigmaVec_OOS,nrow=length(smallSigmaVec_OOS))
  as.numeric((t(wts) %*% smallSigmaVec %*% t(smallSigmaVec) %*% wts) / (t(wts)%*%bigSigma%*%wts))
}

minimumConcentration <- function(wts, smallSigmaVec)
{
  #smallSigmaVec_matrix = as.matrix(smallSigmaVec,nrow=length(smallSigmaVec))
  as.numeric(  (t(wts) %*% diag(smallSigmaVec*smallSigmaVec) %*% wts) / (t(wts)%*%smallSigmaVec%*%t(smallSigmaVec)%*%wts))
}


fReturnCalc = function(wts, allDataWithZerosLocal)
{
  sum(  unlist(wts) *  apply(allDataWithZerosLocal,2,function(x)  {   sum(x)   })   )
}

fVolCalc = function(wts, bigSigmaLocal)
{
  sqrt(as.numeric(t(wts)%*%bigSigmaLocal%*%wts))
}


PMD_OOS = data.table("constraint" = sumabsvs, 
                     "MDPRatio_IS" = apply(wstars_PMD,2, function(wts) {  MDPRatio(wts, bigSigma_IS, as.matrix(smallSigmaVec_IS,nrow=length(smallSigmaVec_IS)))  }),
                     "MinConcRatio_IS" = apply(wstars_PMD,2, function(wts) {  minimumConcentration(wts, smallSigmaVec_IS)  }),
                     "MDPRatio_OOS" = apply(wstars_PMD,2, function(wts) {  MDPRatio(wts, bigSigma_OOS, as.matrix(smallSigmaVec_OOS,nrow=length(smallSigmaVec_OOS)))  }),
                     "MinConcRatio_OOS" = apply(wstars_PMD,2, function(wts) {  minimumConcentration(wts, smallSigmaVec_OOS)  }),
                     "returns_IS" = apply(wstars_PMD,2, function(wts) {  fReturnCalc(wts, allDataWithZeros_IS)  }),
                     "vols_IS" = apply(wstars_PMD,2, function(wts) {  fVolCalc(wts, bigSigma_IS)  }),
                     "returns_OOS" = apply(wstars_PMD,2, function(wts) {  fReturnCalc(wts, allDataWithZeros_OOS)  }),
                     "vols_OOS" = apply(wstars_PMD,2, function(wts) {  fVolCalc(wts, bigSigma_OOS)  }),
                     "sharpe_IS" = apply(wstars_PMD,2, function(wts) {  fReturnCalc(wts, allDataWithZeros_IS)/fVolCalc(wts, bigSigma_IS)  }),
                     "sharpe_OOS" = apply(wstars_PMD,2, function(wts) {  fReturnCalc(wts, allDataWithZeros_OOS)/fVolCalc(wts, bigSigma_OOS)  })
)

# fwrite(PMD_OOS,"D:/papers/sparse_rayleigh/nochange_PMD_cac40.csv")

wstar_MDPRatio_OOS = MDPRatio(wstar, bigSigma_OOS, as.matrix(smallSigmaVec_OOS,nrow=length(smallSigmaVec_OOS)))
wstar_closed_MDPRatio_OOS = MDPRatio(wstar_closed, bigSigma_OOS, as.matrix(smallSigmaVec_OOS,nrow=length(smallSigmaVec_OOS)))

wstar_MDPRatio_IS = MDPRatio(wstar, bigSigma_IS, as.matrix(smallSigmaVec_IS,nrow=length(smallSigmaVec_IS)))
wstar_closed_MDPRatio_IS = MDPRatio(wstar_closed, bigSigma_IS, as.matrix(smallSigmaVec_IS,nrow=length(smallSigmaVec_IS)))


require(latex2exp)
require(gridExtra)
require(ggplot2)
require(ggthemes)


p1 <- ggplot(PMD_OOS, aes(x=constraint, y=MDPRatio_IS)) + geom_point() + ylim(0.0, wstar_MDPRatio_IS) + theme_gdocs() + 
  geom_hline(yintercept=wstar_MDPRatio_IS, col="red",size=1.5) + 
  annotate("text", x=1.1, y=wstar_MDPRatio_IS-0.1, label= TeX("w^{*}"), parse=TRUE, col="red",size=8)  + 
  geom_hline(yintercept=wstar_closed_MDPRatio_IS, col="blue",size=1.5) + 
  annotate("text", x=1.1, y=wstar_closed_MDPRatio_IS+0.1, label= TeX("w_{closed}^{*}"), parse=TRUE, col="blue",size=8)  +
  ggtitle("MDP Ratio (eq.4) in sample for different constraint levels of the PMD")
p1
p2 <- ggplot(PMD_OOS, aes(x=constraint, y=MDPRatio_OOS)) + geom_point() + ylim(0.0, 3.5) + theme_gdocs() + 
  geom_hline(yintercept=wstar_MDPRatio_OOS, col="red",size=1.5) +

annotate("text", x=1.1, y=wstar_MDPRatio_OOS+0.1, label= TeX("w^{*}"), parse=TRUE, col="red",size=8)  + 
  geom_hline(yintercept=wstar_closed_MDPRatio_OOS, col="blue",size=1.5) + 
  annotate("text", x=1.1, y=wstar_closed_MDPRatio_OOS+0.1, label= TeX("w_{closed}^{*}"), parse=TRUE, col="blue",size=8)  +
  ggtitle("MDP Ratio (eq.4) out of sample for different constraint levels of the PMD")
grid.arrange(p1,p2,nrow=1)


p3 = ggplot(PMD_OOS, aes(x=constraint, y=log10(MinConcRatio_IS))) + geom_point() +  
  theme_gdocs() +
  ggtitle("log10(Minimum Concentration Ratio) (eq.5) in sample for different constraint levels of the PMD")
p4 = ggplot(PMD_OOS, aes(x=constraint, y=log10(MinConcRatio_OOS))) + geom_point() +  
  theme_gdocs() +
  ggtitle("log10(Minimum Concentration Ratio) (eq.5) out of sample for different constraint levels of the PMD")
grid.arrange(p3,p4,nrow=1)


p5 = ggplot(PMD_OOS, aes(x=constraint, y=returns_IS)) + geom_point() +  theme_gdocs() +
  ggtitle("Returns in sample for different constraint levels of the PMD")
p6 = ggplot(PMD_OOS, aes(x=constraint, y=returns_OOS)) + geom_point() +  theme_gdocs() +
  ggtitle("Returns out of sample for different constraint levels of the PMD")
grid.arrange(p5,p6,nrow=1)


p7 = ggplot(PMD_OOS, aes(x=constraint, y=vols_IS)) + geom_point() +  theme_gdocs() +
  ggtitle("Vols in sample for different constraint levels of the PMD")
p8 = ggplot(PMD_OOS, aes(x=constraint, y=vols_OOS)) + geom_point() +  theme_gdocs() +
  ggtitle("Vols out of sample for different constraint levels of the PMD")
grid.arrange(p7,p8,nrow=1)

p9 = ggplot(PMD_OOS, aes(x=constraint, y=sharpe_IS)) + geom_point() +  theme_gdocs() +
  ggtitle("Sharpe in sample for different constraint levels of the PMD")
p10 = ggplot(PMD_OOS, aes(x=constraint, y=sharpe_OOS)) + geom_point() +  theme_gdocs() +
  ggtitle("Sharpe out of sample for different constraint levels of the PMD")
grid.arrange(p9,p10,nrow=1)
