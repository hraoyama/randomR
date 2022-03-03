rm(list=ls())
require(ggplot2)
require(ggthemes)
require(statnet)
require(coda)
require(igraph)
require(foreach)
require(pracma)
require(rgexf)
require(fields)
require(plsgenomics)
require(data.table)

numOfSampleGraphs=1E3
n=50
pspace=linspace(n^-1.5, 1.1*log(n)/n, n=20)

gs <- foreach(p=pspace) %do%
  {
    graphCollection = foreach(gg = 1:numOfSampleGraphs) %do%
      {
        erdos.renyi.game(n,p,type="gnp", directed=F, loops=F)
      }
    return(graphCollection)
  }

mypg = function(g,lt=NULL)
{
  if(!is.null(lt))
    plot(g, layout=lt(g))
  else
    plot(g)
}

adjMs <- foreach(gc = gs) %do%
  { return ( foreach(g = gc) %do% 
               { as_adjacency_matrix(g, sparse=F)}) }

eigenVls <- lapply(adjMs, function(x) {  lapply(x, function(z) { eigen(z, symmetric=T)$values }) })
eigenVcts <- lapply(adjMs, function(x) { lapply(x, function(z) { eigen(z, symmetric=T)$vectors }) })

# local centrality measures
btwns <- lapply(gs, function(x) {lapply(x, betweenness)})
btwns_std = lapply(btwns, function(xs) {  lapply(xs, function(x) {  x/((n-1)*(n-2)/2)    })   })
eigenCs <- lapply(gs, function(x) {  lapply(x, function(z) { eigen_centrality(z)$vector } )} )  
  
# local clustering measures
localClustering <- lapply(gs, function(x) {
  lapply(x, function(z)
  {
    g1 = transitivity(z, type="local")#
    idxs = which(is.nan(g1))
    if(length(idxs)>0)
    { g1[idxs] = 0.0 }
    return(g1)
  }) })   
    
# average global clustering per probability
glocalClustering = lapply(gs,
  function(x)
  {
    mean(unlist(lapply(x, function(z)
      { k = transitivity(z, type="global")
        return(if(is.nan(k)) 0.0 else k)
    })))
  } )


# ucc measure of individual notdes for each graph with variance for the number of eigenvalues (K) considered
iprSimple = function(inputMatrix) { return(colSums(eigen(inputMatrix)$vectors^4))}
ucc_k = lapply(adjMs,
               function(x)
               {
                 lapply(x, function(z)
                   {
                    eigs = eigen(z)
                    eigVals = eigs$values
                    eigVecs = eigs$vectors
                    eigValsCS = cumsum(eigVals)
                    eigValTimesU_sq = sweep(eigVecs^2,2,eigVals,"*")
                    numerator = apply(eigValTimesU_sq, 1, cumsum)
                    uccOfNodes = sweep(numerator, 1, eigValsCS,"/")
                    return(uccOfNodes)
                 })
               }               )

ipr_K = lapply(ucc_k,
               function(x)
               {
                lapply(x, function(z)
                {
                  rowSums(z^2)
                } )  
               })

iprKClustering = sapply(ipr_K, function(x) { apply(do.call(rbind, x), 2, mean, na.rm=T)}  )

plot(apply(sapply(btwns, sapply, mean), 2, mean)) ; grid()
plot(apply(sapply(btwns_std, sapply, mean), 2, mean)) ; grid()
plot(apply(sapply(eigenCs, sapply, mean), 2, mean)) ; grid()
plot(apply(sapply(localClustering, sapply, mean), 2, mean)); grid()
matrix.heatmap(iprKClustering[1:49,])

clusterCompareDT = data.table( "P" = pspace,
                               "AVG_LOCAL" = apply(sapply(localClustering, sapply, mean),2,mean),
                               "AVG_GLOBAL" = unlist(glocalClustering),
                               "AVG_IPR_1" = iprKClustering[1,],
                               "AVG_IPR_2" = iprKClustering[2,],
                               "AVG_IPR_3" = iprKClustering[3,],
                               "AVG_IPR_5" = iprKClustering[5,],
                               "AVG_IPR_25" = iprKClustering[25,],
                               "AVG_IPR_30" = iprKClustering[30,],
                               "AVG_IPR_35" = iprKClustering[35,],
                               "AVG_IPR_40" = iprKClustering[40,],
                               "AVG_IPR_45" = iprKClustering[45,])

clusterCompareDT = melt(clusterCompareDT, id.vars="P")
ggplot(clusterCompareDT, aes(x=P, y=value, group = variable, colour = variable)) + 
  geom_point() + geom_line(aes(lty=variable)) + theme_fivethirtyeight()
                               





  



