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
require(data.table)
require(igraph)
require(doParallel)

numTimeSeries_Space = round(linspace(2,50,n=10))
lengthOfTimeSeries_Space = round(linspace(20,1000,n=10))

test_frame = data.table(expand.grid(numTimeSeries_Space, lengthOfTimeSeries_Space))
colnames(test_frame) = c("NUM_SERIES", "LENGTH_SERIES")
test_frame = test_frame[LENGTH_SERIES > NUM_SERIES]

numSamples = 1E4
iprSimple = function(inputMatrix) { return (colSums(eigen(inputMatrix)$vectors^4)) }

num_cores = detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)

clusterExport(cl, c("numSamples","test_frame", "iprSimple"))
allIPRks <- foreach(rc = 1:nrow(test_frame), .packages=c("data.table")) %dopar% 
  {
    p = test_frame[rc, NUM_SERIES]
    n = test_frame[rc, LENGTH_SERIES]
    S = diag(p)
    draws = rWishart(numSamples, n, S)/n
    iprks = t(apply(draws, 3, iprSimple))
    return(list("NUM_SERIES"=p,
                "LENGTH_SERIES"=n,
                "IPR_k"=iprks))
  }

IPRKSims <- lapply( lapply (allIPRks, "[[", "IPR_k"), 
                    function(x) { structure(x, class = c("IPRKSim","matrix"))} )

iprkFrame = data.table("NUM_SERIES" = unlist(lapply(allIPRks, "[[", "NUM_SERIES")),
                       "LENGTH_SERIES" = unlist(lapply(allIPRks, "[[", "LENGTH_SERIES")),
                       "IPRKs" = IPRKSims
                       )

plot(density(unlist(iprkFrame[NUM_SERIES == 2 & LENGTH_SERIES == 129, IPRKs])))
plot(density(unlist(iprkFrame[NUM_SERIES == 7 & LENGTH_SERIES == 129, IPRKs])))
plot(density(unlist(iprkFrame[NUM_SERIES == 23 & LENGTH_SERIES == 129, IPRKs])))
plot(density(unlist(iprkFrame[NUM_SERIES == 50 & LENGTH_SERIES == 129, IPRKs])))

plot(density(unlist(iprkFrame[NUM_SERIES == 50 & LENGTH_SERIES == 1000, IPRKs])))
plot(density(unlist(iprkFrame[NUM_SERIES == 7 & LENGTH_SERIES == 1000, IPRKs])))
plot(density(unlist(iprkFrame[NUM_SERIES == 2 & LENGTH_SERIES == 1000, IPRKs])))

plot(density(unlist(iprkFrame[NUM_SERIES == 50 & LENGTH_SERIES == 456, IPRKs])))
plot(density(unlist(iprkFrame[NUM_SERIES == 7 & LENGTH_SERIES == 456, IPRKs])))
plot(density(unlist(iprkFrame[NUM_SERIES == 2 & LENGTH_SERIES == 456, IPRKs])))

# just individual IPR_k where k=1,10
plot(density(unlist(iprkFrame[NUM_SERIES == 18 & LENGTH_SERIES == 1000, IPRKs][[1]][,1] )))
plot(density(unlist(iprkFrame[NUM_SERIES == 18 & LENGTH_SERIES == 1000, IPRKs][[1]][,10] )))
plot(density(unlist(iprkFrame[NUM_SERIES == 18 & LENGTH_SERIES == 1000, IPRKs][[1]][,18] )))

plot(density(unlist(iprkFrame[NUM_SERIES == 2 & LENGTH_SERIES == 1000, IPRKs][[1]][,1] )))
plot(density(unlist(iprkFrame[NUM_SERIES == 2 & LENGTH_SERIES == 1000, IPRKs][[1]][,2] )))







