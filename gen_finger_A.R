id = as.numeric(commandArgs(trailingOnly = T))

if(!dir.exists('A_finger')) dir.create('A_finger')
if(file.exists( file.path('A_finger', paste0('Avec_', id, '.rds')) )) quit()


require(Rcpp)
require(RcppArmadillo)
require(geoR)
library(magrittr)
require(parallel)
require(rgeos)
require(mvtnorm)
require(ggplot2)
require(dplyr)
require(mgcv)
require(geostatsp)
require(Rfast)
require(sp)
require(TSP)
require(tidyverse)
require(igraph)
require(blockCV)
require(sfsmisc)
require(cppRouting)
require(BRISC)
require(MLmetrics)
library(gRbase)
require(viridis)
require(Rfast)
set.seed(434)
load('finger_points.Rdata')




con = function(i, p.grid, thresh = 1e5) {#takes in an integer and returns a vector of that point's adjacency with all other points
  require(sp)
  require(rdist)
  
  n = nrow(p.grid)
  
  dists = rdist::cdist(p.grid[1:i, , drop=F], p.grid[i, ,drop=F])
  within.thresh = which(dists < thresh)
  
  l <- vector("list", length(within.thresh))
  
  for (j in seq_along(within.thresh)) {
    l[[j]] <- Lines(list(Line(as.matrix(rbind(p.grid[i,], p.grid[within.thresh[j],])))), as.character(within.thresh[j]))
  }
  
  out = rep(F, i)
  out[within.thresh] = gCovers(water, SpatialLines(l), byid=T)
  
  return(out)
  
}

out.list = mclapply((1+(id-1)*100):(id*100), FUN = con, mc.cores = detectCores()/2, p.grid = p.coord)
saveRDS(out.list, file.path('A_finger', paste0('Avec_', id, '.rds')))
