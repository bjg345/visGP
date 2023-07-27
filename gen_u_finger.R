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

A = matrix(nrow = nrow(p.coord), ncol=nrow(p.coord))
for(i in 1:250){
  for(j in 1:100){
      A[((i-1)*100+j), 1:((i-1)*100+j)] = readRDS(file.path('A_u', paste0('Avec_', i, '.rds')))[[j]]
    }
}
A[upper.tri(A)] = t(A)[upper.tri(A)]
diag(A) = 0
saveRDS(A, 'A_finger.rds')

grid = points@coords; rownames(grid)=c()

ground = as.vector(sapply(1:250, function(i) readRDS(file.path('ground_u', paste0('ground_', i, '.rds')))))

ground = (ground-mean(ground))/sd(ground) #standardize

dat = SpatialPointsDataFrame(coords=grid, data = data.frame(val=ground))
blocked =spatialBlock(dat, rows=5, cols=5, selection='checkerboard') #blocked test/train desgin

train.id = which(blocked$foldID==1)
test.id = which(blocked$foldID==2)

grid.train = grid[train.id,]
grid.test = grid[test.id,]

ground.train = ground[train.id]
ground.test = ground[test.id]

A.train = A[train.id, train.id]
A.test = A[test.id, train.id]

save.image('finger_data.Rdata')
