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

nu = .5 #exponential
n = 2.5e4

points =SpatialPoints(cbind(runif(n, -6, 6), runif(n, -6, 6)))
p.coord = points@coords

water = Polygon(matrix(c(-6, -6, -6, 6, 6, 6, 6, -6, 0, -6, 0, 2, 0, -6, -6, -6),
                       ncol=2, byrow=T)) #u design
water = Polygons(list(water), 'water')
water = SpatialPolygons(list(water=water))

save.image('u_points.Rdata')



