id = as.numeric(commandArgs(trailingOnly = T))

if(file.exists(file.path('ground_finger', paste0('ground_', id, '.rds')))) quit()

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

if(!dir.exists('ground_finger')) dir.create('ground_finger')

vert = water@polygons$water@Polygons[[1]]@coords
rownames(vert) = 1:nrow(vert)
vert = vert[-nrow(vert),]
angles = vector('numeric', length = nrow(vert))
for(i in 1:nrow(vert)){
  
  vec1 =  vert[ifelse(i==nrow(vert), 1, i+1),] - vert[i,] 
  vec2 = vert[i,] - vert[ifelse(i==1, nrow(vert), i-1),] 
  
  dot = sum(vec1*vec2)
  det = det(rbind(vec1, vec2))
  
  angles[i] = base::atan2(det, dot)
  
}
angles[which(angles==pi)] = 0
reflex = which(angles<=0) #reflex angles of polygon, to find shortest paths (line must go through reflex angles only)


water.dist = function(p1, p2){
  
  points.temp = SpatialPoints(rbind(p1@coords,p2@coords, vert[reflex,]))
  
  f.temp = function(i, j) gCovers(water, list(line = Line(rbind(points.temp[i]@coords, points.temp[j]@coords)) %>% Lines(ID='_')) %>% SpatialLines)
  adj.temp = outer(1:length(points.temp), 1:length(points.temp), Vectorize(f.temp))
  
  costs = as.matrix(dist(points.temp@coords)) + (1-adj.temp)*1e50 #if not adjacent, direct cost is ~infinite
  
  graph.temp = expand.grid(1:nrow(points.temp@coords), 1:nrow(points.temp@coords)) %>% 
    mutate(cost = apply(., 1, function(x) costs[x[1], x[2]])) %>% rename(from=Var1, to=Var2) %>%
    makegraph()
  
  
  return(get_distance_pair(graph.temp, 1,2))
  
}

sources = cbind(c(-5, -3, -1, 1), 5)

fixed = function(s){
  
  dist1 = water.dist(SpatialPoints(matrix(sources[1,], nrow=1)), SpatialPoints(matrix(s, nrow=1)))
  dist2 = water.dist(SpatialPoints(matrix(sources[2,], nrow=1)), SpatialPoints(matrix(s, nrow=1)))
  dist3 = water.dist(SpatialPoints(matrix(sources[3,], nrow=1)), SpatialPoints(matrix(s, nrow=1)))
  dist4 = water.dist(SpatialPoints(matrix(sources[4,], nrow=1)), SpatialPoints(matrix(s, nrow=1)))
  return(dist1^2/3 + 3*sin(dist3) - dist2*dist4 )
}


ground = mclapply((1+(id-1)*100):(id*100), FUN = function(i) fixed(points@coords[i,]), mc.cores = detectCores()/2)  %>% unlist() #ground truth function values

saveRDS(ground, file.path('ground_finger', paste0('ground_', id, '.rds')))
