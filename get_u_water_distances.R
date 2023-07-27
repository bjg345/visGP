require(Rcpp)
require(RcppArmadillo)
require(geoR)
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

load('u_points.Rdata')

water = Polygon(matrix(c(-6, -6, -6, 6, 6, 6, 6, -6, 0, -6, 0, 2, 0, -6, -6, -6),
                       ncol=2, byrow=T))
water = Polygons(list(water), 'water')
water = SpatialPolygons(list(water=water))


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

reflex = which(angles<0)

water.dist = function(p1, p2){
  
  points.temp = SpatialPoints(rbind(p1@coords,p2@coords, vert[reflex,]))
  
  f.temp = function(i, j) gCovers(water, list(line = Line(rbind(points.temp[i]@coords, points.temp[j]@coords)) %>% Lines(ID='_')) %>% SpatialLines)
  adj.temp = outer(1:length(points.temp), 1:length(points.temp), Vectorize(f.temp))
  
  costs = as.matrix(dist(points.temp@coords)) + (1-adj.temp)*1e50
  
  graph.temp = expand.grid(1:nrow(points.temp@coords), 1:nrow(points.temp@coords)) %>% 
    mutate(cost = apply(., 1, function(x) costs[x[1], x[2]])) %>% rename(from=Var1, to=Var2) %>%
    makegraph()
  
  
  return(get_distance_pair(graph.temp, 1,2))
  
}

n= nrow(grid)
u_water_distances = matrix(NA, nrow = n, n)

for(i in 1:n){
	for(j in 1:i){
		u_water_distances[i,j] = water.dist( SpatialPoints(matrix(grid[i,], nrow=1)),  SpatialPoints(matrix(grid[j,], nrow=1)))
	}
}

u_water_distances[upper.tri(u_water_distances)] = u_water_distances[lower.tri(u_water_distances)]

saveRDS(u_water_distances, 'u_water_distances.rds')
