require(Rcpp)
require(RcppArmadillo)
require(geoR)
require(rgeos)
require(mvtnorm)
require(ggplot2)
require(parallel)
require(magrittr)
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

x.min <-  -10
x.max <- 10
y.min <- -10
y.max <- 10

boundary = Polygon(rbind(c(x.min, y.min), c(x.max, y.min), c(x.max, y.max), c(x.min, y.max)))
boundary = Polygons(list(boundary),'boundary')
boundary = SpatialPolygons(list(boundary=boundary)) #spatial boundaries


water = Polygon(matrix(c(-6, -6, -6, 6, -4.5, 6, -4.5, -5, -4, -5, -4, 6, -2.5, 6, -2.5, -5, -2, -5, -2, 6, -0.5, 6, -0.5, -5, 0, -5, 0, 6, 1.5, 6, 1.5,-6,-6, -6),
                       ncol=2, byrow=T)) #finger design
water = Polygons(list(water), 'water')
water = SpatialPolygons(list(water=water))

land = gDifference(boundary, water) 

nu = .5 #exponential
n = 2.5e4

points =SpatialPoints(cbind(runif(n*100, x.min, x.max), runif(n*100, y.min, y.max)))
landed = points[land]
points = sample(gDifference(points, landed), n)
p.coord = points@coords

save.image('finger_points.Rdata')

quit()
p1=Sys.time()
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

out.list = mclapply(1:nrow(p.coord), FUN = con, mc.cores = 24, p.grid = p.coord)
A = matrix(nrow = nrow(p.coord), ncol=nrow(p.coord))
for(i in 1:nrow(A)){
  A[i, 1:i] = out.list[[i]]
}
A[upper.tri(A)] = t(A)[upper.tri(A)]
diag(A) = 0
saveRDS(A, 'A_finger_nothresh.rds')
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

reflex = which(angles<0) #reflex angles of polygon, to find shortest paths (line must go through reflex angles only)


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
  dist3 = water.dist(SpatialPoints(matrix(sources[1,], nrow=1)), SpatialPoints(matrix(s, nrow=1)))
  dist4 = water.dist(SpatialPoints(matrix(sources[2,], nrow=1)), SpatialPoints(matrix(s, nrow=1)))
  return(dist1^2/3 + 3*sin(dist1) - dist3*dist4 ) #arbitrary complicated function of distance
}


grid = points@coords; rownames(grid)=c()

ground = mclapply(1:n, FUN = function(i) fixed(points@coords[i,]), mc.cores = 24) %>% unlist() #ground truth function values
ground = (ground-mean(ground))/sd(ground) #standardize

noise = rnorm(n, sd = .2)

vals = ground + noise


dat = SpatialPointsDataFrame(coords=grid, data = data.frame(val=vals))
blocked =spatialBlock(dat, rows=5, cols=5, selection='checkerboard') #blocked test/train desgin

train.id = which(blocked$foldID==1)
test.id = which(blocked$foldID==2)

grid.train = grid[train.id,]
grid.test = grid[test.id,]

vals.train = vals[train.id]
vals.test = vals[test.id]

ground.train = ground[train.id]
ground.test = ground[test.id]

A.train = A[train.id, train.id]
A.test = A[test.id, train.id]

save.image('finger_points_nothresh.Rdata')

