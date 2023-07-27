library(geostatsp)
library(boraGP)
library(sf)
library(foreach)
library(rgeos)
library(tidyverse)
library(gRbase)
library(igraph)

load('u_data.Rdata')

nu=sigma.sq=tau.sq=1
phi = .1

id = 25

set.seed(id)

if(id <= 1500){
  n = 250
} else if (id <= 3000){
  n = 1200
} else (n = 10000)

if((id %% 1500) <= 500 & (id %% 1500L) > 0L){
  stdev = .1
} else if ((id %% 1500) <= 1000 & (id %% 1500L) > 0L){
  stdev = .25
} else (stdev = 1)

n.samples = 10000
burn = 5000

ind.train = sample.int(10000, .8*n, replace = F)
noise.train = rnorm(.8*n, sd = stdev)

ind.test = sample.int(10000, .2*n, replace = F)
noise.test = rnorm(.2*n, sd = stdev)

grid.train = grid.train[ind.train,]
grid.test = grid.test[ind.test,]

ground.train = ground.train[ind.train]
ground.test = ground.test[ind.test]

bound = as(water, 'SpatialLines')

m = 10
ord = order(grid.train[,1])
barrier_nninfo_all <- barrier_neighbor(coords = grid.train, coords.0 =  grid.test, ord = ord,
                                                n.neighbors = m,
                                                barrier = st_as_sf(bound),
                                                cores = 1,
                                                verbose = T,
                                                debug = list(barrier_n.indx = NULL,
                                                            barrier_dist = NULL,
                                                             barrier_nn.indx.0_list = NULL,
                                                             barrier_dist0 = NULL,
                                                            ref_window = NULL,
                                                             nonref_window = NULL,
                                                             ref_fill = TRUE,
                                                             nonref_fill = TRUE))
saveRDS(barrier_nninfo_all, 'cov_comp_neigbhors.rds')

#BORA
b = list()
b[[1]] = 0
f = list()
f[[1]] = tau.sq + sigma.sq

N = length(ord)
coords = grid.train[ord,]
D.ord = as.matrix(dist(coords))

for(i in 2:nrow(grid.train)){
  ind = barrier_nninfo_all$barrier_n.indx[[i]]
  M = length(ind)
  dists = D.ord[ind, ind]
  CC = matern(dists, param=c(range=1/phi, variance=sigma.sq, shape=nu)) + tau.sq*diag(M)
  C =  matern(barrier_nninfo_all$barrier_dist[[i]], param=c(range=1/phi, variance=sigma.sq, shape=nu))
  b[[i]] = solve(CC) %*% C
  f[[i]] = tau.sq+sigma.sq-sum(b[[i]] * C)
}
                                                             
F = diag(unlist(f))
B = matrix(0, nrow=N, ncol=N)
for(i in 1:N) B[i, barrier_nninfo_all$barrier_n.indx[[i]]] = b[[i]]
    
cov.mat.est.bora = (solve(diag(length(ord))-B) %*% F %*% t(solve(diag(length(ord))-B)))[order(ord), order(ord)]

A = A.train[ind.train, ind.train]

covs.est.wat.bora = cov.mat.est.bora[which(A==1)]
covs.est.land.bora = cov.mat.est.bora[which(A==0 & row( cov.mat.est.bora) != col( cov.mat.est.bora))]



#visGP
D = as.matrix(dist(grid.train))
g.unchord = graph_from_adjacency_matrix(A, 'undirected')
g.chord = is_chordal(g.unchord, newgraph=T)$newgraph
rip = rip(as(g.chord, 'graphNEL'))
cliques = sapply(rip$cliques, as.numeric)
seps = sapply(rip$separators, as.numeric)
clique.prec = function(clique){
  mat = matrix(0, nrow = N, ncol = N)
  if(length(clique)==0) return(mat)
  sub.cov = matern(D[clique,clique], param=c(range=1/phi, variance=sigma.sq, shape=nu)) + tau.sq*diag(length(clique))
  sub.prec = solve(sub.cov)
  mat[clique, clique] = sub.prec
  return(mat)
}

cov.mat.est.visGP = solve(Reduce('+', lapply(cliques, clique.prec)) - Reduce('+', lapply(seps, clique.prec)))
covs.est.wat.visGP = cov.mat.est.visGP[which(A==1)]
covs.est.land.visGP = cov.mat.est.visGP[which(A==0 & row( cov.mat.est.visGP) != col( cov.mat.est.visGP))]


##True
get.center.dist = function(p) sqrt(p[1]^2+(p[2]-2)^2)
get.center.dist2 = function(i,j) {get.center.dist(grid.train[i,])+get.center.dist(grid.train[j,])}
D.center = matrix(NA, nrow=N, ncol=N)
for(i in 1:N){for(j in 1:N){D.center[i,j] = ifelse(i==j, 0, get.center.dist2(i,j))}}
D.geo = D*(A) + D.center*(1-A)

cov.mat = matern(D.geo, param=c(range=1/phi, variance=sigma.sq, shape=nu)) + tau.sq*diag(N)

covs.tru.wat = cov.mat[which(A==1)]
covs.tru.land = cov.mat[which(A==0 & row( cov.mat) != col( cov.mat))]

png(file='cov_comp_points.png')
plot(water, main='Locations', cex.main=1.5, cex.lab=1.5, cex.axis=1.5)
points(grid.train)
dev.off()

png(file='cov_comp_vars_water.png')
plot(covs.tru.wat, covs.est.wat.bora, xlab='Raw variance/covariance values', ylab='Induced variance/covariance values', main='Comparison of raw and modeled\n covariances on domain-connected points', col='blue', cex.main=1.5, cex.lab=1.5, cex.axis=1.5) 
points(cbind(covs.tru.wat, covs.est.wat.visGP), col='red')
legend('bottomright',legend=c("BORA-GP", "visGP"),
       col=c("blue", "red"), pch=1)
abline(a = 0, b = 1, col = "black", lty = 2, lwd=3)
dev.off()

png(file='cov_comp_vars_land.png')
plot(covs.tru.land, covs.est.land.bora, xlab='Raw (geodesic) variance/covariance values', ylab='Induced variance/covariance values', main='Comparison of raw and modeled\n covariances on non-domain-connected points', col='blue', cex.main=1.5, cex.lab=1.5, cex.axis=1.5, ylim=c(0, max(covs.est.land.visGP+.05))) 
points(cbind(covs.tru.land, covs.est.land.visGP), col='red')
legend('bottomright',legend=c("BORA-GP", "visGP"),
       col=c("blue", "red"), pch=1)
abline(a = 0, b = 1, col = "black", lty = 2, lwd=3)
dev.off()

png(file='cov_comp_ord.png')
plot(diag(cov.mat.est.bora)[ord], xlab='Nearest-neighbor ordering', ylab='BORA-GP-induced marginal variance', main='Impact of ordering on variance', cex.main=1.5, cex.lab=1.5, cex.axis=1.5, col='blue')
dev.off()


