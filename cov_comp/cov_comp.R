library(geostatsp)
library(boraGP)
library(sf)
library(foreach)
library(rgeos)
library(tidyverse)
library(gRbase)
library(igraph)
library(sfsmisc)

load('u_data.Rdata')

my_theme <- theme_gray() +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 18), 
    legend.text = element_text(size=18),
    legend.title=element_text(size=18),
    legend.position = 'bottom'
  )
  
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
  C =  matern(D.ord[i, ind], param=c(range=1/phi, variance=sigma.sq, shape=nu))
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
V(g.unchord)$name = 1:N
g.chord = is_chordal(g.unchord, newgraph=T)$newgraph
rip = rip(g.chord)
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


##True geodesic
get.center.dist = function(p) sqrt(p[1]^2+(p[2]-2)^2)
get.center.dist2 = function(i,j) {get.center.dist(grid.train[i,])+get.center.dist(grid.train[j,])}
D.center = matrix(NA, nrow=N, ncol=N)
for(i in 1:N){for(j in 1:N){D.center[i,j] = ifelse(i==j, 0, get.center.dist2(i,j))}}
D.geo = D*(A) + D.center*(1-A)

cov.mat = matern(D.geo, param=c(range=1/phi, variance=sigma.sq, shape=nu)) + tau.sq*diag(N)

covs.tru.wat = cov.mat[which(A==1)]
covs.tru.land = cov.mat[which(A==0 & row( cov.mat) != col( cov.mat))]

# MDS
embedding <- stats::cmdscale(D.geo, k = 3) # mds embed into R^3
D.mds <- as.matrix(dist(embedding))
covs.est.mds <-  matern(D.mds, param=c(range=1/phi, variance=sigma.sq, shape=nu)) + tau.sq*diag(N)
covs.est.wat.mds <- covs.est.mds[which(A == 1)]
covs.est.land.mds <- covs.est.mds[which(A == 0 & row(covs.est.mds) != col(covs.est.mds))]

#plots
png(file='cov_comp_points.png')
plot(water, main='Locations', cex.main=2.75, cex.lab=2.75, cex.axis=2.75)
points(grid.train)
dev.off()

library(ggplot2)

library(ggplot2)

df1 <- data.frame(true = rep(covs.tru.wat, 3), value = c(covs.est.wat.bora, covs.est.wat.mds,  covs.est.wat.visGP), 
                  type = factor(rep(c("BORA-GP","MDS", "visGP"), each = length(covs.tru.wat))))
df1 <- df1[sample.int(nrow(df1)),] #for plotting purposes

p1 = ggplot(df1, aes(x = true, y = value, color = type)) +
  geom_point(alpha = 0.05, size = 3) +
  labs(x = 'Raw variance/covariance values', y = 'Induced variance/covariance values',
       title = 'Comparison of raw and modeled\n covariances on domain-connected points') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = c("BORA-GP" = "blue", "visGP" = "red", "MDS" = "yellow"),
                     name = "Method", 
                     guide = guide_legend(override.aes = list(alpha = 1))) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", size = 1) +
  my_theme
ggsave(filename = 'cov_comp_vars_water.png', p1)

df2 <- data.frame(true = rep(covs.tru.land, 3), value = c(covs.est.land.bora, covs.est.land.mds, covs.est.land.visGP),
                  type = factor(rep(c("BORA-GP", "MDS", "visGP"), each = length(covs.tru.land))))
df2 <- df2[sample.int(nrow(df2)),] #for plotting purposes

p2 = ggplot(df2, aes(x = true, y = value, color = type)) +
  geom_point(alpha = 0.05, size = 3) +
  labs(x = 'Raw (geodesic) variance/covariance values', y = 'Induced variance/covariance values', 
       title = 'Comparison of raw and modeled covariances\non non-domain-connected points') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = c("BORA-GP" = "blue", "visGP" = "red", "MDS" = "yellow"),
                     name = "Method", 
                     guide = guide_legend(override.aes = list(alpha = 1))) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", size = 1) +
  my_theme
ggsave(filename = 'cov_comp_vars_land.png', plot = p2)

png(file='cov_comp_ord.png')
plot(diag(cov.mat.est.bora)[ord], 
     xlab='Nearest-neighbor ordering', 
     ylab='', 
     main='Impact of ordering on variance', 
     cex.main=2.25, 
     cex.lab=2.25, 
     cex.axis=2.25, 
     col='blue')
mtext('BORA-GP-induced marginal variance', side=2, line=2.5, cex=2.25)
dev.off()


