# Compare the log-likelihoods on visGP model achieved by boraGP posterior means vs visGP point estimates

source('../../functions.R')

library(tidyverse)
library(mgcv)
library(fields)
library(rgeos)
library(sf)
library(boraGP)
library(INLA)
library(raster)
library(foreach)
library(geostatsp)
library(igraph)
library(gRbase)

id = as.numeric(commandArgs(trailingOnly = T))
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

out_dir = 'visgp_llcomp'
if(!dir.exists(out_dir)) dir.create(out_dir)

load('../../u_data.Rdata')


id = as.numeric(commandArgs(trailingOnly = T))
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


ind.train = sample.int(10000, .8*n, replace = F)
noise.train = rnorm(.8*n, sd = stdev)

ind.test = sample.int(10000, .2*n, replace = F)
noise.test = rnorm(.2*n, sd = stdev)

grid.train = grid.train[ind.train,]
grid.test = grid.test[ind.test,]

ground.train = ground.train[ind.train]
ground.test = ground.test[ind.test]

vals.train = ground.train + noise.train
vals.test = ground.test + noise.test


nu <- 0.5

A.train = A.train[ind.train, ind.train]

if(n > 2500) A.train = A.train * (as.matrix(dist(grid.train)) < 1)

# function to get loglikelihood from parameters and observations
get_loglik = function(params, locations, vals, adj){
  
  nu = 0.5
  
  D = as.matrix(dist(locations))
 
  g.unchord = graph_from_adjacency_matrix(adj, 'undirected')
  V(g.unchord)$name = 1:(length(vals))
  g.chord = is_chordal(g.unchord, newgraph=T)$newgraph
  
  rip = rip(g.chord)

  cliques = sapply(rip$cliques, as.numeric)
  seps = sapply(rip$separators, as.numeric)
      
  ll = cov.select.loglik(params, nu = nu, y = vals, dists = D, cliques = cliques,
      seps = seps, coords = locations)
      
  return(-1*ll)
  
}


params_visgp = readRDS(file.path('../visgp_fit', paste0('fit_', n, '_', stdev, '_', id, '.rds')))$par       
params_bora = readRDS(file.path('../maxprec_borapoint', paste0('par_', n, '_', stdev, '_', id, '.rds')))

ll_visgp = get_loglik(params_visgp, grid.train, vals.train, A.train)
ll_bora = get_loglik(params_bora, grid.train, vals.train, A.train)

res = c(ll_visgp, ll_bora); names(res) = c('ll_visgp', 'll_bora')

saveRDS(res, file.path(out_dir, paste0('res_', id, '.rds')))
