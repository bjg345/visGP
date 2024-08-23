library(igraph)
library(BRISC)
library(mvtnorm)
if(!dir.exists('visgp_fit_addedge')) dir.create('visgp_fit_addedge')

id = as.numeric(commandArgs(trailingOnly = T))
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

if(file.exists(file.path('visgp_fit_addedge', paste0('fit_', n, '_', stdev, '_', id, '.rds')))) quit()

load('../fork_data.Rdata')
source('../functions.R')

id = as.numeric(commandArgs(trailingOnly = T))
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
set.seed(id)

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

A.train = A.train[ind.train, ind.train]

# test of effect of chordal completion
G = graph_from_adjacency_matrix(A.train, 'undirected')
G.new = is_chordal(G, newgraph=T)$newgraph
edge_diff = ecount(G.new) - ecount(G) 
saveRDS(edge_diff, file.path('visgp_fit_addedge', paste0('edge_diff_', id, '.rds')))
saveRDS(G, file.path('visgp_fit_addedge', paste0('G_', id, '.rds')))
saveRDS(G.new, file.path('visgp_fit_addedge', paste0('G.new_', id, '.rds')))