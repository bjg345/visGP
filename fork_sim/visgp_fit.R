library(BRISC)
library(mvtnorm)
library(igraph)
library(gRbase)

if(!dir.exists('visgp_fit')) dir.create('visgp_fit')

load('../fork_data.Rdata')
source('../functions.R')

id = as.numeric(commandArgs(trailingOnly = TRUE))
if(id <= 1500){
  n = 250
} else if (id <= 3000){
  n = 1200
} else {
  n = 10000
}

if((id %% 1500) <= 500 & (id %% 1500L) > 0L){
  stdev = 0.1
} else if ((id %% 1500) <= 1000 & (id %% 1500L) > 0L){
  stdev = 0.25
} else {
  stdev = 1
}

set.seed(id)

ind.train = sample.int(10000, 0.8*n, replace = FALSE)
noise.train = rnorm(0.8*n, sd = stdev)

ind.test = sample.int(10000, 0.2*n, replace = FALSE)
noise.test = rnorm(0.2*n, sd = stdev)

grid.train = grid.train[ind.train,]
grid.test = grid.test[ind.test,]

ground.train = ground.train[ind.train]
ground.test = ground.test[ind.test]

vals.train = ground.train + noise.train
vals.test = ground.test + noise.test

A.train = A.train[ind.train, ind.train]

if(n > 2500) A.train = A.train * (as.matrix(dist(grid.train)) < 1)

visgp.fit = fit_water(grid.train, A.train, y=vals.train, nu = 0.5)
saveRDS(visgp.fit, file.path('visgp_fit', paste0('fit_', n, '_', stdev, '_', id, '.rds')))