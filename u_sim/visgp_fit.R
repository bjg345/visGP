library(BRISC)
library(mvtnorm)
if(!dir.exists('visgp_fit')) dir.create('visgp_fit')
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

if(file.exists(file.path('visgp_fit', paste0('fit_', n, '_', stdev, '_', id, '.rds')))) quit()

load('../../u_data.Rdata')
source('../../functions.R')

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

ind.train = sample.int(25000, .8*n, replace = F)
noise.train = rnorm(.8*n, sd = stdev)

ind.test = sample(setdiff(1:25000, ind.train), .2*n, replace = F)
noise.test = rnorm(.2*n, sd = stdev)

grid.train = grid[ind.train,]
grid.test = grid[ind.test,]

ground.train = ground[ind.train]
ground.test = ground[ind.test]

vals.train = ground.train + noise.train
vals.test = ground.test + noise.test

A.train = A[ind.train, ind.train]


if(n > 2500) A.train = A.train * (as.matrix(dist(grid.train)) < 1)
visgp.fit = fit_water(grid.train, A.train, y=vals.train, nu = .5)
saveRDS(visgp.fit, file.path('visgp_fit', paste0('fit_', n, '_', stdev, '_', id, '.rds')))

