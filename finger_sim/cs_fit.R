library(BRISC)
library(mvtnorm)
if(!dir.exists('cs_fit')) dir.create('cs_fit')
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

#if(file.exists(file.path('cs_fit', paste0('fit_', n, '_', stdev, '_', id, '.rds')))) quit()

load('../finger_data.Rdata')
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



if(n > 25000) k = 10 else k = NULL
if(n > 25000) n.cliques = 100 else n.cliques = NULL
if(n > 2500) A.train = A.train * (as.matrix(dist(grid.train)) < 1)
t0=Sys.time()
cs.fit = fit.water(grid.train, A.train, y=vals.train, nu = .5, method='cov.select', n.neighbors=k, n.cliques=n.cliques)
Sys.time()-t0
saveRDS(cs.fit, file.path('cs_fit', paste0('fit_', n, '_', stdev, '_', id, 'test', '.rds')))

