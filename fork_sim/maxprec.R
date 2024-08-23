library(BRISC)
library(fields)

if(!dir.exists('maxprec')) dir.create('maxprec')
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

if(file.exists(file.path('maxprec', paste0('vals_', n, '_', stdev, '_', id, '.rds')))) quit()

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
A.test = A.test[ind.test, ind.train]

D.train = as.matrix(dist(grid.train))



maxprec.fit = readRDS(file.path('visgp_fit', paste0('fit_', n, '_', stdev, '_', id, '.rds')))

maxprec.pred = sapply(1:nrow(grid.test), function(i) pred_neighbor(new_loc = grid.test[i,], A=A.train, A_vec = A.test[i, ], 
                                          grid.train, vals.train, maxprec.fit$par, nu=.5, method = 'maxprec', n_neighbors = 10, 
                                           D=D.train, water) )

maxprec.err = unlist(maxprec.pred[1,]) - ground.test

saveRDS(maxprec.fit, file.path('maxprec', paste0('fit_', n, '_', stdev, '_', id, '.rds')))
saveRDS(maxprec.pred, file.path('maxprec', paste0('pred_', n, '_', stdev, '_', id, '.rds')))
saveRDS(maxprec.err, file.path('maxprec', paste0('err_', n, '_', stdev, '_', id, '.rds')))
saveRDS(vals.test, file.path('maxprec', paste0('vals_', n, '_', stdev, '_', id, '.rds')))

