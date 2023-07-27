library(BRISC)

if(!dir.exists('neighbor')) dir.create('neighbor')

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

#if(file.exists(file.path('neighbor', paste0('vals_', n, '_', stdev, '_', id, '.rds')))) quit()

load('../finger_data.Rdata')
source('../functions.R')

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

A.train = A.train[ind.train, ind.train]

neighbor.fit = fit.water(grid.train, A.train, y=vals.train, n.neighbors  = 10, nu = .5, method='neighbor')
 neighbor.pred = sapply(1:nrow(grid.test), function(i) pred_neighbor(new.loc = grid.test[i,], A= A.train, A.vec = NULL, 
                                                                    grid.train, vals.train, neighbor.fit$fit$par, method='euclidean',
                                                                    nu=.5,n.neighbors=10, D=as.matrix(dist(grid.train)), water) )
 neighbor.err = unlist(neighbor.pred[1,]) - ground.test

saveRDS(neighbor.fit, file.path('neighbor', paste0('fit_', n, '_', stdev, '_', id, '.rds')))
saveRDS(neighbor.pred, file.path('neighbor', paste0('pred_', n, '_', stdev, '_', id, '.rds')))
saveRDS(neighbor.err, file.path('neighbor', paste0('err_', n, '_', stdev, '_', id, '.rds')))
saveRDS(vals.test, file.path('neighbor', paste0('vals_', n, '_', stdev, '_', id, '.rds')))

