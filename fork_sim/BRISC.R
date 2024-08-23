library(BRISC)

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


if(!dir.exists('BRISC')) dir.create('BRISC')
if(file.exists( file.path('BRISC', paste0('fit_', n, '_', stdev, '_', id, '.rds')))) quit()
rm(list=ls())
load('../fork_data.Rdata')

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

BRISC.fit = BRISC_estimation(grid.train, y = vals.train, cov.model = 'exponential', n.neighbors=10)
BRISC.pred = BRISC_prediction(BRISC.fit, grid.test)
BRISC.err = BRISC.pred$prediction - ground.test

saveRDS(BRISC.fit, file.path('BRISC', paste0('fit_', n, '_', stdev, '_', id, '.rds')))
saveRDS(BRISC.pred, file.path('BRISC', paste0('pred_', n, '_', stdev, '_', id, '.rds')))
saveRDS(BRISC.err, file.path('BRISC', paste0('err_', n, '_', stdev, '_', id, '.rds')))
saveRDS(vals.test, file.path('BRISC', paste0('vals_', n, '_', stdev, '_', id, '.rds')))
