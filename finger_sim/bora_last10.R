library(tidyverse)
library(mgcv)
library(fields)
library(dplyr)
library(rgeos)
library(sf)
library(boraGP)
library(INLA)
library(raster)
library(foreach)

library(BRISC)

if(!dir.exists('bora_last10')) dir.create('bora_last10')


load('../finger_data.Rdata')

t0=Sys.time()
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

n.samples = 10000
burn = 5000

ind.train = sample.int(10000, .8*n, replace = F)
noise.train = rnorm(.8*n, sd = stdev)

ind.test = sample.int(10000, .2*n, replace = F)
noise.test = rnorm(.2*n, sd = stdev)

grid.train = grid.train[ind.train,]
grid.test = grid.test[ind.test,]

pred = readRDS(file.path('bora', paste0('pred_', n, '_', stdev, '_', id, '.rds')))
vals = readRDS(file.path('bora', paste0('vals_', n, '_', stdev, '_', id, '.rds')))

if(id %in% c(5, 6, 177, 185, 270, 278, 397, 481, 498,525, 554, 
  699, 829,  845, 853, 921, 958, 1009, 1083, 1114, 1216,
 1258, 1311, 1346, 1420, 2369, 2370, 2627)){
   ord = order(grid.train[,2])
} else if(id %in% c(17, 249, 519, 677, 837, 998, 1002, 1067, 1215)){
  ord = order(grid.train[,1]+grid.train[,2])
} else{
  ord = order(grid.train[,1])
}

ind.eval = which(ord > .9*length(ord))

vals = ground.train[ind.train]

fit = readRDS(file.path('bora', paste0('fit_', n, '_', stdev, '_', id, '.rds')))


cov = sapply(ind.eval, function(i) between(vals[i], fit$y.hat.quant[i,2], fit$y.hat.quant[i,3]))
saveRDS(cov, file.path('bora_last10', paste0('cov_', n, '_', stdev, '_', id, '.rds')))

         
print(Sys.time()-t0)

