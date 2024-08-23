# Compare the log-likelihoods on boraGP model achieved by boraGP posterior means vs visGP point estimates

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

out_dir = 'bora_llcomp'
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

vals.train = ground.train + noise.train
vals.test = ground.test + noise.test


nu <- 0.5

bound=as(water, 'SpatialLines')

m = 10
ord = order(grid.train[,1])
barrier_nninfo_all <- barrier_neighbor(coords = grid.train, coords.0 =  grid.test, ord = ord,
                                                n.neighbors = m,
                                                barrier = st_as_sf(bound),
                                                cores = 1,
                                                verbose = T, min.gridsize=.01,
                                                debug = list(barrier_n.indx = NULL,
                                                            barrier_dist = NULL,
                                                             barrier_nn.indx.0_list = NULL,
                                                             barrier_dist0 = NULL,
                                                            ref_window = NULL,
                                                             nonref_window = NULL,
                                                             ref_fill = TRUE,
                                                             nonref_fill = TRUE))


barrier_nn.indx.0_list <- barrier_nninfo_all$barrier_nn.indx.0_list
       barrier_nn.indx.0 <- do.call("rbind", barrier_nn.indx.0_list)
       barrier_n.indx <- barrier_nninfo_all$barrier_n.indx
       barrier_nn.indx <- as.integer(unlist(barrier_n.indx)[-1]-1)
       barrier_nn.indx.lu <- c(cumsum(c(0,0,sapply(barrier_n.indx, length)[-1]))[1:length(ind.train)],
                               c(0,sapply(barrier_n.indx, length)[-1])) %>% as.integer()
       barrier_nninfo <- list(type = "barrier",
                              n.indx = barrier_n.indx,
                              n.neighbors = m, nn.indx = barrier_nn.indx,
                             nn.indx.lu = barrier_nn.indx.lu, ord = ord)
 
# function to get loglikelihood from parameters and observations
get_loglik = function(params, locations, vals, ord){
  
  nu = 0.5
  range = exp(params[1]) %>% as.numeric()
  nugget = exp(params[2]) %>% as.numeric()
  sigma = exp(params[3]) %>% as.numeric()
  mean = params[4] %>% as.numeric()
  
  b = list()
  b[[1]] = 0
  f = list()
  f[[1]] = nugget + sigma^2
  ll = list() # log likelihoods
  
  N = length(ord)
  coords.ord = grid.train[ord,]
  D.ord = as.matrix(dist(coords.ord))
  y.ord = vals[ord] - mean
  
  ll[[1]] = dnorm(y.ord[1], mean = 0, sd = sqrt(f[[1]]), log = T)
  
  for(i in 2:nrow(grid.train)){
    ind = barrier_nninfo_all$barrier_n.indx[[i]]
    M = length(ind)
    dists = D.ord[ind, ind]

    CC = geostatsp::matern(dists, param=c(range=range, variance=sigma^2, shape=nu)) + nugget*diag(M)
    C =  geostatsp::matern(D.ord[i, ind], param=c(range=range, variance=sigma^2, shape=nu))
    b[[i]] = solve(CC) %*% C
    f[[i]] = nugget+sigma^2-sum(b[[i]] * C)
    
    y_neighbors = y.ord[ind]
   
    mu = sum(b[[i]]*y_neighbors)
    ll[[i]] = dnorm(y.ord[i], mean = mu, sd = sqrt(f[[i]]), log = T)
  }

  return(sum(unlist(ll)))
  
}


params_visgp = readRDS(file.path('../visgp_fit', paste0('fit_', n, '_', stdev, '_', id, '.rds')))$par       
params_bora = readRDS(file.path('../maxprec_borapoint', paste0('par_', n, '_', stdev, '_', id, '.rds')))

ll_visgp = get_loglik(params_visgp, grid.train, vals.train, ord)
ll_bora = get_loglik(params_bora, grid.train, vals.train, ord)

res = c(ll_visgp, ll_bora); names(res) = c('ll_visgp', 'll_bora')

saveRDS(res, file.path(out_dir, paste0('res_', id, '.rds')))
