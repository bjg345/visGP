library(tidyverse)
library(mgcv)
library(fields)
library(rgeos)
library(sf)
library(boraGP)
library(INLA)
library(raster)
library(foreach)

library(BRISC)
dir <- 'bora_2'
if(!dir.exists(dir)) dir.create(dir)



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

if(file.exists( file.path(dir, paste0('vals_', n, '_', stdev, '_', id, '.rds')) )) quit()

load('../../fork_data.Rdata')


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

nu <- 0.5
if(n < 2500){
        BRISC.fit = BRISC_estimation(y=y, coords=coords)
        phi = BRISC.fit$Theta[3]
	      tau.sq = BRISC.fit$Theta[2]
	      sigma.sq = sqrt(BRISC.fit$Theta[1])
} else{ # Use defaults to avoid numerical issues of BRISC in large dense sample
    phi <- 10
    tau.sq <- 1
    sigma.sq <- 1
}

starting <- list("phi" = phi, "sigma.sq" = sigma.sq, 
                 "tau.sq" = tau.sq, "nu" = nu)
tuning <- list("phi" = 0.5, "sigma.sq" = 0.5, "tau.sq" = 0.1, "nu" = 0)         # nu is fixed 
min_phi <- 1e-6
priors <- list("phi.Unif" = c(min_phi, 10*phi),                              
               "sigma.sq.IG" = c(2, sigma.sq), "tau.sq.IG" = c(2, tau.sq+.01), 
               "nu.Unif" = c(nu-0.1,nu+0.1))

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
      
      barrier_m.s <- spNNGP(vals.train~ 1, coords = grid.train, starting = starting,
                             method = "response", n.neighbors = m,
                             tuning = tuning, priors = priors, 
                             cov.model = "matern",
                             n.samples = n.samples, n.omp.threads = 1,
                             neighbor.info = barrier_nninfo, verbose = T, return.neighbor.info = (n<1000))
       barrier_p.s <- predict(barrier_m.s, X.0 = matrix(1, nrow = nrow(grid.test), ncol = 1),
                              coords.0 = as.matrix(grid.test),
                              sub.sample = 
                                list(start = burn+1, end = n.samples, thin = 1),
                              nn.indx.0 = barrier_nn.indx.0, 
                             n.omp.threads = 1, verbose = T)
                             
       ystarBRGP <- rowMeans(barrier_p.s$p.y.0)
              yquantBRGP <- apply(barrier_p.s$p.y.0, 1, 
                           function(x) quantile(x, probs = c(0.025, 0.975)))
                           
       saveRDS(barrier_m.s, file.path(dir, paste0('fit_', n, '_', stdev, '_', id, '.rds')))
        saveRDS(list(pred=ystarBRGP, ci = yquantBRGP), file.path(dir, paste0('pred_', n, '_', stdev, '_', id, '.rds')))
          saveRDS(ystarBRGP-vals.test, file.path(dir, paste0('err_', n, '_', stdev, '_', id, '.rds')))
      saveRDS(vals.test, file.path(dir, paste0('vals_', n, '_', stdev, '_', id, '.rds')))


         

