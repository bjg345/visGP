library(sfsmisc)
library(tidyverse)
library(magrittr)
load('../finger_points.Rdata')
z = qnorm(.975)
method.names = c("maxprec", "nearest", "weighted", "neighbor", "euclj", "BRISC")

dat = data.frame(method = c(), n = c(), stdev = c(), quad = c())

ind = 1:1500
getpar = function(j){

        if(j <= 1500){
                    n <- 250
                } else if (j <= 3000){
                    n <- 1200
                  } else (n = 10000)
        
                if((j %% 1500) <= 500){
                    stdev <- .1
                    } else if ((j %% 1500) <= 1000){
                        stdev <- .25
                    } else (stdev <- 1)
                    return(c(n, stdev))
            }
for(method in method.names){
  for(j in ind){
     n=getpar(j)[1]; stdev=getpar(j)[2];  
    set.seed(j)
    
    ind.train = sample.int(10000, .8*n, replace = F)
    noise.train = rnorm(.8*n, sd = stdev)
  
    ind.test = sample.int(10000, .2*n, replace = F)
    
    points = grid.test[ind.test,]
    quad = apply(points, 1, function(p) quadrant(p[1], p[2]))
     
    covered = function(j){
       n=getpar(j)[1]; stdev=getpar(j)[2];   
      pred = readRDS(file.path('maxprec', paste0('pred_', n, '_', stdev, '_', j, '.rds')))
      vals = readRDS(file.path('maxprec', paste0('vals_', n, '_', stdev, '_', j, '.rds')))
      return(sapply(1:length(vals), function(i) vals[i] >= pred[[1,i]] - z*pred[[2,i]] & vals[i] <= pred[[1,i]] + z*pred[[2,i]]))
    }
    
    coverage = covered(j) 
    
    dat = rbind(dat, data.frame(method = method, n = sapply(ind, function(j) getpar(j)[1]), stdev=sapply(ind, function(j) getpar(j)[2]),quad = quad, coverage = coverage))
  }
}


out = dat %>% group_by(n, stdev, method, quad) %>% summarise(coverage = mean(coverage))

print(out, n = nrow(out))
