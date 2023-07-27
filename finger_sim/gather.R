library(magrittr)
library(xtable)
library(tidyverse)
library(R.matlab)
z = qnorm(.975)
method.names = c("maxprec", "nearest", "weighted", "euclid", "glgp", "BRISC", "bora")#, "bora_csfeed")

cl = function(name){
	switch(name, 
	"maxprec" = "CS: Maximum precision",
	"nearest" = "CS: Nearest clique",
	"weighted" = "CS: Precision-weighted",
	"bora" = "BORA GP",
	"bora_csfeed" = "BORA: fed parameters",
	"euclid" = "CS: Standard kriging",
	"glgp" = "GLGP",
	"BRISC" = "Euclidean (BRISC)",
	"euclid_approx" = "CS: approximate")
}

replace_na_percent <- function(df) {
  df[df == "NA%"] <- ""
  return(df)
}

out1 = matrix(NA, nrow=length(method.names), ncol = 2)
rownames(out1) = method.names
colnames(out1) = c('mse', 'coverage')

dat = data.frame(method = c(), n = c(), stdev = c(), mse = c(), coverage = c(), len=c())

ind = c(1:1500)

getpar = function(j){

        if(j <= 1500){
                    n <- 250
                } else if (j <= 3000){
                    n <- 1200
                  } else (n = 10000)
        
                if((j %% 1500) <= 500 & (j %% 1500L) > 0L){
                    stdev <- .1
                    } else if ((j %% 1500) <= 1000 & (j %% 1500L) > 0L){
                        stdev <- .25
                    } else (stdev <- 1)
                    return(c(n, stdev))
            }

  for(method in method.names){
    
    if(method == 'glgp') {readFunc = function(f) readMat(f)[[1]]; ext='.mat'}
    else{readFunc = readRDS; ext='.rds'}	    

    if(method == 'glgp'){
      err2.list = lapply(ind, function(j) {n=getpar(j)[1]; stdev=getpar(j)[2]; 
        ( readFunc(file.path(method, paste0('pred_', n, '_', stdev, '_', j, ext))) - 
        readFunc(file.path(method, paste0('vals_', n, '_', stdev, '_', j, ext))) )^2})
    } else if (method == 'BRISC'){
      err2.list = lapply(ind, function(j) {n=getpar(j)[1]; stdev=getpar(j)[2]; 
        ( readFunc(file.path(method, paste0('pred_', n, '_', stdev, '_', j, ext)))$prediction - 
        readFunc(file.path(method, paste0('vals_', n, '_', stdev, '_', j, ext))) )^2})
    } else if(method == 'bora' | method == 'bora_csfeed'){
	err2.list = lapply(ind, function(j) {n=getpar(j)[1]; stdev=getpar(j)[2];
        ( readFunc(file.path(method, paste0('pred_', n, '_', stdev, '_', j, ext)))$pred %>% unlist -
        readFunc(file.path(method, paste0('vals_', n, '_', stdev, '_', j, ext))) )^2})
	}
	else {
      err2.list = lapply(ind, function(j) {n=getpar(j)[1]; stdev=getpar(j)[2]; 
        ( readFunc(file.path(method, paste0('pred_', n, '_', stdev, '_', j, ext)))[1,] %>% unlist - 
        readFunc(file.path(method, paste0('vals_', n, '_', stdev, '_', j, ext))) )^2})
    }
    mse = sapply(err2.list, mean) 
    
    ci.stat = function(j){
       n=getpar(j)[1]; stdev=getpar(j)[2];   
      pred = readFunc(file.path(method, paste0('pred_', n, '_', stdev, '_', j, ext)))
      vals = readFunc(file.path(method, paste0('vals_', n, '_', stdev, '_', j, ext)))
      if(method == 'glgp') {return(NA)}
      else  if(method == 'bora' | method == 'bora_csfeed') return(list(
      cov=sapply(1:length(vals), function(i) vals[i] >= pred$ci[[1,i]]  & vals[i] <= pred$ci[[2,i]]) %>% mean,
      len=sapply(1:length(vals), function(i)  pred$ci[[2,i]] - pred$ci[[1,i]]) %>% mean))
      
	
      else if(method != 'BRISC'){
        return(list(cov = sapply(1:length(vals), function(i) vals[i] >= pred[[1,i]] - z*pred[[2,i]] & vals[i] <= pred[[1,i]] + z*pred[[2,i]]) %>% mean,
        len = sapply(1:length(vals), function(i) 2*z*pred[[2,i]]) %>% mean))
        }
      else{
        return(list(cov=sapply(1:length(vals), function(i) vals[i] >= pred$prediction.ci[i,1] & vals[i] <= pred$prediction.ci[i,2]) %>% mean,
        len=sapply(1:length(vals), function(i)  pred$prediction.ci[i,2]- pred$prediction.ci[i,1]) %>% mean))
      }
    }
    
    if(method == 'glgp'){
      res =  NA
      coverage = NA
      len = NA
    }
    else{
      res =  sapply(ind, ci.stat)
      coverage = res[1,] %>% unlist
      len = res[2,] %>% unlist
    }
    
    
    dat = rbind(dat, data.frame(method = cl(method), n = sapply(ind, function(j) getpar(j)[1]), stdev=sapply(ind, function(j) getpar(j)[2]), mse = mse, coverage = coverage, len = len))
  }

out1 = dat %>% group_by(n, stdev, method) %>% summarise(MSE = mean(mse), Coverage =  (mean(coverage)*100 )%>% round() %>% format() %>% paste0('%'), CI.length = mean(len)) %>% rename(Nugget.sd=stdev, Method = method)
 out1b = dat %>% group_by(n, stdev, method) %>% summarise(MSE = mean(mse), Coverage =  (mean(coverage)), CI.length = mean(len)) %>% rename(Nugget.sd=stdev, Method = method)

#print(xtable(out, digits=-3), math.style.exponents=T)


out2 = matrix(NA, nrow=length(method.names), ncol = 2)
rownames(out2) = method.names
colnames(out2) = c('mse', 'coverage')

dat = data.frame(method = c(), n = c(), stdev = c(), mse = c(), coverage = c(), len=c())

ind = c(1501:3000)

getpar = function(j){

        if(j <= 1500){
                    n <- 250
                } else if (j <= 3000){
                    n <- 1200
                  } else (n = 10000)
        
                if((j %% 1500) <= 500 & (j %% 1500L) > 0L){
                    stdev <- .1
                    } else if ((j %% 1500) <= 1000 & (j %% 1500L) > 0L){
                        stdev <- .25
                    } else (stdev <- 1)
                    return(c(n, stdev))
            }

  for(method in method.names[-5]){
    
    if(method == 'glgp') {readFunc = function(f) readMat(f)[[1]]; ext='.mat'}
    else{readFunc = readRDS; ext='.rds'}	    

    if(method == 'glgp'){
      err2.list = lapply(ind, function(j) {n=getpar(j)[1]; stdev=getpar(j)[2]; 
        ( readFunc(file.path(method, paste0('pred_', n, '_', stdev, '_', j, ext))) - 
        readFunc(file.path(method, paste0('vals_', n, '_', stdev, '_', j, ext))) )^2})
    } else if (method == 'BRISC'){
      err2.list = lapply(ind, function(j) {n=getpar(j)[1]; stdev=getpar(j)[2]; 
        ( readFunc(file.path(method, paste0('pred_', n, '_', stdev, '_', j, ext)))$prediction - 
        readFunc(file.path(method, paste0('vals_', n, '_', stdev, '_', j, ext))) )^2})
    } else if(method == 'bora' | method == 'bora_csfeed'){
	err2.list = lapply(ind, function(j) {n=getpar(j)[1]; stdev=getpar(j)[2];
        ( readFunc(file.path(method, paste0('pred_', n, '_', stdev, '_', j, ext)))$pred %>% unlist -
        readFunc(file.path(method, paste0('vals_', n, '_', stdev, '_', j, ext))) )^2})
	}
	else {
      err2.list = lapply(ind, function(j) {n=getpar(j)[1]; stdev=getpar(j)[2]; 
        ( readFunc(file.path(method, paste0('pred_', n, '_', stdev, '_', j, ext)))[1,] %>% unlist - 
        readFunc(file.path(method, paste0('vals_', n, '_', stdev, '_', j, ext))) )^2})
    }
    mse = sapply(err2.list, mean) 
    
    ci.stat = function(j){
       n=getpar(j)[1]; stdev=getpar(j)[2];   
      pred = readFunc(file.path(method, paste0('pred_', n, '_', stdev, '_', j, ext)))
      vals = readFunc(file.path(method, paste0('vals_', n, '_', stdev, '_', j, ext)))
      if(method == 'glgp') {return(NA)}
      else if(method == 'bora' | method == 'bora_csfeed') return(list(
      cov=sapply(1:length(vals), function(i) vals[i] >= pred$ci[[1,i]]  & vals[i] <= pred$ci[[2,i]]) %>% mean,
      len=sapply(1:length(vals), function(i)  pred$ci[[2,i]] - pred$ci[[1,i]]) %>% mean))
      
	
      else if(method != 'BRISC'){
        return(list(cov = sapply(1:length(vals), function(i) vals[i] >= pred[[1,i]] - z*pred[[2,i]] & vals[i] <= pred[[1,i]] + z*pred[[2,i]]) %>% mean,
        len = sapply(1:length(vals), function(i) 2*z*pred[[2,i]]) %>% mean))
        }
      else{
        return(list(cov=sapply(1:length(vals), function(i) vals[i] >= pred$prediction.ci[i,1] & vals[i] <= pred$prediction.ci[i,2]) %>% mean,
        len=sapply(1:length(vals), function(i)  pred$prediction.ci[i,2]- pred$prediction.ci[i,1]) %>% mean))
      }
    }
    
    if(method == 'glgp'){
      res =  NA
      coverage = NA
      len = NA
    }
    else{
      res =  sapply(ind, ci.stat)
      coverage = res[1,] %>% unlist
      len = res[2,] %>% unlist
    }
    
    
    dat = rbind(dat, data.frame(method = cl(method), n = sapply(ind, function(j) getpar(j)[1]), stdev=sapply(ind, function(j) getpar(j)[2]), mse = mse, coverage = coverage, len = len))
  }

out2 = dat %>% group_by(n, stdev, method) %>% summarise(MSE = mean(mse), Coverage =  (mean(coverage)*100 )%>% round() %>% format() %>% paste0('%'), CI.length = mean(len)) %>% rename(Nugget.sd=stdev, Method = method)
out2b = dat %>% group_by(n, stdev, method) %>% summarise(MSE = mean(mse), Coverage =  (mean(coverage) ), CI.length = mean(len)) %>% rename(Nugget.sd=stdev, Method = method)


out3 = matrix(NA, nrow=length(method.names), ncol = 2)
rownames(out3) = method.names
colnames(out3) = c('mse', 'coverage')

dat = data.frame(method = c(), n = c(), stdev = c(), mse = c(), coverage = c(), len=c())
 
ind = c(3001:4500)

getpar = function(j){

        if(j <= 1500){
                    n <- 250
                } else if (j <= 3000){
                    n <- 1200
                  } else (n = 10000)
        
                if((j %% 1500) <= 500 & (j %% 1500L) > 0L){
                    stdev <- .1
                    } else if ((j %% 1500) <= 1000 & (j %% 1500L) > 0L){
                        stdev <- .25
                    } else (stdev <- 1)
                    return(c(n, stdev))
            }

  for(method in method.names[-5]){
    
    if(method == 'glgp') {readFunc = function(f) readMat(f)[[1]]; ext='.mat'}
    else{readFunc = readRDS; ext='.rds'}	    

    if(method == 'glgp'){
      err2.list = lapply(ind, function(j) {n=getpar(j)[1]; stdev=getpar(j)[2]; 
        ( readFunc(file.path(method, paste0('pred_', n, '_', stdev, '_', j, ext))) - 
        readFunc(file.path(method, paste0('vals_', n, '_', stdev, '_', j, ext))) )^2})
    } else if (method == 'BRISC'){
      err2.list = lapply(ind, function(j) {n=getpar(j)[1]; stdev=getpar(j)[2]; 
        ( readFunc(file.path(method, paste0('pred_', n, '_', stdev, '_', j, ext)))$prediction - 
        readFunc(file.path(method, paste0('vals_', n, '_', stdev, '_', j, ext))) )^2})
    } else if(method == 'bora' | method == 'bora_csfeed'){
	err2.list = lapply(ind, function(j) {n=getpar(j)[1]; stdev=getpar(j)[2];
        ( readFunc(file.path(method, paste0('pred_', n, '_', stdev, '_', j, ext)))$pred %>% unlist -
        readFunc(file.path(method, paste0('vals_', n, '_', stdev, '_', j, ext))) )^2})
	}
	else {
      err2.list = lapply(ind, function(j) {n=getpar(j)[1]; stdev=getpar(j)[2]; 
        ( readFunc(file.path(method, paste0('pred_', n, '_', stdev, '_', j, ext)))[1,] %>% unlist - 
        readFunc(file.path(method, paste0('vals_', n, '_', stdev, '_', j, ext))) )^2})
    }
    mse = sapply(err2.list, mean) 
    
    ci.stat = function(j){
       n=getpar(j)[1]; stdev=getpar(j)[2];   
      pred = readFunc(file.path(method, paste0('pred_', n, '_', stdev, '_', j, ext)))
      vals = readFunc(file.path(method, paste0('vals_', n, '_', stdev, '_', j, ext)))
      if(method == 'glgp') {return(NA)}
      else if(method == 'bora' | method == 'bora_csfeed') return(list(
      cov=sapply(1:length(vals), function(i) vals[i] >= pred$ci[[1,i]]  & vals[i] <= pred$ci[[2,i]]) %>% mean,
      len=sapply(1:length(vals), function(i)  pred$ci[[2,i]] - pred$ci[[1,i]]) %>% mean))
      
	
      else if(method != 'BRISC'){
        return(list(cov = sapply(1:length(vals), function(i) vals[i] >= pred[[1,i]] - z*pred[[2,i]] & vals[i] <= pred[[1,i]] + z*pred[[2,i]]) %>% mean,
        len = sapply(1:length(vals), function(i) 2*z*pred[[2,i]]) %>% mean))
        }
      else{
        return(list(cov=sapply(1:length(vals), function(i) vals[i] >= pred$prediction.ci[i,1] & vals[i] <= pred$prediction.ci[i,2]) %>% mean,
        len=sapply(1:length(vals), function(i)  pred$prediction.ci[i,2]- pred$prediction.ci[i,1]) %>% mean))
      }
    }
    
    if(method == 'glgp'){
      res =  NA
      coverage = NA
      len = NA
    }
    else{
      res =  sapply(ind, ci.stat)
      coverage = res[1,] %>% unlist
      len = res[2,] %>% unlist
    }
    
    
    dat = rbind(dat, data.frame(method = cl(method), n = sapply(ind, function(j) getpar(j)[1]), stdev=sapply(ind, function(j) getpar(j)[2]), mse = mse, coverage = coverage, len = len))
  }

out3 = dat %>% group_by(n, stdev, method) %>% summarise(MSE = mean(mse), Coverage =  (mean(coverage)*100 )%>% round() %>% format() %>% paste0('%'), CI.length = mean(len)) %>% rename(Nugget.sd=stdev, Method = method)
out3b = dat %>% group_by(n, stdev, method) %>% summarise(MSE = mean(mse), Coverage =  (mean(coverage) ), CI.length = mean(len)) %>% rename(Nugget.sd=stdev, Method = method)

res = replace_na_percent(rbind(out1, out2, out3))
print(xtable(res, digits=-3), math.style.exponents=T, include.rownames=F, file='simtab.txt')
saveRDS(rbind(out1b, out2b, out3b), 'res.rds')

