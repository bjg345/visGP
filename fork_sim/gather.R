library(magrittr)
library(xtable)
library(tidyverse)
library(R.matlab)

z = qnorm(.975)
method.names = c("maxprec", "nearest", "weighted", "euclid", "glgp", "BRISC", "bora")

cl = function(name) {
  switch(name, 
         "maxprec" = "CS: Maximum precision",
         "nearest" = "CS: Nearest clique",
         "weighted" = "CS: Precision-weighted",
         "bora" = "BORA GP",
         "euclid" = "CS: Standard kriging",
         "glgp" = "GLGP",
         "BRISC" = "Euclidean (BRISC)",
         "maxprec_borapoint" = "BORA: Point estimates, max. prec. prediction"
  )
}

replace_na_percent <- function(df) {
  df[df == "NA%"] <- ""
  return(df)
}

getpar = function(j) {
  if (j <= 1500) {
    n <- 250
  } else if (j <= 3000) {
    n <- 1200
  } else {
    n <- 10000
  }
  
  if ((j %% 1500) <= 500 & (j %% 1500L) > 0L) {
    stdev <- 0.1
  } else if ((j %% 1500) <= 1000 & (j %% 1500L) > 0L) {
    stdev <- 0.25
  } else {
    stdev <- 1
  }
  
  return(c(n, stdev))
}

process_data = function(ind, include_glgp = TRUE) {
  dat = data.frame(method = c(), n = c(), stdev = c(), mse = c(), coverage = c(), len = c())
  
  methods_to_process = if (include_glgp) method.names else method.names[-5]
  
  for (method in methods_to_process) {
    if (method == 'glgp') {
      readFunc = function(f) readMat(f)[[1]]
      ext = '.mat'
    } else {
      readFunc = readRDS
      ext = '.rds'
    }
    
    err2.list = lapply(ind, function(j) {
      n = getpar(j)[1]
      stdev = getpar(j)[2]
      
      if (method == 'glgp') {
        (readFunc(file.path(method, paste0('pred_', n, '_', stdev, '_', j, ext))) - 
           readFunc(file.path(method, paste0('vals_', n, '_', stdev, '_', j, ext))))^2
      } else if (method == 'BRISC') {
        (readFunc(file.path(method, paste0('pred_', n, '_', stdev, '_', j, ext)))$prediction - 
           readFunc(file.path(method, paste0('vals_', n, '_', stdev, '_', j, ext))))^2
      } else if (method == 'bora' | method == 'bora_visgpfeed') {
        (readFunc(file.path(method, paste0('pred_', n, '_', stdev, '_', j, ext)))$pred %>% unlist -
           readFunc(file.path(method, paste0('vals_', n, '_', stdev, '_', j, ext))))^2
      } else {
        (readFunc(file.path(method, paste0('pred_', n, '_', stdev, '_', j, ext)))[1,] %>% unlist - 
           readFunc(file.path(method, paste0('vals_', n, '_', stdev, '_', j, ext))))^2
      }
    })
    
    mse = sapply(err2.list, mean)
    
    ci.stat = function(j) {
      n = getpar(j)[1]
      stdev = getpar(j)[2]
      pred = readFunc(file.path(method, paste0('pred_', n, '_', stdev, '_', j, ext)))
      vals = readFunc(file.path(method, paste0('vals_', n, '_', stdev, '_', j, ext)))
      
      if (method == 'glgp') {
        return(NA)
      } else if (method == 'bora' | method == 'bora_visgpfeed') {
        return(list(
          cov = sapply(1:length(vals), function(i) vals[i] >= pred$ci[[1,i]] & vals[i] <= pred$ci[[2,i]]) %>% mean,
          len = sapply(1:length(vals), function(i) pred$ci[[2,i]] - pred$ci[[1,i]]) %>% mean
        ))
      } else if (method != 'BRISC') {
        return(list(
          cov = sapply(1:length(vals), function(i) vals[i] >= pred[[1,i]] - z*pred[[2,i]] & vals[i] <= pred[[1,i]] + z*pred[[2,i]]) %>% mean,
          len = sapply(1:length(vals), function(i) 2*z*pred[[2,i]]) %>% mean
        ))
      } else {
        return(list(
          cov = sapply(1:length(vals), function(i) vals[i] >= pred$prediction.ci[i,1] & vals[i] <= pred$prediction.ci[i,2]) %>% mean,
          len = sapply(1:length(vals), function(i) pred$prediction.ci[i,2] - pred$prediction.ci[i,1]) %>% mean
        ))
      }
    }
    
    if (method == 'glgp') {
      res = NA
      coverage = NA
      len = NA
    } else {
      res = sapply(ind, ci.stat)
      coverage = res[1,] %>% unlist
      len = res[2,] %>% unlist
    }
    
    dat = rbind(dat, data.frame(
      method = cl(method),
      n = sapply(ind, function(j) getpar(j)[1]),
      stdev = sapply(ind, function(j) getpar(j)[2]),
      mse = mse,
      coverage = coverage,
      len = len
    ))
  }
  
  out = dat %>% 
    group_by(n, stdev, method) %>% 
    summarise(
      MSE = mean(mse),
      Coverage = (mean(coverage)*100) %>% round() %>% format() %>% paste0('%'),
      CI.length = mean(len)
    ) %>% 
    rename(Nugget.sd = stdev, Method = method)
  
  outb = dat %>% 
    group_by(n, stdev, method) %>% 
    summarise(
      MSE = mean(mse),
      Coverage = mean(coverage),
      CI.length = mean(len)
    ) %>% 
    rename(Nugget.sd = stdev, Method = method)
  
  return(list(out = out, outb = outb))
}

# Process data for different ranges
result1 = process_data(1:1500)
result2 = process_data(1501:3000, include_glgp = FALSE)
result3 = process_data(3001:4500, include_glgp = FALSE)

# Combine results
res = replace_na_percent(rbind(result1$out, result2$out, result3$out))
res_b = rbind(result1$outb, result2$outb, result3$outb)

# Output results
print(xtable(res, digits=-3), math.style.exponents=T, include.rownames=F, file='simtab.txt')
saveRDS(res_b, 'res.rds')