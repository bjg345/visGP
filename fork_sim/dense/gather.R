library(magrittr)
library(xtable)
library(tidyverse)
library(R.matlab)

# Constants
z = qnorm(.975)
method.names = c("maxprec", "glgp", "bora")

# Function to classify methods
cl <- function(name) {
  switch(name, 
         "maxprec" = "CS: Maximum precision",
         "nearest" = "CS: Nearest clique",
         "weighted" = "CS: Precision-weighted",
         "bora" = "BORA GP",
         "euclid" = "CS: Standard kriging",
         "glgp" = "GLGP",
         "BRISC" = "Euclidean (BRISC)",
         "maxprec_borapoint" = "BORA: Point estimates, max. prec. prediction")
}

# Function to replace NA percentages
replace_na_percent <- function(df) {
  df[df == "NA%"] <- ""
  return(df)
}

# Function to get parameters based on index
getpar <- function(j) {
  n <- ifelse(j <= 1500, 250, ifelse(j <= 3000, 1200, 10000))
  stdev <- ifelse((j %% 1500) <= 500 & (j %% 1500) != 0, 0.1, ifelse((j %% 1500) <= 1000 & (j %% 1500) != 0, 0.25, 1))
  return(c(n, stdev))
}

# Function to compute statistics for each method
compute_stats <- function(method, indices) {
  ext <- if (method == 'glgp') '.mat' else '.rds'
  readFunc <- if (method == 'glgp') function(f) readMat(f)[[1]] else readRDS
  
  results <- data.frame()
  for (j in indices) {
    n <- getpar(j)[1]
    stdev <- getpar(j)[2]
    file_pred <- file.path(method, paste0('pred_', n, '_', stdev, '_', j, ext))
    file_vals <- file.path(method, paste0('vals_', n, '_', stdev, '_', j, ext))
    
    pred <- readFunc(file_pred)
    vals <- readFunc(file_vals)
    
    err2 <- if (method == 'glgp' || method == 'BRISC') {
      (pred - vals) ^ 2
    } else if (method == 'bora') {
      (unlist(pred$pred) - vals) ^ 2
    } else {
      (unlist(pred[1,]) - vals) ^ 2
    }
    
    mse <- mean(err2)
    
    cov <- NA
    len <- NA
    
    if (method == 'glgp') {
      cov <- NA  # Placeholder, as GLGP does not use this calculation
      len <- NA  # Placeholder, as GLGP does not use this calculation
    } else if (method == 'bora' || method == 'bora_visgpfeed') {
      cov <- mean(sapply(1:length(vals), function(i) vals[i] >= pred$ci[[1,i]] & vals[i] <= pred$ci[[2,i]]))
      len <- mean(sapply(1:length(vals), function(i) pred$ci[[2,i]] - pred$ci[[1,i]]))
    } else if (method == 'BRISC') {
      cov <- mean(sapply(1:length(vals), function(i) vals[i] >= pred$prediction.ci[i,1] & vals[i] <= pred$prediction.ci[i,2]))
      len <- mean(sapply(1:length(vals), function(i) pred$prediction.ci[i,2] - pred$prediction.ci[i,1]))
    } else {
      cov <- mean(sapply(1:length(vals), function(i) vals[i] >= pred[[1,i]] - z*pred[[2,i]] & vals[i] <= pred[[1,i]] + z*pred[[2,i]]))
      len <- mean(sapply(1:length(vals), function(i) 2*z*pred[[2,i]]))
    }
    
    method_result <- data.frame(
      method = cl(method), 
      n = n, 
      stdev = stdev, 
      mse = mse, 
      coverage = cov, 
      len = len
    )
    
    results <- rbind(results, method_result)
  }
  return(results)
}


# Main analysis function
analyze_methods <- function(methods, indices) {
  results <- data.frame()
  for (method in methods) {
    method_data <- compute_stats(method, indices)
    results <- rbind(results, method_data)
  }
  
  # Grouping and summarizing the results
  summarized_results <- results %>%
    group_by(n, stdev, method) %>%
    summarise(
      MSE = mean(mse), 
      Coverage = mean(coverage, na.rm = TRUE),
      CI.length = mean(len, na.rm = TRUE)
    ) %>%
    rename(Nugget.sd = stdev, Method = method)
  
  return(summarized_results)
}

# Run analysis
indices1 <- 1:1500
indices2 <- 1501:3000
indices3 <- 3001:4500

out1 <- analyze_methods(method.names, indices1)
out2 <- analyze_methods(method.names[-2], indices2) # Excluding GLGP
out3 <- analyze_methods(method.names[-2], indices3) # Excluding GLGP

res <- replace_na_percent(bind_rows(out1, out2, out3))
print(xtable(res, digits=-3), math.style.exponents=T, include.rownames=F, file='simtab.txt')
saveRDS(res, 'res.rds')
