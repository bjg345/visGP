library(R.matlab)
library(xtable)

# get simulation results 
res <- readRDS('res.rds')

# make table

tab <- data.frame('MSE' = rep(NA, 3), 'Coverage' = rep(NA, 3), 'CI length' = rep(NA, 3))

rownames(tab) <- c('BRISC (Euclidean)', 'visGP (maximum precision)', 'BORA-GP')

for(i in 1:3){
    tab[i, 1] = mean(res$mse_results[,i])
    tab[i, 2] = mean(res$coverage_results[,i])
    tab[i, 3] = mean(res$length_results[,i])
}

formatted_tab <- data.frame(
  Col1 = formatC(tab[, 1], format = "e", digits = 2),
  Col2 = paste0(round(tab[, 2] * 100, 1), "%"),
  Col3 = round(tab[, 3], 3)
)

rownames(formatted_tab) <- rownames(tab); colnames(formatted_tab) <- colnames(tab)

xtable_obj <- xtable(formatted_tab)

print(xtable_obj, file = "tab.txt")
