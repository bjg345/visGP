library(ggplot2)

inds = 1501:2000

files = sapply(inds, function(i) list.files(path = 'visgp_llcomp', pattern = paste0('res_', i, '.rds'), full.names = T))
params = sapply(files, readRDS)

ll_visGP = unlist(params[1,])
ll_bora = unlist(params[2,])

ll_plot = qplot(ll_visGP, ll_bora)
ll_plot = ll_plot + ggtitle('visGP-based log-likelihood comparison for U-shaped\n, n=1200, nugget.sd=.1') +
labs(x = 'LL from visGP estimates', y = 'LL from boraGP estimates') +
geom_abline(slope = 1, intercept = 0, color = "blue", linetype = "dotted")
ggsave('ll_comp_visgp.png', ll_plot)
