library(ggplot2)

inds = 1501:2000

f_bora = sapply(inds, function(i) list.files(path = '../maxprec_borapoint', pattern = paste0('par.*_', i, '.rds'), full.names = T))
params_bora = sapply(f_bora, readRDS)

f_visGP = sapply(inds, function(i) list.files(path = '../visgp_fit', pattern = paste0('fit.*_', i, '.rds'), full.names = T))
params_visGP = sapply(f_visGP, function(f) readRDS(f)$par)

range_bora = exp(params_bora[1,])
nugget_bora = exp(params_bora[2,])
var_bora = exp(2*params_bora[3,])

range_visGP = exp(params_visGP[1,])
nugget_visGP = exp(params_visGP[2,])
var_visGP = exp(2*params_visGP[3,])

range_plot = qplot(range_visGP, range_bora)
range_plot = range_plot + ggtitle('Range comparison for U-shaped, n=1200, nugget.sd=.1')
ggsave('range_comp.png', range_plot)
var_plot = qplot(var_visGP, var_bora)
var_plot = var_plot + ggtitle('Variance comparison for U-shaped, n=1200, nugget.sd=.1')
ggsave('sigma2_comp.png', var_plot)
nug_plot = qplot(nugget_visGP, nugget_bora, alpha = .25)
nug_plot = nug_plot + ggtitle('Nugget comparison for U-shaped, n=1200, nugget.sd=.1') +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "blue") +
  xlim(0, .015) + ylim(0, .015)
ggsave('nugget_comp.png', nug_plot)
rat_plot = qplot(var_visGP/range_visGP, var_bora/range_bora, alpha = .25)
rat_plot = rat_plot + ggtitle('Variance/range comparison for U-shaped, n=1200, nugget.sd=.1') +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "blue") +
  xlim(0, .007) + ylim(0, .016)
ggsave('ratio_comp.png', rat_plot)