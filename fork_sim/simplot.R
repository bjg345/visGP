library(magrittr)
library(ggplot2)
library(gridExtra)
library(tidyverse)

my_theme <- theme_gray() +
  theme(
    plot.title = element_text(size = 24, hjust = 0.5),
    axis.title = element_text(size = 22),
    strip.text = element_text(size = 12),
    axis.text = element_text(size = 20),
    legend.text=element_text(size=20),
    legend.title=element_text(size=20),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
  
res = readRDS('res.rds')
res$n = as.factor(res$n)

res = res %>%
  gather(Measure, Value, MSE, Coverage, CI.length) %>%
  filter(Method %in% c('BORA GP', 'Euclidean (BRISC)', 'GLGP', 'CS: Maximum precision'))
res = res %>% mutate_all(~ifelse(. == "Euclidean (BRISC)", "Euclidean", .))
res = res %>% mutate_all(~ifelse(. == "CS: Maximum precision", "visGP", .))
res = res %>% mutate_all(~ifelse(. == "BORA GP", "BORA-GP", .))


res = transform(res,
  Measure = factor(Measure, levels=c("MSE","Coverage","CI.length")))

# filter for MSE measure only
res_mse = filter(res, Measure == "MSE")

# create a custom labeller function
custom_labeller <- function(value){
  return(paste0("n = ", value))
}

# define color scheme
color_scheme <- c('BORA-GP' = 'blue', 'visGP' = 'red', 'Euclidean' = 'green', 'GLGP' = 'purple')

# create a plot for all four methods
p_all = ggplot(res_mse, aes(x=factor(Nugget.sd), y=Value, fill=Method)) +
  geom_bar(stat='identity', position='dodge') +
  facet_wrap(~n, scales='free', labeller = as_labeller(custom_labeller)) +
  scale_fill_manual(values = color_scheme) +
  xlab('Nugget standard deviation') +
  ylab('MSE') +
  ggtitle('Simulation results in\nFork-shaped domain') +
  my_theme

# filter for visGP and BORA-GP methods only
res_zoom = filter(res_mse, Method %in% c('BORA-GP', 'visGP'))

# create a "zoomed in" plot for visGP and BORA-GP methods
p_zoom = ggplot(res_zoom, aes(x=factor(Nugget.sd), y=Value, fill=Method)) +
  geom_bar(stat='identity', position='dodge') +
  facet_wrap(~n, scales='free', labeller = as_labeller(custom_labeller)) +
  scale_fill_manual(values = color_scheme) +
  xlab('Nugget standard deviation') +
  ylab('MSE') +
  ggtitle('Simulation results in\nFork-shaped domain (two methods)') +
  my_theme 

# save the plots
ggsave('simplot_all.png', p_all, width = 16, height = 6)
ggsave('simplot_zoom.png', p_zoom, width = 16, height = 6)
