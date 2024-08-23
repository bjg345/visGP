library(ggplot2)
library(viridis)

load('fork_data.Rdata')

my_theme <- theme_gray() +
  theme(
    plot.title = element_text(size = 24, hjust = 0.5),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 20),
    legend.title=element_text(size=20),
    legend.text=element_text(size=20)
  )
  
x=water@polygons$water@Polygons[[1]]@coords[,1]
y=water@polygons$water@Polygons[[1]]@coords[,2]

values=ground
group = rep(NA, nrow(grid))

group[test.id] = 'test'
group[train.id] = 'train'

used.id = c(test.id[1:10000], train.id[1:10000])

group = group[used.id]
grid = grid[used.id,]
values = values[used.id]

p1=ggplot(mapping=aes(x=x,y=y))+
  geom_point(mapping=aes(x=grid[,1], y=grid[,2], color=values))+
  geom_polygon(fill=NA, colour='black') + scale_color_viridis()+
 ggtitle('Underlying function values\nin fork-shaped domain') +
 my_theme
ggsave('fork_vals.png', p1)


p2=ggplot(mapping=aes(x=x,y=y)) + 
  geom_point(mapping=aes(x=grid[,1], y=grid[,2], color=group))+geom_polygon(fill=NA, colour='black')+
  ggtitle('Test/train groups in\nfork-shaped domain') +
  my_theme
ggsave('fork_groups.png', p2)
