library(ggplot2)
library(viridis)

load('../finger_data.Rdata')

x=water@polygons$water@Polygons[[1]]@coords[,1]
y=water@polygons$water@Polygons[[1]]@coords[,2]
values=ground
group = rep(NA, nrow(grid))
group[test.id] = 'test'
group[train.id] = 'train'

p1=ggplot(mapping=aes(x=x,y=y))+
  geom_point(mapping=aes(x=grid[,1], y=grid[,2], color=values))+
  geom_polygon(fill=NA, colour='black') + scale_color_viridis()+
 ggtitle('Underlying function values in fork-shaped domain')
ggsave('finger_vals.png', p1)


p2=ggplot(mapping=aes(x=x,y=y)) + 
  geom_point(mapping=aes(x=grid[,1], y=grid[,2], color=group))+geom_polygon(fill=NA, colour='black')+
  ggtitle('Test/train groups in fork-shaped domain')
ggsave('finger_groups.png', p2)

rm(list=ls())

load('../u_data.Rdata')

x=water@polygons$water@Polygons[[1]]@coords[,1]
y=water@polygons$water@Polygons[[1]]@coords[,2]
values=ground
group = rep(NA, nrow(grid))
group[test.id] = 'test'
group[train.id] = 'train'

p3=ggplot(mapping=aes(x=x,y=y))+
  geom_point(mapping=aes(x=grid[,1], y=grid[,2], color=values))+
  geom_polygon(fill=NA, colour='black') + scale_color_viridis()+
  ggtitle('Underlying function values in U-shaped domain')
ggsave('u_vals.png', p3)

p4=ggplot(mapping=aes(x=x,y=y))+
  geom_point(mapping=aes(x=grid[,1], y=grid[,2], color=group))+
  geom_polygon(fill=NA, colour='black') + 
  ggtitle('Test/train groups in U-shaped domain')
ggsave('u_groups.png', p4)
