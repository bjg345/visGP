load('u_data.Rdata')

write.csv(grid.train, file.path('u_sim', 'grid.train.csv'), row.names = F)
write.csv(grid.test, file.path('u_sim', 'grid.test.csv'), row.names = F)
#write.csv(vals.test, file.path('u_sim', 'vals.test.csv'), row.names = F)
write.csv(ground.test, file.path('u_sim', 'ground.test.csv'), row.names = F)
#write.csv(vals.train, file.path('u_sim', 'vals.train.csv'), row.names = F)
write.csv(ground.train, file.path('u_sim', 'ground.train.csv'), row.names = F)

load('finger_data.Rdata')
write.csv(grid.train, file.path('finger_sim', 'grid.train.csv'), row.names = F)
write.csv(grid.test, file.path('finger_sim', 'grid.test.csv'), row.names = F)
#write.csv(vals.test, file.path('finger_sim', 'vals.test.csv'), row.names = F)
write.csv(ground.test, file.path('finger_sim', 'ground.test.csv'), row.names = F)
#write.csv(vals.train, file.path('finger_sim', 'vals.train.csv'), row.names = F)
write.csv(ground.train, file.path('finger_sim', 'ground.train.csv'), row.names = F)

