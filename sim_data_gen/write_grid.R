# Write simulated data to csv files for use by GLGP matlab scripts

load('u_data.Rdata')

write.csv(grid.train, file.path('../u_sim', 'grid.train.csv'), row.names = F)
write.csv(grid.test, file.path('../u_sim', 'grid.test.csv'), row.names = F)
write.csv(ground.test, file.path('../u_sim', 'ground.test.csv'), row.names = F)
write.csv(ground.train, file.path('../u_sim', 'ground.train.csv'), row.names = F)

load('fork_data.Rdata')
write.csv(grid.train, file.path('../fork_sim', 'grid.train.csv'), row.names = F)
write.csv(grid.test, file.path('../fork_sim', 'grid.test.csv'), row.names = F)
write.csv(ground.test, file.path('../fork_sim', 'ground.test.csv'), row.names = F)
write.csv(ground.train, file.path('../fork_sim', 'ground.train.csv'), row.names = F)

