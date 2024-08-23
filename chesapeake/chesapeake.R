library(tidyverse)
library(magrittr)
library(sf)
library(viridis)
library(BRISC)
library(MLmetrics)
library(ggplot2)
library(data.table)
library(xtable)
library(raster)
library(R.matlab)
library(geosphere)
library(boraGP)
library(MBA)
library(foreach)
source('../functions.R')

nu = .5

# Large display for plots
my_theme <- theme_gray() +
  theme(
    plot.title = element_text(size = 24, hjust = 0.5),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 20),
    legend.text=element_text(size=20),
    legend.title=element_text(size=20),
  )

# Read shapefiles using st_read, which creates sf objects
water <- st_read("mygeodata/chesapeake-polygon.shp")
chesa.detail <- st_read("Chesapeake_Bay_Shoreline_Medium_Resolution/Chesapeake_Bay_Shoreline_Medium_Resolution.shp")
projcrs <- st_crs(water)

# Read and process CSV data
ph.dat <- fread('WaterQualityWaterQualityStation.csv', quote="") %>%
  filter(Parameter == "\"PH\"") %>%
  group_by(Station) %>%
  summarise(y = mean(MeasureValue, na.rm=T),
            lat = mean(Latitude), lon = mean(Longitude)) %>% # latitude and longitude are constant for each station
  dplyr::select(lat, lon, y)

# Convert ph.dat to an sf object
ph.sf <- st_as_sf(ph.dat, coords = c('lon', 'lat'), crs = st_crs(4326))  # 4326 is WGS84

# Transform coordinates to match CRS
ph.sf <- st_transform(ph.sf, projcrs)

# make demonstration grid
bbox <- st_bbox(water)
x_range <- bbox$xmax - bbox$xmin
y_range <- bbox$ymax - bbox$ymin
n_points <- 1000
grid_spacing <- sqrt((x_range * y_range) / n_points)
grid <- st_make_grid(water, cellsize = grid_spacing, what = 'centers')
grid <- st_intersection(grid, water)
grid_sf <- st_sf(geometry = grid)


# attach grid to coords
obs_coords <- st_geometry(ph.sf)
obs_coords_matrix <- st_coordinates(ph.sf)
grid_matrix <- st_coordinates(grid_sf)[,1:2]
all_coords <- rbind(obs_coords_matrix, grid_matrix)
dist_matrix <- distm(all_coords) # geodesic distances
mds_coords <- cmdscale(dist_matrix, k = 2)  # project to R^2 for convenience

# normalize coordinates
centroid <- colMeans(mds_coords)
translated_coords <- sweep(mds_coords, 2, centroid, FUN = "-")
scaling_factor <- sqrt(1 / mean(rowSums(translated_coords^2)))
coords.norm <- translated_coords * scaling_factor


# get point assignments
n <- nrow(ph.dat)

grid.id = (n+1):nrow(mds_coords)
grid.coords = coords.norm[grid.id,]

# function to return adjacency of pairs
con <- function(i, p.grid) {
  results <- vector("logical", length = i)

  for (j in 1:i) {
    # Generate great circle points
    gc_points <- gcIntermediate(p.grid[i, ], p.grid[j, ], n = 100, addStartEnd = TRUE)

    # Convert the points to a LineString and then to an sfc object
    gc_line <- st_sfc(st_linestring(gc_points), crs = st_crs(water))

    # Check if the great circle line is completely covered by the water body
    is_covered <- st_covers(st_geometry(water), gc_line, sparse=F)

    results[j] <- is_covered[1,1]
  }

  print(i)
  return(results)
}

# create adjacency matrix
if(file.exists('chesa_adjacency.rds')){
    A = readRDS('chesa_adjacency.rds')
} else{
    out.list = sapply(1:nrow(all_coords), FUN = con, p.grid = all_coords)
    A = matrix(nrow = nrow(all_coords), ncol=nrow(all_coords))
    for(i in 1:nrow(A)){
      A[i, 1:i] = out.list[[i]]
    }
    A[upper.tri(A)] = t(A)[upper.tri(A)]
    diag(A) = 0
    hist(rowSums(A), breaks = 15)
}

saveRDS(ph.dat, 'ph.rds')
saveRDS(A, 'chesa_adjacency.rds')

# raw pH plot
p1=ggplot() + 
  ggtitle("Chesapeake pH") +
  geom_sf(data = chesa.detail, size = .3, color = "black", fill = "cyan1") + 
  geom_sf(mapping = aes(color = y), data = ph.sf) +
  coord_sf() +
  scale_color_viridis(name='pH') + my_theme
ggsave('pH.png', p1)


# buffered boundary plot
p2=ggplot() + 
  ggtitle("Chesapeake pH - buffered boundary") +
  geom_sf(data = water, size = .3, color = "black", fill = "cyan1") + 
  geom_sf(mapping = aes(color = y), data = ph.sf) +
  coord_sf() +
  scale_color_viridis(name='pH') + my_theme
ggsave('PH_buff.png', p2)


# Accuracy experiment with multiple folds 
set.seed(89)
n_folds <- 3
folds <- sample(rep(1:n_folds, length.out = n))

mse_results <- matrix(0, nrow = 3, ncol = 3)
coverage_results <- matrix(0, nrow = 3, ncol = 3) 
length_results <- matrix(0, nrow = 3, ncol = 3)

for(i in 1:n_folds){
  test.id <- which(folds==i)
  train.id <- which(folds!=i)
  
  train.dat <- ph.dat[train.id,]
  test.dat <- ph.dat[test.id,]
  
  # BRISC model to inform starting parameters
  BRISC.mod = BRISC_estimation(coords = coords.norm[train.id,], y = train.dat$y, cov.model = 'exponential',
                               n.neighbors=15)
  BRISC.pred = BRISC_prediction(BRISC.mod, coords.0 = coords.norm[test.id,])
  
  
  # bora setup
  nu <- 0.5
  phi = BRISC.mod$Theta[3]
  tau.sq = BRISC.mod$Theta[2]
  sigma.sq = sqrt(BRISC.mod$Theta[1])
  starting <- list("phi" = phi, "sigma.sq" = sigma.sq, 
                   "tau.sq" = tau.sq, "nu" = nu)
  tuning <- list("phi" = 0.5, "sigma.sq" = 0.5, "tau.sq" = 0.1, "nu" = 0)         # nu is fixed 
  priors <- list("phi.Unif" = c(1e-6, 10*phi),                              
                 "sigma.sq.IG" = c(2, sigma.sq), "tau.sq.IG" = c(2, tau.sq+.01), 
                 "nu.Unif" = c(nu-0.1,nu+0.1))
  
  # Convert the extent of 'water' to a polygon and then to an sf object
  water_bbox <- st_bbox(water)
  water_bbox_polygon <- st_as_sfc(st_bbox(water_bbox), crs = st_crs(water))
  
  # Calculate the difference to find boundary
  bound <- st_difference(water_bbox_polygon, water)
  
  # fit bora model
  m = 15
  n.samples=10000
  burn=5000
  ord = order(all_coords[train.id,1])
  barrier_nninfo_all <- barrier_neighbor(coords = all_coords[train.id,], coords.0 =  all_coords[test.id,], ord = ord,
                                         n.neighbors = m,
                                         barrier = st_as_sf(bound),
                                         cores = 1,
                                         verbose = T,
                                         debug = list(barrier_n.indx = NULL,
                                                      barrier_dist = NULL,
                                                      barrier_nn.indx.0_list = NULL,
                                                      barrier_dist0 = NULL,
                                                      ref_window = NULL,
                                                      nonref_window = NULL,
                                                      ref_fill = TRUE,
                                                      nonref_fill = TRUE))
  
  barrier_nn.indx.0_list <- barrier_nninfo_all$barrier_nn.indx.0_list
  barrier_nn.indx.0 <- do.call("rbind", barrier_nn.indx.0_list)
  barrier_n.indx <- barrier_nninfo_all$barrier_n.indx
  barrier_nn.indx <- as.integer(unlist(barrier_n.indx)[-1]-1)
  barrier_nn.indx.lu <- c(cumsum(c(0,0,sapply(barrier_n.indx, length)[-1]))[1:nrow(train.dat)],
                          c(0,sapply(barrier_n.indx, length)[-1])) %>% as.integer()
  barrier_nninfo <- list(type = "barrier",
                         n.indx = barrier_n.indx,
                         n.neighbors = m, nn.indx = barrier_nn.indx,
                         nn.indx.lu = barrier_nn.indx.lu, ord = ord)
  barrier_m.s <- spNNGP(train.dat$y~ 1, coords = all_coords[train.id,], starting = starting,
                        method = "response", n.neighbors = m,
                        tuning = tuning, priors = priors, 
                        cov.model = "matern",
                        n.samples = n.samples, n.omp.threads = 1,
                        neighbor.info = barrier_nninfo, verbose = T)
  barrier_p.s <- predict(barrier_m.s, X.0 = matrix(1, nrow = nrow(test.dat), ncol = 1),
                         coords.0 = as.matrix(all_coords[test.id,]),
                         sub.sample = 
                           list(start = burn+1, end = n.samples, thin = 1),
                         nn.indx.0 = barrier_nn.indx.0, 
                         n.omp.threads = 1, verbose = T)
  bora.pred.vals <- rowMeans(barrier_p.s$p.y.0)
  yquantBRGP <- apply(barrier_p.s$p.y.0, 1, 
                      function(x) quantile(x, probs = c(0.025, 0.975)))
  
  # fit visGP
  visGP.fit = fit_water(coords.norm[train.id,], A[train.id, train.id], y= ph.dat$y[train.id])
  
  
  visGP.pred = sapply(test.id, function(i)
    pred_neighbor(new_loc = coords.norm[i,], A=A[train.id,train.id], A_vec = A[i, train.id], 
                  coords.norm[train.id,], ph.dat$y[train.id], visGP.fit$par, nu=.5, 
                  method = 'maxprec', n_neighbors = 15,
                  D=as.matrix(dist(coords.norm[train.id, ]))) )
  visGP.pred.vals = visGP.pred[1,] %>% unlist
  
  # Calculate metrics
  BRISC.res = c(MSE(BRISC.pred$prediction, ph.dat$y[test.id]),
                mean(sapply(1:length(test.id), function(i) BRISC.pred$prediction.ci[i,1] <= test.dat$y & BRISC.pred$prediction.ci[i,2] >= test.dat$y)),
                mean(sapply(1:length(test.id), function(i) BRISC.pred$prediction.ci[i,2] - BRISC.pred$prediction.ci[i,1])))
  
  visGP.res = c(MSE(visGP.pred.vals, ph.dat$y[test.id]), 
             mean(sapply(1:length(test.id), function(i) visGP.pred.vals[i]-qnorm(.975)*visGP.pred[[2,i]] <= test.dat$y &  visGP.pred.vals[i]+qnorm(.975)*visGP.pred[[2,i]]  >= test.dat$y)),
             mean(sapply(1:length(test.id), function(i) 2*qnorm(.975)*visGP.pred[[2,i]] )))
  
  BORA.res = c(MSE(bora.pred.vals, ph.dat$y[test.id]),
               mean(sapply(1:length(test.id), function(i) yquantBRGP[1, i] <= test.dat$y & yquantBRGP[2, i]  >= test.dat$y)),
               mean(sapply(1:length(test.id), function(i) yquantBRGP[2, i] - yquantBRGP[1, i] )))
  
  # Store results           
  mse_results[i,] <- c(BRISC.res[1], visGP.res[1], BORA.res[1])
  coverage_results[i,] <- c(BRISC.res[2], visGP.res[2], BORA.res[2])
  length_results[i,] <- c(BRISC.res[3], visGP.res[3], BORA.res[3])
  
  # save data
  write.csv(ph.dat$y[test.id], paste0('test_y_', i, '.csv'), row.names=F)
  write.csv(coords.norm[test.id,], paste0('test_loc_', i, '.csv'), row.names=F)
  write.csv(ph.dat$y[train.id], paste0('train_y_', i, '.csv'), row.names=F)
  write.csv(coords.norm[train.id,], paste0('train_loc_', i, '.csv'), row.names=F)
}

# save results
saveRDS(list(mse_results=mse_results, coverage_results=coverage_results, length_results=length_results), 
        'res.rds')

# full data models for demonstration 
BRISC.mod = BRISC_estimation(coords = coords.norm[1:n,], y = ph.dat$y, cov.model = 'exponential',
                             n.neighbors=15)
visGP.fit = fit_water(coords.norm[1:n,], A[1:n, 1:n], y = ph.dat$y[1:n])

BRISC.grid.pred <- BRISC_prediction(BRISC.mod, grid.coords)
BRISC.pred.vals.grid <- BRISC.grid.pred$prediction

visGP.pred.grid = sapply(grid.id, function(i)
  pred_neighbor(new_loc = coords.norm[i,], A=A[1:n,1:n], A_vec = A[i, 1:n], 
                coords.norm[1:n,], ph.dat$y[1:n], visGP.fit$par, nu=.5, 
                method = 'maxprec', n_neighbors = 15,
                D=as.matrix(dist(coords.norm[1:n, ]))) )
visGP.pred.vals.grid = visGP.pred.grid[1,] %>% unlist

grid.sf <- st_as_sf(data.frame(all_coords[grid.id,]), coords=c('X','Y'), crs=st_crs(4326))
diff.sf <- st_as_sf(grid.sf %>% mutate(diff = visGP.pred.vals.grid - BRISC.pred.vals.grid)
                         , coords = c('lon', 'lat'), crs = 'WGS84')


# Create a dense grid of coordinates
dense_grid <- expand.grid(x = seq(min(all_coords[,1]), max(all_coords[,1]), length.out = 100),
                          y = seq(min(all_coords[,2]), max(all_coords[,2]), length.out = 100))

# Interpolate BRISC predictions on the dense grid
BRISC.dense.pred <- mba.surf(cbind(all_coords[grid.id,1], all_coords[grid.id,2], BRISC.pred.vals.grid),
                             no.X = 100, no.Y = 100, extend = TRUE)$xyz.est$z

# Interpolate visGP predictions on the dense grid  
visGP.dense.pred <- mba.surf(cbind(all_coords[grid.id,1], all_coords[grid.id,2], visGP.pred.vals.grid),
                          no.X = 100, no.Y = 100, extend = TRUE)$xyz.est$z

# Create sf object with dense grid and differences
dense.grid.sf <- st_as_sf(data.frame(dense_grid), coords=c('x','y'), crs=st_crs(4326))
diff.dense.sf <- st_as_sf(dense.grid.sf %>% 
                            mutate(diff = as.vector(visGP.dense.pred - BRISC.dense.pred)),
                          coords = c('x', 'y'), crs = 'WGS84')

# Transform the CRS of chesa.detail to match diff.dense.sf
chesa.detail <- st_transform(chesa.detail, crs = st_crs(diff.dense.sf))

# Intersect the differences with the Chesapeake Bay shoreline
diff.chesa.sf <- st_intersection(diff.dense.sf, st_make_valid(chesa.detail))

# dense grid plot
p8=ggplot() + 
  ggtitle("visGP minus BRISC (Euclidean)") +
  geom_sf(data = diff.chesa.sf, size = .3, color = "black", fill = "cyan1") + 
  geom_sf(mapping = aes(color = diff), data = diff.chesa.sf) +
  coord_sf() +
  scale_color_viridis() + my_theme
ggsave('visGP_BRISC_diff.png', p8)

# Interpolate visGP standard deviations on the dense grid  
visGP.dense.ci <- mba.surf(cbind(all_coords[grid.id,1], all_coords[grid.id,2], 2*qnorm(.975)*visGP.pred.grid[2,] %>% unlist),
                        no.X = 100, no.Y = 100, extend = TRUE)$xyz.est$z

# Create sf object with dense grid and standard deviations
visGP.ci.dense.sf <- st_as_sf(dense.grid.sf %>% 
                          mutate(ci = as.vector(visGP.dense.ci)),
                        coords = c('x', 'y'), crs = 'WGS84')

# Transform the CRS of chesa.detail to match sd.dense.sf
chesa.detail <- st_transform(chesa.detail, crs = st_crs(visGP.ci.dense.sf))

# Intersect the standard deviations with the Chesapeake Bay shoreline
visGP.ci.chesa.sf <- st_intersection(visGP.ci.dense.sf, st_make_valid(chesa.detail))

# dense grid plot
p9 <- ggplot() + 
  ggtitle("visGP predictive interval width") +
  geom_sf(data = visGP.ci.chesa.sf, size = .3, color = "black", fill = "cyan1") + 
  geom_sf(mapping = aes(color = ci), data = visGP.ci.chesa.sf) +
  coord_sf() +
  scale_color_viridis() + 
  labs(color = "width") +  # Change legend title for color
  my_theme



# Interpolate BRISC standard deviations on the dense grid  
BRISC.dense.ci <- mba.surf(cbind(all_coords[grid.id,1], all_coords[grid.id,2], rowDiffs(BRISC.grid.pred$prediction.ci)),
                        no.X = 100, no.Y = 100, extend = TRUE)$xyz.est$z

# Create sf object with dense grid and standard deviations
BRISC.ci.dense.sf <- st_as_sf(dense.grid.sf %>% 
                          mutate(ci = as.vector(BRISC.dense.ci)),
                        coords = c('x', 'y'), crs = 'WGS84')

# Transform the CRS of chesa.detail to match sd.dense.sf
chesa.detail <- st_transform(chesa.detail, crs = st_crs(BRISC.ci.dense.sf))

# Intersect the standard deviations with the Chesapeake Bay shoreline
BRISC.ci.chesa.sf <- st_intersection(BRISC.ci.dense.sf, st_make_valid(chesa.detail))

# dense grid plot
p10 <- ggplot() + 
  ggtitle("BRISC (Euclidean) \npredictive interval width") +
  geom_sf(data = BRISC.ci.chesa.sf, size = .3, color = "black", fill = "cyan1") + 
  geom_sf(mapping = aes(color = ci), data = BRISC.ci.chesa.sf) +
  coord_sf() +
  scale_color_viridis() + my_theme
ggsave('BRISC_pi.png', p10)








