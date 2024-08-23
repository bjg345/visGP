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
library(remotes)
library(MBA)
library(foreach)
source('../functions.R')

library(igraph)
library(gRbase)
library(geostatsp)

nu = .5

# Read shapefiles using st_read, which creates sf objects
water <- st_read("mygeodata/chesapeake-polygon.shp")
chesa.detail <- st_read("Chesapeake_Bay_Shoreline_Medium_Resolution/Chesapeake_Bay_Shoreline_Medium_Resolution.shp")
projcrs <- st_crs(water)

# Read and process CSV data
ph.dat <- fread('WaterQualityWaterQualityStation.csv', quote="") %>%
  filter(Parameter == "\"PH\"") %>%
  group_by(Station) %>%
  summarise(y = mean(MeasureValue, na.rm=T),
            lat = mean(Latitude), lon = mean(Longitude)) %>%
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
#saveRDS(A, 'chesa_adjacency.rds')

# raw pH plot
p1=ggplot() + 
  ggtitle("Chesapeake pH") +
  geom_sf(data = chesa.detail, size = .3, color = "black", fill = "cyan1") + 
  geom_sf(mapping = aes(color = y), data = ph.sf) +
  coord_sf() +
  scale_color_viridis(name='pH',option="magma")
#ggsave('pH.png', p1)


# buffered boundary plot
p2=ggplot() + 
  ggtitle("Chesapeake pH - buffered boundary") +
  geom_sf(data = water, size = .3, color = "black", fill = "cyan1") + 
  geom_sf(mapping = aes(color = y), data = ph.sf) +
  coord_sf() +
  scale_color_viridis(name='pH',option="magma")

# Accuracy experiment with multiple folds 
set.seed(89)
nfolds=n
folds=1:n

mse_results <- matrix(0, nrow = nfolds, ncol = 3)
coverage_results <- matrix(0, nrow = nfolds, ncol = 3) 
length_results <- matrix(0, nrow = nfolds, ncol = 3)

for(i in 1:nfolds){
  print(paste0("Fold: ",i))
  test.id <- which(folds==i)
  train.id <- which(folds!=i)
  
  train.dat <- ph.dat[train.id,]
  test.dat <- ph.dat[test.id,,drop=F]
  
  # BRISC model to inform starting parameters
  BRISC.mod = BRISC_estimation(coords = coords.norm[train.id,], y = train.dat$y, cov.model = 'exponential',
                               n.neighbors=15)
  BRISC.pred = BRISC_prediction(BRISC.mod, coords.0 = coords.norm[test.id,,drop=F])
  
  # fit visGP
  visGP.fit = fit_water(coords.norm[train.id,], A[train.id, train.id], y= ph.dat$y[train.id])
  
  
  visGP.pred = sapply(test.id, function(i)
    pred_neighbor(new_loc = coords.norm[i,], A=A[train.id,train.id], A_vec = A[i, train.id], 
                  coords.norm[train.id,], ph.dat$y[train.id], visGP.fit$par, nu=.5, 
                  method = 'maxprec', n_neighbors = 15,
                  D=as.matrix(dist(coords.norm[train.id, ]))) )
  visGP.pred.vals = visGP.pred[1,] %>% unlist
  
  # Calculate metrics
  BRISC.res = c(MSE(BRISC.pred$prediction, ph.dat$y[test.id,drop=F]),
                mean(sapply(1:length(test.id), function(i) BRISC.pred$prediction.ci[i,1] <= test.dat$y & BRISC.pred$prediction.ci[i,2] >= test.dat$y)),
                mean(sapply(1:length(test.id), function(i) BRISC.pred$prediction.ci[i,2] - BRISC.pred$prediction.ci[i,1])))
  
  visGP.res = c(MSE(visGP.pred.vals, ph.dat$y[test.id]), 
             mean(sapply(1:length(test.id), function(i) visGP.pred.vals[i]-qnorm(.975)*visGP.pred[[2,i]] <= test.dat$y &  visGP.pred.vals[i]+qnorm(.975)*visGP.pred[[2,i]]  >= test.dat$y)),
             mean(sapply(1:length(test.id), function(i) 2*qnorm(.975)*visGP.pred[[2,i]] )))
  
  # Store results           
  mse_results[i,] <- c(BRISC.res[1], visGP.res[1], NA)
  coverage_results[i,] <- c(BRISC.res[2], visGP.res[2], NA)
  length_results[i,] <- c(BRISC.res[3], visGP.res[3], NA)
}

round(100*colMeans(mse_results),2)
round(100*colMeans(coverage_results),1)
round(colMeans(length_results),2)

# save results
saveRDS(list(mse_results=mse_results, coverage_results=coverage_results, length_results=length_results), 
        'loores.rds')
