library(tidyverse)
library(magrittr)
library(sf)
library(blockCV)
library(viridis)
library(BRISC)
library(MLmetrics)
library(rgeos)
library(ggplot2)
library(data.table)
library(xtable)
library(raster)
library(gdistance)
library(R.matlab)
library(boraGP)
library(foreach)
source('../functions.R')

nu = .5

chesa = st_read("mygeodata/chesapeake_3-polygon.shp")
chesa.detail = st_read("Chesapeake_Bay_Shoreline_Medium_Resolution/Chesapeake_Bay_Shoreline_Medium_Resolution.shp")
projcrs <- st_crs(chesa)

#x = extent(chesa)
#buff.rec = 10000
#x@xmin=x@xmin-buff.rec; x@xmax = x@xmax+buff.rec; x@ymin = x@ymin-buff.rec; x@ymax = x@ymax+buff.rec
#ras <- raster(nrow = 1000, ncol = 1000, ext = x)
#values(ras) <- 1
#ras <- mask(ras, chesa)
#tr <- transition(ras, transitionFunction=mean, directions=8) %>% geoCorrection('c')

ph.dat= fread('WaterQualityWaterQualityStation.csv', quote="") %>%
  filter(Parameter =="\"PH\""  ) %>%
  group_by(Station) %>%
  summarise(y = mean(MeasureValue, na.rm=T),
            lat = mean(Latitude), lon = mean(Longitude)) %>% #checked no variation in lat/lon
  dplyr::select(lat, lon, y) 

ph.sf <- st_as_sf(ph.dat, coords = c('lon', 'lat'), crs = 'WGS84')

water <- as_Spatial(chesa)# %>% spTransform(CRS('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'))
coords <- SpatialPoints(ph.dat[, c('lon', 'lat')],
                        proj4string = CRS('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')) %>% 
  spTransform(projcrs$wkt)




#dists <- costDistance(tr, coords@coords, coords@coords)
#diag(dists) = NA


con = function(i, p.grid) {
  require(sp)
  
  n = nrow(p.grid)
  
  l <- vector("list", i)
  
  for (j in 1:i) {
    l[[j]] <- Lines(list(Line(as.matrix(rbind(p.grid[i,], p.grid[j,])))), as.character(j))
  }
  
  out = gCovers(water, 
                SpatialLines(l, proj4string=water@proj4string), byid=T)
  print(i)
  return(out)
  
}

out.list = sapply(1:nrow(ph.dat), FUN = con, p.grid = coords@coords)
A = matrix(nrow = nrow(ph.dat), ncol=nrow(ph.dat))
for(i in 1:nrow(A)){
  A[i, 1:i] = out.list[[i]]
}
A[upper.tri(A)] = t(A)[upper.tri(A)]
diag(A) = 0
hist(rowSums(A), breaks = 15)


saveRDS(ph.dat, 'ph.rds')

n <- nrow(ph.dat)

set.seed(89)
blocked =spatialBlock(coords, rows=10, cols=10, selection='checkerboard') #blocked test/train desgin
blocked$foldID[69]=2
blocked$foldID[136]=1
assignment = vector('character', length=length(blocked$foldID))
for(i in 1:length(assignment)){
	if(blocked$foldID[i] == 1) assignment[i] = 'train'
	else assignment[i] = 'test'
}

p1=ggplot() + 
  ggtitle("Chesapeake pH") +
  geom_sf(data = chesa.detail, size = .3, color = "black", fill = "cyan1") + 
  geom_sf(mapping = aes(color = y), data = ph.sf) +
  coord_sf() +
  scale_color_viridis(name='pH')
ggsave('pH.png', p1)

p0=ggplot() +
  ggtitle("Chesapeake train/test assignment") +
  geom_sf(data = chesa.detail, size = .3, color = "black", fill = "cyan1") +
  geom_sf(mapping = aes(color= assignment), data = ph.sf) +
  scale_color_manual(values = c("test" = "red",
                                "train" = "green")) + 
  coord_sf() 
ggsave('assignment.png', p0)



p2=ggplot() + 
  ggtitle("Chesapeake pH - buffered boundary") +
  geom_sf(data = chesa, size = .3, color = "black", fill = "cyan1") + 
  geom_sf(mapping = aes(color = y), data = ph.sf) +
  coord_sf() +
  scale_color_viridis(name='pH')
ggsave('PH_buff.png', p2)


train.id = which(blocked$foldID==1)
test.id = which(blocked$foldID==2)

train.dat <- ph.dat[train.id,]
test.dat <- ph.dat[test.id,]

#for(i in 1:nrow(A)){
  
 # neigh <- base::sort(dists[i,train.id], na.last = T, index.return=T)$ix[1:3]
 # A[i, train.id[neigh]] <- A[train.id[neigh], i] <- 1
  
  
#}



coords.norm <- scale(coords@coords)
BRISC.mod = BRISC_estimation(coords = coords.norm[train.id,], y = train.dat$y, cov.model = 'exponential',
                             n.neighbors=15)
BRISC.pred = BRISC_prediction(BRISC.mod, coords.0 = coords.norm[test.id,])



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

bound = gDifference(as(extent(water), "SpatialPolygons"), water)


m = 15
n.samples=10000
burn=5000
ord = order(coords.norm[train.id,1])
barrier_nninfo_all <- barrier_neighbor(coords = coords.norm[train.id,], coords.0 =  coords.norm[test.id,], ord = ord,
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
      
      barrier_m.s <- spNNGP(train.dat$y~ 1, coords = coords.norm[train.id,], starting = starting,
                             method = "response", n.neighbors = m,
                             tuning = tuning, priors = priors, 
                             cov.model = "matern",
                             n.samples = n.samples, n.omp.threads = 1,
                             neighbor.info = barrier_nninfo, verbose = T)
       barrier_p.s <- predict(barrier_m.s, X.0 = matrix(1, nrow = nrow(test.dat), ncol = 1),
                              coords.0 = as.matrix(coords.norm[test.id,]),
                              sub.sample = 
                                list(start = burn+1, end = n.samples, thin = 1),
                              nn.indx.0 = barrier_nn.indx.0, 
                             n.omp.threads = 1, verbose = T)
                             
       ystarBRGP <- rowMeans(barrier_p.s$p.y.0)
              yquantBRGP <- apply(barrier_p.s$p.y.0, 1, 
                           function(x) quantile(x, probs = c(0.025, 0.975)))
                           
                           


cs.fit = fit.water(coords.norm[train.id,], A[train.id, train.id], y= ph.dat$y[train.id] ,method = 'cov.select')


cs.pred = sapply(test.id, function(i)
  pred_neighbor(new.loc = coords.norm[i,], A=A[train.id,train.id], A.vec = A[i, train.id], 
                coords.norm[train.id,], ph.dat$y[train.id], cs.fit$par, nu=.5, 
                method = 'maxprec', n.neighbors = 15,
                D=as.matrix(dist(coords.norm[train.id, ]))) )
cs.pred.vals = cs.pred[1,] %>% unlist

bora.pred.vals = ystarBRGP
yquantBRGP <- apply(barrier_p.s$p.y.0, 1, 
                           function(x) quantile(x, probs = c(0.025, 0.975)))
bora.err <- bora.pred.vals - test.dat$y
bora.err.sf <- st_as_sf(ph.dat[test.id,] %>% mutate(err = bora.err)
                         , coords = c('lon', 'lat'), crs = 'WGS84')

BRISC.err <- BRISC.pred$prediction - test.dat$y
BRISC.err.sf <- st_as_sf(ph.dat[test.id,] %>% mutate(err = BRISC.err)
                         , coords = c('lon', 'lat'), crs = 'WGS84')


glgp.pred.vals = readMat('glgp_pred.mat')$vals.pred %>%as.vector
glgp.err <- glgp.pred.vals - test.dat$y
glgp.err.sf <- st_as_sf(ph.dat[test.id,] %>% mutate(err = glgp.err)
                         , coords = c('lon', 'lat'), crs = 'WGS84')

p3=ggplot() + 
  ggtitle("BRISC error") +
  geom_sf(data = chesa.detail, size = .3, color = "black", fill = "cyan1") + 
  geom_sf(mapping = aes(color = err), data = BRISC.err.sf) +
  coord_sf() +
  scale_color_viridis(limits=c(-.3, 1.11))
ggsave('BRISC_error.png', p3)

cs.err <- cs.pred.vals - test.dat$y
cs.err.sf <- st_as_sf(ph.dat[test.id,] %>% mutate(err = cs.err)
                         , coords = c('lon', 'lat'), crs = 'WGS84')
p4=ggplot() + 
  ggtitle("Covariance selection error") +
  geom_sf(data = chesa.detail, size = .3, color = "black", fill = "cyan1") + 
  geom_sf(mapping = aes(color = err), data = cs.err.sf) +
  coord_sf() +
  scale_color_viridis(limits=c(-.3,1.11))
ggsave('cs_error.png', p4)

p5=ggplot() + 
  ggtitle("BORA-GP error") +
  geom_sf(data = chesa.detail, size = .3, color = "black", fill = "cyan1") + 
  geom_sf(mapping = aes(color = err), data = bora.err.sf) +
  coord_sf() +
  scale_color_viridis(limits=c(-.3,1.11))
ggsave('bora_error.png', p5)

p5=ggplot() + 
  ggtitle("GLGP error") +
  geom_sf(data = chesa.detail, size = .3, color = "black", fill = "cyan1") + 
  geom_sf(mapping = aes(color = err), data = glgp.err.sf) +
  coord_sf() +
  scale_color_viridis(limits=c(-.3,1.11))
ggsave('glgp_error.png', p5)

BRISC.res = c(MSE(BRISC.pred$prediction, ph.dat$y[test.id]),
  mean(sapply(1:length(test.id), function(i) BRISC.pred$prediction.ci[i,1] <= test.dat$y & BRISC.pred$prediction.ci[i,2] >= test.dat$y)),
  mean(sapply(1:length(test.id), function(i) BRISC.pred$prediction.ci[i,2] - BRISC.pred$prediction.ci[i,1])))
  
GLGP.res = c(MSE(glgp.pred.vals, ph.dat$y[test.id]), NA, NA)

CS.res = c(MSE(cs.pred.vals, ph.dat$y[test.id]), 
  mean(sapply(1:length(test.id), function(i) cs.pred.vals[i]-qnorm(.975)*cs.pred[[2,i]] <= test.dat$y &  cs.pred.vals[i]+qnorm(.975)*cs.pred[[2,i]]  >= test.dat$y)),
  mean(sapply(1:length(test.id), function(i) 2*qnorm(.975)*cs.pred[[2,i]] )))
  
BORA.res = c(MSE(bora.pred.vals, ph.dat$y[test.id]),
  mean(sapply(1:length(test.id), function(i) yquantBRGP[1, i] <= test.dat$y & yquantBRGP[2, i]  >= test.dat$y)),
  mean(sapply(1:length(test.id), function(i) yquantBRGP[2, i] - yquantBRGP[1, i] )))

df = data.frame(Method = rep(NA, 4), MSE=rep(NA,4), Coverage=rep(NA,4), CI.length=rep(NA,4))
df$Method = c('BRISC', 'Covariance selection (maximum precision)', 'BORA-GP', 'GLGP')
df[1,2:4] = BRISC.res
df[2,2:4] = CS.res
df[3,2:4] = BORA.res
df[4,2:4] = GLGP.res

print(xtable(df, digits=-3), math.style.exponents=T, include.rownames=F)



#write.csv(ph.dat$y[test.id], 'test_y.csv', row.names=F)
#write.csv(coords.norm[test.id,], 'test_loc.csv', row.names=F)
#write.csv(ph.dat$y[train.id], 'train_y.csv', row.names=F)
#write.csv(coords.norm[train.id,], 'train_loc.csv', row.names=F)

