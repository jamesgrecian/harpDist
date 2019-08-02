##############################################################################
### Harp seal INLA model with seasonal spatial field estimated across time ###
### May 2019                                                               ###
##############################################################################

# Following chat with Joe Watson it should be possible to estimate seasonal spatial field across all years combined
# Each year should contribute to the estimation of the seasonal pattern
# Still need to include each season:year combination BUT add extra season index for time mesh
# So sea ice data is extracted to 20 layers BUT those 20 layers are instead indexed by season in the time mesh
# The inla spatial random effect can operate over 4 seasons even with the 20 season x year combination
# Just index season 1 to 4 5 times... so s.group should be 1 to 4 not 1 to 20.

# Load libraries
require(inlabru)
require(INLA)
require(tidyverse)
require(sf)
require(viridis)
source("discrete_gradient.R")

# Load data
dat <- readRDS("harp data/harps2500_indexed.rds")
dat <- dat %>% mutate(x = x/1000,
                      y = y/1000)
# Define albers projection
prj = "+proj=laea +lat_0=75 +lon_0=-25 +x_0=0 +y_0=0 +datum=WGS84 +units=km +no_defs +ellps=WGS84 +towgs84=0,0,0"

# Generate land map
land <- mapr::mapr(dat,
                   prj,
                   buff = 2000)

# Generate simple boundary for inla mesh
b <- mapr::meshr(dat,
                 prj,
                 buff = 500,
                 keep = 0.5,
                 Neumann = F) 

# Mesh boundary parameters
max.edge = max(c(diff(range(dat$x)), diff(range(dat$y))))/15
bound.outer = 1000

# Create mesh
mesh = inla.mesh.2d(boundary = b,
                    loc = cbind(dat$x, dat$y),
                    max.edge = c(1, 5) * max.edge,
                    cutoff = 25,
                    offset = c(max.edge, bound.outer))

# plot to check
ggplot() +
  theme_bw() + ylab("") + xlab("") +
  inlabru::gg(mesh) + 
  geom_sf(aes(), data = st_as_sf(dat, coords = c("lon", "lat"))
          %>% st_set_crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) +
  coord_sf(xlim = c(-4000, 3000), ylim = c(-4000, 3000), crs = prj, expand = F) +
  ggtitle("Haakon barrier model")

# Remove points on land
in.water = over(b, SpatialPoints(cbind(dat$x, dat$y), proj4string = CRS(prj)), returnList=T)[[1]]
print(paste("There are", nrow(dat)-length(in.water), "points on land in the original polygon"))
dat <- dat[in.water,]

# Set up boundary matern model
tl = length(mesh$graph$tv[,1]) # the number of triangles in the mesh

# compute triangle positions
posTri = matrix(0, tl, 2) 
for (t in 1:tl){
  temp = mesh$loc[mesh$graph$tv[t, ], ]
  posTri[t,] = colMeans(temp)[c(1,2)] 
}
posTri = SpatialPoints(posTri, proj4string = CRS(prj))

normal = over(b, posTri, returnList = T) # check which mesh triangles are inside the normal area
normal = unlist(normal)

barrier.triangles = setdiff(1:tl, normal)
poly.barrier = inla.barrier.polygon(mesh, barrier.triangles)

# create the matern object
barrier.model = inla.barrier.pcmatern(mesh,
                                      barrier.triangles = barrier.triangles,
                                      prior.range = c(50, .1),
                                      prior.sigma = c(10, 0.01))


# Harp seal locations need to be put into an inla stack...
# For the lgcp make sure the locations are on mesh nodes
# Then set mesh nodes to 1 where locations are, and zero where not
nv <- mesh$n # number of mesh nodes
n <- nrow(dat) # number of data points

# assigning weight 0 to all the nodes of the mesh that fall on land
ips <- ipoints(mesh)
ips <- st_as_sf(ips)
st_crs(ips) <- prj
ips$weight[sapply(st_intersects(ips, land),function(x){length(x)>0}),] <- 0
table(ips$weight > 0) # check

# Create a 1d time mesh for the time series (something weird with using cyclic)
# This can be seasonal (1-4) but replicated across the 5 year bins...
tmesh <- inla.mesh.1d(loc = 1:4, boundary = "free")
tmesh$loc
(k <- length(tmesh$loc))

# relabel spatial mesh
smesh <- mesh

# Find Voroni polygons for mesh nodes
library(deldir)
dd <- deldir(smesh$loc[, 1], smesh$loc[, 2])
tiles <- tile.list(dd)

# Convert to Spatial Polygons
polys <- SpatialPolygons(lapply(1:length(tiles), function(i) {
  p <- cbind(tiles[[i]]$x, tiles[[i]]$y)
  n <- nrow(p)
  Polygons(list(Polygon(p[c(1:n, 1), ])), i)
}))

# Assign points to each polygon
area <- factor(over(SpatialPoints(cbind(dat$x, dat$y)), polys),
               levels = 1:length(polys))

# Replicate dataframe to match number of time indices
time <- dat$index
table(time)

# Use both space and time index to aggregate data
agg.dat <- as.data.frame(table(area, time))
for(j in 1:2) # set time and area as integer
  agg.dat[[j]] <- as.integer(as.character(agg.dat[[j]])) 
str(agg.dat)

# INLA input dataframe now has mesh nodes replicated 20 times
# All season:year combinations are represented

# Load and stack all 20 raster layers together
ice <- raster::stack(raster::stack("NSIDC Sea Ice/Seasonal_NSIDC_Dec94_Nov99"),
                     raster::stack("NSIDC Sea Ice/Seasonal_NSIDC_Dec99_Nov04"),
                     raster::stack("NSIDC Sea Ice/Seasonal_NSIDC_Dec04_Nov09"),
                     raster::stack("NSIDC Sea Ice/Seasonal_NSIDC_Dec09_Nov14"),
                     raster::stack("NSIDC Sea Ice/Seasonal_NSIDC_Dec14_Nov19"))

# Extract data needed by unique times..
agg.dat$ice <- NA

for (i in sort(unique(time))){
  agg.dat$ice[agg.dat$time == i] <- raster::extract(raster::subset(ice, i),
                                                    ips %>% st_transform("+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +a=6378273 +b=6356889.449 +units=m +no_defs"),
                                                    method = "bilinear")
}

# average across years
av_ice <- raster::stackApply(ice, indices = rep(1:4, times = 5), fun = mean, na.rm = T)

# append seasonal average to agg.dat
unique(agg.dat$time)

for (i in 1:4){
  av_ice_vals <- raster::extract(raster::subset(av_ice, i),
                                 ips %>% st_transform("+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +a=6378273 +b=6356889.449 +units=m +no_defs"),
                                 method = "bilinear")
  agg.dat$av_ice[agg.dat$time %in% seq(i, 20, by = 4)] <- rep(av_ice_vals,
                                                              times = length(unique(agg.dat$time[agg.dat$time %in% seq(i, 20, by = 4)])))
}


# Expected number of points needs to be defined as
# proportional to the area of the polygon * the width of time knot
w.areas <- sapply(1:length(tiles), function(i) {
  p <- cbind(tiles[[i]]$x, tiles[[i]]$y)
  n <- nrow(p)
  pl <- SpatialPolygons(
    list(Polygons(list(Polygon(p[c(1:n, 1),])), i)), proj4string = CRS(prj))
  if (gIntersects(pl, b))
    return(gArea(gIntersection(pl, b)))
  else return(0)
})

# Summary of the polygon areas
summary(w.areas)
sum(w.areas)

# Check with total area
s.area <- gArea(b)
s.area

# Width of each knot
w.t <- diag(inla.mesh.fem(tmesh)$c0)
w.t

# Estimate intensity function
i0 <- n / (gArea(b) * diff(range(tmesh$loc)))
c(i0, log(i0))

# The space-time volume at each polygon and time knot
e0 <- w.areas[agg.dat$area] * (w.t[agg.dat$time])
summary(e0)

# Add season index to the data frame as this is what the time mesh will operate on
agg.dat$season = NA
agg.dat$season[agg.dat$time %in% c(1, 5, 9, 13, 17)] <- 1
agg.dat$season[agg.dat$time %in% c(2, 6, 10, 14, 18)] <- 2
agg.dat$season[agg.dat$time %in% c(3, 7, 11, 15, 19)] <- 3
agg.dat$season[agg.dat$time %in% c(4, 8, 12, 16, 20)] <- 4

# Create projector matrix - base this on mesh not barrier model?
A.st <- inla.spde.make.A(mesh = smesh,
                         loc = smesh$loc[agg.dat$area, ],
                         group = agg.dat$season,
                         mesh.group = tmesh)

# Create space-time index
idx <- inla.spde.make.index(name = 's',
                            n.spde = barrier.model$f$n,
                            n.group = 4)

# Define the data stack for the inla model
stk <- inla.stack(
  data = list(y = agg.dat$Freq, exposure = e0), 
  A = list(A.st, 1), 
  effects = list(idx, list(b0 = rep(1, nrow(agg.dat)),
                           ice_av = agg.dat$av_ice,
                           ice_dev = agg.dat$ice - agg.dat$av_ice)))

# PC prior on temporal correlation
pcrho <- list(theta = list(prior = 'pccor1', param = c(.7, .7))) # order is mu and alpha (1/sd^2)

# Model formula
# Model failed when passing all depth data - so group it instead...
# see ?inla.group for more info
f_1 <- y ~ 0 + b0 +
  f(s, model = barrier.model, group = s.group, control.group = list(model = 'ar1', hyper = pcrho)) +
  f(inla.group(ice_av, n = 25, method = "cut"), model = 'rw2', scale.model = T, hyper = list(theta = list(prior = "pc.prec", param = c(120, 0.01)))) +
  f(inla.group(ice_dev, n = 25, method = "cut"), model = 'rw2', scale.model = T, hyper = list(theta = list(prior = "pc.prec", param = c(120, 0.01))))

start_vals <- readRDS("start_vals.rds")

# Fit the model
m_1 <- inla(f_1,
            family = 'poisson', 
            data = inla.stack.data(stk),
            control.predictor = list(A = inla.stack.A(stk)),
            E = exposure,
            control.inla = list(int.strategy = "eb"), # strategy ='adaptive' fails
            control.compute = list(config = TRUE,
                                   dic = T,
                                   waic = T),
            control.mode = list(theta = start_vals, restart = TRUE), # provide starting values from previous model runs
            verbose = T) # switch on when trialling



