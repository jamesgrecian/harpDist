###############################################
### Harp seal INLA barrier space-time model ###
###############################################

# The data are harp seal telemetry locations collated from studies going back to 1995
# We are interested in:
# 1. Estimating the seasonal north-south migration
# 2. Estimating the link between harp seal distribution and sea ice
# 3. Estimating how the migration may changed due to climate-change induced reductions in ice cover

# I have binned the data into seasons (Dec-Feb; Mar-May; Jun-Aug; Sep-Nov)
# and into year blocks (95-99; 00-04; 05-09; 10-14; 15-19)
# These are combined into an index 1:20 (4 seasons x 5 years)
# I will then be fitting covariates such as sea ice conditions that will vary across these 20 time indices

# In this model:
# 1. share the seasonal space field over all years
# as we don't have all populations tracked in all years

# 2. Seperate the populations
# as they are different sizes and don't overlap

# Load libraries
require(inlabru)
require(INLA)
require(tidyverse)
require(sf)

# You'll need my 'mapr' package from github
# custom functions for mapping and plotting
# devtools::install_github("jamesgrecian/mapr")
require(mapr)
# You might also need the rmapshaper package
# devtools::install_github("ateucher/rmapshaper")  

# Here is a random subsample of 2500 animal locations
dat <- readRDS("harp data/harps2500_indexed.rds")

# To make things a little easier later on convert the projected locations from metres to kilometres
dat <- dat %>% mutate(x = x/1000,
                      y = y/1000)
# Define an Albers projection with km rather than m
prj = "+proj=laea +lat_0=75 +lon_0=-25 +x_0=0 +y_0=0 +datum=WGS84 +units=km +no_defs +ellps=WGS84 +towgs84=0,0,0"

# For visualising the data later on generate a land shapefile using the mapr package
land <- mapr::mapr(dat,
                   prj,
                   buff = 2000)

# Generate simple boundary for the inla mesh
b <- mapr::meshr(dat,
                 prj,
                 buff = 500,
                 keep = 0.5,
                 Neumann = F) 

# Define the parameters of the boundary
max.edge = max(c(diff(range(dat$x)), diff(range(dat$y))))/15
bound.outer = 1000

# Changing cutoff is an effective way to change number of mesh nodes
# Fewer mesh nodes = less computation time...
# Create an inla mesh
mesh = inla.mesh.2d(boundary = b,
                    loc = cbind(dat$x, dat$y),
                    max.edge = c(1, 5) * max.edge,
                    cutoff = 50, # this has normall been 25 to match ice...
                    offset = c(max.edge, bound.outer),
                    crs = CRS(prj))
(nv <- mesh$n) # number of mesh nodes

# Plot to check
p1 <- ggplot() +
  theme_bw() + ylab("") + xlab("") +
  inlabru::gg(mesh) + 
  geom_sf(aes(), data = st_as_sf(dat, coords = c("lon", "lat"))
          %>% st_set_crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) +
  coord_sf(xlim = c(-4000, 3000), ylim = c(-4000, 3000), crs = prj, expand = F) +
  ggtitle("Haakon barrier model")
print(p1)

# Follow Haakon's example script to set up inla model
# Remove points on land
in.water = over(b, SpatialPoints(cbind(dat$x, dat$y), proj4string = CRS(prj)), returnList=T)[[1]]
print(paste("There are", nrow(dat)-length(in.water), "points on land in the original polygon"))
dat <- dat[in.water,]

# Set up boundary matern model
tl = length(mesh$graph$tv[,1]) # the number of triangles in the mesh

# Compute triangle positions
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

# Create a 1d time mesh for the annual cycle
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

#############################
### Set up year replicate ###
#############################

# The dataframe needs to be set up in such a way that all year x season interactions are included
# Make sure there are NAs not 0s when no sampling occured
# Currently we have nothing for year 4 - 2010:2014

# Assign points to each polygon
area <- factor(over(SpatialPoints(cbind(dat$x, dat$y)), polys),
               levels = 1:length(polys))

# Replicate dataframe to match number of time indices
# Add back the year and season index - there will be a better way to do this...
dat$year_i <- NA
dat$year_i[dat$index %in% c(1:4)] <- 1
dat$year_i[dat$index %in% c(5:8)] <- 2
dat$year_i[dat$index %in% c(9:12)] <- 3
dat$year_i[dat$index %in% c(13:16)] <- 4
dat$year_i[dat$index %in% c(17:20)] <- 5

dat$season <- NA
dat$season[dat$index %in% c(1, 5, 9, 13, 17)] <- 1
dat$season[dat$index %in% c(2, 6, 10, 14, 18)] <- 2
dat$season[dat$index %in% c(3, 7, 11, 15, 19)] <- 3
dat$season[dat$index %in% c(4, 8, 12, 16, 20)] <- 4

# Load in the population reference table
pop_ref <- readRDS("populations.rds")

# Match the ids with pop ref and take the corresponding population reference
dat <- dat %>% left_join(pop_ref, by = c("id" = "ref"))
dat <- dat %>% select("id", "date", "lat", "lon", "x", "y", "year", "month", "index", "year_i", "season", "location")

dat <- dat %>% mutate(population = case_when(location == "Newfoundland" ~ 1,
                                             location == "West Ice" ~ 2,
                                             location == "East Ice" ~ 3))

<<<<<<< HEAD
ggplot() + geom_point(aes(x = lon, y = lat, colour = factor(population)), data = dat)

# How to standardise effort by number of tags deployed...?
dat %>% group_by(index, population) %>% summarise(n_distinct(id))
# How to standardise effort by number of locations...?
dat %>% group_by(index, population) %>% summarise(n())
=======
ggplot() + geom_point(aes(x = lon, y = lat, colour = population), data = dat)

>>>>>>> be8b41563369880aac6a1e1fa0371b909eb55aec

# what is missing?
table(dat$index)
table(dat$season)
table(dat$year_i)
table(dat$population)

season <- dat$season
year <- dat$year_i
population <- dat$population

# Use both space and time index to aggregate data
agg.dat <- as.data.frame(table(area, season, population, year))

for(j in 1:4) # set season, year and area as integer
  agg.dat[[j]] <- as.integer(as.character(agg.dat[[j]])) 
str(agg.dat)

# INLA input dataframe now has mesh nodes replicated 48 times (16 * 3) - we have no data for 2010-14
# All season:year combinations are represented
# It is easy to then add the correct space-time covariates

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
e0 <- w.areas[agg.dat$area] #* (w.t[agg.dat$time])
summary(e0)

# Create projector matrix
A.st <- inla.spde.make.A(mesh = smesh,
                         loc = smesh$loc[agg.dat$area, ],
                         group = agg.dat$season,
                         n.group = k,
                         mesh.group = tmesh,
                         repl = agg.dat$population,
                         n.repl = 3)

# Create space-time index
idx <- inla.spde.make.index(name = 's',
                            n.spde = barrier.model$f$n,
                            n.group = k,
                            n.repl = 3)

# Define the data stack for the inla model
stk <- inla.stack(
  data = list(y = agg.dat$Freq, exposure = e0), 
  A = list(A.st, 1), 
  effects = list(idx, list(b0 = rep(1, nrow(agg.dat)))))

# PC prior on temporal correlation
pcrho <- list(theta = list(prior = 'pccor1', param = c(.7, .7))) # order is mu and alpha (1/sd^2)

<<<<<<< HEAD
=======
# Model formula
#form_1 <- y ~ 0 + b0 + f(s, model = barrier.model, group = s.group, control.group = list(model = 'ar1', hyper = pcrho))

# Fit the model
# NB this will take several hours...
#m_1 <- inla(form_1,
#            family = 'poisson',
#            data = inla.stack.data(stk),
#            control.predictor = list(A = inla.stack.A(stk)),
#            E = exposure,
#            control.inla = list(int.strategy = "eb"), # strategy ='adaptive' fails
#            control.compute = list(config = TRUE,
#                                   dic = T,
#                                   waic = T))


# IT WORKS!!! AND IN 2.3 HOURS!!! THATS with half the mesh though... does that matter?

###
### now to include ice...
###

>>>>>>> be8b41563369880aac6a1e1fa0371b909eb55aec
#append the sea ice data to the agg.dat
# the old code was based on the 1:20 index
agg.dat$index <- NA

agg.dat$index[agg.dat$year == 1 & agg.dat$season == 1] <- 1
agg.dat$index[agg.dat$year == 1 & agg.dat$season == 2] <- 2
agg.dat$index[agg.dat$year == 1 & agg.dat$season == 3] <- 3
agg.dat$index[agg.dat$year == 1 & agg.dat$season == 4] <- 4

agg.dat$index[agg.dat$year == 2 & agg.dat$season == 1] <- 5
agg.dat$index[agg.dat$year == 2 & agg.dat$season == 2] <- 6
agg.dat$index[agg.dat$year == 2 & agg.dat$season == 3] <- 7
agg.dat$index[agg.dat$year == 2 & agg.dat$season == 4] <- 8

agg.dat$index[agg.dat$year == 3 & agg.dat$season == 1] <- 9
agg.dat$index[agg.dat$year == 3 & agg.dat$season == 2] <- 10
agg.dat$index[agg.dat$year == 3 & agg.dat$season == 3] <- 11
agg.dat$index[agg.dat$year == 3 & agg.dat$season == 4] <- 12

agg.dat$index[agg.dat$year == 4 & agg.dat$season == 1] <- 13
agg.dat$index[agg.dat$year == 4 & agg.dat$season == 2] <- 14
agg.dat$index[agg.dat$year == 4 & agg.dat$season == 3] <- 15
agg.dat$index[agg.dat$year == 4 & agg.dat$season == 4] <- 16

agg.dat$index[agg.dat$year == 5 & agg.dat$season == 1] <- 17
agg.dat$index[agg.dat$year == 5 & agg.dat$season == 2] <- 18
agg.dat$index[agg.dat$year == 5 & agg.dat$season == 3] <- 19
agg.dat$index[agg.dat$year == 5 & agg.dat$season == 4] <- 20

#stack all 20 raster layers together
ice <- raster::stack(raster::stack("NSIDC Sea Ice/Seasonal_NSIDC_Dec94_Nov99"),
                     raster::stack("NSIDC Sea Ice/Seasonal_NSIDC_Dec99_Nov04"),
                     raster::stack("NSIDC Sea Ice/Seasonal_NSIDC_Dec04_Nov09"),
                     raster::stack("NSIDC Sea Ice/Seasonal_NSIDC_Dec09_Nov14"),
                     raster::stack("NSIDC Sea Ice/Seasonal_NSIDC_Dec14_Nov19"))

# then extract only data needed by unique times..
agg.dat$ice <- NA
for (i in sort(unique(agg.dat$index))){
  agg.dat$ice[agg.dat$index == i] <- raster::extract(raster::subset(ice, i),
                                                     ips %>% st_transform("+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +a=6378273 +b=6356889.449 +units=m +no_defs"),
                                                     method = "bilinear")
}

# average across years
av_ice <- raster::stackApply(ice, indices = rep(1:4, times = 5), fun = mean, na.rm = T)

# append seasonal average to agg.dat
unique(agg.dat$index)

for (i in 1:4){
  av_ice_vals <- raster::extract(raster::subset(av_ice, i),
                                 ips %>% st_transform("+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +a=6378273 +b=6356889.449 +units=m +no_defs"),
                                 method = "bilinear")
  agg.dat$av_ice[agg.dat$index %in% seq(i, 20, by = 4)] <- rep(av_ice_vals,
                                                               times = length(unique(agg.dat$index[agg.dat$index %in% seq(i, 20, by = 4)])))
}

# Define the data stack for the inla model
stk <- inla.stack(
  data = list(y = agg.dat$Freq, exposure = e0), 
  A = list(A.st, 1), 
  effects = list(idx, list(b0 = rep(1, nrow(agg.dat)),
                           ice = agg.dat$ice,
                           ice_av = agg.dat$av_ice,
                           ice_dev = agg.dat$ice - agg.dat$av_ice)))

# PC prior on temporal correlation
pcrho <- list(theta = list(prior = 'pccor1', param = c(.7, .7))) # order is mu and alpha (1/sd^2)

# Model formula
# Model failed when passing all depth data - so group it instead...
# see ?inla.group for more info

# 'null' model of season space use only
f_1 <- y ~ 0 + b0 + f(s, model = barrier.model, group = s.group, control.group = list(model = 'ar1', hyper = pcrho))

# are the seals responding to the ice they see in that year?
# model this as the additive effect of average ice and year specific deviation from average
f_2 <- y ~ 0 + b0 +
  f(s, model = barrier.model, group = s.group, control.group = list(model = 'ar1', hyper = pcrho)) +
  f(inla.group(ice_av, n = 100, method = "cut"), model = 'rw2', scale.model = T, hyper = list(theta = list(prior = "pc.prec", param = c(120, 0.01)))) +
  f(inla.group(ice_dev, n = 100, method = "cut"), model = 'rw2', scale.model = T, hyper = list(theta = list(prior = "pc.prec", param = c(120, 0.01))))

<<<<<<< HEAD
# can we replicate spatial field across populations?]
f_3 <- y ~ 0 + b0 +
  f(s, model = barrier.model, group = s.group, replicate = s.repl, control.group = list(model = 'ar1', hyper = pcrho)) +
  f(inla.group(ice_av, n = 100, method = "cut"), model = 'rw2', scale.model = T, hyper = list(theta = list(prior = "pc.prec", param = c(120, 0.01)))) +
  f(inla.group(ice_dev, n = 100, method = "cut"), model = 'rw2', scale.model = T, hyper = list(theta = list(prior = "pc.prec", param = c(120, 0.01))))

# Is it better just to model the actual ice - not as seasonal variations?
#f_3 <- y ~ 0 + b0 +
#  f(s, model = barrier.model, group = s.group control.group = list(model = 'ar1', hyper = pcrho)) +
#  f(inla.group(ice, n = 100, method = "cut"), model = 'rw2', scale.model = T, hyper = list(theta = list(prior = "pc.prec", param = c(120, 0.01))))
=======
# Is it better just to model the actual ice - not as seasonal variations?
f_3 <- y ~ 0 + b0 +
  f(s, model = barrier.model, group = s.group, control.group = list(model = 'ar1', hyper = pcrho)) +
  f(inla.group(ice, n = 100, method = "cut"), model = 'rw2', scale.model = T, hyper = list(theta = list(prior = "pc.prec", param = c(120, 0.01))))

>>>>>>> be8b41563369880aac6a1e1fa0371b909eb55aec

# previous knowledge
start_vals <- readRDS("start_vals.rds")

# Fit the model in 4 hours
m_2 <- inla(f_2,
            family = 'poisson', 
            data = inla.stack.data(stk),
            control.predictor = list(A = inla.stack.A(stk)),
            E = exposure,
            control.inla = list(int.strategy = "eb"), # strategy ='adaptive' fails
            control.compute = list(config = TRUE,
                                   dic = T,
                                   waic = T),
            control.mode = list(theta = start_vals, restart = TRUE),
            verbose = T) # switch on when trialling

<<<<<<< HEAD
# Fit the model in 32 hours
=======
# Fit the model in 4 hours
>>>>>>> be8b41563369880aac6a1e1fa0371b909eb55aec
m_3 <- inla(f_3,
            family = 'poisson', 
            data = inla.stack.data(stk),
            control.predictor = list(A = inla.stack.A(stk)),
            E = exposure,
            control.inla = list(int.strategy = "eb"), # strategy ='adaptive' fails
            control.compute = list(config = TRUE,
                                   dic = T,
                                   waic = T),
            #            control.mode = list(theta = start_vals, restart = TRUE),
            verbose = T) # switch on when trialling


require(animation)
saveVideo({
  for (i in 1:4){
    p1 <- ggplot() +
      theme_bw() + ylab("") + xlab("") +
<<<<<<< HEAD
      gg(smesh, col = m_2$summary.random$s$mean[(nv+1):(nv*2)]) +
=======
      gg(mesh, col = m_2$summary.random$s$mean[idx$s.group == i]) +
>>>>>>> be8b41563369880aac6a1e1fa0371b909eb55aec
      scale_fill_viridis("", limits = c(-5, 15), breaks = seq(-5, 15, 5), na.value = "transparent") +
      geom_sf(aes(), fill = "grey", colour = "grey", data = land) +
      coord_sf(xlim = c(-4000, 3000), ylim = c(-4000, 3000), crs = prj, expand = F) +
      ggtitle(i)
    print(p1)
  }
}, movie.name = "INLA_harps_grouped.mp4", interval = 1, ani.width = 750, ani.height = 750, other.opts = "-pix_fmt yuv420p -b:v 1080k")




p1 <- ggplot() +
  theme_bw() + ylab("") + xlab("") +
  gg(mesh, col = m_3$summary.random$s$mean[idx$s.group == 1]) +
  scale_fill_viridis("", limits = c(-5, 25), breaks = seq(-5, 25, 5), na.value = "transparent") +
  geom_sf(aes(), fill = "grey", colour = "grey", data = land) +
  coord_sf(xlim = c(-4000, 3000), ylim = c(-4000, 3000), crs = prj, expand = F) +
  ggtitle(1)
p2 <- ggplot() +
  theme_bw() + ylab("") + xlab("") +
  gg(mesh, col = m_3$summary.random$s$mean[idx$s.group == 2]) +
  scale_fill_viridis("", limits = c(-5, 25), breaks = seq(-5, 25, 5), na.value = "transparent") +
  geom_sf(aes(), fill = "grey", colour = "grey", data = land) +
  coord_sf(xlim = c(-4000, 3000), ylim = c(-4000, 3000), crs = prj, expand = F) +
  ggtitle(2)
p3 <- ggplot() +
  theme_bw() + ylab("") + xlab("") +
  gg(mesh, col = m_3$summary.random$s$mean[idx$s.group == 3]) +
  scale_fill_viridis("", limits = c(-5, 25), breaks = seq(-5, 25, 5), na.value = "transparent") +
  geom_sf(aes(), fill = "grey", colour = "grey", data = land) +
  coord_sf(xlim = c(-4000, 3000), ylim = c(-4000, 3000), crs = prj, expand = F) +
  ggtitle(3)
p4 <- ggplot() +
  theme_bw() + ylab("") + xlab("") +
  gg(mesh, col = m_3$summary.random$s$mean[idx$s.group == 4]) +
  scale_fill_viridis("", limits = c(-5, 25), breaks = seq(-5, 25, 5), na.value = "transparent") +
  geom_sf(aes(), fill = "grey", colour = "grey", data = land) +
  coord_sf(xlim = c(-4000, 3000), ylim = c(-4000, 3000), crs = prj, expand = F) +
  ggtitle(4)

quartz(width = 10, height = 10)
gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 2)
quartz.save("Replicate year spatial field.jpeg",
            type = "jpeg",
            dev = dev.cur(),
            dpi = 500)
dev.off()


pals <- tibble(group = idx$s.group,
               vals = m_2$summary.random$s$mean)
ggplot() +
  geom_sf(aes(colour = pals$vals[pals$group == 3]), data = ips)

ggplot() +  
  gg(mesh, col = pals$vals[pals$group == 1]) +
  facet_wrap(~pals$group)



ggplot() +
  theme_bw() + ylab("") + xlab("") +
  inlabru::gg(mesh) + 
  geom_sf(aes(), data = st_as_sf(dat, coords = c("lon", "lat"))
          %>% st_set_crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) +
  coord_sf(xlim = c(-4000, 3000), ylim = c(-4000, 3000), crs = prj, expand = F) +
  facet_wrap(~index)



require(viridis)

ggplot() +
  theme_bw() + ylab("") + xlab("") +
  gg(mesh) +
  geom_sf(aes(), data = ips)

gg(mesh, col = m_2$summary.random$s$mean[idx$s.group == 1]) + 
  scale_fill_viridis("", limits = c(-5, 15), breaks = seq(-5, 15, 5), na.value = "transparent") +
  geom_sf(aes(), fill = "grey", colour = "grey", data = land) +
  coord_sf(xlim = c(-4000, 3000), ylim = c(-4000, 3000), crs = prj, expand = F)



p1 <- ggplot() +
<<<<<<< HEAD
  geom_ribbon(aes(x = ID, ymin = exp(`0.025quant`), ymax = exp(`0.975quant`)), data = m_2$summary.random[[2]], alpha = 0.3) +
  geom_line(aes(x = ID, y = exp(mean)), data = m_2$summary.random[[2]]) +
  xlab("Sea Ice Concentration (%)")

p2 <- ggplot() +
  geom_ribbon(aes(x = ID, ymin = exp(`0.025quant`), ymax = exp(`0.975quant`)), data = m_3$summary.random[[3]], alpha = 0.3) +
  geom_line(aes(x = ID, y = exp(mean)), data = m_3$summary.random[[3]]) +
=======
  geom_ribbon(aes(x = ID, ymin = exp(`0.025quant`), ymax = exp(`0.975quant`)), data = m_3$summary.random[[2]], alpha = 0.3) +
  geom_line(aes(x = ID, y = exp(mean)), data = m_3$summary.random[[2]]) +
  xlab("Sea Ice Concentration (%)")

p2 <- ggplot() +
  geom_ribbon(aes(x = ID, ymin = exp(`0.025quant`), ymax = exp(`0.975quant`)), data = m_2$summary.random[[3]], alpha = 0.3) +
  geom_line(aes(x = ID, y = exp(mean)), data = m_2$summary.random[[3]]) +
>>>>>>> be8b41563369880aac6a1e1fa0371b909eb55aec
  xlab("Deviation from seasonal average sea ice concentration (%)") +
  coord_cartesian(ylim = c(0, 20))


ggplot() +
  geom_line(aes(x = ID, y = exp(mean)), data = cbind.data.frame(m_2$summary.random[[2]], m_2$summary.random[[3]])) +
  xlab("Deviation from seasonal average sea ice concentration (%)") +
  coord_cartesian(ylim = c(0, 20))

m_2$summary.random[[2]]$mean + m_2$summary.random[[3]]$mean


quartz(width = 10, height = 6)
gridExtra::grid.arrange(p1, p2, ncol = 2)
quartz.save("2 component ice output.jpeg",
            type = "jpeg",
            dev = dev.cur(),
            dpi = 500)
dev.off()


foo <- cbind.data.frame(ips %>% st_coordinates(),
                        as.matrix(A.st) %*% m_3$summary.random$s$mean)
names(foo) <- c("x", "y", "z")
ggplot() + 
  geom_point(aes(x = x, y = y, colour = z), data = foo)
su

# Marginals
foo <- m_2$marginals.hyperpar %>% map_dfr(~ .x %>% as_tibble(), .id = "name")
ggplot() + 
  geom_line(aes(x = x, y = y), data = foo) +
  ylab("Density") + xlab("") +
  facet_wrap(~ name, ncol = 1, scales = "free")




foo <- inla.posterior.sample(n = 1000, result = m_2, use.improved.mean = T)
effect <- "APredictor"
id.effect <- which(m_2$misc$configs$contents$tag == effect)
ind.effect <- m_2$misc$configs$contents$start[id.effect] - 1 + (1:m_2$misc$configs$contents$length[id.effect])
Predictor_sample <- lapply(foo, function(x) x$latent[ind.effect])

# this is a list with number of elements equal to number of iterations
# how do you average across the iterations?

preds <- Predictor_sample %>%
  tibble() %>%
  unnest(.id = "name") %>%
  group_by(name) %>%
  mutate(row = 1:n()) %>%
  ungroup() %>%
  group_by(row) %>%
  summarise(mean = mean(., na.rm = T),
            sd = sd(., na.rm = T)) %>%
  mutate(group = rep(1:16, each = nv))

require(animation)
saveVideo({
  for (i in 1:16){
    p1 <- ggplot() +
      theme_bw() + ylab("") + xlab("") +
      gg(mesh, col = exp(preds$mean[preds$group == i])) + 
      scale_fill_viridis("", na.value = "transparent") +
      #  scale_fill_viridis("", limits = c(-5, 15), breaks = seq(-5, 15, 5), na.value = "transparent") +
      geom_sf(aes(), fill = "grey", colour = "grey", data = land) +
      coord_sf(xlim = c(-4000, 3000), ylim = c(-4000, 3000), crs = prj, expand = F) +
      ggtitle(i)
    print(p1)
  }
}, movie.name = "INLA_harps_predicted.mp4", interval = 1, ani.width = 750, ani.height = 750, other.opts = "-pix_fmt yuv420p -b:v 1080k")


<<<<<<< HEAD

# Use Finn's nice interpolation in the gg function to generate faceted plot for each season x population
foo <- bind_rows(
  gg(mesh, col = m_3$summary.random$s$mean[idx$s.repl == 1 & idx$s.group == 1])[[1]]$data %>%
    as_tibble() %>% mutate(population = "NW Atlantic",
                           season = "Dec-Feb"),
  gg(mesh, col = m_3$summary.random$s$mean[idx$s.repl == 1 & idx$s.group == 2])[[1]]$data %>%
    as_tibble() %>% mutate(population = "NW Atlantic",
                           season = "Mar-May"),
  gg(mesh, col = m_3$summary.random$s$mean[idx$s.repl == 1 & idx$s.group == 3])[[1]]$data %>%
    as_tibble() %>% mutate(population = "NW Atlantic",
                           season = "Jun-Aug"),
  gg(mesh, col = m_3$summary.random$s$mean[idx$s.repl == 1 & idx$s.group == 4])[[1]]$data %>%
    as_tibble() %>% mutate(population = "NW Atlantic",
                           season = "Sep-Nov"),
  gg(mesh, col = m_3$summary.random$s$mean[idx$s.repl == 2 & idx$s.group == 1])[[1]]$data %>%
    as_tibble() %>% mutate(population = "West Ice",
                           season = "Dec-Feb"),
  gg(mesh, col = m_3$summary.random$s$mean[idx$s.repl == 2 & idx$s.group == 2])[[1]]$data %>%
    as_tibble() %>% mutate(population = "West Ice",
                           season = "Mar-May"),
  gg(mesh, col = m_3$summary.random$s$mean[idx$s.repl == 2 & idx$s.group == 3])[[1]]$data %>%
    as_tibble() %>% mutate(population = "West Ice",
                           season = "Jun-Aug"),
  gg(mesh, col = m_3$summary.random$s$mean[idx$s.repl == 2 & idx$s.group == 4])[[1]]$data %>%
    as_tibble() %>% mutate(population = "West Ice",
                           season = "Sep-Nov"),
  gg(mesh, col = m_3$summary.random$s$mean[idx$s.repl == 3 & idx$s.group == 1])[[1]]$data %>%
    as_tibble() %>% mutate(population = "East Ice",
                           season = "Dec-Feb"),
  gg(mesh, col = m_3$summary.random$s$mean[idx$s.repl == 3 & idx$s.group == 2])[[1]]$data %>%
    as_tibble() %>% mutate(population = "East Ice",
                           season = "Mar-May"),
  gg(mesh, col = m_3$summary.random$s$mean[idx$s.repl == 3 & idx$s.group == 3])[[1]]$data %>%
    as_tibble() %>% mutate(population = "East Ice",
                           season = "Jun-Aug"),
  gg(mesh, col = m_3$summary.random$s$mean[idx$s.repl == 3 & idx$s.group == 4])[[1]]$data %>%
    as_tibble() %>% mutate(population = "East Ice",
                           season = "Sep-Nov")
  )

foo$season <- factor(foo$season)
foo$season <- fct_relevel(foo$season, "Dec-Feb", "Mar-May")
foo$population <- factor(foo$population)
foo$population <- fct_relevel(foo$population, "NW Atlantic", "West Ice")
source("discrete_gradient.r")

p <- ggplot() +
  theme_minimal() +
  ylab("") + xlab("") +
  geom_point(aes(x = x, y = y, colour = color), data = foo) +
  geom_tile(aes(x = x, y = y, fill = color), data = foo) +
  scale_colour_discrete_gradient("Posterior mean estimate of spatial field",
                                 colours = viridis::viridis(14),
                                 bins = 14,
                                 limits = c(-7.5, 27.5),
                                 breaks = seq(-5, 25, 5),
                                 guide = guide_colourbar(nbin = 500,
                                                         raster = T,
                                                         frame.colour = "black",
                                                         ticks.colour = "black",
                                                         frame.linewidth = 1,
                                                         barwidth = 1,
                                                         barheight = 25,
                                                         direction = "vertical",
                                                         title.position = "right",
                                                         title.theme = element_text(angle = 90,
                                                                                    hjust = 0.5))) +
  geom_sf(aes(), fill = "grey", colour = "grey", data = land) +
  coord_sf(xlim = c(-4000, 3000), ylim = c(-4000, 3000), crs = prj, expand = F) +
  facet_grid(season ~ population)


quartz(width = 9, height = 12)
print(p)
quartz.save(file = "~/Desktop/posterior mean by season and population.jpeg",
            type = "jpeg",
            dev = dev.cur(),
            dpi = 500)
dev.off()

# Hyperparameter posteriors look better too
hype <- bind_rows(
  m_2$marginals.hyperpar[[1]] %>% as_tibble() %>% mutate(group = names(m_2$marginals.hyperpar)[[1]],
                                                       model = "ignoring population") %>% filter(x > 0) %>% filter(x < 10),
  m_2$marginals.hyperpar[[2]] %>% as_tibble() %>% mutate(group = names(m_2$marginals.hyperpar)[[2]],
                                                       model = "ignoring population") %>% filter(x > 0) %>% filter(x < 10),
  m_2$marginals.hyperpar[[3]] %>% as_tibble() %>% mutate(group = names(m_2$marginals.hyperpar)[[3]],
                                                       model = "ignoring population") %>% filter(x > 0) %>% filter(x < 10),
  m_2$marginals.hyperpar[[4]] %>% as_tibble() %>% mutate(group = names(m_2$marginals.hyperpar)[[4]],
                                                       model = "ignoring population") %>% filter(x > 0) %>% filter(x < 10),
  m_2$marginals.hyperpar[[5]] %>% as_tibble() %>% mutate(group = names(m_2$marginals.hyperpar)[[5]],
                                                       model = "ignoring population") %>% filter(x > 0) %>% filter(x < 10),
  m_3$marginals.hyperpar[[1]] %>% as_tibble() %>% mutate(group = names(m_3$marginals.hyperpar)[[1]],
                                                       model = "population replicates") %>% filter(x > 0) %>% filter(x < 10),
  m_3$marginals.hyperpar[[2]] %>% as_tibble() %>% mutate(group = names(m_3$marginals.hyperpar)[[2]],
                                                       model = "population replicates") %>% filter(x > 0) %>% filter(x < 10),
  m_3$marginals.hyperpar[[3]] %>% as_tibble() %>% mutate(group = names(m_3$marginals.hyperpar)[[3]],
                                                       model = "population replicates") %>% filter(x > 0) %>% filter(x < 10),
  m_3$marginals.hyperpar[[4]] %>% as_tibble() %>% mutate(group = names(m_3$marginals.hyperpar)[[4]],
                                                       model = "population replicates") %>% filter(x > 0) %>% filter(x < 10),
  m_3$marginals.hyperpar[[5]] %>% as_tibble() %>% mutate(group = names(m_3$marginals.hyperpar)[[5]],
                                                       model = "population replicates") %>% filter(x > 0) %>% filter(x < 10)
)
=======
preds$season <- agg.dat$season
preds$season <- agg.dat$season
preds$year <- agg.dat$year


preds %>%
  group_by(season) %>%
  slice(1:2)








Predictor_sample %>%
  tibble()



%>%
  glimpse()

mean(Predictor_sample[[1]][1], Predictor_sample[[2]][1])



Predictor_sample %>%
  purrr::transpose() %>%
  purrr::map_df(~rowMeans(as.data.frame(.x), na.rm = T))

zoo <- plyr::aaply(plyr::laply(Predictor_sample, as.matrix), c(2,3), mean)

purr:map_df()

Predictor_sample %>% tibble() %>% map_df(., function(df){summarise_all(df, mean)})

map_dfr(.)



poo <- 
  
  Predictor_sample %>%
  tibble() %>%
  rowwise() %>%
  summarise(mean(.))
mutate(mean = map_dbl(., mean))

unnest(.id = "name") %>%
  group_by(name) %>%
  summarise(Mean = mean(.))

Predictor_sample %>% map_dfr(~ .x %>% as_tibble(), .id = "name")

>>>>>>> be8b41563369880aac6a1e1fa0371b909eb55aec



ggplot() +
<<<<<<< HEAD
  geom_line(aes(x = x, y = y, group = group), data = hype) +
  facet_wrap( ~ group + model,
              scales = "free",
              ncol = 2)
=======
  theme_bw() + ylab("") + xlab("") +
  gg(mesh, col = Predictor_sample[[1]][1:nv])



+ 
  scale_fill_viridis("", limits = c(-5, 15), breaks = seq(-5, 15, 5), na.value = "transparent") +
  geom_sf(aes(), fill = "grey", colour = "grey", data = land) +
  coord_sf(xlim = c(-4000, 3000), ylim = c(-4000, 3000), crs = prj, expand = F)



>>>>>>> be8b41563369880aac6a1e1fa0371b909eb55aec


