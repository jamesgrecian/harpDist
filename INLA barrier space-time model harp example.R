###############################################
### Harp seal INLA barrier space-time model ###
###############################################

# The data are harp seal telemetry locations collated from studies going back to 1995
# We are interested in:
# 1. Estimating the seasonal north-south migration
# 2. Estimating how that migration may have changed over the last 25 years

# I have binned the data into seasons (Dec-Feb; Mar-May; Jun-Aug; Sep-Nov)
# and into year blocks (95-99; 00-04; 05-09; 10-14; 15-19)
# These are combined into an index 1:20 (4 seasons x 5 years)
# I will then be fitting covariates such as sea ice conditions that will vary across these 20 time indices

# However, I would like to share the seasonal space field over all years
# as we don't have all populations tracked in all years

# I am struggling to understand how to index the spatial field from 1:4 when the data are indexed 1:20
# Here is how I would run the INLA model if I were only looking at the seasonal effect only

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

# Here is a random subsample of 500 animal locations
dat <- readRDS("harp data/harps500_indexed.rds")

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

# Create an inla mesh
mesh = inla.mesh.2d(boundary = b,
                    loc = cbind(dat$x, dat$y),
                    max.edge = c(1, 5) * max.edge,
                    cutoff = 25,
                    offset = c(max.edge, bound.outer))

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
  effects = list(idx, list(b0 = rep(1, nrow(agg.dat)))))

# PC prior on temporal correlation
pcrho <- list(theta = list(prior = 'pccor1', param = c(.7, .7))) # order is mu and alpha (1/sd^2)

# Model formula
form_2 <- y ~ 0 + b0 + f(s, model = barrier.model, group = s.group, control.group = list(model = 'ar1', hyper = pcrho))

# Fit the model
m_2 <- inla(form_2,
            family = 'poisson',
            data = inla.stack.data(stk),
            control.predictor = list(A = inla.stack.A(stk)),
            E = exposure,
            control.inla = list(int.strategy = "eb"), # strategy ='adaptive' fails
            control.compute = list(config = TRUE,
                                   dic = T,
                                   waic = T))
print(m_2$dic$dic)

### This fails but if we run the model on a time mesh 1:20 it will run...

tmesh <- inla.mesh.1d(loc = 1:20, boundary = "free")
tmesh$loc
(k <- length(tmesh$loc))

# Create projector matrix
A.st <- inla.spde.make.A(mesh = smesh,
                         loc = smesh$loc[agg.dat$area, ],
                         group = agg.dat$time,
                         mesh.group = tmesh)

# Create space-time index
idx <- inla.spde.make.index(name = 's',
                            n.spde = barrier.model$f$n,
                            n.group = 20)

# Define the data stack for the inla model
stk <- inla.stack(
  data = list(y = agg.dat$Freq, exposure = e0), 
  A = list(A.st, 1), 
  effects = list(idx, list(b0 = rep(1, nrow(agg.dat)))))

# PC prior on temporal correlation
pcrho <- list(theta = list(prior = 'pccor1', param = c(.7, .7))) # order is mu and alpha (1/sd^2)

# Model formula
form_2 <- y ~ 0 + b0 + f(s, model = barrier.model, group = s.group, control.group = list(model = 'ar1', hyper = pcrho))

# Fit the model
m_2 <- inla(form_2,
            family = 'poisson',
            data = inla.stack.data(stk),
            control.predictor = list(A = inla.stack.A(stk)),
            E = exposure,
            control.inla = list(int.strategy = "eb"), # strategy ='adaptive' fails
            control.compute = list(config = TRUE,
                                   dic = T,
                                   waic = T))
print(m_2$dic$dic)




