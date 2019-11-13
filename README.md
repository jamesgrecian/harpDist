<!-- README.md is generated from README.Rmd. Please edit that file -->
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)

harpDist
========

**harpDist** is a space-time modelling framework to estimate the impact
of climate change on harp seal migratory behaviour.

Given a set of animal telemetry data, fit a spatio-temporal log Gaussian
Cox Process to estimate the distribution and the influence of
environmental covariates.

An example of how to fit a LGCP including the Bakka barrier model is
below:

An example analysis
-------------------

Require packages and download data from repo

``` r

# Load libraries
require(INLA)
require(tidyverse)
require(sf)
require(viridis)
require(raster)
require(mapr)
require(rgeos)
source("R/spde-book-functions.R")
source("R/discrete_gradient.r")

# Here is a random subsample of 2500 animal locations
dat <- readRDS("data/harps2500_indexed.rds")
```

Use the `mapr` package to generate a land shapefile and the boundary for
the INLA mesh

``` r

# Define an Albers equal area projection using kilometres as units
prj = "+proj=laea +lat_0=75 +lon_0=-25 +x_0=0 +y_0=0 +datum=WGS84 +units=km +no_defs +ellps=WGS84 +towgs84=0,0,0"

# Generate a land shapefile for generating a map later on
land <- mapr::mapr(dat,
                   prj,
                   buff = 2000)

# Generate simple boundary for the inla mesh
b <- mapr::meshr(dat,
                 prj,
                 buff = 500,
                 keep = 0.5,
                 Neumann = F) 
```

Use the boundary file to generate the INLA mesh, construct the Voroni
polygons and calculate weights

``` r

# Remove any animal locations that fall on land
in.water = over(b, SpatialPoints(cbind(dat$x, dat$y), proj4string = CRS(prj)), returnList=T)[[1]]
print(paste("There are", nrow(dat)-length(in.water), "points on land in the original polygon"))
#> [1] "There are 62 points on land in the original polygon"
dat <- dat[in.water,]

# Define the parameters of the boundary mesh
max.edge = max(c(diff(range(dat$x)), diff(range(dat$y))))/15
bound.outer = 1000

# Create an inla mesh
# Changing cutoff is an effective way to change number of mesh nodes
# Fewer mesh nodes ~ less computation time
mesh = inla.mesh.2d(boundary = b,
                    loc = cbind(dat$x, dat$y),
                    max.edge = c(1, 5) * max.edge,
                    cutoff = 50, # this has been 25 to match ice...
                    offset = c(max.edge, bound.outer),
                    crs = CRS(prj))

# Construct Voroni polygons
dmesh <- book.mesh.dual(mesh)
proj4string(dmesh) <- prj

# Calculate weights for polygons
# Set weights to 0 when polygon is outside study area
w <- sapply(1:length(dmesh), function(i){
  if(gIntersects(dmesh[i,], b))
    return(gArea(gIntersection(dmesh[i,], b)))
  else return(0)
})
```

First fit a lgcp with matérn covariance spde

``` r

# Set up a matern covariance on the 2D spde
spde = inla.spde2.pcmatern(mesh,
                           prior.range = c(50, .1),
                           prior.sigma = c(10, .01))

# Set up the lgcp data structure
nv <- mesh$n # number of mesh nodes
n <- nrow(dat) # number of data points
y.pp <- rep(0:1, c(nv, n)) # response (0 for all meshnodes, 1 for all data points)
e.pp <- c(w, rep(0, n))

# Create projector matrix
lmat <- inla.spde.make.A(mesh = mesh,
                         loc = cbind(dat$x, dat$y))
imat <- Diagonal(nv, rep(1, nv))
A.pp <- rbind(imat, lmat)

# Define the data stack for the inla model
stk <- inla.stack(data = list(y = y.pp,
                              e = e.pp),
                  A = list(1, A.pp),
                  effects = list(list(b0 = rep(1, nv + n)),
                                 list(s = 1:nv)))

fit <- inla(y ~ 0 + b0 + f(s, model = spde),
            family = "poisson",
            data = inla.stack.data(stk),
            control.predictor = list(A = inla.stack.A(stk)),
            E = inla.stack.data(stk)$e,
            control.inla = list(int.strategy = "eb"), # strategy ='adaptive' fails
            control.compute = list(config = TRUE,
                                   dic = T,
                                   waic = T))
```

This model will ignore coastlines, but marine mammals are unlikely to
travel across land!

Instead fit the barrier matérn model. For details see:
<a href="https://doi.org/10.1016/j.spasta.2019.01.002" class="uri">https://doi.org/10.1016/j.spasta.2019.01.002</a>

``` r

# Set up boundary matern model
tl = length(mesh$graph$tv[,1]) # the number of triangles in the mesh

# Compute triangle positions
posTri = matrix(0, tl, 2) 
for (t in 1:tl){
  temp = mesh$loc[mesh$graph$tv[t, ], ]
  posTri[t,] = colMeans(temp)[c(1,2)] 
}
posTri = SpatialPoints(posTri, proj4string = CRS(prj))

normal = unlist(over(b, posTri, returnList = T)) # check which mesh triangles are inside the normal area
barrier.triangles = setdiff(1:tl, normal)
poly.barrier = inla.barrier.polygon(mesh, barrier.triangles)

# create the matern object
barrier.model = inla.barrier.pcmatern(mesh,
                                      barrier.triangles = barrier.triangles,
                                      prior.range = c(50, .1),
                                      prior.sigma = c(10, 0.01))

# inla call is the same just with barrier model specified in the spde
barrier.fit <- inla(y ~ 0 + b0 + f(s, model = barrier.model),
               family = "poisson",
               data = inla.stack.data(stk),
               control.predictor = list(A = inla.stack.A(stk)),
               E = inla.stack.data(stk)$e,
               control.inla = list(int.strategy = "eb"), # strategy ='adaptive' fails
               control.compute = list(config = TRUE,
                                      dic = T,
                                      waic = T))
#> Warning in inla.model.properties.generic(inla.trim.family(model), (mm[names(mm) == : Model 'rgeneric' in section 'latent' is marked as 'experimental'; changes may appear at any time.
#>   Use this model with extra care!!! Further warnings are disabled.
```

Visualise the posterior spatial fields for the two lgcp models

``` r

# Generate a tidy plot
preds <- dmesh %>% st_as_sf()
preds <- rbind(preds, preds)
preds <- preds %>% mutate(model = rep(c("spde", "barrier"), each = nv),
                          posterior = c(fit$summary.random$s$mean,
                                        barrier.fit$summary.random$s$mean))
preds$model <- factor(preds$model)
preds$model = relevel(preds$model, ref = "spde")

p <- ggplot() +
  theme_minimal() +
  ylab("") + xlab("") +
  geom_sf(aes(colour = posterior,
              fill = posterior), data = preds) +
  scale_colour_discrete_gradient("Posterior mean estimate of spatial field",
                                 colours = viridis::viridis(14),
                                 bins = 14,
                                 limits = c(-5, 12.5),
                                 breaks = seq(-5, 12.5, 2.5),
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
  geom_sf(aes(), colour = "white", fill = NA, data = st_as_sf(b)) +
  facet_wrap(~ model)
print(p)
```

![](README-visualise%20the%20lgcp-1.png)
