<!-- README.md is generated from README.Rmd. Please edit that file -->
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)

harpDist
========

**harpDist** is a framework for fitting spatio-temporal log Gaussian Cox
Processes to animal telemetry data.

Given a set of locations, for example from a tagged marine animal,
`harpDist` will generate a global shapefile for mapping, assist with
`INLA` mesh generation and then fit an LGCP model with or without
physical barriers using `R-INLA`

An example of how to do this is below:

An example analysis
-------------------

Require packages, download data from GitHub repository and define
geographic projection

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

# Define an Albers equal area projection using kilometres as units
prj = "+proj=laea +lat_0=75 +lon_0=-25 +x_0=0 +y_0=0 +datum=WGS84 +units=km +no_defs +ellps=WGS84 +towgs84=0,0,0"
```

Use the `mapr` package
(<a href="https://github.com/jamesgrecian/mapr" class="uri">https://github.com/jamesgrecian/mapr</a>)
to generate a land shapefile and the boundary for the INLA mesh

``` r
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
                    cutoff = 50,
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

Non-barrier model
-----------------

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
            control.inla = list(int.strategy = "eb"),
            control.compute = list(config = TRUE,
                                   dic = T,
                                   waic = T))
```

Bakka barrier model
-------------------

Marine mammals are unlikely to travel across land, so fit a spatial
model that includes coastlines as physical barriers to the estimation.

For details see:
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

Plot model outputs
------------------

Visualise the posterior spatial fields for the two lgcp models

``` r
# Generate a tidy plot
preds <- dmesh %>% st_as_sf()
preds <- rbind(preds, preds) # double length for facet_wrap plotting by model type
preds <- preds %>% mutate(model = rep(c("spde", "barrier"), each = nv),
                          posterior = c(fit$summary.random$s$mean,
                                        barrier.fit$summary.random$s$mean))
preds$model <- factor(preds$model)
preds$model <- relevel(preds$model, ref = "spde")

p <- ggplot() +
  theme_bw() +
  ylab("") + xlab("") +
  geom_sf(aes(colour = posterior,
              fill = posterior), data = preds) +
  scale_colour_discrete_gradient("Posterior mean estimate of spatial field",
                                 colours = viridis::viridis(13),
                                 bins = 13,
                                 limits = c(-6.5, 13.0),
                                 breaks = seq(-6.5, 13.0, 3),
                                 guide = guide_colourbar(nbin = 500,
                                                         raster = T,
                                                         frame.colour = "black",
                                                         ticks.colour = "black",
                                                         frame.linewidth = 1,
                                                         barwidth = 20,
                                                         barheight = 1,
                                                         direction = "horizontal",
                                                         title.position = "top",
                                                         title.theme = element_text(angle = 0,
                                                                                    hjust = 0.5))) +
  theme(legend.position = "bottom") +
  geom_sf(aes(), colour = "white", fill = NA, data = st_as_sf(b)) +
  facet_wrap(~ model)
print(p)
```

![](README-visualise%20the%20lgcp-1.png)
