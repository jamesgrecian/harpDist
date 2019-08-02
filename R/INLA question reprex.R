########################################
### Space-time LGCP indexing in INLA ###
########################################

# Generate a point process that varies in time and space ###
# The data are from three populations and five time periods
# following the SPDE book:
# https://becarioprecario.bitbucket.io/spde-gitbook/ch-stapp.html

# Load libraries
require(tidyverse)
require(INLA)
require(inlabru)
require(rgeos)
require(sf)

# Generate some random points
pp <- spatstat::rpoispp(1200) %>% as_tibble()

# Divide into three populations seperated on the x axis
pp <- pp %>% mutate(population = case_when(x < 1/3 ~ 1,
                                      x > 1/3 & x < 2/3 ~ 2,
                                      x > 2/3  ~ 3))
# Divide into four seasons based on y axis movements 
pp <- pp %>% mutate(season = case_when(y < 1/3 ~ 4,
                                       y > 1/3 & y < 2/3 ~ 1,
                                       y > 2/3  ~ 2))
# Fudge to add autumn season
pp <- rbind(pp,
            pp %>% filter(season == 1) %>% mutate(season = 3,
                                                         x = jitter(x, amount = .05),
                                                         y = jitter(y, amount = .05)))
# Seperate into 5 equal sampling years
pp <- pp %>% group_by(population, season) %>% slice(1:100)
pp <- pp %>% mutate(year = rep(1:5, each = 20))

# Plot the data to illustrate the sampling design
ggplot() +
  geom_point(aes(x = x, y = y, colour = factor(year)), data = pp) +
  facet_wrap(~ population + season)

#############################
### Set up the INLA model ###
#############################
# Define the domain
b <- list(inla.nonconvex.hull(as.matrix(pp[1:2])))

# Create the spatial mesh
s_mesh <- inla.mesh.2d(pp[1:2], boundary = b, max.edge=c(1, 0.5))

# Define spde
spde = inla.spde2.pcmatern(s_mesh, prior.range = c(.05, .01), prior.sigma = c(1, 0.01))

# Plot to check
ggplot() +
  inlabru::gg(s_mesh) +
  geom_point(aes(x = x, y = y), data = pp)

#How many mesh nodes and how many points
nv <- s_mesh$n
n <- nrow(pp)
ips <- ipoints(s_mesh)

# Create a 1d time mesh for the annual cycle
# This can be seasonal (1-4) by ignoring year
t_mesh <- inla.mesh.1d(loc = 1:4, boundary = "free")
t_mesh$loc
(k <- length(t_mesh$loc))

# Find Voroni polygons for mesh nodes
library(deldir)
dd <- deldir(s_mesh$loc[, 1], s_mesh$loc[, 2])
tiles <- tile.list(dd)

# Convert to Spatial Polygons
polys <- SpatialPolygons(lapply(1:length(tiles), function(i) {
  p <- cbind(tiles[[i]]$x, tiles[[i]]$y)
  n <- nrow(p)
  Polygons(list(Polygon(p[c(1:n, 1), ])), i)
}))

# Assign points to each polygon
area <- factor(over(SpatialPoints(cbind(pp$x, pp$y)), polys),
               levels = 1:length(polys))

# Replicate dataframe to match number of time indices
time <- pp$season
table(time)

# Use both space and time index to aggregate data
agg.dat <- as.data.frame(table(area, time))
for(j in 1:2) # set time and area as integer
  agg.dat[[j]] <- as.integer(as.character(agg.dat[[j]])) 
str(agg.dat)

# Expected number of points needs to be defined as
# proportional to the area of the polygon * the width of time knot
domain <- SpatialPolygons(list(Polygons(list(Polygon(as.matrix(b[[1]]$loc))), 1)))

w.areas <- sapply(1:length(tiles), function(i) {
  p <- cbind(tiles[[i]]$x, tiles[[i]]$y)
  n <- nrow(p)
  pl <- SpatialPolygons(
    list(Polygons(list(Polygon(p[c(1:n, 1),])), i)))
  if (gIntersects(pl, domain))
    return(gArea(gIntersection(pl, domain)))
  else return(0)
})

# Summary of the polygon areas
summary(w.areas)
sum(w.areas)

# Check with total area
(s.area <- gArea(domain))

# Width of each knot
(w.t <- diag(inla.mesh.fem(t_mesh)$c0))

# Estimate intensity based on domain and time
i0 <- n / (gArea(domain) * diff(range(t_mesh$loc)))
c(i0, log(i0))

# The space-time volume at each polygon and time knot
e0 <- w.areas[agg.dat$area] * (w.t[agg.dat$time])
summary(e0)

# Create projector matrix
A.st <- inla.spde.make.A(mesh = s_mesh,
                         loc = s_mesh$loc[agg.dat$area, ],
                         group = agg.dat$time,
                         mesh.group = t_mesh)

# Create space-time index
idx <- inla.spde.make.index(name = 's',
                            n.spde = spde$f$n,
                            n.group = k)

# Define the data stack for the inla model
stk <- inla.stack(
  data = list(y = agg.dat$Freq, exposure = e0), 
  A = list(A.st, 1), 
  effects = list(idx, list(b0 = rep(1, nrow(agg.dat)))))

# PC prior on temporal correlation
pcrho <- list(theta = list(prior = 'pccor1', param = c(.3, .7))) # order is mu and alpha (1/sd^2)

# Model formula
form_1 <- y ~ 0 + b0 + f(s, model = spde, group = s.group, control.group = list(model = 'ar1', hyper = pcrho))

# Fit the model
# NB this will take about 8 minutes...
m_1 <- inla(form_1,
            family = 'poisson',
            data = inla.stack.data(stk),
            control.predictor = list(A = inla.stack.A(stk)),
            E = exposure,
            control.inla = list(int.strategy = "eb"),
            control.compute = list(config = TRUE,
                                   dic = T,
                                   waic = T))

p1 <- ggplot() +
  geom_sf(aes(colour = m_1$summary.random$s$mean[idx$s.group == 1]), data = sf::st_as_sf(ips)) +
  scale_color_viridis_c("") + ggtitle("1")
p2 <- ggplot() +
  geom_sf(aes(colour = m_1$summary.random$s$mean[idx$s.group == 2]), data = sf::st_as_sf(ips)) +
  scale_color_viridis_c("") + ggtitle("2")
p3 <- ggplot() +
  geom_sf(aes(colour = m_1$summary.random$s$mean[idx$s.group == 3]), data = sf::st_as_sf(ips)) +
  scale_color_viridis_c("") + ggtitle("3")
p4 <- ggplot() +
  geom_sf(aes(colour = m_1$summary.random$s$mean[idx$s.group == 4]), data = sf::st_as_sf(ips)) +
  scale_color_viridis_c("") + ggtitle("4")

gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 2)

