###############################################
### Subsample harp seal point data for INLA ###
###############################################

#load libraries
require(inlabru)
require(INLA)
require(tidyverse)
require(sf)
#devtools::install_github("jamesgrecian/mapr")
require(mapr)
require(RSRDL)
require(viridis)

#Define albers projection
prj = "+proj=laea +lat_0=75 +lon_0=-25 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"

##########################
### Load in point data ###
##########################

dat = get.SRDLpg(theDB = 'harps', theTable = 'diag', theFields = 'All', theDep = 'All', theRef = 'All',
                 theHost = 'localhost', thePort = 5432, theUser = 'postgres', thePwd = 'admin')

#Check format of dataframe
str(dat); summary(dat)
dat$ref = as.factor(dat$ref)
dat$d.date = as.POSIXct(dat$d.date, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
dat$lq = as.factor(dat$lq)
levels(dat$lq) = c("Z", "B", "A", "0", "1", "2", "3") #relabel from -3 to 3

#Resort by id and date
dat = arrange(dat, ref, d.date)

#Remove duplicates
dat = dat[!duplicated(dat$d.date),]

#Exclude animals that have fewer than 20 locations
dat = dat[dat$ref != "hp1-9276-04",]; dat$ref = factor(dat$ref)
dat = dat[dat$ref != "hp1-9278-04",]; dat$ref = factor(dat$ref)
dat = dat[dat$ref != "hp1c-296-96",]; dat$ref = factor(dat$ref)
dat = dat[dat$ref != "hp2-9325-04",]; dat$ref = factor(dat$ref)
dat = dat[dat$ref != "hp2b-08-99",]; dat$ref = factor(dat$ref)
dat = dat[dat$ref != "hp4-L517-17",]; dat$ref = factor(dat$ref)
dat = dat[dat$ref != "hp1b-01-95",]; dat$ref = factor(dat$ref)
dat = dat[dat$ref != "hp1b-02-95",]; dat$ref = factor(dat$ref)

#Remove animals where tag failed within 1 week
dat = dat[dat$ref != "hp2-9334-04",]; dat$ref = factor(dat$ref)

#Format the data for ssmTMB - same format as bsam
dat2 = data.frame(id = dat$ref, date = dat$d.date, lc = dat$lq, lon = dat$lon, lat = dat$lat)
#Drop the Z locations
dat2 = dat2[!dat2$lc == "Z",]
#Remove all locations below 40 degrees North
dat2 = dat2[dat2$lat > 40,]
#Remove all locations above 85 degrees North
dat2 = dat2[dat2$lat < 85,]
#Remove all locations below -100 degrees West
dat2 = dat2[dat2$lon > -100,]
#Remove all locations above 100 degrees East
dat2 = dat2[dat2$lon < 100,]

# Run the model through foieGras to error correct and regularise
require(foieGras)
fls <- foieGras::fit_ssm(dat2, model = "crw", time.step = 6)
pred <- foieGras::pluck(fls, "fitted") #this should probably be predicted... but will need tweaking for particular individuals...

#Project to something more useful...
locs = SpatialPointsDataFrame(coords = cbind(pred$lon, pred$lat), data = pred[,1:4], proj4string = CRS("+proj=longlat +datum=WGS84"))
locs_laea = spTransform(locs, CRS(prj))
locs_laea = as.data.frame(locs_laea)
names(locs_laea)[5] = "x"
names(locs_laea)[6] = "y"

# Output subsample...
dat3 <- sample_n(locs_laea, 2500, replace = F)
#dat4 <- sample_n(locs_laea, 500, replace = F)

saveRDS(dat3, "harp data/harps2500.rds")


####################################
### Generate land shape and mesh ###
####################################

land <- mapr(dat3, prj, buff = 5e5)
b <- meshr(dat3, prj, buff = 5e5, keep = 0.02) 

#have a look at meshes powerpoint from inlabru course
#boundary = list(area of interest, list(buffer on area of interest, coastline))
mesh = inla.mesh.2d(loc = cbind(dat3$x, dat3$y), boundary = b, max.edge = c(250000, 1e+06), cutoff = 50000, max.n = 500)

ggplot() +
  inlabru::gg(mesh) + 
  geom_sf(aes(), data = land) +
  coord_sf(xlim = c(-4000000, 4000000), ylim = c(-4000000, 4000000), crs = prj, expand = T) +
  geom_sf(aes(), data = st_as_sf(dat3, coords = c("lon", "lat")) %>% st_set_crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))



