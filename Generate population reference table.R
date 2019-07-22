#################################################################################
### Generate reference table assigning each individual to breeding population ###
#################################################################################

# The three harp seal populations have different distributions
# Generate index referenced against tag id that can be used in INLA dataframe

# Load libraries
require(tidyverse)
require(sf)
#devtools::install_github("jamesgrecian/mapr")
require(mapr)
require(RSRDL)
require(viridis)

# Define albers projection
prj = "+proj=laea +lat_0=75 +lon_0=-25 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"

# Load location database to check movements of animals
dat = get.SRDLpg(theDB = 'harps', theTable = 'diag', theFields = 'All', theDep = 'All', theRef = 'All',
                 theHost = 'localhost', thePort = 5432, theUser = 'postgres', thePwd = 'admin')

# Check format of dataframe
str(dat); summary(dat)
dat$ref = as.factor(dat$ref)
dat$d.date = as.POSIXct(dat$d.date, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
dat$lq = as.factor(dat$lq)
levels(dat$lq) = c("Z", "B", "A", "0", "1", "2", "3") #relabel from -3 to 3

# Resort by id and date
dat = arrange(dat, ref, d.date)

# Remove duplicates
dat = dat[!duplicated(dat$d.date),]

# Exclude animals that have fewer than 20 locations
dat = dat[dat$ref != "hp1-9276-04",]; dat$ref = factor(dat$ref)
dat = dat[dat$ref != "hp1-9278-04",]; dat$ref = factor(dat$ref)
dat = dat[dat$ref != "hp1c-296-96",]; dat$ref = factor(dat$ref)
dat = dat[dat$ref != "hp2-9325-04",]; dat$ref = factor(dat$ref)
dat = dat[dat$ref != "hp2b-08-99",]; dat$ref = factor(dat$ref)
dat = dat[dat$ref != "hp4-L517-17",]; dat$ref = factor(dat$ref)
dat = dat[dat$ref != "hp1b-01-95",]; dat$ref = factor(dat$ref)
dat = dat[dat$ref != "hp1b-02-95",]; dat$ref = factor(dat$ref)

# Remove animals where tag failed within 1 week
dat = dat[dat$ref != "hp2-9334-04",]; dat$ref = factor(dat$ref)

# Format the data for ssmTMB - same format as bsam
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

# Generate land shapefile
land <- mapr(dat2, prj, buff = 5e5)


# Load in deployment info
dat_info = get.SRDLpg(theDB = 'harps', theTable = 'deployments', theFields = 'All', theDep = 'All', theRef = 'All',
                      theHost = 'localhost', thePort = 5432, theUser = 'postgres', thePwd = 'admin')

# Strip out columns of interest
populations <- dat_info %>% select(ref, gref, location, home.lon, home.lat)

# Group into three populations 
populations$location[populations$location == "GSL"] <- "Newfoundland"
populations$location[populations$gref == "hp4"] <- "West Ice"
populations$location[populations$location == "Tromso"] <- "West Ice"

# Plot to check
ggplot() +
  geom_sf(aes(), data = land %>% st_transform(crs = "+proj=longlat +datum=WGS84")) +
  geom_point(aes(x = home.lon, y = home.lat, colour = location), data = populations) +
  #  geom_point(aes(x = lon, y = lat), data = subset(dat, ref == "hp3-Notag-16")) +
  #  geom_point(aes(x = lon, y = lat), data = subset(dat, ref == "hp3-Skinny-14")) +
  geom_point(aes(x = lon, y = lat), data = subset(dat2, id == "hp3-Yellow-16"))

# Ouptut table for reference
saveRDS(populations, "~/populations.rds")

