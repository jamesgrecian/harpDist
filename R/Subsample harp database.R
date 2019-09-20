###############################################
### Subsample harp seal point data for INLA ###
###############################################

# Load libraries
require(tidyverse)
require(sf)
#devtools::install_github("jamesgrecian/mapr")
require(mapr)
require(RSRDL)
require(viridis)
require(lubridate)

# Define albers projection
prj = "+proj=laea +lat_0=75 +lon_0=-25 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"

##########################
### Load in point data ###
##########################

# Load in Martin's harp database
dat = get.SRDLpg(theDB = 'harps', theTable = 'diag', theFields = 'All', theDep = 'All', theRef = 'All',
                 theHost = 'localhost', thePort = 5432, theUser = 'postgres', thePwd = 'admin') %>% as_tibble()
dat <- dat %>% dplyr::select("ref", "d.date", "lq", "lon", "lat", "semi.major.axis", "semi.minor.axis", "ellipse.orientation")
dat = arrange(dat, ref, d.date) # resort by id and date
dat = dat[!duplicated(dat$d.date),] # remove duplicates

#Exclude animals that have fewer than 20 locations
#Remove animals where tag failed within 1 week
dat = dat[!dat$ref %in% c("hp1-9276-04",
                          "hp1-9278-04",
                          "hp1c-296-96",
                          "hp2-9325-04",
                          "hp2b-08-99",
                          "hp4-L517-17",
                          "hp1b-01-95",
                          "hp1b-02-95",
                          "hp2-9334-04"),]
names(dat) <- c("id", "date", "lc", "lon", "lat", "smaj", "smin", "eor")

# Load hp5
hp5 = get.SRDLpg(theDB = 'hp5', theTable = 'diag', theFields = 'All', theDep = 'All', theRef = 'All',
                 theHost = 'localhost', thePort = 5432, theUser = 'postgres', thePwd = 'admin') %>% as_tibble()
hp5 <- hp5 %>% dplyr::select("REF", "D_DATE", "LQ", "LON", "LAT", "SEMI_MAJOR_AXIS", "SEMI_MINOR_AXIS", "ELLIPSE_ORIENTATION")
hp5 = hp5[hp5$REF %in% c("hp5-L764-18", "hp5-L766-18"),] 
hp5 <- hp5 %>% rename(id = REF,
                      date = D_DATE,
                      lc = LQ,
                      lon = LON,
                      lat = LAT,
                      smaj = SEMI_MAJOR_AXIS,
                      smin = SEMI_MINOR_AXIS,
                      eor = ELLIPSE_ORIENTATION)
# Load in hp6
hp6 <- Hmisc::mdb.get("~/hp6.mdb",  tables = "diag") %>% as_tibble()
hp6 <- hp6 %>% dplyr::select("REF", "D.DATE", "LQ", "LON", "LAT", "SEMI.MAJOR.AXIS", "SEMI.MINOR.AXIS", "ELLIPSE.ORIENTATION")
hp6 <- hp6 %>% rename(id = REF,
                      date = D.DATE,
                      lc = LQ,
                      lon = LON,
                      lat = LAT,
                      smaj = SEMI.MAJOR.AXIS,
                      smin = SEMI.MINOR.AXIS,
                      eor = ELLIPSE.ORIENTATION)
# Format date time
hp6 <- hp6 %>% mutate(date = mdy_hms(date, tz = "UTC"))
# Filter out data prior to deployment
hp6 <- hp6 %>% filter(date > "2019-03-22 00:00:01")

# Combine three datasets
dat <- rbind(dat, hp5, hp6)
dat <- dat %>% arrange(id, date)

# Recode location class for prefilter algoritm
dat <- dat %>% mutate(lc = recode_factor(lc,
                                         `3` = "3",
                                         `2` = "2",
                                         `1` = "1",
                                         `0` = "0",
                                         `-1` = "A",
                                         `-2` = "B",
                                         `-9` = "Z"))
dat <- dat[!dat$lc == "Z",] # drop the Z locations
dat <- dat[dat$lat > 40,] # drop locations below 30 N
dat <- dat[dat$lat < 85,] # drop locations above 85 N
dat <- dat[dat$lon > -100,] # drop locations below 100 W
dat <- dat[dat$lon < 100,] # drop locations above 100 E

####################################
### Error correct using foieGras ###
####################################

# Run the model through foieGras to error correct and regularise
require(foieGras)
fls <- foieGras::fit_ssm(dat, model = "rw", time.step = 24)
pred <- foieGras::grab(fls, what = "predicted", as_sf = F) #this should probably be predicted... but will need tweaking for particular individuals...

# Subsample for output...
out <- sample_n(pred, 2500, replace = F)
saveRDS(out, "data/harps2500.rds")

# Generate land shapefile for plot
land <- mapr(pred, prj, buff = 5e5)

p1 <- ggplot() +
  geom_sf(aes(), data = land) +
  geom_sf(aes(), data = out %>% st_as_sf(coords = c("lon", "lat")) %>% st_set_crs(4326)) +
  coord_sf(xlim = c(-3000000, 3000000), ylim = c(-3500000, 2500000), crs = prj, expand = T)
print(p1)

# ends
