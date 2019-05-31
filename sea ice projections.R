
require(raster)
require(tidyverse)


# Jan 2006 to Dec 2100
# 1:1140
# Start Dec 2007 - Feb 2008
#

id <- c(seq(from = 12, to = 1130, by = 12),
        seq(from = 13, to = 1130, by = 12),
        seq(from = 14, to = 1130, by = 12))
id <- sort(id)

#hist <- raster::stack("/Users/jamesgrecian/Dropbox/James Grecian - Sea ice and seals/sic_IPSL_1x1regrid_historical.nc")
rcp26 <- raster::stack("/Users/jamesgrecian/Dropbox/James Grecian - Sea ice and seals/sic_IPSL_1x1regrid_rcp26.nc")
rcp85 <- raster::stack("/Users/jamesgrecian/Dropbox/James Grecian - Sea ice and seals/sic_IPSL_1x1regrid_rcp85.nc")

#only intersted in northern hemisphere
CP <- extent(-180, 180, 0, 90)
rcp26 <- crop(rcp26, CP)
rcp85 <- crop(rcp85, CP)

# extract the contour and calculate the mean latitude
out <- NA
for(i in id){
  foo <- rasterToContour(subset(rcp85, i), levels = 15)
  foo <- fortify(foo)
  out[i] <- mean(foo$lat)
}

plot(out, type = "l")

df <- as_tibble(cbind(id = id,
                      lat = out[!is.na(out)],
                      group = rep(1:94, each = 3)))

df <- df %>% group_by(group) %>% 
  summarize(mean = mean(lat, na.rm=TRUE))

df %>% ggplot() + geom_line(aes(x = group, y = mean))


prj = "+proj=laea +lat_0=75 +lon_0=-25 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

rcp26_prj = projectRaster(subset(rcp26, 1140), crs = CRS(prj), method = 'ngb', res = 50000) #bilinear seems to give NAs?!
rcp26_prj <- as.data.frame(rasterToPoints(rcp26_prj))
names(rcp26_prj)[3] <- "z"

# Plot to check
ggplot() +
  geom_raster(aes(x = x, y = y, fill = z), data = rcp26_prj) +
  scale_fill_distiller("% Sea Ice", palette = "Blues", limits = c(0,100)) +
  theme_bw() + ylab("") + xlab("") +
  coord_sf(xlim = c(-4000000, 3000000), ylim = c(-4000000, 3000000), crs = prj, expand = F) +
  
  
  
  #  geom_sf(aes(), fill = "grey", colour = "grey", data = land) +
  geom_sf(aes(), data = st_as_sf(dat, coords = c("lon", "lat")) %>%
            st_set_crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) +
  ggtitle("Neumann model")




