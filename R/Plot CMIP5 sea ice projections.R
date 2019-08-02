###########################################
### How to use the sea ice projections? ###
###########################################

# Pearse sent through some example projections under RCP25 and RCP85
# Group the months into seasons
# Estimate harp seal distributions for mid and late century

# Should this be averaged over
# i.e. 2050-2054
# i.e. 2094-2099

# or just 2050 and 2100???
# How to use this if using mean centred covariates...?

# Start simple...
# Sea ice projections run from Jan 2006 to Dec 2100
# 95 years x 12 months = 1140 layers
# Extract predicted sea ice concentration in:
# Dec 2049 - Feb 2050 = layers 528:530
# Mar 2050 - May 2050 = layers 531:533
# Jun 2050 - Aug 2050 = layers 534:536
# Sep 2050 - Nov 2050 = layers 537:539

# Dec 2099 - Feb 2100 = layers 1128:1130
# Mar 2100 - May 2100 = layers 1131:1133
# Jun 2100 - Aug 2100 = layers 1134:1136
# Sep 2100 - Nov 2100 = layers 1137:1139

 # Load libraries
require(raster)
require(tidyverse)
require(sf)
source("discrete_gradient.R")

# Load in sea ice data from Pearse
rcp26 <- raster::stack("/Users/jamesgrecian/Dropbox/James Grecian - Sea ice and seals/sic_IPSL_1x1regrid_rcp26.nc")
rcp85 <- raster::stack("/Users/jamesgrecian/Dropbox/James Grecian - Sea ice and seals/sic_IPSL_1x1regrid_rcp85.nc")
#only intersted in northern hemisphere and area the harp seals use
CP <- extent(-180, 180, 0, 90)
rcp26 <- crop(rcp26, CP)
rcp85 <- crop(rcp85, CP)

ice_2050_rcp26 <- stack(mean(subset(rcp26, 528:530), na.rm = T),
                        mean(subset(rcp26, 531:533), na.rm = T),
                        mean(subset(rcp26, 534:536), na.rm = T),
                        mean(subset(rcp26, 537:539), na.rm = T))

ice_2100_rcp26 <- stack(mean(subset(rcp26, 1128:1130), na.rm = T),
                        mean(subset(rcp26, 1131:1133), na.rm = T),
                        mean(subset(rcp26, 1134:1136), na.rm = T),
                        mean(subset(rcp26, 1137:1139), na.rm = T))

ice_2050_rcp85 <- stack(mean(subset(rcp85, 528:530), na.rm = T),
                        mean(subset(rcp85, 531:533), na.rm = T),
                        mean(subset(rcp85, 534:536), na.rm = T),
                        mean(subset(rcp85, 537:539), na.rm = T))

ice_2100_rcp85 <- stack(mean(subset(rcp85, 1128:1130), na.rm = T),
                        mean(subset(rcp85, 1131:1133), na.rm = T),
                        mean(subset(rcp85, 1134:1136), na.rm = T),
                        mean(subset(rcp85, 1137:1139), na.rm = T))


#prj = "+proj=laea +lat_0=75 +lon_0=-25 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
# Use NSIDC projection
prj = "+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +a=6378273 +b=6356889.449 +units=m +no_defs"

ice_2050_rcp26_df = projectRaster(ice_2050_rcp26, crs = CRS(prj), method = 'ngb', res = 50000) #bilinear seems to give NAs?!
ice_2050_rcp26_df <- as.data.frame(rasterToPoints(ice_2050_rcp26_df))
names(ice_2050_rcp26_df)[3:6] <- c("Winter", "Spring", "Summer", "Autumn")
ice_2050_rcp26_df <- ice_2050_rcp26_df %>% gather(key = "season", value = "ice", Winter:Autumn)
ice_2050_rcp26_df$season <- factor(ice_2050_rcp26_df$season,
                                   levels = c("Winter", "Spring", "Summer", "Autumn"))
ice_2100_rcp26_df = projectRaster(ice_2100_rcp26, crs = CRS(prj), method = 'ngb', res = 50000) #bilinear seems to give NAs?!
ice_2100_rcp26_df <- as.data.frame(rasterToPoints(ice_2100_rcp26_df))
names(ice_2100_rcp26_df)[3:6] <- c("Winter", "Spring", "Summer", "Autumn")
ice_2100_rcp26_df <- ice_2100_rcp26_df %>% gather(key = "season", value = "ice", Winter:Autumn)
ice_2100_rcp26_df$season <- factor(ice_2100_rcp26_df$season,
                                   levels = c("Winter", "Spring", "Summer", "Autumn"))
ice_2050_rcp26_df$year = 2050  
ice_2100_rcp26_df$year = 2100
ice_rcp26 <- rbind.data.frame(ice_2050_rcp26_df,
                              ice_2100_rcp26_df)
ice_rcp26$scenario = "rcp26" 

ice_2050_rcp85_df = projectRaster(ice_2050_rcp85, crs = CRS(prj), method = 'ngb', res = 50000) #bilinear seems to give NAs?!
ice_2050_rcp85_df <- as.data.frame(rasterToPoints(ice_2050_rcp85_df))
names(ice_2050_rcp85_df)[3:6] <- c("Winter", "Spring", "Summer", "Autumn")
ice_2050_rcp85_df <- ice_2050_rcp85_df %>% gather(key = "season", value = "ice", Winter:Autumn)
ice_2050_rcp85_df$season <- factor(ice_2050_rcp85_df$season,
                                   levels = c("Winter", "Spring", "Summer", "Autumn"))
ice_2100_rcp85_df = projectRaster(ice_2100_rcp85, crs = CRS(prj), method = 'ngb', res = 50000) #bilinear seems to give NAs?!
ice_2100_rcp85_df <- as.data.frame(rasterToPoints(ice_2100_rcp85_df))
names(ice_2100_rcp85_df)[3:6] <- c("Winter", "Spring", "Summer", "Autumn")

ice_2100_rcp85_df <- ice_2100_rcp85_df %>% gather(key = "season", value = "ice", Winter:Autumn)



ice_2100_rcp85_df$season <- factor(ice_2100_rcp85_df$season,
                                   levels = c("Winter", "Spring", "Summer", "Autumn"))
ice_2050_rcp85_df$year = 2050  
ice_2100_rcp85_df$year = 2100
ice_rcp85 <- rbind.data.frame(ice_2050_rcp85_df,
                              ice_2100_rcp85_df)
ice_rcp85$scenario = "rcp85" 

future_ice <- rbind.data.frame(ice_rcp26,
                               ice_rcp85)


# Generate land map
dat <- readRDS("harp data/harps2500_indexed.rds")
land <- mapr::mapr(dat,
                   prj,
                   buff = 5000000)
#ext <- raster::extent(-95, 65, 40, 85) %>% st_bbox() %>% st_as_sfc()
#ext <- ext %>% st_segmentize(1)
#ext <- ext %>% st_set_crs(4326)
#ext <-ext %>% sf::st_transform(crs = prj)

# get the bounding box in transformed coordinates and expand by 10%
#xlim <- st_bbox(ext)[c("xmin", "xmax")]*1.5
#ylim <- st_bbox(ext)[c("ymin", "ymax")]*1.5

# turn into enclosing rectangle
#encl_rect <- 
#  list(
#    cbind(
#      c(xlim[1], xlim[2], xlim[2], xlim[1], xlim[1]), 
#      c(ylim[1], ylim[1], ylim[2], ylim[2], ylim[1])
#    )
#  ) %>%
#  st_polygon() %>%
#  st_sfc(crs = prj)

# calculate the area outside the earth outline as the difference
# between the enclosing rectangle and the earth outline
#cookie <- st_difference(encl_rect, ext)

# Add this to the plot to mask out the areas not of interest
# More effort required to move axis labels...

unique(future_ice$scenario)
# Plot to check

ggplot() +
  theme_minimal() + ylab("") + xlab("") +
  geom_raster(aes(x = x, y = y, fill = ice), data = future_ice) +
  scale_fill_distiller("% Sea Ice", palette = "Blues", limits = c(0,100)) +
  geom_sf(aes(), fill = "grey", colour = "dark grey", data = land) +
  facet_grid(season ~ year)
#  geom_sf(fill = "white", color = "black", data = cookie) +
#  coord_sf(xlim = c(-5500000, 5500000), ylim = c(-3900000, 1600000), crs = prj, expand = F) +
#  facet_wrap(~ scenario + year + season,
#             ncol = 4,
#             strip.position = "left")


p <- ggplot() +
  ylab("") + xlab("") +
  geom_point(aes(x = x, y = y, colour = ice), data = future_ice) +
  geom_raster(aes(x = x, y = y, fill = ice), data = future_ice) +
  scale_colour_discrete_gradient("Sea Ice Concentration (%)",
                                 colours = rev(RColorBrewer::brewer.pal(10, "Blues")),
                                 limits = c(0, 100),
                                 breaks = seq(0, 100, 10),
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
  coord_sf(xlim = c(-3650000, 3750000), ylim = c(-5250000, 3000000), crs = "+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +a=6378273 +b=6356889.449 +units=m +no_defs", expand = F) +
  scale_x_continuous(breaks = seq(-180, 180, 20)) + scale_y_continuous(breaks = seq(40, 90, 10))

quartz(title = "Panel Plot", width = 10, height = 9)
p + facet_grid(season ~ scenario + year)
quartz.save(file = "~/Desktop/CMIP5 multipanel plot.jpeg", type = "jpeg",
            dev  = dev.cur(), dpi = 500)
dev.off()
