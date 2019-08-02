# Load libraries
require(raster)
require(tidyverse)
require(sf)
source("~/Dropbox/git_projects/harpDist/discrete_gradient.R")

# Load in sea ice data from Pearse
rcp26 <- raster::stack("/Users/jamesgrecian/Dropbox/James Grecian - Sea ice and seals/sic_IPSL_1x1regrid_rcp26.nc")
rcp85 <- raster::stack("/Users/jamesgrecian/Dropbox/James Grecian - Sea ice and seals/sic_IPSL_1x1regrid_rcp85.nc")
#only intersted in northern hemisphere and area the harp seals use
CP <- extent(-180, 180, 0, 90)
rcp26 <- crop(rcp26, CP)
rcp85 <- crop(rcp85, CP)

zones <- init(rcp26, v = 'y')
#foo <- zonal(rcp26, zones, 'mean', na.rm = T)

xvals <- matrix(,46, nlayers(rcp26))
for (i in 1:nlayers(rcp26)){
  xvals[,i] = zonal(subset(rcp85, i), zones, 'mean', na.rm = T)[,2]
}
xvals <- xvals %>% as_tibble() 

# Sea ice projections run from Jan 2006 to Dec 2100
require(lubridate)
dates <- seq(ymd("2006-01-01"), ymd("2100-12-1"), by = "1 month")
names(xvals) <- dates
xvals <- xvals %>% gather(., key = year, value = ice, 1:1140)
xvals <- xvals %>% mutate(year =  ymd(year))
xvals <- xvals %>% mutate(latitude = rep(zonal(subset(rcp26, i), zones, 'mean', na.rm = T)[,1], 1140))

ggplot() +
  theme_bw() +
  geom_tile(aes(x = year, y = latitude, fill = ice, colour = ice), data = xvals, alpha = 0.5) +
  geom_contour(aes(x = year, y = latitude, z = ice), breaks = 15, data = xvals, colour = "black") +
  scale_x_date(date_breaks = "10 years", date_labels = "%Y", limits = c(ymd("2010-01-01"), ymd("2100-12-1")), expand = c(0, 0)) +
  scale_y_continuous(limits = c(39, 91), expand = c(0, 0)) +
  #scale_fill_distiller("% Sea Ice", palette = "Blues", limits = c(0,100)) +
  ylab("Latitude") + xlab("Season:Year") +
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
                                                                                    hjust = 0.5)))



