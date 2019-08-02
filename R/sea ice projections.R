
require(raster)
require(tidyverse)

# Jan 2006 to Dec 2100
# 1:1140
# Start Dec 2007 - Feb 2008

# Load in sea ice data from Pearse
#hist <- raster::stack("/Users/jamesgrecian/Dropbox/James Grecian - Sea ice and seals/sic_IPSL_1x1regrid_historical.nc")
rcp26 <- raster::stack("/Users/jamesgrecian/Dropbox/James Grecian - Sea ice and seals/sic_IPSL_1x1regrid_rcp26.nc")
rcp85 <- raster::stack("/Users/jamesgrecian/Dropbox/James Grecian - Sea ice and seals/sic_IPSL_1x1regrid_rcp85.nc")

#only intersted in northern hemisphere and area the harp seals use
CP <- extent(-100, 100, 0, 90)
rcp26 <- crop(rcp26, CP)
rcp85 <- crop(rcp85, CP)

# Create sequence of March indices
id <- seq(from = 3, to = 1140, by = 12)
# extract the RCP26 contour and calculate the mean latitude
rcp26_mar <- NA
for(i in id){
  foo <- rasterToContour(subset(rcp26, i), levels = 15)
  foo <- fortify(foo)
  rcp26_mar[i] <- mean(foo$lat)
}

# extract the RCP85 contour and calculate the mean latitude
rcp85_mar <- NA
for(i in id){
  foo <- rasterToContour(subset(rcp85, i), levels = 15)
  foo <- fortify(foo)
  rcp85_mar[i] <- mean(foo$lat)
}

# Create dataframe
df_mar <- as_tibble(cbind(id = id,
                          year = seq(from = 2006, to = 2100, by = 1),
                          rcp26 = rcp26_mar[!is.na(rcp26_mar)],
                          rcp85 = rcp85_mar[!is.na(rcp85_mar)]))
df_mar <- gather(df_mar, key = "rcp", value = "lat", c("rcp26", "rcp85"))

# Generate a plot
p1 <- ggplot() +
  geom_point(aes(x = year, y = lat), data = df_mar) + 
  geom_smooth(aes(x = year, y = lat), method = mgcv::gam, formula = y ~ s(x, bs = "cs"), data = df_mar) +
  xlab("Year") + ylab("Mean latitude of 15% sea ice contour in March") +
  facet_wrap(~rcp)

quartz(width = 10, height = 6) 
print(p1)
quartz.save(file = "~/March sea ice contour rcp26 rcp85.jpeg", type = "jpeg", device = dev.cur(), dpi = 300)
dev.off()



id <- seq(from = 9, to = 1140, by = 12)
rcp26_sep <- NA
for(i in id){
  foo <- tryCatch(rasterToContour(subset(rcp26, i), levels = 15), error = function(err) NA)
  foo <- tryCatch(fortify(foo), error = function(err) NA)
  rcp26_sep[i] <- tryCatch(mean(foo$lat, na.rm = T), error = function(err) 90)
}
rcp85_sep <- NA
for(i in id){
  foo <- tryCatch(rasterToContour(subset(rcp85, i), levels = 15), error = function(err) NA)
  foo <- tryCatch(fortify(foo), error = function(err) NA)
  rcp85_sep[i] <- tryCatch(mean(foo$lat, na.rm = T), error = function(err) 90)
}

df_sep <- as_tibble(cbind(id = id,
                          year = seq(from = 2006, to = 2100, by = 1),
                          rcp26 = rcp26_sep[!is.na(rcp26_sep)],
                          rcp85 = rcp85_sep[!is.na(rcp85_sep)]))
df_sep <- gather(df_sep, key = "rcp", value = "lat", c("rcp26", "rcp85"))

# Generate a plot
p2 <- ggplot() +
  geom_point(aes(x = year, y = lat), data = df_sep) + 
  geom_smooth(aes(x = year, y = lat), method = mgcv::gam, formula = y ~ s(x, bs = "cs"), data = df_sep) +
  xlab("Year") + ylab("Mean latitude of 15% sea ice contour in September") +
  facet_wrap(~rcp)

quartz(width = 10, height = 6) 
print(p2)
quartz.save(file = "~/September sea ice contour rcp26 rcp85.jpeg", type = "jpeg", device = dev.cur(), dpi = 300)
dev.off()









df_sep26 <- as_tibble(cbind(id = id,
                            rcp26 = rcp26_sep[!is.na(rcp26_sep)],
                            year = seq(from = 2006, to = 2100, by = 1)))

rcp85_sep <- NA
#some years there is no ice...
id <- id[-match(c(777, 873, 993, 1029, 1125, 1137), id)]
for(i in id){
  foo <- rasterToContour(subset(rcp85, i), levels = 15)
  foo <- fortify(foo)
  rcp85_sep[i] <- mean(foo$lat, na.rm = T)
}
df_sep85 <- as_tibble(cbind(id = id,
                            rcp85 = rcp85_sep[!is.na(rcp85_sep)]))
#add the blank years back in
df_sep85 <- df_sep85 %>% add_row(id = c(777, 873, 993, 1029, 1125, 1137))
df_sep85 <- df_sep85 %>% arrange(id)
df_sep85 <- df_sep85 %>% mutate(year = seq(from = 2006, to = 2100, by = 1))

df_sep <- cbind.data.frame(df_sep26, df_sep85)
df_sep <- df_sep[,c(1, 3, 2, 5)]
df_sep <- gather(df_sep, key = "rcp", value = "lat", c("rcp26", "rcp85"))

df_sep %>%
  ggplot() +
  geom_point(aes(x = year, y = lat)) + 
  geom_smooth(aes(x = year, y = lat), method = mgcv::gam, formula = y ~ s(x, bs = "cs")) +
  facet_wrap(~rcp)




foo <- rasterToContour(subset(rcp85, 777), levels = 15)
foo <- fortify(foo)
mean(foo$lat, na.rm = T)


foo <- rasterToContour(subset(rcp85, 849), levels = 15)



#df <- as_tibble(cbind(id = id,
#                      lat = out[!is.na(out)],
#                      group = rep(1:95, each = 3)))

#df <- df %>% group_by(group) %>% 
#  summarize(mean = mean(lat, na.rm=TRUE))

                          




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
  
  
  




# Dec-Feb is 12:14
# Mar-May is 15:17
# Jun-Aug is 18:20
# Sep-Nov is 21:23
#id <- c(seq(from = 15, to = 1130, by = 12),
#        seq(from = 16, to = 1130, by = 12),
#        seq(from = 17, to = 1130, by = 12))
#id <- sort(id)



