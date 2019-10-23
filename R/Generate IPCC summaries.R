#############################################
### Generate seasonal sea ice projections ###
#############################################

require(raster)
source("R/summarise_IPCC.R")

# Load in sea ice data from Pearse
rcp26 <- raster::stack("/Users/jamesgrecian/Dropbox/James Grecian - Sea ice and seals/sic_IPSL_1x1regrid_rcp26.nc")
rcp85 <- raster::stack("/Users/jamesgrecian/Dropbox/James Grecian - Sea ice and seals/sic_IPSL_1x1regrid_rcp85.nc")

# only intersted in northern hemisphere and area the harp seals use
CP <- raster::extent(-180, 180, 0, 90)
rcp26 <- raster::crop(rcp26, CP)
rcp85 <- raster::crop(rcp85, CP)

# Try function..!
rcp26_seasonal <- summarise_IPCC(rcp26)
rcp85_seasonal <- summarise_IPCC(rcp85)

# Write output
raster::writeRaster(rcp26_seasonal, filename = "data/IPSL_rcp26_seasonal", bandorder='BIL', overwrite = T)
raster::writeRaster(rcp85_seasonal, filename = "data/IPSL_rcp85_seasonal", bandorder='BIL', overwrite = T)

