######################################################
### Format sea ice projections for INLA prediction ###
######################################################

# Pearse sent through some example projections under RCP25 and RCP85
# Group the months into seasons
# Estimate harp seal distributions for mid and late century

# Should this be averaged as per the 'real' data
# i.e. 2045-2049
# i.e. 2095-2099
# Also include 2015-19 for comparison with NSIDC

# Start simple...
# Sea ice projections run from Jan 2006 [1] to Dec 2100 [1140]
# 95 years x 12 months = 1140 layers
# i.e Jan = seq(1, 1140, 12)

# Extract predicted sea ice concentration
# Average for each of the seasons across 5 years...

# Dec 2014 - Feb 2015 = layers 108:110
# Dec 2015 - Feb 2016 = layers 120:122
# Dec 2016 - Feb 2017 = layers 132:134
# Dec 2017 - Feb 2018 = layers 144:146
# Dec 2018 - Feb 2019 = layers 156:158

# This can be written as a function to avoid copying and pasting a lot of code...
summarise_IPCC <- function(data){

  # Create date sequence that matches RCP data
  dates <- seq(lubridate::ymd("2006-01-01"), lubridate::ymd("2100-12-31"), by = "month")
  # Create sequence of ids to match with RCP rasters
  id_1 <- 1:length(dates)
  # Create recieving vector for stackApply ids
  # Start at 1 as stackApply behaves strangely if the first number it encounters isn't the lowest in the sequence
  id_2 <- rep(1, times = length(dates))
  
  id_2[id_1[lubridate::month(dates) %in% c(1:3) & lubridate::year(dates) %in% c(2015:2019)]-1] <- 2
  id_2[id_1[lubridate::month(dates) %in% c(4:6) & lubridate::year(dates) %in% c(2015:2019)]-1] <- 3
  id_2[id_1[lubridate::month(dates) %in% c(7:9) & lubridate::year(dates) %in% c(2015:2019)]-1] <- 4
  id_2[id_1[lubridate::month(dates) %in% c(10:12) & lubridate::year(dates) %in% c(2015:2019)]-1] <- 5
  
  id_2[id_1[lubridate::month(dates) %in% c(1:3) & lubridate::year(dates) %in% c(2045:2049)]-1] <- 6
  id_2[id_1[lubridate::month(dates) %in% c(4:6) & lubridate::year(dates) %in% c(2045:2049)]-1] <- 7
  id_2[id_1[lubridate::month(dates) %in% c(7:9) & lubridate::year(dates) %in% c(2045:2049)]-1] <- 8
  id_2[id_1[lubridate::month(dates) %in% c(10:12) & lubridate::year(dates) %in% c(2045:2049)]-1] <- 9
  
  id_2[id_1[lubridate::month(dates) %in% c(1:3) & lubridate::year(dates) %in% c(2095:2099)]-1] <- 10
  id_2[id_1[lubridate::month(dates) %in% c(4:6) & lubridate::year(dates) %in% c(2095:2099)]-1] <- 11
  id_2[id_1[lubridate::month(dates) %in% c(7:9) & lubridate::year(dates) %in% c(2095:2099)]-1] <- 12
  id_2[id_1[lubridate::month(dates) %in% c(10:12) & lubridate::year(dates) %in% c(2095:2099)]-1] <- 13
  

  out <- raster::stackApply(data, indices = id_2, fun = mean, na.rm = T)
  out <- raster::subset(out, 2:13)
  names(out) <- c("Dec14_Feb19", "Mar15_May19", "Jun15_Aug19", "Sep15_Nov19",
                  "Dec44_Feb49", "Mar45_May49", "Jun45_Aug49", "Sep45_Nov49",
                  "Dec94_Feb99", "Mat95_May99", "Jun95_Aug99", "Sep95_Nov99")
  return(out)
}

# The function is based on this longer way of doing it...
# #Subset to extract seasonal mean across 5 years
#rcp26_seasonal <- stack(mean(subset(rcp26, ids[month(dates) %in% c(1:3) & year(dates) %in% c(2015:2019)]-1), na.rm = T),
#                        mean(subset(rcp26, ids[month(dates) %in% c(4:6) & year(dates) %in% c(2015:2019)]-1), na.rm = T),
#                        mean(subset(rcp26, ids[month(dates) %in% c(7:9) & year(dates) %in% c(2015:2019)]-1), na.rm = T),
#                        mean(subset(rcp26, ids[month(dates) %in% c(10:12) & year(dates) %in% c(2015:2019)]-1), na.rm = T),
#                        mean(subset(rcp26, ids[month(dates) %in% c(1:3) & year(dates) %in% c(2045:2049)]-1), na.rm = T),
#                        mean(subset(rcp26, ids[month(dates) %in% c(4:6) & year(dates) %in% c(2045:2049)]-1), na.rm = T),
#                        mean(subset(rcp26, ids[month(dates) %in% c(7:9) & year(dates) %in% c(2045:2049)]-1), na.rm = T),
#                        mean(subset(rcp26, ids[month(dates) %in% c(10:12) & year(dates) %in% c(2045:2049)]-1), na.rm = T),
#                        mean(subset(rcp26, ids[month(dates) %in% c(1:3) & year(dates) %in% c(2095:2099)]-1), na.rm = T),
#                        mean(subset(rcp26, ids[month(dates) %in% c(4:6) & year(dates) %in% c(2095:2099)]-1), na.rm = T),
#                        mean(subset(rcp26, ids[month(dates) %in% c(7:9) & year(dates) %in% c(2095:2099)]-1), na.rm = T),
#                        mean(subset(rcp26, ids[month(dates) %in% c(10:12) & year(dates) %in% c(2095:2099)]-1), na.rm = T))

