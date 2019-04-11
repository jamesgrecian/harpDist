###########################################
### Append time index to harp seal data ###
###########################################

require(dplyr)

# Import harp seal data
dat <- readRDS("harp data/harps500.rds")

# Harp seal data spans 1995 to 2018
# To simplify analysis split harp seals into 4 seasons and 5 year bins
# Group into 20 time points
dat <- dat %>% arrange(date)
dat <- dat %>% mutate(year = lubridate::year(date))
dat <- dat %>% mutate(month = lubridate::month(date))

# Initially tried this with mutate and filter in dplyr
# Old skool subsetting easier...
dat$index <- NA

years <- c(1995:1999)
dat$index[dat$month == 12 & dat$year %in% c(years-1)] <- 1
dat$index[dat$month %in% c(1:2) & dat$year %in% c(years)] <- 1
dat$index[dat$month %in% c(3:5) & dat$year %in% c(years)] <- 2
dat$index[dat$month %in% c(6:8) & dat$year %in% c(years)] <- 3
dat$index[dat$month %in% c(9:11) & dat$year %in% c(years)] <- 4

years <- c(2000:2004)
dat$index[dat$month == 12 & dat$year %in% c(years-1)] <- 5
dat$index[dat$month %in% c(1:2) & dat$year %in% c(years)] <- 5
dat$index[dat$month %in% c(3:5) & dat$year %in% c(years)] <- 6
dat$index[dat$month %in% c(6:8) & dat$year %in% c(years)] <- 7
dat$index[dat$month %in% c(9:11) & dat$year %in% c(years)] <- 8

years <- c(2005:2009)
dat$index[dat$month == 12 & dat$year %in% c(years-1)] <- 9
dat$index[dat$month %in% c(1:2) & dat$year %in% c(years)] <- 9
dat$index[dat$month %in% c(3:5) & dat$year %in% c(years)] <- 10
dat$index[dat$month %in% c(6:8) & dat$year %in% c(years)] <- 11
dat$index[dat$month %in% c(9:11) & dat$year %in% c(years)] <- 12

years <- c(2010:2014)
dat$index[dat$month == 12 & dat$year %in% c(years-1)] <- 13
dat$index[dat$month %in% c(1:2) & dat$year %in% c(years)] <- 13
dat$index[dat$month %in% c(3:5) & dat$year %in% c(years)] <- 14
dat$index[dat$month %in% c(6:8) & dat$year %in% c(years)] <- 15
dat$index[dat$month %in% c(9:11) & dat$year %in% c(years)] <- 16

years <- c(2015:2019)
dat$index[dat$month == 12 & dat$year %in% c(years-1)] <- 17
dat$index[dat$month %in% c(1:2) & dat$year %in% c(years)] <- 17
dat$index[dat$month %in% c(3:5) & dat$year %in% c(years)] <- 18
dat$index[dat$month %in% c(6:8) & dat$year %in% c(years)] <- 19
dat$index[dat$month %in% c(9:11) & dat$year %in% c(years)] <- 20

# Output
saveRDS(dat, "harp data/harps500_indexed.rds")

# ends
