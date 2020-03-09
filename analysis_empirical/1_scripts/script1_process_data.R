setwd("~/Projects/MigrationTrends/scripts")

# Required packages
my_packs <- c(
  
  # For data management and plotting
  'tidyverse','reshape2','viridis',
  
  # For analysis
  'jagsUI',
  
  # For parallel processing
  'parallel','doParallel'
  
)

# if any of them are not installed, install them
if (any(!my_packs %in% installed.packages()[, 'Package'])) {install.packages(my_packs[which(!my_packs %in% installed.packages()[, 'Package'])],dependencies = TRUE)}
lapply(my_packs, require, character.only = TRUE)

rm(list=ls())

#******************************************************************************************************************************************
#******************************************************************************************************************************************
# PART 1: FORMAT MIGRATION COUNT DATA
#******************************************************************************************************************************************
#******************************************************************************************************************************************

#------------------------------------------------------------------------------------------------------------------------------------------
# Read/format data from BSC (Canadian data)
#------------------------------------------------------------------------------------------------------------------------------------------

#Data provided by Danielle Ethier (BSC).  Has been pre-processed/cleaned.  Does not include offsets or effort (i.e., hours nets were open).

dat_can <- rbind(read.csv("../data/migration_counts/CAN/2020_01_03/ACBO.BLPW.2018.csv") %>% add_column(., station = "ACBO"),
                 read.csv("../data/migration_counts/CAN/2020_01_03/BPBO.BLPW.2018.csv") %>% add_column(., station = "BPBO"),
                 read.csv("../data/migration_counts/CAN/2020_01_03/IPBO.BLPW.2018.csv") %>% add_column(., station = "IPBO"),
                 
                 read.csv("../data/migration_counts/CAN/2020_01_03/LMBO.BLPW.2018.csv") %>% add_column(., station = "LMBO"),
                 read.csv("../data/migration_counts/CAN/2020_01_03/LPBO.BLPW.2018.csv") %>% add_column(., station = "LPBO"),
                 read.csv("../data/migration_counts/CAN/2020_01_03/MGBO.BLPW.2018.csv") %>% add_column(., station = "MGBO"),
                 
                 read.csv("../data/migration_counts/CAN/2020_01_03/MNO.BLPW.2018.csv") %>% add_column(., station = "MNO"),
                 read.csv("../data/migration_counts/CAN/2020_01_03/PEPBO.BLPW.2018.csv") %>% add_column(., station = "PEPBO"),
                 read.csv("../data/migration_counts/CAN/2020_01_03/PIBO.BLPW.2018.csv") %>% add_column(., station = "PIBO"),
                 
                 read.csv("../data/migration_counts/CAN/2020_01_03/RUTH.BLPW.2018.csv") %>% add_column(., station = "RUTH"),
                 read.csv("../data/migration_counts/CAN/2020_01_03/TCBO.BLPW.2018.csv") %>% add_column(., station = "TCBO"),
                 read.csv("../data/migration_counts/CAN/2020_01_03/TLBBS.BLPW.2018.csv") %>% add_column(., station = "TLBBS"),
                 read.csv("../data/migration_counts/CAN/2020_01_03/TTPBRS.BLPW.2018.csv") %>% add_column(., station = "TTPBRS")
)

# Create a column to distinguish sub-areas at LPBO (area will be 1 for all other stations)
dat_can$area <- 1
dat_can$area[which(dat_can$SurveyAreaIdentifier == "LPBO2")] <- 2
dat_can$area[which(dat_can$SurveyAreaIdentifier == "LPBO3")] <- 3

# Number of days with data each year at each station
day_range_per_station <- aggregate(doy~YearCollected + SurveyAreaIdentifier + season, data = dat_can, FUN = range)
ndays_per_station <- aggregate(doy~YearCollected + SurveyAreaIdentifier + season, data = dat_can, FUN = length)
years_per_station <- aggregate(YearCollected ~ SurveyAreaIdentifier + season, data = dat_can, FUN = range)

# Ensure that all stations/days/years/seasons are included as data
area_season_combinations_can <- unique(dat_can[,c("SurveyAreaIdentifier","season")])

dat_combined_can = data.frame()
for (i in 1:nrow(area_season_combinations_can)){
  
  dat <- subset(dat_can, SurveyAreaIdentifier == area_season_combinations_can$SurveyAreaIdentifier[i] & season == area_season_combinations_can$season[i])
  if (nrow(dat) == 0) next
  
  
  min_doy <- min(dat$doy)
  max_doy <- max(dat$doy)
  min_year <- min(dat$YearCollected)
  max_year <- max(dat$YearCollected)
  
  # Create a "full" dataframe to store counts on all days (including NAs)
  dat_full <- expand.grid(YearCollected = seq(min_year,max_year),
                          doy = seq(min_doy,max_doy))
  
  # Fill with counts
  dat_full <- merge(dat_full, dat, all.x = TRUE)
  
  # Ensure relevant data is filled in
  dat_full$SurveyAreaIdentifier = dat$SurveyAreaIdentifier[1]
  dat_full$station = dat$station[1]
  dat_full$area = dat$area[1]
  dat_full$season = dat$season[1]
  dat_full$min_doy = min_doy
  dat_full$max_doy = max_doy
  dat_full$min_year = min_year
  dat_full$max_year = max_year
  
  dat_combined_can = rbind(dat_combined_can, dat_full)
}
dat_combined_can$country = "CAN"

#------------------------------------------------------------------------------------------------------------------------------------------
# Read/format data from Ricky Dunn (USA data)
#------------------------------------------------------------------------------------------------------------------------------------------
dat_usa <- rbind(readxl::read_xlsx("../data/migration_counts/USA/Cleaned BLPW - AIMS spring.xlsx") %>% as.data.frame() %>% add_column(., season = "Spring"),
                 
                 readxl::read_xlsx("../data/migration_counts/USA/Cleaned BLPW - BIBS fall.xlsx") %>% as.data.frame() %>% add_column(., season = "Fall"),
                 
                 readxl::read_xlsx("../data/migration_counts/USA/Cleaned BLPW - BSBO fall.xlsx") %>% as.data.frame() %>% add_column(., season = "Fall"),
                 readxl::read_xlsx("../data/migration_counts/USA/Cleaned BLPW - BSBO spring.xlsx") %>% as.data.frame() %>% add_column(., season = "Spring"),
                 
                 readxl::read_xlsx("../data/migration_counts/USA/Cleaned BLPW - FBBS spring.xlsx") %>% as.data.frame() %>% add_column(., season = "Spring"),
                 readxl::read_xlsx("../data/migration_counts/USA/Cleaned BLPW - FBBS fall.xlsx") %>% as.data.frame() %>% add_column(., season = "Fall"),
                 
                 readxl::read_xlsx("../data/migration_counts/USA/Cleaned BLPW - KWRS fall.xlsx") %>% as.data.frame() %>% add_column(., season = "Fall"),
                 
                 readxl::read_xlsx("../data/migration_counts/USA/Cleaned BLPW - MCCS fall.xlsx") %>% as.data.frame() %>% add_column(., season = "Fall"),
                 readxl::read_xlsx("../data/migration_counts/USA/Cleaned BLPW - MCCS spring.xlsx") %>% as.data.frame() %>% add_column(., season = "Spring"))

# datasets with different column names
tmp = readxl::read_xlsx("../data/migration_counts/USA/Cleaned BLPW - PARC fall.xlsx") %>% as.data.frame() %>% add_column(., season = "Fall")
colnames(tmp) <- colnames(dat_usa)
dat_usa <- rbind(dat_usa, tmp)
rm(tmp)    

tmp = readxl::read_xlsx("../data/migration_counts/USA/Cleaned BLPW - CFMS fall.xlsx") %>% as.data.frame() %>% add_column(., season = "Fall")
colnames(tmp) <- colnames(dat_usa)
dat_usa <- rbind(dat_usa, tmp)
rm(tmp) 

# Fill net N/100 net-hr column
dat_usa$`N/100 net-hr` = dat_usa$N/dat_usa$`Net-hrs`*100
dat_usa = na.omit(dat_usa)
summary(dat_usa)

# Change column names
colnames(dat_usa) = c("station","YearCollected","doy","net.hrs","ObservationCount","N.per.100.net.hrs","season")
dat_usa$area = 1

# Load analysis windows for each station and restrict data to those dates
us_windows = read.csv("../data/migration_counts/USA/US_station_windows.csv")
dat_usa2 = data.frame()
for (i in 1:nrow(us_windows)){
  x = us_windows[i,]
  station = us_windows$station[i]
  area = us_windows$area[i]
  season = us_windows$season[i]
  
  # data inside that range
  dat_usa2 = rbind(dat_usa2, subset(dat_usa, station == x$station & area == x$area & season == x$season & doy >= x$start_date & doy <= x$end_date))
}
dat_usa = dat_usa2
dat_usa = merge(dat_usa, us_windows, all.x = TRUE)


area_season_combinations_usa <- unique(dat_usa[,c("station","season")])
dat_combined_usa = data.frame()
for (i in 1:nrow(area_season_combinations_usa)){
  dat <- subset(dat_usa, station == area_season_combinations_usa$station[i] & season == area_season_combinations_usa$season[i])
  if (nrow(dat) == 0) next
  min_doy <- min(dat$doy)
  max_doy <- max(dat$doy)
  min_year <- min(dat$YearCollected)
  max_year <- max(dat$YearCollected)
  
  # Create a "full" dataframe to store counts on all days (including NAs)
  dat_full <- expand.grid(YearCollected = seq(min_year,max_year),
                          doy = seq(min_doy,max_doy))
  
  # Fill with counts
  dat_full <- merge(dat_full, dat, all.x = TRUE)
  
  # Ensure relevant data is filled in
  dat_full$station = dat$station[1]
  dat_full$station = dat$station[1]
  dat_full$area = dat$area[1]
  dat_full$season = dat$season[1]
  dat_full$min_doy = min_doy
  dat_full$max_doy = max_doy
  dat_full$min_year = min_year
  dat_full$max_year = max_year
  
  dat_combined_usa = rbind(dat_combined_usa, dat_full)
}
dat_combined_usa$country = "USA"

#------------------------------------------------------------------------------------------------------------------------------------------
# Combine CAN and USA data into a single dataframe; generate plots
#------------------------------------------------------------------------------------------------------------------------------------------
dat_combined = dplyr::bind_rows(dat_combined_can, dat_combined_usa)

#----------
# For any sites without recorded net hours, fill in with median
dat_combined$net.hrs[which(is.na(dat_combined$net.hrs))] = median(dat_combined$net.hrs, na.rm = TRUE)
#----------

rm(list=setdiff(ls(), c("dat_combined")))

# Limit to data collected after 2005
dat_combined = subset(dat_combined, YearCollected >= 2006)

#----------------------------------
#----------------------------------
# Plots of daily counts
#----------------------------------
#----------------------------------

# CAN Spring
CAN_Spring_plot <- ggplot(data = subset(dat_combined, country == "CAN" & season == "Spring")) +
  geom_point(aes(x = doy, y = ObservationCount), col = "blue")+
  facet_grid(station~YearCollected, scales = "free_y")+
  xlab("Day of Year")+
  ylab("Count")+
  theme_bw()
pdf(file = paste0(file = "../figures/data_summaries/CAN_Spring_plot.pdf"), width = 30,height=6)
print(CAN_Spring_plot )
dev.off()

# CAN Fall
CAN_Fall_plot <- ggplot(data = subset(dat_combined, country == "CAN" & season == "Fall")) +
  geom_point(aes(x = doy, y = ObservationCount), col = "blue")+
  facet_grid(station~YearCollected, scales = "free_y")+
  xlab("Day of Year")+
  ylab("Count")+
  theme_bw()
pdf(file = paste0(file = "../figures/data_summaries/CAN_Fall_plot.pdf"), width = 30,height=6)
print(CAN_Fall_plot )
dev.off()

# USA Spring
USA_Spring_plot <- ggplot(data = subset(dat_combined, country == "USA" & season == "Spring")) +
  geom_point(aes(x = doy, y = N.per.100.net.hrs), col = "blue")+
  facet_grid(station~YearCollected, scales = "free_y")+
  xlab("Day of Year")+
  ylab("Count(N/100 net hours)")+
  theme_bw()
pdf(file = paste0(file = "../figures/data_summaries/USA_Spring_plot.pdf"), width = 30,height=6)
print(USA_Spring_plot )
dev.off()

# USA Fall
USA_Fall_plot <- ggplot(data = subset(dat_combined, country == "USA" & season == "Fall")) +
  geom_point(aes(x = doy, y = N.per.100.net.hrs), col = "blue")+
  facet_grid(station~YearCollected, scales = "free_y")+
  xlab("Day of Year")+
  ylab("Count(N/100 net hours)")+
  theme_bw()
pdf(file = paste0(file = "../figures/data_summaries/USA_Fall_plot.pdf"), width = 30,height=6)
print(USA_Fall_plot )
dev.off()

#--------------------------------------------
# Restrict analysis to sites with monitoring windows longer than 3 weeks
#--------------------------------------------

# Examine analysis windows (date ranges with observed counts)
date_ranges_start <- aggregate(doy ~ YearCollected + season + station, data = na.omit(dat_combined[,c("YearCollected","season","station","doy","ObservationCount")]), FUN = min)
date_ranges_end <- aggregate(doy ~ YearCollected + season + station, data = na.omit(dat_combined[,c("YearCollected","season","station","doy","ObservationCount")]), FUN = max)
colnames(date_ranges_start)[4] <- "start"
colnames(date_ranges_end)[4] <- "end"
date_ranges <- merge(date_ranges_start, date_ranges_end, all = TRUE)
date_ranges$window_length <- date_ranges$end - date_ranges$start

#Median window length
median_window_lengths <- aggregate(window_length ~ station + season, data = date_ranges, FUN = median)
median_window_lengths$station_season <- paste0(median_window_lengths$station,"_",median_window_lengths$season)
longer_than_3wks <- subset(median_window_lengths, window_length >= 21)

dat_combined$station_season = paste0(dat_combined$station,"_",dat_combined$season)
dat_combined <- subset(dat_combined, station_season %in% longer_than_3wks$station_season)
#--------------------------------------------

write.csv(dat_combined, file = "./processed_data/dat_combined.csv", row.names = FALSE)
