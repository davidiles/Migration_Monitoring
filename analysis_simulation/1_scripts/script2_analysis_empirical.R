#setwd("~/Projects/MigrationTrends/scripts")

# Required packages
my_packs <- c(
  
  # For data management and plotting
  'tidyverse','reshape2','viridis',
  
  # For analysis
  'jagsUI',
  
  # For parallel processing
  'parallel','doParallel',
  
  'sf','raster','sp','rgdal',
  'maptools','rgeos'
  
)

# if any of them are not installed, install them
if (any(!my_packs %in% installed.packages()[, 'Package'])) {install.packages(my_packs[which(!my_packs %in% installed.packages()[, 'Package'])],dependencies = TRUE)}
lapply(my_packs, require, character.only = TRUE)

rm(list=ls())

#-------------------------------------------------------------------------------------------------------
# Part 1: Define catchment zones
#-------------------------------------------------------------------------------------------------------

# Generate a plot of BCRs and associated range regions
lcc = "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs"

# Read BCR boundaries
bcr1 <- st_read("../data/boundary_shapefiles/bcr/BCR_Terrestrial_master.shp") %>% 
subset(. , COUNTRY %in% c("CANADA","USA") & PROVINCE_S != "HAWAIIAN ISLANDS") %>% 
st_transform( ., crs = lcc)

# Discrete range region zones
bcr1$region <- "South"

bcr1$region[bcr1$PROVINCE_S %in% c("NUNAVUT","SASKATCHEWAN","MANITOBA","ONTARIO","NORTHWEST TERRITORIES")] <- "Central"
bcr1$region[bcr1$PROVINCE_S %in% c("ALASKA","YUKON","BRITISH COLUMBIA","ALBERTA")] <- "West"
bcr1$region[bcr1$PROVINCE_S == "NORTHWEST TERRITORIES" & (bcr1$BCR == 6 | bcr1$BCR == 4)] <- "West"
bcr1$region[bcr1$region == "South" & bcr1$COUNTRY == "CANADA"] <- "East"
bcr1$region <- factor(bcr1$region, levels = c("West","Central","East","South","Far North"))

bcr1$region[bcr1$BCR == 3] <- "Far North"
bcr1$region[bcr1$BCR %in% c(10,5,11,12,9,13)] <- "South"

col_pal <- RColorBrewer::brewer.pal(length(unique(bcr1$region)),"Set2")
col_pal[4:5] <- "gray90"

range_zones_map <- ggplot() +   theme_bw() +
  geom_sf(data = bcr1, aes(fill = region), col = "gray95")+
  scale_fill_manual(values=col_pal)

tiff(filename = "../figures/range_zones_map.tiff", width = 6, height = 6, unit = "in", res = 300)
print(range_zones_map)
dev.off()

#-------------------------------------------------------------------------------------------------------
# Part 2: Plot available stations on map
#-------------------------------------------------------------------------------------------------------

# Read in data and subset to Fall (eventually combine)
dat_combined <- read.csv("./processed_data/dat_combined.csv")
dat_combined <- subset(dat_combined, season == "Fall")

stations_to_keep <- c(

  #West
  "CFMS" = "West",  # Creamer's Field Migration Station
  "TLBBS"= "West", # Teslin Lake
  "MNO"= "West",   # Mackenzie Nature Observatory
  "LMBO"= "West",  # Last Mountain Bird Observatory
  #"TCBO"= "West",


  #Central
  "LPBO" = "Central",
  "PEPBO" = "Central", # Prince edward point
  "BPBO" = "Central",  # Bruce peninsula
  #"PIBO" = "Central",  # Pelee island bird observatory
  #"BSBO" = "Central",
  #"TTPBRS" = "Central",
  #"RUTH" = "Central",



  #East
  "MGBO" = "East",  # McGill
  #"FBBO" = "East",
  #"PARC" = "East",
  "MCCS" = "East",  # Manomet
  #"BIBS" = "East",  # Block Island
  "KWRS" = "East"  # Kingston

)

#Unassigned stations
unique(dat_combined$station)[which(unique(dat_combined$station) %in% names(stations_to_keep) == FALSE)]
dat_combined <- subset(dat_combined, station %in% names(stations_to_keep) & area == 1)

cmmn_coordinates <- read.csv("../data/locations/CMMN_locations.csv")
colnames(cmmn_coordinates) <- c("results_code","station","name","statprov_code","lat","lon")
cmmn_coordinates$station <- as.character(cmmn_coordinates$station)

# Coordinates of US migration stations
us_coordinates <- rbind(data.frame(station = "MCCS", lat = 41.91, lon = -70.54),
                        data.frame(station = "AIMS", lat = 42.99, lon = -70.61),
                        data.frame(station = "KWRS", lat = 41.48, lon = -71.53),
                        data.frame(station = "BIBS", lat = 41.21, lon = -71.58),
                        data.frame(station = "BSBO", lat = 41.61, lon = -83.19),
                        data.frame(station = "FBBO", lat = 39.20, lon = -76.06),
                        data.frame(station = "PARC", lat = 40.16, lon = -79.27),
                        data.frame(station = "CFMS", lat = 64.86, lon = -147.74))

station_coordinates <- rbind(us_coordinates, cmmn_coordinates[,c("station","lat","lon")])
station_coordinates[station_coordinates == "NULL"] <- NA
station_coordinates <- subset(station_coordinates, station %in% dat_combined$station)
station_sf <-  st_as_sf(na.omit(station_coordinates), coords = c("lon", "lat"),crs = 4326, agr = "constant")

range_zones_map_sites <- range_zones_map + 
  geom_sf_text(data = na.omit(station_sf),
               aes(label = station), size = 2, col = "black") +
  ggtitle("Sites with fall counts")

tiff(filename = "../figures/range_zones_map_sites.tiff", width = 6, height = 6, unit = "in", res = 300)
print(range_zones_map_sites)
dev.off()

#-------------------------------------------------------------------------------------------------------
# Part 3: Read in daily count data
#-------------------------------------------------------------------------------------------------------

dat_combined$station_number <- as.numeric(factor(dat_combined$station, levels = names(stations_to_keep)))

dat_combined <- subset(dat_combined, YearCollected < 2019)

min_date <- min(dat_combined$doy, na.rm = TRUE)
min_year <- min(dat_combined$YearCollected, na.rm = TRUE)

dat_combined$doy_adj <- dat_combined$doy - min_date + 1
dat_combined$year_adj <- dat_combined$YearCollected - min_year + 1

nday <- max(dat_combined$doy_adj)
nstation <- max(dat_combined$station_number)
nyear <- max(dat_combined$year_adj)

# Create a nday x nstation x nyear array to store daily counts and effort offsets
daily.count <- array(NA, dim = c(nday,nstation,nyear))
daily.offsets <- array(NA, dim = c(nday,nstation,nyear))
for (i in 1:nrow(dat_combined)){
  daily.count[dat_combined$doy_adj[i],dat_combined$station_number[i],dat_combined$year_adj[i]] = dat_combined$ObservationCount[i]
  daily.offsets[dat_combined$doy_adj[i],dat_combined$station_number[i],dat_combined$year_adj[i]] = dat_combined$net.hrs[i]
  
}
daily.offsets[is.na(daily.offsets)] <- median(daily.offsets, na.rm=TRUE)


#-------------------------------------------------------------------------------------------------------
# Part 4: Assign stations to catchment zones
#-------------------------------------------------------------------------------------------------------
nregion <- 3

# Station indices
station.summary = unique(dat_combined[,c("station","station_number")]) %>% arrange(. , station_number)
station.indices = station.summary$station_number
names(station.indices) <- station.summary$station

# Western sites
col_pal <- RColorBrewer::brewer.pal(length(unique(bcr1$region)),"Set2")
col_pal[c(2:4,5)] <- "gray90"

west_sites <- ggplot() +   theme_bw() +
  geom_sf(data = bcr1, aes(fill = region), col = "gray95")+
  scale_fill_manual(values=col_pal) +
  
  geom_sf_text(data = na.omit(subset(station_sf, station %in% names(stations_to_keep[stations_to_keep == "West"]))),
               aes(label = station), size = 3, col = "black") +
  ggtitle("Stations sampling western zone")

tiff(filename = "../figures/west_sites.tiff", width = 6, height = 6, unit = "in", res = 300)
print(west_sites)
dev.off()

# Central sites
col_pal <- RColorBrewer::brewer.pal(length(unique(bcr1$region)),"Set2")
col_pal[c(1,3:4,5)] <- "gray90"
central_sites <- ggplot() +   theme_bw() +
  geom_sf(data = bcr1, aes(fill = region), col = "gray95")+
  scale_fill_manual(values=col_pal) +
  geom_sf_text(data = na.omit(subset(station_sf, station %in% names(stations_to_keep[stations_to_keep == "Central"]))),
               aes(label = station), size = 3, col = "black") +
  ggtitle("Stations sampling central zone")

tiff(filename = "../figures/central_sites.tiff", width = 6, height = 6, unit = "in", res = 300)
print(central_sites)
dev.off()

# East sites
col_pal <- RColorBrewer::brewer.pal(length(unique(bcr1$region)),"Set2")
col_pal[c(1:2,4,5)] <- "gray90"
east_sites <- ggplot() +   theme_bw() +
  geom_sf(data = bcr1, aes(fill = region), col = "gray95")+
  scale_fill_manual(values=col_pal) +
  
  geom_sf_text(data = na.omit(subset(station_sf, station %in% names(stations_to_keep[stations_to_keep == "East"]))),
               aes(label = station), size = 3, col = "black") +
  ggtitle("Stations sampling east zone")

tiff(filename = "../figures/east_sites.tiff", width = 6, height = 6, unit = "in", res = 300)
print(east_sites)
dev.off()

# [From region, to station, in year]
N.isotope = array(0, dim = c(nregion,nstation,nyear))

# Region 1: "Western" stations
N.isotope[1,which(station.summary$station == "CFMS"),] = 25
N.isotope[1,which(station.summary$station == "TLBBS"),] = 25
N.isotope[1,which(station.summary$station == "MNO"),] = 25
N.isotope[1,which(station.summary$station == "LMBO"),] = 25
N.isotope[1,which(station.summary$station == "TCBO"),] = 25
N.isotope[1,which(station.summary$station == "LPBO"),] = 25


# Region 2: Central stations
N.isotope[2,which(station.summary$station == "PEPBO"),] = 25
N.isotope[2,which(station.summary$station == "BPBO"),] = 25
N.isotope[2,which(station.summary$station == "PIBO"),] = 25
N.isotope[2,which(station.summary$station == "BSBO"),] = 25
N.isotope[2,which(station.summary$station == "TTPBRS"),] = 25
N.isotope[2,which(station.summary$station == "RUTH"),] = 25
N.isotope[2,which(station.summary$station == "MCCS"),] = 25
N.isotope[2,which(station.summary$station == "KWRS"),] = 25
N.isotope[2,which(station.summary$station == "BIBS"),] = 25

# Region 3: Eastern stations
N.isotope[3,which(station.summary$station == "MGBO"),] = 25
N.isotope[3,which(station.summary$station == "FBBO"),] = 25
N.isotope[3,which(station.summary$station == "PARC"),] = 25

N.station.sampled <- apply(N.isotope, c(2,3), FUN = sum) # Number of birds sampled each year

# Remove all isotope data except for 5th year
N.isotope[,1:nstation,1:4] = NA
N.isotope[,1:nstation,6:nyear] = NA

#-------------------------------------------------------------------------------------------------------
# Part 5: Model
#-------------------------------------------------------------------------------------------------------
sink("cmmn_analysis.jags")
cat("
    
    model {
      
      #---------------------------------------------
      # Model for population dynamics in each region
      #---------------------------------------------
 
      for (r in 1:nregion){
      
        # Regional trends
        trend[r] ~ dnorm(0,0.1) 
        
        
      }
      
      # Temporal variance in trend
      proc.sd ~ dunif(0,2)
      proc.tau <- pow(proc.sd,-2)
      
      
      # True (unobserved) population dynamics in each region
      for (y in 1:nyear){
        for (r in 1:nregion){
        
          # Exponential population model
          proc.noise[r,y] ~ dnorm(0,proc.tau)
          logN[r,y] <- logN0[r] + trend[r] * (y-1) + proc.noise[r,y]
          N[r,y] <- exp(logN[r,y])
        }
      }
      
      #---------------------------------------------
      # Observed numbers of birds at station [s] originating from region [r] in year [y]
      #---------------------------------------------
      
      # Proportion of birds from region [r] that are 
      # captured by station [s] (constant through time)
       for (r in 1:nregion){
         for (s in 1:nstation){
         
           log.rho[r,s] ~ dnorm(0,0.1)
           rho[r,s] <- exp(log.rho[r,s]) * rho.fix[r,s]
           
        }
       }
      
      for (s in 1:nstation){
      
      station.sd[s] ~ dunif(0,2)
      station.tau[s] <- pow(station.sd[s],-2)
      
      
        for (y in 1:nyear){
          for (r in 1:nregion){
          
            # Seasonal total arriving at station [s] from region [r]
            log.mu[r,s,y] <- logN[r,y] + log.rho[r,s]
            mu[r,s,y] <- exp(log.mu[r,s,y])
            
          }
        
          # Total seasonal count at station [s]
          mu.total[s,y] <- sum(mu[1:nregion,s,y]) 
          log.mu.total[s,y] <- log(mu.total[s,y])
          
          # Noise in migration strength year
          log.total[s,y] ~ dnorm(log.mu.total[s,y], station.tau[s])
          total[s,y] <- exp(log.total[s,y])
          
          # Multinomial to describe regional composition in this year, at this station
          N.isotope[1:nregion,s,y] ~ dmulti(mu[,s,y], N.station.sampled[s,y])
          
        }
      }
      
      #---------------------------------------------
      # Estimate totals at each station each year
      #---------------------------------------------
      
      for (s in 1:nstation){
      
        mean.migrate[s] ~ dunif(1,nday)
        sd.migrate[s] ~ dunif(0,20)
      
        daily.noise.sd[s] ~ dunif(0,2)
        daily.noise.tau[s] <- pow(daily.noise.sd[s],-2)
      
        for (y in 1:nyear){
          for (d in 1:nday){
            
            norm.density[d,s,y] <- 1/(sqrt(2*pi)*sd.migrate[s])*
                exp(-((d-mean.migrate[s])^2/(2*sd.migrate[s]^2)))
            
            # Expected count on each day
            expected.count[d,s,y] <- norm.density[d,s,y] * total[s,y] * exp(log.offset[d,s,y])
            
            # Daily observation error
            daily.noise[d,s,y] ~ dnorm(0,daily.noise.tau[s])
            
            lambda[d,s,y] <- exp(log(expected.count[d,s,y]) + daily.noise[d,s,y])
            daily.count[d,s,y] ~ dpois(lambda[d,s,y])
            
          }
        }
      }
      
      #---------------------------------------------
      # Derived quantities
      #---------------------------------------------
      
      # Total birds across all regions (used to determine national trend)
      for (y in 1:nyear){
        N.total[y] <- sum(N[,y])
      }
      
    }
    ",fill = TRUE)
sink()

#-------------------------------------------------------------------------------------------------------
# Part 6: Fit model
#-------------------------------------------------------------------------------------------------------

parameters.to.save = c("trend",
                       "proc.sd",
                       "rho",
                       
                       "total",
                       
                       "mean.migrate",
                       "sd.migrate",
                       "daily.noise.sd",

                       "N",
                       
                       "N.total"
                       )

#----------
# Package data for analysis
#----------
rho.fix <- (apply(N.isotope,c(1,2),sum, na.rm = TRUE) > 0) * 1 #In cases where a transition was never observed, fix it to zero

#----------
# Specify initial values
#----------
rho.daily = daily.count / daily.offsets
apply(rho.daily,2,mean, na.rm = TRUE) * dim(rho.daily)[1]

rho.station.init = apply(rho.daily,2,mean, na.rm = TRUE) * dim(rho.daily)[1]
rho.variable.init <- matrix(rep(rho.station.init,nstation),nrow = nregion,ncol = nstation, byrow = TRUE)

inits <- function() list(trend = rep(0,nregion))

jags.data = list(
  
  # Regional regions and initial abundances
  nregion = nregion,
  logN0 = rep(log(1),nregion),
  
  # Number of stations in dataset
  nstation = nstation,
  
  # Number of years in dataset
  nyear = nyear,
  
  # Number of days in each season
  nday = nday,
  
  # Isotope/catchment data
  N.isotope = N.isotope,
  N.station.sampled = N.station.sampled,
  
  daily.count = daily.count,
  
  log.offset = log(daily.offsets),
  pi = pi,
  
  rho.fix = rho.fix, #Can set certain transitions to zero if necessary
  
  max.rho = rho.variable.init[1,]*10
)

#----------
# Fit model in JAGS
#----------
out <- jags(data = jags.data,
            model.file = "cmmn_analysis.jags",
            parameters.to.save = parameters.to.save,
            inits = inits,
            n.chains = 2,
            n.thin = 5,
            n.iter = 40000,
            n.burnin = 10000)

print(out)
max(unlist(out$Rhat),na.rm=TRUE)

#Residual checks: 

    # lambda vs daily count
    # expected count vs daily count
    # true totals vs observed totals
    # true totals over time
    # logN over time
    # trends

#-------------------------------------------------------------------------------------------------------
# Part 7: Summarize posterior densities
#-------------------------------------------------------------------------------------------------------
station.region <- data.frame(region = stations_to_keep,
                             station = names(stations_to_keep))
station.region$region <- factor(station.region$region, levels = c("West","Central","East"))
pd <- out$sims.list


dim(pd$total)

stotal.med <- melt(apply(pd$total, c(2,3), median)) %>% rename(station_number = Var1, year_adj = Var2, index.med = value)
stotal.05 <- melt(apply(pd$total, c(2,3), function(x) quantile(x, 0.05))) %>% rename(station_number = Var1, year_adj = Var2, index.05 = value)
stotal.95 <- melt(apply(pd$total, c(2,3), function(x) quantile(x, 0.95))) %>% rename(station_number = Var1, year_adj = Var2, index.95 = value)

stotal <- full_join(stotal.med, stotal.05) %>% full_join(stotal.95) %>% full_join(station.summary)
stotal$year <- stotal$year_adj + min_year - 1

stotal <- full_join(stotal, station.region)

col_pal <- RColorBrewer::brewer.pal(length(unique(bcr1$region)),"Set2")

stotal.plot <- ggplot(data = stotal) +
  geom_ribbon(aes(x = year, ymin = log(index.05), ymax = log(index.95), fill = region), alpha = 0.5)+
  geom_line(aes(x = year, y = log(index.med), col = region))+
  theme_bw()+
  facet_wrap(region~station, scales = "free")+
  scale_color_manual(values = col_pal, guide = FALSE)+
  scale_fill_manual(values = col_pal, guide = FALSE)+
  ylab("log index")

print(stotal.plot)


rtotal.med <- melt(apply(pd$N, c(2,3), median)) %>% rename(region = Var1, year = Var2, index.med = value)
rtotal.05 <- melt(apply(pd$N, c(2,3), function(x) quantile(x, 0.05))) %>% rename(region = Var1, year = Var2, index.05 = value)
rtotal.95 <- melt(apply(pd$N, c(2,3), function(x) quantile(x, 0.95))) %>% rename(region = Var1, year = Var2, index.95 = value)

rtotal <- full_join(rtotal.med, rtotal.05) %>% full_join(rtotal.95)

rtotal.plot <- ggplot(data = rtotal) +
  geom_ribbon(aes(x = year, ymin = log(index.05), ymax = log(index.95)), alpha = 0.2)+
  geom_line(aes(x = year, y = log(index.med)))+
  theme_bw()+
  facet_wrap(region~.)

print(rtotal.plot)


#-------------------------------------------------------------------------------------------------------
# Part 8: Extract BAM relative abundance (density) estimates for each region
#-------------------------------------------------------------------------------------------------------

# bam map
GDALinfo("../data/relative_abundance_maps/mosaic-BLPW-run3.tif")
bam1 <- raster("../data/relative_abundance_maps/mosaic-BLPW-run3.tif")

# Combine regional polygons
poly.west <- st_combine(subset(bcr1,region == "West")) %>% as("Spatial")
poly.central <- st_combine(subset(bcr1,region == "Central")) %>% as("Spatial")
poly.east <- st_combine(subset(bcr1,region == "East")) %>% as("Spatial")
poly.north <- st_combine(subset(bcr1,region == "Far North")) %>% as("Spatial")
poly.south <- st_combine(subset(bcr1,region == "South")) %>% as("Spatial")

sum.west <- extract(bam1, poly.west, fun = sum, na.rm = TRUE)
sum.central <- extract(bam1, poly.central, fun = sum, na.rm = TRUE)
sum.east <- extract(bam1, poly.east, fun = sum, na.rm = TRUE)
sum.north <- extract(bam1, poly.north, fun = sum, na.rm = TRUE)
sum.south <- extract(bam1, poly.south, fun = sum, na.rm = TRUE)

#-------------------------------------------------------------------------------------------------------
# Part 9: Adjust regional population trajectories based on regional abundances
#-------------------------------------------------------------------------------------------------------

logN.west <- log(sum.west)
logN.central <- log(sum.central)
logN.east <- log(sum.east)

N.adj <- pd$N

N.adj[,1,] <- N.adj[,1,] * as.numeric(sum.west)
N.adj[,2,] <- N.adj[,2,] * as.numeric(sum.central)
N.adj[,3,] <- N.adj[,3,] * as.numeric(sum.east)


rtotal.adj.med <- melt(apply(N.adj, c(2,3), median)) %>% rename(region = Var1, year = Var2, index.med = value)
rtotal.adj.05 <- melt(apply(N.adj, c(2,3), function(x) quantile(x, 0.05))) %>% rename(region = Var1, year = Var2, index.05 = value)
rtotal.adj.95 <- melt(apply(N.adj, c(2,3), function(x) quantile(x, 0.95))) %>% rename(region = Var1, year = Var2, index.95 = value)

rtotal.adj <- full_join(rtotal.adj.med, rtotal.adj.05) %>% full_join(rtotal.adj.95)

rtotal.adj$region.name <- "West"
rtotal.adj$region.name[rtotal.adj$region == 2] <- "Central"
rtotal.adj$region.name[rtotal.adj$region == 3] <- "East"
rtotal.adj$region.name <- factor(rtotal.adj$region.name, levels = c("West","Central","East"))
rtotal.adj.plot <- ggplot(data = rtotal.adj) +
  geom_ribbon(aes(x = year, ymin = log(index.05), ymax = log(index.95)), alpha = 0.2)+
  geom_line(aes(x = year, y = log(index.med)))+
  theme_bw()+
  facet_wrap(region.name~.)

print(rtotal.adj.plot)


rtotal.adj.plot <- ggplot(data = rtotal.adj) +
  geom_ribbon(aes(x = year, ymin = index.05, ymax = index.95), alpha = 0.2)+
  geom_line(aes(x = year, y = index.med))+
  theme_bw()+
  ylab("Relative Abundance")+
  facet_wrap(region.name~.)

print(rtotal.adj.plot)
#-------------------------------------------------------------------------------------------------------
# Part 10: Calculate national totals
#-------------------------------------------------------------------------------------------------------
N.total <- apply(N.adj,c(1,3),sum)


N.total.df <- data.frame(year_adj = 1:nyear,
                         year = (1:nyear) + min_year - 1,
                         N.total.med = apply(N.total, 2, median), 
                         N.total.05 = apply(N.total, 2, function(x) quantile(x, 0.05)), 
                         N.total.95 = apply(N.total, 2, function(x) quantile(x, 0.95)))

N.total.plot <- ggplot(data = N.total.df) +
  geom_ribbon(aes(x = year, ymin = N.total.05, ymax = N.total.95), alpha = 0.2)+
  geom_line(aes(x = year, y = N.total.med))+
  ylab("Relative Abundance")+
  theme_bw()

print(N.total.plot)


logN.total.plot <- ggplot(data = N.total.df) +
  geom_ribbon(aes(x = year, ymin = log(N.total.05), ymax = log(N.total.95)), alpha = 0.2)+
  geom_line(aes(x = year, y = log(N.total.med)))+
  ylab("log(national index)")+
  xlab("year")+
  theme_bw()

print(logN.total.plot)