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

#---------------------------------------------------------------------------------------------
# Load and plot raw data
#---------------------------------------------------------------------------------------------
dat = read.csv(file = "PEPBO_example.csv")

# Plot the raw daily counts and save as PDF
plot1 <- ggplot(data = dat) +
  geom_point(aes(x = doy, y = ObservationCount), col = "blue")+
  facet_grid(station~YearCollected, scales = "free_y")+
  xlab("Day of Year")+
  ylab("Count")+
  theme_bw()
pdf(file = "plot1.pdf", width = 30,height=4)
print(plot1)
dev.off()

#---------------------------------------------------------------------------------------------
# Bayesian analysis (separate for each station)
#---------------------------------------------------------------------------------------------

# Model with random effects for annual phenology
sink("cmmn_separate.jags")
cat("
    
    model {
      
      #---------------------------------------------
      # Model for population dynamics
      #---------------------------------------------

      # Trend and annual noise
      log.trend ~ dnorm(0,1)
      proc.sd ~ dunif(0,2)
      proc.tau <- pow(proc.sd,-2)
      
      intercept ~ dunif(0,upper_limit)
      log.intercept <- log(intercept)
        
      # Population trend model
      for (y in 1:nyear){
        
        mu[y] <- log.intercept + log.trend*(y-1)
        logN[y] <- mu[y] + noise[y]
        noise[y] ~ dnorm(0,proc.tau)
        N[y] <- exp(logN[y])

      } # close year loop

      #---------------------------------------------
      # Model for daily counts
      #---------------------------------------------

      # Parameters describing the mean date of migration, and variation in that peak date among years
      mean.migrate.HYPERMEAN ~ dunif(1,nday)
      mean.migrate.HYPERSD ~ dunif(0,nday)
      mean.migrate.HYPERTAU <- pow(mean.migrate.HYPERSD, -2)
      
      # Parameter describing the width of the migration window (assume this window is constant)
      sd.migrate ~ dunif(0,nday)
         
      # Magnitude of dailyobservation error
      daily.noise.sd ~ dunif(0,2)
      daily.noise.tau <- pow(daily.noise.sd,-2)
      
        for (y in 1:nyear){
        
        #-------------------------------------------
        # Different models of peak phenology (either variable among years or fixed among years)
        #-------------------------------------------
        # Different migration periods for each year (drawn from a shared distribution)
        #mean.migrate[y] ~ dnorm(mean.migrate.HYPERMEAN, mean.migrate.HYPERTAU)
        
        # Force phenology to be the same in all years
        mean.migrate[y] <- mean.migrate.HYPERMEAN
        #-------------------------------------------
        
        for (d in 1:nday){
              
            # Each day within each year, estimate the proportion of total annual detections.
            norm.density[d,y] <- 1/(sqrt(2*pi)*sd.migrate)*exp(-((d-mean.migrate[y])^2/(2*sd.migrate^2)))
              
            # Expected count on each day (probability density * total abundance)
            expected.count[d,y] <- norm.density[d,y] * N[y]
              
            # Daily observation error (some days have higher/lower counts than expected)
            daily.noise[d,y] ~ dnorm(0, daily.noise.tau)
            log.lambda[d,y] <- log(expected.count[d,y]) + daily.noise[d,y]
              
          } # close day loop
            
      } # close year loop 
      
      
      # Observation model
      for (i in 1:nobs){
        lam[i] <- exp(log.lambda[day[i],year[i]])
        daily.count[i] ~ dpois(lam[i])
      }
      
      
    }
    ",fill = TRUE)
sink()

#---------------------------------------------------------------------------------------------
# Package data for JAGS
#---------------------------------------------------------------------------------------------

# Clean data a bit more before
dat$doy_adjusted = dat$doy - min(dat$doy) + 1
dat$year_adjusted = dat$YearCollected - min(dat$YearCollected) + 1
dat = subset(dat, !is.na(ObservationCount))


jags.data = list(daily.count = dat$ObservationCount,
                 nobs = nrow(dat),
                 
                 day = dat$doy_adjusted,
                 nday = max(dat$doy_adjusted),
                 
                 year = dat$year_adjusted,
                 nyear = max(dat$year_adjusted),
                 
                 upper_limit = max(aggregate(ObservationCount~year_adjusted, data = dat, FUN = sum)$ObservationCount)*10,
                 
                 pi = pi
)

inits <- function() list(log.trend = rnorm(1,0,0.05),
                         proc.sd = runif(1,0.1,0.5))

parameters.to.save = c("log.trend",
                       "proc.sd",
                       
                       "log.intercept",
                       "daily.noise.sd",
                       "mean.migrate.HYPERMEAN",
                       "mean.migrate.HYPERSD",
                       
                       "mean.migrate",
                       "sd.migrate",
                       
                       "expected.count",
                       "log.lambda",
                       "mu",
                       "N"
                       
)

#---------------------------------------------------------------------------------------------
# Fit model in JAGS
#---------------------------------------------------------------------------------------------

out <- jags(data = jags.data,
            model.file = "cmmn_separate.jags",
            parameters.to.save = parameters.to.save,
            inits = inits,
            n.chains = 3,
            n.thin = 5,
            n.iter = 20000,
            n.burnin = 10000)

# Posterior samples
ps <- out$sims.list

# Convergence statistics
max.Rhat = max(unlist(out$Rhat),na.rm=TRUE) # Ideally less than 1.1
#traceplot(out)

#---------------------------------------------------------------------------------------------
# Extract and plot results
#---------------------------------------------------------------------------------------------

# Estimate of temporal trend (95% CRI and median) and associated histogram
quantile(ps$log.trend, c(0.025, 0.5 ,0.975))
hist(ps$log.trend, breaks = 200, col = "gray80")

#-------------------------
# Expected annual totals
#-------------------------
# Store and plot time series of annual index estimates at this station
N.est = data.frame(station = dat$station[1],
                              season = dat$season[1],
                              
                              year = (1:jags.data$nyear) + min(dat$YearCollected) - 1,
                              
                              index.500 = apply(ps$N,2,function(x) quantile(x,0.500)),
                              index.025 = apply(ps$N,2,function(x) quantile(x,0.025)),
                              index.975 = apply(ps$N,2,function(x) quantile(x,0.975)))

plot2 <- ggplot(data = N.est) +
  geom_ribbon(aes(x = year, ymin = index.025, ymax = index.975), alpha = 0.2)+
  geom_line(aes(x = year, y = index.500), col = "black")+
  xlab("Year")+
  ylab("Expected seasonal total")+
  ggtitle(dat$station_season)+
  theme_bw()
pdf(file = "plot2.pdf", width = 6, height= 3)
print(plot2)
dev.off()

plot3 <- ggplot(data = N.est) +
  geom_ribbon(aes(x = year, ymin = log(index.025), ymax = log(index.975)), alpha = 0.2)+
  geom_line(aes(x = year, y = log(index.500)), col = "black")+
  xlab("Year")+
  ylab("Annual index\n(log expected seasonal total)")+
  ggtitle(dat$station_season)+
  theme_bw()
pdf(file = "plot3.pdf", width = 6, height= 3)
print(plot3)
dev.off()

#-------------------------
# Expected daily counts
#-------------------------
# Store and plot daily expected counts in each year
daily.est <- expand.grid(doy_adjusted = 1:dim(ps$expected.count)[2],
                        year_adjusted = 1:dim(ps$expected.count)[3],
                        expected.025 = NA, expected.500 = NA, expected.975 = NA)
daily.est$year <- daily.est$year_adjusted + min(dat$YearCollected) - 1
daily.est$doy <- daily.est$doy_adjusted + min(dat$doy) - 1

for (i in 1:nrow(daily.est)){
  
  #Extract posterior estimates for this day, in this year
  ps.daily <- ps$expected.count[,daily.est$doy_adjusted[i], daily.est$year_adjusted[i]]
  
  daily.est$expected.025[i] <- quantile(ps.daily, 0.025)
  daily.est$expected.500[i] <- quantile(ps.daily, 0.500)
  daily.est$expected.975[i] <- quantile(ps.daily, 0.975)
}

plot4 <- ggplot(data = daily.est) +
  geom_ribbon(aes(x = doy, ymin = expected.025, ymax = expected.975), alpha = 0.2)+
  geom_line(aes(x = doy, y = expected.500), col = "black")+
  xlab("Day of Year")+
  ylab("Daily index\n(expected daily count)")+
  facet_wrap(.~year)+
  ggtitle(dat$station_season)+
  theme_bw()

pdf(file = "plot4.pdf", width = 10, height= 6)
print(plot4)
dev.off()

#-------------------------
# Overlay raw daily counts
#-------------------------

for (i in 1:nrow(daily.est)){
  
  dat_for_this_day <- subset(dat, YearCollected == daily.est$year[i] & doy == daily.est$doy[i])
  if (nrow(dat_for_this_day) != 1) next
  
  daily.est$observed[i] <- dat_for_this_day$ObservationCount
}


plot5 <- ggplot(data = daily.est) +
  geom_ribbon(aes(x = doy, ymin = expected.025, ymax = expected.975), alpha = 0.2)+
  geom_line(aes(x = doy, y = expected.500), col = "black")+
  geom_point(aes(x = doy, y = observed), col = "blue")+
  xlab("Day of Year")+
  ylab("Daily index\n(expected daily count)")+
  facet_wrap(.~year)+
  ggtitle(dat$station_season)+
  theme_bw()

pdf(file = "plot5.pdf", width = 10, height= 6)
print(plot5)
dev.off()


#-------------------------
# Relationship between observed total and estimated annual index
#-------------------------

N.est$observed.total <- aggregate(ObservationCount ~ year_adjusted, data = dat, FUN = sum)$ObservationCount

plot6 <- ggplot(data = N.est) +
  geom_errorbar(aes(x = log(observed.total), ymin = log(index.025), ymax = log(index.975)), width = 0)+
  geom_point(aes(x = log(observed.total), y = log(index.500)), col = "black")+
  xlab("log(total seasonal count)")+
  ylab("Population index")+
  ggtitle(dat$station_season)+
  theme_bw()

pdf(file = "plot6.pdf", width = 5, height= 5)
print(plot6)
dev.off()
