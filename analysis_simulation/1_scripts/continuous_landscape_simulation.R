
rm(list=ls())

# Required packages
my.packs <- c(
  'tidyverse',   # For data manipulation/plotting
  'reshape2',    # For data manipulation
  'jagsUI',      # For Bayesian analysis
  'lme4',        # For Frequentist analysis
  'viridis',     # For plotting
  'cowplot',#,     # For plotting
  #'optimx'       # For fitting models
  'parallel',
  'doParallel'
)

# if any of them are not installed, install them
if (any(!my.packs %in% installed.packages()[, 'Package'])) {install.packages(my.packs[which(!my.packs %in% installed.packages()[, 'Package'])],dependencies = TRUE)}

lapply(my.packs, require, character.only = TRUE) # Load packages

# Describe landscape as aggregate of many cells with independent dynamics
ncell <- 1000
nyear <- 20

#---------------
# Assign discrete regions to each cell
#---------------
nregion <- 3
regions <- rep(1,ncell)
regions[201:600] <- 2
regions[601:1000] <- 3


simulation_results <- data.frame()

#for (iteration in 1:1000){

#----------------------
# Parallel processing setup
numCores <- detectCores() # Detect number of cores on machine
numCores <- numCores - 2  # Reserve 2 cores for other tasks
registerDoParallel(numCores) # Register cores to use for parallel processing
#----------------------

all_results = foreach(iteration = 1:1000, .combine = rbind, .packages = my.packs) %dopar% {
  
  set.seed(iteration)
  
  # define initial abundances in each cell
  x <- seq(0,3,length.out = ncell)
  cell.intercept <- log(20)
  
  # define trend in each cell
  y <- seq(0,5,length.out = ncell)
  #cell.trend <- sin(y) * 0.1 - 0.05
  
  cell.trend <- rep(NA, ncell)
  cell.trend[regions == 1] = rnorm(sum(regions==1),0.05,0)
  cell.trend[regions == 2] = rnorm(sum(regions==2),0,0)
  cell.trend[regions == 3] = rnorm(sum(regions==3),-0.05,0)
  
  
  # matrix to store time series of abundance in each cell
  mumat <- matrix(NA, nrow = nyear, ncol = ncell)
  Nmat <- mumat
  for (year in 1:nyear){
    ye1 <- rnorm(ncell,0,0.2* (regions == 1))
    ye2 <- rnorm(ncell,0,0.2* (regions == 2))
    ye3 <- rnorm(ncell,0,0.2* (regions == 3))
    
    mumat[year,] <- cell.intercept + cell.trend*year + ye1 + ye2 + ye3
    Nmat[year,] <- rpois(ncell, exp(mumat[year,]))
  }
  
  
  # National totals each year
  yvec <- 1:nyear
  Ntotal <- apply(Nmat, 1, sum)
  
  #---------------
  # Migration process
  #---------------
  
  # Define number of stations
  nstation <- 6
  
  # Define probabilities each bird passes each station
  mig.prob <- matrix(NA, nrow = nstation, ncol = ncell)
  
  for (c in 1:ncell){
    
    for (s in 1:nstation){
      
      prob <- (1 - abs(c/ncell - s/nstation))^5 * 0.02
      
      mig.prob[s,c] <- prob #prob bird from cell c migrates past station s
      
    }
    
  }
  
  
  # Actual number of birds recorded at each station in each year from each cell
  obsmat <- array(NA, dim = c(nyear,nstation,ncell))
  
  for (s in 1:nstation){
    for (y in 1:nyear){
      obsmat[y,s,] <- rbinom(length(Nmat[y,]),Nmat[y,],mig.prob[s,])
      
    }
  }
  
  
  # Annual totals at each station
  obsmat.sum <- apply(obsmat, c(1,2), sum)
  obsmat.sum
  
  #---------------
  # Sample 25 birds in year 5 from each station (stable isotopes)
  #---------------
  
  bird_sample <- matrix(NA, nrow = 25, ncol = nstation)
  
  for (s in 1:nstation){
    birds_at_station <- rep(which(obsmat[5,s,] != 0), obsmat[5,s,][which(obsmat[5,s,] != 0)])
    bird_sample[,s] <- sort(sample(birds_at_station, 25))
  }
  
  # Assign each sampled bird to a region
  bird_region <- bird_sample * 0 + 1
  bird_region[bird_sample >= 301 & bird_sample <= 600] = 2
  bird_region[bird_sample >600] = 3
  
  # Recast into appropriate form for analysis (r x s x y array)
  
  N.isotope <- array(NA, dim = c(nregion, nstation, nyear))
  
  for (r in 1:nregion){
    for (s in 1:nstation){
      
      N.isotope[r,s,5] <- sum(bird_region[,s] == r)
      
    }
  }
  
  N.station.sampled <- apply(N.isotope, c(2,3), FUN = sum) # Number of birds sampled each year
  N.station.sampled[is.na(N.station.sampled)] = 25
  
  #---------------
  # Prepare data for analysis
  #---------------
  
  nregion <- 3
  nstation <- 6
  
  sink("regional_analysis.jags")
  cat("
    
    model {
      
      #---------------------------------------------
      # Model for population dynamics in each region
      #---------------------------------------------
 
      for (r in 1:nregion){
      
        # Regional trends
        trend[r] ~ dnorm(0,0.1) 
      
      }
      
      # Temporal variance in regional trend (assumed to be same magnitude of variance in all regions)
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
      
        # Year to year variance in abundance at station
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
          log.Nobs[s,y] ~ dnorm(log.mu.total[s,y], station.tau[s]) # Observations
          
          # Multinomial to describe regional composition in this year, at this station
          N.isotope[1:nregion,s,y] ~ dmulti(mu[,s,y], N.station.sampled[s,y])
          
          
          # Noise in migration strength year
          #log.total[s,y] ~ dnorm(log.mu.total[s,y], station.tau[s])
          #total[s,y] <- exp(log.total[s,y])
          #N.obs[s,y] ~ dnorm(total[s,y], 0.1)
          
        }
      }
      # 
      # #---------------------------------------------
      # # Estimate totals at each station each year
      # #---------------------------------------------
      # 
      # for (s in 1:nstation){
      # 
      #   mean.migrate[s] ~ dunif(1,nday)
      #   sd.migrate[s] ~ dunif(0,20)
      # 
      #   daily.noise.sd[s] ~ dunif(0,2)
      #   daily.noise.tau[s] <- pow(daily.noise.sd[s],-2)
      # 
      #   for (y in 1:nyear){
      #     for (d in 1:nday){
      #       
      #       norm.density[d,s,y] <- 1/(sqrt(2*pi)*sd.migrate[s])*
      #           exp(-((d-mean.migrate[s])^2/(2*sd.migrate[s]^2)))
      #       
      #       # Expected count on each day
      #       expected.count[d,s,y] <- norm.density[d,s,y] * total[s,y] * exp(log.offset[d,s,y])
      #       
      #       # Daily observation error
      #       daily.noise[d,s,y] ~ dnorm(0,daily.noise.tau[s])
      #       
      #       lambda[d,s,y] <- exp(log(expected.count[d,s,y]) + daily.noise[d,s,y])
      #       daily.count[d,s,y] ~ dpois(lambda[d,s,y])
      #       
      #     }
      #   }
      # }
      
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
  
  rho.fix <- (apply(N.isotope,c(1,2),sum, na.rm = TRUE) > 0) * 1 #In cases where a transition was never observed, fix it to zero
  
  jags.data <- list(#N.obs  = t(obsmat.sum),
                    log.Nobs = log(t(obsmat.sum)),
                    nregion = 3,
                    
                    logN0 = rep(log(1),nregion),
                    
                    # Number of stations in dataset
                    nstation = nstation,
                    
                    # Number of years in dataset
                    nyear = nyear,
                    
                    # Isotope/catchment data
                    N.isotope = N.isotope,
                    N.station.sampled = N.station.sampled,
                    
                    pi = pi,
                    
                    rho.fix = rho.fix #Can set certain transitions to zero if necessary
                    
                    
  )
  
  inits <- function() list(trend = rep(0,nregion),
                           proc.sd = runif(1,0,0.1))
  
  parameters.to.save = c("trend",
                         "proc.sd",
                         "station.sd",
                         "rho",
                         "total",
                         "N",
                         "N.total")
  
  #----------
  # Fit model in JAGS
  #----------
  out <- jags(data = jags.data,
              model.file = "regional_analysis.jags",
              parameters.to.save = parameters.to.save,
              inits = inits,
              n.chains = 3,
              n.thin = 5,
              n.iter = 40000,
              n.burnin = 20000)
  
  #print(out)
  
  #----------
  # Output
  #----------
  pd <- out$sims.list
  
  stotal.med <- melt(apply(pd$total, c(2,3), median)) %>% rename(station_number = Var1, year_adj = Var2, index.med = value)
  stotal.05 <- melt(apply(pd$total, c(2,3), function(x) quantile(x, 0.05))) %>% rename(station_number = Var1, year_adj = Var2, index.05 = value)
  stotal.95 <- melt(apply(pd$total, c(2,3), function(x) quantile(x, 0.95))) %>% rename(station_number = Var1, year_adj = Var2, index.95 = value)
  stotal <- full_join(stotal.med, stotal.05) %>% full_join(stotal.95)
  stotal$year <- stotal$year_adj
  
  col_pal <- RColorBrewer::brewer.pal(nregion,"Set2")
  
  stotal.plot <- ggplot(data = stotal) +
    geom_ribbon(aes(x = year, ymin = log(index.05), ymax = log(index.95)), alpha = 0.5)+
    geom_line(aes(x = year, y = log(index.med)))+
    theme_bw()+
    facet_wrap(station_number~.)+
    scale_color_manual(values = col_pal, guide = FALSE)+
    scale_fill_manual(values = col_pal, guide = FALSE)+
    ylab("log index")
  
  
  
  # Regional indices over time (not weighted by relative abundance)
  rtotal.med <- melt(apply(pd$N, c(2,3), median)) %>% rename(region = Var1, year = Var2, index.med = value)
  rtotal.05 <- melt(apply(pd$N, c(2,3), function(x) quantile(x, 0.05))) %>% rename(region = Var1, year = Var2, index.05 = value)
  rtotal.95 <- melt(apply(pd$N, c(2,3), function(x) quantile(x, 0.95))) %>% rename(region = Var1, year = Var2, index.95 = value)
  
  rtotal <- full_join(rtotal.med, rtotal.05) %>% full_join(rtotal.95)
  
  rtotal.plot <- ggplot(data = rtotal) +
    geom_ribbon(aes(x = year, ymin = log(index.05), ymax = log(index.95)), alpha = 0.2)+
    geom_line(aes(x = year, y = log(index.med)))+
    theme_bw()+
    facet_wrap(region~.)
  
  
  #-------------------------------------------------------------------------------------------------------
  # Part 9: Adjust regional population trajectories based on regional abundances
  #-------------------------------------------------------------------------------------------------------
  
  sum.west <- sum(Nmat[1,regions == 1])
  sum.central <- sum(Nmat[1,regions == 2])
  sum.east <- sum(Nmat[1,regions == 3])
  
  N.adj <- pd$N
  
  N.adj[,1,] <- N.adj[,1,] * as.numeric(sum.west)
  N.adj[,2,] <- N.adj[,2,] * as.numeric(sum.central)
  N.adj[,3,] <- N.adj[,3,] * as.numeric(sum.east)
  
  rtotal.adj.med <- melt(apply(N.adj, c(2,3), median)) %>% rename(region = Var1, year = Var2, index.med = value)
  rtotal.adj.05 <- melt(apply(N.adj, c(2,3), function(x) quantile(x, 0.05))) %>% rename(region = Var1, year = Var2, index.05 = value)
  rtotal.adj.95 <- melt(apply(N.adj, c(2,3), function(x) quantile(x, 0.95))) %>% rename(region = Var1, year = Var2, index.95 = value)
  
  # True regional trends
  rtotal.true <- rbind(data.frame(region = 1, year = 1:nyear, N.true = apply(Nmat[,which(regions == 1)],1,sum)),
                       data.frame(region = 2, year = 1:nyear, N.true = apply(Nmat[,which(regions == 2)],1,sum)),
                       data.frame(region = 3, year = 1:nyear, N.true = apply(Nmat[,which(regions == 3)],1,sum)))
  
  
  rtotal.adj <- full_join(rtotal.adj.med, rtotal.adj.05) %>% full_join(rtotal.adj.95) %>% full_join(rtotal.true)
  
  
  rtotal.adj$region.name <- "West"
  rtotal.adj$region.name[rtotal.adj$region == 2] <- "Central"
  rtotal.adj$region.name[rtotal.adj$region == 3] <- "East"
  rtotal.adj$region.name <- factor(rtotal.adj$region.name, levels = c("West","Central","East"))
  
  
  
  
  
  log.rtotal.adj.plot <- ggplot(data = rtotal.adj) +
    geom_ribbon(aes(x = year, ymin = log(index.05), ymax = log(index.95), fill = region.name), alpha = 0.5)+
    geom_line(aes(x = year, y = log(index.med), col = region.name))+
    geom_line(aes(x = year, y = log(N.true)), col = "blue", size = 1.5)+
    scale_color_manual(values = col_pal, guide = FALSE)+
    scale_fill_manual(values = col_pal, guide = FALSE)+
    theme_bw()+
    ylab("log(Regional Index)")+
    facet_wrap(region.name~.)
  
  rtotal.adj.plot <- ggplot(data = rtotal.adj) +
    geom_ribbon(aes(x = year, ymin = index.05, ymax = index.95, fill = region.name), alpha = 0.5)+
    geom_line(aes(x = year, y = index.med, col = region.name))+
    geom_line(aes(x = year, y = N.true), col = "blue", size = 1.5)+
    scale_color_manual(values = col_pal, guide = FALSE)+
    scale_fill_manual(values = col_pal, guide = FALSE)+
    
    theme_bw()+
    ylab("Regional Index")+
    facet_wrap(region.name~.)
  
  regional.index.plot <- plot_grid(log.rtotal.adj.plot,
                                   rtotal.adj.plot, nrow = 2, align = "hv")
  
  #print(regional.index.plot)
  
  #-------------------------------------------------------------------------------------------------------
  # Part 10: Calculate national totals
  #-------------------------------------------------------------------------------------------------------
  
  N.total <- apply(N.adj,c(1,3),sum)
  
  N.total.df <- data.frame(year_adj = 1:nyear,
                           year = 1:nyear,
                           N.total.med = apply(N.total, 2, median), 
                           N.total.05 = apply(N.total, 2, function(x) quantile(x, 0.05)), 
                           N.total.95 = apply(N.total, 2, function(x) quantile(x, 0.95)),
                           N.true = apply(Nmat,1,sum))
  
  
  logN.total.plot <- ggplot(data = N.total.df) +
    geom_ribbon(aes(x = year, ymin = log(N.total.05), ymax = log(N.total.95)), alpha = 0.2)+
    geom_line(aes(x = year, y = log(N.total.med)))+
    geom_line(aes(x = year, y = log(N.true)), col = "blue", size = 1.5)+
    
    ylab("log(National Index)")+
    xlab("year")+
    theme_bw()
  
  N.total.plot <- ggplot(data = N.total.df) +
    geom_ribbon(aes(x = year, ymin = N.total.05, ymax = N.total.95), alpha = 0.2)+
    geom_line(aes(x = year, y = N.total.med))+
    geom_line(aes(x = year, y = N.true), col = "blue", size = 1.5)+
    
    ylab("National Index")+
    theme_bw()
  
  
  national.index.plot <- plot_grid(logN.total.plot,
                                   N.total.plot, nrow = 2, 
                                   align = "hv")
  
  #print(national.index.plot)
  
  #-------------------------------------------------------------------------------------------------------
  # Part 11: Calculate relative bias in regional and national empirical trends
  #-------------------------------------------------------------------------------------------------------
  
  trend.national.true <- (log(Ntotal[nyear]) - log(Ntotal[1])) / (nyear-1)
  trend.national.est <- (log(N.total[,nyear]) - log(N.total[,1])) / (nyear-1)
  trend.national.RelBias <- mean((trend.national.est - trend.national.true)/trend.national.true)
  
  trend.reg1.true <- ( log(apply(Nmat[,which(regions == 1)],1,sum)[nyear]) - log(apply(Nmat[,which(regions == 1)],1,sum)[1])) / (nyear-1)
  trend.reg1.est <- (log(N.adj[,1,nyear]) - log(N.adj[,1,1])) / (nyear-1)
  trend.reg1.RelBias <- mean((trend.reg1.est - trend.reg1.true)/trend.reg1.true)
  
  trend.reg2.true <- ( log(apply(Nmat[,which(regions == 2)],1,sum)[nyear]) - log(apply(Nmat[,which(regions == 2)],1,sum)[1])) / (nyear-1)
  trend.reg2.est <- (log(N.adj[,2,nyear]) - log(N.adj[,2,1])) / (nyear-1)
  trend.reg2.RelBias <- mean((trend.reg2.est - trend.reg2.true)/trend.reg2.true)
  
  trend.reg3.true <- ( log(apply(Nmat[,which(regions == 3)],1,sum)[nyear]) - log(apply(Nmat[,which(regions == 3)],1,sum)[1])) / (nyear-1)
  trend.reg3.est <- (log(N.adj[,3,nyear]) - log(N.adj[,3,1])) / (nyear-1)
  trend.reg3.RelBias <- mean((trend.reg3.est - trend.reg3.true)/trend.reg3.true)
  
  
  # Check if species_results already exists
  if (file.exists("../2_Routput/continuous_landscape.csv")){
    simulation_results = read.csv("../2_Routput/continuous_landscape.csv")
  }
  
  run_results <- data.frame(iteration = iteration,
                            
                            max.Rhat = max(unlist(out$Rhat),na.rm=TRUE),
                            
                            trend.national.true = trend.national.true,
                            trend.national.est = mean(trend.national.est),
                            trend.national.RelBias = trend.national.RelBias,
                            
                            trend.reg1.true = trend.reg1.true,
                            trend.reg1.est = mean(trend.reg1.est),
                            trend.reg1.RelBias = trend.reg1.RelBias,
                            
                            trend.reg2.true = trend.reg2.true,
                            trend.reg2.est = mean(trend.reg2.est),
                            trend.reg2.RelBias = trend.reg2.RelBias,
                            
                            trend.reg3.true = trend.reg3.true,
                            trend.reg3.est = mean(trend.reg3.est),
                            trend.reg3.RelBias = trend.reg3.RelBias)
  
  simulation_results <- rbind(simulation_results, run_results)
  write.csv(simulation_results, file = "../2_Routput/continuous_landscape.csv", row.names = FALSE)
  
}


# Evaluate bias/precision
simulation_results <- read.csv(file = "../2_Routput/continuous_landscape.csv")

# Subset to only "converged" runs
simulation_results <- subset(simulation_results, max.Rhat <= 1.2)

ggplot(data = simulation_results) +
  geom_point(aes(x = trend.reg1.true, y = trend.reg1.est), col = "blue")+
  geom_point(aes(x = trend.reg2.true, y = trend.reg2.est), col = "green")+
  geom_point(aes(x = trend.reg3.true, y = trend.reg3.est), col = "orangered")+
  
  geom_point(aes(x = trend.national.true, y = trend.national.est), col = "black")+
  
  geom_abline(intercept = 0, slope = 1)+
  coord_equal(xlim=c(-0.1,0.1),ylim=c(-0.1,0.1))+
  
  theme_bw()
