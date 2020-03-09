
rm(list=ls())

# Describe landscape as aggregate of many cells with independent dynamics
ncell <- 1000
nyear <- 20

# define initial abundances in each cell
x <- seq(0,3,length.out = ncell)
cell.intercept <- log(20)
plot(cell.intercept~x)

# define trend in each cell
y <- seq(0,5,length.out = ncell)
cell.trend <- sin(y) * 0.1 - 0.05
plot(cell.trend~y)

# define 3 strengths of year effects
ye1.sd <- seq(0.1,0,length.out = ncell)
ye2.sd <- seq(0.1,0.1,length.out = ncell)
ye3.sd <- seq(0,0.2,length.out = ncell)



# matrix to store time series of abundance in each cell
mumat <- matrix(NA, nrow = nyear, ncol = ncell)
Nmat <- mumat
for (year in 1:nyear){
  mumat[year,] <- cell.intercept + cell.trend*year + rnorm(ncell,0,ye1.sd) + rnorm(ncell,0,ye2.sd) + rnorm(ncell,0,ye3.sd)
  Nmat[year,] <- rpois(ncell, exp(mumat[year,]))
}

plot(Nmat[1,])
plot(Nmat[nyear,])

# National totals each year
yvec <- 1:nyear
Ntotal <- apply(Nmat, 1, sum)
plot(Ntotal~yvec, type = "l")



#---------------
# Migration process
#---------------

# Define number of stations
nstation <- 6

# Define probabilities each bird passes each station
mig.prob <- matrix(NA, nrow = nstation, ncol = ncell)

for (c in 1:ncell){
  
  for (s in 1:nstation){
    
    prob <- (1 - abs(c/ncell - s/nstation))^3 * 0.02
    
    mig.prob[s,c] <- prob #prob bird from cell c migrates past station s
    
  }
  
}

plot(mig.prob[3,])


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

matplot(obsmat.sum, type = "l")


#---------------
# Sample 25 birds in year 5 from each station (stable isotopes)
#---------------

bird_sample <- matrix(NA, nrow = 25, ncol = nstation)

for (s in 1:nstation){
  birds_at_station <- rep(which(obsmat[5,s,] != 0), obsmat[5,s,][which(obsmat[5,s,] != 0)])
  bird_sample[,s] <- sort(sample(birds_at_station, 25))
}


#---------------
# Assign discrete regions to each cell
#---------------

regions <- rep(1,ncell)
regions[301:600] <- 2
regions[601:1000] <- 3

# Assign each sampled bird to a region

bird_region <- bird_sample * 0 + 1
bird_region[bird_sample >= 301 & bird_sample <= 600] = 2
bird_region[bird_sample >600] = 3






