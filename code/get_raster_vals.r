## projection matrix coordinates
locs <- expand.grid(proj$x,proj$y)
## initalize lists
lums <- pops <- list()

## read in covariate values worldwide from rasters (assumes a raster for each year)
## NOTE this assumes that the list of rasters in each folder are in the same order as the split data frame
## (i.e., orderd by year) and give the raw values of the raster images extracted at long&lat locations
## this isn't the most efficient way to do this but fine for now
## luminosity
list.rasters <- list.files("luminosity",pattern = ".tif$",full.names = TRUE)
rasters <- sapply(list.rasters,raster::raster)
for (i in 1:length(rasters)){
    lums[[i]] <- extract(rasters[[i]], locs)
}
## functions for intercalibration based on Elvidge2013 for luminosity (to make yearly data comparable)
fun.lums <- list(f10 <- function(x){ifelse(x>0, 2.3430+0.5102*x+0.0010*x^2, 0)},
                 f11 <- function(x){ifelse(x>0, 1.8956+0.7345*x+0.0010*x^2, 0)},
                 f12 <- function(x){ifelse(x>0, 1.8750+0.6203*x+0.0010*x^2, 0)},
                 f13 <- function(x){ifelse(x>0, 1.8750+0.6203*x+0.0010*x^2, 0)},
                 f14 <- function(x){ifelse(x>0, 1.8750+0.6203*x+0.0010*x^2, 0)},
                 f15 <- function(x){ifelse(x>0, 1.8750+0.6203*x+0.0010*x^2, 0)},
                 f16 <- function(x){ifelse(x>0, 1.8750+0.6203*x+0.0010*x^2, 0)},
                 f17 <- function(x){ifelse(x>0, 1.8750+0.6203*x+0.0010*x^2, 0)})
## intercalibration for luminosity
for(i in 1:length(lums)){lums[[i]] <- fun.lums[[i]](lums[[i]])}
## scale
## luminosity means in year order previously calculated
means.lum <- c(22.88600, 32.17235, 22.55995, 27.53695, 25.93906, 24.73974, 24.04713, 19.97450)
sd.lum <- c(15.28341, 20.42197, 17.60979, 18.03055, 18.22354, 18.52160, 18.39245, 18.55277)
for(i in 1:length(lums)){lums[[i]] <- (lums[[i]] - means.lum[i])/sd.lum[i]}
### population density
list.rasters <- list.files("popdensity",pattern = ".tif$",full.names = TRUE)
rasters <- sapply(list.rasters,raster::raster)
for (i in 1:length(rasters)){
    pops[[i]] <- extract(rasters[[i]], locs)
}
## scale
## population density means in year order previously calculated
means.pop <- c(1476.681, 1267.425, 1277.815, 1699.609, 1605.749, 2374.827, 1496.576, 1291.569)
sd.pop <- c(2634.002, 2544.810, 2910.881, 4185.542, 4594.430, 6254.001, 3630.597, 3322.203)
for(i in 1:length(pops)){pops[[i]] <- (pops[[i]] - means.pop[i])/sd.pop[i]}
## time to city covariate
raster <- list.files("tt",pattern = ".tif$",full.names = TRUE)
tt <- raster(raster)
tt <- extract(tt,locs)
## sort out time to city covariate
tt <- ifelse(tt < 0, NA, tt)
## scale
## time to city scaling parameters previously calculated
mean.tt <- 33.14502
sd.tt <- 103.2259
tt <- (tt - mean.tt)/sd.tt

