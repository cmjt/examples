## required libraries
library(lgcpSPDE) ## install from github by running devtools::install_github("cmjt/lgcpSPDE") in R
library(raster) ## to extract covariates for plotting
library(rgeos)
library(rgdal) ## for nearby countries
#############################################
## control coarseness of the projections
dims <- c(500,500)
## full country names for data we are interested in
countries.full <- c("Afghanistan","Iraq","Pakistan","Iran")
# list of spatial polygons of above countries
sps <- sapply(countries.full, function(x) world[world$name == x,])
## tag to indicate if prediction is carried out (i.e., NA values put in place of data for countries of interest)
PRED <- FALSE
## prediction year if prediction is to be carried out
if(PRED){pred.year <- 2017}


##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
## out-of-sample prediction
nearby.countries <- lapply(sps,function(x) as.character(world$name[gTouches(x,world,byid = TRUE)]))
for(i in countries.full){nearby.countries[[i]] <- c(nearby.countries[[i]],i)}
sp.near <- sapply(nearby.countries,function(x) world[world$name %in% x,])

meshs <- list()

mn <- list(2/180,1/180,2/180,1/180)
max <- list(10/180,5/180,10/180,3/180)
names(mn) <- names(max) <- countries.full ## names mesh resolution lists


## create mesh for each "set" of countries projected onto the unit sphere
for(cont in countries.full){
    bdry <- inla.sp2segment(sp.near[[cont]])
    bdry$loc <- inla.mesh.map(bdry$loc, projection = "longlat",inverse = TRUE)
    meshs[[cont]] <- inla.mesh.2d(boundary = bdry, max.edge = c(mn[[cont]],max[[cont]]),cutoff = mn[[cont]]) 
}


##############################################
fits <- list()
##############################################
for(cont in countries.full){
    data.full <- terrorism_aggregate
    data <- data.full[data.full$country %in% c(countries.full[cont],nearby.countries[[cont]]),]
    ## Create a named data frame of covariates
    covariates <- data.frame(population = data$pop, time.to.city = data$tt,
                             luminosity = data$lum)
    ## x, y, z locations (latitude and longitude projected onto the unit sphere)
    locs <- cbind(data$x.coord, data$y.coord, data$z.coord)
    ## create temporal indecies
    time <- data$iyear - min(data$iyear) + 1
    ## fit for out of sample predictions
    ## Put NA values at pred locations if PRED is TRUE
    if(PRED){data$total[data$iyear == pred.year & data$country == countries.full[cont]] <- NA}
    fits.tmp <- geo.fit(mesh = meshs[[cont]], locs = locs, response = data$total,covariates = covariates,
                             control.time = list(model = "rw1",
                                                 param = list(theta = list(prior = "pc.prec",param=c(1,0.01)))),
                             temp = time,family = "poisson", sig0 = 0.2, rho0 = 0.01,Prho = 0.9,
                             control.compute = list(waic = TRUE,config = TRUE,openmp.strategy = "huge"),
                             control.inla = list(int.strategy = "eb",strategy = "gaussian",diagonal = 100))
    cat(cont, " initial model fitted","\n")
    fits[[cont]] <- geo.fit(mesh = meshs[[cont]], locs = locs, response = data$total,covariates = covariates,
                            control.time = list(model = "rw1",
                                                param = list(theta = list(prior = "pc.prec",param=c(1,0.01)))),
                            temp = time,family = "poisson", sig0 = 0.2, rho0 = 0.01,Prho = 0.9,
                            control.compute = list(waic = TRUE,openmp.strategy = "huge"),
                            control.mode = list(result = fits.tmp, restart = TRUE),
                            control.inla = list(diagonal = 100))
    cat(cont, "full model fitted","\n")
}


