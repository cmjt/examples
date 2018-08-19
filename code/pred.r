## required libraries
library(lgcpSPDE) ## install from github by running devtools::install_github("cmjt/lgcpSPDE") in R
#############################################
##### Do you want to fit a "dirty" model or not (run the following line each time you want to change
## change to FALSE if you don't want the quick eb and gaussian inla strategies to be used
quick <- TRUE; if(quick){control.inla <- list(int.strategy = "eb",strategy = "gaussian",diagonal = 100)};if(!quick){ control.inla <- list(diagonal = 100)}
## control coarseness of the projections
dims <- c(2000,2000)
## vector of countries we are interested in
countries <- c("AFG","IRQ","IND","PHL","RUS","LBY","PAK","NGA","IRN","SYR","TUR","YEM","UKR")
## full country names for data
countries.full <- c("Afghanistan","Iraq","India"," Philippines","Russia","Libya",
                    "Pakistan","Nigeria","Iran","Syria","Turkey","Yemen","Ukraine")
# list of spatial polygons of above countries
sps <- sapply(countries, function(x) world[world$sov_a3 == x,])
## create mesh of the world projected onto the unit sphere
bdry <- inla.sp2segment(world)
bdry$loc <- inla.mesh.map(bdry$loc, projection = "longlat",inverse = TRUE)
mesh <- inla.mesh.2d(boundary = bdry, max.edge = c(6,100)/180,cutoff = 6/180) ## plot(mesh) to vizualise


## in-sample predicition (simply just one model worldwide)
##################################################

data <- terrorism_aggregate
## Create a named data frame of covariates
covariates <- data.frame(population = data$pop, time.to.city = data$tt,
                         luminosity = data$lum)
## x, y, z locations (latitude and longitude projected onto the unit sphere)
locs <- cbind(data$x.coord, data$y.coord, data$z.coord)
## create temporal indecies
time <- data$iyear - min(data$iyear) + 1

## First This is a quick fit, set control.inla = list(diagonal = 100) if you want a more robust fit,
## or include the argument control.mode = list(result = fit.sv, restart = TRUE) where fit.sv is this first fit
## this gives starting values for the model (see setting quick above)
## fit for in sample predictions
## openmp.strategy = "huge" for inla parralelization
fit <- geo.fit(mesh = mesh, locs = locs, response = data$total,
               covariates = covariates,
               control.time = list(model = "rw1",
                                   param = list(theta = list(prior = "pc.prec",
                                                             param=c(1,0.01)))),
               temp = time,family = "poisson", sig0 = 0.2, rho0 = 0.01,Prho = 0.9,
               control.compute = list(waic = TRUE,config = TRUE,openmp.strategy = "huge"), 
               control.inla = control.inla)

## calculate RMSE and % bias
fitted <- exp(summary(fit)$fixed[1,1] + summary(fit)$fixed[2,1]*covariates[,1] +
                                                    summary(fit)$fixed[3,1]*covariates[,2] +
                                                                summary(fit)$fixed[4,1]*covariates[,3])
observed <- data$total
## boxplot % bias
boxplot((fitted - observed)/observed,ylab = "% bias")
pbias <- sum(fitted - observed, na.rm = TRUE)/sum(observed,na.rm = TRUE) ##  % bias
pbias
mse <- sum((observed - fitted)^2,na.rm = TRUE)/length(fitted) ## MSE
mse
rmse <- sqrt(mse) ## RMSE
rmse
## out-of-sample prediction
##############################################
## which year we want to predict
pred.year <- 2017
pred.fit <- list()
##################################################
for(i in 1:length(countries)){
    data <- terrorism_aggregate
    ## Create a named data frame of covariates
    covariates <- data.frame(population = data$pop, time.to.city = data$tt,
                             luminosity = data$lum)
    ## x, y, z locations (latitude and longitude projected onto the unit sphere)
    locs <- cbind(data$x.coord, data$y.coord, data$z.coord)
    ## create temporal indecies
    time <- data$iyear - min(data$iyear) + 1
    ## fit for out of sample predictions
    ## Put NA values at pred locations
    data$total[data$iyear == pred.year & data$country == countries.full[i]] <- NA
    pred.fit[[i]] <- geo.fit(mesh = mesh, locs = locs, response = data$total,covariates = covariates,
                             control.time = list(model = "rw1",
                                                 param = list(theta = list(prior = "pc.prec",param=c(1,0.01)))),
                             temp = time,family = "poisson", sig0 = 0.2, rho0 = 0.01,Prho = 0.9,
                             control.compute = list(waic = TRUE,config = TRUE),
                             control.inla = control.inla)
}

## name lists
names(sps) <- names(pred.fit) <- countries.full


## Call plotting file which sould be in the working directory
## In addition, note that this assumes that the list of rasters in each folder are in the same order as the split data frame
## (i.e., orderd by year) and give the raw values of the raster images extracted at long&lat locations

source("plot.r") ## will take a while as extracting raster values
