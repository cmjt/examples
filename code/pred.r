## required libraries
library(lgcpSPDE) ## install from github by running devtools::install_github("cmjt/lgcpSPDE") in R
#############################################
library(foreach) ## for parallelization
library(doParallel)
##### Do you want to fit a "dirty" model or not (run the following line each time you want to change
## change to FALSE if you don't want the quick eb and gaussian inla strategies to be used
quick <- TRUE; if(quick){control.inla <- list(int.strategy = "eb",strategy = "gaussian",diagonal = 100)};if(!quick){ control.inla <- list(diagonal = 100)}
## control coarseness of the projections
dims <- c(3000,3000)
## vector of countries we are interested in
countries <- c("AFG","IRQ","IND","PHL","RUS","LBY","PAK","NGA","IRN","SYR","TUR","YEM","UKR")
# list of spatial polygons of above countries
sps <- sapply(countries, function(x) world[world$sov_a3 == x,])
## create mesh of the world projected onto the unit sphere
bdry <- inla.sp2segment(world)
bdry$loc <- inla.mesh.map(bdry$loc, projection = "longlat",inverse = TRUE)
mesh <- inla.mesh.2d(boundary = bdry, max.edge = c(6,100)/180,cutoff = 6/180) ## plot(mesh) to vizualise

## in-sample predicition (simply just one model worldwide)
##################################################

data <- terrorism_data
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
## fit <- geo.fit(mesh = mesh, locs = locs, response = data$total,
##                covariates = covariates,
##                control.time = list(model = "rw1",
##                                    param = list(theta = list(prior = "pc.prec",
##                                                              param=c(1,0.01)))),
##                temp = time,family = "poisson", sig0 = 0.2, rho0 = 0.01,Prho = 0.9,
##                control.compute = list(waic = TRUE,config = TRUE,openmp.strategy = "huge"), 
##                control.inla = control.inla)

## ## extract "in-sample" fields for the whole world
## fit.fields <- find.fields(fit, mesh = mesh, n.t = length(table(time)),
##                               spatial.polygon = world,dims = dims)
## proj <- inla.mesh.projector(mesh,dims = dims)## set up projection
## ## plot these "in-sample" fields worldwide
## cols <- topo.colors(100) ## colours for plotting
## pdf("insamples_fields_world.pdf", paper='A4r',width = 11,height = 8)
## for(i in 1:length(fit.fields[[1]])){
##     image.plot(proj$x,proj$y,fit.fields[[1]][[i]],
##                axes  = FALSE, xlab = "",ylab = "",col = cols,
##                main = paste("In-sample estimated spatial effect---",names(table(data$iyear))[i],sep = ""))
##     plot(world, add = TRUE)
## }
## dev.off()
    

## out-of-sample prediction
#########################################################
cores <- detectCores() ## to detect the number of cores your computer has
cl <- makeCluster(cores[1]-1) ## so as to not overload your computer use one fewer of its cores
registerDoParallel(cl) ## register that you want to use cores number of cores...
##############################################
## which year we want to predict
pred.year <- 2016
##################################################
pred.fits <- foreach(i = 1:length(countries), .packages = "lgcpSPDE",
                       .errorhandling = "pass") %dopar%
    {
        data <- terrorism_data
        ## Create a named data frame of covariates
        covariates <- data.frame(population = data$pop, time.to.city = data$tt,
                                 luminosity = data$lum)
        ## x, y, z locations (latitude and longitude projected onto the unit sphere)
        locs <- cbind(data$x.coord, data$y.coord, data$z.coord)
        ## create temporal indecies
        time <- data$iyear - min(data$iyear) + 1
        ## fit for out of sample predictions
        ## Put NA values at pred locations
        data$total[data$iyear == pred.year & data$country == countries[i]] <- NA
        pred.fit <- geo.fit(mesh = mesh, locs = locs, response = data$total,covariates = covariates,
                            control.time = list(model = "rw1",
                                                param = list(theta = list(prior = "pc.prec",param=c(1,0.01)))),
                            temp = time,family = "poisson", sig0 = 0.2, rho0 = 0.01,Prho = 0.9,
                            control.compute = list(waic = TRUE,config = TRUE),
                            control.inla = control.inla)
    }

stopCluster(cl) ## stop cluster

## extract list of fields, subset to each country's spatial polygon.

pred.fields <- list()
for(i in 1:length(pred.fits)){
    pred.fields[[i]] <- find.fields(pred.fits[[i]], mesh = mesh, n.t = length(table(time)),
                                    spatial.polygon = sps[[i]],dims = dims)
}


## name lists
names(sps) <- names(pred.fields) <- names(pred.fits) <- countries

### plotting loop for each countries in- and out-of sample prediction
pdf(file = "pnas_predictions.pdf", paper='A4r',width = 11,height = 8)
for(i in names(pred.fields)){
    par(mfrow = c(1,2),mar = c(0,0,2,6))
    tmp.fld <- find.fields(fit, mesh = mesh, n.t = length(table(time)),
                              spatial.polygon = sps[[i]],dims = dims)
    image.plot(proj$x,proj$y,tmp.fld[[1]][[7]],axes  = FALSE, xlab = "",ylab = "",col = cols,
               xlim = sps[[i]]@bbox[1,],ylim = sps[[i]]@bbox[2,])
    title(paste(i,"---",pred.year, "spatial effect for in-sample prediction"),
               cex.main = 0.7)
    plot(sps[[i]], add = TRUE)
    image.plot(proj$x,proj$y,pred.fields[[i]][[1]][[7]],axes  = FALSE, xlab = "",ylab = "",col = cols,
               xlim = sps[[i]]@bbox[1,],ylim = sps[[i]]@bbox[2,])
    title(paste(i,"---",pred.year, "spatial effect for out-of-sample prediction"),
               cex.main = 0.7)
    plot(sps[[i]], add = TRUE)
}
dev.off()
