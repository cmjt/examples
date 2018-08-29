## required libraries
library(lgcpSPDE) ## install from github by running devtools::install_github("cmjt/lgcpSPDE") in R
library(raster) ## to extract covariates for plotting
library(rgeos)
library(rgdal) ## for nearby countries
#############################################
##### Do you want to fit a "dirty" model or not (run the following line each time you want to change
## change to FALSE if you don't want the quick eb and gaussian inla strategies to be used
quick <- TRUE; if(quick){control.inla <- list(int.strategy = "eb",strategy = "gaussian",diagonal = 100)};if(!quick){ control.inla <- list(diagonal = 100)}
## control coarseness of the projections
dims <- c(2000,2000)
## full country names for data we are interested in
countries.full <- c("Afghanistan","Iraq","India","Philippines","Russia","Libya",
                    "Pakistan","Nigeria","Iran","Syria","Turkey","Yemen","Ukraine")
# list of spatial polygons of above countries
sps <- sapply(countries.full, function(x) world[world$name == x,])
## prediction year
pred.year <- 2017


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

mn <- list(2/180,1/180,2/180,0.5/180,2/180,1.5/180,2/180,1/180,1/180,0.5/180,1/180,1/180,0.5/180)
max <- list(10/180,5/180,10/180,3/180,10/180,5/180,10/180,3/180,3/180,2/180,5/180,5/180,3/180)
names(mn) <- names(max) <- countries.full ## names mesh resolution lists


## create mesh for each "set" of countries projected onto the unit sphere
for(cont in countries.full){
    bdry <- inla.sp2segment(sp.near[[cont]])
    bdry$loc <- inla.mesh.map(bdry$loc, projection = "longlat",inverse = TRUE)
    meshs[[cont]] <- inla.mesh.2d(boundary = bdry, max.edge = c(mn[[cont]],max[[cont]]),cutoff = mn[[cont]]) 
}


##############################################
pred.fits <- list()
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
    ## Put NA values at pred locations
    data$total[data$iyear == pred.year & data$country == countries.full[cont]] <- NA
    pred.fitx[[cont]] <- geo.fit(mesh = meshs[[cont]], locs = locs, response = data$total,covariates = covariates,
                            control.time = list(model = "rw1",
                                                param = list(theta = list(prior = "pc.prec",param=c(1,0.01)))),
                            temp = time,family = "poisson", sig0 = 0.2, rho0 = 0.01,Prho = 0.9,
                            control.compute = list(waic = TRUE,config = TRUE),
                            control.inla = control.inla)
    cat(cont, "model fitted","\n")
}


### plotting loop for each country out-of sample prediction
## cols <- topo.colors(100) ## colours for plotting
## pdf(file = "pnas_out-predictions_fields.pdf", paper='A4r',width = 11,height = 8)
## for(i in names(pred.fields)){
##     par(mar = c(0,0,2,6))
##     image.plot(projs[[i]]$x,projs[[i]]$y,pred.fields[[i]][[8]],axes  = FALSE, xlab = "",ylab = "",col = cols,
##                xlim = sps[[i]]@bbox[1,],ylim = sps[[i]]@bbox[2,])
##     title(paste(i,"---",pred.year, "spatial effect for out-of-sample prediction"),
##                cex.main = 0.7)
##     plot(sps[[i]], add = TRUE)
## }
## dev.off()



## ### plotting loop for each country out-of sample prediction on scaled response scale
## pdf(file = "pnas_out-predictions_scaled.pdf", paper='A4r',width = 11,height = 8)
## for(i in names(pred.fields)){
##     par(mar = c(0,0,2,6))
##     image.plot(proj$x,proj$y,resp.sc[[i]],axes  = FALSE, xlab = "",ylab = "",col = cols,
##                xlim = sps[[i]]@bbox[1,],ylim = sps[[i]]@bbox[2,])
##     title(paste(i,"---",pred.year, "scaled out-of-sample prediction"),
##                cex.main = 0.7)
##     plot(sps[[i]], add = TRUE)
## }
## dev.off()
