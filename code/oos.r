## required libraries
library(lgcpSPDE) ## install from github by running devtools::install_github("cmjt/lgcpSPDE") in R
library(raster) ## to extract covariates for plotting
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

## out-of-sample prediction
##############################################
## which year we want to predict
pred.year <- 2017
pred.fit.summary <- pred.fields <- list()
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
    pred.fit.tmp <- geo.fit(mesh = mesh, locs = locs, response = data$total,covariates = covariates,
                            control.time = list(model = "rw1",
                                                param = list(theta = list(prior = "pc.prec",param=c(1,0.01)))),
                            temp = time,family = "poisson", sig0 = 0.2, rho0 = 0.01,Prho = 0.9,
                            control.compute = list(waic = TRUE,config = TRUE),
                            control.inla = control.inla)
    pred.fit.summary[[i]] <- summary(pred.fit.tmp)$fixed ## summaries
    ## extract fields for each country in each year
    pred.fields[[i]] <- find.fields(pred.fit.tmp, mesh = mesh,n.t = length(table(time)),
                                    spatial.polygon = sps[[i]],dims = dims)[[1]]
    cat(countries[[i]], "model fitted","\n")
    
}

## name lists
names(sps) <- names(pred.fit.summary) <- names(pred.fields) <- countries.full


proj <- inla.mesh.projector(mesh,dims = dims) ## set up projection


### plotting loop for each country out-of sample prediction
cols <- topo.colors(100) ## colours for plotting
pdf(file = "pnas_out-predictions_fields.pdf", paper='A4r',width = 11,height = 8)
for(i in names(pred.fields)){
    par(mar = c(0,0,2,6))
    image.plot(proj$x,proj$y,pred.fields[[i]][[8]],axes  = FALSE, xlab = "",ylab = "",col = cols,
               xlim = sps[[i]]@bbox[1,],ylim = sps[[i]]@bbox[2,])
    title(paste(i,"---",pred.year, "spatial effect for out-of-sample prediction"),
               cex.main = 0.7)
    plot(sps[[i]], add = TRUE)
}
dev.off()

## get raster values (will take a while)
source("get_raster_vals.r")

### plotting loop for each country out-of sample prediction on scaled response scale
pdf(file = "pnas_out-predictions_scaled.pdf", paper='A4r',width = 11,height = 8)
for(i in names(pred.fields)){
    par(mar = c(0,0,2,6))
    coefs <- pred.fit.summary[[i]][,1]
    resp <- coefs[1] + coefs[2]*pops[[8]] + coefs[3]*tt + coefs[4]*lums[[8]]
    resp <- matrix(resp,ncol = dims[2],nrow = dims[1]) + pred.fields[[i]][[8]]
    resp <- exp(resp)
    mx <- max(c(resp),na.rm = TRUE)
    resp.sc <- resp/mx
    image.plot(proj$x,proj$y,resp.sc,axes  = FALSE, xlab = "",ylab = "",col = cols,
               xlim = sps[[i]]@bbox[1,],ylim = sps[[i]]@bbox[2,])
    title(paste(i,"---",pred.year, "scaled out-of-sample prediction"),
               cex.main = 0.7)
    plot(sps[[i]], add = TRUE)
}
dev.off()
