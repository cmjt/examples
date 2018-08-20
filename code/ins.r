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


## Code to produce fitted values of models scaled between 0 and 1
## extract "in-sample" fields for the whole world
fit.fields <- find.fields(fit, mesh = mesh, n.t = length(table(time)),
                              spatial.polygon = world,dims = dims)
proj <- inla.mesh.projector(mesh,dims = dims) ## set up projection


## get raster values (will take a while)
source("get_raster_vals.r")
## manually calculate fitted values
coefs <- summary(fit)$fixed[,1] ## coefficients of the scaled covariates

## manually calcluate values for each projection pixel and plot
cols <- topo.colors(100) ## colours for plotting
## plot on link scale
pdf("link_world.pdf", paper='A4r',width = 11,height = 8)
for(i in 1:length(fit.fields[[1]])){
    resp <- coefs[1] + coefs[2]*pops[[i]] + coefs[3]*tt + coefs[4]*lums[[i]]
    resp <- matrix(resp,ncol = dims[2],nrow = dims[1]) + fit.fields[[1]][[i]]
    image.plot(proj$x,proj$y,resp,
               axes  = FALSE, xlab = "",ylab = "",
               main = paste("Predicted response on link scale---",names(table(data$iyear))[i],sep = ""))
    plot(world, add = TRUE)
}
dev.off()
## plot on response scale
pdf("response_world.pdf", paper='A4r',width = 11,height = 8)
for(i in 1:length(fit.fields[[1]])){
    resp <- coefs[1] + coefs[2]*pops[[i]] + coefs[3]*tt + coefs[4]*lums[[i]]
    resp <- matrix(resp,ncol = dims[2],nrow = dims[1]) + fit.fields[[1]][[i]]
    image.plot(proj$x,proj$y,exp(resp),
               axes  = FALSE, xlab = "",ylab = "",
               main = paste("Predicted response on response scale---",names(table(data$iyear))[i],sep = ""))
    plot(world, add = TRUE)
}
dev.off()
## plot on scaled (0--1) response scale
pdf("response_world_scaled.pdf", paper='A4r',width = 11,height = 8)
for(i in 1:length(fit.fields[[1]])){
    resp <- coefs[1] + coefs[2]*pops[[i]] + coefs[3]*tt + coefs[4]*lums[[i]]
    resp <- matrix(resp,ncol = dims[2],nrow = dims[1]) + fit.fields[[1]][[i]]
    resp <- exp(resp)
    mx <- max(c(resp),na.rm = TRUE)
    resp.sc <- resp/mx
    image.plot(proj$x,proj$y,resp.sc,
               axes  = FALSE, xlab = "",ylab = "",
               main = paste("Predicted response on response scale---",names(table(data$iyear))[i],sep = ""))
    plot(world, add = TRUE)
}
dev.off()
### plotting loop for each country in- sample prediction
pdf(file = "pnas_in-predictions_fields.pdf", paper='A4r',width = 11,height = 8)
for(i in names(pred.fields)){
    par(mar = c(0,0,2,6))
    tmp.fld <- find.fields(fit, mesh = mesh, n.t = length(table(time)),
                              spatial.polygon = sps[[i]],dims = dims)
    image.plot(proj$x,proj$y,tmp.fld[[1]][[8]],axes  = FALSE, xlab = "",ylab = "",col = cols,
               xlim = sps[[i]]@bbox[1,],ylim = sps[[i]]@bbox[2,])
    title(paste(i,"---",pred.year, "spatial effect for in-sample prediction"),
               cex.main = 0.7)
    plot(sps[[i]], add = TRUE)
}
dev.off()

### plotting loop for each country in- sample prediction on scaled response scale
pdf(file = "pnas_in-predictions_scaled.pdf", paper='A4r',width = 11,height = 8)
for(i in names(pred.fields)){
    par(mar = c(0,0,2,6))
    tmp.fld <- find.fields(fit, mesh = mesh, n.t = length(table(time)),
                           spatial.polygon = sps[[i]],dims = dims)
    resp <- coefs[1] + coefs[2]*pops[[8]] + coefs[3]*tt + coefs[4]*lums[[8]]
    resp <- matrix(resp,ncol = dims[2],nrow = dims[1]) + tmp.fld[[1]][[8]]
    resp <- exp(resp)
    mx <- max(c(resp),na.rm = TRUE)
    resp.sc <- resp/mx
    image.plot(proj$x,proj$y,resp.sc,axes  = FALSE, xlab = "",ylab = "",col = cols,
               xlim = sps[[i]]@bbox[1,],ylim = sps[[i]]@bbox[2,])
    title(paste(i,"---",pred.year, "spatial effect for in-sample prediction"),
               cex.main = 0.7)
    plot(sps[[i]], add = TRUE)
}
dev.off()
