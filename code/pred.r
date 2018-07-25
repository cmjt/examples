## required libraries
library(lgcpSPDE) ## install from github by running devtools::install_github("cmjt/lgcpSPDE") in R
#############################################
library(foreach) ## for parallelization
library(doParallel)
cores <- detectCores() ## to detect the number of cores your computer has
cl <- makeCluster(cores[1]-1) ## so as to not overload your computer use one fewer of its cores
registerDoParallel(cl) ## register that you want to use cores number of cores...
##############################################

## vector of countries we are interested in
countries <- c("AFG","IRQ","IND","PHL","RUS","LBY","PAK","NGA","IRN","SYR","TUR","YEM","UKR")
## list of spatial polygons of above countries
sps <- sapply(countries, function(x) world[world$sov_a3 == x,])
## list of meshes for each country of interest
meshes <- lapply(sps, function(x) inla.mesh.2d(boundary = inla.sp2segment(x),
                                               max.edge=c(0.4,1),cutoff = 0.4)) ## far too course ATM

## in sample predicition
##################################################
fit.fields <- foreach(i = 1:length(countries),.combine = rbind, .packages = "lgcpSPDE",
                       .errorhandling = "pass") %dopar% {
    ## subset data to country of interest
    data <- subset(terrorism_data, terrorism_data$country == countries[i])
    ## Create a named data frame of covariates
    covariates <- data.frame(population = data$pop, time.to.city = data$tt,
                             luminosity = data$lum)
    ## latitude and longitude locations
    locs <- cbind(data$longitude,data$latitude)
    ## create temporal indecies
    time <- data$iyear - min(data$iyear) + 1
    ## fit for in sample predictions
    fit <- geo.fit(mesh = meshes[[i]], locs = locs, response = data$total,
                   covariates = covariates,
                   control.time = list(model = "rw1",
                                       param = list(theta = list(prior = "pc.prec",
                                                                 param=c(1,0.01)))),
                   temp = time,family = "poisson", sig0 = 7, rho0 = 4.7,
                   control.compute = list(waic = TRUE,config = TRUE),
                   control.inla=list(int.strategy="eb"))
    ## extract "in-sample" fields
    tmp.fit.fields <- find.fields(fit, mesh = meshes[[i]], n.t = length(table(time)),
                                  spatial.polygon = sps[[i]],dims = c(1000,1000))
}

## out-of-sample prediction
## which year we want to predict
pred.year <- 2016
##################################################
pred.fields <- foreach(i = 1:length(countries),.combine = rbind, .packages = "lgcpSPDE",
                       .errorhandling = "pass") %dopar% {
    ## subset data to country of interest
    data <- subset(terrorism_data, terrorism_data$country == countries[i])
    ## Create a named data frame of covariates
    covariates <- data.frame(population = data$pop, time.to.city = data$tt,
                             luminosity = data$lum)
    ## latitude and longitude locations
    locs <- cbind(data$longitude,data$latitude)
    ## create temporal indecies
    time <- data$iyear - min(data$iyear) + 1
    ## fit for out of sample predictions
    ## Put NA values at pred locations
    data$total[data$iyear == pred.year] <- NA
    pred.fit <- geo.fit(mesh = meshes[[i]], locs = locs, response = data$total,covariates = covariates,
                        control.time = list(model = "rw1",
                                            param = list(theta = list(prior = "pc.prec",param=c(1,0.01)))),
                        temp = time,family = "poisson", sig0 = 7, rho0 = 4.7,
                        control.compute = list(waic = TRUE,config = TRUE),
                        control.inla=list(int.strategy="eb"))
    ## extract "out-of-sample" fields
    tmp.pred.fields <- find.fields(pred.fit, mesh = meshes[[i]], n.t = length(table(time)),
                                   spatial.polygon = sps[[i]],dims = c(1000,1000))
}

stopCluster(cl) ## stop cluster

## name lists
names(sps) <- names(meshes) <- names(fit.fields) <- names(pred.fields) <- countries
## There is a problem with Russia (the mesh is far too fine)
fit.fields[which(sapply(fit.fields,class) != "list")] <- NULL; pred.fields[which(sapply(pred.fields,class) != "list")]<- NULL
## Niether Libya nor Syria have any observations in all years
fit.fields[which(sapply(fit.fields,length) != 7)] <- NULL; pred.fields[which(sapply(pred.fields,length) != 7)]<- NULL


### plotting loop for each countries in and out-of sample prediction
pdf(file = "pnas_predictions.pdf", paper='A4r',width = 11,height = 8)
for(i in names(fit.fields)){
    cols = heat.colors(100)
    par(mfrow = c(1,2),mar = c(0,0,0,0))
    tmp.dat <- subset(terrorism_data, terrorism_data$iyear == pred.year &
                                      terrorism_data$country == i)[,c(1,2,4)]
    proj <- inla.mesh.projector(meshes[[i]],dims = c(1000,1000))
    image.plot(proj$x,proj$y,fit.fields[[i]][[7]],axes  = FALSE, xlab = "",ylab = "",col = cols)
    legend("top", bty = "n", legend = paste(i,"---",pred.year,
                                            "spatial effect for in-sample prediction"))
    plot(sps[[i]], add = TRUE)
    image.plot(proj$x,proj$y,pred.fields[[i]][[7]],axes  = FALSE, xlab = "",ylab = "",col = cols)
    legend("top", bty = "n", legend = paste(i,"---",pred.year,
                                            "spatial effect for out-of-sample prediction"))
    plot(sps[[i]], add = TRUE)
}
dev.off()
