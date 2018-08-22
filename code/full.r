## required libraries
library(lgcpSPDE) ## install from github by running devtools::install_github("cmjt/lgcpSPDE") in R
data <- terrorism_aggregate
#############################################
## x, y, z locations (latitude and longitude projected onto the unit sphere)
locs <- cbind(data$x.coord, data$y.coord, data$z.coord)
## create mesh of the world projected onto the unit sphere
bdry <- inla.sp2segment(world)
bdry$loc <- inla.mesh.map(bdry$loc, projection = "longlat",inverse = TRUE)
mesh <- inla.mesh.2d(boundary = bdry, loc = locs, max.edge = c(5,50)/180,cutoff = 5/180 ) ## plot(mesh) to vizualise

##############################################
##############################################

## Create a named data frame of covariates
covariates <- data.frame(population = data$pop, time.to.city = data$tt,
                         luminosity = data$lum)

## create temporal indecies
time <- data$iyear - min(data$iyear) + 1


fit.quick <- geo.fit(mesh = mesh, locs = locs, response = data$total,
                     covariates = covariates,
                     control.time = list(model = "rw1",
                                         param = list(theta = list(prior = "pc.prec",
                                                                   param=c(1,0.01)))),
                     temp = time,family = "poisson", sig0 = 0.2, rho0 = 0.01,Prho = 0.9,
                     control.compute = list(waic = TRUE,config = TRUE,openmp.strategy = "huge"), 
                     control.inla = list(int.strategy = "eb",strategy = "gaussian",diagonal = 100)) ## initial fir for decent starting values

fit <- geo.fit(mesh = mesh, locs = locs, response = data$total,
               covariates = covariates,
               control.time = list(model = "rw1",
                                   param = list(theta = list(prior = "pc.prec",
                                                             param=c(1,0.01)))),
               temp = time,family = "poisson", sig0 = 0.2, rho0 = 0.01,Prho = 0.9,
               control.compute = list(waic = TRUE,config = TRUE,openmp.strategy = "huge"),
               control.mode = list(result = fit.quick, restart = TRUE),
               control.inla = list(diagonal = 100)) ## full fit using starting values for quick fit above
