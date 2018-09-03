## calculate RMSE and % bias
fitted <- exp(summary(fit)$fixed[1,1] + summary(fit)$fixed[2,1]*covariates[,1] +
                                                    summary(fit)$fixed[3,1]*covariates[,2] +
                                                                summary(fit)$fixed[4,1]*covariates[,3])
observed <- data$total
## boxplot % bias
## boxplot((fitted - observed)/observed,ylab = "% bias")
pbias <- 100*sum(fitted - observed, na.rm = TRUE)/sum(observed,na.rm = TRUE) ##  % bias
pbias
mse <- sum((fitted - observed)^2,na.rm = TRUE)/length(fitted) ## MSE
mse
rmse <- sqrt(mse) ## RMSE
rmse


## Code to produce fitted values of models scaled between 0 and 1
## extract "in-sample" fields for the whole world
dims <- c(1000,1000)
fit.fields <- find.fields(fit, mesh = mesh, n.t = length(table(time)),
                              spatial.polygon = world,dims = dims)
proj <- inla.mesh.projector(mesh,dims = dims) ## set up projection


## get raster values for the full model (will take a while)
require(raster)
## make sure your working directory has the appropriatly named covariate raster folders
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

## plot on scaled (0--1) link scale
pdf("link_world_scaled.pdf", paper='A4r',width = 11,height = 8)
for(i in 1:length(fit.fields[[1]])){
    resp <- coefs[1] + coefs[2]*pops[[i]] + coefs[3]*tt + coefs[4]*lums[[i]]
    resp <- matrix(resp,ncol = dims[2],nrow = dims[1]) + fit.fields[[1]][[i]]
    mn <-  min(c(resp),na.rm = TRUE)
    resp.sh <- resp + abs(mn)
    mx <-  max(c(resp.sh),na.rm = TRUE)
    resp.sc <- resp.sh/mx 
    image.plot(proj$x,proj$y,resp.sc,
               axes  = FALSE, xlab = "",ylab = "",
               main = paste("Scaled response on link scale---",names(table(data$iyear))[i],sep = ""))
    plot(world, add = TRUE)
}
dev.off()


#############################
#############################
############################
#### out-of-sample plotting
## need AFG, IRN, IRQ, and PAK models & meshes
library(ggplot2)
library(reshape2)
dims <- c(500,500)
## create a list of projections
projs <- lapply(meshs, inla.mesh.projector, dims = dims)
## extract a list of fields
pred.fields <- list()
for(cont in countries.full){
    data.full <- terrorism_aggregate
    data <- data.full[data.full$country %in% c(countries.full[cont],nearby.countries[[cont]]),]
    ## create temporal indecies
    n.t <- length(table(data$iyear - min(data$iyear) + 1))
    pred.fields[[cont]] <- find.fields(pred.fits[[cont]],mesh = meshs[[cont]],n.t = n.t,
                                       spatial.polygon = sps[[cont]], dims = dims)
}

## make sure your working directory has the appropriatly named covariate raster folders
require(raster)
### plotting loop for each country out-of sample prediction
plts.out <- list()
for(cont in names(pred.fields)){
    proj <- projs[[cont]]; source("get_raster_vals.r")
    coefs <- summary(pred.fits[[cont]])$fixed[,1] 
    resp <- coefs[1] + coefs[2]*pops[[8]] + coefs[3]*tt + coefs[4]*lums[[8]]
    resp <- matrix(resp,ncol = dims[2],nrow = dims[1]) + pred.fields[[cont]][[1]][[8]]
    mn <-  min(c(resp),na.rm = TRUE)
    resp.sh <- resp + abs(mn)
    mx <-  max(c(resp.sh),na.rm = TRUE)
    resp.sc <- resp.sh/mx 
    rs <- apply(which(!is.na(resp.sc),arr.ind = TRUE),2,range)
    z <- resp.sc[rs[1,1]:rs[2,1],rs[1,2]:rs[2,2]]
    dimnames(z) <- list(x = proj$x[rs[1,1]:rs[2,1]],y = proj$y[rs[1,2]:rs[2,2]])
    plot2d_1 <- melt(z,value.name="z")
    gg0 <- ggplot(plot2d_1, aes(x,y,z = z,fill = z))
    plts.out[[cont]] <- gg0 + geom_raster(interpolate=TRUE) +
        geom_polygon(data = fortify(sps[[cont]]),
                     aes(x = long, y = lat, z = NULL, fill = NULL ), alpha = 0, color = "black") +
        theme_bw() +
        theme(axis.title.x = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.title.y = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              panel.border = element_blank(), panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()) +
        scale_fill_continuous(low="blue", high="red",guide = "colorbar",na.value="white",
                              name = "Scaled response")
}


### for the in-sample prediction, plotted for each country only (need fit object)

data <- terrorism_aggregate
#############################################
## x, y, z locations (latitude and longitude projected onto the unit sphere)
locs <- cbind(data$x.coord, data$y.coord, data$z.coord)
## create mesh of the world projected onto the unit sphere
bdry <- inla.sp2segment(world)
bdry$loc <- inla.mesh.map(bdry$loc, projection = "longlat",inverse = TRUE)
mesh <- inla.mesh.2d(boundary = bdry, loc = locs, max.edge = c(5,50)/180,cutoff = 5/180 ) ## plot(mesh) to vizualise
proj <- inla.mesh.projector(mesh,dims = dims) ## set up projection
source("get_raster_vals.r")

## manually calculate fitted values
coefs <- summary(fit)$fixed[,1] ## coefficients of the scaled covariates

countries.full <- c("Afghanistan","Iraq","Pakistan","Iran")
plts.in <- list()
for(cont in countries.full){
    fields <- find.fields(fit, mesh = mesh, n.t = length(table(data$iyear - min(data$iyear) + 1)),
                              spatial.polygon = sps[[cont]],dims = dims)[[1]][[8]]

    resp <- coefs[1] + coefs[2]*pops[[8]] + coefs[3]*tt + coefs[4]*lums[[8]]
    resp <- matrix(resp,ncol = dims[2],nrow = dims[1]) + fields
    mn <-  min(c(resp),na.rm = TRUE)
    resp.sh <- resp + abs(mn)
    mx <-  max(c(resp.sh),na.rm = TRUE)
    resp.sc <- resp.sh/mx 
    rs <- apply(which(!is.na(resp.sc),arr.ind = TRUE),2,range)
    z <- resp.sc[rs[1,1]:rs[2,1],rs[1,2]:rs[2,2]]
    dimnames(z) <- list(x = proj$x[rs[1,1]:rs[2,1]],y = proj$y[rs[1,2]:rs[2,2]])
    plot2d_1 <- melt(z,value.name="z")

    gg0 <- ggplot(plot2d_1, aes(x,y,z = z,fill = z))
    plts.in[[cont]] <- gg0 + geom_raster(interpolate=TRUE) +
        geom_polygon(data = fortify(sps[[cont]]),
                     aes(x = long, y = lat, z = NULL, fill = NULL ), alpha = 0, color = "black") +
        theme_bw() +
        theme(axis.title.x = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.title.y = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              panel.border = element_blank(), panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()) +
        scale_fill_continuous(low="blue", high="red",guide = "colorbar",na.value="white",
                              name = "Scaled response") + ggtitle(cont)
        
}
require(gridExtra)
pdf("~/Desktop/COVARIATE/comparisons.pdf")
grid.arrange(arrangeGrob(plts.in[[1]],plts.in[[2]],plts.in[[3]],plts.in[[4]],top = "In-sample", ncol = 1),
             arrangeGrob(plts.out[[1]],plts.out[[2]], plts.out[[3]],plts.out[[4]],top = "Out-of-sample", ncol = 1),
             ncol = 2)
dev.off()
