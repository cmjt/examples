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
