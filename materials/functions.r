#A Bayesian Approach to Modelling Subnational Spatial Dynamics of Worldwide                   #
#Non-State Terrorism, 2010 - 2014                                                             #
#Year: 2017                                                                                   #
#Institution: University of St Andrews, UK                                                    #
#Journal: JRSS Series A                                                                       #
###############################################################################################
###############################################################################################
#functions: functions to fit spatial and space-time model with covariates
###############################################################################################
###############################################################################################
#' Function to fit a  either a spatial or spatio-temporal model to geo-statistical data with an intercept and covariates
#'
#' @return A \code{inla} result object
#'
#' @param mesh a ``mesh'' object i.e. delauney triangulation of the domain, an
#' object returned by \link{make.mesh}.
#' @param locs a matrix of observation locations, where each row corresponds to the observation. 
#' @param response a vector of response variable, each corresponds to the spatial locations
#' in \code{locs}.
#' @param temp (optional) a numeric vector specifying a temporal index for each observation (starting at 1.....T).
#' @param covariates (optional) a named data.frame of covariates.
#' @param family a character vector specifying the assumed likelihood of the response, by default is "gaussian".
#' @param control.time (optional) supplied if the \code{temp} argumet is given to fit a spatio-temporal model. This argument
#' controls the model and prior put on the hyperparameters of the model for the temporal component of the spatio-temporal
#' model. By default this is \code{list(model = 'ar1', param = list(theta = list(prior='pccor1', param = c(0, 0.9))))}
#' which is a pc.prior put on the rho coefficient of a AR(1) model with P(rho>0)=0.9.
#' Refer to Simpson, martins, and rue for further details *****put in proper refs*****
#' @param control.inla a list which controls the fitting procedures INLA uses see Rue et al. ***ref book***
#' by default this is \code{list(strategy='gaussian',int.strategy = 'eb')} for quick and dirty fitting.
#' @param control.compute a list of fit statistics the user wants INLA to return. By default this
#' is \code{list(dic = TRUE, waic = TRUE,cpo = TRUE, config = TRUE)}.
#' @param non.linear (optional) should be used if the user requires a non-linear covariate to be included in the model
#' Must be supplied as a named list with elements \code{random.effect} a numeric vector of the random effect indecies,
#' and \code{model} the random effect model the user wishes to use for \code{random.effect}
#' @param prediction (optional) should be used if the uses wnts to run a prediction model. Must be supplied as a
#' named list with \code{pred.locs} the locations where predictionas are to be made, and only if a spatio-tempral
#' model is to be fitted \code{pred.temp} the temporal indecies for the predictions (the same length as \code{pred.locs}).
#' @param sig0 by default = 1, typical standard deviation to use pc priors for hyperparams of spde model
#' @param Psig by default = 0.5 prob for sigma of pc prior
#' @param rho0 by default = 0.3, typical range to use pc priors for hyperparams of spde model
#' @param Prho by default = 0.5 prob for rho of pc prior
#' @param verbose Logical if \code{TRUE} model fit is output to screen.
#' @param ... add inla options to speed up computation i.e., by giving starting values from a previos model
#'

geo.fit <- function(mesh = NULL,  locs = NULL, response = NULL, temp = NULL,covariates = NULL, family = "gaussian",
                    control.time = list(model = 'ar1', param = list(theta = list(prior='pccor1', param = c(0, 0.9)))),
                    control.inla = list(strategy='gaussian',int.strategy = 'eb'),
                    control.compute = list(dic = TRUE, waic = TRUE,cpo = TRUE, config = TRUE),
                    non.linear = NULL, prediction = NULL, sig0 = 1,Psig = 0.5, rho0 = 0.3,Prho = 0.5,
                    verbose = FALSE, ...){
    if(is.null(temp)&is.null(covariates)){
        fit <- geo.spatial.fit(mesh = mesh, locs = locs, response = response,covariates = NULL,
                               family = family, control.inla = control.inla, control.compute = control.compute,
                               non.linear = non.linear, prediction = prediction,sig0 = sig0, Psig = Psig,
                               rho0 = rho0, Prho = Prho,
                               verbose = verbose, ...)
    }
    if(is.null(temp)&!is.null(covariates)){
        fit <- geo.spatial.fit(mesh = mesh, locs = locs, response = response,covariates = covariates,
                               family = family, control.inla = control.inla, control.compute = control.compute,
                               non.linear = non.linear, prediction = prediction, sig0 = sig0, Psig = Psig,
                               rho0 = rho0, Prho = Prho,
                               verbose = verbose, ...)
    }
    if(!is.null(temp)&is.null(covariates)){
        fit <- geo.spatial.temporal.fit(mesh = mesh, locs = locs, response = response, temp = temp,
                                        family = family, covariates = NULL,
                                        control.time = control.time, control.inla = control.inla,
                                        control.compute = control.compute,
                                        non.linear = non.linear, prediction = prediction,sig0 = sig0, Psig = Psig,
                                        rho0 = rho0, Prho = Prho, 
                                        verbose = verbose, ...)
    }
    if(!is.null(temp)&!is.null(covariates)){
        fit <- geo.spatial.temporal.fit(mesh = mesh, locs = locs, response = response,covariates = covariates,
                                        temp = temp, family = family, control.time = control.time,
                                        control.inla = control.inla, control.compute = control.compute,
                                        non.linear = non.linear, prediction = prediction, sig0 = sig0,  Psig = Psig,
                                        rho0 = rho0, Prho = Prho,
                                        verbose = verbose, ...)
    }
    return(fit)
}
    




#' spatial only model fitting
#'
geo.spatial.fit <- function(mesh, locs, response, covariates, family, control.inla, control.compute,
                            non.linear, prediction, sig0, Psig, rho0, Prho, verbose, ... ){
    spde <-inla.spde2.matern.new(mesh, prior.pc.rho = c(rho0, Prho), prior.pc.sig = c(sig0, Psig))
    nv <- mesh$n
    n <- nrow(locs)
    if(is.null(covariates)&is.null(non.linear)&is.null(prediction)){
        index <- inla.spde.make.A(mesh = mesh, loc = locs)
        stack <- inla.stack(data=list(y = response),
                            A=list(index, 1), tag='observation',
                            effects=list(field = 1:mesh$n, b0 = rep(1,n)))
        formula <- y ~ 0 + b0 +  f(field,model = spde)
    }
    if(is.null(covariates)&is.null(non.linear)&!is.null(prediction)){
        locs.o <- locs
        locs.p <- prediction[["pred.locs"]]
        index.o <- inla.spde.make.A(mesh = mesh, loc = locs.o)
        index.p <- inla.spde.make.A(mesh = mesh, loc = locs.p)
        stk.obvs <- inla.stack(data=list(y = response),
                               A=list(index.o,1), tag = 'observation',
                               effects=list(field = 1:mesh$n, b0 = rep(1,n)))
        stk.prd <- inla.stack(data=list(y = NA), A = list(index.p),
                              effects=list(field = 1:spde$n.spde), tag='prediction')
        stack<-inla.stack(stk.obvs,stk.prd)
        formula <- y ~ 0 + b0 +  f(field,model = spde)
    }
    if(!is.null(covariates)&is.null(non.linear)&is.null(prediction)){
        index <- inla.spde.make.A(mesh = mesh, loc = locs)
        m <- make.covs(covariates)
        cov.effects <- m[[1]]
        cov.form <- m[[2]]
        stack <- inla.stack(data=list(y = response),
                            A=list(index, 1, 1), tag='observation',
                            effects=list(field = 1:mesh$n, b0 = rep(1,n), cov.effects = cov.effects))
        formula <- paste("y", "~  0 + b0 +", cov.form," + f(field, model=spde)")
    }
    if(!is.null(covariates)&is.null(non.linear)&!is.null(prediction)){
        locs.o <- locs
        locs.p <- prediction[["pred.locs"]]
        index.o <- inla.spde.make.A(mesh = mesh, loc = locs.o)
        index.p <- inla.spde.make.A(mesh = mesh, loc = locs.p)
        m <- make.covs(covariates)
        cov.effects <- m[[1]]
        cov.form <- m[[2]]
        stk.obvs <- inla.stack(data=list(y = response),
                            A=list(index.o, 1, 1), tag='observation',
                            effects=list(field = 1:mesh$n, b0 = rep(1,n), cov.effects = cov.effects))
        stk.prd <- inla.stack(data=list(y = NA), A = list(index.p),
                              effects=list(field = 1:spde$n.spde), tag='prediction')
        stack<-inla.stack(stk.obvs,stk.prd)
        formula <- paste("y", "~  0 + b0 +", cov.form," + f(field, model=spde)")
    }
    if(is.null(covariates)&!is.null(non.linear)&is.null(prediction)){
        u <- non.linear[["random.effect"]]
        u.mod <- non.linear[["model"]]
        index <- inla.spde.make.A(mesh = mesh, loc = locs)
        stack <- inla.stack(data=list(y = response),
                            A=list(index, 1, 1), tag='observation',
                            effects=list(field = 1:mesh$n, b0 = rep(1,n),random = u))
        formula <- y ~ 0 + b0 +  f(random, model = u.mod) + f(field,model = spde)
    }
    if(!is.null(covariates)&!is.null(non.linear)&is.null(prediction)){
        u <- non.linear[["random.effect"]]
        u.mod <- non.linear[["model"]]
        index <- inla.spde.make.A(mesh = mesh, loc = locs)
        m <- make.covs(covariates)
        cov.effects <- m[[1]]
        cov.form <- m[[2]]
        stack <- inla.stack(data=list(y = response),
                            A=list(index, 1, 1, 1), tag='observation',
                            effects=list(field = 1:mesh$n, b0 = rep(1,n),random = u, cov.effects = cov.effects))
        formula <- paste("y", "~  0 + b0 +  f(random, model = u.mod) +", cov.form," + f(field, model=spde)")
    }
    if(is.null(covariates)&!is.null(non.linear)&!is.null(prediction)){
        u <- non.linear[["random.effect"]]
        u.mod <- non.linear[["model"]]
        locs.o <- locs
        locs.p <- prediction[["pred.locs"]]
        index.o <- inla.spde.make.A(mesh = mesh, loc = locs.o)
        index.p <- inla.spde.make.A(mesh = mesh, loc = locs.p)
        stk.obvs <- inla.stack(data=list(y = response),
                               A=list(index.o,1,1), tag = 'observation',
                               effects=list(field = 1:mesh$n, b0 = rep(1,n), random = u))
        stk.prd <- inla.stack(data=list(y = NA), A = list(index.p),
                              effects=list(field = 1:spde$n.spde), tag='prediction')
        stack<-inla.stack(stk.obvs,stk.prd)
        formula <- y ~ 0 + b0  +  f(random, model = u.mod) +  f(field,model = spde)
    }
    if(!is.null(covariates)&!is.null(non.linear)&!is.null(prediction)){
        u <- non.linear[["random.effect"]]
        u.mod <- non.linear[["model"]]
        locs.o <- locs
        locs.p <- prediction[["pred.locs"]]
        index.o <- inla.spde.make.A(mesh = mesh, loc = locs.o)
        index.p <- inla.spde.make.A(mesh = mesh, loc = locs.p)
        m <- make.covs(covariates)
        cov.effects <- m[[1]]
        cov.form <- m[[2]]
        stk.obvs <- inla.stack(data=list(y = response),
                            A=list(index.o, 1, 1,1), tag='observation',
                            effects=list(field = 1:mesh$n, b0 = rep(1,n), cov.effects = cov.effects, random = u))
        stk.prd <- inla.stack(data=list(y = NA), A = list(index.p),
                              effects=list(field = 1:spde$n.spde), tag='prediction')
        stack<-inla.stack(stk.obvs,stk.prd)
        formula <- paste("y", "~  0 + b0  +  f(random, model = u.mod) +", cov.form," + f(field, model=spde)")
    }
    result <- inla(as.formula(formula), family = family,
                   data = inla.stack.data(stack),
                   control.predictor=list(A=inla.stack.A(stack)),
                   control.inla = control.inla,
                   control.compute = control.compute,
                   verbose = verbose,
                   ...)
    return(result)
}






#' spatio temporal model fitting
#'
geo.spatial.temporal.fit <- function(mesh, locs, response, covariates, temp,  family,
                                    control.time, control.inla, control.compute,
                                    non.linear, prediction,sig0,Psig, rho0, Prho, verbose, ... ){
    nv <- mesh$n
    k <- (mesh.t <- inla.mesh.1d(temp))$n
    n <- nrow(locs)
    spde <-inla.spde2.matern.new(mesh, prior.pc.rho = c(rho0, Prho), prior.pc.sig = c(sig0, Psig))
    k <- (mesh.t <- inla.mesh.1d(temp))$n
    field <- inla.spde.make.index('field',n.spde = spde$n.spde, group = temp, n.group = k)
    Ast <- inla.spde.make.A(mesh = mesh, loc = locs ,group = temp, n.group = k)
    Y <- response
    if(is.null(covariates)&is.null(non.linear)&is.null(prediction)){
        stack <- inla.stack(data=list(y = Y),
                            A=list(Ast, 1), tag='observation',
                            effects=list(field = field, b0 = rep(1,n)))
        formula <- y ~ 0 + b0 +  f(field,model = spde, group = field.group, control.group = control.time)
    }
    if(is.null(covariates)&is.null(non.linear)&!is.null(prediction)){
        locs.p <- prediction[["pred.locs"]]
        temp.p <- prediction[["pred.temp"]]
        k.p <- (mesh.t <- inla.mesh.1d(temp.p))$n
        Ast.p <- inla.spde.make.A(mesh = mesh, loc = locs.p ,group = temp.p, n.group = k.p)
        field.p <- inla.spde.make.index('field.p',n.spde = spde$n.spde, group = temp.p, n.group = k.p)
        stk.obvs <- inla.stack(data=list(y = Y),
                               A=list(Ast,1), tag = 'observation',
                               effects=list(field = field, b0 = rep(1,n)))
        stk.prd <- inla.stack(data=list(y = NA), A = list(Ast.p),
                              effects=list(field = field.p), tag='prediction')
        stack<-inla.stack(stk.obvs,stk.prd)
        formula <- y ~ 0 + b0 +  f(field,model = spde, group = field.group, control.group = control.time)
    }
    if(!is.null(covariates)&is.null(non.linear)&is.null(prediction)){
        m <- make.covs(covariates)
        cov.effects <- m[[1]]
        cov.form <- m[[2]]
        stack <- inla.stack(data=list(y = Y),
                            A=list(Ast, 1, 1), tag='observation',
                            effects=list(field = field, b0 = rep(1,n), cov.effects = cov.effects))
        formula <- paste("y", "~  0 + b0 +", cov.form,
                         " + f(field,model = spde, group = field.group, control.group = control.time)")
    }
    if(!is.null(covariates)&is.null(non.linear)&!is.null(prediction)){
        locs.p <- prediction[["pred.locs"]]
        temp.p <- prediction[["pred.temp"]]
        k.p <- (mesh.t <- inla.mesh.1d(temp.p))$n
        Ast.p <- inla.spde.make.A(mesh = mesh, loc = locs.p ,group = temp.p, n.group = k.p)
        field.p <- inla.spde.make.index('field.p',n.spde = spde$n.spde, group = temp.p, n.group = k.p)
        m <- make.covs(covariates)
        cov.effects <- m[[1]]
        cov.form <- m[[2]]
        stk.obvs <- inla.stack(data=list(y = Y),
                            A=list(Ast, 1, 1), tag='observation',
                            effects=list(field = field, b0 = rep(1,n), cov.effects = cov.effects))
        stk.prd <- inla.stack(data=list(y = NA), A = list(Ast.p),
                              effects=list(field = field.p), tag='prediction')
        stack<-inla.stack(stk.obvs,stk.prd)
        formula <- paste("y", "~  0 + b0 +", cov.form," + f(field,model = spde, group = field.group, control.group = control.time)")
    }
    if(is.null(covariates)&!is.null(non.linear)&is.null(prediction)){
        u <- non.linear[["random.effect"]]
        u.mod <- non.linear[["model"]]
        stack <- inla.stack(data=list(y = Y),
                            A=list(Ast, 1, 1), tag='observation',
                            effects=list(field = field, b0 = rep(1,n),random = u))
        formula <- y ~ 0 + b0 +  f(random, model = u.mod) + f(field,model = spde,
                                                              group = field.group, control.group = control.time)
    }
    if(!is.null(covariates)&!is.null(non.linear)&is.null(prediction)){
        u <- non.linear[["random.effect"]]
        u.mod <- non.linear[["model"]]
        m <- make.covs(covariates)
        cov.effects <- m[[1]]
        cov.form <- m[[2]]
        stack <- inla.stack(data=list(y = Y),
                            A=list(Ast, 1, 1, 1), tag='observation',
                            effects=list(field = field, b0 = rep(1,n),random = u, cov.effects = cov.effects))
        formula <- paste("y", "~  0 + b0 +  f(random, model = u.mod) +", cov.form,
                         " + f(field,model = spde, group = field.group, control.group = control.time)")
    }
    if(is.null(covariates)&!is.null(non.linear)&!is.null(prediction)){
        locs.p <- prediction[["pred.locs"]]
        temp.p <- prediction[["pred.temp"]]
        k.p <- (mesh.t <- inla.mesh.1d(temp.p))$n
        Ast.p <- inla.spde.make.A(mesh = mesh, loc = locs.p ,group = temp.p, n.group = k.p)
        field.p <- inla.spde.make.index('field.p',n.spde = spde$n.spde, group = temp.p, n.group = k.p)
        u <- non.linear[["random.effect"]]
        u.mod <- non.linear[["model"]]
        stk.obvs <- inla.stack(data=list(y = Y),
                               A=list(Ast,1,1), tag = 'observation',
                               effects=list(field = field, b0 = rep(1,n), random = u))
        stk.prd <- inla.stack(data=list(y = NA), A = list(Ast.p),
                              effects=list(field = field.p), tag='prediction')
        stack<-inla.stack(stk.obvs,stk.prd)
        formula <- y ~ 0 + b0  +  f(random, model = u.mod) +
            f(field,model = spde, group = field.group, control.group = control.time)
    }
    if(!is.null(covariates)&!is.null(non.linear)&!is.null(prediction)){
        locs.p <- prediction[["pred.locs"]]
        temp.p <- prediction[["pred.temp"]]
        k.p <- (mesh.t <- inla.mesh.1d(temp.p))$n
        Ast.p <- inla.spde.make.A(mesh = mesh, loc = locs.p ,group = temp.p, n.group = k.p)
        field.p <- inla.spde.make.index('field.p',n.spde = spde$n.spde, group = temp.p, n.group = k.p)
        u <- non.linear[["random.effect"]]
        u.mod <- non.linear[["model"]]
        m <- make.covs(covariates)
        cov.effects <- m[[1]]
        cov.form <- m[[2]]
        stk.obvs <- inla.stack(data=list(y = Y),
                            A=list(Ast, 1, 1,1), tag='observation',
                            effects=list(field = field, b0 = rep(1,n), cov.effects = cov.effects, random = u))
        stk.prd <- inla.stack(data=list(y = NA), A = list(Ast.p),
                              effects=list(field = field.p), tag='prediction')
        stack<-inla.stack(stk.obvs,stk.prd)
        formula <- paste("y", "~  0 + b0  +  f(random, model = u.mod) +",
                         cov.form," + f(field,model = spde, group = field.group, control.group = control.time)")
    }
    result <- inla(as.formula(formula), family = family,
                   data = inla.stack.data(stack),
                   control.predictor=list(A=inla.stack.A(stack)),
                   control.inla = control.inla,
                   control.compute = control.compute,
                   verbose = verbose,
                   ...)
    return(result)
}






make.covs <- function(covariates){
    n.covs <- ncol(covariates)
    for(i in 1:n.covs){
        assign(colnames(covariates)[i],covariates[,i],envir = .GlobalEnv)
    }
    cov.effects <- sapply(colnames(covariates),get,simplify = FALSE)
    cov.form <- paste(colnames(covariates),collapse = " + ")
    return(list(cov.effects,cov.form))
}

########################## Joint
#### Haakon's code to use pc priors for parameters of latent field

inla.spde2.matern.new = function(mesh, alpha=2, prior.pc.rho, prior.pc.sig){
    d = INLA:::inla.ifelse(inherits(mesh, "inla.mesh"), 2, 1)
    nu = alpha-d/2
    kappa0 = log(8*nu)/2
    tau0   = 0.5*(lgamma(nu)-lgamma(nu+d/2)-d/2*log(4*pi))-nu*kappa0
    spde   = inla.spde2.matern(mesh = mesh,
                               B.tau   = cbind(tau0,   nu,  -1),
                               B.kappa = cbind(kappa0, -1, 0))
    param = c(prior.pc.rho, prior.pc.sig)
    spde$f$hyper.default$theta1$prior = "pcspdega"
    spde$f$hyper.default$theta1$param = param
    spde$f$hyper.default$theta1$initial = log(prior.pc.rho[1])+1
    spde$f$hyper.default$theta2$initial = log(prior.pc.sig[1])-1
    
  # End and return
    return(invisible(spde))  
}

#' Function to fit a  either a spatial or spatio-temporal joint model to geo-statistical data with an intercept and covariates for each likelihood
#'
#' @return A \code{inla} result object
#'
#' @param mesh a ``mesh'' object i.e. delauney triangulation of the domain, an
#' object returned by \link{make.mesh}.
#' @param locs a list of matrcies. The first element holds observation locations for the first likelihood, where each row corresponds to an observation. The second elemenr holds the observation locations for the second likelihood, each row corresponds to an observation. If no second element is supplied the observation locations for the first likelihood are used.
#' @param response a list (length two) of vectors of each response variable, each corresponds to the respective spatial locations
#' in \code{locs}.
#' @param temp (optional) a list of  numeric vectors specifying the temporal indcies for each response respectively.
#' @param covariates (optional) a list (length 2) each element should contain a named data.frame of covariates. The first corresponding to the first likelihood, the second corresponding the the second likelihood.
#' @param family a character vector of length two specifying the assumed likelihood of each response, by default is c("gaussian","gaussian").
#' @param control.time (optional) supplied if the \code{temp} argumet is given to fit a spatio-temporal model. This argument
#' controls the model and prior put on the hyperparameters of the model for the temporal component of the spatio-temporal
#' model. By default this is \code{list(model = 'ar1', param = list(theta = list(prior='pccor1', param = c(0, 0.9))))}
#' which is a pc.prior put on the rho coefficient of a AR(1) model with P(rho>0)=0.9. Assumed to be shared accross both responses
#' Refer to Simpson, martins, and rue for further details *****put in proper refs*****
#' @param control.inla a list which controls the fitting procedures INLA uses see Rue et al. ***ref book***
#' by default this is \code{list(strategy='gaussian',int.strategy = 'eb')} for quick and dirty fitting.
#' @param hyper prior for the copy parameter by default is a N(0,10) i.e.,  list(theta=list(prior='normal', param=c(0,10)))
#' @param control.compute a list of fit statistics the user wants INLA to return. By default this
#' is \code{list(dic = TRUE, waic = TRUE,cpo = TRUE, config = TRUE)}.
#' @param non.linear (optional) a list of named lists should be used if the user requires a non-linear covariate to be included for each likelihood. (i.e., non.linear = list(list(random.effect = idx.1, model = "iid"),list(random.effect = idx.2, model = "iid")) if the user wnats a iid effect for some idx.1 for the first likelihood and another for idx.2 for the second)
#' Must be supplied as a named list with elements \code{random.effect} a numeric vector of the random effect indecies,
#' and \code{model} the random effect model the user wishes to use for \code{random.effect}
#' @param sig0 by default = 1, typical standard deviation to use pc priors for hyperparams of spde model
#' @param Psig by default = 0.5 prob for sigma of pc prior
#' @param rho0 by default = 0.3, typical range to use pc priors for hyperparams of spde model
#' @param Prho by default = 0.5 prob for rho of pc prior
#' @param verbose Logical if \code{TRUE} model fit is output to screen.
#' @param ... add inla options to speed up computation i.e., by giving starting values from a previos model
#' @export

geo.joint.fit <- function(mesh = NULL,  locs = NULL, response = NULL, temp = NULL,covariates = NULL,
                          family = c("gaussian","gaussian"),
                          control.time = list(model = 'ar1', param = list(theta = list(prior='pccor1', param = c(0, 0.9)))),
                          control.inla = list(strategy='gaussian',int.strategy = 'eb'),
                          hyper = list(theta=list(prior='normal', param=c(0,10))),
                          control.compute = list(dic = TRUE, waic = TRUE,cpo = TRUE, config = TRUE),
                          non.linear = NULL, sig0 = 1,Psig = 0.5, rho0 = 0.3,Prho = 0.5,
                          verbose = FALSE, ...){
    if(is.null(temp)){
        fit <- geo.spatial.j.fit(mesh = mesh, locs = locs, response = response,covariates = covariates,
                                 family = family, control.inla = control.inla, control.compute = control.compute,
                                 hyper = hyper, sig0 = sig0, Psig = Psig, rho0 = rho0, Prho = Prho,
                                 verbose = verbose, ...)
    }
    if(!is.null(temp)&is.null(non.linear)){
        fit <- geo.spatial.j.temporal.fit(mesh = mesh, locs = locs, response = response, temp = temp,
                                          family = family, covariates = covariates,
                                          control.time = control.time, control.inla = control.inla,
                                          hyper = hyper, control.compute = control.compute,
                                          sig0 = sig0, Psig = Psig, rho0 = rho0, Prho = Prho,
                                          verbose = verbose, ...)
    }
    if(!is.null(temp)&!is.null(non.linear)){
        fit <- geo.spatial.j.nl.temporal.fit(mesh = mesh, locs = locs, response = response,covariates = covariates,
                                             temp = temp, family = family, control.time = control.time,
                                             sig0 = sig0, Psig = Psig, rho0 = rho0, Prho = Prho,
                                             control.inla = control.inla, hyper = hyper, control.compute = control.compute,
                                             non.linear = non.linear,  verbose = verbose, ...)
    }
    return(fit)
}


#' spatial only fitting
#' 
geo.spatial.j.fit <- function(mesh, locs, response, covariates, family,  control.inla,
                              hyper, control.compute, sig0, Psig, rho0, Prho, verbose, ...){
    spde <-inla.spde2.matern.new(mesh, prior.pc.rho = c(rho0, Prho), prior.pc.sig = c(sig0, Psig))
    nv <- mesh$n
    response.1 <- response[[1]]
    response.2 <- response[[2]]
    locs.1 <- locs[[1]]
    locs.2 <- locs[[2]]
    Ast1 <- inla.spde.make.A(mesh = mesh, loc = locs.1)
    Ast2 <- inla.spde.make.A(mesh = mesh, loc = locs.2)
    field.1 <- field.2 <-  copy.field <-1:nv
    if(!is.null(covariates)){
        m.1 <- make.covs(covariates[[1]])
        m.2 <- make.covs(covariates[[2]])
        cov.effects.1 <- m.1[[1]]
        cov.form.1 <- m.1[[2]]
        cov.effects.2 <- m.2[[1]]
        cov.form.2 <- m.2[[2]]
        stk.1 <- inla.stack(data=list(y=cbind(response.1,NA)),
                             A=list( Ast1,1,1),
                            effects=list(field.1 = field.1, beta0 = rep(1,nrow(locs.1)),cov.effects = cov.effects.1))
        stk.2 <- inla.stack(data=list(y=cbind(NA,response.2)),
                             A=list( Ast2,Ast2,1,1),
                            effects=list(field.2 = field.2,copy.field = copy.field, alpha0 = rep(1,nrow(locs.2)),cov.effects = cov.effects.2))
        stack <- inla.stack(stk.1,stk.2)
        x = "\"field.1\""
        formula = paste("y", "~  0 + beta0 + alpha0 +", cov.form.1, cov.form.2,
                        " + f(field.1, model=spde)",
                        "+ f(field.2, model=spde)",
                    "+ f(copy.field, copy =", x, ",fixed=FALSE, hyper = hyper )")
    }else{
        stk.1 <- inla.stack(data=list(y=cbind(response.1,NA)),
                             A=list( Ast1,1),
                             effects=list(field.1 = field.1, beta0 = rep(1,nrow(locs.1))))
        stk.2 <- inla.stack(data=list(y=cbind(NA,response.2)),
                             A=list( Ast2,Ast2,1),
                            effects=list(field.2 = field.2, copy.field = copy.field,alpha0 = rep(1,nrow(locs.2))))
        stack <- inla.stack(stk.1,stk.2)
        formula <- y ~ 0 + beta0 + alpha0 + f(field.1, model=spde) +
            f(field.2, model=spde) +
            f(copy.field, copy = "field.1", fixed=FALSE, hyper = hyper )
    }
    result <- inla(as.formula(formula), family = family,
                   data = inla.stack.data(stack),
                   control.predictor=list(A=inla.stack.A(stack)),
                   control.inla = control.inla,
                   control.compute = control.compute,
                   verbose = verbose,
                   ...)
    return(result)

}


#' spatio-temporal model fitting 
geo.spatial.j.temporal.fit <-function(mesh, locs, response, covariates, temp, family, control.time, control.inla,
                              hyper, control.compute, sig0, Psig, rho0, Prho, verbose, ...){
    spde <- inla.spde2.matern.new(mesh, prior.pc.rho = c(rho0, Prho), prior.pc.sig = c(sig0, Psig))
    nv <- mesh$n
    temp <- temp # temporal dimension
    temp.1 <- temp[[1]]
    temp.2 <- temp[[2]]
    k.1 <- (mesh.t1 <- inla.mesh.1d(temp.1))$n
    k.2 <- (mesh.t2 <- inla.mesh.1d(temp.2))$n
    ## create projection matrix for loacations
    response.1 <- response[[1]]
    response.2 <- response[[2]]
    locs.1 <- locs[[1]]
    locs.2 <- locs[[2]]
    Ast1 <- inla.spde.make.A(mesh = mesh, loc = locs.1, group = temp.1, n.group = k.1)
    Ast2 <- inla.spde.make.A(mesh = mesh, loc = locs.2, group = temp.2, n.group = k.2)
    field.1 <- inla.spde.make.index('field.1', n.spde = spde$n.spde, group = temp, n.group = k.1)
    field.2 <- inla.spde.make.index('field.2', n.spde = spde$n.spde, group = temp, n.group = k.2)
    copy.field <- inla.spde.make.index('copy.field', n.spde = spde$n.spde, group = temp.2, n.group = k.2)
    if(!is.null(covariates)){
        m.1 <- make.covs(covariates[[1]])
        m.2 <- make.covs(covariates[[2]])
        cov.effects.1 <- m.1[[1]]
        cov.form.1 <- m.1[[2]]
        cov.effects.2 <- m.2[[1]]
        cov.form.2 <- m.2[[2]]
        stk.1 <- inla.stack(data=list(y=cbind(response.1,NA)),
                            A=list( Ast1,1,1),
                            effects=list(field.1 = field.1, beta0 = rep(1,nrow(locs.1)),cov.effects = cov.effects.1))
        stk.2 <- inla.stack(data=list(y=cbind(NA,response.2)),
                            A=list( Ast2,Ast2,1,1),
                            effects=list(field.2 = field.2,copy.field = copy.field,
                                         alpha0 = rep(1,nrow(locs.2)),cov.effects = cov.effects.2))
        stack <- inla.stack(stk.1,stk.2)
        x = "\"field.1\""
        cov.form <- paste(cov.form.1,"+",cov.form.2)
        formula = paste("y", "~  0 + beta0 + alpha0 +", cov.form,
                        " + f(field.1, model=spde, group = field.1.group, control.group = control.time)",
                        "+ f(field.2, model=spde, group = field.2.group , control.group = control.time)",
                    "+ f(copy.field, copy =", x, ",fixed=FALSE, hyper = hyper )")
    }else{
        stk.1 <- inla.stack(data=list(y=cbind(response.1,NA)),
                            A=list( Ast1,1),
                            effects=list(field.1 = field.1, beta0 = rep(1,nrow(locs.1))))
        stk.2 <- inla.stack(data=list(y=cbind(NA,response.2)),
                            A=list( Ast2,Ast2,1),
                            effects=list(field.2 = field.2, copy.field = copy.field,alpha0 = rep(1,nrow(locs.2))))
        stack <- inla.stack(stk.1,stk.2)
        formula <- y ~ 0 + beta0 + alpha0 + f(field.1, model=spde) +
            f(field.2, model=spde) +
            f(copy.field, copy = "field.1", fixed=FALSE, hyper = hyper )
    }
    result <- inla(as.formula(formula), family = family,
                   data = inla.stack.data(stack),
                   control.predictor=list(A=inla.stack.A(stack)),
                   control.inla = control.inla,
                   control.compute = control.compute,
                   verbose = verbose,
                   ...)
    return(result)
}



geo.spatial.j.nl.temporal.fit <-function(mesh, locs, response, covariates, temp, family, control.time, control.inla,
                              hyper, non.linear, control.compute,sig0, Psig, rho0, Prho,  verbose, ...){
    spde <- inla.spde2.matern.new(mesh, prior.pc.rho = c(rho0, Prho), prior.pc.sig = c(sig0, Psig))
    nv <- mesh$n
    temp <- temp # temporal dimension
    temp.1 <- temp[[1]]
    temp.2 <- temp[[2]]
    k.1 <- (mesh.t1 <- inla.mesh.1d(temp.1))$n
    k.2 <- (mesh.t2 <- inla.mesh.1d(temp.2))$n
    ## create projection matrix for loacations
    response.1 <- response[[1]]
    response.2 <- response[[2]]
    locs.1 <- locs[[1]]
    locs.2 <- locs[[2]]
    u <- non.linear[[1]][["random.effect"]]
    u.mod <- non.linear[[1]][["model"]]
    u.2 <- non.linear[[2]][["random.effect"]]
    u.mod.2 <- non.linear[[2]][["model"]]
    Ast1 <- inla.spde.make.A(mesh = mesh, loc = locs.1, group = temp.1, n.group = k.1)
    Ast2 <- inla.spde.make.A(mesh = mesh, loc = locs.2, group = temp.2, n.group = k.2)
    field.1 <- inla.spde.make.index('field.1', n.spde = spde$n.spde, group = temp, n.group = k.1)
    field.2 <- inla.spde.make.index('field.2', n.spde = spde$n.spde, group = temp, n.group = k.2)
    copy.field <- inla.spde.make.index('copy.field', n.spde = spde$n.spde, group = temp.2, n.group = k.2)
    if(!is.null(covariates)){
        m.1 <- make.covs(covariates[[1]])
        m.2 <- make.covs(covariates[[2]])
        cov.effects.1 <- m.1[[1]]
        cov.form.1 <- m.1[[2]]
        cov.effects.2 <- m.2[[1]]
        cov.form.2 <- m.2[[2]]
        stk.1 <- inla.stack(data=list(y=cbind(response.1,NA)),
                            A=list( Ast1,1,1,1),
                            effects=list(field.1 = field.1, beta0 = rep(1,nrow(locs.1)),cov.effects = cov.effects.1,
                                         random = u))
        stk.2 <- inla.stack(data=list(y=cbind(NA,response.2)),
                            A=list( Ast2,Ast2,1,1,1),
                            effects=list(field.2 = field.2,copy.field = copy.field,
                                         alpha0 = rep(1,nrow(locs.2)),cov.effects = cov.effects.2, random.2 = u.2))
        stack <- inla.stack(stk.1,stk.2)
        x = "\"field.1\""
        cov.form <- paste(cov.form.1,"+",cov.form.2)
        formula = paste("y", "~  0 + beta0 + alpha0 + f(random, model = u.mod) + f(random.2, model = u.mod.2) +", cov.form,
                        " + f(field.1, model=spde, group = field.1.group, control.group = control.time)",
                        "+ f(field.2, model=spde, group = field.2.group , control.group = control.time)",
                    "+ f(copy.field, copy =", x, ",fixed=FALSE, hyper = hyper )")
    }else{
        stk.1 <- inla.stack(data=list(y=cbind(response.1,NA)),
                            A=list( Ast1,1,1),
                            effects=list(field.1 = field.1, beta0 = rep(1,nrow(locs.1)),
                                         random = u))
        stk.2 <- inla.stack(data=list(y=cbind(NA,response.2)),
                            A=list( Ast2,Ast2,1,1),
                            effects=list(field.2 = field.2, copy.field = copy.field,alpha0 = rep(1,nrow(locs.2)),
                                         random.2 = u.2))
        stack <- inla.stack(stk.1,stk.2)
        formula <- y ~ 0 + beta0 + alpha0 + f(field.1, model=spde) +  f(random, model = u.mod) + f(random.2, model = u.mod.2) +
            f(field.2, model=spde) +
            f(copy.field, copy = "field.1", fixed=FALSE, hyper = hyper )
    }
    result <- inla(as.formula(formula), family = family,
                   data = inla.stack.data(stack),
                   control.predictor=list(A=inla.stack.A(stack)),
                   control.inla = control.inla,
                   control.compute = control.compute,
                   verbose = verbose,
                   ...)
    return(result)
}
                         
#' Function that extracts the ''random fields'' of the model fitted.
#'
#' Plots the estimated random fields or parameter posterior densities from an object returned by
#' \link{mark.pp.fit}().
#'
#' @param x A fitted model from \link{mark.pp.fit}().
#' @param mesh the mesh used in the model fit.
#' @param n.t numeric, the number of time points.
#' @param sd Logical, if \code{FALSE} means of random fields aer returned.
#' @param plot Logical, if \code{TRUE} the returned matricies (either SD or Mean of
#' random fields are plotted.
#' @param spatial.polygon Optional, if a spatial polygon of the domain is supplied, only
#' values of the random field within the domain will be returned
#' @param ... additional graphical parameters
#'  @importFrom fields image.plot
#'
#' @export
find.fields <- function(x = NULL, mesh = NULL, n.t = NULL, sd = FALSE, plot = FALSE, spatial.polygon = NULL,dims = c(100,100),...){
    if(is.null(attributes(x)$mesh) & is.null(mesh)){
        stop("no mesh has been supplied")}
    if(!is.null(attributes(x)$mesh)){mesh <- attributes(x)$mesh}else{mesh <- mesh}
    proj <- inla.mesh.projector(mesh,dims = dims)
    if(!is.null(spatial.polygon)) inside <- inwin(proj,as.owin(spatial.polygon))
    spde <-inla.spde2.matern(mesh = mesh, alpha = 2)
    fields <- summary(x)$random.names[summary(x)$random.model=="SPDE2 model" | summary(x)$random.model=="Copy"]
    idx <- which(summary(x)$random.model=="SPDE2 model" | summary(x)$random.model=="Copy")
    n <- length(fields)
    if(!is.null(n.t)){
        t <- n.t
        means <- list()
        for (i in idx[1]:idx[n]){
            means [[i-idx[1]+1]] <- lapply(1:t, function(j) { r <- inla.mesh.project(proj, field = x$summary.random[[i]]$mean[1:spde$n.spde + (j-1)*spde$n.spde]);  if(!is.null(spatial.polygon)) r[!inside] <- NA; return(r)})
        }
        sds <- list()
        for (i in idx[1]:idx[n]){
            sds [[i-idx[1]+1]] <- lapply(1:t, function(j) {r <- inla.mesh.project(proj, field = x$summary.random[[i]]$sd[1:spde$n.spde + (j-1)*spde$n.spde]);if(!is.null(spatial.polygon)) r[!inside] <- NA;  return(r)})
        }
        if(plot){plot.fields( x = x, mesh = mesh, n.t = n.t, sd = sd, spatial.polygon = spatial.polygon,...)}
    }else{
        means <- list()
        for (i in idx[1]:idx[n]){
            means[[i-idx[1]+1]] <- inla.mesh.project(proj,x$summary.random[[i]]$mean)
            if(!is.null(spatial.polygon)) means[[i-idx[1]+1]][!inside] <- NA; 
            }
        sds <- list()
        for (i in idx[1]:idx[n]){
            sds[[i-idx[1]+1]] <- inla.mesh.project(proj,x$summary.random[[i]]$sd)
            if(!is.null(spatial.polygon)) sds[[i-idx[1]+1]][!inside] <- NA; 
            }
        if(plot){plot.fields( x = x, mesh = mesh, n.t = n.t, sd = sd, spatial.polygon = spatial.polygon,...)}
    }
    ifelse(sd,return(sds),return(means))
}
## find which parts of the random field are inside the supplied spatial polygon
inwin<-function(proj, window){
    e<-expand.grid(proj$x,proj$y)
    o<-inside.owin(e[,1],e[,2],window)
    o<-matrix(o,nrow=length(proj$x))
}
#####################

##extract values of random field at specified locations
extract <- function(x = NULL, locs = NULL,mesh = NULL, t = NULL, sd = FALSE){
    fields <- summary(x)$random.names[summary(x)$random.model == 
                                                "SPDE2 model" | summary(x)$random.model == "Copy"]
    idx <- which(summary(x)$random.model == "SPDE2 model" | summary(x)$random.model == 
                                                                      "Copy")
    n <- length(fields)
    A <- inla.mesh.project(mesh = mesh,loc = locs)$A
    if (!is.null(t)) {
        means <- list()
        for (i in idx[1]:idx[n]) {
            means[[i - idx[1] + 1]] <- lapply(1:t, function(j) {
                meshNodes <- x$summary.random[[i]]$mean[1:mesh$n + (i-1)*mesh$n]
                means[[i]] <- drop(A%*%meshNodes)
            })
        }
        sds <- list()
        for (i in idx[1]:idx[n]) {
            sds[[i - idx[1] + 1]] <- lapply(1:t, function(j) {
                meshNodes <- x$summary.random[[i]]$mean[1:mesh$n + (i-1)*mesh$n]
                sds[[i]] <- drop(A%*%meshNodes)
            })
        }
    } else {
        means <- list()
        for (i in idx[1]:idx[n]) {
            meshNodes <- x$summary.random[[i - idx[1] + 1]]$mean
            means[[i - idx[1] + 1]] <-  drop(A%*%meshNodes)
         
        }
        sds <- list()
        for (i in idx[1]:idx[n]) {
            meshNodes <- x$summary.random[[i - idx[1] + 1]]$sd
            sds[[i - idx[1] + 1]] <-  drop(A%*%meshNodes)
         
        }    
    }
    names(sds) <- names(means) <- fields
    ifelse(sd, return(sds), return(means))
}
