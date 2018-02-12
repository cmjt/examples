This supplementary material provides the code used to fit the
spatio-temporal marked log-Gaussian Cox process discussed in submitted
to . There is an onlie version of this tutorial also avaiable
[here](https://github.com/cmjt/examples/blob/master/Ecography_ESD.md).

Please note that the data cannot be supplied along with the
supplementary material due to the protection status of the species.
However, if readers do wish to access the data please contact [Andrea
Soriano Redondo](A.Soriano-Redondo@exeter.ac.uk).

The refered to below is a data frame of 65676 rows and 20 columns. Each
row corresponds to an observation, and the columns picked out for the
analysis correspond to the latitude and longitude of wetland locations,
sum of breeding crane pairs on each wetland, wetland density, wetland
area, a wetland density buffer and urban density buffer.

    ## required packages
    library(INLA)
    ## load data.
    wetland <- read.csv("Wetland.csv")

Code section below extracts wetland coordinates and creates time index
and binomial mark vector.

    ## get coordinates (latitude and longitude of wetland
    ## locations)
    locs <- cbind(wetland[, 14], wetland[, 15])
    plot(locs, main = "Wetland locations in England and Wales", axes = FALSE, 
        xlab = "", ylab = "")

![](Ecography_ESD_files/figure-markdown_strict/data-1.png)

    ## transform time from 2000-2015 to 1-16
    t.index <- (wetland$Year - min(wetland$Year)) + 1
    ## Crane presence as binomial mark (0 absent/ 1 present; sum
    ## of breeding crane pairs on each wetland)
    mark <- ifelse(wetland$SUM_PAIR == 0, 0, 1)
    mark.family <- "binomial"

This section sets up the covariate data frame as required for model
fitting

    ## Setting covariates
    wetland$PA_ratio <- wetland$Perimeter/wetland$Area
    ## wetland density, wetland area, a wetland density buffer and
    ## urban density buffer.
    cov1 <- wetland[, c(4, 12, 20, 21)]
    cov1$Area_sc <- scale(cov1$Area)
    cov1$PA_ratio_sc <- scale(cov1$PA_ratio)
    cov1$Wet_density_nosea_sc <- scale(cov1$Wet_density_buf_NoSea)
    cov1$Urb_density_nosea_sc <- scale(cov1$Urb_density_buf_NoSea)
    ## covariate data frame
    covariates <- data.matrix(cov1[, c(5, 6, 7, 8)])
    # function that takes in a named matrix of covariates with
    # ncol equal to the number of covariates returns a list
    # containing the effects ready to be read by (inla.stack) and
    # a covariate formula (ready to be read by inla)
    make.covs <- function(covariates) {
        n.covs <- ncol(covariates)
        for (i in 1:n.covs) {
            assign(colnames(covariates)[i], covariates[, i], envir = .GlobalEnv)
        }
        cov.effects <- sapply(colnames(covariates), get, simplify = FALSE)
        cov.form <- paste(colnames(covariates), collapse = " + ")
        return(list(cov.effects, cov.form))
    }
    m <- make.covs(covariates)
    cov.effects <- m[[1]]
    cov.form <- m[[2]]

The mesh is a requirement of the SPDE model fot the latent fields.

    ## make mesh for SPDE model
    mesh <- inla.mesh.2d(loc = locs, cutoff = 0.6, max.edge = c(0.2, 
        2))

This section sets priors used for hyperparameters of the model and sets
up latent structures.

    ## prior distributions prior for ar1 parameter
    prior.rho = list(theta = list(prior = "pccor1", param = c(0, 
        0.9)))
    ## prior for interaction parameter
    hyper = list(theta = list(prior = "normal", param = c(0, 10)))
    ## Define SPDE model for latent fields
    spde <- inla.spde2.matern(mesh = mesh, alpha = 2)
    ## number of observations
    n <- nrow(locs)
    ## number of mesh nodes
    nv <- mesh$n
    temp <- t.index  ## temporal dimentson
    k <- (mesh.t <- inla.mesh.1d(temp))$n  ## number of groups
    ## the response for the point pattern locations
    y.pp <- rep(0:1, c(k * nv, n))
    ## create projection matrix for loacations
    Ast <- inla.spde.make.A(mesh = mesh, loc = locs, group = temp, 
        n.group = k)
    ## effect for LGCP used for point pattern
    st.volume <- diag(kronecker(Diagonal(n = k), spde$param.inla$M0))
    expected <- c(st.volume, rep(0, n))
    ## create field indecies
    field.pp <- inla.spde.make.index("field.pp", n.spde = spde$n.spde, 
        group = temp, n.group = k)
    field.mark <- inla.spde.make.index("field.mark", n.spde = spde$n.spde, 
        group = temp, n.group = k)
    copy.field <- inla.spde.make.index("copy.field", n.spde = spde$n.spde, 
        group = temp, n.group = k)
    ## temporal model 'ar1'
    ctr.g <- list(model = "ar1", param = prior.rho)

Here the data stacks, a requirement of , are constructed. They contain
the information (i.e., covariates and fields) for each component of the
joint response model.

    ## Prepare data stacks for use in call to inla point process
    ## stack
    stk.pp <- inla.stack(data = list(y = cbind(y.pp, NA), e = expected), 
        A = list(rBind(Diagonal(n = k * nv), Ast)), effects = list(field.pp = field.pp))
    ## mark stack
    stk.mark <- inla.stack(data = list(y = cbind(NA, mark)), A = list(Ast, 
        Ast, 1), effects = list(field.mark = field.mark, copy.field = copy.field, 
        cov.effects = cov.effects))
    ## combine data stacks
    stack <- inla.stack(stk.pp, stk.mark)

Finally the formula is defined and is called to fit the joibt
spatio-temporalmarked log-Gaussian Coz process model.

    ## formula for model as given by equation X in article
    formula <- y ~ 0 + Area_sc + PA_ratio_sc + Wet_density_nosea_sc + 
        Urb_density_nosea_sc + f(field.pp, model = spde, group = field.pp.group, 
        control.group = ctr.g) + f(field.mark, model = spde, group = field.mark.group, 
        control.group = ctr.g) + f(copy.field, copy = "field.pp", 
        fixed = FALSE, hyper = hyper)
    fit <- inla(formula, family = c("poisson", mark.family), data = inla.stack.data(stack), 
        E = inla.stack.data(stack)$e, control.predictor = list(A = inla.stack.A(stack)), 
        control.inla = list(strategy = "gaussian", int.strategy = "eb"), 
        verbose = TRUE)  ## verbose = TRUE will print out modelling output to console