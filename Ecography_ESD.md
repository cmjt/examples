    ## required packages
    library(INLA)
    ## load data.
    wetland <- read.csv("~/Wetland_breeding_50000_200025_england_10bufNOSEA_mask_urban.csv")

Please note that the data cannot be supplied along with the
supplementary material. is a data frame of x rows and 16 columns. Each
row corresponds to an observation, and the columns correspond to x, y ,z
respectively.

    ## get coordinates
    locs <- cbind(wetland[, 14], wetland[, 15])
    ## transform time from 2000-2015 to 1-16
    t.index <- (wetland$Year - min(wetland$Year)) + 1
    ## Crane presence as binomial mark (0 absent/ 1 present)
    mark <- ifelse(wetland$SUM_PAIR == 0, 0, 1)
    mark.family <- "binomial"

    ## Setting covariates
    wetland$PA_ratio <- wetland$Perimeter/wetland$Area
    cov1 < -wetland[, c(4, 12, 20, 21)]
    cov1$Area_sc <- scale(cov1$Area)
    cov1$PA_ratio_sc <- scale(cov1$PA_ratio)
    cov1$Wet_density_nosea_sc <- scale(cov1$Wet_density_buf_NoSea)
    cov1$Urb_density_nosea_sc <- scale(cov1$Urb_density_buf_NoSea)
    ## covariate data frame
    covariates <- data.matrix(cov1[, c(5, 6, 7, 8)])

    ## make mesh for SPDE model
    mesh.pars <- c(max.edge.min = 0.2, max.edge.max = 2, cutoff = 0.6)
    mesh <- make.mesh(locs = locs, mesh.pars = mesh.pars)

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
