This supplementary material illustrate fitting the model discussed in
*Detecting Subnational Diffusion Processes of Lethal Terrorism: Global
Study, 2010-2016* using a subset of the data discussed in that article.
Please note that this material should be treated as an illustration of
the method only. Due to the computational time required to fit the
models discussed in the article the illustration here uses only a subset
of the data and a coarser mesh to enable readers to run through an
example and amend for their own data.

The package , available [here](https://github.com/cmjt/lgcpSPDE),
contains the data and the functionality to fit the model deailed herein.

    ## to install the lgcpSPDE package from github
    devtools::install_github("cmjt/lgcpSPDE")

    library(lgcpSPDE)
    ## only using a subset of the data for illustration
    data <- terrorism_aggregate[terrorism_aggregate$iyear <= 2011, 
        ]

    ## x, y, z locations (latitude and longitude projected onto
    ## the unit sphere)
    locs <- cbind(data$x.coord, data$y.coord, data$z.coord)
    ## create mesh of the world projected onto the unit sphere
    bdry <- inla.sp2segment(world)
    bdry$loc <- inla.mesh.map(bdry$loc, projection = "longlat", inverse = TRUE)
    ## course mesh for illustration
    mesh <- inla.mesh.2d(boundary = bdry, loc = locs, max.edge = c(5, 
        50)/180, cutoff = 5/180)  ## plot(mesh) to vizualise
    ## Create a named data frame of covariates
    covariates <- data.frame(population = data$pop, time.to.city = data$tt, 
        luminosity = data$lum)
    ## create temporal indecies
    time <- data$iyear - min(data$iyear) + 1

To fit this crude model run the following code.

    fit <- geo.fit(mesh = mesh, locs = locs, response = data$total, 
        covariates = covariates, control.time = list(model = "rw1", 
            param = list(theta = list(prior = "pc.prec", param = c(1, 
                0.01)))), temp = time, family = "poisson", sig0 = 0.2, 
        rho0 = 0.01, Prho = 0.9, control.compute = list(waic = TRUE, 
            config = TRUE, openmp.strategy = "huge"), control.inla = list(int.strategy = "eb", 
            strategy = "gaussian", diagonal = 100))

The function is a wrapper for 's (Finn Lindgren, Havard Rue (2015))
function and has the following main arguments:

-   **mesh**: a \`\`mesh'' object i.e. Delauney triangulation of the
    domain,
-   **locs**: a matrix of observation locations, where each row
    corresponds to the observation.
-   **response**: a vector of response variable, each corresponds to the
    spatial locations in .
-   **temp**: (optional) a numeric vector specifying a temporal index
    for each observation (starting at 1.....T).
-   **covariates**: (optional) a named data.frame of covariates.
-   **family**: a character vector specifying the assumed likelihood of
    the response, by default is "gaussian".
-   **control.time**: (optional) supplied if the argument is given to
    fit a spatio-temporal model. This argument controls the model and
    prior put on the hyperparameters of the model for the temporal
    component of the spatio-temporal model. By default this is which is
    a pc.prior put on the rho coefficient of a AR(1) model with
    P(rho&gt;0)=0.9.
-   **control.inla**: a list which controls the fitting procedures INLA
    uses. By default this is for quick and dirty fitting.
    **control.compute**: a list of fit statistics the user wants INLA to
    return. By default this is .
-   **non.linear**: (optional) should be used if the user requires a
    non-linear covariate to be included in the model Must be supplied as
    a named list with elements a numeric vector of the random effect
    indices, and the random effect model the user wishes to use for

### References

Finn Lindgren, Havard Rue. 2015. “Bayesian Spatial Modelling with
R-INLA.” *Nature Materials* 19 (63). Journal of Statistical Software:
1–25. <http://www.jstatsoft.org/v63/i19/>.
