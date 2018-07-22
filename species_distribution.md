



This supplementary material provides the code used to fit the
spatio-temporal marked log-Gaussian Cox process discussed in *Estimating
species distribution in highly dynamic populations using point process
models* submitted to *Ecography*. There is an online version of this
tutorial also available
[here](https://github.com/cmjt/examples/blob/master/species_distribution.md).

Please note that the data used in the article cannot be supplied along
with the supplementary material due to the protection status of the
species (if readers do wish to access the data please contact [Andrea
Soriano Redondo](A.Soriano-Redondo@exeter.ac.uk)). However, we supply
data simulated from the model to illustrate the methodology detailed in
the article. Alongside this we supply functionality so that the readers
can fit the models discussed; the procedure for this is outlined below.

The functionality is provided in the suite of R functions, **lgcpSPDE**,
available [here](https://github.com/cmjt/lgcpSPDE); to install the
**lgcpSPDE** package from GitHub run the following R code.

    devtools::install_github("cmjt/lgcpSPDE") 

The data
--------

The simulated example data used here, contained in the R object
**cranes**, is included in the **lgcpSPDE** package. This dataframe has
5052 wetland observations and 10 variables:

-   **Wetland\_Identity**, a numeric wetland ID;
-   **Lon**, longitude epicentre location of wetland;
-   **Lat**, latitude epicentre location of wetland;
-   **Area**, Area of wetland in *m*<sup>2</sup>;
-   **Perimiter**, Wetland perimiter in *m*;
-   **Wet\_density\_buf\_NoSea**, Surrounding wetland density;
-   **Urb\_density\_buf\_NoSea**, Surrouning urban density;
-   **Year**, Year of observation (i.e., 2014 or 2015);
-   **mark**, A simulated binary mark indicating presence of a Grus Grus
    breeding pair at the wetland;
-   **PA\_ratio**, Wetland perimiter to area ratio.

To load and inspect the data run the R following code.

    library(lgcpSPDE)
    data(cranes)
    str(cranes) ## to see further details run ?cranes

    ## 'data.frame':    5052 obs. of  10 variables:
    ##  $ Wetland_Identity     : int  35 37 41 67 68 80 97 99 119 120 ...
    ##  $ Lon                  : num  -5.21 -5.19 -5.17 -5.25 -5.26 ...
    ##  $ Lat                  : num  50 50 50 50 50 ...
    ##  $ Area                 : int  349375 201250 405625 206875 301250 608749 628751 229374 193750 104375 ...
    ##  $ Perimeter            : int  9450 6450 8250 5400 7050 8600 9350 8300 5350 1700 ...
    ##  $ Wet_density_buf_NoSea: num  0.01086 0.01064 0.00665 0.00853 0.00643 ...
    ##  $ Urb_density_buf_NoSea: num  0.0191 0.018 0.0155 0.0324 0.0311 ...
    ##  $ Year                 : int  2014 2014 2014 2014 2014 2014 2014 2014 2014 2014 ...
    ##  $ mark                 : int  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ PA_ratio             : num  0.027 0.032 0.0203 0.0261 0.0234 ...

The point pattern (of wetlands) and assocciated mark (binary presence or
absence of a pair of breeding cranes) is shown below.

![Epicentre locations of wetlands in England (i.e., potentially suitable
habitats for Grus Grus). Black plotting characters indicate the wetlands
that were (simulated to be) occupied by a breeding pair in 2014 an/or
2015.](species_distribution_files/figure-markdown_strict/plot-1.png)

Modelling
---------

In order to fit a marked spatio-temporal model to this example data, as
detailed in the article, first construct the required R objects as
below.

    locs <- as.matrix(cranes[,2:3]) ## matrix of wetland epicentre locations
    mark <- cranes$mark ## crane breeding pair presence as binary mark (0 absent/ 1 present)
    table(mark,cranes$Year) ## 10 occupied sites in 2014, 20 in 2015

    ##     
    ## mark 2014 2015
    ##    0 2516 2506
    ##    1   10   20

    mark.family <- "binomial" 
    t.index <- (cranes$Year-min(cranes$Year)) + 1 ## create year index
    ## Create covariate data frame (scaled)
    covariates <- data.frame(Area_sc = scale(cranes$Area), PA_ratio_sc = scale(cranes$PA_ratio),
                              Wet_density_nosea_sc = scale(cranes$Wet_density_buf_NoSea),
                              Urb_density_nosea_sc = scale(cranes$Urb_density_buf_NoSea))
    ## construct the mesh ( course just for illustration purposes)
    mesh <- inla.mesh.2d(loc = locs, cutoff = 0.6, max.edge = c(0.2,2)) 

The triangulation (mesh) is an integral part of the INLA-SPDE
methodology. The mesh provides the basis on which the approximation of
the Gaussian random fields are constructed. Choices regarding the
sparsity of the mesh as well as its shape are not simple to make. The
mesh should (among other things)

-   reflect the spatial structure you wish to capture, and

-   provide a decent representation of the spatial domain.

Futher details on mesh construction can be found in Lindgren, Rue, and
Lindström (2011).

This supplementary material is simply to illustrate how the models
discussed in the article were fitted and to enable readers to go through
the same procedure in a timely manner; therefore, the mesh we construct
(shown below) is badly formed and too sparse---simply aid in computation
time.

![Example Delauney triangulation (mesh) required for the marked
spatio-temporal log-Gaussian Cox process model. Overlain are the wetland
locations, black occupies and grey un-occupied, and a spatial polygon of
the UK (not including Northen
Irelan).](species_distribution_files/figure-markdown_strict/mesh-1.png)

To fit the model the **fit.marked.lgcp()** from **lgcpSPDE** is called.
This function offers a number of options, run

    args(fit.marked.lccp)

in R to see all options.

To fit the model discussed in *Estimating species distribution in highly
dynamic populations using point process models* to this simulated
example data we use the following arguments,

-   **mesh**, the Delauney trangulation discussed above;
-   **locs**, a 2 times n matrix of wetland locations;
-   **t.index**, a vector of year indecies for the n observations;
-   **covariates**, a data frame of the named covariates;
-   **mark**, a vector of the binary mark values;
-   **mark.family**, the distribution that the marks are assumed to be a
    realisation of (i.e., "binomial");
-   **prior.range** & **prior.sigma**, peanalised complexity priors on
    the spatial range and standard deviation (see Appendix A2 for more
    details);
-   **pp.int**, logical indicating if the linear predictor for the point
    pattern component of the model should contain an intercept term. If
    TRUE the parameter *α*<sub>0</sub> is estimated;
-   **mark.int**, logical indicating if the linear predictor for the
    mark component of the model should contain an intercept term. If
    TRUE the parameter *β*<sub>0</sub> is estimated.

The following finction call in R will fit a marked spatio-temporal
log-Gaussian Cox process model, as detailed in the article, to the
example data.

    fit <- fit.marked.lgcp(mesh = mesh, locs = locs, t.index = t.index, 
                           covariates = covariates, mark = mark,
                           mark.family = mark.family,
                           prior.range = c(4,0.5),
                           prior.sigma = c(1,0.05),
                           pp.int = TRUE, mark.int = TRUE)

### Model inference

Use the **summary** utility function to view the posterior estimates of
the model parameters. The fixed effect parameters relate to the
coefficients of the fixed effects. In the model fitted here these
parameters are

-   *α*<sub>0</sub>, the intercept for the point pattern component of
    the model;
-   *β*<sub>0</sub>, the intercept for the mark component of the model;
-   the coefficients of the named covariates of the mark component.

The parameters (hyper-parameters) of the latent fields are the spatial
range and standard deviation of the named field along with the *ρ*
parameter for the assumed AR(1) temporal process. The fields in this
model are **field.pp**, **z**(**s**, *t*) in the article (i.e., reflects
the spatially varying intensity of wetlands), and **field.mark**,
**g**(**s**, *t*) in the artile (i.e., conditional on the spatial
intensity of the wetlands this reflects the spatially varying process of
presence of a breeding crane pairs). As the Gaussian field
**z**(**s**, *t*) is shared between both components of the model the
interaction parameter *β* (named *Beta for copy.field*) is estimated.

    summary(fit)$fixed[,1:5] ## print out some summary information for fixed effects

    ##                          mean     sd 0.025quant 0.5quant 0.975quant
    ## beta0                -13.8249 2.0758   -17.9005 -13.8250    -9.7527
    ## alpha0                 0.8576 1.2571    -1.6106   0.8576     3.3237
    ## Area_sc                0.2284 0.0971     0.0378   0.2284     0.4188
    ## PA_ratio_sc           -2.8761 0.5596    -3.9748  -2.8761    -1.7784
    ## Wet_density_nosea_sc  -0.6579 0.4874    -1.6149  -0.6579     0.2982
    ## Urb_density_nosea_sc  -1.0324 0.7318    -2.4691  -1.0324     0.4031

    summary(fit)$hyperpar[,1:5] ## some summary information for parameters of the latent fields

    ##                           mean     sd 0.025quant 0.5quant 0.975quant
    ## Range for field.pp      3.8064 0.5551     2.9607   3.7187     5.1155
    ## Stdev for field.pp      2.4933 0.2991     2.0279   2.4488     3.1927
    ## GroupRho for field.pp   0.9998 0.0002     0.9993   0.9999     1.0000
    ## Range for field.mark    1.8996 0.5650     1.0554   1.8094     3.2554
    ## Stdev for field.mark    3.7561 0.7497     2.6117   3.6376     5.5364
    ## GroupRho for field.mark 0.9547 0.0498     0.8185   0.9706     0.9980
    ## Beta for copy.field     1.0098 0.3072     0.4169   1.0046     1.6260

    ## extract the random fields from the model object over the n.t = 2 years 
    ## This returns a list of matrices
    fields <- find.fields(fit, mesh, n.t = 2) 

The **fields** object above is a named list of matrices of each Gaussian
field for each year. Below **z**(**s**, *t*) (left) and
**g**(**s**, *t*) (right) are plotted for the first year on the link
scale.

![Posterior means of the random fields on the link scales. Left plot
shows the random field of the point process intensity. Right plot shows
the mark specific random field once the spatial structure of the
wetlands, and fixed effects, have been accounted
for.](species_distribution_files/figure-markdown_strict/plots-1.png)

References
----------

Lindgren, Finn, Håvard Rue, and Johan Lindström. 2011. “An Explicit Link
Between Gaussian Fields and Gaussian Markov Random Fields: The
Stochastic Partial Differential Equation Approach.” *Journal of the
Royal Statistical Society: Series B (Statistical Methodology)* 73 (4).
Wiley Online Library: 423–98.
