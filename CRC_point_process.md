Simulated examples
==================

Parameter estimation for the Thomas, Matérn, and void process is carried
out by the R package **palm** (Stevenson (2017)), which is available on
CRAN. All models are fitted using either the **fit.ns()** function, for
NSPP, or the **fit.void()** function, for the void process.

Below, code examples demonstrate using functionality of the **palm**
(Stevenson (2017)) package to simulate from and fit a two-dimensional
point process. This is done for of each of the types mentioned above.
The simulated data for each process are shown---parents are plotted at
grey crosses and daughters by black dots---along with the fitted (solid
lines) and empirical (dashed lines) Palm intensity functions.

Please note that parameters of the Palm intensity functions may differ,
in name only, from the descriptions in Jones-Todd et al. (2018).

    library(palm)

Neyman-Scott point processes
----------------------------

Cluster processes where unobserved parent points randomly generate
daughters centered at their unobserved location. The spatial
distribution of daughters may eiher follow a bivariate normal
distribution (Thomas) or be uniformally distributed in a circle
(Mat'ern).

### Thomas process

The Palm intensity function of a Thomas process may be characterised by
the parameter vector **θ** = (*D*, *λ*, *σ*). Here *D* is the density of
parents, *λ* gives the expected number of daughters per parent, and *σ*
is the standard deviation of the bivariate normal distribution by which
daughters are scattered around their parents.

    set.seed(1234)
    lims <- rbind(c(0, 1),c(0,1)) ## 2D limits of domain (i.e., unit square)
    thomas <- sim.ns(c(D = 7, lambda = 8, sigma = 0.05), lims = lims) ## simulate a Thomas process
    fit.thomas <- fit.ns(thomas$points,lims = lims, R = 0.5) ## fit a Thomas process
    coef(fit.thomas)

    ##          D     lambda      sigma 
    ## 3.87800024 6.29007765 0.05679637

![](CRC_point_process_files/figure-markdown_strict/plot%20thomas-1.png)

### Matérn process

The Palm intensity function of a Matérn process may be characterised by
the parameter vector **θ** = (*D*, *λ*, *τ*). Here *D* is the density of
parents, *λ* gives the expected number of daughters per parent, and *τ*
is the radius of the sphere, centered at a parent location, within which
daughters are uniformally scattered.

    set.seed(2344)
    lims <- rbind(c(0, 1),c(0,1)) ## 2D limits of domain (i.e., unit square)
    matern <- sim.ns(c(D = 7, lambda = 8, tau = 0.05), lims = lims,disp = "uniform") ## simulate a Matern process
    fit.matern <- fit.ns(matern$points,lims = lims, R = 0.5, disp = "uniform") ## fit a Matern process
    coef(fit.matern)

    ##          D     lambda        tau 
    ## 5.43910354 5.05162971 0.04685375

![](CRC_point_process_files/figure-markdown_strict/plot%20matern-1.png)

Void point process
------------------

The Palm intensity function of a void process may be characterised by
the parameter vector
**θ** = (*D*<sub>*c*</sub>, *D*<sub>*p*</sub>, *τ*). Here
*D*<sub>*c*</sub> is the density of children, *D*<sub>*p*</sub> is the
density of parents, and *τ* is the radius of the voids centered at each
parent.

    set.seed(3454)
    lims <- rbind(c(0, 1),c(0,1)) ## 2D limits of domain (i.e., unit square)
    void <- sim.void(pars = c(Dc = 300, Dp = 10,tau = 0.075),lims = lims) ## simulate a void process
    fit.void <- fit.void(points = void$points, lims = lims, R = 0.5,
                         bounds = list(Dc = c(280,320), Dp = c(8,12), tau = c(0,0.2))) ## fit a void process
    coef(fit.void)

    ##           Dp           Dc          tau 
    ##   8.00503329 280.52032288   0.07164047

![](CRC_point_process_files/figure-markdown_strict/plot%20void-1.png)

CRC data
========

This section discusses fitting models to the colorectal cancer (CRC)
data discussed in Jones-Todd et al. (2018). Please note that the CRC
data cannot be supplied however to illustrate the methods in Jones-Todd
et al. (2018) we simulate data from a Thomas point process model.

![Illustration of one image of a patient's slide which enables the
pinpointing of nuclei. Plot i) is a composite immunofluorescence digital
image (red fluorescence highlights tumour cells and blue fluorescence
highlights all nuclei in the image). Plot ii) is an image analysis mask
overlay from automatic machine learnt segmentation of the digital image:
Plot iii) is the point pattern formed by the nuclei of the tumour
(black) and stroma (grey) cells shown in the previous two
images.](CRC_point_process_files/figure-markdown_strict/cancer.png)

The **palm** functions **fit.ns()** and **fit.void()** can also take
lists for and **lims** arquments. This allows us to estimate parameters
at the patient level for each tumour and stroma pattern in the CRC data.
Below we first simulate thirty NSPPs, fifteen for each of two patients
and illustrate fitting a Thomas process to the data.

    ## simulate fifteen point patterns for two "patients"
    N <- 15 ## number of patterns to simulate
    x <- list(rbind(c(0,1),c(0,1))) ## limits of one pattern
    lims <- rep(x, N) ## repeat these limits for each patters
    ## simulate for "first" patient
    sim <- sim.ns(c(D = 7, lambda = 8, sigma = 0.05), lims = lims)
    sim.points.one <- lapply(sim, function(x) x$points) ## get all daughter points
    str(sim.points.one)

    ## List of 15
    ##  $ : num [1:78, 1:2] 0.731 0.767 0.874 0.878 0.888 ...
    ##  $ : num [1:33, 1:2] 0.529 0.5 0.578 0.603 0.49 ...
    ##  $ : num [1:40, 1:2] 0.444 0.532 0.582 0.466 0.508 ...
    ##  $ : num [1:85, 1:2] 0.1348 0.1348 0.1606 0.151 0.0409 ...
    ##  $ : num [1:51, 1:2] 0.36 0.484 0.377 0.497 0.415 ...
    ##  $ : num [1:80, 1:2] 0.811 0.749 0.749 0.815 0.687 ...
    ##  $ : num [1:46, 1:2] 0.904 0.817 0.926 0.76 0.756 ...
    ##  $ : num [1:76, 1:2] 0.567 0.523 0.499 0.592 0.674 ...
    ##  $ : num [1:60, 1:2] 0.653 0.613 0.769 0.672 0.688 ...
    ##  $ : num [1:42, 1:2] 0.645 0.708 0.688 0.694 0.651 ...
    ##  $ : num [1:58, 1:2] 0.0726 0.026 0.0241 0.2717 0.3281 ...
    ##  $ : num [1:42, 1:2] 0.385 0.5 0.417 0.392 0.395 ...
    ##  $ : num [1:63, 1:2] 0.767 0.739 0.761 0.705 0.771 ...
    ##  $ : num [1:32, 1:2] 0.846 0.811 0.992 0.794 0.881 ...
    ##  $ : num [1:58, 1:2] 0.558 0.493 0.574 0.394 0.513 ...

    ## simulate for "second" patient
    sim <- sim.ns(c(D = 10, lambda = 4, sigma = 0.1), lims = lims)
    sim.points.two <- lapply(sim, function(x) x$points) ## get all daughter points
    str(sim.points.two)

    ## List of 15
    ##  $ : num [1:42, 1:2] 0.5241 0.4017 0.2367 0.0337 0.2431 ...
    ##  $ : num [1:19, 1:2] 0.0642 0.1441 0.371 0.234 0.3733 ...
    ##  $ : num [1:39, 1:2] 0.711 0.934 0.964 0.889 0.832 ...
    ##  $ : num [1:41, 1:2] 0.51 0.305 0.403 0.309 0.343 ...
    ##  $ : num [1:26, 1:2] 0.7556 0.7364 0.8178 0.9225 0.0182 ...
    ##  $ : num [1:34, 1:2] 0.984 0.867 0.659 0.752 0.718 ...
    ##  $ : num [1:45, 1:2] 0.496 0.729 0.471 0.542 0.307 ...
    ##  $ : num [1:38, 1:2] 0.903 0.827 0.703 0.44 0.639 ...
    ##  $ : num [1:34, 1:2] 0.877 0.968 0.898 0.9 0.991 ...
    ##  $ : num [1:12, 1:2] 0.1479 0.0216 0.1411 0.0308 0.1113 ...
    ##  $ : num [1:32, 1:2] 0.113 0.251 0.153 0.264 0.392 ...
    ##  $ : num [1:39, 1:2] 0.705 0.679 0.825 0.577 0.727 ...
    ##  $ : num [1:53, 1:2] 0.218 0.211 0.169 0.271 0.204 ...
    ##  $ : num [1:37, 1:2] 0.1453 0.0975 0.097 0.7379 0.6825 ...
    ##  $ : num [1:46, 1:2] 0.0335 0.1141 0.0647 0.012 0.0679 ...

    ## fit models
    fit.one <- fit.ns(points = sim.points.one, lims = lims, R = 0.5)
    coef(fit.one)

    ##          D     lambda      sigma 
    ## 7.95838999 6.80621184 0.04658654

    fit.two <- fit.ns(points = sim.points.two, lims = lims, R = 0.5)
    coef(fit.two)

    ##          D     lambda      sigma 
    ## 9.66838372 3.57232273 0.09250772

    ## plot Palm intensity for all patterns
    plot(fit.one)

![](CRC_point_process_files/figure-markdown_strict/plot%20palms-1.png)

    plot(fit.two)

![](CRC_point_process_files/figure-markdown_strict/plot%20palms-2.png)

References
==========

Jones-Todd, C. M, P Caie, J Illian, B. C Stevenson, A Savage, D
Harrison, and J Bown. 2018. “Identifying Unusual Structures Inherent in
Point Pattern Data and Its Application in Predicting Cancer Patient
Survival.” *Statistics in Medicine* DOI:10.10002/sim.8046.

Stevenson, B. C. 2017. *Palm: Fitting Point Process Models Using the
Palm Likelihood*. <https://github.com/b-steve/palm>.
