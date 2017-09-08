Simulated examples
==================

Parameter estimation for the Thomas, Matérn, and void process is carried
out by the R package **palm** (B. C. Stevenson (2017)), which is
available on CRAN. All models are fitted using either the **fit.ns()**
function, for NSPP, or the **fit.void()** function, for the void
process.

Below, code examples demonstrate using functionality of the **palm** (B.
C. Stevenson (2017)) package to simulate from and fit a two-dimensional
point process. This is done for of each of the types mentioned above.
The simulated data for each process are shown---parents are plotted at
grey crosses and daughters by black dots---along with the fitted (solid
lines) and empirical (dashed lines) Palm intensity functions.

Please note that parameters of the Palm intensity functions may differ,
in name only, from the descriptions in Jones-Todd et al. (2017).

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

    lims <- rbind(c(0, 1),c(0,1)) ## 2D limits of domain (i.e., unit square)
    thomas <- sim.ns(c(D = 7, lambda = 8, sigma = 0.05), lims = lims) ## simulate a Thomas process
    fit.thomas <- fit.ns(thomas$points,lims = lims, R = 0.5) ## fit a Thomas process
    coef(fit.thomas)

    ##           D      lambda       sigma 
    ## 10.34116799  4.83815432  0.03896552

![](/home/charlotte/Git/examples/CRC_point_process_files/figure-markdown_strict/plot%20thomas-1.png)

### Matérn process

The Palm intensity function of a Matérn process may be characterised by
the parameter vector **θ** = (*D*, *λ*, *τ*). Here *D* is the density of
parents, *λ* gives the expected number of daughters per parent, and *τ*
is the radius of the sphere, centered at a parent location, within which
daughters are uniformally scattered.

    lims <- rbind(c(0, 1),c(0,1)) ## 2D limits of domain (i.e., unit square)
    matern <- sim.ns(c(D = 7, lambda = 8, tau = 0.05), lims = lims,disp = "uniform") ## simulate a Matern process
    fit.matern <- fit.ns(matern$points,lims = lims, R = 0.5, disp = "uniform") ## fit a Matern process
    coef(fit.matern)

    ##          D     lambda        tau 
    ## 3.43483934 6.78890506 0.04931027

![](/home/charlotte/Git/examples/CRC_point_process_files/figure-markdown_strict/plot%20matern-1.png)

Void point process
------------------

The Palm intensity function of a void process may be characterised by
the parameter vector
**θ** = (*D*<sub>*c*</sub>, *D*<sub>*p*</sub>, *τ*). Here
*D*<sub>*c*</sub> is the density of children, *D*<sub>*p*</sub> is the
density of parents, and *τ* is the radius of the voids centered at each
parent.

    lims <- rbind(c(0, 1),c(0,1)) ## 2D limits of domain (i.e., unit square)
    void <- sim.void(pars = c(Dc = 300, Dp = 10,tau = 0.075),lims = lims) ## simulate a void process
    fit.void <- fit.void(points = void$points, lims = lims, R = 0.5,
                         bounds = list(Dc = c(280,320), Dp = c(8,12), tau = c(0,0.2))) ## fit a void process
    coef(fit.void)

    ##           Dp           Dc          tau 
    ##  12.00000000 305.35979603   0.06703315

![](/home/charlotte/Git/examples/CRC_point_process_files/figure-markdown_strict/plot%20void-1.png)

Variance estimation
-------------------

Variance estimation is achieved via a parametric bootstrap, which can be
carried out through the use of the **boot.palm()** function. For
example, running **boot.palm(fit.thomas,N = 1000)** will perform 1000
bootstrap resamples of the fitted thomas process in order to estimate
standard errors of the parameters.

CRC data
========

References
==========

Jones-Todd, C. M, P Caie, J Illian, B. C Stevenson, Savage A, D
Harrison, and J Bown. 2017. “Identifying Unusual Structures Inherent in
Point Pattern Data and Its Application in Predicting Cancer Patient
Survival.” *arXiv Preprint arXiv:1705.05938*.

Stevenson, Ben C. 2017. *Palm: Fitting Point Process Models Using the
Palm Likelihood*. <https://github.com/b-steve/palm>.
