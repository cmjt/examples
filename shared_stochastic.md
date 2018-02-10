Applications of marked LGCPs
----------------------------

All funtions used below are available in the `lgcpSPDE` package
[here](https://github.com/cmjt/lgcpSPDE), which can be installed by
using `devtools::install_github("cmjt/lgcpSPDE")`.

### Spatial dynamics of terrorism activity worldwide

The terrorism data is available within the `lgcpSPDE` package and can be
found by running `data(terrorism)`. In addition a
`SpatialPolygonsDataFrame` of the world used in the construction of the
mesh can be found by running `data(world)`.

Terrorism activity, 2010--2016, is overlaid onto a map of the world
below.

![](shared_stochastic_files/figure-markdown_strict/terrorism%20data%7D-1.png)

    head(terrorism)

    ##     nkill  latitude longitude iyear   x.coord     y.coord    z.coord fatal
    ## 1       4 35.077911 62.873525  2010 0.3731415  0.72835230 0.57468979     1
    ## 14      6 32.983300 67.966600  2010 0.3146843  0.77756565 0.54439456     1
    ## 20      1  5.116561  7.367093  2010 0.9877932  0.12771508 0.08918219     1
    ## 86      0 13.833110 20.834835  2010 0.9075024  0.34535935 0.23909462     0
    ## 113     0  5.316667 -4.033333  2010 0.9932317 -0.07003421 0.09266023     0
    ## 120     0  5.419303 -4.021559  2010 0.9930789 -0.06981835 0.09444371     0

    locs <- as.matrix(terrorism[,5:7])
    t.index <- (terrorism$iyear - min(terrorism$iyear)) + 1
    mark <- terrorism$fatal ## 1 = fatalities, 2 = no fatalities

    mesh <- make.mesh(locs = locs,mesh.pars = c(max = 3, min = 0.03,cutoff = 0.03),
                      sphere = TRUE,
                      spatial.polygon = world) ##over sphere with land constraint
    fit <- fit.marked.lgcp(mesh = mesh, locs = locs,t.index = t.index, 
                       mark = mark, mark.family = "binomial",
                       verbose=TRUE, 
                       hyper = list(theta = list(prior = "normal", param = c(0, 10))),
                       link = 2) ## use quick strategy first

    T.fit.imp <- fit.marked.lgcp(mesh = mesh, locs = locs, t.index = t.index, 
                                 mark = mark, mark.family = "binomial",
                                 verbose=TRUE, 
                                 hyper = list(theta = list(prior = "normal", param = c(0, 10))),
                                 control.mode = list(result = fit,restart = TRUE),
                                 control.inla = list(strategy='laplace'),
                                 link = 2)## improve using starting values from first model

### Analysis of eye movement data

The eye movement data is not supplied within the `lgcpSPDE` package as
we do not have permission to distribute. The data refer to the data
collected by Nathalie Guyader & Anne Guérin (GIPSA-lab, Grenoble,
France), and made available through the [RSS 2015
Challenge](https://rsschallenge.wordpress.com/the-2015-challenge/).
However, in this section we will describe the data and detail the
fitting of a spatiotemporal LGCP.

    head(eye_movement)

    ##          fx     fy duration fixation subject condition image    format
    ## 4273  -86.1  -11.0      150        1      n1    normal     1 landscape
    ## 4274 -146.3 -140.9      347        2      n1    normal     1 landscape
    ## 4275  283.0   56.1      163        3      n1    normal     1 landscape
    ## 4276   16.6  228.9      207        4      n1    normal     1 landscape
    ## 4277   80.4  228.6      130        5      n1    normal     1 landscape
    ## 4278 -207.4  222.2      294        6      n1    normal     1 landscape
    ##             sl        ta
    ## 4273        NA        NA
    ## 4274 143.17140 204.86458
    ## 4275 472.34256  65.35025
    ## 4276 317.53551 302.96940
    ## 4277  63.80071  90.26941
    ## 4278 287.87115 268.72609

    mesh <- inla.mesh.2d(loc = as.matrix(eye_movement[,1:2]),max.edge = 50,cutoff = 30)
    plot(mesh)

![](shared_stochastic_files/figure-markdown_strict/eye%20movement-1.png)

    t.index <- as.numeric(cut(eye_movement[,4],7))
    locs <- as.matrix(eye_movement[,1:2])
    mark <- eye_movement[,9]

    fit <- fit.marked.lgcp(mesh = mesh, locs = locs, t.index = t.index,
                           mark = mark, mark.family = "gamma",
                           prior.range = c(500,0.5),
                           prior.sigma = c(5,0.05),
                           verbose = TRUE, 
                           hyper = list(theta = list(prior = "normal", param = c(0, 10))),
                           link = 2)
