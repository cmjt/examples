This supplementary material illustrate the fitting of the models
discussed in *Detecting Subnational Diffusion Processes of Lethal
Terrorism: Global Study, 2010-2016* using a subset of the data discussed
in that article. Please note that this material should be treated as an
illustration of the method only. Due to the computational time required
to fit the models discussed in the article the illustration here uses
only a subset of the data and a coarser mesh to enable reader's to run
through an example and amend for their own data.

Fitting code available either in the package available
[here](https://github.com/cmjt/lgcpSPDE) or functions are provided in
the supplementary
[functions.r](https://github.com/cmjt/examples/blob/master/materials/functions.r)
file (as per [Python, A., et al,
2018](https://rss.onlinelibrary.wiley.com/doi/abs/10.1111/rssa.12384)).

Hotspot evaluation code was written by [A.
Python](andre.python@bdi.ox.ac.uk).

    ## load required packages
    require(INLA)  ## model fitting
    require(spatstat)  ## functions for spatial data
    require(sp)
    require(GISTools)
    require(raster)
    require(rgeos)
    require(maptools)

    ## example data supplied with this material
    load("pnas_example.RData")
    ## fiting functions supplied with this material
    source("functions.r")

Model fitting example
---------------------

    head(example.df)  ## data for Afghanistan and Pakistan

    ##           x         y         z country country.idx latitude longitude
    ## 1 0.3463239 0.8395703 0.4185469     PAK           2     24.5      67.5
    ## 2 0.3364243 0.8437153 0.4182859     PAK           2     24.5      68.5
    ## 3 0.3136284 0.8523466 0.4185002     PAK           2     24.5      70.0
    ## 4 0.4237748 0.7994831 0.4257250     PAK           2     25.0      62.0
    ## 5 0.4202658 0.8015146 0.4253835     PAK           2     25.0      62.5
    ## 6 0.3544914 0.8352477 0.4203536     PAK           2     25.0      67.0
    ##   year year.idx population  time.to.city luminosity total
    ## 1 2014        1 -0.3558181  0.0005527337 -0.2549383     1
    ## 2 2015        2 -0.3304876 -0.0636544096 -0.9392858     1
    ## 3 2015        2 -0.3613648  0.5665875776 -0.2874247     1
    ## 4 2014        1 -0.3671028  1.5319930208 -1.3627242     1
    ## 5 2015        2 -0.3664247  1.1245760187 -0.7493187     3
    ## 6 2014        1  3.8020456 -0.4888692221  1.1150337   277

    plot(sp)  ## plot Afghanistan and Pakistan boundry
    points(example.df$longitude, example.df$latitude, pch = 20)  ## centroid points

![](diffusion_files/figure-markdown_strict/show%20data-1.png)

    ## make mesh; this is an example the mesh only
    bdry <- inla.sp2segment(sp)
    bdry$loc <- inla.mesh.map(bdry$loc, projection = "longlat", inverse = TRUE)  ## project to unit spere representing Earth
    mesh <- inla.mesh.2d(boundary = bdry, max.edge = c(4, 7)/300, 
        cutoff = 4/300)  ## WARNING MESH FAR TOO COARSE. FOR ILLUSTRATION ONLY
    plot(mesh, asp = 1.75)  ## aspect ratio contrls tilt of the plot

![](diffusion_files/figure-markdown_strict/fit%20model-1.png)

    ## set up model info
    locs <- as.matrix(example.df[, 1:3])
    response <- example.df$total
    covariates <- data.frame(population = example.df$population, 
        time.to.city = example.df$time.to.city, luminosity = example.df$luminosity)
    t.idx <- example.df$year.idx
    ## the model fitted here ignores country levelrandom effects
    ## to speed up compuation time call fitting function
    fit <- geo.fit(mesh = mesh, locs = locs, response = response, 
        covariates = covariates, control.time = list(model = "rw1", 
            param = list(theta = list(prior = "pc.prec", param = c(1, 
                0.01)))), temp = t.idx, verbose = FALSE, family = "poisson")
    ## WARNING NOT A GOOD FIT DUE TO COARSENESS OF MESH. FOR
    ## ILLUSTRATION ONLY

Finding "Hotspots"
------------------

#### Find the random fields

    k <- 2  ## number of years
    resol <- c(1440, 720)  ## set the resolution of the grid to detect hotspot and diffusion
    fields <- find.fields(fit, mesh, n.t = k, spatial.polygon = sp, 
        dims = resol)
    sdfields <- find.fields(fit, sd = TRUE, mesh, n.t = k, spatial.polygon = sp, 
        dims = resol)
    proj <- inla.mesh.projector(mesh, projection = "longlat", dims = resol)

#### Fitted values

    ## Crude way to extract fitted values and upper/lower limits
    surf <- list()
    xmean <- list()
    xsd <- list()
    bdq975 <- list()
    bdq25 <- list()
    bdq975surf <- list()
    bdq25surf <- list()
    for (j in 1:k) {
        xmean[[j]] <- fields[[1]][[j]]
        xsd[[j]] <- sdfields[[1]][[j]]
        surf[[j]] <- poisson(link = "log")$linkinv(fit$summary.fix[2, 
            1] * covariates$population + fit$summary.fix[3, 1] * 
            covariates$time.to.city + fit$summary.fix[4, 1] * covariates$luminosity + 
            fields[[1]][[j]] + fit$summary.fix[1, 1])
        bdq975[[j]] <- xmean[[j]] + 1.96 * xsd[[j]]
        bdq25[[j]] <- xmean[[j]] - 1.96 * xsd[[j]]
        bdq975surf[[j]] <- poisson(link = "log")$linkinv(fit$summary.fix[2, 
            1] * covariates$population + fit$summary.fix[3, 1] * 
            covariates$time.to.city + fit$summary.fix[4, 1] * covariates$luminosity + 
            bdq975[[j]] + fit$summary.fix[1, 1])
        bdq25surf[[j]] <- poisson(link = "log")$linkinv(fit$summary.fix[2, 
            1] * covariates$population + fit$summary.fix[3, 1] * 
            covariates$time.to.city + fit$summary.fix[4, 1] * covariates$luminosity + 
            bdq25[[j]] + fit$summary.fix[1, 1])
    }

    min <- min(unlist(lapply(surf, min, na.rm = TRUE)))
    min

    ## [1] 0.0008019408

    max <- max(unlist(lapply(surf, max, na.rm = TRUE)))
    max

    ## [1] 42147.27

    logmin <- log(quantile(response, probs = 0.95))
    logmin

    ##      95% 
    ## 3.401197

    logmax <- log(max)
    logmax

    ## [1] 10.64893

    brk <- seq(logmin, logmax, by = 0.1)
    m <- length(brk) - 1
    targetcy <- sp[(sp$iso_a3 == "AFG" | sp$iso_a3 == "PAK"), ]

#### Construct diffusion process

    # make diffusion processs
    hot <- list()
    for (j in 1:k) {
        hot[[j]] <- ifelse(log(bdq25surf[[j]]) > logmin, 1, NA)  ## put NA everywhere except where there are hotspots (useful for next steps)
    }
    ## create vectors with projection values (x,y) and hotspot
    ## values (z)
    hotdf <- list()
    hotr <- list()
    for (j in 1:k) {
        hotdf[[j]] <- list()
        hotdf[[j]]$x <- proj$x  ## longitude
        hotdf[[j]]$y <- proj$y  ## latitude
        hotdf[[j]]$z <- hot[[j]]
        ## convert vectors to raster
        hotr[[j]] <- raster(hotdf[[j]])
    }
    ## create dataframe from raster
    hotspotdf <- list()
    for (j in 1:k) {
        hotspotdf[[j]] <- as.data.frame(hotr[[j]], xy = TRUE)  ## add xy option to extract lon and lat
        names(hotspotdf[[j]]) <- c("lon", "lat", "hot")
    }
    ## delete hotpsot values if NA
    for (j in 1:k) {
        hotspotdf[[j]] <- hotspotdf[[j]][complete.cases(hotspotdf[[j]]), 
            ]
    }
    ## create polygon of hotspots from raster cells
    hotrpoly <- list()
    for (j in 1:k) {
        hotrpoly[[j]] <- rasterToPolygons(hotr[[j]], fun = NULL, 
            n = 16, na.rm = TRUE, digits = 5, dissolve = TRUE)
        hotrpoly[[j]] <- unionSpatialPolygons(hotrpoly[[j]], rep(1, 
            length(hotrpoly[[j]])))  ## put adjacent polygons together
        hotrpoly[[j]] <- disaggregate(hotrpoly[[j]])  ## disaggregate into multiple polygons
    }
    ## optional: remove too small hotspots
    tau <- 0.25  ## (0.25 correspons to PRIO-grid cell area)
    areas <- list()
    maxareas <- list()
    index <- list()
    for (j in 1:length(hotrpoly)) {
        ## calculate areas of all polygons included in hotspots for
        ## each year
        areas[[j]] <- sapply(slot(hotrpoly[[j]], "polygons"), function(x) sapply(slot(x, 
            "Polygons"), slot, "area"))
        maxareas[[j]] <- sapply(areas[[j]], max)
        ## in case of possible divided polygons keep the one with
        ## highest area (is required for following steps)
        index[[j]] <- which(maxareas[[j]] > tau)  ## find position where area >tau
        hotrpoly[[j]] <- hotrpoly[[j]][index[[j]], ]  ## keep polygon with area>tau
    }
    ## Create buffer
    hotbuff <- list()
    hotadj <- list()
    diffnet <- list()  ## the final result clipped within sp (if option below not done)
    for (j in 1:k) {
        hotbuff[[j]] <- createSPComment(hotrpoly[[j]])
        hotbuff[[j]] <- gBuffer(hotbuff[[j]], byid = FALSE, id = NULL, 
            width = 0.4545455, quadsegs = 1, capStyle = "SQUARE")
        hotbuff[[j]] <- gDifference(hotbuff[[j]], hotrpoly[[j]])
        proj4string(hotbuff[[j]]) <- CRS(proj4string(sp))  ## gives reference
        hotadj[[j]] <- gIntersection(hotbuff[[j]], sp)
        hotadj[[j]] <- createSPComment(hotadj[[j]])  ## problem with holes...
        hotadj[[j]] <- disaggregate(hotadj[[j]])  ## from one to multiple polygons
        proj4string(hotadj[[j]]) <- CRS(proj4string(sp))  ## gives reference
        diffnet[[j]] <- hotadj[[j]]
        diffnet[[j]] <- disaggregate(diffnet[[j]])
    }
    ## END NEIGHBOURHOOD CHOICE intersect grid with diffusion
    ## areas
    diffgrid <- list()
    grid <- list()
    gridpolygon <- list()
    drops <- c("layer")  ##  list of col names
    for (j in 1:k) {
        grid[[j]] <- raster(extent(diffnet[[j]]))
        res(grid[[j]]) <- 0.5
        proj4string(grid[[j]]) <- proj4string(sp)
        gridpolygon[[j]] <- rasterToPolygons(grid[[j]])
        diffgrid[[j]] <- intersect(diffnet[[j]], gridpolygon[[j]])
        diffgrid[[j]]@data$ID <- seq.int(nrow(diffgrid[[j]]@data))
        diffgrid[[j]] <- diffgrid[[j]][, !(names(diffgrid[[j]]) %in% 
            drops)]  ## remove unuseful variable
    }
    ## diffusion
    surflowdf <- list()
    surflowr <- list()
    for (j in 1:k) {
        surflowdf[[j]] <- list()
        surflowdf[[j]]$x <- proj$x
        surflowdf[[j]]$y <- proj$y
        surflowdf[[j]]$z <- bdq25surf[[j]]  ## using lower bound CI
        ## convert dataframe to raster
        surflowr[[j]] <- raster(surflowdf[[j]])
        proj4string(surflowr[[j]]) <- CRS(proj4string(sp))  ## required 
        surflowr[[j]] <- crop(surflowr[[j]], sp)
    }
    surfhighdf <- list()
    surfhighr <- list()
    for (j in 1:k) {
        surfhighdf[[j]] <- list()
        surfhighdf[[j]]$x <- proj$x
        surfhighdf[[j]]$y <- proj$y
        surfhighdf[[j]]$z <- bdq975surf[[j]]  ## using lower bound CI
        ## convert dataframe to raster
        surfhighr[[j]] <- raster(surfhighdf[[j]])
        proj4string(surfhighr[[j]]) <- CRS(proj4string(sp))  ## required 
        surfhighr[[j]] <- crop(surfhighr[[j]], sp)
    }
    surfmeandf <- list()
    surfmeanr <- list()
    for (j in 1:k) {
        surfmeandf[[j]] <- list()
        surfmeandf[[j]]$x <- proj$x
        surfmeandf[[j]]$y <- proj$y
        surfmeandf[[j]]$z <- surf[[j]]  ## using lower bound CI 
        ## convert dataframe to raster
        surfmeanr[[j]] <- raster(surfmeandf[[j]])
        proj4string(surfmeanr[[j]]) <- CRS(proj4string(sp))  ## required 
        surfmeanr[[j]] <- crop(surfmeanr[[j]], sp)
    }
    ## get centroid of diffusion area polygons
    diffxy <- list()
    for (j in 1:k) {
        diffxy[[j]] <- gCentroid(diffgrid[[j]])
        diffxy[[j]] <- coordinates(diffxy[[j]])
    }
    rm(extract)
    meandiff <- list()
    highdiff <- list()
    lowdiff <- list()
    escpos <- list()
    escneg <- list()
    escdf <- list()
    ################## OPTIONAL threshold#####################
    thres <- 0
    ################## OPTIONAL threshold#####################
    for (j in 1:(k - 1)) {
        meandiff[[j]] <- unlist(extract(surfmeanr[[j + 1]], diffgrid[[j]], 
            fun = mean, na.rm = T)) - unlist(extract(surfmeanr[[j]], 
            diffgrid[[j]], fun = mean, na.rm = T))
        highdiff[[j]] <- unlist(extract(surfhighr[[j + 1]], diffgrid[[j]], 
            fun = mean, na.rm = T)) - unlist(extract(surfhighr[[j]], 
            diffgrid[[j]], fun = mean, na.rm = T))
        lowdiff[[j]] <- unlist(extract(surflowr[[j + 1]], diffgrid[[j]], 
            fun = mean, na.rm = T)) - unlist(extract(surflowr[[j]], 
            diffgrid[[j]], fun = mean, na.rm = T))
        escpos[[j]] <- ifelse(meandiff[[j]] > thres & highdiff[[j]] > 
            thres & lowdiff[[j]] > thres, 1, 0)
        escneg[[j]] <- ifelse(meandiff[[j]] < -1 * thres & highdiff[[j]] < 
            -1 * thres & lowdiff[[j]] < -1 * thres, 1, 0)
    }
    for (j in 1:(k - 1)) {
        ## Join mean values to polygon data
        escpos[[j]][is.na(escpos[[j]])] <- 0
        escneg[[j]][is.na(escneg[[j]])] <- 0
        diffgrid[[j]]@data$meandiff <- meandiff[[j]]
        diffgrid[[j]]@data$highdiff <- highdiff[[j]]
        diffgrid[[j]]@data$lowdiff <- lowdiff[[j]]
        diffgrid[[j]]@data$escpos <- escpos[[j]]
        diffgrid[[j]]@data$escneg <- escneg[[j]]
    }
    ## Crop diffusion process
    diffcrop <- list()
    for (j in 1:(k - 1)) {
        diffcrop[[j]] <- crop(diffgrid[[j]], extent(40, 75, 25, 40))
    }
    ## plot
    h <- 1200
    w <- 1.2 * h

#### Plot hotspot and diffusion process

    ## define colour palette
    coldis <- "green4"  ## dark (dissipation)
    colhot <- "red"  ## red (hotspot)
    coldif <- "green1"  ## light (diffusion)
    colnull <- "grey25"  ## light (diffusion)
    plot(diffcrop[[1]], col = ifelse(diffcrop[[1]]$escpos > 0, coldif, 
        ifelse(diffcrop[[1]]$escneg > 0, coldis, colnull)), border = "NA", 
        main = "", xlab = "", ylab = "", axes = F, xlim = c(60, 80), 
        ylim = c(20, 40))
    plot(hotrpoly[[1]], density = 10, angle = 45, col = colhot, border = "black", 
        add = TRUE)
    plot(sp, yaxt = "n", lwd = 3, add = TRUE, border = "grey25")
    plot(targetcy, yaxt = "n", lwd = 3, add = TRUE, border = "grey25")
    text(targetcy, targetcy@data$sov_a3, cex = 1.5)
    box()
    axis(1, col.axis = "black", at = seq(60, 80, by = 10), labels = seq(60, 
        80, by = 10), las = 2)
    axis(2, col.axis = "black", at = seq(20, 40, by = 5), labels = seq(20, 
        40, by = 5), las = 0)
    title(main = 2014, line = 0.5, cex.main = 1.5)
    mtext("Latitude", side = 2, line = 2, at = 30, cex = 1)
    mtext("Longitude", side = 1, line = 2, at = 70, cex = 1)
    legend("bottomright", title = "", inset = 0.005, bty = "n", c("hotspot", 
        "diffusion", "dissipation", "no process"), fill = c(colhot, 
        coldif, coldis, colnull), density = c(10, NA, NA, NA), angle = c(45, 
        NA, NA, NA, NA))

![](diffusion_files/figure-markdown_strict/plot%20hot-1.png)
