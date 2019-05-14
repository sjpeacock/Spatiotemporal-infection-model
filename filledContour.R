#My own filled contour function that does not use layout

filledContour<-function (x = seq(0, 1, length.out = nrow(z)), y = seq(0, 1,
    length.out = ncol(z)), z, xlim = range(x, finite = TRUE),
    ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE),
    levels = pretty(zlim, nlevels), nlevels = 20, color.palette = cm.colors,
    col = color.palette(length(levels) - 1), plot.title, plot.axes,
    key.title, key.axes, asp = NA, xaxs = "i", yaxs = "i", las = 1,
    axes = TRUE, frame.plot = axes, ...){
    mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar

    par(las = 1)

    mar <- mar.orig
    par(mar = mar)
    plot.new()
    plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, asp = asp)
    if (!is.matrix(z) || nrow(z) <= 1L || ncol(z) <= 1L)
      stop("no proper 'z' matrix specified")
    if (!is.double(z))
        storage.mode(z) <- "double"
    .filled.contour(as.double(x), as.double(y), z, as.double(levels), 
    								col = col)
    if (missing(plot.axes)) {
        if (axes) {
            title(main = "", xlab = "", ylab = "")
            Axis(x, side = 1)
            Axis(y, side = 2)
        }
    }
    else plot.axes
    if (frame.plot)
        box()
    if (missing(plot.title))
        title(...)
    else plot.title
    invisible()
}

plot.scale.bar<-function(zlim = range(z, finite = TRUE), 
												 levels = pretty(zlim, nlevels), nlevels = 20, color.palette = cm.colors, 
												 col = color.palette(length(levels) - 1), plot.title, plot.axes, 
												 key.title, key.axes, asp = NA, xaxs = "i", yaxs = "i", las = 1, 
												 axes = TRUE, frame.plot = axes, ...){
	
	# mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
	# # on.exit(par(par.orig))
	# w <- (3 + mar.orig[2L]) * par("csi") * 2.54
	# par(las = las)
	# mar <- mar.orig
	# mar[4L] <- mar[2L]
	# mar[2L] <- 1
	# par(mar = mar)
	par(mar=c(4,1,2,3))
	plot.new()
	plot.window(xlim = c(0, 1), ylim = range(levels), xaxs = "i", 
							yaxs = "i")
	rect(0, levels[-length(levels)], 1, levels[-1L], col = col)
	if (missing(key.axes)) {
		if (axes) 
			axis(4)
	}
	# else key.axes
	# box()
	# if (!missing(key.title)) 
	# key.title
	# mar <- mar.orig
	# mar[4L] <- 1
	# par(mar = mar)
	
}
