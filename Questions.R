# Questions.R
############################################################################
# SET UP
############################################################################

setwd("C:\\Users\\D059331\\Desktop\\DM GIC\\data\\")
require('raster')
require(sp)
require(rgdal)
require(randomForest)

raster.scale <- function(rst, norm="min-max") {
	require(raster)
	nbands <- nlayers(rst)
	if (norm == "min-max") {
		Min <- vector(mode="numeric", length=nbands)
		Max <- vector(mode="numeric", length=nbands)
		for (k in 1:nbands) {
			Min[k] <- rst[[k]]@data@min
			Max[k] <- rst[[k]]@data@max
		}
		return(scale(rst, center=Min, scale=Max-Min))
	}
	# by default mean-std normalization
	return(scale(rst))
}


# Compute normalized difference (e.g. NDVI)
normdiff <- function(stk, red=3, nir=4, stack.it=FALSE, normalize=FALSE, norm="min-max") {
	require(raster)
	dims <- dim(stk)
	max.bands = max(red,nir)
	if (dims[3] < max.bands)
		stop("Raster file does not have enough bands\n")
	nirBand <- getValues(stk[[nir]])
	redBand <- getValues(stk[[red]])
	ndviBand <- (nirBand - redBand) / (nirBand + redBand)
	if (normalize == TRUE) {
		if (norm == "min-max") {
			center <- min(ndviBand)
			ampl <- max(ndviBand)-center
		} else {
			center <- mean(ndviBand)
			ampl <- sd(ndviBand)
		}
		ndviBand <- scale(ndviBand, center=center, scale=ampl)
	}
	ndvi <- raster(nrows=dims[1], ncols=dims[2], crs=crs(stk), ext=extent(stk), resolution=res(stk), vals=ndviBand)
	if (stack.it == FALSE)
		return(ndvi)
	return(stack(stk, ndvi))
}

# Moving function
# This function computes the moving (average, standard-deviation, etc) of a raster dataset
mov.fun <- function(rst, window.size, fun, normalize=FALSE, verbose=FALSE) {
	require(raster)
	if (!is.integer(window.size))
		window.size <- abs(ceiling(window.size))
	# need to be odd
	if (window.size%%2 == 0)
		window.size <- window.size + 1
	w <- matrix(1, nrow=window.size, ncol=window.size)
	n <- nlayers(rst)
	stk <- focal(rst[[1]], w, fun, na.rm=TRUE)
	if (verbose == TRUE)
		cat(">> Band ", 1," of ", n, " done.\n")
	for (k in 2:n) {
		tmp <- focal(rst[[k]], w, fun, na.rm=TRUE)
		stk <- stack(stk, tmp)
		if (verbose == TRUE)
			cat(">> Band ", k, " of ", n, " done.\n")
	}
	if (normalize == TRUE)
		stk <- raster.scale(stk)
	return(stk)
}

# function to calculate and plot kmeans for a raster object
performKMeans <- function(inputRaster, noClusters) {
	raster_df <- as.data.frame(inputRaster)
	clustering <- kmeans(raster_df, noClusters, iter.max = 100, nstart = 10)
	clusterRaster <- raster(inputRaster)
	clusterRaster <- setValues(clusterRaster, clustering$cluster)
	plot(clusterRaster)
	return(clusterRaster)
}

calculateError <- function(prediction, actual) {
	trainDiff <- prediction - actual
	trainDiffCount <- 0
	for(i in 1:length(prediction)) {
		if(trainDiff[i] != 0) {
			trainDiffCount <- trainDiffCount + 1
		}
	}
	return(trainDiffCount  / length(trainDiff))
}

rasterJ<-brick("img\\J_04AUG14112729-M2AS-000000137917_01_P001_etrs89.TIF")


############################################################################
# Why does the second layer contain NA's after data frame conversion?
############################################################################
# plot map
plotRGB(rasterJ, 3,2,1)

# choose an sector for clustering
ext <- drawExtent()
sector <- crop(rasterJ, ext)
plotRGB(sector, 3, 2, 1)

# check for NA's
summary(sector)

# convert to dataframe and check for NA's
sector_df <- as.data.frame(sector)
summary(sector_df)

############################################################################
# Raster A and E seem to be in a strange color format.
# How can we work with them?
############################################################################
rasterA <- brick("img\\A_05SEP22114039-M2AS-005509561050_01_P001_etrs89.TIF")
summary(rasterA) # looks good
plotRGB(rasterA, 3, 2, 1) # Nothing
plot(rasterA) # reverse NA
rasterAScaled <- raster.scale(rasterA)
plot(rasterAScaled) # Nothing new
ndviA <- normdiff(rasterA)
plot(ndviA) # this works!

############################################################################
# Hints on unsupervised clustering?
############################################################################
