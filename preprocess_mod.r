# Basic functions to pre-process satellite images
# Contents:
# raster.scale - normalizes the dataset
# normdiff - computes normalized difference (e.g. NDVI, NDBI, NDWI)
# mov.fun - computes the moving function (= mean, std, etc) with a given window size (e.g. 3x3, 5x5, etc)
# Author: J. D. Silva (Feb 2012)

# Normalize a raster dataset

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
		ndviBand <- raster.scale(ndviBand, norm=norm)
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


##Kristen's stuff
setwd("C:\\Users\\m2015382\\Documents\\data\\img")

##the tif file from Joel is a stack because it has 4 bands, or layers
rasterA<-stack("A_05SEP22114039-M2AS-005509561050_01_P001_etrs89.TIF")



##Scaling data with Joel´s function
Ascale<-raster.scale(rasterA)
Ascale

#caculating NDVI 
ndvi=normdiff(Ascale)
ndvi

getValues(Ascale, row=10)

plot(ndvi)
plot(rasterA)
#to just plot one band - not working though...
plot(rasterA[[1]])

#playing with summary
summary(ndvi)
summary(rasterA)
summary(rasterA[[2]])
summary(Ascale)


#trying to create table from ndvi raster, not getting expected results, even though it 
#plots properly above
ndvitable<-getValues(ndvi)
head(ndvitable, n=10)
ndvitable


