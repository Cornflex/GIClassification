#setwd("C:\\Users\\D059331\\Desktop\\DM GIC\\data\\")
##Kristen's wd
setwd('C:\\Users\\Owner\\Documents\\Portugal\\Sat_Image_Classification\\data\\data')
require('raster')
require(sp)
require(rgdal)
require(randomForest)
require(rmarkdown)
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



##the tif file from Joel is a stack because it has 4 bands, or layers
rasterJ<-brick("img\\J_04AUG14112729-M2AS-000000137917_01_P001_etrs89.TIF")
#rasterE<-stack("E_04SEP24113435-M2AS-000000152724_01_P001_etrs89.TIF")
#rasterA<-stack("A_05SEP22114039-M2AS-005509561050_01_P001_etrs89.TIF")



# plot map
plotRGB(rasterJ, 3,2,1)

# choose a sector for clustering
ext <- drawExtent()
sector <- crop(rasterJ, ext)
plotRGB(sector, 3, 2, 1)

# add normdiff layer for better vegetation recognition
sector[[5]] <- normdiff(sector)
# plot each layer separately
plot(sector)
plotRGB(sector, 3,2,1)

# check for NA's
summary(sector)

# convert to dataframe and check for NA's
#sector_df <- as.data.frame(sector)
#summary(sector_df)

# for some reason, now there are NA's in 2nd layer
# lets remove it
# plot(sector)
# sector <- brick(
# 	sector[[1]],
# 	sector[[3]],
# 	sector[[4]],
# 	sector[[5]]
# )
# sector_df <- as.data.frame(sector)
# summary(sector_df)

# use only nir and ndvi layer for k-means
sector_mod <- brick(sector[[4]],sector[[5]])

# perform kmeans and plot result
performKMeans(sector_mod, 6)
performKMeans(sector,6)

#Using focal function to calculate mean, stddev values 
#for a neighborhood of cells
means<-mov.fun(sector,3,mean, normalize=FALSE, verbose=TRUE)
sds<-mov.fun(sector,3,sd, normalize=FALSE, verbose=TRUE)
#convert output of above to brick and summarize to check for na values
stk<-stack(means,sds)
sectorMovings<-brick(stk)
summary(sectorMovings)
#We see that there are equal NA values in each layer
#We get dimension information to determine number of cells in total in each layer
sectorMovings
#Plot the layers to see if the NA values are noticeable (grouped together for example)

#We see that NA values are a very small percentage of total cells (less than one percent) of all cells
plot(sectorMovings)
#they are not visible in the plots

#use trim to remove NA values
#from documentation: Trim (shrink) a Raster* object by removing outer rows and columns that all have the same value (e.g. NA). 
sectorMovingsTrim<-trim(sectorMovings, padding=0, values=NA)
#Check output for NA values
summary(sectorMovingsTrim)
#No NA values in any layers
sectorMovingsTrim
#Plote output of trim to see if any issues
plot(sectorMovingsTrim)
#It looks good. The moving function may have created NA values around edges, maybe where a full
#three by three grid could not be made.

#Resize extent of sector brick to match the mean trimmed so they can be stacked
newExtent <- extent(sectorMovingsTrim)
sectorModCrop <- crop(sector_mod,newExtent)

# #Can trim solve the issue of layer 2 having NA values when it converts to data frame?
# sector_df <- as.data.frame(sector[[2]])
# summary(sector_df)
# sector_trim<-trim(sector_df, padding=0, values=NA)
# summary(sector_trim)
# plot(sector[[2]])
# #Nope, doesn't remove and NA values, this means the NA values are not around the edges
# #What happeens if we normalize the layer?
# sector_norm<-raster.scale(sector[[2]])
# summary(sector_norm)

#test kmeans on just the movings stack, which is ten layers, means and sds of the original five layers
performKMeans(sectorMovingsTrim,6)
#

#combine movings stack with original 5 layer stack to get 15 layers
stk <- stack(sectorMovingsTrim,sectorModCrop)
sector_all<-brick(stk)
sector_all
#Normalize the stack
sector_all_norm<-raster.scale(sector_all)
sector_all_norm
#Perform kmeans on this normalized stack
performKMeans(sector_all_norm,6)
#does not look good

#try kmeans on just the normalized mean stack, which had good results when 
#not normalized
performKMeans(sectorMeansTrim,6)
sectorMeansTrimNorm<-raster.scale(sectorMeansTrim)
sectorMeansTrimNorm
performKMeans(sectorMeansTrim,6)

#cant try kmeans on non normalized stack of ten
performKMeans(sector_all,6)

#ok, go back and do the felix number of bands
#add these to 



  ### supervised


rasterJ
trainShapes <- readOGR(dsn="C:\\Users\\D059331\\Desktop\\DM GIC\\data\\shp\\trn", layer="J_treino_QB_point")

nrow(trainShapes)
plotRGB(rasterJ, 3, 2, 1)


rasterAtTrain <- extract(rasterJ,trainShapes)
train <- data.frame(rasterAtTrain, trainShapes)
trainFactor <- as.factor(train[,7])
#levels(trainFactor)
#[1] "4"  "7"  "12" "21" "22" "31" "32" "34" "35" "36"
# 10 levels
colors <- c("#CC0000", "#FFD700", "#556B2F", "#008B8B", "#191970", "#8A2BE2", "#D8BFD8", "#8B4513", "#000000", "#FF6347")
plot(trainShapes["Label"], add=TRUE, col=colors[trainFactor])
# light blue = water, red = city

forest <- randomForest(
	x=train[,1:4],
	y=trainFactor
)
plot(forest)
varImpPlot(forest, type=1)
#predict sector
plotRGB(sector, 3,2,1)
#prediction <- predict(rasterJ, forest, ext=ext, filename="ext1.tif", type="response", index=1, na.rm=TRUE, progress="window", overwrite=TRUE)
prediction <- predict(rasterJ, forest, ext=ext, type="response", index=1, na.rm=TRUE, progress="window", overwrite=TRUE)
prediction <- raster("RfClassPred.tif")
plotRGB(sector,3,2,1)
plot(prediction, add=TRUE, col=colors[trainFactor])


plotRGB(rasterJ, 3,2,1)
plot(trainShapes, add=TRUE)
head(trainShapes)
head(trainShapes@data)
head(v)


## with complete training data
trainShapesTotal <- readOGR(dsn="C:\\Users\\D059331\\Desktop\\DM GIC\\data\\shp\\trn", layer="J_treino_QB_Tot_point")
plotRGB(rasterJ, 3, 2, 1)
plot(trainShapesTotal["Label"], add=TRUE, col=color)
rasterAtTrainTotal <- extract(rasterJ,trainShapesTotal)
trainTotal <- data.frame(rasterAtTrainTotal, trainShapesTotal)
trainTotalFactor <- as.factor(trainTotal[,7])
#levels(trainTotalFactor)
#[1] "4"  "7"  "12" "21" "22" "31" "32" "34" "35" "36"
# 10 levels
colors <- c("#CC0000", "#FFD700", "#556B2F", "#008B8B", "#191970", "#8A2BE2", "#D8BFD8", "#8B4513", "#000000", "#FF6347")
plot(trainShapesTotal["Label"], add=TRUE, col=colors[trainTotalFactor])
forest <- randomForest(
	x=trainTotal[,1:4],
	y=)
)


# measure accuracy to determine optimal model fit
relevantTrainShapes <- crop(trainShapes, ext)
predictionMatch <- extract(prediction, relevantTrainShapes)
trainError <- calculateError(predictionMatch, relevantTrainShapes$Label)