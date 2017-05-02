#!/usr/bin/env Rscript

##########################################################################################
#                        MaxEntAnalyzer Rscript v1.0, April 2017                         #
#  Copyright (c)2017 Justinc C. Bagley, Virginia Commonwealth University, Richmond, VA,  #
#  USA; Universidade de Brasília, Brasília, DF, Brazil. See README and license on GitHub #
#  (https://github.com/justincbagley) for further information. Last update: April 27,    #
#  2017. For questions, please email jcbagley@vcu.edu.                                   #
#
#  This script was developed by JCB based on a MaxEnt script that was originally written #
#  by Pietro Mello and Leonardo Gonçalves Tedeschi. JCB translated this code to English, #
#  added new code for ENMeval, and changed the analysis of the threshold value.          #
##########################################################################################

###### INSTRUCTIONS.
##--Make sure the order of coordinates is given as Longitude, followed by Latitude. Mixing
##--this up will render results (if any) useless. Also, divide test and training data. A
##--buffer will be required around occurrence points for species/lineages with less than
##--10 occurrences. Finally, as a general rule of thumb, only use test and training data
##--with bootstrapping when there are at least 15 occurrence points for the species/lineage.
#
##--This script assumes that bioclimatic data layers have been downloaded and prepped and
##--are in current working dir. Layers also assumed to have 0.3 arc-min (=1 degree) resolution.

setwd("./")

#rm(list=ls())	## Uncomment this line to clear the workspace when this script is sourced.

##--Load needed library, R code, or package stuff. Install package if not present.
#source('MaxEntAnalyzer.R', chdir = TRUE)
packages <- c('sp', 'grid', 'lattice', 'ape', 'raster', 'maptools', 'phyloclim', 'dismo','rJava','rgdal')
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
    install.packages(setdiff(packages, rownames(installed.packages())))
}

library(maptools) # Loads dependencies sp, grid, and lattice.
library(phyloclim) # Loads dependencies ape and raster.
library(dismo)
library(rgdal)
library(rJava)
Maxent

###### DEFINE EXTENT OF STUDY AREA TO BE MODELED, AND LOAD GIS DATA LAYERS.
##--Here, we define the geographical extent of the study area, assuming that the extent is 
##--the same for all models run in this instance of running the script. Change as needed if 
##--extent varies across species/lineages, or other units (ext = degrees, ext2 = UTM meters).
ext <- extent(-97,-77,6,25)
ext2 <- extent(-1500000,300000,500000,25000000)

##--Create object with Lambert Azimuthal Equal Area projection (lat and long from Hijmans 
##--and Graham, 2007):
crs.laea <- CRS("+proj=laea +lat_0=0 +lon_0=-80 +ellps=WGS84 +units=m +no_defs")

bioclim.vars <- list.files(path = paste(system.file(package = "dismo"),'/ex/wc0.5', sep = ''), pattern = 'bil', full.names = TRUE)
bioclim.vars <- stack(bioclim.vars)
bioclim.vars.crop <- crop(bioclim.vars, ext)
bioclim.vars.LAEA <- projectRaster(bioclim.vars.crop, crs=crs.laea)
plot(bioclim.vars.LAEA$bio_1)

##--Load paleoenvironmental data layers for the LGM:
paleo.LGM <- list.files(path = paste(system.file(package = "dismo"),'/ex/LGM', sep = ''), pattern = 'tif', full.names = TRUE)
paleo.LGM
paleo.LGM <- stack(paleo.LGM)
paleo.LGM.crop <- crop(paleo.LGM, ext)
paleo.LGM.LAEA <- projectRaster(paleo.LGM.crop, crs=crs.laea)

##--Load paleoenvironmental data layers for the LIG:
paleo.LIG <- list.files(path = paste(system.file(package = "dismo"),'/ex/LIG', sep = ''), pattern = 'bil', full.names = TRUE)
paleo.LIG
paleo.LIG <- stack(paleo.LIG)
paleo.LIG.crop <- crop(paleo.LIG, ext)
paleo.LIG.LAEA <- projectRaster(paleo.LIG.crop, crs=crs.laea)


###### RESCALE THE DATA LAYERS AS NEEDED.
##--Rescale all of the rasters. Here we'll rescale only the LGM variables so that they have
##--the same cell sizes as those of the LIG.
paleo.LGM.LAEA <- resample(paleo.LGM.LAEA, paleo.LIG.LAEA, method = "bilinear")


###### PLOT SHAPE FILES OF LAYERS ON LOCAL MACHINE.
## e.g. 'shape' <- readOGR(dsn= ..., layer= ...)
CA <- readOGR(dsn="C:/Users/leonardo/Desktop/Maxent", layer="ac_pol", proj4string(crs.laea))
CA2 <- spTransform(CA, CRS=crs.laea)
bioclim.vars.CA <- mask(bioclim.vars.LAEA,CA2)
col <- c("BIO1","BIO12","BIO16","BIO17","BIO18","BIO19","BIO2","BIO3","BIO4","BIO5","BIO6","BIO7","BIO8","BIO9")
names(bioclim.vars.CA) <- col
bioclim.vars.CA <- subset(bioclim.vars.CA, col)
paleo.LGM.CA <- mask(paleo.LGM.LAEA,CA2)
paleo.LIG.CA <- mask(paleo.LIG.LAEA,CA2)
plot(bioclim.vars.CA$BIO1)
plot(paleo.LGM.CA$BIO1)
plot(bioclim.LIG.CA$BIO1)



#@@@@@@@@@@@@@@@@@@@@@@@@@@@ SINGLE SPECIES/LINEAGE ANALYSES @@@@@@@@@@@@@@@@@@@@@@@@@@@@#


##########################################################################################
#                          ***  Alfaro cultratus ANALYSIS  ***                           #
##########################################################################################
### I. Load species occurrences input data file:
sp1 <- read.table("Alfaro cultratus.txt", h=T)
head(sp1)
summary(sp1)
tail(sp1)
coordinates(sp1) <- ~lon+lat
proj4string(sp1) <- CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs")
sp1 <- spTransform(sp1, CRS=crs.laea)
class(sp1)
head(sp1)
plot(sp1)

### II. Clip the layers and plot the occurrence points for this species:
##--Geographically resize/cut the map, where o = east, l = west, s = south, n = north; shape =
##--map that will be utilized.
plot(CA2,xlim=c(-1500000,500000),ylim=c(500000,3000000), main= "", bty="o", axes=TRUE, col="lightgreen")
title (main = "Presence points for Alfaro cultratus")

##--Plot the observed points:
points(sp1$lon, sp1$lat, col='blue', pch=20, cex=0.75)

##--If you want, you can highlight the points by making them red in color:
points(sp1$lon, sp1$lat, col='red', cex=0.75)

### III. Prepare testing and training data for MaxEnt:
##--Here, we divide the data points into testing and training data, based on Heuberty (1994):
group <- kfold(sp1, k=5)
pres_train <- sp1[group != 1, ]		## Training data
pres_test <- sp1[group == 1, ]		## Testing data

## Como ja feito anteriormente, para definir o tamanho da area de predicao
##para melhora a velocidade de processamento ext = extent(W,E,S,N)coordenadas dos vertices do retangulo

###  IV. RUN MAXENT. ### 
system.file("java", package="dismo")
jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='')
if (file.exists(jar)) {
	mx.AC <- maxent(bioclim.vars.CA, pres_train, overwrite=TRUE, args=c("-J","-r"), 
	path= "C:/Users/leonardo/Desktop/Maxent/modelresult/Alfaro cultratus/present")
	## Here, -J calls Jackknife method, -r asks if you want to 'overwrite' previous/existing results files
	plot(mx.AC, xlab="Percentage of contribution", ylab="Abiotic variable")
}

##--Background data points:
bg <- randomPoints(bioclim.vars.CA, n=10000, ext=CA2, extf = 1.25)

if (file.exists(jar)) {
	pvtest <- data.frame(extract(bioclim.vars.CA, pres_test))
	avtest <- data.frame(extract(bioclim.vars.CA, bg))
	testp <- predict(mx.AC, pvtest) 
	testa <- predict(mx.AC, avtest) 
	e <- evaluate(p=testp, a=testa, bioclim.vars.CA)
	e
	px.AC = predict(mx.AC, bioclim.vars.CA, ext=ext2, overwrite=TRUE, progress='',  
	filename="C:/Users/leonardo/Desktop/Maxent/modelresult/Alfaro cultratus/present")
#
	plot(px.AC, main='Maxent, raw values,Alfaro cultratus', bty="o")
}

if (file.exists(jar)) {
	## SSSmax was the method that varies the resampling of pseudo-absence (points) less, 
	## yielding a more constant threshold (Liu et al.,2013)
	threshold.AC <- e@t[which.max(e@TPR + e@TNR)]
	pre_abs.AC <- px.AC > threshold.AC
	plot(pre_abs.AC, main="Presence/absence-SSSmax, Alfaro cultratus")
	plot(CA, add=TRUE, border='dark grey')
	points(pres_train, pch='+')} else { plot(1)
}

##--Plot the ROC curve: 
threshold(e)
plot(e, 'ROC', main= "ROC plot", sub= "Alfaro cultratus", bty="o")

## LGM ##
px.LGM.AC <- predict(mx.AC, paleo.LGM.CA, overwrite=TRUE, type= "prop", ext=ext2, progress='', 
filename="C:/Users/leonardo/Desktop/Maxent/modelresult/Alfaro cultratus/LGM")
#
plot(px.LGM.AC, main='Maxent, raw values, Alfaro cultratus', bty="o")
scalebar(1500, xy = c(0,-500000), type = 'bar', divs = 2, below = "km")

if (file.exists(jar)) {
	## SSSmax was the method that varies the resampling of pseudo-absence (points) less, 
	## yielding a more constant threshold (Liu et al.,2013)
	e <- evaluate(p=testp, a=testa, paleo.LGM.CA)
	threshold.LGM.AC <- e@t[which.max(e@TPR + e@TNR)]
	pre_abs.LGM.AC <- px.LGM.AC > threshold.LGM.AC
	plot(pre_abs.LGM.AC, main="Presence/absence-SSSmax, Alfaro cultratus")
	plot(CA, add=TRUE, border='dark grey')
	points(pres_train, pch='+')} else { plot(1)
}

##--Plot the ROC curve: 
threshold(e)
plot(e, 'ROC', main= "ROC plot", sub= "Alfaro cultratus", bty="o")


## LIG ##
px.LIG.AC = predict(mx.AC, paleo.LIG.CA, overwrite=TRUE, type= "prop", ext=ext2, progress='', filename="C:/Users/leonardo/Desktop/Maxent/modelresult/Alfaro cultratus/LIG")
plot(px.LIG.AC, main='Maxent, raw values, Alfaro cultratus', bty="o")


if (file.exists(jar)) {
	## SSSmax was the method that varies the resampling of pseudo-absence (points) less, 
	## yielding a more constant threshold (Liu et al.,2013)
	threshold.LIG.AC <- e@t[which.max(e@TPR + e@TNR)]
	pre_abs.LIG.AC <- px.LIG.AC > threshold.LIG.AC
	plot(pre_abs.LIG.AC, main="Presence/absence-SSSmax, Alfaro cultratus")
	plot(CA, add=TRUE, border='dark grey')
	points(pres_train, pch='+')} else { plot(1)
}

######
# Models
######
models.AC <- stack(px.AC,px.LGM.AC,px.LIG.AC)
plot(models.AC)

######
# Stability areas
######
##--Superimpose the three models from presence-absence data, to create plots of stability
##--areas for Alfaro cultratus, across the three time slices modeled (present, LGM, LIG).
presence.AC <- stack(pre_abs.AC, pre_abs.LGM.AC, pre_abs.LIG.AC)
stability.AC <- sum(presence.AC)
plot(stability.AC)


##########################################################################################
#                          ***  Amatitlania spp. ANALYSIS  ***                           #
##########################################################################################
### I. Load species occurrences input data file:
sp2 <- read.table("Amatitlania spp.txt", h=T)
head(sp2)
summary(sp2)
tail(sp2)
coordinates(sp2) <- ~lon+lat
proj4string(sp2) <- CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs")
sp2 <- spTransform(sp2, CRS=crs.laea)
class(sp2)
head(sp2)
plot(sp2)

### II. Clip the layers and plot the occurrence points for this species:
##--Geographically resize/cut the map, where o = east, l = west, s = south, n = north; shape =
##--map that will be utilized.
plot(CA2,xlim=c(-1500000,500000),ylim=c(500000,3000000), main= "", bty="o", axes=TRUE, col="lightgreen")
title (main = "Presence points for Amatitlania spp.")

##--Plot the observed points:
points(sp2$lon, sp2$lat, col='blue', pch=20, cex=0.75)

##--If you want, you can highlight the points by making them red in color:
points(sp2$lon, sp2$lat, col='red', cex=0.75)

### III. Prepare testing and training data for MaxEnt:
##--Here, we divide the data points into testing and training data, based on Heuberty (1994):
group <- kfold(sp2, k=5)
pres_train <- sp2[group != 1, ]		## Training data
pres_test <- sp2[group == 1, ]		## Testing data

## Como ja feito anteriormente, para definir o tamanho da area de predicao
##para melhora a velocidade de processamento ext = extent(W,E,S,N)coordenadas dos vertices do retangulo

###  IV. RUN MAXENT. ### 
system.file("java", package="dismo")
jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='')
if (file.exists(jar)) {
	mx.AM <- maxent(bioclim.vars.CA, pres_train, overwrite=TRUE, args=c("-J","-r"), path= "C:/Users/leonardo/Desktop/Maxent/modelresult/Amatitlania spp/present")
	## Here, -J calls Jackknife method, -r asks if you want to 'overwrite' previous/existing results files
	plot(mx.AM, xlab="Percentage of contribution", ylab="Abiotic variable")
}

##--Background data points:
bg <- randomPoints(bioclim.vars.CA, n=10000, ext=CA2, extf = 1.25)

if (file.exists(jar)) {
	pvtest <- data.frame(extract(bioclim.vars.CA, pres_test))
	avtest <- data.frame(extract(bioclim.vars.CA, bg))
	testp <- predict(mx.AM, pvtest) 
	testa <- predict(mx.AM, avtest) 
	e <- evaluate(p=testp, a=testa, bioclim.vars.CA)
	e
	px.AM = predict(mx.AM, bioclim.vars.CA, ext=ext2, progress='', overwrite=TRUE, 
	filename="C:/Users/leonardo/Desktop/Maxent/modelresult/Amatitlania spp/present")
#
	plot(px.AM, main='Maxent, raw values,Amatitlania spp', bty="o")
}

if (file.exists(jar)) {
	## SSSmax was the method that varies the resampling of pseudo-absence (points) less, 
	## yielding a more constant threshold (Liu et al.,2013)
	threshold.AM <- e@t[which.max(e@TPR + e@TNR)]
	pre_abs.AM <- px.AM > threshold.AM
	plot(pre_abs.AM, main="Presence/absence-SSSmax, Amatitlania spp.")
	plot(CA2, add=TRUE, border='dark grey')
	points(pres_train, pch='+')} else { plot(1)
}

# para plotar ROC 
threshold(e)
plot(e, 'ROC', main= "ROC plot", sub= "Amatitlania spp",bty="o")

LGM
px.LGM.AM <- predict(mx.AM, paleo.LGM.CA, overwrite=TRUE, type= "prop", ext=ext2, progress='', filename="C:/Users/leonardo/Desktop/Maxent/modelresult/Amatitlania spp/LGM")
plot(px.LGM.AM, main='Maxent, raw values,Amatitlania spp', bty="o")

if (file.exists(jar)) {
	## SSSmax was the method that varies the resampling of pseudo-absence (points) less, 
	## yielding a more constant threshold (Liu et al.,2013)
	threshold.LGM.AM <- e@t[which.max(e@TPR + e@TNR)]
	pre_abs.LGM.AM <- px.LGM.AM > threshold.LGM.AM
	plot(pre_abs.LGM.AM, main="Presence/absence-SSSmax, Amatitlania spp.")
	plot(CA, add=TRUE, border='dark grey')
	points(pres_train, pch='+')} else { plot(1)
}

LIG
px.LIG.AM = predict(mx.AM, paleo.LIG.CA, overwrite=TRUE, ext=ext2, progress='', filename="C:/Users/leonardo/Desktop/Maxent/modelresult/Amatitlania spp/LIG")
plot(px.LIG.AM, main='Maxent, raw values,Amatitlania spp', bty="o")


if (file.exists(jar)) {
	## SSSmax was the method that varies the resampling of pseudo-absence (points) less, 
	## yielding a more constant threshold (Liu et al.,2013)
	threshold.LIG.AM <- e@t[which.max(e@TPR + e@TNR)]
	pre_abs.LIG.AM <- px.LIG.AM > threshold.LIG.AM
	plot(pre_abs.LIG.AM, main="Presence/absence-SSSmax, Amatitlania spp.")
	plot(CA, add=TRUE, border='dark grey')
	points(pres_train, pch='+')} else { plot(1)
}

######
# Models
######
models.AM <- stack(px.AM, px.LGM.AM, px.LIG.AM)
plot(models.AM)

######
# Stability areas
######
##--Superimpose the three models from presence-absence data, to create plots of stability
##--areas for Amatitlania spp., across the three time slices modeled (present, LGM, LIG).
presence.AM <- stack(pre_abs.AM, pre_abs.LGM.AM, pre_abs.LIG.AM)
stability.AM <- sum(presence.AM)
plot(stability.AM)


##########################################################################################
#                            ***  Astyanax spp. ANALYSIS  ***                            #
##########################################################################################
### I. Load species occurrences input data file:
sp3 <- read.table("Astyanax spp.txt", h=T)
head(sp3)
summary(sp3)
tail(sp3)
coordinates(sp3) <- ~lon+lat
proj4string(sp3) <- CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs")
sp3 <- spTransform(sp3, CRS=crs.laea)
class(sp3)
head(sp3)
plot(sp3)

### II. Clip the layers and plot the occurrence points for this species:
##--Geographically resize/cut the map, where o = east, l = west, s = south, n = north; shape =
##--map that will be utilized.
plot(CA2,xlim=c(-1500000,500000),ylim=c(500000,3000000), main= "", bty="o", axes=TRUE, col="lightgreen")
title (main = "Presence points for Astyanax spp.")

##--Plot the observed points:
points(sp3$lon, sp3$lat, col='blue', pch=20, cex=0.75)

##--If you want, you can highlight the points by making them red in color:
points(sp3$lon, sp3$lat, col='red', cex=0.75)

### III. Prepare testing and training data for MaxEnt:
##--Here, we divide the data points into testing and training data, based on Heuberty (1994):
group <- kfold(sp3, k=5)
pres_train <- sp3[group != 1, ]		## Training data
pres_test <- sp3[group == 1, ]		## Testing data

## Como ja feito anteriormente, para definir o tamanho da area de predicao
##para melhora a velocidade de processamento ext = extent(W,E,S,N)coordenadas dos vertices do retangulo

###  IV. RUN MAXENT. ### 
system.file("java", package="dismo")
jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='')
if (file.exists(jar)) {
	mx.AS <- maxent(bioclim.vars.CA, pres_train, overwrite=TRUE, args=c("-J","-r"), path= "C:/Users/leonardo/Desktop/Maxent/modelresult/Astyanax spp/present")
	## Here, -J calls Jackknife method, -r asks if you want to 'overwrite' previous/existing results files
	plot(mx.AS, xlab="Percentage of contribution", ylab="Abiotic variable")
}

##--Background data points:
bg <- randomPoints(bioclim.vars.CA, n=10000, ext=CA2, extf = 1.25)

if (file.exists(jar)) {
	pvtest <- data.frame(extract(bioclim.vars.CA, pres_test))
	avtest <- data.frame(extract(bioclim.vars.CA, bg))
	testp <- predict(mx.AS, pvtest) 
	testa <- predict(mx.AS, avtest) 
	e <- evaluate(p=testp, a=testa, bioclim.vars.CA)
	e
	px.AS = predict(mx.AS, bioclim.vars.CA, ext=ext2, overwrite=TRUE, progress='', 
	filename="C:/Users/leonardo/Desktop/Maxent/modelresult/Astyanax spp/present")
#
	plot(px.AS, main='Maxent, raw values,Astyanax spp', bty="o")
}

if (file.exists(jar)) {
	## SSSmax was the method that varies the resampling of pseudo-absence (points) less, 
	## yielding a more constant threshold (Liu et al.,2013)
	threshold.AS <- e@t[which.max(e@TPR + e@TNR)]
	pre_abs.AS <- px.AS > threshold.AS
	plot(pre_abs.AS, main="Presence/absence-SSSmax, Astyanax spp")
	plot(CA, add=TRUE, border='dark grey')
	points(pres_train, pch='+')} else { plot(1)
}

# para plotar ROC 
threshold(e)
plot(e, 'ROC', main= "ROC plot", sub= "Astyanax spp",bty="o")

LGM
px.LGM.AS <- predict(mx.AS, paleo.LGM.CA, overwrite=TRUE, type= "prop", ext=ext2, progress='', filename="C:/Users/leonardo/Desktop/Maxent/modelresult/Astyanax spp/LGM")
plot(px.LGM.AS, main='Maxent, raw values,Astyanax spp', bty="o")

if (file.exists(jar)) {
	## SSSmax was the method that varies the resampling of pseudo-absence (points) less, 
	## yielding a more constant threshold (Liu et al.,2013)
	threshold.LGM.AS <- e@t[which.max(e@TPR + e@TNR)]
	pre_abs.LGM.AS <- px.LGM.AS > threshold.LGM.AS
	plot(pre_abs.LGM.AS, main="Presence/absence-SSSmax, Astyanax spp")
	plot(CA, add=TRUE, border='dark grey')
	points(pres_train, pch='+')} else { plot(1)
}

LIG
px.LIG.AS = predict(mx.AS, paleo.LIG.CA, overwrite=TRUE, ext=ext2, progress='', filename="C:/Users/leonardo/Desktop/Maxent/modelresult/Astyanax spp/LIG")
plot(px.LIG.AS, main='Maxent, raw values,Astyanax spp', bty="o")


if (file.exists(jar)) {
	## SSSmax was the method that varies the resampling of pseudo-absence (points) less, 
	## yielding a more constant threshold (Liu et al.,2013)
	threshold.LIG.AS <- e@t[which.max(e@TPR + e@TNR)]
	pre_abs.LIG.AS <- px.LIG.AS > threshold.LIG.AS
	plot(pre_abs.LIG.AS, main="Presence/absence-SSSmax, Astyanax spp")
	plot(CA, add=TRUE, border='dark grey')
	points(pres_train, pch='+')} else { plot(1)
}

######
# Models
######
models.AS <- stack(px.AS,px.LGM.AS,px.LIG.AS)
plot(models.AS)

######
# Stability areas
######
##--Superimpose the three models from presence-absence data, to create plots of stability
##--areas for Amatitlania spp., across the three time slices modeled (present, LGM, LIG).
presence.AS <- stack(pre_abs.AS, pre_abs.LGM.AS, pre_abs.LIG.AS)
stability.AS <- sum(presence.AS)
plot(stability.AS)


##########################################################################################
#                         ***  Phallichthys amates ANALYSIS  ***                         #
##########################################################################################
### I. Load species occurrences input data file:
sp4 <- read.table("Phallichthys amates.txt", h=T)
head(sp4)
summary(sp4)
tail(sp4)
coordinates(sp4) <- ~lon+lat
proj4string(sp4) <- CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs")
sp4 <- spTransform(sp4, CRS=crs.laea)
class(sp4)
head(sp4)
plot(sp4)

### II. Clip the layers and plot the occurrence points for this species:
##--Geographically resize/cut the map, where o = east, l = west, s = south, n = north; shape =
##--map that will be utilized.
plot(CA2,xlim=c(-1500000,500000),ylim=c(500000,3000000), main= "", bty="o", axes=TRUE, col="lightgreen")
title (main = "Presence points for Phallichthys amates")

##--Plot the observed points:
points(sp4$lon, sp4$lat, col='blue', pch=20, cex=0.75)

##--If you want, you can highlight the points by making them red in color:
points(sp4$lon, sp4$lat, col='red', cex=0.75)

### III. Prepare testing and training data for MaxEnt:
##--Here, we divide the data points into testing and training data, based on Heuberty (1994):
group <- kfold(sp4, k=5)
pres_train <- sp4[group != 1, ]		## Training data
pres_test <- sp4[group == 1, ]		## Testing data

## Como ja feito anteriormente, para definir o tamanho da area de predicao
##para melhora a velocidade de processamento ext = extent(W,E,S,N)coordenadas dos vertices do retangulo

###  IV. RUN MAXENT. ### 
system.file("java", package="dismo")
jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='')
if (file.exists(jar)) {
	mx.PA <- maxent(bioclim.vars.CA, pres_train, overwrite=TRUE, args=c("-J","-r"), path= "C:/Users/leonardo/Desktop/Maxent/modelresult/Phallichthys amates/present")
	## Here, -J calls Jackknife method, -r asks if you want to 'overwrite' previous/existing results files
	plot(mx.PA, xlab="Percentage of contribution", ylab="Abiotic variable")
}

##--Background data points:
bg <- randomPoints(bioclim.vars.CA, n=10000, ext=CA2, extf = 1.25)

if (file.exists(jar)) {
	pvtest <- data.frame(extract(bioclim.vars.CA, pres_test))
	avtest <- data.frame(extract(bioclim.vars.CA, bg))
	testp <- predict(mx.PA, pvtest) 
	testa <- predict(mx.PA, avtest) 
	e <- evaluate(p=testp, a=testa, bioclim.vars.CA)
	e
	px.PA = predict(mx.PA, bioclim.vars.CA, ext=ext2, overwrite=TRUE, progress='',  
	filename="C:/Users/leonardo/Desktop/Maxent/modelresult/Phallichthys amates/present")
#
	plot(px.PA, main='Maxent, raw values,Phallichthys amates', bty="o")
}

if (file.exists(jar)) {
	## SSSmax was the method that varies the resampling of pseudo-absence (points) less, 
	## yielding a more constant threshold (Liu et al.,2013)
	threshold.PA <- e@t[which.max(e@TPR + e@TNR)]
	pre_abs.PA <- px.PA > threshold.PA
	plot(pre_abs.PA, main="Presence/absence-SSSmax, Phallichthys amates")
	plot(CA, add=TRUE, border='dark grey')
	points(pres_train, pch='+')} else { plot(1)
}

# para plotar ROC 
threshold(e)
plot(e, 'ROC', main= "ROC plot", sub= "Phallichthys amates",bty="o")

LGM
px.LGM.PA <- predict(mx.PA, paleo.LGM.CA, overwrite=TRUE, type= "prop", ext=ext2, progress='', filename="C:/Users/leonardo/Desktop/Maxent/modelresult/Phallichthys amates/LGM")
plot(px.LGM.PA, main='Maxent, raw values,Phallichthys amates', bty="o")

if (file.exists(jar)) {
	## SSSmax was the method that varies the resampling of pseudo-absence (points) less, 
	## yielding a more constant threshold (Liu et al.,2013)
	threshold.LGM.PA <- e@t[which.max(e@TPR + e@TNR)]
	pre_abs.LGM.PA <- px.LGM.PA > threshold.LGM.PA
	plot(pre_abs.LGM.PA, main="Presence/absence-SSSmax, Phallichthys amates")
	plot(CA, add=TRUE, border='dark grey')
	points(pres_train, pch='+')} else { plot(1)
}

LIG
px.LIG.PA = predict(mx.PA, paleo.LIG.CA, overwrite=TRUE, ext=ext2, progress='', filename="C:/Users/leonardo/Desktop/Maxent/modelresult/Phallichthys amates/LIG")
plot(px.LIG.PA, main='Maxent, raw values,Phallichthys amates', bty="o")
scalebar(500000, xy = c(-1500000,950000), type = 'bar', divs = 2, below = "km")


if (file.exists(jar)) {
	## SSSmax was the method that varies the resampling of pseudo-absence (points) less, 
	## yielding a more constant threshold (Liu et al.,2013)
	threshold.LIG.PA <- e@t[which.max(e@TPR + e@TNR)]
	pre_abs.LIG.PA <- px.LIG.PA > threshold.LIG.PA
	plot(pre_abs.LIG.PA, main="Presence/absence-SSSmax, Phallichthys amates")
	plot(CA, add=TRUE, border='dark grey')
	points(pres_train, pch='+')} else { plot(1)
}

######
# Models
######
models.PA <- stack(px.PA, px.LGM.PA, px.LIG.PA)
plot(models.PA)
scalebar(500000, xy = c(-1500000,950000), type = 'bar', divs = 2, below = "km")

######
# Stability areas
######
##--Superimpose the three models from presence-absence data, to create plots of stability
##--areas for Amatitlania spp., across the three time slices modeled (present, LGM, LIG).
presence.PA <- stack(pre_abs.PA, pre_abs.LGM.PA, pre_abs.LIG.PA)
stability.PA <- sum(presence.PA)
plot(stability.PA)


##########################################################################################
#                          ***  Poecilia mexicana ANALYSIS  ***                          #
##########################################################################################
### I. Load species occurrences input data file:
sp5 <- read.table("Poecilia mexicana.txt", h=T)
head(sp5)
summary(sp5)
tail(sp5)
coordinates(sp5) <- ~lon+lat
proj4string(sp5) <- CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs")
sp5 <- spTransform(sp5, CRS=crs.laea)
class(sp5)
head(sp5)
plot(sp5)

### II. Clip the layers and plot the occurrence points for this species:
##--Geographically resize/cut the map, where o = east, l = west, s = south, n = north; shape =
##--map that will be utilized.
plot(CA2,xlim=c(-1500000,500000),ylim=c(500000,3000000), main= "", bty="o", axes=TRUE, col="lightgreen")
title (main = "Presence points for Poecilia mexicana")

##--Plot the observed points:
points(sp5$lon, sp5$lat, col='blue', pch=20, cex=0.75)

##--If you want, you can highlight the points by making them red in color:
points(sp5$lon, sp5$lat, col='red', cex=0.75)

### III. Prepare testing and training data for MaxEnt:
##--Here, we divide the data points into testing and training data, based on Heuberty (1994):
group <- kfold(sp5, k=5)
pres_train <- sp5[group != 1, ]		## Training data
pres_test <- sp5[group == 1, ]		## Testing data

## Como ja feito anteriormente, para definir o tamanho da area de predicao
##para melhora a velocidade de processamento ext = extent(W,E,S,N)coordenadas dos vertices do retangulo

###  IV. RUN MAXENT. ### 
system.file("java", package="dismo")
jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='')
if (file.exists(jar)) {
	mx.PM <- maxent(bioclim.vars.CA, pres_train, overwrite=TRUE, args=c("-J","-r"), path= "C:/Users/leonardo/Desktop/Maxent/modelresult/Poecilia mexicana/present")
	## Here, -J calls Jackknife method, -r asks if you want to 'overwrite' previous/existing results files
	plot(mx.PM, xlab="Percentage of contribution", ylab="Abiotic variable")
}

##--Background data points:
bg <- randomPoints(bioclim.vars.CA, n=10000, ext=CA2, extf = 1.25)

if (file.exists(jar)) {
	pvtest <- data.frame(extract(bioclim.vars.CA, pres_test))
	avtest <- data.frame(extract(bioclim.vars.CA, bg))
	testp <- predict(mx.PM, pvtest) 
	testa <- predict(mx.PM, avtest) 
	e <- evaluate(p=testp, a=testa, bioclim.vars.CA)
	e
	px.PM = predict(mx.PM, bioclim.vars.CA, ext=ext2, overwrite=TRUE, progress='',  
	filename="C:/Users/leonardo/Desktop/Maxent/modelresult/Poecilia mexicana/present")
#
	plot(px.PM, main='Maxent, raw values,Poecilia mexicana', bty="o")
}

if (file.exists(jar)) {
	## SSSmax was the method that varies the resampling of pseudo-absence (points) less, 
	## yielding a more constant threshold (Liu et al.,2013)
	threshold.PM <- e@t[which.max(e@TPR + e@TNR)]
	pre_abs.PM <- px.PM > threshold.PM
	plot(pre_abs.PM, main="Presence/absence-SSSmax, Poecilia mexicana")
	plot(CA, add=TRUE, border='dark grey')
	points(pres_train, pch='+')} else { plot(1)
}

# para plotar ROC 
threshold(e)
plot(e, 'ROC', main= "ROC plot", sub= "Poecilia mexicana",bty="o")

LGM
px.LGM.PM <- predict(mx.PM, paleo.LGM.CA, overwrite=TRUE, type= "prop", ext=ext2, progress='', filename="C:/Users/leonardo/Desktop/Maxent/modelresult/Poecilia mexicana/LGM")
plot(px.LGM.PM, main='Maxent, raw values,Poecilia mexicana', bty="o")

if (file.exists(jar)) {
	## SSSmax was the method that varies the resampling of pseudo-absence (points) less, 
	## yielding a more constant threshold (Liu et al.,2013)
	threshold.LGM.PM <- e@t[which.max(e@TPR + e@TNR)]
	pre_abs.LGM.PM <- px.LGM.PM > threshold.LGM.PM
	plot(pre_abs.LGM.PM, main="Presence/absence-SSSmax, Poecilia mexicana")
	plot(CA, add=TRUE, border='dark grey')
	points(pres_train, pch='+')} else { plot(1)
}

LIG
px.LIG.PM = predict(mx.PM, paleo.LIG.CA, overwrite=TRUE, ext=ext2, progress='', filename="C:/Users/leonardo/Desktop/Maxent/modelresult/Poecilia mexicana/LIG")
plot(px.LIG.PM, main='Maxent, raw values,Poecilia mexicana', bty="o")


if (file.exists(jar)) {
	## SSSmax was the method that varies the resampling of pseudo-absence (points) less, 
	## yielding a more constant threshold (Liu et al.,2013)
	threshold.LIG.PM <- e@t[which.max(e@TPR + e@TNR)]
	pre_abs.LIG.PM <- px.LIG.PM > threshold.LIG.PM
	plot(pre_abs.LIG.PM, main="Presence/absence-SSSmax, Poecilia mexicana")
	plot(CA, add=TRUE, border='dark grey')
	points(pres_train, pch='+')} else { plot(1)
}

######
# Models
######
models.PM <- stack(px.PM, px.LGM.PM, px.LIG.PM)
plot(models.PM)

######
# Stability areas
######
##--Superimpose the three models from presence-absence data, to create plots of stability
##--areas for Amatitlania spp., across the three time slices modeled (present, LGM, LIG).
pres.PM <- stack(pre_abs.PM, pre_abs.LGM.PM, pre_abs.LIG.PM)
estab.PM <- sum(pres.PM)
plot(estab.PM)


##########################################################################################
#                       ***  Priapichthys annectens ANALYSIS  ***                        #
##########################################################################################
### I. Load species occurrences input data file:
sp6<-read.table("Priapichthys annectens.txt", h=T)
head(sp6)
summary(sp6)
tail(sp6)
coordinates(sp6) <- ~lon+lat
proj4string(sp6) <- CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs")
sp6 <- spTransform(sp6, CRS=crs.laea)
class(sp6)
head(sp6)
plot(sp6)

### II. Clip the layers and plot the occurrence points for this species:
##--Geographically resize/cut the map, where o = east, l = west, s = south, n = north; shape =
##--map that will be utilized.
plot(CA2,xlim=c(-1500000,500000),ylim=c(500000,3000000), main= "", bty="o", axes=TRUE, col="lightgreen")
title (main = "Presence points for Priapichthys annectens")

##--Plot the observed points:
points(sp6$lon, sp6$lat, col='blue', pch=20, cex=0.75)

##--If you want, you can highlight the points by making them red in color:
points(sp6$lon, sp6$lat, col='red', cex=0.75)

### III. Prepare testing and training data for MaxEnt:
##--Here, we divide the data points into testing and training data, based on Heuberty (1994):
group <- kfold(sp6, k=5)
pres_train <- sp6[group != 1, ]		## Training data
pres_test <- sp6[group == 1, ]		## Testing data

## Como ja feito anteriormente, para definir o tamanho da area de predicao
##para melhora a velocidade de processamento ext = extent(W,E,S,N)coordenadas dos vertices do retangulo

###  IV. RUN MAXENT. ### 
system.file("java", package="dismo")
jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='')
if (file.exists(jar)) {
	mx.PR <- maxent(bioclim.vars.CA, pres_train, overwrite=TRUE, args=c("-J","-r"), path= "C:/Users/leonardo/Desktop/Maxent/modelresult/Priapichthys annectens/present")
	## Here, -J calls Jackknife method, -r asks if you want to 'overwrite' previous/existing results files
	plot(mx.PR, xlab="Percentage of contribution", ylab="Abiotic variable")
}

##--Background data points:
bg <- randomPoints(bioclim.vars.CA, n=10000, ext=CA2, extf = 1.25)

if (file.exists(jar)) {
	pvtest <- data.frame(extract(bioclim.vars.CA, pres_test))
	avtest <- data.frame(extract(bioclim.vars.CA, bg))
	testp <- predict(mx.PR, pvtest) 
	testa <- predict(mx.PR, avtest) 
	e <- evaluate(p=testp, a=testa, bioclim.vars.CA)
	e
	px.PR = predict(mx.PR, bioclim.vars.CA, ext=ext2, overwrite=TRUE, progress='',  
	filename="C:/Users/leonardo/Desktop/Maxent/modelresult/Priapichthys annectens/present")
#
	plot(px.PR, main='Maxent, raw values,Priapichthys annectens', bty="o")
}

if (file.exists(jar)) {
	## SSSmax was the method that varies the resampling of pseudo-absence (points) less, 
	## yielding a more constant threshold (Liu et al.,2013)
	threshold.PR <- e@t[which.max(e@TPR + e@TNR)]
	pre_abs.PR <- px.PR > threshold.PR
	plot(pre_abs.PR, main="Presence/absence-SSSmax, Priapichthys annectens")
	plot(CA, add=TRUE, border='dark grey')
	points(pres_train, pch='+')} else { plot(1)
}

# para plotar ROC 
threshold(e)
plot(e, 'ROC', main= "ROC plot", sub= "Priapichthys annectens",bty="o")

LGM
px.LGM.PR <- predict(mx.PR, paleo.LGM.CA, overwrite=TRUE, type= "prop", ext=ext2, progress='', filename="C:/Users/leonardo/Desktop/Maxent/modelresult/Priapichthys annectens/LGM")
plot(px.LGM.PR, main='Maxent, raw values,Priapichthys annectens', bty="o")

if (file.exists(jar)) {
	## SSSmax was the method that varies the resampling of pseudo-absence (points) less, 
	## yielding a more constant threshold (Liu et al.,2013)
	threshold.LGM.PR <- e@t[which.max(e@TPR + e@TNR)]
	pre_abs.LGM.PR <- px.LGM.PR > threshold.LGM.PR
	plot(pre_abs.LGM.PR, main="Presence/absence-SSSmax, Priapichthys annectens")
	plot(CA, add=TRUE, border='dark grey')
	points(pres_train, pch='+')} else { plot(1)
}

LIG
px.LIG.PR = predict(mx.PR, paleo.LIG.CA, overwrite=TRUE, ext=ext2, progress='', filename="C:/Users/leonardo/Desktop/Maxent/modelresult/Priapichthys annectens/LIG")
plot(px.LIG.PR, main='Maxent, raw values,Priapichthys annectens', bty="o")


if (file.exists(jar)) {
	## SSSmax was the method that varies the resampling of pseudo-absence (points) less, 
	## yielding a more constant threshold (Liu et al.,2013)
	threshold.LIG.PR <- e@t[which.max(e@TPR + e@TNR)]
	pre_abs.LIG.PR <- px.LIG.PR > threshold.LIG.PR
	plot(pre_abs.LIG.PR, main="Presence/absence-SSSmax, Priapichthys annectens")
	plot(CA, add=TRUE, border='dark grey')
	points(pres_train, pch='+')} else { plot(1)
}

######
# Models
######
models.PR <- stack(px.PR, px.LGM.PR, px.LIG.PR)
plot(models.PR)

######
# Stability areas
######
##--Superimpose the three models from presence-absence data, to create plots of stability
##--areas for Amatitlania spp., across the three time slices modeled (present, LGM, LIG).
presence.PR <- stack(pre_abs.PR, pre_abs.LGM.PR, pre_abs.LIG.PR)
stability.PR <- sum(presence.PR)
plot(stability.PR)


##########################################################################################
#                        ***  Roeboides bouchellei ANALYSIS  ***                         #
##########################################################################################
### I. Load species occurrences input data file:
sp7<-read.table("Roeboides bouchellei.txt",h=T)
head(sp7)
summary(sp7)
tail(sp7)
coordinates(sp7) <- ~lon+lat
proj4string(sp7) <- CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs")
sp7 <- spTransform(sp7, CRS=crs.laea)
class(sp7)
head(sp7)
plot(sp7)

### II. Clip the layers and plot the occurrence points for this species:
##--Geographically resize/cut the map, where o = east, l = west, s = south, n = north; shape =
##--map that will be utilized.
plot(CA2,xlim=c(-1500000,500000),ylim=c(500000,3000000), main= "", bty="o", axes=TRUE, col="lightgreen")
title (main = "Presence points for Roeboides bouchellei")

##--Plot the observed points:
points(sp7$lon, sp7$lat, col='blue', pch=20, cex=0.75)

##--If you want, you can highlight the points by making them red in color:
points(sp7$lon, sp7$lat, col='red', cex=0.75)

### III. Prepare testing and training data for MaxEnt:
##--Here, we divide the data points into testing and training data, based on Heuberty (1994):
group <- kfold(sp7, k=5)
pres_train <- sp7[group != 1, ]		## Training data
pres_test <- sp7[group == 1, ]		## Testing data

## Como ja feito anteriormente, para definir o tamanho da area de predicao
##para melhora a velocidade de processamento ext = extent(W,E,S,N)coordenadas dos vertices do retangulo

###  IV. RUN MAXENT. ### 
system.file("java", package="dismo")
jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='')
if (file.exists(jar)) {
	mx.RB <- maxent(bioclim.vars.CA, pres_train, overwrite=TRUE, args=c("-J","-r"), path= "C:/Users/leonardo/Desktop/Maxent/modelresult/Roeboides bouchellei/present")
	## Here, -J calls Jackknife method, -r asks if you want to 'overwrite' previous/existing results files
	plot(mx, xlab="Percentage of contribution", ylab="Abiotic variable")
}

##--Background data points:
bg <- randomPoints(bioclim.vars.CA, n=10000, ext=CA2, extf = 1.25)

if (file.exists(jar)) {
	pvtest <- data.frame(extract(bioclim.vars.CA, pres_test))
	avtest <- data.frame(extract(bioclim.vars.CA, bg))
	testp <- predict(mx.RB, pvtest) 
	testa <- predict(mx.RB, avtest) 
	e <- evaluate(p=testp, a=testa, bioclim.vars.CA)
	e
	px.RB = predict(mx.RB, bioclim.vars.CA, ext=ext2, overwrite=TRUE, progress='',  
	filename="C:/Users/leonardo/Desktop/Maxent/modelresult/Roeboides bouchellei/present")
#
	plot(px.RB, main='Maxent, raw values,Roeboides bouchellei', bty="o")
}

if (file.exists(jar)) {
	## SSSmax was the method that varies the resampling of pseudo-absence (points) less, 
	## yielding a more constant threshold (Liu et al.,2013)
	threshold.RB <- e@t[which.max(e@TPR + e@TNR)]
	pre_abs.RB <- px.RB > threshold.RB
	plot(pre_abs.RB, main="Presence/absence-SSSmax, Roeboides bouchellei")
	plot(CA, add=TRUE, border='dark grey')
	points(pres_train, pch='+')} else { plot(1)
}

# para plotar ROC 
threshold(e)
plot(e, 'ROC', main= "ROC plot", sub= "Roeboides bouchellei",bty="o")

LGM
px.LGM.RB <- predict(mx.RB, paleo.LGM.CA, overwrite=TRUE, type= "prop", ext=ext2, progress='', filename="C:/Users/leonardo/Desktop/Maxent/modelresult/Roeboides bouchellei/LGM")
plot(px.LGM.RB, main='Maxent, raw values,Roeboides bouchellei', bty="o")

if (file.exists(jar)) {
	## SSSmax was the method that varies the resampling of pseudo-absence (points) less, 
	## yielding a more constant threshold (Liu et al.,2013)
	threshold.LGM.RB <- e@t[which.max(e@TPR + e@TNR)]
	pre_abs.LGM.RB <- px.LGM.RB > threshold.LGM.RB
	plot(pre_abs.LGM.RB, main="Presence/absence-SSSmax, Roeboides bouchellei")
	plot(CA, add=TRUE, border='dark grey')
	points(pres_train, pch='+')} else { plot(1)
}

LIG
px.LIG.RB = predict(mx.RB, paleo.LIG.CA, overwrite=TRUE, ext=ext2, progress='', filename="C:/Users/leonardo/Desktop/Maxent/modelresult/Roeboides bouchellei/LIG")
plot(px.LIG.RB, main='Maxent, raw values,Roeboides bouchellei', bty="o")


if (file.exists(jar)) {
	## SSSmax was the method that varies the resampling of pseudo-absence (points) less, 
	## yielding a more constant threshold (Liu et al.,2013)
	threshold.LIG.RB <- e@t[which.max(e@TPR + e@TNR)]
	pre_abs.LIG.RB <- px.LIG.RB > threshold.LIG.RB
	plot(pre_abs.LIG.RB, main="Presence/absence-SSSmax, Roeboides bouchellei")
	plot(CA, add=TRUE, border='dark grey')
	points(pres_train, pch='+')} else { plot(1)
}

######
# Models
######
models.RB <- stack(px.RB, px.LGM.RB, px.LIG.RB)
plot(models.RB)

######
# Stability areas
######
##--Superimpose the three models from presence-absence data, to create plots of stability
##--areas for Amatitlania spp., across the three time slices modeled (present, LGM, LIG).
presence.RB <- stack(pre_abs.RB, pre_abs.LGM.RB, pre_abs.LIG.RB)
stability.RB <- sum(presence.RB)
plot(stability.RB)

##########################################################################################
#                       ***  Xenophallus umbratilis ANALYSIS  ***                        #
##########################################################################################
### I. Load species occurrences input data file:
sp8<-read.table("Xenophallus umbratilis.txt", h=T)
head(sp8)
summary(sp8)
tail(sp8)
coordinates(sp8) <- ~lon+lat
proj4string(sp8) <- CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs")
sp8 <- spTransform(sp8, CRS=crs.laea)
class(sp8)
head(sp8)
plot(sp8)

### II. Clip the layers and plot the occurrence points for this species:
##--Geographically resize/cut the map, where o = east, l = west, s = south, n = north; shape =
##--map that will be utilized.
plot(CA2,xlim=c(-1500000,500000),ylim=c(500000,3000000), main= "", bty="o", axes=TRUE, col="lightgreen")
title (main = "Presence points for Xenophallus umbratilis")

##--Plot the observed data points:
points(sp8$lon, sp8$lat, col='blue', pch=20, cex=0.75)

##--If you want, you can highlight the points by making them red in color:
points(sp8$lon, sp8$lat, col='red', cex=0.75)

### III. Prepare testing and training data for MaxEnt:
##--Here, we divide the data points into testing and training data, based on Heuberty (1994):
group <- kfold(sp8, k=5)
pres_train <- sp8[group != 1, ]		## Training data
pres_test <- sp8[group == 1, ]		## Testing data

## Como ja feito anteriormente, para definir o tamanho da area de predicao
##para melhora a velocidade de processamento ext = extent(W,E,S,N)coordenadas dos vertices do retangulo

###  IV. RUN MAXENT. ### 
system.file("java", package="dismo")
jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='')
if (file.exists(jar)) {
	mx.XU <- maxent(bioclim.vars.CA, pres_train, overwrite=TRUE, args=c("-J","-r"), path= "C:/Users/leonardo/Desktop/Maxent/modelresult/Xenophallus umbratilis/present")
	## Here, -J calls Jackknife method, -r asks if you want to 'overwrite' previous/existing results files
	plot(mx.XU, xlab="Percentage of contribution", ylab="Abiotic variable")
}

##--Background data points:
bg <- randomPoints(bioclim.vars.CA, n=10000, ext=CA2, extf = 1.25)

if (file.exists(jar)) {
pvtest <- data.frame(extract(bioclim.vars.CA, pres_test))
avtest <- data.frame(extract(bioclim.vars.CA, bg))
testp <- predict(mx.XU, pvtest) 
testa <- predict(mx.XU, avtest) 
e <- evaluate(p=testp, a=testa, bioclim.vars.CA)
e
px.XU = predict(mx.XU, bioclim.vars.CA, overwrite=TRUE, ext=ext2, progress='',  filename="C:/Users/leonardo/Desktop/Maxent/modelresult/Xenophallus umbratilis/present")
plot(px.XU, main='Maxent, raw values,Xenophallus umbratilis', bty="o")
}

if (file.exists(jar)) {
# SSSmax foi o metodo que varia menos a reamostragem das pseudo-ausencia, sendo um threshold mais constante (Liu et al.,2013)
threshold.XU <- e@t[which.max(e@TPR + e@TNR)]
pre_abs.XU <- px.XU > threshold.XU
plot(pre_abs.XU, main="Presence/absence-SSSmax, Xenophallus umbratilis")
plot(CA, add=TRUE, border='dark grey')
points(pres_train, pch='+')
} else {
plot(1)
}

# para plotar ROC 
threshold(e)
plot(e, 'ROC', main= "ROC plot", sub= "Xenophallus umbratilis",bty="o")

LGM
px.LGM.XU <- predict(mx.XU, paleo.LGM.CA, overwrite=TRUE, type= "prop", ext=ext2, progress='', filename="C:/Users/leonardo/Desktop/Maxent/modelresult/Xenophallus umbratilis/LGM")
plot(px.LGM.XU, main='Maxent, raw values,Xenophallus umbratilis', bty="o")

if (file.exists(jar)) {
# SSSmax foi o metodo que varia menos a reamostragem das pseudo-ausencia, sendo um threshold mais constante (Liu et al.,2013)
e <- evaluate(p=testp, a=testa, paleo.LGM.CA)
threshold.LGM.XU <- e@t[which.max(e@TPR + e@TNR)]
pre_abs.LGM.XU <- px.LGM.XU > threshold.LGM.XU
plot(pre_abs.LGM.XU, main="Presence/absence-SSSmax, Xenophallus umbratilis")
plot(CA, add=TRUE, border='dark grey')
points(pres_train, pch='+')
} else {
plot(1)
}

LIG
px.LIG.XU = predict(mx.XU, paleo.LIG.CA, overwrite=TRUE, ext=ext2, progress='', filename="C:/Users/leonardo/Desktop/Maxent/modelresult/Xenophallus umbratilis/LIG")
plot(px.LIG.XU, main='Maxent, raw values,Xenophallus umbratilis', bty="o")


if (file.exists(jar)) {
# SSSmax foi o metodo que varia menos a reamostragem das pseudo-ausencia, sendo um threshold mais constante (Liu et al.,2013)
e <- evaluate(p=testp, a=testa, paleo.LIG.CA)
threshold.LIG.XU <- e@t[which.max(e@TPR + e@TNR)]
pre_abs.LIG.XU <- px.LIG.XU > threshold.LIG.XU
plot(pre_abs.LIG.XU, main="Presence/absence-SSSmax, Xenophallus umbratilis")
plot(CA, add=TRUE, border='dark grey')
points(pres_train, pch='+')
} else {
plot(1)
}

######
# Models
######
models.XU <- stack(px.XU, px.LGM.XU, px.LIG.XU)
plot(models.XU)

######
# Stability areas
######
##--Superimpose the three models from presence-absence data, to create plots of stability
##--areas for Amatitlania spp., across the three time slices modeled (present, LGM, LIG).
presence.XU <- stack(pre_abs.XU, pre_abs.LGM.XU, pre_abs.LIG.XU)
stability.XU <- sum(presence.XU)
plot(stability.XU)



#@@@@@@@@@@@@@@@@@@@@@@@@@@ COMBINED SPECIES/LINEAGE ANALYSES @@@@@@@@@@@@@@@@@@@@@@@@@@@#

######
#Todas as especies
#areas de estabilidade
######
#soprepor os 3 modelos de presenÁa e ausencia

modelos <- stack(estab.AC, estab.AM, estab.AS, estab.PA, estab.PM, estab.PR, estab.RB, estab.XU)
estab <- sum(modelos)
modelos1 <- stack(estab.AC, estab.AM, estab.AS)
modelos2 <- stack(estab.PA, estab.PM, estab.PR)
modelos3 <- stack(estab.RB, estab.XU, estab)
plot(modelos1)
plot(modelos2)
plot(modelos3)
plot(estab.AC)
plot(estab.AM)
plot(estab.AS)
plot(estab.PA)
plot(estab.PM)
plot(estab.PR)
plot(estab.RB)
plot(estab.XU)
plot(estab)



####################
Testes
####################

#####################################################

# carregar arquivo de dados
specie<-read.table("species.txt",h=T)
head(specie)
summary(specie)
tail(specie)
coordinates(specie) <- ~lon+lat
proj4string(specie) <- CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs")
specie <- spTransform(specie, CRS=crs.laea)
class(specie)
head(specie)
plot(specie)


#################################################

# Recorte geografico do mapa, no qual o= oeste, l=leste, s= sul, n= norte
# shape= mapa a ser utilizado
plot(CA2,xlim=c(-1500000,500000),ylim=c(500000,3000000), main= "", bty="o", axes=TRUE,col="lightgreen")
title (main = "Pontos de presenca para especie'x'")

# plotar pontos
points(specie$lon, specie$lat, col='blue', pch=20, cex=0.75)

# se quiser pode destacar ponto com contorno
points(specie$lon, specie$lat, col='red', cex=0.75)

######################################################
#*********************MODELAGEM**********************# 
######################################################

# Mantem-se apenas lon e lat
especie <- specie [,2:3]


# Para dividir em test e training data, baseado em Heuberty(1994)
# dividir em dados de treino e teste
group <- kfold(especie, k=5)
# treino
pres_train <- specie[group != 1, ]
# teste
pres_test <- specie[group == 1, ]

## Como ja feito anteriormente, para definir o tamanho da area de predicao
##para melhora a velocidade de processamento ext = extent(W,E,S,N)coordenadas dos vertices do retangulo

############################################################################
# para Maxent, com condicional que o programa esta corretamente instalado###
############################################################################

system.file("java", package="dismo")
jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='')
if (file.exists(jar)) {
mx <- maxent(bioclim.vars.CA, pres_train, args=c("-J","-r"), path= "C:/Users/leonardo/Desktop/Maxent/especie'x'")
# Nos quais "-J" indica o uso de Jackknife method, "r" para perguntar se deseja fazer "overwrite" dos dados, caso ja existam
plot(mx, xlab="Percentage of contribution", ylab="Abiotic variable")
}

# background data
bg <- randomPoints(bioclim.vars.CA, n=10000, ext=CA2, extf = 1.25)

if (file.exists(jar)) {
pvtest <- data.frame(extract(bioclim.vars.CA, pres_test))
avtest <- data.frame(extract(bioclim.vars.CA, bg))
testp <- predict(mx, pvtest) 
testa <- predict(mx, avtest) 
e <- evaluate(p=testp, a=testa, bioclim.vars.CA)
e
px = predict(mx, bioclim.vars.CA, ext=ext2, progress='',  filename="C:/Users/leonardo/Desktop/Maxent/especie'x'")
plot(px, main='Maxent, raw values,especie', bty="o")
}

if (file.exists(jar)) {
# SSSmax foi o metodo que varia menos a reamostragem das pseudo-ausencia, sendo um threshold mais constante (Liu et al.,2013)
threshold <- e@t[which.max(e@TPR + e@TNR)]
pre_abs <- px > threshold
plot(pre_abs, main="Presence/absence-SSSmax, especie")
plot(CA2, add=TRUE, border='dark grey')
points(pres_train, pch='+')
} else {
plot(1)
}

# para plotar ROC 
threshold(e)
plot(e, 'ROC', main= "ROC plot", sub= "especie'x'",bty="o")

LGM
if (file.exists(jar)) {
mx.LGM <- BIOMOD_Modeling(paleo.LGM.CA, models = c('MAXENT'),
pres_train, args=c("-J","-r"), path= "C:/Users/leonardo/Desktop/Maxent/especie'x'/LGM")

mx.LGM <- maxent(paleo.LGM.CA, pres_train, args=c("-J","-r"), path= "C:/Users/leonardo/Desktop/Maxent/especie'x'/LGM")
# Nos quais "-J" indica o uso de Jackknife method, "r" para perguntar se deseja fazer "overwrite" dos dados, caso ja existam
mx.LGM.proj <- BIOMOD_Projection(modeling.output = mx.LGM, new.env = bioclim.vars.CA,
proj.name = 'geral.LGM', selected.models = 'all', binary.meth = 'TSS', 
compress = FALSE, build.clamping.mask = TRUE)
plot(mx, xlab="Percentage of contribution", ylab="Abiotic variable")
}

# background data
bg <- randomPoints(paleo.LGM.CA, n=10000, ext=CA2, extf = 1.25)

if (file.exists(jar)) {
pvtest <- data.frame(extract(paleo.LGM.CA, pres_test))
avtest <- data.frame(extract(paleo.LGM.CA, bg))
testp <- predict(mx.LGM, pvtest) 
testa <- predict(mx.LGM, avtest) 
e <- evaluate(p=testp, a=testa, paleo.LGM.CA)
e
px.LGM = predict(mx.LGM, paleo.LGM.CA, verwrite=TRUE, ext=bioclim.vars.CA, progress='', filename="C:/Users/leonardo/Desktop/Maxent/especie'x'/LGM")
plot(px.LGM, main='Maxent, raw values,especie', bty="o")
}

if (file.exists(jar)) {
# SSSmax foi o metodo que varia menos a reamostragem das pseudo-ausencia, sendo um threshold mais constante (Liu et al.,2013)
threshold.LGM <- e@t[which.max(e@TPR + e@TNR)]
pre_abs.LGM <- px.LGM > threshold.LGM
plot(pre_abs.LGM, main="Presence/absence-SSSmax, especie'x'")
plot(CA, add=TRUE, border='dark grey')
points(pres_train, pch='+')
} else {
plot(1)
}

# para plotar ROC 
threshold(e)
plot(e, 'ROC', main= "ROC plot", sub= "especie'x'",bty="o")

LIG
if (file.exists(jar)) {
mx.LIG <- maxent(paleo.LIG.CA, pres_train, args=c("-J","-r"), path= "C:/Users/leonardo/Desktop/Maxent/especie'x'")
# Nos quais "-J" indica o uso de Jackknife method, "r" para perguntar se deseja fazer "overwrite" dos dados, caso ja existam
plot(mx, xlab="Percentage of contribution", ylab="Abiotic variable")
}

# background data
bg <- randomPoints(paleo.LIG.CA, n=10000, ext=CA, extf = 1.25)

if (file.exists(jar)) {
pvtest <- data.frame(extract(paleo.LIG.CA, pres_test))
avtest <- data.frame(extract(paleo.LIG.CA, bg))
testp <- predict(mx.LIG, pvtest) 
testa <- predict(mx.LIG, avtest) 
e <- evaluate(p=testp, a=testa, paleo.LIG.CA)
e
px.LIG = predict(mx.LIG, paleo.LIG.CA, overwrite=TRUE, ext=bioclim.vars.CA, progress='', filename="C:/Users/leonardo/Desktop/Maxent/especie'x'/LIG")
plot(px.LIG, main='Maxent, raw values,especie', bty="o")
}

if (file.exists(jar)) {
# SSSmax foi o metodo que varia menos a reamostragem das pseudo-ausencia, sendo um threshold mais constante (Liu et al.,2013)
threshold.LIG <- e@t[which.max(e@TPR + e@TNR)]
pre_abs.LIG <- px.LIG > threshold.LIG
plot(pre_abs.LIG, main="Presence/absence-SSSmax, especie'x'")
plot(CA, add=TRUE, border='dark grey')
points(pres_train, pch='+')
} else {
plot(1)
}

# para plotar ROC 
threshold(e)
plot(e, 'ROC', main= "ROC plot", sub= "especie'x'",bty="o")

######################################### END ############################################
