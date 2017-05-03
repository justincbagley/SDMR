#!/usr/bin/env Rscript

##########################################################################################
#                         MaxEntAnalyzer Rscript v1.0, May 2017                          #
#  Copyright (c)2017 Justinc C. Bagley, Virginia Commonwealth University, Richmond, VA,  #
#  USA; Universidade de Brasília, Brasília, DF, Brazil. See README and license on GitHub #
#  (https://github.com/justincbagley) for further information. Last update: May 3, 2017. #
#  For questions, please email jcbagley@vcu.edu.                                         #
#                                                                                        #
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

setwd("~/PATH/TO/ANALYSIS/FOLDER/")

#rm(list=ls())	## Not run. Uncomment this line to clear the workspace when this script is sourced.

##--Load needed library, R code, or package stuff. Install package if not present.
#source('MaxEntAnalyzer.R', chdir = TRUE)
packages <- c('sp', 'grid', 'lattice', 'ape', 'raster', 'maptools', 'phyloclim', 'dismo', 'rJava', 'rgdal', 'spocc', 'ENMeval')
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
    install.packages(setdiff(packages, rownames(installed.packages())))
}

library(sp); library(grid); library(lattice); library(maptools); library(phyloclim) # Loads dependencies ape and raster (but masks part of raster).
library(dismo); library(rgdal); library(rJava);
library(spocc); library(ENMeval);

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
## JCB: I'm commenting this out, because I think it should be as on newline below:  plot(bioclim.vars.LAEA$bio_1)
plot(bioclim.vars.LAEA[[1]], main=names(bioclim.vars.LAEA )[1])

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
plot(CA2, xlim=c(-1500000,500000), ylim=c(500000,3000000), main= "", bty="o", axes=TRUE, col="lightgreen")
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
	path= "~/PATH/TO/Alfaro\_cultratus/present/")
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
	filename="~/PATH/TO/Alfaro\_cultratus/present")
#
	plot(px.AC, main='Maxent, raw values, Alfaro cultratus', bty="o")
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
filename="~/PATH/TO/Alfaro\_cultratus/LGM")
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
px.LIG.AC = predict(mx.AC, paleo.LIG.CA, overwrite=TRUE, type= "prop", ext=ext2, progress='', filename="~/PATH/TO/Alfaro\_cultratus/LIG")
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






[ ADD MORE SPECIES MODELING BLOCKS HERE BY MODIFYING THE EXAMPLE BLOCK ABOVE ] 



#@@@@@@@@@@@@@@@@@@@@@@@@@@ COMBINED SPECIES/LINEAGE ANALYSES @@@@@@@@@@@@@@@@@@@@@@@@@@@#

##########################################################################################
#                     ***  Map stability areas for ALL species  ***                      #
##########################################################################################
##--Superimposes the 3 presence-absence models for each species and plots a stability map.
models <- stack(estab.AC, estab.AM, estab.AS, estab.PA, estab.PM, estab.PR, estab.RB, estab.XU)
stability <- sum(models)
row1_3mods <- stack(estab.AC, estab.AM, estab.AS)
row2_3mods <- stack(estab.PA, estab.PM, estab.PR)
row3_2mod_stab <- stack(estab.RB, estab.XU, stability)
plot(row1_3mods)
plot(row2_3mods)
plot(row3_2mod_stab)
#
plot(estab.AC)
plot(estab.AM)
plot(estab.AS)
plot(estab.PA)
plot(estab.PM)
plot(estab.PR)
plot(estab.RB)
plot(estab.XU)
plot(stability)


## NOTE: Add graphics.off(), pdf(), and dev.off() code to previous section saving graphical output to file.




######################################### END ############################################
