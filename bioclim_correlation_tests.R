
install.packages("ENMtools")
require(devtools)
install_github("dyerlab/gstudio")

install.packages(c("raster", "sp", "gdistance"), dependencies = TRUE)
install.packages("rgdal")
require(gstudio)
require(raster)
require(rgdal)
install.packages("dismo")
require(dismo)

install.packages ("ENMeval")

setwd("~/Desktop/modeling_materials_bio1-19_30s_bil")


#### correlation test ####

bilFile <- list.files (pattern=".bil")


for(file in bilFile){
  print(file)
  r<- raster(file, native=TRUE)
  x <- extract(r, Data_A[,c('Longitude', 'Latitude')])
  Data_A[[file]] <- x
}



setwd("~/Desktop/modeling_materials")
Data_Asty <- read.csv(file="72_Astyanax_genus_coords_for_Maxent_noRedundant.csv")
Data_A <- Data_Asty[,2:3]
setwd("~/Desktop/modeling_materials/cropped_bils_Asty")

Data_Xeno <- read.csv(file="Xu_58_coords.csv")
Data_X <- Data_Xeno[,2:3]
setwd("~/Desktop/modeling_materials/Xeno_bils")

Data_Pha <- read.csv(file="Pha_79_coords.csv")
Data_P <- Data_Pha[,2:3]
setwd("~/Desktop/modeling_materials/Ro_and_Pha_bils")
                        
clim <- list.files(pattern = ".bil",full.names = TRUE)    

Data_X
Data_P

for(file in clim){
  print(file)
  r<- raster(file, native=TRUE)
  x <- extract(r, Data_P[,c('Longitude', 'Latitude')])
  Data_P[[file]] <- x
}


install.packages("corrplot")
library(corrplot)
                        
Corplot1 <- cor(Data_P, use="pairwise.complete.obs")
Corplot1
setwd("~/Desktop/modeling_materials")
write.csv(Corplot1, file="Pha_bioclim_correlations.csv")



setwd("~/Desktop/modeling_materials/bio1-19_30s_bil")
ex <- extent(-87,-81, 8, 13) # Xenophallus
ex <- extent (-95,-76, 7, 17) # Pha_79_coords AND  Roehoides_bouchellei 

bio1 <- raster("bio_1.bil")
bio1c <- crop(x=bio1, y=ex)
plot(bio1c)
setwd("~/Desktop/modeling_materials/Xeno_bils")
setwd("~/Desktop/modeling_materials/Ro_and_Pha_bils")
writeRaster(x=bio1c, file="bio1_pha.bil", overwrite=TRUE)

bio2 <- raster("bio_2.bil")
bio2c <- crop(x=bio2, y=ex)
plot(bio2c)
writeRaster(x=bio2c, file="bio2_pha.bil", overwrite=TRUE)

bio3 <- raster("bio_3.bil")
bio3c <- crop(x=bio3, y=ex)
writeRaster(x=bio3c, file="bio3_pha.bil", overwrite=TRUE)

bio4 <- raster("bio_4.bil")
bio4c <- crop(x=bio4, y=ex)
writeRaster(x=bio4c, file="bio4_pha.bil", overwrite=TRUE)

bio5 <- raster("bio_5.bil")
bio5c <- crop(x=bio5, y=ex)
plot(bio5c)
writeRaster(x=bio5c, file="bio5_pha.bil", overwrite=TRUE)

bio6 <- raster("bio_6.bil")
bio6c <- crop(x=bio6, y=ex)
plot(bio6c)
writeRaster(x=bio6c, file="bio6_pha.bil", overwrite=TRUE)

bio7 <- raster("bio_7.bil")
bio7c <- crop(x=bio7, y=ex)
plot(bio7c)
writeRaster(x=bio7c, file="bio7_pha.bil", overwrite=TRUE)

bio8 <- raster("bio_8.bil")
bio8c <- crop(x=bio8, y=ex)
writeRaster(x=bio8c, file="bio8_pha.bil", overwrite=TRUE)

bio9 <- raster("bio_9.bil")
bio9c <- crop(x=bio9, y=ex)
writeRaster(x=bio9c, file="bio9_pha.bil", overwrite=TRUE)

bio10 <-raster("bio_10.bil")
bio10c <- crop(x=bio10, y=ex)
writeRaster(x=bio10c, file="bio10_pha.bil", overwrite=TRUE)

bio11 <-raster("bio_11.bil")
bio11c <- crop(x=bio11, y=ex)
writeRaster(x=bio11c, file="bio11_pha.bil", overwrite=TRUE)

bio12 <-raster("bio_12.bil")
bio12c <- crop(x=bio12, y=ex)
writeRaster(x=bio12c, file="bio12_pha.bil", overwrite=TRUE)

bio13 <-raster("bio_13.bil")
bio13c <- crop(x=bio13, y=ex)
writeRaster(x=bio13c, file="bio13_pha.bil", overwrite=TRUE)

bio14 <-raster("bio_14.bil")
bio14c <- crop(x=bio14, y=ex)
writeRaster(x=bio14c, file="bio14_pha.bil", overwrite=TRUE)

bio15 <-raster("bio_15.bil")
bio15c <- crop(x=bio15, y=ex)
writeRaster(x=bio15c, file="bio15_pha.bil", overwrite=TRUE)

bio16 <-raster("bio_16.bil")
bio16c <- crop(x=bio16, y=ex)
writeRaster(x=bio16c, file="bio16_pha.bil", overwrite=TRUE)

bio17 <-raster("bio_17.bil")
bio17c <- crop(x=bio17, y=ex)
writeRaster(x=bio17c, file="bio17_pha.bil", overwrite=TRUE)

bio18 <-raster("bio_18.bil")
bio18c <- crop(x=bio18, y=ex)
writeRaster(x=bio18c, file="bio18_pha.bil", overwrite=TRUE)

bio19 <-raster("bio_19.bil")
bio19c <- crop(x=bio19, y=ex)
writeRaster(x=bio19c, file="bio19_pha.bil", overwrite=TRUE)

setwd("~/Desktop/modeling_materials/alt_23")
alt23 <- raster("alt_23.bil")
alt23c <- crop(x=alt23, y=ex)
writeRaster(x=alt23c, file="alt23_pha.bil", overwrite=TRUE)


#### From bil to .grd ....check the extents ####
bio1 <- raster("bio_1.bil")
ex <- extent(-97,-77,6,25)
bio1c <- crop(x=bio1, y=ex)
plot(bio1c)
writeRaster(x=bio1c, file="bio1c.grd", overwrite=TRUE)

bio2 <- raster("bio_2.bil")
ex <- extent(-97,-77,6,25)
bio2c <- crop(x=bio2, y=ex)
plot(bio2c)
writeRaster(x=bio2c, file="bio2c.grd", overwrite=TRUE)

bio3 <- raster("bio_3.bil")
ex <- extent(-97,-77,6,25)
bio3c <- crop(x=bio3, y=ex)
plot(bio3c)
writeRaster(x=bio3c, file="bio3c.grd", overwrite=TRUE)

bio4 <- raster("bio_4.bil")
bio4c <- crop(x=bio4, y=ex)
writeRaster(x=bio4c, file="bio4c.grd", overwrite=TRUE)

bio5 <- raster("bio_5.bil")
ex <- extent(-97,-77,6,25)
bio5c <- crop(x=bio5, y=ex)
plot(bio5c)
writeRaster(x=bio5c, file="bio5c.grd", overwrite=TRUE)

bio6 <- raster("bio_6.bil")
ex <- extent(-97,-77,6,25)
bio6c <- crop(x=bio6, y=ex)
plot(bio6c)
writeRaster(x=bio6c, file="bio6c.grd", overwrite=TRUE)

bio7 <- raster("bio_7.bil")
ex <- extent(-97,-77,6,25)
bio7c <- crop(x=bio7, y=ex)
plot(bio7c)
writeRaster(x=bio7c, file="bio7c.grd", overwrite=TRUE)

bio8 <- raster("bio_8.bil")
ex <- extent(-97,-77,6,25)
bio8c <- crop(x=bio8, y=ex)
plot(bio8c)
writeRaster(x=bio8c, file="bio8c.grd", overwrite=TRUE)

bio9 <- raster("bio_9.bil")
bio9c <- crop(x=bio9, y=ex)
writeRaster(x=bio9c, file="bio9c.grd", overwrite=TRUE)

bio10 <-raster("bio_10.bil")
bio10c <- crop(x=bio10, y=ex)
writeRaster(x=bio10c, file="bio10c.grd", overwrite=TRUE)

bio11 <-raster("bio_11.bil")
bio11c <- crop(x=bio11, y=ex)
writeRaster(x=bio11c, file="bio11c.grd", overwrite=TRUE)

bio12 <-raster("bio_12.bil")
bio12c <- crop(x=bio12, y=ex)
writeRaster(x=bio12c, file="bio12c.grd", overwrite=TRUE)

bio13 <-raster("bio_13.bil")
bio13c <- crop(x=bio13, y=ex)
writeRaster(x=bio13c, file="bio13c.grd", overwrite=TRUE)

bio14 <-raster("bio_14.bil")
bio14c <- crop(x=bio14, y=ex)
writeRaster(x=bio14c, file="bio14c.grd", overwrite=TRUE)

bio15 <-raster("bio_15.bil")
bio15c <- crop(x=bio15, y=ex)
writeRaster(x=bio15c, file="bio15c.grd", overwrite=TRUE)

bio16 <-raster("bio_16.bil")
bio16c <- crop(x=bio16, y=ex)
writeRaster(x=bio16c, file="bio16c.grd", overwrite=TRUE)

bio17 <-raster("bio_17.bil")
bio17c <- crop(x=bio17, y=ex)
writeRaster(x=bio17c, file="bio17c.grd", overwrite=TRUE)

bio18 <-raster("bio_18.bil")
bio18c <- crop(x=bio18, y=ex)
writeRaster(x=bio18c, file="bio18c.grd", overwrite=TRUE)

bio19 <-raster("bio_19.bil")
bio19c <- crop(x=bio19, y=ex)
writeRaster(x=bio19c, file="bio19c.grd", overwrite=TRUE)





##### convert extent to UTM _ from Justin's script #####
crs.laea <- CRS("+proj=laea +lat_0=0 +lon_0=-80 +ellps=WGS84 +units=m +no_defs")

bioclim.vars <- list.files(path = paste(system.file(package = "dismo"),'/ex/wc0.5', sep = ''), pattern = 'grd', full.names = TRUE)  
## NOTE: You may need to change this from pattern='bil' to pattern='grd' to get it to work.
bioclim.vars2 <- stack(bioclim.vars)
bioclim.vars.crop <- crop(bioclim.vars, ext)
bioclim.vars.LAEA <- projectRaster(bioclim.vars.crop, crs=crs.laea)
## JCB: I'm commenting this out, because I think it should be as on newline below:  plot(bioclim.vars.LAEA$bio_1)
plot(bioclim.vars.LAEA[[1]], main=names(bioclim.vars.LAEA )[1])


#### convert to UTM_ my script #####

setwd("~/Desktop/bio1-19_30s_grd")
bio1 <- raster("bio1c.grd")
bio2 <- raster("bio2c.grd")
bio3 <- raster("bio3c.grd")
bio4 <- raster("bio4c.grd")
bio5 <- raster("bio5c.grd")
bio6 <- raster("bio6c.grd")
bio7 <- raster("bio7c.grd")
bio8 <- raster("bio8c.grd")
bio9 <- raster("bio9c.grd")
bio10 <- raster("bio10c.grd")
bio11 <- raster("bio11c.grd")
bio12 <- raster("bio12c.grd")
bio13 <- raster("bio13c.grd")
bio14 <- raster("bio14c.grd")
bio15 <- raster("bio15c.grd")
bio16 <- raster("bio16c.grd")
bio17 <- raster("bio17c.grd")
bio18 <- raster("bio18c.grd")
bio19<- raster("bio19c.grd")


bio1p <- projectRaster(bio1c, crs=crs.laea)
writeRaster(bio1p, "bio1_utm.bil")
plot(bio1p)
points(Data_A$Longitude, Data_A$Latitude, pch=15, col="dodgerblue4", cex=20, add=TRUE)

bio2p <- projectRaster(bio2c, crs=crs.laea)
writeRaster(bio2p, "bio2_utm.bil")

bio3p <- projectRaster(bio3c, crs=crs.laea)
writeRaster(bio3p, "bio3_utm.bil")

bio4p <- projectRaster(bio4c, crs=crs.laea)
writeRaster(bio4p, "bio4_utm.bil")
bio5p <- projectRaster(bio5c, crs=crs.laea)
writeRaster(bio5p, "bio5_utm.bil")
bio6p <- projectRaster(bio6c, crs=crs.laea)
writeRaster(bio6p, "bio6_utm.bil")
bio7p <- projectRaster(bio7c, crs=crs.laea)
writeRaster(bio7p, "bio7_utm.bil")
bio8p <- projectRaster(bio8c, crs=crs.laea)
writeRaster(bio8p, "bio8_utm.bil")
bio9p <- projectRaster(bio9c, crs=crs.laea)
writeRaster(bio9p, "bio9_utm.bil")
bio10p <- projectRaster(bio10c, crs=crs.laea)
writeRaster(bio10p, "bio10_utm.bil")
bio11p <- projectRaster(bio11c, crs=crs.laea)
writeRaster(bio11p, "bio11_utm.bil")
bio12p <- projectRaster(bio12c, crs=crs.laea)
writeRaster(bio12p, "bio12_utm.bil")
bio13p <- projectRaster(bio13c, crs=crs.laea)
writeRaster(bio13p, "bio13_utm.bil")
bio14p <- projectRaster(bio14c, crs=crs.laea)
writeRaster(bio14p, "bio14_utm.bil")
bio15p <- projectRaster(bio15c, crs=crs.laea)
writeRaster(bio15p, "bio15_utm.bil")
bio16p <- projectRaster(bio16c, crs=crs.laea)
writeRaster(bio16p, "bio16_utm.bil")
bio17p <- projectRaster(bio17c, crs=crs.laea)
writeRaster(bio17p, "bio17_utm.bil")
bio18p <- projectRaster(bio18c, crs=crs.laea)
writeRaster(bio18p, "bio18_utm.bil")
bio19p <- projectRaster(bio19c, crs=crs.laea)
writeRaster(bio19p, "bio19_utm.bil")


plot(bio18c)
plot(bio19c)
crs.laea <- CRS("+proj=laea +lat_0=0 +lon_0=-80 +ellps=WGS84 +units=m +no_defs")
bioStack <- stack(bio1, bio2, bio3, bio4, bio5, bio6, bio7, bio8, bio9, bio10, bio11, bio12, bio13, bio14, bio15, bio16, bio17, bio18, bio19)
bioclim.vars.LAEA <- projectRaster(bioStack, crs=crs.laea)


plot(bioclim.vars.LAEA[[1]], main=names(bioclim.vars.LAEA )[1])

#### didn't use this... ext2 <- extent(-1500000,300000,500000,25000000) ####

