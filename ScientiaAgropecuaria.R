####################################
########### LIBRERIAS ##############
####################################
library(raster)
library(cluster)
library(FactoMineR)
library(agricolae)
library(psych)
library(DMwR)
library(fpc)
library(StatMatch)
##LECTURA DE COVARIABLES

startdir <- getwd()
setwd("F:/roya/ScientiaAgropecuaria2020/SAGAGIS") 
files <- list.files(pattern="sdat$"); files
stack1 <- list()
for(i in 1:length(files)) {
  stack1[[i]] <- raster(files[i])}
covariables <- do.call(stack, stack1) ### JO!
setwd(startdir)
covariables <- dropLayer(covariables,i=1)
names(covariables)
rm(files,i,stack1,startdir)

##LECTURA DE DATOS DE ROYA
roya <- read.delim("clipboard",header = TRUE)
roya <- roya[,c(1,7,8,15,26)]

##Extraer covariables a datos
roya.sp <- roya
coordinates(roya.sp) <- ~X+Y
crs(roya.sp) <- "+proj=utm +zone=18 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
roya2 <-extract(covariables, roya.sp)
roya <- cbind(roya,roya2)
rm(roya2)

