library(googlesheets)
library(raster)
library(rasterVis)
library(spdep)
library(dplyr)
library(maptools)
library(agricolae)
##################################
#Lectura de datos con googlesheets
##################################

gs_auth()
my_sheets <- gs_ls()
gap <- gs_title("Roya_Data")
roya <- gap %>%
  gs_read(ws = "Results")
roya <- as.data.frame(roya)
rm(my_sheets,gap)
roya$S1 <-as.factor(roya$S1)
roya$S2 <-as.factor(roya$S2)
roya$S3 <-as.factor(roya$S3)
roya$DS <-as.factor(roya$DS)
roya$PS <-as.factor(roya$PS)
roya$PC <-as.factor(roya$PC)
roya$FU <-as.factor(roya$FU)
roya$DH <-as.factor(roya$DH)
roya$FE <-as.factor(roya$FE)


roya$HE2[roya$HE2 == 0] <- 1
roya$HT2[roya$HT2 == 0] <- 1
#########################
#Regresion Poisson simple
#########################
md1 <- glm( formula =HE2  ~ Z, family = poisson(link = log),data=roya,na.action = na.omit)
summary(md1)
source('http://www.poleto.com/funcoes/envel.pois.txt')
source('http://www.poleto.com/funcoes/diag.pois.txt')
envel.pois(md1)
diag.pois(md1)

phi<-md1$deviance/md1$df.residual
phi

#Regresion Poisson con offset
md2 <- glm(formula = HE2 ~ Z +offset(log(HT2)), family = poisson(link = log),data=roya,na.action = na.omit)
summary(md2)
envel.pois(md2)
diag.pois(md2,iden=c(2,2,2,0,3,3,3,0))

phi<-md2$deviance/md2$df.residual
phi

#Regresion binomial negativa
library(MASS)
md3=glm.nb(HE2 ~ Z,link=log,data=roya)
summary(md3)
source('http://www.poleto.com/funcoes/diag.nb.txt')
diag.nb(md3)
source('http://www.poleto.com/funcoes/envel.nb.txt')
envel.nb(md3)


#####################
## Analisis exploratorio espacial
#####################
library(raster)
library(sp)
library(rgdal)
#Asignar coordenadas y residuos
roya.sp <- roya
roya.sp[,"Residuos"] <- residuals(md2)
coordinates(roya.sp) <- ~X+Y
crs(roya.sp) <- "+proj=utm +zone=18 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

#Verificar la autocorrelacion espacial
#1. matriz de pesos
sw <- knearneigh(roya.sp, longlat = FALSE)
#2. Crear objeto nb
sw <- knn2nb(sw,sym = TRUE)
is.symmetric.nb(sw)
#3. Crear lista de vecinos
sw <- nb2listw(sw,zero.policy = TRUE)
#4. Comprobar autocorrelacion espacial
#Ho: El coeficiennte es 0 - no hay autocorrelacion espacial
#Ha: El coeficiente es diferente de 0
moran.test(x = roya$HE2,listw = sw,randomisation=TRUE,adjust.n=TRUE)
#4.1 Comprobar autocorrelacion espacial "no para metrico"
moran.mc(x = roya$HE2,listw = sw,nsim=100000,adjust.n=TRUE)

#####################
## Lectura de covariables
#####################
dem <- raster(file.choose())
dem <-projectRaster(dem,crs="+proj=utm +zone=18 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
plot(dem); plot(roya.sp,add=TRUE)
aoi <- drawExtent()
dem <- crop(dem,aoi)

incidence <- predict(md2,dem)
