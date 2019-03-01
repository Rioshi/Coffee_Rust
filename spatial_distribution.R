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
cof <- gap %>%
  gs_read(ws = "Results")
cof <- as.data.frame(cof)
rm(my_sheets,gap)
cof$S1 <-as.factor(cof$S1)
cof$S2 <-as.factor(cof$S2)
cof$S3 <-as.factor(cof$S3)
cof$DS <-as.factor(cof$DS)
cof$PS <-as.factor(cof$PS)
cof$PC <-as.factor(cof$PC)
cof$FU <-as.factor(cof$FU)
cof$DH <-as.factor(cof$DH)
cof$FE <-as.factor(cof$FE)

##################################
#Analisis exploratorio de datos
##################################

summary(cof)
boxplot(cof$Incidencia3)
hist(cof$Incidencia3)
