####################################
########### LIBRERIAS ##############
####################################
library(cluster)
library(FactoMineR)
library(factoextra)
library(agricolae)
library(psych)
library(DMwR)
library(fpc)
library(StatMatch)
library(betareg)
library(caret)
library(raster)
library(gstat)
library(sp)
library(rgdal)
library(spdep)
library(ggplot2)


########################
##LECTURA DE COVARIABLES
########################

#startdir <- getwd()
#setwd("F:/roya/ScientiaAgropecuaria2020/SAGAGIS") 
#files <- list.files(pattern="sdat$"); files
#stack1 <- list()
#for(i in 1:length(files)) {
#  stack1[[i]] <- raster(files[i])}
#covariables <- do.call(stack, stack1) ### JO!
#setwd(startdir)
#names(covariables)
#rm(files,i,stack1,startdir)
#saveRDS(covariables,file="F:/roya/ScientiaAgropecuaria2020/Rdata/covariables.rds")

covariables <- readRDS("F:/roya/ScientiaAgropecuaria2020/Rdata/covariables.rds")
covariables <- dropLayer(covariables,i=c(1))
nn <- names(covariables)
nn[c(1,2,6,7,8,17,18)] <- c("Asp","Elv","Ppt","LST_B10","LST_B11","Tm","Ws")
names(covariables) <- nn
rm(nn)

##########################
##LECTURA DE DATOS DE ROYA
##########################

#roya <- read.delim("clipboard",header = TRUE)
#saveRDS(roya,file="F:/roya/ScientiaAgropecuaria2020/Rdata/roya.rds")

roya <- readRDS("F:/roya/ScientiaAgropecuaria2020/Rdata/roya.rds")
roya[,"ITI"] <- roya$HE1/roya$HT1 #Incidencia tercio inferior
roya[,"ITM"] <- roya$HE2/roya$HT2 #Incidencia tercio medio
roya[,"ITS"] <- roya$HE3/roya$HT3 #Incidencia tercio superior
roya[,"IT"]  <- (roya$HE1+roya$HE2+roya$HE3)/(roya$HT1+roya$HT2+roya$HT3) #Incidencia total de la planta

roya$ITI <- replace(x=roya$ITI,is.na(roya$ITI),0)
roya$ITM <- replace(x=roya$ITM,is.na(roya$ITM),0)
roya$ITS <- replace(x=roya$ITS,is.na(roya$ITS),0)
roya$IT <- replace(x=roya$IT,is.na(roya$IT),0)

roya[,"ST"] <- apply(X = roya[,c("S1","S2","S3")],MARGIN = 1,FUN=max) #Severidad maxima de la planta

nn <- names(roya) ; nn <- nn[c(1,7,8,12,15,18,19:24,28:length(nn))]
roya <- roya[,nn] ; rm(nn)

roya$S1 <- ordered(roya$S1) ; roya$S2 <- ordered(roya$S2); roya$S3 <- ordered(roya$S3) ; 
roya$S2 <- ordered(roya$S2) ; roya$ST <- ordered(roya$ST)

##Transformar en base a (Smithson and Verkuilen 2006) para no tener 0 y 1
#roya$ITI <- replace(x=roya$ITI,roya$ITI==0,0.01) 
roya$ITI <- (roya$ITI*(length(roya$ITI)-1)+0.5)/length(roya$ITI)
roya$ITM <- (roya$ITM*(length(roya$ITM)-1)+0.5)/length(roya$ITM)
roya$ITS <- (roya$ITS*(length(roya$ITS)-1)+0.5)/length(roya$ITS)
roya$IT <- (roya$IT*(length(roya$IT)-1)+0.5)/length(roya$IT)

##########################################
## PREPROCESAMIENTO DE COVARIABLES
##########################################

#Ver el numero de predictores
nlayers(covariables)
covars <- na.omit(as.data.frame(covariables))
covars.sp <- na.omit(as.data.frame(covariables,xy=TRUE))

#Analisis de Componentes Principales
res.PCA <- PCA(na.omit(covars),scale.unit = TRUE,graph = FALSE)
summary(res.PCA)

#Para revisar los autovalores por dimension
eig.val <- get_eigenvalue(res.PCA)
eig.val

jpeg("F:/roya/ScientiaAgropecuaria2020/Figuras/eigenvalue.jpeg", width = 25, height = 17, units = 'cm', res = 600, pointsize = 50)
fviz_eig(res.PCA,choice = "eigenvalue",barcolor="black",linecolor = "red", 
         font.x=12,font.y=12,font.xtickslab=12,font.ytickslab=12) +
  labs(title ="", x = "Componentes", y = "Autovalor")
dev.off()

#Para revisar el aporte de las variables
var.val <- get_pca_var(res.PCA)
var.val$cor
write.csv(var.val$cor,file="F:/roya/ScientiaAgropecuaria2020/Rdata/PCA_corr.csv")

jpeg("F:/roya/ScientiaAgropecuaria2020/Figuras/biplot.jpeg", width = 25, height = 17, units = 'cm', res = 600, pointsize = 50)
fviz_pca_var(res.PCA, col.var = "black", repel = TRUE,
             font.x=12,font.y=12,font.xtickslab=12,font.ytickslab=12) +
  labs(title = "")
dev.off()

#Para revisar el aporte de las observaciones
ind.val <- get_pca_ind(res.PCA)
head(ind.val$coord) #estos son los Scores (nuevos valores de cada variable)
fviz_pca_ind(res.PCA)

#CREAR VARIABLES PCA 4
new.roya.pca <- ind.val$coord[,1:4]
new.roya.pca <- cbind(covars.sp[,1:2],new.roya.pca)
head(new.roya.pca)
coordinates(new.roya.pca) <- ~x+y
crs(new.roya.pca) <- crs(covariables)
gridded(new.roya.pca) <- TRUE
D1 <- raster(new.roya.pca,"Dim.1")
D2 <- raster(new.roya.pca,"Dim.2")
D3 <- raster(new.roya.pca,"Dim.3")
D4 <- raster(new.roya.pca,"Dim.4")
roya.PCA <- stack(D1,D2,D3,D4)


rm(eig.val,ind.val,numCores,var.val,covars,covars.sp,new.roya.pca,
   D1,D2,D3,D4,res.PCA)

rm(roya.PCA)
#############################
##EXTRAER COVARIABLES A ROYA
#############################
roya.sp <- roya
coordinates(roya.sp) <- ~X+Y
crs(roya.sp) <- crs(covariables)
roya2 <-extract(covariables, roya.sp)
roya3 <-extract(roya.PCA, roya.sp)
roya <- cbind(roya,roya2)
roya <- cbind(roya,roya3)
rm(roya2,roya3,roya.sp)

##############################
#### ANALISIS EXPLORATORIO ###
##############################

#Descripcion variables de roya
psych::describe(roya)
dt <- table(roya$ST) ; dt; prop.table(dt)
rm(dt)

#Autocorrelacion espacial global
roya.sp <- roya ; coordinates(roya.sp) <- ~X+Y ; crs(roya.sp) <- crs(covariables)
sw <- knearneigh(roya.sp, longlat = FALSE)
sw <- knn2nb(sw,sym = TRUE)
sw <- nb2listw(sw,zero.policy = TRUE)
moran.mc(x = roya$ITS,listw = sw,nsim=100000,adjust.n=TRUE)
joincount.mc(fx=roya$ST,listw = sw,nsim=100000)
rm(sw)

#Comparacion severidad
#ST <- read.delim("clipboard",header = TRUE)
#saveRDS(ST,file="F:/roya/ScientiaAgropecuaria2020/Rdata/SVfrecuency.rds")
ST <- readRDS("F:/roya/ScientiaAgropecuaria2020/Rdata/SVfrecuency.rds")
ggplot(data=ST,aes(x=Severidad,y=Fr,fill=Nivel)) +
  geom_bar("identity")

jpeg("F:/roya/ScientiaAgropecuaria2020/Figuras/ST_fr.jpeg", width = 25, height = 17, units = 'cm', res = 600, pointsize = 50)
ggplot(data=ST, aes(x=Severidad, y=Fr, fill=as.factor(Nivel))) +
  geom_bar(stat="identity",color="black") +
  scale_fill_brewer(palette="YlOrRd") +
  xlab("") + ylab("Frecuencia relativa (%)") +
  theme(legend.position="right",legend.title = element_blank(),legend.key.size = unit(0.5,"line"))
dev.off()

rm(ST)

################################################################################################
################################################################################################
################################################################################################

###############################
## MODELAMIENTO DE INCIDENCIA 
###############################
predictores <- c("NDVI","Elv","LST_B11","Asp","Slope")
md.full <- paste("IT ~ ", paste(predictores,collapse = " + ")) 
pairs(as.formula(md.full),data=roya)

md <- lm(IT ~ NDVI+Slope+LST_B11+poly(Elv,2)+poly(Asp,2),data=roya)
summary(md)
lmtest:: bptest(md)
car::vif(md)
car::durbinWatsonTest(md)


md0 <- betareg(IT ~ NDVI+Slope+LST_B11+Elv+Asp,data=roya,link="logit")
md1 <- betareg(ITM ~ NDVI+Slope+LST_B11+poly(Elv,2)+poly(Asp,2),data=roya,link="logit")
summary(md1) ; mean((roya$IT - predict(md1,roya))^2)

plot(IT ~ Elv, data = roya)
md_elv <- lm(IT~poly(Elv,2),data=roya)
lines(x=roya$Elv,y=predict(md1,roya))


#plot(md1, which = 1:4, type = "pearson")
#plot(md1, which = 5, type = "deviance", sub.caption = "")

set.seed(666)
md2 <- breg.cv(roya,fm.es = IT ~ NDVI+Slope+LST_B11+poly(Elv,2)+poly(Asp,2),respuesta = 16,k=10,t=1000,family="logit")
md2

library(ggplot2)
library(gridExtra)
p1 <- ggplot(data = roya, aes(Elv, residuals(md0,"pearson"))) +
  geom_point() + geom_smooth(color = "firebrick") + geom_hline(yintercept = 0) +
  theme_bw() + xlab("Elevación (m.s.n.m.)") + ylab("Residuales") + 
  ggtitle("Predictor Lineal")

p2 <- ggplot(data = roya, aes(Elv, residuals(md1,"pearson"))) +
  geom_point() + geom_smooth(color = "firebrick") + geom_hline(yintercept = 0) +
  theme_bw() + xlab("Elevación (m.s.n.m.)") + ylab("Residuales") +
  ggtitle("Predictor Cuadrático")

p3 <- ggplot(data = roya, aes(Asp, residuals(md0,"pearson"))) +
  geom_point() + geom_smooth(color = "firebrick") + geom_hline(yintercept = 0) +
  theme_bw() + xlab("Orientación (°)") + ylab("Residuales") 

p4 <- ggplot(data = roya, aes(Asp, residuals(md1,"pearson"))) +
  geom_point() + geom_smooth(color = "firebrick") + geom_hline(yintercept = 0) +
  theme_bw() + xlab("Orientación (°)") + ylab("Residuales")

jpeg("F:/roya/ScientiaAgropecuaria2020/Figuras/cuadratico.jpeg", width = 25, height = 17, units = 'cm', res = 600, pointsize = 50)
grid.arrange(p1,p2,p3,p4)
dev.off()

p6 <- ggplot(data = roya, aes(as.numeric(row.names(roya)), residuals(md1,type="deviance"),label=row.names(roya))) +
  geom_point() + geom_smooth(color = "firebrick") + geom_hline(yintercept = 0) +
  geom_hline(yintercept=c(-2,2), linetype="dashed", color = "red") +
  geom_text(aes(label=row.names(roya)),hjust=0, vjust=0) +
  theme_bw()
p6


#Diagnostico Inicial
df <- roya
md0 <- betareg(IT~1,data=df,link="log") ; AIC(md0)
md1 <- betareg(md.full,data=df,link = "log") ; AIC(md1)
plot(md1,1:6)

#Valores Leverage > 2*p/n son influyentes
HV <- hatvalues(md1)
Lev <- which(HV > 2*length(md1)/nrow(df) )

rm(md.full)

#Prediccion Incidencia
covars <- na.omit(as.data.frame(covariables,xy=TRUE))
covars[,"IT.pred"] <- predict(md1,covars,type="response")
covars[,"Var"] <- predict(md1,covars,type="variance")
covars.sp <- covars
coordinates(covars.sp) <- ~x+y ; crs(covars.sp) <- crs(covariables)
gridded(covars.sp) <- TRUE
IT <- raster(covars.sp,"IT.pred")
VAR <- raster(covars.sp,"Var")


aoi <- raster("F:/roya/AgriScientia/RAIFA_2019/DEM.tif")
aoi <- raster::resample(aoi,IT,method='bilinear')
IT <- mask(IT,aoi) ; VAR <- mask(IT,aoi)
library(rasterVis)
p1 <- rasterVis::levelplot(IT*100, par.settings = RdBuTheme,margin=FALSE,
                           colorkey=list(space="bottom")) +
  layer(sp.points(roya.sp, size=3,pch=20,col="black"))

################################################################################################
################################################################################################
################################################################################################

###############################
## MODELAMIENTO DE LA SEVERIDAD
###############################
library(VGAM)
library(ordinal)
mod0 <- vglm(ST ~ 1, family=cumulative(parallel=TRUE), data = roya)
mod1 <- vglm(ST ~ NDVI+Slope+LST_B11+poly(Elv,2)+poly(Asp,2), family=cumulative(parallel=TRUE), data = roya) 
summary(mod1) ;  1 - (logLik(mod1)/logLik(mod0))
AIC(mod1)

#Prediccion de variable ordinal
res.ord <- predict(mod1,covars,type="response")
#res.ord <- levels(roya$ST)[max.col(res.ord)] obtiene el vector directo de clases
res.ord <- max.col(res.ord)
res.ord <- cbind(covars[,1:2],res.ord)
coordinates(res.ord) <- ~x+y ; crs(res.ord) <- crs(covariables)
gridded(res.ord) <- TRUE
ST <- raster(res.ord,"res.ord")
ST <- mask(ST,aoi)
ST <- ratify(ST)
p2 <- rasterVis::levelplot(ST,att='ID',col.regions=brewer.pal(n = 4, name = "BrBG"),
                           colorkey=list(space="bottom"))+
  layer(sp.points(roya.sp, size=3,pch=20,col="black"))

wombo <- c(p1, p2, layout = c(2, 1), merge.legends = TRUE)

jpeg("F:/roya/ScientiaAgropecuaria2020/Figuras/Infeccion.jpeg", width = 25, height = 17, units = 'cm', res = 600, pointsize = 50)
print(wombo)
dev.off()

#Incertidumbre
ST_clases <- max.col(predict(mod1,roya,type="response"))
ST_clases <- ordered(ST_clases,levels = levels(roya$ST))
CM <- confusionMatrix(data=ST_clases,reference = roya$ST)
CM$overall

library(irr)
dt <- data.frame(a=ST_clases,b=roya$ST)
kp <- kappa2(dt, "equal")
kp$value


################################################################################################
################################################################################################
################################################################################################
### FUNCIONES DE APOYO ##


##############################################
#K-fold Cross validation for beta regression
##############################################
#data = data.frame con los datos
#fm.es = Estructura del modelo
#respuesta = indice de la variable respuesta en data
#k y t = fold y repeticiones
breg.cv <- function(data,fm.es,respuesta,k,t,...){
  require("caret")
  n <- nrow(data)
  particion <- createMultiFolds(y=data[,respuesta],k = k,times = t)
  modelos <- list()
  MSE_te <- rep(0,times=length(particion))
  pseR2 <- rep(0,times=length(particion))
  phi <- rep(0,times=length(particion))
  nk <- rep(0,times=length(particion))
  for (j in 1:length(particion)) {
    nk[j] <- n - length(particion[[j]])
  }
  for (i in 1:length(particion)) {
    modelos[[i]] <- betareg(formula = fm.es,data = data[particion[[i]],],link="logit")
    MSE_te[i] <- mean((data[-particion[[i]],respuesta] - predict(modelos[[i]],data[-particion[[i]],]))^2)
    pseR2[i] <- modelos[[i]]$pseudo.r.squared
    phi[i] <- modelos[[i]]$coefficients$precision
  }
  MSE_te <- matrix(MSE_te,nrow = t,ncol = k)
  pseR2 <- matrix(pseR2,nrow = t,ncol = k)
  phi <- matrix(phi,nrow = t,ncol = k)
  nk <- matrix(nk,nrow=t,ncol=k)
  MSE_CV <- rowSums(MSE_te*nk)/n
  pseR2_CV <- rowSums(pseR2*nk)/n
  phi_CV <- rowSums(phi*nk)/n
  return(list("MSE"=mean(MSE_CV),"MSE_sd"=sd(MSE_CV),"pR2"=mean(pseR2_CV),"pR2_sd"=sd(pseR2_CV),
              "phi"=mean(phi_CV),"phi_sd"=sd(phi_CV)))
}



breg.cv(roya,fm.es = md.full,respuesta = 16,k=5,t=10,family="log")


##############################################
#Cross validation ordinal
##############################################
Ord.cv <- function(data,fm.es,respuesta,k,t,...){
  require("caret")
  n <- nrow(data)
  particion <- createMultiFolds(y=data[,respuesta],k = k,times = t)
  modelos <- list()
  Acc <- rep(0,times=length(particion))
  Acc_pv <- rep(0,times=length(particion))
  Kap <- rep(0,times=length(particion))
  Kap_pv <- rep(0,times=length(particion))
  nk <- rep(0,times=length(particion))
  for (j in 1:length(particion)) {
    nk[j] <- n - length(particion[[j]])
  }
  for (i in 1:length(particion)) {
    modelos[[i]] <- vglm(fm.es, family=cumulative(parallel=TRUE), data = data[particion[[i]],])
    ST_clase <- max.col(predict(modelos[[i]],data[-particion[[i]],],type="response"))
    ST_clase <- ordered(ST_clases,levels = levels(data[,respuesta]))
    MC <- confusionMatrix(data=ST_clases,reference = data[-particion[[i]],respuesta])
    Acc[i] <- MC$overall[1]
    Acc_pv[i] <- MC$overall[6]
    Kap[i] <- kappa2(data.frame(a=ST_clases,b=data[-particion[[i]],respuesta]),"equal")$value
    kap_pv[i] <- kappa2(data.frame(a=ST_clases,b=data[-particion[[i]],respuesta]),"equal")$p.value

  }
  Acc <- matrix(Acc,nrow = t,ncol = k)
  Acc_pv <- matrix(Acc_pv,nrow = t,ncol = k)
  Kap <- matrix(Kap,nrow = t,ncol = k)
  Kap_pv <- matrix(Kap_pv,nrow = t,ncol = k)
  nk <- matrix(nk,nrow=t,ncol=k)
  Acc_CV <- rowSums(Acc*nk)/n
  Acc_pv_CV <- rowSums(Acc_pv*nk)/n
  Kap_CV <- rowSums(Kap*nk)/n
  Kap_pv_CV <- rowSums(Kap_pv*nk)/n
  return(list("Accuracy"=mean(Acc_CV),"Acc_sd"=sd(Acc_CV), "pvalor" = mean(Acc_pv_CV), "sd" = sd(Acc_pv_CV),
              "Kappa"=mean(Kap_CV),"KP_sd"=sd(Kap_CV), "pvalor" = mean(Kap_pv_CV), "sd" = sd(Kap_pv_CV)))

}

Ord.cv(roya,fm.es = ST ~ NDVI+Slope+LST_B11+poly(Elv,2)+poly(Asp,2),respuesta = 17,k=10,t=1)



##################################################
#### Beast Model subset for beta regression 
##################################################

#Crear una lista TRUE/FALSE con el numero de predictores
ve <- c(TRUE,FALSE)
lista <- rep(list(ve),length(predictores))
names(lista) <- predictores

#Crear una matriz TRUE/FALSE de las combinaciones 
require(utils)
regMat <- expand.grid(lista) ; rm(ve,lista)

#Crear la formula de todos los modelos
Mod.full <- apply(regMat, 1, function(x) as.formula(
  paste(c("IT ~ 1", predictores[x]),
        collapse=" + ")) )
length(Mod.full) #Cuantos hay

#Crear lista de modelos para cada enlace
modelos_logit <- list()
modelos_probit <- list()
modelos_cloglog <- list()
modelos_log <- list()
modelos_loglog <- list()

for (i in 1:length(Mod.full)) {
  modelos_logit[[i]] <- betareg(formula = Mod.full[[i]],data = roya,link = "logit")
  modelos_probit[[i]] <- betareg(formula = Mod.full[[i]],data = roya,link = "probit")
  modelos_cloglog[[i]] <- betareg(formula = Mod.full[[i]],data = roya,link = "cloglog")
  modelos_log[[i]] <- betareg(formula = Mod.full[[i]],data = roya,link = "log")
  modelos_loglog[[i]] <- betareg(formula = Mod.full[[i]],data = roya,link = "loglog")
}

#Base de datos de pseudoR2
for (i in 1:ncol(regMat)) {
  regMat[,i] <- as.numeric(regMat[,i])
}
regMat[,"N_predic"] <- base::rowSums(regMat)

psR <- matrix(0,nrow=length(Mod.full),ncol = 5)
for (i in 1:length(Mod.full)) {
  psR[i,1] <- modelos_logit[[i]]$pseudo.r.squared
  psR[i,2] <- modelos_probit[[i]]$pseudo.r.squared
  psR[i,3] <- modelos_cloglog[[i]]$pseudo.r.squared
  psR[i,4] <- modelos_log[[i]]$pseudo.r.squared
  psR[i,5] <- modelos_loglog[[i]]$pseudo.r.squared
}

regMat <- cbind(regMat,psR)
names(regMat) <- c("NDVI","Elv","LST_B11","Asp","Slope","N_pred","Logit","Probit","cloglog","Log","LogLog")
rm(psR)


############################################
### STEPWISE FORWARD BETAREG ###############
############################################
predictores <- names(roya[,18:ncol(roya)])

BRstep <- function(data,respuesta,predictores){
  M0 <- betareg(as.formula(paste(respuesta,"~ 1")),data=roya)
  M1 <- list()
  Akaike <- rep(0,length(predictores))
  for (i in 1:length(predictores)) {
    M1[[i]] <- betareg(as.formula(paste(respuesta," ~ ",predictores[i])),data=roya)
    Akaike[i] <- AIC(M1[[i]]) 
  }
  M1 <- M1[[which.min(Akaike)]]
  return(M1)
}

BRstep <- function(data,respuesta,predictores){
  #Crear todos los modelos
  ve <- c(TRUE,FALSE) ; lista <- rep(list(ve),length(predictores)) ;names(lista) <- predictores
  regMat <- expand.grid(lista)
  ModelosTotal <- apply(regMat, 1, function(x) as.formula(paste(c(paste(respuesta,"~ 1"), predictores[x]),collapse=" + ")))
  N_predi <- rowSums(regMat) #numero de predictores del modelo
  #Paso uno Modelo nulo
  M0 <- betareg(ModelosTotal[[which(N_predi==0)]],data=data)
  return(M0)
}


