roya <- read.delim("clipboard",header=TRUE)

#Regresion Poisson simple
md1 <- glm( formula = X.HE1 ~ Z, family = poisson(link = log),data=roya)
summary(md1)
source('http://www.poleto.com/funcoes/envel.pois.txt')
source('http://www.poleto.com/funcoes/diag.pois.txt')
envel.pois(md1)
diag.pois(md1)

phi<-md1$deviance/md1$df.residual
phi

#Regresion Poisson con offset
md2 <- glm(formula = X.HE1 ~ Z +offset(log(X.HT1)), family = poisson(link = log),data=roya)
summary(md2)
envel.pois(md2)
diag.pois(md2,iden=c(2,2,2,0,3,3,3,0))

phi<-md2$deviance/md2$df.residual
phi

#Regresion binomial negativa
library(MASS)
md3=glm.nb(X.HE1 ~ Z,link=log,data=roya)
summary(md3)
source('http://www.poleto.com/funcoes/diag.nb.txt')
diag.nb(md3)
source('http://www.poleto.com/funcoes/envel.nb.txt')
envel.nb(md3)
