#CARGA DE LIBRERIAS
library(corrplot)
library(sqldf)
library(bnlearn)
library(visNetwork)

#Se carga libreria para carga de datos desde github
library(RCurl)

#Carga de la base de datos
url<-"https://raw.githubusercontent.com/Ricardo27cruz27/Redes_Bayesianas/master/diabetes.csv"
url_csv<-getURL(url)
datos<-read.csv(text=url_csv,header = F)
#nombres de las variables
names(datos)<-c("Embarazo",
                "Glucosa",
                "Presion",
                "Tricep",
                "Insulina",
                "IMC",
                "Pedigree",
                "Edad",
                "Diabetes")
attach(datos)


#Análisis exploratorio
summary(datos[,1])
sort(unique(datos[,1]))
#hist(datos[,1])

summary(datos[,2])
sort(unique(datos[,2]))
#hist(datos[,2])

#NA'S
datos.limp<-datos
datos.limp[(which(Glucosa==0)),2]=NA
datos.limp[(which(Presion==0)),3]=NA
datos.limp[(which(Tricep==0)),4]=NA
datos.limp[(which(Insulina==0)),5]=NA
datos.limp[(which(IMC==0)),6]=NA


#Categorización.
#Embarazos
datos.limp[(which(datos.limp[,1]<=1)),1]=0
datos.limp[(which(datos.limp[,1]>1 & datos.limp[,1]<=4)),1]=1
datos.limp[(which(datos.limp[,1]>4)),1]=2

#Glucosa
datos.limp[(which(datos.limp[,2]<140)),2]=0
#datos.limp[(which(datos.limp[,2]>70 & datos.limp[,2]<140)),2]=1
datos.limp[(which(datos.limp[,2]>=140)),2]=1

#Presion
datos.limp[(which(datos.limp[,3]<=60)),3]=0
datos.limp[(which(datos.limp[,3]>60 & datos.limp[,3]<90)),3]=1
datos.limp[(which(datos.limp[,3]>=90)),3]=2

#Tricep
datos.limp[(which(datos.limp[,4]<=25)),4]=0
datos.limp[(which(datos.limp[,4]>25)),4]=1

#Insulina
datos.limp[(which(datos.limp[,5]<=120)),5]=0
datos.limp[(which(datos.limp[,5]>120)),5]=1

#IMC
datos.limp[(which(datos.limp[,6]<=25)),6]=0
datos.limp[(which(datos.limp[,6]>25 & datos.limp[,6]<=30)),6]=1
datos.limp[(which(datos.limp[,6]>30)),6]=2

#Pedigree
datos.limp[(which(datos.limp[,7]>quantile(Pedigree,.75))),7]=3
datos.limp[(which(datos.limp[,7]>quantile(Pedigree,.5) &
                    datos.limp[,7]<=quantile(Pedigree,.75))),7]=2
datos.limp[(which(datos.limp[,7]>quantile(Pedigree,.25) &
                    datos.limp[,7]<=quantile(Pedigree,.5))),7]=1
datos.limp[(which(datos.limp[,7]<=quantile(Pedigree,.25))),7]=0

#Edad
datos.limp[(which(datos.limp[,8]<=41)),8]=0
#datos.limp[(which(datos.limp[,8]>41 & datos.limp[,8]<=61)),8]=1
datos.limp[(which(datos.limp[,8]>41)),8]=1

#Correlación
correlacion<-cor(as.matrix(datos.limp),use = "pairwise.complete.obs")
corrplot.mixed(correlacion, upper="square",lower="number")
#corrplot(correlacion)

#D como el conjunto de datos
D<-datos.limp
#como factores
D[names(datos)]<-lapply(D[names(datos)],as.factor)
sapply(D, class)


#Conjunto de entrenamiento y prueba
set.seed(27)
muestra<-sample(1:nrow(D),size = 614)
D.train<-D[muestra,]
D.test<-D[-muestra,]

nrow(D.train); nrow(D.test)

#Aprendizaje estructural
#Generacion de un modelo a través de la 
#muestra de entrenamiento

#Funcion para plotear redes:
plot.network <- function(structure, ht = "400px",title,subtitle){
  nodes.uniq <- unique(c(structure$arcs[,1], structure$arcs[,2]))
  nodes <- data.frame(id = nodes.uniq,
                      label = nodes.uniq,
                      color = "darkturquoise",
                      shadow = TRUE,
                      shape="circle")
  
  edges <- data.frame(from = structure$arcs[,1],
                      to = structure$arcs[,2],
                      arrows = "to",
                      smooth = TRUE,
                      shadow = TRUE,
                      color = "black")
  
  return(visNetwork(nodes, edges, height = ht, width = "100%",main=title,submain=subtitle
  ) %>% visLayout(randomSeed = 27) 
  )
}

#RB sin ninguna restricción
RB.sinrest<-structural.em(D.train,return.all = TRUE,maximize.args = list(restart=100))
plot.network(RB.sinrest$dag,title = "",subtitle = "")

#funcion para fuerza de los arcos
fuerza.bootsrap<-function(RB,args=list()){
  arcs = boot.strength(RB$imputed, algorithm = "hc",algorithm.args = args)
  print(arcs[(arcs$strength > 0.5) & (arcs$direction >= 0.5), ])
  #prueba<-(averaged.network(arcs,threshold = .5))
  return(arcs)
}
boot.RB.sinrestricciones<-fuerza.bootsrap(RB.sinrest)
boot.RB.sinrestricciones[which(boot.RB.sinrestricciones$from=="Diabetes"),
                         ]

#RED CON RESTRICCIONES
## ARCOS PROHIBIDOS:
blacklist<-data.frame(from=names(D)[-8],to=rep("Edad",8))
whitelist4=data.frame(from=c("IMC","Glucosa","Pedigree","IMC","Embarazo"),
                      to=c("Diabetes","Diabetes","Diabetes","Tricep","Tricep"))
RB.b1w4<-structural.em(D.train,
                       maximize.args = list(blacklist=blacklist,restart=100,whitelist=whitelist4),
                       return.all = TRUE)
plot.network(RB.b1w4$dag,title = "",subtitle = "")
#bootsrap:
boot.RB.b1w4<-fuerza.bootsrap(RB.b1w4,list(whitelist=whitelist4,blacklist=blacklist))


#Tablas de probabilidad condicionada:
RB.b1w4$fitted$Pedigree
RB.b1w4$fitted$Edad
RB.b1w4$fitted$Presion
RB.b1w4$fitted$Embarazo
RB.b1w4$fitted$Glucosa
RB.b1w4$fitted$IMC
RB.b1w4$fitted$Tricep
RB.b1w4$fitted$Insulina
RB.b1w4$fitted$Diabetes



#INFERENCIA

#carga de librerias
library(BiocManager)
library(caret)
BiocManager::install("RBGL")
library(gRbase)
library(gRain)

#Redes como grain
gr.RB1<-as.grain(RB.sinrest$fitted)
gr.RB10<-as.grain(RB.b1w4$fitted)

#redes a probar
redes<-list(gr.RB1,gr.RB10)
#Funcion de propagacion
prop<-function(red){
  propagacion=red
  prediccion.train<-predict(propagacion,response = "Diabetes",newdata = D.train,predictors = names(D)[-9],type="distribution")
  pred.train<-as.integer(prediccion.train$pred$Diabetes)
  ac.train<-sum(Diab.train==pred.train)/nrow(D.train)
  
  prediccion.test<-predict(propagacion,response = "Diabetes",newdata = D.test,predictors = names(D)[-9],type="distribution")
  pred.test<-as.integer(prediccion.test$pred$Diabetes)
  ac.test<-sum(Diab.test==pred.test)/nrow(D.test)
  return(list(prediccion.train$pred$Diabetes,prediccion.test$pred$Diabetes))
}

ACC.PRU<-lapply(redes, prop)

#Probabilidades 
ACC.PRU[[2]][1]

#red
rn<-2
#punto de corte
co<-.5
r.train<-as.integer(ACC.PRU[[rn]][[1]][,2]>co)
#matriz de confusion
cm.train<-table(r.train,D.train$Diabetes)
cm.train<-confusionMatrix(cm.train,positive = "1")
cm.train$overall[1]
cm.train$byClass[1]

#gráfico
his.train<-ACC.PRU[[2]][[1]][,2]
ceros.train<-his.train[which(Diab.train==0)]
unos.train<-his.train[which(Diab.train==1)]
dc.train<-density(ceros.train,from = 0,to = 1)
du.train<-density(unos.train,from = 0,to = 1)
x11()
plot(dc.train,xlim=c(0,1),ylim=c(0,2.5),main="",xlab="p()",ylab="",col="darkblue")
par(new=T)
plot(du.train,xlim=c(0,1),ylim=c(0,2.5),main="",xlab="p()",ylab="",col="red")
abline(v=.5,col="green")


#Cambio del punto de corte
his.train<-ACC.PRU[[2]][[1]][,2]
ceros.train<-his.train[which(Diab.train==0)]
unos.train<-his.train[which(Diab.train==1)]
dc.train<-density(ceros.train,from = 0,to = 1)
du.train<-density(unos.train,from = 0,to = 1)
m<-max(du.train$y[which(du.train$x<.35)])
cutoff<-du.train$x[which(du.train$y==m)]
x11()
plot(dc.train,xlim=c(0,1),ylim=c(0,2.5),main="",xlab="p()",ylab="",col="darkblue")
par(new=T)
plot(du.train,xlim=c(0,1),ylim=c(0,2.5),main="",xlab="p()",ylab="",col="red")
abline(v=cutoff,col="green")

his.test<-ACC.PRU[[2]][[2]][,2]
ceros.test<-his.test[which(Diab.test==0)]
unos.test<-his.test[which(Diab.test==1)]
dc.test<-density(ceros.test)
du.test<-density(unos.test)
x11()
plot(dc.test,xlim=c(0,1),ylim=c(0,2.5),main="",xlab="p()",ylab="",col="darkblue")
par(new=T)
plot(du.test,xlim=c(0,1),ylim=c(0,2.5),main="",xlab="p()",ylab="",col="red")
abline(v=cutoff,col="green")

#precision y sensibilidad con diferente punto de corte
rn<-2
co<-cutoff
r.train<-as.integer(ACC.PRU[[rn]][[1]][,2]>co)
r.test<-as.integer(ACC.PRU[[rn]][[2]][,2]>co)

cm.train<-table(r.train,D.train$Diabetes)
cm.train<-confusionMatrix(cm.train,positive = "1")

cm.test<-table(r.test,D.test$Diabetes)
cm.test<-confusionMatrix(cm.test,positive = "1")

cm.train
cm.test

#validación cruzada
estructura<-RB.b1w4$dag
#tamaños de los folds
tamaños<-c(154,153,154,153,154)
#entrenamiento de cada iteración
sets.test<-list(B1=D[0:cumsum(tamaños)[1],],
                B2=D[(cumsum(tamaños)[1]+1):cumsum(tamaños)[2],],
                B3=D[(cumsum(tamaños)[2]+1):cumsum(tamaños)[3],],
                B4=D[(cumsum(tamaños)[3]+1):cumsum(tamaños)[4],],
                B5=D[(cumsum(tamaños)[4]+1):cumsum(tamaños)[5],])
#prueba de cada iteración
sets.train<-list(B1=D[-(0:cumsum(tamaños)[1]),],
                 B2=D[-((cumsum(tamaños)[1]+1):cumsum(tamaños)[2]),],
                 B3=D[-((cumsum(tamaños)[2]+1):cumsum(tamaños)[3]),],
                 B4=D[-((cumsum(tamaños)[3]+1):cumsum(tamaños)[4]),],
                 B5=D[-((cumsum(tamaños)[4]+1):cumsum(tamaños)[5]),])
#Función de cada iteración en el cv
prueba.cv<-function(red,entrenamiento,pruebas){
  propagacion=red
  prediccion.train<-predict(propagacion,response = "Diabetes",newdata = entrenamiento,predictors = names(D)[-9],type="distribution")
  #pred.train<-as.integer(prediccion.train$pred$Diabetes)
  #ac.train<-sum(Diab.train==pred.train)/nrow(entrenamiento)
  
  prediccion.test<-predict(propagacion,response = "Diabetes",newdata = pruebas,predictors = names(D)[-9],type="distribution")
  #pred.test<-as.integer(prediccion.test$pred$Diabetes)
  #ac.test<-sum(Diab.test==pred.test)/nrow(pruebas)
  return(list(prediccion.train$pred$Diabetes,prediccion.test$pred$Diabetes))
}
#matriz de resumen del cv
Resumen<-matrix(0,5,5)
colnames(Resumen)<-c("Punto de corte",
                     "Acc train",
                     "Sensit train",
                     "Acc test",
                     "Sensit test")
#CV
for(i in 1:5){
  red.ajustada<-bn.fit(x = estructura,
                       data =sets.train[[i]])
  grain.ajustada<-as.grain(red.ajustada)
  probability<-prueba.cv(grain.ajustada,
                         sets.train[[i]],
                         sets.test[[i]])
  his.train<-probability[[1]][,2]
  ceros.train<-his.train[which(sets.train[[i]][,9]==0)]
  unos.train<-his.train[which(sets.train[[i]][,9]==1)]
  dc.train<-density(ceros.train,from = 0,to=1)
  du.train<-density(unos.train,from = 0,to=1)
  m<-max(du.train$y[which(du.train$x<.35)])
  cutoff<-du.train$x[which(du.train$y==m)]
  #m<-which.min(abs(du.train$y-dc.train$y))
  #cutoff<-du.train$x[m[1]]
  x11()
  plot(dc.train,xlim=c(0,1),ylim=c(0,2.5),main=paste("train",i),xlab="p()",ylab="",col="darkblue")
  par(new=T)
  plot(du.train,xlim=c(0,1),ylim=c(0,2.5),main=paste("train",i),xlab="p()",ylab="",col="red")
  abline(v=cutoff,col="green")
  
  his.test<-probability[[2]][,2]
  ceros.test<-his.test[which(sets.test[[i]][,9]==0)]
  unos.test<-his.test[which(sets.test[[i]][,9]==1)]
  dc.test<-density(ceros.test)
  du.test<-density(unos.test)
  x11()
  plot(dc.test,xlim=c(0,1),ylim=c(0,2.5),main=paste("test",i),xlab="p()",ylab="",col="darkblue")
  par(new=T)
  plot(du.test,xlim=c(0,1),ylim=c(0,2.5),main=paste("test",i),xlab="p()",ylab="",col="red")
  abline(v=cutoff,col="green")
  #cortes<-c(cortes,cutoff)
  
  r.train<-as.integer(probability[[1]][,2]>cutoff)
  r.test<-as.integer(probability[[2]][,2]>cutoff)
  
  cm.train<-table(r.train,sets.train[[i]][,9])
  cm.train<-confusionMatrix(cm.train,positive = "1")
  
  cm.test<-table(r.test,sets.test[[i]][,9])
  cm.test<-confusionMatrix(cm.test,positive = "1")
  
  Resumen[i,1]<-cutoff
  Resumen[i,2]<-cm.train$overall[1]
  Resumen[i,3]<-cm.train$byClass[1]
  Resumen[i,4]<-cm.test$overall[1]
  Resumen[i,5]<-cm.test$byClass[1]
}

#Resumen
row.names(Resumen)<-c("fold1","fold2","fold3","fold4","fold5")
Resumen
colMeans(Resumen)




#consultas:
diab1<-D[which(D$Diabetes==1),c("Diabetes","Presion")]
pred<-predict(gr.RB10,
              response = "Presion",
              newdata = diab1,
              predictors = "Diabetes",
              type="distribution")
pred$pred


querygrain(gr.RB10, 
           nodes = c("Diabetes", "IMC"),
           type = "joint")
gr.RB10.evi<-setFinding(gr.RB10,
                        nodes=c("Embarazo","Edad"),
                        states=c("1","1"))
querygrain(gr.RB10.evi, 
           nodes = c("Diabetes", "IMC"),
           type = "joint")



