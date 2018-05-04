library(tidyverse)
library(plot3D)

options(stringsAsFactors=FALSE)

setwd('~/gdrive/RNAEditing/')
df1<-read.table('editing_levels_counts_and_area.tsv',header=FALSE)
names(df1)<-c("EditingRate","Counts","Area")

png(file="CountsAreaEditing3d.png",height=6,width=6,unit="in",res=300)
scatter3D(x=df1$Counts,y=df1$Area,z=df1$EditingRate,
         xlab="Counts",ylab="Area",zlab="Editing Rate",
         phi=25,theta=45)
dev.off()

lm1<-lm(EditingRate~Counts+Area,data=df1)
summary(lm1)
anova(lm1)

gridLines<-30
xPred<-seq(min(df1$Counts),max(df1$Counts),length.out=gridLines)
yPred<-seq(min(df1$Area),max(df1$Area),length.out=gridLines)
xyPred<-expand.grid(Counts=xPred,Area=yPred)
zPred<-matrix(predict(lm1,newdata=xyPred),nrow=gridLines,ncol=gridLines)

scatter3D(x=df1$Counts,y=df1$Area,z=df1$EditingRate,
          xlab="Counts",ylab="Area",zlab="Editing Rate",
          phi=25,theta=45,surf=list(x=xPred,y=yPred,z=zPred,col="grey"))
surf3D(x=xPred,y=yPred,z=zPred,add=TRUE)
