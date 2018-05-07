library(tidyverse)
library(plot3D)

options(stringsAsFactors=FALSE)

setwd('~/gdrive/RNAEditing/')

############# Import data #############
df1a<-read.table('GM12787_polya.tsv',header=FALSE)
df1a$cellLine<-"GM12787"
df1a$type<-"polyA"
df1b<-read.table('K562_polya.tsv',header=FALSE)
df1b$cellLine<-"K562"
df1b$type<-"polyA"
df1c<-read.table('GM12787_totalrna.tsv',header=FALSE)
df1c$cellLine<-"GM12787"
df1c$type<-"totalRNA"

df1<-rbind(df1a,df1b,df1c)
names(df1)[match(names(df1),c("V1","V2","V3"),nomatch=FALSE)]<-c("Coverage","Counts","Area")

############# CAE plots #############
CAEGrid<-expand.grid(phi=seq(0,2*pi,pi/4)*180/pi,theta=seq(0,2*pi,pi/4)*180/pi)
for(i in 1:nrow(CAEGrid)){
  file<-paste0("plots/CountsAreaEditing3d_",CAEGrid$phi[i],"_",CAEGrid$theta[i],".png")
  png(file=file,height=5,width=6,unit="in",res=300)
  scatter3D(z=df1$Counts,y=df1$Area,x=df1$Coverage,
            zlab="Counts",ylab="Area",xlab="Coverage",
            phi=CAEGrid$phi[i],theta=CAEGrid$theta[i])
  dev.off()
}

lm1<-lm(Coverage~Counts+Area,data=df1)
summary(lm1)
anova(lm1)

gridLines<-30
xPred<-seq(min(df1$Counts),max(df1$Counts),length.out=gridLines)
yPred<-seq(min(df1$Area),max(df1$Area),length.out=gridLines)
xyPred<-expand.grid(Counts=xPred,Area=yPred)
zPred<-matrix(predict(lm1,newdata=xyPred),nrow=gridLines,ncol=gridLines)

scatter3D(x=df1$Counts,y=df1$Area,z=df1$Coverage,
          xlab="Counts",ylab="Area",zlab="Editing Rate",
          phi=25,theta=45,surf=list(x=xPred,y=yPred,z=zPred,col="grey"))
surf3D(x=xPred,y=yPred,z=zPred,add=TRUE)
