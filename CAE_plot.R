library(tidyverse)
library(scatterplot3d)
library(rgl)
library(lme4)

options(stringsAsFactors=FALSE)
oldPar<-par()
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

# Bind and replace names:
df1<-rbind(df1a,df1b,df1c)
names(df1)[match(names(df1),c("V1","V2","V3"),nomatch=FALSE)]<-c("Coverage","Counts","Area")

# Add color:
pal1<-RColorBrewer::brewer.pal(3,"Set1") # Red, blue, green
df1$color<-pal1[1]
df1$color[df1$cellLine=="K562"]<-pal1[2]
df1$color[df1$type=="totalRNA"]<-pal1[3]

############# CAE plots #############
png(file="plots/CovAreaCount_Color.png",height=5,width=6,unit="in",res=300)
par(las=2)
idk<-scatterplot3d(x=df1$Coverage,y=df1$Area,z=df1$Counts,color=df1$color,
              xlab="Coverage",ylab="Area",zlab="Counts",angle=40,
              cex.axis=.7)
dev.off()

scatterplot3d(x=df1$Coverage,y=df1$Area,z=df1$Counts/df1$Area,color=df1$color,
                   xlab="Coverage",ylab="Area",zlab="Rate",angle=50,
                   cex.axis=.7)

plot3d(x=df1$Coverage,y=df1$Area,z=df1$Counts, col=df1$color,type = "p")
plot3d(x=df1$Coverage,y=log(df1$Area),z=df1$Counts, col=df1$color,type = "p")

############# Model #############
df1$areaSq<-df1$Area**2
df1$coverageSq<-df1$Coverage**2
df1$areaSqR<-sqrt(df1$Area)
df1$coverageSqR<-sqrt(df1$Coverage)
df1$areaCub<-df1$Area**3

lmBig<-lm(Counts~Area+Coverage+areaSq+areaSqR+coverageSq+coverageSqR,data=df1)
summary(lmBig)

lmBig2<-update(lmBig,.~.-areaSqR-coverageSq-coverageSqR)
summary(lmBig2)

lm1<-lm(Counts~Area*Coverage+areaSq,data=df1)
summary(lm1)

lm2<-lm(Counts~Area+areaSq+areaCub,data=df1)
summary(lm2)

lm3<-update(lm2,.~.-areaCub)
anova(lm3,lm2)

lm3me<-lmer(Counts~Area+areaSq+areaCub+(1|cellLine),data=df1)

df1$pred<-predict(lm3)

plot3d(x=df1$Coverage,y=df1$Area,z=df1$Counts, col=df1$color,type = "p")
points3d(x=df1$Coverage,y=df1$Area,z=df1$pred)

df1$normCounts<-df1$Counts-df1$pred
plot3d(x=df1$Coverage,y=df1$Area,z=df1$normCounts,col=df1$color,type = "p")

plot(x=df1$Area,y=df1$normCounts/df1$Area,col=df1$color)
plot(x=df1$Area,y=df1$Counts/df1$Area,col=df1$color)
