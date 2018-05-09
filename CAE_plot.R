library(tidyverse)
library(scatterplot3d)
library(rgl)
library(lme4)

options(stringsAsFactors=FALSE)
oldPar<-par()
setwd('~/gdrive/RNAEditing/')

############# Import data #############
files<-list.files(path="split_by_file",recursive=TRUE)

temp2<-list()
for(i in 1:length(files)){
  temp<-read.table(file=paste0("split_by_file/",files[i]),header=FALSE)
  splitFile<-unlist(str_split(files[i],"/"))
  temp$pheno<-splitFile[1]
  names(temp)[match(names(temp),c("V1","V2","V3"),nomatch=FALSE)]<-
    c("Coverage","Counts","Area")
  temp$id<-unlist(str_split(splitFile[2],"\\."))[1]
  temp$opt<-unlist(str_split(unlist(str_split(splitFile[2],"_"))[3],"\\."))[1]
  temp2[[i]]<-temp
}
df1<-do.call("rbind",temp2)

table(df1$opt)

# Cell line:
df1$cellLine<-sapply(str_split(df1$pheno,"_"),function(x) x[1])
df1$type<-sapply(str_split(df1$pheno,"_"),function(x) x[2])

# Add color:
pal1<-RColorBrewer::brewer.pal(4,"Set1") # Red, blue, green
df1$color<-pal1[1]
df1$color[df1$pheno=="GM12787_total"]<-pal1[2]
df1$color[df1$pheno=="K562_polya"]<-pal1[3]
df1$color[df1$pheno=="K562_total"]<-pal1[4]

# First analysis:
df0<-df1
df1<-df0[df0$opt=="c5e3",]

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

par3d(windowRect = c(0,50,612,612))
parent <- currentSubscene3d()
plot3d(x=df1$Coverage,y=df1$Area,z=df1$Counts, col=df1$color,type = "p")
legend3d("topright",legend=names(table(df1$pheno)),pch=16,
         col=pal1,cex=1,inset=c(0.02))

############# Model #############
df1$Area<-log(df1$Area)
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

lm3me<-lmer(Counts~Area+areaSq+areaCub+(1|id),data=df1)
coef1<-as.matrix(coef(lm3me)$id[1,])
df1$pred<-mean(coef(lm3me)$id[,1])+coef1[,'Area']*df1$Area+coef1[,'areaSq']*df1$areaSq+
  coef1[,'areaCub']*df1$areaCub

plot3d(x=df1$Coverage,y=df1$Area,z=df1$Counts, col=df1$color,type = "p")
points3d(x=df1$Coverage,y=df1$Area,z=df1$pred)

############# Statistical test #############
lm4me<-update(lm3me,.~.+cellLine)
summary(lm4me)
