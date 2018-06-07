############# Prereqs #############
library(tidyverse)
library(lme4)
library(lmerTest)
library(zoib)

options(stringsAsFactors=FALSE)
oldPar<-par()
setwd('~/gdrive/RNAEditing/')

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

# Cell line:
df1$cellLine<-sapply(str_split(df1$pheno,"_"),function(x) x[1])
df1$type<-sapply(str_split(df1$pheno,"_"),function(x) x[2])

############# Model #############
df1$Area<-log(df1$Area)
df1$areaSq<-df1$Area**2
df1$coverageSq<-df1$Coverage**2
df1$areaSqR<-sqrt(df1$Area)
df1$coverageSqR<-sqrt(df1$Coverage)
df1$areaCub<-df1$Area**3

df1<-df1[df1$opt=="c5e3",]
lm3me<-lmer(Counts~Area+areaSq+areaCub+Coverage+Coverage:Area+Area*type+(1|id),data=df1)

coef1<-as.matrix(coef(lm3me)$id[1,])
df1$pred<-mean(coef(lm3me)$id[,1])+coef1[,'Area']*df1$Area+coef1[,'areaSq']*df1$areaSq+
  coef1[,'areaCub']*df1$areaCub+coef1[,'Coverage']*df1$Coverage+
  coef1[,'Area:Coverage']*df1$Area*df1$Coverage+
  ifelse(df1$type=="total",1,0)*coef1[,'typetotal']+
  ifelse(df1$type=="total",1,0)*df1$Area*coef1[,'Area:typetotal']

############# Statistical test #############
df1$norm<-df1$pred-df1$Counts
lm4me<-update(lm3me,.~.+cellLine)
summary(lm4me)
anova(lm4me)

idNorm<-df1 %>% group_by(id) %>% summarize(mean(norm))

############# Sites #############
sites<-read.table('all_huvec_site_stats.srt.ann.tsv')
names(sites)<-c("Sample","Group","Island_ID","Location","Transition_Type",
                "Edited_Bases", "Total_Mapped_Bases")

sites$prop<-sites$Edited_Bases/sites$Total_Mapped_Bases
sites$prop[is.nan(sites$prop)]<-NA

############# Sites #############
sites<-sites %>% group_by(Location) %>% 
  mutate(locN=sum(ifelse(Total_Mapped_Bases>0,1L,0L)))
sites<-sites %>% group_by(Island_ID) %>% 
  mutate(islandN=sum(ifelse(Total_Mapped_Bases>0,1L,0L)))

hist(sites$prop)
