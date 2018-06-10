############# Prereqs #############
library(tidyverse)
library(lme4)
library(lmerTest)
library(betareg)

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

############# Site imputation #############
# Counts: 
sites<-sites %>% group_by(Location,Group) %>% 
  mutate(GroupLocN=sum(ifelse(Total_Mapped_Bases>0,1L,0L)))
sites<-sites %>% group_by(Location) %>% 
  mutate(minGroupLocN=min(GroupLocN),locN=sum(ifelse(Total_Mapped_Bases>0L,1L,0L)))

sites<-sites %>% group_by(Island_ID,Group) %>% 
  mutate(groupIslandN=sum(ifelse(Total_Mapped_Bases>0L,1L,0L)))
sites<-sites %>% group_by(Location) %>% 
  mutate(minGroupIslandN=min(groupIslandN),
         groupIslandN=sum(ifelse(Total_Mapped_Bases>0L,1L,0L)))

# Minimum value imputation function 
minImpFun<-function(x){
  if(all(is.na(x))) return(NA)
  else{
    x2<-x[complete.cases(x)]
    x2<-x2[x2>0]
    if(length(x2)>0) min(x2)
    else return(0)
  }
}
sites<-sites %>% group_by(Location) %>% mutate(minSiteProp=minImpFun(prop))
sites<-sites %>% group_by(Island_ID) %>% mutate(minIslandProp=minImpFun(prop))

# Test logical:
sites$siteTest<-sites$minGroupLocN>=2L & sites$minSiteProp>0
sites$islandTest<-sites$minGroupIslandN>=2L & sites$minIslandProp>0

# Add pseudo-props & Logits:
sites$prop2<-sites$prop
sites$prop2[is.na(sites$prop) & !is.na(sites$minSiteProp)]<-
  sites$minSiteProp[is.na(sites$prop) & !is.na(sites$minSiteProp)]
sites$prop2<-(sites$prop2*99+.5)/100

########### Beta reg model ###########
# Order by location
sites$chr<-str_split(sites$Location,":",simplify=TRUE)[,1]
sites$loc<-as.integer(str_split(sites$Location,":",simplify=TRUE)[,2])
sites<-sites %>% arrange(chr,loc)
sites$location2<-factor(sites$Location,levels=unique(sites$Location))

sites2<-sites %>% filter(siteTest)

locs<-unique(sites2$Location)
sites3<-sites2 %>% filter(Location %in% locs[1:500])
br1<-betareg(prop2~Group+Location|Total_Mapped_Bases,data=sites3)
