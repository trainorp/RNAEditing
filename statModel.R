########### Prereqs ###########
library(tidyverse)

options(stringsAsFactors=FALSE)
setwd("/home/patrick/gdrive/RNAEditing/")

########### Process ###########
# Arguments (for when this is made into a function):
grp1<-"K562"
grp2<-"GM12787"

# Import tsv's:
alt<-read.table('alt_matrix.chr1.tsv',header=TRUE)
ref<-read.table('ref_matrix.chr1.tsv',header=TRUE)
tot<-read.table('total_matrix.chr1.tsv',header=TRUE)
ecpms<-read.table('ecpms.tsv',header=TRUE)

# Make into one dataframe:
alt<-alt %>% gather(key="sample",value="alt",-(chr:pos))
ref<-ref %>% gather(key="sample",value="ref",-(chr:pos))
df1<-alt %>% full_join(ref)
df1$tot<-df1$alt+df1$ref
df1$prop<-df1$alt/df1$tot
df1$psuAlt<-(df1$alt+.1*df1$tot)/(1.1*df1$tot)
df1$psuRef<-(df1$ref+.1*df1$tot)/(1.1*df1$tot)
df1$lo<-log(df1$psuAlt/df1$psuRef)

# Add group and replicate:
df1$group<-str_split(df1$sample,"_",simplify=TRUE)[,1]
df1$rep<-str_split(df1$sample,"_",simplify=TRUE)[,2]

# Join ecpms:
df1<-df1 %>% left_join(ecpms,by=c("sample"="Sample_Name"))

# Change position to character
df1$pos<-as.character(df1$pos)

# As factor:
df1$group<-factor(df1$group,levels=c("K562","GM12787"))
df1$grp<-ifelse(df1$group==grp1,"grp1","grp2")

########### Exploratory analysis ###########
df1$nonZero<-as.integer(df1$tot>0)
xtabs(~nonZero+group,data=df1)
nonZeroDf<-as.data.frame(prop.table(xtabs(~nonZero+sample,data=df1)))
nonZeroDf<-nonZeroDf %>% filter(nonZero==TRUE) %>% select(-nonZero)

# Min information filter:
df1<-df1 %>% left_join(df1 %>% group_by(pos,group) %>% 
                    summarize(nonZeroCount=sum(nonZero)) %>% 
                    group_by(pos) %>% 
                    summarize(noInfo=sum(as.integer(nonZeroCount<3))>0))

ggplot(df1 %>% filter(!noInfo),aes(group=group,x=lo,y=..density..,fill=group))+
  geom_histogram(binwidth=.05,position="identity",alpha=.5)+
  theme_bw()

ggplot(df1 %>% filter(!noInfo),aes(group=group,x=lo,y=..density..,fill=group))+
  geom_histogram(binwidth=.1,position="identity",alpha=.5)+
  theme_bw()

########### Model ###########
lm1<-lm(lo~EPCM,data=df1[!df1$noInfo,])
coefLm1<-coef(lm1)
df1$off<-coefLm1['EPCM']*df1$EPCM

df2<-df1 %>% select(chr,pos,noInfo) %>% distinct()
df2$propGrp1<-df2$propGrp2<-df2$pVal<-NA
for(i in 1:nrow(df2)){
  if(!df2$noInfo[i]){
    df3<-df1 %>% filter(pos==df2$pos[i])
    df3$w<-df3$tot/sum(df3$tot)
    df4<-df3 %>% group_by(grp) %>% summarize(mean=mean(prop,na.rm=TRUE))
    lm3<-lm(lo~off+group,weights=w,data=df3)
    df2$pVal[i]<-anova(lm3)['group','Pr(>F)']
    df2$propGrp1[i]<-df4$mean[df4$grp=="grp1"]
    df2$propGrp2[i]<-df4$mean[df4$grp=="grp2"]
  }
}
