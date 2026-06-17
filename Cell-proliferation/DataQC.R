setwd("/Users/maria/Documents/work/Weis_lab/Acropora larvae infections/analysis_Nov2024")

#larv=read.csv("data/LarvEduSizeSym_wMeta_final.csv")
larv=read.csv("data/LarvEduSizeSym_wMeta_6Dec2024.csv")
larv[,c("symClade","symTreat","timepoint")]=lapply(larv[,c("symClade","symTreat","timepoint")],as.factor)
larv$symClade=factor(larv$symClade,levels=c('apo','A','B','C','D'))
larv$symTreat=factor(larv$symTreat,levels=c('apo','25','32'))
larv$sym.size=larv$sym.count/larv$size
larv$Edu.size=larv$Edu/larv$size

library(dplyr)
sample_size = larv %>% group_by(symTreat,symClade,timepoint) %>% summarise(n=n())
larv[!is.na(larv$Edu.size),] %>% group_by(symTreat,symClade,timepoint) %>% summarise(n=n())
larv[!is.na(larv$sym.count),] %>% group_by(symTreat,symClade,timepoint) %>% summarise(n=n())
larv[duplicated(larv$StackID),]
# no duplicates
# filter out damaged and meta samples
larv2=larv[!larv$filter=='Y',]
sample_size2 = larv2 %>% group_by(symTreat,symClade,timepoint) %>% summarise(n=n())
sample_size3 = cbind(sample_size,'n2'=sample_size2$n)
larv2[!is.na(larv$Edu.size),] %>% group_by(symTreat,symClade,timepoint) %>% summarise(n=n())
larv2[!is.na(larv$sym.count),] %>% group_by(symTreat,symClade,timepoint) %>% summarise(n=n())

larv=larv2
#### check distruibutions ####
trait=larv$Edu.size #change to trait of interest
# first check normality, outliers, etc
hist(trait)
boxplot(trait)
outliers=boxplot(trait)$out
larv[trait %in% outliers,] 
# size: 14J (some damage but size still accurate), 32G
# Edu: 9 outliers
# Edu/size: 6 outliers (most look same as Edu) 
# sym: 20 outliers
# sym.size: 21 outliers
library(e1071)
skewness(na.omit(trait))
# size: 0.2378297
# Edu: 1.171549
# Edu.size: 1.718824
# sym: 2.237483
# sym.size: 2.44853
shapiro.test(na.omit(trait))
# size: normal
# Edu: not normal
# Edu.size: not normal
# sym: not normal
# sym.size: not normal
library(car)
qqPlot(na.omit(trait))

# checked outliers, no obvious issues
# size is normal, but both Edu and sym counts not normal
# could normalize but also need to do batch correction and feel like that should be done first

#### plots #####
library(ggplot2)
# checking for size -- can only use ones where we could measure full area
ggplot(larv[!larv$symClade %in% c('C') & larv$offscreen %in% c('no'),],aes(x=timepoint,y=size,color=symTreat))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(color=symTreat),position=position_jitterdodge(),size=2)+
  scale_color_manual(values=c("#5BBCD6","#F98400"))+
  facet_wrap(~symClade,nrow=1)+
  ylab('area')+
  theme_bw(base_size=20)

# edu
ggplot(larv[!larv$symClade %in% c('C'),],aes(x=timepoint,y=Edu,color=symTreat))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(color=symTreat),position=position_jitterdodge(),size=2)+
  scale_color_manual(values=c("#5BBCD6","#F98400"))+
  facet_wrap(~symClade,nrow=1)+
  ylab('Edu/larvae')+
  theme_bw(base_size=20)

# syms per larvae
ggplot(larv[!larv$symClade %in% c('C'),],aes(x=timepoint,y=sym.count,color=symTreat))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(color=symTreat),position=position_jitterdodge(),size=2)+
  scale_color_manual(values=c("#5BBCD6","#F98400"))+
  facet_wrap(~symClade,nrow=1)+
  ylab('syms/larvae')+
  theme_bw(base_size=20)

# syms per area
ggplot(larv[!larv$symClade %in% c('C'),],aes(x=timepoint,y=sym.size,color=symTreat))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(color=symTreat),position=position_jitterdodge(),size=2)+
  scale_color_manual(values=c("#5BBCD6","#F98400"))+
  facet_wrap(~symClade,nrow=1)+
  ylab('syms/area')+
  theme_bw(base_size=20)

##### now need to normalize for batch effect between new and old confocal ####
# size should be consistent except that some were cutoff in the old images
# Edu definitely needs to be corrected for, will need to check bias in syms
# try batch effect correction following https://cran.r-project.org/web/packages/batchtma/vignettes/batchtma.html
#remotes::install_github("stopsack/batchtma")
library(batchtma)
library(tidyverse)
larv %>% plot_batch(marker = Edu.size, batch = confocal, color = symTreat)
larv %>% plot_batch(marker = sym.count, batch = confocal, color = symTreat) #syms don't have batch effect, yay!
#method = simple calculates the mean for each batch and subtracts the difference between this mean and the grand mean, such that all batches end up having a mean equivalent to the grand mean. 
#Differences in variance between batches will remain, if they exist.
larv %>% 
  adjust_batch(markers = Edu.size, batch = confocal, 
               method = simple) %>%
  plot_batch(marker = Edu.size_adj2, batch = confocal, color = symTreat)
# method = quantnorm performs quantile normalization: values are ranked within each batch, and then each rank is assigned the mean per rank across batches. Quantile normalization ensures that all batches have near-identical biomarker distributions. 
# However, quantile normalization does not allow for accounting for confounders.
larv %>% 
  adjust_batch(markers = Edu.size, batch = confocal, 
               method = quantnorm) %>%
  plot_batch(marker = Edu.size_adj6, batch = confocal, color = symTreat)
# this accounts for mean and variance
# make new columns with batch corrected values
larv_adj <- larv %>% 
  adjust_batch(markers = Edu.size, batch = confocal, 
               method = simple) %>%
  adjust_batch(markers = Edu.size, batch = confocal, 
               method = quantnorm)
# rename columns to batch correction applied
colnames(larv_adj)[18:19]=c('Edu.size.simple','Edu.size.quantnorm')

# let's see how this affects results
ggplot(larv_adj[!larv_adj$symClade %in% c('C'),],aes(x=timepoint,y=Edu.size.simple,color=symTreat))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(color=symTreat),position=position_jitterdodge(),size=2)+
  scale_color_manual(values=c("#5BBCD6","#F98400"))+
  facet_wrap(~symClade,nrow=1)+
  theme_bw(base_size=20)
ggplot(larv_adj[!larv_adj$symClade %in% c('C'),],aes(x=timepoint,y=Edu.size.quantnorm,color=symTreat))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(color=symTreat),position=position_jitterdodge(),size=2)+
  scale_color_manual(values=c("#5BBCD6","#F98400"))+
  facet_wrap(~symClade,nrow=1)+
  theme_bw(base_size=20)
# both normalizations look like same story! that is great
# a little concerned about smaller sample sizes after quality control

# what if we just batch correct edu (not normalized value)
larv %>% 
  adjust_batch(markers = Edu, batch = confocal, 
               method = simple) %>%
  plot_batch(marker = Edu_adj2, batch = confocal, color = symTreat)
larv %>% 
  adjust_batch(markers = Edu, batch = confocal, 
               method = quantnorm) %>%
  plot_batch(marker = Edu_adj6, batch = confocal, color = symTreat)
# make new columns with batch corrected values
larv_adj2 <- larv_adj %>% 
  adjust_batch(markers = Edu, batch = confocal, 
               method = simple) %>%
  adjust_batch(markers = Edu, batch = confocal, 
               method = quantnorm)
# rename columns to batch correction applied
colnames(larv_adj2)[20:21]=c('Edu.simple','Edu.quantnorm')
# let's see how this affects results
ggplot(larv_adj2[!larv_adj2$symClade %in% c('C'),],aes(x=timepoint,y=Edu.simple/size,color=symTreat))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(color=symTreat),position=position_jitterdodge(),size=2)+
  scale_color_manual(values=c("#5BBCD6","#F98400"))+
  facet_wrap(~symClade,nrow=1)+
  theme_bw(base_size=20)
ggplot(larv_adj2[!larv_adj2$symClade %in% c('C'),],aes(x=timepoint,y=Edu.quantnorm/size,color=symTreat))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(color=symTreat),position=position_jitterdodge(),size=2)+
  scale_color_manual(values=c("#5BBCD6","#F98400"))+
  facet_wrap(~symClade,nrow=1)+
  theme_bw(base_size=20)


#### need to normalize symbiont counts ####
hist(larv$sym.count)
hist(log10(larv$sym.count))
larv_adj$sym.count.log=log10(larv_adj$sym.count)
larv_adj$sym.count.log[larv_adj$sym.count.log=='-Inf']<-NA # need to set apos to NA because log(0) is -Inf
skewness(na.omit(larv$sym.count))
skewness(na.omit(larv_adj$sym.count.log))
shapiro.test(na.omit(larv$sym.count))
shapiro.test(na.omit(larv_adj$sym.count.log))
qqPlot(larv_adj$sym.count.log)
# not exactly normal but much closer, probs good enough for mixed models

# visualize log transform
ggplot(larv_adj[!larv_adj$symClade %in% c('C','apo'),],aes(x=timepoint,y=sym.count.log,color=symTreat))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(color=symTreat),position=position_jitterdodge(),size=2)+
  scale_color_manual(values=c("#5BBCD6","#F98400"))+
  facet_wrap(~symClade,nrow=1)+
  theme_bw(base_size=20)
ggplot(larv_adj[!larv_adj$symClade %in% c('C','apo'),],aes(x=timepoint,y=sym.count.log/size,color=symTreat))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(color=symTreat),position=position_jitterdodge(),size=2)+
  scale_color_manual(values=c("#5BBCD6","#F98400"))+
  facet_wrap(~symClade,nrow=1)+
  theme_bw(base_size=20)

larv_adj2$sym.count.log=log10(larv_adj2$sym.count)
larv_adj2$sym.count.log[larv_adj2$sym.count.log=='-Inf']<-NA

# may also want Edu.size normalization for mixed model
hist(larv$Edu.size)
hist(log10(larv$Edu.size)) # skews other direction
hist(sqrt(larv$Edu.size))
skewness(na.omit(larv$Edu.size))
skewness(na.omit(log10(larv$Edu.size)))
skewness(na.omit(sqrt(larv$Edu.size))) #sqrt is least skewed
shapiro.test(na.omit(larv$Edu.size))
shapiro.test(na.omit(log10(larv$Edu.size)))
shapiro.test(na.omit(sqrt(larv$Edu.size))) # normal!
qqPlot(larv$Edu.size)
qqPlot(log10(larv$Edu.size))
qqPlot(sqrt(larv$Edu.size))

larv_adj2$Edu.size.sqrt=sqrt(larv_adj2$Edu.size)

write.csv(larv_adj2,"data/LarvAdj4Stats_EduSizeSym.csv")


##################
#### Settlers ####
##################
#set=read.csv('data/SetEduSizeSym_wMeta_final.csv')
set=read.csv('data/SetEduSizeSym_wMeta_QC_21Jan2025.csv')
# sym_score with N indicates those should be filtered out, but just for sym score (I think, will need to check Edus)
set[,c("symClade","symTreat")]=lapply(set[,c("symClade","symTreat")],as.factor)
set$symClade=factor(set$symClade,levels=c('apo','A','B','C','D'))
set$symTreat=factor(set$symTreat,levels=c('apo','25','32'))
setQCbyJQ=read.csv('data/RecruitEduQC_JQ3Feb2025.csv')
setQC=merge(set,setQCbyJQ[,c('File','damage','upside.down','thresholding','usable','add.notes')],by='File',all.x=T)
write.csv(setQC,'data/SetEduSizeSym_wMeta_QC_15Jan2026.csv')

set=setQC
set %>% group_by(symClade,symTreat) %>% summarise(n())
set[!set$sym_score=='N',] %>% group_by(symClade,symTreat) %>% summarise(n())
count(set[set$damage=='y',]) # 32
count(set[set$upside.down=='y',]) # 22
count(set[set$thresholding=='b',]) # 9
count(set[set$thresholding=='b'&set$upside.down=='y'&set$damage=='y',]) # 2
set[!set$Edu_score=='not usable',] %>% group_by(symClade,symTreat) %>% summarise(n())
set[!set$upside.down=='y',] %>% group_by(symClade,symTreat) %>% summarise(n())


#### check distruibutions ####
trait=set$sym_area/set$size #change to trait of interest
# first check normality, outliers, etc
hist(trait)
boxplot(trait)
outliers=boxplot(trait)$out
set[trait %in% outliers,] 
# size: no outliers
# Edu: no outliers
# Edu/size: 1 outlier, 120B  
# sym: 5 outliers
# sym.size: 6 outliers
library(e1071)
skewness(na.omit(trait))
# size: 0.2161633
# Edu: 0.30389
# Edu.size: 0.699814
# sym: 1.907911
# sym.size: 1.766217
shapiro.test(na.omit(trait))
# size: normal
# Edu: not normal, but doesn't look so bad
# Edu.size: not normal
# sym: not normal, very skewed
# sym.size: not normal
library(car)
qqPlot(na.omit(trait))

# probably need to transform Edu and sym in sets
# let's try same transformations as larvae
hist(sqrt(set$Edu_area[!set$usable=='N']/set$size[!set$usable=='N'])) # better
hist(log10(set$Edu_area[!set$usable=='N']/set$size[!set$usable=='N'])) # more skewed

hist(set$sym_area[!set$sym_score=='N'&!set$symClade=='apo'])
hist(log10(set$sym_area[!set$sym_score=='N'&!set$symClade=='apo']/set$size[!set$sym_score=='N'&!set$symClade=='apo'])) # worse
hist(log(set$sym_area[!set$sym_score=='N'&!set$symClade=='apo']/set$size[!set$sym_score=='N'&!set$symClade=='apo'])) 
hist(sqrt(set$sym_area[!set$sym_score=='N'&!set$symClade=='apo']/set$size[!set$sym_score=='N'&!set$symClade=='apo'])) # maybe better but still not normal
hist((set$sym_area[!set$sym_score=='N'&!set$symClade=='apo']/set$size[!set$sym_score=='N'&!set$symClade=='apo'])^(1/3)) #this is the best!

ggplot(set[!set$symClade=='apo' & !set$sym_score=='N',],aes(x=symTreat,y=(sym_area/size)^(1/3),color=symTreat))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(color=symTreat),position=position_jitterdodge(),size=2)+
  scale_color_manual(values=c("#5BBCD6","#F98400"))+
  facet_wrap(~symClade,nrow=1)+
  theme_bw(base_size=20)

ggplot(set[!set$symClade=='apo'&!set$Edu_score=='not usable',],aes(x=symTreat,y=sqrt(Edu_area/size),color=symTreat))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(color=symTreat),position=position_jitterdodge(),size=2)+
  scale_color_manual(values=c("#5BBCD6","#F98400"))+
  facet_wrap(~symClade,nrow=1)+
  theme_bw(base_size=20)

ggplot(set[!set$symClade %in% c('apo','C') & !set$sym_score=='N',],aes(x=sqrt(Edu_area/size),y=log10(sym_area/size),color=symClade,shape=symTreat,linetype=symTreat))+
  geom_point(size=2)+
  scale_color_manual(values=c("#e67e22","#f4d03f","#52be80"))+
  geom_smooth(method='lm', formula= y~x,se=F,size=1.5)+
  theme_bw(base_size=20)+
  ylab(expression(log[10]("syms/area")))+
  xlab(expression(sqrt("Edu/area")))


