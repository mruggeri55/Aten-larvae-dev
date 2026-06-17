setwd('/Users/maria/Library/CloudStorage/GoogleDrive-mariaruggeri55@gmail.com/Shared drives/Okinawa coral larvae EdU/glycan expt')

df=read.csv('MQYres - 24Feb2025.csv')
#df=read.csv('EQYmes - 21Feb2025.csv')

df$spp=gsub('[[:digit:]]+|-', '', df$sampleID)
df$culture=sapply(strsplit(df$sampleID,'-'),FUN='[',1)
df$rep=sapply(strsplit(df$sampleID,'-'),FUN='[',2)

library(ggplot2)
ggplot(df,aes(x=date,y=Fv.Fm,color=treatment))+
  geom_boxplot()+
  geom_point(position=position_jitterdodge())

ggplot(df,aes(x=date,y=Fv.Fm,color=treatment))+
  geom_boxplot()+
  geom_point(position=position_jitterdodge())+
  facet_wrap(~spp)

ggplot(df,aes(x=treatment,y=Fv.Fm,color=spp))+
  geom_boxplot()+
  geom_point(position=position_jitterdodge())+
  facet_wrap(~date)


## plot avgs
library(dplyr)
sumdf = df %>% group_by(date,spp,treatment) %>% summarise(avg=mean(Fv.Fm),err=sd(Fv.Fm))
write.csv(sumdf, file="/Users/maria/Documents/work/Weis_lab/Acropora larvae infections/analysis_Nov2024/data/FvFm_sum_2Apr2026.csv")

ggplot(sumdf, aes(x=date,y=avg,color=spp,group=interaction(spp,treatment)))+
  geom_line(aes(linetype=treatment),position=position_dodge(0.4))+
  geom_point(position=position_dodge(0.4))+
  geom_errorbar(aes(ymin=avg-err,ymax=avg+err),position=position_dodge(0.4))

ggplot(sumdf[sumdf$date == '2/20/25',], aes(x=treatment,y=avg,color=spp))+
  geom_point(position=position_dodge(0.4))+
  geom_line(position=position_dodge(0.4),aes(group=spp))+
  geom_errorbar(aes(ymin=avg-err,ymax=avg+err),position=position_dodge(0.4))+
  theme_bw(base_size = 20)+
  ylab('Fv/Fm')
FvFmP=ggplot(sumdf[sumdf$date == '2/20/25',], aes(x=spp,y=avg,color=treatment))+
  geom_point(position=position_dodge(1),size=2)+
  geom_errorbar(aes(ymin=avg-err,ymax=avg+err),position=position_dodge(1),size=1,width=0.5)+
  scale_color_manual(values=c("#5BBCD6","#F98400"))+
  theme_bw(base_size = 20)+
  ylab('Fv/Fm')


#### make new variable for change in Fv/Fm from initial control value

baseline = sumdf[sumdf$treatment %in% c('C','c') & sumdf$date=='2/17/25',]
baseline$baseline = baseline$avg

delFvFm = merge(df,baseline[c('spp','baseline')],by='spp',all.x=T)
delFvFm$delFvFm = delFvFm$Fv.Fm - delFvFm$baseline

ggplot(delFvFm,aes(x=date,y=delFvFm,color=treatment))+
  geom_boxplot()+
  geom_point(position=position_jitterdodge())+
  facet_wrap(~spp)
ggplot(delFvFm,aes(x=treatment,y=delFvFm,color=spp))+
  geom_boxplot()+
  geom_point(position=position_jitterdodge())+
  facet_wrap(~date)



######################################
#### read in symbiont density data #####
#####################################
dens=read.csv('SymQuantFinal_21Feb2025 - Sheet1.csv')
dens_df=dens[,c('Sample','treatment','density')]
dens_df$spp=gsub('[[:digit:]]+|-', '', dens_df$Sample)
dens_df$spp[dens_df$Sample=='BI - 3']<-'B'
dens_df$spp[dens_df$spp=='C ']<-'C'

ggplot(dens_df, aes(x=treatment,y=density,fill=spp))+
  geom_boxplot(outliers=F,position=position_dodge(0.4))+
  geom_point(position=position_dodge(0.4))

dens_df$del_dens=dens_df$density-375000
ggplot(dens_df, aes(x=treatment,y=del_dens,fill=spp))+
  geom_boxplot(outliers=F,position=position_dodge(0.4))+
  geom_point(position=position_dodge(0.4))
ggplot(dens_df, aes(x=spp,y=del_dens,fill=treatment))+
  geom_boxplot(outliers=F,position=position_dodge(1))+
  geom_point(position=position_dodge(1))

sum_delTP = dens_df %>% group_by(spp,treatment) %>% summarise(avg=mean(del_dens),err=sd(del_dens))
delTPp=ggplot(sum_delTP, aes(x=spp,y=avg,color=treatment))+
  geom_point(position=position_dodge(1),size=2)+
  geom_errorbar(aes(ymin=avg-err,ymax=avg+err),position=position_dodge(1),size=1,width=0.5)+
  scale_color_manual(values=c("#5BBCD6","#F98400"))+
  theme_bw(base_size = 20)+
  ylab(expression(atop(Delta * density, "(final - initial)")))+
  scale_x_discrete(labels=c("A"="S.mic","B"="B.min","C"="C.gor","D"="D.tren"))+
  xlab('symbiont species')



# need to compare to control vals I think
sum_dens = dens_df %>% group_by(spp,treatment) %>% summarise(avg=mean(density),err=sd(density))
baseline = sum_dens[sum_dens$treatment %in% c('C','c'),]
baseline$baseline = baseline$avg

delC_dens = merge(dens_df,baseline[c('spp','baseline')],by='spp',all.x=T)
delC_dens$delC_dens = delC_dens$density - delC_dens$baseline

ggplot(delC_dens,aes(x=spp,y=delC_dens,color=treatment))+
  geom_boxplot(outliers=F,position=position_dodge(1))+
  geom_point(position=position_dodge(1))
ggplot(delC_dens,aes(x=treatment,y=delC_dens,color=spp))+
  geom_boxplot(outliers=F,position=position_dodge(1))+
  geom_point(position=position_dodge(1))

sum_del_dens = delC_dens %>% group_by(spp,treatment) %>% 
  summarise(avg=mean(delC_dens),err=sd(delC_dens))

ggplot(sum_del_dens, aes(x=treatment,y=avg,color=spp))+
  geom_point(position=position_dodge(1))+
  geom_line(position=position_dodge(1),aes(group=spp))+
  geom_errorbar(aes(ymin=avg-err,ymax=avg+err),position=position_dodge(1))+
  theme_bw(base_size = 20)

densP=ggplot(sum_del_dens, aes(x=spp,y=avg,color=treatment))+
  geom_point(position=position_dodge(1),size=2)+
  geom_errorbar(aes(ymin=avg-err,ymax=avg+err),position=position_dodge(1),size=1,width=0.5)+
  scale_color_manual(values=c("#5BBCD6","#F98400"))+
  theme_bw(base_size = 20)+
  ylab(expression(paste(Delta,'density')))


# make combined density and FvFm plot
library(cowplot)
legend <- get_legend(FvFmP)
# del density being compared to control avgs
pgrid=plot_grid(FvFmP+theme(axis.title.x = element_blank(),legend.position = 'none',plot.margin = unit(c(0.5,0.5,0,0.5), "cm"),axis.ticks.x=element_blank(),axis.text.x=element_blank()),
             densP+theme(axis.title.x = element_blank(),legend.position = 'none',plot.margin = unit(c(0,0.5,0.5,0.5), "cm")),
          ncol=1,align='v',rel_heights = c(0.85,1))

plot_grid(pgrid,legend,ncol=2,rel_widths = c(1, .3))

# del density being comparing to initial
pgrid2=plot_grid(FvFmP+theme(axis.title.x = element_blank(),legend.position = 'none',plot.margin = unit(c(0.5,0.5,0,0.5), "cm"),axis.ticks.x=element_blank(),axis.text.x=element_blank()),
                delTPp+theme(axis.title.x = element_blank(),legend.position = 'none',plot.margin = unit(c(0,0.5,0.5,0.5), "cm")),
                ncol=1,align='v',rel_heights = c(0.85,1))

plot_grid(pgrid2,legend,ncol=2,rel_widths = c(1, .3))




#################
##### STATS #####
#################
library(lme4)
library(lmerTest)
library(emmeans)
library(sjPlot)
library(piecewiseSEM)
library(car)

# check for normality
hist(df$Fv.Fm)
boxplot(df$Fv.Fm)
library(e1071)
skewness(df$Fv.Fm)
shapiro.test(df$Fv.Fm)
qqPlot(df$Fv.Fm)

library(lubridate)
df$date_time=mdy(df$date)
df$treatment=as.factor(df$treatment)
df$spp=as.factor(df$spp)

allMix=lmer(Fv.Fm ~ treatment*spp + date_time + (1|rep:spp), data=df)
qqPlot(residuals(allMix))
plot(allMix,which=1) 
summary(allMix)

rsquared(allMix)
rand(allMix)
allRes=anova(allMix,ddf='Kenward-Roger') # everything is significant!!!
joint_tests(allMix) # to get sum p vals by treatment, spp, date, and interaction

pairs(emmeans(allMix,c('spp','treatment')),simple = 'each')
# temp has significant effect on all spp but D!
# also all species different in control except A-C
# all different in the heat

# this does not give me pairwise comparisons of interaction terms though
# I think the following will work
contrast(emmeans(allMix,c('spp','treatment')),interaction = 'pairwise')

### now just stats for last timepoint
finalLM=lm(Fv.Fm ~ treatment*spp, data=df[df$date=='2/20/25',])
qqPlot(residuals(finalLM))
plot(finalLM,which=1) 
summary(finalLM)

rsquared(finalLM)
joint_tests(finalLM) # to get sum p vals by treatment, spp, and interaction

emmeans(finalLM,c('spp','treatment'))
pairs(emmeans(finalLM,c('spp','treatment')),simple = 'each')
# temp has significant effect on all spp but D!
# also all species different in control except A-C
# all different in the heat
contrast(emmeans(finalLM,c('spp','treatment')),interaction = 'pairwise')


### density stats
# density compared to control mean (within spp)
hist(delC_dens$delC_dens)
boxplot(delC_dens$delC_dens)
skewness(delC_dens$delC_dens)
shapiro.test(delC_dens$delC_dens)
qqPlot(delC_dens$delC_dens)

# looks normal, moving on to stats
delC_dens$treatment=as.factor(delC_dens$treatment)
delC_dens$spp=as.factor(delC_dens$spp)
densLM=lm(delC_dens ~ treatment*spp, data=delC_dens)
qqPlot(residuals(densLM))
plot(densLM,which=1) 
summary(densLM)

# only treatment significant, and not best R2... only 0.3657
# not surprised no spp effect because normalized to baseline
# but also no interactive effect
# must be too much variation
# but let's try not normalized densities just to check

dens_df$treatment=as.factor(dens_df$treatment)
dens_df$spp=as.factor(dens_df$spp)
dens_rawLM=lm(density ~ treatment*spp, data=dens_df)
qqPlot(residuals(dens_rawLM))
plot(dens_rawLM,which=1) 
summary(dens_rawLM)
# much better R2, but still no interaction
# mainly driven by baseline differences between spp and treatment

# let's try final-initial (although should really change anything because all the same starting)
dens_delTP_LM=lm(del_dens ~ treatment*spp, data=dens_df)
qqPlot(residuals(dens_delTP_LM))
plot(dens_delTP_LM,which=1) 
summary(dens_delTP_LM)
emmeans(dens_delTP_LM,c('spp','treatment'))
pairs(emmeans(dens_delTP_LM,c('spp','treatment')),simple = 'each')
# only A significant change in density due to treatment
