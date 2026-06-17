setwd("/Users/maria/Documents/work/Weis_lab/Acropora larvae infections/analysis_Nov2024")

library(dplyr)
library(lme4)
library(lmerTest)
library(emmeans)
library(sjPlot)
library(piecewiseSEM)
library(car)
library(stargazer)

larv=read.csv("data/LarvAdj4Stats_EduSizeSym.csv",row.names = 1)
larv[,c("symClade","symTreat","timepoint","confocal")]=lapply(larv[,c("symClade","symTreat","timepoint","confocal")],as.factor)
larv$symClade=factor(larv$symClade,levels=c('apo','A','B','C','D'))
larv$symTreat=factor(larv$symTreat,levels=c('apo','25','32'))
str(larv)

##### larvae Edu #####
# first let's use raw Edu by size values and control for confocal as random effect
EduMix=lmer(sqrt(Edu.size) ~ symTreat*symClade*timepoint + (1|confocal), data=larv[!larv$symClade %in% c('C','apo'),])
qqPlot(residuals(EduMix))
plot(EduMix,which=1) # heteroskedastic, but log and sqrt transformation fixes this!
summary(EduMix)
rsquared(EduMix)
rand(EduMix) # significant
tab_model(EduMix, file = "stats/EduMixLarv_results_table.html")
EduRes=anova(EduMix,ddf='Kenward-Roger') # everything is significant!!!
#write.csv(EduRes,'stats/LarvEduLMMres.csv')
EduRand=as.data.frame(rand(EduMix))
tab_dfs(list(EduRes,EduRand),show.rownames=T,digits=4)

emmeans(EduMix,'symTreat')
emmeans(EduMix,'symClade')
pairs(emmeans(EduMix,'symClade'))
emmeans(EduMix,'timepoint')
EduTab=as.data.frame(pairs(emmeans(EduMix,c('symClade','symTreat'))))
write.csv(EduTab,"stats/Larv_EduMixRes_pw_comps.csv")

# check if Edu is significantly lower in apo compared to D at day 9
# subset 32 treatments at day 9 and apos
sub32d9=larv[larv$symTreat=='32'&larv$timepoint=='9',]
sub32d9wApo=rbind(sub32d9,larv[larv$symClade=='apo'&larv$timepoint=='9',])
summary(aov(sqrt(Edu.size)~symClade,sub32d9wApo))
TukeyHSD(aov(sqrt(Edu.size)~symClade,sub32d9wApo))
emmeans(aov(sqrt(Edu.size)~symClade,sub32d9wApo),'symClade')
Edu_apoVspp32Res=as.data.frame(TukeyHSD(aov(sqrt(Edu.size)~symClade,sub32d9wApo))$symClade)
write.csv(Edu_apoVspp32Res,'stats/Larv_Edu_apoVspp32Res.csv')

# yes Edu significantly reduced in 32D compared to apo on day9, also lower than 32A and 32B
ggplot(sub32d9wApo,aes(x=symClade,y=sqrt(Edu.size)))+geom_boxplot()



# Edu=lm(Edu.size.simple ~ symTreat*symClade*timepoint, data=larv[!larv$symClade %in% c('C','apo'),])
# summary(Edu)
# Edu_res=as.data.frame(pairs(emmeans(Edu,c('symTreat','symClade','timepoint'))))
# Edu_res[Edu_res$p.value<0.05,]
# 
# A=lm(Edu.size.quantnorm ~ symTreat*timepoint, data=larv[larv$symClade %in% c('A'),])
# B=lm(Edu.size.quantnorm ~ symTreat*timepoint, data=larv[larv$symClade %in% c('B'),])
# D=lm(Edu.size.quantnorm ~ symTreat*timepoint, data=larv[larv$symClade %in% c('D'),])
# summary(A)
# summary(B)
# summary(D)
# so no interaction between sym treatment and timepoint for A
# but significant interaction for B & D! Aka sym treatment affects host proliferation


##### larvae sym ######
symMix=lmer(sym.count.log/size ~ symTreat*symClade*timepoint + (1|confocal), data=larv[!larv$symClade %in% c('C','apo'),])
# error is singular, apparently confocal explains no variance
ggplot(larv,aes(x=confocal,y=sym.count.log/size))+geom_boxplot()+geom_point()
# sym counts do not appear to be affected by confocal
qqPlot(residuals(symMix))
plot(symMix,which=1) 
summary(symMix)
rsquared(symMix)
rand(symMix) # confocal explains no variance
anova(symMix,ddf='Kenward-Roger') 
tab_model(symMix, file = "stats/symMixLarv_wInt_results_table.html")

# no effect of interactions between timepoint and symtreat or symclade, let's try removing interaction to see if it helps singularity
symMix2=lmer(sym.count.log/size ~ symTreat*symClade + timepoint + (1|confocal), data=larv[!larv$symClade %in% c('C','apo'),])
rand(symMix2) 
# still singular, confocal explaining no variance
# let's try linear regression
symLM=lm(sym.count.log/size ~ symTreat*symClade*timepoint, data=larv[!larv$symClade %in% c('C','apo'),])
summary(symLM)
symLM2=lm(sym.count.log/size ~ symTreat + symClade + timepoint, data=larv[!larv$symClade %in% c('C','apo'),])
summary(symLM2)
tab_model(symLM2, file = "stats/symLMLarv_noInt_results_table.html")

pairs(emmeans(symLM2,c('symClade')))
ggplot(larv[!larv$symClade %in% c('C','apo'),],aes(x=symClade,y=sym.count.log/size))+
  geom_boxplot(outliers = F)+
  geom_point(aes(color=symTreat))
emmeans(symLM2,c('symClade'))

emmeans(symLM2,c('symTreat','symClade'))
larv[!larv$symClade %in% c('C','apo'),] %>% group_by(symClade,symTreat) %>% summarise(mean=mean(sym.count.log/size,na.rm=T))
pairs(emmeans(symLM2,c('symClade','timepoint')))
symLM2pw_allfactors=as.data.frame(pairs(emmeans(symLM2,c('symTreat','symClade','timepoint'))))
write.csv(symLM2pw_allfactors,"stats/Larv_SymLMres_pw_comps_allfactors.csv")

symLM2pw=as.data.frame(pairs(emmeans(symLM2,c('symClade','symTreat'))))
write.csv(symLM2pw,"stats/Larv_SymLMres_pw_comps.csv")

#### larvae host-sym correlation ####
larv$sqrt_Edu_size=sqrt(larv$Edu.size)
cormod=lm(sym.count.log/size~sqrt_Edu_size*symClade*symTreat,data=larv[!larv$symClade %in% c('C','apo'),])
summary(cormod)
anova(cormod)

test(emtrends(cormod, ~ symClade,var = "sqrt_Edu_size"))
test(emtrends(cormod, ~ symTreat,var = "sqrt_Edu_size"))
test(emtrends(cormod, ~ symClade * symTreat,var = "sqrt_Edu_size"))

exp_cormod=as.data.frame(anova(cormod))
exp_slopes=as.data.frame(test(emtrends(cormod, ~ symClade * symTreat,var = "sqrt_Edu_size")))
write.csv(exp_cormod,file='stats/symdensVSedu_correlation_anova.csv')
write.csv(exp_slopes,file='stats/symdensVSedu_correlation_slopes.csv')


-----------------------
########## set #########
-----------------------
set=read.csv('data/SetEduSizeSym_wMeta_QC_15Jan2026.csv')
# sym_score with N indicates those should be filtered out
# usable column is for Edu, N need to be filtered out
set[,c("symClade","symTreat")]=lapply(set[,c("symClade","symTreat")],as.factor)
set$symClade=factor(set$symClade,levels=c('apo','A','B','C','D'))
set$symTreat=factor(set$symTreat,levels=c('apo','25','32'))
set %>% group_by(symClade,symTreat) %>% summarise(n())
set[!set$sym_score=='N'&!set$usable=='N',] %>% group_by(symClade,symTreat) %>% summarise(n())
str(set)

###### Edu
EduMixSet=lmer(sqrt(Edu_area/size) ~ symTreat*symClade + (1|confocal), 
               data=set[!set$symClade %in% c('apo','C')&!set$usable=='N',])
qqPlot(residuals(EduMixSet))
plot(EduMixSet,which=1) 
summary(EduMixSet)
rsquared(EduMixSet) # not great
rand(EduMixSet) # confocal not significant
anova(EduMixSet,ddf='Kenward-Roger') # sym treat and sym clade ns, but interaction is
tab_model(EduMixSet, file = "stats/EduMixSet_results_table.html")

pairs(emmeans(EduMixSet,c('symClade','symTreat')))
# interaction driven by differing responses of A and D to pretreatment
# A pretreatment at 32 increased host proliferation, whereas D pretreatment decreased host proliferation

# what about size? might be driving some of these patterns
SizeMixSet=lmer(size ~ symTreat*symClade + (1|confocal), data=set[!set$symClade %in% c('apo','C')&!set$usable=='N',])
qqPlot(residuals(SizeMixSet))
plot(SizeMixSet,which=1) 
summary(SizeMixSet)
rsquared(SizeMixSet) # not great
rand(SizeMixSet) # confocal not significant
anova(SizeMixSet,ddf='Kenward-Roger') # sym treat and sym clade ns, but interaction is
tab_model(SizeMixSet, file = "stats/SizeMixSet_results_table.html")

pairs(emmeans(SizeMixSet,c('symClade','symTreat')))
# D and A significantly different size at 25

t.test(data=set[set$symClade=='A',],size~symTreat)
# with t test, size different between A at 25 and 32 so could be driving diffs in Edu
t.test(data=set[set$symClade=='A',],Edu_area~symTreat)
# but unnormalized Edu also sig different, so both probably factoring in

##### sym
SymMixSet=lmer((sym_area/size)^(1/3) ~ symTreat*symClade + (1|confocal), data=set[!set$symClade %in% c('apo','C')&!set$sym_score=="N",])
qqPlot(residuals(SymMixSet))
plot(SymMixSet,which=1) 
summary(SymMixSet)
rsquared(SymMixSet) # TERRIBLE
rand(SymMixSet) # confocal not significant
anova(SymMixSet,ddf='Kenward-Roger') # nothing significant and model is singular
tab_model(SymMixSet, file = "stats/SymMixSet_results_table.html")

# try lm removing confocal
SymLMSet=lm((sym_area/size)^(1/3) ~ symTreat*symClade, data=set[!set$symClade %in% c('apo','C')&!set$sym_score=="N",])
qqPlot(residuals(SymLMSet))
plot(SymLMSet,which=1) 
summary(SymLMSet) #nothing significant, even without interaction


#### Hoechst to size correlation ####
HS=read.csv("/Users/maria/Documents/work/Weis_lab/Acropora larvae infections/analysis_Nov2024/data/size-hoechst correlation - Sheet1.csv")
ggplot(HS,aes(x=hoechst,y=size))+geom_point(size=2)+geom_smooth(method='lm', formula= y~x)+
  ylab("total area")+xlab("Hoechst area")+
  theme_bw(base_size=20)
summary(lm(size~hoechst,data=HS))

#### infection rate ####
inf=read.csv("/Users/maria/Documents/work/Weis_lab/Acropora larvae infections/analysis_Nov2024/data/inf_rate.csv",row.names = 1)
str(inf)
inf[,c('symClade','hostStage','symTreat')]=lapply(inf[,c('symClade','hostStage','symTreat')],factor)

# note this is binary (sym/tot) so will analyze differently
summary(glm(prop_sym~symTreat*symClade,data=inf,weights=tot,family="binomial"))
summary(glm(prop_sym~symTreat*symClade,
            data=inf[!inf$hostStage=='set'&!inf$symClade=='apo',],
            weights=tot,family="binomial"))
summary(glm(prop_sym~symTreat*symClade+timepoint,
            data=inf[!inf$hostStage=='set'&!inf$symClade=='apo',],
            weights=tot,family="binomial"))
summary(glm(prop_sym~symTreat*symClade+hostStage,
            data=inf[!inf$symClade=='apo',],
            weights=tot,family="binomial"))

# no effect of thermal treatment
# effect of symbiont -- C different from others

GLM=glm(prop_sym~symTreat+symClade+hostStage+timepoint,
        data=inf[!inf$symClade=='apo',],
        weights=tot,family="binomial")
joint_tests(GLM)

# most recent analysis below
# might be better to test the effect of life stage only using 9dpi
summary(glm(prop_sym~symTreat+symClade+hostStage,
            data=inf[!inf$symClade=='apo'&inf$timepoint=='9',],
            weights=tot,family="binomial"))
# just looking at TP9, there is an effect of thermal treatment and life stage on inf rate
# effect of sym spp but again driven by C
# try without C to see if still an effect of thermal treatment
summary(glm(prop_sym~symTreat+symClade+hostStage,
            data=inf[!inf$symClade %in% c('apo','C') & inf$timepoint=='9',],
            weights=tot,family="binomial"))
# no effect of treatment after removing C so previous model was driven by absence of C at 32
# just C
summary(glm(prop_sym~symTreat+hostStage,
            data=inf[inf$symClade %in% c('C'),],
            weights=tot,family="binomial"))
# significant effect of treatment and life stage
# and time without life stage
summary(glm(prop_sym~symTreat+symClade+timepoint,
            data=inf[!inf$symClade=='apo'&inf$hostStage=='larv',],
            weights=tot,family="binomial"))
mod=glm(prop_sym~symTreat+symClade+timepoint,
    data=inf[!inf$symClade=='apo'&inf$hostStage=='larv',],
    weights=tot,family="binomial")
emmeans(mod,pairwise~symClade,type = "response")
tab_model(mod, file = "stats/Inf_results_table_larvae.html")
# so just looking at larvae, no effect of treatment, effect of sym spp driven by C

inf %>% group_by(symClade,hostStage) %>% summarise(mean(prop_sym))
inf[!inf$symClade=='apo',] %>% group_by(hostStage) %>% summarise(mean(prop_sym))
