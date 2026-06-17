setwd("/Users/maria/Documents/work/Weis_lab/Acropora larvae infections/analysis_Nov2024")

##### LARVAE ######
larv=read.csv("data/LarvAdj4Stats_EduSizeSym.csv",row.names = 1)
larv[,c("symClade","symTreat","timepoint","confocal")]=lapply(larv[,c("symClade","symTreat","timepoint","confocal")],as.factor)
larv$symClade=factor(larv$symClade,levels=c('apo','A','B','C','D'))
larv$symTreat=factor(larv$symTreat,levels=c('apo','25','32'))
str(larv)

library(ggplot2)
##### plot Edu
p=ggplot(larv[!larv$symClade %in% c('C'),],aes(x=timepoint,y=sqrt(Edu.size),color=symTreat))+
  geom_boxplot(outlier.shape = NA,size=1)+
  geom_point(aes(color=symTreat),position=position_jitterdodge(),size=2.5)+
  scale_color_manual(values=c("#5BBCD6","#F98400"))+
  facet_wrap(~symClade,nrow=1)+
  theme_bw(base_size=20)+
  ylab(expression(sqrt("Edu/size")))+
  guides(color=guide_legend(title="symbiont\npre-treatment"))

# save legend separately
library(cowplot)
legend <- get_legend(p)
library(ggpubr)
as_ggplot(legend)

# save without legend for better sizing
p+theme(legend.position = 'none')

#### plot syms
p2=ggplot(larv[!larv$symClade %in% c('C'),],aes(x=timepoint,y=sym.count.log/size,color=symTreat))+
  geom_boxplot(outlier.shape = NA,size=1)+
  geom_point(aes(color=symTreat),position=position_jitterdodge(),size=2.5)+
  scale_color_manual(values=c("#5BBCD6","#F98400"))+
  facet_wrap(~symClade,nrow=1)+
  theme_bw(base_size=20)+
  ylab(expression(log[10]("syms/area")))+
  guides(color=guide_legend(title="symbiont\npre-treatment"))+
  scale_y_continuous(labels = ~ sprintf(fmt = "%0.00e", .))

p2+theme(legend.position = 'none')


### plot Edu against sym density
ggplot(larv[!larv$symClade %in% c('C'),],aes(x=sqrt(Edu.size),y=sym.count.log/size,color=symTreat))+
  geom_point()+
  scale_color_manual(values=c("#5BBCD6","#F98400"))+
  geom_smooth(method='lm', formula= y~x,se=F)

ggplot(larv[!larv$symClade %in% c('C','apo'),],aes(x=sqrt(Edu.size),y=sym.count.log/size,color=symClade,shape=symTreat,linetype=symTreat))+
  geom_point(size=2)+
  scale_color_manual(values = c("B"="#0571B0","A"="#4DAC26","C"="gold","D"="#E66101"),
                     labels = c("B"="B.min","A"="S.mic","C"="C.gor","D"="D.tren")) +
  geom_smooth(method='lm', formula= y~x,se=F,size=1.5)+
  theme_bw(base_size=20)+
  ylab(expression(log[10]("syms/area")))+
  xlab(expression(sqrt("Edu/area")))+
  labs(color="Symbiont spp",linetype="pre-treatment",shape="pre-treatment")


#### SETTLERS ######
set=read.csv('data/SetEduSizeSym_wMeta_final.csv')
# sym_score with N indicates those should be filtered out, but just for sym score (I think, will need to check Edus)
set[,c("symClade","symTreat")]=lapply(set[,c("symClade","symTreat")],as.factor)
set$symClade=factor(set$symClade,levels=c('apo','A','B','C','D'))
set$symTreat=factor(set$symTreat,levels=c('apo','25','32'))

ggplot(set,
       aes(x=symTreat,y=sqrt(Edu_area/size),color=symTreat))+
  geom_boxplot(outlier.shape = NA,size=1)+
  geom_point(aes(color=symTreat),position=position_jitterdodge(),size=2.5)+
  scale_color_manual(values=c("#5BBCD6","#F98400"))+
  facet_wrap(~symClade,nrow=1)+
  theme_bw(base_size=20)+
  ylab(expression(sqrt("Edu/size")))+
  guides(color=guide_legend(title="symbiont\npre-treatment"))+
  theme(axis.title.x = element_blank())
ggplot(set[!set$symClade=='apo' & !set$sym_score=='N',],
       aes(x=symTreat,y=log10(sym_area/size),color=symTreat))+
  geom_boxplot(outlier.shape = NA,size=1)+
  geom_point(aes(color=symTreat),position=position_jitterdodge(),size=2.5)+
  scale_color_manual(values=c("#5BBCD6","#F98400"))+
  facet_wrap(~symClade,nrow=1)+
  theme_bw(base_size=20)+
  ylab(expression(log10("symbiont area")))+
  guides(color=guide_legend(title="symbiont\npre-treatment"))+
  theme(axis.title.x = element_blank())

ggplot(set[!set$symClade=='apo' & !set$sym_score=='N',],
       aes(x=sqrt(Edu_area/size),y=log10(sym_area/size),color=symClade,shape=symTreat,linetype=symTreat))+
  geom_point(size=2)+
  geom_smooth(method='lm', formula= y~x,se=F,size=1.5)+
  theme_bw(base_size=18)+
  ggtitle('recruits')+
  scale_color_manual(values=c("#e67e22","#f4d03f","skyblue3","#52be80"))

##### uptake rate ####
inf_rate=read.csv("/Users/maria/Documents/work/Weis_lab/Acropora larvae infections/analysis_Nov2024/data/inf_rate_sum_2Apr2026.csv",
                  row.names = 1)
inf_rate$symTreat=as.factor(inf_rate$symTreat)

infP=ggplot(inf_rate[!inf_rate$symClade=='apo' & !inf_rate$hostStage=='set',],aes(x=symClade,y=prop_sym,fill=symTreat))+
  geom_bar(stat="identity",position=position_dodge())+
  scale_fill_manual(values = c('deepskyblue','darkorange'),
                    name = "symbiont\npre-treatment",
                    labels = c("25"="25°C","32"="32°C"))+
  facet_wrap(~timepoint, labeller=labeller(timepoint=c("3"="3 dpi","9"="9 dpi")))+
  theme_bw(base_size=16)+
  ylab('symbiont\nuptake rate')+
  xlab('species')+
  scale_x_discrete(labels=c("A"="S.mic","B"="B.min","C"="C.gor","D"="D.tren"))+
  theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=0.5),legend.title = element_text(size = 14))

##### Fv/Fm #######
FvFm=read.csv("/Users/maria/Documents/work/Weis_lab/Acropora larvae infections/analysis_Nov2024/data/FvFm_sum_2Apr2026.csv")

FvFmP=ggplot(FvFm[FvFm$date == '2/20/25',], aes(x=spp,y=avg,color=treatment))+
  geom_point(position=position_dodge(1),size=2)+
  geom_errorbar(aes(ymin=avg-err,ymax=avg+err),position=position_dodge(1),size=1,width=0.5)+
  scale_color_manual(values=c("#5BBCD6","#F98400"))+
  theme_bw(base_size = 16)+
  ylab('Fv/Fm')+
  xlab('species')+
  scale_x_discrete(labels=c("A"="S.mic","B"="B.min","C"="C.gor","D"="D.tren"))+
  theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=0.5))+
  ylim(0.5,0.77)

# combine infection rate and FvFm plots
library(cowplot)
legend <- get_legend(infP)
plot_grid(FvFmP+theme(legend.position = 'none',plot.title.position = "plot", plot.title = element_text(hjust = 0, size=14))+ggtitle("A) Pre-inoculation"),
          infP+theme(legend.position = 'none',plot.title.position = "plot",plot.title = element_text(hjust = 0, size=14))+ggtitle("B) Post-inoculation"),
          legend,
          nrow=1,rel_widths = c(1.6, 2,0.8))


