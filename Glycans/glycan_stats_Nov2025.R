setwd("/Users/maria/Documents/work/Weis_lab/Acropora larvae infections")

##### stats for glycan data
# Nov 4 2025
# using limma for diffferential abundance of individual glycans
# using permanova and diversity stats for glycan composition from relative abundance

# read in relevant files
common_corrected=read.csv("glycans/all_data/MR_norm/common_33glycans_PQN_log2_batch_4Nov2025.csv",row.names = 1)
common_corrected_noC=read.csv("glycans/all_data/MR_norm/common_33glycans_noC_PQN_log2_batch_11Nov2025.csv",row.names = 1)
batch2=read.csv("glycans/all_data/MR_norm/batch2only_PQN_log2_4Nov2025.csv",row.names = 1)
batch2_noC1=read.csv("glycans/all_data/MR_norm/batch2only_noC1_PQN_log2_25Nov2025.csv",row.names = 1)
meta_all=read.csv("glycans/all_data/glycan_meta.csv")
glycan_meta=read.csv('glycans/all_data/glycan_info.csv',row.names = 1)

#### limma ####
# change these settings
prefix = "batch2_noC1" # set prefix for naming output
df = batch2_noC1 # set to df of interest
#meta  = meta_all # for common_corrected
#meta = meta_all[meta_all$batch=="B",] # for batch2
#meta = meta_all[!meta_all$spp=="Cgor",] # for noC
meta = meta_all[meta_all$sample %in% colnames(df),]
##############

# then run the following
meta$SppTreat=paste(meta$spp,meta$treat,sep='_')
# check sample order aligns
colnames(df)==meta$sample

# then limma for stats
library(limma)
mod=model.matrix(~ 0 + spp + treat,data=meta)
colnames(mod)
fit <- lmFit(df, mod)

spp_levels <- unique(meta$spp)
pairs <- combn(spp_levels, 2, simplify = FALSE)
contrast_strings <- sapply(pairs, function(x)
  paste0("spp", x[2], " - spp", x[1])
)
spp.contrast.matrix <- makeContrasts(contrasts = contrast_strings, levels = mod)
fit2 <- contrasts.fit(fit, spp.contrast.matrix)
fit2 <- eBayes(fit2)

library(dplyr)
pvals = as.data.frame(fit2$p.value) %>% filter(if_any(everything(), ~ .x < 0.05))
colSums(pvals<0.05)
binary <- ifelse(pvals > 0.05, 1, 0)
library(UpSetR)
upset(as.data.frame(binary),sets=colnames(binary))

# get toptable with adjusted pvalues for each contrast
results_list <- lapply(contrast_strings, function(cn) {
  topTable(fit2, coef = cn, number = Inf, adjust.method = "BH")
})
names(results_list) <- contrast_strings
results_all <- bind_rows(lapply(names(results_list), function(cn) {
  df <- results_list[[cn]]
  df$Contrast <- cn
  df
}), .id = NULL)
results_all$Glycan=sapply(strsplit(rownames(results_all),"[...]"),FUN='[',1)
results_all[results_all$adj.P.Val<0.05,]

saveRDS(results_list,file=paste("glycans/all_data/stats/limma_out_",prefix,".RDS",sep=''))

# check for treatment effect
fit3 <- contrasts.fit(fit, coefficients='treatH')
fit3 <- eBayes(fit3)
treat_res=topTable(fit3, coef = 'treatH', number = Inf, adjust.method = "BH")
# nothing significant
write.csv(treat_res,paste("glycans/all_data/stats/limma_out_treat",prefix,".csv",sep=''))

##### compare output from datasets #####
common_out=readRDS('glycans/all_data/stats/limma_out_common_corrected.RDS')
noC_out=readRDS('glycans/all_data/stats/limma_out_noC.RDS')
batch2_out=readRDS('glycans/all_data/stats/limma_out_batch2.RDS')
batch2_noC1=readRDS('glycans/all_data/stats/limma_out_batch2_noC1.RDS')

gCommon=rownames(common_out$`sppSmic - sppDtren`[common_out$`sppSmic - sppDtren`$adj.P.Val<0.05,])
gNoC=rownames(noC_out$`sppSmic - sppDtren`[noC_out$`sppSmic - sppDtren`$adj.P.Val<0.05,])
gBatch2=rownames(batch2_out$`sppSmic - sppDtren`[batch2_out$`sppSmic - sppDtren`$adj.P.Val<0.05,])

library(gplots)
list=list('common'=gCommon,'noC'=gNoC,'batch2'=gBatch2)
venn(list)

# just for one df looking at overlap in contrasts
DtrenBmin=rownames(noC_out$`sppDtren - sppBmin`[noC_out$`sppDtren - sppBmin`$adj.P.Val<0.05,])
SmicBmin=rownames(noC_out$`sppSmic - sppBmin`[noC_out$`sppSmic - sppBmin`$adj.P.Val<0.05,])
SmicDtren=rownames(noC_out$`sppSmic - sppDtren`[noC_out$`sppSmic - sppDtren`$adj.P.Val<0.05,])
list=list('DtrenBmin'=DtrenBmin,'SmicBmin'=SmicBmin,'SmicDtren'=SmicDtren)
venn(list)

# count number of differenitally abundant glycans for each species comparison
l=c()
for (i in names(batch2_noC1)){
  df=batch2_noC1[[i]]
  l=c(l,length(df$logFC[df$adj.P.Val<0.05]))
}
names(l)<-names(batch2_noC1)

#### glycan abundance plots ####
# convert tables to long format and add glycan id
df=batch2_noC1
prefix="batch2_noC1"
abund=read.csv("glycans/all_data/MR_norm/batch2only_noC1_PQN_log2_25Nov2025.csv",row.names = 1)

long_df <- bind_rows(df, .id = "comp")
long_df$Glycan=sapply(strsplit(rownames(long_df),'[...]'),FUN='[',1)
# make list of significant glycans and subset those from abundance df
#sig_glycans=unique(c(DtrenBmin,SmicBmin,SmicDtren)) 
sig_glycans=unique(long_df$Glycan[long_df$adj.P.Val<0.05])
#common_corrected_noC=read.csv("glycans/all_data/MR_norm/common_33glycans_noC_PQN_log2_batch_11Nov2025.csv",row.names = 1)
#abund=read.csv("glycans/all_data/MR_norm/batch2only_noC1_PQN_log2_25Nov2025.csv",row.names = 1)
#noC_abund_sub_sig=common_corrected_noC[rownames(common_corrected_noC) %in% sig_glycans,]
abund_sub_sig=abund[rownames(abund) %in% sig_glycans,]

# convert to long format and add in sample and glycan metadata
library(reshape2)
#noC_abund_sub_sig$Glycan=rownames(noC_abund_sub_sig)
abund_sub_sig$Glycan=rownames(abund_sub_sig)
#long=melt(noC_abund_sub_sig,id.vars='Glycan',variable.name="symSppTreatRep",value.name = "log2abund")
long=melt(abund_sub_sig,id.vars='Glycan',variable.name="symSppTreatRep",value.name = "log2abund")
long$symSpp=sapply(strsplit(as.character(long$symSppTreatRep),'_'),FUN='[',1)
long$symTreat=sapply(strsplit(as.character(long$symSppTreatRep),'_'),FUN='[',2)
long$SppTreat=paste(long$symSpp,long$symTreat)
long_merge=merge(long,glycan_meta,by='Glycan',all.x=T)

plist=vector("list", length(unique(long_merge$Glycan.Types)))
for (glycan in unique(long_merge$Glycan.Types)){
p<-ggplot(long_merge[long_merge$Glycan.Types == glycan,],
       aes(x=symSpp,y=log2abund,fill=symSpp))+
  geom_boxplot(outliers=F)+
  geom_point()+
  scale_x_discrete(labels=c('Bmin','Dtren','Smic',"Cgor"))+ #remove Cgor if doing subset
  facet_wrap(~Glycan)+
  theme(legend.position = 'none')+
  theme_bw()+
  ggtitle(glycan)
plist[[glycan]] <- p
}
plist

plist.family=vector("list", length(unique(long_merge$Glycan.family)))
for (glycan in unique(long_merge$Glycan.family)){
  p<-ggplot(long_merge[long_merge$Glycan.family == glycan,],
            aes(x=symSpp,y=log2abund,fill=symSpp))+
    geom_boxplot(outliers=F)+
    geom_point()+
    scale_x_discrete(labels=c('Bmin','Dtren','Smic','Cgor'))+ #remove Cgor if doing subset
    facet_wrap(~Glycan)+
    theme(legend.position = 'none')+
    theme_bw()+
    ggtitle(glycan)
  plist.family[[glycan]] <- p
}
plist.family


#### heatmap all sig glycans
library(pheatmap)
#centered_matrix <- as.data.frame(t(apply(noC_abund_sub_sig[,!colnames(noC_abund_sub_sig)=="Glycan"], 1, scale, center = TRUE, scale = FALSE)))
#colnames(centered_matrix)<-colnames(noC_abund_sub_sig[,!colnames(noC_abund_sub_sig)=="Glycan"])
centered_matrix <- as.data.frame(t(apply(abund_sub_sig[,!colnames(abund_sub_sig)=="Glycan"], 1, scale, center = TRUE, scale = FALSE)))
colnames(centered_matrix)<-colnames(abund_sub_sig[,!colnames(abund_sub_sig)=="Glycan"])

#col <- colorRampPalette(c("blue", "white", "red"),bias=1.05)(50)
col <- colorRampPalette(c("#045A8D", "white", "#E07700"),bias=1.05)(50)
Glycan.Type.df=glycan_meta[glycan_meta$Glycan %in% sig_glycans,]
#Glycan.Type.df$Glycan==noC_abund_sub_sig$Glycan
Glycan.Type.df$Glycan==abund_sub_sig$Glycan
#rownames(Glycan.Type.df)<-rownames(noC_abund_sub_sig)
rownames(Glycan.Type.df)<-rownames(abund_sub_sig)
meta = meta_all[meta_all$sample %in% colnames(abund_sub_sig),]
meta$sample==colnames(centered_matrix)
rownames(meta)<-meta$sample

ann_colors <- list(
  spp = c(Cgor = "#7F3C8D", Smic = "forestgreen",Dtren = "firebrick3",Bmin = "#F2B701"),
  treat = c(C = "#5BBCD6", H = "#F98400"),
  # Glycan.Types = c("Paucimannose"        = "#5EBF5A",  # bright green
  #                  "High Mannose"        = "#2E8B3E",  # deeper green
  #                  
  #                  # Base categories with increased hue separation
  #                  "Fucosylated"         = "#F2A529",  # vivid amber
  #                  "Sialylated"          = "#2BA8E8",  # cyan-leaning blue
  #                  "NeuGc"               = "#C556C9",  # saturated magenta-violet
  #                  
  #                  # Hybrids: LCH blends of new bases
  #                  "Sialofucosylated"    = "#6BB9A7",  # teal-turquoise blend (amber ↔ cyan)
  #                  "Fucosylated & NeuGc" = "#E27A84",  # coral-rose blend (amber ↔ magenta)
  #                  
  #                  "Others"              = "#8E8E8E")
  Glycan.Types = c("High Mannose"="#0571B0","Paucimannose"="lightblue","Sialylated"="#E66101","Sialofucosylated"= "firebrick3", "Fucosylated"="#4DAC26",
                   "Fucosylated & NeuGc" = "darkgreen", "NeuGc" = "#7F3C8D","Others"="gold")
)

x11()
pheatmap(centered_matrix,
         color=col,
         #annotation_col=meta[,c("spp","treat","batch")],
         annotation_row=Glycan.Type.df[,c("Glycan.Types"),drop = FALSE],
         annotation_colors = ann_colors,
         labels_col = c('Bmin_C2','Bmin_C3','Bmin_H2','Bmin_H3','Dtren_C2','Dtren_C3','Dtren_H2','Dtren_H3',
                        'Smic_C2','Smic_C3','Smic_H2','Smic_H3','Cgor_C2','Cgor_C3','Cgor_H2','Cgor_H3'),
         show_rownames = F,
         cluster_cols = F,
         gaps_col = c(4,8,12),
         border_color = 'black',
         annotation_names_row = F,
         fontsize = 16)
dev.copy2pdf(file = paste("glycans/all_data/plots/heatmap_",length(sig_glycans),"glycans_diffAbundBySpp_",prefix,".pdf",sep=''))
dev.off()


##### glycan composition #####
rel_all=read.csv("~/Documents/work/Weis_lab/Acropora larvae infections/glycans/all_data/MR_norm/all_87glycans_rel_abund_4Nov2025.csv", row.names = 1,check.names=F)
#rel_all=read.csv('/Users/maria/Documents/work/Weis_lab/Acropora larvae infections/glycans/all_data/MR_norm/common_33glycans_rel_abund_4Nov2025.csv',
#                   row.names=1,check.names=F)

t_rel_all=as.data.frame(t(rel_all))

meta_all=read.csv("glycans/all_data/glycan_meta.csv", row.names=1, check.names = FALSE, stringsAsFactors = FALSE)
meta_all$SppTreat=paste(meta_all$spp, meta_all$treat)

glycan_meta=read.csv('glycans/all_data/glycan_info.csv',row.names = 1)
rownames(glycan_meta)=glycan_meta$Glycan
glycan_sub=glycan_meta[,c("Glycan.family","Glycan.Types")]

library(phyloseq)
OTU=otu_table(as.matrix(rel_all),taxa_are_rows = F) # otu table will be our glycan relative abundance
TAX=tax_table(as.matrix(glycan_sub)) # taxa are glycans
sampledata = sample_data(meta_all)
ps.rel=phyloseq(OTU,TAX,sampledata)
#saveRDS(ps.rel,'/Users/maria/Documents/work/Weis_lab/Acropora larvae infections/glycans/all_data/MR_norm/ps.rel.RDS')

# now using phyloseq to summarise across glycan ID levels
type.rel <- tax_glom(ps.rel, taxrank = "Glycan.Types" )
fam.rel <- tax_glom(ps.rel, taxrank = "Glycan.family")

#### beta diversity 
library(vegan)
bet.ps.spp <- betadisper(vegdist(ps.rel@otu_table),ps.rel@sam_data$spp)
bet.ps.treat <- betadisper(vegdist(ps.rel@otu_table),ps.rel@sam_data$treat)

#bet.ps.type <- betadisper(vegdist(type.rel@otu_table),type.rel@sam_data$spp)
#bet.ps.fam <- betadisper(vegdist(fam.rel@otu_table),fam.rel@sam_data$spp)
# note default vegdist is bray-curtis

anova(bet.ps.spp)
anova(bet.ps.treat)

#bet.ps.HSD.prof <- TukeyHSD(bet.ps.prof)
#plot(bet.ps.HSD.prof,las=1)

#### adonis tests 
# do species have different composition? YES
samdf=data.frame(ps.rel@sam_data)
adonis2(ps.rel@otu_table ~ spp, data=samdf, permutations=999)
adonis2(type.rel@otu_table ~ spp, data=samdf, permutations=999)
adonis2(fam.rel@otu_table ~ spp, data=samdf, permutations=999)
# do treatments have different composition? NO
adonis2(ps.rel@otu_table ~ treat, data=samdf, permutations=999)
adonis2(type.rel@otu_table ~ treat, data=samdf, permutations=999)
adonis2(fam.rel@otu_table ~ treat, data=samdf, permutations=999)

library(funfuns)
pw_ad_glycan=pairwise.adonis(ps.rel@otu_table, factors=samdf$spp, permutations=999) 
pw_ad_type=pairwise.adonis(type.rel@otu_table, factors=samdf$spp, permutations=999) 
pw_ad_fam=pairwise.adonis(fam.rel@otu_table, factors=samdf$spp, permutations=999) 
# at the individual glycan level all are different from one another
# but at higher order levels (glycan type and family), Bmin different from everyone and Dtren V Smic
# meaning C is similar to Dtren and Smic in terms of general composition, but early infection rate is very different

##### alpha diversity
# need to use raw counts for this not relative abundance
# will not work yet
alpha.ps <- estimate_richness(ps, measures = c("Shannon", "Simpson"))
alpha.ps$spp <- sample_data(ps)$spp

#### NMDS plot
ord = ordinate(ps.rel, "NMDS", "bray")
plot_ordination(ps.rel, ord, type="samples", color="spp", shape="treat")
# batch effect in the composition data
# this is because there are glycans observed in one batch but not the other
# only analyze glycans observed in both batches -- still have batch effect
# subset batch 2 only
ps_batch2 <- subset_samples(ps.rel, batch == "B" & sample_names(ps.rel) != "Cladocopium.goreaui_25C_1")
batch2_filtered <- prune_taxa(taxa_sums(ps_batch2) > 0, ps_batch2)
ord = ordinate(batch2_filtered, "NMDS", "bray")
p <- plot_ordination(batch2_filtered, ord, type="samples", color="spp", shape="treat")
p + 
  geom_point(size=2) + 
  stat_ellipse(aes(group=spp)) +
  scale_color_manual(values = c("Bmin"="#0571B0","Smic"="#4DAC26","Cgor"="gold","Dtren"="#E66101")) +
  theme_bw()
