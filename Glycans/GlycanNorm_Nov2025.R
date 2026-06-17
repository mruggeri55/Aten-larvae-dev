setwd("/Users/maria/Documents/work/Weis_lab/Acropora larvae infections")

#### Glycan normalization
# Nov 2025
# will perform PQN normalization across species because does not seem to be major differences and then batch correct
# will also perform PQN on only batch 2 for sanity checking
# will generate relative abundance for both datasets from abs counts

# read in absolute abundance and metadata files
Abs=read.csv("glycans/all_data/Abs_GlycanQuant_4Sep2025.csv")
meta=read.csv("glycans/all_data/glycan_meta.csv")

# transpose absolute abundance df and save glycan metadata
tAbs=as.data.frame(t(Abs[,!colnames(Abs) %in% c("Glycan.Types","Glycan.family","Glycan","Hex.Nac","Hex","Fuc","NeuAc","NeuGc","X")]))
colnames(tAbs)=Abs$Glycan
glycan_meta=Abs[,colnames(Abs) %in% c("Glycan.Types","Glycan.family","Glycan","Hex.Nac","Hex","Fuc","NeuAc","NeuGc")]
write.csv(glycan_meta,'glycans/all_data/glycan_info.csv')

# check if sample order matched between metadata and count table
meta$sample==rownames(tAbs)

##### make list of glycans exclusively in batch A or B 
uniq_A=names(tAbs)[colSums(tAbs[rownames(tAbs) %in% meta$sample[meta$batch=='B'],]) == 0]
uniq_B=names(tAbs)[colSums(tAbs[rownames(tAbs) %in% meta$sample[meta$batch=='A'],]) == 0]

#### define function for half min gap filling ####
impute_half_min <- function(data) {
  # Apply the transformation column-wise
  result <- apply(data, 2, function(x) {
    min_val <- min(x[x > 0], na.rm = TRUE)  # smallest non-zero
    x[x == 0] <- min_val / 2                # replace zeros
    return(x)
  })
  
  # Convert back to data frame (if you passed a data frame)
  if (is.data.frame(data)) {
    result <- as.data.frame(result)
  }
  
  return(result)
}

##### PQN normalization for batch 2 only ######
# subset batch 2 and remove non-informative glycans
batch2=tAbs[rownames(tAbs) %in% meta$sample[meta$batch=='B'],!colnames(tAbs) %in% uniq_A]
batch2=batch2[!rownames(batch2) %in% "Cladocopium.goreaui_25C_1",] # to sub without C1 which is technically its own batch
# fill in zeros with half min 
impute_batch2 <- impute_half_min(batch2)
# then PQN normalize
library(Rcpm)
batch2_pqn <- as.data.frame(pqn(impute_batch2, n = "median"))
#write.csv(batch2_pqn,"glycans/all_data/MR_norm/batch2only_PQN_4Nov2025.csv")
#write.csv(t(log2(batch2_pqn)),"glycans/all_data/MR_norm/batch2only_PQN_log2_4Nov2025.csv")
write.csv(batch2_pqn,"glycans/all_data/MR_norm/batch2only_noC1_PQN_25Nov2025.csv")
write.csv(t(log2(batch2_pqn)),"glycans/all_data/MR_norm/batch2only_noC1_PQN_log2_25Nov2025.csv")

#### PQN normalization for common glycans and batch correction ####
# gap fill and pqn
impute_common=impute_half_min(tAbs[,!colnames(tAbs) %in% c(uniq_A,uniq_B)])
common_pqn=as.data.frame(pqn(impute_common, n="median"))
write.csv(common_pqn,"glycans/all_data/MR_norm/common_33glycans_PQN_4Nov2025.csv")

# log transform and batch correction
batch=meta$batch
mod=model.matrix(~as.factor(spp)+as.factor(treat),data=meta)
library(limma)
common_corrected <- as.data.frame(removeBatchEffect(t(log2(common_pqn)), batch = batch, design = mod))
write.csv(common_corrected,"glycans/all_data/MR_norm/common_33glycans_PQN_log2_batch_4Nov2025.csv")

#### normalization and batch correction without Cs #####
# Cs are exclusive to batch B
impute_common_noC=impute_half_min(tAbs[!grepl("Cladocopium",rownames(tAbs)),!colnames(tAbs) %in% c(uniq_A,uniq_B)])
common_noC_pqn=as.data.frame(pqn(impute_common_noC, n="median"))
batch=meta$batch[!meta$spp=="Cgor"]
mod=model.matrix(~as.factor(spp)+as.factor(treat),data=meta[!meta$spp=="Cgor",])
library(limma)
common_corrected_noC <- as.data.frame(removeBatchEffect(t(log2(common_noC_pqn)), batch = batch, design = mod))
write.csv(common_corrected_noC,"glycans/all_data/MR_norm/common_33glycans_noC_PQN_log2_batch_11Nov2025.csv")


#### relative abundance #####
rel_all=sweep(tAbs, 1, rowSums(tAbs), FUN = "/") 
rel_common=sweep(tAbs[,!colnames(tAbs) %in% c(uniq_A,uniq_B)], 1, rowSums(tAbs[,!colnames(tAbs) %in% c(uniq_A,uniq_B)]), FUN = "/") 
# note these subsets don't actually change relative abundance but generating them for ease later

write.csv(rel_all,"glycans/all_data/MR_norm/all_87glycans_rel_abund_4Nov2025.csv")
write.csv(rel_common,"glycans/all_data/MR_norm/common_33glycans_rel_abund_4Nov2025.csv")
