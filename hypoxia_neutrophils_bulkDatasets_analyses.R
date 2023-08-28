
# Reproducibility script for Figure 4j,k, Suppl. Fig. 9a,b,c. 
# The normalized gene expression matrices and clinical tables for each public RNA-seq dataset used are not uploaded to the repository for size constraints, but they can be retrieve from the respective publications cited in the paper.
# if you wish to just reproduce the figues, comment-out the lines from 12 to 426.

library(singscore)
library(survival)

DataDir = "data/hypoxia_neutrophils_bulkDatasets_analyses/"
OutDir = "results/hypoxia_neutrophils_bulkDatasets_analyses/"

# ### Functions
# recalib = function(values,thMax,thMin,recalib_percent){
#    unit = recalib_percent/100
#    newv = ((values-(min(c(thMin,thMax))) )/abs(thMax-thMin))*(1/unit)
#    return(newv)
# }

# ### Hypoxia signature
# gavish_mps = paste0(DataDir,"gavish_mp.txt") # meta-programs from Gavish et al. 
# aa = read.table(file = gavish_mps, header = T, row.names = NULL, stringsAsFactors = F, sep = "\t", quote = '')
# hy = list(gavish_hypoxia = sort(as.character(aa$MP6.Hypoxia)),
#           TAN = c("CCR3","CCRL2","DDIT3","FLOT1","HIF1A","IRAK2","MAFF","MAP1LC3B2","MCOLN1","NBN","NOD2","PI3","PLAU","PPIF","TGM3","TOM1","UBR5-AS1","ZNF267"), 
#           NAN = c("AGO4","ARG1","CYP4F3","ERGIC1","FLOT2","FRAT2","LRP10","MGAM","MMP25","MSRB1","NDEL1","NFE2","PADI4","PBX2","PHOSPHO1","RASGRP4","REPS2","SULT1B1", "TSEN34","XKR8")) # https://www.sciencedirect.com/science/article/pii/S1535610822004998?via%3Dihub#app2
# hy$gavish_hypoxia = hy$gavish_hypoxia[!(hy$gavish_hypoxia=="DDIT3")] # the only one in common, to avoid spurious correlations
# hy$TAN = hy$TAN[!(hy$TAN=="DDIT3")] # the only one in common, to avoid spurious correlations
# hy$TRN = c(hy$TAN,hy$NAN)
# recalib_percent = 10
# colorz = c("dodgerblue4","firebrick")
# all_datasets = c("TCGA","Chen","Micke","Yokota","Beg","Shedden","Pintilie")
# cor_summary = data.frame(row.names = all_datasets, dataset = all_datasets, cor_tan = NA, pval_tan = NA, cor_nan = NA, pval_nan = NA, dataset_label = NA, Npatients = NA)
# cox_MultiVar_df = data.frame(row.names = all_datasets, dataset = all_datasets, HR_gavish_hypoxia = NA, HR_lower95_gavish_hypoxia = NA, HR_upper95_gavish_hypoxia = NA, pval_gavish_hypoxia = NA,
#   dataset_label = NA, Npatients = NA, Covariates = NA, stringsAsFactors = F)


# ### TCGA
# dataset = "TCGA"
# load("/mnt/ndata/daniele/lung_multiregion/rna-seq/Data/TCGA_expression_gdc_TPM/TCGA_LUAD_ge_TPM_GeneSymbols.RData")
# cn = as.character(colnames(ge))
# cn = cn[substr(cn,14,15) %in% c("01")] # only primary
# ge = ge[,cn]
# cn_short = substr(cn,1,12)
# colnames(ge) = cn_short
# logtpm = log2(ge+1)
# rankData = rankGenes(logtpm)
# scoredf = simpleScore(rankData, upSet = hy[["gavish_hypoxia"]])
# scores = data.frame( row.names = colnames(ge), Patient = colnames(ge), gavish_hypoxia = scoredf[colnames(ge),"TotalScore"], Dataset = paste0(dataset,"-LUAD (N=",nrow(scoredf),")"))
# for (fea in c( "TRN","TAN","NAN" ))
# {
#   scoredf = simpleScore(rankData, upSet = hy[[fea]])
#   scores[,fea] = scoredf[rownames(scores),"TotalScore"]
# }
# cor_summary[dataset,"Npatients"] = nrow(scores)
# cor_summary[dataset,"dataset_label"] = scores[1,"Dataset"]
# cor_summary[dataset,"cor_tan"] = cor(scores$gavish_hypoxia, scores$TAN, method = "pearson")
# cor_summary[dataset,"pval_tan"] = cor.test(scores$gavish_hypoxia, scores$TAN, method = "pearson")$p.value
# cor_summary[dataset,"cor_nan"] = cor(scores$gavish_hypoxia, scores$NAN, method = "pearson")
# cor_summary[dataset,"pval_nan"] = cor.test(scores$gavish_hypoxia, scores$NAN, method = "pearson")$p.value
# all_scores = scores
# load(paste0("/mnt/ndata/daniele/lung_multiregion/rna-seq/Processed/PatternSignatures/LtoS_SST/ClassicalSignature_29Oct_3PatientsSplitted/TCGA_gs_clin.RData"))
# this_Clin = cbind(scores,gs_clin[rownames(scores),c( "vital_status_num","Times","gender","age_at_initial_pathologic_diagnosis","stage_collapsed" )])
# cox_MultiVar_df[dataset,"Npatients"] = nrow(this_Clin)
# cox_MultiVar_df[dataset,"dataset_label"] = paste0(dataset,"-LUAD (N=",nrow(this_Clin),")")
# for (feature in c("TAN","NAN","gavish_hypoxia" ) )
# {
#   this_Clin[,feature] = recalib(values = this_Clin[,feature], thMax = 0.5, thMin = -0.5 , recalib_percent = recalib_percent)
#   # Splitting half
#   this_Clin$enrich = (this_Clin[,feature] > median(this_Clin[,feature]))
#   Surv_df = data.frame(Patient = this_Clin$Patient, vital_status = as.numeric(this_Clin$vital_status_num), Times = this_Clin$Times, Signature_value = this_Clin[,feature],
#            Signature = factor(this_Clin$enrich, levels = c(F,T) ), sex = this_Clin$gender, age = as.numeric(this_Clin$age_at_initial_pathologic_diagnosis), Stage = as.numeric(this_Clin$stage_collapsed) )
#   Surv_df$SurvObj = with(Surv_df, Surv(Times, as.numeric(as.character(vital_status))))
#   # Cox regressions
#   coxr = coxph(as.formula(paste0("SurvObj ~ Signature_value + sex + age + Stage")), data = Surv_df) #  Age + Sex +
#   a = (summary(coxr))
#   cox_MultiVar_df[dataset,paste0("HR_",feature)] = a$coefficients["Signature_value","exp(coef)"]
#   cox_MultiVar_df[dataset,paste0("HR_lower95_",feature)] = a$conf.int["Signature_value","lower .95"]
#   cox_MultiVar_df[dataset,paste0("HR_upper95_",feature)] = a$conf.int["Signature_value","upper .95"]
#   cox_MultiVar_df[dataset,paste0("pval_",feature)] = a$coefficients["Signature_value","Pr(>|z|)"]
#   cox_MultiVar_df[dataset,"Covariates"] = "sex, age, stage"
# }

# ### Chen
# dataset = "Chen"
# load(file = paste0("/mnt/ndata/daniele/lung_multiregion/Data/","Chen2020/ge_Chen2020_normalized_GeneSymbols.RData"))
# ge = ge[!is.na(rowSums(ge)),]
# logtpm = ge
# rankData = rankGenes(logtpm)
# scoredf = simpleScore(rankData, upSet = hy[["gavish_hypoxia"]])
# scores = data.frame( row.names = colnames(ge), Patient = colnames(ge), gavish_hypoxia = scoredf[colnames(ge),"TotalScore"], Dataset = paste0(dataset," et al. (N=",nrow(scoredf),")"))
# for (fea in c( "TRN","TAN","NAN" ))
# {
#   scoredf = simpleScore(rankData, upSet = hy[[fea]])
#   scores[,fea] = scoredf[rownames(scores),"TotalScore"]
# }
# cor_summary[dataset,"Npatients"] = nrow(scores)
# cor_summary[dataset,"dataset_label"] = scores[1,"Dataset"]
# cor_summary[dataset,"cor_tan"] = cor(scores$gavish_hypoxia, scores$TAN, method = "pearson")
# cor_summary[dataset,"pval_tan"] = cor.test(scores$gavish_hypoxia, scores$TAN, method = "pearson")$p.value
# cor_summary[dataset,"cor_nan"] = cor(scores$gavish_hypoxia, scores$NAN, method = "pearson")
# cor_summary[dataset,"pval_nan"] = cor.test(scores$gavish_hypoxia, scores$NAN, method = "pearson")$p.value
# all_scores = rbind(all_scores,scores)
# load(paste0("/mnt/ndata/daniele/lung_multiregion/rna-seq/Processed/PatternSignatures/LtoS_SST/ClassicalSignature_29Oct_3PatientsSplitted/Chen_gs_clin.RData"))
# this_Clin = gs_clin
# this_Clin$Times = as.numeric(this_Clin$OS.Month)
# this_Clin$vital_status_num = NA
# this_Clin[this_Clin$OS.Status=="Dead","vital_status_num"] = 1
# this_Clin[this_Clin$OS.Status=="Alive","vital_status_num"] = 0
# this_Clin$gender = this_Clin$Gender
# this_Clin$age_at_initial_pathologic_diagnosis = this_Clin$Age
# this_Clin = cbind(scores,this_Clin[rownames(scores),c( "vital_status_num","Times","gender","age_at_initial_pathologic_diagnosis","stage_collapsed" )])
# cox_MultiVar_df[dataset,"Npatients"] = nrow(this_Clin)
# cox_MultiVar_df[dataset,"dataset_label"] = scores[1,"Dataset"]
# for (feature in c("TAN","NAN","gavish_hypoxia" ) )
# {
#   this_Clin[,feature] = recalib(values = this_Clin[,feature], thMax = 0.5, thMin = -0.5 , recalib_percent = recalib_percent)
#   # Splitting half
#   this_Clin$enrich = (this_Clin[,feature] > median(this_Clin[,feature]))
#   Surv_df = data.frame(Patient = this_Clin$Patient, vital_status = as.numeric(this_Clin$vital_status_num), Times = this_Clin$Times, Signature_value = this_Clin[,feature],
#            Signature = factor(this_Clin$enrich, levels = c(F,T) ), sex = this_Clin$gender, age = as.numeric(this_Clin$age_at_initial_pathologic_diagnosis), Stage = as.numeric(this_Clin$stage_collapsed) )
#   Surv_df$SurvObj = with(Surv_df, Surv(Times, as.numeric(as.character(vital_status))))
#   # Cox regressions
#   coxr = coxph(as.formula(paste0("SurvObj ~ Signature_value + sex + age + Stage")), data = Surv_df) #  Age + Sex +
#   a = (summary(coxr))
#   cox_MultiVar_df[dataset,paste0("HR_",feature)] = a$coefficients["Signature_value","exp(coef)"]
#   cox_MultiVar_df[dataset,paste0("HR_lower95_",feature)] = a$conf.int["Signature_value","lower .95"]
#   cox_MultiVar_df[dataset,paste0("HR_upper95_",feature)] = a$conf.int["Signature_value","upper .95"]
#   cox_MultiVar_df[dataset,paste0("pval_",feature)] = a$coefficients["Signature_value","Pr(>|z|)"]
#   cox_MultiVar_df[dataset,"Covariates"] = "sex, age, stage"
# }


# # Beg
# dataset = "Beg"
# PcDir = "/mnt/ndata/daniele/lung_multiregion/Data/Beg2016/"
# load(file = paste0(PcDir,"Beg_LUAD_ge_microarray_GeneSymbols.RData"))
# logtpm = log2(ge)
# rankData = rankGenes(logtpm)
# scoredf = simpleScore(rankData, upSet = hy[["gavish_hypoxia"]])
# scores = data.frame( row.names = colnames(ge), Patient = colnames(ge), gavish_hypoxia = scoredf[colnames(ge),"TotalScore"], Dataset = paste0("Schabath et al. (N=",nrow(scoredf),")"))
# for (fea in c( "TRN","TAN","NAN" ))
# {
#   scoredf = simpleScore(rankData, upSet = hy[[fea]])
#   scores[,fea] = scoredf[rownames(scores),"TotalScore"]
# }
# cor_summary[dataset,"Npatients"] = nrow(scores)
# cor_summary[dataset,"dataset_label"] = scores[1,"Dataset"]
# cor_summary[dataset,"cor_tan"] = cor(scores$gavish_hypoxia, scores$TAN, method = "pearson")
# cor_summary[dataset,"pval_tan"] = cor.test(scores$gavish_hypoxia, scores$TAN, method = "pearson")$p.value
# cor_summary[dataset,"cor_nan"] = cor(scores$gavish_hypoxia, scores$NAN, method = "pearson")
# cor_summary[dataset,"pval_nan"] = cor.test(scores$gavish_hypoxia, scores$NAN, method = "pearson")$p.value
# all_scores = rbind(all_scores,scores)
# load(file = paste0(PcDir, "Clin_Beg.RData"))
# this_Clin = merge(scores, Clin[,c("Patient","vital_status_num", "times", "Stage", "Age","Gender","Smoking","KRAS_status","EGFR_status","STK11_status","TP53_status")], by = "Patient")
# rownames(this_Clin) = this_Clin$Patient
# cox_MultiVar_df[dataset,"Npatients"] = nrow(this_Clin)
# cox_MultiVar_df[dataset,"dataset_label"] = scores[1,"Dataset"]
# # Harmonization
# this_Clin$age_at_initial_pathologic_diagnosis = as.numeric(this_Clin$Age)
# this_Clin$gender = this_Clin$Gender
# this_Clin$Times = as.numeric(this_Clin$times)
# this_Clin$stage_collapsed = as.numeric(this_Clin$Stage)
# for (feature in c("TAN","NAN","gavish_hypoxia" ) )
# {
#   this_Clin[,feature] = recalib(values = this_Clin[,feature], thMax = 0.5, thMin = -0.5 , recalib_percent = recalib_percent)
#   # Splitting half
#   this_Clin$enrich = (this_Clin[,feature] > median(this_Clin[,feature]))
#   Surv_df = data.frame(Patient = this_Clin$Patient, vital_status = as.numeric(this_Clin$vital_status_num), Times = this_Clin$Times, Signature_value = this_Clin[,feature],
#            Signature = factor(this_Clin$enrich, levels = c(F,T) ), sex = this_Clin$gender, age = as.numeric(this_Clin$age_at_initial_pathologic_diagnosis), Stage = as.numeric(this_Clin$stage_collapsed) )
#   Surv_df$SurvObj = with(Surv_df, Surv(Times, as.numeric(as.character(vital_status))))
#   # Cox regressions
#   coxr = coxph(as.formula(paste0("SurvObj ~ Signature_value + sex + age + Stage")), data = Surv_df) #  Age + Sex +
#   a = (summary(coxr))
#   cox_MultiVar_df[dataset,paste0("HR_",feature)] = a$coefficients["Signature_value","exp(coef)"]
#   cox_MultiVar_df[dataset,paste0("HR_lower95_",feature)] = a$conf.int["Signature_value","lower .95"]
#   cox_MultiVar_df[dataset,paste0("HR_upper95_",feature)] = a$conf.int["Signature_value","upper .95"]
#   cox_MultiVar_df[dataset,paste0("pval_",feature)] = a$coefficients["Signature_value","Pr(>|z|)"]
#   cox_MultiVar_df[dataset,"Covariates"] = "sex, age, stage"
# }

# ### Micke
# dataset = "Micke"
# PcDir = "/mnt/ndata/daniele/lung_multiregion/Data/Micke2018/"
# load(file = paste0("/mnt/ndata/daniele/lung_multiregion/Data/Micke2018/", "Clin_Micke.RData")) # Clin
# load(file = paste0(PcDir,"Micke_LUAD_ge_FPKM_GeneSymbols.RData"))
# logtpm = log2(ge+1)
# rankData = rankGenes(logtpm)
# scoredf = simpleScore(rankData, upSet = hy[["gavish_hypoxia"]])
# scores = data.frame( row.names = colnames(ge), Patient = colnames(ge), gavish_hypoxia = scoredf[colnames(ge),"TotalScore"], Dataset = paste0("Mezheyeuski et al. (N=",nrow(scoredf),")"))
# for (fea in c( "TRN","TAN","NAN" ))
# {
#   scoredf = simpleScore(rankData, upSet = hy[[fea]])
#   scores[,fea] = scoredf[rownames(scores),"TotalScore"]
# }
# cor_summary[dataset,"Npatients"] = nrow(scores)
# cor_summary[dataset,"dataset_label"] = scores[1,"Dataset"]
# cor_summary[dataset,"cor_tan"] = cor(scores$gavish_hypoxia, scores$TAN, method = "pearson")
# cor_summary[dataset,"pval_tan"] = cor.test(scores$gavish_hypoxia, scores$TAN, method = "pearson")$p.value
# cor_summary[dataset,"cor_nan"] = cor(scores$gavish_hypoxia, scores$NAN, method = "pearson")
# cor_summary[dataset,"pval_nan"] = cor.test(scores$gavish_hypoxia, scores$NAN, method = "pearson")$p.value
# all_scores = rbind(all_scores,scores)
# load(file = paste0("/mnt/ndata/daniele/lung_multiregion/Data/Micke2018/", "Clin_Micke.RData")) # Clin
# this_Clin = merge(scores, Clin[,c("Patient","vital_status_num", "times", "Stage", "Age","Gender")], by = "Patient")
# rownames(this_Clin) = this_Clin$Patient
# cox_MultiVar_df[dataset,"Npatients"] = nrow(this_Clin)
# cox_MultiVar_df[dataset,"dataset_label"] = scores[1,"Dataset"]
# # Harmonization
# this_Clin$age_at_initial_pathologic_diagnosis = as.numeric(this_Clin$Age)
# this_Clin$gender = this_Clin$Gender
# this_Clin$Times = as.numeric(this_Clin$times)
# this_Clin$stage_collapsed = as.numeric(this_Clin$Stage)
# for (feature in c("TAN","NAN","gavish_hypoxia" ) )
# {
#   this_Clin[,feature] = recalib(values = this_Clin[,feature], thMax = 0.5, thMin = -0.5 , recalib_percent = recalib_percent)
#   # Splitting half
#   this_Clin$enrich = (this_Clin[,feature] > median(this_Clin[,feature]))
#   Surv_df = data.frame(Patient = this_Clin$Patient, vital_status = as.numeric(this_Clin$vital_status_num), Times = this_Clin$Times, Signature_value = this_Clin[,feature],
#            Signature = factor(this_Clin$enrich, levels = c(F,T) ), sex = this_Clin$gender, age = as.numeric(this_Clin$age_at_initial_pathologic_diagnosis), Stage = as.numeric(this_Clin$stage_collapsed) )
#   Surv_df$SurvObj = with(Surv_df, Surv(Times, as.numeric(as.character(vital_status))))
#   # Cox regressions
#   coxr = coxph(as.formula(paste0("SurvObj ~ Signature_value + sex + age + Stage")), data = Surv_df) #  Age + Sex +
#   a = (summary(coxr))
#   cox_MultiVar_df[dataset,paste0("HR_",feature)] = a$coefficients["Signature_value","exp(coef)"]
#   cox_MultiVar_df[dataset,paste0("HR_lower95_",feature)] = a$conf.int["Signature_value","lower .95"]
#   cox_MultiVar_df[dataset,paste0("HR_upper95_",feature)] = a$conf.int["Signature_value","upper .95"]
#   cox_MultiVar_df[dataset,paste0("pval_",feature)] = a$coefficients["Signature_value","Pr(>|z|)"]
#   cox_MultiVar_df[dataset,"Covariates"] = "sex, age, stage"
# }



# ### Pintilie
# dataset = "Pintilie"
# PcDir = "/mnt/ndata/daniele/lung_multiregion/Data/Pintilie2013/"
# load(file = paste0(PcDir,"Pintilie_LUAD_ge_logcounts_GeneSymbols.RData"))
# logtpm = ge
# rankData = rankGenes(logtpm)
# scoredf = simpleScore(rankData, upSet = hy[["gavish_hypoxia"]])
# scores = data.frame( row.names = colnames(ge), Patient = colnames(ge), gavish_hypoxia = scoredf[colnames(ge),"TotalScore"], Dataset = paste0("Der et al. (N=",nrow(scoredf),")"))
# for (fea in c( "TRN","TAN","NAN" ))
# {
#   scoredf = simpleScore(rankData, upSet = hy[[fea]])
#   scores[,fea] = scoredf[rownames(scores),"TotalScore"]
# }
# cor_summary[dataset,"Npatients"] = nrow(scores)
# cor_summary[dataset,"dataset_label"] = scores[1,"Dataset"]
# cor_summary[dataset,"cor_tan"] = cor(scores$gavish_hypoxia, scores$TAN, method = "pearson")
# cor_summary[dataset,"pval_tan"] = cor.test(scores$gavish_hypoxia, scores$TAN, method = "pearson")$p.value
# cor_summary[dataset,"cor_nan"] = cor(scores$gavish_hypoxia, scores$NAN, method = "pearson")
# cor_summary[dataset,"pval_nan"] = cor.test(scores$gavish_hypoxia, scores$NAN, method = "pearson")$p.value
# all_scores = rbind(all_scores,scores)
# load(file = paste0(PcDir, "Clin_Pintilie.RData"))
# this_Clin = merge(scores, Clin[,c("Patient","vital_status_num", "times", "Stage", "Age","Gender","Smoking")], by = "Patient")
# rownames(this_Clin) = this_Clin$Patient
# cox_MultiVar_df[dataset,"Npatients"] = nrow(this_Clin)
# cox_MultiVar_df[dataset,"dataset_label"] = scores[1,"Dataset"]
# # Harmonization
# this_Clin$age_at_initial_pathologic_diagnosis = as.numeric(this_Clin$Age)
# this_Clin$gender = this_Clin$Gender
# this_Clin$Times = as.numeric(this_Clin$times)
# this_Clin$stage_collapsed = as.numeric(this_Clin$Stage)
# for (feature in c("TAN","NAN","gavish_hypoxia" ) )
# {
#   this_Clin[,feature] = recalib(values = this_Clin[,feature], thMax = 0.5, thMin = -0.5 , recalib_percent = recalib_percent)
#   # Splitting half
#   this_Clin$enrich = (this_Clin[,feature] > median(this_Clin[,feature]))
#   Surv_df = data.frame(Patient = this_Clin$Patient, vital_status = as.numeric(this_Clin$vital_status_num), Times = this_Clin$Times, Signature_value = this_Clin[,feature],
#            Signature = factor(this_Clin$enrich, levels = c(F,T) ), sex = this_Clin$gender, age = as.numeric(this_Clin$age_at_initial_pathologic_diagnosis), Stage = as.numeric(this_Clin$stage_collapsed) )
#   Surv_df$SurvObj = with(Surv_df, Surv(Times, as.numeric(as.character(vital_status))))
#   # Cox regressions
#   coxr = coxph(as.formula(paste0("SurvObj ~ Signature_value + sex + age + Stage")), data = Surv_df) #  Age + Sex +
#   a = (summary(coxr))
#   cox_MultiVar_df[dataset,paste0("HR_",feature)] = a$coefficients["Signature_value","exp(coef)"]
#   cox_MultiVar_df[dataset,paste0("HR_lower95_",feature)] = a$conf.int["Signature_value","lower .95"]
#   cox_MultiVar_df[dataset,paste0("HR_upper95_",feature)] = a$conf.int["Signature_value","upper .95"]
#   cox_MultiVar_df[dataset,paste0("pval_",feature)] = a$coefficients["Signature_value","Pr(>|z|)"]
#   cox_MultiVar_df[dataset,"Covariates"] = "sex, age, stage"
# }



# # Shedden
# dataset = "Shedden"
# PcDir = "/mnt/ndata/daniele/lung_multiregion/Data/Shedden2008/"
# load(file = paste0(PcDir,"Shedden_LUAD_ge_rawcounts_GeneSymbols.RData"))
# logtpm = log2(ge)
# rankData = rankGenes(logtpm)
# scoredf = simpleScore(rankData, upSet = hy[["gavish_hypoxia"]])
# scores = data.frame( row.names = colnames(ge), Patient = colnames(ge), gavish_hypoxia = scoredf[colnames(ge),"TotalScore"], Dataset = paste0(dataset," et al. (N=",nrow(scoredf),")"))
# for (fea in c( "TRN","TAN","NAN" ))
# {
#   scoredf = simpleScore(rankData, upSet = hy[[fea]])
#   scores[,fea] = scoredf[rownames(scores),"TotalScore"]
# }
# cor_summary[dataset,"Npatients"] = nrow(scores)
# cor_summary[dataset,"dataset_label"] = scores[1,"Dataset"]
# cor_summary[dataset,"cor_tan"] = cor(scores$gavish_hypoxia, scores$TAN, method = "pearson")
# cor_summary[dataset,"pval_tan"] = cor.test(scores$gavish_hypoxia, scores$TAN, method = "pearson")$p.value
# cor_summary[dataset,"cor_nan"] = cor(scores$gavish_hypoxia, scores$NAN, method = "pearson")
# cor_summary[dataset,"pval_nan"] = cor.test(scores$gavish_hypoxia, scores$NAN, method = "pearson")$p.value
# all_scores = rbind(all_scores,scores)
# load(file = paste0(PcDir, "Clin_Shedden.RData"))
# this_Clin = merge(scores, Clin[,c("Patient","vital_status_num","times","Stage","Age","Gender", "Adjuvant_chemo", "Adjuvant_rt", "ever_smoked", "grade")], by = "Patient")
# rownames(this_Clin) = this_Clin$Patient
# cox_MultiVar_df[dataset,"Npatients"] = nrow(this_Clin)
# cox_MultiVar_df[dataset,"dataset_label"] = scores[1,"Dataset"]
# # Harmonization
# this_Clin$age_at_initial_pathologic_diagnosis = as.numeric(this_Clin$Age)
# this_Clin$gender = this_Clin$Gender
# this_Clin$Times = as.numeric(this_Clin$times)
# this_Clin$stage_collapsed = as.numeric(this_Clin$Stage)
# for (feature in c("TAN","NAN","gavish_hypoxia" ) )
# {
#   this_Clin[,feature] = recalib(values = this_Clin[,feature], thMax = 0.5, thMin = -0.5 , recalib_percent = recalib_percent)
#   # Splitting half
#   this_Clin$enrich = (this_Clin[,feature] > median(this_Clin[,feature]))
#   Surv_df = data.frame(Patient = this_Clin$Patient, vital_status = as.numeric(this_Clin$vital_status_num), Times = this_Clin$Times, Signature_value = this_Clin[,feature],
#            Signature = factor(this_Clin$enrich, levels = c(F,T) ), sex = this_Clin$gender, age = as.numeric(this_Clin$age_at_initial_pathologic_diagnosis), Stage = as.numeric(this_Clin$stage_collapsed) )
#   Surv_df$SurvObj = with(Surv_df, Surv(Times, as.numeric(as.character(vital_status))))
#   # Cox regressions
#   coxr = coxph(as.formula(paste0("SurvObj ~ Signature_value + sex + age + Stage")), data = Surv_df) #  Age + Sex +
#   a = (summary(coxr))
#   cox_MultiVar_df[dataset,paste0("HR_",feature)] = a$coefficients["Signature_value","exp(coef)"]
#   cox_MultiVar_df[dataset,paste0("HR_lower95_",feature)] = a$conf.int["Signature_value","lower .95"]
#   cox_MultiVar_df[dataset,paste0("HR_upper95_",feature)] = a$conf.int["Signature_value","upper .95"]
#   cox_MultiVar_df[dataset,paste0("pval_",feature)] = a$coefficients["Signature_value","Pr(>|z|)"]
#   cox_MultiVar_df[dataset,"Covariates"] = "sex, age, stage"
# }


# # Yokota
# dataset = "Yokota"
# PcDir = "/mnt/ndata/daniele/lung_multiregion/Data/Yokota2012/"
# load(file = paste0(PcDir,"Yokota_LUAD_ge_rawcounts_GeneSymbols.RData"))
# logtpm = log2(ge)
# rankData = rankGenes(logtpm)
# scoredf = simpleScore(rankData, upSet = hy[["gavish_hypoxia"]])
# scores = data.frame( row.names = colnames(ge), Patient = colnames(ge), gavish_hypoxia = scoredf[colnames(ge),"TotalScore"], Dataset = paste0("Okayama et al. (N=",nrow(scoredf),")"))
# for (fea in c( "TRN","TAN","NAN" ))
# {
#   scoredf = simpleScore(rankData, upSet = hy[[fea]])
#   scores[,fea] = scoredf[rownames(scores),"TotalScore"]
# }
# cor_summary[dataset,"Npatients"] = nrow(scores)
# cor_summary[dataset,"dataset_label"] = scores[1,"Dataset"]
# cor_summary[dataset,"cor_tan"] = cor(scores$gavish_hypoxia, scores$TAN, method = "pearson")
# cor_summary[dataset,"pval_tan"] = cor.test(scores$gavish_hypoxia, scores$TAN, method = "pearson")$p.value
# cor_summary[dataset,"cor_nan"] = cor(scores$gavish_hypoxia, scores$NAN, method = "pearson")
# cor_summary[dataset,"pval_nan"] = cor.test(scores$gavish_hypoxia, scores$NAN, method = "pearson")$p.value
# all_scores = rbind(all_scores,scores)
# load(file = paste0(PcDir, "Clin_Yokota.RData"))
# this_Clin = merge(scores, Clin[,c("Patient","vital_status_num", "times", "Stage", "Age","Gender")], by = "Patient")
# rownames(this_Clin) = this_Clin$Patient
# cox_MultiVar_df[dataset,"Npatients"] = nrow(this_Clin)
# cox_MultiVar_df[dataset,"dataset_label"] = scores[1,"Dataset"]
# # Harmonization
# this_Clin$age_at_initial_pathologic_diagnosis = as.numeric(this_Clin$Age)
# this_Clin$gender = this_Clin$Gender
# this_Clin$Times = as.numeric(this_Clin$times)
# this_Clin$stage_collapsed = as.numeric(this_Clin$Stage)
# for (feature in c("TAN","NAN","gavish_hypoxia" ) )
# {
#   this_Clin[,feature] = recalib(values = this_Clin[,feature], thMax = 0.5, thMin = -0.5 , recalib_percent = recalib_percent)
#   # Splitting half
#   this_Clin$enrich = (this_Clin[,feature] > median(this_Clin[,feature]))
#   Surv_df = data.frame(Patient = this_Clin$Patient, vital_status = as.numeric(this_Clin$vital_status_num), Times = this_Clin$Times, Signature_value = this_Clin[,feature],
#            Signature = factor(this_Clin$enrich, levels = c(F,T) ), sex = this_Clin$gender, age = as.numeric(this_Clin$age_at_initial_pathologic_diagnosis), Stage = as.numeric(this_Clin$stage_collapsed) )
#   Surv_df$SurvObj = with(Surv_df, Surv(Times, as.numeric(as.character(vital_status))))
#   # Cox regressions
#   coxr = coxph(as.formula(paste0("SurvObj ~ Signature_value + sex + age + Stage")), data = Surv_df) #  Age + Sex +
#   a = (summary(coxr))
#   cox_MultiVar_df[dataset,paste0("HR_",feature)] = a$coefficients["Signature_value","exp(coef)"]
#   cox_MultiVar_df[dataset,paste0("HR_lower95_",feature)] = a$conf.int["Signature_value","lower .95"]
#   cox_MultiVar_df[dataset,paste0("HR_upper95_",feature)] = a$conf.int["Signature_value","upper .95"]
#   cox_MultiVar_df[dataset,paste0("pval_",feature)] = a$coefficients["Signature_value","Pr(>|z|)"]
#   cox_MultiVar_df[dataset,"Covariates"] = "sex, age, stage"
# }



# # Ding
# dataset = "Ding"
# a=load(file = paste0("/mnt/ndata/daniele/lung_multiregion/rna-seq/Processed/PatternSignatures/LtoSsignature_onDing/MicropapillaryMax10/","Ding_logtpm_AllGenes.RData"))
# rankData = rankGenes(logtpm)
# scoredf = simpleScore(rankData, upSet = hy[["gavish_hypoxia"]])
# scores = data.frame( row.names = colnames(logtpm), Patient = colnames(logtpm), gavish_hypoxia = scoredf[colnames(logtpm),"TotalScore"], Dataset = paste0(dataset," et al. (N=",nrow(scoredf),")"))
# for (fea in c( "TRN","TAN","NAN" ))
# {
#   scoredf = simpleScore(rankData, upSet = hy[[fea]])
#   scores[,fea] = scoredf[rownames(scores),"TotalScore"]
# }
# cor_summary[dataset,"Npatients"] = nrow(scores)
# cor_summary[dataset,"dataset_label"] = scores[1,"Dataset"]
# cor_summary[dataset,"cor_tan"] = cor(scores$gavish_hypoxia, scores$TAN, method = "pearson")
# cor_summary[dataset,"pval_tan"] = cor.test(scores$gavish_hypoxia, scores$TAN, method = "pearson")$p.value
# cor_summary[dataset,"cor_nan"] = cor(scores$gavish_hypoxia, scores$NAN, method = "pearson")
# cor_summary[dataset,"pval_nan"] = cor.test(scores$gavish_hypoxia, scores$NAN, method = "pearson")$p.value
# all_scores = rbind(all_scores,scores)

# # LuMu
# dataset = "LuMu"
# load( file = "/mnt/ndata/daniele/lung_multiregion/rna-seq/Processed/ExprTables/lung_multiregion_log2TPMcombat_table_GeneSymbols.RData" )
# logtpm = ge
# rankData = rankGenes(logtpm)
# scoredf = simpleScore(rankData, upSet = hy[["gavish_hypoxia"]])
# scores = data.frame( row.names = colnames(logtpm), Patient = colnames(logtpm), gavish_hypoxia = scoredf[colnames(logtpm),"TotalScore"], Dataset = paste0("Tavernari et al. (N=",nrow(scoredf),")"))
# for (fea in c( "TRN","TAN","NAN" ))
# {
#   scoredf = simpleScore(rankData, upSet = hy[[fea]])
#   scores[,fea] = scoredf[rownames(scores),"TotalScore"]
# }
# cor_summary[dataset,"Npatients"] = nrow(scores)
# cor_summary[dataset,"dataset_label"] = scores[1,"Dataset"]
# cor_summary[dataset,"cor_tan"] = cor(scores$gavish_hypoxia, scores$TAN, method = "pearson")
# cor_summary[dataset,"pval_tan"] = cor.test(scores$gavish_hypoxia, scores$TAN, method = "pearson")$p.value
# cor_summary[dataset,"cor_nan"] = cor(scores$gavish_hypoxia, scores$NAN, method = "pearson")
# cor_summary[dataset,"pval_nan"] = cor.test(scores$gavish_hypoxia, scores$NAN, method = "pearson")$p.value
# all_scores = rbind(all_scores,scores)

# all_scores$Dataset = factor(all_scores$Dataset, levels = c( "TCGA-LUAD (N=513)","Shedden et al. (N=443)","Schabath et al. (N=398)","Okayama et al. (N=226)","Der et al. (N=181)","Chen et al. (N=169)","Mezheyeuski et al. (N=106)","Ding et al. (N=38)","Tavernari et al. (N=29)" ))
# cor_summary_melted = data.frame( dataset_label = c(cor_summary$dataset_label,cor_summary$dataset_label), cor = c(cor_summary$cor_tan,cor_summary$cor_nan), side = rep(c("TAN","NAN" ), each = nrow(cor_summary)) ,stringsAsFactors = F)
# cor_summary_melted$dataset_label = gsub("\\(.*","",cor_summary_melted$dataset_label)
# cor_summary_melted$dataset_label = paste0(substr(cor_summary_melted$dataset_label,1,nchar(cor_summary_melted$dataset_label)-1),"\n","(N=",c(cor_summary$Npatients,cor_summary$Npatients),")")
# cor_summary_melted$dataset_label = factor(cor_summary_melted$dataset_label,levels = c("TCGA-LUAD\n(N=513)","Shedden et al.\n(N=443)","Schabath et al.\n(N=398)","Okayama et al.\n(N=226)","Der et al.\n(N=181)","Chen et al.\n(N=169)","Mezheyeuski et al.\n(N=106)","Ding et al.\n(N=38)","Tavernari et al.\n(N=29)"))
# save(all_scores,file = paste0(OutDir,"AllDatasets_vs_hypoxia_scatterplots_table.RData"))
# save(cor_summary_melted,file = paste0(OutDir,"AllDatasets_TanNan_cor_vs_hypoxia_Barplot_table.RData"))
# save(cox_MultiVar_df,file=paste0(OutDir,"cox_MultiVar_df.RData") )







load(file = paste0(OutDir,"AllDatasets_TanNan_cor_vs_hypoxia_Barplot_table.RData"))
load(file = paste0(OutDir,"AllDatasets_vs_hypoxia_scatterplots_table.RData"))
load(file=paste0(OutDir,"cox_MultiVar_df.RData") )

library(ggplot2)
library(ggpubr)
library(ggplot2); library(gridExtra)

pdf( paste0(OutDir,"SupplFig9a_AllDatasets_TAN_vs_hypoxia.pdf"), width = 6.5, height = 6.5, useDingbats = F)
ggplot(all_scores, aes(gavish_hypoxia, TAN)) +
  geom_point(size=0.1) + stat_cor(p.digits = 1, r.accuracy = 0.01, label.sep = '\n', label.x.npc = "left", label.y.npc = "top", output.type = "text")+
  geom_smooth(method="lm", color = "firebrick", ) + xlab( "Hypoxia score" ) + ylab("TAN score") + theme_classic() +
  facet_wrap(~ Dataset, scales = "free")
dev.off()

pdf( paste0(OutDir,"SupplFig9b_AllDatasets_NAN_vs_hypoxia.pdf"), width = 6.5, height = 6.5, useDingbats = F)
ggplot(all_scores, aes(gavish_hypoxia, NAN)) +
  geom_point(size=0.1) + stat_cor(p.digits = 1, r.accuracy = 0.01, label.sep = '\n', label.x.npc = "left", label.y.npc = "top", output.type = "text")+
  geom_smooth(method="lm", color = "firebrick", ) + xlab( "Hypoxia score" ) + ylab("NAN score") + theme_classic() +
  facet_wrap(~ Dataset, scales = "free")
dev.off()

# Barplot
pdf( paste0(OutDir,"Fig4j_AllDatasets_TanNan_cor_vs_hypoxia_Barplot.pdf"), width = 9, height = 4, useDingbats = F)
ggplot(cor_summary_melted, aes(x=dataset_label, y=cor,fill=side)) +
  geom_bar(stat="identity", position=position_dodge(), color = 'black') + scale_fill_manual(values=c("#084887","#E69F00"))+ 
  xlab( "" ) + ylab("Pearson R with hypoxia score") + labs(fill="") + theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

library(forestplot)
clip = c(0.2,10)
load(file=paste0(OutDir,"cox_MultiVar_df.RData") ) # cox_MultiVar_df
cox_MultiVar_df = cox_MultiVar_df
   tabletext = cbind(
     c("Dataset", substr(cox_MultiVar_df$dataset_label,1,nchar(cox_MultiVar_df$dataset_label)-8), "Summary"),
     c("# patients", cox_MultiVar_df$Npatients, sum(cox_MultiVar_df$Npatients)),
     c("HR\n(TAN)", signif(cox_MultiVar_df$HR_TAN,3), signif(weighted.mean(x = cox_MultiVar_df$HR_TAN, w = cox_MultiVar_df$Npatients),3)),
     c("p-value\n(TAN)", formatC(signif(cox_MultiVar_df$pval_TAN,2)), NA),
     c("HR\n(NAN)", signif(cox_MultiVar_df$HR_NAN,3), signif(weighted.mean(x = cox_MultiVar_df$HR_NAN, w = cox_MultiVar_df$Npatients),3)),
     c("p-value\n(NAN)", formatC(signif(cox_MultiVar_df$pval_NAN,2)), NA),
     c("HR\n(hypoxia)", signif(cox_MultiVar_df$HR_gavish_hypoxia,3), signif(weighted.mean(x = cox_MultiVar_df$HR_gavish_hypoxia, w = cox_MultiVar_df$Npatients),3)),
     c("p-value\n(hypoxia)", formatC(signif(cox_MultiVar_df$pval_gavish_hypoxia,2)), NA),
     rep("      ",length(c("Dataset", rownames(cox_MultiVar_df), "Summary"))))

   pdf(paste0(OutDir,"Fig4k_SupplFig9c_cox_MultiVar_df_ForestPlot_AllDatasets_TAN_NAN_Hypoxia.pdf"),12,10, useDingbats=F,onefile=F)
   forestplot(tabletext, 
            legend = c("TAN","NAN","hypoxia"),
            align = c("l","l","l"),
            fn.ci_norm = c(fpDrawCircleCI,fpDrawCircleCI,fpDrawCircleCI),
           hrzl_lines = gpar(col="#444444"),
           colgap = unit(5,"mm"),
            boxsize = 1.5*cbind(c(NA,cox_MultiVar_df$Npatients/max(cox_MultiVar_df$Npatients),NA)/4,c(NA,cox_MultiVar_df$Npatients/max(cox_MultiVar_df$Npatients),NA)/4),
            mean = cbind(c(NA,cox_MultiVar_df$HR_TAN,NA),c(NA,cox_MultiVar_df$HR_NAN,NA),c(NA,cox_MultiVar_df$HR_gavish_hypoxia,NA)), 
            lower = cbind(c(NA,cox_MultiVar_df$HR_lower95_TAN,NA),c(NA,cox_MultiVar_df$HR_lower95_NAN,NA),c(NA,cox_MultiVar_df$HR_lower95_gavish_hypoxia,NA)), 
            upper = cbind(c(NA,cox_MultiVar_df$HR_upper95_TAN,NA),c(NA,cox_MultiVar_df$HR_upper95_NAN,NA),c(NA,cox_MultiVar_df$HR_upper95_gavish_hypoxia,NA)),
           is.summary=c(TRUE,rep(FALSE,length(cox_MultiVar_df$Npatients)),TRUE),
           xlog=T,
           clip = clip,
           xlab = "Hazard ratio",
           col=fpColors(box=c("#909CC2","#084887","#F58A07"),line=c("#909CC2","#084887","#F58A07")))
   dev.off()

