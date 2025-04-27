library(openxlsx)
library(maftools)
library(survival)
library(survminer)

#data_input----
mutation <- read.table("./msk_impact_2017/data_mutations.txt", sep = "\t", quote = '', header = T, row.names = NULL, check.names = F)
patient <- read.table("./msk_impact_2017/data_clinical_patient.txt", sep = "\t", quote = '', header = T, row.names = NULL, check.names = F)
sample <- read.table("./msk_impact_2017/data_clinical_sample.txt", sep = "\t", quote = '', header = T, row.names = NULL, check.names = F)
 
clin_merge <- merge(sample, patient,by= "PATIENT_ID")
clin_merge <- clin_merge[clin_merge$OS_STATUS %in% c("0:LIVING","1:DECEASED"),]
clin_merge <- clin_merge[clin_merge$CANCER_TYPE %in% c("Non-Small Cell Lung Cancer"),]

mut_filter <- mutation[mutation$Tumor_Sample_Barcode %in% clin_merge$SAMPLE_ID,]
clin_filter <- clin_merge[clin_merge$SAMPLE_ID %in% mut_filter$Tumor_Sample_Barcode,]
# 
# write.xlsx(mut_filter, "msk_impact_mutation.xlsx")
# write.xlsx(clin_merge, "msk_impact_clinical.xlsx")

#all_mutation_gene&pathway----
mut_maf <- read.maf("msk_impact_mutation.maf")
# oncoplot(mut_maf, top = 10000, removeNonMutated = F, writeMatrix = T)

onco_mat <- read.table("msk_impact_onco_matrix.txt", sep = "\t", quote = '', check.names = F, header = T, row.names = 1)
onco_mat <- ifelse(onco_mat != "", 1,0)
allmut_onco <- as.data.frame(t(onco_mat))

#pathway
pathdb <- system.file("extdata", "oncogenic_sig_patwhays.tsv", package = "maftools")
pathdb = data.table::fread(input = pathdb)
pathdb_RTKRAS <- pathdb[pathdb$Pathway == "RTK-RAS",]
# oncoplot(mut_maf, genes = pathdb_RTKRAS$Gene, fontSize = 0.7 ,showTumorSampleBarcodes = F,writeMatrix=T,sortByMutation = T, removeNonMutated = F)

onco_rtkras <- read.table("msk_impact_rtkras_onco_matrix.txt", sep = '\t', quote = '', check.names = F, row.names = 1, header = T)
onco_rtkras[,1320:1618] <- c("")
onco_rtkras <- ifelse(onco_rtkras != "", 1,0)

rtkres_onco <- as.data.frame(t(onco_rtkras))
rtkres_onco$sum <- rowSums(rtkres_onco[,1:ncol(rtkres_onco)])
rtkres_onco$RTK_RAS <- ifelse(rtkres_onco$sum != 0, 1, 0) 
rtkres_onco <- rtkres_onco[match(clin_filter$SAMPLE_ID, rownames(rtkres_onco)),]

#RS_model_validate----
coef <- read.table("lasso_coef.txt", sep = '\t', quote = '', header = T, row.names = NULL,check.names = F)

allmut_onco <- allmut_onco[match(clin_filter$SAMPLE_ID, rownames(allmut_onco)),]
identical(clin_filter$SAMPLE_ID, rownames(allmut_onco))

clin_filter$ARID1B <- allmut_onco$ARID1B
clin_filter$SETD2 <- allmut_onco$SETD2
clin_filter$RTKRAS <- rtkres_onco$RTK_RAS
clin_filter$RISKSCORE <- clin_filter$ARID1B*coef[1,2] + clin_filter$SETD2*coef[6,2]+clin_filter$RTKRAS*coef[7,2]
clin_filter$RS <- ifelse(clin_filter$RISKSCORE > median(clin_filter$RISKSCORE),1,0)

clin_filter$OS_STATUS <- ifelse(clin_filter$OS_STATUS == "0:LIVING",0,1)
clin_filter$OS_MONTHS <- as.numeric(clin_filter$OS_MONTHS)

###
luad_filter <- clin_filter[clin_filter$ONCOTREE_CODE == "LUAD",]
lusc_filter <- clin_filter[clin_filter$ONCOTREE_CODE == "LUSC",]
fit2 <- survfit(Surv(OS_MONTHS, OS_STATUS) ~ RS, data = luad_filter)
p <- ggsurvplot(fit2, data = luad_filter,font.title=c(20, "bold"),conf.int=F, pval=TRUE, legend.title ="",
                ylab=c("Overall Survival"),  surv.median.line = "hv")
res_cox<-coxph(Surv(OS_MONTHS, OS_STATUS) ~ RS, data = luad_filter)
p$plot = p$plot + 
  ggplot2::annotate("text",x = 20, y = 0.15,
                    label = paste("HR :",round(summary(res_cox)$conf.int[1],2))) + 
  ggplot2::annotate("text",x = 20, y = 0.10,
                    label = paste("(","95%CI:",round(summary(res_cox)$conf.int[3],2),"-",round(summary(res_cox)$conf.int[4],2),")",sep = ""))
pdf("msk_impact_luad_survival.pdf", width = 5, height = 5, onefile = F)
p
dev.off()

fit2 <- survfit(Surv(OS_MONTHS, OS_STATUS) ~ RS, data = lusc_filter)
p <- ggsurvplot(fit2, data = lusc_filter,font.title=c(20, "bold"),conf.int=F, pval=TRUE, legend.title ="",
                ylab=c("Overall Survival"),  surv.median.line = "hv")
res_cox<-coxph(Surv(OS_MONTHS, OS_STATUS) ~ RS, data = lusc_filter)
p$plot = p$plot + 
  ggplot2::annotate("text",x = 20, y = 0.15,
                    label = paste("HR :",round(summary(res_cox)$conf.int[1],2))) + 
  ggplot2::annotate("text",x = 20, y = 0.10,
                    label = paste("(","95%CI:",round(summary(res_cox)$conf.int[3],2),"-",round(summary(res_cox)$conf.int[4],2),")",sep = ""))
pdf("msk_impact_lusc_survival.pdf", width = 5, height = 5, onefile = F)
p
dev.off()

Biopsy_filter <- luad_filter[luad_filter$SPECIMEN_TYPE == "Biopsy",]
Cytology_filter <- luad_filter[luad_filter$SPECIMEN_TYPE == "Cytology",]
Resection_filter <- luad_filter[luad_filter$SPECIMEN_TYPE == "Resection",]

fit2 <- survfit(Surv(OS_MONTHS, OS_STATUS) ~ RS, data = Biopsy_filter)
p <- ggsurvplot(fit2, data = Biopsy_filter,font.title=c(20, "bold"),conf.int=F, pval=TRUE, legend.title ="",
                ylab=c("Overall Survival"),  surv.median.line = "hv")
res_cox<-coxph(Surv(OS_MONTHS, OS_STATUS) ~ RS, data = Biopsy_filter)
p$plot = p$plot + 
  ggplot2::annotate("text",x = 20, y = 0.15,
                    label = paste("HR :",round(summary(res_cox)$conf.int[1],2))) + 
  ggplot2::annotate("text",x = 20, y = 0.10,
                    label = paste("(","95%CI:",round(summary(res_cox)$conf.int[3],2),"-",round(summary(res_cox)$conf.int[4],2),")",sep = ""))
pdf("msk_impact_Biopsy_survival.pdf", width = 5, height = 5, onefile = F)
p
dev.off()

fit2 <- survfit(Surv(OS_MONTHS, OS_STATUS) ~ RS, data = Cytology_filter)
p <- ggsurvplot(fit2, data = Cytology_filter,font.title=c(20, "bold"),conf.int=F, pval=TRUE, legend.title ="",
                ylab=c("Overall Survival"),  surv.median.line = "hv")
res_cox<-coxph(Surv(OS_MONTHS, OS_STATUS) ~ RS, data = Cytology_filter)
p$plot = p$plot + 
  ggplot2::annotate("text",x = 20, y = 0.15,
                    label = paste("HR :",round(summary(res_cox)$conf.int[1],2))) + 
  ggplot2::annotate("text",x = 20, y = 0.10,
                    label = paste("(","95%CI:",round(summary(res_cox)$conf.int[3],2),"-",round(summary(res_cox)$conf.int[4],2),")",sep = ""))
pdf("msk_impact_Cytology_survival.pdf", width = 5, height = 5, onefile = F)
p
dev.off()

fit2 <- survfit(Surv(OS_MONTHS, OS_STATUS) ~ RS, data = Resection_filter)
p <- ggsurvplot(fit2, data = Resection_filter,font.title=c(20, "bold"),conf.int=F, pval=TRUE, legend.title ="",
                ylab=c("Overall Survival"),  surv.median.line = "hv")
res_cox<-coxph(Surv(OS_MONTHS, OS_STATUS) ~ RS, data = Resection_filter)
p$plot = p$plot + 
  ggplot2::annotate("text",x = 20, y = 0.15,
                    label = paste("HR :",round(summary(res_cox)$conf.int[1],2))) + 
  ggplot2::annotate("text",x = 20, y = 0.10,
                    label = paste("(","95%CI:",round(summary(res_cox)$conf.int[3],2),"-",round(summary(res_cox)$conf.int[4],2),")",sep = ""))
pdf("msk_impact_Resection_survival.pdf", width = 5, height = 5, onefile = F)
p
dev.off()

####
new_filter <- clin_filter[clin_filter$ONCOTREE_CODE == "LUAD" &
                            clin_filter$SAMPLE_COLLECTION_SOURCE == "In-House" &
                            clin_filter$TUMOR_PURITY > 10,]
new_filter <- new_filter[!is.na(new_filter$PATIENT_ID),]
fit2 <- survfit(Surv(OS_MONTHS, OS_STATUS) ~ RS, data = new_filter)
p <- ggsurvplot(fit2, data = new_filter,font.title=c(20, "bold"),conf.int=F, pval=TRUE, 
                legend.title ="",palette = c("#000080","#B22222"),legend.labs = c("LRS","HRS"),
                ylab=c("Overall Survival"),  surv.median.line = "hv")
res_cox<-coxph(Surv(OS_MONTHS, OS_STATUS) ~ RS, data = new_filter)
p$plot = p$plot + 
  ggplot2::annotate("text",x = 20, y = 0.15,
                    label = paste("HR :",round(summary(res_cox)$conf.int[1],2))) + 
  ggplot2::annotate("text",x = 20, y = 0.10,
                    label = paste("(","95%CI:",round(summary(res_cox)$conf.int[3],2),"-",round(summary(res_cox)$conf.int[4],2),")",sep = ""))
pdf("new_msk_impact_luad_inhouse_tumor_survival.pdf", width = 5, height = 5, onefile = F)
p
dev.off()

Biopsy_filter <- new_filter[new_filter$SPECIMEN_TYPE == "Biopsy",]
Cytology_filter <- new_filter[new_filter$SPECIMEN_TYPE == "Cytology",]
Resection_filter <- new_filter[new_filter$SPECIMEN_TYPE == "Resection",]

fit2 <- survfit(Surv(OS_MONTHS, OS_STATUS) ~ RS, data = Biopsy_filter)
p <- ggsurvplot(fit2, data = Biopsy_filter,font.title=c(20, "bold"),conf.int=F, pval=TRUE, legend.title ="", 
                palette = c("#000080","#B22222"),legend.labs = c("LRS","HRS"),
                ylab=c("Overall Survival"),  surv.median.line = "hv")
res_cox<-coxph(Surv(OS_MONTHS, OS_STATUS) ~ RS, data = Biopsy_filter)
p$plot = p$plot + 
  ggplot2::annotate("text",x = 20, y = 0.15,
                    label = paste("HR :",round(summary(res_cox)$conf.int[1],2))) + 
  ggplot2::annotate("text",x = 20, y = 0.10,
                    label = paste("(","95%CI:",round(summary(res_cox)$conf.int[3],2),"-",round(summary(res_cox)$conf.int[4],2),")",sep = ""))
pdf("new_msk_impact_Biopsy_survival.pdf", width = 5, height = 5, onefile = F)
p
dev.off()

fit2 <- survfit(Surv(OS_MONTHS, OS_STATUS) ~ RS, data = Cytology_filter)
p <- ggsurvplot(fit2, data = Cytology_filter,font.title=c(20, "bold"),conf.int=F, pval=TRUE, legend.title ="", 
                palette = c("#000080","#B22222"),legend.labs = c("LRS","HRS"),
                ylab=c("Overall Survival"),  surv.median.line = "hv")
res_cox<-coxph(Surv(OS_MONTHS, OS_STATUS) ~ RS, data = Cytology_filter)
p$plot = p$plot + 
  ggplot2::annotate("text",x = 20, y = 0.15,
                    label = paste("HR :",round(summary(res_cox)$conf.int[1],2))) + 
  ggplot2::annotate("text",x = 20, y = 0.10,
                    label = paste("(","95%CI:",round(summary(res_cox)$conf.int[3],2),"-",round(summary(res_cox)$conf.int[4],2),")",sep = ""))
pdf("new_msk_impact_Cytology_survival.pdf", width = 5, height = 5, onefile = F)
p
dev.off()

fit2 <- survfit(Surv(OS_MONTHS, OS_STATUS) ~ RS, data = Resection_filter)
p <- ggsurvplot(fit2, data = Resection_filter,font.title=c(20, "bold"),conf.int=F, pval=TRUE, legend.title ="", 
                palette = c("#000080","#B22222"),legend.labs = c("LRS","HRS"),
                ylab=c("Overall Survival"),  surv.median.line = "hv")
res_cox<-coxph(Surv(OS_MONTHS, OS_STATUS) ~ RS, data = Resection_filter)
p$plot = p$plot + 
  ggplot2::annotate("text",x = 20, y = 0.15,
                    label = paste("HR :",round(summary(res_cox)$conf.int[1],2))) + 
  ggplot2::annotate("text",x = 20, y = 0.10,
                    label = paste("(","95%CI:",round(summary(res_cox)$conf.int[3],2),"-",round(summary(res_cox)$conf.int[4],2),")",sep = ""))
pdf("new_msk_impact_Resection_survival.pdf", width = 5, height = 5, onefile = F)
p
dev.off()