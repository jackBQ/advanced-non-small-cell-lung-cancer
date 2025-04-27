library(maftools)
library(stringr)
library(survival)
library(survminer)
library(pROC)
library(survivalROC)
library(openxlsx)

#####SEGMENT-Gistic
# a <-  list.files("./wes_segment")
# cnv_list = list()
# for (i in a){
#   # print(i)
#   cnv <- read.table(i, sep = '\t', quote = '', row.names = NULL, header = T, check.names = F,fill = T, fileEncoding = "utf-8")
#   cnv_filter <- cnv[,c("chromosome","start","end","probes","log2")]
#   cnv_filter$Sample <- str_sub(i,1,str_locate(i,"\\_T"),-1)
#   cnv_list[[i]] = as.data.frame(t(cnv_filter))
# }
# 
# cnv_sample <- as.data.frame(cnv_list)
# cnv_sample <- as.data.frame(t(cnv_sample))
# colnames(cnv_sample) <- c("Chromosome","Start","End", "Num_probes","Segment_Mean", "Sample")
# cnv_sample$Chromosome <- str_sub(cnv_sample$Chromosome, 4,-1)
# write.table(cnv_sample[,c(6,1,2,3,4,5)], "cnv_segment.txt", sep = "\t", quote = F, row.names = F)

# gistic2 -b ./wes_segment -seg ./wes_segment/cnv_segment.txt -refgene ./refgenefiles/hg19.mat -conf 0.90

#####gistic_result_summary----
clin <- read.xlsx("clin_final.xlsx", sheet = 1)
all_lesions <- read.table("all_lesions.conf_95.txt", sep = "\t", quote = '', check.names = F, header = T, row.names = NULL)
all_lesions <- all_lesions[1:119,]
rownames(all_lesions) <- paste0(str_sub(all_lesions$`Unique Name`,1,3),str_sub(all_lesions$`Unique Name`,-3,-1),"_",str_sub(all_lesions$Descriptor))

lesions <- as.data.frame(t(all_lesions[,10:94]))
lesions$Sample_ID <- rownames(lesions)
clin_lesions <- merge(lesions, clin[,c("Sample_ID","OS_STATUS", "OS_DAYS","PFS_STATUS","PFS_DAYS")], by = "Sample_ID")

#OS_Univarite&Multivarte_cox----
source("Cox.R")

clin_lesions <- merge(lesions, clin[,c("Sample_ID","OS_STATUS", "OS_DAYS","PFS_STATUS","PFS_DAYS")], by = "Sample_ID")
clin_lesions$OS_STATUS <- ifelse(clin_lesions$OS_STATUS == "Yes", 1,0)
clin_lesions$PFS_STATUS <- ifelse(clin_lesions$PFS_STATUS == "Yes",1,0)
clin_lesions$OS_DAYS <- as.numeric(clin_lesions$OS_DAYS)
clin_lesions$PFS_DAYS <- as.numeric(clin_lesions$PFS_DAYS)

for (i in 2:(ncol(clin_lesions)-4)){
  clin_lesions[,i]<-ifelse(clin_lesions[,i]>0,1,0)}

cox_clin<-clin_lesions[,c(1,121,122)]  ###OS_stats os_month
dc<-clin_lesions[,c(2:120)]

rownames(dc)<-clin_lesions[,1]
dc<-t(dc)

colnames(cox_clin)<-c("sample","Events","Survival")
cox_uni<-as.data.frame(getUniOrMultiCOXAnalysis(dc,cox_clin,method="Univariate"))

sig_uni <- cox_uni[cox_uni$`Pr(>|Z|)[pval]` < 0.05,]
dc_mul <- clin_lesions[,colnames(clin_lesions) %in% sig_uni$geneName]
rownames(dc_mul)<-clin_lesions[,1]
dc_mul<-t(dc_mul)
cox_mul<-getUniOrMultiCOXAnalysis(dc_mul,cox_clin,method="Multivariate")

# write.table(cox_uni,"os_cox_univariate.txt",sep="\t",quote=F,row.names = F)
# write.table(cox_mul,"os_cox_Multivariate.txt",sep="\t",quote=F,row.names = F)

#PFS_Univarite&Multivarte_cox----
cox_clin<-clin_lesions[,c(1,123,124)]  ###PFS_stats PFS_month
dc<-clin_lesions[,c(2:120)]

rownames(dc)<-clin_lesions[,1]
dc<-t(dc)

colnames(cox_clin)<-c("sample","Events","Survival")
cox_uni<-as.data.frame(getUniOrMultiCOXAnalysis(dc,cox_clin,method="Univariate"))

sig_uni <- cox_uni[cox_uni$`Pr(>|Z|)[pval]` < 0.05,]
dc_mul <- clin_lesions[,colnames(clin_lesions) %in% sig_uni$geneName]
rownames(dc_mul)<-clin_lesions[,1]
dc_mul<-t(dc_mul)
cox_mul<-getUniOrMultiCOXAnalysis(dc_mul,cox_clin,method="Multivariate")

# write.table(cox_uni,"pfs_cox_univariate.txt",sep="\t",quote=F,row.names = F)
# write.table(cox_mul,"pfs_cox_Multivariate.txt",sep="\t",quote=F,row.names = F)

#CNA-PeakNum_survival----
clin_lesions$cnv_num <- rowSums(clin_lesions[,c(2:120)])
clin_lesions$cnv_burden <- ifelse(clin_lesions$cnv_num >= median(clin_lesions$cnv_num), "CNA-peakNum High","CNA-peakNum Low")
clin_lesions$cnv_group <- ifelse(clin_lesions$cnv_burden == "CNA-peakNum High", 1,0)
pdf("new_os_cnv_burden_survival_median.pdf", width = 5, height = 5, onefile = F)
fit2 <- survfit(Surv(OS_DAYS, OS_STATUS) ~ cnv_burden, data = clin_lesions)
p <- ggsurvplot(fit2, data=clin_lesions,font.title=c(20, "bold"),conf.int=F, pval=TRUE, legend.title ="",ylab=c("Overall Survival"),surv.median.line = "hv",
                legend.labs = c("CNA-peakNum High","CNA-peakNum Low"),palette = c("#CD5C5C","#2E8B57"))
res_cox<-coxph(Surv(OS_DAYS, OS_STATUS) ~ cnv_group, data = clin_lesions)
p$plot = p$plot + 
  ggplot2::annotate("text",x = 80, y = 0.10,
                    label = paste("Hazard ratio = ",round(summary(res_cox)$conf.int[1],2))) + 
  ggplot2::annotate("text",x = 80, y = 0.03,
                    label = paste("95%CI:",round(summary(res_cox)$conf.int[3],2)," to ",round(summary(res_cox)$conf.int[4],2),sep = ""))
dev.off()

pdf("new_pfs_cnv_burden_survival_median.pdf", width = 5, height = 5, onefile = F)
fit2 <- survfit(Surv(PFS_DAYS, PFS_STATUS) ~ cnv_burden, data = clin_lesions)
p <- ggsurvplot(fit2, data=clin_lesions,font.title=c(20, "bold"),conf.int=F, pval=TRUE, legend.title ="",ylab=c("Progression Free Survival"),surv.median.line = "hv",
                legend.labs = c("CNA-peakNum High","CNA-peakNum Low"),palette = c("#CD5C5C","#2E8B57"))

res_cox<-coxph(Surv(PFS_DAYS, PFS_STATUS) ~ cnv_group, data = clin_lesions)
p$plot = p$plot + 
  ggplot2::annotate("text",x = 80, y = 0.10,
                    label = paste("Hazard ratio = ",round(summary(res_cox)$conf.int[1],2))) + 
  ggplot2::annotate("text",x = 80, y = 0.03,
                    label = paste("95%CI:",round(summary(res_cox)$conf.int[3],2)," to ",round(summary(res_cox)$conf.int[4],2),sep = ""))
dev.off()

data <- clin_lesions[,c(1,125,126)]
colnames(data)[3] <- c("cnv_burden")
write.xlsx(data,"new_cnv_burden_85samples.xlsx", rowNames = F)
