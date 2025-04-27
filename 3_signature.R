#突变特征----
library(BSgenome.Hsapiens.UCSC.hg19, quietly = TRUE)
library(openxlsx)
library(maftools)
maf <- read.maf(maf = "somatic_mutation_nonsys.maf")
sample_info<-read.xlsx("clin_samples.xlsx")

laml.tnm = trinucleotideMatrix(maf = maf, prefix = NULL, 
                               add = TRUE, 
                               ref_genome = "BSgenome.Hsapiens.UCSC.hg19")#'chr'
pdf("APOBEC富集与非富集的区别.pdf")
plotApobecDiff(tnm = laml.tnm, maf = maf, pVal = 1)
dev.off()
library(NMF)
laml.sign<-estimateSignatures(mat = laml.tnm, nTry = 10,pConstant = 0.001)
#绘制elbow曲线图，根据上述结果可视化并确定最佳signatures数量。Cophenetic相关性显著下降的值;
pdf("elbow_curve.pdf")
plotCophenetic(res = laml.sign)
dev.off()
laml.sig = extractSignatures(mat = laml.tnm, n = 2,pConstant = 0.001)

#Compate against original 30 signatures
laml.og30.cosm = compareSignatures(nmfRes = laml.sig, 
                                   sig_db = "legacy")#"SBS"


library('pheatmap')
pdf("mutation_2signature_heatmap.pdf")
pheatmap::pheatmap(mat = laml.og30.cosm$cosine_similarities,
                   cluster_rows = FALSE,
                   main = "cosine similarity against validated signatures")
dev.off()
pdf("mutation_2signature.pdf")
plotSignatures(nmfRes = laml.sig, title_size = 1.2, sig_db = "legacy")
dev.off()

signature_data<-as.data.frame(t(laml.sig$contributions))
signature_data$signature_group<-max.col(signature_data[,1:2])
signature_data$signature_group<-ifelse(signature_data$signature_group == 1,"APOBEC","Smoking")
#signature_data$signature_group<-ifelse(signature_data$signature_group == "Signature_1","APOBEC","Smoking")
signature_data$Sample_ID<-rownames(signature_data)
signature_data<-merge(signature_data,sample_info[,which(colnames(sample_info) %in% c("Sample_ID","OS_STATUS","OS_DAYS","PFS_STATUS","PFS_DAYS"))],by="Sample_ID")
signature_data$OS_STATUS<-gsub("Yes",1,signature_data$OS_STATUS)
signature_data$OS_STATUS<-gsub("No",0,signature_data$OS_STATUS)
signature_data$OS_STATUS<-as.numeric(signature_data$OS_STATUS)
signature_data$PFS_STATUS<-gsub("No",0,signature_data$PFS_STATUS)
signature_data$PFS_STATUS<-gsub("Yes",1,signature_data$PFS_STATUS)
signature_data$PFS_STATUS<-as.numeric(signature_data$PFS_STATUS)

library(survival)
library(survminer)
pdf("2signature_OS_output.pdf",onefile = F)
signature_fit <- survfit(Surv(OS_DAYS, OS_STATUS) ~ signature_group, data = signature_data)
ggsurvplot(signature_fit,
           pval = TRUE, #conf.int = TRUE,
           xlab = "survival time(days)",
           ylab = "OS",
           #palette = c("#E4F6A8","#A8CFE6"),
           #legend = c(0.8,0.8), # 指定图例位置
           legend.title = "", # 设置图例标题
           legend.labs = c("APOBEC","Smoking"),#"POLE",,"Unknown"
           risk.table = TRUE,
)
dev.off()

pdf("2signature_PFS_output.pdf",onefile = F)
signature_fit <- survfit(Surv(PFS_DAYS, PFS_STATUS) ~ signature_group, data = signature_data)
ggsurvplot(signature_fit,
           pval = TRUE, #conf.int = TRUE,
           xlab = "survival time(days)",
           ylab = "PFS",
           #palette = c("#E4F6A8","#A8CFE6"),
           legend = c(0.8,0.8), # 指定图例位置
           legend.title = "", # 设置图例标题
           legend.labs = c("APOBEC","Smoking"),
           risk.table = TRUE,
)
dev.off()

laml.sig_4 = extractSignatures(mat = laml.tnm, n = 4,pConstant = 0.001)

#Compate against original 30 signatures
laml.og30.cosm_4 = compareSignatures(nmfRes = laml.sig_4, 
                                     sig_db = "legacy")#"SBS"


pdf("mutation_4signature_heatmap.pdf")
pheatmap::pheatmap(mat = laml.og30.cosm_4$cosine_similarities,
                   cluster_rows = FALSE,
                   main = "cosine similarity against validated signatures")
dev.off()
pdf("mutation_4signature.pdf")
plotSignatures(nmfRes = laml.sig_4, title_size = 1.2, sig_db = "legacy")
dev.off()

signature_data_4<-as.data.frame(t(laml.sig_4$contributions))
signature_data_4$signature_group<-max.col(signature_data_4[,1:4])
signature_data_4$signature_group<-ifelse(signature_data_4$signature_group == 1,"APOBEC",
                                         ifelse(signature_data_4$signature_group == 2,"Smoking",
                                                ifelse(signature_data_4$signature_group == 3,"Smoking","POLE")))
signature_data_4$Sample_ID<-rownames(signature_data_4)
signature_data_4<-merge(signature_data_4,sample_info[,which(colnames(sample_info) %in% c("Sample_ID","OS_STATUS","OS_DAYS","PFS_STATUS","PFS_DAYS"))],by="Sample_ID")
signature_data_4$OS_STATUS<-gsub("Yes",1,signature_data_4$OS_STATUS)
signature_data_4$OS_STATUS<-gsub("No",0,signature_data_4$OS_STATUS)
signature_data_4$OS_STATUS<-as.numeric(signature_data_4$OS_STATUS)
signature_data_4$PFS_STATUS<-gsub("No",0,signature_data_4$PFS_STATUS)
signature_data_4$PFS_STATUS<-gsub("Yes",1,signature_data_4$PFS_STATUS)
signature_data_4$PFS_STATUS<-as.numeric(signature_data_4$PFS_STATUS)

signature_data_4$signature_group1<-ifelse(signature_data_4$signature_group =="POLE","POLE","non-POLE")
pdf("POLE_4signature_OS_output.pdf",onefile = F)
signature_fit <- survfit(Surv(OS_DAYS, OS_STATUS) ~ signature_group1, data = signature_data_4)
ggsurvplot(signature_fit,
           pval = TRUE, #conf.int = TRUE,
           xlab = "survival time(days)",
           ylab = "OS",
           #palette = c("#E4F6A8","#A8CFE6"),
           #legend = c(0.8,0.8), # 指定图例位置
           legend.title = "", # 设置图例标题
           #legend.labs = c("non-POLE","POLE"),#,"Smoking","Unknown"
           risk.table = TRUE,
)
dev.off()

pdf("POLE_4signature_PFS_output.pdf",onefile = F)
signature_fit <- survfit(Surv(PFS_DAYS, PFS_STATUS) ~ signature_group1, data = signature_data_4)
ggsurvplot(signature_fit,
           pval = TRUE, #conf.int = TRUE,
           xlab = "survival time(days)",
           ylab = "PFS",
           #palette = c("#E4F6A8","#A8CFE6"),
           legend = c(0.8,0.8), # 指定图例位置
           legend.title = "", # 设置图例标题
           legend.labs = c("non-POLE","POLE"),#"APOBEC","POLE","Smoking","Unknown"
           risk.table = TRUE,
)
dev.off()

####sig_waterfall
load("signature.RData")
var_data<-read.table("new_var_filter.txt",header = T,quote = "",sep = "\t")
signature_data<-merge(signature_data,var_data[,c("Sample_ID","RECIST", "TMB_group","new_PDL1")], by = "Sample_ID")
colnames(signature_data)[14:15]<-c("TMB", "PDL1")
signature_data<-signature_data[match(colnames(top_bar_data),signature_data$Sample_ID),]

signature_data$Smoking<-ifelse(signature_data$Smoking == "Yes","Smoking","Non-smoking")
signature_data$Drinking<-ifelse(signature_data$Drinking == "Yes","Drinking","Non_drinking")
signature_data$Age<-ifelse(signature_data$Age == 1,"<=65",">65")
signature_data$Gender<-ifelse(signature_data$Gender == 1,"Male","Female")
signature_data$TMB <- ifelse(signature_data$TMB == 0, "TMB-L", "TMB-H")
signature_data$RECIST <- ifelse(signature_data$RECIST == "待评估","Nx",signature_data$RECIST)


identical(signature_data$Sample_ID,colnames(top_bar_data))
colann=HeatmapAnnotation(df=signature_data[,c(6:11,13:15)],
                         col=list(Age=c("<=65"="#E1FFE0",">65"="#6E7B8B"),
                                  Gender=c("Female"="#FFE4E1","Male"="#6CA6CD"),
                                  Smoking=c("Non-smoking"="#FD9F73","Smoking"="#76E6B4"),
                                  Drinking=c("Non_drinking"="#CBEED8","Drinking"="#3EAF37"),
                                  HLA_group=c("heterozygous"="#F3E8D8","homozygous"="#FFB889"),
                                  HED=c("1" ="#BDDAEB","0"="#4189C1"),
                                  RECIST=c("CR"="#006DB0" ,"PR"="#9370DB","PD"="#40b8bb","SD"="#F8D568","Nx"="#F9CFF9"),
                                  TMB=c("TMB-H"="#DC9445","TMB-L"="#358DB9"),
                                  PDL1=c(">=50%"="#3DA873", "1%-49%"="#2F509E","<1%"="#8C57A2")),
                         gp = gpar(col = "white"),
                         bar=anno_barplot(as.data.frame(t(top_bar_data)),
                                          gp=gpar(fill =vc_cols,col = vc_cols),
                                          height = unit(4, "cm")))

col_dend = as.dendrogram(hclust(dist(signature_data[,2:5])))
col_dend = color_branches(col_dend, k = 4,col = c("#405EA2","#AECCE4","#B07F79","#EEBC67")) # `color_branches()` returns a dendrogram object

pdf("NEW_signature_heatmap.pdf",width = 15,height = 5)
Heatmap(t(signature_data[,2:5]),rect_gp = gpar(type="none"),
        cluster_rows = F,
        show_row_dend = F,column_dend_height = unit(3, "cm"),
        cluster_columns =  col_dend, #column_dend_gp = gpar(col = c("#405EA2","#AECCE4","#B07F79","#EEBC67")),
        show_row_names = F,show_column_names = F,
        top_annotation = colann)
dev.off()
Heatmap(laml.sig$contributions)

cut<-cutree(hclust(dist(signature_data[,2:5])),4)
# cut
# 1  2  3  4 
# 69 33 18 42 
signature_data$signature_group<-cut
signature_data$signature_group<-gsub(1,"Smoking",signature_data$signature_group)
signature_data$signature_group<-gsub(2,"Unknown",signature_data$signature_group)
signature_data$signature_group<-gsub(3,"APOBEC",signature_data$signature_group)
signature_data$signature_group<-gsub(4,"POLE",signature_data$signature_group)

chisq.test(signature_data$Age,signature_data$signature_group)#p-value =0.3341
fisher.test(signature_data$Age,signature_data$signature_group)#p-value = 0.3366

chisq.test(signature_data$Gender,signature_data$signature_group)#p-value = 0.01967可能不准;每组样本量小于5
fisher.test(signature_data$Gender,signature_data$signature_group)#p-value = 0.02672

chisq.test(signature_data$Smoking,signature_data$signature_group)#p-value =0.07399可能不准
fisher.test(signature_data$Smoking,signature_data$signature_group)#p-value = 0.0686

chisq.test(signature_data$Drinking,signature_data$signature_group)#p-value =0.8788
fisher.test(signature_data$Drinking,signature_data$signature_group)#p-value = 0.8562

chisq.test(signature_data$HLA_group,signature_data$signature_group)#p-value =0.2478可能不准
fisher.test(signature_data$HLA_group,signature_data$signature_group)#p-value = 0.2474

chisq.test(signature_data$HED,signature_data$signature_group)#p-value =0.5578
fisher.test(signature_data$HED,signature_data$signature_group)#p-value = 0.5651

chisq.test(signature_data$RECIST,signature_data$signature_group)#p-value =0.7321
fisher.test(signature_data$RECIST,signature_data$signature_group)

chisq.test(signature_data$TMB,signature_data$signature_group)#p-value =0.01057
fisher.test(signature_data$TMB,signature_data$signature_group)#p-value = 0.009664

chisq.test(signature_data$PDL1,signature_data$signature_group)#p-value = 0.3962
fisher.test(signature_data$PDL1,signature_data$signature_group)#p-value = 0.4272