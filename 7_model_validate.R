library(openxlsx)
library(timeROC)
library(survival)
library(survminer)
library(RColorBrewer)

#data_input----
lasso_coef_pre <- read.table("var_lasso_group.txt", sep = '\t', quote = '', check.names = F, header = T, row.names = NULL)

wsi <- read.table("wsi_filter_var.txt", row.names = NULL, sep = '\t', quote = '', check.names = F, header = T)
stk11_comut <- read.table("sub2_mutation_data_phen_output.txt",  sep = '\t', quote = '', check.names = F, header = T, row.names = NULL)
cnv_burden <- read.xlsx("new_cnv_burden.xlsx", sheet = 1)

lasso_pre_filter <- lasso_coef_pre[,c("Sample_ID", "TMB_group", "new_PDL1", "Mean_HE", "Score", "RS", "OS_STATUS", "OS_DAYS", "PFS_STATUS", "PFS_DAYS")]
colnames(lasso_pre_filter)[5:10] <- c("pre_Score", "pre_RS", "pre_OS_STATUS", "pre_OS_DAYS", "pre_PFS_STATUS", "pre_PFS_DAYS")

wsi_filter <- wsi[,1:6]
stk11_comut <- stk11_comut[,c("Sample_ID","mutation_group")]
cnv_burden <- cnv_burden[match(stk11_comut$Sample_ID, cnv_burden$Sample_ID),]
cnv_burden$Sample_ID <- stk11_comut$Sample_ID
cnv_burden <- cnv_burden[,c("Sample_ID", "cnv_burden")]
cnv_burden$cnv_burden <- ifelse(cnv_burden$cnv_burden == "CNA-peakNum High", 0,1)

ALL_pre=list(stk11_comut, cnv_burden, wsi_filter, lasso_pre_filter)

multimerge<-function(dat=list(),...){
  if(length(dat)<2)return(as.data.frame(dat))
  mergedat<-dat[[1]]
  dat[[1]]<-NULL
  for(i in dat){
    mergedat<-merge(all=TRUE,mergedat,i,...)
  }
  return(mergedat)
}

pre_merge_mat <- multimerge(ALL_pre)
pre_merge_mat$mutation_group <- ifelse(pre_merge_mat$mutation_group == "SM", 2,
                                       ifelse(pre_merge_mat$mutation_group == "WT", 1,0))
pre_merge_mat$new_PDL1 <- ifelse(pre_merge_mat$new_PDL1 == "<1%", 2, 
                                 ifelse(pre_merge_mat$new_PDL1 == "1%-49%", 1,0))
pre_merge_mat$Neoplastic.cells <- ifelse(pre_merge_mat$Neoplastic.cells > median(pre_merge_mat$Neoplastic.cells), 1,0)
pre_merge_mat$Inflammatory <- ifelse(pre_merge_mat$Inflammatory > median(pre_merge_mat$Inflammatory), 0,1)
pre_merge_mat$Soft.tissue.cells <- ifelse(pre_merge_mat$Soft.tissue.cells > median(pre_merge_mat$Soft.tissue.cells), 0,1)
pre_merge_mat$Dead.Cells <- ifelse(pre_merge_mat$Dead.Cells > median(pre_merge_mat$Dead.Cells), 0,1)
pre_merge_mat$Epithelial_value <- pre_merge_mat$Epithelial
pre_merge_mat$Epithelial <- ifelse(pre_merge_mat$Epithelial > median(pre_merge_mat$Epithelial), 0,1)
pre_merge_mat$TMB_group <- ifelse(pre_merge_mat$TMB_group == 0, 1,0)

#roc_merge----
pre_merge_mat$RS <- ifelse(pre_merge_mat$pre_RS == "HRS", 1,0)
for (m in 2:11){
  feature <- colnames(pre_merge_mat)[m]
  pre_merge_mat[,paste0(feature,"_RS")] <- pre_merge_mat[,m]+pre_merge_mat$pre_Score
}
pre_roc_list <- list()
pre_roc_ci <- list()
for (i in c(2:12,20:29)){
  ROC <- timeROC(T = pre_merge_mat$pre_PFS_DAYS, # 将结局时间转换为月
                 delta = as.numeric(pre_merge_mat$pre_PFS_STATUS), #生存结局
                 marker = pre_merge_mat[,i], #预测变量 
                 cause = 1, #阳性结局赋值
                 weighting = "marginal", #权重计算方法，marginal是默认值，采用km计算删失分布
                 times = c(30*6, 30*8, 30*12, 30*16), #预测1、3、5、10年时间
                 iid = T) #只有marginal可以计算置信区间 
  feature <- colnames(pre_merge_mat)[i]
  pre_roc_list[[feature]] = ROC$AUC
  pre_roc_ci[[feature]] = confint(ROC, level = 0.95)$CI_AUC
}

pre_auc_list <- as.data.frame(t(as.data.frame(pre_roc_list)))

pre_auc_list[,"Max_auc"] <- apply(pre_auc_list,1,max)
rownames(pre_auc_list)[1:11] <- c("STK11_coMut", "CNV_burden", "Neoplastic_Cells", "Inflammatory","Soft_tissue_cells","Dead_Cells", "Epithelial", "TMB", "PDL1", "HED","RS")
pre_auc_list$Feature <- rownames(pre_auc_list)
# write.xlsx(pre_auc_list, "pre_auc_list_162.xlsx", rowNames = F)

###
names(pre_roc_ci)[1:11] <-  c("STK11_coMut", "CNV_burden", "Neoplastic_Cells", "Inflammatory","Soft_tissue_cells","Dead_Cells", "Epithelial", "TMB", "PDL1", "HED","RS")
pre_roc_ci <- as.data.frame(t(as.data.frame(pre_roc_ci)))
pre_roc_ci <- pre_roc_ci/100

for (i in rownames(pre_auc_list)){
  pre_auc_list[i,"CI.5_180"] <- pre_roc_ci[paste0(i,".2.5."),1]
  pre_auc_list[i,"CI.95_180"] <- pre_roc_ci[paste0(i,".97.5."),1]
  
  pre_auc_list[i,"CI.5_240"] <- pre_roc_ci[paste0(i,".2.5."),2]
  pre_auc_list[i,"CI.95_240"] <- pre_roc_ci[paste0(i,".97.5."),2] 
  
  pre_auc_list[i,"CI.5_360"] <- pre_roc_ci[paste0(i,".2.5."),3]
  pre_auc_list[i,"CI.95_360"] <- pre_roc_ci[paste0(i,".97.5."),3] 
  
  pre_auc_list[i,"CI.5_480"] <- pre_roc_ci[paste0(i,".2.5."),4]
  pre_auc_list[i,"CI.95_480"] <- pre_roc_ci[paste0(i,".97.5."),4] 
  
  pre_auc_list[i,"CI.5_max"] <- ifelse(pre_auc_list[i,"Max_auc"] == pre_auc_list[i,"t=180"],pre_roc_ci[paste0(i,".2.5."),1],
                                       ifelse(pre_auc_list[i,"Max_auc"] == pre_auc_list[i,"t=240"],pre_roc_ci[paste0(i,".2.5."),2],
                                              ifelse(pre_auc_list[i,"Max_auc"] == pre_auc_list[i,"t=360"],pre_roc_ci[paste0(i,".2.5."),3],pre_roc_ci[paste0(i,".2.5."),4])))
  pre_auc_list[i,"CI.95_max"] <- ifelse(pre_auc_list[i,"Max_auc"] == pre_auc_list[i,"t=180"],pre_roc_ci[paste0(i,".97.5."),1],
                                        ifelse(pre_auc_list[i,"Max_auc"] == pre_auc_list[i,"t=240"],pre_roc_ci[paste0(i,".97.5."),2],
                                               ifelse(pre_auc_list[i,"Max_auc"] == pre_auc_list[i,"t=360"],pre_roc_ci[paste0(i,".97.5."),3],pre_roc_ci[paste0(i,".97.5."),4])))
}

pre_auc_list_single <- pre_auc_list[1:11,]
pre_auc_list_single$Feature <- factor(pre_auc_list_single$Feature, levels = c("Neoplastic_Cells", "Inflammatory","Soft_tissue_cells","Dead_Cells", "Epithelial", "STK11_coMut", "CNV_burden", "TMB", "PDL1", "HED","RS"))
pre_auc_list_mul <- pre_auc_list[11:nrow(pre_auc_list),]
pre_auc_list_mul$Feature <- c("RS", "RS_STK11_coMut", "RS_CNV_burden", "RS_Neoplastic_Cells", "RS_Inflammatory","RS_Soft_tissue_cells","RS_Dead_Cells", "RS_Epithelial", "RS_TMB", "RS_PDL1", "RS_HED")
pre_auc_list_mul$Feature <- factor(pre_auc_list_mul$Feature,levels = c("RS", "RS_STK11_coMut", "RS_CNV_burden", "RS_Neoplastic_Cells", "RS_Inflammatory","RS_Soft_tissue_cells","RS_Dead_Cells", "RS_Epithelial", "RS_TMB", "RS_PDL1", "RS_HED"))

pdf("480_muti_feature_auc_barplot_single.pdf", width = 6, height = 5, onefile = F)
ggplot(pre_auc_list_single, aes(Feature, `t=480`, fill = Feature,label = round(`t=480`,3)))+
  geom_bar(stat="identity",width=0.5)+
  geom_errorbar(data=pre_auc_list_single,mapping=aes(x = Feature,ymin = CI.5_480, ymax = CI.95_480), width = 0.2, color = 'black', size=0.5)+
  theme_bw()+xlab(NULL)+ylab("AUC")+
  scale_fill_manual(values=c("#72484F", "#1E1E4C", "#434532", "#3B858D", "#7F624D", "#406C9D", "#6A5096", "#163224", "#5B7F5A",
                             "#36233A", "#3A4490"))+
  geom_text(size = 4, position = position_stack(vjust = 1.05))  +
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1))+
  stat_summary(fun.data = 'mean_sd', geom = "errorbar", colour = "black",
               width = 0.25,position = position_dodge( .5))+
  guides(fill = FALSE)
dev.off()

pdf("480_muti_feature_auc_barplot_cor.pdf", width = 6, height = 5, onefile = F)
ggplot(pre_auc_list_mul, aes(Feature, `t=480`, fill = Feature,label = round(`t=480`,3)))+
  geom_bar(stat="identity",width=0.5)+
  geom_errorbar(data=pre_auc_list_mul,mapping=aes(x = Feature,ymin = CI.5_480, ymax = CI.95_480), width = 0.2, color = 'black', size=0.5)+
  theme_bw()+xlab(NULL)+ylab("AUC")+
  scale_fill_manual(values=c("#3A4490", "#996273", "#41839D", "#72674D", "#3D7553", "#7A5CAB", "#536AD0", "#3A4E6B", "#85B17E", "#66B3FF","#B378B5"))+
  geom_text(size = 4, position = position_stack(vjust = 1.05))  +
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1))+
  guides(fill = FALSE)
dev.off()

#AUC_compare----
ROC_RS <- timeROC(T = pre_merge_mat$pre_PFS_DAYS, # 将结局时间转换为月
                  delta = as.numeric(pre_merge_mat$pre_PFS_STATUS), #生存结局
                  marker = pre_merge_mat[,"pre_Score"], #预测变量 
                  cause = 1, #阳性结局赋值
                  weighting = "marginal", #权重计算方法，marginal是默认值，采用km计算删失分布
                  times = c(30*6, 30*8, 30*12, 30*16), #预测1、3、5、10年时间
                  iid = T) #只有marginal可以计算置信区间 

ROC_Epithelial_RS <- timeROC(T = pre_merge_mat$pre_PFS_DAYS, # 将结局时间转换为月
                             delta = as.numeric(pre_merge_mat$pre_PFS_STATUS), #生存结局
                             marker = pre_merge_mat[,"Epithelial_RS"], #预测变量 
                             cause = 1, #阳性结局赋值
                             weighting = "marginal", #权重计算方法，marginal是默认值，采用km计算删失分布
                             times = c(30*6, 30*8, 30*12, 30*16), #预测1、3、5、10年时间
                             iid = T) #只有marginal可以计算置信区间 

ROC_Epithelial <- timeROC(T = pre_merge_mat$pre_PFS_DAYS, # 将结局时间转换为月
                          delta = as.numeric(pre_merge_mat$pre_PFS_STATUS), #生存结局
                          marker = pre_merge_mat[,"Epithelial"], #预测变量 
                          cause = 1, #阳性结局赋值
                          weighting = "marginal", #权重计算方法，marginal是默认值，采用km计算删失分布
                          times = c(30*6, 30*8, 30*12, 30*16), #预测1、3、5、10年时间
                          iid = T) #只有marginal可以计算置信区间 
timeROC::compare(ROC_RS,ROC_Epithelial_RS,adjusted=TRUE)
timeROC::compare(ROC_Epithelial,ROC_Epithelial_RS,adjusted=TRUE)

#uni_forestplot----
pre_merge_mat$pre_RS <- ifelse(pre_merge_mat$pre_RS == "HRS",1,0)
source("Cox.R")

cox_clin<-pre_merge_mat[,c("Sample_ID","pre_PFS_STATUS","pre_PFS_DAYS")] 
dc<-pre_merge_mat[,c(2:11,13,20:29)]

rownames(dc)<-pre_merge_mat[,1]
dc<-t(dc)

colnames(cox_clin)<-c("sample","Events","Survival")
cox_uni<-as.data.frame(getUniOrMultiCOXAnalysis(dc,cox_clin,method="Univariate"))
cox_uni$Log_likelihood <- log(as.numeric(cox_uni$`Likelihood ratio test`))
write.xlsx(cox_uni, "new_cox_uni_model_LLR.xlsx")

for (i in c(3,6,7,8)){cox_uni[,i] <- as.numeric(cox_uni[,i])}
select <- cox_uni[cox_uni$geneName %in% c("TMB_group","new_PDL1","Epithelial_RS"),c(1,3,7,8,6)]
select$CI <- paste0(round(select$lower.95,3),"-",round(select$upper.95,3))
select$`Pr(>|Z|)[pval]` <- round(select$`Pr(>|Z|)[pval]`,3)
select<-rbind(c("Features", NA, NA, NA,"p","HR(95%CI)"),select)
for(i in 2:4) {select[, i] = as.numeric(select[, i])}
pdf("feature_PFS_uni_forestplot.pdf",onefile = FALSE,width = 8,height = 5)
forestplot(select[,c(1,6,5)],  #1,6,5列显示为变量 HR(CI) p数值形式
           mean=select[,2],   #表格第3列为HR，要变成森林图的小方块
           lower=select[,3],  #表格第7列为5%CI，
           upper=select[,4],  #表格第8列为95%CI，它俩要化作线段，穿过方块
           zero=1,            #零线或参考线为HR=1即x轴的垂直线
           boxsize=0.2,       #设置小黑块的大小
           graph.pos=4,#"right",#森林图放在最右侧
           hrzl_lines=list("1" = gpar(lty=1,lwd=2),
                           "2" = gpar(lty=2),
                           "5"= gpar(lwd=2,lty=1)),
           graphwidth = unit(.25,"npc"),
           # xlab="HR (95% CI)",
           #xticks=c(0,1,2,4,6,8,10,) ,
           #----------------#线条粗细（x轴、置信区间）
           lwd.zero=1,
           lwd.ci=1.5,
           lwd.xaxis=1, 
           lty.ci=1,
           ci.vertices =T,
           ci.vertices.height=0.2, 
           clip=c(0.1,10),
           #----------------#行间距、字间距/box形状  
           ineheight=unit(8, 'mm'), 
           line.margin=unit(8, 'mm'),
           colgap=unit(6, 'mm'),
           col=fpColors(zero = "gray",
                        box = '#3CB371', 
                        lines = '#1a37d9'),
           fn.ci_norm="fpDrawCircleCI", 
           title="Univariate cox forest plot") 
dev.off()

###
dc<-pre_merge_mat[,c("TMB_group","new_PDL1","Epithelial_RS")]
rownames(dc)<-pre_merge_mat[,1]
dc<-t(dc)

cox_mul<-getUniOrMultiCOXAnalysis(dc,cox_clin,method="Multivariate")
cox_mul <- cox_mul[,c(1,3,7,8,6)]
cox_mul$CI <- paste0(round(cox_mul$lower.95,3),"-",round(cox_mul$upper.95,3))
cox_mul$`Pr(>|Z|)[pval]` <- round(cox_mul$`Pr(>|Z|)[pval]`,3)
cox_mul<-rbind(c("Features", NA, NA, NA,"p","HR(95%CI)"),cox_mul)
for(i in 2:4) {cox_mul[, i] = as.numeric(cox_mul[, i])}
pdf("feature_PFS_mul_forestplot.pdf",onefile = FALSE,width = 8,height = 5)
forestplot(cox_mul[,c(1,6,5)],  #1,6,5列显示为变量 HR(CI) p数值形式
           mean=cox_mul[,2],   #表格第3列为HR，要变成森林图的小方块
           lower=cox_mul[,3],  #表格第7列为5%CI，
           upper=cox_mul[,4],  #表格第8列为95%CI，它俩要化作线段，穿过方块
           zero=1,            #零线或参考线为HR=1即x轴的垂直线
           boxsize=0.2,       #设置小黑块的大小
           graph.pos=4,#"right",#森林图放在最右侧
           hrzl_lines=list("1" = gpar(lty=1,lwd=2),
                           "2" = gpar(lty=2),
                           "5"= gpar(lwd=2,lty=1)),
           graphwidth = unit(.25,"npc"),
           # xlab="HR (95% CI)",
           #xticks=c(0,1,2,4,6,8,10,) ,
           #----------------#线条粗细（x轴、置信区间）
           lwd.zero=1,
           lwd.ci=1.5,
           lwd.xaxis=1, 
           lty.ci=1,
           ci.vertices =T,
           ci.vertices.height=0.2, 
           clip=c(0.1,10),
           #----------------#行间距、字间距/box形状  
           ineheight=unit(8, 'mm'), 
           line.margin=unit(8, 'mm'),
           colgap=unit(6, 'mm'),
           col=fpColors(zero = "gray",
                        box = '#3CB371', 
                        lines = '#1a37d9'),
           fn.ci_norm="fpDrawCircleCI", 
           title="Multivariate cox forest plot") 
dev.off()