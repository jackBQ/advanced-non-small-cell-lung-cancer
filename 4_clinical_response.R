library(openxlsx)
library(survival)
library(survminer)
library(plyr)
library(maftools)
library(tidyr)

#clin_input----
clin_final <- read.xlsx("clin_pdl1_response.xlsx", sheet = 1)
table(clin_final$RECIST)
# Nx     PD     PR     SD 待评估 
# 9     12    107     32      2

clin_final$Response <- ifelse(clin_final$RECIST == "PR", "Response",
                              ifelse(clin_final$RECIST == "SD"|clin_final$RECIST == "PD", "nonResponse", NA))
table(clin_final$Response)
# nonResponse    Response 
# 44         107 

#PDL1----
clin_final$OS_STATUS <- ifelse(clin_final$OS_STATUS == "Yes", 1, 0)
clin_final$PFS_STATUS <- ifelse(clin_final$PFS_STATUS == "Yes", 1, 0)
clin_final$TMB_group <- ifelse(clin_final$Report_TMB >= 10, "TMB-H","TMB-L")
clin_final$new_PDL1 <- factor(clin_final$new_PDL1, levels = c("<1%", "1%-49%", ">=50%"))
clin_final$PDL1_0513 <- ifelse(clin_final$new_PDL1 == "<1%", 0,
                               ifelse(clin_final$new_PDL1 == "1%-49%", 1,2))
fit2 <- survfit(Surv(OS_DAYS, OS_STATUS) ~ new_PDL1, data = clin_final)
p <- ggsurvplot(fit2, data = clin_final,font.title=c(20, "bold"),conf.int=F, pval=TRUE, legend.title ="",ylab=c("Overall Survival"),surv.median.line = "hv",
                legend.labs = c("<1%", "1%-49%", "≥50%"), palette = c("#8c57a2b2","#2f509eb2", "#3da873b2"))
# clin_final$new_PDL1_group <- ifelse(clin_final$new_PDL1 == ">=50%",0,ifelse(clin_final$new_PDL1 == "1%-49%", 1,2))
res_cox<-coxph(Surv(OS_DAYS, OS_STATUS) ~PDL1_0513, data=clin_final)
p$plot = p$plot + 
  ggplot2::annotate("text",x = 100, y = 0.12,
                    label = paste("Hazard ratio = ",round(summary(res_cox)$conf.int[1],2))) + 
  ggplot2::annotate("text",x = 100, y = 0.05,
                    label = paste("95%CI:",round(summary(res_cox)$conf.int[3],2),"-",round(summary(res_cox)$conf.int[4],2),sep = ""))
pdf("PDL1_OS_survival.pdf", width = 5, height = 5, onefile = F)
p
dev.off()

fit2 <- survfit(Surv(PFS_DAYS, PFS_STATUS) ~ new_PDL1, data = clin_final)
p <- ggsurvplot(fit2, data = clin_final,font.title=c(20, "bold"),conf.int=F, pval=TRUE, legend.title ="",ylab=c("Progression Free Survival"),surv.median.line = "hv",
                legend.labs = c("<1%", "1%-49%", "≥50%"), palette = c("#8c57a2b2","#2f509eb2", "#3da873b2"))
res_cox<-coxph(Surv(PFS_DAYS, PFS_STATUS) ~PDL1_0513, data=clin_final)
p$plot = p$plot + 
  ggplot2::annotate("text",x = 100, y = 0.12,
                    label = paste("Hazard ratio = ",round(summary(res_cox)$conf.int[1],2))) + 
  ggplot2::annotate("text",x = 100, y = 0.05,
                    label = paste("95%CI:",round(summary(res_cox)$conf.int[3],2),"-",round(summary(res_cox)$conf.int[4],2),sep = ""))
pdf("PDL1_PFS_survival.pdf", width = 5, height = 5, onefile = F)
p
dev.off()

#TMB----
clin_final$TMB_group <- ifelse(clin_final$TMB_group == "TMB-H", 1, 0)
fit2 <- survfit(Surv(OS_DAYS, OS_STATUS) ~ TMB_group, data = clin_final)
p <- ggsurvplot(fit2, data = clin_final,font.title=c(20, "bold"),conf.int=F, pval=TRUE, legend.title ="",ylab=c("Overall Survival"),surv.median.line = "hv",
                legend.labs = c("TMB-L","TMB-H"), palette = c("#358DB9D2","#dc9445b2"))
res_cox<-coxph(Surv(OS_DAYS, OS_STATUS) ~TMB_group, data=clin_final)
p$plot = p$plot + 
  ggplot2::annotate("text",x = 100, y = 0.12,
                    label = paste("Hazard ratio = ",round(summary(res_cox)$conf.int[1],2))) + 
  ggplot2::annotate("text",x = 100, y = 0.05,
                    label = paste("95%CI:",round(summary(res_cox)$conf.int[3],2),"-",round(summary(res_cox)$conf.int[4],2),sep = ""))
pdf("TMB_OS_survival.pdf", width = 5, height = 5, onefile = F)
p
dev.off()

fit2 <- survfit(Surv(PFS_DAYS, PFS_STATUS) ~ TMB_group, data = clin_final)
p <- ggsurvplot(fit2, data = clin_final,font.title=c(20, "bold"),conf.int=F, pval=TRUE, legend.title ="",ylab=c("Progression Free Survival"), surv.median.line = "hv",               
                legend.labs = c("TMB-L","TMB-H"), palette = c("#358DB9D2","#dc9445b2"))
res_cox<-coxph(Surv(PFS_DAYS, PFS_STATUS) ~TMB_group, data=clin_final)
p$plot = p$plot + 
  ggplot2::annotate("text",x = 100, y = 0.12,
                    label = paste("Hazard ratio = ",round(summary(res_cox)$conf.int[1],2))) + 
  ggplot2::annotate("text",x = 100, y = 0.05,
                    label = paste("95%CI:",round(summary(res_cox)$conf.int[3],2),"-",round(summary(res_cox)$conf.int[4],2),sep = ""))
pdf("TMB_PFS_survival.pdf", width = 5, height = 5, onefile = F)
p
dev.off()

#**Percent_Weight----
a <- data.frame(table(clin_final$Response,clin_final$new_PDL1))
a <- ddply(a,.(Var1),transform,percent=Freq/sum(Freq)*100)
a$label = paste0(sprintf("%.1f",a$percent),"%")
a$Var2 <- factor(a$Var2, levels = c("<1%", "1%-49%", ">=50%"))
pvalue <- fisher.test(table(clin_final$Response,clin_final$new_PDL1))$p.value

pdf("PDL1_Response_percent.pdf", width = 4, height = 5, onefile = NULL)
ggplot(a,aes(Var1,percent,fill=Var2, label = label))+
  geom_bar(stat="identity",position = position_stack(), width = 0.8)+
  scale_fill_manual(values = c("#8c57a2b2","#2f509eb2", "#3da873b2"),label=c("<1%", "1%-49%", "≥50%"))+
  scale_y_continuous(labels = scales::percent_format(scale = 1))+ #百分比y轴
  labs(x=NULL,y="Percent Weight",
       fill="")+
  geom_text(size = 4, position = position_stack(vjust = 0.5))  +
  annotate(geom = "text",
           cex=6,
           x=1.5, y=105, # 根据自己的数据调节p value的位置
           label=paste0("P ", ifelse(pvalue<0.001, "< 0.001", paste0("= ",round(pvalue,3)))), # 添加P值
           color="black")+
  theme_classic()+
  theme(legend.position = "top",
        legend.text = element_text(size=12),
        axis.text = element_text(size=12),
        axis.title = element_text(size=12))
dev.off()

a <- data.frame(table(clin_final$Response,clin_final$TMB_group))
a <- ddply(a,.(Var1),transform,percent=Freq/sum(Freq)*100)
a$label = paste0(sprintf("%.1f",a$percent),"%")
a$Var2 <- factor(a$Var2, levels = c("TMB-H","TMB-L"))
pvalue <- fisher.test(table(clin_final$Response,clin_final$TMB_group))$p.value

pdf("TMB_Response_percent.pdf", width = 4, height = 5, onefile = NULL)
ggplot(a,aes(Var1,percent,fill=Var2, label = label))+
  geom_bar(stat="identity",position = position_stack(), width = 0.8)+
  scale_fill_manual(values = c("#dc9445b2","#358DB9D2"),label=c("TMB-H","TMB-L"))+
  scale_y_continuous(labels = scales::percent_format(scale = 1))+ #百分比y轴
  labs(x=NULL,y="Percent Weight",
       fill="")+
  geom_text(size = 4, position = position_stack(vjust = 0.5))  +
  annotate(geom = "text",
           cex=6,
           x=1.5, y=105, # 根据自己的数据调节p value的位置
           label=paste0("P ", ifelse(pvalue<0.001, "< 0.001", paste0("= ",round(pvalue,3)))), # 添加P值
           color="black")+
  theme_classic()+
  theme(legend.position = "top",
        legend.text = element_text(size=12),
        axis.text = element_text(size=12),
        axis.title = element_text(size=12))
dev.off()

#RTKRAS_analysis----
rtkras_mut <- read.xlsx("summary_var.xlsx", sheet = 1)
clin_final <- merge(clin_final, rtkras_mut[,c("Sample_ID","RTK-RAS_pathway")], by = "Sample_ID")
colnames(clin_final)[209] <- c("RTKRAS")
# clin_final$RTKRAS <- factor(clin_final$RTKRAS, levels = c(1,0))
clin_rtkras <- clin_final[!is.na(clin_final$Response),]

fit2 <- survfit(Surv(OS_DAYS, OS_STATUS) ~ RTKRAS, data = clin_final)
p <- ggsurvplot(fit2, data = clin_final,font.title=c(20, "bold"),conf.int=F, pval=TRUE, legend.title ="",ylab=c("Overall Survival"),surv.median.line = "hv",
                legend.labs = c("RTKRAS_WT","RTKRAS_MT"),palette = c("#00828099","#bb002199"))
res_cox<-coxph(Surv(OS_DAYS, OS_STATUS) ~RTKRAS, data=clin_final)
p$plot = p$plot + 
  ggplot2::annotate("text",x = 100, y = 0.12,
                    label = paste("Hazard ratio = ",round(summary(res_cox)$conf.int[1],2))) + 
  ggplot2::annotate("text",x = 100, y = 0.05,
                    label = paste("95%CI:",round(summary(res_cox)$conf.int[3],2),"-",round(summary(res_cox)$conf.int[4],2),sep = ""))
pdf("RTKRAS_OS_survival.pdf", width = 5, height = 5, onefile = F)
p
dev.off()

fit2 <- survfit(Surv(PFS_DAYS, PFS_STATUS) ~ RTKRAS, data = clin_final)
p <- ggsurvplot(fit2, data = clin_final,font.title=c(20, "bold"),conf.int=F, pval=TRUE, legend.title ="",ylab=c("Progression Free Survival"),  surv.median.line = "hv",              
                legend.labs = c("RTKRAS_WT","RTKRAS_MT"),palette = c("#00828099","#bb002199"))
res_cox<-coxph(Surv(PFS_DAYS, PFS_STATUS) ~RTKRAS, data=clin_final)
p$plot = p$plot + 
  ggplot2::annotate("text",x = 100, y = 0.12,
                    label = paste("Hazard ratio = ",round(summary(res_cox)$conf.int[1],2))) + 
  ggplot2::annotate("text",x = 100, y = 0.05,
                    label = paste("95%CI:",round(summary(res_cox)$conf.int[3],2),"-",round(summary(res_cox)$conf.int[4],2),sep = ""))
pdf("RTKRAS_PFS_survival.pdf", width = 5, height = 5, onefile = F)
p
dev.off()

#####
a <- data.frame(table(clin_rtkras$Response, clin_rtkras$RTKRAS))
a <- ddply(a,.(Var1),transform,percent=Freq/sum(Freq)*100)
a$label = paste0(sprintf("%.1f",a$percent),"%")
# a$Var2 <- factor(a$Var2, levels = c(,0))
pvalue <- fisher.test(table(clin_rtkras$Response, clin_rtkras$RTKRAS))$p.value

pdf("Response_RTKRAS_percent.pdf", width = 5, height = 5, onefile = NULL)
ggplot(a,aes(Var1,percent,fill=Var2, label = label))+
  geom_bar(stat="identity",position = position_stack())+
  scale_fill_manual(values = c("#bb002199","#00828099"),label=c("MT","WT"))+
  scale_y_continuous(labels = scales::percent_format(scale = 1))+ #百分比y轴
  labs(x=NULL,y="Percent Weight",
       fill="RTK-RAS")+
  geom_text(size = 4, position = position_stack(vjust = 0.5))  +
  annotate(geom = "text",
           cex=4,
           x=1.5, y=105, # 根据自己的数据调节p value的位置
           label=paste0("P ", ifelse(pvalue<0.001, "< 0.001", paste0("= ",round(pvalue,3)))), # 添加P值
           color="black")+
  theme_classic()+
  theme(legend.position = "top",
        legend.text = element_text(size=12),
        axis.text = element_text(size=12),
        axis.title = element_text(size=12))
dev.off()

#####
mut_filter <- read.table("somatic_mutation_nonsys.maf", sep = "\t", quote = '', header = T, row.names = NULL, check.names = F)
response_maf <- read.maf(mut_filter[mut_filter$Tumor_Sample_Barcode %in% clin_rtkras[clin_rtkras$Response == "Response",]$Sample_ID,])
non_response_maf <- read.maf(mut_filter[mut_filter$Tumor_Sample_Barcode %in% clin_rtkras[clin_rtkras$Response == "nonResponse",]$Sample_ID,])

pathway_nonPD<-OncogenicPathways(maf = non_response_maf)
pathway_PD<-OncogenicPathways(maf = response_maf)
forest_data<-merge(as.data.frame(pathway_nonPD[,c(1,5)]),
                   as.data.frame(pathway_PD[,c(1,5)]),by="Pathway",all.x=T)

forest_data[is.na(forest_data)]=0
colnames(forest_data)[2:3]<-c("nonPD_MUT","PD_MUT")

nonPD.sampleSize<-as.numeric(non_response_maf@summary[3, summary])
PD.sampleSize<-as.numeric(response_maf@summary[3, summary])
forest_data$nonPD_Wild<- 44 - forest_data$nonPD_MUT
forest_data$PD_Wild<- 107 - forest_data$PD_MUT
rownames(forest_data)<-forest_data$Pathway
forest_data<-forest_data[,-1]
forest_data<-forest_data[,c(2,4,1,3)]
#*fisher检验-----
fisher_in_bulk <- function(pathway_list, mutation_info){
  fisher_result <- function(pathway){
    gene_data<-matrix(unlist(mutation_info[pathway,]), byrow = TRUE, nrow = 2)
    # 在创建matrix之前我们需要unlist函数 去解列 data.frame 为vector
    xf=fisher.test(gene_data, conf.int = TRUE, conf.level = 0.95)
    fisher_report <- data.frame('Pathway'=pathway,
                                'Response'=gene_data[1,1],
                                'nonResponse'=gene_data[2,1],
                                'OR'=round(xf$estimate,3),
                                '5%CI'=round(xf$conf.int[1],3),
                                '95%CI'=round(xf$conf.int[2],3),
                                'pvalue'=round(xf$p.value,3),
                                'OR(CI)'=paste0(round(xf$estimate,3),
                                                "(",round(xf$conf.int[1],3),
                                                "-",round(xf$conf.int[2],3),")"))
    fisher_report  
  }
  fisher_result_list <- lapply(pathway_list, fisher_result)
  do.call(rbind, fisher_result_list)
}
pathway_fisher <- fisher_in_bulk(pathway_list = rownames(forest_data),
                                 mutation_info = forest_data)
pathway_barpolt <- gather(pathway_fisher[,c(1:3)], Response, MT, -Pathway)

pdf("Pathway_Response.pdf",width = 8, height = 6, onefile = F)
ggplot(data = pathway_barpolt,aes(x=Pathway,y=MT,fill=Response,label = MT))+
  geom_bar(stat = "identity",position = "stack")+
  theme_bw()+theme(panel.grid = element_blank(),legend.position = "top")+
  scale_fill_manual(values = c("Response" = "#A8450C93", "nonResponse" = "#16545993"))+
  annotate("text",x=1,y=135,label="0.700",size=4)+
  annotate("text",x=2,y=135,label="1.000",size=4)+
  annotate("text",x=3,y=135,label="0.373",size=4)+
  annotate("text",x=4,y=135,label="0.718",size=4)+
  annotate("text",x=5,y=135,label="0.213",size=4)+
  annotate("text",x=6,y=135,label="0.716",size=4)+
  annotate("text",x=7,y=135,label="0.034",size=4)+
  annotate("text",x=8,y=135,label="0.282",size=4)+
  annotate("text",x=9,y=135,label="0.630",size=4)+
  annotate("text",x=10,y=135,label="0.534",size=4)+
  geom_text(size = 4, position = position_stack(vjust = 0.5))  
dev.off()

#####
#    HRS = "#B22222"   LRS = "#000080"
var <- clin_final[,c("Sample_ID","Report_TMB","new_PDL1","ARID1A","RTKRAS","RECIST")]
date <- read.xlsx("new_最后随访-更新疗效与生存 发liu 2024-3-22.xlsx", sheet = 1)
date_filter <- date[date$患者姓名 %in% clin_pre$患者姓名,]
data_filter <- merge(date_filter, clin_pre[,c("Sample_ID","患者姓名")], by = "患者姓名")

var_date <- merge(var, data_filter[,c("Sample_ID", "首次PR/CR日期")], by = "Sample_ID")
var_date$RTKRAS <- ifelse(var_date$RTKRAS == 1, "MT", "WT")
write.xlsx(var_date, "select_clin_matrix.xlsx", rowNames = F)

#####
lasso <- read.table("var_lasso_group.txt", sep = '\t', quote = '', check.names = F, header = T, row.names = NULL)
identical(lasso$Sample_ID, clin_final$Sample_ID)
clin_final$RS <- lasso$RS
clin_forest <- clin_final[,c("Sample_ID","PFS_STATUS","PFS_DAYS","Age","Gender","Smoking","Drinking","Differentiated","Pathological_subtype","Stage","new_PDL1","TMB_group", "PD1.PDL1","RECIST","RS")]
colnames(clin_forest) <- c("Sample_ID","PFS_STATUS","PFS_DAYS","Age","Sex","Smoking","Drinking","Differentiated","Pathological_subtype","Stage","PDL1","TMB", "First_line_Therapy","BOR","RS")
clin_forest <- clin_forest[!is.na(clin_forest$PFS_STATUS),]
clin_forest$BOR <- ifelse(clin_forest$BOR == "CR"|clin_forest$BOR == "PR", "Response", 
                          ifelse(clin_forest$BOR == "SD"|clin_forest$BOR == "PD","nonResponse", NA))

clin_forest$BOR <- factor(clin_forest$BOR, levels = c("Response", "nonResponse"))
clin_forest$Differentiated <- factor(clin_forest$Differentiated, levels = c("Moderately", "Low", "High")) 
clin_forest[is.na(clin_forest$First_line_Therapy),]$First_line_Therapy <- c("Chemotherapy+checkpoint inhibitor anticancer drugs")
clin_forest$First_line_Therapy <- ifelse(clin_forest$First_line_Therapy == "PD1", "Chemotherapy+PD-1 inhibitors",
                                         ifelse(clin_forest$First_line_Therapy == "PDL1", "Chemotherapy+PDL1 inhibitors", "Chemotherapy+checkpoint inhibitor anticancer drugs")) 
clin_forest$RS <- factor(clin_forest$RS, levels = c("LRS","HRS"))
clin_forest$First_line_Therapy <- factor(clin_forest$First_line_Therapy,levels = c("Chemotherapy+PD-1 inhibitors", "Chemotherapy+PDL1 inhibitors", "Chemotherapy+checkpoint inhibitor anticancer drugs"))

###uni_cox
y = Surv(time = clin_forest$PFS_DAYS, event = clin_forest$PFS_STATUS == 1)
Uni_cox_model<- function(x){
  FML <- as.formula(paste0 ("y~",x))
  feature = x
  cox<- coxph(FML,data=clin_forest)
  cox1<-summary(cox)
  HR <- round(cox1$coefficients[,2],3)
  PValue <- round(cox1$coefficients[,5],3)
  CI5 <-round(cox1$conf.int[,3],3)
  CI95 <-round(cox1$conf.int[,4],3)
  subchar <- rownames(cox1$coefficients)
  Uni_cox_model<- data.frame("Feature" = x,
                             'Characteristics' = subchar,
                             'HR' = HR,
                             'Lower' = CI5,
                             'Upper' = CI95,
                             'Pvalue' = PValue)
  return(Uni_cox_model)
}  
variable.names<- colnames(clin_forest)[c(4:ncol(clin_forest))]
Uni_cox <- lapply(variable.names, Uni_cox_model)
Uni_cox<- ldply(Uni_cox,data.frame)
Uni_cox$HR.CI <- paste0(Uni_cox$HR,"(",Uni_cox$Lower,"-",Uni_cox$Upper, ")")
for (i in 1:nrow(Uni_cox)){
  Uni_cox[i,"Characteristics"] <- str_sub(Uni_cox[i,"Characteristics"], nchar(Uni_cox[i, "Feature"])+1, -1)
}
result <- Uni_cox

ins_1 <- function(x) {c(x, rep(NA, ncol(result)-1))}
ins_2 <- function(x) {c(NA,x, rep(NA, ncol(result)-3), "Reference")}

##2-2：插入空行，形成一个新表
for(i in 6:7) {result[, i] = as.character(result[, i])}

result[1,2] <- c("Age")
result$Feature <- NA
result<-rbind(c("Feature","Characteristics", NA, NA, NA, "p", "HR(95%CI)"),
              result[1, ],
              ins_1("Sex"),
              ins_2("Female"),  
              result[2, ], 
              ins_1("Smoking"),
              ins_2("No"),
              result[3, ],
              ins_1("Drinking"),
              ins_2("No"),
              result[4, ],
              ins_1("Differentiated"),
              ins_2("Moderately"),
              result[5:6, ],
              ins_1("Pathological_subtype"),
              ins_2("LUAD"),
              result[7:8,],
              ins_1("Stage"),
              ins_2("III"),
              result[9,],
              ins_1("PDL1"),
              ins_2("<1%"),
              result[10:11,],
              ins_1("TMB"),
              ins_2("TMB-H"),
              result[12,],
              ins_1("First_line_Therapy"),
              ins_2("Chemotherapy+PD-1 inhibitors"),
              result[13:14,],
              ins_1("Response"),
              ins_2("Response"),
              result[15,],
              ins_1("RS"),
              ins_2("LRS"),
              result[16,],
              c(NA, NA, NA, NA, NA,NA))
# for(i in 2:4) {result[, i] = as.numeric(result[, i])}
colnames(result) <- c("Feature","Characteristics", NA, NA, NA, "p", "HR(95%CI)")

pdf("new_pfs_uni_clin_forest.pdf", width = 11, height = 8, onefile = F)
forestplot(
  result[,c(1,2,7,6)],        # 需要显示在森林图中的列
  mean = as.numeric(result[, 3]),             # 均值列（HR），它将显示为森林图的小方块或其他形状哈哈哈哈哈
  lower = as.numeric(result[, 4]),            # 95%置信区间的下限数据列
  upper = as.numeric(result[, 5]),            # 95%置信区间的上限数据列
  zero=1,            #零线或参考线为HR=1即x轴的垂直线
  boxsize=0.3,       #设置小黑块的大小
  lineheight=unit(6,'mm'),graph.pos=5,#"right",#森林图放在最右侧
  #graphwidth = unit(.25,"npc"),
  #----------------#线条粗细（x轴、置信区间）
  lwd.zero=1,#线宽度
  clip=c(-1,10), # 设置x轴的范围，若置信区间落在设定的范围外，则用箭头表示
  lwd.ci=2,lwd.xaxis=1, 
  lty.ci=1,# HR线（穿过box的直线）的线型，默认为1（实直线）
  #箱线图两端添加小竖线，高度
  ci.vertices =T,ci.vertices.height=0.2, 
  #----------------#行间距、字间距 
  ineheight=unit(8, 'mm'), line.margin=unit(8, 'mm'),colgap=unit(6, 'mm'),
  col=fpColors(zero = "grey", box = 'green', lines = 'blue'),
  fn.ci_norm="fpDrawCircleCI")
dev.off()

###mul_cox
mul_cox_model<- as.formula(paste0 ("y~",paste0(Uni_cox$Feature[Uni_cox$Pvalue<=1],collapse = "+")))
mul_cox<-coxph(mul_cox_model,data=clin_forest)
cox4<-summary(mul_cox)
mul_HR<- round(cox4$coefficients[,2],3)
mul_PValue<- round(cox4$coefficients[,5],3)
mul_CI1<-round(cox4$conf.int[,3],3)
mul_CI2<-round(cox4$conf.int[,4],3)
mut_HRCI <- paste0(round(cox4$coefficients[,2],3),"(",round(cox4$conf.int[,3],3),"-",round(cox4$conf.int[,4],3),")")
mul_CI<-paste(mul_CI1,'-',mul_CI2)
mul_cox1<- data.frame("P"=mul_PValue, "HR.CI" = mut_HRCI, "HR"=mul_HR,"lower"=mul_CI1, "upper" = mul_CI2 )
mul_cox1 <- mul_cox1[-14,]
mul_cox1$Characteristics <- c("Age", "Male", "Yes", "Yes", "Low", "High", "LUSC", "Other", "IV", "1%-49%", ">=50%", "TMB-L", "Chemotherapy+PDL1 inhibitors", "nonResponse", "HRS")
mul_cox1$Feature <- NA
mul_cox1 <- mul_cox1[,c(7,6,3,4,5,1,2)]

mul_cox1<-rbind(c("Feature","Characteristics", NA, NA, NA, "p", "HR(95%CI)"),
                mul_cox1[1, ],
                ins_1("Sex"),
                ins_2("Female"),  
                mul_cox1[2, ], 
                ins_1("Smoking"),
                ins_2("No"),
                mul_cox1[3, ],
                ins_1("Drinking"),
                ins_2("No"),
                mul_cox1[4, ],
                ins_1("Differentiated"),
                ins_2("Moderately"),
                mul_cox1[5:6, ],
                ins_1("Pathological_subtype"),
                ins_2("LUAD"),
                mul_cox1[7:8,],
                ins_1("Stage"),
                ins_2("III"),
                mul_cox1[9,],
                ins_1("PDL1"),
                ins_2("<1%"),
                mul_cox1[10:11,],
                ins_1("TMB"),
                ins_2("TMB-H"),
                mul_cox1[12,],
                ins_1("First_line_Therapy"),
                ins_2("Chemotherapy+PD-1 inhibitors"),
                mul_cox1[13,],
                ins_1("Response"),
                ins_2("Response"),
                mul_cox1[14,],
                ins_1("RS"),
                ins_2("LRS"),
                mul_cox1[15,],
                c(NA, NA, NA, NA, NA,NA))

pdf("new_pfs_mul_clin_forest.pdf", width = 9, height = 8, onefile = F)
forestplot(
  mul_cox1[,c(1,2,7,6)],        # 需要显示在森林图中的列
  mean = as.numeric(mul_cox1[, 3]),             # 均值列（HR），它将显示为森林图的小方块或其他形状哈哈哈哈哈
  lower = as.numeric(mul_cox1[, 4]),            # 95%置信区间的下限数据列
  upper = as.numeric(mul_cox1[, 5]),            # 95%置信区间的上限数据列
  zero=1,            #零线或参考线为HR=1即x轴的垂直线
  boxsize=0.3,       #设置小黑块的大小
  lineheight=unit(6,'mm'),graph.pos=5,#"right",#森林图放在最右侧
  #graphwidth = unit(.25,"npc"),
  #----------------#线条粗细（x轴、置信区间）
  lwd.zero=1,#线宽度
  clip=c(-1,10), # 设置x轴的范围，若置信区间落在设定的范围外，则用箭头表示
  lwd.ci=2,lwd.xaxis=1, 
  lty.ci=1,# HR线（穿过box的直线）的线型，默认为1（实直线）
  #箱线图两端添加小竖线，高度
  ci.vertices =T,ci.vertices.height=0.2, 
  #----------------#行间距、字间距 
  ineheight=unit(8, 'mm'), line.margin=unit(8, 'mm'),colgap=unit(6, 'mm'),
  col=fpColors(zero = "grey", box = 'green', lines = 'blue'),
  fn.ci_norm="fpDrawCircleCI")
dev.off()