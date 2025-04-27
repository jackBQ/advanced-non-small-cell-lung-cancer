library(openxlsx)
sample_info<-read.xlsx("new_clin.xlsx")
var_data<-read.table("new_var_filter.txt",header = T,quote = "",sep = "\t")
clin<-sample_info[,c("Sample_ID","new_PDL1","OS_DAYS",
                     "OS_STATUS","PFS_DAYS","PFS_STATUS")]
clin<-merge(clin,var_data[,c("Sample_ID","KDM6A","TMB_group","RECIST")])
clin$TMB_group<-ifelse(clin$TMB_group == 1,"TMB-H","TMB-L")
clin$new_PDL1<-gsub("<1%|1%-49%","<50%",clin$new_PDL1)
clin$KDM6A<-ifelse(clin$KDM6A == 1,"MT","WT")
clin$OS_STATUS<-gsub("Yes",1,clin$OS_STATUS)
clin$OS_STATUS<-gsub("No",0,clin$OS_STATUS)
clin$OS_STATUS<-as.numeric(clin$OS_STATUS)
clin$PFS_STATUS<-gsub("No",0,clin$PFS_STATUS)
clin$PFS_STATUS<-gsub("Yes",1,clin$PFS_STATUS)
clin$PFS_STATUS<-as.numeric(clin$PFS_STATUS)

pvalue <- fisher.test(table(clin$TMB_group,clin$KDM6A))$p.value#方差检验
TMB_data<-data.frame(table(clin$TMB_group,clin$KDM6A))
TMB_data<- ddply(TMB_data,.(Var2),transform,percent=Freq/sum(Freq)*100) 
TMB_data$label = paste0(sprintf("%.1f", TMB_data$percent), "%")
colnames(TMB_data)[1]<-"TMB"
TMB_plot<-ggplot(TMB_data , aes(x=Var2,y=percent,fill=TMB)) + 
  geom_bar(position = position_fill(),width = 0.8,stat = "identity") + 
  labs(y = 'Percent Weight',x='',fill="") +#设置y轴名为'Percent'
  scale_fill_manual(values = c("TMB-H"="#DC9445","TMB-L"="#358DB9"))+
  scale_y_continuous(labels = scales::percent_format(scale = 100))+ #百分比y轴
  geom_text(aes(label=label),position = position_fill(vjust = 0.5),color="black")+
  annotate(geom = "text",
           cex=6,
           x=1.5, y=1.05, # 根据自己的数据调节p value的位置
           label=paste0("P ", ifelse(pvalue<0.001, "< 0.001", paste0("= ",round(pvalue,3)))), # 添加P值
           color="black")+
  theme_classic()+
  theme(legend.position = "top",
        legend.text = element_text(size=12),
        axis.text = element_text(size=12),
        axis.title = element_text(size=12))
#ggsave("TMB_plot_output.pdf",TMB_plot,height=5,width=4.5)

pvalue <- fisher.test(table(clin$new_PDL1,clin$KDM6A))$p.value
PDL1_data<-data.frame(table(clin$new_PDL1,clin$KDM6A))
PDL1_data<- ddply(PDL1_data,.(Var2),transform,percent=Freq/sum(Freq)*100) 
PDL1_data$label = paste0(sprintf("%.1f", PDL1_data$percent), "%")
colnames(PDL1_data)[1]<-"PDL1"
PDL1_plot<-ggplot(PDL1_data , aes(x=Var2,y=percent,fill=PDL1)) + 
  geom_bar(position = position_fill(),width = 0.8,stat = "identity") + 
  labs(y = 'Percent Weight',x='',fill="") +#设置y轴名为'Percent'
  scale_fill_manual(values = c(">=50%"="#3DA873", "<50%"="#8C57A2"))+
  scale_y_continuous(labels = scales::percent_format(scale = 100))+ #百分比y轴
  geom_text(aes(label=label),position = position_fill(vjust = 0.5),color="black")+
  annotate(geom = "text",
           cex=6,
           x=1.5, y=1.05, # 根据自己的数据调节p value的位置
           label=paste0("P ", ifelse(pvalue<0.001, "< 0.001", paste0("= ",round(pvalue,3)))), # 添加P值
           color="black")+
  theme_classic()+
  theme(legend.position = "top",
        legend.text = element_text(size=12),
        axis.text = element_text(size=12),
        axis.title = element_text(size=12))
#ggsave("PDL1_plot_output.pdf",PDL1_plot,height=5,width=4.5)

library(survival)
library(survminer)
###OS-----
clin1<-clin
clin1$KDM6A<-ifelse(clin1$KDM6A == "MT",1,0)
KDM6A_OS_cox<-summary(coxph(Surv(OS_DAYS, OS_STATUS) ~ KDM6A, data = clin1))
KDM6A_OS_kmfit<-survfit(Surv(OS_DAYS, OS_STATUS) ~ KDM6A, data = clin)
p<-ggsurvplot(KDM6A_OS_kmfit,
              pval = TRUE, #conf.int = TRUE,
              xlab = "Time (Days)",#ggtheme =theme_bw(),
              ylab = "Overall Survival",#ggtheme =theme_bw(),
              palette = c("MT" = "#BB0021" ,"WT" = "#008280"),
              legend = "top",#c(0.8,0.8), # 指定图例位置
              legend.title = "", # 设置图例标题
              legend.labs = sort(unique(clin$KDM6A)),
              risk.table = F,
              surv.median.line = "hv"
)
p$plot<-p$plot+
  annotate("text", x = 0, y = 0.08, 
           label = paste0("Hazard ratio = ",sprintf("%.3f",KDM6A_OS_cox$conf.int[1]),
                          "\n95% CI:",sprintf("%.3f",KDM6A_OS_cox$conf.int[3]),
                          " to ",sprintf("%.3f",KDM6A_OS_cox$conf.int[4]),
                          "\nMedian OS:",surv_median(KDM6A_OS_kmfit)$median[2]," Days vs ",
                          surv_median(KDM6A_OS_kmfit)$median[1]," Days"),   ###添加P和HR 95%CI
           size = 5, color = "black", hjust = 0)+
  theme(text = element_text(size = 15))
pdf("KDM6A_mutation_OS_output_xin.pdf",onefile = F)    
p
dev.off()

###PFS----
KDM6A_PFS_cox<-summary(coxph(Surv(PFS_DAYS, PFS_STATUS) ~ KDM6A, data = clin1))
KDM6A_PFS_kmfit<-survfit(Surv(PFS_DAYS, PFS_STATUS) ~ KDM6A, data = clin)
p<-ggsurvplot(KDM6A_PFS_kmfit,
              pval = TRUE, #conf.int = TRUE,
              xlab = "Time (Days)",
              ylab = "Progression Free Survival",#ggtheme =theme_bw(),
              palette = c("MT" = "#BB0021" ,"WT" = "#008280"),
              legend = "top",#c(0.8,0.8), # 指定图例位置
              legend.title = "", # 设置图例标题
              legend.labs = sort(unique(clin$KDM6A)),
              risk.table = F,
              surv.median.line = "hv",
)
p$plot<-p$plot+
  annotate("text", x = 0, y = 0.08, 
           label = paste0("Hazard ratio = ",sprintf("%.3f",KDM6A_PFS_cox$conf.int[1]),
                          "\n95% CI:",sprintf("%.3f",KDM6A_PFS_cox$conf.int[3]),
                          " to ",sprintf("%.3f",KDM6A_PFS_cox$conf.int[4]),
                          "\nMedian PFS:",surv_median(KDM6A_PFS_kmfit)$median[2]," Days vs ",
                          surv_median(KDM6A_PFS_kmfit)$median[1]," Days"),   ###添加P和HR 95%CI
           size = 5, color = "black", hjust = 0)+
  theme(text = element_text(size = 15))
pdf("KDM6A_mutation_PFS_output_xin.pdf",onefile = F)    
p
dev.off()

 
#MSK---
KM_data<-read.table("Overall_cbio .txt",quote = "",sep = "\t",fill = T)
KM_data<-KM_data[-c(1:2,14:15),]
colnames(KM_data)<-KM_data[1,]
KM_data<-KM_data[-1,]
KM_data$KDM6A<-c(rep("MT",10),rep("WT",728))
KM_data$OS<-ifelse(KM_data$Status == "deceased",1,0)
KM_data$OS_time<-as.numeric(KM_data$`Time (months)`)

KM_data1<-KM_data
KM_data1$KDM6A<-ifelse(KM_data1$KDM6A == "MT",1,0)
KDM6A_MSK_cox<-summary(coxph(Surv(OS_time, OS) ~ KDM6A, data = KM_data1))
KDM6A_MSK_kmfit<-survfit(Surv(OS_time, OS) ~ KDM6A, data = KM_data)
p<-ggsurvplot(KDM6A_MSK_kmfit,
              pval = TRUE, #conf.int = TRUE,
              xlab = "Time (Months)",
              ylab = "Overall Survival",#ggtheme =theme_bw(),
              palette = c("MT" = "#BB0021" ,"WT" = "#008280"),
              legend = "top",#c(0.8,0.8), # 指定图例位置
              legend.title = "", # 设置图例标题
              legend.labs = sort(unique(clin$KDM6A)),
              risk.table = F,
              surv.median.line = "hv",
)
p$plot<-p$plot+
  annotate("text", x = 0, y = 0.05, 
           label = paste0("Hazard ratio = ",sprintf("%.3f",KDM6A_MSK_cox$conf.int[1]),
                          "\n95% CI:",sprintf("%.3f",KDM6A_MSK_cox$conf.int[3]),
                          " to ",sprintf("%.3f",KDM6A_MSK_cox$conf.int[4]),
                          "\nMedian OS:",sprintf("%.1f",surv_median(KDM6A_MSK_kmfit)$median[2])," Months vs ",
                          sprintf("%.1f",surv_median(KDM6A_MSK_kmfit)$median[1])," Months"),   ###添加P和HR 95%CI
           size = 5, color = "black", hjust = 0)+
  theme(text = element_text(size = 15))
pdf("MSK_KDM6A_mutation_OS_output_xin1.pdf",onefile = F)    
p
dev.off()

#TMB
KDM6A_mutation<-read.table("nsclc_ctdx_msk_2022/data_mutations.txt",header = T,
                           quote = "",sep = "\t",comment.char = "#")
KDM6A_sample<-unique(KDM6A_mutation[which(KDM6A_mutation$Hugo_Symbol == "KDM6A"),"Tumor_Sample_Barcode"])
KDM6A_sample<-substr(KDM6A_sample,1,9)
KDM6A_TMB<-read.table("nsclc_ctdx_msk_2022/data_clinical_sample.txt",header = T,
                      quote = "",sep = "\t",comment.char = "#")
KDM6A_TMB$group<-ifelse(KDM6A_TMB$PATIENT_ID %in% KDM6A_sample,"MT","WT")
KDM6A_TMB1<-KDM6A_TMB#[which(KDM6A_TMB$PATIENT_ID %in% KM_data$`Case ID`),]
my_comparsion<-list(c("MT","WT"))
KDM6A_TMB_plot<-ggplot(KDM6A_TMB1,aes(x=group,y=TMB_NONSYNONYMOUS,fill = factor(group)))+
  geom_boxplot(width = 0.4,notch = F,position=position_dodge(0.8)) +
  geom_point(alpha=1, shape=21, size = 2.5,
             position=position_jitterdodge(dodge.width = 0.8))+ #添加数据点
  labs( y= 'TMB',x='KDM6A')+
  scale_fill_manual(values = c("MT" = "#BB0021" ,"WT" = "#008280"))+
  theme_classic()+
  theme(panel.grid=element_blank())+
  stat_compare_means(aes(group=group), method = "wilcox.test",label = "p.signif", label.x = 1.5,comparisons = my_comparsion)+
  guides(fill="none")
ggsave("MSK_KDM6A_boxplot_output_xin.pdf",KDM6A_TMB_plot,width = 5,height = 5)
