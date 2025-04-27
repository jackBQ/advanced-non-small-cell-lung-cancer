rm(list = ls())
#ssgsea----
load("ssgsea_LateNSCLC.Rdata")
LUAD_group<-read.table("TCGA_LUAD_cutoff_group.txt",header = T)
LUSC_group<-read.table("TCGA_LUSC_cutoff_group.txt",header = T)
NSCLC_group<-read.table("TCGA_NSCLC_cutoff_group.txt",header = T)
  
ssgsea_scores<-as.data.frame(t(ssgsea))
ssgsea_scores$PATIENT_ID<-substr(rownames(ssgsea_scores),1,15)
LUAD_DSS_ssgsea<-merge(ssgsea_scores,LUAD_group[,c("PATIENT_ID","DSS_group")])
LUSC_DSS_ssgsea<-merge(ssgsea_scores,LUSC_group[,c("PATIENT_ID","DSS_group")])
NSCLC_DSS_ssgsea<-merge(ssgsea_scores,NSCLC_group[,c("PATIENT_ID","DSS_group")])

ssgsea_scores1<-rbind(LUAD_DSS_ssgsea,LUSC_DSS_ssgsea,NSCLC_DSS_ssgsea)

library(reshape2)
ssgsea_DSS<-melt(ssgsea_scores1,id.vars=c("DSS_group"),
                 measure.vars=c(colnames(ssgsea_scores1)[2:(ncol(ssgsea_scores1)-1)]))

library(ggplot2)
library(ggpubr)
new_comparisons <- list(c("LUAD_HRS","LUAD_LRS"),c("LUSC_HRS","LUSC_LRS"),c("NSCLC_HRS","NSCLC_LRS"))
for (i in 2:(dim(ssgsea_scores1)[2]-1))  { 
  print(paste0("run to cells:",colnames(ssgsea_scores1)[i]," i=",i))
  pdf(file = paste0("ssgsea/",colnames(ssgsea_scores1)[i],"_ssgsea_output.pdf"),onefile=F)
  a<-ggplot(ssgsea_DSS[which(ssgsea_DSS$variable == colnames(ssgsea_scores1)[i]),],aes(x=DSS_group,y=value,fill = factor(DSS_group)))+
    geom_violin(aes(fill=DSS_group),position=position_dodge(0.8))+
    geom_boxplot(width = 0.2,notch = F,position=position_dodge(0.8)) +
    labs( y= colnames(ssgsea_scores1)[i],x='')+
    scale_fill_manual(values = c("NSCLC_HRS"="#F77300","NSCLC_LRS"="#00AF86",
                                 "LUAD_HRS"="#FF615E","LUAD_LRS"="#00B8CF",
                                 "LUSC_HRS"="#F7A034","LUSC_LRS"="#8661C1"
    ))+
    theme_classic()+
    theme(panel.grid=element_blank())+
    stat_compare_means(aes(group=DSS_group),label = "p.signif", label.x = 1.5,comparisons = new_comparisons)+
    guides(fill="none")
  print(a)
  dev.off()
}

LUAD_DFI_ssgsea<-merge(ssgsea_scores,LUAD_group[,c("PATIENT_ID","DFI_group")])
LUSC_DFI_ssgsea<-merge(ssgsea_scores,LUSC_group[,c("PATIENT_ID","DFI_group")])
NSCLC_DFI_ssgsea<-merge(ssgsea_scores,NSCLC_group[,c("PATIENT_ID","DFI_group")])

ssgsea_scores2<-rbind(LUAD_DFI_ssgsea,LUSC_DFI_ssgsea,NSCLC_DFI_ssgsea)

library(reshape2)
ssgsea_DFI<-melt(ssgsea_scores2,id.vars=c("DFI_group"),
                 measure.vars=c(colnames(ssgsea_scores2)[2:(ncol(ssgsea_scores2)-1)]))

library(ggplot2)
library(ggpubr)
new_comparisons <- list(c("LUAD_HRS","LUAD_LRS"),c("LUSC_HRS","LUSC_LRS"),c("NSCLC_HRS","NSCLC_LRS"))
for (i in 2:(dim(ssgsea_scores2)[2]-1))  { 
  print(paste0("run to cells:",colnames(ssgsea_scores2)[i]," i=",i))
  pdf(file = paste0("ssgsea/",colnames(ssgsea_scores2)[i],"_ssgsea_output.pdf"),onefile=F)
  a<-ggplot(ssgsea_DFI[which(ssgsea_DFI$variable == colnames(ssgsea_scores2)[i]),],aes(x=DFI_group,y=value,fill = factor(DFI_group)))+
    geom_violin(aes(fill=DFI_group),position=position_dodge(0.8))+
    geom_boxplot(width = 0.2,notch = F,position=position_dodge(0.8)) +
    labs( y= colnames(ssgsea_scores2)[i],x='')+
    scale_fill_manual(values = c("NSCLC_HRS"="#F77300","NSCLC_LRS"="#00AF86",
                                 "LUAD_HRS"="#FF615E","LUAD_LRS"="#00B8CF",
                                 "LUSC_HRS"="#F7A034","LUSC_LRS"="#8661C1"
    ))+
    theme_classic()+
    theme(panel.grid=element_blank())+
    stat_compare_means(aes(group=DFI_group),label = "p.signif", label.x = 1.5,comparisons = new_comparisons)+
    guides(fill="none")
  print(a)
  dev.off()
}


#免疫检查点-------
Immunoinhibitors_genes<-read.table("Immunoinhibitors.txt",
                                   header = T)
Immunostimulators_genes<-read.table("Immunostimulators.txt",
                                    header = T)
HLA_genes<-read.table("HLA基因.txt",header=F)
exp_NSCLC_TPM<-read.table("NSCLC_exp_TPM_975_output.txt",
                          header = T,row.names = 1,check.names = F,quote = "",
                          sep = "\t")
exp_NSCLC_TPM_log<-log2(exp_NSCLC_TPM+1)
exp_NSCLC_TPM_log<-as.data.frame(t(exp_NSCLC_TPM_log))
exp_NSCLC_TPM_log$PATIENT_ID<-substr(rownames(exp_NSCLC_TPM_log),1,15)
Immunoinhibitors_exp<-exp_NSCLC_TPM_log[,which(colnames(exp_NSCLC_TPM_log) %in% c(Immunoinhibitors_genes[,1],"PATIENT_ID"))]
Immunostimulators_exp<-exp_NSCLC_TPM_log[,which(colnames(exp_NSCLC_TPM_log) %in% c(Immunostimulators_genes[,1],"PATIENT_ID"))]
HLA_exp<-exp_NSCLC_TPM_log[,which(colnames(exp_NSCLC_TPM_log) %in% c(HLA_genes[,1],"PATIENT_ID"))]

#Immunoinhibitors----
LUAD_DSS_Immunoinhibitors<-merge(Immunoinhibitors_exp,LUAD_group[,c("PATIENT_ID","DSS_group")])
LUSC_DSS_Immunoinhibitors<-merge(Immunoinhibitors_exp,LUSC_group[,c("PATIENT_ID","DSS_group")])
NSCLC_DSS_Immunoinhibitors<-merge(Immunoinhibitors_exp,NSCLC_group[,c("PATIENT_ID","DSS_group")])

Immunoinhibitors_scores1<-rbind(LUAD_DSS_Immunoinhibitors,LUSC_DSS_Immunoinhibitors,NSCLC_DSS_Immunoinhibitors)

library(reshape2)
Immunoinhibitors_data1<-melt(Immunoinhibitors_scores1,id.vars=c("DSS_group"),
                             measure.vars=c(colnames(Immunoinhibitors_scores1)[2:(ncol(Immunoinhibitors_scores1)-1)]))

library(ggplot2)
library(ggpubr)
#new_comparisons <- list(c("LUAD_HRS","LUAD_LRS"),c("LUSC_HRS","LUSC_LRS"),c("NSCLC_HRS","NSCLC_LRS"))
for (i in 2:(dim(Immunoinhibitors_scores1)[2]-1))  { 
  print(paste0("run to cells:",colnames(Immunoinhibitors_scores1)[i]," i=",i))
  pdf(file = paste0("Immunoinhibitors/",colnames(Immunoinhibitors_scores1)[i],"_Immunoinhibitors_output.pdf"),onefile=F)
  a<-ggplot(Immunoinhibitors_data1[which(Immunoinhibitors_data1$variable == colnames(Immunoinhibitors_scores1)[i]),],aes(x=DSS_group,y=value,fill = factor(DSS_group)))+
    geom_violin(aes(fill=DSS_group),position=position_dodge(0.8))+
    geom_boxplot(width = 0.2,notch = F,position=position_dodge(0.8)) +
    labs( y= colnames(Immunoinhibitors_scores1)[i],x='')+
    scale_fill_manual(values = c("NSCLC_HRS"="#F77300","NSCLC_LRS"="#00AF86",
                                 "LUAD_HRS"="#FF615E","LUAD_LRS"="#00B8CF",
                                 "LUSC_HRS"="#F7A034","LUSC_LRS"="#8661C1"
    ))+
    theme_classic()+
    theme(panel.grid=element_blank())+
    stat_compare_means(aes(group=DSS_group),label = "p.signif", label.x = 1.5,comparisons = new_comparisons)+
    guides(fill="none")
  print(a)
  dev.off()
}

LUAD_DFI_Immunoinhibitors<-merge(Immunoinhibitors_exp,LUAD_group[,c("PATIENT_ID","DFI_group")])
LUSC_DFI_Immunoinhibitors<-merge(Immunoinhibitors_exp,LUSC_group[,c("PATIENT_ID","DFI_group")])
NSCLC_DFI_Immunoinhibitors<-merge(Immunoinhibitors_exp,NSCLC_group[,c("PATIENT_ID","DFI_group")])

Immunoinhibitors_scores2<-rbind(LUAD_DFI_Immunoinhibitors,LUSC_DFI_Immunoinhibitors,NSCLC_DFI_Immunoinhibitors)

library(reshape2)
Immunoinhibitors_DFI<-melt(Immunoinhibitors_scores2,id.vars=c("DFI_group"),
                           measure.vars=c(colnames(Immunoinhibitors_scores2)[2:(ncol(Immunoinhibitors_scores2)-1)]))

library(ggplot2)
library(ggpubr)
#new_comparisons <- list(c("LUAD_HRS","LUAD_LRS"),c("LUSC_HRS","LUSC_LRS"),c("NSCLC_HRS","NSCLC_LRS"))
for (i in 2:(dim(Immunoinhibitors_scores2)[2]-1))  { 
  print(paste0("run to cells:",colnames(Immunoinhibitors_scores2)[i]," i=",i))
  pdf(file = paste0("Immunoinhibitors/",colnames(Immunoinhibitors_scores2)[i],"_Immunoinhibitors_output.pdf"),onefile=F)
  a<-ggplot(Immunoinhibitors_DFI[which(Immunoinhibitors_DFI$variable == colnames(Immunoinhibitors_scores2)[i]),],aes(x=DFI_group,y=value,fill = factor(DFI_group)))+
    geom_violin(aes(fill=DFI_group),position=position_dodge(0.8))+
    geom_boxplot(width = 0.2,notch = F,position=position_dodge(0.8)) +
    labs( y= colnames(Immunoinhibitors_scores2)[i],x='')+
    scale_fill_manual(values = c("NSCLC_HRS"="#F77300","NSCLC_LRS"="#00AF86",
                                 "LUAD_HRS"="#FF615E","LUAD_LRS"="#00B8CF",
                                 "LUSC_HRS"="#F7A034","LUSC_LRS"="#8661C1"
    ))+
    theme_classic()+
    theme(panel.grid=element_blank())+
    stat_compare_means(aes(group=DFI_group),label = "p.signif", label.x = 1.5,comparisons = new_comparisons)+
    guides(fill="none")
  print(a)
  dev.off()
}


#Immunostimulators----
LUAD_DSS_Immunostimulators<-merge(Immunostimulators_exp,LUAD_group[,c("PATIENT_ID","DSS_group")])
LUSC_DSS_Immunostimulators<-merge(Immunostimulators_exp,LUSC_group[,c("PATIENT_ID","DSS_group")])
NSCLC_DSS_Immunostimulators<-merge(Immunostimulators_exp,NSCLC_group[,c("PATIENT_ID","DSS_group")])

Immunostimulators_scores1<-rbind(LUAD_DSS_Immunostimulators,LUSC_DSS_Immunostimulators,NSCLC_DSS_Immunostimulators)

library(reshape2)
Immunostimulators_DSS<-melt(Immunostimulators_scores1,id.vars=c("DSS_group"),
                            measure.vars=c(colnames(Immunostimulators_scores1)[2:(ncol(Immunostimulators_scores1)-1)]))

library(ggplot2)
library(ggpubr)
#new_comparisons <- list(c("LUAD_HRS","LUAD_LRS"),c("LUSC_HRS","LUSC_LRS"),c("NSCLC_HRS","NSCLC_LRS"))
for (i in 2:(dim(Immunostimulators_scores1)[2]-1))  { 
  print(paste0("run to cells:",colnames(Immunostimulators_scores1)[i]," i=",i))
  pdf(file = paste0("Immunostimulators/",colnames(Immunostimulators_scores1)[i],"_Immunostimulators_output.pdf"),onefile=F)
  a<-ggplot(Immunostimulators_DSS[which(Immunostimulators_DSS$variable == colnames(Immunostimulators_scores1)[i]),],aes(x=DSS_group,y=value,fill = factor(DSS_group)))+
    geom_violin(aes(fill=DSS_group),position=position_dodge(0.8))+
    geom_boxplot(width = 0.2,notch = F,position=position_dodge(0.8)) +
    labs( y= colnames(Immunostimulators_scores1)[i],x='')+
    scale_fill_manual(values = c("NSCLC_HRS"="#F77300","NSCLC_LRS"="#00AF86",
                                 "LUAD_HRS"="#FF615E","LUAD_LRS"="#00B8CF",
                                 "LUSC_HRS"="#F7A034","LUSC_LRS"="#8661C1"
    ))+
    theme_classic()+
    theme(panel.grid=element_blank())+
    stat_compare_means(aes(group=DSS_group),label = "p.signif", label.x = 1.5,comparisons = new_comparisons)+
    guides(fill="none")
  print(a)
  dev.off()
}

LUAD_DFI_Immunostimulators<-merge(Immunostimulators_exp,LUAD_group[,c("PATIENT_ID","DFI_group")])
LUSC_DFI_Immunostimulators<-merge(Immunostimulators_exp,LUSC_group[,c("PATIENT_ID","DFI_group")])
NSCLC_DFI_Immunostimulators<-merge(Immunostimulators_exp,NSCLC_group[,c("PATIENT_ID","DFI_group")])

Immunostimulators_scores2<-rbind(LUAD_DFI_Immunostimulators,LUSC_DFI_Immunostimulators,NSCLC_DFI_Immunostimulators)

library(reshape2)
Immunostimulators_DFI<-melt(Immunostimulators_scores2,id.vars=c("DFI_group"),
                            measure.vars=c(colnames(Immunostimulators_scores2)[2:(ncol(Immunostimulators_scores2)-1)]))

library(ggplot2)
library(ggpubr)
#new_comparisons <- list(c("LUAD_HRS","LUAD_LRS"),c("LUSC_HRS","LUSC_LRS"),c("NSCLC_HRS","NSCLC_LRS"))
for (i in 2:(dim(Immunostimulators_scores2)[2]-1))  { 
  print(paste0("run to cells:",colnames(Immunostimulators_scores2)[i]," i=",i))
  pdf(file = paste0("Immunostimulators/",colnames(Immunostimulators_scores2)[i],"_Immunostimulators_output.pdf"),onefile=F)
  a<-ggplot(Immunostimulators_DFI[which(Immunostimulators_DFI$variable == colnames(Immunostimulators_scores2)[i]),],aes(x=DFI_group,y=value,fill = factor(DFI_group)))+
    geom_violin(aes(fill=DFI_group),position=position_dodge(0.8))+
    geom_boxplot(width = 0.2,notch = F,position=position_dodge(0.8)) +
    labs( y= colnames(Immunostimulators_scores2)[i],x='')+
    scale_fill_manual(values = c("NSCLC_HRS"="#F77300","NSCLC_LRS"="#00AF86",
                                 "LUAD_HRS"="#FF615E","LUAD_LRS"="#00B8CF",
                                 "LUSC_HRS"="#F7A034","LUSC_LRS"="#8661C1"
    ))+
    theme_classic()+
    theme(panel.grid=element_blank())+
    stat_compare_means(aes(group=DFI_group),label = "p.signif", label.x = 1.5,comparisons = new_comparisons)+
    guides(fill="none")
  print(a)
  dev.off()
}


#HLA----
LUAD_DSS_HLA<-merge(HLA_exp,LUAD_group[,c("PATIENT_ID","DSS_group")])
LUSC_DSS_HLA<-merge(HLA_exp,LUSC_group[,c("PATIENT_ID","DSS_group")])
NSCLC_DSS_HLA<-merge(HLA_exp,NSCLC_group[,c("PATIENT_ID","DSS_group")])

HLA_scores1<-rbind(LUAD_DSS_HLA,LUSC_DSS_HLA,NSCLC_DSS_HLA)

library(reshape2)
HLA_DSS<-melt(HLA_scores1,id.vars=c("DSS_group"),
              measure.vars=c(colnames(HLA_scores1)[2:(ncol(HLA_scores1)-1)]))

library(ggplot2)
library(ggpubr)
#new_comparisons <- list(c("LUAD_HRS","LUAD_LRS"),c("LUSC_HRS","LUSC_LRS"),c("NSCLC_HRS","NSCLC_LRS"))
for (i in 2:(dim(HLA_scores1)[2]-1))  { 
  print(paste0("run to cells:",colnames(HLA_scores1)[i]," i=",i))
  pdf(file = paste0("HLA/",colnames(HLA_scores1)[i],"_HLA_output.pdf"),onefile=F)
  a<-ggplot(HLA_DSS[which(HLA_DSS$variable == colnames(HLA_scores1)[i]),],aes(x=DSS_group,y=value,fill = factor(DSS_group)))+
    geom_violin(aes(fill=DSS_group),position=position_dodge(0.8))+
    geom_boxplot(width = 0.2,notch = F,position=position_dodge(0.8)) +
    labs( y= colnames(HLA_scores1)[i],x='')+
    scale_fill_manual(values = c("NSCLC_HRS"="#F77300","NSCLC_LRS"="#00AF86",
                                 "LUAD_HRS"="#FF615E","LUAD_LRS"="#00B8CF",
                                 "LUSC_HRS"="#F7A034","LUSC_LRS"="#8661C1"
    ))+
    theme_classic()+
    theme(panel.grid=element_blank())+
    stat_compare_means(aes(group=DSS_group),label = "p.signif", label.x = 1.5,comparisons = new_comparisons)+
    guides(fill="none")
  print(a)
  dev.off()
}

LUAD_DFI_HLA<-merge(HLA_exp,LUAD_group[,c("PATIENT_ID","DFI_group")])
LUSC_DFI_HLA<-merge(HLA_exp,LUSC_group[,c("PATIENT_ID","DFI_group")])
NSCLC_DFI_HLA<-merge(HLA_exp,NSCLC_group[,c("PATIENT_ID","DFI_group")])

HLA_scores2<-rbind(LUAD_DFI_HLA,LUSC_DFI_HLA,NSCLC_DFI_HLA)

library(reshape2)
HLA_DFI<-melt(HLA_scores2,id.vars=c("DFI_group"),
              measure.vars=c(colnames(HLA_scores2)[2:(ncol(HLA_scores2)-1)]))

library(ggplot2)
library(ggpubr)
#new_comparisons <- list(c("LUAD_HRS","LUAD_LRS"),c("LUSC_HRS","LUSC_LRS"),c("NSCLC_HRS","NSCLC_LRS"))
for (i in 2:(dim(HLA_scores2)[2]-1))  { 
  print(paste0("run to cells:",colnames(HLA_scores2)[i]," i=",i))
  pdf(file = paste0("HLA/",colnames(HLA_scores2)[i],"_HLA_output.pdf"),onefile=F)
  a<-ggplot(HLA_DFI[which(HLA_DFI$variable == colnames(HLA_scores2)[i]),],aes(x=DFI_group,y=value,fill = factor(DFI_group)))+
    geom_violin(aes(fill=DFI_group),position=position_dodge(0.8))+
    geom_boxplot(width = 0.2,notch = F,position=position_dodge(0.8)) +
    labs( y= colnames(HLA_scores2)[i],x='')+
    scale_fill_manual(values = c("NSCLC_HRS"="#F77300","NSCLC_LRS"="#00AF86",
                                 "LUAD_HRS"="#FF615E","LUAD_LRS"="#00B8CF",
                                 "LUSC_HRS"="#F7A034","LUSC_LRS"="#8661C1"
    ))+
    theme_classic()+
    theme(panel.grid=element_blank())+
    stat_compare_means(aes(group=DFI_group),label = "p.signif", label.x = 1.5,comparisons = new_comparisons)+
    guides(fill="none")
  print(a)
  dev.off()
}