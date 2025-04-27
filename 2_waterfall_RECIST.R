library(openxlsx)
library(maftools)

maf <- read.maf(maf = "somatic_mutation_nonsys.maf")
clin<-sample_info[,c("Sample_ID","WES","DNA_Fusion","RNA_Fusion","Smoking",
                     "Drinking","new_PDL1","Pathological_subtype",
                     "T","M","N","Stage","RECIST",
                     "OS_STATUS","PFS_STATUS")]
clin<-merge(clin,var_data[,c("Sample_ID","Age","Gender",
                             "Mean_HE","TMB_group","HLA_group")])
clin$OS_STATUS<-gsub("Yes","Death",clin$OS_STATUS)
clin$OS_STATUS<-gsub("No","Alive",clin$OS_STATUS)
clin$DNA_Fusion<-gsub("Fusion","DNA_Fusion",clin$DNA_Fusion)
clin$RNA_Fusion<-gsub("Fusion","RNA_Fusion",clin$RNA_Fusion)
clin$Smoking<-ifelse(clin$Smoking == "Yes","Smoking","Non-smoking")
clin$Drinking<-ifelse(clin$Drinking == "Yes","Drinking","Non_drinking")
clin$Age<-ifelse(clin$Age == 1,"<=65",">65")
clin$Gender<-ifelse(clin$Gender == 1,"Male","Female")
colnames(clin)<-gsub("Mean_HE","HED",colnames(clin))
clin$TMB_group<-ifelse(clin$TMB_group == 1,"TMB-H","TMB-L")
clin$HLA_group<-ifelse(clin$HLA_group ==0 ,"homozygous","heterozygous")
clin$RECIST<-gsub("Nx","NX",clin$RECIST)
clin$HED<-ifelse(clin$HED == 0,"HED-L","HED-H")

#*PR/CR-PD/SD样本分子图谱----
#突变图谱------
# To add annotation on right y axis: https://github.com/tidyverse/ggplot2/issues/3171
guide_axis_label_trans <- function(label_trans = identity, ...) {
  axis_guide <- guide_axis(...)
  axis_guide$label_trans <- rlang::as_function(label_trans)
  class(axis_guide) <- c("guide_axis_trans", class(axis_guide))
  axis_guide
}
guide_train.guide_axis_trans <- function(x, ...) {
  trained <- NextMethod()
  trained$key$.label <- x$label_trans(trained$key$.label)
  trained
}

#自定义函数，用于创建图例
create_legends <- function(labels_list, colors_list, titles_list, orientation = "horizontal", size = 15) {
  num_legends <- length(labels_list)
  legend_plots <- list()
  # 计算每个分组中类型数的最大值
  max_labels <- sapply(labels_list, function(x) length(x))
  max_labels <- max(max_labels)
  
  for (i in 1:num_legends) {
    data <- data.frame(Label = labels_list[[i]], Color = colors_list[[i]])
    if (orientation == "horizontal") {
      p <- ggplot(data, aes(x = Label, y = 0, fill = Color)) +
        geom_tile(width=0.9*length(labels_list[[i]])/max_labels,height=0.9*length(labels_list[[i]])/max_labels) +
        coord_equal()+
        scale_fill_identity() +
        scale_y_discrete(position = "right")+
        theme_void() +
        theme(axis.text.y.right = element_text(size = 8),
              legend.position = "none",
              plot.title = element_text(hjust = 0.5)) +
        ggtitle(titles_list[[i]])+
        guides(fill = guide_legend(title = titles_list[[i]], label.position = "bottom", keywidth = unit(size, "points")))
    } else {
      p <- ggplot(data, aes(x = 0, y = Label, fill = Color)) +
        geom_tile(width=0.9*length(labels_list[[i]])/max_labels,height=0.9*length(labels_list[[i]])/max_labels) +
        coord_equal()+
        scale_fill_identity() +
        scale_y_discrete(position = "right")+
        theme_void() +
        theme(axis.text.y.right = element_text(size = 8),
              legend.position = "none",
              plot.title = element_text(hjust = 0.5)) +
        ggtitle(titles_list[[i]])+
        guides(fill = guide_legend(title = titles_list[[i]], label.position = "right", keyheight = unit(size, "points")))
    }
    
    legend_plots[[i]] <- p
  }
  
  return(legend_plots)
}



genes = getGeneSummary(x = maf)[1:30, Hugo_Symbol]  # Top30 genes
vc_cols = c("#3578AD","#FFFF33","#8F4B99","#040000","#FF7E00","#4CA74A","#D6231F","#F780BF","#EDEDED")#
names(vc_cols) = c('Missense_Mutation', 'Frame_Shift_Del','Nonsense_Mutation',
                   'Multi_Hit','Frame_Shift_Ins', 'In_Frame_Ins','Splice_Site',
                   'In_Frame_Del',"Non_mut")
type_levels_2<-c('Missense_Mutation', 'Frame_Shift_Del','Nonsense_Mutation',
                 'Multi_Hit','Frame_Shift_Ins', 'In_Frame_Ins','Splice_Site',
                 'In_Frame_Del',"Non_mut")

col=c("WES"="#8DB8FF",
      "DNA_Fusion"="#79BDA4",
      "RNA_Fusion"="#DDA0DD",
      "CR"="#006DB0" ,"PR"="#9370DB","PD"="#40b8bb","SD"="#F8D568","NX"="#F9CFF9","待评估"="#F0945C",
      "Death"="#FF91AF","Alive"="#9E74A5",
      "Yes"="#827717","No"="#277FB0",
      "Female"="#FFE4E1","Male"="#6CA6CD",
      "<=65"="#E1FFE0",">65"="#6E7B8B",
      "Non-smoking"="#FD9F73","Smoking"="#76E6B4",
      "Non_drinking"="#CBEED8","Drinking"="#3EAF37",
      "LUAD"="#D8B365","LUSC"="#5BB5AC","Other"="#DE526C",
      "T1"="#FAE5B8","T2"="#E4CD87","T3"="#EFBC91","T4"="#E79A90","Tx"="#C69287",
      "M0"="#B8DBB3","M1"="#72B063", "Mx"= "#2E6D4E",
      "N0"="#C6DBEF","N1"="#6BAED6","N2"="#4583B6","N3"="#2171B5","Nx"="#08519C",
      "III"="#257D8B","IV"="#EAA558",
      ">=50%"="#BFDFD2", "1%-49%"="#94CFD6","<1%"="#68BED9",
      "TMB-H"="#DEAA87","TMB-L"="#E9DDAF",
      #"WES-TMB-H"="#93AC93","WES-TMB-L"="#6F8A91",
      "heterozygous"="#F3E8D8","homozygous"="#FFB889",
      "HED-H" ="#BDDAEB","HED-L"="#4189C1")#HED



#**PR maf----------
PR_maf = read.maf(maf = "PR_somatic_mutation.maf")

#step1: data preparation were based on maftools package
PR_om = maftools:::createOncoMatrix(m = PR_maf,g=PR_maf@data$Hugo_Symbol)
PR_numMat = PR_om$numericMatrix  # gene/sample ~ frequency
PR_numMat<-PR_numMat[which(rownames(PR_numMat) %in% genes),]
PR_mat_origin = PR_om$oncoMatrix # gene/sample ~ variant type
PR_mat_origin<-PR_mat_origin[which(rownames(PR_mat_origin) %in% genes),]
#step2  main
library(dplyr)
library(ggplot2)
library(stringr)
PR_totSamps = as.numeric(PR_maf@summary[3, summary])
PR_percent_alt = paste0(round(100 * (apply(PR_numMat, 1, function(x) length(x[x != 0]))/PR_totSamps)), "%")
PR_p_main = PR_mat_origin[rev(rownames(PR_mat_origin)),] %>%
  reshape2::melt() %>%
  dplyr::rename(Gene=Var1, Sample=Var2, Type=value) %>%
  dplyr::mutate(Type2 = ifelse(Type=="","Non_mut",Type)) %>%
  dplyr::mutate(Type2 = factor(Type2, levels = type_levels_2)) %>%
  ggplot(aes(x=Sample, y=Gene, fill=Type2)) +
  geom_tile(colour="white") +
  scale_fill_manual(name = NULL,values = vc_cols,
                    breaks = rev(setdiff(type_levels_2,"Non_mut"))) +#guide=guide_legend(byrow = T,nrow = 1)
  theme_void() +
  theme(legend.position = "bottom",#"none"#修改图例位置
        legend.text = element_text(size = 12),
        legend.key.size = unit(12,"pt")) +#"none"
  theme(axis.text.x = element_blank(),#element_text(angle=90),
        axis.text.y.left = element_text(hjust = 0.95),
        axis.text.y.right = element_text(hjust = 0.05),
        plot.margin = unit(c(0, 0, 3, 0), "lines")) +
  scale_y_discrete(labels = rev(PR_percent_alt))+
  guides(y.sec = guide_axis_label_trans(~str_pad(rev(rownames(PR_mat_origin)),5,side = "right")))

#step3:TMB plot
PR_samp_sum = getSampleSummary(x = PR_maf) %>%
  as.data.frame() %>%
  dplyr::select(!total) %>%
  tibble::column_to_rownames("Tumor_Sample_Barcode")
PR_top_bar_data = t(PR_samp_sum[colnames(PR_numMat),, drop = FALSE])
dim(PR_top_bar_data)
# [1]   7 16

PR_p_top = PR_top_bar_data %>%
  reshape2::melt() %>%
  dplyr::rename(Type=Var1, Sample=Var2, Freq=value) %>%
  dplyr::mutate(Type=factor(Type, levels = type_levels_2)) %>%
  ggplot(aes(x=Sample,y=Freq,fill=Type)) +
  geom_col(position="stack") +
  scale_fill_manual(values = vc_cols) +
  theme_void() +
  theme(legend.position = "none") +
  theme(plot.margin = unit(c(1, 0, 1, 2), "lines"),
        axis.line.y = element_line(color="black",linewidth=0.6),
        axis.ticks.length.y.left = unit(6, "pt"),
        axis.text.y.left = element_text(vjust = 0.5, hjust = 0.5),
        axis.ticks.y.left = element_line(linewidth=0.6),
        axis.title.y.left = element_text(angle = 90, size = 13,vjust = 0.5, hjust = 0.5,
                                         margin = margin(r = -6, unit = "pt"))
  ) +
  scale_y_continuous(expand=c(0,0), breaks = c(0, max(colSums(PR_top_bar_data)))) +
  ylab("Total mutation")

###step4: the clinical plot
library(ComplexHeatmap)
PR_clin<-clin[which(clin$RECIST %in% c("PR","CR")),]

PR_clin<-PR_clin[match(colnames(PR_mat_origin),PR_clin$Sample_ID),]
rownames(PR_clin)<-PR_clin$Sample_ID
PR_clin<-PR_clin[,-1]

PR_clin1<-as.matrix(PR_clin[rownames(PR_clin),]) %>%
  reshape2::melt() %>%
  dplyr::rename(patient=Var1, variable=Var2, Type=value) 

PR_clin1$patient=factor(PR_clin1$patient,levels = unique(PR_clin1$patient))
PR_clin1$variable=factor(PR_clin1$variable,levels = rev(unique(PR_clin1$variable)))

PR_p_middle<-PR_clin1%>%
  ggplot(aes(x=patient,y=variable))+
  geom_tile(aes(fill=Type),color="white",size=1)+ #color和size分别指定方块边线的颜色和粗细
  theme_void() +
  scale_fill_manual(values = col,na.value = "white")+ #指定自定义的颜色
  theme(legend.position = "none",
        axis.text.y.right = element_text(hjust = 0.05),
        axis.text.x.bottom = element_blank(),
        axis.text.y.left = element_blank(), 
        axis.ticks = element_blank() )+#不显示坐标轴刻度
  guides(y.sec = guide_axis_label_trans(~str_pad(rev(unique(PR_clin1$variable)),6,side = "right")))
#step5 merge
library(aplot)
options(aplot_guides = "keep")
PR_merge=PR_p_middle %>%
  insert_bottom(PR_p_main,height = 2) %>%
  insert_top(PR_p_top,height = 0.3)


#**PDSD maf-------
PDSD_maf = read.maf(maf = "PDSD_somatic_mutation.maf")

#step1: data preparation were based on maftools package
PDSD_om = maftools:::createOncoMatrix(m = PDSD_maf,g=PDSD_maf@data$Hugo_Symbol)
PDSD_numMat = PDSD_om$numericMatrix  # gene/sample ~ frequency
PDSD_numMat<-PDSD_numMat[which(rownames(PDSD_numMat) %in% genes),]
PDSD_numMat<-PDSD_numMat[match(rownames(PR_numMat),rownames(PDSD_numMat)),]
rownames(PDSD_numMat)<-rownames(PR_numMat)
PDSD_numMat[is.na(PDSD_numMat)]=0

PDSD_mat_origin = PDSD_om$oncoMatrix # gene/sample ~ variant type
PDSD_mat_origin<-PDSD_mat_origin[which(rownames(PDSD_mat_origin) %in% genes),]
PDSD_mat_origin<-PDSD_mat_origin[match(rownames(PR_mat_origin),rownames(PDSD_mat_origin)),]
rownames(PDSD_mat_origin)<-rownames(PR_mat_origin)
PDSD_mat_origin[is.na(PDSD_mat_origin)] <- ""

#step2  main
PDSD_totSamps = as.numeric(PDSD_maf@summary[3, summary])
PDSD_percent_alt = paste0(round(100 * (apply(PDSD_numMat, 1, function(x) length(x[x != 0]))/PDSD_totSamps)), "%")
PDSD_p_main = PDSD_mat_origin[rev(rownames(PDSD_mat_origin)),] %>%
  reshape2::melt() %>%
  dplyr::rename(Gene=Var1, Sample=Var2, Type=value) %>%
  dplyr::mutate(Type2 = ifelse(Type=="","Non_mut",Type)) %>%
  dplyr::mutate(Type2 = factor(Type2, levels = type_levels_2)) %>%
  ggplot(aes(x=Sample, y=Gene, fill=Type2)) +
  geom_tile(colour="white") +
  scale_fill_manual(name = NULL,values = vc_cols,
                    breaks = rev(setdiff(type_levels_2,"Non_mut"))) +
  theme_void() +
  theme(legend.position = "bottom",#"none"#修改图例位置
        legend.text = element_text(size = 12),
        legend.key.size = unit(12,"pt")) +
  theme(axis.text.x = element_blank(),#element_text(angle=90),
        axis.text.y.left = element_text(hjust = 0.95),
        axis.text.y.right = element_text(hjust = 0.05),
        plot.margin = unit(c(0, 0, 3, 0), "lines")) +
  #scale_y_discrete(labels = rev(PDSD_percent_alt))+
  guides(y.sec = guide_axis_label_trans(~str_pad(rev(PDSD_percent_alt),5,side = "right")))

#step3:TMB plot
PDSD_samp_sum = getSampleSummary(x = PDSD_maf) %>%
  as.data.frame() %>%
  dplyr::select(!total) %>%
  tibble::column_to_rownames("Tumor_Sample_Barcode")
PDSD_top_bar_data = t(PDSD_samp_sum[colnames(PDSD_numMat),, drop = FALSE])
dim(PDSD_top_bar_data)
# [1]   7 16

PDSD_p_top = PDSD_top_bar_data %>%
  reshape2::melt() %>%
  dplyr::rename(Type=Var1, Sample=Var2, Freq=value) %>%
  dplyr::mutate(Type=factor(Type, levels = type_levels_2)) %>%
  ggplot(aes(x=Sample,y=Freq,fill=Type)) +
  geom_col(position="stack") +
  scale_fill_manual(values = vc_cols) +
  theme_void() +
  theme(legend.position = "none") +
  theme(plot.margin = unit(c(1, 0, 1, 2), "lines"),
        axis.line.y = element_line(color="black",linewidth=0.6),
        axis.ticks.length.y.left = unit(6, "pt"),
        axis.text.y.left = element_text(vjust = 0.5, hjust = 0.5),
        axis.ticks.y.left = element_line(linewidth=0.6),
        axis.title.y.left = element_text(angle = 90, size = 13,vjust = 0.5, hjust = 0.5,
                                         margin = margin(r = -6, unit = "pt"))
  ) +
  scale_y_continuous(expand=c(0,0), breaks = c(0, max(colSums(PDSD_top_bar_data)))) +
  ylab("Total mutation")

###step4: the clinical plot
PDSD_clin<-clin[which(clin$RECIST %in% c("PD","SD")),]

PDSD_clin<-PDSD_clin[match(colnames(PDSD_mat_origin),PDSD_clin$Sample_ID),]
rownames(PDSD_clin)<-PDSD_clin$Sample_ID
PDSD_clin<-PDSD_clin[,-1]

PDSD_clin1<-as.matrix(PDSD_clin[rownames(PDSD_clin),]) %>%
  reshape2::melt() %>%
  dplyr::rename(patient=Var1, variable=Var2, Type=value) 

PDSD_clin1$patient=factor(PDSD_clin1$patient,levels = unique(PDSD_clin1$patient))
PDSD_clin1$variable=factor(PDSD_clin1$variable,levels = rev(unique(PDSD_clin1$variable)))

PDSD_p_middle<-PDSD_clin1%>%
  ggplot(aes(x=patient,y=variable))+
  geom_tile(aes(fill=Type),color="white",size=1)+ #color和size分别指定方块边线的颜色和粗细
  theme_void() +
  scale_fill_manual(values = col,na.value = "white")+ #指定自定义的颜色
  theme(legend.position = "none",
        axis.text.y.right = element_text(hjust = 0.05),
        axis.text.x.bottom = element_blank(),
        axis.text.y.left = element_blank(), 
        axis.ticks = element_blank() )+#不显示坐标轴刻度
  guides(y.sec = guide_axis_label_trans(~str_pad(rev(unique(PDSD_clin1$variable)),6,side = "right")))
#step5 merge
options(aplot_guides = "keep")
PDSD_merge=PDSD_p_middle %>%
  insert_bottom(PDSD_p_main,height = 2) %>%
  insert_top(PDSD_p_top,height = 0.3)

# 数据
generate_labels_list <- function(clin,vc_cols) {
  labels_list <- list()
  # 遍历数据框的每一列
  if(names(vc_cols)[1] == "Mutation"){
    labels_list[[1]] <- "Mutation"
    for (i in 1:ncol(clin)) {
      unique_values <- na.omit(unique(clin[[i]]))
      labels_list[[i+1]] <- unique_values
    }
  }else{
    for (i in 1:ncol(clin)) {
      unique_values <- na.omit(unique(clin[[i]]))
      labels_list[[i]] <- unique_values
    }
  }
  return(labels_list)
}
generate_colors_list <- function(clin,col,vc_cols) {
  colors_list <- list()
  # 遍历数据框的每一列
  if(names(vc_cols)[1] == "Mutation"){
    colors_list[[1]] <- vc_cols["Mutation"]
    for (i in 1:ncol(clin)) {
      unique_values <- na.omit(unique(clin[[i]]))
      colors_list[[i+1]] <- col[unique_values]
    }
  }else{
    for (i in 1:ncol(clin)) {
      unique_values <- na.omit(unique(clin[[i]]))
      colors_list[[i]] <- col[unique_values]
    }
  }
  return(colors_list)
}
titles_list<-as.list(c(colnames(PDSD_clin)))
labels_list<-generate_labels_list(PDSD_clin,vc_cols)
colors_list<-generate_colors_list(PDSD_clin,col,vc_cols)

legend_plots<-create_legends(labels_list, colors_list, titles_list,
                             orientation = "vertical", size = 2)


library(patchwork)
library(gridExtra)
# library(showtext)
# font_add("myfont","C:/Windows/Fonts/simhei.ttf")
# font_families()
# showtext_auto()
pdf("RECIST_waterfall_output_1.pdf",onefile = F,width = 16,height = 10.5)
p1<-(plot_list(PR_merge)+plot_list(PDSD_merge))+plot_layout(widths =  c(2,1.4))
p1/grid.arrange(grobs = legend_plots, ncol = 10)+
  plot_layout(heights = c(1,0.2))
dev.off()
