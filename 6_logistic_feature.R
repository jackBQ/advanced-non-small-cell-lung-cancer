setwd("D:/工作目录/TMB/24_npj_reviewer")
library(maftools)
library(plyr)
load("20240125.RData")

###
mut_filter <- read.table("somatic_mutation_nonsys.maf", sep = "\t", quote = '', header = T, row.names = NULL, check.names = F)
var_v1 <- read.xlsx("new_var_filter.xlsx", sheet = 1)
mut_filter <- mut_filter[mut_filter$Tumor_Sample_Barcode %in% var_v1[var_v1$RECIST %in% c("PD","SD","PR"),]$Sample_ID,]

mut_maf <- read.maf(mut_filter)
# oncoplot(mut_maf,top = 10000, removeNonMutated = F, writeMatrix = T)
allgene_mut <- read.table("allgene_mut_onco_matrix.txt", sep = '\t', quote = '', check.names = F, header = T, row.names = 1)
allgene_mut <- ifelse(allgene_mut == "",0,1)
genemut_mat <- as.data.frame(t(as.data.frame(allgene_mut)))
var_v1 <- var_v1[match(rownames(genemut_mat), var_v1$Sample_ID),]
sample_info <- sample_info[match(rownames(genemut_mat), sample_info$Sample_ID),]
identical(rownames(genemut_mat), sample_info$Sample_ID)
identical(rownames(genemut_mat), var_v1$Sample_ID)

genemut_mat$RECIST <- var_v1$RECIST
genemut_mat$TMB <- as.numeric(sample_info$Report_TMB)
genemut_mat$PDL1 <- var_v1$new_PDL1
genemut_mat$PDL1 <- ifelse(genemut_mat$PDL1 == "<1%",0,
                           ifelse(genemut_mat$PDL1 == "1%-49%",1,2))
genemut_mat$RECIST <- ifelse(genemut_mat$RECIST == "PR", 0,1)
data <- genemut_mat
colnames(data)[216] <- c("HLA_B")
colnames(data)[414] <- c("H3_4")
colnames(data)[454] <- c("NKX2_1")

###logistic----
Uni_glm_model=function(x){ #写循环函数
  print(x)
  FML=as.formula(paste0("RECIST==1~",x)) #构筑分析，注意这里的Group为你准备数据中第二列的列名，不然运行不成功的。
  glm1<- glm(FML,family = binomial,data = data) #单因素分析
  glm2=summary(glm1) #处理分析结果
  OR=round(exp(coef(glm1)),2) #提取风险率
  SE=glm2$coefficients[,2]
  CI5=round(exp(coef(glm1)-1.96*SE),2) #计算风险率范围
  CI95=round(exp(coef(glm1)+1.96*SE),2)
  CI=paste0(CI5,"-",CI95)
  P=signif(glm2$coefficients[,4],4) #提取P值，保留4位小数
  Uni_glm_model <- data.frame("characteristics"=x,
                              "OR"=OR,
                              "CI.UP"=CI95,
                              "CI.LOW"=CI5,
                              "p"=P)[-1,]
  return(Uni_glm_model)
}

variable.names=colnames(data)[c(1:518,520,521)] #把要分析的基因放进去，注意这里的“2:5”表示你要分析基因的位置，从“Group”那一列为第一列。

Uni_glm=lapply(variable.names,Uni_glm_model) #应用函数
Uni_glm=ldply(Uni_glm,data.frame) #整理成数据框
Uni_glm$adj.p <- round(p.adjust(Uni_glm$p, method = "BH"),3)
write.xlsx(Uni_glm,"Uni_glm_allmut_gene.xlsx") #保存结果

###mul
model <- glm(RECIST ~ ARID1A+LZTR1+SERPINB3+PDL1, 
             data = data, family = binomial)
summary(model)

###
sample_info <- sample_info[match(signature_data$Sample_ID, sample_info$Sample_ID),]
identical(signature_data$Sample_ID, sample_info$Sample_ID)
identical(signature_data$Sample_ID, var_v1$Sample_ID)
signature_data$TMB <- sample_info$Report_TMB
signature_data$PDL1 <- var_v1$new_PDL1
signature_data$RECIST <- var_v1$RECIST
signature_data$Response <- ifelse(signature_data$RECIST == "PR", "Response", ifelse(signature_data$RECIST %in% c("PD","SD"), "nonResponse", NA))
table(signature_data$Response)

data <- signature_data[!is.na(signature_data$Response),c("Sample_ID", "TMB","PDL1", "Response")]
data$TMB <- as.numeric(data$TMB)
data$PDL1 <- ifelse(data$PDL1 == "<1%",0, ifelse(data$PDL1 == "1%-49%",1,2))
data$Response <- ifelse(data$Response == "Response", 0,1)
# data$PDL1 <- factor(data$PDL1, levels = c(0,1,2))

###
model <- glm(Response ~ TMB, data = data, family = binomial)
summary(model)$coef

model <- glm(Response ~ PDL1, data = data, family = binomial)
summary(model)$coef

data %>% 
  mutate(prob = ifelse(Response == 1, 1, 0)) %>%
  ggplot(aes(TMB, prob)) +
  geom_point(alpha = 0.2) +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  labs(title = "Logistic Regression Model", 
    x = "TMB",
    y = "Probability of being diabete-pos")

model <- glm( Response ~ TMB+PDL1, data = data, family = binomial)
summary(model)$coef
