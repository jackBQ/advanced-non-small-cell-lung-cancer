DDR_pathway<-c("MMR","BER","NER","HRR","NHEJ","CPF","FA","TLS")
DDR_pathway<-paste0(DDR_pathway,"_pathway",sep = "")
pathway_onco<-read.table("var_lasso_group.txt",header=T,
                         quote="",sep="\t",check.name=F)
pathway_onco<-pathway_onco[,which(colnames(pathway_onco) %in% c(DDR_pathway,"Sample_ID","OS_STATUS","OS_DAYS",
                                                                "PFS_STATUS","PFS_DAYS","Response","RS"))]
rownames(pathway_onco)<-pathway_onco$Sample_ID
pathway_onco<-pathway_onco[,-1]

#DDR 与 Response----
res_onco <- as.data.frame(t(pathway_onco[which(pathway_onco$Response == "Response"),1:7]))
nres_onco <- as.data.frame(t(pathway_onco[which(pathway_onco$Response == "nonResponse"),1:7]))

for (i in (1:nrow(res_onco))){
  res_onco[i,"res_mut_num"] <- length(which(res_onco[i,c(1:ncol(res_onco))] == "1") )
  res_onco[i,"res_non_mut_num"] <- length(which(res_onco[i,c(1:ncol(res_onco))] == "0") )
}
for (i in (1:nrow(nres_onco))){
  nres_onco[i,"nres_mut_num"] <- length(which(nres_onco[i,c(1:ncol(nres_onco))] == "1") )
  nres_onco[i,"nres_non_mut_num"] <- length(which(nres_onco[i,c(1:ncol(nres_onco))] == "0") )
}

res_onco$gene <- rownames(res_onco)
nres_onco$gene <- rownames(nres_onco)
res_diff <- merge(res_onco[,c("res_mut_num","res_non_mut_num","gene")], nres_onco[,c("nres_mut_num","nres_non_mut_num","gene")], by.x = "gene", by.y = "gene")

out <- data.frame()
for (i in 1:nrow(res_diff )){
  t <- fisher.test(matrix(as.vector(t(res_diff[i, 2:5])), ncol=2))
  
  d <- res_diff[i, ]
  d$p.value <- t$p.value
  d$OR <- t$estimate[[1]]
  d$OR.lower95 <- t$conf.int[1]
  d$OR.upper95 <- t$conf.int[2]
  out <- rbind(out, d)
}
Response_fisher <- out
write.table(Response_fisher,"Fisher_Response_DDR_result.txt", sep = '\t', quote = F, row.names = F)

#####barplot
library(tidyr)
library(ggplot2)
pathway_response_barpolt <- gather(Response_fisher[,c(1,2,4)], Response, MT, -gene)
pathway_response_barpolt[pathway_response_barpolt$Response == "nres_mut_num",]$Response <- c("nonResponse")
pathway_response_barpolt[pathway_response_barpolt$Response == "res_mut_num",]$Response <- c("Response")

pdf("Pathway_Response_DDR.pdf",width = 8, height = 6, onefile = F)
ggplot(data = pathway_response_barpolt,aes(x=gene,y=MT,fill=Response))+
  geom_bar(stat = "identity",position = "stack")+
  theme_bw()+
  theme(panel.grid = element_blank(),legend.position = "top")+
  scale_fill_manual(values = c("Response" = "#DAA520", "nonResponse" = "#9370DB"))+
  annotate(geom = "text",
           cex=4,
           x=1:7, y=140, # 根据自己的数据调节p value的位置
           label=paste0("P ", 
                        ifelse(Response_fisher$p.value<0.001, 
                               "< 0.001", 
                               paste0("= ",round(Response_fisher$p.value,3)))), # 添加P值
           color="black")+
  geom_text(aes(label = MT), position = position_stack(vjust = 0.5), color = "black", size = 4)+
  labs(x='Pathways',y='MT count')
dev.off()



#cox分析-----
#*PFS-----
library(survival)
gene_Uni_cox_model<- function(x,data){
  FML <- as.formula(paste0 ("Surv(time=PFS_DAYS,event=PFS_STATUS)~",x))
  cox<- coxph(FML,data)
  cox1<-summary(cox)
  HR <- round(cox1$coefficients[,2],2)
  PValue <- round(cox1$coefficients[,5],3)
  CI5 <-round(cox1$conf.int[,3],2)
  CI95 <-round(cox1$conf.int[,4],2)
  Uni_cox_model<- data.frame('Characteristics' = x,
                             'HR' = HR,
                             'CI5' = CI5,
                             'CI95' = CI95,
                             'p' = PValue)
  return(Uni_cox_model)}  
#**3.将想要进行的单因素回归变量输入模型
#**4.输出结果
gene_Uni_cox <- lapply(colnames(pathway_onco)[1:7],data=pathway_onco, gene_Uni_cox_model)
library(plyr)
gene_Uni_cox<- ldply(gene_Uni_cox,data.frame)
#**5.优化表格，这里举例HR+95% CI+P 风格
gene_Uni_cox$CI<-paste(gene_Uni_cox$CI5,'-',gene_Uni_cox$CI95)
gene_Uni_cox<-cbind(gene_Uni_cox,res_diff[,c("res_mut_num","nres_mut_num")])
write.table(gene_Uni_cox,"DDR_PFS_Uni_cox_output.txt",sep="\t",row.names = F,col.names = T,quote=F)

library(forestplot)
gene_Uni_cox1<-gene_Uni_cox
#**2. 给参考变量插入空行
ins <- function(x) {c(x, rep(NA, ncol(gene_Uni_cox1)-1))}
for(i in 5:6) {gene_Uni_cox1[, i] = as.character(gene_Uni_cox1[, i])}
gene_Uni_cox1<-rbind(c("Genes", NA, NA, NA,"p","HR(95%CI)","Response","nonResponse"),
                     gene_Uni_cox1)
for(i in 2:4) {gene_Uni_cox1[, i] = as.numeric(gene_Uni_cox1[, i])}
pdf("DDR_PFS_uni_forestplot.pdf",onefile = FALSE,width = 8,height = 5)
forestplot(gene_Uni_cox1[,c(1,7,8,6,5)],  #1,6,5列显示为变量 HR(CI) p数值形式
           mean=gene_Uni_cox1[,2],   #表格第3列为HR，要变成森林图的小方块
           lower=gene_Uni_cox1[,3],  #表格第7列为5%CI，
           upper=gene_Uni_cox1[,4],  #表格第8列为95%CI，它俩要化作线段，穿过方块
           zero=1,            #零线或参考线为HR=1即x轴的垂直线
           boxsize=0.2,       #设置小黑块的大小
           graph.pos=4,#"right",#森林图放在最右侧
           hrzl_lines=list("1" = gpar(lty=1,lwd=2),
                           "2" = gpar(lty=2),
                           "9"= gpar(lwd=2,lty=1)),
           graphwidth = unit(.25,"npc"),
           #xlab="HR (95% CI)",
           #xticks=c(0,1,2,4,6,8,10,) ,
           #----------------#线条粗细（x轴、置信区间）
           lwd.zero=1,
           lwd.ci=1.5,
           lwd.xaxis=1, 
           lty.ci=1,
           ci.vertices =T,
           ci.vertices.height=0.2, 
           clip=c(0.1,5),
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

#DDR 与 model group----
LRS_onco <- as.data.frame(t(pathway_onco[which(pathway_onco$RS == "LRS"),1:7]))
HRS_onco <- as.data.frame(t(pathway_onco[which(pathway_onco$RS == "HRS"),1:7]))

for (i in (1:nrow(LRS_onco))){
  LRS_onco[i,"LRS_mut_num"] <- length(which(LRS_onco[i,c(1:ncol(LRS_onco))] == "1") )
  LRS_onco[i,"LRS_non_mut_num"] <- length(which(LRS_onco[i,c(1:ncol(LRS_onco))] == "0") )
}
for (i in (1:nrow(HRS_onco))){
  HRS_onco[i,"HRS_mut_num"] <- length(which(HRS_onco[i,c(1:ncol(HRS_onco))] == "1") )
  HRS_onco[i,"HRS_non_mut_num"] <- length(which(HRS_onco[i,c(1:ncol(HRS_onco))] == "0") )
}

LRS_onco$gene <- rownames(LRS_onco)
HRS_onco$gene <- rownames(HRS_onco)
RS_diff <- merge(LRS_onco[,c("LRS_mut_num","LRS_non_mut_num","gene")], HRS_onco[,c("HRS_mut_num","HRS_non_mut_num","gene")], by.x = "gene", by.y = "gene")

RS_out <- data.frame()
for (i in 1:nrow(RS_diff )){
  RS_t <- fisher.test(matrix(as.vector(t(RS_diff[i, 2:5])), ncol=2))
  
  RS_d <- RS_diff[i, ]
  RS_d$p.value <- RS_t$p.value
  RS_d$OR <- RS_t$estimate[[1]]
  RS_d$OR.lower95 <- RS_t$conf.int[1]
  RS_d$OR.upper95 <- RS_t$conf.int[2]
  RS_out <- rbind(RS_out, RS_d)
}
RS_fisher <- RS_out
write.table(RS_fisher,"Fisher_RS_DDR_result.txt", sep = '\t', quote = F, row.names = F)

#####barplot
library(tidyr)
library(ggplot2)
pathway_RS_barpolt <- gather(RS_fisher[,c(1,2,4)], RS, MT, -gene)
pathway_RS_barpolt[pathway_RS_barpolt$RS == "HRS_mut_num",]$RS <- c("HRS")
pathway_RS_barpolt[pathway_RS_barpolt$RS == "LRS_mut_num",]$RS <- c("LRS")
  
pdf("Pathway_RS_DDR.pdf",width = 8, height = 6, onefile = F)
ggplot(data = pathway_RS_barpolt,aes(x=gene,y=MT,fill=RS))+
  geom_bar(stat = "identity",position = "stack")+
  theme_bw()+
  theme(panel.grid = element_blank(),legend.position = "top")+
  scale_fill_manual(values = c("HRS" = "#CC0000", "LRS" = "#2f5688"))+
  annotate(geom = "text",
           cex=4,
           x=1:7, y=140, # 根据自己的数据调节p value的位置
           label=paste0("P ", 
                        ifelse(RS_fisher$p.value<0.001, 
                               "< 0.001", 
                               paste0("= ",round(RS_fisher$p.value,3)))), # 添加P值
           color="black")+
  geom_text(aes(label = MT), position = position_stack(vjust = 0.5), color = "black", size = 4)+
  labs(x='Pathways',y='MT count')
dev.off()
