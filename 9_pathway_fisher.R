pathway_onco<-read.table("var_lasso_group.txt",header=T,quote="",sep="\t",check.name=F)
pathway<-colnames(pathway_onco)[210:230]
pathway_onco<-pathway_onco[,which(colnames(pathway_onco) %in% c(pathway,"Sample_ID","OS_STATUS","OS_DAYS",
                                                                "PFS_STATUS","PFS_DAYS","Response","RS"))]
rownames(pathway_onco)<-pathway_onco$Sample_ID
pathway_onco<-pathway_onco[,-1]

#fisher.test----
res_onco <- as.data.frame(t(pathway_onco[which(pathway_onco$Response == "Response"),1:21]))
nres_onco <- as.data.frame(t(pathway_onco[which(pathway_onco$Response == "nonResponse"),1:21]))

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
write.table(Response_fisher,"fisher_Response_allpathway_result.txt", sep = '\t', quote = F, row.names = F)

#pathway_barplot_view----
library(tidyr)
library(ggplot2)
pathway_response_barpolt <- gather(Response_fisher[,c(1,2,4)], Response, MT, -gene)
pathway_response_barpolt[pathway_response_barpolt$Response == "nres_mut_num",]$Response <- c("nonResponse")
pathway_response_barpolt[pathway_response_barpolt$Response == "res_mut_num",]$Response <- c("Response")

pdf("Pathway_Response_allPathway.pdf",width = 15, height = 6, onefile = F)
ggplot(data = pathway_response_barpolt,aes(x=gene,y=MT,fill=Response))+
  geom_bar(stat = "identity",position = "stack")+
  theme_bw()+
  theme(panel.grid = element_blank(),legend.position = "top")+
  scale_fill_manual(values = c("Response" = "#DAA520", "nonResponse" = "#9370DB"))+
  annotate(geom = "text",
           cex=4,
           x=1:21, y=140, 
           label=paste0("P ", 
                        ifelse(Response_fisher$p.value<0.001, 
                               "< 0.001", 
                               paste0("= ",round(Response_fisher$p.value,3)))), # 添加P值
           color="black")+
  geom_text(aes(label = MT), position = position_stack(vjust = 0.5), color = "black", size = 4)+
  labs(x='Pathways',y='MT count')
dev.off()

###data_input
pathway_gene <- read.table("new_pathways.txt",header=T,quote="",sep="\t",check.name=F)
mut_maf <- read.maf(maf = "somatic_mutation_nonsys.maf")
summary_pathway <- pathways(mut_maf,pathways = pathway_gene)
write.xlsx(summary_pathway,"allPathway_summary.xlsx")

#logistic----
pathway_onco$Response <- ifelse(pathway_onco$Response == "Response",0,1)
Uni_glm_model=function(x){ 
  print(x)
  FML=as.formula(paste0("Response==1~",x)) 
  glm1<- glm(FML,family = binomial,data = data) 
  glm2=summary(glm1) 
  OR=round(exp(coef(glm1)),2) 
  SE=glm2$coefficients[,2]
  CI5=round(exp(coef(glm1)-1.96*SE),2) 
  CI95=round(exp(coef(glm1)+1.96*SE),2)
  CI=paste0(CI5,"-",CI95)
  P=signif(glm2$coefficients[,4],4) 
  Uni_glm_model <- data.frame("characteristics"=x,
                              "OR"=OR,
                              "CI.UP"=CI95,
                              "CI.LOW"=CI5,
                              "p"=P)[-1,]
  return(Uni_glm_model)
}

variable.names=colnames(pathway_onco)[c(1:21)]

Uni_glm=lapply(variable.names,Uni_glm_model) 
Uni_glm=ldply(Uni_glm,data.frame) 
Uni_glm$adj.p <- round(p.adjust(Uni_glm$p, method = "BH"),3)
write.xlsx(Uni_glm,"glm_allmut_pathway.xlsx")
