library(openxlsx)
library(lars) 
library(glmnet) 
library(ggplot2)
library(survival)
library(survminer)
library(timeROC)
library(tidyr)
library(plyr)
library(ggsci)

#####
var_v1 <- read.xlsx("new_var_filter.xlsx", sheet = 1)
var_v2 <- read.table("newresponse_var_filter.txt", sep = '\t', quote = '', check.names = F, header = T, row.names = NULL)

#Univariate_cox----
source("Cox.R")

cox_clin<-var_v1[,c("Sample_ID","PFS_STATUS","PFS_DAYS")] 
dc<-var_v1[,c(2:230)]

rownames(dc)<-var_v1[,1]
dc<-t(dc)

colnames(cox_clin)<-c("sample","Events","Survival")
cox_uni<-as.data.frame(getUniOrMultiCOXAnalysis(dc,cox_clin,method="Univariate"))
# write.table(cox_uni, "pfs_cox_univariate_allvars.txt",sep="\t",quote=F,row.names = F)

#Lasso----
sig_feature <- cox_uni[cox_uni$`Pr(>|Z|)[pval]` < 0.05,]
sig_feature <- sig_feature[!is.na(sig_feature$`Pr(>|Z|)[pval]`),]
rownames(var_v1) <- var_v1$Sample_ID
sig_mut <- var_v1[,colnames(var_v1) %in% sig_feature$geneName]
table(rownames(sig_mut) == rownames(var_v1))
sig_mut$PFS_STATUS <- var_v1$PFS_STATUS
sig_mut$PFS_DAYS <- var_v1$PFS_DAYS
sig_mut <- sig_mut[!is.na(sig_mut$PFS_STATUS),]

phe <- sig_mut[,c("PFS_STATUS","PFS_DAYS")]
rownames(phe) <- rownames(sig_mut)
mutSet <- sig_mut[,1:(ncol(sig_mut)-2)]
mutSet <- t(mutSet)
x=t(mutSet)
y=phe$PFS_STATUS

set.seed(4)
model_lasso <- glmnet(x, y, family="binomial", nlambda=50, alpha=1)#拉手回归模型
print(model_lasso)

# pdf(file="uni_lambda_pfs.pdf",width=12,height = 10)
plot(model_lasso, xvar="lambda", label=TRUE)
# dev.off()

# pdf(file = "uni_cvfit_pfs.pdf",width = 12,height = 10)
cv_fit <- cv.glmnet(x=as(x,"dgCMatrix"), y=y, alpha = 1, nfolds = 10)
plot(cv_fit)
# dev.off()

x <- coef(model_lasso) 
tmp <- as.data.frame(as.matrix(x)) 
tmp$coef <- row.names(tmp) 
tmp <- reshape::melt(tmp, id = "coef") 
tmp$variable <- as.numeric(gsub("s", "", tmp$variable)) 
tmp$coef <- gsub('_','-',tmp$coef) 
tmp$lambda <- model_lasso$lambda[tmp$variable+1] 
# extract the lambda values 
tmp$norm <- apply(abs(x[-1,]), 2, sum)[tmp$variable+1] 
# compute L1 norm

pdf("lambda.pdf",15,10)
ggplot(tmp,aes(log(lambda),value,color = coef)) + 
  geom_vline(xintercept = log(cv_fit$lambda.min),
             size=0.8,color='grey60',
             alpha=0.8,linetype=2)+
  geom_line(size=1) + 
  xlab("Lambda (log scale)") + 
  ylab('Coefficients')+ 
  theme_bw(base_rect_size = 2)+ 
  scale_color_manual(values= c(pal_npg()(10),
                               pal_d3()(10),
                               pal_jco()(3),
                               pal_lancet()(10),
                               pal_aaas()(10),
                               pal_simpsons()(10),
                               pal_gsea()(10),
                               pal_jama()(10)))+ 
  scale_x_continuous(expand = c(0.01,0.01))+ 
  scale_y_continuous(expand = c(0.01,0.01))+ 
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size=15,color='black'), 
        axis.text = element_text(size=12,color='black'), 
        legend.title = element_blank(), 
        legend.text = element_text(size=12,color='black'), 
        legend.position = 'right')+ 
  annotate('text',x = -3.3,y=0.26,
           label='Optimal Lambda = 0.0351',
           color='black')+ 
  guides(col=guide_legend(ncol = 2))
dev.off()

####RS_model----
lasso_coef <- as.data.frame(as.matrix(l.coef1))
lasso_coef$gene <- rownames(lasso_coef)
lasso_coef <- lasso_coef[which(lasso_coef$s1 != 0),]
lasso_coef <- lasso_coef[-1,]
# write.table(lasso_coef[,c(2,1)],"lasso_coef.txt", sep = '\t', quote = F, row.names = F)

var_v1$Score <- var_v1[,lasso_coef[1,2]]*lasso_coef[1,1]+
  var_v1[,lasso_coef[2,2]]*lasso_coef[2,1]+
  var_v1[,lasso_coef[3,2]]*lasso_coef[3,1]+
  var_v1[,lasso_coef[4,2]]*lasso_coef[4,1]+
  var_v1[,lasso_coef[5,2]]*lasso_coef[5,1]+
  var_v1[,lasso_coef[6,2]]*lasso_coef[6,1]+
  var_v1[,lasso_coef[7,2]]*lasso_coef[7,1]

var_v1$RS <- ifelse(var_v1$Score > median(var_v1$Score), "HRS", "LRS")
# write.table(var_v1, "var_lasso_group.txt", sep = '\t', quote = F, row.names = F)

tmb_var <- var_v1[var_v1$TMB_group == 1,]

var_v1$rs_group <- ifelse(var_v1$RS == "HRS", 1, 0)
fit2 <- survfit(Surv(OS_DAYS, OS_STATUS) ~ RS, data = var_v1)
p <- ggsurvplot(fit2, data = var_v1,font.title=c(20, "bold"),conf.int=F, pval=TRUE, legend.title ="",ylab=c("Overall Survival"),surv.median.line = "hv",
                palette = c("#B22222", "#000080"),legend.labs = c("HRS","LRS"))
res_cox<-coxph(Surv(OS_DAYS, OS_STATUS) ~ rs_group, data = var_v1)
p$plot = p$plot + 
  ggplot2::annotate("text",x = 20, y = 0.15,
                    label = paste("HR :",round(summary(res_cox)$conf.int[1],2))) + 
  ggplot2::annotate("text",x = 20, y = 0.10,
                    label = paste("(","95%CI:",round(summary(res_cox)$conf.int[3],2),"-",round(summary(res_cox)$conf.int[4],2),")",sep = ""))
surv_median(fit2)
pdf("RS_OS_survival_0507.pdf", width = 5, height = 5,onefile = F)
p
dev.off()

fit2 <- survfit(Surv(PFS_DAYS, PFS_STATUS) ~ RS, data = var_v1)
p <- ggsurvplot(fit2, data = var_v1,font.title=c(20, "bold"),conf.int=F, pval=TRUE, legend.title ="",ylab=c("Progression Free Survival"),surv.median.line = "hv",
                palette = c("#B22222", "#000080"),legend.labs = c("HRS","LRS"))
res_cox<-coxph(Surv(PFS_DAYS, PFS_STATUS) ~ rs_group, data = var_v1)
p$plot = p$plot + 
  ggplot2::annotate("text",x = 20, y = 0.15,
                    label = paste("HR :",round(summary(res_cox)$conf.int[1],2))) + 
  ggplot2::annotate("text",x = 20, y = 0.10,
                    label = paste("(","95%CI:",round(summary(res_cox)$conf.int[3],2),"-",round(summary(res_cox)$conf.int[4],2),")",sep = ""))
surv_median(fit2)
pdf("RS_PFS_survival_0507.pdf", width = 5, height = 5,onefile = F)
p
dev.off()

#####TIMEroc
ROC <- timeROC(T = var_v1$PFS_DAYS, #将结局时间转换为月
               delta = as.numeric(var_v1$PFS_STATUS), #生存结局
               marker = var_v1$Score, #预测变量 
               cause = 1, #阳性结局赋值
               weighting = "marginal", #权重计算方法，marginal是默认值，采用km计算删失分布
               times = c(30*6, 30*8, 30*12, 30*16, 30*24), #预测1、3、5、10年时间
               iid = T) #只有marginal可以计算置信区间 
confint(ROC)$CI_AUC

# pdf("Lasso_RS_ROC.pdf", width = 7, height = 7)
{
  plot(ROC, time=30*6, col = "DarkCyan", add = F, title = F, lwd = 2)
  plot(ROC, time=30*8, col = "ForestGreen", add = T, lwd = 2)
  plot(ROC, time=30*12, col = "FireBrick", add = T, lwd = 2)
  plot(ROC, time=30*16, col = "Chocolate", add = T, lwd = 2)
  plot(ROC, time=30*24, col = "Violet", add = T, lwd = 2)
  legend("bottomright", lty = 1, cex = 1.0,
         col = c("DarkCyan", "ForestGreen", "FireBrick","Chocolate","Violet"),
         legend = c("AUC of 6 months survival:0.706(0.606-0.806)", 
                    "AUC of 8 months survival:0.791(0.708-0.875)",
                    "AUC of 12 months survival:0.709(0.611-0.808)",
                    "AUC of 16 months survival:0.734(0.641-0.826)",
                    "AUC of 24 months survival:0.777(0.710-0.844)"))
}
# dev.off()
