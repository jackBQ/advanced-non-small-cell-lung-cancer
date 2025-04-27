library(maftools)
library(openxlsx)
library(survival)
library(survminer)
library(tidyr)
library(plyr)
library(ComplexHeatmap)
library(circlize)

#####
mut_maf <- read.maf(maf = "somatic_mutation_nonsys.maf")
# oncoplot(maf=mut_maf, borderCol=NULL, top=20, writeMatrix = T, removeNonMutated = F)

top20_onco <- read.table("top20_onco_matrix.txt", sep = '\t', quote = '', check.names = F, header = T, row.names = 1)
top20_onco[,160:162] <- c("")
top20_onco <- ifelse(top20_onco != "",1,0)
top20_matrix <- as.data.frame(t(top20_onco))
top20_matrix$Sample_ID <- rownames(top20_matrix)

clin <- read.xlsx("clin_final.xlsx", sheet = 1)
#####pathway
pathway_gene <- read.table("new_pathways.txt", sep = "\t", quote = '', check.names = F, header = T, row.names = NULL)
# oncoplot(maf=mut_maf, borderCol=NULL, top=10000, writeMatrix = T)

onco_matrix <- read.table("allgenes_onco_matrix.txt", sep = "\t", quote = '', check.names = F, header = T ,row.names = 1)

pathway_mut <- onco_matrix[rownames(onco_matrix) %in% pathway_gene$Gene,]
pathway_mut <- as.data.frame(ifelse(pathway_mut != "",1,0))
pathway_mut$Gene <- rownames(pathway_mut)
pathway_merge <- merge(pathway_mut, pathway_gene, by = "Gene")

pathway_matrix <- data.frame(matrix(ncol = ncol(pathway_merge)-3, nrow = length(table(pathway_merge$Pathway))))
for (i in 2:(ncol(pathway_merge)-2)){
  pathway_matrix[,i-1] <- aggregate(pathway_merge[,i]~Pathway,data=pathway_merge,FUN=sum)[,2]
}
rownames(pathway_matrix) <- aggregate(pathway_merge[,i]~Pathway,data=pathway_merge,FUN=sum)[,1]
colnames(pathway_matrix) <- colnames(pathway_merge)[2:163] #num_sample
pathway_matrix <- as.data.frame(ifelse(pathway_matrix == 0, 0, 1))
pathway_matrix <- as.data.frame(t(pathway_matrix))
colnames(pathway_matrix)[1:21] <- paste0(colnames(pathway_matrix)[1:21],"_pathway")

#####
pathway_matrix$Sample_ID <- rownames(pathway_matrix)
pathway_clin <- merge(pathway_matrix, clin[,c("Sample_ID","OS_STATUS","OS_DAYS","PFS_STATUS","PFS_DAYS")], by = "Sample_ID")
pathway_clin$OS_STATUS <- ifelse(pathway_clin$OS_STATUS == "Yes", 1,0)
pathway_clin$PFS_STATUS <- ifelse(pathway_clin$PFS_STATUS == "Yes", 1,0)
pathway_clin$OS_DAYS <- as.numeric(pathway_clin$OS_DAYS)
pathway_clin$PFS_DAYS <- as.numeric(pathway_clin$PFS_DAYS)

#####
top20_pathway <- merge(top20_matrix, pathway_matrix, by = "Sample_ID")
top20_pathway_clin <- merge(top20_pathway, clin[,c("Sample_ID", "OS_STATUS","OS_DAYS")], by = "Sample_ID")

###cox
source("Cox.R")

top20_pathway_clin$OS_STATUS <- ifelse(top20_pathway_clin$OS_STATUS == "Yes", 1,0)
top20_pathway_clin$OS_DAYS <- as.numeric(top20_pathway_clin$OS_DAYS)

cox_clin<-top20_pathway_clin[,c("Sample_ID","OS_STATUS","OS_DAYS")]  ###OS_stats os_month
dc<-top20_pathway_clin[,c(2:(ncol(top20_pathway_clin)-2))]

rownames(dc)<-top20_pathway_clin[,1]
dc<-t(dc)

colnames(cox_clin)<-c("sample","Events","Survival")
cox_uni<-as.data.frame(getUniOrMultiCOXAnalysis(dc,cox_clin,method="Univariate"))
write.table(cox_uni,"new_os_cox_univariate_top20_pathway.txt",sep="\t",quote=F,row.names = F)

#####co_mut
comut_input <- top20_pathway
rownames(comut_input) <- comut_input$Sample_ID
comut_input <- comut_input[,-1]
value_group <- combn(colnames(comut_input),2)
for (i in 1:ncol(value_group)){
  value_1 <- value_group[1,i]
  value_2 <- value_group[2,i]
  comut_input[,paste0(value_1,":",value_2)] <- ifelse(comut_input[,value_1] == 1 & comut_input[,value_2] == 1, "1","0")}

comut_matrix <- comut_input[,42:ncol(comut_input)]
mut_summary <- list()
for (j in colnames(comut_matrix)){
  mut_table <- as.data.frame(table(comut_matrix[,j]))
  if (mut_table[1,2] > 142){
    next
  } else {
    mut_summary[[j]] <- mut_table}}
comut_matrix$Sample_ID <- rownames(comut_matrix)
comut_matrix_clin <- merge(comut_matrix, clin[,c("Sample_ID", "OS_STATUS","OS_DAYS")], by = "Sample_ID")

###cox
comut_matrix_clin$OS_STATUS <- ifelse(comut_matrix_clin$OS_STATUS == "Yes", 1,0)
comut_matrix_clin$OS_DAYS <- as.numeric(comut_matrix_clin$OS_DAYS)

cox_clin<-comut_matrix_clin[,c("Sample_ID","OS_STATUS","OS_DAYS")]  ###OS_stats os_month
dc<-comut_matrix_clin[,c(2:(ncol(comut_matrix_clin)-2))]

rownames(dc)<-comut_matrix_clin[,1]
dc<-t(dc)

colnames(cox_clin)<-c("sample","Events","Survival")
cox_uni<-as.data.frame(getUniOrMultiCOXAnalysis(dc,cox_clin,method="Univariate"))
write.table(cox_uni,"new_os_cox_univariate_comut.txt",sep="\t",quote=F,row.names = F)

cox_uni_filter <- cox_uni[cox_uni$geneName %in% names(mut_summary),]
write.table(cox_uni_filter,"new_os_cox_univariate_comut_filter_panel.txt",sep="\t",quote=F,row.names = F)

#####
cox_uni_sig <- cox_uni[cox_uni$`Pr(>|Z|)[pval]` < 0.1,]

comut_matrix_clin_pfs <- merge(comut_matrix_clin, clin[,c("Sample_ID","PFS_STATUS","PFS_DAYS")], by = "Sample_ID")
comut_matrix_clin_pfs$PFS_STATUS <- ifelse(comut_matrix_clin_pfs$PFS_STATUS == "Yes",1, 0)
comut_matrix_clin_pfs$PFS_DAYS <- as.numeric(comut_matrix_clin_pfs$PFS_DAYS)

pdf("ARID1B_DDRpathway_survival_pfs_0510.pdf", width = 6, height = 6, onefile = F)
fit=survfit(Surv(as.numeric(PFS_DAYS),PFS_STATUS )~ comut_matrix_clin_pfs$`ARID1B:DDR_pathway`, data = comut_matrix_clin_pfs)
p <- ggsurvplot(fit, data = comut_matrix_clin_pfs,font.title=c(20, "bold"),conf.int=F, pval=TRUE, legend.title ="",ylab=c("Progression Free Survival"),surv.median.line = "hv",
                legend.labs = c("ARID1B_DDRpathway_WT","ARID1B_DDRpathway_MT"),palette = c("#00828099","#bb002199"),)
res_cox<-coxph(Surv(as.numeric(PFS_DAYS),PFS_STATUS )~ comut_matrix_clin_pfs$`ARID1B:DDR_pathway`, data = comut_matrix_clin_pfs)
p$plot = p$plot + 
  ggplot2::annotate("text",x = 100, y = 0.15,
                    label = paste("HR :",round(summary(res_cox)$conf.int[1],2))) + 
  ggplot2::annotate("text",x = 100, y = 0.10,
                    label = paste("(","95%CI:",round(summary(res_cox)$conf.int[3],2),"-",round(summary(res_cox)$conf.int[4],2),")",sep = ""))
dev.off()


pdf("ARID1B_DDRpathway_survival_os_0511.pdf", width = 6, height = 6, onefile = F)
fit=survfit(Surv(as.numeric(OS_DAYS),OS_STATUS )~ comut_matrix_clin_pfs$`ARID1B:DDR_pathway`, data = comut_matrix_clin_pfs)
p <- ggsurvplot(fit, data = comut_matrix_clin_pfs,font.title=c(20, "bold"),conf.int=F, pval=TRUE, legend.title ="",ylab=c("Overall Survival"),surv.median.line = "hv",
                legend.labs = c("ARID1B_DDRpathway_WT","ARID1B_DDRpathway_MT"),palette = c("#00828099","#bb002199"),)
res_cox<-coxph(Surv(as.numeric(OS_DAYS),OS_STATUS )~ comut_matrix_clin_pfs$`ARID1B:DDR_pathway`, data = comut_matrix_clin_pfs)
p$plot = p$plot + 
  ggplot2::annotate("text",x = 100, y = 0.15,
                    label = paste("HR :",round(summary(res_cox)$conf.int[1],2))) + 
  ggplot2::annotate("text",x = 100, y = 0.10,
                    label = paste("(","95%CI:",round(summary(res_cox)$conf.int[3],2),"-",round(summary(res_cox)$conf.int[4],2),")",sep = ""))
dev.off()

#####
uni_result <- cox_uni
for (n in 1:nrow(uni_result)){
  group <- strsplit(uni_result[n,"geneName"], ":", fixed= T)[1]
  uni_result[n, "Value_1"] <- group[[1]][1]
  uni_result[n, "Value_2"] <- group[[1]][2]}

cox_uni_p <- uni_result[,c("Value_1","Value_2","Pr(>|Z|)[pval]")]
colnames(cox_uni_p)[3] <- c("P_value")
cox_uni_hr <- uni_result[,c("Value_1","Value_2","exp(coef)[HR]")]
colnames(cox_uni_hr)[3] <- c("HR")
 
#####spread
uni_p_spread <- spread(cox_uni_p, Value_2, P_value)
rownames(uni_p_spread) <- uni_p_spread$Value_1
uni_p_spread <- uni_p_spread[,-1]
uni_hr_spread <- spread(cox_uni_hr, Value_2, HR)
rownames(uni_hr_spread) <- uni_hr_spread$Value_1
uni_hr_spread <- uni_hr_spread[,-1]
uni_p_spread["WNT_pathway",] <- NA
uni_hr_spread["WNT_pathway",] <- NA
uni_p_spread[,"TP53"] <- NA
uni_hr_spread[,"TP53"] <- NA

uni_p_spread <- uni_p_spread[match(colnames(uni_p_spread), rownames(uni_p_spread)),]
uni_hr_spread <- uni_hr_spread[match(colnames(uni_hr_spread), rownames(uni_hr_spread)),]


for (i in 1:ncol(uni_p_spread)){
  for (j in 1:nrow(uni_p_spread)){
    if (is.na(uni_p_spread[i,j]) == FALSE){
      uni_p_spread[j,i] <- uni_p_spread[i,j]
    } else {
      uni_p_spread[i,j] <- uni_p_spread[j,i]}}}
uni_p_spread[!upper.tri(uni_p_spread, diag = TRUE)] <- 0

for (i in 1:ncol(uni_hr_spread)){
  for (j in 1:nrow(uni_hr_spread)){
    if (is.na(uni_hr_spread[i,j]) == FALSE){
      uni_hr_spread[j,i] <- uni_hr_spread[i,j]
    } else {
      uni_hr_spread[i,j] <- uni_hr_spread[j,i]}}}
uni_hr_spread[!lower.tri(uni_hr_spread, diag = TRUE)] <- 0

uni_merge <- data.frame()
for (m in 1:ncol(uni_p_spread)){
  for (n in 1:nrow(uni_p_spread)){
    uni_merge[m,n] <- uni_hr_spread[m,n]
    uni_merge[n,m] <- uni_p_spread[n,m]}}
rownames(uni_merge) <- rownames(uni_p_spread)
colnames(uni_merge) <- colnames(uni_p_spread)
for (i in 1:ncol(uni_merge)) {uni_merge[,i] <- as.numeric(uni_merge[,i])}
uni_merge <- as.matrix(uni_merge)

uni_p_spread <- as.data.frame(uni_p_spread)
for (i in 1:ncol(uni_p_spread)) {uni_p_spread[,i] <- as.numeric(uni_p_spread[,i])}
uni_p_spread <- as.matrix(uni_p_spread)
uni_p_spread[uni_p_spread == 0] <- NA
uni_hr_spread <- as.data.frame(uni_hr_spread)
for (i in 1:ncol(uni_hr_spread)) {uni_hr_spread[,i] <- as.numeric(uni_hr_spread[,i])}
uni_hr_spread <- as.matrix(uni_hr_spread)
uni_hr_spread[uni_hr_spread == 0] <- NA

#####mut_freq
top20_pathway_freq <- top20_pathway_clin[,1:(ncol(top20_pathway_clin)-2)]
rownames(top20_pathway_freq) <- top20_pathway_freq$Sample_ID
top20_pathway_freq <- as.data.frame(t(top20_pathway_freq[,-1]))
top20_pathway_freq$Num <- rowSums(top20_pathway_freq[,1:ncol(top20_pathway_freq)])
top20_pathway_freq$Freq <- top20_pathway_freq$Num/162

#####
single_uni <- read.table("new_os_cox_univariate_top20_pathway.txt", sep = '\t', quote = '', check.names = F, header = T, row.names = NULL)
single_value <- single_uni[,c("geneName","exp(coef)[HR]","Pr(>|Z|)[pval]")]
rownames(single_value) <- single_value$geneName
single_value <- single_value[,-1]
colnames(single_value) <- c("HR","P")

single_value <- single_value[match(rownames(uni_merge), rownames(single_value)),]
top20_pathway_freq <- top20_pathway_freq[match(rownames(single_value), rownames(top20_pathway_freq)),]
identical(rownames(single_value), rownames(top20_pathway_freq))
single_value$Mut_Freq <- top20_pathway_freq$Freq
single_value <- as.matrix(single_value)

#####complexheatmap
col1 = colorRamp2(c( 0,1,2), c("white","#FF6347","#1E90FF"))
col2 = colorRamp2(c( 0,0.5, 1), c( "purple", "orange","white"))
col3 = colorRamp2(c( 0,0.5, 1), c("white", "#40E0D0", "#DB7093"))

ht1 <- Heatmap(uni_merge, rect_gp = gpar(type = "none"), show_heatmap_legend = FALSE,border = "black",
               cluster_rows = FALSE, cluster_columns = FALSE,row_names_gp = gpar(fontsize = 10),column_names_gp = gpar(fontsize = 10),
               row_names_side = "left",column_names_side = 'top',
               
               layer_fun = function(j, i, x, y, w, h, fill) {
                 l = i > j
                 grid.rect(x[l], y[l], w[l], h[l], 
                           gp = gpar(fill = col1(pindex(uni_merge, i[l], j[l])), col = "grey90"))
                 l = i < j
                 grid.rect(x[l], y[l], w[l], h[l], 
                           gp = gpar(fill = col2(pindex(uni_merge, i[l], j[l])), col = "grey90"))
               },
               cell_fun = function(j, i, x, y, width, height, fill) {
                 grid.text(sprintf("%.2f", uni_merge[i, j]), x, y, gp = gpar(fontsize = 10))
               })
draw(ht1, heatmap_legend_list = list(
  Legend(title = "Hazard Ratio", col_fun = col1),
  Legend(title = "P", col_fun = col2)))

ht2 <- Heatmap(single_value[,"HR"], col = col1, name = "HR", cluster_rows = F, 
               row_names_gp = gpar(fontsize = 10),column_names_gp = gpar(fontsize = 10),
               row_names_side = "left",column_names_side = 'top',
               cell_fun = function(j, i, x, y, width, height, fill) {
                 grid.text(sprintf("%.2f", single_value[i, "HR"]), x, y, gp = gpar(fontsize = 10))})  
ht3 <- Heatmap(single_value[,"P"], col = col2, name = "P", cluster_rows = F, 
               row_names_gp = gpar(fontsize = 10),column_names_gp = gpar(fontsize = 10),
               row_names_side = "left",column_names_side = 'top',
               cell_fun = function(j, i, x, y, width, height, fill) {
                 grid.text(sprintf("%.2f", single_value[i, "P"]), x, y, gp = gpar(fontsize = 10))})  
ht4 <- Heatmap(single_value[,"Mut_Freq"], col = col3, name = "Mut_freq", cluster_rows = F, 
               row_names_gp = gpar(fontsize = 10),column_names_gp = gpar(fontsize = 10),
               row_names_side = "left",column_names_side = 'top',
               cell_fun = function(j, i, x, y, width, height, fill) {
                 grid.text(sprintf("%.2f", single_value[i, "Mut_Freq"]), x, y, gp = gpar(fontsize = 10))})  
pdf("new_top20_pathway_uni_panel.pdf", width = 20, height = 15, onefile = F)
draw(ht4 + ht2 + ht3 + ht1)
dev.off()
