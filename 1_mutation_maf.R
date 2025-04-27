###mutation
library(maftools)

a <-  list.files("./")
mutation_list = list()
for (i in a){
  # print(i)
  mutation <- read.table(i, sep = '\t', quote = '', row.names = NULL, header = T, check.names = F,fill = T, fileEncoding = "utf-8")
  mutation <- mutation[which(mutation$FILTER_flag == "True" & mutation$QC_flag == "True" & mutation$MAF_flag == "True" & mutation$Black_flag == "False"),]
  mutation <- mutation[which(mutation$Black_flag == "False"),]
  mutation <- mutation[,c("Chr","Start","End","Ref","Alt","Func.refGene", "Gene.refGene","ExonicFunc.refGene","AAChange.refGene", "Sample", "GeneDetail.refGene")]
  mutation <- as.data.frame(t(mutation))
  mutation_list[[i]] = mutation
}
mutation_sample <- as.data.frame(mutation_list)
mutation_sample <- as.data.frame(t(mutation_sample))
# write.table(mutation_sample,"somatic_data.txt",sep="\t",quote=F,row.names = F)

somatic_mutation <- read.table("somatic_data.txt", sep = "\t",quote = '', header = T)
mutation_maf = annovarToMaf(annovar= "somatic_data.txt", Center = NULL, refBuild = "hg19", tsbCol = NULL, table = "refGene", ens2hugo = TRUE, basename = NULL, sep = "\t", MAFobj = FALSE, sampleAnno = NULL)
mutation_maf$Tumor_Sample_Barcode <- mutation_maf$Sample
mutation_maf <- mutation_maf[,-18]
# write.table(mutation_maf,file="somatic_mutation.maf",sep="\t",quote=F,row.names = FALSE)
length(table(mutation_maf$Tumor_Sample_Barcode))

###oncoplot
library(openxlsx)
library(maftools)
#基因突变临床数据-----
input<-read.table("./somatic_maf_input.txt",sep = "\t",quote = "",header = T)
nonsys<-c("Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins","Missense_Mutation","Nonsense_Mutation","Splice_Site")
input<-input[which(input$Variant_Classification %in% nonsys),]

sample_info<-read.xlsx("clin_samples.xlsx")

write.table(input[which(input$Tumor_Sample_Barcode %in% sample_info[which(sample_info$RECIST == "PR"),1]),],
            file="PR_somatic_mutation.maf",sep="\t",quote=F,row.names = FALSE)
write.table(input[which(input$Tumor_Sample_Barcode %in% sample_info[which(sample_info$RECIST %in% c("PD","SD")),1]),],
            file="PDSD_somatic_mutation.maf",sep="\t",quote=F,row.names = FALSE)
#write.table(input,file="somatic_mutation_nonsys.maf",sep="\t",quote=F,row.names = FALSE)

