

library("survival")
 library("survminer")

getPrintFile<-function(coxph.fit,gene,GeneCoefResult){
  result<-summary(coxph.fit);#txt<-paste(fileName,"_cox.txt",sep="")  # print(result)
  header=t(c("geneName","coef","exp(coef)[HR]","se(coef)","Z[coef/se]","Pr(>|Z|)[pval]","lower.95","upper.95","Likelihood ratio test","Wald test","Log rank test"))
  count<-rownames(result[[7]])
  
  if(length(count)==1){
    ##单因素输出
    aa<-result$coefficients
    bb<-result$conf.int
    temp<-t(c(gene,aa,bb[c(3,4)],result[[9]][3],result[[12]][3],result[[10]][3]))
    if(is.null(GeneCoefResult)==TRUE){
      GeneCoefResult<-rbind(GeneCoefResult,temp)
      colnames(GeneCoefResult)<-header
    }else{
      GeneCoefResult<-rbind(GeneCoefResult,temp)########合并所有基因的风险得分结果
    }
    
    
  }else{
    ##多因素输出
    aa<-result$coefficients
    bb<-result$conf.int
    n<-dim(aa)[1]
    GeneCoefResult<-data.frame(gene[1:n],aa,bb[,c(3,4)],rep(result[[9]][3],times=n),rep(result[[12]][3],times=n),rep(result[[10]][3],times=n))
    colnames(GeneCoefResult)<-header
    
    
  }
  
  
  return(GeneCoefResult)
}
####################COX survival analysis######################
getUniOrMultiCOXAnalysis<-function(subprof,clin,method){
  ###将表达谱和临床数据合并到一个文件###
  subprof<-t(subprof);
  samples<-rownames(subprof);
  mm<-match(samples,as.character(clin[,1]))
  newsubprof<-subprof[which(is.na(mm)==FALSE), ,drop=FALSE];
  n1<-dim(newsubprof)[2]
  clinInfor<-clin[mm[which(is.na(mm)==FALSE)], ,drop=FALSE];
  n2<-dim(clinInfor)[2]
  cox.data<-cbind(newsubprof,clinInfor[,-1])[,c((n1+1):(n1+2),(1:n1))]
  numdata<-apply(unname(cox.data),c(1,2),as.numeric)
  cox.data<-apply(cox.data,c(1,2),as.numeric)
  cox.data<-as.data.frame(cox.data)
  ###单因素或多因素COX###
  gname<-colnames(cox.data)[-c(1,2)]
  if(method=="Univariate"){
    GeneCoefResult<-c()
    for(i in 1:length(gname)){
      print(i)
      x<-numdata[,2+i];
      if(sd(na.omit(x))==0){next
      }else{
        
        
        coxph.fit<-coxph(Surv(Survival,Events)~x,data = cox.data) 
        GeneCoefResult<-getPrintFile(coxph.fit,gname[i],GeneCoefResult);
        #detach(cox.data)
      }
    }
  }else if(method=="Multivariate"){
    GeneCoefResult<-c();
    data<-numdata;attach(cox.data);nn<-dim(numdata)[2];
    coxph.fit<-coxph(Surv(Survival,Events)~data[,3:nn]);
    GeneCoefResult<-getPrintFile(coxph.fit,gname,GeneCoefResult);
    detach(cox.data)
  } 
  #UpGeneCoefResult <- GeneCoefResult [which(GeneCoefResult[,3]>1), ,drop=FALSE]
  #LoGeneCoefResult <- GeneCoefResult [which(GeneCoefResult[,3]<1), , drop=FALSE]
  return(GeneCoefResult)
}



# getCutoffValue<-function(coxRes,subprof){
#   coef<-coxRes[,c("geneName","coef"),drop=FALSE]
#   inter<-intersect(gsub("gene.","",as.character(coef[,1])),as.character(rownames((subprof))))
#   #nn<-match(gsub("gene.","",as.character(coef[,1])),as.character(subprof[,1]))
#   nn<-match(as.character(inter),gsub("gene.","",as.character(coef[,1])))
#   nn<-nn[which(is.na(nn)==FALSE)]
#   #nnn<-match(as.character(subprof[,1]),gsub("gene.","",as.character(coef[,1])))
#   nnn<-match(as.character(inter),as.character(rownames((subprof))))
#   nnn<-nnn[which(is.na(nnn)==FALSE)]
#   coefExp<-c();numprof<-apply(unname(subprof),c(1,2),as.numeric)
#   for(i in 1:dim(numprof)[2]){
#     coefExp<-c(coefExp,sum(as.numeric(as.character(coef[nn,2]))*as.numeric(as.character(numprof[nnn,i]))))
#   }
#   
#   #median<-mean(coefExp)
#   median<-quantile(coefExp,  probs = c(0.5))
#   return(median)     
# }
##########################################
getCutoffValue<-function(coxRes,subprof){
  coef<-coxRes[,c("geneName","coef"),drop=FALSE]
  inter<-intersect(gsub("gene.","",as.character(coef[,1])),as.character(rownames((subprof))))
  #nn<-match(gsub("gene.","",as.character(coef[,1])),as.character(subprof[,1]))
  nn<-match(as.character(inter),gsub("gene.","",as.character(coef[,1])))
  nn<-nn[which(is.na(nn)==FALSE)]
  #nnn<-match(as.character(subprof[,1]),gsub("gene.","",as.character(coef[,1])))
  nnn<-match(as.character(inter),as.character(rownames((subprof))))
  nnn<-nnn[which(is.na(nnn)==FALSE)]
  coefExp<-c();numprof<-apply(unname(subprof),c(1,2),as.numeric)
  for(i in 1:dim(numprof)[2]){
    coefExp<-c(coefExp,sum(as.numeric(as.character(coef[nn,2]))*as.numeric(as.character(numprof[nnn,i]))))
  }
  
 median<-median(na.omit(coefExp))
# median<-mean(coefExp)
  #median<-quantile(coefExp,  probs = c(0.5))
  return(median)     
}
#####

getCoefExpCluster<-function(coxRes,subprof2,subprof){
  coef<-coxRes[,c("geneName","coef"), drop=FALSE ]
  inter<-intersect(gsub("gene.","",as.character(coef[,1])),as.character(rownames((subprof2))))
  #nn<-match(gsub("gene.","",as.character(coef[,1])),as.character(subprof2[,1]))
  nn<-match(as.character(inter),gsub("gene.","",as.character(coef[,1])))
  nn<-nn[which(is.na(nn)==FALSE)]
  #nnn<-match(as.character(subprof2[,1]),gsub("gene.","",as.character(coef[,1])))
  nnn<-match(as.character(inter),as.character(rownames((subprof2))))
  nnn<-nnn[which(is.na(nnn)==FALSE)]
  coefExp2<-c();numprof<-apply(unname(subprof2),c(1,2),as.numeric)
  for(i in 1:dim(numprof)[2]){
    coefExp2<-c(coefExp2,sum(as.numeric(as.character(coef[nn,2]))*as.numeric(as.character(numprof[nnn,i]))))
  }
  cutoff<-getCutoffValue(coxRes,subprof)
  riskscore<-coefExp2
  coefExp22<-coefExp2
  coefExp22[which(coefExp2>=cutoff)]<-1
  coefExp22[which(coefExp2<cutoff)]<-0
  return(cbind(colnames(subprof2),coefExp22,riskscore))
}


getKMSurAnalyPvalue<-function(subprof2,clin2,label2){
  
  subprof2<-t(subprof2);#colnames(subprof2)<-paste("gene.",as.character(subprof2[1,]),sep="");subprof2<-subprof2[-1,]
  colnames(subprof2)<-paste("gene.",as.character(colnames((subprof2))),sep="");#subprof2<-subprof2[-1,]
  samples<-rownames(subprof2);
  inter<-intersect(as.character(clin2[,1]),samples)
  mm<-match(as.character(inter),as.character(clin2[,1]))
  nn<-match(as.character(inter),samples)
  
  newsubprof<-subprof2[nn[which(is.na(nn)==FALSE)],,drop=FALSE];n1<-dim(newsubprof)[2]
  clinInfor<-clin2[mm[which(is.na(mm)==FALSE)],];n2<-dim(clinInfor)[2]
  KMdata<-data.frame(newsubprof,clinInfor[,-1])[,c("Survival","Events")]
  interT<-intersect(label2[,1],rownames(KMdata))
  mm2<-match(as.character(interT),rownames(KMdata))
  mm2<-mm2[which(is.na(mm2)==FALSE)]
  nn2<-match(as.character(interT),label2[,1])
  nn2<-nn2[which(is.na(nn2)==FALSE)]
  colnames(label2)[1]<-"Sample"
  KMdata<-cbind(label2[nn2,],KMdata[mm2,])[,c(1,4:5,2)]
  rownames(KMdata)<-NULL;
  #colnames(KMdata)<-c("Sample","Survival","Events","coefExp")
  ####
  ####write.table(KMdata,"Survival_ClusterInfo.txt",col.names=T,row.names=F,sep="\t",quote=F)
  ####
  attach(KMdata);
  coefExp<-as.numeric(KMdata[,4])
  fit<-survfit(Surv(Survival,as.numeric(Events))~coefExp)
  
  plot(fit,  mark.time=TRUE,lty = 1,lwd=2,xlab="Survival Time (months)",
       ylab="Survival Probability",col=c(rgb(31,144,72,maxColorValue=255),
                                         rgb(240,24,40,maxColorValue=255)))
  
  dif<-survdiff(Surv(Survival,as.numeric(Events))~coefExp)
  dif.pval<-pchisq(dif$chisq,1,lower.tail = F)
  dif.pval<-round(dif.pval,digits=100)
  detach(KMdata);
  return(list(dif.pval,KMdata))
  
}
