#=============得到两组之间的p值和fold change==============
getFC_P<-function(gseMat,sampleCategory){
  tumor<-data.frame(gseMat[,sampleCategory=="1"])
  print(cbind(colnames(tumor),"非野生型样本"))
  notumor<-data.frame(gseMat[,sampleCategory=="0"])
  print(cbind(colnames(notumor),"野生型样本"))
  FC_P<-NULL
  matFC<-NULL
  matP<-NULL
  matP.adjust<-NULL
  FCupOrDown<-NULL
  #if(tPMethod=="two.side")
  for(i in seq(1,nrow(gseMat))){
    p<-t.test(notumor[i,],tumor[i,])
    matP<-c(matP,p$p.value)
    #fdr修正
    matFC<-c(matFC,2^(sum(tumor[i,])/ncol(tumor[i,])-sum(notumor[i,])/ncol(notumor[i,])))
    if(matFC[i]>1) FCupOrDown=c(FCupOrDown,1) else FCupOrDown=c(FCupOrDown,-1)
    #====================进度条===================
    if(i%%round(nrow(gseMat)/100)==0){
      print(paste0(round(i/round(nrow(gseMat)/100)),"%的表达谱分析已完成，第",i,"条!!!!"))
    }
    if(i==nrow(gseMat)){
     print(paste0("100%的表达谱分析已完成，第",i,"条!!!!"))
    }
    #====================进度条===================
  }
  matP<-p.adjust(matP,method = "fdr",n=length(matP))
  FC_P<-data.frame(geneSYMOL=rownames(gseMat),foldChange=round(matFC,5),
                   P_value=matP,FCupOrDown)
  return(FC_P)
}
#=============limma包得到两组之间的p值和fold change==============
getFC_P_limma<-function(gseMat,sampleCategory){
  checkLimma()
  design <- model.matrix(~0+factor(sampleCategory))
  rownames(design)<-colnames(gseMat)
  fit=lmFit(gseMat,design)
  fit=eBayes(fit)
  gseMat<-topTable(fit,coef=1,adjust='BH',n=Inf)
  FC_P<-data.frame(geneSYMOL=rownames(gseMat),foldChange=2^gseMat$logFC,
                   P_value=gseMat$P.Value)
  return(FC_P)
}

#=====================获取up and down标签==========
getYouUpDown<-function(matFC_P,fc_yours,p_yours){
  threshold<-NULL
  for(i in seq(1,nrow(matFC_P))){
    if(matFC_P$foldChange[i]>fc_yours & matFC_P$P_value[i]<p_yours){
      threshold<-c(threshold,"up")
    }else if(matFC_P$foldChange[i]<1/fc_yours & matFC_P$P_value[i]<p_yours){
      threshold<-c(threshold,"down")
    }else{
      threshold<-c(threshold,"nosig")
    }
  }
  return(threshold)
}
#==================火山图标注个数设置==============
labelHotspot<-function(volcanoMat,topNUM){
  volcanoMat_up<-volcanoMat[volcanoMat$Threshold=="up",]
  volcanoMat_up_gene<-tail(volcanoMat_up[order(as.numeric(as.vector(volcanoMat_up$YNegalog10P))),],topNUM)
  volcanoMat_down<-volcanoMat[volcanoMat$Threshold=="down",]
  volcanoMat_down_gene<-tail(volcanoMat_down[order(as.numeric(as.vector(volcanoMat_down$YNegalog10P))),],topNUM)
  volcanoMat$labelGHY<-""
  topNUMgene<-c(as.vector(volcanoMat_up_gene$geneSYMOL),as.vector(volcanoMat_down_gene$geneSYMOL))
  volcanoMat$labelGHY[match(topNUMgene,volcanoMat$geneSYMOL)]<-topNUMgene
  return(volcanoMat$labelGHY)
  }