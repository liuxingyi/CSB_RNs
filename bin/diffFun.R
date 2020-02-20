writeLogs<-function(GSEGetted,GSEplatform,sampleCategory,p_yours,fc_yours){
  MessageOur<-paste0("\n=========================\n","gseID:",GSEGetted,
                     "\nGSEplatform:",GSEplatform,"\nsampleCategory:",paste0(sampleCategory,collapse =" ")[1],
                     "\np_value:",p_yours,"\nFC:",fc_yours)
  write.table(MessageOur,file=paste0("../logs/",GSEGetted,"finished.log"),append = T,col.names = F)
  write.table(MessageOur,file=paste0("../result/",GSEGetted,"/parameterIntroFile.txt"),append = F,col.names = F)
}
#============================将箱式图拉到一个水平===========================
normalMedian<-function(gseMat){
  meanMedian<-0
  newGseMat<-NULL
  for(i in seq(1,ncol(gseMat))){
    meanMedian<-meanMedian+median(gseMat[,i])
  }
  meanMedian<-meanMedian/ncol(gseMat)
  for(i in seq(1,ncol(gseMat))){
    newGseMat<-cbind(newGseMat,gseMat[,i]-(median(gseMat[,i])-meanMedian))
  }
  newGseMat<-data.frame(newGseMat)
  rownames(newGseMat)<-rownames(gseMat)
  return(newGseMat)
}
#===========================注释注释包不能注释的注释的平台==================
preDealDowInfor<-function(annofileDir,GSEplatform,idcol,symbolcol){
  annoData<-read.table(annofileDir,header=T,fill=T,sep = "\t",comment.char = "#")
  annoData<-data.frame(probe_id=annoData$ID,symbol=annoData$GENE_SYMBOL,label=annoData$GENE_SYMBOL=="")
  annoData<-annoData[annoData$label==FALSE,c("probe_id","symbol")]
  write.table(annoData,paste0("../config_file/",GSEplatform,".idSymbol"),col.names = T,quote = F,row.names = F)
}
#============================id symbol对应文件注释=============================
fileannotaMat<-function(geneMat,idsymbolDir){
  ids<-read.table(idsymbolDir,header=T,fill=F,sep = " ")
  affySetExprs<-merge(ids,geneMat,by.x="probe_id",by.y="ID_REF")
  
  exprSet<-affySetExprs[affySetExprs$probe_id %in% ids$probe_id,]
  #check no annnotation about probes
  ids<-ids[ids$probe_id %in% exprSet$probe_id,]
  #acoording to order of befor dataset,sorting probe_id
  ids<-ids[match(exprSet$probe_id,ids$probe_id),]
  rownames(exprSet)<-exprSet[,1]
  tmp<-by(exprSet[,3:ncol(exprSet)],ids$symbol,function(x)rownames(x)[which(rowMeans(x)%in%mean_one(rowMeans(x)))])
  probes<-as.character(tmp)
  print("已经去除某个gene多探针的情况(多探针情况选择最接近均值的探针)")
  print(paste0("未处理前的探针情况是:",length(exprSet$symbol),"个"))
  print("未处理前的探针情况----基因含有探针的分布状况:")
  print(table(data.frame(table(exprSet$symbol))$Freq))
  print("处理去重后的探针情况是:")
  print(table(data.frame(table(probes))$Freq))
  expeSet<-exprSet[rownames(exprSet) %in% probes,]
  ids<-ids[ids$probe_id %in% probes,]
  rownames(expeSet)<-ids[match(expeSet$probe_id,ids$probe_id),2]
  expeSet<-expeSet[,-1]
  expeSet<-expeSet[,-1]
  #clear name space
  return(expeSet)
}