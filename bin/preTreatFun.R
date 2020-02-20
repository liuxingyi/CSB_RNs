#=============================================================
#输入GSExxxxx,返回下载matrix数据的ftp号
GSEToUrl<-function(GSEGetted){
  if(nchar(GSEGetted)==9){
    dirNnn<-paste(substring(GSEGetted,1,6),"nnn/",sep="")
  }else{
    dirNnn<-paste(substring(GSEGetted,1,5),"nnn/",sep="")
  }
  
  prefixUrl<-"ftp://ftp.ncbi.nlm.nih.gov/geo/series/"
  matName<-paste(GSEGetted,"_series_matrix.txt.gz",sep="")
  url<-paste(prefixUrl,dirNnn,GSEGetted,"/matrix/",matName,sep="")
  return(url)
}
#============================================================
#选中位数（必须存在的一个）
mean_one<-function(number_array){
  number_array_distance<-abs(number_array-mean(number_array))
  return(number_array[number_array_distance%in%min(number_array_distance)][1])
}
#============================================================
#输入未注释的matrix矩阵，输出注释的matrix文件
annotaMat<-function(geneMat,gplPlatform){
  geneMat<-na.omit(geneMat)
  matDb<-read.table("../config_file/chip_platform",header=T)
  annotation_name<-sub("\\s+","",as.character(matDb[matDb[,1]==gplPlatform,3]))
  if(length(annotation_name) != 0){
    annotation_db<-paste(annotation_name,".db",sep="")
    if (!requireNamespace(annotation_db)){
      BiocManager::install(annotation_db)
    }
    library(annotation_db,character.only = T)
    annotation_db_SYMBOL<-paste(annotation_name,"SYMBOL",sep="")
    annotation_db_cmd<-paste("ids<-toTable(",annotation_db_SYMBOL,")")
    eval(parse(text = annotation_db_cmd))
    
    affySetExprs<-merge(ids,geneMat,by.x="probe_id",by.y="ID_REF")
    
    exprSet<-affySetExprs[affySetExprs$probe_id %in% ids$probe_id,]
    #check no annnotation about probes
    ids<-ids[ids$probe_id %in% exprSet$probe_id,]
    #acoording to order of befor dataset,sorting probe_id
    ids<-ids[match(exprSet$probe_id,ids$probe_id),]
    rownames(exprSet)<-exprSet[,1]
    tmp<-by(exprSet[,3:ncol(exprSet)],exprSet$symbol,function(x)rownames(x)[which(rowMeans(x)%in%mean_one(rowMeans(x)))])
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
}
