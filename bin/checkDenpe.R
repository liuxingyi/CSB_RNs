checkDenpe<-function(GSEGetted){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  if (!file.exists("../result")){
    dir.create("../result")
    print(paste0("../result不存在，已创建！"))
  }
  if (!file.exists("../logs")){
    dir.create("../logs")
    print(paste0("../logs不存在，已创建！"))
  }
  if (!file.exists(paste0("../result/",GSEGetted))){
    dir.create(paste0("../result/",GSEGetted))
    print(paste0("../result/",GSEGetted,"不存在已创建！"))
  }
  if (!file.exists(paste0("../logs/",GSEGetted))){
    dir.create(paste0("../logs/",GSEGetted))
    print(paste0("../logs/",GSEGetted,"不存在已创建！"))
  }
  fileDirs<-c(paste0("../result/",GSEGetted,"/"),paste0("../logs/",GSEGetted,"/"))
  return(fileDirs)
}
#======================检查heatmap====================
checkHeatmap<-function(){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  if (!requireNamespace(c("ALL"), quietly = TRUE))
      BiocManager::install(c("ALL"))
  if (!requireNamespace(c("gplots"), quietly = TRUE))
    BiocManager::install(c("gplots"))
  if (!requireNamespace(c("pheatmap"), quietly = TRUE))
    BiocManager::install(c("pheatmap"))
  library("ALL")
  library("gplots")
  library("pheatmap")
  
}
#======================检查volcano====================
checkvolcano<-function(){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  if (!requireNamespace(c("colorspace"), quietly = TRUE))
    BiocManager::install(c("colorspace"))
  if (!requireNamespace(c("ggplot2"), quietly = TRUE))
    BiocManager::install(c("ggplot2"))
	  if (!requireNamespace(c("ggrepel"), quietly = TRUE))
    BiocManager::install(c("ggrepel"))
  library("ggplot2")
  library("ggrepel")
}
#====================GO富集分析======================
checkGO<-function(){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  if (!requireNamespace(c("AnnotationHub"), quietly = TRUE))
    BiocManager::install(c("AnnotationHub"))
  if (!requireNamespace(c("org.Hs.eg.db"), quietly = TRUE))
    BiocManager::install(c("org.Hs.eg.db"))
  if (!requireNamespace(c("graph"), quietly = TRUE))
    BiocManager::install(c("graph"))
  if (!requireNamespace(c("Rgraphviz"), quietly = TRUE))
    BiocManager::install(c("Rgraphviz"))
  if (!requireNamespace(c("topGO"), quietly = TRUE))
    BiocManager::install(c("topGO"))
  if (!requireNamespace(c("Rgraphviz"), quietly = TRUE))
    BiocManager::install(c("Rgraphviz"))
  if (!requireNamespace(c("clusterProfiler"), quietly = TRUE))
    BiocManager::install(c("clusterProfiler"))
  suppressMessages(library("AnnotationHub"))
  library("org.Hs.eg.db")
  library("clusterProfiler")
  library("graph")
  library("Rgraphviz")
  library("topGO")
}
#====================kegg、go分析=========================
checkOntDir<-function(GSEGetted,ont){
  if (!file.exists(paste0("../result/",GSEGetted,"/",ont,"/"))){
    dir.create(paste0("../result/",GSEGetted,"/",ont,"/"))
    print(paste0("../result/",GSEGetted,"/",ont,"/","不存在已创建！"))
  }
  return(paste0("../result/",GSEGetted,"/",ont,"/"))
}
#==================limma分析==============================
checkLimma<-function(){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  if (!requireNamespace("limma", quietly = TRUE))
    BiocManager::install(c("limma"))
  suppressMessages(library("limma"))
}
