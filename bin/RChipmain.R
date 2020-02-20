rm(list=ls())
options(warn =-1)
source("preTreatFun.R",encoding="utf-8")
source("checkDenpe.R",encoding="utf-8")
source("getDiffFun.R",encoding="utf-8")
source("diffFun.R",encoding="utf-8")
source("svmrfeFeatureRanking.R",encoding="utf-8")

GSEGetted<-"GSE28735"
GSEplatform<-"GPL6244"
sampleCategory<-c(rep(c(1,0),times=45))
p_yours<-0.01
fc_yours<-1.5



fileDirs<-checkDenpe(GSEGetted)
savedGSEFile<-paste("../result/",GSEGetted,"/",GSEGetted,"_series_matrix.txt.gz",sep="")


#=================Read the file=====================
if(file.exists(savedGSEFile)){
  print("The file has been downloaded")
}else{
  GSEUrl<-GSEToUrl(GSEGetted)
  download.file(GSEUrl,destfile=paste(savedGSEFile,sep=""))
}

data<-read.table(paste(savedGSEFile,sep=""),header=T,comment.char = "!",fill=T)
gseMat<-annotaMat(data,GSEplatform)

#rownames(data)<-data[,1]
#data<-data[-1,-1]
#=========================================================
#Download the GEO database comment file
##annofileDir<-"../config_file/GPL22763.anno"
##preDealDowInfor(annofileDir,GSEplatform,"ID","GENE_SYMBOL")
##gseMat<-fileannotaMat(data,paste0("../config_file/",GSEplatform,".idSymbol"))

#gseMat<-log2(gseMat)
if(max(gseMat)>=20){
  gseMat<-log2(gseMat)
  gseMat<-t(scale(t(gseMat),center = TRUE, scale = TRUE))
}
#write.table(gseMat[,sampleCategory==0 | sampleCategory==1],paste0(fileDirs[1],GSEGetted,"_annoMat.txt"))

#================boxplot=======================
boxplot(gseMat,col="gray",xlab="samples",ylab="express label",main="expression quentity boxplt")
pdf(paste0(fileDirs[1],GSEGetted,"_quentity_boxplt.pdf"),width = 8, height = 5)
boxplot(gseMat,col="gray",xlab="samples",ylab="express label",main="expression quentity boxplt")
dev.off()


matFC_P<-getFC_P(gseMat,sampleCategory)
write.table(matFC_P,file = paste0(fileDirs[1],GSEGetted,"matFC_P.txt"),col.names = T,row.names = F,quote = F)


#=================Volcano Plot====================
checkvolcano()
matFC_P<-matFC_P[matFC_P$P_value!=0,]
volcanoThreshold<-getYouUpDown(matFC_P,fc_yours,p_yours)
windowsFonts(TNR = windowsFont("Times New Roman"))
#testMat<-cbind(volcanoThreshold,matFC_P)
#table(testMat[testMat$volcanoThreshold=="up",]$P_value>0.05)
volcanoXLog2FC<-log2(matFC_P$foldChange)
volcanoYNegalog10P<--log10(matFC_P$P_value)
volcanoMat<-data.frame(cbind(geneSYMOL=as.vector(matFC_P$geneSYMOL),Threshold=volcanoThreshold,
        XLog2FC=volcanoXLog2FC,
        YNegalog10P=volcanoYNegalog10P))

volcanoMatLabel<-labelHotspot(volcanoMat,10)
volcanoplot<-ggplot(data=volcanoMat,aes(x=volcanoXLog2FC,
        y=volcanoYNegalog10P,colour=volcanoThreshold),family="TNR")+
        geom_point(size=2.0,na.rm = TRUE)+
        scale_color_manual(values=c("#2f5688","#BBBBBB","#CC0000"))+
        geom_vline(xintercept = c(-log2(fc_yours),log2(fc_yours)),lty=4,col="black",lwd=0.5)+
        geom_hline(yintercept=-log10(p_yours),lty=4,col="black",lwd=0.5)+
        labs(x="log2(fold change)",y="-log10(FDR)",title=paste0(GSEGetted," Volcano Plot (PDR<",p_yours,"and FC >",fc_yours,")"))+
        xlim(min(volcanoXLog2FC)-1,max(volcanoXLog2FC)+1)+
        theme(plot.title=element_text(hjust = 0.5,vjust = 2),axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),axis.text = element_text(size=16))+
        geom_text_repel(label=volcanoMatLabel,colour="black",
                  size=4,na.rm = TRUE,nudge_x = 0.2,)

print(volcanoplot)
pdf(paste0(fileDirs[1],GSEGetted,"_volcan.pdf"),width=9,height=7)
print(volcanoplot)
dev.off()
#===============Written to the file======================
matP_yours<-matFC_P[matFC_P$P_value<p_yours,]
write.table(gseMat[rownames(gseMat)%in%matP_yours$geneSYMOL,],file = paste0(fileDirs[1],"matP_",p_yours,".txt"),col.names = T,row.names = F,quote = F)
matFC_yours<-matFC_P[matFC_P$foldChange>fc_yours | matFC_P$foldChange<1/fc_yours,]
write.table(gseMat[rownames(gseMat)%in%matFC_yours$geneSYMOL,],file = paste0(fileDirs[1],"matFC_",fc_yours,".txt"),col.names = T,row.names = F,quote = F)
matFC_P_yours<-matFC_P[(matFC_P$foldChange>fc_yours | matFC_P$foldChange<1/fc_yours)&matFC_P$P_value<p_yours,]
write.table(gseMat[rownames(gseMat)%in%matFC_P_yours$geneSYMOL,],file = paste0(fileDirs[1],"matFC_P_",p_yours,"_",fc_yours,".txt"),col.names = T,row.names = F,quote = F)
#GENE OF NUMBERS
dim(gseMat)
#DEGs genes
diffgenep_fc<-matFC_P[(matFC_P$foldChange>fc_yours | matFC_P$foldChange<1/fc_yours)&matFC_P$P_value<p_yours,]
dim(diffgenep_fc)
write.table(diffgenep_fc[,1],paste0(fileDirs[1],GSEGetted,"_diffgene_",p_yours,"_",fc_yours,".txt"),row.names = F,col.names=F,quote=F)
#raise gene
dim(matFC_P[(matFC_P$foldChange>fc_yours)&matFC_P$P_value<p_yours,])
write.table(matFC_P[(matFC_P$foldChange>fc_yours)&matFC_P$P_value<p_yours,"geneSYMOL"],paste0(fileDirs[1],GSEGetted,"_raisegene_",p_yours,"_",fc_yours,".txt"),row.names = F,col.names=F,quote=F)
#reduce gene
dim(matFC_P[(matFC_P$foldChange<1/fc_yours)&matFC_P$P_value<p_yours,])
write.table(matFC_P[(matFC_P$foldChange<1/fc_yours)&matFC_P$P_value<p_yours,"geneSYMOL"],paste0(fileDirs[1],GSEGetted,"_downgene_",p_yours,"_",fc_yours,".txt"),row.names = F,col.names=F,quote=F)
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
matDiffgenes<-as.vector(matFC_P_yours$geneSYMOL)
matDiffMat<-as.matrix(gseMat[rownames(gseMat)%in%matDiffgenes,])
ColSideColor<-sampleCategory
ColSideColor[ColSideColor==0]<-"blue"
ColSideColor[ColSideColor==1]<-"red"
#==================heatmap======================
checkHeatmap()
heatmapDiffGene<-volcanoMatLabel[volcanoMatLabel!=""]
volcanoMatP_FC<-matP_yours[matP_yours$geneSYMOL%in%heatmapDiffGene,]
write.table(volcanoMatP_FC,paste0(fileDirs[1],GSEGetted,"_top20FC.txt"),row.names = F,col.names=F,quote=F)

pdf(paste0(fileDirs[1],GSEGetted,"_blue_red2.pdf"),width = 8, height = 8)
  heatmap.2(matDiffMat, col =greenred(100), scale = "column",
            labRow=F,key=TRUE, symkey=FALSE, density.info="none", 
            trace="none", cexRow=0.5,ColSideColors = ColSideColor
            ,main = paste0(GSEGetted," heatPlot (PDR<",p_yours,"and FC >",fc_yours,")"))
dev.off()
pdf(paste0(fileDirs[1],GSEGetted,"_green_orange.pdf"),width = 8, height = 8)
  annotation_ROW=data.frame(CellType=cbind(ColSideColor))
  rownames(annotation_ROW)<-colnames(matDiffMat)
  pheatmap(matDiffMat,annotation_col=annotation_ROW,show_rownames = F,fontsize = 14
           ,main = paste0(GSEGetted," heatPlot (PDR<",p_yours,"and FC >",fc_yours,")"))
dev.off()

#================GO Enrichment============================
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
fileDirs[1]<-paste0("../result/",GSEGetted,"/")
outLabel<-"BP"
#One of "BP", "MF", and "CC" subontologies, or "ALL" for all three
checkGO()
checkOntDir(GSEGetted,outLabel)
k<-keys(org.Hs.eg.db,keytype = "SYMBOL")
listSymolEntraz<-select(org.Hs.eg.db,key=k,columns=c("ENTREZID","SYMBOL","ENSEMBL"),keytype="SYMBOL")
ourGeneEntrazId<-listSymolEntraz[listSymolEntraz$SYMBOL%in%matFC_P_yours$geneSYMOL,]$ENTREZID
ego <- enrichGO(OrgDb="org.Hs.eg.db", gene = ourGeneEntrazId, ont = outLabel, pvalueCutoff = 0.05, readable= TRUE)
write.table(data.frame(ego),file=paste0(fileDirs[1],outLabel,"/","GO_term.txt"))
#===============KEGG Enrichment=======================
ekk <- enrichKEGG(gene= ourGeneEntrazId,organism  = 'hsa', qvalueCutoff = 0.05)
write.table(data.frame(ekk),file=paste0(fileDirs[1],outLabel,"/","KEGG_pathway.txt"))
#==============Written to the log============================
writeLogs(GSEGetted,GSEplatform,sampleCategory,p_yours,fc_yours)

BiocManager::install("e1071")
library("e1071")
pTable<-gseMat
y<-sampleCategory
pTable<-pTable[rownames(pTable)%in%matFC_P_yours$geneSYMOL,]
pTable<-t(pTable)
#1 sick,-1 no sick 
featureRankedList = svmrfeFeatureRanking(pTable,y)
print(featureRankedList[1:10])
geneOrder<-data.frame(cbind(rank=featureRankedList,ID=colnames(pTable)))
geneOrder<-geneOrder[order(featureRankedList),]
geneOrder<-cbind(geneOrder,score=round(((1+length(geneOrder[,1]))-as.numeric(as.vector(geneOrder$rank)))/length(geneOrder[,1]),5))
write.table(geneOrder,paste0("../result/",GSEGetted,"/","SVMRFEgeneOrder.txt"),col.names = T,row.names = F,quote = F)

other_score<-read.table(paste0("../result/",GSEGetted,"/","string_interactions.tsv default node.txt"),header = T,sep="\t")
other_score<-data.frame(AverageShortestPathLength=other_score$AverageShortestPathLength,Degree=other_score$Degree,name=other_score$name)
lastScore<-merge(other_score,geneOrder[,c(2,3)],by.x="name",by.y="ID",all.x=F,all.y=F)
resultmat<-data.frame(lastScore,ONEscore=lastScore[,3]*lastScore[,4]/lastScore[,2])
resultmat<-resultmat[order(resultmat$ONEscore,decreasing = T),]
write.table(resultmat,paste0("../result/",GSEGetted,"/","resultmat.txt"),col.names = T,row.names = F,quote = F)
