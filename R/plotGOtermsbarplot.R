#读入mRNA、miRNA的数据表达谱
library("edgeR")
library("targetFinder")
library(AnnotationHub)	#library导入需要使用的数据包
library(org.Hs.eg.db)   #人类注释数据库
library(clusterProfiler)
library(dplyr)
library(ggplot2)
data <- read.csv("mRNACount.csv",header = T,row.names = 1)
miRNAdata <- read.csv("miRNACount.csv",header = T,row.names = 1)
lncRNAdata <- read.csv("lncRNACount.csv",header = T,row.names = 1)
data <- round(data)
lncRNAdata <- round(lncRNAdata)
dataCleaned <- data[!rowSums(data)==0,]
miRNAdata <- miRNAdata[!rowSums(miRNAdata)==0,]
lncRNAdata <- lncRNAdata[!rowSums(lncRNAdata)==0,]
#读入样本注释信息
clinical <- read.table("clinacal.txt",header = T,row.names = 1)
rownames(clinical) <- gsub("-",".",rownames(clinical))
group <- as.factor(clinical[,1])
dataCleaned <- dataCleaned[,match(rownames(clinical),substr(gsub("\\.","-",colnames(dataCleaned)),start = 1,stop = nchar(colnames(dataCleaned))-6))]
colnames(dataCleaned) <- substr(colnames(dataCleaned),start = 1,stop = nchar(colnames(dataCleaned))-6)
miRNAdata <- miRNAdata[,match(rownames(clinical),substr(gsub("\\.","-",colnames(miRNAdata)),start = 1,stop = nchar(colnames(miRNAdata))-6))]
#starBase数据库miRNA和mRNA
starBaseMiRnaMRNA <- read.table("ENCORI_hg19_CLIP-seq_all_allGene.txt",skip = 3,header=T)
starBaselncRnaMRNA <- read.table("ENCORI_hg19_CLIP-seq_all_alllncRNA.txt",skip = 3,header=T,sep = "\t")
starBaseMiRnaMRNA <- starBaseMiRnaMRNA[rowSums(starBaseMiRnaMRNA[,seq(15,21)])>4 & starBaseMiRnaMRNA[,22]>6 &starBaseMiRnaMRNA[,12]>6,]
starBaselncRnaMRNA <- starBaselncRnaMRNA[starBaselncRnaMRNA[,17]>6 & starBaselncRnaMRNA[,10]>6,]
#额外的编号对应信息
emsembleAnno <- read.table("./hg19.fa",fill = T,sep = ";")
ENSGENST <- cbind(substr(as.character(emsembleAnno[,1]),start = 9,stop = 23),substr(as.character(emsembleAnno[,2]),start =16,stop = 33))

#差异分析
edgeRT <- function(x,group,min.count = 10, min.total.count = 15){
    y <- DGEList(counts=x,group=group)
    keep <- filterByExpr(y,min.count = 10, min.total.count = 15)
    y <- y[keep,,keep.lib.sizes=FALSE]
    y <- calcNormFactors(y)
    design <- model.matrix(~group)
    y <- estimateDisp(y,design)
    fit <- glmQLFit(y,design)
    qlf <- glmQLFTest(fit,coef=2)
    qlfTable <- cbind(qlf$table,FDR.BH=p.adjust(qlf$table[,4],method = "BH"))
    return(qlfTable)
}
dataCleaned <- dataCleaned[,match(colnames(dataCleaned),rownames(clinical))]
question1mRNA <- dataCleaned[,c("N.He.IR.1","N.He.IR.2","N.He.IR.3","N.He.Ctrl.1","N.He.Ctrl.2","N.He.Ctrl.3")]
question1mRNAclinical <- clinical[c("N.He.IR.1","N.He.IR.2","N.He.IR.3","N.He.Ctrl.1","N.He.Ctrl.2","N.He.Ctrl.3"),]
question1mRNAGroup <- as.factor(question1mRNAclinical[,2])
mRNAqlfTable <- edgeRT(x=question1mRNA,group=question1mRNAGroup)
mRNAqlfTableDiff <-mRNAqlfTable[mRNAqlfTable[,4]<0.05 &(mRNAqlfTable[,1]>log(3/2) | mRNAqlfTable[,1]<log(2/3)),]
k<-keys(org.Hs.eg.db,keytype = "ENSEMBL")
listSymolEntraz <- bitr(k, fromType="ENSEMBL", toType=c("ENTREZID","SYMBOL"), OrgDb="org.Hs.eg.db");
ourGeneEntrazId<-listSymolEntraz[listSymolEntraz$ENSEMBL%in%rownames(mRNAqlfTableDiff),]$SYMBOL
plotGOtermsbarplot(ourGeneEntrazId)

miRNAqlfTable <- edgeRT(x=miRNAdata,group=group)
miRNAqlfTableDiff <-miRNAqlfTable[(miRNAqlfTable[,5]<0.01),]
lncRNAqlfTable <- edgeRT(x=lncRNAdata,group=group)
lncRNAqlfTableDiff <-lncRNAqlfTable[(lncRNAqlfTable[,5]<0.01),]


unique(rownames(miRNAqlfTableDiff))
rownames(mRNAqlfTableDiff)
lncRNAqlfTableDiffSelected  <- unique(ENSGENST[ENSGENST[,2] %in% rownames(lncRNAqlfTableDiff),1])
#获取miRNA和mRNA网络
ourMiM <- NULL
for(i in seq(1,nrow(starBaseMiRnaMRNA))){
    mRNA.tmp <- as.character(starBaseMiRnaMRNA[i,2]) %in% rownames(miRNAqlfTableDiff)
    miRNA.tmp <- as.character(starBaseMiRnaMRNA[i,3]) %in% rownames(mRNAqlfTableDiff)
    if(!(mRNA.tmp & miRNA.tmp))next
    ourMiM <- rbind(ourMiM,starBaseMiRnaMRNA[i,])
    print(paste0(i,"|",nrow(starBaseMiRnaMRNA)))
}
ourLnM <- NULL
for(i in seq(1,nrow(starBaselncRnaMRNA))){
    lncRNA.tmp <-  as.character(starBaselncRnaMRNA[i,3]) %in% lncRNAqlfTableDiffSelected
    miRNA.tmp <- as.character(starBaselncRnaMRNA[i,2]) %in% rownames(miRNAqlfTableDiff)
    if(!(mRNA.tmp & miRNA.tmp))next
    ourLnM <- rbind(ourLnM,starBaselncRnaMRNA[i,])
    print(paste0(i,"|",nrow(starBaselncRnaMRNA)))
}
ResTable <- merge(ourMiM,ourLnM,by="miRNAname")
ResTable <- ResTable[,c("miRNAname","geneID.x","geneName.x","geneType.x","geneID.y","geneName.y","geneType.y")]
ResTable <- unique(ResTable[ResTable[,7]=="lincRNA",])

unique(paste0(ResTable$miRNAname,ResTable$geneID.x,ResTable$geneName.x,ResTable$geneType.x,ResTable$geneID.y,ResTable$geneName.y,ResTable$geneType.y))
write.table(as.character(unique(ResTable$geneName.x)),file="ceRNAgene.txt",quote = F,row.names = F,col.names = F)
oxygen_module <- read.table("oxygen_module.txt")
write.table(as.character(oxygen_module[oxygen_module[,3]=="turquoise",2]),file="turquoise_gene.txt",quote = F,row.names = F,col.names = F)
write.table(rbind(cbind(as.character(ResTable$miRNAname),geneSymbol=as.character(ResTable$geneID.x),as.character(ResTable$geneType.x)),cbind(as.character(ResTable$miRNAname),geneSymbol=as.character(ResTable$geneID.y),as.character(ResTable$geneType.y))),file="ceRNAnet.txt",quote = F,row.names = F,col.names = F)

APPSelected <-ResTable[ResTable$geneName.x=="APP" | ResTable$geneName.x=="FOXO3" | ResTable$geneName.x=="FBN1",]
APPSelectedW <- rbind(cbind(as.character(APPSelected$miRNAname),geneSymbol=as.character(APPSelected$geneName.x),as.character(APPSelected$geneType.x)),cbind(as.character(APPSelected$miRNAname),geneSymbol=as.character(APPSelected$geneName.y),as.character(APPSelected$geneType.y)))
APPSelectedWV <- strsplit(unique(paste0(APPSelectedW[,1],":",APPSelectedW[,2],":",APPSelectedW[,3])),split = ":")
APPSelectedWVC <- NULL
for(i in seq(1,length(APPSelectedWV))){
    APPSelectedWVC <- rbind(APPSelectedWVC,APPSelectedWV[[i]])
}
write.table(APPSelectedWVC,file="APPceRNAnet.txt",quote = F,row.names = F,col.names = F)


gene159 <-read.table("107.txt",header=F)
plotGOtermsbarplot(ourGeneEntrazId)
k<-keys(org.Hs.eg.db,keytype = "SYMBOL")
listSymolEntraz <- bitr(k, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db");
ourGeneEntrazId<-listSymolEntraz[listSymolEntraz$ENSEMBL%in%rownames(mRNAqlfTableDiff),]$SYMBOL
library(clusterProfiler)
ekk <- enrichKEGG(gene= ourGeneEntrazId,organism  = 'hsa', pvalueCutoff = 0.05)
dotplot(ekk,font.size=8)
ego <- enrichGO(OrgDb="org.Hs.eg.db", gene = ourGeneEntrazId, ont = "BP", pvalueCutoff = 0.01, readable= TRUE) #GO富集分析
dotplot(ego,showCategory=10,title="Enrichment GO Top10")

BPdata <- read.table("BPenrichment.Process.tsv",header=F,fill =T,sep = "\t",quote = "")
BPdataV <- NULL
for(i in seq(1,10)){
    BPdataV <- c(BPdataV,unlist(strsplit(as.character(BPdata[i,7]),split = "[,]")))
}
write.table(data.frame(table(BPdataV))[data.frame(table(BPdataV))[,2]>5,],file="BPdataV.txt",row.names = F,col.names = T,quote = F)
