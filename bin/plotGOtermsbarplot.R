rm(list=ls())
source("preTreatFun.R",encoding="utf-8")
source("checkDenpe.R",encoding="utf-8")
source("getDiffFun.R",encoding="utf-8")
source("diffFun.R",encoding="utf-8")

matFC_P_yours<-read.table("../记录和报告/ALLnetworkGene.txt",header=F,stringsAsFactors =F)
matFC_P_yours<-data.frame(cbind(geneSYMOL=matFC_P_yours[,1]))
fileDirs<-c("../记录和报告/")
#One of "BP", "MF", and "CC" subontologies, or "ALL" for all three
checkGO()
k<-keys(org.Hs.eg.db,keytype = "SYMBOL")
listSymolEntraz<-select(org.Hs.eg.db,key=k,columns=c("ENTREZID","SYMBOL","ENSEMBL"),keytype="SYMBOL")
ourGeneEntrazId<-listSymolEntraz[listSymolEntraz$SYMBOL%in%matFC_P_yours$geneSYMOL,]$ENTREZID
#BP分析
egoBP<- enrichGO(OrgDb="org.Hs.eg.db", gene = ourGeneEntrazId, ont = "BP", 
                pvalueCutoff = 0.05, readable= TRUE) #GO富集分析
write.table(data.frame(egoBP),"BPGO.txt")
egoBPData<-data.frame(egoBP)
#MF分析
egoMF <- enrichGO(OrgDb="org.Hs.eg.db", gene = ourGeneEntrazId, ont = "MF", 
                pvalueCutoff = 0.05, readable= TRUE) #GO富集分析
egoMFData<-data.frame(egoMF)
#CC分析
egoCC <- enrichGO(OrgDb="org.Hs.eg.db", gene = ourGeneEntrazId, ont = "CC", 
                pvalueCutoff = 0.05, readable= TRUE) #GO富集分析
egoCCData<-data.frame(egoCC)
#选择前10个
display_number<-c(rep(10,3))
ego_result_BP<-egoBPData[seq(1,display_number[1]),]
ego_result_CC<-egoCCData[seq(1,display_number[2]),]
ego_result_MF<-egoMFData[seq(1,display_number[3]),]


go_enrich_df <- data.frame(ID=c(ego_result_BP$ID, ego_result_CC$ID, ego_result_MF$ID),
                           Description=c(ego_result_BP$Description, ego_result_CC$Description,
                           ego_result_MF$Description),GeneNumber=c(ego_result_BP$Count, 
                           ego_result_CC$Count, ego_result_MF$Count),
                           FDR=signif(as.numeric(as.vector(c(ego_result_BP$p.adjust, ego_result_CC$p.adjust,ego_result_MF$p.adjust))),3),
                           type=factor(c(rep("biological process", display_number[1]), 
                           rep("cellular component", display_number[2]),
                           rep("molecular function", display_number[3])), 
                           levels=c("molecular function", "cellular component", "biological process")))

## shorten the names of GO terms
shorten_names <- function(x, n_word=6, n_char=40){
  
  if (length(strsplit(x, " ")[[1]]) > n_word || (nchar(x) > 40))
  {
    if (nchar(x) > 40) x <- substr(x, 1, 40)
    x <- paste(paste(strsplit(x, " ")[[1]][1:min(length(strsplit(x," ")[[1]]), n_word)],
                     collapse=" "), "...", sep="")
    return(x)
  } 
  else
  {
    return(x)
  }
}

labels=as.vector(sapply(
  as.character(go_enrich_df$Description[as.numeric(go_enrich_df$Description)]),
  shorten_names))
names(labels)<-as.factor(rev(seq(1,length(labels))))

CPCOLS <- c("#8DA1CB", "#FD8D62", "#66C3A5")
library(ggplot2)
#========================================================
#以FDR值为x轴，标记gene数量为标签
p <- ggplot(data=go_enrich_df, aes(x=factor(Description,level=rev(Description)), y=-log10(FDR), fill=type)) +
  geom_bar(stat="identity", width=0.8) + coord_flip() +
  scale_fill_manual(values = CPCOLS) + theme_bw() + 
  scale_x_discrete(labels=labels) +
  xlab("GO term") +ylim(0,max(-log10(go_enrich_df$FDR))+5)+
  theme(axis.text=element_text(face = "bold", color="gray40"),
        axis.text.y = element_text(face = "bold", color="gray40",size = 12),
        axis.text.x = element_text(face = "bold", color="gray40",size = 12),
        axis.title.x = element_text(face = "bold",size = 14),
        axis.title.y = element_text(face = "bold",size = 14),
        axis.title  = element_text(face = "bold",size = 14),
        title = element_text(size=14,face="bold"),
        legend.text =element_text(face = "bold",size = 10),
        legend.title=element_blank()) +
  geom_text(aes(label=go_enrich_df$GeneNumber),hjust=0)
#========================================================
#以gene数量为x轴，标记FDR值为标签
p <- ggplot(data=go_enrich_df, aes(x=factor(Description,level=rev(Description)), y=GeneNumber, fill=type)) +
  geom_bar(stat="identity", width=0.8) + coord_flip() +
  scale_fill_manual(values = CPCOLS) + theme_bw() + 
  scale_x_discrete(labels=labels) +
  xlab("GO term") +ylim(0,max(go_enrich_df$GeneNumber)+5)+
  theme(axis.text=element_text(face = "bold", color="gray40"),
        axis.text.y = element_text(face = "bold", color="gray40",size = 12),
        axis.text.x = element_text(face = "bold", color="gray40",size = 12),
        axis.title.x = element_text(face = "bold",size = 14),
        axis.title.y = element_text(face = "bold",size = 14),
        axis.title  = element_text(face = "bold",size = 14),
        title = element_text(size=14,face="bold"),
        legend.text =element_text(face = "bold",size = 10),
        legend.title=element_blank()) +
  geom_text(aes(label=go_enrich_df$FDR),hjust=0)
#========================================================
#以gene数量为x轴，标记FDR值为标签,排序从大到小
go_enrich_df<-go_enrich_df[order(go_enrich_df$type,go_enrich_df$GeneNumber,decreasing = T),]
p <- ggplot(data=go_enrich_df, aes(x=factor(Description,level=rev(Description)), y=GeneNumber, fill=type)) +
  geom_bar(stat="identity", width=0.8) + coord_flip() +
  scale_fill_manual(values = CPCOLS) + theme_bw() + 
  scale_x_discrete(labels=labels) +
  xlab("GO term") +ylim(0,max(go_enrich_df$GeneNumber)+5)+
  theme(axis.text=element_text(face = "bold", color="gray40"),
        axis.text.y = element_text(face = "bold", color="gray40",size = 12),
        axis.text.x = element_text(face = "bold", color="gray40",size = 12),
        axis.title.x = element_text(face = "bold",size = 14),
        axis.title.y = element_text(face = "bold",size = 14),
        axis.title  = element_text(face = "bold",size = 14),
        title = element_text(size=14,face="bold"),
        legend.text =element_text(face = "bold",size = 10),
        legend.title=element_blank()) +
  geom_text(aes(label=go_enrich_df$FDR),hjust=0)
#======================================================
#绘制plot图
install.packages("GOplot")
library(GOplot)
#饼状图
#install.packages("ggpubr")
library(ggpubr)
pie1Data<-go_enrich_df[go_enrich_df$type=="molecular function",]
pie1Data<-pie1Data[order(pie1Data$GeneNumber),]
pie1 <- ggplot(data=pie1Data,aes(x = "", y = factor(GeneNumber,levels = seq(1,length(GeneNumber))), fill = Description)) + 
  geom_bar(stat = "identity", width = 1) +ylab("")+ xlab("")+ ## width >= 1 时中心的杂点将消失
  coord_polar(theta = "y")+ theme(axis.text.x = element_blank())+
  theme(axis.text=element_text(face = "bold", color="gray40"),
        axis.text.y = element_text(face = "bold", color="gray40",size = 12),
        axis.title.x = element_text(face = "bold",size = 14),
        axis.title.y = element_text(face = "bold",size = 14),
        axis.title  = element_text(face = "bold",size = 14),
        title = element_text(size=14,face="bold"),
        legend.text =element_text(face = "bold",size = 10),
        axis.ticks = element_blank(),axis.text.x = element_blank(),
        panel.grid=element_blank(),panel.border=element_blank())

pie2Data<-go_enrich_df[go_enrich_df$type=="cellular component",]
pie2Data<-pie2Data[order(pie2Data$GeneNumber),]
pie2 <- ggplot(data=pie2Data,aes(x = "", y = GeneNumber, fill = Description)) + 
  geom_bar(stat = "identity", width = 1) + ylab("")+ xlab("")+  ## width >= 1 时中心的杂点将消失
  coord_polar(theta = "y")+ theme(axis.text.x = element_blank())+
  theme(axis.text=element_text(face = "bold", color="gray40"),
        axis.text.y = element_text(face = "bold", color="gray40",size = 12),
        axis.title.x = element_text(face = "bold",size = 14),
        axis.title.y = element_text(face = "bold",size = 14),
        axis.title  = element_text(face = "bold",size = 14),
        title = element_text(size=14,face="bold"),
        legend.text =element_text(face = "bold",size = 10),
        axis.ticks = element_blank(),axis.text.x = element_blank(),
        panel.grid=element_blank(),panel.border=element_blank())

pie3Data<-go_enrich_df[go_enrich_df$type=="biological process",]
pie3Data<-pie3Data[order(pie3Data$GeneNumber),]
pie3 <- ggplot(data=pie3Data,aes(x = "", y = GeneNumber, fill = Description)) + 
  geom_bar(stat = "identity", width = 1) + ylab("")+ xlab("")+ ## width >= 1 时中心的杂点将消失
  coord_polar(theta = "y")+ theme(axis.text.x = element_blank())+
  theme(axis.text=element_text(face = "bold", color="gray40"),
        axis.text.y = element_text(face = "bold", color="gray40",size = 12),
        axis.title.x = element_text(face = "bold",size = 14),
        axis.title.y = element_text(face = "bold",size = 14),
        axis.title  = element_text(face = "bold",size = 14),
        title = element_text(size=14,face="bold"),
        legend.text =element_text(face = "bold",size = 10),
        axis.ticks = element_blank(),axis.text.x = element_blank(),
        panel.grid=element_blank(),panel.border=element_blank())

ggarrange(pie1,pie2,pie3,ncol=1,nrow=3,labels=c("MF","CC","BP"),align = "hv")

