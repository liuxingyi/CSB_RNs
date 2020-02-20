install.packages("VennDiagram")
library(VennDiagram)
library(ggpubr)
#=================绘制火山???====================
checkvolcano()
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
matFC_P<-matFC_P[matFC_P$P_value!=0,]
volcanoThresholdGSE71989<-getYouUpDown(matFC_P,fc_yours,p_yours)
windowsFonts(TNR = windowsFont("Times New Roman"))
#testMat<-cbind(volcanoThreshold,matFC_P)
#table(testMat[testMat$volcanoThreshold=="up",]$P_value>0.05)
volcanoXLog2FCGSE71989<-log2(matFC_P$foldChange)
volcanoYNegalog10PGSE71989<--log10(matFC_P$P_value)
volcanoMatGSE71989<-data.frame(cbind(geneSYMOL=as.vector(matFC_P$geneSYMOL),Threshold=volcanoThresholdGSE71989,
                             XLog2FC=volcanoXLog2FCGSE71989,
                             YNegalog10P=volcanoYNegalog10PGSE71989))

diffGeneGSE71989<-volcanoMatGSE71989[volcanoMatGSE71989$Threshold=="up" | volcanoMatGSE71989$Threshold=="down","geneSYMOL"]

volcanoMatLabelGSE71989<-labelHotspot(volcanoMatGSE71989,10)
volcanoplotGSE71989<-ggplot(data=volcanoMatGSE71989,aes(x=volcanoXLog2FCGSE71989,
                                        y=volcanoYNegalog10PGSE71989,colour=volcanoThresholdGSE71989),family="TNR")+
  geom_point(size=2.0,na.rm = TRUE)+
  scale_color_manual(values=c("#2f5688","#BBBBBB","#CC0000"))+
  geom_vline(xintercept = c(-log2(fc_yours),log2(fc_yours)),lty=4,col="black",lwd=0.5)+
  geom_hline(yintercept=-log10(p_yours),lty=4,col="black",lwd=0.5)+
  labs(x="log2(fold change)",y="-log10(FDR)",title=paste0("GSE71989 Volcano Plot (PDR<",p_yours,"and FC >",fc_yours,")"))+
  xlim(min(volcanoXLog2FCGSE71989)-1,max(volcanoXLog2FCGSE71989)+1)+
  theme(plot.title=element_text(hjust = 0.5,vjust = 2),axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),axis.text = element_text(size=16),
        legend.position='none')+
        geom_text_repel(label=volcanoMatLabelGSE71989,colour="black",
                  size=4,na.rm = TRUE,nudge_x = 0.2,)
print(volcanoplotGSE71989)
#====================================================================
checkvolcano()
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
matFC_P<-matFC_P[matFC_P$P_value!=0,]
volcanoThresholdGSE28735<-getYouUpDown(matFC_P,fc_yours,p_yours)
windowsFonts(TNR = windowsFont("Times New Roman"))
#testMat<-cbind(volcanoThreshold,matFC_P)
#table(testMat[testMat$volcanoThreshold=="up",]$P_value>0.05)
volcanoXLog2FCGSE28735<-log2(matFC_P$foldChange)
volcanoYNegalog10PGSE28735<--log10(matFC_P$P_value)
volcanoMatGSE28735<-data.frame(cbind(geneSYMOL=as.vector(matFC_P$geneSYMOL),Threshold=volcanoThresholdGSE28735,
                                     XLog2FC=volcanoXLog2FCGSE28735,
                                     YNegalog10P=volcanoYNegalog10PGSE28735))
diffGeneGSE28735<-volcanoMatGSE28735[volcanoMatGSE28735$Threshold=="up" | volcanoMatGSE28735$Threshold=="down","geneSYMOL"]

volcanoMatLabelGSE28735<-labelHotspot(volcanoMatGSE28735,10)
volcanoplotGSE28735<-ggplot(data=volcanoMatGSE28735,aes(x=volcanoXLog2FCGSE28735,
                                                        y=volcanoYNegalog10PGSE28735,colour=volcanoThresholdGSE28735),family="TNR")+
  geom_point(size=2.0,na.rm = TRUE)+
  scale_color_manual(values=c("#2f5688","#BBBBBB","#CC0000"))+
  geom_vline(xintercept = c(-log2(fc_yours),log2(fc_yours)),lty=4,col="black",lwd=0.5)+
  geom_hline(yintercept=-log10(p_yours),lty=4,col="black",lwd=0.5)+
  labs(x="log2(fold change)",y="",title=paste0("GSE28735 Volcano Plot (PDR<",p_yours,"and FC >",fc_yours,")"))+
  xlim(min(volcanoXLog2FCGSE28735)-1,max(volcanoXLog2FCGSE28735)+1)+
  theme(plot.title=element_text(hjust = 0.5,vjust = 2),axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),axis.text = element_text(size=16))+
  geom_text_repel(label=volcanoMatLabelGSE28735,colour="black",
                  size=4,na.rm = TRUE,nudge_x = 0.2,)
print(volcanoplotGSE28735)
#=====================================================
checkvolcano()
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
matFC_P<-matFC_P[matFC_P$P_value!=0,]
volcanoThresholdGSE15471<-getYouUpDown(matFC_P,fc_yours,p_yours)
windowsFonts(TNR = windowsFont("Times New Roman"))
#testMat<-cbind(volcanoThreshold,matFC_P)
#table(testMat[testMat$volcanoThreshold=="up",]$P_value>0.05)
volcanoXLog2FCGSE15471<-log2(matFC_P$foldChange)
volcanoYNegalog10PGSE15471<--log10(matFC_P$P_value)
volcanoMatGSE15471<-data.frame(cbind(geneSYMOL=as.vector(matFC_P$geneSYMOL),Threshold=volcanoThresholdGSE15471,
                                     XLog2FC=volcanoXLog2FCGSE15471,
                                     YNegalog10P=volcanoYNegalog10PGSE15471))

diffGeneGSE15471<-volcanoMatGSE15471[volcanoMatGSE15471$Threshold=="up" | volcanoMatGSE15471$Threshold=="down","geneSYMOL"]

classification<-volcanoThresholdGSE15471
volcanoMatLabelGSE15471<-labelHotspot(volcanoMatGSE15471,10)
volcanoplotGSE15471<-ggplot(data=volcanoMatGSE15471,aes(x=volcanoXLog2FCGSE15471,
                                                        y=volcanoYNegalog10PGSE15471,colour=classification),family="TNR")+
  geom_point(size=2.0,na.rm = TRUE)+
  scale_color_manual(values=c("#2f5688","#BBBBBB","#CC0000"))+
  geom_vline(xintercept = c(-log2(fc_yours),log2(fc_yours)),lty=4,col="black",lwd=0.5)+
  geom_hline(yintercept=-log10(p_yours),lty=4,col="black",lwd=0.5)+
  labs(x="",y="-log10(FDR)",title=paste0(GSEGetted," Volcano Plot (PDR<",p_yours,"and FC >",fc_yours,")"))+
  xlim(min(volcanoXLog2FCGSE15471)-1,max(volcanoXLog2FCGSE15471)+1)+
  theme(plot.title=element_text(hjust = 0.5,vjust = 2),axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),axis.text = element_text(size=16),
        legend.position='none')+
        geom_text_repel(label=volcanoMatLabelGSE15471,colour="black",
                  size=4,na.rm = TRUE,nudge_x = 0.2,)
print(volcanoplotGSE15471)

T<-venn.diagram(list(GSE15471=volcanoplotGSE15471,GSE71989=volcanoplotGSE71989,
                GSE28735=volcanoplotGSE28735),filename=NULL
                ,lwd=1,lty=2,col=c('red','green','blue')
                ,fill=c('red','green','blue')
                ,cat.col=c('red','green','blue')
                ,reverse=TRUE)
grid.draw(T)
venn(list(GSE15471=volcanoplotGSE15471,GSE71989=volcanoplotGSE71989,
          GSE28735=volcanoplotGSE28735))
dev.off()
ggarrange(volcanoplotGSE15471,volcanoplotGSE28735,volcanoplotGSE71989,
          nrow=2,ncol=2,heights = c(1,1,1,1),widths=c(1,1,1,1),common.legend = TRUE,hjust=-2)

