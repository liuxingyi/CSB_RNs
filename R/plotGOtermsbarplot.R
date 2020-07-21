#' Draw a bar chart for GO enrichment analysis
#'
#' @description According to the input geneSymol collection,
#'  convert to entrazid and draw three types of bar charts.
#' @param geneSet as.vector,gene dataset.
#' @import Rgraphviz
#' @importFrom clusterProfiler enrichGO
#' @export plotGOtermsbarplot
#' @author Xingyi Liu


plotGOtermsbarplot <-function(geneSet,
                              label.min.n.words = 15,
                              label.min.n.char = 90,
                              display.number = 10,
                              Species.label="hsa",
                              enrich.pvalue = 0.05,
                              plot.type = 2,
                              coord_flip =TRUE,
                              CPCOLS = c("#2f5688","#0F1F95","#CC0000","red"),
                              xlab = "GO term",
                              ylab = "-log10(FDR)",
                              title = "GO enrichment plot",
                              face = "bold",
                              axis.text.color = "gray40",
                              axis.text.font.size = 12,
                              axis.title.font.size = 14,
                              legend.text.font.size = 10,
                              prefixLable="enrichment"){
    if (!requireNamespace("org.Hs.eg.db")) {
        BiocManager::install("org.Hs.eg.db")
    }
    suppressMessages(library("org.Hs.eg.db"))
    k <- keys(org.Hs.eg.db,keytype = "SYMBOL")
    listSymolEntraz <- bitr(k, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db");
    ourGeneEntrazId<-listSymolEntraz[listSymolEntraz$SYMBOL%in%geneSet,]$ENTREZID
    if(Species.label == "hsa") OrgDb="org.Hs.eg.db"
    egoBP<- enrichGO(OrgDb = OrgDb,
                     gene = ourGeneEntrazId,
                     ont = "BP",
                     pvalueCutoff =  enrich.pvalue,
                     readable= TRUE)
    egoBPData<-data.frame(egoBP)
    egoMF <- enrichGO(OrgDb = OrgDb,
                      gene = ourGeneEntrazId,
                      ont = "MF",
                      pvalueCutoff =  enrich.pvalue,
                      readable= TRUE)
    egoMFData<-data.frame(egoMF)
    egoCC <- enrichGO(OrgDb = OrgDb,
                      gene = ourGeneEntrazId,
                      ont = "CC",
                      pvalueCutoff =  enrich.pvalue,
                      readable= TRUE)
    egoCCData<-data.frame(egoCC)
    ekk <- enrichKEGG(gene= ourGeneEntrazId,
                      organism  = Species.label,
                      pvalueCutoff = enrich.pvalue)
    KEGGData <- data.frame(ekk)
    
    display_number<-c(rep(display.number,4))
    ego_result_BP<-egoBPData[seq(1,display_number[1]),]
    ego_result_CC<-egoCCData[seq(1,display_number[2]),]
    ego_result_MF<-egoMFData[seq(1,display_number[3]),]
    ego_result_KEGG<-KEGGData[seq(1,display_number[4]),]

    go_enrich_df_writeted <- data.frame(ID=c(ego_result_BP$ID, ego_result_CC$ID, ego_result_MF$ID,ego_result_KEGG$ID),
                               Description=c(ego_result_BP$Description, ego_result_CC$Description,ego_result_MF$Description,ego_result_KEGG$Description),
                               GeneNumber=c(ego_result_BP$Count, ego_result_CC$Count, ego_result_MF$Count,as.numeric(unlist(strsplit(ego_result_KEGG$GeneRatio,split = "[//]"))[seq(1,2*display.number,2)])),
                               FDR=signif(as.numeric(as.vector(c(ego_result_BP$p.adjust, ego_result_CC$p.adjust,ego_result_MF$p.adjust,ego_result_KEGG$p.adjust))),3),
                               type=factor(c(rep("biological process", display_number[1]), rep("cellular component", display_number[2]),rep("molecular function", display_number[3]),rep("kegg pathway",display_number[4])),
                                           levels=c("molecular function", "cellular component", "biological process","kegg pathway")))
    write.table(go_enrich_df_writeted,file=paste0(prefixLable,".txt"),col.names = F,row.names = F,quote = F)
    go_enrich_df<-go_enrich_df_writeted

    if(plot.type == 1)
    {
      
      go_enrich_df <- rbind(go_enrich_df[seq(1,display_number[1]),][order(go_enrich_df[seq(1,display_number[1]),]$FDR,decreasing = F),],
                            go_enrich_df[seq(display_number[1]+1,2*display_number[2]),][order(go_enrich_df[seq(display_number[1]+1,2*display_number[2]),]$FDR,decreasing = F),],
                            go_enrich_df[seq(2*display_number[1]+1,3*display_number[3]),][order(go_enrich_df[seq(2*display_number[1]+1,3*display_number[3]),]$FDR,decreasing = F),],
                            go_enrich_df[seq(3*display_number[1]+1,4*display_number[4]),][order(go_enrich_df[seq(3*display_number[1]+1,4*display_number[4]),]$FDR,decreasing = F),])
      labelSet <- as.character(go_enrich_df$Description)
      newlabelSet <- NULL
      for(i in seq(1,length(labelSet))){
        newlabelSet <- c(newlabelSet,shorten_names(labelSet[i],n_word=label.min.n.words, n_char=label.min.n.char))
      }
      labelSet <- newlabelSet
      labels <- as.vector(labelSet)
      go_enrich_df$Description <- labels
      go_enrich_df <- go_enrich_df[!is.na(go_enrich_df$GeneNumber),]
      labels <- go_enrich_df$Description
      #names(labels)<-as.factor(rev(seq(1,length(labels))))
      p <- ggplot(data=go_enrich_df, aes(x=factor(Description,level=rev(Description)), y=-log10(FDR), fill=type)) +
        geom_bar(stat="identity", width=0.8) +
        scale_fill_manual(values = CPCOLS) + theme_bw() +
        scale_x_discrete(labels=labels) +
        xlab(xlab) + ylab(ylab) +ylim(0,max(-log10(go_enrich_df$FDR))+10)+labs(title=title)+
        theme(axis.text=element_text(face = face, color = axis.text.color),
              axis.text.y = element_text(face = face, color=axis.text.color,size = axis.text.font.size),
              axis.text.x = element_text(face = face, color=axis.text.color,size = axis.text.font.size),
              axis.title.x = element_text(face = face,size = axis.title.font.size),
              axis.title.y = element_text(face = face,size = axis.title.font.size),
              axis.title  = element_text(face = face,size = axis.title.font.size),
              title = element_text(size= axis.title.font.size,face=face),
              legend.text =element_text(face = face,size = legend.text.font.size),
              legend.title=element_blank()) +
        geom_text(aes(label=GeneNumber),hjust=0)
    }else if(plot.type == 2){
      if(ylab == "-log10(FDR)")ylab="the number of genes"
      go_enrich_df <- rbind(go_enrich_df[seq(1,display_number[1]),][order(go_enrich_df[seq(1,display_number[1]),]$GeneNumber,decreasing = T),],
                            go_enrich_df[seq(display_number[1]+1,2*display_number[2]),][order(go_enrich_df[seq(display_number[1]+1,2*display_number[2]),]$GeneNumber,decreasing = T),],
                            go_enrich_df[seq(2*display_number[1]+1,3*display_number[3]),][order(go_enrich_df[seq(2*display_number[1]+1,3*display_number[3]),]$GeneNumber,decreasing = T),],
                            go_enrich_df[seq(3*display_number[1]+1,4*display_number[4]),][order(go_enrich_df[seq(3*display_number[1]+1,4*display_number[4]),]$GeneNumber,decreasing = T),])
      labelSet <- as.character(go_enrich_df$Description)
      newlabelSet <- NULL
      for(i in seq(1,length(labelSet))){
        newlabelSet <- c(newlabelSet,shorten_names(labelSet[i],n_word=label.min.n.words, n_char=label.min.n.char))
      }
      labelSet <- newlabelSet
      labels <- as.vector(labelSet)
      go_enrich_df$Description <- labels
      go_enrich_df <- go_enrich_df[!is.na(go_enrich_df$GeneNumber),]
      labels <- go_enrich_df$Description
      p <- ggplot(data=go_enrich_df, aes(x=factor(Description,level=rev(Description)), y=GeneNumber, fill=type)) +
        geom_bar(stat="identity", width=0.8)  +
        scale_fill_manual(values = CPCOLS) + theme_bw() +
        scale_x_discrete(labels=labels) +
        xlab(xlab) + ylab(ylab) +ylim(0,max(go_enrich_df$GeneNumber)+10)+labs(title=title)+
        theme(axis.text=element_text(face = face, color = axis.text.color),
              axis.text.y = element_text(face = face, color= axis.text.color,size = axis.text.font.size),
              axis.text.x = element_text(face = face, color= axis.text.color,size = axis.text.font.size),
              axis.title.x = element_text(face = face,size = axis.title.font.size),
              axis.title.y = element_text(face = face,size = axis.title.font.size),
              axis.title  = element_text(face = face,size = axis.title.font.size),
              title = element_text(size=14,face= face),
              legend.text =element_text(face = face,size = legend.text.font.size),
              legend.title=element_blank()) +
        geom_text(aes(label=FDR),hjust=0)
    }else {
      stop("plot.type Error!!!")
    }
    if(coord_flip ==TRUE) p<-p+coord_flip()
    return(p)
    }

shorten_names <- function(x, n_word=3, n_char=10){
    if(is.na(x))return(x)
    if (length(strsplit(x, " ")[[1]]) > n_word || (nchar(x) > n_char))
    {
        if (nchar(x) > n_char) x <- substr(x, 1, n_char)
        x <- paste(paste(strsplit(x, " ")[[1]][1:min(length(strsplit(x," ")[[1]]), n_word)],
                         collapse=" "), "...", sep="")
        return(x)
    }
    else
    {
        return(x)
    }
}
