#' Draw a bar chart for GO enrichment analysis
#'
#' @description According to the input geneSymol collection,
#'  convert to entrazid and draw three types of bar charts.
#' @param geneSet as.vector,gene dataset.
#' @import Rgraphviz
#' @importFrom clusterProfiler enrichGO
#' @export plotGOtermsbarplot
#' @author Xingyi Liu

plotGOtermsbarplot <- function(geneSet,
                       label.min.n.words = 10,
                       label.min.n.char = 190,
                       display.number = 10,
                       plot.type = 2,
                       enrich.pvalue = 0.05,
                       xlab = "GO term",
                       ylab = "-log10(FDR)",
                       title = "GO enrichment plot",
                       CPCOLS = c("#2f5688","#0F1F95","#F2D06D","#CC0000"),
                       face = "bold",
                       axis.text.color = "black",
                       axis.text.font.size = 14,
                       axis.title.font.size = 14,
                       legend.text.font.size = 14,
                       prefixLable="enrichment"){
    
    if (!requireNamespace("org.Hs.eg.db")) {
        BiocManager::install("org.Hs.eg.db")
    }
    suppressMessages(library("org.Hs.eg.db"))
    k <- keys(org.Hs.eg.db,keytype = "SYMBOL")
    listSymolEntraz <- bitr(k, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db");
    ourGeneEntrazId<-listSymolEntraz[listSymolEntraz$SYMBOL%in%geneSet,]$ENTREZID
    OrgDb="org.Hs.eg.db"
    egoBP<- enrichGO(OrgDb = OrgDb,
                     gene = ourGeneEntrazId,
                     ont = "BP",
                     pvalueCutoff =  enrich.pvalue,
                     readable= TRUE)
    egoBPData<-data.frame(egoBP)
    print("molecular function enrichment finished!!")
    egoMF <- enrichGO(OrgDb = OrgDb,
                      gene = ourGeneEntrazId,
                      ont = "MF",
                      pvalueCutoff =  enrich.pvalue,
                      readable= TRUE)
    egoMFData<-data.frame(egoMF)
    print("cellular component enrichment finished!!")
    egoCC <- enrichGO(OrgDb = OrgDb,
                      gene = ourGeneEntrazId,
                      ont = "CC",
                      pvalueCutoff =  enrich.pvalue,
                      readable= TRUE)
    egoCCData<-data.frame(egoCC)
    print("biological process enrichment finished!!")
    ekk <- enrichKEGG(gene= ourGeneEntrazId,
                      organism  = "hsa",
                      pvalueCutoff = enrich.pvalue)
    KEGGData <- data.frame(ekk)
    print("kegg pathway enrichment finished!!")
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
                                                    levels=c("biological process", "cellular component", "molecular function","kegg pathway")),
                                        color=factor(c(rep(CPCOLS[1],display_number[1]), rep(CPCOLS[2], display_number[2]),rep(CPCOLS[3], display_number[3]),rep(CPCOLS[4],display_number[4])),
                                                          levels=CPCOLS))
    print(paste0("BP CC MF kegg were writed to",prefixLable,".txt"))
    write.table(rbind(ego_result_BP,ego_result_CC,ego_result_MF,ego_result_KEGG),file=paste0(prefixLable,".txt"),col.names = F,row.names = F,quote = F)
    #delete all NA lines
    go_enrich_df <- go_enrich_df_writeted[!is.na(go_enrich_df_writeted$Description),]
    #shorten GO or Kegg terms.
    labelSet <- as.character(go_enrich_df$Description)
    newlabelSet <- NULL
    for(i in seq(1,length(labelSet))){
        newlabelSet <- c(newlabelSet,shorten_names(labelSet[i],n_word=label.min.n.words, n_char=label.min.n.char))
    }
    go_enrich_df$Description <- as.character(newlabelSet)
    go_enrich_df$color <- as.factor(as.character(go_enrich_df$color))
    
    if(plot.type == 1)
    {
        go_enrich_df_new <- NULL
        for(i in CPCOLS){
            go_enrich_df_line <- go_enrich_df[go_enrich_df$color==i,]
            go_enrich_df_line <-go_enrich_df_line[order(go_enrich_df_line$FDR,decreasing = F),]
            go_enrich_df_new <- rbind(go_enrich_df_new,go_enrich_df_line)
        }
        go_enrich_df <- go_enrich_df_new
        levels(go_enrich_df$color) <- CPCOLS[CPCOLS%in%levels(go_enrich_df$color)]
        p <- ggplot(data=go_enrich_df, aes(x=factor(Description,level=rev(Description)), y=-log10(FDR), fill=type)) +
            geom_bar(stat="identity", width=0.8) +
            scale_fill_manual(values = levels(go_enrich_df$color)) + theme_bw() +
            scale_x_discrete(labels=rev(go_enrich_df$Description)) +
            xlab(xlab) + ylab(ylab) +ylim(0,max(-log10(go_enrich_df$FDR))*1.2)+labs(title=title)+
            theme(axis.text=element_text(face = face, color = axis.text.color),
                  axis.text.y = element_text(face = face, color=axis.text.color,size = axis.text.font.size),
                  axis.text.x = element_text(face = face, color=axis.text.color,size = axis.text.font.size),
                  axis.title.x = element_text(face = face,size = axis.title.font.size),
                  axis.title.y = element_text(face = face,size = axis.title.font.size),
                  axis.title  = element_text(face = face,size = axis.title.font.size),
                  title = element_text(size= axis.title.font.size,face=face),
                  legend.text =element_text(face = face,size = legend.text.font.size),
                  legend.title=element_blank()) +
            geom_text(aes(label=GeneNumber),hjust=0) +
            coord_flip()
    }else if(plot.type == 2){
        go_enrich_df_new <- NULL
        for(i in CPCOLS){
            go_enrich_df_line <- go_enrich_df[go_enrich_df$color==i,]
            go_enrich_df_line <-go_enrich_df_line[order(go_enrich_df_line$GeneNumber,decreasing = T),]
            go_enrich_df_new <- rbind(go_enrich_df_new,go_enrich_df_line)
        }
        go_enrich_df <- go_enrich_df_new
        levels(go_enrich_df$color) <- CPCOLS[CPCOLS%in%levels(go_enrich_df$color)]
        p <- ggplot(data=go_enrich_df, aes(x=factor(Description,level=rev(Description)), y=GeneNumber, fill=type)) +
            geom_bar(stat="identity", width=0.8)  +
            scale_fill_manual(values = levels(go_enrich_df$color)) + theme_bw() +
            scale_x_discrete(labels=rev(go_enrich_df$Description)) +
            xlab(xlab) + ylab(ylab) +ylim(0,max(go_enrich_df$GeneNumber)*1.15)+labs(title=title)+
            theme(axis.text=element_text(face = face, color = axis.text.color),
                  axis.text.y = element_text(face = face, color= axis.text.color,size = axis.text.font.size),
                  axis.text.x = element_text(face = face, color= axis.text.color,size = axis.text.font.size),
                  axis.title.x = element_text(face = face,size = axis.title.font.size),
                  axis.title.y = element_text(face = face,size = axis.title.font.size),
                  axis.title  = element_text(face = face,size = axis.title.font.size),
                  title = element_text(size=14,face= face),
                  legend.text =element_text(face = face,size = legend.text.font.size),
                  legend.title=element_blank()) +
            geom_text(aes(label=FDR),hjust=0) +
            coord_flip()
    }else {
        stop("plot.type Error!!!")
    }
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

prexShort <- function(x){
    x <- as.character(go_enrich_df$type)
    x[x == "biological process"] = "BP"
    x[x == "cellular component"] = "CC"
    x[x == "molecular function"] = "MF"
    x[x == "kegg pathway"] = "KEGG"
    return(x)
}
