# CSB_RNs
Procedure flow for drug target discovery
## example codes
```
profilerFrame <- get.GEO.data(GEO.acession = "28735",platfrom = "GPL6244")
boxplot(profilerFrame,col="gray",xlab="samples",ylab="express label",main="expression quentity boxplt")
profilerFramecleaned<-geneentropyfilter(profilerFrame)
profileData <- analysisIt(profilerFramecleaned,sample.label = c(rep(c(1,0),times=45)))
volcanoPlot(profileData)
geneSet <- rownames(profileData@data)[profileData@selected.gene]
plotGOtermsbarplot(geneSet = geneSet)
```
