### 1) EDIT THE FOLLOWING THREE LINES TO PROVIDE PATHS TO THE fineRADstructure OUTPUT 
setwd("~/Dropbox/Postdoc_Milan/LK_Joan/intersect_Tasos_Joan_vcfs/fineradstructure/") ## The directory where the files are located
chunkfile<-"populations.haps_noinbred_def_def_chunks.out" ## RADpainter output file
mcmcfile<-"populations.haps_noinbred_def_def_chunks.mcmc.xml" ## finestructure mcmc file
treefile<-"populations.haps_noinbred_def_def_chunks.mcmcTree.xml" ## finestructure tree file
### 2) EDIT THIS PATH TO WHERE YOU WANT THE PLOTS:
plotsFolder <- "~/Dropbox/Postdoc_Milan/LK_Joan/intersect_Tasos_Joan_vcfs/fineradstructure/"
### 3) SET VALUES FOR THESE VARIABLES: "analysisName" will be included in output plots
analysisName <- "stacks_unfiltered_noinbred_def_def";  maxIndv <- 10000; maxPop<-10000


### 4) EDIT THE PATH TO YOUR COPY of FinestructureLibrary.R
source("~/Dropbox/Tesi/Data/Durham/fineRADStructure/FinestructureLibrary.R", chdir = TRUE) # read in the R functions, which also calls the needed packages

### 5) EXECUTE THE CODE ABOVE AND THE REST OF THE CODE BELOW
## make some colours
some.colors<-MakeColorYRP() # these are yellow-red-purple
some.colorsEnd<-MakeColorYRP(final=c(0.2,0.2,0.2)) # as above, but with a dark grey final for capped values
library(scales)
#some.colors <- viridis_pal(option="D", direction = -1)(61)
#some.colorsEnd <- append(some.colors, "#333333")
#nb.cols <- 61
#some.colors <- colorRampPalette(brewer.pal(8, "YlOrRd"))(nb.cols)
#some.colorsEnd <- append(some.colors, "#333333")
###### READ IN THE CHUNKCOUNT FILE
dataraw<-as.matrix(read.table(chunkfile,row.names=1,header=T,skip=1)) # read in the pairwise coincidence 
###### READ IN THE MCMC FILES
mcmcxml<-xmlTreeParse(mcmcfile) ## read into xml format
mcmcdata<-as.data.frame.myres(mcmcxml) ## convert this into a data frame
###### READ IN THE TREE FILES
treexml<-xmlTreeParse(treefile) ## read the tree as xml format
ttree<-extractTree(treexml) ## extract the tree into ape's phylo format

## Reduce the amount of significant digits printed in the posteror assignment probabilities (numbers shown in the tree):
ttree$node.label[ttree$node.label!=""] <-format(as.numeric(ttree$node.label[ttree$node.label!=""]),digits=2)
# convert to dendrogram format
tdend<-myapetodend(ttree,factor=1)
## Now we work on the MAP state
mapstate<-extractValue(treexml,"Pop") # map state as a finestructure clustering
mapstatelist<-popAsList(mapstate) # .. and as a list of individuals in populations
popnames<-lapply(mapstatelist,NameSummary) # population names IN A REVERSIBLE FORMAT (I.E LOSSLESS)
## NOTE: if your population labels don't correspond to the format we used (NAME<number>) YOU MAY HAVE TROUBLE HERE. YOU MAY NEED TO RENAME THEM INTO THIS FORM AND DEFINE YOUR POPULATION NAMES IN popnamesplot BELOW
popnamesplot<-lapply(mapstatelist,NameMoreSummary) # a nicer summary of the populations
names(popnames)<-popnamesplot # for nicety only
names(popnamesplot)<-popnamesplot # for nicety only
popdend<-makemydend(tdend,mapstatelist) # use NameSummary to make popdend
popdend<-fixMidpointMembers(popdend) # needed for obscure dendrogram reasons
popdendclear<-makemydend(tdend,mapstatelist,"NameMoreSummary")# use NameMoreSummary to make popdend
popdendclear<-fixMidpointMembers(popdendclear) # needed for obscure dendrogram reasons


########################
## Plot 1: COANCESTRY MATRIX
fullorder<-labels(tdend) # the order according to the tree
datamatrix<-dataraw[fullorder,fullorder] # reorder the data matrix

tmpmat<-datamatrix 
tmpmat[tmpmat>maxIndv]<-maxIndv #  # cap the heatmap
pdf(file=paste(plotsFolder,analysisName,"-SimpleCoancestry.pdf",sep=""),height=25,width=25)
plotFinestructure(tmpmat,dimnames(tmpmat)[[1]],dend=tdend,cols=some.colorsEnd,cex.axis=1.1,edgePar=list(p.lwd=0,t.srt=90,t.off=-0.1,t.cex=1.2))
dev.off()

########################
## Plot 2: POPULATIONS AND COANCESTRY AVERAGES
popmeanmatrix<-getPopMeanMatrix(datamatrix,mapstatelist)

tmpmat<-popmeanmatrix
tmpmat[tmpmat>maxPop]<-maxPop # cap the heatmap
pdf(file=paste(plotsFolder,analysisName,"-PopAveragedCoancestry.pdf",sep=""),height=20,width=20)
plotFinestructure(tmpmat,dimnames(tmpmat)[[1]],dend=tdend,cols=some.colorsEnd,cex.axis=1.1,edgePar=list(p.lwd=0,t.srt=90,t.off=-0.1,t.cex=1.2))
dev.off()
