## Load libraries
library(tidyverse) # version tidyverse_2.0.0

## PCA without inbred individuals, no sex chromosomes and LD pruned with plink
setwd("/Users/apple/Dropbox/Postdoc_Milan/LK_Joan/PCA/")

eigenvec_table <- read.table('LK_intersect_norelated_wCro_nosexchrom.eigenvec', header = FALSE)

## Add populations
eigenvec_table$Populations <- as.character(c("ESN","ESN","ESN","GRC","GRC","ESN","GRC","GRC","ISR",
                                             "TUR","ESN","KAZ","ISR","ISR","TUR","TUR","CRO","ISR","CRO","ESN",
                                             "SIC","GRC","RUS","RUS","RUS","MOS","MOS","MON","ISR",
                                             "ISR","ESS","ESS","ESS","ESN","ESN","ITN","ITN","ITN",
                                             "SIC","RUS","MON","MON","ITN","ITS","ITS","ITS","ITS",
                                             "SIC","RUS","ISR","ISR","ITS","ITS","SIC","KAZ","RUS",
                                             "MON","CRO","ESS","ESS","ESS","ESS","ITN","ITN","ITN","ITS",
                                             "ITS","GRL","GRL","GRL","GRG","GRG","KAZ","MOS","MOS",
                                             "MON","MON","SIC","SIC","GRL","GRL","GRL","GRG","KAZ"))
## Add population order
pop_order <- c("ESN","ESS","SIC","ITS","ITN","CRO","GRG","GRC","GRL","TUR","ISR","KAZ","RUS","MON","MOS")
## Add colours
palette<-(c("#db0000","#dc2906","#dd4a11","#e16d19","#e68f21","#edb42a","#f5d933","#f5d933",
            "#e4e940","#cad34f","#b0bc66","#98a881","#7f929e","#677dbc","#5068dc"))
names(palette)<-pop_order

## Data wrangling
label <- eigenvec_table$V1
eigenvec_table <- eigenvec_table[-2]
eigenvec_table <- eigenvec_table[-1]

head(eigenvec_table)
for (i in 1:10){
  colnames(eigenvec_table)[i]<-paste0("PC",i)
}

eigenvec_table$Populations <- factor(eigenvec_table$Populations, levels=pop_order)

## Load eigenvalues and calculate percentage of the variance accounted for each PC
eigenval <- read.table('LK_intersect_norelated_wCro_nosexchrom.eigenval', header = F)
percentage <- round(eigenval$V1/sum(eigenval$V1)*100,2)
percentage <- paste0(colnames(eigenvec_table)[1:10]," (",paste(as.character(percentage),"%)"))
percentage

## Plot PC1 vs PC2
pca12<- ggplot(eigenvec_table,aes(x=PC1,y=PC2)) +
  geom_point(aes(colour=Populations), size=4) + 
  scale_color_manual(values=palette) +
  theme_bw(base_family = "Arial") +
  theme(panel.grid = element_blank()) +
  theme(axis.text = element_text(color="black", size=16),
        axis.title = element_text(color="black", size=18),
        legend.title = element_text(color="black", size=18),
        legend.text = element_text(color="black", size=16)) +
  xlab(percentage[1]) +
  ylab(percentage[2])
ggsave("LK_PC1-PC2_wallCro.pdf", pca12, device=cairo_pdf, units="cm", width=20, height=15, limitsize=FALSE)

## Plot PC1 vs PC3
pca13<- ggplot(eigenvec_table,aes(x=PC1,y=PC3)) +
  geom_point(aes(colour=Populations), size=2) + 
  scale_color_manual(values=palette) +
  theme_bw() +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  xlab(percentage[1]) +
  ylab(percentage[3])
ggsave("LK_PC1-PC3_wallCro.pdf", pca13, device="pdf", units="cm", width=20, height=15, limitsize=FALSE)