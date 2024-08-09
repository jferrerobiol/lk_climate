# Lesser kestrel local adaptation and responses to climate change

### Load required libraries
library(readxl)
library(raster)
library(rgdal)
library(tidyverse)
library(pegas)
library(LEA)
library(rnaturalearth)
library(rnaturalearthdata)
library(RColorBrewer)
library(ggpubr)
library(vegan)
library(qvalue)
library(robust)
library(WMDB)
library(ggVennDiagram)
library(cowplot)
library(corrplot)
library(rgeos)
library(ggrepel)
library(sf)
library(foreach)
library(factoextra)
library(vcfR)
library(ggsci)
library(gradientForest)
library(fields)
library(elevatr)
library(reshape2)
library(geosphere)
library(gdm)
library(parallel)
library(doParallel)
library(Hmisc)
library(gplots)
library(viridis)
library(ggspatial)
library(adegenet)

## 1. Loading and formatting data
### Genetic data
#### Loading genotype data
setwd("~/Dropbox/Postdoc_Milan/LK_Joan/intersect_Tasos_Joan_vcfs/vcfs/")
Genotypes <- read.table("./RDA_input.raw", header = T, sep=" ")
Genotypes <- Genotypes[,-(3:6)]
colnames(Genotypes)[1:2] <- c("ind_id","pop_id")

#### Estimating population allele frequencies
AllFreq <- aggregate(Genotypes, by = list(Genotypes$pop_id), function(x) mean(x, na.rm = T)/2)
row.names(AllFreq) <- as.character(AllFreq$Group.1)
AllFreq$ind_id <- NULL
AllFreq$pop_id <- NULL
AllFreq$Group.1 <- NULL
#### There are no columns with missing data for none of the pops, so no need to filter by missing data and impute missing genotypes

#### Filtering on MAF
freq_mean <- colMeans(AllFreq)
AllFreq <- AllFreq[,-which(freq_mean>=0.95 | freq_mean<=0.05)]

#### Ordering loci based on their scaffold
AllFreq <- AllFreq[,order(colnames(AllFreq))]

#### Load population coordinates
setwd("~/Dropbox/Postdoc_Milan/LK_Joan/GEA/")
coord <- as.data.frame(read_excel("genetics_colonies_coords.xlsx"))
rownames(coord) <- c("ISR","ESS","SIC","TUR","ESN","GRC","GRG","GRL","ITS","MOS","CRO","ITN","KAZ","MON","RUS")

#### Load the rasters where the remove.NAs.stack function from Capblancq & Forester 2021 has been applied
ras_current <- stack("~/Dropbox/Postdoc_Milan/LK_Joan/GEA/ras_current.grd")
ras_2040_moderate <- stack("~/Dropbox/Postdoc_Milan/LK_Joan/GEA/chelsa_future_bioclim/ras_2040_moderate.grd")
ras_2070_moderate <- stack("~/Dropbox/Postdoc_Milan/LK_Joan/GEA/chelsa_future_bioclim/ras_2070_moderate.grd")
ras_2040_extreme <- stack("~/Dropbox/Postdoc_Milan/LK_Joan/GEA/chelsa_future_bioclim/ras_2040_extreme.grd")
ras_2070_extreme <- stack("~/Dropbox/Postdoc_Milan/LK_Joan/GEA/chelsa_future_bioclim/ras_2070_extreme.grd")

#### Extracting environmental values for each source population from the rasters
Env <- data.frame(raster::extract(ras_current, coord[,2:3]))

#### Standardization of the variables
Env <- scale(Env, center=TRUE, scale=TRUE) # center=TRUE, scale=TRUE are the defaults for scale()

#### Recovering scaling coefficients
scale_env <- attr(Env, 'scaled:scale')
center_env <- attr(Env, 'scaled:center')

#### Climatic table
Env <- as.data.frame(Env)
row.names(Env) <- rownames(coord)
head(Env)

### Inferring population structure
#### Get mean position in PC1 and 2 from the PCA
setwd("/Users/apple/Dropbox/Postdoc_Milan/LK_Joan/PCA/")

eigenvec_table <- read.table('LK_intersect_LDpruned_norelated_wCro.eigenvec', header = FALSE)

#### Add populations
eigenvec_table$Populations <- as.character(c("ESN","ESN","ESN","GRC","GRC","ESN","GRC","GRC","ISR",
                                             "TUR","ESN","KAZ","ISR","ISR","TUR","TUR","CRO","ISR","CRO","ESN",
                                             "SIC","GRC","RUS","RUS","RUS","MOS","MOS","MON","ISR",
                                             "ISR","ESS","ESS","ESS","ESN","ESN","ITN","ITN","ITN",
                                             "SIC","RUS","MON","MON","ITN","ITS","ITS","ITS","ITS",
                                             "SIC","RUS","ISR","ISR","ITS","ITS","SIC","KAZ","RUS",
                                             "MON","CRO","ESS","ESS","ESS","ESS","ITN","ITN","ITN","ITS",
                                             "ITS","GRL","GRL","GRL","GRG","GRG","KAZ","MOS","MOS",
                                             "MON","MON","SIC","SIC","GRL","GRL","GRL","GRG","KAZ"))

eigenvec_table <- eigenvec_table[,c(2:4,23)]
for (i in 2:3){
  colnames(eigenvec_table)[i]<-paste0("PC",i-1)
}
head(eigenvec_table)

#### Calculate means
eigenvec_pop <- aggregate(eigenvec_table, by = list(eigenvec_table$Populations), function(x) mean(x, na.rm = T))
head(eigenvec_pop)
rownames(eigenvec_pop) <- eigenvec_pop$Group.1
colnames(eigenvec_pop)[1] <- "pop"
eigenvec_pop <- eigenvec_pop[,c(1,3,4)]
head(eigenvec_pop)

#### Plot population structure
pop_order <- c("ESN","ESS","SIC","ITS","ITN","CRO","GRG","GRC","GRL","TUR","ISR","KAZ","RUS","MON","MOS")
palette<-(c("#db0000","#dc2906","#dd4a11","#e16d19","#e68f21","#edb42a","#f5d933","#f5d933",
            "#e4e940","#cad34f","#b0bc66","#98a881","#7f929e","#677dbc","#5068dc"))
names(palette)<-pop_order

pca12<- ggplot(eigenvec_pop,aes(x=PC1,y=PC2)) +
  geom_point(aes(colour=pop), size=2) + 
  scale_color_manual(values=palette) +
  geom_label_repel(aes(label = pop), box.padding   = 0.35, point.padding = 0.5, segment.color = 'grey50') +
  theme_bw() +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0)

pca12

### Table gathering all variables and removing populations with only 3 individuals
Variables <- merge(eigenvec_pop, coord, by=0, all=TRUE)
rownames(Variables) <- Variables$Row.names
Variables <- merge(Variables[-1], Env, by=0, all=TRUE)
Variables <- Variables[-1]
Variables <- Variables[Variables$pop!="TUR" & Variables$pop!="GRG" & Variables$pop!="CRO",]
AllFreq <- AllFreq[rownames(AllFreq)!="TUR" & rownames(AllFreq)!="GRG" & rownames(AllFreq)!="CRO",]
head(Variables)

### Load functions to use later
#### Function to identify outlier SNPs from Forester 2018
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

## 2. GEA analyses using RDA with PCs of bioclimatic variables in order to maximise the amount of variance in allele frequency associated with adaptation to climate
### Perform PCA of the bioclimatic variables to summarise as much environmental information as possible
Env <- Env[rownames(Env)!="TUR" & rownames(Env)!="GRG" & rownames(Env)!="CRO",]
pca_bioclim <- prcomp(Env, scale=F)
fviz_eig(pca_bioclim) # Retain 3 PCs which explain ~90% of the variance
eigenval_pca_bioclim <- pca_bioclim$sdev*pca_bioclim$sdev
var_PCs <- eigenval_pca_bioclim/sum(eigenval_pca_bioclim)

#### Plot populations along PC1 and PC2 to check if populations separate well along the PC space
fviz_pca_ind(pca_bioclim,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

#### Plot variables to see variable contributions and correlation among variables
fviz_pca_var(pca_bioclim,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

res.var <- get_pca_var(pca_bioclim) # Results for Variables
res.var$contrib # Contributions to the PCs

### Merge the first 3 PCs to the Variables dataset
rownames(Variables) <- Variables$pop
PCs <- pca_bioclim$x[,c(1:3)]
colnames(PCs) <- c("bio_PC1","bio_PC2","bio_PC3")
Variables <- merge(Variables, PCs, by=0, all=TRUE)
Variables <- Variables[,-1]

### Run RDA
RDA_env_all <- rda(AllFreq ~ bio_PC1 + bio_PC2 + bio_PC3,  Variables)
RDA_env_all
RsquareAdj(RDA_env_all)

### Assess number of RDA axes to include
screeplot(RDA_env_all, main="Eigenvalues of constrained axes")
summary(eigenvals(RDA_env_all, model = "constrained"))
signif.full <- anova.cca(RDA_env_all, parallel=getOption("mc.cores"))
signif.full
signif.axis <- anova.cca(RDA_env_all, by="axis", parallel=getOption("mc.cores"))
signif.axis # 1 significant axis, so we will focus on SNPs associated with RDA1

### Plot the RDA using ggplot
#### Formatting table for ggplot
locus_scores <- scores(RDA_env_all, choices=c(1:2), display="species", scaling="none") # vegan references "species", here these are the loci
TAB_loci <- data.frame(names = row.names(locus_scores), locus_scores)
pop_scores <- scores(RDA_env_all, choices=c(1:2), display="sites", scaling="none") # vegan references "sites", here these are the populations
TAB_pop <- data.frame(names = row.names(pop_scores), pop_scores)
TAB_var <- as.data.frame(scores(RDA_env_all, choices=c(1,2), display="bp")) # pull the biplot scores

#### Plot
biplot_pop_SNPs_PCs <- ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = TAB_loci, aes(x=RDA1*20, y=RDA2*20), size = 1.4, colour = "gray90") +
  geom_point(data = TAB_pop, aes(x=RDA1*20, y=RDA2*20, colour=names), size = 2) +
  scale_color_manual(values=palette) +
  geom_text_repel(data = TAB_pop, aes(x=RDA1*20, y=RDA2*20, label = names), size=3, family = "Arial") +
  geom_segment(data = TAB_var, aes(xend=RDA1*10, yend=RDA2*10, x=0, y=0), colour="#0868ac", size=0.25, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_text(data = TAB_var, aes(x=5*RDA1, y=5*RDA2-0.5, label = c("bio_PC1","bio_PC2","bio_PC3")), col="#0868ac", size = 3, family = "Helvetica") +
  xlab("RDA 1") + ylab("RDA 2") +
  facet_wrap(~"RDA space") +
  guides(color=guide_legend(title="Population")) +
  theme_bw(base_size = 11, base_family = "Arial") +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))

biplot_pop_SNPs_PCs
ggsave("~/Dropbox/Postdoc_Milan/LK_Joan/GEA/plots/RDA_biplot_pops_SNPs.pdf", biplot_pop_SNPs_PCs, device=cairo_pdf, units="cm", width=18, height=15, limitsize=FALSE)

## 3. Prepare list of candidate SNPs
### Identify candidate SNPs based on SNP loadings on RDA1
load.rda1 <- scores(RDA_env_all, choices=1, display="species")
hist(load.rda1[,1], main="Loadings on RDA1")

ncand <- outliers(load.rda1[,1],3) #338 candidates
ncand_PCs <- names(ncand)

### For each candidate SNP, calculate coefficient of correlation with environment vs coefficient of correlation with structure and vs coefficient of correlation with longitude as a proxy of IBD
all_cand_snps <- data.frame(SNP=character(), R2_adj_clim=numeric(), 
                            R2_adj_struct=numeric(), R2_adj_lon=numeric(), 
                            ratio_clim_struct=numeric(), logratio_clim_struct=numeric(),
                            ratio_clim_lon=numeric(), logratio_clim_lon=numeric(),
                            stringsAsFactors=FALSE)
for (i in 1:length(ncand)) {
  cand_snp <- as.data.frame(AllFreq[,colnames(AllFreq)==names(ncand)[i]])
  rownames(cand_snp) <- rownames(AllFreq)
  colnames(cand_snp) <- "Allele_Freq"
  cand_snp <- cbind(cand_snp,Variables[,c(2,3,5,26,27,28)])
  m_struct <- glm(Allele_Freq ~ PC1 + PC2, data=cand_snp, family="gaussian")
  r_square_adj_struct <- RsquareAdj(m_struct)$adj.r.squared
  r_square_adj_struct <- ifelse(r_square_adj_struct > 0, r_square_adj_struct, 0)
  m_lon <- glm(Allele_Freq ~ lon_pop, data=cand_snp, family="gaussian")
  r_square_adj_lon <- RsquareAdj(m_lon)$adj.r.squared
  r_square_adj_lon <- ifelse(r_square_adj_lon > 0, r_square_adj_lon, 0)
  m_clim <- glm(Allele_Freq ~ bio_PC1 + bio_PC2 + bio_PC3, data=cand_snp, family="gaussian")
  r_square_adj_clim <- RsquareAdj(m_clim)$adj.r.squared
  r_square_adj_clim <- ifelse(r_square_adj_clim > 0, r_square_adj_clim, 0)
  ratio_clim_struct <- r_square_adj_clim/r_square_adj_struct
  log_ratio_clim_struct <- log(ratio_clim_struct)
  ratio_clim_lon <- r_square_adj_clim/r_square_adj_lon
  log_ratio_clim_lon <- log(ratio_clim_lon)
  cand_snp_row <- c(names(ncand)[i], r_square_adj_clim, r_square_adj_struct, 
                    r_square_adj_lon, ratio_clim_struct, log_ratio_clim_struct,
                    ratio_clim_lon, log_ratio_clim_lon)
  all_cand_snps <- rbind(all_cand_snps, cand_snp_row)
}
colnames(all_cand_snps) <- c("SNP", "R2_adj_clim", "R2_adj_struct", "R2_adj_lon", 
                             "ratio_clim_struct","logratio_clim_struct",
                             "ratio_clim_lon","logratio_clim_lon")
all_cand_snps$R2_adj_clim <- as.numeric(all_cand_snps$R2_adj_clim)
all_cand_snps$R2_adj_struct <- as.numeric(all_cand_snps$R2_adj_struct)
all_cand_snps$R2_adj_lon <- as.numeric(all_cand_snps$R2_adj_lon)
all_cand_snps$ratio_clim_struct <- as.numeric(all_cand_snps$ratio_clim_struct)
all_cand_snps$logratio_clim_struct <- as.numeric(all_cand_snps$logratio_clim_struct)
all_cand_snps$ratio_clim_lon <- as.numeric(all_cand_snps$ratio_clim_lon)
all_cand_snps$logratio_clim_lon <- as.numeric(all_cand_snps$logratio_clim_lon)

### Select candidate SNPs that have a higher coefficient of correlation with environment than coefficient of correlation with structure and with longitude and have coefficient of correlation with environment > 0.5 as there are probably meaningful. When more than one SNP per contig exist, take the SNP with the highest R2_adj for climate
sel_cand_snps <- all_cand_snps[all_cand_snps$logratio_clim_struct > 0 &
                                 all_cand_snps$logratio_clim_lon > 0 &
                                 all_cand_snps$R2_adj_clim >= 0.5,]
sel_cand_snps$radtag <- gsub('X([0-9]+).*', '\\1', sel_cand_snps$SNP)
sel_cand_snps <- sel_cand_snps[order(sel_cand_snps$radtag, -sel_cand_snps$R2_adj_clim),] #sort by radtag and reverse sort by R2_adj_clim
sel_cand_snps <- sel_cand_snps[!duplicated(sel_cand_snps$radtag), ] #take the first row within each radtag
sel_cand_snps_RDA <- sel_cand_snps$SNP

### Load list of SNPs that overlapped between PCAdapt and Outflank analyses
#### Remove duplicates to leave only one SNP per locus
fst_outliers <- read.csv("~/Dropbox/Postdoc_Milan/LK_Joan/GEA/outlier_lists_Tasos/top.candicates_intersect84.csv")
fst_outliers <- fst_outliers[fst_outliers$X!=22442,]
fst_outliers$radtag <- gsub('.*_([0-9]+):.*','\\1',fst_outliers$LocusName)
fst_outliers <- fst_outliers[order(fst_outliers$radtag, fst_outliers$LocusName),] #sort by radtag
fst_outliers <- fst_outliers[!duplicated(fst_outliers$radtag), ] #take the first row within each radtag

#### Check overlap between Fst outliers and RDA outliers
snps_RDA <- gsub('X([0-9]+).([0-9]+).*','\\1:\\2', sel_cand_snps_RDA)
snps_fst <- gsub('.*_([0-9]+:[0-9]+):.*','\\1', fst_outliers$LocusName)

overlap_RDA_fst <- list(RDA=snps_RDA, fst=snps_fst)
ggVennDiagram(overlap_RDA_fst, category.names = c("RDA", "Fst"), lty="solid", size=0.2) + 
  scale_fill_gradient2(low = "white", high = 'gray40') + scale_color_manual(values = c("grey", "grey", "grey", "grey")) + guides(fill = "none") + theme(text = element_text(size=16, family = "Helvetica"))

### Generate list containing RDA and Fst outliers
snps_RDA_fst <- unique(c(snps_RDA, snps_fst))

### Retain SNPs that are within genes
#### Read SNPEff annotated VCF file
vcf <- read.vcfR("~/Dropbox/Postdoc_Milan/LK_Joan/GEA/intersect_Tasos_Joan.nomono.nosexchr_ncbi_annotated.vcf", verbose = FALSE )
vcf_info <- vcf@fix[,c(1:3)]
vcf_annot <- data.frame(mut_type=character(), SNPEff_type=character(), gene_ID=character(), gene_type=character(), substitution=character(), 
                        stringsAsFactors=FALSE)
for (i in 1:length(vcf@fix[,8])) {
  vcf_annot <- rbind(vcf_annot, strsplit(vcf@fix[i,8], split='\\|')$INFO[c(2:4,8,10)])
}
colnames(vcf_annot) <- c("mut_type", "SNPEff_type", "gene_ID", "gene_type", "substitution")
vcf_info <- cbind(vcf_info, vcf_annot)
vcf_info$radtag_pos <- gsub('(.*):[+-]','\\1', vcf_info$ID)

#### Select rows corresponding to the list of candidate SNPs from RDA and Fst
vcf_info_cand_snps <- vcf_info[vcf_info$radtag_pos %in% snps_RDA_fst,]

#### Retain only SNPs within genes and remove SNPs in uncharacterised genes
vcf_info_cand_snps_ingenes <- vcf_info_cand_snps[vcf_info_cand_snps$gene_type=="protein_coding",]
vcf_info_cand_snps_ingenes <- vcf_info_cand_snps_ingenes[!startsWith(vcf_info_cand_snps_ingenes$gene_ID,"LOC121"),]

#### Write table with gene information for set of candidate genes
write.table(vcf_info_cand_snps_ingenes, file= "~/Dropbox/Postdoc_Milan/LK_Joan/GEA/SNPS_RDA_12_Fst_gene_info_last.tsv", sep='\t', quote=F, row.names=F)

#### List of final candidate genes
snps_final_set <- vcf_info_cand_snps_ingenes$radtag_pos

### PCA of climate-associated SNPs
vcf_clim_SNPs <- read.vcfR("~/Dropbox/Postdoc_Milan/LK_Joan/GEA/climate_associated_SNPs.vcf")

#### Load popmap and check that samples in the pops and vcf datasets match
popmap <- read.delim("~/Dropbox/Postdoc_Milan/LK_Joan/intersect_Tasos_Joan_vcfs/popmap.txt")
popmap <- popmap[popmap$ind %in% colnames(vcf@gt),]
colnames(vcf_clim_SNPs@gt)[-1] == popmap$ind

#### convert to genlight
gen<-vcfR2genlight(vcf_clim_SNPs)

#### perform PCA
pca<-glPca(gen, nf=30)
var_frac <- pca$eig/sum(pca$eig)

#### pull pca scores out of df
pca.scores <- as.data.frame(pca$scores)
pca.scores$pop <- popmap$pop
pca.scores$pop <- factor(pca.scores$pop, levels=c("ESN","ESS","SIC","ITS","ITN","CRO","GRG","GRC","GRL","TUR","ISR","KAZ","RUS","MON","MOS"))
mypal<-(c("#db0000","#dc2906","#dd4a11","#e16d19","#e68f21","#edb42a","#f5d933","#f5d933",
          "#e4e940","#cad34f","#b0bc66","#98a881","#7f929e","#677dbc","#5068dc"))

#### ggplot color by pop
p_pca_clim_assoc_snps <- ggplot(pca.scores,aes(x=PC1, y=PC2, colour=pop)) +
  geom_point(size=2) +
  scale_color_manual(values=mypal) +
  theme_bw(base_family = "Arial") +
  theme(panel.grid = element_blank()) +
  theme(axis.text = element_text(color="black", size=14),
        axis.title = element_text(size=16)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  coord_fixed(ratio = 1) +
  xlab(paste0("PC1 ","(",round(var_frac[1]*100,2)," %)")) +
  ylab(paste0("PC2 ","(",round(var_frac[2]*100,2)," %)"))
p_pca_clim_assoc_snps
ggsave("~/Dropbox/Postdoc_Milan/LK_Joan/GEA/plots/pca_clim_assoc_snps_PC1-PC2.pdf", p_pca_clim_assoc_snps, device=cairo_pdf, units="cm", width=18, height=15, limitsize=FALSE)

#### ggplot color by cluster
mypal<-(c("#ff6d00","#ff6d00","#ff6d00","#ff6d00","#ff6d00","#ff6d00","#ff6d00","#ff6d00",
          "#ff6d00","#ff6d00","#ff6d00","#5c7ec0","#5c7ec0","#5c7ec0","#5c7ec0"))

p_pca_clim_assoc_snps <- ggplot(pca.scores,aes(x=PC1, y=PC2, colour=pop)) +
  geom_point(size=2) +
  scale_color_manual(values=mypal) +
  theme_bw(base_family = "Arial") +
  theme(panel.grid = element_blank()) +
  theme(axis.text = element_text(color="black", size=14),
        axis.title = element_text(size=16)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  coord_fixed(ratio = 0.5) +
  xlab(paste0("PC1 ","(",round(var_frac[1]*100,2)," %)")) +
  ylab(paste0("PC2 ","(",round(var_frac[2]*100,2)," %)"))
p_pca_clim_assoc_snps
ggsave("~/Dropbox/Postdoc_Milan/LK_Joan/GEA/plots/pca_clim_assoc_snps_PC1-PC2_per_cluster.pdf", p_pca_clim_assoc_snps, device=cairo_pdf, units="cm", width=18, height=15, limitsize=FALSE)

## 4. Allele turnover across the landscape
### Create a list of AllFreq headers that match the candidate SNPs
selection <- vector()
for (i in 1:length(snps_final_set)) {
  selection <- c(selection, colnames(AllFreq)[grep(gsub('([0-9]{4,6})\\:([0-9]{1,3})','X\\1.\\2..', snps_final_set[i]), colnames(AllFreq))])
}

selection <- selection[-c(5,6)] # Remove extra grepped elements

### Prepare the dataframe
gf_data <- data.frame(raster::extract(ras_current, coord[,2:3]))
row.names(gf_data) <- rownames(coord)
gf_data <- merge(gf_data, AllFreq, by=0)
rownames(gf_data) <- gf_data$Row.names

#### Separate the data between candidate and rest of SNPs
pot_outliers <- c(gsub("(X[0-9]+.)[0-9]+.*","\\1",ncand_PCs), 
                  gsub("([0-9]+)","X\\1.",fst_outliers$radtag)) #list of radtags containing potential outliers for either of the RDAs or the Fst scans
pot_outliers <- unique(pot_outliers)
pot_outliers_n <- colnames(gf_data)[grep(paste(pot_outliers, collapse="|"), colnames(gf_data))] #grep the radtags containing at least one potential outlier

gf_data_candidate <- gf_data[,selection]
gf_data_background <- gf_data[,!(colnames(gf_data) %in% pot_outliers_n)]
gf_data_background <- gf_data_background[,-c(1:20)]

#### Create a table with the bioclimatic information of populations
bioclim_present <- gf_data[,c(grep("bio",names(gf_data)))]

#### Define the bioclimatic variables of interest
bioclimatic <- paste("bio_",1:19,sep = "")

#### set the importance of the permutation distribution of each variable.
maxLevel <- log2(0.368*nrow(gf_data_candidate)/2)

### Run the gradient forest function for each set of SNPs
gf_candidate <- gradientForest(cbind(bioclim_present[,bioclimatic], gf_data_candidate),
                               predictor.vars=colnames(bioclim_present[,bioclimatic]),
                               response.vars=colnames(gf_data_candidate), ntree=500,
                               maxLevel=maxLevel, trace=T, corr.threshold=0.50)
gf_background <- gradientForest(cbind(bioclim_present[,bioclimatic], gf_data_background),
                                predictor.vars=colnames(bioclim_present[,bioclimatic]),
                                response.vars=colnames(gf_data_background), ntree=500,
                                maxLevel=maxLevel, trace=T, corr.threshold=0.50)
gf_runs <- list(gf_background=gf_background,
                gf_candidate=gf_candidate) # combine the GF models into a list and save it

### Estimate the importance of each variable in the model
bio_cand <- gf_candidate$overall.imp[order(gf_candidate$overall.imp,decreasing = T)]
bio_most_imp_cand <- names(bio_cand[1])
bio_cand_d <- as.data.frame(bio_cand)
bio_cand_d$bio_var <- row.names(bio_cand_d)
bio_cand_d$bio_var <- factor(bio_cand_d$bio_var, levels=bio_cand_d$bio_var[order(bio_cand_d$bio_cand)])

mycolours <- colorRampPalette(c("#e4ece2","#bae0c1","#95d5b5","#73c9ba","#54adbe","#3877b2","#1f38a7"))(19)
axis_col <- c("#3877b2","#3877b2","#E65A4C","#E65A4C",rep("#3877b2",4),rep("#E65A4C",7),"#3877b2","#E65A4C","#E65A4C","#3877b2")
#### Barplot of variable weighted importance
p_bio_imp <- ggplot(bio_cand_d) +
  theme_minimal() +
  theme(panel.grid=element_blank(), axis.ticks.x = element_line(size=0.5), 
        text = element_text(colour = "black", family = "Arial"), axis.text.y = element_text(colour = axis_col)) +
  geom_col(aes(x=bio_cand, y=bio_var, fill=bio_var)) +
  scale_fill_manual(values=mycolours, guide="none") +
  scale_y_discrete(labels=c("BIO18","BIO14","BIO5","BIO10","BIO16","BIO17","BIO13","BIO15","BIO2","BIO3","BIO9"
                            ,"BIO1","BIO8","BIO11","BIO6","BIO12","BIO4","BIO7","BIO19")) +
  geom_segment(aes(x = 0, y = 0, xend = 0.003, yend = 0), size=0.5) +
  xlab(expression(R^2~" weighted importance")) + ylab("")

p_bio_imp
ggsave("~/Dropbox/Postdoc_Milan/LK_Joan/GEA/plots/GF_bioclim_imp.pdf", p_bio_imp, device=cairo_pdf, units="cm", width=6, height=11, limitsize=FALSE)

### Allele turnover functions across the landscape
#### for BIO19
temp_cand_overall_bio19 <- cumimp(gf_candidate,predictor= "bio_19",
                                  type=c("Overall"),standardize = T) # all candidate SNPs
temp_cand_SNP_bio19 <- cumimp(gf_candidate,predictor = "bio_19",
                              type=c("Species"),standardize = T) #each individual candidate allele
temp_ref_overall_bio19 <- cumimp(gf_background,predictor = "bio_19",
                                 type=c("Overall"),standardize = T) #all neutral SNPs
temp_ref_SNP_bio19 <- cumimp(gf_background,predictor = "bio_19",
                             type=c("Species"),standardize = T) #each individidual neutral SNPs

#### Incorporate the populations into the allele turnover functions
pop_turn <- predict(gf_candidate,bioclim_present[,grep("bio",names(bioclim_present))])
temp_bio19 <- data.frame(bio=bioclim_present[,"bio_19"],imp=pop_turn[,"bio_19"]) # get the x (biovalue) and y (predicted cumulative importance) values of each population

#### for BIO7
temp_cand_overall_bio7 <- cumimp(gf_candidate,predictor= "bio_7",
                                 type=c("Overall"),standardize = T) # all candidate SNPs
temp_cand_SNP_bio7 <- cumimp(gf_candidate,predictor = "bio_7",
                             type=c("Species"),standardize = T) #each individual candidate allele
temp_ref_overall_bio7 <- cumimp(gf_background,predictor = "bio_7",
                                type=c("Overall"),standardize = T) #all neutral SNPs
temp_ref_SNP_bio7 <- cumimp(gf_background,predictor = "bio_7",
                            type=c("Species"),standardize = T) #each individidual neutral SNPs

#### Incorporate the populations into the allele turnover functions
temp_bio7 <- data.frame(bio=bioclim_present[,"bio_7"],imp=pop_turn[,"bio_7"]) # get the x (biovalue) and y (predicted cumulative importance) values of each population

#### Plot showing populations
pop_col <- c("#db0000","#dc2906","#f5d933","#e4e940","#b0bc66","#e68f21","#e16d19","#98a881","#677dbc","#5068dc","#7f929e","#dd4a11")
par(mfrow=c(1,2))
par(mai=c(0.9,0.8,0.4,0))
plot(temp_cand_overall_bio19,type="n",xlim=c(0,310), ylim=c(0,0.16),mgp=c(2,0.6,0), ylab="Cumulative importance",xlab="Precipitation of coldest quarter (mm)")
for(j in 1:length(temp_cand_SNP_bio19)){
  lines(temp_cand_SNP_bio19[[j]],col=adjustcolor(MaizePal::maize_pal("MaizAzul")[5],alpha.f=0.6))
}
lines(temp_cand_overall_bio19,col=MaizePal::maize_pal("MaizAzul")[2],lwd=4)
lines(temp_ref_overall_bio19,col=gray(.10),lwd=4)
points(temp_bio19$bio,temp_bio19$imp,pch=21,bg=pop_col,cex=1.5)

par(mai=c(0.9,0.1,0.4,0.6),tcl=-0.2)
plot(temp_cand_overall_bio7,type="n",xlim=c(23,56), ylim=c(0,0.16),mgp=c(2,0.6,0), ylab="Cumulative importance",xlab="Temperature annual range (ºC)")
for(j in 1:length(temp_cand_SNP_bio7)){
  lines(temp_cand_SNP_bio7[[j]],col=adjustcolor(MaizePal::maize_pal("RubyGold")[5],alpha.f=0.6))
}
lines(temp_cand_overall_bio7,col=MaizePal::maize_pal("RubyGold")[2],lwd=4)
lines(temp_ref_overall_bio7,col=gray(.10),lwd=4)
points(temp_bio7$bio,temp_bio7$imp,pch=21,bg=pop_col,cex=1.5)

## 5. Adaptive landscape: projecting adaptive gradient across space
### Calculate correlations among bioclimatic variables to choose one among those correlated (|r| > 0.7) and following the GF weighted importance
cor_var <- cor(Variables[,bioclimatic])
corrplot(cor_var, method = 'number')  # 

### Perform an adaptively enriched RDA with the 65 candidate SNPs and the most important variables from the GF analysis that are not correlated to other more important variables
RDA_outliers <- rda(AllFreq[,selection] ~ bio_19 + bio_8 + bio_2 + bio_15 + bio_10,  Variables)

#### Plot the RDA biplot
TAB_loci <- as.data.frame(scores(RDA_outliers, choices=c(1:2), display="species", scaling="none"))
TAB_var <- as.data.frame(scores(RDA_outliers, choices=c(1:2), display="bp"))
p_enriched_RDA <- ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = TAB_loci, aes(x=RDA1*3, y=RDA2*3), colour = "#EB8055FF", size = 2, alpha = 0.8) + #"#F9A242FF"
  geom_segment(data = TAB_var, aes(xend=RDA1, yend=RDA2, x=0, y=0), colour="black", size=0.15, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_text(data = TAB_var, aes(x=1.1*RDA1, y=1.1*RDA2, label = row.names(TAB_var)), size = 2.5, family = "Arial") +
  xlab("RDA 1 (67%)") + ylab("RDA 2 (21%)") +
  facet_wrap(~"Adaptively enriched RDA space") +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11, base_family = "Helvetica") +
  theme(panel.grid = element_blank(), plot.background = element_blank(), panel.background = element_blank(), strip.text = element_text(size=11))

p_enriched_RDA
ggsave("~/Dropbox/Postdoc_Milan/LK_Joan/GEA/plots/cand_SNPs_RDA_space.pdf", p_enriched_RDA, device=cairo_pdf, units="cm", width=15, height=15, limitsize=FALSE)

#### Upload the map and species range shape file
admin <- ne_countries(scale = "medium", returnclass = "sf")
range <- readOGR("~/Dropbox/Postdoc_Milan/LK_Joan/GEA/current_SDM_shapefile/lower_res/current_suitable_breeding_range_buff005degrees.shp") 
crs(range) <- '+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs'

#### Load the function to predict the adaptive index across the landscape
source("~/Dropbox/Postdoc_Milan/LK_Joan/GEA/RDA_Capblanq_example/RDA-landscape-genomics/src/adaptive_index.R")

#### Running the function for all the climatic pixels of Lesser kestrel distribution range
res_RDA_proj_current <- adaptive_index(RDA = RDA_outliers, K = 1, env_pres = ras_current, range = range, method = "loadings", scale_env = scale_env, center_env = center_env)

#### Vectorization of the climatic rasters for ggplot
RDA_proj <- list(res_RDA_proj_current$RDA1)
RDA_proj <- lapply(RDA_proj, function(x) rasterToPoints(x))
for(i in 1:length(RDA_proj)){
  RDA_proj[[i]][,3] <- (RDA_proj[[i]][,3]-min(RDA_proj[[i]][,3]))/(max(RDA_proj[[i]][,3])-min(RDA_proj[[i]][,3]))
}

#### Adaptive genetic turnover projected across Lesser kestrel range for RDA1
TAB_RDA <- as.data.frame(do.call(rbind, RDA_proj[1:2]))
colnames(TAB_RDA)[3] <- "value"
TAB_RDA$variable <- "RDA1"
TAB_RDA$variable <- as.factor(TAB_RDA$variable)

p_adaptive_landscape <- ggplot(data = TAB_RDA) + 
  geom_sf(data = admin, fill=gray(.8), colour=gray(.8), size=0.05) +
  geom_tile(aes(x = x, y = y, fill = cut(value, breaks=seq(0, 1, length.out=10), include.lowest = T))) + 
  scale_fill_viridis_d(alpha = 0.8, direction = -1, option = "A", labels = c("Negative scores","","","","Intermediate scores","","","","Positive scores")) +
  coord_sf(xlim = c(150, -20), ylim = c(20, 65), expand = F) +
  scale_x_continuous(breaks=seq(-20,150,20)) +
  xlab("Longitude") + ylab("Latitude") +
  guides(fill=guide_legend(title="Adaptive index")) +
  theme_bw(base_size = 11, base_family = "Arial") +
  theme(panel.grid = element_blank(), plot.background = element_blank(), panel.background = element_blank(), strip.text = element_text(size=11))

p_adaptive_landscape
ggsave("~/Dropbox/Postdoc_Milan/LK_Joan/GEA/plots/adaptive_landscape_bio-19-8-2-10-15.pdf", p_adaptive_landscape, device = cairo_pdf, units="cm", width=18, height=6, limitsize=FALSE)

## 6. Manhattan plot of SNPs associated with RDA1
### Get info to plot Manhattan
lims <- mean(load.rda1) + c(-1, 1) * 3 * sd(load.rda1)
SNP_ID <- gsub('X([0-9]+).([0-9]+).*','\\1:\\2', colnames(AllFreq))
Outliers <- rep("Neutral", length(colnames(AllFreq)))
Outliers[colnames(AllFreq)%in%selection] <- "Outlier"
Outliers <- factor(Outliers, levels = c("Neutral", "Outlier"))

TAB_manhatan <- data.frame(radtag_pos = SNP_ID, 
                           loadings = load.rda1, 
                           Outliers = Outliers)
TAB_manhatan <- TAB_manhatan %>% left_join(vcf_info)
TAB_manhatan <- TAB_manhatan[,c(1:5,7,9)]
chrom_names <- read.delim("~/Dropbox/Postdoc_Milan/LK_Joan/GEA/correspondence_chrom_names.txt", header=F)
colnames(chrom_names) <- c("chromosome", "CHROM")
TAB_manhatan <- TAB_manhatan %>% left_join(chrom_names)
TAB_manhatan$RDA1_pol <- ifelse(TAB_manhatan$RDA1<0,-TAB_manhatan$RDA1,TAB_manhatan$RDA1)
Outliers_pol <- rep("Neutral", length(colnames(AllFreq)))
Outliers_pol[colnames(AllFreq)%in%ncand_PCs] <- "Outlier"
Outliers_pol[colnames(AllFreq)%in%sel_cand_snps_RDA] <- "Climate"
Outliers_pol[colnames(AllFreq)%in%selection] <- "Selection"
Outliers_pol <- factor(Outliers_pol, levels = c("Neutral","Outlier","Climate","Selection"))
TAB_manhatan$Outliers_pol <- Outliers_pol
TAB_manhatan$POS <- as.numeric(TAB_manhatan$POS)
TAB_manhatan$chromosome <- as.numeric(TAB_manhatan$chromosome)
TAB_manhatan <- TAB_manhatan[order(TAB_manhatan$chromosome,TAB_manhatan$POS),]

TAB_manhatan$BPcum <- NA
TAB_manhatan$BPcum[1] <- as.numeric(TAB_manhatan$POS[1])
acum = 0

for (i in 2:length(TAB_manhatan$POS)) {
  if(TAB_manhatan$chromosome[i] == TAB_manhatan$chromosome[i-1]) {
    TAB_manhatan$BPcum[i] <- as.numeric(TAB_manhatan$POS[i]) + acum
  }
  if(TAB_manhatan$chromosome[i] != TAB_manhatan$chromosome[i-1]) {
    acum = acum + as.numeric(TAB_manhatan$POS[i-1])
    TAB_manhatan$BPcum[i] <- as.numeric(TAB_manhatan$POS[i]) + acum
  }
}

mean_chrom <- vector()
for(i in levels(as.factor(TAB_manhatan$chromosome))){
  mean_chrom[i] <- (min(TAB_manhatan$BPcum[TAB_manhatan$chromosome==i]) + max(TAB_manhatan$BPcum[TAB_manhatan$chromosome==i]))/2
}
mean_chrom

# Manhattan plot
p_manhattan_climate <- ggplot(TAB_manhatan[TAB_manhatan$Outliers=="Neutral",], aes(x=BPcum, y=RDA1)) +
  # Show all points
  geom_point(aes(color=as.factor(chromosome)), size=1) +
  scale_color_manual(values = rep(c("gray80", "gray50"), 10 )) +
  geom_point(data = TAB_manhatan[TAB_manhatan$Outliers=="Outlier",], colour="#f74343FF", size=1) +
  geom_hline(yintercept=lims, linetype="dashed", color = "black", size=0.6) +
  # custom X axis:
  scale_x_continuous(label = names(mean_chrom), breaks = mean_chrom) +
  scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  xlab("Chromosome") + ylab("Loadings on RDA1") +
  # Custom the theme:
  theme_classic() +
  theme(legend.position="none", axis.text = element_text(color="black", size=14),
        axis.title = element_text(size=16))
p_manhattan_climate

# Polarised Manhattan plot
p_manhattan_climate_pol <- ggplot(TAB_manhatan[TAB_manhatan$Outliers=="Neutral",], aes(x=BPcum, y=RDA1_pol)) +
  # Show all points
  geom_point(aes(color=as.factor(chromosome)), size=1) +
  scale_color_manual(values = rep(c("gray80", "gray50"), 10 )) +
  geom_point(data = TAB_manhatan[TAB_manhatan$Outliers_pol=="Outlier",], colour="#fddea0", size=1) +
  geom_point(data = TAB_manhatan[TAB_manhatan$Outliers_pol=="Climate",], colour="#f97b5d", size=1) +
  geom_point(data = TAB_manhatan[TAB_manhatan$Outliers_pol=="Selection",], colour="#bd3977", size=1) +
  # custom X axis:
  scale_x_continuous(label = names(mean_chrom), breaks = mean_chrom) +
  scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  xlab("Chromosome") + ylab("Loadings on RDA1") +
  # Custom the theme:
  theme_classic() +
  theme(legend.position="none", axis.text = element_text(color="black", size=14),
        axis.title = element_text(size=16))
p_manhattan_climate_pol

ggsave("~/Dropbox/Postdoc_Milan/LK_Joan/GEA/plots/manhattan_climate.pdf", p_manhattan_climate_pol, device=cairo_pdf, units="cm", width=20, height=7, limitsize=FALSE)

### 7. Spearman correlation between candidate SNPs and bioclimatic variables
cand_genes_function <- read_excel("~/Dropbox/Postdoc_Milan/LK_Joan/GEA/Local_adapt_SNPS_RDA_12_Fst_gene_info_last.xlsx")
records_cand <- gsub("^X([0-9]+)\\.([0-9]+).*","\\1:\\2",colnames(gf_data_candidate)) %in% cand_genes_function$radtag_pos
cand_snps_data <- gf_data_candidate[,records_cand]
cand_snps_data <- merge(cand_snps_data,Env, by=0)
cand_snps_data <- cand_snps_data[-1]
cand_snps_data <- as.matrix(cand_snps_data)
cand_snps_env_var_corr <- rcorr(cand_snps_data, type=c("spearman"))
corr_coef_cand_snps_env_var <- cand_snps_env_var_corr$r[c(62:80),c(1:61)]
corr_coef_cand_snps_env_var <- abs(corr_coef_cand_snps_env_var)

#### Get gene names to print them in the heatmap
d_heatmap <- data.frame(colnames_heatmap_dataset=colnames(corr_coef_cand_snps_env_var),
                        radtag_pos=gsub("^X([0-9]+)\\.([0-9]+).*","\\1:\\2",colnames(corr_coef_cand_snps_env_var)))
d_gene_ID_radtagpos <- cand_genes_function %>% select(radtag_pos,gene_ID)
d_heatmap <- d_heatmap %>% left_join(d_gene_ID_radtagpos)
colnames(corr_coef_cand_snps_env_var) <- d_heatmap$gene_ID

#### Plot heatmap of correlations between SNPs and climatic variables
colours <- colorRampPalette(brewer.pal(n=9,"YlOrRd"))
some.colors <- colours(100)

heatmap.2(corr_coef_cand_snps_env_var, trace="none", col=some.colors)

## 8. Plot heatmap of candidate genes literature associations
cand_genes_function$Loc_adapt <- ifelse(cand_genes_function$Local_adaptation %in% "x","yes","no")
cand_genes_function$Urb_adapt_domest <- ifelse(cand_genes_function$Adaptationurban_Domestication_selection_domestic_animals %in% "x","yes","no")
cand_genes_function$Phen_tr <- ifelse(cand_genes_function$Phenotypic_traits_relevant_local_adaptation %in% "x","yes","no")
cand_genes_function$str_res <- ifelse(cand_genes_function$Stress_response %in% "x","yes","no")
cand_genes_function <- mutate(cand_genes_function,
                              num = case_when(
                                Loc_adapt == "yes" &  Phen_tr == "yes" &  Urb_adapt_domest == "yes" &  str_res == "yes" ~ "four", 
                                Loc_adapt == "yes" &  Phen_tr == "yes" &  Urb_adapt_domest == "yes" &  str_res == "no" ~ "three_1", 
                                Loc_adapt == "yes" &  Phen_tr == "yes" &  Urb_adapt_domest == "no" &  str_res == "yes" ~ "three_2", 
                                Loc_adapt == "yes" &  Phen_tr == "no" &  Urb_adapt_domest == "yes" &  str_res == "yes" ~ "three_3", 
                                Loc_adapt == "no" &  Phen_tr == "yes" &  Urb_adapt_domest == "yes" &  str_res == "yes" ~ "three_4", 
                                Loc_adapt == "yes" &  Phen_tr == "yes" &  Urb_adapt_domest == "no" &  str_res == "no" ~ "two_1", 
                                Loc_adapt == "yes" &  Phen_tr == "no" &  Urb_adapt_domest == "yes" &  str_res == "no" ~ "two_2", 
                                Loc_adapt == "yes" &  Phen_tr == "no" &  Urb_adapt_domest == "no" &  str_res == "yes" ~ "two_3", 
                                Loc_adapt == "no" &  Phen_tr == "yes" &  Urb_adapt_domest == "yes" &  str_res == "no" ~ "two_4", 
                                Loc_adapt == "no" &  Phen_tr == "yes" &  Urb_adapt_domest == "no" &  str_res == "yes" ~ "two_5", 
                                Loc_adapt == "no" &  Phen_tr == "no" &  Urb_adapt_domest == "yes" &  str_res == "yes" ~ "two_6", 
                                Loc_adapt == "yes" &  Phen_tr == "no" &  Urb_adapt_domest == "no" &  str_res == "no" ~ "one_1", 
                                Loc_adapt == "no" &  Phen_tr == "yes" &  Urb_adapt_domest == "no" &  str_res == "no" ~ "one_2", 
                                Loc_adapt == "no" &  Phen_tr == "no" &  Urb_adapt_domest == "yes" &  str_res == "no" ~ "one_3", 
                                Loc_adapt == "no" &  Phen_tr == "no" &  Urb_adapt_domest == "no" &  str_res == "yes" ~ "one_4", 
                                Loc_adapt == "no" &  Phen_tr == "no" &  Urb_adapt_domest == "no" &  str_res == "no" ~ "zero"))
cand_genes_function_sel <- cand_genes_function[,c(7,28:32)]
cand_genes_function_long <- melt(cand_genes_function_sel,id.vars = c("gene_ID","num"))
cand_genes_function_long$variable <- factor(cand_genes_function_long$variable, levels=c("Loc_adapt","Phen_tr","Urb_adapt_domest","str_res"))
cand_genes_function_long$gene_ID <- factor(cand_genes_function_long$gene_ID)
cand_genes_function_long <- mutate(cand_genes_function_long,
                                   num_num = case_when(
                                     num=="four" ~ 1,
                                     num=="three_1" ~ 2,
                                     num=="three_2" ~ 3,
                                     num=="three_3" ~ 4,
                                     num=="three_4" ~ 5,
                                     num=="two_1" ~ 6,
                                     num=="two_2" ~ 7,
                                     num=="two_3" ~ 8,
                                     num=="two_4" ~ 9,
                                     num=="two_5" ~ 10,
                                     num=="two_6" ~ 12,
                                     num=="one_1" ~ 13,
                                     num=="one_2" ~ 14,
                                     num=="one_3" ~ 15,
                                     num=="one_4" ~ 16,
                                     num=="zero" ~ 17
                                   ))
cand_genes_function_long$gene_ID <- factor(cand_genes_function_long$gene_ID, levels=c("ZNF280D","DPP4","IRAG2","ZFAND3","OSBPL6","GALNTL6","RALYL","HERC1","TACC3","LRP4","LDB3","FAR1","ALS2","MYBBP1A","RORA","RNF43","CNTNAP5","GPR22","EXOC3L4","PPM1H","SLC2A9","CDK6","DYNC1H1","CASK","ITPK1","SPEN","MYBPC3","FSTL4","CPNE4","INTS12","EMP2","ANKS6","CFAP61","LMBR1","MAEA","NEO1","EIF4ENIF1","ADORA2A","KCNMB2","ZBTB20","MSL2","EVI5","PAK3","RSPH14","SLC5A1","PTPN5","UBE2K","AATF","MGRN1","PRTG","EPB41L2","ADGRB3","TECRL","CDH8","CACNA1C","BRSK2","RSF1","FASN")) # order in the sasme way than the heatmap

p_cand_genes_heatmap <- ggplot(cand_genes_function_long, aes(x = variable, y = gene_ID, fill = value)) +
  theme_minimal(base_size = 11, base_family = "Arial") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_tile(color="white") +
  scale_fill_manual(values=c("white", "#607D8B"), guide="none") +
  xlab("") + ylab("Gene ID")
p_cand_genes_heatmap

ggsave("~/Dropbox/Postdoc_Milan/LK_Joan/GEA/plots/candidate_genes_heatmap.pdf", p_cand_genes_heatmap, device = cairo_pdf, units="cm", width=5, height=18, limitsize=FALSE)

## 9. Calculate local genetic offsets with RDA
### Function to predict genomic offset from a RDA model
source("~/Dropbox/Postdoc_Milan/LK_Joan/GEA/RDA_Capblanq_example/RDA-landscape-genomics/src/genomic_offset.R")

### Running the function for the 3 future scenarios
res_RDA_proj2040mod <- genomic_offset(RDA_outliers, K = 1, env_pres = ras_current, env_fut = ras_2040_moderate, range = range, method = "loadings", scale_env = scale_env, center_env = center_env)
res_RDA_proj2070mod <- genomic_offset(RDA_outliers, K = 1, env_pres = ras_current, env_fut = ras_2070_moderate, range = range, method = "loadings", scale_env = scale_env, center_env = center_env)
res_RDA_proj2040ext <- genomic_offset(RDA_outliers, K = 1, env_pres = ras_current, env_fut = ras_2040_extreme, range = range, method = "loadings", scale_env = scale_env, center_env = center_env)
res_RDA_proj2070ext <- genomic_offset(RDA_outliers, K = 1, env_pres = ras_current, env_fut = ras_2070_extreme, range = range, method = "loadings", scale_env = scale_env, center_env = center_env)

## Tables containing global genetic offsets predicted for 2041-2070 under moderate and extreme scenarios and for 2071-2100 under extreme scenario
RDA_proj_offset_2040mod <- data.frame(rbind(rasterToPoints(res_RDA_proj2040mod$Proj_offset_global), 
                                            rasterToPoints(res_RDA_proj2040mod$Proj_offset_global)))
RDA_proj_offset_2040mod$Global_offset_b <- ifelse(RDA_proj_offset_2040mod$Global_offset<2,RDA_proj_offset_2040mod$Global_offset,2.1)
RDA_proj_offset_2040ext <- data.frame(rbind(rasterToPoints(res_RDA_proj2040ext$Proj_offset_global), 
                                            rasterToPoints(res_RDA_proj2040ext$Proj_offset_global)))
RDA_proj_offset_2040ext$Global_offset_b <- ifelse(RDA_proj_offset_2040ext$Global_offset<2,RDA_proj_offset_2040ext$Global_offset,2.1)
RDA_proj_offset_2070ext <- data.frame(rbind(rasterToPoints(res_RDA_proj2070ext$Proj_offset_global), 
                                            rasterToPoints(res_RDA_proj2070ext$Proj_offset_global)))
RDA_proj_offset_2070ext$Global_offset_b <- ifelse(RDA_proj_offset_2070ext$Global_offset<2,RDA_proj_offset_2070ext$Global_offset,2.1)

## Projecting genomic offset on a map
colfunc <- colorRampPalette(c("#FFFFBF","#FEE08B","#FDAE61","#F46D43","#D53E4F","#9E0142"))
colors <- colfunc(5)

p_genom_offset_2040_moderate <- ggplot(data = RDA_proj_offset_2040mod) + 
  geom_sf(data = admin, fill=gray(.8), colour=gray(.8), size=0.05) +
  geom_tile(aes(x = x, y = y, fill = cut(Global_offset_b, breaks=seq(0, 2.5, by = 0.5), include.lowest = T)), alpha = 1) + 
  scale_fill_manual(values = colors, labels = c("0-0.5","0.5-1.0","1.0-1.5","1.5-2.0","> 2"),
                    guide = guide_legend(title="Genetic offset", title.position = "top",
                                         title.hjust = 0.5, ncol = 1, label.position="right"), na.translate = F) +
  coord_sf(xlim = c(150, -20), ylim = c(20, 65), expand = F) +
  xlab("Longitude") + ylab("Latitude") +
  theme_bw(base_size = 11, base_family = "Arial") +
  theme(panel.grid = element_blank(), plot.background = element_blank(), panel.background = element_blank(), strip.text = element_text(size=11))

p_genom_offset_2040_moderate
ggsave("~/Dropbox/Postdoc_Milan/LK_Joan/GEA/plots/genom_offset_2040_moderate_bio19-8-2-15-10.pdf", p_genom_offset_2040_moderate, device = cairo_pdf, units="cm", width=18, height=6, limitsize=FALSE)

p_genom_offset_2040_extreme <- ggplot(data = RDA_proj_offset_2040ext) + 
  geom_sf(data = admin, fill=gray(.8), colour=gray(.8), size=0.05) +
  geom_tile(aes(x = x, y = y, fill = cut(Global_offset_b, breaks=seq(0, 2.5, by = 0.5), include.lowest = T)), alpha = 1) + 
  scale_fill_manual(values = colors, labels = c("0-0.5","0.5-1.0","1.0-1.5","1.5-2.0","> 2"),
                    guide = guide_legend(title="Genetic offset", title.position = "top",
                                         title.hjust = 0.5, ncol = 1, label.position="right"), na.translate = F) +
  coord_sf(xlim = c(150, -20), ylim = c(20, 65), expand = F) +
  xlab("Longitude") + ylab("Latitude") +
  theme_bw(base_size = 11, base_family = "Arial") +
  theme(panel.grid = element_blank(), plot.background = element_blank(), panel.background = element_blank(), strip.text = element_text(size=11))

p_genom_offset_2040_extreme
ggsave("~/Dropbox/Postdoc_Milan/LK_Joan/GEA/plots/genom_offset_2040_extreme_bio19-8-2-15-10.pdf", p_genom_offset_2040_extreme, device = cairo_pdf, units="cm", width=18, height=6, limitsize=FALSE)

p_genom_offset_2070_extreme <- ggplot(data = RDA_proj_offset_2070ext) + 
  geom_sf(data = admin, fill=gray(.8), colour=gray(.8), size=0.05) +
  geom_tile(aes(x = x, y = y, fill = cut(Global_offset_b, breaks=seq(0, 2.5, by = 0.5), include.lowest = T)), alpha = 1) + 
  scale_fill_manual(values = colors, labels = c("0-0.5","0.5-1.0","1.0-1.5","1.5-2.0","> 2"),
                    guide = guide_legend(title="Genetic offset", title.position = "top",
                                         title.hjust = 0.5, ncol = 1, label.position="right"), na.translate = F) +
  coord_sf(xlim = c(150, -20), ylim = c(20, 65), expand = F) +
  xlab("Longitude") + ylab("Latitude") +
  theme_bw(base_size = 11, base_family = "Arial") +
  theme(panel.grid = element_blank(), plot.background = element_blank(), panel.background = element_blank(), strip.text = element_text(size=11))

p_genom_offset_2070_extreme
ggsave("~/Dropbox/Postdoc_Milan/LK_Joan/GEA/plots/genom_offset_2070_extreme_bio19-8-2-15-10.pdf", p_genom_offset_2070_extreme, device = cairo_pdf, units="cm", width=18, height=6, limitsize=FALSE)

## 10. Calculate genomic offsets for the two main groups
### Load breeding data for each of the groups
east <- read.csv("~/Dropbox/Postdoc_Milan/LK_Joan/GEA/old_Mattia/eastern_data_bioclim_5km_buffer_no_dupl.csv")
east <- east[,c(5,12,13)]
east$group <- "eastern"
west <- read.csv("~/Dropbox/Postdoc_Milan/LK_Joan/GEA/old_Mattia/western_data_bioclim_5km_buffer_no_dupl.csv")
west <- west[,c(5,13,14)]
west$group <- "western"
presence <- rbind(east,west)
colnames(presence) <- c("country","x","y","group")

### Find nearest neighbour
#### Reorganise genomic offset data frame
genom_offset <- cbind(RDA_proj_offset_2040mod,RDA_proj_offset_2040ext,RDA_proj_offset_2070ext)
genom_offset <- genom_offset[,c(1:3,7,11)]
colnames(genom_offset) <- c("x","y","Global_offset_2040_mod","Global_offset_2040_ext","Global_offset_2070_ext")
#### Mask raster to range shapefile
env_pres <- raster::mask(ras_current, range)
#### Extract raster cell numbers with presence data
presence_sp <- data.frame(raster::extract(env_pres, presence[,2:3], cellnumbers=T))
#### Convert raster cell numbers to lat lon and add to presence_sp data frame
presence_sp <- cbind(presence_sp,as.data.frame(xyFromCell(env_pres, presence_sp$cells)))
presence_sp <- presence_sp[,c(21,22)]
#### Merge genomic offset and presence data
presence_correspond <- cbind(presence_sp,presence)
colnames(presence_correspond) <- c("x_ras","y_ras","country","x","y","group")
genom_offset_presence <- merge(presence_correspond,genom_offset, by.x=c("x_ras", "y_ras"), by.y=c("x","y"))
genom_offset_presence$group <- factor(genom_offset_presence$group, levels=c("western","eastern"))
genom_offset_presence <- unique(genom_offset_presence)

### Plot violin plots of genomic offset for eastern and western clade
p_gen_off_groups_2040 <- ggplot(data=genom_offset_presence, aes(x=group,y=Global_offset_2040_ext)) +
  theme_bw(base_size = 11, base_family = "Arial") +
  theme(panel.grid = element_blank()) +
  geom_violin(aes(colour=group),width=0.3) +
  geom_boxplot(aes(colour=group),width=0.2) +
  geom_jitter(aes(colour=group),alpha=0.5,width=0.05) +
  scale_fill_manual(values=c(MaizePal::maize_pal("RubyGold")[2],MaizePal::maize_pal("MaizAzul")[2]), guide="none") +
  scale_colour_manual(values=c(MaizePal::maize_pal("RubyGold")[2],MaizePal::maize_pal("MaizAzul")[2]), guide="none") +
  ylim(0,5.5) +
  ylab("Genetic offset")
p_gen_off_groups_2040

### Plot density plots of genomic offset for eastern and western clade
median_west <- median(genom_offset_presence$Global_offset_2040_mod[genom_offset_presence$group=="western"])
median_east <- median(genom_offset_presence$Global_offset_2040_mod[genom_offset_presence$group=="eastern"])

p_gen_off_groups_2040_mod <- ggplot(data=genom_offset_presence, aes(y=Global_offset_2040_mod,fill=group)) +
  theme_bw(base_size = 11, base_family = "Arial") +
  theme(panel.grid = element_blank()) +
  geom_density(alpha=0.6, adjust=2, colour=NA) +
  geom_hline(yintercept=median_west, colour="#ff6d00") +
  geom_hline(yintercept=median_east, colour="#5c7ec0") +
  scale_fill_manual(values=c("#ff6d00","#5c7ec0"), guide="none") +
  ylim(0,5.5) +
  ylab("Genetic offset")
p_gen_off_groups_2040_mod

ggsave("~/Dropbox/Postdoc_Milan/LK_Joan/GEA/plots/local_offset_densities_groups_2040_mod_bio19-8-2-15-10.pdf", p_gen_off_groups_2040_mod, device = cairo_pdf, units="cm", width=4, height=12, limitsize=FALSE)

median_west <- median(genom_offset_presence$Global_offset_2040_ext[genom_offset_presence$group=="western"])
median_east <- median(genom_offset_presence$Global_offset_2040_ext[genom_offset_presence$group=="eastern"])

p_gen_off_groups_2040_ext <- ggplot(data=genom_offset_presence, aes(y=Global_offset_2040_ext,fill=group)) +
  theme_bw(base_size = 11, base_family = "Arial") +
  theme(panel.grid = element_blank()) +
  geom_density(alpha=0.6, adjust=2, colour=NA) +
  geom_hline(yintercept=median_west, colour="#ff6d00") +
  geom_hline(yintercept=median_east, colour="#5c7ec0") +
  scale_fill_manual(values=c("#ff6d00","#5c7ec0"), guide="none") +
  ylim(0,5.5) +
  ylab("Genetic offset")
p_gen_off_groups_2040_ext

ggsave("~/Dropbox/Postdoc_Milan/LK_Joan/GEA/plots/local_offset_densities_groups_2040_bio19-8-2-15-10.pdf", p_gen_off_groups_2040, device = cairo_pdf, units="cm", width=4, height=12, limitsize=FALSE)

median_west <- median(genom_offset_presence$Global_offset_2070_ext[genom_offset_presence$group=="western"])
median_east <- median(genom_offset_presence$Global_offset_2070_ext[genom_offset_presence$group=="eastern"])

p_gen_off_groups_2070_ext <- ggplot(data=genom_offset_presence, aes(y=Global_offset_2070_ext,fill=group)) +
  theme_bw(base_size = 11, base_family = "Arial") +
  theme(panel.grid = element_blank()) +
  geom_density(alpha=0.6, adjust=2, colour=NA) +
  geom_hline(yintercept=median_west, colour="#ff6d00") +
  geom_hline(yintercept=median_east, colour="#5c7ec0") +
  scale_fill_manual(values=c("#ff6d00","#5c7ec0"), guide="none") +
  ylim(0,5.5) +
  ylab("Genetic offset")
p_gen_off_groups_2070_ext

ggsave("~/Dropbox/Postdoc_Milan/LK_Joan/GEA/plots/local_offset_densities_groups_2070_bio19-8-2-15-10.pdf", p_gen_off_groups_2070, device = cairo_pdf, units="cm", width=4, height=12, limitsize=FALSE)

### Statistical differences between western and eastern group in genetic offsets
t.test(genom_offset_presence$Global_offset_2040[genom_offset_presence$group=="western"],
       genom_offset_presence$Global_offset_2040[genom_offset_presence$group=="eastern"]) #p-value = 4.918e-14
wilcox.test(genom_offset_presence$Global_offset_2040 ~ genom_offset_presence$group) #p-value < 2.2e-16
mean(genom_offset_presence$Global_offset_2040[genom_offset_presence$group=="western"]) #0.6014195
mean(genom_offset_presence$Global_offset_2040[genom_offset_presence$group=="eastern"]) #1.119877

t.test(genom_offset_presence$Global_offset_2070[genom_offset_presence$group=="western"],
       genom_offset_presence$Global_offset_2070[genom_offset_presence$group=="eastern"]) #p-value < 2.2e-16
wilcox.test(genom_offset_presence$Global_offset_2070 ~ genom_offset_presence$group) #p-value < 2.2e-16
mean(genom_offset_presence$Global_offset_2070[genom_offset_presence$group=="western"]) #0.6027215
mean(genom_offset_presence$Global_offset_2070[genom_offset_presence$group=="eastern"]) #1.394648
