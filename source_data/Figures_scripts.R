######################################################
## Scripts to produce figures for LK genomics paper ##
######################################################
rm(list=ls())

library(tidyverse)
library(scales)
library(raster)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggnewscale)
library(ragg)
library(patchwork)
library(svglite)

# Figure 1: sampling sites and occurrence data
## Load range data
range <- sf::st_read("species_22696357.shp") %>% # Read shapefile of Birdlife distribution. This file is not provided as source data but can be requested from Birdlife International (https://datazone.birdlife.org/contact-us/request-our-data)
  filter(legend == "Extant (non breeding)" | legend == "Extant (breeding)") %>% # Filter breeding and non-breeding distributions
  mutate_if(is.character, as.factor) # Convert character variables to factors

## Load sampling data
sampling <- read.csv("summary_sampling.csv")
sampling$Acronym <- factor(sampling$Acronym,levels = c("POR","ESN","ESS","SIC","ITS","ITN","CRO","GRG","GRC","GRL","TUR","ISR","KAZ","RUS","MON","MOS"))
### Add population order
pop_order <- c("POR","ESN","ESS","SIC","ITS","ITN","CRO","GRG","GRC","GRL","TUR","ISR","KAZ","RUS","MON","MOS")
### Add colours
palette<-(c("#a80000","#db0000","#dc2906","#dd4a11","#e16d19","#e68f21","#edb42a","#f5d933","#f5d933",
            "#e4e940","#cad34f","#b0bc66","#98a881","#7f929e","#677dbc","#5068dc"))
names(palette)<-pop_order

## Load occurrence data
occurrence <- read.csv("occurrence_data.csv")
occurrence$ESU_season <- paste(occurrence$ESU, occurrence$season, sep="_")

## Load country polygons
admin <- ne_countries(scale = "medium", returnclass = "sf")

## Plot sampling
p_sampling <- ggplot() + 
  theme_bw(base_family = "Arial") +
  theme(panel.border = element_rect(color = "black", fill = NA, size=1),
        panel.grid = element_blank(),
        legend.position='none',
        axis.text = element_text(color="black", size=12),
        axis.title = element_text(size=14)) +
  geom_sf(data = range, aes(fill=legend), lwd=0, alpha=1) +
  geom_sf(data = admin, fill=NA, lwd=0.25, colour=gray(.5)) +
  scale_fill_manual(values=c("#C5AA83","#859EC6"), guide="none") +
  scale_color_manual(values=c("#C5AA83","#859EC6"), guide="none") +
  new_scale_fill() +
  geom_point(shape=21, data = sampling, aes(x = Longitude, y = Latitude, fill = Acronym), colour="black", size=3.5) +
  scale_fill_manual(values=palette, guide="none") +
  coord_sf(xlim = c(150, -20), ylim = c(25, 60), expand = F) +
  xlab("Longitude") + ylab("Latitude")
p_sampling

## Plot occurrence
p_occurrence <- ggplot() + 
  theme_bw(base_family = "Arial") +
  theme(panel.border = element_rect(color = "black", fill = NA, size=1),
        panel.grid = element_blank(),
        legend.position='none',
        axis.text = element_text(color="black", size=12),
        axis.title = element_text(size=14)) +
  geom_sf(data = range, aes(fill=legend), lwd=0, alpha=1) +
  geom_sf(data = admin, fill=NA, lwd=0.25, colour=gray(.5)) +
  scale_fill_manual(values=c("#C5AA83","#859EC6"), guide="none") +
  scale_color_manual(values=c("#C5AA83","#859EC6"), guide="none") +
  new_scale_fill() +
  geom_point(shape=21, data = occurrence, aes(x = longitude, y = latitude, fill = ESU_season), colour="black") +
  scale_fill_manual(values=c("#5c7ec0","#BDCBE5","#ff6d00","#ffc499"), guide="none") +
  coord_sf(xlim = c(150, -20), ylim = c(-38, 65), expand = F) +
  xlab("Longitude") + ylab("Latitude")
p_occurrence

## Arrange plots and save
p_fig1 <- p_sampling / p_occurrence + plot_layout(guides = "collect")
ggsave("sampling_sites_occurrence_data.svg", p_fig1, device="svg", units="cm", width=20, height=20, limitsize=FALSE)
ggsave("sampling_sites_occurrence_data.pdf", p_fig1, device = ragg::agg_pdf, width = 20, height = 20, units = "cm", limitsize = FALSE)

# Figure 2: Population structure and gene flow
## Figure 2a: Plot admixture results: script based on https://luisdva.github.io/rstats/model-cluster-plots/ 

### Load libs
lib<-c("ggplot2","gridExtra","grid","dplyr","stringi","forcats","ggrepel","purrr","cowplot","tidyr","plyr")
lapply(lib,library,character.only=T)

### Read in individual to population map and add populations 
inds <- read.table("LK_intersect_norelated_wCro_nosexchrom.fam",header=F)
inds <- inds$V1
inds <- as.data.frame(inds)
colnames(inds)[1] <- "Sample"
inds$Populations <- as.character(c("ESN","ESN","ESN","GRC","GRC","ESN","GRC","GRC","ISR",
                                   "TUR","ESN","KAZ","ISR","ISR","TUR","TUR","CRO","ISR","CRO","ESN",
                                   "SIC","GRC","RUS","RUS","RUS","MOS","MOS","MON","ISR",
                                   "ISR","ESS","ESS","ESS","ESN","ESN","ITN","ITN","ITN",
                                   "SIC","RUS","MON","MON","ITN","ITS","ITS","ITS","ITS",
                                   "SIC","RUS","ISR","ISR","ITS","ITS","SIC","KAZ","RUS",
                                   "MON","CRO","ESS","ESS","ESS","ESS","ITN","ITN","ITN","ITS",
                                   "ITS","GRL","GRL","GRL","GRG","GRG","KAZ","MOS","MOS",
                                   "MON","MON","SIC","SIC","GRL","GRL","GRL","GRG","KAZ"))
### Add population order
pop_order <- c("ESN","ESS","SIC","ITS","ITN","CRO","GRG","GRC","GRL","TUR","ISR","KAZ","RUS","MON","MOS")
inds$Populations <- factor(inds$Populations, levels=pop_order)

### Read in K values
K2 <- read.table("LK_intersect_norelated_wCro_nosexchrom.2.Q",col.names=c("1","2"),check.names = F)

### Add sample and population information and select K values
K2$sampleID <- inds$Sample
K2$loc <- inds$Populations
K2 <- K2 %>% gather(key=popGroup,value=prob,1:2) %>%
  dplyr::select(sampleID,popGroup,prob,loc)
K2$loc <- factor(K2$loc, levels=pop_order)

### Set up for plotting
palette1 <- c("#ff6d00","#5c7ec0")

### plot:
K2_plot <- ggplot(data=K2, aes(factor(sampleID), prob, fill = factor(popGroup), colour = factor(popGroup))) +
  geom_col(size = 1) +
  ylab(expression(italic(K)~"= 2"))+
  facet_grid(~loc, switch = "x", scales = "free", space = "free") +
  theme_minimal() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  scale_colour_manual(values=palette1)+
  scale_fill_manual(values=palette1)+
  theme(panel.spacing.x = unit(0.1, "lines"),
        axis.text.x = element_blank(),
        axis.text.y=element_text(family="Arial",size=14),
        axis.title.x=element_blank(),
        axis.title.y=element_text(family="Arial",size=16),
        panel.grid = element_blank(),
        legend.position="none",
        strip.background = element_blank(),
        strip.text.x=element_text(family="Arial",angle=90,size=12))

K2_plot

ggsave("admixture_K2.pdf", K2_plot, width = 20, height = 5, device=cairo_pdf, units="cm", limitsize = FALSE)

## Figure 2b: Plot PCA
### Read eigenvectors
eigenvec_table <- read.table('LK_intersect_norelated_wCro_nosexchrom.eigenvec', header = FALSE)

### Add populations
eigenvec_table$Populations <- as.character(c("ESN","ESN","ESN","GRC","GRC","ESN","GRC","GRC","ISR",
                                             "TUR","ESN","KAZ","ISR","ISR","TUR","TUR","CRO","ISR","CRO","ESN",
                                             "SIC","GRC","RUS","RUS","RUS","MOS","MOS","MON","ISR",
                                             "ISR","ESS","ESS","ESS","ESN","ESN","ITN","ITN","ITN",
                                             "SIC","RUS","MON","MON","ITN","ITS","ITS","ITS","ITS",
                                             "SIC","RUS","ISR","ISR","ITS","ITS","SIC","KAZ","RUS",
                                             "MON","CRO","ESS","ESS","ESS","ESS","ITN","ITN","ITN","ITS",
                                             "ITS","GRL","GRL","GRL","GRG","GRG","KAZ","MOS","MOS",
                                             "MON","MON","SIC","SIC","GRL","GRL","GRL","GRG","KAZ"))
### Add population order
pop_order <- c("ESN","ESS","SIC","ITS","ITN","CRO","GRG","GRC","GRL","TUR","ISR","KAZ","RUS","MON","MOS")

### Add colours
palette<-(c("#db0000","#dc2906","#dd4a11","#e16d19","#e68f21","#edb42a","#f5d933","#f5d933",
            "#e4e940","#cad34f","#b0bc66","#98a881","#7f929e","#677dbc","#5068dc"))
names(palette)<-pop_order

label <- eigenvec_table$V1
eigenvec_table <- eigenvec_table[-2]
eigenvec_table <- eigenvec_table[-1]

head(eigenvec_table)
for (i in 1:10){
  colnames(eigenvec_table)[i]<-paste0("PC",i)
}
eigenval <- read.table('LK_intersect_norelated_wCro_nosexchrom.eigenval', header = F)
percentage <- round(eigenval$V1/sum(eigenval$V1)*100,2)
percentage <- paste0(colnames(eigenvec_table)[1:10]," (",paste(as.character(percentage),"%)"))
percentage
eigenvec_table$Populations <- factor(eigenvec_table$Populations, levels=pop_order)

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

## Figure 2c: Plot EEMS
### Load libraries
library(rEEMSplots)
library(rgdal)
library(rworldmap)
library(rworldxtra)

### Load results
eems_results <- file.path("eems_results")

### Plot results
eems.plots(mcmcpath=eems_results, plotpath=paste0(eems_results,"final"),
           projection.in = "+proj=longlat +datum=WGS84",
           projection.out = "+proj=merc +datum=WGS84",
           longlat=T, out.png = F, add.grid = F, add.outline = T, lwd.outline = 1, add.demes = T,
           add.map=T, lwd.map=0.5,min.cex.demes = 0.5, max.cex.demes = 1.5)

eems_polygon <- read.delim("LK_minconvexpoly.txt")
world <- ne_countries(scale = "medium", returnclass = "sf")
ggplot(data = world) +
  geom_sf() +
  geom_point(data = eems_polygon, aes(x = longitude, y = latitude), size = 4, 
             shape = 23, fill = "darkred") +
  coord_sf(xlim = c(-10, 130), ylim = c(25, 65), expand = FALSE)

## Figure 2d: Plot Treemix population tree
source("treemix_data/plotting_funcs_treemix.R")
plot_tree("treemix_data/final.0.0")

## Figure 2e: Plot Treemix residuals
plot_resid("treemix_data/final.0.0", "treemix_data/poporder_treemix.txt")

# Figure 3: Ecological differentiation
## Figure 3a: Plot current SDMs

### Load stacks from SDMs
ras_east_breed_past <- stack("SDM_rasters_current_past/stack_breeding_east_chelsa_past_recl.grd")
ras_west_breed_past <- stack("SDM_rasters_current_past/stack_breeding_west_chelsa_past_recl.grd")
ras_east_wint_past <- stack("SDM_rasters_current_past/stack_wintering_east_chelsa_past_recl.grd")
ras_west_wint_past <- stack("SDM_rasters_current_past/stack_wintering_west_chelsa_past_recl.grd")

# Vectorization of the climatic rasters for ggplot
ras_east_breed_past <- as(ras_east_breed_past, "SpatialPixelsDataFrame")
ras_east_breed_past <- as.data.frame(ras_east_breed_past)
ras_west_breed_past <- as(ras_west_breed_past, "SpatialPixelsDataFrame")
ras_west_breed_past <- as.data.frame(ras_west_breed_past)
ras_east_wint_past <- as(ras_east_wint_past, "SpatialPixelsDataFrame")
ras_east_wint_past <- as.data.frame(ras_east_wint_past)
ras_west_wint_past <- as(ras_west_wint_past, "SpatialPixelsDataFrame")
ras_west_wint_past <- as.data.frame(ras_west_wint_past)

### Plot and save
p_current_SDM <- ggplot() + 
  geom_sf(data = admin, fill="black", lwd=0, colour="black") +
  geom_tile(data = subset(ras_east_breed_past, p0_breeding_east_prediction_crop==1), aes(x = x, y = y, fill = MH_breeding_est_prediction), fill = "#5c7ec0") + 
  geom_tile(data = subset(ras_east_wint_past, p0_wintering_east_prediction_crop==1), aes(x = x, y = y, fill = MH_wintering_est_prediction), fill = "#BDCBE5") + 
  geom_tile(data = subset(ras_west_breed_past, p0_breeding_west_prediction_crop==1), aes(x = x, y = y, fill = MH_breeding_ovest_prediction), fill = "#ff6d00") + 
  geom_tile(data = subset(ras_west_wint_past, p0_wintering_west_prediction_crop==1), aes(x = x, y = y, fill = MH_wintering_west_prediction), fill = "#ffc499") +
  coord_sf(xlim = c(150, -20), ylim = c(-38, 65), expand = F) +
  theme_void()
p_current_SDM

ggsave("current_SDM.pdf", p_current_SDM, device="pdf", units="cm", width=50, height=35, limitsize=FALSE)

## Figure 3b: Plot habitat use and selection
### See script https://github.com/jferrerobiol/lk_climate/blob/main/ecological_differentiation/habitat_use_and_preference.R

## Figure 3c and d: Plot climatic niche
### See script https://github.com/jferrerobiol/lk_climate/blob/main/ecological_differentiation/ecospat_niche_modelling.R

# Figure 4: Climate-associated genetic variation
## Figure 4a: IBD and IBC
### See script https://github.com/jferrerobiol/lk_climate/blob/main/climate_associated_genetic_variation/IBD_IBC.R

## Figure 4b: PCA climate-associated SNPs
### See script https://github.com/jferrerobiol/lk_climate/blob/main/climate_associated_genetic_variation/GEA_genetic_offsets.R

## Figure 4c: Gradient forest
### See script https://github.com/jferrerobiol/lk_climate/blob/main/climate_associated_genetic_variation/GEA_genetic_offsets.R

## Figure 4d: Heatmap climate-associated SNPs
### See script https://github.com/jferrerobiol/lk_climate/blob/main/climate_associated_genetic_variation/GEA_genetic_offsets.R

# Figure 5: Past impacts climate
## Figure 5a-f: Past fluctuations

### Load TºC data
temp <- read.table("epica_dome_c_temp_800kyr.txt", header=TRUE)
temp$Age <- temp$Age + 71

mysmoothmodel <- stats::loess(Temperature ~ Age, data = temp, span=0.05)
predict(mysmoothmodel, newdata = data.frame(Age = c(130000,20000,18000,16000,14000,12000,10000,8000,6000,4000,2000))) # I use this script to extract TºC anomalies

### Define mu per generation and site, generation time and read main results and bootstraps
mu <- 3.3e-9
gen <- 2 # used in BSP and MSMC2

### Load MSMC2 data
lk <- read.table("Eur_6hap.final.txt", header=TRUE)

boot_1 <- read.table("Eur_3ind_bootstrap_1/Eur_6hap.final.txt", header=TRUE)
boot_2 <- read.table("Eur_3ind_bootstrap_2/Eur_6hap.final.txt", header=TRUE)
boot_3 <- read.table("Eur_3ind_bootstrap_3/Eur_6hap.final.txt", header=TRUE)
boot_4 <- read.table("Eur_3ind_bootstrap_4/Eur_6hap.final.txt", header=TRUE)
boot_5 <- read.table("Eur_3ind_bootstrap_5/Eur_6hap.final.txt", header=TRUE)
boot_6 <- read.table("Eur_3ind_bootstrap_6/Eur_6hap.final.txt", header=TRUE)
boot_7 <- read.table("Eur_3ind_bootstrap_7/Eur_6hap.final.txt", header=TRUE)
boot_8 <- read.table("Eur_3ind_bootstrap_8/Eur_6hap.final.txt", header=TRUE)
boot_9 <- read.table("Eur_3ind_bootstrap_9/Eur_6hap.final.txt", header=TRUE)
boot_10 <- read.table("Eur_3ind_bootstrap_10/Eur_6hap.final.txt", header=TRUE)
boot_11 <- read.table("Eur_3ind_bootstrap_11/Eur_6hap.final.txt", header=TRUE)
boot_12 <- read.table("Eur_3ind_bootstrap_12/Eur_6hap.final.txt", header=TRUE)
boot_13 <- read.table("Eur_3ind_bootstrap_13/Eur_6hap.final.txt", header=TRUE)
boot_14 <- read.table("Eur_3ind_bootstrap_14/Eur_6hap.final.txt", header=TRUE)
boot_15 <- read.table("Eur_3ind_bootstrap_15/Eur_6hap.final.txt", header=TRUE)
boot_16 <- read.table("Eur_3ind_bootstrap_16/Eur_6hap.final.txt", header=TRUE)
boot_17 <- read.table("Eur_3ind_bootstrap_17/Eur_6hap.final.txt", header=TRUE)
boot_18 <- read.table("Eur_3ind_bootstrap_18/Eur_6hap.final.txt", header=TRUE)
boot_19 <- read.table("Eur_3ind_bootstrap_19/Eur_6hap.final.txt", header=TRUE)
boot_20 <- read.table("Eur_3ind_bootstrap_20/Eur_6hap.final.txt", header=TRUE)
boot_21 <- read.table("Eur_3ind_bootstrap_21/Eur_6hap.final.txt", header=TRUE)
boot_22 <- read.table("Eur_3ind_bootstrap_22/Eur_6hap.final.txt", header=TRUE)
boot_23 <- read.table("Eur_3ind_bootstrap_23/Eur_6hap.final.txt", header=TRUE)
boot_24 <- read.table("Eur_3ind_bootstrap_24/Eur_6hap.final.txt", header=TRUE)
boot_25 <- read.table("Eur_3ind_bootstrap_25/Eur_6hap.final.txt", header=TRUE)
boot_26 <- read.table("Eur_3ind_bootstrap_26/Eur_6hap.final.txt", header=TRUE)
boot_27 <- read.table("Eur_3ind_bootstrap_27/Eur_6hap.final.txt", header=TRUE)
boot_28 <- read.table("Eur_3ind_bootstrap_28/Eur_6hap.final.txt", header=TRUE)
boot_30 <- read.table("Eur_3ind_bootstrap_30/Eur_6hap.final.txt", header=TRUE)
boot_31 <- read.table("Eur_3ind_bootstrap_31/Eur_6hap.final.txt", header=TRUE)
boot_32 <- read.table("Eur_3ind_bootstrap_32/Eur_6hap.final.txt", header=TRUE)
boot_33 <- read.table("Eur_3ind_bootstrap_33/Eur_6hap.final.txt", header=TRUE)
boot_34 <- read.table("Eur_3ind_bootstrap_34/Eur_6hap.final.txt", header=TRUE)
boot_35 <- read.table("Eur_3ind_bootstrap_35/Eur_6hap.final.txt", header=TRUE)
boot_36 <- read.table("Eur_3ind_bootstrap_36/Eur_6hap.final.txt", header=TRUE)
boot_38 <- read.table("Eur_3ind_bootstrap_38/Eur_6hap.final.txt", header=TRUE)
boot_39 <- read.table("Eur_3ind_bootstrap_39/Eur_6hap.final.txt", header=TRUE)
boot_40 <- read.table("Eur_3ind_bootstrap_40/Eur_6hap.final.txt", header=TRUE)
boot_41 <- read.table("Eur_3ind_bootstrap_41/Eur_6hap.final.txt", header=TRUE)
boot_42 <- read.table("Eur_3ind_bootstrap_42/Eur_6hap.final.txt", header=TRUE)
boot_43 <- read.table("Eur_3ind_bootstrap_43/Eur_6hap.final.txt", header=TRUE)
boot_44 <- read.table("Eur_3ind_bootstrap_44/Eur_6hap.final.txt", header=TRUE)
boot_45 <- read.table("Eur_3ind_bootstrap_45/Eur_6hap.final.txt", header=TRUE)
boot_46 <- read.table("Eur_3ind_bootstrap_46/Eur_6hap.final.txt", header=TRUE)
boot_47 <- read.table("Eur_3ind_bootstrap_47/Eur_6hap.final.txt", header=TRUE)
boot_48 <- read.table("Eur_3ind_bootstrap_48/Eur_6hap.final.txt", header=TRUE)
boot_49 <- read.table("Eur_3ind_bootstrap_49/Eur_6hap.final.txt", header=TRUE)
boot_50 <- read.table("Eur_3ind_bootstrap_50/Eur_6hap.final.txt", header=TRUE)
boot_51 <- read.table("Eur_3ind_bootstrap_51/Eur_6hap.final.txt", header=TRUE)
boot_52 <- read.table("Eur_3ind_bootstrap_52/Eur_6hap.final.txt", header=TRUE)
boot_53 <- read.table("Eur_3ind_bootstrap_53/Eur_6hap.final.txt", header=TRUE)
boot_54 <- read.table("Eur_3ind_bootstrap_54/Eur_6hap.final.txt", header=TRUE)
boot_55 <- read.table("Eur_3ind_bootstrap_55/Eur_6hap.final.txt", header=TRUE)
boot_56 <- read.table("Eur_3ind_bootstrap_56/Eur_6hap.final.txt", header=TRUE)
boot_57 <- read.table("Eur_3ind_bootstrap_57/Eur_6hap.final.txt", header=TRUE)
boot_58 <- read.table("Eur_3ind_bootstrap_58/Eur_6hap.final.txt", header=TRUE)
boot_59 <- read.table("Eur_3ind_bootstrap_59/Eur_6hap.final.txt", header=TRUE)
boot_60 <- read.table("Eur_3ind_bootstrap_60/Eur_6hap.final.txt", header=TRUE)
boot_61 <- read.table("Eur_3ind_bootstrap_61/Eur_6hap.final.txt", header=TRUE)
boot_62 <- read.table("Eur_3ind_bootstrap_62/Eur_6hap.final.txt", header=TRUE)
boot_63 <- read.table("Eur_3ind_bootstrap_63/Eur_6hap.final.txt", header=TRUE)
boot_64 <- read.table("Eur_3ind_bootstrap_64/Eur_6hap.final.txt", header=TRUE)
boot_65 <- read.table("Eur_3ind_bootstrap_65/Eur_6hap.final.txt", header=TRUE)
boot_66 <- read.table("Eur_3ind_bootstrap_66/Eur_6hap.final.txt", header=TRUE)
boot_67 <- read.table("Eur_3ind_bootstrap_67/Eur_6hap.final.txt", header=TRUE)
boot_68 <- read.table("Eur_3ind_bootstrap_68/Eur_6hap.final.txt", header=TRUE)
boot_69 <- read.table("Eur_3ind_bootstrap_69/Eur_6hap.final.txt", header=TRUE)
boot_70 <- read.table("Eur_3ind_bootstrap_70/Eur_6hap.final.txt", header=TRUE)
boot_71 <- read.table("Eur_3ind_bootstrap_71/Eur_6hap.final.txt", header=TRUE)
boot_72 <- read.table("Eur_3ind_bootstrap_72/Eur_6hap.final.txt", header=TRUE)
boot_73 <- read.table("Eur_3ind_bootstrap_73/Eur_6hap.final.txt", header=TRUE)
boot_74 <- read.table("Eur_3ind_bootstrap_74/Eur_6hap.final.txt", header=TRUE)
boot_75 <- read.table("Eur_3ind_bootstrap_75/Eur_6hap.final.txt", header=TRUE)
boot_76 <- read.table("Eur_3ind_bootstrap_76/Eur_6hap.final.txt", header=TRUE)
boot_77 <- read.table("Eur_3ind_bootstrap_77/Eur_6hap.final.txt", header=TRUE)
boot_78 <- read.table("Eur_3ind_bootstrap_78/Eur_6hap.final.txt", header=TRUE)
boot_79 <- read.table("Eur_3ind_bootstrap_79/Eur_6hap.final.txt", header=TRUE)
boot_80 <- read.table("Eur_3ind_bootstrap_80/Eur_6hap.final.txt", header=TRUE)
boot_81 <- read.table("Eur_3ind_bootstrap_81/Eur_6hap.final.txt", header=TRUE)
boot_82 <- read.table("Eur_3ind_bootstrap_82/Eur_6hap.final.txt", header=TRUE)
boot_83 <- read.table("Eur_3ind_bootstrap_83/Eur_6hap.final.txt", header=TRUE)
boot_84 <- read.table("Eur_3ind_bootstrap_84/Eur_6hap.final.txt", header=TRUE)
boot_85 <- read.table("Eur_3ind_bootstrap_85/Eur_6hap.final.txt", header=TRUE)
boot_86 <- read.table("Eur_3ind_bootstrap_86/Eur_6hap.final.txt", header=TRUE)
boot_87 <- read.table("Eur_3ind_bootstrap_87/Eur_6hap.final.txt", header=TRUE)
boot_88 <- read.table("Eur_3ind_bootstrap_88/Eur_6hap.final.txt", header=TRUE)
boot_89 <- read.table("Eur_3ind_bootstrap_89/Eur_6hap.final.txt", header=TRUE)
boot_90 <- read.table("Eur_3ind_bootstrap_90/Eur_6hap.final.txt", header=TRUE)
boot_91 <- read.table("Eur_3ind_bootstrap_91/Eur_6hap.final.txt", header=TRUE)
boot_92 <- read.table("Eur_3ind_bootstrap_92/Eur_6hap.final.txt", header=TRUE)
boot_93 <- read.table("Eur_3ind_bootstrap_93/Eur_6hap.final.txt", header=TRUE)
boot_94 <- read.table("Eur_3ind_bootstrap_94/Eur_6hap.final.txt", header=TRUE)
boot_95 <- read.table("Eur_3ind_bootstrap_95/Eur_6hap.final.txt", header=TRUE)
boot_96 <- read.table("Eur_3ind_bootstrap_96/Eur_6hap.final.txt", header=TRUE)
boot_97 <- read.table("Eur_3ind_bootstrap_97/Eur_6hap.final.txt", header=TRUE)
boot_98 <- read.table("Eur_3ind_bootstrap_98/Eur_6hap.final.txt", header=TRUE)
boot_99 <- read.table("Eur_3ind_bootstrap_99/Eur_6hap.final.txt", header=TRUE)

### Merge datasets and plot with ggplot2
lk$analysis <- "main"
boot_1$analysis <- "boot_1"
boot_2$analysis <- "boot_2"
boot_3$analysis <- "boot_3"
boot_4$analysis <- "boot_4"
boot_5$analysis <- "boot_5"
boot_6$analysis <- "boot_6"
boot_7$analysis <- "boot_7"
boot_8$analysis <- "boot_8"
boot_9$analysis <- "boot_9"
boot_10$analysis <- "boot_10"
boot_11$analysis <- "boot_11"
boot_12$analysis <- "boot_12"
boot_13$analysis <- "boot_13"
boot_14$analysis <- "boot_14"
boot_15$analysis <- "boot_15"
boot_16$analysis <- "boot_16"
boot_17$analysis <- "boot_17"
boot_18$analysis <- "boot_18"
boot_19$analysis <- "boot_19"
boot_20$analysis <- "boot_20"
boot_21$analysis <- "boot_21"
boot_22$analysis <- "boot_22"
boot_23$analysis <- "boot_23"
boot_24$analysis <- "boot_24"
boot_25$analysis <- "boot_25"
boot_26$analysis <- "boot_26"
boot_27$analysis <- "boot_27"
boot_28$analysis <- "boot_28"
boot_30$analysis <- "boot_30"
boot_31$analysis <- "boot_31"
boot_32$analysis <- "boot_32"
boot_33$analysis <- "boot_33"
boot_34$analysis <- "boot_34"
boot_35$analysis <- "boot_35"
boot_36$analysis <- "boot_36"
boot_38$analysis <- "boot_38"
boot_39$analysis <- "boot_39"
boot_40$analysis <- "boot_40"
boot_41$analysis <- "boot_41"
boot_42$analysis <- "boot_42"
boot_43$analysis <- "boot_43"
boot_44$analysis <- "boot_44"
boot_45$analysis <- "boot_45"
boot_46$analysis <- "boot_46"
boot_47$analysis <- "boot_47"
boot_48$analysis <- "boot_48"
boot_49$analysis <- "boot_49"
boot_50$analysis <- "boot_50"
boot_51$analysis <- "boot_51"
boot_52$analysis <- "boot_51"
boot_53$analysis <- "boot_53"
boot_54$analysis <- "boot_54"
boot_55$analysis <- "boot_55"
boot_56$analysis <- "boot_56"
boot_57$analysis <- "boot_57"
boot_58$analysis <- "boot_58"
boot_59$analysis <- "boot_59"
boot_60$analysis <- "boot_60"
boot_61$analysis <- "boot_61"
boot_62$analysis <- "boot_62"
boot_63$analysis <- "boot_63"
boot_64$analysis <- "boot_64"
boot_65$analysis <- "boot_65"
boot_66$analysis <- "boot_66"
boot_67$analysis <- "boot_67"
boot_68$analysis <- "boot_68"
boot_69$analysis <- "boot_69"
boot_70$analysis <- "boot_70"
boot_71$analysis <- "boot_71"
boot_72$analysis <- "boot_72"
boot_73$analysis <- "boot_73"
boot_74$analysis <- "boot_74"
boot_75$analysis <- "boot_75"
boot_76$analysis <- "boot_76"
boot_77$analysis <- "boot_77"
boot_78$analysis <- "boot_78"
boot_79$analysis <- "boot_79"
boot_80$analysis <- "boot_80"
boot_81$analysis <- "boot_81"
boot_82$analysis <- "boot_82"
boot_83$analysis <- "boot_83"
boot_84$analysis <- "boot_84"
boot_85$analysis <- "boot_85"
boot_86$analysis <- "boot_86"
boot_87$analysis <- "boot_87"
boot_88$analysis <- "boot_88"
boot_89$analysis <- "boot_89"
boot_90$analysis <- "boot_90"
boot_91$analysis <- "boot_91"
boot_92$analysis <- "boot_92"
boot_93$analysis <- "boot_93"
boot_94$analysis <- "boot_94"
boot_95$analysis <- "boot_95"
boot_96$analysis <- "boot_96"
boot_97$analysis <- "boot_97"
boot_98$analysis <- "boot_98"
boot_99$analysis <- "boot_99"

msmc <- bind_rows(lk,boot_1,boot_2,boot_3,boot_4,boot_5,boot_6,
                  boot_7,boot_8,boot_9,boot_10,boot_11,boot_12,boot_13,
                  boot_14,boot_15,boot_16,boot_17,boot_18,boot_19,
                  boot_20,boot_21,boot_22,boot_23,boot_24,boot_25,
                  boot_26,boot_27,boot_28,boot_30,boot_31,boot_32,
                  boot_33,boot_34,boot_35,boot_36,boot_38,boot_39,
                  boot_40,boot_41,boot_42,boot_43,boot_44,boot_45,
                  boot_46,boot_47,boot_48,boot_49,boot_50,boot_51,
                  boot_52,boot_53,boot_54,boot_55,boot_56,boot_57,
                  boot_58,boot_59,boot_60,boot_61,boot_62,boot_63,
                  boot_64,boot_65,boot_66,boot_67,boot_68,boot_69,
                  boot_70,boot_71,boot_72,boot_73,boot_74,boot_75,
                  boot_76,boot_77,boot_78,boot_79,boot_80,boot_81,
                  boot_82,boot_83,boot_84,boot_85,boot_86,boot_87,
                  boot_88,boot_89,boot_90,boot_91,boot_92,boot_93,
                  boot_94,boot_95,boot_96,boot_97,boot_98,boot_99)

msmc$analysis <- factor(msmc$analysis, levels=c("main","boot_1","boot_2","boot_3","boot_4",
                                                "boot_5","boot_6","boot_7","boot_8","boot_9",
                                                "boot_10","boot_11","boot_12","boot_13","boot_14",
                                                "boot_15","boot_16","boot_17","boot_18","boot_19",
                                                "boot_20","boot_21","boot_22","boot_23","boot_24",
                                                "boot_25","boot_26","boot_27","boot_28","boot_30",
                                                "boot_31","boot_32","boot_33","boot_34","boot_35",
                                                "boot_36","boot_38","boot_39","boot_40","boot_41",
                                                "boot_42","boot_43","boot_44","boot_45","boot_46",
                                                "boot_47","boot_48","boot_49","boot_50","boot_51",
                                                "boot_52","boot_53","boot_54","boot_55","boot_56",
                                                "boot_57","boot_58","boot_59","boot_60","boot_61",
                                                "boot_62","boot_63","boot_64","boot_65","boot_66",
                                                "boot_67","boot_68","boot_69","boot_70","boot_71",
                                                "boot_72","boot_73","boot_74","boot_75","boot_76",
                                                "boot_77","boot_78","boot_79","boot_80","boot_81",
                                                "boot_82","boot_83","boot_84","boot_85","boot_86",
                                                "boot_87","boot_88","boot_89","boot_90","boot_91",
                                                "boot_92","boot_93","boot_94","boot_95","boot_96",
                                                "boot_97","boot_98","boot_99"))

analysis_alpha <- c(1,rep(0.06, 97))

msmc_main <- msmc[msmc$analysis=="main",]
data_demo_msmc <- data.frame(time=msmc_main$left_time_boundary/mu*gen, Ne=(1/msmc_main$lambda)/(2*mu)) # I get MSMC2 Ne values from here

### Load BSP data
bsp_EUR <- read.delim("BSP_EUR_samples.txt")
bsp_ASI <- read.delim("BSP_ASI_samples.txt")

### Load suitable habitat data
suit <- read_excel("../corr_demo_suithab_temp/demo_suithab_temp.xlsx")

### Plot TºC data
temp_plot <- ggplot(data=temp, aes(x=Age, y=Temperature)) +
  theme_bw(base_family = "Arial") +
  theme(panel.grid = element_blank(),
        legend.position='none',
        axis.text = element_text(color="black", size=14),
        axis.title = element_text(size=16)) +
  scale_x_log10(limits=c(10^3.2,10^5.2),
                breaks = c(10000,100000),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_continuous(limits = c(-11,5.5), expand = c(0,0), breaks = seq(-9,4,3)) +
  geom_rect(aes(xmin=19000,xmax=33000,ymin=-Inf,ymax=Inf),fill="#f2f2f2") +
  geom_vline(xintercept=c(seq(2000,18000,4000),130000), col="black", size=0.75, alpha=0.5) +
  geom_line(col="black", size=0.75, alpha=0.5) +
  geom_smooth(method="loess", span=0.05, se=FALSE, col="black") +
  xlab("") + 
  ylab("Temperature anomaly (ºC)") +
  annotation_logticks(sides="b") 

temp_plot

### Plot dyabc best divergence scenario
dyabc_plot <- ggplot() +
  theme_bw(base_family = "Arial") +
  theme(panel.grid = element_blank(),
        legend.position='none',
        axis.text = element_text(color="black", size=14),
        axis.title = element_text(size=16)) +
  scale_x_log10(limits=c(10^3.2,10^5.2),
                breaks = c(10000,100000),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_continuous(limits = c(0,4), expand = c(0,0)) +
  geom_rect(aes(xmin=19000,xmax=33000,ymin=-Inf,ymax=Inf),fill="#f2f2f2") +
  geom_vline(xintercept=c(seq(2000,18000,4000),130000), col="black", size=0.75, alpha=0.5) +
  geom_segment(aes(x=40400,xend=130000,y=2,yend=2),col="black",size=1) +
  geom_segment(aes(x=0,xend=40400,y=3.5,yend=2),col="#D04335",size=1) +
  geom_segment(aes(x=0,xend=40400,y=0.5,yend=2),col="#4671B1",size=1) +
  geom_segment(aes(x=0,xend=12100,y=1.5,yend=2),col="#B4BB73",size=1) +
  geom_segment(aes(x=0,xend=4900,y=2.5,yend=2.35),col="#EBDC5C",size=1) +
  geom_segment(aes(x=4900,xend=4900,y=2.9,yend=2.35),col="#D04335",size=1) +
  geom_segment(aes(x=4900,xend=4900,y=1.8,yend=2.35),col="#B4BB73",size=1) +
  geom_segment(aes(x=12100,xend=12100,y=2.52,yend=2),col="#D04335",size=1) +
  geom_segment(aes(x=12100,xend=12100,y=1.48,yend=2),col="#4671B1",size=1) +
  geom_segment(aes(x=10^3.2,xend=9000,y=3.7,yend=3.7),col="#c67149",size=0.5) +
  geom_segment(aes(x=4000,xend=23100,y=3.2,yend=3.2),col="#c67149",size=0.5) +
  geom_segment(aes(x=15800,xend=68700,y=2.7,yend=2.7),col="#c67149",size=0.5) +
  geom_point(aes(x=4900,y=3.7),col="#c67149",size=3) +
  geom_point(aes(x=12100,y=3.2),col="#c67149",size=3) +
  geom_point(aes(x=40400,y=2.7),col="#c67149",size=3) +
  xlab("") + 
  ylab("") +
  annotation_logticks(sides="b") 

dyabc_plot

### Plot MSMC2 demographic trends
msmc_plot <- ggplot() +
  theme_bw(base_family = "Arial") +
  theme(panel.grid = element_blank(),
        legend.position='none',
        axis.text = element_text(color="black", size=14),
        axis.title = element_text(size=16)) +
  scale_x_log10(limits=c(10^3.2,10^5.2),
                breaks = c(10000,100000),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_log10(limits=c(10^5.3,10^6.5),
                breaks = c(10000,100000,1000000),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  geom_rect(aes(xmin=19000,xmax=33000,ymin=-Inf,ymax=Inf),fill="#f2f2f2") +
  geom_vline(xintercept=c(seq(2000,18000,4000),130000), col="black", size=0.75, alpha=0.5) +
  geom_step(data=msmc[msmc$left_time_boundary/mu*gen>3000,], aes(x=left_time_boundary/mu*gen, y=(1/lambda)/(2*mu), alpha=analysis), size=0.75, col="#ff6d00") +
  scale_alpha_manual(values=analysis_alpha, guide="none") +
  xlab("") +
  ylab(expression(paste(italic("N")["e"]," MSMC2"))) +
  annotation_logticks(sides="bl") 

msmc_plot

### Plot BSP demographic trends
bsp_plot <- ggplot() +
  theme_bw(base_family = "Arial") +
  theme(panel.grid = element_blank(),
        legend.position='none',
        axis.text = element_text(color="black", size=14),
        axis.title = element_text(size=16)) +
  scale_x_log10(limits=c(10^3.2,10^5.2),
                breaks = c(10000,100000),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_log10(limits=c(10^3.3,10^6),
                breaks = c(10000,100000,1000000),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  geom_rect(aes(xmin=19000,xmax=33000,ymin=-Inf,ymax=Inf),fill="#f2f2f2") +
  geom_vline(xintercept=c(seq(2000,18000,4000),130000), col="black", size=0.75, alpha=0.5) +
  geom_line(data=bsp_EUR, aes(x=Time, y=Ne), size=0.75, col="#ff6d00") +
  geom_ribbon(data=bsp_EUR, aes(x=Time, y=Ne, ymin=LowerNe, ymax=UpperNe), size=0.75, fill="#ff6d00", alpha=0.1) +
  geom_line(data=bsp_ASI, aes(x=Time, y=Ne), size=0.75, col="#5c7ec0") +
  geom_ribbon(data=bsp_ASI, aes(x=Time, y=Ne, ymin=LowerNe, ymax=UpperNe), size=0.75, fill="#5c7ec0", alpha=0.1) +
  xlab("") +
  ylab(expression(paste(italic("N")["e"]," BSP"))) +
  annotation_logticks(sides="bl") 

bsp_plot

### Plot breeding suitable habitat trends
breed_suit_plot <- ggplot() +
  theme_bw(base_family = "Arial") +
  theme(panel.grid = element_blank(),
        legend.position='none',
        axis.text = element_text(color="black", size=14),
        axis.title = element_text(size=16)) +
  scale_x_log10(limits=c(10^3.2,10^5.2),
                breaks = c(10000,100000),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_log10(limits = c(1,18), expand = c(0,0), breaks = c(2,5,10,15),
                labels = c(2,5,10,15)) +
  geom_rect(aes(xmin=19000,xmax=33000,ymin=-Inf,ymax=Inf),fill="#f2f2f2") +
  geom_vline(xintercept=c(seq(2000,18000,4000),130000), col="black", size=0.75, alpha=0.5) +
  geom_line(data=suit[suit$time<130000,], aes(x=time, y=western_breed_suit_hab/10^6), size=0.75, col="#ff6d00") +
  geom_line(data=suit[suit$time<130000,], aes(x=time, y=eastern_breed_suit_hab/10^6), size=0.75, col="#5c7ec0") +
  geom_line(data=suit[suit$time>18000,], aes(x=time, y=western_breed_suit_hab/10^6), size=0.75, col="#ff6d00", linetype=3) +
  geom_line(data=suit[suit$time>18000,], aes(x=time, y=eastern_breed_suit_hab/10^6), size=0.75, col="#5c7ec0", linetype=3) +
  xlab("") + 
  ylab(expression(paste("Suitable breeding habitat (million k", m^{2}, ")"))) +
  annotation_logticks(sides="bl") 

breed_suit_plot

### Plot non-breeding suitable habitat trends
nonbreed_suit_plot <- ggplot() +
  theme_bw(base_family = "Arial") +
  theme(panel.grid = element_blank(),
        legend.position='none',
        axis.text = element_text(color="black", size=14),
        axis.title = element_text(size=16)) +
  scale_x_log10(limits=c(10^3.2,10^5.2),
                breaks = c(10000,100000),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_log10(limits = c(0.8,6), expand = c(0,0), breaks = c(1:5),
                labels = c(1:5)) +
  geom_rect(aes(xmin=19000,xmax=33000,ymin=-Inf,ymax=Inf),fill="#f2f2f2") +
  geom_vline(xintercept=c(seq(2000,18000,4000),130000), col="black", size=0.75, alpha=0.5) +
  geom_line(data=suit[suit$time<130000,], aes(x=time, y=western_wint_suit_hab/10^6), size=0.75, col="#ff6d00") +
  geom_line(data=suit[suit$time<130000,], aes(x=time, y=eastern_wint_suit_hab/10^6), size=0.75, col="#5c7ec0") +
  geom_line(data=suit[suit$time>18000,], aes(x=time, y=western_wint_suit_hab/10^6), size=0.75, col="#ff6d00", linetype=3) +
  geom_line(data=suit[suit$time>18000,], aes(x=time, y=eastern_wint_suit_hab/10^6), size=0.75, col="#5c7ec0", linetype=3) +
  xlab("Years ago") + 
  ylab(expression(paste("Suitable non-breeding habitat (million k", m^{2}, ")"))) +
  annotation_logticks(sides="b") 

nonbreed_suit_plot

ggsave("temp_demo_suithab_split.pdf", grobz.plot, device=cairo_pdf, units="cm", width=20, height=25, limitsize=FALSE)

grobz <- lapply(list(temp_plot, dyabc_plot, msmc_plot, bsp_plot, breed_suit_plot, nonbreed_suit_plot), ggplotGrob)
grobz.plot <- arrangeGrob(grobs = list(rbind(grobz[[1]], grobz[[2]], grobz[[3]], grobz[[4]], grobz[[5]], grobz[[6]], size = "last")), ncol = 1)

ggsave("temp_abc_demo_suithab_split.pdf", grobz.plot, device=cairo_pdf, units="cm", width=20, height=42, limitsize=FALSE)

## Figure 5g: Past SDMs

### Load SDM rasters for LIG and future
ras_east_breed <- stack("stack_breeding_east_crop_recl_2071_2100.grd")
ras_east_wint <- stack("stack_wintering_east_crop_recl_2071_2100.grd")
ras_west_breed <- stack("stack_breeding_west_crop_recl_2071_2100.grd")
ras_west_wint <- stack("stack_wintering_west_crop_recl_2071_2100.grd")

#### Vectorization of the climatic rasters for ggplot
ras_east_breed <- as(ras_east_breed, "SpatialPixelsDataFrame")
ras_east_breed <- as.data.frame(ras_east_breed)
ras_east_wint <- as(ras_east_wint, "SpatialPixelsDataFrame")
ras_east_wint <- as.data.frame(ras_east_wint)
ras_west_breed <- as(ras_west_breed, "SpatialPixelsDataFrame")
ras_west_breed <- as.data.frame(ras_west_breed)
ras_west_wint <- as(ras_west_wint, "SpatialPixelsDataFrame")
ras_west_wint <- as.data.frame(ras_west_wint)

### LIG (~130 kya)
p_LIG_SDM <- ggplot() + 
  geom_sf(data = admin, fill="black", lwd=0, colour="black") +
  geom_tile(data = subset(ras_east_breed, LIG_breeding_est_prediction==1), aes(x = x, y = y, fill = LIG_breeding_est_prediction), fill = "#5c7ec0") + 
  geom_tile(data = subset(ras_east_wint, LIG_wintering_est_prediction==1), aes(x = x, y = y, fill = LIG_wintering_est_prediction), fill = "#BDCBE5") + 
  geom_tile(data = subset(ras_west_breed, LIG_breeding_ovest_prediction==1), aes(x = x, y = y, fill = LIG_breeding_ovest_prediction), fill = "#ff6d00") + 
  geom_tile(data = subset(ras_west_wint, LIG_wintering_west_prediction==1), aes(x = x, y = y, fill = LIG_wintering_west_prediction), fill = "#ffc499") +
  coord_sf(xlim = c(150, -20), ylim = c(-38, 65), expand = F) +
  theme_void()
p_LIG_SDM

# -18 kya
p_minuseighteen_SDM <- ggplot() + 
  geom_sf(data = admin, fill="black", lwd=0, colour="black") +
  geom_tile(data = subset(ras_east_breed_past, p180_breeding_east_prediction_crop==1), aes(x = x, y = y, fill = MH_breeding_est_prediction), fill = "#5c7ec0") + 
  geom_tile(data = subset(ras_east_wint_past, p180_wintering_east_prediction_crop==1), aes(x = x, y = y, fill = MH_wintering_est_prediction), fill = "#BDCBE5") + 
  geom_tile(data = subset(ras_west_breed_past, p180_breeding_west_prediction_crop==1), aes(x = x, y = y, fill = MH_breeding_ovest_prediction), fill = "#ff6d00") + 
  geom_tile(data = subset(ras_west_wint_past, p180_wintering_west_prediction_crop==1), aes(x = x, y = y, fill = MH_wintering_west_prediction), fill = "#ffc499") +
  coord_sf(xlim = c(150, -20), ylim = c(-38, 65), expand = F) +
  theme_void()
p_minuseighteen_SDM

ggsave("-18_SDM.pdf", p_minuseighteen_SDM, device="pdf", units="cm", width=50, height=35, limitsize=FALSE)

# -14 kya
p_minusfourteen_SDM <- ggplot() + 
  geom_sf(data = admin, fill="black", lwd=0, colour="black") +
  geom_tile(data = subset(ras_east_breed_past, p140_breeding_east_prediction_crop==1), aes(x = x, y = y, fill = MH_breeding_est_prediction), fill = "#5c7ec0") + 
  geom_tile(data = subset(ras_east_wint_past, p140_wintering_east_prediction_crop==1), aes(x = x, y = y, fill = MH_wintering_est_prediction), fill = "#BDCBE5") + 
  geom_tile(data = subset(ras_west_breed_past, p140_breeding_west_prediction_crop==1), aes(x = x, y = y, fill = MH_breeding_ovest_prediction), fill = "#ff6d00") + 
  geom_tile(data = subset(ras_west_wint_past, p140_wintering_west_prediction_crop==1), aes(x = x, y = y, fill = MH_wintering_west_prediction), fill = "#ffc499") +
  coord_sf(xlim = c(150, -20), ylim = c(-38, 65), expand = F) +
  theme_void()
p_minusfourteen_SDM

ggsave("-14_SDM.pdf", p_minusfourteen_SDM, device="pdf", units="cm", width=50, height=35, limitsize=FALSE)

# -10 kya
p_minusten_SDM <- ggplot() + 
  geom_sf(data = admin, fill="black", lwd=0, colour="black") +
  geom_tile(data = subset(ras_east_breed_past, p100_breeding_east_prediction_crop==1), aes(x = x, y = y, fill = MH_breeding_est_prediction), fill = "#5c7ec0") + 
  geom_tile(data = subset(ras_east_wint_past, p100_wintering_east_prediction_crop==1), aes(x = x, y = y, fill = MH_wintering_est_prediction), fill = "#BDCBE5") + 
  geom_tile(data = subset(ras_west_breed_past, p100_breeding_west_prediction_crop==1), aes(x = x, y = y, fill = MH_breeding_ovest_prediction), fill = "#ff6d00") + 
  geom_tile(data = subset(ras_west_wint_past, p100_wintering_west_prediction_crop==1), aes(x = x, y = y, fill = MH_wintering_west_prediction), fill = "#ffc499") +
  coord_sf(xlim = c(150, -20), ylim = c(-38, 65), expand = F) +
  theme_void()
p_minusten_SDM

ggsave("-10_SDM.pdf", p_minusten_SDM, device="pdf", units="cm", width=50, height=35, limitsize=FALSE)

# -6 kya
p_minussix_SDM <- ggplot() + 
  geom_sf(data = admin, fill="black", lwd=0, colour="black") +
  geom_tile(data = subset(ras_east_breed_past, p60_breeding_east_prediction_crop==1), aes(x = x, y = y, fill = MH_breeding_est_prediction), fill = "#5c7ec0") + 
  geom_tile(data = subset(ras_east_wint_past, p60_wintering_east_prediction_crop==1), aes(x = x, y = y, fill = MH_wintering_est_prediction), fill = "#BDCBE5") + 
  geom_tile(data = subset(ras_west_breed_past, p60_breeding_west_prediction_crop==1), aes(x = x, y = y, fill = MH_breeding_ovest_prediction), fill = "#ff6d00") + 
  geom_tile(data = subset(ras_west_wint_past, p60_wintering_west_prediction_crop==1), aes(x = x, y = y, fill = MH_wintering_west_prediction), fill = "#ffc499") +
  coord_sf(xlim = c(150, -20), ylim = c(-38, 65), expand = F) +
  theme_void()
p_minussix_SDM

ggsave("-6_SDM.pdf", p_minussix_SDM, device="pdf", units="cm", width=50, height=35, limitsize=FALSE)

# -2 kya
p_minustwo_SDM <- ggplot() + 
  geom_sf(data = admin, fill="black", lwd=0, colour="black") +
  geom_tile(data = subset(ras_east_breed_past, p20_breeding_east_prediction_crop==1), aes(x = x, y = y, fill = MH_breeding_est_prediction), fill = "#5c7ec0") + 
  geom_tile(data = subset(ras_east_wint_past, p20_wintering_east_prediction_crop==1), aes(x = x, y = y, fill = MH_wintering_est_prediction), fill = "#BDCBE5") + 
  geom_tile(data = subset(ras_west_breed_past, p20_breeding_west_prediction_crop==1), aes(x = x, y = y, fill = MH_breeding_ovest_prediction), fill = "#ff6d00") + 
  geom_tile(data = subset(ras_west_wint_past, p20_wintering_west_prediction_crop==1), aes(x = x, y = y, fill = MH_wintering_west_prediction), fill = "#ffc499") +
  coord_sf(xlim = c(150, -20), ylim = c(-38, 65), expand = F) +
  theme_void()
p_minustwo_SDM

ggsave("-2_SDM.pdf", p_minustwo_SDM, device="pdf", units="cm", width=50, height=35, limitsize=FALSE)

# Figure 6: Future distribution forecasts
## Figure 6a: Future SDMs

### Calculate centroids of distributions
centroid_current_east <- colMeans(xyFromCell(ras_east_breed_past, which(ras_east_breed_past$p0_breeding_east_prediction_crop[]==1)))
centroid_2040_east <- colMeans(xyFromCell(ras_east_breed, which(ras_east_breed$ukesm1_585_breeding_east_prediction[]==1)))
centroid_2070_east <- colMeans(xyFromCell(ras_east_breed, which(ras_east_breed$ukesm1_585_breeding_east_prediction_2071_2100[]==1)))

centroid_current_west <- colMeans(xyFromCell(ras_west_breed_past, which(ras_west_breed_past$p0_breeding_west_prediction_crop[]==1)))
centroid_2040_west <- colMeans(xyFromCell(ras_west_breed, which(ras_west_breed$ukesm1_585_breeding_west_prediction[]==1)))
centroid_2070_west <- colMeans(xyFromCell(ras_west_breed, which(ras_west_breed$ukesm1_585_breeding_west_prediction_2071_2100[]==1)))

### Calculate distance between east and west centroids
dist_current <- as.numeric(earth.dist(rbind(t(as.data.frame(centroid_current_east)),t(as.data.frame(centroid_current_west))))) #4952 km
dist_2040 <- as.numeric(earth.dist(rbind(t(as.data.frame(centroid_2040_east)),t(as.data.frame(centroid_2040_west))))) #5508 km
dist_2070 <- as.numeric(earth.dist(rbind(t(as.data.frame(centroid_2070_east)),t(as.data.frame(centroid_2070_west))))) #5252 km

### Plot SDMs for current, 2040 and 2070 with distance between centroids
p_current_SDM <- ggplot() + 
  geom_sf(data = admin, fill="black", lwd=0, colour="black") +
  geom_tile(data = subset(ras_east_breed_past, p0_breeding_east_prediction_crop==1), aes(x = x, y = y, fill = MH_breeding_est_prediction), fill = "#5c7ec0") + 
  geom_tile(data = subset(ras_east_wint_past, p0_wintering_east_prediction_crop==1), aes(x = x, y = y, fill = MH_wintering_est_prediction), fill = "#BDCBE5") + 
  geom_tile(data = subset(ras_west_breed_past, p0_breeding_west_prediction_crop==1), aes(x = x, y = y, fill = MH_breeding_ovest_prediction), fill = "#ff6d00") + 
  geom_tile(data = subset(ras_west_wint_past, p0_wintering_west_prediction_crop==1), aes(x = x, y = y, fill = MH_wintering_west_prediction), fill = "#ffc499") +
  geom_segment(aes(x=centroid_current_east[1], y=centroid_current_east[2], xend=centroid_current_west[1], yend=centroid_current_west[2]), colour="#a4e036", size=1, alpha=0.5) +
  annotate("point", x=centroid_current_east[1], y=centroid_current_east[2], col="#a4e036", size=2) +
  annotate("point", x=centroid_current_west[1], y=centroid_current_west[2], col="#a4e036", size=2) +
  coord_sf(xlim = c(150, -20), ylim = c(-38, 65), expand = F) +
  theme_void()
p_current_SDM

p_2040_SDM <- ggplot() + 
  geom_sf(data = admin, fill="black", lwd=0, colour="black") +
  geom_tile(data = subset(ras_east_breed, ukesm1_585_breeding_east_prediction==1), aes(x = x, y = y, fill = ukesm1_585_breeding_east_prediction), fill = "#5c7ec0") + 
  geom_tile(data = subset(ras_east_wint, ukesm1_585_wintering_east_prediction==1), aes(x = x, y = y, fill = ukesm1_585_wintering_east_prediction), fill = "#BDCBE5") + 
  geom_tile(data = subset(ras_west_breed, ukesm1_585_breeding_west_prediction==1), aes(x = x, y = y, fill = ukesm1_585_breeding_west_prediction), fill = "#ff6d00") + 
  geom_tile(data = subset(ras_west_wint, ukesm1_585_wintering_west_prediction==1), aes(x = x, y = y, fill = ukesm1_585_wintering_west_prediction), fill = "#ffc499") +
  geom_segment(aes(x=centroid_2040_east[1], y=centroid_2040_east[2], xend=centroid_2040_west[1], yend=centroid_2040_west[2]), colour="#a4e036", size=1, alpha=0.5) +
  annotate("point", x=centroid_2040_east[1], y=centroid_2040_east[2], col="#a4e036", size=2) +
  annotate("point", x=centroid_2040_west[1], y=centroid_2040_west[2], col="#a4e036", size=2) +
  coord_sf(xlim = c(150, -20), ylim = c(-38, 65), expand = F) +
  theme_void()
p_2040_SDM

p_2070_SDM <- ggplot() + 
  geom_sf(data = admin, fill="black", lwd=0, colour="black") +
  geom_tile(data = subset(ras_east_breed, ukesm1_585_breeding_east_prediction_2071_2100==1), aes(x = x, y = y, fill = ukesm1_585_breeding_east_prediction_2071_2100), fill = "#5c7ec0") + 
  geom_tile(data = subset(ras_east_wint, ukesm1_585_wintering_east_prediction_2071_2100==1), aes(x = x, y = y, fill = ukesm1_585_wintering_east_prediction_2071_2100), fill = "#BDCBE5") + 
  geom_tile(data = subset(ras_west_breed, ukesm1_585_breeding_west_prediction_2071_2100==1), aes(x = x, y = y, fill = ukesm1_585_breeding_west_prediction_2071_2100), fill = "#ff6d00") + 
  geom_tile(data = subset(ras_west_wint, ukesm1_585_wintering_west_prediction_2071_2100==1), aes(x = x, y = y, fill = ukesm1_585_wintering_west_prediction_2071_2100), fill = "#ffc499") +
  geom_segment(aes(x=centroid_2070_east[1], y=centroid_2070_east[2], xend=centroid_2070_west[1], yend=centroid_2070_west[2]), colour="#a4e036", size=1, alpha=0.5) +
  annotate("point", x=centroid_2070_east[1], y=centroid_2070_east[2], col="#a4e036", size=2) +
  annotate("point", x=centroid_2070_west[1], y=centroid_2070_west[2], col="#a4e036", size=2) +
  coord_sf(xlim = c(150, -20), ylim = c(-38, 65), expand = F) +
  theme_void()
p_2070_SDM

dist_centroids <- plot_grid(p_current_SDM, p_2040_SDM, p_2070_SDM, ncol=3, nrow=1)

ggsave("future_SDMs_dist_centroids.pdf", dist_centroids, device="pdf", units="cm", width=15, height=7.5, limitsize=FALSE)

## Figure 6b: Future range trends
### See data source range_size_scenario_statistics_R.xlsx

# Figure 7: Genetic offsets
## See script https://github.com/jferrerobiol/lk_climate/blob/main/climate_associated_genetic_variation/GEA_genetic_offsets.R

# Figure S3: PCA downsampling
## 20000 SNPs
eigenvec_table <- read.table('LK_intersect_norelated_wCro_nosexchrom_down20000.eigenvec', header = FALSE)

### Add populations
pops <- as.character(c("ESN","ESN","ESN","GRC","GRC","ESN","GRC","GRC","ISR",
                       "TUR","ESN","KAZ","ISR","ISR","TUR","TUR","CRO","ISR","CRO","ESN",
                       "SIC","GRC","RUS","RUS","RUS","MOS","MOS","MON","ISR",
                       "ISR","ESS","ESS","ESS","ESN","ESN","ITN","ITN","ITN",
                       "SIC","RUS","MON","MON","ITN","ITS","ITS","ITS","ITS",
                       "SIC","RUS","ISR","ISR","ITS","ITS","SIC","KAZ","RUS",
                       "MON","CRO","ESS","ESS","ESS","ESS","ITN","ITN","ITN","ITS",
                       "ITS","GRL","GRL","GRL","GRG","GRG","KAZ","MOS","MOS",
                       "MON","MON","SIC","SIC","GRL","GRL","GRL","GRG","KAZ"))
eigenvec_table$Populations <- pops
### Add population order
pop_order <- c("ESN","ESS","SIC","ITS","ITN","CRO","GRG","GRC","GRL","TUR","ISR","KAZ","RUS","MON","MOS")
### Add colours
palette<-(c("#db0000","#dc2906","#dd4a11","#e16d19","#e68f21","#edb42a","#f5d933","#f5d933",
            "#e4e940","#cad34f","#b0bc66","#98a881","#7f929e","#677dbc","#5068dc"))
names(palette)<-pop_order

label <- eigenvec_table$V1
eigenvec_table <- eigenvec_table[-2]
eigenvec_table <- eigenvec_table[-1]

head(eigenvec_table)
for (i in 1:10){
  colnames(eigenvec_table)[i]<-paste0("PC",i)
}
eigenval <- read.table('LK_intersect_norelated_wCro_nosexchrom_down20000.eigenval', header = F)
percentage <- round(eigenval$V1/sum(eigenval$V1)*100,2)
percentage <- paste0(colnames(eigenvec_table)[1:10]," (",paste(as.character(percentage),"%)"))
percentage
eigenvec_table$Populations <- factor(eigenvec_table$Populations, levels=pop_order)

pca12_down20000<- ggplot(eigenvec_table,aes(x=PC1,y=PC2)) +
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

pca12_down20000

## 10000 SNPs
eigenvec_table <- read.table('LK_intersect_norelated_wCro_nosexchrom_down10000.eigenvec', header = FALSE)

### Add populations
eigenvec_table$Populations <- pops
### Add population order
pop_order <- c("ESN","ESS","SIC","ITS","ITN","CRO","GRG","GRC","GRL","TUR","ISR","KAZ","RUS","MON","MOS")
### Add colours
palette<-(c("#db0000","#dc2906","#dd4a11","#e16d19","#e68f21","#edb42a","#f5d933","#f5d933",
            "#e4e940","#cad34f","#b0bc66","#98a881","#7f929e","#677dbc","#5068dc"))
names(palette)<-pop_order

label <- eigenvec_table$V1
eigenvec_table <- eigenvec_table[-2]
eigenvec_table <- eigenvec_table[-1]

head(eigenvec_table)
for (i in 1:10){
  colnames(eigenvec_table)[i]<-paste0("PC",i)
}
eigenval <- read.table('LK_intersect_norelated_wCro_nosexchrom_down10000.eigenval', header = F)
percentage <- round(eigenval$V1/sum(eigenval$V1)*100,2)
percentage <- paste0(colnames(eigenvec_table)[1:10]," (",paste(as.character(percentage),"%)"))
percentage
eigenvec_table$Populations <- factor(eigenvec_table$Populations, levels=pop_order)

pca12_down10000<- ggplot(eigenvec_table,aes(x=PC1,y=PC2)) +
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

pca12_down10000

## 5000 SNPs
eigenvec_table <- read.table('LK_intersect_norelated_wCro_nosexchrom_down5000.eigenvec', header = FALSE)

### Add populations
eigenvec_table$Populations <- pops
### Add population order
pop_order <- c("ESN","ESS","SIC","ITS","ITN","CRO","GRG","GRC","GRL","TUR","ISR","KAZ","RUS","MON","MOS")
### Add colours
palette<-(c("#db0000","#dc2906","#dd4a11","#e16d19","#e68f21","#edb42a","#f5d933","#f5d933",
            "#e4e940","#cad34f","#b0bc66","#98a881","#7f929e","#677dbc","#5068dc"))
names(palette)<-pop_order

label <- eigenvec_table$V1
eigenvec_table <- eigenvec_table[-2]
eigenvec_table <- eigenvec_table[-1]

head(eigenvec_table)
for (i in 1:10){
  colnames(eigenvec_table)[i]<-paste0("PC",i)
}
eigenval <- read.table('LK_intersect_norelated_wCro_nosexchrom_down5000.eigenval', header = F)
percentage <- round(eigenval$V1/sum(eigenval$V1)*100,2)
percentage <- paste0(colnames(eigenvec_table)[1:10]," (",paste(as.character(percentage),"%)"))
percentage
eigenvec_table$Populations <- factor(eigenvec_table$Populations, levels=pop_order)

pca12_down5000<- ggplot(eigenvec_table,aes(x=PC1,y=PC2)) +
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

pca12_down5000

## 2000 SNPs
eigenvec_table <- read.table('LK_intersect_norelated_wCro_nosexchrom_down2000.eigenvec', header = FALSE)

### Add populations
eigenvec_table$Populations <- pops
### Add population order
pop_order <- c("ESN","ESS","SIC","ITS","ITN","CRO","GRG","GRC","GRL","TUR","ISR","KAZ","RUS","MON","MOS")
### Add colours
palette<-(c("#db0000","#dc2906","#dd4a11","#e16d19","#e68f21","#edb42a","#f5d933","#f5d933",
            "#e4e940","#cad34f","#b0bc66","#98a881","#7f929e","#677dbc","#5068dc"))
names(palette)<-pop_order

label <- eigenvec_table$V1
eigenvec_table <- eigenvec_table[-2]
eigenvec_table <- eigenvec_table[-1]

head(eigenvec_table)
for (i in 1:10){
  colnames(eigenvec_table)[i]<-paste0("PC",i)
}
eigenval <- read.table('LK_intersect_norelated_wCro_nosexchrom_down2000.eigenval', header = F)
percentage <- round(eigenval$V1/sum(eigenval$V1)*100,2)
percentage <- paste0(colnames(eigenvec_table)[1:10]," (",paste(as.character(percentage),"%)"))
percentage
eigenvec_table$Populations <- factor(eigenvec_table$Populations, levels=pop_order)

pca12_down2000<- ggplot(eigenvec_table,aes(x=PC1,y=PC2)) +
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

pca12_down2000

## 500 SNPs
eigenvec_table <- read.table('LK_intersect_norelated_wCro_nosexchrom_down500.eigenvec', header = FALSE)

### Add populations
eigenvec_table$Populations <- pops
### Add population order
pop_order <- c("ESN","ESS","SIC","ITS","ITN","CRO","GRG","GRC","GRL","TUR","ISR","KAZ","RUS","MON","MOS")
### Add colours
palette<-(c("#db0000","#dc2906","#dd4a11","#e16d19","#e68f21","#edb42a","#f5d933","#f5d933",
            "#e4e940","#cad34f","#b0bc66","#98a881","#7f929e","#677dbc","#5068dc"))
names(palette)<-pop_order

label <- eigenvec_table$V1
eigenvec_table <- eigenvec_table[-2]
eigenvec_table <- eigenvec_table[-1]

head(eigenvec_table)
for (i in 1:10){
  colnames(eigenvec_table)[i]<-paste0("PC",i)
}
eigenval <- read.table('LK_intersect_norelated_wCro_nosexchrom_down500.eigenval', header = F)
percentage <- round(eigenval$V1/sum(eigenval$V1)*100,2)
percentage <- paste0(colnames(eigenvec_table)[1:10]," (",paste(as.character(percentage),"%)"))
percentage
eigenvec_table$Populations <- factor(eigenvec_table$Populations, levels=pop_order)

pca12_down500<- ggplot(eigenvec_table,aes(x=PC1,y=PC2)) +
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

pca12_down500

pca_downsampled <- ggarrange(pca12, pca12_down20000, pca12_down10000, pca12_down5000,
                             pca12_down2000, pca12_down500, ncol=3, nrow=2,
                             labels= c("(a)","(b)","(c)","(d)","(e)","(f)"))

ggsave("LK_PC1-PC2_wallCro_downsampling.pdf", pca_downsampled, device=cairo_pdf, units="cm", width=60, height=30, limitsize=FALSE)

# Figure S4: Phist heatmap
library(RColorBrewer)
library(otuSummary)

fst <- read.delim("phist_pop_pairs_stacks_populations_75.tsv", header=T, row.names=1)

## Remove the comparisons with populations with 3 or less samples
fst <- fst[c(1:3,5,7:14), c(1:3,5,7:14)]

fst_long1 <- matrixConvert(fst, colname = c("pop1", "pop2", "Phist"))
fst_long2 <- matrixConvert(fst, colname = c("pop2", "pop1", "Phist"))
fst_long2 <- fst_long2[,c(2,1,3)]
fst_long <- rbind(fst_long1,fst_long2)
fst_long$pop1 <- factor(fst_long$pop1, levels = c("ESN","ESS","SIC","ITS","ITN","CRO","GRG","GRC","GRL","TUR","ISR","KAZ","RUS","MON","MOS"))
fst_long$pop2 <- factor(fst_long$pop2, levels = c("ESN","ESS","SIC","ITS","ITN","CRO","GRG","GRC","GRL","TUR","ISR","KAZ","RUS","MON","MOS"))
fst_long$Phist <- ifelse(fst_long$Phist<0,0,fst_long$Phist)

#Plot heatmap
heatmap_Phist <- ggplot(fst_long, aes(pop1, pop2)) +
  theme_bw() +
  geom_tile(aes(fill = Phist), color='white') +
  scale_fill_distiller(palette = "YlOrRd", direction = 1) +
  theme(axis.line=element_blank(),
        axis.ticks=element_blank(),
        panel.border=element_blank(),
        panel.grid=element_blank())

heatmap_Phist

ggsave("LK_phist_heatmap.pdf", heatmap_Phist, device="pdf", units="cm", width=15, height=12, limitsize=FALSE)

# Figure S5: Plot PC1 and PC3
pca13<- ggplot(eigenvec_table,aes(x=PC1,y=PC3)) +
  geom_point(aes(colour=Populations), size=2) + 
  scale_color_manual(values=palette) +
  theme_bw() +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  xlab(percentage[1]) +
  ylab(percentage[3])
ggsave("LK_PC1-PC3_wallCro.pdf", pca13, device="pdf", units="cm", width=20, height=15, limitsize=FALSE)

# Figure S6: Plot fineRADStructure
chunkfile<-"populations.haps_noinbred_def_def_chunks.out" ## RADpainter output file
mcmcfile<-"populations.haps_noinbred_def_def_chunks.mcmc.xml" ## finestructure mcmc file
treefile<-"populations.haps_noinbred_def_def_chunks.mcmcTree.xml" ## finestructure tree file
### 2) EDIT THIS PATH TO WHERE YOU WANT THE PLOTS:
plotsFolder <- "fineradstructure/"
### 3) SET VALUES FOR THESE VARIABLES: "analysisName" will be included in output plots
analysisName <- "stacks_unfiltered_noinbred_def_def";  maxIndv <- 10000; maxPop<-10000


### 4) EDIT THE PATH TO YOUR COPY of FinestructureLibrary.R
source("FinestructureLibrary.R", chdir = TRUE) # read in the R functions, which also calls the needed packages

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

# Figure S8: Plot individual heterozygosity
## Load individual heterozygosities file
ind_het_100 <- read.delim("intersect_Tasos_Joan.nomono.nosexchr_ncbi_annotated_nomissing.het")
popmap <- read.delim("popmap.txt")
ind_het_100 <- ind_het_100 %>% left_join(popmap, by=c("INDV"="ind"))
ind_het_100$pop <- factor(ind_het_100$pop, levels=c("ESN","ESS","SIC","ITS","ITN","CRO","GRG","GRC","GRL","TUR","ISR","KAZ","RUS","MON","MOS"))
mypal<-(c("#db0000","#dc2906","#dd4a11","#e16d19","#e68f21","#edb42a","#f5d933","#f5d933",
          "#e4e940","#cad34f","#b0bc66","#98a881","#7f929e","#677dbc","#5068dc"))

## Plot individual heterozygosities
p_ind_het_100 <- ggplot(ind_het_100, aes(x=pop, y=(N_SITES-O.HOM.)/N_SITES)) +
  theme_bw(base_family = "Arial") +
  theme(panel.grid = element_blank()) +
  theme(axis.text = element_text(color="black", size=16),
        axis.title = element_text(size=18),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_point(aes(col=pop), size=2) +
  stat_summary(fun="mean", geom="segment", mapping=aes(xend=..x.. - 0.25, yend=..y..), size=1) +
  stat_summary(fun="mean", geom="segment", mapping=aes(xend=..x.. + 0.25, yend=..y..), size=1) +
  scale_color_manual(values=mypal, guide="none") +
  annotate("rect", xmin=0, xmax=2.5, ymin=-Inf, ymax=Inf, alpha=0.3, fill="#DC2906") +
  annotate("rect", xmin=2.5, xmax=10.5, ymin=-Inf, ymax=Inf, alpha=0.3, fill="#EADE32") +
  annotate("rect", xmin=10.5, xmax=11.5, ymin=-Inf, ymax=Inf, alpha=0.3, fill="#B0BC66") +
  annotate("rect", xmin=11.5, xmax=Inf, ymin=-Inf, ymax=Inf, alpha=0.3, fill="#5C7EC0") +
  labs(x="", y="Proportion of heterozygote sites")

p_ind_het_100
ggsave("LK_ind_hets_100.pdf", p_ind_het_100, device=cairo_pdf, units="cm", width=16.479, height=7.2, limitsize=FALSE)

# Figure S9: genome-wide heterozygosity
## Load libraries
library(ggridges)
library(ggpubr)

# Read individual heterozigosities and plot
#paternal
pat_het <- read.delim("paternal.merged.dups.bam.1MB_windows_meet_filters.nhetSNPs.bed.gz", header=F)
names(pat_het)[c(seq(1,5))] <- c("chr", "window_start","window_end","sites_called","nsnps")
pat_het$window_size <- pat_het$window_end-pat_het$window_start
pat_het$chr <- factor(pat_het$chr, levels=c("SUPER_1","SUPER_2","SUPER_3","SUPER_4","SUPER_5","SUPER_6",
                                            "SUPER_7","SUPER_8","SUPER_9","SUPER_10","SUPER_11","SUPER_12",
                                            "SUPER_13","SUPER_14","SUPER_15","SUPER_16","SUPER_17","SUPER_18",
                                            "SUPER_19","SUPER_20","SUPER_21","SUPER_22"))
pat_het <- pat_het %>% arrange(chr, window_start)
pat_het$xaxis <- seq(1,length(pat_het$chr))
pat_het <- pat_het[pat_het$window_size>500000,]
pat_het <- pat_het[pat_het$sites_called > pat_het$window_size*0.8,]
pat_het$pi <- pat_het$nsnps/pat_het$sites_called

loc_chr <- pat_het %>% group_by(chr) %>%
  summarise(loc_chr=mean(xaxis))

dens_pat <- ggplot(pat_het) +
  theme_linedraw(base_family = "Arial") +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size=1),
        legend.position='none',
        axis.text = element_text(color="black", size=14),
        axis.title = element_text(size=16)) +
  geom_histogram(binwidth=0.0003, aes(x=pi), col="#ff6d00", fill="#ff6d00", alpha=0.3) +
  xlim(0,0.01) + ylim(0,130) +
  labs(x=expression(paste("Observed ",pi)),y="# of windows") + 
  theme(legend.position = "none")

col <- rep(c("#ff6d00","#ffc499"),11)
distr_pat <- ggplot(pat_het, aes(x=xaxis, ymin=0, ymax=pi, col=chr)) +
  theme_linedraw(base_family = "Arial") +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size=1),
        legend.position='none',
        axis.text = element_text(color="black", size=14),
        axis.title = element_text(size=16)) +
  scale_x_continuous(breaks=loc_chr$loc_chr, labels=rownames(loc_chr)) +
  geom_linerange() +
  scale_color_manual(values=col) +
  ylim(0,0.01) +
  labs(x="Chromosome",y=expression(paste("Observed ",pi)))

mean(pat_het$pi)

#maternal
mat_het <- read.delim("maternal.merged.dups.bam.1MB_windows_meet_filters.nhetSNPs.bed.gz", header=F)
names(mat_het)[c(seq(1,5))] <- c("chr", "window_start","window_end","sites_called","nsnps")
mat_het$window_size <- mat_het$window_end-mat_het$window_start
mat_het$chr <- factor(mat_het$chr, levels=c("SUPER_1","SUPER_2","SUPER_3","SUPER_4","SUPER_5","SUPER_6",
                                            "SUPER_7","SUPER_8","SUPER_9","SUPER_10","SUPER_11","SUPER_12",
                                            "SUPER_13","SUPER_14","SUPER_15","SUPER_16","SUPER_17","SUPER_18",
                                            "SUPER_19","SUPER_20","SUPER_21","SUPER_22"))
mat_het <- mat_het %>% arrange(chr, window_start)
mat_het$xaxis <- seq(1,length(mat_het$chr))
mat_het <- mat_het[mat_het$window_size>500000,]
mat_het <- mat_het[mat_het$sites_called > mat_het$window_size*0.8,]
mat_het$pi <- mat_het$nsnps/mat_het$sites_called

dens_mat <- ggplot(mat_het) +
  theme_linedraw(base_family = "Arial") +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size=1),
        legend.position='none',
        axis.text = element_text(color="black", size=14),
        axis.title = element_text(size=16)) +
  geom_histogram(binwidth=0.0003, aes(x=pi), col="#ff6d00", fill="#ff6d00", alpha=0.3) +
  xlim(0,0.01) + ylim(0,130) +
  labs(x=expression(paste("Observed ",pi)),y="# of windows") + 
  theme(legend.position = "none")

col <- rep(c("#ff6d00","#ffc499"),11)
distr_mat <- ggplot(mat_het, aes(x=xaxis, ymin=0, ymax=pi, col=chr)) +
  theme_linedraw(base_family = "Arial") +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size=1),
        legend.position='none',
        axis.text = element_text(color="black", size=14),
        axis.title = element_text(size=16)) +
  scale_x_continuous(breaks=loc_chr$loc_chr, labels=rownames(loc_chr)) +
  geom_linerange() +
  scale_color_manual(values=col) +
  ylim(0,0.01) +
  labs(x="Chromosome",y=expression(paste("Observed ",pi)))

mean(mat_het$pi)

#LK_83M
LK_83M_het <- read.delim("LK_83M.merged.dups.bam.1MB_windows_meet_filters.nhetSNPs.bed.gz", header=F)
names(LK_83M_het)[c(seq(1,5))] <- c("chr", "window_start","window_end","sites_called","nsnps")
LK_83M_het$window_size <- LK_83M_het$window_end-LK_83M_het$window_start
LK_83M_het$chr <- factor(LK_83M_het$chr, levels=c("SUPER_1","SUPER_2","SUPER_3","SUPER_4","SUPER_5","SUPER_6",
                                                  "SUPER_7","SUPER_8","SUPER_9","SUPER_10","SUPER_11","SUPER_12",
                                                  "SUPER_13","SUPER_14","SUPER_15","SUPER_16","SUPER_17","SUPER_18",
                                                  "SUPER_19","SUPER_20","SUPER_21","SUPER_22"))
LK_83M_het <- LK_83M_het %>% arrange(chr, window_start)
LK_83M_het$xaxis <- seq(1,length(LK_83M_het$chr))
LK_83M_het <- LK_83M_het[LK_83M_het$window_size>500000,]
LK_83M_het <- LK_83M_het[LK_83M_het$sites_called > LK_83M_het$window_size*0.8,]
LK_83M_het$pi <- LK_83M_het$nsnps/LK_83M_het$sites_called

dens_LK_83M <- ggplot(LK_83M_het) +
  theme_linedraw(base_family = "Arial") +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size=1),
        legend.position='none',
        axis.text = element_text(color="black", size=14),
        axis.title = element_text(size=16)) +
  geom_histogram(binwidth=0.0003, aes(x=pi), col="#ff6d00", fill="#ff6d00", alpha=0.3) +
  xlim(0,0.01) + ylim(0,130) +
  labs(x=expression(paste("Observed ",pi)),y="# of windows") + 
  theme(legend.position = "none")

col <- rep(c("#ff6d00","#ffc499"),11)
distr_LK_83M <- ggplot(LK_83M_het, aes(x=xaxis, ymin=0, ymax=pi, col=chr)) +
  theme_linedraw(base_family = "Arial") +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size=1),
        legend.position='none',
        axis.text = element_text(color="black", size=14),
        axis.title = element_text(size=16)) +
  scale_x_continuous(breaks=loc_chr$loc_chr, labels=rownames(loc_chr)) +
  geom_linerange() +
  scale_color_manual(values=col) +
  ylim(0,0.01) +
  labs(x="Chromosome",y=expression(paste("Observed ",pi)))

mean(LK_83M_het$pi)

het_Eur <- ggarrange(dens_pat, dens_mat, dens_LK_83M, distr_pat, distr_mat, distr_LK_83M,
                     ncol=3, nrow=2, align = "v")
ggsave("heterozygosities_Eur.pdf", het_Eur, device=cairo_pdf, units="cm", width=40, height=10, limitsize=FALSE)

#LK_F7
LK_F7_het <- read.delim("LK_F7.merged.dups.bam.1MB_windows_meet_filters.nhetSNPs.bed.gz", header=F)
names(LK_F7_het)[c(seq(1,5))] <- c("chr", "window_start","window_end","sites_called","nsnps")
LK_F7_het$window_size <- LK_F7_het$window_end-LK_F7_het$window_start
LK_F7_het$chr <- factor(LK_F7_het$chr, levels=c("SUPER_1","SUPER_2","SUPER_3","SUPER_4","SUPER_5","SUPER_6",
                                                "SUPER_7","SUPER_8","SUPER_9","SUPER_10","SUPER_11","SUPER_12",
                                                "SUPER_13","SUPER_14","SUPER_15","SUPER_16","SUPER_17","SUPER_18",
                                                "SUPER_19","SUPER_20","SUPER_21","SUPER_22"))
LK_F7_het <- LK_F7_het %>% arrange(chr, window_start)
LK_F7_het$xaxis <- seq(1,length(LK_F7_het$chr))
LK_F7_het <- LK_F7_het[LK_F7_het$window_size>500000,]
LK_F7_het <- LK_F7_het[LK_F7_het$sites_called > LK_F7_het$window_size*0.8,]
LK_F7_het$pi <- LK_F7_het$nsnps/LK_F7_het$sites_called

dens_LK_F7 <- ggplot(LK_F7_het) +
  theme_linedraw(base_family = "Arial") +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size=1),
        legend.position='none',
        axis.text = element_text(color="black", size=14),
        axis.title = element_text(size=16)) +
  geom_histogram(binwidth=0.0003, aes(x=pi), col="#5c7ec0", fill="#5c7ec0", alpha=0.3) +
  xlim(0,0.01) + ylim(0,130) +
  labs(x=expression(paste("Observed ",pi)),y="# of windows") + 
  theme(legend.position = "none")

col <- rep(c("#5c7ec0","#BDCBE5"),11)
distr_LK_F7 <- ggplot(LK_F7_het, aes(x=xaxis, ymin=0, ymax=pi, col=chr)) +
  theme_linedraw(base_family = "Arial") +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size=1),
        legend.position='none',
        axis.text = element_text(color="black", size=14),
        axis.title = element_text(size=16)) +
  scale_x_continuous(breaks=loc_chr$loc_chr, labels=rownames(loc_chr)) +
  geom_linerange() +
  scale_color_manual(values=col) +
  ylim(0,0.01) +
  labs(x="Chromosome",y=expression(paste("Observed ",pi)))

mean(LK_F7_het$pi)

#LK_F8
LK_F8_het <- read.delim("LK_F8.merged.dups.bam.1MB_windows_meet_filters.nhetSNPs.bed.gz", header=F)
names(LK_F8_het)[c(seq(1,5))] <- c("chr", "window_start","window_end","sites_called","nsnps")
LK_F8_het$window_size <- LK_F8_het$window_end-LK_F8_het$window_start
LK_F8_het$chr <- factor(LK_F8_het$chr, levels=c("SUPER_1","SUPER_2","SUPER_3","SUPER_4","SUPER_5","SUPER_6",
                                                "SUPER_7","SUPER_8","SUPER_9","SUPER_10","SUPER_11","SUPER_12",
                                                "SUPER_13","SUPER_14","SUPER_15","SUPER_16","SUPER_17","SUPER_18",
                                                "SUPER_19","SUPER_20","SUPER_21","SUPER_22"))
LK_F8_het <- LK_F8_het %>% arrange(chr, window_start)
LK_F8_het$xaxis <- seq(1,length(LK_F8_het$chr))
LK_F8_het <- LK_F8_het[LK_F8_het$window_size>500000,]
LK_F8_het <- LK_F8_het[LK_F8_het$sites_called > LK_F8_het$window_size*0.8,]
LK_F8_het$pi <- LK_F8_het$nsnps/LK_F8_het$sites_called

dens_LK_F8 <- ggplot(LK_F8_het) +
  theme_linedraw(base_family = "Arial") +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size=1),
        legend.position='none',
        axis.text = element_text(color="black", size=14),
        axis.title = element_text(size=16)) +
  geom_histogram(binwidth=0.0003, aes(x=pi), col="#5c7ec0", fill="#5c7ec0", alpha=0.3) +
  xlim(0,0.01) + ylim(0,130) +
  labs(x=expression(paste("Observed ",pi)),y="# of windows") + 
  theme(legend.position = "none")

col <- rep(c("#5c7ec0","#BDCBE5"),11)
distr_LK_F8 <- ggplot(LK_F8_het, aes(x=xaxis, ymin=0, ymax=pi, col=chr)) +
  theme_linedraw(base_family = "Arial") +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size=1),
        legend.position='none',
        axis.text = element_text(color="black", size=14),
        axis.title = element_text(size=16)) +
  scale_x_continuous(breaks=loc_chr$loc_chr, labels=rownames(loc_chr)) +
  geom_linerange() +
  scale_color_manual(values=col) +
  ylim(0,0.01) +
  labs(x="Chromosome",y=expression(paste("Observed ",pi)))

mean(LK_F8_het$pi)

#LK_M2
LK_M2_het <- read.delim("LK_M2.merged.dups.bam.1MB_windows_meet_filters.nhetSNPs.bed.gz", header=F)
names(LK_M2_het)[c(seq(1,5))] <- c("chr", "window_start","window_end","sites_called","nsnps")
LK_M2_het$window_size <- LK_M2_het$window_end-LK_M2_het$window_start
LK_M2_het$chr <- factor(LK_M2_het$chr, levels=c("SUPER_1","SUPER_2","SUPER_3","SUPER_4","SUPER_5","SUPER_6",
                                                "SUPER_7","SUPER_8","SUPER_9","SUPER_10","SUPER_11","SUPER_12",
                                                "SUPER_13","SUPER_14","SUPER_15","SUPER_16","SUPER_17","SUPER_18",
                                                "SUPER_19","SUPER_20","SUPER_21","SUPER_22"))
LK_M2_het <- LK_M2_het %>% arrange(chr, window_start)
LK_M2_het$xaxis <- seq(1,length(LK_M2_het$chr))
LK_M2_het <- LK_M2_het[LK_M2_het$window_size>500000,]
LK_M2_het <- LK_M2_het[LK_M2_het$sites_called > LK_M2_het$window_size*0.8,]
LK_M2_het$pi <- LK_M2_het$nsnps/LK_M2_het$sites_called

dens_LK_M2 <- ggplot(LK_M2_het) +
  theme_linedraw(base_family = "Arial") +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size=1),
        legend.position='none',
        axis.text = element_text(color="black", size=14),
        axis.title = element_text(size=16)) +
  geom_histogram(binwidth=0.0003, aes(x=pi), col="#5c7ec0", fill="#5c7ec0", alpha=0.3) +
  xlim(0,0.01) + ylim(0,130) +
  labs(x=expression(paste("Observed ",pi)),y="# of windows") + 
  theme(legend.position = "none")

col <- rep(c("#5c7ec0","#BDCBE5"),11)
distr_LK_M2 <- ggplot(LK_M2_het, aes(x=xaxis, ymin=0, ymax=pi, col=chr)) +
  theme_linedraw(base_family = "Arial") +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size=1),
        legend.position='none',
        axis.text = element_text(color="black", size=14),
        axis.title = element_text(size=16)) +
  scale_x_continuous(breaks=loc_chr$loc_chr, labels=rownames(loc_chr)) +
  geom_linerange() +
  scale_color_manual(values=col) +
  ylim(0,0.01) +
  labs(x="Chromosome",y=expression(paste("Observed ",pi)))

mean(LK_M2_het$pi)

het_Asi <- ggarrange(dens_LK_F7, dens_LK_F8, dens_LK_M2, distr_LK_F7, distr_LK_F8, distr_LK_M2,
                     ncol=3, nrow=2, align = "v")
ggsave("heterozygosities_Asi.pdf", het_Asi, device=cairo_pdf, units="cm", width=40, height=10, limitsize=FALSE)

## Mean datasets
## Eur
Eur_het <- pat_het %>%
  full_join(mat_het, by=c("chr", "window_start", "window_end", "xaxis")) %>%
  full_join(LK_83M_het, by=c("chr", "window_start", "window_end", "xaxis"))

Eur_het$mean_pi <- (Eur_het$pi.x + Eur_het$pi.y + Eur_het$pi)/3

dens_Eur <- ggplot(Eur_het) +
  theme_linedraw(base_family = "Arial") +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size=1),
        legend.position='none',
        axis.text = element_text(color="black", size=14),
        axis.title = element_text(size=16)) +
  geom_histogram(binwidth=0.0003, aes(x=mean_pi), col="#ff6d00", fill="#ff6d00", alpha=0.3) +
  xlim(0,0.01) + ylim(0,140) +
  labs(x=expression(paste("Observed ",pi)),y="# of windows") + 
  theme(legend.position = "none")

col <- rep(c("#ff6d00","#ffc499"),11)
distr_Eur <- ggplot(Eur_het, aes(x=xaxis, ymin=0, ymax=mean_pi, col=chr)) +
  theme_linedraw(base_family = "Arial") +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size=1),
        legend.position='none',
        axis.text = element_text(color="black", size=14),
        axis.title = element_text(size=16)) +
  scale_x_continuous(breaks=loc_chr$loc_chr, labels=rownames(loc_chr)) +
  geom_linerange() +
  scale_color_manual(values=col) +
  ylim(0,0.01) +
  labs(x="Chromosome",y=expression(paste("Observed ",pi)))

mean(na.omit(Eur_het$mean_pi))

## Asi
Asi_het <- LK_F8_het %>%
  full_join(LK_F7_het, by=c("chr", "window_start", "window_end", "xaxis")) %>%
  full_join(LK_M2_het, by=c("chr", "window_start", "window_end", "xaxis"))

Asi_het$mean_pi <- (Asi_het$pi.x + Asi_het$pi.y + Asi_het$pi.x)/3

dens_Asi <- ggplot(Asi_het) +
  theme_linedraw(base_family = "Arial") +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size=1),
        legend.position='none',
        axis.text = element_text(color="black", size=14),
        axis.title = element_text(size=16)) +
  geom_histogram(binwidth=0.0003, aes(x=mean_pi), col="#5c7ec0", fill="#5c7ec0", alpha=0.3) +
  xlim(0,0.01) + ylim(0,140) +
  labs(x=expression(paste("Observed ",pi)),y="# of windows") + 
  theme(legend.position = "none")

col <- rep(c("#5c7ec0","#BDCBE5"),11)
distr_Asi <- ggplot(Asi_het, aes(x=xaxis, ymin=0, ymax=mean_pi, col=chr)) +
  theme_linedraw(base_family = "Arial") +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size=1),
        legend.position='none',
        axis.text = element_text(color="black", size=14),
        axis.title = element_text(size=16)) +
  scale_x_continuous(breaks=loc_chr$loc_chr, labels=rownames(loc_chr)) +
  geom_linerange() +
  scale_color_manual(values=col) +
  ylim(0,0.01) +
  labs(x="Chromosome",y=expression(paste("Observed ",pi)))

mean(na.omit(Asi_het$mean_pi))

het_Eur_Asi <- ggarrange(distr_Eur, dens_Eur, distr_Asi, dens_Asi,
                         ncol=2, nrow=2, align = "v")

ggsave("heterozygosities_Eur_Asi.pdf", het_Eur_Asi, device=cairo_pdf, units="cm", width=26, height=10, limitsize=FALSE)

# Figure S10: NeighbourNet
## Load libraries
library(vcfR)
library(ape)
library(adegenet)
library(phangorn)

## Load vcf file and convert to nexus distance file
vcf <- read.vcfR("populations.snps.wFTin.vcf.gz")
dnabin <- vcfR2DNAbin(vcf,extract.indels=T,consensus=T,extract.haps=F,unphased_as_NA=F)

dist <- dist.dna(dnabin, model = "TN93", variance = FALSE,
                 gamma = FALSE, pairwise.deletion = TRUE,
                 base.freq = NULL, as.matrix = TRUE)

write.nexus.dist(dist, file = "LK_stacks_woutgroup_dist.nex", append = FALSE, upper = FALSE,
                 diag = TRUE, digits = getOption("digits"))
## The resulting file is used as input in SplitsTree5 to generate the NeighbourNet

# Figure S11: Admixture cross-validation

cverror <- read.delim("../github/lk_climate/source_data/data/cross_validation.txt", header=F)

plot(cverror, type="b", pch = 16, cex = 1, lwd = 1)

# Figure S13: Elevational distribution

## Load libraries
library(elevatr)
library(gghalves)
library(ggdist)

## Load data for each of the groups
east_breed <- read.csv("eastern_data_bioclim_5km_buffer_no_dupl.csv")
east_breed <- east_breed %>% dplyr::select(X1,Y1)
east_breed$group <- "eastern_breeding"
west_breed <- read.csv("western_data_bioclim_5km_buffer_no_dupl.csv")
west_breed <- west_breed %>% dplyr::select(X1,Y1)
west_breed$group <- "western_breeding"
east_nonbreed <- read.csv("wintering_eastern_bioclim_5km_buffer_no_dupl.csv")
east_nonbreed <- east_nonbreed %>% dplyr::select(X1,Y1)
east_nonbreed$group <- "eastern_nonbreeding"
west_nonbreed <- read.csv("wintering_western_bioclim_5km_buffer_no_dupl.csv")
west_nonbreed <- west_nonbreed %>% dplyr::select(X1,Y1)
west_nonbreed$group <- "western_nonbreeding"

presence <- rbind(east_breed,west_breed,east_nonbreed,west_nonbreed)
colnames(presence) <- c("x","y","group")
presence$group <- factor(presence$group, levels=c("western_breeding","eastern_breeding","western_nonbreeding","eastern_nonbreeding"),
                         labels=c("Western breeding","Eastern breeding","Western non-breeding","Eastern non-breeding"))

### Retrieve altitude data
df_elev_epqs <- get_elev_point(presence[,c(1,2)], prj = "EPSG:4326", src = "aws")
presence$elevation <- df_elev_epqs@data$elevation

p_elevation <- ggplot(presence, aes(x=group, y=elevation, colour=group, fill=group)) + 
  theme_bw(base_family = "Arial") +
  theme(panel.grid = element_blank(),
        legend.position='none',
        axis.text = element_text(color="black", size=12),
        axis.title = element_text(size=14)) +
  ggdist::stat_halfeye(adjust = .5, width = .3, .width = 0, justification = -.3, point_colour = NA, alpha = .5) + 
  geom_boxplot(width = .1, outlier.shape = NA, alpha = .5) +
  gghalves::geom_half_point(side = "l", range_scale = .4, alpha = .5) +
  scale_colour_manual(values=c("#ff6d00","#5c7ec0","#ffc499","#bdcbe5")) +
  scale_fill_manual(values=c("#ff6d00","#5c7ec0","#ffc499","#bdcbe5")) +
  xlab("") + ylab("Elevation (m)")
p_elevation
ggsave("LK_elevation_per_group.pdf", p_elevation, device=cairo_pdf, units="cm", width=24, height=12, limitsize=FALSE)

# Figure S14: Occurrence data and predicted distribution ranges
occurrence <- read.csv("occurrence_data.csv")
occurrence$ESU_season <- paste(occurrence$ESU, occurrence$season, sep="_")

p_occurrence_sdm <- ggplot() + 
  theme_bw(base_family = "Arial") +
  theme(panel.grid = element_blank(),
        legend.position='none',
        axis.text = element_text(color="black", size=12),
        axis.title = element_text(size=14)) +
  geom_sf(data = admin, fill="#9e9a95", lwd=0, colour="#9e9a95") +
  geom_tile(data = subset(ras_east_breed, breeding_prediction==1), aes(x = x, y = y, fill = breeding_prediction), fill = "#5c7ec0") + 
  geom_tile(data = subset(ras_east_wint, wintering_prediction==1), aes(x = x, y = y, fill = wintering_prediction), fill = "#BDCBE5") + 
  geom_tile(data = subset(ras_west_breed, breeding_prediction==1), aes(x = x, y = y, fill = breeding_prediction), fill = "#ff6d00") + 
  geom_tile(data = subset(ras_west_wint, wintering_prediction==1), aes(x = x, y = y, fill = wintering_prediction), fill = "#ffc499") +
  geom_point(shape=21, data = occurrence, aes(x = longitude, y = latitude, fill = ESU_season), colour="black") +
  scale_fill_manual(values=c("#5c7ec0","#BDCBE5","#ff6d00","#ffc499"), guide="none") +
  coord_sf(xlim = c(150, -20), ylim = c(-38, 65), expand = F) +
  xlab("Longitude") + ylab("Latitude")
p_occurrence_sdm

ggsave("occurrence_current_SDM.pdf", p_occurrence_sdm, device=cairo_pdf, units="cm", width=20, height=14, limitsize=FALSE)

# Figure S15: Manhattan plot climate-associated SNPs
## See script https://github.com/jferrerobiol/lk_climate/blob/main/climate_associated_genetic_variation/GEA_genetic_offsets.R

# Figure S16: Gradient forest variable importance
## See script https://github.com/jferrerobiol/lk_climate/blob/main/climate_associated_genetic_variation/GEA_genetic_offsets.R

# Figure S18: Latitudinal trends
### See data source range_size_scenario_statistics_R.xlsx

# Figure S19: Uncropped SDMs

## Load SDM rasters
uncropped.list <- list.files(path="SDMs_uncropped", pattern =".tif", full.names=TRUE)
ras_uncropped <- stack(uncropped.list)
ras_uncropped <- as(ras_uncropped, "SpatialPixelsDataFrame")
ras_uncropped <- as.data.frame(ras_uncropped)

ras_uncropped_current_east <- stack("SDMs_uncropped/wintering_east_chelsa_pres_recl_uncropped_land.tif")
ras_uncropped_current_east <- as(ras_uncropped_current_east, "SpatialPixelsDataFrame")
ras_uncropped_current_east <- as.data.frame(ras_uncropped_current_east)

ras_uncropped_current_west <- stack("SDMs_uncropped/wintering_west_chelsa_pres_recl_uncropped_land.tif")
ras_uncropped_current_west <- as(ras_uncropped_current_west, "SpatialPixelsDataFrame")
ras_uncropped_current_west <- as.data.frame(ras_uncropped_current_west)

p_western_present_uncropped <- ggplot() + 
  geom_sf(data = admin, fill="black", lwd=0, colour="black") +
  geom_tile(data = subset(ras_uncropped_current_west, layer==1), aes(x = x, y = y), fill = "#ffc499") + 
  coord_sf(xlim = c(150, -20), ylim = c(-38, 65), expand = F) +
  theme_void()
p_western_present_uncropped

p_eastern_present_uncropped <- ggplot() + 
  geom_sf(data = admin, fill="black", lwd=0, colour="black") +
  geom_tile(data = subset(ras_uncropped_current_east, layer==1), aes(x = x, y = y), fill = "#BDCBE5") + 
  coord_sf(xlim = c(150, -20), ylim = c(-38, 65), expand = F) +
  theme_void()
p_western_present_uncropped

p_western_2040_uncropped <- ggplot() + 
  geom_sf(data = admin, fill="black", lwd=0, colour="black") +
  geom_tile(data = subset(ras_uncropped, ukesm1_585_wintering_west_prediction==1), aes(x = x, y = y), fill = "#ffc499") + 
  coord_sf(xlim = c(150, -20), ylim = c(-38, 65), expand = F) +
  theme_void()
p_western_2040_uncropped

p_eastern_2040_uncropped <- ggplot() + 
  geom_sf(data = admin, fill="black", lwd=0, colour="black") +
  geom_tile(data = subset(ras_uncropped, ukesm1_585_wintering_east_prediction==1), aes(x = x, y = y), fill = "#BDCBE5") + 
  coord_sf(xlim = c(150, -20), ylim = c(-38, 65), expand = F) +
  theme_void()
p_eastern_2040_uncropped

p_western_2070_uncropped <- ggplot() + 
  geom_sf(data = admin, fill="black", lwd=0, colour="black") +
  geom_tile(data = subset(ras_uncropped, ukesm1_585_wintering_west_prediction_2071_2100==1), aes(x = x, y = y), fill = "#ffc499") + 
  coord_sf(xlim = c(150, -20), ylim = c(-38, 65), expand = F) +
  theme_void()
p_western_2070_uncropped

p_eastern_2070_uncropped <- ggplot() + 
  geom_sf(data = admin, fill="black", lwd=0, colour="black") +
  geom_tile(data = subset(ras_uncropped, ukesm1_585_wintering_east_prediction_2071_2100==1), aes(x = x, y = y), fill = "#BDCBE5") + 
  coord_sf(xlim = c(150, -20), ylim = c(-38, 65), expand = F) +
  theme_void()
p_eastern_2070_uncropped

uncropped_nonbreeding_present_future <- plot_grid(p_western_present_uncropped, p_western_2040_uncropped, p_western_2070_uncropped,
                                                  p_eastern_present_uncropped, p_eastern_2040_uncropped, p_eastern_2070_uncropped,
                                                  ncol=3, nrow=2, labels = c("(a)","(b)","(c)","(d)","(e)","(f)"),
                                                  label_size = 10, label_y = 1)

# Figure S20: Genetic offsets under moderate warming climate
## See script https://github.com/jferrerobiol/lk_climate/blob/main/climate_associated_genetic_variation/GEA_genetic_offsets.R
