library(geosphere)
library(sf)
library(marmap)
library(fossil)
library(reshape2)
library(dendextend)
library(adegenet)
library(metagMisc)
library(tidyverse)
library(ggpubr)
library(usedist)
library(vegan)
library(factoextra)
library(ecodist)

###################################
## Script to compute IBD and IBE ##
###################################

setwd("~/Dropbox/Postdoc_Milan/LK_Joan/IBD_analyses")

# IBD
# Load colony coordinates
bred_sites <- read.delim("loc_coordinates_bo.txt", dec = ",")
bred_sites_coord <- bred_sites[,c(3,2)]
rownames(bred_sites_coord) <- bred_sites$pop
colnames(bred_sites_coord) <- c("x","y")

dgeo <- round(earth.dist(bred_sites_coord), 0)
dgeo <- dist_setNames(dgeo, rownames(bred_sites_coord))

# Perform Mantel tests between genetic and geographic distances
## Load fst table
fst <- read.delim("phist_pop_pairs_stacks_populations_75.tsv", header=T, row.names=1)

## Remove the comparisons with populations with 3 or less samples
fst <- fst[c(1:3,5,7:14), c(1:3,5,7:14)]

## Convert to distance matrix
dfst <- as.dist(fst)

## Order fst and distance matrices the same way
dfst <- sort_dist_mat(dfst, by_rows = TRUE, by_cols = TRUE)
dgeo <- sort_dist_mat(dgeo, by_rows = TRUE, by_cols = TRUE)
dgen <- dfst/(1-dfst)

## Perform mantel test in vegan
ibd <- vegan::mantel(dgen, dgeo, method = "pearson", permutations = 9999, na.rm = TRUE)

## Perform mantel test for both ESUs separately
dgeo_western <- dist_subset(dgeo, c("ESN","ESS","GRC","GRL","ISR","ITN","ITS","SIC"))
dgen_western <- dist_subset(dgen, c("ESN","ESS","GRC","GRL","ISR","ITN","ITS","SIC"))

dgeo_eastern <- dist_subset(dgeo, c("KAZ","MON","MOS","RUS"))
dgen_eastern <- dist_subset(dgen, c("KAZ","MON","MOS","RUS"))

ibd_western <- vegan::mantel(dgen_western, dgeo_western, method = "pearson", permutations = 9999, na.rm = TRUE)
ibd_eastern <- vegan::mantel(dgen_eastern, dgeo_eastern, method = "pearson", permutations = 9999, na.rm = TRUE)

## Plot geographic distance against genetic distance
d_dgen <- dist2list(dgen, tri = TRUE)
d_dgeo <- dist2list(dgeo, tri = TRUE)
d <- merge(d_dgen, d_dgeo, by=c("col","row"))
colnames(d) <- c("p1","p2","dgen","dgeo")
EUR <- c("ESS","ESN","SIC","ITS","ITN","GRL","GRC","ISR")
ASI <- c("KAZ","MON","MOS","RUS")

d <- mutate(d, comparison = case_when(p1 %in% EUR & p2 %in% EUR ~ "within_gen_cluster",
                                      p1 %in% ASI & p2 %in% ASI ~ "within_gen_cluster",
                                      p1 %in% EUR & p2 %in% ASI ~ "between_gen_cluster",
                                      p1 %in% ASI & p2 %in% EUR ~ "between_gen_cluster"))
d$comparison <- as.factor(d$comparison)

p_IBD <- ggplot(d, aes(x=dgeo, y=dgen)) +
  theme_bw(base_family = "Arial") +
  theme(panel.grid = element_blank(),
        legend.position='none',
        axis.text = element_text(color="black", size=14),
        axis.title = element_text(size=16)) +
  stat_density_2d(aes(fill = after_stat(density)), geom = "raster", contour = FALSE) +
  scale_fill_distiller(palette= "Spectral", direction=-1) +
  geom_point(aes(shape=comparison), size = 3) +
  geom_smooth(method="lm", se=FALSE, size=0.5, col="black") +
  annotate("text", label=expression(paste(italic(r), " = 0.8091")), x = 1700, y = 0.048, size=6) +
  annotate("text", label=expression(paste(italic(P), " = 0.0003")), x = 1700, y = 0.044, size=6) +
  labs(y=expression(paste("Genetic distance (",F[ST],"/(1-",F[ST],"))")), x="Geographic distance (km)")

p_IBD

ggsave("LK_IBD.pdf", p_IBD, device=cairo_pdf, units="cm", width=15, height=15, limitsize=FALSE)

# IBC (isolation-by-climate)
## Calculate climatic distance using PCA
### Load the raster of current climatic data
ras_current <- stack("~/Dropbox/Postdoc_Milan/LK_Joan/GEA/ras_current.grd")

### Extracting environmental values for each source population from the rasters
Env <- data.frame(raster::extract(ras_current, coord[,2:3]))

### Standardization of the variables
Env <- scale(Env, center=TRUE, scale=TRUE) # center=TRUE, scale=TRUE are the defaults for scale()

### Recovering scaling coefficients
scale_env <- attr(Env, 'scaled:scale')
center_env <- attr(Env, 'scaled:center')

### Climatic table
Env <- as.data.frame(Env)
row.names(Env) <- rownames(coord)
head(Env)

res.pca <- prcomp(Env, scale = FALSE)
fviz_eig(res.pca)
res.ind <- get_pca_ind(res.pca)
pca_clim <- as.data.frame(res.ind$coord)
dclim <- dist(pca_clim)
dclim <- sort_dist_mat(dclim, by_rows = TRUE, by_cols = TRUE)

## Perform mantel test in vegan
ibc <- vegan::mantel(dgen, dclim, method = "pearson", permutations = 9999, na.rm = TRUE)

## Plot climatic distance against genetic distance
d_dgen <- dist2list(dgen, tri = TRUE)
d_dclim <- dist2list(dclim, tri = TRUE)
d <- merge(d_dgen, d_dclim, by=c("col","row"))
colnames(d) <- c("p1","p2","dgen","dclim")
EUR <- c("ESS","ESN","SIC","ITS","ITN","GRL","GRC","ISR")
ASI <- c("KAZ","MON","MOS","RUS")

d <- mutate(d, comparison = case_when(p1 %in% EUR & p2 %in% EUR ~ "within_gen_cluster",
                                      p1 %in% ASI & p2 %in% ASI ~ "within_gen_cluster",
                                      p1 %in% EUR & p2 %in% ASI ~ "between_gen_cluster",
                                      p1 %in% ASI & p2 %in% EUR ~ "between_gen_cluster"))
d$comparison <- as.factor(d$comparison)

p_IBC <- ggplot(d, aes(x=dclim, y=dgen)) +
  theme_bw(base_family = "Arial") +
  theme(panel.grid = element_blank(),
        legend.position='none',
        axis.text = element_text(color="black", size=14),
        axis.title = element_text(size=16)) +
  stat_density_2d(aes(fill = after_stat(density)), geom = "raster", contour = FALSE) +
  scale_fill_distiller(palette= "Spectral", direction=-1) +
  geom_point(aes(shape=comparison), size = 3) +
  geom_smooth(method="lm", se=FALSE, size=0.5, col="black") +
  annotate("text", label=expression(paste(italic(r), " = 0.6751")), x = 3, y = 0.048, size=6) +
  annotate("text", label=expression(paste(italic(P), " = 0.0011")), x = 3, y = 0.044, size=6) +
  labs(y=expression(paste("Genetic distance (",F[ST],"/(1-",F[ST],"))")), x="Climatic distance")

p_IBC

ggsave("LK_IBC.pdf", p_IBC, device=cairo_pdf, units="cm", width=15, height=15, limitsize=FALSE)

p_IBD_IBC <- ggarrange(p_IBD, p_IBC, nrow=1)

ggsave("LK_IBD-IBC.pdf", p_IBD_IBC, device=cairo_pdf, units="cm", width=30, height=15, limitsize=FALSE)

# Run partial mantel test with geographic distance and climatic distance
mantel.partial(dgen,dclim,dgeo, method = "pearson", permutations = 9999, na.rm = TRUE) # Climate distance is not significant when controlling for geographic distance
mantel.partial(dclim,dgen,dgeo, method = "pearson", permutations = 9999, na.rm = TRUE) # Climate distance is not significant when controlling for geographic distance

# Run MRM as partial Mantel tests have been heavily criticised
MRM(dgen ~ dgeo + dclim, nperm=9999)
MRM(dgen ~ dclim + dgeo, nperm=9999)

# Infer covariance between dgeo and dclim
cov_geo_clim <- vegan::mantel(dgeo, dclim, method = "pearson", permutations = 9999, na.rm = TRUE)
cov_clim_geo <- vegan::mantel(dclim, dgeo, method = "pearson", permutations = 9999, na.rm = TRUE)

# Calculate correlations between PC1 from population structure, PC1 on bioclimatic variables and longitude
cor.test(Variables$bio_PC1, Variables$y) #r=0.877
cor.test(Variables$bio_PC1, Variables$PC1) #r=0.9213558