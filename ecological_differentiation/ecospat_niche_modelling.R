## Load libraries
library(raster)
library(sf)
library(rgdal)
library(dplyr)
library(ggplot2)
library(ecospat)
library(ggnewscale)
library(ade4)
library(tidyterra)

## Load occurrence and background data
setwd("~/Dropbox/Postdoc_Milan/LK_Joan/niche_modelling")

### Occurrence data
#### Load raster of Chelsa data
current.list <- list.files(path="./chelsa_present", pattern =".tif", full.names=TRUE)
ras <- stack(current.list)
names(ras) <- c("bio01","bio02","bio03","bio04","bio05","bio06","bio07","bio08","bio09","bio10","bio11","bio12",
                "bio13","bio14","bio15","bio16","bio17","bio18","bio19")

#### Load occurrence data and split it into ESU/season
occurrence <- read.csv("../SDMs/occurrence_data/Table_S5_occurrence_data.csv")

west_breed <- occurrence %>% filter(ESU=="Western") %>% filter(season=="breeding")
colnames(west_breed)[c(1,2)] <- c("x","y") 
coordinates(west_breed) <- c("x","y")
west_breed <- cbind(west_breed, extract(ras, west_breed))

east_breed <- occurrence %>% filter(ESU=="Eastern") %>% filter(season=="breeding")
colnames(east_breed)[c(1,2)] <- c("x","y") 
coordinates(east_breed) <- c("x","y")
east_breed <- cbind(east_breed, extract(ras, east_breed))

west_nonbreed <- occurrence %>% filter(ESU=="Western") %>% filter(season=="non-breeding")
colnames(west_nonbreed)[c(1,2)] <- c("x","y") 
coordinates(west_nonbreed) <- c("x","y")
west_nonbreed <- cbind(west_nonbreed, extract(ras, west_nonbreed))

east_nonbreed <- occurrence %>% filter(ESU=="Eastern") %>% filter(season=="non-breeding")
colnames(east_nonbreed)[c(1,2)] <- c("x","y") 
coordinates(east_nonbreed) <- c("x","y")
east_nonbreed <- cbind(east_nonbreed, extract(ras, east_nonbreed))

#### Load background data
background_compl <- read.csv("background_complessivo.csv")
background_compl <- background_compl[,c(2,3)] # keep only coordinates
background_compl$background <- "general_background"
colnames(background_compl)[c(1,2)] <- c("x","y") 
coordinates(background_compl) <- c("x","y") # convert to spatialpointsdataframe
background_compl <- cbind(background_compl, extract(ras, background_compl)) # extract bioclimatic variables from Chelsa raster 

background_W_breed <- read.csv("background_entro_1500_breeding_west.csv")
background_W_breed <- background_W_breed[,c(2,3)]
colnames(background_W_breed)[c(1,2)] <- c("x","y")
background_W_breed$ESU <- "western"
background_W_breed$season <- "breeding"
coordinates(background_W_breed) <- c("x","y")
background_W_breed <- cbind(background_W_breed, extract(ras, background_W_breed))

background_E_breed <- read.csv("background_entro_1500_breeding_east.csv")
background_E_breed <- background_E_breed[,c(3,4)]
colnames(background_E_breed)[c(1,2)] <- c("x","y")
background_E_breed$ESU <- "eastern"
background_E_breed$season <- "breeding"
coordinates(background_E_breed) <- c("x","y")
background_E_breed <- cbind(background_E_breed, extract(ras, background_E_breed))

background_W_nonbreed <- read.csv("background_entro_1500_wintering_west.csv")
background_W_nonbreed <- background_W_nonbreed[,c(2,3)]
colnames(background_W_nonbreed)[c(1,2)] <- c("x","y")
background_W_nonbreed$ESU <- "western"
background_W_nonbreed$season <- "non-breeding"
coordinates(background_W_nonbreed) <- c("x","y")
background_W_nonbreed <- cbind(background_W_nonbreed, extract(ras, background_W_nonbreed))

background_E_nonbreed <- read.csv("background_entro_1500_wintering_east.csv")
background_E_nonbreed <- background_E_nonbreed[,c(2,3)]
colnames(background_E_nonbreed)[c(1,2)] <- c("x","y")
background_E_nonbreed$ESU <- "eastern"
background_E_nonbreed$season <- "non-breeding"
coordinates(background_E_nonbreed) <- c("x","y")
background_E_nonbreed <- cbind(background_E_nonbreed, extract(ras, background_E_nonbreed))

## Produce global environmental background data
### For breeding season
background_all_breed <- rbind(background_W_breed, background_E_breed)

### For non-breeding season
background_all_nonbreed <- rbind(background_W_nonbreed, background_E_nonbreed)

## Perform a Principal Component Analysis (PCA) of the background environmental data
### Prepare the breeding data
bioclim_background_all_breed <- background_all_breed[,c(3:21)]
bioclim_background_all_breed_mat <- as.matrix(bioclim_background_all_breed@data)
bioclim_background_all_breed_mat <- bioclim_background_all_breed_mat[complete.cases(bioclim_background_all_breed_mat), ]

### Prepare the non-breeding data
bioclim_background_all_nonbreed <- background_all_nonbreed[,c(3:21)]
bioclim_background_all_nonbreed_mat <- as.matrix(bioclim_background_all_nonbreed@data)
bioclim_background_all_nonbreed_mat <- bioclim_background_all_nonbreed_mat[complete.cases(bioclim_background_all_nonbreed_mat), ]

### Perform the PCA
#### Breeding data
pca.clim_breed <- dudi.pca(bioclim_background_all_breed_mat, center = TRUE,
                     scale = TRUE, scannf = FALSE, nf = 2)

global.scores_lk_breed <- pca.clim_breed$li

#### Non-breeding data
pca.clim_nonbreed <- dudi.pca(bioclim_background_all_nonbreed_mat, center = TRUE,
                           scale = TRUE, scannf = FALSE, nf = 2)

global.scores_lk_nonbreed <- pca.clim_nonbreed$li

### Map occurrence data to PCA space
#### Breeding data
westLS.scores_breed <- suprow(pca.clim_breed,
         data.frame(west_breed)[, colnames(bioclim_background_all_breed_mat)])$li   
eastLS.scores_breed <-suprow(pca.clim_breed,
         data.frame(east_breed)[, colnames(bioclim_background_all_breed_mat)])$li

#### Non-breeding data
westLS.scores_nonbreed <- suprow(pca.clim_nonbreed,
                              data.frame(west_nonbreed)[, colnames(bioclim_background_all_nonbreed_mat)])$li   
eastLS.scores_nonbreed <-suprow(pca.clim_nonbreed,
                             data.frame(east_nonbreed)[, colnames(bioclim_background_all_nonbreed_mat)])$li

### Map Western and Eastern background data to PCA space
#### Breeding data
bioclim_background_W_breed <- background_W_breed[,c(3:21)]
bioclim_background_W_breed_mat <- as.matrix(bioclim_background_W_breed@data)
bioclim_background_W_breed_mat <- bioclim_background_W_breed_mat[complete.cases(bioclim_background_W_breed_mat), ]

bioclim_background_E_breed <- background_E_breed[,c(3:21)]
bioclim_background_E_breed_mat <- as.matrix(bioclim_background_E_breed@data)
bioclim_background_E_breed_mat <- bioclim_background_E_breed_mat[complete.cases(bioclim_background_E_breed_mat), ]

westEnv.scores_breed <- suprow(pca.clim_breed, bioclim_background_W_breed_mat)$li
eastEnv.scores_breed <- suprow(pca.clim_breed, bioclim_background_E_breed_mat)$li

#### Non-breeding data
bioclim_background_W_nonbreed <- background_W_nonbreed[,c(3:21)]
bioclim_background_W_nonbreed_mat <- as.matrix(bioclim_background_W_nonbreed@data)
bioclim_background_W_nonbreed_mat <- bioclim_background_W_nonbreed_mat[complete.cases(bioclim_background_W_nonbreed_mat), ]

bioclim_background_E_nonbreed <- background_E_nonbreed[,c(3:21)]
bioclim_background_E_nonbreed_mat <- as.matrix(bioclim_background_E_nonbreed@data)
bioclim_background_E_nonbreed_mat <- bioclim_background_E_nonbreed_mat[complete.cases(bioclim_background_E_nonbreed_mat), ]

westEnv.scores_nonbreed <- suprow(pca.clim_nonbreed, bioclim_background_W_nonbreed_mat)$li
eastEnv.scores_nonbreed <- suprow(pca.clim_nonbreed, bioclim_background_E_nonbreed_mat)$li

## Calculate the Occurrence Density Grid for western and eastern ESUs
### Breeding
westGrid_breed <- ecospat.grid.clim.dyn(global.scores_lk_breed,
                                    westEnv.scores_breed,
                                    westLS.scores_breed)

eastGrid_breed <- ecospat.grid.clim.dyn(global.scores_lk_breed,
                                      eastEnv.scores_breed, 
                                      eastLS.scores_breed)

### Non-breeding
westGrid_nonbreed <- ecospat.grid.clim.dyn(global.scores_lk_nonbreed,
                                        westEnv.scores_nonbreed,
                                        westLS.scores_nonbreed)

eastGrid_nonbreed <- ecospat.grid.clim.dyn(global.scores_lk_nonbreed,
                                        eastEnv.scores_nonbreed, 
                                        eastLS.scores_nonbreed)

## Compare niches of both ESUs
## Default ecospat plot
### Breeding
ecospat.plot.niche.dyn(westGrid_breed, eastGrid_breed, quant = 0.1, interest = 2, name.axis1 = "PC1 (44.8 %)",
                       name.axis2 = "PC2 (25.5 %)", col.unf = "#ff6d00", col.exp = "#5c7ec0",
                       col.stab = "black", colZ1 = "#ff6d00", colZ2 = "#5c7ec0", transparency=50)

### Non-breeding
ecospat.plot.niche.dyn(westGrid_nonbreed, eastGrid_nonbreed, quant = 0.1, interest = 1, name.axis1 = "PC1 (51.8 %)",
                       name.axis2 = "PC2 (19.4 %)", col.unf = "#ff6d00", col.exp = "#5c7ec0",
                       col.stab = "black", colZ1 = "#ff6d00", colZ2 = "#5c7ec0", transparency=50)

## Prepare data for ggplot2
### Breeding
#### Western ESU
west_grid_zcor_breed <- raster(westGrid_breed$z.uncor)
west_grid_zcor_breed <- as(west_grid_zcor_breed, "SpatialPixelsDataFrame")
west_grid_zcor_breed <- as.data.frame(west_grid_zcor_breed)
west_grid_zcor_breed$val <- west_grid_zcor_breed$lyr.1/sum(west_grid_zcor_breed$lyr.1) # Normalise density

west_grid_Z_background_breed <- raster(westGrid_breed$Z)
west_grid_Z_background_breed <- as(west_grid_Z_background_breed, "SpatialPixelsDataFrame")
west_grid_Z_background_breed <- as.data.frame(west_grid_Z_background_breed)
west_background_breed_quantiles <- quantile(west_grid_Z_background_breed$lyr.1[west_grid_Z_background_breed$lyr.1>0], c(0,0.25)) 

#### Eastern ESU
east_grid_zcor_breed <- raster(eastGrid_breed$z.uncor)
east_grid_zcor_breed <- as(east_grid_zcor_breed, "SpatialPixelsDataFrame")
east_grid_zcor_breed <- as.data.frame(east_grid_zcor_breed)
east_grid_zcor_breed$val <- east_grid_zcor_breed$lyr.1/sum(east_grid_zcor_breed$lyr.1) # Normalise density

east_grid_Z_background_breed <- raster(eastGrid_breed$Z)
east_grid_Z_background_breed <- as(east_grid_Z_background_breed, "SpatialPixelsDataFrame")
east_grid_Z_background_breed <- as.data.frame(east_grid_Z_background_breed)
east_background_breed_quantiles <- quantile(east_grid_Z_background_breed$lyr.1[east_grid_Z_background_breed$lyr.1>0], c(0,0.25)) 

### Non-breeding
#### Western ESU
west_grid_zcor_nonbreed <- raster(westGrid_nonbreed$z.uncor)
west_grid_zcor_nonbreed <- as(west_grid_zcor_nonbreed, "SpatialPixelsDataFrame")
west_grid_zcor_nonbreed <- as.data.frame(west_grid_zcor_nonbreed)
west_grid_zcor_nonbreed$val <- west_grid_zcor_nonbreed$lyr.1/sum(west_grid_zcor_nonbreed$lyr.1) # Normalise density

west_grid_Z_background_nonbreed <- raster(westGrid_nonbreed$Z)
west_grid_Z_background_nonbreed <- as(west_grid_Z_background_nonbreed, "SpatialPixelsDataFrame")
west_grid_Z_background_nonbreed <- as.data.frame(west_grid_Z_background_nonbreed)
west_background_nonbreed_quantiles <- quantile(west_grid_Z_background_nonbreed$lyr.1[west_grid_Z_background_nonbreed$lyr.1>0], c(0,0.25)) 

#### Eastern ESU
east_grid_zcor_nonbreed <- raster(eastGrid_nonbreed$z.uncor)
east_grid_zcor_nonbreed <- as(east_grid_zcor_nonbreed, "SpatialPixelsDataFrame")
east_grid_zcor_nonbreed <- as.data.frame(east_grid_zcor_nonbreed)
east_grid_zcor_nonbreed$val <- east_grid_zcor_nonbreed$lyr.1/sum(east_grid_zcor_nonbreed$lyr.1) # Normalise density

east_grid_Z_background_nonbreed <- raster(eastGrid_nonbreed$Z)
east_grid_Z_background_nonbreed <- as(east_grid_Z_background_nonbreed, "SpatialPixelsDataFrame")
east_grid_Z_background_nonbreed <- as.data.frame(east_grid_Z_background_nonbreed)
east_background_nonbreed_quantiles <- quantile(east_grid_Z_background_nonbreed$lyr.1[east_grid_Z_background_nonbreed$lyr.1>0], c(0,0.25)) 

## Compare niches of both ESUs using ggplot
### Breeding
p_ecospat_niche_breed <- ggplot() +
  theme_bw(base_family = "Arial") +
  theme(panel.grid = element_blank(),
        legend.position='none',
        axis.text = element_text(color="black", size=14),
        axis.title = element_text(size=16)) +
  geom_tile(data=subset(east_grid_zcor_breed, lyr.1 > 0), aes(x=x, y=y, fill=val), alpha=0.8) +
  scale_fill_gradient(low="#cdd8ec", high="#395894", guide="none") +
  new_scale_fill() +
  geom_tile(data=subset(west_grid_zcor_breed, lyr.1 > 0), aes(x=x, y=y, fill=val), alpha=0.8) +
  scale_fill_gradient(low="#ffddc4", high="#d85c00", guide="none") +
  geom_spatraster_contour(data=eastGrid_breed$Z, breaks=c(east_background_breed_quantiles[2]), colour="#5c7ec0", linetype=2) +
  geom_spatraster_contour(data=eastGrid_breed$Z, breaks=c(east_background_breed_quantiles[1]), colour="#5c7ec0", linetype=1) +
  geom_spatraster_contour(data=westGrid_breed$Z, breaks=c(west_background_breed_quantiles[2]), colour="#ff6d00", linetype=2) +
  geom_spatraster_contour(data=westGrid_breed$Z, breaks=c(west_background_breed_quantiles[1]), colour="#ff6d00", linetype=1) +
  coord_fixed(xlim = c(min(west_grid_Z_background_breed$x),12),
              ylim = c(min(west_grid_Z_background_breed$y),18), expand = F) +
  xlab("PC1 (44.8 %)") + ylab("PC2 (25.5 %)")

p_ecospat_niche_breed  

### Non-breeding
p_ecospat_niche_nonbreed <- ggplot() +
  theme_bw(base_family = "Arial") +
  theme(panel.grid = element_blank(),
        legend.position='none',
        axis.text = element_text(color="black", size=14),
        axis.title = element_text(size=16)) +
  geom_tile(data=subset(east_grid_zcor_nonbreed, lyr.1 > 0), aes(x=x, y=y, fill=val), alpha=0.8) +
  scale_fill_gradient(low="#cdd8ec", high="#395894", guide="none") +
  new_scale_fill() +
  geom_tile(data=subset(west_grid_zcor_nonbreed, lyr.1 > 0), aes(x=x, y=y, fill=val), alpha=0.8) +
  scale_fill_gradient(low="#ffddc4", high="#d85c00", guide="none") +
  geom_spatraster_contour(data=eastGrid_nonbreed$Z, breaks=c(east_background_nonbreed_quantiles[2]), colour="#5c7ec0", linetype=2) +
  geom_spatraster_contour(data=eastGrid_nonbreed$Z, breaks=c(east_background_nonbreed_quantiles[1]), colour="#5c7ec0", linetype=1) +
  geom_spatraster_contour(data=westGrid_nonbreed$Z, breaks=c(west_background_nonbreed_quantiles[2]), colour="#ff6d00", linetype=2) +
  geom_spatraster_contour(data=westGrid_nonbreed$Z, breaks=c(west_background_nonbreed_quantiles[1]), colour="#ff6d00", linetype=1) +
  coord_fixed(xlim = c(min(west_grid_Z_background_nonbreed$x), max(west_grid_Z_background_nonbreed$x)),
              ylim = c(min(west_grid_Z_background_nonbreed$y), max(west_grid_Z_background_nonbreed$y)), expand = F) +
  xlab("PC1 (51.8 %)") + ylab("PC2 (19.4 %)")

p_ecospat_niche_nonbreed
ggsave("ecospat_niche_nonbreed.pdf", p_ecospat_niche_nonbreed, device=cairo_pdf, units="cm", width=8, height=10, limitsize=FALSE)


## Plot bioclimatic variables contribution
### Breeding
ecospat.plot.contrib(contrib=pca.clim_breed$co, eigen=pca.clim_breed$eig)
pca.clim_breed$co

### Non-breeding
ecospat.plot.contrib(contrib=pca.clim_nonbreed$co, eigen=pca.clim_nonbreed$eig)
pca.clim_nonbreed$co

## Calculate niche overlap
### Breeding
ecospat.niche.overlap(westGrid_breed, eastGrid_breed, cor=T)

#### Perform the Niche Equivalency Test
eq.test_breed <- ecospat.niche.equivalency.test(westGrid_breed, eastGrid_breed, overlap.alternative = "lower",
                                                expansion.alternative = "higher", stability.alternative = "lower",
                                                unfilling.alternative = "higher", rep = 1000, ncores = 4)

#### Perform the Niche Similarity Test
sim.test_breed <- ecospat.niche.similarity.test(westGrid_breed, eastGrid_breed, rep = 1000, rand.type = 1, ncores = 4)

#### Plot Equivalency and Similarity Test
par(mfrow=c(1,2))
ecospat.plot.overlap.test(eq.test_breed, "D", "Equivalency") 
ecospat.plot.overlap.test(sim.test_breed, "D", "Similarity")

### Non-breeding
ecospat.niche.overlap(westGrid_nonbreed, eastGrid_nonbreed, cor=T)

#### Perform the Niche Equivalency Test
eq.test_nonbreed <- ecospat.niche.equivalency.test(westGrid_nonbreed, eastGrid_nonbreed, overlap.alternative = "lower",
                                                   expansion.alternative = "higher", stability.alternative = "lower",
                                                   unfilling.alternative = "higher", rep = 1000, ncores = 4)

#### Perform the Niche Similarity Test
sim.test_nonbreed <- ecospat.niche.similarity.test(westGrid_nonbreed, eastGrid_nonbreed, rep = 1000, rand.type = 1, ncores = 4)

#### Plot Equivalency and Similarity Test
par(mfrow=c(1,2))
ecospat.plot.overlap.test(eq.test_nonbreed, "D", "Equivalency") 
ecospat.plot.overlap.test(sim.test_nonbreed, "D", "Similarity")

save.image("~/Dropbox/Postdoc_Milan/LK_Joan/niche_modelling/LK_ecospat_niche_modelling.RData")
###########################
########## Extra ##########
###########################

# gridding the native niche
grid.clim.t.west <- ecospat.grid.clim.dyn(glob = background_all_breed@data[,3],
                                         glob1 = data.frame(background_W_breed@data[,3]),
                                         data.frame(west_breed@data)[,3], R = 1000, th.sp = 0)

# gridding the invasive niche
grid.clim.t.east <- ecospat.grid.clim.dyn (glob = background_all_breed@data[,3], 
                                          glob1 = data.frame(background_E_breed@data[,3]), 
                                          data.frame(east_breed@data)[,3], R = 1000, th.sp = 0)

t.dyn <- ecospat.niche.dyn.index (grid.clim.t.west, grid.clim.t.east, intersection=0.1)

ecospat.plot.niche.dyn(grid.clim.t.west, grid.clim.t.east, quant=0.1, interest=2,
                       name.axis1="Annual mean temperature")

# showing the shift of the niche centroid along the temperature gradient (compared to the shift of the available climate in the study area)
ecospat.shift.centroids(data.frame(ocEUR)[,4],
                        data.frame(ocAUS)[,4],
                        data.frame(eurEnvM)[,1],
                        data.frame(ausEnvM)[,1])
