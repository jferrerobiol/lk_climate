##############################################
## Script to calculate land cover use by LK ##
##############################################

## The first part is done in the cluster as it needs a lot of memory
## Load libraries
library(raster)
library(terra)
library(tidyverse)
library(tidyr)

## Load tifs of land cover use
setwd("/home/users/joan.ferrer/GRILLAIO2/GRILLAIO/LK_genomics/land_use")

### Breeding
ras.list_breed <- list.files(path="./breeding", pattern =".tif", full.names=TRUE)
for (i in 1:length(ras.list_breed)) {
  ras <- stack(ras.list_breed[i])
  assign(paste0("ras_breed_", i), ras)
}

## Merge tiffs of land cover use to build a single raster
ras_breed <- raster::mosaic(ras_breed_1,ras_breed_2,ras_breed_3,ras_breed_4,ras_breed_5,ras_breed_6,ras_breed_7,
                            ras_breed_8,ras_breed_9,ras_breed_10,ras_breed_11,ras_breed_12,ras_breed_13,ras_breed_14,
                            ras_breed_15,ras_breed_16,ras_breed_17,ras_breed_18,ras_breed_19,ras_breed_20,ras_breed_21,
                            ras_breed_22,ras_breed_23,fun="mean")
print("merging of tiffs done")
## Convert to spatraster
ras_breed_spatraster <- terra::rast(ras_breed)
## Split land use categories to individual layers
ras_breed_single_layers <- terra::segregate(ras_breed_spatraster)
print("split land use to individual layers done")
## Merge land use categories
closed_forest <- ras_breed_single_layers$`111` + ras_breed_single_layers$`113` + ras_breed_single_layers$`112` + ras_breed_single_layers$`114` + ras_breed_single_layers$`115`+ ras_breed_single_layers$`116`
open_forest <- ras_breed_single_layers$`121` + ras_breed_single_layers$`123` + ras_breed_single_layers$`122` + ras_breed_single_layers$`124` + ras_breed_single_layers$`125`+ ras_breed_single_layers$`126`
shrubs <- ras_breed_single_layers$`20`
herb_veg <- ras_breed_single_layers$`30`
herb_wet <- ras_breed_single_layers$`90`
moss_lich <- ras_breed_single_layers$`100`
bare <- ras_breed_single_layers$`60`
cropland <- ras_breed_single_layers$`40`
urban <- ras_breed_single_layers$`50`
snow_ice <- ras_breed_single_layers$`70`
perm_water <- ras_breed_single_layers$`80`
open_sea <- ras_breed_single_layers$`200`
## Create raster of reclassified land use categories
ras_breed_combined_layers <- c(closed_forest, open_forest, shrubs, herb_veg, herb_wet, moss_lich, bare, cropland, urban, snow_ice, perm_water, open_sea)
print("raster of reclassified land use categories done")
## Load raster at 2.5 arcmin resolution
ras_to_right_res <- raster::stack("/home/users/joan.ferrer/GRILLAIO2/GRILLAIO/LK_genomics/land_use/stack_breeding_east_chelsa_past_recl.grd")
## Convert it to spatraster
ras_to_right_res <- terra::rast(ras_to_right_res)
## Change resolution of land use data raster to 2.5 arcmin
ras_breed_combined_layers_rightres <- terra::resample(ras_breed_combined_layers, ras_to_right_res)
print("resolution changed to 2.5 arcmin")
## Write the land use raster at 2.5 arcmin resolution 
terra::writeRaster(ras_breed_combined_layers_rightres, "/home/users/joan.ferrer/GRILLAIO2/GRILLAIO/LK_genomics/land_use/ras_breed_combined_layers_rightres.grd")
## Convert to dataframe
d_ras_breed_combined_layers_rightres <- terra::as.data.frame(ras_breed_combined_layers_rightres, xy=T, cells=T)
print("raster converted to data frame")
## Write dataframe
write.table(d_ras_breed_combined_layers_rightres, "/home/users/joan.ferrer/GRILLAIO2/GRILLAIO/LK_genomics/land_use/d_ras_breed_combined_layers_rightres.tsv", sep ="\t", quote=F, row.names=F)
print("finished")

### Wintering
ras.list_wint <- list.files(path="./wintering", pattern =".tif", full.names=TRUE)
for (i in 1:length(ras.list_wint)) {
  ras <- stack(ras.list_wint[i])
  assign(paste0("ras_wint_", i), ras)
}

## Merge tiffs of land cover use to build a single raster
ras_wint <- raster::mosaic(ras_wint_1,ras_wint_2,ras_wint_3,ras_wint_4,ras_wint_5,ras_wint_6,ras_wint_7,
                           ras_wint_8,ras_wint_9,ras_wint_10,ras_wint_11,ras_wint_12,ras_wint_13,ras_wint_14,
                           fun="mean")
print("merging of tiffs done")
## Convert to spatraster
ras_wint_spatraster <- terra::rast(ras_wint)
## Split land use categories to individual layers
ras_wint_single_layers <- terra::segregate(ras_wint_spatraster)
print("split land use to individual layers done")
#terra::writeRaster(ras_wint_single_layers, "/home/users/joan.ferrer/GRILLAIO2/GRILLAIO/LK_genomics/land_use/ras_wint_ind_layers_highres.grd")
## Merge land use categories
closed_forest <- ras_wint_single_layers$`111` + ras_wint_single_layers$`112` + ras_wint_single_layers$`114` + ras_wint_single_layers$`115`+ ras_wint_single_layers$`116`
open_forest <- ras_wint_single_layers$`121` + ras_wint_single_layers$`122` + ras_wint_single_layers$`124` + ras_wint_single_layers$`125`+ ras_wint_single_layers$`126`
shrubs <- ras_wint_single_layers$`20`
herb_veg <- ras_wint_single_layers$`30`
herb_wet <- ras_wint_single_layers$`90`
moss_lich <- ras_wint_single_layers$`100`
bare <- ras_wint_single_layers$`60`
cropland <- ras_wint_single_layers$`40`
urban <- ras_wint_single_layers$`50`
snow_ice <- ras_wint_single_layers$`70`
perm_water <- ras_wint_single_layers$`80`
open_sea <- ras_wint_single_layers$`200`
## Create raster of reclassified land use categories
ras_wint_combined_layers <- c(closed_forest, open_forest, shrubs, herb_veg, herb_wet, moss_lich, bare, cropland, urban, snow_ice, perm_water, open_sea)
print("raster of reclassified land use categories done")
## Load raster at 2.5 arcmin resolution
ras_to_right_res <- raster::stack("/home/users/joan.ferrer/GRILLAIO2/GRILLAIO/LK_genomics/land_use/stack_breeding_east_chelsa_past_recl.grd")
## Convert it to spatraster
ras_to_right_res <- terra::rast(ras_to_right_res)
## Change resolution of land use data raster to 2.5 arcmin
ras_wint_combined_layers_rightres <- terra::resample(ras_wint_combined_layers, ras_to_right_res)
print("resolution changed to 2.5 arcmin")
## Write the land use raster at 2.5 arcmin resolution 
terra::writeRaster(ras_wint_combined_layers_rightres, "/home/users/joan.ferrer/GRILLAIO2/GRILLAIO/LK_genomics/land_use/ras_wint_combined_layers_rightres.grd")
## Convert to dataframe
d_ras_wint_combined_layers_rightres <- terra::as.data.frame(ras_wint_combined_layers_rightres, xy=T, cells=T)
print("raster converted to data frame")
## Write dataframe
write.table(d_ras_wint_combined_layers_rightres, "/home/users/joan.ferrer/GRILLAIO2/GRILLAIO/LK_genomics/land_use/d_ras_wint_combined_layers_rightres.tsv", sep ="\t", quote=F, row.names=F)
print("finished")

## Load rasters to extract data from occurrence and background
ras_breed_combined_layers_rightres <- terra::rast("ras_breed_combined_layers_rightres.grd")
ras_wint_combined_layers_rightres <- terra::rast("ras_wint_combined_layers_rightres.grd")

## Load presence data for both clusters and both seasons
occurrence <- read.csv("Table_S5_occurrence_data.csv")

west_breed <- occurrence %>% filter(ESU=="Western") %>% filter(season=="breeding")
colnames(west_breed)[c(1,2)] <- c("x","y") 

east_breed <- occurrence %>% filter(ESU=="Eastern") %>% filter(season=="breeding")
colnames(east_breed)[c(1,2)] <- c("x","y") 

west_nonbreed <- occurrence %>% filter(ESU=="Western") %>% filter(season=="non-breeding")
colnames(west_nonbreed)[c(1,2)] <- c("x","y") 

east_nonbreed <- occurrence %>% filter(ESU=="Eastern") %>% filter(season=="non-breeding")
colnames(east_nonbreed)[c(1,2)] <- c("x","y") 

## Load background data for both clusters and both seasons
background_W_breed <- read.csv("background_entro_1500_breeding_west.csv")
background_W_breed <- background_W_breed[,c(2,3)]
colnames(background_W_breed)[c(1,2)] <- c("x","y")
background_W_breed$ESU <- "Western"
background_W_breed$season <- "breeding"

background_E_breed <- read.csv("background_entro_1500_breeding_east.csv")
background_E_breed <- background_E_breed[,c(3,4)]
colnames(background_E_breed)[c(1,2)] <- c("x","y")
background_E_breed$ESU <- "Eastern"
background_E_breed$season <- "breeding"

background_W_nonbreed <- read.csv("background_entro_1500_wintering_west.csv")
background_W_nonbreed <- background_W_nonbreed[,c(2,3)]
colnames(background_W_nonbreed)[c(1,2)] <- c("x","y")
background_W_nonbreed$ESU <- "Western"
background_W_nonbreed$season <- "non-breeding"

background_E_nonbreed <- read.csv("background_entro_1500_wintering_east.csv")
background_E_nonbreed <- background_E_nonbreed[,c(2,3)]
colnames(background_E_nonbreed)[c(1,2)] <- c("x","y")
background_E_nonbreed$ESU <- "Eastern"
background_E_nonbreed$season <- "non-breeding"

## Extract raster cells for each background point
### For breeding Europe
background_breed_Eur <- terra::extract(ras_breed_combined_layers_rightres, background_W_breed[,1:2], xy=T, cells=T)
colnames(background_breed_Eur)[c(2:13)] <- c("closed_forest","open_forest","shrubs","herb_veg","herb_wet","moss_lich","bare","cropland","urban","snow_ice","perm_water","open_sea")
background_breed_Eur <- na.omit(background_breed_Eur)
#### show totals
totals_background_breed_Eur_abs <- c(sum(background_breed_Eur$closed_forest)/length(background_breed_Eur$urban),
                                     sum(background_breed_Eur$open_forest)/length(background_breed_Eur$urban),
                                     sum(background_breed_Eur$shrubs)/length(background_breed_Eur$urban),
                                     sum(background_breed_Eur$herb_veg)/length(background_breed_Eur$urban),
                                     sum(background_breed_Eur$herb_wet)/length(background_breed_Eur$urban),
                                     sum(background_breed_Eur$moss_lich)/length(background_breed_Eur$urban),
                                     sum(background_breed_Eur$bare)/length(background_breed_Eur$urban),
                                     sum(background_breed_Eur$cropland)/length(background_breed_Eur$urban),
                                     sum(background_breed_Eur$urban)/length(background_breed_Eur$urban),
                                     sum(background_breed_Eur$snow_ice)/length(background_breed_Eur$urban),
                                     sum(background_breed_Eur$perm_water)/length(background_breed_Eur$urban),
                                     sum(background_breed_Eur$open_sea)/length(background_breed_Eur$urban))
names(totals_background_breed_Eur_abs) <- c("background_closed_forest","background_open_forest","background_shrubs","background_herb_veg","background_herb_wet",
                                            "background_moss_lich","background_bare","background_cropland","background_urban","background_snow_ice",
                                            "background_perm_water","background_open_sea")

### For breeding Asia
background_breed_Asi <- terra::extract(ras_breed_combined_layers_rightres, background_E_breed[,1:2], xy=T, cells=T)
colnames(background_breed_Asi)[c(2:13)] <- c("closed_forest","open_forest","shrubs","herb_veg","herb_wet","moss_lich","bare","cropland","urban","snow_ice","perm_water","open_sea")
background_breed_Asi <- na.omit(background_breed_Asi)
#### show totals
totals_background_breed_Asi_abs <- c(sum(background_breed_Asi$closed_forest)/length(background_breed_Asi$urban),
                                     sum(background_breed_Asi$open_forest)/length(background_breed_Asi$urban),
                                     sum(background_breed_Asi$shrubs)/length(background_breed_Asi$urban),
                                     sum(background_breed_Asi$herb_veg)/length(background_breed_Asi$urban),
                                     sum(background_breed_Asi$herb_wet)/length(background_breed_Asi$urban),
                                     sum(background_breed_Asi$moss_lich)/length(background_breed_Asi$urban),
                                     sum(background_breed_Asi$bare)/length(background_breed_Asi$urban),
                                     sum(background_breed_Asi$cropland)/length(background_breed_Asi$urban),
                                     sum(background_breed_Asi$urban)/length(background_breed_Asi$urban),
                                     sum(background_breed_Asi$snow_ice)/length(background_breed_Asi$urban),
                                     sum(background_breed_Asi$perm_water)/length(background_breed_Asi$urban),
                                     sum(background_breed_Asi$open_sea)/length(background_breed_Asi$urban))
names(totals_background_breed_Asi_abs) <- c("background_closed_forest","background_open_forest","background_shrubs","background_herb_veg","background_herb_wet",
                                            "background_moss_lich","background_bare","background_cropland","background_urban","background_snow_ice",
                                            "background_perm_water","background_open_sea")

### For non-breeding Europe
background_winter_Eur <- terra::extract(ras_wint_combined_layers_rightres, background_W_nonbreed[,1:2], xy=T, cells=T)
colnames(background_winter_Eur)[c(2:13)] <- c("closed_forest","open_forest","shrubs","herb_veg","herb_wet","moss_lich","bare","cropland","urban","snow_ice","perm_water","open_sea")
background_winter_Eur <- na.omit(background_winter_Eur)
#### show totals
totals_background_winter_Eur_abs <- c(sum(background_winter_Eur$closed_forest)/length(background_winter_Eur$urban),
                                      sum(background_winter_Eur$open_forest)/length(background_winter_Eur$urban),
                                      sum(background_winter_Eur$shrubs)/length(background_winter_Eur$urban),
                                      sum(background_winter_Eur$herb_veg)/length(background_winter_Eur$urban),
                                      sum(background_winter_Eur$herb_wet)/length(background_winter_Eur$urban),
                                      sum(background_winter_Eur$moss_lich)/length(background_winter_Eur$urban),
                                      sum(background_winter_Eur$bare)/length(background_winter_Eur$urban),
                                      sum(background_winter_Eur$cropland)/length(background_winter_Eur$urban),
                                      sum(background_winter_Eur$urban)/length(background_winter_Eur$urban),
                                      sum(background_winter_Eur$snow_ice)/length(background_winter_Eur$urban),
                                      sum(background_winter_Eur$perm_water)/length(background_winter_Eur$urban),
                                      sum(background_winter_Eur$open_sea)/length(background_winter_Eur$urban))
names(totals_background_winter_Eur_abs) <- c("background_closed_forest","background_open_forest","background_shrubs","background_herb_veg","background_herb_wet",
                                             "background_moss_lich","background_bare","background_cropland","background_urban","background_snow_ice",
                                             "background_perm_water","background_open_sea")

### For non-breeding Asia
background_winter_Asi <- terra::extract(ras_wint_combined_layers_rightres, background_E_nonbreed[,1:2], xy=T, cells=T)
colnames(background_winter_Asi)[c(2:13)] <- c("closed_forest","open_forest","shrubs","herb_veg","herb_wet","moss_lich","bare","cropland","urban","snow_ice","perm_water","open_sea")
background_winter_Asi <- na.omit(background_winter_Asi)
#### show totals
totals_background_winter_Asi_abs <- c(sum(background_winter_Asi$closed_forest)/length(background_winter_Asi$urban),
                                      sum(background_winter_Asi$open_forest)/length(background_winter_Asi$urban),
                                      sum(background_winter_Asi$shrubs)/length(background_winter_Asi$urban),
                                      sum(background_winter_Asi$herb_veg)/length(background_winter_Asi$urban),
                                      sum(background_winter_Asi$herb_wet)/length(background_winter_Asi$urban),
                                      sum(background_winter_Asi$moss_lich)/length(background_winter_Asi$urban),
                                      sum(background_winter_Asi$bare)/length(background_winter_Asi$urban),
                                      sum(background_winter_Asi$cropland)/length(background_winter_Asi$urban),
                                      sum(background_winter_Asi$urban)/length(background_winter_Asi$urban),
                                      sum(background_winter_Asi$snow_ice)/length(background_winter_Asi$urban),
                                      sum(background_winter_Asi$perm_water)/length(background_winter_Asi$urban),
                                      sum(background_winter_Asi$open_sea)/length(background_winter_Asi$urban))
names(totals_background_winter_Asi_abs) <- c("background_closed_forest","background_open_forest","background_shrubs","background_herb_veg","background_herb_wet",
                                             "background_moss_lich","background_bare","background_cropland","background_urban","background_snow_ice",
                                             "background_perm_water","background_open_sea")

## Extract raster cells for each occurrence point
### For breeding Europe
presence_breed_Eur <- terra::extract(ras_breed_combined_layers_rightres, west_breed[,1:2], xy=T, cells=T)
colnames(presence_breed_Eur)[c(2:13)] <- c("closed_forest","open_forest","shrubs","herb_veg","herb_wet","moss_lich","bare","cropland","urban","snow_ice","perm_water","open_sea")
presence_breed_Eur <- presence_breed_Eur %>% mutate(background_closed_forest=totals_background_breed_Eur_abs[1]) %>%
  mutate(background_open_forest=totals_background_breed_Eur_abs[2]) %>%
  mutate(background_shrubs=totals_background_breed_Eur_abs[3]) %>%
  mutate(background_herb_veg=totals_background_breed_Eur_abs[4]) %>%
  mutate(background_herb_wet=totals_background_breed_Eur_abs[5]) %>%
  mutate(background_moss_lich=totals_background_breed_Eur_abs[6]) %>%
  mutate(background_bare=totals_background_breed_Eur_abs[7]) %>%
  mutate(background_cropland=totals_background_breed_Eur_abs[8]) %>%
  mutate(background_urban=totals_background_breed_Eur_abs[9]) %>%
  mutate(background_snow_ice=totals_background_breed_Eur_abs[10]) %>%
  mutate(background_perm_water=totals_background_breed_Eur_abs[11]) %>%
  mutate(background_open_sea=totals_background_breed_Eur_abs[12])
write.table(presence_breed_Eur, "/home/users/joan.ferrer/GRILLAIO2/GRILLAIO/LK_genomics/land_use/land_use_breed_Western.csv", sep =",", quote=F, row.names=F)

#### show totals
totals_presence_breed_Eur_abs <- c(sum(presence_breed_Eur$closed_forest)/length(presence_breed_Eur$urban),
                                   sum(presence_breed_Eur$open_forest)/length(presence_breed_Eur$urban),
                                   sum(presence_breed_Eur$shrubs)/length(presence_breed_Eur$urban),
                                   sum(presence_breed_Eur$herb_veg)/length(presence_breed_Eur$urban),
                                   sum(presence_breed_Eur$herb_wet)/length(presence_breed_Eur$urban),
                                   sum(presence_breed_Eur$moss_lich)/length(presence_breed_Eur$urban),
                                   sum(presence_breed_Eur$bare)/length(presence_breed_Eur$urban),
                                   sum(presence_breed_Eur$cropland)/length(presence_breed_Eur$urban),
                                   sum(presence_breed_Eur$urban)/length(presence_breed_Eur$urban),
                                   sum(presence_breed_Eur$snow_ice)/length(presence_breed_Eur$urban),
                                   sum(presence_breed_Eur$perm_water)/length(presence_breed_Eur$urban),
                                   sum(presence_breed_Eur$open_sea)/length(presence_breed_Eur$urban))
names(totals_presence_breed_Eur_abs) <- c("closed_forest","open_forest","shrubs","herb_veg","herb_wet","moss_lich","bare","cropland","urban","snow_ice","perm_water","open_sea")

### For breeding Asia
presence_breed_Asi <- terra::extract(ras_breed_combined_layers_rightres, east_breed[,1:2], xy=T, cells=T)
#presence_breed_Asi <- extract(ras_breed, east_breed[,1:2], cellnumbers=T, buffer=5000)
#presence_breed_Asi <- as.data.frame(do.call(rbind, presence_breed_Asi))
#presence_breed_Asi$value <- as.factor(presence_breed_Asi$value)
colnames(presence_breed_Asi)[c(2:13)] <- c("closed_forest","open_forest","shrubs","herb_veg","herb_wet","moss_lich","bare","cropland","urban","snow_ice","perm_water","open_sea")
presence_breed_Asi <- presence_breed_Asi %>% mutate(background_closed_forest=totals_background_breed_Asi_abs[1]) %>%
  mutate(background_open_forest=totals_background_breed_Asi_abs[2]) %>%
  mutate(background_shrubs=totals_background_breed_Asi_abs[3]) %>%
  mutate(background_herb_veg=totals_background_breed_Asi_abs[4]) %>%
  mutate(background_herb_wet=totals_background_breed_Asi_abs[5]) %>%
  mutate(background_moss_lich=totals_background_breed_Asi_abs[6]) %>%
  mutate(background_bare=totals_background_breed_Asi_abs[7]) %>%
  mutate(background_cropland=totals_background_breed_Asi_abs[8]) %>%
  mutate(background_urban=totals_background_breed_Asi_abs[9]) %>%
  mutate(background_snow_ice=totals_background_breed_Asi_abs[10]) %>%
  mutate(background_perm_water=totals_background_breed_Asi_abs[11]) %>%
  mutate(background_open_sea=totals_background_breed_Asi_abs[12])
write.table(presence_breed_Asi, "/home/users/joan.ferrer/GRILLAIO2/GRILLAIO/LK_genomics/land_use/land_use_breed_Eastern.csv", sep =",", quote=F, row.names=F)

#### show totals
totals_presence_breed_Asi_abs <- c(sum(presence_breed_Asi$closed_forest)/length(presence_breed_Asi$urban),
                                   sum(presence_breed_Asi$open_forest)/length(presence_breed_Asi$urban),
                                   sum(presence_breed_Asi$shrubs)/length(presence_breed_Asi$urban),
                                   sum(presence_breed_Asi$herb_veg)/length(presence_breed_Asi$urban),
                                   sum(presence_breed_Asi$herb_wet)/length(presence_breed_Asi$urban),
                                   sum(presence_breed_Asi$moss_lich)/length(presence_breed_Asi$urban),
                                   sum(presence_breed_Asi$bare)/length(presence_breed_Asi$urban),
                                   sum(presence_breed_Asi$cropland)/length(presence_breed_Asi$urban),
                                   sum(presence_breed_Asi$urban)/length(presence_breed_Asi$urban),
                                   sum(presence_breed_Asi$snow_ice)/length(presence_breed_Asi$urban),
                                   sum(presence_breed_Asi$perm_water)/length(presence_breed_Asi$urban),
                                   sum(presence_breed_Asi$open_sea)/length(presence_breed_Asi$urban))
names(totals_presence_breed_Asi_abs) <- c("closed_forest","open_forest","shrubs","herb_veg","herb_wet","moss_lich","bare","cropland","urban","snow_ice","perm_water","open_sea")

### For wintering Europe
presence_winter_Eur <- terra::extract(ras_wint_combined_layers_rightres, west_nonbreed[,1:2], xy=T, cells=T)
colnames(presence_winter_Eur)[c(2:13)] <- c("closed_forest","open_forest","shrubs","herb_veg","herb_wet","moss_lich","bare","cropland","urban","snow_ice","perm_water","open_sea")
presence_winter_Eur <- presence_winter_Eur %>% mutate(background_closed_forest=totals_background_winter_Eur_abs[1]) %>%
  mutate(background_open_forest=totals_background_winter_Eur_abs[2]) %>%
  mutate(background_shrubs=totals_background_winter_Eur_abs[3]) %>%
  mutate(background_herb_veg=totals_background_winter_Eur_abs[4]) %>%
  mutate(background_herb_wet=totals_background_winter_Eur_abs[5]) %>%
  mutate(background_moss_lich=totals_background_winter_Eur_abs[6]) %>%
  mutate(background_bare=totals_background_winter_Eur_abs[7]) %>%
  mutate(background_cropland=totals_background_winter_Eur_abs[8]) %>%
  mutate(background_urban=totals_background_winter_Eur_abs[9]) %>%
  mutate(background_snow_ice=totals_background_winter_Eur_abs[10]) %>%
  mutate(background_perm_water=totals_background_winter_Eur_abs[11]) %>%
  mutate(background_open_sea=totals_background_winter_Eur_abs[12])
write.table(presence_winter_Eur, "/home/users/joan.ferrer/GRILLAIO2/GRILLAIO/LK_genomics/land_use/land_use_winter_Western.csv", sep =",", quote=F, row.names=F)

#### show totals
totals_presence_winter_Eur_abs <- c(sum(presence_winter_Eur$closed_forest)/length(presence_winter_Eur$urban),
                                   sum(presence_winter_Eur$open_forest)/length(presence_winter_Eur$urban),
                                   sum(presence_winter_Eur$shrubs)/length(presence_winter_Eur$urban),
                                   sum(presence_winter_Eur$herb_veg)/length(presence_winter_Eur$urban),
                                   sum(presence_winter_Eur$herb_wet)/length(presence_winter_Eur$urban),
                                   sum(presence_winter_Eur$moss_lich)/length(presence_winter_Eur$urban),
                                   sum(presence_winter_Eur$bare)/length(presence_winter_Eur$urban),
                                   sum(presence_winter_Eur$cropland)/length(presence_winter_Eur$urban),
                                   sum(presence_winter_Eur$urban)/length(presence_winter_Eur$urban),
                                   sum(presence_winter_Eur$snow_ice)/length(presence_winter_Eur$urban),
                                   sum(presence_winter_Eur$perm_water)/length(presence_winter_Eur$urban),
                                   sum(presence_winter_Eur$open_sea)/length(presence_winter_Eur$urban))
names(totals_presence_winter_Eur_abs) <- c("closed_forest","open_forest","shrubs","herb_veg","herb_wet","moss_lich","bare","cropland","urban","snow_ice","perm_water","open_sea")

### For wintering Asia
presence_winter_Asi <- terra::extract(ras_wint_combined_layers_rightres, east_nonbreed[,1:2], xy=T, cells=T)
colnames(presence_winter_Asi)[c(2:13)] <- c("closed_forest","open_forest","shrubs","herb_veg","herb_wet","moss_lich","bare","cropland","urban","snow_ice","perm_water","open_sea")
presence_winter_Asi <- presence_winter_Asi %>% mutate(background_closed_forest=totals_background_winter_Asi_abs[1]) %>%
  mutate(background_open_forest=totals_background_winter_Asi_abs[2]) %>%
  mutate(background_shrubs=totals_background_winter_Asi_abs[3]) %>%
  mutate(background_herb_veg=totals_background_winter_Asi_abs[4]) %>%
  mutate(background_herb_wet=totals_background_winter_Asi_abs[5]) %>%
  mutate(background_moss_lich=totals_background_winter_Asi_abs[6]) %>%
  mutate(background_bare=totals_background_winter_Asi_abs[7]) %>%
  mutate(background_cropland=totals_background_winter_Asi_abs[8]) %>%
  mutate(background_urban=totals_background_winter_Asi_abs[9]) %>%
  mutate(background_snow_ice=totals_background_winter_Asi_abs[10]) %>%
  mutate(background_perm_water=totals_background_winter_Asi_abs[11]) %>%
  mutate(background_open_sea=totals_background_winter_Asi_abs[12])
write.table(presence_winter_Asi, "/home/users/joan.ferrer/GRILLAIO2/GRILLAIO/LK_genomics/land_use/land_use_winter_Eastern.csv", sep =",", quote=F, row.names=F)

#### show totals
totals_presence_winter_Asi_abs <- c(sum(presence_winter_Asi$closed_forest)/length(presence_winter_Asi$urban),
                                    sum(presence_winter_Asi$open_forest)/length(presence_winter_Asi$urban),
                                    sum(presence_winter_Asi$shrubs)/length(presence_winter_Asi$urban),
                                    sum(presence_winter_Asi$herb_veg)/length(presence_winter_Asi$urban),
                                    sum(presence_winter_Asi$herb_wet)/length(presence_winter_Asi$urban),
                                    sum(presence_winter_Asi$moss_lich)/length(presence_winter_Asi$urban),
                                    sum(presence_winter_Asi$bare)/length(presence_winter_Asi$urban),
                                    sum(presence_winter_Asi$cropland)/length(presence_winter_Asi$urban),
                                    sum(presence_winter_Asi$urban)/length(presence_winter_Asi$urban),
                                    sum(presence_winter_Asi$snow_ice)/length(presence_winter_Asi$urban),
                                    sum(presence_winter_Asi$perm_water)/length(presence_winter_Asi$urban),
                                    sum(presence_winter_Asi$open_sea)/length(presence_winter_Asi$urban))
names(totals_presence_winter_Asi_abs) <- c("closed_forest","open_forest","shrubs","herb_veg","herb_wet","moss_lich","bare","cropland","urban","snow_ice","perm_water","open_sea")

## Prepare dataset for plotting
### Occurrence data
d_land_use_percent_wide <- totals_presence_breed_Eur_abs %>% bind_rows(totals_presence_breed_Asi_abs) %>%
  bind_rows(totals_presence_winter_Eur_abs) %>% bind_rows(totals_presence_winter_Asi_abs) %>%
  mutate(ESU_season=c("Western_breed","Eastern_breed","Western_winter","Eastern_winter"))

d_land_use_percent_long <- d_land_use_percent_wide %>% pivot_longer(cols = !ESU_season)
write.table(d_land_use_percent_long, "/home/users/joan.ferrer/GRILLAIO2/GRILLAIO/LK_genomics/land_use/d_occurrence_land_use_percent_long.csv", sep =",", quote=F, row.names=F)

### Background data
d_land_use_background_percent_wide <- totals_background_breed_Eur_abs %>% bind_rows(totals_background_breed_Asi_abs) %>%
  bind_rows(totals_background_winter_Eur_abs) %>% bind_rows(totals_background_winter_Asi_abs) %>%
  mutate(ESU_season=c("Western_breed","Eastern_breed","Western_winter","Eastern_winter"))

d_land_use_background_percent_long <- d_land_use_background_percent_wide %>% pivot_longer(cols = !ESU_season)
write.table(d_land_use_background_percent_long, "/home/users/joan.ferrer/GRILLAIO2/GRILLAIO/LK_genomics/land_use/d_background_land_use_percent_long.csv", sep =",", quote=F, row.names=F)

## From here I do it again in my local laptop
## Plot occurrence and background land use categories with barplot
### Occurrence
setwd("~/Dropbox/Postdoc_Milan/LK_Joan/land_cover_use/")

d_land_use_percent_long <- read.csv("d_occurrence_land_use_percent_long.csv")

d_land_use_percent_long$ESU_season <- factor(d_land_use_percent_long$ESU_season,
                                             levels=c("Western_breed","Eastern_breed","Western_winter","Eastern_winter"))
d_land_use_percent_long$cat <- ifelse(d_land_use_percent_long$category=="shrubs" | d_land_use_percent_long$category=="bare" |
                                        d_land_use_percent_long$category=="cropland" | d_land_use_percent_long$category=="herb_veg" |
                                        d_land_use_percent_long$category=="herb_wet" | d_land_use_percent_long$category=="urban" |
                                        d_land_use_percent_long$category=="open_forest" | d_land_use_percent_long$category=="closed_forest",
                                      d_land_use_percent_long$category, "other")
d_land_use_percent_long$cat <- recode_factor(d_land_use_percent_long$cat, herb_wet="herb_veg", open_forest="forest", closed_forest="forest")
d_land_use_percent_long$cat <- factor(d_land_use_percent_long$cat,
                                             levels=c("forest","shrubs","herb_veg","bare","cropland","urban","other"))

#colours <- c("#007800","#706629","#FFCC33","#FFFF4C","#0096a0","#fae6a0","#b4b4b4","#f096ff","#FF0000","#f0f0f0","#0032C8","#000080")
colours <- c("#007800","#FFCC33","#FFFF4C","#b4b4b4","#f096ff","#FF0000","#000000")

p_land_use <- ggplot(d_land_use_percent_long, aes(x=ESU_season, y=proportion*100, fill=cat)) + 
  theme_bw(base_family = "Arial") +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color="black", size=14),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title = element_text(size=16)) +
  geom_col() +
  scale_fill_manual(values=colours, name="Land cover") +
  scale_x_discrete(labels=c("Breeding Western","Breeding Eastern","Non-breeding Western","Non-breeding Eastern")) +
  xlab("") + ylab("Proportion (%)")

p_land_use
ggsave("LK_land_use_occurrence.pdf", p_land_use, device=cairo_pdf, units="cm", width=15, height=15, limitsize=FALSE)

### Background
d_land_use_background_percent_long <- read.csv("d_background_land_use_percent_long.csv")

d_land_use_background_percent_long$ESU_season <- factor(d_land_use_background_percent_long$ESU_season,
                                             levels=c("Western_breed","Eastern_breed","Western_winter","Eastern_winter"))
d_land_use_background_percent_long$cat <- ifelse(d_land_use_background_percent_long$category=="background_shrubs" | d_land_use_background_percent_long$category=="background_bare" | d_land_use_background_percent_long$category=="background_cropland" |
                                                   d_land_use_background_percent_long$category=="background_herb_veg" | d_land_use_background_percent_long$category=="background_herb_wet" | d_land_use_background_percent_long$category=="background_urban" |
                                                   d_land_use_background_percent_long$category=="background_open_forest" | d_land_use_background_percent_long$category=="background_closed_forest",
                                                   d_land_use_background_percent_long$category, "background_other")
d_land_use_background_percent_long$cat <- recode_factor(d_land_use_background_percent_long$cat, background_herb_wet="background_herb_veg", background_open_forest="background_forest", background_closed_forest="background_forest")
d_land_use_background_percent_long$cat <- factor(d_land_use_background_percent_long$cat,
                                      levels=c("background_forest","background_shrubs","background_herb_veg","background_bare","background_cropland","background_urban","background_other"))

p_land_use_background <- ggplot(d_land_use_background_percent_long, aes(x=ESU_season, y=proportion*100, fill=cat)) + 
  theme_bw(base_family = "Arial") +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color="black", size=14),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title = element_text(size=16)) +
  geom_col() +
  scale_fill_manual(values=colours, name="Land cover") +
  scale_x_discrete(labels=c("Breeding Western","Breeding Eastern","Non-breeding Western","Non-breeding Eastern")) +
  xlab("") + ylab("Proportion (%)")

p_land_use_background
ggsave("LK_land_use_background.pdf", p_land_use_background, device=cairo_pdf, units="cm", width=15, height=15, limitsize=FALSE)

## Heatmap of habitat preference f values for the most relevant habitats
d_hab_pref <- read.csv("habpref_R.csv", header=T, stringsAsFactors = T)
d_hab_pref$ESU_season <- factor(d_hab_pref$ESU_season,levels=c("Western_breeding","Eastern_breeding","Western_nonbreeding","Eastern_nonbreeding"))
d_hab_pref$habitat <- factor(d_hab_pref$habitat,levels=c("other","urban","cropland","bare","herb_veg","shrubs","forest"))

p_hab_pref <- ggplot(d_hab_pref, aes(x=ESU_season, y=habitat, fill=f)) +
  theme_bw(base_family = "Arial") +
  theme(panel.grid = element_blank(),
          axis.text = element_text(color="black", size=14),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.title = element_text(size=16)) +
  geom_tile() +
# scale_fill_viridis_c(option = "turbo") +
  scale_fill_distiller(palette="RdBu", limits=c(0,1)) +
  xlab("") + ylab("")

p_hab_pref
ggsave("LK_land_use_selection_test.pdf", p_hab_pref, device=cairo_pdf, units="cm", width=15, height=10, limitsize=FALSE)
