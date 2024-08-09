##################################################
## Scripts to produce figures LK genomics paper ##
##################################################

library(tidyverse)
library(scales)
library(cowplot)
library(raster)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(fossil)
library(gtable)
library(gridExtra)
library(readxl)
library(scales)
library(ggpubr)
library(ggnewscale)
library(ragg)

setwd("~/Dropbox/Postdoc_Milan/LK_Joan/Figures/")

# Plot present SDM
ras_east_breed <- stack("~/Dropbox/Postdoc_Milan/LK_Joan/SDMs/all_SDMs/stack_breeding_east_crop_recl_2071_2100.grd")
ras_east_wint <- stack("~/Dropbox/Postdoc_Milan/LK_Joan/SDMs/all_SDMs/stack_wintering_east_crop_recl_2071_2100.grd")
ras_west_breed <- stack("~/Dropbox/Postdoc_Milan/LK_Joan/SDMs/all_SDMs/stack_breeding_west_crop_recl_2071_2100.grd")
ras_west_wint <- stack("~/Dropbox/Postdoc_Milan/LK_Joan/SDMs/all_SDMs/stack_wintering_west_crop_recl_2071_2100.grd")

admin <- ne_countries(scale = "medium", returnclass = "sf")

#### Vectorization of the climatic rasters for ggplot
ras_east_breed <- as(ras_east_breed, "SpatialPixelsDataFrame")
ras_east_breed <- as.data.frame(ras_east_breed)

ras_east_wint <- as(ras_east_wint, "SpatialPixelsDataFrame")
ras_east_wint <- as.data.frame(ras_east_wint)

ras_west_breed <- as(ras_west_breed, "SpatialPixelsDataFrame")
ras_west_breed <- as.data.frame(ras_west_breed)

ras_west_wint <- as(ras_west_wint, "SpatialPixelsDataFrame")
ras_west_wint <- as.data.frame(ras_west_wint)

p_current_SDM <- ggplot() + 
  geom_sf(data = admin, fill="black", lwd=0, colour="black") +
  geom_tile(data = subset(ras_east_breed, breeding_prediction==1), aes(x = x, y = y, fill = breeding_prediction), fill = "#5c7ec0") + 
  geom_tile(data = subset(ras_east_wint, wintering_prediction==1), aes(x = x, y = y, fill = wintering_prediction), fill = "#BDCBE5") + 
  geom_tile(data = subset(ras_west_breed, breeding_prediction==1), aes(x = x, y = y, fill = breeding_prediction), fill = "#ff6d00") + 
  geom_tile(data = subset(ras_west_wint, wintering_prediction==1), aes(x = x, y = y, fill = wintering_prediction), fill = "#ffc499") +
  coord_sf(xlim = c(150, -20), ylim = c(-38, 65), expand = F) +
  theme_void()
p_current_SDM

setwd("~/Dropbox/Postdoc_Milan/LK_Joan/Figures/")
ggsave("current_SDM.pdf", p_current_SDM, device="pdf", units="cm", width=5, height=3.5, limitsize=FALSE)

# Past SDMs
# LIG
p_LIG_SDM <- ggplot() + 
  geom_sf(data = admin, fill="black", lwd=0, colour="black") +
  geom_tile(data = subset(ras_east_breed, LIG_breeding_est_prediction==1), aes(x = x, y = y, fill = LIG_breeding_est_prediction), fill = "#5c7ec0") + 
  geom_tile(data = subset(ras_east_wint, LIG_wintering_est_prediction==1), aes(x = x, y = y, fill = LIG_wintering_est_prediction), fill = "#BDCBE5") + 
  geom_tile(data = subset(ras_west_breed, LIG_breeding_ovest_prediction==1), aes(x = x, y = y, fill = LIG_breeding_ovest_prediction), fill = "#ff6d00") + 
  geom_tile(data = subset(ras_west_wint, LIG_wintering_west_prediction==1), aes(x = x, y = y, fill = LIG_wintering_west_prediction), fill = "#ffc499") +
  coord_sf(xlim = c(150, -20), ylim = c(-38, 65), expand = F) +
  theme_void()
p_LIG_SDM

setwd("~/Dropbox/Postdoc_Milan/LK_Joan/Figures/")
ggsave("LIG_SDM.pdf", p_LIG_SDM, device="pdf", units="cm", width=50, height=35, limitsize=FALSE)

# Load stacks from CHELSA past high temporal resolution
ras_east_breed_past <- stack("~/Dropbox/Postdoc_Milan/LK_Joan/SDMs/CHELSA_past_models/stack_breeding_east_chelsa_past_recl.grd")
ras_west_breed_past <- stack("~/Dropbox/Postdoc_Milan/LK_Joan/SDMs/CHELSA_past_models/stack_breeding_west_chelsa_past_recl.grd")
ras_east_wint_past <- stack("~/Dropbox/Postdoc_Milan/LK_Joan/SDMs/CHELSA_past_models/stack_wintering_east_chelsa_past_recl.grd")
ras_west_wint_past <- stack("~/Dropbox/Postdoc_Milan/LK_Joan/SDMs/CHELSA_past_models/stack_wintering_west_chelsa_past_recl.grd")

# Vectorization of the climatic rasters for ggplot
ras_east_breed_past <- as(ras_east_breed_past, "SpatialPixelsDataFrame")
ras_east_breed_past <- as.data.frame(ras_east_breed_past)

ras_west_breed_past <- as(ras_west_breed_past, "SpatialPixelsDataFrame")
ras_west_breed_past <- as.data.frame(ras_west_breed_past)

ras_east_wint_past <- as(ras_east_wint_past, "SpatialPixelsDataFrame")
ras_east_wint_past <- as.data.frame(ras_east_wint_past)

ras_west_wint_past <- as(ras_west_wint_past, "SpatialPixelsDataFrame")
ras_west_wint_past <- as.data.frame(ras_west_wint_past)

# -20 kya
p_minustwenty_SDM <- ggplot() + 
  geom_sf(data = admin, fill="black", lwd=0, colour="black") +
  geom_tile(data = subset(ras_east_breed_past, p200_breeding_east_prediction_crop==1), aes(x = x, y = y, fill = MH_breeding_est_prediction), fill = "#5c7ec0") + 
  geom_tile(data = subset(ras_east_wint_past, p200_wintering_east_prediction_crop==1), aes(x = x, y = y, fill = MH_wintering_est_prediction), fill = "#BDCBE5") + 
  geom_tile(data = subset(ras_west_breed_past, p200_breeding_west_prediction_crop==1), aes(x = x, y = y, fill = MH_breeding_ovest_prediction), fill = "#ff6d00") + 
  geom_tile(data = subset(ras_west_wint_past, p200_wintering_west_prediction_crop==1), aes(x = x, y = y, fill = MH_wintering_west_prediction), fill = "#ffc499") +
  coord_sf(xlim = c(150, -20), ylim = c(-38, 65), expand = F) +
  theme_void()
p_minustwenty_SDM

setwd("~/Dropbox/Postdoc_Milan/LK_Joan/Figures/")
ggsave("-20_SDM.pdf", p_minustwenty_SDM, device="pdf", units="cm", width=50, height=35, limitsize=FALSE)

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

setwd("~/Dropbox/Postdoc_Milan/LK_Joan/Figures/")
ggsave("-18_SDM.pdf", p_minuseighteen_SDM, device="pdf", units="cm", width=50, height=35, limitsize=FALSE)

# -16 kya
p_minussixteen_SDM <- ggplot() + 
  geom_sf(data = admin, fill="black", lwd=0, colour="black") +
  geom_tile(data = subset(ras_east_breed_past, p160_breeding_east_prediction_crop==1), aes(x = x, y = y, fill = MH_breeding_est_prediction), fill = "#5c7ec0") + 
  geom_tile(data = subset(ras_east_wint_past, p160_wintering_east_prediction_crop==1), aes(x = x, y = y, fill = MH_wintering_est_prediction), fill = "#BDCBE5") + 
  geom_tile(data = subset(ras_west_breed_past, p160_breeding_west_prediction_crop==1), aes(x = x, y = y, fill = MH_breeding_ovest_prediction), fill = "#ff6d00") + 
  geom_tile(data = subset(ras_west_wint_past, p160_wintering_west_prediction_crop==1), aes(x = x, y = y, fill = MH_wintering_west_prediction), fill = "#ffc499") +
  coord_sf(xlim = c(150, -20), ylim = c(-38, 65), expand = F) +
  theme_void()
p_minussixteen_SDM

setwd("~/Dropbox/Postdoc_Milan/LK_Joan/Figures/")
ggsave("-16_SDM.pdf", p_minussixteen_SDM, device="pdf", units="cm", width=50, height=35, limitsize=FALSE)

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

setwd("~/Dropbox/Postdoc_Milan/LK_Joan/Figures/")
ggsave("-14_SDM.pdf", p_minusfourteen_SDM, device="pdf", units="cm", width=50, height=35, limitsize=FALSE)

# -12 kya
p_minustwelve_SDM <- ggplot() + 
  geom_sf(data = admin, fill="black", lwd=0, colour="black") +
  geom_tile(data = subset(ras_east_breed_past, p120_breeding_east_prediction_crop==1), aes(x = x, y = y, fill = MH_breeding_est_prediction), fill = "#5c7ec0") + 
  geom_tile(data = subset(ras_east_wint_past, p120_wintering_east_prediction_crop==1), aes(x = x, y = y, fill = MH_wintering_est_prediction), fill = "#BDCBE5") + 
  geom_tile(data = subset(ras_west_breed_past, p120_breeding_west_prediction_crop==1), aes(x = x, y = y, fill = MH_breeding_ovest_prediction), fill = "#ff6d00") + 
  geom_tile(data = subset(ras_west_wint_past, p120_wintering_west_prediction_crop==1), aes(x = x, y = y, fill = MH_wintering_west_prediction), fill = "#ffc499") +
  coord_sf(xlim = c(150, -20), ylim = c(-38, 65), expand = F) +
  theme_void()
p_minustwelve_SDM

setwd("~/Dropbox/Postdoc_Milan/LK_Joan/Figures/")
ggsave("-12_SDM.pdf", p_minustwelve_SDM, device="pdf", units="cm", width=50, height=35, limitsize=FALSE)

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

setwd("~/Dropbox/Postdoc_Milan/LK_Joan/Figures/")
ggsave("-10_SDM.pdf", p_minusten_SDM, device="pdf", units="cm", width=50, height=35, limitsize=FALSE)

# -8 kya
p_minuseight_SDM <- ggplot() + 
  geom_sf(data = admin, fill="black", lwd=0, colour="black") +
  geom_tile(data = subset(ras_east_breed_past, p80_breeding_east_prediction_crop==1), aes(x = x, y = y, fill = MH_breeding_est_prediction), fill = "#5c7ec0") + 
  geom_tile(data = subset(ras_east_wint_past, p80_wintering_east_prediction_crop==1), aes(x = x, y = y, fill = MH_wintering_est_prediction), fill = "#BDCBE5") + 
  geom_tile(data = subset(ras_west_breed_past, p80_breeding_west_prediction_crop==1), aes(x = x, y = y, fill = MH_breeding_ovest_prediction), fill = "#ff6d00") + 
  geom_tile(data = subset(ras_west_wint_past, p80_wintering_west_prediction_crop==1), aes(x = x, y = y, fill = MH_wintering_west_prediction), fill = "#ffc499") +
  coord_sf(xlim = c(150, -20), ylim = c(-38, 65), expand = F) +
  theme_void()
p_minuseight_SDM

setwd("~/Dropbox/Postdoc_Milan/LK_Joan/Figures/")
ggsave("-8_SDM.pdf", p_minuseight_SDM, device="pdf", units="cm", width=50, height=35, limitsize=FALSE)

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

setwd("~/Dropbox/Postdoc_Milan/LK_Joan/Figures/")
ggsave("-6_SDM.pdf", p_minussix_SDM, device="pdf", units="cm", width=50, height=35, limitsize=FALSE)

# -4 kya
p_minusfour_SDM <- ggplot() + 
  geom_sf(data = admin, fill="black", lwd=0, colour="black") +
  geom_tile(data = subset(ras_east_breed_past, p40_breeding_east_prediction_crop==1), aes(x = x, y = y, fill = MH_breeding_est_prediction), fill = "#5c7ec0") + 
  geom_tile(data = subset(ras_east_wint_past, p40_wintering_east_prediction_crop==1), aes(x = x, y = y, fill = MH_wintering_est_prediction), fill = "#BDCBE5") + 
  geom_tile(data = subset(ras_west_breed_past, p40_breeding_west_prediction_crop==1), aes(x = x, y = y, fill = MH_breeding_ovest_prediction), fill = "#ff6d00") + 
  geom_tile(data = subset(ras_west_wint_past, p40_wintering_west_prediction_crop==1), aes(x = x, y = y, fill = MH_wintering_west_prediction), fill = "#ffc499") +
  coord_sf(xlim = c(150, -20), ylim = c(-38, 65), expand = F) +
  theme_void()
p_minusfour_SDM

setwd("~/Dropbox/Postdoc_Milan/LK_Joan/Figures/")
ggsave("-4_SDM.pdf", p_minusfour_SDM, device="pdf", units="cm", width=50, height=35, limitsize=FALSE)

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

setwd("~/Dropbox/Postdoc_Milan/LK_Joan/Figures/")
ggsave("-2_SDM.pdf", p_minustwo_SDM, device="pdf", units="cm", width=50, height=35, limitsize=FALSE)

# 0 kya
p_zero_SDM <- ggplot() + 
  geom_sf(data = admin, fill="black", lwd=0, colour="black") +
  geom_tile(data = subset(ras_east_breed_past, p0_breeding_east_prediction_crop==1), aes(x = x, y = y, fill = MH_breeding_est_prediction), fill = "#5c7ec0") + 
  geom_tile(data = subset(ras_east_wint_past, p0_wintering_east_prediction_crop==1), aes(x = x, y = y, fill = MH_wintering_est_prediction), fill = "#BDCBE5") + 
  geom_tile(data = subset(ras_west_breed_past, p0_breeding_west_prediction_crop==1), aes(x = x, y = y, fill = MH_breeding_ovest_prediction), fill = "#ff6d00") + 
  geom_tile(data = subset(ras_west_wint_past, p0_wintering_west_prediction_crop==1), aes(x = x, y = y, fill = MH_wintering_west_prediction), fill = "#ffc499") +
  coord_sf(xlim = c(150, -20), ylim = c(-38, 65), expand = F) +
  theme_void()
p_zero_SDM

setwd("~/Dropbox/Postdoc_Milan/LK_Joan/Figures/")
ggsave("0_SDM.pdf", p_zero_SDM, device="pdf", units="cm", width=50, height=35, limitsize=FALSE)

# Plot SDMs showing distance between centroids of breeding distribution between 2 clusters
# Calculate centroids (to do that, load rasters again, then convert them as data frames for plotting)
centroid_current_east <- colMeans(xyFromCell(ras_east_breed_past, which(ras_east_breed_past$p0_breeding_east_prediction_crop[]==1)))
centroid_2040_east <- colMeans(xyFromCell(ras_east_breed, which(ras_east_breed$ukesm1_585_breeding_east_prediction[]==1)))
centroid_2070_east <- colMeans(xyFromCell(ras_east_breed, which(ras_east_breed$ukesm1_585_breeding_east_prediction_2071_2100[]==1)))

centroid_current_west <- colMeans(xyFromCell(ras_west_breed_past, which(ras_west_breed_past$p0_breeding_west_prediction_crop[]==1)))
centroid_2040_west <- colMeans(xyFromCell(ras_west_breed, which(ras_west_breed$ukesm1_585_breeding_west_prediction[]==1)))
centroid_2070_west <- colMeans(xyFromCell(ras_west_breed, which(ras_west_breed$ukesm1_585_breeding_west_prediction_2071_2100[]==1)))

# Calculate distance between east and west centroids
dist_current <- as.numeric(earth.dist(rbind(t(as.data.frame(centroid_current_east)),t(as.data.frame(centroid_current_west))))) #4952 km
dist_2040 <- as.numeric(earth.dist(rbind(t(as.data.frame(centroid_2040_east)),t(as.data.frame(centroid_2040_west))))) #5508 km
dist_2070 <- as.numeric(earth.dist(rbind(t(as.data.frame(centroid_2070_east)),t(as.data.frame(centroid_2070_west))))) #5252 km

# Plot SDMs for current, 2040 and 2070 with distance between centroids
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

setwd("~/Dropbox/Postdoc_Milan/LK_Joan/Figures/")
ggsave("future_SDMs_dist_centroids.pdf", dist_centroids, device="pdf", units="cm", width=15, height=7.5, limitsize=FALSE)

## Figure S1: presence data together with the SDM predicted ranges to show how well SDMs represent current known range
occurrence <- read.csv("../SDMs/occurrence_data/Table_S5_occurrence_data.csv")
occurrence$ESU_season <- paste(occurrence$ESU, occurrence$season, sep="_")

p_occurrence_sdm <- ggplot() + 
  geom_sf(data = admin, fill="#9e9a95", lwd=0, colour="#9e9a95") +
  geom_tile(data = subset(ras_east_breed, breeding_prediction==1), aes(x = x, y = y, fill = breeding_prediction), fill = "#5c7ec0") + 
  geom_tile(data = subset(ras_east_wint, wintering_prediction==1), aes(x = x, y = y, fill = wintering_prediction), fill = "#BDCBE5") + 
  geom_tile(data = subset(ras_west_breed, breeding_prediction==1), aes(x = x, y = y, fill = breeding_prediction), fill = "#ff6d00") + 
  geom_tile(data = subset(ras_west_wint, wintering_prediction==1), aes(x = x, y = y, fill = wintering_prediction), fill = "#ffc499") +
  geom_point(shape=21, data = occurrence, aes(x = longitude, y = latitude, fill = ESU_season), colour="black") +
  scale_fill_manual(values=c("#5c7ec0","#BDCBE5","#ff6d00","#ffc499"), guide="none") +
  coord_sf(xlim = c(150, -20), ylim = c(-38, 65), expand = F) +
  theme_void()
p_occurrence_sdm

ggsave("occurrence_current_SDM.pdf", p_occurrence_sdm, device="pdf", units="cm", width=20, height=14, limitsize=FALSE)

## Figure 1: sampling sites and occurrence data
range <- sf::st_read("~/Dropbox/Postdoc_Milan/LK_Joan/presence_data/Lesser_Kestrel/species_22696357.shp") %>% # Read shapefile of Birdlife distribution
  filter(legend == "Extant (non breeding)" | legend == "Extant (breeding)") %>% # Filter breeding and non-breeding distributions
  mutate_if(is.character, as.factor) # Convert character variables to factors

p_occurrence <- ggplot() + 
  geom_sf(data = range, aes(fill=legend), lwd=0, alpha=1) +
  geom_sf(data = admin, fill=NA, lwd=0.25, colour=gray(.5)) +
  scale_fill_manual(values=c("#C5AA83","#859EC6"), guide="none") +
  scale_color_manual(values=c("#C5AA83","#859EC6"), guide="none") +
  new_scale_fill() +
  geom_point(shape=21, data = occurrence, aes(x = longitude, y = latitude, fill = ESU_season), colour="black") +
  scale_fill_manual(values=c("#5c7ec0","#BDCBE5","#ff6d00","#ffc499"), guide="none") +
  coord_sf(xlim = c(150, -20), ylim = c(-38, 65), expand = F) +
  theme_void() +
  theme(panel.border = element_rect(color = "black", fill = NA, size=1))
p_occurrence

sampling <- read.csv("~/Google Drive/LK_manuscript/Figures_tables/tables/Table_S1_summary_sampling - Sheet1.csv")
sampling$Acronym <- factor(sampling$Acronym,levels = c("POR","ESN","ESS","SIC","ITS","ITN","CRO","GRG","GRC","GRL","TUR","ISR","KAZ","RUS","MON","MOS"))
## Add population order
pop_order <- c("POR","ESN","ESS","SIC","ITS","ITN","CRO","GRG","GRC","GRL","TUR","ISR","KAZ","RUS","MON","MOS")
## Add colours
palette<-(c("#a80000","#db0000","#dc2906","#dd4a11","#e16d19","#e68f21","#edb42a","#f5d933","#f5d933",
            "#e4e940","#cad34f","#b0bc66","#98a881","#7f929e","#677dbc","#5068dc"))
names(palette)<-pop_order

p_sampling <- ggplot() + 
  geom_sf(data = range, aes(fill=legend), lwd=0, alpha=1) +
  geom_sf(data = admin, fill=NA, lwd=0.25, colour=gray(.5)) +
  scale_fill_manual(values=c("#C5AA83","#859EC6"), guide="none") +
  scale_color_manual(values=c("#C5AA83","#859EC6"), guide="none") +
  new_scale_fill() +
  geom_point(shape=21, data = sampling, aes(x = Longitude, y = Latitude, fill = Acronym), colour="black", size=3.5) +
  scale_fill_manual(values=palette, guide="none") +
  coord_sf(xlim = c(150, -20), ylim = c(25, 60), expand = F) +
  theme_void() +
  theme(panel.border = element_rect(color = "black", fill = NA, size=1))
p_sampling

p_fig1 <- ggarrange(p_sampling, p_occurrence, ncol=1, align = "v")
ggsave("sampling_sites_occurrence_data.pdf", p_fig1, device="pdf", units="cm", width=20, height=20, limitsize=FALSE)
