## Script to get elevation of all presence points and plot elevation of Asian vs. European ##
library(elevatr)
library(tidyverse)
library(gghalves)
library(ggdist)

### Load data for each of the groups
east_breed <- read.csv("~/Dropbox/Postdoc_Milan/LK_Joan/GEA/old_Mattia/eastern_data_bioclim_5km_buffer_no_dupl.csv")
east_breed <- east_breed %>% select(X1,Y1)
east_breed$group <- "eastern_breeding"
west_breed <- read.csv("~/Dropbox/Postdoc_Milan/LK_Joan/GEA/old_Mattia/western_data_bioclim_5km_buffer_no_dupl.csv")
west_breed <- west_breed %>% select(X1,Y1)
west_breed$group <- "western_breeding"
east_nonbreed <- read.csv("~/Dropbox/Postdoc_Milan/LK_Joan/GEA/old_Mattia/wintering_eastern_bioclim_5km_buffer_no_dupl.csv")
east_nonbreed <- east_nonbreed %>% select(X1,Y1)
east_nonbreed$group <- "eastern_nonbreeding"
west_nonbreed <- read.csv("~/Dropbox/Postdoc_Milan/LK_Joan/GEA/old_Mattia/wintering_western_bioclim_5km_buffer_no_dupl.csv")
west_nonbreed <- west_nonbreed %>% select(X1,Y1)
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
ggsave("~/Dropbox/Postdoc_Milan/LK_Joan/presence_data/LK_elevation_per_group.pdf", p_elevation, device=cairo_pdf, units="cm", width=24, height=12, limitsize=FALSE)

# Check preesence data
ggplot(data = presence) + 
  geom_sf(data = admin, fill=gray(.8), colour=gray(.8), size=0.05) +
  geom_point(aes(x=x, y=y,col=group)) + 
  coord_sf(xlim = c(130, -30), ylim = c(-45, 65), expand = F) +
  scale_x_continuous(breaks=seq(-20,150,20)) +
  xlab("Longitude") + ylab("Latitude") +
  guides(fill=guide_legend(title="Adaptive index")) +
  theme_bw(base_size = 11, base_family = "Arial") +
  theme(panel.grid = element_blank(), plot.background = element_blank(), panel.background = element_blank(), strip.text = element_text(size=11))

t.test(presence$elevation[presence$group=="Western breeding"],
       presence$elevation[presence$group=="Eastern breeding"]) #p-value < 2.2e-16
wilcox.test(presence$elevation[presence$group=="Western breeding" | presence$group=="Eastern breeding"] ~ presence$group[presence$group=="Western breeding" | presence$group=="Eastern breeding"]) #p-value < 2.2e-16
mean(presence$elevation[presence$group=="Western breeding"]) #455.08 m
mean(presence$elevation[presence$group=="Eastern breeding"]) #890.43 m

t.test(presence$elevation[presence$group=="Western non-breeding"],
       presence$elevation[presence$group=="Eastern non-breeding"]) #p-value < 2.2e-16
wilcox.test(presence$elevation[presence$group=="Western non-breeding" | presence$group=="Eastern non-breeding"] ~ presence$group[presence$group=="Western non-breeding" | presence$group=="Eastern non-breeding"]) #p-value < 2.2e-16
mean(presence$elevation[presence$group=="Western non-breeding"]) #454.61 m
mean(presence$elevation[presence$group=="Eastern non-breeding"]) #890.91 m
