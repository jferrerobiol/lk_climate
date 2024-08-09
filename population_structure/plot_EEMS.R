library(rEEMSplots)
library(rgdal)
library(rworldmap)
library(rworldxtra)

#example
eems_results <- file.path("~/Dropbox/Postdoc_Milan/LK_Joan/EEMS/stacks_all/chain1/")

eems.plots(mcmcpath=eems_results, plotpath=paste0(eems_results,"final"),
           projection.in = "+proj=longlat +datum=WGS84",
           projection.out = "+proj=merc +datum=WGS84",
           longlat=T, out.png = F, add.grid = F, add.outline = T, lwd.outline = 1, add.demes = T,
           add.map=T, lwd.map=0.5,min.cex.demes = 0.5, max.cex.demes = 1.5)

#new polygon
eems_polygon <- read.delim("~/Dropbox/Postdoc_Milan/LK_Joan/EEMS/input/LK_minconvexpoly.txt")
world <- ne_countries(scale = "medium", returnclass = "sf")
ggplot(data = world) +
  geom_sf() +
  geom_point(data = eems_polygon, aes(x = longitude, y = latitude), size = 4, 
             shape = 23, fill = "darkred") +
  coord_sf(xlim = c(-10, 130), ylim = c(25, 65), expand = FALSE)
