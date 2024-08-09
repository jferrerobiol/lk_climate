### plot admixture results

## Script based on https://luisdva.github.io/rstats/model-cluster-plots/ 

## Load libs
lib<-c("ggplot2","gridExtra","grid","dplyr","stringi","forcats","ggrepel","purrr","cowplot","tidyr","plyr")
lapply(lib,library,character.only=T)

## setwd
setwd("~/Dropbox/Postdoc_Milan/LK_Joan/admixture/")

## read in individual to population map 
## Add populations
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
## Add population order
pop_order <- c("ESN","ESS","SIC","ITS","ITN","CRO","GRG","GRC","GRL","TUR","ISR","KAZ","RUS","MON","MOS")
inds$Populations <- factor(inds$Populations, levels=pop_order)

## read in K values
K2 <- read.table("LK_intersect_norelated_wCro_nosexchrom.2.Q",col.names=c("1","2"),check.names = F)
K3 <- read.table("LK_intersect_norelated_wCro_nosexchrom.3.Q",col.names=c("1","2","3"),check.names = F)
K4 <- read.table("LK_intersect_norelated_wCro_nosexchrom.4.Q",col.names=c("1","2","3","4"),check.names = F)

## add sample and population information and select K values
K2$sampleID <- inds$Sample
K2$loc <- inds$Populations
K2 <- K2 %>% gather(key=popGroup,value=prob,1:2) %>%
  dplyr::select(sampleID,popGroup,prob,loc)
K2$loc <- factor(K2$loc, levels=pop_order)

K3$sampleID <- inds$Sample
K3$loc <- inds$Populations
K3 <- K3 %>% gather(key=popGroup,value=prob,1:3) %>%
  dplyr::select(sampleID,popGroup,prob,loc)
K3$loc <- factor(K3$loc, levels=pop_order)

K4$sampleID <- inds$Sample
K4$loc <- inds$Populations
K4 <- K4 %>% gather(key=popGroup,value=prob,1:4) %>%
  dplyr::select(sampleID,popGroup,prob,loc)
K4$loc <- factor(K4$loc, levels=pop_order)

## set up for plotting
## set colour palettes

## use palette1 for K2
palette1 <- c("#ff6d00","#5c7ec0")
## use palette2 for K3
palette2 <- c("#5c7ec0","#40A4D0","#ff6d00")
## use palette3 for K4
palette3 <- c("#40A4D0","#D43F3A","#8CB04F","#9632B8")

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

K3_plot <- ggplot(data=K3, aes(factor(sampleID), prob, fill = factor(popGroup), colour = factor(popGroup))) +
  geom_col(size = 1) +
  ylab(expression(italic(K)~"= 3"))+
  facet_grid(~loc, switch = "x", scales = "free", space = "free") +
  theme_minimal() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  scale_colour_manual(values=palette2)+
  scale_fill_manual(values=palette2)+
  theme(panel.spacing.x = unit(0.1, "lines"),
        axis.text.x = element_blank(),
        axis.text.y=element_text(family="Arial",size=14),
        axis.title.x=element_blank(),
        axis.title.y=element_text(family="Arial",size=16),
        panel.grid = element_blank(),
        legend.position="none",
        strip.background = element_blank(),
        strip.text.x=element_text(family="Arial",angle=90,size=12))

K3_plot

K4_plot <- ggplot(data=K4, aes(factor(sampleID), prob, fill = factor(popGroup), colour = factor(popGroup))) +
  geom_col(size = 1) +
  ylab(expression(italic(K)~"= 4"))+
  facet_grid(~loc, switch = "x", scales = "free", space = "free") +
  theme_minimal() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  scale_colour_manual(values=palette3)+
  scale_fill_manual(values=palette3)+
  theme(panel.spacing.x = unit(0.1, "lines"),
        axis.text.x = element_blank(),
        axis.text.y=element_text(family="Arial",size=14),
        axis.title.x=element_blank(),
        axis.title.y=element_text(family="Arial",size=16),
        panel.grid = element_blank(),
        legend.position="none",
        strip.background = element_blank(),
        strip.text.x=element_text(family="Arial",angle=90,size=12))

K4_plot

## put plots together and save them
all_K <- plot_grid(K2_plot,K3_plot,K4_plot,ncol=1)
ggsave("admixture_K2-4.pdf", all_K, width = 20, height = 25, device=cairo_pdf, units="cm", limitsize = FALSE)

## load CV error table
cverror <- read.delim("CV_error_maxmiss1.tsv")

CV_error <- ddply(cverror, c("K"), summarise,
                  N    = length(CV_error),
                  mean = mean(CV_error),
                  sd   = sd(CV_error),
                  se   = sd / sqrt(N)
)

CV_error

CV_error$K <- factor(CV_error$K, levels = CV_error$K)

## plot CV error
cverror_plot <- ggplot(CV_error, aes(x=K, y=mean)) +
  theme_bw(base_family = "Arial") +
  theme(panel.grid = element_blank()) +
  theme(axis.text = element_text(color="black", size=14),
        axis.title = element_text(size=16)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.3) +
  geom_line(group=1) +
  geom_point(size=2)+
  labs(x=expression(italic(K)), y="CV error")

cverror_plot
ggsave("admixture_CVerror_maxmiss1.pdf", cverror_plot, width = 20, height = 10, device=cairo_pdf, units="cm", limitsize = FALSE)