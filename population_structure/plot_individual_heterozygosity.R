## First I calculate the proportion of polymorphic sites for each individual with the script:
# vcftools --vcf ../GEA/intersect_Tasos_Joan.nomono.nosexchr_ncbi_annotated.vcf --max-missing 1 --het --out intersect_Tasos_Joan.nomono.nosexchr_ncbi_annotated_nomissing

## Load libraries
library(tidyverse)
library(ggsci)
library(ggpubr)
library(scales)

## Load individual heterozygosities file
setwd("~/Dropbox/Postdoc_Milan/LK_Joan/ind_het/")
ind_het_100 <- read.delim("intersect_Tasos_Joan.nomono.nosexchr_ncbi_annotated_nomissing.het")
popmap <- read.delim("../intersect_Tasos_Joan_vcfs/popmap.txt")
ind_het_100 <- ind_het_100 %>% left_join(popmap, by=c("INDV"="ind"))
ind_het_100$pop <- factor(ind_het_100$pop, levels=c("ESN","ESS","SIC","ITS","ITN","CRO","GRG","GRC","GRL","TUR","ISR","KAZ","RUS","MON","MOS"))
mypal<-(c("#db0000","#dc2906","#dd4a11","#e16d19","#e68f21","#edb42a","#f5d933","#f5d933",
            "#e4e940","#cad34f","#b0bc66","#98a881","#7f929e","#677dbc","#5068dc"))

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
