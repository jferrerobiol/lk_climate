library(vcfR)
library(ape)
library(adegenet)
library(phangorn)

setwd("~/Dropbox/Postdoc_Milan/LK_Joan/intersect_Tasos_Joan_vcfs/vcfs/")
vcf <- read.vcfR("intersect_Tasos_Joan.nomono.nosexchr_ncbi_annotated.vcf")
dnabin <- vcfR2DNAbin(vcf,extract.indels=T,consensus=T,extract.haps=F,unphased_as_NA=F)

dist <- dist.dna(dnabin, model = "TN93", variance = FALSE,
                 gamma = FALSE, pairwise.deletion = TRUE,
                 base.freq = NULL, as.matrix = TRUE)

setwd("~/Dropbox/Postdoc_Milan/LK_Joan/splitstree/")

write.nexus.dist(dist, file = "LK_intersect_dist.nex", append = FALSE, upper = FALSE,
                 diag = TRUE, digits = getOption("digits"))
