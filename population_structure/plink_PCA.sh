## Perform linkage pruning
plink --vcf intersect_Tasos_Joan.nomono.nosexchr_ncbi_annotated.vcf --double-id --allow-extra-chr --indep-pairwise 50 10 0.1 --out LK_intersect_norelated_wCro_nosexchrom

## Perform PCA
plink --vcf intersect_Tasos_Joan.nomono.nosexchr_ncbi_annotated.vcf --double-id --allow-extra-chr --extract LK_intersect_norelated_wCro_nosexchrom.prune.in --make-bed --pca --out LK_intersect_norelated_wCro_nosexchrom
