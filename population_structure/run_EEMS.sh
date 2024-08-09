## Run plink
plink --vcf intersect_Tasos_Joan.nomono.nosexchr_ncbi_annotated.vcf --double-id --allow-extra-chr --make-bed --out LK

#Before using plink I modify the vcf file leaving only numbers in sample names, removing SUPER_Z and SUPER_W and removing SUPER_ from the name of chromosomes

## Run bed2diffs making sure you run it from the soft folder
./bed2diffs_v1 --bfile /Users/apple/Dropbox/Postdoc_Milan/LK_Joan/EEMS/stacks_all/LK

## Run runeems_snps
/Users/apple/soft/eems/runeems_snps/src/runeems_snps --params params-chain1.ini --seed 1234
