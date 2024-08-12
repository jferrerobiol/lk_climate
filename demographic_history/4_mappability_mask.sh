#!/bin/bash
#SBATCH -o /home/users/joan.ferrer/GRILLAIO/LK_genomics/MSMC2/scripts/logs/4_mappability_mask.out
#SBATCH -J 4_mappability_mask
#SBATCH --get-user-env
#SBATCH -p light
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mail-type=end
#SBATCH --mail-user=joan.ferrer@unimi.it
#SBATCH --account=grillaio
#SBATCH --time=0-05:00:00

MASTER=/home/users/joan.ferrer/GRILLAIO/LK_genomics/MSMC2
REF=/home/users/joan.ferrer/GRILLAIO/genomic_data/bFalNau1/reference

## Run genmap with kmer size = 150 and allowing 2 mismatches and on 30 threads ##
genmap index -F $REF/bFalNau1.pri.cur.20200818.fasta -I $MASTER/map_mask
genmap map -K 150 -E 2 -T 30 -I $MASTER/map_mask -O $MASTER/map_mask -t -w -bg

## Generate a bed file containing only areas with mappability > 0.5 ##
cat $MASTER/map_mask/bFalNau1.pri.cur.20200818.genmap.bedgraph | awk '$4 >= 0.5 {print $0}' > $MASTER/map_mask/genmap_150_output_more0.5.bed

## Subtract repetitive areas from mappability mask
sort -k1,1 -k2,2n $MASTER/map_mask/genmap_150_output_more0.5.bed > $MASTER/map_mask/genmap_150_output_more0.5.sorted.bed
subtractBed -a $MASTER/map_mask/genmap_150_output_more0.5.sorted.bed -b $REF/masking/bFalNau1.pri.cur.20200818_REPMASK_REPEATS.bed > $MASTER/map_mask/mapp_mask_repeats_subtract.bed

## Calculate percentage of the genome covered by the mappability mask ##
cat $MASTER/map_mask/mapp_mask_repeats_subtract.bed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
# Divide this number by total length of assembly and multiply by 100 to obtain percentage

## Divide the mappability mask per scaffold ##
chroms=($(cat $REF/*fai | cut -f1 | grep -vE 'unloc|cafol|_Z|_W|arrow' | sort -V | uniq))
for chrom in ${chroms[@]}; do cat $MASTER/map_mask/mapp_mask_repeats_subtract.bed | grep $chrom | gzip > $MASTER/bams/$chrom\_mapp_mask.bed.gz; done
