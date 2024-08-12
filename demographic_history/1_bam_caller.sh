#!/bin/bash
#SBATCH -o /home/users/joan.ferrer/GRILLAIO/LK_genomics/MSMC2/scripts/logs/1_bam_caller.out
#SBATCH -J 1_bam_caller
#SBATCH --get-user-env
#SBATCH -p light
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mail-type=end
#SBATCH --mail-user=joan.ferrer@unimi.it
#SBATCH --account=grillaio
#SBATCH --time=2-00:00:00
#SBATCH --array=0-5

## Load modules
module load genetics/broadinstitute

## Set your master path
MASTER=/home/users/joan.ferrer/GRILLAIO/LK_genomics/MSMC2/bams
REF=/home/users/joan.ferrer/GRILLAIO/genomic_data/bFalNau1/reference
bamcaller=/home/users/joan.ferrer/GRILLAIO/soft/msmc-tools

## Read bams and chromosomes ##
bams=($(ls $MASTER/*bam))
chroms=($(cat $REF/*fai | cut -f1 | grep -vE 'unloc|cafol|_Z|_W|arrow' | sort -V))

## Estimate the average coverage from chromosome SUPER_20

depth=$(samtools depth -r SUPER_11 ${bams[(($SLURM_ARRAY_TASK_ID))]} | awk '{sum += $3} END {print sum / NR}')

## Call SNP and generate mask and vcf files

for chrom in ${chroms[@]}; do 
    samtools mpileup -B -q 20 -Q 20 -C 50 -g -r $chrom -f $REF/*fasta ${bams[(($SLURM_ARRAY_TASK_ID))]} | bcftools call -c -V indels |\
    $bamcaller/bamCaller.py $depth ${bams[(($SLURM_ARRAY_TASK_ID))]}\.$chrom.mask.bed.gz | gzip -c > ${bams[(($SLURM_ARRAY_TASK_ID))]}\.$chrom.vcf.gz
done
