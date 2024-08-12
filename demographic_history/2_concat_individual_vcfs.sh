#!/bin/bash
#SBATCH -o /home/users/joan.ferrer/GRILLAIO/LK_genomics/MSMC2/scripts/logs/2_concat_individual_vcfs
#SBATCH -J 2_concat_individual_vcfs
#SBATCH --get-user-env
#SBATCH -p light
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mail-type=end
#SBATCH --mail-user=joan.ferrer@unimi.it
#SBATCH --account=grillaio
#SBATCH --time=2-00:00:00
#SBATCH --array=0-5

## Set your master path
MASTER=/home/users/joan.ferrer/GRILLAIO/LK_genomics/MSMC2/bams
REF=/home/users/joan.ferrer/GRILLAIO/genomic_data/bFalNau1/reference

## Read chromosomes
chroms=($(cat $REF/*fai | cut -f1 | grep -vE 'unloc|cafol|_Z|_W|arrow' | sort -V | uniq))

## Read individuals
IND=($(ls $MASTER/*bam | sed 's/.merged.dups.bam//g'))

## List vcf files
vcfs=($(for chr in ${chroms[@]}; do echo ${IND[(($SLURM_ARRAY_TASK_ID))]}\.merged.dups.bam.$chr\.vcf.gz; done))

## Concatenate vcf files per individual
bcftools concat -Oz -o ${IND[(($SLURM_ARRAY_TASK_ID))]}\.concat.vcf.gz ${vcfs[@]}

## Sort and index vcf files
bcftools sort -Oz -o ${IND[(($SLURM_ARRAY_TASK_ID))]}\.concat.sorted.vcf.gz ${IND[(($SLURM_ARRAY_TASK_ID))]}\.concat.vcf.gz
tabix ${IND[(($SLURM_ARRAY_TASK_ID))]}\.concat.vcf.gz
