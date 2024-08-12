#!/bin/bash
#SBATCH -o /home/users/joan.ferrer/GRILLAIO/LK_genomics/MSMC2/scripts/logs/7_split_vcfs_per_individual.out
#SBATCH -J 7_split_vcfs_per_individual
#SBATCH --get-user-env
#SBATCH -p light
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mail-type=end
#SBATCH --mail-user=joan.ferrer@unimi.it
#SBATCH --account=grillaio
#SBATCH --time=2-00:00:00
#SBATCH --array=0-21%4

## Set your master path
MASTER=/home/users/joan.ferrer/GRILLAIO/LK_genomics/MSMC2/bams
REF=/home/users/joan.ferrer/GRILLAIO/genomic_data/bFalNau1/reference

## Read chromosomes
chroms=($(cat $REF/*fai | cut -f1 | grep -vE 'unloc|cafol|_Z|_W|arrow' | sort -V | uniq))

## Get individuals
inds=($(zcat $MASTER/${chroms[(($SLURM_ARRAY_TASK_ID))]}\_whatshap_shapeit4_phased.vcf.gz | grep -v "^##" | grep "^#" | cut -f10- | tr '\t' '\n'))

## Split vcfs per individual
for ind in ${inds[@]}; do
  bcftools view -Oz -s $ind -o $MASTER/$ind\_${chroms[(($SLURM_ARRAY_TASK_ID))]}\_whatshap_shapeit4_phased.vcf.gz $MASTER/${chroms[(($SLURM_ARRAY_TASK_ID))]}\_whatshap_shapeit4_phased.vcf.gz
done
