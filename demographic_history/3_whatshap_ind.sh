#!/bin/bash
#SBATCH -o /home/users/joan.ferrer/GRILLAIO/LK_genomics/MSMC2/scripts/logs/3_whatshap_ind.out
#SBATCH -J 3_whatshap_ind
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

## Read individuals
IND=($(ls $MASTER/*bam | sed 's/.merged.dups.bam//g'))

## Phase
whatshap phase -o ${IND[(($SLURM_ARRAY_TASK_ID))]}\.whatshap_phased.vcf.gz \
  --reference=/home/users/joan.ferrer/GRILLAIO/genomic_data/bFalNau1/reference/bFalNau1.pri.cur.20200818.fasta \
  ${IND[(($SLURM_ARRAY_TASK_ID))]}\.concat.vcf.gz ${IND[(($SLURM_ARRAY_TASK_ID))]}\.merged.dups.bam

## Index vcf files
tabix ${IND[(($SLURM_ARRAY_TASK_ID))]}\.whatshap_phased.vcf.gz
