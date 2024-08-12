#!/bin/bash
#SBATCH -o /home/users/joan.ferrer/GRILLAIO/LK_genomics/MSMC2/scripts/logs/5_merge_whatshap_individual_vcfs.out
#SBATCH -J 5_merge_whatshap_individual_vcfs
#SBATCH --get-user-env
#SBATCH -p light
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mail-type=end
#SBATCH --mail-user=joan.ferrer@unimi.it
#SBATCH --account=grillaio
#SBATCH --time=0-03:00:00

## Set your master path
MASTER=/home/users/joan.ferrer/GRILLAIO/LK_genomics/MSMC2/bams

## Read individuals
VCFs=($(ls $MASTER/*whatshap_phased.vcf.gz))

## Merge individual vcfs using bcftools merge
bcftools merge -Oz -o $MASTER/merged_whatshap_phased.vcf.gz --threads 12 -0 ${VCFs[@]}
tabix $MASTER/merged_whatshap_phased.vcf.gz
