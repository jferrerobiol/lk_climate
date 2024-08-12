#!/bin/bash
#SBATCH -o /home/users/joan.ferrer/GRILLAIO/LK_genomics/MSMC2/scripts/logs/6_run_shapeit4.out
#SBATCH -J 6_run_shapeit4
#SBATCH --get-user-env
#SBATCH -p light
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mail-type=end
#SBATCH --mail-user=joan.ferrer@unimi.it
#SBATCH --account=grillaio
#SBATCH --time=2-00:00:00
#SBATCH --array=0-21

## Set your master path
MASTER=/home/users/joan.ferrer/GRILLAIO/LK_genomics/MSMC2/bams
REF=/home/users/joan.ferrer/GRILLAIO/genomic_data/bFalNau1/reference

## Read chromosomes
chroms=($(cat $REF/*fai | cut -f1 | grep -vE 'unloc|cafol|_Z|_W|arrow' | sort -V | uniq))

## Run shapeit using phase-informative read information and using a constant recombination rate of 1cM/Mb
shapeit4 --input $MASTER/merged_whatshap_phased.vcf.gz --region ${chroms[(($SLURM_ARRAY_TASK_ID))]} --use-PS 0.0001 --thread 4 --output $MASTER/${chroms[(($SLURM_ARRAY_TASK_ID))]}\_whatshap_shapeit4_phased.vcf.gz
tabix $MASTER/${chroms[(($SLURM_ARRAY_TASK_ID))]}\_whatshap_shapeit4_phased.vcf.gz
