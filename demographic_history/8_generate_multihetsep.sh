#!/bin/bash
#SBATCH -o /home/users/joan.ferrer/GRILLAIO/LK_genomics/MSMC2/scripts/logs/8_generate_multihetsep.out
#SBATCH -J 8_generate_multihetsep
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
OUTDIR=/home/users/joan.ferrer/GRILLAIO/LK_genomics/multihetsep
chroms=($(cat $REF/*fai | cut -f1 | grep -vE 'unloc|cafol|_Z|_W|arrow' | sort -V | uniq))
SOFT=/home/users/joan.ferrer/GRILLAIO/soft/msmc-tools

$SOFT/generate_multihetsep.py --chr ${chroms[(($SLURM_ARRAY_TASK_ID))]} \
    --mask $MASTER/maternal.merged.dups.bam.${chroms[(($SLURM_ARRAY_TASK_ID))]}.mask.bed.gz --mask $MASTER/paternal.merged.dups.bam.${chroms[(($SLURM_ARRAY_TASK_ID))]}.mask.bed.gz \
    --mask $MASTER/LK_83M.merged.dups.bam.${chroms[(($SLURM_ARRAY_TASK_ID))]}.mask.bed.gz --mask $MASTER/LK_F7.merged.dups.bam.${chroms[(($SLURM_ARRAY_TASK_ID))]}.mask.bed.gz \
    --mask $MASTER/LK_F8.merged.dups.bam.${chroms[(($SLURM_ARRAY_TASK_ID))]}.mask.bed.gz --mask $MASTER/LK_M2.merged.dups.bam.${chroms[(($SLURM_ARRAY_TASK_ID))]}.mask.bed.gz \
    --mask $MASTER/${chroms[(($SLURM_ARRAY_TASK_ID))]}\_mapp_mask.bed.gz \
    $MASTER/maternal_${chroms[(($SLURM_ARRAY_TASK_ID))]}_whatshap_shapeit4_phased.vcf.gz $MASTER/paternal_${chroms[(($SLURM_ARRAY_TASK_ID))]}_whatshap_shapeit4_phased.vcf.gz \
    $MASTER/83M_LK_Italy_${chroms[(($SLURM_ARRAY_TASK_ID))]}_whatshap_shapeit4_phased.vcf.gz $MASTER/F7_LK_Tuva_${chroms[(($SLURM_ARRAY_TASK_ID))]}_whatshap_shapeit4_phased.vcf.gz \
    $MASTER/F8_LK_Tuva_${chroms[(($SLURM_ARRAY_TASK_ID))]}_whatshap_shapeit4_phased.vcf.gz $MASTER/M2_LK_Tuva_${chroms[(($SLURM_ARRAY_TASK_ID))]}_whatshap_shapeit4_phased.vcf.gz \
    > $OUTDIR/Eur_Asi_${chroms[(($SLURM_ARRAY_TASK_ID))]}.multihetsep.txt
