#!/bin/bash
#SBATCH -o /home/users/joan.ferrer/GRILLAIO/LK_genomics/MSMC2/scripts/logs/9a_run_msmc2_time_patterning_Eur.out
#SBATCH -J 9a_run_msmc2_time_patterning_Eur
#SBATCH --get-user-env
#SBATCH -p bigmem
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=11
#SBATCH --mem=500G
#SBATCH --mail-type=end
#SBATCH --mail-user=joan.ferrer@unimi.it
#SBATCH --account=grillaio
#SBATCH --time=2-00:00:00

## Set your master path
INDIR=/home/users/joan.ferrer/GRILLAIO/LK_genomics/MSMC2/multihetsep
OUTDIR=/home/users/joan.ferrer/GRILLAIO/LK_genomics/MSMC2/msmc2_output
SOFT=/home/users/joan.ferrer/GRILLAIO/soft/msmc-tools

msmc2 -t 11 -p 1*4+30*2+1*4+1*6+1*10 -s -I 0,1,2,3,4,5 -o $OUTDIR/Eur_6hap $INDIR/Eur_Asi_*.multihetsep.txt
