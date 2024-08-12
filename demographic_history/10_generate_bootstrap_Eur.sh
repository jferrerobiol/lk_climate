#!/bin/bash
#SBATCH -o /home/users/joan.ferrer/GRILLAIO/LK_genomics/MSMC2/scripts/logs/11_generate_bootstrap_Eur.out
#SBATCH -J 11_generate_bootstrap_Eur
#SBATCH --get-user-env
#SBATCH -p light
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --mail-type=end
#SBATCH --mail-user=joan.ferrer@unimi.it
#SBATCH --account=grillaio
#SBATCH --time=0-02:00:00

## Set your master path
INDIR=/home/users/joan.ferrer/GRILLAIO/LK_genomics/MSMC2/multihetsep
SOFT=/home/users/joan.ferrer/GRILLAIO/soft/msmc-tools

python3 $SOFT/multihetsep_bootstrap.py -n 100 --chunks_per_chromosome 10 --nr_chromosomes 22 $INDIR/bootstrap_dir $INDIR/Eur_Asi_SUPER_*.multihetsep.txt
