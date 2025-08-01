#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mem-per-cpu=5G
#SBATCH --job-name=PDS_Lap_simulations
#SBATCH --time=5:00:00
#SBATCH --partition=scavenge
#SBATCH --requeue
#SBATCH --array=1-500

module purge
module load R/4.3
module load R-bundle-CRAN/2023.12-foss-2022b
export OMP_NUM_THREADS=1
module list

params=$(sed -n ${SLURM_ARRAY_TASK_ID}p ./repeat.txt)
Rscript PDS-Lap_simulations.R $params 