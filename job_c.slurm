#!/bin/bash
#SBATCH -A ctbp-levine
#SBATCH --job-name=smSimplex
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --mail-user=yz47@rice.edu
#SBATCH --mail-type=END
#SBATCH --export=ALL

ulimit -c unlimited
./cpu
