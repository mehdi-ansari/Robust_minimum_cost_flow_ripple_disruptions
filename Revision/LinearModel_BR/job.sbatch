#!/bin/bash
#SBATCH -p batch
#SBATCH -t 6:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mail-user=meansar@okstate.edu
#SBATCH --mail-type=end

module load cplex/12.10

./globalplanner