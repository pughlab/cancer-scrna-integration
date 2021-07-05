#!/bin/bash
#SBATCH -t 6:00:00
#SBATCH --mem=60G
#SBATCH -p himem
#SBATCH -c 1
#SBATCH -N 1

module load R/4.0.0

# where $1 is the Rscript
# ... $2 = kweight
# ... $3 = distpct
# ... $4 = numvargenes

Rscript $1 --kweight $2 --distpct $3 --vargenes $4
