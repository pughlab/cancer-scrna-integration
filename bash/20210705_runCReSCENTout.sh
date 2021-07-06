#!/bin/bash
#SBATCH -t 15:00:00
#SBATCH --mem=30G
#SBATCH -p all
#SBATCH -c 10
#SBATCH -N 1

module load R/4.0.0

# where $1 is the seurat file

Rscript crescent_files.r --seurat $1 --inputPath /cluster/projects/pughlab/projects/cancer_scrna_integration/integration
