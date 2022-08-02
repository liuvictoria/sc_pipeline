#!/bin/bash
#SBATCH --job-name=panT_aggregation
#SBATCH --mail-type=END
#SBATCH --mail-user=victoria.liu@jax.org
#SBATCH -p compute
#SBATCH -q batch
#SBATCH -t 48:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4

module load singularity

singularity run ~/seurat_latest.sif R -f panT_aggregation.R
