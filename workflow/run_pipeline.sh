#!/usr/bin/bash

#SBATCH --job-name HAPLO_WF
#SBATCH --output %j_HAPLO_WF.log
#SBATCH --partition batch
#SBATCH --cpus-per-task 1
#SBATCH --mem 8GB
#SBATCH --time 10-00:00:00

export MAMBA_ROOT_PREFIX=/shared/Software/micromamba
eval "$(micromamba shell hook --shell bash)"
micromamba activate snakemake

#snakemake --jobs 30 --profile slurm/   #--unlock
snakemake --jobs 30  --slurm  --cluster-config slurm/config.yaml
1
#micromamba deactivate

#sbatch run_pipeline.sh
