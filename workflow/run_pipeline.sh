#!/usr/bin/bash

#SBATCH --job-name HAPLO_WF
#SBATCH --output %j_HAPLO_WF.log
#SBATCH --partition batch
#SBATCH --cpus-per-task 1
#SBATCH --mem 128GB
#SBATCH --time 8:00:00

export MAMBA_ROOT_PREFIX=/shared/Software/micromamba
eval "$(micromamba shell hook --shell bash)"
micromamba activate snakemake

#snakemake --jobs 8  --reason --until get_dosage --default-resource mem_gb=8192  --latency-wait 10  --keep-going  --cluster 'sbatch  -p fast -c 8 --mem-per-cpu=8GB'
snakemake --jobs 30 --profile slurm/   #--unlock

#micromamba deactivate

#sbatch run_pipeline.txt -p fast
