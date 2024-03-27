#!/usr/bin/bash

#SBATCH --job-name HAPLO_WF
#SBATCH --output %j_HAPLO_WF.log
#SBATCH --partition fast
#SBATCH --cpus-per-task 1
#SBATCH --mem-per-cpu 8GB


export MAMBA_ROOT_PREFIX=/shared/Software/micromamba
eval "$(micromamba shell hook --shell bash)"
micromamba activate snakemake

#snakemake --jobs 8  --reason --until get_dosage --default-resource mem_gb=8192  --latency-wait 10  --keep-going  --cluster 'sbatch  -p fast -c 8 --mem-per-cpu=8GB'
snakemake --jobs 1 --profile slurm/

conda deactivate

#sbatch run_pipeline.txt -p fast
