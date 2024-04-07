#!/usr/bin/bash

#SBATCH --job-name Haplo_WF
#SBATCH --output %j_report_wf_PIPK.log
#SBATCH --partition batch
#SBATCH --cpus-per-task 1
#SBATCH --mem 8GB
#SBATCH --time 10-00:00:00

export MAMBA_ROOT_PREFIX=/shared/Software/micromamba
eval "$(micromamba shell hook --shell bash)"
micromamba activate snakemake

#snakemake --unlock
#snakemake -np
#snakemake --jobs 30 --profile slurm/
#snakemake --jobs 50  --slurm  --cluster-config slurm/config.yaml
snakemake --jobs 10  --slurm  --cluster-config slurm/config.yaml --snakefile make_report.smk  
#snakemake --jobs 30  --slurm  --cluster-config slurm/config.yaml  --until   

micromamba deactivate
#sbatch run_pipeline.sh
