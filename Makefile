TARGETS=dependencies dag run unlock

CONDA_ENV_DIR=$(shell dirname ${CONDA_EXE})
HN=$(shell hostname | sed "s/[0-9]//g")


ifeq ($(HN),$(filter $(HN),cnode gnode hnode))
  CONDA_ENV_NAME=/exchange/healthds/software/envs/snakemake
else
  CONDA_ENV_NAME=snakemake
endif

all:
	@echo "Try one of: ${TARGETS}"

dag:
	source $(CONDA_ENV_DIR)/activate $(CONDA_ENV_NAME) && \
	snakemake --dag | dot -Tsvg > dag.svg

dependencies:
	mamba env update -n snakemake --file environment.yml

dev-dependencies: dependencies
	mamba env update -n snakemake --file environment_dev.yml

dry-run:
	source $(CONDA_ENV_DIR)/activate $(CONDA_ENV_NAME) && \
	snakemake --sdm conda --dry-run --profile slurm --snakefile workflow/Snakefile

pre-commit:
	source $(CONDA_ENV_DIR)/activate $(CONDA_ENV_NAME) && \
	if [ ! -f .git/hooks/pre-commit ]; then pre-commit install; fi && \
	pre-commit run --all-files

run:
	source $(CONDA_ENV_DIR)/activate $(CONDA_ENV_NAME) && \
	snakemake --profile slurm --snakefile workflow/Snakefile

rerun:
	source $(CONDA_ENV_DIR)/activate $(CONDA_ENV_NAME) && \
	snakemake --profile slurm --snakefile workflow/Snakefile --rerun-incomplete

unlock:
	source $(CONDA_ENV_DIR)/activate $(CONDA_ENV_NAME) && \
	snakemake --unlock

dockerfile_:
	source $(CONDA_ENV_DIR)/activate $(CONDA_ENV_NAME) && \
	snakemake --containerize --snakefile workflow/Snakefile > Dockerfile

test:
	source $(CONDA_ENV_DIR)/activate $(CONDA_ENV_NAME) && \
	snakemake --profile slurm --snakefile workflow/test.smk

dry-test:
	source $(CONDA_ENV_DIR)/activate $(CONDA_ENV_NAME) && \
	snakemake --sdm conda --dry-run --profile slurm --snakefile workflow/test.smk
