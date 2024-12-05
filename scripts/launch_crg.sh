#!/usr/bin/env bash
#SBATCH -c 1
#SBATCH --mem=2G
#SBATCH --qos=pipelines
set -e
set -u

unset SLURM_JOB_ID
snakemake --profile profiles/crg "${@}"
