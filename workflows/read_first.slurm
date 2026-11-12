#!/usr/bin/env bash
#SBATCH --job-name=hgtlong
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=24:00:00
set -euo pipefail
: "${PROJECT:?Set PROJECT to an absolute project directory before sbatch}"
: "${SLURM_JOB_ID:?Submit this script through Slurm}"
cd "$PROJECT"
command -v hgtlong
command -v minimap2
command -v samtools
mkdir -p run
hgtlong build-panel --genomes input/MAGs --outdir run/panel
# Reserve two of the 16 CPUs for BAM writing and orchestration.
minimap2 -t 14 -ax map-ont -p 0 -N 100 -Y -K 100M   run/panel/panel.fa input/reads.fastq.gz   | samtools view -b -o run/mappings.bam -
samtools collate -Ou run/mappings.bam run/collate_tmp   | samtools view -h -   | hgtlong extract --sam - --panel run/panel/panel.tsv       --workers 8 --out run/evidence.tsv.gz
length_args=()
if [[ -n "${MIN_HGT_LENGTH:-}" ]]; then length_args=(--min-hgt-length "$MIN_HGT_LENGTH"); fi
hgtlong call "${length_args[@]}" --evidence run/evidence.tsv.gz --outdir run/calls
