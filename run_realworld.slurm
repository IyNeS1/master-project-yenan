#!/bin/bash
#SBATCH --job-name=realworld_pipeline
#SBATCH --partition=sapphire
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=slurm_%j.out
#SBATCH --error=slurm_%j.err
#SBATCH --account=punim2974

module load Nextflow/25.10.4

export PYTHONPATH=""

cd /data/projects/punim2974/master-project-yenan

nextflow run main.nf \
    --bam /data/projects/punim2974/master-project-yenan/results_hek293t/basecall/calls.bam \
    --control_pod5 /data/projects/punim2974/data/pod5_ivt \
    --reference /data/projects/punim2974/data/references/GRCh38.primary_assembly.genome.fa \
    --ground_truth /data/projects/punim2974/data/references/GSM6432590_293T-mRNA-1_35bp_m2.totalm6A.FDR.csv \
    --outdir results_hek293t \
    -resume