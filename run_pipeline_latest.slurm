#!/bin/bash
#SBATCH --job-name=nf_launch
#SBATCH --partition=sapphire
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=48:00:00
#SBATCH --output=slurm_%j.out
#SBATCH --error=slurm_%j.err
#SBATCH --account=punim2974

module load Nextflow/25.10.4

cd /data/projects/punim2974/master-project-yenan

nextflow run main.nf \
    --pod5 /data/projects/punim2974/data/pod5/m6A_DRACH.pod5 \
    --control_pod5 /data/projects/punim2974/data/pod5/control_DRACH.pod5 \
    --reference /data/projects/punim2974/data/references/drach_context_strands.fa \
    --ground_truth /data/projects/punim2974/data/references/drach_context_m6A_sites.bed \
    --control_ground_truth /data/projects/punim2974/data/references/drach_context_A_sites.bed \
    --outdir results_drach_final \
    -resume