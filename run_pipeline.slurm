#!/bin/bash
#SBATCH --job-name=nanopore_pipeline
#SBATCH --partition=gpu-a100
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --gpus=1
#SBATCH --time=24:00:00
#SBATCH --output=slurm_%j.out
#SBATCH --error=slurm_%j.err
#SBATCH --account=punim2974

# Load Nextflow and conda
module load Nextflow/25.10.4
module load Mamba/23.1.0-4

# Activate conda environment
source activate nanopore_env

# Clear system Python paths to avoid conflicts with conda
export PYTHONPATH=""

# Navigate to project directory
cd /data/projects/punim2974/master-project-yenan

# Add dorado and CUDA libraries to PATH without loading GCC/OpenMPI modules
# This avoids Python/numpy conflicts with conda environment
export PATH=/apps/easybuild-2022/easybuild/software/MPI/GCC/11.3.0/OpenMPI/4.1.4/dorado/1.3.1-CUDA-12.2.0/bin:$PATH
export LD_LIBRARY_PATH=/apps/easybuild-2022/easybuild/software/MPI/GCC/11.3.0/OpenMPI/4.1.4/dorado/1.3.1-CUDA-12.2.0/lib:/apps/easybuild-2022/easybuild/software/Core/CUDA/12.2.0/lib:/apps/easybuild-2022/easybuild/software/Core/CUDA/12.2.0/nvvm/lib64:/apps/easybuild-2022/easybuild/software/Core/cuDNN/8.9.3.28-CUDA-12.2.0/lib:/apps/easybuild-2022/easybuild/software/Core/CUDA/12.2.0/extras/CUPTI/lib64:$LD_LIBRARY_PATH

# Navigate to project directory
cd /data/projects/punim2974/master-project-yenan

# Run pipeline (resume from last successful step)
nextflow run main.nf \
    --pod5 /data/projects/punim2974/data/pod5/m6A_rep1_full.pod5 \
    --control_pod5 /data/projects/punim2974/data/pod5/control_rep1.pod5 \
    --reference /data/projects/punim2974/data/references/sampled_context_strands.fa \
    --ground_truth /data/projects/punim2974/data/references/sampled_context_m6A_sites.bed \
    --control_ground_truth /data/projects/punim2974/data/references/sampled_context_A_sites.bed \
    --outdir results_full
