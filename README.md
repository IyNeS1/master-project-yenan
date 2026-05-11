# Nanopore RNA Modification (m6A) Detection Pipeline

A Nextflow DSL2 pipeline for detecting m6A RNA modifications from Oxford Nanopore direct RNA sequencing data, with benchmarking against ONT's official validation dataset.

## Pipeline Overview

| Step | Process | Tool | Description |
|------|---------|------|-------------|
| 1 | Basecalling | Dorado (GPU) | Raw signal → base sequences with m6A modification tags |
| 2 | Read QC | NanoPlot | Quality assessment of basecalled reads |
| 3 | Sequence Alignment | minimap2 + samtools | Map reads to reference, preserving mod tags (MM/ML) |
| 4 | Modification Detection | modkit pileup | Aggregate per-site m6A modification frequencies |
| 5a | Calibration Sweep | calibration_sweep.py | Site-level evaluation across frequency thresholds |
| 5b | Modkit Validate | modkit validate | Per-read balanced accuracy (ONT benchmark method) |

The pipeline processes both a **modified sample** (m6A) and an optional **control sample** (canonical bases) in parallel. When a control sample is provided, `modkit validate` produces a balanced accuracy score directly comparable to ONT's official benchmarks.

## Requirements

- Nextflow ≥ 25.10.4
- Dorado ≥ 1.3.1 (requires NVIDIA GPU)
- Conda environment `nanopore_env` with: NanoPlot, minimap2, samtools, modkit, pandas, matplotlib

## Test Data

From ONT's [RNA Modified Base Best Practices and Benchmarking](https://epi2me.nanoporetech.com/rna-mod-validation-data/) dataset:

```bash
# Modified sample (m6A) — subset
aws s3 cp --no-sign-request \
  s3://ont-open-data/rna-modbase-validation_2025.03/subset/m6A_rep1.pod5 ./

# Control sample — subset
aws s3 cp --no-sign-request \
  s3://ont-open-data/rna-modbase-validation_2025.03/subset/control_rep1.pod5 ./

# Modified sample (m6A) — full dataset (64GB)
aws s3 cp --no-sign-request \
  s3://ont-open-data/rna-modbase-validation_2025.03/full/m6A_rep1.pod5 ./

# Reference and ground truth
aws s3 cp --no-sign-request \
  s3://ont-open-data/rna-modbase-validation_2025.03/references/sampled_context_strands.fa ./
aws s3 cp --no-sign-request \
  s3://ont-open-data/rna-modbase-validation_2025.03/references/sampled_context_m6A_sites.bed ./
aws s3 cp --no-sign-request \
  s3://ont-open-data/rna-modbase-validation_2025.03/references/sampled_context_A_sites.bed ./
```

## How to Run

### Full pipeline (from POD5, with control)

```bash
nextflow run main.nf \
  --pod5 /path/to/m6A_rep1.pod5 \
  --control_pod5 /path/to/control_rep1.pod5 \
  --reference /path/to/sampled_context_strands.fa \
  --ground_truth /path/to/sampled_context_m6A_sites.bed \
  --control_ground_truth /path/to/sampled_context_A_sites.bed \
  --outdir results
```

### From pre-basecalled BAM (skip basecalling)

```bash
nextflow run main.nf \
  --bam /path/to/m6A_rep1.bam \
  --control_bam /path/to/control_rep1.bam \
  --reference /path/to/sampled_context_strands.fa \
  --ground_truth /path/to/sampled_context_m6A_sites.bed \
  --control_ground_truth /path/to/sampled_context_A_sites.bed \
  --outdir results
```

### Without control (calibration sweep only, no validate)

```bash
nextflow run main.nf \
  --pod5 /path/to/m6A_rep1.pod5 \
  --reference /path/to/sampled_context_strands.fa \
  --ground_truth /path/to/sampled_context_m6A_sites.bed \
  --outdir results
```

## Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--pod5` | m6A sample POD5 file | null |
| `--bam` | Pre-basecalled m6A BAM (skip basecalling) | null |
| `--control_pod5` | Control sample POD5 file | null |
| `--control_bam` | Pre-basecalled control BAM | null |
| `--reference` | Reference FASTA | null |
| `--ground_truth` | m6A ground truth BED | null |
| `--control_ground_truth` | Canonical A ground truth BED | null |
| `--outdir` | Output directory | results |
| `--dorado_model` | Dorado basecalling model | rna004_130bps_sup@v5.1.0 |

## Results

### modkit validate (per-read balanced accuracy)

| Dataset | Filtered Accuracy | False Positive Rate | False Negative Rate |
|---------|-------------------|---------------------|---------------------|
| Subset (m6A + control) | 97.13% | 0.39% | 5.63% |
| Full (m6A + control) | 95.96% | 0.42% | 8.10% |
| ONT official benchmark | 97.21% | ~0.38% | ~0.98% |

### Configuration

| Component | This Pipeline | ONT Official |
|-----------|---------------|--------------|
| Dorado | v1.3.1 | v0.9.1 |
| Model | rna004_130bps_sup@v5.1.0 | SUP v5.1 |
| modkit | validate | validate |
| GPU | NVIDIA A100 80GB | Tesla V100 16GB |

## Project Structure

```
.
├── main.nf                 # Nextflow pipeline
├── nextflow.config         # Nextflow configuration
├── run_pipeline.slurm      # SLURM job submission script
├── bin/
│   ├── calibration_sweep.py    # Threshold sweep evaluation
│   └── performance_estimation.py
├── results/                # Subset results
├── results_full/           # Full dataset results
└── figures/                # Calibration sweep plots
```
