# Nanopore RNA Modification (m6A) Detection Pipeline

A Nextflow DSL2 pipeline for detecting N6-methyladenosine (m6A) RNA modifications from Oxford Nanopore Technologies (ONT) direct RNA sequencing data (RNA004 chemistry). The pipeline processes raw POD5 files or pre-basecalled BAMs and produces standardised bedMethyl output for downstream differential modification analysis.

## Pipeline Overview

The pipeline consists of 9 modular processes shared across analysis and benchmarking workflows. Each process is defined once in `modules/` and reused via Nextflow `include` aliases.

| Step | Process | Tool | Description |
|------|---------|------|-------------|
| 1 | BASECALL | Dorado (GPU) | Raw signal to base sequences with m6A modification tags |
| 2 | READ_QC | NanoPlot | Quality assessment of basecalled reads |
| 3 | ALIGN | minimap2 + samtools | Splice-aware alignment preserving MM/ML tags |
| 4 | MODKIT_PILEUP | modkit pileup | Aggregate per-site m6A modification frequencies into bedMethyl |
| 5 | CALIBRATION_SWEEP | calibration_sweep.py | Site-level precision/recall/F1 across frequency thresholds |
| 6 | CALIBRATION_SWEEP_IVT | calibration_sweep_ivt.py | False positive analysis for IVT negative control |
| 7 | MODKIT_VALIDATE | modkit validate | Per-read balanced accuracy (ONT benchmark method) |
| 8 | ACCURACY_COMPARISON | glori_correlation.py | Spearman/Pearson correlation against GLORI ground truth |
| 9 | SUMMARISE | summarise.nf | Aggregate benchmark results across configurations |

The pipeline processes both a modified sample (WT) and an optional control sample (IVT) in parallel using the same process definitions. Steps 5-8 are conditional on ground truth availability.

## Requirements

- Nextflow >= 25.10.4
- One of: Docker, Singularity/Apptainer, or Conda
- Dorado >= 1.3.1 (requires NVIDIA GPU, only needed when starting from POD5)

### Containerised dependencies (recommended)

Four Docker images are available at `ghcr.io/iynes1/`:

| Image | Version | Processes |
|-------|---------|-----------|
| minimap2-samtools | 1.0 | ALIGN |
| modkit | 1.0 (v0.6.1) | MODKIT_PILEUP, MODKIT_VALIDATE |
| nanoplot | 1.0 | READ_QC |
| python-scripts | 1.0 | CALIBRATION_SWEEP, ACCURACY_COMPARISON, SUMMARISE |

## Quick Start

### From POD5 (full pipeline)

```bash
nextflow run main.nf -profile singularity \
  --pod5 /path/to/sample.pod5 \
  --reference /path/to/reference.fa \
  --outdir results
```

### From pre-basecalled BAM (skip basecalling)

```bash
nextflow run main.nf -profile singularity \
  --bam /path/to/calls.bam \
  --reference /path/to/reference.fa \
  --outdir results
```

### With IVT control and ground truth

```bash
nextflow run main.nf -profile singularity \
  --pod5 /path/to/wt.pod5 \
  --control_pod5 /path/to/ivt.pod5 \
  --reference /path/to/reference.fa \
  --ground_truth /path/to/ground_truth.bed \
  --control_ground_truth /path/to/control_ground_truth.bed \
  --outdir results
```

### Batch benchmarking (multiple models x k-values)

```bash
nextflow run benchmark_from_bam.nf -profile singularity \
  --bam_dir /path/to/basecall_results/ \
  --reference /path/to/reference.fa \
  --ground_truth /path/to/ground_truth.csv \
  --kvalues "13,14,16" \
  --outdir results_benchmark
```

### Spartan HPC (University of Melbourne)

```bash
nextflow run main.nf -profile singularity,spartan \
  --bam /path/to/calls.bam \
  --reference /path/to/reference.fa \
  --outdir results
```

## Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--pod5` | Sample POD5 file or directory | null |
| `--bam` | Pre-basecalled BAM (skips basecalling) | null |
| `--control_pod5` | IVT control POD5 | null |
| `--control_bam` | IVT control BAM | null |
| `--reference` | Reference genome FASTA | null |
| `--ground_truth` | Ground truth BED or GLORI CSV | null |
| `--control_ground_truth` | Canonical A ground truth for IVT | null |
| `--outdir` | Output directory | results |
| `--dorado_model` | Dorado basecalling model | rna004_130bps_sup@v5.1.0 |
| `--kvalue` | minimap2 k-mer size | 14 |
| `--min_coverage` | Minimum read coverage per site | 50 |
| `--mod_threshold` | Modification frequency threshold (%) | 15 |

The default k-mer size of 14 balances memory usage and speed. Accuracy is not significantly affected by k-value choice (F1 difference <= 0.002 across k13, k14, k16). Users with sufficient memory may use `--kvalue 16` for faster alignment.

## Output

MODKIT_PILEUP produces two bedMethyl files per sample:

- `modifications.bed` — filtered (m6A only, coverage >= min_coverage, frequency > mod_threshold)
- `all_modifications.bed` — unfiltered, for custom downstream analysis

## Benchmarking Results

### ONT Validation Dataset (DRACH context, 20 ground truth sites)

| Metric | Value |
|--------|-------|
| Precision | 0.952 |
| Recall | 1.000 |
| F1 Score | 0.976 |
| Per-read accuracy | ~97% |

### Real-World: HEK293T (Zou et al. 2025, GLORI ground truth)

Best configuration: sup@v5.2.0, k14, threshold 15%.

| Metric | Value |
|--------|-------|
| F1 Score | 0.747 |
| Spearman rho | 0.8995 |
| Pearson r | 0.9203 |

Key findings from systematic benchmarking (18 configurations):
- sup models outperform hac by ~4% F1
- k-value has negligible effect on accuracy (F1 difference <= 0.002)
- sup@v5.2.0 achieves the highest F1 across all k-values
- k12 is computationally infeasible (SIGKILL)

## Key Design Decisions

**MM/ML tag preservation**: Dorado encodes modification calls in MM/ML BAM tags. These must be preserved during BAM-FASTQ-BAM conversion using `samtools fastq -T MM,ML` and `minimap2 -y`. Without both flags, modification data is silently lost.

**Coverage filter**: A minimum of 50 reads per site removes false positives from high modification frequencies at low-coverage regions.

**Splice-aware alignment**: Uses `minimap2 -ax splice -uf` for reads spanning exon-exon junctions.

## Project Structure

```
main.nf                      Main pipeline entry point
benchmark_from_bam.nf         Batch benchmarking (from BAM)
benchmark.nf                  Legacy benchmarking (from POD5)
nextflow.config               Universal config and profiles
environment.yml               Conda fallback
modules/                      Shared process definitions (9 modules)
conf/                         Environment-specific configs
  spartan.config              Spartan HPC (SLURM, Dorado paths)
  benchmark.config            ONT validation data paths
  benchmark_realworld_drach.config  Real-world data paths
bin/                          Custom analysis scripts
containers/                   Dockerfiles for each tool
.github/workflows/            CI for container builds
figures/                      Thesis and report figures
```

## References

- Zou, Y., et al. A comparative evaluation of computational models for RNA modification detection using nanopore sequencing with RNA004 chemistry. *Brief Bioinform* 2025;26(4).
- Di Tommaso, P., et al. Nextflow enables reproducible computational workflows. *Nature Biotechnology* 2017;35(4):316-319.
