# Nanopore RNA Modification (m6A) Detection Pipeline

A Nextflow DSL2 pipeline for detecting N6-methyladenosine (m6A) RNA modifications from Oxford Nanopore Technologies (ONT) direct RNA sequencing data (RNA004 chemistry). The pipeline processes raw POD5 files and produces standardised bedMethyl output for downstream differential modification abundance analysis.

## Pipeline Overview

| Step | Process | Tool | Description |
|------|---------|------|-------------|
| 1 | BASECALL | Dorado (GPU) | Raw signal to base sequences with m6A modification tags |
| 2 | BASECALL_CONTROL | Dorado (GPU) | Basecall control sample for comparative validation |
| 3 | READ_QC | NanoPlot | Quality assessment of basecalled reads |
| 4 | ALIGN | minimap2 + samtools | Splice-aware alignment to reference genome, preserving MM/ML tags |
| 5 | ALIGN_CONTROL | minimap2 + samtools | Align control sample reads |
| 6 | MODKIT_PILEUP | modkit pileup | Aggregate per-site m6A modification frequencies into bedMethyl |
| 7 | CALIBRATION_SWEEP | calibration_sweep.py | Site-level precision/recall/F1 across frequency thresholds |
| 8 | MODKIT_VALIDATE | modkit validate | Per-read balanced accuracy (ONT benchmark method) |
| 9 | ACCURACY_COMPARISON | glori_correlation.py | Compare pipeline output against GLORI ground truth |

Steps 7, 8, and 9 are conditional: CALIBRATION_SWEEP and MODKIT_VALIDATE run when `--ground_truth` is provided. ACCURACY_COMPARISON runs when a GLORI ground truth CSV is provided.

The pipeline processes both a modified sample (m6A) and an optional control sample (canonical bases) in parallel. When a control sample is provided, `modkit validate` produces a balanced accuracy score directly comparable to ONT's official benchmarks.

## Requirements

- Nextflow >= 25.10.4
- Dorado >= 1.3.1 (requires NVIDIA GPU)
- Conda environment with: NanoPlot, minimap2, samtools, modkit, pandas, matplotlib, scipy

## Benchmarking Results

### ONT Validation Dataset (DRACH context, 20 ground truth sites)

| Metric | Value |
|--------|-------|
| Precision | 0.952 |
| Recall | 1.000 |
| F1 Score | 0.976 |
| Per-read accuracy (modkit validate) | ~97% |

Default filtering: modification type = 'a' (m6A), coverage >= 50, frequency > 80%.

### Real-World Validation: HEK293T (Zou et al. 2025)

Data source: ENA PRJEB80229 (ERR13767652), 5,581,139 reads processed.

Compared against GLORI ground truth (GSM6432590/GSM6432591):

| Metric | Value |
|--------|-------|
| Spearman rho | 0.8995 |
| Pearson r | 0.9203 |
| Pipeline m6A sites detected | 5,389 (coverage >= 50, frequency > 80%) |

The Spearman correlation matches the value reported by Zou et al. (~0.89 for Dorado).

## How to Run

### ONT benchmark dataset

```bash
nextflow run main.nf \
  --pod5 /path/to/m6A_DRACH.pod5 \
  --control_pod5 /path/to/control_DRACH.pod5 \
  --reference /path/to/drach_context_strands.fa \
  --ground_truth /path/to/drach_context_m6A_sites.bed \
  --control_ground_truth /path/to/drach_context_A_sites.bed \
  --outdir results_drach_final \
  -resume
```

### Real-world data (e.g. HEK293T with GLORI ground truth)

```bash
nextflow run main.nf \
  --pod5 /path/to/pod5_directory \
  --reference /path/to/GRCh38.primary_assembly.genome.fa \
  --ground_truth /path/to/GLORI_totalm6A.FDR.csv \
  --outdir results_hek293t \
  -resume
```

### Without ground truth (modification detection only)

```bash
nextflow run main.nf \
  --pod5 /path/to/sample.pod5 \
  --reference /path/to/reference.fa \
  --outdir results
```

## Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--pod5` | Sample POD5 file or directory | null |
| `--control_pod5` | Control sample POD5 file | null |
| `--reference` | Reference FASTA | null |
| `--ground_truth` | Ground truth BED or GLORI CSV | null |
| `--control_ground_truth` | Canonical A ground truth BED | null |
| `--outdir` | Output directory | results |
| `--dorado_model` | Dorado basecalling model | rna004_130bps_sup@v5.1.0 |

## Output

MODKIT_PILEUP produces two bedMethyl files:

- `modifications.bed`: filtered output (m6A only, coverage >= 50, frequency > 80%)
- `all_modifications.bed`: unfiltered output for custom downstream analysis

The 80% modification frequency threshold is a recommended default. Users can adjust this value for their specific use case.

## Project Structure

```
.
├── main.nf                    # Nextflow pipeline definition
├── nextflow.config            # Nextflow configuration (SLURM executor, GPU/CPU partitions)
├── environment.yml            # Conda environment specification
├── run_benchmark.slurm        # SLURM launch script for ONT benchmark dataset
├── run_realworld.slurm        # SLURM launch script for real-world data
├── bin/
│   ├── calibration_sweep.py       # Threshold sweep evaluation (precision/recall/F1)
│   ├── glori_correlation.py       # GLORI ground truth accuracy comparison
│   └── performance_estimation.py  # Site-level performance metrics
└── figures/
    ├── calibration_sweep_combined_drach.png
    ├── calibration_sweep_separate_drach.png
    └── pipeline_workflow.svg
```

## Key Design Decisions

**Splice-aware alignment**: The pipeline uses `minimap2 -ax splice -uf -k14` for alignment to full genome references. This is necessary for real-world RNA-seq data where reads span exon-exon junctions.

**MM/ML tag preservation**: Dorado encodes per-read modification information in MM and ML tags in BAM files. These tags must be explicitly preserved during BAM to FASTQ to BAM conversion using `samtools fastq -T MM,ML` and `minimap2 -y`. Without both flags, modification data is silently lost.

**Coverage filter**: A minimum coverage of 50 reads per site removes false positives caused by high modification frequencies at low-coverage regions. This filter was validated through systematic calibration sweeps.

## References

- Zou, Y., et al. A comparative evaluation of computational models for RNA modification detection using nanopore sequencing with RNA004 chemistry. *Brief Bioinform* 2025;26(4).
- Di Tommaso, P., et al. Nextflow enables reproducible computational workflows. *Nature Biotechnology* 2017;35(4):316-319.