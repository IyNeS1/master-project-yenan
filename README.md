# Nanopore RNA Modification (m6A) Detection Pipeline

A Nextflow DSL2 pipeline for detecting N6-methyladenosine (m6A) RNA modifications from Oxford Nanopore Technologies (ONT) direct RNA sequencing data (RNA004 chemistry). The pipeline processes raw POD5 files or pre-basecalled BAMs and produces standardised bedMethyl output for downstream differential modification analysis.

Current version: **v0.1.0** (alpha, pre-release)

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

The pipeline processes both a modified sample (WT) and an optional control sample (IVT) in parallel using the same process definitions. Steps 5-8 are conditional on ground truth availability. Step 8 runs only when the ground truth is a GLORI CSV file.

## Requirements

- Nextflow >= 25.10.4
- Docker or Singularity/Apptainer
- NVIDIA GPU, only needed when starting from POD5 (basecalling with Dorado)
- Pre-downloaded Dorado models, only needed when starting from POD5

Windows users: run the pipeline inside WSL2 (Ubuntu) with Docker Desktop and its WSL integration enabled.

### Containerised dependencies

All tools run in Docker images hosted at `ghcr.io/iynes1/`:

| Image | Version | Processes |
|-------|---------|-----------|
| dorado | 1.0 (Dorado 1.3.1) | BASECALL |
| minimap2-samtools | 1.0 | ALIGN |
| modkit | 1.0 (v0.6.1) | MODKIT_PILEUP, MODKIT_VALIDATE |
| nanoplot | 1.0 | READ_QC |
| python-scripts | 1.0 | CALIBRATION_SWEEP, CALIBRATION_SWEEP_IVT, ACCURACY_COMPARISON, SUMMARISE |

The Singularity profile passes the host GPU into the Dorado container with `--nv`, and the Docker profile with `--gpus all`.

## Quick Start

### Test run

A small test dataset (ONT DRACH context validation data) is downloaded automatically from ONT's public open data bucket. The test starts from pre-basecalled BAM, so no GPU is needed.

```bash
git clone https://github.com/IyNeS1/master-project-yenan.git
cd master-project-yenan
nextflow run main.nf -profile test,docker
```

Use `-profile test,singularity` on HPC systems without Docker. Results are written to `results_test/`.

Expected results (threshold 80%, 20 ground truth sites):

| Metric | Value |
|--------|-------|
| Precision | 0.909 |
| Recall | 1.000 |
| F1 Score | 0.952 |
| Per-read accuracy (modkit validate) | 97.57% |

The test BAMs were basecalled by ONT with the all-context `m6A_inosine` model, so the results differ slightly from the benchmark below, which uses the `m6A_DRACH` model.

### From POD5 (full pipeline)

```bash
nextflow run main.nf -profile singularity \
  --pod5 /path/to/sample.pod5 \
  --reference /path/to/reference.fa \
  --dorado_models_dir /path/to/dorado_models \
  --outdir results
```

Dorado models can be downloaded in advance with `dorado download --model <model_name> --models-directory /path/to/dorado_models`, so compute nodes do not need internet access.

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
  --dorado_models_dir /path/to/dorado_models \
  --outdir results
```

### Using a locally installed Dorado

If Dorado is already installed on your system, for example as an HPC module, pass its `bin` directory with `--dorado_path`. The Dorado container is then disabled for the basecalling step:

```bash
nextflow run main.nf -profile singularity,spartan \
  --pod5 /path/to/sample.pod5 \
  --reference /path/to/reference.fa \
  --dorado_models_dir /path/to/dorado_models \
  --dorado_path /path/to/dorado/bin \
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
  --pod5 /path/to/sample.pod5 \
  --reference /path/to/reference.fa \
  --dorado_models_dir /path/to/dorado_models \
  --outdir results
```

For small validation datasets, add `-c conf/gpu_short_test.config` to request short GPU jobs, which are usually scheduled much sooner.

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
| `--dorado_models_dir` | Directory with pre-downloaded Dorado models (required with POD5 input) | null |
| `--dorado_path` | `bin` directory of a local Dorado; if set, the Dorado container is not used | null |
| `--dorado_container` | Dorado container image | ghcr.io/iynes1/dorado:1.0 |
| `--kvalue` | minimap2 k-mer size | 14 |
| `--min_coverage` | Minimum read coverage per site | 50 |
| `--mod_threshold` | Modification frequency threshold (%) | 15 |

The default k-mer size of 14 balances memory usage and speed. Accuracy is not significantly affected by k-value choice (F1 difference <= 0.002 across k13, k14, k16). Users with sufficient memory may use `--kvalue 16` for faster alignment.

The default modification threshold of 15% was optimised on real-world HEK293T data. For synthetic datasets with near-complete modification, such as the ONT validation data, a threshold of 80% is used.

## Output

MODKIT_PILEUP produces two bedMethyl files per sample:

- `modifications.bed`: filtered (m6A only, coverage >= min_coverage, frequency > mod_threshold)
- `all_modifications.bed`: unfiltered, for custom downstream analysis

## Benchmarking Results

### ONT Validation Dataset (DRACH context, 20 ground truth sites)

End-to-end run from POD5 with the containerised Dorado (sup@v5.1.0, m6A_DRACH), threshold 80%:

| Metric | Value |
|--------|-------|
| Precision | 0.952 |
| Recall | 1.000 |
| F1 Score | 0.976 |
| Per-read accuracy (modkit validate) | 97.89% |

Results are identical when using a locally installed Dorado.

### Real-World: HEK293T (Zou et al. 2025, GLORI ground truth)

Best configuration: sup@v5.2.0, k14, threshold 15%.

| Metric | Value |
|--------|-------|
| F1 Score | 0.747 |
| Spearman rho | 0.9294 |
| Pearson r | 0.9519 |

Key findings from systematic benchmarking (18 configurations):
- sup models outperform hac by ~4% F1
- k-value has negligible effect on accuracy (F1 difference <= 0.002)
- sup@v5.2.0 achieves the highest F1 across all k-values
- k12 is computationally infeasible (SIGKILL)

## Key Design Decisions

**MM/ML tag preservation**: Dorado encodes modification calls in MM/ML BAM tags. These must be preserved during BAM-FASTQ-BAM conversion using `samtools fastq -T MM,ML` and `minimap2 -y`. Without both flags, modification data is silently lost.

**Coverage filter**: A minimum of 50 reads per site removes false positives from high modification frequencies at low-coverage regions.

**Splice-aware alignment**: Uses `minimap2 -ax splice -uf` for reads spanning exon-exon junctions.

**Container or local Dorado**: Dorado runs in a container by default for portability. A local installation can be selected with `--dorado_path`, and all system-specific paths stay in the environment config (for example `conf/spartan.config`).

## Project Structure

```
main.nf                      Main pipeline entry point
benchmark_from_bam.nf         Batch benchmarking (from BAM)
benchmark.nf                  Legacy benchmarking (from POD5)
nextflow.config               Universal config, profiles and manifest
environment.yml               Conda environment (legacy, not used by profiles)
modules/                      Shared process definitions (9 modules)
conf/                         Environment-specific configs
  test.config                 Test profile (ONT open data)
  spartan.config              Spartan HPC (SLURM, GPU, local Dorado paths)
  gpu_short_test.config       Short GPU requests for small datasets on Spartan
  benchmark.config            ONT validation data paths
  benchmark_realworld_drach.config  Real-world data paths
bin/                          Custom analysis scripts
containers/                   Dockerfiles for each tool
.github/workflows/            CI for container builds
figures/                      Thesis and report figures
```

## Data

The test profile and ONT validation results use the ONT RNA modified base validation dataset (`s3://ont-open-data/rna-modbase-validation_2025.03`), described at https://epi2me.nanoporetech.com/rna-mod-validation-data/.

## References

- Zou, Y., et al. A comparative evaluation of computational models for RNA modification detection using nanopore sequencing with RNA004 chemistry. *Brief Bioinform* 2025;26(4).
- Di Tommaso, P., et al. Nextflow enables reproducible computational workflows. *Nature Biotechnology* 2017;35(4):316-319.