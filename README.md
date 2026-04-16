# Nanopore RNA Modification Pipeline

A Nextflow pipeline for detecting RNA modifications (m6A) from Oxford Nanopore sequencing data.

## Pipeline overview

| Step | Process                                 | Tool                                                               |
| ---- | --------------------------------------- | ------------------------------------------------------------------ |
| 1    | Raw data QC                             | _(planned: LongReadSum / ToulligQC)_                               |
| 2    | Basecalling                             | _(planned: Dorado - not yet implemented, Intel Mac not supported)_ |
| 3    | Read QC                                 | NanoPlot 1.46.2                                                    |
| 4    | Sequence alignment                      | Minimap2 2.30 + Samtools 1.23                                      |
| 5/6  | Modification calling + bedMethyl output | modkit                                                             |

## Current status

Steps 3-6 are functional. Steps 1-2 are planned but not yet implemented.
Pre-basecalled BAM files from ONT open data are used as pipeline input for now.

## Requirements

- Nextflow 25.10.4
- Conda
- conda environment: `nanopore_env` with NanoPlot, Minimap2, Samtools, modkit, pandas

> **Note:** Before running, update `nextflow.config` with your own conda environment path:
>
> ```
> process.conda = "/path/to/your/conda/envs/nanopore_env"
> ```

## Test data

From ONT's RNA Modified Base Best Practices and Benchmarking dataset:

- https://epi2me.nanoporetech.com/rna-mod-validation-data/

Download with:

```bash
aws s3 cp --no-sign-request \
  s3://ont-open-data/rna-modbase-validation_2025.03/basecalls/m6A_rep1.bam ./

aws s3 cp --no-sign-request \
  s3://ont-open-data/rna-modbase-validation_2025.03/basecalls/control_rep1.bam ./

aws s3 cp --no-sign-request \
  s3://ont-open-data/rna-modbase-validation_2025.03/references/sampled_context_strands.fa ./

aws s3 cp --no-sign-request \
  s3://ont-open-data/rna-modbase-validation_2025.03/references/sampled_context_m6A_sites.bed ./
```

## How to run

**Pipeline:**

```bash
nextflow run main.nf \
  --bam /path/to/your.bam \
  --reference /path/to/reference.fa \
  --outdir /path/to/results
```

**Performance estimation:**

```bash
python performance_estimation.py \
  --m6a /path/to/results/modkit/modifications.bed \
  --ground_truth /path/to/ground_truth.bed
```

## Performance estimation results (test dataset)

| Metric    | Value |
| --------- | ----- |
| TP        | 29    |
| FP        | 5     |
| FN        | 1     |
| Precision | 0.853 |
| Recall    | 0.967 |
| F1 Score  | 0.906 |

m6A sites were filtered by modification frequency > 80% and modification type 'a' only.

## Calibration sweep

To find the optimal frequency threshold, run the calibration sweep script:

```bash
python calibration_sweep.py \
  --m6a /path/to/modkit_output/m6A_default.bed \
  --ground_truth /path/to/references/ground_truth.bed
```

This will test frequency thresholds from 70% to 100% at 1% intervals and generate two plots:

- `calibration_sweep_combined.png` — Precision, Recall and F1 Score on one graph
- `calibration_sweep_separate.png` — separate plots for each metric

## Calibration sweep results (test dataset)

Best threshold: **80%**

| Threshold | Precision | Recall | F1 Score     |
| --------- | --------- | ------ | ------------ |
| 70%       | 0.770     | 1.000  | 0.870        |
| 75%       | 0.806     | 0.967  | 0.879        |
| 80%       | 0.853     | 0.967  | **0.906** ✅ |
| 85%       | 0.853     | 0.967  | 0.906        |
| 90%       | 0.857     | 0.800  | 0.828        |
