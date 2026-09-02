# Project Structure

## Pipeline Files

- `main.nf` — Main pipeline: basecalling, alignment, modkit pileup, calibration sweep, modkit validate, GLORI correlation. Supports --pod5 or --bam input for WT, --control_pod5 or --control_bam for IVT control. All processes imported from modules/.
- `nextflow.config` — Universal configuration: default parameters, container definitions, profiles (singularity, docker, spartan).
- `benchmark_from_bam.nf` — Benchmarking pipeline starting from pre-basecalled BAMs. Runs multiple model x k-value combinations. Reuses modules from modules/.
- `benchmark.nf` — (Legacy) End-to-end benchmarking from POD5. Not yet refactored to use modules.
- `benchmark.config` — (Legacy) Config for benchmark.nf with ONT DRACH context validation data paths.
- `benchmark_realworld_drach.config` — (Legacy) Config for benchmark_from_bam.nf with real-world HEK293T data paths.
- `environment.yml` — Conda environment definition (backup for systems without container support).

## modules/

Shared Nextflow process definitions, imported by main.nf and benchmark_from_bam.nf.

- `basecall.nf` — Dorado basecalling with integrated modification detection (GPU)
- `read_qc.nf` — Read quality control with NanoPlot
- `align.nf` — Splice-aware alignment with minimap2, k-value passed via channel
- `modkit_pileup.nf` — Modification aggregation with modkit pileup, outputs filtered and unfiltered bedMethyl
- `calibration_sweep.nf` — Threshold sweep for WT data with DRACH motif and GLORI overlap filtering
- `calibration_sweep_ivt.nf` — Threshold sweep for IVT negative control false positive analysis
- `modkit_validate.nf` — Per-read accuracy assessment using modkit validate
- `accuracy_comparison.nf` — Spearman/Pearson correlation against GLORI ground truth
- `summarise.nf` — Aggregates calibration sweep CSVs into benchmark summary

## conf/

Environment-specific configuration files, loaded via -profile.

- `spartan.config` — University of Melbourne Spartan HPC: SLURM executor, resource allocation, Dorado module paths

## bin/

Custom scripts called by pipeline processes.

- `calibration_sweep.py` — Main calibration sweep script for WT data
- `calibration_sweep_ivt.py` — Calibration sweep for IVT negative control data
- `fp_trend_plot_thesis.py` — Generates FP trend plots for thesis figures
- `glori_correlation.py` — Spearman/Pearson correlation analysis between pipeline output and GLORI ground truth
- `test_calibration_sweep.py` — Unit tests for calibration sweep logic

## containers/

Docker/Apptainer containerisation files.

- `minimap2-samtools/` — Dockerfile for minimap2 2.30 + samtools 1.23.1
- `modkit/` — Dockerfile for ont-modkit 0.6.1
- `nanoplot/` — Dockerfile for NanoPlot (plotly<6 pinned)
- `python-scripts/` — Dockerfile for pandas 3.0.2, matplotlib, pysam, scipy
- `sif/` — Apptainer SIF image files for Spartan HPC

Images pushed to ghcr.io/iynes1/. GitHub Actions workflow in .github/ triggers builds on containers/** changes.

## Results Directories

- `results_drach_final/` — Final results from main pipeline on ONT DRACH context validation data (F1 = 0.976)
- `results_benchmark_realworld_drach/` — Real-world benchmarking results (HEK293T WT, Zou et al. 2025, GLORI ground truth). 18 model x k-value combinations.
- `results_benchmark/` — Early benchmarking results from ONT validation data
- `results_stacked_bar/` — Runtime and memory benchmarking data. Contains timing_log.tsv (24 rows).
- `basecall_results/` — Dorado basecalling output BAMs for each model (hac/sup x v5.1.0/v5.2.0)

## figures/

Thesis and report figures.

- `drach_context/` — Figures from DRACH context validation dataset
- `hek293t/` — Figures from real-world HEK293T dataset
- `sampled_context/` — Figures from sampled context validation dataset (early development)
- `fp_trend_ivt_thesis.pdf/.png` — IVT FP trend plot for thesis
- `fp_trend_ivt_vs_wt.png` — IVT vs WT false positive comparison
- `pipeline_workflow.svg` — Pipeline workflow diagram

## slurm_scripts/

SLURM job submission scripts for Spartan HPC.

## slurm_logs/

Retained SLURM job output logs for final runs only.

## Other

- `.github/` — GitHub Actions workflows (docker-build.yml)
- `software/` — Local tool installations: dorado_models/, dorado-0.9.1-linux-x64/
- `README.md` — Project readme
- `.gitignore` — Excludes results directories, slurm logs, .nextflow cache
