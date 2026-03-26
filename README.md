# master-project

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
- conda environment: `nanopore_env` with NanoPlot, Minimap2, Samtools, modkit

## Test data

From ONT's RNA Modified Base Best Practices and Benchmarking dataset:

- https://epi2me.nanoporetech.com/rna-mod-validation-data/

Download with:

```bash
aws s3 cp --no-sign-request \
  s3://ont-open-data/rna-modbase-validation_2025.03/basecalls/m6A_rep1.bam ./

aws s3 cp --no-sign-request \
  s3://ont-open-data/rna-modbase-validation_2025.03/references/sampled_context_strands.fa ./
```

## How to run

```bash
nextflow run main.nf \
  --bam /path/to/your.bam \
  --reference /path/to/reference.fa \
  --outdir /path/to/results
```

## Output

- `results/qc/` - NanoPlot QC report
- `results/aligned/` - aligned BAM file
- `results/modkit/` - modifications.bed (bedMethyl format)
