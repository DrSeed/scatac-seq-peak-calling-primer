# scatac-seq-peak-calling-primer

Single-cell ATAC-seq promises chromatin accessibility at single-cell resolution. What it delivers is a brutally sparse matrix where most entries are zero and you spend half your time fighting the noise. This repo is a teaching scaffold that walks you from a fragments file to a cell-by-peak matrix, then quality-controls the wreckage.

## Why This Exists

Most tutorials hand you a pre-built matrix and pretend the hard part is clustering. The hard part is everything before that: defining peaks across cells that each have a handful of reads, deciding which barcodes are real cells, and not deceiving yourself about what the binary accessibility signal actually means.

This project demonstrates the core moves in your own words, with runnable Python that fakes a small fragments dataset so you can see the pipeline behave without downloading 40 GB.

## The Pipeline in Plain Terms

1. Read a fragments file (chrom, start, end, barcode, count).
2. Filter barcodes by fragment count to separate cells from empty droplets.
3. Build candidate peaks by binning the genome and counting overlapping fragments.
4. Construct a sparse cell-by-peak matrix.
5. Compute per-cell QC: total fragments, fraction in peaks (FRiP), peak counts.

## When NOT to Use This

| Situation | Use this scaffold? | Reach for instead |
|-----------|--------------------|-------------------|
| Learning the data flow | Yes | - |
| Production atlas-scale data | No | ArchR, SnapATAC2, signac |
| You need TF motif footprinting | No | chromVAR, TOBIAS |
| You only have a count matrix already | No | scanpy / muon |

## The Uncomfortable Truth

A cell-by-peak matrix from scATAC-seq is roughly 1-3% non-zero. Almost every clever downstream method exists to cope with that emptiness, not to extract hidden richness. If your biological conclusion survives only because of an aggressive imputation step, it was never really there. Treat binarisation as the honest default and be suspicious of anything that smells like wishful smoothing.

## Quick Start

```bash
python scatac_pipeline.py --simulate --out results/
```

This simulates fragments, builds the matrix, and writes QC to `results/`.

## Further Reading

Inspired by Ming 'Tommy' Tang, "Computational Analyses and Challenges of Single-cell ATAC-seq" (https://divingintogeneticsandgenomics.com/publication/2025-12-01-scatac-seq-review/).
