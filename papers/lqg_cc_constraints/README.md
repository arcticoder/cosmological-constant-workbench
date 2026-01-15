# LQG CC Constraints (Short Paper)

This folder contains a short manuscript and a reproducible data/figure pipeline supporting an empirical constraint rememberable as a “no-go style” bound for the current LQG CC predictor formulation in a bounded parameter domain.

## Reproduce the scan data

From the repo root:

- `PYTHONPATH=src python papers/lqg_cc_constraints/generate_scan_data.py`

This writes:

- `papers/lqg_cc_constraints/data/scan_results.tsv`

## Build the PDF

- `bash papers/lqg_cc_constraints/build_paper.sh`

This regenerates the TSV and compiles `paper.tex` using `latexmk` (preferred) or `pdflatex`.

## Files

- `paper.tex`: manuscript (pgfplots histogram reads `data/scan_results.tsv`)
- `generate_scan_data.py`: bounded parameter scan → stable TSV export
- `data/scan_results.tsv`: scan outputs (108 points)
