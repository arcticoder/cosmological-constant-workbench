# Cosmological Constant Workbench

[![GitHub](https://img.shields.io/badge/GitHub-DawsonInstitute%2Fcosmological--constant--workbench-blue)](https://github.com/DawsonInstitute/cosmological-constant-workbench)

A reproducible, test-first workbench for the cosmological constant problem:

- Reproduce the baseline discrepancy between naive QFT vacuum energy estimates and the observed dark-energy density.
- Provide small, explicit calculators (with units) and a CLI that outputs JSON.
- Create a place to evaluate candidate mechanisms *by their assumptions* and *their observable consequences*, not by curve-fitting.
- Integrate external LQG-inspired predictors and perform bounded parameter scans to establish empirical constraints (e.g., 150+ order mismatch demonstrated in [papers/lqg_cc_constraints](papers/lqg_cc_constraints)).

## Quick start

```bash
cd cosmological-constant-workbench
python -m venv .venv && source .venv/bin/activate
pip install -e .

# Baseline: observed ρ_Λ, Λ, and naive QFT cutoffs
ccw-baseline --h0 67.4 --omega-lambda 0.6889 --cutoff planck
ccw-baseline --cutoff electroweak
ccw-baseline --cutoff 1TeV

# Sweep toy mechanism parameters (writes JSON + CSV)
ccw-sweep --mechanism cpl --z 0,1,2 --grid-json '[{"w0":-1.0,"wa":0.2},{"w0":-0.9,"wa":0.0}]'
```

## What this repo is (and isn’t)

- This is **not** a claim of a solved cosmological constant problem.
- It is a **rigorous, incremental** environment to (1) reproduce the problem, (2) encode assumptions in code, and (3) compare mechanisms against constraints.

## Repository layout

- `src/ccw/` — library code
- `scripts/` — runnable scripts (thin wrappers)
- `tests/` — unit tests (pytest)
- `docs/` — roadmap, notes, and provenance

## Next steps

See `docs/TASKS.md` for the working task list.

## Provenance

See `docs/PROVENANCE.md` for equation provenance and scope.
