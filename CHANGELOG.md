## Unreleased

- Refactor: migrated from legacy scripts and `setup.py` packaging to a modern `src/` layout with `pyproject.toml` and `uv.lock`.
- CLI: introduced short commands (`add-names`, `heatmap`, `louvain`, `louvain-signed`, `convert-10x`) while keeping the original long command names for compatibility.
- Output conventions: when a plain output prefix is provided (e.g. `Sample`), EEISP now writes outputs under `data/output/` by default.
- Optional features: added a 10X CellRanger converter (`eeisp-convert-10x`) behind an optional `scanpy` extra (`eeisp[scanpy]`).
- Cleanup: removed legacy/duplicated scripts and old packaging artifacts.
- Tests/CI: added a sample-based pytest integration test and a GitHub Actions workflow to run it.

### Benefits

- Reproducible installs via `uv sync` and `uv tool install`.
- Cleaner importable modules (easier to test, extend, and reuse from Python).
- Shorter commands for day-to-day usage, especially with `uv run`.
- Default `data/input` / `data/output` layout matches typical data-pipeline conventions.

## 0.6.0 (2024-07-06)
- Added `LouvainSigned.py` and `network.py`

## 0.5.0 (2023-10-21)
- Added this ChangeLog
- Added a function to plot the histogram of CDI and EEI

## 0.4.1 (2022-07-06)
- Dramatically improved EEISP. It now uses much less memory and is faster than the previous version.
- Use multiple CPUs to compute the joint_matrix to reduce computation time.

## 0.3.0 (2021-09-25)
- First release
- Released on PyPI
- Updated README


