# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Running Scripts

Always use the venv Python interpreter — it has Biopython installed:

```bash
# Run any script
.venv/Scripts/python.exe <script>.py

# Local_Alignment .py has a space in its filename — always quote it
.venv/Scripts/python.exe "Local_Alignment .py"
```

If on Windows and seeing `UnicodeEncodeError`, prefix with:
```bash
set PYTHONIOENCODING=utf-8
```

## Environment

- Python venv at `.venv/` — Biopython, numpy, pandas installed there
- Windows cp1252 terminal encoding — avoid non-ASCII characters (arrows, box-drawing, etc.) in `print()` strings
- `Data_Retrieval.py` requires an active internet connection (NCBI + UniProt); all other scripts are fully offline

## Architecture

This is a DCIT 411 bioinformatics coursework project structured as five tasks:

**Task 2 — Data retrieval** (`Data_Retrieval.py`): Fetches sequences from NCBI Entrez (`Bio.Entrez.efetch`) and UniProt (`Bio.ExPASy`). Exports `human_hbb_seq` and `mouse_hbb_seq` which `different_parameters.py` imports at runtime (with a stub fallback if the import fails).

**Task 3 — Pairwise alignment** (`global_alignment_demo.py`, `Local_Alignment .py`, `Pairwise_sequence.py`, `different_parameters.py`): All use `Bio.Align.PairwiseAligner`. DNA uses numeric match/mismatch scores; protein uses BLOSUM62 from `Bio.Align.substitution_matrices`. `different_parameters.py` also contains legacy `Bio.Align.Applications` wrappers for ClustalW/MUSCLE/MAFFT — these are kept for reference alongside the newer subprocess-based wrappers in `MSA.py`.

**Task 4 — MSA** (`MSA.py`): Detects external tools via `shutil.which()` and tries ClustalW → MUSCLE → MAFFT in order. If none are found, falls back to a pure-Python center-star progressive aligner (`fallback_msa`) that uses `PairwiseAligner` + BLOSUM62. Requires `hemoglobin.fasta` to be present. Outputs `fallback_aln.fasta`.

**Task 5 — Advanced topics** (`advanced_topics.py`): Fully standalone. Three demos in sequence:
- PSSM built manually from a `MultipleSeqAlignment` (log-odds against uniform background); `Bio.Align.AlignInfo.SummaryInfo` is deprecated — the code uses a manual implementation instead
- Kabsch structural superimposition using numpy SVD (falls back to centroid-only without numpy)
- Consensus generation with Shannon entropy per column; `SummaryInfo.dumb_consensus` is called first but the code always runs the manual fallback regardless

## Key Gotchas

- `Bio.Align.AlignInfo.SummaryInfo` is deprecated in current Biopython — do not add new usage of it; use manual column-iteration instead
- `different_parameters.py` has duplicate MSA runner functions near the bottom — these are legacy artifacts, not a bug
- `Local_Alignment .py` has a literal space in its filename — preserve it, it matches the existing file
