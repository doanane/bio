# DCIT 411: Bioinformatics — Sequence Alignment with Biopython

A complete implementation of pairwise and multiple sequence alignment
techniques using Biopython, covering dynamic programming algorithms,
substitution matrices, MSA tools, and advanced profile/structural alignment.

---

## Table of Contents

1. [Project Overview](#project-overview)
2. [Repository Structure](#repository-structure)
3. [Requirements](#requirements)
4. [Setup and Installation](#setup-and-installation)
5. [File-by-File Guide](#file-by-file-guide)
6. [Data Files](#data-files)
7. [How to Run Each Script](#how-to-run-each-script)
8. [Optional: External MSA Tools](#optional-external-msa-tools)
9. [Expected Output Summary](#expected-output-summary)
10. [Troubleshooting](#troubleshooting)

---

## Project Overview

This project implements and demonstrates sequence alignment in Python
using Biopython, structured around five tasks:

| Task | Description |
|------|-------------|
| 1 | Literature review (background knowledge — no code) |
| 2 | Data retrieval from NCBI and UniProt, sequence preprocessing |
| 3 | Pairwise alignment: Needleman-Wunsch (global) and Smith-Waterman (local) |
| 4 | Multiple Sequence Alignment (MSA): ClustalW, MUSCLE, MAFFT, built-in fallback |
| 5 | Advanced topics: PSSM/HMM profile alignment, structural alignment, consensus generation |

All scripts run fully offline. External MSA tools (ClustalW, MUSCLE, MAFFT)
are optional — every script has a built-in fallback.

---

## Repository Structure

```
bio/
|
|-- Data_Retrieval.py          # Task 2  — fetch & preprocess sequences from NCBI / UniProt
|-- global_alignment_demo.py   # Task 3a — Needleman-Wunsch global alignment
|-- Local_Alignment .py        # Task 3a — Smith-Waterman local alignment
|-- Pairwise_sequence.py       # Task 3b — BLOSUM62 substitution matrix configuration
|-- different_parameters.py    # Task 3c — gap penalty experiments, parameter comparison
|-- MSA.py                     # Task 4  — MSA with ClustalW/MUSCLE/MAFFT + built-in fallback
|-- advanced_topics.py         # Task 5  — PSSM/HMM profile alignment, structural alignment,
|                              #           consensus sequence generation
|
|-- hemoglobin.fasta           # Input FASTA: 5 hemoglobin alpha-chain sequences (Human,
|                              #   Mouse, Rabbit, Dog, Horse) used by MSA.py
|-- NC_005816.fna              # FASTA genome: Yersinia pestis plasmid pPCP1
|-- NC_005816.gb               # GenBank record for NC_005816
|-- mafft_aln.fasta            # Example MAFFT output alignment (pre-generated)
|-- ls_orchid.fasta            # Orchid 18S rRNA sequences (used in tutorial notebooks)
|-- ls_orchid.gbk              # GenBank records for orchid sequences
|-- example.fastq              # Example FASTQ file for SeqIO tutorial
|
|-- tutorial1.ipynb            # Biopython basics: Seq, SeqRecord, SeqIO
|-- 02-seq_objects.ipynb       # Seq object methods and properties
|-- 03-seq_annot.ipynb         # Sequence annotations and SeqRecord features
|-- 04-seqio.ipynb             # Reading/writing sequence file formats
|-- class_quiz.ipynb           # In-class exercises
|-- Untitled.ipynb             # Scratch notebook
|
|-- .venv/                     # Python virtual environment (not committed to git)
|-- __pycache__/               # Python bytecode cache (auto-generated)
```

---

## Requirements

### Python version

Python 3.9 or higher is required. Python 3.11+ is recommended.

### Python packages

| Package | Version | Purpose |
|---------|---------|---------|
| `biopython` | >= 1.80 | Core library for all sequence operations |
| `numpy` | >= 1.23 | Used by advanced_topics.py for SVD in Kabsch structural alignment |
| `pandas` | >= 1.5 (optional) | Used by different_parameters.py for tabular output |

Install all at once (see Setup below).

### External tools (all optional)

These are third-party alignment programs. The code detects them
automatically via `shutil.which()`. If none are found, MSA.py falls
back to its built-in progressive aligner — no error is raised.

| Tool | Download | Used in |
|------|----------|---------|
| ClustalW2 | http://www.clustal.org/clustal2/ | MSA.py, different_parameters.py |
| MUSCLE | https://drive5.com/muscle/ | MSA.py, different_parameters.py |
| MAFFT | https://mafft.cbrc.jp/alignment/software/ | MSA.py, different_parameters.py |

If installed, add the executable to your system PATH so the scripts
can find it automatically.

---

## Setup and Installation

### Step 1 — Clone the repository

```bash
git clone <repository-url>
cd bio
```

### Step 2 — Create a virtual environment

```bash
# Windows
python -m venv .venv
.venv\Scripts\activate

# macOS / Linux
python3 -m venv .venv
source .venv/bin/activate
```

### Step 3 — Install Python dependencies

```bash
pip install biopython numpy pandas
```

To verify Biopython installed correctly:

```bash
python -c "import Bio; print(Bio.__version__)"
```

### Step 4 — (Optional) Install external MSA tools

Only needed if you want ClustalW, MUSCLE, or MAFFT comparisons in MSA.py.
Without them the scripts still run using the built-in fallback aligner.

**Windows — MAFFT example:**
Download the Windows installer from https://mafft.cbrc.jp/alignment/software/
and add the install directory to your PATH environment variable.

**Linux (Ubuntu/Debian):**
```bash
sudo apt-get install clustalw muscle mafft
```

**macOS (Homebrew):**
```bash
brew install clustal-w muscle mafft
```

---

## File-by-File Guide

### `Data_Retrieval.py` — Task 2: Data Collection and Preprocessing

Connects to NCBI Entrez and UniProt (ExPASy) at runtime to download
real biological sequences.

**What it does:**
- Fetches the human beta-actin mRNA record (`NM_001101.5`) from NCBI
  in FASTA format using `Bio.Entrez.efetch`.
- Fetches the human hemoglobin beta protein (`P68871`) from UniProt
  via `Bio.ExPASy.get_sprot_raw` and `Bio.SwissProt.read`.
- Implements `preprocess_sequence()` which cleans a SeqRecord:
  strips non-standard characters from DNA sequences (keeps only A/T/C/G),
  and passes protein sequences through unchanged.
- Prints the first 50 bases of the processed actin sequence.

**Note:** requires an active internet connection. The email address
in `Entrez.email` is set as required by NCBI's usage policy.

---

### `global_alignment_demo.py` — Task 3a: Needleman-Wunsch Global Alignment

Implements and demonstrates the Needleman-Wunsch algorithm through
Biopython's `Bio.Align.PairwiseAligner` in `mode="global"`.

**What it does:**
- `global_alignment_demo(seq1, seq2, seq_type)` — configures an aligner
  for the given sequence type:
  - DNA: `match=2.0`, `mismatch=-1.0`, `open=-2.0`, `extend=-0.5`
  - Protein: BLOSUM62 substitution matrix, `open=-10.0`, `extend=-0.5`
- Prints the alignment score, the formatted optimal alignment, and
  the percent sequence identity.
- `calculate_identity(alignment)` — counts matching non-gap pairs
  and returns the identity as a percentage.
- Runs a demo on two short human/mouse actin sequences (`ATGCGTACGTTAGC`
  vs `ATGCGTTCGTTAGC`).

---

### `Local_Alignment .py` — Task 3a: Smith-Waterman Local Alignment

Implements local alignment using `Bio.Align.PairwiseAligner` in
`mode="local"` (Smith-Waterman).

**What it does:**
- `local_alignment_demo(seq1, seq2, seq_type)` — same configuration
  options as the global aligner but mode is local.
- Finds the highest-scoring local alignment between two sequences,
  which may cover only a sub-region of each sequence.
- Prints the aligned region, its score, and the start/end coordinates
  of the aligned region within each original sequence.
- Runs a demo on human and mouse hemoglobin alpha N-terminal peptides
  (`MVLSPADKTNVKAAWGKVG` vs `MVLSGEDKSNVKAAWGKVG`) to find the
  conserved core.

**Key difference from global:** local alignment skips poorly matching
ends and finds the best matching sub-region — useful for finding
conserved domains in otherwise divergent proteins.

---

### `Pairwise_sequence.py` — Task 3b: BLOSUM62 Substitution Matrix

A focused demonstration of configuring a pairwise aligner with the
BLOSUM62 substitution matrix.

**What it does:**
- Creates a `PairwiseAligner` and loads `BLOSUM62` from
  `Bio.Align.substitution_matrices`.
- Applies the matrix alongside explicit gap penalties:
  `open=-2.0`, `extend=-0.5`.
- Prints all aligner settings to show exactly how the substitution
  matrix interacts with the rest of the scoring model.

**Why BLOSUM62?** It is derived from conserved blocks in known protein
alignments and is the standard default for protein similarity searches
(e.g. BLAST). It scores conservative substitutions (e.g. I/V) higher
than radical ones (e.g. W/G).

---

### `different_parameters.py` — Task 3c: Gap Penalty Experiments

Systematically tests how different gap penalty values affect alignment
quality, and includes wrapper functions for external MSA tools.

**What it does:**

*Part 1 — Gap penalty experiment:*
- `calculate_identity(alignment)` — computes percent identity over the
  full alignment length (including gaps in the denominator).
- `parameter_experiment(seq1, seq2, seq_type)` — iterates over a grid
  of 3 gap-open penalties (`-5`, `-10`, `-15`) x 3 gap-extend penalties
  (`-0.5`, `-1.0`, `-2.0`) = 9 combinations.
- For each combination it reports: alignment score, percent identity,
  and total gap count.
- Output is printed as a pandas DataFrame if pandas is installed,
  or as plain dicts if not.
- When run directly, imports `human_hbb_seq` and `mouse_hbb_seq` from
  `Data_Retrieval.py` (falls back to short stubs if unavailable).

*Part 2 — MSA tool wrappers (legacy / reference):*
- `run_clustalw()`, `run_muscle()`, `run_mafft()` — wrapper functions
  using `Bio.Align.Applications` command-line wrappers (older API).
  These are reference implementations kept for comparison; MSA.py
  uses updated subprocess-based wrappers.

---

### `MSA.py` — Task 4: Multiple Sequence Alignment

The main MSA module. Tries external tools in order, falls back to a
built-in aligner, then visualizes and reports statistics.

**Functions:**

`run_clustalw(input_file, output_file)` — calls `clustalw2` via
subprocess, reads the `.aln` output using `AlignIO.read(..., "clustal")`.

`run_muscle(input_file, output_file)` — calls `muscle -in ... -out ...`,
reads the FASTA output.

`run_mafft(input_file, output_file)` — calls `mafft --auto ...`,
captures stdout to file, reads FASTA output.

`compare_msa_methods(input_fasta)` — detects installed tools with
`shutil.which()`, runs each one, records runtime and alignment length.
Returns a results dict for all three tools.

`fallback_msa(input_fasta, output_file)` — a pure-Python progressive
center-star aligner. It aligns every sequence against the first
(reference) sequence using `Bio.Align.PairwiseAligner` (with BLOSUM62
for proteins, match/mismatch scores for DNA), then merges the resulting
column-aligned sequences into a single MSA. Does not require any
external executable. Saves result to `fallback_aln.fasta`.

`visualize_msa(alignment, label)` — prints the alignment in block
format with a position ruler and a conservation line:
`*` = all sequences identical at this column,
`.` = no gaps but residues vary,
`(blank)` = one or more sequences have a gap.

`compute_msa_stats(alignment)` — returns a dict with: number of
sequences, alignment length, count of fully conserved columns,
count of gap-containing columns, and percent conserved.

**`__main__` block flow:**
1. Reads `hemoglobin.fasta`.
2. Calls `compare_msa_methods()` — tries ClustalW, MUSCLE, MAFFT.
3. If no external tool succeeded, calls `fallback_msa()`.
4. Prints a comparison summary table (tool, status, runtime, length).
5. For each successful alignment, prints the visualization and statistics.

---

### `advanced_topics.py` — Task 5: Advanced Alignment Topics

Covers three advanced sub-tasks in one standalone script. Requires no
external tools or internet connection.

#### Part A — Profile-based Alignment (PSSM / HMM)

`build_pssm(alignment, pseudocount=0.5)` — constructs a
Position-Specific Scoring Matrix from a `MultipleSeqAlignment`.
For each column it counts residue frequencies (with pseudocounts to
avoid log(0)), then computes log-odds scores against a uniform
background of 1/20 per amino acid. Returns a list of dicts, one
per alignment column.

`score_sequence_against_profile(sequence, pssm)` — sums the PSSM
log-odds score for every residue in a query sequence. This is the
exact scoring model used by PSI-BLAST and HMMER profile search.
High scores indicate the query fits the conserved pattern; low
scores indicate divergence.

`print_pssm(pssm, max_cols)` — prints a 20-row (amino acid) x
N-column table of log-odds scores.

`demo_profile_alignment()` — builds a 5-species hemoglobin alpha
PSSM (Human, Mouse, Rabbit, Pig, Dog), prints it, then scores
two queries: a gorilla sequence (closely related, high score ~29.8)
and a random protein (dissimilar, low score ~-2.0). Explains the
link to PSI-BLAST and profile HMMs (HMMER).

#### Part B — Structural Alignment (3-D Coordinates)

`compute_rmsd(coords1, coords2)` — computes the Root Mean Square
Deviation between two equal-length sets of 3-D coordinates:
`RMSD = sqrt(mean(sum|r1-r2|^2))`.

`kabsch_superimpose(P, Q)` — implements the Kabsch algorithm in pure
Python (with numpy SVD if available, centroid-only fallback without).
Centres both sets on their centroids, builds the covariance matrix H,
applies SVD to find the optimal rotation matrix R, rotates P onto Q,
and returns the rotated coordinates along with the post-superimposition
RMSD.

`build_mock_atoms(coords)` — creates lightweight atom-like objects
that satisfy Bio.PDB Superimposer's interface, allowing the demo to
call `Bio.PDB.Superimposer.set_atoms()` with plain coordinate data.

`demo_structural_alignment()` — uses two simulated C-alpha helix
fragments (8 atoms each), computes RMSD before superimposition (~0.40 A),
then uses Bio.PDB Superimposer if available (falls back to Kabsch),
and reports RMSD after (~0.03 A).

#### Part C — Consensus Sequence Generation

`compute_conservation(alignment)` — iterates over every column of a
MultipleSeqAlignment and returns a list of tuples:
`(consensus_residue, frequency, shannon_entropy)`. Gaps are excluded
from frequency calculations. Entropy of 0 means fully conserved;
higher entropy means more variation.

`generate_consensus(alignment, threshold=0.5)` — produces a consensus
string. Positions where the most frequent residue exceeds the threshold
get that residue; others get `?` (ambiguous).

`print_conservation_table(alignment, conservation)` — tabular report
showing column index, consensus residue, frequency, Shannon entropy
in bits, and a three-tier conservation label:
`***` highly conserved (>=80%), `**` moderate (50-80%), `*` variable (<50%).

`demo_consensus_generation()` — aligns 5 full-length hemoglobin alpha
sequences, tries Biopython's `SummaryInfo.dumb_consensus` first, then
always computes the manual consensus. Prints the first 20 columns of
the conservation table and a summary showing ~82% of positions are
highly conserved across species.

---

## Data Files

### `hemoglobin.fasta`

Five full-length hemoglobin alpha-chain protein sequences in FASTA format.
Used as the primary input for `MSA.py`.

| ID | Species | Length |
|----|---------|--------|
| Human_HBA | Homo sapiens | 141 aa |
| Mouse_HBA | Mus musculus | 141 aa |
| Rabbit_HBA | Oryctolagus cuniculus | 141 aa |
| Dog_HBA | Canis lupus familiaris | 141 aa |
| Horse_HBA | Equus caballus | 141 aa |

### `NC_005816.fna` / `NC_005816.gb`

Complete genomic sequence and GenBank annotation for *Yersinia pestis*
plasmid pPCP1. Used in tutorial notebooks to demonstrate `SeqIO.read`
with GenBank and FASTA formats.

### `ls_orchid.fasta` / `ls_orchid.gbk`

A set of orchid 18S rRNA sequences, a classic Biopython tutorial dataset.
Used in `tutorial1.ipynb` and `04-seqio.ipynb` to demonstrate parsing
and iterating over multi-record files.

### `mafft_aln.fasta`

A pre-generated MAFFT alignment output. Can be inspected without running
MAFFT:

```python
from Bio import AlignIO
aln = AlignIO.read("mafft_aln.fasta", "fasta")
print(aln)
```

### `example.fastq`

A short FASTQ file for demonstrating quality-score parsing in
`04-seqio.ipynb`.

---

## How to Run Each Script

Make sure your virtual environment is activated first:

```bash
# Windows
.venv\Scripts\activate

# macOS / Linux
source .venv/bin/activate
```

### Task 2 — Data Retrieval (requires internet)

```bash
python Data_Retrieval.py
```

Downloads human beta-actin mRNA from NCBI and hemoglobin beta from
UniProt, then prints the sequences and their preprocessed lengths.

### Task 3a — Global Alignment (Needleman-Wunsch)

```bash
python global_alignment_demo.py
```

Aligns `ATGCGTACGTTAGC` vs `ATGCGTTCGTTAGC` and prints the score,
formatted alignment, and percent identity.

### Task 3a — Local Alignment (Smith-Waterman)

```bash
python "Local_Alignment .py"
```

Note the space in the filename — quotes are required on the command line.
Aligns two hemoglobin alpha peptides and prints the best local region.

### Task 3b — BLOSUM62 Substitution Matrix

```bash
python Pairwise_sequence.py
```

Configures and prints all aligner settings with BLOSUM62.

### Task 3c — Gap Penalty Experiments

```bash
python different_parameters.py
```

Prints a 9-row table comparing alignment score, identity, and gap count
across different gap open/extend penalty combinations.

### Task 4 — Multiple Sequence Alignment

```bash
python MSA.py
```

Attempts ClustalW, MUSCLE, MAFFT in order. Falls back to the built-in
progressive aligner if none are found. Prints a summary table,
a visualized alignment with conservation line, and per-alignment statistics.

### Task 5 — Advanced Topics

```bash
python advanced_topics.py
```

Runs all three demos in sequence:
- Part A: PSSM profile building and sequence scoring (~29.8 for gorilla, ~-2.0 for random)
- Part B: Structural superimposition with RMSD (~0.40 A before, ~0.03 A after)
- Part C: Consensus generation — ~82% positions highly conserved in hemoglobin alpha

### Jupyter Notebooks

```bash
jupyter notebook
```

Then open any `.ipynb` file from the browser interface.

| Notebook | Content |
|----------|---------|
| `tutorial1.ipynb` | Biopython basics, Seq/SeqRecord usage |
| `02-seq_objects.ipynb` | Seq methods: complement, transcription, translation |
| `03-seq_annot.ipynb` | SeqRecord features and annotations |
| `04-seqio.ipynb` | Reading/writing FASTA, GenBank, FASTQ files |
| `class_quiz.ipynb` | Practice exercises |

---

## Optional: External MSA Tools

The scripts detect external tools automatically — no configuration needed.
Just ensure the executable is on your system PATH.

To check if a tool is available:

```bash
# Windows
where clustalw2
where muscle
where mafft

# macOS / Linux
which clustalw2
which muscle
which mafft
```

When any tool is found, `MSA.py` runs it and shows its runtime and
alignment length alongside the built-in result for direct comparison.

---

## Expected Output Summary

| Script | Key output |
|--------|-----------|
| `Data_Retrieval.py` | NCBI / UniProt sequences + preprocessed length |
| `global_alignment_demo.py` | Score, formatted block alignment, % identity |
| `Local_Alignment .py` | Best local region, score, start/end coordinates |
| `Pairwise_sequence.py` | Aligner configuration printout with BLOSUM62 |
| `different_parameters.py` | 9-row table of gap penalty effects |
| `MSA.py` | Tool comparison table, block MSA with ruler, statistics |
| `advanced_topics.py` | PSSM table, profile scores, RMSD values, conservation table |

---

## Troubleshooting

**`UnicodeEncodeError` on Windows**

Set the console encoding before running:

```bash
set PYTHONIOENCODING=utf-8
python <script>.py
```

**`ModuleNotFoundError: No module named 'Bio'`**

The virtual environment is not activated, or Biopython was not installed:

```bash
.venv\Scripts\activate       # Windows
pip install biopython
```

**`Data_Retrieval.py` hangs or raises `urllib.error.URLError`**

This script fetches live data from NCBI and UniProt. Check your network
connection. NCBI may throttle repeated requests — wait a few seconds
and retry.

**`MSA.py` shows all tools as NOT FOUND**

This is expected if ClustalW, MUSCLE, and MAFFT are not installed.
The script automatically runs the built-in progressive aligner and
produces a valid alignment — no action needed.

**`Local_Alignment .py` — command not found**

The filename contains a space. Always quote it on the command line:

```bash
python "Local_Alignment .py"
```

**Biopython deprecation warning about `SummaryInfo`**

This warning from `advanced_topics.py` is harmless. The code catches
the condition and falls back to a manual implementation automatically.

---

*DCIT 411 Bioinformatics — Sequence Alignment with Biopython*
*Submission deadline: February 27, 2026*
