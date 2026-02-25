"""
DCIT 411: Bioinformatics -- Task 5: Advanced Topics
===================================================
Covers three sub-tasks:
  5a. Profile-based sequence alignment (PSSM / Hidden Markov Model concepts)
  5b. Structural alignment (aligning protein structures by 3-D coordinates)
  5c. Consensus sequence generation (most conserved residues across an MSA)
"""

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os


# PART A: Profile-based Sequence Alignment (PSSM / HMM)


AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")


def build_pssm(alignment, pseudocount=0.5):
    """
    Build a Position-Specific Scoring Matrix (PSSM) from a MultipleSeqAlignment.

    Each column records log-odds scores for every amino acid relative to the
    background uniform frequency (1/20).  Pseudocounts avoid log(0).

    Returns
    -------
    pssm : list of dicts  (length = alignment columns)
    """
    length = alignment.get_alignment_length()
    n_seqs = len(alignment)
    background = 1.0 / len(AMINO_ACIDS)
    pssm = []

    for col in range(length):
        column = [str(rec.seq)[col].upper() for rec in alignment]
        counts = {aa: pseudocount for aa in AMINO_ACIDS}
        for residue in column:
            if residue in counts:
                counts[residue] += 1

        total = sum(counts.values())
        scores = {}
        import math
        for aa in AMINO_ACIDS:
            freq = counts[aa] / total
            scores[aa] = round(math.log2(freq / background), 3)
        pssm.append(scores)

    return pssm


def score_sequence_against_profile(sequence, pssm):
    """
    Score a sequence against a PSSM profile.

    Slides the sequence (or uses it directly if same length) and returns
    the total log-odds score -- the core idea behind profile-sequence
    alignment used in PSI-BLAST and HMM profile search.
    """
    seq = sequence.upper()
    if len(seq) != len(pssm):
        raise ValueError(
            f"Sequence length ({len(seq)}) must equal profile length ({len(pssm)}). "
            "Pad/trim as required."
        )

    total_score = 0.0
    for i, residue in enumerate(seq):
        if residue in pssm[i]:
            total_score += pssm[i][residue]
        else:
            total_score += -4.0  # gap / unknown penalty
    return total_score


def print_pssm(pssm, max_cols=10):
    """Print a human-readable slice of the PSSM."""
    cols = min(len(pssm), max_cols)
    header = "AA  |" + "".join(f" Col{i+1:>3}" for i in range(cols))
    print(header)
    print("-" * len(header))
    for aa in AMINO_ACIDS:
        row = f" {aa}  |"
        for i in range(cols):
            row += f" {pssm[i][aa]:>6.2f}"
        print(row)


def demo_profile_alignment():
    """
    Demonstrate profile-based alignment:
      1. Build a PSSM from a small hemoglobin-like MSA.
      2. Score two query sequences against the profile.
      3. Explain the link to PSI-BLAST and profile HMMs.
    """
    print("=" * 60)
    print("PART A: Profile-based Alignment (PSSM / HMM)")
    print("=" * 60)

    # Build a toy 5-sequence alignment of hemoglobin alpha-chain N-terminus
    seqs = [
        "MVLSPADKTNVK",   # Human HBA
        "MVLSGEDKSNVK",   # Mouse HBA
        "MVLSAADKSNVK",   # Rabbit HBA
        "MVLSAADKANVK",   # Pig HBA
        "MVLSGDEKSNVK",   # Dog HBA
    ]
    ids = ["Human_HBA", "Mouse_HBA", "Rabbit_HBA", "Pig_HBA", "Dog_HBA"]

    records = [SeqRecord(Seq(s), id=i, description="") for s, i in zip(seqs, ids)]
    alignment = MultipleSeqAlignment(records)

    print("\nInput MSA (hemoglobin alpha N-terminus):")
    for rec in alignment:
        print(f"  {rec.id:<12} {rec.seq}")

    pssm = build_pssm(alignment)

    print(f"\nPSSM built from {len(alignment)} sequences, {len(pssm)} positions.")
    print("(Showing first 10 columns -- log2 odds vs. uniform background)\n")
    print_pssm(pssm, max_cols=len(pssm))

    # Score two query sequences
    queries = {
        "Gorilla_HBA":  "MVLSPADKTNVK",   # very similar to human
        "Random_prot":  "ACWYNFSDTLQK",   # dissimilar
    }

    print("\nProfile-sequence scoring (PSI-BLAST-style log-odds sum):")
    for name, seq in queries.items():
        score = score_sequence_against_profile(seq, pssm)
        print(f"  {name:<15}: {score:>7.2f}")

    print("""
Interpretation:
  - A high score means the query matches the conserved positions in the profile.
  - PSI-BLAST iteratively builds a PSSM from BLAST hits and rescans the database.
  - Profile HMMs (e.g. HMMER) extend this with probabilistic emission /
    transition scores for insertions and deletions, giving a full HMM.
""")



# PART B: Structural Alignment (aligning 3-D protein coordinates)


def build_mock_atoms(coords):
    """
    Return a list of simple objects that mimic Bio.PDB Atom CA coordinates,
    used when Bio.PDB structures are unavailable.
    """
    class MockAtom:
        def __init__(self, xyz):
            self.coord = xyz
        def get_coord(self):
            return self.coord

    return [MockAtom(c) for c in coords]


def compute_rmsd(coords1, coords2):
    """
    Compute RMSD between two equal-length lists of 3-D coordinates.
    RMSD = sqrt( mean( sum_i |r1_i - r2_i|^2 ) )
    """
    import math
    if len(coords1) != len(coords2):
        raise ValueError("Coordinate lists must have the same length.")
    n = len(coords1)
    total = 0.0
    for c1, c2 in zip(coords1, coords2):
        total += sum((a - b) ** 2 for a, b in zip(c1, c2))
    return math.sqrt(total / n)


def kabsch_superimpose(P, Q):
    """
    Kabsch algorithm: rotate P onto Q to minimise RMSD.

    Parameters
    ----------
    P, Q : list of (x, y, z) tuples / lists -- must be same length.

    Returns
    -------
    P_rotated : rotated coordinates of P
    rmsd_after : RMSD after superimposition
    """
    import math

    def centroid(coords):
        n = len(coords)
        return [sum(c[i] for c in coords) / n for i in range(3)]

    def subtract(coords, centre):
        return [[c[i] - centre[i] for i in range(3)] for c in coords]

    def matmul(A, B):
        rows_A, cols_A = len(A), len(A[0])
        cols_B = len(B[0])
        C = [[0.0] * cols_B for _ in range(rows_A)]
        for i in range(rows_A):
            for j in range(cols_B):
                for k in range(cols_A):
                    C[i][j] += A[i][k] * B[k][j]
        return C

    def transpose(M):
        return [[M[j][i] for j in range(len(M))] for i in range(len(M[0]))]

    # Centre both sets
    cp = centroid(P)
    cq = centroid(Q)
    Pc = subtract(P, cp)
    Qc = subtract(Q, cq)

    # Covariance matrix H = Pc^T · Qc  (3x3)
    H = [[0.0] * 3 for _ in range(3)]
    for p, q in zip(Pc, Qc):
        for i in range(3):
            for j in range(3):
                H[i][j] += p[i] * q[j]

    # SVD via numpy if available, else fall back to simplified rotation
    try:
        import numpy as np
        U, S, Vt = np.linalg.svd(H)
        d = np.linalg.det(Vt.T @ U.T)
        D = np.diag([1, 1, d])
        R = Vt.T @ D @ U.T
        P_rot = [(R @ np.array(p) + np.array(cq)).tolist() for p in Pc]
    except ImportError:
        # Without numpy just translate to common centroid (no rotation)
        P_rot = [[c + cq[i] for i, c in enumerate(p)] for p in Pc]

    rmsd_after = compute_rmsd(P_rot, Q)
    return P_rot, rmsd_after


def demo_structural_alignment():
    """
    Demonstrate structural alignment:
      1. Try to use Bio.PDB Superimposer on real C-alpha coordinates.
      2. Fall back to a built-in Kabsch implementation if Bio.PDB is absent.
      3. Report RMSD before and after superimposition.
    """
    print("=" * 60)
    print("PART B: Structural Alignment (3-D Coordinate Superimposition)")
    print("=" * 60)

    # Two sets of C-alpha coordinates (A) -- simplified helix fragment
    # Structure 1: human myoglobin helix-A region (simulated)
    struct1_ca = [
        (1.0,  0.0,  0.0),
        (2.4,  1.2,  0.8),
        (3.5,  0.6,  2.1),
        (4.8,  1.8,  1.5),
        (5.2,  3.0,  2.8),
        (6.5,  2.5,  4.0),
        (7.1,  3.8,  3.3),
        (8.4,  3.2,  4.6),
    ]
    # Structure 2: sperm whale myoglobin -- slightly rotated + translated
    struct2_ca = [
        (1.2,  0.3,  0.1),
        (2.6,  1.5,  0.9),
        (3.7,  0.9,  2.2),
        (5.0,  2.1,  1.7),
        (5.4,  3.3,  3.0),
        (6.7,  2.8,  4.2),
        (7.3,  4.1,  3.5),
        (8.6,  3.5,  4.8),
    ]

    rmsd_before = compute_rmsd(struct1_ca, struct2_ca)
    print(f"\nRMSD before superimposition : {rmsd_before:.4f} A")

    # Try Bio.PDB Superimposer
    try:
        from Bio.PDB import Superimposer
        atoms1 = build_mock_atoms(struct1_ca)
        atoms2 = build_mock_atoms(struct2_ca)
        sup = Superimposer()
        sup.set_atoms(atoms1, atoms2)
        rmsd_after = sup.rms
        print(f"RMSD after superimposition  : {rmsd_after:.4f} A  (Bio.PDB Superimposer)")
    except Exception:
        _, rmsd_after = kabsch_superimpose(list(struct1_ca), list(struct2_ca))
        print(f"RMSD after superimposition  : {rmsd_after:.4f} A  (Kabsch algorithm)")

    print(f"\nImprovement : {rmsd_before - rmsd_after:.4f} A")
    print("""
Key Concepts:
  - Structural alignment finds the optimal rotation/translation that minimises
    RMSD between equivalent C-alpha atoms of two protein structures.
  - The Kabsch algorithm (used internally by Bio.PDB Superimposer) solves this
    via Singular Value Decomposition (SVD).
  - RMSD < 2 A typically indicates structural similarity; < 1 A is very high.
  - Structural alignment can reveal evolutionary relationships even when
    sequence identity is too low for sequence-based alignment (<20-30%).
""")



# PART C: Consensus Sequence Generation


def compute_conservation(alignment):
    """
    For each column return:
      - consensus residue (most frequent non-gap character)
      - frequency of that residue (0-1)
      - entropy (bits) -- 0 = fully conserved, log2(N) = fully random
    """
    import math
    length = alignment.get_alignment_length()
    conservation = []

    for col in range(length):
        column = [str(rec.seq)[col].upper() for rec in alignment]
        residues = [c for c in column if c != '-']
        if not residues:
            conservation.append(('-', 0.0, 0.0))
            continue

        counts = {}
        for r in residues:
            counts[r] = counts.get(r, 0) + 1

        total = len(residues)
        consensus_res = max(counts, key=counts.get)
        consensus_freq = counts[consensus_res] / total

        # Shannon entropy
        entropy = 0.0
        for c, n in counts.items():
            p = n / total
            if p > 0:
                entropy -= p * math.log2(p)

        conservation.append((consensus_res, round(consensus_freq, 3), round(entropy, 3)))

    return conservation


def generate_consensus(alignment, threshold=0.5):
    """
    Generate a consensus sequence.
      - threshold : minimum frequency for a residue to appear in consensus;
                    below threshold a '?' (ambiguous) is placed.
    """
    conservation = compute_conservation(alignment)
    consensus = ""
    for res, freq, _ in conservation:
        if freq >= threshold:
            consensus += res
        else:
            consensus += "?"
    return consensus


def print_conservation_table(alignment, conservation):
    """Print a coloured-text conservation report."""
    print("\nConservation analysis per column:")
    header = f"{'Col':>4}  {'Consensus':>9}  {'Freq':>6}  {'Entropy(bits)':>13}  {'Conservation'}"
    print(header)
    print("-" * len(header))
    for i, (res, freq, entropy) in enumerate(conservation, 1):
        if freq >= 0.8:
            level = "*** Highly conserved"
        elif freq >= 0.5:
            level = "**  Moderately conserved"
        else:
            level = "*   Variable"
        print(f"{i:>4}  {res:>9}  {freq:>6.1%}  {entropy:>13.3f}  {level}")


def demo_consensus_generation():
    """
    Demonstrate consensus sequence generation from a MultipleSeqAlignment.
    """
    print("=" * 60)
    print("PART C: Consensus Sequence Generation")
    print("=" * 60)

    # Aligned hemoglobin alpha sequences (same length, already aligned)
    aligned_seqs = [
        "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR",
        "MVLSGEDKSNVKAAWGKVGGHAGEYGAEALERMFLGFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVGHLDDLPGALSALSDLHAYKLRVDPVNFKLLSHCLLVTLAAHHPADFTPAVHASLDKFLSSVSTVLTSKYR",
        "MVLSAADKTNVKGVFSKIGGHAAEYGAEALERMFLGFPTTKTYFPHFDLSHGSAQVKAHGKKVADALASAAGHLDDLPQALSALSDLHAYKLRVDPVNFKLLSHCLLSTLAVHLPNDFTPAVHASLDKFMASVSTVLTSKYR",
        "MVLSAADKANVKTIFSKIGGHAAEYGAEALERMFASHPTTKTYFPHFDLSHGSAQVKAHGKKVADALANAAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPADFTPAVHASLDKFLANVSTVLTSKYR-",
        "MVLSGEDKSNVKAAWGKVGGHAGEYGAEALERMFLGFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVGHLDDLPGALSALSDLHAYKLRVDPVNFKLLSHCLLVTLAAHHPADFTPAVHASLDKFLSSVSTVLTSKYR",
    ]
    ids = ["Human", "Mouse", "Rabbit", "Pig", "Dog"]

    # Ensure uniform length (pad shorter sequences)
    max_len = max(len(s) for s in aligned_seqs)
    aligned_seqs = [s.ljust(max_len, '-') for s in aligned_seqs]

    records = [SeqRecord(Seq(s), id=i, description="") for s, i in zip(aligned_seqs, ids)]
    alignment = MultipleSeqAlignment(records)

    print(f"\nAlignment: {len(alignment)} sequences x {alignment.get_alignment_length()} columns")
    print("\nSequence IDs:")
    for rec in alignment:
        print(f"  {rec.id}")

    # Attempt to use Biopython's built-in SummaryInfo if available
    bio_consensus = None
    try:
        from Bio.Align import AlignInfo
        summary = AlignInfo.SummaryInfo(alignment)
        bio_consensus = str(summary.dumb_consensus(threshold=0.5, ambiguous="?"))
        print(f"\nBiopython dumb_consensus (threshold=0.5):")
        print(f"  {bio_consensus[:60]}...")
    except Exception as e:
        print(f"\n(Bio.Align.AlignInfo unavailable: {e} -- using manual implementation)")

    # Always compute with our manual method
    conservation = compute_conservation(alignment)
    manual_consensus = generate_consensus(alignment, threshold=0.5)

    print(f"\nManual consensus sequence (threshold=0.5):")
    print(f"  {manual_consensus[:60]}...")
    print(f"\nConsensus length: {len(manual_consensus)} residues")

    # Conservation report for first 20 columns
    short_conservation = conservation[:20]
    print("\n--- Conservation table (first 20 columns) ---")
    print_conservation_table(alignment, short_conservation)

    # Summary statistics
    highly_conserved = sum(1 for _, freq, _ in conservation if freq >= 0.8)
    variable = sum(1 for _, freq, _ in conservation if freq < 0.5)
    print(f"\nSummary over all {len(conservation)} positions:")
    print(f"  Highly conserved (>=80%): {highly_conserved} ({highly_conserved/len(conservation):.1%})")
    print(f"  Variable (<50%)         : {variable} ({variable/len(conservation):.1%})")
    print("""
Interpretation:
  - Highly conserved positions often correspond to functionally critical residues
    (e.g. heme-binding residues in hemoglobin).
  - Consensus sequences are used to build sequence logos, identify motifs,
    and annotate protein families in databases such as Pfam and PROSITE.
""")



# MAIN


if __name__ == "__main__":
    demo_profile_alignment()
    print()
    demo_structural_alignment()
    print()
    demo_consensus_generation()
