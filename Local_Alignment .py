from Bio import Align
from Bio.Align import substitution_matrices


def local_alignment_demo(seq1, seq2, seq_type="protein"):
    """
    Perform local alignment using Smith-Waterman algorithm
    """
    aligner = Align.PairwiseAligner()
    aligner.mode = "local"
    
    if seq_type == "protein":
        aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
        aligner.open_gap_score = -10.0
        aligner.extend_gap_score = -0.5
    else:
        aligner.match_score = 2.0
        aligner.mismatch_score = -1.0
        aligner.open_gap_score = -2.0
        aligner.extend_gap_score = -0.5
    
    score = aligner.score(seq1, seq2)
    print(f"Local alignment score: {score}")
    
    alignments = aligner.align(seq1, seq2)
    
    if len(alignments) > 0:
        alignment = alignments[0]
        print("\nBest local alignment:")
        print(alignment)
        
        # Show coordinates of aligned region
        coordinates = alignment.coordinates
        print(f"\nAligned region coordinates:")
        print(f"Sequence A: {coordinates[0][0]}-{coordinates[0][-1]}")
        print(f"Sequence B: {coordinates[1][0]}-{coordinates[1][-1]}")
        
        return alignment
    else:
        print("No local alignment found")
        return None

# Test with human and mouse hemoglobin alpha (should find conserved regions)
# Example sequences for demo (define these to avoid undefined variable errors)
# These are short illustrative peptide sequences — replace with real sequences as needed
human_hba_seq = "MVLSPADKTNVKAAWGKVG"
mouse_hba_seq = "MVLSGEDKSNVKAAWGKVG"

alignment = local_alignment_demo(human_hba_seq, mouse_hba_seq, "protein")