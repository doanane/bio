from Bio import Align
from Bio.Align import substitution_matrices


def global_alignment_demo(seq1, seq2, seq_type="dna"):
    """
    Perform global alignment using Needleman-Wunsch algorithm
    """
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    
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
    print(f"Alignment score: {score}")
    
    
    alignments = aligner.align(seq1, seq2)
    
    
    alignment = alignments[0]
    print("\nOptimal global alignment:")
    print(alignment)
    
    
    identity = calculate_identity(alignment)
    print(f"Sequence identity: {identity:.2f}%")
    
    return alignment

def calculate_identity(alignment):
    """Calculate percent identity from alignment"""
    seq1_aligned = alignment[0]
    seq2_aligned = alignment[1]
    
    matches = 0
    total = 0
    
    for a, b in zip(seq1_aligned, seq2_aligned):
        if a != '-' and b != '-':  
            total += 1
            if a == b:
                matches += 1
    
    if total == 0:
        return 0.0
    return (matches / total) * 100


human_actin_seq = "ATGCGTACGTTAGC"
mouse_actin_seq = "ATGCGTTCGTTAGC"

alignment = global_alignment_demo(human_actin_seq, mouse_actin_seq, "dna")