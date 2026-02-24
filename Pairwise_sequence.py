from Bio import Align
from Bio.Align import substitution_matrices


aligner = Align.PairwiseAligner()


aligner.mode = "global"  
aligner.match_score = 2.0
aligner.mismatch_score = -1.0
aligner.open_gap_score = -2.0
aligner.extend_gap_score = -0.5


blosum62 = substitution_matrices.load("BLOSUM62")
aligner.substitution_matrix = blosum62
print("Aligner configured with BLOSUM62 substitution matrix.") 
print(f"Match score: {aligner.match_score}")
print(f"Mismatch score: {aligner.mismatch_score}")
print(f"Gap open score: {aligner.open_gap_score}")
print(f"Gap extend score: {aligner.extend_gap_score}")
print(f"Alignment mode: {aligner.mode}")