from Bio import Align
from Bio.Align import substitution_matrices


try:
    import pandas as pd # type: ignore
except Exception:
    pd = None


def calculate_identity(alignment):
    """
    Calculate percent identity from a PairwiseAlignment object.
    Identity is computed as matches / alignment length (including gaps).
    """
    
    a = str(alignment[0])
    b = str(alignment[1])
    if len(a) != len(b):
        
        length = max(len(a), len(b))
        a = a.ljust(length, "-")
        b = b.ljust(length, "-")
    else:
        length = len(a)

    matches = 0
    for x, y in zip(a, b):
        if x == y and x != '-':
            matches += 1

    try:
        return matches / length if length > 0 else 0.0
    except Exception:
        return 0.0

def parameter_experiment(seq1, seq2, seq_type="protein"):
    """
    Test different gap penalty combinations
    """
    results = []
    
    
    gap_open_penalties = [-5, -10, -15]
    gap_extend_penalties = [-0.5, -1, -2]
    
    for open_score in gap_open_penalties:
        for extend_score in gap_extend_penalties:
            aligner = Align.PairwiseAligner()
            aligner.mode = "global"
            
            if seq_type == "protein":
                aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
            else:
                aligner.match_score = 2.0
                aligner.mismatch_score = -1.0
            
            aligner.open_gap_score = open_score
            aligner.extend_gap_score = extend_score
            
            score = aligner.score(seq1, seq2)
            alignments = aligner.align(seq1, seq2)
            
            if len(alignments) > 0:
                alignment = alignments[0]
                identity = calculate_identity(alignment)
                
                results.append({
                    'open': open_score,
                    'extend': extend_score,
                    'score': score,
                    'identity': identity,
                    'gaps': str(alignment[0]).count('-') + str(alignment[1]).count('-')
                })
    
    
    if pd is not None:
        df = pd.DataFrame(results)
        print("Effect of gap penalties on alignment:")
        print(df.to_string(index=False))
        return df
    else:
        print("Effect of gap penalties on alignment:")
        for r in results:
            print(r)
        return results

if __name__ == "__main__":
    
    try:
        from Data_Retrieval import human_hbb_seq
    except Exception:
        human_hbb_seq = (
            "VLSPADKTNVKAAW"  
        )

    try:
        from Data_Retrieval import mouse_hbb_seq
    except Exception:
        mouse_hbb_seq = (
            "VLSPADKTNVKAAW"  
        )

    param_results = parameter_experiment(human_hbb_seq[:50], mouse_hbb_seq[:50], "protein")

    
    if pd is None:
        for r in param_results.to_dict(orient='records') if hasattr(param_results, 'to_dict') else param_results:
            print(r)


from Bio import AlignIO
from Bio.Align.Applications import ClustalwCommandline, MuscleCommandline
import subprocess
import time
import os

def run_clustalw(input_file, output_file="clustalw_aln.aln"):
    """
    Run ClustalW multiple sequence alignment
    """
    cline = ClustalwCommandline("clustalw2", infile=input_file, outfile=output_file)
    
    start_time = time.time()
    stdout, stderr = cline()
    runtime = time.time() - start_time
    
    # Read alignment
    alignment = AlignIO.read(output_file, "clustal")
    
    return alignment, runtime

def run_muscle(input_file, output_file="muscle_aln.fasta"):
    """
    Run MUSCLE multiple sequence alignment
    """
    cline = MuscleCommandline("muscle", input=input_file, out=output_file)
    
    start_time = time.time()
    stdout, stderr = cline()
    runtime = time.time() - start_time
    
    # Read alignment
    alignment = AlignIO.read(output_file, "fasta")
    
    return alignment, runtime

def run_mafft(input_file, output_file="mafft_aln.fasta"):
    """
    Run MAFFT multiple sequence alignment
    """
    cmd = ["mafft", "--auto", input_file]
    
    start_time = time.time()
    with open(output_file, 'w') as out_f:
        subprocess.run(cmd, stdout=out_f, stderr=subprocess.PIPE, text=True)
    runtime = time.time() - start_time
    
    # Read alignment
    alignment = AlignIO.read(output_file, "fasta")
    
    return alignment, runtime