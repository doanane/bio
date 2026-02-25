from Bio import AlignIO
import subprocess
import time
import os

def run_clustalw(input_file, output_file="clustalw_aln.aln"):
    """
    Run ClustalW multiple sequence alignment
    """
    
    cmd = ["clustalw2", f"-INFILE={input_file}", f"-OUTFILE={output_file}"]
    start_time = time.time()
    proc = subprocess.run(cmd, capture_output=True, text=True)
    runtime = time.time() - start_time

    if proc.returncode != 0:
        raise RuntimeError(f"clustalw2 failed: {proc.stderr}")

    
    alignment = AlignIO.read(output_file, "clustal")

    return alignment, runtime

def run_muscle(input_file, output_file="muscle_aln.fasta"):
    """
    Run MUSCLE multiple sequence alignment
    """
    
    cmd = ["muscle", "-in", input_file, "-out", output_file]
    start_time = time.time()
    proc = subprocess.run(cmd, capture_output=True, text=True)
    runtime = time.time() - start_time

    if proc.returncode != 0:
        raise RuntimeError(f"muscle failed: {proc.stderr}")

    
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
    
    
    alignment = AlignIO.read(output_file, "fasta")
    
    return alignment, runtime

def write_fasta_for_msa(sequences_dict, filename):
    """
    Write sequences to FASTA file for MSA tools
    """
    with open(filename, 'w') as f:
        for seq_id, sequence in sequences_dict.items():
            f.write(f">{seq_id}\n")
            
            for i in range(0, len(sequence), 60):
                f.write(f"{sequence[i:i+60]}\n")

if __name__ == "__main__":
    
    try:
        from Data_Retrieval import human_hbb_seq, human_hba_seq
    except Exception:
        human_hbb_seq = "VLSPADKTNVKAAW"  
        human_hba_seq = "MVLSGEDKSNVKAAWGK"  

    try:
        from Local_Alignment import mouse_hba_seq
    except Exception:
        mouse_hba_seq = "MVLSGEDKSNVKAAWGK"  

    
    try:
        from Data_Retrieval import mouse_hbb_seq
    except Exception:
        mouse_hbb_seq = "VLSPADKTNVKAAW"

    
    hemoglobin_seqs = {
        "Human_HBB": human_hbb_seq,
        "Mouse_HBB": mouse_hbb_seq,
        "Human_HBA": human_hba_seq,
        "Mouse_HBA": mouse_hba_seq
    }

    out_file = "hemoglobin.fasta"
    write_fasta_for_msa(hemoglobin_seqs, out_file)

    
    print(f"Wrote {len(hemoglobin_seqs)} sequences to {out_file}")