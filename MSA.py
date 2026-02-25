from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import Align
from Bio import SeqIO
import subprocess
import time
import os
import shutil


def run_clustalw(input_file, output_file="clustalw_aln.aln"):
    cmd = ["clustalw2", f"-INFILE={input_file}", f"-OUTFILE={output_file}"]
    start_time = time.time()
    proc = subprocess.run(cmd, capture_output=True, text=True)
    runtime = time.time() - start_time

    if proc.returncode != 0:
        raise RuntimeError(f"clustalw2 failed: {proc.stderr}")

    alignment = AlignIO.read(output_file, "clustal")
    return alignment, runtime


def run_muscle(input_file, output_file="muscle_aln.fasta"):
    cmd = ["muscle", "-in", input_file, "-out", output_file]
    start_time = time.time()
    proc = subprocess.run(cmd, capture_output=True, text=True)
    runtime = time.time() - start_time

    if proc.returncode != 0:
        raise RuntimeError(f"muscle failed: {proc.stderr}")

    alignment = AlignIO.read(output_file, "fasta")
    return alignment, runtime


def run_mafft(input_file, output_file="mafft_aln.fasta"):
    cmd = ["mafft", "--auto", input_file]
    start_time = time.time()
    with open(output_file, 'w') as out_f:
        proc = subprocess.run(cmd, stdout=out_f, stderr=subprocess.PIPE, text=True)
    runtime = time.time() - start_time

    if proc.returncode != 0:
        raise RuntimeError(f"mafft failed: {proc.stderr}")

    alignment = AlignIO.read(output_file, "fasta")
    return alignment, runtime


def compare_msa_methods(input_fasta):
    """Compare performance of different MSA tools with basic error handling."""
    results = {}

    # Detect available executables on PATH
    available = {
        'ClustalW': shutil.which('clustalw2') or shutil.which('clustalw'),
        'MUSCLE': shutil.which('muscle'),
        'MAFFT': shutil.which('mafft')
    }

    print("Tool availability:")
    for tool, path in available.items():
        print(f" - {tool}: {'found at ' + path if path else 'NOT FOUND'}")

    # ClustalW
    if available['ClustalW']:
        try:
            print("Running ClustalW...")
            clustalw_aln, clustalw_time = run_clustalw(input_fasta)
            results['ClustalW'] = {
                'alignment': clustalw_aln,
                'runtime': clustalw_time,
                'length': clustalw_aln.get_alignment_length()
            }
            print(f"  Runtime: {clustalw_time:.2f}s")
            print(f"  Alignment length: {clustalw_aln.get_alignment_length()}")
        except Exception as e:
            results['ClustalW'] = {'error': str(e)}
            print("  ClustalW failed:", e)
    else:
        results['ClustalW'] = {'error': 'executable not found'}
        print("Skipping ClustalW (not installed)")

    # MUSCLE
    if available['MUSCLE']:
        try:
            print("Running MUSCLE...")
            muscle_aln, muscle_time = run_muscle(input_fasta)
            results['MUSCLE'] = {
                'alignment': muscle_aln,
                'runtime': muscle_time,
                'length': muscle_aln.get_alignment_length()
            }
            print(f"  Runtime: {muscle_time:.2f}s")
            print(f"  Alignment length: {muscle_aln.get_alignment_length()}")
        except Exception as e:
            results['MUSCLE'] = {'error': str(e)}
            print("  MUSCLE failed:", e)
    else:
        results['MUSCLE'] = {'error': 'executable not found'}
        print("Skipping MUSCLE (not installed)")

    # MAFFT
    if available['MAFFT']:
        try:
            print("Running MAFFT...")
            mafft_aln, mafft_time = run_mafft(input_fasta)
            results['MAFFT'] = {
                'alignment': mafft_aln,
                'runtime': mafft_time,
                'length': mafft_aln.get_alignment_length()
            }
            print(f"  Runtime: {mafft_time:.2f}s")
            print(f"  Alignment length: {mafft_aln.get_alignment_length()}")
        except Exception as e:
            results['MAFFT'] = {'error': str(e)}
            print("  MAFFT failed:", e)
    else:
        results['MAFFT'] = {'error': 'executable not found'}
        print("Skipping MAFFT (not installed)")

    return results


def fallback_msa(input_fasta, output_file="fallback_aln.fasta"):
    """Create a simple center-star MSA by aligning all sequences to the first sequence.
    This is a heuristic fallback when external tools are not available.
    """
    records = list(SeqIO.parse(input_fasta, "fasta"))
    if len(records) == 0:
        raise RuntimeError("No sequences found in input fasta")
    if len(records) == 1:
        SeqIO.write(records, output_file, "fasta")
        msa = MultipleSeqAlignment([SeqRecord(Seq(str(records[0].seq)), id=records[0].id)])
        return msa, 0.0

    ref = str(records[0].seq)
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"

    aligned_ref = None
    aligned_seqs = []  # list of strings, first will be ref
    ids = []

    for idx, rec in enumerate(records):
        seq = str(rec.seq)
        ids.append(rec.id)
        if idx == 0:
            # reference
            aligned_ref = ref
            aligned_seqs.append(ref)
            continue

        # choose substitution matrix heuristically by alphabet
        if set(seq.upper()).issubset(set('ACDEFGHIKLMNPQRSTVWY')) and set(ref.upper()).issubset(set('ACDEFGHIKLMNPQRSTVWY')):
            try:
                aligner.substitution_matrix = Align.substitution_matrices.load("BLOSUM62")
            except Exception:
                pass
        else:
            aligner.match_score = 2.0
            aligner.mismatch_score = -1.0

        alns = aligner.align(ref, seq)
        aln = alns[0]
        a_ref = str(aln[0])
        a_seq = str(aln[1])

        if len(aligned_seqs) == 1 and aligned_seqs[0] == ref:
            # first pair, initialize aligned_ref and aligned_seqs
            aligned_ref = a_ref
            aligned_seqs = [a_ref, a_seq]
        else:
            # merge current aligned_ref with new a_ref
            old_ref = aligned_ref
            old_seqs = aligned_seqs
            new_ref = a_ref
            new_seq = a_seq

            i = j = 0
            merged_ref = []
            new_aligned = ['' for _ in range(len(old_seqs) + 1)]
            while i < len(old_ref) or j < len(new_ref):
                c1 = old_ref[i] if i < len(old_ref) else None
                c2 = new_ref[j] if j < len(new_ref) else None

                if c1 is None:
                    # append remaining c2
                    merged_ref.append(c2)
                    for k, s in enumerate(old_seqs):
                        new_aligned[k] += '-'
                    new_aligned[-1] += new_seq[j]
                    j += 1
                elif c2 is None:
                    merged_ref.append(c1)
                    for k, s in enumerate(old_seqs):
                        new_aligned[k] += old_seqs[k][i]
                    new_aligned[-1] += '-'
                    i += 1
                elif c1 == c2:
                    merged_ref.append(c1)
                    for k in range(len(old_seqs)):
                        new_aligned[k] += old_seqs[k][i]
                    new_aligned[-1] += new_seq[j]
                    i += 1
                    j += 1
                elif c1 == '-':
                    merged_ref.append('-')
                    for k in range(len(old_seqs)):
                        new_aligned[k] += old_seqs[k][i]
                    new_aligned[-1] += '-'
                    i += 1
                elif c2 == '-':
                    merged_ref.append('-')
                    for k in range(len(old_seqs)):
                        new_aligned[k] += '-'
                    new_aligned[-1] += new_seq[j]
                    j += 1
                else:
                    # different non-gap chars, consume both
                    merged_ref.append(c1)
                    for k in range(len(old_seqs)):
                        new_aligned[k] += old_seqs[k][i]
                    new_aligned[-1] += new_seq[j]
                    i += 1
                    j += 1

            aligned_ref = ''.join(merged_ref)
            aligned_seqs = new_aligned

    # Build SeqRecords
    seq_records = []
    for idx, s in enumerate(aligned_seqs):
        seq_records.append(SeqRecord(Seq(s), id=ids[idx]))

    SeqIO.write(seq_records, output_file, "fasta")
    msa = MultipleSeqAlignment(seq_records)
    return msa, 0.0


if __name__ == "__main__":
    
    fasta = "hemoglobin.fasta"
    if not os.path.exists(fasta):
        
        try:
            import importlib
            importlib.import_module('MSA')
        except Exception:
            print(f"{fasta} not found. Please create it before running alignment.")

    msa_results = compare_msa_methods(fasta)

    
    print("\nMSA Summary:")
    for tool, res in msa_results.items():
        if 'error' in res:
            print(f"- {tool}: ERROR - {res['error']}")
        else:
            print(f"- {tool}: runtime={res['runtime']:.2f}s, length={res['length']}")