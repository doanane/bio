from Bio import Entrez, SeqIO
from Bio import ExPASy
from Bio import SwissProt


Entrez.email = "anane365221@gmail.com"


handle = Entrez.efetch(db="nucleotide", id="NM_001101.5", 
                       rettype="fasta", retmode="text")
human_actin = SeqIO.read(handle, "fasta")
handle.close()


handle = ExPASy.get_sprot_raw("P68871")
human_hbb_record = SwissProt.read(handle)
human_hbb_seq = human_hbb_record.sequence
handle.close()
print("Human beta-actin DNA sequence:")
print(human_actin.seq)
print("\nHuman hemoglobin beta amino acid sequence:")
print(human_hbb_seq)


def preprocess_sequence(record):
    """Clean and prepare sequence for alignment"""
    seq = str(record.seq).upper()
    
    
    if set(seq).issubset({'A', 'T', 'C', 'G'}):
        
        seq = ''.join([c for c in seq if c in 'ATCG'])
    
    
    elif set(seq).issubset(set('ACDEFGHIKLMNPQRSTVWY')):
        
        pass
    
    
    return seq


human_actin_seq = preprocess_sequence(human_actin)
print(f"Human actin length: {len(human_actin_seq)}")
print(f"First 50 bases: {human_actin_seq[:50]}")