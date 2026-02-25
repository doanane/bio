from Bio import Align, SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def progressive_msa(input_fasta, output_file="fallback_aln.fasta"):
    """
    Progressive multiple sequence alignment using center-star method.
    Aligns all sequences to the first sequence and merges progressively.
    """
    # Read sequences
    records = list(SeqIO.parse(input_fasta, "fasta"))
    if len(records) == 0:
        raise RuntimeError("No sequences found")
    
    # Handle single sequence case
    if len(records) == 1:
        msa = MultipleSeqAlignment([records[0]])
        return msa
    
    # Use first sequence as reference (center)
    reference = str(records[0].seq)
    aligned_sequences = [reference]  # List of aligned strings
    sequence_ids = [records[0].id]
    
    # Create aligner
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.open_gap_score = -10.0
    aligner.extend_gap_score = -0.5
    
    # Determine sequence type for appropriate scoring
    if is_protein_sequence(reference):
        try:
            aligner.substitution_matrix = Align.substitution_matrices.load("BLOSUM62")
        except:
            # Fallback to simple scoring
            aligner.match_score = 2.0
            aligner.mismatch_score = -1.0
    else:
        aligner.match_score = 2.0
        aligner.mismatch_score = -1.0
    
    # Align each sequence to reference
    pairwise_alignments = []
    for i, record in enumerate(records[1:], 1):
        query = str(record.seq)
        alignments = aligner.align(reference, query)
        best_alignment = alignments[0]
        
        pairwise_alignments.append({
            'id': record.id,
            'ref_aligned': str(best_alignment[0]),
            'query_aligned': str(best_alignment[1])
        })
    
    # Merge all alignments progressively
    final_alignment = merge_alignments(aligned_sequences[0], 
                                      [p['ref_aligned'] for p in pairwise_alignments],
                                      [p['query_aligned'] for p in pairwise_alignments])
    
    # Create final MSA object
    seq_records = []
    for i, seq_str in enumerate(final_alignment['sequences']):
        seq_records.append(SeqRecord(Seq(seq_str), 
                                     id=sequence_ids[i] if i < len(sequence_ids) else pairwise_alignments[i-1]['id']))
    
    msa = MultipleSeqAlignment(seq_records)
    return msa

def is_protein_sequence(sequence):
    """Detect if sequence is protein based on alphabet"""
    protein_letters = set('ACDEFGHIKLMNPQRSTVWY')
    seq_letters = set(sequence.upper())
    return seq_letters.issubset(protein_letters)

def merge_alignments(reference_seq, ref_alignments, query_alignments):
    """
    Merge multiple pairwise alignments into a single MSA.
    
    Args:
        reference_seq: Original reference sequence (unaligned)
        ref_alignments: List of reference sequences after pairwise alignment
        query_alignments: List of query sequences after pairwise alignment
    """
    # Start with first pairwise alignment
    aligned_seqs = [ref_alignments[0], query_alignments[0]]
    
    # Progressively add remaining sequences
    for i in range(1, len(ref_alignments)):
        new_ref = ref_alignments[i]
        new_query = query_alignments[i]
        
        # Merge existing alignment with new pairwise alignment
        merged = []
        for existing_seq in aligned_seqs:
            # Insert gaps to match new_ref's gap pattern
            merged_seq = insert_gaps(existing_seq, new_ref)
            merged.append(merged_seq)
        
        # Add the new query sequence with gaps matching new_ref
        merged.append(new_query)
        aligned_seqs = merged
    
    return {'sequences': aligned_seqs}

def insert_gaps(sequence, gap_pattern):
    """
    Insert gaps into sequence to match the gap pattern of gap_pattern.
    This ensures all sequences have the same length in the final MSA.
    """
    result = []
    pattern_idx = 0
    seq_idx = 0
    
    while pattern_idx < len(gap_pattern):
        if gap_pattern[pattern_idx] == '-':
            result.append('-')
        else:
            if seq_idx < len(sequence):
                result.append(sequence[seq_idx])
                seq_idx += 1
            else:
                result.append('-')
        pattern_idx += 1
    
    # Add any remaining sequence characters
    while seq_idx < len(sequence):
        result.append(sequence[seq_idx])
        seq_idx += 1
    
    return ''.join(result)