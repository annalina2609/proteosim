# Enzyme cleavage patterns is a variable containing a sample of cleavage patterns

import re

enzyme_cleavage_patterns = {
    'LysC': r'(?<=K)',
    'LysN': r'(?=K)',
    'ArgC': r'(?<=R)',
    'Trypsin': r'(?<=[KR])(?!P)',
}

enzyme_cleavage_patterns

# Function for digestion of a single protein sequence

def digest_protein_sequence(protein_seq, cleave_pattern, min_pep_len=5, max_pep_len=30):
    """
    Digest a protein sequence using a given enzymatic cleavage pattern (regex).

    The protein sequence is split into peptide fragments according to the
    specified cleavage pattern. The resulting peptides are filtered so that
    only fragments with a length suitable for detection by mass spectrometry
    (MS) are returned.

    Parameters
    ----------
    protein_seq = str
        Amino acid sequence of a protein as a string.
    cleave_pattern = regex pattern
        Regular expression describing the enzyme-specific cleavage rule.
        You might chose from the `enzyme_cleavage_patterns` dictionary.
    min_pep_len = int
        Minimum peptide length to keep (default is 5).
    max_pep_len = int
        Maximum peptide length to keep (default is 30).

    Returns
    -------
    filtered_peptides = list of str
        list of peptides as strings (only those peptide fragments that are detactable)
    """
    peptides = re.split(cleave_pattern, protein_seq)

    filtered_peptides = [p for p in peptides if min_pep_len <= len(p) <= max_pep_len]
    
    return filtered_peptides

# Function for digestion of a collection of protein sequences

def digest_protein_collection(protein_map, cleave_pattern, min_pep_len=5, max_pep_len=30):
    """
    Digest a collection of protein sequences using a given enzymatic cleavage pattern (regex).

    The protein sequences are split into peptide fragments according to the
    specified cleavage pattern. The resulting peptides are filtered so that
    only fragments with a length suitable for detection by mass spectrometry
    (MS) are returned.

    Parameters
    ----------
    protein_map = dict
        dictionary of proteins and their according amino acide sequence
    cleave_pattern = regex pattern
        Regular expression describing the enzyme-specific cleavage rule.
        You might chose from the `enzyme_cleavage_patterns` dictionary.
    min_pep_len = int
        Minimum peptide length to keep (default is 5).
    max_pep_len = int
        Maximum peptide length to keep (default is 30).

    Returns
    -------
    digested_proteins = dict
        dictionary containing protein IDs as keys and corresponding lists of digested peptides as strings (only those peptide fragments that are detactable)
    """
    digested_proteins = {}

    for protein_id, sequence in protein_map.items():
        
        protein_fragemnts = re.split(cleave_pattern, sequence)
        filtered_fragments = [p for p in protein_fragemnts if min_pep_len <= len(p) <= max_pep_len]

        digested_proteins[protein_id] = filtered_fragments
    
    return digested_proteins

# Function to compute sequence coverage 

def compute_sequence_coverage(protein_seq, peptides):
    """
    Function reports coverage percentage for any protein/peptide combination. 

    Parameters
    ----------
    protein_seq : str
        amino acide sequence of protein as string
    peptides: list
        list of peptide sequences as strings
    Returns
    -------
    covered_percentage: float
        percentage of protein sequence covered by peptide fragments 
    """
    covered = [False] * len(protein_seq)

    for pep in peptides: 
        start = 0
        while True: 
            idx = protein_seq.find(pep, start)
            if idx == -1:
                break
            for i in range(idx, idx+len(pep)):
                covered[i] = True
            start = start + 1
        
    covered_percentage = sum (covered) / len(protein_seq) * 100

    return covered_percentage

