def read_fasta(filepath):
    """
    Reads fasta files and creates a dictionary called protein_map, which contains all protein IDs as keys and their corresponding amino acid sequences 

    Parameters
    ----------
    filepath: str

    Returns
    -------
    protein_map: dict
        Dictionary mapping protein ids to protein amino acid sequences
    """
    fasta_path = filepath
    
    protein_map = {}
    current_id = None
    current_sequence = []
    
    with open(fasta_path, 'r', encoding='utf-8') as fasta_handle:
        for line in fasta_handle:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith('>'):
                if current_id is not None:
                    protein_map[current_id] = ''.join(current_sequence)
                    current_sequence = []
                current_id = stripped.split('|')[1]
            else:
                    current_sequence.append(stripped)
    
    if current_id is not None:
        protein_map[current_id] = ''.join(current_sequence)
    
    return protein_map
    

