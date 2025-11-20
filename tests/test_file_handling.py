from proteosim.file_handling import read_fasta

def test_read_fasta():
    tmp_fasta_path = '../data/dummy_proteins.fasta'
    protein_map = read_fasta(tmp_fasta_path)
    
    assert protein_map['P11802'] == "MAT"
    assert protein_map['A0A087WTH1'] == "MVI"