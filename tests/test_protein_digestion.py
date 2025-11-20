# Test for Digest Protein Collection

from proteosim.protein_digestion import digest_protein_collection
from proteosim.protein_digestion import enzyme_cleavage_patterns

def test_digest_protein_collection():
    dummy_proteins = {
        'DUMMY1': 'AGHIKOPIROLHGIEFTKL',
        'DUMMY2': 'AGHJOIKOPIROLHGIEFTKLIJEOP'   # 'AGHIK' & 'OLHGIEFTK'
    }

    cleave_pattern = enzyme_cleavage_patterns['Trypsin']

    test_digested_peptides_map = digest_protein_collection(
    dummy_proteins,
    cleave_pattern= cleave_pattern,
    min_pep_len=5,
    max_pep_len=30,
)
    assert test_digested_peptides_map['DUMMY1'] == ['AGHIK','OLHGIEFTK']
    assert test_digested_peptides_map['DUMMY2'] == ['AGHJOIK','OLHGIEFTK', 'LIJEOP']

# Test for Coverage Function

from proteosim.protein_digestion import compute_sequence_coverage

def test_compute_sequence_coverage():
    dummy_prot_seq = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ' #Alphabet d.h. 26 Buchstaben
    dummy_peps = ['ABCD','HIJKLMNOP'] # 13 

    coverage = compute_sequence_coverage(dummy_prot_seq, dummy_peps)
    assert coverage == 50

    