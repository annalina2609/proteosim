from proteosim.mass_spectra_simulation import calculate_mol_mass
def test_calculate_mol_mass():
    aa_mass_dict = {
    'A': 71.08, 'R': 156.19, 'N': 114.10, 'D': 115.09,
    'C': 103.15, 'E': 129.12, 'Q': 128.13, 'G': 57.05,
    'H': 137.14, 'I': 113.16, 'L': 113.16, 'K': 128.17,
    'M': 131.19, 'F': 147.18, 'P': 97.12, 'S': 87.08,
    'T': 101.11, 'W': 186.21, 'Y': 163.18, 'V': 99.13,
    }
    test_peptide = 'ACH'
    expected_mass = {'peptide': 311.37 }
    calculated_mass = calculate_mol_mass(test_peptide, aa_mass_dict)
    assert calculated_mass == expected_mass

from proteosim.mass_spectra_simulation import calculate_mol_mass_collection
def test_calculate_mol_mass_collection():
    aa_mass_dict ={
    'A': 71.08, 'R': 156.19, 'N': 114.10, 'D': 115.09,
    'C': 103.15, 'E': 129.12, 'Q': 128.13, 'G': 57.05,
    'H': 137.14, 'I': 113.16, 'L': 113.16, 'K': 128.17,
    'M': 131.19, 'F': 147.18, 'P': 97.12, 'S': 87.08,
    'T': 101.11, 'W': 186.21, 'Y': 163.18, 'V': 99.13,
    }
    peptides = ['AA', 'PEP']
    expected = {
        'AA': 142.16,
        'PEP': 323.36
        }
    actual = calculate_mol_mass_collection(peptides,aa_mass_dict)

    assert actual == expected

from proteosim.mass_spectra_simulation import calculate_mz_collection
def test_calculate_mz_collection():
    peptide_mass_map = {
        'AA': 142.16,
        'PP': 194.24}
    actual = calculate_mz_collection(peptide_mass_map, charge=2)
    expected = {
        'AA': 72.087,
        'PP': 98.127}
    assert actual == expected

from proteosim.mass_spectra_simulation import fragment_peptide
def test_fragment_peptide():
    peptide = 'ABC'
    expected = ['A', 'AB', 'C', 'BC']
    actual = fragment_peptide(peptide)

    assert actual == expected