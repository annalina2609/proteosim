# Function to calculate molecular mass of a single peptide from amino acid sequence

def calculate_mol_mass(peptide_seq, amino_acid_mass_dict):
    """
    Function takes a single peptide sequence and 
    returns a dictionary `{peptide: mass}` by summing the amino-acid masses using 
    a given dictionary containing amino acid masses.

    Parameters
    ----------
    peptide_seq: str
        Amino acid sequence of peptide.
    amino_acid_mass_dict: dict
        Dictionary containg all amino acids as keys and their corresponding masses. 

    Returns
    -------
    mol_mass_of_peptide: dict
        Dictionary containg mass of peptide.

    """
    mol_mass_of_peptide = {}
    mol_mass = 0.0

    for a in peptide_seq:
        mol_mass += amino_acid_mass_dict[a]
    
    mol_mass_of_peptide['peptide'] = mol_mass

    return mol_mass_of_peptide

# Function to calculate molecular mass of list of peptides from amino acid sequence

def calculate_mol_mass_collection(peptides, amino_acid_mass_dict):
    """
    Function takes a list of peptide sequences and 
    returns a dictionary `{peptide: mass}` by summing the amino-acid masses using 
    a given dictionary containing amino acid masses.

    Parameters
    ----------
    peptides: list
        List of Amino acid sequence of peptides as strings.
    amino_acid_mass_dict: dict
        Dictionary containg all amino acids as keys and their corresponding masses. 

    Returns
    -------
    mol_mass_of_peptides: dict
        Dictionary mapping masses of peptides to peptide sequence.

    """
    mol_mass_of_peptides = {}

    for p in peptides: 
        peptide_seq = p

        mol_mass_of_peptide = {}
        mol_mass = 0.0

        for a in peptide_seq:
            mol_mass += amino_acid_mass_dict[a]
            mol_mass = round(mol_mass, 2)
    
        mol_mass_of_peptides[p] = mol_mass
    
    return mol_mass_of_peptides

# Function to calculate m/z from neutral mass values of peptides

def calculate_mz_collection(peptide_mass_map, charge=2, proton_mass=1.007):
    """
    Function computes m/z values from neutral mass values of peptides after hypothetical electrospray ionization.

    Parameters
    ----------
    peptide_mass_map: dict
         Dictionary mapping peptide masses to peptide sequences.
    charge: int
        Number of charges per peptide / charge state of peptides (default value = 2)
    proton_mass: float
        Mass of a proton (default = 1.007).

    Returns
    -------
    peptide_mz_map: dict
        Dictionary mapping peptide m/z values to peptide sequences.

    """
    peptide_mz_map = {}

    for pep,mass in peptide_mass_map.items():
        mz_pep = (mass + (charge * proton_mass))/ charge
        mz_pep = round(mz_pep, 3)

        peptide_mz_map[pep] = mz_pep
    
    return peptide_mz_map

# Function plotting a simulated MS Spectrum 

import numpy as np
import matplotlib.pyplot as plt

def plot_spectrum(mz_values, random_count_range=(0, 30000), seed=42, title = "MS Spectrum", width = 15):
    """
    Generate and plot a simulated MS1 spectrum with random intensities.

    For each given m/z value, a random intensity is drawn from the specified
    range and visualized as a bar plot (spectrum).

    Parameters
    ----------
    mz_values : iterable of float
        m/z values (x-axis) for which random intensities will be generated.
    random_count_range : tuple of int, optional
        (low, high) range from which random integer intensities are sampled.
        Default is (0, 30000).
    seed : int, optional
        Random seed for reproducibility of the simulated intensities.
        Default is 42.
    title: str
        Title of the Plot.
    width: int
        Width of Spectrum Lines (default is set to 15)

    Returns
    -------
    None
        The function creates a bar plot of the simulated spectrum and shows it.
    """

    mz_intensities_map = {}
    np.random.seed(seed)

    for m in mz_values:
        mz_intensities_map[m] = np.random.randint(*random_count_range)
    
    plt.bar(mz_values, mz_intensities_map.values(), color= "pink", width = width)
    plt.title(title)
    plt.xlabel("m/z values")
    plt.ylabel("Intensities")
    plt.show

# Function fragmenting a given list of peptides

def fragment_peptide(peptide):
    """
    Generate the b- and y-ion fragment sequences for a given peptide.

    For the input peptide sequence, this function constructs the b-ion series
    as N-terminal prefixes and the y-ion series as C-terminal suffixes, and
    returns the combined list of all fragment sequences.

    Parameters
    ----------
    peptide : str
        Amino acid sequence of the peptide (e.g. "PEPTIDE").

    Returns
    -------
    combined: list of str
        Combined list of b- and y-ion sequences for the peptide.
    """

    b_ions = [peptide[:a] for a in range(1, len(peptide))]
    y_ions = [peptide[-a:] for a in range(1, len(peptide))]

    combined = b_ions + y_ions

    return combined

# Dictionary mapping all amino acide masses [dalton]

amino_acid_mass_dalton = {
    'A': 71.08, 'R': 156.19, 'N': 114.10, 'D': 115.09,
    'C': 103.15, 'E': 129.12, 'Q': 128.13, 'G': 57.05,
    'H': 137.14, 'I': 113.16, 'L': 113.16, 'K': 128.17,
    'M': 131.19, 'F': 147.18, 'P': 97.12, 'S': 87.08,
    'T': 101.11, 'W': 186.21, 'Y': 163.18, 'V': 99.13,
    }