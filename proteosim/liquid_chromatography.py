# Fuction for prediction of Liquid Chromatography Retention times from AS-Sequence

from pyteomics import achrom

def predict_lc_retention_times(peptides):
    """
    This function takes a list of peptide sequences and returns a dictionary mapping 
    each peptide to its predicted retetion time by using achrom (based on experimental data from Guo et al.).

    Parameters
    ----------
    peptides: list 
        List of strings, containing peptide amino acide sequences.

    Returns
    -------
    relative_RT_peptides: dict
        Dictionary mapping each peptide to its predicted retention time.

    """
    relative_RT_peptides = {}

    for p in peptides: 
        relative_RT = achrom.calculate_RT( p, achrom.RCs_guo_ph7_0)
        relative_RT = round(float(relative_RT), 2)
        relative_RT_peptides[p] = relative_RT
    
    return relative_RT_peptides

# Function to plot LC Retention Times
import matplotlib.pyplot as plt

def plot_retention_time(retention_times, resolution=30):
    """
    Function plots a histogram of given retention times. 

    Parameters
    ----------
    retention_times: list or dict
        List or Dictionary of proteins retention times.
    resolution: int
        Number of histogram bins with default value 30.

    Returns
    -------
    hist_RT: histogram
        Histogram of retetion times

    """
    plt.figure()
    plt.hist(retention_times,resolution, color ="pink", edgecolor = "gray", linewidth = 1.5) 

    plt.title("Histogram of Retention Times")
    plt.xlabel("Relative retention times")
    plt.ylabel("Intensities")
    plt.show()

# Function filters peptide sequences for peptides within the selected retention time window 

def select_retention_time_window(peptide_rt_map, lower_ret_time, upper_ret_time):
    """
    Function scans dictionary with retention times and returns new dictionary only containg
    those peptides with retention times inbetween the given lower and upper bounds.

    Parameters
    ----------
    peptide_rt_map: dict
        Dictionary containing peptide sequences as keys and their corresponding retention times.
    lower_ret_time: float
        Lower bound for retention times accepted by MS.
    upper_ret_time: float
        Upper bound for retention times accepted by MS.

    Returns
    -------
    filtered_RT_peptides: dict
        Dictionary containing only peptides with retention times within the given boundaries. 

    """
    filtered_RT_peptides = {}

    for pep, rt in peptide_rt_map.items():
            if lower_ret_time <= rt <= upper_ret_time:
                  filtered_RT_peptides[pep] = rt
    
    return filtered_RT_peptides