import pandas as pd
import numpy as np
import argparse
import glob
from typing import List
from multiprocessing import Pool

def parse_args():
    """
    Parses the command line arguments using ArgumentParser.
    
    Returns:
        args (Namespace): An object containing the parsed command line arguments.
    """
    parser = argparse.ArgumentParser(description="Argument Parser for Optimize Code")
    
    parser.add_argument('--number_bins', type=int, help='Number of bins')
    parser.add_argument('--cpus', type=int, required=False, default=1, help='Number of CPUs (optional, default=1)')
    
    return parser.parse_args()

def get_bin(z: float, bin_centers: np.ndarray) -> float:
    """
    Finds the closest bin center to a given value.

    Args:
        z (float): The value for which to find the closest bin center.
        bin_centers (ndarray): An array of bin centers.

    Returns:
        float: The bin center that is closest to the given value.
    """
    distances = np.abs(bin_centers - z)
    return bin_centers[np.argmin(distances)]

def make_bin_centers(box_height: float, number_of_bins: int) -> np.ndarray:
    """
    Calculate the bin centers for a given box height and number of bins.
    
    Parameters:
        box_height (float): The height of the box.
        number_of_bins (int): The number of bins.
        
    Returns:
        numpy.ndarray: An array of bin centers.
    """
    half_box_height = 0.5 * box_height
    bin_size = box_height / number_of_bins
    
    z_min = -half_box_height
    z_max = half_box_height
    
    start = z_min + 0.5 * bin_size
    end = z_max - 0.5 * bin_size
    
    bin_centers = np.linspace(start, end, number_of_bins)
    
    return bin_centers

def get_number_molecules_profile_results(results_list: List[pd.DataFrame], bin_centers: np.ndarray) -> pd.DataFrame:
    """
    Calculate the number of molecules for each bin in the given results list.

    Args:
        results_list (list): A list of pandas DataFrames containing the results.
        bin_centers (array-like): An array or list of bin centers.

    Returns:
        pandas.DataFrame: The number of molecules for each bin.
    """

    df_number_molecules = (
        pd.concat(
            [result[["number_molecules"]].rename(columns={"number_molecules":i}) for i, result in enumerate(results_list)],
            axis=1
        )
        .fillna(0)
        .T
    )

    return df_number_molecules


def get_orietantion_profile_results(results_list: List[pd.DataFrame], bin_centers: np.ndarray) -> pd.DataFrame:
    """
    Calculate the orientation profile results.

    Parameters:
    - results_list: A list of dataframes containing the results.
    - bin_centers: An array or list of bin centers.

    Returns:
    - The transposed orientation profile dataframe.
    """

    df_orientation = (
        pd.concat(
            [result[["S"]].rename(columns={"S":i}) for i, result in enumerate(results_list)],
            axis=1
        )
        .fillna(0)
        .T
    )
    
    return df_orientation

def get_profile_summary(profiles):
    """
    Calculates the mean, standard deviation, and number of unique values for a given dataframe of profiles.

    Parameters:
    - profiles (pandas.DataFrame): A dataframe containing profiles data

    Returns:
    - statistics (pandas.DataFrame): A dataframe with columns 'index', 'mean', 'std', and 'nunique'
                                     representing the statistical summary of the profiles data
    """
    midpoint = profiles.shape[0] // 2
    profiles_production = profiles.iloc[midpoint:]
    
    statistics = (
        profiles_production
        .agg(['mean', 'std', 'nunique'])
        .T
        .reset_index(names="Bin")
    )

    statistics_mirror = (
        statistics
        .copy()
        .assign(
            Bin = lambda df: (
                -1*df["Bin"]
            )
        )
    )

    statistics_final = (
        pd.concat([statistics, statistics_mirror], ignore_index=True)
        .sort_values("Bin")
    )
    
    return statistics_final

def get_profiles(filename: str) -> None:
    """
    Processes the given file and generates various profile summaries.

    Parameters:
        filename (str): The name of the file to process.

    Returns:
        None

    Raises:
        Any exceptions that occur during the processing.

    """
    try:
        base_name = filename.split('raw_profile')[0]
        
        with open(filename) as raw_profile:
            aux_df = pd.DataFrame(columns=['z', 'S'])
            results_list = []
            index_counter = 0
            bin_centers = None
            box_z_size = None
            
            for line_number, line in enumerate(raw_profile):
                splited_line = line.split()

                if line_number == 0:
                    if len(splited_line) == 4:
                        box_x_size, box_y_size, box_z_size, molar_mass = map(float, splited_line)
                    else:    
                        box_x_size, box_y_size, box_z_size, _, molar_mass = map(float, splited_line)
                    bin_centers = make_bin_centers(box_z_size, NUMBER_BINS)
                    bin_size = box_z_size / NUMBER_BINS
                    bin_volume = bin_size * box_x_size * box_y_size * 1.0E-30

                if len(splited_line) >= 4 :
                    if len(aux_df) > 0:
                        results_inter = (
                            aux_df
                            .apply(lambda x: x.astype(float))
                            .assign(
                                bin = lambda df: (
                                    df["z"]
                                    .apply(lambda z: get_bin(z, bin_centers))
                                )
                            )
                        )
                        
                        results = (
                            results_inter
                            .groupby("bin", as_index=False)
                            .agg(
                                S = ("S", "mean"),
                                number_molecules = ("z", "count")
                            )
                            .assign(
                                abs_bin = lambda df: (
                                    np.abs(df["bin"])
                                )
                            )
                            .groupby("abs_bin")
                            .agg(
                                S = ("S", "mean"),
                                number_molecules = ("number_molecules", "mean")
                            )
                        )
                        
                        results_list.append(results.copy())
                        aux_df = pd.DataFrame(columns=['z', 'S'])
                        index_counter = 0
                elif len(splited_line) == 2:
                    aux_df.loc[index_counter] = splited_line
                    index_counter += 1
        
        number_molecules_profiles = get_number_molecules_profile_results(results_list, bin_centers)
        orietantion_profiles = get_orietantion_profile_results(results_list, bin_centers)

        number_molecules_summary = get_profile_summary(number_molecules_profiles)

        molar_density_summary = number_molecules_summary.copy()
        molar_density_summary['mean'] /= bin_volume * AVOGADRO_NUMBER

        mass_density_summary = molar_density_summary.copy()
        mass_density_summary['mean'] *= molar_mass / 1000

        orietantion_summary = get_profile_summary(orietantion_profiles)

        number_molecules_summary.to_csv(base_name + 'number_molecules_profile.csv', index=False)
        molar_density_summary.to_csv(base_name + 'molar_density_profile.csv', index=False)
        mass_density_summary.to_csv(base_name + 'mass_density_profile.csv', index=False)
        orietantion_summary.to_csv(base_name + 'orientation_profile.csv', index=False)
    
    except Exception as e:
        print("An error occurred:", str(e))


if __name__ == '__main__':
    NUMBER_BINS = int(parse_args().number_bins)
    CPUS = int(parse_args().cpus)
    AVOGADRO_NUMBER = 6.023E23
    filename_list = glob.glob('*raw_profile*')

    with Pool(processes=CPUS) as pool:
        pool.map(get_profiles, filename_list)
