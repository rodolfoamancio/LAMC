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
                    # inicialização
                    if len(splited_line) == 4:
                        box_x_size, box_y_size, box_z_size, molar_mass = map(float, splited_line)
                    else:    
                        box_x_size, box_y_size, box_z_size, _, molar_mass = map(float, splited_line)
                    bin_centers = make_bin_centers(box_z_size, NUMBER_BINS)
                    bin_size = box_z_size / NUMBER_BINS
                    bin_volume = bin_size * box_x_size * box_y_size
                elif len(splited_line) == 2:
                    # adição novo dado
                    aux_df.loc[index_counter] = splited_line
                    index_counter += 1
                elif len(splited_line) == 0:
                    # consolidação de resultados
                    results_inter = (
                        aux_df
                        .apply(lambda x: x.astype(float))
                        .assign(
                            bin = lambda d: d.z.apply(lambda x: get_bin(x, bin_centers))
                        )
                        .groupby("bin", as_index=False)
                        .agg(
                            s = ("S", "mean"),
                            n = ("z", "count")
                        )
                        .merge(
                            pd.DataFrame({"bin":bin_centers}),
                            on="bin",
                            how="right"
                        )
                        .fillna(0)
                        .assign(
                            abs_bin = lambda d: np.round(d.bin.abs(), 4)
                        )
                        .groupby("abs_bin", as_index=False)
                        .agg(
                            s = ("s", "mean"),
                            n = ("n", "mean")
                        )
                    )

                    results_list.append(results_inter.copy())
                    aux_df = pd.DataFrame(columns=['z', 'S'])
                    index_counter = 0
                    
        
        production_list = results_list[int(len(results_list)*0.5):]

        N = 5
        size = len(production_list)//N

        chunks = [pd.concat(production_list[size*i:size*(i+1) if i < N - 1 else None]) for i in range(N)]

        chunks_avg = [
            (
                chunk
                .groupby("abs_bin", as_index=False)
                .agg(
                    s=("s", "mean"),
                    n=("n", "mean")
                )
            )
            for chunk in chunks
        ]

        df_summary = (
            pd.concat(chunks_avg, ignore_index=True)
            .groupby("abs_bin", as_index=False)
            .agg(
                s_avg = ("s", "mean"),
                s_std = ("s", "std"),
                n_avg = ("n", "mean"),
                n_std = ("n", "std")
            )
        )

        df_mirror = (
            df_summary
            .query("abs_bin != 0")
            .assign(abs_bin = lambda d: -1*d.abs_bin)
        )

        df_complete = (
            pd.concat([df_summary, df_mirror], ignore_index=True)
            .rename(columns={"abs_bin":"bin"})
            .assign(
                # molecules/A³
                d_n_avg = lambda d: d.n_avg/bin_volume,
                d_n_std = lambda d: d.n_std/bin_volume,
                # mol/m³
                d_mol_avg = lambda d: d.d_n_avg*1E30/AVOGADRO_NUMBER,
                d_mol_std = lambda d: d.d_n_std*1E30/AVOGADRO_NUMBER,
                # kg/m³
                d_mass_avg = lambda d: d.d_mol_avg*molar_mass/1000,
                d_mass_std = lambda d: d.d_mol_std*molar_mass/1000,
            )
        )

        number_molecules_summary = df_complete[["bin", "d_n_avg", "d_n_std"]].rename(columns={"d_n_avg":"avg", "d_n_std":"std"})
        molar_density_summary = df_complete[["bin", "d_mol_avg", "d_mol_std"]].rename(columns={"d_mol_avg":"avg", "d_mol_std":"std"})
        mass_density_summary = df_complete[["bin", "d_mass_avg", "d_mass_std"]].rename(columns={"d_mass_avg":"avg", "d_mass_std":"std"})
        orietantion_summary = df_complete[["bin", "s_avg", "s_std"]].rename(columns={"s_avg":"avg", "s_std":"std"})

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
