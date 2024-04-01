import pandas as pd
import numpy as np
import glob
import argparse
from multiprocessing import Pool
from sklearn.linear_model import LinearRegression
from typing import Tuple

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cpus", type=int, required=False, default=1)
    return parser.parse_args()

def get_block_std(series: pd.Series, number_blocks: int) -> float:
    total_size = len(series)
    block_size = int(total_size / number_blocks)
    block_means = []
    
    for i in range(number_blocks):
        subset = series[i * block_size:(i + 1) * block_size]
        block_mean = subset.mean()
        block_means.append(block_mean)
    
    std_dev = np.std(block_means)
    
    return std_dev

def get_statistical_inefficiency_curve(series: pd.Series) -> Tuple[np.array]:
    total_size = len(series)
    std_total = series.std()
    n_blocks_array = np.linspace(2, int(total_size/10), 10, dtype=int)
    stat_ineff_list = []
    inverse_block_size_list = []
    for n_blocks in n_blocks_array:
        block_size = int(total_size/n_blocks)
        std = get_block_std(series, n_blocks)
        stat_ineff = block_size*(std**2)/(std_total**2)
        inverse_block_size_list.append(n_blocks/total_size)
        stat_ineff_list.append(stat_ineff)
    return np.array(inverse_block_size_list).reshape(-1,1), np.array(stat_ineff_list).reshape(-1,1)

def get_series_statistics(series: pd.Series) -> pd.DataFrame:
    if series.nunique() > 2:
        x, y = get_statistical_inefficiency_curve(series)
        regression = LinearRegression().fit(x, y)
        estimated_statistical_inefficiency = regression.intercept_[0]
        total_steps = len(series)
        std_global = series.std()
        std_corrected = np.sqrt((std_global**2) * estimated_statistical_inefficiency / total_steps)
    else:
        std_corrected = 0
        estimated_statistical_inefficiency = "NA"
    
    mean = series.mean()

    results = pd.DataFrame({
            "Property":series.name,
            "Mean":mean,
            "STD":std_corrected,
            "Statistical inefficiency":estimated_statistical_inefficiency
        },
        index=[0]
    )
    
    return results

def get_summary(production_data):
    
    properties = [
        "temperature",
        "volume",
        "number_molecules", 
        "density_number", 
        "density_mass", 
        "density_mol", 
        "potential_total", 
        "potential_bonded", 
        "potential_nonbonded", 
        "potential_lrc", 
        "potential_walls", 
        "potential_perturbed",
        "potential_perturbed_squared",
        "weight_ghost_molecule", 
        "pressure_total", 
        "pressure_excess", 
        "pressure_ideal", 
        "pressure_lrc", 
        "weight_ideal_chain",
        "x_size",
        "y_size", 
        "z_size"
    ]
    
    production_data = production_data.apply(pd.to_numeric, errors="coerce")
    dfs = []
    for column in properties:
        if column in production_data.columns:
            dfs.append(get_series_statistics(production_data[column]))

    summary_df = pd.concat(dfs)

    return summary_df

def get_potential_perturbed_covariance(series: pd.Series) -> float:
    series_squared = series**2
    global_cov = series.cov(series_squared)
    num_points = len(series)
    n_blocks_array = np.linspace(2, int(num_points/10), 10, dtype=int)
    stat_ineff_list = []
    inverse_block_size_list = []
    for n_blocks in n_blocks_array:
        block_size = int(num_points/n_blocks)
        list_a = []
        list_b = []
        for i in range(n_blocks):
            list_a.append(series[i*block_size:(i+1)*block_size].mean())
            list_b.append(series_squared[i*block_size:(i+1)*block_size].mean())
        cov = pd.Series(list_a).cov(pd.Series(list_b))
        stat_ineff = block_size*(cov**2)/(global_cov**2)
        inverse_block_size_list.append(n_blocks/num_points)
        stat_ineff_list.append(stat_ineff)

    x = np.array(inverse_block_size_list).reshape(-1, 1)
    y = np.array(stat_ineff_list).reshape(-1, 1)

    model = LinearRegression().fit(x, y)

    final_statistical_inefficienty = model.intercept_[0]

    cov_corrected = final_statistical_inefficienty*global_cov/num_points
    return cov_corrected

def get_properties_summary(filename):
    system_name = filename.split("properties")[0]

    raw_data = pd.read_csv(filename, sep="\s+")
    production_data = (
        raw_data
        .apply(pd.to_numeric, errors="coerce")
        .query("production == 1")
        .reset_index(drop=True)
    )

    production_data["potential_perturbed_squared"] = (
        production_data
        .potential_perturbed
        .map(lambda x: x**2)
    )

    results_summary = get_summary(production_data).set_index("Property")

    beta = 1/(BOLTZMANN_CONSTANT*results_summary.loc["temperature", "Mean"])
    n = results_summary.loc["number_molecules", "Mean"]
    a1 = 0
    a2 = 0
    a1_std = 0
    a2_std = 0

    mean_u1 = results_summary.loc["potential_perturbed", "Mean"]
    mean_u1_sqr = results_summary.loc["potential_perturbed_squared", "Mean"]
    std_u1 = results_summary.loc["potential_perturbed", "STD"]
    std_u1_sqr = results_summary.loc["potential_perturbed_squared", "STD"]
    cov_u1_u1_sqr = get_potential_perturbed_covariance(production_data.potential_perturbed)

    a1 = (beta/n)*mean_u1
    a2 = -0.5*((beta**2)/n)*(mean_u1_sqr - (mean_u1**2))

    a1_std = (beta/n)*std_u1
    a2_std = 0.5*((beta**2)/n)*np.sqrt(
        (std_u1_sqr**2) 
        + 4*((mean_u1**2)*(std_u1**2))
        - 4*(mean_u1*cov_u1_u1_sqr)
    )

    perturbation_results = {
        "a1": {"Mean": a1, "STD": a1_std, "Statistical inefficiency": "NA"},
        "a2": {"Mean": a2, "STD": a2_std, "Statistical inefficiency": "NA"}
    }

    displacement_acceptance = production_data["displacement_acceptance"].iloc[-1]

    displacement_results = {
        "displacement_acceptance": {"Mean": displacement_acceptance, "STD": 0, "Statistical inefficiency": "NA"}
    }

    results_concat = pd.concat([
        results_summary,
        pd.DataFrame.from_dict(perturbation_results, orient="index"),
        pd.DataFrame.from_dict(displacement_results, orient="index")
    ])

    results_concat.to_csv(system_name+"summary.csv")

if __name__ == "__main__":
    BOLTZMANN_CONSTANT = 1.38064852E-23
    CPUS = int(argparser().cpus)
    filename_list = glob.glob("*properties*")
    with Pool(processes=CPUS) as pool:
        pool.map(get_properties_summary, filename_list)