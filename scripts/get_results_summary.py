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

def get_series_statistics(series: pd.Series) -> pd.DataFrame:
    if series.nunique() > 2:
        series = series.dropna()
        total_size = len(series)
        total_var = np.var(series, ddof=1)

        n_block_list = np.arange(2, total_size//10, 1)
        size_list = total_size//n_block_list
        inv_size_list = 1/size_list

        s_list = []

        for i, nb in enumerate(n_block_list):
            chunks = [series.iloc[j*size_list[i]:(j+1)*size_list[i]] for j in range(nb)]
            chunks_avg = [np.mean(chunk) for chunk in chunks]
            chunks_var = np.var(chunks_avg, ddof=1)
            s_list.append(size_list[i]*chunks_var/total_var)

        x = inv_size_list.reshape(-1, 1)
        y = np.array(s_list).reshape(-1, 1)

        regression = LinearRegression().fit(x, y)
        estimated_statistical_inefficiency = regression.intercept_[0]
        std_corrected = np.sqrt(estimated_statistical_inefficiency*total_var/total_size)
    else:
        std_corrected = None
        estimated_statistical_inefficiency = None
    
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
        "pressure_total", 
        "pressure_excess", 
        "pressure_ideal", 
        "pressure_lrc", 
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
    s_u1 = results_summary.loc["potential_perturbed", "Statistical inefficiency"]
    std_u1 = results_summary.loc["potential_perturbed", "STD"]

    s_u1_sqr = results_summary.loc["potential_perturbed_squared", "Statistical inefficiency"]
    std_u1_sqr = results_summary.loc["potential_perturbed_squared", "STD"]

    total_cov = production_data.potential_perturbed.cov(production_data.potential_perturbed_squared)

    cov_u1_u1_sqr = np.sqrt(s_u1*s_u1_sqr)*total_cov/len(production_data)

    results_summary.loc["pot_pert_pert_sqr_cov", ["Mean", "STD", "Statistical inefficiency"]] = [cov_u1_u1_sqr, cov_u1_u1_sqr, None]

    a1 = (beta/n)*mean_u1
    a2 = -0.5*((beta**2)/n)*(production_data.potential_perturbed.var(ddof=0))

    a1_std = (beta/n)*std_u1

    a2_std = 0.5*((beta**2)/n)*np.sqrt(
        (std_u1_sqr**2) 
        + 4*((mean_u1**2)*(std_u1**2))
        - 4*(mean_u1*cov_u1_u1_sqr)
    )

    perturbation_results = {
        "a1": {"Mean": a1, "STD": a1_std, "Statistical inefficiency": None},
        "a2": {"Mean": a2, "STD": a2_std, "Statistical inefficiency": None}
    }

    displacement_acceptance = production_data["displacement_acceptance"].iloc[-1]

    displacement_results = {
        "displacement_acceptance": {"Mean": displacement_acceptance, "STD": 0, "Statistical inefficiency": None}
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