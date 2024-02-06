from hashlib import new
from math import prod
from types import new_class
import pandas as pd
import numpy as np
import glob
import argparse
from multiprocessing import Pool
from sklearn.linear_model import LinearRegression

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cpus', type=int, required=False, default=1)
    return parser.parse_args()

def get_block_std(series, number_blocks):
    total_size = len(series)
    block_size = int(total_size / number_blocks)
    block_means = []
    
    for i in range(number_blocks):
        subset = series[i * block_size:(i + 1) * block_size]
        block_mean = subset.mean()
        block_means.append(block_mean)
    
    std_dev = np.std(block_means)
    
    return std_dev

def get_statistical_inefficiency_curve(series):
    total_size = len(series)
    std_total = series.std()
    n_blocks_array = np.linspace(2, int(total_size/10), 50, dtype=int)
    std_list = []
    stat_ineff_list = []
    inverse_block_size_list = []
    for n_blocks in n_blocks_array:
        block_size = int(total_size/n_blocks)
        std = get_block_std(series, n_blocks)
        stat_ineff = block_size*(std**2)/(std_total**2)
        inverse_block_size_list.append(n_blocks/total_size)
        std_list.append(std)
        stat_ineff_list.append(stat_ineff)
    return np.array(inverse_block_size_list).reshape(-1,1), np.array(stat_ineff_list).reshape(-1,1)

def get_series_statistics(series):
    if series.nunique() > 2:
        x, y = get_statistical_inefficiency_curve(series)
        regression = LinearRegression().fit(x, y)
        estimated_statistical_inefficiency = regression.intercept_[0]
        total_steps = len(series)
        std_global = series.std()
        std_corrected = np.sqrt((std_global**2) * estimated_statistical_inefficiency / total_steps)
    else:
        std_corrected = 0
        estimated_statistical_inefficiency = 'NA'
    
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
        'temperature',
        'volume',
        'number_molecules', 
        'density_number', 
        'density_mass', 
        'density_mol', 
        'potential_total', 
        'potential_bonded', 
        'potential_nonbonded', 
        'potential_lrc', 
        'potential_walls', 
        'potential_perturbed',
        'weight_ghost_molecule', 
        'pressure_total', 
        'pressure_excess', 
        'pressure_ideal', 
        'pressure_lrc', 
        'weight_ideal_chain',
        'x_size',
        'y_size', 
        'z_size'
    ]
    
    production_data = production_data.apply(pd.to_numeric, errors='coerce')
    dfs = []
    for column in properties:
        if column in production_data.columns:
            dfs.append(get_series_statistics(production_data[column]))

    summary_df = pd.concat(dfs)

    return summary_df

def get_properties_summary(filename):
    try:
        system_name = filename.split('properties')[0]

        raw_data = pd.read_csv(filename, sep='\s+')
        production_data = raw_data.loc[int(raw_data.shape[0]/2):].reset_index(drop=True)
        

        object_columns = production_data.select_dtypes(include=['object']).columns
        production_data[object_columns] = production_data[object_columns].apply(pd.to_numeric, errors='coerce')

        
        results_summary = get_summary(production_data).set_index('Property')

        if 'weight_ghost_molecule' in production_data.columns:
            mean_weight_ghost_molecule = results_summary.loc['weight_ghost_molecule', 'Mean']
            std_weight_ghost_molecule = results_summary.loc['weight_ghost_molecule', 'STD']
            mean_weight_ideal_chain = production_data['weight_ideal_chain'].mean()
            std_weight_ideal_chain = production_data['std_weight_ideal_chain'].mean()
            pressure_ideal_gas = results_summary.loc['pressure_ideal', 'Mean']

            std_fugacity = np.sqrt(
                (pressure_ideal_gas**2) *
                ((std_weight_ideal_chain**2) / (mean_weight_ghost_molecule**2) +
                 (mean_weight_ideal_chain**2) * (std_weight_ghost_molecule**2) / (mean_weight_ghost_molecule**4))
            )
            fugacity = pressure_ideal_gas * mean_weight_ideal_chain / mean_weight_ghost_molecule
            mean_pressure = results_summary.loc['pressure_total', 'Mean']
            std_pressure = results_summary.loc['pressure_total', 'STD']
            std_phi = np.sqrt(
                (std_fugacity**2) / (mean_pressure**2) +
                (fugacity**2) * (std_pressure**2) / (mean_pressure**4)
            )
            phi = fugacity / results_summary.loc['pressure_total', 'Mean']

            fugacity_results = {
                'phi': {'Mean': phi, 'STD': std_phi, 'Statistical inefficiency': 'NA'},
                'fugacity': {'Mean': fugacity, 'STD': std_fugacity, 'Statistical inefficiency': 'NA'}
            }
            results_summary = pd.concat([results_summary, pd.DataFrame.from_dict(fugacity_results, orient='index')])

        beta = 0
        a1 = 0
        a2 = 0

        try:
            beta = 1 / (BOLTZMANN_CONSTANT * production_data["temperature"].mean())
            a1 = (production_data["potential_perturbed"].mean() * beta / production_data["number_molecules"].mean())
            a2 = (-0.5*production_data["potential_perturbed"].var(ddof=0) * (beta**2) / production_data["number_molecules"].mean())
        except:
            pass

        perturbation_results = {
            'a1': {'Mean': a1, 'STD': 0, 'Statistical inefficiency': 'NA'},
            'a2': {'Mean': a2, 'STD': 0, 'Statistical inefficiency': 'NA'}
        }

        displacement_acceptance = 0
        try:
            displacement_acceptance = production_data["displacement_acceptance"].iloc[-1]
        except:
            pass
            
        displacement_results = {
            'displacement_acceptance': {'Mean': displacement_acceptance, 'STD': 0, 'Statistical inefficiency': 'NA'}
        }

        results_concat = pd.concat([
            results_summary,
            pd.DataFrame.from_dict(perturbation_results, orient='index'),
            pd.DataFrame.from_dict(displacement_results, orient='index')
        ])

        results_concat.to_csv(system_name+'summary.csv')

    except Exception as e:
        print(f"An error occurred: {str(e)}")

if __name__ == "__main__":
    BOLTZMANN_CONSTANT = 1.38064852E-23
    CPUS = int(argparser().cpus)
    filename_list = glob.glob('*properties*')
    with Pool(processes=CPUS) as pool:
        pool.map(get_properties_summary, filename_list)