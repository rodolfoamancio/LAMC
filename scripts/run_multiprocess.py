import os
import subprocess
from multiprocessing import Pool
import glob
import argparse

def run_mc_program(input_file, lamc_dir):
    """
    Runs the MC program with the given input file and LAMC directory.
    
    Args:
        input_file (str): The path of the input file.
        lamc_dir (str): The directory path where LAMC is located.
    """
    subprocess.call([lamc_dir+"/src/LAMC.exe", "./"+input_file])

if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--base_name', type=str)
    parser.add_argument('--cpus', type=int)
    args = parser.parse_args()

    # Get LAMC directory path from environment variable
    lamc_dir = os.getenv('LAMC_DIR')

    # Get list of files matching the base name pattern
    file_list = glob.glob("*" + args.base_name + "*")

    # Execute MC programs in parallel using multiprocessing.Pool
    with Pool(processes=args.cpus) as pool:
        pool.starmap(run_mc_program, [(file, lamc_dir) for file in file_list])