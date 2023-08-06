import pandas as pd
import glob
import argparse

"""
In summary, this code reads multiple CSV files that match the 'summary' pattern in the current directory. 
It extracts relevant data from each file, adds a 'filename' column, and concatenates the processed data 
into a single dataframe. Finally, it saves the dataframe as a CSV file with the specified output name.
"""

# Importing necessary libraries

def get_relevant_data(data):
    # Function to extract relevant data from the input dataframe
    # It selects columns 'Property', 'Mean' and 'STD', assigns new column names, and concatenates them
    mean_properties = (
        data
        .copy()
        [['Property', 'Mean']]
        .assign(Property=lambda x: 'mean_' + x['Property'])
        .rename(columns={'Mean':'value'})
    )
    # Selects columns 'Property' and 'Mean' from the input dataframe
    # Assigns a new column name by adding 'mean_' prefix to each value in the 'Property' column
    # Renames the 'Mean' column as 'value'
    
    std_properties = (
        data
        .copy()
        [['Property', 'STD']]
        .assign(Property=lambda x: 'std_' + x['Property'])
        .rename(columns={'STD':'value'})
    )
    # Selects columns 'Property' and 'STD' from the input dataframe
    # Assigns a new column name by adding 'std_' prefix to each value in the 'Property' column
    # Renames the 'STD' column as 'value'
    
    full_properties = pd.concat([mean_properties, std_properties])
    # Concatenates the mean_properties and std_properties dataframes vertically
    
    full_properties.set_index('Property', inplace=True)
    # Sets the 'Property' column as the index of the dataframe
    
    return full_properties.T
    # Transposes the dataframe and returns it

def main():
    # Main function to process the data
    
    output_name = parse_args().output_name
    # Parses the command line arguments and retrieves the value for 'output_name'
    
    filename_list = glob.glob('*summary*')
    # Retrieves a list of file names matching the '*summary*' pattern using glob library
    
    dfs = []
    # List to store the processed dataframes for each file
    
    for filename in filename_list:
        raw_data = pd.read_csv(filename, header=None, names=["Property", "Mean", "STD", "Statistical inefficiency"], skiprows=[0])
        selected_data = get_relevant_data(raw_data)
        # Calls the get_relevant_data function to extract relevant data from the dataframe
        
        selected_data['filename'] = filename
        # Adds a new column 'filename' and assigns the current filename as its value
        
        dfs.append(selected_data)
        # Appends the processed dataframe to the list of dfs
            
    
    new_df = pd.concat(dfs, ignore_index=True)
    # Concatenates all the processed dataframes into a single dataframe
    
    new_df.to_csv(output_name)
    # Saves the dataframe as a CSV file with the specified output name


def parse_args():
    # Function to parse command line arguments
    
    parser = argparse.ArgumentParser()
    # Creates an argument parser object
    
    parser.add_argument('--output_name', type=str, help='Name for output file')
    # Defines a command line argument '--output_name' of type string with a help message
    
    return parser.parse_args()
    # Parses the command line arguments and returns the namespace

if __name__ == '__main__':
    main()
