import os
import pickle
import re

import pandas as pd
from tabulate import tabulate

# Read enivornment voriable
PROJ_DIR = os.environ["PROJ_DIR"]
os.chdir(PROJ_DIR)


import re
import os
import pandas as pd
import pickle

def parse_design_id(design_id: str) -> tuple:
    """
    Parse the design ID from the directory name.
    
    Args:
        design_id (str): Directory name.
    
    Returns:
        tuple: Tuple containing design ID, T, and p.
    """
    match = re.match(r"(\w+)_T(\d+)_p(\d+)", design_id)
    if match:
        design_id = match.group(1)
        T, p = int(match.group(2)), int(match.group(3))
        return design_id, T, p

    raise ValueError(f"Invalid design_id: {design_id}")

def get_csv_files(directory):
    """
    Get a list of CSV files from a directory.

    Args:
        directory (str): Path to the directory.
    
    Returns:
        list: List of CSV file names.
    """
    return [file for file in os.listdir(directory) if file.endswith('.csv')]

def parse_csv_file(file):
    """
    Parse the CSV file name.

    Args:
        file (str): CSV file name.
    
    Returns:
        tuple: Tuple containing the model name and item. 
               If item is not present, returns None.
    """
    splits = file.split(".")[0].split("__")
    return splits[0], splits[1] if len(splits) > 1 else None

def create_data_from_csv_files(directory):
    """
    Create a dictionary from CSV files in a directory.

    Args:
        directory (str): Path to the directory.
    
    Returns:
        dict: Dictionary containing data from CSV files.
    """
    data = {}
    csv_files = get_csv_files(directory)

    for csv_file in csv_files:
        model_name, item = parse_csv_file(csv_file)
        if item is None:
            continue

        df = pd.read_csv(os.path.join(directory, csv_file), header=0)
        if model_name not in data:
            data[model_name] = {}
        data[model_name][item] = df

    return data

def create_full_data_dictionary(design: str = "designB", dump: bool = False):
    """
    Create a full data dictionary from a directory structure.

    Args:
        design (str, optional): Design name. Defaults to "designB".
        dump (bool, optional): If True, dumps the data dictionary to a pickle file. 
                               Defaults to False.
    
    Returns:
        dict: Dictionary containing the full data.
    """
    fit_dir = os.path.join(os.getcwd(), "out/simulation/fit/")
    if not os.path.exists(fit_dir):
        raise FileNotFoundError(f"Directory not found: {fit_dir}")
    
    data = {}

    for design_dir in os.listdir(fit_dir):
        design_id, T, p = parse_design_id(design_dir)

        if design_id not in data:
            data[design_id] = {}

        uuid_dir = os.path.join(fit_dir, design_dir)

        for uuid in os.listdir(uuid_dir):
            data_dir = os.path.join(uuid_dir, uuid)
            data[design_id][uuid] = create_data_from_csv_files(data_dir)

    if dump:
        with open(f"out/simulation/fit/{design}_data.pkl", "wb") as f:
            pickle.dump(data, f)

    return data

def aggregate_to_dataframe(data: dict):
    


def write_table_to_latex(df: pd.DataFrame):
    # Define the columns
    columns = [
        "T",
        "p",
        "GF-SPLASH1",
        "GF-SPLASH2",
        "GF-SPLASH3",
        "F-SPLASH",
        "SSF-SPLASH",
        "SPLASH1",
        "SPLASH2",
        "PVAR",
    ]

    # Create an empty DataFrame with the specified columns
    df = pd.DataFrame(columns=columns)

    # Add rows to the DataFrame
    # Here you will replace this with your actual data
    for i in range(10):  # replace 10 with the number of rows you want
        df.loc[i] = [0] * len(columns)  # replace [0]*len(columns) with your actual data

    # Convert the DataFrame to a LaTeX table
    latex_table = tabulate(df, tablefmt="latex_raw", headers="keys", showindex=False)

    # Write the LaTeX table to a file
    with open("table.tex", "w") as f:
        f.write(latex_table)


def main():
    return create_full_data_dictionary()


if __name__ == "__main__":
    data = main()
