import os
import pickle
import re

import pandas as pd
from tabulate import tabulate

# Read enivornment voriable
PROJ_DIR = os.environ["PROJ_DIR"]
os.chdir(PROJ_DIR)


def parse_design_id(design_id: str) -> tuple:
    match = re.match(r"(\w+)_T(\d+)_p(\d+)", design_id)
    if match:
        design_id = match.group(1)
        T, p = int(match.group(2)), int(match.group(3))
        return design_id, T, p

    raise ValueError(f"Invalid design_id: {design_id}")


# %%
def create_full_data_dictionary(design: str = "designB", dump: bool = False):
    # Get the path to the out/simulation/fit/ directory
    fit_dir = os.path.join(os.getcwd(), "out/simulation/fit/")
    data = {}

    # Iterate over the directories in the out/simulation/fit/ directory
    for design_dir in os.listdir(fit_dir):
        design_id, T, p = parse_design_id(design_dir)
        data[design_dir] = data.get(design_id, {})

        # Get the UUIDs of the subdirectories in the design_dir directory
        uuids = os.listdir(os.path.join(fit_dir, design_dir))

        # Iterate over the UUIDs
        for uuid in uuids:
            # Create a dictionary for the uuid
            data[design_dir][uuid] = data[design_dir].get(uuid, {})

            # Get the path to the csv files in the uuid directory
            uuid_dir = os.path.join(fit_dir, design_dir, uuid)
            csv_files = os.listdir(uuid_dir)

            # Iterate over the csv files
            for csv_file in csv_files:
                splits = csv_files.split(".")[0].split(
                    "__"
                )  # Remove .csv from filename and then split on __
                if len(splits) < 2:
                    continue

                model_name = splits[0]
                item = splits[1]

                # Read the csv file into a pandas dataframe
                df = pd.read_csv(
                    os.path.join(fit_dir, design_dir, uuid, csv_file), header=0
                )

                # Add the dataframe to the data dictionary
                data[design_dir][uuid][model_name] = data[design_dir][uuid].get(
                    model_name, {}
                )
                data[design_dir][uuid][model_name][item] = df

    # Dump the dictionary to a pickle file
    if dump:
        with open(f"out/simulation/fit/{design}_data.pkl", "wb") as f:
            pickle.dump(data, f)

    return data


def aggregate_to_dataframe(d: dict):
    pass


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
