import os

import pandas as pd
from tabulate import tabulate

# Read enivornment voriable
PROJ_DIR = os.environ["PROJ_DIR"]
os.chdir(PROJ_DIR)


# %%
def create_full_data_dictionary(design: str = "designB"):
    # Get the path to the out/simulation/fit/ directory
    fit_dir = os.path.join(os.getcwd(), "out/simulation/fit/")
    data = {}

    # Iterate over the directories in the out/simulation/fit/ directory
    for design_dir in os.listdir(fit_dir):
        if design_dir.startswith("designB"):
            # Get the UUIDs of the subdirectories in the design_dir directory
            uuids = os.listdir(os.path.join(fit_dir, design_dir))

            # Iterate over the UUIDs
            for uuid in uuids:
                # Get the path to the csv files in the uuid directory
                uuid_dir = os.path.join(fit_dir, design_dir, uuid)
                csv_files = os.listdir(uuid_dir)

                # Iterate over the csv files
                for csv_file in csv_files:
                    splits = csv_file.split("_")
                    print(splits)
                    
                    # If the csv file starts with the prefix "splash_", then it is a model output
                    if csv_file.startswith("splash_"):
                        # Read the csv file into a pandas dataframe
                        df = pd.read_csv(
                            os.path.join(fit_dir, design_dir, uuid, csv_file)
                        )

                        # Add the dataframe to the data dictionary
                        data[design_dir] = data.get(design_dir, {})
                        data[design_dir][uuid] = df

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
    create_full_data_dictionary()


if __name__ == "__main__":
    main()
