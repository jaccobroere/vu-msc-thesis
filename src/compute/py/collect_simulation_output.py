import os
import pickle
import re

import numpy as np
import pandas as pd
from tabulate import tabulate

# Read enivornment voriable
PROJ_DIR = os.environ["PROJ_DIR"]
os.chdir(PROJ_DIR)


def is_uuid(uuid: str) -> bool:
    """
    Check whether a string is a UUID.

    Args:
        uuid (str): String to check.

    Returns:
        bool: True if the string is a UUID, False otherwise.
    """
    return re.match(
        r"[a-f0-9]{8}-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]" r"{4}-[a-f0-9]{12}", uuid
    )


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

    print(f"WARNING: Could not parse design ID: {design_id}")
    return None, None, None


def get_csv_files(directory):
    """
    Get a list of CSV files from a directory.

    Args:
        directory (str): Path to the directory.

    Returns:
        list: List of CSV file names.
    """
    return [file for file in os.listdir(directory) if file.endswith(".csv")]


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
        if item in [None, "y_pred"]:
            continue

        df = pd.read_csv(os.path.join(directory, csv_file), header=0)
        if model_name not in data:
            data[model_name] = {}
        data[model_name][item] = df

    return data


def create_full_data_dictionary(design: str = "designB", dump: bool = False):
    """
        Create a full data dictionary from a directory structure.

        The returned dictionary looks like this:

        {
        'design_id1': {
            'uuid1': {
                'model_name1': {
                    'item1': DataFrame1,
                    'item2': DataFrame2,
                    ...
                },
                'model_name2': {
                    'item1': DataFrame3,
                    'item2': DataFrame4,
                    ...
                },
                ...
            },
            'uuid2': {
                'model_name1': {
                    'item1': DataFrame5,
                    'item2': DataFrame6,
                    ...
                },
                'model_name2': {
                    'item1': DataFrame7,
                    'item2': DataFrame8,
                    ...
                },
                ...
            },
            ...
        },
        'design_id2': {
            ...
        },
        ...
    }


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

        # Check whether the design ID matches the specified design and whether it is a design at all
        if design_id is None or design_id != design:
            continue

        # Create dictionary if not present
        data[design_dir] = data.get(design_dir, {})

        uuid_dir = os.path.join(fit_dir, design_dir)

        # Walk through the UUID directories and create data from CSV files
        for uuid in os.listdir(uuid_dir):
            if not is_uuid(uuid):
                continue

            data_dir = os.path.join(uuid_dir, uuid)
            data[design_dir][uuid] = create_data_from_csv_files(data_dir)

    if dump:
        with open(f"out/simulation/fit/{design}_data.pkl", "wb") as f:
            pickle.dump(data, f)

    return data


def collect_rmsfe_data(data: dict):
    columns = [
        "design_id",
        "fsplash",
        "ssfsplash",
        "gfsplash_a05",
        "gfsplash_sym_a0",
        "gfsplash_sym_a05",
        "splash_a0",
        "splash_a05",
        "pvar",
    ]

    df = pd.DataFrame(columns=columns)

    for design_id in data:
        for uuid in data[design_id]:
            # Skip if there is no data for a certain UUID
            if len(data[design_id][uuid]) == 0:
                continue

            # Initialize the dictionary to save RMSFE values in
            temp_dict = {key: [0] for key in columns}
            for model_name, item in data[design_id][uuid].items():
                # Get the RMSFE value
                temp_dict[model_name] = item["rmsfe"].to_numpy().flatten()
                # Assign the design_id to a column to group by later

            temp_dict["design_id"] = design_id
            temp_df = pd.DataFrame(temp_dict)
            df = pd.concat([df, temp_df], ignore_index=True)

    return df.groupby("design_id").mean()


def read_AB_true(design_id: str, uuid: str):
    A_true = pd.read_csv(
        os.path.join("out", "simulation", "fit", design_id, uuid, "A_true.csv"),
        header=0,
    )
    B_true = pd.read_csv(
        os.path.join("out", "simulation", "fit", design_id, uuid, "B_true.csv"),
        header=0,
    )
    return A_true.to_numpy(), B_true.to_numpy()


def collect_estimation_error_data(data: dict):
    columns = [
        "design_id",
        "fsplash",
        "ssfsplash",
        "gfsplash_a05",
        "gfsplash_sym_a0",
        "gfsplash_sym_a05",
        "splash_a0",
        "splash_a05",
        "pvar",
    ]

    df_A = pd.DataFrame(columns=columns)
    df_B = pd.DataFrame(columns=columns)

    for design_id in data:
        for uuid in data[design_id]:
            # Skip if there is no data for a certain UUID
            if len(data[design_id][uuid]) == 0:
                continue

            # Initialize the dictionary to save EE values in
            temp_dict_A = {key: [0] for key in columns}
            temp_dict_B = {key: [0] for key in columns}
            A_true, B_true = read_AB_true(design_id, uuid)

            for model_name, item in data[design_id][uuid].items():
                temp_dict_A[model_name] = [np.nan]
                temp_dict_B[model_name] = [np.nan]
                if model_name == "pvar":
                    continue
                # Get the EEA value as the spectral norm of the differences
                temp_dict_A[model_name] = np.linalg.norm(
                    A_true - item["estimate_A"].to_numpy(), ord=2
                )
                temp_dict_B[model_name] = np.linalg.norm(
                    B_true - item["estimate_B"].to_numpy(), ord=2
                )
                # Assign the design_id to a column to group by later

            temp_dict_A["design_id"] = design_id
            temp_dict_B["design_id"] = design_id
            temp_df_A = pd.DataFrame(temp_dict_A)
            temp_df_B = pd.DataFrame(temp_dict_B)
            df_A = pd.concat([df_A, temp_df_A], ignore_index=True)
            df_B = pd.concat([df_B, temp_df_B], ignore_index=True)

    return df_A.groupby("design_id").mean(), df_B.groupby("design_id").mean()


def write_table_to_latex(df: pd.DataFrame, filename: str):
    # Define the columns
    columns = [
        r"F-SPLASH($\lambda$)",
        r"SSF-SPLASH($\alpha=0.5$, $\lambda$)",
        r"GF-SPLASH($\alpha=0.5$, $\lambda$, $\sigma=0$)",
        r"GF-SPLASH($\alpha=0$, $\lambda$, $\sigma=1$)",
        r"GF-SPLASH($\alpha=0.5$, $\lambda$, $\sigma=1$)",
        r"SPLASH($0$, $\lambda$)",
        r"SPLASH($0.5$, $\lambda$)",
        r"PVAR($\lambda$)",
    ]

    # Create an empty DataFrame with the specified columns
    df.columns = columns

    # Add the T and p values as columns and set as index
    df["design_id"] = df.index
    df["T"] = df["design_id"].apply(lambda x: parse_design_id(x)[1])
    df["p"] = df["design_id"].apply(lambda x: parse_design_id(x)[2])
    df.drop("design_id", axis=1, inplace=True)
    df.sort_values(["p", "T"], inplace=True)
    df = df.loc[:, ["p", "T"] + columns]
    # df.set_index(["T", "p"], inplace=True)

    # Convert the DataFrame to a LaTeX table
    latex_table = tabulate(
        df, tablefmt="latex_raw", headers="keys", showindex=False, numalign="center"
    )

    # Write the LaTeX table to a file
    with open(
        os.path.join(
            "out",
            "tables",
            filename,
        ),
        "w",
    ) as f:
        f.write(latex_table)

    return latex_table


def main(design: str = "designB"):
    data = create_full_data_dictionary(design)
    df_rmsfe = collect_rmsfe_data(data)
    df_A, df_B = collect_estimation_error_data(data)
    latex_table_rmsfe = write_table_to_latex(df_rmsfe, f"{design}_rmsfe.tex")
    latex_table_A = write_table_to_latex(df_A, f"{design}_EEA.tex")
    latex_table_B = write_table_to_latex(df_B, f"{design}_EEB.tex")

    return latex_table_rmsfe, latex_table_A, latex_table_B


if __name__ == "__main__":
    design = "designB"
    main(design)
