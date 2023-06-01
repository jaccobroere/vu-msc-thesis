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
        r"[a-f0-9]{8}-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]" r"{4}-[a-f0-9]{12}",
        uuid.lower(),
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
        # if item in [None, "y_pred"]:
        if item in [None]:
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
    fit_dir = os.path.join(os.getcwd(), "out/simulation/final/")
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
        with open(f"out/simulation/final/{design}_data.pkl", "wb") as f:
            pickle.dump(data, f)

    return data


def calc_rmsfe(y_true, y_pred, y_hat_true):
    """
    Calculate the RMSFE.

    Args:
        y_true (np.array): True values.
        y_pred (np.array): Predicted values.
        y_hat_true (np.array): True values of the y_hat model.

    Returns:
        float: RMSFE.
    """
    return np.sqrt(np.sum((y_true - y_pred) ** 2) / np.sum((y_true - y_hat_true) ** 2))


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

    df_rmsfe_h1 = pd.DataFrame(columns=columns)
    df_rmsfe_long = pd.DataFrame(columns=columns)

    for design_id in data:
        for uuid in data[design_id]:
            # Skip if there is no data for a certain UUID
            if len(data[design_id][uuid]) == 0:
                continue

            y_true = pd.read_csv(
                os.path.join(
                    "out", "simulation", "final", design_id, uuid, "y_true.csv"
                ),
                header=0,
            ).to_numpy()
            y_test = y_true[:, int(np.floor(y_true.shape[1] * 0.8)) :]
            y_hat_true = pd.read_csv(
                os.path.join(
                    "out", "simulation", "final", design_id, uuid, "y_hat_true.csv"
                ),
                header=0,
            ).to_numpy()
            # Initialize the dictionary to save RMSFE values in
            temp_dict_h1 = {key: [0] for key in columns}
            temp_dict_long = {key: [0] for key in columns}
            for model_name, item in data[design_id][uuid].items():
                # Get the RMSFE value
                temp_dict_h1[model_name] = item["rmsfe"].to_numpy().flatten()
                # Assign the design_id to a column to group by later
                y_pred = item["y_pred"].to_numpy()
                temp_dict_long[model_name] = calc_rmsfe(
                    y_test, y_pred, y_hat_true
                ).flatten()

            temp_dict_h1["design_id"] = design_id
            temp_dict_long["design_id"] = design_id
            temp_df_rmsfe_h1 = pd.DataFrame(temp_dict_h1)
            temp_df_rmsfe_long = pd.DataFrame(temp_dict_long)
            df_rmsfe_h1 = pd.concat([df_rmsfe_h1, temp_df_rmsfe_h1], ignore_index=True)
            df_rmsfe_long = pd.concat(
                [df_rmsfe_long, temp_df_rmsfe_long], ignore_index=True
            )

    return (
        df_rmsfe_h1.groupby("design_id").mean(),
        df_rmsfe_long.groupby("design_id").mean(),
    )


def read_AB_true(design_id: str, uuid: str):
    A_true = pd.read_csv(
        os.path.join("out", "simulation", "final", design_id, uuid, "A_true.csv"),
        header=0,
    )
    B_true = pd.read_csv(
        os.path.join("out", "simulation", "final", design_id, uuid, "B_true.csv"),
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


def write_table_to_latex(df: pd.DataFrame, filename: str = None):
    temp = df.copy()

    # Define the columns
    columns = [
        r"F-SPL($\lambda$)",
        r"SSF-SPL($0.5, \lambda$)",
        r"GF-SPL($0.5, 0, \lambda$)",
        r"GF-SPL($0, 1, \lambda$)",
        r"GF-SPL($0.5, 1, \lambda$)",
        r"SPLASH($0, \lambda$)",
        r"SPLASH($0.5, \lambda$)",
        r"PVAR($\lambda$)",
    ]

    # Create an empty DataFrame with the specified columns
    temp.columns = columns

    # Find the element that is the mininum of each row
    argmin = temp.loc[:, columns].idxmin(axis=1)

    # Add the T and p values as columns and set as index
    temp["design_id"] = temp.index
    temp["T"] = temp["design_id"].apply(lambda x: parse_design_id(x)[1])
    temp["p"] = temp["design_id"].apply(lambda x: parse_design_id(x)[2])
    temp.drop("design_id", axis=1, inplace=True)
    temp.sort_values(["p", "T"], inplace=True)
    temp = temp.loc[:, ["p", "T"] + columns]

    # Round to three decimals and convert to string type
    temp = temp.round(3).applymap("{:.3f}".format)

    # Replace nan values with a dash for formatting
    temp.replace("nan", "-", inplace=True)

    # Make the argmin value bold
    for i, row in temp.iterrows():
        temp.loc[i, argmin[i]] = r"\textbf{" + temp.loc[i, argmin[i]] + "}"

    # Convert the DataFrame to a LaTeX table
    latex_table = tabulate(
        temp, tablefmt="latex_raw", headers="keys", showindex=False, numalign="center"
    )

    # Write the LaTeX table to a file
    if filename is not None:
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


def combine_tables(tables, design="designB", filename: str = None):
    spl = design[:-1], design[-1]
    design_cap, letter = spl[0].capitalize(), spl[1]
    header = f"""\\begin{{landscape}}
    \\bgroup
    \\def\\arraystretch{{1.3}}
    \\begin{{table}}[p]
    \\footnotesize
    \\centering
    \\caption{{Simulation results for {design_cap} {letter}}}
    \\label{{tab:results_{design}}}
    \\begin{{tabular}}{{cccccccccc}}    
    \\hline \\hline
    $p$  &  $T$   &  F-SPL($\lambda$)  & SSF-SPL($0.5, \lambda$)  &  GF-SPL($0.5, 0, \lambda$)  &  GF-SPL($0, 1, \lambda$)  &  GF-SPL($0.5, 1, \lambda$)  &  SPLASH($0, \lambda$)  &  SPLASH($0.5, \lambda$)  &  PVAR($\lambda$)  \\\\
    \\hline
    """

    subheaders = [
        f"\\multicolumn{{10}}{{l}}{{\\textbf{{RMSFE}}}} \\\\",
        f"\t\\multicolumn{{10}}{{l}}{{$\\mathbf{{EE_A}}$}} \\\\",
        f"\t\\multicolumn{{10}}{{l}}{{$\\mathbf{{EE_B}}$}} \\\\",
    ]

    for i, table in enumerate(tables):
        lines = table.split("\n")
        header += subheaders[i] + "\n"
        header += "\t" + r"\hline" + "\n"

        for j, line in enumerate(lines[4:-2]):
            if j == len(lines[4:-2]) - 1:
                header += "\t" + line + "\n"
            else:
                header += "\t" + line + r" \hdashline" + "\n"

        header += "\t" + r"\hline" + "\n"

    tail = f"""\t\\hline
    \\multicolumn{{10}}{{l}}{{\\textbf{{Note:}} Numbers in \\textbf{bold} inidicate best performance for each combination of $p$ and $T$}} \\\\
    \\multicolumn{{10}}{{l}}{{\\textbf{{Note:}} Simulation results are based on $N_\\text{{sim}} = 250$ Monte Carlo simulations}}
    \\end{{tabular}}
    \\end{{table}}
    \\egroup
\\end{{landscape}}
    """
    header += tail

    if filename is not None:
        with open(
            os.path.join(
                "out",
                "tables",
                filename,
            ),
            "w",
        ) as f:
            f.write(header)

    return header


def main(design: str = "designB"):
    data = create_full_data_dictionary(design)
    df_rmsfe_h1, df_rmsfe_long = collect_rmsfe_data(data)
    df_A, df_B = collect_estimation_error_data(data)
    latex_table_rmsfe = write_table_to_latex(df_rmsfe_h1, f"{design}_rmsfe.tex")
    latex_table_rmsfe_long = write_table_to_latex(
        df_rmsfe_long, f"{design}_rmsfe_long.tex"
    )
    latex_table_A = write_table_to_latex(df_A, f"{design}_EEA.tex")
    latex_table_B = write_table_to_latex(df_B, f"{design}_EEB.tex")
    full_table = combine_tables(
        [latex_table_rmsfe, latex_table_A, latex_table_B],
        design,
        f"results_{design}.tex",
    )

    return (
        full_table,
        latex_table_rmsfe,
        latex_table_rmsfe_long,
        latex_table_A,
        latex_table_B,
    )


if __name__ == "__main__":
    designs = ["designA", "designB", "designC"]
    for design in designs:
        main(design)
