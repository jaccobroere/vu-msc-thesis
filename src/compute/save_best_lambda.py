import os
import re
import sys

import numpy as np
import pandas as pd


def parse_float(s):
    pattern = r"\d+(?:\.\d+)?(?:[eE][-+]?\d+)?"  # regular expression pattern to match scientific notation and floating-point numbers
    match = re.search(pattern, s)  # search for the pattern in the string
    if match:
        return float(
            match.group()
        )  # extract the matched string and convert it to a float
    else:
        return None  # return None if no match is found


# get project directory path
proj_dir = os.environ["PROJ_DIR"]
os.chdir(proj_dir)

# read CLI arguments
path_prefix = "designA_T500_p25"  # sys.argv[1]
sim_id = "1343"  # sys.argv[2]

# set up directory paths
data_dir = os.path.join(proj_dir, "data", "simulation")
out_dir = os.path.join(proj_dir, "out")
lambdas_dir = os.path.join(out_dir, "simulation", "lambdas")
sim_id_dir = os.path.join(lambdas_dir, sim_id)

# read in the lambda dataframes
grid_lam_reg_a0 = pd.read_csv(os.path.join(sim_id_dir, "reg_a0.csv"), delimiter=",")
grid_lam_reg_a05 = pd.read_csv(os.path.join(sim_id_dir, "reg_a05.csv"), delimiter=",")
grid_lam_sym_a0 = pd.read_csv(os.path.join(sim_id_dir, "sym_a0.csv"), delimiter=",")
grid_lam_sym_a05 = pd.read_csv(os.path.join(sim_id_dir, "sym_a05.csv"), delimiter=",")

# calculate the average for each column
avg_lam_reg_a0 = np.mean(grid_lam_reg_a0, axis=0)
avg_lam_reg_a05 = np.mean(grid_lam_reg_a05, axis=0)
avg_lam_sym_a0 = np.mean(grid_lam_sym_a0, axis=0)
avg_lam_sym_a05 = np.mean(grid_lam_sym_a05, axis=0)


# calculate what lambda was the best
best_lam_reg_a0 = parse_float(avg_lam_reg_a0.index[np.argmin(avg_lam_reg_a0)])
best_lam_reg_a05 = parse_float(avg_lam_reg_a05.index[np.argmin(avg_lam_reg_a05)])
best_lam_sym_a0 = parse_float(avg_lam_sym_a0.index[np.argmin(avg_lam_sym_a0)])
best_lam_sym_a05 = parse_float(avg_lam_sym_a05.index[np.argmin(avg_lam_sym_a05)])


# save the grid of best lambdas selected lambda values
best_lambdas = pd.DataFrame(
    {
        "model": [
            "best_lam_reg_a0",
            "best_lam_reg_a05",
            "best_lam_sym_a0",
            "best_lam_sym_a05",
        ],
        "lambda": [
            best_lam_reg_a0,
            best_lam_reg_a05,
            best_lam_sym_a0,
            best_lam_sym_a05,
        ],
    }
)
# %%
best_lambdas.to_csv(
    os.path.join(out_dir, sim_id_dir, f"{path_prefix}_best_lambdas.csv"),
    index=False,
    header=True,
)
