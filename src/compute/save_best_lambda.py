# This script walks the out/simlulation/lambdas directory and
# saves the best lambda value for each simulation design, which resides in a subfolder

import os
import sys

import numpy as np
import pandas as pd
from utils import parse_float


def collect_dataframes(lambdas_dir, model_name):
    res = pd.DataFrame()

    for uuid_dir in os.listdir(lambdas_dir):
        path = os.path.join(lambdas_dir, uuid_dir)
        if os.path.isdir(path):
            df = pd.read_csv(
                os.path.join(path, f"{model_name}.csv"),
                delimiter=",",
                header=0,
            )
            res = pd.concat([res, df], axis=0, ignore_index=True)

    res.columns = df.columns

    return res


def save_best_lam(lambdas_dir):
    # Read in the dataframes
    grid_lam_reg_a0 = collect_dataframes(lambdas_dir, "reg_a0")
    grid_lam_reg_a05 = collect_dataframes(lambdas_dir, "reg_a05")
    grid_lam_sym_a0 = collect_dataframes(lambdas_dir, "sym_a0")
    grid_lam_sym_a05 = collect_dataframes(lambdas_dir, "sym_a05")
    grid_lam_spl_a0 = collect_dataframes(lambdas_dir, "spl_a0")
    grid_lam_spl_a05 = collect_dataframes(lambdas_dir, "spl_a05")

    # calculate the average for each column
    avg_lam_reg_a0 = np.mean(grid_lam_reg_a0, axis=0)
    avg_lam_reg_a05 = np.mean(grid_lam_reg_a05, axis=0)
    avg_lam_sym_a0 = np.mean(grid_lam_sym_a0, axis=0)
    avg_lam_sym_a05 = np.mean(grid_lam_sym_a05, axis=0)
    avg_lam_spl_a0 = np.mean(grid_lam_spl_a0, axis=0)
    avg_lam_spl_a05 = np.mean(grid_lam_spl_a05, axis=0)

    # calculate what lambda was the best
    best_lam_reg_a0 = parse_float(avg_lam_reg_a0.index[np.argmin(avg_lam_reg_a0)])
    best_lam_reg_a05 = parse_float(avg_lam_reg_a05.index[np.argmin(avg_lam_reg_a05)])
    best_lam_sym_a0 = parse_float(avg_lam_sym_a0.index[np.argmin(avg_lam_sym_a0)])
    best_lam_sym_a05 = parse_float(avg_lam_sym_a05.index[np.argmin(avg_lam_sym_a05)])
    best_lam_spl_a0 = parse_float(avg_lam_spl_a0.index[np.argmin(avg_lam_spl_a0)])
    best_lam_spl_a05 = parse_float(avg_lam_spl_a05.index[np.argmin(avg_lam_spl_a05)])

    # save the grid of best lambdas selected lambda values
    best_lambdas = pd.DataFrame(
        {
            "model": [
                "best_lam_reg_a0",
                "best_lam_reg_a05",
                "best_lam_sym_a0",
                "best_lam_sym_a05",
                "best_lam_spl_a0",
                "best_lam_spl_a05",
            ],
            "lambda": [
                best_lam_reg_a0,
                best_lam_reg_a05,
                best_lam_sym_a0,
                best_lam_sym_a05,
                best_lam_spl_a0,
                best_lam_spl_a05,
            ],
        }
    )

    best_lambdas.to_csv(
        os.path.join(lambdas_dir, "best_lambdas.csv"),
        index=False,
        header=True,
    )


if __name__ == "__main__":
    # Read CLI arguments
    sim_design_id = sys.argv[1]

    # get project directory path
    proj_dir = os.environ["PROJ_DIR"]
    os.chdir(proj_dir)

    # set up directory paths
    lambdas_dir = os.path.join(
        "out", "simulation", "lambdas", f"{sim_design_id}_detlam"
    )

    save_best_lam(lambdas_dir)
