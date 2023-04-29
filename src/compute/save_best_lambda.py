# This script walks the out/simlulation/lambdas directory and
# saves the best lambda value for each simulation design, which resides in a subfolder

import os

import numpy as np
import pandas as pd
from utils import parse_float


def save_best_lam(subdir):
    sim_id_dir = os.path.join(lambdas_dir, subdir)

    # read in the lambda dataframes
    grid_lam_reg_a0 = pd.read_csv(os.path.join(sim_id_dir, "reg_a0.csv"), delimiter=",")
    grid_lam_reg_a05 = pd.read_csv(
        os.path.join(sim_id_dir, "reg_a05.csv"), delimiter=","
    )
    grid_lam_sym_a0 = pd.read_csv(os.path.join(sim_id_dir, "sym_a0.csv"), delimiter=",")
    grid_lam_sym_a05 = pd.read_csv(
        os.path.join(sim_id_dir, "sym_a05.csv"), delimiter=","
    )
    grid_lam_spl_a0 = pd.read_csv(os.path.join(sim_id_dir, "spl_a0.csv"), delimiter=",")
    grid_lam_spl_a05 = pd.read_csv(
        os.path.join(sim_id_dir, "spl_a05.csv"), delimiter=","
    )

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
        os.path.join(out_dir, sim_id_dir, f"{subdir}_best_lambdas.csv"),
        index=False,
        header=True,
    )


if __name__ == "__main__":
    # get project directory path
    proj_dir = os.environ["PROJ_DIR"]
    os.chdir(proj_dir)

    # set up directory paths
    data_dir = os.path.join(proj_dir, "data", "simulation")
    out_dir = os.path.join(proj_dir, "out")
    lambdas_dir = os.path.join(out_dir, "simulation", "lambdas")

    # Save the best lambda value for each simulation design
    for subdir in os.listdir(lambdas_dir):
        if os.path.isdir(os.path.join(lambdas_dir, subdir)):
            save_best_lam(subdir)
