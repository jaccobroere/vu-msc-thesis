//
//  gen_lasso.h
//  gen_lasso
//
//  Created by Yunzhang Zhu on 10/10/14.
//  Copyright (c) 2014 OSU. All rights reserved.
//

#ifndef __gen_lasso__gen_lasso__
#define __gen_lasso__gen_lasso__

#include <stdio.h>
#include "dmatrix.h"
#include "util.h"

/*
 * structure for problem data
 */
enum penalty_types
{
    general_L1,
    graph_L1,
    graph_sparse_L1
};

enum varying_rho_types
{
    adaptive,       // change rho adaptively using the varying penalty schemes.
    non_adaptive,   // change rho for different A's
    constant        // doesn't change rho, user specify a value and use it for all A's
    
};

enum algorithm_type
{
    special_admm,
    standard_admm,
    linearized_admm
};

typedef struct
{
    int p, n, m, lambda_graph_grid, lambda_sparse_grid,lambda_pairs_grid;
    double *Y, *X, *lambda_graph, *lambda_sparse,*lambda_pairs,
    *fun_vals;
    dmatrix_csr* A;
    penalty_types pen_types;
    varying_rho_types varying_rho_scheme;
} linreg_data;

typedef struct
{
    double *theta, *eps_abs, *eps_rel, mu, tau,rho;
    int *max_num_iter;
    bool diag_mat_given, general_A, reporting;
    algorithm_type algorithm;
} parameters;

typedef struct
{
    int npm, j,k, *degree,varying_rho_points;
    bool skinny;
    double  prires, dualres, eps_pri, eps_dual, *prires_vec,
    *dualres_vec, *beta, *beta_old, *alpha, *alpha1,
    *alpha2, *Xtalpha3, *alpha_bar, *Abeta, *alpha_barA, *alpha2A,
    *diag_mat, *diag_mat_inv, *L, *XTY, *temp, *Xbeta;
    dmatrix_csc* A_csc;
    bool* graph;
} tmpvars;


void linreg_path(int* beta_dim, int* A_row_dim,
                 int* sample_size, penalty_types* pen_types,
                 varying_rho_types* varying_rho_scheme,
                 double* Y,
                 double* X, dmatrix_csr* A,
                 double* theta,
                 double* lambda_graph, int* lambda_graph_grid,
                 double* lambda_sparse, int* lambda_sparse_grid,
                 return_data* out, int* num_iter,
                 algorithm_type* algorithm,
                 bool reporting, double* rho,
                 double* eps_abs, double* eps_rel);

void linreg_path_v2(int* beta_dim, int* A_row_dim,
                    int* sample_size,
                    varying_rho_types* varying_rho_scheme,
                    double* Y,
                    double* X, dmatrix_csr* A,
                    double* theta,
                    double* lambda_pairs, int* num_of_pairs,
                    return_data* out, int* num_iter,
                    algorithm_type* algorithm,
                    double* diag_mat, bool diag_mat_given,
                    bool reporting, double* rho,
                    double* eps_abs, double* eps_rel);


#endif /* defined(__gen_lasso__gen_lasso__) */
