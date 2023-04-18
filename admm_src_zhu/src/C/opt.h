#include "dmatrix.h"

typedef struct
{
    int p, n, m, lambda_grid;
    double *Y, *X, *lambda,
    *fun_vals;
    dmatrix_csr* A;
} linreg_data;

typedef struct
{
    int *chol_num, *admm_iter_num;
    double *fun_val, *fun_path, *beta_path;
    
} return_data;

typedef struct
{
    double *theta, *rho, *eps_abs, *eps_rel, mu, tau;
    int *max_num_iter;
    bool variant, varying_rho, general_A, reporting;
} parameters;

typedef struct
{
    int npm, k;
    bool skinny;
    double  prires, dualres, eps_pri, eps_dual, *prires_vec,
    *dualres_vec, *beta, *beta_old, *alpha, *alpha1,
    *alpha2, *alpha_bar, *Abeta, *alpha_barA, *alpha2A,
    *diag_mat, *D, *diag_mat_inv, *L, *XTY, *temp;
    dmatrix_csc* A_csc;
    bool* graph;
} tmpvars;



void filtering_old(int* beta_dim, int* A_row_dim,
               double* Y, int* A_idx,
               int* A_jdx, int* A_rdx, int* nz,
               double* A_val, double* beta,
               double* alpha, double* theta,
               double* lambda, int* lambda_grid,
               double* rho, int* num_iter,double* eps);

void filtering(int* beta_dim, int* A_row_dim,
               double* Y, int* A_idx,
               int* A_jdx, int* A_rdx, int* nz,
               double* A_val, double* beta,
               double* theta,
               double* lambda, int* lambda_grid,
               double* rho, int* num_iter,
               double* eps, double* eps_abs,
               double* eps_rel, bool general_A,
               bool varying_rho);


void linreg_path(int* beta_dim, int* A_row_dim,
                     int* sample_size, double* Y,
                     double* X, dmatrix_csr* A,
                     double* theta,
                     double* lambda, int* lambda_grid,
                     return_data* out,
                     double* rho, int* num_iter,
                     bool variant, bool varying_rho,
                     bool general_A, bool reporting,
                     double* eps_abs, double* eps_rel);
