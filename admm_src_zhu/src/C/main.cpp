//
//  main.c
//  general_lasso_xcode
//
//  Created by Yunzhang Zhu on 8/30/14.
//  Copyright (c) 2014 OSU. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "dmatrix.h"
#include "gen_lasso.h"
#include "def.h"
#include "util.h"
#include "filtering.h"



void test_fuse_simple(){
    // a fused filtering example //
    srand((unsigned)time(NULL));
    int p = 500;
    int m = p-1;
    int nz = 2*m;
    double eps = 1e-3;
    int* idx = (int*) malloc(sizeof(int)*nz);
    int* jdx = (int*) malloc(sizeof(int)*nz);
    double* val = (double*) malloc(sizeof(double)*nz);
    int* rdx = (int*) malloc(sizeof(int)*(m+1));
    for (int i = 0; i < nz; i++) {
        idx[i] = i / 2;
        jdx[i] = (i+1) / 2;
        val[i] = ((i % 2) == 0)? 1:-1;
    }
    for (int i = 0; i < m+1; i++) {
        rdx[i] = 2*i;
    }
    dmatrix_csr* M;
    dmat_new_sparse_csr(&M, m, p, nz, idx, jdx, rdx, val);
    double* beta = (double*) malloc(sizeof(double)*p);
    double* alpha = (double*) malloc(sizeof(double)*p);
    dmat_vset(p, .0, beta);
    dmat_vset(p, .0, alpha);
    double* Y = (double*) malloc(sizeof(double)*p);
    dmat_vset(p/2, 1.0, Y);
    dmat_vset(p/2, 2.0, Y+p/2);
    srand(1);
    for (int i = 0; i < p; i++) {
        Y[i] = Y[i] + (((double)rand()/(double)RAND_MAX) - .5) / 5;
    }
    double theta = 1.0;
    int lambda_grid = 10;
    double lambda[10] = {.2,.5,1.0,2.0,5.0,8,10,12,20,50};
    // int lambda_grid = 1;
    // double lambda[6] = {20.0};
    double rho = 1.0;
    int num_iter = 1e5;
    // filtering_old(&p, &m, Y, idx, jdx, rdx, &nz, val, beta, alpha, &theta, lambda, &lambda_grid, &rho, &num_iter,&eps);
    
    bool general_A = false, varying_rho = true, variant = false, reporting=false;
    
    return_data* out = (return_data*) malloc(sizeof(return_data));
    out->beta_path = (double*) malloc(sizeof(double)*p*lambda_grid);
    out->fun_val = (double*) malloc(sizeof(double)*lambda_grid);
    out->chol_num = (int*) malloc(sizeof(int)*lambda_grid);
    out->admm_iter_num = (int*) malloc(sizeof(int)*lambda_grid);
    if (reporting) out->fun_path = (double*) malloc(sizeof(double)*num_iter*lambda_grid);
    
    filtering(&p, &m, Y, idx, jdx, rdx, &nz, val, &theta, lambda,
              &lambda_grid, &rho, &num_iter,&eps,variant, general_A,reporting,varying_rho,out);
    free(M); free(Y); free(beta); free(alpha);
    free(idx); free(jdx); free(val);
}
/*
void test_linreg(){
    srand((unsigned)time(NULL));
    int p = 500;
    int m = p-1;
    int n = 200;
    int nz = 2*m;
    double eps_abs = 1e-5;
    double eps_rel = 1e-4;
    double* X = (double*) malloc(sizeof(double)*n*p);
    srand(1);
    for (int i = 0; i < n*p; i++) {
        X[i] = 2*((double)rand()/(double)RAND_MAX)-1.0;
    }
    int* idx = (int*) malloc(sizeof(int)*nz);
    int* jdx = (int*) malloc(sizeof(int)*nz);
    double* val = (double*) malloc(sizeof(double)*nz);
    int* rdx = (int*) malloc(sizeof(int)*(m+1));
    for (int i = 0; i < nz; i++) {
        idx[i] = i / 2;
        jdx[i] = (i+1) / 2;
        val[i] = ((i % 2) == 0)? 1:-1;
    }
    for (int i = 0; i < m+1; i++) {
        rdx[i] = 2*i;
    }
    dmatrix_csr* M;
    dmat_new_sparse_csr(&M, m, p, nz, idx, jdx, rdx, val);
    double* beta = (double*) malloc(sizeof(double)*p);
    double* beta_initial = (double*) malloc(sizeof(double)*p);
    double* alpha = (double*) malloc(sizeof(double)*p);
    dmat_vset(p/2, 1.0, beta);
    dmat_vset(p/2, .0, beta+(p/2));
    double* Y = (double*) malloc(sizeof(double)*n);
    dmat_yAx(n, p, X, beta, Y);
    for (int i = 0; i < n; i++) {
        Y[i] = Y[i] + (((double)rand()/(double)RAND_MAX) - .5) / 5;
    }
    double theta = 1.0;
    // int lambda_grid = 10;
    // double lambda[10] = {3.0,2.5,2.0,1.5,1.3,1.0,.8,.5,.3,.2};
    int lambda_grid = 7;
    double lambda[7] = {.1,.5,2.0,5.0,10.0,20.0,50.0};
    double rho = 1.0;
    int max_num_iter = 1000;
    double* beta_path = (double*) malloc(sizeof(double)*p*lambda_grid);
    bool variant = true;
    bool varying_rho = true;
    bool general_A = false;
    bool reporting = true;
    return_data* out = (return_data*) malloc(sizeof(return_data));
    out->beta_path = (double*) malloc(sizeof(double)*p*lambda_grid);
    out->fun_val = (double*) malloc(sizeof(double)*lambda_grid);
    out->chol_num = (int*) malloc(sizeof(int)*lambda_grid);
    out->admm_iter_num = (int*) malloc(sizeof(int)*lambda_grid);
    if (reporting) out->fun_path = (double*) malloc(sizeof(double)*max_num_iter*lambda_grid);
    
    linreg_path(&p, &m, &n, Y, X, M,
                &theta, lambda, &lambda_grid,
                out, &rho,
                &max_num_iter, variant, varying_rho,
                general_A, reporting,
                &eps_abs, &eps_rel);
    // linreg_single(&p, &m, &n, Y, X, M, beta_path, &theta, lambda, &rho, &max_num_iter, variant, varying_rho, general_A, &eps_abs, &eps_rel);
    
    if (out->fun_path) free(out->fun_path);
    free(out->admm_iter_num); free(out->chol_num);
    free(out->fun_val); free(out->beta_path);
    free(M);
    free(out);
    free(Y);
    free(beta);
    free(alpha);
    free(beta_initial);
    free(idx); free(jdx); free(val); free(X); free(beta_path);
}
*/
    /* input is row-oriented matrix */
void test_blas(){
    srand((unsigned)time(NULL));
    int p = 2;
    int n = 4;
    double* temp = (double*) malloc(sizeof(double)*max(n, p));
    double* X = (double*) malloc(sizeof(double)*n*p);
    srand(1);
    for (int i = 0; i < n*p; i++) {
        X[i] = 2*((double)rand()/(double)RAND_MAX)-1.0;
    }
    double* beta = (double*) malloc(sizeof(double)*p);
    double* beta_lse = (double*) malloc(sizeof(double)*p);
    dmat_vset(p/2, 1.0, beta);
    dmat_vset(p/2, .0, beta+(p/2));
    srand(1);
    double* Y = (double*) malloc(sizeof(double)*n);
    dmat_yAx(n, p, X, beta, Y);
    for (int i = 0; i < n; i++) {
        Y[i] += (((double)rand()/(double)RAND_MAX) - .5) / 5;
    }
    double* XtX = (double*) malloc(sizeof(double)*p*p);
    double* L = (double*) malloc(sizeof(double)*p*p);
    dmat_B_ATA(n, p, X, XtX);
    double* XtY = (double*) malloc(sizeof(double)*p);
    /*
    double* B = (double*) malloc(sizeof(double)*n*n);
    double* D = (double*) malloc(sizeof(double)*p);
    for (int i = 0; i < p; i++) {
        D[i] = (double)i+1;
    }
    print_dmatrix(X, 1, n*p);
    dmat_B_ADAT(n, p, X, D, B);
    print_dmatrix(B, 1, n*n);
    */
    
    dmat_vcopy(p*p, XtX, L);
    dmat_potrf(L,p);
    dmat_yATx(n, p, X, Y, XtY);
    dmat_vcopy(p, XtY, beta_lse);
    dmat_potrs(L, p, beta_lse);
    dmat_yAx(n, p, X, beta_lse, temp);
    print_dmatrix(X, n, p);
    print_dmatrix(temp, 1, n);
    dmat_yATx(n, p, X, temp, beta_lse);
    print_dmatrix(beta_lse, 1, p);
}

void test_UF(){
    int p = 101;
    int m = 100;
    double* beta = (double*) malloc(sizeof(double)*p);
    double* beta_bar = (double*) malloc(sizeof(double)*p);
    for (int i = 0;  i < p; i++) {
        beta[i] = i+1.0;
    }
    int* idx = (int*) malloc(sizeof(int)*m);
    int* jdx = (int*) malloc(sizeof(int)*m);
    for (int i = 0; i < m/2; i++) {
        idx[i] = i;
        jdx[i] = i+1;
    }
    int* graph = (int*) malloc(sizeof(int)*m);
    for (int i = 0; i < m; i++) {
        graph[i] = 1;
    }
    
    UF TMP(p,beta);
    for (int i = 0; i < m; i++) {
        if (graph[i]) TMP.unite(idx[i], jdx[i]);
    }
    
    for (int i = 0; i < p; i++) {
        beta_bar[i] = TMP.mean[TMP.root(TMP.id[i])];
    }
    
    print_dmatrix(beta_bar, 1, p);
    free(idx); free(jdx);
}

void test_sparse_graph_linreg()
{
    srand((unsigned)time(NULL));
    int p = 200;
    int m = p-1;
    int n = 100;
    int nz = 2*m;
    double eps_abs = 1e-5;
    double eps_rel = 1e-4;
    double* X = (double*) malloc(sizeof(double)*n*p);
    srand(1);
    for (int i = 0; i < n*p; i++) {
        X[i] = 2*((double)rand()/(double)RAND_MAX)-1.0;
    }
    int* idx = (int*) malloc(sizeof(int)*nz);
    int* jdx = (int*) malloc(sizeof(int)*nz);
    double* val = (double*) malloc(sizeof(double)*nz);
    int* rdx = (int*) malloc(sizeof(int)*(m+1));
    for (int i = 0; i < nz; i++) {
        idx[i] = i / 2;
        jdx[i] = (i+1) / 2;
        val[i] = ((i % 2) == 0)? 1:-1;
    }
    for (int i = 0; i < m+1; i++) {
        rdx[i] = 2*i;
    }
    
    dmatrix_csr* M;
    dmat_new_sparse_csr(&M, m, p, nz, idx, jdx, rdx, val);
    double* beta = (double*) malloc(sizeof(double)*p);
    double* beta_initial = (double*) malloc(sizeof(double)*p);
    double* alpha = (double*) malloc(sizeof(double)*p);
    dmat_vset(p/2, 1.0, beta);
    dmat_vset(p/2, .0, beta+(p/2));
    double* Y = (double*) malloc(sizeof(double)*n);
    dmat_yAx(n, p, X, beta, Y);
    for (int i = 0; i < n; i++) {
        Y[i] = Y[i] + (((double)rand()/(double)RAND_MAX) - .5) / 5;
    }
    double theta = 1.0;
//    int lambda_graph_grid = 1;
//    int lambda_sparse_grid = 1;
//    double lambda_graph[1] = {.1};
//    double lambda_sparse[1] = {1};
    
    int lambda_graph_grid = 5;
    int lambda_sparse_grid = 1;
    double lambda_graph[5] = {1e-2,.1,.2,.4,1.0};
    double lambda_sparse[1] = {1e-2};
    
    double lambda_pairs[10] = {1e-2,1e-2,.1,1e-2,.2,1e-2,.4,1e-2,1.0,1e-2};
    int num_of_lam_pairs = 5;
    
    int max_num_iter = 1e4;
    double* beta_path = (double*) malloc(sizeof(double)*p*lambda_graph_grid*lambda_sparse_grid);
    bool reporting = false;
    return_data* out = (return_data*) malloc(sizeof(return_data));
    out->beta_path = (double*) malloc(sizeof(double)*p*lambda_graph_grid*lambda_sparse_grid);
    out->fun_val = (double*) malloc(sizeof(double)*lambda_graph_grid*lambda_sparse_grid);
    out->chol_num = (int*) malloc(sizeof(int)*lambda_graph_grid*lambda_sparse_grid);
    out->admm_iter_num = (int*) malloc(sizeof(int)*lambda_graph_grid*lambda_sparse_grid);
    if (reporting) out->fun_path = (double*) malloc(sizeof(double)*max_num_iter*lambda_graph_grid*lambda_sparse_grid);
    penalty_types pen_types = graph_sparse_L1;
    varying_rho_types varying_rho_scheme = adaptive;
    algorithm_type algorithm = special_admm;
    double rho = 1.0;
    dmat_vset(p*lambda_graph_grid*lambda_sparse_grid,.0,out->beta_path);
    linreg_path(&p, &m, &n, &pen_types, &varying_rho_scheme, Y, X, M, &theta, lambda_graph, &lambda_graph_grid, lambda_sparse,
                             &lambda_sparse_grid, out, &max_num_iter, &algorithm, reporting, &rho, &eps_abs, &eps_rel);
    // print_dmatrix(out->beta_path, 1, p);
    print_dmatrix(out->fun_val, 1, lambda_graph_grid*lambda_sparse_grid);
    print_imatrix(out->admm_iter_num, 1, lambda_graph_grid*lambda_sparse_grid);
    rho = 1.0;
    double diag_mat = 100.0; // set a number greater than the largest eigenvalue of AtA //
    bool diag_mat_given = true;
    dmat_vset(p*lambda_graph_grid*lambda_sparse_grid,.0,out->beta_path);
    linreg_path_v2(&p, &m, &n, &varying_rho_scheme, Y, X, M, &theta, lambda_pairs, &num_of_lam_pairs, out, &max_num_iter, &algorithm, &diag_mat, diag_mat_given, reporting, &rho, &eps_abs, &eps_rel);
    
    // print_dmatrix(out->beta_path, 1, p);
    print_dmatrix(out->fun_val, 1, lambda_graph_grid*lambda_sparse_grid);
    print_imatrix(out->admm_iter_num, 1, lambda_graph_grid*lambda_sparse_grid);
    
    
    rho = 1.0;
    diag_mat = 100.0; // set a number greater than the largest eigenvalue of AtA //
    diag_mat_given = false;
    algorithm = standard_admm;
    dmat_vset(p*lambda_graph_grid*lambda_sparse_grid,.0,out->beta_path);
    linreg_path_v2(&p, &m, &n, &varying_rho_scheme, Y, X, M, &theta, lambda_pairs, &num_of_lam_pairs, out, &max_num_iter, &algorithm, &diag_mat, diag_mat_given, reporting, &rho, &eps_abs, &eps_rel);
    
    // print_dmatrix(out->beta_path, 1, p);
    print_dmatrix(out->fun_val, 1, lambda_graph_grid*lambda_sparse_grid);
    print_imatrix(out->admm_iter_num, 1, lambda_graph_grid*lambda_sparse_grid);
    
    // if (out->fun_path) free(out->fun_path);
    free(out->admm_iter_num); free(out->chol_num);
    free(out->fun_val); free(out->beta_path);
    free(M);
    free(out);
    free(Y);
    free(beta);
    free(alpha);
    free(beta_initial);
    free(idx); free(jdx); free(val); free(X); free(beta_path);
}


int main(int argc, const char * argv[])
{
    test_sparse_graph_linreg();
    // test_blas();
    // test_fuse_simple();
    // test_linreg();
    // test_UF();
    return 0;
}


