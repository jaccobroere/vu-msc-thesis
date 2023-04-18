//
//  gen_lasso.cpp
//  gen_lasso
//
//  Created by Yunzhang Zhu on 10/10/14.
//  Copyright (c) 2014 OSU. All rights reserved.
//
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>


#include "gen_lasso.h"
#include "dmatrix.h"
#include "def.h"
#include "util.h"

/*
 void graph_update(int m, double* Abeta,
 double* alpha2, double* lambda,
 double* rho,bool* graph){
 double tmp1 = lambda[0] / rho[0];
 double tmp;
 for (int i = 0; i < m; i++){
 tmp = Abeta[i]+alpha2[i];
 graph[i] = (abs(tmp) <= tmp1);
 }
 }*/

void graph_update(linreg_data* data, parameters* para, tmpvars* tmp){
    double a, tmp0, tmp1, tmp2, rho_inv; int i;
    tmp2 = data->lambda_graph[0] / para->rho;
    rho_inv = 1.0 / para->rho;
    if (data->pen_types == graph_sparse_L1)
    {
        tmp1 = data->lambda_sparse[0] / para->rho;
        for (i = 0; i < data->m; i++) {
            a = tmp->Abeta[i] + rho_inv*tmp->alpha2[i];
            tmp0 = (i < data->p)? tmp1 : tmp2;
            tmp->graph[i] = (abs(a) <= tmp0);
        }
    } else
    {
        for (i = 0; i < data->m; i++) {
            a = tmp->Abeta[i] + rho_inv*tmp->alpha2[i];
            tmp->graph[i] = (abs(a) <= tmp2);
        }
    }
}

/* calculate function value */
double fun_val(linreg_data* data, parameters* para, tmpvars* tmp)
{
    // dmat_yAx(data->n, data->p, data->X, tmp->beta, tmp->temp);
    dmat_waxpby(data->n, 1.0, data->Y, -1.0, tmp->Xbeta, tmp->temp);
    double fun_val = dmat_norm2(data->n, tmp->temp);
    fun_val = .5*fun_val*fun_val + data->lambda_sparse[0]*dmat_norm1(data->p, tmp->Abeta)
    + data->lambda_graph[0]*dmat_norm1(data->m - data->p, tmp->Abeta+data->p);
    return (fun_val);
}


void ini_linreg_prob(int* p, int* m,
                     int* n, penalty_types* pen_types,
                     varying_rho_types* varying_rho_scheme,
                     double* Y,
                     double* X, dmatrix_csr* A,
                     double* theta,
                     double* lambda_graph, int* lambda_graph_grid,
                     double* lambda_sparse, int* lambda_sparse_grid,
                     int* num_iter,
                     algorithm_type* algorithm,
                     bool reporting, double* rho,
                     double* eps_abs, double* eps_rel,
                     linreg_data** data_, parameters** para_, tmpvars** tmp_)
{
    linreg_data* data   = (linreg_data*) malloc(sizeof(linreg_data));
    parameters* para    = (parameters*) malloc(sizeof(parameters));
    tmpvars* tmp        = (tmpvars*) malloc(sizeof(tmpvars));
    int i;
    
    data->p             = p[0];
    data->m             = (pen_types[0] == graph_sparse_L1)? (m[0]+p[0]):(m[0]);
    data->n             = n[0];
    data->X             = X;
    data->Y             = Y;
    data->A             = A;
    data->lambda_graph        = lambda_graph;
    data->lambda_graph_grid   = lambda_graph_grid[0];
    data->lambda_sparse       = lambda_sparse;
    data->lambda_sparse_grid  = lambda_sparse_grid[0];
    data->pen_types     = pen_types[0];
    data->varying_rho_scheme = varying_rho_scheme[0];
    // data->constraints   = NULL;
    
    para->theta         = theta;
    para->max_num_iter  = num_iter;
    para->eps_abs       = eps_abs;
    para->eps_rel       = eps_rel;
    para->algorithm     = algorithm[0];
    para->reporting     = reporting;
    para->mu            = .1;
    para->tau           = 2.0;
    if (varying_rho_scheme[0] == constant) { para->rho = rho[0]; }
    if (varying_rho_scheme[0] == adaptive)
    {
        if (data->pen_types == graph_sparse_L1)
        { para->rho = max(1e-10,.5*(lambda_sparse[0]+lambda_graph[0])); }
        else { para->rho = data->lambda_graph[0]; }
    }
    tmp->npm            = max(data->n,max(data->m,data->p));
    tmp->skinny         = false;
    tmp->k              = 0;
    tmp->j              = 0;
    tmp->prires         = .0;
    tmp->dualres        = .0;
    tmp->eps_pri        = .0;
    tmp->eps_dual       = .0;
    tmp->beta           = (double*) malloc(sizeof(double)*data->p);
    tmp->dualres_vec    = (double*) malloc(sizeof(double)*data->p);
    tmp->beta_old       = (double*) malloc(sizeof(double)*data->p);
    tmp->alpha_barA     = (double*) malloc(sizeof(double)*data->p);
    tmp->alpha2A        = (double*) malloc(sizeof(double)*data->p);
    tmp->XTY            = (double*) malloc(sizeof(double)*data->p);
    tmp->alpha          = (double*) malloc(sizeof(double)*data->m);
    tmp->alpha1         = (double*) malloc(sizeof(double)*data->m);
    tmp->alpha2         = (double*) malloc(sizeof(double)*data->m);
    tmp->alpha_bar      = (double*) malloc(sizeof(double)*data->m);
    tmp->graph          = (bool*  ) malloc(sizeof(bool  )*data->m);
    tmp->prires_vec     = (double*) malloc(sizeof(double)*data->m);
    tmp->Abeta          = (double*) malloc(sizeof(double)*data->m);
    tmp->temp           = (double*) malloc(sizeof(double)*tmp->npm);
    tmp->Xbeta          = (double*) malloc(sizeof(double)*data->n);
    tmp->degree         = NULL;
    tmp->diag_mat       = NULL;
    tmp->diag_mat_inv   = NULL;
    tmp->L              = NULL;
    tmp->A_csc          = NULL;
    tmp->varying_rho_points = min(n[0], p[0]);
    tmp->Xtalpha3       = NULL;
    
    dmat_yATx(data->n, data->p, data->X, data->Y, tmp->XTY);
    dmat_vset(data->p, .0, tmp->beta);
    dmat_vset(data->m, .0, tmp->alpha);
    dmat_vset(data->m, .0, tmp->alpha1);
    dmat_vset(data->m, .0, tmp->alpha2);
    dmat_vset(data->m, .0, tmp->alpha_bar);
    
    
    //    if (algorithm[0] == linearized_admm) {
    //        tmp->skinny = (data->n >= data->p);
    //        tmp->degree = (int*) malloc(sizeof(double)*data->p);
    //        dmat_iset(data->p, 0, tmp->degree);
    //        tmp->diag_mat = (double*) malloc(sizeof(double)*data->p);
    //        tmp->diag_mat_inv = (double*) malloc(sizeof(double)*data->p);
    //
    //        if (data->pen_types == general_L1) {
    //            double max_sig_A;
    //            dmatcsr_to_csc(data->A, &(tmp->A_csc));
    //            get_diag_general(data->A, tmp->A_csc, &max_sig_A);
    //            dmat_vset(data->p, max_sig_A, tmp->diag_mat);
    //        } else {
    //            for (i = 0; i < data->A->nz; i++) { tmp->degree[data->A->jdx[i]] += 1; }
    //        }
    //        if (tmp->skinny) { tmp->L = (double*) malloc(sizeof(double)*data->p*data->p); }
    //        else { tmp->L = (double*) malloc(sizeof(double)*data->n*data->n); }
    //    } else if (algorithm[0] == standard_admm) { // standard ADMM //
    //        tmp->L = (double*) malloc(sizeof(double)*data->p*data->p);
    //        dmatcsr_to_csc(data->A, &(tmp->A_csc));
    //    }
    
    if (algorithm[0] == special_admm) {
        tmp->Xtalpha3 = (double*) malloc(sizeof(double)*p[0]);
        dmat_vset(p[0], .0, tmp->Xtalpha3);
    }
    if (algorithm[0] == standard_admm) {
        tmp->L = (double*) malloc(sizeof(double)*data->p*data->p);
        dmatcsr_to_csc(data->A, &(tmp->A_csc));
    } else {
        tmp->skinny = (data->n >= data->p);
        tmp->degree = (int*) malloc(sizeof(double)*data->p);
        dmat_iset(data->p, 0, tmp->degree);
        tmp->diag_mat = (double*) malloc(sizeof(double)*data->p);
        tmp->diag_mat_inv = (double*) malloc(sizeof(double)*data->p);
        
        if (data->pen_types == general_L1) {
            double max_sig_A;
            dmatcsr_to_csc(data->A, &(tmp->A_csc));
            get_diag_general(data->A, tmp->A_csc, &max_sig_A);
            dmat_vset(data->p, max_sig_A, tmp->diag_mat);
        } else {
            for (i = 0; i < data->A->nz; i++) { tmp->degree[data->A->jdx[i]] += 1; }
        }
        if (tmp->skinny) { tmp->L = (double*) malloc(sizeof(double)*data->p*data->p); }
        else { tmp->L = (double*) malloc(sizeof(double)*data->n*data->n); }
    }
    
    *data_  = data;
    *para_  = para;
    *tmp_   = tmp;
}


// initialization for linreg_path_v2(...) //
void ini_linreg_prob_v2(int* p, int* m,
                        int* n,
                        varying_rho_types* varying_rho_scheme,
                        double* Y,
                        double* X, dmatrix_csr* A,
                        double* theta,
                        double* lambda_pairs, int* num_of_pairs,
                        int* num_iter,
                        algorithm_type* algorithm,
                        double* diag_mat, bool diag_mat_given,
                        bool reporting, double* rho,
                        double* eps_abs, double* eps_rel,
                        linreg_data** data_, parameters** para_, tmpvars** tmp_)
{
    linreg_data* data   = (linreg_data*) malloc(sizeof(linreg_data));
    parameters* para    = (parameters*) malloc(sizeof(parameters));
    tmpvars* tmp        = (tmpvars*) malloc(sizeof(tmpvars));
    int i;
    
    data->p             = p[0];
    data->m             = m[0]+p[0];
    data->n             = n[0];
    data->X             = X;
    data->Y             = Y;
    data->A             = A;
    data->lambda_pairs  = lambda_pairs;
    data->lambda_pairs_grid = num_of_pairs[0];
    data->pen_types     = graph_sparse_L1;
    data->varying_rho_scheme = varying_rho_scheme[0];
    // data->constraints   = NULL;
    
    para->theta         = theta;
    para->max_num_iter  = num_iter;
    para->eps_abs       = eps_abs;
    para->eps_rel       = eps_rel;
    para->algorithm     = algorithm[0];
    para->diag_mat_given= diag_mat_given;
    para->reporting     = reporting;
    para->mu            = .1;
    para->tau           = 2.0;
    if (varying_rho_scheme[0] == constant) { para->rho = rho[0]; }
    if (varying_rho_scheme[0] == adaptive)
    {
        para->rho = max(1e-10,.5*(lambda_pairs[0]+lambda_pairs[1]));
    }
    tmp->npm            = max(data->n,max(data->m,data->p));
    tmp->skinny         = false;
    tmp->k              = 0;
    tmp->j              = 0;
    tmp->prires         = .0;
    tmp->dualres        = .0;
    tmp->eps_pri        = .0;
    tmp->eps_dual       = .0;
    tmp->beta           = (double*) malloc(sizeof(double)*data->p);
    tmp->dualres_vec    = (double*) malloc(sizeof(double)*data->p);
    tmp->beta_old       = (double*) malloc(sizeof(double)*data->p);
    tmp->alpha_barA     = (double*) malloc(sizeof(double)*data->p);
    tmp->alpha2A        = (double*) malloc(sizeof(double)*data->p);
    tmp->XTY            = (double*) malloc(sizeof(double)*data->p);
    tmp->alpha          = (double*) malloc(sizeof(double)*data->m);
    tmp->alpha1         = (double*) malloc(sizeof(double)*data->m);
    tmp->alpha2         = (double*) malloc(sizeof(double)*data->m);
    tmp->alpha_bar      = (double*) malloc(sizeof(double)*data->m);
    tmp->graph          = (bool*  ) malloc(sizeof(bool  )*data->m);
    tmp->prires_vec     = (double*) malloc(sizeof(double)*data->m);
    tmp->Abeta          = (double*) malloc(sizeof(double)*data->m);
    tmp->temp           = (double*) malloc(sizeof(double)*tmp->npm);
    tmp->Xbeta          = (double*) malloc(sizeof(double)*data->n);
    tmp->degree         = NULL;
    tmp->diag_mat       = NULL;
    tmp->diag_mat_inv   = NULL;
    tmp->L              = NULL;
    tmp->A_csc          = NULL;
    tmp->varying_rho_points = min(n[0], p[0]);
    tmp->Xtalpha3       = NULL;
    
    dmat_yATx(data->n, data->p, data->X, data->Y, tmp->XTY);
    dmat_vset(data->p, .0, tmp->beta);
    dmat_vset(data->m, .0, tmp->alpha);
    dmat_vset(data->m, .0, tmp->alpha1);
    dmat_vset(data->m, .0, tmp->alpha2);
    dmat_vset(data->m, .0, tmp->alpha_bar);
    
    
    //    if (algorithm[0] == linearized_admm) {
    //        tmp->skinny = (data->n >= data->p);
    //        tmp->degree = (int*) malloc(sizeof(double)*data->p);
    //        dmat_iset(data->p, 0, tmp->degree);
    //        tmp->diag_mat = (double*) malloc(sizeof(double)*data->p);
    //        tmp->diag_mat_inv = (double*) malloc(sizeof(double)*data->p);
    //
    //        if (data->pen_types == general_L1) {
    //            double max_sig_A;
    //            dmatcsr_to_csc(data->A, &(tmp->A_csc));
    //            get_diag_general(data->A, tmp->A_csc, &max_sig_A);
    //            dmat_vset(data->p, max_sig_A, tmp->diag_mat);
    //        } else {
    //            for (i = 0; i < data->A->nz; i++) { tmp->degree[data->A->jdx[i]] += 1; }
    //        }
    //        if (tmp->skinny) { tmp->L = (double*) malloc(sizeof(double)*data->p*data->p); }
    //        else { tmp->L = (double*) malloc(sizeof(double)*data->n*data->n); }
    //    } else if (algorithm[0] == standard_admm) { // standard ADMM //
    //        tmp->L = (double*) malloc(sizeof(double)*data->p*data->p);
    //        dmatcsr_to_csc(data->A, &(tmp->A_csc));
    //    }
    
    if (algorithm[0] == special_admm) {
        tmp->Xtalpha3 = (double*) malloc(sizeof(double)*p[0]);
        dmat_vset(p[0], .0, tmp->Xtalpha3);
    }
    if (algorithm[0] == standard_admm) {
        tmp->L = (double*) malloc(sizeof(double)*data->p*data->p);
        dmatcsr_to_csc(data->A, &(tmp->A_csc));
    } else {
        tmp->skinny = (data->n >= data->p);
        tmp->degree = (int*) malloc(sizeof(double)*data->p);
        dmat_iset(data->p, 0, tmp->degree);
        tmp->diag_mat = (double*) malloc(sizeof(double)*data->p);
        tmp->diag_mat_inv = (double*) malloc(sizeof(double)*data->p);
        
        
        if (diag_mat_given) // given D //
        { dmat_vset(data->p, diag_mat[0], tmp->diag_mat); }
        else if (data->pen_types == general_L1) // for general A, set D to be max. singular value //
        {
            double max_sig_A;
            dmatcsr_to_csc(data->A, &(tmp->A_csc));
            get_diag_general(data->A, tmp->A_csc, &max_sig_A);
            dmat_vset(data->p, max_sig_A, tmp->diag_mat);
        } else // for oriented incidence matrix A, set a special D //
        {
            for (i = 0; i < data->A->nz; i++) { tmp->degree[data->A->jdx[i]] += 1; }
        }
        if (tmp->skinny) { tmp->L = (double*) malloc(sizeof(double)*data->p*data->p); }
        else { tmp->L = (double*) malloc(sizeof(double)*data->n*data->n); }
    }
    
    *data_  = data;
    *para_  = para;
    *tmp_   = tmp;
}


void free_tmpvars_linreg(tmpvars* tmp)
{
    if (tmp) {
        if (tmp->prires_vec)    { free(tmp->prires_vec);    }
        if (tmp->dualres_vec)   { free(tmp->dualres_vec);   }
        if (tmp->beta)          { free(tmp->beta);          }
        if (tmp->beta_old)      { free(tmp->beta_old);      }
        if (tmp->alpha)         { free(tmp->alpha);         }
        if (tmp->alpha1)        { free(tmp->alpha1);        }
        if (tmp->alpha2)        { free(tmp->alpha2);        }
        if (tmp->Xtalpha3)      { free(tmp->Xtalpha3);      }
        if (tmp->Abeta)         { free(tmp->Abeta);         }
        if (tmp->alpha_barA)    { free(tmp->alpha_barA);    }
        if (tmp->alpha2A)       { free(tmp->alpha2A);       }
        if (tmp->diag_mat)      { free(tmp->diag_mat);      }
        if (tmp->degree)        { free(tmp->degree);        }
        if (tmp->diag_mat_inv)  { free(tmp->diag_mat_inv);  }
        if (tmp->L)             { free(tmp->L);             }
        if (tmp->XTY)           { free(tmp->XTY);           }
        if (tmp->temp)          { free(tmp->temp);          }
        if (tmp->graph)         { free(tmp->graph);         }
        if (tmp->A_csc)         { dmat_csc_free(tmp->A_csc);}
        if (tmp->Xbeta)         { free(tmp->Xbeta);         }
        free(tmp);
    }
    
}

void factorize(int n, int p, int m, penalty_types pen_types,
               algorithm_type algorithm, bool skinny, int* degree,
               double* diag_mat, double* diag_mat_inv,
               dmatrix_csr* A, dmatrix_csc* A_csc,
               double* X, double* alpha, double* beta_old,
               double* rho, double* lambda_graph,
               double* lambda_sparse, double* L)
{
    int i,k;
    //    if (variant) {
    //        // diag_mat = rho*D //
    //        if (pen_types == general_L1) {
    //            for (i = 0; i < p; i++) { diag_mat[i] = rho[0]*diag_mat[i]; }
    //        }
    //        if (pen_types == graph_L1) {
    //            for (i = 0; i < p; i++) { diag_mat[i] = 2*rho[0]*degree[i]; }
    //        }
    //        if (pen_types == graph_sparse_L1) {
    //            for (i = 0; i < p; i++)
    //            { diag_mat[i] = rho[0]*(2*degree[i]+1); }
    //        }
    //        dmat_yinvx(p, diag_mat, diag_mat_inv);
    //        if (skinny) {
    //            dmat_B_ATA(n, p, X, L);
    //            for (i = 0; i < p; i++) {
    //                L[i*p+i] += diag_mat[i];
    //            }
    //            /* L = chol(rho * D + Xt*X) */
    //            dmat_potrf(L,p);
    //        } else {
    //            dmat_B_ADAT(n, p, X, diag_mat_inv, L);
    //            for (i = 0; i < n; i++) {
    //                L[i*n+i] += 1.0;
    //            }
    //            /* L = chol(I + rho^{-1} X D^{-1} Xt) */
    //            dmat_potrf(L,n);
    //        }
    //    } else /* standard ADMM */{
    //        dmat_B_ATA(n, p, X, L);
    //        for (i = 0; i < p; i++) {
    //            dmat_vset(m, .0, alpha); // store i-th column of A
    //            for (k = A_csc->cdx[i]; k < A_csc->cdx[i+1]; k++) {
    //                alpha[A_csc->idx[k]] = A_csc->val[k];
    //            }
    //            dmat_ySTx(A, alpha, beta_old);
    //            dmat_waxpby(p, 1.0, L+i*p, rho[0], beta_old, L+i*p);
    //            if (pen_types == graph_sparse_L1) { L[i*p+i] += rho[0]; }
    //        }
    //        /* L = chol(rho * At*A + Xt*X) */
    //        dmat_potrf(L,p);
    //    }
    
    if (algorithm == standard_admm) {
        dmat_B_ATA(n, p, X, L);
        for (i = 0; i < p; i++) {
            dmat_vset(m, .0, alpha); // store i-th column of A
            for (k = A_csc->cdx[i]; k < A_csc->cdx[i+1]; k++) {
                alpha[A_csc->idx[k]] = A_csc->val[k];
            }
            dmat_ySTx(A, alpha, beta_old);
            dmat_waxpby(p, 1.0, L+i*p, rho[0], beta_old, L+i*p);
            if (pen_types == graph_sparse_L1) { L[i*p+i] += rho[0]; }
        }
        /* L = chol(rho * At*A + Xt*X) */
        dmat_potrf(L,p);
    } else {
        if (pen_types == general_L1) {
            for (i = 0; i < p; i++) { diag_mat[i] = diag_mat[i]; }
        }
        if (pen_types == graph_L1) {
            for (i = 0; i < p; i++) { diag_mat[i] = 2*degree[i]; }
        }
        if (pen_types == graph_sparse_L1) {
            for (i = 0; i < p; i++)
            { diag_mat[i] = 2*degree[i]+1; }
        }
        if (algorithm == linearized_admm) { for (i = 0; i < p; i++) { diag_mat[i] *= rho[0]; } }
        dmat_yinvx(p, diag_mat, diag_mat_inv);
        if (skinny) {
            dmat_B_ATA(n, p, X, L);
            for (i = 0; i < p; i++) {
                L[i*p+i] += diag_mat[i];
            }
            /* L = chol(rho * D + Xt*X) */
            dmat_potrf(L,p);
        } else {
            dmat_B_ADAT(n, p, X, diag_mat_inv, L);
            for (i = 0; i < n; i++) {
                L[i*n+i] += 1.0;
            }
            /* L = chol(I + rho^{-1} X D^{-1} Xt) */
            dmat_potrf(L,n);
        }
    }
}

void factorize_v2(int n, int p, int m, penalty_types pen_types,
               algorithm_type algorithm, bool diag_mat_given,
                  bool skinny, int* degree,
               double* diag_mat, double* diag_mat_inv,
               dmatrix_csr* A, dmatrix_csc* A_csc,
               double* X, double* alpha, double* beta_old,
               double* rho, double* lambda_graph,
               double* lambda_sparse, double* L)
{
    int i,k;
    
    if (algorithm == standard_admm) {
        dmat_B_ATA(n, p, X, L);
        for (i = 0; i < p; i++) {
            dmat_vset(m, .0, alpha); // store i-th column of A
            for (k = A_csc->cdx[i]; k < A_csc->cdx[i+1]; k++) {
                alpha[A_csc->idx[k]] = A_csc->val[k];
            }
            dmat_ySTx(A, alpha, beta_old);
            dmat_waxpby(p, 1.0, L+i*p, rho[0], beta_old, L+i*p);
            if (pen_types == graph_sparse_L1) { L[i*p+i] += rho[0]; }
        }
        /* L = chol(rho * At*A + Xt*X) */
        dmat_potrf(L,p);
    } else {
        if (!diag_mat_given) {
            if (pen_types == general_L1) {
                for (i = 0; i < p; i++) { diag_mat[i] = diag_mat[i]; }
            }
            if (pen_types == graph_L1) {
                for (i = 0; i < p; i++) { diag_mat[i] = 2*degree[i]; }
            }
            if (pen_types == graph_sparse_L1) {
                for (i = 0; i < p; i++)
                { diag_mat[i] = 2*degree[i]+1; }
            }
            if (algorithm == linearized_admm) { for (i = 0; i < p; i++) { diag_mat[i] *= rho[0]; } }
        }
        dmat_yinvx(p, diag_mat, diag_mat_inv);
        if (skinny) {
            dmat_B_ATA(n, p, X, L);
            for (i = 0; i < p; i++) {
                L[i*p+i] += diag_mat[i];
            }
            /* L = chol(rho * D + Xt*X) */
            dmat_potrf(L,p);
        } else {
            dmat_B_ADAT(n, p, X, diag_mat_inv, L);
            for (i = 0; i < n; i++) {
                L[i*n+i] += 1.0;
            }
            /* L = chol(I + rho^{-1} X D^{-1} Xt) */
            dmat_potrf(L,n);
        }
    }
}


void beta_update_linreg(int p, int n, double* XTY, double* X,
                        double* alpha_bar, double* alpha,
                        double* alpha_barA, double* beta,
                        double* diag_mat, double* diag_mat_inv,
                        double* rho, dmatrix_csr* A,
                        double* L, double* tmp,
                        double* dualres_vec, double* Xtalpha3,
                        algorithm_type algorithm, bool skinny,
                        penalty_types pen_types){
    if (pen_types == graph_sparse_L1)
    {
        dmat_ySTx(A, alpha_bar+p, alpha_barA); // At * alpha_bar2
        dmat_waxpby(p, 1.0, alpha_bar, 1.0, alpha_barA, alpha_barA); // alpha_bar1 + At * alpha_bar2
        
    } else { dmat_ySTx(A, alpha_bar, alpha_barA); } // At * alpha
    
    dmat_waxpby(p, 1.0, XTY, -1.0, alpha_barA, tmp);
    if (algorithm == special_admm) {
        dmat_waxpby(p, 1.0, tmp, rho[0] - 1.0, Xtalpha3, tmp);
    }
    
    if (algorithm == standard_admm) {
        if (pen_types == graph_sparse_L1) { dmat_waxpby(p, rho[0], beta, 1.0, tmp, tmp); } // tmp = rho*beta + XtY - alpha_barA
        dmat_ySx(A, beta, alpha);
        dmat_ySTx(A, alpha, beta);
        dmat_waxpby(p, rho[0], beta, 1.0, tmp, beta); // beta = rho*AtA*beta + XtY - alpha_barA //
        dmat_potrs(L, p, beta);
    } else {
        dmat_elemprod(p, diag_mat, beta, beta);
        if (algorithm == special_admm) {
            dmat_waxpby(p, 1.0, beta, 1/rho[0], tmp, beta);
        } else {
            dmat_waxpby(p, 1.0, beta, 1.0, tmp, beta);
        }
        if (skinny) {
            dmat_potrs(L, p, beta);
        } else {
            dmat_elemprod(p, diag_mat_inv, beta, beta); // beta = beta ./ (rho*D)
            dmat_vcopy(p, beta, dualres_vec); // cache = beta
            dmat_yAx(n, p, X, beta, tmp); // tmp = X * beta
            dmat_potrs(L, n, tmp); // tmp = (I + rho^{-1} X D^{-1} Xt)^{-1} * tmp
            dmat_yATx(n, p, X, tmp, beta); // beta = X^T * beta
            dmat_elemprod(p, diag_mat_inv, beta, beta); // beta = beta ./ (rho*D)
            dmat_waxpby(p, 1.0, dualres_vec, -1.0, beta, beta); // beta = cache - beta
        }
    }
}

void alpha_update_linreg(int m, int n, double* X,
                         dmatrix_csr* A,
                         double* beta, double* Abeta,
                         double* alpha1, double* alpha2,
                         double* rho, double* lambda_sparse,
                         double* lambda_graph,
                         penalty_types pen_types,
                         algorithm_type algorithm,
                         double* Xtalpha3, double* XTY,
                         double* temp, double* Xbeta)
{
    double tmp; int i, p = A->n;
    dmat_yAx(n, p, X, beta, Xbeta);
    dmat_vcopy(m,alpha2,alpha1);
    if (pen_types == graph_sparse_L1)
    { // with sparse penalty //
        dmat_vcopy(p, beta, Abeta);
        dmat_ySx(A, beta, Abeta+p);
        double lambda;
        for (i = 0; i < m; i++) {
            lambda = (i < p)? lambda_sparse[0] : lambda_graph[0];
            tmp = rho[0]*Abeta[i]+alpha2[i];
            alpha2[i] = (abs(tmp) <= lambda)? tmp : (lambda*sign(tmp));
        }
    } else // without sparse penalty //
    {
        dmat_ySx(A, beta, Abeta);
        for (i = p; i < m; i++){
            tmp = rho[0]*Abeta[i]+alpha2[i];
            alpha2[i] = (abs(tmp) <= lambda_graph[0])? tmp : (lambda_graph[0]*sign(tmp));
        }
    }
    
    if (algorithm == special_admm) {
        dmat_yATx(n, p, X, Xbeta, temp);
        dmat_waxpby(p, 1.0/(1.0+rho[0]), Xtalpha3, rho[0]/(1.0+rho[0]), temp, Xtalpha3);
    }
}

void stopping_para_update(int m, int n, int p, double* beta,
                          double* alpha1, double* alpha2,
                          double* prires_vec, double* rho,
                          double* dualres_vec, double* beta_old,
                          double* diag_mat, dmatrix_csr* A,
                          double* alpha_barA, double* alpha2A,
                          double* alpha, double* eps_abs,
                          double* eps_rel, double* Abeta,
                          penalty_types pen_types,
                          tmpvars* tmp)
{
    // calculate primal residuals //
    dmat_waxpby(m, 1.0, alpha2, -1.0, alpha1, prires_vec);
    tmp->prires = dmat_norm2(m, prires_vec) / rho[0];
    
    // calculate dual residuals //
    dmat_waxpby(p, 1.0, beta, -1.0, beta_old, alpha2A);
    dmat_ySx(A, alpha2A, alpha);
    dmat_ySTx(A, alpha, dualres_vec);
    if (pen_types == graph_sparse_L1) { dmat_waxpby(p, 1.0, dualres_vec, 1.0, alpha2A, dualres_vec); }
    dmat_waxpby(p, 1.0, alpha_barA, rho[0], dualres_vec, dualres_vec);
    if (pen_types == graph_sparse_L1) {
        dmat_ySTx(A, alpha2+p, alpha2A);
        dmat_waxpby(p, 1.0, alpha2, 1.0, alpha2A, alpha2A);
    } else { dmat_ySTx(A, alpha2, alpha2A); }
    dmat_waxpby(p, -1.0, alpha2A, 1.0, dualres_vec, dualres_vec);
    tmp->dualres = dmat_norm2(p, dualres_vec);
    
    // calculate pri_abs and dual_abs //
    tmp->eps_pri = sqrt(m)*eps_abs[0]  + eps_rel[0] * dmat_norm2(m, Abeta);
    tmp->eps_dual = sqrt(p)*eps_abs[0] + dmat_norm2(p, alpha2A)*eps_rel[0];
}




void linreg_sub(linreg_data* data, return_data* out, parameters* para, tmpvars* tmp)
{
    int i, s = 1;
    // printf("%3s %10s %10s %10s %10s %10s\n", "#", "r norm", "eps_pri", "s norm", "eps_dual", "objective");
    
    if (data->varying_rho_scheme == non_adaptive) // if varying_rho_scheme = non_adaptive, set rho = lambda, //
    {
        if (data->pen_types == graph_sparse_L1)
        { para->rho = max(1e-10,.5*(data->lambda_sparse[0]+data->lambda_graph[0])); }
        else { para->rho = data->lambda_graph[0]; }
        
        factorize(data->n, data->p, data->m, data->pen_types, para->algorithm,
                  tmp->skinny, tmp->degree, tmp->diag_mat, tmp->diag_mat_inv,
                  data->A, tmp->A_csc, data->X, tmp->alpha, tmp->beta_old,
                  &para->rho, data->lambda_graph, data->lambda_sparse,tmp->L);
        (out->chol_num[data->lambda_sparse_grid*tmp->k+tmp->j])++;
    }
    // printf("rho is: %f. \n", para->rho);
    for (i = 0; i < para->max_num_iter[0]; i++) {
        dmat_vcopy(data->p, tmp->beta, tmp->beta_old); // save beta
        // update beta //
        beta_update_linreg(data->p, data->n,
                           tmp->XTY, data->X,
                           tmp->alpha_bar,tmp->alpha,
                           tmp->alpha_barA, tmp->beta,
                           tmp->diag_mat, tmp->diag_mat_inv,
                           &para->rho, data->A, tmp->L,
                           tmp->temp, tmp->dualres_vec,
                           tmp->Xtalpha3,
                           para->algorithm, tmp->skinny,
                           data->pen_types);
        
        // update alpha //
        alpha_update_linreg(data->m,data->n, data->X,
                            data->A,tmp->beta,
                            tmp->Abeta, tmp->alpha1, tmp->alpha2,
                            &para->rho,data->lambda_sparse,
                            data->lambda_graph,data->pen_types,
                            para->algorithm, tmp->Xtalpha3,
                            tmp->XTY, tmp->temp,tmp->Xbeta);
        
        // alpha_bar = alpha2 + theta*(alpha2 - alpha1)
        dmat_waxpby(data->m, 1.0+para->theta[0], tmp->alpha2,
                    -para->theta[0], tmp->alpha1, tmp->alpha_bar);
        
        stopping_para_update(data->m, data->n, data->p, tmp->beta,
                             tmp->alpha1, tmp->alpha2,
                             tmp->prires_vec, &para->rho,
                             tmp->dualres_vec, tmp->beta_old,
                             tmp->diag_mat, data->A,
                             tmp->alpha_barA, tmp->alpha2A,
                             tmp->alpha, para->eps_abs,
                             para->eps_rel, tmp->Abeta,
                             data->pen_types, tmp);
        
        //        printf("%4d %10.6f %10.6f %10.6f %10.6f %10.6f\n",
        //               i, tmp->prires, tmp->eps_pri, tmp->dualres,
        //               tmp->eps_dual, fun_val(data,para,tmp));
        
        if (data->varying_rho_scheme == adaptive && (i == ((s+1)*s*tmp->varying_rho_points/2)))
            // only change rho if i = .5(s+1)s min(n,p) //
        {
            s++;
            if ((tmp->prires / tmp->eps_pri) < para->mu*(tmp->dualres / tmp->eps_dual)) {
                (out->chol_num[data->lambda_sparse_grid*tmp->k+tmp->j])++;
                para->rho /= para->tau;
                // printf("rho / 2 \n");
                if (para->algorithm != special_admm) {
                    factorize(data->n, data->p, data->m, data->pen_types, para->algorithm,
                              tmp->skinny, tmp->degree, tmp->diag_mat, tmp->diag_mat_inv,
                              data->A, tmp->A_csc, data->X, tmp->alpha, tmp->beta_old,
                              &para->rho, data->lambda_graph, data->lambda_sparse,tmp->L);
                }
            }
            if ((tmp->dualres / tmp->eps_dual) < para->mu*(tmp->prires / tmp->eps_pri)) {
                (out->chol_num[data->lambda_sparse_grid*tmp->k+tmp->j])++;
                para->rho *= para->tau;
                // printf("rho * 2 \n");
                if (para->algorithm != special_admm) {
                    factorize(data->n, data->p, data->m, data->pen_types, para->algorithm,
                              tmp->skinny, tmp->degree, tmp->diag_mat, tmp->diag_mat_inv,
                              data->A, tmp->A_csc, data->X, tmp->alpha, tmp->beta_old,
                              &para->rho, data->lambda_graph, data->lambda_sparse,tmp->L);
                }
            }
        }
        
        if (para->reporting) // report function values and primal residuals //
        {
            out->fun_path[para->max_num_iter[0]*(data->lambda_sparse_grid*tmp->k+tmp->j)+i] = fun_val(data, para, tmp);
            continue; // when reporting, we let it run until reaching max iter. //
        }
        
        if ((tmp->prires < tmp->eps_pri) && (tmp->dualres < tmp->eps_dual)) // check convergence //
        {
            out->fun_val[data->lambda_sparse_grid*tmp->k+tmp->j] = fun_val(data,para,tmp);
            out->admm_iter_num[data->lambda_sparse_grid*tmp->k+tmp->j] = i;
//                        printf("#iter:%d, #chol:%d, rho:%f, snorm:%10.6f, fun_val:%10.6f. \n",
//                               i, out->chol_num[data->lambda_sparse_grid*tmp->k+tmp->j],
//                               para->rho,tmp->prires,out->fun_val[data->lambda_sparse_grid*tmp->k+tmp->j]);
            break;
        }
        
        if (i == (para->max_num_iter[0]-1)) {
            out->fun_val[data->lambda_sparse_grid*tmp->k+tmp->j] = fun_val(data,para,tmp);
            out->admm_iter_num[data->lambda_sparse_grid*tmp->k+tmp->j] = i;
//                        printf("#iter:>%d, #chol:%d, rho:%f, snorm:%10.6f, fun_val:%10.6f. \n",
//                               i, out->chol_num[data->lambda_sparse_grid*tmp->k+tmp->j],
//                               para->rho,tmp->prires,out->fun_val[data->lambda_sparse_grid*tmp->k+tmp->j]);
        }
    }
}



void linreg_sub_v2(linreg_data* data, return_data* out, parameters* para, tmpvars* tmp)
{
    int i, s = 1;
    // printf("%3s %10s %10s %10s %10s %10s\n", "#", "r norm", "eps_pri", "s norm", "eps_dual", "objective");
    
    if (data->varying_rho_scheme == non_adaptive) // if varying_rho_scheme = non_adaptive, set rho = lambda, //
    {
        if (data->pen_types == graph_sparse_L1)
        { para->rho = max(1e-10,.5*(data->lambda_sparse[0]+data->lambda_graph[0])); }
        else { para->rho = data->lambda_graph[0]; }
        
        factorize_v2(data->n, data->p, data->m, data->pen_types, para->algorithm,
                  para->diag_mat_given,
                  tmp->skinny, tmp->degree, tmp->diag_mat, tmp->diag_mat_inv,
                  data->A, tmp->A_csc, data->X, tmp->alpha, tmp->beta_old,
                  &para->rho, data->lambda_graph, data->lambda_sparse,tmp->L);
        (out->chol_num[tmp->k])++;
    }
    // printf("rho is: %f. \n", para->rho);
    for (i = 0; i < para->max_num_iter[0]; i++) {
        dmat_vcopy(data->p, tmp->beta, tmp->beta_old); // save beta
        // update beta //
        beta_update_linreg(data->p, data->n,
                           tmp->XTY, data->X,
                           tmp->alpha_bar,tmp->alpha,
                           tmp->alpha_barA, tmp->beta,
                           tmp->diag_mat, tmp->diag_mat_inv,
                           &para->rho, data->A, tmp->L,
                           tmp->temp, tmp->dualres_vec,
                           tmp->Xtalpha3,
                           para->algorithm, tmp->skinny,
                           data->pen_types);
        
        // update alpha //
        alpha_update_linreg(data->m,data->n, data->X,
                            data->A,tmp->beta,
                            tmp->Abeta, tmp->alpha1, tmp->alpha2,
                            &para->rho,data->lambda_sparse,
                            data->lambda_graph,data->pen_types,
                            para->algorithm, tmp->Xtalpha3,
                            tmp->XTY, tmp->temp,tmp->Xbeta);
        
        // alpha_bar = alpha2 + theta*(alpha2 - alpha1)
        dmat_waxpby(data->m, 1.0+para->theta[0], tmp->alpha2,
                    -para->theta[0], tmp->alpha1, tmp->alpha_bar);
        
        stopping_para_update(data->m, data->n, data->p, tmp->beta,
                             tmp->alpha1, tmp->alpha2,
                             tmp->prires_vec, &para->rho,
                             tmp->dualres_vec, tmp->beta_old,
                             tmp->diag_mat, data->A,
                             tmp->alpha_barA, tmp->alpha2A,
                             tmp->alpha, para->eps_abs,
                             para->eps_rel, tmp->Abeta,
                             data->pen_types, tmp);
        
        //        printf("%4d %10.6f %10.6f %10.6f %10.6f %10.6f\n",
        //               i, tmp->prires, tmp->eps_pri, tmp->dualres,
        //               tmp->eps_dual, fun_val(data,para,tmp));
        
        if (data->varying_rho_scheme == adaptive && (i == ((s+1)*s*tmp->varying_rho_points/2)))
            // only change rho if i = .5(s+1)s min(n,p) //
        {
            s++;
            if ((tmp->prires / tmp->eps_pri) < para->mu*(tmp->dualres / tmp->eps_dual)) {
                (out->chol_num[tmp->k])++;
                para->rho /= para->tau;
                // printf("rho / 2 \n");
                if (para->algorithm != special_admm) {
                    factorize_v2(data->n, data->p, data->m, data->pen_types, para->algorithm,
                                 para->diag_mat_given,
                              tmp->skinny, tmp->degree, tmp->diag_mat, tmp->diag_mat_inv,
                              data->A, tmp->A_csc, data->X, tmp->alpha, tmp->beta_old,
                              &para->rho, data->lambda_graph, data->lambda_sparse,tmp->L);
                }
            }
            if ((tmp->dualres / tmp->eps_dual) < para->mu*(tmp->prires / tmp->eps_pri)) {
                (out->chol_num[tmp->k])++;
                para->rho *= para->tau;
                // printf("rho * 2 \n");
                if (para->algorithm != special_admm) {
                    factorize_v2(data->n, data->p, data->m, data->pen_types, para->algorithm,
                                 para->diag_mat_given,
                              tmp->skinny, tmp->degree, tmp->diag_mat, tmp->diag_mat_inv,
                              data->A, tmp->A_csc, data->X, tmp->alpha, tmp->beta_old,
                              &para->rho, data->lambda_graph, data->lambda_sparse,tmp->L);
                }
            }
        }
        
        if (para->reporting) // report function values and primal residuals //
        {
            out->fun_path[para->max_num_iter[0]*tmp->k+i] = fun_val(data, para, tmp);
            continue; // when reporting, we let it run until reaching max iter. //
        }
        
        if ((tmp->prires < tmp->eps_pri) && (tmp->dualres < tmp->eps_dual)) // check convergence //
        {
            out->fun_val[tmp->k] = fun_val(data,para,tmp);
            out->admm_iter_num[tmp->k] = i;
//            printf("#iter:%d, #chol:%d, rho:%f, snorm:%10.6f, fun_val:%10.6f. \n",
//                               i, out->chol_num[tmp->k],
//                               para->rho,tmp->prires,out->fun_val[tmp->k]);
            break;
        }
        
        if (i == (para->max_num_iter[0]-1)) {
            out->fun_val[tmp->k] = fun_val(data,para,tmp);
            out->admm_iter_num[tmp->k] = i;
//            printf("#iter:>%d, #chol:%d, rho:%f, snorm:%10.6f, fun_val:%10.6f. \n",
//                               i, out->chol_num[tmp->k],
//                               para->rho,tmp->prires,out->fun_val[tmp->k]);
        }
    }
}


// path following for lambda1 and lambda2 //
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
                 double* eps_abs, double* eps_rel)
{
    int j,k;
    linreg_data* data;
    parameters* para;
    tmpvars* tmp;
    ini_linreg_prob(beta_dim, A_row_dim, sample_size, pen_types,
                    varying_rho_scheme, Y, X, A, theta,
                    lambda_graph, lambda_graph_grid, lambda_sparse,
                    lambda_sparse_grid, num_iter, algorithm, reporting,
                    rho, eps_abs, eps_rel, &data, &para, &tmp);
    
    if (data->varying_rho_scheme != non_adaptive)
    {
        factorize(data->n, data->p, data->m, data->pen_types, para->algorithm,
                  tmp->skinny, tmp->degree, tmp->diag_mat, tmp->diag_mat_inv,
                  data->A, tmp->A_csc, data->X, tmp->alpha, tmp->beta_old,
                  &para->rho, data->lambda_graph, data->lambda_sparse, tmp->L);
        out->chol_num[0]++;
    }
    
    if (pen_types[0] == graph_sparse_L1) {
        dmat_iset(lambda_graph_grid[0]*lambda_sparse_grid[0], 0, out->chol_num);
        if (data->varying_rho_scheme == adaptive) { (out->chol_num[0])++; }
        for (k = 0; k < lambda_graph_grid[0]; k++) {
            data->lambda_graph = lambda_graph + k;
            tmp->k = k;
            for (j = 0; j < lambda_sparse_grid[0]; j++) {
                // printf("\nlambda_sparse, lambda_graph are : %f, %f. \n", lambda_sparse[j],lambda_graph[k]);
                data->lambda_sparse = lambda_sparse + j;
                tmp->j = j;
                
                linreg_sub(data, out, para, tmp);
                
                if (!reporting) {
                    graph_update(data,para,tmp);
                    UF TMP(data->p,tmp->beta);
                    for (int i = data->p; i < data->m; i++)
                    { if (tmp->graph[i])  TMP.unite(A->jdx[2*(i-data->p)], A->jdx[2*(i-data->p)+1]);}
                    for (int i = 0; i < data->p; i++)
                    { tmp->beta[i] = TMP.mean[TMP.root(TMP.id[i])]; }
                    for (int i = 0; i < data->p; i++) { if (tmp->graph[i]) tmp->beta[i] = .0; }
                }
                
                // print_dmatrix(tmp->beta, 1, data->p);
                /* save current solution */
                dmat_vcopy(data->p,tmp->beta,out->beta_path + j*data->p + k*lambda_sparse_grid[0]*data->p);
            }
        }
    }
    
    if (pen_types[0] == graph_L1) {
        dmat_iset(lambda_graph_grid[0], 0, out->chol_num);
        if (data->varying_rho_scheme == adaptive) { (out->chol_num[0])++; }
        for (int k = 0; k < lambda_graph_grid[0]; k++) {
            // printf("\nlambda is : %f. \n", lambda_graph[k]);
            data->lambda_graph = lambda_graph + k;
            tmp->k = k;
            linreg_sub(data, out, para, tmp);
            
            // graph_update(data->m,tmp->Abeta,tmp->alpha2,data->lambda_graph,para->rho,tmp->graph);
            graph_update(data,para,tmp);
            UF TMP(data->p,tmp->beta);
            for (int i =data->p; i < data->m; i++)
            { if (tmp->graph[i])  TMP.unite(A->jdx[2*(i-data->p)], A->jdx[2*(i-data->p)+1]);}
            for (int i = 0; i < data->p; i++)
            { tmp->beta[i] = TMP.mean[TMP.root(TMP.id[i])]; }
            
            /* save current solution */
            dmat_vcopy(data->p,tmp->beta,out->beta_path+k*data->p);
        }
    }
    
    free_tmpvars_linreg(tmp);
    free(data);
    free(para);
}


// path following for a sequence of pairs of (lambda1, lambda2) //
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
                    double* eps_abs, double* eps_rel)
{
    int k;
    linreg_data* data;
    parameters* para;
    tmpvars* tmp;
    ini_linreg_prob_v2(beta_dim, A_row_dim, sample_size,
                       varying_rho_scheme, Y, X, A, theta,
                       lambda_pairs, num_of_pairs, num_iter,
                       algorithm, diag_mat, diag_mat_given,reporting,
                       rho, eps_abs, eps_rel, &data, &para, &tmp);
    
    if (data->varying_rho_scheme != non_adaptive)
    {
        factorize_v2(data->n, data->p, data->m, data->pen_types, para->algorithm,
                     para->diag_mat_given,
                  tmp->skinny, tmp->degree, tmp->diag_mat, tmp->diag_mat_inv,
                  data->A, tmp->A_csc, data->X, tmp->alpha, tmp->beta_old,
                  &para->rho, data->lambda_graph, data->lambda_sparse, tmp->L);
        out->chol_num[0]++;
    }
    
    
    dmat_iset(num_of_pairs[0], 0, out->chol_num);
    if (data->varying_rho_scheme == adaptive) { (out->chol_num[0])++; }
    for (k = 0; k < num_of_pairs[0]; k++) {
        data->lambda_graph = lambda_pairs + 2*k;
        data->lambda_sparse = lambda_pairs + 2*k + 1;
        tmp->k = k;
        
        linreg_sub_v2(data, out, para, tmp);
        
        if (!reporting) {
            graph_update(data,para,tmp);
            UF TMP(data->p,tmp->beta);
            for (int i = data->p; i < data->m; i++)
            { if (tmp->graph[i])  TMP.unite(A->jdx[2*(i-data->p)], A->jdx[2*(i-data->p)+1]);}
            for (int i = 0; i < data->p; i++)
            { tmp->beta[i] = TMP.mean[TMP.root(TMP.id[i])]; }
            for (int i = 0; i < data->p; i++) { if (tmp->graph[i]) tmp->beta[i] = .0; }
        }
        
        // print_dmatrix(tmp->beta, 1, data->p);
        /* save current solution */
        dmat_vcopy(data->p,tmp->beta,out->beta_path + k*data->p);
    }
    
    free_tmpvars_linreg(tmp);
    free(data);
    free(para);
}





