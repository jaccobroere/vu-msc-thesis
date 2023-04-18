//
//  filtering.cpp
//  gen_lasso
//
//  Created by Yunzhang Zhu on 10/24/14.
//  Copyright (c) 2014 OSU. All rights reserved.
//

#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#include "filtering.h"
#include "dmatrix.h"
#include "def.h"
#include "util.h"

void factorize_filtering(int p, int m, double* rho, double* tmp1, double* tmp2,
                         dmatrix_csr* A, dmatrix_csc* A_csc, double* L){
    int i,k;
    dmat_vset(p*p, .0, L);
    for (i = 0; i < p; i++) { L[i+i*p] = 1.0; }
    for (i = 0; i < p; i++) {
        dmat_vset(m, .0, tmp1); // store i-th column of A
        for (k = A_csc->cdx[i]; k < A_csc->cdx[i+1]; k++) {
            tmp1[A_csc->idx[k]] = A_csc->val[k];
        }
        dmat_ySTx(A, tmp1, tmp2);
        dmat_waxpby(p, 1.0, L+i*p, rho[0], tmp2, L+i*p);
    }
    /* L = chol(I + rho * At*A) */
    dmat_potrf(L,p);
}

bool stopping_simple(double* eps, int p, int m,
                     double* beta, double* lambda,
                     double* Y, double* Abeta,
                     double* alpha2A, double* rho,
                     dmatrix_csr* A,double* alpha2,double* fun_val)
{
    double p_star, d_star, tmp;
    dmat_waxpby(p, 1.0, beta, -1.0, Y, alpha2A);
    tmp = dmat_norm2(p, alpha2A);
    p_star = .5*tmp*tmp + lambda[0]*dmat_norm1(m, Abeta);
    
    dmat_ySTx(A, alpha2, alpha2A);
    tmp = dmat_norm2(p, alpha2A);
    d_star = dmat_dot(p, Y, alpha2A) - .5*tmp*tmp;
    
    tmp = p_star - d_star;
    fun_val[0] = p_star;
    // printf("primal obj is: %10.8f \n", p_star);
    // printf("duality gap is: %10.8f \n", tmp);
    return ((abs(tmp)/p) < eps[0]);
    // return (abs(tmp/(max(abs(p_star),abs(d_star)))) < eps[0]);
}

// solves .5 \|Y - beta\|_2^2 + lambda \|A beta\|_1 //
void filtering(int* beta_dim, int* A_row_dim,
               double* Y, int* A_idx,
               int* A_jdx, int* A_rdx, int* nz,
               double* A_val,
               double* theta,
               double* lambda, int* lambda_grid,
               double* rho, int* num_iter,
               double* eps, bool variant,
               bool general_A, bool reporting,
               bool varying_rho, return_data* out){
    int p = beta_dim[0];
    int m = A_row_dim[0];
    int MAX_NUM_ITER = num_iter[0];
    rho[0] = lambda[0];
    bool converge;
    int s, varying_rho_points = 5e2, i, j, k;
    double tmp, lamk, prires, dualres, eps_pri, eps_dual,mu=.1,tau=2.0,fun_val;
    dmatrix_csr* A;
    dmat_new_sparse_csr(&A, m, p, nz[0], A_idx, A_jdx, A_rdx, A_val);
    double* D = (double*) malloc(sizeof(double)*p);
    dmat_vset(p, .0, D);
    double* beta = (double*) malloc(sizeof(double)*p);
    double* alpha1 = (double*) malloc(sizeof(double)*m);
    double* alpha2 = (double*) malloc(sizeof(double)*m);
    double* alpha_bar = (double*) malloc(sizeof(double)*max(m,p));
    double* Abeta = (double*) malloc(sizeof(double)*m);
    double* Abeta_old = (double*) malloc(sizeof(double)*max(m,p));
    double* alpha_barA = (double*) malloc(sizeof(double)*p);
    double* alpha2A = (double*) malloc(sizeof(double)*p);
    double* rhoDpI_inv = (double*) malloc(sizeof(double)*p);
    double* rhoDpI_inv_Y = (double*) malloc(sizeof(double)*p);
    double* allones = (double*) malloc(sizeof(double)*p);
    dmat_vset(p, 1.0, allones);
    double* L = NULL;
    dmatrix_csc* A_csc = NULL;
    if (!variant) {
        L = (double*) malloc(sizeof(double)*p*p);
        dmatcsr_to_csc(A, &A_csc);
    }
    dmat_vset(m, .0, alpha2);
    dmat_vset(m, .0, Abeta);
    dmat_vset(m, .0, alpha_bar);
    dmat_vcopy(p, Y, beta);
    if (general_A) {
        double max_sig_A;
        dmatcsr_to_csc(A, &A_csc);
        get_diag_general(A, A_csc, &max_sig_A);
        dmat_vset(p, max_sig_A, D);
    } else {
        for (int i = 0; i < (A->nz); i++) { D[A->jdx[i]] += 2.0; }
    }
    
    if (varying_rho) {
        if (variant) {
            dmat_waxpby(p, rho[0], D, 1.0, allones, Abeta_old);
            dmat_yinvx(p, Abeta_old, rhoDpI_inv);
            dmat_elemprod(p, rhoDpI_inv, Y, rhoDpI_inv_Y);
        } else {
            factorize_filtering(p, m, rho, alpha_bar, Abeta_old, A, A_csc, L);
        }
    }
    for (k = 0; k < lambda_grid[0]; k++) {
        if (!varying_rho) // non-adaptive strategy, rho = lambda[k]
        {
            (out->chol_num[k])++;
            rho[0] = lambda[k];
            if (variant) // rhoDpI_inv_Y = (rho*D + I)^{-1} Y //
            {
                dmat_waxpby(p, rho[0], D, 1.0, allones, Abeta_old);
                dmat_yinvx(p, Abeta_old, rhoDpI_inv);
                dmat_elemprod(p, rhoDpI_inv, Y, rhoDpI_inv_Y);
            } else // rhoDpI_inv_Y = (rho*AtA + I)^{-1} Y //
            {
                factorize_filtering(p, m, rho, alpha_bar, Abeta_old, A, A_csc, L);
                dmat_vcopy(p, Y, rhoDpI_inv_Y);
                dmat_potrs(L, p, rhoDpI_inv_Y);
            }
        }
        s = 1;
        lamk = lambda[k];
        
        for (i = 0; i < MAX_NUM_ITER; i++) {
            /* update beta */
            dmat_ySTx(A, alpha_bar, alpha_barA);
            dmat_waxpby(p, 1.0, alpha_barA, 1.0, beta, Abeta_old);
            dmat_waxpby(p, 1.0, beta, 1.0, rhoDpI_inv_Y, beta);
            if (variant) { dmat_elemprod(p, rhoDpI_inv, Abeta_old, Abeta_old); }
            else         { dmat_potrs(L, p, Abeta_old);                        }
            dmat_waxpby(p, -1.0, Abeta_old, 1.0, beta, beta);
            
            
            /* update alpha */
            dmat_vcopy(m, Abeta, Abeta_old);
            dmat_ySx(A, beta, Abeta);
            dmat_vcopy(m,alpha2,alpha1);
            for (j = 0; j < m; j++){
                tmp = rho[0]*Abeta[j]+alpha2[j];
                alpha2[j] = (abs(tmp) <= lamk)? tmp : (lamk*sign(tmp));
            }
            
            /* update the parameters */
            // update_parameter(theta, rho);
            
            // stopping rule //
            converge = stopping_simple(eps,p,m,beta,&lamk,Y,Abeta,alpha2A,rho,A,alpha2,&fun_val);
            out->admm_iter_num[k] = i;
            if (converge){
                out->fun_val[k] = fun_val;
                printf("#iter: %d, obj: %10.8f, rho: %f\n", i,fun_val,rho[0]);
                break;
            }
            if (i == (MAX_NUM_ITER - 1)) { out->fun_val[k] = fun_val; }
            if  (reporting) { out->fun_path[k*MAX_NUM_ITER + i] = fun_val; }
            
            /*
            dmat_waxpby(m, 1.0, alpha2, -1.0, alpha1, alpha_bar);
            prires = dmat_norm2(m, alpha_bar) / rho[0];
            
            dmat_waxpby(m, rho[0], Abeta, -rho[0], Abeta_old, Abeta_old);
            dmat_ySTx(A, Abeta_old, alpha_bar);
            // printf("first part: %10.8f\n",rho[0]*dmat_norm2(p, alpha_bar));
            print_dmatrix(alpha_bar, 1, p);
            
            dmat_waxpby(p, 1.0, alpha2A, -1.0, alpha_barA, Abeta_old);
            // printf("second part: %10.8f\n",dmat_norm2(p, Abeta_old));
            print_dmatrix(Abeta_old, 1, p);
            // print_dmatrix(beta, 1, p); print_dmatrix(alpha1, 1, m); print_dmatrix(alpha2, 1, m);
            
            dmat_waxpby(p, 1.0, alpha_barA, 1.0, alpha_bar, alpha_bar);
            dmat_waxpby(p, -1.0, alpha2A, 1.0, alpha_bar, alpha_bar);
            dualres = dmat_norm2(p, alpha_bar);
            
            eps_pri = .1*sqrt(m)*eps[0] + eps[0]*dmat_norm2(m, Abeta);
            eps_dual = .1*sqrt(p)*eps[0] + eps[0]*dmat_norm2(p, alpha2A);
            printf("%5d %10.8f %10.8f %10.8f %10.8f %10.5f \n", i, prires, eps_pri, dualres, eps_dual,rho[0]);
            */
            
            // varying rho strategy //
            if (varying_rho && (i == ((s+1)*s*varying_rho_points/2)))
            {
                s++;
                dmat_waxpby(m, 1.0, alpha2, -1.0, alpha1, alpha_bar);
                prires = dmat_norm2(m, alpha_bar) / rho[0];
                
                dmat_waxpby(m, 1.0, Abeta, -1.0, Abeta_old, Abeta_old);
                dmat_ySTx(A, Abeta_old, alpha_bar);
                dmat_waxpby(p, 1.0, alpha_barA, rho[0], alpha_bar, alpha_bar);
                dmat_waxpby(p, -1.0, alpha2A, 1.0, alpha_bar, alpha_bar);
                dualres = dmat_norm2(p, alpha_bar);
                
                eps_pri = .1*sqrt(m)*eps[0] + eps[0]*dmat_norm2(m, Abeta);
                eps_dual = .1*sqrt(p)*eps[0] + eps[0]*dmat_norm2(p, alpha2A);
                printf("%5d %10.8f %10.8f %10.8f %10.8f %10.5f \n", i, prires, eps_pri, dualres, eps_dual,rho[0]);
                
                if ((prires/eps_pri) < (mu*dualres/eps_dual)) {
                    (out->chol_num[k])++;
                    rho[0] /= tau;
                    if (variant) // rhoDpI_inv_Y = (rho*D + I)^{-1} Y //
                    {
                        dmat_waxpby(p, rho[0], D, 1.0, allones, Abeta_old);
                        dmat_yinvx(p, Abeta_old, rhoDpI_inv);
                        dmat_elemprod(p, rhoDpI_inv, Y, rhoDpI_inv_Y);
                    } else // rhoDpI_inv_Y = (rho*AtA + I)^{-1} Y //
                    {
                        factorize_filtering(p, m, rho, alpha_bar, Abeta_old, A, A_csc, L);
                        dmat_vcopy(p, Y, rhoDpI_inv_Y);
                        dmat_potrs(L, p, rhoDpI_inv_Y);
                    }
                }
                if ((dualres/eps_dual) < (mu*prires/eps_pri)) {
                    (out->chol_num[k])++;
                    rho[0] *= tau;
                    if (variant) // rhoDpI_inv_Y = (rho*D + I)^{-1} Y //
                    {
                        dmat_waxpby(p, rho[0], D, 1.0, allones, Abeta_old);
                        dmat_yinvx(p, Abeta_old, rhoDpI_inv);
                        dmat_elemprod(p, rhoDpI_inv, Y, rhoDpI_inv_Y);
                    } else // rhoDpI_inv_Y = (rho*AtA + I)^{-1} Y //
                    {
                        factorize_filtering(p, m, rho, alpha_bar, Abeta_old, A, A_csc, L);
                        dmat_vcopy(p, Y, rhoDpI_inv_Y);
                        dmat_potrs(L, p, rhoDpI_inv_Y);
                    }
                }
            }
            /* alpha_bar = alpha2 + theta*(alpha2 - alpha1) */
            dmat_waxpby(m, 1.0+theta[0], alpha2, -theta[0], alpha1, alpha_bar);
        }
        
        dmat_vcopy(p, beta, out->beta_path + k*p);
        
    }
    free(D); free(alpha1); free(alpha2); free(alpha_bar);
    free(rhoDpI_inv); free(rhoDpI_inv_Y); free(beta);
    free(Abeta); free(alpha_barA); free(alpha2A); free(A);
    if (!variant) { dmat_csc_free(A_csc); free(L); }
}