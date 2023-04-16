#include <math.h>
#include <vector>
#include <algorithm>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>


#include "def.h"
#include "dmatrix.h"
#include "opt.h"
#include "util.h"

// return tau, an estimate of the maximum singular value of A //
void get_diag_general(dmatrix_csr* A, dmatrix_csc* A_csc, double* tau)
{
    int p               = A->n;
    int m               = A->m;
    int* rdx = A->rdx;
    int* cdx = A_csc->cdx;
    double* val = A->val;
    double* val_csc = A_csc->val;
    double max_Acol = .0;
    double max_Atrow = .0;
    double* row_norm = (double*) malloc(sizeof(double)*m);
    double* col_norm  = (double*) malloc(sizeof(double)*p);
    double* tmp = (double*) malloc(sizeof(double)*max(m, p));
    for (int i = 0; i < m; i++) { row_norm[i] = dmat_norm1(rdx[i+1]-rdx[i], val+rdx[i]); }
    for (int j = 0; j < p; j++) { col_norm[j] = dmat_norm1(cdx[j+1]-cdx[j], val_csc+cdx[j]); }
    dmat_yabsSx(A, col_norm, tmp);
    for (int i = 0; i < m; i++) { max_Acol = max(max_Acol, tmp[i]); }
    dmat_yabsSTx(A, row_norm, tmp);
    for (int j = 0; j < p; j++) { max_Atrow = max(max_Atrow, tmp[j]); }
    tau[0] = min(max_Acol,max_Atrow);
    printf("max singular value is: %f\n", tau[0]);
}

void ini_linreg_prob(int* p, int* m,
                     int* n, double* Y,
                     double* X, dmatrix_csr* A,
                     double* theta,
                     double* lambda, int* lambda_grid,
                     double* rho, int* num_iter,
                     bool variant, bool varying_rho,
                     bool general_A, bool reporting,
                     double* eps_abs, double* eps_rel,
                     linreg_data** data_, parameters** para_, tmpvars** tmp_)
{
    linreg_data* data   = (linreg_data*) malloc(sizeof(linreg_data));
    parameters* para    = (parameters*) malloc(sizeof(parameters));
    tmpvars* tmp        = (tmpvars*) malloc(sizeof(tmpvars));
    data->p             = p[0];
    data->m             = m[0];
    data->n             = n[0];
    data->X             = X;
    data->Y             = Y;
    data->A             = A;
    data->lambda        = lambda;
    data->lambda_grid   = lambda_grid[0];
    // data->constraints   = NULL;
    
    para->theta         = theta;
    para->rho           = rho;
    para->max_num_iter  = num_iter;
    para->eps_abs       = eps_abs;
    para->eps_rel       = eps_rel;
    para->variant       = variant;
    para->varying_rho   = varying_rho;
    para->general_A     = general_A;
    para->reporting     = reporting;
    para->mu            = .1;
    para->tau           = 2.0;
    
    tmp->npm            = max(data->n,max(data->p,data->m));
    tmp->skinny         = false;
    tmp->k              = 0;
    tmp->prires         = .0;
    tmp->dualres        = .0;
    tmp->eps_pri        = .0;
    tmp->eps_dual       = .0;
    tmp->beta           = (double*) malloc(sizeof(double)*data->p);
    tmp->alpha          = (double*) malloc(sizeof(double)*data->m);
    tmp->alpha1         = (double*) malloc(sizeof(double)*data->m);
    tmp->alpha2         = (double*) malloc(sizeof(double)*data->m);
    tmp->alpha_bar      = (double*) malloc(sizeof(double)*data->m);
    tmp->graph          = (bool*  ) malloc(sizeof(bool  )*data->m);
    tmp->temp           = (double*) malloc(sizeof(double)*tmp->npm);
    tmp->beta_old       = (double*) malloc(sizeof(double)*data->p);
    tmp->prires_vec     = (double*) malloc(sizeof(double)*data->m);
    tmp->dualres_vec    = (double*) malloc(sizeof(double)*data->p);
    tmp->Abeta          = (double*) malloc(sizeof(double)*data->m);
    tmp->alpha_barA     = (double*) malloc(sizeof(double)*data->p);
    tmp->alpha2A        = (double*) malloc(sizeof(double)*data->p);
    tmp->XTY            = (double*) malloc(sizeof(double)*data->p);
    // tmp->AtAbeta        = (double*) malloc(sizeof(double)*data->p);
    // tmp->AtAbeta_old    = (double*) malloc(sizeof(double)*data->p);
    tmp->D              = NULL;
    tmp->diag_mat       = NULL;
    tmp->diag_mat_inv   = NULL;
    tmp->L              = NULL;
    tmp->A_csc          = NULL;
    
    dmat_yATx(data->n, data->p, data->X, data->Y, tmp->XTY);
    dmat_vset(data->p, .0, tmp->beta);
    dmat_vset(data->m, .0, tmp->alpha);
    dmat_vset(data->m, .0, tmp->alpha1);
    dmat_vset(data->m, .0, tmp->alpha2);
    dmat_vset(data->m, .0, tmp->alpha_bar);
    
    /*
     if (reporting) {
        // data->constraints = (double*) malloc(sizeof(double)*lambda_grid[0]*num_iter[0]);
        data->fun_vals    = (double*) malloc(sizeof(double)*lambda_grid[0]*num_iter[0]);
     }
    */
    
    if (variant) {
        tmp->skinny = (data->n >= data->p);
        tmp->D = (double*) malloc(sizeof(double)*data->p);
        dmat_vset(data->p, .0, tmp->D);
        tmp->diag_mat = (double*) malloc(sizeof(double)*data->p);
        tmp->diag_mat_inv = (double*) malloc(sizeof(double)*data->p);
        if (general_A) {
            double max_sig_A;
            dmatcsr_to_csc(data->A, &(tmp->A_csc));
            get_diag_general(data->A, tmp->A_csc, &max_sig_A);
            dmat_vset(data->p, max_sig_A, tmp->D);
        } else {
            for (int i = 0; i < (data->A->nz); i++) { tmp->D[data->A->jdx[i]] += 2.0; }
        }
        if (tmp->skinny) { tmp->L = (double*) malloc(sizeof(double)*data->p*data->p); }
        else { tmp->L = (double*) malloc(sizeof(double)*data->n*data->n); }
    } else {
        tmp->L = (double*) malloc(sizeof(double)*data->p*data->p);
        dmatcsr_to_csc(data->A, &(tmp->A_csc));
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
        if (tmp->Abeta)         { free(tmp->Abeta);         }
        if (tmp->alpha_barA)    { free(tmp->alpha_barA);    }
        if (tmp->alpha2A)       { free(tmp->alpha2A);       }
        if (tmp->diag_mat)      { free(tmp->diag_mat);      }
        if (tmp->D)             { free(tmp->D);             }
        if (tmp->diag_mat_inv)  { free(tmp->diag_mat_inv);  }
        if (tmp->L)             { free(tmp->L);             }
        if (tmp->XTY)           { free(tmp->XTY);           }
        if (tmp->temp)          { free(tmp->temp);          }
        if (tmp->graph)         { free(tmp->graph);         }
        if (tmp->A_csc)         { dmat_csc_free(tmp->A_csc);}
        //  if (tmp->AtAbeta)       { free(tmp->AtAbeta);       }
        //  if (tmp->AtAbeta_old)   { free(tmp->AtAbeta_old);   }
        free(tmp);
    }
    
}


bool stopping_simple(double* eps, int p, int m,
                     double* beta, double* lambda,
                     double* Y, double* Abeta,
                     double* alpha2A, double* rho,
                     dmatrix_csr* A,double* alpha2){
    double p_star, d_star, tmp;
    dmat_waxpby(p, 1.0, beta, -1.0, Y, alpha2A);
    tmp = dmat_norm2(p, alpha2A);
    p_star = .5*tmp*tmp + lambda[0]*dmat_norm1(m, Abeta);
    
    dmat_ySTx(A, alpha2, alpha2A);
    tmp = dmat_norm2(p, alpha2A);
    d_star = dmat_dot(p, Y, alpha2A) - .5*tmp*tmp;
    
    tmp = p_star - d_star;
    // printf("duality gap is: %10.8f \n", tmp);
    return (abs(tmp/(max(abs(p_star),abs(d_star)))) < eps[0]);
}

void get_diag(dmatrix_csr* A, double* diag_mat){
    for (int i = 0; i < (A->nz); i++) {
        diag_mat[A->jdx[i]] += 2;
    }
}

void get_diag_linreg(dmatrix_csr* A, double* diag_mat,double* rho){
    for (int i = 0; i < (A->nz); i++) {
        diag_mat[A->jdx[i]] += 2;
    }
    dmat_waxpby(A->n, rho[0], diag_mat, .0, diag_mat, diag_mat);
}


void get_diag_identity(dmatrix_csr* A, double* diag_mat){
    get_diag(A, diag_mat);
    double max = .0;
    double tmp = .0;
    for (int i = 0; i < A->n; i++) {
        tmp = diag_mat[i];
        max = (max > tmp)? max:tmp;
    }
    for (int i = 0; i < A->n; i++) {
        diag_mat[i] = max;
    }
}

void alpha_update_l1(int m, double* Abeta,
                     double* alpha1, double* alpha2,
                     double* lambda, double* rho){
    dmat_vcopy(m,alpha2,alpha1);
    double tmp1 = lambda[0] / rho[0];
    for (int i = 0; i < m; i++){
        double tmp = Abeta[i]+alpha2[i];
        alpha2[i] = (abs(tmp) <= tmp1)? tmp : (tmp1*sign(tmp));
    }
}

void alpha_update_linreg(int m, dmatrix_csr* A,
                         double* beta, double* Abeta,
                         double* alpha1, double* alpha2,
                         double* lambda, double* rho){
    dmat_ySx(A, beta, Abeta);
    dmat_vcopy(m,alpha2,alpha1);
    double tmp1 = lambda[0] / rho[0];
    for (int i = 0; i < m; i++){
        double tmp = Abeta[i]+alpha2[i];
        alpha2[i] = (abs(tmp) <= tmp1)? tmp : (tmp1*sign(tmp));
    }
}

void graph_update(int m, double* Abeta,
                  double* alpha2, double* lambda,
                  double* rho,bool* graph){
    double tmp1 = lambda[0] / rho[0];
    double tmp;
    for (int i = 0; i < m; i++){
        tmp = Abeta[i]+alpha2[i];
        graph[i] = (abs(tmp) <= tmp1);
    }
}

void beta_update_simple(int p, double* Y,
                        double* alpha_barA, double* beta,
                        double* diag_mat, double* rho){
    double rho_inv = 1.0 / rho[0];
    double tmp;
    for (int i = 0; i < p; i++) {
        tmp = 1.0 / (diag_mat[i] + rho_inv);
        beta[i] = tmp*(diag_mat[i]*beta[i] + rho_inv*Y[i] - alpha_barA[i]);
    }
}

/* solves .5 \|Y - beta\|_2^2 + lambda \sum_{(i,i')\in E} |beta_i-beta_i'|,
 where A is a m by p matrix */
void filtering_old(int* beta_dim, int* A_row_dim,
                   double* Y, int* A_idx,
                   int* A_jdx, int* A_rdx, int* nz,
                   double* A_val, double* beta,
                   double* alpha, double* theta,
                   double* lambda, int* lambda_grid,
                   double* rho, int* num_iter,
                   double* eps){
    int p = beta_dim[0];
    int m = A_row_dim[0];
    int MAX_NUM_ITER = num_iter[0];
    dmatrix_csr* A;
    dmat_new_sparse_csr(&A, m, p, nz[0], A_idx, A_jdx, A_rdx, A_val);
    double* diag_mat = (double*) malloc(sizeof(double)*p);
    dmat_vset(p, .0, diag_mat);
    double* alpha1 = (double*) malloc(sizeof(double)*m);
    double* alpha2 = (double*) malloc(sizeof(double)*m);
    double* alpha_bar = (double*) malloc(sizeof(double)*m);
    double* Abeta = (double*) malloc(sizeof(double)*m);
    double* alpha_barA = (double*) malloc(sizeof(double)*p);
    double* alpha2A = (double*) malloc(sizeof(double)*p);
    dmat_vcopy(m,alpha,alpha1); dmat_vcopy(m,alpha,alpha2);
    dmat_waxpby(m, 1.0+theta[0], alpha2, -theta[0],alpha1, alpha_bar);
    get_diag(A, diag_mat);
    int total_num_iters = 0;
    for (int k = 0; k < lambda_grid[0]; k++) {
        // optional update alpha first //
        // dmat_ySx(A, beta, Abeta);
        // alpha_update_l1(m, Abeta, alpha1, alpha2, lambda+k, rho);
        for (int i = 0; i < MAX_NUM_ITER; i++) {
            /* update beta */
            dmat_ySTx(A, alpha_bar, alpha_barA);
            beta_update_simple(p, Y, alpha_barA, beta, diag_mat, rho);
            /* update alpha*/
            dmat_ySx(A, beta, Abeta);
            alpha_update_l1(m, Abeta, alpha1, alpha2, lambda+k, rho);
            
            /* update the parameters */
            // update_parameter(theta, rho);
            
            /* alpha_bar = alpha2 + theta*(alpha2 - alpha1) */
            dmat_waxpby(m, 1.0+theta[0], alpha2, -theta[0], alpha1, alpha_bar);
            // stopping rule //
            dmat_ySTx(A, alpha2, alpha2A);
            if (stopping_simple(eps,p,m,beta,lambda+k,Y,Abeta,alpha2A,rho,A,alpha2)){
                total_num_iters += i;
                printf("converge at %d iterations!\n", i);
                break;
            }
        }
    }
    printf("total number of iterations is %d\n", total_num_iters);
    free(diag_mat); free(alpha1); free(alpha2); free(alpha_bar);
    free(Abeta); free(alpha_barA); free(alpha2A); free(A);
}


// solves .5 \|Y - beta\|_2^2 + lambda \|A beta\|_1 //
void filtering(int* beta_dim, int* A_row_dim,
               double* Y, int* A_idx,
               int* A_jdx, int* A_rdx, int* nz,
               double* A_val, double* beta,
               double* theta,
               double* lambda, int* lambda_grid,
               double* rho, int* num_iter,
               double* eps, double* eps_abs,
               double* eps_rel, bool general_A,
               bool varying_rho){
    int p = beta_dim[0];
    int m = A_row_dim[0];
    int MAX_NUM_ITER = num_iter[0];
    double tmp, lamk, prires, dualres, eps_pri, eps_dual,mu=.1,tau=2.0;
    dmatrix_csr* A;
    dmatrix_csc* A_csc = NULL;
    dmat_new_sparse_csr(&A, m, p, nz[0], A_idx, A_jdx, A_rdx, A_val);
    double* D = (double*) malloc(sizeof(double)*p);
    dmat_vset(p, .0, D);
    double* alpha1 = (double*) malloc(sizeof(double)*m);
    double* alpha2 = (double*) malloc(sizeof(double)*m);
    double* alpha_bar = (double*) malloc(sizeof(double)*max(m,p));
    double* Abeta = (double*) malloc(sizeof(double)*m);
    double* Abeta_old = (double*) malloc(sizeof(double)*max(m,p));
    double* alpha_barA = (double*) malloc(sizeof(double)*p);
    double* alpha2A = (double*) malloc(sizeof(double)*p);
    double* rhoD = (double*) malloc(sizeof(double)*p);
    double* rhoDpI_inv = (double*) malloc(sizeof(double)*p);
    double* allones = (double*) malloc(sizeof(double)*p);
    dmat_vset(m, .0, alpha2);
    dmat_vset(m, .0, Abeta);
    dmat_vset(m, .0, alpha_bar);
    dmat_vset(p, 1.0, allones);
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
        dmat_waxpby(p, rho[0], D, .0, rhoD, rhoD);
        dmat_waxpby(p, 1.0, rhoD, 1.0, allones, Abeta_old);
        dmat_yinvx(p, Abeta_old, rhoDpI_inv);
    }
    for (int k = 0; k < lambda_grid[0]; k++) {
        if (!varying_rho) {
            rho[0] = lambda[k];
            dmat_waxpby(p, rho[0], D, .0, rhoD, rhoD);
            dmat_waxpby(p, 1.0, rhoD, 1.0, allones, Abeta_old);
            dmat_yinvx(p, Abeta_old, rhoDpI_inv);
        }
        lamk = lambda[k];
        for (int i = 0; i < MAX_NUM_ITER; i++) {
            /* update beta */
            dmat_elemprod(p, rhoD, beta, beta);
            dmat_waxpby(p, 1.0, Y, 1.0, beta, beta);
            dmat_ySTx(A, alpha_bar, alpha_barA);
            dmat_waxpby(p, -1.0, alpha_barA, 1.0, beta, beta);
            dmat_elemprod(p, rhoDpI_inv, beta, beta);
            
            /* update alpha*/
            dmat_vcopy(m, Abeta, Abeta_old);
            dmat_ySx(A, beta, Abeta);
            dmat_vcopy(m,alpha2,alpha1);
            for (int j = 0; j < m; j++){
                tmp = rho[0]*Abeta[j]+alpha2[j];
                alpha2[j] = (abs(tmp) <= lamk)? tmp : (lamk*sign(tmp));
            }
            
            /* update the parameters */
            // update_parameter(theta, rho);
            
            // stopping rule //
            if (stopping_simple(eps,p,m,beta,&lamk,Y,Abeta,alpha2A,rho,A,alpha2)){
                printf("converge at %d iterations!\n", i);
                break;
            } else {
                if (varying_rho && (i % p==(p-1))) /* not stable!! */
                {
                    dmat_waxpby(m, 1.0, alpha2, -1.0, alpha1, alpha_bar);
                    prires = dmat_norm2(m, alpha_bar);
                    
                    dmat_waxpby(m, 1.0, Abeta, -1.0, Abeta_old, Abeta_old);
                    dmat_ySTx(A, Abeta_old, alpha_bar);
                    dmat_waxpby(p, 1.0, alpha_barA, rho[0], alpha_bar, alpha_bar);
                    dmat_waxpby(p, -1.0, alpha2A, 1.0, alpha_bar, alpha_bar);
                    dualres = dmat_norm2(p, alpha_bar);
                    
                    eps_pri = sqrt(m)*eps_abs[0] + eps_rel[0]*dmat_norm2(m, Abeta);
                    eps_dual = sqrt(p)*eps_abs[0] + eps_rel[0]*dmat_norm2(p, alpha2A);
                    
                    if (prires < (mu*dualres)) {
                        rho[0] /= tau;
                        dmat_waxpby(p, rho[0], D, .0, rhoD, rhoD);
                        dmat_waxpby(p, 1.0, rhoD, 1.0, allones, Abeta_old);
                        dmat_yinvx(p, Abeta_old, rhoDpI_inv);
                    }
                    if (dualres < (mu*prires)) {
                        rho[0] *= tau;
                        dmat_waxpby(p, rho[0], D, .0, rhoD, rhoD);
                        dmat_waxpby(p, 1.0, rhoD, 1.0, allones, Abeta_old);
                        dmat_yinvx(p, Abeta_old, rhoDpI_inv);
                    }
                    // printf("%3d %10.5f %10.5f %10.5f %10.5f %10.5f \n", i, prires, eps_pri, dualres, eps_dual,rho[0]);
                }
            }
            /* alpha_bar = alpha2 + theta*(alpha2 - alpha1) */
            dmat_waxpby(m, 1.0+theta[0], alpha2, -theta[0], alpha1, alpha_bar);
        }
        
    }
    free(D); free(alpha1); free(alpha2); free(alpha_bar);
    free(Abeta); free(alpha_barA); free(alpha2A); free(A);
}

void allocate_vars(int n,int p,int m, dmatrix_csr* A, double** D_,
                   double** diag_mat_,double** diag_mat_inv_,
                   double** L_,dmatrix_csc** A_csc_,bool* skinny,
                   bool variant){
    double *diag_mat=NULL, *diag_mat_inv=NULL, *D=NULL, *L=NULL;
    dmatrix_csc* A_csc=NULL;
    if (variant) {
        skinny[0] = (n >= p);
        D = (double*) malloc(sizeof(double)*p);
        dmat_vset(p, .0, D);
        diag_mat = (double*) malloc(sizeof(double)*p);
        diag_mat_inv = (double*) malloc(sizeof(double)*p);
        for (int i = 0; i < (A->nz); i++) { D[A->jdx[i]] += 2.0; }
        if (skinny[0]) { L = (double*) malloc(sizeof(double)*p*p); }
        else { L = (double*) malloc(sizeof(double)*n*n); }
    } else {
        L = (double*) malloc(sizeof(double)*p*p);
        dmatcsr_to_csc(A, &A_csc);
    }
    *D_ = D;
    *diag_mat_ = diag_mat;
    *diag_mat_inv_ = diag_mat_inv;
    *L_ = L;
    *A_csc_ = A_csc;
}

void factorize(int n, int p, int m,
               bool variant, bool skinny, double* D,
               double* diag_mat, double* diag_mat_inv,
               dmatrix_csr* A, dmatrix_csc* A_csc,
               double* X, double* alpha, double* beta_old,
               double* rho, double* L){
    if (variant) {
        // dmat_vset(p, .0, diag_mat);
        // get_diag_linreg(A, diag_mat,rho); /* diag_mat = rho*D */
        dmat_waxpby(p, rho[0], D, .0, D, diag_mat); /* diag_mat = rho*D */
        dmat_yinvx(p, diag_mat, diag_mat_inv);
        if (skinny) {
            dmat_B_ATA(n, p, X, L);
            for (int i = 0; i < p; i++) {
                L[i*p+i] += diag_mat[i];
            }
            /* L = chol(rho * D + Xt*X) */
            dmat_potrf(L,p);
        } else {
            dmat_B_ADAT(n, p, X, diag_mat_inv, L);
            for (int i = 0; i < n; i++) {
                L[i*n+i] += 1.0;
            }
            /* L = chol(I + rho^{-1} X D^{-1} Xt) */
            dmat_potrf(L,n);
        }
    } else /* standard ADMM */{
        dmat_B_ATA(n, p, X, L);
        for (int i = 0; i < p; i++) {
            dmat_vset(m, .0, alpha); // store i-th column of A
            for (int k = A_csc->cdx[i]; k < A_csc->cdx[i+1]; k++) {
                alpha[A_csc->idx[k]] = A_csc->val[k];
            }
            dmat_ySTx(A, alpha, beta_old);
            dmat_waxpby(p, 1.0, L+i*p, rho[0], beta_old, L+i*p);
        }
        /* L = chol(rho * At*A + Xt*X) */
        dmat_potrf(L,p);
    }
}

void beta_update_linreg(int p, int n, double* XTY, double* X,
                        double* alpha_bar, double* alpha,
                        double* alpha_barA, double* beta,
                        double* diag_mat, double* diag_mat_inv,
                        double* rho, dmatrix_csr* A,
                        double* L, double* tmp,
                        double* dualres_vec,
                        bool variant,bool skinny){
    dmat_ySTx(A, alpha_bar, alpha_barA);
    dmat_waxpby(p, 1.0, XTY, -rho[0], alpha_barA, tmp);
    if (variant) {
        dmat_elemprod(p, diag_mat, beta, beta);
        dmat_waxpby(p, 1.0, beta, 1.0, tmp, beta);
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
    } else {
        dmat_ySx(A, beta, alpha);
        dmat_ySTx(A, alpha, beta);
        dmat_waxpby(p, rho[0], beta, 1.0, tmp, beta);
        dmat_potrs(L, p, beta);
    }
}



/* calculate function value */
double fun_val(linreg_data* data, parameters* para, tmpvars* tmp)
{
    dmat_yAx(data->n, data->p, data->X, tmp->beta, tmp->temp);
    dmat_waxpby(data->n, 1.0, data->Y, -1.0, tmp->temp, tmp->temp);
    double fun_val = dmat_norm2(data->n, tmp->temp);
    return (.5*fun_val*fun_val + data->lambda[0] * dmat_norm1(data->m, tmp->Abeta));
}

void stopping_para_update(int m, int n, int p, double* beta,
                          double* alpha1, double* alpha2,
                          double* prires_vec, double* rho,
                          double* dualres_vec, double* beta_old,
                          double* diag_mat, dmatrix_csr* A,
                          double* alpha_barA, double* alpha2A,
                          double* alpha, double* eps_abs,
                          double* eps_rel, double* Abeta,
                          bool variant, tmpvars* tmp)
{
    // calculate primal and dual residuals
    dmat_waxpby(m, 1.0, alpha2, -1.0, alpha1, prires_vec);
    tmp->prires = dmat_norm2(m, prires_vec);
    
    dmat_waxpby(p, 1.0, beta, -1.0, beta_old, dualres_vec);
    
    
    /*
    // use separate stopping criterions for special and standard admm //
    if (variant) {
        dmat_elemprod(p, diag_mat, dualres_vec, dualres_vec);
        dmat_waxpby(p, .0, dualres_vec, 1.0/rho[0], dualres_vec, dualres_vec);
    } else {
        dmat_ySx(A, dualres_vec, alpha);
        dmat_ySTx(A, alpha, dualres_vec);
    }
    dmat_waxpby(p, 1.0, alpha_barA, 1.0, dualres_vec, dualres_vec);
    dmat_ySTx(A, alpha2, alpha2A);
    dmat_waxpby(p, -1.0, alpha2A, 1.0, dualres_vec, dualres_vec);
    tmp->dualres = rho[0]*dmat_norm2(p, dualres_vec);
    
    // calculate eps_pri and eps_dual
    if (variant) {
        dmat_elemprod(p, diag_mat, beta, dualres_vec);
        tmp->eps_pri = sqrt(m+p)*eps_abs[0] + eps_rel[0] * sqrt(dmat_dot(p, dualres_vec, beta) / rho[0]);
    } else {
        tmp->eps_pri = sqrt(m)*eps_abs[0]   + eps_rel[0] * dmat_norm2(m, Abeta);
    }
    tmp->eps_dual = sqrt(p)*eps_abs[0] + rho[0]*dmat_norm2(p, alpha2A)*eps_rel[0];
    */
    
    
    // use the standard admm's stopping criterion //
    dmat_ySx(A, dualres_vec, alpha);
    dmat_ySTx(A, alpha, dualres_vec);
    dmat_waxpby(p, 1.0, alpha_barA, 1.0, dualres_vec, dualres_vec);
    dmat_ySTx(A, alpha2, alpha2A);
    dmat_waxpby(p, -1.0, alpha2A, 1.0, dualres_vec, dualres_vec);
    tmp->dualres = rho[0]*dmat_norm2(p, dualres_vec);
    
    tmp->eps_pri = sqrt(m)*eps_abs[0]   + eps_rel[0] * dmat_norm2(m, Abeta);
    tmp->eps_dual = sqrt(p)*eps_abs[0] + rho[0]*dmat_norm2(p, alpha2A)*eps_rel[0];
}

/*
 void stopping_para_update_new(int m, int n, int p, double* beta,
 double* alpha1, double* alpha2,
 double* prires_vec, double* rho,
 double* dualres_vec, double* beta_old,
 double* diag_mat, dmatrix_csr* A,
 double* alpha_barA, double* alpha2A,
 double* alpha, double* eps_abs,
 double* eps_rel, double* Abeta,
 bool variant, tmpvars* tmp)
 {
 // calculate primal and dual residuals
 dmat_waxpby(m, 1.0, alpha2, -1.0, alpha1, prires_vec);
 tmp->prires = dmat_norm2(m, prires_vec);
 tmp->eps_pri = eps_rel[0]*dmat_norm2(m, alpha2)+sqrt(m)*eps_abs[0];
 
 
 dmat_waxpby(p, 1.0, beta, -1.0, beta_old, dualres_vec);
 dmat_elemprod(p, diag_mat, dualres_vec, dualres_vec);
 dmat_waxpby(p, .0, dualres_vec, 1.0/rho[0], dualres_vec, dualres_vec);
 dmat_waxpby(p, 1.0, alpha_barA, 1.0, dualres_vec, dualres_vec);
 dmat_ySTx(A, alpha2, alpha2A);
 dmat_waxpby(p, -1.0, alpha2A, 1.0, dualres_vec, dualres_vec);
 tmp->dualres = dmat_norm2(p, dualres_vec);
 
 dmat_elemprod(p, diag_mat, beta, dualres_vec);
 dmat_waxpby(p, -rho[0], alpha2A, 1.0, dualres_vec, dualres_vec);
 dmat_ySTx(A, alpha1, alpha2A);
 dmat_waxpby(p, rho[0], alpha2A, 1.0, dualres_vec, dualres_vec);
 tmp->eps_dual = eps_rel[0]*dmat_norm2(p, dualres_vec)/rho[0]+sqrt(p)*eps_abs[0];
 }
 */

/*
 void stopping_para_update_old(linreg_data* data, parameters* para, tmpvars* tmp)
 {
 // calculate primal and dual residuals
 dmat_waxpby(data->m, 1.0, tmp->alpha2, -1.0, tmp->alpha1, tmp->prires_vec);
 tmp->prires = dmat_norm2(data->m, tmp->prires_vec);
 
 dmat_waxpby(data->p, 1.0, tmp->beta, -1.0, tmp->beta_old, tmp->dualres_vec);
 if (para->variant) {
 dmat_elemprod(data->p, tmp->diag_mat, tmp->dualres_vec, tmp->dualres_vec);
 dmat_waxpby(data->p, .0, tmp->dualres_vec,
 1.0/para->rho[0], tmp->dualres_vec, tmp->dualres_vec);
 } else {
 dmat_ySx(data->A, tmp->dualres_vec, tmp->alpha);
 dmat_ySTx(data->A, tmp->alpha, tmp->dualres_vec);
 }
 dmat_waxpby(data->p, 1.0, tmp->alpha_barA, 1.0, tmp->dualres_vec, tmp->dualres_vec);
 dmat_ySTx(data->A, tmp->alpha2, tmp->alpha2A);
 dmat_waxpby(data->p, -1.0, tmp->alpha2A, 1.0, tmp->dualres_vec, tmp->dualres_vec);
 tmp->dualres = para->rho[0]*dmat_norm2(data->p, tmp->dualres_vec);
 
 // calculate eps_pri and eps_dual
 if (para->variant) {
 dmat_elemprod(data->p, tmp->diag_mat, tmp->beta, tmp->dualres_vec);
 tmp->eps_pri = sqrt(data->m+data->p)*para->eps_abs[0] +
 para->eps_rel[0] * sqrt(dmat_dot(data->p, tmp->dualres_vec, tmp->beta) / para->rho[0]);
 } else {
 tmp->eps_pri = para->eps_rel[0] * dmat_norm2(data->m, tmp->Abeta)
 + sqrt(data->m)*para->eps_abs[0];
 }
 tmp->eps_dual = sqrt(data->p)*para->eps_abs[0] + para->rho[0]*dmat_norm2(data->p, tmp->alpha2A)*para->eps_rel[0];
 }
 */

void linreg_sub(linreg_data* data, return_data* out, parameters* para, tmpvars* tmp)
{
    // printf("%3s %10s %10s %10s %10s %10s\n", "#", "r norm", "eps_pri", "s norm", "eps_dual", "objective");
    // tmp->num_of_chol = 0;
    out->chol_num[tmp->k] = 0;
    if (!para->varying_rho) // if not varying rho, set rho = lambda, //
        // then factor once at each iteration   //
    {
        // tmp->num_of_chol++;
        (out->chol_num[tmp->k])++;
        para->rho[0] = max(1e-10,data->lambda[0]);
        factorize(data->n, data->p, data->m, para->variant,
                  tmp->skinny, tmp->D, tmp->diag_mat, tmp->diag_mat_inv,
                  data->A, tmp->A_csc, data->X, tmp->alpha, tmp->beta_old,
                  para->rho, tmp->L);
    }
    for (int i = 0; i < para->max_num_iter[0]; i++) {
        dmat_vcopy(data->p, tmp->beta, tmp->beta_old); // save beta
        
        /* update beta */
        beta_update_linreg(data->p, data->n,
                           tmp->XTY, data->X,
                           tmp->alpha_bar,tmp->alpha,
                           tmp->alpha_barA, tmp->beta,
                           tmp->diag_mat, tmp->diag_mat_inv,
                           para->rho, data->A, tmp->L,
                           tmp->temp, tmp->dualres_vec,
                           para->variant, tmp->skinny);
        
        
        /* update alpha */
        alpha_update_linreg(data->m,data->A,tmp->beta,
                            tmp->Abeta, tmp->alpha1, tmp->alpha2,
                            data->lambda, para->rho);
        
        /* alpha_bar = alpha2 + theta*(alpha2 - alpha1) */
        dmat_waxpby(data->m, 1.0+para->theta[0], tmp->alpha2,
                    -para->theta[0], tmp->alpha1, tmp->alpha_bar);
        
        // stopping_para_update(data,para,tmp);
        stopping_para_update(data->m, data->n, data->p, tmp->beta,
                             tmp->alpha1, tmp->alpha2,
                             tmp->prires_vec, para->rho,
                             tmp->dualres_vec, tmp->beta_old,
                             tmp->diag_mat, data->A,
                             tmp->alpha_barA, tmp->alpha2A,
                             tmp->alpha, para->eps_abs,
                             para->eps_rel, tmp->Abeta,
                             para->variant, tmp);
        
        // printf("%3d %10.4f %10.4f %10.4f %10.4f %10.4f\n", i, tmp->prires, tmp->eps_pri, tmp->dualres, tmp->eps_dual, fun_val(data,para,tmp));
        
        // varying rho strategy //
        if (para->varying_rho && (i % data->n == data->n-1))
        {
            /*
            if (tmp->prires < para->mu*tmp->dualres) {
                (out->chol_num[tmp->k])++;
                para->rho[0] /= para->tau;
                factorize(data->n, data->p, data->m, para->variant,
                          tmp->skinny, tmp->D, tmp->diag_mat, tmp->diag_mat_inv,
                          data->A, tmp->A_csc, data->X, tmp->alpha, tmp->beta_old,
                          para->rho, tmp->L);
                dmat_waxpby(data->m, para->tau, tmp->alpha1, .0, tmp->alpha1, tmp->alpha1);
                dmat_waxpby(data->m, para->tau, tmp->alpha2, .0, tmp->alpha2, tmp->alpha2);
                dmat_waxpby(data->m, para->tau, tmp->alpha_bar, .0, tmp->alpha_bar, tmp->alpha_bar);
            }
            if (tmp->dualres < para->mu*tmp->prires) {
                (out->chol_num[tmp->k])++;
                para->rho[0] *= para->tau;
                factorize(data->n, data->p, data->m, para->variant,
                          tmp->skinny, tmp->D, tmp->diag_mat, tmp->diag_mat_inv,
                          data->A, tmp->A_csc, data->X, tmp->alpha, tmp->beta_old,
                          para->rho, tmp->L);
                dmat_waxpby(data->m, 1.0/para->tau, tmp->alpha1, .0, tmp->alpha1, tmp->alpha1);
                dmat_waxpby(data->m, 1.0/para->tau, tmp->alpha2, .0, tmp->alpha2, tmp->alpha2);
                dmat_waxpby(data->m, 1.0/para->tau, tmp->alpha_bar, .0, tmp->alpha_bar, tmp->alpha_bar);
            } */
            
            if ((tmp->prires / tmp->eps_pri) < para->mu*(tmp->dualres / tmp->eps_dual)) {
                (out->chol_num[tmp->k])++;
                para->rho[0] /= para->tau;
                factorize(data->n, data->p, data->m, para->variant,
                          tmp->skinny, tmp->D, tmp->diag_mat, tmp->diag_mat_inv,
                          data->A, tmp->A_csc, data->X, tmp->alpha, tmp->beta_old,
                          para->rho, tmp->L);
                dmat_waxpby(data->m, para->tau, tmp->alpha1, .0, tmp->alpha1, tmp->alpha1);
                dmat_waxpby(data->m, para->tau, tmp->alpha2, .0, tmp->alpha2, tmp->alpha2);
                dmat_waxpby(data->m, para->tau, tmp->alpha_bar, .0, tmp->alpha_bar, tmp->alpha_bar);
            }
            if ((tmp->dualres / tmp->eps_dual) < para->mu*(tmp->prires / tmp->eps_pri)) {
                (out->chol_num[tmp->k])++;
                para->rho[0] *= para->tau;
                factorize(data->n, data->p, data->m, para->variant,
                          tmp->skinny, tmp->D, tmp->diag_mat, tmp->diag_mat_inv,
                          data->A, tmp->A_csc, data->X, tmp->alpha, tmp->beta_old,
                          para->rho, tmp->L);
                dmat_waxpby(data->m, 1.0/para->tau, tmp->alpha1, .0, tmp->alpha1, tmp->alpha1);
                dmat_waxpby(data->m, 1.0/para->tau, tmp->alpha2, .0, tmp->alpha2, tmp->alpha2);
                dmat_waxpby(data->m, 1.0/para->tau, tmp->alpha_bar, .0, tmp->alpha_bar, tmp->alpha_bar);
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
            printf("#iter:%d, #chol:%d, rho:%f, snorm:%10.6f, fun_val:%10.6f. \n", i, out->chol_num[tmp->k], para->rho[0],tmp->prires,out->fun_val[tmp->k]);
            break;
        }
        
        if (i == (para->max_num_iter[0]-1)) {
            out->fun_val[tmp->k] = fun_val(data,para,tmp);
            out->admm_iter_num[tmp->k] = i;
        }
        
    }
}

/*
void linreg_single(int* beta_dim, int* A_row_dim,
                   int* sample_size, double* Y,
                   double* X, dmatrix_csr* A,
                   double* beta_sol, double* theta,
                   double* lambda,
                   double* rho, int* num_iter,
                   bool variant, bool varying_rho, bool general_A,
                   double* eps_abs, double* eps_rel)
{
    linreg_data* data;
    parameters* para;
    tmpvars* tmp;
    ini_linreg_prob(beta_dim, A_row_dim, sample_size, Y, X, A, theta,
                    lambda, rho, num_iter, variant,
                    varying_rho, general_A, eps_abs, eps_rel, &data, &para, &tmp);
    
    if (varying_rho) // if varying rho, then factorize once at the beginning //
    {
        factorize(data->n, data->p, data->m, para->variant,
                  tmp->skinny, tmp->D, tmp->diag_mat, tmp->diag_mat_inv,
                  data->A, tmp->A_csc, data->X, tmp->alpha, tmp->beta_old,
                  para->rho, tmp->L);
    }
    
    linreg_sub(data, para, tmp);
    
    graph_update(data->m,tmp->Abeta,tmp->alpha2,data->lambda,para->rho,tmp->graph);
    UF TMP(data->p,tmp->beta);
    for (int i = 0; i < data->m; i++)
    { if (tmp->graph[i]) TMP.unite(A->idx[i], data->A->jdx[i]); }
    for (int i = 0; i < data->p; i++)
    { tmp->beta[i] = TMP.mean[TMP.root(TMP.id[i])]; }
    
    // print_dmatrix(beta,1,p);
    // save  solution
    dmat_vcopy(data->p,tmp->beta,beta_sol);
}
*/


/* solves .5 \|Y - X beta\|_2^2 + lambda \sum_{(i,i')\in E} |beta_i-beta_i'|, */
void linreg_path(int* beta_dim, int* A_row_dim,
                 int* sample_size, double* Y,
                 double* X, dmatrix_csr* A,
                     double* theta,
                 double* lambda, int* lambda_grid,
                     return_data* out,
                 double* rho, int* num_iter,
                 bool variant, bool varying_rho,
                 bool general_A, bool reporting,
                 double* eps_abs, double* eps_rel)
{
    linreg_data* data;
    parameters* para;
    tmpvars* tmp;
    ini_linreg_prob(beta_dim, A_row_dim, sample_size, Y, X, A, theta,
                    lambda, lambda_grid, rho, num_iter, variant,
                    varying_rho, general_A, reporting,
                    eps_abs, eps_rel, &data, &para, &tmp);
    
    if (varying_rho) // if varying rho, then factorize once at the beginning //
    {
        factorize(data->n, data->p, data->m, para->variant,
                  tmp->skinny, tmp->D, tmp->diag_mat, tmp->diag_mat_inv,
                  data->A, tmp->A_csc, data->X, tmp->alpha, tmp->beta_old,
                  para->rho, tmp->L);
    }
    
    for (int k = 0; k < lambda_grid[0]; k++) {
        printf("\nlambda is : %f. \n", lambda[k]);
        data->lambda = lambda + k;
        tmp->k = k;
        linreg_sub(data, out, para, tmp);
        
        
        if (!general_A)
        {
            graph_update(data->m,tmp->Abeta,tmp->alpha2,data->lambda,para->rho,tmp->graph);
            UF TMP(data->p,tmp->beta);
            for (int i = 0; i < data->m; i++)
            { if (tmp->graph[i])  TMP.unite(A->jdx[2*i], A->jdx[2*i+1]);}
            for (int i = 0; i < data->p; i++)
            { tmp->beta[i] = TMP.mean[TMP.root(TMP.id[i])]; }
        }
        /* save current solution */
        dmat_vcopy(data->p,tmp->beta,out->beta_path+k*data->p);
    }
    
    free_tmpvars_linreg(tmp);
    free(data);
    free(para);
}
/* solves .5 \|Y - X beta\|_2^2 + lambda1 \|beta\|_1 + lambda2 \sum_{(i,i')\in E} |beta_i-beta_i'|, */
void linreg_sparse_path(){
    // to be implemented ! //
}


