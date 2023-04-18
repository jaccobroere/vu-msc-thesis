#include <R.h>
#include <Rinternals.h>
#include <stdio.h>
#include <stdlib.h>

#include "gen_lasso.h"
#include "filtering.h"
#include "dmatrix.h"
#include "util.h"


extern "C"{
    
    
    SEXP filtering(SEXP Y_,
                SEXP val_, SEXP idx_,
                SEXP jdx_,
                SEXP Lambda_,
                SEXP m_, SEXP rho_,
                SEXP eps_, SEXP maxiter_,
                SEXP admm_type_, SEXP adaptive_varying_rho_,
                SEXP general_A_, SEXP reporting_)
    {
        double *Y; Y = REAL(Y_);
        double* Lambda = REAL(Lambda_);
        int lambda_grid = length(Lambda_);
        int p = length(Y_);
        int m = *INTEGER(m_);
        double eps = *REAL(eps_);
        double rho     = *REAL(rho_);
        int max_num_iter = *INTEGER(maxiter_);
        double theta = 1.0;
        bool admm_type = *LOGICAL(admm_type_);
        bool adaptive_varying_rho = *LOGICAL(adaptive_varying_rho_);
        bool general_A = *LOGICAL(general_A_);
        bool reporting = *LOGICAL(reporting_);
        
        // bool admm_type = true;
        dmatrix* A_src;
        dmatrix_csr* A;
        A_src = (dmatrix*) malloc(sizeof(dmatrix));
        A_src->m = m;
        A_src->n = p;
        A_src->nz = length(val_);
        A_src->val = REAL(val_);
        
        A_src->idx = (int*) malloc(sizeof(int)*(A_src->nz));
        A_src->jdx = (int*) malloc(sizeof(int)*(A_src->nz));
        // C uses zero-based indexing while R uses one-based indexing //
        for (int i = 0; i < A_src->nz; i++ ) {
            A_src->jdx[i] = INTEGER(jdx_)[i] - 1;
            A_src->idx[i] = INTEGER(idx_)[i] - 1;
        }
        
        dmat_to_csr(A_src, &A);
        free(A_src->idx); free(A_src->jdx); free(A_src);
        
        // construct output //
        return_data* out = (return_data*) malloc(sizeof(return_data));
        
        SEXP Rbeta_path = PROTECT(allocVector(REALSXP, p*lambda_grid));
        SEXP Rchol_num = PROTECT(allocVector(INTSXP, lambda_grid));
        SEXP Rfun_val = PROTECT(allocVector(REALSXP, lambda_grid));
        SEXP Radmm_iter_num = PROTECT(allocVector(INTSXP, lambda_grid));
        SEXP Rfun_path  = PROTECT(allocVector(REALSXP, (reporting)? (max_num_iter*lambda_grid):0));
        
        out->beta_path = REAL(Rbeta_path);
        out->fun_path  = REAL(Rfun_path);
        out->fun_val   = REAL(Rfun_val);
        out->chol_num  = INTEGER(Rchol_num);
        out->admm_iter_num = INTEGER(Radmm_iter_num);
        
        filtering(&p, &m, Y, A->idx, A->jdx, A->rdx, &A->nz,
                  A->val,  &theta, Lambda, &lambda_grid,
                  &rho, &max_num_iter, &eps, admm_type,
                  general_A, reporting, adaptive_varying_rho, out);
        
        // free(A->idx); free(A->jdx); free(A->rdx);
        free(A);
        
        SEXP res = PROTECT(allocVector(VECSXP,5));
        SET_VECTOR_ELT(res, 0, Rbeta_path);
        SET_VECTOR_ELT(res, 1, Rfun_path);
        SET_VECTOR_ELT(res, 2, Rchol_num);
        SET_VECTOR_ELT(res, 3, Rfun_val);
        SET_VECTOR_ELT(res, 4, Radmm_iter_num);
        SEXP names = PROTECT(allocVector(STRSXP,5));
        SET_STRING_ELT(names, 0, mkChar("beta_path"));
        SET_STRING_ELT(names, 1, mkChar("fun_path"));
        SET_STRING_ELT(names, 2, mkChar("chol_num"));
        SET_STRING_ELT(names, 3, mkChar("fun_val"));
        SET_STRING_ELT(names, 4, mkChar("admm_iter_num"));
        setAttrib(res, R_NamesSymbol, names);
        
        UNPROTECT(7);
        return res;
    }
    
    
    SEXP linreg(SEXP Y_, SEXP X_,
                SEXP val_, SEXP idx_,
                SEXP jdx_, SEXP Lambda_graph_,
                    SEXP Lambda_sparse_,
                SEXP p_, SEXP m_, SEXP rho_,
                SEXP eps_abs_, SEXP eps_rel_,
                SEXP maxiter_,
                SEXP special_admm_,
                SEXP linearized_admm_,
                SEXP standard_admm_,
                SEXP adaptive_varying_rho_,
                SEXP nonadaptive_varying_rho_,
                SEXP nonvarying_rho_,
                SEXP general_A_, SEXP graph_only_,
                    SEXP graph_sparse_,
                    SEXP reporting_)
    {
        double eps_abs = *REAL(eps_abs_);
        double eps_rel = *REAL(eps_rel_);
        double rho     = *REAL(rho_);
        int max_num_iter = *INTEGER(maxiter_);
        double theta = 1.0;
        
        bool special = *LOGICAL(special_admm_);
        bool linearized = *LOGICAL(linearized_admm_);
        bool standard = *LOGICAL(standard_admm_);
        bool adaptive_varying_rho = *LOGICAL(adaptive_varying_rho_);
        bool nonadaptive_varying_rho = *LOGICAL(nonadaptive_varying_rho_);
        bool nonvarying_rho = *LOGICAL(nonvarying_rho_);
        bool general_A = *LOGICAL(general_A_);
        bool graph_only = *LOGICAL(graph_only_);
        bool graph_sparse = *LOGICAL(graph_sparse_);
        
        
        algorithm_type algorithm;
        if (special)    { algorithm = special_admm;    }
        if (linearized) { algorithm = linearized_admm; }
        if (standard)   { algorithm = standard_admm;   }
        penalty_types pen_types;
        if (general_A) pen_types = general_L1;
        if (graph_only) pen_types = graph_L1;
        if (graph_sparse) pen_types = graph_sparse_L1;
        varying_rho_types varying_rho_scheme;
        if (nonvarying_rho) varying_rho_scheme = constant;
        if (nonadaptive_varying_rho) varying_rho_scheme = non_adaptive;
        if (adaptive_varying_rho) varying_rho_scheme = adaptive;
        
        
        bool reporting = *LOGICAL(reporting_);
        // bool admm_type = true;
        int n = length(Y_);
        int p = *INTEGER(p_);
        int m = *INTEGER(m_);
        double *Y, *X;
        Y = REAL(Y_);
        X = REAL(X_);
		dmatrix* A_src; 
        dmatrix_csr* A;
		A_src = (dmatrix*) malloc(sizeof(dmatrix));
        A_src->m = m;
        A_src->n = p;
        A_src->nz = length(val_);
        A_src->val = REAL(val_);
        
        A_src->idx = (int*) malloc(sizeof(int)*(A_src->nz));
        A_src->jdx = (int*) malloc(sizeof(int)*(A_src->nz));
        // C uses zero-based indexing while R uses one-based indexing //
        for (int i = 0; i < A_src->nz; i++ ) {
			A_src->jdx[i] = INTEGER(jdx_)[i] - 1;
			A_src->idx[i] = INTEGER(idx_)[i] - 1;
		}
        
        dmat_to_csr(A_src, &A);
        free(A_src->idx); free(A_src->jdx); free(A_src);
        
        double* Lambda_graph = REAL(Lambda_graph_);
        int lambda_graph_grid = length(Lambda_graph_);
        double* Lambda_sparse = REAL(Lambda_sparse_);
        int lambda_sparse_grid = length(Lambda_sparse_);
        int lambda_grid = lambda_sparse_grid*lambda_graph_grid;
        
        // construct output //
        return_data* out = (return_data*) malloc(sizeof(return_data));
        
        SEXP Rbeta_path = PROTECT(allocVector(REALSXP, p*lambda_grid));
        SEXP Rchol_num = PROTECT(allocVector(INTSXP, lambda_grid));
        SEXP Rfun_val = PROTECT(allocVector(REALSXP, lambda_grid));
        SEXP Radmm_iter_num = PROTECT(allocVector(INTSXP, lambda_grid));
        SEXP Rfun_path  = PROTECT(allocVector(REALSXP, (reporting)? (max_num_iter*lambda_grid):0));
        
        out->beta_path = REAL(Rbeta_path);
        out->fun_path  = REAL(Rfun_path);
        out->fun_val   = REAL(Rfun_val);
        out->chol_num  = INTEGER(Rchol_num);
        out->admm_iter_num = INTEGER(Radmm_iter_num);
        
        linreg_path(&p, &m, &n, &pen_types, &varying_rho_scheme,
                    Y, X, A, &theta, Lambda_graph, &lambda_graph_grid,
                    Lambda_sparse, &lambda_sparse_grid,
                    out, &max_num_iter, &algorithm, reporting,
                     &rho, &eps_abs, &eps_rel);
        
        // free(A->idx); free(A->jdx); free(A->rdx);
        free(A);
        
        SEXP res = PROTECT(allocVector(VECSXP,5));
        SET_VECTOR_ELT(res, 0, Rbeta_path);
        SET_VECTOR_ELT(res, 1, Rfun_path);
        SET_VECTOR_ELT(res, 2, Rchol_num);
        SET_VECTOR_ELT(res, 3, Rfun_val);
        SET_VECTOR_ELT(res, 4, Radmm_iter_num);
        SEXP names = PROTECT(allocVector(STRSXP,5));
        SET_STRING_ELT(names, 0, mkChar("beta_path"));
        SET_STRING_ELT(names, 1, mkChar("fun_path"));
        SET_STRING_ELT(names, 2, mkChar("chol_num"));
        SET_STRING_ELT(names, 3, mkChar("fun_val"));
        SET_STRING_ELT(names, 4, mkChar("admm_iter_num"));
        setAttrib(res, R_NamesSymbol, names);
        
        UNPROTECT(7);
        return res;
    }
    
    
    SEXP linreg_v2(SEXP Y_, SEXP X_,
                SEXP val_, SEXP idx_,
                SEXP jdx_, SEXP Lambda_pairs_,
                SEXP p_, SEXP m_, SEXP rho_,
                SEXP eps_abs_, SEXP eps_rel_,
                SEXP maxiter_,
                   SEXP diag_mat_,
                   SEXP diag_mat_given_,
                SEXP special_admm_,
                SEXP linearized_admm_,
                SEXP standard_admm_,
                SEXP adaptive_varying_rho_,
                SEXP nonadaptive_varying_rho_,
                SEXP nonvarying_rho_,
                SEXP general_A_, SEXP graph_only_,
                SEXP graph_sparse_,
                SEXP reporting_)
    {
        double eps_abs = *REAL(eps_abs_);
        double eps_rel = *REAL(eps_rel_);
        double rho     = *REAL(rho_);
        double diag_mat = *REAL(diag_mat_);
        int max_num_iter = *INTEGER(maxiter_);
        double theta = 1.0;
        
        bool special = *LOGICAL(special_admm_);
        bool linearized = *LOGICAL(linearized_admm_);
        bool standard = *LOGICAL(standard_admm_);
        bool adaptive_varying_rho = *LOGICAL(adaptive_varying_rho_);
        bool nonadaptive_varying_rho = *LOGICAL(nonadaptive_varying_rho_);
        bool nonvarying_rho = *LOGICAL(nonvarying_rho_);
        bool general_A = *LOGICAL(general_A_);
        bool graph_only = *LOGICAL(graph_only_);
        bool graph_sparse = *LOGICAL(graph_sparse_);
        bool diag_mat_given = *LOGICAL(diag_mat_given_);
        
        
        algorithm_type algorithm;
        if (special)    { algorithm = special_admm;    }
        if (linearized) { algorithm = linearized_admm; }
        if (standard)   { algorithm = standard_admm;   }
        penalty_types pen_types;
        if (general_A) pen_types = general_L1;
        if (graph_only) pen_types = graph_L1;
        if (graph_sparse) pen_types = graph_sparse_L1;
        varying_rho_types varying_rho_scheme;
        if (nonvarying_rho) varying_rho_scheme = constant;
        if (nonadaptive_varying_rho) varying_rho_scheme = non_adaptive;
        if (adaptive_varying_rho) varying_rho_scheme = adaptive;
        
        
        bool reporting = *LOGICAL(reporting_);
        // bool admm_type = true;
        int n = length(Y_);
        int p = *INTEGER(p_);
        int m = *INTEGER(m_);
        double *Y, *X;
        Y = REAL(Y_);
        X = REAL(X_);
        dmatrix* A_src;
        dmatrix_csr* A;
        A_src = (dmatrix*) malloc(sizeof(dmatrix));
        A_src->m = m;
        A_src->n = p;
        A_src->nz = length(val_);
        A_src->val = REAL(val_);
        
        A_src->idx = (int*) malloc(sizeof(int)*(A_src->nz));
        A_src->jdx = (int*) malloc(sizeof(int)*(A_src->nz));
        // C uses zero-based indexing while R uses one-based indexing //
        for (int i = 0; i < A_src->nz; i++ ) {
            A_src->jdx[i] = INTEGER(jdx_)[i] - 1;
            A_src->idx[i] = INTEGER(idx_)[i] - 1;
        }
        
        dmat_to_csr(A_src, &A);
        free(A_src->idx); free(A_src->jdx); free(A_src);
        
        double* Lambda_pairs = REAL(Lambda_pairs_);
        int lambda_grid = length(Lambda_pairs_) / 2;
        // int lambda_grid = lambda_sparse_grid*lambda_graph_grid;
        
        // construct output //
        return_data* out = (return_data*) malloc(sizeof(return_data));
        
        SEXP Rbeta_path = PROTECT(allocVector(REALSXP, p*lambda_grid));
        SEXP Rchol_num = PROTECT(allocVector(INTSXP, lambda_grid));
        SEXP Rfun_val = PROTECT(allocVector(REALSXP, lambda_grid));
        SEXP Radmm_iter_num = PROTECT(allocVector(INTSXP, lambda_grid));
        SEXP Rfun_path  = PROTECT(allocVector(REALSXP, (reporting)? (max_num_iter*lambda_grid):0));
        
        out->beta_path = REAL(Rbeta_path);
        out->fun_path  = REAL(Rfun_path);
        out->fun_val   = REAL(Rfun_val);
        out->chol_num  = INTEGER(Rchol_num);
        out->admm_iter_num = INTEGER(Radmm_iter_num);
        
        linreg_path_v2(&p, &m, &n, &varying_rho_scheme,
                    Y, X, A, &theta, Lambda_pairs, &lambda_grid,
                    out, &max_num_iter, &algorithm, &diag_mat,
                       diag_mat_given, reporting,
                    &rho, &eps_abs, &eps_rel);
        
        // free(A->idx); free(A->jdx); free(A->rdx);
        free(A);
        
        SEXP res = PROTECT(allocVector(VECSXP,5));
        SET_VECTOR_ELT(res, 0, Rbeta_path);
        SET_VECTOR_ELT(res, 1, Rfun_path);
        SET_VECTOR_ELT(res, 2, Rchol_num);
        SET_VECTOR_ELT(res, 3, Rfun_val);
        SET_VECTOR_ELT(res, 4, Radmm_iter_num);
        SEXP names = PROTECT(allocVector(STRSXP,5));
        SET_STRING_ELT(names, 0, mkChar("beta_path"));
        SET_STRING_ELT(names, 1, mkChar("fun_path"));
        SET_STRING_ELT(names, 2, mkChar("chol_num"));
        SET_STRING_ELT(names, 3, mkChar("fun_val"));
        SET_STRING_ELT(names, 4, mkChar("admm_iter_num"));
        setAttrib(res, R_NamesSymbol, names);
        
        UNPROTECT(7);
        return res;
    }
    
}