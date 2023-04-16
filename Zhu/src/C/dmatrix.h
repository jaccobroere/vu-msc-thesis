#ifndef DMATRIX_H
#define DMATRIX_H
/** \file   dmatrix.h
 *  \brief  Header file for matrix and vector manipulation functions.
 */

#ifdef __cplusplus
extern "C" {
#endif
    
    
    /** \struct dmatrix
     *  \brief  Structure for both dense and sparse matrix.
     *  - if the matrix is dense,
     *      then nz is set to -1 and idx, jdx are set to NULL.
     *  - if the matrix is sparse,
     *      then nz is set to the number of non-zero elements
     *      and the coordiante is saved at idx and jdx.
     *
     *  dmatrix stores a matrix in compressed sparse row (CSR) format.
     *  For computational efficiency, it also store the row indices
     *  as well as row start indices.
     *
     *  NOTE: As for sparse matrix, the matrix indices are not modified
     *  throughout the whole program. It is possible since
     *  matrix-matrix multiplication (except diagonal matrix) is
     *  not performed in PCG mode.
     */
    typedef struct
    {
    int     m;      /**< number of rows */
    int     n;      /**< number of columns */
    int     nz;     /**< number of non-zero entries*/
    /**< nz >= 0: sparse, nz == -1: dense */
    double  *val;   /**< entry values */
    
    /* fields for sparse matrix */
    
    int     *jdx;   /**< column indices    (for both csr and coord) */
    int     *idx;   /**< row indices       (for coordinate) */
    
    } dmatrix;
    
    
    // compressed sparse row representation //
    typedef struct
    {
    int     m;      /**< number of rows */
    int     n;      /**< number of columns */
    int     nz;     /**< number of non-zero entries*/
    /**< nz >= 0: sparse, nz == -1: dense */
    double  *val;   /**< entry values */
    
    /* fields for sparse matrix */
    
    int     *jdx;   /**< column indices    (for both csr and coord) */
    int     *idx;   /**< row indices       (for coordinate) */
    int     *rdx;   /**< row start indices (for csr) */
    
    } dmatrix_csr;
    
    
        // compressed sparse column representation //
    typedef struct
    {
    int     m;      /**< number of rows */
    int     n;      /**< number of columns */
    int     nz;     /**< number of non-zero entries*/
    /**< nz >= 0: sparse, nz == -1: dense */
    double  *val;   /**< entry values */
    
    /* fields for sparse matrix */
    
    int     *jdx;   /**< column indices    (for both csc and coord) */
    int     *idx;   /**< row indices       (for coordinate) */
    int     *cdx;   /**< column start indices (for csc) */
    
    } dmatrix_csc;
    
    int dmat_norm0(const int n, const double* x);
    double dmat_norm0_surrogate(const int n, const double* x, const double* Thred);
    double dmat_norm1(const int n, const double *x);
    double dmat_norm2(const int n, const double *x);
    double dmat_norminf(const int n, const double *x);
    
    double dmat_dot(const int n, const double *x, const double *y);
    
    void dmat_vset(int n, const double val, double *dst);
    void dmat_iset(int n, const int val, int *dst);
    void dmat_vcopy(const int n, const double *src, double *dst);
    void dmat_icopy(const int n, const int *src, int *dst);
    void dmat_yexpx(const int n, const double *x, double *y);
    void dmat_ysqrtx(const int n, const double *x, double *y);
    void dmat_yinvx(const int n, const double *x, double *y);
    
    void dmat_waxpby(int n, double alpha, const double *x, double beta,
                     const double *y, double *w);
    
    void dmat_elemprod(const int n, const double *x, const double *y, double *z);
    void dmat_ideleprod(const int n, const int *x, double *y);
    void dmat_elemdivi(const int n, const double *x, const double *y, double *z);
    double dmat_xAx(int n, double *tmp, const double *A, const double *x);
    
        // dense matrix-matrix(vector) multiplication //
    void dmat_yAx(int m,int n, const double *A, const double *x, double *y);
    void dmat_yATx(int m, int n, const double *A, const double *x, double *y);
    void dmat_C_ATB(int m, int n1, int n2, double* A, double* B, double* C);
    void dmat_B_ATA(int m, int n, double *A, double *B);
    void dmat_B_AAT(int m, int n, double *A, double *B);
    void dmat_B_ATDA(int m, int n, double *A, double* D, double *B);
    void dmat_B_ADAT(int m, int n, double *A, double* D, double *B);
    
        //sparse matrix vector multiplication //
    void dmat_ySx(const dmatrix_csr *A, const double *x, double *y);
    void dmat_ySTx(const dmatrix_csr *A, const double *x, double *y);
    void dmat_yabsSx(const dmatrix_csr *A, const double *x, double *y);
    void dmat_yabsSTx(const dmatrix_csr *A, const double *x, double *y);
    
    void dmat_diagadd(dmatrix *M, const double *d);
    void dmat_colavg(const dmatrix *M, double *y);
    
    void dmat_potrs(double *A, const int m, double *b);
    void dmat_posv(double *A, const int m, double *b);
    void dmat_potrf(double *A, const int m);
    void eigen_decomp(int n, double* X, double *eigvec, double *eigval);
    double dmat_det(int n, const double *A);
    
    void dmat_to_csr(const dmatrix* M, dmatrix_csr** dst);
    void dmat_to_csc(const dmatrix* M, dmatrix_csc** dst);
    void dmatcsr_to_csc(const dmatrix_csr* M, dmatrix_csc** dst);
    void dmatcsr_rank1_update(dmatrix_csr *data_csr, const double *a, const double* b, double w);
    void dmatcsc_rank1_update(dmatrix_csc *data_csc, const double *a, const double* b, double w);
    void dmat_new_sparse(dmatrix** M, const int m, const int n, const int nz, int *idx, int *jdx, double* val);
    void dmat_new_sparse_csr(dmatrix_csr** M, const int m, const int n, const int nz, int *idx, int *jdx, int *rdx, double* val);
    void dmat_free(dmatrix* M);
    void dmat_csr_free(dmatrix_csr* M);
    void dmat_csc_free(dmatrix_csc* M);
    void dmat_build_idx(dmatrix_csr *M);
    
    /* print for debugging */
    void dmat_vprint(const int n, const double *v);
    void print_dmatrix(const double* matrix,int m,int n);
    void print_imatrix(const int* matrix,int m,int n);
    void dmat_ATx(int *graph, int NumEdges, int p, double *x, double *y);
#ifdef __cplusplus
}
#endif

#endif /* DMATRIX_H */
