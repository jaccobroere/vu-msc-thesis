/** \file   dmatrix.c
 *  \brief  Source file for matrix and vector manipulation functions.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

#ifndef Rpackage
#include "blas.h"
#include "lapack.h"
#else
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#endif

#include "def.h"
#include "dmatrix.h"

/****************************************************************************/
/*                                                                          */
/*                           VECTOR OPERATIONS                              */
/*                                                                          */
/****************************************************************************/

int dmat_norm0(const int n, const double* x)
{
    int i, nz = 0;
    for (i = 0; i < n; i++)
        if (abs(x[i]) > 1e-10) nz++;
    
    return nz;
}

double dmat_norm0_surrogate(const int n, const double* x, const double* Thred)
{
    double thred = Thred[0];
    double part = 0.0;
    int i;
    for (i = 0; i < n; i++) {
        part += min(abs(x[i]), thred);
    }
    return part;
}


/** \brief \f$ \|x\|_1 \f$
 *
 *  Returns 1-norm of a vector x.
 *
 *  @param  n   length of a vector x.
 *  @param  x   pointer to a vector x.
 *  @return     result.
 */


double dmat_norm1(const int n, const double *x)
{
    return F77_CALL(dasum)(&n, x, &ione);
}


/** \brief \f$ \|x\|_2 \f$
 *
 *  Returns 2-norm of a vector x.
 *
 *  @param  n   length of a vector x.
 *  @param  x   pointer to a vector x.
 *  @return     result.
 */
double dmat_norm2(const int n, const double *x)
{
    return F77_CALL(dnrm2)(&n, x, &ione);
}


/** \brief \f$ \|x\|_{\infty} \f$
 *
 *  Returns infinity-norm of a vector x
 *
 *  @param  n   length of a vector x
 *  @param  x   pointer to a vector x
 *  @return     result.
 */
double dmat_norminf(const int n, const double *x)
{
    return fabs(x[F77_CALL(idamax)(&n, x, &ione)-1]);
}


/** \brief \f$ x^Ty \f$
 *
 *  Returns dot product of a vector x and y.
 *
 *  @param  n   length of a vector x.
 *  @param  x   pointer to a vector x.
 *  @param  y   pointer to a vector y.
 *  @return     result.
 */
double dmat_dot(const int n, const double *x, const double *y)
{
    return F77_CALL(ddot)(&n, x, &ione, y, &ione);
}


/** \brief \f$ \mbox{dst}_i \leftarrow \mbox{val} \f$
 *
 *  Sets all the elements of a vector with a constant value.
 *
 *  @param  n   length of a vector.
 *  @param  val constant value to set.
 *  @param  dst pointer to a vector.
 */
void dmat_vset(int n, const double val, double *dst)
{
    while (n-- != 0)
        *dst++ = val;
}

void dmat_iset(int n, const int val, int *dst)
{
    while (n-- != 0)
        *dst++ = val;
}


/** \brief \f$ \mbox{dst} \leftarrow \mbox{src} \f$
 *
 *  Copies a vector.
 *
 *  @param  n   length of vectors.
 *  @param  src pointer to a source vector.
 *  @param  dst pointer to a destination vector.
 */
void dmat_vcopy(const int n, const double *src, double *dst)
{
    F77_CALL(dcopy)(&n, src, &ione, dst, &ione);
}

void dmat_icopy(const int n, const int *src, int *dst)
{
    memcpy(dst, src, sizeof(int)*n);
}


/** \brief \f$ y_i = \exp(x_i) \f$
 *
 *  Computes elementwise exp() of a vector.
 *
 *  @param  n   length of vectors.
 *  @param  x   pointer to a source vector.
 *  @param  y   pointer to a destination vector.
 */
void dmat_yexpx(const int n, const double *x, double *y)
{
#ifdef MKL
    vdExp(n, x, y);
#else
#ifdef ACML
    vrda_exp(n, x, y);
#else
    const double *xi;
    double *yi;
    xi = x+n-1;
    yi = y+n-1;
    do {
        *yi-- = exp(*xi--);
    } while (xi >= x);
#endif
#endif
}


/** \brief \f$ y_i = x_i^{1/2} \f$
 *
 *  Computes elementwise sqrt() of a vector.
 *
 *  @param  n   length of vectors.
 *  @param  x   pointer to a source vector.
 *  @param  y   pointer to a destination vector.
 */
void dmat_ysqrtx(const int n, const double *x, double *y)
{
#ifdef MKL
    vdSqrt(n, x, y);
#else
    const double *xi;
    double *yi;
    xi = x+n-1;
    yi = y+n-1;
    do {
        *yi-- = sqrt(*xi--);
    } while (xi >= x);
#endif
}


/** \brief \f$ y_i = 1/x_i \f$
 *
 *  Computes elementwise inv() of a vector.
 *
 *  @param  n   length of vectors.
 *  @param  x   pointer to a source vector.
 *  @param  y   pointer to a destination vector.
 */
void dmat_yinvx(const int n, const double *x, double *y)
{
#ifdef MKL
    vdInv(n, x, y);
#else
    const double *xi;
    double *yi;
    xi = x+n-1;
    yi = y+n-1;
    do {
        *yi-- = 1/(*xi--);
    } while (xi >= x);
#endif
}


/** \brief \f$ w = \alpha*x + \beta*y \f$
 *
 *  Computes weighted vector sum.
 *  - w = -x
 *  - w = alpha*x
 *  - w = -x + y
 *  - w = alpha*x + y
 *  - w = -x - y
 *  - w = alpha*x - y
 *  - w = alpha*x + beta*y
 *
 *  @param  n       length of vectors.
 *  @param  alpha   constant
 *  @param  x       pointer to a vector.
 *  @param  beta    constant
 *  @param  y       pointer to a vector.
 *  @param  w       pointer to a result vector.
 */
void dmat_waxpby(int n, double alpha, const double *x, double beta,
                 const double *y, double *w)

{
#if 1
    if (w != x && w != y)
    {
        dmat_vset(n, 0, w);
        F77_CALL(daxpy)(&n, &alpha, x, &ione, w, &ione);
        F77_CALL(daxpy)(&n, &beta , y, &ione, w, &ione);
    }
    else if (w == x && w == y)
    {
        double tmp;
        tmp = alpha+beta;
        F77_CALL(dscal)(&n, &tmp, w, &ione);
    }
    else if (w == x /*&& w != y */)
    {
        if (alpha != 1.0) F77_CALL(dscal)(&n, &alpha, w, &ione);
        if (beta  != 0.0) F77_CALL(daxpy)(&n, &beta , y, &ione, w, &ione);
    }
    else /* if (w == y && w != x ) */
    {
        if (beta  != 1.0) F77_CALL(dscal)(&n, &beta , w, &ione);
        if (alpha != 0.0) F77_CALL(daxpy)(&n, &alpha, x, &ione, w, &ione);
    }
#else
    int i;
    
    if (beta == 0.0)
    {
        if (alpha == -1.0)
        {
            for (i = 0; i < n; i++)
                w[i] = -x[i];
        }
        else
        {
            for (i = 0; i < n; i++)
                w[i] = alpha*x[i];
        }
    }
    else if (beta == 1.0)
    {
        if (alpha == -1.0)
        {
            for (i = 0; i < n; i++)
                w[i] = -x[i] + y[i];
        }
        else
        {
            for (i = 0; i < n; i++)
                w[i] = alpha*x[i] + y[i];
        }
    }
    else if (beta == -1.0)
    {
        if (alpha == -1.0)
        {
            for (i = 0; i < n; i++)
                w[i] = -x[i] - y[i];
        }
        else
        {
            for (i = 0; i < n; i++)
                w[i] = alpha*x[i] - y[i];
        }
    }
    else
    {
        for (i = 0; i < n; i++)
            w[i] = alpha*x[i] + beta*y[i];
    }
#endif
}


/**  \brief \f$ z_i = x_i*y_i \f$
 *
 *  Computes elementwise product of vectors.
 *
 *  NOTE: x = x.*y is not allowed, i.e., z should not be the same with x.
 *
 *  @param  n   length of a vector x.
 *  @param  x   pointer to a vector x.
 *  @param  y   pointer to a vector y.
 *  @param  z   pointer to a result vector z.
 */
void dmat_elemprod(const int n, const double *x, const double *y, double *z)
{
    if (y != z) {
        F77_CALL(dcopy)(&n, y, &ione, z, &ione);
    }
    F77_CALL(dtbmv)("U", "N", "N", &n, &izero, x, &ione, z, &ione);
}

/*
 x is a zero-one vector
 calculate y = y.*x
 */

void dmat_ideleprod(const int n, const int *x, double *y)
{
    int i;
    for (i = 0; i < n; i++) {
        if (x[i] == 0) {
            y[i] = 0;
        }
    }
}

/** \brief \f$ z_i = x_i/y_i \f$
 *
 *  Computes elementwise division of vectors.
 *
 *  NOTE: y = x./y is not allowed, i.e., w should not be the same with y.
 *
 *  @param  n   length of a vector x.
 *  @param  x   pointer to a vector x.
 *  @param  y   pointer to a vector y.
 *  @param  z   pointer to a result vector z.
 */
void dmat_elemdivi(const int n, const double *x, const double *y, double *z)
{
#ifdef MKL
    vdDiv(n, x, y, z);
#else
    if (x != z) F77_CALL(dcopy)(&n, x, &ione, z, &ione);
    F77_CALL(dtbsv)("U", "N", "N", &n, &izero, y, &ione, z, &ione);
#endif
}


/****************************************************************************/
/*                                                                          */
/*                           MATRIX OPERATIONS                              */
/*                                                                          */
/****************************************************************************/

/*
 This function calculates the eigenvalues and eigenvectors of
 the n*n symmetric matrix X.
 The matrices have to be in Fortran vector format.
 The eigenvectors will be put columnwise in the n*n matrix eigvec,
 where the corresponding eigenvalues will be put in the vector
 eigval (length n of course). Only the lower triangle of the matrix
 X is used. The content of X is not changed.
 
 This function first queries the Lapack routines for optimal workspace
 sizes. These memoryblocks are then allocated and the decomposition is
 calculated using the Lapack function "dsyevr". The allocated memory
 is then freed.
 */


/** \brief \f$ y = Ax \f$
 *
 *  Computes matrix-vector product.
 *
 *  @param  A   pointer to a matrix.
 *  @param  x   pointer to a vector x.
 *  @param  y   pointer to a result vector y.
 */
void dmat_ySx(const dmatrix_csr *A, const double *x, double *y)
{
    int k, l;
    int m, n, nz;
    int *jdx, *rdx;
    double *val;
    double tmp;
    
    m   = A->m;
    n   = A->n;
    nz  = A->nz;
    val = A->val;
    jdx = A->jdx;
    rdx = A->rdx;
    
    if (nz >= 0)    /* sparse matrix */
    {
        for (k = 0; k < m ; k++)
        {
            tmp = 0.0;
            for (l = rdx[k]; l < rdx[k+1]; l++)
            {
                tmp += val[l]*x[jdx[l]];
            }
            y[k] = tmp;
        }
    }
    else            /* dense matrix */
    {
        F77_CALL(dgemv)("T",&n,&m,&done,val,&n,x,&ione,&dzero,y,&ione);
    }
}


/** \brief \f$ y = A^Tx \f$
 *
 *  Computes transposed matrixr-vector product.
 *
 *  @param  A   pointer to a matrix.
 *  @param  x   pointer to a vector x.
 *  @param  y   pointer to a result vector y.
 */
void dmat_ySTx(const dmatrix_csr *A, const double *x, double *y)
{
    int k, l;
    int m, n, nz;
    int *jdx, *rdx;
    double *val;
    double tmp;
    
    m   = A->m;
    n   = A->n;
    nz  = A->nz;
    val = A->val;
    jdx = A->jdx;
    rdx = A->rdx;
    
    dmat_vset(n, 0, y);
    if (nz >= 0)    /* sparse matrix */
    {
        for (k = 0; k < m ; k++)
        {
            tmp = x[k];
            for (l = rdx[k]; l < rdx[k+1]; l++)
            {
                y[jdx[l]] += val[l]*tmp;
            }
        }
    }
    else            /* dense matrix */
    {
        F77_CALL(dgemv)("N",&n,&m,&done,val,&n,x,&ione,&dzero,y,&ione);
    }
}


/** \brief \f$ y = |A|x \f$
 *
 *  Computes matrix-vector product.
 *
 *  @param  A   pointer to a matrix.
 *  @param  x   pointer to a vector x.
 *  @param  y   pointer to a result vector y.
 */
void dmat_yabsSx(const dmatrix_csr *A, const double *x, double *y)
{
    int k, l, m;
    int *jdx, *rdx;
    double *val;
    double tmp;
    
    m   = A->m;
    val = A->val;
    jdx = A->jdx;
    rdx = A->rdx;
    for (k = 0; k < m ; k++)
    {
        tmp = 0.0;
        for (l = rdx[k]; l < rdx[k+1]; l++)
        {
            tmp += abs(val[l])*x[jdx[l]];
        }
        y[k] = tmp;
    }
}


/** \brief \f$ y = |A^T|x \f$
 *
 *  Computes transposed matrixr-vector product.
 *
 *  @param  A   pointer to a matrix.
 *  @param  x   pointer to a vector x.
 *  @param  y   pointer to a result vector y.
 */
void dmat_yabsSTx(const dmatrix_csr *A, const double *x, double *y)
{
    int k, l;
    int m, n;
    int *jdx, *rdx;
    double *val;
    double tmp;
    
    m   = A->m;
    n   = A->n;
    val = A->val;
    jdx = A->jdx;
    rdx = A->rdx;
    
    dmat_vset(n, 0, y);
    for (k = 0; k < m ; k++)
    {
        tmp = x[k];
        for (l = rdx[k]; l < rdx[k+1]; l++)
        {
            y[jdx[l]] += abs(val[l])*tmp;
        }
    }
}




/** \brief \f$ M = M+\mbox{diag}(d) \f$
 *
 *  Adds diagonal entries to a matrix.
 *
 *  @param  M       pointer to a matrix.
 *  @param  d       pointer to a diagonal vector.
 */
void dmat_diagadd(dmatrix *M, const double *d)
{
    int m, n, nz;
    m  = M->m;
    n  = M->n;
    nz = M->nz;
    if (nz >= 0)
    {
        /* TO BE IMPLEMENTED */
        printf("dmat_diagadd: sparse case should be implemented\n");
        exit(0);
    }
    else
    {
        int i;
        double *val;
        val = M->val;
        for (i = 0; i < min(m,n); i++)
        {
            val[(n+1)*i] += d[i];
        }
    }
}




/** \brief \f$ y_j = (1/m)\sum_{i} M_{ij} \f$
 *
 *  Computes column average of a matrix.
 *
 *  @param  M       pointer to a matrix.
 *  @param  y       pointer to a result vector (column average).
 */
void dmat_colavg(const dmatrix *M, double *y)
{
    int i, j, m, n, nz;
    int *jdx;
    double *val;
    n   = M->n;
    m   = M->m;
    nz  = M->nz;
    val = M->val;
    jdx = M->jdx;
    
    if (nz >= 0)
    {
        double tmp;
        dmat_vset(n, 0, y);
        for (i = 0; i < nz; i++)
        {
            y[jdx[i]] += val[i];
        }
        tmp = 1.0/m;
        F77_CALL(dscal)(&n, &tmp, y, &ione);
    }
    else    /* dense matrix */
    {
        for (j = 0; j < n; j++) /* col */
        {
            double sum = 0.0;
            for (i = 0; i < m; i++)
                sum += val[i*n+j];
            
            y[j] = sum/m;
        }
    }
}



/** \brief \f$ Ax = b \f$
 *
 *  Computes the solution of a linear system Ax = b via Cholesky method.
 *  Wrapper function to dpotrs.
 *  That is, A should be a Cholesky factor computed by potrf.
 *
 *  @param  A       pointer to a matrix (Cholesky factor).
 *  @param  b       pointer to a vector.
 */
void dmat_potrs(double *A, const int m, double *b)
{
    int one = 1;
    int info;
    F77_CALL(dpotrs)("L", &m, &one, A, &m, b, &m, &info);
}

/** \brief \f$ Ax = b \f$
 *
 *  Computes the solution of a linear system Ax = b via Cholesky method.
 *  Wrapper function to dports.
 *  That is, the cholesky factor is stored in A after computation.
 *
 *  @param  A       pointer to a matrix.
 *  @param  b       pointer to a vector.
 */
void dmat_posv(double *A, const int m, double *b)
{
    int one = 1;
    int info;
    F77_CALL(dposv)("L", &m, &one, A, &m, b, &m, &info);
}

/*
 computes the Cholesky factorization of a real symmetric
 positive definite matrix A.
 */
void dmat_potrf(double *A, const int m){
    int info;
    F77_CALL(dpotrf)("L", &m, A, &m, &info);
}

void eigen_decomp(int n, double* X, double *eigvec, double *eigval) {
    
    double *WORK;
    double abstol, WORKopt, vl, vu;
    int *IWORK;
    int numeig, sizeWORK, sizeIWORK, IWORKopt, il, iu,info;
    vl = 0.0;
    vu = 0.0;
    il = 0;
    iu = 0;
    /*  The support of the eigenvectors. We will not use this but the routine needs it  */
    int ISUPPZ[2*n];
    abstol = -1.0; // default tolerance
    
    /*  Query the Lapack routine for optimal sizes for workspace arrays  */
    sizeWORK = -1;
    sizeIWORK = -1;
    F77_CALL(dsyevr)("V", "A", "L", &n, X, &n, &vl,&vu,&il,&iu, &abstol, &numeig, eigval, eigvec, &n, ISUPPZ, &WORKopt, &sizeWORK, &IWORKopt, &sizeIWORK,&info);
    sizeWORK = (int)WORKopt;
    sizeIWORK = IWORKopt;
    WORK = (double*)malloc (sizeWORK*sizeof(double));
    IWORK = (int*)malloc (sizeIWORK*sizeof(int));
    /*  Now calculate the eigenvalues and vectors using optimal workspaces  */
    F77_CALL(dsyevr)("V", "A", "L", &n, X, &n, &vl,&vu,&il,&iu, &abstol, &numeig, eigval, eigvec, &n, ISUPPZ, WORK, &sizeWORK, IWORK, &sizeIWORK,&info);
    /*  Cleanup  */
    free((void*)(WORK)); free((void*)IWORK);
}


/** \brief \f$ y = x'Ax \f$
 *
 *  Computes quadratic forms
 *
 *  @param  A   pointer to a matrix.
 *  @param  x   pointer to a vector x.
 */
double dmat_xAx(int n, double *tmp, const double *A, const double *x)
{
    F77_CALL(dgemv)("T",&n,&n,&done,A,&n,x,&ione,&dzero,tmp,&ione);
    return dmat_dot(n, x, tmp);
}


/** \brief \f$ y = Ax \f$
 *
 *  Computes matrix-vector product.
 *
 *  @param  A   pointer to m by n matrix.
 *  @param  x   pointer to a vector x.
 *  @param  y   pointer to a result vector y.
 */
void dmat_yAx(int m, int n, const double *A, const double *x, double *y)
{
    F77_CALL(dgemv)("T",&n,&m,&done,A,&n,x,&ione,&dzero,y,&ione);
}

/** \brief \f$ y = A^Tx \f$
 *
 *  Computes transposed matrixr-vector product.
 *
 *  @param  A   pointer to m by n matrix.
 *  @param  x   pointer to a vector x.
 *  @param  y   pointer to a result vector y.
 */
void dmat_yATx(int m, int n, const double *A, const double *x, double *y)
{
    F77_CALL(dgemv)("N",&n,&m,&done,A,&n,x,&ione,&dzero,y,&ione);
}

/** \brief \f$ C = A^TB \f$
 *
 *  Computes dense matrix-transpose-matrix product.
 *
 *  @param  A       pointer to m by n1 matrix.
 *  @param  B       pointer to m by n2 matrix.
 *  @param  C       pointer to n1 by n2 matrix.
 */
void dmat_C_ATB(int m, int n1, int n2, double* A, double* B, double* C)
{
    F77_CALL(dgemm)("N","T",&n2,&n1,&m,&done,B,&n2,A,&n1,&dzero,C,&n2);
}

/** \brief \f$ B = A^TA \f$
 *
 *  Computes dense matrix-transpose-matrix product; return a
 *  row-oriented upper-diag symmetric matrix.
 *
 *  @param  A       pointer to m by n matrix.
 *  @param  B       pointer to a result matrix.
 */
void dmat_B_ATA(int m, int n, double *A, double *B)
{
    F77_CALL(dsyrk)("L","N",&n,&m,&done,A,&n,&dzero,B,&n);
}

/** \brief \f$ B = AA^T \f$
 *
 *  Computes dense matrix-transpose-matrix product; return a
 *  row-oriented upper-diag symmetric matrix.
 *
 *  @param  A       pointer to m by n matrix.
 *  @param  B       pointer to a result matrix.
 */
void dmat_B_AAT(int m, int n, double *A, double *B)
{
    F77_CALL(dsyrk)("L","T",&m,&n,&done,A,&n,&dzero,B,&m);
}


/** \brief \f$ B = A^TDA \f$
 *
 *  Computes dense matrix-transpose-matrix product.
 *
 *  @param  A       pointer to m by n matrix.
 *  @param  D       pointer to a m dimensional vector.
 *  @param  B       pointer to a result matrix.
 */
void dmat_B_ATDA(int m, int n, double *A, double* D, double *B)
{
    int i;
    double* Tmp = (double*) malloc(sizeof(double)*m*n);
    // scale rows of A by D and stores the scaled matrix in Tmp //
    for (i = 0; i < m; i++) {
        dmat_waxpby(n, sqrt(D[i]), A+i*n, 0.0, Tmp, Tmp+i*n);
    }
    dmat_B_ATA(m, n, Tmp, B);
    free(Tmp);
}

/** \brief \f$ B = ADA^T \f$
 *
 *  Computes dense matrix-transpose-matrix product.
 *
 *  @param  A       pointer to m by n matrix.
 *  @param  D       pointer to a n dimensional vector.
 *  @param  B       pointer to a result matrix.
 */
void dmat_B_ADAT(int m, int n, double *A, double* D, double *B)
{
    int i;
    double* Tmp = (double*) malloc(sizeof(double)*m*n);
    double* sqrtD = Tmp+(m-1)*n;
    dmat_ysqrtx(n, D, sqrtD);
    // scale cols of A by D and stores the scaled matrix in Tmp //
    for (i = 0; i < m; i++) {
        dmat_elemprod(n, A+i*n, sqrtD, Tmp+i*n);
    }
    dmat_B_AAT(m, n, Tmp, B);
    free(Tmp);
}


/** \brief \f$ y = det(A) \f$
 *
 *  Computes determinant of positive definite matrix
 *
 *  @param  A   pointer to a matrix.
 *  @param  n   dimension of A
 */
double dmat_det(int n, const double *A)
{
    int i;
    double* A_copy = (double*) malloc(sizeof(double)*(n+1)*n/2);
    for (i = 0; i < n; i++) {
        dmat_vcopy(i+1, A+i*n, A_copy+(i+1)*i/2);
    }
    int ipiv[n];
    double b[n];
    dmat_vset(n, 0, b);
    F77_CALL(dspsv)("U", &n, &ione, A_copy, ipiv, b, &n, &izero);
    if (izero != 0) {
        printf("determinant failed \n");
    }
    
    /*
     ** compute the determinant det = det(A)
     ** if ipiv[i] > 0, then D(i,i) is a 1x1 block diagonal
     ** if ipiv[i] = ipiv[i-1] < 0, then D(i-1,i-1),
     ** D(i-1,i), and D(i,i) form a 2x2 block diagonal
     */
    double det = 1;
    for(i = 0; i < n; i++) {
        if (ipiv[i]>0) {
            det *= A_copy[i+i*(i+1)/2];
        }
        if ((i>0) && (ipiv[i]<0) && (ipiv[i-1] == ipiv[i])) {
            det *= (A_copy[i+i*(i+1)/2]*A_copy[i-1+i*(i-1)/2]-A_copy[i-1+i*(i+1)/2]*A_copy[i-1+i*(i+1)/2]);
        }
    }
    free(A_copy);
    return det;
}

// calculate y = A^T x, where A(m,i) = 1, A(m, j) = -1
// with (i,j) being the m-th edge in the graph //
void dmat_ATx(int *graph, int NumEdges, int p, double *x, double *y){
    dmat_vset(p, 0, y);
    int m;
    for (m = 0; m < NumEdges; m++) {
        //int i = graph[2*m];
        //int j = graph[2*m+1];
        y[graph[2*m]] += x[m];
        y[graph[2*m+1]] -= x[m];
    }
}



/** \brief Build row index from csr info. (jdx and rdx)
 *
 *  -
 *
 *  @param  M       pointer to a matrix.
 */
void dmat_build_idx(dmatrix_csr *M)
{
    int i, j, m;
    int *idx, *rdx;
    
    m = M->m;
    rdx = M->rdx;
    idx = M->idx;
    
    for (i = 0; i < m; i++)
    {
        for (j = rdx[i]; j < rdx[i+1]; j++) *idx++ = i;
    }
}

void dmat_to_csr(const dmatrix* M, dmatrix_csr** dst){
    // duplicate to dmatrix_csr //
    int m, n, nz;
    int *idx, *jdx;
    double *val;
    n   = M->n;
    m   = M->m;
    nz  = M->nz;
    val = M->val;
    idx = M->idx;
    jdx = M->jdx;
    dmatrix_csr *mcp;
    mcp = (dmatrix_csr*) malloc(sizeof(dmatrix_csr));
    mcp->m = m;
    mcp->n = n;
    mcp->nz = nz;
    if (nz >= 0) {
        mcp->val = (double*) malloc(sizeof(double)*nz);
        mcp->idx = (int*) malloc(sizeof(int)*nz);
        mcp->jdx = (int*) malloc(sizeof(int)*nz);
        mcp->rdx = (int*) malloc(sizeof(int)*(m+1));
    }
    else {
        mcp->val = (double*) malloc(sizeof(double)*(m*n));
        mcp->idx = NULL;
        mcp->jdx = NULL;
        mcp->rdx = NULL;
    }
    // conversion of sparse to compressed sparse row //
    int *row_size = (int*) malloc(sizeof(int)*m);
    dmat_iset(m, 0, row_size);
    int i,k;
    for (i = 0; i < nz; i++) {
        if (idx[i] >= m) {
            printf("maximum dimention went wrong! \n");
        }
        row_size[idx[i]]++;
    }
    dmat_iset(m+1, 0, mcp->rdx);
    for (k = 0; k < m; k++) {
        mcp->rdx[k+1] = mcp->rdx[k] + row_size[k];
    }
    dmat_iset(m, 0, row_size);
    for (i = 0; i < nz; i++) {
        if ((mcp->rdx[idx[i]]+row_size[idx[i]])>=nz) {
            printf("out of bounds in dmat_to_csr!\n");
        }
        mcp->val[mcp->rdx[idx[i]]+row_size[idx[i]]] = val[i];
        mcp->idx[mcp->rdx[idx[i]]+row_size[idx[i]]] = idx[i];
        mcp->jdx[mcp->rdx[idx[i]]+row_size[idx[i]]] = jdx[i];
        row_size[idx[i]]++;
    }
    free(row_size); row_size = NULL;
    *dst = mcp;
}

void dmat_to_csc(const dmatrix* M, dmatrix_csc** dst){
    // duplicate to dmatrix_csr //
    int m, n, nz;
    int *idx, *jdx;
    double *val;
    n   = M->n;
    m   = M->m;
    nz  = M->nz;
    val = M->val;
    idx = M->idx;
    jdx = M->jdx;
    dmatrix_csc *mcp;
    mcp = (dmatrix_csc*) malloc(sizeof(dmatrix_csc));
    mcp->m = m;
    mcp->n = n;
    mcp->nz = nz;
    if (nz >= 0) {
        mcp->val = (double*) malloc(sizeof(double)*nz);
        mcp->idx = (int*) malloc(sizeof(int)*nz);
        mcp->jdx = (int*) malloc(sizeof(int)*nz);
        mcp->cdx = (int*) malloc(sizeof(int)*(n+1));
    }
    else {
        mcp->val = (double*) malloc(sizeof(double)*(m*n));
        mcp->idx = NULL;
        mcp->jdx = NULL;
        mcp->cdx = NULL;
    }
    // conversion of sparse to compressed sparse row //
    
    int *col_size = (int*) malloc(sizeof(int)*n);
    dmat_iset(n, 0, col_size);
    int i,k;
    for (i = 0; i < nz; i++) {
        if (jdx[i] >= n) {
            printf("maximum dimention went wrong! \n");
        }
        col_size[jdx[i]]++;
    }
    dmat_iset(n+1, 0, mcp->cdx);
    for (k = 0; k < n; k++) {
        mcp->cdx[k+1] = mcp->cdx[k] + col_size[k];
    }
    dmat_iset(n, 0, col_size);
    for (i = 0; i < nz; i++) {
        if ((mcp->cdx[jdx[i]]+col_size[jdx[i]])>=nz) {
            printf("out of bounds in dmat_to_csr!\n");
        }
        mcp->val[mcp->cdx[jdx[i]]+col_size[jdx[i]]] = val[i];
        mcp->idx[mcp->cdx[jdx[i]]+col_size[jdx[i]]] = idx[i];
        mcp->jdx[mcp->cdx[jdx[i]]+col_size[jdx[i]]] = jdx[i];
        col_size[jdx[i]]++;
    }
    free(col_size); col_size = NULL;
    *dst = mcp;
}

void dmatcsr_to_csc(const dmatrix_csr* M, dmatrix_csc** dst){
    // duplicate to dmatrix_csr //
    int m, n, nz;
    int *idx, *jdx;
    double *val;
    n   = M->n;
    m   = M->m;
    nz  = M->nz;
    val = M->val;
    idx = M->idx;
    jdx = M->jdx;
    dmatrix_csc *mcp;
    mcp = (dmatrix_csc*) malloc(sizeof(dmatrix_csc));
    mcp->m = m;
    mcp->n = n;
    mcp->nz = nz;
    if (nz >= 0) {
        mcp->val = (double*) malloc(sizeof(double)*nz);
        mcp->idx = (int*) malloc(sizeof(int)*nz);
        mcp->jdx = (int*) malloc(sizeof(int)*nz);
        mcp->cdx = (int*) malloc(sizeof(int)*(n+1));
    }
    else {
        mcp->val = (double*) malloc(sizeof(double)*(m*n));
        mcp->idx = NULL;
        mcp->jdx = NULL;
        mcp->cdx = NULL;
    }
    // conversion of sparse to compressed sparse row //
    
    int *col_size = (int*) malloc(sizeof(int)*n);
    dmat_iset(n, 0, col_size);
    int i,k;
    for (i = 0; i < nz; i++) {
        if (jdx[i] >= n) {
            printf("maximum dimention went wrong! \n");
        }
        col_size[jdx[i]]++;
    }
    dmat_iset(n+1, 0, mcp->cdx);
    for (k = 0; k < n; k++) {
        mcp->cdx[k+1] = mcp->cdx[k] + col_size[k];
    }
    dmat_iset(n, 0, col_size);
    for (i = 0; i < nz; i++) {
        if ((mcp->cdx[jdx[i]]+col_size[jdx[i]])>=nz) {
            printf("out of bounds in dmat_to_csr!\n");
        }
        mcp->val[mcp->cdx[jdx[i]]+col_size[jdx[i]]] = val[i];
        mcp->idx[mcp->cdx[jdx[i]]+col_size[jdx[i]]] = idx[i];
        mcp->jdx[mcp->cdx[jdx[i]]+col_size[jdx[i]]] = jdx[i];
        col_size[jdx[i]]++;
    }
    free(col_size); col_size = NULL;
    *dst = mcp;
}


// sparse matrix (row-oriented) rank-1 update: mat = mat + w * (a b^T) //
void dmatcsr_rank1_update(dmatrix_csr *data_csr, const double *a, const double* b, double w){
    int* jdx = data_csr->jdx;
    int* rdx = data_csr->rdx;
    int m = data_csr->m;
    int i, k;
    for (k = 0; k < m; k++) {
        int sample_size = rdx[k+1] - rdx[k];
        if (sample_size == 0) continue;
        for (i = 0; i < sample_size; i++)
            data_csr->val[rdx[k] + i] += (w * a[k] * b[jdx[rdx[k] + i]]);
    }
}
// sparse matrix (col-oriented) rank-1 update: mat = mat + w * (a b^T) //
void dmatcsc_rank1_update(dmatrix_csc *data_csc, const double *a, const double* b, double w){
    int* idx = data_csc->idx;
    int* cdx = data_csc->cdx;
    int n = data_csc->n;
    int i, k;
    for (k = 0; k < n; k++) {
        int sample_size = cdx[k+1] - cdx[k];
        if (sample_size == 0) continue;
        for (i = 0; i < sample_size; i++)
            data_csc->val[cdx[k] + i] += (w * b[k] * a[idx[cdx[k] + i]]);
    }
}

void dmat_new_sparse(dmatrix** M, const int m, const int n, const int nz, int *idx, int *jdx, double* val){
    dmatrix *dmat;
    dmat = (dmatrix*) malloc(sizeof(dmatrix));
    dmat->m = m;
    dmat->n = n;
    dmat->nz = nz;
    dmat->val = val;
    dmat->idx = idx;
    dmat->jdx = jdx;
    *M = dmat;
}

void dmat_new_sparse_csr(dmatrix_csr** M, const int m, const int n, const int nz, int *idx, int *jdx, int *rdx, double* val) {
    dmatrix_csr *dmat;
    dmat = (dmatrix_csr*) malloc(sizeof(dmatrix_csr));
    dmat->m = m;
    dmat->n = n;
    dmat->nz = nz;
    dmat->val = val;
    dmat->idx = idx;
    dmat->jdx = jdx;
    dmat->rdx = rdx;
    *M = dmat;
}


void dmat_csr_free(dmatrix_csr* M){
    if (M) {
        if (M->val) free(M->val);
        if (M->idx) free(M->idx);
        if (M->jdx) free(M->jdx);
        if (M->rdx) free(M->rdx);
        free(M);
    }
}

void dmat_csc_free(dmatrix_csc* M){
    if (M) {
        if (M->val) free(M->val);
        if (M->idx) free(M->idx);
        if (M->jdx) free(M->jdx);
        if (M->cdx) free(M->cdx);
        free(M);
    }
}

void dmat_free(dmatrix* M)
{
    if (M){
        if (M->val) free(M->val);
        if (M->idx) free(M->idx);
        if (M->jdx) free(M->jdx);
        free(M);
    }
}



// for debugging //
void print_dmatrix(const double* matrix,int m,int n){
    int i,j;
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            printf("%5.8f, ", matrix[i*n+j]);
        }
        printf("\n");
    }
    printf("\n");
}

void print_imatrix(const int* matrix,int m,int n){
    int i,j;
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            printf("%d ", matrix[i*n+j]);
        }
        printf("\n");
    }
    printf("\n");
}




// for debugging //



