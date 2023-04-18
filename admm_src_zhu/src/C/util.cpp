//
//  util.cpp
//  gen_lasso
//
//  Created by Yunzhang Zhu on 9/16/14.
//  Copyright (c) 2014 OSU. All rights reserved.
//
#include <stdlib.h>
#include <stdio.h>

#include "util.h"
#include "def.h"

// UF class function definitions //
UF::UF(int N, double* beta)   {
    cnt = N;
    id = new int[N];
    sz = new int[N];
    mean = new double[N];
    for(int i=0; i<N; i++)	{
        id[i] = i;
        sz[i] = 1;
        mean[i] = beta[i];
    }
    
}
UF::~UF()	{
    delete [] id;
    delete [] sz;
    delete [] mean;
}

// Return the root of p.
int UF::root(int p)	{
    int i = p;
    while (i != id[i]) {
        // id[i] = id[id[i]];
        i = id[i];
    }
    return i;
}

// union p and q
void UF::unite(int p, int q)	{
    int i = root(p);
    int j = root(q);
    if (i == j) return;
    if (sz[i] > sz[j]) {
        id[j] = i;
        mean[i] = (mean[i]*sz[i] + mean[j]*sz[j])/(sz[i]+sz[j]);
        sz[i] += sz[j];
    }
    else {
        id[i] = j;
        mean[j] = (mean[i]*sz[i] + mean[j]*sz[j])/(sz[i]+sz[j]);
        sz[j] += sz[i];
    }
    cnt--;
}

// Are objects p and q in the same set?
bool UF::connected(int p, int q)    {
    return root(p) == root(q);
}

// Return the number of disjoint sets.
int UF::count() {
    return cnt;
}


void get_diag_general(dmatrix_csr* A, dmatrix_csc* A_csc, double* tau)
{
    int p = A->n;
    int m = A->m;
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










