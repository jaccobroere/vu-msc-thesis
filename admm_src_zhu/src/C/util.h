#ifndef UTIL_H
#define UTIL_H
//
//  util.h
//  gen_lasso
//
//  Created by Yunzhang Zhu on 9/16/14.
//  Copyright (c) 2014 OSU. All rights reserved.
//

#include "dmatrix.h"

class UF {
public:
    int *id, cnt, *sz;
    double* mean;
    UF(int,double*);
    ~UF();
    int root(int);
    void unite(int, int);
    bool connected(int, int);
    int count(void);
};

void get_diag_general(dmatrix_csr* A, dmatrix_csc* A_csc, double* tau);

typedef struct
{
    int *chol_num, *admm_iter_num;
    double *fun_val, *fun_path, *beta_path;
    
} return_data;

#endif