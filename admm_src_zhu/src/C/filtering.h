//
//  filtering.h
//  gen_lasso
//
//  Created by Yunzhang Zhu on 10/24/14.
//  Copyright (c) 2014 OSU. All rights reserved.
//

#ifndef __gen_lasso__filtering__
#define __gen_lasso__filtering__

#include <stdio.h>

#include "util.h"

void filtering(int* beta_dim, int* A_row_dim,
               double* Y, int* A_idx,
               int* A_jdx, int* A_rdx, int* nz,
               double* A_val, 
               double* theta,
               double* lambda, int* lambda_grid,
               double* rho, int* num_iter,
               double* eps, bool variant,
               bool general_A, bool reporting,
               bool varying_rho,return_data* out);

#endif /* defined(__gen_lasso__filtering__) */
