/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 et: */
/*
   IGraph library.
   Copyright (C) 2007-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard street, Cambridge, MA 02139 USA

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#ifndef BLAS_INTERNAL_H
#define BLAS_INTERNAL_H

/* Note: only files calling the BLAS routines directly need to
   include this header.
*/

#include "config.h"
#include "igraph_types.h"

#ifndef INTERNAL_BLAS
#define igraphdaxpy_ daxpy
#define igraphdger_ dger
#define igraphdcopy_ dcopy
#define igraphdscal_ dscal
#define igraphdswap_ dswap
#define igraphdgemm_ dgemm
#define igraphdgemv_ dgemv
#define igraphddot_ ddot
#define igraphdnrm2_ dnrm2
#define igraphlsame_ lsame
#define igraphdrot_ drot
#define igraphidamax_ idamax
#define igraphdtrmm_ dtrmm
#define igraphdasum_ dasum
#define igraphdtrsm_ dtrsm
#define igraphdtrsv_ dtrsv
#define igraphdnrm2_ dnrm2
#endif

int igraphdgemv_(char *trans, int *m, int *n, igraph_real_t *alpha,
                 igraph_real_t *a, int *lda, igraph_real_t *x, int *incx,
                 igraph_real_t *beta, igraph_real_t *y, int *incy);

int igraphdgemm_(char *transa, char *transb, int *m, int *n, int *k,
                 double *alpha, double *a, int *lda, double *b, int *ldb,
                 double *beta, double *c__, int *ldc);

double igraphdnrm2_(int *n, double *x, int *incx);

#endif
