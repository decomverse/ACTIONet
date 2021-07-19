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

#ifndef LAPACK_INTERNAL_H
#define LAPACK_INTERNAL_H

/* Note: only files calling the LAPACK routines directly need to
   include this header.
*/

#include "igraph_types.h"
#include "config.h"

#ifndef INTERNAL_LAPACK
    #define igraphdgeevx_   dgeevx
    #define igraphdgeev_    dgeev
    #define igraphdgebak_   dgebak
    #define igraphxerbla_   xerbla
    #define igraphdgebal_   dgebal
    #define igraphdisnan_   disnan
    #define igraphdlaisnan_ dlaisnan
    #define igraphdgehrd_   dgehrd
    #define igraphdgehd2_   dgehd2
    #define igraphiladlc_   iladlc
    #define igraphiladlr_   iladlr
    #define igraphdlapy2_   dlapy2
    #define igraphdlahr2_   dlahr2
    #define igraphdlacpy_   dlacpy
    #define igraphilaenv_   ilaenv
    #define igraphieeeck_   ieeeck
    #define igraphiparmq_   iparmq
    #define igraphdhseqr_   dhseqr
    #define igraphdlahqr_   dlahqr
    #define igraphdlabad_   dlabad
    #define igraphdlanv2_   dlanv2
    #define igraphdlaqr0_   dlaqr0
    #define igraphdlaqr3_   dlaqr3
    #define igraphdlaqr4_   dlaqr4
    #define igraphdlaqr2_   dlaqr2
    #define igraphdlaset_   dlaset
    #define igraphdormhr_   dormhr
    #define igraphdormqr_   dormqr
    #define igraphdtrexc_   dtrexc
    #define igraphdlaexc_   dlaexc
    #define igraphdlange_   dlange
    #define igraphdlassq_   dlassq
    #define igraphdlasy2_   dlasy2
    #define igraphdlaqr5_   dlaqr5
    #define igraphdlaqr1_   dlaqr1
    #define igraphdlascl_   dlascl
    #define igraphdorghr_   dorghr
    #define igraphdorgqr_   dorgqr
    #define igraphdorg2r_   dorg2r
    #define igraphdtrevc_   dtrevc
    #define igraphdlaln2_   dlaln2
    #define igraphdladiv_   dladiv
    #define igraphdsyevr_   dsyevr
    #define igraphdsyrk_    dsyrk
    #define igraphdlansy_   dlansy
    #define igraphdormtr_   dormtr
    #define igraphdormql_   dormql
    #define igraphdstebz_   dstebz
    #define igraphdlaebz_   dlaebz
    #define igraphdstein_   dstein
    #define igraphdlagtf_   dlagtf
    #define igraphdlagts_   dlagts
    #define igraphdlarnv_   dlarnv
    #define igraphdlaruv_   dlaruv
    #define igraphdstemr_   dstemr
    #define igraphdlae2_    dlae2
    #define igraphdlaev2_   dlaev2
    #define igraphdlanst_   dlanst
    #define igraphdlarrc_   dlarrc
    #define igraphdlarre_   dlarre
    #define igraphdlarra_   dlarra
    #define igraphdlarrb_   dlarrb
    #define igraphdlaneg_   dlaneg
    #define igraphdlarrd_   dlarrd
    #define igraphdlarrk_   dlarrk
    #define igraphdlasq2_   dlasq2
    #define igraphdlasq3_   dlasq3
    #define igraphdlasq4_   dlasq4
    #define igraphdlasq5_   dlasq5
    #define igraphdlasq6_   dlasq6
    #define igraphdlasrt_   dlasrt
    #define igraphdlarrj_   dlarrj
    #define igraphdlarrr_   dlarrr
    #define igraphdlarrv_   dlarrv
    #define igraphdlar1v_   dlar1v
    #define igraphdlarrf_   dlarrf
    #define igraphdpotrf_   dpotrf
    #define igraphdsterf_   dsterf
    #define igraphdsytrd_   dsytrd
    #define igraphdlatrd_   dlatrd
    #define igraphdsytd2_   dsytd2
    #define igraphdlanhs_   dlanhs
    #define igraphdgeqr2_   dgeqr2
    #define igraphdtrsen_   dtrsen
    #define igraphdlacn2_   dlacn2
    #define igraphdtrsyl_   dtrsyl
    #define igraphdlasr_    dlasr
    #define igraphdsteqr_   dsteqr
    #define igraphdgesv_    dgesv
    #define igraphdgetrf_   dgetrf
    #define igraphdgetf2_   dgetf2
    #define igraphdlaswp_   dlaswp
    #define igraphdgetrs_   dgetrs
    #define igraphlen_trim_ len_trim
    #define igraph_dlamc1_  dlamc1
    #define igraph_dlamc2_  dlamc2
    #define igraph_dlamc3_  dlamc3
    #define igraph_dlamc4_  dlamc4
    #define igraph_dlamc5_  dlamc5
    #define igraphddot_     ddot
#endif

int igraphdgetrf_(int *m, int *n, igraph_real_t *a, int *lda, int *ipiv,
                  int *info);
int igraphdgetrs_(char *trans, int *n, int *nrhs, igraph_real_t *a,
                  int *lda, int *ipiv, igraph_real_t *b, int *ldb,
                  int *info);
int igraphdgesv_(int *n, int *nrhs, igraph_real_t *a, int *lda,
                 int *ipiv, igraph_real_t *b, int *ldb, int *info);

igraph_real_t igraphdlapy2_(igraph_real_t *x, igraph_real_t *y);

int igraphdsyevr_(char *jobz, char *range, char *uplo, int *n,
                  igraph_real_t *a, int *lda, igraph_real_t *vl,
                  igraph_real_t *vu, int * il, int *iu,
                  igraph_real_t *abstol, int *m, igraph_real_t *w,
                  igraph_real_t *z, int *ldz, int *isuppz,
                  igraph_real_t *work, int *lwork, int *iwork,
                  int *liwork, int *info);

int igraphdgeev_(char *jobvl, char *jobvr, int *n, igraph_real_t *a,
                 int *lda, igraph_real_t *wr, igraph_real_t *wi,
                 igraph_real_t *vl, int *ldvl, igraph_real_t *vr, int *ldvr,
                 igraph_real_t *work, int *lwork, int *info);

int igraphdgeevx_(char *balanc, char *jobvl, char *jobvr, char *sense,
                  int *n, igraph_real_t *a, int *lda, igraph_real_t *wr,
                  igraph_real_t *wi, igraph_real_t *vl, int *ldvl,
                  igraph_real_t *vr, int *ldvr, int *ilo, int *ihi,
                  igraph_real_t *scale, igraph_real_t *abnrm,
                  igraph_real_t *rconde, igraph_real_t *rcondv,
                  igraph_real_t *work, int *lwork, int *iwork, int *info);

int igraphdgehrd_(int *n, int *ilo, int *ihi, igraph_real_t *A, int *lda,
                  igraph_real_t *tau, igraph_real_t *work, int *lwork,
                  int *info);

igraph_real_t igraphddot_(int *n, igraph_real_t *dx, int *incx,
                          igraph_real_t *dy, int *incy);

#endif
