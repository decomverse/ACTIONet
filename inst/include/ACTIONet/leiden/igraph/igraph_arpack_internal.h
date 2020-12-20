/* -*- mode: C -*-  */
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

#ifndef ARPACK_INTERNAL_H
#define ARPACK_INTERNAL_H

/* Note: only files calling the arpack routines directly need to
   include this header.
*/

#include "config.h"
#include "igraph_types.h"

#ifndef INTERNAL_ARPACK
#define igraphdsaupd_ dsaupd_
#define igraphdseupd_ dseupd_
#define igraphdsaup2_ dsaup2_
#define igraphdstats_ dstats_
#define igraphdsesrt_ dsesrt_
#define igraphdsortr_ dsortr_
#define igraphdsortc_ dsortc_
#define igraphdgetv0_ dgetv0_
#define igraphdsaitr_ dsaitr_
#define igraphdsapps_ dsapps_
#define igraphdsconv_ dsconv_
#define igraphdseigt_ dseigt_
#define igraphdsgets_ dsgets_
#define igraphdstqrb_ dstqrb_
#define igraphdmout_ dmout_
#define igraphivout_ ivout_
#define igraphsecond_ second_
#define igraphdvout_ dvout_
#define igraphdnaitr_ dnaitr_
#define igraphdnapps_ dnapps_
#define igraphdnaup2_ dnaup2_
#define igraphdnaupd_ dnaupd_
#define igraphdnconv_ dnconv_
#define igraphdlabad_ dlabad_
#define igraphdlanhs_ dlanhs_
#define igraphdsortc_ dsortc_
#define igraphdneigh_ dneigh_
#define igraphdngets_ dngets_
#define igraphdstatn_ dstatn_
#define igraphdlaqrb_ dlaqrb_

#define igraphdsaupd_ dsaupd_
#define igraphdseupd_ dseupd_
#define igraphdnaupd_ dnaupd_
#define igraphdneupd_ dneupd_
#endif

#ifndef INTERNAL_LAPACK
#define igraphdlarnv_ dlarnv
#define igraphdlascl_ dlascl
#define igraphdlartg_ dlartg
#define igraphdlaset_ dlaset
#define igraphdlae2_ dlae2
#define igraphdlaev2_ dlaev2
#define igraphdlasr_ dlasr
#define igraphdlasrt_ dlasrt
#define igraphdgeqr2_ dgeqr2
#define igraphdlacpy_ dlacpy
#define igraphdsteqr_ dsteqr
#define igraphdlanst_ dlanst
#define igraphdlapy2_ dlapy2
#define igraphdlamch_ dlamch
#define igraphdlaruv_ dlaruv
#define igraphdlassq_ dlassq
#define igraphdlamc2_ dlamc2
#define igraphdlamc1_ dlamc1
#define igraphdlamc2_ dlamc2
#define igraphdlamc3_ dlamc3
#define igraphdlamc4_ dlamc4
#define igraphdlamc5_ dlamc5
#define igraphdlabad_ dlabad
#define igraphdlanhs_ dlanhs
#define igraphdtrevc_ dtrevc
#define igraphdlanv2_ dlanv2
#define igraphdlaln2_ dlaln2
#define igraphdladiv_ dladiv
#define igraphdtrsen_ dtrsen
#define igraphdlahqr_ dlahqr
#define igraphdtrsen_ dtrsen
#define igraphdlacon_ dlacon
#define igraphdtrsyl_ dtrsyl
#define igraphdtrexc_ dtrexc
#define igraphdlange_ dlange
#define igraphdlaexc_ dlaexc
#define igraphdlasy2_ dlasy2
#endif

#if 0 /* internal f2c functions always used */
#define igraphd_sign d_sign
#define igraphetime_ etime_
#define igraphpow_dd pow_dd
#define igraphpow_di pow_di
#define igraphs_cmp s_cmp
#define igraphs_copy s_copy
#define igraphd_lg10_ d_lg10_
#define igraphi_dnnt_ i_dnnt_
#endif

#ifdef HAVE_GFORTRAN

int igraphdsaupd_(int *ido, char *bmat, int *n, char *which, int *nev,
                  igraph_real_t *tol, igraph_real_t *resid, int *ncv,
                  igraph_real_t *v, int *ldv, int *iparam, int *ipntr,
                  igraph_real_t *workd, igraph_real_t *workl, int *lworkl,
                  int *info, int bmat_len, int which_len);

int igraphdseupd_(int *rvec, char *howmny, int *select, igraph_real_t *d,
                  igraph_real_t *z, int *ldz, igraph_real_t *sigma, char *bmat,
                  int *n, char *which, int *nev, igraph_real_t *tol,
                  igraph_real_t *resid, int *ncv, igraph_real_t *v, int *ldv,
                  int *iparam, int *ipntr, igraph_real_t *workd,
                  igraph_real_t *workl, int *lworkl, int *info, int howmny_len,
                  int bmat_len, int which_len);

int igraphdnaupd_(int *ido, char *bmat, int *n, char *which, int *nev,
                  igraph_real_t *tol, igraph_real_t *resid, int *ncv,
                  igraph_real_t *v, int *ldv, int *iparam, int *ipntr,
                  igraph_real_t *workd, igraph_real_t *workl, int *lworkl,
                  int *info, int bmat_len, int which_len);

int igraphdneupd_(int *rvec, char *howmny, int *select, igraph_real_t *dr,
                  igraph_real_t *di, igraph_real_t *z, int *ldz,
                  igraph_real_t *sigmar, igraph_real_t *sigmai,
                  igraph_real_t *workev, char *bmat, int *n, char *which,
                  int *nev, igraph_real_t *tol, igraph_real_t *resid, int *ncv,
                  igraph_real_t *v, int *ldv, int *iparam, int *ipntr,
                  igraph_real_t *workd, igraph_real_t *workl, int *lworkl,
                  int *info, int howmny_len, int bmat_len, int which_len);

int igraphdsortr_(char *which, int *apply, int *n, igraph_real_t *x1,
                  igraph_real_t *x2, int which_len);

int igraphdsortc_(char *which, int *apply, int *n, igraph_real_t *xreal,
                  igraph_real_t *ximag, igraph_real_t *y, int which_len);

#else

int igraphdsaupd_(int *ido, char *bmat, int *n, char *which, int *nev,
                  igraph_real_t *tol, igraph_real_t *resid, int *ncv,
                  igraph_real_t *v, int *ldv, int *iparam, int *ipntr,
                  igraph_real_t *workd, igraph_real_t *workl, int *lworkl,
                  int *info);

int igraphdseupd_(int *rvec, char *howmny, int *select, igraph_real_t *d,
                  igraph_real_t *z, int *ldz, igraph_real_t *sigma, char *bmat,
                  int *n, char *which, int *nev, igraph_real_t *tol,
                  igraph_real_t *resid, int *ncv, igraph_real_t *v, int *ldv,
                  int *iparam, int *ipntr, igraph_real_t *workd,
                  igraph_real_t *workl, int *lworkl, int *info);

int igraphdnaupd_(int *ido, char *bmat, int *n, char *which, int *nev,
                  igraph_real_t *tol, igraph_real_t *resid, int *ncv,
                  igraph_real_t *v, int *ldv, int *iparam, int *ipntr,
                  igraph_real_t *workd, igraph_real_t *workl, int *lworkl,
                  int *info);

int igraphdneupd_(int *rvec, char *howmny, int *select, igraph_real_t *dr,
                  igraph_real_t *di, igraph_real_t *z, int *ldz,
                  igraph_real_t *sigmar, igraph_real_t *sigmai,
                  igraph_real_t *workev, char *bmat, int *n, char *which,
                  int *nev, igraph_real_t *tol, igraph_real_t *resid, int *ncv,
                  igraph_real_t *v, int *ldv, int *iparam, int *ipntr,
                  igraph_real_t *workd, igraph_real_t *workl, int *lworkl,
                  int *info);

int igraphdsortr_(char *which, int *apply, int *n, igraph_real_t *x1,
                  igraph_real_t *x2);

int igraphdsortc_(char *which, int *apply, int *n, igraph_real_t *xreal,
                  igraph_real_t *ximag, igraph_real_t *y);

#endif

#endif /* ARPACK_INTERNAL_H */
