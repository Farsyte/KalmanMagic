#ifndef __WORKERS_H__
#define __WORKERS_H__

/** \file
 * \brief Prototypes for worker functions.
 *
 * This file provides C Function Declarations
 * for the low level worker functions.
 */

#include "cfort.h"
#include "fprotos.h"

extern void     hh_tri (
    int m,
    int n,
    double *a,
    int d
);

extern void     potter_data (
    int n,
    int m,
    double x[],
    double s[],
    double h[],
    double z[]
);

extern void     udut_fact (
    int n,
    double p[],
    double u[],
    double d[]
);

extern void     udut_data (
    int n,
    int m,
    double x[],
    double u[],
    double d[],
    double h[],
    double r[],
    double z[]
);

extern void     udut_time (
    int n,
    int m,
    double *x,
    double *u,
    double *d,
    double *f,
    double *g
);

extern void     srif_data (
    int n,
    int m,
    double im[],
    double iv[],
    double h[],
    double z[]
);

extern void     srif_read (
    int n,
    double x[],
    double p[],
    double s[],
    double im[],
    double iv[]
);

extern void     srif_time (
    int n,
    int m,
    double im[],
    double iv[],
    double fi[],
    double g[],
    double wm[],
    double wv[]
);

#endif/*__WORKERS_H__*/
