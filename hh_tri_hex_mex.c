#include "workers.h"
#include "mex.h"
#include "getmba.h"

#include <stdio.h>

/** Matlab MEX Gateway Routine
 */
void
mexFunction (
    int nlhs,
    mxArray * plhs[],
    int nrhs,
    mxArray const *prhs[]
)
{
    integer         d1,
                    d2,
                    d3;
    doublereal     *Aa,
                   *Ab,
                   *Ac,
                   *Ad,
                   *Ae,
                   *Af;

    setbuf (stderr, 0);

    if (nrhs != 6)
        mexErrMsgTxt ("input is (Aa,Ab,Ac,Ad,Ae,Af)");
    if (nlhs != 6)
        mexErrMsgTxt ("output is (Aa,Ab,Ac,Ad,Ae,Af)");

    plhs[0] = mxDuplicateArray (prhs[0]);
    plhs[1] = mxDuplicateArray (prhs[1]);
    plhs[2] = mxDuplicateArray (prhs[2]);
    plhs[3] = mxDuplicateArray (prhs[3]);
    plhs[4] = mxDuplicateArray (prhs[4]);
    plhs[5] = mxDuplicateArray (prhs[5]);

    d1 = mxGetN (plhs[0]);
    d2 = mxGetN (plhs[1]);
    d3 = mxGetN (plhs[2]);

    GETMBA (Aa, plhs[0], d1, d1);
    GETMBA (Ab, plhs[1], d1, d2);
    GETMBA (Ac, plhs[2], d1, d3);
    GETMBA (Ad, plhs[3], d2, d1);
    GETMBA (Ae, plhs[4], d2, d2);
    GETMBA (Af, plhs[5], d2, d3);

    hh_tri_hex__ (&d1, &d2, &d3, Aa, Ab, Ac, Ad, Ae, Af);
}
