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
    integer         d,
                    m,
                    n;
    doublereal     *Aa,
                   *Ab,
                   *Ac,
                   *Ad;

    setbuf (stderr, 0);

    if (nrhs != 4)
        mexErrMsgTxt ("input is (Aa,Ab,Ac,Ad)");
    if (nlhs != 4)
        mexErrMsgTxt ("output is (Aa,Ab,Ac,Ad)");

    plhs[0] = mxDuplicateArray (prhs[0]);
    plhs[1] = mxDuplicateArray (prhs[1]);
    plhs[2] = mxDuplicateArray (prhs[2]);
    plhs[3] = mxDuplicateArray (prhs[3]);

    d = mxGetM (plhs[0]);
    m = mxGetM (plhs[3]);
    n = mxGetN (plhs[3]);

    GETMBA (Aa, plhs[0], d, d);
    GETMBA (Ab, plhs[1], d, n);
    GETMBA (Ac, plhs[2], m, d);
    GETMBA (Ad, plhs[3], m, n);

    hh_tri_quad__ (&d, &m, &n, Aa, Ab, Ac, Ad);
}
