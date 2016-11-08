/** \file
 * \brief Square Root Information Filter Estimate Production
 */

#include "workers.h"
#include <math.h>

/** Square Root Information Filter Estimate Production
 *
 * This function examines the state of the SRIF and
 * provides the current state estimate and covariance;
 * as a byproduct, also provides S (where S S' = P).
 *
 * The system is described by
 * .     Iv = Im x + v
 *
 * \note Im is upper triangular,
 * which simplifies this tremendously.
 *
 * ACTUAL DIMENSIONS: X(n),P(n,n),S(n,n),IM(n,n),IV(n)
 */
void
srif_read (
    int n,
    double x[],
    double p[],
    double s[],
    double im[],
    double iv[]
)
{
    int             i,
                    j,
                    k;
    int             ni,
                    nj,
                    nk;
    double          sum;

/*
 *     X = Im \ IV
 */
    nj = n * n;
    for (j = n - 1; j >= 0; --j)
    {
        nj -= n;
        sum = 0;
        nk = nj;
        for (k = j + 1; k < n; ++k)
        {
            nk += n;
            sum += im[j + nk] * x[k];
        }
        x[j] = (iv[j] - sum) / im[j + nj];
    }
/*
 *     S = inv(Im)
 */
    s[0] = 1 / im[0];
    nj = 0;
    for (j = 1; j < n; ++j)
    {
        nj += n;
        s[j + nj] = 1 / im[j + nj];
        nk = 0;
        for (k = 0; k < j; ++k)
        {
            sum = 0;
            ni = nk;
            for (i = k; i < j; ++i)
            {
                sum -= s[k + ni] * im[i + nj];
                ni += n;
            }
            s[k + nj] = sum * s[j + nj];
            nk += n;
        }
    }
/*
 *     P = S * S'
 */
    nj = 0;
    for (j = 0; j < n; ++j)
    {
        ni = 0;
        for (i = 0; i < j; ++i)
        {
            sum = 0;
            nk = nj;
            for (k = j; k < n; ++k)
            {
                sum += s[i + nk] * s[j + nk];
                nk += n;
            }
            p[i + nj] = p[j + ni] = sum;
            ni += n;
        }
        sum = 0;
        nk = nj;
        for (k = j; k < n; ++k)
        {
            sum += s[j + nk] * s[j + nk];
            nk += n;
        }
        p[j + nj] = sum;
        nj += n;
    }
}
