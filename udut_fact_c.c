/** \file
 * \brief Implementation of U D U^t Factorization
 */

#include "workers.h"
#include <math.h>

/** Upper Triangular udut Factorizaton
 *
 * P(N,N)   matrix to be factored
 * U(N,N)    where to store U factor
 * D(N)      where to store D factor
 *
 * Given a positive definite matrix P,
 * construct the diagonal matrix D
 * and the unit upper triangular matrix U
 * such that P = U D U^t
 *
 * The content of P is destroyed.
 *
 * IMPLEMENTATION VARIATIONS:
 *
 * The D vector could be laid out along the
 * diagonal of the U matrix.
 *
 * Alternately, the U matrix can be compressed
 * into efficient vector storage, where only
 * the elements above the diagonal are stored,
 * column by column.
 */

void
udut_fact (
    int n,
    double p[],
    double u[],
    double d[]
)
{
    double          alpha;
    double          beta;
    int             i;
    int             j;
    int             k;
    int             nj;
    int             nk;

    nj = n * n;
    for (j = n - 1; j > 0; --j)
    {
        nj -= n;
        d[j] = p[j + nj];
        alpha = 1 / d[j];
        nk = 0;
        for (k = 0; k < j; ++k)
        {
            beta = p[k + nj];
            u[k + nj] = alpha * beta;
            for (i = 0; i <= k; ++i)
            {
                p[i + nk] -= beta * u[i + nj];
            }
            nk += n;
        }

        u[j + nj] = 1;
        for (k = j + 1; k < n; ++k)
        {
            u[k + nj] = 0;
        }
    }
    d[0] = p[0];
    u[0] = 1;
}
