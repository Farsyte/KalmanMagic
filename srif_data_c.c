/** \file
 * \brief SRIF Data Processing using Householder Tranfsormations
 */

#include "workers.h"
#include <math.h>

/** SRIF Data Processing using Householder Tranfsormations
 *
 * Observations are merged into the SRIF state by rewriting
 * the problem in an interesting form, then using Householder
 * transformations on a matrix, and drawing the results out
 * of the transformed matrix.
 *
 * We arrange our equations in this form:
 *
 *     .      ||  |  Im |       |  Iv | ||2
 *     J(x) = ||  |     | X0 -  |     | ||
 *     .      ||  |  H  |       |  z  | ||
 *
 *
 * where the state of the sytem "x" minimizes J. We want
 * of course to transform this to give us a new Im,Iv to
 * replace the old one, which includes the H,z information.
 *
 * Applying an orthogonal transformation T to the vector being
 * minimized will not change its length, and J will be minimized by
 * the same X vector:
 *
 *     .      ||  |  Im |       |  Iv | ||2
 *     J(x) = || T|     | X0 - T|     | ||
 *     .      ||  |  H  |       |  z  | ||
 *
 * We carefully select the T matrix to triangularize the left hand
 * block, leaving us with an equation in this form:
 *
 *     .      ||  |  Im |       |  Iv | ||2
 *     J(x) = ||  |     | X0 -  |     | ||
 *     .      ||  |  0  |       |  e  | ||
 *
 * Apologies for re-using Im and Iv here, expanding the names to
 * differentiate between prior estimates and updated estimates makes
 * it pretty unreadable. Note that the new Im is upper triangular;
 * the new Im and Iv are appropriate for use going forward to
 * represent a refined estimate of the state X, and the E value left
 * in the lower part of the right box represents a cumulative amount
 * of known error that can not be reduced out of J by any X.
 *
 * Note that we never actually construct T. Our only interest in it
 * is the assertion that it is an orthogonal transform!
 *
 * This function operates "in place" on all data items, updating Im
 * and Iv in place, and destroying H and Z.
 *
 * @param n     dimension of state vector
 * @param m     dimension of sense vector
 * @param im    estimate information matrix
 * @param iv    estimate information vector
 * @param h     sensor coefficients
 * @param z     sensor readings
 *
 * Runtime is O(N * N * (N+M))
 *
 * ACTUAL DIMENSIONS: IM(n,n),IV(n),H(m,n),Z(m)
 */

void
srif_data (
    int n,
    int m,
    double im[],
    double iv[],
    double h[],
    double z[]
)
{
    int             d,
                    r,
                    c;
    int             md,
                    nd;
    int             mc,
                    nc;
    double          s,
                    beta,
                    gamma,
                    imd,
                    imd2,
                    utop;

    md = 0;
    nd = 0;

    for (d = 0; d < n; ++d)
    {

        imd = im[d + nd];
        imd2 = imd * imd;

        s = imd2;
        for (r = d + 1; r < n; ++r)
            s += im[r + nd] * im[r + nd];
        for (r = 0; r < m; ++r)
            s += h[r + md] * h[r + md];

        if (s > imd2)
        {

            /*
             * WARNING: SQRT IS EXPENSIVE
             */

            s = sqrt (s);
            if (imd > 0)
                s = -s;

            utop = imd - s;
            beta = 1 / (s * utop);

            im[d + nd] = s;

            mc = md;
            nc = nd;
            for (c = d + 1; c < n; ++c)
            {
                mc += m;
                nc += n;

                gamma = utop * im[d + nc];
                for (r = d + 1; r < n; ++r)
                    gamma += im[r + nd] * im[r + nc];
                for (r = 0; r < m; ++r)
                    gamma += h[r + md] * h[r + mc];
                gamma *= beta;

                im[d + nc] += gamma * utop;
                for (r = d + 1; r < n; ++r)
                    im[r + nc] += gamma * im[r + nd];
                for (r = 0; r < m; ++r)
                    h[r + mc] += gamma * h[r + md];
            }

            gamma = utop * iv[d];
            for (r = d + 1; r < n; ++r)
                gamma += im[r + nd] * iv[r];
            for (r = 0; r < m; ++r)
                gamma += h[r + md] * z[r];
            gamma = beta * gamma;

            iv[d] += gamma * utop;
            for (r = d + 1; r < n; ++r)
                iv[r] += gamma * im[r + nd];
            for (r = 0; r < m; ++r)
                z[r] += gamma * h[r + md];
        }

        for (r = d + 1; r < n; ++r)
            im[r + nd] = 0;
        for (r = 0; r < m; ++r)
            h[r + md] = 0;

        md += m;
        nd += n;
    }
}
