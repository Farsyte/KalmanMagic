/** \file
 * \brief Time update of U-D factors, using MWG-S algorithm.
 */

#include "workers.h"
#include <alloca.h>

/** Time update of U-D factors, using MWG-S algorithm.
 *
 * \param[in]    n      = state dimension
 * \param[in]    m      = sensor dimension
 * \param[inout] x      = state estimate vector
 * \param[inout] u      = covariance, U matrix part
 * \param[inout] d      = covariance, D vector part
 * \param[in]    f      = system update model
 * \param[in]    g      = 
 *
 * IMPLEMENTATION VARIATIONS:
 *      The D vector could be laid out along the
 *      diagonal of the U matrix.
 *
 *      Alternately, the U matrix can be compressed
 *      into efficient vector storage, where only
 *      the elements above the diagonal are stored,
 *      column by column.
 */
void
udut_time (
    int n,
    int m,
    double x[],
    double u[],
    double d[],
    double f[],
    double g[]
)
{
    int             i;
    int             j;
    int             k;
    int             ni;
    int             nj;
    int             nk;
    int             nn;
    double         *gj = alloca (m * sizeof *gj);
    double         *vj = alloca (n * sizeof *vj);
    double         *tx = alloca (n * sizeof *tx);
    double         *dv = alloca (n * sizeof *dv);
    double          sum;
    double          dinv;
    double          sigma;

    nn = n * n;

    /*
     * Calculate U := F*U
     *     using TX to hold U(:,J)
     *     while it is being rewritten.
     */
    nj = nn;
    for (j = n - 1; j >= 0; --j)
    {
        nj -= n;
        for (i = 0; i < j; ++i)
        {
            tx[i] = u[i + nj];
        }
        for (i = 0; i < n; ++i)
        {
            sum = f[i + nj];
            nk = 0;
            for (k = 0; k < j; ++k)
            {
                sum += f[i + nk] * tx[k];
                nk += n;
            }
            u[i + nj] = sum;
        }
    }

    /*
     * Calculate X := F*X
     *     using TX to hold the old X.
     */

    for (j = 0; j < n; ++j)
    {
        tx[j] = x[j];
    }

    for (j = 0; j < n; ++j)
    {
        sum = 0;
        nk = 0;
        for (k = 0; k < n; ++k)
        {
            sum += f[j + nk] * tx[k];
            nk += n;
        }
        x[j] = sum;
    }

    /*
     * Modified Weighted Gram-Schmidt Orthogonalization
     */

    nj = nn;
    for (j = n - 1; j > 0; --j)
    {
        nj -= n;

        /*
         * Ujj := Vj D Vj + Gj Q Gj
         */

        sum = 0;
        nk = 0;
        for (k = 0; k < n; ++k)
        {
            vj[k] = u[j + nk];
            dv[k] = d[k] * vj[k];
            sum += vj[k] * dv[k];
            nk += n;
        }
        nk = 0;
        for (k = 0; k < m; ++k)
        {
            gj[k] = g[j + nk];
            sum += gj[k] * gj[k];
            nk += n;
        }
        u[j + nj] = sum;
        dinv = 1 / sum;

        /*
         * DINV is 1/Dj
         */

        nk = 0;
        for (k = 0; k < j; ++k)
        {

            /*
             * Ujk := ( Uk D Vj + Gk Q Gj ) / Ujj
             */

            sum = 0;
            ni = 0;
            for (i = 0; i < n; ++i)
            {
                sum += u[k + ni] * dv[i];
                ni += n;
            }
            ni = 0;
            for (i = 0; i < m; ++i)
            {
                sum += g[k + ni] * gj[i];
                ni += n;
            }
            sigma = sum * dinv;
            u[j + nk] = sigma;

            /*
             * Uk := Uk - Vj Ujk
             */

            ni = 0;
            for (i = 0; i < n; ++i)
            {
                u[k + ni] -= sigma * vj[i];
                ni += n;
            }

            /*
             * Gk := Gk - Gj Ujk
             *
             *     This destroys the G parameter, but skipping
             *     this step introduces errors.
             */

            ni = 0;
            for (i = 0; i < m; ++i)
            {
                g[k + ni] -= sigma * gj[i];
                ni += n;
            }
            nk += n;
        }
    }

    /*
     * Ujj := Vj D Vj + Gj Q Gj   (specific to j=0)
     */

    sum = 0;
    nk = 0;
    for (k = 0; k < n; ++k)
    {
        sum += u[nk] * u[nk] * d[k];
        nk += n;
    }
    nk = 0;
    for (k = 0; k < m; ++k)
    {
        sum += g[nk] * g[nk];
        nk += n;
    }

    /*
     * WAIT! Instead of stuffing SUM back into U(1,1) ...
     *
     *     Move all results into their final locations:
     *     D data is found along the diagonal of the U storage, and
     *     U data is found in the lower triangular part.
     *     Need to move D data out, copy U data into place,
     *     then set U diagonal elements to 1 and elements below
     *     the diagonal to 0.
     *
     *     NOTE: be careful to copy the resuts out of the lower
     *     triangular area BEFORE clearing that area. Do not fall into
     *     the trap of reordering the loops to optimally pass over the
     *     U memory.
     */
    d[0] = sum;
    u[0] = 1;
    nj = 0;
    for (j = 1; j < n; ++j)
    {
        nj += n;
        d[j] = u[j + nj];
        ni = 0;
        for (i = 0; i < j; ++i)
        {
            u[i + nj] = u[j + ni];
            u[j + ni] = 0;
            ni += n;
        }
        u[j + nj] = 1;
    }
}
