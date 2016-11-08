/** \file
 * \brief Implementation of the U-D Update
 */
#include <alloca.h>

/** U-D Update
 *
 * \param[in]    n      = state dimension
 * \param[in]    m      = sensor dimension
 * \param[inout] x      = state estimate vector
 * \param[inout] u      = covariance, U matrix part
 * \param[inout] d      = covariance, D vector part
 * \param[in]    h      = sensor coefficient matrix
 * \param[in]    r      = sensor error variance vector
 * \param[in]    z      = sensor observation vector
 *
 * \note
 *      This procedure requires
 *      that the R matrix is diagonal.
 *
 * The covariance is P = U D U'
 * where D is diagonal, and U is unit upper triangular.
 *
 * Elements on and below the diagonal of U are neither
 * referenced nor updated.
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
udut_data (
    int n,
    int m,
    double x[],
    double u[],
    double d[],
    double h[],
    double r[],
    double z[]
)
{
    int             i,
                    j,
                    k,
                    l,
                    mk,
                    nk,
                    ml;
    double          alpha,
                    beta,
                    gamma,
                    lamda;
    double          Hjk;
    double         *b = alloca (n * sizeof *b);

    for (j = 0; j < m; ++j)
    {
        for (k = 0; k < n; ++k)
            z[j] -= h[j + m * k] * x[k];

        mk = m * n;
        nk = n * n;
        for (k = n - 1; k > 0; --k)
        {
            mk -= m;
            nk -= n;
            Hjk = h[j + mk];
            ml = 0;
            for (l = 0; l < k; ++l)
            {
                Hjk += u[l + nk] * h[j + ml];
                ml += m;
            }
            h[j + mk] = Hjk;
            b[k] = d[k] * Hjk;
        }
        b[0] = d[0] * h[j];

/*
 *     z[j] := z[j] - H'X (the computed residual),
 *     B = D U' H, and H := U'H have all been computed.
 */
        alpha = r[j] + b[0] * h[j];
        gamma = 1 / alpha;
        d[0] = r[j] * gamma * d[0];

        mk = 0;
        nk = 0;
        for (k = 1; k < n; ++k)
        {
            mk += m;
            nk += n;
            beta = alpha;
            alpha += b[k] * h[j + mk];
            lamda = -h[j + mk] * gamma;
            gamma = 1 / alpha;
            d[k] *= beta * gamma;
            for (i = 0; i < k; ++i)
            {
                beta = u[i + nk];
                u[i + nk] = beta + b[i] * lamda;
                b[i] += b[k] * beta;
            }
        }

        z[j] = z[j] * gamma;

        for (k = 0; k < n; ++k)
            x[k] += b[k] * z[j];
    }
}
