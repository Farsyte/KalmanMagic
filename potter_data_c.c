/**  \file
 * \brief Data Processing using the Potter Mechanization
 */

#include <alloca.h>
#include <math.h>

/** Data Processing using the Potter Mechanization
 *
 * This function updates the estimate of the state x and the state
 * error covariance square root estimate S based on an observation
 * vector z derived from the state via the observation coefficient
 * matrix H plus a noise vector.
 *
 * The observation coefficient vector and observation value must be
 * normalized so that the variance of the observation noise is one.
 *
 * \param n state vector dimension
 * \param m sensor vector dimension
 * \param x (N   input)  state estimate (in-place update)
 * \param s (NxN input)  error covariance (in-place update)
 * \param h (MxN input)  normalized observation coefficients
 * \param z (M   input)  normalized noisy observation
 */
void
potter_data (
    int n,
    int m,
    double x[],
    double s[],
    double h[],
    double z[]
)
{
    int             i,
                    j,
                    k;
    double          t,
                    d,
                    a,
                    vi;
    double         *v = alloca (n * sizeof *v);

    for (j = 0; j < m; ++j)
    {

        t = 1;
        d = z[j];

        for (i = 0; i < n; ++i)
        {
            vi = 0;
            for (k = 0; k < n; ++k)
                vi += s[k + n * i] * h[j + m * k];
            v[i] = vi;
            d -= h[j + m * i] * x[i];
            t += vi * vi;
        }

        t = 1 / t;
        d *= t;

        t /= (1 + sqrt (t));

        for (i = 0; i < n; ++i)
        {
            a = 0;
            for (k = 0; k < n; ++k)
                a += s[i + n * k] * v[k];
            x[i] += a * d;
            a *= t;
            for (k = 0; k < n; ++k)
                s[i + n * k] -= a * v[k];
        }
    }
}
