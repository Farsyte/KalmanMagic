/** \file
 * \brief Square Root Information Filtering: Time Processing
 */
#include "workers.h"
#include <alloca.h>
#include <math.h>

/** Square Root Information Filtering: Time Processing
 *
 * We are going to use a trick here, similar to what we just did
 * in the SRIF Data processing above: transform the whole problem
 * into a single matrix, and triangularize the matrix to get the
 * results that we are seeking.
 *
 * The equation of interest is as follows:
 *
 * .    new X = F X + G W
 *
 * where
 * .    F   is the state update transform,
 * .    G   is the process noise mapping matrix,
 * .    W   is the process noise vector from N(0,I)
 *
 * We are going to use G to represent our process noise sources:
 * it is merely the mapping of a collection of independent
 * gaussian noise sources with zero mean and unit variance
 * (represented by the W vector in the equation) into the state
 * space. The covariance matrix for the process noise Q is G*G'.
 *
 * We could factor Q into GG' but in practice, it is as easy to
 * initially construct G for the system as it is to construct Q
 * (in fact, when I look back at my constructions, I tended to
 * construct Q by first building a matrix that could have been
 * calleed G).
 *
 * We have the following equations coming into the procedure:
 *
 * .    Iv0 = Im0 X0 +   Vx,   Vx from N(0,I)
 * .     X1 =  F  X0 + G W0
 * and describing "W0" we have:
 * .    Wv0 = Wm0 W0 +   Vw,   Vw from N(0,I)
 *
 * ( Note that this actualy gives us a lot more freedom in
 * . what is going on with W0 than we need; in practice,
 * . we can select G so that W0 is from N(0,I) but that
 * . does not actually simplify the math that is coming. )
 *
 * Similar to before, we want an (Iv1,Im1) that will
 * determine the unknown X1 state, as always, minimizing
 * the error in a least-squares sense.
 *
 * First, express it all in terms of X1. This will of course
 * require us to have on hand Fi, the inverse of F ...
 *
 * Solving the second equation for X0 in terms of X1,
 * .     X0 = Fi (X1 - G w)
 * then substituting that back into the first,
 * .    Iv0 = Im0 Fi (X1 - G W0) + v
 *
 * Combining the equation describing W0 with
 * the new equation describing X1 and W0,
 *
 * .    Wv0 =                    Wm0   W0 + Vw
 * .    Iv0 = (Im0 Fi) X1 - (Im0 Fi G) W0 + Vx
 *
 * rewriting this back into our matrix forms and turning it
 * into a performance functional,
 *
 * .          ||  |  Wm0         0    |   | W0 |    | Wv0 | ||2
 * J(W0,X1) = ||  |                   | * |    | -  |     | ||
 * .          ||  | -Im0*Fi*G  Im0*Fi |   | X1 |    | Iv0 | ||
 *
 * Now hit the stuff inside the squared-two-norm with a nice
 * Householder diagonalization transform ... and we get:
 *
 * .          ||  |  Wm0         0    |   | W0 |    | Wv0 | ||2
 * J(W0,X1) = || T|                   | * |    | - T|     | ||
 * .          ||  | -Im0*Fi*G  Im0*Fi |   | X1 |    | Iv0 | ||
 *
 * .          ||  |  ^T1        ^T2   |   | W0 |    | ^T3 | ||2
 * .        = ||  |                   | * |    | -  |     | ||
 * .          ||  |   0         Im1   |   | X1 |    | Iv1 | ||
 *
 * which separates into minimizing
 * .    ^T1 W0 + ^T2 X1 - ^T3
 * and
 * .             Im1 X1 - Iv1
 *
 * Note that this last line is precisely the form we want. All we
 * lack is the proof that J(W0,X1) is minimized by X1 when the
 * last line is minimized; Bierman goes through the magic that
 * proves this result. Essentially, ^T1 is invertable diagonal, so
 * the top part can always minimize to zero when W0 = (^T2 X1 -
 * ^T3) / ^T1 which means that the value of X1 is determined only
 * by the last line.
 *
 * So once again, we marshal our data and apply Householder to the
 * result, giving us a nice upper triangular Im1 matrix and an Iv1
 * vector that we can use as the filter state going forward.
 *
 * REAL DIMENSIONS: IM(N,N),IV(N),FI(N,N),G(N,M),WM(N,N),WV(M)
 */
void
srif_time (
    int n,
    int m,
    double im[],
    double iv[],
    double fi[],
    double g[],
    double wm[],
    double wv[]
)
{
    int             i,
                    j,
                    k,
                    l;
    double          sum;
    double          oldd,
                    oldd2,
                    s2,
                    s,
                    topu;
    double          beta,
                    gamma;
    double         *wx = alloca (m * n * sizeof *wx);
    double         *iw = alloca (n * m * sizeof *iw);
    double         *tn = alloca (n * sizeof *tn);

/*
 *     Im = Im*Fi
 */
    for (i = 0; i < n; ++i)
    {
        for (j = 0; j < n; ++j)
        {
            sum = 0;
            for (k = 0; k < n; ++k)
            {
                sum = sum + im[i + n * k] * fi[k + n * j];
            }
            tn[j] = sum;
        }
        for (j = 0; j < n; ++j)
        {
            im[i + n * j] = tn[j];
        }
    }

/*
 *     iw = -Im*G
 */
    for (j = 0; j < m; ++j)
    {
        for (i = 0; i < n; ++i)
        {
            sum = 0;
            for (k = 0; k < n; ++k)
            {
                sum = sum + im[i + n * k] * g[k + n * j];
            }
            iw[i + n * j] = -sum;
        }
    }
/*
 *     wx = zeros(m,n)
 */
    for (j = 0; j < n; ++j)
    {
        for (i = 0; i < m; ++i)
        {
            wx[i + m * j] = 0;
        }
    }
/*
 *     Matrix to be triangularized is
 *
 *     .     | Wm  wx  Wv |  (m)
 *     .     |            |
 *     .     | iw  Im  Iv |  (n)
 *     .       (m) (n) (1)
 *
 *     "wx" is zero when we start.
 *     "iw" is zero when we finish.
 *
 *     This function partially triangularizes this matrix by walking
 *     down the diagonal of Wm then Im, selecting a householder
 *     transformation that clears elements below that diagonal
 *     element, and applying that transformation to the remainder of
 *     the matrix. This will result in Wm triangular, iw zero, and Im
 *     triangular.
 *
 *     The Householder transform can be viewed geometrically as
 *     producing the mirror image of a vector Y using a plane that is
 *     defined by the transforming vector U.
 *
 *     At each iteration, we use the transform that takes the column
 *     vector from the selected element to the bottom of the matrix,
 *     and transforms it so that all elements but the first are zero.
 *
 *     The sign of "s" is selected to avoid subtracting two large
 *     numbers with a small difference, and more importantly to avoid
 *     inverting a very tiny number dominated by numerical errors.
 *
 *     Runtime is O(m*(m+n+1)*(m+n))
 */
    for (l = 0; l < m; ++l)
    {

        oldd = wm[l + n * l];
        oldd2 = oldd * oldd;

        s2 = oldd2;
        for (i = l + 1; i < m; ++i)
        {
            s2 = s2 + wm[i + n * l] * wm[i + n * l];
        }
        for (i = 0; i < n; ++i)
        {
            s2 = s2 + iw[i + n * l] * iw[i + n * l];
        }

        if (s2 > oldd2)
        {

            /*
             * WARNING: SQRT IS EXPENSIVE
             */
            s = sqrt (s2);

            if (oldd > 0)
                s = -s;

            topu = oldd - s;
            beta = 1 / (s * topu);

            for (j = l + 1; j < m; ++j)
            {

                gamma = topu * wm[l + n * j];
                for (i = l + 1; i < m; ++i)
                {
                    gamma = gamma + wm[i + n * l] * wm[i + n * j];
                }
                for (i = 0; i < n; ++i)
                {
                    gamma = gamma + iw[i + n * l] * iw[i + n * j];
                }
                gamma = gamma * beta;

                wm[l + n * j] = wm[l + n * j] + gamma * topu;
                for (i = l + 1; i < m; ++i)
                {
                    wm[i + n * j] =
                        wm[i + n * j] + gamma * wm[i + n * l];
                }
                for (i = 0; i < n; ++i)
                {
                    iw[i + n * j] =
                        iw[i + n * j] + gamma * iw[i + n * l];
                }
            }

            for (j = 0; j < n; ++j)
            {
                gamma = topu * wx[l + m * j];
                for (i = l + 1; i < m; ++i)
                {
                    gamma = gamma + wm[i + n * l] * wx[i + m * j];
                }
                for (i = 0; i < n; ++i)
                {
                    gamma = gamma + iw[i + n * l] * im[i + n * j];
                }
                gamma = gamma * beta;

                wx[l + m * j] = wx[l + m * j] + gamma * topu;
                for (i = l + 1; i < m; ++i)
                {
                    wx[i + m * j] =
                        wx[i + m * j] + gamma * wm[i + n * l];
                }
                for (i = 0; i < n; ++i)
                {
                    im[i + n * j] =
                        im[i + n * j] + gamma * iw[i + n * l];
                }
            }

            gamma = topu * wv[l];
            for (i = l + 1; i < m; ++i)
            {
                gamma = gamma + wm[i + n * l] * wv[i];
            }
            for (i = 0; i < n; ++i)
            {
                gamma = gamma + iw[i + n * l] * iv[i];
            }
            gamma = gamma * beta;

            wv[l] = wv[l] + gamma * topu;
            for (i = l + 1; i < m; ++i)
            {
                wv[i] = wv[i] + gamma * wm[i + n * l];
            }
            for (i = 0; i < n; ++i)
            {
                iv[i] = iv[i] + gamma * iw[i + n * l];
            }

            wm[l + n * l] = s;
        }
        for (i = l + 1; i < m; ++i)
        {
            wm[i + n * l] = 0;
        }
        for (i = 0; i < n; ++i)
        {
            iw[i + n * l] = 0;
        }
    }

    for (l = 0; l < n - 1; ++l)
    {

        oldd = im[l + n * l];
        oldd2 = oldd * oldd;
        s2 = oldd2;
        for (i = l + 1; i < m; ++i)
        {
            s2 = s2 + im[i + n * l] * im[i + n * l];
        }
        if (s2 > oldd2)
        {

            /*
             * WARNING: SQRT IS EXPENSIVE
             */
            s = sqrt (s2);

            if (oldd > 0)
                s = -s;

            topu = oldd - s;
            beta = 1 / (s * topu);

            for (j = l + 1; j < n; ++j)
            {
                gamma = topu * im[l + n * j];
                for (i = l + 1; i < m; ++i)
                {
                    gamma = gamma + im[i + n * l] * im[i + n * j];
                }
                gamma = gamma * beta;

                im[l + n * j] = im[l + n * j] + gamma * topu;
                for (i = l + 1; i < m; ++i)
                {
                    im[i + n * j] =
                        im[i + n * j] + gamma * im[i + n * l];
                }
            }

            gamma = topu * iv[l];
            for (i = l + 1; i < m; ++i)
            {
                gamma = gamma + im[i + n * l] * iv[i];
            }
            gamma = gamma * beta;

            iv[l] = iv[l] + gamma * topu;
            for (i = l + 1; i < m; ++i)
            {
                iv[i] = iv[i] + gamma * im[i + n * l];
            }

            im[l + n * l] = s;
        }
        for (i = l + 1; i < m; ++i)
        {
            im[i + n * l] = 0;
        }
    }

}
