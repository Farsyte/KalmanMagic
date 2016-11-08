/**  \file
 * \brief Householder Transformation Implementation
 */

#include "workers.h"
#include <math.h>

/** Triangularize (part of) a matrix using Householder Transforms
 *
 * For each of the first "Nd" diagonal elements of A, a Householder
 * transform that produces zeros in all entries below that element is
 * applied to the entire submatrix from that element to the lower right
 * corner of the matrix, inclusive.
 *
 * The Householder transform can be viewed geometrically as producing
 * the mirror image of a vector V using a plane that is defined by the
 * transforming vector U.
 *
 * At each iteration, we use the transform that examines another column
 * of the matrix, and transforms it so that all elements below the
 * diagonal are zero.
 *
 * If all elements below the diagonal are already zero, then no
 * transform is required, of couse.
 *
 * The sign of "s" is selected to avoid subtracting two large numbers
 * with a small difference, and more importantly to avoid inverting a
 * very tiny number dominated by numerical errors.
 *
 * Note that this function never explicitly stores the T matrix and the
 * U vector is overwritten by result data. Frequently the actual T
 * matrix is not needed. If the cumulative T matrix is required, append
 * an identity matrix to the right side of the input; the T matrix will
 * appear in that area.
 *
 * Runtime is O(dmn), with O(Nd) square roots and O(Nd) scalar
 * inversions.
 */

void
hh_tri (
    int Nr,
    int Nc,
    double *A,
    int Nd
)
{
    int             i,
                    k,
                    l;
    int             kNr,
                    lNr;
    double          s,
                    beta,
                    gamma;

    lNr = 0;
    for (l = 0; l < Nd; ++l)
    {

        s = 0;
        for (i = l; i < Nr; ++i)
            s += A[i + lNr] * A[i + lNr];
        s = sqrt (s);

        if (A[l + lNr] > 0)
            s = -s;

        A[l + lNr] = A[l + lNr] - s;
        beta = 1 / (s * A[l + lNr]);

        kNr = lNr + Nr;
        for (k = l + 1; k < Nc; ++k)
        {
            gamma = 0;
            for (i = l; i < Nr; ++i)
                gamma += A[i + lNr] * A[i + kNr];
            gamma *= beta;
            for (i = l; i < Nr; ++i)
                A[i + kNr] = A[i + kNr] + gamma * A[i + lNr];

            kNr += Nr;
        }

        A[l + lNr] = s;
        for (i = l + 1; i < Nr; ++i)
            A[i + lNr] = 0;

        lNr += Nr;
    }
}
