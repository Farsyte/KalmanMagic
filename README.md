# KalmanMagic

## It's an old code, but it checks out.

Many years ago, I was involved with a work project that wanted to use
a Kalman Filter to estimate the system state for some simple test
vehicles, based on prior estimates, the known behavior of the systems,
and incoming noisy sensor observations. The initial implementation
using the AutoBayes package made use of a ``Bierman Update'' to
process the sensor data, which led me to a deeper study of the Kalman
filter and related methods. This project attempts to capture in
distilled form the highlights of what I have learned.

In order to organize my thinking, I wrote up my notes in TeX and set
up some prototype implementations in Matlab (we were doing model based
design using Simulink), Fortran (we still had some code in Fortran
that was slightly faster than its C equivalent), and C.

Recently, one of my good friends came to me for advice on Kalman
Filters. My notes here are all based on what was old math at the time
I wrote it -- there have been many advances in Kalman Filters beyond
what I needed to learn -- but this project should still provide a
reasonable foundation.

Additionally, where the methods described here produce acceptable
results with acceptable performance, more recent (and fancier and more
complicated and potentially more stable and potentially fster)
developments may not be needed.

In the interests of making this Readme useful, I include some of the
introductory text from the final PDF document. Markdown is fairly
limited, so refer to the actual PDF (or the TeX sources) for the
actual equations.

## Linear Discrete System Model

All filters considered in this project
operate on a linear discrete model
of the system (plus noise), as expressed
by these equations:

        z = H x + N(0,R)
        y = F x + N(0,Q)

Some of the data processing algorithms require
that the observations are normalized so that
the sensor noise on each channel is independent,
gaussian, zero mean and unit variance.
A system with a general R matrix can be normalized
by factoring R into "Sr Sr^t" then inverting Sr and
generating the normalized H matrix and $z$ vector.

        R = Sr Sr^t
        z = H x + Sr N(0,I)
        Hbar = Sr^{-1} H
        zbar = Sr^{-1} z

## Kalman Filter Basics

Standard Kalman Filters maintain an estimate x of the state of the
system, and a matrix P containing the covariance of the error of that
estimate. P=0 would represent absolute certainty that x was the exact
state of the system. When even one element of the state is completely
unknown, infinite (or IEEE NaN) values appear in P.

## Information Filter Basics

Information filters turn the Kalman filter's system state
representation on its head by storing an information matrix Im
which is the inverse of the P matrix, and storing a vector Iv
which is the product of Im and the estimated state x.
When Im is invertable,
an estimate xest of the system state
and an estimate Pest of the covaraince of the error
can be generated:

        xest = Im^{-1} Iv
        Pest = Im^{-1}

## Potter's Square Root Covariance filter

Normal Kalman filters (as described above) work directly with the
covariance of the data, and the dynamic range in magnitudes involved
can cause numerical issues. Some of these issues can be addressed by
instead working directly with a square-root form of the data. Since
the covariance matrix P is positive-definite, we can calculate another
matrix S such that S S^t = P (call this the square root matrix).
Potter's mechanization of the Kalman filter uses the S matrix directly
as part of the filter state, exhibiting improved numerical properties.

This algorithm operates on one scalar observation at a time, and
requires using the normalized Hbar matrix and zbar vectors
to eliminate the R matrix as discussed earlier. Hbarj is used
to refer to the j-th row of the normalized coefficient matrix, and
zbarj refers to the j-th entry in the normalized observation vector.

## Using the U-D Factorization

Bierman presents a significant improvement over Potter
by factoring the covariance matrix P = U D U^t
where U is a unit upper triangular matrix and D
is a diagonal matrix.
This method has benefits similar to Potter's formulation but increases
efficiency by working with triangular matrices
and reducing the number of scalar ``square root'' operations.

## Square Root Information Filters

This is a variation of Information Filters
along the same direction as
the Potter variation of the Kalman Filter.
The SRIF maintains an upper triangular matrix Im
and a vector Iv as the filter state.
The relationship between the filter state
and the system state (and our certainty)
is best captured by the ``Data Equation''

        Iv = Im x + N(0,I)

or more informally,
``Iv is Im times x,
plus or minus one-ish.''

## Appendix: Agee-Turner PD Factorzation Update

This procedure is used to update a positive definite matrix
when it is being incremented by another matrix of a
very specific form, as follows:

        Ubar Dbar Ubar^t = U D U^t + c a a^t

where:

- U and Ubar are unit upper triangular,
- D and Dbar are diagonal,
- c is a positive scalar,
- and a is a vector.

The purpose of the update is to calculate
Ubar and Dbar
when given U, D, c and a.

## Appendix: Gram-Schmidt Orthogonalization

This procedure is used to generate a collection
of vectors that are orthogonal, when weighted by
a given vector, from an input collection of
vectors that are assumed to be linearly independent.

## Appendix: Householder Triangularization

Several of the algorithms discussed in this project
require the partial triangularization of a matrix,
using orthogonal transformations;
they need the transformed matrix as a result
but do not need the transformation itself.
