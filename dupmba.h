#ifndef __DUPMBA_H__
#define __DUPMBA_H__

/** \file
 * \brief Provide MEX Array Duplication macro for mexFunction
 */

#include <alloca.h>

/** Duplicate a MEX array
 *
 * This is useful when the implementation wants a parameter
 * to contain the input data, but will destroy it. We do not
 * want to scribble on the array Matlab hands us as an input.
 *
 * This is SPECIFICALLY ONLY to be used from the mexFunction
 * itself, as failures are reported via mexErrMsgTxt calls
 * which do not operate properly from anywhere else.
 *
 * The storage allocated goes away when mexFunction returns.
 *
 * @param dst where to store new pointer
 * @param mxa matlab array to duplicate
 * @param d1 first dimension of array
 * @param d2 second dimension of array
 */
#define DUPMBA(dst,mxa,d1,d2)                                   \
    do {                                                        \
        long eix, noe = d1 * d2;                                \
        double *src;                                            \
        GETMBA(src,mxa,d1,d2);                                  \
        dst = alloca(sizeof *dst * noe);                        \
        if (!dst)                                               \
            mexErrMsgTxt(                                       \
                "unable to allocate output space for " #dst);   \
        for (eix = 0; eix < noe; ++eix)                         \
            dst[eix] = src[eix];                                \
    } while (0)

#endif/*__DUPMBA_H__*/
