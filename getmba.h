#ifndef __GETMBA_H__
#define __GETMBA_H__

/** \file
 * \brief Provide MEX parameter access macro
 */

#define F1      " must be numeric."
#define F2      " must be double precision floating point."
#define F3      " must not be complex."
#define F4      " must not be sparse."
#define F5      " must not have so many dimensions."
#define F6      "Unable to locate Matlab storage for "
#define F7      "wrong total number of elements in "

/** Get MEX Parameter Base Address
 *
 * This macro assures that the (mxa) parameter is numeric, double
 * precision, not complex, not sparse, has (d1) rows and (d2) columns,
 * has two dimensions, and (d1*d2) total elements. It then sets (dst)
 * to the address of the first element, and checks that it is not
 * null. The mexErrMsgTxt() function is called if any tested condition
 * is not met.
 *
 * \param[out] dst = variable to hold result
 * \param[in]  mxa = Matlab Array to examine
 * \param[in]  d1  = required number of rows
 * \param[in]  d2  = required number of columns
 */
#define GETMBA(dst,mxa,d1,d2)                                          \
        do {                                                           \
            int         pd = mxGetNumberOfDimensions(mxa);             \
            long        pm = mxGetM(mxa);                              \
            long        pn = mxGetN(mxa);                              \
            long        pc = mxGetNumberOfElements(mxa);               \
            if (!mxIsNumeric(mxa))        mexErrMsgTxt(#dst F1);       \
            if (!mxIsDouble(mxa))         mexErrMsgTxt(#dst F2);       \
            if (mxIsComplex(mxa))         mexErrMsgTxt(#dst F3);       \
            if (mxIsSparse(mxa))          mexErrMsgTxt(#dst F4);       \
            if (d1 != pm) {                                            \
                if (d1 == 1)                                           \
                    mexErrMsgTxt(#dst " must have one row");           \
                else                                                   \
                    mexErrMsgTxt(#dst " must have " #d1 " rows");      \
            }                                                          \
            if (d2 != pn) {                                            \
                if (d2 == 1)                                           \
                    mexErrMsgTxt(#dst " must have one column");        \
                else                                                   \
                    mexErrMsgTxt(#dst " must have " #d2 " columns");   \
            }                                                          \
            if (2 != pd)                mexErrMsgTxt(#dst F5);         \
            if (pc != (d1 * d2))        mexErrMsgTxt(#dst F7);         \
                                                                       \
            dst = mxGetPr(mxa);                                        \
            if (!dst)                   mexErrMsgTxt(F6 #dst);         \
        } while (0)

#endif/*__GETMBA_H__*/
