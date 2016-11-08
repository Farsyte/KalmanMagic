#ifndef __CFORT_H__
#define __CFORT_H__

/** \file
 * \brief Data types to use in C code that calls Fortran code
 */

/** Integer Data Type
 * In this package, C code calling Fortran code
 * should use the "integer" data type to match
 * Fortran INTEGER data.
 */
typedef long    integer;

/** Floating Point Data Type
 * In this package, C code calling Fortran code
 * should use the "doublereal" data type to match
 * Fortran DOUBLE PRECISION data.
 */
typedef double  doublereal;

#endif/*__CFORT_H__*/
