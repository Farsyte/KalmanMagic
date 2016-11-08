#!/bin/sh -

if which f2c >/dev/null 2>/dev/null
then

#
# NOTE: "f2c" will generate a number of complaints:
#
#        hh_tri:
#        hh_tri_hex:
#        hh_tri_quad:
#        potter_data:
#     Error on line 427: Declaration error for v: adjustable dimension on non-argument
#     Error on line 427: wr_ardecls:  nonconstant array size
#        srif_data:
#        srif_read:
#        srif_time:
#     Error on line 926: Declaration error for iw: adjustable dimension on non-argument
#     Error on line 926: Declaration error for tn: adjustable dimension on non-argument
#     Error on line 926: Declaration error for wx: adjustable dimension on non-argument
#     Error on line 926: wr_ardecls:  nonconstant array size
#     Error on line 926: wr_ardecls:  nonconstant array size
#     Error on line 926: wr_ardecls:  nonconstant array size
#     Error on line 926: wr_ardecls:  nonconstant array size
#     Error on line 926: wr_ardecls:  nonconstant array size
#        udut_data:
#     Error on line 1001: Declaration error for b: adjustable dimension on non-argument
#     Error on line 1001: wr_ardecls:  nonconstant array size
#        udut_fact:
#        udut_time:
#     Error on line 1208: Declaration error for gj: adjustable dimension on non-argument
#     Error on line 1208: Declaration error for vj: adjustable dimension on non-argument
#     Error on line 1208: Declaration error for tx: adjustable dimension on non-argument
#     Error on line 1208: Declaration error for dvj: adjustable dimension on non-argument
#     Error on line 1208: wr_ardecls:  nonconstant array size
#     Error on line 1208: wr_ardecls:  nonconstant array size
#     Error on line 1208: wr_ardecls:  nonconstant array size
#     Error on line 1208: wr_ardecls:  nonconstant array size
#
# We can safely ignore these, as we are only going to be
# looking at the prototypes for the functions.
#
# If our target platform demanded strict Fortran 77 then
# these variable sized arrays, mostly working storage, would
# have to be provided by the caller.
#
# And we would get complaints when we *COMPILED* them, so we
# would not be looking at the fprotos.sh error list anyway.
#

    cat fprotos.h.head
    cat "$@" < /dev/null | f2c -P -!c 2>/dev/null
    cat fprotos.h.tail

elif [ -f fprotos.h.ref ]; then

    cat fprotos.h.ref

else

    cat fprotos.h.head
    printf '#warning "Prototypes for Fortran implementations are not available: f2c not installed and fprotos.h.ref is missing."\n'
    cat fprotos.h.tail

fi
