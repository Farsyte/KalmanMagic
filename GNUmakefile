# Note: this build process generates prototypes for
# the Fortran code using F2C which means that it is
# compatible with g77 but not gfortran.
#
# Shifting from g77 to gfortran will require some
# changes to how the function prototypes are built.
#
# Also ...
#
# MEX creates, uses, and destroys an object file
# called "mexversion.o" -- so attempting to run more
# than one MEX command at the same time in the same
# directory can (and will) fail.
#
# Arrange for MEX to use different subdirectories
# for each execution.
#

P               := KalmanMagic

Q		= @		# GMake prefix: do not echo this command.
I		= -		# GMake prefix: continue even if this fails.
C		= $I $Q		# Use Both Prefixes for "just do this" commands.
N		= $C true			# Common "do nothing" command
R		= $C /bin/rm -vfr		# Common "remove" command

N		= -@true	# silently do nothing.

CC		= gcc		# select a C compiler (gcc recommended)
CSTD		= --std=c99	# select --std=c89 or --std=c99

FC		= gfortran	# select a Fortran compiler (g77 or gfortran).
FSTD		= --std=f95	# select -fno-f77 or --std=f95

# USE_F		= -DUSE_FORTRAN
USE_F		= -UUSE_FORTRAN

WARN		= -W -Wall -Wextra -pedantic
FNAMES		= -fsecond-underscore
OPT		= -O
DBG		= -g
CFLAGS		= ${OPT} ${DBG} ${CSTD} ${WARN} -fPIC ${USE_F}
FFLAGS		= ${OPT} ${DBG} ${FSTD} ${WARN} -fPIC ${FNAMES}

MEX		= mex		# Matlab MEX compilation
MEXFLAGS	=		# Matlab MEX flags

# Specify final result file names.

WORKLIB		= libwork.a

# Locate source files.

ALLTEX		= ${wildcard *.tex}
HEADERS		= ${wildcard *.h}
CWORKERS	= ${wildcard *_c.c}
FWORKERS	= ${wildcard *_f.f}
MEXC		= ${wildcard *_mex.c}
# MEXC		= ${shell grep -l mexFunction ${CMEXFUNC}}
HCF		= $(HEADERS) $(MEXC) $(CWORKERS) $(FWORKERS)

# Lists of generated files.

WORKERS		= ${CWORKERS:%.c=%.o} \
		  ${FWORKERS:%.f=%.o}
MEXDIR		= ${MEXC:%.c=%.d}
MEXOUT		= ${MEXC:%.c=%.${MEXEXT}}

# Select proper extention for MEX binaries,
# which changes from platform to platform.

MEXEXT		= ${shell mexext}

# Stock Targets

default:			; $N
prep:				; $N
build:				; $N
post:				; $N
publish:			; $N
clean::				; $N
cleaner::			; $N

rebuild_world:
	$(MAKE) clean
	$(MAKE) cleaner
	$(MAKE) prep
	$(MAKE) build
	$(MAKE) post
	$(MAKE) publish


# Stock Target Heirarchy

build:		prep
post:		build
default:	post
cleaner::	clean

# This project requires a header file containing
# C prototypes for the Fortran functions. Build it.
#
# Assumes we have "f2c" which may no longer exist.

fprotos.h:			; /bin/sh fprotos.sh ${FWORKERS} | indent -st > $@
prep:		fprotos.h
cleaner::			; $R fprotos.h

# Rules to compile the Workers.

${WORKERS}:	${HEADERS} ${FWPROTOS} fprotos.h
prep:		${WORKERS}
clean::				; $R ${WORKERS}

# Combine the workers into a library.

${WORKLIB}:	${WORKERS}	; ${AR} ${ARFLAGS} $@ ${WORKERS}
prep:		${WORKLIB}
clean::				; $R ${WORKLIB}

# Execute MEX processing (in a subdirectory).

${MEXDIR}: 	${UTIL} ${WORKLIB}
%.d: %.c
	-@mkdir -p $@
	${MEX} ${MEXFLAGS} $< ${WORKLIB} -outdir $@
build:		${MEXDIR}
clean::				; $R ${MEXDIR}

# Retrieve the MEX object from the subdirectory.

%.${MEXEXT}: %.d
	cp $</$@ $@
build:		${MEXOUT}
clean::				; $R ${MEXOUT}


# Various source code statistics.
# Too bad PMCCABE does not grok Matlab code.

stats:		$(HCF)		; ./stats.sh $(HCF) > $@
post:		stats
cleaner::			; $R stats

# Use Doxygen to extract some documentation.
doxy:           Doxyfile        ; doxygen
cleaner::                       ; $R html latex


# Publish the TEX document as a DVI file.
# XXX: modern latex installations go directly to PDF.
# $P.dvi:         ${ALLTEX}	; ./tex_dvi.sh $P $P.tex ${filter-out $P.tex, ${ALLTEX}}
# dvi:				$P.dvi
# cleaner::			; $R $P.dvi
# cleaner::			; $R $P.dvi.d

# Publish the TEX document as a PDF file.

$P.pdf:         ${ALLTEX}	; ./tex_pdf.sh $P $P.tex ${filter-out $P.tex, ${ALLTEX}}
pdf:				$P.pdf
cleaner::			; $R $P.pdf
cleaner::			; $R $P.pdf.d

PFILES		= ${shell cat $P.list}
$P.tbz:	$P.list ${PFILES}		; tar cfj $@ ${PFILES}

publish:	$P.tbz

cleaner::					; $R $P.tbz
cleaner::					; $R $P.aux
cleaner::					; $R $P.log
cleaner::					; $R $P.toc
cleaner::					; $R $P.db
