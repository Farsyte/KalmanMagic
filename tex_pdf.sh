#!/bin/sh -
#
# Translate a TeX document to PDF format.
#
# usage:
#	tex_pdf Name file1.tex file2.tex ...

set -e

base="$1"
outf="${base}.pdf"
outd="${outf}.d"
shift

mkdir -p "${outd}"
ln -f "$@" "${outd}"
cd "${outd}"

    echo "Running pdflatex...."
pdflatex "$@"			>/dev/null

if [ -f ${base}.idx ]
then
    echo "Running makeindex...."
    makeindex ${base}.idx	>/dev/null

    echo "Rerunning pdflatex...."
    pdflatex "$@"			>/dev/null
fi

pdflatex_count=5
while egrep -qs 'Rerun (Pdflatex|to get cross-references right)' ${base}.log && [ $pdflatex_count -gt 0 ]
do
    echo "Rerunning pdflatex...."
    pdflatex "$@"			>/dev/null
    pdflatex_count=`expr $pdflatex_count - 1`
done

ln -f "${outf}" ../"${outf}"
echo "Generated ${outf}"
