#!/bin/sh -
#
# Translate a TeX document to DVI format.
#
# usage:
#	tex_dvi Name file1.tex file2.tex ...

set -e

base="$1"
outf="${base}.dvi"
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

latex_count=5
while egrep -qs 'Rerun (LaTeX|to get cross-references right)' ${base}.log && [ $latex_count -gt 0 ]
do
    echo "Rerunning latex...."
    pdflatex "$@"			>/dev/null
    latex_count=`expr $latex_count - 1`
done

ln -f "${outf}" ../"${outf}"
echo "Generated ${outf}"
