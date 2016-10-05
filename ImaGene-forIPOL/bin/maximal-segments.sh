#!/bin/bash

usage() {
    echo "$0 -x <x> -y <y> <prefix> <shape> ..."
    echo "    Analyse the maximal segments on a given shape <shape> and output XFIG and EPS files."
    echo "    - <prefix>: all the files will be prefixed with it."
    echo "    - <shape>: is some -circle <R>, ..."
    echo "    - you may add VIEW parameters as: "
    echo " VIEW=\"-view 0 0 32 16\"; $0 ..."
    echo "      to indicate which part you observe."
}

if test $# -le 6; then 
    usage
    exit 0;
fi

SPACE="$1 $2 $3 $4"
let xmax="$2"
let ymax="$4"
shift 4
PREFIX=$1
shift 1

# GENERATE FIRST QUADRANT
if test -z "$VIEW"; then
    VIEW="-view 0 0 $xmax $ymax"
fi

# Generate shape geometry
gencontour ${SPACE} $* -dFC > ${PREFIX}2.fc
cat ${PREFIX}2.fc | freeman2freeman -ul 2 2 > ${PREFIX}.fc
cat ${PREFIX}.fc | freeman2pgm ${SPACE} > ${PREFIX}.pgm
cat ${PREFIX}.fc | freeman2fig ${SPACE} ${VIEW} -header -unitcm 0.5 -contour 0 3 -pointel_grid 0 1 -pixels 1 0 -pgm_file ${PREFIX}.pgm -ms 0 0 -ms_display COLOR 58 0.5 > ${PREFIX}.fig
fig2dev -L eps ${PREFIX}.fig ${PREFIX}.eps
fig2dev -L png ${PREFIX}.fig ${PREFIX}.png

