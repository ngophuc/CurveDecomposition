#!/bin/bash

usage() {
    echo "$0 <prefix> <shape> ..."
    echo "    Compare the respective errors of CC and GMC curvature estimation on a given shape <shape> by calling error_analysis."
    echo "    - <prefix>: all the files will be prefixed with it."
    echo "    - <shape>: is some -circle <R>, ..."
    echo "    - you may add -x <x> -y <y> to specify a bigger discrete space."
}

if test $# -le 2; then 
    usage
    exit 0;
fi

PREFIX=$1
shift 1

# Generate shape geometry
gencontour $* -dSG > ${PREFIX}.geom

# generate Freeman chain code
gencontour $* -dFC > ${PREFIX}.fc

# generate PGM
gencontour $* -dPGM > ${PREFIX}.pgm
display ${PREFIX}.pgm &

###############################################################################
# GMC curvature estimator.
###############################################################################

# estimate curvature and length with GMC and l-MST
EPS="0.00000001 -1"
echo "# GMC epsilon is ${EPS}"
cat ${PREFIX}.fc | curvature_gmc -minimizer RLX 0 -eps ${EPS} > ${PREFIX}-gmc.geom

###############################################################################
# CC curvature estimator.
###############################################################################

# estimate curvature and length with GMC and l-MST
cat ${PREFIX}.fc | curvature_cc > ${PREFIX}-cc.geom


join ${PREFIX}-gmc.geom ${PREFIX}.geom > ${PREFIX}-gmc-ext.geom
join ${PREFIX}-cc.geom ${PREFIX}.geom > ${PREFIX}-cc-ext.geom

echo > ${PREFIX}.err-gmc
echo "############################################" >> ${PREFIX}.err-gmc
echo "######  GMC ERROR ANALYSIS             #####" >> ${PREFIX}.err-gmc
echo "############################################" >> ${PREFIX}.err-gmc

cat ${PREFIX}-gmc-ext.geom | error_analysis -expected_idx 10 -estimated_idx 2 -measure_idx 16 >> ${PREFIX}.err-gmc

echo > ${PREFIX}.err-cc
echo "############################################" >> ${PREFIX}.err-cc
echo "######   CC ERROR ANALYSIS             #####" >> ${PREFIX}.err-cc
echo "############################################" >> ${PREFIX}.err-cc

cat ${PREFIX}-cc-ext.geom | error_analysis -expected_idx 10 -estimated_idx 2 -measure_idx 16 >> ${PREFIX}.err-cc

cat ${PREFIX}.err-gmc
cat ${PREFIX}.err-cc

cat > ${PREFIX}-curv-gmc-cc.gnuplot <<EOF
set terminal postscript epsf 14
set output "${PREFIX}-curv-gmc-cc.eps"
plot "${PREFIX}-gmc.geom" using 1:2 title "Estimated curvature by GMC" with lines, "${PREFIX}-cc.geom" using 1:2 title "Estimated curvature by CC" with lines, "${PREFIX}.geom" using 1:4 title "Expected curvature" with lines
set output
set terminal x11
replot
EOF

cat ${PREFIX}-curv-gmc-cc.gnuplot | gnuplot -persist
