#!/bin/bash

usage() {
    echo "$0 <prefix> <shape> ..."
    echo "    Compare the respective qualities of CC and GMC curvature estimation on a given shape <shape> by reintegrating this estimated geometric information into a curve in the plane."
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

# Generate curve by integrating geometric information (curvature and length)
cat ${PREFIX}.geom | integrate_shape -curv_idx 4 -cabs_idx 9 -t0 1.57079632679489661922 \
    | sdp2pts -centering > ${PREFIX}-int.pts

# generate "true shape"
gencontour $* -dS 2000  \
    | sdp2pts -centering > ${PREFIX}-true.pts

# generate "discrete shape"
gencontour $* -dFC | freeman2sdp \
    | sdp2pts -centering > ${PREFIX}-digital.pts

# generate Freeman chain code
gencontour $* -dFC > ${PREFIX}.fc


###############################################################################
# GMC curvature estimator.
###############################################################################

# estimate curvature and length with GMC and l-MST
cat ${PREFIX}.fc | curvature_gmc -minimizer RLX 0 -eps 0.0000001 -1 > ${PREFIX}-gmc.geom

# Generate curve by integrating geometric information (GMC curvature and l-MST length)
cat ${PREFIX}-gmc.geom | integrate_shape -curv_idx 2 -cabs_idx 4 -t0 1.57079632679489661922 \
    | sdp2pts -centering > ${PREFIX}-gmc-intlmst.pts

# Generate curve by integrating geometric information (GMC curvature and GMC-theta length)
cat ${PREFIX}-gmc.geom | integrate_shape -curv_idx 2 -cabs_idx 7 -t0 1.57079632679489661922 \
    | sdp2pts -centering > ${PREFIX}-gmc-inttheta.pts

# Generate curve by integrating geometric information (GMC curvature and GMC-theta length)
join ${PREFIX}-gmc.geom ${PREFIX}.geom > ${PREFIX}-gmc-ext.geom
cat ${PREFIX}-gmc-ext.geom | integrate_shape -curv_idx 2 -cabs_idx 15 -t0 1.57079632679489661922 \
    | sdp2pts -centering > ${PREFIX}-gmc-inttrue.pts


###############################################################################
# CC curvature estimator.
###############################################################################

# estimate curvature and length with GMC and l-MST
cat ${PREFIX}.fc | curvature_cc > ${PREFIX}-cc.geom

# Generate curve by integrating geometric information (CC curvature and l-MST length)
cat ${PREFIX}-cc.geom | integrate_shape -curv_idx 2 -cabs_idx 4 -t0 1.57079632679489661922 \
    | sdp2pts -centering > ${PREFIX}-cc-intlmst.pts

# Generate curve by integrating geometric information (GMC curvature and GMC-theta length)
cat ${PREFIX}-cc.geom | integrate_shape -curv_idx 2 -cabs_idx 7 -t0 1.57079632679489661922 \
    | sdp2pts -centering > ${PREFIX}-cc-inttheta.pts

# Generate curve by integrating geometric information (GMC curvature and GMC-theta length)
join ${PREFIX}-cc.geom ${PREFIX}.geom > ${PREFIX}-cc-ext.geom
cat ${PREFIX}-cc-ext.geom | integrate_shape -curv_idx 2 -cabs_idx 15 -t0 1.57079632679489661922 \
    | sdp2pts -centering > ${PREFIX}-cc-inttrue.pts



cat > ${PREFIX}-gmc.gnuplot <<EOF
set terminal postscript epsf 14
set output "${PREFIX}-gmc.eps"
plot "${PREFIX}-gmc-intlmst.pts" using 1:2 title "Curve by integration of estimated GMC curvature and estimated l-MST length" with lines, "${PREFIX}-gmc-inttheta.pts" using 1:2 title "Curve by integration of estimated GMC curvature and estimated GMC-theta length" with lines, "${PREFIX}-gmc-inttrue.pts" using 1:2 title "Curve by integration of estimated GMC curvature and true arclength" with lines, "${PREFIX}-int.pts" using 1:2 title "Curve by integration of curvature and length" with lines, "${PREFIX}-digital.pts" using 1:2 title "Digital curve" with lines, "${PREFIX}-true.pts" using 1:2 title "True curve" with lines 
set output
set terminal x11
replot
EOF

cat ${PREFIX}-gmc.gnuplot | gnuplot -persist


cat > ${PREFIX}-cc.gnuplot <<EOF
set terminal postscript epsf 14
set output "${PREFIX}-cc.eps"
plot "${PREFIX}-cc-intlmst.pts" using 1:2 title "Curve by integration of estimated CC curvature and estimated l-MST length" with lines, "${PREFIX}-cc-inttheta.pts" using 1:2 title "Curve by integration of estimated CC curvature and estimated CC-theta length" with lines, "${PREFIX}-cc-inttrue.pts" using 1:2 title "Curve by integration of estimated CC curvature and true arclength" with lines, "${PREFIX}-int.pts" using 1:2 title "Curve by integration of curvature and length" with lines, "${PREFIX}-digital.pts" using 1:2 title "Digital curve" with lines, "${PREFIX}-true.pts" using 1:2 title "True curve" with lines 
set output
set terminal x11
replot
EOF

cat ${PREFIX}-cc.gnuplot | gnuplot -persist

cat > ${PREFIX}-gmc-cc.gnuplot <<EOF
set terminal postscript epsf 14
set output "${PREFIX}-gmc-cc.eps"
plot "${PREFIX}-gmc-intlmst.pts" using 1:2 title "Curve by integration of estimated GMC curvature and estimated l-MST length" with lines, "${PREFIX}-cc-inttrue.pts" using 1:2 title "Curve by integration of estimated CC curvature and true arclength" with lines, "${PREFIX}-digital.pts" using 1:2 title "Digital curve" with lines
set output
set terminal x11
replot
EOF

cat ${PREFIX}-gmc-cc.gnuplot | gnuplot -persist
