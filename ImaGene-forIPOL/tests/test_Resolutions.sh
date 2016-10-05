#!/bin/sh

usage() {
    echo "$0 <samplingSizeMax> <samplingSizeAnalys>  <threshold curve/Flat> <shape> [position analysis [MEAN/MAX/MIN] ]..."
    echo "    Analyse the shape segments on a given shape <shape> and output XFIG and EPS files."
}

if test $# -le 3; then 
    usage
    exit 0;
fi

SSM="$1"
SSA="$2" 
TCF="$3"

SHAPE="$4"

POSITION=0


if test $# -ge 5; then
    POSITION="$5"
fi
STAT="MEAN"
if test $# -eq 6; then
    STAT="$6"
fi


test_Resolutions -samplingSizeMax $SSM -afficheStat $POSITION  -affichePointIndex $POSITION -afficheContourSRC -enteteXFIG -afficheContourSampled  -samplingSizeStartAnalyse ${SSA}  -afficheSlope $STAT -afficheFlat $TCF 13  -afficheCurve $TCF  0.0  11  -afficheBruit   -agrandissementEPS 5 -afficheRegLineaire  < $SHAPE > tmp.fig; fig2dev -Leps tmp.fig tmp.eps    ;



cat > pltMSL2 <<EOF
set terminal postscript epsf color 22
set output 'msl2.eps';
set logscale xy;
plot   "msl.txt" using (\$1):(\$5) title "moy" w p, "mslMoy.txt" using (\$1):(\$2) title "moyenne stat (moy)" w l , "mslSel.txt" using  (\$1):(\$5) title "selection " w l ;
set output
set terminal x11
replot
EOF

cat pltMSL2 | gnuplot -persist


cat > pltSlopes <<EOF
set terminal postscript epsf color 22
set output 'slopes.eps';
plot "slopes.txt" using (\$1):(\$2) title "slope" w l, "slopes.txt" using (\$1):(\$3) title "b" w l;
set output
set terminal x11
replot
EOF


cat pltSlopes | gnuplot -persist





gv tmp.eps