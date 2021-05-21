#!/bin/sh

if test $# -lt 1; 
then
    echo "$0 <file1>"
    echo "\t - view the generated DMLP for that file."
    exit 0
fi

j=$1
echo 'plot "XXX.pts" using 1:2:($4-$1):($5-$2) title "Contour XXX" with vectors head filled lt 1,"XXX-dmlp-ccw.pts" using 1:2:($4-$1):($5-$2) title "DMLP CCW" with vectors head filled lt 2, "XXX-dmlp-ccw.pts" using 1:2:($3/4):($3/4) notitle with boxxyerrorbars 2,"XXX-dmlp-cw.pts" using 1:2:($4-$1):($5-$2) title "DMLP CW" with vectors head filled lt 3, "XXX-dmlp-cw.pts" using 1:2:($3/4):($3/4) notitle with boxxyerrorbars 3'| sed 's@XXX@'$j'@g' | gnuplot -persist
