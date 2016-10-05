#!/bin/sh

if test $# -lt 1; 
then
    echo "$0 <file1.fc> [<file2.fc> ...]"
    echo "\t - tests that DMLP is correct for all given input digital shapes."
    exit 0
fi

for i in $*; do
    echo $i
    j=`echo $i | sed 's/.fc//'`
    ./tests/test_DMLP-debug -input "$i" -orientation CCW -dPTS 2> /dev/null > "$j.pts"
    ./tests/test_DMLP-debug -input "$i" -orientation CCW -trace 1 -debug -mlp -dPTS 2> "$j-dmlp-ccw.txt" > "$j-dmlp-ccw.pts"
    ./tests/test_DMLP-debug -input "$i" -orientation CW -trace 1 -debug -mlp -dPTS 2> "$j-dmlp-cw.txt" > "$j-dmlp-cw.pts"
    # close contours
    cat "$j.pts" | head -3 | tail -1 >> "$j.pts"
    cat "$j-dmlp-ccw.pts" | head -3 | tail -1 >> "$j-dmlp-ccw.pts"
    cat "$j-dmlp-cw.pts" | head -3 | tail -1 >> "$j-dmlp-cw.pts"
#    echo 'plot "XXX.pts" using 1:2 title "Contour XXX" w l 1, "XXX-dmlp-ccw.pts" using 1:2 title "DMLP CCW" w lp 2, "XXX-dmlp-ccw.pts" using 1:2:($3/4):($3/4) notitle with boxxyerrorbars 2,"XXX-dmlp-cw.pts" using 1:2 title "DMLP CW" w lp 3, "XXX-dmlp-cw.pts" using 1:2:($3/4):($3/4) notitle with boxxyerrorbars 3'| sed 's/XXX/'$j'/g' | gnuplot -persist
    echo 'plot "XXX.pts" using 1:2:($4-$1):($5-$2) title "Contour XXX" with vectors head filled lt 1,"XXX-dmlp-ccw.pts" using 1:2:($4-$1):($5-$2) title "DMLP CCW" with vectors head filled lt 2, "XXX-dmlp-ccw.pts" using 1:2:($3/4):($3/4) notitle with boxxyerrorbars 2,"XXX-dmlp-cw.pts" using 1:2:($4-$1):($5-$2) title "DMLP CW" with vectors head filled lt 3, "XXX-dmlp-cw.pts" using 1:2:($3/4):($3/4) notitle with boxxyerrorbars 3'| sed 's/XXX/'$j'/g' | gnuplot -persist
done

