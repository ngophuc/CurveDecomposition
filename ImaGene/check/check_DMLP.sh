#!/bin/sh

if test $# -lt 1; 
then
    echo "$0 <file1.fc> [<file2.fc> ...]"
    echo "\t - checks that DMLP is correct for all given input digital shapes."
    exit 0
fi

DMLPBIN=../build/tests/test_DMLP
let n=0;
let nok=0;
for i in $*; do
    echo "- $i ";
    j=`basename $i | sed 's/.fc//'`
    ${DMLPBIN} -input "$i" -orientation CCW -dPTS 2> /dev/null > "output/$j.pts"
    ${DMLPBIN} -input "$i" -orientation CCW -mlp -dPTS 2> "output/$j-dmlp-ccw.txt" > "output/$j-dmlp-ccw.pts"
    ${DMLPBIN} -input "$i" -orientation CW -mlp -dPTS 2> "output/$j-dmlp-cw.txt" > "output/$j-dmlp-cw.pts"
    # close contours
    cat "output/$j-dmlp-ccw.pts" | head -3 | tail -1 >> "output/$j-dmlp-ccw.pts"
    cat "output/$j-dmlp-cw.pts" | head -3 | tail -1 >> "output/$j-dmlp-cw.pts"
    diff -q "output/$j-dmlp-ccw.pts" "expected/$j-dmlp-ccw.pts"
    if test $? -eq 0; then
	echo "--- $i cw : OK";
	let nok=nok+1
    else echo "--- $i cw : ERREUR";
    fi
    let n=n+1
    diff -q "output/$j-dmlp-ccw.pts" "expected/$j-dmlp-ccw.pts"
    if test $? -eq 0; then
	echo "--- $i ccw: OK";
	let nok=nok+1
    else echo "--- $i ccw: ERREUR";
    fi
    let n=n+1
done
echo "- OK : $nok/$n"
