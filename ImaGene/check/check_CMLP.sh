#!/bin/bash

CMLPBIN='../build/bin/freeman2polygon'
file_prefix='.tmp-'

usage(){
    echo "Usage: $0 <file1.fc> ..."
}


check_cmlp() {
    if test ! -r $1; then
	echo "File $1 is not readable."
	exit 1
    fi

    fc=`cat $1 | grep -v '^#'`
    x=`echo $fc | cut -d ' ' -f 1`
    y=`echo $fc | cut -d ' ' -f 2`
    chain=`echo $fc | cut -d ' ' -f 3`
    len=`echo $chain | awk '{ print(length($0)); }'`

# echo "x=$x"
# echo "y=$y"
# echo "c=$chain"
# echo "l=$len"

    let i=1
    let nbok=0
    let nbko=0
    err=0.0
    max_err=0.0
    cmlp=`${CMLPBIN} -input $1 -lenCMLP -auto_center 2> /dev/null | grep -v '^#'`
    while test $i -le $len; do
	start_c=`echo $chain | awk '{ print(substr($0,'$i')); }'`
	let j=i-1
	end_c=`echo $chain | awk '{ print(substr($0,0,'$j')); }'`
	circc=`echo "${start_c}${end_c}"`
	filename="${file_prefix}$i"
	let i=i+1
	echo "$x $y $circc" > ${filename}
	cmlplen=`${CMLPBIN} -input ${filename} -lenCMLP -auto_center 2> /dev/null | grep -v '^#'`
#     echo "-- $i: $cmlp $cmlplen"
	if test "$cmlp" = "$cmlplen"; then
	    let nbok=nbok+1
	else
	    let nbko=nbko+1
	fi
	max_err=`echo "define abs(x) {if (x>=0) return x else return -x } if (abs($cmlp-$cmlplen)>$max_err) abs($cmlp-$cmlplen) else $max_err" | bc -l`
	err=`echo "define abs(x) {if (x>=0) return x else return -x } $err+abs($cmlp-$cmlplen)" | bc -l`
    done
    let nb=nbok+nbko
    echo "-- $1 --- $nbok/$nb --- len_cmlp=$cmlp --- total_err=$err --- max_err=$max_err"
}

if test $# = 0; then
    usage
    exit 0
fi
if test ! -x $CMLPBIN; then
    echo "File $CMLPBIN is not executable."
    exit 1
fi

for fc in $*; do
    check_cmlp $fc
done