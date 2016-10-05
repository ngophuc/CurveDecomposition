#!/bin/bash

usage() {
    echo "$0 <curvature_visualc parameters>"
    echo "    Creates the simplified polygon that is constructed from the multiscale digital visual curvature (see [Liu, Latecki,Liu IJCV 2008]). The user must provide a Freemanchain in the standard input. Outputs a list of points that is a subset of the input digital points."
    echo "    - <curvature_visualc parameters>: all the parameters are handled to the program curvature_visualc."
    echo "    - Example: -radius 5 -window 4 -scale 0.01"
}

if test $# -le 1; then 
    usage
    exit 0;
fi

curvature_visualc $* \
    | grep -v '#' \
    | awk '{ if ($4 != 0) printf("%d %d\n", $6, $7); }'