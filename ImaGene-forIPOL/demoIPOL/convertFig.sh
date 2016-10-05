#! /bin/sh
PATH=/usr/bin/:/usr/local/bin:/opt/local/bin


fig2dev -L eps $1 resu.eps
convert -density 300 resu.eps resu.png

