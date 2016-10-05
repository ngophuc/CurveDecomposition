#! /bin/sh
PATH=/usr/bin/:/usr/local/bin:/opt/local/bin
ARG=$*
$(pwd)/../../bin/meaningfulScaleEstim ${ARG} > noiseLevels.txt 2>> info.txt

 


