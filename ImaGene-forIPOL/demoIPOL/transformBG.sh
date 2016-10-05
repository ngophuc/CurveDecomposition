#! /bin/sh
PATH=/usr/bin/:/usr/local/bin:/opt/local/bin

convert +contrast +contrast +contrast +contrast +contrast -modulate 160,100,100  -type grayscale -depth 8  input_0.png input_0BG.png                                 

                                      
#OLD does not work on previous versions:
#convert  -brightness-contrast 40x-40 input_0.png input_0BG.png



