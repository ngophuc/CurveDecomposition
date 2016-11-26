# Curve Decomposition Algorithm

Command lines for reproducing the results in <a href="https://hal.inria.fr/hal-01375089">A discrete approach for decomposing noisy digital contours into arcs and segments</a>, submitted to DGMM4CV'16.

Results:
Fig 9:
- Line 1: File <a href="https://github.com/ngophuc/CurveDecomposition/blob/master/Samples/pentagonNoise4.sdp">pentagonNoise4.sdp</a>
<pre class="code highlight js-syntax-highlight plaintext">
<code>./testContourDecom -i ../Samples/pentagonNoise4.sdp -o ../Results/pentagonNoise4 -d ../ImaGene-forIPOL &#x000A; --maxScale 10 --samplingStep 1.0 --alphaMax 0.78 --thickness 0.2 --isseTol 4.0 --nbPointCircle 3</code>
</pre>&#x000A;&#x000A;
- Line 2: <p>File <a href="https://github.com/ngophuc/CurveDecomposition/blob/master/Samples/CircleRectNew.sdp">CircleRect.sdp</a>
<code>./testContourDecom -i ../Samples/CircleRectNew.sdp -o ../Results/CircleRect -d ../ImaGene-forIPOL --maxScale 5 --samplingStep 0.2 --alphaMax 0.78 --thickness 0.2 --isseTol 2.0 --nbPointCircle 3</code>
</pre>&#x000A;&#x000A;
- Line 3:

Fig 10
