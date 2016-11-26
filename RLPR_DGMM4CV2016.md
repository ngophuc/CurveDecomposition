# Curve Decomposition Algorithm

Command lines for reproducing the results in <a href="https://hal.inria.fr/hal-01375089">A discrete approach for decomposing noisy digital contours into arcs and segments</a>, submitted to DGMM4CV'16.

Results:
Fig 9:
- Line 1: <p>File <a href="https://github.com/ngophuc/CurveDecomposition/blob/master/Samples/pentagonNoise4.sdp">pentagonNoise4.sdp</a>: </p>&#x000A;&#x000A;
<pre class="code highlight js-syntax-highlight plaintext">
<code>./testContourDecom -i ../Samples/pentagonNoise4.sdp -o ../Results/pentagonNoise4 -d ../ImaGene-forIPOL &#x000A; --maxScale 10 --samplingStep 1.0 --alphaMax 0.78 --thickness 0.2 --isseTol 4.0 --nbPointCircle 3</code>
</pre>&#x000A;&#x000A;
- Line 2:
- Line 3:
Fig 10