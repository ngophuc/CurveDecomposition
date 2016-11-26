# Curve Decomposition Algorithm

Command lines for reproducing the results in <a href="https://hal.inria.fr/hal-01375089">A discrete approach for decomposing noisy digital contours into arcs and segments</a>, submitted to DGMM4CV'16.

Results:

Fig 9:
- Line 1: File <a href="https://github.com/ngophuc/CurveDecomposition/blob/master/Samples/pentagonNoise4.sdp">pentagonNoise4.sdp</a>
<pre class="code highlight js-syntax-highlight plaintext">
<code>./testContourDecom -i ../Samples/pentagonNoise4.sdp -o ../Results/pentagonNoise4 -d ../ImaGene-forIPOL --maxScale 10 --samplingStep 1.0 --alphaMax 0.78 --thickness 0.2 --isseTol 4.0 --nbPointCircle 3</code>
</pre>
- Line 2: File <a href="https://github.com/ngophuc/CurveDecomposition/blob/master/Samples/CircleRectBis.sdp">CircleRectBis.sdp</a>
<pre class="code highlight js-syntax-highlight plaintext">
<code>./testContourDecom -i ../Samples/CircleRectBis.sdp -o ../Results/CircleRect -d ../ImaGene-forIPOL --maxScale 5 --samplingStep 0.2 --alphaMax 0.78 --thickness 0.2 --isseTol 2.0 --nbPointCircle 3</code>
</pre>
- Line 3: File <a href="https://github.com/ngophuc/CurveDecomposition/blob/master/Samples/flower100Noise4.sdp">flower100Noise4.sdp</a>
<pre class="code highlight js-syntax-highlight plaintext">
<code>./testContourDecom -i ../Samples/flower100Noise4.sdp -o ../Results/flower100Noise4 -d ../ImaGene-forIPOL --maxScale 10 --samplingStep 0.2 --alphaMax 0.78 --thickness 0.2 --isseTol 4.0 --nbPointCircle 3</code>
</pre>

Fig 10
- Line 1: File <a href="https://github.com/ngophuc/CurveDecomposition/blob/master/Samples/exemple13_KNoise_03.sdp">exemple13_KNoise_03.sdp</a>
<pre class="code highlight js-syntax-highlight plaintext">
<code>./testContourDecom -i ../Samples/exemple13_KNoise_03.sdp -o ../Results/exemple13_KNoise_03 -d ../ImaGene-forIPOL --maxScale 10 --samplingStep 1.0 --alphaMax 0.78 --thickness 0.2 --isseTol 4.0 --nbPointCircle 3</code>
</pre>
- Line 2: File <a href="https://github.com/ngophuc/CurveDecomposition/blob/master/Samples/symbol004_KNoise_05.sdp">symbol004_KNoise_05.sdp</a>: </p>
<pre class="code highlight js-syntax-highlight plaintext">
<code>./testContourDecom -i ../Samples/symbol004_KNoise_05.sdp -o ../Results/symbol004_KNoise_05 -d ../ImaGene-forIPOL --maxScale 10 --samplingStep 1.0 --alphaMax 0.78 --thickness 0.2 --isseTol 4.0 --nbPointCircle 3</code>
</pre>
- Line 3:
- Line 4:
- Line 5: File <a href="https://github.com/ngophuc/CurveDecomposition/blob/master/Samples/example7_KNoise_07.sdp">example7_KNoise_07.sdp</a>: </p>
<pre class="code highlight js-syntax-highlight plaintext">
<code>./testContourDecom -i ../Samples/example7_KNoise_07.sdp -o ../Results/example7_KNoise_07 -d ../ImaGene-forIPOL &#x000A; --maxScale 10 --samplingStep 1.0 --alphaMax 0.78 --thickness 0.2 --isseTol 4.0 --nbPointCircle 3</code>
</pre>

Fig 11

Fig 12
- Line 1:
- Line 2:
Fig 13
