# Curve Decomposition Algorithm

Command lines for reproducing the results in <a href="https://hal.inria.fr/hal-01375089">A discrete approach for decomposing noisy digital contours into arcs and segments</a>, submitted to DGMM4CV'16.

Results:

Fig 9:
- Line 1: File <a href="https://github.com/ngophuc/CurveDecomposition/blob/master/Samples/pentagonNoise4.sdp">pentagonNoise4.sdp</a> <pre class="code highlight js-syntax-highlight plaintext">
<code>./testContourDecom -i ../Samples/pentagonNoise4.sdp -o ../Results/pentagonNoise4 -d ../ImaGene-forIPOL --maxScale 10 --samplingStep 1.0 --alphaMax 0.78 --thickness 0.2 --isseTol 4.0 --nbPointCircle 3</code>
</pre>
- Line 2: File <a href="https://github.com/ngophuc/CurveDecomposition/blob/master/Samples/RectCircle.sdp">RectCircle.sdp</a>
<pre class="code highlight js-syntax-highlight plaintext">
<code>./testContourDecom -i ../Samples/RectCircle.sdp -o ../Results/RectCircle -d ../ImaGene-forIPOL --maxScale 5 --samplingStep 0.2 --alphaMax 0.78 --thickness 0.2 --isseTol 2.0 --nbPointCircle 3</code>
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
- Line 2: File <a href="https://github.com/ngophuc/CurveDecomposition/blob/master/Samples/symbol004_KNoise_05.sdp">symbol004_KNoise_05.sdp</a>
<pre class="code highlight js-syntax-highlight plaintext">
<code>./testContourDecom -i ../Samples/symbol004_KNoise_05.sdp -o ../Results/symbol004_KNoise_05 -d ../ImaGene-forIPOL --maxScale 10 --samplingStep 1.0 --alphaMax 0.78 --thickness 0.2 --isseTol 4.0 --nbPointCircle 3</code>
</pre>
- Line 3: File <a href="https://github.com/ngophuc/CurveDecomposition/blob/master/Samples/symbol066_KNoise_05.sdp">symbol066_KNoise_05.sdp</a>
<pre class="code highlight js-syntax-highlight plaintext">
<code>./testContourDecom -i ../Samples/symbol066_KNoise_05.sdp -o ../Results/symbol066_KNoise_05 -d ../ImaGene-forIPOL --maxScale 15 --samplingStep 0.5 --alphaMax 0.78 --thickness 0.2 --isseTol 3.0 --nbPointCircle 3</code>
</pre>
- Line 4: File <a href="https://github.com/ngophuc/CurveDecomposition/blob/master/Samples/symbol119_KNoise_07.sdp">symbol119_KNoise_07.sdp</a>
<pre class="code highlight js-syntax-highlight plaintext">
<code>./testContourDecom -i ../Samples/symbol119_KNoise_07.sdp -o ../Results/symbol119_KNoise_07 -d ../ImaGene-forIPOL --maxScale 15 --samplingStep 1 --alphaMax 0.78 --thickness 0.2 --isseTol 3.0 --nbPointCircle 3</code>
</pre>
- Line 5: File <a href="https://github.com/ngophuc/CurveDecomposition/blob/master/Samples/example7_KNoise_07.sdp">example7_KNoise_07.sdp</a>
<pre class="code highlight js-syntax-highlight plaintext">
<code>./testContourDecom -i ../Samples/example7_KNoise_07.sdp -o ../Results/example7_KNoise_07 -d ../ImaGene-forIPOL --maxScale 10 --samplingStep 1.0 --alphaMax 0.78 --thickness 0.2 --isseTol 4.0 --nbPointCircle 3</code>
</pre>

Fig 11: File <a href="https://github.com/ngophuc/CurveDecomposition/blob/master/Samples/letterR.sdp">letterR.sdp</a>
<pre class="code highlight js-syntax-highlight plaintext">
<code>./testContourDecom -i ../Samples/letterR.sdp -o ../Results/letterR -d ../ImaGene-forIPOL --maxScale 20 --samplingStep 3 --alphaMax 1.0 --thickness 0.2 --isseTol 3.0 --nbPointCircle 3</code>
</pre>

Fig 12
- Line 1: File <a href="https://github.com/ngophuc/CurveDecomposition/blob/master/Samples/techImage_1.sdp">techImage_1.sdp</a>
<pre class="code highlight js-syntax-highlight plaintext">
<code>./testContourDecom -i ../Samples/techImage_1.sdp -o ../Results/techImage_1 -d ../ImaGene-forIPOL --maxScale 10 --samplingStep 1.0 --alphaMax 0.78 --thickness 0.2 --isseTol 15.0 --nbPointCircle 3</code>
</pre>
- Line 2: File <a href="https://github.com/ngophuc/CurveDecomposition/blob/master/Samples/techImage_2.sdp">techImage_2.sdp</a>
<pre class="code highlight js-syntax-highlight plaintext">
<code>./testContourDecom -i ../Samples/techImage_2.sdp -o ../Results/techImage_2 -d ../ImaGene-forIPOL --maxScale 10 --samplingStep 0.2 --alphaMax 0.78 --thickness 0.2 --isseTol 4.0 --nbPointCircle 3
</code>
</pre>

Fig 13: File <a href="https://github.com/ngophuc/CurveDecomposition/blob/master/Samples/ContourFlower.sdp">ContourFlower.sdp</a>
<pre class="code highlight js-syntax-highlight plaintext">
<code>./testContourDecom -i ../Samples/ContourFlower.sdp -o ../Results/ContourFlower -d ../ImaGene-forIPOL --maxScale 10 --samplingStep 0.2 --alphaMax 0.78 --thickness 0.2 --isseTol 4.0 --nbPointCircle 3</code>
</pre>
