# CurveDecomposition

Source of code associated to submitted RRPR paper.


To install the program see <a href="https://github.com/ngophuc/CurveDecomposition/blob/master/INSTALL.txt">INSTALL.txt</a> file


If you want you can also directly test the programm online:

http://ipol-geometry.loria.fr/~phuc/ipol_demo/RRPR_demo/


* [![Build Status](https://travis-ci.org/ngophuc/CurveDecomposition.svg?branch=master)](https://travis-ci.org/ngophuc/CurveDecomposition)

# Examples

<p>File <a href="https://github.com/ngophuc/CurveDecomposition/blob/master/Samples/circle50.sdp">circle50.sdp</a>: </p>&#x000A;&#x000A;
<pre class="code highlight js-syntax-highlight plaintext">
<code>./testContourDecom -i ../Samples/circle50.sdp -o ../Results/circle50 -d ../ImaGene-forIPOL &#x000A; --samplingStep 1.0 --maxScale 10 -a 0.78 -t 0.2 -n 3 -s 4.0</code>
</pre>&#x000A;&#x000A;
<p>
	<table cellpadding="5">
		<tr>
		<td align="center" valign="center">
			<a href="https://github.com/ngophuc/CurveDecomposition/blob/master/Results/circle50.pdf">
				<img width="300" src="https://github.com/ngophuc/CurveDecomposition/blob/master/Results/circle50.png" alt="Input curve" />
			</a>	
		<br />
		Input curve
		</td>
		<td align="center" valign="center">
			<a href="https://github.com/ngophuc/CurveDecomposition/blob/master/Results/circle50_OnlyArcSeg.pdf">
				<img width="300" src="https://github.com/ngophuc/CurveDecomposition/blob/master/Results/circle50_OnlyArcSeg.png" alt="Decomposition result" />
			</a>
		<br />
		Decomposition result
		</td>
		</tr>
	</table>
</p>