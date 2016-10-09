# CurveDecomposition

Source of code associated to submitted RRPR paper.


To install the program see <a href="https://github.com/ngophuc/CurveDecomposition/blob/master/INSTALL.txt">INSTALL.txt</a> file


If you want you can also directly test the programm online:

http://ipol-geometry.loria.fr/~phuc/ipol_demo/RRPR_demo/


* [![Build Status](https://travis-ci.org/ngophuc/CurveDecomposition.svg?branch=master)](https://travis-ci.org/ngophuc/CurveDecomposition)

# Examples

<p>File <a href="https://github.com/ngophuc/CurveDecomposition/blob/master/Samples/circle50.sdp">circle50.sdp</a>: </p>&#x000A;&#x000A;
<pre class="code highlight js-syntax-highlight plaintext">
	<code>./testTriangle -i ../Samples/flower-30-8-3.sdp  -a 30  -c  -o flowerCDTC30.eps --holes ../Samples/flower-30-8-3Holes.sdp&#x000A;ps2pdf -dEPSCrop flowerCDTC30.eps flowerCDTC30.pdf&#x000A;</code>
</pre>&#x000A;&#x000A;
<p>
	<table cellpadding="5">
		<tr>
		<td align="center" valign="center">
			<a href="https://github.com/ngophuc/CurveDecomposition/blob/master/Samples/circle50.png">
				<img width="300" src="https://github.com/ngophuc/CurveDecomposition/blob/master/Samples/circle50.png" alt="Input curve" />
			</a>	
		<br />
		Input curve
		</td>

		<td align="center" valign="center">
			<a href="https://github.com/ngophuc/CurveDecomposition/blob/master/Samples/circle50.png">
				<img width="300" src="https://github.com/ngophuc/CurveDecomposition/blob/master/Samples/circle50.png" alt="Dominant points" />
			</a>
		<br />
		Dominant point detected 
		</td>

		<td align="center" valign="center">
			<a href="https://github.com/ngophuc/CurveDecomposition/blob/master/Samples/circle50.png">
				<img width="300" src="https://github.com/ngophuc/CurveDecomposition/blob/master/Samples/circle50.png" alt="Decomposition result" />
			</a>
		<br />
		Decomposition result
		</td>
		</tr>
	</table>
</p>