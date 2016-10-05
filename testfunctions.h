#ifndef TESTFUNCTION
#define TESTFUNCTION

#include <iostream>
#include <exception>
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/boards/Board2D.h"
#include "DGtal/io/colormaps/HueShadeColorMap.h"

#include "DGtal/geometry/curves/ArithmeticalDSSComputer.h"
#include "DGtal/geometry/curves/SaturatedSegmentation.h"
#include <DGtal/geometry/curves/AlphaThickSegmentComputer.h>

#include "myfunctions.h"

using namespace std;
using namespace DGtal;
using namespace Z2i;

#define FILENAMESIZE 200

/*** Adaptative tangent cover decomposition ***/
vector<vector<AlphaThickSegmentComputer2DD> > adaptiveTangentCoverDecomposition(const vector<vector<RealPoint> >& aContour, const vector<vector<double> >& vecMT, string filename, unsigned int w=200, unsigned int h=200);
/*** Adaptative tangent cover decomposition ***/

/**** Dominant points detections ****/
vector<vector<RealPoint> > dominantPointDetection(const vector<vector<AlphaThickSegmentComputer2DD> > &vecTangentCover, const vector<vector<RealPoint> > &aContour, string filename, unsigned int w=200, unsigned int h=200);
/**** Dominant points detections ****/

/**** Dominant points selection ****/
vector<vector<RealPoint> > dominantPointSimplification(const vector<vector<RealPoint> > &DP, const vector<vector<int> >& indexDP, const vector<vector<RealPoint> > &aContour, string filename, unsigned int w=200, unsigned int h=200);
/**** Dominant points selection ****/

/**** Tangent space transformation ****/
vector<vector<RealPoint> > tangentSpaceTransform(const vector<vector<RealPoint> > &DP, const vector<vector<RealPoint> > &aContour);
/**** Tangent space transformation ****/

/****Decomposition of Curve into Segments and Arcs ***/
vector<vector<int> > arcSegmentDecomposition(const vector<vector<RealPoint> > &aContour, const vector<vector<int> > &indexDP, const vector<vector<RealPoint> > &MP, double alphaMax, double thickness, double iseTol, double angleTol, int nbPointCir, vector<vector<double> > &segments, vector<vector<double> > &arcs);
/****Decomposition of Curve into Segments and Arcs ***/

/**** Draw decomposition from seg of seg and arcs ****/
void drawDecomposition(const vector<vector<RealPoint> >& aContour, const vector<vector<RealPoint> >& DP, const vector<vector<int> > isolatedPoint, const vector<vector<double> >& segments, const vector<vector<double> >& arcs, string filename, string onlyfilename, unsigned int w=200, unsigned int h=200);
/**** Draw decomposition from seg of seg and arcs ****/
#endif // TESTFUNCTION

