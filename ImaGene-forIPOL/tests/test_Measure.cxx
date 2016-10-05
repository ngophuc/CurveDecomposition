/**
 * @file testr_measure.cpp
 * @author David Coeurjolly (\c david.coeurjolly@liris.cnrs.fr )
 *
 *
 * @date 2010/03/03
 *
 * This file is part of the ImaGene library
 */

/**
 * Description of test_trace' <p>
 * Aim: simple test of \ref MeasureOfStraighLines
 */

#include <cstdio>
#include <cmath>
#include <iostream>
#include <fstream>
#include "ImaGene/mathutils/MeasureOfStraightLines.h"
#include "ImaGene/dgeometry2d/C4CSegment.h"

using namespace ImaGene;
using namespace std;



/**
 * Compute the measure of the unit square [0,1]x[0,1]
 *
 * Expected value : sqrt{2}/2
 **/
void testUnitSquare()
{
  vector<double> a;
  vector<double> b;
  
  a.push_back(0); b.push_back(0);
  a.push_back(1); b.push_back(0);
  a.push_back(1); b.push_back(1);
  a.push_back(0); b.push_back(1);
  
  MeasureOfStraightLines measure;

  cout << "Measure of the Straight of Lines of the unit square = " << measure.computeMeasure(a,b)<< std::endl;
  cout <<"Expected value = 0.707107 (sqrt(2)/2)"<<endl;
}


/**
 * Compute the measure of the unit square [0,1]x[0,1]
 *
 * Expected value : sqrt{2}/2
 **/
void testUnitSquareCentroid()
{
  vector<double> a;
  vector<double> b;
  
  a.push_back(0); b.push_back(0);
  a.push_back(1); b.push_back(0);
  a.push_back(1); b.push_back(1);
  a.push_back(0); b.push_back(1);
  
  MeasureOfStraightLines measure;

  cout << "Centroid measure of the unit square = (" << measure.computeCentroidA(a,b)
	       << ","<<measure.computeCentroidB(a,b)<<")"<<std::endl;
  cout <<"Expected value = (0.4142,0.5)"<<endl;
}


/**
* Test C4CSegment link
*
**/
void testC4CSegement()
{
  C4CSegment c;

  c.init();
  //Small segmnet "--|--|"
  c.extendsFront(2);
  c.extendsFront(2);
  c.extendsFront(3);
  c.extendsFront(2);
  c.extendsFront(2);
  c.extendsFront(3);

  c.selfDisplay(cout);
  
  cout << "C4CSegment Measure = " << c.getMeasure()<<std::endl;
  cout << "C4CSegment slope = " << c.getTangentMeasure()<<std::endl;
}

int main(int argc, char **argv)
{
  
  testUnitSquare();
  testUnitSquareCentroid();
  testC4CSegement();
  return 0;
}

