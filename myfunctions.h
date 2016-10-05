#ifndef MYFUNCTIONS
#define MYFUNCTIONS

#include "DGtal/helpers/StdDefs.h"
#include "DGtal/base/Common.h"
#include "DGtal/io/boards/Board2D.h"
#include <stdio.h>
#include <numeric>
#include <algorithm>
#include <list>
#include <vector>
#include <iostream>
#include <fstream>
#include <exception>

using namespace std;
using namespace DGtal;
using namespace Z2i;

typedef std::vector<Point> Range;
typedef std::vector<Point>::const_iterator ConstIterator;
typedef std::vector<RealPoint>::const_iterator ConstIteratorD;
typedef ArithmeticalDSS<ConstIterator,int,4> DSS4;

typedef  AlphaThickSegmentComputer< Point > AlphaThickSegmentComputer2D;
typedef  AlphaThickSegmentComputer<RealPoint> AlphaThickSegmentComputer2DD;

struct MyDrawStyleCustomColor : public DrawableWithBoard2D
{
  Color myPenColor;
  Color myFillColor;
  MyDrawStyleCustomColor( const Color & penColor,
        const Color & fillColor )
    : myPenColor( penColor ), myFillColor( fillColor )
  {}
  virtual void setStyle( Board2D & aboard) const
  {
    aboard.setFillColor( myFillColor);
    aboard.setPenColor( myPenColor );
  }
};

typedef PointVector<2, double> RealVector;
struct Value
{
  RealVector first;
  double second;
  Value () : second ( 0. ) {}
  Value & operator+= ( const Value & ch )
  {
    this->first += ch.first;
    this->second += ch.second;
    return *this;
  }
};

double sort_increase(double a, double b);
vector<int> absSortIndex(vector<double> const& values, bool isIncrease=true);

vector<vector<RealPoint> > getContours(const vector<vector<RealPoint> > &vecPoint);

vector<double> readMeanindfulThicknessFile(const char* filename);

void writeFile(const vector<RealPoint>& v, const char* filename);

Point getStartPoint(const AlphaThickSegmentComputer2D s);
Point getEndPoint(const AlphaThickSegmentComputer2D s);
RealPoint getStartPoint(const AlphaThickSegmentComputer2DD s);
RealPoint getEndPoint(const AlphaThickSegmentComputer2DD s);

double getSlope(double Nx,double Ny);//(Nx,Ny) : normal vector of the line

double acuteAngle(RealPoint bp1, RealPoint bp2, RealPoint bp3); //acute angle between three points
double signedAngle(RealPoint v1, RealPoint v2); //angle between two vectors
double relativeAngle(RealPoint v1, RealPoint v2); //relative angle between two vectors
double relativeAngle(RealPoint bp1, RealPoint bp2, RealPoint bp3); //relative angle between three points

double distancePointSegment(RealPoint p, RealPoint s1, RealPoint s2);
double distancePoints(RealPoint p1, RealPoint p2);
double distancePointCircle(Point p, RealPoint center, double radius);
double signedDistancePointCircle(Point p, RealPoint center, double radius);

int findElement(const vector<RealPoint>& vec, RealPoint p);
int findElement(const vector<RealPoint>& vec, RealPoint p, int start);
int findStartElement(const vector<RealPoint>& vec, const AlphaThickSegmentComputer2DD s);
int findEndElement(const vector<RealPoint>& vec, const AlphaThickSegmentComputer2DD s);

double error_CR(const vector<RealPoint>& contour, const vector<RealPoint>& DP);
double error_ISE(const vector<RealPoint>& contour, const vector<RealPoint>& DP, const vector<int>& indexDP, bool isClosed=false);
double error_L_infini(const vector<RealPoint>& contour, const vector<RealPoint>& DP, const vector<int>& indexDP, bool isClosed=false);
double error_FOM(const vector<RealPoint>& contour, const vector<RealPoint>& DP, const vector<int>& indexDP, bool isClosed=false);
double error_FOM_M(const vector<RealPoint>& contour, const vector<RealPoint>& DP, const vector<int>& indexDP, bool isClosed=false);
double error_FOM_ND(const vector<RealPoint>& contour, const vector<RealPoint>& DP, const vector<int>& indexDP, bool isClosed=false);
void error_All(const vector<RealPoint>& contour, const vector<RealPoint>& DP, const vector<int>& indexDP, bool isClosed=false);

int isLeft(RealPoint p1, RealPoint p2, RealPoint p3);

RealPoint determineCenter(RealPoint p1, RealPoint p2, RealPoint p3);
double determineRadius(RealPoint p1, RealPoint p2, RealPoint p3);

RealPoint determineCenter(Point p1, Point p2, Point p3);
RealPoint determineCenter(Point p1, Point p2, double radius, bool negatif);
double determineRadius(Point p1, Point p2, Point p3);
double determineRadius(RealPoint centre, Point p);
double arcLength(Point p1, Point p2, Point p3);

double iseContourCircle(const vector<RealPoint>& contour, Point p1, Point p2, RealPoint center, double radius);//ise of points between two points p1,p2 w.r.t the circle (center,radius)
double iseContourCircle(const vector<RealPoint>& contour, int indexP1, int indexP2, RealPoint center, double radius);//ise of points between two points of index idp1,idp2 w.r.t the circle (center,radius)
double lmaxContourCircle(const vector<RealPoint>& contour, Point p1, Point p2, RealPoint center, double radius);//lmax
double lmaxContourCircle(const vector<RealPoint>& contour, int indexP1, int indexP2, RealPoint center, double radius);//lmax

double iseContourSegment(const vector<RealPoint>& contour, Point p1, Point p2);//ise of points between two points p1,p2 w.r.t the circle (center,radius)
double iseContourSegment(const vector<RealPoint>& contour, int indexP1, int indexP2);//ise of points between two points of index idp1,idp2 w.r.t the circle (center,radius)
double lmaxContourSegment(const vector<RealPoint>& contour, Point p1, Point p2);//lmax
double lmaxContourSegment(const vector<RealPoint>& contour, int indexP1, int indexP2);//lmax

#endif // MYFUNCTIONS

