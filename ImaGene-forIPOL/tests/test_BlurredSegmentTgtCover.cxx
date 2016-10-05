///////////////////////////////////////////////////////////////////////////////
// Test .
///////////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include <iostream>
#include <iomanip>
#include "ImaGene/mathutils/Mathutils.h"
#include "ImaGene/base/Proxy.h"
#include "ImaGene/base/Vector.h"
#include "ImaGene/base/Vector2i.h"
#include "ImaGene/base/Arguments.h"
#include "ImaGene/base/StandardArguments.h"
#include "ImaGene/base/Shapes.h"
#include "ImaGene/base/VectorUtils.h"
#include "ImaGene/timetools/Clock.h"
#include "ImaGene/dgeometry2d/C4CSegment.h"
#include "ImaGene/dgeometry2d/C4CIterator.h"
#include "ImaGene/dgeometry2d/C4CGeometry.h"
#include "ImaGene/dgeometry2d/C4CSegmentPencil.h"
#include "ImaGene/dgeometry2d/C4CTangentialCover.h"
#include "ImaGene/digitalnD/KnCharSet.h"
#include "ImaGene/digitalnD/BelAdjacency.h"
#include "ImaGene/digitalnD/ObjectBoundary.h"
#include "ImaGene/digitalnD/KnShapes.h"
#include "ImaGene/digitalnD/Embedder.h"
#include "ImaGene/digitalnD/GridEmbedder.h"
#include "ImaGene/digitalnD/Frame2D.h"
#include "ImaGene/digitalnD/KnSpaceScanner.h"
#include "ImaGene/helper/ShapeHelper.h"
#include "ImaGene/helper/C4CTangentialCoverGeometry.h"
#include "ImaGene/helper/GlobalC4CGeometry.h"
#include "ImaGene/helper/DrawingXFIG.h"
#include "ImaGene/mathutils/Statistics.h"
#include "ImaGene/base/Signal.h"
#include "ImaGene/dgeometry2d/BlurredSegmentTgtCover.h"
#include "ImaGene/helper/ScaleProfile.h"


using namespace std;
using namespace ImaGene;

static Arguments args;
static int fillIntensity = 20;
static const int nbIterationSpikes = 10;

// /**
//  * @return the four points of each tangential cover of a contour. 
//  * @param vectPointsContour is a set of the pixels of a contour.
//  * @param i is the index of the tangential cover.
//  */

// vector<Vector2D> ComputeRealBounds( const std::vector<Vector2D> & vectPointsContour,
//                                     const vector<BlurredSegmentTgtCover::Segment> & m_Segments, 
//                                     int thick, 
//                                     int i )
// {
//   float x,y;
//   std::vector<Vector2D> vectPointsBound;
//   Vector2D fp ( vectPointsContour[m_Segments[i].m_first].x() , vectPointsContour[m_Segments[i].m_first].y() );
//   Vector2D lp ( vectPointsContour[m_Segments[i].m_last].x() , vectPointsContour[m_Segments[i].m_last].y() );
//   cerr << "first_point = " << m_Segments[i].m_first << "   last_point = " << m_Segments[i].m_last << endl;
//   cerr << "fp(" << fp.x() << "," << fp.y() << ")" << endl;
//   cerr << "lp(" << lp.x() << "," << lp.y() << ")" << endl;
//   int a1 = m_Segments[i].m_a;
//   int b1 = m_Segments[i].m_b;
//   int mu1 = m_Segments[i].m_mu;
//   int omega1 = m_Segments[i].m_omega;
// //  cout << "mu = " << mu1 << "  omega = " << omega1 << endl;
//   int mu_sup = mu1 + omega1;
//   //cout << "a1 = " << a1 << "  b1 = " << b1 << "   mu1 = " << mu1 << "   mu_sup = " << mu_sup << endl;
//   int a2 = - b1;
//   int b2 = a1;
//   int mfp = a2 * fp.x() - b2 * fp.y();
//   int mlp = a2 * lp.x() - b2 * lp.y();
//   //cout << "a2 = " << a2 << "  b2 = " << b2 << "   mfp = " << mfp << "   mlp = " << mlp << endl;
//   // coordinate of P0.
//   cout << fixed << setprecision (2);
//   x = (float)( mu1 * b2 - b1 * mfp ) / ( a1 * b2 - a2 * b1 ); 
//   y = (float)( b1 != 0 ? ( a1 * x - mu1 ) / b1 : fp.y() );
//   //cout << "x0 = " << x << "  y0 = " << y << endl; 
//   Vector2D P0( x , y );
//   vectPointsBound.push_back( P0 );
//   // coordinate of P1.
//   x = (float)( mu_sup * b2 - b1 * mfp ) / ( a1 * b2 - a2 * b1 ); 
//   y = (float)( b1 != 0 ? ( a1 * x - mu_sup ) / b1 : fp.y() );
//   Vector2D P1( x , y );
//   vectPointsBound.push_back( P1 );
//   // coordinate of P2.
//   x = (float)( mu_sup * b2 - b1 * mlp ) / ( a1 * b2 - a2 * b1 ); 
//   y = (float)( b1 != 0 ? ( a1 * x - mu_sup ) / b1 : lp.y() );
//   Vector2D P2( x , y );
//   vectPointsBound.push_back( P2 );
//   // coordinate of P3.
//   x = (float)( mu1 * b2 - b1 * mlp ) / ( a1 * b2 - a2 * b1 ); 
//   y = (float)( b1 != 0 ? ( a1 * x - mu1 ) / b1 : lp.y() );
//   Vector2D P3( x , y );
//   vectPointsBound.push_back( P3 );
//   cerr << "i = " << i << endl;
//   cerr << "P0(" << P0.x() << "," << P0.y() << ")" << endl;
//   cerr << "P1(" << P1.x() << "," << P1.y() << ")" << endl;
//   cerr << "P2(" << P2.x() << "," << P2.y() << ")" << endl;
//   cerr << "P3(" << P3.x() << "," << P3.y() << ")" << endl;
//   return vectPointsBound;
// }

// vector<Vector2D> ComputeBasicBounds( const std::vector<Vector2D> & vectPointsContour,
//                                      const vector<BlurredSegmentTgtCover::Segment> & m_Segments, 
//                                      int thick, 
//                                      int i )
// {
//   vector<Vector2D> VectPointsBound = ComputeRealBounds(vectPointsContour,m_Segments,thick,i);
//   float x,y;
//   x = VectPointsBound[0].x(); y = VectPointsBound[0].y();
//   Vector2D P0(x,y);
//   x = VectPointsBound[1].x(); y = VectPointsBound[1].y();
//   Vector2D P1(x,y);
//   x = VectPointsBound[2].x(); y = VectPointsBound[2].y();
//   Vector2D P2(x,y);
//   x = VectPointsBound[3].x(); y = VectPointsBound[3].y();
//   Vector2D P3(x,y);
//   int a1 = m_Segments[i].m_a;
//   int b1 = m_Segments[i].m_b;
//   int mu1 = m_Segments[i].m_mu;
//   int omega1 = m_Segments[i].m_omega;
//   int mu_sup = mu1 + omega1;
//   float aD1 = (float)( P0.y() - P1.y() )/ ( P0.x() - P1.x() );
//   float bD1 = P0.y() - aD1 * P0.x();
//   cerr << "D1(" << aD1 << "," << bD1 << ")" << endl;
//   float aD2 = (float)( P3.y() - P2.y() )/ ( P3.x() - P2.x() );
//   float bD2 = P3.y() - aD2 * P3.x();
//   cerr << "D2(" << aD2 << "," << bD2 << ")" << endl;
//   int first = m_Segments[i].m_first;
//   int last = m_Segments[i].m_last;
//   float maxdistance = 0.0, r1, r2, r;
// //  cout << "i = " << i << "  first = " << first << "  last = " << last << endl;
//   Vector2D fp ( vectPointsContour[m_Segments[i].m_first].x() , vectPointsContour[m_Segments[i].m_first].y() );
//   Vector2D lp ( vectPointsContour[m_Segments[i].m_last].x() , vectPointsContour[m_Segments[i].m_last].y() );
//   int nb = vectPointsContour.size();
// /*  if( last - first < 0 ) last = last + first;
//   for( int p = first+1 ; p <= last ; ++p)
//   {
//     //cout << "p = " << p << "  first = " << first << "  last = " << last << endl;
//     Vector2i nextpoint ( vectPointsContour[p%nb].x() , vectPointsContour[p%nb].y() );
// //    cout << "Nextpoint = (" << nextpoint.x() << "," << nextpoint.y() << ")" << endl;
//     r1 = aD1 * nextpoint.x() + bD1;
// //    cout << "r1 = " << r1 << "   nextpoint.y() = " << nextpoint.y() << endl;
    
//    // if( bD1 > bD2 ? nextpoint.y() < r1 : nextpoint.y() > r1 ) { break;}//if( p > first+(last-first)/3 ) break;}
//    // else
//    // {
//       float distance = (float) ( abs( aD1 * nextpoint.x() - nextpoint.y() + bD1 ) ) / ( sqrt( aD1 * aD1 + 1 ) );
//       //cout << "distance = " << distance << endl;
//       bool b = bD1 > bD2 ? true : false;
//       //cerr << "(bD1,bD2,bD1>bD2)=(" << bD1 << ";" << bD2 << ";" << b << ")" << endl;
//       if( bD1 < bD2 ? nextpoint.y() > r1 : nextpoint.y() < r1  && maxdistance <= distance )
//       { 
//         maxdistance = distance;
//         fp = nextpoint;
//   //      cout << "Correctfp = (" << fp.x() << "," << fp.y() << ")" << endl;
//       }
//    // }
//   }
//   maxdistance = 0.0;
//   lp = lp;
//   for( int p = last-1 ; p >= first ; --p)
//   {
// //    cout << "p = " << p << endl;
//     if( p < 0 ) p = p + nb;
//     Vector2i previouspoint ( vectPointsContour[p%nb].x() , vectPointsContour[p%nb].y() );
//     //cout << "Previouspoint = (" << previouspoint.x() << "," << previouspoint.y() << ")" << endl;
//     r2 = aD2 * previouspoint.x() + bD2;
//     //cerr << "   r2 = " << r2 << "   Previouspoint.y() = " << previouspoint.y() << endl;
//     //if( bD1 > bD2 ? previouspoint.y() > r2 : previouspoint.y() < r2 ) { break;}//if( p < last-(last-first)/3 ) break;}
//     //else
//     //{
//       float distance = (float) ( abs( aD2 * previouspoint.x() - previouspoint.y() + bD2 ) ) / ( sqrt( aD2 * aD2 + 1 ) );
// //      cout << "distance = " << distance << endl;
//       if( bD1 < bD2 ? previouspoint.y() < r2 : previouspoint.y() > r2 && maxdistance <= distance )
//       { 
//         maxdistance = distance;
//         lp = previouspoint;
// //        cout << "Correctlp = (" << lp.x() << "," << lp.y() << ")" << endl;
//       }
//    // }
//   }*/
//     BlurredSegmentTgtCover tgc;
//     std::deque<int> D = tgc.VectDeque[i];
//     for( int p = 0 ; p < nb ; ++p)
//     {
//       Vector2D nextpoint ( vectPointsContour[ D[p] ].x() , vectPointsContour[ D[p] ].y() );
//       cout << "Nextpoint = (" << nextpoint.x() << "," << nextpoint.y() << ")" << endl;
//       r1 = aD1 * nextpoint.x() + bD1;
//       float distance = (float) ( abs( aD1 * nextpoint.x() - nextpoint.y() + bD1 ) ) / ( sqrt( aD1 * aD1 + 1 ) );
//       if( bD1 < bD2 ? nextpoint.y() > r1 : nextpoint.y() < r1  && maxdistance <= distance )
//       { 
//         maxdistance = distance;
//         fp = nextpoint;
//       }
//     }
//     maxdistance = 0.0;
//     lp = lp;
//     for( int p = 0 ; p < nb ; ++p)
//     {
//        Vector2D previouspoint ( vectPointsContour[ D[p] ].x() , vectPointsContour[ D[p] ].y() );
//        r2 = aD2 * previouspoint.x() + bD2;
//        float distance = (float) ( abs( aD2 * previouspoint.x() - previouspoint.y() + bD2 ) ) / ( sqrt( aD2 * aD2 + 1 ) );
//        if( bD1 < bD2 ? previouspoint.y() < r2 : previouspoint.y() > r2 && maxdistance <= distance )
//        { 
//           maxdistance = distance;
//           lp = previouspoint;
//        }
//     }
//   /////////////////////////////////////////////////////////////////////////////////////
//   std::vector<Vector2D> vectPointsBound;
//   float a2 = - b1;
//   float b2 = a1;
//   float mfp = a2 * fp.x() - b2 * fp.y();
//   float mlp = a2 * lp.x() - b2 * lp.y();
//   //cout << "a2 = " << a2 << "  b2 = " << b2 << "   mfp = " << mfp << "   mlp = " << mlp << endl;
//   // coordinate of P0.
//   cout << fixed << setprecision (2);
//   x = (float)( mu1 * b2 - b1 * mfp ) / ( a1 * b2 - a2 * b1 ); 
//   y = (float)( b1 != 0 ? ( a1 * x - mu1 ) / b1 : fp.y() );
//   //cout << "x0 = " << x << "  y0 = " << y << endl; 
//   Vector2D NewP0( x , y );
//   vectPointsBound.push_back( NewP0 );
//   // coordinate of P1.
//   x = (float)( mu_sup * b2 - b1 * mfp ) / ( a1 * b2 - a2 * b1 ); 
//   y = (float)( b1 != 0 ? ( a1 * x - mu_sup ) / b1 : fp.y() );
//   Vector2D NewP1( x , y );
//   vectPointsBound.push_back( NewP1 );
//   // coordinate of P2.
//   x = (float)( mu_sup * b2 - b1 * mlp ) / ( a1 * b2 - a2 * b1 ); 
//   y = (float)( b1 != 0 ? ( a1 * x - mu_sup ) / b1 : lp.y() );
//   Vector2D NewP2( x , y );
//   vectPointsBound.push_back( NewP2 );
//   // coordinate of P3.
//   x = (float)( mu1 * b2 - b1 * mlp ) / ( a1 * b2 - a2 * b1 ); 
//   y = (float)( b1 != 0 ? ( a1 * x - mu1 ) / b1 : lp.y() );
//   Vector2D NewP3( x , y );
//   vectPointsBound.push_back( NewP3 );
//   cerr << "NewP0(" << NewP0.x() << "," << NewP0.y() << ")" << endl;
//   cerr << "NewP1(" << NewP1.x() << "," << NewP1.y() << ")" << endl;
//   cerr << "NewP2(" << NewP2.x() << "," << NewP2.y() << ")" << endl;
//   cerr << "NewP3(" << NewP3.x() << "," << NewP3.y() << ")" << endl;
//   /////////////////////////////////////////////////////////////////////////////////////
//   return vectPointsBound;
// }
// //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
// void
// transformFreemanChain(FreemanChain & fc, vector<uint> & c2trans , const FreemanChain &fcSrc,  int samplingSize, 
// 		      int xIni, int yIni){
//   FreemanChain c2;
//   vector<uint> c2subc;
//   vector<uint> subc2c;
//   FreemanChain::subsample( c2, c2subc, subc2c, fcSrc, 
// 			   samplingSize, samplingSize,
// 			   xIni, yIni);
  
//   int tailleChaine = fcSrc.chain.size();
//   c2trans.clear();

//   for(int k=0; k<tailleChaine; k++){
//     c2trans.push_back(c2subc.at(k));
//   }

//   for(int k=0; k<tailleChaine; k++){
//     c2trans.push_back(c2subc.at(k));
//   }

//   for(int i=0; i<nbIterationSpikes; i++){
//     FreemanChain c3;
//     FreemanChain c4;
//     vector<uint> o2i;
//     vector<uint> i2o;     
//     vector<uint> o2i2;
//     vector<uint> i2o2;     
//     if ( ! FreemanChain::cleanOuterSpikes( c3, o2i, i2o, c2, true  ) )
//       cerr << "Contour with no interior !" << endl;
//     if ( ! FreemanChain::cleanOuterSpikes( c4, o2i2, i2o2, c3, false  ) )
//     cerr << "Contour with no interior !" << endl;
//     for(int k=0; k<tailleChaine; k++){
//       int newIndex = c2trans.at(k);      
//       newIndex = o2i.at(newIndex%o2i.size());
//       newIndex = o2i2.at(newIndex%o2i2.size());
//       c2trans.at(k)=newIndex%tailleChaine;
//     }
    
//     c2.chain=c4.chain;    
//   }
//   fc.chain= c2.chain;  
// }

// ////////////////////////////////////////////////////////////////////////////////////////////////
// ////////////////////////////////////////////////////////////////////////////////////////////////

// class XFIGFilter {

// private:
//   ostream & m_out;

// public:
//   const int RESOLUTION;
//   const double EPS_MAGNIFICATION;
//   const double ONE_IN_CM;
//   int h;
//   int v;
//   double xmin;
//   double xmax;
//   double ymin;
//   double ymax;
//   int color_number;

//   XFIGFilter( ostream & out, int res = 1200, double eps_mag = 100.0, 
// 	      double one_in_cm = 0.5, int h = 1, int v = 1 )
//     : m_out( out ), RESOLUTION( res ), EPS_MAGNIFICATION( eps_mag ),
//       ONE_IN_CM( one_in_cm ), h (h ), v (v )
//   {
//     xmin = 0.0;
//     ymin = 0.0;
//     xmax = 29.7;
//     ymax = 21.0;
//     color_number = 32;
//   }

//   void setView( double x1, double y1, double x2, double y2 )
//   {
//     xmin = x1;
//     xmax = x2;
//     ymin = y1;
//     ymax = y2;
//   }

//   void outputHeader()
//   {
//     m_out << "#FIG 3.2" << endl
// 	  << "Landscape" << endl
// 	  << "Center" << endl
// 	  << "Inches" << endl
// 	  << "A4" << endl
// 	  << EPS_MAGNIFICATION << endl
// 	  << "Single" << endl
// 	  << "-1" << endl 
//       // color number for transparent color for GIF
//       // export. -3=background, -2=None, -1=Default, 0-31 for standard
//       // colors or 32- for user colors)
// 	  << RESOLUTION <<" 2" << endl;
//   }

//   int toX( double x )
//   {
//     return (int) floor( RESOLUTION * x * ONE_IN_CM * h / 2.54 );
//   }

//   int toY( double y )
//   {
//     return (int) floor( RESOLUTION * y * ONE_IN_CM * v / 2.54 );
//   }

//   void fitInView( deque< pair<int,int> > & inside_idx,
// 		  const vector<double> & x,
// 		  const vector<double> & y )
//   {
//     bool inside = false;
//     int first = 0;
//     for ( int i = 0; i < (int) x.size(); ++i )
//       {
// 	if ( ( x[ i ] >= xmin ) && ( x[ i ] <= xmax ) &&
// 	     ( y[ i ] >= ymin ) && ( y[ i ] <= ymax ) )
// 	  {
// 	    if ( ! inside ) 
// 	      {
// 		inside = true;
// 		first = i;
// 	      }
// 	  }
// 	else
// 	  {
// 	    if ( inside )
// 	      {
// 		inside_idx.push_back( make_pair( first, i ) );
// 		inside = false;
// 	      }
// 	  }
//       }
//     if ( inside )
//       inside_idx.push_back( make_pair( first, x.size() ) );
//   }

//  void outputPolyline( const vector<double> & x,
// 		       const vector<double> & y,
// 		       int first,
// 		       int after_last,
// 		       int color = 0,
// 		       int thickness = 1,
// 		       int line_style = 0,
// 		       int depth = 50 )
//   {
//     //param 6 :fill color
//     //param 9: fill style
//     //param 4: epaisseur
//     //param 5: couleur trait
//     m_out << "2 1 " << line_style << " "
// 	  << thickness << " " 
// 	  << color << " "
// 	  << "1" << " " // fill color
// 	  << depth << " " // depth
// 	  << "-1 -1 " // pen_style, area_fill (enumeration type, -1 = no fill)
// 	  << "3.14 " // style_val (1/80 inch, spec for dash/dotted lines)
// 	  << "0 0 " // join_style, cap_style (only used for POLYLINE)
// 	  << "-1 0 0 " // radius (1/80 inch, radius of arc-boxes), forward_arrow (0: off, 1: on), backward_arrow (0: off, 1: on)
// 	  << after_last - first << " " << endl;
//     //	      << x.size() << " " << endl;
//     for ( int i = first; i < after_last; ++i )
//     {
//       //m_out << " x y : " << x[ i ] << " " << y[ i ] << endl;
//       m_out << toX( x[ i ] ) << " " << toY( y[ i ] ) << endl;
//     }
//   }
// /**
//    * @param line_style -1 = Default, 0 = Solid, 1 = Dashed, 2 =
//    * Dotted, 3 = Dash-dotted, 4 = Dash-double-dotted, 5 =
//    * Dash-triple-dotted
//    */
//   void outputPolylineInView( const vector<double> & x,
// 			     const vector<double> & y,
// 			     int color = 0,
// 			     int thickness = 1,
// 			     int line_style = 0,
// 			     int depth = 50 )
//   {
//     //deque< pair<int,int> > inside_idx;
//     //fitInView( inside_idx, x, y );
//     //for ( deque< pair<int,int> >::const_iterator it = inside_idx.begin();
// 	//  it != inside_idx.end();
// 	//  ++it )
//       //{
// 	outputPolyline( x, y, 0, x.size(),
// 			color, thickness, line_style, depth );
//      // }
//   }


// void getCoords( vector<double> & x,
// 		vector<double> & y,
// 		const FreemanChain & c,
// 		int start, 
// 		int end )
// {
//   int pos = 0;
//   FreemanChain::const_iterator it = c.begin();
//   for ( ;
// 	it != c.end();
// 	++it )
//     {
//       if ( ( start <= pos )
// 	   && ( ( end <= 0 ) || ( pos < end ) ) )
// 	{
// 	  Vector2i xy( *it );
// 	  x.push_back( xy.x() );
// 	  y.push_back( xy.y() );
// 	  // cout << xy.x() << " " << xy.y() << endl;
// 	}
//       ++pos;
//     }
//   if ( ( start <= pos )
//        && ( ( end <= 0 ) || ( pos < end ) ) )
//     {
//       Vector2i xy( *it );
//       x.push_back( xy.x() );
//       y.push_back( xy.y() );
//     }
  
// }


// void outputLine( double x1, double y1, double x2, double y2,
// 		   int color = 0,
// 		   int thickness = 1,
// 		   int line_style = 0,
// 		   int depth = 50 )
//   {
//     vector<double> x( 2 );
//     vector<double> y( 2 );
//     x[ 0 ] = x1;
//     y[ 0 ] = y1;
//     x[ 1 ] = x2;
//     y[ 1 ] = y2;
//     outputPolyline( x, y, 0, 2, color, thickness, line_style, depth );
//   }
 
// void outputGrid( double x1, double y1, double x2, double y2,
// 		   uint div_x, uint div_y,
// 		   int color = 0,
// 		   int thickness = 1,
// 		   int line_style = 0,
// 		   int depth = 60 )
//   {
//     double y1p = y1 > ymin ? y1 : ymin;
//     double y2p = y2 < ymax ? y2 : ymax;
//     if ( y1p <= y2p )
//       for ( uint i = 0; i <= div_x; ++i )
// 	{
// 	  double x = x1 + ((x2 -x1)*i)/(double)div_x + 0.5;
// 	  if ( ( x >= xmin ) && ( x <= xmax ) )
// 	    outputLine( x, y1p, x, y2p, color, thickness, line_style, depth );
// 	}

//     double x1p = x1 > xmin ? x1 : xmin;
//     double x2p = x2 < xmax ? x2 : xmax;
//     if ( x1p <= x2p )
//       for ( uint j = 0; j <= div_y; ++j )
// 	{
// 	  double y = y1 + ((y2 -y1)*j)/(double)div_y + 0.5;
// 	  if ( ( y >= ymin ) && ( y <= ymax ) )
// 	    outputLine( x1p, y, x2p, y, color, thickness, line_style, depth );
// 	}
//   }


// void DisplayContourXFIG( FreemanChain & fc, 
//                          int color, int thickness, 
//                          bool filled, bool closeCnt, int zoom )
// {
//   //param 6 :fill color
//   //param 9: fill style
//   //param 4: epaisseur
//   //param 5: couleur trait
//   cout <<setiosflags(ios_base::fixed);
//   cout<< setprecision(0);
  
//   cout << "\n 2 " << ((closeCnt)? 1 : 0 )<< " 0 "  << thickness << " " << color  
//        << "  10 19 -1 "<<((filled)? 0:-1)  <<" 0.000 0 0 -1 0 0 " << fc.chain.size()+((closeCnt)? 1 : 0 ) 
//        << " "  <<endl ;
//   for (FreemanChain::const_iterator it = fc.begin();it != fc.end(); ++it )
//    cout << ((*it).x())*RESOLUTION*zoom << " " << ((*it).y())*RESOLUTION*zoom << " " ;        
//   if(closeCnt)
//    cout << (((*(fc.begin())).x()))*RESOLUTION*zoom << " " << ((*(fc.begin())).y())*RESOLUTION*zoom << " " ;        
// }


// void displayPixel( int pixelx, int pixely, int color, int colorFill, 
//                    float tailleh, float taillev, int depth )
// {
//   cout <<setiosflags(ios_base::fixed);
//   cout<< setprecision(1);
//   Vector2D p0 (pixelx+tailleh/2.0, pixely+taillev/2.0);
//   Vector2D p1 (pixelx+tailleh/2.0, pixely-taillev/2.0);
//   Vector2D p2 (pixelx-tailleh/2.0, pixely-taillev/2.0);
//   Vector2D p3 (pixelx-tailleh/2.0, pixely+taillev/2.0);
// //  cout << "taille = " << tailleh << endl;
// //  cout << "Pixel(" << tailleh/2 << "," << taillev/2 << ")" << endl;
// //  cout << "P0(" << p0.x() << "," << p0.y() << ")" << endl;
//   cout << setprecision(0) << "\n 2 3 0 1 " << color  << " " << colorFill 
//        << " " << depth<<    " -1 20 0.000 0 0 -1 0 0 5 " <<endl ;  
  
//   cout << p0.x()*RESOLUTION << " " << p0.y()*RESOLUTION << " " << p1.x()*RESOLUTION 
//        << " " << p1.y()*RESOLUTION << " " << p2.x()*RESOLUTION << " " 
//        << p2.y()*RESOLUTION << " "<< p3.x()*RESOLUTION << " " << p3.y()*RESOLUTION 
//        << " "  << p0.x()*RESOLUTION << " " << p0.y()*RESOLUTION << endl;
// }


// void displayContourPixelsXFIG( FreemanChain & fc, 
//                                int color, int colorfill, 
//                                int depth, int h, int v )
// {
//   FreemanChain::const_iterator it = fc.begin();
//   int nb = fc.chain.size();
//   float valx = fc.x0;
//   float valy = fc.y0; 
//   displayPixel(valx, valy, color, colorfill, h, v, depth);
//   for( int i = 1 ; i < nb ; ++i )
//   {  
//     if( fc.code( i ) == 0 ) valx += h;
//     else if( fc.code( i ) == 1 ) valy -= v;
//     else if( fc.code( i ) == 2 ) valx -= h;
//     else valy += v;
//     displayPixel(valx, valy, color, colorfill, h, v, depth);
//     it.nextInLoop();
    
//   }
// }

// void displayCroix( float pixelx, float pixely, int color, 
//                    int width, int heigth, int depth )
// { 
//   cout <<setiosflags(ios_base::fixed);
//   cout<< setprecision(0);
//   cout<< "\n 2 1 0 5 " << color << " " << color << " " << depth << " -1 -1 0.000 0 0 -1 0 0 2"<< endl;
//   cout << (pixelx-width/2.0)*RESOLUTION << " " << (pixely-heigth/2.0)*RESOLUTION << " " 
//        << (pixelx+width/2.0)*RESOLUTION << " " << (pixely+heigth/2.0)*RESOLUTION << endl; 
//   cout<< "\n 2 1 0 5 " << color << " " << color << " " << depth << " -1 -1 0.000 0 0 -1 0 0 2"<< endl;
//   cout << (pixelx+width/2.0)*RESOLUTION << " " << (pixely-heigth/2.0)*RESOLUTION << " " 
//        << (pixelx-width/2.0)*RESOLUTION << " " << (pixely+heigth/2.0)*RESOLUTION << endl;   
// }

// void drawLine( ostream & os, 
//                const Vector2D & point1, 
//                const Vector2D & point2,
//                int color, 
//                int linewidth, 
//                int depth )
// { 
//    os <<setiosflags(ios_base::fixed);
//    os << setprecision(0);
//    os << "\n 2 1 0 " << linewidth << " " <<color << " 7 " << depth << " -1 -1 0.000 0 0 -1 0 0 2"
//       << endl;
//    os << (point1.x())* RESOLUTION << " " << (point1.y())* RESOLUTION << " " 
//       << (point2.x())* RESOLUTION << " " << (point2.y())* RESOLUTION 
//       << endl; 
// }

// /**
//  *  Display the tangential cover of a contour between i and j whose i and j 
//  *  are the indexes of the tangential cover.
//  *  @param vectPointsContour contains all the points of a contour.
//  *  @param thick the diagonal thickness of a blurred segment.
//  */
// void DisplayTangentialCover( ostream & os, 
//                              const vector<Vector2D> & vectPointsContour, 
//                              const vector<BlurredSegmentTgtCover::Segment> & m_Segments, 
//                              int color, 
//                              int depth,
//                              int linewidth, 
//                              int thick,
//                              int i,
//                              int j )
// {
//   int nbSegment = m_Segments.size();
//   //cout << "Nb of segment = " << nbSegment << endl;
//   vector<Vector2D> vectPointsTC;
//   if( i >= nbSegment )
//      i = nbSegment - 1;
//   if( j >= nbSegment )
//      j = nbSegment - 1;
//   for( int k = i ; k <= j ; ++k )
//   {
//     vectPointsTC = ComputeBasicBounds( vectPointsContour, m_Segments , thick , k );
//     Vector2D p0 ( vectPointsTC[0].x() , vectPointsTC[0].y() );
//     Vector2D p1 ( vectPointsTC[1].x() , vectPointsTC[1].y() );
//     Vector2D p2 ( vectPointsTC[2].x() , vectPointsTC[2].y() );
//     Vector2D p3 ( vectPointsTC[3].x() , vectPointsTC[3].y() );
//     cerr << "k = " << k 
//          << "   p0(" << p0.x() << "," << p0.y() << ")"  
//          << "   p1(" << p1.x() << "," << p1.y() << ")" 
//          << "   p2(" << p2.x() << "," << p2.y() << ")" 
//          << "   p3(" << p3.x() << "," << p3.y() << ")" << endl;
//     cerr << "==================================================" << endl;
//     drawLine( cout, p0, p1, color, linewidth, depth);
//     drawLine( cout, p1, p2, color, linewidth, depth);
//     drawLine( cout, p2, p3, color, linewidth, depth);
//     drawLine( cout, p3, p0, color, linewidth, depth);
//   }
// }

// };

// ///////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////
// //
// // M A I N
// //
// ///////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////

// int
// main( int argc, char** argv ) 
// {
//   // -------------------------------------------------------------------------
//   // Prepare arguments.
//   //ShapeHelper::addStarShapedArgs( args );
//   args.addBooleanOption( "-header", "-header: add XFIG header." );
//   args.addBooleanOption( "-auto_center", "-auto_center: automatically centers the shape in the PGM image." );
//   args.addOption( "-NoiseLevelForEachPixelWithIntegerThick", "-NoiseLevelForEachPixelWithIntegerThick <thick> <min_width> <max_slope> <min_slope> : ","10","1","-1","0" );
//   args.addOption( "-NoiseLevelForEachPixelWithNonIntegerThick", "-NoiseLevelForEachPixelWithNonIntegerThick <thick> <min_width> <max_slope> <min_slope> : ","10","1","-1","0" );
//   args.addOption( "-DistributionOfSlopesWithIntegerThick", "-DistributionOfSlopesWithIntegerThick <thick> <flatside> : Distribution of the slopes of the multi-thickness lengths of blurred segments for different kinds of shapes generated with several integer thickness","10","0" );
//   args.addOption( "-DistributionOfSlopesWithNonIntegerThick", "-DistributionOfSlopesWithNonIntegerThick <thick> <flatside> : Distribution of the slopes of the multi-thickness lengths of blurred segments for different kinds of shapes generated with several non integer thickness","10","0" );
//   args.addOption( "-unitcm", "-unitcm <u>: tells the length in cm of a unit step.", "0.5" );
//   args.addOption( "-subsampling", "-subsampling <h> <v>: the discretization step or resolution, the closer to 0, the finer.", "1", "1" );
//   args.addOption( "-NumberTCofEachPoint", "-NumberTCofEachPoint <thick>: computes the numbers of blurred segments of each point in the contour", "1" );
//   args.addOption( "-TimeExtractSegment", "-TimeExtractSegment <thick>: the diagonal thickness of the blurred segment.", "1" );
//   args.addOption( "-connexity", "-connexity <connexity>: connexity of shape.", "0" );
//   args.addOption( "-blurred", "-blurred <blurred>: 0 <-> regular segment, 1 <-> blur segments.", "0" );
//   args.addOption( "-DisplayContourXFIG", "-DisplayContourXFIG <color> <depth>: display a contour XFIG.", "1", "1" );
//   args.addOption( "-DisplayContourXFIGPixels", "-DisplayContourXFIGPixels <color> <depth>: display a pixels contour XFIG.", "1","1");
//   args.addOption( "-DisplayTangentialCover", "-DisplayTangentialCover <first_TC> <last_TC> <depth>: display the Tangential Cover of a contour with thickness 'thick'.", "0", "0", "1");
//   args.addOption( "-AverageLengthOfMSForPointWithIntegerThick", "-AverageLengthOfMSForPointWithIntegerThick <thick> <indexpoint> : computes the average length of maximal segments for a point indexpoint at Integer thickness thick","10", "0" );
//   args.addOption( "-AverageLengthOfMSForPointWithNonIntegerThick", "-AverageLengthOfMSForPointWithNonIntegerThick <thick> <indexpoint> : computes the average length of maximal segments for a point indexpoint at non integer thickness thick","10", "0" );
//   args.addOption( "-DisplayNoiseLevel", "-DisplayNoiseLevel <color> <depth> : display a number inside each pixel to specify its Noise Levels.", "1","1");
//   args.addOption( "-DisplayPoint", "-DisplayPoint<p1><p2><p3><color> show three points like crosses.","0","1","3","4" );
//   args.addOption( "-ShowNoiseDetection", "-ShowNoiseDetection <color><fillcolor><depth> : display a pixel of size equal to the noise level.", "1","1","1");
// // -------------------------------------------------------------------------
//   // Read Freeman chain and creates space/contour.
//   Clock::startClock();
//   if ( ( argc <= 1 ) 
//        || ! args.readArguments( argc, argv ) ) 
//     {
//       cerr << args.usage( "test_TCover", 
// 			  "Tests the tangential cover of a digital curve based on the blurred segment with different thicknesses.",
// 			  "" )      
// 	   << endl;
//       return 1;
//     }  
//     FreemanChain c;
//   FreemanChain::read( cin, c );

//   if ( ! cin.good() )
//     {
//       cerr << "Error reading Freeman chain code." << endl; 
//       return 2;
//     }
  
//   cerr<< "here" <<endl;  
//   KnSpace* ks;
//   if ( args.check( "-auto_center" ) )
//     ks = ShapeHelper::makeSpaceFromFreemanChain( c );
//   else 
//   {
//      // Build space.
//      uint d = StandardArguments::dim( args );
//      if ( d != 2 )
//      {
// 	cerr << "Dimension is 2." << endl;
// 	return 2;
//      }
//      Kn_size sizes[ d ];
//      StandardArguments::fillSizes( args, sizes );
//      ks = new KnSpace( 2, sizes );
//   }

//   cerr << "here 3 "<< endl;


  
//   int h,v;
//   FreemanChain fc;
//   vector<uint> c2trans;
//   if ( args.check( "-subsampling" ) )
//   {
//     h= args.getOption( "-subsampling" )->getIntValue( 0 );
//     v= args.getOption( "-subsampling" )->getIntValue( 1 ); 
//     transformFreemanChain( fc, c2trans , c,  h, 0, 0 );
//     cerr << "FreemanChain = " << fc.chain << endl;               
//   }
  

//   cerr<< "here2" <<endl;
  
//   int connexity;  // 0 <-> 4-connex
//   if ( args.check( "-connexity" ) )
//   {
//     connexity = args.getOption( "-connexity" ) -> getIntValue( 0 );  
//   }

//   int blur;  // 0 <-> 4-connex
//   if ( args.check( "-blurred" ) )
//   {
//     blur = args.getOption( "-blurred" ) -> getIntValue( 0 );  
//   }

//   BlurredSegmentTgtCover tgc, tgc1, tgc2;
//   vector<Vector2D> vectPointsContour, vectPointsContour1;
//   std::vector<Vector2D> pointsbound;
//   float thick;

//   if ( args.check( "-TimeExtractSegment" ) )
//   {
//     float thickmax = args.getOption( "-TimeExtractSegment" ) -> getIntValue( 0 );
//     BlurredSegmentTgtCover tgc;
//     tgc.init( c );
//     //vectPointsContour = tgc.getPointsContour();
//     tgc.ExtractSegment( thickmax, connexity, blur );
//     //vectm_Segment_thick.push_back( tgc.getSegmentContour() );
//     //vectvectLengthMS.push_back( tgc.getLengthMS( vectPointsContour) );
//   }
//   long t = Clock::stopClock();
//   cerr << "# MultiThicknesses in " << t << " ms." << endl;  

//   if ( args.check( "-NumberTCofEachPoint" ) )
//   {
//     thick = args.getOption( "-NumberTCofEachPoint" ) -> getIntValue( 0 );
//     tgc.init( c );
//     vectPointsContour = tgc.getPointsContour();
//     //cout << "Number of points = " << vectPointsContour.size() << endl;
//     tgc.init();
//     tgc.ExtractSegment( thick, connexity, blur );
//     //tgc.PrintTgCover();
//     pointsbound = tgc.getPointsBoundTC();
//     std::vector<BlurredSegmentTgtCover::Segment> m_Segment_thick;
//     m_Segment_thick = tgc.getSegmentContour();
//     tgc.NumberTCofEachPoint( vectPointsContour, m_Segment_thick );
//   }

//   std::vector<float> NoiseLevel;
//   tgc.init( c );
//   vectPointsContour = tgc.getPointsContour();
//   int nb = vectPointsContour.size();
//   cerr << "Number of points = " << nb << endl;
//   std::vector<float> VectSlope( nb ); 

//   if ( args.check( "-NoiseLevelForEachPixelWithIntegerThick" ) )
//   {
//     int thickmax = args.getOption( "-NoiseLevelForEachPixelWithIntegerThick" )-> getIntValue( 0 );
//     int min_width = args.getOption( "-NoiseLevelForEachPixelWithIntegerThick" )-> getIntValue( 1 );
//     double max_slope = args.getOption( "-NoiseLevelForEachPixelWithIntegerThick" )-> getDoubleValue( 2 );
//     double min_slope = args.getOption( "-NoiseLevelForEachPixelWithIntegerThick" )-> getDoubleValue( 3 );
//     std::vector< vector<BlurredSegmentTgtCover::Segment> > vectm_Segment_thick;
//     vector< vector<float> > vectvectLengthMS;
//     for( uint i = 1 ; i <= thickmax ; ++i ) 
//     {
//       BlurredSegmentTgtCover tgc;
//       tgc.init( c );
//       vectPointsContour = tgc.getPointsContour();
//       tgc.ExtractSegment( i, connexity, blur );
//       vectm_Segment_thick.push_back( tgc.getSegmentContour() );
//       vectvectLengthMS.push_back( tgc.getLengthMS( vectPointsContour) );
//     }      
//     for( int j = 0 ; j < vectPointsContour.size() ; ++j ) 
//     {
//       int nb = 0; 
//       ScaleProfile sp;
//       sp.init( thickmax );
//       //sp.initNonInteger( thickmax );
//       vector<float> vectAverageLengthOfMS;
//       for( uint i = 1 ; i <= thickmax ; ++i ) 
//       {
//         BlurredSegmentTgtCover tgc1;
//         vector<BlurredSegmentTgtCover::Segment> m_Segment_thick = vectm_Segment_thick[nb];
//         vector<float> vectLengthMS = vectvectLengthMS[nb];
//         float AverageLengthOfMS = tgc1.getAverageLengthOfMS( j, m_Segment_thick, vectPointsContour, vectLengthMS );
//         //cout << "j = " << j << "  k = " << k << "  AverageLengthOfMS = " << AverageLengthOfMS << endl;
//         vectAverageLengthOfMS.push_back( AverageLengthOfMS );
//         sp.addValue( nb, AverageLengthOfMS );
//         ++nb;
//        }

//       //BK Change 18/02/2011
//       //VectSlope[j] = sp.SlopeRegressionLinear( vectAverageLengthOfMS, thickmax );
//       VectSlope[j] = sp.slope();// vectAverageLengthOfMS, thickmax );
//        // cout << "j = " << j << "   VectSlope[j] = " << VectSlope[j] << endl;
//        std::vector<double> X;
//        std::vector<double> Y;
//        sp.getProfile(X,Y);
//        float noiselevel = sp.noiseLevel( min_width, max_slope, min_slope );
//        NoiseLevel.push_back( noiselevel );
//        //cout << "j = " << j << " Noise Level =  " << noiselevel << endl;
//     }
//   } 

//   if ( args.check( "-NoiseLevelForEachPixelWithNonIntegerThick" ) )
//   {
//     int thickmax = args.getOption( "-NoiseLevelForEachPixelWithNonIntegerThick" )-> getIntValue( 0 );
//     int min_width = args.getOption( "-NoiseLevelForEachPixelWithNonIntegerThick" )-> getIntValue( 1 );
//     double max_slope = args.getOption( "-NoiseLevelForEachPixelWithNonIntegerThick" )-> getDoubleValue( 2 );
//     double min_slope = args.getOption( "-NoiseLevelForEachPixelWithNonIntegerThick" )-> getDoubleValue( 3 );
//     std::vector< vector<BlurredSegmentTgtCover::Segment> > vectm_Segment_thick;
//     vector< vector<float> > vectvectLengthMS;
//     for( uint i = 1 ; i <= thickmax ; ++i ) 
//     {
//       float k = i;
//       for( int l = 1 ; l <= 5 ; ++l )
//       {       
//         BlurredSegmentTgtCover tgc;
//         tgc.init( c );
//         vectPointsContour = tgc.getPointsContour();
//         tgc.ExtractSegment( k, connexity, blur );
//         vectm_Segment_thick.push_back( tgc.getSegmentContour() );
//         vectvectLengthMS.push_back( tgc.getLengthMS( vectPointsContour) );
//         k += 0.2;
//       }
//     }      
//     for( int j = 0 ; j < vectPointsContour.size() ; ++j ) 
//     {
//       int nb = 0; 
//       ScaleProfile sp;
//       sp.init(thickmax);//initNonInteger( thickmax );
//       vector<float> vectAverageLengthOfMS;
//       for( uint i = 1 ; i <= thickmax ; ++i ) 
//       {
//          float k = i;
//          for( int l = 1 ; l <= 5 ; ++l )
//          {       
//            BlurredSegmentTgtCover tgc1;
//            vector<BlurredSegmentTgtCover::Segment> m_Segment_thick = vectm_Segment_thick[nb];
//            vector<float> vectLengthMS = vectvectLengthMS[nb];
//            float AverageLengthOfMS = tgc1.getAverageLengthOfMS( j, m_Segment_thick, vectPointsContour, vectLengthMS );
//            //cout << "j = " << j << "  k = " << k << "  AverageLengthOfMS = " << AverageLengthOfMS << endl;
//            vectAverageLengthOfMS.push_back( AverageLengthOfMS );
//            sp.addValue( nb, AverageLengthOfMS );
//            ++nb;
//            k += 0.2;
//          }
//        }
//       //BK Change 18/02/2011
//       //VectSlope[j] = sp.SlopeRegressionLinear( vectAverageLengthOfMS, thickmax );
//        //cout << "j = " << j << "   VectSlope[j] = " << VectSlope[j] << endl;
//        std::vector<double> X;
//        std::vector<double> Y;
//        sp.getProfile(X,Y);
//        float noiselevel = sp.noiseLevel( min_width, max_slope, min_slope );
//        NoiseLevel.push_back( noiselevel );
//        //cout << "j = " << j << " Noise Level =  " << noiselevel << endl;
//      }
//    } 


//   if ( args.check( "-DistributionOfSlopesWithIntegerThick" ) )
//   {
//     int thickmax = args.getOption( "-DistributionOfSlopesWithIntegerThick" )-> getIntValue( 0 );
//     int flatside = args.getOption( "-DistributionOfSlopesWithIntegerThick" )-> getIntValue( 1 );
//     vector< float > vectSlope( nb );
//     ScaleProfile sp;
//     for( int i = 0 ; i < nb ; ++i ) // nb is the number of pixels in the contour.
//     {
//       cerr << "i = " << i << endl;
//       vector<float> vectAverageLengthOfMS;
//       for( int j = 1 ; j <= thickmax ; ++j )
//       {
//           cerr << "j = " << j << endl;
//           BlurredSegmentTgtCover tgc;
//           tgc.init( c );
//           vectPointsContour = tgc.getPointsContour();
//           tgc.ExtractSegment( j, connexity, blur );  
//           std::vector<BlurredSegmentTgtCover::Segment> m_Segment_thick;
//           m_Segment_thick = tgc.getSegmentContour();
//           vector<float> vectLengthMS = tgc.getLengthMS( vectPointsContour );
//           float AverageLengthOfMS = tgc.getAverageLengthOfMS( i, m_Segment_thick, vectPointsContour, vectLengthMS );
//           //cout << "j = " << j << "   AverageLengthOfMS = " << AverageLengthOfMS << endl;
//           vectAverageLengthOfMS.push_back( AverageLengthOfMS );
//       }
//       //BK Change 18/02/2011
//       //vectSlope[i] = sp.SlopeRegressionLinear( vectAverageLengthOfMS, thickmax );
//       vectSlope[i] = sp.slope();

//       cerr << "Slope = " << vectSlope[i] << endl;
//     }
//     sp.DistributionOfSlopes( vectSlope, nb, thickmax, flatside );
//   } 

//   if ( args.check( "-DistributionOfSlopesWithNonIntegerThick" ) )
//   {
//     int thickmax = args.getOption( "-DistributionOfSlopesWithNonIntegerThick" )-> getIntValue( 0 );
//     int flatside = args.getOption( "-DistributionOfSlopesWithIntegerThick" )-> getIntValue( 1 );
//     vector< float > vectSlope( nb );
//     ScaleProfile sp;
//     for( int i = 0 ; i < nb ; ++i ) // nb is the number of pixels in the contour.
//     {
//       cerr << "i = " << i << endl;
//       vector<float> vectAverageLengthOfMS;
//       for( int j = 1 ; j <= thickmax ; ++j )
//       {
//          float k = j;
//          for( int l = 1 ; l <= 5 ; ++l )
//          {
//            //cerr << "j = " << j << endl;
//            BlurredSegmentTgtCover tgc;
//            tgc.init( c );
//            vectPointsContour = tgc.getPointsContour();
//            tgc.ExtractSegment( j, connexity, blur );  
//            std::vector<BlurredSegmentTgtCover::Segment> m_Segment_thick;
//            m_Segment_thick = tgc.getSegmentContour();
//            vector<float> vectLengthMS = tgc.getLengthMS( vectPointsContour );
//            float AverageLengthOfMS = tgc.getAverageLengthOfMS( i, m_Segment_thick, vectPointsContour, vectLengthMS );
//            //cout << "k = " << k << "   AverageLengthOfMS = " << AverageLengthOfMS << endl;
//            vectAverageLengthOfMS.push_back( AverageLengthOfMS );
//            k += 0.2;
//          }
//       }
//       //BK Change 18/02/2011
//       //vectSlope[i] = sp.SlopeRegressionLinear( vectAverageLengthOfMS, thickmax );
//       vectSlope[i] = sp.slope();
      
//       cerr << "Slope = " << vectSlope[i] << endl;
//     }
//     sp.DistributionOfSlopes( vectSlope, nb, thickmax, flatside );
//   } 

//   if ( args.check( "-AverageLengthOfMSForPointWithIntegerThick" ) )
//   {
//     int thickmax = args.getOption( "-AverageLengthOfMSForPointWithIntegerThick" )-> getIntValue( 0 );
//     int indexpoint = args.getOption( "-AverageLengthOfMSForPointWithIntegerThick" )-> getIntValue( 1 );
//     vector<float> vectAverageLengthOfMS;
//     float AverageLength;
//     BlurredSegmentTgtCover tgc;
//     for( int i = 1 ; i <= thickmax ; ++i )
//     {
//         BlurredSegmentTgtCover tgc;
//         tgc.init( c );
//         vectPointsContour = tgc.getPointsContour();
//         tgc.ExtractSegment( i, connexity, blur );
//         std::vector<BlurredSegmentTgtCover::Segment> m_Segment_thick;
//         m_Segment_thick = tgc.getSegmentContour();
//         vector<float> vectLengthMS = tgc.getLengthMS( vectPointsContour );
//         float AverageLengthOfMS = tgc.getAverageLengthOfMS( indexpoint, m_Segment_thick, vectPointsContour, vectLengthMS );
//         vectAverageLengthOfMS.push_back( AverageLengthOfMS );
//         cout << i << "   " << AverageLengthOfMS << endl;
//     }
//   } 

//   if ( args.check( "-AverageLengthOfMSForPointWithNonIntegerThick" ) )
//   {
//     int thickmax = args.getOption( "-AverageLengthOfMSForPointWithNonIntegerThick" )-> getIntValue( 0 );
//     int indexpoint = args.getOption( "-AverageLengthOfMSForPointWithNonIntegerThick" )-> getIntValue( 1 );
//     vector<float> vectAverageLengthOfMS;
//     float AverageLength;
//     BlurredSegmentTgtCover tgc;
//     for( int i = 1 ; i <= thickmax ; ++i )
//     {
//       float k = i;
//       for( int j = 1 ; j <= 5 ; ++j )
//       {
//         BlurredSegmentTgtCover tgc;
//         tgc.init( c );
//         vectPointsContour = tgc.getPointsContour();
//         tgc.ExtractSegment( i, connexity, blur );
//         std::vector<BlurredSegmentTgtCover::Segment> m_Segment_thick;
//         m_Segment_thick = tgc.getSegmentContour();
//         vector<float> vectLengthMS = tgc.getLengthMS( vectPointsContour );
//         float AverageLengthOfMS = tgc.getAverageLengthOfMS( indexpoint, m_Segment_thick, vectPointsContour, vectLengthMS );
//         vectAverageLengthOfMS.push_back( AverageLengthOfMS );
//         cout << k << "   " << AverageLengthOfMS << endl;
//         k += 0.2;
//       }
//     }
//     //float k = 1;
//     //for( int i = 0 ; i < vectAverageLengthOfMS.size()-1 ; ++i )
//     //{
//     //  float slope = ( log( vectAverageLengthOfMS[i+1]/(k+0.2) ) - log( vectAverageLengthOfMS[i]/k ) ) / ( log( k+0.2 ) - log( k ) );
//     //  //cout << k << "   " << vectAverageLengthOfMS[i] << "   " << slope << endl;
//     //  k += 0.2;
//     //}
//   } 

//   XFIGFilter filter( cout, 1200, 10.0, 
// 		     args.getOption( "-unitcm" )->getFloatValue( 0 ), h, v );
  

//   if ( args.check( "-header" ) )
//     //filter.outputHeader();
//     DrawingXFIG::includeXFIGHeader(cout,120,10);
  
//   if ( args.check( "-DisplayContourXFIG" ) )
//   {
//     int color = args.getOption( "-DisplayContourXFIG" )->getIntValue( 0 );
//     int thickness = args.getOption( "-DisplayContourXFIG" )->getIntValue( 1 );
//     DrawingXFIG::drawContour(cout, vectPointsContour, color, color, thickness, true, false, 5, 1, 0, 0);
//   }
 
//   if ( args.check( "-ShowNoiseDetection" ) )
//   {
//     int color = args.getOption( "-ShowNoiseDetection" )->getIntValue( 0 );
//     int fillcolor = args.getOption( "-ShowNoiseDetection" )->getIntValue( 1 );
//     int depth = args.getOption( "-ShowNoiseDetection" )->getIntValue( 2 );
//     FreemanChain::const_iterator it = c.begin();
//     int nb = c.chain.size();
//     float valx = c.x0;
//     float valy = c.y0; 
//     for( int i = 0 ; i < nb ; ++i )
//     {  
//       float size = NoiseLevel[i]-0.2;
//       if( c.code( i ) == 0 ) valx += 1;
//       else if( c.code( i ) == 1 ) valy += 1;
//       else if( c.code( i ) == 2 ) valx -= 1;
//       else valy -= 1;
//       filter.displayPixel(valx, valy, color, fillcolor, size, size, depth);
//     }
//   }
 
//   if ( args.check( "-DisplayNoiseLevel" ) )
//   {
//     int color = args.getOption( "-DisplayNoiseLevel" )->getIntValue( 0 );
//     int depth = args.getOption( "-DisplayNoiseLevel" )->getIntValue( 1 );
//     //BK Change 18/02/2011
//     //DrawingXFIG::drawTextInsidePixel(cout, vectPointsContour, NoiseLevel, 0, 0, 0, 50, color, depth);
//   }
  
//   if ( args.check( "-DisplayContourXFIGPixels" ) )
//   {
//     int color = args.getOption( "-DisplayContourXFIGPixels" )->getIntValue( 0 );
//     int depth = args.getOption( "-DisplayContourXFIGPixels" )->getIntValue( 1 );
//     vector<Vector2i> v; 
//     for(int i=0; i< vectPointsContour.size(); i++){
//       Vector2i v ((int)((vectPointsContour.at(i)).x()),(int)((vectPointsContour.at(i)).y()));
//       //v.push_back(Vector2i();
//     }
//     DrawingXFIG::drawContourPixels(cout, v, color, 1, 0, 0, depth);
//   }

//   if ( args.check( "-DisplayTangentialCover" ) )
//   {
//     int first_TC = args.getOption( "-DisplayTangentialCover" )->getIntValue( 0 );
//     int last_TC = args.getOption( "-DisplayTangentialCover" )->getIntValue( 1 );
//     //int color = args.getOption( "-DisplayTangentialCover" )->getIntValue( 2 );
//     //int linewidth = args.getOption( "-DisplayTangentialCover" )->getIntValue( 3 );
//     int depth = args.getOption( "-DisplayTangentialCover" )->getIntValue( 2 );
//     std::vector<BlurredSegmentTgtCover::Segment> m_Segments;
//     m_Segments = tgc.getSegmentContour();
//     cerr << "Number of TC = " << m_Segments.size() << endl;
//     filter.DisplayTangentialCover( cout, vectPointsContour, m_Segments, 1, depth, 3, thick, first_TC, last_TC );
//   }

//   if ( args.check( "-DisplayPoint" ) )
//   {
//     int p1 = args.getOption( "-DisplayPoint" )->getIntValue( 0 );
//     int p2 = args.getOption( "-DisplayPoint" )->getIntValue( 1 );
//     int p3 = args.getOption( "-DisplayPoint" )->getIntValue( 2 );
//     int color = args.getOption( "-DisplayPoint" )->getIntValue( 3 );
//     DrawingXFIG::drawCross(cout, vectPointsContour.at(p1), color+3, 2, 20, 0 );
//     DrawingXFIG::drawCross(cout, vectPointsContour.at(p2), color+1, 2, 20, 0 );
//     DrawingXFIG::drawCross(cout, vectPointsContour.at(p3), color, 2, 20, 0 );
//   }
  
//   return 0;
// }






///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//
// M A I N
//
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

int
main( int argc, char** argv ) 
{
  


}


