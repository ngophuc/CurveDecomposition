///////////////////////////////////////////////////////////////////////////////
// Reads a freeman contour and a curve given as a set of points, in
// order to sample the curve according to the sampling of the freeman
// contour.
///////////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <string>

#include "ImaGene/base/Arguments.h"
#include "ImaGene/base/StandardArguments.h"
#include "ImaGene/base/Proxy.h"
#include "ImaGene/base/Vector.h"
#include "ImaGene/base/Shapes.h"
#include "ImaGene/dgeometry2d/FreemanChain.h"
#include "ImaGene/dgeometry2d/C4CIteratorOnFreemanChain.h"

using namespace std;
using namespace ImaGene;


static Arguments args;

struct Point2D {
  double x[ 2 ];

  static double sqr( double x )
  {
    return x*x;
  }

  double dist2( const Point2D & other ) const
  {
    return sqrt( sqr( x[ 0 ] - other.x[ 0 ] )
		 + sqr( x[ 1 ] - other.x[ 1 ] ) );
  }
  double norm2() const
  {
    return sqrt( sqr( x[ 0 ] ) + sqr( x[ 1 ] ) );
  }
  Point2D & operator+=( const Point2D & other )
  {
    x[ 0 ] += other.x[ 0 ];
    x[ 1 ] += other.x[ 1 ];
    return *this;
  } 
  Point2D & operator-=( const Point2D & other )
  {
    x[ 0 ] -= other.x[ 0 ];
    x[ 1 ] -= other.x[ 1 ];
    return *this;
  } 
  
  Point2D & operator/=( double d )
  {
    x[ 0 ] /= d;
    x[ 1 ] /= d;
    return *this;
  } 

};

ostream & operator<<( ostream & out, const Point2D & p )
{
  out << p.x[ 0 ] << " " << p.x[ 1 ];
  return out;
}


class SampledCurve {
private:
  vector<Point2D> m_points;

public:
  SampledCurve()
  {
    m_points.clear();
  }
  
  void add( const Point2D & p )
  {
    m_points.push_back( p );
  }
  
  uint nbPoints() const
  {
    return m_points.size();
  }

  uint indexClosest( const Point2D & p ) const
  {
    uint j = 0;
    double min_d = p.dist2( m_points[ j ] );

    for ( uint i = 1; i < m_points.size(); ++i )
      {
	double d = p.dist2( m_points[ i ] );
	if ( d < min_d ) 
	  {
	    j = i;
	    min_d = d;
	  }
      }
    return j;
  }

  const Point2D & x( uint i ) const
  {
    return m_points[ i ];
  }

  Point2D smooth_x( uint i, uint k ) const
  {
    if ( k == 0 )
      return x( i );
    uint ni = ( i + 1 ) % nbPoints();
    uint pi = ( i + nbPoints() - 1 ) % nbPoints();
    Point2D mid( smooth_x( i, k - 1 ) );
    mid += mid;
    Point2D prev( smooth_x( ni, k - 1 ) );
    Point2D next( smooth_x( pi, k - 1 ) );
    mid += prev;
    mid += next;
    mid /= 4.0;
    return mid;
  }

  Point2D dx_di( uint i, uint k ) const
  {
    uint ni = ( i + 1 ) % nbPoints();
    uint pi = ( i + nbPoints() - 1 ) % nbPoints();
    Point2D next( smooth_x( ni, k ) );
    Point2D prev( smooth_x( pi, k ) );
    next -= prev;
    next /= 2.0;
    return next;
  }

  Point2D d2x_di2( uint i, uint k ) const
  {
    uint ni = ( i + 1 ) % nbPoints();
    uint pi = ( i + nbPoints() - 1 ) % nbPoints();
    Point2D next( dx_di( ni, k ) );
    Point2D prev( dx_di( pi, k ) );
    next -= prev;
    next /= 2.0;
    return next;
  }

  double ds( uint i, uint j, uint k ) const
  {
    double d = 0.0;
    Point2D cur( smooth_x( i, k ) );
    while ( i != j )
      {
	i = ( i + 1 ) % nbPoints();
	Point2D next( smooth_x( i, k ) );
	d += cur.dist2( next );
	cur = next;
      }
    return d;
  }

  Point2D dx_ds( uint i, uint k ) const
  {
    uint ni = ( i + 1 ) % nbPoints();
    uint pi = ( i + nbPoints() - 1 ) % nbPoints();
    Point2D next( smooth_x( ni, k ) );
    Point2D prev( smooth_x( pi, k ) );
    double d = ds( pi, ni, k );
    next -= prev;
    next /= d;
    return next;
  }

  Point2D d2x_ds2( uint i, uint k ) const
  {
    uint ni = ( i + 1 ) % nbPoints();
    uint pi = ( i + nbPoints() - 1 ) % nbPoints();
    Point2D next( dx_ds( ni, k ) );
    Point2D prev( dx_ds( pi, k ) );
    double d = ds( pi, ni, k );
    next -= prev;
    next /= d;
    return next;
  }

};

istream & operator>>( istream & in, Point2D & p )
{
  in >> p.x[ 0 ];
  in >> p.x[ 1 ];
}

istream & operator>>( istream & in, SampledCurve & curve )
{
  string str;
  while ( true )
    {
      getline( in, str );
      if ( ! in.good() ) return in;
      if ( ( str.size() > 0 ) && ( str[ 0 ] != '#' ) )
	{
	  istringstream str_in( str );
	  Point2D p;
	  str_in >> p;
	  curve.add( p );
	}
    }
  return in;
}


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
  // -------------------------------------------------------------------------
  // Prepare arguments.
  //StandardArguments::addDigitalArgs( args, 2, false, false );
  StandardArguments::addIOArgs( args, true, true );
  args.addOption( "-curve_file", "-curve_file <filename>: the name of the file listing the coordinates of the points of the curve.", "" );
  args.addOption( "-step", "-step <dh>: the grid step of the digital space.", "1.0" );
  args.addOption( "-origin", "-origin <x0> <y0>: the position of the digitak point (0,0) in the Euclidean plane.", "0.0", "0.0" );
  args.addOption( "-smooth", "-smooth <k>: the smoothing factor, 0 none, 1=[1 2 1], 2=[1 4 6 4 1], etc.", "0" );
  args.addOption( "-orientation", "-orientation <CW/CCW>: the orientation of the contour, CW: interior is to the right, CCW: interior is to the left.", "CW" );
  args.addBooleanOption( "-dG", "-dG: displays the geometry of the curve." );
  args.addBooleanOption( "-dPG", "-dPG: displays the geometry of the curve at pointels of the digital contour." );
  args.addBooleanOption( "-dSG", "-dSG: displays the geometry of the curve at surfels (middle) of the digital contour." );
 
  if ( ( argc <= 0 ) 
       || ! args.readArguments( argc, argv ) ) 
    {
      cerr << args.usage( "sample_curve", 
			  "Reads a digital contour (freeman chain) and a curve given as a set of points, in order to sample the curve according to the sampling of the freeman contour. The freeman contour is given in a file or standard input. The curve is a file composed of two columns of floating points, the coordinates of the points of the curve. The digital contour should have been created in the same frame as the curve. The user may specified a grid step for the digital space, as well as a translation vector.",
			  "" )
	   << endl;
      return 1;
    }

  // -------------------------------------------------------------------------
  // Read Freeman chain and creates space/contour.
  FreemanChain c; 
  istream & in_str = StandardArguments::openInput( args );
  FreemanChain::read( in_str, c );
  if ( ! in_str.good() )
    {
      cerr << "Error reading Freeman chain code." << endl;
      return 2;
    }
  ostream & out = StandardArguments::openOutput( args );
  double dh = args.getOption( "-step" )->getDoubleValue( 0 );
  double x0 = args.getOption( "-origin" )->getDoubleValue( 0 );
  double y0 = args.getOption( "-origin" )->getDoubleValue( 1 );
  int k = args.getOption( "-smooth" )->getIntValue( 0 );
  bool cw = args.getOption( "-orientation" )->getValue( 0 ) == "CW";
  string curve_file =  args.getOption( "-curve_file" )->getValue( 0 );
  if ( curve_file == "" )
    {
      cerr << "Curve filename is the empty string." << endl;
      return 2;
    }

  ifstream in_curve( curve_file.c_str() );
  if ( ! in_curve.good() )
    {
      cerr << "Error reading curve file <" << curve_file << ">" << endl;
      return 2;
    }
  SampledCurve curve;
  in_curve >> curve;

  cerr << "--- curve has " << curve.nbPoints() << " points." << endl;

  if ( args.check( "-dG" ) )
    {
      out << "# -dG: display geometry of curve. (sx,sy) means smoothed version of (x,y)." << endl;
      out << "# k=" << k << " smoothing iterations." << endl;
      out << "# idx x y sx sy dsx/ds dsy/ds d2sx/ds d2sy/ds |kappa|" << endl;
      for ( uint i = 0; i < curve.nbPoints(); ++i )
	{
	  Point2D d1 = curve.dx_ds( i, k );
	  Point2D d2 = curve.d2x_ds2( i, k );
	  double abs_kappa = d2.norm2();
	  double kappa = 
	    ( d1.x[ 0 ] * d2.x[ 1 ] - d1.x[ 1 ] * d2.x[ 0 ] ) >= 0.0 
	    ? ( cw ? -abs_kappa : abs_kappa ) 
	    : ( cw ? abs_kappa : -abs_kappa );
	  out << i << " " 
	      << curve.x( i ) << " " 
	      << curve.smooth_x( i, k ) << " " 
	      << d1 << " " << d2 << " " << kappa << " " 
	      << endl;
	}
    }
  if ( args.check( "-dPG" ) )
    {
      out << "# -dPG: display geometry of curve at pointels. (sx,sy) means smoothed version of (x,y)." << endl;
      out << "# k=" << k << " smoothing iterations." << endl;
      out << "# p_idx dx dy c_idx x y sx sy dsx/ds dsy/ds d2sx/ds d2sy/ds |kappa|" << endl;
      
      uint idx = 0;
      for ( FreemanChain::const_iterator it = c.begin();
	    it != c.end();
	    ++it, ++idx )
	{
	  Point2D p;
	  p.x[ 0 ] = x0 + dh * (*it).x();
	  p.x[ 1 ] = y0 + dh * (*it).y();
	  uint i = curve.indexClosest( p );
	  Point2D d1 = curve.dx_ds( i, k );
	  Point2D d2 = curve.d2x_ds2( i, k );
	  double abs_kappa = d2.norm2();
	  double kappa = 
	    ( d1.x[ 0 ] * d2.x[ 1 ] - d1.x[ 1 ] * d2.x[ 0 ] ) >= 0.0 
	    ? ( cw ? -abs_kappa : abs_kappa ) 
	    : ( cw ? abs_kappa : -abs_kappa );
	  out << idx << " " << p << " "
	      << i << " "
	      << curve.x( i ) << " " 
	      << curve.smooth_x( i, k ) << " " 
	      << d1 << " " << d2 << " " << kappa << " " 
	      << endl;
	}
    }
  if ( args.check( "-dSG" ) )
    {
      out << "# -dSG: display geometry of curve at surfels. (sx,sy) means smoothed version of (x,y)." << endl;
      out << "# k=" << k << " smoothing iterations." << endl;
      int nb_ccw = c.isClosed();
      out << "# nb_ccw=" << nb_ccw << endl;
      out << "# s_idx dx dy c_idx x y sx sy dsx/ds dsy/ds d2sx/ds d2sy/ds |kappa|" << endl;
      
      uint idx = 0;
      FreemanChain::const_iterator it = c.begin();
      Point2D prev;
      prev.x[ 0 ] = x0 + dh * (*it).x();
      prev.x[ 1 ] = y0 + dh * (*it).y();
      ++it;
      for ( ;
	    it != c.end();
	    ++it, ++idx )
	{
	  Point2D p;
	  p.x[ 0 ] = x0 + dh * (*it).x();
	  p.x[ 1 ] = y0 + dh * (*it).y();
	  prev += p;
	  prev /= 2.0;
	  uint i = curve.indexClosest( prev );
	  Point2D d1 = curve.dx_ds( i, k );
	  Point2D d2 = curve.d2x_ds2( i, k );
	  double abs_kappa = d2.norm2();
	  double kappa = 
	    ( d1.x[ 0 ] * d2.x[ 1 ] - d1.x[ 1 ] * d2.x[ 0 ] ) >= 0.0 
	    ? ( cw ? -abs_kappa : abs_kappa ) 
	    : ( cw ? abs_kappa : -abs_kappa );
	  out << idx << " " << prev << " "
	      << i << " "
	      << curve.x( i ) << " " 
	      << curve.smooth_x( i, k ) << " " 
	      << d1 << " " << d2 << " " << kappa << " " 
	      << endl;
	  prev = p;
	}
      it = c.begin();
      Point2D p;
      p.x[ 0 ] = x0 + dh * (*it).x();
      p.x[ 1 ] = y0 + dh * (*it).y();
      prev += p;
      prev /= 2.0;
      uint i = curve.indexClosest( prev );
      Point2D d1 = curve.dx_ds( i, k );
      Point2D d2 = curve.d2x_ds2( i, k );
      double abs_kappa = d2.norm2();
      double kappa = 
	( d1.x[ 0 ] * d2.x[ 1 ] - d1.x[ 1 ] * d2.x[ 0 ] ) >= 0.0 
	? ( cw ? -abs_kappa : abs_kappa ) 
	: ( cw ? abs_kappa : -abs_kappa );
      out << idx << " " << prev << " "
	  << i << " "
	  << curve.x( i ) << " " 
	  << curve.smooth_x( i, k ) << " " 
	  << d1 << " " << d2 << " " << kappa << " " 
	  << endl;

    }
  return 0;

}
