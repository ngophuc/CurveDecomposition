///////////////////////////////////////////////////////////////////////////////
// Converts a Freeman chaincode into an XFIG.
///////////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <deque>

#include "ImaGene/base/Arguments.h"
#include "ImaGene/base/StandardArguments.h"
#include "ImaGene/base/Proxy.h"
#include "ImaGene/base/Vector.h"
#include "ImaGene/base/Shapes.h"
#include "ImaGene/timetools/Clock.h"
#include "ImaGene/dgeometry2d/C4CTangentialCover.h"
#include "ImaGene/dgeometry2d/C4CIteratorOnFreemanChain.h"
#include "ImaGene/digitalnD/C4CIteratorOnFreemanChainSurface.h"
#include "ImaGene/digitalnD/KnSpace.h"
#include "ImaGene/digitalnD/KnCharSet.h"
#include "ImaGene/digitalnD/KnShapes.h"
#include "ImaGene/digitalnD/KnSpaceScanner.h"
#include "ImaGene/digitalnD/C4CIteratorOnSurface.h"
#include "ImaGene/digitalnD/BelAdjacency.h"
#include "ImaGene/digitalnD/ObjectBoundary.h"
#include "ImaGene/digitalnD/Frame2D.h"
#include "ImaGene/helper/C4CTangentialCoverGeometry.h"
#include "ImaGene/helper/GlobalC4CGeometry.h"
#include "ImaGene/helper/ContourHelper.h"
#include "ImaGene/helper/ShapeHelper.h"
#include "ImaGene/helper/CharacteristicPolygon.h"

using namespace std;
using namespace ImaGene;


static Arguments args;

void buildTangentialCover( C4CTangentialCover & tcover, 
			   C4CIterator & cp, uint max_size )
{
  Clock::startClock();
  tcover.init( cp, max_size );
  long t = Clock::stopClock();
  cerr << "# Tangential cover in " << t << " ms."
       << " nbsurf=" << tcover.nbSurfels() 
       << " nbms=" << tcover.nbMaximalSegments() << endl;
}

class XFIGFilter {

private:
  ostream & m_out;

public:
  const int RESOLUTION;
  const double EPS_MAGNIFICATION;
  const double ONE_IN_CM;
  double xmin;
  double xmax;
  double ymin;
  double ymax;
  int color_number;

  XFIGFilter( ostream & out, int res = 1200, double eps_mag = 100.0, 
	      double one_in_cm = 0.5 )
    : m_out( out ), RESOLUTION( res ), EPS_MAGNIFICATION( eps_mag ),
      ONE_IN_CM( one_in_cm )
  {
    xmin = 0.0;
    ymin = 0.0;
    xmax = 29.7;
    ymax = 21.0;
    color_number = 32;
  }

  void setView( double x1, double y1, double x2, double y2 )
  {
    xmin = x1;
    xmax = x2;
    ymin = y1;
    ymax = y2;
  }

  void outputHeader()
  {
    m_out << "#FIG 3.2" << endl
	  << "Landscape" << endl
	  << "Center" << endl
	  << "Inches" << endl
	  << "A4" << endl
	  << EPS_MAGNIFICATION << endl
	  << "Single" << endl
	  << "-1" << endl 
      // color number for transparent color for GIF
      // export. -3=background, -2=None, -1=Default, 0-31 for standard
      // colors or 32- for user colors)
	  << RESOLUTION <<" 2" << endl;
  }
  int toX( double x )
  {
    return (int) floor( RESOLUTION*x * ONE_IN_CM / 2.54 );
  }
  int toY( double y )
  {
    return (int) floor( RESOLUTION*y * ONE_IN_CM / 2.54 );
  }

  void fitInView( deque< pair<int,int> > & inside_idx,
		  const vector<double> & x,
		  const vector<double> & y )
  {
    bool inside = false;
    int first = 0;
    for ( int i = 0; i < (int) x.size(); ++i )
      {
	if ( ( x[ i ] >= xmin ) && ( x[ i ] <= xmax ) &&
	     ( y[ i ] >= ymin ) && ( y[ i ] <= ymax ) )
	  {
	    if ( ! inside ) 
	      {
		inside = true;
		first = i;
	      }
	  }
	else
	  {
	    if ( inside )
	      {
		inside_idx.push_back( make_pair( first, i ) );
		inside = false;
	      }
	  }
      }
    if ( inside )
      inside_idx.push_back( make_pair( first, x.size() ) );
  }

  /**
   * @param line_style -1 = Default, 0 = Solid, 1 = Dashed, 2 =
   * Dotted, 3 = Dash-dotted, 4 = Dash-double-dotted, 5 =
   * Dash-triple-dotted
   */
  void outputPolyline( const vector<double> & x,
		       const vector<double> & y,
		       int first,
		       int after_last,
		       int color = 0,
		       int thickness = 1,
		       int line_style = 0,
		       int depth = 50 )
  {
    //param 6 :fill color
    //param 9: fill style
    //param 4: epaisseur
    //param 5: couleur trait
    m_out << "2 1 " << line_style << " "
	  << thickness << " " 
	  << color << " "
	  << "1" << " " // fill color
	  << depth << " " // depth
	  << "-1 -1 " // pen_style, area_fill (enumeration type, -1 = no fill)
	  << "3.14 " // style_val (1/80 inch, spec for dash/dotted lines)
	  << "0 0 " // join_style, cap_style (only used for POLYLINE)
	  << "-1 0 0 " // radius (1/80 inch, radius of arc-boxes), forward_arrow (0: off, 1: on), backward_arrow (0: off, 1: on)
	  << after_last - first << " " << endl;
    //	      << x.size() << " " << endl;
    for ( int i = first; i < after_last; ++i )
      m_out << toX( x[ i ] ) << " " << toY( y[ i ] ) << endl;
  }

  /**
   * @param line_style -1 = Default, 0 = Solid, 1 = Dashed, 2 =
   * Dotted, 3 = Dash-dotted, 4 = Dash-double-dotted, 5 =
   * Dash-triple-dotted
   */
  void outputPolylineInView( const vector<double> & x,
			     const vector<double> & y,
			     int color = 0,
			     int thickness = 1,
			     int line_style = 0,
			     int depth = 50 )
  {
    deque< pair<int,int> > inside_idx;

    fitInView( inside_idx, x, y );
    
    for ( deque< pair<int,int> >::const_iterator it = inside_idx.begin();
	  it != inside_idx.end();
	  ++it )
      {
	outputPolyline( x, y, it->first, it->second,
			color, thickness, line_style, depth );
      }
  }

  void outputBox( double x1, double y1, double x2, double y2,
		  int color = 0,
		  int thickness = 1,
		  int line_style = 0,
		  int depth = 50,
		  int fill_color = 0,
		  int fill_style = 0 )
  {
    double x1p = x1 > xmin ? x1 : xmin;
    double y1p = y1 > ymin ? y1 : ymin;
    double x2p = x2 < xmax ? x2 : xmax;
    double y2p = y2 < ymax ? y2 : ymax;
    if ( ( x1p > x2p ) || ( y1p > y2p ) )
      return;
    
    m_out << "2 2 " << line_style << " "
	  << thickness << " " << color << " "
	  << fill_color << " " // fill_color
	  << depth << " " // depth
	  << "-1 " // pen_style
	  << fill_style << " " // (enumeration type, -1 = no fill)
	  << "0.000 " // style_val (1/80 inch, spec for dash/dotted lines)
	  << "0 0 " // join_style, cap_style (only used for POLYLINE)
	  << "-1 0 0 " // radius (1/80 inch, radius of arc-boxes), forward_arrow (0: off, 1: on), backward_arrow (0: off, 1: on)
	  << "5" << " " << endl;
    //	      << x.size() << " " << endl;
    m_out <<  toX( x1p ) << " " << toY( y1p ) << " "
	  <<  toX( x2p ) << " " << toY( y1p ) << " "
	  <<  toX( x2p ) << " " << toY( y2p ) << " "
	  <<  toX( x1p ) << " " << toY( y2p ) << " "
	  <<  toX( x1p ) << " " << toY( y1p ) << " " << endl;
  }

  void outputTriangle( double x1, double y1, 
		       double x2, double y2,		       
		       double x3, double y3,
		       int color = 0,
		       int thickness = 1,
		       int line_style = 0,
		       int depth = 50,
		       int fill_color = 0,
		       int fill_style = 0 )
  {
    if ( ( x1 < xmin ) || ( x2 > xmax ) || ( x3 > xmax )
	 || ( y1 < ymin ) || ( y2 > ymax ) || ( y3 > ymax ) )
      return;
    
    m_out << "2 3 " << line_style << " "
	  << thickness << " " << color << " "
	  << fill_color << " " // fill_color
	  << depth << " " // depth
	  << "-1 " // pen_style
	  << fill_style << " " // (enumeration type, -1 = no fill)
	  << "0.000 " // style_val (1/80 inch, spec for dash/dotted lines)
	  << "0 0 " // join_style, cap_style (only used for POLYLINE)
	  << "-1 0 0 " // radius (1/80 inch, radius of arc-boxes), forward_arrow (0: off, 1: on), backward_arrow (0: off, 1: on)
	  << "4" << " " << endl;
    //	      << x.size() << " " << endl;
    m_out <<  toX( x1 ) << " " << toY( y1 ) << " "
	  <<  toX( x2 ) << " " << toY( y2 ) << " "
	  <<  toX( x3 ) << " " << toY( y3 ) << " "
	  <<  toX( x1 ) << " " << toY( y1 ) << endl;
  }

  void outputEllipse( double x1, double y1, double x2, double y2,
		      int color = 0,
		      int thickness = 1,
		      int line_style = 0,
		      int depth = 50,
		      int fill_color = 0,
		      int fill_style = 0 )
  {
    if ( ( x1 < xmin ) || ( x2 > xmax )
	 || ( y1 < ymin ) || ( y2 > ymax ) )
      return;
    
    m_out << "1 2 " << line_style << " "
	  << thickness << " " << color << " "
	  << fill_color << " " // fill_color
	  << depth << " " // depth
	  << "-1 " // pen_style
	  << fill_style << " " // (enumeration type, -1 = no fill)
	  << "0.000 " // style_val (1/80 inch, spec for dash/dotted lines)
	  << "1 " // direction (always 1)
	  << "0.0 " // angle wrt x-axis 
	  << toX( ( x1 + x2 ) / 2.0 ) << " "
	  << toY( ( y1 + y2 ) / 2.0 ) << " "
	  << toX( ( x2 - x1 ) / 2.0 ) << " "
	  << toY( ( y2 - y1 ) / 2.0 ) << " ";
// 	  << "0 0 " // center_x center_y 
// 	  << "0 0 "; // radius_x radius_y 

// 	  << "0 0 " // join_style, cap_style (only used for POLYLINE)
// 	  << "-1 0 0 " // radius (1/80 inch, radius of arc-boxes), forward_arrow (0: off, 1: on), backward_arrow (0: off, 1: on)
// 	  << "5" << " " << endl;
    //	      << x.size() << " " << endl;
    m_out <<  toX( x1 ) << " " << toY( y1 ) << " "
	  <<  toX( x2 ) << " " << toY( y2 ) << endl;
  }

  void outputShadedTriangle( double x1, double y1, 
			     double x2, double y2,		       
			     double x3, double y3,
			     int color = 0,
			     int thickness = 1,
			     int fill_color = 0,
			     double shade_level = 0.20,
			     int depth = 50 )
  {
    int fill_style = (int) floor( shade_level * 20.0 );
    outputTriangle( x1, y1, x2, y2, x3, y3, 
		    color, thickness, 0, depth,
		    fill_color, fill_style );
  }

  
  void outputBoxedPoint( double i, double j,
			 float size = 0.35,
			 int color = 0,
			 int thickness = 0,
			 int fill_color = 0,
			 double shade_level = 0.20,
			 int depth = 60 )
  {
    double x1 = i;
    double y1 = j;
    int fill_style = (int) floor( shade_level * 20.0 );
    outputBox( x1 - size, y1 - size, x1 + size, y1 + size,
	       color, thickness, 0, depth, fill_color, fill_style );
  }

  void outputCircularPoint( double i, double j,
			    float size = 0.15,
			    int color = 0,
			    int thickness = 1,
			    int fill_color = 0,
			    double shade_level = 0.10,
			    int depth = 60 )
  {
    double x1 = i;
    double y1 = j;
    int fill_style = (int) floor( shade_level * 20.0 );
    outputEllipse( x1 - size, y1 - size, x1 + size, y1 + size,
		   color, thickness, 0, depth, fill_color, fill_style );
  }


  void outputBoxedPixel( int i, int j, bool inside = true )
  {
    if ( inside )
      outputBoxedPoint( i + 0.5, j + 0.5, 0.35, 0, 1, 0, 0.2, 60 );
    else
      outputBoxedPoint( i + 0.5, j + 0.5, 0.35, 0, 1, 0, 0.0, 60 );
  }
  void outputCircularPixel( int i, int j, bool inside = true )
  {
    if ( inside )
      outputCircularPoint( i + 0.5, j + 0.5, 0.15, 0, 1, 0, 0.6, 60 );
    else
      outputCircularPoint( i + 0.5, j + 0.5, 0.15, 0, 1, 0, 0.0, 60 );
  }

  void outputBoxedPointel( int i, int j, bool inside = true )
  {
    if ( inside )
      outputBoxedPoint( i, j, 0.2, 0, 1, 0, 0.8, 30 );
    else
      outputBoxedPoint( i, j, 0.2, 0, 1, 0, 0.0, 30 );
  }
  void outputCircularPointel( int i, int j, bool inside = true )
  {
    if ( inside )
      outputCircularPoint( i, j, 0.2, 0, 1, 0, 0.8, 30 );
    else
      outputCircularPoint( i, j, 0.2, 0, 1, 0, 0.0, 30 );
  }

  void outputLine( double x1, double y1, double x2, double y2,
		   int color = 0,
		   int thickness = 1,
		   int line_style = 0,
		   int depth = 50 )
  {
    vector<double> x( 2 );
    vector<double> y( 2 );
    x[ 0 ] = x1;
    y[ 0 ] = y1;
    x[ 1 ] = x2;
    y[ 1 ] = y2;
    outputPolyline( x, y, 0, 2, color, thickness, line_style, depth );
  }

  void outputGrid( double x1, double y1, double x2, double y2,
		   uint div_x, uint div_y,
		   int color = 0,
		   int thickness = 1,
		   int line_style = 0,
		   int depth = 60 )
  {
    double y1p = y1 > ymin ? y1 : ymin;
    double y2p = y2 < ymax ? y2 : ymax;
    if ( y1p <= y2p )
      for ( uint i = 0; i <= div_x; ++i )
	{
	  double x = x1 + ((x2 -x1)*i)/(double)div_x;
	  if ( ( x >= xmin ) && ( x <= xmax ) )
	    outputLine( x, y1p, x, y2p, color, thickness, line_style, depth );
	}

    double x1p = x1 > xmin ? x1 : xmin;
    double x2p = x2 < xmax ? x2 : xmax;
    if ( x1p <= x2p )
      for ( uint j = 0; j <= div_y; ++j )
	{
	  double y = y1 + ((y2 -y1)*j)/(double)div_y;
	  if ( ( y >= ymin ) && ( y <= ymax ) )
	    outputLine( x1p, y, x2p, y, color, thickness, line_style, depth );
	}
  }

  int outputNewColor( string rgb )
  {
    m_out << "0 " << color_number << " " << rgb << endl;
    return color_number++;
  }

  void outputText( double x, double y, string txt,
		   float angle,
		   int alignment = 0,
		   int ps_font = 0, // times roman
		   int font_size = 12,
		   int color = 0,
		   int depth = 50 )
  {
    if ( ( x < xmin ) || ( x > xmax )
	 || ( y < ymin ) || ( y > ymax ) )
      return;

    m_out << "4 " << alignment << " "
	  << color << " " 
	  << depth << " "
	  << "0 " // pen_style (not used)
	  << ps_font << " " // PostScript font number
	  << font_size << " "
	  << angle << " "
	  << "2 " // Font flag = PostScript font
	  << "0.0 0.0 " // height and length
	  << toX( x ) << " " << toY( y ) << " " 
	  << txt.c_str() << "\\001" << endl;

  }

  /**
   * @param mode 0: pixel interieur, 1: pixel exterieur, 2: pointel
   * boite gris sombre, 3: pointel boite vide, 4: pointel boite noire,
   * 5: pointel boite U, 6: pointel boite L, 7: pointel circulaire
   * plein, 8: pointel circulaire vide, 9: cross, 10: bracket [, 11: bracket ]
   */
  void outputElement( double i, double j, uint mode, int color = 0 )
  {
    double rpix = 0.35;
    double rptl = 0.2;

    switch ( mode ) {
    case 0:
      outputBoxedPoint( i, j, rpix, color, 1, color, 0.2, 60 ); 
      break;
    case 1:
      outputBoxedPoint( i, j, rpix, color, 1, color, 0, 60 ); 
      break;
    case 2:
      outputBoxedPoint( i, j, rptl, color, 1, color, 0.6, 30 );
      break;
    case 3:
      outputBoxedPoint( i, j, rptl, color, 1, color, 0.0, 30 );
      break;
    case 4:
      outputBoxedPoint( i, j, rptl, color, 1, color, 1.0, 30 );
      break;
    case 5:
      outputShadedTriangle( i - rptl , j - rptl, 
			    i + rptl, j - rptl,
			    i - rptl, j + rptl, 
			    color, 1, color, 0.0, 30 );
      outputShadedTriangle( i + rptl , j - rptl, 
			    i - rptl, j + rptl,
			    i + rptl, j + rptl, 
			    color, 1, color, 0.6, 30 );
      break;
    case 6:
      outputShadedTriangle( i - rptl , j - rptl, 
			    i + rptl, j - rptl,
			    i - rptl, j + rptl,
			    color, 1, color, 0.6, 30 );
      outputShadedTriangle( i + rptl , j - rptl, 
			    i - rptl, j + rptl,
			    i + rptl, j + rptl, 
			    color, 1, color, 0.0, 30 );
      break;
    case 7:
      outputCircularPoint( i, j, rptl, 
			   color, 1, color, 0.8, 30 );
      break;
    case 8:
      outputCircularPoint( i, j, rptl, 
			   color, 1, color, 0.0, 30 );
      break;
    case 9: // cross X
      outputLine( i - rptl, j - rptl, 
		  i + rptl, j + rptl, 
		  color, 2, color, 20 );
      outputLine( i + rptl, j - rptl, 
		  i - rptl, j + rptl, 
		  color, 2, color, 20 );
      break;
    case 10: // bracket [
      outputLine( i, j - rptl, 
		  i, j + rptl, color, 2, color, 20 );
      outputLine( i, j - rptl, 
		  i + rptl, j - rptl, color, 2, color, 20 );
      outputLine( i, j + rptl, 
		  i + rptl, j + rptl, color, 2, color, 20 );
      break;
    case 11: // bracket ]
      outputLine( i, j - rptl, 
		  i, j + rptl, color, 2, color, 20 );
      outputLine( i, j - rptl, 
		  i - rptl, j - rptl, color, 2, color, 20 );
      outputLine( i, j + rptl, 
		  i - rptl, j + rptl, color, 2, color, 20 );
      break;
    }
  }  
};


void getCoords( vector<double> & x,
		vector<double> & y,
		const FreemanChain & c,
		int start, 
		int end )
{
  int pos = 0;
  FreemanChain::const_iterator it = c.begin();
  for ( ;
	it != c.end();
	++it )
    {
      if ( ( start <= pos )
	   && ( ( end <= 0 ) || ( pos < end ) ) )
	{
	  Vector2i xy( *it );
	  x.push_back( xy.x() );
	  y.push_back( xy.y() );
	  // cout << xy.x() << " " << xy.y() << endl;
	}
      ++pos;
    }
  if ( ( start <= pos )
       && ( ( end <= 0 ) || ( pos < end ) ) )
    {
      Vector2i xy( *it );
      x.push_back( xy.x() );
      y.push_back( xy.y() );
    }
  
}


void outputDSS( XFIGFilter & filter, 
		const C4CSegment & dss, 
		const Frame2D & frame,
		int color,
		uint ms_mode = 0,
		double ms_coef_mode = 0.0 )
{
  Vector2i u( frame.transformPoint( dss.u() ) );
  Vector2i up( frame.transformPoint( dss.up() ) );
  Vector2i l( frame.transformPoint( dss.l() ) );
  Vector2i lp( frame.transformPoint( dss.lp() ) );
  if ( ms_mode & 32 )
    {
      filter.outputElement( u.x(), u.y(), 5, color );
      filter.outputElement( up.x(), up.y(), 5, color );
    }
  if ( ms_mode & 16 )
    {
      filter.outputElement( l.x(), l.y(), 6, color );
      filter.outputElement( lp.x(), lp.y(), 6, color );
    }
  Vector2i c( frame.transformPoint( dss.c_n() ) );
  Vector2i cp( frame.transformPoint( dss.cp_n() ) );
  if ( ms_mode & 64 )
    {
      filter.outputElement( c.x(), c.y(), 9, color );
      filter.outputElement( cp.x(), cp.y(), 9, color );
    }
  Vector2D c1 = frame.transformPoint( dss.project( dss.c_n(), dss.u() ) );
  Vector2D c2 = frame.transformPoint( dss.project( dss.c_n(), dss.l() ) );
  Vector2D cp1 = frame.transformPoint( dss.project( dss.cp_n(), dss.u() ) );
  Vector2D cp2 = frame.transformPoint( dss.project( dss.cp_n(), dss.l() ) );
  vector<double> x( 5 );
  vector<double> y( 5 );
  x[ 0 ] = x[ 4 ] = c1.x();
  y[ 0 ] = y[ 4 ] = c1.y();
  x[ 1 ] = c2.x(); y[ 1 ] = c2.y();
  x[ 2 ] = cp2.x(); y[ 2 ] = cp2.y();
  x[ 3 ] = cp1.x(); y[ 3 ] = cp1.y();

  vector<uint> z;
  Mathutils::cfrac( z, dss.a() < 0 ? -dss.a() : dss.a(), dss.b() );

  int line_style = 0; // solid
  if ( ms_mode & 4 )
    line_style = ( z.size() % 2 ) == 0 
      ? 0 // even
      : 1; // odd

  int thickness = ( ms_mode % 8 ) != 0 ? 2 : 1;
  filter.outputPolylineInView( x, y, color, thickness, line_style, 45 );

  ostringstream out_str;
  if ( ms_mode & 1 )
    out_str << "(" << dss.a() << "," << dss.b() << ";" << dss.mu() << ")"; 
  if ( ms_mode & 2 )
    {
      out_str << "[" << z[ 0 ]; 
      for ( uint i = 1; i < z.size(); ++i )
	out_str << "," << z[ i ];
      out_str << "]";
    }
  // cerr << out_str.str() << endl;
  if ( out_str.str() != "" )
    {
      float angle = atan2( (double) dss.a(), (double) dss.b() );
      double xbase = ( c1.x() + cp1.x() ) / 2.0;
      double ybase = ( c1.y() + cp1.y() ) / 2.0;
      xbase += ms_coef_mode * ( c2.x() - c1.x() );
      ybase += ms_coef_mode * ( c2.y() - c1.y() );

      filter.outputText( xbase, ybase,
			 out_str.str(),
			 -frame.angleToX( angle ),
			 1, 0, 12, color, 20 );
    }

}

void outputMS( XFIGFilter & filter, 
	       KnSpace* ks,
	       const C4CTangentialCover & tcover,
	       uint idx_ms,
	       int color,
	       uint ms_mode,
	       double ms_coef_mode )
{
  Frame2D frame;
  frame.init( ks, 0, 1 );
  Proxy<C4CIteratorOnSurface> it 
    ( (C4CIteratorOnSurface*) tcover.getSurfelFront( idx_ms ) );
  Kn_sid sbel = it->current();
  frame.setSurfelFrame( sbel, it->trackDir() );
  outputDSS( filter, tcover.getMaximalSegment( idx_ms ).dss, 
	     frame, color, ms_mode,
	     ( idx_ms % 2 ) == 0 ? - ms_coef_mode : 1.5 + ms_coef_mode );

}


/**
 * Computes the characteristic polygon of the digital contour from its
 * tangential cover.
 *
 */
void 
computeCharacteristicPolygon
( vector<double> & x,
  vector<double> & y,
  //  KnSpace* ks,
  const C4CTangentialCover & tcover,
  C4CTangentialCoverGeometry & tcover_geometry )
{
  CharacteristicPolygon P;
  P.init( tcover, tcover_geometry );
  for ( uint i = 0; i < P.edges().size(); ++i )
    {
      Vector2D P1 = P.edges()[ i ].getP1();
      Vector2D P2 = P.edges()[ i ].getP2();
      x.push_back( P1.x() );
      y.push_back( P1.y() );
//       x.push_back( P2.x() );
//       y.push_back( P2.y() );
    }
  x.push_back( x[ 0 ] );
  y.push_back( y[ 0 ] );
//   Mathutils::AngleComputer ac;
//   uint nb_ms = tcover.nbMaximalSegments();
//   for ( uint i = 0; i < nb_ms; ++i )
//     {
//       uint i_prev = ( i + nb_ms - 1 ) % nb_ms;
//       const C4CTangentialCoverGeometry::MaximalSegmentGeometry & msg_prev 
// 	= tcover_geometry.geometry( i_prev );

//       const C4CTangentialCover::MaximalSegment & ms
// 	= tcover.getMaximalSegment( i );
//       const C4CTangentialCoverGeometry::MaximalSegmentGeometry & msg
// 	= tcover_geometry.geometry( i );
//       const Frame2D & frame
// 	= tcover_geometry.sgeometry( ms.front_surfel_idx ).frame;

//       uint i_next = ( i + 1 ) % nb_ms;
//       const C4CTangentialCoverGeometry::MaximalSegmentGeometry & msg_next
// 	= tcover_geometry.geometry( i_next );

//       if ( ac.less( msg_prev.angle_to_x, msg.angle_to_x ) )
// 	{ // Convex to back
// 	  cerr << "convex back" << endl;
// 	  Vector2D Up( ms.dss.u().x(), ms.dss.u().y() - 0.5 );
// 	  if ( ms.dss.a() < 0 ) Up.x() -= 0.5; 
// 	  else Up.x() += 0.5; 
// 	  Vector2D U1 = frame.transformPoint( Up );
// 	  x.push_back( U1.x() );
// 	  y.push_back( U1.y() );
// 	}
//       else
// 	{ // Concave to back
// 	  cerr << "concave back" << endl;
// 	  Vector2D Lp( ms.dss.l().x(), ms.dss.l().y() + 0.5 );
// 	  if ( ms.dss.a() > 0 ) Lp.x() -= 0.5; 
// 	  else Lp.x() += 0.5; 
// 	  Vector2D L1 = frame.transformPoint( Lp );
// 	  x.push_back( L1.x() );
// 	  y.push_back( L1.y() );
// 	}

//       if ( ac.less( msg.angle_to_x, msg_next.angle_to_x ) )
// 	{ // Convex to front
// 	  cerr << "convex front" << endl;
// 	  Vector2D Up( ms.dss.up().x(), ms.dss.up().y() - 0.5 );
// 	  if ( ms.dss.a() <= 0 ) Up.x() -= 0.5; 
// 	  else Up.x() += 0.5; 
// 	  Vector2D U1 = frame.transformPoint( Up );
// 	  x.push_back( U1.x() );
// 	  y.push_back( U1.y() );
// 	}
//       else
// 	{ // Concave to front
// 	  cerr << "concave front" << endl;
// 	  Vector2D Lp( ms.dss.lp().x(), ms.dss.lp().y() + 0.5 );
// 	  if ( ms.dss.a() >= 0 ) Lp.x() -= 0.5; 
// 	  else Lp.x() += 0.5; 
// 	  Vector2D L1 = frame.transformPoint( Lp );
// 	  x.push_back( L1.x() );
// 	  y.push_back( L1.y() );
// 	}

//     }
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
  StandardArguments::addDigitalArgs( args, 2, false, false );
  args.addOption( "-step", "-step <h>: the discretization step or resolution, the closer to 0, the finer.", "1.0" );
  args.addBooleanOption( "-header", "-header: add XFIG header." );
  args.addOption( "-contour", "-contour <color> <thickness>: displays the contour as a polygonal line", "0", "3" );
  args.addOption( "-restrict", "-restrict <begin> <past_end>: specifies which part of the contour is to be displayed.", "0", "0" );
  args.addOption( "-ms", "-ms <begin> <past_end>: specifies which maximal segments are to be displayed.", "0", "0" );
  args.addOption( "-ms_display", "-ms_display <BW|COLOR> <bits 6:b 5:u 4:l 3:t 2:p 1:z 0:s> <shift>: specifies how maximal segments are displayed. b: begin and end of MS, u: upper leaning points, l: lower leaning points, t: thick box, p: display parity, z: display frac, s: display slope ; shift= shift for displaying text", "BW", "6", "0.5" );
  args.addOption( "-pixel_grid", "-pixel_grid <color> <thickness>: displays an embedding grid.", "1", "1" );
  args.addOption( "-pointel_grid", "-pointel_grid <color> <thickness>: displays an embedding grid.", "2", "1" );
  args.addOption( "-view", "-view <x1> <y1> <x2> <y2>: specifies the subset of Z2 to display", "0", "0", "64", "64" );
  args.addOption( "-unitcm", "-unitcm <u>: tells the length in cm of a unit step.", "0.5" );
  args.addOption( "-pgm_file", "-pgm_file <fname>: specifies a PGM file to display a set of pixels as background (def. is 128).", "128" );
  args.addOption( "-threshold", "-threshold <val>: threshold value for binarizing PGM gray values (def. is 128).", "128" );
  args.addOption( "-pixels", "-pixels <number WB> <0/1>: displays white (bit 1)/black (bit 0) pixels with 0: boxes, 1: disks", "1", "0" );
  args.addBooleanOption( "-auto_center", "-auto_center: automatically centers the contour in the space." );
  args.addOption( "-cpolygon", "-cpolygon <color> <thickness>: compute characteristic polygon.", "2", "2" );
  args.addBooleanOption( "-convexity", "-convexity: outputs the convexity of each surfel as given by the characteristic polygon." );
  if ( ( argc <= 0 ) 
       || ! args.readArguments( argc, argv ) ) 
    {
      cerr << args.usage( "freeman2fig", 
			  "Visualizes a digital contour (given as a freeman chaincode on the standard input) by outputing an XFIG file (or part of file).",
			  "" )
	   << endl;
      return 1;
    }

  // -------------------------------------------------------------------------
  // Read Freeman chain and creates space/contour.

  FreemanChain c;
  FreemanChain::read( cin, c );
  if ( ! cin.good() )
    {
      cerr << "Error reading Freeman chain code." << endl;
      return 2;
    }

  KnSpace* ks;
  if ( args.check( "-auto_center" ) )
    ks = ShapeHelper::makeSpaceFromFreemanChain( c );
  else 
    {
      // Build space.
      uint d = StandardArguments::dim( args );
      if ( d != 2 )
	{
	  cerr << "Dimension is 2." << endl;
	  return 2;
	}
      Kn_size sizes[ d ];
      StandardArguments::fillSizes( args, sizes );
      ks = new KnSpace( 2, sizes );
    }

  // -------------------------------------------------------------------------
  // Prepare XFIG filter.

  XFIGFilter filter( cout, 1200, 100.0, 
		     args.getOption( "-unitcm" )->getFloatValue( 0 ) );

  if ( args.check( "-header" ) )
    filter.outputHeader();

  bool color_mode = 
    args.getOption( "-ms_display" )->getValue( 0 ) == "COLOR";
  uint ms_mode = 
    (uint) args.getOption( "-ms_display" )->getIntValue( 1 );
  double ms_coef_mode = 
    args.getOption( "-ms_display" )->getDoubleValue( 2 );
  int c0 = 0;
  int c1 = 0;
  if ( color_mode )
    {
      c0 = filter.outputNewColor( "#000080" );
      filter.outputNewColor( "#0000ff" );
      filter.outputNewColor( "#8080ff" );
      filter.outputNewColor( "#ff80ff" );
      filter.outputNewColor( "#ff00ff" );
      filter.outputNewColor( "#ff0080" );
      filter.outputNewColor( "#ff0000" );
      filter.outputNewColor( "#800000" );
      filter.outputNewColor( "#008000" );
      c1 = filter.outputNewColor( "#00ff00" );
    }

  double xmin = args.getOption( "-view" )->getFloatValue( 0 );
  double ymin = args.getOption( "-view" )->getFloatValue( 1 );
  double xmax = args.getOption( "-view" )->getFloatValue( 2 );
  double ymax = args.getOption( "-view" )->getFloatValue( 3 );
  filter.setView( xmin, ymin, xmax, ymax );

  // -------------------------------------------------------------------------
  

  if ( args.check( "-contour" ) )
    {
      int b = args.getOption( "-restrict" )->getIntValue( 0 );
      int e = args.getOption( "-restrict" )->getIntValue( 1 );
      vector<double> x;
      vector<double> y;
      getCoords( x, y, c, b, e );
      filter.outputPolylineInView( x, y, 
				   args.getOption( "-contour" )->getIntValue( 0 ),
				   args.getOption( "-contour" )->getIntValue( 1 ),
				   0, 40 );
    }
  if ( args.check( "-pointel_grid" ) )
    {
      uint div_x = (uint) floor( xmax - xmin );
      uint div_y = (uint) floor( ymax - ymin );
      filter.outputGrid( xmin, ymin, xmax, ymax, div_x, div_y,
			 args.getOption( "-pointel_grid" )->getIntValue( 0 ),
			 args.getOption( "-pointel_grid" )->getIntValue( 1 ),
			 0, 60 ); 
    }
  if ( args.check( "-pixel_grid" ) )
    {
      uint div_x = (uint) floor( xmax - xmin ) - 1;
      uint div_y = (uint) floor( ymax - ymin ) - 1;
      filter.outputGrid( xmin + 0.5, ymin + 0.5, xmax - 0.5, ymax - 0.5,
			 div_x, div_y,
			 args.getOption( "-pixel_grid" )->getIntValue( 0 ),
			 args.getOption( "-pixel_grid" )->getIntValue( 1 ),
			 1, 70 ); 
    }

  uint pixel_bw_mode = args.getOption( "-pixels" )->getIntValue( 0 );
  uint pixel_style_mode = args.getOption( "-pixels" )->getIntValue( 1 );

  if ( args.check( "-pgm_file" ) )
    {
      KnSpace* ks;
      KnCharSet* voxset;
      string fname = args.getOption( "-pgm_file" )->getValue( 0 );
      ifstream fin( fname.c_str() );
      if ( ! fin.good() )
	{
	  cerr << "Error opening PGM file: " 
	       << fname
	       << endl;
	  return 2;
	}
      uint threshold = (uint) args.getOption( "-threshold" )->getIntValue( 0 );
      if ( ! ShapeHelper::importFromPGM( fin, ks, voxset, threshold ) )
	{
	  cerr << "Error reading PGM file: "  << fname
	       << endl;
	  return 2;
	}
      KnSpaceScanner scan2( *ks, 
			    ks->ufirstCell( ks->dim() ),
			    ks->ulastCell( ks->dim() ) );
      
      Kn_uid p = scan2.begin();
      for ( Kn_uid last_y = scan2.last( p, 1 );
	    p <= last_y; 
	    p += scan2.gotonext( 1 ) )
	{
	  for ( Kn_uid last_x = scan2.last( p, 0 ); 
		p <= last_x; 
		p++ ) // NB: 'scan.gotonext( 0 )' == 1;
	    { //... whatever
	      uint i = ks->udecodeCoord( p, 0 );
	      uint j = ks->udecodeCoord( p, 1 );
	      bool vox = (*voxset)[ p ];
	      if ( pixel_style_mode == 0 ) // boxes
		{
		  if ( vox && ( pixel_bw_mode & 1 ) )
		    filter.outputBoxedPixel( i, j, vox );
		  else if ( ( ! vox ) && ( pixel_bw_mode & 2 ) )
		    filter.outputBoxedPixel( i, j, vox );
		}
	      else // disks
		{
		  if ( vox && ( pixel_bw_mode & 1 ) )
		    filter.outputCircularPixel( i, j, vox );
		  else if ( ( ! vox ) && ( pixel_bw_mode & 2 ) )
		    filter.outputCircularPixel( i, j, vox );
		}
	    }
	}
      delete voxset;
      delete ks;
      
      fin.close();
    }


  // KnSpace* ks = ShapeHelper::makeSpaceFromFreemanChain( c );
  C4CIteratorOnFreemanChain itfc;
  itfc.init( c.begin(), true );
  C4CIteratorOnFreemanChainSurface itfcs;
  itfcs.init( ks, itfc );

  // double dh = (double) args.getOption( "-step" )->getDoubleValue( 0 );

  // build tangential cover
  Clock::startClock();
  C4CTangentialCover tcover;
  buildTangentialCover( tcover, itfcs, 0 );
  C4CTangentialCoverGeometry tcover_geometry;
  Frame2D frame;
  frame.init( ks, 0, 1 );
  C4CTangentialCoverGeometry::MSGeometryComputerByLP geocomp;
  tcover_geometry.init( tcover, itfcs, geocomp, frame );
  long t2 = Clock::stopClock();
  cerr << "# Tangential cover & geometry in " << t2 << " ms." << endl;


  if ( args.check( "-ms" ) )
    {
      uint min = args.getOption( "-ms" )->getIntValue( 0 );
      uint max = args.getOption( "-ms" )->getIntValue( 1 );
      if ( ( max == 0 )
	   || ( max > tcover.nbMaximalSegments() ) )
	max = tcover.nbMaximalSegments();
      if ( max >= min )
	{
	  int c = c0;
	  for ( uint idx = min; idx < max; ++idx )
	    {
	      outputMS( filter, ks, tcover, idx, c, ms_mode, ms_coef_mode );
	      ++c;
	      if ( c > c1 ) c = c0;
	    }
	}
    }

  if ( args.check( "-cpolygon" ) )
    {
      uint color = args.getOption( "-cpolygon" )->getIntValue( 0 );
      uint thickness = args.getOption( "-cpolygon" )->getIntValue( 1 );

      vector<double> x;
      vector<double> y;
      computeCharacteristicPolygon( x, y, tcover, tcover_geometry );

      filter.outputPolylineInView( x, y,
				   color, thickness,
				   0, 60 );
      for ( uint i = 0; i < x.size(); ++i )
	{
	  filter.outputElement( x[ i ], y[ i ], 9, color );
	}
    }

  if ( args.check( "-convexity" ) )
    {
      int b = args.getOption( "-restrict" )->getIntValue( 0 );
      int e = args.getOption( "-restrict" )->getIntValue( 1 );
      vector<double> x;
      vector<double> y;
      getCoords( x, y, c, b, e );
      CharacteristicPolygon P;
      P.init( tcover, tcover_geometry );
      // P.extractEdges();
      for ( uint i = 0; i < ( x.size() - 1 ); ++i )
	{
	  filter.outputLine( x[ i ], y[ i ], x[ i + 1 ], y[ i + 1 ],
			     (int) ( P.shape()[ i ] - '0' ),
			     args.getOption( "-contour" )->getIntValue( 1 ),
			     0, 40 );
	}
      for ( uint i = 0; i < P.edges().size(); ++i )
	{
	  uint i1 = P.edges()[ i ].idx1;
	  uint i2 = P.edges()[ i ].idx2;
	  filter.outputLine( x[ i1 ], y[ i1 ], x[ i2 ], y[ i2 ],
			     4, 2,
			     0, 20 );
	}
    }


  delete ks;
  return 0;
}
