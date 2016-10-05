///////////////////////////////////////////////////////////////////////////////
// Generates a new polygon from a freeman contour.
///////////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <deque>
#include <vector>
#include <sstream>
#include <string>


#include "ImaGene/base/Arguments.h"
#include "ImaGene/base/StandardArguments.h"
#include "ImaGene/base/Proxy.h"
#include "ImaGene/base/Vector.h"
#include "ImaGene/base/Shapes.h"
#include "ImaGene/timetools/Clock.h"
#include "ImaGene/dgeometry2d/FreemanChain.h"
#include "ImaGene/dgeometry2d/C4CIteratorOnFreemanChain.h"
#include "ImaGene/digitalnD/GridEmbedder.h"
#include "ImaGene/digitalnD/KnSpace.h"
#include "ImaGene/digitalnD/KnSpaceScanner.h"
#include "ImaGene/digitalnD/KnCharSet.h"
#include "ImaGene/digitalnD/KnRCellSet.h"
#include "ImaGene/digitalnD/KnShapes.h"
#include "ImaGene/digitalnD/C4CIteratorOnBdry.h"
#include "ImaGene/digitalnD/C4CIteratorOnFreemanChainSurface.h"
#include "ImaGene/digitalnD/BelAdjacency.h"
#include "ImaGene/digitalnD/ObjectBoundary.h"
#include "ImaGene/digitalnD/Frame2D.h"
#include "ImaGene/helper/ContourHelper.h"
#include "ImaGene/helper/ShapeHelper.h"
#include "ImaGene/helper/CharacteristicPolygon.h"
#include "ImaGene/helper/ChangeQuadrant.h"

using namespace std;
using namespace ImaGene;


static Arguments args;

class XFIGFilter {

private:
  ostream & m_out;

public:
  const int RESOLUTION;
  const double EPS_MAGNIFICATION;
  const double ONE_IN_CM;
  int h;
  int v;
  double xmin;
  double xmax;
  double ymin;
  double ymax;
  int color_number;

  XFIGFilter( ostream & out, int res = 1200, double eps_mag = 100.0, 
	      double one_in_cm = 0.5, int h = 1, int v = 1 )
    : m_out( out ), RESOLUTION( res ), EPS_MAGNIFICATION( eps_mag ),
      ONE_IN_CM( one_in_cm ), h (h ), v (v )
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
    return (int) floor( RESOLUTION * x * ONE_IN_CM * h / 2.54 );
  }
  int toY( double y )
  {
    return (int) floor( RESOLUTION * y * ONE_IN_CM * v / 2.54 );
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
    {
      //m_out << " x y : " << x[ i ] << " " << y[ i ] << endl;
      m_out << toX( x[ i ] ) << " " << toY( y[ i ] ) << endl;
    }
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
    //deque< pair<int,int> > inside_idx;
    //fitInView( inside_idx, x, y );
    //for ( deque< pair<int,int> >::const_iterator it = inside_idx.begin();
	//  it != inside_idx.end();
	//  ++it )
      //{
	outputPolyline( x, y, 0, x.size(),
			color, thickness, line_style, depth );
     // }
  }


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
	  double x = x1 + ((x2 -x1)*i)/(double)div_x + 0.5;
	  if ( ( x >= xmin ) && ( x <= xmax ) )
	    outputLine( x, y1p, x, y2p, color, thickness, line_style, depth );
	}

    double x1p = x1 > xmin ? x1 : xmin;
    double x2p = x2 < xmax ? x2 : xmax;
    if ( x1p <= x2p )
      for ( uint j = 0; j <= div_y; ++j )
	{
	  double y = y1 + ((y2 -y1)*j)/(double)div_y + 0.5;
	  if ( ( y >= ymin ) && ( y <= ymax ) )
	    outputLine( x1p, y, x2p, y, color, thickness, line_style, depth );
	}
  }


void DisplayContourXFIG( FreemanChain & fc, 
                         int color, int thickness, 
                         bool filled, bool closeCnt, int zoom )
{
  //param 6 :fill color
  //param 9: fill style
  //param 4: epaisseur
  //param 5: couleur trait
  cout <<setiosflags(ios_base::fixed);
  cout<< setprecision(0);
  
  cout << "\n 2 " << ((closeCnt)? 1 : 0 )<< " 0 "  << thickness << " " << color  
       << "  10 19 -1 "<<((filled)? 0:-1)  <<" 0.000 0 0 -1 0 0 " << fc.chain.size()+((closeCnt)? 1 : 0 ) 
       << " "  <<endl ;
  for (FreemanChain::const_iterator it = fc.begin();it != fc.end(); ++it )
   cout << ((*it).x())*RESOLUTION*zoom << " " << ((*it).y())*RESOLUTION*zoom << " " ;        
  if(closeCnt)
   cout << (((*(fc.begin())).x()))*RESOLUTION*zoom << " " << ((*(fc.begin())).y())*RESOLUTION*zoom << " " ;        
}


void displayPixel( float pixelx, float pixely, int color, int colorFill, 
                   double tailleh, int taillev, int depth )
{
  Vector2D p0 (pixelx+tailleh/2.0, pixely+taillev/2.0);
  Vector2D p1 (pixelx+tailleh/2.0, pixely-taillev/2.0);
  Vector2D p2 (pixelx-tailleh/2.0, pixely-taillev/2.0);
  Vector2D p3 (pixelx-tailleh/2.0, pixely+taillev/2.0);
  
  cout << setprecision(15) << "\n 2 3 0 1 " << color  << " " << colorFill 
       << " " << depth<<    " -1 20 0.000 0 0 -1 0 0 5 " <<endl ;  
  
  cout << p0.x()*RESOLUTION << " " << p0.y()*RESOLUTION << " " << p1.x()*RESOLUTION 
       << " " << p1.y()*RESOLUTION << " " << p2.x()*RESOLUTION << " " 
       << p2.y()*RESOLUTION << " "<< p3.x()*RESOLUTION << " " << p3.y()*RESOLUTION 
       << " "  << p0.x()*RESOLUTION << " " << p0.y()*RESOLUTION << endl;
}


void displayContourPixelsXFIG( FreemanChain & fc, 
                               int color, int colorfill, 
                               int depth, int h, int v )
{
  FreemanChain::const_iterator it = fc.begin();
  int nb = fc.chain.size();
  float valx = fc.x0;
  float valy = fc.y0; 
  displayPixel(valx, valy, color, colorfill, h, v, depth);
  for( int i = 1 ; i < nb ; ++i )
  {  
    if( fc.code( i ) == 0 ) valx += h;
    else if( fc.code( i ) == 1 ) valy -= v;
    else if( fc.code( i ) == 2 ) valx -= h;
    else valy += v;
    displayPixel(valx, valy, color, colorfill, h, v, depth);
    it.nextInLoop();
    
  }
}

void displayCroix( float pixelx, float pixely, int color, 
                   int width, int heigth, int depth )
{ 
  cout <<setiosflags(ios_base::fixed);
  cout<< setprecision(0);
  cout<< "\n 2 1 0 5 " << color << " " << color << " " << depth << " -1 -1 0.000 0 0 -1 0 0 2"<< endl;
  cout << (pixelx-width/2.0)*RESOLUTION << " " << (pixely-heigth/2.0)*RESOLUTION << " " 
       << (pixelx+width/2.0)*RESOLUTION << " " << (pixely+heigth/2.0)*RESOLUTION << endl; 
  cout<< "\n 2 1 0 5 " << color << " " << color << " " << depth << " -1 -1 0.000 0 0 -1 0 0 2"<< endl;
  cout << (pixelx+width/2.0)*RESOLUTION << " " << (pixely-heigth/2.0)*RESOLUTION << " " 
       << (pixelx-width/2.0)*RESOLUTION << " " << (pixely+heigth/2.0)*RESOLUTION << endl;   
}


void displayText( float pixelx, float pixely, string txt, 
                  float angle, int color, int font_size )
{ 
  cout <<setiosflags(ios_base::fixed);
  cout<< setprecision(0);
  cout << "\n 4 1 " << color << " " << "1 0 0 " << font_size << " " << angle << " 2 0.0 0.0 " << " "
       << pixelx * RESOLUTION << " " << pixely * RESOLUTION << " " << txt.c_str() << "\\001" << endl;  
}

};


///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////


struct Point2i {
  float x;
  float y;
};

int norm( const Point2i & A, const Point2i & B )
{
  float dx = A.x - B.x;
  float dy = A.y - B.y;
  return ( dx >= 0 ? dx : -dx )
    + ( dy >= 0 ? dy : -dy );
}

ostream & operator<<( ostream & out, const Point2i & p )
{
  out << "(" << p.x << "," << p.y << ")";
}

class Set {

public:
  virtual bool isInside( float x, float y ) const = 0;
  bool isInside( const Point2i & p ) const
  {
    return this->isInside( p.x, p.y );
  }
};

class DSL : public Set {

public:
  using Set::isInside;

  static int nb_requetes;
  int a;
  int b;
  float mu;
  int nb_point_teste;
  float remainder( const Point2i & p ) const {
    return remainder( p.x, p.y );
  }
  float remainder( float x, float y ) const {
    return a*x - b*y;
  }
  virtual bool isInside( float x, float y ) const {
    nb_requetes++;
    float r = remainder( x, y );
    return ( r >= mu ) && ( r < ( mu + a + b ) );
  }
  bool isWU( float x, float y ) {
    float r = remainder( x, y );
    return ( r == ( mu - 1 ) );
  }
  bool isWL( float x, float y ) {
    float r = remainder( x, y );
    return ( r == ( mu + a + b ) );
  }
};

int DSL::nb_requetes = 0;



typedef std::vector<int> vint;

void display( ostream & out, 
	      const vint & u, const vint & p, const vint & q )
{
  int c = u.size() - 1;
  out << p[ c ] << "/" << q[ c ] << "=[" << u[ 0 ];
  for ( int i = 1; i <= c; i++ )
    out << ',' << u[ i ];
  out << "]";
}

/**
 * Computes the bezout coefficient of (a,b) a*u-b*v = 1;
 */
void Bezout( const vint & a, const vint & b, int & u, int & v )
{
  int c = a.size() - 1;
  if ( c == 0 ) 
    { 
      u = b[ c ] - 0; 
      v = a[ c ] - 1;
    }
  else if ( ( c & 0x1 ) == 0 ) 
    {
      u = b[ c ] - b[ c - 1 ];
      v = a[ c ] - a[ c - 1 ];
    }
  else
    {
      u = b[ c - 1 ];
      v = a[ c - 1 ];
    }
}

bool checkBezout( const vint & a, const vint & b, int & u, int & v )
{
  int c = a.size() - 1;
  return ( a.back() * u - b.back() * v ) == 1;
}

/**
 * Updates the slope coefficients in case of weak upper or lower
 * leaning point.
 */
void updateSlope( bool wu, int delta, vint & u, vint & p, vint & q )
{
  int c = u.size() - 1;
  bool simple = wu == ( ( c & 0x1 ) == 0 );
  int p_cm1 = ( c > 0 ) ? p[ c - 1 ] : 1; // defines p_-1 
  int q_cm1 = ( c > 0 ) ? q[ c - 1 ] : 0; // defines q_-1 
  if ( simple )
    { // ( WU and even ) or ( WL and odd )
      if ( delta == 1 )
	{
	  u[ c ]++; 
	  p[ c ] += p_cm1;
	  q[ c ] += q_cm1;
	}
      else 
	{
	  u.push_back( delta );
	  p.push_back( p[ c ] * delta + p_cm1 );
	  q.push_back( q[ c ] * delta + q_cm1 );
	}
    }
  else
    { // ( WU and odd ) or ( WL and even )
      u[ c ]--; 
      p[ c ] -= p_cm1;
      q[ c ] -= q_cm1;
      u.push_back( 1 );
      p.push_back( p[ c ] + p_cm1 );
      q.push_back( q[ c ] + q_cm1 );
      ++c;
      if ( delta == 1 )
	{
	  u[ c ]++; 
	  p[ c ] += p[ c - 1 ]; // p_cm1;
	  q[ c ] += q[ c - 1 ]; // q_cm1;
	}
      else
	{
	  u.push_back( delta );
	  p.push_back( p[ c ] * delta + p[ c - 1 ] );
	  q.push_back( q[ c ] * delta + q[ c - 1 ] );
	}
    }
}

class DSS {

public:
  DSL myD;
  Point2i U1;
  Point2i U2;
  Point2i L1;
  Point2i L2;

  Point2i X0;
  Point2i XF;
  vint u;
  vint p;
  vint q;

  void init()
  {
    //myD.a = 0;
    //myD.b = 1;
    u.clear();
    p.clear();
    q.clear();
    u.push_back( 0 );
    p.push_back( 0 );
    q.push_back( 1 );
  }

  int sum() const
  {
    int acc = 0;
    for ( int i = 0; i < u.size(); ++i )
      acc += u[ i ];
    return acc;
  }

/**
   * Determines the characteristics of the DSS when it is known that
   * the set S under consideration is a piece of DSL in the first quadrant.
   * Stops as soon as the slope is the same as the expected DSL.
   * #nb_requete_isInside <= 2*sum_coefs_u + 2
   */
bool smartRec( const DSL & S, const Point2i & x0, const Point2i & xf)
{
    init();
    int nb_points_testes = 0;
    myD.a = S.a;
    myD.b = S.b;
    myD.mu = myD.remainder( x0 );
    U1 = x0;
    L1 = x0;
    U2.x = x0.x + 1;
    U2.y = x0.y;
    L2 = U2;
    if ( ! ( S.isInside( U1 ) && S.isInside( U2 ) ) )
      return false;
    X0 = x0;
    XF = xf;
    Point2i MU;
    Point2i ML;
    bool inside = true;
    while ( inside && ( ( p.back() != S.a ) || ( q.back() != S.b ) ) )
    {
	// Compute Bezout coefficients.
	int ap; int bp;
	Bezout( p, q, bp, ap );
	int a = p.back();
	int b = q.back();
        //cout << "pq.back(): " << a << " " << b << endl;
        //cout << "- z=";
	//display( cout, u, p, q );
	//cout << " Bezout=(" << bp << "," << ap << ") "
	//     << ap << "/" << bp << " "
	//     << ( checkBezout( p, q, bp, ap ) ? "OK" : "ERREUR" )
	//     << endl;

	// MU = U1 + k(b,a) + (b',a')
	MU = U1;  
	MU.x += b - bp; 
	MU.y += a - ap;
	// ML = L1 + (k+1)(b,a) - (b',a')
	ML = L1;
	ML.x += bp;
	ML.y += ap;
	int delta = 1;
	int loop = 0; // 1: WU, 2: WL, 3: out
	while ( loop == 0 )
	  {
	    MU.x += b; MU.y += a;
	    ML.x += b; ML.y += a;

	    //cout << "  (d=" << delta << ") MU" << MU << endl;
	    //cout << "  (d=" << delta << ") ML" << ML << endl;
            if ( ( MU.y <= xf.y ) && ( ML.x <= xf.x ) ) nb_points_testes += 2;
            else if ( ( MU.y <= xf.y ) || ( ML.x <= xf.x ) ) nb_points_testes += 1;
             
	    if ( ( MU.y <= xf.y ) && S.isInside( MU ) )
	    { loop = 1; break; }
	    else if ( ( ML.x <= xf.x ) && S.isInside( ML ) )
	    { loop = 2; break; }
	    else if ( ( MU.y > xf.y ) || ( ML.x > xf.x ) )
	    { loop = 3; break; }
	    else 
	      ++delta;  
	  }
	if ( loop == 1 )
	  { // increase slope with WU
	    //cout << "  => MU" << MU << endl;
	    updateSlope( true, delta, u, p, q );
	    U2 = MU;
	    L1.x = ML.x - bp;
	    L1.y = ML.y - ap;
	    if ( norm( x0, L1 ) >= norm( x0, U2 ) )
	      { 
		L1.x -= b;
		L1.y -= a;
	      }
	    L2 = L1;
	    //cout << "     U1" << U1 << " L1" << L1 << endl;
	  }
	else if ( loop == 2 )
	  { // decrease slope with WL
	    //cout << "  => ML" << ML << endl;
	    updateSlope( false, delta, u, p, q );
	    L2 = ML;
	    U1.x = MU.x - b + bp;
	    U1.y = MU.y - a + ap;
	    if ( norm( X0, U1 ) >= norm( X0, L2 ) )
	      { 
		U1.x -= b;
		U1.y -= a;
	      }
	    U2 = U1;
	    //cout << "     U1" << U1 << " L1" << L1 << endl;
	  }
	if ( loop == 3 ) inside = false; 
      }
    myD.a = p.back();
    myD.b = q.back();
    myD.mu = myD.remainder( U1 );
    myD.nb_point_teste = nb_points_testes;
    return true;
}

};

///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////

/**
 * @return the pgcd of a*h and b*v.
 */
int get_pgcd(int ah, int bv)
{
  int pgcd = 0;
  if( ah == 0 || bv == 0 )
  {
    return 1;
  }
  while(1)
  {
    pgcd = ah % bv;
    if(pgcd == 0)
    {
      pgcd = bv;
      break;
    }
    ah = bv;
    bv = pgcd;
  }
  return pgcd;
}

/**
 * Compute the new characteristics (new_a,new_b,new_mu) of segment that is the   
 * subsampling of the segment of characteristics (a,b,mu).
 */
void CalculCharacteristicSegmenthv( ImaGene::ChangeQuadrant::segmentMLP & OriginalSegmentMLPhv, ImaGene::ChangeQuadrant::segmentMLP & OriginalSegmentMLP, int h, int v )
{
  int a = OriginalSegmentMLP.a;
  int b = OriginalSegmentMLP.b;
  int mu = OriginalSegmentMLP.mu;
  int quadrant = OriginalSegmentMLP.Quadrant;
  int SI, SS, gcd; 
  int a1 = ( a >= 0 ? a : -a );
  int b1 = ( b >= 0 ? b : -b );
  if( a1 != 0 && b1 != 0 )
  {
    gcd = get_pgcd( a1 * h , b1 * v );
    OriginalSegmentMLPhv.a = ( a * h ) / gcd;
    OriginalSegmentMLPhv.b = ( b * v ) / gcd;
  }
  else
  {
    if( a1 == 0 )
    {
      OriginalSegmentMLPhv.a = 0;
      OriginalSegmentMLPhv.b = b;
      gcd = b1 * v;
    }
    else
    {
      OriginalSegmentMLPhv.a = a;
      OriginalSegmentMLPhv.b = 0;
      gcd = a1 * h;
    }
  }
  int new_a = OriginalSegmentMLPhv.a;
  int new_b = OriginalSegmentMLPhv.b;
  int p = ( new_a >=0 ? new_a : - new_a ) + ( new_b >=0 ? new_b : - new_b ) ;
  int p1 = new_a >=0 ? new_a : - new_a;
  int p2 = new_b >=0 ? new_b : - new_b;
  int Q1, Q2, Q3, R1, R2, R3;
  int val1, val2, val3;
  val1 = a1 + b1 - 1;
  val2 = ( quadrant == 1 || quadrant == 3 ) ? mu + 2 * a1 + 2 * b1 - 1 : mu + 2 * a1 + b1 - 1;
  val3 = ( quadrant == 1 || quadrant == 3 ) ? 2 * mu + 3 * a1 + 3 * b1 - 1 : 2 * mu + 3 * a1 + b1 - 1;
  if( val1 >= 0 || val1 % gcd == 0 )
  {
    Q1 = val1 / gcd;
    R1 = val1 % gcd;
  }
  else
  {
    Q1 = val1 / gcd - 1;
    R1 = val1 - Q1 * gcd;
  }
  if( val2 >= 0 || val2 % gcd == 0 )
  {
    Q2 = val2 / gcd;
    R2 = val2 % gcd;
  }
  else
  {
    Q2 = val2 / gcd - 1;
    R2 = val2 - Q2 * gcd;
  }
  if( val3 >= 0 || val3 % gcd == 0 )
  {
    Q3 = val3 / gcd;
    R3 = val3 % gcd;
  }
  else
  {
    Q3 = val3 / gcd - 1;
    R3 = val3 - Q3 * gcd;
  }
  if( R2 <= R1 ) SI = 0;
   else SI = 1;
  if( R3 <= R2 ) SS = 0;
   else SS = 1; 
  int mu1 = ( quadrant == 1 || quadrant == 3 ) ? - p + Q2 - Q1 + SI : - p1 + Q2 - Q1 + SI;
  int mu2 = ( quadrant == 1 || quadrant == 3 ) ? Q3 - Q2 + SS : p2 + Q3 - Q2 + SS;
  OriginalSegmentMLPhv.mu = ( mu1 < mu2 ) ? mu1 : mu2;
  OriginalSegmentMLPhv.mu_sup = ( mu1 < mu2 ) ? mu2 - 1 : mu1 - 1;
}

/**
 * Compute the FreemanChain of the segment of characteristics (a,b,mu) between 
 * the first point fp and the last point lp.
 */
void InnerFreemanChain( FreemanChain & InnerFreemanChainSegment, 
                        Vector2i & fp, Vector2i & lp, 
                        int a, int b, int mu )
{
  int xf = fp.x();
  int yf = fp.y();
  int xl = lp.x();
  int yl = lp.y();
  int x, y;
  char code;
  if( a >= 0 && b >= 0 )
  {
    while( xf != xl || yf != yl )
    {
      x = xf + 1; y = yf;
      int r = a * x - b * y;
      if( r >= mu && r < mu + a + b ) code = '0';
      else { x = xf; y = yf + 1; code = '1';  }
      xf = x; yf = y;
      InnerFreemanChainSegment.chain += code;
      //if( first ){ first_code = code == 0 ? 0 : 1; first = false; }
    }
  }
  else
  {
    if( a >= 0 && b < 0 )
    {
      while( xf != xl || yf != yl )
      {
        x = xf - 1; y = yf;
        int r = a * x - b * y;
        if( r >= mu && r < mu + a - b ) code = '2';
        else { x = xf; y = yf + 1; code = '1'; }
        xf = x; yf = y;
        InnerFreemanChainSegment.chain += code;
        //if( first ){ first_code = code == '2' ? 2 : 1; first = false; }
      }
    }
    else
    {
      if( a < 0 && b < 0 )
      {
        while( xf != xl || yf != yl )
        {
          x = xf - 1; y = yf;
          int r = - a * x + b * y;
          if( r >= mu && r < mu - a - b ) code = '2';
          else { x = xf; y = yf - 1; code = '3'; }
          xf = x; yf = y;
          InnerFreemanChainSegment.chain += code;
          //if( first ){ first_code = code == 2 ? 2 : 3; first = false; }
        }
      }
      else
      {
        while( xf != xl || yf != yl )
        {
          x = xf + 1; y = yf;
          int r = - a * x + b * y;
          if( r >= mu && r < mu - a + b ) code = '0';
          else { x = xf; y = yf - 1; code = '3'; }
          xf = x; yf = y;
          InnerFreemanChainSegment.chain += code;
          //if( first ){ first_code = code == 0 ? 0 : 3; first = false; }
        }
      }
    }
  }
}

/**
 * Compute the new coefficients (new_a,new_b) of the slope  where (new_a,new_b) = 1.
 */
void getslope( Vector2D & fp, Vector2D & lp, int & new_a, int & new_b, int & quadrant )
{   
  int a,b;
  a = lp.y() - fp.y();
  b = lp.x() - fp.x();
  if( a >= 0 && b >= 0 ) quadrant = 0;
  else if( a >= 0 && b < 0 ) quadrant = 1;
  else if( a < 0 && b < 0 ) quadrant = 2;
  else quadrant = 3;
  int pgcd;
  if( a != 0 && b != 0 )
  {
    pgcd = get_pgcd(( a > 0 ) ? a : -a,( b > 0 ) ? b : -b);
  }
  else
  {
    if( a == 0 ) pgcd = ( b > 0 ) ? b : -b;
    else pgcd = ( a > 0 ) ? a : -a;
  }
  new_a = a / pgcd;
  new_b = b / pgcd;
}

/**
 * @return the first value of delta where v divides delta*a
 * and h divides delta*b.
 */
int Deltamin( int a, int b, int h, int v )
{
  int g = get_pgcd( h , v);
  int g1 = get_pgcd( a , v);
  int g2 = get_pgcd( b , h);
  int g3 = get_pgcd( h , v / g1 );
  int g4 = get_pgcd( v , h / g2 );
  int deltamin;
  if( h == v ) deltamin = v;
  else if( g == 1 ) deltamin = ( h * v ) / ( g1 * g2 );
  else deltamin = ( h * v * g ) / ( g1 * g2 * g3 * g4 );
  
  return deltamin;
}

/**
 * @return the number of patterns of slope (a,b).
 */
int NumberofPatterns( ImaGene::ChangeQuadrant::segmentMLP & OriginalSegmentMLP )
{
  int numberofpatterns;
  int a = OriginalSegmentMLP.a;
  int b = OriginalSegmentMLP.b;
  int mu = OriginalSegmentMLP.mu;
  Vector2i first_point;
  first_point.x() = OriginalSegmentMLP.first_point.x();
  first_point.y() = OriginalSegmentMLP.first_point.y();
  Vector2i last_point;
  last_point.x() = OriginalSegmentMLP.last_point.x();
  last_point.y() = OriginalSegmentMLP.last_point.y();
  int diffx = last_point.x() - first_point.x();
  int diffy = last_point.y() - first_point.y();
  if( b != 0 ) numberofpatterns = ( diffx > 0 ? diffx : - diffx ) / ( b > 0 ? b : - b );
  else numberofpatterns = ( diffy > 0 ? diffy : - diffy ) / ( a > 0 ? a : - a );

  return numberofpatterns;
}

/**
 * @return the first code of inner FreemanChain.
 */
void getfirstcode( int & first_code, const FreemanChain & c, 
                   bool convex_back, bool convex_front, 
                   int indexfc1, int a, int b )
{
  int code_c;
  if( ( ( convex_back && !convex_front ) || ( !convex_back && convex_front ) ) && ( a == 0 || b == 0 ) ) 
    first_code = c.code( indexfc1 + 1 );
  else first_code = c.code( indexfc1 );
}

/**
 * After extracting edges of the contour, computes the characteristics of each 
 * edge by the tiling (h,v).
 */
void TransformPolygon2Polygonhv( vector<ImaGene::ChangeQuadrant::segmentMLP> & VectOriginalSegmentMLPhv,                vector<ImaGene::ChangeQuadrant::segmentMLP> & VectGlobalSegmentMLP, 
vector<ImaGene::CharacteristicPolygon::Edge> & VectSegmentMLP, 
const FreemanChain & c, int h, int v, int x0, int y0 )
{
  if ( ( h == 0 ) || ( v == 0 ) ) return;
  //cout << "h: " << h << "   " << "v: " << v << "   FreemanChain: " << c.chain << endl;  
  //cout << "x0 y0: " << c.x0 << "  " << c.y0 << endl;
  uint nb = c.chain.size();
  //cout << "Number of points of the initial contour: " << nb << endl;
  //cout << "=======================================================================================" << endl;
  int number_SegmentMLP = VectSegmentMLP.size(); 
  bool first_segment = true, first_step = true;
  int nbpointtested;
  int numbertotalpointtested = 0;
  int nb_equal_XY = 0;
  int nb_deltamin = 0;
  FreemanChain chv;
  for( int i = 0; i < number_SegmentMLP ; ++i )
  {
    ImaGene::ChangeQuadrant::segmentMLP GlobalSegmentMLP;
    int indexfc1 = VectSegmentMLP[i].idx1;
    int indexfc2 = VectSegmentMLP[i].idx2;
    GlobalSegmentMLP.convex_back = VectSegmentMLP[i].convex_back;
    GlobalSegmentMLP.convex_front = VectSegmentMLP[i].convex_front;
    Vector2D fp = VectSegmentMLP[i].getP1();
    Vector2D lp = VectSegmentMLP[i].getP2();
    //cout << "indexfc1: " << indexfc1 << "   " << "indexfc2: " << indexfc2 << endl;
    //cout << "convex_back: " << GlobalSegmentMLP.convex_back << "   " << "convex_front: " 
    //     << GlobalSegmentMLP.convex_front << "   " << "  fp: " << fp.x() << " ; " << fp.y() 
    //     << "  lp: " << lp.x() << " ; " << lp.y() << endl;
    GlobalSegmentMLP.first_point.x() = fp.x() - x0;
    GlobalSegmentMLP.first_point.y() = fp.y() - y0;
    GlobalSegmentMLP.last_point.x() = lp.x() - x0;
    GlobalSegmentMLP.last_point.y() = lp.y() - y0;
    if( ( fp.x()-x0 == lp.x()-x0 ) && ( fp.y()-y0 == lp.y()-y0 ) ) continue;
    int Quadrant;
    int a; int b; 
    int first_code;
    getslope( fp, lp, a, b, Quadrant );
    GlobalSegmentMLP.Quadrant = Quadrant;
    getfirstcode( first_code, c, GlobalSegmentMLP.convex_back, GlobalSegmentMLP.convex_front, indexfc1, a, b );
    GlobalSegmentMLP.first_code = first_code;
    GlobalSegmentMLP.a = a;
    GlobalSegmentMLP.b = b;
    //cout << "I- Global Quadrant:  First_point_segment:( " << GlobalSegmentMLP.first_point.x() << "," 
    //     << GlobalSegmentMLP.first_point.y() << " )   Last_point_segment:( " << GlobalSegmentMLP.last_point.x() 
    //     << "," << GlobalSegmentMLP.last_point.y() 
    //     << " )   (a,b):( " << GlobalSegmentMLP.a << "," << GlobalSegmentMLP.b  
    //    << " )    Quadrant of this FreemanChain: " << Quadrant
    //     << endl;
    ChangeQuadrant segment;
    ImaGene::ChangeQuadrant::segmentMLP OriginalSegmentMLP;
    segment.TransformGlobalQuad2OriginalQuad( GlobalSegmentMLP, OriginalSegmentMLP, h, v );
    FreemanChain InnerFreemanChainSegment;
    InnerFreemanChain( InnerFreemanChainSegment, OriginalSegmentMLP.first_point, OriginalSegmentMLP.last_point,       OriginalSegmentMLP.a, OriginalSegmentMLP.b, OriginalSegmentMLP.mu);
    OriginalSegmentMLP.fc.chain = InnerFreemanChainSegment.chain;
    OriginalSegmentMLP.first_code = GlobalSegmentMLP.first_code;
    GlobalSegmentMLP.fc.chain = OriginalSegmentMLP.fc.chain;
    VectGlobalSegmentMLP.push_back( GlobalSegmentMLP );
    //cout << "II- Original Quadrant:  First_point_segment:( " << OriginalSegmentMLP.first_point.x() << "," 
    //     << OriginalSegmentMLP.first_point.y() 
    //     << " )Last_point_segment:( " << OriginalSegmentMLP.last_point.x() << "," << OriginalSegmentMLP.last_point.y() 
    //     << " )   Slope (a,b,mu,mu_sup) of this FreemanChain: (" << OriginalSegmentMLP.a << ";" 
    //     << OriginalSegmentMLP.b << ";" <<  OriginalSegmentMLP.mu << ";" << OriginalSegmentMLP.mu_sup << ")" 
    //     << "     FreemanChain in the Original Quadrant: " << OriginalSegmentMLP.fc.chain 
    //     << endl;

    // Multiscale of this segment.
    int FX = ( fp.x() - x0 ) / h;
    int FY = ( fp.y() - y0 ) / h;
    int LX = ( lp.x() - x0 ) / v;
    int LY = ( lp.y() - y0 ) / v;
    if( ( FX == LX ) && ( FY == LY ) ) 
    {
      ++nb_equal_XY; 
      continue;  
    }
    ImaGene::ChangeQuadrant::segmentMLP OriginalSegmentMLPhv;
    CalculCharacteristicSegmenthv( OriginalSegmentMLPhv, OriginalSegmentMLP, h, v);
    int delta, deltamin;
    deltamin = Deltamin( a, b, h, v );
    delta = NumberofPatterns( OriginalSegmentMLP );
    DSS dss;
    DSL D;
    ImaGene::ChangeQuadrant::segmentMLP OriginalDSSMLPhv;
    segment.TransformFCOriginalQuad2OriginalQuadhv ( OriginalSegmentMLP, OriginalSegmentMLPhv, h, v);
    if( delta >= 2 * deltamin )
    {
      ++nb_deltamin;
      nbpointtested = 0;
      OriginalDSSMLPhv.nbpointtested = nbpointtested;
      int quadrant = OriginalSegmentMLP.Quadrant;
      OriginalDSSMLPhv.a = OriginalSegmentMLPhv.a;
      OriginalDSSMLPhv.b = OriginalSegmentMLPhv.b;
      OriginalDSSMLPhv.mu = OriginalSegmentMLPhv.mu;
      OriginalDSSMLPhv.fc.chain = OriginalSegmentMLPhv.fc.chain;
      //cout << "III- Original Quadrant (h,v): " << " FreemanChain in the Original Quadrant (h,v): " 
      //     << OriginalDSSMLPhv.fc.chain 
      //     << "    S'_hv(DSS.a,DSS.b,DSS.mu) = ( " << OriginalDSSMLPhv.a << " , " 
      //     << OriginalDSSMLPhv.b << " , " << OriginalDSSMLPhv.mu << " ) " << "  nb_points_tested:  " 
      //     << OriginalDSSMLPhv.nbpointtested
      //     << endl;
    }
    else 
    {
      //segment.TransformFCOriginalQuad2OriginalQuadhv ( OriginalSegmentMLP, OriginalSegmentMLPhv, h, v);
      Vector2i fp_hv;
      Vector2i lp_hv;
      //cout << "III- Original Quadrant (h,v):  First_point_segment_hv:(" << OriginalSegmentMLPhv.first_point.x() << "," 
      //     << OriginalSegmentMLPhv.first_point.y() << " )   Last_point_segment (h,v):( " 
      //     << OriginalSegmentMLPhv.last_point.x() << "," << OriginalSegmentMLPhv.last_point.y() 
      //     <<  " )   (a,b,mu,mu_sup): ( " << OriginalSegmentMLPhv.a << "," << OriginalSegmentMLPhv.b << ";" 
      //     << OriginalSegmentMLPhv.mu << ";" << OriginalSegmentMLPhv.mu_sup <<  " )"
      //     << "     FreemanChain in the Original Quadrant hv: " << OriginalSegmentMLPhv.fc.chain 
      //     << endl;
      ImaGene::ChangeQuadrant::segmentMLP FirstSegmentMLPhv;
      segment.TransformOriginalQuad2FirstQuad( FirstSegmentMLPhv, OriginalSegmentMLPhv );
      //cout << "IV- First Quadrant (h,v):  First_point_segment_hv:(" << FirstSegmentMLPhv.first_point.x() << "," 
      //     << FirstSegmentMLPhv.first_point.y() << " )   Last_point_segment_hv:( " 
      //     << FirstSegmentMLPhv.last_point.x() << "," << FirstSegmentMLPhv.last_point.y() 
      //     <<  " )   (a,b,mu,mu_sup): ( " << FirstSegmentMLPhv.a << "," << FirstSegmentMLPhv.b << ";" 
      //     << FirstSegmentMLPhv.mu << ";" << FirstSegmentMLPhv.mu_sup << " )"
      //     << endl;
      Point2i P0, PF;
      P0.x = FirstSegmentMLPhv.first_point.x();
      P0.y = FirstSegmentMLPhv.first_point.y();
      PF.x = FirstSegmentMLPhv.last_point.x();
      PF.y = FirstSegmentMLPhv.last_point.y();
      D.a = FirstSegmentMLPhv.a;
      D.b = FirstSegmentMLPhv.b;
      D.mu = FirstSegmentMLPhv.mu;   
      dss.smartRec( D, P0, PF );
      nbpointtested = dss.myD.nb_point_teste;
      ImaGene::ChangeQuadrant::segmentMLP FirstDSSMLPhv;
      FirstDSSMLPhv.a = dss.myD.a;
      FirstDSSMLPhv.b = dss.myD.b;
      FirstDSSMLPhv.mu = dss.myD.mu;
      FirstDSSMLPhv.first_point.x() = P0.x;
      FirstDSSMLPhv.first_point.y() = P0.y;
      //cout << "V- First Quadrant (h,v): First_point_segment_hv:( " << P0.x << " , " << P0.y << " ) " 
      //     << "  Last_point_segment_hv:( " << PF.x << " , " << PF.y << " ) "
      //     << "  S'_hv(DSS.a,DSS.b,DSS.mu) = ( " << FirstDSSMLPhv.a << " , " << FirstDSSMLPhv.b << " , " 
      //     << FirstDSSMLPhv.mu << " ) " << "  nb_points_tested:  " << nbpointtested
      //     << endl;
      OriginalDSSMLPhv.nbpointtested = nbpointtested;
      segment.TransformFirstQuad2OriginalQuad( FirstDSSMLPhv, OriginalDSSMLPhv, OriginalSegmentMLPhv, GlobalSegmentMLP );
      OriginalDSSMLPhv.fc.chain = OriginalSegmentMLPhv.fc.chain;
      //cout << "VI- Original Quadrant (h,v): " << " FreemanChain in the Original Quadrant (h,v): " 
      //     << OriginalDSSMLPhv.fc.chain 
      //     << "    S'_hv(DSS.a,DSS.b,DSS.mu) = ( " << OriginalDSSMLPhv.a << " , " 
      //     << OriginalDSSMLPhv.b << " , " << OriginalDSSMLPhv.mu << " ) " << "  nb_points_tested:  " 
      //     << OriginalDSSMLPhv.nbpointtested
      //     << endl;
    //  if( OriginalDSSMLPhv.fc.chain != "" ) numbertotalpointtested += OriginalDSSMLPhv.nbpointtested;
    }
    if( OriginalDSSMLPhv.fc.chain != "" ) chv.chain += OriginalDSSMLPhv.fc.chain;
    VectOriginalSegmentMLPhv.push_back( OriginalDSSMLPhv );
    //cout << " ============================================================================================ " << endl;
  }
  //cout << "New Freeman Chain: " << chv.chain << endl;
  //cout << " Size of New Freeman Chain: " << chv.chain.size() << endl; 
  // cout << "numbertotalpointtested:  " << numbertotalpointtested << "   Nbsegmentshv: " << number_SegmentMLP << endl;
  // cout << "nb_equal_XY: " << nb_equal_XY << endl;
  // cout << "nb_deltamin: " << nb_deltamin << endl;
}


/**
 * Given a Freeman chain [fc] coding a 4-connected pixel loop, computes
 * its subsampling by the transformation:
 * X = ( x - x0 ) div h, 
 * Y = ( y - y0 ) div v.
 *
 * @param subc (output) the subsampled Freeman chain code (may contain spikes)
 * 
 * @param first_point the first point of the Freeman Chain.
 *
 * @param fc the input chain code.
 *
 * @param h the subsampling along x
 * @param v the subsampling along y
 * @param x0 the x-origin of the frame (X,Y) in (x,y)
 * @param y0 the y-origin of the frame (X,Y) in (x,y)
 */
void NewFreemanChainhv( FreemanChain & subc,
                        Vector2D first_point,
		        const FreemanChain & fc,
		        uint h, uint v,
		        int x0, int y0 )
{
  if ( ( h == 0 ) || ( v == 0 ) ) return;
  uint nb = fc.chain.size();
  Vector2D fxy( first_point );
  Vector2D I_fXY, I_nXY;
  Vector2D fXY( ( fxy.x() - x0 ) / h, ( fxy.y() - y0 ) / v );
  I_fXY.x() = (int) fXY.x();  
  I_fXY.y() = (int) fXY.y();
  Vector2D nxy;
  int valx, valy;
  for ( uint i = 0; i < nb; ++i )
  {
    if( fc.code( i ) == 0 ) { valx =  fxy.x() + 1 ; valy = fxy.y(); }
    else if( fc.code( i ) == 1 ) { valx =  fxy.x() ; valy = fxy.y() - 1; }
    else if( fc.code( i ) == 2 ) { valx =  fxy.x() - 1 ; valy = fxy.y(); }
    else { valx =  fxy.x() ; valy = fxy.y() + 1; }
    nxy.x() = valx;
    nxy.y() = valy;
    fxy.x() = nxy.x();
    fxy.y() = nxy.y();
    Vector2D nXY( ( nxy.x() - x0 ) / h, ( nxy.y() - y0 ) / v );
    I_nXY.x() = (int) nXY.x();  
    I_nXY.y() = (int) nXY.y();
    if ( ( I_nXY.x() != I_fXY.x() ) || ( I_nXY.y() != I_fXY.y() ) )
    {
      char code;
      if ( I_nXY.x() > I_fXY.x() )       code = '0';
      else if ( I_nXY.x() < I_fXY.x() )  code = '2';
      else if ( I_nXY.y() > I_fXY.y() )  code = '1';
      else                               code = '3';
      I_fXY = I_nXY;
      subc.chain += code;
    }    
  }
}
//**********************************************************************************************************************//
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
  args.addBooleanOption( "-header", "-header: add XFIG header." );
  args.addOption( "-ColorSize", "-ColorSize <color> <thickness>: displays the contour as a polygonal line", "0", "3" );
  args.addOption("-agrandissementEPS", "-agrandissmentEPS <coef> : agrandissement utilisé par eps (defaut=1)", "1");
  args.addOption( "-DisplayContourSurfelSampled", "-DisplayContourSurfelSampled <color> <thickness> ", "0", "3" );
  args.addOption( "-SegmentContourSampled", "-SegmentContourSampled <color> <thickness> ", "0", "3" );
  args.addOption( "-displayPixelSegment", "-displayPixelSegment <color> <fillcolor> <depth>", "0", "0", "3" );
  args.addOption( "-displayPixelSegmenthv", "-displayPixelSegmenthv <color> <depth>", "0", "3" );
  args.addOption( "-view", "-view <x1> <y1> <x2> <y2>: specifies the subset of Z2 to display", "0", "0", "64", "64" );
  args.addOption( "-unitcm", "-unitcm <u>: tells the length in cm of a unit step.", "0.5" );
  args.addOption( "-restrict", "-restrict <begin> <past_end>: specifies which part of the contour is to be displayed.", "0", "0" );
  args.addOption( "-displayabmu", "-displayabmu <COLOR> <font_size>", "0", "15" );
  args.addOption( "-pixel_grid", "-pixel_grid <color> <thickness>: displays an embedding grid.", "1", "1" );
  args.addOption( "-displayContourPixels", "-displayContourPixels <color> <colorfill> <depth>: displays a contour pixels.", "1", "1", "2" );
  args.addOption( "-displayContourPixelshv", "-displayContourPixelshv <color> <colorfill> <depth>: displays a contour pixels hv.", "1", "1", "1" );
  args.addOption( "-pointel_grid", "-pointel_grid <color> <thickness>: displays an embedding grid.", "2", "1" );
  args.addBooleanOption( "-auto_center", "-auto_center: automatically centers the shape in the PGM image." );
  args.addBooleanOption( "-Digitization_hv", "-Digitization_hv");
  args.addOption( "-hv", "-hv <h><v>: h and v are elements of the sampling.","1","1");
  if ( ( argc <= 0 ) 
       || ! args.readArguments( argc, argv ) ) 
    {
      cerr << args.usage( "test_Multiscale_segment", 
			  "Reconstructs a new polygon from a Freeman chaincode by the tiling (h,v).",
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
  int h = 1;
  int v = 1;
  if ( args.check( "-hv" ) )
  {
    h= args.getOption( "-hv" )->getIntValue( 0 );
    v= args.getOption( "-hv" )->getIntValue( 1 );      
  }
  C4CIteratorOnFreemanChain itfc;
  itfc.init( c.begin(), true );
  C4CIteratorOnFreemanChainSurface itfcs;
  itfcs.init( ks, itfc );
  Clock::startClock();
  C4CTangentialCover tcover;
  tcover.init( itfcs, 0 );
  C4CTangentialCoverGeometry tcover_geometry;
  Frame2D frame;
  frame.init( ks, 0, 1 );
  C4CTangentialCoverGeometry::MSGeometryComputerByLP geocomp;
  tcover_geometry.init( tcover, itfcs, geocomp, frame );
  long t2 = Clock::stopClock();
  cerr << "# Tangential cover & geometry in " << t2 << " ms." << endl;
  CharacteristicPolygon P;
  P.init( tcover, tcover_geometry );
  vector<ImaGene::CharacteristicPolygon::Edge> VectSegmentMLP;
  vector<ImaGene::ChangeQuadrant::segmentMLP> VectOriginalSegmentMLPhv;
  vector<ImaGene::ChangeQuadrant::segmentMLP> VectGlobalSegmentMLP;
  int x0 = 0;
  int y0 = 0;
  Clock::startClock();
  //VectSegmentMLP = P.edges();
  //for( int k = 1; k < 15 ; ++k )
  //for( int x0 = 0; x0 < k ; ++x0 )
  //for( int y0 = 0; y0 < k ; ++y0 )
  //FreemanChain inner_chain;
  //vector<uint> outer2inner;
  //vector<uint> inner2outer;
  //FreemanChain New_outer_chain;
  //Vector2D first_point( c.x0 , c.y0 );
  //for( int i = 0 ; i < 1000 ; ++i ) 
  //{ 
    VectSegmentMLP = P.edges();
    TransformPolygon2Polygonhv( VectOriginalSegmentMLPhv, VectGlobalSegmentMLP, VectSegmentMLP, c, h, v, x0, y0);
    
    //FreemanChain::innerContour( inner_chain, outer2inner, inner2outer, c, true);
    //NewFreemanChainhv( New_outer_chain, first_point, inner_chain, h, v, x0, y0 );
  //}
  long t3 = Clock::stopClock();
  cerr << "# Multiscale in " << t3 << " ms." << endl;
  //******************************************************************************************************//
  //                                         Prepare XFIG filter                                          //
  //*****************************************************************************************************//

  XFIGFilter filter( cout, 1200, 100.0, 
		     args.getOption( "-unitcm" )->getFloatValue( 0 ), h, v );

  if ( args.check( "-header" ) )
    filter.outputHeader();

  double xmin = args.getOption( "-view" )->getFloatValue( 0 );
  double ymin = args.getOption( "-view" )->getFloatValue( 1 );
  double xmax = args.getOption( "-view" )->getFloatValue( 2 );
  double ymax = args.getOption( "-view" )->getFloatValue( 3 );
  filter.setView( xmin, ymin, xmax, ymax );

  if ( args.check( "-displayContourPixels" ) )
  {
    int color = args.getOption( "-displayContourPixels" )->getIntValue( 0 );
    int colorfill = args.getOption( "-displayContourPixels" )->getIntValue( 1 );
    int thickness = args.getOption( "-displayContourPixels" )->getIntValue( 2 );
    int nb = VectGlobalSegmentMLP.size();
    FreemanChain fc;
    string s = " ";
    for( int i = 0; i < nb ; ++i )
      s += VectGlobalSegmentMLP[i].fc.chain;
    fc.chain = s;
    fc.x0 = VectGlobalSegmentMLP[0].first_point.x();
    fc.y0 = VectGlobalSegmentMLP[0].first_point.y();
    filter.displayContourPixelsXFIG(fc, color, colorfill, thickness, 1, 1);        
  }
   
  if ( args.check( "-displayContourPixelshv" ) )
  {
    int color = args.getOption( "-displayContourPixelshv" )->getIntValue( 0 );
    int colorfill = args.getOption( "-displayContourPixelshv" )->getIntValue( 1 );
    int thickness = args.getOption( "-displayContourPixelshv" )->getIntValue( 2 );
    int nb = VectOriginalSegmentMLPhv.size();
    FreemanChain fc;
    string s = " ";
    for( int i = 0 ; i < nb ; ++i )
      s += VectOriginalSegmentMLPhv[i].fc.chain;
    fc.chain = s;
    int x0 = VectGlobalSegmentMLP[0].first_point.x();
    int y0 = VectGlobalSegmentMLP[0].first_point.y();
    fc.x0 = x0 - x0 % h + (float)h / 2 - 0.5;
    fc.y0 = y0 + y0 % v - ( v - 1 ) - 0.5 + (float) v / 2;
    filter.displayContourPixelsXFIG(fc, color, colorfill, thickness, h, v);        
  }

  if ( args.check( "-displayPixelSegment" ) )
  {
    int color = args.getOption( "-displayPixelSegment" )->getIntValue( 0 );
    int fillcolor = args.getOption( "-displayPixelSegment" )->getIntValue( 1 );
    int depth = args.getOption( "-displayPixelSegment" )->getIntValue( 2 );
    int nb = VectOriginalSegmentMLPhv.size();
    float pixelx, pixely;
    pixelx = VectGlobalSegmentMLP[0].first_point.x();
    pixely = VectGlobalSegmentMLP[0].first_point.y();
    filter.displayPixel( pixelx, pixely, color, fillcolor, 1, 1, depth );
    for( int i = 0; i < nb ; ++i )
    {
      FreemanChain fc;
      fc = VectGlobalSegmentMLP[i].fc;
      int nbfc = fc.chain.size();
      for( int j = 0; j < nbfc ; ++j )
      {
        if( fc.code( j ) == 0 ) pixelx += 1;
        else  if( fc.code( j ) == 2 ) pixelx -= 1;
        else  if( fc.code( j ) == 1 ) pixely -= 1;
        else pixely += 1;
      }
      filter.displayPixel( pixelx, pixely, color, fillcolor, 1, 1, depth );
    }  
  }

  if ( args.check( "-displayPixelSegmenthv" ) )
  {
    int color = args.getOption( "-displayPixelSegmenthv" )->getIntValue( 0 );
    int depth = args.getOption( "-displayPixelSegmenthv" )->getIntValue( 1 );
    int nb = VectOriginalSegmentMLPhv.size();
    int x0 = VectGlobalSegmentMLP[0].first_point.x();
    int y0 = VectGlobalSegmentMLP[0].first_point.y();
    float pixelx = x0 - x0 % h + (float) h / 2 - 0.5;
    float pixely = y0 + y0 % v - ( v - 1 ) - 0.5 + (float) v / 2;
    filter.displayCroix( pixelx, pixely, color, h, v, depth ); 
    for( int i = 0; i < nb ; ++i )
    {
      FreemanChain fc;
      fc = VectOriginalSegmentMLPhv[i].fc;
      int nbfc = fc.chain.size();
      for( int j = 0; j < nbfc ; ++j )
      {
        if( fc.code( j ) == 0 ) pixelx += h;
        else  if( fc.code( j ) == 2 ) pixelx -= h;
        else  if( fc.code( j ) == 1 ) pixely -= v;
        else pixely += v;
      }
      filter.displayCroix( pixelx, pixely, color, h, v, depth );        
    }  
  }

//  if ( args.check( "-Digitization_hv" ) )
//  {
    //FreemanChain New_outer_chain;
    //Vector2D first_point( 0, 0 );
    //NewFreemanChainhv( New_outer_chain, first_point, c, h, v, x0, y0 );
    
    //FreemanChain inner_chain;
    //vector<uint> outer2inner;
    //vector<uint> inner2outer;
    //FreemanChain::innerContour( inner_chain, outer2inner, inner2outer, c, true);
    //FreemanChain New_outer_chain;
    //Vector2D first_point( VectGlobalSegmentMLP[0].first_point.x() , VectGlobalSegmentMLP[0].first_point.y() );
    //NewFreemanChainhv( New_outer_chain, first_point, inner_chain, h, v, x0, y0 );
    //FreemanChain New_inner_chain;
    //vector<uint> ic2nic;
    //vector<uint> nic2ic;
    //FreemanChain::subsample( New_inner_chain, ic2nic, nic2ic, inner_chain, h, v, 0, 0);
    //FreemanChain New_inner_chain1;
    //vector<uint> New_outer2inner;
    //vector<uint> New_inner2outer;
    //FreemanChain::New_Contour( New_inner_chain1, New_inner2outer, New_outer2inner, New_inner_chain, true  );
    //cout << " New_inner_chain1: " << New_inner_chain1.chain << endl;
    //FreemanChain New_outer_chain;
    //vector<uint> nic2noc;
    //vector<uint> noc2nic;
    //FreemanChain::outerContour( New_outer_chain, nic2noc, noc2nic, New_inner_chain1, true);
    //cout << " New_outer_chain: " << New_outer_chain.chain << endl;
 
  //  int nb = New_outer_chain.chain.size();
  //  cout << "nb: " << nb << endl;
  //  vector<double> Cx;
  //  vector<double> Cy;
  //  FreemanChain::const_iterator it = c.begin();
  //  Vector2i fxy( it.get() );
 //   Cx.push_back( fxy.x() );
 //   Cy.push_back( fxy.y() );
 //   int valx = fxy.x();
 //   int valy = fxy.y();
 //   for( int i = 0 ; i < nb ; ++i )
 //   {
 //     if( New_outer_chain.code ( i ) == 0 ) { valx = valx + 1; valy = valy; }
 //     else if( New_outer_chain.code ( i ) == 1 ) { valx = valx; valy = valy + 1; }
 //     else if( New_outer_chain.code ( i ) == 2 ) { valx = valx - 1; valy = valy; }
 //     else { valx = valx; valy = valy - 1; }
 //     Cx.push_back( valx );
 //     Cy.push_back( valy );
 //   }
 //   Cx.push_back( Cx.at(0) );
 //   Cy.push_back( Cy.at(0) );
 //   filter.outputPolylineInView( Cx, Cy, 
 //  	                 	   args.getOption( "-ContourSampled" )->getIntValue( 0 ),
//				   args.getOption( "-ContourSampled" )->getIntValue( 1 ),
//				   0, 40 );
//  }

if ( args.check( "-DisplayContourSurfelSampled" ) )
{
  int color = args.getOption( "-DisplayContourSurfelSampled" )->getIntValue( 0 );
  int thickness = args.getOption( "-DisplayContourSurfelSampled" )->getIntValue( 1 );
  int nb = VectOriginalSegmentMLPhv.size();
  FreemanChain fc;
  int convex_back, convex_front;
  Vector2i first_point;
  Vector2D xy;
  vector<double> CS_x;
  vector<double> CS_y;
  int x0 = VectGlobalSegmentMLP[0].first_point.x();
  int y0 = VectGlobalSegmentMLP[0].first_point.y();
  for( int i = 0 ; i < nb ; ++i )
  {
    convex_back = VectOriginalSegmentMLPhv[i].convex_back;
    convex_front = VectOriginalSegmentMLPhv[i].convex_front;
    fc.chain = VectOriginalSegmentMLPhv[i].fc.chain;
    if(fc.chain == "" ) continue;
    if( convex_back ) { xy.x() = x0 - 0.5; xy.y() = y0 + 0.5; }
    else { xy.x() = x0 - 0.5; xy.y() = y0 - 0.5; }
    CS_x.push_back( xy.x() );
    CS_y.push_back( xy.y() );
    fc.chain = VectOriginalSegmentMLPhv[i].fc.chain;
    cout << "freeman: " << fc.chain << endl;
    cout << "CSx: " << xy.x() << "  " << "CSy: " << xy.y() << endl;
    int nb_fc = fc.chain.size();
    for( int j = 0; j < nb_fc; ++j )
    {
      if( fc.code ( j ) == 0 ) { xy.x() += h; x0 += h; }
      else if( fc.code ( j ) == 1 ) { xy.y() -= v; y0 -= v; }
      else if( fc.code ( j ) == 2 ) { xy.x() -= h; x0 -= h; }
      else { xy.y() += v; y0 += v; }
      cout << "fc.code: " << fc.code(j) << "  xy: " << xy.x() << "  " << xy.y() << endl;
      //cout << "x0y0: " << x0 << "  " << y0 << endl;
      CS_x.push_back( xy.x() );
      CS_y.push_back( xy.y() );
    }
    int a = VectOriginalSegmentMLPhv[i].a;
    int b = VectOriginalSegmentMLPhv[i].b;
    if( convex_back && convex_front && ( a == 0 || b == 0 ) ) 
    {
      if( fc.code ( nb_fc - 1 ) == 0 ) { xy.x() += h; x0 += h; }
      else if( fc.code ( nb_fc - 1  ) == 1 ) { xy.y() -= v; y0 -= v; }
      else if( fc.code ( nb_fc - 1  ) == 2 ) { xy.x() -= h; x0 -= h; }
      else { xy.y() += v; y0 += v; }
      CS_x.push_back( xy.x() );
      CS_y.push_back( xy.y() );
    }
  } 
  CS_x.push_back( CS_x.at(0) );
  CS_y.push_back( CS_y.at(0) );
  filter.outputPolylineInView( CS_x, CS_y, color, thickness, 0, 40 );
}

if( args.check("-displayabmu") )
{  
  int color = args.getOption( "-displayabmu" )->getIntValue( 0 );
  int font_size = args.getOption( "-displayabmu" )->getIntValue( 1 );
  int nb = VectOriginalSegmentMLPhv.size();
  ostringstream out_str;
  float angle;
  int x0 = VectGlobalSegmentMLP[0].first_point.x();
  int y0 = VectGlobalSegmentMLP[0].first_point.y();
  float pixelxi = x0 - x0 % h + (float) h / 2 - 0.5;
  float pixelyi = y0 + y0 % v - ( v - 1 ) - 0.5 + (float) v / 2;
  float pixelxf = pixelxi;
  float pixelyf = pixelyi;
  for( int i = 0 ; i < nb ; ++i )
  {
    FreemanChain fc;
    fc = VectOriginalSegmentMLPhv[i].fc;
    if ( fc.chain != "")
    { 
      ostringstream out_str;
      out_str << "(" << VectOriginalSegmentMLPhv[i].a << "," << VectOriginalSegmentMLPhv[i].b << "," 
              << VectOriginalSegmentMLPhv[i].mu << "," << VectOriginalSegmentMLPhv[i].nbpointtested << ") " << endl; 
      double a = (double) VectOriginalSegmentMLPhv[i].a;
      double b = (double) VectOriginalSegmentMLPhv[i].b;
      if( b == 0 )  angle = atan2( b , a );
      else angle = atan2( a , b );
      int nbfc = fc.chain.size();
      for( int j = 0; j < nbfc ; ++j )
      {
        if( fc.code( j ) == 0 ) pixelxf += h;
        else  if( fc.code( j ) == 2 ) pixelxf -= h;
        else  if( fc.code( j ) == 1 ) pixelyf -= v;
        else pixelyf += v;
      }
      float pixelx = ( pixelxi + pixelxf ) / 2.0;
      float pixely = ( pixelyi + pixelyf ) / 2.0;
      if( ( a > 0 && b < 0 ) || ( a < 0 && b == 0 ) ) { pixelx += h / 2 + 1;}
      else if( a > 0 && b > 0 ){ pixelx += h / 2 + 1;  pixely += v; }
      else if( a == 0 && b < 0 ) { pixely -= v; }
      else if( a == 0 && b > 0 ) { pixely += v; }
      else if( ( a < 0 && b > 0 ) ) { pixelx -= h; pixely += 1.8; }
      else if( ( a < 0 && b < 0 ) ) { pixelx -= h/2; pixely -= 1; }
      filter.displayText(pixelx, pixely, out_str.str(), -frame.angleToX( angle ), color, font_size);
      pixelxi = pixelxf;
      pixelyi = pixelyf;
    }
  }
}

//if ( args.check( "-pointel_grid" ) )
//{
//  uint div_x = (uint) floor( xmax - xmin );
//  uint div_y = (uint) floor( ymax - ymin );
//  filter.outputGrid( xmin, ymin, xmax, ymax, div_x, div_y,
//                	 args.getOption( "-pointel_grid" )->getIntValue( 0 ),
//			 args.getOption( "-pointel_grid" )->getIntValue( 1 ),
//			 0, 60 ); 
//}
//
//if ( args.check( "-pixel_grid" ) )
//{
//  uint div_x = (uint) floor( xmax - xmin ) - 1;
//  uint div_y = (uint) floor( ymax - ymin ) - 1;
//  filter.outputGrid( xmin, ymin, xmax, ymax,
//             		 div_x, div_y,
//			 args.getOption( "-pixel_grid" )->getIntValue( 0 ),
//			 args.getOption( "-pixel_grid" )->getIntValue( 1 ),
//			 1, 70 ); 
//}
}
