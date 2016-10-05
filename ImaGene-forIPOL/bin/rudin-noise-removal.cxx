///////////////////////////////////////////////////////////////////////////////
// Generates contours from pgm image
///////////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <iomanip>

#include "ImaGene/base/Arguments.h"
#include "ImaGene/base/StandardArguments.h"

using namespace std;
using namespace ImaGene;


// double abs( double x )
// {
//   return x >= 0.0 ? x : -x ;
// }

double min( double x, double y )
{
  return x < y ? x : y ;
}
double max( double x, double y )
{
  return x > y ? x : y ;
}
double sqr( double x )
{
  return x * x;
}

double minmod( double a, double b )
{
  return 
    ( ( ( a > 0.0 ) ? 1.0 : ( ( a < 0.0 ) ? -1.0 : 0.0 ) )
      + ( ( b > 0.0 ) ? 1.0 : ( ( b < 0.0 ) ? -1.0 : 0.0 ) ) )
    * min( abs( a ), abs( b ) )  
    * 0.5; 
}

template <class PixelType>
class Image2D 
{
protected:
  uint m_w;
  uint m_h;
  PixelType* m_data;
public:
  Image2D() : m_data( 0 ) {}
  Image2D( uint w, uint h ) 
  : m_w( w ), m_h( h )
  {
    m_data = new PixelType[ m_w * m_h ];
  }
  Image2D( const Image2D & other ) 
    : m_w( other.m_w ), m_h( other.m_h )
  {
    m_data = new PixelType[ m_w * m_h ];
    copy( other.m_data, other.m_data + m_w * m_h, this->m_data );
  }
  Image2D<PixelType>& operator=( const Image2D<PixelType> & other ) 
  {
    if ( this != &other )
      {
	resize( other.m_w, other.m_h );
	copy( other.m_data, other.m_data + m_w * m_h, this->m_data );
      }
    return *this;

  }
  ~Image2D() 
  {
    if ( m_data != 0 ) delete[] m_data;
  }
  
  void exportPGM( Image2D<unsigned char> & other )
  {
    other.resize( w(), h() );
    for ( uint y = 0; y < h(); ++y )
      for ( uint x = 0; x < w(); ++x )
	{
	  PixelType v = getPixel( x, y );
	  if ( v < (PixelType) 0 ) v = (PixelType) 0;
	  else if ( v > (PixelType) 255 ) v = (PixelType) 255;
	  other.setPixel( x, y, (unsigned char) v ); 
	}
  }

  uint w() const { return m_w; }
  uint h() const { return m_h; }
  void resize( uint w, uint h )
  {
    if ( ( m_data != 0 )
	 && ( m_w == w )
	 && ( m_h == h ) ) return;
    if ( m_data != 0 ) delete[] m_data;
    m_w = w;
    m_h = h;
    m_data = new PixelType[ m_w * m_h ];
    
  }
  void setPixel( uint x, uint y, PixelType v )
  {
    m_data[ x + y*m_w ] = v;
  }
  PixelType getPixel( uint x, uint y ) const
  {
    return m_data[ x + y*m_w ];
  }

  void reset( PixelType v )
  {
    PixelType* src_b = m_data;
    PixelType* src_e = m_data + w() + h();
    while ( src_b != src_e )
      {
	*src_b = v;
	++src_b;
      }
  }

private:
  void copy( PixelType* src_b, PixelType* src_e, PixelType* dst )
  {
    while ( src_b != src_e )
      {
	*dst = *src_b;
	++dst; ++src_b;
      }
  }
};


class Grille2D : public Image2D<double>
{
public:
  Grille2D( uint w, uint h )
    : Image2D<double>( w, h )
  {}

  Grille2D( const Grille2D & img )
    : Image2D<double>( img )
  {}

  Grille2D( const Image2D<double> & img )
    : Image2D<double>( img )
  {}

  Grille2D( const Image2D<unsigned char> & img )
    : Image2D<double>( img.w(), img.h() )
  {
    for ( uint y = 0; y < img.h(); ++y )
      for ( uint x = 0; x < img.w(); ++x )
	setPixel( x, y, (double) img.getPixel( x, y ) );

  }
  Grille2D& operator=( const Grille2D & other )
  {
    Image2D<double>::operator=( other );
    return *this;
  }
  double u( uint x, uint y ) const
  {
    return getPixel( x, y );
  }
  
  double dxp( uint x, uint y ) const
  {
    return ( x < ( w() - 1 ) )
      ? u( x + 1, y ) - u( x, y )
      : 0.0;
  }
  double dxm( uint x, uint y ) const
  {
    return ( x > 0 )
      ? u( x, y ) - u( x - 1, y )
      : 0.0;
  }
  double dyp( uint x, uint y ) const
  {
    return ( y < ( h() - 1 ) )
      ? u( x, y + 1 ) - u( x, y )
      : 0.0;
  }
  double dym( uint x, uint y ) const
  {
    return ( y > 0 )
      ? u( x, y ) - u( x, y - 1 )
      : 0.0;
  }

  double normGradpp( uint x, uint y ) const
  {
    return sqrt( sqr( dxp( x, y ) ) + sqr( dyp( x, y ) ) );
  }

  double dx0( uint x, uint y ) const
  {
    return ( ( x > 0 ) && ( x < ( w() - 1 ) ) )
      ? ( u( x + 1, y ) - u( x - 1, y ) ) * 0.5
      : 0.0;
  }

  double dy0( uint x, uint y ) const
  {
    return ( ( y > 0 ) && ( y < ( h() - 1 ) ) )
      ? ( u( x, y + 1 ) - u( x, y - 1 ) ) * 0.5
      : 0.0;
  }
  
};



/**
 * Imports a PGM file from the stream
 * [in]. Creates an image.
 *
 * @param in (read) the input stream.
 * @param refImage (output) the created image.
 * @return 'true' if everything went well, otherwise return 'false'.
 */
template <class Image2D>
bool 
readPGM( std::istream & in, 
	 Image2D & refImage )
{
  string str;
  getline( in, str );
  if ( ! in.good() ) return false;
  if ( str != "P5" ) return false;
  do
    {
      getline( in, str );
      if ( ! in.good() ) return false;
    }
  while ( str[ 0 ] == '#' || str=="");
  istringstream str_in( str );
  uint sizes[ 2 ];
  str_in >> sizes[ 0 ] >> sizes[ 1 ];
  
  getline( in, str );
  istringstream str2_in( str );
  int max_value;
  str2_in >> max_value;

  if ( ! in.good() ) return false;
  cerr << "# PGM " << sizes[ 0 ] << " " << sizes[ 1 ] 
       << " " << max_value << " from <" << str << ">" << endl;

  refImage.resize( sizes[ 0 ], sizes[ 1 ] );

  uint nb_read = 0;
  in >> noskipws;
  for ( uint y = 0; y < sizes[ 1 ]; ++y )
    for ( uint x = 0; x < sizes[ 0 ]; ++x )
      { //... whatever
	unsigned char c; 
	in >> c;
	if ( in.good() ) ++nb_read;
	refImage.setPixel( x, y, c );
    }
  if ( in.fail() || in.bad() )
    {
      cerr << "# nbread=" << nb_read << endl;
      return false;
    }
  in >> skipws;
  return true;
}


/**
 * Exports a PGM file in the stream
 * [out], given some image.
 *
 * @param out (written) the output stream.
 * @param refImage the image to save.
 * @return 'true' if everything went well, otherwise return 'false'.
 */
template <class Image2D>
bool 
writePGM( std::ostream & out, 
	  const Image2D & image )
{
  string str;
  out << "P5" << endl;
  out << "# CREATOR: rudin-noise-removal " 
      << "(jacques-olivier.lachaud@univ-savoie.fr)." << endl;
  out << image.w() << " " << image.h() << endl;
  out << 255 << endl;

  uint nb_written = 0;
  for ( uint y = 0; y < image.h(); ++y )
    for ( uint x = 0; x < image.w(); ++x )
      { //... whatever
	unsigned char c = image.getPixel( x, y ); 
	out << c;
	if ( out.good() ) ++nb_written;
    }
  if ( out.fail() || out.bad() )
    {
      cerr << "# nbwritten=" << nb_written << endl;
      return false;
    }
  out << endl;
  return true;
}

double
randUniforme (double min, double max){
  double r;
  r = (double) rand() / (double) RAND_MAX;
  return ( r*(max-min) + min);
}

double
randomGaussien (double moyenne, double ecartType){
  double r1 = randUniforme (0.0, 1.0);
  double r2 = randUniforme (0.0, 1.0);
  double r = sqrt (-2.0 * log(r1)) * cos (2.0 * M_PI * r2);
  return (moyenne + ecartType * r);
}


void perturbate( Grille2D & u, double sigma )
{
  for ( uint j = 0; j < u.h(); ++j )
    for ( uint i = 0; i < u.w(); ++i )
      u.setPixel( i, j, u.u( i, j ) + randomGaussien( 0.0, sigma ) );
}

void rudinByVese( const Grille2D & u0, Grille2D & uf,
		  double dt, double dh, uint nb, double lambda )
{
  cerr << "--- Rudin et al. By Vese (" << u0.w() << "," << u0.h() << ")" 
       << " dt=" << dt << " dh=" << dh << " lambda=" << lambda << endl;
  
  Grille2D un( u0 );
  uint w = u0.w();
  uint h = u0.h();
  uint wm = w - 1;
  uint hm = h - 1;
  Grille2D Gxu( w, h );
  Grille2D Gyu( w, h );
  Gxu.reset( 0.0 );
  Gyu.reset( 0.0 );
  double eps = 1e-6;
  double eps2 = sqr( eps );
  Grille2D unn( u0 );
  for ( uint n = 0; n < nb; ++n )
    {
      cerr << "--- n = " << n;
      cerr << ", lambda = " << lambda;
      // Boundary conditions
      for ( uint j = 1; j < hm; ++j )
	{
	  un.setPixel( 0, j, un.u( 1, j ) );
	  un.setPixel( wm, j, un.u( wm - 1, j ) );
	}
      for ( uint i = 1; i < wm; ++i )
	{
	  un.setPixel( i, 0, un.u( i, 1 ) );
	  un.setPixel( i, hm, un.u( i, hm - 1 ) );
	}
      un.setPixel( 0, 0, un.u( 1, 1 ) );
      un.setPixel( wm, 0, un.u( wm - 1, 1 ) );
      un.setPixel( 0, hm, un.u( 1, hm - 1 ) );
      un.setPixel( wm, hm, un.u( wm - 1, hm - 1 ) );
      // Compute u_n+1
      double v;
      double dv;
      double d = 1 / ( 2.0 * lambda * sqr( dh ) );
      for ( uint j = 1; j < hm; ++j )
	for ( uint i = 1; i < wm; ++i )
	  {
	    double c1 = 1.0 /
	      sqrt( eps2 + sqr( un.dxp( i, j ) / dh ) 
		    + sqr( un.dy0( i, j ) / dh ) );
	    double c2 = 1.0 /
	      sqrt( eps2 + sqr( un.dxm( i, j ) / dh ) 
		    + sqr( un.dy0( i - 1, j ) / dh ) );
	    double c3 = 1.0 /
	      sqrt( eps2 + sqr( un.dx0( i, j ) / dh ) 
		    + sqr( un.dyp( i, j ) / dh ) );
	    double c4 = 1.0 /
	      sqrt( eps2 + sqr( un.dx0( i, j - 1 ) / dh ) 
		    + sqr( un.dym( i, j ) / dh ) );
	    v = 1.0 / ( 1.0 + d * ( c1 + c2 + c3 + c4 ) )
	      * ( u0.u( i, j ) 
		  + d * ( c1 * un.u( i + 1, j ) + c2 * un.u( i - 1, j ) 
			  + c3 * un.u( i, j + 1 ) + c4 * un.u( i, j - 1 ) ) ); 
	    unn.setPixel( i, j, v );
	  }
      double sum_u = 0.0;
      double sum_u2 = 0.0;
      for ( uint j = 0; j < h; ++j )
	for ( uint i = 0; i < w; ++i )
	  {
	    sum_u += un.u( i, j );
	    sum_u2 += sqr( un.u( i, j ) );
	  }
      cerr << ", E[u]=" << ( sum_u / (double) (w*h) )
	   << ", sigma[u]=" << sqrt( - sqr( sum_u / (double) (w*h) ) 
				     + sum_u2 / (double) (w*h) ) << endl;
      un = unn;
    }
  uf = un;
}

/**
 * Rudin, Osher and Fatami algorithm for denoising images. It is based
 * on the optimization of a functional which is a fitting in the
 * L1-sense constrained with some smoothness. The scheme itself is not
 * the original one (which does not work properly) but is based on
 * some scheme proposed by Vese. Strangely, the [lambda] coefficient
 * is now a given user-defined paramter while it was estimated from
 * the data in the original ROF. However, this trick makes the
 * denoising work well.
 *
 * @param u0 the input image
 * @param uf (returns) the denoised image
 * @param dh the grid step which is 1/max(width,height)
 * @param nb the number of iterations
 *
 * @param lambda the scale coefficient (small (< 0.5): smooth a lot
 * and removes noise, high (>2): keeps features but also noise).
 */
void rudinByVese2( const Grille2D & u0, Grille2D & uf,
		   double dh, uint nb, double lambda )
{
  cerr << "--- Rudin et al. By Vese 2 (" << u0.w() << "," << u0.h() << ")" 
       << " dh=" << dh << " lambda=" << lambda << endl;
  cerr << "    [";
  
  Grille2D un( u0 );
  uint w = u0.w();
  uint h = u0.h();
  uint wm = w - 1;
  uint hm = h - 1;
  double eps = 1e-6;
  double eps2 = sqr( eps );
  double eps2dh2 = sqr( eps * dh );
  double d = 1 / ( 2.0 * lambda * sqr( dh ) );
  Grille2D unn( u0 );
  for ( uint n = 0; n < nb; ++n )
    {
      if ( n % ( nb / 20 ) == 0 )
	cerr << "." << flush;
      // Boundary conditions
      for ( uint j = 1; j < hm; ++j )
	{
	  un.setPixel( 0, j, un.u( 1, j ) );
	  un.setPixel( wm, j, un.u( wm - 1, j ) );
	}
      for ( uint i = 1; i < wm; ++i )
	{
	  un.setPixel( i, 0, un.u( i, 1 ) );
	  un.setPixel( i, hm, un.u( i, hm - 1 ) );
	}
      un.setPixel( 0, 0, un.u( 1, 1 ) );
      un.setPixel( wm, 0, un.u( wm - 1, 1 ) );
      un.setPixel( 0, hm, un.u( 1, hm - 1 ) );
      un.setPixel( wm, hm, un.u( wm - 1, hm - 1 ) );
      // Compute u_n+1
      double v;
      double dv;
      for ( uint j = 1; j < hm; ++j )
	for ( uint i = 1; i < wm; ++i )
	  {
	    double c1 = dh /
	      sqrt( eps2dh2 + sqr( un.dxp( i, j ) ) 
		    + sqr( un.dy0( i, j ) ) );
	    double c2 = dh /
	      sqrt( eps2dh2 + sqr( un.dxm( i, j ) ) 
		    + sqr( un.dy0( i - 1, j ) ) );
	    double c3 = dh /
	      sqrt( eps2dh2 + sqr( un.dx0( i, j ) ) 
		    + sqr( un.dyp( i, j ) ) );
	    double c4 = dh /
	      sqrt( eps2dh2 + sqr( un.dx0( i, j - 1 ) ) 
		    + sqr( un.dym( i, j ) ) );
	    v = 1.0 / ( 1.0 + d * ( c1 + c2 + c3 + c4 ) )
	      * ( u0.u( i, j ) 
		  + d * ( c1 * un.u( i + 1, j ) + c2 * un.u( i - 1, j ) 
			  + c3 * un.u( i, j + 1 ) + c4 * un.u( i, j - 1 ) ) ); 
	    unn.setPixel( i, j, v );
	  }
      // double sum_u = 0.0;
      // double sum_u2 = 0.0;
      // for ( uint j = 0; j < h; ++j )
      // 	for ( uint i = 0; i < w; ++i )
      // 	  {
      // 	    sum_u += un.u( i, j );
      // 	    sum_u2 += sqr( un.u( i, j ) );
      // 	  }
      // cerr << ", E[u]=" << ( sum_u / (double) (w*h) )
      // 	   << ", sigma[u]=" << sqrt( - sqr( sum_u / (double) (w*h) ) 
      // 				     + sum_u2 / (double) (w*h) ) << endl;
      un = unn;
    }
  uf = un;
  cerr << "]" << endl;
}

void rudin( const Grille2D & u0, Grille2D & uf,
	    double dt, double dh, uint nb, double sigma )
{
  cerr << "--- Rudin et al. (" << u0.w() << "," << u0.h() << ")" 
       << " dt=" << dt << " dh=" << dh << " sigma=" << sigma << endl;
  
  Grille2D un( u0 );
  perturbate( un, sigma );
  uint w = u0.w();
  uint h = u0.h();
  uint wm = w - 1;
  uint hm = h - 1;
  Grille2D Gxu( w, h );
  Grille2D Gyu( w, h );
  Gxu.reset( 0.0 );
  Gyu.reset( 0.0 );
  double eps = 1e-6;
  for ( uint n = 0; n < nb; ++n )
    {
      cerr << "--- n = " << n;
      // compute lambda
      double lambda = 0.0;
      for ( uint j = 1; j < hm; ++j )
	for ( uint i = 1; i < wm; ++i )
	  {
	    double ng = un.normGradpp( i, j );  
	    // cerr << " " << ng ;
	    if ( ng > 0.0 )
	      lambda +=  ng - ( ( u0.dxp( i, j ) * un.dxp( i, j )
				  + u0.dyp( i, j ) * un.dyp( i, j ) ) 
				/ ng );
	  }
      lambda *= - dh / ( 2.0 * sqr( sigma ) );

      cerr << ", lambda = " << lambda;
      // Compute Gxu and Gyu
      double v;
      double dv;
      for ( uint j = 1; j < hm; ++j )
	for ( uint i = 1; i < wm; ++i )
	  {
	    dv = sqrt( sqr( un.dxp( i, j ) ) 
		      + sqr( minmod( un.dyp( i,j ), un.dym( i,j ) ) ) );
	    v = ( dv > eps ) ? ( un.dxp( i, j ) / dv ) : 0.0;
	    Gxu.setPixel( i, j, v );

	    dv = sqrt( sqr( un.dyp( i, j ) ) 
		       + sqr( minmod( un.dxp( i,j ), un.dxm( i,j ) ) ) );
	    v = ( dv > eps ) ? ( un.dyp( i, j ) / dv ) : 0.0;
	    Gyu.setPixel( i, j, v );
	  }

      // Boundary conditions
      for ( uint j = 1; j < hm; ++j )
	{
	  un.setPixel( 0, j, un.u( 1, j ) );
	  un.setPixel( wm, j, un.u( wm - 1, j ) );
	}
      for ( uint i = 1; i < wm; ++i )
	{
	  un.setPixel( i, 0, un.u( i, 1) );
	  un.setPixel( i, hm, un.u( i, hm - 1 ) );
	}
      // Compute u_n+1
      double diff = 0.0;
      double dt_dh = dt / dh;
      for ( uint j = 1; j < hm; ++j )
	for ( uint i = 1; i < wm; ++i )
	  {
	    v = dt_dh * ( Gxu.dxm( i, j ) + Gyu.dym( i, j ) )
	      - dt * lambda * ( un.u( i, j ) - u0.u( i, j ) );
	    diff += abs( v );
	    un.setPixel( i, j, un.u( i, j ) + v );
	  }      
      double sum_u = 0.0;
      double sum_u2 = 0.0;
      for ( uint j = 0; j < h; ++j )
	for ( uint i = 0; i < w; ++i )
	  {
	    sum_u += un.u( i, j );
	    sum_u2 += sqr( un.u( i, j ) );
	  }
      cerr << ", diff=" << ( diff / (double) (w*h) );
      cerr << ", E[u]=" << ( sum_u / (double) (w*h) )
	   << ", sigma[u]=" << sqrt( - sqr( sum_u / (double) (w*h) ) 
				     + sum_u2 / (double) (w*h) ) << endl;
    }
  uf = un;
}


static Arguments args;

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
  StandardArguments::addIOArgs( args, true, true );
  args.addOption( "-sigma", "-sigma <val>: standard deviation of input image noise.", "1.0" );
  args.addOption( "-c", "-c <val>: the constant in the CFL criterion.", "0.5" );
  args.addOption( "-nb", "-nb <n>: the number of iterations.", "10" );
  args.addOption( "-scheme", "-scheme <n>: 0: Rudin, 1: Vese, 2: Vese opt.", "2" );

  if ( ( argc <= 0 ) 
       || ! args.readArguments( argc, argv ) ) 
    {
      cerr << args.usage( "rudin-noise-removal", 
			  "Remove noise from a PGM image given on the standard input and writes the restored image on the standard output as a PGM image.",
			  "" )
	   << endl;
      return 1;
    }

  double sigma = args.getOption( "-sigma" )->getDoubleValue( 0 );
  double c = args.getOption( "-c" )->getDoubleValue( 0 );
  uint nb = args.getOption( "-nb" )->getIntValue( 0 );

  istream & in_str = StandardArguments::openInput( args );
  Image2D<unsigned char> img;
  if ( readPGM( in_str, img ) )
    {
      cerr << "Success loading image." << endl;
    }
  else 
    cerr << "Error loading image." << endl;


  Grille2D u0( img );
  Grille2D uf( img.w(), img.h() );
  double dh = 1.0 / (double) max( img.w(), img.h() );
  double dt = c * sqr( dh );
  int scheme = args.getOption( "-scheme" )->getIntValue( 0 );
  if ( scheme == 2 )
    rudinByVese2( u0, uf, dh, nb, sigma );
  else if ( scheme == 1 )
    rudinByVese( u0, uf, dt, dh, nb, sigma );
  else
    rudin( u0, uf, dt, dh, nb, sigma );

  uf.exportPGM( img );
  
//   for ( uint y = 0; y < img.h(); ++y )
//     for ( uint x = 0; x < img.w(); ++x )
//       img.setPixel( x, y, 255 - img.getPixel( x, y ) );



  ostream & out_str = StandardArguments::openOutput( args );
  if ( writePGM( out_str, img ) )
    {
      cerr << "Success writing image." << endl;
    }
  else 
    cerr << "Error writing image." << endl;
  
}




 
