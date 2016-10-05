///////////////////////////////////////////////////////////////////////////////
// Generates polygons from digital contours
///////////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <vector>

#include "ImaGene/base/Arguments.h"
#include "ImaGene/base/StandardArguments.h"
#include "ImaGene/base/Vector2i.h"
#include "ImaGene/base/Vector.h"
#include "ImaGene/mathutils/Line2D.h"
#include "ImaGene/timetools/Clock.h"
#include "ImaGene/dgeometry2d/FreemanChain.h"
#include "ImaGene/dgeometry2d/C4CTangentialCover.h"
#include "ImaGene/dgeometry2d/C4CIteratorOnFreemanChain.h"
#include "ImaGene/digitalnD/C4CIteratorOnFreemanChainSurface.h"
#include "ImaGene/helper/C4CTangentialCoverGeometry.h"
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

void computeCoordinates( vector<double> & x,
			 vector<double> & y, 
			 const CharacteristicPolygon & P )
{
  for ( uint i = 0; i < P.edges().size(); ++i )
    {
      Vector2D P1 = P.edges()[ i ].getP1();
      Vector2D P2 = P.edges()[ i ].getP2();
      if ( ( i == 0 )
	   || ( x.back() != P1.x() ) 
	   || ( y.back() != P1.y() ) )
	{
	  x.push_back( P1.x() );
	  y.push_back( P1.y() );
	}
    }
  if ( ( x.back() != x[ 0 ] ) 
       || ( y.back() != y[ 0 ] ) )
    {
      x.push_back( x[ 0 ] );
      y.push_back( y[ 0 ] );
    }
}

void computeMedianPolygon( vector<double> & x,
			   vector<double> & y, 
			   const CharacteristicPolygon & P )
{
  uint nbe = P.edges().size();
  // cerr << "# Compute median polygon with nbe=" << nbe << endl;
  Vector2D pi;
  Vector2D pj;
  // Start edge
  deque<Line2D> F;
  deque<Line2D> Fe;
  Fe.push_front( P.edges()[ 0 ].medianLine() );
  for ( uint i = 1; i < nbe; ++i )
    {
      Line2D cur_line = P.edges()[ i ].medianLine();
      if ( cur_line.intersectionPoint( Fe.back(), pi ) )
	Fe.push_back( cur_line );
    }
  while ( ! Fe.back().intersectionPoint( Fe.front(), pi ) )
    Fe.pop_back();
  nbe = Fe.size();
  // cerr << "# Only " << nbe << " different directions." << endl;

  F.push_front( Fe.front() );
  Fe.push_back( Fe.front() );
  Fe.pop_front();
  uint i = 1;
  while ( i <= nbe )
    {
      Line2D & prev_line = F.back();
      if ( ! prev_line.intersectionPoint( Fe.front(), pi ) )
	cerr << "# Erreur 1 !" << endl;
      Line2D cur_line = Fe.front();
      Fe.push_back( Fe.front() );
      Fe.pop_front();
      ++i;
      float t = prev_line.position( pi );
      while ( true )
	{
	  if ( ! prev_line.intersectionPoint( Fe.front(), pj ) )
	    break;
	  float t2 = prev_line.position( pj );
	  if ( t <= t2 ) break;
	  t = t2;
	  pi = pj;
	  cur_line = Fe.front();
	  Fe.push_back( Fe.front() );
	  Fe.pop_front();
	  ++i;
	}
      // pi is the intersection
      // cout << pi.x() << " " << pi.y() << endl;
      x.push_back( pi.x() );
      y.push_back( pi.y() );
      F.push_back( cur_line );
    }
  if ( ( x.back() != x[ 0 ] ) 
       || ( y.back() != y[ 0 ] ) )
    {
      x.push_back( x[ 0 ] );
      y.push_back( y[ 0 ] );
    }
//   uint i = 1;
//   while ( i < nbe )
//     {
//       cerr << "# edge " << i << endl;
//       Line2D & prev_line = F.back();
//       Vector2D pi;
//       Vector2D pj;
//       while ( ! prev_line.intersectionPoint
// 	      ( P.edges()[ i % nbe ].medianLine(), pi ) )
// 	{
// 	  cerr << P.edges()[ i % nbe ].medianLine().direction() << endl;
// 	  ++i;
// 	}
//       float t = prev_line.position( pi );
//       uint j = i+1;
//       while ( true )
// 	{
// 	  while ( ! prev_line.intersectionPoint
// 		  ( P.edges()[ j % nbe ].medianLine(), pj ) )
// 	    ++j;
// 	  float t2 = prev_line.position( pj );
// 	  if ( t <= t2 ) break;
// 	  t = t2;
// 	  pi = pj;
// 	  i = j;
// 	  ++j;
// 	}
//       // pi is the intersection
//       cout << pi.x() << " " << pi.y() << endl;
//       x.push_back( pi.x() );
//       y.push_back( pi.y() );
//       F.push_back( P.edges()[ i % nbe ].medianLine() );
//     }
}

void computeIntermediatePolygon( vector<double> & x,
				 vector<double> & y, 
				 const CharacteristicPolygon & P,
				 float alpha )
{
  uint nbe = P.edges().size();
  // cerr << "# Compute median polygon with nbe=" << nbe << endl;
  Vector2D pi;
  Vector2D pj;
  // Start edge
  deque<Line2D> F;
  deque<Line2D> Fe;
  Fe.push_front( P.edges()[ 0 ].intermediateLine( alpha ) );
  for ( uint i = 1; i < nbe; ++i )
    {
      Line2D cur_line = P.edges()[ i ].intermediateLine( alpha );
      if ( cur_line.intersectionPoint( Fe.back(), pi ) )
	Fe.push_back( cur_line );
    }
  while ( ! Fe.back().intersectionPoint( Fe.front(), pi ) )
    Fe.pop_back();
  nbe = Fe.size();
  // cerr << "# Only " << nbe << " different directions." << endl;

  F.push_front( Fe.front() );
  Fe.push_back( Fe.front() );
  Fe.pop_front();
  uint i = 1;
  while ( i <= nbe )
    {
      Line2D & prev_line = F.back();
      if ( ! prev_line.intersectionPoint( Fe.front(), pi ) )
	cerr << "# Erreur 1 !" << endl;
      Line2D cur_line = Fe.front();
      Fe.push_back( Fe.front() );
      Fe.pop_front();
      ++i;
      float t = prev_line.position( pi );
      while ( true )
	{
	  if ( ! prev_line.intersectionPoint( Fe.front(), pj ) )
	    break;
	  float t2 = prev_line.position( pj );
	  if ( t <= t2 ) break;
	  t = t2;
	  pi = pj;
	  cur_line = Fe.front();
	  Fe.push_back( Fe.front() );
	  Fe.pop_front();
	  ++i;
	}
      // pi is the intersection
      // cout << pi.x() << " " << pi.y() << endl;
      x.push_back( pi.x() );
      y.push_back( pi.y() );
      F.push_back( cur_line );
    }
  if ( ( x.back() != x[ 0 ] ) 
       || ( y.back() != y[ 0 ] ) )
    {
      x.push_back( x[ 0 ] );
      y.push_back( y[ 0 ] );
    }
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
  StandardArguments::addIOArgs( args, true, false );
  args.addOption( "-step", "-step <h>: the discretization step or resolution, the closer to 0, the finer.", "1.0" );
  args.addBooleanOption( "-dCP", "-dCP: displays characteristic polygon as a set of vertices in Z2+(1/2,1/2)." );
  args.addBooleanOption( "-dCMLP", "-dCMLP: displays combinatorial minimal length polygon as a set of vertices in Z2+(1/2,1/2)." );
  args.addBooleanOption( "-dCMLP2", "-dCMLP2: displays combinatorial minimal length polygon as a set of vertices in R2. Coordinates are correct wrt input digital contour in the half-integer plane." );
  args.addBooleanOption( "-dCPh", "-dCPh: displays characteristic polygon as a set of vertices embedded in h(Z+1/2) x h(Z+1/2)." );
  args.addBooleanOption( "-dMP", "-dMP: displays median polygon as a set of vertices in R2." );
  args.addOption( "-dIP", "-dIP <alpha>: displays intermediate polygon of coefficient <alpha> as a set of vertices in R2 (alpha = 0 => IP=CP, alpha = 0.5 => IP=MP)", "0.5" );
  args.addBooleanOption( "-auto_center", "-auto_center: automatically centers the contour in the space." );
  args.addBooleanOption( "-lenCP", "-lenCP: displays Euclidean length of characteristic polygon." );
  args.addBooleanOption( "-lenCMLP", "-lenCMLP: displays Euclidean length of combinatorial MLP." );
  args.addBooleanOption( "-lenMP", "-lenMP: displays Euclidean length of median polygon." );
  args.addOption( "-lenIP", "-lenIP <alpha>: displays Euclidean length of intermediate polygon of coefficient alpha.", "0.5" );
  args.addOption( "-repeat", "-repeat <N>: repeat computations N times. Useful for timings.", "1" );
  args.addBooleanOption( "-v", "-v: verbose mode." );
  args.addBooleanOption( "-timeCMLP", "-timeCMLP: displays the computation time of ``N`` CMLP.");
  args.addBooleanOption( "-timeCP", "-timeCP: displays the computation time of ``N`` CP.");
  
  if ( ( argc <= 0 ) 
       || ! args.readArguments( argc, argv ) ) 
    {
      cerr << args.usage( "freeman2polygon", 
			  "Converts a contour defined as a Freeman chain code to a polygon. It writes it on the standard output.",
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

  double dh = (double) args.getOption( "-step" )->getDoubleValue( 0 );

  // -------------------------------------------------------------------------
  // Build tangential cover and its geometry.

  C4CIteratorOnFreemanChain itfc;
  itfc.init( c.begin(), true );
  C4CIteratorOnFreemanChainSurface itfcs;
  itfcs.init( ks, itfc );

  // cerr << *ks << endl;
  // C4CIteratorOnFreemanChainSurface itfcs2( itfcs );
  // do {
  //   ks->displayKn_sid( itfcs2.current(), cerr );
  //   cerr << endl;
  //   itfcs2.next();
  // } while ( ! itfcs2.equals( itfcs ) );


  vector<double> x;
  vector<double> y;
  vector<int> vx;
  vector<int> vy;
  vector<uint> vi;
  vector<bool> vt;

  CharacteristicPolygon P;
  if (   args.check( "-lenMP" ) || args.check( "-dMP" ) 
      || args.check( "-lenCP" ) || args.check( "-dCPh" ) 
      || args.check( "-dCP"   ) || args.check( "-timeCP" ) ){
    C4CTangentialCover tcover;
    buildTangentialCover( tcover, itfcs, 0 );
    C4CTangentialCoverGeometry tcover_geometry;
    Frame2D frame;
    frame.init( ks, 0, 1 );
    C4CTangentialCoverGeometry::MSGeometryComputerByLP geocomp;
    tcover_geometry.init( tcover, itfcs, geocomp, frame );


    // -------------------------------------------------------------------------
    // Computes characteristic polygon.

    uint N = args.getOption( "-repeat" )->getIntValue( 0 );
    double test = 0.0;
    Clock::startClock();
    for (int i=0; i<N; i++) {
      P.init( tcover, tcover_geometry );
      computeCoordinates( x, y, P );
      test = 0.0;
      for ( uint i = 0; i < (x.size()-1); ++i )
        test += abs(x[i+1] - x[i]) + abs(y[i+1] - y[i]);
      x.clear();
      y.clear();
    }
    long t2 = Clock::stopClock();
    cerr << "# Tangential cover & geometry, characteristic polygon in " << t2 << " ms." << endl;
    if (args.check( "-timeCP")) {
      cerr << "# CP   test value : " << (int) test << "." <<  test  - (int) test  << endl;
      cout << t2 ;
    }
    if (args.check( "-timeCMLP"))
      cout << " ";
    else
      cout << endl;
  }

  cout << setprecision(15);
  if ( args.check( "-dCP" ) )
    {
      computeCoordinates( x, y, P );
      if ( args.check( "-v" ) )
	{
	  cout << "# Half-integer characteristic polygon " << endl;
	  cout << "# nbv=" << x.size() - 1 << endl;
	}
      for ( uint i = 0; i < x.size(); ++i )
	cout << x[ i ] << " " << y[ i ] << endl;
    }

  if ( args.check( "-dCMLP" ) )
    {
      FreemanChain::computeMLP( vx, vy, vi, c );
      double len = 0.0;
      uint n = vx.size();
      for ( uint i = 0; i < vx.size(); ++i )
	len += sqrt( (vx[(i+1)%n]-vx[i])*(vx[(i+1)%n]-vx[i])
		     + (vy[(i+1)%n]-vy[i])*(vy[(i+1)%n]-vy[i]) );
      if ( args.check( "-v" ) )
	{
	  cout << "# integer CMLP " << endl;
	  cout << "# nbv=" << vx.size() << endl;
	}
      for ( uint i = 0; i < vx.size(); ++i )
	cout << vx[ i ] << " " << vy[ i ] << " " << vi[ i ] << endl;
      cout << "# len=" << len << endl;
    }

  if ( args.check( "-dCMLP2" ) )
    {
      int nb_ccw_loops = c.isClosed();
      bool cw = nb_ccw_loops < 0;
      Vector2i twice_dv;
      FreemanChain::computeMLP( vx, vy, vi, vt, twice_dv, c, cw );
      double len = 0.0;
      uint n = vx.size();
      for ( uint i = 0; i < vx.size(); ++i )
	len += sqrt( (vx[(i+1)%n]-vx[i])*(vx[(i+1)%n]-vx[i])
		     + (vy[(i+1)%n]-vy[i])*(vy[(i+1)%n]-vy[i]) );
      if ( args.check( "-v" ) )
	{
	  cout << "# CMLP " << endl;
	  cout << "# nbv=" << vx.size() << endl;
	}
      for ( uint i = 0; i < vx.size(); ++i )
	cout << ( vx[ i ] + 0.5 * twice_dv.x() ) 
	     << " " << ( vy[ i ] + 0.5 * twice_dv.y() )
	     << " " << vi[ i ] << " " << vt[ i ] << endl;
      cout << "# len=" << len << endl;
    }

  if ( args.check( "-dCPh" ) )
    {
      computeCoordinates( x, y, P );
      if ( args.check( "-v" ) )
	{
	  cout << "# Embedded characteristic polygon " << endl;
	  cout << "# nbv=" << x.size() - 1 << endl;
	  cout << "# dh=" << dh << endl;
	}
      for ( uint i = 0; i < x.size(); ++i )
	cout << x[ i ]*dh << " " << y[ i ]*dh << endl;
    }

  if ( args.check( "-lenCP" ) )
  {
    computeCoordinates( x, y, P );
    double len = 0.0;
    for ( uint i = 0; i < (x.size()-1); ++i )
      len += sqrt( (x[i+1]-x[i])*(x[i+1]-x[i])
          + (y[i+1]-y[i])*(y[i+1]-y[i]) );
    if ( args.check( "-v" ) )
    {
      cout << "# Eucliden length of characteristic polygon " << endl;
      cout << "# nbv=" << x.size() - 1 << endl;
      cout << "# dh=" << dh << endl;
    }
    len *= dh;
    cout << len << endl;
  }

  if ( args.check( "-lenCMLP" ) )
  {
    uint N = args.getOption( "-repeat" )->getIntValue( 0 );
    Clock::startClock();
    double len;
    for ( uint i = 0; i < N; ++i )
      len = FreemanChain::lengthMLP( c );
    long t = Clock::stopClock();
    if ( args.check( "-v" ) )
    {
      cout << "# Euclidean length of combinatorial MLP" << endl;
      cout << "# dh=" << dh << endl;
      cout << "# CMLP computed " << N << " times in " << t << " ms." << endl;
    }
    len *= dh;
    cout << len << endl;
  }

  if ( args.check( "-timeCMLP" ) )
  {
    uint N = args.getOption( "-repeat" )->getIntValue( 0 );
    double test;
    Clock::startClock();
    for ( uint i = 0; i < N; ++i )
      test = FreemanChain::testCMLP( c );
    long t = Clock::stopClock();
    cerr << "# CMLP test value : " << (int) test << "." <<  test  - (int) test  << endl;
    cout << t << endl;
  }

  if ( args.check( "-dMP" ) )
  {
      computeMedianPolygon( x, y, P );
      if ( args.check( "-v" ) )
	{
	  cout << "# median polygon " << endl;
	  cout << "# nbv=" << x.size() - 1 << endl;
	}
      for ( uint i = 0; i < x.size(); ++i )
	cout << x[ i ] << " " << y[ i ] << endl;
    }

  if ( args.check( "-dIP" ) )
    {
      float alpha = args.getOption( "-dIP" )->getFloatValue( 0 );
      computeIntermediatePolygon( x, y, P, alpha );
      if ( args.check( "-v" ) )
	{
	  cout << "# intermediate polygon alpha=" << alpha << endl;
	  cout << "# nbv=" << x.size() - 1 << endl;
	}
      for ( uint i = 0; i < x.size(); ++i )
	cout << x[ i ] << " " << y[ i ] << endl;
    }

  if ( args.check( "-lenMP" ) )
    {
      computeMedianPolygon( x, y, P );
      double len = 0.0;
      for ( uint i = 0; i < (x.size()-1); ++i )
	len += sqrt( (x[i+1]-x[i])*(x[i+1]-x[i])
		     + (y[i+1]-y[i])*(y[i+1]-y[i]) );
      len *= dh;
      if ( args.check( "-v" ) )
	{
	  cout << "# Eucliden length of median polygon " << endl;
	  cout << "# nbv=" << x.size() - 1 << endl;
	  cout << "# dh=" << dh << endl;
	}
      cout << len << endl;
    }
  if ( args.check( "-lenIP" ) )
    {
      float alpha = args.getOption( "-lenIP" )->getFloatValue( 0 );
      computeIntermediatePolygon( x, y, P, alpha );
      double len = 0.0;
      for ( uint i = 0; i < (x.size()-1); ++i )
	len += sqrt( (x[i+1]-x[i])*(x[i+1]-x[i])
		     + (y[i+1]-y[i])*(y[i+1]-y[i]) );
      len *= dh;
      if ( args.check( "-v" ) )
	{
	  cout << "# Eucliden length of intermediate polygon alpha=" << alpha
	       << endl;
	  cout << "# nbv=" << x.size() - 1 << endl;
	  cout << "# dh=" << dh << endl;
	}
      cout << len << " " << alpha << endl;
    }


  return 0;
}
