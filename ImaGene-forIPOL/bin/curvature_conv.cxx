///////////////////////////////////////////////////////////////////////////////
// Computes curvature with estimation of circumcircle
///////////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>

#include "ImaGene/base/Arguments.h"
#include "ImaGene/base/Proxy.h"
#include "ImaGene/base/Vector.h"
#include "ImaGene/base/VectorUtils.h"
#include "ImaGene/base/Shapes.h"
#include "ImaGene/base/Signal.h"
#include "ImaGene/timetools/Clock.h"
#include "ImaGene/dgeometry2d/C4CTangentialCover.h"
#include "ImaGene/dgeometry2d/C4CIteratorOnFreemanChain.h"
#include "ImaGene/dgeometry2d/EuclideanGeometry.h"
#include "ImaGene/digitalnD/C4CIteratorOnFreemanChainSurface.h"
#include "ImaGene/digitalnD/KnSpace.h"
#include "ImaGene/digitalnD/KnCharSet.h"
#include "ImaGene/digitalnD/KnShapes.h"
#include "ImaGene/digitalnD/C4CIteratorOnSurface.h"
#include "ImaGene/digitalnD/BelAdjacency.h"
#include "ImaGene/digitalnD/ObjectBoundary.h"
#include "ImaGene/digitalnD/Frame2D.h"
#include "ImaGene/helper/C4CTangentialCoverGeometry.h"
#include "ImaGene/helper/GlobalC4CGeometry.h"
#include "ImaGene/helper/ContourHelper.h"
#include "ImaGene/helper/ShapeHelper.h"

using namespace std;
using namespace ImaGene;


static Arguments args;

void buildTangentialCover( C4CTangentialCover & tcover, 
			   C4CIterator & cp, uint max_size )
{
  tcover.init( cp, max_size );
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

  args.addOption( "-step", "-step <h>: the discretization step or resolution, the closer to 0, the finer.", "1.0" );
  args.addOption( "-mask_size", "-mask_size <n>: the size of the mask is 2*n+1. The size of the mask is related to the digitization step. If none is given, the mask size is ceil( 1/2*(h^(-4/3)) ).", "1" );
  
  if ( ( argc <= 0 ) 
       || ! args.readArguments( argc, argv ) ) 
    {
      cerr << args.usage( "curvature_conv", 
			  "Estimates the curvature along a digital contour given as a Freeman chaincode on the standard input. The estimator is based on successive convolution by (derivative of) binomial kernels [Brunet, Malgouyres 2008].   The size of the mask is related to the digitization step. If none is given, the mask size is ceil( 1/2*(h^(-4/3)) ).",
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

  KnSpace* ks = ShapeHelper::makeSpaceFromFreemanChain( c );
  C4CIteratorOnFreemanChain itfc;
  itfc.init( c.begin(), true );
  C4CIteratorOnFreemanChainSurface itfcs;
  itfcs.init( ks, itfc );
  double dh = (double) args.getOption( "-step" )->getDoubleValue( 0 );

  bool open;
  const int nb_bels = (int) C4CIterator::size( itfcs, open );
  cerr << "# Size boundary Nb=" << nb_bels << " open=" << open << endl;

  // ----------------------------------------------------------------------
  // The size of the mask is related to the digitization step. If none
  // is given, the mask size is ceil( 1/2*(h^(-4/3)) ).
  Clock::startClock();

  uint n = args.check( "-mask_size" ) 
    ? args.getOption( "-mask_size" )->getIntValue( 0 )
    : (uint) ceil( 0.5 / pow( dh, 4.0/3.0 ) ); 
  
  Signal<double> X( nb_bels, 0, true, 0.0 );
  Signal<double> Y( nb_bels, 0, true, 0.0 );

  // First, define the curve with the coordinates of the surfels.
  Proxy<C4CIteratorOnSurface> it( (C4CIteratorOnSurface*) itfcs.clone() );
  Frame2D frame;
  frame.init( ks, 0, 1 );
  for ( int i = 0; i < nb_bels; ++i )
    {
      Kn_sid sbel = it->current();
      frame.setSurfelFrame( sbel, it->trackDir() );
      Vector2D p( frame.transformPoint( Vector2D( 0.5, 0.0 ) ) );
      X[ i ] = p.x();
      Y[ i ] = p.y();
      // cerr << p.x() << " " << p.y() << endl;
      it->next();
    }

  Signal<double> G( Signal<double>::H2() );
  G.multiply( 0.25 );
  Signal<double> H2n = G;
  for ( uint i = 1; i < n; ++i )
    H2n = H2n * G; 

  Signal<double> HX( X * H2n );
  Signal<double> HY( Y * H2n );
  Signal<double> DX( HX * Signal<double>::Delta() );
  Signal<double> DY( HY * Signal<double>::Delta() );
  Signal<double> DDX( DX * Signal<double>::Delta() );
  Signal<double> DDY( DY * Signal<double>::Delta() );

  long t = Clock::stopClock();
  cerr << "# Curvature by convolution in " << t << " ms." << endl;
  cerr << "# TIME=" << t << endl;

  cout << setprecision(15);
  cout << "#######################################################" << endl;
  cout << "# Curvature by convolutions." << endl
       << "# Size boundary Nb=" << nb_bels << endl
       << "# dh=" << dh << endl
       << "# num curv_est tgt_angle_est cabs dabs el l xc yc xc_est yc_est dx_est dy_est" << endl
       << "# (el: elementary length of a surfel estimated from curvature)" << endl
       << "# (l: sum of elementary lengthes from initial surfel)" << endl;


  double dabs = 0.5;
  double cabs = 0.0;
  for ( int i = 0; i < nb_bels; ++i )
    {
      double tgt_angle_to_x = atan2( DY[ i ], DX[ i ] );
      double l_surfel = sqrt( Mathutils::sqr( HX[ i + 1 ] - HX[ i ] )
			      + Mathutils::sqr( HY[ i + 1 ] - HY[ i ] ) );
      double denom = pow( Mathutils::sqr( DX[ i ] )
			  + Mathutils::sqr( DY[ i ] ), 1.5 );
      double curv = ( denom != 0.0 )
	? ( DDX[ i ] * DY[ i ] - DDY[ i ] * DX[ i ] ) / denom
	: 0.0;
      cout << i << " " 
	   << curv / dh << " "
	   << tgt_angle_to_x << " "
	   << cabs << " "
	   << dabs << " "
	   << l_surfel << " "
	   << cabs << " "
	   << X[ i ] << " "
	   << Y[ i ] << " "
	   << HX[ i ] << " "
	   << HY[ i ] << " "
	   << DX[ i ] << " "
	   << DY[ i ] << " "
	   << endl;
      cabs += l_surfel;
      dabs += 0.5;
    }

  // Frees some stuff.
  delete ks;
  return 0;
}
