///////////////////////////////////////////////////////////////////////////////
// Computes curvature with global optimization method.
///////////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>

#include "ImaGene/base/Arguments.h"
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



C4CELength* computeGlobalGeometry( GlobalC4CGeometry & global_geom,
			    KnSpace* ks, 
			    C4CIteratorOnSurface & cp,
			    LinearMinimizer* lm )
{
  Clock::startClock();
  C4CTangentialCover tcover;
  buildTangentialCover( tcover, cp, 0 );
  long t = Clock::stopClock();
  cerr << "# Tangential cover in " << t << " ms." << endl;

  C4CTangentialCoverGeometry tcover_geometry;
  Frame2D frame;
  frame.init( ks, 0, 1 );
  C4CTangentialCoverGeometry::MSGeometryComputerByLP geocomp;

  Clock::startClock();
  tcover_geometry.init( tcover, cp, geocomp, frame );
  long t2 = Clock::stopClock();
  cerr << "# Tangential cover geometry in " << t2 << " ms." << endl;

  Clock::startClock();
  TriangleFunction lambda;
  C4CELength* el
    ( tcover_geometry.computeELengthByLambdaMS( tcover,
						tcover_geometry, 
						lambda ) );
  long t3 = Clock::stopClock();
  cerr << "# curvilinear abscissa estimation by l-MST in " << t3 << " ms."
       << endl;

  Clock::startClock();
  global_geom.init( tcover, tcover_geometry, *el, lm );

//   cerr << "LENGTH=" << el->length( 0, 0 ) << endl;
//   for ( uint i = 0; i < el->nbSurfels(); i++ )
//     cerr << "EL[" << i << "]=" << el->elength( i ) << endl;

  global_geom.computeContourGeometry
    ( args.getOption( "-eps" )->getDoubleValue( 0 ),
      args.getOption( "-eps" )->getDoubleValue( 1 ) );

  long t4 = Clock::stopClock();
  cerr << "# GMC optimization in " << t4 << " ms." << endl;
  cerr << "# TIME=" << ( t + t2 + t3 + t4 ) << endl;

  return el;

}


void displayGMCestimation( ostream & out,
			   KnSpace* ks, 
			   C4CIteratorOnSurface & cp,
			   LinearMinimizer* lm,
			   double dh )
{
  bool open;
  uint nb_surfels = C4CIterator::size( cp, open );

  out << setprecision(15);
  out << "#######################################################" << endl;
  out << "# displayGMCestimation" << endl
      << "# displays the estimated geometry of the real shape." << endl 
      << "# Size boundary Nb=" << nb_surfels << " open=" << open << endl
      << "# dh=" << dh << endl
      << "# num est_curv est_tgt_angle cabs dabs el_gmc l_gmc" << endl
      << "# (cabs: curvilinear abscissa by summation of cos/sin of l-MST tangent)" << endl
      << "# (el: surfel elementary length by cos/sin of the GMC tangents )" << endl
      << "# (l: sum of el from beginning)" << endl;


  GlobalC4CGeometry global_geom;
  C4CELength* el = computeGlobalGeometry( global_geom, ks, cp, lm );
  GlobalC4CGeometry::LocalGeometry lgeo;
  Proxy<C4CIteratorOnSurface> it_copy( (C4CIteratorOnSurface*) cp.clone() );

  double dabs = 0.5;
  double l = 0.0;
  double l_lmst = 0.0;
  for ( uint i = 0; i < nb_surfels; ++i )
    {
      Kn_sid bel = it_copy->current();
      Vector v( ks->svectorBasis( bel ) );
      global_geom.geometryFromDiscreteAbscissa( (double) dabs, lgeo );
      double el_int = fabs( v.ro( 0 ) * cos( lgeo.angle_to_x ) 
			+ v.ro( 1 ) * sin( lgeo.angle_to_x ) );
      double el_lmst = el->elength( i );
      out << i 
	  << " " << lgeo.curvature / dh
	  << " " << lgeo.angle_to_x
	  << " " << l_lmst * dh
	  << " " << dabs
	  << " " << el_int * dh
	  << " " << l *dh
	  << endl;
      l += el_int;
      l_lmst += el_lmst;
      dabs += 1.0;
      it_copy->next();
    }
  delete el;
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
  args.addOption( "-minimizer", "-minimizer <type> <a>: choose the minimizer for global geometry computation, <type>=STD mix gradient/optimal, <type>=GD gradient descent, <type>=AGD adaptive step gradient descent, <type>=RLX Relaxation.", "RLX", "0.5" );
  args.addOption( "-eps", "-eps <max> <sum>: specifies max and sum epsilon to stop optimization in global geometry computation.", "0.0000001", "-1.0" );
  
  if ( ( argc <= 0 ) 
       || ! args.readArguments( argc, argv ) ) 
    {
      cerr << args.usage( "curvature_gmc", 
			  "Estimates the curvature along a digital contour given as a Freeman chaincode on the standard input.",
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

  // -------------------------------------------------------------------------
  // Choose Minimizer.
  LinearMinimizer* lm = 0;
  if ( args.getOption( "-minimizer" )->getValue( 0 ) == "GD" )
    lm = new LinearMinimizerByGradientDescent
      ( args.getOption( "-minimizer" )->getDoubleValue( 1 ) );
  else if ( args.getOption( "-minimizer" )->getValue( 0 ) == "AGD" )
    lm = new LinearMinimizerByAdaptiveStepGradientDescent
      ( args.getOption( "-minimizer" )->getDoubleValue( 1 ) );
  else if ( args.getOption( "-minimizer" )->getValue( 0 ) == "RLX" )
    lm = new LinearMinimizerByRelaxation;
  else
    lm = new LinearMinimizer;

  if ( lm != 0 ) 
    {
      cout << "# " << *lm << endl
	   << "# max_eps=" <<  args.getOption( "-eps" )->getDoubleValue( 0 )
	   << " sum_eps=" <<  args.getOption( "-eps" )->getDoubleValue( 1 )
	   << endl;
    }


  // -------------------------------------------------------------------------
  // Estimate curvature
  displayGMCestimation( cout, ks, itfcs, lm, dh ); 

  if ( lm != 0 ) delete lm;
  if ( ks != 0 ) delete ks;
  return 0;
}
