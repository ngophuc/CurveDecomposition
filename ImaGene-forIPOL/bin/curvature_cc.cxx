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
  
  if ( ( argc <= 0 ) 
       || ! args.readArguments( argc, argv ) ) 
    {
      cerr << args.usage( "curvature_cc", 
			  "Estimates the curvature along a digital contour given as a Freeman chaincode on the standard input. On each surfel, the estimator computes a circle that is circumscribed to the surfel center, the longest front and back standard digital segments. [Coeurjolly et al, 2002]",
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
  // Estimate curvature
  Mathutils::AngleComputer ac;
  C4CTangentialCover tcover;
  Clock::startClock();
  buildTangentialCover( tcover, itfcs, 0 );
  long t = Clock::stopClock();
  cerr << "# Tangential cover in " << t << " ms."
       << " nbsurf=" << tcover.nbSurfels() 
       << " nbms=" << tcover.nbMaximalSegments() << endl;

  Clock::startClock();
  C4CTangentialCoverGeometry tcover_geometry;
  Frame2D frame;
  frame.init( ks, 0, 1 );
  C4CTangentialCoverGeometry::MSGeometryComputerByLP geocomp;
  tcover_geometry.init( tcover, itfcs, geocomp, frame );
  long t2 = Clock::stopClock();
  cerr << "# Tangential cover geometry in " << t2 << " ms." << endl;

  Clock::startClock();
  TriangleFunction lambda;
  Proxy<C4CELength> el
    ( tcover_geometry.computeELengthByLambdaMS( tcover,
						tcover_geometry, 
						lambda ) );
  long t3 = Clock::stopClock();
  cerr << "# curvilinear abscissa estimation by l-MST in " << t3 << " ms."
       << endl;


  cout << setprecision(15);
  cout << "#######################################################" << endl;
  cout << "# displayCCestimation" << endl
       << "# displays the estimated geometry of the real shape." << endl 
       << "# Size boundary Nb=" << tcover.nbSurfels() << endl
       << "# dh=" << dh << endl
       << "# num est_curv est_tgt_angle cabs dabs el l" << endl
       << "# (el: elementary length of a surfel estimated from curvature)" << endl
       << "# (l: sum of elementary lengthes from initial surfel)" << endl;
  
  Clock::startClock();

  double dabs = 0.5;
  double cabs = 0.0;
  C4CTangentialCover::SurfelMaximalSegments sms 
    = tcover.beginSMS( 0 );
  for ( uint idx = 0; idx < tcover.nbSurfels(); ++idx )
    {
      // left/back maximal segment.
      const C4CSegment & left_dss 
	= tcover.getMaximalSegment( sms.begin_ms ).dss;
      uint left_front_surfel
	= tcover.getMaximalSegment( sms.begin_ms ).front_surfel_idx;
      const Frame2D & left_frame 
	= tcover_geometry.sgeometry( left_front_surfel ).frame;

      // right/front maximal segment
      const C4CSegment & right_dss 
	= tcover.getMaximalSegment( sms.end_ms ).dss;
      uint right_front_surfel
	= tcover.getMaximalSegment( sms.end_ms ).front_surfel_idx;
      const Frame2D & right_frame 
	= tcover_geometry.sgeometry( right_front_surfel ).frame;

      // local frame
      const Frame2D & surfel_frame 
	= tcover_geometry.sgeometry( idx ).frame;

      // Extracts three points
      Vector2i Ci( left_frame.transformPoint( left_dss.c_n() ) );
      Vector2i Bi( right_frame.transformPoint( right_dss.cp_n() ) );
      Vector2D A( surfel_frame.transformPoint( Vector2D( 0.5, 0.0 ) ) );
      Vector2D C( Ci.x(), Ci.y() );
      Vector2D B( Bi.x(), Bi.y() );

      // The curvature is estimated by computing the circumcircle to ABC
      double curv = 
	(double) EuclideanGeometry::curvatureCircumscribedCircle
	( A.x(), A.y(), B.x(), B.y(), C.x(), C.y() );

      Vector2D AB( B ); AB -= A;
      Vector2D CA( A ); CA -= C;
      Vector2D BC( C ); BC -= B;
      double radius;
      Vector2D tgtA;
      if ( curv < 0.0 )
	{
	  curv = 0.0;
	  radius = 0.0;
	  tgtA = BC;
	}
      else
	{
	  radius = 1 / curv;
	  Vector2D M( B ); M += C;
	  M /= 2.0;
	  Vector2D CM( M ); CM -= C;
	  double dCM = VectorUtils::norm( CM );
	  VectorUtils::normalize( CM );
	  Vector2D n( CM.y(), -CM.x() );
	  double dMP = radius*radius - dCM*dCM;
	  Vector2D P( M );
	  n *= sqrt( dMP );
	  P += n;
	  Vector2D PA( A ); PA -= P;
	  VectorUtils::normalize( PA );
	  tgtA = Vector2D( PA.y(), -PA.x() );
	}

      if ( VectorUtils::det( CA, AB ) < 0 ) curv = -curv;

      double tgt_angle_to_x = ac.cast( atan2( tgtA.y(), tgtA.x() ) );
      double l_surfel = el->elength( idx );

      cout << idx << " " 
	   << curv / dh << " "
	   << tgt_angle_to_x << " "
	   << cabs << " "
	   << dabs << " "
	   << l_surfel << " "
	   << cabs << " "
	   << endl;

      cabs += l_surfel;
      dabs += 0.5;
      tcover.nextSMS( sms );
    }
  long t4 = Clock::stopClock();
  cerr << "# Curvature by circumcircle in " << t4 << " ms." << endl;
  cerr << "# TIME=" << ( t + t2 + t3 + t4 ) << endl;
  if ( ks != 0 ) delete ks;
  return 0;
}

