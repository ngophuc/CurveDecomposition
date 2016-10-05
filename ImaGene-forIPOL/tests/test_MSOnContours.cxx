///////////////////////////////////////////////////////////////////////////////
// Tests maximal segments on open and closed contours.
///////////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <vector>
#include "ImaGene/base/Vector2i.h"
#include "ImaGene/base/Arguments.h"
#include "ImaGene/base/StandardArguments.h"
#include "ImaGene/mathutils/Mathutils.h"
#include "ImaGene/mathutils/Statistics.h"
#include "ImaGene/timetools/Clock.h"
#include "ImaGene/dgeometry2d/C4CTangentialCover.h"
#include "ImaGene/dgeometry2d/FreemanChain.h"
#include "ImaGene/dgeometry2d/C4CIteratorOnFreemanChain.h"
#include "ImaGene/digitalnD/C4CIteratorOnFreemanChainSurface.h"
#include "ImaGene/helper/ShapeHelper.h"
#include "ImaGene/helper/DrawingXFIG.h"

using namespace std;
using namespace ImaGene;

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
  StandardArguments::addIOArgs( args, true, false );
  args.addBooleanOption( "-headerXFIG", "-headerXFIG: add xfig header." );
  args.addOption( "-drawContour", "-drawContour <color> <linewidth>: draw the input contour with the specified parameters in XFIG format.", "0", "1" );
  args.addOption( "-drawMS", "-drawMS <color> <linewidth> <ms_mode>: draw the maximal segments of the input contour with the specified parameters in XFIG format.", "0", "1", "0" );
  if ( ( argc <= 0 ) 
       || ! args.readArguments( argc, argv ) ) 
    {
      cerr << args.usage( "test_MSOnContours", 
			  "Tests maximal segments on open and closed contours."
			  , "" ) << endl;
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

  int nb_turn_ccw = c.isClosed();
  KnSpace* ks = ShapeHelper::makeSpaceFromFreemanChain( c );
  C4CIteratorOnFreemanChain itfc;
  itfc.init( c.begin(), nb_turn_ccw != 0 );
  C4CIteratorOnFreemanChainSurface itfcs;
  itfcs.init( ks, itfc );

  Clock::startClock();
  C4CTangentialCover tcover;
  tcover.init( itfcs, 0 );
  long t = Clock::stopClock();
  cerr << "# Tangential cover in " << t << " ms." << endl;

  if ( args.check( "-headerXFIG" ) )
    DrawingXFIG::includeXFIGHeader( cout, 1200, 1 );

  if ( args.check( "-drawContour" ) )
    {
      int color = args.getOption( "-drawContour" )->getIntValue( 0 );
      int linewidth = args.getOption( "-drawContour" )->getIntValue( 1 );
      DrawingXFIG::drawContour( cout, c, color, linewidth, 0, 0, 50);
    }

  if ( args.check( "-drawMS" ) )
    {
      int color = args.getOption( "-drawMS" )->getIntValue( 0 );
      int linewidth = args.getOption( "-drawMS" )->getIntValue( 1 );
      int ms_mode = args.getOption( "-drawMS" )->getIntValue( 2 );
      Frame2D frame;
      frame.init( ks, 0, 1 );
      for ( uint i = 0; i < tcover.nbMaximalSegments(); ++i )
	{
	  const C4CTangentialCover::MaximalSegment & ms 
	    = tcover.getMaximalSegment( i );
	  C4CIteratorOnSurface* itfront 
	    = (C4CIteratorOnSurface*) tcover.getSurfelFront( i );
	  Kn_sid surfel = itfront->current();
	  frame.setSurfelFrame( surfel, 
			        ( ks->sorthDir( surfel ) + 1 ) % 2 );
	  DrawingXFIG::drawDSS( cout, ms.dss, frame, color, 1, 0.0, 0.0, 
				linewidth, ms_mode, 0.0 );
	  delete itfront;
	}
    }
  cerr << "size = " << tcover.nbSurfels() << endl
       << "open = " << ( tcover.isContourOpen() ? "true" : "false" ) << endl
       << "nbms = " << tcover.nbMaximalSegments() << endl;
  
}

