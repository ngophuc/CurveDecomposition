///////////////////////////////////////////////////////////////////////////////
// Generates a pgm image from a freeman contour.
///////////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <deque>

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
  args.addBooleanOption( "-auto_center", "-auto_center: automatically centers the shape in the PGM image." );
  
  if ( ( argc <= 0 ) 
       || ! args.readArguments( argc, argv ) ) 
    {
      cerr << args.usage( "test_Reconstruction", 
			  "Reconstructs a polygon/spline from a Freeman chaincode.",
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

  uint d = StandardArguments::dim( args );
  if ( d != 2 )
    {
      cerr << "Dimension is 2." << endl;
      return 2;
    }
  Kn_size sizes[ d ];
  StandardArguments::fillSizes( args, sizes );
  KnSpace* ks;
  if ( args.check( "-auto_center" ) )
    ks = ShapeHelper::makeSpaceFromFreemanChain( c );
  else 
    ks = new KnSpace( 2, sizes );

  // KnRCellSet contour = ShapeHelper::makeContourFromFreemanChain( ks, c, true );

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
  //  P.extractEdges();
}
