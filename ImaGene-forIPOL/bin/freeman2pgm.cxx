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
  StandardArguments::addDigitalArgs( args, 2, false, false );
  args.addBooleanOption( "-auto_center", "-auto_center: automatically centers the shape in the PGM image." );
  args.addBooleanOption( "-inverse", "-inverse: inverse the orientation of the contour so that the image is the negative of the expected image." );
  
  if ( ( argc <= 0 ) 
       || ! args.readArguments( argc, argv ) ) 
    {
      cerr << args.usage( "freeman2pgm", 
			  "Fills a contour given as a Freeman chain code and output a PGM image. By default, the user should provide the size of the digital space. Otherwise, with -auto_center, the program computes a big enough image for the contour.",
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
  KnRCellSet contour = ShapeHelper::makeContourFromFreemanChain( ks, c, true );

  // -------------------------------------------------------------------------
  // Compute PGM image

  KnCharSet image = KnShapes::ucomputeInterior( *ks, contour );

  bool inverse = args.check( "-inverse" );


  if ( inverse )
    ShapeHelper::exportToPGM( cout, ks, ~image );
  else
    ShapeHelper::exportToPGM( cout, ks, image );

  delete ks;
  return 0;
}
