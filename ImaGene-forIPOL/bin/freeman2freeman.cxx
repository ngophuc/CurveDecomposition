///////////////////////////////////////////////////////////////////////////////
// Transforms a Freeman chaincode
///////////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <vector>

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

  args.addBooleanOption( "-ulzero", "-ulzero: upper left point becomes (0,0)." );
  args.addOption( "-ul", "-ul <x0> <y0>: upper left point becomes (x0,y0).", "1", "1" );
  args.addOption( "-inner_contour", "-inner_contour <s>: extracts the 4-connected contour just inside the given boundary. s=CCW contour turns counterclockwise with the inside to its left, s=CW contour turns clockwise with the inside to its right.", "CCW" );
  args.addOption( "-outer_contour", "-outer_contour <s>: same as -inner_contour with CW and CCW reversed.", "CW" );
  args.addOption( "-clean_outer_spikes", "-clean_outer_spikes <s>: removes spikes pointing outside along the given contour. s=CCW contour turns counterclockwise with the inside to its left, s=CW contour turns clockwise with the inside to its right.", "CCW" );
  args.addOption( "-clean_inner_spikes", "-clean_inner_spikes <s>: same as -clean_outer_spikes with CW and CCW reversed.", "CW" );

  args.addOption( "-subsample", "-subsample <h> <v> <x0> <y0>: creates a (h,v)-subsampled freeman chain by the transformation X=(x-x0) div h, Y=(y-y0) div v.", "2", "2", "0", "0" );

  if ( ( argc <= 0 ) 
       || ! args.readArguments( argc, argv ) ) 
    {
      cerr << args.usage( "freeman2freeman", 
			  "Transforms a digital contour (given as a freeman chaincode on the standard input) into another digital contour according to the options.",
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

  int min_x;
  int min_y;
  int max_x;
  int max_y;
  c.computeBoundingBox( min_x, min_y, max_x, max_y );
  
  cerr << "# SDP contour" << endl
       << "# size=" << c.chain.size() << endl
       << "# bbox=" << min_x << " " << min_y
       << " " << max_x << " " << max_y << endl;

  if ( args.check( "-ulzero" ) )
    {
      c.x0 -= min_x;
      c.y0 -= min_y;
    }
  if ( args.check( "-ul" ) )
    {
      c.x0 -= min_x - args.getOption( "-ul" )->getIntValue( 0 );
      c.y0 -= min_y - args.getOption( "-ul" )->getIntValue( 1 );
    }

  if ( args.check( "-subsample" ) )
    {
      FreemanChain c2;
      vector<uint> c2subc;
      vector<uint> subc2c;
      FreemanChain::subsample( c2, c2subc, subc2c, c, 
			       args.getOption( "-subsample" )->getIntValue( 0 ),
			       args.getOption( "-subsample" )->getIntValue( 1 ),
			       args.getOption( "-subsample" )->getIntValue( 2 ),
			       args.getOption( "-subsample" )->getIntValue( 3 )
			       );
      c.chain = c2.chain;
    }

  if ( args.check( "-inner_contour" ) || args.check( "-outer_contour" ))
    {
      bool ccw = args.check( "-inner_contour" ) ?
	args.getOption( "-inner_contour" )->getValue( 0 ) == "CCW" :
	args.getOption( "-outer_contour" )->getValue( 0 ) == "CW";
      FreemanChain c2;
      vector<uint> o2i;
      vector<uint> i2o;
      FreemanChain::innerContour( c2, o2i, i2o, c, ccw );
      c.chain = c2.chain;
      c.x0 = c2.x0;
      c.y0 = c2.y0;
    }

  if ( args.check( "-clean_outer_spikes" ) || args.check( "-clean_inner_spikes" ))
    {
      bool ccw = args.check( "-clean_outer_spikes" ) ?
	args.getOption( "-clean_outer_spikes" )->getValue( 0 ) == "CCW" :
	args.getOption( "-clean_inner_spikes" )->getValue( 0 ) == "CW";
      FreemanChain c2;
      vector<uint> o2i;
      vector<uint> i2o;
      if ( ! FreemanChain::cleanOuterSpikes( c2, o2i, i2o, c, ccw ) )
	cerr << "Contour with no interior !" << endl;
      c.chain = c2.chain;
    }

  FreemanChain::write( cout, c );
  return 0;
}
