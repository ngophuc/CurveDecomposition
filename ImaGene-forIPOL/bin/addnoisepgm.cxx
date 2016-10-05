///////////////////////////////////////////////////////////////////////////////
// Add noise to a pgm image
///////////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>

#include "ImaGene/base/Arguments.h"
#include "ImaGene/base/Vector.h"
#include "ImaGene/base/Shapes.h"
#include "ImaGene/digitalnD/GridEmbedder.h"
#include "ImaGene/digitalnD/KnSpace.h"
#include "ImaGene/digitalnD/KnCharSet.h"
#include "ImaGene/digitalnD/KnShapes.h"
#include "ImaGene/digitalnD/C4CIteratorOnSurface.h"
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
  ShapeHelper::addNoiseArgs( args );
  // -------------------------------------------------------------------------
  // Prepare arguments.
  args.addOption( "-threshold", "-threshold <val>: threshold value for binarizing PGM gray values (def. is 128).", "128" );
  
  if ( ( argc <= 0 ) 
       || ! args.readArguments( argc, argv ) ) 
    {
      cerr << args.usage( "addnoisepgm", 
			  "Add (binary) noise to the binarization of input image and outputs the image on the standard output.",
			  "" )
	   << endl;
      return 1;
    }

  // reading input image.
  KnSpace* ks;
  KnCharSet* voxset;
  uint threshold = (uint) args.getOption( "-threshold" )->getIntValue( 0 );
  if ( ! ShapeHelper::importFromPGM( cin, ks, voxset, threshold ) )
    {
      cerr << "Error reading PGM file." << endl;
      return 2;
    }
  KnCharSet voxset2( *voxset );

  // -------------------------------------------------------------------------
  // Change shape.
  // KnCharSet aux_voxset( *voxset );
  Kn_uid in, out;
  if ( args.check( "-noisify" ) )
    {
      Kn_uid bel = 0;
      ShapeHelper::addNoiseFromArgs( args, ks, voxset2, bel, in, out );

//       if ( args.check( "-inner_object" ) )
// 	{
// 	  // NB: change contour and object.
// 	  bel = ShapeHelper::findInnerObject( ks, *voxset, in, aux_voxset );
// 	  delete voxset;
// 	  voxset = &aux_voxset;
// 	}
//       else if ( args.check( "-outer_object" ) )
// 	{
// 	  // NB: change contour and object.
// 	  bel = ShapeHelper::findOuterObject( ks, *voxset, out, 
// 					      aux_voxset, in );
// 	  delete voxset;
// 	  voxset = &aux_voxset;
// 	}
//       else delete voxset;
    }
  //  else delete voxset;
  
  //ShapeHelper::exportToPGM( cout, ks, aux_voxset );
  ShapeHelper::exportToPGM( cout, ks, voxset2 );
  delete voxset;
  delete ks;

  return 0;
}
