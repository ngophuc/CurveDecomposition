
#include "ImaGene/base/Arguments.h"
#include "ImaGene/base/StandardArguments.h"
#include "ImaGene/dgeometry2d/C4CTangentialCover.h"
#include "ImaGene/digitalnD/KnSpace.h"
#include "ImaGene/digitalnD/KnCharSet.h"
#include "ImaGene/digitalnD/KnShapes.h"
#include "ImaGene/digitalnD/C4CIteratorOnSurface.h"
#include "ImaGene/digitalnD/BelAdjacency.h"
#include "ImaGene/digitalnD/ObjectBoundary.h"
#include "ImaGene/digitalnD/KnRCellSet.h"
#include "ImaGene/digitalnD/KnRCellVector.h"
#include "ImaGene/helper/C4CTangentialCoverGeometry.h"
#include "ImaGene/helper/CharacteristicPolygon.h"
#include "ImaGene/helper/ShapeHelper.h"

#include <iostream>
#include <fstream>



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
     
  StandardArguments::addDigitalArgs( args, 3, false, false );
  ShapeHelper::addSimple3DShapesArgs( args );
 
  
  if ( ( argc <= 1 ) 
       || ! args.readArguments( argc, argv ) ) 
    {
      cerr << args.usage( "gen3dShapes ", 
			  "Display shapes as sets of voxels",
			  "" )
	   << endl;
      return 1;
    }


  // -------------------------------------------------------------------------
  // Build space.
  uint d = StandardArguments::dim( args );
  if ( d != 3 )
    {
      cerr << "Dimension should be 3." << endl;
      return 2;
    }
  Kn_size sizes[ d ];
  StandardArguments::fillSizes( args, sizes );
  KnSpace ks( d, sizes );
  cerr << "--- Space: " << ks << endl;

  KnCharSet voxset = ShapeHelper::makeSimple3DShapesFromArgs( args, ks  );

  cerr << "voxset size=" << voxset.nbElements()<< endl;
  ShapeHelper::exportToPGM3d (cout, &ks, voxset);
    
  
  
}
  
  
