///////////////////////////////////////////////////////////////////////////////
// Example 3 bis: Exporting a 3d shapes to a pgm3d file
///////////////////////////////////////////////////////////////////////////////

#include <fstream>
#include <iostream>
#include "ImaGene/digitalnD/KnShapes.h"
#include "ImaGene/digitalnD/KnCharSet.h"
#include "ImaGene/digitalnD/BelAdjacency.h"
#include "ImaGene/digitalnD/C4CIteratorOnBdry.h"
#include "ImaGene/helper/ShapeHelper.h"


using namespace std;
using namespace ImaGene;

int
main( int argc, char** argv ) 
{
// Creates a 3D space of size 32x32x32
  const uint D = 3;
  Kn_size sizes[ D ] = { 32, 32, 32 }; 
  KnSpace ks( D, sizes );
  // Creates a voxel/spel at coordinate 16, 16, 16
  Kn_size xyz[ D ] = { 16*2+1, 16*2+1, 16*2+1 }; 
  Kn_uid center = ks.ukcode( xyz ); // code the cell
  // Creates a ball centered on 'center' and of radius 13.0
  KnCharSet ball = KnShapes::umakeVolumicSphere( ks, center, 13.0 );
  
  // Exports as PGM image on standard output.
  ShapeHelper::exportToPGM3d( cout, &ks, ball ); 
  
  return 0;
}
