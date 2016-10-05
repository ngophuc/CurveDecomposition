///////////////////////////////////////////////////////////////////////////////
// Example 3: creates a PGM image from a Freeman chain code (2D space).
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include "ImaGene/dgeometry2d/FreemanChain.h"
#include "ImaGene/digitalnD/K2Space.h"
#include "ImaGene/digitalnD/KnCharSet.h"
#include "ImaGene/digitalnD/KnRCellSet.h"
#include "ImaGene/digitalnD/KnShapes.h"
#include "ImaGene/helper/ShapeHelper.h"

using namespace std;
using namespace ImaGene;

int
main( int argc, char** argv ) 
{
  // creates 2D space.
  K2Space ks( 128, 128 ); 
  // creates some 4-connected Freeman chain code. 
  FreemanChain c;
  c.x0 = 104;
  c.y0 = 64;
  c.chain  = "1212121221222122212222212222212110110110110110110111101111112112";
  c.chain += "1222222232322323232323232323233232222112121121211212121212121212";
  c.chain += "2122223232323333333333033303330333033332222222223222222322232232";
  c.chain += "2323333330303030030003000300000300000303323323323323323323333233";
  c.chain += "3333033030000000101001010101010101010110100003303033030330303030";
  c.chain += "3030303003000010101011111111112111211121112111100000000010000001";
  c.chain += "0001001001011111";
  // Builds the set of oriented 1-cells, which is the 4-connected
  // interpixel contour.
  KnRCellSet contour = ShapeHelper::makeContourFromFreemanChain( &ks, c, true );
  // Extracts pixels that are interior to this contour.
  KnCharSet image = KnShapes::ucomputeInterior( ks, contour );
  // Complements this set (i.e. exterior pixels).
  image = ~image;
  // Exports as PGM image on standard output.
  ShapeHelper::exportToPGM( cout, &ks, image );
  return 0;
}
