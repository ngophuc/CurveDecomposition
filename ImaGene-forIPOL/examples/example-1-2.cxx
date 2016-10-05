///////////////////////////////////////////////////////////////////////////////
// Example 2: creates a sequence of nD digital spaces, a digital
// sphere within, then extracts its boundary
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include "ImaGene/digitalnD/KnSpace.h"
#include "ImaGene/digitalnD/KnCharSet.h"
#include "ImaGene/digitalnD/KnRCellSet.h"
#include "ImaGene/digitalnD/KnShapes.h"

using namespace std;
using namespace ImaGene;

int main( int argc, char** argv )
{
  // Creates spaces of size 8^n, with increasing n.
  Kn_size sizes[] = { 8, 8, 8, 8, 8, 8, 8, 8 }; 
  for ( uint D = 2; D <= 8; ++ D )
    {
      cout << " - Dimension is " << D << "." << endl;
      KnSpace ks( D, sizes );
      if ( ! ks.OK() ) break;
      // Creates a voxel/spel at coordinate 3, ..., 3
      Kn_size xyz[ D ];
      for ( uint i = 0; i < D; ++i ) xyz[ i ] = 3*2+1;
      Kn_uid center = ks.ukcode( xyz ); // code the cell
      // Creates a ball centered on 'center' and of radius 1.0
      KnCharSet ball = KnShapes::umakeVolumicSphere( ks, center, 1.0 );
      cout << " --- ball has " << ball.nbElements() << " voxels." << endl;
      // Extracts one bel on the ball boundary.
      Kn_sid bel = KnShapes::sfindClosestBel( ks, center, ball );
      // Extracts the ball boundary, ie a sphere.
      KnRCellSet sphere = KnShapes::strackBoundary( ks, ball, bel );
      cout << " --- sphere has " << sphere.nbElements() << " surfels." << endl;
    }
}
