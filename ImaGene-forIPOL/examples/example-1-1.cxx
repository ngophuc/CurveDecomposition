///////////////////////////////////////////////////////////////////////////////
// Example 1: creates a 3D digital space, a digital sphere and extracts
// its boundary
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
  // Creates a 3D space of size 32x32x32
  const uint D = 3;
  Kn_size sizes[ D ] = { 32, 32, 32 }; 
  KnSpace ks( D, sizes );
  // Creates a voxel/spel at coordinate 12, 13, 14
  Kn_size xyz[ D ] = { 12*2+1, 13*2+1, 14*2+1 }; 
//   Kn_size xyz[ D ] = { 13*2+1, 3*2+1, 3*2+1 }; 
  Kn_uid center = ks.ukcode( xyz ); // code the cell
  cout<<"code du voxel (12,13,14) : "<<center<<endl;
  // Creates a ball centered on 'center' and of radius 3.0
    KnCharSet ball = KnShapes::umakeVolumicSphere( ks, center, 4.0 );
//   KnCharSet ball = KnShapes::umakeVolumicSphere( ks, center, 1.0 );
  cout << " --- ball has " << ball.nbElements() << " voxels." << endl;
  // Extracts one bel on the ball boundary.
  Kn_sid bel = KnShapes::sfindClosestBel( ks, center, ball );
  // Extracts the ball boundary, ie a sphere.
  KnRCellSet sphere = KnShapes::strackBoundary( ks, ball, bel );
  cout << " --- sphere has " << sphere.nbElements() << " surfels." << endl;
}
