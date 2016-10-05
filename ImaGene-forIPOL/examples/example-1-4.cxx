///////////////////////////////////////////////////////////////////////////////
// Example 4: extracts maximal segments of a digital contour and
// displays their slopes.
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include "ImaGene/dgeometry2d/FreemanChain.h"
#include "ImaGene/dgeometry2d/C4CTangentialCover.h"
#include "ImaGene/digitalnD/K2Space.h"
#include "ImaGene/digitalnD/Frame2D.h"
#include "ImaGene/dgeometry2d/C4CIteratorOnFreemanChain.h"
#include "ImaGene/digitalnD/C4CIteratorOnFreemanChainSurface.h"

using namespace std;
using namespace ImaGene;

int
main( int argc, char** argv ) 
{
  // creates 2D space.
  K2Space ks( 128, 128 ); 
  // creates some 4-connected Freeman chain code, ie a digital contour. 
  FreemanChain c;
  c.x0 = 36;
  c.y0 = 32;
  c.chain  = "1211221232233230330030100110";
  // Creates an iterator on the Freeman chain that tells if the path
  // is going straight, or turning left or right.
  C4CIteratorOnFreemanChain itfc;
  itfc.init( c.begin(), true );
  // This customized iterator knows also what is the cell code for the
  // current step of the Freeman chain.
  C4CIteratorOnFreemanChainSurface itfcs;
  itfcs.init( &ks, itfc );
  // This computes in linear time all the maximal segments of the
  // digital contour.
  C4CTangentialCover tcover;
  tcover.init( itfcs, 0 );

  // The following loop will display slopes of each maximal segment as
  // a 2D vector.
  // needs a Frame to embed the digital contour.
  Frame2D frame;
  frame.init( &ks, 0, 1 );
  uint i = 0;
  for ( uint m = 0; m < tcover.nbMaximalSegments(); ++m )
    {
      // Obtain the current maximal segment.
      const C4CTangentialCover::MaximalSegment & ms 
	= tcover.getMaximalSegment( m );
      // Place iterator at the front of the maximal segment.
      while ( i != ms.front_surfel_idx ) 
	{
	  i = (i+1) % tcover.nbSurfels();
	  itfcs.next();
	}
      // Gets the cell at front of iterator. It defines the local
      // frame of the digital straight segment.
      Kn_sid front_surfel = itfcs.current();
      frame.setSurfelFrame( front_surfel, ks.stanDir( front_surfel ) );
      // Transform the tangent vector into the global reference frame.
      const C4CSegment & segment = ms.dss; 
      Vector2i tgt = frame.transform( segment.getTangent() );
      cout << "[" << m <<  "] tgt=(" << tgt.x() << "," << tgt.y() << ")"
	   << endl;
    }
  return 0;
}
