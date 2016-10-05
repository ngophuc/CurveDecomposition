///////////////////////////////////////////////////////////////////////////////
// Example 5: computes the lambda-MST tangent estimation on each
// surfel of a digital contour.
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include "ImaGene/base/Proxy.h"
#include "ImaGene/dgeometry2d/FreemanChain.h"
#include "ImaGene/dgeometry2d/C4CSegment.h"
#include "ImaGene/dgeometry2d/C4CGeometry.h"
#include "ImaGene/dgeometry2d/C4CSegmentPencil.h"
#include "ImaGene/digitalnD/K2Space.h"
#include "ImaGene/digitalnD/Frame2D.h"
#include "ImaGene/dgeometry2d/C4CIteratorOnFreemanChain.h"
#include "ImaGene/digitalnD/C4CIteratorOnFreemanChainSurface.h"

using namespace std;
using namespace ImaGene;

// Displays the lambda-MST tangent estimation at each surfel of a
// digital contour given as an iterator.
void lambdaMSTEstimator( K2Space & ks, 
			 const C4CIteratorOnSurface & input_it )
{
  // Functions used for combining all the angle of the maximal
  // segments inside the pencil.
  TriangleFunction l;
  DTriangleFunction lp;

  // Memorizes starting point.
  Proxy<C4CIteratorOnSurface> cur_it
    ( (C4CIteratorOnSurface*) input_it.clone() );
  // The following loop will display the estimated tangent direction
  // at each surfel. It needs a Frame to embed the digital contour.
  Frame2D frame;
  frame.init( &ks, 0, 1 );
  // A pencil has experimentally no more than 7 segments.
  const uint m = 20;
  C4CSegment segments[ m ];
  uint idx = 0;
  do
    {
      // Builds the pencil of maximal segments.
      uint j = 0;
      uint k;
      if ( C4CGeometry::maximalSegments( *cur_it, segments, j, k, m ) )
	{
	  // All geometric computations are made in the local frame of
	  // the current boundary element.
	  C4CSegmentPencil pencil( segments, j, k, m, l, lp );
	  float theta = pencil.angleToX( Vector2D( 0.5, 0.0 ) );
	  // Cast angle in the global frame.
	  Kn_sid surfel = cur_it->current();
	  frame.setSurfelFrame( surfel, ks.stanDir( surfel ) );
	  cout << idx << " " << frame.angleToX( theta ) << endl;
	}
      // Go to next element.
      ++idx;
      if ( cur_it->next() == 0) break;
    }
  while ( ! cur_it->equals( input_it ) );
}

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
  // Displays the tangent estimation at each surfel.
  lambdaMSTEstimator( ks, itfcs );
}
