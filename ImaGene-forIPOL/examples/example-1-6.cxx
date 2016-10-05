///////////////////////////////////////////////////////////////////////////////
// Example 6: Extract a contour from a grayscale PGM image.
///////////////////////////////////////////////////////////////////////////////

#include <fstream>
#include <iostream>
#include "ImaGene/digitalnD/Frame2D.h"
#include "ImaGene/digitalnD/KnCharSet.h"
#include "ImaGene/digitalnD/BelAdjacency.h"
#include "ImaGene/digitalnD/C4CIteratorOnBdry.h"
#include "ImaGene/helper/ShapeHelper.h"
#include "ImaGene/helper/DrawingXFIG.h"

using namespace std;
using namespace ImaGene;

int
main( int argc, char** argv ) 
{
  KnSpace* ks;
  KnCharSet* voxset;
  
  // Import the set of pixels extracted from a threshold set to 128.
  ShapeHelper::importFromPGM( cin, ks, voxset, 128 );

  
  // Get the sets of signed surfel cells of the voxset bondary.  
  KnRCellSet bdry = KnShapes::smakeBoundary( *ks, *voxset );
  
  // Get the first element of the sets. 
  Kn_sid bel = *(bdry.begin());
  uint k = *( ks->sbegin_dirs( bel ) );
  
  // Display boundary with C4Citerator:
  BelAdjacency badj( *ks, true );  
  C4CIteratorOnBdry c4c_it( badj, bel, k, *voxset );
  Kn_sid belInit = c4c_it.current();
  bel = belInit;
  vector<Vector2D> vectContour;
  
  // Display using XFIG tools
  DrawingXFIG::includeXFIGHeader( cout, 1024, 12 );  

  do {
    if ( c4c_it.next() == 0 ) break;
    bel = c4c_it.current();
    Frame2D frame;
    frame.init( ks, 0, 1 );
    frame.setSurfelFrame( bel, c4c_it.trackDir() );
    Vector2D surfel( frame.transformPoint( Vector2D( 0.0, 0.0 ) ) );    
    Vector2D surfelCenter( frame.transformPoint( Vector2D( 0.5, 0.0 ) ) );        
    vectContour.push_back(surfel);
    DrawingXFIG::drawCircle( cout, surfelCenter, 4, 0.1, 0 );
  } while(bel!=belInit);
  
   
  DrawingXFIG::drawContour( cout, vectContour, 1, 0, 20,true ,false, 10, 1.0, 0.0, 0.0 );
  DrawingXFIG::drawImage( cout, "contourSrc.jpeg", 32, 32, 40 );
	
  return 0;
}
