#include <iostream>
#include "ImaGene/dgeometry2d/FreemanChain.h"
#include "ImaGene/dgeometry2d/FreemanChainTransform.h"
#include "ImaGene/helper/MultiscaleProfile.h"

using namespace ImaGene;
using namespace std;

int
main( int argc, char** argv ) 
{
  // Setting the default parameter:
  uint maxScale = 10;
  double maxSlope = 0.0;
  uint minSize = 1;
  uint nbIterationSpikeDetection = 5;
  // -------------------------------------------------------------------------
  // Read Freeman chain and creates space/contour.  
  FreemanChain fc; 
  FreemanChain::read(cin, fc );
  cerr << "Reading Freeman chain [done]" << endl;
  
  // -------------------------------------------------------------------------
  // Computing the multiscale profile
  FreemanChainSubsample fcsub( 1, 1, 0, 0 );
  FreemanChainCleanSpikesCCW fccs( nbIterationSpikeDetection );
  FreemanChainCompose fcomp( fccs, fcsub );
  FreemanChainTransform* ptr_fct = &fcomp;
  FreemanChainSubsample* ptr_fcsub = &fcsub;
  
  MultiscaleProfile MP;
  MP.chooseSubsampler( *ptr_fct, *ptr_fcsub );
  MP.init( fc, maxScale);
  
  //-------------------------------------------------------------------------
  // Displaying the noise by iterating on the freeman chain.
  int i=0;
  for (FreemanChain::const_iterator it = fc.begin() ; it!=fc.end(); ++it){
    uint noiseLevel =  MP.noiseLevel(i, minSize, maxSlope);
    cout << "noise level point index" << i << ": " << noiseLevel << endl;
    i++;
  }

  return 0;
}

