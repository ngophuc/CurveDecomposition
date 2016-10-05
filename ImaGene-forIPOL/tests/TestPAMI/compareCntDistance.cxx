/////////////////////////////////////////////////////////////////////////////
// Display the Noise level obtained after multi scales analysis 
//////////////////////////////////////////////////////////////////////////////
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <algorithm>
#include "ImaGene/base/Proxy.h"
#include "ImaGene/dgeometry2d/C4CGeometry.h"
#include "ImaGene/dgeometry2d/C4CSegmentPencil.h"
#include "ImaGene/base/Vector2i.h"
#include "ImaGene/base/Arguments.h"
#include "ImaGene/base/StandardArguments.h"
#include "ImaGene/mathutils/Mathutils.h"
#include "ImaGene/mathutils/Statistics.h"
#include "ImaGene/timetools/Clock.h"
#include "ImaGene/dgeometry2d/C4CTangentialCover.h"
#include "ImaGene/helper/C4CTangentialCoverGeometry.h"
#include "ImaGene/dgeometry2d/FreemanChain.h"
#include "ImaGene/dgeometry2d/FreemanChainTransform.h"
#include "ImaGene/dgeometry2d/C4CIteratorOnFreemanChain.h"
#include "ImaGene/digitalnD/C4CIteratorOnFreemanChainSurface.h"
#include "ImaGene/digitalnD/KnSpace.h"
#include "ImaGene/digitalnD/C4CIteratorOnSurface.h"
#include "ImaGene/helper/ShapeHelper.h"
#include "ImaGene/helper/MultiscaleProfile.h"
#include "ImaGene/helper/DrawingXFIG.h"
#include "ImaGene/helper/CurveVariationsHelper.h"
#include "ImaGene/mathutils/SimpleLinearRegression.h"



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
  StandardArguments::addIOArgs( args, true, false );
 
  
 
 
  args.addOption("-cntA", "-cntA <contourA.fc>" , " ");
  args.addOption("-cntB", "-cntB <contourB.fc>" , " ");
  args.addOption("-gridSize", "-gridSize <g> " , "1.0");
  
  
  
  
  
  if ( ( argc <= 1 )|| ! args.readArguments( argc, argv ))
    {
      cerr << args.usage( "compareCntDistance", 
			  "compareCntDistance -cntA <contourA.fc> -cntB <contourB.fc> : return the sum of the totgal distance of the contour B to the contour B."
			  ,"" ) << endl;
      return 1;
    }
  
  
  FreemanChain fcA;
  FreemanChain fcB;
  
  string cntA = args.getOption("-cntA")->getValue(0);
  string cntB = args.getOption("-cntB")->getValue(0);

  
  
  ifstream in;
  in.open(cntA.c_str(), ifstream::in);
  FreemanChain::read( in, fcA );
  in.close();
  in.open(cntB.c_str(), ifstream::in);
  FreemanChain::read( in, fcB );
  
  
  KnSpace* ks = ShapeHelper::makeSpaceFromFreemanChain(fcB);
  KnSpace* ksA = ShapeHelper::makeSpaceFromFreemanChain(fcA);
  KnRCellSet contourA = ShapeHelper::makeContourFromFreemanChain(ks, fcA, true );
  KnRUCellVector<int>* dmap  = KnShapes::computeBdryCityBlockDistanceMap( *ksA, contourA );
  
  KnRCellSet contourB = ShapeHelper::makeContourFromFreemanChain(ks, fcB, true );

  int nb=0;
  int sommeE=0;
  for ( KnRCellSet::cell_iterator p = contourB.begin();
	p != contourB.end();
	++p )
    {
      
      Kn_sid s = *p;
      uint i = ks->sorthDir( s );
      bool d = ks->sdirect( s, i );
      Kn_uid inner_spel = ks->unsigns( ks->sincident( s, i, d ) );
      sommeE+= (abs((*dmap)[ inner_spel ])-1);
      
  
    }


  double gridSize  = args.getOption("-gridSize")->getFloatValue(0);
  cerr << "somme E " << sommeE << " size " << fcB.chain.size()<< " "  << endl;
  cerr << "erreur totale = " << ((gridSize*sommeE)/(double)(fcB.chain.size())) << endl;
  

  
  
//   ImagePGM image1 =  importPGM(str1.c_str());
//   ImagePGM image2 =  importPGM(str2.c_str());
  
  
//   uint numDifference = 0;
  
//   uint width = image1.width;
//   uint height = image1.height;
  

  

  
//   for(int i=0; i<height; i++){
//     for(int j=0; j<width; j++){
//       uint indice = i*width+j;
//       if(image1.tabImage[indice]!=image2.tabImage[indice]){
// 	numDifference++;	
//       }
       
//     }
//   }
  
  

//   cout << "Total number of differences between " << str1 << "  and " << str2 << ": " << numDifference << endl;
  
  
}







