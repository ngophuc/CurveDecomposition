///////////////////////////////////////////////////////////////////////////////
// Test the length variation of maximal segments on digital contour
///////////////////////////////////////////////////////////////////////////////
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
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



using namespace std;
using namespace ImaGene;

static uint x0Decal =0;
static uint y0Decal =0;


static Arguments args;

Vector2D  getRegLinear(const  vector<Vector2D> &vectorPoints );


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
  args.addOption( "-samplingSizeMax", "-samplingSizeMax <n>: choose how many scales are computed.", "10"  );
  args.addOption("-idx", "-idx <index> set the index of the point for which the profile is extracter", "0"); 
  args.addOption("-affProfileRegs", "-affProfileRegs <filename> extract slopes of the linear regression of all contour points and save it in <filename> "," ");

  

  
  if ( ( argc <= 0 ) 
       || ! args.readArguments( argc, argv ) ) 
    {
      cerr << args.usage( "test_Profiles", 
			  "test_Profiles -idx <index>  shows the  profile for the point of index idx of contour given in std input."
			  ,"" ) << endl;
      return 1;
    }
  



  
  // -------------------------------------------------------------------------
  // Read Freeman chain and creates space/contour.  
  
  FreemanChain c; 
  istream & in_str = StandardArguments::openInput( args );
  FreemanChain::read( in_str, c );
  if ( ! in_str.good() )
    {
      cerr << "Error reading Freeman chain code." << endl;
      return 2;
    }


  // -------------------------------------------------------------------------
  // Read some arguments.
  uint samplingSizeMax = args.getOption( "-samplingSizeMax" )->getIntValue( 0 );
  uint index = args.getOption("-idx")->getIntValue(0);
  
  int nbIterationSpikes = 10;
  
  FreemanChainSubsample fcsub( 1, 1, 0, 0 );
  FreemanChainCleanSpikesCCW fccs( nbIterationSpikes );
  FreemanChainCompose fcomp( fccs, fcsub );
  FreemanChainTransform* ptr_fct = &fcomp;
  FreemanChainSubsample* ptr_fcsub = &fcsub;
  
  MultiscaleProfile MP;
  MP.chooseSubsampler( *ptr_fct, *ptr_fcsub );
  MP.init( c, samplingSizeMax );
  



  if(args.check("-affProfileRegs")){
    fstream strProfilesReg;
    const string name = args.getOption("-affProfileRegs")->getValue(0);
    strProfilesReg.open(name.c_str(), ios::out);    
    
    for(int i =0; i<c.chain.size();i++){
      vector<double> vectx; 
      vector<double> vecty; 
      MP.profile(vectx, vecty, i);
      vector<Vector2D> vPoints;
      for( int k=0; k<vectx.size(); k++){
	Vector2D v(vectx.at(k), vecty.at(k)); 
	vPoints.push_back(v);
      }
      Vector2D reg = getRegLinear(vPoints );
      strProfilesReg<< i << " " << reg.x() << endl; 
    }
    
    strProfilesReg.close();
  }



  vector<double> vx; 
  vector<double> vy; 
  MP.profile(vx, vy, index);
    
  vector<double> vxMedian; 
  vector<double> vyMedian; 
  MP.profileFromMedian(vxMedian, vyMedian, index);
  


  fstream strProfiles;
  strProfiles.open("profiles.txt", ios::out);

  
  for(int i=0; i< vx.size(); i++){
    strProfiles << vx.at(i) << " " << vy.at(i) << " " << vyMedian.at(i) << endl;    
  }
  
  DrawingXFIG::includeXFIGHeader (cout, 1024,2);
  DrawingXFIG::drawContour (cout, c, 0, 10, 0, 0, 10);

  FreemanChain::const_iterator it = c.begin();
  for (int i=0; i< index; i++){
    ++it;
  }  
  Vector2i pt ((*it).x(), (*it).y());
  DrawingXFIG::drawCross (cout,pt ,1, 3, 10, 40);
  

  
  exit(0);

}


















// Return the two coefficients of the line obtained from linear regression
Vector2D 
getRegLinear(const  vector<Vector2D> &vectorPoints ){
  int size = vectorPoints.size();
  Statistics statReg(2, true);
  double coVariance  =0.0;
  
  for(int k=0; k< size;  k++){           
    Vector2D point = vectorPoints.at(k);
    statReg.addValue(0,point.x());      	
    statReg.addValue(1,point.y());      	
    coVariance+=(point.x()*point.y());
  }
  statReg.terminate();
  int nbSamples = statReg.samples(0); 
  double sampleMoy = statReg.mean(0);
  double tailleSegMaxMoy = statReg.mean(1);
  coVariance = coVariance/nbSamples;
  coVariance = coVariance - sampleMoy*tailleSegMaxMoy;
  double slope = coVariance/statReg.unbiasedVariance(0);
  double b = statReg.mean(1)-slope*statReg.mean(0);
  return Vector2D(slope, b);
} 

