/////////////////////////////////////////////////////////////////////////////
// Display meaningful scale level obtained after multi scales analysis 
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

static uint x0Decal =0;
static uint y0Decal =0;


static Arguments args;
static const int RESOLUTION=1200;
static int samplingSizeMax = 20;



Vector2i getPointFromFreemanChain(const FreemanChain &fc, uint pos);
uint estimMaxSamplingSize(FreemanChain fc);



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
  args.addOption("-setSamplingSizeMax", "-setSamplingSizeMax <max_scale>: set the maximal scale used for contour analysis ", "20" );
  args.addBooleanOption("-processAllContours", "-processAllContours: process all contours (by default process onmly the first point).");
  
  //affichage du bruit
  args.addBooleanOption("-printNoiseLevel", "-printNoiseLevel: displays noise level for each surfel.");
  args.addOption("-setFileNameFigure", "-setFileNameFigure <name> set the <name> for the output of the -drawContour and displayBoxes", "" );
  args.addOption("-setFileNameNoiseLevel", "-setFileNameNoiseLevel <name> set the <name> for the output of the -drawContour and displayBoxes", "" );
  
  args.addBooleanOption( "-enteteXFIG", "-enteteXFIG: add xfig entete" );
  args.addOption( "-drawContourSRC", "-drawContourSRC <color> <linewidth>: draw the original source contour.", "1", "1");
  args.addBooleanOption("-drawXFIGNoiseLevel", "-drawXFIGNoiseLevel" );
  

  args.addBooleanOption("-estimSamplingSizeMax", "-estimSamplingSizeMax: return an estimation of the maximal possible sampling size " );
  
  
  // def de bruit:
  args.addOption( "-meaningfulScale", "-meaningfulScale <min_size> <max_slope>: specifies parameters for defining meaningful scales: minimum size of the interval of scales and maximum slopes between consecutive samples within. default parameter minSize= 1 and minSlopes 0.0", "1", "0.0" );  

  //Affichage de stats:
  args.addOption("-affBoxesStat", "-affBoxesStat <scale> parameter", "1.0");
  
 
  

  


  // pour afficher l'image source.
  args.addOption("-setPosImage", "-setPosImage <x> <y>", "0", "0");
  args.addOption("-afficheImage", "-afficheImage <strinf nomImage> <largeur> <hauteur>","0", "0", "0");
  
  
  



  
  if ( ( argc <= 0 ) 
       || ! args.readArguments( argc, argv ) ) 
    {
      cerr << args.usage( "meaningfulScaleEstim", 
			  "meaningfulScaleEstim  display the first level of scale without noise."
			  ,"" ) << endl;
      return 1;
    }

    
  static int    agrandissementEPS = 30;




    



  


    
  
  
  uint mscales_min_size = args.getOption( "-meaningfulScale" )->getIntValue( 0 );
  double mscales_max_slope = args.getOption( "-meaningfulScale" )->getDoubleValue( 1 );
    

  
  // -------------------------------------------------------------------------
  // Read Freeman chain and creates space/contour.  
  
  std::vector<FreemanChain> vectFC;
  bool containsFC=true;
  while( containsFC){
    FreemanChain fc; 
    istream & in_str = StandardArguments::openInput( args );
    FreemanChain::read( in_str, fc );
    containsFC = in_str.good();

    if ( ! in_str.good() )
      {
	if(vectFC.size()==0){
	  cerr << "Error reading Freeman chain code." << endl;
	  return 2;
	}else{
	  continue;
	}
      }
    containsFC=containsFC && args.check("-processAllContours");
    vectFC.push_back(fc);
  }
  if(args.check("-processAllContours")){
    cerr << "Read " << vectFC.size() << " Freemanchains" << std::endl;
  }
  


  
  ofstream ofFig;
  ofstream ofNoise;

  if(args.check("-setFileNameFigure")){
    string name = args.getOption("-setFileNameFigure")->getValue(0);
    ofFig.open(name.c_str(), ios_base::out);
  }

  if(args.check("-enteteXFIG")){
    args.check("-setFileNameFigure")? ofFig : cout << "#FIG 3.2 \nLandscape \nCenter \nInches \nLetter  \n" <<agrandissementEPS<< "\nSingle \n1" 
						   << " \n"<<RESOLUTION<<" 1" << endl;

    ofFig << "#FIG 3.2 \nLandscape \nCenter \nInches \nLetter  \n" <<agrandissementEPS<< "\nSingle \n1" 
	  << " \n"<<RESOLUTION<<" 1" << endl;
  }
  if(args.check("-afficheImage")){    
    string nomImage = args.getOption( "-afficheImage" )->getValue( 0 );
    int largeur = args.getOption( "-afficheImage" )->getIntValue( 1 );
    int hauteur = args.getOption( "-afficheImage" )->getIntValue( 2 );
    if(args.check("-setPosImage")){
      int posX = args.getOption("-setPosImage")->getIntValue(0);
      int posY = args.getOption("-setPosImage")->getIntValue(1);
      DrawingXFIG::drawImage(args.check("-setFileNameFigure")? ofFig :cout, nomImage, posX, posY, largeur, hauteur, 100);
      
    }else{
      DrawingXFIG::drawImage(args.check("-setFileNameFigure")? ofFig :cout, nomImage, 0, 0, largeur, hauteur, 100);
    }        
  }
  

  long time = Clock::stopClock();

  Clock::startClock();




  for(int k =0; k<vectFC.size(); k++){
    samplingSizeMax=20;
    if(args.check("-processAllContours")){
      cerr << "Processing contour " << k << endl;    
    }
    FreemanChain fc = vectFC.at(k);
  
    if(args.check("-setSamplingSizeMax")){    
      samplingSizeMax = args.getOption("-setSamplingSizeMax")->getIntValue(0);    
    } 
  
    uint samplingSizeMaxEstim = estimMaxSamplingSize(fc);
  
    if(samplingSizeMaxEstim<samplingSizeMax)
      samplingSizeMax= samplingSizeMaxEstim;

  

    if(args.check("-setFileNameNoiseLevel")){
      string name = args.getOption("-setFileNameNoiseLevel")->getValue(0);
      cerr << "name " << name <<endl;
      ofNoise.open(name.c_str(), ios_base::out);
    }

  

  
    if(args.check("-drawContourSRC")){
      uint color = args.getOption("-drawContourSRC")->getIntValue(0);
      uint linewidth = args.getOption("-drawContourSRC")->getIntValue(1);
      DrawingXFIG::setFillIntensity(100);
      DrawingXFIG::drawContour(args.check("-setFileNameFigure")? ofFig :cout, fc, color, linewidth,0,0 , 2);    
    }
  



    //Computing the noise level for each pixel:    
  
    int nbIterationSpikes = 5;
  

  
    FreemanChainSubsample fcsub( 1, 1, 0, 0 );
    FreemanChainCleanSpikesCCW fccs( nbIterationSpikes );
    FreemanChainCompose fcomp( fccs, fcsub );
    FreemanChainTransform* ptr_fct = &fcomp;
    FreemanChainSubsample* ptr_fcsub = &fcsub;
  
    MultiscaleProfile MP;
    MP.chooseSubsampler( *ptr_fct, *ptr_fcsub );
    MP.init( fc, samplingSizeMax );

    Clock::startClock();
    cerr<< "Multi-scale computed in :" << time << " ms" << endl;
    cerr << "Contour size: " << fc.chain.size() << " surfels" << endl;
    cerr << "Sampling size max used: " << samplingSizeMax << endl ;
  
    if(args.check("-printNoiseLevel")){
      ((args.check("-setFileNameNoiseLevel"))? ofNoise :cout) << "# -printNoiseLevel: displays noise level for each surfel." 
							      << endl
							      << "# idx noiselvl code x y" << endl;  
    }

  
    Statistics statBoxes (1, false);
  
    double h= args.getOption("-affBoxesStat")->getFloatValue(0);
  

    uint i =0;
    for (FreemanChain::const_iterator it = fc.begin() ; it!=fc.end(); ++it){
      Vector2i pt ((*it).x(), (*it).y());
    
      uint noiseLevel =0;
      uint idx = it.getPosition();
      uint code = it.getCode();
    
      Vector2i xy( *it );
      noiseLevel =  MP.noiseLevel(i, mscales_min_size, mscales_max_slope);

    
    
      if(args.check("-affBoxesStat")){
	statBoxes.addValue(0, noiseLevel);
      }
    
    
      if(args.check("-drawXFIGNoiseLevel")){
	DrawingXFIG::setFillIntensity(38);
	if(noiseLevel==0){
	  DrawingXFIG::drawPixel (args.check("-setFileNameFigure")? ofFig :cout, pt, 2, 2, samplingSizeMax, 50);
	}else{
	  DrawingXFIG::drawPixel (args.check("-setFileNameFigure")? ofFig :cout, pt, 1, 1, noiseLevel, 50);
	}
      }
      if(args.check("-printNoiseLevel")){
	(args.check("-setFileNameNoiseLevel")? ofNoise :cout) << idx << " " << noiseLevel << " " << code
							      << " " << xy.x() << " " << xy.y() << endl;
      } 
    
      i++;
    }   


    statBoxes.terminate();
  

    if(args.check("-affBoxesStat")){
      cerr << "# Noise estimation stats: Mean:" << statBoxes.mean(0)*h << " variance: "<< statBoxes.variance(0)*h << " Min: " <<  statBoxes.min(0)*h << " Max: " <<  statBoxes.max(0) << endl ;
    
    }
  
    if(args.check("-estimSamplingSizeMax")){
      cerr << "Possible sampling size max: " << estimMaxSamplingSize(fc) << endl;    
    }

  }
  
  






  long time2 = Clock::stopClock();
  cerr << "total time =" << time+time2 <<endl;
}
























uint 
estimMaxSamplingSize(FreemanChain fc){
  int minX = 0;
  int minY = 0;
  int maxX = 0;
  int maxY = 0;
  
  fc.computeBoundingBox(minX, minY, maxX, maxY);
  int largeur = maxX-minX;
  int hauteur = maxY-minY;
  
  return min( largeur/4, hauteur/4);
  
}
