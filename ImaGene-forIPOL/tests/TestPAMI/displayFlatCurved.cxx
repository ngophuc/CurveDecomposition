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

static uint x0Decal =0;
static uint y0Decal =0;


static Arguments args;
static const int RESOLUTION=1200;
static int samplingSizeMax = 20;



Vector2i getPointFromFreemanChain(const FreemanChain &fc, uint pos);

double getSlopeFromMeaningfulScales(const MultiscaleProfile &mp, uint index, uint min_size, uint max_slope);




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
  args.addBooleanOption( "-enteteXFIG", "-enteteXFIG: add xfig entete" );
  args.addOption( "-drawContourSRC", "-drawContourSRC <color> <linewidth>: draw the original source contour.", "1", "10");
  args.addOption("-setSamplingSizeMax", "-setSamplingSizeMax <max_scale>: set the maximal scale used for contour analysis ", "20" );

  
  // def de bruit:
  args.addOption( "-meaningfulScales", "-meaningfulScales <min_size> <max_slope>: specifies parameters for defining meaningful scales: minimum size of the interval of scales and maximum slopes between consecutivesamples within.", "1", "-0.2" );  
 
  
  
  // options d'affichage
    args.addOption("-setPosImage", "-setPosImage <x> <y>", "0", "0");
      
    args.addBooleanOption("-displayNoisePPM", "-displayNoisePPM <image.pgm>  affiche le bruit en couleur dans l'image source."  );
    args.addBooleanOption("-initDisplayNoise", "-initDisplayNoise" );
    args.addOption("-afficheImage", "-afficheImage <strinf nomImage> <largeur> <hauteur>","0", "0", "0");
    
    args.addOption("-displayFlat", "-displayFlat <threshold> <boxWidth> <boxColor>  diplay with color box the flat areas defined from a threshold on the slope of the linear regression on meaningful scale", "-0.6", "1",  "3" );

    
    args.addOption("-saveFlatCurvedInfo", "-saveFlatCurvedInfo <filename> save contour information on flat/curved areas (default: saved in flatCurvedArea.txt) ", 
		   "flatCurvedArea.txt");
    
  if ( ( argc <= 0 ) 
       || ! args.readArguments( argc, argv ) ) 
    {
      cerr << args.usage( "displayFlatCurved", 
			  "displayFlatCurved  display the flat/curved areas."
			  ,"-saveFlatCurvedInfo" ) << endl;
      return 1;
    }
  
    
    static int    agrandissementEPS = 1;

  if(args.check("-enteteXFIG")){
    cout << "#FIG 3.2 \nLandscape \nCenter \nInches \nLetter  \n" <<agrandissementEPS<< "\nSingle \n1" 
	 << " \n"<<RESOLUTION<<" 1" << endl;
  }
  


  
  uint mscales_min_size = args.getOption( "-meaningfulScales" )->getIntValue( 0 );
  double mscales_max_slope = args.getOption( "-meaningfulScales" )->getDoubleValue( 1 );
  bool standard_scale = args.check( "-standardScale" );
  bool standard_scaleMax = args.check( "-standardScaleMax" );
  
  uint n = args.getOption( "-standardScale" )->getIntValue( 0 );    
  double alpha = args.getOption( "-standardScale" )->getDoubleValue( 1 );    

  if(standard_scaleMax){
    n = args.getOption( "-standardScaleMax" )->getIntValue( 0 );    
    alpha = args.getOption( "-standardScaleMax" )->getDoubleValue( 1 );    
  }
  
  

    

  
  // -------------------------------------------------------------------------
  // Read Freeman chain and creates space/contour.  
  FreemanChain fc; 
  istream & in_str = StandardArguments::openInput( args );
  FreemanChain::read( in_str, fc );
  if ( ! in_str.good() )
    {
      cerr << "Error reading Freeman chain code." << endl;
      return 2;
    }
  



  
    

  

  if(args.check("-setSamplingSizeMax")){
     samplingSizeMax = args.getOption("-setSamplingSizeMax")->getIntValue(0);    
  } 








  
  //Computing the noise level for each pixel:    
  
  int nbIterationSpikes = 15;


  Clock::startClock();
  cerr << "starting computing multi resolution and analysis defined on " << fc.chain.size() << endl;
  FreemanChainSubsample fcsub( 1, 1, 0, 0 );
  FreemanChainCleanSpikesCCW fccs( nbIterationSpikes );
  FreemanChainCompose fcomp( fccs, fcsub );
  FreemanChainTransform* ptr_fct = &fcomp;
  FreemanChainSubsample* ptr_fcsub = &fcsub;
  
  MultiscaleProfile MP;
  MP.chooseSubsampler( *ptr_fct, *ptr_fcsub );
  MP.init( fc, samplingSizeMax );
  long time = Clock::stopClock();
  cerr<< "Multi-scale computed in :" << time << endl;
    

  
  ofstream ofNoise;
 
  bool savingFaltCurved = args.check("-saveFlatCurvedInfo");
  string saveName = args.getOption("-saveFlatCurvedInfo")->getValue(0);
  if(savingFaltCurved){
    ofNoise.open(saveName.c_str(), ios_base::out);
    ofNoise<< "# x y noiseLevel  0|1   0=FLAT 1=CURVED " << endl;
  }
  
  
  uint i =0;
  for (FreemanChain::const_iterator it = fc.begin() ; it!=fc.end(); ++it){
    Vector2i pt ((*it).x(), (*it).y());
    
    uint noiseLevel =0;
    
    
    
      noiseLevel =  MP.noiseLevel(i, mscales_min_size, mscales_max_slope);
    
    if(args.check("-displayFlat")){
      
      double threshold = args.getOption("-displayFlat")->getFloatValue(0);
      double slope = getSlopeFromMeaningfulScales(MP, i, mscales_min_size, mscales_max_slope);
      int boxSize =  args.getOption("-displayFlat")->getIntValue(1);
      int boxColor =  args.getOption("-displayFlat")->getIntValue(2);

      
      if(savingFaltCurved){
	ofNoise << pt.x() << " " << pt.y() << " " ;
      }

      bool isCurved=false;
      if(slope<threshold){
	//DrawingXFIG::setFillIntensity(30);
	//	DrawingXFIG::drawPixel (cout, pt, 1, 1, 3, 10);
	if(savingFaltCurved){
	  ofNoise<<0 ;
	}
      }	else if(slope<=0.0){
	DrawingXFIG::setFillIntensity(20);
	DrawingXFIG::drawPixel(cout, pt, boxColor, boxColor, boxSize, 30);
	isCurved=true;
	if(savingFaltCurved){
	  ofNoise<<1 ;
	}
      }else{
	DrawingXFIG::setFillIntensity(20);
	DrawingXFIG::drawPixel(cout, pt, 2, 2, boxSize, 30);
	if(savingFaltCurved){
	  ofNoise<<2 ;
	}
      }
      
      
      
    }
    if(savingFaltCurved){
      ofNoise << " " << noiseLevel << endl;
      
    }

    
    DrawingXFIG::setFillIntensity(38);
    //    DrawingXFIG::drawPixel (cout, pt, 1, 1, noiseLevel, 50);
    i++;
  } 
  
  ofNoise.close();

  if(args.check("-drawContourSRC")){
    uint color = args.getOption("-drawContourSRC")->getIntValue(0);
    uint linewidth = args.getOption("-drawContourSRC")->getIntValue(1);
    DrawingXFIG::setFillIntensity(100);
    DrawingXFIG::drawContour(cout, fc, color, linewidth,0,0 , 2);
    
  }
  
  if(args.check("-afficheImage")){    
    string nomImage = args.getOption( "-afficheImage" )->getValue( 0 );
    int largeur = args.getOption( "-afficheImage" )->getIntValue( 1 );
    int hauteur = args.getOption( "-afficheImage" )->getIntValue( 2 );
    if(args.check("-setPosImage")){
      int posX = args.getOption("-setPosImage")->getIntValue(0);
      int posY = args.getOption("-setPosImage")->getIntValue(1);
      DrawingXFIG::drawImage(cout, nomImage, posX, posY, largeur, hauteur, 100);
      
    }else{
      DrawingXFIG::drawImage(cout, nomImage, 0, 0, largeur, hauteur, 100);
    }    
    
  }
    

  
    

}










Vector2i 
getPointFromFreemanChain(const FreemanChain &fc, uint pos){
  int tailleFc = fc.chain.size();
  FreemanChain::const_iterator it = fc.begin();  
  // recherche du point à l'indice pos
  int i=0;
  while(i!=pos && i< tailleFc){
    it.nextInLoop();
    i++;
  }
  return *it;
}










double
getSlopeFromMeaningfulScales(const MultiscaleProfile &mp, uint index, uint min_size, uint max_slope){
  vector< std::pair< uint, uint > > vectIntervals;
  mp.meaningfulScales(vectIntervals,  index, min_size, max_slope);
  vector<double> vx;
  vector<double> vy;
  mp.profile(vx, vy, index);
  
  
  uint mfsDeb=0;
  uint mfsFin=mp.all_stats.size();
  if(vectIntervals.size()!=0){
    mfsDeb=vectIntervals.at(0).first-1;
    mfsFin=vectIntervals.at(0).second-1;        
  }else{
    return 100.0;
  }

  SimpleLinearRegression SLR;
  for(int i=mfsDeb; i<=mfsFin; i++){
    SLR.addSample(vx.at(i), vy.at(i));
  }
  if ( ! SLR.computeRegression() ){
    Vector2D pt1(vx.at(mfsDeb),vy.at(mfsDeb));
    Vector2D pt2(vx.at(mfsFin),vy.at(mfsFin));
    double slope = (pt2.y() - pt1.y())/(pt2.x() -pt1.x());
    cerr << "PROBLEM in SLR" << endl;
    cerr << "size: " << SLR.m_n << "slope " << SLR.slope() << "mfs:"<<  mfsDeb << " " << mfsFin << "slope R:" << slope<<  endl;
    return slope;
  }

  return SLR.slope();
}


























