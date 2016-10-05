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



struct ImagePGM{
  uint width;
  uint height;
  uint maxValue;
  uint *tabImage;
};


struct ImagePPM{
  uint width;
  uint height;
  uint maxValue;
  uint * tabImageR;
  uint * tabImageG;
  uint * tabImageB;
};



void exportPGM(string  name, ImagePGM image);
ImagePGM importPGM(string name);


void exportPPM(string name, ImagePPM image);
ImagePPM importPPM( string name);



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
  args.addOption( "-drawContourSubsampled", "-drawContourSubsampled <color> <linewidth>: draw the original source contour.", "1", "10");
  args.addOption( "-drawContourSubsampledK", "-drawContourSubsampledK <color> <linewidth> <Level>: draw the original source contour subsampled at the resolution K.", "1", "10", "1" );
  args.addOption( "-drawContourSRCOnly", "-drawContourSRCOnly <color> <linewidth>: draw the original source contour.", "1", "10");
  args.addOption("-setSamplingSizeMax", "-setSamplingSizeMax <max_scale>: set the maximal scale used for contour analysis ", "20" );

  
  // def de bruit:
  args.addOption( "-meaningfulScales", "-meaningfulScales <min_size> <max_slope>: specifies parameters for defining meaningful scales: minimum size of the interval of scales and maximum slopes between consecutivesamples within.", "1", "-0.2" );  
  
  

  
  args.addOption("-maxValidScale", "-maxValidScale <max Slope> show the noise defined as largest possible interval containing slope less than <maxSlope>.", "0.0"); 
  args.addOption("-maxMeaningfulScales" ,"-maxMeaningfulScales <minSize> <maxSlope> : return the meaningfulScale with longest maximal segments ","1", "-0.2" );


  // options d'affichage
  args.addOption("-drawProfil", "-drawProfil <index> : draw the multiscale profil of the point <index> and save it in profil.txt", "0");
  args.addOption("-drawDetailedProfil", "-drawDetailedProfil <index> : draw all the points of the profiles.", "0" );
  args.addOption("-plotReg", "-plotReg <index> <minSize>  <alpha> : plot linear regression in regProfil.txt ", "0", "3", "0.1");

    args.addOption("-setPosImage", "-setPosImage <x> <y>", "0", "0");
  //args.addBooleanOption( "-singlePixelMode", "-singlePixelMode draw the noise level of a pixel with its representant in the subsampled contour" );
  
    
    args.addOption("-getSize", "-getSize <name> ", "");
    
    
    args.addBooleanOption("-displayNoisePPM", "-displayNoisePPM <image.pgm>  affiche le bruit en couleur dans l'image source."  );
    args.addBooleanOption("-initDisplayNoise", "-initDisplayNoise" );
    args.addOption("-afficheImage", "-afficheImage <strinf nomImage> <largeur> <hauteur>","0", "0", "0");
    
  

  
  if ( ( argc <= 0 ) 
       || ! args.readArguments( argc, argv ) ) 
    {
      cerr << args.usage( "displayNoise", 
			  "displayNoise  display the first level of scale without noise."
			  ,"" ) << endl;
      return 1;
    }
  
    
    static int    agrandissementEPS = 30;

  if(args.check("-enteteXFIG")){
    cout << "#FIG 3.2 \nLandscape \nCenter \nInches \nLetter  \n" <<agrandissementEPS<< "\nSingle \n1" 
	 << " \n"<<RESOLUTION<<" 1" << endl;
  }
  


  
  uint mscales_min_size = args.getOption( "-meaningfulScales" )->getIntValue( 0 );
  double mscales_max_slope = args.getOption( "-meaningfulScales" )->getDoubleValue( 1 );
  

  

  
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
  



  
    
  if(args.check("-drawContourSRC")){
    uint color = args.getOption("-drawContourSRC")->getIntValue(0);
    uint linewidth = args.getOption("-drawContourSRC")->getIntValue(1);
    DrawingXFIG::setFillIntensity(100);
    DrawingXFIG::drawContour(cout, fc, color, linewidth,0,0 , 2);
    
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
    


  if(args.check("-drawProfil")){
    uint index = args.getOption("-drawProfil")->getIntValue(0);
    vector <double> vx ;
    vector <double> vy ;
    MP.profile(vx, vy, index);
    fstream fstrProfiles;
    fstrProfiles.open("profil.txt", ios::out);
    for(int i=0; i< vx.size(); i++){
      fstrProfiles << exp(vx.at(i)) << " " << exp(vy.at(i)) << endl;
    }
    FreemanChain::const_iterator it = fc.begin();
    for (int i=0; i< index; i++){
      ++it;
    }  
    Vector2i pt ((*it).x(), (*it).y());
    DrawingXFIG::drawCross (cout,pt ,4, 10, 10, 40);
    
    uint error =  MP.maxMeaningfulScale(index , 1, 0.0);
    cerr << "error scale for maxMeaningfulScales 1 0.0=" <<  error << endl;
    
    
  }


  if(args.check("-drawDetailedProfil")){
    uint index = args.getOption("-drawDetailedProfil")->getIntValue(0);
    vector <double> vx ;
    vector <double> vy ;
    vector <uint> nb;

    MP.detailedProfile(vx, vy, nb, index);
    
    vector <double> mx ;
    vector <double> my ;
    
    MP.profile(mx, my, index);
    
    fstream fstrProfiles;
    fstrProfiles.open("detailedProfil.txt", ios::out);
    uint pos= 0;
    for(int i=0; i< nb.size(); i++){
      
      uint nbCnt = nb.at(i);
      for (int j=0; j < nbCnt; j++){
	fstrProfiles << exp(vx.at(pos)) << " " << exp(vy.at(pos)) << endl;
	pos++;
      }
      
    }
    
    FreemanChain::const_iterator it = fc.begin();
    for (int i=0; i< index; i++){
      ++it;
    }  
    Vector2i pt ((*it).x(), (*it).y());
    DrawingXFIG::drawCross (cout,pt ,4, 10, 10, 40);
    
  }



  if(args.check("-plotReg")){
    uint index = args.getOption("-plotReg")->getIntValue(0);
    uint minSize = args.getOption("-plotReg")->getIntValue(1);    
    double alpha = args.getOption("-plotReg")->getDoubleValue(2);
    
    fstream fstrProfiles;
    fstrProfiles.open("regProfil.txt", ios::out);
    
    

  }  
  

  if(args.check("-drawContourSubsampled")){
    uint color = args.getOption("-drawContourSubsampled")->getIntValue(0);
    uint linewidth = args.getOption("-drawContourSubsampled")->getIntValue(1);
    MultiscaleFreemanChain::SubsampledChainKey key( samplingSizeMax, samplingSizeMax, 3, 3 );
    const MultiscaleFreemanChain::SubsampledChain* subchain = MP.get( key );
    DrawingXFIG::    drawContourPixelsCorner(cout, subchain->subc, 7,samplingSizeMax , 0 ,0 , 100);
    //
    if(args.check("-drawContourSRCOnly")){
      return 0;      
    }
  }
  if(args.check("-drawContourSubsampledK")){
    uint color = args.getOption("-drawContourSubsampledK")->getIntValue(0);
    uint linewidth = args.getOption("-drawContourSubsampledK")->getIntValue(1);
    uint resol = args.getOption("-drawContourSubsampledK")->getIntValue(2);
    MultiscaleFreemanChain::SubsampledChainKey key( resol, resol, 0, 0 );
    const MultiscaleFreemanChain::SubsampledChain* subchain = MP.get( key );
    DrawingXFIG::drawContourPixelsCorner(cout, subchain->subc, 7,resol ,0 ,0 , 100);
    //
    if(args.check("-drawContourSRCOnly")){
      return 0;      
    }
  }



  
  uint i =0;
  for (FreemanChain::const_iterator it = fc.begin() ; it!=fc.end(); ++it){
    Vector2i pt ((*it).x(), (*it).y());
    
    uint noiseLevel =0;
    
    if(args.check("-meaningfulScales")){
      noiseLevel =  MP.noiseLevel(i, mscales_min_size, mscales_max_slope);
    }else if(args.check("-maxValidScale")){
      double maxSlope = args.getOption("-maxValidScale")->getDoubleValue(0);
      noiseLevel = MP.maximalValidScale(i, maxSlope);
    }else if(args.check("-maxMeaningfulScales")){
      uint minSize = args.getOption( "-maxMeaningfulScales" )->getIntValue( 0 );
      double maxSlope = args.getOption( "-maxMeaningfulScales" )->getDoubleValue( 1 );
      noiseLevel = MP.maxMeaningfulScale (i,  minSize, maxSlope );
    }    
    
    
    DrawingXFIG::setFillIntensity(38);
    DrawingXFIG::drawPixel (cout, pt, 1, 1, noiseLevel, 50);
    i++;
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
    


  if(args.check("-displayNoisePPM")){    
    ImagePPM imageBruit  = importPPM("resuBruit.ppm");
    ImagePPM imageCompteur  = importPPM("resuBruitCompteur.ppm");
    
    uint maxBruitCpt = 0;
    uint maxBruit = 0;
    uint i =0;
    if(args.check("-initDisplayNoise")){
      for(uint i=0; i< imageBruit.height; i++){
	for(uint j=0; j < imageBruit.width; j++){
	  imageBruit.tabImageR[i*imageBruit.width+j]=0;
	  imageBruit.tabImageG[i*imageBruit.width+j]=0;
	  imageBruit.tabImageB[i*imageBruit.width+j]=0;
	  imageCompteur.tabImageR[i*imageBruit.width+j]=0;
	  imageCompteur.tabImageG[i*imageBruit.width+j]=0;
	  imageCompteur.tabImageB[i*imageBruit.width+j]=0;	  
	}      
      }
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























//---------------------------
// Sauvegarde format PGM/PPM 




ImagePGM importPGM(string name){

    
  ifstream in;
  in.open(name.c_str(), ifstream::in);
    
  
  ImagePGM imageImported;
  string str;
  getline( in, str );
  if ( ! in.good() ) return imageImported;
  if ( str != "P5" ) return imageImported;
  do
    {
      getline( in, str );
      if ( ! in.good() ) return imageImported;
    }
  while ( str[ 0 ] == '#' );
  istringstream str_in( str );
  str_in >> imageImported.width >> imageImported.height;
  in >> noskipws;
  getline( in, str );
  istringstream str2_in( str );
  int max_value;
  
  str2_in >> max_value;
  imageImported.maxValue=max_value;
  imageImported.tabImage = new uint[imageImported.width*imageImported.height];  
  
  for(int i=0; i<imageImported.width*imageImported.height; i++){
    unsigned char c; 
    in >> c;
    imageImported.tabImage[i] = (uint)c;    
  }  
  return imageImported;
}







ImagePPM importPPM(string name){

  ifstream in;
  in.open(name.c_str(), ifstream::in);

  ImagePPM imageImported;
  string str;
  getline( in, str );
  if ( ! in.good() ) return imageImported;
  if ( str != "P6" ) return imageImported;
  do
    {
      getline( in, str );
      if ( ! in.good() ) return imageImported;
    }
  while ( str[ 0 ] == '#' );
  istringstream str_in( str );
  str_in >> imageImported.width >> imageImported.height;
  in >> noskipws;
  getline( in, str );
  istringstream str2_in( str );
  int max_value;
  
  str2_in >> max_value;
  imageImported.maxValue=max_value;
  imageImported.tabImageR = new uint[imageImported.width*imageImported.height];  
  imageImported.tabImageG = new uint[imageImported.width*imageImported.height];  
  imageImported.tabImageB = new uint[imageImported.width*imageImported.height];  
  
  for(int i=0; i<imageImported.width*imageImported.height; i++){
    unsigned char c; 
    in >> c;
    imageImported.tabImageR[i] = (uint)c;    
    in >> c;
    imageImported.tabImageG[i] = (uint)c;    
    in >> c;
    imageImported.tabImageB[i] = (uint)c;    
  }  
  in.close();
  return imageImported;
}




void exportPGM(string name, ImagePGM image){
  ofstream out;
  out.open(name.c_str(), ios_base::out);

  out << "P5" << endl
      << "# CREATOR: addGaussianNoisePGM " 
      << "(kerautret@loria.fr)" << endl;
  out << image.width << " " << image.height << endl
      << image.maxValue << endl;
  
  
  for(int i=0; i<image.width*image.height; i++){
    out << (unsigned char)(image.tabImage[i]) ;    
  }
  
  out << endl;
  out.close();
}





void exportPPM(string name, ImagePPM image){
  ofstream out;
  out.open(name.c_str(), ios_base::out);
    
  out << "P6" << endl
      << "# CREATOR: displayNoise " 
      << "(kerautret@loria.fr)" << endl;
  out << image.width << " " << image.height << endl
      << image.maxValue << endl;
  
  
  for(int i=0; i<image.width*image.height; i++){
    out << (unsigned char)(image.tabImageR[i])<< (unsigned char)(image.tabImageG[i])
	<< (unsigned char)(image.tabImageB[i]);    
  }
  out << endl;
  out.close();
}










