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




void exportPGM(const char* name, ImagePGM image);
ImagePGM importPGM(const char * name);


void exportPPM(const char* name, ImagePPM image);
ImagePPM importPPM( const  char* name);


Vector2i getOptimalPosition(const ImagePGM &imgRef, const ImagePGM &img, int backgroundValue, int dx, int dy);



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
 
  
  
  args.addOption("-srcA", "-srcA <sourceImage>" , " ");
  args.addOption("-srcB", "-srcB <sourceImage>" , " ");
  args.addBooleanOption("-quiet", "-quiet");
  
  
  if ( ( argc <= 1 )|| ! args.readArguments( argc, argv ))
    {
      cerr << args.usage( "ajustImages", 
			  "ajustImages -srcA <image1.ppm> -srcB <image2.ppm> : return ."
			  ,"" ) << endl;
      return 1;
    }

  
  string str1 = args.getOption("-srcA")->getValue(0);
  string str2 = args.getOption("-srcB")->getValue(0);
  
  ImagePGM image1 =  importPGM(str1.c_str());
  ImagePGM image2 =  importPGM(str2.c_str());
  
  
  
  Vector2i ptOptimal = getOptimalPosition(image1, image2, image1.maxValue, 20, 20);
  if(!args.check("-quiet")){
    cerr << "Position optimal: " << ptOptimal.x() << " " << ptOptimal.y() << endl;
  }

  
  
  // image 1 : modèle
  
  ImagePGM imgAjusted ;
  imgAjusted.width = image1.width;
  imgAjusted.height = image1.height;  
  imgAjusted.tabImage = new uint[image1.width*image1.height];
  imgAjusted.maxValue = image1.maxValue;
  for(int i=0; i< imgAjusted.height; i++){
    for(int j=0; j< imgAjusted.width; j++){
      int posNewX = j - ptOptimal.x();
      int posNewY = i - ptOptimal.y();
      int valRef = (posNewX<0 || posNewX>= image2.width || posNewY<0 || posNewY>= image2.height )? image1.maxValue : 
	(int)(image2.tabImage[posNewY*image2.width +posNewX]);      
      imgAjusted.tabImage[i*imgAjusted.width+j]=valRef;
    }    
  }
  
  exportPGM("imgAjusted.pgm", imgAjusted);
  
  
}
  
  

  

	        








//---------------------------
// Sauvegarde format PGM/PPM 




ImagePGM importPGM(const  char * name){

    
  ifstream in;
  in.open(name, ifstream::in);
    
  
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







ImagePPM importPPM(const  char* name){

  ifstream in;
  in.open(name, ifstream::in);

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
  return imageImported;
}




void exportPGM(const char* name, ImagePGM image){
  ofstream out;
  out.open(name, ios_base::out);

  out << "P5" << endl
      << "# CREATOR: addGaussianNoisePGM " 
      << "(kerautret@loria.fr)" << endl;
  out << image.width << " " << image.height << endl
      << image.maxValue << endl;
  
  
  for(int i=0; i<image.width*image.height; i++){
    out << (unsigned char)(image.tabImage[i]) ;    
  }
  
  out << endl;
}





void exportPPM(const char* name, ImagePPM image){

  ofstream out;
  out.open(name, ios_base::out);
    
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
}





 
// Recherche de la position la meilleure minimisant la différence entre les deux images.
// img1: image de ref
// img2: image 

Vector2i
getOptimalPosition(const ImagePGM &imgRef, const ImagePGM &img, int backgroundValue, int dx, int dy){
  Vector2i resu(0,0);
  int widthRef = (int) imgRef.width;
  int heightRef = (int) imgRef.height;

  
  int widthImg = (int) img.width;
  int heightImg = (int) img.height;
  bool first=true;
  int erreurMin;
  
  // pTestX,Y postion de l'image à re positionner
  for(int pTestY=-dy; pTestY + heightImg < heightRef+dy; pTestY++ ){
    for(int pTestX=-dy; pTestX +widthImg < widthRef+dx; pTestX++ ){
      int erreurTot=0;

      // Evaluation de l'erreur:
      for(int i =0; i < heightImg; i++){
	for(int j =0; j < widthImg; j++){
	  int yref= pTestY + i;
	  int xref= pTestX + j;	  
	  int indice = i*widthImg+j;
	  int valImg = (int)( img.tabImage[indice]);
	  int valRef = (xref<0 || xref>= widthRef || yref<0 || yref>= heightRef )? backgroundValue : 
	    (int)(imgRef.tabImage[yref*widthRef +xref]);
	  
	  if(valImg!=valRef)
	    erreurTot++;
	}
      }
      if(first) {
	first=false;
	erreurMin=erreurTot;	
      }
      if(erreurTot<=erreurMin){
	erreurMin = erreurTot;
	resu.x() = pTestX;
	resu.y() = pTestY;	
      }
	
    }
  }
  if(args.check("-quiet")){
    cout <<erreurMin;
  }else{
    cout << "erreur OPT = " << erreurMin << endl;
  } 
  return resu;
}





