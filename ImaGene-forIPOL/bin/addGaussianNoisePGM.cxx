///////////////////////////////////////////////////////////////////////////////
// Add Gaussian noise to a pgm image
///////////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <sstream>

#include "ImaGene/base/Arguments.h"
#include "ImaGene/base/Vector.h"
#include "ImaGene/base/Shapes.h"
#include "ImaGene/digitalnD/GridEmbedder.h"
#include "ImaGene/digitalnD/KnSpace.h"
#include "ImaGene/digitalnD/KnCharSet.h"
#include "ImaGene/digitalnD/KnShapes.h"
#include "ImaGene/digitalnD/C4CIteratorOnSurface.h"
#include "ImaGene/digitalnD/BelAdjacency.h"
#include "ImaGene/digitalnD/ObjectBoundary.h"
#include "ImaGene/digitalnD/Frame2D.h"
#include "ImaGene/helper/ContourHelper.h"
#include "ImaGene/helper/ShapeHelper.h"

using namespace std;
using namespace ImaGene;


static Arguments args;



struct ImagePGM{
  int width;
  int height;
  int maxValue;
  int *tabImage;
};


void exportPGM(ostream & out, ImagePGM image);
ImagePGM importPGM( std::istream & in);
double randUniforme (double min, double max);
double randomGaussien (double moyenne, double ecartType);




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
  srand ( time(NULL) );

  args.addOption( "-stdDev", "-stdDev <val>: Gaussian noise with std Dev <val>", "10" );
  args.addOption( "-stdDevSector", "-stdDevSector <std Dev sector1> <std Dev sector2> <std Dev sector3> <std Dev sector4>: Gaussian noise with std Dev <val>", "10", "10", "10", "10" );
  
  if ( ( argc <= 0 ) 
       || ! args.readArguments( argc, argv ) ) 
    {
      cerr << args.usage( "addGaussianNoisePGM", 
			  "Add Gaussian noise to the grey level of the image given as pgm in standard input to the standard output. The noise distribution is centered to 0.",
			  "-stdDev -stdDevSector" )
	   << endl;
      return 1;
    }


  if(args.check( "-stdDev" )){
    double deviation =  args.getOption( "-stdDev" )->getFloatValue( 0 );    
    ImagePGM image = importPGM(cin);
    for(int i=0; i<image.height*image.width; i++){
      image.tabImage[i] = image.tabImage[i] + randomGaussien (0, deviation);
      if(image.tabImage[i]<0)
      image.tabImage[i]=0;
      if(image.tabImage[i]>image.maxValue)
	image.tabImage[i]=image.maxValue;
    }
    exportPGM(cout,image);
    delete image.tabImage;
  
  }

  if(args.check( "-stdDevSector" )){
    double dev1 =  args.getOption( "-stdDevSector" )->getFloatValue( 0 );    
    double dev2 =  args.getOption( "-stdDevSector" )->getFloatValue( 1 );    
    double dev3 =  args.getOption( "-stdDevSector" )->getFloatValue( 2 );    
    double dev4 =  args.getOption( "-stdDevSector" )->getFloatValue( 3 );    
    
    double deviation=0.0;
    ImagePGM image = importPGM(cin);
    for(int i=0; i<image.height; i++){
      for(int j=0; j<image.width; j++){
	if((i<image.width/2) && (j<image.height/2)){
	  deviation=dev1;	  
	}else if((i>=image.width/2)&(j<image.height/2)){
	  deviation=dev2;	  
	}else if((i<image.width/2) && (j>=image.height/2)){
	  deviation=dev4;	  
	}else if((i>=image.width/2)&(j>=image.height/2)){
	  deviation=dev3;	  
	}
	
	image.tabImage[i*image.width+j] = image.tabImage[i*image.width+j] + randomGaussien (0, deviation);
	if(image.tabImage[i*image.width+j]<0)
	  image.tabImage[i*image.width+j]=0;
	if(image.tabImage[i*image.width+j]>image.maxValue)
	  image.tabImage[i*image.width+j]=image.maxValue;
      }
    }
    exportPGM(cout,image);
    delete image.tabImage;
  
  }
  
  

  
  return 0;
}
 






ImagePGM importPGM( std::istream & in){
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
  while ( str[ 0 ] == '#' || str=="" );
  istringstream str_in( str );
  str_in >> imageImported.width >> imageImported.height;
  in >> noskipws;
  getline( in, str );
  istringstream str2_in( str );
  int max_value;
  
  str2_in >> max_value;
  imageImported.maxValue=max_value;
  imageImported.tabImage = new int[imageImported.width*imageImported.height];  
  
  for(int i=0; i<imageImported.width*imageImported.height; i++){
    unsigned char c; 
    in >> c;
    imageImported.tabImage[i] = (int)c;    
  }  
  return imageImported;
}






void exportPGM(ostream & out, ImagePGM image){
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






double
randUniforme (double min, double max){
  double r;
  r = (double) rand() / (double) RAND_MAX;
  return ( r*(max-min) + min);
}




double
randomGaussien (double moyenne, double ecartType){
  double r1 = randUniforme (0.0, 1.0);
  double r2 = randUniforme (0.0, 1.0);
  double r = sqrt (-2.0 * log(r1)) * cos (2.0 * M_PI * r2);
  return (moyenne + ecartType * r);
}



