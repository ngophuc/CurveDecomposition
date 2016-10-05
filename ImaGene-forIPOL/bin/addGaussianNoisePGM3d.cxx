///////////////////////////////////////////////////////////////////////////////
// Add Gaussian noise to a pgm image
///////////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <sstream>
#include "ImaGene/base/StandardArguments.h"
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
#include "ImaGene/digitalnD/KnSpaceScanner.h"
#include "ImaGene/digitalnD/KnRCellVector.h"

using namespace std;
using namespace ImaGene;


static Arguments args;





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
  
  StandardArguments::addIOArgs( args, true, true );
  args.addOption( "-stdDevSector", "-stdDevSector <std Dev sector1> <std Dev sector2> <std Dev sector3> <std Dev sector4>: \
		   Gaussian noise with std Dev <val> (default 10, 10, 10 ,10)", "10", "10", "10", "10" );
  args.addOption( "-stdDev", "-stdDev <std Dev> Gaussian noise with std Dev <val> (default 10)", "10" );
  args.addBooleanOption("-binary", "-binary used to consider binary image (0-1) by transforming with range 0-255"); 
  
    if ( ( argc <= 0 ) 
       || ! args.readArguments( argc, argv ) ) 
    {
      cerr << args.usage( "addGaussianNoisePGM3d", 
			  "Add Gaussian noise to the grey level of the image given as pgm in standard input to the standard output. The noise distribution is centered to 0.",
			  "-stdDev -stdDevSector -binary" )
	   << endl;
      return 1;
    }

  


  KnSpace *ks;
  KnRUCellVector <int> *vectorVal3d;


  istream & in_str = StandardArguments::openInput( args );
  ShapeHelper::importFromPGM3d(in_str, ks, vectorVal3d, 3);
  uint largeur= ks->size(0);
  uint hauteur= ks->size(1);
  uint prof = ks->size(2);
    
  
  


    
    double dev1 =  args.getOption( "-stdDevSector" )->getFloatValue( 0 );    
    double dev2 =  args.getOption( "-stdDevSector" )->getFloatValue( 1 );    
    double dev3 =  args.getOption( "-stdDevSector" )->getFloatValue( 2 );    
    double dev4 =  args.getOption( "-stdDevSector" )->getFloatValue( 3 );    
    
    double deviation=args.getOption( "-stdDev" )->getFloatValue( 0 );
    

    int x,y,z;
    KnSpaceScanner scan2( *ks, 
			  ks->ufirstCell( ks->dim() ),
			  ks->ulastCell( ks->dim() ) );
    
    Kn_uid p = scan2.begin();
    for ( Kn_uid last_z = scan2.last( p, 2 ), z=0;
	  p <= last_z; 
	  z++,p += scan2.gotonext( 2 ) )
      {  
	for (Kn_uid last_y = scan2.last( p, 1 ),y=0;
	     p <= last_y; 
	     y++, p += scan2.gotonext( 1 ) )
	  {
	    for ( Kn_uid last_x = scan2.last( p, 0 ),x=0; 
		  p <= last_x; 
		  x++, p++ ) // NB: 'scan.gotonext( 0 )' == 1;
	      {
		
		
		if(args.check( "-stdDevSector" )){
		  if((x<largeur/2) && (y<hauteur/2)&&(z<prof/2)){
		    deviation=dev1;	  
		  }else if((x>=largeur/2)&(y<hauteur/2)&&(z<prof/2)){
		    deviation=dev2;	  
		  }else if((y<largeur/2) && (y>=hauteur/2)&&(z<prof/2)){
		    deviation=dev4;	  
		  }else if((x>=largeur/2)&(y>=hauteur/2)&&(z<prof/2)){
		    deviation=dev3;	  
		  }else if((x<largeur/2) && (y<hauteur/2)&&(z>=prof/2)){
		    deviation=dev4;	  
		  }else if((x>=largeur/2)&(y<hauteur/2)&&(z>=prof/2)){
		    deviation=dev3;	  
		  }else if((y<largeur/2) && (y>=hauteur/2)&&(z>=prof/2)){
		    deviation=dev1;	  
		  }else if((x>=largeur/2)&(y>=hauteur/2)&&(z>=prof/2)){
		    deviation=dev2;	  
		  }
		}
		
		int valImage = (*vectorVal3d)[p];
		if (args.check("-binary")){
		  valImage = (valImage==1)? 255:0;
		}
		int val = valImage+randomGaussien (0, deviation);
		if(val>255)
		  (*vectorVal3d)[p]= 255;
		else if(val<0)
		  (*vectorVal3d)[p]= 0;		
		else
		  (*vectorVal3d)[p]=val;
	      }
	  }
      }
    
    
    ShapeHelper::exportToPGM3d(cout,  ks, *vectorVal3d);
  
    
  
  
  
  return 0; 
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



