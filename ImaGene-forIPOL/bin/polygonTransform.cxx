///////////////////////////////////////////////////////////////////////////////
// Add Gaussian noise to a pgm image
///////////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
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

using namespace std;
using namespace ImaGene;


static Arguments args;
#include "ImaGene/base/StandardArguments.h"



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
  StandardArguments::addIOArgs( args, true, false );

  
  //args.addOption( "-stdDev", "-stdDev <val>: Gaussian noise with std Dev <val>", "10" );

  args.addOption("-addGaussianNoise", "-addGaussianNoise <stdDev> add gaussian noise by moving  each vertex of the polygon with gaussian noise of <stdDev>", "5.0");
  args.addOption("-addGaussianNoiseSect", "-addGaussianNoiseSect <stdDev1> <stdDev2> <stdDev3> <stdDev4> add gaussian noise by moving  each vertex of the polygon with gaussian noise of <stdDev1>, <stdDev2>, <stdDev3>, <stdDev4> according the vertex position ", "0.0", "0.2", "0.4", "0.6" );
  
  args.addOption("-overSampling", "-overSampling <nb> add nb iso points between each vertex", "1" );
  args.addOption("-reduceSampling", 
		 "-reduceSampling <samplingRate> reduce the number of polygon point by removing all point of index i for which i%samplingRate!=0. ", "5");
  args.addOption("-reduceSamplingSect", 
		 "-reduceSamplingSect <samplingRate1> <samplingRate2> <samplingRate3> <samplingRate4> reduce the number of polygon point by removing all point of index i for which i%samplingRate!=0. ", "1", "4", "8", "16");
  
  args.addOption("-changeScaleFactor", "-changeScaleFactor <scale> ","1.0"); 
  
  args.addBooleanOption("-minEuclideanSpanningTree", "-minEuclideanSpanningTree deduce the contour as the largest path of the Euclidean Spanning Tree of the initial points set." );
  args.addBooleanOption("-transformLoop", "-transformLoop " );
  args.addBooleanOption("-removeLoop", "-removeLoop " );
  args.addOption("-savePointSet", "savePointSet <filename>", " ");
  


  


  
  if ( ( argc <= 0 ) 
       || ! args.readArguments( argc, argv ) ) 
    {
      cerr << args.usage( "polygonTransform", 
			  "Transform the polygon given as input ( X Y field)  with several possible transformations .",
			  "" )
	   << endl;
      return 1;
    }

  


  istream & in_str = StandardArguments::openInput( args );
  vector<Vector2D> polygon = ContourHelper::getPolygonFromStream(in_str);


  
  if(args.check("-reduceSampling")){
    uint samplingRate = args.getOption("-reduceSampling")->getIntValue(0);
    vector<Vector2D> newContour;
    for(int i =0; i< polygon.size(); i++){
      if(i%samplingRate==0){
	newContour.push_back(polygon.at(i));
      }
    }
    polygon=newContour;
  }


  if(args.check("-reduceSamplingSect")){
    
    uint samplingRate1 = args.getOption("-reduceSamplingSect")->getIntValue(0);
    uint samplingRate2 = args.getOption("-reduceSamplingSect")->getIntValue(1);
    uint samplingRate3 = args.getOption("-reduceSamplingSect")->getIntValue(2);
    uint samplingRate4 = args.getOption("-reduceSamplingSect")->getIntValue(3);
    Vector2D ptMin, ptMax;
    ContourHelper::computeContourBoundingBox(polygon, ptMin, ptMax);
    uint samplingRate;
    vector<Vector2D> newContour;
    for(int i =0; i< polygon.size(); i++){
      Vector2D p = polygon.at(i);
      if((p.x()<(ptMin.x()+(ptMax.x()-ptMin.x())/2.0))){
	if((p.y()<(ptMin.y()+(ptMax.y()-ptMin.y())/2.0))){
	  samplingRate=samplingRate1;
	}else{
	  samplingRate=samplingRate2;
	}
      }else{
	if((p.y()<(ptMin.y()+(ptMax.y()-ptMin.y())/2.0))){
	  samplingRate=samplingRate4;
	}else{
	  samplingRate=samplingRate3;
	}
      }
      
      if(i%samplingRate==0){
	newContour.push_back(polygon.at(i));
      }
    }
    polygon=newContour;
  }

  

  
  if(args.check("-overSampling")){
    
    uint nbAdd = args.getOption("-overSampling")->getIntValue(0);    
    vector<Vector2D> newContour;
    for(int i =0; i< polygon.size(); i++){
      Vector2D ptA = polygon.at(i);
      Vector2D ptB = polygon.at((i+1)%polygon.size());
      newContour.push_back(ptA);
      for (int k=0; k<nbAdd; k++){
	Vector2D ptNew((ptA.x()/(nbAdd+1.0))*(k+1)+(ptB.x()/(nbAdd+1.0))*(nbAdd-k),
		       (ptA.y()/(nbAdd+1.0))*(k+1)+(ptB.y()/(nbAdd+1.0))*(nbAdd-k) );
	newContour.push_back(ptNew);
      }
      newContour.push_back(ptB);
    }
    polygon=newContour;
  }
  

  
  


  if(args.check("-addGaussianNoise")){
    double deviation = args.getOption("-addGaussianNoise")->getFloatValue(0);
    vector<Vector2D> newContour;
    for(int i =0; i< polygon.size(); i++){
      Vector2D pt = polygon.at(i);
      newContour.push_back(Vector2D(pt.x()+randomGaussien (0, deviation),pt.y()+randomGaussien (0, deviation) ));
      
    }
    polygon=newContour;
  }
  




  if(args.check("-addGaussianNoiseSect")){
    
    double deviation1 = args.getOption("-addGaussianNoiseSect")->getFloatValue(0);
    double deviation2 = args.getOption("-addGaussianNoiseSect")->getFloatValue(1);
    double deviation3 = args.getOption("-addGaussianNoiseSect")->getFloatValue(2);
    double deviation4 = args.getOption("-addGaussianNoiseSect")->getFloatValue(3);
    double deviation;
    Vector2D ptMin, ptMax;
    ContourHelper::computeContourBoundingBox(polygon, ptMin, ptMax);
    uint samplingRate;
    vector<Vector2D> newContour;
    for(int i =0; i< polygon.size(); i++){
      Vector2D p = polygon.at(i);
      if((p.x()<(ptMin.x()+(ptMax.x()-ptMin.x())/2.0))){
	if((p.y()<(ptMin.y()+(ptMax.y()-ptMin.y())/2.0))){
	  deviation=deviation1;
	}else{
	  deviation=deviation2;
	}
      }else{
	if((p.y()<(ptMin.y()+(ptMax.y()-ptMin.y())/2.0))){
	  deviation=deviation4;
	}else{
	  deviation=deviation3;
	}
      }
      
      newContour.push_back(Vector2D(p.x()+randomGaussien (0, deviation),p.y()+randomGaussien (0, deviation) ));
      
    }
    polygon=newContour;
  }

  
  
  
  if(args.check("-savePointSet")){
    string filename = args.getOption("-savePointSet")->getValue(0);
    ofstream outfile;
    outfile.open(filename.c_str());
    outfile<< "# Output from polygonTransform (unordered noisy point set) " <<endl;
    for(int i=0;i< polygon.size(); i++){
      outfile << polygon.at(i).x() << " " << polygon.at(i).y() << endl;
    }
    outfile.close();
  }

  
  
  
  if(args.check("-transformLoop")||args.check("-removeLoop")){
    polygon=ContourHelper::transformLoop(polygon,args.check("-removeLoop"),  0.1);
  }

  
  if(args.check("-minEuclideanSpanningTree")){
    ContourHelper::ValuedEdgeSet g= ContourHelper::getMinEuclideanSpanningTree(polygon);
     vector<uint> longPath =ContourHelper::getLongestPathFromTree(g, polygon.size());
     vector<Vector2D> newContour;
     for(int i=0; i< longPath.size()-1; i++){
       newContour.push_back(polygon.at(longPath.at(i)));
     }
     polygon = newContour;
  }



  
  double scale = args.getOption("-changeScaleFactor")->getFloatValue(0);

  for(int i =0; i< polygon.size(); i++){
    cout << polygon.at(i).x()*scale << " " << polygon.at(i).y()*scale << endl;
  }

  
  


  
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



