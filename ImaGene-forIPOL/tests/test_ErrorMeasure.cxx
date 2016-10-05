//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <algorithm>
#include "ImaGene/base/Vector2i.h"
#include "ImaGene/base/Arguments.h"
#include "ImaGene/base/StandardArguments.h"
#include "ImaGene/mathutils/Mathutils.h"
#include "ImaGene/mathutils/Statistics.h"
#include "ImaGene/timetools/Clock.h"
#include "ImaGene/dgeometry2d/FreemanChain.h"
#include "ImaGene/dgeometry2d/FreemanChainTransform.h"
#include "ImaGene/dgeometry2d/C4CIteratorOnFreemanChain.h"
#include "ImaGene/digitalnD/C4CIteratorOnFreemanChainSurface.h"
#include "ImaGene/helper/DrawingXFIG.h"

#include "ImaGene/helper/ContourHelper.h"


using namespace std;
using namespace ImaGene;

static Arguments args;


vector<double> lambdaMSTEstimator (FreemanChain & fc);



void 
computeErrorStat(vector<Vector2D> polygon, 
		 vector<Vector2i> contour, vector<double> vectTangentesSRC,
		 vector<int> vectAssociations,
		 Statistics & stat );

template<class T> 
vector<int> 
getAssociations(const vector<Vector2i> &contour, const vector<T> &polygon);


#define INFTY 10000001
 

double d_euc ( unsigned int i, unsigned int j, vector<Vector2D> &P, vector<Vector2D> &Q )
 {
   return sqrt ( ( P[i].x() - Q[j].x() ) * ( P[i].x() - Q[j].x() ) +
   ( P[i].y() - Q[j].y() ) * ( P[i].y() - Q[j].y() ) );
 }
 
double c ( unsigned int i, unsigned int j,
 	 vector<Vector2D> &P, vector<Vector2D> &Q,
 	 double *ca )
{
  if ( ca[i+j*P.size() ] > -1 )
    return ca[i+j*P.size() ];
  else
    if ( ( i ==0 ) && ( j==0 ) )
      ca[0] = d_euc( 0,0,P,Q );
    else
      if ( ( i>0 ) && ( j==0 ) )
	ca[i+j*P.size() ] = max ( c ( i-1,0,P,Q,ca ),d_euc ( i,0,P,Q ) );
      else
	if ( ( i==0 ) && ( j>0 ) )
	  ca[i+j*P.size() ] = max ( c ( 0,j-1,P,Q,ca ),d_euc( 0,j,P,Q ) );
	else
	  if ( ( i>0 ) && ( j>0 ) )
	    ca[i+j*P.size() ] = max ( min ( c ( i-1,j,P,Q,ca ),
				min ( c ( i-1,j-1,P,Q,ca ),c ( i,j-1,P,Q,ca ) ) ),d_euc ( i,j,P,Q ) );
	  else
	    ca[i+j*P.size() ] = INFTY;
 
   return ca[ i+j*P.size() ];
 }
 
double discreteFrechet ( vector<Vector2D> &P, vector<Vector2D> &Q )
{
  double *ca;

  ca = (double *)malloc(sizeof(double)*P.size() *Q.size() );

  if (ca==NULL)
  {
    cout<< "Allocation error!"<<endl;
    exit(2);
  }
    
  for ( unsigned int i=0; i < P.size(); i++ )
    for ( unsigned int j=0; j < Q.size(); j++ )
      ca[i + j*P.size() ] = -1;
 
 return c ( P.size()-1,Q.size()-1,P,Q,ca );

}
 

 

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
 
  args.addOption("-polygon", "-polygon <filename> <posX> <posY> the polygon file containing X Y coordinates at positions posX , posY (converted to CCW is necessary)", " ", "0", "1");
 args.addOption("-noisyContour", "-noisyContour <noisyContour.fc> freeman chain of the noisy contour (should be CCW)", " ");
 args.addOption("-sourceContour", "-sourceContour <sourceContour.fc> freeman chain of the source contour (should be CCW)", " ");
 
 args.addOption("-sourceContourGeom", "-sourceContourGeom <fichierGeom.txt> <pos x> <pos y> <pos tgt> specify the contour and normal for the geometric file the position <position> ", " ", "1", "2", "3");
 
 if ( ( argc <= 1 ) || ! args.readArguments( argc, argv ))
 {
 cerr << args.usage( "test_ErrorMeasure",
 "Error measure testing..."
 	 ,"-polygon -noisyContour -sourceContour" ) << endl;
 return 1;
 }
 

 
 // Recup polygon...
 string filePoly = args.getOption("-polygon")->getValue(0);
 int posX = args.getOption("-polygon")->getIntValue(1);
 int posY = args.getOption("-polygon")->getIntValue(2);
 fstream fstPoly;
 fstPoly.open (filePoly.c_str(), ios::in);
 vector<Vector2D> vectPolygon = ContourHelper::getPolygonFromStream(fstPoly, posX, posY);
 cerr << "Polygon CCW " << ((ContourHelper::isCCW(vectPolygon))? "[OK]" :"[False]" )<< endl;
 if(ContourHelper::transformCCW(vectPolygon))
 cerr << "change to CCW ..."<< endl;
 
 vector<double> vectTgt;
 
 // Recup contour initial + tangentes
 vector<Vector2i> vectContourSRC;
 FreemanChain fcSrc;
 if(!args.check("-sourceContourGeom")){
   string fileSrc = args.getOption("-sourceContour")->getValue(0);
   fstream fstSrc;
   fstSrc.open(fileSrc.c_str(), ios::in);
   FreemanChain::read(fstSrc, fcSrc);
   if ( ! fstSrc.good() ){
     cerr << "Error reading Freeman chain code." << endl;
     return 2;
   }
   FreemanChain::getContourPoints(fcSrc, vectContourSRC);
   if(ContourHelper::isCCW(vectContourSRC)){
 cerr << "Initial contour counter CCW [OK]" <<endl;
   }else{
     cerr << "Initial contour should be CCW."<< endl;
     exit(2);
 }
   // Estimation des tangentes par L-MST
   vectTgt = lambdaMSTEstimator(fcSrc);
 }else{
   string fileSrcGeom = args.getOption("-sourceContourGeom")->getValue(0);
   uint indx= args.getOption("-sourceContourGeom")->getIntValue(1);
   uint indy= args.getOption("-sourceContourGeom")->getIntValue(2);
   uint indTgt= args.getOption("-sourceContourGeom")->getIntValue(3);
   fstream fstGeom;
   cerr << "ind x " << indx << endl;
   cerr << "ind y " << indy << endl;
   fstGeom.open (fileSrcGeom.c_str(), ios::in);
   vector<double> vx = ContourHelper::extractField<double>(fstGeom,indx);
   fstGeom.close();
   fstGeom.open (fileSrcGeom.c_str(), ios::in);
   vector<double> vy = ContourHelper::extractField<double>(fstGeom,indy);
   fstGeom.close();
   cerr << "size=" << vx.size() << endl;
   cerr << "size=" << vy.size() << endl;
   for (int i=0;i< vx.size();i++){
     cerr <<"i = " << i<< endl;
     vectContourSRC.push_back(Vector2i(vx.at(i), vy.at(i)));
   }
   
   cerr << "ok" << vectContourSRC.size()<<endl;
   fstGeom.open (fileSrcGeom.c_str(), ios::in);
   vectTgt=ContourHelper::extractField<double>(fstGeom,indTgt);
 }
 
 // Recup contour bruité
 vector<Vector2i> vectContourNoisy;
 FreemanChain fcNoisy;
 string fileNoisy = args.getOption("-noisyContour")->getValue(0);
 fstream fstNoisy;
 fstNoisy.open(fileNoisy.c_str(), ios::in);
 FreemanChain::read(fstNoisy, fcNoisy);
 if ( ! fstNoisy.good() ){
   cerr << "Error reading Freeman chain code." << endl;
   return 2;
 }
 FreemanChain::getContourPoints(fcNoisy, vectContourNoisy);
 if(ContourHelper::isCCW(vectContourNoisy)){
   cerr << "Initial contourNoisy counter CCW [OK]" <<endl;
 }else{
   cerr << "Initial contour should be CCW."<< endl;
 }
 
 
 // Test affichage des contours
 DrawingXFIG::includeXFIGHeader(cout, 1024, 2);
 DrawingXFIG::drawContour (cout, vectContourSRC, 1, 1, 20, true, false, 30, 1.0, 0.0,0.0);
 // DrawingXFIG::drawContour (cout, vectContourNoisy, 2, 2, 20, true, false, 55, 1.0, 0.0,0.0);
 DrawingXFIG::drawContour (cout, vectPolygon, 0, 3, 20, true, false, 35, 1.0, 0.0,0.0);

 
 //Testing contours associations:
 vector<int> vectAsso = getAssociations(vectContourSRC, vectPolygon);
 cerr << "vectPolygone" << vectPolygon.size()<< endl;
 cerr << "getting associations [OK] " <<vectAsso.size()<< endl;
 
 // for(int i=0;i< vectPolygon.size();i++){
 //  Vector2D pol = vectPolygon.at(i);
 //  Vector2i pCnt = vectContourSRC.at(vectAsso.at(i));
 //  DrawingXFIG::drawLine(cout, Vector2i(pol.x()+1, pol.y()+1), pCnt, 3, 10, 0);
 //}


 Statistics st (2, false);
 computeErrorStat(vectPolygon, vectContourSRC, vectTgt,vectAsso, st);
 st.terminate();

  ///Frechet
  /// TODO: Ugly Vector2i->Vector2D conversion
  vector<Vector2D> vectContourSRC_float(vectContourSRC.size());
  vector<Vector2D> vectContourNoisy_float(vectContourNoisy.size());
  for(unsigned int i=0; i<vectContourSRC.size(); i++)
    vectContourSRC_float.push_back(Vector2D(vectContourSRC[i].x(),vectContourSRC[i].y()));
  
  
  for(unsigned int i=0; i<vectContourNoisy.size(); i++)
    vectContourNoisy_float.push_back(Vector2D(vectContourNoisy[i].x(),vectContourNoisy[i].y()));
  
  

  cerr<< "Total Points:" << st.samples(0) << endl<< 
    "Tangent error: mean: "<< st.mean(0) << " max: "<< st.max(0)<< " variance: " << st.variance(0)<<endl;
  cerr << "Distance error: mean: "<< st.mean(1) << " max: "<< st.max(1)<< " variance: " << st.variance(1)<<endl;
  

  cerr<< " discreteFrechet distance (SRC-Noisy) "<< discreteFrechet(vectContourSRC_float, vectContourNoisy_float)<<endl;
  cerr<< " discreteFrechet distance (Poly-SRC) "<< discreteFrechet(vectPolygon, vectContourSRC_float)<<endl;
  cerr<< " discreteFrechet distance (Poly-Noisy) "<< discreteFrechet(vectPolygon, vectContourNoisy_float)<<endl;

 return 0;

}





// Compute Error from tangent and distance to initial discrete contour

void 
computeErrorStat(vector<Vector2D> polygon, 
		 vector<Vector2i> contour, vector<double> vectTangentesSRC,
		 vector<int> vectAssociations,
		 Statistics & stat ){
  

  cerr << "Calcul de l'erreur sur les normales...";
  for (int i=0; i < polygon.size(); i++){
    Vector2i pt1 ((int) (polygon.at(i%polygon.size()).x()),
		  (int) (polygon.at(i%polygon.size()).y()));
    Vector2i pt2 ((int) (polygon.at((i+1)%polygon.size()).x()),
		  (int) (polygon.at((i+1)%polygon.size()).y()));
    int indicePt1Cnt=vectAssociations.at((i)%polygon.size());
    int indicePt2Cnt=vectAssociations.at((i+1)%polygon.size());
    
    Vector2D tangente ( (pt2.x() -pt1.x()), (pt2.y() - pt1.y()));
    double norm= sqrt(tangente.x()*tangente.x()+tangente.y()*tangente.y());
    tangente.x()/=norm;
    tangente.y()/=norm;
    for(int j = indicePt1Cnt; j!=indicePt2Cnt;j=(j+1)%contour.size()){
      double thetaMST= vectTangentesSRC.at(j%(vectTangentesSRC.size()));
      Vector2D tangenteMST (cos(thetaMST), sin(thetaMST));
      Vector2D p1 (contour.at(j%contour.size()).x(),contour.at(j%contour.size()).y());
      double xInter, yInter;
      if(pt2.y()== pt1.y()){
	xInter= p1.x();
	yInter = pt1.y();
      }else if(pt2.x()== pt1.x()){
	xInter= pt1.x();
	yInter = p1.y();
      }else{
	double a = tangente.y()/tangente.x();
	double b = pt1.y()-a*pt1.x();
	double ap = -1.0/a;
	double bp = p1.y()-ap*p1.x();
	xInter = -(b-bp)/(a-ap);
	yInter = a*xInter + b;
      }
      Vector2D p2 (xInter, yInter);
      DrawingXFIG::drawLine(cout, p1, p2, 3, 10, 0);
      
      double distance = sqrt((p1.x()-p2.x())*(p1.x()-p2.x())+(p1.y()-p2.y())*(p1.y()-p2.y()));
      stat.addValue(1, distance);
      double dx= tangenteMST.x() - tangente.x();
      double dy= tangenteMST.y() - tangente.y();
      double distanceTgt = sqrt(dx*dx+dy*dy);
      DrawingXFIG::drawCircle(cout, p1, 4, distanceTgt*1.0, 48);
      stat.addValue(0,distanceTgt);
    }
  }
  cerr<< "[OK]" <<endl;
}








vector<double>
lambdaMSTEstimator (FreemanChain & fc){
 
 KnSpace* ks = ShapeHelper::makeSpaceFromFreemanChain( fc );
 vector<double> vectResult;
 
 // Computes maximal segments.
 C4CIteratorOnFreemanChain itfc;
 itfc.init( fc.begin(), true );
 C4CIteratorOnFreemanChainSurface itfcs;
 
 itfcs.init( ks, itfc );
 C4CTangentialCover tcover;
 bool isInit = tcover.init( itfcs,0 );
 uint surfaceSize = tcover.nbSurfels(); 

 C4CTangentialCoverGeometry tcovergeom ;
 tcovergeom.init(ks, 0,1, tcover, itfcs );

 uint idx = 0;
 Mathutils::ModuloComputer msmc( tcover.nbMaximalSegments());
 Mathutils::ModuloComputer sc( surfaceSize );

 TriangleFunction lambda;
 
 // First, we look for the first index.
 uint idx_source = 0;

 C4CTangentialCover::SurfelMaximalSegments sms = tcover.beginSMS( 0 );

 for ( uint idx= 0; idx < tcover.nbSurfels(); idx++ ){
 double tangente = C4CTangentialCoverGeometry::angleByLambdaMS (tcover, tcovergeom, lambda, sms);
 vectResult.push_back( tangente );
 tcover.nextSMS(sms );
 }
 if ( ks != 0 ) delete ks; 
 return vectResult;
 
}









// Associating each point of the polygon to the nearest point of the
// contour starting from the last point association. The initial point
// association is given by the nearest point of the whole contour.

template<class T> 
vector<int> 
getAssociations(const vector<Vector2i> &contour, const vector<T> &polygon){
  
  vector<int> vectAssociations;
   
  //initial association:
  T ptPol0 = polygon.at(0);
  Vector2D v (contour.at(0).x()-ptPol0.x(), contour.at(0).y()-ptPol0.y());
  double minDist = v.x()*v.x()+v.y()*v.y();
  int indiceMin=0;
  for(int j=1; j<contour.size(); j++){
    Vector2D vDist (contour.at(j).x()-ptPol0.x(), contour.at(j).y()-ptPol0.y() );
    double norm= vDist.x()*vDist.x()+vDist.y()*vDist.y();
    if(norm<minDist){
      indiceMin=j;	
      minDist=norm;
    }
  }
  vectAssociations.push_back(indiceMin);
  
  // Association of each polygon point to the nearest point 
  Vector2i pt0= contour.at(indiceMin);
  int posLastAssociation=indiceMin;
  int nbPtContourDone =0;
  T ptPol;
  int posPol=1;
  while(posPol!=polygon.size()){
    ptPol=polygon.at(posPol);
    Vector2i p=contour.at((posLastAssociation+1)%contour.size());
    Vector2D v (p.x()-ptPol.x(), p.y()-ptPol.y());
    double minDist = v.x()*v.x()+v.y()*v.y();
    int indiceMin=(posLastAssociation+1)%contour.size();
    for(int i=0; i< contour.size()-nbPtContourDone; i++){
      int posContour=(posLastAssociation+i)%contour.size();
      Vector2i p = contour.at(posContour);
      Vector2D vDist (p.x()-ptPol.x(), p.y()-ptPol.y());
      double norm= vDist.x()*vDist.x()+vDist.y()*vDist.y();
      if(norm<minDist){
	indiceMin=posContour;
	minDist=norm;
      }
    } 
    
    vectAssociations.push_back(indiceMin);
    if(indiceMin>=posLastAssociation){
      nbPtContourDone += indiceMin-posLastAssociation; 
    }else{
      nbPtContourDone += indiceMin+contour.size()-posLastAssociation; 
    }
    posLastAssociation= indiceMin;
    posPol++;
  }
  
  return vectAssociations;
}



