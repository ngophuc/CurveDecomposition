///////////////////////////////////////////////////////////////////////////////
// Test the length variation of maximal segments on digital contour
///////////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <vector>
#include "ImaGene/base/Vector2i.h"
#include "ImaGene/base/Arguments.h"
#include "ImaGene/base/StandardArguments.h"
#include "ImaGene/mathutils/Mathutils.h"
#include "ImaGene/mathutils/Statistics.h"
#include "ImaGene/timetools/Clock.h"
#include "ImaGene/dgeometry2d/C4CTangentialCover.h"
#include "ImaGene/dgeometry2d/FreemanChain.h"
#include "ImaGene/dgeometry2d/C4CIteratorOnFreemanChain.h"
#include "ImaGene/digitalnD/C4CIteratorOnFreemanChainSurface.h"
#include "ImaGene/digitalnD/KnSpace.h"
#include "ImaGene/digitalnD/C4CIteratorOnSurface.h"
#include "ImaGene/helper/ShapeHelper.h"


using namespace std;
using namespace ImaGene;


static Arguments args;






struct FreemanAndIndex{
  FreemanChain *fc;
  vector<uint> fc2trans;
};



struct IndexedZone{
  uint begin;
  uint end;
};






static const int nbIterationSpikes = 20;
static const int RESOLUTION=1200;
static int agrandissementEPS=1;





/**
 * @param fc the freeman chain code of the (closed) contour.
 * 
 * @param c2trans the mapping from the original contour to the contour
 * [fc], which is a subsampled version.
 */
bool
gatherStatsMSForOneScale( Statistics & stats_scale,
			  const FreemanChain & fc,
			  uint resolution );

/**
 * @param fc (maybe updated) the freeman chain code of the (closed)
 * contour. The chain code may be translated.
 * 
 * @return a dyn. alloc. statistics object storing the length of
 * maximal segments for each surfel. Contains as many statistics
 * variables as the number of surfels of [fc].
 */
Statistics*
fillStatsMS( FreemanChain & fc );


IndexedZone
getIndexZone( const vector<Vector2D> & vectPoints, 
	      double minSlopes );

void 
transformFreemanChain( FreemanChain & fc, 
		       vector<uint> & c2trans , 
		       const FreemanChain &fcSrc,  
		       int samplingSize, int xIni, int yIni,
		       uint nbIterationSpikes = 10 );

Statistics*
getStatMSFromFreeman( const FreemanChain & fc, 
		      int idx_surfel, 
		      int resolution, int x0, int y0 );

Statistics*
getStatMSFromFreeman( FreemanChain & fc, 
		      int idx_surfel );

vector< vector <FreemanAndIndex> > 
getSampledContoursDecal( const FreemanChain & fc,
			 int resolutionMax);

vector<Vector2D> 
getVectorSlopes( const FreemanChain & fc, 
		 int samplingSizeMax, 
		 vector<IndexedZone> & vectResolution, 
		 int samplingSizeStartAnalyse, 
		 double penteMin );

Vector2D 
getRegLinear( const vector<Vector2D> & vectorPoints );





IndexedZone getIndexZone(const vector<Vector2D> &vectPoints, double minSlopes );


Statistics * getStatMSFromFreeman(const FreemanChain & fc, int idx_surfel, int resolution, int x0, int y0 );







void afficheLine(Vector2i point1, Vector2i point2,int color);
void afficheContourXFIGPixels(FreemanChain & fc, int color, int zoom, int x0, int y0, int prof);
void afficheContourXFIG(FreemanChain & fc, int color, int epaisseurTrait, bool filled, bool closeCnt, int zoom, int x0, int y0);
void affichePoints(Vector2i point, int color);
void afficheCroix(Vector2i point, int color, int largeur);
void affichePoint(Vector2i point, int color, double taillePoints, int profondeur);
void affichePixel(Vector2i pixel, int color, int colorFill, double taille, int profondeur);
void afficheImage(string nomImage, int largeur, int hauteur);
void afficheContourPointSelect(const FreemanChain & fc, vector<IndexedZone> vectResolution );
void  afficheContourPointSelect(const FreemanChain & fc, vector<Vector2D> vectPentes, double thresholdMin, double thresholdMax, int color, int colorFill );
void afficheContourPointSelectResolution(const FreemanChain & fc, vector<IndexedZone> vectResolution );
void extraitPointSelect(const FreemanChain & fc, vector<IndexedZone> vectResolution );


void afficheContourValide(const FreemanChain & fc, vector<IndexedZone> vectResolution, int resolution, bool included );




void affichePixelCenter(Vector2i pixel,  double taille, iostream &iostr);





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
  args.addOption( "-samplingSizeMax", "-samplingSizeMax <n> ", "3"  );

  
  args.addBooleanOption( "-enteteXFIG", "-enteteXFIG: add xfig entete" );
  args.addOption( "-afficheContourSRC", "-afficheContourSRC <color> <epaisseur>  ", "0", "1" );
  args.addBooleanOption( "-afficheContourSampled", "-afficheContourSampled " );
  args.addBooleanOption( "-afficheContourSampledSAnalyse", "-afficheContourSampledSAnalyse " );
  
  args.addBooleanOption( "-afficheRegLineaire", "-afficheRegLineaire <index>");  
  args.addBooleanOption( "-afficheBruit", "-afficheBruit");  
  args.addBooleanOption("-afficheContourEpaisseur", "-afficheContourEpaisseur " );
  args.addOption("-afficheImage", "-afficheImage <string nomImage> <largeur> <hauteur>: image format is GIF, JPEG","0", "0", "0");
  


  args.addBooleanOption("-extraitCarre", "-extraitCarre (format x1 y1 lx1 lx1)");  
  
  
  args.addOption("-samplingSizeStartAnalyse", "-samplingSizeStartAnalyse <size>: set the sampling size for the begining of the analyse slope of the segment length  (defaut=1) ", "1"  );
  args.addOption( "-afficheSlope", "-afficheSlope  <MIN|MEAN|MAX>: affiche pour chaque point la pente de l'evolution de la longueur des segments max (en logscale), 0=Min 1= Mean, 2= Max moyennÃ© sur une fenÃªtre de taille width ajustÃ© en fonction de la taille premier contour", "MEAN" );
  
  args.addOption("-agrandissementEPS", "-agrandissementEPS <coef> : agrandissement utilisé par eps (defaut=1)", "1");
  args.addOption("-affichePointIndex", "-affichePointIndex <index>: affiche le point d'indexe initiale <index> ", "0");
  args.addOption("-affichePointsIndex", "-affichePointsIndex <index1> <index2> <index3>: affiche les  points d'indexes initiales <index I> (Croix, rond moyen, rond Gros) ", "0", "1", "2");
  args.addOption( "-afficheStat", "-afficheStat <n>: affiche stat des tailles de segments max pour le surfel <n>.", "0" );
  args.addOption( "-afficheContourValide", "-afficheContourValide <n> : affiche les points valides (n est supérieur à a si [a,b] est le premier intervalle contenant des pentes < slopeMax (slopeMax =0 par defaut)) pour la résolution <n>", "1" );
  args.addOption( "-afficheContourValideIncluded", "-afficheContourValideIncluded <n>: affiche les points valides (n est inclus [a,b] cad le premier intervalle contenant des pentes < slopeMax ) pour la résolution <n> et selon l'angle theta max (slopeMax =0 par defaut)", "1" );
  args.addOption( "-setSlopeMax"," -setSlopeMax  <theta>: défini la pente maxi  (utilisé par afficheContourValide, afficheContourValideIncluded and afficheSlopes) ", "0" );
  
  args.addOption("-afficheZonesPlates", "-afficheZonesPlates <seuilMin><seuilMax>: affiche les zones plates dont la pente est inférieur à seuil Min  ", "-0.75", "-0.5"  );
  args.addOption("-afficheZonesCourbes", "-afficheZonesCourbes <seuilMin> <seuilMax>: affiche les zones plates dont la pente est supérieur à seuilMin  (seuil max defaut=0)", "-0.5", "0.0"  );
  args.addOption("-afficheCurve", "-afficheCurve <seuilMin> <seuilMax>: draw in color areas for which slope is sup as  ", "0.75","4", "4"   );
  args.addOption("-afficheContourSampledK", "-afficheContourSampledK <K> <color><profondeur XFIG> affiche le contour ramené à la résolution K  ", "1","50", "1"   );

  

  args.addOption("-afficheAssociation" , "-afficheAssociation <K> affiche les association entre le contour d'origine et le niveau K", "2"); 

  
  if ( ( argc <= 0 ) 
       || ! args.readArguments( argc, argv ) ) 
    {
      cerr << args.usage( "analysemslength", 
			  "Analyse for each pixel of the contour the length of the maximal segments. Exemple: testResolutions -samplingSizeMax 5 -afficheStat 0 -afficheContourSRC -enteteXFIG -afficheContourSampled  -samplingSizeStartAnalyse 1 -affichePointIndex 0 -afficheSlope  MEAN -afficheFlat -0.5 2  -afficheCurve -0.5  0.0  3  -afficheBruit   -agrandissementEPS 5 -afficheRegLineaire -extraitCarre  < formeTest8.fc > tmp.fig; fig2dev -Leps tmp.fig tmp.eps  ; gnuplot pltMSL2 ; gnuplot pltSlope "
			  ,"" ) << endl;
      return 1;
    }
  
  if(args.check("-agrandissementEPS")){
    agrandissementEPS = args.getOption( "-agrandissementEPS" )->getIntValue( 0 );
  }

  
  if(args.check("-enteteXFIG")){
    cout << "#FIG 3.2 \nLandscape \nCenter \nInches \nLetter  \n" <<agrandissementEPS<< "\nSingle \n1" 
	 << " \n"<<RESOLUTION<<" 1" << endl;
  }
  


  if(args.check("-afficheImage")){    
    string nomImage = args.getOption( "-afficheImage" )->getValue( 0 );
    int largeur = args.getOption( "-afficheImage" )->getIntValue( 1 );
    int hauteur = args.getOption( "-afficheImage" )->getIntValue( 2 );

    afficheImage(nomImage, largeur, hauteur);
  }
    
  


  vector<Vector2D> vectSlopes;
  vector<IndexedZone> vectResolution;





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

  
  uint samplingSizeMax = args.getOption( "-samplingSizeMax" )->getIntValue( 0 );
  
  //Vector2D point0(0.0,0.0);
  //  Vector2D 
  FreemanChain::const_iterator it = c.begin();
  Vector2i point0 ((*it).x(),(*it).y() );
  cerr << "first point = [" << point0.x() << ", "  << point0.y() << "]" << endl;
  
  Vector2i point0NK(0,0);
  
  

  
  vector<uint> c2trans ;
  FreemanChain fcNew;
  int samplingSizeStartAnalyse =1;
  double slopeMax = 0.0;
  
  
  transformFreemanChain(fcNew, c2trans,c, samplingSizeMax,0,0);  


  

  if(args.check("-afficheContourSampled")){  
    FreemanChain::const_iterator ite = fcNew.begin();
    for(int i=0;i <c2trans[0]; i++) {    
      ++ite;
    }         
    point0NK = Vector2i((*ite).x(),(*ite).y());
    ///Ramène dans les coordonnées de la résolution initiale.

    int decalX0 =  (point0.x()%samplingSizeMax) - (samplingSizeMax/2.0);
    int decalY0 =  (point0.y()%samplingSizeMax) - (samplingSizeMax/2.0);
    
    point0NK.x() = (point0NK.x())*samplingSizeMax+decalX0;
    point0NK.y() = (point0NK.y())*samplingSizeMax+decalY0;

    //    afficheContourXFIG(fcNew, 0, 50, false, true,samplingSizeMax, dx,dy);
    afficheContourXFIGPixels(fcNew, 30, samplingSizeMax,point0.x()-point0NK.x() ,point0.y()-point0NK.y(), 50);
  }  
  



  if(args.check("-afficheContourSampledK")){  
    
    vector<uint> c2transK ;
    FreemanChain fcNewK;
    int  k = args.getOption( "-afficheContourSampledK" )->getIntValue( 0 );
    int  color = args.getOption( "-afficheContourSampledK" )->getIntValue( 1 );
    int  profXFIG = args.getOption( "-afficheContourSampledK" )->getIntValue( 2 );
    transformFreemanChain(fcNewK, c2transK,c, k,0,0);  
    
    FreemanChain::const_iterator ite = fcNewK.begin();
    for(int i=0;i <c2transK[0]; i++) {    
      ++ite;
    }         
    point0NK = Vector2i((*ite).x(),(*ite).y());
    ///Ramène dans les coordonnées de la résolution initiale.

    int decalX0 =  (point0.x()%k) - (k/2.0);
    int decalY0 =  (point0.y()%k) - (k/2.0);
    
    point0NK.x() = (point0NK.x())*k+decalX0;
    point0NK.y() = (point0NK.y())*k+decalY0;

    //    afficheContourXFIG(fcNew, 0, 50, false, true,samplingSizeMax, dx,dy);
    afficheContourXFIGPixels(fcNewK, color, k,point0.x()-point0NK.x() ,point0.y()-point0NK.y(), profXFIG);
  }  


  if(args.check("-afficheAssociation")){      
    vector<uint> c2transK ;
    FreemanChain fcNewK;
    int  k = args.getOption( "-afficheAssociation" )->getIntValue( 0 );
    transformFreemanChain(fcNewK, c2transK,c, k,0,0);  
    
    
    FreemanChain::const_iterator ite = fcNewK.begin();
    for(int i=0;i <c2transK[0]; i++) {    
      ++ite;
    }
        
    point0NK = Vector2i((*ite).x(),(*ite).y());
    ///Ramène dans les coordonnées de la résolution initiale.

    int decalX0 =  (point0.x()%k) - (k/2.0);
    int decalY0 =  (point0.y()%k) - (k/2.0);
    
    point0NK.x() = (point0NK.x())*k+decalX0;
    point0NK.y() = (point0NK.y())*k+decalY0;

    int x0= point0.x()-point0NK.x();
    int y0= point0.y()-point0NK.y();
    
    
        
    FreemanChain::const_iterator iterC = c.begin();
    int posC=0;
    while(iterC!=c.end()){
      Vector2i pt = *iterC;
      FreemanChain::const_iterator iterFcNew = fcNewK.begin();
      for(int i=0;i <c2transK[posC]; i++) {    
	++iterFcNew;
      }
      Vector2i ptNew = *iterFcNew;
      ptNew.x()=ptNew.x()*k+x0;
      ptNew.y()=ptNew.y()*k+y0;
      //if(posC%2==0){
      afficheLine(pt, ptNew,4);
      affichePoint(pt,0,0.1,10);
      affichePoint(ptNew,1,0.1,10);
      //}
      ++iterC;
      posC++;
    }


    

//     //Vector2i pt ((*it).x()*zoom+x0, (*it).y()*zoom+y0);
//     //affichePixel(pt, 0,color, zoom, 59);

    
//     Vector2i pt (ptRk0.x()*resolution+x0, ptRk0.y()*resolution+y0);


      
  }  
  
  
  

  
  if(args.check("-afficheContourSRC")){  
    int color = args.getOption( "-afficheContourSRC" )->getIntValue( 0 );
    int epaisseur = args.getOption( "-afficheContourSRC" )->getIntValue( 1 );
    afficheContourXFIG(c,color, epaisseur, false, true, 1,0,0);
    //afficheContourXFIGPixels(c,0);
    
  }  


  
  if(args.check("-samplingSizeStartAnalyse")){
    samplingSizeStartAnalyse = args.getOption( "-samplingSizeStartAnalyse" )->getIntValue( 0 );
  }
  

  if(args.check("-affichePointsIndex")){

    int index1 = args.getOption( "-affichePointsIndex" )->getIntValue( 0 );
    int index2 = args.getOption( "-affichePointsIndex" )->getIntValue( 1 );
    int index3 = args.getOption( "-affichePointsIndex" )->getIntValue( 2 );
    
    FreemanChain::const_iterator ite = c.begin();
    for(int i=0;i <index1; i++) {    
      ++ite;
    }     
    afficheCroix(*ite,4, 4);    
    
    ite = c.begin();
    for(int i=0;i <index2; i++) {    
      ++ite;
    }     
    affichePoint(*ite,4,4, 20);    
    ite = c.begin();
    for(int i=0;i <index3; i++) {    
      ++ite;
    }     
    affichePoint(*ite,4, 8, 20);    
  }    
  
  




 
  if(args.check("-affichePointIndex")){
    int index = args.getOption( "-affichePointIndex" )->getIntValue( 0 );
    
    FreemanChain::const_iterator ite = c.begin();
    for(int i=0;i <index; i++) {    
      ++ite;
    }     
    afficheCroix(*ite,4, 4);    
    ite = fcNew.begin();    
    if(args.check("-afficheContourSampled")){
      uint newIndex = c2trans.at(index);    
      for(int i=0;i <newIndex; i++) {    
	++ite;
      }  
      int dx= point0.x() - point0NK.x();
      int dy= point0.y() - point0NK.y();

      Vector2i pRK((*ite).x()*samplingSizeMax+dx,(*ite).y()*samplingSizeMax+dy); 
      afficheCroix(pRK,1, 10);    
    }
  }    

  
  
  transformFreemanChain(fcNew, c2trans,c, samplingSizeStartAnalyse,0,0);  
  if(args.check("-afficheContourSampledSAnalyse")){  
    afficheContourXFIG(fcNew, 2, 20, false, true,1,0,0);
  } 

  
  if(args.check("-setSlopeMax")){  
    slopeMax = args.getOption( "-setSlopeMax" )->getFloatValue( 0 );
  }  
  


  if(args.check("-afficheStat")){
    int index = args.getOption( "-afficheStat" )->getIntValue( 0 );
    fstream msl, mslSel;
    
    msl.open("msl.txt",ios::out);
    mslSel.open("mslSel.txt", ios::out);
    
    
    fstream mslMoy;
    mslMoy.open ("mslMoy.txt", ios::out);
    Statistics statTotal(1, false);

    bool* statOK= new bool [samplingSizeMax]; 
    vector<Statistics*> all_stats;
    for( int k = samplingSizeStartAnalyse; k <= samplingSizeMax; k++ ) {
      all_stats.push_back( new Statistics( c.chain.size(), false ) );
      statOK[k-samplingSizeStartAnalyse] = gatherStatsMSForOneScale( *( all_stats.back() ), c, k );
      if(statOK[k-samplingSizeStartAnalyse]){
	mslMoy << k << " " << all_stats[k-samplingSizeStartAnalyse]->mean(index) << endl;
      }
    }
  



    
    //for(int k=samplingSizeStartAnalyse; k<= samplingSizeMax; k++){      
    //   mslMoy << k << " " << statDecalages.mean(0) << endl;       
    //}
    //cerr << "resolution associé à la taille des DSS max = " << statTotal.maxIndice(0)+1 << endl;
  }
  
  
  
  if(args.check("-afficheRegLineaire")){
    int index = args.getOption( "-afficheStat" )->getIntValue( 0 );    
    vector<Vector2D> vectPoints; 
    fstream streamReg;
    streamReg.open ("regLin.txt", ios::out);    
    
    for(int k=1; k<= samplingSizeMax; k++){           
      Statistics *stats = getStatMSFromFreeman(c, index, k, 0,0 );     
      if(stats->samples(0)!=0){
	vectPoints.push_back(Vector2D(log(k), log(stats->mean(0))));		
	streamReg << k << "  " << stats->mean(0) << endl;
      }                  
      
    }
    
    
    Vector2D regLin = getRegLinear(vectPoints );
    cerr << "slope = "<< regLin.x() << " b = " << exp(regLin.y()) << endl; 
  }


  
   if(args.check("-afficheSlope")){
     vectSlopes = getVectorSlopes(c, samplingSizeMax, vectResolution, samplingSizeStartAnalyse, slopeMax);    
     fstream streamSlope;
     streamSlope.open ("slopes.txt", ios::out);    
     vector<Vector2D>::iterator it= vectSlopes.begin();
     int i=0;
     while(it!=vectSlopes.end()){
       streamSlope << i << " " << (*it).x() << " " << (*it).y() << endl;
       ++it;
       ++i;
     }  
   }
   

   
   if(args.check("-afficheContourValide")){
     int resolution = args.getOption( "-afficheContourValide" )->getIntValue( 0 );    
     vectSlopes = getVectorSlopes(c, samplingSizeMax, vectResolution, samplingSizeStartAnalyse, slopeMax);    
     afficheContourValide(c, vectResolution,  resolution, false );     
   }

   if(args.check("-afficheContourValide")){
     int resolution = args.getOption( "-afficheContourValide" )->getIntValue( 0 );    
     vectSlopes = getVectorSlopes(c, samplingSizeMax, vectResolution, samplingSizeStartAnalyse, slopeMax);    
     afficheContourValide(c, vectResolution,  resolution, false );     
   }


   if(args.check("-afficheContourValideIncluded")){
     int resolution = args.getOption( "-afficheContourValideIncluded" )->getIntValue( 0 );    
     vectSlopes = getVectorSlopes(c, samplingSizeMax, vectResolution, samplingSizeStartAnalyse, slopeMax);    
     afficheContourValide(c, vectResolution,  resolution, true );     
   }

  
  if(args.check("-afficheBruit")){
    if(vectResolution.size()==0){
      vectSlopes =  getVectorSlopes(c, samplingSizeMax, vectResolution, samplingSizeStartAnalyse, slopeMax);    
    }
    afficheContourPointSelect(c, vectResolution);    
  }



  //###
  if(args.check("-extraitCarre")){
    if(vectResolution.size()==0){
      vectSlopes =  getVectorSlopes(c, samplingSizeMax, vectResolution, samplingSizeStartAnalyse, slopeMax);    
    }
    //affichePixelCenter(Vector2i pixel,  double taille);
    extraitPointSelect(c,  vectResolution );
  }


  
  if(args.check("-afficheContourEpaisseur")){
    if(vectResolution.size()==0){
      vectSlopes = getVectorSlopes(c, samplingSizeMax, vectResolution, samplingSizeStartAnalyse, slopeMax);    
    }
    //    afficheContourPointSelect(c, vectResolution);    
    afficheContourPointSelectResolution(c, vectResolution);
  }
  
  
  if(args.check("-afficheZonesPlates")){
    double seuilMin = args.getOption( "-afficheZonesPlates" )->getFloatValue( 0 );    
    double seuilMax = args.getOption( "-afficheZonesPlates" )->getFloatValue( 1 );    
     if(vectResolution.size()==0){
       vectSlopes = getVectorSlopes(c, samplingSizeMax, vectResolution, samplingSizeStartAnalyse, slopeMax);    
     }
     afficheContourPointSelect(c, vectSlopes, seuilMin, seuilMax, 1,18 );
     //afficheContourValide(c, vectResolution,  resolution, false );     
   }
  if(args.check("-afficheZonesCourbes")){
    double seuilMin = args.getOption( "-afficheZonesCourbes" )->getFloatValue( 0 );    
    double seuilMax = args.getOption( "-afficheZonesCourbes" )->getFloatValue( 1 );    
     if(vectResolution.size()==0){
       vectSlopes = getVectorSlopes(c, samplingSizeMax, vectResolution, samplingSizeStartAnalyse, slopeMax);    
     }
     afficheContourPointSelect(c, vectSlopes, seuilMin,seuilMax , 3,3  );
     //afficheContourValide(c, vectResolution,  resolution, false );     
   }



  

}




// Fin du main.
/////////////////////////////////////////////////////////////////////////////















vector<Vector2D>
getVectorSlopes(const FreemanChain & c, int samplingSizeMax, vector<IndexedZone> &vectResolution, 
		int samplingSizeStartAnalyse, double penteMin){  
  vector<Vector2D> vectSlopeAndOrigin;
  

  bool* statOK= new bool [samplingSizeMax]; 
  vector<Statistics*> all_stats;
  for( int k = samplingSizeStartAnalyse; k <= samplingSizeMax; k++ ) {
    all_stats.push_back( new Statistics( c.chain.size(), false ) );
    statOK[k-samplingSizeStartAnalyse] = gatherStatsMSForOneScale( *( all_stats.back() ), c, k );
  }
  
  
    
  for(int i=0; i<c.chain.size();i++){    
    vector<Vector2D> vectPoint;
    vector<Vector2D> vectPointAll;
    vector<double> vectSize;
    Vector2D ptPrev;
    bool foundFirstResolution = false;
    
    cerr << "[" << i << "]";
    
    int resolutionDebut; 
    int resolutionFin; 
    
    // Filling of the vectors containing all the points (vectPointsofPlot)
    // of the plot associated to lengths evolution acording the
    // resolution.
    
    vector<Vector2D> vectPointsofPlot;
    
    for(int k=samplingSizeStartAnalyse; k<= samplingSizeMax; k++){                 
      if(statOK[k-samplingSizeStartAnalyse]){
	Vector2D pt(log(k), log(all_stats[k-samplingSizeStartAnalyse]->mean(i)));
        vectPointsofPlot.push_back(pt);
	//cerr<< "taille Samples " << all_stats[k-samplingSizeStartAnalyse]->samples(i) << ", " << pt.y() << " ] " ;
      }else{
	cerr << "invalid contour " << endl; ;
      }
    }

    vector<Vector2D> vectPointsSelected;
    IndexedZone iz = getIndexZone(vectPointsofPlot, penteMin);
    int firstResolution;
    cerr << "[ "<< iz.begin << "," << iz.end << "] "; 
    if(iz.begin !=iz.end ){           
      for(int l=0; l < vectPointsofPlot.size(); l++){
	Vector2D pt = vectPointsofPlot.at(l);
	int position = round(exp(pt.x()));
	if(position >=iz.begin && position <=iz.end){
	  vectPointsSelected.push_back(vectPointsofPlot.at(l));
	}	
      }

    }
    
    cerr << " sel size[" << vectPointsSelected.size() << "] total sel [" << vectPointsofPlot.size() <<"]" ;
    //    iz.begin-=samplingSizeStartAnalyse;
    //iz.end-=samplingSizeStartAnalyse;
    
    vectResolution.push_back(iz);  
    Vector2D regLin = getRegLinear((vectPointsSelected.size()==0)?vectPointsofPlot:vectPointsSelected);
    vectSlopeAndOrigin.push_back(regLin);    
    cerr << endl;
  }   

  delete statOK;
  return vectSlopeAndOrigin;
  
}



/*
 * Return the index asociated to the first longest area with slope less or equal to @minSlopes 
 * If no slopes are found the indexedZone.begin and indexedZone.end are set to 0.  
 */

IndexedZone
getIndexZone(const vector<Vector2D> &vectPoints, double minSlopes ){
  int tailleVectP =  vectPoints.size();
  IndexedZone indexResult;
  indexResult.begin =1;
  indexResult.end =1;
  
  Vector2D previousPoint = vectPoints.at(0);
  Vector2D currentPoint;
  bool firstFound = false;
  for(int i=1; i< tailleVectP; i++){
    currentPoint = vectPoints.at(i); 
    double dy = currentPoint.y() - previousPoint.y();
    double dx = currentPoint.x() - previousPoint.x();
    if((dy/dx)<=minSlopes){	    
      if(!firstFound){
	indexResult.begin=round(exp(previousPoint.x()));
	indexResult.end=round(exp(currentPoint.x()));
	firstFound=true;
      }else{
	indexResult.end=round(exp(currentPoint.x()));
      }
    }else{
      if(firstFound){
	break;
      }
    }
    previousPoint = currentPoint;
  }
  return indexResult;
}

































Statistics * 
getStatMSFromFreeman(const FreemanChain & fc, int idx_surfel, int resolution, int x0, int y0 ) {
  FreemanChain fcNew;  
  vector<uint> c2trans;
  transformFreemanChain(fcNew, c2trans, fc, resolution, x0, y0);  
  return getStatMSFromFreeman(fcNew, c2trans.at(idx_surfel));
}







Statistics * 
getStatMSFromFreeman( FreemanChain & fc, int idx_surfel ) {
  int width =0;
  KnSpace* ks = ShapeHelper::makeSpaceFromFreemanChain( fc );
  
  C4CIteratorOnFreemanChain itfc;
  itfc.init( fc.begin(), true );
  C4CIteratorOnFreemanChainSurface itfcs;
  itfcs.init( ks, itfc );  
  C4CTangentialCover tcover;
  bool isInit = tcover.init( itfcs,0 );
  uint surfaceSize =   tcover.nbSurfels();
  
  uint idx = 0;
  Mathutils::ModuloComputer msmc( tcover.nbMaximalSegments());
  Mathutils::ModuloComputer sc( surfaceSize );
  Statistics* stats = new Statistics( 1, false );  

  //  cerr << "avant recalage: nb max segments " <<tcover.nbMaximalSegments()  << "nb  surf " << surfaceSize << "itfcs size: "
  //    << C4CIterator::size(itfcs) <<endl;
  if(tcover.nbMaximalSegments()<=1){
    stats->terminate();      
    return stats;
  }

  
  //recalage sur l'indice initiale
  C4CTangentialCover::SurfelMaximalSegments sms = tcover.beginSMS(0);
   while ( ( idx < sc.cast(idx_surfel) )
	   && ( tcover.nextSMS( sms ) ) ) ++idx;
      
   
   for ( uint idx_ms = sms.begin_ms; 
	 idx_ms != sms.end_ms; 
	 msmc.increment( idx_ms ) )
     {
       const C4CTangentialCover::MaximalSegment & ms =
	 tcover.getMaximalSegment( idx_ms );
       stats->addValue( 0, (double) ms.dss.size() - 1.0 );
     }

  stats->terminate();      
  return stats;
}











vector< vector <FreemanAndIndex> >
getSampledContoursDecal(const FreemanChain &fc, int resolutionMax){
  vector< vector<FreemanAndIndex> >  vectFaI;   
  
  for(int i=1; i<=resolutionMax; i++){
    vector<FreemanAndIndex> vectContourResoX;
    for(int x0=0; x0<i; x0++){  
      for(int y0=0; y0<i; y0++){  
	FreemanChain *fcNew = new FreemanChain();  
	vector<uint> c2trans;
	transformFreemanChain(*fcNew, c2trans, fc, i, x0, y0);
	FreemanAndIndex fai;
	fai.fc= fcNew;
	fai.fc2trans =c2trans; 
	vectContourResoX.push_back(fai);    
      }
    }
    
    vectFaI.push_back( vectContourResoX);
  }  
  return vectFaI;
}














//////////////////////////////////////////////////////////////////////////////////////
// Fonctions relatives à l'affichage
//////////////////////////////////////////////////////////////////////////////////////








void afficheContourXFIGPixels(FreemanChain & fc, int color, int zoom, int x0, int y0, int prof){
  for (FreemanChain::const_iterator it = fc.begin();it != fc.end(); ++it ){  
    Vector2i pt ((*it).x()*zoom+x0, (*it).y()*zoom+y0);
    affichePixel(pt, 0,color, zoom, prof);
  }
}



void afficheContourXFIG(FreemanChain & fc, int color, int epaisseurTrait, bool filled, bool closeCnt, int zoom, int x0, int y0){
  //param 6 :fill color
  //param 9: fill style
  //param 4: epaisseur
  //param 5: couleur trait
  cout <<setiosflags(ios_base::fixed);
  cout<< setprecision(0);
  
  cout << "\n 2 " << ((closeCnt)? 1 : 0 )<< " 0 "  << epaisseurTrait << " " << color  << 
    "  10 19 -1 "<<((filled)? 0:-1)  <<" 0.000 0 0 -1 0 0 " << fc.chain.size()+((closeCnt)? 1 : 0 ) << " "  <<endl ;
  for (FreemanChain::const_iterator it = fc.begin();it != fc.end(); ++it ){  
    cout << ((*it).x())*RESOLUTION*zoom+x0*RESOLUTION << " " << ((*it).y())*RESOLUTION*zoom+y0*RESOLUTION << " " ;        
  }
  if(closeCnt){
    cout << (((*(fc.begin())).x()))*RESOLUTION*zoom+x0*RESOLUTION << " " << ((*(fc.begin())).y())*RESOLUTION*zoom+y0*RESOLUTION << " " ;        
  }

}


void afficheCroix(Vector2i point, int color, int largeur){ 
  cout <<setiosflags(ios_base::fixed);
  cout<< setprecision(0);
  cout<< "\n 2 1 0 100 " <<color << " 7 10 -1 -1 0.000 0 0 -1 0 0 2"<<endl;
  cout << (point.x()-largeur)*RESOLUTION << " " << (point.y()-largeur)*RESOLUTION << " " 
       << (point.x()+largeur)*RESOLUTION << " " << (point.y()+largeur)*RESOLUTION << endl; 
  cout<< "\n 2 1 0 100 " <<color << " 7 10 -1 -1 0.000 0 0 -1 0 0 2"<<endl;
  cout << (point.x()+largeur)*RESOLUTION << " " << (point.y()-largeur)*RESOLUTION << " " 
       << (point.x()-largeur)*RESOLUTION << " " << (point.y()+largeur)*RESOLUTION << endl;   
}




void afficheLine(Vector2i point1, Vector2i point2,int color){ 
   cout <<setiosflags(ios_base::fixed);
   cout<< setprecision(0);
   cout<< "\n 2 1 0 1 " <<color << " 7 20 -1 -1 0.000 0 0 -1 0 0 2"<<endl;
   cout << (point1.x())*RESOLUTION << " " << (point1.y())*RESOLUTION << " " 
	<< (point2.x())*RESOLUTION << " " << (point2.y())*RESOLUTION << endl; 
}



void affichePoint(Vector2i point, int color, double taillePoints ,int profondeur){ 
  cout <<setiosflags(ios_base::fixed);
  cout<< setprecision(0);
  int taille = (int)(800.0*taillePoints);
  cout << "1 4 0 1 "<< color << " " << color <<" " << profondeur<< " -1 20 0.000 1 0.0000 "<< point.x()*RESOLUTION << " "<< point.y()*RESOLUTION 
       << " " << taille << " " << taille << " " << point.x()*RESOLUTION - 540 << " " << point.y()*RESOLUTION << " " << 
    point.x()*RESOLUTION << " "  <<point.y()*RESOLUTION << endl;
}





void affichePixel(Vector2i pixel, int color, int colorFill, double taille, int profondeur){
  Vector2D p0 (pixel.x()+taille/2.0, pixel.y()+taille/2.0);
  Vector2D p1 (pixel.x()+taille/2.0, pixel.y()-taille/2.0);
  Vector2D p2 (pixel.x()-taille/2.0, pixel.y()-taille/2.0);
  Vector2D p3 (pixel.x()-taille/2.0, pixel.y()+taille/2.0);
  
  cout << "\n 2 3 0 1 " << color  << " " << colorFill << 
    " " << profondeur<<    " -1 20 0.000 0 0 -1 0 0 5 " <<endl ;  
  
  cout << p0.x()*RESOLUTION << " " << p0.y()*RESOLUTION << " " << p1.x()*RESOLUTION << " " << p1.y()*RESOLUTION << " " << p2.x()*RESOLUTION << " " << p2.y()*RESOLUTION << " "<< 
    p3.x()*RESOLUTION << " " << p3.y()*RESOLUTION << " "  << p0.x()*RESOLUTION << " " << p0.y()*RESOLUTION << endl;
  
}






// affichage des points à partir du centre (pour David Coeurjolly)

void affichePixelCenter(Vector2i pixel,  double taille, iostream &iostr){
  iostr << pixel.x() << " " << pixel.y() << " " << taille << " " << taille << endl ;
}







void 
afficheContourPointSelect(const FreemanChain & fc, vector<Vector2D> vectPentes, double thresholdMin,
			  double thresholdMax, int color, int colorFill ){
  int i=0;
  for (FreemanChain::const_iterator it = fc.begin();it != fc.end(); ++it ){      
    Vector2i point = (*it);   
    double  pente = vectPentes.at(i).x();
    if(thresholdMin < pente && pente < thresholdMax){
      //  affichePoint(point, color, 1, 50);      
      affichePixel(point, colorFill, colorFill,4, 50);      
    }    
    i++;
  } 
  
}






void afficheContourPointSelect(const FreemanChain & fc, vector<IndexedZone> vectResolution ){
  int i=0;
  for (FreemanChain::const_iterator it = fc.begin();it != fc.end(); ++it ){      
    Vector2i point = (*it);
    int  resolution = vectResolution.at(i).begin;
    //    if(resolution!=1){
      affichePixel(point,1,11, resolution, 50);      
      // }
    i++;
  }   
}






void extraitPointSelect(const FreemanChain & fc, vector<IndexedZone> vectResolution ){
  int i=0;
  fstream fstr;
    
  fstr.open("square.txt",ios::out);
  


  for (FreemanChain::const_iterator it = fc.begin();it != fc.end(); ++it ){      
    Vector2i point = (*it);
    int  resolution = vectResolution.at(i).begin;
    //    if(resolution!=1){
    affichePixelCenter(point, resolution, fstr);      
    // }
    
    i++;
  }   
}




void afficheContourPointSelectResolution(const FreemanChain & fc, vector<IndexedZone> vectResolution ){
  int i=0;
  FreemanChain::const_iterator iii = fc.begin();
  Vector2i point0 ((*iii).x(),(*iii).y() );

  
  for (FreemanChain::const_iterator it = fc.begin();it != fc.end(); ++it ){      
    Vector2i point = (*it);
    int  resolution = vectResolution.at(i).begin;
    //    if(resolution!=1){

    vector<uint> c2tr ;
    FreemanChain fcN;
    transformFreemanChain(fcN, c2tr,fc, resolution,0,0);  
    
    FreemanChain::const_iterator itera = fcN.begin();
    for(int j=0;j <c2tr[i]; j++) {    
      ++itera;
    }         
    FreemanChain::const_iterator itera0 = fcN.begin();
    for(int j=0;j <c2tr[0]; j++) {    
      ++itera0;
    }         
    Vector2i pt0k =  *(itera0);
    
    
    Vector2i ptRk =  *itera;
    Vector2i ptRk0 =  *itera;
    int decalX0 =  (point0.x()%resolution) - (resolution/2.0);
    int decalY0 =  (point0.y()%resolution) - (resolution/2.0);
    
    ptRk.x() = (ptRk.x())*resolution+decalX0;
    ptRk.y() = (ptRk.y())*resolution+decalY0;

    int x0= point0.x()-pt0k.x()*resolution;
    int y0= point0.y()-pt0k.y()*resolution;
    //    afficheContourXFIG(fcNew, 0, 50, false, true,samplingSizeMax, dx,dy);
    //    afficheContourXFIGPixels(fcNew, 30, samplingSizeMax,point0.x()-point0NK.x() ,point0.y()-point0NK.y());
    

    //Vector2i pt ((*it).x()*zoom+x0, (*it).y()*zoom+y0);
    //affichePixel(pt, 0,color, zoom, 59);

    
    Vector2i pt (ptRk0.x()*resolution+x0, ptRk0.y()*resolution+y0);
    affichePixel(pt, 1,11, resolution, 59);
    
    i++;
  }   
}




void afficheContourValide(const FreemanChain & fc, vector<IndexedZone> vectResolution, int resolution, bool included ){
  int i=0;
  for (FreemanChain::const_iterator it = fc.begin();it != fc.end(); ++it ){      
    Vector2i point = (*it);
    int  resDeb = vectResolution.at(i).begin;
    int  resFin = vectResolution.at(i).end;
    if(resDeb<=resolution   && resDeb!=resFin && (included? (resFin>=resolution) :true)){  
      //affichePixel(point,11,11, 1, 50);      //affichePoint(point,4, 2, 50);      
    }else{
      affichePixel(point,4,27, resolution+2, 49);      //affichePoint(point,1, 1, 49);      
    }
    i++;
  }   
}





void  afficheImage(string nomImage, int largeur, int hauteur){
    cout << "2 5 0 1 0 -1 155 -1 -1 0.000 0 0 -1 0 0 5"<< endl<< " 0 " <<nomImage << endl << " 0 0 " << largeur*RESOLUTION  << " 0 " << largeur*RESOLUTION << " " << hauteur*RESOLUTION << " 0 " << hauteur*RESOLUTION << " 0 0 " << endl ; 
  

}








//---------------------------------------------------------------
//---------------------------------------------------------------
// Code reprise Jaco




/**
 * Given the freeman chain [fcSrc], computes a subsampled freeman
 * chain [fc] with multiscale factor (samplingSize,samplingSize) and
 * origin shift of (xIni,yIni). [c2trans] gives the index
 * transformation from the source contour to the destination contour.
 */
void
transformFreemanChain( FreemanChain & fc, 
		       vector<uint> & c2trans, 
		       const FreemanChain & fcSrc, 
		       int samplingSize,
		       int xIni, int yIni,
		       uint nbIterationSpikes ) {
  FreemanChain c2;
  vector<uint> c2subc;
  vector<uint> subc2c;
  int h = samplingSize;
  int v = samplingSize;
  FreemanChain::subsample( c2, c2subc, subc2c, fcSrc, 
			   h, v, 
			   xIni, yIni);
  
  int tailleChaine = fcSrc.chain.size();
  c2trans.clear();

  for(int k=0; k<tailleChaine; k++){
    c2trans.push_back(c2subc.at(k));
  }

  for(uint i=0; i<nbIterationSpikes; i++){
    FreemanChain c3;
    FreemanChain c4;
    vector<uint> o2i;
    vector<uint> i2o;     
    vector<uint> o2i2;
    vector<uint> i2o2;     
    if ( ! FreemanChain::cleanOuterSpikes( c3, o2i, i2o, c2, true  ) )
      cerr << "Contour with no interior !" << endl;
    if ( ! FreemanChain::cleanOuterSpikes( c4, o2i2, i2o2, c3, false  ) )
    cerr << "Contour with no interior !" << endl;
    for(int k=0; k<tailleChaine; k++){
      int newIndex = c2trans.at(k);      
      newIndex = o2i.at(newIndex%o2i.size());
      newIndex = o2i2.at(newIndex%o2i2.size());
      c2trans.at(k)=newIndex%tailleChaine;
    }
    
    c2.chain=c4.chain;    
  }
  fc.chain= c2.chain;  
}

/**
 * @param fc the freeman chain code of the (closed) contour.
 * 
 * @param c2trans the mapping from the original contour to the contour
 * [fc], which is a subsampled version.
 */
bool
gatherStatsMSForOneScale( Statistics & stats_scale,
			  const FreemanChain & fc,
			  uint resolution )
{
  int k = (int) resolution;
  uint src_size = fc.chain.size();
  if ( stats_scale.nb() != fc.chain.size() )
    cerr << "[gatherStatsMSForOneScale] statistics and freeman chain mismatch."
	 << endl;
  
  // Computes all possible shifts for more robust multiscale analysis.
  cerr << "+---- computing scale " << k << " " << flush;
  bool auMoins1 = false;
  for(int x0 = 0; x0 < k; x0++ ) {
    for(int y0 = 0; y0 < k; y0++ ) {	  
      FreemanChain fcNew;  
      vector<uint> c2trans;
      transformFreemanChain(fcNew, c2trans, fc, resolution, x0, y0);
      uint size = fcNew.chain.size();
      // Computes ms length statistics for one shift. 
      Statistics* stats1 = fillStatsMS( fcNew );
      // Relates these statistics to surfels on the original contour.
      if(stats1->samples(0)!=0){
	auMoins1=true;
	for ( uint i = 0; i < src_size; ++i ){
	  stats_scale.addValue( i, stats1->mean( c2trans[ i ] ) );	  
	}
      }
      delete stats1;
      cerr << "." << flush;
    }
  }
  cerr << " ended." << endl;
  stats_scale.terminate();
  return auMoins1;
}


/**
 * @param fc (maybe updated) the freeman chain code of the (closed)
 * contour. The chain code may be translated.
 * 
 * @return a dyn. alloc. statistics object storing the length of
 * maximal segments for each surfel. Contains as many statistics
 * variables as the number of surfels of [fc].
 */
Statistics*
fillStatsMS( FreemanChain & fc )
{
  int width =0;
  KnSpace* ks = ShapeHelper::makeSpaceFromFreemanChain( fc );
  
  // Computes maximal segments.
  C4CIteratorOnFreemanChain itfc;
  itfc.init( fc.begin(), true );
  C4CIteratorOnFreemanChainSurface itfcs;
  itfcs.init( ks, itfc );  
  C4CTangentialCover tcover;
  bool isInit = tcover.init( itfcs,0 );
  uint surfaceSize =   tcover.nbSurfels();
  Statistics* stats = new Statistics( surfaceSize, false );  
  

  uint idx = 0;
  Mathutils::ModuloComputer msmc( tcover.nbMaximalSegments());
  Mathutils::ModuloComputer sc( surfaceSize );
  // cerr << "[fillStatsMS] nbsurf=" << surfaceSize
  //      <<  " nbms=" << tcover.nbMaximalSegments() << endl;
  if(tcover.nbMaximalSegments()<=1){
    if ( ks != 0 ) delete ks;
    return stats;
  }
  
  // First, we look for the first index.
  uint idx_source = 0;
  C4CTangentialCover::SurfelMaximalSegments sms = 
    tcover.beginSMS( idx_source );
  do
    {
      // Compute statistics for this surfel.
      for ( uint idx_ms = sms.begin_ms; idx_ms != sms.end_ms; 
	    msmc.increment( idx_ms ) )
	{
	  const C4CTangentialCover::MaximalSegment & ms =
	    tcover.getMaximalSegment( idx_ms );
	  stats->addValue( sms.idx_surfel, (double) ms.dss.size() - 1.0 );
	}
      if ( ! tcover.nextSMS( sms ) ) break;
    }
  while ( sms.idx_surfel != idx_source );
  stats->terminate();      
  if ( ks != 0 ) delete ks;
  return stats;
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







