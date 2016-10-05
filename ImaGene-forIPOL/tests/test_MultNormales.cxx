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

struct FreemanAndIndex{
  FreemanChain *fc;
  vector<uint> fc2trans;
};



struct IndexedZone{
  uint begin;
  uint end;
};



struct ResolTangente{
  uint x0;
  uint y0;
  vector<double> tabTangente;
  bool isDefined;
};






/**
 * @param fc the freeman chain code of the (closed) contour.
 * 
 * @param c2trans the mapping from the original contour to the contour
 * [fc], which is a subsampled version.
 */
void
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




/**
 * @param fc (maybe updated) the freeman chain code of the (closed)
 * contour. The chain code may be translated.
 * 
 * @return a vector storing the normal
 * orientation for each surfel. 
 */

Statistics* 
fillStatsLambdaMST( FreemanChain & fc );



IndexedZone
getIndexZone( const vector<Vector2D> & vectPoints, 
	      double minSlopes );

void 
transformFreemanChain( FreemanChain & fc, 
		       vector<uint> & c2trans , 
		       const FreemanChain &fcSrc,  
		       int samplingSize, int xIni, int yIni,
		       uint nbIterationSpikes = 10 );
void 
transformFreemanChainTesting( FreemanChain & fc, 
			      vector<uint> & c2trans , 
			      const FreemanChain &fcSrc,  
			      int samplingSize, int xIni, int yIni,
			      uint nbIterationSpikes = 10 );
/**
 * Given the interpixel freeman chain [fcSrc], computes a subsampled
 * freeman chain [fc] with multiscale factor
 * (samplingSize,samplingSize) and origin shift of
 * (xIni,yIni). [c2trans] gives the index transformation from the
 * source contour to the destination contour.
 */
void
transformIPFreemanChainTesting( FreemanChain & fc, 
				vector<uint> & c2trans, 
				const FreemanChain & fcSrc, 
				int samplingSize,
				int xIni, int yIni,
				uint nbIterationSpikes );

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
		 int samplingSize, 
		 vector<IndexedZone> & vectResolution, 
		 int samplingSizeStartAnalyse, 
		 double penteMin );

Vector2D 
getRegLinear( const vector<Vector2D> & vectorPoints );




vector<uint> 
getNumContourMaximalLength(const FreemanChain &fc, uint resolution);


vector<ResolTangente> 
getAllTangentsLambdaMST( FreemanChain & fc, uint resolution );





vector<double>
getLambdaMSTFromScale(FreemanChain & fc, uint resolution, string mode );


vector<double> lambdaMSTEstimator (  FreemanChain & fc);



double  computeMeanAngle(const vector<double> &vectToMean);

vector<double>  lambdaMSTEstimator( const FreemanChain &fc,  FreemanChain &fcNew, const vector<uint> &fcNew2fc);



vector< vector<double> >  getLambdaMSTFromAllScales(FreemanChain & fc, uint scaleMax, string mode );




vector <uint>  getRegularizedNoiseLevel(const vector<uint> &vectNoiseLevel, uint tailleMasque);



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
  args.addOption( "-samplingSizeMax", "-samplingSizeMax <n>: choose how many scales are computed.", "3"  );
  args.addOption( "-nbIterationSpikes", "-nbIterationSpikes <n>: useful for very noisy contours, where their subsampling might be ill-formed. 10 is good for very noisy, 3 is enough for smooth.", "10"  );
  args.addOption( "-counting", "-counting <STD|CURVED|FLAT>: choose in which bins slopes are enumerated. STD: covers all slopes from -1.16 to 0.16. CURVED: covers mainly around -0.5:-0.333 with 5 bins. FLAT: covers mainly around -1:-0.75 with 5 bins.", "STD"  );
  args.addBooleanOption( "-global_stats", "-global_stats: displays mean, std dev, variance, unbiased variance, min and max of slopes."  );
  args.addOption( "-checkm", "-checkm <m> <x0> <y0> <T/N/I>: check computation of subsampled contour. T=Testing N=Normal I=IPTesting", "2", "0", "0", "T" );
  
  

  // Options pour test des normales:


  args.addOption("-extraitLambdaMST", "-extraitLambdaMST <resol> <TYPE> extrait les normales par Lambda MST à partir du niveau de resolution <resol> en fonction des différents contours obtenus par différents décallages (TYPE=<MEAN> moyenne des tangentes sur tous les contours, TYPE=<OPT> tangente du contour sur un point ayant une taille de segments max min maximal, TYPE=<SINGLE> un seul contour (le premier) ." , "1","MEAN" );
  

  args.addOption("-multiScalesLambdaMST", "-multiScalesLambdaMST <resol> <TYPE> extrait les normales par Lambda MST à partir du niveau de resolution <resol> en fonction des différents contours obtenus par différents décallages (TYPE=<MEAN> moyenne des tangentes sur tous les contours, TYPE=<OPT> tangente du contour sur un point ayant une taille de segments max min maximal, TYPE=<SINGLE> un seul contour (le premier) ." , "1","MEAN" );

  args.addOption("-setContourDecal", "-setContourDecal <x0> <y0> sélectionne le contour avec le décalage x0 et y0", "0", "0");
  args.addOption("-setRegul", "-setRegul <taille>  transforme la courbe des niveaux de bruit tel que la distance verticale entre deux points soit inf ou égal à taille ", "1");
  
  args.addOption("-drawNormals", "-drawNormals <freq> affiche le contour et les normales avec une frequence <freq>", "5");

  

  

  if ( ( argc <= 0 ) 
       || ! args.readArguments( argc, argv ) ) 
    {
      cerr << args.usage( "test_Multiscale", 
			  "Analyse some multiscale properties of a shape with respect to its expected asymptotic properties"
			  ,"" ) << endl;
      return 1;
    }

  // -------------------------------------------------------------------------
  // Read some arguments.
  uint samplingSizeMax = args.getOption( "-samplingSizeMax" )->getIntValue( 0 );
  uint nbIterationSpikes = args.getOption( "-nbIterationSpikes" )->getIntValue( 0 );
  int samplingSizeStartAnalyse =1;



  if(args.check("-setContourDecal")){
    x0Decal=args.getOption("-setContourDecal")->getIntValue(0);
    y0Decal=args.getOption("-setContourDecal")->getIntValue(1);
    
    
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







  // Test du calcul des tangentes 
  
  if(args.check("-extraitLambdaMST")){
    uint samplingSize = args.getOption("-extraitLambdaMST")->getIntValue(0);
    string type = args.getOption("-extraitLambdaMST")->getValue(1);
    
    vector<double> vectTangentes = getLambdaMSTFromScale(c, samplingSize, type); 
    
    if(args.check("-drawNormals")){

      DrawingXFIG::includeXFIGHeader(cout, 1024, 2);
      //DrawingXFIG::drawContourPixels(cout, c, 2, 1, 0, 0, 50);
      DrawingXFIG::drawContour(cout, c, 0, 1, 0, 0, 50);
	  vector<Vector2i> vectorPoints;
	  FreemanChain::getContourPoints(c,vectorPoints) ;
	  double tailleVector = 10.0;
	  uint mod = args.getOption("-drawNormals")->getIntValue(0);
	  for(int i=0; i < c.chain.size(); i++){
	    if(i%mod==0){
	      double tangente= vectTangentes.at(i);
	      Vector2i pt1 = vectorPoints.at(i);
	      Vector2i pt2 (pt1.x()+cos(tangente+M_PI/2.0)*tailleVector, pt1.y()+sin(tangente+M_PI/2.0)*tailleVector);
	      //DrawingXFIG::drawLine(cout, pt1, pt2, 2, 10, 40);
	      Vector2i pt3 (pt1.x()+cos(tangente-M_PI/2.0)*tailleVector, pt1.y()+sin(tangente-M_PI/2.0)*tailleVector);
	      DrawingXFIG::drawLine(cout, pt2, pt3, 2, 10, 40);
	    }
	  }
	  
    }else{
      for(int i=0; i < vectTangentes.size(); i++){
      cout << i << " " << vectTangentes.at(i) << endl;
    }
    }
    
  }



  // Test du calcul des tangentes 
  
  if(args.check("-multiScalesLambdaMST")){
    bool regul = args.check("-setRegul");
    
    
    uint samplingSizeMax = args.getOption("-multiScalesLambdaMST")->getIntValue(0);    
    FreemanChainSubsample fcsub( 1, 1, 0, 0 );
    FreemanChainCleanSpikesCCW fccs( nbIterationSpikes );
    FreemanChainCompose fcomp( fccs, fcsub );
    FreemanChainTransform* ptr_fct = &fcomp;
    FreemanChainSubsample* ptr_fcsub = &fcsub;
    
    MultiscaleProfile MP;
    MP.chooseSubsampler( *ptr_fct, *ptr_fcsub );
    MP.init( c, samplingSizeMax );

    string type = args.getOption("-multiScalesLambdaMST")->getValue(1);
    cerr << "computing all normals for all resolutions ..." ;
    vector <vector<double> > vectTangentes = getLambdaMSTFromAllScales(c, samplingSizeMax, type); 
    cerr << " [OK] " <<endl;
    
    vector<double> vectNoiseLevel ;
    for(int i=0; i < c.chain.size(); i++){
      uint noiseLevel = MP.noiseLevel(i,1, 0.0);
      vectNoiseLevel.push_back((double)noiseLevel);
    }
    
    vector<double> vectNoiseLevelRegul;
    if(regul){
      uint maxVariation  =  args.getOption("-setRegul")->getIntValue(0);
      cerr << "regullllllllll "<< maxVariation << endl;
      vectNoiseLevelRegul = CurveVariationsHelper::getCurveWithLowVariationsFromMax(vectNoiseLevel, maxVariation, true);
    }else
      vectNoiseLevelRegul= vectNoiseLevel;
    
    if(!args.check("-drawNormals")) {
      for(int i=0; i < c.chain.size(); i++){
	int noiseLevelRegul = ((int) vectNoiseLevelRegul.at(i))-1;
	int noiseLevel = ((int) vectNoiseLevel.at(i))-1;
	cout << i << " " << vectTangentes.at(noiseLevelRegul).at(i) << " " << noiseLevelRegul <<  " " << noiseLevel << endl;
      }
    }else{
      DrawingXFIG::includeXFIGHeader(cout, 1024, 2);
      //DrawingXFIG::drawContourPixels(cout, c, 2, 1, 0, 0, 50);
      DrawingXFIG::drawContour(cout, c, 0, 1, 0, 0, 50);
      vector<Vector2i> vectorPoints;
      FreemanChain::getContourPoints(c,vectorPoints) ;
      double tailleVector = 10.0;
      uint mod = args.getOption("-drawNormals")-> getIntValue(0);
      for(int i=0; i < c.chain.size(); i++){
	if(i%mod==0){
	  int noiseLevel = ((int)vectNoiseLevelRegul.at(i))-1;
	  double tangente = vectTangentes.at(noiseLevel).at(i);
	  Vector2i pt1 = vectorPoints.at(i);
	  Vector2i pt2 (pt1.x()+cos(tangente+M_PI/2.0)*tailleVector, pt1.y()+sin(tangente+M_PI/2.0)*tailleVector);
	  //DrawingXFIG::drawLine(cout, pt1, pt2, 2, 10, 40);
	  Vector2i pt3 (pt1.x()+cos(tangente-M_PI/2.0)*tailleVector, pt1.y()+sin(tangente-M_PI/2.0)*tailleVector);
	  DrawingXFIG::drawLine(cout, pt2, pt3, 2, 10, 40);
	  //  DrawingXFIG::drawLine(cout, pt3, pt1, 1, 10, 40);
	}
      }
      
    }
    

        
  }






  if ( args.check( "-checkm" ) )
    {
      FreemanChain subc;
      vector<uint> c2trans;
      int samplingSize = args.getOption( "-checkm" )->getIntValue( 0 );
      int xIni = args.getOption( "-checkm" )->getIntValue( 1 );
      int yIni = args.getOption( "-checkm" )->getIntValue( 2 );
      if ( args.getOption( "-checkm" )->getValue( 3 ) == "T" )
	transformFreemanChainTesting( subc, c2trans, c, 
				      samplingSize, xIni, yIni, 
				      nbIterationSpikes );
      else if ( args.getOption( "-checkm" )->getValue( 3 ) == "I" )
	transformIPFreemanChainTesting( subc, c2trans, c, 
					samplingSize, xIni, yIni, 
					nbIterationSpikes );
      else
	transformFreemanChain( subc, c2trans, c, 
			       samplingSize, xIni, yIni, 
			       nbIterationSpikes );
      cout << subc.x0 << " " << subc.y0 << " " << subc.chain << endl;
      return 0;
    }

  // -------------------------------------------------------------------------
  // Computes multiscale.
  uint src_size = c.chain.size();
  double slopeMax = 0.0;
  
  vector<Statistics*> all_stats;
  cerr << "+-- computing all scales from " << samplingSizeStartAnalyse
       << " to " << samplingSizeMax << "." << endl;
  for( int k = samplingSizeStartAnalyse; k <= samplingSizeMax; k++ ) {
    all_stats.push_back( new Statistics( src_size, false ) );
    gatherStatsMSForOneScale( *( all_stats.back() ), c, k );
  }
  
  // -------------------------------------------------------------------------
  // Verifies asymptotic properties.
  double delta = 1.0 / 24.0;
  double first = -1.0 - delta * 4.0;
  double last = 0.0 + delta * 4.0;
  if (args.getOption( "-counting" )->getValue( 0 ) == "CURVED" )
    {
      delta = 1.0 / 6.0;
      first = -1.0 / 2.0 - delta * 2.0;
      last = -1.0 / 3.0 + delta * 2.0;
    }
  else if (args.getOption( "-counting" )->getValue( 0 ) == "FLAT" )
    {
      delta = 1.0 / 4.0;
      first = -1.0 - delta * 1.0;
      last = -1.0 / 2.0 + delta * 1.0;
    }

  int n = (int) ( round( ( last - first ) / delta ) ) + 2;
  Statistics slopes( n, false );
  Statistics all_slopes( 1, false );
  for ( uint idx = 0; idx < src_size; ++idx )
    {
      vector<Vector2D> vectPoints; 
      for ( uint k = 0; k < all_stats.size(); ++k )
	{
	  if ( all_stats[ k ]->samples( idx ) !=0 ) 
	    vectPoints.push_back( Vector2D( log( k+1 ), 
					    log( all_stats[ k ]->mean( idx ) )
					    ) );
	}
      Vector2D regLin = getRegLinear(vectPoints );
      // cerr << "[" << idx << "] slope = "<< regLin.x() 
      // 	   << " b = " << exp(regLin.y()) << endl; 
      // Counts the surfel in the appropriate bin.
      double s = regLin.x();
      all_slopes.addValue( 0, s );
      uint bin = 0;
      if ( s >= last ) bin = n - 1;
      else if ( s >= first )
	bin = 1+(uint) ( floor( ( s - first ) / delta ) );
      slopes.addValue( bin, s );
    }
  for ( uint k = 0; k < all_stats.size(); ++k )
    delete all_stats[ k ];

  // Display number of elements in bins
  cout << "# test_Multiscale" << endl
       << "# Displays statistics relating multiscale to asymptotic properties."
       << endl;
  slopes.terminate();
  all_slopes.terminate();
  if ( args.check( "-counting" ) ) 
    {
      cerr << "# i nb percent moy interval_begin" << endl;
      for ( int i = 0; i < n; ++i ) 
	{
	  cout << i << " " << slopes.samples( i ) 
	       << " " << ( 100.0 * slopes.samples( i ) / src_size )
	       << " " << slopes.mean( i )
	       << " " << ( first + (i-1)*delta )
	       << endl;
	}
    }
  if ( args.check( "-global_stats" ) ) 
    {
      cerr << "# mean stddev var uvar min max" << endl;
      cerr << all_slopes.mean( 0 ) 
	   << " " << sqrt( all_slopes.variance( 0 ) )
	   << " " << all_slopes.variance( 0 )
	   << " " << all_slopes.unbiasedVariance( 0 )
	   << " " << all_slopes.min( 0 )
	   << " " << all_slopes.max( 0 ) 
	   << endl;
    }



   






  return 0;
}












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
    c2.x0 = c4.x0;
    c2.y0 = c4.y0;
    c2.chain=c4.chain;    
  }
  fc.x0 = c2.x0;
  fc.y0 = c2.y0;
  fc.chain= c2.chain;  
}



/**
 * Given the interpixel freeman chain [fcSrc], computes a subsampled
 * freeman chain [fc] with multiscale factor
 * (samplingSize,samplingSize) and origin shift of
 * (xIni,yIni). [c2trans] gives the index transformation from the
 * source contour to the destination contour.
 */
void
transformIPFreemanChainTesting( FreemanChain & fc, 
				vector<uint> & c2trans, 
				const FreemanChain & fcSrc, 
				int samplingSize,
				int xIni, int yIni,
				uint nbIterationSpikes ) 
{
  FreemanChainInnerCCW fcinner;
  FreemanChainOuterCCW fcouter;
  FreemanChainSubsample fcsub( samplingSize, samplingSize, xIni, yIni );
  FreemanChainCleanSpikesCCW fccs( nbIterationSpikes );
  FreemanChainCompose fcomp1( fcsub, fcinner );
  FreemanChainCompose fcomp2( fccs, fcomp1 );
  FreemanChainCompose fcomp3( fcouter, fcomp2 );
  FreemanChainCompose fcomp4( fccs, fcomp3 );
  vector<uint> trans2c;
  if ( ! fcomp4.apply( fc, c2trans, trans2c, fcSrc ) )
    cerr << "[transformIPFreemanChainTesting] Failure." << endl;

}

/**
 * Given the freeman chain [fcSrc], computes a subsampled freeman
 * chain [fc] with multiscale factor (samplingSize,samplingSize) and
 * origin shift of (xIni,yIni). [c2trans] gives the index
 * transformation from the source contour to the destination contour.
 */
void
transformFreemanChainTesting( FreemanChain & fc, 
			      vector<uint> & c2trans, 
			      const FreemanChain & fcSrc, 
			      int samplingSize,
			      int xIni, int yIni,
			      uint nbIterationSpikes ) 
{
  FreemanChainSubsample fcsub( samplingSize, samplingSize, xIni, yIni );
  FreemanChainCleanSpikesCCW fccs( nbIterationSpikes );
  FreemanChainCompose fcomp( fccs, fcsub );
  vector<uint> trans2c;
    if ( ! fcomp.apply( fc, c2trans, trans2c, fcSrc ) )
    cerr << "[transformFreemanChainTesting] Failure." << endl;


  // if ( fcSrc.chain.size() == 0 ) return;
  // // cout << "[transformFreemanChainTesting] 1" << endl;
  // FreemanChain c2;
  // vector<uint> c2subc;
  // vector<uint> subc2c;
  // int h = samplingSize;
  // int v = samplingSize;
  // FreemanChain::subsample( c2, c2subc, subc2c, fcSrc, 
  // 			   h, v, 
  // 			   xIni, yIni);
  // // cout << "[transformFreemanChainTesting] 2 " 
  // //      << " src_c=" << fcSrc.chain.size()
  // //      << " c2subc=" << c2subc.size()
  // //      << endl;
  
  // int tailleChaine = fcSrc.chain.size();
  // c2trans.clear();

  // for(int k=0; k<tailleChaine; k++){
  //   c2trans.push_back(c2subc.at(k));
  // }

  // for(uint i=0; i<nbIterationSpikes; i++){
  //   // cout << "[transformFreemanChainTesting] 3" << endl;
  //   FreemanChain c3;
  //   FreemanChain c4;
  //   vector<uint> o2i;
  //   vector<uint> i2o;     
  //   vector<uint> o2i2;
  //   vector<uint> i2o2;
  //   FreemanChain::cleanOuterSpikes( c3, o2i, i2o, c2, true  );
  //   FreemanChain::cleanOuterSpikes( c4, o2i2, i2o2, c3, false  );
  //   // if ( ! FreemanChain::cleanOuterSpikes( c3, o2i, i2o, c2, true  ) )
  //   //   {
  //   // 	cerr << "[FreemanChain::cleanOuterSpikes( c3, o2i, i2o, c2, true  )]"
  //   // 	     << " Contour with no interior !" << endl
  //   // 	     << " c_in =" << c2.chain << endl
  //   // 	     << " c_out=" << c3.chain << endl;
  //   //   }
  //   // if ( ! FreemanChain::cleanOuterSpikes( c4, o2i2, i2o2, c3, false  ) )
  //   //   {
  //   // 	cerr << "[FreemanChain::cleanOuterSpikes( c3, o2i, i2o, c2, false  )]"
  //   // 	     << " Contour with no interior !" << endl
  //   // 	     << " c_in =" << c3.chain << endl
  //   // 	     << " c_out=" << c4.chain << endl;
  //   //   }
  //   for ( int k = 0; k < tailleChaine; k++ ) {
  //     int newIndex = c2trans.at( k );      
  //     // JOL : take into account empty contours.
  //     newIndex = o2i.size() != 0 ? o2i.at( newIndex % o2i.size() ) : 0;
  //     newIndex = o2i2.size() != 0 ? o2i2.at( newIndex % o2i2.size() ) : 0;
  //     // newIndex = o2i2.at( newIndex % o2i2.size() );
  //     c2trans.at( k ) = newIndex % tailleChaine;
  //   }
  //   c2.x0 = c4.x0;
  //   c2.y0 = c4.y0;
  //   c2.chain=c4.chain;    
  // }
  // fc.x0 = c2.x0;
  // fc.y0 = c2.y0;
  // fc.chain= c2.chain;  
}

/**
 * @param fc the freeman chain code of the (closed) contour.
 * 
 * @param c2trans the mapping from the original contour to the contour
 * [fc], which is a subsampled version.
 */
void
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
  
  for(int x0 = 0; x0 < k; x0++ ) {
    for(int y0 = 0; y0 < k; y0++ ) {	  
      FreemanChain fcNew;  
      vector<uint> c2trans;
      //%%%%
      transformFreemanChainTesting( fcNew, c2trans ,fc ,resolution, x0, y0,10 );
      //transformFreemanChain(fcNew, c2trans, fc, resolution, x0, y0);
      uint size = fcNew.chain.size();
      // Computes ms length statistics for one shift. 
      Statistics* stats1 = fillStatsMS( fcNew );
      // Relates these statistics to surfels on the original contour.
      for ( uint i = 0; i < src_size; ++i ){
	stats_scale.addValue( i ,stats1->min( c2trans[ i ] ) );
      }
      delete stats1;
      
      cerr << "." << flush;
    }
  }
  cerr << endl;
  stats_scale.terminate();
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













/**
 * Pour une résolution donnée, détermine pour chaque pixel le numéro
 * du contour présentant en ce point le plus grand segment max en
 * moyenne.
 *
 **/


vector<uint> 
getNumContourMaximalLength(const FreemanChain &fc, uint resolution){
  vector<uint> vectNumContour;  
  Statistics *statScale =  new Statistics(fc.chain.size(), true);

  gatherStatsMSForOneScale(*statScale, fc, resolution);

  for(int i=0; i< fc.chain.size(); i++){
    vectNumContour.push_back(statScale->maxIndice(i));    
    //vectNumContour.push_back(0);    
  }
  delete statScale;
  
  return vectNumContour;
}










/**
 * @param fc (maybe updated) the freeman chain code of the (closed)
 * contour. The chain code may be translated.
 * 
 * @return a vector storing the normal
 * orientation for each surfel. 
 */


vector<double>
getLambdaMSTFromScale(FreemanChain & fc, uint resolution, string mode ){
  

  // Recupération du numéro du contour optimal pour chaque point du contour
  vector<uint> vectNumContourForSurfel = getNumContourMaximalLength(fc, resolution);
  
  uint tailleContour = fc.chain.size();
  
  vector<double> vectLambdaMST;

  
  vector<ResolTangente> resTgt = getAllTangentsLambdaMST(fc, resolution);

  
  for(int i=0; i< tailleContour; i++){
    
    uint numContour =  vectNumContourForSurfel.at(i);

    double tangente;
    if(mode=="OPT")
      tangente = resTgt.at(numContour).tabTangente.at(i);
    else if (mode == "MEAN"){
      vector<double> vectToMean ;
      for(int k =0 ; k< resolution*resolution; k++){
	if(resTgt.at(k).isDefined){
	  vectToMean.push_back(resTgt.at(k).tabTangente.at(i));
	}
      }
      tangente = computeMeanAngle(vectToMean);
    }
    else if (mode=="SINGLE"){
      tangente = resTgt.at(x0Decal*resolution+y0Decal).tabTangente.at(i);
    }else if (mode =="MEDIAN"){
      vector<double> vectToMed ;
      for(int k =0 ; k< resolution*resolution; k++){
	if(resTgt.at(k).isDefined){
	  vectToMed.push_back(resTgt.at(k).tabTangente.at(i)+2.0*M_PI);
	}
      }
      nth_element (vectToMed.begin(), vectToMed.begin()+(vectToMed.size()/2), vectToMed.end());
      
      tangente = Mathutils::AngleComputer::cast((float)*(vectToMed.begin()+(vectToMed.size()/2)));
      
    }
      
    

    vectLambdaMST.push_back(tangente);        
  }
  

  return vectLambdaMST;
}






/**
 *  Calcul les normales pour toutes les résolutions des contours
 *   
 *
 **/




vector< vector<double> > 
getLambdaMSTFromAllScales(FreemanChain & fc, uint scaleMax, string mode ){
  
  vector< vector <double> > vectResult;
  for(int i=1; i< scaleMax; i++){
    cerr << "Computing nomals scale " << i << endl;
    vectResult.push_back( getLambdaMSTFromScale(fc, i, mode)); 
    cerr << "Computing nomals scale " << i << " OK" <<endl;
  }  
  return vectResult;

}








/**
 * @param fc (maybe updated) the freeman chain code of the (closed)
 * contour. The chain code may be translated.
 * 
 * @return a vector storing the normal
 * orientation for each surfel. 
 */


vector<ResolTangente> 
getAllTangentsLambdaMST( FreemanChain & fc, uint resolution ){
  int tailleChaine = fc.chain.size();
  vector<ResolTangente> vectResult ;
  
  for(int x0 = 0; x0 < resolution; x0++ ) {
    for(int y0 = 0; y0 < resolution; y0++ ) {	  
      FreemanChain fcNew;  
      vector<uint> c2transf, trans2cf;
      transformFreemanChainTesting( fcNew, c2transf ,fc ,resolution, x0, y0,10 );
      //FreemanChain::subsample(fcNew, c2transf, trans2cf, fc, resolution, resolution, x0, y0);
      //transformFreemanChain(fcNew, c2trans, fc, resolution, x0, y0);
      
      ResolTangente resolTgt;
      resolTgt.x0 = x0;
      resolTgt.y0 = y0;
      resolTgt.tabTangente = lambdaMSTEstimator(fc, fcNew, c2transf);
      resolTgt.isDefined = (resolTgt.tabTangente.size()!=0);
      if(!resolTgt.isDefined)
	cerr << "isNotDefined......" << endl;
      vectResult.push_back(resolTgt); 
      
    }
  }
  
  return vectResult;
}





double 
computeMeanAngle(const vector<double> &vectToMean){
  double angleResult = 0.0;
  double sizeVectToMean =  (double)vectToMean.size();
  for(uint i = 0; i< vectToMean.size();i++){
    angleResult+=Mathutils::AngleComputer::deviation(vectToMean.at(i),vectToMean.at(0));
  }
  
  return Mathutils::AngleComputer::cast( angleResult/sizeVectToMean + vectToMean.at(0));
}






vector<double> 
lambdaMSTEstimator(const  FreemanChain &fc,  FreemanChain & fcNew, const vector<uint> &fc2fcNew){
  uint tailleSRC = fc.chain.size();
  vector<double> vectPosition ; 
  
  for(int i=0;i < tailleSRC; i++){
    vectPosition.push_back(0.5);
  }
  
  //vector<double> tabTangent  = lambdaMSTEstimator(fcNew, vectPosition);
  vector<double> tabTangent  = lambdaMSTEstimator(fcNew);
  vector<double> tabResult;  
  for(int i=0; i< tailleSRC; i++){
    tabResult.push_back(tabTangent.at(fc2fcNew.at(i)));
  }
  return tabResult;
}












/**
 * @param fc (maybe updated) the freeman chain code of the (closed)
 * contour. The chain code may be translated.
 * 
 * @return a vector storing the normal
 * orientation for each surfel. 
 */


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
  uint surfaceSize =   tcover.nbSurfels();
  

  C4CTangentialCoverGeometry tcovergeom ;
  tcovergeom.init(ks, 0,1, tcover, itfcs );  

  uint idx = 0;
  Mathutils::ModuloComputer msmc( tcover.nbMaximalSegments());
  Mathutils::ModuloComputer sc( surfaceSize );

  if(tcover.nbMaximalSegments()<=1){
    cerr << "pb nb max segments <1!!!!" << endl;
    if ( ks != 0 ) delete ks;
    return vectResult;
  }
  
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











// Essaie pour la régulation des niveaux de bruit


vector <uint> 
getRegularizedNoiseLevel(const vector<uint> &vectNoiseLevel, uint tailleMasque){
  
  vector<uint> vectResult;
  uint taille = vectNoiseLevel.size();
  for(uint i=0; i<taille; i++){
    vector<uint> vToMean; 
    
    vToMean.push_back(vectNoiseLevel.at(i));
    for(uint k=1; k<= tailleMasque; k++){
      vToMean.push_back(vectNoiseLevel.at((i-k+taille)%taille));
      vToMean.push_back(vectNoiseLevel.at((i+k+taille)%taille));
    }
    nth_element (vToMean.begin(), vToMean.end()-1, vToMean.end());
    uint meanVal =  *(vToMean.end()-1);
    vectResult.push_back(meanVal);    
  }
  cerr << "comp mean ok" << endl;
  
  return vectResult;
}


