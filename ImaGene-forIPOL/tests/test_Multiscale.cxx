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
#include "ImaGene/dgeometry2d/FreemanChainTransform.h"
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
  cerr << "+---- computing scale " << k << " " << flush;
  for(int x0 = 0; x0 < k; x0++ ) {
    for(int y0 = 0; y0 < k; y0++ ) {	  
      FreemanChain fcNew;  
      vector<uint> c2trans;
      transformFreemanChain(fcNew, c2trans, fc, resolution, x0, y0);
      uint size = fcNew.chain.size();
      // Computes ms length statistics for one shift. 
      Statistics* stats1 = fillStatsMS( fcNew );
      // Relates these statistics to surfels on the original contour.
      for ( uint i = 0; i < src_size; ++i )
	stats_scale.addValue( i, stats1->mean( c2trans[ i ] ) );

      delete stats1;
      cerr << "." << flush;
    }
  }
  cerr << " ended." << endl;
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

