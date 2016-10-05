///////////////////////////////////////////////////////////////////////////////
// Extracts the noise level and meaningful scale of a digital contour.
///////////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include "ImaGene/helper/DrawingXFIG.h"
#include "ImaGene/base/Vector2i.h"
#include "ImaGene/base/Arguments.h"
#include "ImaGene/base/StandardArguments.h"
#include "ImaGene/mathutils/G.h"
#include "ImaGene/mathutils/Mathutils.h"
#include "ImaGene/mathutils/Statistics.h"
#include "ImaGene/timetools/Clock.h"
#include "ImaGene/dgeometry2d/C4CTangentialCover.h"
#include "ImaGene/dgeometry2d/FreemanChain.h"
#include "ImaGene/dgeometry2d/C4CIteratorOnFreemanChain.h"
#include "ImaGene/dgeometry2d/C4CGeometry.h"
#include "ImaGene/dgeometry2d/C4CSegmentPencil.h"
#include "ImaGene/dgeometry2d/MultiscaleFreemanChain.h"
#include "ImaGene/digitalnD/C4CIteratorOnFreemanChainSurface.h"
#include "ImaGene/digitalnD/KnSpace.h"
#include "ImaGene/digitalnD/K2Space.h"
#include "ImaGene/digitalnD/KnShapes.h"
#include "ImaGene/digitalnD/Frame2D.h"
#include "ImaGene/digitalnD/C4CIteratorOnSurface.h"
#include "ImaGene/helper/ShapeHelper.h"
#include "ImaGene/helper/MultiscaleProfile.h"


using namespace std;
using namespace ImaGene;


static Arguments args;

Vector2D 
getRegLinear( const vector<Vector2D> & vectorPoints );

void
testLinearRegression( double p );

void
convolveByGaussianKernels( vector<double> & sx,
			   vector<double> & sy,
			   const vector<double> & x,
			   const vector<double> & y,
			   const vector<uint> & nl,
			   const vector<G> & g_fcts );




class LambdaMSTEstimator
{
  static const uint m_nb_ms = 20;
  KnSpace* m_ks;
  mutable Frame2D m_frame;
  TriangleFunction m_l;
  DTriangleFunction m_lp;
  C4CSegment* m_segments;
  uint m_x;
  uint m_y;

public:
  LambdaMSTEstimator()
    : m_ks( 0 ), m_segments( 0 )
  {
    m_segments = new C4CSegment[ m_nb_ms ];
  }
  LambdaMSTEstimator( const LambdaMSTEstimator & other )
    : m_ks( other.m_ks ), m_segments( 0 )
  {
    m_segments = new C4CSegment[ m_nb_ms ];
    init( *m_ks, other.m_x, other.m_y );
  }

  LambdaMSTEstimator &
  operator=( const LambdaMSTEstimator & other )
  {
    if ( this != &other )
      {
	m_ks = other.m_ks;
	init( *m_ks, other.m_x, other.m_y );
      }
    return *this;
  }

  ~LambdaMSTEstimator()
  {
    if ( m_segments != 0 ) delete[] m_segments;
  }

  const KnSpace* space() const
  { 
    return m_ks;
  }


  void init( KnSpace & ks, uint x = 0, uint y = 1 )
  {
    m_ks = &ks;
    m_frame.init( m_ks, x, y );
    m_x = x;
    m_y = y;
  }

  uint otherDir( uint d ) const
  {
    return ( d == m_x ) ? m_y : m_x;
  }

  /**
   * Tangent estimation at surfel position [shift_it], the [shift]
   * specifies the position within.
   *
   * @param it any iterator.
   * @param shift a value between 0.0 and 1.0
   */
  double estimation( const C4CIteratorOnSurface & input_it, 
		     double shift = 0.5 ) const
  {
    // Builds the pencil of maximal segments.
    uint j = 0;
    uint k;
    if ( C4CGeometry::maximalSegments( input_it, m_segments, j, k, m_nb_ms ) )
      {
	// All geometric computations are made in the local frame of
	// the current boundary element.
	C4CSegmentPencil pencil( m_segments, j, k, m_nb_ms, m_l, m_lp );
	double theta = pencil.angleToX( Vector2D( shift, 0.0 ) );
	// Cast angle in the global frame.
	Kn_sid surfel = input_it.current();
	// cerr << "surfel=";
	// m_ks->displayKn_sid( surfel, std::cerr ); 
	// cerr << endl;
	m_frame.setSurfelFrame( surfel, 
				otherDir( m_ks->sorthDir( surfel ) ) );
	return m_frame.angleToX( theta );
      }
    return -1.0;
  }

};  

void 
createLambdaMSTEstimators
( map<MultiscaleFreemanChain::SubsampledChainKey, LambdaMSTEstimator>
  & map_estimators,
  MultiscaleProfile & MP,
  uint samplingSizeMax )
{
  for ( uint k = 1; k <= samplingSizeMax; ++k )
    for ( int x0 = 0; x0 < k; ++x0 )
      for ( int y0 = 0; y0 < k; ++y0 )
	{
	  MultiscaleFreemanChain::SubsampledChainKey key( k, k, x0, y0 );
	  MultiscaleFreemanChain::SubsampledChain* subchain = MP.get( key );
	  if ( subchain != 0 )
	    {
	      LambdaMSTEstimator LMST;
	      // FreemanChain copyc;
	      // copyc.chain = subchain->subc.chain;
	      // copyc.x0 = subchain->subc.x0;
	      // copyc.y0 = subchain->subc.y0;
	      KnSpace* ks 
		= ShapeHelper::makeSpaceFromFreemanChain( subchain->subc );
	      LMST.init( *ks );
	      map_estimators[ key ] = LMST;
	      // cerr << "key=(" << k << "," << k 
	      // 	   << "," << x0 << "," << y0 << ") LMST=[ks=("
	      // 	   << ks->size( 0 ) << "," << ks->size( 1 ) << ") chain="
	      // 	   << subchain->subc.chain << "]" << endl;
	    }
	}
}

double averageEstimation
( map<MultiscaleFreemanChain::SubsampledChainKey, LambdaMSTEstimator>
  & map_estimators,
  const MultiscaleProfile & MP,
  uint idx, uint res )
{
  Mathutils::AngleComputer ac;
  C4CIteratorOnFreemanChain itfc;
  C4CIteratorOnFreemanChainSurface itfcs;
  double t0 = -1.0;
  double average = 0.0;
  uint nb = 0;
  for ( int x0 = 0; x0 < res; ++x0 )
    for ( int y0 = 0; y0 < res; ++y0 )
      {
	MultiscaleFreemanChain::SubsampledChainKey key( res, res, x0, y0 );
	const MultiscaleFreemanChain::SubsampledChain* subchain
	  = MP.get( key );
	if ( subchain != 0 )
	  {
	    const LambdaMSTEstimator & LMLP
	      = map_estimators[ key ];
	    pair<uint,double> surfel = subchain->subsurfel( idx );
	    uint surf_idx = surfel.first;
	    double shift = surfel.second;
	    FreemanChain::const_iterator p ( subchain->subc,
					     surf_idx );
	    itfc.init( p, true );
	    itfcs.init( LMLP.space(), itfc );
	    double t = LMLP.estimation( itfcs, shift );
	    ++nb;
	    if ( t0 < 0.0 ) 
	      t0 = t;
	    else
	      {
		average += ac.deviation( t, t0 );
	      }
	  }
      }
  average = ac.cast( ( average / nb ) + t0 );
  return average;
}

void 
multiresMedianPGM( const MultiscaleProfile & MP,
		   uint mscales_min_size, 
		   double mscales_max_slope );

void
multiresMedianPGMOnOriginalImage
( const MultiscaleProfile & MP,
  const FreemanChain & fc,
  uint mscales_min_size, 
  double mscales_max_slope,
  KnSpace* ks,
  KnRUCellVector<int> & vectorVal);

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
  args.addOption( "-samplingSizeMax", "-samplingSizeMax <n>: choose how many scales are computed.", "10"  );
  args.addOption( "-nbIterationSpikes", "-nbIterationSpikes <n>: useful for very noisy contours, where their subsampling might be ill-formed. 10 is good for very noisy, 3 is enough for smooth.", "10"  );
  args.addOption( "-profile", "-profile <idx>: draws the profile of the surfel with index <idx> in the input freeman chain.", "0"  );
  args.addOption( "-detailedProfile", "-detailedProfile <idx>: draws the detailed profile of the surfel with index <idx> in the input freeman chain.", "0"  );
  args.addOption( "-affMapping", "-affMapping <r> <x0> <y0>: outputs the mapping between finest resolution and the given subsampling.", "1", "0", "0" );
  args.addOption( "-affTgt", "-affTgt <r> <x0> <y0>: outputs the tangents estimated at the given subsampling.", "1", "0", "0" );
  args.addOption( "-meaningfulScales", "-meaningfulScales <min_size> <max_slope>: specifies parameters for defining meaningful scales: minimum size of the interval of scales and maximum slopes between consecutivesamples within.", "1", "-0.2" );
  args.addBooleanOption( "-affNoiseLevel", "-affNoiseLevel: gives for each surfel its noise level (as specified by -meaningfulScales)." );
  args.addOption( "-multiresTgt", "-multiresTgt <SINGLE|MEAN|MAX>: computes the tangent direction according to the noise level (as specified by -meaningfulScales). For each surfel, SINGLE means selecting the subsampling at (0,0), MEAN means selecting the subsampling with longest mean of maximal segments, MAX means selecting the subsampling with the longest maximal segment.", "MEAN" );
  args.addBooleanOption( "-averageMultiresTgt", "-averageMultiresTgt: computes the tangent direction according to the noise level (as specified by -meaningfulScales). For each surfel, all results at a given scale are averaged." );
  args.addOption( "-averageTgt", "-averageTgt <scale>: computes the tangent direction. For each surfel, all results at the given scale are averaged.", "1" );
  args.addBooleanOption( "-smoothContour", "-smoothContour: computes a smooth version of the input digital contour according to the noise level (as specified by -meaningfulScales). NB: not a very good idea since noise level is related to tangent estimation and not position estimation." );
  args.addOption( "-gsmooth", "-gsmooth <sigma>: computes a smooth version of the input digital contour by convolution with a Gaussian kernel of deviation <sigma>.", "1.0" );
  args.addBooleanOption( "-multiresMedianPGM", "-multiresMedianPGM: computes the median filtering of the PGM image directed by the noise level of the contour." );
  args.addOption( "-multiresMedianPGMOnOriginalImage", "-multiresMedianPGMOnOriginalImage <filename.pgm>: computes the median filtering of the specified PGM image directed by the noise level of the contour.", "toto.pgm" );
 
  // Rajout BK ExpePAMI

  args.addOption("-setCstNoise","-setCstNoise <level> used to set a constant noise level used by -multiresMedianPGM", "1" );
  

  if ( ( argc <= 0 ) 
       || ! args.readArguments( argc, argv ) ) 
    {
      cerr << args.usage( "scale_analysis", 
			  "Analyse some multiscale properties of a shape given as a Freemanc chain with respect to its expected asymptotic properties"
			  ,"" ) << endl;
      return 1;
    }
  
  // -------------------------------------------------------------------------
  // Read some arguments.
  uint samplingSizeMax = args.getOption( "-samplingSizeMax" )->getIntValue( 0 );
  uint nbIterationSpikes = args.getOption( "-nbIterationSpikes" )->getIntValue( 0 );
  uint mscales_min_size = args.getOption( "-meaningfulScales" )->getIntValue( 0 );
  double mscales_max_slope = args.getOption( "-meaningfulScales" )->getDoubleValue( 1 );

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
  FreemanChainSubsample fcsub( 1, 1, 0, 0 );
  FreemanChainCleanSpikesCCW fccs( nbIterationSpikes );
  FreemanChainCompose fcomp( fccs, fcsub );
  FreemanChainTransform* ptr_fct = &fcomp;
  FreemanChainSubsample* ptr_fcsub = &fcsub;

  MultiscaleProfile MP;
  MP.chooseSubsampler( *ptr_fct, *ptr_fcsub );
  MP.init( c, samplingSizeMax );

  if ( args.check( "-affMapping" ) )
    {
      uint r = args.getOption( "-affMapping" )->getIntValue( 0 );
      int x0 = args.getOption( "-affMapping" )->getIntValue( 1 );
      int y0 = args.getOption( "-affMapping" )->getIntValue( 2 );
      cout << "# -affMapping: mapping for subsampling ("
	   << r << "," << r << "," << x0 << "," << y0 << ")" << endl
	   << "# pointel_idx surfel_idx sub_pointel_idx sub_surfel_idx"
	   << endl;
      MultiscaleFreemanChain::SubsampledChainKey key( r, r, x0, y0 );
      const MultiscaleFreemanChain::SubsampledChain* subchain = MP.get( key );
      if ( subchain != 0 )
	{
	  for ( uint i = 0; i < subchain->c2subc.size(); ++i )
	    {
	      cout << i 
		   << " " << (i+0.5) 
		   << " " << subchain->subpointel( i )
		   << " " << ( subchain->subsurfel( i ).first
			       + subchain->subsurfel( i ).second ) 
		   << endl;
	    }
	}
    }
  if ( args.check( "-profile" ) )
    {
      uint idx = args.getOption( "-profile" )->getIntValue( 0 );
      vector<double> x;
      vector<double> y;
      MP.profile( x, y, idx );
      cout << "# -profile: Profile of surfel idx=" << idx << endl
	   << "# scale x y slope*x+shift" << endl;
      Statistics stats( 2, true );
      stats.addValues( 0, x.begin(), x.end() );
      stats.addValues( 1, y.begin(), y.end() );
      stats.terminate();
      std::pair< double, double > reg 
	= stats.linearRegression( 0, 1 );
      cout << "# slope=" << reg.first 
	   << " shift=" << reg.second << endl;
      uint s = x.size();
      for ( uint i = 0; i != s; ++i )
	{      
       	  cout << i 
	       << " " << x[ i ] << " " << y[ i ]
 	       << " " << ( x[ i ] * reg.first + reg.second )
 	       << endl;
	}
      std::vector< std::pair< uint, uint > > scales;
      MP.meaningfulScales( scales, idx, 
			   mscales_min_size, mscales_max_slope );
      cout << "# " << scales.size() << " meaningful scales" << endl
	   << "#";
      for ( uint i = 0; i != scales.size(); ++i )
	cout << " (" << scales[ i ].first 
	     << "," << scales[ i ].second
	     << ")";
      cout << endl;
    }

  if ( args.check( "-detailedProfile" ) )
    {
      uint idx = args.getOption( "-detailedProfile" )->getIntValue( 0 );
      vector<double> x;
      vector<double> y;
      vector<uint> nb;
      MP.detailedProfile( x, y, nb, idx );
      cout << "# -profile: Profile of surfel idx=" << idx << endl
	   << "# scale x y slope*x+shift" << endl;
      Statistics stats( 2, true );
      stats.addValues( 0, x.begin(), x.end() );
      stats.addValues( 1, y.begin(), y.end() );
      stats.terminate();
      std::pair< double, double > reg 
	= stats.linearRegression( 0, 1 );
      cout << "# slope=" << reg.first 
	   << " shift=" << reg.second << endl;
      uint s = x.size();
      for ( uint i = 0; i != s; ++i )
	{      
       	  cout << i 
	       << " " << x[ i ] << " " << y[ i ]
 	       << " " << ( x[ i ] * reg.first + reg.second )
 	       << endl;
	}
    }

  if ( args.check( "-affNoiseLevel" ) )
    {
      cout << "# -affNoiseLevel: displays noise level for each surfel." 
	   << endl
	   << "# idx noiselvl code x y" << endl;
      for ( FreemanChain::const_iterator it = c.begin();
	    it != c.end();
	    ++it )
	{
	  uint idx = it.getPosition();
	  uint code = it.getCode();
	  Vector2i xy( *it );
	  uint nl = MP.noiseLevel( idx, mscales_min_size, mscales_max_slope );
	  cout << idx << " " << nl << " " << code
	       << " " << xy.x() << " " << xy.y() << endl;
	}
    }
  // testLinearRegression( 1.5 );
  // testLinearRegression( 1.8 );
  // testLinearRegression( -0.3 );

  if ( args.check( "-affTgt" ) )
    {
      uint r = args.getOption( "-affTgt" )->getIntValue( 0 );
      int x0 = args.getOption( "-affTgt" )->getIntValue( 1 );
      int y0 = args.getOption( "-affTgt" )->getIntValue( 2 );
      cout << "# -affTgt: tangent estimation for subsampling ("
	   << r << "," << r << "," << x0 << "," << y0 << ")" << endl
	   << "# idx subc_idx tgt_angle"
	   << endl;
      MultiscaleFreemanChain::SubsampledChainKey key( r, r, x0, y0 );
      const MultiscaleFreemanChain::SubsampledChain* subchain = MP.get( key );
      if ( subchain != 0 )
	{
	  LambdaMSTEstimator LMST;
	  FreemanChain copyc;
	  copyc.chain = subchain->subc.chain;
	  copyc.x0 = subchain->subc.x0;
	  copyc.y0 = subchain->subc.y0;
	  KnSpace* ks 
	    = ShapeHelper::makeSpaceFromFreemanChain( copyc );
	  LMST.init( *ks );
	  C4CIteratorOnFreemanChain itfc;
	  C4CIteratorOnFreemanChainSurface itfcs;
	  for ( uint i = 0; i < subchain->c2subc.size(); ++i )
	    {
	      pair<uint,double> surfel = subchain->subsurfel( i );
	      uint surf_idx = surfel.first;
	      double shift = surfel.second;
	      FreemanChain::const_iterator p ( copyc,
					       surf_idx );
	      itfc.init( p, true );
	      itfcs.init( ks, itfc );
	      cout << i 
		   << " " << ( surf_idx + shift )
		   << " " << LMST.estimation( itfcs, shift )
		   << endl;
	    }
	  delete ks;
	}
    }
  if ( args.check( "-multiresTgt" ) )
    {
      uint multires_mode = 
	( args.getOption( "-multiresTgt" )->getValue( 0 ) == "MEAN" ) ? 0 
	: ( args.getOption( "-multiresTgt" )->getValue( 0 ) == "MAX" ) ? 1 
	: 2;
      map<MultiscaleFreemanChain::SubsampledChainKey, LambdaMSTEstimator>
	map_estimators;      
      createLambdaMSTEstimators( map_estimators, MP, samplingSizeMax );
      C4CIteratorOnFreemanChain itfc;
      C4CIteratorOnFreemanChainSurface itfcs;
      for ( FreemanChain::const_iterator it = c.begin();
	    it != c.end();
	    ++it )
	{
	  uint idx = it.getPosition();
	  uint code = it.getCode();
	  Vector2i xy( *it );
	  uint nl = MP.noiseLevel( idx, mscales_min_size, mscales_max_slope );
	  MultiscaleFreemanChain::SubsampledChainKey key
	    = ( multires_mode == 0 )
	    ? MP.longestMeanAtScale( idx, nl )
	    : ( multires_mode == 1 ) 
	    ? MP.longestMaxAtScale( idx, nl )
	    : MultiscaleFreemanChain::SubsampledChainKey( nl, nl, 0, 0 );
	  // if ( nl != 0 )
	  //   cerr << "nl=" << nl << " key=" 
	  // 	 << key.h << "," << key.v << "," << key.x0 << ","
	  // 	 << key.y0 << endl;
	  const MultiscaleFreemanChain::SubsampledChain* subchain
	    = MP.get( key );
	  const LambdaMSTEstimator & LMLP
	    = map_estimators[ key ];
	  // cerr << "MULTIRESTGT key=(" << key.h << "," << key.v 
	  //      << "," << key.x0 << "," << key.y0 << ") LMST=[ks=("
	  //      << LMLP.space()->size( 0 ) << "," << LMLP.space()->size( 1 ) << ") chain="
	  //      << subchain->subc.chain << "]" << endl;

	  pair<uint,double> surfel = subchain->subsurfel( idx );
	  uint surf_idx = surfel.first;
	  double shift = surfel.second;
	  FreemanChain::const_iterator p ( subchain->subc,
					   surf_idx );
	  itfc.init( p, true );
	  itfcs.init( LMLP.space(), itfc );

	  cout << idx 
	       << " " << ( surf_idx + shift )
	       << " " << LMLP.estimation( itfcs, shift )
	       << " " << nl << " " << code
	       << " " << xy.x() << " " << xy.y() << endl;
	}
      // to do : delete KnSpaces
    }

  if ( args.check( "-averageMultiresTgt" ) )
    {
      map<MultiscaleFreemanChain::SubsampledChainKey, LambdaMSTEstimator>
	map_estimators;      
      createLambdaMSTEstimators( map_estimators, MP, samplingSizeMax );
      for ( FreemanChain::const_iterator it = c.begin();
	    it != c.end();
	    ++it )
	{
	  uint idx = it.getPosition();
	  uint code = it.getCode();
	  Vector2i xy( *it );
	  uint nl = MP.noiseLevel( idx, mscales_min_size, mscales_max_slope );
	  double tgt = averageEstimation( map_estimators, MP, idx, nl );

	  cout << idx 
	       << " " << (idx+0.5)
	       << " " << tgt
	       << " " << nl << " " << code
	       << " " << xy.x() << " " << xy.y() << endl;
	}
      // to do : delete KnSpaces
    }

  if ( args.check( "-averageTgt" ) )
    {
      uint scale = args.getOption( "-averageTgt" )->getIntValue( 0 );
      map<MultiscaleFreemanChain::SubsampledChainKey, LambdaMSTEstimator>
	map_estimators;      
      createLambdaMSTEstimators( map_estimators, MP, samplingSizeMax );
      for ( FreemanChain::const_iterator it = c.begin();
	    it != c.end();
	    ++it )
	{
	  uint idx = it.getPosition();
	  uint code = it.getCode();
	  Vector2i xy( *it );
	  double tgt = averageEstimation( map_estimators, MP, idx, scale );

	  cout << idx 
	       << " " << (idx+0.5)
	       << " " << tgt
	       << " " << scale << " " << code
	       << " " << xy.x() << " " << xy.y() << endl;
	}
      // to do : delete KnSpaces
    }

  if ( args.check( "-gsmooth" ) )
    {
      double s = args.getOption( "-gsmooth" )->getDoubleValue( 0 );
      cout << "# -gsmoothContour: displays a smoother version of the contour with gaussian kernel" 
	   << endl << "# " << s << endl
	   << "# idx sx sy noiselvl x y" << endl;
      vector<double> x;
      vector<double> y;
      vector<double> sx;
      vector<double> sy;
      vector<uint> nl;
      for ( FreemanChain::const_iterator it = c.begin();
	    it != c.end();
	    ++it )
	{
	  uint idx = it.getPosition();
	  uint code = it.getCode();
	  Vector2i xy( *it );
	  x.push_back( (double) xy.x() );
	  y.push_back( (double) xy.y() );
	  nl.push_back( 0 );
	}
      vector<G> g_fcts( 1 );
      g_fcts[ 0 ].setSigma( s );
      convolveByGaussianKernels( sx, sy, x, y, nl, g_fcts );
      for ( uint i = 0; i < x.size(); ++i )
	cout << i << " " << sx[ i ] << " " << sy[ i ]
	     << " " << nl[ i ] << " " << x[ i ] << " " << y[ i ] << endl;
    

    }

  if ( args.check( "-multiresMedianPGM" ) )
    {
      multiresMedianPGM( MP, mscales_min_size, mscales_max_slope );
    }
  if ( args.check( "-multiresMedianPGMOnOriginalImage" ) )
    {
      string filename = args.getOption("-multiresMedianPGMOnOriginalImage")
	->getValue(0);
      ifstream in_str;
      in_str.open(filename.c_str(), ios::in);
      if ( in_str.good() )
	{
	  KnSpace* ks = 0;
	  KnRUCellVector<int>* vec = 0;
	  bool ok = ShapeHelper::importFromPGM( in_str, ks, vec );
	  if ( ok )
	    multiresMedianPGMOnOriginalImage( MP, 
					      c,
					      mscales_min_size,
					      mscales_max_slope,
					      ks, 
					      *vec );
	  else
	    cerr << "-multiresMedianPGMOnOriginalImage: error reading " << filename << endl;
	  if ( vec != 0 ) delete vec;
	  if ( ks != 0 ) delete ks;
	}
      else
	cerr << "-multiresMedianPGMOnOriginalImage: error opening " << filename << endl;
    }

}


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

void
testLinearRegression( double p )
{
  std::vector<Vector2D> pts;
  Statistics tstat( 2, true );
  double x = 10.0;
  for ( uint i = 1; i < 1000; ++i )
    {
      double y = pow( x, p );
      pts.push_back( Vector2D( log( x ), log( y ) ) );
      tstat.addValue( 0, pts.back().x() );
      tstat.addValue( 1, pts.back().y() );
      x += 5.0;
    } 
  tstat.terminate();
  Vector2D reg = getRegLinear( pts );
  std::pair<double,double> sreg = tstat.linearRegression( 0, 1 );
  cout << "+--- test linear regression " << endl
       << " exp=" << p
       << " slo=" << reg.x()
       << " shi=" << reg.y()
       << " sslo=" << sreg.first
       << " sshi=" << sreg.second
       << endl;
}

double 
convolution( const vector<double> & x, const G & g,
	     uint t )
{
  uint nx = x.size();
  double acc = x[ t ] * g( 0.0 );
  double acc2 = g( 0.0 );
  // cerr << "t=" << t << " g(t)=" << g(0.0);
  for ( uint i = 1; i < (nx/2); ++i )
    {
      double pos = (double) i;
      double gpos = g( pos );
      uint k1 = ( t + i ) % nx;
      uint k2 = ( t + nx - (int) i ) % nx;
      // ( t - k ) >= (double) nx
      // 	? t - k - nx
      // 	: ( ( t < k )
      // 	    ? t - k + nx
      // 	    : t - k );
      acc += x[ k1 ] * gpos;
      acc += x[ k2 ] * gpos;
      acc2 += gpos;
      acc2 += gpos;
      // cerr << " k1=" << k1 << " g(k1)=" <<  g( - (double) i )
      // 	   << " k2=" << k2 << " g(k2)=" <<  g( (double) i );
    }
  // cerr << endl;
  return acc / acc2;
}

void
convolveByGaussianKernels( vector<double> & sx,
			   vector<double> & sy,
			   const vector<double> & x,
			   const vector<double> & y,
			   const vector<uint> & nl,
			   const vector<G> & g_fcts )
{
  for ( uint t = 0; t < x.size(); ++t )
    {
      sx.push_back( convolution( x, g_fcts[ nl[ t ] ], t ) );
      sy.push_back( convolution( y, g_fcts[ nl[ t ] ], t ) );
    }
}



void 
multiresMedianPGMOnOriginalImage
( const MultiscaleProfile & MP,
  const FreemanChain & c,
  uint mscales_min_size, 
  double mscales_max_slope,
  KnSpace* ks,
  KnRUCellVector<int> & vectorVal)
{
  MultiscaleFreemanChain::SubsampledChainKey key( 1, 1, 0, 0 );
  const MultiscaleFreemanChain::SubsampledChain* subchain = MP.get( key );
  // JOL: problem with coordinates.
//   FreemanChain c;
//   c.chain = subchain->subc.chain;
//   c.x0 = subchain->subc.x0;
//   c.y0 = subchain->subc.y0;

//  cerr << c.chain << " " << c.x0 << " " << c.y0 << endl;
//   K2Space* ks = ShapeHelper::makeSpaceFromFreemanChain( c );
//   KnRCellSet contour 
//     = ShapeHelper::makeContourFromFreemanChain( ks, c, true );
//   KnCharSet image = KnShapes::ucomputeInterior( *ks, contour );
  Kn_size w = ks->size( 0 );
  Kn_size h = ks->size( 1 );
  Kn_uid aspel = ks->ufirstCell( 2 );
  int input[ w ][ h ];
  int nlmask[ w ][ h ];
  int output[ w ][ h ];

  // Rajout BK Expe PAMI
  int maxScale = MP.all_stats.size();
  
  bool cstNoise = args.check("-setCstNoise");
  int levelNoise;
  int minSize;
  double alpha;
  

  if(cstNoise){
    levelNoise = args.getOption("-setCstNoise")->getIntValue(0);
  }
  //---

  for ( Kn_size y = 0; y < h; ++y )
    for ( Kn_size x = 0; x < w; ++x )
      {
	Kn_size tbl[ 2 ];
	tbl[ 0 ] = x;
	tbl[ 1 ] = y;
	Kn_uid pixel = ks->ucode( tbl, aspel );
	input[ x ][ y ] = vectorVal[ pixel ];
	nlmask[ x ][ y ] = 0;
      }
  for ( FreemanChain::const_iterator it = c.begin(); it != c.end(); ++it )
    {
      Vector2i p = *it;
      int nlv=0;
            
      nlv = (int) MP.noiseLevel( it.getPosition(), 
				 mscales_min_size, mscales_max_slope);
      
      
      // Rajout pour effectuer les tests dans PAMI:
      // Meaningfull scales par défaut, StdScale ou niveau constant.
      
      if(cstNoise){
	nlv= levelNoise;	
      }
      //cerr << p.x() << "," << p.y() << "=" << nlv << endl;
      
      //nlv = nlv > 0 ? nlv - 1 : nlv;
      nlv = nlv > 0 ? nlv - 1 : maxScale;
      for ( int k = -nlv; k <= nlv; ++k )
	for ( int l = -nlv; l <= nlv; ++l )
	  {
	    if ( ( ( p.x() + k ) >= 0 ) && ( ( p.y() + l ) >= 0 )
		 && ( ( p.x() + k ) < w ) && ( ( p.y() + l ) < h ) )
	      {
		int absk = ( k < 0) ? -k : k;
		int absl = ( l < 0) ? -l : l;
		int d = absk < absl ? absl : absk;
		if ( nlmask[ p.x() + k ][ p.y() + l ] < ( nlv - d ) )
		  nlmask[ p.x() + k ][ p.y() + l ] = nlv - d;
	      }
	  }
    }
  for ( uint y = 0; y < h; ++y )
    for ( uint x = 0; x < w; ++x )
      {
	vector<int> val;
	int nlv = nlmask[ x ][ y ];
	int max = (nlv+1) / 2;
	int min = -max;
	for ( int k = min; k <= max; ++k )
	  for ( int l = min; l <= max; ++l )
	    {
	      if ( ( ( x + k ) >= 0 ) && ( ( y + l ) >= 0 )
		   && ( ( x + k ) < w ) && ( ( y + l ) < h ) )
		val.push_back( input[ x + k ][ y + l ] );
	    }
	uint med = val.size() / 2;
	nth_element( val.begin(), val.begin() + med, val.end() );
	output[ x ][ y ] = val[ med ];
      }

  for ( Kn_size y = 0; y < h; ++y )
    for ( Kn_size x = 0; x < w; ++x )
      {
	Kn_size tbl[ 2 ];
	tbl[ 0 ] = x;
	tbl[ 1 ] = y;
	Kn_uid pixel = ks->ucode( tbl, aspel );
	//Kn_uid pixel = ks->upixel( x, y );
	vectorVal[ pixel ] = output[ x ][ y ];
      }
  ShapeHelper::exportToPGM( cout, ks, vectorVal );

}

void 
multiresMedianPGM( const MultiscaleProfile & MP,
		   uint mscales_min_size, 
		   double mscales_max_slope )
{
  MultiscaleFreemanChain::SubsampledChainKey key( 1, 1, 0, 0 );
  const MultiscaleFreemanChain::SubsampledChain* subchain = MP.get( key );
  FreemanChain c;
  c.chain = subchain->subc.chain;
  c.x0 = 0;
  c.y0 = 0;
  
  K2Space* ks = ShapeHelper::makeSpaceFromFreemanChain( c );
  KnRCellSet contour 
    = ShapeHelper::makeContourFromFreemanChain( ks, c, true );
  KnCharSet image = KnShapes::ucomputeInterior( *ks, contour );
  Kn_size w = ks->width();
  Kn_size h = ks->height();
  bool input[ w ][ h ];
  int nlmask[ w ][ h ];
  bool output[ w ][ h ];

  // Rajout BK Expe PAMI
  int maxScale = MP.all_stats.size();
  
  bool cstNoise = args.check("-setCstNoise");
  int levelNoise;
  int minSize;
  double alpha;
  
  if(cstNoise){
    levelNoise = args.getOption("-setCstNoise")->getIntValue(0);
  }
  //---

  for ( Kn_size y = 0; y < h; ++y )
    for ( Kn_size x = 0; x < w; ++x )
      {
	Kn_uid pixel = ks->upixel( x, y );
	input[ x ][ y ] = image[ pixel ];
	nlmask[ x ][ y ] = 0;
      }
  for ( FreemanChain::const_iterator it = c.begin(); it != c.end(); ++it )
    {
      Vector2i p = *it;
      int nlv=0;
            
      nlv = (int) MP.noiseLevel( it.getPosition(), 
				 mscales_min_size, mscales_max_slope);
      
      
      if(cstNoise){
	nlv= levelNoise;	
      }
      
      //nlv = nlv > 0 ? nlv - 1 : nlv;
      nlv = nlv > 0 ? nlv - 1 : maxScale;
      for ( int k = -nlv; k <= nlv; ++k )
	for ( int l = -nlv; l <= nlv; ++l )
	  {
	    if ( ( ( p.x() + k ) >= 0 ) && ( ( p.y() + l ) >= 0 )
		 && ( ( p.x() + k ) < w ) && ( ( p.y() + l ) < h ) )
	      {
		int absk = ( k < 0) ? -k : k;
		int absl = ( l < 0) ? -l : l;
		int d = absk < absl ? absl : absk;
		if ( nlmask[ p.x() + k ][ p.y() + l ] < ( nlv - d ) )
		  nlmask[ p.x() + k ][ p.y() + l ] = nlv - d;
	      }
	  }
    }
  for ( uint y = 0; y < h; ++y )
    for ( uint x = 0; x < w; ++x )
      {
	vector<int> val;
	int nlv = nlmask[ x ][ y ];
	int max = (nlv+1) / 2;
	int min = -max;
	for ( int k = min; k <= max; ++k )
	  for ( int l = min; l <= max; ++l )
	    {
	      if ( ( ( x + k ) >= 0 ) && ( ( y + l ) >= 0 )
		   && ( ( x + k ) < w ) && ( ( y + l ) < h ) )
		val.push_back( input[ x + k ][ y + l ]
			       ? 1 : 0 );
	    }
	uint med = val.size() / 2;
	nth_element( val.begin(), val.begin() + med, val.end() );
	output[ x ][ y ] = val[ med ];
      }

  for ( Kn_size y = 0; y < h; ++y )
    for ( Kn_size x = 0; x < w; ++x )
      {
	Kn_uid pixel = ks->upixel( x, y );
	image[ pixel ] = output[ x ][ y ];
      }
  ShapeHelper::exportToPGM( cout, ks, image );
  delete ks;
  

}



