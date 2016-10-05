///////////////////////////////////////////////////////////////////////////////
// Test the length variation of maximal segments on images in order to
// detect its meaningful scales.
///////////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <vector>
#include "ImaGene/base/Vector2i.h"
#include "ImaGene/base/Arguments.h"
#include "ImaGene/base/StandardArguments.h"
#include "ImaGene/mathutils/G.h"
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
#include "ImaGene/helper/DigitalElevationMap.h"
#include "ImaGene/helper/DigitalElevationMapStats.h"
#include "ImaGene/helper/ImageMultiscaleAnalysis.h"
#include "ImaGene/image/Image2D.h"
#include "ImaGene/image/PGMFilter.h"
#include "ImaGene/image/Image2DUtils.h"


using namespace std;
using namespace ImaGene;


static Arguments args;

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
  StandardArguments::addIOArgs( args, true, true );
  args.addOption( "-gaussian_smooth", "-gaussian_smooth <sigma>: smoothes the input image by convolution with the gaussian function of variance sigma.", "0.15915494309189533577"  );
  args.addOption( "-subsample", "-subsample <h> <v> <x0> <y0>: subsample the image.", "1", "1", "0", "0" );
  args.addOption( "-quantification", "-quantification <z>: quantification when subsampling.", "1" );
  args.addOption( "-meaningfulScales", "-meaningfulScales <min_size> <max_slope> <min_slope>: specifies parameters for defining meaningful scales: minimum size of the interval of scales and maximum/minimum slopes between consecutive samples within.", "1", "-0.2", "-10.0" );
  args.addBooleanOption( "-dPGM", "-dPGM: outputs the computed image." );
  args.addBooleanOption( "-dEMap", "-dEmap: outputs the elevation map as a 3D object in PGM3D." );
  args.addBooleanOption( "-MSX", "-MSX: displays the number of maximal segments along X." );
  args.addOption( "-IMA", "-IMA <max>: computes the image multiscale analysis till scale <max>.", "10" );

  if ( ( argc <= 0 ) 
       || ! args.readArguments( argc, argv ) ) 
    {
      cerr << args.usage( "test_ImageScaleAnalysis", 
			  "Test the length variation of maximal segments on images in order to detect its meaningful scales."
			  ,"" ) << endl;
      return 1;
    }
  
  // -------------------------------------------------------------------------
  // Read some arguments.

  uint mscales_min_size 
    = args.getOption( "-meaningfulScales" )->getIntValue( 0 );
  double mscales_max_slope 
    = args.getOption( "-meaningfulScales" )->getDoubleValue( 1 );
  double mscales_min_slope
    = args.getOption( "-meaningfulScales" )->getDoubleValue( 2 );

  // -------------------------------------------------------------------------
  // Read Freeman chain and creates space/contour.  
  
  istream & in_str = StandardArguments::openInput( args );
  ByteImage2D input_image;
  ByteImage2D output_image;
  Image2DUtils<ByteImage2D> image_utils;
  PGMFilter<ByteImage2D> pgm_filter;
//   Image2D< uint8 > image;
//   PGMFilter< Image2D< uint8 > >  pgm_filter;
  bool success = pgm_filter.read( input_image, in_str );
  if ( ! success )
    {
      cerr << "[test_ImageScaleAnalysis] Error reading PGM image: " 
	   << pgm_filter.errorMessage() << endl;
      return 2;
    }
  input_image.copy( output_image );
  if ( args.check( "-gaussian_smooth" ) )
    {
      double sigma = args.getOption( "-gaussian_smooth" )->getDoubleValue( 0 );
      ByteImage2D smoothed_image;
      image_utils.gaussianSmoothing( smoothed_image, input_image, sigma );
      smoothed_image.copy( output_image );
    }
  if ( args.check( "-subsample" ) )
    {
      uint dx = args.getOption( "-subsample" )->getIntValue( 0 );
      uint dy = args.getOption( "-subsample" )->getIntValue( 1 );
      int x0 = args.getOption( "-subsample" )->getIntValue( 2 );
      int y0 = args.getOption( "-subsample" )->getIntValue( 3 );
      uint z = args.getOption( "-quantification" )->getIntValue( 0 );
      ByteImage2D subsampled_image;
      Pixel origin( x0, y0 );
      image_utils.subsample( subsampled_image, output_image, 
			     dx, dy, origin, z );
      subsampled_image.copy( output_image );
    }
  if ( args.check( "-dEMap" ) )
    {
      DigitalElevationMap emap;
      emap.init( output_image );
      ostream & out_str = StandardArguments::openOutput( args );
      ShapeHelper::exportToPGM3d( out_str,
				  emap.space(),
				  emap.object() );
    }
  if ( args.check( "-MSX" ) )
    {
      // DigitalElevationMap emap;
      // emap.init( output_image );
      // Pixel p = emap.lowest();
      // Pixel q = emap.highest();
      // Pixel z = p;
      // for ( Pixel z = p; z.y < q.y; ++z.y )
      // 	{
      // 	  C4CIteratorOnSurface* iter = emap.createHIterator( z );
      // 	  C4CTangentialCover tcover;
      // 	  tcover.init( *iter, 0 );
      // 	  cerr << "Row " << z.y << ": " << tcover.nbMaximalSegments() << endl;
      // 	  delete iter;
      // 	}
      DigitalElevationMapStats emap;
      emap.init( output_image );
      emap.computeStats();
      Kn_size sizes[ 2 ];
      Pixel p = emap.elevation().lowest();
      Pixel q = emap.elevation().highest();
      for ( Pixel z(p.x + 1, p.y + 1); z.y < q.y; ++z.y )
	for ( z.x = p.x + 1; z.x < q.x; ++z.x )
	  {
	    sizes[ 0 ] = ( (uint) z.x - p.x ) * 2;
	    sizes[ 1 ] = ( (uint) z.y - p.y ) * 2 + 1;
	    Kn_uid linelx = emap.space2D()->ukcode( sizes );
	    sizes[ 0 ] = ( (uint) z.x - p.x ) * 2 + 1;
	    sizes[ 1 ] = ( (uint) z.y - p.y ) * 2;
	    Kn_uid linely = emap.space2D()->ukcode( sizes );
	    Statistic<float> s( emap.stat( linelx ) );
	    s += emap.stat( linely );
	    output_image.set( z, (int) s.mean() );
	  }
    }
  if ( args.check( "-IMA" ) )
    {
      uint max_scale = args.getOption( "-IMA" )->getIntValue( 0 );
      ImageMultiscaleAnalysis ima;
      ima.init( output_image, max_scale );
      Kn_size sizes[ 2 ];
      Pixel p = output_image.lowest();
      Pixel q = output_image.highest();
      ScaleProfile sp;
      for ( Pixel z(p.x + 1, p.y + 1); z.y < q.y; ++z.y )
	for ( z.x = p.x + 1; z.x < q.x; ++z.x )
	  {
	    sizes[ 0 ] = ( (uint) z.x - p.x ) * 2;
	    sizes[ 1 ] = ( (uint) z.y - p.y ) * 2 + 1;
	    ima.getScaleProfile( sp, sizes[ 0 ], sizes[ 1 ] );
	    uint nlx = sp.noiseLevel( mscales_min_size, 
				      mscales_max_slope, mscales_min_slope );
	    sizes[ 0 ] = ( (uint) z.x - p.x ) * 2 + 1;
	    sizes[ 1 ] = ( (uint) z.y - p.y ) * 2;
	    ima.getScaleProfile( sp, sizes[ 0 ], sizes[ 1 ] );
	    uint nly = sp.noiseLevel( mscales_min_size, 
				      mscales_max_slope, mscales_min_slope );

	    if ( ( nlx == 0 ) || ( nly == 0 ) )
	      output_image.set( z, 255 );
	    else
	      output_image.set( z, nlx + nly );
	  }      
    }
  if ( args.check( "-dPGM" ) )
    {
      ostream & out_str = StandardArguments::openOutput( args );
      bool success = pgm_filter.write( output_image, out_str );
      if ( ! success )
	{
	  cerr << "[test_ImageScaleAnalysis] Error writing PGM image: " 
	       << pgm_filter.errorMessage() << endl;
	  return 2;
	}
    }

  return 0;
}
