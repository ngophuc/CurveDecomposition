///////////////////////////////////////////////////////////////////////////////
// Test module for tangent computation along contours in 2D.
///////////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <iostream>
#include <iomanip>

//#include "LinAlg/LinAlg2D/Vector2D.hpp"
#include "ImaGene/mathutils/Mathutils.h"
#include "ImaGene/base/Proxy.h"
#include "ImaGene/base/Vector.h"
#include "ImaGene/base/Arguments.h"
#include "ImaGene/base/StandardArguments.h"
#include "ImaGene/base/Shapes.h"
#include "ImaGene/base/VectorUtils.h"
#include "ImaGene/timetools/Clock.h"
#include "ImaGene/dgeometry2d/C4CSegment.h"
#include "ImaGene/dgeometry2d/C4CIterator.h"
#include "ImaGene/dgeometry2d/C4CGeometry.h"
#include "ImaGene/dgeometry2d/C4CSegmentPencil.h"
#include "ImaGene/dgeometry2d/C4CTangentialCover.h"
#include "ImaGene/digitalnD/KnCharSet.h"
#include "ImaGene/digitalnD/BelAdjacency.h"
#include "ImaGene/digitalnD/ObjectBoundary.h"
#include "ImaGene/digitalnD/KnShapes.h"
#include "ImaGene/digitalnD/Embedder.h"
#include "ImaGene/digitalnD/GridEmbedder.h"
#include "ImaGene/digitalnD/Frame2D.h"
#include "ImaGene/digitalnD/KnSpaceScanner.h"
#include "ImaGene/helper/ShapeHelper.h"
#include "ImaGene/helper/C4CTangentialCoverGeometry.h"
#include "ImaGene/helper/GlobalC4CGeometry.h"
#include "ImaGene/mathutils/Statistics.h"
#include "ImaGene/base/Signal.h"


using namespace std;
using namespace ImaGene;

static Arguments args;



StarShaped* 
shapeFromArgs( KnSpace & ks, const Embedder & embedder,
	       KnCharSet & voxset,
	       Kn_sid & bel, uint & nb_bels )
{
  StarShaped* shape = ShapeHelper::makeStarShapedFromArgs( args );
  Vector2D dir( 1.0, 0.0 );
  cerr << "# --- digitizing shape " << endl;
  voxset = ShapeHelper::makeShape( ks, embedder, *shape );
  cerr << "# --- find start bel " << endl;
  bel = ShapeHelper::getBelOnStarShapedBoundary( ks, embedder, *shape, dir );
  cerr << "# --- compute contour size " << endl;
  nb_bels = KnShapes::sgetContourSize( ks, voxset, bel,
				       *( ks.sbegin_dirs( bel ) ) );
  
			  //cerr << "Set=" << voxset << endl;
  //cerr << "Bel=";
  //ks.displayKn_sid( bel, cerr );
  //cerr << endl;
  return shape;
}




///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//
// T E S T S
//
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

void testModuloComputer()
{
  Mathutils::ModuloComputer mc( 11 );
  uint j = 1;
  uint i = 0;
  do
    {
      cout << "[" << i << "<" << j << "(" << mc.k << ")]==" 
	   << mc.less( i, j ) << "  "
	   << "[" << i << ">" << j << "(" << mc.k << ")]==" 
	   << mc.less( j, i ) << endl;
      mc.increment( i );
    }
  while ( i != 0 );


  j = 3;
  i = 0;
  do
    {
      cout << "[" << j << "-" << i << "(" << mc.k << ")]==" 
	   << mc.posDiff( j, i )
	   << endl;
      mc.increment( i );
      mc.increment( j );
    }
  while ( i != 0 );
}



void buildTangentialCover( C4CTangentialCover & tcover, 
			   C4CIterator & cp, uint max_size )
{
  Clock::startClock();
  tcover.init( cp, max_size );
  long t = Clock::stopClock();
  cerr << "Tangential cover in " << t << " ms."
       << " nbsurf=" << tcover.nbSurfels() 
       << " nbms=" << tcover.nbMaximalSegments() << endl;
}


void testTangentialCoverGeometry( KnSpace* ks, C4CIteratorOnSurface & cp )
{
  C4CTangentialCover tcover;
  buildTangentialCover( tcover, cp, 0 );
  C4CTangentialCoverGeometry tcover_geometry;
  tcover_geometry.init( ks, 0, 1, tcover, cp );
  Mathutils::AngleComputer ac;

  float sum_diff = 0.0f;
  for ( uint i = 0; i < tcover_geometry.nbMaximalSegments(); ++i)
    {
      const C4CTangentialCoverGeometry::MaximalSegmentGeometry & msgeom
	= tcover_geometry.geometry( i );
      cout << "Geometry [" << i << "]"
	   << " angle=" << msgeom.angle_to_x
	   << " min=" << msgeom.min_angle
	   << " max=" << msgeom.max_angle
	   << endl;
      sum_diff += ac.posDiff( msgeom.max_angle, msgeom.min_angle );
    }
  cout << "Moyenne diff=" << sum_diff / tcover_geometry.nbMaximalSegments()
       << endl;
}

void displayCommonPartGeometry( KnSpace* ks, C4CIteratorOnSurface & cp,
				LinearMinimizer* lm )
{
  // build tangential cover
  C4CTangentialCover tcover;
  buildTangentialCover( tcover, cp, 0 );

  // build tangential cover geometry, using a surfel frame and
  // constraints given by leaning points.
  C4CTangentialCoverGeometry tcover_geometry;
  Kn_sid sbel = cp.current();
  Frame2D frame;
  frame.init( ks, 0, 1 );
  frame.setSurfelFrame( sbel, cp.trackDir() );
  C4CTangentialCoverGeometry::MSGeometryComputerByLP geocomp;
  tcover_geometry.init( tcover, cp, geocomp, frame );

  // compute curvature by global geometry. The length of each common
  // part is estimated with normal integration, normal estimated with
  // lambda-MS tangent estimator.
  GlobalC4CGeometry global_geom;
  TriangleFunction lambda;
  C4CELength* el = tcover_geometry.computeELengthByLambdaMS
	( tcover, tcover_geometry, lambda );

  global_geom.init( tcover, tcover_geometry, *el, lm );
  global_geom.computeContourGeometry
    ( args.getOption( "-eps" )->getDoubleValue( 0 ),
      args.getOption( "-eps" )->getDoubleValue( 1 ) );
      
  GlobalC4CGeometry::displayGeometry( cout, global_geom );
  delete el;
}

// OLD VERSION (JOL 2007/03/09)
// void displayCommonPartGeometry( KnSpace* ks, C4CIteratorOnSurface & cp )
// {
//   C4CTangentialCover tcover;
//   buildTangentialCover( tcover, cp, 0 );
//   C4CTangentialCoverGeometry tcover_geometry;
//   tcover_geometry.init( ks, 0, 1, tcover, cp );
//   GlobalC4CGeometry global_geom;
//   global_geom.init( tcover, tcover_geometry );
//   global_geom.computeContourGeometry( 0.00001 );
//   GlobalC4CGeometry::displayGeometry( cout, global_geom );
// }

void displayTangentialCoverGeometry( KnSpace* ks, C4CIteratorOnSurface & cp )
{
  C4CTangentialCover tcover;
  buildTangentialCover( tcover, cp, 0 );
  C4CTangentialCoverGeometry tcover_geometry;
  tcover_geometry.init( ks, 0, 1, tcover, cp );
  C4CTangentialCoverGeometry::displayGeometry( cout, tcover, tcover_geometry );
}

void displayShape( KnSpace* ks, C4CIteratorOnSurface & cp,
		   StarShaped* shape, GridEmbedder & embedder,
		   float dh )
{
  cout << "# Display shape information" << endl
       << "# surf_idx xc yc t angle_to_x curvature" << endl;
  
  Kn_sid sbel = cp.current();
  Kn_sid bel = sbel;
  float th_a_t_x;
  uint n = 0;
  do
    {
      Vector xbel( 2 );
      embedder.sembed( bel, xbel );
      float theta = shape->parameter( xbel );   // surfel's centroid angle to origin
      Vector2D t = shape->tangent( theta );     // theoretical tangent computed at surfel's centroid angle
      float curv = shape->curvature( theta );
      
      float c_t_a_x = VectorUtils::dotProduct( Vector2D( 1.0, 0.0 ) , t );
      c_t_a_x = c_t_a_x > 1.0 ? 1.0 : ( c_t_a_x < -1.0 ? -1.0 : c_t_a_x );
      float s_t_a_x = VectorUtils::det( Vector2D( 1.0, 0.0 ) , t );
      if ( s_t_a_x >= 0 ) th_a_t_x = acos( c_t_a_x );
      else th_a_t_x = 2 * M_PI - acos( c_t_a_x );

      cout << n << " " << xbel.ro( 0 ) << " " << xbel.ro( 1 )
	   << " " << theta << " " << th_a_t_x << " " << curv
	   << endl;
      cp.next();
      bel = cp.current();
      ++n;
      
    }
  while ( bel != sbel );
}


// Pour Bertrand Kerautret.
void computeCurvature( KnSpace* ks, const string & fname,
		       GridEmbedder & embedder,
		       float dh,
		       LinearMinimizer* lm )
{
  // TODO: read file fname

  // Cree l'ensemble des pixels interieurs.
  KnCharSet pixel_set = KnCharSet::create( *ks, ks->dim(), false, 0 );
  uint tx = ks->size( 0 );
  uint ty = ks->size( 1 );
  Kn_size coords[ 2 ];
  bool first_bel = true;
  Kn_sid bel = 0;
  for ( uint y = 0; y < ty; ++y )
    for ( uint x = 0; x < tx; ++x )
      {
	coords[ 0 ] = 2*x+1;
	coords[ 1 ] = 2*y+1;
	Kn_sid p = ks->ukcode( coords );
	// inside test
	if ( ( x > 10 ) && ( x < 20 ) && ( y > 10 ) && ( y < 25 ) )
	  {
	    pixel_set += p;
	    if ( first_bel )
	      {
		bel = ks->sincident( ks->signsPos( p ), 0, false );
		first_bel = false;
	      }
	  }
      }

  // Create boundary.
  BelAdjacency badj( *ks, true );
  ObjectBoundary bdry( badj, pixel_set );
  uint dir = *( ks->sbegin_dirs( bel ) );
  C4CIteratorOnSurface* cp = bdry.newC4CIterator( bel, dir );

  // Calcule l'ensemble des segments maximaux de la courbe.
  C4CTangentialCover tcover;
  buildTangentialCover( tcover, *cp, 0 );
  C4CTangentialCoverGeometry tcover_geometry;
  tcover_geometry.init( ks, 0, 1, tcover, *cp );
  GlobalC4CGeometry global_geom;
  TriangleFunction lambda;
  C4CELength* el = tcover_geometry.computeELengthByLambdaMS
	( tcover, tcover_geometry, lambda );
  global_geom.init( tcover, tcover_geometry, *el, lm );
  global_geom.computeContourGeometry
    ( args.getOption( "-eps" )->getDoubleValue( 0 ),
      args.getOption( "-eps" )->getDoubleValue( 1 ) );

  cout << "# Test global geometry" << endl
       << "# sidx dabs xc yc est_t est_c" << endl
       << "# Modulo surfel is " << tcover.nbSurfels() << endl
       << "# dh is " << dh << endl;

  Kn_sid sbel = cp->current();
  double dabs = 0.5;
  uint i = 0;
  do
    {
      Vector xbel( 2 );
      embedder.sembed( bel, xbel );
      GlobalC4CGeometry::LocalGeometry lgeo;
      global_geom.geometryFromDiscreteAbscissa( (float) dabs, lgeo );
      
      cout << i << " " << dabs << " "
	   << " " << xbel.ro( 0 ) << " " << xbel.ro( 1 )
	   << " " << lgeo.angle_to_x
	   << " " << lgeo.curvature / dh
	   << endl;

      cp->next();
      bel = cp->current();
      dabs += 1.0;
      if ( dabs >= (double) tcover.nbSurfels() ) 
	dabs -= (double) tcover.nbSurfels();
      ++i;
    }
  while ( bel != sbel );
  delete cp;
}




void testGlobalGeometry( KnSpace* ks, C4CIteratorOnSurface & cp,
			 StarShaped* shape, GridEmbedder & embedder,
			 double dh,
			 LinearMinimizer* lm )
{
  C4CTangentialCover tcover;
  buildTangentialCover( tcover, cp, 0 );
  Proxy<C4CIterator> bs( tcover.beginSurfel() );
  Proxy<C4CIterator> es( tcover.endSurfel() );
  if ( bs->equals( *es ) && ! tcover.isContourOpen() )
    cout << "# Contour is closed : ok" << endl;
  else if ( ! bs->equals( *es ) && tcover.isContourOpen() )
    cout << "# Contour is open : ok" << endl;
  else cout << "# Contour is " << tcover.isContourOpen() << " : ko" << endl;

  C4CTangentialCoverGeometry tcover_geometry;
  // tcover_geometry.init( ks, 0, 1, tcover, cp );
  Kn_sid sbel = cp.current();
  Frame2D frame;
  frame.init( ks, 0, 1 );
  frame.setSurfelFrame( sbel, cp.trackDir() );
  C4CTangentialCoverGeometry::MSGeometryComputerByLP geocomp;
  tcover_geometry.init( tcover, cp, geocomp, frame );

  GlobalC4CGeometry global_geom;
  //  global_geom.init( tcover, tcover_geometry );
  TriangleFunction lambda;
  Proxy<C4CELength> el ( tcover_geometry.computeELengthByLambdaMS
	( tcover, tcover_geometry, lambda ) );
  global_geom.init( tcover, tcover_geometry, *el, lm );

//   cerr << "LENGTH=" << el->length( 0, 0 ) << endl;
//   for ( uint i = 0; i < el->nbSurfels(); i++ )
//     cerr << "EL[" << i << "]=" << el->elength( i ) << endl;

  global_geom.computeContourGeometry
    ( args.getOption( "-eps" )->getDoubleValue( 0 ),
      args.getOption( "-eps" )->getDoubleValue( 1 ) );


  cout << setprecision(15);
  cout << "#######################################################" << endl;
  cout << "# testGlobalGeometry" << endl
       << "# sidx dabs theta xc yc t est_t c est_c cabs length" << endl
       << "# Modulo surfel is " << tcover.nbSurfels() << endl
       << "# dh is " << dh << endl;

  Kn_sid bel = sbel;
  double cabs = 0.0; //el->elength( 0 ) / 2.0;
  double dabs = 0.5;
  float th_a_t_x;
  uint i = 0;
  do
    {
      Vector xbel( 2 );
      embedder.sembed( bel, xbel );
      float theta = shape->parameter( xbel );   // surfel's centroid angle to origin
      Vector2D t = shape->tangent( theta );     // theoretical tangent computed at surfel's centroid angle
      float curv = shape->curvature( theta );

      float c_t_a_x = VectorUtils::dotProduct( Vector2D( 1.0, 0.0 ) , t );
      c_t_a_x = c_t_a_x > 1.0 ? 1.0 : ( c_t_a_x < -1.0 ? -1.0 : c_t_a_x );
      float s_t_a_x = VectorUtils::det( Vector2D( 1.0, 0.0 ) , t );
      if ( s_t_a_x >= 0 ) th_a_t_x = acos( c_t_a_x );
      else th_a_t_x = 2 * M_PI - acos( c_t_a_x );
      
      GlobalC4CGeometry::LocalGeometry lgeo;
      global_geom.geometryFromDiscreteAbscissa( (double) dabs, lgeo );
      
      cout << i << " " << dabs << " " << theta 
	   << " " << xbel.ro( 0 ) << " " << xbel.ro( 1 )
	   << " " << th_a_t_x << " " << lgeo.angle_to_x
	   << " " << curv << " " << lgeo.curvature / dh
	   << " " << cabs << " " << el->elength( i )
	   << endl;

      cp.next();
      bel = cp.current();
      dabs += 1.0;
      cabs += el->elength( i );
      if ( dabs >= (double) tcover.nbSurfels() ) 
	dabs -= (double) tcover.nbSurfels();
      ++i;
    }
  while ( bel != sbel );
}

void testGlobalGeometryMinimizer( KnSpace* ks, C4CIteratorOnSurface & cp,
				  StarShaped* shape, GridEmbedder & embedder,
				  double dh,
				  LinearMinimizer* lm )
{
  C4CTangentialCover tcover;
  buildTangentialCover( tcover, cp, 0 );
  Proxy<C4CIterator> bs( tcover.beginSurfel() );
  Proxy<C4CIterator> es( tcover.endSurfel() );

  C4CTangentialCoverGeometry tcover_geometry;
  Kn_sid sbel = cp.current();
  Frame2D frame;
  frame.init( ks, 0, 1 );
  frame.setSurfelFrame( sbel, cp.trackDir() );
  C4CTangentialCoverGeometry::MSGeometryComputerByLP geocomp;
  tcover_geometry.init( tcover, cp, geocomp, frame );

  GlobalC4CGeometry global_geom;
  TriangleFunction lambda;
  Proxy<C4CELength> el ( tcover_geometry.computeELengthByLambdaMS
	( tcover, tcover_geometry, lambda ) );
  global_geom.init( tcover, tcover_geometry, *el, lm );

//   cerr << "LENGTH=" << el->length( 0, 0 ) << endl;
//   for ( uint i = 0; i < el->nbSurfels(); i++ )
//     cerr << "EL[" << i << "]=" << el->elength( i ) << endl;

  
  Clock::startClock();
  global_geom.computeContourGeometry
    ( args.getOption( "-eps" )->getDoubleValue( 0 ),
      args.getOption( "-eps" )->getDoubleValue( 1 ) );
  long t = Clock::stopClock();
  cout << "# GlobalGeometry in " << t << " ms." << endl;

  cout << setprecision(15);
  cout << "#######################################################" << endl;
  cout << "# testGlobalGeometryMinimizer" << endl
       << "# nbsurfels nbms nbparams " << endl
       << "# " << tcover.nbSurfels() << " " << tcover.nbMaximalSegments() 
       << " " << global_geom.nbCommonParts() << endl
       << "# Modulo surfel is " << tcover.nbSurfels() << endl
       << "# dh is " << dh << endl
       << "# average shape curvature is TK_avg." << endl
       << "# eps_sigma eps_mu TK_avg TK_M TK_m K_avg K_M K_m (K_M-K_m)/TK_avg"
       << "# meanABSERROR maxABSERROR" << endl;

  Statistics stats( 3 ); // 0 true curvature, 1 estimated curvature, 2 abs diff

  Kn_sid bel = sbel;
  double dabs = 0.5;
  do
    {
      Vector xbel( 2 );
      embedder.sembed( bel, xbel );
      float theta = shape->parameter( xbel );   // surfel's centroid angle to origin
      Vector2D t = shape->tangent( theta );     // theoretical tangent computed at surfel's centroid angle
      double curv = (double) shape->curvature( theta );

      GlobalC4CGeometry::LocalGeometry lgeo;
      global_geom.geometryFromDiscreteAbscissa( (double) dabs, lgeo );

      stats.addValue( 0, curv );
      stats.addValue( 1, lgeo.curvature / dh );
      stats.addValue( 2, fabs( curv - lgeo.curvature / dh ) );

      dabs += 1.0;
      if ( dabs >= (double) tcover.nbSurfels() ) 
	dabs -= (double) tcover.nbSurfels();
      
      cp.next();
      bel = cp.current();
    }
  while ( bel != sbel );
  stats.terminate();
  cout << args.getOption( "-eps" )->getDoubleValue( 1 )
       << " " << args.getOption( "-eps" )->getDoubleValue( 0 )
       << " " << stats.mean( 0 )
       << " " << stats.max( 0 )
       << " " << stats.min( 0 )
       << " " << stats.mean( 1 )
       << " " << stats.max( 1 )
       << " " << stats.min( 1 )
       << " " << ( stats.max( 1 ) - stats.min( 1 ) )/stats.mean( 0 )
       << " " << stats.mean( 2 )
       << " " << sqrt( stats.unbiasedVariance( 2 ) )
       << " " << stats.max( 2 )
       << endl;

  // TOFINISH
}


void testSizeCommonParts( C4CIterator & c )
{
  C4CTangentialCover tcover;
  buildTangentialCover( tcover, c, 0 );
  Mathutils::ModuloComputer mcsu( tcover.nbSurfels() );
  Statistics stats( 1 );
  C4CTangentialCover::CommonPart cp = tcover.beginCP();
  uint first_idx = cp.back_surfel_idx;
  do 
    {
      uint s = mcsu.posDiff( cp.after_front_surfel_idx, cp.back_surfel_idx );
      stats.addValue( 0, s );
      tcover.nextCP( cp );
    }
  while ( cp.back_surfel_idx != first_idx );
  
  stats.terminate();
  cout << "#######################################################" << endl;
  cout << "# testSizeCommonParts" << endl
       << "# nbSurfels nbCP Stats_LongCP( moy min max var uvar )" << endl;
  cout << tcover.nbSurfels()
       << " " << stats.samples( 0 ) 
       << " " << stats.mean( 0 )
       << " " << stats.min( 0 )
       << " " << stats.max( 0 )
       << " " << stats.variance( 0 )
       << " " << stats.unbiasedVariance( 0 )
       << endl;
}


void testNbMaximalSegments( C4CIterator & cp )
{
  C4CTangentialCover tcover;
  buildTangentialCover( tcover, cp, 0 );
  Mathutils::ModuloComputer mcsu( tcover.nbSurfels() );
  Statistics stats( 1 );
  
  for ( uint i = 0; i < tcover.nbMaximalSegments(); ++i )
    {
      const C4CTangentialCover::MaximalSegment & ms 
	= tcover.getMaximalSegment( i );
      stats.addValue( 0, mcsu.posDiff( ms.front_surfel_idx + 1,
				       ms.back_surfel_idx ) );
//       cout << ms.dss 
// 	   << " " << ms.back_surfel_idx 
// 	   << " " << ms.front_surfel_idx 
// 	   << " " << ms.front_surfel_idx - ms.back_surfel_idx 
// 	   << endl;
    }

  stats.terminate();
  cout << "#######################################################" << endl;
  cout << "# testNbMaximalSegments" << endl
       << "# nbSurfels nbMS Stats_LongMS( moy min max var uvar )" << endl;
  cout << tcover.nbSurfels()
       << " " << stats.samples( 0 ) 
       << " " << stats.mean( 0 )
       << " " << stats.min( 0 )
       << " " << stats.max( 0 )
       << " " << stats.variance( 0 )
       << " " << stats.unbiasedVariance( 0 )
       << endl;
}


void testNbMSPerSurfel( C4CIterator & cp )
{
  C4CTangentialCover tcover;
  buildTangentialCover( tcover, cp, 0 );

  Statistics stats( 1 );
  Mathutils::ModuloComputer mcms( tcover.nbMaximalSegments() );
  
  C4CTangentialCover::SurfelMaximalSegments sms
    = tcover.beginSMS( 0 );
  for ( uint i = 0; i < tcover.nbSurfels(); ++i )
    {
      uint k = sms.begin_ms; //< index of the first MS containing the surfel
      uint l = sms.end_ms;   //< index of the MS just after the last MS containing the surfel
      stats.addValue( 0, mcms.posDiff( l, k ) );
//       cout << i
//  	   << " " << k
//  	   << " " << l
//  	   << " " << l - k
//  	   << endl;
      tcover.nextSMS( sms );
    }

  stats.terminate();
  cout << "#######################################################" << endl;
  cout << "# testNbMSPerSurfel" << endl
       << "# nbMS nbSurfels Stats_NbMSPerSurfel( moy min max var uvar )"
       << endl;
  cout << tcover.nbMaximalSegments()
       << " " << stats.samples( 0 ) 
       << " " << stats.mean( 0 )
       << " " << stats.min( 0 )
       << " " << stats.max( 0 )
       << " " << stats.variance( 0 )
       << " " << stats.unbiasedVariance( 0 )
       << endl;
}


void
displayBoundary( KnSpace* ks, C4CIteratorOnSurface* cp, bool loop )
{
  cout << "#######################################################" << endl;
  cout << "# displayBoundary" << endl
       << "# displays the points on the boundary as couples." << endl
       << "# x y" << endl;
  Proxy<C4CIteratorOnSurface> it( (C4CIteratorOnSurface*) cp->clone() );
  Proxy<C4CIteratorOnSurface> stop( (C4CIteratorOnSurface*) cp->clone() );
  Frame2D frame;
  frame.init( ks, 0, 1 );
  do
    {
      Kn_sid sbel = it->current();
      frame.setSurfelFrame( sbel, it->trackDir() );
      Vector2i p( frame.transformPoint( Vector2i( 0, 0 ) ) );
      cout << p.x() << " " << p.y() << endl;
      it->next();
    }
  while ( ! it->equals( *stop ) );
  if ( loop )
    {
      Kn_sid sbel = it->current();
      frame.setSurfelFrame( sbel, it->trackDir() );
      Vector2i p( frame.transformPoint( Vector2i( 0, 0 ) ) );
      cout << p.x() << " " << p.y() << endl;
    }
}

void testSurfelFrame( KnSpace* ks, C4CIteratorOnSurface* cp )
{
  cout << "#######################################################" << endl;
  cout << "# testSurfelFrame" << endl
       << "# for each surfel, extract tangent and transforms it in" << endl
       << "# the standard frame." << endl
       << "# tgt.x tgt.y tgt2d.x tgt2d.y theta T[tgt2d].x T[tgt2d].y T[theta]"
       << "# T[(0,0)].x T[(0,0)].y T[c_n()].x T[c_n()].y T[cp_n()].x T[cp_n()].y" 
       << endl;
  Frame2D frame;
  frame.init( ks, 0, 1 );
  Kn_sid sbel = cp->current();
  Kn_sid bel = sbel;
  do
    {
      Proxy<C4CIterator> cpfwd( cp->clone() );
      Proxy<C4CIterator> cpbwd( cp->clone() );
      C4CSegment segment =
	C4CGeometry::maximalFrontTangent( *cpfwd, *cpbwd, 0 );
      Vector2i tgt = segment.getTangent();
      Vector2D tgt2d( tgt.x(), tgt.y() );
      float theta = 
	(float) atan2( (double) tgt2d.ro( 1 ), (double) tgt2d.ro( 0 )  );
      frame.setSurfelFrame( bel, *( ks->sbegin_dirs( bel ) ) );
      cout << tgt.x() << " " << tgt.y()
	   << " " << tgt2d.x() << " " << tgt2d.y()
	   << " " << theta
	   << " " << frame.transform( tgt2d ).x()
	   << " " << frame.transform( tgt2d ).y()
	   << " " << frame.angleToX( theta )
	   << " " << frame.transformPoint( Vector2D( 0, 0 ) ).x()
	   << " " << frame.transformPoint( Vector2D( 0, 0 ) ).y()
	   << " " << frame.transformPoint( segment.c_n() ).x()
	   << " " << frame.transformPoint( segment.c_n() ).y()
	   << " " << frame.transformPoint( segment.cp_n() ).x()
	   << " " << frame.transformPoint( segment.cp_n() ).y()
	   << endl;
      cp->next();
      bel = cp->current();
    }
  while ( bel != sbel );
}


void
testSignal( KnSpace* ks, C4CIteratorOnSurface* cp, double dh, uint nb_bels )
{
  uint n = (uint) ceil( 0.5 / pow( dh, 4.0/3.0 ) ); 
  cout << "#######################################################" << endl;
  cout << "# testSignal n=" << n << endl
       << "# x y hx hy dx dy i" << endl;

  Signal<double> X( nb_bels, 0, true, 0.0 );
  Signal<double> Y( nb_bels, 0, true, 0.0 );
 
  Proxy<C4CIteratorOnSurface> it( (C4CIteratorOnSurface*) cp->clone() );
  Proxy<C4CIteratorOnSurface> stop( (C4CIteratorOnSurface*) cp->clone() );
  Frame2D frame;
  frame.init( ks, 0, 1 );
  int i = 0;
  do
    {
      Kn_sid sbel = it->current();
      frame.setSurfelFrame( sbel, it->trackDir() );
      Vector2i p( frame.transformPoint( Vector2i( 0, 0 ) ) );
      X[ i ] = p.x();
      Y[ i ] = p.y();
      // cout << p.x() << " " << p.y() << endl;
      it->next();
      ++i;
    }
  while ( ! it->equals( *stop ) );

  //       Signal<double> Hn( Signal<double>::H2n( i ) );
  //       cout << "H_" << 2*i << " = ";
  //       for ( int j = -10; j <= 10; j++ )
  // 	cout << " " << Hn[ j ];
  //       cout << endl;
  Clock::startClock();

  // Second technique
  Signal<double> HX( X );
  Signal<double> HY( Y );
  Signal<double> DX( X * Signal<double>::Delta() );
  Signal<double> DY( Y * Signal<double>::Delta() );
  Signal<double> G( Signal<double>::H2() );
  G.multiply( 0.25 );
  while ( n > 0 ) {
    HX = HX * G;
    HY = HY * G;
    DX = DX * G;
    DY = DY * G;
    --n;
  }
  double s = 1.0;

  // First technique.
//   Signal<double> Dn( Signal<double>::D2n( n ) );
//   Signal<double> HX( X * Signal<double>::H2n( n ) );
//   Signal<double> HY( Y * Signal<double>::H2n( n ) );
//   Signal<double> DX( X * Signal<double>::D2n( n ) );
//   Signal<double> DY( Y * Signal<double>::D2n( n ) );
//   cout << "# D_" << 2*n << " = ";
//   for ( int j = -10; j <= 10; j++ )
//     cout << " " << Dn[ j ];
//   cout << endl;
//  double s = 1.0 / pow( 4.0, (int) n );


  long t = Clock::stopClock();
  cerr << "# Convolution time = " << t << " (ms)." << endl;

  for ( uint i = 0; i < nb_bels; i++ )
    cout << X[ i ] << " " << Y[ i ] << " "
	 << HX[ i ] * s << " " << HY[ i ] * s << " "
	 << DX[ i ] * s << " " << DY[ i ] * s << " "
	 << i
	 << endl;
}



void
testDistanceTransform( KnSpace* ks, KnCharSet voxset, Kn_sid bel )
{
  KnRCellSet bdry = KnShapes::strackBoundary( *ks, voxset, bel );
  KnRUCellVector<int>* dmap 
    = KnShapes::computeBdryCityBlockDistanceMap( *ks, bdry );
  KnSpaceScanner scanner( *ks, 
			  ks->ufirstCell( ks->dim() ),
			  ks->ulastCell( ks->dim() ) );
  Kn_uid p = scanner.begin();
  do 
    { // ... whatever [p] is the current cell
      cout << (*dmap)[ p ] << " ";
    }
  while ( scanner.increment( p ) ); 
  cout << endl;
  int m = * ( min_element( dmap->begin(), dmap->end() ) );
  int M = * ( max_element( dmap->begin(), dmap->end() ) );
  cout << "min=" << m << " max=" << M << endl;
  delete dmap;

  Kn_uid in, out;
  vector<double> p_in( 3 );
  p_in[0] = 0.0;
  p_in[1] = 0.5;
  p_in[2] = 0.15;
  vector<double> p_out( 4 );
  p_out[0] = 0.0;
  p_out[1] = 0.3;
  p_out[2] = 0.1;
  p_out[3] = 0.05;
  KnCharSet nobj
    = KnShapes::noisifyObject( *ks, voxset, bdry, p_in, p_out, in, out );
  
  KnSpaceScanner scan2( *ks, 
			ks->ufirstCell( ks->dim() ),
			ks->ulastCell( ks->dim() ) );

  Kn_uid last_y, last_x;
  p = scan2.begin();
  for ( last_y = scan2.last( p, 1 );
	p <= last_y; 
	p += scan2.gotonext( 1 ) )
    {
      for ( last_x = scan2.last( p, 0 ); 
	    p <= last_x; 
	    p++ ) // NB: 'scan.gotonext( 0 )' == 1;
	{ //... whatever
	  if ( voxset[ p ] == nobj[ p ] )
	    cout << ( nobj[ p ] ? 'A' : '.' );
	  else
	    cout << ( nobj[ p ] ? '1' : '0' );
	}
      cout << endl;
    }
}


void
noisify( KnSpace* ks, KnCharSet & voxset, Kn_sid bel, 
	 Kn_uid & in, Kn_uid & out )
{
  KnRCellSet bdry = bel != 0 
    ? KnShapes::strackBoundary( *ks, voxset, bel )
    : KnShapes::smakeBoundary( *ks, voxset );
  
  KnRUCellVector<int>* dmap 
    = KnShapes::computeBdryCityBlockDistanceMap( *ks, bdry );
  int m = * ( min_element( dmap->begin(), dmap->end() ) );
  int M = * ( max_element( dmap->begin(), dmap->end() ) );

  string law = args.getOption( "-noisify" )->getValue( 0 );
  double a = (double) args.getOption( "-noisify" )->getDoubleValue( 1 );
  double b = (double) args.getOption( "-noisify" )->getDoubleValue( 2 );
  
  vector<double> p_in( M + 1 );
  vector<double> p_out( -m + 1 );
  if ( law == "POWER" )
    {
      p_in[ 0 ] = 1.0;
      p_out[ 0 ] = 1.0;
      double bb = a/b;
      for ( uint k = 1; k < p_in.size(); ++k )
	{
	  p_in[ k ] = bb;
	  bb /= b;
	}
      bb = a/b;
      for ( uint k = 1; k < p_out.size(); ++k )
	{
	  p_out[ k ] = bb;
	  bb /= b;
	}
    }

  voxset = KnShapes::noisifyObject( *ks, voxset, bdry, p_in, p_out, in, out );
}


Kn_sid
findInnerObject( KnSpace* ks, KnCharSet voxset, Kn_uid in,
		 KnCharSet & main_inner_comp )
{
  main_inner_comp 
    = KnShapes::uexpandSeedWithinBounds( *ks, in, ~voxset );
  // cerr << main_inner_comp << endl;
  
  return KnShapes::sfindFurthestBel( *ks, in, main_inner_comp );
}

  
Kn_sid
findOuterObject( KnSpace* ks, KnCharSet voxset, Kn_uid out, 
		 KnCharSet & main_outer_comp, Kn_uid in  )
{
  main_outer_comp 
    = KnShapes::uexpandSeedWithinBounds( *ks, out, voxset );
  // cerr << main_outer_comp << endl;
  
  return ks->sopp( KnShapes::sfindClosestBel( *ks, in, ~main_outer_comp ) );
}

void
textview( KnSpace* ks, KnCharSet voxset, char out, char in )
{
  KnSpaceScanner scan2( *ks, 
			ks->ufirstCell( ks->dim() ),
			ks->ulastCell( ks->dim() ) );

  Kn_uid last_y, last_x;
  Kn_uid p = scan2.begin();
  for ( last_y = scan2.last( p, 1 );
	p <= last_y; 
	p += scan2.gotonext( 1 ) )
    {
      for ( last_x = scan2.last( p, 0 ); 
	    p <= last_x; 
	    p++ ) // NB: 'scan.gotonext( 0 )' == 1;
	{ //... whatever
	  cout << ( voxset[ p ] ? in : out );
	}
      cout << endl;
    }
}


void
outputPGM( KnSpace* ks, KnCharSet voxset )
{
  cout << "P5" << endl
       << "# CREATOR: test_TangentialCover " 
       << "(lachaud@labri.fr)" << endl;
  cout << ks->size( 0 ) << " " << ks->size( 1 ) << endl
       << "255" << endl;
  
  KnSpaceScanner scan2( *ks, 
			ks->ufirstCell( ks->dim() ),
			ks->ulastCell( ks->dim() ) );

  Kn_uid last_y, last_x;
  Kn_uid p = scan2.begin();
  for ( last_y = scan2.last( p, 1 );
	p <= last_y; 
	p += scan2.gotonext( 1 ) )
    {
      for ( last_x = scan2.last( p, 0 ); 
	    p <= last_x; 
	    p++ ) // NB: 'scan.gotonext( 0 )' == 1;
	{ //... whatever
	  cout << ( voxset[ p ] ? (char) 0 : (char) 255 );
	}
    }
  cout << endl;
}


void viewMaximalSegments( KnSpace* ks, C4CIterator & cp )
{
  C4CTangentialCover tcover;
  buildTangentialCover( tcover, cp, 0 );
  uint si = 0; // index of surfel
  Mathutils::ModuloComputer mcsu( tcover.nbSurfels() );
  Frame2D frame;
  frame.init( ks, 0, 1 );

  cout << "# viewMaximalSegments" << endl
       << "# provides information on maximal segment in the global frame." 
       << endl
       << "# msi a b n nbU nbL xU yU xU' yU' xL yL xL' yL'" << endl;
  uint nbULU = 0;
  for (uint msi = 0; msi < tcover.nbMaximalSegments(); ++msi )
    {
      const C4CTangentialCover::MaximalSegment & ms 
	= tcover.getMaximalSegment( msi );
      while ( si != ms.front_surfel_idx ) 
	{
	  cp.next();
	  si = mcsu.next( si );
	}
      C4CIteratorOnSurface* cpi = dynamic_cast<C4CIteratorOnSurface*>( &cp );
      Kn_sid bel = cpi->current();
      uint t = cpi->trackDir();
      frame.setSurfelFrame( bel, t );
      const C4CSegment & segment = ms.dss;
      Vector2i p = frame.transform( segment.getTangent() );
      Vector2i u = frame.transformPoint( segment.u() );
      Vector2i up = frame.transformPoint( segment.up() );
      Vector2i l = frame.transformPoint( segment.l() );
      Vector2i lp = frame.transformPoint( segment.lp() );
      Vector2i c = frame.transformPoint( segment.c_n() );
      Vector2i cp = frame.transformPoint( segment.cp_n() );
      cout << msi << " " << p.y() << " " << p.x()
	   << " " << "TODO"
	   << " " << segment.nbU()
	   << " " << segment.nbL()
	   << " " << u.x() << " " << u.y()
	   << " " << up.x() << " " << up.y()
	   << " " << l.x() << " " << l.y()
	   << " " << lp.x() << " " << lp.y()
	   << " " << c.x() << " " << c.y()
	   << " " << cp.x() << " " << cp.y()
	   << endl;
      nbULU += segment.isULU() ? 1 : 0;
    }
  cout << "# NBMS NBULU NBLUL" << endl
       << "# " << tcover.nbMaximalSegments() << " " << nbULU 
       << " " << tcover.nbMaximalSegments() - nbULU << endl;
}




void
displayBoundaryWord( KnSpace* ks, C4CIteratorOnSurface* cp )
{
  cout << "#######################################################" << endl;
  cout << "# displayBoundaryWord" << endl
       << "# displays the contour as a word with freeman codes." << endl
       << "# x y" << endl;
  Proxy<C4CIteratorOnSurface> it( (C4CIteratorOnSurface*) cp->clone() );
  Proxy<C4CIteratorOnSurface> stop( (C4CIteratorOnSurface*) cp->clone() );
  Frame2D frame;
  frame.init( ks, 0, 1 );
  int nb = 0;
  do
    {
      Kn_sid sbel = it->current();
      frame.setSurfelFrame( sbel, it->trackDir() );
      Vector2i p1( frame.transformPoint( Vector2i( 0, 0 ) ) );
      Vector2i p2( frame.transformPoint( Vector2i( 1, 0 ) ) );
      int dx = p2.x() - p1.x();
      int dy = p2.y() - p1.y();
      char code = '0' + (char) ( dx != 0 ? (1 - dx) : (2 - dy) );
      cout << code; //<< "(" << dx << "," << dy << ")";
      it->next();
    }
  while ( ! it->equals( *stop ) );
  cout << endl << "# nb=" << nb << endl;

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
  // -------------------------------------------------------------------------
  // Prepare arguments.
  StandardArguments::addDigitalArgs( args, 2, false, false );
  ShapeHelper::addStarShapedArgs( args );

//   args.addOption( "-step", "-step <h>: the discretization step or resolution, the closer to 0, the finer.", "1.0" );
//   args.addOption( "-circle", "-circle <R>: the test shape is a circle of radius R", "10.0" );
//   args.addOption( "-flower", "-flower <R> <r> <k> <phi>: the test shape is a flower with k extremeties with mean radius R and variability of radius r, phi is the phase of the flower", "3", "10.0", "5.0", "0.0" );
//   args.addOption( "-accflower", "-accflower <R> <r> <k>: the test shape is a phase accelerating flower with k extremeties with mean radius R and variability of radius r", "4.0", "2.0", "4" );
//   args.addOption( "-rsquare", "-rsquare <R> <r>: the test shape is rounded square of big radius R and small corner radius r", "4.0", "1.0" );
//   args.addBooleanOption( "-fuzzy_uniform", "-fuzzy_uniform: possible uniform random shift of the shape center." );


  args.addOption( "-file", "-file <fname>: read the file <fname> to define a set of interior pixels", "" );
  args.addOption( "-noisify", "-noisify <type> <a> <b>: noisify the shape within study according to the distance k to the boundary ; with the law <type>=POWER such that P_change(k)=a/b^k", "POWER", "1.0", "2.0" );

  args.addBooleanOption( "-tSF", "-tSF: testSurfelFrame." );
  args.addBooleanOption( "-tGG", "-tGG: testGlobalGeometry." );
  args.addBooleanOption( "-tGGM", "-tGGM: testGlobalGeometryMinimizer." );
  args.addBooleanOption( "-tSCP", "-tSCP: testSizeCommonParts." );
  args.addBooleanOption( "-tNMS", "-tNMS: testNbMaximalSegments." );
  args.addBooleanOption( "-tNMSPS", "-tNMSPS: testNbMSPerSurfel." );
  args.addBooleanOption( "-dCPG", "-dCPG: displayCommonPartGeometry. Allows to display common parts, result of global geometry optimization, and estimated curvature. See file global_geom.txt" );
  args.addOption( "-dB", "-dB <loop>: display boundary. loop=1 : rewrite the first point.", "1" );
  args.addBooleanOption( "-tSignal", "-tSignal: testSignal. Tangent estimator of Fourey-Malgouyres based on discrete convolutions." );
  args.addBooleanOption( "-tDT", "-tDT: testDistanceTransform. Computes the city-block distance transform to an object boundary." );
  args.addBooleanOption( "-textview", "-textview: display the shape on the standard output." );
  args.addBooleanOption( "-pgm", "-pgm: output shape in pgm format on the standard output stream." );
  args.addBooleanOption( "-inner_object", "-inner_object: extract inner contour (useful with noisify)." );
  args.addBooleanOption( "-outer_object", "-outer_object: extract outer contour (useful with noisify)." );
  args.addBooleanOption( "-vMS", "-vMS: viewMaximalSegment: provides lots of information on each MS in the global frame." );
  args.addBooleanOption( "-dTCG", "-dTCG: displayTangentialCoverGeometry: provides information on each MS related to global curvature estimator." );
  args.addBooleanOption( "-dBW", "-dBW: display boundary word with freeman codes." );
  args.addOption( "-minimizer", "-minimizer <type> <a>: choose the minimizer for global geometry computation, <type>=STD mix gradient/optimal, <type>=GD gradient descent, <type>=AGD adaptive step gradient descent, <type>=RLX Relaxation.", "STD", "0.5" );
  args.addOption( "-eps", "-eps <max> <sum>: specifies max and sum epsilon to stop optimization in global geometry computation.", "0.0000001", "-1.0" );

  if ( ( argc <= 1 ) 
       || ! args.readArguments( argc, argv ) ) 
    {
      cerr << args.usage( "test_TangentialCover", 
			  "Tests the tangential cover of a digital curve.",
			  "" )
	   << endl;
      // -x -y -circle -flower -accflower -rsquare -step -fuzzy_uniform -center -file -minimizer -eps -tSF -tGG -tSCP -tNMS -tNSMPS -dCPG -tSignal -tDT -noisify -textview -pgm -inner_object -outer_object -vMS -dTCG -dB -dBW      
      return 1;
    }

  // -------------------------------------------------------------------------
  // Display command line.
//   cout << "#";
//   for ( int i = 0; i < argc; ++i )
//     cout << " " << argv[ i ];
//   cout << endl;

  // -------------------------------------------------------------------------
  // Build space.
  uint d = StandardArguments::dim( args );
  if ( d != 2 )
    {
      cerr << "Dimension is 2." << endl;
      return 2;
    }
  Kn_size sizes[ d ];
  StandardArguments::fillSizes( args, sizes );
  KnSpace ks( 2, sizes );
  Kn_uid vcenter = ks.uspel( ks.ukcode( sizes ) ); // center
  Vector xcenter = ks.ucentroid( vcenter );
  
  // -------------------------------------------------------------------------
  // Get some parameters.
  double dh = (double) args.getOption( "-step" )->getDoubleValue( 0 );
  GridEmbedder embedder;
  embedder.init( &ks );
  embedder.setCenter( xcenter );
  embedder.setScale( dh );

  // -------------------------------------------------------------------------
  // Make shape.
  uint nb_bels = 0;
  Kn_sid bel = 0;
  KnCharSet voxset = KnCharSet::create( ks, ks.dim(), false, 0 );
  KnCharSet aux_voxset( voxset );
  
  StarShaped* shape = shapeFromArgs( ks, embedder, voxset, bel, nb_bels );
  
  // -------------------------------------------------------------------------
  // Change shape.

  Kn_uid in, out;
  if ( args.check( "-noisify" ) )
    {
      noisify( &ks, voxset, bel, in, out );
      bel = findInnerObject( &ks, voxset, in, aux_voxset );
    }      
  bool interior = true;
  if ( args.check( "-inner_object" ) )
    {
      // NB: change contour and object.
      bel = findInnerObject( &ks, voxset, in, aux_voxset );
      voxset = aux_voxset;
      interior = true;
    }
  if ( args.check( "-outer_object" ) )
    {
      // NB: change contour and object.
      bel = findOuterObject( &ks, voxset, out, aux_voxset, in );
      voxset = aux_voxset;
      interior = false;
    }

  cerr << "# Bel:";
  ks.displayKn_sid( bel, cerr );
  cerr << endl;
  
  if ( args.check( "-textview" ) )
    textview( &ks, voxset, '0', '1' );
  if ( args.check( "-pgm" ) )
    outputPGM( &ks, voxset );

  // -------------------------------------------------------------------------
  // Choose Minimizer.
  LinearMinimizer* lm = 0;
  if ( args.getOption( "-minimizer" )->getValue( 0 ) == "GD" )
    lm = new LinearMinimizerByGradientDescent
      ( args.getOption( "-minimizer" )->getDoubleValue( 1 ) );
  else if ( args.getOption( "-minimizer" )->getValue( 0 ) == "AGD" )
    lm = new LinearMinimizerByAdaptiveStepGradientDescent
      ( args.getOption( "-minimizer" )->getDoubleValue( 1 ) );
  else if ( args.getOption( "-minimizer" )->getValue( 0 ) == "RLX" )
    lm = new LinearMinimizerByRelaxation;
  else
    lm = new LinearMinimizer;

  if ( lm != 0 ) 
    {
      cout << "# " << *lm << endl
	   << "# max_eps=" <<  args.getOption( "-eps" )->getDoubleValue( 0 )
	   << " sum_eps=" <<  args.getOption( "-eps" )->getDoubleValue( 1 )
	   << endl;
    }

  // -------------------------------------------------------------------------
  // Create boundary.
  BelAdjacency badj( ks, interior );
  ObjectBoundary bdry( badj, voxset );
  uint dir = *( ks.sbegin_dirs( bel ) );
  C4CIteratorOnSurface* cp = bdry.newC4CIterator( bel, dir );

  bool open;
  uint nb = C4CIterator::size( *cp, open );
  cerr << "# Size boundary Nb=" << nb << " open=" << open << endl;
  
  // testModuloComputer();
  if ( args.check( "-tSF" ) ) 
    testSurfelFrame( &ks, cp );
  if ( args.check( "-tSCP" ) ) 
    testSizeCommonParts( *cp );
  if ( args.check( "-tNMS" ) ) 
    testNbMaximalSegments( *cp );
  if ( args.check( "-tNMSPS" ) ) 
    testNbMSPerSurfel( *cp );
  if ( args.check( "-tGG" ) ) 
    testGlobalGeometry( &ks, *cp, shape, embedder, dh, lm );
  if ( args.check( "-tGGM" ) ) 
    testGlobalGeometryMinimizer( &ks, *cp, shape, embedder, dh, lm );
  if ( args.check( "-dCPG" ) ) 
    displayCommonPartGeometry( &ks, *cp, lm );
  if ( args.check( "-tSignal" ) ) 
    testSignal( &ks, cp, dh, nb_bels );
  if ( args.check( "-tDT" ) ) 
    testDistanceTransform( &ks, voxset, bel );

  
  if ( args.check( "-dB" ) ) 
    displayBoundary( &ks, cp, 
		     (bool) args.getOption( "-dB" )->getIntValue( 0 ) );

  if ( args.check( "-dBW" ) ) 
    displayBoundaryWord( &ks, cp );
  
  if ( args.check( "-vMS" ) ) 
    viewMaximalSegments( &ks, *cp );

  // testTangentialCoverGeometry( &ks, *cp );
//   computeCurvature( &ks, 
// 		    args.getOption( "-file" )->getValue( 0 ), 
// 		    embedder, 
// 		    dh );

  // displayCommonPartGeometry( &ks, *cp );

  if ( args.check( "-dTCG" ) ) 
    displayTangentialCoverGeometry( &ks, *cp );

  // displayShape( &ks, *cp, shape, embedder, dh );
  
  delete cp;
  delete shape;
  return 0;
}
