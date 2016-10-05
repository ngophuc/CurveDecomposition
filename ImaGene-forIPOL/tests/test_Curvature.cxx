///////////////////////////////////////////////////////////////////////////////
// Test module for digitalnD and Kn_uid, Kn_sid
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cmath>
#include <cstdlib>

//#include "LinAlg/LinAlg2D/Vector2D.hpp"
#include "ImaGene/base/Vector.h"
#include "ImaGene/base/Arguments.h"
#include "ImaGene/base/StandardArguments.h"
#include "ImaGene/timetools/Clock.h"
#include "ImaGene/dgeometry2d/C4CSegment.h"
#include "ImaGene/dgeometry2d/C4CIterator.h"
#include "ImaGene/dgeometry2d/C4CGeometry.h"
#include "ImaGene/dgeometry2d/EuclideanGeometry.h"
#include "ImaGene/digitalnD/KnCharSet.h"
#include "ImaGene/digitalnD/KnRCellSet.h"
#include "ImaGene/digitalnD/BelAdjacency.h"
#include "ImaGene/digitalnD/ObjectBoundary.h"
#include "ImaGene/digitalnD/KnShapes.h"


using namespace std;
using namespace ImaGene;

static Arguments args;


void computeCurvatureByWindowedAngle( KnSpace & ks, 
				      KnCharSet voxset,
				      KnRCellSet dsurf,
				      Kn_uid center,
				      float radius )
{
  KnRCellSet::cell_iterator p = dsurf.begin();
  Kn_sid bel = *p;

  float cx = ks.ucentroid( center, 0 );
  float cy = ks.ucentroid( center, 1 );
  
  // Get tangent plane.
  BelAdjacency badj( ks, true );
  ObjectBoundary bdry( badj, voxset );
  uint track_dir = *( ks.sbegin_dirs( bel ) );
  C4CIteratorOnSurface* cp = bdry.newC4CIterator( bel, track_dir );

  float esp = 0.0;
  float esp2 = 0.0;
  float integrated_esp = 0.0;
  uint nb = 0;
  float l = 0.0;
  Kn_sid starting_bel = bel; 
  do 
    {
      // Computes sym tgt centered.
      C4CIteratorOnSurface* cpback = 
	static_cast<C4CIteratorOnSurface*>( cp->clone() );
      C4CIteratorOnSurface* cpfront = 
	static_cast<C4CIteratorOnSurface*>( cp->clone() );
      C4CSegment symtgt_c;
      symtgt_c.init();
      uint nb_c = 0;
      while ( symtgt_c.extendsUnsecure( cpback->previous(), cpfront->next() ) )
	++nb_c;
      //cpback->next();
      //cpfront->previous();
      float ds_c = symtgt_c.getTangent().averagedLength();
      Vector2D C( ks.scentroid( cp->current() ) );
      Vector2D R( ks.scentroid( cpfront->current() ) );
      Vector2D L( ks.scentroid( cpback->current() ) );

      // Computes sym tgt left
      C4CIterator* leftfront = cpback->clone();
      C4CIterator* leftback = cpback->clone();
      Vector2i tgt_l = 
	C4CGeometry::symmetricTangent( *leftfront, *leftback, 0 )
	.getLine().tangent();
      float ds_l = tgt_l.averagedLength();
      delete leftback;
      delete leftfront;

      // Computes sym tgt right
      C4CIterator* rightfront = cpfront->clone();
      C4CIterator* rightback = cpfront->clone();
      Vector2i tgt_r = 
	C4CGeometry::symmetricTangent( *rightfront, *rightback, 0 )
	.getLine().tangent();
      float ds_r = tgt_r.averagedLength();
      delete rightback;
      delete rightfront;

      while ( cpback->current() != cp->current() )
	{
	  C4CSegment::slidesForward( cpback->next(), tgt_l );
	}
      while ( cpfront->current() != cp->current() )
	{
	  C4CSegment::slidesBackward( cpfront->previous(), tgt_r );
	}
      

      delete cpfront;
      delete cpback;

      // Computes geodesic distance.
//       float hh = 
// 	sqrt( ( C.x() - L.x() ) * ( C.x() - L.x() )
// 	      + ( C.y() - L.y() ) * ( C.y() - L.y() ) )
// 	+ sqrt( ( C.x() - R.x() ) * ( C.x() - R.x() )
// 		+ ( C.y() - R.y() ) * ( C.y() - R.y() ) );
      float hh = ( nb_c + 1 )* ( ds_l + 2*ds_c + ds_r ) / 2;

      float normc_p = symtgt_c.getTangent().norm();
      float xc_p = symtgt_c.getTangent().x() / normc_p;
      float yc_p = symtgt_c.getTangent().y() / normc_p;

      float norml_p = tgt_l.norm();
      float xl_p = tgt_l.x() / norml_p;
      float yl_p = tgt_l.y() / norml_p;

      float normr_p = tgt_r.norm();
      float xr_p = tgt_r.x() / normr_p;
      float yr_p = tgt_r.y() / normr_p;

      float xc_pp = ( xr_p - xl_p ) / ( hh );
      float yc_pp = ( yr_p - yl_p ) / ( hh );

      float curv = fabs( xc_p * yc_pp - xc_pp * yc_p );

      float bx = ks.scentroid( bel, 0 );
      float by = ks.scentroid( bel, 1 );
      if ( args.check( "-per_surfel" ) )
	cout << "K_av_symtgt( " << bx - cx << " , " << by -cy << " )= " 
	     << curv << endl;

      esp += curv;
      esp2 += curv*curv;
      integrated_esp += curv; // * tline.tangent().averagedLength();
      // Go to next one.
      cp->next();
      bel = cp->current();
      ++nb;
    }
  while ( bel != starting_bel );

  delete cp;
  
  cout << "Exp= " << 1/radius << "  Moy= " << esp/nb 
       << "  IntK= " << integrated_esp / l
       << "  var= " << sqrt( esp2/nb - (esp/nb)*(esp/nb) )
       << "  nb= " << nb << "  L/2piR= " << l / ( 2 * M_PI * radius ) 
       << endl;
  
}



void computeCurvatureByAngleVariation( KnSpace & ks, 
				       KnCharSet voxset,
				       KnRCellSet dsurf,
				       Kn_uid center,
				       float radius )
{
  KnRCellSet::cell_iterator p = dsurf.begin();
  Kn_sid bel = *p;

  float cx = ks.ucentroid( center, 0 );
  float cy = ks.ucentroid( center, 1 );
  
  // Get tangent plane.
  BelAdjacency badj( ks, true );
  ObjectBoundary bdry( badj, voxset );
  uint track_dir = *( ks.sbegin_dirs( bel ) );
  C4CIteratorOnSurface* cp = bdry.newC4CIterator( bel, track_dir );

  float esp = 0.0;
  float esp2 = 0.0;
  float integrated_esp = 0.0;
  uint nb = 0;
  float l = 0.0;
  Kn_sid starting_bel = bel; 
  do 
    {
      C4CIterator* cpback = cp->clone();
      C4CIterator* cpfront = cp->clone();
      DLine tline = C4CGeometry::symmetricTangent( *cpfront, *cpback, 0 )
	.getLine();
      l += tline.tangent().averagedLength();
      delete cpfront;
      delete cpback;
      C4CIterator* cpcopy = cp->clone();
      float curv = C4CGeometry::curvatureBySymAngleVariation( *cpcopy, 0 );
      delete cpcopy;
//       C4CIteratorOnSurface* cfront = cp->clone();
//       C4CIteratorOnSurface* cback = cp->clone();
//       float curv = C4CGeometry::curvatureByAngleVariation( *cfront, *cback, 0 );
//       delete cback;
//       delete cfront;
      float bx = ks.scentroid( bel, 0 );
      float by = ks.scentroid( bel, 1 );
      if ( args.check( "-per_surfel" ) )
	cout << "K_av_symtgt( " << bx - cx << " , " << by -cy << " )= " 
	     << curv << endl;

      esp += curv;
      esp2 += curv*curv;
      integrated_esp += curv * tline.tangent().averagedLength();
      // Go to next one.
      cp->next();
      bel = cp->current();
      ++nb;
    }
  while ( bel != starting_bel );

  delete cp;
  
  cout << "Exp= " << 1/radius << "  Moy= " << esp/nb 
       << "  IntK= " << integrated_esp / l
       << "  var= " << sqrt( esp2/nb - (esp/nb)*(esp/nb) )
       << "  nb= " << nb << "  L/2piR= " << l / ( 2 * M_PI * radius ) 
       << endl;
  
}

void computeCurvatureByCircumscribedCircle( KnSpace & ks, 
					    KnCharSet voxset,
					    KnRCellSet dsurf,
					    Kn_uid center,
					    float radius )
{
  KnRCellSet::cell_iterator p = dsurf.begin();
  Kn_sid bel = *p;

  float cx = ks.ucentroid( center, 0 );
  float cy = ks.ucentroid( center, 1 );
  
  // Get tangent plane.
  BelAdjacency badj( ks, true );
  ObjectBoundary bdry( badj, voxset );
  uint track_dir = *( ks.sbegin_dirs( bel ) );
  C4CIteratorOnSurface* cp = bdry.newC4CIterator( bel, track_dir );

  float esp = 0.0;
  float esp2 = 0.0;
  float integrated_esp = 0.0;
  uint nb = 0;
  float l = 0.0;
  Kn_sid starting_bel = bel; 
  do 
    {
      C4CIteratorOnSurface* cpback = 
	static_cast<C4CIteratorOnSurface*>( cp->clone() );
      C4CIteratorOnSurface* cpfront = 
	static_cast<C4CIteratorOnSurface*>( cp->clone() );
      C4CSegment poshalftgt;
      poshalftgt.init();
      C4CSegment neghalftgt;
      neghalftgt.init();
      C4CGeometry::longestPositiveSegment( poshalftgt, *cpfront, 0 );
      C4CGeometry::longestNegativeSegment( neghalftgt, *cpback, 0 );
      Vector2D A( ks.scentroid( cp->current() ) );
      Vector2D B( ks.scentroid( cpfront->current() ) );
      Vector2D C( ks.scentroid( cpback->current() ) );
      delete cpfront;
      delete cpback;
      //DLine tline = C4CGeometry::symmetricTangent( *cpfront, *cpback, 0 );
      //l += tline.tangent().averagedLength();

      float curv = EuclideanGeometry::curvatureCircumscribedCircle
	( A.x(), A.y(), B.x(), B.y(), C.x(), C.y() );
      if ( curv < 0.0 ) curv = 0.0;
      float bx = ks.scentroid( bel, 0 );
      float by = ks.scentroid( bel, 1 );
      if ( args.check( "-per_surfel" ) )
	cout << "K_circum( " << A.x() - cx << " , " << A.y() - cy << " )= " 
	     << curv << endl;

      esp += curv;
      esp2 += curv*curv;
      integrated_esp += curv; // * tline.tangent().averagedLength();
      // Go to next one.
      cp->next();
      bel = cp->current();
      ++nb;
    }
  while ( bel != starting_bel );

  delete cp;
  
  cout << "Exp= " << 1/radius << "  Moy= " << esp/nb 
       << "  IntK= " << integrated_esp / l
       << "  var= " << sqrt( esp2/nb - (esp/nb)*(esp/nb) )
       << "  nb= " << nb << "  L/2piR= " << l / ( 2 * M_PI * radius ) 
       << endl;
  
}


int
main( int argc, char** argv ) 
{
  StandardArguments::addDigitalArgs( args, 2, true, false );

  args.addBooleanOption( "-per_surfel", 
			 "-per_surfel: asks to display the curvature for each surfel" );
  args.addBooleanOption( "-curv_symtgt_angle", 
			 "-curv_symtgt_angle: choose the angle variation of symmetric tangents as curvature estimator" );
  args.addBooleanOption( "-curv_circumcircle", "-curv_circumcircle: choose the curvature of the circumcircle to endpoints of halftangents as curvature estimator" );

  args.addBooleanOption( "-curv_windowed_angle", 
			 "-curv_windowed_angle: computes curvature by a variation angle defined over an adaptative window" );

  args.addOption( "-fuzzy", "-fuzzy <variance[=1.0]>: possible random shift of the circle center.", "1.0" );
  
  if ( ( argc <= 1 ) 
       || ! args.readArguments( argc, argv ) ) 
    {
      cout << args.usage( "test_Curvature", 
			  "Tests several discrete curvature estimators and output some statistics.",
			  "-h -d -x -y -ri -curv_symtgt_angle -curv_circumcircle -curved_windowed_angle -per_surfel -fuzzy" ) 
	   << endl;
      return 1;
    }

  float r1;

  // Build space.
  uint d = StandardArguments::dim( args );
  if ( d != 2 )
    {
      cerr << "Dimension is 2." << endl;
      return 2;
    }
  Kn_size sizes[ d ];
  StandardArguments::fillSizes( args, sizes );
  // cerr << "[ " << sizes[ 0 ] << " x " << sizes[ 1 ] << " ]" << endl;
  KnSpace ks( d, sizes );

  // center
  Kn_uid clow = ks.uspel( ks.ukcode( sizes ) );
  cout << "--- Space: " << ks << endl;
  float radius_first = args.getOption( "-ri" )->getFloatValue( 0 ); 
  float radius_last = args.getOption( "-ri" )->getFloatValue( 1 ); 
  float radius_incr = args.getOption( "-ri" )->getFloatValue( 2 ); 
  for ( r1 = radius_first; 
	r1 <= radius_last; 
	r1 += radius_incr )
    {
      cout << "--- creating sphere (r=" << r1 << ") -----" << endl;
      Clock::startClock();
      KnCharSet sph1 = KnShapes::umakeVolumicSphere( ks, clow, r1 );
      KnCharSet voxset = sph1;
      long ti2 = Clock::stopClock();
      cout << "in " << ti2 << " ms." << endl;
      cout << "r= " << r1 << " voxset   = " << voxset.nbElements() << " spels." << endl;
      
      //   uint belx[ 100 ];
      //   belx[ 0 ] = ( sizes[ 0 ] & ~1 ) + 2 * r + 2;
      //   for ( i = 1; i < ks.dim(); i++ )
      //     belx[ i ] = ( sizes[ i ] & ~1 ) + 1;
      //   Kn_sid bel = ks.negkcode( belx ); // starting bel for tracking
      cout << "--- boundary of sphere -----" << endl;
      Clock::startClock();
      KnRCellSet digsurf = KnShapes::smakeBoundary( ks, voxset );
      long ti3 = Clock::stopClock();
      cout << "in " << ti3 << " ms." << endl;
      cout << "digsurf   = " << digsurf.nbElements() << " surfels." << endl;
      
      if ( args.check( "-curv_symtgt_angle" ) )
 	computeCurvatureByAngleVariation( ks, voxset, digsurf, clow, r1 );
      if ( args.check( "-curv_circumcircle" ) )
 	computeCurvatureByCircumscribedCircle( ks, voxset, digsurf, clow, r1 );
      if ( args.check( "-curv_windowed_angle" ) )
 	computeCurvatureByWindowedAngle( ks, voxset, digsurf, clow, r1 );
    }
  return 0;
}
  
