///////////////////////////////////////////////////////////////////////////////
// Test module for digitalnD and Kn_uid, Kn_sid
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cmath>

#include "ImaGene/base/Arguments.h"
#include "ImaGene/base/StandardArguments.h"
#include "ImaGene/digitalnD/KnSpace.h"
#include "ImaGene/digitalnD/KnCharSet.h"
#include "ImaGene/digitalnD/KnShapes.h"
#include "ImaGene/digitalnD/K3Isosurface.h"
#include "ImaGene/digitalnD/K3Isosurface_k6_l18.h"
#include "ImaGene/digitalnD/BelAdjacency.h"
#include "ImaGene/digitalnD/ObjectBoundary.h"
#include "ImaGene/digitalnD/KnRCellVector.h"
#include "ImaGene/digitalnD/DigitalSurfaceGeometry.h"
#include "ImaGene/digitalnD/GeometryComputerByContour4.h"
#include "ImaGene/digitalnD/MeanCurvatureOnSurfaceComputerByDG.h"
#include "ImaGene/digitalnD/AreaOnSurfaceComputerByEuler.h"
#include "ImaGene/timetools/Clock.h"

using namespace std;
using namespace ImaGene;

static Arguments args;

/**
 * Add to [result] a SoNormal and a SoNormalBinding starting at index [index].
 * The normals are computed using a discrete normal estimator.
 */
void computeDiscreteNormals( const KnSpace & ks,
			     const KnCharSet & voxset,
			     const KnRCellSet & dig_surf )
{
  cout << "Computing normals..."; cout.flush();
  uint nb_vtx = dig_surf.nbElements();
  BelAdjacency badj( ks, true );
  ObjectBoundary bdry( badj, voxset );
  DigitalSurfaceGeometry ctxt_geometry;
  ctxt_geometry.setSurface( &bdry );
  GeometryComputerByContour4 concrete_geometry;
  ctxt_geometry.tangent_computer = &concrete_geometry;
  ctxt_geometry.normal_computer = &concrete_geometry;
  Clock::startClock();
  
  KnRCellSet::cell_iterator ip = dig_surf.begin();
  KnRCellSet::cell_iterator ip_end = dig_surf.end();
  Vector n( ks.dim() );
  while ( ip != ip_end )
    {
      // normals
      ctxt_geometry.normal( *ip, n );
      ++ip;
    }
  long ti5 = Clock::stopClock();
  cout << "in " << ti5 << " ms." << endl;
}

/**
 * Compute deviation of discrete normals compared to  perfect sphere normals.
 */
void computeDiscreteNormalsDeviation( const KnSpace & ks,
				      const KnCharSet & voxset,
				      const KnRCellSet & dig_surf )
{
  cout << "Computing normal deviation to perfect sphere normals ..."; 
  cout.flush();
  uint nb_vtx = dig_surf.nbElements();
  BelAdjacency badj( ks, true );
  ObjectBoundary bdry( badj, voxset );
  DigitalSurfaceGeometry ctxt_geometry;
  ctxt_geometry.setSurface( &bdry );
  GeometryComputerByContour4 concrete_geometry;
  ctxt_geometry.tangent_computer = &concrete_geometry;
  ctxt_geometry.normal_computer = &concrete_geometry;

  Clock::startClock();
  
  Kn_sid p;
  KnRCellSet::cell_iterator ip = dig_surf.begin();
  KnRCellSet::cell_iterator ip_end = dig_surf.end();

  Vector n_discrete( ks.dim() );
  Vector n_sphere( ks.dim() );
  Kn_size sizes[ ks.dim() ];
  uint i = 0;
  for ( i = 0; i < ks.dim(); i++ )
    sizes[ i ] = ks.size( i );
  Kn_uid center = ks.uspel( ks.ukcode( sizes ) );
  Vector vcenter = ks.ucentroid( center );
  double* angles = new double[ nb_vtx ];
  double medium = 0.0;
  i = 0;
  uint j;
  while ( ip != ip_end )
    {
      // normals
      p = *ip;
      ctxt_geometry.normal( p, n_discrete );
      ks.scentroid( p, n_sphere );
      n_sphere -= vcenter;
      double norm = 0.0;
      for ( j = 0; j < ks.dim(); j++ )
	norm += n_sphere.ro( j ) * n_sphere.ro( j );
      n_sphere /= sqrt( norm );
      // scalar product to mesure angle.
      norm = 0.0;
      for ( j = 0; j < ks.dim(); j++ )
	norm += - n_sphere.ro( j ) * n_discrete.ro( j );
      if ( norm >= 1.0 ) norm = 1.0;
      medium += angles[ i++ ] = acos( norm );
      //cout << " " << angles[ i - 1 ]* 180 / M_PI ;
      ++ip;
    }
  long ti5 = Clock::stopClock();
  cout << "in " << ti5 << " ms." << endl;
  medium /= nb_vtx;
  double acc = 0.0;
  for ( i = 0; i < nb_vtx; i++ )
    acc += ( angles[ i ] - medium ) * ( angles[ i ] - medium );
  double acc2 = 0.0;
  for ( i = 0; i < nb_vtx; i++ )
    acc2 += angles[ i ] * angles[ i ];
  cout << "Medium angle deviation = " << medium * 180 / M_PI 
       << "  standard dev = " << sqrt( acc / nb_vtx ) * 180 / M_PI 
       << " standard dev from 0 = "
       << sqrt( acc2 / nb_vtx ) * 180 / M_PI << endl;
  delete[] angles;
}

/**
 * Add to [result] a SoBaseColor and Normal and a SoMaterialBinding starting at index [index].
 * The normals are computed using a discrete normal estimator.
 */
float computeAreaMeasure( const KnSpace & ks,
			  const KnCharSet & voxset,
			  const KnRCellSet & dig_surf,
			  bool euler_integration )
{
  cout << "Computing areas "; cout.flush();
  uint nb_vtx = dig_surf.nbElements();
  BelAdjacency badj( ks, true );
  ObjectBoundary bdry( badj, voxset );
  DigitalSurfaceGeometry ctxt_geometry;
  ctxt_geometry.setSurface( &bdry );
  GeometryComputerByContour4 concrete_geometry;
  AreaOnSurfaceComputerByEuler area_by_euler;
  ctxt_geometry.tangent_computer = &concrete_geometry;
  ctxt_geometry.normal_computer = &concrete_geometry;
  if ( euler_integration )
    ctxt_geometry.area_computer = &area_by_euler;
  else
    ctxt_geometry.area_computer = &concrete_geometry;

  Clock::startClock();
  float warea = 0.0;
  Kn_sid p;
  KnRCellSet::cell_iterator ip = dig_surf.begin();
  KnRCellSet::cell_iterator ip_end = dig_surf.end();
  while ( ip != ip_end )
    {
      // normals
      p = *ip;
      // ks.displayKn_sid( p, cerr ); cerr << endl;
      float area = ctxt_geometry.area( p ); 
      warea += area;
      ++ip;
    }
  long ti5 = Clock::stopClock();
  cout << "in " << ti5 << " ms." << endl;
  return warea;
}


/**
 * Add to [result] a SoBaseColor and Normal and a SoMaterialBinding starting at index [index].
 * The normals are computed using a discrete normal estimator.
 */
void computeTangentPlanes( const KnSpace & ks,
			   const KnCharSet & voxset,
			   const KnRCellSet & dig_surf )
{
  cout << "Computing tangent planes "; cout.flush();
  uint nb_vtx = dig_surf.nbElements();
  BelAdjacency badj( ks, true );
  ObjectBoundary bdry( badj, voxset );
  DigitalSurfaceGeometry ctxt_geometry;
  ctxt_geometry.setSurface( &bdry );
  GeometryComputerByContour4 concrete_geometry;
  AreaOnSurfaceComputerByEuler area_by_euler;
  ctxt_geometry.tangent_computer = &concrete_geometry;
  ctxt_geometry.normal_computer = &concrete_geometry;
  ctxt_geometry.area_computer = &concrete_geometry;

  Clock::startClock();
  float x_low, x_up; 
  Kn_sid p;
  KnRCellSet::cell_iterator ip = dig_surf.begin();
  KnRCellSet::cell_iterator ip_end = dig_surf.end();
  Vector vn( ks.dim() );
  while ( ip != ip_end )
    {
      // normals
      p = *ip;
      //ctxt_geometry.normal( p, vn );
      ctxt_geometry.tangentPlane( p, vn, x_low, x_up ); 
      cout << " (" << x_low << "," << x_up 
	   << "," << ks.sdirect( p, ks.sorthDir( p ) ) << ")";
      ++ip;
    }
  cout << endl;
  long ti5 = Clock::stopClock();
  cout << "in " << ti5 << " ms." << endl;
}

/**
 * Add to [result] a SoBaseColor and Normal and a SoMaterialBinding starting at index [index].
 */
void computeMCurvatureMeasure( const KnSpace & ks,
			       const KnCharSet & voxset,
			       const KnRCellSet & dig_surf,
			       float sigma = 1.0 )
{
  cout << "Computing mean curvatures "; cout.flush();
  uint nb_vtx = dig_surf.nbElements();
  BelAdjacency badj( ks, true );
  ObjectBoundary bdry( badj, voxset );
  DigitalSurfaceGeometry ctxt_geometry;
  ctxt_geometry.setSurface( &bdry );
  GeometryComputerByContour4 concrete_geometry;
  MeanCurvatureOnSurfaceComputerByDG concrete_curv_geometry;
  concrete_curv_geometry.setParameters( sigma, 5000, 0.01, 0.01 );
  
  ctxt_geometry.tangent_computer = &concrete_geometry;
  ctxt_geometry.normal_computer = &concrete_geometry;
  ctxt_geometry.area_computer = &concrete_geometry;
  ctxt_geometry.length_computer = &concrete_geometry;
  ctxt_geometry.mean_curvature_computer = &concrete_curv_geometry;
  
  Clock::startClock();
  KnRCellSet::cell_iterator ip = dig_surf.begin();
  KnRCellSet::cell_iterator ip_end = dig_surf.end();
  while ( ip != ip_end )
    {
      // normals
      float mcurv = ctxt_geometry.meanCurvature( *ip ); 
      ++ip;
    }
  long ti5 = Clock::stopClock();
  cout << "in " << ti5 << " ms." << endl;
}

void displayArea( uint d, float r, float carea )
{
  if ( d == 2 )
    {
      cout << "S1: area=" << carea << " area/(2*PI*r)="
	   << carea / ( 2 * M_PI * r ) << endl;
      cout << "C1: area=" << carea << " area/(4*(2r)="
	   << carea / ( 8 * r ) << endl;
    }
  else if ( d == 3 )
    {
      cout << "S2: area=" << carea << " area/(4*PI*r^2)="
	   << carea / ( 4 * M_PI * r * r ) << endl;
      cout << "C2: area=" << carea << " area/(6*(4r^2)="
	   << carea / ( 24 * r * r ) << endl;
    }
  else if ( d == 4 )
    {
      cout << "S3: area=" << carea << " area/(2*PI^2*r^3)="
	   << carea / ( 2 * M_PI * M_PI * r * r * r ) << endl;
      cout << "C3: area=" << carea << " area/(8*(8r^3)="
	   << carea / ( 64 * r * r * r ) << endl;
    }
}

void outputArea( ostream & out, 
		 uint d, 
		 float r,
		 float carea )
{
  out << d << " " << r << " " << carea << " ";
  if ( d == 2 )
    out << carea / ( 2 * M_PI * r );
  else if ( d == 3 )
    out << carea / ( 4 * M_PI * r * r );
  else if ( d == 4 )
    out << carea / ( 6 * M_PI * r * r * r );
  out << endl;
}


int testKnMC( int argc, char** argv )
{
  StandardArguments::addDigitalArgs( args, 8, true, true );
  args.addBooleanOption( "-dnormals" ,"-dnormals: compute discrete normals." );
  args.addBooleanOption( "-planes" ,"-planes: compute planes." );
  args.addBooleanOption( "-area" ,"-planes: compute area." );
  args.addBooleanOption( "-area_euler" ,"-planes: compute area (Euler method)." );
  args.addBooleanOption( "-ndeviation" ,"-ndeviation: compute normal deviation." );
  args.addOption( "-mcurv" ,"-mcurv <sigma>: compute curvature averaged with s<sigma>.", "1.0" );

  if ( ( argc <= 1 ) 
       || ! args.readArguments( argc, argv ) ) 
    {
      cerr << args.usage( "test_nD", 
			  "Tests some nD features.",
			  "" )
	   << endl;
      return 1;
    }

  KnTypes::displayBaseTypesSize();
  uint d = StandardArguments::dim( args );
  Kn_size sizes[ d ];
  StandardArguments::fillSizes( args, sizes );
  KnSpace ks( d, sizes );

  float r = args.getOption( "-r" )->getFloatValue( 0 );
  float r1 = args.getOption( "-ri" )->getFloatValue( 0 );
  float r2 = args.getOption( "-ri" )->getFloatValue( 1 );
  float r_incr = args.getOption( "-ri" )->getFloatValue( 2 );
  if ( args.check( "-r" ) )
    {
      r1 = r2 = r;
      r_incr = 1.0;
    }
  bool dnormals = args.check( "-dnormals" );
  bool planes = args.check( "-planes" );
  bool area = args.check( "-area" );
  bool area_euler = args.check( "-area_euler" );
  bool ndeviation = args.check( "-ndeviation" );
  bool mcurv = args.check( "-mcurv" );
  float sigma = mcurv 
    ? args.getOption( "-mcurv" )->getFloatValue( 0 ) 
    : 1.0;

  for ( ; r1 <= r2; r1 += r_incr )
    {
      r = r1;
      // center
      Kn_uid clow = ks.uspel( ks.ukcode( sizes ) );
      cout << "--- Space: " << ks << endl;
      cout << "--- creating sphere (r=" << r << ") -----" << endl;
      Clock::startClock();
      
      KnCharSet sph1 = KnShapes::umakeVolumicSphere( ks, clow, r );
      
      //    cout << "--- creating parallepiped (r=" << r << ") -----" << endl;
      Kn_uid low_uid = KnShapes::ugetCubeLowerBound( ks, clow, r + 1 );
      Kn_uid up_uid = KnShapes::ugetCubeUpperBound( ks, clow, r + 1 );
      //    Kn_uid low_uid1 = KnShapes::ugetCubeLowerBound( ks, clow, r - 3 );
      //    Kn_uid up_uid1 = KnShapes::ugetCubeUpperBound( ks, clow, r - 3 );
      
      //    KnCharSet par1 = KnShapes::umakeVolumicParallelepiped( ks,
      //    							 low_uid1,
      //    							 up_uid1 );
      //    cout << "--- parallepiped - sphere -----" << endl;
      //    KnCharSet voxset = par1 - sph1;
      
      cout << "--- sphere -----" << endl;
      KnCharSet voxset = sph1;
      long ti2 = Clock::stopClock();
      cout << "in " << ti2 << " ms." << endl;
      cout << "voxset   = " << voxset.nbElements() << " spels." << endl;
      // for ( KnCharSet::cell_iterator p = voxset.begin(); 
      // 	    p != voxset.end(); ++p )
      // 	{
      // 	  ks.displayKn_uid( *p, cout );
      // 	  cout << endl;
      // 	}
      // cout << endl;
      
      //   uint belx[ 100 ];
      //   belx[ 0 ] = ( sizes[ 0 ] & ~1 ) + 2 * r + 2;
      //   for ( i = 1; i < ks.dim(); i++ )
      //     belx[ i ] = ( sizes[ i ] & ~1 ) + 1;
      //   Kn_sid bel = ks.negkcode( belx ); // starting bel for tracking
      cout << "--- boundary of sphere -----" << endl;
      Clock::startClock();
      KnRCellSet digsurf = KnShapes::smakeBoundary( ks, voxset, 
						    low_uid, up_uid );
      long ti3 = Clock::stopClock();
      cout << "in " << ti3 << " ms." << endl;
      cout << "digsurf   = " << digsurf.nbElements() << " surfels." << endl;
      
      if ( dnormals ) 
	computeDiscreteNormals( ks, voxset, digsurf );
      if ( area )
	{
	  cout << "--- area by enumeration ---------------------" << endl;
	  float carea = computeAreaMeasure( ks, voxset, digsurf, false );
	  displayArea( d, r, carea);
	  outputArea( cerr, d, r, carea );
	}
      if ( area_euler )
	{
	  cout << "--- area by euler integration ------------------" << endl;
	  float carea = computeAreaMeasure( ks, voxset, digsurf, true );
	  displayArea( d, r, carea);
	  outputArea( cerr, d, r, carea );
	}
      if ( mcurv )
	computeMCurvatureMeasure( ks, voxset, digsurf, sigma );
      if ( ndeviation )
	computeDiscreteNormalsDeviation( ks, voxset, digsurf );
      if ( planes )
	computeTangentPlanes( ks, voxset, digsurf );
    }
}

int main( int argc, char** argv )
{
  return testKnMC( argc, argv );
}

