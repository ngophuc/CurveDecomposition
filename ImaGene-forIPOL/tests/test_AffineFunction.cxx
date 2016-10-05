///////////////////////////////////////////////////////////////////////////////
// Test module for affine functions
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cmath>
#include <cstdlib>

#include "ImaGene/base/Arguments.h"
#include "ImaGene/mathutils/PWFAffineFunction.h"

using namespace std;
using namespace ImaGene;
//using namespace imagelib;

static Arguments args;

int
main( int argc, char** argv ) 
{
  float m,x;
  int limit;
  
  cout << "--------------- Uniform laws -------------------------" << endl;
  
  PWFAffineFunction f1 = PWFAffineFunction::uniformDP( 0.0, 1.0 );
  PWFAffineFunction f2 = PWFAffineFunction::uniformDP( 0.5, 1.5 );
  PWFAffineFunction f1p2 = PWFAffineFunction::add( f1, f2 );
  f1p2 = PWFAffineFunction::multiply( f1p2, 0.5 );
  
  cout << f1 << endl;
  cout << "f1: ";
  for ( float x = -1.0; x < 2.0; x += 0.1 )
    cout << "(" << x << " " << f1.data().value( x ) << ") ";
  cout << endl;
  cout << "f1(1)=" << f1.data().value(1.0) << endl;
  cout << "E(f1)=" << f1.data().centroid() << endl;
  cout << "D(f1)=" << sqrt( f1.data().variance() ) << endl;
  m = f1.data().firstMax( x, limit );
  cout << "Max(f1)=" << m << " at " << x << " with limit=" << limit << endl;
  
  cout << f2 << endl;
  cout << "f2: ";
  for ( float x = -1.0; x < 2.0; x += 0.1 )
    cout << "(" << x << " " << f2.data().value( x ) << ") ";
  cout << endl;
  cout << "f2(1.5)=" << f2.data().value(1.5) << endl;
  cout << "E(f2)=" << f2.data().centroid() << endl;
  cout << "D(f2)=" << sqrt( f2.data().variance() ) << endl;
  m = f2.data().firstMax( x, limit );
  cout << "Max(f2)=" << m << " at " << x << " with limit=" << limit << endl;

  cout << f1p2 << endl;
  cout << "f1 + f2: ";
  for ( float x = -1.0; x < 2.0; x += 0.1 )
    cout << "(" << x << " " << f1p2.data().value( x ) << ") ";
  cout << endl;
  cout << "E(f1 + f2)=" << f1p2.data().centroid() << endl;
  cout << "D(f1 + f2)=" << sqrt( f1p2.data().variance() ) << endl;
  m = f1p2.data().firstMax( x, limit );
  cout << "Max(f1 + f2)=" << m << " at " << x << " with limit=" << limit << endl;


  cout << "--------------- Triangular laws -------------------------" << endl;
  
  PWFAffineFunction tf1 = PWFAffineFunction::triangularDP( 0.0, 1.0 );
  PWFAffineFunction tf2 = PWFAffineFunction::triangularDP( 0.6, 1.4 );
  PWFAffineFunction tf1p2 = PWFAffineFunction::add( tf1, tf2 );
  tf1p2 = PWFAffineFunction::multiply( tf1p2, 0.5 );
  
  cout << tf1 << endl;
  cout << "tf1: ";
  for ( float x = -1.0; x < 2.0; x += 0.1 )
    cout << "(" << x << " " << tf1.data().value( x ) << ") ";
  cout << endl;
  cout << "tf1(1)=" << tf1.data().value(1.0) << endl;
  cout << "E(tf1)=" << tf1.data().centroid() << endl;
  cout << "D(tf1)=" << sqrt( tf1.data().variance() ) << endl;
  
  cout << tf2 << endl;
  cout << "tf2: ";
  for ( float x = -1.0; x < 2.0; x += 0.1 )
    cout << "(" << x << " " << tf2.data().value( x ) << ") ";
  cout << endl;
  cout << "tf2(1.5)=" << tf2.data().value(1.5) << endl;
  cout << "E(tf2)=" << tf2.data().centroid() << endl;
  cout << "D(tf2)=" << sqrt( tf2.data().variance() ) << endl;

  cout << tf1p2 << endl;
  cout << "tf1 + tf2: ";
  for ( float x = -1.0; x < 2.0; x += 0.1 )
    cout << "(" << x << " " << tf1p2.data().value( x ) << ") ";
  cout << endl;
  cout << "E(tf1 + tf2)=" << tf1p2.data().centroid() << endl;
  cout << "D(tf1 + tf2)=" << sqrt( tf1p2.data().variance() ) << endl;


  cout << "--------------- Trapezoidal laws -------------------------" << endl;
  
  PWFAffineFunction pf1 = PWFAffineFunction::trapezoidalDP( 0.0, 1.0, 0.2 );
  PWFAffineFunction pf2 = PWFAffineFunction::trapezoidalDP( 0.6, 1.4, 0.4 );
  PWFAffineFunction pf1p2 = 0.5 * ( pf1 + pf2 );
//   PWFAffineFunction::add( pf1, pf2 );
//   pf1p2 = PWFAffineFunction::multiply( pf1p2, 0.5 );
  
  cout << pf1 << endl;
  cout << "pf1: ";
  for ( float x = -1.0; x < 2.0; x += 0.1 )
    cout << "(" << x << " " << pf1.data().value( x ) << ") ";
  cout << endl;
  cout << "pf1(1)=" << pf1.data().value(1.0) << endl;
  cout << "E(pf1)=" << pf1.data().centroid() << endl;
  cout << "D(pf1)=" << sqrt( pf1.data().variance() ) << endl;
  
  cout << pf2 << endl;
  cout << "pf2: ";
  for ( float x = -1.0; x < 2.0; x += 0.1 )
    cout << "(" << x << " " << pf2.data().value( x ) << ") ";
  cout << endl;
  cout << "pf2(1.5)=" << pf2.data().value(1.5) << endl;
  cout << "E(pf2)=" << pf2.data().centroid() << endl;
  cout << "D(pf2)=" << sqrt( pf2.data().variance() ) << endl;

  cout << pf1p2 << endl;
  cout << "pf1 + pf2: ";
  for ( float x = -1.0; x < 2.0; x += 0.1 )
    cout << "(" << x << " " << pf1p2.data().value( x ) << ") ";
  cout << endl;
  cout << "E(pf1 + pf2)=" << pf1p2.data().centroid() << endl;
  cout << "D(pf1 + pf2)=" << sqrt( pf1p2.data().variance() ) << endl;
}

