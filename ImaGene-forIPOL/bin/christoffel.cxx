///////////////////////////////////////////////////////////////////////////////
// Creates a Freeman chaincode corresponding to a Christoffel word.
///////////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <cstdlib>
#include <iostream>

#include "ImaGene/base/Arguments.h"
#include "ImaGene/mathutils/Mathutils.h"
#include "ImaGene/dgeometry2d/FreemanChain.h"


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

string christoffel( vector<uint> & z, uint n, char zero, char one )
{
  string s = "";
  if ( n == 0 ) s += zero;
  else if ( n == 1 ) 
    {
      for ( uint i = 0; i < z[ n ]; ++i )
	s += zero;
      s += one;
    }
  else if ( n % 2 == 0 ) 
    {
      s += christoffel( z, n - 2, zero, one );
      for ( uint i = 0; i < z[ n ]; ++i )
	s += christoffel( z, n - 1, zero, one );
    }
  else
    {
      for ( uint i = 0; i < z[ n ]; ++i )
	s += christoffel( z, n - 1, zero, one );
      s += christoffel( z, n - 2, zero, one );
    }
  return s;
}

int
main( int argc, char** argv ) 
{
  // -------------------------------------------------------------------------
  // Prepare arguments.
  args.addOption( "-slope", " -slope <a> <b>: specifies the slope a/b of the Christoffel word.", "0", "1" );  
  args.addBooleanOption( "-close", " -close: close the Freeman chain to get a closed 4-connected contour." );

  if ( ( argc <= 0 ) 
       || ! args.readArguments( argc, argv ) ) 
    {
      cerr << args.usage( "christoffel", 
			  "Creates the Freeman chain code corresponding to the Christoffel word of the given slope a/b.",
			  "" )
	   << endl;
      return 1;
    }

  int a = args.getOption( "-slope" )->getIntValue( 0 );
  int b = args.getOption( "-slope" )->getIntValue( 1 );
  bool close = args.check( "-close" );
  char zero = '0';
  char one = close ? '3' : '1';
  vector<uint> z;
  Mathutils::cfrac( z, a, b );
  string s = ( a <= b ) 
    ? christoffel( z, z.size() - 1, zero, one )
    : christoffel( z, z.size() - 1, one, zero );
  cout << "# christoffel -slope " << a << " " << b << endl;
  if ( close )
    {
      cout << "1 " << a+1 << " " << s << zero;
      for ( uint i = 0; i <= a; ++i ) cout << "1";
      for ( uint i = 0; i <= b; ++i ) cout << "2";
      cout << "3" << endl;
    }
  else
    cout << s << endl;
}
