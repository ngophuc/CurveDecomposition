///////////////////////////////////////////////////////////////////////////////
// Test class OrderedAlphabet.
///////////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>

#include "ImaGene/base/Arguments.h"
#include "ImaGene/base/StandardArguments.h"
#include "ImaGene/base/OrderedAlphabet.h"

using namespace std;
using namespace ImaGene;


static Arguments args;

int main( int argc, char** argv )
{
  // -------------------------------------------------------------------------
  // Prepare arguments.
  StandardArguments::addDebugArgs( args, true, true, true );
  StandardArguments::addIOArgs( args, true, false );
  args.addBooleanOption( "-test_FLF", "-test_FLF: test first Lyndon Factor." );
  
  if ( ( argc <= 0 ) 
       || ! args.readArguments( argc, argv ) ) 
    {
      cerr << args.usage( "test_OrderedAlphabet", 
			  "Test class OrderedAlphabet.",
			  "" )
	   << endl;
      return 1;
    }
  int err = 0;
  err += OrderedAlphabet::testClass( args );
  // if ( args.check( "-test_slinel" ) )
  //   {
  //     err += test_K2Space_slinel( K2, contour, c );
  //   }
  return err;
}
