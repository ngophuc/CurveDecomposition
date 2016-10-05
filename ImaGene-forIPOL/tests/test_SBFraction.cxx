///////////////////////////////////////////////////////////////////////////////
// Test fractions implemented via Stern-Brocot tree.
///////////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <cstdlib>
#include <iostream>
#include "ImaGene/base/Vector2i.h"
#include "ImaGene/base/Arguments.h"
#include "ImaGene/base/StandardArguments.h"
#include "ImaGene/mathutils/SternBrocot.h"
#include "ImaGene/timetools/Clock.h"


using namespace std;
using namespace ImaGene;


static Arguments args;

int
main( int argc, char** argv ) 
{
  // -------------------------------------------------------------------------
  // Prepare arguments.
  args.addOption( "-fraction", "-fraction <p> <q>: display various development of this fraction.", "3", "8" );
  args.addOption( "-benchmark", "-benchmark <nb> <p> <q>: computes <nb> fractions a/b randomly chosen such that 0 < a <= p and 0 < b <= q.", "100000", "1000", "1000" );
  if ( ( argc <= 0 ) 
       || ! args.readArguments( argc, argv ) ) 
    {
      cerr << args.usage( "test_Math", 
			  "Test module SternBrocot fractions.",
			  "" ) << endl;
      return 1;
    }
  
  typedef SternBrocot::SBFraction Fraction;

  if ( args.check( "-fraction" ) )
    {
      unsigned int p = args.getOption( "-fraction" )->getIntValue( 0 );
      unsigned int q = args.getOption( "-fraction" )->getIntValue( 1 );
      Fraction f = SternBrocot::fraction( p, q );
      cout << "---- Reduced fractions of f ---" << endl;
      cout << "f=" << f << endl;
      int i = 0;
      Fraction g = f;
      while ( g.k() > 0 )
	{
	  g = f.reduced( i );
	  cout << "f_" << f.k()-i << "=" << g << endl;
	  i++;
	}
      cout << "---- Splits de f ---" << endl;
      Fraction f1, f2;
      unsigned int n1, n2;
      f.getSplit( f1, f2 );
      cout << f1 << " . " << f2 << endl;
      f.getSplitBerstel( f1, n1, f2, n2 );
      cout << "(" << f1 << ")^" << n1 << " . " 
	   << "(" << f2 << ")^" << n2 << endl;
    }
  if ( args.check( "-benchmark" ) )
    {
      cout << "---- benchmarking fractions ---" << endl;
      Clock::startClock();
      unsigned int nb =  args.getOption( "-benchmark" )->getIntValue( 0 );
      unsigned int p = args.getOption( "-benchmark" )->getIntValue( 1 );
      unsigned int q = args.getOption( "-benchmark" )->getIntValue( 2 );

      unsigned long long profs = 0;
      for ( unsigned int i = 0; i < nb; ++i )
	{
	  int a = ( random() % p ) + 1;
	  int b = ( random() % q ) + 1;
	  int g = Mathutils::gcd( a, b );
	  Fraction f = SternBrocot::fraction( a/g, b/g );
	  profs += (unsigned long long) f.k();
	} 
      long t = Clock::stopClock();
      double moy = (double) ( ( (double) profs ) / nb );
      cout << "     p_moy=" << moy << " in " << t << " ms." << endl;
    }
  return 0;
}
