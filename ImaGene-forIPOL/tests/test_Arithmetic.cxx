///////////////////////////////////////////////////////////////////////////////
// Test module for some arithmetic properties
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cmath>
#include <cstdlib>

#include "ImaGene/base/Arguments.h"
#include "ImaGene/timetools/Clock.h"
#include "ImaGene/mathutils/Mathutils.h"
#include "ImaGene/mathutils/Statistics.h"

using namespace std;
using namespace ImaGene;

static Arguments args;

/**
   The greatest common denominator (GCD) is the largest positive integer
   that divides into both numbers without a remainder.
   Examples: GCD(256,64)=64, GCD(12,8)=4, GCD(5,3)=1 
*/
int GCD(int x, int y)
{
 int g;

 /* Work with absolute values (positive integers) */
 if (x < 0) x = -x;
 if (y < 0) y = -y;

 if ((x+y) == 0) return 0; /* Error, both parameters zero */

 g = y;
 /* Iterate until x = 0 */
 while(x > 0)
 {
  g = x;
  x = y % x;
  y = g;
 }

 return g;
}

void statsPartialQuotients( int n )
{
  // 0: average delta
  // 1: average p
  // 2: average q
  // 3: average u
  // 4: average sum_u
  Clock::startClock();
  cerr << "--- statsPartialQuotients( " << n << " ) in ";
  cerr.flush();

  Statistics stats( 6 );
  vector<uint> u;
  u.reserve( n );
  for ( int q = 1; q <= n; ++q ) 
    {
      if ( ( q % 1000 ) == 0 ) cerr << "." << flush;
      for ( int p = 1; p <= q; ++p )
	{
	  int g = GCD( p, q );
	  int a = p / g;
	  int b = q / g;
	  Mathutils::cfrac( u, a, b );
	  uint k = u.size();
	  uint sum = 0;
	  for ( uint i = 0; i < k; ++i )
	    {
	      uint v = u[ i ];
	      sum += v;
	    }
	  u.clear();
	  stats.addValue( 0, (double) g );
	  stats.addValue( 1, (double) a );
	  stats.addValue( 2, (double) b );
	  stats.addValue( 3, (double) sum / (double) k );
	  stats.addValue( 4, (double) sum );
	  stats.addValue( 5, (double) k );
	}
    }
  stats.terminate();

  long t = Clock::stopClock();
  cerr << t << " ms." << endl;

  cout << "# Arithmetic stats for 1 <= p <= q <= n" << endl
       << "# g = gcd(p,q), [u_0, ..., u_k ] = cfrac( p/q ), s = sum(u)" << endl
       << "# n E(g) E(p/g) E(q/g) E(u) E(s) E(k) Dev(g) Dev(p) Dev(q) Dev(u) Dev(s) Dev(k)" << endl
       << n 
       << " " << stats.mean( 0 )
       << " " << stats.mean( 1 )
       << " " << stats.mean( 2 )
       << " " << stats.mean( 3 )
       << " " << stats.mean( 4 )
       << " " << stats.mean( 5 )
       << " " << sqrt( stats.variance( 0 ) )
       << " " << sqrt( stats.variance( 1 ) )
       << " " << sqrt( stats.variance( 2 ) )
       << " " << sqrt( stats.variance( 3 ) )
       << " " << sqrt( stats.variance( 4 ) )
       << " " << sqrt( stats.variance( 5 ) )
       << endl;
}

void statsPartialQuotientsExact( int n )
{
  // 0: average delta
  // 1: average p
  // 2: average q
  // 3: average u
  // 4: average sum_u
  Clock::startClock();
  cerr << "--- statsPartialQuotientsExact( " << n << " ) in ";
  cerr.flush();

  Statistics stats( 6 );
  vector<uint> u;
  u.reserve( n );
  for ( int p = 1; p <= n / 2; ++p ) 
    {
      if ( ( p % 1000 ) == 0 ) cerr << "." << flush;
      int q = n - p;
      int g = GCD( p, q );
      int a = p / g;
      int b = q / g;
      Mathutils::cfrac( u, a, b );
      uint k = u.size();
      uint sum = 0;
      for ( uint i = 0; i < k; ++i )
	{
	  uint v = u[ i ];
	  sum += v;
	}
      u.clear();
      stats.addValue( 0, (double) g );
      stats.addValue( 1, (double) a );
      stats.addValue( 2, (double) b );
      stats.addValue( 3, (double) sum / (double) k );
      stats.addValue( 4, (double) sum );
      stats.addValue( 5, (double) k );
    }
  stats.terminate();

  long t = Clock::stopClock();
  cerr << t << " ms." << endl;

  cout << "# Arithmetic stats for 1 <= p <= q <= n" << endl
       << "# g = gcd(p,q), [u_0, ..., u_k ] = cfrac( p/q ), s = sum(u)" << endl
       << "# n E(g) E(p/g) E(q/g) E(u) E(s) E(k) Dev(g) Dev(p) Dev(q) Dev(u) Dev(s) Dev(k)" << endl
       << n 
       << " " << stats.mean( 0 )
       << " " << stats.mean( 1 )
       << " " << stats.mean( 2 )
       << " " << stats.mean( 3 )
       << " " << stats.mean( 4 )
       << " " << stats.mean( 5 )
       << " " << sqrt( stats.variance( 0 ) )
       << " " << sqrt( stats.variance( 1 ) )
       << " " << sqrt( stats.variance( 2 ) )
       << " " << sqrt( stats.variance( 3 ) )
       << " " << sqrt( stats.variance( 4 ) )
       << " " << sqrt( stats.variance( 5 ) )
       << endl;
}

int 
main( int argc, char** argv ) 
{
  args.addOption( "-stats_partial_quotients", "-stats_partial_quotients <n>: computes some statistics on partial quotients of fractions with integer in 1..n.", "100" );
  args.addOption( "-stats_partial_quotients_exact", "-stats_partial_quotients_exact <n>: computes some statistics on partial quotients of fractions a/b a+b=n.", "100" );
  args.addOption( "-exp_stats_partial_quotients", "-exp_stats_partial_quotients <n>: computes some statistics on partial quotients of fractions with integer in 1..k, with k growing exponentially from 1 to n.", "100" );
  args.addOption( "-lin_stats_partial_quotients", "-lin_stats_partial_quotients <n>: computes some statistics on partial quotients of fractions with integer in 1..k, with k growing linearly from 1 to n.", "100" );
  
  if ( ( argc <= 1 ) 
       || ! args.readArguments( argc, argv ) ) 
    {
      cerr << args.usage( "test_TangentialCover", 
			  "Computes some statistics related to arithmetic properties.",
			  "" )
	   << endl;
      return 1;
    }
  
  if ( args.check( "-stats_partial_quotients" ) ) 
    {
      int n = args.getOption( "-stats_partial_quotients" )->getIntValue( 0 );
      statsPartialQuotients( n );
    }

  if ( args.check( "-stats_partial_quotients_exact" ) ) 
    {
      int n = args.getOption( "-stats_partial_quotients_exact" )->getIntValue( 0 );
      statsPartialQuotientsExact( n );
    }

  if ( args.check( "-exp_stats_partial_quotients" ) ) 
    {
      int n = args.getOption( "-exp_stats_partial_quotients" )->getIntValue( 0 );
      int i = 1;
      while ( i <= n ) 
	{
	  statsPartialQuotients( i );
	  i = i * 2;
	}
    }

  if ( args.check( "-lin_stats_partial_quotients" ) ) 
    {
      int n = args.getOption( "-lin_stats_partial_quotients" )->getIntValue( 0 );
      for ( int i = 1; i <= n; i++ )
	{
	  statsPartialQuotients( i );
	}
    }
    

  return 0;
}
