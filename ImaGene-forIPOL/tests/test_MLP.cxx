///////////////////////////////////////////////////////////////////////////////
// Test MLP
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <iomanip>
#include "ImaGene/dgeometry2d/FreemanChain.h"

using namespace std;
using namespace ImaGene;

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
  FreemanChain c;
  c.chain = "01010303032323212121";
  double len = FreemanChain::lengthMLP( c );
  cout << setprecision(15) << len << endl;

  c.chain = "10103030323232121210";
  len = FreemanChain::lengthMLP( c );
  cout << setprecision(15) << len << endl;

  c.chain = "01030303232321212101";
  len = FreemanChain::lengthMLP( c );
  cout << setprecision(15) << len << endl;

  c.chain = "10303032323212121010";
  len = FreemanChain::lengthMLP( c );
  cout << setprecision(15) << len << endl;


  c.chain = "00333332222211100110";
  len = FreemanChain::lengthMLP( c );
  cout << setprecision(15) << len << endl;

  c.chain = "00033333222221110011";
  len = FreemanChain::lengthMLP( c );
  cout << setprecision(15) << len << endl;
}
