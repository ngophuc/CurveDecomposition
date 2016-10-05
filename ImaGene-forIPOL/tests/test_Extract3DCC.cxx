#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>

#include "ImaGene/helper/ShapeHelper.h"
#include "ImaGene/digitalnD/KnShapes.h"


using namespace ImaGene;
using namespace std;

static Arguments args;


int
main( int argc, char** argv ) 
{ 
  args.addOption("-threshold", "-threshold <min> <max> default: 0-128","0", "128" );
  args.addOption("-min_size", "-min_size <min> default 100 ", "100");
  
  if ( ( argc <= 1 ) || ! args.readArguments( argc, argv )){
    cerr << args.usage( "test_Extract3DCC",
			"Extract 3D Connected Compoment from a 3D pgm file and threshold values ",
			"-threshold -min_size " ) << endl;
    return 1;
  }
  
  KnSpace *ks;
  KnCharSet *voxset;
  
  int thresholdMin = args.getOption("-threshold")->getIntValue(0);
  int thresholdMax = args.getOption("-threshold")->getIntValue(1);  
  uint min_size = args.getOption("-min_size")->getIntValue(0);
  
  ShapeHelper::importFromPGM3d(cin, ks, voxset, thresholdMin, thresholdMax);    
  
  KnRCellSet bdry = KnShapes::smakeBoundary( *ks, *voxset );  
  KnRCellSet markedSurfelsSet = KnRCellSet::create(*ks, 2, true,0);
  KnRUCellVector<int> image3DLabel(*ks, 3);
  KnCharSet bounds = ~(*voxset);
  
  int nbDetected =0;
  int nbExtracted=0;
  KnRCellSet::cell_iterator boundaryIter = bdry.begin();
  while( boundaryIter!=bdry.end()){
    Kn_sid s = *boundaryIter;    
    
    uint orth_dir = ks->sorthDir( s );
    uint orth_direct = ks->sdirect( s, orth_dir );
    Kn_uid intVox = ks->unsigns( ks->sincident( s, orth_dir, orth_direct ) );
    
    if(!markedSurfelsSet[s]){
      KnRCellSet setConnectedBoundary = KnShapes::strackBoundary (*ks, *voxset, s);          
      cerr << "[Found] tracking connected surface  " << nbDetected;
      for ( KnRCellSet::cell_iterator cellBd_it = setConnectedBoundary.begin();
	    cellBd_it != setConnectedBoundary.end();
	    ++cellBd_it ){
	Kn_sid bel = *cellBd_it;
	markedSurfelsSet+=bel;	
      }
      
      KnCharSet interior = KnShapes::uexpandSeedWithinBounds(*ks, intVox, bounds);
      uint nb=interior.nbElements();
      cerr << " [nb el = " << nb << "] "; 
      if(nb>min_size && image3DLabel[intVox]==0){
	for ( KnCharSet::cell_iterator cellInt_it = interior.begin();
	      cellInt_it != interior.end();
	      ++cellInt_it ){
	  Kn_uid bel = *cellInt_it;
	  image3DLabel[bel]=nbExtracted+1;	 
	}
	cerr << " [Saved] "<< endl;	  
	nbExtracted++;
      }else
	  cerr << endl;
      ++nbDetected;
    }
    ++boundaryIter;
  }
  cerr << nbExtracted << " Connected Components extracted from a total of " << 
    nbDetected << " Connected Components" << endl;
  
  ShapeHelper::exportToPGM3d(cout, ks, image3DLabel);
  
}

