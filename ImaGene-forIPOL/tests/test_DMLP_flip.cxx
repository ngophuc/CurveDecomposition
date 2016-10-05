#include <cassert>
#include <iostream>
#include "ImaGene/dgeometry2d/FreemanChain.h"
#include "ImaGene/dgeometry2d/GridCurve.h"
#include "ImaGene/digitalnD/KnSpace.h"
#include "ImaGene/digitalnD/K2Space.h"
#include "ImaGene/digitalnD/KnShapes.h"
#include "ImaGene/helper/ShapeHelper.h"


using namespace ImaGene;





void displayDMLP(GridCurve& dmlp)
{
  std::cout<<dmlp.getFreemanChain();
  /*  GridCurve::Visitor itm = dmlp.begin();
  std::cout<<"Base: "<<(*itm).first<<": ";
  do
    {
       std::cout<<(int)(*itm).second<<" ";
       ++itm;
    }
    while (itm!=dmlp.begin());*/
}

bool check_visitor_coordinates( GridCurve & dmlp )
{
  GridCurve::Visitor visit_begin = dmlp.begin(); 
  GridCurve::Visitor visit = visit_begin;
  GridCurve::Visitor visit_end = dmlp.end();
  FreemanChain fc = dmlp.getFreemanChain();
  FreemanChain::const_iterator it = fc.begin();
  FreemanChain::const_iterator it_end = fc.end();
  bool ok = true;
  while ( it != it_end )
    {
      GridCurve::Visitor::Value visit_value = *visit;
      Vector2i it_v = *it;
      unsigned int it_step = it.getCode();
      std::cerr << "[check] visitor: "
                << visit_value.first << " " << (int) visit_value.second 
                << " fc: " << it_v << " " << (int) it_step;
      if ( ( visit_value.first != it_v ) || ( visit_value.second != it_step ) )
        {
          std::cerr << " ERR " << *(visit.edgeIterator());
          ok = false;
        }
      else
        std::cerr << " OK  " << *(visit.edgeIterator());;
      std::cerr << std::endl;
      ++it; ++visit;
    }
  return ok;
}


double diffLength( GridCurve & dmlp, double expectedLength) 
{
  double l = 0.0;
  for (GridCurve::EdgeListIterator it = dmlp.beginEdge(), itEnd= dmlp.endEdge(); it != itEnd; ++it) 
  {
    if ( ! it->isQuadrantChange() ) {
      l += it->euclideanLength();
    } 
  }
  return (l > expectedLength) ? l - expectedLength : expectedLength - l;
}

// Compute a digital space big enough for containing the two given
// freeman chaincodes.
K2Space* 
makeSpaceFromTwoFreemanChains( FreemanChain & c1, 
			       FreemanChain & c2, Vector2i & delta )
{
  Kn_ssize min_x,min_x2;
  Kn_ssize min_y,min_y2;
  Kn_ssize max_x,max_x2;
  Kn_ssize max_y,max_y2;
  c1.computeBoundingBox( min_x, min_y, max_x, max_y );
  c2.computeBoundingBox( min_x2, min_y2, max_x2, max_y2 );
  min_x = min_x2 < min_x ? min_x2 : min_x;
  min_y = min_y2 < min_y ? min_y2 : min_y;
  max_x = max_x2 > max_x ? max_x2 : max_x;
  max_y = max_y2 > max_y ? max_y2 : max_y;

  Kn_size width = (Kn_size) ( 3 + max_x - min_x );
  Kn_size height = (Kn_size) ( 3 + max_y - min_y );
  Kn_size i = 1;
  while ( i < ( width + 2 ) ) i <<=1;
  Kn_size j = 1;
  while ( j < ( height + 2 ) ) j <<=1;
  
  delta = Vector2i( 1+( (Kn_ssize) i - max_x - min_x ) / 2,
		    1+( (Kn_ssize) j - max_y - min_y ) / 2 );
  c1.x0 += 1+( (Kn_ssize) i - max_x - min_x ) / 2;
  c1.y0 += 1+( (Kn_ssize) j - max_y - min_y ) / 2;
  c2.x0 += 1+( (Kn_ssize) i - max_x - min_x ) / 2;
  c2.y0 += 1+( (Kn_ssize) j - max_y - min_y ) / 2;
  Kn_size sizes[ 2 ];
  sizes[ 0 ] = i;
  sizes[ 1 ] = j;
  return new K2Space( i, j ); 
  // return new KnSpace( 2, sizes );
}

void displayPixels( ostream & out, 
		    const K2Space & ks, 
		    KnCharSet s, Vector2i delta )
{
  for ( KnCharSet::cell_iterator it = s.begin(), it_end = s.end();
	it != it_end; ++it )
    out << " (" << ( ks.ux( *it ) - delta.x() ) 
	<< "," << ( ks.uy( *it ) -delta.y() ) << ")";
  out << std::endl;
}


void cleanFreemanChain ( FreemanChain & clean, const FreemanChain & c ) 
{
  //cerr << "[DIRTY FREEMANCHAIN] " << c.chain << endl;
  std::vector<unsigned int> v1, v2;
  bool ccw = c.isClosed() > 0;
  FreemanChain tmp;
  FreemanChain::cleanOuterSpikes( tmp, v1, v2, c, ccw);
  FreemanChain::cleanOuterSpikes( clean, v1, v2, tmp, ! ccw);
  //cerr << "[CLEAN FREEMANCHAIN] " << clean.chain << endl;
}



KnRCellSet inverseCells( const KnSpace & ks, KnRCellSet cells )
{
    KnRCellSet inverse = KnRCellSet::create( ks, cells.dim(), true, 0 );
    for ( KnRCellSet::cell_iterator it = cells.begin(), it_end = cells.end(); 
          it != it_end; ++it )
    {
      inverse += ks.sopp( *it );
    }
    return inverse;
}

// Compute the interior pixels to a freeman chain code.
// @returns 0 if c1 == c2, 
//          1 if c1 has exactly one more element than c2
//         -1 if c2 has exactly one more element than c1
//          2 otherwise.
//
// if returns 1 or -1, diff is the corresponding pixel.
int checkDifference ( FreemanChain c1, FreemanChain c2, Vector2i & diff  )
{ 
  Vector2i delta;
  //K2Space* ks = makeSpaceFromTwoFreemanChains( c1, c2, delta );
  K2Space* ks = new K2Space( 256, 256 ) ;
  delta = Vector2i( 0, 0 );
  KnRCellSet contour1 = ShapeHelper::makeContourFromFreemanChain( ks, c1, true );
  if ( c1.isClosed() < 0 ) 
    contour1 = inverseCells( *ks, contour1 ); 
  KnCharSet image1 = KnShapes::ucomputeInterior( *ks, contour1, false );
  KnRCellSet contour2 = ShapeHelper::makeContourFromFreemanChain( ks, c2, true );
  if ( c2.isClosed() < 0 ) 
    contour2 = inverseCells( *ks, contour2 ); 
  KnCharSet image2 = KnShapes::ucomputeInterior( *ks, contour2, false );
  KnCharSet minus12 = image1 - image2;
  KnCharSet minus21 = image2 - image1;
  bool ok = true;

  //////////////////////
//  cerr << "[CHECKDIFFERENCE minus12 ]" ;
//  displayPixels( std::cerr, *ks, minus12, delta );
//  cerr << "[CHECKDIFFERENCE minus21 ]" ;
//  displayPixels( std::cerr, *ks, minus21, delta );
  //////////////////////

  if ( minus12.nbElements() == 1 )
    {
      diff = Vector2i ( ks->ux( *minus12.begin() ) - delta.x() , 
                        ks->uy( *minus12.begin() ) - delta.y() );
      return (minus21.nbElements() == 0) ? 1 : 2 ;
    }

  if ( minus21.nbElements() == 1 )
    {
      diff = Vector2i ( ks->ux( *minus21.begin() ) - delta.x() , 
                        ks->uy( *minus21.begin() ) - delta.y() );
      return (minus12.nbElements() == 0) ? -1 : 2 ;
    }
  
  return ( (minus12.nbElements() == 0) && ( minus21.nbElements() == 0) ) ? 0 : 2;
}



bool testFlipPosition( FreemanChain ch, Vector2i pos, unsigned char step, bool inside );

bool checkFlipPosition( FreemanChain ch )
{
  GridCurve dmlp(ch, "INTERPIXEL");
  #ifdef EXTRA_OUTPUT_DMLP
  std::cout << "original_shape_from_constructor = " << dmlp << endl;
  #endif
  //GridCurve dmlp;
  //dmlp.initFromInterpixelFreemanChain( ch );
  //#ifdef EXTRA_OUTPUT_DMLP
  //std::cout << "original_shape_from_cmlp = " << dmlp << endl;
  //#endif

  GridCurve::Visitor itmlp = dmlp.begin();
  GridCurve::Visitor itBegin = itmlp;
  #ifdef EXTRA_OUTPUT_DMLP
  std::cout << "l_in = []\n";
  std::cout << "l_out = []\n";
  #endif
  bool ok;
  do 
    {
      ok = false;
      Vector2i pos = (*itmlp).first;
      GridCurve::Visitor::Value v = *itmlp;
      Vector2i in_pixel = itmlp.insidePixel();
      Vector2i out_pixel = itmlp.outsidePixel();
      GridCurve::EdgeListIterator e = itmlp.edgeIterator();
      #ifdef EXTRA_OUTPUT_DMLP
      std::cerr << "+- Flip at pointel=" 
        << " step=" << (char) ('0' + v.second)
        << " ("  <<  pos.x() << "," << pos.y() << ") "
        << " inpixel=("  <<  in_pixel.x() << "," << in_pixel.y() << ") "
        << " outpixel=("  <<  out_pixel.x() << "," << out_pixel.y() << ") "
        << " edge=" << *e 
        << std::endl;
      #endif

      dmlp.flipAndSimplify(itmlp,true,true);
      #ifdef EXTRA_OUTPUT_DMLP
      std::cout << "l_in.append(" << dmlp << ")\n";
      #endif
      Vector2i diff;
      FreemanChain dmlp_ch_flip;
      cleanFreemanChain( dmlp_ch_flip, dmlp.getFreemanChain() );
      int after_inside = checkDifference( ch, dmlp_ch_flip, diff );
      #ifdef EXTRA_OUTPUT_DMLP
      std::cerr << "  +- inside flip: ";
      #endif
      if ( after_inside * after_inside == 1 ) {
          // after_inside is 1 or -1.
          if ( diff == in_pixel ) {
              #ifdef EXTRA_OUTPUT_DMLP
              std::cerr << "OK" << endl;
              #endif
              ok = true;
          } else {
              std::cerr << "error : one pixel flipped at : " << diff << std::endl;
          } 
      } else {
          std::cerr << "error : " << (after_inside == 0 ? "No pixel flipped" : 
                                      "More then one (?) pixels flipped") << std::endl; 
          testFlipPosition ( ch, pos, v.second,true );
      }
      if (!ok) break;
      ok = false;

      dmlp.undoFlip();
      dmlp.flipAndSimplify(itmlp,false,true);
      #ifdef EXTRA_OUTPUT_DMLP
      std::cout << "l_out.append(" << dmlp << ")\n";
      #endif

      cleanFreemanChain( dmlp_ch_flip, dmlp.getFreemanChain() );
      int after_outside = checkDifference( ch, dmlp_ch_flip, diff );
      #ifdef EXTRA_OUTPUT_DMLP
      std::cerr << "  +- outside flip: ";
       #endif
      if ( after_outside * after_outside == 1 ) {
          // after_outside is 1 or -1.
          if ( diff == out_pixel ) {
              #ifdef EXTRA_OUTPUT_DMLP
              std::cerr << "OK" << endl;
              #endif
              ok = true;
          } else {
              std::cerr << "error : one pixel flipped at : " << diff << std::endl;
          } 
      } else {
          std::cerr << "error : " << (after_outside == 0 ? "No pixel flipped" : 
                                      "More then one (?) pixels flipped") << std::endl; 
          testFlipPosition ( ch, pos, v.second, false );
      }
      dmlp.undoFlip();
      ++itmlp;
    } 
  while ( ok && ( itmlp != itBegin ) );
  return ok;
}

// Compute the interior pixels to a freeman chain code.
bool compareInteriors( FreemanChain c1, FreemanChain c2 )
{
  Vector2i delta;
  K2Space* ks = makeSpaceFromTwoFreemanChains( c1, c2, delta );
  KnRCellSet contour1 = ShapeHelper::makeContourFromFreemanChain( ks, c1, true );
  KnCharSet image1 = KnShapes::ucomputeInterior( *ks, contour1, false );
  KnRCellSet contour2 = ShapeHelper::makeContourFromFreemanChain( ks, c2, true );
  KnCharSet image2 = KnShapes::ucomputeInterior( *ks, contour2, false );
  KnCharSet minus12 = image1 - image2;
  KnCharSet minus21 = image2 - image1;
  bool ok = true;
  if ( minus12.nbElements() != 0 )
    {
      std::cerr << "  [C1-C2 != 0]";
      displayPixels( std::cerr, *ks, minus12, delta );
      ok = false;
    }
  if ( minus21.nbElements() != 0 )
    {
      std::cerr << "  [C2-C1 != 0]";
      displayPixels( std::cerr, *ks, minus21, delta );
      ok = false;
    }
  return ok;
}

bool testFlipPosition( FreemanChain ch,
               Vector2i pos, unsigned char step, bool inside )
{
  GridCurve dmlp(ch, "INTERPIXEL");
  GridCurve::Visitor itmlp=dmlp.findPointel( pos.x(), pos.y(), step, dmlp.begin() );
  GridCurve::Visitor::Value v = *itmlp;
  Vector2i in_pixel = itmlp.insidePixel();
  Vector2i out_pixel = itmlp.outsidePixel();
  #ifdef EXTRA_OUTPUT_DMLP
  std::cerr << "+- Flip at pointel="
        << " step=" << (char) ('0' + v.second)
        << " ("  <<  pos.x() << "," << pos.y() << ") "
        << " inpixel=("  <<  in_pixel.x() << "," << in_pixel.y() << ") "
        << " outpixel=("  <<  out_pixel.x() << "," << out_pixel.y() << ") "

        << std::endl;
  std::cerr << "+- before flip: " << std::endl;
  #endif
  FreemanChain dmlp_ch ;
  cleanFreemanChain( dmlp_ch, dmlp.getFreemanChain() );

  #ifdef EXTRA_OUTPUT_DMLP
  std::cerr << "   +- original : " << "(" << ch.x0 << "," << ch.y0 << ") " << ch.chain << std::endl;
  std::cerr << "   +- dmlped   : " << "(" << dmlp_ch.x0 << "," << dmlp_ch.y0 << ") " << dmlp_ch.chain << std::endl;
  #endif
  bool ok_before = compareInteriors( ch, dmlp_ch );
  #ifdef EXTRA_OUTPUT_DMLP
  std::cerr << "   +- " << ( ok_before ? "OK" : "KO" ) << std::endl;
  #endif
  dmlp.flipAndSimplify(itmlp,inside,true);
  FreemanChain dmlp_ch_flip;
  cleanFreemanChain( dmlp_ch_flip, dmlp.getFreemanChain() ) ;
  bool ok_after = compareInteriors( ch, dmlp_ch_flip );

  #ifdef EXTRA_OUTPUT_DMLP
  std::cerr << "+- after flip: " << std::endl;
  std::cerr << "   +- " << ( ok_after ? "KO" : "OK ?" ) << std::endl;
  std::cerr << " DMLP = "  << dmlp << std::endl;
  #endif

  return ok_before;
} 

bool testDMLP(FreemanChain ch) {
  unsigned int nb = 0;
  unsigned int nbok = 0;


  GridCurve dmlp(ch, "INTERPIXEL");
  std::cout << "Input=" << ch << std::endl;
  std::cout << "DMLPc=" << dmlp.getFreemanChain() << std::endl;
  std::cout << "Input=" << ch.chain << std::endl;
  std::cout << "DMLPc=" << dmlp.getFreemanChain().chain << std::endl;
  std::cout << "DMLP =" << dmlp << std::endl;
  std::cout << "DMLPc=" << dmlp.getFreemanChain().chain << std::endl;
  for ( GridCurve::EdgeListIterator it = dmlp.beginEdge(),
	it_end = dmlp.endEdge(); it != it_end; ++it )
    {
      std::vector<unsigned char> steps;
      it->pushBackFreemanCode( steps );
      std::cout << *it << " ab=" << it->a << it->b << " => ";
      for ( unsigned int i = 0; i < steps.size(); ++i )
	std::cout << (char) ('0' + steps[ i ] );
      std::cout << std::endl;
    }
  check_visitor_coordinates( dmlp );

  //assert( dmlp.getFreemanChain() == ch );
  // std::cout << chaine << std::endl
  //           << dmlp_test.getFreemanChain() << std::endl;
  /*
  // Code to test if starting from ch and flip pointel 4,1 give ch2.
    displayDMLP(dmlp);  
    GridCurve::Visitor itmlp=dmlp.findPointel(4,1);
    dmlp.flipAndSimplify(itmlp,true,true);
    displayDMLP(dmlp);
    return EXIT_SUCCESS;
  */

  srand(time(NULL));



  std::cout << "//////////////////////////" << std::endl;
  
  //  GridCurve::Visitor itmlp=dmlp.findPointel(4,1);
  double diff = diffLength(dmlp, FreemanChain::lengthMLP( dmlp.getFreemanChain() ) );
  ++nb; nbok += ( diff < 0.0000001 ) ? 1 : 0;
  std::cout<<"\t: "
    <<" DMLP="<< dmlp.getLength() << ";  "
    <<" CMLP="<< FreemanChain::lengthMLP(ch);
  displayDMLP(dmlp);std::cout<< std::endl;
  std::cout << "DMLP =" << dmlp << std::endl;
  for (int i=0; i<ch.size()+1; ++i)
  {
    GridCurve::Visitor itmlp=dmlp.findPointel(0,0);
    //for (int j=0,k=rand()%100;j<k;++j)
    for (int j=0,k=i;j<k;++j)
      ++itmlp;
    for (int inside=0; inside < 2; inside++) {
    //bool inside=(rand()%2?true:false);

      std::cout << "\n\n" << std::endl;
      dmlp.flipAndSimplify(itmlp,inside,true);      
      std::cout<<"Run "<<i<<"\t: "
        <<" DMLP="<< dmlp.getLength() << ";  "
        <<" CMLP="<< FreemanChain::lengthMLP(dmlp.getFreemanChain())
        << std::endl;
      std::cout << "DMLP =" << dmlp << std::endl;
      double diff = diffLength(dmlp, FreemanChain::lengthMLP( dmlp.getFreemanChain() ));
      ++nb; nbok += ( diff  < 0.0000001 ) ? 1 : 0;
      if ( diff >= 0.0000001 ) {
        std::cout << "[BUG]\n buggedgc = " << dmlp << std::endl; 
        dmlp.undoFlip();
        std::cout << "[BUG]\n original = " << dmlp << std::endl; 
      } else {
        dmlp.undoFlip();
      }
      std::cout<<"\t: "
        <<" DMLP="<< dmlp.getLength() << ";  "
        <<" CPLM="<< FreemanChain::lengthMLP(ch);
      displayDMLP(dmlp);std::cout<< std::endl;
    }
  }
  std::cout << "(" << nbok << "/" << nb << ") test succesful" << std::endl;
  return nbok == nb;
}

int main()
{
  // CCW
  //  FreemanChain ch("2222222222222222222323232323232333333333333333303030303030303030000000000000001010101010101111111111111111111121212121", 27, 31);
  // CW
  // FreemanChain ch("0000000000000000000303030303030333333333333333323232323232323232222222222222221212121212121111111111111111111101010101", 27, 31 );
  // OK
  // CW and CCW seems to work fine.
  // FreemanChain ch("00001112222333",0,0);
  // FreemanChain ch("0000333322221111",0,0);
  // FreemanChain ch("00030333232221211101",0,0);
  // FreemanChain ch("00010111212223233303",0,0);
  // FreemanChain ch("2222222222222222222323232323232333333333333333303030303030303030000000000000001010101010101111111111111111111121212121", 27, 31);
  FreemanChain ch("2222222222222222222121212121212111111111111111101010101010101010000000000000003030303030303333333333333333333323232323", 27, 31);
  // FreemanChain ch("0001210100033332222221",0,0);
  // FreemanChain ch("00112233",0,0);
  // FreemanChain ch("03321210",0,0);
  // FreemanChain ch("01123230",0,0);
  // FreemanChain ch("011233",0,0);
  // FreemanChain ch("033211",0,0);

  // KO
  // Too small.
  // FreemanChain ch("0123",0,0);

  // ------------ check_visitor_coordinates.---------------------
  // OK 
  // FreemanChain ch("0001210100033332222221",0,0);
  //FreemanChain ch("00012221000100033332222221",0,0);

  // From Guillaume :
  // J'ai essayé de tracer un peu, apparemment dans
  // GridCurve::preSimplifyAt il y a un cas ou isMergeable retourne
  // vrai mais ou l'objet n'est pas modifié par le merge.
  // flipAndSimplify((13,10)) pourPolygon( [ Segment((1,2),tilde=false,nb=2),
  // Segment((2,1),tilde=false,nb=3), Segment((1,0),tilde=false,nb=4),
  // Operateur('plus'), Segment((1,3),tilde=false,nb=1),
  // Segment((1,1),tilde=false,nb=2), Segment((4,1),tilde=false,nb=1),
  // Segment((1,0),tilde=false,nb=2), Operateur('plus'),
  // Segment((1,4),tilde=false,nb=1), Segment((1,2),tilde=false,nb=1),
  // Segment((1,1),tilde=false,nb=1), Segment((3,1),tilde=false,nb=1),
  // Segment((1,0),tilde=false,nb=2), Operateur('plus'),
  // Segment((1,2),tilde=false,nb=1), Segment((3,2),tilde=true,nb=1),
  // Segment((2,3),tilde=true,nb=1), Segment((2,1),tilde=false,nb=2),
  // Segment((1,0),tilde=false,nb=1), Operateur('plus') ] ) 
  //FreemanChain ch("01101100100100100000300030303333033323333233232223222122121212212211211211",0,0);



  // ------
  //FreemanChain ch2("00001211222333",0,0);
  std::vector<FreemanChain> theChains;
  //theChains.push_back( FreemanChain("110301030103010303332222222211", 10, 12 ));
  //theChains.push_back( FreemanChain("01101100100100100000300030303333033323333233232223222122121212212211211211",20,20));
  //theChains.push_back( FreemanChain("1000010000100010000100010110111011011100033333333333333332222222222222222222222222",20,20));
  theChains.push_back( FreemanChain("00001112222333",20,20));
  theChains.push_back( FreemanChain("22221110000333",20,20));
  theChains.push_back( FreemanChain("00010300111111223223223333",20,20));
  theChains.push_back( FreemanChain("0000333322221111",20,20));
  theChains.push_back( FreemanChain("00030333232221211101",20,20));
  theChains.push_back( FreemanChain("00010111212223233303",20,20));
  theChains.push_back( FreemanChain("2222222222222222222323232323232333333333333333303030303030303030000000000000001010101010101111111111111111111121212121", 27, 31));
  theChains.push_back( FreemanChain("2222222222222222222121212121212111111111111111101010101010101010000000000000003030303030303333333333333333333323232323", 27, 31));
  theChains.push_back( FreemanChain("000121010003333322222211",20,20));
  theChains.push_back( FreemanChain("00112233",20,20));
  theChains.push_back( FreemanChain("03321210",20,20));
  theChains.push_back( FreemanChain("01123230",20,20));
  theChains.push_back( FreemanChain("011233",20,20));
  // this one requires 1024x1024 Knspace
  // theChains.push_back( FreemanChain("12121212212121221212122122122121221221222122122122212212221222122212222122212222122221222212222212222221222221222222212222222212222222221222222222222212222222222222222222222222222222222222222222222222222222222222222222222222222222221222121111010110101010101010101010101010101010101010101010100101010101010101010101010101010101010101010110101010101010101010101101010101010101101010101101010101101010101101010110101101010110101101011010110110101101011011011011010110110110110111011011011101110110111101110111110111111111101211111111211121121121212121212212212222122222222222232222232223222322232232223223223223223223232232232322323223232232322323232232323223232323223232323232322323232323232323232323232322323232332323232323232323232323232323233232323232323233232323232332323232332323232332323232332323233232323323232332323232332323233232323233232323232232222121112111112111111112111111111112111111111111211111111112111111111211111112111111121111111211111211111121111121111121111211112111121111211112111211112111211121112111211211121112112112111211211211211211211211211212112112121121212112121212121212121212121221212212212212222212223222223223232232323232323233232332332332332333233323332333233333233333233333332333333333323333333333333333333333333333333330333333333333033333333330333333303333333033333303333330333330333330333330333303333303333033330333303333033330333303333033330333303333303333333323222221221221221221212212122121221212212122121221212212122121221221212212122121221221212212212122122121221221221212212212212212212212212212212212212212212212212221221221222122122212212221222122122212221222212221222212222122222122222122222221222222222222222222223222223222232232232232323232323233233323333333333330333033303303303303303303033030330303033030303030303303030303030303030303003030303030303003030303003030303003030300303003030030303003030030030300303003003030030030300300303003003003003003030030030030030030030030030030030030003003003003003003000300300300300300030030030030030030030300303303233232232232232223222322232223222322232223222322232223222322232223222322232223222322322232223223222322322232232223223223223222322322322322322322322322322322322322322323223223223232232322322323223232232322323232232323223232323232232323232323232323232323233232323233232332332323332332333323333333333033303303030303030300303000300300030000300000030000000000000000000000000010000000010000001000001000001000001000010001000010001000010001000100010001000100010010001000100100010010001001001000100100100100100010010010010010010010010010010010010010010010010010010010010010010010010001000100030303333233332333233323332333233233323323332333233233323332332333233323332333233323332333233323332333233332333233332333323333233332333323333233333233333233333323333332333333323333333333323333333333333333333333333333303333333303333303333033303330330330330303303030303030300300030000000000010001001001001010010101010010101010101011010101010110101011010110101101011011010110110110110110110110110110110110110111011011011101101110111011011101110111011101110111011110111011101111011110111101110111110111101111011111011110111110111110111110111110111101111101111011101101101000303030303033030330303303303033033030330330303303303033033030330330303303033033030330303303033030330303303033030330303033030303303030330303033030303033030303033030303030303033030303030303030303030303030303030303030303003030303030030303030030300303003030030300300300300300300300030000300003000000000100001001001010101010110110111011110111111111111111112111111121111211112111121112111211121121112112111211211211121121121121121121121211211211211212112112112121121211211212112121121211212112121121212112121121212112121211212121121212121121212121121212121121212121211212121212112121212121121212121121212112121121111100100000000003000000003000000003000000030000000300000000300000000030000000000030000000000000003000000000000000000000000000000000000000001000000000000100000000010000001000000100000100001000010001000100010001001001000100101001001010010101010101010110101101110111111111211121121121212112", 805, 512 ) );

  bool ok = true;
  std::vector<FreemanChain>::const_iterator it = theChains.begin();
  for (std::vector<FreemanChain>::const_iterator itEnd = theChains.end(); it != itEnd; ++it) 
    {
      #ifdef EXTRA_OUTPUT_DMLP
      cerr << "==========================================================" << endl;
      cerr << "==========================================================" << endl;
      cerr << "==========================================================" << endl;
      cerr << "==========================================================" << endl;
      cerr << "\n" << *it << std::endl;
      #endif
      ok = ok && checkFlipPosition( *it );
      if (!ok) break;
    }
  if (ok) 
    { 
      cout << "All test passed" << endl;
    } 
  else 
    {
      cout << "Problem with FreemanChain :\n";
      cout << (*it).chain << endl;
    }
    return (ok) ? EXIT_SUCCESS : 1;
}
