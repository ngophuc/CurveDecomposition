///////////////////////////////////////////////////////////////////////////////
// Test dynamic minimum length polygon
///////////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <vector>
#include <list>
#include "ImaGene/base/Arguments.h"
#include "ImaGene/base/StandardArguments.h"
#include "ImaGene/dgeometry2d/CurveCode.h"
#include "ImaGene/dgeometry2d/FreemanChain.h"
#include "ImaGene/dgeometry2d/GridCurve.h"

using namespace std;
using namespace ImaGene;


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//
// M A I N
//
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
static Arguments args;

int
main( int argc, char** argv ) 
{
  // -------------------------------------------------------------------------
  // Prepare arguments.
  // StandardArguments::addDigitalArgs( args, 2, false, false );
  StandardArguments::addIOArgs( args, true, false );
  args.addOption( "-orientation", "-orientation <CW,CCW,COMP>: gives the orientation of the input Freeman chain. COMP means that the orientation is given by the number of cw or ccw loops made by the contour. Default is COMP.", "COMP" );
  args.addOption( "-contour_type", "-contour_type <INNER,OUTER,INTERPIXEL,RPR>: indicates whether the input freeman chain codes the inner border of a digital object made of inside pixels (INNER), its outer border made of outside pixels (OUTER) or its interpixel digital contour (INTERPIXEL), RPR: the file contains a RPR. Default is INNER.", "INNER" );
  args.addOption( "-splitEdge", "-splitEdge <U|D|R|L|pos> <edge_index>: splits the edge of given index and according to the specified type Up, Down, Right, Left, or splits this edge at position pos.", "U", "0" );
  args.addOption( "-flipEdge", "-flipEdge <pos> <edge_index> <0/1>: flips the edge of given index at position pos, 0: flip outside pixel in, 1: flip inside pixel out.", "0", "0", "0" );
  args.addOption( "-mergeEdge", "-mergeEdge <edge_index>: tries to merge the edge of given index with its successor.", "0" );
  args.addOption( "-simplifyEdge", "-simplifyEdge <edge_index>: tries to simplify the edge of given index with its successor.", "0" );
  args.addOption( "-simplifyZone", "-simplifyZone <edge_index_begin> <edge_index_end>: tries to simplify the specified zone.", "0", "1" );
  args.addOption( "-flip", "-flip <position> <inside=0/1> <simplify=0/1>: flips the curve at the visitor position [pos], 0: flip outside pixel in, 1: flip inside pixel out.", "0", "0", "1" );
  args.addBooleanOption( "-undoFlip", "-undoFlip: undo the flipEdge." );
  args.addBooleanOption( "-dFC", "-dFC: displays the freeman chain code associated with the current grid curve." );
  args.addBooleanOption( "-dV", "-dV: displays the visitor movement in the current grid curve." );

  if ( ( argc <= 0 ) 
       || ! args.readArguments( argc, argv ) ) 
    {
      cerr << args.usage( "test_GridCurve", 
			  "Reads a freeman chaincode. Tests dynamic minimum length polygon.",
			  "" )
	   << endl;
      return 1;
    }

  // -------------------------------------------------------------------------
  // Define curve.
  GridCurve curve;
  ReversiblePolygonalRepresentation rpr;
  // Read Freeman chain and creates space/contour.
  istream & in_str = StandardArguments::openInput( args );
  string contour_type = args.getOption( "-contour_type" )->getValue( 0 );
  if ( contour_type == "RPR" )
    {
      bool ok = rpr.read( in_str );
      if ( ! ok ) 
	{
	  cerr << "Error reading RPR." << endl;
	  return 3;
	}
      curve.initFromReversiblePolygonalRepresentation( rpr );
    }
  else
    {
      FreemanChain c;
      FreemanChain::read( in_str, c );
      if ( ! in_str.good() )
	{
	  cerr << "Error reading Freeman chain code." << endl;
	  return 2;
	}
      // Define orientation.
      string str_orientation = args.getOption( "-orientation" )->getValue( 0 );
      bool cw = true;
      int nb_ccw_loops = c.isClosed();
      cerr << "Number of CCW loops = " << nb_ccw_loops << endl;
      if ( str_orientation == "CW" )       cw = true;
      else if ( str_orientation == "CCW" ) cw = false;
      else                                 cw = nb_ccw_loops < 0;

      if ( contour_type == "INNER" )
	curve.initFromFreemanChain( c.begin(), c.end(), true, cw );
      else if ( contour_type == "OUTER" )
	curve.initFromFreemanChain( c.begin(), c.end(), false, cw );
      else if ( contour_type == "INTERPIXEL" )
	{
	  CurveCode cc;
	  cc.init( Vector2i( c.x0, c.y0 ), cw, c.chain );
	  std::cerr << cc << std::endl;
	  
	  // Builds the naive RPR.
	  // rpr.initFromInterpixelCurveCode( cc, true );

	  rpr.cw = cw;
	  std::vector<uint> vi;
	  bool closed =
	    FreemanChain::computeMLP( rpr.vx, rpr.vy, vi, rpr.vt,
	  			      rpr.twice_dv, c, cw )
	    != c.end();
	  if ( closed )
	    {
	      rpr.vx.push_back( rpr.vx[ 0 ] );
	      rpr.vy.push_back( rpr.vy[ 0 ] );
	      rpr.vt.push_back( rpr.vt[ 0 ] );
	    }
	  cerr << "Freeman chain is " << ( closed ? "closed" : "open" ) << "."
	       << endl;
	  for ( unsigned int i = 0; i < rpr.vx.size(); ++i )
	    cerr << "(" << rpr.vx[i] << "," << rpr.vy[i] << ","
		 << rpr.vt[i] << ") ";
	  cerr  << endl;
	  //rpr.purge();
	  curve.initFromReversiblePolygonalRepresentation( rpr );
	}
      unsigned nb_found = 0;
      unsigned nb = 0;
      cerr << "------ Visitor list -------"  << endl;
      GridCurve::Visitor visitor = curve.begin();
      GridCurve::Visitor visitor_begin = visitor;
      do { 
	cerr << (*visitor).first.x() << ' '
	     << (*visitor).first.y() << ' ' 
	     << (char) ('0'+(*visitor).second) << endl;
	++visitor;
      } while ( visitor != visitor_begin );
      cerr << "------ Curve list -------"  << endl;
      for ( FreemanChain::const_iterator it = c.begin(), it_end = c.end();
	    it != it_end; ++it )
	{
	  cerr << (*it).x() << ' ' << (*it).y() << ' ' << it.getCode() << endl;
	  // visitor = 
	  //   curve.findPointel( (*it).x(), (*it).y(), it.getCode(),
	  // 		       curve.begin() );
	  // if ( visitor != curve.end() ) 
	  //   {
	  //     GridCurve::Visitor::Value val = *visitor;
	  //     if ( ( val.first == *it ) && ( val.second == it.getCode() ) )
	  // 	++nb_found;
	  //     else
	  // 	cerr << "FindPointel: not found " 
	  // 	     << *it << " " << it.getCode() << endl;
	  //   }
	  // ++nb;
	}
      cerr << "GridCurve::findPointel: found " << nb_found << "/" << nb << endl;
    }
  // cout << "freeman   = " << c << endl;
  // cout << "gridcurve = " << curve << endl;
  // curve.verboseSelfDisplay( cout );
  // cout << endl;
  if ( args.check( "-splitEdge" ) )
    {
      int edge_index = args.getOption( "-splitEdge" )->getIntValue( 1 );
      string split_type =  args.getOption( "-splitEdge" )->getValue( 0 );
      typedef GridCurve::Edge Edge;
      typedef GridCurve::EdgeListIterator EdgeListIterator;
      EdgeListIterator it = curve.beginEdge();
      while ( ( edge_index > 0 ) && ( it != curve.endEdge() ) )
	{
	  ++it; --edge_index;
	}
      if ( it != curve.endEdge() )
	{
	  EdgeListIterator it_dummy;
	  cerr << "edge   = " << *it << endl;
	  if ( split_type == "U" )
	    curve.splitUp( it, it_dummy );
	  else if ( split_type == "D" )
	    curve.splitDown( it, it_dummy );
	  else if ( split_type == "L" )
	    curve.splitLeft( it, it_dummy );
	  else if ( split_type == "R" )
	    curve.splitRight( it, it_dummy );
	  else
	    {
	      unsigned int position = atoi( split_type.c_str() );
	      curve.splitEdgeAt( it, position, true ); // split inside
	    }
	}
    }
  if ( args.check( "-flipEdge" ) )
    {
      unsigned int position = args.getOption( "-flipEdge" )->getIntValue( 0 );
      int edge_index = args.getOption( "-flipEdge" )->getIntValue( 1 );
      bool inside = args.getOption( "-flipEdge" )->getIntValue( 2 );
      typedef GridCurve::Edge Edge;
      typedef GridCurve::Visitor Visitor;
      typedef GridCurve::EdgeListIterator EdgeListIterator;
      // find edge.
      EdgeListIterator it = curve.beginEdge();
      while ( ( edge_index > 0 ) && ( it != curve.endEdge() ) )
	{
	  ++it; --edge_index;
	}
      // determine position within edge.
      Visitor v( curve, it, position );
      Visitor b( v ); 
      Visitor e( v );
      // do the flip.
      curve.flip( v, inside, b, e, true );
    }
  if ( args.check( "-mergeEdge" ) )
    {
      int edge_index = args.getOption( "-mergeEdge" )->getIntValue( 0 );
      typedef GridCurve::Edge Edge;
      typedef GridCurve::EdgeListIterator EdgeListIterator;
      EdgeListIterator it = curve.beginEdge();
      while ( ( edge_index > 0 ) && ( it != curve.endEdge() ) )
	{
	  ++it; --edge_index;
	}
      if ( it != curve.endEdge() )
	{
	  const Edge & edge = *it;
	  EdgeListIterator it_next = curve.nextEdge( it );
	  cerr << "edge   = " << edge << "," << *it_next << endl;
	  bool merge = curve.merge( it );
	  if ( merge )
	    cerr << "merge  = " << *it << endl;
	  else 
	    cerr << "no merge." << endl;
	}
    }
  if ( args.check( "-simplifyEdge" ) )
    {
      int edge_index = args.getOption( "-simplifyEdge" )->getIntValue( 0 );
      typedef GridCurve::Edge Edge;
      typedef GridCurve::EdgeListIterator EdgeListIterator;
      EdgeListIterator it = curve.beginEdge();
      while ( ( edge_index > 0 ) && ( it != curve.endEdge() ) )
	{
	  ++it; --edge_index;
	}
      if ( it != curve.endEdge() )
	{
	  const Edge & edge = *it;
	  EdgeListIterator it_next = curve.nextEdge( it );
	  cerr << "edge   = " << edge << "," << *it_next << endl;
	  EdgeListIterator it_bf, it_al;
	  EdgeListIterator min, max;
	  bool simplify = curve.preSimplifyAt( it, it_bf, it_al, min, max );
	  if ( simplify )
	    {
	      cerr << "simplify  = ";
	      for ( it_bf = curve.nextEdge( it_bf ); 
		    it_bf != it_al; 
		    it_bf = curve.nextEdge( it_bf ) )
		cerr << "," << *it_bf;
	      cerr << endl;
	    }
	  else 
	    cerr << "no simplify." << endl;
	}
    }
  typedef GridCurve::Edge Edge;
  typedef GridCurve::EdgeListIterator EdgeListIterator;
  if ( args.check( "-simplifyZone" ) )
    {
      int edge_index_b = args.getOption( "-simplifyZone" )->getIntValue( 0 );
      int edge_index_e = args.getOption( "-simplifyZone" )->getIntValue( 1 );
      EdgeListIterator it = curve.beginEdge();
      EdgeListIterator itb = curve.endEdge();
      EdgeListIterator ite = curve.endEdge();
      EdgeListIterator itmid = curve.beginEdge();
      for ( int i = 0; it != curve.endEdge(); ++it, ++i )
	{
	  if ( i == edge_index_b ) itb = it;
	  if ( i == edge_index_e ) ite = it;
	  if ( i & 1 ) ++itmid;
	}
      if ( ite == curve.endEdge() )
	{ // loop
	  bool modif = true;
	  while ( modif )
	    {
	      cerr << "LOOP: first simplification" << endl;
	      ite = curve.previousEdge( curve.beginEdge() );
	      modif = curve.simplifyZone( itb, ite );
	      cerr << "LOOP: second simplification" << endl;
	      itb = curve.beginEdge();
	      ite = curve.previousEdge( curve.beginEdge() );
	      bool modif2 = curve.simplifyZone( ite, itb );
	      modif = modif || modif2;
	    }
	}
      else
	{
	  cerr << "LOOP: unique simplification" << endl;
	  curve.simplifyZone( itb, ite );
	}
    }

  if ( args.check( "-flip" ) )
    {
      unsigned int visitor_pos = args.getOption( "-flip" )->getIntValue( 0 );
      bool inside = args.getOption( "-flip" )->getIntValue( 1 );
      bool simplify = args.getOption( "-flip" )->getIntValue( 2 );
      typedef GridCurve::Edge Edge;
      typedef GridCurve::Visitor Visitor;
      typedef GridCurve::EdgeListIterator EdgeListIterator;
      Visitor v = curve.begin();
      Visitor endv = curve.end();
      for ( unsigned int i = 0; i < visitor_pos; ++i )
	{
	  ++v;
	  // if ( v == endv ) v = curve.begin();
	}
      Visitor b( v ); 
      Visitor e( v );
      curve.flip( v, inside, b, e, true );
      EdgeListIterator begin_mod =  b.edgeIterator();
      // ( b == endv ) ? curve.previousEdge( curve.endEdge() ) : b.edgeIterator();
      EdgeListIterator end_mod = e.edgeIterator();
      //( e == endv ) ? curve.nextEdge( curve.endEdge() ) : e.edgeIterator();
      if ( simplify )
	curve.simplifyZone( begin_mod, end_mod );
    }


  cout << "mylist.append( " << curve << " )" << endl;

  if ( args.check( "-dFC" ) )
    {
      vector<unsigned char> v;
      EdgeListIterator itb = curve.beginEdge();
      EdgeListIterator it = itb;
      EdgeListIterator ite = curve.endEdge();
      for ( ; it != ite; ++it )
	{
	  it->pushBackFreemanCode( v );
	}
      GridCurve::simplifyFreemanCode( v );
      cout << itb->pointel.x() << " " << itb->pointel.y() << " ";
      for ( unsigned int i = 0; i < v.size(); ++i )
	cout << (char) ( '0' + v[ i ] );
      cout << endl;
    }

  if ( args.check( "-dV" ) )
    {
      for ( GridCurve::Visitor it = curve.begin(), it_end = curve.end();
	    it != it_end; ++it )
	{
	  GridCurve::Visitor::Value v = *it;
	  cout << (char) ( '0' + v.second ) << " " << v.first.x()
	       << " " << v.first.y() << endl;
	}
    }

  if ( args.check( "-undoFlip" ) )
    {
      curve.undoFlip();
      cout << "undocurve = " << curve << endl;
    }
  return 0;
}
