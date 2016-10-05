///////////////////////////////////////////////////////////////////////////////
// Computes incidence matrices of a cubical complex from a binary image
///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <vector>
#include "ImaGene/base/Arguments.h"
#include "ImaGene/base/StandardArguments.h"
#include "ImaGene/base/BasicTypes.h"
#include "ImaGene/base/HashTable.h"
#include "ImaGene/timetools/Clock.h"
#include "ImaGene/digitalnD/KnSpace.h"
#include "ImaGene/digitalnD/KnShapes.h"
#include "ImaGene/digitalnD/KnSpaceCoordScanner.h"
#include "ImaGene/digitalnD/CubicalComplex.h"

using namespace std;
using namespace ImaGene;

namespace ImaGene
{
  class IncMatrix 
  {
  private:
    uint m_rows;
    uint m_columns;
    HashTable<int> m_values;
    
    public:
    IncMatrix()
    {
    }
    ~IncMatrix()
    {
    }
    void init( uint nb_rows, uint nb_columns )
    {
      m_values.init( nb_rows * nb_columns, nb_rows + nb_columns );
      m_rows = nb_rows;
      m_columns = nb_columns;
    }
    void setValue( uint row, uint column, int value )
    {
      uint key = row * m_columns + column;
      m_values.put( key, value );
    }
    int getValue( uint row, uint column ) const
    {
      int value;
      uint key = row * m_columns + column;
      if ( m_values.get( key, value ) )
	return value;
      else
	return 0;
    }
    

  };
  
  

  
  class IMCubicalComplex
  {
  private:
    
    /**
     * The digital space.
     */
    const KnSpace & m_ks;
    
    /**
     * The dimension of the maximal cells.
     */
    uint m_n;

  public:

    /**
     * For each dimension 0..n, the set of (unsigned cells).
     * @TODO for now public: Just for quick test. 
     */
    std::vector<Kn_uid>* m_cells;

    /**
     * For each dimension 0..n-1, the set of (unsigned cells).
     * @TODO for now public: Just for quick test. 
     */
    IncMatrix* m_imatrices;
    
  public:
    
    /**
     * Constructs an empty cubical complex.
     * 
     * @param ks the digital space in which it is embedded.
     */
    IMCubicalComplex( const KnSpace & ks )
      : m_ks( ks ), m_n( 0 ), m_cells( 0 ), m_imatrices( 0 )
    {
    }
    
    /**
     * Destructor.
     */
    ~IMCubicalComplex()
    {
      reset();
    }

    /**
     * Empties the cubical complex.
     */
    void reset()
    {
      m_n = 0;
      if ( m_cells != 0 )
	delete[] m_cells;
      if ( m_imatrices != 0 )
	delete[] m_imatrices;
    }
    
    /** 
     * Constructs an homogeneous unrestricted cubical complex from the
     * given cells.  All the cells are stored in the array
     * [m_cells]. Incidences are stored in the array of "matrices"
     * [m_imatrices].
     * 
     * @param object a set of unsigned cells of dimension [n].
     * @param n the dimension of the cells
     */
    void uinit( KnCharSet object, uint n );

    /** 
     * Constructs an homogeneous unrestricted cubical complex from the
     * given cells.  All the cells are stored in the array
     * [m_cells]. Incidences are stored in the array of "matrices"
     * [m_imatrices].
     * 
     * @param object a set of signed cells of dimension [n].
     * @param n the dimension of the cells
     */
    void sinit( KnRCellSet object, uint n );

    /** 
     * Displays the incidence matrix of dimension [m] as a set of m-1 chains.
     * 
     * @param m the dimension of the incidence matrix.
     */
    void display( uint m ) const;
    
    
  };
  
}


/** 
 * Constructs an homogeneous unrestricted cubical complex from the
 * given cells.  All the cells are stored in the array
 * [m_cells]. Incidences are stored in the array of "matrices"
 * [m_imatrices].
 * 
 * @param object a set of signed cells of dimension [n].
 * @param n the dimension of the cells
 */
void
ImaGene::IMCubicalComplex::sinit( KnRCellSet object, uint n )
{
  KnCharSet uobject = KnCharSet::create( m_ks, n, false, 0 );
  for ( KnRCellSet::cell_iterator p = object.begin(); 
	p != object.end();
	++p )
    uobject += m_ks.unsigns( *p );
  uinit( uobject, n );
}



/** 
 * Constructs an homogeneous unrestricted cubical complex from the
 * given cells.  All the cells are stored in the array
 * [m_cells]. Incidences are stored in the array of "matrices"
 * [m_imatrices].
 * 
 * @param object a set of unsigned cells of dimension [n].
 * @param n the dimension of the cells
 */
void
ImaGene::IMCubicalComplex::uinit( KnCharSet object, uint n )
{
  reset();
  m_n = n;
  m_cells = new std::vector<Kn_uid>[ n + 1 ];
  m_imatrices = new IncMatrix[ n ];

  KnCharSet all_cells = KnCharSet::create( m_ks, false, 0 );
  // 1st pass: determine all maximal cells
  KnCharSet::cell_iterator p = object.begin();
  KnCharSet::cell_iterator p_end = object.end();
  for ( ; p != p_end; ++p )
    {
      Kn_uid c = *p;       // current cell
      if ( m_ks.udim( c ) == n )
	{
	  m_cells[ n ].push_back( c );
	  all_cells += c;
	}
    }
  // 2nd pass: for each dimension, compute incident cells then incidence.
  for ( uint m = n; m != 0; --m )
    {
      uint j = 0;
      HashTable<uint> index;
      index.init( m_ks.ulast(), m_cells[ m ].size() * m * 2 );
      for ( uint i = 0; i < m_cells[ m ].size(); ++i )
	{
	  Kn_uid c = m_cells[ m ][ i ];
	  for ( KnSpace::dir_iterator q = m_ks.ubegin_dirs( c );
		! q.end();
		++q )
	    {
	      Kn_uid s1 = m_ks.uincident( c, *q, false );
	      Kn_uid s2 = m_ks.uincident( c, *q, true );
	      if ( ! all_cells[ s1 ] )
		{
		  m_cells[ m - 1 ].push_back( s1 );
		  index[ s1 ] = j++;
		  all_cells += s1;
		}
	      if ( ! all_cells[ s2 ] )
		{
		  m_cells[ m - 1 ].push_back( s2 );
		  index[ s2 ] = j++;
		  all_cells += s2;
		}
	    }
	}
      // Now compute incidence.
      m_imatrices[ m - 1 ].init( m_cells[ m ].size(), // rows
				 m_cells[ m - 1 ].size() // columns
				 );
      for ( uint i = 0; i < m_cells[ m ].size(); ++i )
	{
	  Kn_uid c = m_cells[ m ][ i ];
	  Kn_sid sc = m_ks.signsPos( c );
	  for ( KnSpace::dir_iterator q = m_ks.ubegin_dirs( c );
		! q.end();
		++q )
	    {
	      Kn_sid ss1 = m_ks.sincident( sc, *q, false );
	      Kn_sid ss2 = m_ks.sincident( sc, *q, true );
	      Kn_uid s1 = m_ks.unsigns( ss1 );
	      Kn_uid s2 = m_ks.unsigns( ss2 );
	      Kn_sign sign_ss1 = m_ks.decodeSign( ss1 );
	      Kn_sign sign_ss2 = m_ks.decodeSign( ss2 );
	      m_imatrices[ m - 1 ].setValue
		( i, index[ s1 ], ( sign_ss1 == KnTypes::POS ) ? 1 : -1 );
	      m_imatrices[ m - 1 ].setValue
		( i, index[ s2 ], ( sign_ss2 == KnTypes::POS ) ? 1 : -1 );
	    }
	}
    }
  
}


/** 
 * Displays the incidence matrix of dimension [m] as a set of m-1 chains.
 * 
 * @param m the dimension of the incidence matrix.
 */
void
ImaGene::IMCubicalComplex::display( uint m ) const
{
  for ( uint i = 0; i < m_cells[ m+1 ].size(); ++i )
    {
      cout << "Bord( [" << i << "|" << m_cells[ m+1 ][ i ] << "] ) = ";
      for ( uint j = 0; j < m_cells[ m ].size(); ++j )
	{
	  int inc = m_imatrices[ m ].getValue( i, j );
	  if ( inc == 1 ) cout << "+";
	  else if ( inc == -1 ) cout << "-";
	  if ( inc != 0 )
	    cout << "[" << j << "|" << m_cells[ m ][ j ] << "] ";
	}
      cout << endl;
    }
}




KnCharSet
makeSphere( KnSpace & ks, 
	    float x0, float y0, float z0, 
	    float rsphere )
{
  Kn_size center[ 3 ];
  center[ 0 ] = 2*( (Kn_size) floor( x0 + 0.5) ) + 1;
  center[ 1 ] = 2*( (Kn_size) floor( y0 + 0.5) ) + 1;
  center[ 2 ] = 2*( (Kn_size) floor( z0 + 0.5) ) + 1;
  Kn_uid clow = ks.uspel( ks.ukcode( center ) );
  
  cout << "--- creating sphere (r=" << rsphere << ")... ";
  cout.flush();
  Clock::startClock();
  KnCharSet sph1 = KnShapes::umakeVolumicSphere( ks, clow, rsphere );
  long ti2 = Clock::stopClock();
  cout << "in " << ti2 << " ms." << endl;
  return sph1;
}


KnCharSet
makeCubeMinusSphere( KnSpace & ks, 
		     float x0, float y0, float z0, 
		     int rc, float rsphere )
{
  Kn_size center[ 3 ];
  center[ 0 ] = 2*( (Kn_size) floor( x0 + 0.5) ) + 1;
  center[ 1 ] = 2*( (Kn_size) floor( y0 + 0.5) ) + 1;
  center[ 2 ] = 2*( (Kn_size) floor( z0 + 0.5) ) + 1;
  Kn_uid ccenter = ks.uspel( ks.ukcode( center ) );
  
  cout << "--- creating sphere (r=" << rsphere << ")... ";
  cout.flush();
  Clock::startClock();
  Kn_uid low_uid1 = KnShapes::ugetCubeLowerBound( ks, ccenter, rc );
  Kn_uid up_uid1 = KnShapes::ugetCubeUpperBound( ks, ccenter, rc );
  
  KnCharSet sph1 = KnShapes::umakeVolumicSphere( ks, ccenter, rsphere );
  cout << "parallepiped (r=" << rc << ")... " << endl;
  KnCharSet par1 = KnShapes::umakeVolumicParallelepiped( ks,
   							 low_uid1,
   							 up_uid1 );
  cout << " parallepiped - sphere... " << endl;
  KnCharSet voxset = par1 - sph1;
  long ti2 = Clock::stopClock();
  cout << "in " << ti2 << " ms." << endl;
  return voxset;
}


KnCharSet
makeObject( KnSpace & ks )
{
  KnCharSet voxset = KnCharSet::create( ks, ks.dim(), false, 0 );
  // first and last pixel/voxel/spel of the image.
  Kn_uid first = ks.ufirstCell( ks.dim() );
  Kn_uid last = ks.ulastCell( ks.dim() );
  KnSpaceCoordScanner scanner( ks, first, last );
  Kn_uid p = scanner.lower_bound;
  Vector vp = ks.ucentroid( p );
  do 
    { // ... whatever [p] is the current cell
      // ... and [vp] its centroid
      if ( ( vp.ro( 0 ) >= 5.0 ) && ( vp.ro( 0 ) <= 8.0 ) 
	   && ( vp.ro( 1 ) >= 3.0 ) && ( vp.ro( 1 ) <= 4.0 ) 
	   && ( ( ks.dim() <= 2 )
		|| ( ( vp.ro( 2 ) >= 4.0 ) && ( vp.ro( 2 ) <= 6.0 ) )
		) )
	voxset += p;
    }
  while ( scanner.increment( p, vp ) ); 
  cout << "OBJECT done" << endl;
  return voxset;
}

/**
 * Doit sans doute etre modifie pour lister differement les faces (dim=2).
 */
void
outputGeometry( const KnSpace & ks, Kn_uid cell, ostream & out );


/**
 * Doit sans doute etre modifie pour lister differement les faces (dim=2).
 */
void
outputGeometryOfFace( const KnSpace & ks, Kn_uid cell, ostream & out )
{
  uint dimcell = ks.udim( cell );
  if ( dimcell == 2 )
    {
      KnSpace::dir_iterator q = ks.ubegin_dirs( cell );
      uint d1 = *q; ++q;
      uint d2 = *q;
      Kn_uid edge1 = ks.uincident( cell, d1, true );
      Kn_uid edge2 = ks.uincident( cell, d1, false );
      Kn_uid vertex1 = ks.uincident( edge1, d2, true );
      Kn_uid vertex2 = ks.uincident( edge1, d2, false );
      Kn_uid vertex3 = ks.uincident( edge2, d2, true );
      Kn_uid vertex4 = ks.uincident( edge2, d2, false );
      outputGeometry( ks, vertex1, out );
      outputGeometry( ks, vertex2, out );
      outputGeometry( ks, vertex4, out );
      outputGeometry( ks, vertex3, out );
    }
}


/**
 * Doit sans doute etre modifie pour lister differement les faces (dim=2).
 */
void
outputGeometry( const KnSpace & ks, Kn_uid cell, ostream & out )
{
  uint dimcell = ks.udim( cell );
  out << "( ";
  if ( dimcell == 0 )
    {
      // lists coordinates
      for ( uint i = 0; i < ks.dim(); ++i )
	out << ks.udecodeCoord( cell, i ) << " ";
    }
  else if ( dimcell >= 1 )
    {
      for ( KnSpace::dir_iterator q = ks.ubegin_dirs( cell );
	    ! q.end();
	    ++q )
	{
	  outputGeometry( ks, 
			  ks.uincident( cell, *q, false ),
			  out );
	  out << " ";
	  outputGeometry( ks, 
			  ks.uincident( cell, *q, true ),
			  out );
	  out << " ";
	}
    }
  out << ")";
}




// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

int main( int argc, char** argv )
{
  Arguments args;
  args.addOption( "-sphere", "-sphere <r>: the shape is a sphere of radius r.", "10" );
  args.addOption( "-cms", "-cms <rc> <rs>: the shape is a cube of edge rc minus a sphere of radius rs.", "10", "13.0" );
  args.addBooleanOption( "-object", "-object: the shape is an object as defined in the code." );
  args.addBooleanOption( "-cubic", "-cubic: tests classes (H)CubicalComplex and CubicalCell." );
  if ( ( argc <= 1 ) 
       || ! args.readArguments( argc, argv ) ) 
    {
      cout << args.usage( "test_IMCubicalComplex", 
			  "Computes incidence matrices for a given cubical complex.",
			  "-d -h -x -y -z -sphere -cms -object -cubic" ) 
	   << endl;
      return 1;
    }

  // dimension, sizes, radius
  uint d = StandardArguments::dim( args );
  if ( d > 3 )
    {
      cerr << "Dimension is 2 or 3." << endl;
      return 2;
    }
  
  Kn_size sizes[ d ];
  StandardArguments::fillSizes( args, sizes );

  // option "-sphere"
  bool sphere = args.check( "-sphere" );
  // option "-cms"
  bool cms = args.check( "-cms" );
  // option "-object"
  bool object = args.check( "-object" );

  // Defines space.
  KnSpace ks( d, sizes ); // Z3
    
  // center
  cout << "--- Space: " << ks << endl;
  KnCharSet voxset = KnCharSet::create( ks, d, false, 0 );
  if ( sphere )
    {
      float rsphere = args.getOption( "-sphere" )->getFloatValue( 0 );
      voxset = makeSphere( ks, 
			   sizes[ 0 ] / 2.0,
			   sizes[ 1 ] / 2.0,
			   sizes[ 2 ] / 2.0,
			   rsphere );
    }
  else if ( cms )
    {
      int rc = args.getOption( "-cms" )->getIntValue( 0 );
      float rsphere = args.getOption( "-cms" )->getFloatValue( 1 );
      voxset = makeCubeMinusSphere( ks, 
				    sizes[ 0 ] / 2.0,
				    sizes[ 1 ] / 2.0,
				    sizes[ 2 ] / 2.0,
				    rc, rsphere );
    }
  else if ( object )
    voxset = makeObject( ks );
  

  if ( args.check( "-cubic" ) )
    {
      CubicalCell::testCubicalCell();
      cout << endl;
      CubicalComplex::testCubicalComplex();
      HierarchicalCubicalComplex::test();
      return 0;
    }

  // Computes incidence.
  IMCubicalComplex k( ks );

  cout << " ---- VOLUME ---------" << endl;
  k.uinit( voxset, d );
  for ( uint m = 0; m < d; ++m )
    {
      cout << "..... dim = " << m << " .........." << endl;
      k.display( m );
    }

  cout << " ---- GEOMETRY ---------" << endl;
  for ( uint m = 0; m < d; ++m )
    {
      cout << "..... dim = " << m << " .........." << endl;
      for ( uint i = 0; i < k.m_cells[ m ].size(); ++i )
	{
	  cout << "Geom[" << i << "|" << k.m_cells[ m ][ i ] << "] = ";
	  outputGeometry( ks, k.m_cells[ m ][ i ], cout );
	  cout << endl;
	}
    }

  cout << " ---- SURFACE ---------" << endl;
  KnRCellSet surfelset = KnShapes::smakeBoundary( ks, voxset );
  k.sinit( surfelset, d - 1 );
  for ( uint m = 0; m < ( d - 1 ); ++m )
    {
      cout << "..... dim = " << m << " .........." << endl;
      k.display( m );
    }

  return 0;
}
