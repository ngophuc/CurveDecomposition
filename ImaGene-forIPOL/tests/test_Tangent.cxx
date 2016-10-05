///////////////////////////////////////////////////////////////////////////////
// Test module for tangent computation along contours in 2D.
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cmath>
#include <cstdlib>

//#include "LinAlg/LinAlg2D/Vector2D.hpp"
//#include "ImageLib/Gauss/G.hpp"
//#include "MINIWIN/RGBWindow.h"
// #include "ImaGeneUtils/dbgvisu2d/K2SpaceViewer.h"


#include "ImaGene/base/Proxy.h"
#include "ImaGene/base/Vector.h"
#include "ImaGene/base/VectorUtils.h"
#include "ImaGene/base/Arguments.h"
#include "ImaGene/base/StandardArguments.h"
#include "ImaGene/base/Shapes.h"
#include "ImaGene/timetools/Clock.h"
#include "ImaGene/mathutils/Mathutils.h"
#include "ImaGene/mathutils/G.h"
#include "ImaGene/mathutils/Functions.h"
#include "ImaGene/mathutils/Polynomial.h"
#include "ImaGene/dgeometry2d/C4CSegment.h"
#include "ImaGene/dgeometry2d/C4CIterator.h"
#include "ImaGene/dgeometry2d/C4CGeometry.h"
#include "ImaGene/dgeometry2d/C4CSegmentPencil.h"
#include "ImaGene/dgeometry2d/C4CTangentialCover.h"
#include "ImaGene/dgeometry2d/EuclideanGeometry.h"
#include "ImaGene/digitalnD/C4CIteratorOnBdry.h"
#include "ImaGene/digitalnD/KnCharSet.h"
#include "ImaGene/digitalnD/KnRCellSet.h"
#include "ImaGene/digitalnD/BelAdjacency.h"
#include "ImaGene/digitalnD/ObjectBoundary.h"
#include "ImaGene/digitalnD/KnShapes.h"
#include "ImaGene/digitalnD/Embedder.h"
#include "ImaGene/digitalnD/GridEmbedder.h"
#include "ImaGene/digitalnD/KnSpaceScanner.h"
#include "ImaGene/helper/C4CTangentialCoverGeometry.h"
#include "ImaGene/helper/GlobalC4CGeometry.h"
#include "ImaGene/helper/ShapeHelper.h"

using namespace std;
using namespace ImaGene;

static Arguments args;


float ang_from_v2d ( Vector2D & v2d ) {

  if ( v2d.x() == 0.0 ) 
    {
      if ( v2d.y() >0 )
	return M_PI/2.0; 
      else 
	return 3*M_PI/2.0; 
    }
  else if (  ( v2d.x() > 0.0 ) && (   v2d.y() >= 0.0 ) )
    return atan(v2d.y()/v2d.x());
  else if (  ( v2d.x() > 0.0 ) && (   v2d.y() <= 0.0 ) )
    return 2*M_PI + atan(v2d.y()/v2d.x());
  else if (  ( v2d.x() < 0.0 ) && (   v2d.y() >= 0.0 ) )
    return atan(v2d.y()/v2d.x()) + M_PI;
  else if (  ( v2d.x() < 0.0 ) && (   v2d.y() <= 0.0 ) )
    return  atan(v2d.y()/v2d.x()) + M_PI;
}

int get_quandrant_from_angle (float ang)
{
  while (ang > 2.0*M_PI) {
    ang -= 2.0*M_PI; }
  while (ang < 0) {
    ang += 2.0*M_PI; }
  if ( ang < M_PI/2.0 ) 
    return 0; 
  else if ( ang < M_PI ) 
    return 1; 
  else if ( ang < 3.0*M_PI/2.0 ) 
    return 2;
  else 
    return 3;
}

StarShaped* 
shapeFromArgs( KnSpace & ks, const Embedder & embedder,
	       KnCharSet & voxset,
	       Kn_sid & bel, uint & nb_bels )
{
  StarShaped* shape = ShapeHelper::makeStarShapedFromArgs( args );
  Vector2D dir( 1.0, 0.0 );
  cerr << "# --- digitizing shape " << endl;
  voxset = ShapeHelper::makeShape( ks, embedder, *shape );
  cerr << "# --- find start bel " << endl;
  bel = ShapeHelper::getBelOnStarShapedBoundary( ks, embedder, *shape, dir );
  cerr << "# --- compute contour size " << endl;
  nb_bels = KnShapes::sgetContourSize( ks, voxset, bel,
				       *( ks.sbegin_dirs( bel ) ) );
  
  //cerr << "Set=" << voxset << endl;
  //cerr << "Bel=";
  //ks.displayKn_sid( bel, cerr );
  //cerr << endl;
  return shape;

}


///////////////////////////////////////////////////////////////////////////////
// TANGENT COMPUTATION
///////////////////////////////////////////////////////////////////////////////

class TangentComputer
{
public:
  uint size;
  uint type;
  float grid_step;
  TangentComputer() 
    : size( 0 ), type(0), grid_step (0.0f)
  {
  }

  virtual ~TangentComputer () 
  {}
  
//   virtual bool init (const C4CIteratorOnSurface* iter)
//   {
//     return true;
//   }

  virtual void reset ()
  {};

  virtual void setType (uint t)
  {
    type = t;
  };

  virtual Vector2D computeTangent( KnSpace & ks,
				   const C4CIteratorOnSurface* iter,
				   const Embedder & embedder, 
				   float & window,
				   float & ds,
				   Vector2D & ext_front,
				   Vector2D & ext_back ) = 0;
};

// class designed to time test purposes only
class TestTimeTangentComputer : public TangentComputer
{
public:
  Vector2D computeTangent( KnSpace & ks,
			   const C4CIteratorOnSurface* iter,
			   const Embedder & embedder, 
			   float & window,
			   float & ds,
			   Vector2D & ext_front,
			   Vector2D & ext_back )
  {
    uint j =0;
    for (int i = 0; i < (int) size; ++i) {
      j=i;
    }
    return Vector2D(cos (1.0/j), sin(1.0/j));
  }
};


// Symmetric discrete tangent computer.
class SDTangentComputer : public TangentComputer
{
public:
  Vector2D computeTangent( KnSpace & ks,
			   const C4CIteratorOnSurface* iter,
			   const Embedder & embedder, 
			   float & window,
			   float & ds,
			   Vector2D & ext_front,
			   Vector2D & ext_back );
};

Vector2D
SDTangentComputer::computeTangent( KnSpace & ks,
				   const C4CIteratorOnSurface* iter,
				   const Embedder & embedder, 
				   float & window,
				   float & ds,
				   Vector2D & ext_front,
				   Vector2D & ext_back )
{
  // Computes symmetric tangent at surfel '*iter'.
  C4CIteratorOnSurface* iterfront = 
    dynamic_cast<C4CIteratorOnSurface*>( iter->clone() );
  C4CIteratorOnSurface* iterback = 
    dynamic_cast<C4CIteratorOnSurface*>( iter->clone() );

  C4CSegment segment = C4CGeometry::symmetricTangent( *iterfront, *iterback, 
						      this->size );
  Vector2i tgt_l = segment.getLine().tangent();
  Kn_sid bel = iter->current();
  uint orth = ks.sorthDir( bel );
  uint track = 1 - orth;
  bool orth_direct = ks.sdirect( bel, orth );
  bool track_direct = ks.sdirect( bel, track );
  float orth_sign = orth_direct ? 1.0 : -1.0;
  float track_sign = track_direct ? 1.0 : -1.0;
  Vector2D tgt;
  if ( orth == 1 )
    tgt = Vector2D( track_sign * tgt_l.x(), - orth_sign * tgt_l.y() );
      else
    tgt = Vector2D( - orth_sign * tgt_l.y(), track_sign * tgt_l.x() );
  VectorUtils::normalize( tgt );
  Vector2D pt_r;
  Vector2D pt_l;
  embedder.sembed( iterfront->current(), pt_r );
  embedder.sembed( iterback->current(), pt_l );
  ext_front = pt_r;
  ext_back = pt_l;

  pt_r -= pt_l;
  window = VectorUtils::norm( pt_r );
  //window = VectorUtils::distance( pt_r, pt_l );

  delete iterback;
  delete iterfront;

  // averaged length.
  ds = 1 / ( fabs( tgt.ro( 0 ) ) + fabs( tgt.ro( 1 ) ) );

  return tgt;
}


// Extended discrete tangent computer.
class EDTangentComputer : public TangentComputer
{
public:
  Vector2D computeTangent( KnSpace & ks,
			   const C4CIteratorOnSurface* iter,
			   const Embedder & embedder, 
			   float & window,
			   float & ds,
			   Vector2D & ext_front,
			   Vector2D & ext_back );
};

Vector2D
EDTangentComputer::computeTangent( KnSpace & ks,
				   const C4CIteratorOnSurface* iter,
				   const Embedder & embedder, 
				   float & window,
				   float & ds,
				   Vector2D & ext_front,
				   Vector2D & ext_back )
{
  // Computes symmetric tangent at surfel '*iter'.
  C4CIteratorOnSurface* iterfront = 
    dynamic_cast<C4CIteratorOnSurface*>( iter->clone() );
  C4CIteratorOnSurface* iterback = 
    dynamic_cast<C4CIteratorOnSurface*>( iter->clone() );
  C4CSegment segment = 
    C4CGeometry::extendedTangent( *iterfront, *iterback, 
				  this->size );
  Vector2i tgt_l = segment.getTangent();

  Kn_sid bel = iter->current();
  uint orth = ks.sorthDir( bel );
  uint track = 1 - orth;
  bool orth_direct = ks.sdirect( bel, orth );
  bool track_direct = ks.sdirect( bel, track );
  float orth_sign = orth_direct ? 1.0 : -1.0;
  float track_sign = track_direct ? 1.0 : -1.0;
  Vector2D tgt;
  if ( orth == 1 )
    tgt = Vector2D( track_sign * tgt_l.x(), - orth_sign * tgt_l.y() );
      else
    tgt = Vector2D( - orth_sign * tgt_l.y(), track_sign * tgt_l.x() );
  VectorUtils::normalize( tgt );
  Vector2D pt_r;
  Vector2D pt_l;
  embedder.sembed( iterfront->current(), pt_r );
  embedder.sembed( iterback->current(), pt_l );
  ext_front = pt_r;
  ext_back = pt_l;

  // Limited window for TDE.
  Vector2D pt;
  embedder.sembed( iter->current(), pt );
  pt_l -= pt;
  float windowl = VectorUtils::norm( pt_l );
  pt_r -= pt;
  float windowr = VectorUtils::norm( pt_r );
  window = windowr <= windowl ? 2*windowr : 2*windowl;
//   pt_r -= pt_l;
//   window = norm( pt_r );

  delete iterback;
  delete iterfront;

  // averaged length.
  ds = 1 / ( fabs( tgt.ro( 0 ) ) + fabs( tgt.ro( 1 ) ) );

  return tgt;
}


// Half-tangents discrete tangent computer.
class HTangentComputer : public TangentComputer
{
public:
  Vector2D computeTangent( KnSpace & ks,
			   const C4CIteratorOnSurface* iter,
			   const Embedder & embedder, 
			   float & window,
			   float & ds,
			   Vector2D & ext_front,
			   Vector2D & ext_back );
};

Vector2D
HTangentComputer::computeTangent( KnSpace & ks,
				  const C4CIteratorOnSurface* iter,
				  const Embedder & embedder, 
				  float & window,
				  float & ds,
				  Vector2D & ext_front,
				  Vector2D & ext_back )
{
  // Computes symmetric tangent at surfel '*iter'.
  C4CIteratorOnSurface* iterfront = 
    dynamic_cast<C4CIteratorOnSurface*>( iter->clone() );
  C4CIteratorOnSurface* iterback = 
    dynamic_cast<C4CIteratorOnSurface*>( iter->clone() );
  C4CSegment segment_pos;
  segment_pos.init();
  C4CGeometry::longestPositiveSegment( segment_pos, *iterfront, 
				       this->size );
  C4CSegment segment_neg;
  segment_neg.init();
  C4CGeometry::longestNegativeSegment( segment_neg, *iterback, 
				       this->size );
  DLine med = DLine::medianLine( segment_pos.getLine(),
				 segment_neg.getLine() );
  Vector2i tgt_l = med.tangent();

  Kn_sid bel = iter->current();
  uint orth = ks.sorthDir( bel );
  uint track = 1 - orth;
  bool orth_direct = ks.sdirect( bel, orth );
  bool track_direct = ks.sdirect( bel, track );
  float orth_sign = orth_direct ? 1.0 : -1.0;
  float track_sign = track_direct ? 1.0 : -1.0;
  Vector2D tgt;
  if ( orth == 1 )
    tgt = Vector2D( track_sign * tgt_l.x(), - orth_sign * tgt_l.y() );
      else
    tgt = Vector2D( - orth_sign * tgt_l.y(), track_sign * tgt_l.x() );
  VectorUtils::normalize( tgt );
  Vector2D pt_r;
  Vector2D pt_l;
  embedder.sembed( iterfront->current(), pt_r );
  embedder.sembed( iterback->current(), pt_l );
  ext_front = pt_r;
  ext_back = pt_l;

  pt_r -= pt_l;
  window = VectorUtils::norm( pt_r );
  //window = VectorUtils::distance( pt_r, pt_l );

  delete iterback;
  delete iterfront;

  // averaged length.
  ds = 1 / ( fabs( tgt.ro( 0 ) ) + fabs( tgt.ro( 1 ) ) );

  return tgt;
}


// Tangents with maximal segments.
class MSTangentComputer : public TangentComputer
{
public:
  enum InterpolationFct {
    TRIANGLE, BELLSHAPE2 
  };
  

  ~MSTangentComputer();
  MSTangentComputer( InterpolationFct type_fct, bool using_analog );

  Vector2D computeTangent( KnSpace & ks,
			   const C4CIteratorOnSurface* iter,
			   const Embedder & embedder, 
			   float & window,
			   float & ds,
			   Vector2D & ext_front,
			   Vector2D & ext_back );

  /*
   * Not very nice... but only for prototyping.
   */
  static const uint m_max = 100;
  C4CSegment m_segments[ m_max ];
  R2RFunction* m_lambda;
  R2RFunction* m_lambda_p;
  bool m_using_analog;
  
};

MSTangentComputer::~MSTangentComputer()
{
  if ( m_lambda != 0 ) delete m_lambda;
  if ( m_lambda_p != 0 ) delete m_lambda_p;
}

MSTangentComputer::MSTangentComputer( InterpolationFct type_fct,
				      bool using_analog )
  : m_using_analog( using_analog )
{
  if ( type_fct == TRIANGLE )
    {
      m_lambda = new TriangleFunction;
      m_lambda_p = new DTriangleFunction;
    }
  else
    {
      m_lambda = new Polynomial( Polynomial::bellShape2( 1.0f ) );
      m_lambda_p = new Polynomial( Polynomial::bellShape2( 1.0f )
				   .derivative() );
    }
}



// Tangents with maximal segments.
Vector2D
MSTangentComputer::computeTangent( KnSpace & ks,
				   const C4CIteratorOnSurface* iter,
				   const Embedder & embedder, 
				   float & window,
				   float & ds,
				   Vector2D & ext_front,
				   Vector2D & ext_back )
{
  C4CIteratorOnSurface* it = 
    dynamic_cast<C4CIteratorOnSurface*>( iter->clone() );
  uint j = 0;
  uint k;
  uint m = m_max;
  if ( ! C4CGeometry::maximalSegments( *it, m_segments, j, k, m ) )
      cerr << "[AnalogComputerWithSegmentPencil::computeAnalog]"
	   << " Not enough segments." << endl;
  C4CSegmentPencil pencil( m_segments, j, k, m, *m_lambda, *m_lambda_p );
  float theta;
  if ( m_using_analog )
    {
      Vector2D dx = pencil.continuousAnalogDerivative( Vector2D(0.5f, 0.0f) );
      theta = atan( dx.ro( 1 ) / dx.ro( 0 ) );
    }
  else
    theta = pencil.angleToX( Vector2D( 0.5f, 0.0f ) );

  Kn_sid bel = iter->current();
  uint orth = ks.sorthDir( bel );
  uint track = 1 - orth;
  bool orth_direct = ks.sdirect( bel, orth );
  bool track_direct = ks.sdirect( bel, track );
  float orth_sign = orth_direct ? 1.0 : -1.0;
  float track_sign = track_direct ? 1.0 : -1.0;

  // Compute unitary tangent vector from angle to X-axis.
  Vector2D tgt;
  if ( orth == 1 )
    tgt = Vector2D( track_sign * cos( theta ), - orth_sign * sin( theta ) );
      else
    tgt = Vector2D( - orth_sign * sin( theta ), track_sign * cos( theta ) );

  // Embedding of leftmost and rightmost points are not valid.
  Vector2D pt_r;
  Vector2D pt_l;
  embedder.sembed( it->current(), pt_r );
  embedder.sembed( it->current(), pt_l );
  ext_front = pt_r;
  ext_back = pt_l;

  pt_r -= pt_l;
  window = VectorUtils::norm( pt_r );
  //window = VectorUtils::distance( pt_r, pt_l );

  // averaged length.
  ds = 1 / ( fabs( tgt.ro( 0 ) ) + fabs( tgt.ro( 1 ) ) );

  delete it;

  return tgt;
}

// Tangents based on maximal segments using pre-computation of 
// all the maximal segments prior to the computation of tangents.
// most of the code is based on MSTangentComputer
// restricted to analog=false

class PreMSTangentComputer : public TangentComputer
{
public:
  enum InterpolationFct {
    TRIANGLE, BELLSHAPE2 
  };
  

  ~PreMSTangentComputer();
  PreMSTangentComputer( InterpolationFct type_fct);

  bool init (KnSpace & ks, const C4CIteratorOnSurface* iter);

  void reset ();

  Vector2D computeTangent( KnSpace & ks,
			   const C4CIteratorOnSurface* iter,
			   const Embedder & embedder, 
			   float & window,
			   float & ds,
			   Vector2D & ext_front,
			   Vector2D & ext_back );

  /*
   * Not very nice... but only for prototyping.
   */
  uint m_cur_surf_idx;
  bool m_is_first_call_done;
  C4CTangentialCover m_t_cover; 
  C4CTangentialCoverGeometry m_t_cover_geometry;
  C4CTangentialCover::SurfelMaximalSegments m_pencil_on_cur_idx;  

  R2RFunction* m_lambda;
};

PreMSTangentComputer::~PreMSTangentComputer()
{
  if ( m_lambda != 0 ) delete m_lambda;
}

PreMSTangentComputer::PreMSTangentComputer( InterpolationFct type_fct )
					    
  :  m_cur_surf_idx (0),  m_is_first_call_done (false)
{
  if ( type_fct == TRIANGLE )
    {
      m_lambda = new TriangleFunction;
    }
  else
    {
      m_lambda = new Polynomial( Polynomial::bellShape2( 1.0f ) );
    }
}

void 
PreMSTangentComputer::reset ()
{
  m_is_first_call_done = false;
  m_cur_surf_idx = 0;
}

bool 
PreMSTangentComputer::init( KnSpace & ks, const C4CIteratorOnSurface* iter ) 
{
  C4CIteratorOnBdry* it_clone = dynamic_cast<C4CIteratorOnBdry* >(iter->clone());

  m_t_cover.init( *it_clone , 0 ); // REQUIRES ImaGene1.5 spec

  C4CIteratorOnBdry* beg_surfel = dynamic_cast<C4CIteratorOnBdry*> (m_t_cover.beginSurfel());

  m_t_cover_geometry.init( &ks, (uint) 0, (uint) 1, m_t_cover, *beg_surfel);  

  bool init = m_t_cover.OK() && m_t_cover_geometry.OK(); 
  
  delete beg_surfel;
  delete it_clone;
  
  return init;
}

Vector2D
PreMSTangentComputer::computeTangent( KnSpace & ks,
				   const C4CIteratorOnSurface* iter,
				   const Embedder & embedder, 
				   float & window,
				   float & ds,
				   Vector2D & ext_front,
				   Vector2D & ext_back )
{
  if (m_is_first_call_done)
    {
      m_t_cover.nextSMS ( m_pencil_on_cur_idx ); // modify m_pencil_on_cur_idx ...
      m_cur_surf_idx++;
    }
  else {
    if (init(ks, iter))
      ;
    else 
      {
	cerr << " Error during PreMSTangentComputer intialisation ... exiting" << endl;
	exit (1);
      }
    m_is_first_call_done = true;
    m_pencil_on_cur_idx = m_t_cover.beginSMS ( m_cur_surf_idx ); // seems to be in log !! :S
  }

  float theta = ImaGene::C4CTangentialCoverGeometry::angleByLambdaMS( m_t_cover,
								      m_t_cover_geometry,
								      *m_lambda,
								      m_pencil_on_cur_idx);

  // Compute unitary tangent vector from angle to X-axis.
  Vector2D tgt;
  tgt = Vector2D( cos( theta ), sin( theta ) );
  VectorUtils::normalize(tgt);
  return tgt;
}


// Tangents based on median filtering

class Matas95TangentComputer : public  TangentComputer
{
protected:
  static const uint m_type_std = 0;
  static const uint m_type_left = 1;
  static const uint m_type_right = 2;
  
public:
    
  Vector2D computeTangent ( KnSpace & ks,
			    const C4CIteratorOnSurface* iter,
			    const Embedder & embedder, 
			    float & window,
			    float & ds,
			    Vector2D & ext_front,
			    Vector2D & ext_back)

  {
    Vector2D v2d_central_point = ks.scentroid( iter->current() );
    Vector2D v2d_current_direction;
    int wsize = static_cast<int>(size);
    
    Proxy<C4CIteratorOnBdry> it_back( dynamic_cast<C4CIteratorOnBdry*>( iter->clone()));
    Proxy<C4CIteratorOnBdry> it_fwd( dynamic_cast<C4CIteratorOnBdry*>( iter->clone()));
    
    Vector2D a_v2d_in_win[2*wsize];
    
    int j=0;    
    for (int i = 0; i < wsize ; i++) {
      it_back->previous();
      v2d_current_direction = v2d_central_point;
      v2d_current_direction -= ks.scentroid( it_back->current() );
      VectorUtils::normalize (v2d_current_direction);
      a_v2d_in_win[wsize-1-i] = v2d_current_direction;

      it_fwd->next();
      v2d_current_direction = ks.scentroid( it_fwd->current() );
      v2d_current_direction -= v2d_central_point; 
      VectorUtils::normalize (v2d_current_direction);
      a_v2d_in_win[wsize+i] = v2d_current_direction;
    }
    
    // get orientations 
    // on ramene les angles entre [0..2*M_PI]
    // TODO should be a float Vector2D::getAngle const () function...

    bool has_quadrant[4];
    has_quadrant[0]=has_quadrant[1]=has_quadrant[2]=has_quadrant[3]=false;

    float a_ang_in_win[2*wsize];
    for (int i = 0; i < 2*wsize; i++) {

      if ( a_v2d_in_win[i].x() == 0.0 ) 
	{
	  if ( a_v2d_in_win[i].y() >0 )
	    a_ang_in_win[i] = M_PI/2.0; 
	  else 
	    a_ang_in_win[i] = 3*M_PI/2.0; 
	}
      else if (  ( a_v2d_in_win[i].x() > 0.0 ) && (   a_v2d_in_win[i].y() >= 0.0 ) ) {
	a_ang_in_win[i] = atan(a_v2d_in_win[i].y()/a_v2d_in_win[i].x());
	has_quadrant[0]  = true;
      }
      else if (  ( a_v2d_in_win[i].x() > 0.0 ) && (   a_v2d_in_win[i].y() <= 0.0 ) ) {
	a_ang_in_win[i] = 2*M_PI + atan(a_v2d_in_win[i].y()/a_v2d_in_win[i].x());
	has_quadrant[3] = true;
      }
      else if (  ( a_v2d_in_win[i].x() < 0.0 ) && (   a_v2d_in_win[i].y() >= 0.0 ) ) {
	a_ang_in_win[i] = atan(a_v2d_in_win[i].y()/a_v2d_in_win[i].x()) + M_PI;
	has_quadrant[1] = true;
      }
      else if (  ( a_v2d_in_win[i].x() < 0.0 ) && (   a_v2d_in_win[i].y() <= 0.0 ) ) {
	a_ang_in_win[i] = atan(a_v2d_in_win[i].y()/a_v2d_in_win[i].x()) + M_PI;
	has_quadrant[2] = true;
      }
    }

    if ( has_quadrant[0] && has_quadrant[1] && has_quadrant[2] && has_quadrant[3] ) {
      std::cerr << "Comparaison impossible dans 'matas95tangentcomputer' : les 4 quadrants sont utilisés...\n Exiting ..." << std::endl;
      exit(1);
    }
    
    //     if ( has_quadrant[0] && has_quadrant[3] ) {
    //       cerr << "Avant modifications : " << endl;
    //       for (int i = 0; i < 2*wsize; i++) {
    // 	cerr << a_ang_in_win[i] << " "; 
    //       }
    //       cerr << endl;
    //     }

    // on ramène les angles en ajoutant 2 pi à ceux qui ont moins de pi
    // traite les cas pi/2 -- 2pi
    // traite les cas 2pi -- pi/2 -- pi

    if ( has_quadrant[0] && has_quadrant[3] ) {
      for (int i = 0; i < 2*wsize; i++) {
	//	if ( a_ang_in_win[i] < 3.0*M_PI/2.0 ) {
	if ( get_quandrant_from_angle(a_ang_in_win[i]) <= 1  ) {
	  a_ang_in_win[i] += 2*M_PI;
	}
      }
    }
        
    // sort the orientations (could be faster ....)
    // could just sort from wsize-1..2*wsize

    float temp;
    // Vector2D temp_vec;
    for (int i =0; i < 2*wsize-1; i++) {
	   for (int j =0; j < 2*wsize-i-1; j++) {
	if ( a_ang_in_win[j] > a_ang_in_win[j+1]) {
	  temp = a_ang_in_win[j+1];
	  a_ang_in_win[j+1] = a_ang_in_win[j];
	  a_ang_in_win[j] = temp;
	  
	  //   temp_vec = a_v2d_in_win[j+1];
	  // 	  a_v2d_in_win[j+1] = a_v2d_in_win[j];
	  // 	  a_v2d_in_win[j] = temp_vec;
	}
      }
    }

    //     if ( has_quadrant[0] && has_quadrant[3] ) {
    //       for (int i = 0; i < 2*wsize; i++) {
    // 	cerr << a_ang_in_win[i] << " "; 
    //       }
    //       cerr << endl;
    //     }
    
    float angle_res;
    
    //     if ( get_quandrant_from_angle(a_ang_in_win[wsize-1]) == get_quandrant_from_angle(a_ang_in_win[wsize])  )
    //       angle_res = (a_ang_in_win[wsize-1] + a_ang_in_win[wsize])/2.0;
    //     else 
    
    //     if ( ( ( get_quandrant_from_angle( a_ang_in_win[wsize-1] ) == 3 ) && 
    // 	   ( get_quandrant_from_angle( a_ang_in_win[wsize]   ) == 0 ) ) || 
    // 	 ( ( get_quandrant_from_angle( a_ang_in_win[wsize-1] ) == 0 ) && 
    // 	   ( get_quandrant_from_angle( a_ang_in_win[wsize]   ) == 3 ) ) )
    //       angle_res = (a_ang_in_win[wsize-1] + a_ang_in_win[wsize])/2.0 - M_PI;
    //     else 
    
    if ( type == m_type_std )
      angle_res = (a_ang_in_win[wsize-1] + a_ang_in_win[wsize])/2.0;
    else if ( type == m_type_left )
      angle_res = a_ang_in_win[wsize-1];
    else if ( type == m_type_right )
      angle_res = a_ang_in_win[wsize];

    //     = (a_ang_in_win[wsize-1] + a_ang_in_win[wsize]);
    //     if (angle_res > 4*M_PI) { // added first and 4 th quadrant
    //       angle_res -= 2*M_PI;
    //       if (angle_res > 4*M_PI) { // added 2 first quadrant with comparison to a  4th quadrant angle
    // 	angle_res -= 2*M_PI;
    //       }
    //     }
    //     angle_res /= 2.0;


    // not shure this is needed...
    //     while (angle_res > 2.0*M_PI) {
    //       angle_res -= 2.0*M_PI;
    //     }
    
    //     cerr << " Autour de : " 
    // 	 << "(" << v2d_central_point.x() << "," << v2d_central_point.y() << ")"
    // 	 << " median : " << angle_res 
    // 	 << endl;

    v2d_current_direction = Vector2D( cos(angle_res),
				      sin(angle_res));
    VectorUtils::normalize( v2d_current_direction);
    return v2d_current_direction;
  }
};

// gaussian derivative for tangent computation

class GaussianDerivativeTangentComputer : public TangentComputer
{
  
public:
    
  Vector2D computeTangent ( KnSpace & ks,
			    const C4CIteratorOnSurface* iter,
			    const Embedder & embedder, 
			    float & window,
			    float & ds,
			    Vector2D & ext_front,
			    Vector2D & ext_back)

  {
    Vector2D central_point = ks.scentroid( iter->current() );
    Vector2D current_direction;

    C4CIterator* iter_clone = iter->clone();
    C4CIteratorOnBdry *iter_copy;
    iter_copy = dynamic_cast<C4CIteratorOnBdry*>(iter_clone);

    Vector2D mp_CoordsInWin[2*((int) size)+1];
    
    // surfel index k
    for (int i = 0; i < ((int) size) ; i++) 
      iter_copy->previous();
    // surfel indexed k - size

    for (int i =0; i < 2*((int) size) +1 ; i++, iter_copy->next() ) {
      mp_CoordsInWin[i] = ks.scentroid( iter_copy->current() );
    }
    // surfel index k - size + 2 *size +1 = k + size + 1

    for (int i = 0; i < ((int) size) +1; i++) iter_copy->previous();

    float sigma = ((float) size)/3.0;
    
    Vector2D res(0.0,0.0);
    int k,anti_k;
    for (int i =0; i < 2*((int) size)+1; i++) {
      k = i - ((int) size);
      anti_k = 2*((int) size) -i ;
      res.x() += -k*expf( -0.5*(k/sigma)*(k/sigma))/ (sigma*sigma*sqrt(2*M_PI))*mp_CoordsInWin[anti_k].x();
      res.y() += -k*expf( -0.5*(k/sigma)*(k/sigma))/ (sigma*sigma*sqrt(2*M_PI))*mp_CoordsInWin[anti_k].y();
    }    
    
    VectorUtils::normalize(res);
    delete iter_clone;
    return res;
  }
};


// hybrid gaussian derivative for tangent computation
// uses maximal segments to choose the  window size

class HybridGaussianDerivativeTangentComputer : public TangentComputer
{
  
public:
    
  Vector2D computeTangent ( KnSpace & ks,
			    const C4CIteratorOnSurface* iter,
			    const Embedder & embedder, 
			    float & window,
			    float & ds,
			    Vector2D & ext_front,
			    Vector2D & ext_back)

  {
    Vector2D central_point = ks.scentroid( iter->current() );
    Vector2D current_direction;

    C4CIterator* iter_clone = iter->clone();
    C4CIteratorOnBdry *iter_copy;
    iter_copy = dynamic_cast<C4CIteratorOnBdry*>(iter_clone);

    // get size of computation window with maximal segments (front and back)
    C4CIterator* iter1 = iter_copy->clone();
    C4CIterator* iter2 = iter_copy->clone();
    C4CIterator* iter3 = iter_copy->clone();
    C4CIterator* iter4 = iter_copy->clone();
    C4CIteratorOnBdry* iterFwdFrontSMax = dynamic_cast<C4CIteratorOnBdry*>(iter1);
    C4CIteratorOnBdry* iterBwdFrontSMax = dynamic_cast<C4CIteratorOnBdry*>(iter2);
    C4CIteratorOnBdry* iterFwdBackSMax = dynamic_cast<C4CIteratorOnBdry*>(iter3);
    C4CIteratorOnBdry* iterBwdBackSMax = dynamic_cast<C4CIteratorOnBdry*>(iter4);
    

    C4CGeometry::maximalFrontTangent( *iterFwdFrontSMax, *iterBwdFrontSMax, 0);
    C4CGeometry::maximalBackTangent( *iterFwdBackSMax, *iterBwdBackSMax, 0);
    
    // on prends taille = max iter
    int size_Forward =0;
    int size_Backward =0;
    int hybrid_size = 0;
    while (! iter_copy->equals(*iterFwdFrontSMax)) {
      iterFwdFrontSMax->previous(); 
      size_Forward++;
    }
    
    while (! iter_copy->equals(*iterBwdBackSMax)) {
      iterBwdBackSMax->next(); 
      size_Backward++;
    }
    
    if (size_Backward > size_Forward) {
      hybrid_size = size_Backward;
    }
    else {
      hybrid_size = size_Forward;
    }
    
    //usual stuff now

    Vector2D mp_CoordsInWin[2*((int) hybrid_size)+1];
    
    // surfel index k
    for (int i = 0; i < ((int) hybrid_size) ; i++) 
      iter_copy->previous();
    // surfel indexed k - hybrid_size

    for (int i =0; i < 2*((int) hybrid_size) +1 ; i++, iter_copy->next() ) {
      mp_CoordsInWin[i] = ks.scentroid( iter_copy->current() );
    }
    // surfel index k - hybrid_size + 2 *hybrid_size +1 = k + hybrid_size + 1

    for (int i = 0; i < ((int) hybrid_size) +1; i++) iter_copy->previous();

    float sigma = ((float) hybrid_size)/3.0;
    
    Vector2D res(0.0,0.0);
    int k,anti_k;
    for (int i =0; i < 2*((int) hybrid_size)+1; i++) {
      k = i - ((int) hybrid_size);
      anti_k = 2*((int) hybrid_size) -i ;
      res.x() += -k*expf( -0.5*(k/sigma)*(k/sigma))/ (sigma*sigma*sqrt(2*M_PI))*mp_CoordsInWin[anti_k].x();
      res.y() += -k*expf( -0.5*(k/sigma)*(k/sigma))/ (sigma*sigma*sqrt(2*M_PI))*mp_CoordsInWin[anti_k].y();
    }    
    
    VectorUtils::normalize(res);
    delete iter_clone;
    delete iter1;
    delete iter2;
    delete iter3;
    delete iter4;
    return res;
  }
};

// 2nd hybrid gaussian derivative
// winsize = (L_MS) *(1/h)^(1/6)

class Hybrid2GaussianDerivativeTangentComputer : public TangentComputer
{
  
public:
    
  Vector2D computeTangent ( KnSpace & ks,
			    const C4CIteratorOnSurface* iter,
			    const Embedder & embedder, 
			    float & window,
			    float & ds,
			    Vector2D & ext_front,
			    Vector2D & ext_back)

  {
    Vector2D central_point = ks.scentroid( iter->current() );
    Vector2D current_direction;

    C4CIterator* iter_clone = iter->clone();
    C4CIteratorOnBdry *iter_copy;
    iter_copy = dynamic_cast<C4CIteratorOnBdry*>(iter_clone);

    // get size of computation window with maximal segments (front and back)
    C4CIterator* iter1 = iter_copy->clone();
    C4CIterator* iter2 = iter_copy->clone();
    C4CIterator* iter3 = iter_copy->clone();
    C4CIterator* iter4 = iter_copy->clone();
    C4CIteratorOnBdry* iterFwdFrontSMax = dynamic_cast<C4CIteratorOnBdry*>(iter1);
    C4CIteratorOnBdry* iterBwdFrontSMax = dynamic_cast<C4CIteratorOnBdry*>(iter2);
    C4CIteratorOnBdry* iterFwdBackSMax = dynamic_cast<C4CIteratorOnBdry*>(iter3);
    C4CIteratorOnBdry* iterBwdBackSMax = dynamic_cast<C4CIteratorOnBdry*>(iter4);
    

    C4CGeometry::maximalFrontTangent( *iterFwdFrontSMax, *iterBwdFrontSMax, 0);
    C4CGeometry::maximalBackTangent( *iterFwdBackSMax, *iterBwdBackSMax, 0);
    
    // on prends taille = max iter
    int size_Forward =0;
    int size_Backward =0;
    int hybrid_size = 0;
    while (! iter_copy->equals(*iterFwdFrontSMax)) {
      iterFwdFrontSMax->previous(); 
      size_Forward++;
    }
    
    while (! iter_copy->equals(*iterBwdBackSMax)) {
      iterBwdBackSMax->next(); 
      size_Backward++;
    }
    
    if (size_Backward > size_Forward) {
      hybrid_size = size_Backward;
    }
    else {
      hybrid_size = size_Forward;
    }
    
    //change power
    //     cerr <<  "Old_value " << hybrid_size << " ";

    hybrid_size *= (int) powf ( 1.0/grid_step, 1.0/6.0);

    //    cerr <<  "New_value " << hybrid_size << endl;

    //usual stuff now

    Vector2D mp_CoordsInWin[2*((int) hybrid_size)+1];
    
    // surfel index k
    for (int i = 0; i < ((int) hybrid_size) ; i++) 
      iter_copy->previous();
    // surfel indexed k - hybrid_size

    for (int i =0; i < 2*((int) hybrid_size) +1 ; i++, iter_copy->next() ) {
      mp_CoordsInWin[i] = ks.scentroid( iter_copy->current() );
    }
    // surfel index k - hybrid_size + 2 *hybrid_size +1 = k + hybrid_size + 1

    for (int i = 0; i < ((int) hybrid_size) +1; i++) iter_copy->previous();

    float sigma = ((float) hybrid_size)/3.0;
    
    Vector2D res(0.0,0.0);
    int k,anti_k;
    for (int i =0; i < 2*((int) hybrid_size)+1; i++) {
      k = i - ((int) hybrid_size);
      anti_k = 2*((int) hybrid_size) -i ;
      res.x() += -k*expf( -0.5*(k/sigma)*(k/sigma))/ (sigma*sigma*sqrt(2*M_PI))*mp_CoordsInWin[anti_k].x();
      res.y() += -k*expf( -0.5*(k/sigma)*(k/sigma))/ (sigma*sigma*sqrt(2*M_PI))*mp_CoordsInWin[anti_k].y();
    }    
    
    VectorUtils::normalize(res);
    delete iter_clone;
    delete iter1;
    delete iter2;
    delete iter3;
    delete iter4;
    return res;
  }
};

// 3rd hybrid gaussian derivative using a global fixed length
// as average (L_MS)
//



class Hybrid3GaussianDerivativeTangentComputer : public TangentComputer
{
public:

  bool init (KnSpace & ks, const C4CIteratorOnSurface* iter);

  void reset ();

  Vector2D computeTangent( KnSpace & ks,
			   const C4CIteratorOnSurface* iter,
			   const Embedder & embedder, 
			   float & window,
			   float & ds,
			   Vector2D & ext_front,
			   Vector2D & ext_back );

  /*
   * Not very nice... but only for prototyping.
   */
  bool m_is_first_call_done;
  uint m_hybrid_size;
};

void 
Hybrid3GaussianDerivativeTangentComputer::reset ()
{
  m_is_first_call_done = false;
  m_hybrid_size = 0;
}

bool 
Hybrid3GaussianDerivativeTangentComputer::init( KnSpace & ks, const C4CIteratorOnSurface* iter ) 
{
  C4CIteratorOnBdry* it_clone = dynamic_cast<C4CIteratorOnBdry* >(iter->clone());
  C4CTangentialCover m_t_cover; 
  m_t_cover.init( *it_clone , 0 ); // REQUIRES ImaGene1.5 spec
  
  uint nb_ms = (int) m_t_cover.nbMaximalSegments();
  uint sum_ms_length = 0;
  for ( uint i =0; i < nb_ms; i++){
    sum_ms_length += m_t_cover.getMaximalSegment(i).dss.size();
  }
  m_hybrid_size = sum_ms_length/nb_ms;
  
  bool init = m_t_cover.OK();
  
  delete it_clone;
  
  return init;
}

Vector2D
Hybrid3GaussianDerivativeTangentComputer::computeTangent( KnSpace & ks,
							  const C4CIteratorOnSurface* iter,
							  const Embedder & embedder, 
							  float & window,
							  float & ds,
							  Vector2D & ext_front,
							  Vector2D & ext_back )
{
  if (m_is_first_call_done)
    {

    }
  else {
    if (init(ks, iter))
      ;
    else 
      {
	cerr << " Error during PreMSTangentComputer intialisation ... exiting" << endl;
	exit (1);
      }
    m_is_first_call_done = true;
  }

  Vector2D central_point = ks.scentroid( iter->current() );
  Vector2D current_direction;
  
  C4CIterator* iter_clone = iter->clone();
  C4CIteratorOnBdry *iter_copy;
  iter_copy = dynamic_cast<C4CIteratorOnBdry*>(iter_clone);
  
  // on prends taille = max iter
  
  //usual stuff now
  
  Vector2D mp_CoordsInWin[2*((int) m_hybrid_size)+1];
  
  // surfel index k
  for (int i = 0; i < ((int) m_hybrid_size) ; i++) 
    iter_copy->previous();
  // surfel indexed k - m_hybrid_size
  
  for (int i =0; i < 2*((int) m_hybrid_size) +1 ; i++, iter_copy->next() ) {
      mp_CoordsInWin[i] = ks.scentroid( iter_copy->current() );
    }
    // surfel index k - m_hybrid_size + 2 *m_hybrid_size +1 = k + m_hybrid_size + 1

    for (int i = 0; i < ((int) m_hybrid_size) +1; i++) iter_copy->previous();

    float sigma = ((float) m_hybrid_size)/3.0;
    
    Vector2D res(0.0,0.0);
    int k,anti_k;
    for (int i =0; i < 2*((int) m_hybrid_size)+1; i++) {
      k = i - ((int) m_hybrid_size);
      anti_k = 2*((int) m_hybrid_size) -i ;
      res.x() += -k*expf( -0.5*(k/sigma)*(k/sigma))/ (sigma*sigma*sqrt(2*M_PI))*mp_CoordsInWin[anti_k].x();
      res.y() += -k*expf( -0.5*(k/sigma)*(k/sigma))/ (sigma*sigma*sqrt(2*M_PI))*mp_CoordsInWin[anti_k].y();
    }    
    
    VectorUtils::normalize(res);
    delete iter_clone;
    return res;
}






// Implicit parabola fitting LEWINNER 2004 :  Implicit parabola fitting

// this method considers y=f(x) functions 
// orientation is given by the displacement of the border points of the computation window

class ImplicitParabolaFittingTangentComputer : public TangentComputer
{
public:

  Vector2D computeTangent( KnSpace & ks,
			  const C4CIteratorOnSurface* iter,
			  const Embedder & embedder, 
			  float & window,
			  float & ds,
			  Vector2D & ext_front,
			  Vector2D & ext_back)
  {
    C4CIterator* iter_clone = iter->clone();
    C4CIteratorOnBdry *iter_copy;
    iter_copy = dynamic_cast<C4CIteratorOnBdry*>(iter_clone);

    float a, b, c, g, h;
    a=b=c=g=h=0.0;
    float dx, dy;
    dx=dy=0.0;
    Vector2D centroids_coordinate[2*((int) size)+1];
    Vector2D central_point = ks.scentroid( iter_copy->current() );
    Vector2D current_direction;

    // Get to the begining of the computation window
    for (int i = 0; i < ((int) size); i++) 
      iter_copy->previous();

    // get the centroids of the surfels
    for (int i = 0; i < 2*((int) size) +1; i++ ,
	   iter_copy->next()) {
      centroids_coordinate[i] = ks.scentroid( iter_copy->current() );
      //   std::cerr << " (" 
      // 		<< " x = " << centroids_coordinate[i].x() 
      // 		<< " y  = " << centroids_coordinate[i].y()
      // 		<< " )" ;
    }

    // Center the samples wrt to the original point
    Vector2D center = centroids_coordinate[size]; 
    for (int i = 0; i < 2*((int) size) +1; i++) {
      centroids_coordinate[i] -= center;
    }
    
    // get the  variations & orientation
    float x_min, x_max, y_min, y_max;
    x_min = x_max = centroids_coordinate[size].x();
    y_min = y_max = centroids_coordinate[size].y();
    for (int i = 0; i < 2*((int) size) +1; i++) {
      if (centroids_coordinate[i].x() < x_min)
	x_min = centroids_coordinate[i].x();
      else if (centroids_coordinate[i].x() > x_max)
	x_max = centroids_coordinate[i].x();
      if (centroids_coordinate[i].y() < y_min)
	y_min = centroids_coordinate[i].y();
      else if (centroids_coordinate[i].y() > y_max)
	y_max = centroids_coordinate[i].y();
    }
    dx = x_max - x_min; 
    dy = y_max - y_min;

    //  std::cerr << " (" 
    // 	      << " x = " << centroids_coordinate[size].x() 
    // 	      << ", y  = " << centroids_coordinate[size].y()
    // 	      << " ;  dx = " << dx
    // 	      << " , dy = " << dy;
    
    Vector2D xtrem = centroids_coordinate[2*size];
    xtrem  -= centroids_coordinate[0];
    //     std::cerr << " ; " 
    // 	      << " x(2s)-x(0) = " << xtrem.x() 
    // 	      << ",  y(2s)-y(0)  = " << xtrem.y(); 
    
    
    Vector2D v;

    if ( dy  < dx ) {
      for (int i = 0; i < 2*((int) size) +1; i++) {
	a += centroids_coordinate[i].x()*centroids_coordinate[i].x();
	b += 0.5 * centroids_coordinate[i].x()*centroids_coordinate[i].x()*centroids_coordinate[i].x();
	c += 0.25 * centroids_coordinate[i].x()*centroids_coordinate[i].x()*centroids_coordinate[i].x()*centroids_coordinate[i].x();
	g += centroids_coordinate[i].x()*centroids_coordinate[i].y();
	h += 0.5 * centroids_coordinate[i].x()*centroids_coordinate[i].x()*centroids_coordinate[i].y();
      }
      
      // f'' = (ah - bg) / ( ac - bb)
      if (xtrem.x() < 0.0)
	v  = Vector2D(-1.0,-1.0*(c*g-b*h)/(a*c - b*b));
      else 
	v  = Vector2D(1.0,(c*g-b*h)/(a*c - b*b));
    }
    else {
      for (int i = 0; i < 2*((int) size) +1; i++) {
	a += centroids_coordinate[i].y()*centroids_coordinate[i].y();
	b += 0.5 * centroids_coordinate[i].y()*centroids_coordinate[i].y()*centroids_coordinate[i].y();
	c += 0.25 * centroids_coordinate[i].y()*centroids_coordinate[i].y()*centroids_coordinate[i].y()*centroids_coordinate[i].y();
	g += centroids_coordinate[i].y()*centroids_coordinate[i].x();
	h += 0.5 * centroids_coordinate[i].y()*centroids_coordinate[i].y()*centroids_coordinate[i].x();
      }
      
      // f'' = (ah - bg) / ( ac - bb)
      if (xtrem.y() < 0.0)
	v = Vector2D(-1.0*(c*g-b*h)/(a*c - b*b),-1.0 );
      else 
	v = Vector2D((c*g-b*h)/(a*c - b*b),1.0 );
    }
    //    std::cerr << " )" << endl;
    VectorUtils::normalize(v);
    
    delete iter_clone;
    return v;
  }

};


// INCREMENTAL Implicit parabola fitting 
// this method considers y=f(x) functions 
// orientation is given by the displacement of the border points of the computation window

// WARNING
// NOT WORKING PERFECTLY YET !!

class PreImplicitParabolaFittingTangentComputer : public TangentComputer
{
public:

  bool init (KnSpace & ks,const C4CIteratorOnSurface* iter);
  void reset ();
  

  Vector2D computeTangent( KnSpace & ks,
			   const C4CIteratorOnSurface* iter,
			   const Embedder & embedder, 
			   float & window,
			   float & ds,
			   Vector2D & ext_front,
			   Vector2D & ext_back);
  
  
  uint m_cur_surf_idx;
  bool m_is_first_call_done;
  C4CIteratorOnBdry* m_it_beg_window;
  C4CIteratorOnBdry* m_it_end_window;

  float m_a1,m_b1,m_c1,m_g1,m_h1; // for m_dy <  m_dx
  float m_a2,m_b2,m_c2,m_g2,m_h2; // for m_dy >= m_dx


  float m_dx,m_dy;
  float m_x_min, m_x_max, m_y_min, m_y_max;
  Vector2D m_central_point;
  Vector2D m_beg_point;
  Vector2D m_end_point;
  
  Vector2D m_xtrem;

};

bool 
PreImplicitParabolaFittingTangentComputer::init (KnSpace & ks, const C4CIteratorOnSurface* iter)
{
  m_it_beg_window = dynamic_cast<C4CIteratorOnBdry*>(iter->clone());
  m_it_end_window = dynamic_cast<C4CIteratorOnBdry*>(iter->clone());
  m_a1=m_b1=m_c1=m_g1=m_h1=m_a2=m_b2=m_c2=m_g2=m_h2=0.0;
  m_dx=m_dy=0.0;
  m_central_point = ks.scentroid( iter->current() );
  
  C4CIteratorOnBdry *iter_copy = dynamic_cast<C4CIteratorOnBdry*>(iter->clone());

  Vector2D centroids_coordinate[2*((int) size)+1];
  Vector2D central_point = ks.scentroid( iter_copy->current() );
  Vector2D current_direction;
  
  for (int i = 0; i < ((int) size); i++) {
    m_it_beg_window->previous();
    m_it_end_window->next();   
    iter_copy->previous();
  }
  
  // get the centroids of the surfels
  for (int i = 0; i < 2*((int) size) +1; i++ ,
	 iter_copy->next()) {
    centroids_coordinate[i] = ks.scentroid( iter_copy->current() );
  }

//   // Center the samples wrt to the original point
//   Vector2D center = centroids_coordinate[size]; 
//   for (int i = 0; i < 2*((int) size) +1; i++) {
//     centroids_coordinate[i] -= center;
//   }
    
  // get the  variations & orientation

  m_x_min = m_x_max = centroids_coordinate[size].x();
  m_y_min = m_y_max = centroids_coordinate[size].y();
  for (int i = 0; i < 2*((int) size) +1; i++) {
    if (centroids_coordinate[i].x() < m_x_min)
      m_x_min = centroids_coordinate[i].x();
    else if (centroids_coordinate[i].x() > m_x_max)
      m_x_max = centroids_coordinate[i].x();
    if (centroids_coordinate[i].y() < m_y_min)
      m_y_min = centroids_coordinate[i].y();
    else if (centroids_coordinate[i].y() > m_y_max)
      m_y_max = centroids_coordinate[i].y();
  }
  m_dx = m_x_max - m_x_min; 
  m_dy = m_y_max - m_y_min;

  m_xtrem = centroids_coordinate[2*size];
  m_xtrem  -= centroids_coordinate[0];

  
  for (int i = 0; i < 2*((int) size) +1; i++) {
    // for m_dy < m_dx
    m_a1 += centroids_coordinate[i].x()*centroids_coordinate[i].x();
    m_b1 += 0.5 * centroids_coordinate[i].x()*centroids_coordinate[i].x()*centroids_coordinate[i].x();
    m_c1 += 0.25 * centroids_coordinate[i].x()*centroids_coordinate[i].x()*centroids_coordinate[i].x()*centroids_coordinate[i].x();
    m_g1 += centroids_coordinate[i].x()*centroids_coordinate[i].y();
    m_h1 += 0.5 * centroids_coordinate[i].x()*centroids_coordinate[i].x()*centroids_coordinate[i].y();
    // for m_dy >= m_dx
    m_a2 += centroids_coordinate[i].y()*centroids_coordinate[i].y();
    m_b2 += 0.5 * centroids_coordinate[i].y()*centroids_coordinate[i].y()*centroids_coordinate[i].y();
    m_c2 += 0.25 * centroids_coordinate[i].y()*centroids_coordinate[i].y()*centroids_coordinate[i].y()*centroids_coordinate[i].y();
    m_g2 += centroids_coordinate[i].y()*centroids_coordinate[i].x();
    m_h2 += 0.5 * centroids_coordinate[i].y()*centroids_coordinate[i].y()*centroids_coordinate[i].x();
  }

  delete iter_copy;

  

  return true;
}

void 
PreImplicitParabolaFittingTangentComputer::reset ()
{
  m_is_first_call_done = false;
  if ( m_it_end_window != 0 )
    delete m_it_end_window; 
  if ( m_it_beg_window != 0 )
    delete m_it_beg_window; 
}


Vector2D 
PreImplicitParabolaFittingTangentComputer::computeTangent( KnSpace & ks,
							   const C4CIteratorOnSurface* iter,
							   const Embedder & embedder, 
							   float & window,
							   float & ds,
							   Vector2D & ext_front,
							   Vector2D & ext_back)
{
  Vector2D v;
  
  if (m_is_first_call_done)
    {
      // modify values: increment both iterator
      // change values a,b,c,g,h
      // remove last value
      // add new one
      m_central_point = ks.scentroid( iter->current() );
      
      // remove  old values
      
      Vector2D old_beg_win = ks.scentroid( m_it_beg_window->current());
      float old_beg_win_x = old_beg_win.x();
      float old_beg_win_y = old_beg_win.y();

      // for m_dy < m_dx
      m_a1 -= old_beg_win_x*old_beg_win_x; 
      m_b1 -= 0.5 * old_beg_win_x*old_beg_win_x*old_beg_win_x; 
      m_c1 -= 0.25 * old_beg_win_x*old_beg_win_x*old_beg_win_x*old_beg_win_x; 
      m_g1 -= old_beg_win_x * old_beg_win_y; 
      m_h1 -= 0.5 * old_beg_win_x*old_beg_win_x*old_beg_win_y; 
      // for m_dy >= m_dx
      m_a2 -= old_beg_win_y*old_beg_win_y; 
      m_b2 -= 0.5 * old_beg_win_y*old_beg_win_x*old_beg_win_y; 
      m_c2 -= 0.25 * old_beg_win_y*old_beg_win_y*old_beg_win_y*old_beg_win_y; 
      m_g2 -= old_beg_win_y * old_beg_win_x; 
      m_h2 -= 0.5 * old_beg_win_y*old_beg_win_y*old_beg_win_x; 
      
      // update beg iterator
      m_it_beg_window->next();

      // update end iterator
      m_it_end_window->next();
      
      // add new values
      Vector2D new_end_win = ks.scentroid( m_it_end_window->current());
      float new_end_win_x = new_end_win.x();
      float new_end_win_y = new_end_win.y();

      // for m_dy < m_dx
      m_a1 += new_end_win_x*new_end_win_x; 
      m_b1 += 0.5 * new_end_win_x*new_end_win_x*new_end_win_x; 
      m_c1 += 0.25 * new_end_win_x*new_end_win_x*new_end_win_x*new_end_win_x; 
      m_g1 += new_end_win_x * new_end_win_y; 
      m_h1 += 0.5 * new_end_win_x*new_end_win_x*new_end_win_y; 
      // for m_dy >= m_dx
      m_a2 += new_end_win_y*new_end_win_y; 
      m_b2 += 0.5 * new_end_win_y*new_end_win_x*new_end_win_y; 
      m_c2 += 0.25 * new_end_win_y*new_end_win_y*new_end_win_y*new_end_win_y; 
      m_g2 += new_end_win_y * new_end_win_x; 
      m_h2 += 0.5 * new_end_win_y*new_end_win_y*new_end_win_x; 
      
      


    }
  else {
    if (init(ks,iter))
      m_is_first_call_done = true;
    else 
      {
	cerr << " Error during PreImplicitParabolaFittingTangentComputer ... exiting" << endl;
	exit (1);
      }
  }

  //compute value
  
//   if ( m_dy  < m_dx ) {    
//     // f'' = (ah - bg) / ( ac - bb)
//     if (m_xtrem.x() < 0.0)
//       v  = Vector2D(-1.0,-1.0*(m_c1*m_g1-m_b1*m_h1)/(m_a1*m_c1 - m_b1*m_b1) -
// 		    1.0*(m_a1*m_h1-m_b1*m_g1)/(m_a1*m_c1 - m_b1*m_b1)*m_central_point.x()  );
//     else 
      v  = Vector2D(1.0,(m_c1*m_g1-m_b1*m_h1)/(m_a1*m_c1 - m_b1*m_b1) + 
		    1.0*(m_a1*m_h1-m_b1*m_g1)/(m_a1*m_c1 - m_b1*m_b1)*m_central_point.x()  );
//   }
//   else {
//     // f'' = (ah - bg) / ( ac - bb)
//     if (m_xtrem.y() < 0.0)
//       v = Vector2D(-1.0*(m_c2*m_g2-m_b2*m_h2)/(m_a2*m_c2 - m_b2*m_b2)   - 
// 		   1.0*(m_a1*m_h1-m_b1*m_g1)/(m_a1*m_c1 - m_b1*m_b1)*m_central_point.y()
// 		   ,-1.0 );
//     else 
//       v = Vector2D((m_c2*m_g2-m_b2*m_h2)/(m_a2*m_c2 - m_b2*m_b2) + 
// 		   1.0*(m_a1*m_h1-m_b1*m_g1)/(m_a1*m_c1 - m_b1*m_b1)*m_central_point.y() 
// 		   ,1.0 );
//   }
   
    
//     // get the  variations & orientation
//     float x_min, x_max, y_min, y_max;
//     x_min = x_max = centroids_coordinate[size].x();
//     y_min = y_max = centroids_coordinate[size].y();
//     for (int i = 0; i < 2*((int) size) +1; i++) {
//       if (centroids_coordinate[i].x() < x_min)
// 	x_min = centroids_coordinate[i].x();
//       else if (centroids_coordinate[i].x() > x_max)
// 	x_max = centroids_coordinate[i].x();
//       if (centroids_coordinate[i].y() < y_min)
// 	y_min = centroids_coordinate[i].y();
//       else if (centroids_coordinate[i].y() > y_max)
// 	y_max = centroids_coordinate[i].y();
//     }
//     dx = x_max - x_min; 
//     dy = y_max - y_min;

//     //  std::cerr << " (" 
//     // 	      << " x = " << centroids_coordinate[size].x() 
//     // 	      << ", y  = " << centroids_coordinate[size].y()
//     // 	      << " ;  dx = " << dx
//     // 	      << " , dy = " << dy;
    
//     Vector2D xtrem = centroids_coordinate[2*size];
//     xtrem  -= centroids_coordinate[0];
//     //     std::cerr << " ; " 
//     // 	      << " x(2s)-x(0) = " << xtrem.x() 
//     // 	      << ",  y(2s)-y(0)  = " << xtrem.y(); 
    
    
  VectorUtils::normalize(v);
  
  return v;
}


// Explicit parabola fitting :

// we fit P_3 (X) = a_0 + a_1 X + a_2 X^2   over 2M+1 points

// Solution is computed using the normal equations:
// AX=B        =>           tAAX=tAB                 =>  X=(tAA)^(-1) AB   
//
// B = ( y_-M )   X = ( a_0 )  A = ( 1  x_M   (x_M)^2  )
//     (  ... )       ( a_1 )      ( 1  ...     ...    )
//     ( y_M  )       ( a_2 )      ( 1  x_-M  (x_-M)^2 ) 
//
// B is a 2M+1 vector and A is a (2M+1,3) matrix, X is a 3 vector
//
//  tAA = ( a   b   c)   det(tAA) = acf+2bec - c^3 - ae^2 - fb^2
//        ( b   c   e)
//        ( c   e   f)
//
//  Com(tAA) = ( cf-e^2          ce -bf           be -c^2 ) = tCom(tAA) = ( a2  b2  c2 )
//             ( ce-bf           af-c^2           bc -ae )                ( b2  d2  e2 )
//             ( be-c^2          bc -ae           ac-b^2 )                ( c2  e2  f2 )
//
//  tAB = ( g ) = ( sum y_i         )
//        ( h )   ( sum x_i y_i     )
//        ( k )   ( sum (x_i)^2 y_i )
//
//  det(tAA) X = ( g a2 + h b2 + k c2 )
//               ( g b2 + h d2 + k e2 )
//               ( g c2 + h e2 + k f2 )
//
// this method considers y=f(x) functions 
// orientation is given by the displacement of the border points of the computation window

class ExplicitParabolaFittingTangentComputer : public TangentComputer
{
public:

  Vector2D computeTangent( KnSpace & ks,
			  const C4CIteratorOnSurface* iter,
			  const Embedder & embedder, 
			  float & window,
			  float & ds,
			  Vector2D & ext_front,
			  Vector2D & ext_back)
  {
    C4CIterator* iter_clone = iter->clone();
    C4CIteratorOnBdry *iter_copy;
    iter_copy = dynamic_cast<C4CIteratorOnBdry*>(iter_clone);
    
    // coefficients
    float a, b, c, e, f; // tAA
    float a2, b2, c2, d2, e2, f2; // com(tAA)
    float g, h, k; // tAB
    a=b=c=e=f=a2=b2=c2=d2=e2=f2=g=h=k=0.0;
    float dx, dy;
    dx=dy=0.0;
    Vector2D centroids_coordinate[2*((int) size)+1];
    Vector2D central_point = ks.scentroid( iter_copy->current() );
    Vector2D current_direction;

    // Get to the begining of the computation window
    for (int i = 0; i < ((int) size); i++) 
      iter_copy->previous();

    // get the centroids of the surfels
    for (int i = 0; i < 2*((int) size) +1; i++ ,
	   iter_copy->next()) {
      centroids_coordinate[i] = ks.scentroid( iter_copy->current() );
      //   std::cerr << " (" 
      // 		<< " x = " << centroids_coordinate[i].x() 
      // 		<< " y  = " << centroids_coordinate[i].y()
      // 		<< " )" ;
    }

    // Center the samples wrt to the original point
    Vector2D center = centroids_coordinate[size]; 
    for (int i = 0; i < 2*((int) size) +1; i++) {
      centroids_coordinate[i] -= center;
    }
    
    // get the  variations & orientation
    float x_min, x_max, y_min, y_max;
    x_min = x_max = centroids_coordinate[size].x();
    y_min = y_max = centroids_coordinate[size].y();
    for (int i = 0; i < 2*((int) size) +1; i++) {
      if (centroids_coordinate[i].x() < x_min)
	x_min = centroids_coordinate[i].x();
      else if (centroids_coordinate[i].x() > x_max)
	x_max = centroids_coordinate[i].x();
      if (centroids_coordinate[i].y() < y_min)
	y_min = centroids_coordinate[i].y();
      else if (centroids_coordinate[i].y() > y_max)
	y_max = centroids_coordinate[i].y();
    }
    dx = x_max - x_min; 
    dy = y_max - y_min;

    //  std::cerr << " (" 
    // 	      << " x = " << centroids_coordinate[size].x() 
    // 	      << ", y  = " << centroids_coordinate[size].y()
    // 	      << " ;  dx = " << dx
    // 	      << " , dy = " << dy;
    
    Vector2D xtrem = centroids_coordinate[2*size];
    xtrem  -= centroids_coordinate[0];
    //     std::cerr << " ; " 
    // 	      << " x(2s)-x(0) = " << xtrem.x() 
    // 	      << ",  y(2s)-y(0)  = " << xtrem.y(); 
    
    
    Vector2D v;

    if ( dy  < dx ) {
      a = 2*((int) size) +1;
      for (int i = 0; i < 2*((int) size) +1; i++) {
	b += centroids_coordinate[i].x();
	c += centroids_coordinate[i].x()*centroids_coordinate[i].x();
	e += centroids_coordinate[i].x()*centroids_coordinate[i].x()*centroids_coordinate[i].x();
	f += centroids_coordinate[i].x()*centroids_coordinate[i].x()*centroids_coordinate[i].x()*centroids_coordinate[i].x();
	g += centroids_coordinate[i].y();
	h += centroids_coordinate[i].x()*centroids_coordinate[i].y();
	k += centroids_coordinate[i].x()*centroids_coordinate[i].x()*centroids_coordinate[i].y();
      }
      float det_tAA = a*c*f + 2*b*e*c - c*c*c -a*e*e -f*b*b;
      b2 = c*e-b*f;
      d2 = a*f -c*c;
      e2 = b*c -a*e;
      
      if (xtrem.x() < 0.0)
	v  = Vector2D(-1.0,-1.0*(b2*g + d2*h + e2*k)/det_tAA);
      else 
	v  = Vector2D(1.0,(b2*g + d2*h + e2*k)/det_tAA);
    }
    else {
      a = 2*((int) size) +1;
      for (int i = 0; i < 2*((int) size) +1; i++) {
	b += centroids_coordinate[i].y();
	c += centroids_coordinate[i].y()*centroids_coordinate[i].y();
	e += centroids_coordinate[i].y()*centroids_coordinate[i].y()*centroids_coordinate[i].y();
	f += centroids_coordinate[i].y()*centroids_coordinate[i].y()*centroids_coordinate[i].y()*centroids_coordinate[i].y();
	g += centroids_coordinate[i].x();
	h += centroids_coordinate[i].y()*centroids_coordinate[i].x();
	k += centroids_coordinate[i].y()*centroids_coordinate[i].y()*centroids_coordinate[i].x();
      }
      float det_tAA = a*c*f + 2*b*e*c - c*c*c -a*e*e -f*b*b;
      b2 = c*e-b*f;
      d2 = a*f -c*c;
      e2 = b*c -a*e;
      
      if (xtrem.y() < 0.0)
	v  = Vector2D(-1.0*(b2*g + d2*h + e2*k)/det_tAA,-1.0);
      else 
	v  = Vector2D((b2*g + d2*h + e2*k)/det_tAA, 1.0);
    }
    //    std::cerr << " )" << endl;
    VectorUtils::normalize(v);
    
    delete iter_clone;
    return v;
  }
  
};

// Implicit parabola fitting on each coordinates
// LEWINNER 2004 :  Implicit parabola fitting with independent coordinates 
// we now apply the preceeding methods on each component of the curve that is implicit parabola fitting on 

// 1   2   3   4   5                             1   2   3   4   5
//                     on the one hand, and on                       on the other hand
// x1  x2  x3  x4  x5                            y1  y2  y3  y4  y5

// we then combine both to get the derivative of the digital curve at the considered point.

// WARNING THIS IS USING A TRIVIAL PARAMETRISATION OF THE CURVE
// + no need to check for bigger variation and orientation, implicit with the parametrisation

class SeparatedImplicitParabolaFittingTangentComputer : public TangentComputer
{
public:

  Vector2D computeTangent( KnSpace & ks,
			  const C4CIteratorOnSurface* iter,
			  const Embedder & embedder, 
			  float & window,
			  float & ds,
			  Vector2D & ext_front,
			  Vector2D & ext_back)
  {
    C4CIterator* iter_clone = iter->clone();
    C4CIteratorOnBdry *iter_copy;
    iter_copy = dynamic_cast<C4CIteratorOnBdry*>(iter_clone);

    float a, b, c, g, h;
    a=b=c=g=h=0.0;
    float dx, dy;
    dx=dy=0.0;
    Vector2D centroids_coordinate[2*((int) size)+1];
    Vector2D central_point = ks.scentroid( iter_copy->current() );
    Vector2D current_direction;

    // Get to the begining of the computation window
    for (int i = 0; i < ((int) size); i++) 
      iter_copy->previous();

    // get the centroids of the surfels
    for (int i = 0; i < 2*((int) size) +1; i++ ,
	   iter_copy->next()) {
      centroids_coordinate[i] = ks.scentroid( iter_copy->current() );
      //   std::cerr << " (" 
      // 		<< " x = " << centroids_coordinate[i].x() 
      // 		<< " y  = " << centroids_coordinate[i].y()
      // 		<< " )" ;
    }

    // center the samples wrt to the original point
    Vector2D center = centroids_coordinate[size]; 
    for (int i = 0; i < 2*((int) size) +1; i++) {
      centroids_coordinate[i] -= center;
    }
       
    Vector2D v;
    float der_x, der_y;
    for (int i = - ((int) size) ; i < ((int) size) +1; i++) {
      a += (float) i*i;
      b += 0.5 * i*i*i;
      c += 0.25 * i*i*i*i;
      g += i*centroids_coordinate[i+((int) size)].x();
      h += 0.5 * i*i*centroids_coordinate[i+((int) size)].x();
    }

      der_x = (c*g-b*h)/(a*c - b*b);

    a=b=c=g=h=0.0;

    for (int i =- ((int) size) ; i < ((int) size) +1; i++) {
      a += (float) i*i;
      b += 0.5 * i*i*i;
      c += 0.25 * i*i*i*i;
      g += i*centroids_coordinate[i+((int) size)].y();
      h += 0.5 * i*i*centroids_coordinate[i+((int) size)].y();
    }
    
    der_y = (c*g-b*h)/(a*c - b*b);
    
    v  = Vector2D ( der_x,der_y);
    VectorUtils::normalize(v);

    delete iter_clone;
    return v;
  }

};

// SEIPF : explicit parabola fittig on each coordinate

class SeparatedExplicitParabolaFittingTangentComputer : public TangentComputer
{
public:

  Vector2D computeTangent( KnSpace & ks,
			  const C4CIteratorOnSurface* iter,
			  const Embedder & embedder, 
			  float & window,
			  float & ds,
			  Vector2D & ext_front,
			  Vector2D & ext_back)
  {
    C4CIterator* iter_clone = iter->clone();
    C4CIteratorOnBdry *iter_copy;
    iter_copy = dynamic_cast<C4CIteratorOnBdry*>(iter_clone);
    
    // coefficients
    float a, b, c, e, f; // tAA
    float a2, b2, c2, d2, e2, f2; // com(tAA)
    float g, h, k; // tAB
    a=b=c=e=f=a2=b2=c2=d2=e2=f2=g=h=k=0.0;
    float dx, dy;
    dx=dy=0.0;
    Vector2D centroids_coordinate[2*((int) size)+1];
    Vector2D central_point = ks.scentroid( iter_copy->current() );
    Vector2D current_direction;

    // Get to the begining of the computation window
    for (int i = 0; i < ((int) size); i++) 
      iter_copy->previous();

    // get the centroids of the surfels
    for (int i = 0; i < 2*((int) size) +1; i++ ,
	   iter_copy->next()) {
      centroids_coordinate[i] = ks.scentroid( iter_copy->current() );
    }

    // Center the samples wrt to the original point
    Vector2D center = centroids_coordinate[size]; 
    for (int i = 0; i < 2*((int) size) +1; i++) {
      centroids_coordinate[i] -= center;
    }
    
    Vector2D v;
    float der_x, der_y;
    
    der_x = (c*g-b*h)/(a*c - b*b);
    
    a = 2*((int) size) +1;
    for (int i = - ((int) size) ; i < ((int) size) +1; i++) {
      b += (float) i;
      c += (float) i*i;
      e += (float) i*i*i;
      f += (float) i*i*i*i;
      g += centroids_coordinate[i+((int) size)].x();
      h += i*centroids_coordinate[i+((int) size)].x();
      k += i*i*centroids_coordinate[i+((int) size)].x();
    }
    float det_tAA = a*c*f + 2*b*e*c - c*c*c -a*e*e -f*b*b;
    b2 = c*e-b*f;
    d2 = a*f -c*c;
    e2 = b*c -a*e;
    
    der_x = (b2*g + d2*h + e2*k)/det_tAA;
    
    g=h=k=0.0;
    
    for (int i = - ((int) size) ; i < ((int) size) +1; i++) {
      g += centroids_coordinate[i+((int) size)].y();
      h += i*centroids_coordinate[i+((int) size)].y();
      k += i*i*centroids_coordinate[i+((int) size)].y();
    }
    
    der_y = (b2*g + d2*h + e2*k)/det_tAA;
    
    v  = Vector2D(der_x, der_y);
    
    VectorUtils::normalize(v);
    
    delete iter_clone;
    return v;
  }
  
};


// hybrid Independent coordinates method  for tangent computation
// uses maximal segments to choose the  window size
// WARNING THIS IS USING A TRIVIAL PARAMETRISATION OF THE CURVE

class HybridSeparatedImplicitParabolaFittingTangentComputer : public TangentComputer
{
public:

  Vector2D computeTangent( KnSpace & ks,
			  const C4CIteratorOnSurface* iter,
			  const Embedder & embedder, 
			  float & window,
			  float & ds,
			  Vector2D & ext_front,
			  Vector2D & ext_back)
  {
    // C4CIterator* iter_copy = iter.clone();
    C4CIterator* iter_clone = iter->clone();
    C4CIteratorOnBdry *iter_copy;
    iter_copy = dynamic_cast<C4CIteratorOnBdry*>(iter_clone);

    // get size of computation window with maximal segments (front and back)
    C4CIterator* iter1 = iter_copy->clone();
    C4CIterator* iter2 = iter_copy->clone();
    C4CIterator* iter3 = iter_copy->clone();
    C4CIterator* iter4 = iter_copy->clone();
    C4CIteratorOnBdry* iterFwdFrontSMax = dynamic_cast<C4CIteratorOnBdry*>(iter1);
    C4CIteratorOnBdry* iterBwdFrontSMax = dynamic_cast<C4CIteratorOnBdry*>(iter2);
    C4CIteratorOnBdry* iterFwdBackSMax = dynamic_cast<C4CIteratorOnBdry*>(iter3);
    C4CIteratorOnBdry* iterBwdBackSMax = dynamic_cast<C4CIteratorOnBdry*>(iter4);
    
    C4CGeometry::maximalFrontTangent( *iterFwdFrontSMax, *iterBwdFrontSMax, 0);
    C4CGeometry::maximalBackTangent( *iterFwdBackSMax, *iterBwdBackSMax, 0);
    
    // on prends taille = max iter
    int size_Forward =0;
    int size_Backward =0;
    int hybrid_size = 0;
    while (! iter_copy->equals(*iterFwdFrontSMax)) {
      iterFwdFrontSMax->previous(); 
      size_Forward++;
    }
    
    while (! iter_copy->equals(*iterBwdBackSMax)) {
      iterBwdBackSMax->next(); 
      size_Backward++;
    }
    
    if (size_Backward > size_Forward) {
      hybrid_size = size_Backward;
    }
    else {
      hybrid_size = size_Forward;
    }

    // usual stuff now

    float a, b, c, g, h;
    a=b=c=g=h=0.0;
    float dx, dy;
    dx=dy=0.0;
    Vector2D centroids_coordinate[2*((int) hybrid_size)+1];
    Vector2D central_point = ks.scentroid( iter_copy->current() );
    Vector2D current_direction;

    for (int i = 0; i < ((int) hybrid_size); i++) 
      iter_copy->previous();

    // We are now at the begining of the computation window
    // get the centroids of the surfels
    
    for (int i = 0; i < 2*((int) hybrid_size) +1; i++ ,
	   iter_copy->next()) {
      centroids_coordinate[i] = ks.scentroid( iter_copy->current() );
    }

    // we center the samples wrt to the original point
    
    Vector2D center = centroids_coordinate[hybrid_size]; 
     for (int i = 0; i < 2*((int) hybrid_size) +1; i++) {
       centroids_coordinate[i] -= center;
     }
       
     // WARNING THIS IS USING A TRIVIAL PARAMETRISATION OF THE CURVE    

    Vector2D v;
    float der_x, der_y;
    for (int i = - ((int) hybrid_size) ; i < ((int) hybrid_size) +1; i++) {
      a += (float) i*i;
      b += 0.5 * i*i*i;
      c += 0.25 * i*i*i*i;
      g += i*centroids_coordinate[i+((int) hybrid_size)].x();
      h += 0.5 * i*i*centroids_coordinate[i+((int) hybrid_size)].x();
    }
 
    der_x = (c*g-b*h)/(a*c - b*b);
    a=b=c=g=h=0.0;

    // WARNING THIS IS USING A TRIVIAL PARAMETRISATION OF THE CURVE

    for (int i =- ((int) hybrid_size) ; i < ((int) hybrid_size) +1; i++) {
      a += (float) i*i;
      b += 0.5 * i*i*i;
      c += 0.25 * i*i*i*i;
      g += i*centroids_coordinate[i+((int) hybrid_size)].y();
      h += 0.5 * i*i*centroids_coordinate[i+((int) hybrid_size)].y();
    }

    der_y = (c*g-b*h)/(a*c - b*b);

    v  = Vector2D ( der_x,der_y);

    VectorUtils::normalize (v);

    delete iter_clone;
    delete iter1;
    delete iter2;
    delete iter3;
    delete iter4;
    return v;
  }

};

//  Hybrid 3 Separated implicit parabalo fitting
//
//   size is determined the average length of maximal segments
//
//

class Hybrid3SeparatedImplicitParabolaFittingTangentComputer : public TangentComputer
{
public:

  bool init (KnSpace & ks, const C4CIteratorOnSurface* iter);

  void reset ();

  Vector2D computeTangent( KnSpace & ks,
			   const C4CIteratorOnSurface* iter,
			   const Embedder & embedder, 
			   float & window,
			   float & ds,
			   Vector2D & ext_front,
			   Vector2D & ext_back );

  /*
   * Not very nice... but only for prototyping.
   */
  bool m_is_first_call_done;
  uint m_hybrid_size;
};

void 
Hybrid3SeparatedImplicitParabolaFittingTangentComputer::reset ()
{
  m_is_first_call_done = false;
  m_hybrid_size = 0;
}

bool 
Hybrid3SeparatedImplicitParabolaFittingTangentComputer::init( KnSpace & ks, const C4CIteratorOnSurface* iter ) 
{
  C4CIteratorOnBdry* it_clone = dynamic_cast<C4CIteratorOnBdry* >(iter->clone());
  C4CTangentialCover m_t_cover; 
  m_t_cover.init( *it_clone , 0 ); // REQUIRES ImaGene1.5 spec
  
  uint nb_ms = (int) m_t_cover.nbMaximalSegments();
  uint sum_ms_length = 0;
  for ( uint i =0; i < nb_ms; i++){
    sum_ms_length += m_t_cover.getMaximalSegment(i).dss.size();
  }
  m_hybrid_size = sum_ms_length/nb_ms;
  
  bool init = m_t_cover.OK();
  
  delete it_clone;
  
  return init;
}

Vector2D
Hybrid3SeparatedImplicitParabolaFittingTangentComputer::computeTangent( KnSpace & ks,
							  const C4CIteratorOnSurface* iter,
							  const Embedder & embedder, 
							  float & window,
							  float & ds,
							  Vector2D & ext_front,
							  Vector2D & ext_back )
{
  if (m_is_first_call_done)
    {
    }
  else {
    if (init(ks, iter))
      ;
    else 
      {
	cerr << " Error during PreMSTangentComputer intialisation ... exiting" << endl;
	exit (1);
      }
    m_is_first_call_done = true;
  }
  
  C4CIterator* iter_clone = iter->clone();
  C4CIteratorOnBdry *iter_copy;
  iter_copy = dynamic_cast<C4CIteratorOnBdry*>(iter_clone);
  
  // usual stuff now
  
  float a, b, c, g, h;
  a=b=c=g=h=0.0;
  float dx, dy;
  dx=dy=0.0;
  Vector2D centroids_coordinate[2*((int) m_hybrid_size)+1];
  Vector2D central_point = ks.scentroid( iter_copy->current() );
  Vector2D current_direction;
  
  for (int i = 0; i < ((int) m_hybrid_size); i++) 
    iter_copy->previous();
  
  // We are now at the begining of the computation window
  // get the centroids of the surfels
  
  for (int i = 0; i < 2*((int) m_hybrid_size) +1; i++ ,
	 iter_copy->next()) {
    centroids_coordinate[i] = ks.scentroid( iter_copy->current() );
  }
  
  // we center the samples wrt to the original point
  
  Vector2D center = centroids_coordinate[m_hybrid_size]; 
  for (int i = 0; i < 2*((int) m_hybrid_size) +1; i++) {
    centroids_coordinate[i] -= center;
  }
  
  // WARNING THIS IS USING A TRIVIAL PARAMETRISATION OF THE CURVE    
  
  Vector2D v;
  float der_x, der_y;
  for (int i = - ((int) m_hybrid_size) ; i < ((int) m_hybrid_size) +1; i++) {
    a += (float) i*i;
    b += 0.5 * i*i*i;
    c += 0.25 * i*i*i*i;
    g += i*centroids_coordinate[i+((int) m_hybrid_size)].x();
    h += 0.5 * i*i*centroids_coordinate[i+((int) m_hybrid_size)].x();
  }
  
  der_x = (c*g-b*h)/(a*c - b*b);
  a=b=c=g=h=0.0;
  
  // WARNING THIS IS USING A TRIVIAL PARAMETRISATION OF THE CURVE
  
  for (int i =- ((int) m_hybrid_size) ; i < ((int) m_hybrid_size) +1; i++) {
    a += (float) i*i;
    b += 0.5 * i*i*i;
    c += 0.25 * i*i*i*i;
    g += i*centroids_coordinate[i+((int) m_hybrid_size)].y();
    h += 0.5 * i*i*centroids_coordinate[i+((int) m_hybrid_size)].y();
  }
  
  der_y = (c*g-b*h)/(a*c - b*b);
  
  v  = Vector2D ( der_x,der_y);
  
  VectorUtils::normalize (v);
  
  delete iter_clone;
  return v;
}


// Regression linéaire sans ponderation 
// sous la forme y=f(x)

class RegLinTangentComputer : public TangentComputer
{
public:

  Vector2D computeTangent( KnSpace & ks,
			  const C4CIteratorOnSurface* iter,
			  const Embedder & embedder, 
			  float & window,
			  float & ds,
			  Vector2D & ext_front,
			  Vector2D & ext_back)
  {
    C4CIterator* iter_clone = iter->clone();
    C4CIteratorOnBdry *iter_copy;
    iter_copy = dynamic_cast<C4CIteratorOnBdry*>(iter_clone);

    float a, b, c, d;
    a=b=c=d=0.0;
    float dx, dy;
    dx=dy=0.0;
    Vector2D centroids_coordinate[2*size+1];
    Vector2D central_point = ks.scentroid( iter_copy->current() );
    Vector2D current_direction;

    // Back to the begining of the computation window
    for (int i = 0; i < ((int) size); i++) 
      iter_copy->previous();
    // get the centroids of the surfels
    for (int i = 0; i < 2*((int) size) +1; i++ ,
	   iter_copy->next()) {
      centroids_coordinate[i] = ks.scentroid( iter_copy->current() );
    }

    // we center the samples wrt to the original point ?? usefulness ??
    
    Vector2D center = centroids_coordinate[size]; 
     for (int i = 0; i < 2*((int) size) +1; i++) {
       centroids_coordinate[i] -= center;
     }
    
     // get the  variations
     // bad !! should be average std. dev. imho
     
    float x_min, x_max, y_min, y_max;
    x_min = x_max = centroids_coordinate[size].x();
    y_min = y_max = centroids_coordinate[size].y();
    for (int i = 0; i < 2*((int) size) +1; i++) {
      if (centroids_coordinate[i].x() < x_min)
	x_min = centroids_coordinate[i].x();
      else if (centroids_coordinate[i].x() > x_max)
	x_max = centroids_coordinate[i].x();
      if (centroids_coordinate[i].y() < y_min)
	y_min = centroids_coordinate[i].y();
      else if (centroids_coordinate[i].y() > y_max)
	y_max = centroids_coordinate[i].y();
    }
    dx = x_max - x_min; 
    dy = y_max - y_min;

    Vector2D xtrem = centroids_coordinate[2*size];
    xtrem  -= centroids_coordinate[0];
    
    Vector2D v;

    if ( dy  < dx ) {
      for (int i = 0; i < 2*((int) size) +1; i++) {
	a += centroids_coordinate[i].x();
	b += centroids_coordinate[i].x()*centroids_coordinate[i].x();
	c += centroids_coordinate[i].y();
	d += centroids_coordinate[i].x()*centroids_coordinate[i].y();
      }
      
      if (xtrem.x() < 0.0)
	v  = Vector2D(-1.0,-1.0*( (2*size +1)*d  - a*c)/((2*size +1)*b - a*a));
      else 
	v  = Vector2D(1.0,( (2*size +1)*d  - a*c)/((2*size +1)*b - a*a));
    }
    else {
      for (int i = 0; i < 2*((int) size) +1; i++) {
	a += centroids_coordinate[i].y();
	b += centroids_coordinate[i].y()*centroids_coordinate[i].y();
	c += centroids_coordinate[i].x();
	d += centroids_coordinate[i].y()*centroids_coordinate[i].x();
      }
      
      if (xtrem.y() < 0.0)
	v = Vector2D(-1.0*( (2*size +1)*d  - a*c)/((2*size +1)*b - a*a),-1.0 );
      else 
	v = Vector2D( ( (2*size +1)*d  - a*c)/((2*size +1)*b - a*a),1.0 );
      
    }
      VectorUtils::normalize(v);
      delete iter_copy;
      return v;
  }

};



///////////////////////////////////////////////////////////////////////////////
// DATAS
///////////////////////////////////////////////////////////////////////////////

struct PerSurfelData 
{
  Kn_sid surfel;
  float theta;
  float xc, yc;
  Vector2D est_tgt;
  float window;
  float angle_to_x;
  float th_angle_to_x;
  float angle_deviation;
  float curvature;
  float est_ds;
  float est_curv;
  Vector2D ext_front;
  Vector2D ext_back;
  // François's changes
  float th_max_angle_to_x; 
  float th_min_angle_to_x;
  float th_moy_angle_to_x; 
  // end françois's changes
};


/**
 * this function completes the data associated with each surfel. It
 * uses 'sd.est_tgt', 'sd.theta' and 'sd.surfel' to compute the other
 * information. 'sd.window' is left unchanged.
 */
void 
completeSurfelData( KnSpace & ks,
		    PerSurfelData & sd,
		    const StarShaped & starshape, 
		    const Embedder & embedder )
{
  Vector xbel( 2 );
  embedder.sembed( sd.surfel, xbel );
  float theta = starshape.parameter( xbel );   // surfel's centroid angle to origin
  Vector2D t = starshape.tangent( theta );     // theoretical tangent computed at surfel's centroid angle
  float curv = starshape.curvature( theta );

  uint i = *(ks.sbegin_dirs( sd.surfel ) );
  bool i_pos = ks.sdirect( sd.surfel, i );
  Kn_sid pt2 = ks.sincident( sd.surfel, i, i_pos );
  Kn_sid pt1 = ks.sincident( sd.surfel, i, ! i_pos );
  embedder.sembed( pt2, xbel );
  float theta2 = starshape.parameter( xbel ); // incident-vertex's angle to origin (positive orientation)
  embedder.sembed( pt1, xbel );
  float theta1 = starshape.parameter( xbel ); // incident-vertex's angle to origin (negative orientation)
  
  embedder.sembed( sd.surfel, xbel );         // We set xbel back to what it was (ie the surfel of interest)
  
  // françois 's changes
  // theoretical(ly computed) min, max and moy angle to x-axis along surfel's projection

  if ( fabs( theta - theta2) > M_PI )  // shift so that theta and thetha2 share the same quadrant
      {
	if ( theta - theta2 > 0.0 )
	  theta2 += 2.0 * M_PI; 
	else 
	  theta2 -= 2.0 * M_PI;
      }
  

  if ( fabs( theta - theta1) > M_PI )  // shift so that theta and thetha1 share the same quadrant
      {
	if ( theta - theta1 > 0.0 )
	  theta1 += 2.0 * M_PI; 
	else 
	  theta1 -= 2.0 * M_PI;
      }
  
  float tmin_a_to_x,tmax_a_to_x; 
  float tmoy_a_to_x=0.0f; 
  float pas_de_calcul = 10.0f; // magic number : # of subdivisions along the surfel
  float c_t_a_x,s_t_a_x,th_a_t_x;
  float norme_euc_local_tan;
  float somme_n_e_l_t=0.0f;
  Vector2D local_tan;
  float local_angle;
  float first_angle;
  // we compute the tangent's centroid theoretical angle to x-axis
  c_t_a_x = VectorUtils::dotProduct( Vector2D( 1.0, 0.0 ) , t );
  c_t_a_x = c_t_a_x > 1.0 ? 1.0 : ( c_t_a_x < -1.0 ? -1.0 : c_t_a_x );
  s_t_a_x = VectorUtils::det( Vector2D( 1.0, 0.0 ) , t );
  if ( s_t_a_x >= 0 ) th_a_t_x = acos( c_t_a_x );
  else th_a_t_x = 2 * M_PI - acos( c_t_a_x );
  
  // reference value for further quadrant comparisons
  first_angle = th_a_t_x;       
  
  tmin_a_to_x = tmax_a_to_x = th_a_t_x;   
  for (int i =0; i < pas_de_calcul; ++i) { // != L788
    local_angle = theta1 + ((float) i)*(theta2 - theta1)/pas_de_calcul;
    local_tan = starshape.tangent( local_angle);
    norme_euc_local_tan = sqrt(local_tan.x()*local_tan.x() + local_tan.y()*local_tan.y());
    somme_n_e_l_t += norme_euc_local_tan;
    
    c_t_a_x = VectorUtils::dotProduct( Vector2D( 1.0, 0.0 ) , local_tan );
    c_t_a_x = c_t_a_x > 1.0 ? 1.0 : ( c_t_a_x < -1.0 ? -1.0 : c_t_a_x );
    s_t_a_x = VectorUtils::det( Vector2D( 1.0, 0.0 ) , local_tan );
    if ( s_t_a_x >= 0 ) th_a_t_x = acos( c_t_a_x );
    else th_a_t_x = 2 * M_PI - acos( c_t_a_x );
    
    // we shift the tangent's th.a.to.x-axis so that it shares the quadrant of the reference value
    if ( fabs( first_angle - th_a_t_x) > M_PI )
      {
	if ( first_angle - th_a_t_x > 0.0 )
	  th_a_t_x += 2.0 * M_PI; 
	else 
	  th_a_t_x -= 2.0 * M_PI;
      }
    
    tmoy_a_to_x += norme_euc_local_tan * th_a_t_x;
    if (tmin_a_to_x > th_a_t_x) tmin_a_to_x = th_a_t_x; 
    if (tmax_a_to_x < th_a_t_x) tmax_a_to_x = th_a_t_x;
  }
  tmoy_a_to_x /= somme_n_e_l_t; 
  // end françois's changes
    
//   cerr << "completing surfel=" << sd.surfel 
//        << " tgt=" << sd.est_tgt.ro( 0 ) << " " << sd.est_tgt.ro( 1 )
//        << endl;
  
  // Angle deviation between estimated tangent and estimated tangent.
  float angle_deviation;
  float cos_angle = VectorUtils::dotProduct( t, sd.est_tgt );
  cos_angle = cos_angle > 1.0 ? 1.0 
    : ( cos_angle < -1.0 ? -1.0 : cos_angle );
  float sin_angle = VectorUtils::det( t, sd.est_tgt );
  if ( sin_angle >= 0 ) angle_deviation = acos( cos_angle );
  else angle_deviation = - acos( cos_angle );
  
  // estimated angle to x-axis.
  float a_to_x;
  float cos_a = VectorUtils::dotProduct( Vector2D( 1.0, 0.0 ) , sd.est_tgt );
  cos_a = cos_a > 1.0 ? 1.0 : ( cos_a < -1.0 ? -1.0 : cos_a );
  float sin_a = VectorUtils::det( Vector2D( 1.0, 0.0 ) , sd.est_tgt );
  if ( sin_a >= 0 ) a_to_x = acos( cos_a );
  else a_to_x = 2 * M_PI - acos( cos_a );

  // again shift to fit the quadrant
  if ( fabs( first_angle - a_to_x) > M_PI ) 
    {
      if ( first_angle - a_to_x > 0.0 )
	a_to_x += 2.0 * M_PI; 
      else 
	a_to_x -= 2.0 * M_PI;
    }
  
  // theoretical angle to x-axis.
  float ta_to_x;
  float cos_ta = VectorUtils::dotProduct( Vector2D( 1.0, 0.0 ) , t );
  cos_ta = cos_ta > 1.0 ? 1.0 : ( cos_ta < -1.0 ? -1.0 : cos_ta );
  float sin_ta = VectorUtils::det( Vector2D( 1.0, 0.0 ) , t );
  if ( sin_ta >= 0 ) ta_to_x = acos( cos_ta );
  else ta_to_x = 2 * M_PI - acos( cos_ta );

  // again shift to fit the quadrant
  if ( fabs( first_angle - ta_to_x) > M_PI )
    {
      if ( first_angle - ta_to_x > 0.0 )
	ta_to_x += 2.0 * M_PI; 
      else 
	ta_to_x -= 2.0 * M_PI;
    }
  
  //      sd.surfel = bel;
  sd.theta = theta;
  sd.xc = xbel.ro( 0 );
  sd.yc = xbel.ro( 1 );
  sd.angle_deviation = angle_deviation;
  sd.angle_to_x = a_to_x;
  sd.th_angle_to_x = ta_to_x;
  sd.curvature = curv;
  //    cerr << "X=( " << xbel.ro( 0 ) << " " << xbel.ro( 1 )
  //         << " |" << ks.sorthDir( bel ) << "): "
  //         << "t= " << theta << " "
  //         << "tgt=( " << t.ro( 0 ) << " " << t.ro( 1 ) << " ) "
  //         << "curv= " << curv << " "
  //         << "est_tgt=( " << estimated_tgt.ro( 0 ) << " " << estimated_tgt.ro( 1 ) << " ) "
  //         << "angle_dev= " << angle_deviation << " "
  //         << "wdw= " << window << " "
  //         << "angle_to_x= " << a_to_x << " "
  //         << "th_angle_to_x= " << ta_to_x << " "
  //         << endl;
  
  // françois's changes
  sd.th_max_angle_to_x = tmax_a_to_x; 
  sd.th_min_angle_to_x = tmin_a_to_x; 
  sd.th_moy_angle_to_x = tmoy_a_to_x;
  // end françois's changes
}



struct PerShapeData
{
  uint nb_surfels;
  PerSurfelData* data;
  uint nb_filled;
  
  ~PerShapeData()
  {
    if ( data != 0 ) delete[] data;
  }
  
  PerShapeData()
    : nb_surfels( 0 ), data( 0 ), nb_filled( 0 )
  {
  }
  
  PerShapeData( uint nb )
  {
    data = new PerSurfelData[ nb ];
    nb_surfels = nb;
    nb_filled = 0;
  }

  PerShapeData( const PerShapeData& other )
  {
    data = new PerSurfelData[ other.nb_filled ];
    nb_surfels = other.nb_filled;
    for ( uint i = 0; i < other.nb_filled; ++i )
      data [ i ] = other.data[ i ];
    nb_filled = other.nb_filled;
  }
 
  PerShapeData& operator=( const PerShapeData& other )
  {
    if ( this != &other )
      {
	if ( nb_surfels < other.nb_filled )
	  {
	    if ( data != 0 ) delete[] data;
	    data = new PerSurfelData[ other.nb_filled ];
	    nb_surfels = other.nb_filled;
	  }
	for ( uint i = 0; i < other.nb_filled; ++i )
	  data[ i ] = other.data[ i ];
	nb_filled = other.nb_filled;
      }
    return *this;
  }
  
  
};


///////////////////////////////////////////////////////////////////////////////
// VARIABLES TO ANALYSE
///////////////////////////////////////////////////////////////////////////////

class PerSurfelValue 
{
public:
  virtual ~PerSurfelValue() {}
  
  virtual float value( const PerSurfelData & d ) const = 0;
};

class WindowValue : public PerSurfelValue
{
public:
  virtual float value( const PerSurfelData & d ) const
  {
    return d.window;
  }
  
};

class AngleToXValue : public PerSurfelValue
{
public:
  virtual float value( const PerSurfelData & d ) const
  {
    return d.angle_to_x;
  }
  
};


class AngleDeviationValue : public PerSurfelValue
{
public:
  virtual float value( const PerSurfelData & d ) const
  {
    return d.angle_deviation;
  }
  
};

class SquareAngleDeviationValue : public PerSurfelValue
{
public:
  virtual float value( const PerSurfelData & d ) const
  {
    return d.angle_deviation * d.angle_deviation;
  }
  
};

class CurvatureValue : public PerSurfelValue
{
public:
  virtual float value( const PerSurfelData & d ) const
  {
    return d.curvature;
  }
  
};

class AbsCurvatureValue : public PerSurfelValue
{
public:
  virtual float value ( const PerSurfelData & d ) const
  {
    return d.curvature >0 ? d.curvature : -d.curvature;
  }
};

class SquareCurvatureDeviationValue : public PerSurfelValue
{
public:
  virtual float value( const PerSurfelData & d ) const
  {
    return ( d.curvature - d.est_curv ) * ( d.curvature - d.est_curv );
  }
  
};

class EstimatedCurvatureValue : public PerSurfelValue
{
public:
  virtual float value( const PerSurfelData & d ) const
  {
    return d.est_curv;
  }
  
};

class RValue : public PerSurfelValue
{
public:
  virtual float value( const PerSurfelData & d ) const
  {
    return sqrt( d.xc * d.xc + d.yc * d.yc );
  }
  
};

class XCValue : public PerSurfelValue
{
public:
  virtual float value( const PerSurfelData & d ) const
  {
    return d.xc;
  }
  
};

class YCValue : public PerSurfelValue
{
public:
  virtual float value( const PerSurfelData & d ) const
  {
    return d.yc;
  }
  
};


class EstimatedNormalXValue : public PerSurfelValue
{
public:
  virtual float value( const PerSurfelData & d ) const
  {
    return -d.est_tgt.ro( 1 );
  }
  
};

class EstimatedNormalYValue : public PerSurfelValue
{
public:
  virtual float value( const PerSurfelData & d ) const
  {
    return d.est_tgt.ro( 0 );
  }
  
};

// françois's changes
class MinAngleToXValue : public PerSurfelValue
{
public:
  virtual float value (const PerSurfelData & d) const
  {
    return d.th_min_angle_to_x; 
  }

};

class MaxAngleToXValue : public PerSurfelValue
{
public:
  virtual float value (const PerSurfelData & d) const
  {
    return d.th_max_angle_to_x; 
  }

};

class MoyAngleToXValue : public PerSurfelValue
{
public:
  virtual float value (const PerSurfelData & d) const
  {
    return d.th_moy_angle_to_x; 
  }

};

class DevEstToMoyAngleValue : public PerSurfelValue 
{
public : 
  virtual float value (const PerSurfelData & d) const
  {
    float a1 = fabs(d.th_moy_angle_to_x - d.angle_to_x);
    return  a1;
  }
  
};

class SqDevEstToMoyAngleValue : public PerSurfelValue
{
public : 
  virtual float value (const PerSurfelData & d) const
  {
    float s1 = fabs(d.th_moy_angle_to_x - d.angle_to_x);
    
    return s1*s1; 
  }

};

class MaxDevEstToAngleValue : public PerSurfelValue
{
public : 
  virtual float value (const PerSurfelData & d) const
  {
    float a1 = fabs(d.th_min_angle_to_x - d.angle_to_x);
    float a2 = fabs(d.th_max_angle_to_x - d.angle_to_x);
    return (a1 >= a2 ? a1 : a2 );
  }

};

// end françois's changes



///////////////////////////////////////////////////////////////////////////////
// EXPERIMENTS
///////////////////////////////////////////////////////////////////////////////


class PerSurfelExperiment 
{
public:
  virtual ~PerSurfelExperiment() {}
  virtual void eval( KnSpace & ks, 
		     const C4CIteratorOnSurface* iter ) = 0;

};

class StarShapedPerSurfelExperiment : public PerSurfelExperiment
{
protected:
  const StarShaped & starshape;
  const Embedder & embedder;

  StarShapedPerSurfelExperiment( const StarShaped & ss,
				 const Embedder & emb )
    : starshape( ss ), embedder( emb )
  {
  }
  
};

class TangentComputationExperiment : public StarShapedPerSurfelExperiment
{
protected:
  TangentComputer & tangent_computer;
  PerShapeData shape_data;
  uint current_surfel;
  
public:
  TangentComputationExperiment( const StarShaped & ss,
				const Embedder & emb,
				TangentComputer & tgt_computer,
				uint nb_surfels )
    : StarShapedPerSurfelExperiment( ss, emb ),
      tangent_computer( tgt_computer ), 
      shape_data( nb_surfels ),
      current_surfel( 0 )
  {}

  const PerShapeData& getShapeData() const
  {
    return shape_data;
  }
  
  void setShapeData( const PerShapeData& other )
  {
    shape_data = other;
  }
  
  virtual void eval( KnSpace & ks, 
		     const C4CIteratorOnSurface* iter );

  void completeAllData( KnSpace & ks );

};


void 
TangentComputationExperiment::completeAllData( KnSpace & ks ) 
{
  for ( uint i = 0; i < shape_data.nb_filled; ++i )
    {
      completeSurfelData( ks, shape_data.data[ i ], starshape, embedder );
    } // for ( uint i = 0; i < shape_data.nb_filled; ++i )
}

// IMHO C LA QUI FAUT FAIRE DES CHANGEMENTS

void
TangentComputationExperiment::eval( KnSpace & ks, 
				    const C4CIteratorOnSurface* iter )
{
//   Vector xbel( 2 );
   Kn_sid bel = iter->current();
//   embedder.sembed( bel, xbel );
//   float theta = starshape.parameter( xbel );
//   Vector2D t = starshape.tangent( theta );
//   float curv = starshape.curvature( theta );
  float window = -1.0;
  float ds = 1.0;
  PerSurfelData sd;

  //   C4CIterator* iter_clone = iter->clone();
  //   C4CIteratorOnBdry *iter_copy;
  //   iter_copy = dynamic_cast<C4CIteratorOnBdry*>(iter_clone);
  //   delete iter_clone;
  
  Vector2D estimated_tgt 
    = tangent_computer.computeTangent( ks,
				       iter,
				       embedder,
				       window,
				       ds,
				       sd.ext_front, sd.ext_back );
//   // Angle deviation between estimated tangent and estimated tangent.
//   float angle_deviation;
//   float cos_angle = VectorUtils::dotProduct( t, estimated_tgt );
//   cos_angle = cos_angle > 1.0 ? 1.0 : cos_angle;
//   float sin_angle = VectorUtils::det( t, estimated_tgt );
//   if ( sin_angle >= 0 ) angle_deviation = acos( cos_angle );
//   else angle_deviation = - acos( cos_angle );

//   // estimated angle to x-axis.
//   float a_to_x;
//   float cos_a = VectorUtils::dotProduct( Vector2D( 1.0, 0.0 ) , estimated_tgt );
//   cos_a = cos_a > 1.0 ? 1.0 : ( cos_a < -1.0 ? -1.0 : cos_a );
//   float sin_a = VectorUtils::det( Vector2D( 1.0, 0.0 ) , estimated_tgt );
//   if ( sin_a >= 0 ) a_to_x = acos( cos_a );
//   else a_to_x = 2 * M_PI - acos( cos_a );

//   // theoretical angle to x-axis.
//   float ta_to_x;
//   float cos_ta = VectorUtils::dotProduct( Vector2D( 1.0, 0.0 ) , t );
//   cos_ta = cos_ta > 1.0 ? 1.0 : ( cos_ta < -1.0 ? -1.0 : cos_ta );
//   float sin_ta = VectorUtils::det( Vector2D( 1.0, 0.0 ) , t );
//   if ( sin_ta >= 0 ) ta_to_x = acos( cos_ta );
//   else ta_to_x = 2 * M_PI - acos( cos_ta );

  sd.surfel = bel;
//   sd.theta = theta;
//   sd.xc = xbel.ro( 0 );
//   sd.yc = xbel.ro( 1 );
  sd.est_tgt = estimated_tgt;
  sd.window = window;
  sd.est_ds = ds;
//   sd.angle_deviation = angle_deviation;
//   sd.angle_to_x = a_to_x;
//   sd.th_angle_to_x = ta_to_x;
//   sd.curvature = curv;
  shape_data.data[ current_surfel ] = sd;
  shape_data.nb_filled++;
  
//    cerr << "X=( " << xbel.ro( 0 ) << " " << xbel.ro( 1 )
//         << " |" << ks.sorthDir( bel ) << "): "
//         << "t= " << theta << " "
//         << "tgt=( " << t.ro( 0 ) << " " << t.ro( 1 ) << " ) "
//         << "curv= " << curv << " "
//         << "est_tgt=( " << estimated_tgt.ro( 0 ) << " " << estimated_tgt.ro( 1 ) << " ) "
//         << "angle_dev= " << angle_deviation << " "
//         << "wdw= " << window << " "
//         << "angle_to_x= " << a_to_x << " "
//         << "th_angle_to_x= " << ta_to_x << " "
//         << endl;

  current_surfel++;
}





///////////////////////////////////////////////////////////////////////////////
// HELPER CLASSES FOR STATISTICS
///////////////////////////////////////////////////////////////////////////////


class Statistics
{
  uint m_nb;
  uint* m_samples;
  float* m_exp;
  float* m_exp2;
  float* m_var;
  float* m_unbiased_var;
  float* m_max;
  float* m_min;

public:
  Statistics( uint size )
  {
    m_samples = 0;
    m_exp = 0;
    m_exp2 = 0;
    m_var = 0;
    m_unbiased_var = 0;
    m_max = 0;
    m_min = 0;
    init( size );
  }

  ~Statistics()
  {
    erase();
  }

  uint nb() const
  {
    return m_nb;
  }

  uint samples( uint k ) const
  {
    return m_samples[ k ];
  }

  float mean( uint k ) const
  {
    return m_exp[ k ];
  }

  float variance( uint k ) const
  {
    return m_var[ k ];
  }

  float unbiasedVariance( uint k ) const
  {
    return m_unbiased_var[ k ];
  }

  float max( uint k ) const
  {
    return m_max[ k ];
  }

  float min( uint k ) const
  {
    return m_min[ k ];
  }
  

  void addValue( uint k, float v )
  {
    m_samples[ k ] += 1;
    m_exp[ k ] += v;
    m_exp2[ k ] += v*v;
    if ( m_samples[ k ] == 1 )
      {
	m_max[ k ] = v;
	m_min[ k ] = v;
      }
    else if ( v > m_max[ k ] )
      m_max[ k ] = v;
    else if ( v < m_min[ k ] )
      m_min[ k ] = v;
  }
  
  void terminate()
  {
    for ( uint k = 0; k < m_nb; ++k )
      {
	m_exp[ k ] /= m_samples[ k ];
	m_exp2[ k ] /= m_samples[ k ];
	m_var[ k ] = m_exp2[ k ] - m_exp[ k ] * m_exp[ k ];
	m_unbiased_var[ k ] = m_samples[ k ] * m_var[ k ] 
	  / ( m_samples[ k ] - 1 );
      }
  }

  void init( uint size )
  {
    erase();
    m_nb = size;
    m_samples = new uint[ size ];
    m_exp = new float[ size ];
    m_exp2 = new float[ size ];
    m_var = new float[ size ];
    m_unbiased_var = new float[ size ];
    m_max = new float[ size ];
    m_min = new float[ size ];
    clear();
  }

  void clear()
  {
    if ( m_nb == 0 ) return;
    for ( uint i = 0; i < m_nb; ++ i )
      {
	m_samples[ i ] = 0;
	m_exp[ i ] = 0.0;
	m_exp2[ i ] = 0.0;
	m_var[ i ] = 0.0;
	m_unbiased_var[ i ] = 0.0;
	m_max[ i ] = 0.0;
	m_min[ i ] = 0.0;
      }
  }
  
  void erase() 
  {
    if ( m_samples != 0 ) delete[] m_samples;
    if ( m_exp != 0 ) delete[] m_exp;
    if ( m_exp2 != 0 ) delete[] m_exp2;
    if ( m_var != 0 ) delete[] m_var;
    if ( m_unbiased_var != 0 ) delete[] m_unbiased_var;
    if ( m_max != 0 ) delete[] m_max;
    if ( m_min != 0 ) delete[] m_min;
    m_samples = 0;
    m_exp = 0;
    m_exp2 = 0;
    m_var = 0;
    m_unbiased_var = 0;
    m_nb = 0;
    m_min = 0;
    m_max = 0;
  }
  
};



class PerThetaStatistic
{

  uint m_nb_trials;
  PerShapeData** m_all_datas;

public:  
  PerThetaStatistic( uint nb_trials )
  {
    m_all_datas = new PerShapeData*[ nb_trials ];
    m_nb_trials = nb_trials;
    for ( uint i = 0; i < nb_trials; ++i )
      m_all_datas[ i ] = new PerShapeData( 1024 );
  }

  ~PerThetaStatistic()
  {
    if ( m_all_datas != 0 )
      {
	for ( uint i = 0; i < m_nb_trials; ++i )
	  delete m_all_datas[ i ];
	delete[] m_all_datas;
      }
  }

  void assign( uint trial, const PerShapeData & data )
  {
    if ( trial < m_nb_trials )
      m_all_datas[ trial ]->operator=( data );
  }
  
  const PerShapeData & getData( uint trial ) const
  {
    return *( m_all_datas[ trial ] );
  }

  void statsForAValue( float theta0, float theta1, 
		       Statistics & stats,
		       const PerSurfelValue & variable )
  {
    stats.clear();
    uint nb_theta = stats.nb();
    
    // For all experiments
    for ( uint i = 0; i < m_nb_trials; ++i )
      {
	// For all theta angle on the shape.
	for ( uint j = 0; j < m_all_datas[ i ]->nb_filled; ++j )
	  {
	    PerSurfelData & surfdata = m_all_datas[ i ]->data[ j ];
	    if ( ( theta0 <= surfdata.theta ) && ( surfdata.theta < theta1 ) )
	      {
		uint k = (uint) floor( ( surfdata.theta - theta0 ) * nb_theta 
				       / ( theta1 - theta0 ) );
		stats.addValue( k, variable.value( surfdata ) );
	      }
	  }
      }
    stats.terminate();
  }
  
};



///////////////////////////////////////////////////////////////////////////////
// BLURRING
///////////////////////////////////////////////////////////////////////////////

/**
 * Abstract class to represent an indexed discrete signal.
 */
class DiscreteSignal 
{
public:
  virtual ~DiscreteSignal() {}
  virtual double value( int idx ) const = 0;
  virtual double distanceToNext( int idx ) const = 0;
  virtual double distanceToPrevious( int idx ) const = 0;
  
};

/**
 * Abstract class to represent an indexed weighted discrete signal.
 */
class WeightedDiscreteSignal : public DiscreteSignal
{
public:
  virtual ~WeightedDiscreteSignal() {}
  virtual double weight( int input_idx ) const = 0;

};


double
gaussianConvolution( const DiscreteSignal & dfct, const G & gfct, 
		     int idx, float band )
{
  // integration around 0.
  double support = 0.5 * ( dfct.distanceToNext( idx ) 
			  + dfct.distanceToPrevious( idx ) );
  double int_gauss = support * gfct( 0.0 );
  double acc = int_gauss;
  double value = int_gauss * dfct.value( idx );

  // integration on the positive half-space.
  //  cerr << "[band=" << band << " ";
  //cerr << "[" << "(" << idx << " " << dfct.value( idx ) << ") ";
  double x = 0.0;
  int pos_idx = idx;
  while ( true )
    {
      support = dfct.distanceToNext( pos_idx );
      x += support;
      // cerr << "(" << x << " " << pos_idx << ") ";
      if ( x > band ) break;
      //cerr << "(" << pos_idx << " " << dfct.value( pos_idx ) << ") ";
      ++pos_idx;
      support += dfct.distanceToNext( pos_idx );
      support *= 0.5;
      float int_gauss = support * gfct( -x );
      acc += int_gauss;
      value += int_gauss * dfct.value( pos_idx );
    }

  // integration on the negative half-space.
  x = 0.0;
  int neg_idx = idx;
  while ( true )
    {
      support = dfct.distanceToPrevious( neg_idx );
      x -= support;
      // cerr << "(" << x << " " << neg_idx << ") ";
      if ( x < -band ) break;
      //cerr << "(" << neg_idx << " " << dfct.value( neg_idx ) << ") ";
      --neg_idx;
      support += dfct.distanceToPrevious( neg_idx );
      support *= 0.5;
      float int_gauss = support * gfct( -x );
      acc += int_gauss;
      value += int_gauss * dfct.value( neg_idx );
    }
  //  cerr << "val=" << value << " acc=" << acc << endl;
  //cerr << "V=" << ( value / acc ) << " C=" << dfct.value( idx ) << endl;
  return value / acc;
}

double
weightedGaussianConvolution( const WeightedDiscreteSignal & dfct,
			     const G & gfct, 
			     int idx, float band )
{
  // integration around 0.
  double support = 0.5 * ( dfct.distanceToNext( idx ) 
			  + dfct.distanceToPrevious( idx ) );
  double int_gauss = support * gfct( 0.0 ) * dfct.weight( idx );
  double acc = int_gauss;
  double value = int_gauss * dfct.value( idx );

  // integration on the positive half-space.
  //  cerr << "[band=" << band << " ";
  //cerr << "[" << "(" << idx << " " << dfct.value( idx ) << ") ";
  double x = 0.0;
  int pos_idx = idx;
  while ( true )
    {
      support = dfct.distanceToNext( pos_idx );
      x += support;
      // cerr << "(" << x << " " << pos_idx << ") ";
      if ( x > band ) break;
      //cerr << "(" << pos_idx << " " << dfct.value( pos_idx ) << ") ";
      ++pos_idx;
      support += dfct.distanceToNext( pos_idx );
      support *= 0.5;
      float int_gauss = support * gfct( -x ) * dfct.weight( pos_idx );
      acc += int_gauss;
      value += int_gauss * dfct.value( pos_idx );
    }

  // integration on the negative half-space.
  x = 0.0;
  int neg_idx = idx;
  while ( true )
    {
      support = dfct.distanceToPrevious( neg_idx );
      x -= support;
      // cerr << "(" << x << " " << neg_idx << ") ";
      if ( x < -band ) break;
      //cerr << "(" << neg_idx << " " << dfct.value( neg_idx ) << ") ";
      --neg_idx;
      support += dfct.distanceToPrevious( neg_idx );
      support *= 0.5;
      float int_gauss = support * gfct( -x ) * dfct.weight( neg_idx );
      acc += int_gauss;
      value += int_gauss * dfct.value( neg_idx );
    }
  //  cerr << "val=" << value << " acc=" << acc << endl;
  //cerr << "V=" << ( value / acc ) << " C=" << dfct.value( idx ) << endl;
  return value / acc;
}


/**
 * Bridge class to see the shape data set as a discrete signal
 * sampling the angle between the contour and the x-axis.
 */ 
class ShapeDataDiscreteSignal : public WeightedDiscreteSignal
{
private:
  const PerShapeData* m_shape_data;
  uint m_base_idx;
  double m_base_angle_to_x;
  double m_dh;
  
public:
  ShapeDataDiscreteSignal() {}
  ~ShapeDataDiscreteSignal() {}
  
  void init( const PerShapeData* shape_data,
	     uint base_idx,
	     float dh );

  uint adjustIndex( int idx ) const;

  virtual double value( int idx ) const;
  virtual double distanceToNext( int idx ) const;
  virtual double distanceToPrevious( int idx ) const;
  virtual double weight( int input_idx ) const;
};


void
ShapeDataDiscreteSignal::init( const PerShapeData* shape_data,
			       uint base_idx,
			       float dh )
{
  m_shape_data = shape_data;
  m_base_idx = base_idx;
  m_base_angle_to_x = (double) m_shape_data->data[ m_base_idx ].angle_to_x;
  m_dh = (double) dh;
}

uint 
ShapeDataDiscreteSignal::adjustIndex( int input_idx ) const
{
  // Adjust index.
  while ( input_idx >= (int) m_shape_data->nb_filled ) 
    input_idx -= m_shape_data->nb_filled;
  while ( input_idx < 0 ) 
    input_idx += m_shape_data->nb_filled;
  return (uint) input_idx;
}

double 
ShapeDataDiscreteSignal::value( int input_idx ) const
{
  uint idx = adjustIndex( input_idx );

  // Adjust angle to relative base.
  double angle_to_x = (double) m_shape_data->data[ idx ].angle_to_x;
  while ( fabs( angle_to_x + 2*M_PI - m_base_angle_to_x ) 
	  < fabs( angle_to_x  - m_base_angle_to_x ) )
    angle_to_x += 2*M_PI;
  while ( fabs( angle_to_x - 2*M_PI - m_base_angle_to_x ) 
	  < fabs( angle_to_x  - m_base_angle_to_x ) )
    angle_to_x -= 2*M_PI;
  return angle_to_x;
}

double 
ShapeDataDiscreteSignal::distanceToNext( int input_idx ) const
{
  uint idx = adjustIndex( input_idx );
  uint idx1 = adjustIndex( input_idx + 1 );
  return 0.5 * m_dh * ( m_shape_data->data[ idx ].est_ds 
			+ m_shape_data->data[ idx1 ].est_ds );
}

double 
ShapeDataDiscreteSignal::distanceToPrevious( int input_idx ) const
{
  uint idx = adjustIndex( input_idx );
  uint idx1 = adjustIndex( input_idx - 1 );
  return 0.5 * m_dh * ( m_shape_data->data[ idx ].est_ds 
			+ m_shape_data->data[ idx1 ].est_ds );
}


double
ShapeDataDiscreteSignal::weight( int input_idx ) const
{
  uint idx = adjustIndex( input_idx );
  double weight = (double) m_shape_data->data[ idx ].window;
  return weight;
}



class BlurEstimatedTangent 
{

public:
  virtual ~BlurEstimatedTangent() {}
  
  virtual void blur( const PerShapeData & input, PerShapeData & output,
		     float dh ) = 0;
  
};

class FixedWindowBlurring : public BlurEstimatedTangent
{
  G m_gauss_fct;
  float m_sigma;
  float m_band;
  
public:
  virtual ~FixedWindowBlurring() {}
  FixedWindowBlurring( float sigma, float band );
  virtual void blur( const PerShapeData & input, PerShapeData & output,
		     float dh );

};

FixedWindowBlurring::FixedWindowBlurring( float sigma, float band )
  : m_gauss_fct( sigma ), m_sigma( sigma ), m_band( band )
{}

void 
FixedWindowBlurring::blur( const PerShapeData & input, PerShapeData & output,
			   float dh )
{
  ShapeDataDiscreteSignal shape_signal;
  // For all surfels
  for ( uint i = 0; i < input.nb_filled; ++i )
    {
      shape_signal.init( &input, i, dh );
      double blurred_angle_to_x 
	= gaussianConvolution( shape_signal, m_gauss_fct, i, m_band );
//       cerr << "Fixed [" << i << "]=" << blurred_angle_to_x 
// 	   << ", old:" << input.data[ i ].angle_to_x
// 	   << ", sigma:" << m_sigma << endl;
//       cerr << "[" << i << "]=" << blurred_angle_to_x 
// 	   << ", " << input.data[ i ].angle_to_x << endl;
//      output.data[ i ].angle_to_x = (float) blurred_angle_to_x;
      Vector2D new_est_tgt( cos( blurred_angle_to_x ),
			    sin( blurred_angle_to_x ) );
      output.data[ i ].est_tgt = new_est_tgt;
    }
}

class VariableWindowBlurring : public BlurEstimatedTangent
{
  G m_gauss_fct;
  float m_coef_sigma;
  float m_coef_band;
  
public:
  virtual ~VariableWindowBlurring() {}
  VariableWindowBlurring( float coef_sigma, float coef_band );
  virtual void blur( const PerShapeData & input, PerShapeData & output,
		     float dh );

};

VariableWindowBlurring::VariableWindowBlurring( float coef_sigma, 
						float coef_band )
  : m_gauss_fct(), 
    m_coef_sigma( coef_sigma ), m_coef_band( coef_band )
{}

void 
VariableWindowBlurring::blur( const PerShapeData & input,
			      PerShapeData & output,
			      float dh )
{
  ShapeDataDiscreteSignal shape_signal;
  // For all surfels
  for ( uint i = 0; i < input.nb_filled; ++i )
    {
      shape_signal.init( &input, i, dh );
      float sigma = m_coef_sigma * input.data[ i ].window;
      m_gauss_fct.setSigma( sigma );
      double blurred_angle_to_x 
	= gaussianConvolution( shape_signal, m_gauss_fct, i, 
			       m_coef_band * sigma );
//       cerr << "Var [" << i << "]=" << blurred_angle_to_x 
// 	   << ", old:" << input.data[ i ].angle_to_x
// 	   << ", sigma:" << sigma << endl;
//       output.data[ i ].angle_to_x = (float) blurred_angle_to_x;
      Vector2D new_est_tgt( cos( blurred_angle_to_x ),
			    sin( blurred_angle_to_x ) );
      output.data[ i ].est_tgt = new_est_tgt;
    }
}

class WeightedVariableWindowBlurring : public BlurEstimatedTangent
{
  G m_gauss_fct;
  float m_coef_sigma;
  float m_coef_band;
  
public:
  virtual ~WeightedVariableWindowBlurring() {}
  WeightedVariableWindowBlurring( float coef_sigma, float coef_band );
  virtual void blur( const PerShapeData & input, PerShapeData & output,
		     float dh );

};

WeightedVariableWindowBlurring::WeightedVariableWindowBlurring
( float coef_sigma, 
  float coef_band )
  : m_gauss_fct(), 
    m_coef_sigma( coef_sigma ), m_coef_band( coef_band )
{}

void 
WeightedVariableWindowBlurring::blur( const PerShapeData & input,
				      PerShapeData & output,
				      float dh )
{
  ShapeDataDiscreteSignal shape_signal;
  // For all surfels
  for ( uint i = 0; i < input.nb_filled; ++i )
    {
      shape_signal.init( &input, i, dh );
      float sigma = m_coef_sigma * input.data[ i ].window;
      m_gauss_fct.setSigma( sigma );
      double blurred_angle_to_x 
	= weightedGaussianConvolution( shape_signal, m_gauss_fct, i, 
				       m_coef_band * sigma );
//       cerr << "Var [" << i << "]=" << blurred_angle_to_x 
// 	   << ", old:" << input.data[ i ].angle_to_x
// 	   << ", sigma:" << sigma << endl;
//       output.data[ i ].angle_to_x = (float) blurred_angle_to_x;
      Vector2D new_est_tgt( cos( blurred_angle_to_x ),
			    sin( blurred_angle_to_x ) );
      output.data[ i ].est_tgt = new_est_tgt;
    }
}


///////////////////////////////////////////////////////////////////////////////
// CURVATURE COMPUTATION
///////////////////////////////////////////////////////////////////////////////

class CurvatureComputer
{
public:
  virtual ~CurvatureComputer() {}
  
  void computeCurvatures( const PerShapeData & shape_data,
				    double dh ) const;
  virtual double curvature( const PerShapeData & shape_data,
			    uint surfel_idx,
			    double dh ) const = 0;
};

void 
CurvatureComputer::computeCurvatures( const PerShapeData & shape_data,
				      double dh ) const
{
  for ( uint i = 0; i < shape_data.nb_filled; ++i )
    {
      double curv = curvature( shape_data, i, dh );
      shape_data.data[ i ].est_curv = curv;
    }
}


class CurvatureByDifferentiation : public CurvatureComputer
{
public:
  virtual ~CurvatureByDifferentiation() {}
  CurvatureByDifferentiation() {}
  
  virtual double curvature( const PerShapeData & shape_data,
			    uint surfel_idx,
			    double dh ) const;
};

double
CurvatureByDifferentiation::curvature
( const PerShapeData & shape_data,
  uint surfel_idx,
  double dh ) const
{
  uint next_idx = ( surfel_idx + 1 ) % shape_data.nb_filled;
  uint prev_idx = ( surfel_idx + shape_data.nb_filled - 1 ) 
    % shape_data.nb_filled;
  double dl = shape_data.data[ surfel_idx ].est_ds 
    + 0.5 * ( shape_data.data[ next_idx ].est_ds 
	      + shape_data.data[ prev_idx ].est_ds );

  float base_angle = shape_data.data[ surfel_idx ].angle_to_x;
  float next_angle = shape_data.data[ next_idx ].angle_to_x;
  while ( fabs( next_angle + 2*M_PI - base_angle ) 
	  < fabs( next_angle - base_angle ) )
    next_angle += 2*M_PI;
  while ( fabs( next_angle - 2*M_PI - base_angle ) 
	  < fabs( next_angle - base_angle ) )
    next_angle -= 2*M_PI;
  float prev_angle = shape_data.data[ prev_idx ].angle_to_x;
  while ( fabs( prev_angle + 2*M_PI - base_angle ) 
	  < fabs( prev_angle - base_angle ) )
    prev_angle += 2*M_PI;
  while ( fabs( prev_angle - 2*M_PI - base_angle ) 
	  < fabs( prev_angle - base_angle ) )
    prev_angle -= 2*M_PI;

  double curv = (double) (( next_angle - prev_angle ) / ( dh * dl ));
  return curv;
}

class CurvatureByWindowedDifferentiation : public CurvatureComputer
{
  float m_coef_wdw;
  
public:
  virtual ~CurvatureByWindowedDifferentiation() {}
  CurvatureByWindowedDifferentiation( float coef_wdw )
  : m_coef_wdw( coef_wdw )
  {}
  
  virtual double curvature( const PerShapeData & shape_data,
			    uint surfel_idx,
			    double dh ) const;
};

double
CurvatureByWindowedDifferentiation::curvature
( const PerShapeData & shape_data,
  uint surfel_idx,
  double dh ) const
{
  uint next_idx = surfel_idx;
  uint prev_idx = surfel_idx;
  float window = m_coef_wdw * shape_data.data[ surfel_idx ].window;
  float dist = 0.0;
  float geodist = 0.0;
  while ( geodist <= window )
    {
      Vector2D pos_oldnext ( shape_data.data[ next_idx ].xc,
			     shape_data.data[ next_idx ].yc );
      Vector2D pos_oldprev( shape_data.data[ prev_idx ].xc,
			    shape_data.data[ prev_idx ].yc );
      next_idx = ( next_idx + 1 ) % shape_data.nb_filled;
      prev_idx = ( prev_idx + shape_data.nb_filled - 1 ) 
	% shape_data.nb_filled;
      Vector2D pos_next( shape_data.data[ next_idx ].xc,
			 shape_data.data[ next_idx ].yc );
      Vector2D pos_prev( shape_data.data[ prev_idx ].xc,
			 shape_data.data[ prev_idx ].yc );

      pos_oldnext -= pos_next;
      pos_oldprev -= pos_prev;
      geodist += VectorUtils::norm( pos_oldnext ) + VectorUtils::norm( pos_oldprev );
      pos_next -= pos_prev;
      dist = VectorUtils::norm( pos_next );
      //cerr << next_idx << " " << prev_idx << " " << dist << endl;
    }

  float base_angle = shape_data.data[ surfel_idx ].angle_to_x;
  float next_angle = shape_data.data[ next_idx ].angle_to_x;
  while ( fabs( next_angle + 2*M_PI - base_angle ) 
	  < fabs( next_angle - base_angle ) )
    next_angle += 2*M_PI;
  while ( fabs( next_angle - 2*M_PI - base_angle ) 
	  < fabs( next_angle - base_angle ) )
    next_angle -= 2*M_PI;
  float prev_angle = shape_data.data[ prev_idx ].angle_to_x;
  while ( fabs( prev_angle + 2*M_PI - base_angle ) 
	  < fabs( prev_angle - base_angle ) )
    prev_angle += 2*M_PI;
  while ( fabs( prev_angle - 2*M_PI - base_angle ) 
	  < fabs( prev_angle - base_angle ) )
    prev_angle -= 2*M_PI;
  
  double curv = (double) ( ( next_angle - prev_angle ) / geodist );
  return curv;
}


class CurvatureByCircumscribedCircle : public CurvatureComputer
{
public:
  virtual ~CurvatureByCircumscribedCircle() {}
  CurvatureByCircumscribedCircle() {}
  
  virtual double curvature( const PerShapeData & shape_data,
			    uint surfel_idx,
			    double dh ) const;
};

double
CurvatureByCircumscribedCircle::curvature
( const PerShapeData & shape_data,
  uint surfel_idx,
  double dh ) const
{
  Vector2D A( shape_data.data[ surfel_idx ].xc,
	    shape_data.data[ surfel_idx ].yc );
  Vector2D B( shape_data.data[ surfel_idx ].ext_front );
  Vector2D C( shape_data.data[ surfel_idx ].ext_back );
  
  double curv = 
    (double) EuclideanGeometry::curvatureCircumscribedCircle
    ( A.x(), A.y(), B.x(), B.y(), C.x(), C.y() );
  if ( curv < 0.0 ) curv = 0.0;
  Vector2D AB( B ); AB -= A;
  Vector2D CA( A ); CA -= C;
  if ( VectorUtils::det( CA, AB ) < 0 ) curv = -curv;
  if ( surfel_idx == 1 )
    {
      cerr << "[DBG] curv=" << curv 
	   << " A=(" << A.ro( 0 ) << " " << A.ro( 1 ) << ")"
	   << " B=(" << B.ro( 0 ) << " " << B.ro( 1 ) << ")"
	   << " C=(" << C.ro( 0 ) << " " << C.ro( 1 ) << ")"
	   << endl;
    }
  
  return curv;
}


///////////////////////////////////////////////////////////////////////////////
// INTEGRATION

void
integrateCurve( float theta0, float theta1,
		const Statistics & curv, 
		const Statistics & normal_x,
		const Statistics & normal_y,
		float r0, float tgt_x0, float tgt_y0 )
{
  uint nb_thetas = curv.nb();
  float dtheta = ( theta1 - theta0 ) / nb_thetas;
  float theta = theta0;
  float tx = tgt_x0;
  float ty = tgt_y0;
  float r = r0;
  for ( uint k = 0; k < nb_thetas; ++k )
    {
      float norm = sqrt( tx * tx + ty * ty );
      tx /= norm;
      ty /= norm;
      float c = curv.mean( k );
      float nx = normal_x.mean( k );
      float ny = normal_y.mean( k );
      float x = r * cos( theta );
      float y = r * sin( theta );
      float h = -r + r * ( tx * sin( theta ) - ty * cos( theta ) )
	/ ( tx * sin( theta + dtheta ) - ty * cos( theta + dtheta ) );
      float ds = sqrt( h * h + r * r * dtheta * dtheta );
      cout << k << " " << theta << " " << h << " " << ds << " "
	   << x << " " << y << " " << tx << " " << ty 
	   << " " << nx << " " << ny  << " " << c << endl;
      tx += c * nx * ds;
      ty += c * ny * ds;
      r += h;
      theta += dtheta;
    }
}

void
integrateCurve2( float theta0, float theta1,
		 const Statistics & true_x, 
		 const Statistics & true_y, 
		 const Statistics & curv, 
		 const Statistics & normal_x,
		 const Statistics & normal_y,
		 float r0, float tgt_x0, float tgt_y0 )
{
  uint nb_thetas = curv.nb();
  float dtheta = ( theta1 - theta0 ) / nb_thetas;
  float theta = theta0;
  float tx = tgt_x0;
  float ty = tgt_y0;
  float r = r0;
  float x = r * cos( theta );
  float y = r * sin( theta );
  float norm = sqrt( tx * tx + ty * ty );
  tx /= norm;
  ty /= norm;
  for ( uint k = 0; k < nb_thetas; ++k )
    {
      float c = curv.mean( k );
      float nx = normal_x.mean( k );
      float ny = normal_y.mean( k );
      float ds = sqrt( ( true_x.mean( k ) - true_x. mean( k + 1 ) )
		       * ( true_x.mean( k ) - true_x. mean( k + 1 ) )
		       + ( true_y.mean( k ) - true_y. mean( k + 1 ) )
		       * ( true_y.mean( k ) - true_y. mean( k + 1 ) ) );
      cout << k << " " << theta << " " << "(x,y)" << " " << ds << " "
	   << x << " " << y << " " << tx << " " << ty 
	   << " " << nx << " " << ny  << " " << c << endl;
      tx += c * nx * ds;
      ty += c * ny * ds;
      norm = sqrt( tx * tx + ty * ty );
      tx /= norm;
      ty /= norm;
      x += tx * ds;
      y += ty * ds;
      theta += dtheta;
    }
}



void 
forallSurfels( KnSpace & ks, KnCharSet voxset, Kn_sid starting_bel,
	       PerSurfelExperiment & experiment)
{
  // Get tangent plane.
  BelAdjacency badj( ks, true );
  ObjectBoundary bdry( badj, voxset );
  uint track_dir = *( ks.sbegin_dirs( starting_bel ) );
  
  // C4CIteratorOnSurface* cp = bdry.newC4CIterator( starting_bel, track_dir );
  C4CIteratorOnBdry* cp = new C4CIteratorOnBdry( badj, starting_bel, track_dir, voxset );

  Kn_sid bel = starting_bel;
  Vector xbel( 2 );
  uint nb = 0;
  
  do 
    {
      
//       if ( K2SpaceViewer::viewer() != 0 )
// 	{
// 	  K2SpaceViewer::viewer()->addSCell( bel,
// 					     Color( 1.0, 0.0, 0.0 ) );
// 	}

      // Make experiment.
      experiment.eval( ks, cp );
      
      // Go to next one.
      cp->next();
      bel = cp->current();
      ++nb;
    }
  while ( bel != starting_bel );

  
}


//void experiment1( KnSpace & ks, const Embedder & embedder )


int
main( int argc, char** argv ) 
{
  StandardArguments::addDigitalArgs( args, 2 );
  ShapeHelper::addStarShapedArgs( args );
  //  args.addOption( "-view", "-view <w> <h> <z>: open a window of size <w>x<h> to display some results with grid size of <z> pixels.", "512", "512", "8" );
  // args.addOption( "-step", "-step <h>: the discretization step or resolution, the closer to 0, the finer.", "1.0" );
  // args.addOption( "-circle", "-circle <R>: the test shape is a circle of radius R", "10.0" );
  // args.addOption( "-flower", "-flower <R> <r> <k> <phi>: the test shape is a flower with k extremeties with mean radius R and variability of radius r, phi is the phase of the flower", "3", "10.0", "5.0", "0.0" );
  // args.addOption( "-accflower", "-accflower <R> <r> <k>: the test shape is a phase accelerating flower with k extremeties with mean radius R and variability of radius r", "4.0", "2.0", "4" );
  // args.addOption( "-rsquare", "-rsquare <R> <r>: the test shape is rounded square of big radius R and small corner radius r", "4.0", "1.0" );
    // args.addBooleanOption( "-fuzzy_uniform", "-fuzzy_uniform: possible uniform random shift of the shape center." );

  args.addOption( "-nbtrials", "-nbtrials <N>: the number of experiments over which the result is averaged", "1" );

  args.addOption( "-tgt_computer", "-tgt_computer {TDE|TDS|HT|MS|MST|MSA|MSTA|MATAS95|GAUSSDERIV|IPF|REGLIN} <n::integer> <m::integer>: the algorithm used to compute the tangent around a contour point, <n> is window size, <m> is a user parameter, \n\t\tTDS is symmetric discrete tangent, \n\t\tTDE is extended discrete tangent, \n\t\tHT is median of half-tangents, \n\t\tMS is interpolation of maximal segments with bellshape fct, \n\t\tMST is interpolation of maximal segments with triangle fct, \n\t\tMSA is same as MS but use point interpolation, \n\t\tMSTA is same as MST but use point interpolation, \n\t\tMATAS95 uses a median filtering technique (n: window 2n+1) (m:0 is std method, 1 and 2 are biaised), \n\t\tGAUSSDERIV uses gaussian derivative, \n\t\tHGAUSSDERIV uses an hybrid maxDSS gaussian derivate, the size of the computation is the max of (|c-b|,|f-c|) were c is the point of interest, b the BACK of the BACK MS and f the FRONT of the FRONT MS,  \n\t\tH2GAUSSDERIV uses an hybrid maxDSS gaussian derivate, the size of the computation is the max of (|c-b|,|f-c|)^(3/2), \n\t\tH3GAUSSDERIV uses an hybrid maxDSS gaussian derivate, the size of the computation window is the average length of max dss,\n\t\tIPF stands for implicit parabola fitting, \n\t\tEPF stands for explicit parabola fitting, \n\t\tSIPF stands for separated implicit parabola fitting, ie implicit parabola fittin on each coordinates, \n\t\tSEPF stands for separated explicit parabola fitting, ie explicit parabola fitting on each coordinates, \n\t\tHSIPF stands for hybrid separated implicit parabola fitting, ie implicit parabola fittin on each coordinates with the window size determined with MS \n\t\tH3SIPF uses an hybrid maxDSS SIPF, the size of the computation window is the average length of max dss,\n\t\tREGLIN stands for linear regression; \n\t\t<n> when 0, these tangents are not limited (except for MATAS95), otherwise they are limited to 2n+1 surfels \n\t\tTEST stands a test with constant value, used for time purposes, makes a loop of size <n>", "TDS", "0", "0" );

  args.addOption( "-thetas", "-thetas <t0> <t1> <nb>: the computations are gathered into <nb> different angle theta intervals from <t0> to <t1>", "0", "6.28318530717958647688", "144" );
  
  args.addBooleanOption( "-v_r", "-v_r: computes the distance to the origin (for displaying the starshape with gnuplot." );

  args.addBooleanOption( "-v_angle_dev", "-v_angle_dev: computes the statistic of angle deviation between true tangent and estimated tangent." );

  args.addBooleanOption( "-v_sqr_angle_dev", "-v_sqr_angle_dev: computes the statistic of angle deviation between true tangent and estimated tangent." );

  args.addBooleanOption( "-v_window", "-v_window: computes the statistic of tangent window." );

  args.addBooleanOption( "-v_angle_to_x", "-v_angle_to_x: computes the statistic of the angle between the estimated tangent and the x-axis." );

  // françois's changes
  args.addBooleanOption( "-v_min_angle_to_x", "-v_min_angle_to_x: computes the statistics of the minimal angle between the theoretical tangent and the x-axis along the surfel's projection on the bound of the shape.");
  
  args.addBooleanOption( "-v_max_angle_to_x", "-v_max_angle_to_x: computes the statistics of the maximal angle between the theoretical tangent and the x-axis along the surfel's projection on the bound of the shape.");
  
  args.addBooleanOption( "-v_moy_angle_to_x", "-v_moy_angle_to_x: computes the statistics of the average angle between the theoretical tangent and the x-axis along the surfel's projection on the bound of the shape.");
  
  args.addBooleanOption( "-v_dev_moy_angle", "-v_dev_moy_angle: computes the statistics of the deviation between estimated angle to the x-axis and average angle between the theoretical tangent and the x-axis along the surfel's projection on the bound of the shape.");

  args.addBooleanOption( "-v_sq_dev_moy_angle", "-v_sq_dev_moy_angle: computes the statistics of the squared deviation between estimated angle to the x-axis and average angle between the theoretical tangent and the x-axis along the surfel's projection on the bound of the shape.");

  args.addBooleanOption( "-v_dev_max_angle", "-v_dev_max_angle: computes the statistics of the maximum deviation between estimated angle to the x-axis and  angle between the theoretical tangent and the x-axis along the surfel's projection on the bound of the shape.");

  args.addBooleanOption( "-v_abs_curvature", "-v_abs_curvature: computes the absolute value of the shape curvature (not really a statistic, but useful for comparison purposes." );

  args.addBooleanOption( "-time", "-time: compute time spent on each xp ." );
  
// end françois's changes

  args.addBooleanOption( "-v_curvature", "-v_curvature: computes the shape curvature (not really a statistic, but useful for comparison purposes." );

  args.addBooleanOption( "-v_estcurv", "-v_estcurv: estimates the shape curvature." );


  args.addBooleanOption( "-v_sqr_curv_dev", "-v_sqr_curv_dev: computes the statistic of squared curvature deviation." );

  args.addOption( "-curv_computer", "-curv_computer {NO|DIFF|WDIFF|CCIRCLE} <n>: the algorithm used to compute the discrete curvature around a contour point, {DIFF} is simple differentiation and <n> is reserved, {WDIFF} is window differentiation and <n> is the window multiplier, {CCIRCLE} is circumscribed circle and <n> is reserved.", "NO", "1.0" );


  args.addOption( "-blur", "-blur {no|fixedw|variablew|weightedw} <c_sigma> <c_band>: blur the tangent estimation according to the parameters: {no} does nothing, {fixed} blur with a gaussian of fixed window size s=c_sigma*dh and the computation is limited to the band <c_band>*s (use c_band>=4), {variablew} blur with a gaussian of variable window size s=c_sigma*window and the computation is limited to the band <c_band>*s (use c_band>=4), {weightedw} same as 'variablew' except that convolution is weighted by the local window size.", "no", "1.0", "4.0" );

  args.addOption( "-integrate_curve", "-integrate_curve <r> <tx> <ty>: estimates the curve from its estimated curvatures and normal vectors, <r> is the initial distance to the origin, (<tx>,<ty>) is the first tangent vector", "4.0", "1.0", "0.0" );

  args.addOption( "-integrate_curve2", "-integrate_curve2 <r> <tx> <ty>: estimates the curve from its estimated curvatures and normal vectors, <r> is the initial distance to the origin, (<tx>,<ty>) is the first tangent vector", "4.0", "1.0", "0.0" );

  args.addOption( "-integrate_truecurve2", "-integrate_curve2 <r> <tx> <ty>: estimates the curve from its true curvatures and normal vectors, <r> is the initial distance to the origin, (<tx>,<ty>) is the first tangent vector", "4.0", "1.0", "0.0" );


//   args.addBooleanOption( "-per_surfel", 
// 			 "-per_surfel: asks to display the curvature for each surfel" );
//   args.addBooleanOption( "-curv_windowed_angle", 
// 			 "-curv_windowed_angle: computes curvature by a variation angle defined over an adaptative window" );

  
  if ( ( argc <= 1 ) 
       || ! args.readArguments( argc, argv ) ) 
    {
      cerr << args.usage( "test_Tangent", 
			  "Tests several discrete tangent estimators and output some statistics.",
			  "-h -x -y -fuzzy_uniform -nbtrials -circle -flower -accflower -rsquare -step -thetas -v_angle_dev -v_sqr_angle_dev -v_window -v_angle_to_x -v_min_angle_to_x -v_max_angle_to_x -v_moy_angle_to_x -v_dev_moy_angle -v_sq_dev_moy_angle -v_dev_max_angle -v_abs_curvature-v_curvature -tgt_computer -blur -v_estcurv -v_sqr_curv_dev -v_r -curv_computer -integrate_curve -integrate_curve2 -integrate_truecurve2 -time" ) 
	   << endl;
      return 1;
    }


  // Display command line.
  cout << "#";
  for ( int i = 0; i < argc; ++i )
    cout << " " << argv[ i ];
  cout << endl;

  // -------------------------------------------------------------------------
  // Build space.
  uint d = StandardArguments::dim( args );
  if ( d != 2 )
    {
      cerr << "Dimension is 2." << endl;
      return 2;
    }
  Kn_size sizes[ d ];
  StandardArguments::fillSizes( args, sizes );
  KnSpace ks( 2, sizes );
  Kn_uid vcenter = ks.uspel( ks.ukcode( sizes ) ); // center
  Vector xcenter = ks.ucentroid( vcenter );
  
  cerr << "--- Space: " << ks << endl;

  // -------------------------------------------------------------------------
  // Take care of visualization.
  //  SWindow* win = 0;
//   if ( args.check( "-view" ) )
//     {
//       K2SpaceViewer* visu2d = new K2SpaceViewer;
//       win = new RGBWindow( args.getOption( "-view" )->getIntValue( 0 ),
// 			   args.getOption( "-view" )->getIntValue( 1 ) );
      
//       visu2d->init( win, &ks );
//       visu2d->setSizes( args.getOption( "-view" )->getIntValue( 2 ),
// 			args.getOption( "-view" )->getIntValue( 2 ),
// 			2, 2 );
//       visu2d->clear();
//       K2SpaceViewer::setViewer( visu2d );
//     }
  

  // -------------------------------------------------------------------------
  // Get some parameters.
  float dh = args.getOption( "-step" )->getFloatValue( 0 );
  uint nb_trials = args.getOption( "-nbtrials" )->getIntValue( 0 );

  cerr << "Embedder" << endl;
  GridEmbedder embedder;
  embedder.init( &ks );
  embedder.setCenter( xcenter );
  embedder.setScale( dh );

  // Tangent computers.
  TangentComputer* tgt_computer = 0;
  std::string tgt_calc_str = args.getOption( "-tgt_computer" )->getValue( 0 );
  uint tgt_calc_winsize = args.getOption( "-tgt_computer" )->getIntValue( 1 );
  uint tgt_calc_type = args.getOption( "-tgt_computer" )->getIntValue( 2 );
  if ( tgt_calc_str == "TEST" )
    tgt_computer = new TestTimeTangentComputer;
  else if ( tgt_calc_str == "TDE" )
    tgt_computer = new EDTangentComputer;
  else if ( tgt_calc_str == "TDS" )
    tgt_computer = new SDTangentComputer;
  else if ( tgt_calc_str == "HT" )
    tgt_computer = new HTangentComputer;
  else if ( tgt_calc_str == "MS" )
    tgt_computer = new MSTangentComputer( MSTangentComputer::BELLSHAPE2, 
					  false );
  else if ( tgt_calc_str == "MST" )
    tgt_computer = new MSTangentComputer( MSTangentComputer::TRIANGLE, 
					  false );
  else if ( tgt_calc_str == "IMST" )
    tgt_computer = new PreMSTangentComputer( PreMSTangentComputer::TRIANGLE);
  else if ( tgt_calc_str == "MSA" )
    tgt_computer = new MSTangentComputer( MSTangentComputer::BELLSHAPE2,
					  true );
  else if ( tgt_calc_str == "MSTA" )
    tgt_computer = new MSTangentComputer( MSTangentComputer::TRIANGLE,
					  true );
  else if ( tgt_calc_str == "MATAS95" ) {
    tgt_computer = new Matas95TangentComputer;
    tgt_computer->setType(tgt_calc_type);
  }

  else if ( tgt_calc_str == "GAUSSDERIV" )
    tgt_computer = new GaussianDerivativeTangentComputer;
  else if ( tgt_calc_str == "HGAUSSDERIV" )
    tgt_computer = new HybridGaussianDerivativeTangentComputer;
  else if ( tgt_calc_str == "H2GAUSSDERIV" )
    tgt_computer = new Hybrid2GaussianDerivativeTangentComputer;
  else if ( tgt_calc_str == "H3GAUSSDERIV" )
    tgt_computer = new Hybrid3GaussianDerivativeTangentComputer;
  else if ( tgt_calc_str == "IPF" )
    tgt_computer = new ImplicitParabolaFittingTangentComputer;
  else if ( tgt_calc_str == "PREIPF" )
    tgt_computer = new PreImplicitParabolaFittingTangentComputer;
  else if ( tgt_calc_str == "EPF" )
    tgt_computer = new ExplicitParabolaFittingTangentComputer;
  else if ( tgt_calc_str == "SIPF" )
    tgt_computer = new SeparatedImplicitParabolaFittingTangentComputer;
  else if ( tgt_calc_str == "SEPF" )
    tgt_computer = new SeparatedExplicitParabolaFittingTangentComputer;
  else if ( tgt_calc_str == "HSIPF" )
    tgt_computer = new HybridSeparatedImplicitParabolaFittingTangentComputer;
  else if ( tgt_calc_str == "H3SIPF" )
    tgt_computer = new Hybrid3SeparatedImplicitParabolaFittingTangentComputer;
  else if ( tgt_calc_str == "REGLIN" )
    tgt_computer = new RegLinTangentComputer;
  else
    tgt_computer = new SDTangentComputer;
  tgt_computer->size = 
    (uint) args.getOption( "-tgt_computer" )->getIntValue( 1 );
  tgt_computer->grid_step  = 
    (float) args.getOption( "-step")->getFloatValue( 0 );

  // ajouter le cas ou MATAS95 est choisi pour mettre tgt_computer->size a 3...
  
  // Blurring techniques
  BlurEstimatedTangent* blurrer = 0;
  std::string blur_str = args.getOption( "-blur" )->getValue( 0 );
  float coef_sigma = args.getOption( "-blur" )->getFloatValue( 1 );
  float coef_band = args.getOption( "-blur" )->getFloatValue( 2 );
  // {no|fixedw|variablew|weightedw}
  if ( blur_str == "fixedw" )
    blurrer = new FixedWindowBlurring( coef_sigma * dh, 
				       coef_sigma * dh * coef_band );
  else if ( blur_str == "variablew" )
    blurrer = new VariableWindowBlurring( coef_sigma, 
					  coef_band );
  else if ( blur_str == "weightedw" )
    blurrer = new WeightedVariableWindowBlurring( coef_sigma, 
						  coef_band );

  // Curvature computers.
  CurvatureComputer* curv_computer = 0;
  std::string curv_calc_str = args.getOption( "-curv_computer" )->getValue( 0 );
  if ( curv_calc_str == "DIFF" )
    curv_computer = new CurvatureByDifferentiation;
  else if ( curv_calc_str == "WDIFF" )
    curv_computer = new CurvatureByWindowedDifferentiation( args.getOption( "-curv_computer" )->getFloatValue( 1 ) );
  else if ( curv_calc_str == "CCIRCLE" )
    curv_computer = new CurvatureByCircumscribedCircle;

  // Starting series of experiments.
  PerThetaStatistic stat_storer( nb_trials );

  long my_time;
  bool compute_xp_time;
    
  if (args.check( "-time"))
    { 
      my_time = 0; 
      compute_xp_time = true;
    }
  else 
    compute_xp_time = false;
  
  for ( uint n = 0; n < nb_trials; ++n )
    {
      cerr <<" ---------------- Trial " << n << " ---------------" << endl;

      // Make shape.
      uint nb_bels = 0;
      Kn_sid bel = 0;
      KnCharSet voxset = KnCharSet::create( ks, ks.dim(), false, 0 );
      StarShaped* shape = shapeFromArgs( ks, embedder, voxset, bel, nb_bels );

//       if ( K2SpaceViewer::viewer() != 0 )
// 	{
// 	  K2SpaceViewer::viewer()->clear();
// 	  K2SpaceViewer::viewer()->addKnCharSet( voxset, false, 
// 						 Color( 0.3, 0.3, 0.3 ) );
// 	}
      


      // Prepare and perform experiment.
      TangentComputationExperiment experiment( *shape, embedder, 
					       *tgt_computer, 
					       nb_bels );

      if (compute_xp_time)
	Clock::startClock();      

      forallSurfels( ks, voxset, bel, experiment);
      
      if (compute_xp_time)
	my_time += Clock::stopClock();

      // reset estimator 
      tgt_computer->reset();
      
      experiment.completeAllData( ks );

      // Blur if needed.
      if ( blurrer != 0 )
	{
	  PerShapeData blurred_shape_data = experiment.getShapeData();
	  blurrer->blur( experiment.getShapeData(), 
			 blurred_shape_data,
			 dh );
	  experiment.setShapeData( blurred_shape_data );
	  experiment.completeAllData( ks );
	}

      // Estimates curvatures if needed.
      if ( curv_computer != 0 )
	{
	  cerr << "Curvatures" << endl;
	  PerShapeData curv_shape_data = experiment.getShapeData();
	  curv_computer->computeCurvatures( curv_shape_data, dh );
	  experiment.setShapeData( curv_shape_data );
	}
      
      // Memorize computations.
      stat_storer.assign( n, experiment.getShapeData() );
      
      // Free some stuff
      if ( shape != 0 ) delete shape;

//       if ( K2SpaceViewer::viewer() != 0 )
// 	{
// 	  K2SpaceViewer::viewer()->view();
// 	}
    } // for ( uint n = 0; n < nb_trials; ++n )

  // Get variable over which to perform stat.
  PerSurfelValue* variable = 0;
  if ( args.check( "-v_angle_dev" ) )
    variable = new AngleDeviationValue;
  else if ( args.check( "-v_sqr_angle_dev" ) )
    variable = new SquareAngleDeviationValue;
  else if ( args.check( "-v_window" ) )
    variable = new WindowValue;
  else if ( args.check( "-v_angle_to_x" ) )
    variable = new AngleToXValue;
  else if ( args.check( "-v_curvature" ) )
    variable = new CurvatureValue;
  else if ( args.check( "-v_abs_curvature" ) )
    variable = new AbsCurvatureValue;
  else if ( args.check( "-v_estcurv" ) )
    variable = new EstimatedCurvatureValue;
  else if ( args.check( "-v_sqr_curv_dev" ) )
    variable = new SquareCurvatureDeviationValue;
  else if ( args.check( "-v_r" ) )
    variable = new RValue;
  // françois's changes
  else if ( args.check( "-v_min_angle_to_x" ) )
    variable = new MinAngleToXValue; 
  else if ( args.check( "-v_max_angle_to_x" ) )
    variable = new MaxAngleToXValue; 
  else if ( args.check( "-v_moy_angle_to_x" ) )
    variable = new MoyAngleToXValue; 
  else if ( args.check( "-v_dev_moy_angle" ) )
    variable = new DevEstToMoyAngleValue;
  else if ( args.check( "-v_sq_dev_moy_angle" ) )
    variable = new SqDevEstToMoyAngleValue;
  else if ( args.check( "-v_dev_max_angle" ) )
    variable = new MaxDevEstToAngleValue;
  // end françois's changes

  float theta0 = args.getOption( "-thetas" )->getFloatValue( 0 );
  float theta1 = args.getOption( "-thetas" )->getFloatValue( 1 );
  uint nb_thetas = args.getOption( "-thetas" )->getIntValue( 2 );
  //  cerr << "Variable=" << variable << endl;
  
  if ( variable != 0 )
    if ( nb_thetas != 0 )
      {
	cout << "# Per theta [tmin:tmax] range data" << endl;
	cout << "# idx tmin #X E(X) var(X) uvar(X) max(X) min(X)" << endl;
	Statistics stat_per_theta( nb_thetas );
	stat_storer.statsForAValue( theta0, theta1, stat_per_theta, *variable );
	for ( uint k = 0; k < nb_thetas; ++k )
	  cout << k << " " 
	       << ( theta0 + k * ( theta1 - theta0 ) / nb_thetas ) << " "
	       << stat_per_theta.samples( k ) << " "
	       << stat_per_theta.mean( k ) << " "
	       << stat_per_theta.variance( k ) << " "
	       << stat_per_theta.unbiasedVariance( k ) << " "
	       << stat_per_theta.max( k ) << " "
	       << stat_per_theta.min( k ) << " "
	       << endl;
      }
    else // nb_thetas = 0
      {
	cout << "# Per surfel data" << endl;
	cout << "# idx t 1 value true_value xc yc" << endl;
	const PerShapeData & shapedata = stat_storer.getData( 0 );
	for ( uint k = 0; k < shapedata.nb_filled; ++k )
	  {
	    PerSurfelData & surfdata = shapedata.data[ k ];
	    cout << k << " "
		 << surfdata.theta << " "
		 << "1 "
		 << variable->value( surfdata ) << " "
		 << surfdata.th_angle_to_x << " "
		 << surfdata.xc << " "
		 << surfdata.yc << " "
		 << endl;
	  }
      }
  
  if ( args.check( "-integrate_curve" ) )
    {
      Statistics curv( nb_thetas );
      Statistics normal_x( nb_thetas );
      Statistics normal_y( nb_thetas );
      stat_storer.statsForAValue( theta0, theta1, curv, 
				  EstimatedCurvatureValue() );
      stat_storer.statsForAValue( theta0, theta1, normal_x, 
				  EstimatedNormalXValue() );
      stat_storer.statsForAValue( theta0, theta1, normal_y, 
				  EstimatedNormalYValue() );
      integrateCurve( theta0, theta1,
		      curv, normal_x, normal_y,
		      args.getOption( "-integrate_curve" )->getFloatValue( 0 ),
		      args.getOption( "-integrate_curve" )->getFloatValue( 1 ),
		      args.getOption( "-integrate_curve" )->getFloatValue( 2 )
		      );
    }
  if ( args.check( "-integrate_curve2" ) )
    {
      Statistics x( nb_thetas );
      Statistics y( nb_thetas );
      Statistics curv( nb_thetas );
      Statistics normal_x( nb_thetas );
      Statistics normal_y( nb_thetas );
      stat_storer.statsForAValue( theta0, theta1, x, 
				  XCValue() );
      stat_storer.statsForAValue( theta0, theta1, y, 
				  YCValue() );
      stat_storer.statsForAValue( theta0, theta1, curv, 
				  EstimatedCurvatureValue() );
      stat_storer.statsForAValue( theta0, theta1, normal_x, 
				  EstimatedNormalXValue() );
      stat_storer.statsForAValue( theta0, theta1, normal_y, 
				  EstimatedNormalYValue() );
      integrateCurve2
	( theta0, theta1,
	  x, y, curv, normal_x, normal_y,
	  args.getOption( "-integrate_curve2" )->getFloatValue( 0 ),
	  args.getOption( "-integrate_curve2" )->getFloatValue( 1 ),
	  args.getOption( "-integrate_curve2" )->getFloatValue( 2 )
	  );
    }
  if ( args.check( "-integrate_truecurve2" ) )
    {
      Statistics x( nb_thetas );
      Statistics y( nb_thetas );
      Statistics curv( nb_thetas );
      Statistics normal_x( nb_thetas );
      Statistics normal_y( nb_thetas );
      stat_storer.statsForAValue( theta0, theta1, x, 
				  XCValue() );
      stat_storer.statsForAValue( theta0, theta1, y, 
				  YCValue() );
      stat_storer.statsForAValue( theta0, theta1, curv, 
				  CurvatureValue() );
      stat_storer.statsForAValue( theta0, theta1, normal_x, 
				  EstimatedNormalXValue() );
      stat_storer.statsForAValue( theta0, theta1, normal_y, 
				  EstimatedNormalYValue() );
      integrateCurve2
	( theta0, theta1,
	  x, y, curv, normal_x, normal_y,
	  args.getOption( "-integrate_truecurve2" )->getFloatValue( 0 ),
	  args.getOption( "-integrate_truecurve2" )->getFloatValue( 1 ),
	  args.getOption( "-integrate_truecurve2" )->getFloatValue( 2 )
	  );
    }

  if ( args.check( "-time" ) )
    cout << my_time << " ms spent on " << nb_trials << " experiment(s)" << endl; 

  if ( variable != 0 ) delete variable;
  if ( tgt_computer != 0 ) delete tgt_computer;
  if ( blurrer != 0 ) delete blurrer;
  
  return 0;
  
}
