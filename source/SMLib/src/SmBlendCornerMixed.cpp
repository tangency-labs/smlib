// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

/**********************************************************************//**
* FILE NAME --- SmBlendCornerMixed.cpp
* PURPOSE: SolidCraft change. The blend corner's Mixed event
*   (SmBlendCorner.h): a convex blend B1 and a concave blend B2 whose edges
*   continue each other G1 at the vertex V, on a common face T (a plane), the
*   third edge E at V, not blended, a straight tangent seam between their
*   other faces F1 and F2 (the test peg's front edge and its concave arc,
*   where the front face runs tangent into the cylinder).
*
* On T the blends' rails cross at X. A convex and a concave blend meeting on
* a face meet at their rails' crossing with no spring curve between them:
* a patch G1 to both along their end sections there cannot be C2 (their
* normal curvatures across, +1/r and -1/r, would ask its shape operator to
* be unsymmetric), and it dips below the faces whatever its form. Each blend
* is instead trimmed along the other's rail continued past X: B1 by the
* surface S2 swept by B2's rail on T along T's normal, from X to A on its
* rail on F1, and B2 by S1, swept by B1's rail on T, from X to B on its rail
* on F2. They leave X along the other's rail, and the corner's shape
* operator there is then symmetric: B1's curvature across its trimline acts
* along T's direction square to B1's axis, B2's along the radius, and the
* trimlines run square to those.
*
* On the wall F1 + F2 (G1 across E) the corner runs along B1's rail from A
* to P, where it meets E, up E's line past V to Q, where B2's rail on F2
* crosses it, and along B2's rail back to B: every corner of that spring has
* the wall's own shape operator on both sides. The corner splits along the
* arc X-V in T, tangent at V to the edges' common direction: the face L
* (X, A, P, V) below T with B1, the face U (X, V, Q, B) above it with B2.
**********************************************************************/

#include "StdAfx.h"

#ifndef __SMBLENDCORNER_H__
#include <SmBlendCorner.h>
#endif

#ifndef __SMFILLETSOLVER_H__
#include <SmFilletSolver.h>
#endif

#include <SmFilletGeom.h>
#include <SmFilletExecutive.h>
#include <SmSolutionArray.h>
#include <SmBSplineCurve.h>
#include <SmBSplineSurface.h>
#include <SmEdgeuse.h>
#include <SmEdge.h>
#include <SmFace.h>
#include <SmVertex.h>
#include <SmAxis2Placement.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <vector>


namespace
{

// The Mixed corner's parts at its vertex, read off the ring (pAllEUs): the convex and the concave
// blended edgeuses, the seam, the faces T, F1 (B1's other), F2 (B2's other), and the blends' solvers.
struct SmMixedRing
{
  SmEdgeuse      * pConvex = NULL, * pConcave = NULL, * pSeam = NULL;
  SmFace         * pT = NULL, * pF1 = NULL, * pF2 = NULL;
  SmFilletSolver * apFS[2] = { NULL, NULL };   // B1, B2
};

SmStatus sm_MixedRing(const SmVertex * /*cpVertex*/, const SmTArray<SmEdgeuse*> & crSolverEUs, const SmTArray<SmEdgeuse*> & crAllEUs,
                      SmFilletExecutive * pExec, SmMixedRing & rRing, SmBoolean & rbIsMixed)
{
  rbIsMixed = FALSE;
  if ( crSolverEUs.GetSize() != 2 || crAllEUs.GetSize() != 3 )
    { return SM_SUCCESS; }
  SmEdgeuse * pEU0 = crSolverEUs[0];
  SmEdgeuse * pEU1 = crSolverEUs[1];
  NER( pEU0 );
  NER( pEU1 );
  SmBoolean abConcave[2] = { FALSE, FALSE };
  SER( SmFilletSolver::IsConcaveEdge( *pEU0->GetEdge(), abConcave[0] ));
  SER( SmFilletSolver::IsConcaveEdge( *pEU1->GetEdge(), abConcave[1] ));
  if ( abConcave[0] == abConcave[1] )
    { return SM_SUCCESS; }
  rRing.pConvex  = abConcave[0] ? pEU1 : pEU0;
  rRing.pConcave = abConcave[0] ? pEU0 : pEU1;
  for (ULONG ii=0; ii<3; ii++)
    {
      if ( crAllEUs[ii] != pEU0 && crAllEUs[ii] != pEU1 )
        { rRing.pSeam = crAllEUs[ii]; }
    }
  NER( rRing.pSeam );
  rRing.pT = sm_GetFaceByTwoEdgeuses( rRing.pConvex, rRing.pConcave, rRing.pSeam );
  if ( rRing.pT == NULL )
    { return SM_SUCCESS; }
  SmFace * apSeam[2] = { rRing.pSeam->GetFace(), rRing.pSeam->GetRadial() != NULL ? rRing.pSeam->GetRadial()->GetFace() : NULL };
  for (ULONG ii=0; ii<2; ii++)
    {
      SmEdgeuse * pB = ii == 0 ? rRing.pConvex : rRing.pConcave;
      SmFace    * pOther = NULL;
      for (ULONG k=0; k<2; k++)
        {
          if ( apSeam[k] != NULL && apSeam[k] != rRing.pT && sm_GetFaceByTwoEdgeuses( pB, rRing.pSeam, NULL ) == apSeam[k] )
            { pOther = apSeam[k]; }
        }
      if ( pOther == NULL )
        { return SM_SUCCESS; }
      ( ii == 0 ? rRing.pF1 : rRing.pF2 ) = pOther;
      rRing.apFS[ii] = pExec->GetFilletSolverOfEdgeuse( pB );
      if ( rRing.apFS[ii] == NULL )
        { rRing.apFS[ii] = pExec->GetFilletSolverOfEdgeuse( pB->GetMate() ); }
      if ( rRing.apFS[ii] == NULL )
        { return SM_SUCCESS; }
    }
  if ( rRing.pF1 == rRing.pF2 )
    { return SM_SUCCESS; }
  rbIsMixed = TRUE;
  return SM_SUCCESS;
}

// The unit tangent of an edge leaving a vertex.
SmStatus sm_MixedLeaving(const SmEdge * pEdge, const SmVertex * cpVertex, SmVector3d & rDir)
{
  const SmBoolean bStart = ( pEdge->GetStartVertex() == cpVertex ) ? TRUE : FALSE;
  SER( SM_CONST_CAST(SmEdge*, pEdge)->GetEndTangent( bStart, rDir ));
  SER( rDir.Unitize() );
  // into the edge: towards its point a little way in from the vertex
  NER( pEdge->GetCurve() );
  const SmExtent1d sIvl = pEdge->GetInterval();
  SmPoint3d sIn;
  SER( pEdge->GetCurve()->EvaluatePoint( sIvl.Evaluate( bStart ? 1e-3 : 1.0 - 1e-3 ), sIn ));
  if ( rDir.Dot( sIn - cpVertex->GetPoint() ) < 0.0 )
    { rDir = -rDir; }
  return SM_SUCCESS;
}

// The parameter of the point of crCurve nearest crPnt.
SmStatus sm_MixedParam(const SmCurve & crCurve, const SmPoint3d & crPnt, double dTol, double & rdT)
{
  SmSolutionArray sSolutions;
  SER( crCurve.GlobalPointSolve( crCurve.GetNaturalInterval(), SM_SO_MINIMIZE, crPnt, dTol, NULL, NULL, SM_SR_SINGLE, sSolutions ));
  if ( sSolutions.GetSize() < 1 )
    { SER( SM_ERR ); }
  rdT = sSolutions[0].m_vStart[0];
  return SM_SUCCESS;
}

// The crossing of a curve with a surface nearest crNear.
SmStatus sm_MixedCurveSurface(const SmCurve & crCurve, const SmSurface & crSurface, const SmPoint3d & crNear, double dTol,
                              SmPoint3d & rPnt, double & rdT)
{
  SmSolutionArray sSolutions;
  SER( crSurface.GlobalCurveIntersect( crSurface.GetNaturalUVDomain(), crCurve, crCurve.GetNaturalInterval(), dTol, sSolutions ));
  double dBest = SM_BIG_DOUBLE;
  for (ULONG ii=0; ii<sSolutions.GetSize(); ii++)
    {
      SmPoint3d sPnt;
      SER( crCurve.EvaluatePoint( sSolutions[ii].m_vStart.m_adParameters[0], sPnt ));
      if ( sPnt.DistanceBetween( crNear ) < dBest )
        {
          dBest = sPnt.DistanceBetween( crNear );
          rPnt  = sPnt;
          rdT   = sSolutions[ii].m_vStart.m_adParameters[0];
        }
    }
  if ( dBest == SM_BIG_DOUBLE )
    { SER_MSG( SM_ERR, _T("Mixed corner: a rail does not meet the other blend's rail surface") ); }
  return SM_SUCCESS;
}

// The surface swept by a rail on the plane T along T's normal, reaching dReach either side.
SmStatus sm_MixedRailWall(const SmContext & crContext, const SmCurve & crRail, const SmVector3d & crNormal, double dReach,
                          SmBSplineSurface *& rpWall)
{
  rpWall = NULL;
  SmCurve * pCopy = NULL;
  SER( crRail.Copy( crContext, pCopy ));
  NER( pCopy );
  SmObjDelete sClean( pCopy );
  SmBSplineCurve * pB0 = SM_CAST_PTR( SmBSplineCurve, pCopy );
  NER( pB0 );
  // the rail continued past both ends as itself (its end spans' own polynomials: a rational arc
  // continues as its circle), so the other blend's rail end, under this one's end, is inside it
  SmBSplineCurve * pB1 = NULL, * pB = NULL;
  SER( pB0->CreateExtendedCurve( crContext, 0.25 * dReach, 1, SM_CT_CINFINITY, pB1 ));
  NER( pB1 );
  SmObjDelete sClean1( pB1 );
  SER( pB1->CreateExtendedCurve( crContext, 0.25 * dReach, 2, SM_CT_CINFINITY, pB ));
  NER( pB );
  SmObjDelete sClean2( pB );
  const SmVector3d sShift = crNormal * ( -dReach );
  SER( pB->Transform( SmAxis2Placement( SmPoint3d( sShift.x, sShift.y, sShift.z ), SmVector3d( 1, 0, 0 ), SmVector3d( 0, 1, 0 ))));
  SER( SmBSplineSurface::CreateLinearSweep( crContext, *pB, crNormal * ( 2.0 * dReach ), rpWall ));
  NER( rpWall );
  return SM_SUCCESS;
}


// ---------------------------------------------------------------------------------------------
// Exact piecewise-Bezier algebra in homogeneous coordinates (x w, y w, z w, w), for the corner's
// faces. Every operation here is exact (to rounding): knot insertion, de Casteljau subdivision,
// degree elevation, Bernstein products, sums, and the composition of a rational surface with a
// curve in its parameters. A curve is a list of Bezier segments of one degree over its breaks.

typedef std::array<double,4> SmMxV4;

SmMxV4 sm_MxAdd(const SmMxV4 & a, const SmMxV4 & b) { return { a[0] + b[0], a[1] + b[1], a[2] + b[2], a[3] + b[3] }; }
SmMxV4 sm_MxSub(const SmMxV4 & a, const SmMxV4 & b) { return { a[0] - b[0], a[1] - b[1], a[2] - b[2], a[3] - b[3] }; }
SmMxV4 sm_MxMul(const SmMxV4 & a, double s)        { return { a[0] * s, a[1] * s, a[2] * s, a[3] * s }; }
SmMxV4 sm_MxLerp(const SmMxV4 & a, const SmMxV4 & b, double t) { return sm_MxAdd( sm_MxMul( a, 1.0 - t ), sm_MxMul( b, t )); }
SmMxV4 sm_MxPoint(const SmPoint3d & p, double w = 1.0) { return { p.x * w, p.y * w, p.z * w, w }; }
SmMxV4 sm_MxDir(const SmVector3d & d) { return { d.x, d.y, d.z, 0.0 }; }
SmPoint3d sm_MxProject(const SmMxV4 & a) { return SmPoint3d( a[0] / a[3], a[1] / a[3], a[2] / a[3] ); }

double sm_MxChoose(int n, int k)
{
  double d = 1.0;
  for (int i=1; i<=k; i++) { d = d * double( n - k + i ) / double( i ); }
  return d;
}

typedef std::vector<SmMxV4> SmMxBez; // degree = size - 1, on [0, 1]

struct SmMxCrv
{
  std::vector<double>  brk; // segment s on [brk[s], brk[s+1]]
  std::vector<SmMxBez> seg;
  int Deg() const { return seg.empty() ? -1 : (int)seg[0].size() - 1; }
};

SmMxV4 sm_MxBezEval(const SmMxBez & b, double t)
{
  SmMxBez c = b;
  for (size_t r=1; r<c.size(); r++)
    for (size_t i=0; i+r<c.size(); i++)
      c[i] = sm_MxLerp( c[i], c[i+1], t );
  return c[0];
}

void sm_MxBezSplit(const SmMxBez & b, double t, SmMxBez & rL, SmMxBez & rR)
{
  const size_t n = b.size();
  SmMxBez c = b;
  rL.assign( n, SmMxV4{} );
  rR.assign( n, SmMxV4{} );
  rL[0] = c[0];
  rR[n-1] = c[n-1];
  for (size_t r=1; r<n; r++)
    {
      for (size_t i=0; i+r<n; i++)
        c[i] = sm_MxLerp( c[i], c[i+1], t );
      rL[r] = c[0];
      rR[n-1-r] = c[n-1-r];
    }
}

SmMxBez sm_MxBezElevate(const SmMxBez & b, int lTo)
{
  SmMxBez c = b;
  while ( (int)c.size() - 1 < lTo )
    {
      const int n = (int)c.size() - 1;
      SmMxBez e( n + 2 );
      e[0] = c[0];
      e[n+1] = c[n];
      for (int i=1; i<=n; i++)
        e[i] = sm_MxAdd( sm_MxMul( c[i-1], double( i ) / double( n + 1 )), sm_MxMul( c[i], 1.0 - double( i ) / double( n + 1 )));
      c = e;
    }
  return c;
}

// d/dt of a segment over [a, b] (a Bezier of one degree less)
SmMxBez sm_MxBezDeriv(const SmMxBez & b, double dLen)
{
  const int n = (int)b.size() - 1;
  if ( n == 0 )
    return SmMxBez( 1, SmMxV4{} );
  SmMxBez d( n );
  for (int i=0; i<n; i++)
    d[i] = sm_MxMul( sm_MxSub( b[i+1], b[i] ), double( n ) / dLen );
  return d;
}

// the product of a scalar Bezier (its first coordinates) and a Bezier
SmMxBez sm_MxBezTimes(const SmMxBez & s, const SmMxBez & p)
{
  const int m = (int)s.size() - 1, n = (int)p.size() - 1;
  SmMxBez r( m + n + 1, SmMxV4{} );
  for (int i=0; i<=m; i++)
    for (int j=0; j<=n; j++)
      r[i+j] = sm_MxAdd( r[i+j], sm_MxMul( p[j], s[i][0] * sm_MxChoose( m, i ) * sm_MxChoose( n, j )));
  for (int k=0; k<=m+n; k++)
    r[k] = sm_MxMul( r[k], 1.0 / sm_MxChoose( m + n, k ));
  return r;
}

SmMxV4 sm_MxEval(const SmMxCrv & c, double t)
{
  size_t s = 0;
  while ( s + 1 < c.seg.size() && t > c.brk[s+1] ) s++;
  const double a = c.brk[s], b = c.brk[s+1];
  return sm_MxBezEval( c.seg[s], ( t - a ) / ( b - a ));
}

SmMxCrv sm_MxElevate(const SmMxCrv & c, int lTo)
{
  SmMxCrv r = c;
  for (size_t s=0; s<r.seg.size(); s++) r.seg[s] = sm_MxBezElevate( r.seg[s], lTo );
  return r;
}

// c's segments split at every break of crBrk inside them
SmMxCrv sm_MxRefine(const SmMxCrv & c, const std::vector<double> & crBrk)
{
  SmMxCrv r;
  r.brk.push_back( c.brk.front() );
  for (size_t s=0; s<c.seg.size(); s++)
    {
      const double a = c.brk[s], b = c.brk[s+1];
      SmMxBez rest = c.seg[s];
      double from = a;
      for (double t : crBrk)
        {
          if ( t <= from + 1e-14 * ( 1.0 + fabs( t )) || t >= b - 1e-14 * ( 1.0 + fabs( t )))
            continue;
          SmMxBez L, R;
          sm_MxBezSplit( rest, ( t - from ) / ( b - from ), L, R );
          r.seg.push_back( L );
          r.brk.push_back( t );
          rest = R;
          from = t;
        }
      r.seg.push_back( rest );
      r.brk.push_back( b );
    }
  return r;
}

std::vector<double> sm_MxUnion(const std::vector<double> & a, const std::vector<double> & b)
{
  std::vector<double> u = a;
  u.insert( u.end(), b.begin(), b.end() );
  std::sort( u.begin(), u.end() );
  std::vector<double> r;
  for (double t : u)
    if ( r.empty() || t > r.back() + 1e-13 * ( 1.0 + fabs( t )))
      r.push_back( t );
  return r;
}

// a and b on common breaks and degree
void sm_MxCommon(SmMxCrv & a, SmMxCrv & b)
{
  const std::vector<double> u = sm_MxUnion( a.brk, b.brk );
  a = sm_MxRefine( a, u );
  b = sm_MxRefine( b, u );
  const int d = std::max( a.Deg(), b.Deg() );
  a = sm_MxElevate( a, d );
  b = sm_MxElevate( b, d );
  // the breaks as one list (rounding apart only)
  b.brk = a.brk;
}

SmMxCrv sm_MxSum(SmMxCrv a, SmMxCrv b, double sa = 1.0, double sb = 1.0)
{
  sm_MxCommon( a, b );
  for (size_t s=0; s<a.seg.size(); s++)
    for (size_t i=0; i<a.seg[s].size(); i++)
      a.seg[s][i] = sm_MxAdd( sm_MxMul( a.seg[s][i], sa ), sm_MxMul( b.seg[s][i], sb ));
  return a;
}

// s (a scalar curve) times p
SmMxCrv sm_MxTimes(SmMxCrv s, SmMxCrv p)
{
  const std::vector<double> u = sm_MxUnion( s.brk, p.brk );
  s = sm_MxRefine( s, u );
  p = sm_MxRefine( p, u );
  SmMxCrv r;
  r.brk = s.brk;
  for (size_t k=0; k<s.seg.size(); k++) r.seg.push_back( sm_MxBezTimes( s.seg[k], p.seg[k] ));
  return r;
}

SmMxCrv sm_MxDeriv(const SmMxCrv & c)
{
  SmMxCrv r;
  r.brk = c.brk;
  for (size_t s=0; s<c.seg.size(); s++) r.seg.push_back( sm_MxBezDeriv( c.seg[s], c.brk[s+1] - c.brk[s] ));
  return r;
}

// c with every control point mapped by f (a linear map of R^4)
template <class F> SmMxCrv sm_MxMap(const SmMxCrv & c, F f)
{
  SmMxCrv r = c;
  for (auto & b : r.seg) for (auto & p : b) p = f( p );
  return r;
}

// the scalar cubic Hermite on [0, 1]: values v0, v1, derivatives d0, d1
SmMxCrv sm_MxHermite(double v0, double v1, double d0, double d1)
{
  SmMxCrv r;
  r.brk = { 0.0, 1.0 };
  r.seg.push_back( SmMxBez{ SmMxV4{ v0, 0, 0, 0 }, SmMxV4{ v0 + d0 / 3.0, 0, 0, 0 }, SmMxV4{ v1 - d1 / 3.0, 0, 0, 0 }, SmMxV4{ v1, 0, 0, 0 } } );
  return r;
}

SmMxCrv sm_MxConst(const SmMxV4 & v)
{
  SmMxCrv r;
  r.brk = { 0.0, 1.0 };
  r.seg.push_back( SmMxBez{ v } );
  return r;
}

// c on [0, 1] (linearly), or reversed
SmMxCrv sm_MxUnit(const SmMxCrv & c, SmBoolean bReverse)
{
  SmMxCrv r;
  const double a = c.brk.front(), b = c.brk.back();
  const size_t n = c.seg.size();
  for (size_t i=0; i<=n; i++)
    {
      const double t = ( c.brk[bReverse ? n - i : i] - a ) / ( b - a );
      r.brk.push_back( bReverse ? 1.0 - t : t );
    }
  for (size_t i=0; i<n; i++)
    {
      SmMxBez s = c.seg[bReverse ? n - 1 - i : i];
      if ( bReverse ) std::reverse( s.begin(), s.end() );
      r.seg.push_back( s );
    }
  r.brk.front() = 0.0;
  r.brk.back()  = 1.0;
  return r;
}

// c scaled by the cubic factor rho that makes its weight 1 with zero derivative at both ends (the
// same curve, the same parameterization): its corners' homogeneous tangents then have no weight
// part, and a face's weights stay near 1. A polynomial curve is left as it is; a rho not positive
// on the whole interval (seen on its Bernstein coefficients) is taken linear, its end weights 1.
SmMxCrv sm_MxUnitWeights(const SmMxCrv & c)
{
  const double a = c.brk.front(), b = c.brk.back(), L = b - a;
  SmMxCrv d = c;
  for (size_t s=0; s<d.seg.size(); s++) d.seg[s] = sm_MxBezDeriv( c.seg[s], c.brk[s+1] - c.brk[s] );
  const double w0 = sm_MxEval( c, a )[3], w1 = sm_MxEval( c, b )[3];
  const double dw0 = sm_MxEval( d, a )[3], dw1 = sm_MxEval( d, b )[3];
  if ( fabs( w0 - 1.0 ) < 1e-15 && fabs( w1 - 1.0 ) < 1e-15 && fabs( dw0 ) < 1e-15 && fabs( dw1 ) < 1e-15 )
    return c;
  // rho on [a, b] as a cubic: values 1/w, derivatives -w'/w^2 (per unit of t)
  const double r0 = 1.0 / w0, r1 = 1.0 / w1, s0 = -dw0 / ( w0 * w0 ) * L, s1 = -dw1 / ( w1 * w1 ) * L;
  const double q[4] = { r0, r0 + s0 / 3.0, r1 - s1 / 3.0, r1 };
  SmMxCrv l;
  l.brk = { a, b };
  // (a curve of degree above 10 takes the linear factor: its ribbon, three degrees above it with the
  // cubic, keeps within the B-spline degree limit, 17)
  if ( q[1] > 0.0 && q[2] > 0.0 && c.Deg() <= 10 )
    l.seg.push_back( SmMxBez{ SmMxV4{ q[0], 0, 0, 0 }, SmMxV4{ q[1], 0, 0, 0 }, SmMxV4{ q[2], 0, 0, 0 }, SmMxV4{ q[3], 0, 0, 0 } } );
  else
    l.seg.push_back( SmMxBez{ SmMxV4{ r0, 0, 0, 0 }, SmMxV4{ r1, 0, 0, 0 } } );
  return sm_MxTimes( l, c );
}

// the B-spline curve's segments by knot insertion (homogeneous; a 2D curve's z is 0)
SmStatus sm_MxFromBSpline(const SmBSplineCurve & crCurve, SmMxCrv & rOut)
{
  ULONG lDim = 0, lDeg = 0;
  SmTArray<SmPoint3d> sPts;
  SmBSplineCurveForm eForm;
  SmTArray<ULONG> sMults;
  SmTArray<double> sKnots, sW;
  SmKnotType eKT;
  SER( crCurve.GetCanonical( lDim, lDeg, sPts, eForm, sMults, sKnots, eKT, sW ));
  const int p = (int)lDeg;
  std::vector<double> U;
  for (ULONG i=0; i<sKnots.GetSize(); i++)
    for (ULONG k=0; k<sMults[i]; k++) U.push_back( sKnots[i] );
  std::vector<SmMxV4> P;
  for (ULONG i=0; i<sPts.GetSize(); i++)
    P.push_back( sm_MxPoint( sPts[i], sW.GetSize() > 0 ? sW[i] : 1.0 ));
  if ( U.size() != P.size() + p + 1 )
    { SER( SM_ERR ); }
  // insert every interior knot up to multiplicity p (Boehm)
  for (ULONG ki=1; ki+1<sKnots.GetSize(); ki++)
    {
      const double t = sKnots[ki];
      for (ULONG m=sMults[ki]; m<(ULONG)p; m++)
        {
          size_t k = 0;
          while ( k + 1 < U.size() && U[k+1] <= t ) k++;
          // t in [U[k], U[k+1]), s = current multiplicity
          int s = 0;
          for (size_t q=0; q<U.size(); q++) if ( U[q] == t ) s++;
          std::vector<SmMxV4> Q( P.size() + 1 );
          for (size_t i=0; i<=k-p; i++) Q[i] = P[i];
          for (size_t i=k-s+1; i<=P.size(); i++) Q[i] = P[i-1];
          for (size_t i=k-p+1; i<=k-s; i++)
            {
              const double a = ( t - U[i] ) / ( U[i+p] - U[i] );
              Q[i] = sm_MxLerp( P[i-1], P[i], a );
            }
          U.insert( U.begin() + k + 1, t );
          P = Q;
        }
    }
  rOut.brk.clear();
  rOut.seg.clear();
  for (ULONG i=0; i<sKnots.GetSize(); i++) rOut.brk.push_back( sKnots[i] );
  for (size_t s=0; s+1<rOut.brk.size(); s++)
    rOut.seg.push_back( SmMxBez( P.begin() + s * p, P.begin() + s * p + p + 1 ));
  return SM_SUCCESS;
}

// A surface's Bezier patches: patch (i, j) over [ub[i], ub[i+1]] x [vb[j], vb[j+1]], control
// point (a, b) at a * (dv + 1) + b.
struct SmMxSrf
{
  int du = 0, dv = 0;
  std::vector<double> ub, vb;
  std::vector<std::vector<SmMxV4>> P; // index i * (vb.size() - 1) + j
  const std::vector<SmMxV4> & At(size_t i, size_t j) const { return P[i * ( vb.size() - 1 ) + j]; }
};

SmStatus sm_MxSrfFromBSpline(const SmBSplineSurface & crSurf, SmMxSrf & rOut)
{
  ULONG lDu = 0, lDv = 0;
  SmTArray<SmPoint3d> sPts;
  SmBSplineSurfaceForm eForm;
  SmTArray<ULONG> sUM, sVM;
  SmTArray<double> sUK, sVK, sW;
  SmKnotType eKT;
  SER( crSurf.GetCanonical( lDu, lDv, sPts, eForm, sUM, sVM, sUK, sVK, eKT, sW ));
  ULONG lNu = 0, lNv = 0;
  for (ULONG i=0; i<sUM.GetSize(); i++) lNu += sUM[i];
  for (ULONG i=0; i<sVM.GetSize(); i++) lNv += sVM[i];
  lNu -= lDu + 1;
  lNv -= lDv + 1;
  if ( lNu * lNv != sPts.GetSize() )
    { SER( SM_ERR ); }
  // each row (u fixed index) a curve in v, made piecewise Bezier; then each column in u
  auto toCrv = [&](const std::vector<SmMxV4> & crP, ULONG lDeg, const SmTArray<ULONG> & crM, const SmTArray<double> & crK, SmMxCrv & rC) -> SmStatus
    {
      // build a temporary B-spline in 4D via the Boehm path of sm_MxFromBSpline, on raw data
      const int p = (int)lDeg;
      std::vector<double> U;
      for (ULONG i=0; i<crK.GetSize(); i++)
        for (ULONG k=0; k<crM[i]; k++) U.push_back( crK[i] );
      std::vector<SmMxV4> P = crP;
      for (ULONG ki=1; ki+1<crK.GetSize(); ki++)
        {
          const double t = crK[ki];
          for (ULONG m=crM[ki]; m<(ULONG)p; m++)
            {
              size_t k = 0;
              while ( k + 1 < U.size() && U[k+1] <= t ) k++;
              int s = 0;
              for (size_t q=0; q<U.size(); q++) if ( U[q] == t ) s++;
              std::vector<SmMxV4> Q( P.size() + 1 );
              for (size_t i=0; i<=k-p; i++) Q[i] = P[i];
              for (size_t i=k-s+1; i<=P.size(); i++) Q[i] = P[i-1];
              for (size_t i=k-p+1; i<=k-s; i++)
                Q[i] = sm_MxLerp( P[i-1], P[i], ( t - U[i] ) / ( U[i+p] - U[i] ));
              U.insert( U.begin() + k + 1, t );
              P = Q;
            }
        }
      rC.brk.clear();
      rC.seg.clear();
      for (ULONG i=0; i<crK.GetSize(); i++) rC.brk.push_back( crK[i] );
      for (size_t s=0; s+1<rC.brk.size(); s++)
        rC.seg.push_back( SmMxBez( P.begin() + s * p, P.begin() + s * p + p + 1 ));
      return SM_SUCCESS;
    };
  std::vector<SmMxCrv> sRows( lNu );
  for (ULONG i=0; i<lNu; i++)
    {
      std::vector<SmMxV4> row;
      for (ULONG j=0; j<lNv; j++)
        row.push_back( sm_MxPoint( sPts[i * lNv + j], sW.GetSize() > 0 ? sW[i * lNv + j] : 1.0 ));
      SER( toCrv( row, lDv, sVM, sVK, sRows[i] ));
    }
  const size_t nvs = sRows[0].seg.size();
  rOut.du = (int)lDu;
  rOut.dv = (int)lDv;
  rOut.vb = sRows[0].brk;
  rOut.ub.clear();
  for (ULONG i=0; i<sUK.GetSize(); i++) rOut.ub.push_back( sUK[i] );
  const size_t nus = rOut.ub.size() - 1;
  rOut.P.assign( nus * nvs, std::vector<SmMxV4>( ( lDu + 1 ) * ( lDv + 1 )));
  for (size_t js=0; js<nvs; js++)
    for (ULONG b=0; b<=lDv; b++)
      {
        std::vector<SmMxV4> col;
        for (ULONG i=0; i<lNu; i++) col.push_back( sRows[i].seg[js][b] );
        SmMxCrv c;
        SER( toCrv( col, lDu, sUM, sUK, c ));
        for (size_t is=0; is<nus; is++)
          for (ULONG a=0; a<=lDu; a++)
            rOut.P[is * nvs + js][a * ( lDv + 1 ) + b] = c.seg[is][a];
      }
  return SM_SUCCESS;
}

// the parameters in (0, 1) where a scalar Bezier crosses dLevel (subdivision, then bisection to
// rounding)
void sm_MxRoots(const std::vector<double> & crB, double dLevel, double a, double b, int depth, std::vector<double> & rRoots)
{
  double lo = SM_BIG_DOUBLE, hi = -SM_BIG_DOUBLE;
  for (double c : crB) { lo = std::min( lo, c - dLevel ); hi = std::max( hi, c - dLevel ); }
  if ( lo > 0.0 || hi < 0.0 )
    return;
  if ( depth > 60 || b - a < 1e-15 )
    {
      rRoots.push_back( 0.5 * ( a + b ));
      return;
    }
  const double fa = crB.front() - dLevel, fb = crB.back() - dLevel;
  // monotone control polygon with a sign change: bisect on the Bezier
  SmBoolean bMono = TRUE;
  for (size_t i=1; i<crB.size() && bMono; i++)
    if ( ( crB[i] - crB[i-1] ) * ( fb - fa ) < 0.0 ) bMono = FALSE;
  if ( bMono && fa * fb < 0.0 )
    {
      SmMxBez z;
      for (double c : crB) z.push_back( SmMxV4{ c - dLevel, 0, 0, 0 } );
      double l = 0.0, h = 1.0;
      for (int it=0; it<200; it++)
        {
          const double m = 0.5 * ( l + h );
          if ( ( sm_MxBezEval( z, m )[0] < 0.0 ) == ( fa < 0.0 )) l = m; else h = m;
          if ( h - l < 1e-17 ) break;
        }
      rRoots.push_back( a + ( b - a ) * 0.5 * ( l + h ));
      return;
    }
  if ( fa == 0.0 || fb == 0.0 )
    return; // ends at a break: no split inside
  SmMxBez z, L, R;
  for (double c : crB) z.push_back( SmMxV4{ c, 0, 0, 0 } );
  sm_MxBezSplit( z, 0.5, L, R );
  std::vector<double> l, r;
  for (auto & q : L) l.push_back( q[0] );
  for (auto & q : R) r.push_back( q[0] );
  sm_MxRoots( l, dLevel, a, 0.5 * ( a + b ), depth + 1, rRoots );
  sm_MxRoots( r, dLevel, 0.5 * ( a + b ), b, depth + 1, rRoots );
}

// S(pc(t)): the surface's homogeneous form along a polynomial curve in its parameters (x, y),
// exact: pc is split where it crosses the patches' breaks, and on each piece the patch's
// Bernstein basis in (s(t), t(t)) is expanded by products, of degree (du + dv) times pc's.
SmStatus sm_MxCompose(const SmMxSrf & S, const SmMxCrv & crPc, SmMxCrv & rOut)
{
  // the breaks of pc where it crosses a patch boundary
  std::vector<double> sCut;
  for (size_t s=0; s<crPc.seg.size(); s++)
    {
      const double a = crPc.brk[s], b = crPc.brk[s+1];
      for (int k=0; k<2; k++)
        {
          std::vector<double> c;
          for (auto & q : crPc.seg[s]) c.push_back( q[k] / q[3] );
          const std::vector<double> & lv = k == 0 ? S.ub : S.vb;
          for (size_t i=1; i+1<lv.size(); i++)
            {
              std::vector<double> roots;
              sm_MxRoots( c, lv[i], 0.0, 1.0, 0, roots );
              for (double r : roots) sCut.push_back( a + ( b - a ) * r );
            }
        }
    }
  SmMxCrv pc = sm_MxRefine( crPc, sCut );
  rOut.brk = pc.brk;
  rOut.seg.clear();
  for (size_t s=0; s<pc.seg.size(); s++)
    {
      const SmMxBez & q = pc.seg[s];
      const SmMxV4 m = sm_MxBezEval( q, 0.5 );
      const double x = m[0] / m[3], y = m[1] / m[3];
      size_t i = 0, j = 0;
      while ( i + 2 < S.ub.size() && x > S.ub[i+1] ) i++;
      while ( j + 2 < S.vb.size() && y > S.vb[j+1] ) j++;
      const double u0 = S.ub[i], du = S.ub[i+1] - S.ub[i], v0 = S.vb[j], dv = S.vb[j+1] - S.vb[j];
      // homogeneous local coordinates: s = X / W, t = Y / W, and W
      SmMxBez X, Y, W, XW, YW; // scalars in [0]
      for (auto & c : q)
        {
          X.push_back( SmMxV4{ ( c[0] - u0 * c[3] ) / du, 0, 0, 0 } );
          Y.push_back( SmMxV4{ ( c[1] - v0 * c[3] ) / dv, 0, 0, 0 } );
          W.push_back( SmMxV4{ c[3], 0, 0, 0 } );
          XW.push_back( SmMxV4{ c[3] - ( c[0] - u0 * c[3] ) / du, 0, 0, 0 } );
          YW.push_back( SmMxV4{ c[3] - ( c[1] - v0 * c[3] ) / dv, 0, 0, 0 } );
        }
      auto power = [](const SmMxBez & b, int k) { SmMxBez r{ SmMxV4{ 1, 0, 0, 0 } }; for (int iPow=0; iPow<k; iPow++) r = sm_MxBezTimes( b, r ); return r; };
      std::vector<SmMxBez> Bu( S.du + 1 ), Bv( S.dv + 1 );
      for (int a=0; a<=S.du; a++)
        {
          SmMxBez t = sm_MxBezTimes( power( X, a ), power( XW, S.du - a ));
          for (auto & c : t) c = sm_MxMul( c, sm_MxChoose( S.du, a ));
          Bu[a] = t;
        }
      for (int b=0; b<=S.dv; b++)
        {
          SmMxBez t = sm_MxBezTimes( power( Y, b ), power( YW, S.dv - b ));
          for (auto & c : t) c = sm_MxMul( c, sm_MxChoose( S.dv, b ));
          Bv[b] = t;
        }
      const std::vector<SmMxV4> & P = S.At( i, j );
      SmMxBez R;
      for (int a=0; a<=S.du; a++)
        for (int b=0; b<=S.dv; b++)
          {
            SmMxBez basis = sm_MxBezTimes( Bu[a], Bv[b] );
            SmMxBez term;
            for (auto & c : basis) term.push_back( sm_MxMul( P[a * ( S.dv + 1 ) + b], c[0] ));
            if ( R.empty() ) R = term; else for (size_t k=0; k<R.size(); k++) R[k] = sm_MxAdd( R[k], term[k] );
          }
      rOut.seg.push_back( R );
    }
  return SM_SUCCESS;
}

// the least-norm solution of A x = b (rows of A given), by the eigen decomposition of A^T A
// (Jacobi); rdResidual |A x - b|
void sm_MxSolve(const std::vector<std::vector<double>> & A, const std::vector<double> & b, std::vector<double> & x, double & rdResidual)
{
  const size_t m = A.size(), n = A.empty() ? 0 : A[0].size();
  std::vector<std::vector<double>> N( n, std::vector<double>( n, 0.0 )), V( n, std::vector<double>( n, 0.0 ));
  std::vector<double> r( n, 0.0 );
  for (size_t i=0; i<n; i++)
    {
      V[i][i] = 1.0;
      for (size_t j=0; j<n; j++)
        for (size_t k=0; k<m; k++) N[i][j] += A[k][i] * A[k][j];
      for (size_t k=0; k<m; k++) r[i] += A[k][i] * b[k];
    }
  for (int sweep=0; sweep<100; sweep++)
    {
      double off = 0.0;
      for (size_t p=0; p<n; p++) for (size_t q=p+1; q<n; q++) off += N[p][q] * N[p][q];
      if ( off < 1e-30 ) break;
      for (size_t p=0; p<n; p++)
        for (size_t q=p+1; q<n; q++)
          {
            if ( fabs( N[p][q] ) < 1e-300 ) continue;
            const double th = 0.5 * atan2( 2.0 * N[p][q], N[q][q] - N[p][p] );
            const double c = cos( th ), s = sin( th );
            for (size_t k=0; k<n; k++)
              {
                const double a = N[k][p], d = N[k][q];
                N[k][p] = c * a - s * d;
                N[k][q] = s * a + c * d;
              }
            for (size_t k=0; k<n; k++)
              {
                const double a = N[p][k], d = N[q][k];
                N[p][k] = c * a - s * d;
                N[q][k] = s * a + c * d;
              }
            for (size_t k=0; k<n; k++)
              {
                const double a = V[k][p], d = V[k][q];
                V[k][p] = c * a - s * d;
                V[k][q] = s * a + c * d;
              }
          }
    }
  double emax = 0.0;
  for (size_t i=0; i<n; i++) emax = std::max( emax, fabs( N[i][i] ));
  x.assign( n, 0.0 );
  for (size_t i=0; i<n; i++)
    {
      if ( N[i][i] <= 1e-14 * emax ) continue;
      double c = 0.0;
      for (size_t k=0; k<n; k++) c += V[k][i] * r[k];
      c /= N[i][i];
      for (size_t k=0; k<n; k++) x[k] += V[k][i] * c;
    }
  rdResidual = 0.0;
  for (size_t k=0; k<m; k++)
    {
      double e = -b[k];
      for (size_t i=0; i<n; i++) e += A[k][i] * x[i];
      rdResidual += e * e;
    }
  rdResidual = sqrt( rdResidual );
}

// ---------------------------------------------------------------------------------------------
// The corner's faces: four-sided patches, each the bicubically blended Coons patch (in R^4) of
// its four sides and their cross-derivative ribbons. A side is a curve c on [0, 1] with end
// weights 1, and its ribbon g = alpha c + beta c' + gamma D, D a direction field whose span with
// the curve's point and tangent projects onto the neighbour's tangent plane along it: G1 to the
// neighbour exactly, whatever the cubic Hermite scalars alpha, beta, gamma. Their end values come
// from the adjacent sides' tangents at the corners, their end derivatives from the corners'
// twists, solved (least-norm) where both ribbons meet; each solve's residual is checked.

struct SmMxSide
{
  SmMxCrv c;                 // the side, on [0, 1], end weights 1
  SmMxCrv D;                 // its transversal direction field
  SmBoolean bShared = FALSE; // D is the field M = A + k B shared along X-V (k solved with the twists)
  double v[2][3] = {};       // alpha, beta, gamma at u = 0, 1
  double d[2][3] = {};       // their derivatives
  SmMxV4 Dv[2], Dd[2], Dk[2]; // D and D' at the ends (Dk: D''s part per unit of k)
};

// the side's D at its ends (an unshared side's own)
void sm_MxSideEnds(SmMxSide & rS)
{
  for (int e=0; e<2; e++)
    {
      rS.Dv[e] = sm_MxEval( rS.D, e );
      SmMxCrv d = rS.D;
      for (size_t s=0; s<d.seg.size(); s++) d.seg[s] = sm_MxBezDeriv( rS.D.seg[s], rS.D.brk[s+1] - rS.D.brk[s] );
      rS.Dd[e] = sm_MxEval( d, e );
      rS.Dk[e] = SmMxV4{};
    }
}

SmMxV4 sm_MxCrvDerivAt(const SmMxCrv & c, double t, int k)
{
  SmMxCrv r = c;
  for (int i=0; i<k; i++) r = sm_MxDeriv( r );
  return sm_MxEval( r, t );
}

// alpha c + beta c' + gamma D at the side's end e = target: the first-order corner condition
SmStatus sm_MxFirstOrder(SmMxSide & rS, int e, const SmMxV4 & crTarget, double dScale)
{
  const double t = e;
  const SmMxV4 c = sm_MxEval( rS.c, t ), c1 = sm_MxCrvDerivAt( rS.c, t, 1 ), D = rS.Dv[e];
  std::vector<std::vector<double>> A( 4, std::vector<double>( 3 ));
  std::vector<double> b( 4 ), x;
  for (int k=0; k<4; k++) { A[k][0] = c[k]; A[k][1] = c1[k]; A[k][2] = D[k]; b[k] = crTarget[k]; }
  double dRes = 0.0;
  sm_MxSolve( A, b, x, dRes );
  if ( dRes > 1e-9 * dScale )
    { SER_MSG( SM_ERR, _T("Mixed corner: a side's tangent leaves its neighbour's tangent plane at a corner") ); }
  for (int k=0; k<3; k++) rS.v[e][k] = x[k];
  return SM_SUCCESS;
}

// One ribbon end in a twist equation: its side, which end, and the sign it enters with.
struct SmMxTwistTerm
{
  SmMxSide * pS;
  int        e;
  double     sign;
};

// the twists at corners: sum over each equation's two terms of sign * g'(e) = 0, the unknowns each
// term's (alpha', beta', gamma') and, where its D is M = A + k B, the scalar k (shared)
SmStatus sm_MxTwists(std::vector<std::array<SmMxTwistTerm,2>> & rEqs, double dScale, double * pK)
{
  // unknowns: 3 per term, then k when asked
  std::vector<SmMxTwistTerm*> sTerms;
  for (auto & q : rEqs) { sTerms.push_back( &q[0] ); sTerms.push_back( &q[1] ); }
  const size_t n = 3 * sTerms.size() + ( pK != NULL ? 1 : 0 );
  std::vector<std::vector<double>> A;
  std::vector<double> b;
  for (size_t qi=0; qi<rEqs.size(); qi++)
    {
      std::vector<std::vector<double>> rows( 4, std::vector<double>( n, 0.0 ));
      std::vector<double> rhs( 4, 0.0 );
      for (int h=0; h<2; h++)
        {
          SmMxTwistTerm & T = rEqs[qi][h];
          const size_t idx = 3 * ( 2 * qi + h );
          const double t = T.e;
          const SmMxV4 c = sm_MxEval( T.pS->c, t ), c1 = sm_MxCrvDerivAt( T.pS->c, t, 1 ), c2 = sm_MxCrvDerivAt( T.pS->c, t, 2 );
          const SmMxV4 D = T.pS->Dv[T.e];
          const double al = T.pS->v[T.e][0], be = T.pS->v[T.e][1], ga = T.pS->v[T.e][2];
          SmMxV4 known = sm_MxAdd( sm_MxAdd( sm_MxMul( c1, al ), sm_MxMul( c2, be )), sm_MxMul( T.pS->Dd[T.e], ga ));
          for (int k=0; k<4; k++)
            {
              rows[k][idx+0] += T.sign * c[k];
              rows[k][idx+1] += T.sign * c1[k];
              rows[k][idx+2] += T.sign * D[k];
              if ( pK != NULL )
                rows[k][n-1] += T.sign * ga * T.pS->Dk[T.e][k];
              rhs[k] -= T.sign * known[k];
            }
        }
      for (int k=0; k<4; k++) { A.push_back( rows[k] ); b.push_back( rhs[k] ); }
    }
  std::vector<double> x;
  double dRes = 0.0;
  sm_MxSolve( A, b, x, dRes );
  if ( dRes > 1e-8 * dScale )
    { SER_MSG( SM_ERR, _T("Mixed corner: the corner's twists are not compatible") ); }
  for (size_t i=0; i<sTerms.size(); i++)
    for (int k=0; k<3; k++) sTerms[i]->pS->d[sTerms[i]->e][k] = x[3*i+k];
  if ( pK != NULL )
    *pK = x[n-1];
  return SM_SUCCESS;
}

SmMxCrv sm_MxRibbon(const SmMxSide & crS)
{
  const SmMxCrv a = sm_MxHermite( crS.v[0][0], crS.v[1][0], crS.d[0][0], crS.d[1][0] );
  const SmMxCrv b = sm_MxHermite( crS.v[0][1], crS.v[1][1], crS.d[0][1], crS.d[1][1] );
  const SmMxCrv g = sm_MxHermite( crS.v[0][2], crS.v[1][2], crS.d[0][2], crS.d[1][2] );
  SmMxCrv r = sm_MxTimes( a, crS.c );
  r = sm_MxSum( r, sm_MxTimes( b, sm_MxDeriv( crS.c )));
  r = sm_MxSum( r, sm_MxTimes( g, crS.D ));
  return r;
}

// surface refinement and elevation, direction by direction (each row/column a curve)
SmMxSrf sm_MxSrfMapU(const SmMxSrf & S, std::function<SmMxCrv(const SmMxCrv &)> f)
{
  const size_t nus = S.ub.size() - 1, nvs = S.vb.size() - 1;
  SmMxSrf R;
  R.dv = S.dv;
  R.vb = S.vb;
  for (size_t j=0; j<nvs; j++)
    for (int b=0; b<=S.dv; b++)
      {
        SmMxCrv c;
        c.brk = S.ub;
        for (size_t i=0; i<nus; i++)
          {
            SmMxBez z;
            for (int a=0; a<=S.du; a++) z.push_back( S.At( i, j )[a * ( S.dv + 1 ) + b] );
            c.seg.push_back( z );
          }
        const SmMxCrv m = f( c );
        if ( R.P.empty() )
          {
            R.du = m.Deg();
            R.ub = m.brk;
            R.P.assign( ( m.brk.size() - 1 ) * nvs, std::vector<SmMxV4>( ( R.du + 1 ) * ( R.dv + 1 )));
          }
        for (size_t i=0; i<m.seg.size(); i++)
          for (int a=0; a<=R.du; a++)
            R.P[i * nvs + j][a * ( R.dv + 1 ) + b] = m.seg[i][a];
      }
  return R;
}

SmMxSrf sm_MxSrfTranspose(const SmMxSrf & S)
{
  SmMxSrf R;
  R.du = S.dv; R.dv = S.du; R.ub = S.vb; R.vb = S.ub;
  const size_t nus = S.ub.size() - 1, nvs = S.vb.size() - 1;
  R.P.assign( nus * nvs, std::vector<SmMxV4>( ( R.du + 1 ) * ( R.dv + 1 )));
  for (size_t i=0; i<nus; i++)
    for (size_t j=0; j<nvs; j++)
      for (int a=0; a<=S.du; a++)
        for (int b=0; b<=S.dv; b++)
          R.P[j * nus + i][b * ( R.dv + 1 ) + a] = S.At( i, j )[a * ( S.dv + 1 ) + b];
  return R;
}

SmMxSrf sm_MxSrfMapV(const SmMxSrf & S, std::function<SmMxCrv(const SmMxCrv &)> f)
{
  return sm_MxSrfTranspose( sm_MxSrfMapU( sm_MxSrfTranspose( S ), f ));
}

SmMxSrf sm_MxSrfCommon(const SmMxSrf & S, const std::vector<double> & ub, int du, const std::vector<double> & vb, int dv)
{
  SmMxSrf R = sm_MxSrfMapU( S, [&](const SmMxCrv & c) { return sm_MxElevate( sm_MxRefine( c, ub ), du ); } );
  return sm_MxSrfMapV( R, [&](const SmMxCrv & c) { return sm_MxElevate( sm_MxRefine( c, vb ), dv ); } );
}

// sum over k of crU[k](u) times the cubic Bernstein row crH[k](v)
SmMxSrf sm_MxTensorU(std::vector<SmMxCrv> crU, const double crH[4][4])
{
  for (size_t k=1; k<crU.size(); k++) sm_MxCommon( crU[0], crU[k] );
  for (size_t k=1; k<crU.size(); k++) sm_MxCommon( crU[0], crU[k] );
  SmMxSrf S;
  S.du = crU[0].Deg();
  S.dv = 3;
  S.ub = crU[0].brk;
  S.vb = { 0.0, 1.0 };
  for (size_t i=0; i<crU[0].seg.size(); i++)
    {
      std::vector<SmMxV4> P( ( S.du + 1 ) * 4, SmMxV4{} );
      for (int a=0; a<=S.du; a++)
        for (int b=0; b<4; b++)
          for (size_t k=0; k<crU.size(); k++)
            P[a * 4 + b] = sm_MxAdd( P[a * 4 + b], sm_MxMul( crU[k].seg[i][a], crH[k][b] ));
      S.P.push_back( P );
    }
  return S;
}

SmMxSrf sm_MxSrfSum(SmMxSrf A, SmMxSrf B, double sb)
{
  const std::vector<double> ub = sm_MxUnion( A.ub, B.ub ), vb = sm_MxUnion( A.vb, B.vb );
  const int du = std::max( A.du, B.du ), dv = std::max( A.dv, B.dv );
  A = sm_MxSrfCommon( A, ub, du, vb, dv );
  B = sm_MxSrfCommon( B, ub, du, vb, dv );
  for (size_t p=0; p<A.P.size(); p++)
    for (size_t k=0; k<A.P[p].size(); k++)
      A.P[p][k] = sm_MxAdd( A.P[p][k], sm_MxMul( B.P[p][k], sb ));
  return A;
}

// the patch of four sides (v = 0, v = 1, u = 0, u = 1) and their ribbons
SmMxSrf sm_MxCoons(const SmMxSide & c0, const SmMxSide & c1, const SmMxSide & d0, const SmMxSide & d1)
{
  // cubic Hermite blends in Bernstein form: H0, H1 (values), H~0, H~1 (derivatives)
  const double H[4][4] = { { 1, 1, 0, 0 }, { 0, 0, 1, 1 }, { 0, 1.0 / 3.0, 0, 0 }, { 0, 0, -1.0 / 3.0, 0 } };
  const SmMxSrf P1 = sm_MxTensorU( { c0.c, c1.c, sm_MxRibbon( c0 ), sm_MxRibbon( c1 ) }, H );
  const SmMxSrf P2 = sm_MxSrfTranspose( sm_MxTensorU( { d0.c, d1.c, sm_MxRibbon( d0 ), sm_MxRibbon( d1 ) }, H ));
  // the corners' data: positions, S_u, S_v, S_uv
  SmMxV4 Sp[2][2], Su[2][2], Sv[2][2], Suv[2][2];
  const SmMxSide * aC[2] = { &c0, &c1 };
  const SmMxSide * aD[2] = { &d0, &d1 };
  for (int i=0; i<2; i++)
    for (int j=0; j<2; j++)
      {
        Sp[i][j] = sm_MxEval( aC[j]->c, i );
        Su[i][j] = sm_MxCrvDerivAt( aC[j]->c, i, 1 );
        Sv[i][j] = sm_MxCrvDerivAt( aD[i]->c, j, 1 );
        Suv[i][j] = sm_MxCrvDerivAt( sm_MxRibbon( *aC[j] ), i, 1 );
      }
  SmMxSrf P12;
  P12.du = P12.dv = 3;
  P12.ub = P12.vb = { 0.0, 1.0 };
  std::vector<SmMxV4> P( 16, SmMxV4{} );
  for (int a=0; a<4; a++)
    for (int b=0; b<4; b++)
      for (int i=0; i<2; i++)
        for (int j=0; j<2; j++)
          {
            P[a*4+b] = sm_MxAdd( P[a*4+b], sm_MxMul( Sp[i][j],  H[i][a]   * H[j][b] ));
            P[a*4+b] = sm_MxAdd( P[a*4+b], sm_MxMul( Su[i][j],  H[2+i][a] * H[j][b] ));
            P[a*4+b] = sm_MxAdd( P[a*4+b], sm_MxMul( Sv[i][j],  H[i][a]   * H[2+j][b] ));
            P[a*4+b] = sm_MxAdd( P[a*4+b], sm_MxMul( Suv[i][j], H[2+i][a] * H[2+j][b] ));
          }
  P12.P.push_back( P );
  return sm_MxSrfSum( sm_MxSrfSum( P1, P2, 1.0 ), P12, -1.0 );
}

// the bilinearly blended Coons patch of four sides (G0)
SmMxSrf sm_MxCoonsBilinear(const SmMxSide & c0, const SmMxSide & c1, const SmMxSide & d0, const SmMxSide & d1)
{
  const double H[2][4] = { { 1, 0, 0, 0 }, { 0, 1, 0, 0 } };
  auto tensor = [](std::vector<SmMxCrv> crU) // crU[k](u) times (1 - v), v
    {
      sm_MxCommon( crU[0], crU[1] );
      SmMxSrf S;
      S.du = crU[0].Deg();
      S.dv = 1;
      S.ub = crU[0].brk;
      S.vb = { 0.0, 1.0 };
      for (size_t i=0; i<crU[0].seg.size(); i++)
        {
          std::vector<SmMxV4> P( ( S.du + 1 ) * 2 );
          for (int a=0; a<=S.du; a++) { P[a * 2] = crU[0].seg[i][a]; P[a * 2 + 1] = crU[1].seg[i][a]; }
          S.P.push_back( P );
        }
      return S;
    };
  (void)H;
  const SmMxSrf P1 = tensor( { c0.c, c1.c } );
  const SmMxSrf P2 = sm_MxSrfTranspose( tensor( { d0.c, d1.c } ));
  SmMxSrf P12;
  P12.du = P12.dv = 1;
  P12.ub = P12.vb = { 0.0, 1.0 };
  P12.P.push_back( { sm_MxEval( c0.c, 0.0 ), sm_MxEval( c1.c, 0.0 ), sm_MxEval( c0.c, 1.0 ), sm_MxEval( c1.c, 1.0 ) } );
  return sm_MxSrfSum( sm_MxSrfSum( P1, P2, 1.0 ), P12, -1.0 );
}

SmStatus sm_MxSrfToBSpline(const SmContext & crContext, const SmMxSrf & S, SmBSplineSurface *& rpOut)
{
  rpOut = NULL;
  const size_t nus = S.ub.size() - 1, nvs = S.vb.size() - 1;
  const size_t nU = nus * S.du + 1, nV = nvs * S.dv + 1;
  SmTArray<SmPoint3d> sPts;
  SmTArray<double> sW;
  for (size_t I=0; I<nU; I++)
    for (size_t J=0; J<nV; J++)
      {
        const size_t i = std::min( I / S.du, nus - 1 ), j = std::min( J / S.dv, nvs - 1 );
        const SmMxV4 & q = S.At( i, j )[( I - i * S.du ) * ( S.dv + 1 ) + ( J - j * S.dv )];
        if ( !( q[3] > 0.0 ))
          { SER_MSG( SM_ERR, _T("Mixed corner: a face's weights are not all positive") ); }
        sPts.Add( sm_MxProject( q ));
        sW.Add( q[3] );
      }
  SmTArray<ULONG> sUM, sVM;
  SmTArray<double> sUK, sVK;
  for (size_t i=0; i<S.ub.size(); i++) { sUK.Add( S.ub[i] ); sUM.Add( ( i == 0 || i + 1 == S.ub.size() ) ? S.du + 1 : S.du ); }
  for (size_t j=0; j<S.vb.size(); j++) { sVK.Add( S.vb[j] ); sVM.Add( ( j == 0 || j + 1 == S.vb.size() ) ? S.dv + 1 : S.dv ); }
  SER( SmBSplineSurface::CreateCanonical( crContext, S.du, S.dv, sPts, SM_SF_UNSPECIFIED, sUM, sVM, sUK, sVK, SM_KT_UNSPECIFIED, &sW, NULL, rpOut ));
  NER( rpOut );
  return SM_SUCCESS;
}

// the piece of a curve between parameters a and b, on [0, 1] from a to b
SmMxCrv sm_MxPiece(const SmMxCrv & c, double a, double b)
{
  const double lo = std::min( a, b ), hi = std::max( a, b );
  SmMxCrv r = sm_MxRefine( c, { lo, hi } ), p;
  for (size_t s=0; s<r.seg.size(); s++)
    if ( r.brk[s] >= lo - 1e-13 * ( 1.0 + fabs( lo )) && r.brk[s+1] <= hi + 1e-13 * ( 1.0 + fabs( hi )))
      {
        if ( p.brk.empty() ) p.brk.push_back( r.brk[s] );
        p.seg.push_back( r.seg[s] );
        p.brk.push_back( r.brk[s+1] );
      }
  return sm_MxUnit( p, a > b ? TRUE : FALSE );
}

// a polynomial uv curve (B-spline) from a 2D piecewise Bezier curve
SmStatus sm_MxToUV(const SmContext & crContext, const SmMxCrv & c, SmBSplineCurve *& rpOut)
{
  rpOut = NULL;
  const int p = c.Deg();
  SmTArray<SmPoint3d> sPts;
  for (size_t s=0; s<c.seg.size(); s++)
    for (int i=( s == 0 ? 0 : 1 ); i<=p; i++)
      sPts.Add( SmPoint3d( c.seg[s][i][0] / c.seg[s][i][3], c.seg[s][i][1] / c.seg[s][i][3], 0.0 ));
  SmTArray<ULONG> sMults;
  SmTArray<double> sKnots;
  for (size_t i=0; i<c.brk.size(); i++)
    {
      sKnots.Add( c.brk[i] );
      sMults.Add( ( i == 0 || i + 1 == c.brk.size() ) ? p + 1 : p );
    }
  SER( SmBSplineCurve::CreateCanonical( crContext, 2, p, sPts, SM_CF_UNSPECIFIED, sMults, sKnots, SM_KT_UNSPECIFIED, NULL, NULL, rpOut ));
  NER( rpOut );
  return SM_SUCCESS;
}

} // namespace

/*******************************************************************//**
PURPOSE: Whether a filleted vertex is a Mixed corner (see the file's head):
    three edges, two blended of opposite convexity continuing each other G1
    there, the third a straight tangent seam between their other faces, the
    blends' common face a plane.
***********************************************************************/
SmStatus SmBlendCorner::ClassifyMixed
  (const SmVertex             * cpVertex,
   const SmTArray<SmEdgeuse*> & crSolverEUs,
   const SmTArray<SmEdgeuse*> & crAllEUs,
   SmFilletExecutive          * pExec,
   double                       dTangencyTolDeg,
   SmBlendCornerEvent         & reEvent)
{
  reEvent = SM_BCE_NONE;
  NER( cpVertex );
  NER( pExec );
  SmMixedRing sRing;
  SmBoolean   bMixed = FALSE;
  SER( sm_MixedRing( cpVertex, crSolverEUs, crAllEUs, pExec, sRing, bMixed ));
  if ( !bMixed )
    { return SM_SUCCESS; }
  if (   !sm_TestG1Continuity( cpVertex, sRing.pConvex->GetEdge(), sRing.pConcave->GetEdge(), dTangencyTolDeg, NULL, 0 )
      || !sm_CheckTangentialEdgeuseAtVertex( sRing.pSeam, cpVertex, dTangencyTolDeg ))
    { return SM_SUCCESS; }
  const double dTol = smos_Min( sRing.apFS[0]->GetThisApproxTol3d(), sRing.apFS[1]->GetThisApproxTol3d() );
  SmPoint3d  sLinePt;
  SmVector3d sAxis;
  if (   sRing.pT->GetSurface() == NULL || !sRing.pT->GetSurface()->IsPlanar( dTol )
      || sRing.pSeam->GetEdge()->GetCurve() == NULL || !sRing.pSeam->GetEdge()->GetCurve()->IsLine( 10, dTol, sLinePt, sAxis ))
    { return SM_SUCCESS; }
  reEvent = SM_BCE_MIXED;
  return SM_SUCCESS;

} // end SmBlendCorner::ClassifyMixed

/*******************************************************************//**
PURPOSE: The Mixed corner's topology (see the file's head): vertices X (the
    rails' crossing on T), A (B1's rail end on F1), P (on the seam), V,
    Q (on F2), B (B2's rail end on F2); edges X-A (B1's trimline), A-P
    (on F1), P-V (the seam), V-X (the inflection, in T), V-Q, Q-B (on F2),
    B-X (B2's trimline).
***********************************************************************/
SmStatus SmBlendCorner::MakeMixedTopology()
{
  NER( m_pSolverEUs );
  NER( m_pAllEUs );
  SmMixedRing sRing;
  SmBoolean   bMixed = FALSE;
  SER( sm_MixedRing( m_cpVertex, *m_pSolverEUs, *m_pAllEUs, GetFilletExecutive(), sRing, bMixed ));
  if ( !bMixed )
    { SER( SM_ERR ); }
  SmFilletGeom * apFG[2];
  for (ULONG ii=0; ii<2; ii++)
    {
      apFG[ii] = sRing.apFS[ii]->GetFirstFilletGeom();
      NER( apFG[ii] );
    }
  const ULONG alOnT[2]    = { sRing.apFS[0]->FindIndexOfRailOnFace( sRing.pT ),  sRing.apFS[1]->FindIndexOfRailOnFace( sRing.pT ) };
  const ULONG alOnWall[2] = { sRing.apFS[0]->FindIndexOfRailOnFace( sRing.pF1 ), sRing.apFS[1]->FindIndexOfRailOnFace( sRing.pF2 ) };
  if ( alOnT[0] == alOnWall[0] || alOnT[1] == alOnWall[1] )
    { SER( SM_ERR ); }

  // the blends run on past the vertex (no blending corner): their rails on the wall then cross the
  // seam's line, at P and Q
  m_bBlending = FALSE;
  m_bBevel    = FALSE;

  // vertices: X, A, P, V, Q, B
  SmFilletVertex * apV[6];
  for (ULONG ii=0; ii<6; ii++)
    {
      apV[ii] = new (m_pPseudoBrep) SmFilletVertex( this );
      NER( apV[ii] );
      m_vVertices.Add( apV[ii] );
    }
  SmFilletVertex * pX = apV[0], * pA = apV[1], * pP = apV[2], * pV = apV[3], * pQ = apV[4], * pB = apV[5];
  pX->SetPointClass( SM_PC_FACE, sRing.pT );
  pX->SetFilletVertexType( SM_FV_RAIL_X_RAIL );
  pA->SetPointClass( SM_PC_FACE, sRing.pF1 );
  pP->SetPointClass( SM_PC_EDGEUSE, sRing.pSeam );
  pV->SetPointClass( SM_PC_VERTEX, (SmObject*)m_cpVertex );
  pV->SetFilletVertexType( SM_FV_ON_VERTEX );
  pQ->SetPointClass( SM_PC_FACE, sRing.pF2 );
  pB->SetPointClass( SM_PC_FACE, sRing.pF2 );
  SmFilletEdgeuse * pNewEU = NULL;
  SER( apFG[0]->MakeRailEdgeuse( m_pPseudoBrep, alOnT[0], pX, pNewEU ));
  SER( apFG[1]->MakeRailEdgeuse( m_pPseudoBrep, alOnT[1], pX, pNewEU ));
  SER( apFG[0]->MakeRailEdgeuse( m_pPseudoBrep, alOnWall[0], pA, pNewEU ));
  SER( apFG[1]->MakeRailEdgeuse( m_pPseudoBrep, alOnWall[1], pB, pNewEU ));

  // edges: X-A, A-P, P-V, V-X, V-Q, Q-B, X-B (each trimline from X, as its face's iso runs)
  const ULONG alFrom[7] = { 0, 1, 2, 3, 3, 4, 0 };
  const ULONG alTo[7]   = { 1, 2, 3, 0, 4, 5, 5 };
  SmFace * apOn[7]      = { NULL, sRing.pF1, sRing.pF1, NULL, sRing.pF2, sRing.pF2, NULL };
  for (ULONG ii=0; ii<7; ii++)
    {
      SmFilletEdge * pEdge = NULL;
      SER( SmFilletEdge::MakeFilletEdge( m_pPseudoBrep, apV[alFrom[ii]], apV[alTo[ii]], pEdge, this ));
      NER( pEdge );
      pEdge->SetFilletEdgeType( SM_FE_PRECOMPUTED );
      if ( apOn[ii] != NULL )
        { pEdge->SetOriginalFace( apOn[ii] ); }
      if ( ii == 0 || ii == 6 )
        {
          SmFilletEdgeuse * pPrimEU = (SmFilletEdgeuse*)pEdge->GetPrimaryEdgeuse();
          NER( pPrimEU );
          pPrimEU->SetFilletGeom( apFG[ii == 0 ? 0 : 1] );
        }
      m_vEdges.Add( pEdge );
    }
  return SM_SUCCESS;

} // end SmBlendCorner::MakeMixedTopology

/*******************************************************************//**
PURPOSE: The Mixed corner's vertices once both blends are traced (X, the
    rails' crossing, and V by their own kinds): A and B, where each blend's
    rail on the wall meets the other's rail swept along T's normal; P, where
    B1's rail on F1 meets the seam; Q, where B2's rail on F2 crosses the
    seam's line past V.
***********************************************************************/
SmStatus SmBlendCorner::CalcMixedVertGeom()
{
  if ( m_vVertices.GetSize() != 6 || m_vEdges.GetSize() != 7 )
    { SER( SM_ERR ); }
  SmMixedRing sRing;
  SmBoolean   bMixed = FALSE;
  SER( sm_MixedRing( m_cpVertex, *m_pSolverEUs, *m_pAllEUs, GetFilletExecutive(), sRing, bMixed ));
  if ( !bMixed )
    { SER( SM_ERR ); }
  SmFilletGeom * apFG[2] = { sRing.apFS[0]->GetFirstFilletGeom(), sRing.apFS[1]->GetFirstFilletGeom() };
  if ( apFG[0] == NULL || apFG[1] == NULL || apFG[0]->GetFilletSurface() == NULL || apFG[1]->GetFilletSurface() == NULL )
    { return SM_SUCCESS; } // not traced yet: the corner's second pass places them
  SmFilletVertex * pX = m_vVertices[0], * pA = m_vVertices[1], * pP = m_vVertices[2], * pV = m_vVertices[3],
                 * pQ = m_vVertices[4], * pB = m_vVertices[5];
  if ( !pX->IsProcessed() || !pV->IsProcessed() )
    { SER( SM_ERR ); }
  const double dTol = GetThisApproxTol3d();

  SmPoint3d  sOnT;
  SmVector3d sNormal;
  if ( !sRing.pT->GetSurface()->IsPlanar( dTol, &sOnT, &sNormal ))
    { SER( SM_ERR ); }
  SER( sNormal.Unitize() );
  const SmPoint3d sV = pV->GetPoint();

  const ULONG alOnT[2]    = { sRing.apFS[0]->FindIndexOfRailOnFace( sRing.pT ),  sRing.apFS[1]->FindIndexOfRailOnFace( sRing.pT ) };
  const ULONG alOnWall[2] = { sRing.apFS[0]->FindIndexOfRailOnFace( sRing.pF1 ), sRing.apFS[1]->FindIndexOfRailOnFace( sRing.pF2 ) };
  SmCurve * apRailT[2]    = { apFG[0]->GetRail( alOnT[0] )->GetCurve(),    apFG[1]->GetRail( alOnT[1] )->GetCurve() };
  SmCurve * apRailWall[2] = { apFG[0]->GetRail( alOnWall[0] )->GetCurve(), apFG[1]->GetRail( alOnWall[1] )->GetCurve() };
  NER( apRailT[0] ); NER( apRailT[1] ); NER( apRailWall[0] ); NER( apRailWall[1] );

  const double dReach = 4.0 * ( sV.DistanceBetween( pX->GetPoint() ) + 1.0 );

  // A and B: each wall rail meets the other blend's rail on T swept along T's normal
  for (ULONG ii=0; ii<2; ii++)
    {
      SmBSplineSurface * pWall = NULL;
      SER( sm_MixedRailWall( m_crContext, *apRailT[1-ii], sNormal, dReach, pWall ));
      SmObjDelete sClean( pWall );
      SmPoint3d sPnt;
      double    dT = 0.0;
      SER( sm_MixedCurveSurface( *apRailWall[ii], *pWall, pX->GetPoint(), dTol, sPnt, dT ));
      ( ii == 0 ? pA : pB )->SetPoint( sPnt );
      ( ii == 0 ? pA : pB )->SetStatus( SM_FIL_PROCESSED );
    }

  // P: B1's wall rail crosses the seam; Q: B2's wall rail crosses the seam's line past V
  SmEdge  * pSeam = sRing.pSeam->GetEdge();
  SmVector3d sAlong;
  SER( sm_MixedLeaving( pSeam, m_cpVertex, sAlong ));
  {
    SmBSplineCurve * pLine = NULL;
    SER( SmBSplineCurve::CreateLineSegment( m_crContext, 3, sV + sAlong * dReach, sV - sAlong * dReach, pLine ));
    SmObjDelete sClean( pLine );
    for (ULONG ii=0; ii<2; ii++)
      {
        SmSolutionArray sSolutions;
        SER( apRailWall[ii]->GlobalCurveIntersect( apRailWall[ii]->GetNaturalInterval(), *pLine, pLine->GetNaturalInterval(), dTol, sSolutions ));
        double    dBest = SM_BIG_DOUBLE;
        SmPoint3d sBest;
        for (ULONG k=0; k<sSolutions.GetSize(); k++)
          {
            SmPoint3d sPnt;
            SER( pLine->EvaluatePoint( sSolutions[k].m_vStart.m_adParameters[1], sPnt ));
            // P on the seam's side of V (along it), Q past V
            if ( ( ii == 0 ) != ( ( sPnt - sV ).Dot( sAlong ) > 0.0 ))
              { continue; }
            if ( sPnt.DistanceBetween( sV ) < dBest )
              { dBest = sPnt.DistanceBetween( sV ); sBest = sPnt; }
          }
        if ( dBest == SM_BIG_DOUBLE )
          { SER_MSG( SM_ERR, _T("Mixed corner: a wall rail does not cross the seam's line") ); }
        ( ii == 0 ? pP : pQ )->SetPoint( sBest );
        ( ii == 0 ? pP : pQ )->SetStatus( SM_FIL_PROCESSED );
      }
  }
  return SM_SUCCESS;

} // end SmBlendCorner::CalcMixedVertGeom

/*******************************************************************//**
PURPOSE: The Mixed corner's edges (see the file's head) and its two faces'
    surfaces, L (X, A, P, V) and U (X, V, Q, B).

NOTES: The trimlines: each blend B_i meets the other's rail wall in a
    certified intersection; its uv curve on B_i, held at X and at its rail end
    exactly on the rails' iso lines, is the trimline's definition, and the
    edge's curve is B_i over it (sm_MxCompose: exactly on B_i). The rail
    pieces A-P and Q-B are the rails' own, the lines P-V and V-Q the seam's,
    the inflection V-X the arc in T tangent at V to the edges' direction.
    Each face is the Coons patch of its sides (sm_MxCoons), every ribbon G1
    to its neighbour by its D:
    - on B1 (a round or chamfer along a straight edge, ruled along it): the
      edge's direction;
    - on B2 (about the concave edge's axis, a surface of revolution): the
      circumferential field, Omega x (P - O);
    - on F1 along A-P: F1's direction square to the rail;
    - on the seam's line (P-V, V-Q): the edges' common direction at V (the
      wall's tangent plane is constant along the line);
    - on F2 along Q-B (a parallel of F2, coaxial): F2's meridian there;
    - along V-X, shared by L and U: M, from T's direction square to the arc
      at X to the seam's at V, its end derivatives solved with the twists at
      X and V of both faces at once.
***********************************************************************/
SmStatus SmBlendCorner::CalcMixedEdgeGeom()
{
  SmMixedRing sRing;
  SmBoolean   bMixed = FALSE;
  SER( sm_MixedRing( m_cpVertex, *m_pSolverEUs, *m_pAllEUs, GetFilletExecutive(), sRing, bMixed ));
  if ( !bMixed )
    { SER( SM_ERR ); }
  for (ULONG ii=0; ii<m_vVertices.GetSize(); ii++)
    {
      if ( !m_vVertices[ii]->IsProcessed() )
        { SER( SM_ERR ); }
    }
  SmFilletGeom * apFG[2] = { sRing.apFS[0]->GetFirstFilletGeom(), sRing.apFS[1]->GetFirstFilletGeom() };
  NER( apFG[0] );
  NER( apFG[1] );
  const double dTol = GetThisApproxTol3d();
  // chamfers meet their faces G0: the corner's faces are then the bilinear Coons patches of their sides
  const SmBoolean bChamfer = (   sRing.apFS[0]->GetFilletSurfaceGenerator()->GetFilletSurfaceGeneratorType() == SM_FSG_LINEAR
                              || sRing.apFS[1]->GetFilletSurfaceGenerator()->GetFilletSurfaceGeneratorType() == SM_FSG_LINEAR ) ? TRUE : FALSE;
  SmPoint3d  sOnT;
  SmVector3d sNormal;
  if ( !sRing.pT->GetSurface()->IsPlanar( dTol, &sOnT, &sNormal ))
    { SER( SM_ERR ); }
  SER( sNormal.Unitize() );
  SmFilletVertex * pX = m_vVertices[0], * pA = m_vVertices[1], * pP = m_vVertices[2], * pV = m_vVertices[3],
                 * pQ = m_vVertices[4], * pB = m_vVertices[5];
  const SmPoint3d sX = pX->GetPoint(), sV = pV->GetPoint();
  const double dSize  = 1.0 + sX.GetMaxDimension() + sV.DistanceBetween( sX );
  const double dReach = 4.0 * ( sV.DistanceBetween( sX ) + 1.0 );
  const ULONG alOnT[2]    = { sRing.apFS[0]->FindIndexOfRailOnFace( sRing.pT ),  sRing.apFS[1]->FindIndexOfRailOnFace( sRing.pT ) };
  const ULONG alOnWall[2] = { sRing.apFS[0]->FindIndexOfRailOnFace( sRing.pF1 ), sRing.apFS[1]->FindIndexOfRailOnFace( sRing.pF2 ) };
  SmCurve * apRailT[2]    = { apFG[0]->GetRail( alOnT[0] )->GetCurve(),    apFG[1]->GetRail( alOnT[1] )->GetCurve() };
  SmCurve * apRailWall[2] = { apFG[0]->GetRail( alOnWall[0] )->GetCurve(), apFG[1]->GetRail( alOnWall[1] )->GetCurve() };
  NER( apRailT[0] ); NER( apRailT[1] ); NER( apRailWall[0] ); NER( apRailWall[1] );

  // the edges' directions at V, the seam's, the concave edge's axis (Omega through O)
  SmVector3d sChain, sSeamDir;
  SER( sm_MixedLeaving( sRing.pConvex->GetEdge(), m_cpVertex, sChain ));
  SER( sm_MixedLeaving( sRing.pSeam->GetEdge(), m_cpVertex, sSeamDir ));
  sChain = -sChain; // along the convex edge into V, on along the concave one
  SmVector3d sConvexDir = sChain;
  SmPoint3d  sO;
  SmVector3d sOmega = sNormal;
  {
    // the concave edge's circle: its centre from three of its points
    SmCurve * pArc = sRing.pConcave->GetEdge()->GetCurve();
    NER( pArc );
    const SmExtent1d sIvl = sRing.pConcave->GetEdge()->GetInterval();
    SmPoint3d p0, p1, p2;
    SER( pArc->EvaluatePoint( sIvl.Evaluate( 0.1 ), p0 ));
    SER( pArc->EvaluatePoint( sIvl.Evaluate( 0.5 ), p1 ));
    SER( pArc->EvaluatePoint( sIvl.Evaluate( 0.9 ), p2 ));
    const SmVector3d a = p1 - p0, b = p2 - p0;
    const SmVector3d n = a * b;
    const double dN = n.Dot( n );
    if ( dN <= SM_EFF_ZERO )
      { SER( SM_ERR ); }
    // the circumcentre: p0 + (|a|^2 (b x n) + |b|^2 (n x a)) / (2 |n|^2)
    sO = p0 + ( ( b * n ) * a.Dot( a ) + ( n * a ) * b.Dot( b )) * ( 0.5 / dN );
  }

  // the outward side of T at X: B1 (convex) lies below it
  double dOut = 1.0;

  // the trimlines: uv on B_i from X to its rail end, then B_i over it
  SmMxCrv asTrim[2], asPc[2];
  for (ULONG ii=0; ii<2; ii++)
    {
      SmBSplineSurface * pWall = NULL;
      SER( sm_MixedRailWall( m_crContext, *apRailT[1-ii], sNormal, dReach, pWall ));
      SmObjDelete sCleanWall( pWall );
      SmFilletVertex * pTo = ii == 0 ? pA : pB;
      SmVector3d sDir = pTo->GetPoint() - sX;
      SmBSplineCurve * p3d = NULL, * pUV1 = NULL, * pUV2 = NULL;
      SER( sm_SrfSrfIntersection( m_crContext, pX, pTo, &sDir, NULL, apFG[ii]->GetFilletSurface(), pWall,
                                  dTol, sRing.apFS[ii]->GetThisAngTolRad(), p3d, pUV1, pUV2, FALSE, NULL, 0 ));
      SmObjDelete sC1( p3d ), sC2( pUV1 ), sC3( pUV2 );
      NER( pUV1 );
      SmBSplineSurface * pRound = SM_CAST_PTR( SmBSplineSurface, apFG[ii]->GetFilletSurface() );
      NER( pRound );
      SmMxCrv pc;
      SER( sm_MxFromBSpline( *pUV1, pc ));
      // which end is at X
      SmPoint3d sEnd0;
      {
        const SmMxV4 q = sm_MxEval( pc, pc.brk.front() );
        SER( pRound->EvaluatePoint( SmPoint2d( q[0] / q[3], q[1] / q[3] ), sEnd0 ));
      }
      if ( sEnd0.DistanceBetween( sX ) > sEnd0.DistanceBetween( pTo->GetPoint() ))
        { pc = sm_MxUnit( pc, TRUE ); }
      else
        { pc = sm_MxUnit( pc, FALSE ); }
      // its ends on the rails' iso lines: the rails are the round's v = min, max lines, along u
      const SmExtent2d sDom = pRound->GetNaturalUVDomain();
      for (int e=0; e<2; e++)
        {
          SmCurve * pRail = e == 0 ? apRailT[ii] : apRailWall[ii];
          const SmPoint3d sAt = e == 0 ? sX : pTo->GetPoint();
          double dT = 0.0;
          SER( sm_MixedParam( *pRail, sAt, dTol, dT ));
          SmPoint3d sR, sLo, sHi;
          SER( pRail->EvaluatePoint( dT, sR ));
          SER( pRound->EvaluatePoint( SmPoint2d( dT, sDom.GetMin().y ), sLo ));
          SER( pRound->EvaluatePoint( SmPoint2d( dT, sDom.GetMax().y ), sHi ));
          const double dIso = sLo.DistanceBetween( sR ) <= sHi.DistanceBetween( sR ) ? sDom.GetMin().y : sDom.GetMax().y;
          SmMxV4 & q = e == 0 ? pc.seg.front().front() : pc.seg.back().back();
          q = SmMxV4{ dT, dIso, 0.0, 1.0 };
        }
      SmMxSrf S;
      SER( sm_MxSrfFromBSpline( *pRound, S ));
      SmMxCrv c;
      SER( sm_MxCompose( S, pc, c ));
      asTrim[ii] = sm_MxUnitWeights( c );
      asPc[ii] = pc;
      if ( ii == 0 )
        {
          const SmMxV4 m = sm_MxEval( asTrim[0], 0.5 );
          dOut = ( ( sm_MxProject( m ) - sX ).Dot( sNormal ) < 0.0 ) ? 1.0 : -1.0;
        }
    }
  const SmVector3d sOut = sNormal * dOut;
  // the rail ends exactly where the trimlines end
  pA->SetPoint( sm_MxProject( sm_MxEval( asTrim[0], 1.0 )));
  pB->SetPoint( sm_MxProject( sm_MxEval( asTrim[1], 1.0 )));

  // the rail pieces A-P (B1's), B-Q (B2's), the lines V-P, V-Q, and the arc X-V
  SmMxCrv sAP, sBQ, sVP, sVQ, sXV, sXVq;
  {
    SmMxCrv r0, r1;
    SER( sm_MxFromBSpline( *SM_CAST_PTR( SmBSplineCurve, apRailWall[0] ), r0 ));
    SER( sm_MxFromBSpline( *SM_CAST_PTR( SmBSplineCurve, apRailWall[1] ), r1 ));
    double dA = 0.0, dP = 0.0, dB = 0.0, dQ = 0.0;
    SER( sm_MixedParam( *apRailWall[0], pA->GetPoint(), dTol, dA ));
    SER( sm_MixedParam( *apRailWall[0], pP->GetPoint(), dTol, dP ));
    SER( sm_MixedParam( *apRailWall[1], pB->GetPoint(), dTol, dB ));
    SER( sm_MixedParam( *apRailWall[1], pQ->GetPoint(), dTol, dQ ));
    sAP = sm_MxUnitWeights( sm_MxPiece( r0, dA, dP ));
    sBQ = sm_MxUnitWeights( sm_MxPiece( r1, dB, dQ ));
    // their ends exactly the vertices (A and B the trimlines' ends)
    sAP.seg.front().front() = sm_MxPoint( pA->GetPoint() );
    sBQ.seg.front().front() = sm_MxPoint( pB->GetPoint() );
    pP->SetPoint( sm_MxProject( sAP.seg.back().back() ));
    pQ->SetPoint( sm_MxProject( sBQ.seg.back().back() ));
  }
  auto line = [](const SmPoint3d & a, const SmPoint3d & b) { SmMxCrv c; c.brk = { 0.0, 1.0 }; c.seg.push_back( SmMxBez{ sm_MxPoint( a ), sm_MxPoint( b ) } ); return c; };
  sVP = line( sV, pP->GetPoint() );
  sVQ = line( sV, pQ->GetPoint() );
  {
    // V-X leaves X along an asymptotic direction of the corner's shape operator there: it lies in T,
    // so each face's normal curvature along it at X is 0, and the two blends' fix the operator (II
    // along each trimline's direction e_i is B_i's; across, B1, ruled along D, gives
    // II(e1, e2) = lambda II(e1, e1) for e2 = lambda e1 + mu D). It arrives at V along the edges.
    // The curve is the conic in T with those end tangents (its middle control point where they
    // meet), weighted as a circle's arc of that turn. Chamfers (G0) take the wedge's bisector.
    SmVector3d asE[2];
    double     adK[2];
    for (int ii=0; ii<2; ii++)
      {
        const SmMxV4 c = sm_MxEval( asTrim[ii], 0.0 ), c1 = sm_MxCrvDerivAt( asTrim[ii], 0.0, 1 ), c2 = sm_MxCrvDerivAt( asTrim[ii], 0.0, 2 );
        const SmVector3d P( c[0] / c[3], c[1] / c[3], c[2] / c[3] );
        const SmVector3d P1( ( c1[0] - P.x * c1[3] ) / c[3], ( c1[1] - P.y * c1[3] ) / c[3], ( c1[2] - P.z * c1[3] ) / c[3] );
        const SmVector3d P2( ( c2[0] - 2.0 * P1.x * c1[3] - P.x * c2[3] ) / c[3], ( c2[1] - 2.0 * P1.y * c1[3] - P.y * c2[3] ) / c[3],
                             ( c2[2] - 2.0 * P1.z * c1[3] - P.z * c2[3] ) / c[3] );
        const double dSpeed2 = P1.Dot( P1 );
        if ( dSpeed2 <= SM_EFF_ZERO )
          { SER( SM_ERR ); }
        adK[ii] = P2.Dot( sOut ) / dSpeed2;
        asE[ii] = P1 - sOut * P1.Dot( sOut );
        SER( asE[ii].Unitize() );
      }
    SmVector3d sDirX = asE[0] + asE[1];
    if ( !bChamfer )
      {
        // e2 = lambda e1 + mu D (in T)
        SmVector3d sD = sConvexDir - sOut * sConvexDir.Dot( sOut );
        SER( sD.Unitize() );
        const double d11 = 1.0, d12 = asE[0].Dot( sD ), d22 = 1.0;
        const double r1 = asE[1].Dot( asE[0] ), r2 = asE[1].Dot( sD );
        const double dDet = d11 * d22 - d12 * d12;
        if ( smos_Fabs( dDet ) <= SM_EFF_ZERO )
          { SER( SM_ERR ); }
        const double dLambda = ( r1 * d22 - r2 * d12 ) / dDet;
        const double II11 = adK[0], II22 = adK[1], II12 = dLambda * adK[0];
        // II22 k^2 + 2 II12 k + II11 = 0, k > 0: a = e1 + k e2
        const double dDisc = II12 * II12 - II11 * II22;
        double dKBest = -1.0;
        if ( dDisc >= 0.0 && smos_Fabs( II22 ) > SM_EFF_ZERO )
          {
            for (int sgn=-1; sgn<=1; sgn+=2)
              {
                const double k = ( -II12 + sgn * sqrt( dDisc )) / II22;
                if ( k > 0.0 && ( dKBest < 0.0 || k < dKBest ))
                  { dKBest = k; }
              }
          }
        if ( dKBest <= 0.0 )
          { SER_MSG( SM_ERR, _T("Mixed corner: the corner's shape at the rails' crossing has no asymptotic direction between the trimlines") ); }
        sDirX = asE[0] + asE[1] * dKBest;
      }
    SER( sDirX.Unitize() );
    SmVector3d sDirV = sChain - sOut * sChain.Dot( sOut ); // arriving at V
    SER( sDirV.Unitize() );
    // X + s a = V - t c
    const SmVector3d sD = sV - sX;
    const double a11 = sDirX.Dot( sDirX ), a12 = sDirX.Dot( sDirV ), a22 = sDirV.Dot( sDirV );
    const double b1 = sDirX.Dot( sD ), b2 = sDirV.Dot( sD );
    const double dDet = a11 * a22 - a12 * a12;
    if ( smos_Fabs( dDet ) <= SM_EFF_ZERO )
      { SER( SM_ERR ); }
    const double dS = ( b1 * a22 - b2 * a12 ) / dDet, dT = ( a11 * b2 - a12 * b1 ) / dDet;
    if ( !( dS > 0.0 ) || !( dT > 0.0 ))
      { SER_MSG( SM_ERR, _T("Mixed corner: the inflection's end tangents do not meet between X and V") ); }
    const SmPoint3d sI = sX + sDirX * dS;
    const double w = cos( 0.5 * acos( smos_Max( -1.0, smos_Min( 1.0, sDirX.Dot( sDirV )))));
    sXV.brk = { 0.0, 1.0 };
    sXV.seg.push_back( SmMxBez{ sm_MxPoint( sX ), sm_MxPoint( sI, w ), sm_MxPoint( sV ) } );
    sXVq = sXV;
    sXV = sm_MxUnitWeights( sXV );
  }

  // the transversal fields
  auto circumferential = [&](const SmMxV4 & c) { const SmVector3d r( c[0] - c[3] * sO.x, c[1] - c[3] * sO.y, c[2] - c[3] * sO.z ); return sm_MxDir( sOmega * r ); };
  SmVector3d sMeridian;
  double dMa = 0.0, dMb = 0.0, dRho = 0.0;
  {
    // F2's meridian at Q, as a Omega + b rho^
    const SmPoint3d sQ = pQ->GetPoint();
    SmVector3d sRad = ( sQ - sO ) - sOmega * ( sQ - sO ).Dot( sOmega );
    dRho = sRad.Length();
    SER( sRad.Unitize() );
    SmPoint2d sUV;
    SmBoolean bFound = FALSE;
    SmSurface * pF2 = sRing.pF2->GetSurface();
    SmSolutionArray sSol;
    SER( pF2->GlobalPointSolve( pF2->GetNaturalUVDomain(), SM_SO_MINIMIZE, sQ, SmZoneTol3d( dTol ), NULL, SM_SR_SINGLE, sSol ));
    if ( sSol.GetSize() < 1 )
      { SER( SM_ERR ); }
    SmVector3d sN;
    SER( pF2->EvaluateNormal( SmPoint2d( sSol[0].m_vStart[0], sSol[0].m_vStart[1] ), TRUE, TRUE, sN ));
    (void)bFound; (void)sUV;
    sMeridian = sN * ( sOmega * sRad );
    SER( sMeridian.Unitize() );
    dMa = sMeridian.Dot( sOmega );
    dMb = sMeridian.Dot( sRad );
  }
  auto meridian = [&](const SmMxV4 & c)
    {
      SmVector3d r( c[0] - c[3] * sO.x, c[1] - c[3] * sO.y, c[2] - c[3] * sO.z );
      r = r - sOmega * r.Dot( sOmega );
      return sm_MxDir( sOmega * ( dMa * c[3] ) + r * ( dMb / dRho ));
    };
  SmVector3d sF1N;
  {
    SmSurface * pF1 = sRing.pF1->GetSurface();
    SmSolutionArray sSol;
    SER( pF1->GlobalPointSolve( pF1->GetNaturalUVDomain(), SM_SO_MINIMIZE, pA->GetPoint(), SmZoneTol3d( dTol ), NULL, SM_SR_SINGLE, sSol ));
    if ( sSol.GetSize() < 1 )
      { SER( SM_ERR ); }
    SER( pF1->EvaluateNormal( SmPoint2d( sSol[0].m_vStart[0], sSol[0].m_vStart[1] ), TRUE, TRUE, sF1N ));
  }

  // the sides of L (X: (0,0), A: (1,0), P: (1,1), V: (0,1)) and U (X: (0,0), V: (1,0), Q: (1,1), B: (0,1))
  SmMxSide L[4], U[4]; // c0, c1, d0, d1
  L[0].c = asTrim[0];         L[0].D = sm_MxConst( sm_MxDir( sConvexDir ));
  L[1].c = sVP;               L[1].D = sm_MxConst( sm_MxDir( sChain ));
  L[2].c = sXV;               L[2].bShared = TRUE;
  L[3].c = sAP;
  U[0].c = sXV;               U[0].bShared = TRUE;
  U[1].c = sBQ;               U[1].D = sm_MxMap( sBQ, meridian );
  U[2].c = asTrim[1];         U[2].D = sm_MxMap( asTrim[1], circumferential );
  U[3].c = sVQ;               U[3].D = sm_MxConst( sm_MxDir( sChain ));
  {
    // A-P's direction (its tangent at A)
    const SmMxV4 t = sm_MxCrvDerivAt( sAP, 0.0, 1 ), p = sm_MxEval( sAP, 0.0 );
    const SmVector3d d( t[0] - p[0] * t[3], t[1] - p[1] * t[3], t[2] - p[2] * t[3] );
    L[3].D = sm_MxConst( sm_MxDir( sF1N * d ));
  }
  // M along V-X, pointing into U: A + k B with
  //   A = d (n_T x T~), d the curve's distance from the wall's plane at V (into the material), T~ its
  //     homogeneous tangent: horizontal, zero at V together with its derivative (d(V) = 0 and the
  //     curve arrives at V along the wall);
  //   B = t w^2 n_T: zero at X, out of T at V.
  // So the faces' tangent plane turns from T at X to the wall at V as the room to the wall closes, and
  // only k, the rate at X, is left: the curvature compatibility there fixes it (solved with the twists).
  SmMxSide * apShared[2] = { &L[2], &U[0] };
  SmMxCrv sMA, sMB;
  {
    // the wall's plane at V: its normal, into the material (the curve lies on that side)
    SmVector3d sNW = sNormal * sChain;
    if ( sNW.Dot( sX - sV ) < 0.0 ) sNW = -sNW;
    SER( sNW.Unitize() );
    const SmMxCrv & Q = sXVq;
    SmMxCrv dist = sm_MxMap( Q, [&](const SmMxV4 & q) { return SmMxV4{ sNW.x * ( q[0] - q[3] * sV.x ) + sNW.y * ( q[1] - q[3] * sV.y ) + sNW.z * ( q[2] - q[3] * sV.z ), 0, 0, 0 }; } );
    SmMxCrv w = sm_MxMap( Q, [](const SmMxV4 & q) { return SmMxV4{ q[3], 0, 0, 0 }; } );
    SmMxCrv dQ = sm_MxDeriv( Q );
    SmMxCrv dw = sm_MxMap( dQ, [](const SmMxV4 & q) { return SmMxV4{ q[3], 0, 0, 0 }; } );
    SmMxCrv T = sm_MxSum( sm_MxTimes( w, dQ ), sm_MxTimes( dw, Q ), 1.0, -1.0 );
    SmMxCrv nT = sm_MxMap( T, [&](const SmMxV4 & q) { const SmVector3d v = sNormal * SmVector3d( q[0], q[1], q[2] ); return SmMxV4{ v.x, v.y, v.z, 0.0 }; } );
    sMA = sm_MxTimes( dist, nT );
    SmMxCrv lin;
    lin.brk = { 0.0, 1.0 };
    lin.seg.push_back( SmMxBez{ SmMxV4{ 0, 0, 0, 0 }, SmMxV4{ 1, 0, 0, 0 } } );
    SmMxCrv sc = sm_MxTimes( lin, sm_MxTimes( w, w ));
    sMB = sm_MxMap( sc, [&](const SmMxV4 & q) { return SmMxV4{ q[0] * sOut.x, q[0] * sOut.y, q[0] * sOut.z, 0.0 }; } );
    // A into U at X (towards U's trimline)
    const SmMxV4 t = sm_MxCrvDerivAt( asTrim[1], 0.0, 1 ), p = sm_MxEval( asTrim[1], 0.0 );
    const SmVector3d sE2( t[0] - p[0] * t[3], t[1] - p[1] * t[3], t[2] - p[2] * t[3] );
    const SmMxV4 a0 = sm_MxEval( sMA, 0.0 );
    if ( SmVector3d( a0[0], a0[1], a0[2] ).Dot( sE2 ) < 0.0 )
      { sMA = sm_MxMap( sMA, [](const SmMxV4 & q) { return sm_MxMul( q, -1.0 ); } ); }
    for (SmMxSide * pS : apShared)
      {
        pS->Dv[0] = sm_MxEval( sMA, 0.0 );
        pS->Dd[0] = sm_MxCrvDerivAt( sMA, 0.0, 1 );
        pS->Dk[0] = sm_MxCrvDerivAt( sMB, 0.0, 1 );
        pS->Dv[1] = sm_MxEval( sMB, 1.0 );      // gamma at V is taken per unit of k
        pS->Dd[1] = sm_MxCrvDerivAt( sMB, 1.0, 1 );
        pS->Dk[1] = SmMxV4{};
      }
  }
  for (int f=0; f<2; f++)
    for (int k=0; k<4; k++)
      if ( !( f == 0 ? L : U )[k].bShared )
        sm_MxSideEnds( ( f == 0 ? L : U )[k] );

  if ( !bChamfer )
  {
  // first order: each ribbon's end = the adjacent side's tangent there
  auto tangentAt = [&](const SmMxSide & s, int e) { return sm_MxCrvDerivAt( s.c, e, 1 ); };
  for (int f=0; f<2; f++)
    {
      SmMxSide * S = f == 0 ? L : U;
      SER( sm_MxFirstOrder( S[0], 0, tangentAt( S[2], 0 ), dSize ));
      SER( sm_MxFirstOrder( S[0], 1, tangentAt( S[3], 0 ), dSize ));
      SER( sm_MxFirstOrder( S[1], 0, tangentAt( S[2], 1 ), dSize ));
      SER( sm_MxFirstOrder( S[1], 1, tangentAt( S[3], 1 ), dSize ));
      SER( sm_MxFirstOrder( S[2], 0, tangentAt( S[0], 0 ), dSize ));
      SER( sm_MxFirstOrder( S[2], 1, tangentAt( S[1], 0 ), dSize ));
      SER( sm_MxFirstOrder( S[3], 0, tangentAt( S[0], 1 ), dSize ));
      SER( sm_MxFirstOrder( S[3], 1, tangentAt( S[1], 1 ), dSize ));
    }
  // twists: g0'(0) = h0'(0), g0'(1) = h1'(0), g1'(0) = h0'(1), g1'(1) = h1'(1); at X and V both
  // faces at once, with k
  double dK = 0.0;
  {
    std::vector<std::array<SmMxTwistTerm,2>> eXV = { { SmMxTwistTerm{ &L[0], 0, 1.0 }, SmMxTwistTerm{ &L[2], 0, -1.0 } },
                                                      { SmMxTwistTerm{ &U[0], 0, 1.0 }, SmMxTwistTerm{ &U[2], 0, -1.0 } },
                                                      { SmMxTwistTerm{ &L[1], 0, 1.0 }, SmMxTwistTerm{ &L[2], 1, -1.0 } },
                                                      { SmMxTwistTerm{ &U[0], 1, 1.0 }, SmMxTwistTerm{ &U[3], 0, -1.0 } } };
    SER( sm_MxTwists( eXV, dSize, &dK ));
    std::vector<std::array<SmMxTwistTerm,2>> eA = { { SmMxTwistTerm{ &L[0], 1, 1.0 }, SmMxTwistTerm{ &L[3], 0, -1.0 } } };
    SER( sm_MxTwists( eA, dSize, NULL ));
    std::vector<std::array<SmMxTwistTerm,2>> eP = { { SmMxTwistTerm{ &L[1], 1, 1.0 }, SmMxTwistTerm{ &L[3], 1, -1.0 } } };
    SER( sm_MxTwists( eP, dSize, NULL ));
    std::vector<std::array<SmMxTwistTerm,2>> eB = { { SmMxTwistTerm{ &U[1], 0, 1.0 }, SmMxTwistTerm{ &U[2], 1, -1.0 } } };
    SER( sm_MxTwists( eB, dSize, NULL ));
    std::vector<std::array<SmMxTwistTerm,2>> eQ = { { SmMxTwistTerm{ &U[1], 1, 1.0 }, SmMxTwistTerm{ &U[3], 1, -1.0 } } };
    SER( sm_MxTwists( eQ, dSize, NULL ));
  }
  if ( !( dK > 0.0 ))
    { SER_MSG( SM_ERR, _T("Mixed corner: the faces' tangent plane would turn the wrong way along the inflection") ); }
  for (SmMxSide * pS : apShared)
    {
      pS->v[1][2] /= dK;
      pS->d[1][2] /= dK;
    }
  L[2].D = U[0].D = sm_MxSum( sMA, sMB, 1.0, dK );
  }
  else
  {
  L[2].D = U[0].D = sm_MxSum( sMA, sMB, 1.0, 1.0 );
  }

  // the faces, each turned so its normal leaves the material (out of T at X)
  DeleteMixedSurfaces();
  SmBSplineCurve * apSide[2][4] = { { NULL, NULL, NULL, NULL }, { NULL, NULL, NULL, NULL } }; // c0, c1, d0, d1 as iso curves
  SmObjDelete asCleanSide[2][4];
  for (int f=0; f<2; f++)
    {
      SmMxSide * S = f == 0 ? L : U;
      SmMxSrf srf = bChamfer ? sm_MxCoonsBilinear( S[0], S[1], S[2], S[3] ) : sm_MxCoons( S[0], S[1], S[2], S[3] );
      // its normal at X: S_u x S_v = c0' x d0'
      const SmPoint3d sP0 = sm_MxProject( sm_MxEval( S[0].c, 0.0 ));
      auto dir = [&](const SmMxCrv & c) { const SmMxV4 t = sm_MxCrvDerivAt( c, 0.0, 1 ), p = sm_MxEval( c, 0.0 ); return SmVector3d( ( t[0] - sP0.x * t[3] ) / p[3], ( t[1] - sP0.y * t[3] ) / p[3], ( t[2] - sP0.z * t[3] ) / p[3] ); };
      const SmBoolean bTurn = ( ( dir( S[0].c ) * dir( S[2].c )).Dot( sOut ) < 0.0 ) ? TRUE : FALSE;
      if ( bTurn )
        { srf = sm_MxSrfTranspose( srf ); }
      SmBSplineSurface * pPatch = NULL;
      SER( sm_MxSrfToBSpline( m_crContext, srf, pPatch ));
      m_apMixedSurface[f] = pPatch;
      // its loop: v = 0 (SAME), u = 1 (SAME), v = 1 (OPPOSITE), u = 0 (OPPOSITE), each the patch's iso
      const SmExtent2d sDom( 0.0, 0.0, 1.0, 1.0 );
      const SmSurfParamType aeDir[4]   = { SM_SP_V, SM_SP_U, SM_SP_V, SM_SP_U };
      const double          adAt[4]    = { 0.0, 1.0, 1.0, 0.0 };
      const SmOrientType    aeOrient[4] = { SM_OT_SAME, SM_OT_SAME, SM_OT_OPPOSITE, SM_OT_OPPOSITE };
      const SmPoint2d asFrom[4] = { SmPoint2d( 0, 0 ), SmPoint2d( 1, 0 ), SmPoint2d( 0, 1 ), SmPoint2d( 0, 0 ) };
      const SmPoint2d asTo[4]   = { SmPoint2d( 1, 0 ), SmPoint2d( 1, 1 ), SmPoint2d( 1, 1 ), SmPoint2d( 0, 1 ) };
      // which logical side each loop curve is: untransposed v=0 c0, u=1 d1, v=1 c1, u=0 d0
      const int aiSide[4]  = { 0, 3, 1, 2 };
      const int aiSideT[4] = { 2, 1, 3, 0 };
      for (int k=0; k<4; k++)
        {
          SmBSplineCurve * pIso = NULL;
          double dGap = 0.0, dDummy = 0.0;
          SER( pPatch->CreateIsoParametricCurve( m_crContext, aeDir[k], adAt[k], dDummy, pIso, &sDom, &dGap ));
          NER( pIso );
          SmBSplineCurve * pUV = NULL;
          SER( SmBSplineCurve::CreateLineSegment( m_crContext, 2, SmPoint3d( asFrom[k].x, asFrom[k].y, 0.0 ), SmPoint3d( asTo[k].x, asTo[k].y, 0.0 ), pUV ));
          m_vMixed3d[f].Add( pIso );
          m_vMixedUV[f].Add( pUV );
          m_vMixedOrients[f].Add( aeOrient[k] );
          SmCurve * pCopy = NULL;
          SER( pIso->Copy( m_crContext, pCopy ));
          const int iSide = bTurn ? aiSideT[k] : aiSide[k];
          apSide[f][iSide] = SM_CAST_PTR( SmBSplineCurve, pCopy );
          asCleanSide[f][iSide].SetObj( pCopy );
        }
    }

  // the edges' curves, from the faces' isos: X-A L c0, A-P L d1, P-V L c1 reversed, V-X L d0
  // reversed, V-Q U d1, Q-B U c1 reversed, X-B U d0
  const int aiFace[7] = { 0, 0, 0, 0, 1, 1, 1 };
  const int aiSide[7] = { 0, 3, 1, 2, 3, 1, 2 };
  const SmBoolean abRev[7] = { FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE };
  for (int k=0; k<7; k++)
    {
      SmCurve * pC = NULL;
      SER( apSide[aiFace[k]][aiSide[k]]->Copy( m_crContext, pC ));
      SmBSplineCurve * pCurve = SM_CAST_PTR( SmBSplineCurve, pC );
      if ( abRev[k] )
        {
          SmExtent1d sNew;
          SER( pCurve->ReverseParameterization( pCurve->GetNaturalInterval(), sNew ));
        }
      SmFilletEdge * pEdge = m_vEdges[k];
      pEdge->SetCurve( pCurve, TRUE );
      pCurve->SetOwner( pEdge );
    }
  // the trimlines' uv curves on the blends (X-A on B1, X-B on B2, both from X)
  for (int ii=0; ii<2; ii++)
    {
      SmBSplineCurve * pUV = NULL;
      SER( sm_MxToUV( m_crContext, asPc[ii], pUV ));
      SmFilletEdge * pEdge = m_vEdges[ii == 0 ? 0 : 6];
      SmFilletEdgeuse * pPrimEU = (SmFilletEdgeuse*)pEdge->GetPrimaryEdgeuse();
      pPrimEU->SetUVCurve( pUV );
      apFG[ii]->AddSideFilletEdgeuse( pPrimEU );
    }
  SER( SmFilletCorner::CalcCornerEdgeGeom() );
  return SM_SUCCESS;

} // end SmBlendCorner::CalcMixedEdgeGeom

void SmBlendCorner::DeleteMixedSurfaces()
{
  for (ULONG lFace=0; lFace<2; lFace++)
    {
      for (ULONG ii=0; ii<m_vMixed3d[lFace].GetSize(); ii++) { delete m_vMixed3d[lFace][ii]; }
      for (ULONG ii=0; ii<m_vMixedUV[lFace].GetSize(); ii++) { delete m_vMixedUV[lFace][ii]; }
      m_vMixed3d[lFace].ReSet();
      m_vMixedUV[lFace].ReSet();
      m_vMixedOrients[lFace].ReSet();
      if ( m_apMixedSurface[lFace] != NULL ) { delete m_apMixedSurface[lFace]; m_apMixedSurface[lFace] = NULL; }
    }
}

/*******************************************************************//**
PURPOSE: The Mixed corner's two faces, each its patch bounded by its
    natural loop. Their edges and the blends' and faces' meet in the stitch
    and the merge.
***********************************************************************/
SmStatus SmBlendCorner::CreateMixedPatches()
{
  SmBrep * pFilletBrep = m_pExecutive->GetFilletBrep();
  NER( pFilletBrep );
  for (ULONG lFace=0; lFace<2; lFace++)
    {
      NER( m_apMixedSurface[lFace] );
      SmTArray<ULONG> sCurveLoops;
      sCurveLoops.Add( m_vMixed3d[lFace].GetSize() );
      SmTArray<SmPoint3d> sLoopPoints;
      SmRegion  * pNewRegion = NULL;
      SmShell   * pNewShell  = NULL;
      SmFace    * pNewFace   = NULL;
      SmSurface * pSurface   = m_apMixedSurface[lFace];
      SmTArray<SmCurve*>        s3d( m_vMixed3d[lFace] );
      SmTArray<SmBSplineCurve*> sUV( m_vMixedUV[lFace] );
      m_apMixedSurface[lFace] = NULL;
      m_vMixed3d[lFace].ReSet();
      m_vMixedUV[lFace].ReSet();
      SER( pFilletBrep->MakeFaceWithCurves( NULL, sCurveLoops, &s3d, &sUV, m_vMixedOrients[lFace], sLoopPoints, pSurface,
                                            pSurface->GetNaturalUVDomain(), SM_OT_SAME, pNewRegion, pNewShell, pNewFace ));
      NER( pNewFace );
      m_vSurfaces.Add( pSurface );
    }
  return SM_SUCCESS;

} // end SmBlendCorner::CreateMixedPatches
