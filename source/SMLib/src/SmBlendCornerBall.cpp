// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

/**********************************************************************//**
* FILE NAME --- SmBlendCornerBall.cpp
* PURPOSE: SolidCraft change. The blend corner's Ball event
*   (SmBlendCorner.h): three constant rounds of one radius r and one
*   convexity at a vertex of three faces, each a plane or a natural quadric
*   (a sphere, a cylinder or a cone).
*
* The ball of radius r touching the three faces F_i is centred at C, the
* point at distance r from each on the side the rounds' balls roll: where
* the three faces' offsets by r meet (the triple offset point). It touches
* F_i at P_i = C - s_i r n_i, n_i the offset's normal at C (the face's at
* P_i). The round between F_i and F_j ends on its section through C, the
* circle of its ball there: its ball's centre runs square to both contact
* normals, so the section lies in the plane of n_i and n_j through C,
* centred at C, radius r, from P_i to P_j, whatever the faces. So all three
* end sections are great circles of the one sphere about C, and each round
* is tangent to that sphere along its end section (both are tangent to the
* ball there). The corner's face is the sphere's triangle P_0 P_1 P_2
* between the three sections: exact, G1 to every round.
*
* C itself, where two faces or more are planes: the line their offsets share
* meets the third's offset, a plane or a quadric (a quadratic along the
* line), in closed form; else Newton's method on the offsets' implicit equations,
* the root taken only where an interval (Krawczyk) test proves it the one
* root within rounding of the part's size of the rounds' own ball.
*
* The face is one rational patch of degree (4, 2) exactly on the sphere
* whose natural boundary is the triangle (sm_BallTrianglePatch), at every
* corner, square or oblique: its loop is that natural boundary, the pole at
* one contact left out, the three sections.
**********************************************************************/

#include "StdAfx.h"

#ifndef __SMBLENDCORNER_H__
#include <SmBlendCorner.h>
#endif

#ifndef __SMFILLETSOLVER_H__
#include <SmFilletSolver.h>
#endif

#include <SmAxis2Placement.h>
#include <SmFilletExact.h>
#include <cfloat>
#include <SmBallPatch.h>
#include <SmOffsetSurface.h>

namespace
{

// The faces of a corner ring of three edgeuses: face i between edgeuses i - 1 and i (as the
// corner's vertex i). NULL entries where two edgeuses share no face.
void sm_BallRingFaces(const SmTArray<SmEdgeuse*> & crAllEUs, SmFace * apFace[3])
{
  for (ULONG ii=0; ii<3; ii++)
    { apFace[ii] = sm_GetFaceByTwoEdgeuses( crAllEUs[(ii+2)%3], crAllEUs[ii], NULL ); }
}

// The quadric (a plane or a natural quadric) a ring face lies on; SM_BQ_NONE where it lies on none.
SmStatus sm_BallFaceQuadric(const SmFace * pFace, double dTol, SmBlendQuadric & rQuadric)
{
  rQuadric = SmBlendQuadric();
  NER( pFace );
  NER( pFace->GetSurface() );
  SER( SmBlendQuadric::FromSurface( pFace->GetSurface(), pFace->GetUVDomain().GetMid(), dTol, rQuadric ));
  if (   rQuadric.m_eKind != SmBlendQuadric::SM_BQ_PLANE && rQuadric.m_eKind != SmBlendQuadric::SM_BQ_SPHERE
      && rQuadric.m_eKind != SmBlendQuadric::SM_BQ_CYLINDER && rQuadric.m_eKind != SmBlendQuadric::SM_BQ_CONE )
    { rQuadric = SmBlendQuadric(); }
  return SM_SUCCESS;
}

// A quadric's implicit form f(x) = x.A.x + 2 b.x + c about the origin crAt (x = point - crAt).
struct SmBallImplicit
{
  double a[3][3], b[3], c;

  void Set(const SmBlendQuadric & crQuadric, const SmPoint3d & crAt)
  {
    SmBlendQuadric sLocal( crQuadric );
    sLocal.m_sPoint = SmPoint3d( crQuadric.m_sPoint.x - crAt.x, crQuadric.m_sPoint.y - crAt.y, crQuadric.m_sPoint.z - crAt.z );
    sLocal.Implicit( a, b, c );
  }
  double Value(const double x[3], double * pdBound = NULL) const
  {
    double d = c, dAbs = smos_Fabs( c );
    for (int i=0; i<3; i++)
      {
        d += 2.0 * b[i] * x[i];
        dAbs += smos_Fabs( 2.0 * b[i] * x[i] );
        for (int j=0; j<3; j++) { d += x[i] * a[i][j] * x[j]; dAbs += smos_Fabs( x[i] * a[i][j] * x[j] ); }
      }
    if ( pdBound != NULL ) { *pdBound = dAbs; }
    return d;
  }
  void Gradient(const double x[3], double g[3]) const
  {
    for (int i=0; i<3; i++)
      { g[i] = 2.0 * ( a[i][0] * x[0] + a[i][1] * x[1] + a[i][2] * x[2] + b[i] ); }
  }
};

SmBoolean sm_BallInverse(const double m[3][3], double inv[3][3])
{
  const double dDet = m[0][0] * ( m[1][1] * m[2][2] - m[1][2] * m[2][1] ) - m[0][1] * ( m[1][0] * m[2][2] - m[1][2] * m[2][0] )
                    + m[0][2] * ( m[1][0] * m[2][1] - m[1][1] * m[2][0] );
  double dScale = 0.0;
  for (int i=0; i<3; i++) for (int j=0; j<3; j++) { dScale = smos_Max( dScale, smos_Fabs( m[i][j] )); }
  if ( !( smos_Fabs( dDet ) > SM_EFF_ZERO_SQRT * dScale * dScale * dScale ))
    { return FALSE; }
  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
      {
        const int i1 = (j+1)%3, i2 = (j+2)%3, j1 = (i+1)%3, j2 = (i+2)%3;
        inv[i][j] = ( m[i1][j1] * m[i2][j2] - m[i1][j2] * m[i2][j1] ) / dDet;
      }
  return TRUE;
}

/*******************************************************************//**
PURPOSE: The one point where three quadrics' implicit equations vanish near
    the origin (their offsets about the rounds' ball): Newton's method from
    the origin, then Krawczyk's test on the box of half-width dHalf about
    the iterate, which proves it holds exactly one root (the gradients are
    affine, so the Jacobian over the box is the iterate's plus 2 A times the
    box, exactly). Rounding is bounded by the sum of the terms' magnitudes.
    rbProved FALSE where Newton fails or the test does not prove the box.
***********************************************************************/
void sm_BallCertifiedRoot(const SmBallImplicit aF[3], double dHalf, double x[3], SmBoolean & rbProved)
{
  rbProved = FALSE;
  x[0] = x[1] = x[2] = 0.0;
  double J[3][3], Y[3][3], f[3];
  for (int it=0; it<50; it++)
    {
      for (int i=0; i<3; i++) { f[i] = aF[i].Value( x ); aF[i].Gradient( x, J[i] ); }
      if ( !sm_BallInverse( J, Y ))
        { return; }
      double dStep = 0.0;
      for (int i=0; i<3; i++)
        {
          const double d = Y[i][0] * f[0] + Y[i][1] * f[1] + Y[i][2] * f[2];
          x[i] -= d;
          dStep = smos_Max( dStep, smos_Fabs( d ));
        }
      if ( dStep <= 0.01 * dHalf && it > 0 )
        { break; }
    }
  // Krawczyk: K = x - Y f(x) + (I - Y J(X)) (X - x) inside X proves exactly one root in X. With
  // J(X) = J(x) + dJ, |dJ_mk| <= 2 dHalf sum_l |A_m,kl|, component i of K - x is within
  // |(Y f)_i| + dHalf sum_k ( |(I - Y J(x))_ik| + sum_m |Y_im| |dJ_mk| ), each term padded for
  // rounding (f by its terms' magnitudes).
  double adBound[3];
  for (int i=0; i<3; i++) { f[i] = aF[i].Value( x, &adBound[i] ); aF[i].Gradient( x, J[i] ); }
  if ( !sm_BallInverse( J, Y ))
    { return; }
  const double dEps = 16.0 * DBL_EPSILON;
  double adSpread[3][3];
  for (int m=0; m<3; m++)
    for (int k=0; k<3; k++)
      { adSpread[m][k] = 2.0 * dHalf * ( smos_Fabs( aF[m].a[k][0] ) + smos_Fabs( aF[m].a[k][1] ) + smos_Fabs( aF[m].a[k][2] )); }
  for (int i=0; i<3; i++)
    {
      double dRad = 0.0;
      for (int k=0; k<3; k++)
        {
          dRad += smos_Fabs( Y[i][k] ) * ( smos_Fabs( f[k] ) + dEps * adBound[k] );
          double dIYJ = ( i == k ) ? 1.0 : 0.0, dScale = 1.0;
          for (int m=0; m<3; m++)
            {
              dIYJ -= Y[i][m] * J[m][k];
              dScale += smos_Fabs( Y[i][m] * J[m][k] );
              dRad += smos_Fabs( Y[i][m] ) * adSpread[m][k] * dHalf;
            }
          dRad += ( smos_Fabs( dIYJ ) + dEps * dScale ) * dHalf;
        }
      if ( !( dRad < dHalf ))
        { return; }
    }
  rbProved = TRUE;
}

// A polynomial in (s, t) of degree (m, n) in the Bernstein basis: c[i][j] the coefficient of
// B^m_i(s) B^n_j(t). Degrees up to (4, 2), the ball patch's.
struct SmBallPoly
{
  ULONG  m = 0, n = 0;
  double c[5][3] = {};
};

double sm_Choose(ULONG n, ULONG k)
{
  double d = 1.0;
  for (ULONG i=1; i<=k; i++) { d = d * double( n - k + i ) / double( i ); }
  return d;
}

// The product of two polynomials in the Bernstein basis, of the sum of their degrees.
SmBallPoly sm_BallTimes(const SmBallPoly & crA, const SmBallPoly & crB)
{
  SmBallPoly sP;
  sP.m = crA.m + crB.m;
  sP.n = crA.n + crB.n;
  for (ULONG i=0; i<=sP.m; i++)
    for (ULONG j=0; j<=sP.n; j++)
      {
        double d = 0.0;
        for (ULONG a=0; a<=crA.m; a++)
          for (ULONG b=0; b<=crA.n; b++)
            {
              if ( a > i || b > j || i - a > crB.m || j - b > crB.n )
                { continue; }
              d += crA.c[a][b] * crB.c[i-a][j-b] * sm_Choose( crA.m, a ) * sm_Choose( crB.m, i - a ) * sm_Choose( crA.n, b ) * sm_Choose( crB.n, j - b );
            }
        sP.c[i][j] = d / ( sm_Choose( sP.m, i ) * sm_Choose( sP.n, j ));
      }
  return sP;
}

// dA crA + dB crB, of one degree.
SmBallPoly sm_BallSum(double dA, const SmBallPoly & crA, double dB, const SmBallPoly & crB)
{
  SmBallPoly sP = crA;
  for (ULONG i=0; i<=sP.m; i++)
    for (ULONG j=0; j<=sP.n; j++)
      { sP.c[i][j] = dA * crA.c[i][j] + dB * crB.c[i][j]; }
  return sP;
}

/*******************************************************************//**
PURPOSE: The ball's triangle between the contacts as one rational patch on
    the sphere, where no contact is a pole with the opposite side its
    equator (see the file's head).

NOTES: The sphere's stereographic projection from the antipode of the
    contact a_k (the pole, at the plane's origin) takes the sides through it
    to segments of lines through the origin, from 0 to A and to B, and the
    side between the others to a circle's arc from A to B: a rational
    quadratic c(s). The plane's sector t c(s), s, t in [0, 1], taken back to
    the sphere by the inverse projection (of degree 2), is a rational patch
    of degree (4, 2) on the sphere: exactly the triangle, t = 0 the pole,
    s = 0 and s = 1 its sides through the pole, t = 1 the third.
    Homogeneously, with (U, V, T) = (t X(s), t Y(s), W(s)):
      S = C + r (2 U T e1 + 2 V T e2 + (U^2 + V^2 - T^2) e3) / (U^2 + V^2 + T^2),
    e3 the projection's centre, -a_k. With bReverseS, s runs the other way
    (and the normal is turned).
***********************************************************************/
SmStatus sm_BallTrianglePatch
  (const SmContext  & crContext,
   const SmPoint3d  & crCentre,
   double             dRadius,
   const SmVector3d   asDir[3],   // unit: the contacts less the centre, over the radius
   ULONG              lPole,
   SmBoolean          bReverseS,
   SmBSplineSurface *& rpPatch,
   ULONG            & rlFrom,     // the contact at s = 0, t = 1
   ULONG            & rlTo)       // and at s = 1
{
  rpPatch = NULL;
  rlFrom = (lPole+1)%3;
  rlTo   = (lPole+2)%3;
  SmVector3d sE3 = -asDir[lPole];
  SmVector3d sE1 = asDir[rlFrom] - sE3 * asDir[rlFrom].Dot( sE3 );
  SER( sE1.Unitize() );
  SmVector3d sE2 = sE3 * sE1;
  // the projection: p -> (p.e1, p.e2) / (1 - p.e3)
  double adA[2], adB[2], adM[2];
  SmVector3d sMid = asDir[rlFrom] + asDir[rlTo];
  SER( sMid.Unitize() );
  const SmVector3d * apP[3] = { &asDir[rlFrom], &asDir[rlTo], &sMid };
  double * apOut[3] = { adA, adB, adM };
  for (ULONG ii=0; ii<3; ii++)
    {
      const double dDen = 1.0 - apP[ii]->Dot( sE3 );
      if ( dDen <= SM_EFF_ZERO )
        { SER( SM_ERR ); }
      apOut[ii][0] = apP[ii]->Dot( sE1 ) / dDen;
      apOut[ii][1] = apP[ii]->Dot( sE2 ) / dDen;
    }
  // the third side's great circle (normal m) projects to the circle about -(m1, m2) / m3,
  // of radius^2 1 + (m1^2 + m2^2) / m3^2
  SmVector3d sM = asDir[rlFrom] * asDir[rlTo];
  SER( sM.Unitize() );
  const double dM3 = sM.Dot( sE3 );
  if ( smos_Fabs( dM3 ) <= SM_EFF_ZERO )
    { SER( SM_ERR ); } // the third side runs through the pole: no triangle
  const double adO[2] = { -sM.Dot( sE1 ) / dM3, -sM.Dot( sE2 ) / dM3 };
  const double dR = sqrt( 1.0 + ( adO[0] * adO[0] + adO[1] * adO[1] ));
  double adUA[2] = { ( adA[0] - adO[0] ) / dR, ( adA[1] - adO[1] ) / dR };
  double adUB[2] = { ( adB[0] - adO[0] ) / dR, ( adB[1] - adO[1] ) / dR };
  const double dCos = adUA[0] * adUB[0] + adUA[1] * adUB[1];
  // the arc from A to B through the middle's image, less than a half turn: its middle control
  // point where the tangents at A and B cross, its weight the cosine of half its angle
  if (   1.0 + dCos <= SM_EFF_ZERO_SQRT
      || ( adM[0] - adO[0] ) * ( adUA[0] + adUB[0] ) + ( adM[1] - adO[1] ) * ( adUA[1] + adUB[1] ) <= 0.0 )
    { SER( SM_ERR ); }
  const double dW = sqrt( 0.5 * ( 1.0 + dCos ));
  const double adD[2] = { adO[0] + dR * ( adUA[0] + adUB[0] ) / ( 1.0 + dCos ), adO[1] + dR * ( adUA[1] + adUB[1] ) / ( 1.0 + dCos ) };

  // (U, V, T) of degree (2, 1)
  SmBallPoly sU, sV, sT;
  sU.m = sV.m = sT.m = 2;
  sU.n = sV.n = sT.n = 1;
  const double adX[3] = { adA[0], dW * adD[0], adB[0] };
  const double adY[3] = { adA[1], dW * adD[1], adB[1] };
  const double adWt[3] = { 1.0, dW, 1.0 };
  for (ULONG i=0; i<3; i++)
    {
      sU.c[i][1] = adX[i];
      sV.c[i][1] = adY[i];
      sT.c[i][0] = sT.c[i][1] = adWt[i];
    }
  const SmBallPoly sUT = sm_BallTimes( sU, sT ), sVT = sm_BallTimes( sV, sT );
  const SmBallPoly sSq = sm_BallSum( 1.0, sm_BallTimes( sU, sU ), 1.0, sm_BallTimes( sV, sV ));
  const SmBallPoly sTT = sm_BallTimes( sT, sT );
  const SmBallPoly sH3 = sm_BallSum( 1.0, sSq, -1.0, sTT ), sH4 = sm_BallSum( 1.0, sSq, 1.0, sTT );

  // control points, s (u) outer, t (v) inner
  SmTArray<SmPoint3d> sPoints;
  SmTArray<double>    sWeights;
  for (ULONG ii=0; ii<=4; ii++)
    {
      const ULONG i = bReverseS ? 4 - ii : ii;
      for (ULONG j=0; j<=2; j++)
        {
          const double dWij = sH4.c[i][j];
          if ( dWij <= SM_EFF_ZERO )
            { SER( SM_ERR ); }
          const SmVector3d sOff = ( sE1 * ( 2.0 * sUT.c[i][j] ) + sE2 * ( 2.0 * sVT.c[i][j] ) + sE3 * sH3.c[i][j] ) * ( dRadius / dWij );
          sPoints.Add( crCentre + sOff );
          sWeights.Add( dWij );
        }
    }
  if ( bReverseS )
    { ULONG lSwap = rlFrom; rlFrom = rlTo; rlTo = lSwap; }
  SmTArray<ULONG> sUMults, sVMults;
  sUMults.Add( 5 ); sUMults.Add( 5 );
  sVMults.Add( 3 ); sVMults.Add( 3 );
  SmTArray<double> sKnots;
  sKnots.Add( 0.0 ); sKnots.Add( 1.0 );
  SER( SmBSplineSurface::CreateCanonical( crContext, 4, 2, sPoints, SM_SF_UNSPECIFIED, sUMults, sVMults, sKnots, sKnots,
                                          SM_KT_PIECEWISE_BEZIER_KNOTS, &sWeights, NULL, rpPatch ));
  NER( rpPatch );
  return SM_SUCCESS;

} // end sm_BallTrianglePatch

} // namespace

/*******************************************************************//**
PURPOSE: Whether a filleted vertex is a Ball: every edge there blended (three
    of them), each a constant circular round of the same radius, all convex
    or all concave, the three faces planes meeting at a point, and the
    caller not setting the vertex back (a setback corner).
***********************************************************************/
SmStatus SmBlendCorner::ClassifyBall
  (const SmVertex             * cpVertex,
   const SmTArray<SmEdgeuse*> & crSolverEUs,
   const SmTArray<SmEdgeuse*> & crAllEUs,
   SmFilletExecutive          * pExec,
   SmBlendCornerEvent         & reEvent)
{
  reEvent = SM_BCE_NONE;
  NER( cpVertex );
  NER( pExec );
  if ( crSolverEUs.GetSize() != 3 || crAllEUs.GetSize() != 3 || pExec->IsSetBackVertex( cpVertex ))
    { return SM_SUCCESS; }

  // constant circular rounds of one radius, of one convexity
  double    dRadius = -1.0;
  double    dTol    = SM_BIG_DOUBLE;
  SmBoolean bConcave0 = FALSE;
  for (ULONG ii=0; ii<3; ii++)
    {
      SmEdgeuse * pEU = crAllEUs[ii];
      NER( pEU );
      SmFilletSolver * pFS = pExec->GetFilletSolverOfEdgeuse( pEU );
      if ( pFS == NULL )
        { pFS = pExec->GetFilletSolverOfEdgeuse( pEU->GetMate() ); }
      if ( pFS == NULL )
        { return SM_SUCCESS; }
      SmFilletSurfaceGenerator * pFSG = pFS->GetFilletSurfaceGenerator();
      const SmOffsetSurface * pOff0 = pFS->GetSurface(0);
      const SmOffsetSurface * pOff1 = pFS->GetSurface(1);
      if (   pFSG == NULL || pOff0 == NULL || pOff1 == NULL || pFSG->GetFilletSurfaceGeneratorType() != SM_FSG_CIRCULAR
          || pFS->GetSolverType() != SM_FS_CONST_RADIUS || pFS->OffsetRadiiCanChange() )
        { return SM_SUCCESS; }
      dTol = smos_Min( dTol, pFS->GetThisApproxTol3d() );
      const double dR0 = smos_Fabs( pOff0->GetOffsetDistance() ), dR1 = smos_Fabs( pOff1->GetOffsetDistance() );
      if ( smos_Fabs( dR0 - dR1 ) > dTol || ( dRadius >= 0.0 && smos_Fabs( dR0 - dRadius ) > dTol ))
        { return SM_SUCCESS; }
      dRadius = dR0;
      SmBoolean bConcave = FALSE;
      SER( SmFilletSolver::IsConcaveEdge( *pEU->GetEdge(), bConcave ));
      if ( ii == 0 )
        { bConcave0 = bConcave; }
      else if ( bConcave != bConcave0 )
        { return SM_SUCCESS; }
    }
  if ( !( dRadius > dTol ))
    { return SM_SUCCESS; }

  // three faces, each a plane or a natural quadric (a cylinder or sphere wider than the ball),
  // their normals at the vertex independent
  SmFace * apFace[3];
  sm_BallRingFaces( crAllEUs, apFace );
  SmVector3d asNormal[3];
  for (ULONG ii=0; ii<3; ii++)
    {
      if ( apFace[ii] == NULL )
        { return SM_SUCCESS; }
      SmBlendQuadric sQuadric;
      SER( sm_BallFaceQuadric( apFace[ii], dTol, sQuadric ));
      if ( sQuadric.m_eKind == SmBlendQuadric::SM_BQ_NONE )
        { return SM_SUCCESS; }
      asNormal[ii] = sQuadric.Normal( cpVertex->GetPoint() );
      if ( asNormal[ii].Unitize( FALSE ) != SM_SUCCESS )
        { return SM_SUCCESS; } // a cone's apex
      if (   ( sQuadric.m_eKind == SmBlendQuadric::SM_BQ_SPHERE || sQuadric.m_eKind == SmBlendQuadric::SM_BQ_CYLINDER )
          && !( sQuadric.m_dRadius > dRadius + dTol ))
        { return SM_SUCCESS; } // a ball as wide as the cylinder or sphere: its offset is its axis or centre
    }
  if ( apFace[0] == apFace[1] || apFace[1] == apFace[2] || apFace[2] == apFace[0] || smos_Fabs( asNormal[0].Dot( asNormal[1] * asNormal[2] )) <= SM_EFF_ZERO_SQRT )
    { return SM_SUCCESS; }
  reEvent = SM_BCE_BALL;
  return SM_SUCCESS;

} // end SmBlendCorner::ClassifyBall

/*******************************************************************//**
PURPOSE: The Ball's topology: vertex i where the rails of rounds i - 1 and i
    cross on their common face (SM_FV_RAIL_X_RAIL, the ball's contact with
    that face), and edge i from vertex i to vertex i + 1, round i's end
    section (SM_FE_CROSS_SECTION).
***********************************************************************/
SmStatus SmBlendCorner::MakeBallTopology()
{
  NER( m_pAllEUs );
  if ( !IsBall() || m_pAllEUs->GetSize() != 3 )
    { SER( SM_ERR ); }
  SmFace * apFace[3];
  sm_BallRingFaces( *m_pAllEUs, apFace );
  SmFilletGeom * apFG[3];
  SmFilletSolver * apFS[3];
  for (ULONG ii=0; ii<3; ii++)
    {
      SmEdgeuse * pEU = (*m_pAllEUs)[ii];
      apFS[ii] = GetFilletExecutive()->GetFilletSolverOfEdgeuse( pEU );
      if ( apFS[ii] == NULL )
        { apFS[ii] = GetFilletExecutive()->GetFilletSolverOfEdgeuse( pEU->GetMate() ); }
      NER( apFS[ii] );
      apFG[ii] = apFS[ii]->GetFirstFilletGeom();
      NER( apFG[ii] );
      NER( apFace[ii] );
    }

  // the contacts: rounds i - 1 and i's rails on face i
  SmFilletEdgeuse * pNewEU = NULL;
  for (ULONG ii=0; ii<3; ii++)
    {
      SmFilletVertex * pFV = new (m_pPseudoBrep) SmFilletVertex( this );
      NER( pFV );
      pFV->SetPointClass( SM_PC_FACE, apFace[ii] );
      pFV->SetFilletVertexType( SM_FV_RAIL_X_RAIL );
      const ULONG lPrev = (ii+2)%3;
      SER( apFG[lPrev]->MakeRailEdgeuse( m_pPseudoBrep, apFS[lPrev]->FindIndexOfRailOnFace( apFace[ii] ), pFV, pNewEU ));
      SER( apFG[ii]->MakeRailEdgeuse( m_pPseudoBrep, apFS[ii]->FindIndexOfRailOnFace( apFace[ii] ), pFV, pNewEU ));
      m_vVertices.Add( pFV );
    }

  // each round ends on its section through the ball's centre; the corner's face is a patch of
  // its own
  m_bBlending = TRUE;
  m_bBevel    = FALSE;
  for (ULONG ii=0; ii<3; ii++)
    {
      SmFilletEdge * pEdge = NULL;
      SER( SmFilletEdge::MakeFilletEdge( m_pPseudoBrep, m_vVertices[ii], m_vVertices[(ii+1)%3], pEdge, this ));
      NER( pEdge );
      pEdge->SetFilletEdgeType( SM_FE_CROSS_SECTION );
      SmFilletEdgeuse * pPrimEU = (SmFilletEdgeuse*)pEdge->GetPrimaryEdgeuse();
      NER( pPrimEU );
      pPrimEU->SetFilletGeom( apFG[ii] );
      m_vEdges.Add( pEdge );
    }
  return SM_SUCCESS;

} // end SmBlendCorner::MakeBallTopology

/*******************************************************************//**
PURPOSE: The Ball's centre and contacts (see the file's head), from the
    three faces and the processed contacts. The centre is the faces' triple
    offset point, solved from them; the rounds' balls at the contacts are
    checked to be it.
***********************************************************************/
SmStatus SmBlendCorner::BallCentre(SmPoint3d & rCentre, double & rdRadius, SmPoint3d asContact[3]) const
{
  NER( m_pAllEUs );
  if ( m_vVertices.GetSize() != 3 || m_vEdges.GetSize() != 3 || m_pAllEUs->GetSize() != 3 )
    { SER( SM_ERR ); }
  const double dTol = GetThisApproxTol3d();

  // a ball of the rounds at the contacts: round 0's offset surfaces at vertex 0
  SmFilletVertex * pV0 = m_vVertices[0];
  NER( pV0 );
  SmFilletEdgeuse * pSection0 = (SmFilletEdgeuse*)m_vEdges[0]->GetPrimaryEdgeuse();
  NER( pSection0 );
  SmFilletGeom * pFG0 = pSection0->GetFilletGeom();
  NER( pFG0 );
  SmFilletVertexuse * pVU = pV0->GetVUAtRailEnd( pFG0 );
  NER( pVU );
  SmTsectPnt & rTsect = pVU->GetTsectPnt();
  SmPoint3d sOn[2];
  for (ULONG k=0; k<2; k++)
    {
      SmOffsetSurface * pOffset = pFG0->GetOffsetSurface( k );
      NER( pOffset );
      SER( pOffset->EvaluatePoint( rTsect.UVPos( k ), sOn[k] ));
    }
  if ( sOn[0].DistanceBetween( sOn[1] ) > dTol )
    { SER( SM_ERR ); } // the section's offset points are its ball's centre
  const SmPoint3d sBall = sOn[0];
  rdRadius = smos_Fabs( pFG0->GetOffsetSurface( 0 )->GetOffsetDistance() );

  // each face's quadric offset by r on the ball's side (the side of it the rounds' ball is on,
  // seen from vertex i, the contact on it)
  SmFace * apFace[3];
  sm_BallRingFaces( *m_pAllEUs, apFace );
  SmBlendQuadric asOffset[3];
  double         adSide[3];
  ULONG          lPlanes = 0;
  for (ULONG ii=0; ii<3; ii++)
    {
      SmBlendQuadric sFace;
      SER( sm_BallFaceQuadric( apFace[ii], dTol, sFace ));
      if ( sFace.m_eKind == SmBlendQuadric::SM_BQ_NONE )
        { SER( SM_ERR ); } // ClassifyBall took the faces to be planes or quadrics
      const SmPoint3d sVertex = m_vVertices[ii]->GetPoint();
      adSide[ii] = ( sFace.Normal( sVertex ).Dot( sBall - sVertex ) >= 0.0 ) ? 1.0 : -1.0;
      asOffset[ii] = sFace.Offset( adSide[ii] * rdRadius, dTol );
      if ( asOffset[ii].m_eKind != sFace.m_eKind || asOffset[ii].m_dSide != sFace.m_dSide )
        { SER( SM_ERR ); } // ClassifyBall took the quadrics wider than the ball
      if ( asOffset[ii].m_eKind == SmBlendQuadric::SM_BQ_PLANE )
        { lPlanes++; }
    }

  // C: where the offsets meet, about the rounds' ball (each quadric's implicit form there)
  SmBallImplicit asF[3];
  ULONG alOrder[3] = { 0, 1, 2 }; // the planes first
  for (ULONG ii=0, lAt=0; ii<3; ii++)
    { if ( asOffset[ii].m_eKind == SmBlendQuadric::SM_BQ_PLANE ) { alOrder[lAt++] = ii; } }
  for (ULONG ii=0, lAt=lPlanes; ii<3; ii++)
    { if ( asOffset[ii].m_eKind != SmBlendQuadric::SM_BQ_PLANE ) { alOrder[lAt++] = ii; } }
  for (ULONG ii=0; ii<3; ii++)
    { asF[ii].Set( asOffset[alOrder[ii]], sBall ); }
  double x[3] = { 0.0, 0.0, 0.0 };
  if ( lPlanes >= 2 )
    {
      // the two offset planes' line, n_0 . x = -c_0 / 2 and n_1 . x = -c_1 / 2 (2 b = n)
      const SmVector3d sN0( 2.0 * asF[0].b[0], 2.0 * asF[0].b[1], 2.0 * asF[0].b[2] ), sN1( 2.0 * asF[1].b[0], 2.0 * asF[1].b[1], 2.0 * asF[1].b[2] );
      const SmVector3d sD = sN0 * sN1;
      const double dDD = sD.Dot( sD );
      if ( dDD <= SM_EFF_ZERO )
        { SER( SM_ERR ); }
      // its point nearest the origin: (h_0 (n_1 x d) + h_1 (d x n_0)) / |d|^2, h_i = -c_i
      const SmVector3d sP0 = ( ( sN1 * sD ) * ( -asF[0].c ) + ( sD * sN0 ) * ( -asF[1].c )) * ( 1.0 / dDD );
      const double p0[3] = { sP0.x, sP0.y, sP0.z }, d[3] = { sD.x, sD.y, sD.z };
      if ( lPlanes == 3 )
        {
          // the third plane: n_2 . (p0 + t d) = -c_2
          const double dDen = 2.0 * ( asF[2].b[0] * d[0] + asF[2].b[1] * d[1] + asF[2].b[2] * d[2] );
          if ( smos_Fabs( dDen ) <= SM_EFF_ZERO )
            { SER( SM_ERR ); }
          const double t = -asF[2].Value( p0 ) / dDen;
          for (int k=0; k<3; k++) { x[k] = p0[k] + t * d[k]; }
        }
      else
        {
          // the quadric along the line: q2 t^2 + 2 q1 t + q0 = 0, the root nearest the rounds' ball
          double dA = 0.0, dB = 0.0;
          for (int k=0; k<3; k++)
            {
              double dAd = 0.0, dAp = 0.0;
              for (int m=0; m<3; m++) { dAd += asF[2].a[k][m] * d[m]; dAp += asF[2].a[k][m] * p0[m]; }
              dA += d[k] * dAd;
              dB += d[k] * ( dAp + asF[2].b[k] );
            }
          const double dC = asF[2].Value( p0 );
          const double dDisc = dB * dB - dA * dC;
          if ( dDisc < 0.0 || ( smos_Fabs( dA ) <= SM_EFF_ZERO * dDD && smos_Fabs( dB ) <= SM_EFF_ZERO * dDD ))
            { SER( SM_ERR ); } // the line misses the offset quadric
          double adT[2];
          ULONG  lRoots = 0;
          if ( smos_Fabs( dA ) <= SM_EFF_ZERO * dDD )
            { adT[lRoots++] = -dC / ( 2.0 * dB ); }
          else
            {
              // without cancellation: q = -(b + sign(b) sqrt(disc)), t = q / a and c / q
              const double dQ = -( dB + ( dB >= 0.0 ? 1.0 : -1.0 ) * sqrt( dDisc ));
              adT[lRoots++] = dQ / dA;
              if ( dQ != 0.0 ) { adT[lRoots++] = dC / dQ; }
            }
          // the root nearest the origin (the rounds' ball): |p0 + t d|, p0 square to d
          double dBest = adT[0];
          for (ULONG k=1; k<lRoots; k++) { if ( smos_Fabs( adT[k] ) < smos_Fabs( dBest )) { dBest = adT[k]; } }
          for (int k=0; k<3; k++) { x[k] = p0[k] + dBest * d[k]; }
        }
    }
  else
    {
      // no closed form: the root certified within rounding of the part's size
      double dSize = 0.0;
      for (ULONG ii=0; ii<3; ii++)
        {
          SmExtent3d sFaceBox;
          SER( apFace[ii]->CalculateBoundingBox( sFaceBox ));
          dSize = smos_Max( dSize, ( sFaceBox.GetMax() - sFaceBox.GetMin() ).Length() );
        }
      SmBoolean bProved = FALSE;
      sm_BallCertifiedRoot( asF, 1e-12 * smos_Max( dSize, 1.0 ), x, bProved );
      if ( !bProved )
        { SER( SM_ERR ); }
    }
  rCentre = SmPoint3d( sBall.x + x[0], sBall.y + x[1], sBall.z + x[2] );
  if ( rCentre.DistanceBetween( sBall ) > dTol )
    { SER( SM_ERR ); } // the rounds' ball is the offsets' meeting point
  // the contacts: the centre moved back along the offset's normal there (the face's at the contact)
  for (ULONG ii=0; ii<3; ii++)
    {
      SmVector3d sN = asOffset[ii].Normal( rCentre );
      SER( sN.Unitize() );
      asContact[ii] = rCentre - sN * ( adSide[ii] * rdRadius );
      if ( asContact[ii].DistanceBetween( m_vVertices[ii]->GetPoint() ) > dTol )
        { SER( SM_ERR ); } // vertex i is the ball's contact with face i
    }
  return SM_SUCCESS;

} // end SmBlendCorner::BallCentre

/*******************************************************************//**
PURPOSE: The Ball's contacts as they are in closed form, not as the rails'
    crossings were solved: each vertex at its contact, and each round's
    section there (its vertexuse's tsect point, which the round's trace ends
    on) at the ball's centre on both its offset surfaces. The exact round
    (SmFilletSolver::CalcLineSpineFilletGeom) runs its axis and span to that
    section, so its end circle is the ball's.
***********************************************************************/
SmStatus SmBlendCorner::PlaceBallContacts()
{
  SmPoint3d sCentre, asContact[3];
  double    dRadius = 0.0;
  SER( BallCentre( sCentre, dRadius, asContact ));
  const double dTol = GetThisApproxTol3d();
  for (ULONG ii=0; ii<3; ii++)
    {
      SmFilletVertex * pFV = m_vVertices[ii];
      pFV->SetPoint( asContact[ii] );
      // every round's section there: the tsect point of each vertexuse of a round's rail
      SmTArray<SmVertexuse*> sVUs;
      pFV->GetVertexuses( sVUs );
      for (ULONG u=0; u<sVUs.GetSize(); u++)
        {
          SmFilletVertexuse * pVU = (SmFilletVertexuse*)sVUs[u];
          SmFilletEdgeuse   * pEU = pVU != NULL ? (SmFilletEdgeuse*)pVU->GetEdgeuse() : NULL;
          SmFilletGeom      * pFG = pEU != NULL ? pEU->GetFilletGeom() : NULL;
          if ( pFG == NULL || ((SmFilletEdge*)pEU->GetEdge())->GetFilletEdgeType() != SM_FE_RAIL )
            { continue; }
          SmTsectPnt & rTsect = pVU->GetTsectPnt();
          for (ULONG k=0; k<2; k++)
            {
              SmOffsetSurface * pOffset = pFG->GetOffsetSurface( k );
              NER( pOffset );
              SmBoolean  bFound = FALSE;
              SmSolution sSolution;
              SER( pOffset->LocalPointSolve( pOffset->GetNaturalUVDomain(), SM_SO_MINIMIZE, sCentre, rTsect.UVPos( k ), bFound, sSolution ));
              const SmPoint2d sUV( sSolution.m_vStart[0], sSolution.m_vStart[1] );
              SmPoint3d sAt;
              SER( pOffset->EvaluatePoint( sUV, sAt ));
              if ( !bFound || sAt.DistanceBetween( sCentre ) > dTol )
                { SER( SM_ERR ); }
              rTsect.UVPos( k ) = sUV;
            }
        }
    }
  return SM_SUCCESS;

} // end SmBlendCorner::PlaceBallContacts

SmStatus SmBlendCorner::BuildBallSurface()
{
  DeleteBallSurface();
  SmPoint3d sCentre, asContact[3];
  double    dRadius = 0.0;
  SER( BallCentre( sCentre, dRadius, asContact ));
  const double dTol = GetThisApproxTol3d();

  // The face's normal is its surface's (SmBrep::OrientTrimmedSurfaces turns a face the other way
  // by swapping its surface's parameters): out of the material, towards the centre where the
  // rounds add material (concave), away from it where they remove it.
  SmBoolean bAddsMaterial = FALSE;
  SER( SmFilletSolver::IsConcaveEdge( *(*m_pAllEUs)[0]->GetEdge(), bAddsMaterial ));

  // the face: one rational patch on the sphere, the pole at contact 0, its normal turned as the
  // face's
  const ULONG lPole = 0;
  ULONG lFrom = 0, lTo = 0;
  SmSurface * pSurface = NULL;
  SmVector3d asDir[3];
  for (ULONG ii=0; ii<3; ii++)
    {
      asDir[ii] = asContact[ii] - sCentre;
      SER( asDir[ii].Unitize() );
    }
  for (ULONG lTurn=0; lTurn<2 && pSurface == NULL; lTurn++)
    {
      SmBSplineSurface * pPatch = NULL;
      SER( sm_BallTrianglePatch( m_crContext, sCentre, dRadius, asDir, lPole, lTurn == 1 ? TRUE : FALSE, pPatch, lFrom, lTo ));
      SmObjDelete sClean( pPatch );
      const SmPoint2d sMid = pPatch->GetNaturalUVDomain().Evaluate( 0.5, 0.5 );
      SmPoint3d  sAt;
      SmVector3d sNormal;
      SER( pPatch->EvaluatePoint( sMid, sAt ));
      SER( pPatch->EvaluateNormal( sMid, TRUE, TRUE, sNormal ));
      if ( ( sNormal.Dot( sAt - sCentre ) < 0.0 ) == ( bAddsMaterial ? TRUE : FALSE ))
        {
          // the face's surface carries its sphere as its definition (its offset is exact)
          SmBallPatch * pBall = NULL;
          SER( SmBallPatch::Create( m_crContext, *pPatch, sCentre, dRadius, pBall ));
          m_pBallSurface = pSurface = pBall;
        }
    }
  NER( pSurface );

  // its corners are the contacts: the pole at v = vmin, the others at v = vmax
  const SmExtent2d sNurbDomain = pSurface->GetNaturalUVDomain();
  const SmPoint2d asUV[3] = { SmPoint2d( sNurbDomain.GetMin().x, sNurbDomain.GetMax().y ), SmPoint2d( sNurbDomain.GetMax().x, sNurbDomain.GetMax().y ),
                              SmPoint2d( sNurbDomain.GetMin().x, sNurbDomain.GetMin().y ) };
  const ULONG alAt[3] = { lFrom, lTo, lPole };
  for (ULONG ii=0; ii<3; ii++)
    {
      SmPoint3d sAt;
      SER( pSurface->EvaluatePoint( asUV[ii], sAt ));
      if ( sAt.DistanceBetween( m_vVertices[alAt[ii]]->GetPoint() ) > dTol )
        { SER( SM_ERR ); }
    }

  SER( pSurface->CreateNaturalUVTrimCurves( m_crContext, sNurbDomain, SM_SP_V, TRUE, m_vBall3d, m_vBallUV, m_vBallOrients ));
  if ( m_vBall3d.GetSize() != 3 || m_vBallUV.GetSize() != 3 || m_vBallOrients.GetSize() != 3 )
    { SER( SM_ERR ); } // the three sections; the side at the pole a point
  return SM_SUCCESS;

} // end SmBlendCorner::BuildBallSurface

void SmBlendCorner::DeleteBallSurface()
{
  for (ULONG ii=0; ii<m_vBall3d.GetSize(); ii++) { delete m_vBall3d[ii]; }
  for (ULONG ii=0; ii<m_vBallUV.GetSize(); ii++) { delete m_vBallUV[ii]; }
  m_vBall3d.ReSet();
  m_vBallUV.ReSet();
  m_vBallOrients.ReSet();
  if ( m_pBallSurface != NULL ) { delete m_pBallSurface; m_pBallSurface = NULL; }
}

/*******************************************************************//**
PURPOSE: The Ball's face: the sphere bounded by its natural loop. Its edges
    and the rounds' meet in the stitch and the merge.
***********************************************************************/
SmStatus SmBlendCorner::CreateBallPatch()
{
  NER( m_pBallSurface );
  SmBrep * pFilletBrep = m_pExecutive->GetFilletBrep();
  NER( pFilletBrep );
  SmTArray<ULONG> sCurveLoops;
  sCurveLoops.Add( m_vBall3d.GetSize() );
  SmTArray<SmPoint3d> sLoopPoints;
  SmRegion  * pNewRegion = NULL;
  SmShell   * pNewShell  = NULL;
  SmFace    * pNewFace   = NULL;
  SmSurface * pSurface   = m_pBallSurface;
  // the face takes the surface and the curves
  SmTArray<SmCurve*>        s3d( m_vBall3d );
  SmTArray<SmBSplineCurve*> sUV( m_vBallUV );
  m_pBallSurface = NULL;
  m_vBall3d.ReSet();
  m_vBallUV.ReSet();
  SER( pFilletBrep->MakeFaceWithCurves( NULL, sCurveLoops, &s3d, &sUV, m_vBallOrients, sLoopPoints, pSurface,
                                        pSurface->GetNaturalUVDomain(), SM_OT_SAME, pNewRegion, pNewShell, pNewFace ));
  NER( pNewFace );
  m_vSurfaces.Add( pSurface ); // the corner's patch face (SmFilletCorner::GetPatchFaces, the stitch)
  return SM_SUCCESS;

} // end SmBlendCorner::CreateBallPatch
