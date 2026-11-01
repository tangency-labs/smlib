// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

/**********************************************************************//**
* FILE NAME --- SmBlendCornerBall.cpp
* PURPOSE: SolidCraft change. The blend corner's Ball event
*   (SmBlendCorner.h): three constant rounds of one radius r and one
*   convexity at a vertex of three planes.
*
* The ball of radius r touching the three planes F_i (unit normals n_i) is
* centred at C, the point at distance r from each on the side the rounds'
* balls roll (the triple offset point): n_i . C = n_i . p_i + s_i r. It
* touches F_i at P_i = C - s_i r n_i. The round between F_i and F_j ends on
* its section through C, the circle of its ball there: centred at C, radius
* r, in the plane of n_i and n_j through C, from P_i to P_j. So all three end
* sections are great circles of the one sphere about C, and each round is
* tangent to that sphere along its end section (both are tangent to the
* ball there). The corner's face is the sphere's triangle P_0 P_1 P_2
* between the three sections: exact, G1 to every round.
*
* Where one plane is square to the other two (n_k . n_i = n_k . n_j = 0) the
* face is the SmSphere's natural patch: the pole at that plane's contact
* P_k, so that the sections through P_k are meridians and the section
* P_i P_j is the equator:
*   S(u, v) = C + r (cos v (cos u X + sin u Y) + sin v Z),
* Z = (P_k - C) / r, X = (P_i - C) / r, Y = Z x X, u in [0, theta] (theta the
* angle from X to P_j - C), v in [0, 90] degrees. Elsewhere no side of the
* triangle is the equator of a pole at another vertex, and a great circle
* off the meridians and the equator has no rational curve in the SmSphere's
* NURB parameters: the face is one rational patch of degree (4, 2) on the
* sphere whose natural boundary is the triangle (sm_BallTrianglePatch). Its
* loop, either way, is that patch's natural boundary, the pole left out: the
* three sections.
**********************************************************************/

#include "StdAfx.h"

#ifndef __SMBLENDCORNER_H__
#include <SmBlendCorner.h>
#endif

#ifndef __SMFILLETSOLVER_H__
#include <SmFilletSolver.h>
#endif

#include <SmAxis2Placement.h>
#include <SmOffsetSurface.h>
#include <SmSphere.h>

namespace
{

// The faces of a corner ring of three edgeuses: face i between edgeuses i - 1 and i (as the
// corner's vertex i). NULL entries where two edgeuses share no face.
void sm_BallRingFaces(const SmTArray<SmEdgeuse*> & crAllEUs, SmFace * apFace[3])
{
  for (ULONG ii=0; ii<3; ii++)
    { apFace[ii] = sm_GetFaceByTwoEdgeuses( crAllEUs[(ii+2)%3], crAllEUs[ii], NULL ); }
}

// The index k of the plane square to the other two (|n_k . n_i| and |n_k . n_j| at most dSinTol),
// or 3 when none is.
ULONG sm_BallPole(const SmVector3d asNormal[3], double dSinTol)
{
  for (ULONG kk=0; kk<3; kk++)
    {
      if (   smos_Fabs( asNormal[kk].Dot( asNormal[(kk+1)%3] )) <= dSinTol
          && smos_Fabs( asNormal[kk].Dot( asNormal[(kk+2)%3] )) <= dSinTol )
        { return kk; }
    }
  return 3;
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

  // three planes
  SmFace * apFace[3];
  sm_BallRingFaces( crAllEUs, apFace );
  SmVector3d asNormal[3];
  for (ULONG ii=0; ii<3; ii++)
    {
      if ( apFace[ii] == NULL || apFace[ii]->GetSurface() == NULL || !apFace[ii]->GetSurface()->IsPlanar( dTol, NULL, &asNormal[ii] ))
        { return SM_SUCCESS; }
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
PURPOSE: The Ball's sphere and its loop (see the file's head), from the
    three planes and the processed contacts. The centre is the planes'
    triple offset point, solved from them; the rounds' balls at the
    contacts are checked to be it.
***********************************************************************/
SmStatus SmBlendCorner::BuildBallSurface()
{
  DeleteBallSurface();
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
  const double dRadius = smos_Fabs( pFG0->GetOffsetSurface( 0 )->GetOffsetDistance() );

  // the planes, and the centre r from each on the ball's side: n_i . C = n_i . p_i + s_i r
  SmFace * apFace[3];
  sm_BallRingFaces( *m_pAllEUs, apFace );
  SmVector3d asNormal[3];
  double     adRhs[3], adSide[3];
  for (ULONG ii=0; ii<3; ii++)
    {
      NER( apFace[ii] );
      NER( apFace[ii]->GetSurface() );
      SmPoint3d sOnPlane;
      if ( !apFace[ii]->GetSurface()->IsPlanar( dTol, &sOnPlane, &asNormal[ii] ))
        { SER( SM_ERR ); } // ClassifyBall took the faces to be planes
      SER( asNormal[ii].Unitize() );
      adSide[ii] = ( asNormal[ii].Dot( sBall - sOnPlane ) >= 0.0 ) ? 1.0 : -1.0;
      adRhs[ii]  = asNormal[ii].x * sOnPlane.x + asNormal[ii].y * sOnPlane.y + asNormal[ii].z * sOnPlane.z + adSide[ii] * dRadius;
    }
  // Cramer's rule: C = (rhs_0 n_1 x n_2 + rhs_1 n_2 x n_0 + rhs_2 n_0 x n_1) / (n_0 . n_1 x n_2)
  const double dDet = asNormal[0].Dot( asNormal[1] * asNormal[2] );
  if ( smos_Fabs( dDet ) <= SM_EFF_ZERO )
    { SER( SM_ERR ); }
  const SmVector3d sC = ( asNormal[1] * asNormal[2] ) * ( adRhs[0] / dDet ) + ( asNormal[2] * asNormal[0] ) * ( adRhs[1] / dDet )
                       + ( asNormal[0] * asNormal[1] ) * ( adRhs[2] / dDet );
  const SmPoint3d sCentre( sC.x, sC.y, sC.z );
  if ( sCentre.DistanceBetween( sBall ) > dTol )
    { SER( SM_ERR ); } // the rounds' ball is the planes' triple offset point
  SmPoint3d asContact[3];
  for (ULONG ii=0; ii<3; ii++)
    {
      asContact[ii] = sCentre - asNormal[ii] * ( adSide[ii] * dRadius );
      if ( asContact[ii].DistanceBetween( m_vVertices[ii]->GetPoint() ) > dTol )
        { SER( SM_ERR ); } // vertex i is the ball's contact with face i
      // the contact in closed form, not as the rails' crossing was solved: the rounds' end
      // sections are made through it (SmFilletEdge::CalcCrossSection)
      m_vVertices[ii]->SetPoint( asContact[ii] );
    }

  // The face's normal is its surface's (SmBrep::OrientTrimmedSurfaces turns a face the other way
  // by swapping its surface's parameters): out of the material, towards the centre where the
  // rounds add material (concave), away from it where they remove it.
  SmBoolean bAddsMaterial = FALSE;
  SER( SmFilletSolver::IsConcaveEdge( *(*m_pAllEUs)[0]->GetEdge(), bAddsMaterial ));

  // the corners of the face's natural patch: the contacts at the pole and at the ends of the
  // side opposite it
  ULONG lPole = sm_BallPole( asNormal, dTol / dRadius ), lFrom = 0, lTo = 0;
  SmSurface * pSurface = NULL;
  if ( lPole != 3 )
    {
      // the sphere, the pole's contact on the plane square to the other two, the equator through
      // the others
      lFrom = (lPole+1)%3;
      lTo   = (lPole+2)%3;
      SmVector3d sZ = asContact[lPole] - sCentre;
      SER( sZ.Unitize() );
      SmVector3d sX = asContact[lFrom] - sCentre;
      sX = sX - sZ * sX.Dot( sZ );
      SER( sX.Unitize() );
      SmVector3d sY = sZ * sX;
      SmVector3d sTo = asContact[lTo] - sCentre;
      if ( sTo.Dot( sY ) < 0.0 )
        {
          // the equator runs from the other contact
          ULONG lSwap = lFrom; lFrom = lTo; lTo = lSwap;
          sX = asContact[lFrom] - sCentre;
          sX = sX - sZ * sX.Dot( sZ );
          SER( sX.Unitize() );
          sY = sZ * sX;
          sTo = asContact[lTo] - sCentre;
        }
      const double dTheta = SM_RAD2DEG( atan2( sTo.Dot( sY ), sTo.Dot( sX )));
      if ( dTheta <= 0.0 || dTheta >= 180.0 )
        { SER( SM_ERR ); }
      // The sphere's normal points away from its centre unless its NURB parameters are swapped.
      SmExtent2d sDomain;
      SER( sDomain.SetMinMax( 0.0, 0.0, dTheta, 90.0 ));
      SmSphere * pSphere = new (m_crContext) SmSphere( sCentre, sX, sY, sDomain, dRadius, bAddsMaterial, FALSE, NULL, &m_crContext );
      NER( pSphere );
      m_pBallSurface = pSurface = pSphere;
      SER( pSphere->MakeNurb() );
    }
  else
    {
      // one rational patch on the sphere, the pole at contact 0, its normal turned as the face's
      lPole = 0;
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
              sClean.Clear();
              m_pBallSurface = pSurface = pPatch;
            }
        }
      NER( pSurface );
    }

  // its corners are the contacts: the pole at v = vmax (the sphere) or v = vmin (the patch)
  const SmExtent2d sNurbDomain = pSurface->GetNaturalUVDomain();
  SmPoint3d sAt;
  if ( lPole != 3 && pSurface->GetType() == SmSphere_TYPE )
    {
      SmSphere * pSphere = SM_CAST_PTR( SmSphere, pSurface );
      NER( pSphere );
      const SmExtent2d sStep = pSphere->GetSTEPUVDomain();
      const SmPoint2d asStep[3] = { SmPoint2d( sStep.GetMin().x, sStep.GetMin().y ), SmPoint2d( sStep.GetMax().x, sStep.GetMin().y ),
                                    SmPoint2d( sStep.GetMin().x, sStep.GetMax().y ) };
      const ULONG alAt[3] = { lFrom, lTo, lPole };
      for (ULONG ii=0; ii<3; ii++)
        {
          SER( pSphere->EvaluateSTEPPoint( asStep[ii], sAt ));
          if ( sAt.DistanceBetween( m_vVertices[alAt[ii]]->GetPoint() ) > dTol )
            { SER( SM_ERR ); }
        }
    }
  else
    {
      const SmPoint2d asUV[3] = { SmPoint2d( sNurbDomain.GetMin().x, sNurbDomain.GetMax().y ), SmPoint2d( sNurbDomain.GetMax().x, sNurbDomain.GetMax().y ),
                                  SmPoint2d( sNurbDomain.GetMin().x, sNurbDomain.GetMin().y ) };
      const ULONG alAt[3] = { lFrom, lTo, lPole };
      for (ULONG ii=0; ii<3; ii++)
        {
          SER( pSurface->EvaluatePoint( asUV[ii], sAt ));
          if ( sAt.DistanceBetween( m_vVertices[alAt[ii]]->GetPoint() ) > dTol )
            { SER( SM_ERR ); }
        }
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
