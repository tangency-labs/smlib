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
* Its face is the sphere's natural patch: the pole at the contact P_k of the
* plane square to the other two (n_k . n_i = n_k . n_j = 0), so that the
* sections through P_k are meridians and the section P_i P_j is the equator:
*   S(u, v) = C + r (cos v (cos u X + sin u Y) + sin v Z),
* Z = (P_k - C) / r, X = (P_i - C) / r, Y = Z x X, u in [0, theta] (theta the
* angle from X to P_j - C), v in [0, 90] degrees. Its loop is that patch's
* natural boundary, the pole left out: two meridians and the equator arc.
* Where no plane is square to the other two, no side of the triangle is a
* meridian and the equator at once: a great circle off them has no rational
* curve in the sphere's NURB parameters. That corner is not a Ball yet.
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

} // namespace

/*******************************************************************//**
PURPOSE: Whether a filleted vertex is a Ball: every edge there blended (three
    of them), each a constant circular round of the same radius, all convex
    or all concave, the three faces planes, one of them square to the other
    two, and the caller not setting the vertex back (a setback corner).
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

  // three planes, one square to the other two: off square by an angle whose sine is s, a section
  // through its pole leaves the meridian's plane by up to s r, held to the fillet's tolerance
  SmFace * apFace[3];
  sm_BallRingFaces( crAllEUs, apFace );
  SmVector3d asNormal[3];
  for (ULONG ii=0; ii<3; ii++)
    {
      if ( apFace[ii] == NULL || apFace[ii]->GetSurface() == NULL || !apFace[ii]->GetSurface()->IsPlanar( dTol, NULL, &asNormal[ii] ))
        { return SM_SUCCESS; }
    }
  if ( apFace[0] == apFace[1] || apFace[1] == apFace[2] || apFace[2] == apFace[0] || sm_BallPole( asNormal, dTol / dRadius ) == 3 )
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
    { SER( SM_ERR ); } // ClassifyBall found one

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
