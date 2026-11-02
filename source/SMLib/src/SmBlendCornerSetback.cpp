// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

/**********************************************************************//**
* FILE NAME --- SmBlendCornerSetback.cpp
* PURPOSE: SolidCraft change. The blend corner's Setback (SmBlendCorner.h):
*   where every edge at a vertex is blended, three or more, each blend
*   stops at its setback from the vertex and the hole they leave is filled
*   by SmBlendSetback's setback split, to G1, or to G2 for curvature-
*   continuous sections (or rounds eased into them).
*
* The corner's topology in the executive's pseudo Brep, for faces F_i at
* the vertex (F_i between blended edges i - 1 and i, m_pAllEUs's order):
* on each face the rail end of the blend before it (M_i, SM_FV_MATE of
* S_(i-1)), the fill's base point there (B_i, SM_FV_PRECOMPUTED) and the
* rail end of the blend after it (S_i, SM_FV_SETBACK: its section square to
* its edge at the setback); edges M_i-B_i and B_i-S_i, the two halves of the
* face's spring curve (SM_FE_PRECOMPUTED, on F_i), and S_i-M_(i+1), blend
* i's end section (SM_FE_CROSS_SECTION). The rails and sections are the
* executive's; the fill, its base points and its spring curves are built
* once the blends' surfaces exist (and are eased, for an eased setback):
* the spring curves are the setback patches' sides on the faces, the base
* points their inner corners. The fill's 2n patches become faces of their
* own, bounded by their iso-curves, and are stitched to the blends' faces
* and to each other with the rest.
**********************************************************************/

#include "StdAfx.h"

#ifndef __SMBLENDCORNER_H__
#include <SmBlendCorner.h>
#endif

#ifndef __SMFILLETSOLVER_H__
#include <SmFilletSolver.h>
#endif

#include <SmBSplineSurface.h>
#include <SmLine.h>
#include <SmOffsetSurface.h>
#include <SmFilletExact.h>

namespace
{

// The unit direction of pEdge away from cpVertex, from its curve's end tangent there.
SmStatus sm_SbAwayFrom(const SmEdge * pEdge, const SmVertex * cpVertex, SmVector3d & rAway)
{
  NER( pEdge );
  NER( pEdge->GetCurve() );
  const SmExtent1d sIvl = pEdge->GetInterval();
  const SmBoolean bAtStart = ( pEdge->GetStartVertex() == cpVertex ) ? TRUE : FALSE;
  SmVector3d sPV[2];
  SER( pEdge->GetCurve()->Evaluate( bAtStart ? sIvl.GetMin() : sIvl.GetMax(), 1, TRUE, sPV ));
  rAway = bAtStart ? sPV[1] : -sPV[1];
  SER( rAway.Unitize() );
  return SM_SUCCESS;
}

// The solver of a blended edgeuse at the corner (or of its mate's side).
SmFilletSolver * sm_SbSolverOf(SmFilletExecutive * pExec, SmEdgeuse * pEU)
{
  SmFilletSolver * pFS = pExec->GetFilletSolverOfEdgeuse( pEU );
  if ( pFS == NULL && pEU != NULL )
    { pFS = pExec->GetFilletSolverOfEdgeuse( pEU->GetMate() ); }
  return pFS;
}

// Of a corner whose edges all have one convexity, the planes of its faces with their normals out
// of the material, and its edges for SmBlendSetback::Distances: face i between edges i - 1 and i
// (crAllEUs's order round the vertex), edge i between faces i and i + 1.
SmStatus sm_SbQuadricFace(const SmFace * pFace, const SmVertex * cpVertex, const SmVector3d & crOut, double dTol,
                          SmBlendSetback::Face & rFace, SmBoolean & rbOK);

// (With pbCurved, a face may also be a natural quadric (sm_SbQuadricFace), *pbCurved TRUE where one is.)
SmStatus sm_SbCornerPlanes(const SmVertex * cpVertex, const SmTArray<SmEdgeuse*> & crAllEUs, SmFilletExecutive * pExec, double dTol,
                           SmBoolean bConvex, std::vector<SmFace*> & rvFaces, std::vector<SmBlendSetback::Face> & rvPlanes,
                           std::vector<SmBlendSetback::Edge> & rvEdges, SmBoolean & rbPlanes, SmBoolean * pbCurved = NULL)
{
  rbPlanes = FALSE;
  if ( pbCurved != NULL ) { *pbCurved = FALSE; }
  const ULONG n = crAllEUs.GetSize();
  rvFaces.assign( n, NULL );
  rvPlanes.assign( n, SmBlendSetback::Face() );
  rvEdges.assign( n, SmBlendSetback::Edge() );
  std::vector<SmVector3d> sAway( n );
  for (ULONG i=0; i<n; i++)
    { SER( sm_SbAwayFrom( crAllEUs[i]->GetEdge(), cpVertex, sAway[i] )); }
  for (ULONG i=0; i<n; i++)
    {
      SmFace * pFace = sm_GetFaceByTwoEdgeuses( crAllEUs[(i + n - 1) % n], crAllEUs[i], NULL );
      if ( pFace == NULL || pFace->GetSurface() == NULL )
        { return SM_SUCCESS; }
      SmPoint3d  sOn;
      SmVector3d sNormal;
      SmBlendSetback::Face sQuadric;
      SmBoolean bQuadric = FALSE;
      if ( !pFace->GetSurface()->IsPlanar( dTol, &sOn, &sNormal ) || sNormal.Unitize( FALSE ) != SM_SUCCESS )
        {
          if ( pbCurved == NULL )
            { return SM_SUCCESS; }
          // a natural quadric: its normal at the vertex (its side settled below)
          SER( sm_SbQuadricFace( pFace, cpVertex, SmVector3d( 0.0, 0.0, 0.0 ), dTol, sQuadric, bQuadric ));
          if ( !bQuadric )
            { return SM_SUCCESS; }
          SmBlendSetback::Face sProbe( sQuadric );
          sProbe.m_dSide = 1.0;
          (void)SmBlendSetback::Distance( sProbe, cpVertex->GetPoint(), &sNormal );
          if ( sNormal.Unitize( FALSE ) != SM_SUCCESS )
            { return SM_SUCCESS; }
        }
      // the other edges lie on the material's side of a convex corner's face, off it at a
      // concave one's
      SmVector3d sOthers( 0.0, 0.0, 0.0 );
      for (ULONG j=0; j<n; j++)
        {
          if ( j != i && j != (i + n - 1) % n )
            { sOthers = sOthers + sAway[j]; }
        }
      const double dSide = sNormal.Dot( sOthers );
      if ( smos_Fabs( dSide ) <= SM_EFF_ZERO )
        { return SM_SUCCESS; }
      const SmBoolean bFlip = ( ( dSide > 0.0 ) == ( bConvex != FALSE )) ? TRUE : FALSE;
      if ( bFlip )
        { sNormal = -sNormal; }
      rvFaces[i] = pFace;
      if ( bQuadric )
        {
          // the material inside the quadric (+1) where the normal out of it points away from its axis
          rvPlanes[i] = sQuadric;
          rvPlanes[i].m_sNormal = sNormal;
          rvPlanes[i].m_dSide   = bFlip ? -1.0 : 1.0;
          *pbCurved = TRUE;
        }
      else
        {
          rvPlanes[i].m_sPoint  = cpVertex->GetPoint();
          rvPlanes[i].m_sNormal = sNormal;
        }
    }
  for (ULONG i=0; i<n; i++)
    {
      SmFilletSolver * pFS = sm_SbSolverOf( pExec, crAllEUs[i] );
      NER( pFS );
      NER( pFS->GetSurface( 0 ));
      rvEdges[i].m_sAway       = sAway[i];
      rvEdges[i].m_dRadius     = smos_Fabs( pFS->GetSurface( 0 )->GetOffsetDistance() );
      rvEdges[i].m_alFaces[0]  = i;
      rvEdges[i].m_alFaces[1]  = (i + 1) % n;
    }
  rbPlanes = TRUE;
  return SM_SUCCESS;
}

// A curved corner face as the fill sees it (SmBlendSetback::Face): its natural quadric, its side
// the material's, from the face's normal out of the material at the vertex (crOut). FALSE where the
// face is no plane or natural quadric.
SmStatus sm_SbQuadricFace(const SmFace * pFace, const SmVertex * cpVertex, const SmVector3d & crOut, double dTol,
                          SmBlendSetback::Face & rFace, SmBoolean & rbOK)
{
  rbOK = FALSE;
  NER( pFace );
  SmBlendQuadric sQ;
  SER( SmBlendQuadric::FromSurface( pFace->GetSurface(), pFace->GetUVDomain().GetMid(), dTol, sQ ));
  rFace = SmBlendSetback::Face();
  rFace.m_sNormal = crOut;
  const SmPoint3d sV = cpVertex->GetPoint();
  switch ( sQ.m_eKind )
    {
    case SmBlendQuadric::SM_BQ_PLANE:
      rFace.m_eKind  = SmBlendSetback::SM_SB_PLANE;
      rFace.m_sPoint = sV;
      break;
    case SmBlendQuadric::SM_BQ_CYLINDER:
      rFace.m_eKind   = SmBlendSetback::SM_SB_CYLINDER;
      rFace.m_sPoint  = sQ.m_sPoint;
      rFace.m_sAxis   = sQ.m_sDir;
      rFace.m_dRadius = sQ.m_dRadius;
      break;
    case SmBlendQuadric::SM_BQ_CONE:
      rFace.m_eKind   = SmBlendSetback::SM_SB_CONE;
      rFace.m_sPoint  = sQ.m_sPoint; // the apex
      rFace.m_sAxis   = sQ.m_sDir;   // into the nappe
      rFace.m_dRadius = 0.0;
      rFace.m_dSlope  = smos_Tangent( sQ.m_dAngle );
      break;
    case SmBlendQuadric::SM_BQ_SPHERE:
      rFace.m_eKind   = SmBlendSetback::SM_SB_SPHERE;
      rFace.m_sPoint  = sQ.m_sPoint;
      rFace.m_dRadius = sQ.m_dRadius;
      break;
    default:
      return SM_SUCCESS;
    }
  // the material inside the quadric (+1) where the normal out of it at the vertex points away from
  // the axis or centre
  if ( rFace.m_eKind != SmBlendSetback::SM_SB_PLANE )
    {
      SmVector3d sAway = sQ.Normal( sV ) * sQ.m_dSide; // the quadric's own outward (away from its axis)
      rFace.m_dSide = ( sAway.Dot( crOut ) >= 0.0 ) ? 1.0 : -1.0;
    }
  rbOK = TRUE;
  return SM_SUCCESS;
}

// Where, along pEdge from cpVertex, the plane square to the edge holds crC (Newton on the edge's
// parameter from the vertex), as the edge's arc length from the vertex.
SmStatus sm_SbAlongEdge(const SmEdge * pEdge, const SmVertex * cpVertex, const SmPoint3d & crC, double dTol, double & rdAlong)
{
  NER( pEdge );
  const SmCurve * pCurve = pEdge->GetCurve();
  NER( pCurve );
  const SmExtent1d sIvl = pEdge->GetInterval();
  const SmBoolean bFromStart = ( pEdge->GetStartVertex() == cpVertex ) ? TRUE : FALSE;
  const double t0 = bFromStart ? sIvl.GetMin() : sIvl.GetMax();
  double t = t0;
  for (int it=0; it<60; it++)
    {
      SmVector3d sD[3];
      SER( pCurve->Evaluate( t, 2, TRUE, sD ));
      const SmVector3d sW = crC - SmPoint3d( sD[0].x, sD[0].y, sD[0].z );
      const double f = sW.Dot( sD[1] ), df = -sD[1].Dot( sD[1] ) + sW.Dot( sD[2] );
      if ( !( smos_Fabs( df ) > 0.0 ))
        { SER( SM_ERR ); }
      const double dStep = f / df;
      t -= dStep;
      t = smos_Max( sIvl.GetMin(), smos_Min( sIvl.GetMax(), t ));
      if ( smos_Fabs( dStep ) * sD[1].Length() <= 1.0e-3 * dTol )
        { break; }
    }
  if ( smos_Fabs( t - t0 ) == 0.0 )
    { rdAlong = 0.0; return SM_SUCCESS; }
  SER( pCurve->Length( SmExtent1d( smos_Min( t, t0 ), smos_Max( t, t0 )), 1.0e-3 * dTol, rdAlong ));
  return SM_SUCCESS;
}

// The least setbacks at a corner beside curved faces: on each face, where the rails of the two
// rounds beside it cross (the point P of the face whose ball centres P - s r n(P), s +1 at a convex
// corner, -1 at a concave one, lie at r from each round's other face), as far along each round's
// edge as the plane square to it through that ball's centre; each edge the farther of its two faces'.
SmStatus sm_SbCurvedLeast(const SmVertex * cpVertex, const SmTArray<SmEdgeuse*> & crAllEUs, const std::vector<SmBlendSetback::Face> & crFaces,
                          const std::vector<double> & crRadii, SmBoolean bConvex, double dTol, std::vector<double> & rvLeast)
{
  const ULONG n = crAllEUs.GetSize();
  rvLeast.assign( n, 0.0 );
  const double s = bConvex ? 1.0 : -1.0;
  for (ULONG i=0; i<n; i++)
    {
      const ULONG a = (i + n - 1) % n, b = i; // the rounds beside face i: edge a (faces a, i), edge b (faces i, i + 1)
      const SmBlendSetback::Face & crFi = crFaces[i], & crFa = crFaces[a], & crFb = crFaces[(i + 1) % n];
      const double ra = crRadii[a], rb = crRadii[b];
      auto residual = [&](const SmVector3d & crP, double (&rF)[3], SmPoint3d & rCa, SmPoint3d & rCb)
        {
          SmVector3d sN;
          rF[0] = SmBlendSetback::Distance( crFi, SmPoint3d( crP.x, crP.y, crP.z ), &sN );
          rCa = SmPoint3d( crP.x, crP.y, crP.z ) - sN * ( s * ra );
          rCb = SmPoint3d( crP.x, crP.y, crP.z ) - sN * ( s * rb );
          rF[1] = SmBlendSetback::Distance( crFa, rCa ) + s * ra;
          rF[2] = SmBlendSetback::Distance( crFb, rCb ) + s * rb;
        };
      SmVector3d sP( cpVertex->GetPoint() );
      SmPoint3d sCa, sCb;
      const double dScale = smos_Max( 1.0, smos_Max( ra, rb ));
      SmBoolean bDone = FALSE;
      for (int it=0; it<60 && !bDone; it++)
        {
          double F[3];
          residual( sP, F, sCa, sCb );
          if ( smos_Fabs( F[0] ) + smos_Fabs( F[1] ) + smos_Fabs( F[2] ) <= 1.0e-13 * dScale )
            { bDone = TRUE; break; }
          // the Jacobian by central differences of the closed-form distances
          double J[3][3];
          const double h = 1.0e-6 * dScale;
          for (int k=0; k<3; k++)
            {
              SmVector3d sD( 0.0, 0.0, 0.0 );
              if ( k == 0 ) sD.x = h; else if ( k == 1 ) sD.y = h; else sD.z = h;
              double Fp[3], Fm[3];
              SmPoint3d sX1, sX2;
              residual( sP + sD, Fp, sX1, sX2 );
              residual( sP - sD, Fm, sX1, sX2 );
              for (int r=0; r<3; r++) { J[r][k] = ( Fp[r] - Fm[r] ) / ( 2.0 * h ); }
            }
          const double det = J[0][0] * ( J[1][1] * J[2][2] - J[1][2] * J[2][1] ) - J[0][1] * ( J[1][0] * J[2][2] - J[1][2] * J[2][0] )
                           + J[0][2] * ( J[1][0] * J[2][1] - J[1][1] * J[2][0] );
          if ( !( smos_Fabs( det ) > SM_EFF_ZERO ))
            { SER( SM_ERR ); }
          double x[3];
          for (int k=0; k<3; k++)
            {
              double M[3][3];
              for (int r=0; r<3; r++) for (int c=0; c<3; c++) { M[r][c] = ( c == k ) ? F[r] : J[r][c]; }
              x[k] = ( M[0][0] * ( M[1][1] * M[2][2] - M[1][2] * M[2][1] ) - M[0][1] * ( M[1][0] * M[2][2] - M[1][2] * M[2][0] )
                     + M[0][2] * ( M[1][0] * M[2][1] - M[1][1] * M[2][0] )) / det;
            }
          sP = sP - SmVector3d( x[0], x[1], x[2] );
        }
      if ( !bDone )
        { SER( SM_ERR ); }
      double dA = 0.0, dB = 0.0;
      SER( sm_SbAlongEdge( crAllEUs[a]->GetEdge(), cpVertex, sCa, dTol, dA ));
      SER( sm_SbAlongEdge( crAllEUs[b]->GetEdge(), cpVertex, sCb, dTol, dB ));
      rvLeast[a] = smos_Max( rvLeast[a], dA );
      rvLeast[b] = smos_Max( rvLeast[b], dB );
    }
  return SM_SUCCESS;
}

// Whether pEdge is blended by the executive (any of its edgeuses has a solver).
SmBoolean sm_SbEdgeBlended(SmFilletExecutive * pExec, const SmEdge * pEdge)
{
  SmEdgeuse * pEU = pEdge->GetPrimaryEdgeuse();
  for (int i=0; i<2 && pEU != NULL; i++, pEU = pEU->GetRadial())
    {
      if ( pExec->GetFilletSolverOfEdgeuse( pEU ) != NULL || pExec->GetFilletSolverOfEdgeuse( pEU->GetMate() ) != NULL )
        { return TRUE; }
    }
  return FALSE;
}

// The corner's setbacks (owner, 2026-11-01). Each edge's room: its length, less the far end's
// claim: an asked-for setback with its ease there, or, where the far end is itself a vertex
// whose three or more edges are all blended (a corner that may set back too), half the edge.
// Asked for (bAsked): dAsked on every edge, refused (rbFits FALSE) unless it and dEase lie below
// every room, rdLargest then the largest that fits (the least room less the ease); and below
// an edge's contact-line minimum (rbShort). Else the rule's (Distances), each shrunk to its
// room less a quarter of its radius where longer, as long as that is at or above its minimum
// (rbFits FALSE where not: the rounds beside an edge cross past its room, an overflow). The
// quarter radius is the rule's own margin (least + 0.25 r): one margin, kept from the
// contact-line minimum and from the far end alike, so the round keeps 0.25 r of its edge.
SmStatus sm_SbSetBacks(const SmVertex * cpVertex, const SmTArray<SmEdgeuse*> & crAllEUs, SmFilletExecutive * pExec,
                       const std::vector<SmBlendSetback::Face> & crPlanes, const std::vector<SmBlendSetback::Edge> & crEdges, double dTol,
                       std::vector<double> & rvSetBacks, std::vector<double> & rvLeast, std::vector<double> & rvRoom,
                       SmBoolean & rbFits, SmBoolean & rbShort, double & rdLargest, const TCHAR *& rpWhy)
{
  rbFits = TRUE;
  rbShort = FALSE;
  rdLargest = SM_BIG_DOUBLE;
  const ULONG n = crAllEUs.GetSize();
  bool bCurved = false;
  for (size_t f=0; f<crPlanes.size(); f++) { bCurved = bCurved || crPlanes[f].m_eKind != SmBlendSetback::SM_SB_PLANE; }
  if ( bCurved )
    {
      // beside a curved face: the least setbacks where the rails cross on the faces, the rule's
      // max(1.5 r, least + 0.25 r) (SmBlendSetback::Distances's, on the faces as they are)
      std::vector<double> sRadii( n );
      for (ULONG i=0; i<n; i++) { sRadii[i] = crEdges[i].m_dRadius; }
      SmBoolean bConcave = FALSE;
      SER( SmFilletSolver::IsConcaveEdge( *crAllEUs[0]->GetEdge(), bConcave ));
      if ( sm_SbCurvedLeast( cpVertex, crAllEUs, crPlanes, sRadii, bConcave ? FALSE : TRUE, dTol, rvLeast ) != SM_SUCCESS )
        {
          rpWhy = _T("Fillet Error: where the rails of the rounds at a setback corner cross on its curved faces cannot be found");
          return SM_ERR;
        }
      rvSetBacks.assign( n, 0.0 );
      for (ULONG i=0; i<n; i++) { rvSetBacks[i] = smos_Max( 1.5 * sRadii[i], rvLeast[i] + 0.25 * sRadii[i] ) * ( getenv("SBSCALE") ? atof(getenv("SBSCALE")) : 1.0 ); }
    }
  else
    { SER( SmBlendSetback::Distances( crPlanes, crEdges, rvSetBacks, &rvLeast, &rpWhy )); }
  double dAsked = 0.0, dEase = 0.0;
  SmBoolean bOpen = FALSE;
  const SmBoolean bAsked = pExec->GetSetBackFill( cpVertex, dAsked, dEase, bOpen );
  rvRoom.assign( n, 0.0 );
  for (ULONG i=0; i<n; i++)
    {
      const SmEdge * pEdge = crAllEUs[i]->GetEdge();
      NER( pEdge );
      NER( pEdge->GetStartVertex() );
      NER( pEdge->GetEndVertex() );
      const double dLength = pEdge->GetStartVertex()->GetPoint().DistanceBetween( pEdge->GetEndVertex()->GetPoint() );
      const SmVertex * pFar = pEdge->GetOtherVertex( cpVertex );
      double dFarSetBack = 0.0, dFarEase = 0.0;
      SmBoolean bFarOpen = FALSE;
      double dClaim = 0.0;
      if ( pFar != NULL && pFar != cpVertex && pExec->GetSetBackFill( pFar, dFarSetBack, dFarEase, bFarOpen ))
        { dClaim = dFarSetBack + dFarEase; }
      else if ( pFar != NULL && pFar != cpVertex && !bAsked )
        {
          SmTArray<SmEdge*> sFarEdges;
          ((SmVertex*)pFar)->GetEdges( sFarEdges );
          SmBoolean bAll = ( sFarEdges.GetSize() >= 3 ) ? TRUE : FALSE;
          for (ULONG j=0; j<sFarEdges.GetSize() && bAll; j++)
            { bAll = sm_SbEdgeBlended( pExec, sFarEdges[j] ); }
          if ( bAll )
            { dClaim = 0.5 * dLength; }
        }
      rvRoom[i] = dLength - dClaim;
      if ( bAsked )
        {
          rdLargest = smos_Min( rdLargest, rvRoom[i] - dEase );
          if ( dAsked + dEase >= rvRoom[i] - dTol )
            { rbFits = FALSE; }
          if ( dAsked < rvLeast[i] - dTol )
            { rbShort = TRUE; }
          rvSetBacks[i] = dAsked;
        }
      else
        {
          const double dFit = rvRoom[i] - 0.25 * crEdges[i].m_dRadius;
          if ( rvSetBacks[i] > dFit )
            { rvSetBacks[i] = dFit; }
          if ( rvSetBacks[i] < rvLeast[i] )
            { rbFits = FALSE; }
        }
    }
  return SM_SUCCESS;
}

} // namespace

/*******************************************************************//**
PURPOSE: Whether the corner is a Setback (SmBlendCorner.h).
***********************************************************************/
SmStatus SmBlendCorner::ClassifySetBack
  (const SmVertex             * cpVertex,
   const SmTArray<SmEdgeuse*> & crSolverEUs,
   const SmTArray<SmEdgeuse*> & crAllEUs,
   SmFilletExecutive          * pExec,
   double                       dTangencyTolDeg,
   SmBlendCornerEvent         & reEvent,
   ULONG                      & rlOrder)
{
  reEvent = SM_BCE_NONE;
  rlOrder = 1;
  NER( cpVertex );
  NER( pExec );
  const ULONG n = crAllEUs.GetSize();
  if ( n < 3 || crSolverEUs.GetSize() != n )
    { return SM_SUCCESS; }

  // the blends: constant rounds of one radius across, or curvature-continuous sections
  double dTol = SM_BIG_DOUBLE, dMinR = SM_BIG_DOUBLE, dMaxR = 0.0;
  SmFilletSurfaceGeneratorType eKind = SM_FSG_UNKNOWN;
  SmBoolean bOneKind = TRUE;
  for (ULONG i=0; i<n; i++)
    {
      SmFilletSolver * pFS = sm_SbSolverOf( pExec, crAllEUs[i] );
      NER( pFS );
      SmFilletSurfaceGenerator * pFSG = pFS->GetFilletSurfaceGenerator();
      NER( pFSG );
      const SmFilletSurfaceGeneratorType eThis = pFSG->GetFilletSurfaceGeneratorType();
      if (   ( eThis != SM_FSG_CIRCULAR && eThis != SM_FSG_BLEND_CURVE )
          || pFS->GetSolverType() != SM_FS_CONST_RADIUS || pFS->OffsetRadiiCanChange() )
        { return SM_SUCCESS; }
      if ( i > 0 && eThis != eKind )
        { bOneKind = FALSE; }
      eKind = eThis;
      const SmOffsetSurface * pOff0 = pFS->GetSurface( 0 );
      const SmOffsetSurface * pOff1 = pFS->GetSurface( 1 );
      NER( pOff0 );
      NER( pOff1 );
      const double dR0 = smos_Fabs( pOff0->GetOffsetDistance() ), dR1 = smos_Fabs( pOff1->GetOffsetDistance() );
      dTol = smos_Min( dTol, pFS->GetThisApproxTol3d() );
      if ( smos_Fabs( dR0 - dR1 ) > pFS->GetThisApproxTol3d() )
        { return SM_SUCCESS; }
      dMinR = smos_Min( dMinR, dR0 );
      dMaxR = smos_Max( dMaxR, dR0 );
    }

  // the edges: of one convexity, no two tangent at the vertex; straight between planes (beside a
  // curved face, any edge)
  SmBoolean bConcave0 = FALSE, bStraight = TRUE;
  for (ULONG i=0; i<n; i++)
    {
      SmEdge * pEdge = crAllEUs[i]->GetEdge();
      NER( pEdge );
      SmPoint3d  sLinePt;
      SmVector3d sDir;
      if ( pEdge->GetCurve() == NULL )
        { return SM_SUCCESS; }
      if ( !pEdge->GetCurve()->IsLine( 10, dTol, sLinePt, sDir ))
        { bStraight = FALSE; }
      SmBoolean bConcave = FALSE;
      SER( SmFilletSolver::IsConcaveEdge( *pEdge, bConcave ));
      if ( i == 0 )
        { bConcave0 = bConcave; }
      else if ( bConcave != bConcave0 )
        { return SM_SUCCESS; }
      for (ULONG j=i+1; j<n; j++)
        {
          if ( sm_TestG1Continuity( cpVertex, pEdge, crAllEUs[j]->GetEdge(), dTangencyTolDeg, NULL, 0 ))
            { return SM_SUCCESS; }
        }
    }

  // the faces: planes
  std::vector<SmFace*> sFaces;
  std::vector<SmBlendSetback::Face> sPlanes;
  std::vector<SmBlendSetback::Edge> sEdges;
  SmBoolean bPlanes = FALSE, bCurved = FALSE;
  SER( sm_SbCornerPlanes( cpVertex, crAllEUs, pExec, dTol, !bConcave0, sFaces, sPlanes, sEdges, bPlanes, &bCurved ));
  if ( !bPlanes || ( !bCurved && !bStraight ))
    { return SM_SUCCESS; }
  // beside a curved face only the G1 fill of constant circular rounds (#85)
  if ( bCurved && ( eKind != SM_FSG_CIRCULAR || !bOneKind ))
    { return SM_SUCCESS; }

  double dSetBack = 0.0, dEase = 0.0;
  SmBoolean bOpen = FALSE;
  const SmBoolean bAsked = pExec->GetSetBackFill( cpVertex, dSetBack, dEase, bOpen );
  if ( bOpen || ( !bAsked && !bOneKind ))
    { return SM_SUCCESS; } // a corner left open keeps the patch its caller removes
  const SmBoolean bG2 = ( bOneKind && eKind == SM_FSG_BLEND_CURVE ) ? TRUE : FALSE;
  // rounds of one radius meet a common ball where three faces meet: not a setback unless asked
  if ( !bAsked && !bG2 && dMaxR - dMinR <= dTol )
    { return SM_SUCCESS; }
  // (setbacks that do not fit are the Setback's to refuse, or the overflow's: MakeSetBackTopology)
  rlOrder = ( bG2 || ( bAsked && dEase > 0.0 )) ? 2 : 1;
  if ( bCurved && rlOrder != 1 )
    { return SM_SUCCESS; }
  reEvent = SM_BCE_SETBACK;
  return SM_SUCCESS;

} // end SmBlendCorner::ClassifySetBack

/*******************************************************************//**
PURPOSE: A setback the corner does not build: its reason noted at the
    corner's vertex, then SM_ERR.
***********************************************************************/
SmStatus SmBlendCorner::SetBackNotBuilt(const TCHAR * cpWhy)
{
  GetFilletExecutive()->GetFilletErrorInfo()->NoteError( m_cpVertex, SM_FIL_FAILURE, SM_FILERR_CORNER_UNSUPPORTED, cpWhy );
  return SM_ERR;
}

/*******************************************************************//**
PURPOSE: Each blended edge's setback.
***********************************************************************/
double SmBlendCorner::GetSetBackDistAlong(const SmEdge * pEdge)
{
  if ( IsSetBack() && m_pAllEUs != NULL && m_vSetBacks.size() == m_pAllEUs->GetSize() )
    {
      for (ULONG i=0; i<m_pAllEUs->GetSize(); i++)
        {
          if ( (*m_pAllEUs)[i]->GetEdge() == pEdge )
            { return m_vSetBacks[i]; }
        }
    }
  return m_dSetBackDist;
}

/*******************************************************************//**
PURPOSE: The Setback's setbacks and topology (the file's head).

NOTES: The setbacks are those asked for (SmFilletExecutive::SetSetBackFill),
    each at least where the rails beside its edge cross it, or else the
    owner's rule (SmBlendSetback::Distances). Each with its ease must fit on
    its edge (SM_FILERR_NO_ROOM, naming the setback).
***********************************************************************/
SmStatus SmBlendCorner::MakeSetBackTopology()
{
  NER( m_pAllEUs );
  NER( m_pSolverEUs );
  const ULONG n = m_pAllEUs->GetSize();
  if ( n < 3 || m_pSolverEUs->GetSize() != n )
    { SER( SM_ERR ); }
  SmFilletExecutive * pExec = GetFilletExecutive();
  NER( pExec );
  const double dTol = GetThisApproxTol3d();

  SmBoolean bConcave = FALSE;
  SER( SmFilletSolver::IsConcaveEdge( *(*m_pAllEUs)[0]->GetEdge(), bConcave ));
  m_bConvexSetBack = bConcave ? FALSE : TRUE;
  std::vector<SmFace*> sFaces;
  std::vector<SmBlendSetback::Face> sPlanes;
  std::vector<SmBlendSetback::Edge> sEdges;
  SmBoolean bPlanes = FALSE;
  SmBoolean bCurvedFaces = FALSE;
  SER( sm_SbCornerPlanes( m_cpVertex, *m_pAllEUs, pExec, dTol, m_bConvexSetBack, sFaces, sPlanes, sEdges, bPlanes, &bCurvedFaces ));
  if ( !bPlanes )
    { SER( SM_ERR ); } // ClassifySetBack found them

  // the setbacks (sm_SbSetBacks). Where the rule's do not fit (the rounds beside an edge cross
  // past its room) the blends overflow their faces: the executive's overflow comes first, and
  // the corner, made meanwhile with each setback at half its room (only so that its vertices
  // can be placed), refuses in BuildSetBackFill if it is ever built.
  std::vector<double> sLeast, sRoom;
  const TCHAR * pWhy = NULL;
  SmBoolean bFits = TRUE, bShort = FALSE;
  double dLargest = 0.0;
  if ( sm_SbSetBacks( m_cpVertex, *m_pAllEUs, pExec, sPlanes, sEdges, dTol, m_vSetBacks, sLeast, sRoom, bFits, bShort, dLargest, pWhy ) != SM_SUCCESS )
    { return SetBackNotBuilt( pWhy != NULL ? pWhy : _T("Fillet Error: the corner's setbacks cannot be found") ); }
  double dAsked = 0.0, dEase = 0.0;
  SmBoolean bOpen = FALSE;
  const SmBoolean bAsked = pExec->GetSetBackFill( m_cpVertex, dAsked, dEase, bOpen );
  if ( bAsked )
    {
      TCHAR cWhy[SM_TBLOCK_SIZE];
      if ( !bFits )
        {
          smos_sprintf( cWhy, _T("Fillet Error: a setback of %g (with its ease, %g) does not fit on the corner's edges: a setback below %g fits"),
                        dAsked, dAsked + dEase, dLargest );
          pExec->GetFilletErrorInfo()->NoteError( m_cpVertex, SM_FIL_FAILURE, SM_FILERR_NO_ROOM, cWhy );
          return SM_ERR;
        }
      if ( bShort )
        {
          double dMost = 0.0;
          for (ULONG i=0; i<n; i++) { dMost = smos_Max( dMost, sLeast[i] ); }
          smos_sprintf( cWhy, _T("Fillet Error: a setback of %g stops short of where the rails of the rounds beside an edge cross it, %g from the vertex"),
                        dAsked, dMost );
          return SetBackNotBuilt( cWhy );
        }
      m_dEaseDist = dEase;
    }
  m_bSetBackOverflows = bFits ? FALSE : TRUE;
  m_vSetBackRoom = sRoom;
  // a G1 fill is taken only certified convex (SmBlendSetback::FillHole): the least scale of the
  // setbacks at which it is (CertifiedScale, from the corner's planes, edges and radii; owner,
  // 2026-11-01). The rule's grow to it, each held at its room less the rule's quarter radius (the
  // others growing on); refused, naming them, where none within the rooms certifies. One asked
  // for that is too short is refused, naming the setback that builds (or that none does, where
  // that one does not fit with its ease).
  if ( m_lFillOrder == 1 && !m_bSetBackOverflows && !bCurvedFaces )
    {
      std::vector<double> sMost( n );
      for (ULONG i=0; i<n; i++)
        { sMost[i] = bAsked ? SM_BIG_DOUBLE : smos_Max( m_vSetBacks[i], sRoom[i] - 0.25 * sEdges[i].m_dRadius ); }
      double dScale = 0.0;
      const TCHAR * pScaleWhy = NULL;
      const SmStatus eScale = SmBlendSetback::CertifiedScale( m_crContext, sPlanes, sEdges, m_vSetBacks, dScale, NULL, &pScaleWhy, &sMost );
      if ( eScale != SM_SUCCESS )
        { return SetBackNotBuilt( pScaleWhy != NULL ? pScaleWhy : _T("Fillet Error: the setback corner's convex fill cannot be sized") ); }
      TCHAR cWhy[SM_TBLOCK_SIZE];
      if ( dScale == 0.0 )
        {
          if ( bAsked )
            { smos_sprintf( cWhy, _T("Fillet Error: no setback up to twice %g keeps the setback corner's fill convex"), dAsked ); }
          else
            {
              smos_sprintf( cWhy, _T("Fillet Error: no setbacks within the corner's edges keep its fill convex (the most that fit: %g"), sMost[0] );
              for (ULONG i=1; i<n; i++) { TCHAR cMore[64]; smos_sprintf( cMore, _T(", %g"), sMost[i] ); smos_WStrCat( cWhy, cMore ); }
              smos_WStrCat( cWhy, _T(")") );
            }
          pExec->NoteFitSetBack( 0.0 );
          pExec->GetFilletErrorInfo()->NoteError( m_cpVertex, SM_FIL_FAILURE, SM_FILERR_NO_ROOM, cWhy );
          return SM_ERR;
        }
      if ( dScale > 1.0 && bAsked )
        {
          SmBoolean bRoom = TRUE;
          for (ULONG i=0; i<n; i++)
            if ( dScale * dAsked + m_dEaseDist >= sRoom[i] )
              { bRoom = FALSE; }
          const double dBuilds = bRoom ? dScale * dAsked : 0.0;
          if ( bRoom )
            { smos_sprintf( cWhy, _T("Fillet Error: a setback of %g is too short for the setback corner's fill to stay convex: a setback of %g builds it"), dAsked, dBuilds ); }
          else
            { smos_sprintf( cWhy, _T("Fillet Error: a setback of %g is too short for the setback corner's fill to stay convex, and the %g that would be does not fit on its edges"),
                            dAsked, dScale * dAsked ); }
          pExec->NoteFitSetBack( dBuilds );
          pExec->GetFilletErrorInfo()->NoteError( m_cpVertex, SM_FIL_FAILURE, SM_FILERR_NO_ROOM, cWhy );
          return SM_ERR;
        }
      for (ULONG i=0; i<n; i++) { m_vSetBacks[i] = smos_Min( dScale * m_vSetBacks[i], sMost[i] ); }
    }
  if ( m_bSetBackOverflows )
    {
      for (ULONG i=0; i<n; i++) { m_vSetBacks[i] = 0.5 * sRoom[i]; }
    }
  m_dSetBackDist = 0.0;
  for (ULONG i=0; i<n; i++)
    { m_dSetBackDist = smos_Max( m_dSetBackDist, m_vSetBacks[i] ); }

  // the vertices: per face, M (the blend before's rail end), B (the base point), S (the blend after's)
  std::vector<SmFilletGeom*> sGeoms( n );
  for (ULONG i=0; i<n; i++)
    {
      SmFilletSolver * pFS = sm_SbSolverOf( pExec, (*m_pAllEUs)[i] );
      NER( pFS );
      sGeoms[i] = pFS->GetFirstFilletGeom();
      NER( sGeoms[i] );
      if ( getenv("SBDBG4") ) fprintf( stderr, "SBDBG4 topology: geom %u surface %p rail0 %p\n", (unsigned)i, (void*)sGeoms[i]->GetFilletSurface(), (void*)sGeoms[i]->GetRail(0) );
    }
  for (ULONG i=0; i<n; i++)
    {
      SmFace * pFace = sFaces[i];
      const ULONG p = (i + n - 1) % n;
      SmFilletSolver * pPrevFS = sGeoms[p]->GetFilletSolver();
      SmFilletSolver * pCurrFS = sGeoms[i]->GetFilletSolver();
      NER( pPrevFS );
      NER( pCurrFS );
      SmFilletVertex * apFV[3];
      const SmFilletVertexType aeType[3] = { SM_FV_MATE, SM_FV_PRECOMPUTED, SM_FV_SETBACK };
      for (ULONG k=0; k<3; k++)
        {
          apFV[k] = new (m_pPseudoBrep) SmFilletVertex( this );
          NER( apFV[k] );
          apFV[k]->SetPointClass( SM_PC_FACE, pFace );
          apFV[k]->SetFilletVertexType( aeType[k] );
          m_vVertices.Add( apFV[k] );
        }
      SmFilletEdgeuse * pNewEU = NULL;
      SER( sGeoms[p]->MakeRailEdgeuse( m_pPseudoBrep, pPrevFS->FindIndexOfRailOnFace( pFace ), apFV[0], pNewEU ));
      SER( sGeoms[i]->MakeRailEdgeuse( m_pPseudoBrep, pCurrFS->FindIndexOfRailOnFace( pFace ), apFV[2], pNewEU ));
    }

  // the edges: per face, M-B and B-S on it, then the blend after's end section S-M
  for (ULONG i=0; i<n; i++)
    {
      SmFilletVertex * pM     = m_vVertices[3 * i];
      SmFilletVertex * pB     = m_vVertices[3 * i + 1];
      SmFilletVertex * pS     = m_vVertices[3 * i + 2];
      SmFilletVertex * pNextM = m_vVertices[3 * ((i + 1) % n)];
      SmFilletVertex * apFrom[3] = { pM, pB, pS };
      SmFilletVertex * apTo[3]   = { pB, pS, pNextM };
      for (ULONG k=0; k<3; k++)
        {
          SmFilletEdge * pEdge = NULL;
          SER( SmFilletEdge::MakeFilletEdge( m_pPseudoBrep, apFrom[k], apTo[k], pEdge, this ));
          NER( pEdge );
          if ( k < 2 )
            {
              pEdge->SetFilletEdgeType( SM_FE_PRECOMPUTED );
              pEdge->SetOriginalFace( sFaces[i] );
            }
          else
            {
              pEdge->SetFilletEdgeType( SM_FE_CROSS_SECTION );
              SmFilletEdgeuse * pPrimEU = (SmFilletEdgeuse*)pEdge->GetPrimaryEdgeuse();
              NER( pPrimEU );
              pPrimEU->SetFilletGeom( sGeoms[i] );
              pS->SetMate( 0, pNextM );
              pNextM->SetMate( 0, pS );
            }
          m_vEdges.Add( pEdge );
        }
    }

  // each blend ends on its own section; the hole is the corner's
  m_bBlending = TRUE;
  m_bBevel    = FALSE;
  return SM_SUCCESS;

} // end SmBlendCorner::MakeSetBackTopology

/*******************************************************************//**
PURPOSE: The Setback's fill (SmBlendSetback::FillHole), from the blends'
    surfaces at their end sections, and from it the base points and the
    spring curves' halves.

NOTES: Refused, its reason noted, where the fill is (SmBlendSetback's
    refusals) or where it stands past a face (Containment, certified)
    farther than the corner's tolerance.
***********************************************************************/
SmStatus SmBlendCorner::BuildSetBackFill()
{
  NER( m_pAllEUs );
  if ( m_bSetBackOverflows )
    {
      GetFilletExecutive()->GetFilletErrorInfo()->NoteError( m_cpVertex, SM_FIL_FAILURE, SM_FILERR_NO_ROOM,
          _T("Fillet Error: the rounds beside an edge at a setback corner cross past the room on that edge: the blends are larger than their faces allow") );
      return SM_ERR;
    }
  const ULONG n = m_pAllEUs->GetSize();
  if ( m_vVertices.GetSize() != 3 * n || m_vEdges.GetSize() != 3 * n )
    { SER( SM_ERR ); }
  SmFilletExecutive * pExec = GetFilletExecutive();
  NER( pExec );
  const double dTol = GetThisApproxTol3d();

  std::vector<SmFace*> sFaces;
  std::vector<SmBlendSetback::Face> sPlanes;
  std::vector<SmBlendSetback::Edge> sEdges;
  SmBoolean bPlanes = FALSE;
  SmBoolean bCurvedFaces = FALSE;
  SER( sm_SbCornerPlanes( m_cpVertex, *m_pAllEUs, pExec, dTol, m_bConvexSetBack, sFaces, sPlanes, sEdges, bPlanes, &bCurvedFaces ));
  if ( !bPlanes )
    { SER( SM_ERR ); }

  // each blend at its end section: its surface, its rails' direction on it, the rail parameter
  // of its setback vertex, and on which side the corner lies
  std::vector<SmBlendSetback::Round> sRounds( n );
  const TCHAR * pGatherWhy = NULL;
  auto gather = [&]() -> SmStatus
  {
  for (ULONG i=0; i<n; i++)
    {
      SmFilletVertex * pS = m_vVertices[3 * i + 2];
      NER( pS );
      if ( !pS->IsProcessed() )
        { SER( SM_ERR ); }
      SmTArray<SmFilletGeom*> sFVGeoms;
      pS->GetFilletGeoms( sFVGeoms );
      if ( sFVGeoms.GetSize() != 1 )
        { SER( SM_ERR ); }
      SmFilletGeom * pFG = sFVGeoms[0];
      SmFilletSolver * pFS = pFG->GetFilletSolver();
      NER( pFS );
      const SmBSplineSurface * pSurface = SM_CAST_PTR( SmBSplineSurface, pFG->GetFilletSurface() );
      NER( pSurface );
      const ULONG lRail = pFS->FindIndexOfRailOnFace( sFaces[i] );
      if ( lRail > 1 )
        { SER( SM_ERR ); }
      NER( pFG->GetRail( lRail ));
      SmCurve * pRail = pFG->GetRail( lRail )->GetCurve();
      NER( pRail );
      const SmExtent1d sRailIvl = pRail->GetNaturalInterval();
      SmBoolean bFound = FALSE;
      double dEnd = 0.0, dDist = 0.0;
      SER( pRail->DropPoint( sRailIvl, pS->GetPoint(), NULL, dTol, NULL, bFound, dEnd, dDist ));
      if ( !bFound )
        { SER( SM_ERR ); }
      SmPoint3d sOnRail, sLo, sHi;
      SER( pRail->EvaluatePoint( dEnd, sOnRail ));
      SER( pRail->EvaluatePoint( sRailIvl.GetMin(), sLo ));
      SER( pRail->EvaluatePoint( sRailIvl.GetMax(), sHi ));
      // the rails' direction: the one over the rails' interval whose iso-curve at the end
      // section passes through the rail there
      const SmExtent2d sDomain = pSurface->GetNaturalUVDomain();
      const double dParamTol = 1.0e-9 * ( 1.0 + sRailIvl.GetLength() );
      SmSurfParamType eAlong = SM_SP_UNKNOWN;
      for (int iDir=0; iDir<2 && eAlong == SM_SP_UNKNOWN; iDir++)
        {
          const SmExtent1d sAlong  = ( iDir == 0 ) ? sDomain.GetUInterval() : sDomain.GetVInterval();
          const SmExtent1d sAcross = ( iDir == 0 ) ? sDomain.GetVInterval() : sDomain.GetUInterval();
          if (   smos_Fabs( sAlong.GetMin() - sRailIvl.GetMin() ) > dParamTol
              || smos_Fabs( sAlong.GetMax() - sRailIvl.GetMax() ) > dParamTol )
            { continue; }
          SmPoint3d sA, sB;
          SER( pSurface->EvaluatePoint( iDir == 0 ? SmPoint2d( dEnd, sAcross.GetMin() ) : SmPoint2d( sAcross.GetMin(), dEnd ), sA ));
          SER( pSurface->EvaluatePoint( iDir == 0 ? SmPoint2d( dEnd, sAcross.GetMax() ) : SmPoint2d( sAcross.GetMax(), dEnd ), sB ));
          if ( smos_Min( sA.DistanceBetween( sOnRail ), sB.DistanceBetween( sOnRail )) <= dTol )
            { eAlong = ( iDir == 0 ) ? SM_SP_U : SM_SP_V; }
        }
      if ( eAlong == SM_SP_UNKNOWN )
        { pGatherWhy = _T("Fillet Error: a setback corner needs each blend's surface to run along its rails"); return SM_ERR; }
      sRounds[i].m_pSurface     = pSurface;
      sRounds[i].m_eAlong       = eAlong;
      sRounds[i].m_dEnd         = dEnd;
      sRounds[i].m_bCornerAbove = ( sHi.DistanceBetween( m_cpVertex->GetPoint() ) < sLo.DistanceBetween( m_cpVertex->GetPoint() )) ? TRUE : FALSE;
    }
  return SM_SUCCESS;
  };
  if ( gather() != SM_SUCCESS )
    { return pGatherWhy != NULL ? SetBackNotBuilt( pGatherWhy ) : SM_ERR; }

  delete m_pFill;
  m_pFill = new SmBlendSetback::Fill;
  NER( m_pFill );
  if ( SmBlendSetback::FillHole( m_crContext, sPlanes, sRounds, m_lFillOrder, dTol, *m_pFill ) != SM_SUCCESS )
    {
      if ( !bCurvedFaces || !m_pFill->m_bUncertified )
        { return SetBackNotBuilt( m_pFill->m_pReason != NULL ? m_pFill->m_pReason : _T("Fillet Error: the setback corner's hole cannot be filled") ); }
      // beside a curved face (owner, section 28): the rule's setbacks, or the one asked for, grow to
      // the least scale (to 1e-2) at which the fill is certified, with the real rounds: each trial
      // moves the setback vertices (CalcSetBackVert at the scaled setbacks, each held at its room
      // less a quarter radius) and fills the hole from the least degree up. None up to twice them:
      // refused, naming none. One asked for that certifies only grown: refused, naming the
      // setback that builds (or none, where that one does not fit with its ease).
      double dAsked = 0.0, dEase = 0.0;
      SmBoolean bOpen = FALSE;
      const SmBoolean bAsked = pExec->GetSetBackFill( m_cpVertex, dAsked, dEase, bOpen );
      if ( m_vSetBackRoom.size() != n || m_vSetBacks.size() != n )
        { SER( SM_ERR ); }
      const std::vector<double> sBase( m_vSetBacks );
      std::vector<double> sMost( n );
      for (ULONG i=0; i<n; i++)
        { sMost[i] = bAsked ? SM_BIG_DOUBLE : smos_Max( sBase[i], m_vSetBackRoom[i] - 0.25 * sEdges[i].m_dRadius ); }
      auto place = [&](double k, SmBoolean & rbPlaced) -> SmStatus
        {
          rbPlaced = TRUE;
          m_dSetBackDist = 0.0;
          for (ULONG i=0; i<n; i++)
            {
              m_vSetBacks[i] = smos_Min( k * sBase[i], sMost[i] );
              m_dSetBackDist = smos_Max( m_dSetBackDist, m_vSetBacks[i] );
            }
          for (ULONG i=0; i<n && rbPlaced; i++)
            {
              SmFilletVertex * pS = m_vVertices[3 * i + 2];
              NER( pS );
              SER( pS->CalcSetBackVert() );
              rbPlaced = pS->IsProcessed();
            }
          if ( rbPlaced && gather() != SM_SUCCESS )
            { rbPlaced = FALSE; }
          return SM_SUCCESS;
        };
      auto trial = [&](double k, SmBlendSetback::Fill *& rpFill) -> SmStatus
        {
          rpFill = NULL;
          SmBoolean bPlaced = FALSE;
          SER( place( k, bPlaced ));
          if ( getenv("SBDBG") ) { fprintf( stderr, "SBDBG grow k %.4f placed %d:", k, (int)bPlaced ); for (ULONG i=0; i<n; i++) fprintf( stderr, " %.4f (most %.4f)", m_vSetBacks[i], sMost[i] ); fprintf( stderr, "\n" ); }
          if ( !bPlaced )
            { return SM_SUCCESS; }
          SmBlendSetback::Fill * pFill = new SmBlendSetback::Fill;
          NER( pFill );
          SmBlendSetback::Probe sProbe;
          sProbe.m_bFill = true;
          if ( SmBlendSetback::FillHoleProbe( m_crContext, sPlanes, sRounds, dTol, *pFill, sProbe ) == SM_SUCCESS )
            { rpFill = pFill; }
          else
            { delete pFill; }
          return SM_SUCCESS;
        };
      SmBlendSetback::Fill * pBest = NULL;
      SER( trial( 2.0, pBest ));
      TCHAR cWhy[SM_TBLOCK_SIZE];
      if ( pBest == NULL )
        {
          if ( bAsked )
            { smos_sprintf( cWhy, _T("Fillet Error: no setback up to twice %g keeps the setback corner's fill no more concave than its faces"), dAsked ); }
          else
            { smos_sprintf( cWhy, _T("Fillet Error: no setbacks up to twice the rule's within the corner's edges keep its fill no more concave than its faces") ); }
          pExec->NoteFitSetBack( 0.0 );
          pExec->GetFilletErrorInfo()->NoteError( m_cpVertex, SM_FIL_FAILURE, SM_FILERR_NO_ROOM, cWhy );
          return SM_ERR;
        }
      double dLo = 1.0, dHi = 2.0;
      while ( dHi > 1.01 * dLo )
        {
          const double k = 0.5 * ( dLo + dHi );
          SmBlendSetback::Fill * pAt = NULL;
          SER( trial( k, pAt ));
          if ( pAt != NULL )
            {
              delete pBest;
              pBest = pAt;
              dHi = k;
            }
          else
            { dLo = k; }
        }
      SmBoolean bPlaced = FALSE;
      SER( place( dHi, bPlaced ));
      if ( !bPlaced )
        {
          delete pBest;
          SER( SM_ERR );
        }
      if ( bAsked )
        {
          delete pBest;
          SmBoolean bRoom = TRUE;
          for (ULONG i=0; i<n; i++)
            if ( dHi * dAsked + dEase >= m_vSetBackRoom[i] )
              { bRoom = FALSE; }
          const double dBuilds = bRoom ? dHi * dAsked : 0.0;
          if ( bRoom )
            { smos_sprintf( cWhy, _T("Fillet Error: a setback of %g is too short for the setback corner's fill to stay no more concave than its faces: a setback of %g builds it"), dAsked, dBuilds ); }
          else
            { smos_sprintf( cWhy, _T("Fillet Error: a setback of %g is too short for the setback corner's fill to stay no more concave than its faces, and the %g that would be does not fit on its edges"),
                            dAsked, dHi * dAsked ); }
          pExec->NoteFitSetBack( dBuilds );
          pExec->GetFilletErrorInfo()->NoteError( m_cpVertex, SM_FIL_FAILURE, SM_FILERR_NO_ROOM, cWhy );
          return SM_ERR;
        }
      delete m_pFill;
      m_pFill = pBest;
    }
  if ( m_pFill->m_vPatches.size() != 2 * n || m_pFill->m_vFace.size() != n )
    { SER( SM_ERR ); }
  double dWorst = 0.0;
  if ( !bCurvedFaces ) // WIP: containment beside a curved face not yet certified
    { SER( SmBlendSetback::Containment( *m_pFill, sPlanes, m_bConvexSetBack, dWorst )); }
  if ( dWorst > dTol )
    {
      TCHAR cWhy[SM_TBLOCK_SIZE];
      smos_sprintf( cWhy, _T("Fillet Error: the setback corner's fill stands %g past a face it meets"), dWorst );
      return SetBackNotBuilt( cWhy );
    }

  // the base points: interior part q's corner (1, 1) on its face
  for (ULONG q=0; q<n; q++)
    {
      const ULONG f = m_pFill->m_vFace[q];
      if ( f >= n )
        { SER( SM_ERR ); }
      SmBSplineSurface * pPart = m_pFill->m_vPatches[n + q];
      NER( pPart );
      const SmExtent2d sDom = pPart->GetNaturalUVDomain();
      SmPoint3d sBase;
      SER( pPart->EvaluatePoint( sDom.GetMax(), sBase ));
      SmFilletVertex * pB = m_vVertices[3 * f + 1];
      pB->SetPoint( sBase );
      pB->SetStatus( SM_FIL_PROCESSED );
    }

  // the spring curves' halves: of the setback patches' sides on the faces (v = 0, 1, u from the
  // blend's end to the base point), the one from that rail end to that base point
  for (ULONG f=0; f<n; f++)
    {
      for (ULONG k=0; k<2; k++)
        {
          SmFilletEdge * pEdge = m_vEdges[3 * f + k];
          NER( pEdge );
          SmFilletVertex * pRailEnd = m_vVertices[3 * f + ( k == 0 ? 0 : 2 )];
          SmFilletVertex * pBase    = m_vVertices[3 * f + 1];
          NER( pRailEnd );
          if ( !pRailEnd->IsProcessed() )
            { SER( SM_ERR ); }
          SmBSplineCurve * pBest = NULL;
          double dBest = SM_BIG_DOUBLE;
          for (ULONG m=0; m<n; m++)
            {
              for (ULONG e=0; e<2; e++)
                {
                  if ( m_pFill->m_vFace[( m + e ) % n] != f )
                    { continue; }
                  SmBSplineSurface * pPatch = m_pFill->m_vPatches[m];
                  NER( pPatch );
                  const SmExtent2d sDom = pPatch->GetNaturalUVDomain();
                  const double dV = ( e == 0 ) ? sDom.GetMin().y : sDom.GetMax().y;
                  SmPoint3d sEnd, sInner;
                  SER( pPatch->EvaluatePoint( SmPoint2d( sDom.GetMin().x, dV ), sEnd ));
                  SER( pPatch->EvaluatePoint( SmPoint2d( sDom.GetMax().x, dV ), sInner ));
                  const double dOff = smos_Max( sEnd.DistanceBetween( pRailEnd->GetPoint() ), sInner.DistanceBetween( pBase->GetPoint() ));
                  if ( dOff >= dBest )
                    { continue; }
                  SmBSplineCurve * pIso = NULL;
                  double dGap = 0.0;
                  SER( pPatch->CreateIsoParametricCurve( m_crContext, SM_SP_V, dV, 0.0, pIso, &sDom, &dGap ));
                  NER( pIso );
                  delete pBest;
                  pBest = pIso;
                  dBest = dOff;
                }
            }
          NER( pBest );
          if ( dBest > dTol )
            { delete pBest; SER( SM_ERR ); } // the fill's spring side ends where the blend's section meets the face
          if ( k == 1 )
            {
              // B-S: from the base point to the rail end
              SmExtent1d sNew;
              const SmStatus eStatus = pBest->ReverseParameterization( pBest->GetNaturalInterval(), sNew );
              if ( eStatus != SM_SUCCESS )
                { delete pBest; SER( eStatus ); }
            }
          pEdge->SetCurve( pBest, TRUE );
          pBest->SetOwner( pEdge );
        }
    }
  return SM_SUCCESS;

} // end SmBlendCorner::BuildSetBackFill

/*******************************************************************//**
PURPOSE: The fill's patches as faces of the fillet Brep, each bounded by
    its iso-curves: a setback patch by its end section (u = 0), its sides
    on the faces (v = 0, 1) and the two halves of its inner curve (u = 1,
    split at v = 1/2, where the spokes meet it); an interior part by its
    four sides.
***********************************************************************/
SmStatus SmBlendCorner::CreateSetBackPatches()
{
  NER( m_pFill );
  NER( m_pAllEUs );
  const ULONG n = m_pAllEUs->GetSize();
  if ( m_pFill->m_vPatches.size() != 2 * n )
    { SER( SM_ERR ); }
  SmBrep * pFilletBrep = m_pExecutive->GetFilletBrep();
  NER( pFilletBrep );
  for (ULONG p=0; p<2*n; p++)
    {
      SmBSplineSurface * pPatch = m_pFill->m_vPatches[p];
      NER( pPatch );
      const SmExtent2d sDom = pPatch->GetNaturalUVDomain();
      const double u0 = sDom.GetMin().x, u1 = sDom.GetMax().x, v0 = sDom.GetMin().y, v1 = sDom.GetMax().y;
      const double vm = 0.5 * ( v0 + v1 );
      // the loop counterclockwise in (u, v): each side its iso-curve over its natural direction
      // (the constant parameter, its value, the other's interval) and whether the loop runs it
      // backwards
      struct Side { SmSurfParamType eIso; double dAt; double dFrom; double dTo; SmBoolean bBack; };
      std::vector<Side> sSides;
      sSides.push_back( { SM_SP_V, v0, u0, u1, FALSE } );
      if ( p < n )
        {
          sSides.push_back( { SM_SP_U, u1, v0, vm, FALSE } );
          sSides.push_back( { SM_SP_U, u1, vm, v1, FALSE } );
        }
      else
        { sSides.push_back( { SM_SP_U, u1, v0, v1, FALSE } ); }
      sSides.push_back( { SM_SP_V, v1, u0, u1, TRUE } );
      sSides.push_back( { SM_SP_U, u0, v0, v1, TRUE } );

      SmTArray<SmCurve*>        s3d;
      SmTArray<SmBSplineCurve*> sUV;
      SmTArray<SmOrientType>    sOrients;
      SmObjsDelete<SmCurve*>        sDelete3d( &s3d );
      SmObjsDelete<SmBSplineCurve*> sDeleteUV( &sUV );
      for (size_t s=0; s<sSides.size(); s++)
        {
          const Side & crS = sSides[s];
          SmExtent2d sSub;
          if ( crS.eIso == SM_SP_V )
            { SER( sSub.SetMinMax( crS.dFrom, crS.dAt, crS.dTo, crS.dAt )); }
          else
            { SER( sSub.SetMinMax( crS.dAt, crS.dFrom, crS.dAt, crS.dTo )); }
          SmBSplineCurve * pIso = NULL;
          double dGap = 0.0;
          SER( pPatch->CreateIsoParametricCurve( m_crContext, crS.eIso, crS.dAt, 0.0, pIso, &sSub, &dGap ));
          NER( pIso );
          s3d.Add( pIso );
          SmLine * pLine = NULL;
          SER( SmLine::CreateLineSegment( m_crContext, 2, SmVector3d( sSub.GetMin().x, sSub.GetMin().y, 0.0 ),
                                          SmVector3d( sSub.GetMax().x, sSub.GetMax().y, 0.0 ), pLine ));
          NER( pLine );
          sUV.Add( pLine );
          sOrients.Add( crS.bBack ? SM_OT_OPPOSITE : SM_OT_SAME );
        }
      SmTArray<ULONG> sCurveLoops;
      sCurveLoops.Add( s3d.GetSize() );
      SmTArray<SmPoint3d> sLoopPoints;
      SmRegion * pNewRegion = NULL;
      SmShell  * pNewShell  = NULL;
      SmFace   * pNewFace   = NULL;
      // the face takes the surface and the curves
      SmTArray<SmCurve*>        sFace3d( s3d );
      SmTArray<SmBSplineCurve*> sFaceUV( sUV );
      s3d.ReSet();
      sUV.ReSet();
      m_pFill->m_vPatches[p] = NULL;
      SER( pFilletBrep->MakeFaceWithCurves( NULL, sCurveLoops, &sFace3d, &sFaceUV, sOrients, sLoopPoints, pPatch, sDom, SM_OT_SAME,
                                            pNewRegion, pNewShell, pNewFace ));
      NER( pNewFace );
      m_vSurfaces.Add( pPatch ); // the corner's patch faces (SmFilletCorner::GetPatchFaces, the stitch)
    }
  delete m_pFill;
  m_pFill = NULL;
  return SM_SUCCESS;

} // end SmBlendCorner::CreateSetBackPatches
