// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

/**********************************************************************//**
* FILE NAME --- SmFilletOverflow.cpp
* PURPOSE: SolidCraft change. A convex blend past the faces it rolls on: a
*   blend whose rails leave their faces across cliffs (an opening's rim, the
*   body's outline) keeps its surface, built on the faces continued, and the
*   sharp edge moves: the blend is trimmed in place where its sections meet
*   the faces beyond, along the crease each makes with it.
*
* The executive computes the fillet as usual on the support faces' surfaces
* (CreateFilletBrep, m_bCliffOverflow): m_pFilletBrep is the sheet of fillet
* and corner faces, its boundary on the surfaces of the faces the blends roll
* on and end on, or past them where the blend runs on. TrimAtCreases
*   0. cuts the fillet off at the ends its Nx1 corners left uncut, by the end
*      face's surface continued (CutFilletAtOpenEnd: the end's own event);
*   1. finds each sheet face's creases: its surface sectioned by the surface of
*      each face of the body it may cross (exactly by a plane: a line on a
*      cylinder, a conic on a cone, a circle on a torus; a certified surface
*      intersection otherwise), each kept where it lies on that face;
*   2. splits the sheet face along them and keeps the pieces inside the body;
*   3. continues a support face on its plane to a kept rail that runs past the
*      face's boundary into the material (a chamfer's leg past its wall's foot:
*      the wall carries on down to the leg; an end rib under a lip: the end face
*      carries on under the lip);
*   4. imprints the kept sheet's boundary (its rails, creases and end ribs) on
*      the body's faces, deletes the faces it encloses round the blended edges,
*      and sews the sheet in.
* Nothing is subtracted: the body's faces are only trimmed, so the blend never
* bridges an opening.
**********************************************************************/

#include "StdAfx.h"

#include <SmFilletExecutive.h>
#include <SmFilletGeom.h>
#include <SmFilletSolver.h>
#include <SmFilletStandardSolver.h>
#include <SmBrep.h>
#include <SmFace.h>
#include <SmEdge.h>
#include <SmVertex.h>
#include <SmRegion.h>
#include <SmOffsetSurface.h>
#include <SmSurface.h>
#include <SmCurve.h>
#include <SmCurveClass.h>
#include <SmPointClass.h>
#include <SmExtent3d.h>
#include <SmTol.h>
#include <SmPlane.h>
#include <SmBSplineCurve.h>

#include <algorithm>
#include <map>
#include <set>
#include <vector>

namespace
{

// The plane a planar surface lies on, by a point and its unit normal: an SmPlane, or any other
// surface SMLib finds planar (a loft's planar wall); rbPlanar FALSE when it is not planar.
SmStatus sm_PlaneOf(const SmSurface & crSurface, SmPoint3d & rAt, SmVector3d & rNormal, SmBoolean & rbPlanar)
{
  rbPlanar = FALSE;
  if ( !crSurface.IsPlanar( SM_EFF_ZERO, &rAt, &rNormal ))
    { return SM_SUCCESS; }
  if ( rNormal.Unitize() != SM_SUCCESS )
    { SER( SM_ERR ); }
  rbPlanar = TRUE;
  return SM_SUCCESS;
}

// The distance from crPnt to crFace's surface: to its plane when planar, else to the surface over
// its own domain.
SmStatus sm_DistanceToSurface(const SmFace & crFace, const SmPoint3d & crPnt, double & rdDist)
{
  rdDist = SM_BIG_DOUBLE;
  const SmSurface * pSurface = crFace.GetSurface();   NER(pSurface);
  SmPoint3d  sAt;
  SmVector3d sNormal;
  SmBoolean  bPlanar = FALSE;
  SER( sm_PlaneOf( *pSurface, sAt, sNormal, bPlanar ));
  if ( bPlanar )
    {
      rdDist = smos_Fabs( ( crPnt - sAt ).Dot( sNormal ));
      return SM_SUCCESS;
    }
  SmBoolean bDropped = FALSE, bMulti = FALSE;
  SmPoint2d sUV;
  double    dGap = SM_BIG_DOUBLE;
  SER( pSurface->DropPoint( crPnt, pSurface->GetNaturalUVDomain(), NULL, bDropped, sUV, dGap, bMulti ));
  if ( bDropped )
    { rdDist = dGap; }
  return SM_SUCCESS;
}

// Whether the edge lies on crFace's surface: its ends and three points between within dTol.
SmStatus sm_EdgeOnSurface(const SmEdge & crEdge, const SmFace & crFace, double dTol, SmBoolean & rbOn)
{
  rbOn = TRUE;
  const SmCurve * pCurve = crEdge.GetCurve();   NER(pCurve);
  const SmExtent1d sIvl = crEdge.GetInterval();
  for (ULONG i=0; i<=4 && rbOn; i++)
    {
      SmPoint3d sPnt;
      SER( pCurve->EvaluatePoint( sIvl.Evaluate( i / 4.0 ), sPnt ));
      double dDist = SM_BIG_DOUBLE;
      SER( sm_DistanceToSurface( crFace, sPnt, dDist ));
      rbOn = ( dDist <= dTol ) ? TRUE : FALSE;
    }
  return SM_SUCCESS;
}

// Owns the curves an intersection hands back.
struct sm_Curves
{
  SmTArray<SmCurve*> m_v3D, m_vUV1, m_vUV2;
  ~sm_Curves()
    {
      std::set<SmCurve*> sAll;
      for (SmTArray<SmCurve*> * pA : { &m_v3D, &m_vUV1, &m_vUV2 })
        {
          for (ULONG i=0; i<pA->GetSize(); i++)
            { sAll.insert( (*pA)[i] ); }
        }
      for (SmCurve * pC : sAll)
        { delete pC; }
    }
};

// The curves where crSheet's surface meets crFace's: a planar face's plane sections it exactly
// (SmSurface::CreatePlanarSectionCurves: lines, conics and circles on the analytic blends), any
// other surface meets it by SMLib's surface intersection over both faces' domains, its touching
// points left out. Into rCurves (3D).
SmStatus sm_Creases(const SmContext & crContext, const SmFace & crSheet, const SmFace & crFace, double dTol, sm_Curves & rCurves)
{
  const SmSurface * pSheet = crSheet.GetSurface();   NER(pSheet);
  const SmSurface * pOther = crFace.GetSurface();    NER(pOther);
  SmApproxTol3d sApprox( dTol );
  SmPoint3d  sAt;
  SmVector3d sNormal;
  SmBoolean  bPlanar = FALSE;
  SER( sm_PlaneOf( *pOther, sAt, sNormal, bPlanar ));
  if ( bPlanar )
    {
      SER( pSheet->CreatePlanarSectionCurves( crContext, crSheet.GetUVDomain(), sAt, sNormal, &sApprox, NULL, &rCurves.m_v3D, &rCurves.m_vUV1 ));
      return SM_SUCCESS;
    }
  SmBoolean abEdges[2] = { TRUE, TRUE };
  SmTArray<SmTsectCurveType> sTypes;
  sm_Curves sAll;
  SER( pSheet->GlobalSurfaceIntersect( crContext, crSheet.GetUVDomain(), *pOther, crFace.GetUVDomain(), abEdges, &sApprox, NULL,
                                       &sAll.m_v3D, &sAll.m_vUV1, &sAll.m_vUV2, &sTypes, NULL ));
  for (ULONG i=0; i<sAll.m_v3D.GetSize(); i++)
    {
      if ( i < sTypes.GetSize() && sTypes[i] == SM_TC_TOUCHING )
        { continue; }
      SmCurve * pCopy = NULL;
      SER( sAll.m_v3D[i]->Copy( crContext, pCopy ));
      NER(pCopy);
      rCurves.m_v3D.Add( pCopy );
    }
  return SM_SUCCESS;
}

// The pieces of crCurve that lie inside crFace (its runs classified in the face), as trimmed
// copies into rPieces. Runs along the face's boundary, or outside it, are left out.
SmStatus sm_PiecesInFace(const SmContext & crContext, const SmCurve & crCurve, const SmFace & crFace, double dTol,
                         SmTArray<SmCurve*> & rPieces)
{
  SmCurveClassification sClass( &crCurve, crCurve.GetNaturalInterval(), NULL, dTol );
  SER( crFace.CurveOnClassify( TRUE, sClass, TRUE ));
  for (ULONG j=0; j<sClass.GetSize(); j++)
    {
      if ( sClass[j].m_vMid.GetPointClass() != SM_PC_FACE )
        { continue; }
      SmExtent1d sRun = sClass.GetSubInterval( j );
      SmPoint3d sA, sB;
      SER( crCurve.EvaluatePoint( sRun.GetMin(), sA ));
      SER( crCurve.EvaluatePoint( sRun.GetMax(), sB ));
      if ( sA.DistanceBetween( sB ) <= dTol && !crCurve.IsClosed( crCurve.GetNaturalInterval(), dTol ) )
        { continue; }
      SmCurve * pPiece = NULL;
      SER( crCurve.Copy( crContext, pPiece ));
      NER(pPiece);
      if ( !( sRun == pPiece->GetNaturalInterval() ))
        {
          SmStatus eStat = pPiece->Trim( sRun );
          if ( eStat != SM_SUCCESS )
            { delete pPiece; SER( eStat ); }
        }
      rPieces.Add( pPiece );
    }
  return SM_SUCCESS;
}

// Whether crPnt lies in the material of crBrep (a solid region, not the infinite one or a void).
SmStatus sm_InMaterial(const SmBrep & crBrep, const SmPoint3d & crPnt, SmBoolean & rbIn, SmBoolean & rbOnBoundary)
{
  rbIn = FALSE;
  rbOnBoundary = FALSE;
  SmPointClassification sClass;
  SER( crBrep.Point3DClassify( crPnt, SmTol::GetZoneTol3d( &crBrep ), TRUE, sClass ));
  if ( sClass.GetPointClass() == SM_PC_REGION )
    {
      const SmRegion * pRegion = SM_CAST_PTR( SmRegion, sClass.GetObject() );
      rbIn = ( pRegion != NULL && !pRegion->IsInfiniteRegion() && !pRegion->IsVoid() ) ? TRUE : FALSE;
    }
  else if ( sClass.GetPointClass() != SM_PC_UNKNOWN )
    { rbOnBoundary = TRUE; }
  return SM_SUCCESS;
}

// A curve of a loop, owned until a face takes it, and how it runs.
struct sm_LoopCurve
{
  SmCurve    * m_pCurve = NULL;
  SmOrientType m_eOrient = SM_OT_SAME;
};

// A support face continued to the rails that run past its foot into the material: the loop on
// its plane, owned, and the face's edges along that loop (the material it stops at).
struct sm_Extension
{
  SmPoint3d                 m_sAt;
  SmVector3d                m_sNormal;
  std::vector<sm_LoopCurve> m_vLoop;
  std::vector<SmEdge*>      m_vAlong;
  std::vector<std::pair<SmPoint3d, SmPoint3d>> m_vAlongPieces;  // pieces of the face's edges it runs along, by their ends
  sm_Extension() = default;
  sm_Extension(const sm_Extension &) = delete;
  sm_Extension & operator=(const sm_Extension &) = delete;
  sm_Extension(sm_Extension && rOther) noexcept
    : m_sAt( rOther.m_sAt ), m_sNormal( rOther.m_sNormal ), m_vLoop( std::move( rOther.m_vLoop )), m_vAlong( std::move( rOther.m_vAlong )),
      m_vAlongPieces( std::move( rOther.m_vAlongPieces ))
    { rOther.m_vLoop.clear(); }
  ~sm_Extension() { for (sm_LoopCurve & rC : m_vLoop) { delete rC.m_pCurve; } }
};

// The two ends of a curve over its natural interval.
SmStatus sm_CurveEnds(const SmCurve & crCurve, SmPoint3d & rStart, SmPoint3d & rEnd)
{
  SmExtent1d sIvl = crCurve.GetNaturalInterval();
  SER( crCurve.EvaluatePoint( sIvl.GetMin(), rStart ));
  SER( crCurve.EvaluatePoint( sIvl.GetMax(), rEnd ));
  return SM_SUCCESS;
}

// Adds pCurve (taken over) to rLoop, run from crFrom to crTo.
SmStatus sm_AddCurve(std::vector<sm_LoopCurve> & rLoop, SmCurve * pCurve, const SmPoint3d & crFrom, const SmPoint3d & crTo)
{
  NER(pCurve);
  sm_LoopCurve sC;
  sC.m_pCurve = pCurve;
  rLoop.push_back( sC );
  SmPoint3d sStart, sEnd;
  SER( sm_CurveEnds( *pCurve, sStart, sEnd ));
  double dSame     = sStart.DistanceBetween( crFrom ) + sEnd.DistanceBetween( crTo );
  double dOpposite = sStart.DistanceBetween( crTo )   + sEnd.DistanceBetween( crFrom );
  rLoop.back().m_eOrient = ( dSame <= dOpposite ) ? SM_OT_SAME : SM_OT_OPPOSITE;
  return SM_SUCCESS;
}

// A copy of the edge's curve over the edge's interval.
SmStatus sm_EdgeCurve(const SmContext & crContext, const SmEdge & crEdge, SmCurve *& rpCurve)
{
  rpCurve = NULL;
  const SmCurve * pCurve = crEdge.GetCurve();   NER(pCurve);
  SER( pCurve->Copy( crContext, rpCurve ));     NER(rpCurve);
  SmExtent1d sIvl = crEdge.GetInterval();
  if ( !( sIvl == rpCurve->GetNaturalInterval() ))
    { SER( rpCurve->Trim( sIvl )); }
  return SM_SUCCESS;
}

// The side edge of crFace whose line holds crPnt (a straight edge of the face that is not
// blended, crPnt within dTol of the line it lies on), and its end nearer crPnt. NULL when none.
SmStatus sm_SideLineThrough(const SmFace & crFace, const SmPoint3d & crPnt, const std::set<SmEdge*> & crBlended, double dTol,
                            SmEdge *& rpSide, SmVertex *& rpFoot)
{
  rpSide = NULL;
  rpFoot = NULL;
  SmTArray<SmEdge*> sEdges;
  crFace.GetEdges( sEdges );
  for (ULONG e=0; e<sEdges.GetSize(); e++)
    {
      SmEdge * pEdge = sEdges[e];
      if ( crBlended.count( pEdge ) || pEdge->GetCurve() == NULL || !pEdge->GetCurve()->IsLinear( SM_EFF_ZERO ))
        { continue; }
      const SmPoint3d sA = pEdge->GetStartVertex()->GetPoint(), sB = pEdge->GetEndVertex()->GetPoint();
      SmVector3d sDir = sB - sA;
      const double dLen = sDir.Length();
      if ( dLen <= dTol )
        { continue; }
      sDir = sDir / dLen;
      const SmVector3d sTo = crPnt - sA;
      const SmVector3d sOff = sTo - sDir * sTo.Dot( sDir );
      if ( sOff.Length() > dTol )
        { continue; }
      if ( rpSide != NULL )
        { SER( SM_ERR ); }  // two side lines through one point
      rpSide = pEdge;
      rpFoot = ( crPnt.DistanceBetween( sA ) <= crPnt.DistanceBetween( sB )) ? pEdge->GetStartVertex() : pEdge->GetEndVertex();
    }
  return SM_SUCCESS;
}

// A run of a kept rail off its support face, in the material: its curve (owned), and at each end
// the sheet vertex it ends at (where the rail ends there) or the face's edge it crosses there.
struct sm_OffRun
{
  SmCurve  * m_pCurve = NULL;
  SmPoint3d  m_asEnd[2];
  SmVertex * m_apVertex[2] = { NULL, NULL };
  SmEdge   * m_apCross[2]  = { NULL, NULL };
  sm_OffRun() = default;
  sm_OffRun(const sm_OffRun &) = delete;
  sm_OffRun & operator=(const sm_OffRun &) = delete;
  sm_OffRun(sm_OffRun && rOther) noexcept
    : m_pCurve( rOther.m_pCurve )
    {
      for (ULONG k=0; k<2; k++) { m_asEnd[k] = rOther.m_asEnd[k]; m_apVertex[k] = rOther.m_apVertex[k]; m_apCross[k] = rOther.m_apCross[k]; }
      rOther.m_pCurve = NULL;
    }
  ~sm_OffRun() { delete m_pCurve; }
};

// A copy of crEdge's curve from crFrom to crTo (both on it).
SmStatus sm_EdgePiece(const SmContext & crContext, const SmEdge & crEdge, const SmPoint3d & crFrom, const SmPoint3d & crTo, double dTol,
                      SmCurve *& rpPiece)
{
  rpPiece = NULL;
  SmCurve * pCurve = NULL;
  SER( sm_EdgeCurve( crContext, crEdge, pCurve ));
  SmObjDelete sClean( pCurve );
  double adT[2] = { 0.0, 0.0 };
  const SmPoint3d asAt[2] = { crFrom, crTo };
  for (ULONG k=0; k<2; k++)
    {
      SmBoolean bOn = FALSE;
      double    dDist = SM_BIG_DOUBLE;
      SER( pCurve->DropPoint( pCurve->GetNaturalInterval(), asAt[k], NULL, SM_BIG_DOUBLE, NULL, bOn, adT[k], dDist ));
      if ( !bOn || dDist > 10.0 * dTol )
        { SER( SM_ERR ); }
    }
  SmExtent1d sIvl;
  SER( sIvl.SetMinMax( smos_Min( adT[0], adT[1] ), smos_Max( adT[0], adT[1] )));
  if ( sIvl.GetLength() <= SM_EFF_ZERO )
    { return SM_SUCCESS; }
  if ( !( sIvl == pCurve->GetNaturalInterval() ))
    { SER( pCurve->Trim( sIvl )); }
  sClean.Clear();
  rpPiece = pCurve;
  return SM_SUCCESS;
}

// The pieces of crEdge from crFrom to crTo (both on it), each way along it that does not leave
// it: the one between them, and round through its vertex as well where it closes on itself.
SmStatus sm_EdgeArcs(const SmContext & crContext, const SmEdge & crEdge, const SmPoint3d & crFrom, const SmPoint3d & crTo, double dTol,
                     std::vector<std::vector<SmCurve*>> & rArcs)
{
  SmCurve * pCurve = NULL;
  SER( sm_EdgeCurve( crContext, crEdge, pCurve ));
  SmObjDelete sClean( pCurve );
  const SmExtent1d sAll = pCurve->GetNaturalInterval();
  double adT[2] = { 0.0, 0.0 };
  const SmPoint3d asAt[2] = { crFrom, crTo };
  for (ULONG k=0; k<2; k++)
    {
      SmBoolean bOn = FALSE;
      double    dDist = SM_BIG_DOUBLE;
      SER( pCurve->DropPoint( sAll, asAt[k], NULL, SM_BIG_DOUBLE, NULL, bOn, adT[k], dDist ));
      if ( !bOn || dDist > 10.0 * dTol )
        { SER( SM_ERR ); }
    }
  const double dLo = smos_Min( adT[0], adT[1] ), dHi = smos_Max( adT[0], adT[1] );
  auto fPiece = [&]( double dA, double dB, std::vector<SmCurve*> & rArc ) -> SmStatus
    {
      if ( dB - dA <= SM_EFF_ZERO * smos_Max( 1.0, sAll.GetLength() ))
        { return SM_SUCCESS; }
      SmCurve * pPiece = NULL;
      SER( pCurve->Copy( crContext, pPiece ));
      NER(pPiece);
      SmExtent1d sIvl( dA, dB );
      if ( !( sIvl == pPiece->GetNaturalInterval() ))
        {
          SmStatus eStat = pPiece->Trim( sIvl );
          if ( eStat != SM_SUCCESS )
            { delete pPiece; SER( eStat ); }
        }
      rArc.push_back( pPiece );
      return SM_SUCCESS;
    };
  // as the arc runs from crFrom: from adT[0] to adT[1]
  std::vector<SmCurve*> vBetween;
  SER( fPiece( dLo, dHi, vBetween ));
  if ( adT[0] > adT[1] )
    { std::reverse( vBetween.begin(), vBetween.end() ); }
  rArcs.push_back( vBetween );
  if ( crEdge.GetStartVertex() == crEdge.GetEndVertex() )
    {
      std::vector<SmCurve*> vRound, vToEnd, vFromStart;
      SER( fPiece( dHi, sAll.GetMax(), vToEnd ));
      SER( fPiece( sAll.GetMin(), dLo, vFromStart ));
      // from adT[0] round through the vertex to adT[1]
      if ( adT[0] >= adT[1] ) { vRound = vToEnd; vRound.insert( vRound.end(), vFromStart.begin(), vFromStart.end() ); }
      else                    { vRound = vFromStart; vRound.insert( vRound.end(), vToEnd.begin(), vToEnd.end() ); }
      rArcs.push_back( vRound );
    }
  return SM_SUCCESS;
}

// A point of a loop on the plane (crAt, crX, crY): its sampled polygon, in order.
SmStatus sm_LoopPolygon(const std::vector<sm_LoopCurve> & crLoop, const SmPoint3d & crAt, const SmVector3d & crX, const SmVector3d & crY,
                        std::vector<SmPoint2d> & rUV)
{
  for (const sm_LoopCurve & rC : crLoop)
    {
      SmExtent1d sIvl = rC.m_pCurve->GetNaturalInterval();
      for (ULONG i=0; i<16; i++)
        {
          const double dT = ( rC.m_eOrient == SM_OT_SAME ) ? i / 16.0 : 1.0 - i / 16.0;
          SmPoint3d sPnt;
          SER( rC.m_pCurve->EvaluatePoint( sIvl.Evaluate( dT ), sPnt ));
          const SmVector3d sD = sPnt - crAt;
          rUV.push_back( SmPoint2d( sD.Dot( crX ), sD.Dot( crY )));
        }
    }
  return SM_SUCCESS;
}

// Axes in the plane square to crNormal.
SmStatus sm_PlaneAxes(const SmVector3d & crNormal, SmVector3d & rX, SmVector3d & rY)
{
  SmVector3d sAxis( 1.0, 0.0, 0.0 );
  if ( smos_Fabs( crNormal.y ) < smos_Fabs( crNormal.x ) && smos_Fabs( crNormal.y ) <= smos_Fabs( crNormal.z ))
    { sAxis = SmVector3d( 0.0, 1.0, 0.0 ); }
  else if ( smos_Fabs( crNormal.z ) < smos_Fabs( crNormal.x ))
    { sAxis = SmVector3d( 0.0, 0.0, 1.0 ); }
  rX = crNormal * sAxis;
  if ( rX.Unitize() != SM_SUCCESS )
    { SER( SM_ERR ); }
  rY = crNormal * rX;
  return SM_SUCCESS;
}

/* The extension of the planar support face crFace to the kept rails on its plane that run off it
   into the material: a chamfer's leg past its wall's foot (the wall carries on down to the leg),
   an end rib under a lip standing on the end face (the end face carries on under the lip). Each
   chain of such runs ends where it crosses the face's boundary, or at a vertex of the sheet on the
   line of a side edge of the face (the walls beside it meet it there), joined to that edge's end
   along its line. The loop runs along the chain and back along the face's boundary, off the
   blended and side edges, the way round whose inside, beside the chain, lies in the material.
   rcWhy is set where the extension is not one this builds. */
SmStatus sm_Extend(const SmContext & crContext, const SmBrep & crTarget, const SmFace & crFace, const std::vector<SmEdge*> & crRails,
                   const std::set<SmEdge*> & crBlended, double dTol, std::vector<sm_Extension> & rExtensions, const TCHAR *& rcWhy)
{
  rcWhy = NULL;
  // the runs of the rails off the face, in the material
  std::vector<sm_OffRun> vRuns;
  for (SmEdge * pRail : crRails)
    {
      SmCurve * pCurve = NULL;
      SER( sm_EdgeCurve( crContext, *pRail, pCurve ));
      SmObjDelete sClean( pCurve );
      const SmExtent1d sAll = pCurve->GetNaturalInterval();
      SmCurveClassification sClass( pCurve, sAll, NULL, dTol );
      SER( crFace.CurveOnClassify( TRUE, sClass, TRUE ));
      for (ULONG j=0; j<sClass.GetSize(); j++)
        {
          if ( sClass[j].m_vMid.GetPointClass() != SM_PC_UNKNOWN )
            { continue; }
          SmExtent1d sRun = sClass.GetSubInterval( j );
          SmPoint3d sA, sB, sMid;
          SER( pCurve->EvaluatePoint( sRun.GetMin(), sA ));
          SER( pCurve->EvaluatePoint( sRun.GetMax(), sB ));
          if ( sA.DistanceBetween( sB ) <= dTol )
            { continue; }
          SER( pCurve->EvaluatePoint( sRun.Evaluate( 0.5 ), sMid ));
          SmBoolean bIn = FALSE, bOnBoundary = FALSE;
          SER( sm_InMaterial( crTarget, sMid, bIn, bOnBoundary ));
          if ( !bIn )
            { continue; }  // in the open: cut off along a crease elsewhere
          sm_OffRun sRunOff;
          SER( pCurve->Copy( crContext, sRunOff.m_pCurve ));
          NER(sRunOff.m_pCurve);
          if ( !( sRun == sRunOff.m_pCurve->GetNaturalInterval() ))
            { SER( sRunOff.m_pCurve->Trim( sRun )); }
          sRunOff.m_asEnd[0] = sA;
          sRunOff.m_asEnd[1] = sB;
          const SmPointClassification * apPC[2] = { &sClass[j].m_vStart, &sClass[j].m_vEnd };
          const double adT[2] = { sRun.GetMin(), sRun.GetMax() };
          for (ULONG k=0; k<2; k++)
            {
              SmVertex * pRailEnd = ( smos_Fabs( adT[k] - sAll.GetMin() ) <= SM_EFF_ZERO_SQRT * smos_Max( 1.0, sAll.GetLength() ))
                                    ? pRail->GetStartVertex()
                                    : ( smos_Fabs( adT[k] - sAll.GetMax() ) <= SM_EFF_ZERO_SQRT * smos_Max( 1.0, sAll.GetLength() ))
                                      ? pRail->GetEndVertex() : NULL;
              if ( pRailEnd != NULL && pRailEnd->GetPoint().DistanceBetween( sRunOff.m_asEnd[k] ) <= 10.0 * dTol )
                { sRunOff.m_apVertex[k] = pRailEnd; }
              else if ( apPC[k]->GetPointClass() == SM_PC_EDGE )
                { sRunOff.m_apCross[k] = SM_CAST_PTR( SmEdge, apPC[k]->GetObject() ); }
              if ( sRunOff.m_apVertex[k] == NULL && sRunOff.m_apCross[k] == NULL )
                {
                  rcWhy = _T("Fillet Error: a blend's rail runs off a face into the material across a vertex of the face");
                  return SM_ERR;
                }
            }
          vRuns.push_back( std::move( sRunOff ));
        }
    }
  if ( vRuns.empty() )
    { return SM_SUCCESS; }

  SmPoint3d  sAt;
  SmVector3d sNormal;
  SmBoolean  bPlanar = FALSE;
  SER( sm_PlaneOf( *crFace.GetSurface(), sAt, sNormal, bPlanar ));
  if ( !bPlanar )
    {
      rcWhy = _T("Fillet Error: a blend's rail runs past the boundary of a face that is not planar into the material");
      return SM_ERR;
    }
  SmVector3d sX, sY;
  SER( sm_PlaneAxes( sNormal, sX, sY ));
  SmTArray<SmEdge*> sFaceEdges;
  crFace.GetEdges( sFaceEdges );

  // the chains of runs, joined at the sheet's vertices; each run once
  std::vector<bool> vUsed( vRuns.size(), false );
  for (size_t r0=0; r0<vRuns.size(); r0++)
    {
      if ( vUsed[r0] )
        { continue; }
      // the chain as (run, reversed) pairs, its ends the first run's start and the last run's end
      std::vector<std::pair<size_t, bool>> vChain( 1, std::make_pair( r0, false ));
      vUsed[r0] = true;
      for (ULONG lSide=0; lSide<2; lSide++)
        {
          for (;;)
            {
              const std::pair<size_t, bool> & rAt = ( lSide == 0 ) ? vChain.back() : vChain.front();
              const sm_OffRun & rRun = vRuns[rAt.first];
              const ULONG lEnd = ( lSide == 0 ) ? ( rAt.second ? 0 : 1 ) : ( rAt.second ? 1 : 0 );
              SmVertex * pV = rRun.m_apVertex[lEnd];
              if ( pV == NULL )
                { break; }
              size_t lNext = vRuns.size();
              bool   bRev  = false;
              for (size_t q=0; q<vRuns.size(); q++)
                {
                  if ( vUsed[q] )
                    { continue; }
                  if ( vRuns[q].m_apVertex[0] == pV ) { lNext = q; bRev = ( lSide != 0 ); break; }
                  if ( vRuns[q].m_apVertex[1] == pV ) { lNext = q; bRev = ( lSide == 0 ); break; }
                }
              if ( lNext == vRuns.size() )
                { break; }
              vUsed[lNext] = true;
              if ( lSide == 0 ) { vChain.push_back( std::make_pair( lNext, bRev )); }
              else              { vChain.insert( vChain.begin(), std::make_pair( lNext, bRev )); }
            }
        }
      // the chain's two ends: [0] its start, [1] its end
      SmPoint3d  asEnd[2];
      SmVertex * apV[2];
      SmEdge   * apCross[2];
      {
        const sm_OffRun & rFirst = vRuns[vChain.front().first];
        const sm_OffRun & rLast  = vRuns[vChain.back().first];
        const ULONG lF = vChain.front().second ? 1 : 0, lL = vChain.back().second ? 0 : 1;
        asEnd[0] = rFirst.m_asEnd[lF]; apV[0] = rFirst.m_apVertex[lF]; apCross[0] = rFirst.m_apCross[lF];
        asEnd[1] = rLast.m_asEnd[lL];  apV[1] = rLast.m_apVertex[lL];  apCross[1] = rLast.m_apCross[lL];
      }
      if ( apV[0] != NULL && apV[0] == apV[1] )
        {
          rcWhy = _T("Fillet Error: a blend's rail runs all round past a face's boundary into the material");
          return SM_ERR;
        }
      // where each end meets the face's boundary: its crossing, or its side line's end on the face
      SmEdge   * apSide[2] = { NULL, NULL };
      SmVertex * apFoot[2] = { NULL, NULL };
      for (ULONG k=0; k<2; k++)
        {
          if ( apCross[k] != NULL )
            { continue; }
          SER( sm_SideLineThrough( crFace, asEnd[k], crBlended, 10.0 * dTol, apSide[k], apFoot[k] ));
          if ( apSide[k] == NULL )
            {
              rcWhy = _T("Fillet Error: a blend's rail past a face's boundary ends off the lines of the face's side edges");
              return SM_ERR;
            }
        }
      // The ways round the face's boundary from the chain's end to its start: from a foot vertex,
      // or from a crossing to either end of the edge crossed, on along edges that are neither
      // blended nor side edges, to the start's foot or crossed edge.
      struct sm_Way { std::vector<SmEdge*> m_vEdges; std::vector<SmVertex*> m_vVerts; };
      std::vector<sm_Way> vWays;
      std::vector<SmVertex*> vStarts;
      if ( apFoot[1] != NULL ) { vStarts.push_back( apFoot[1] ); }
      else { vStarts.push_back( apCross[1]->GetStartVertex() ); vStarts.push_back( apCross[1]->GetEndVertex() ); }
      for (SmVertex * pStart : vStarts)
        {
          // a start edge's own direction, from a foot: each edge at it
          std::vector<SmEdge*> vFirst;
          if ( apFoot[1] != NULL )
            {
              for (ULONG e=0; e<sFaceEdges.GetSize(); e++)
                {
                  SmEdge * pE = sFaceEdges[e];
                  if ( pE != apSide[0] && pE != apSide[1] && !crBlended.count( pE ) && ( pE->GetStartVertex() == pStart || pE->GetEndVertex() == pStart ))
                    { vFirst.push_back( pE ); }
                }
            }
          else
            { vFirst.push_back( NULL ); }  // leaving the crossed edge at pStart
          for (SmEdge * pFirst : vFirst)
            {
              sm_Way sWay;
              sWay.m_vVerts.push_back( pStart );
              SmVertex * pAt = pStart;
              SmEdge   * pFrom = ( pFirst == NULL ) ? apCross[1] : NULL;
              SmEdge   * pNext = pFirst;
              SmBoolean  bDone = FALSE;
              for (ULONG lSteps=0; lSteps<=sFaceEdges.GetSize() && !bDone; lSteps++)
                {
                  // arrived: at the start's foot, or at an end of the edge the start crosses
                  if ( pNext == NULL && (( apFoot[0] != NULL && pAt == apFoot[0] ) ||
                                         ( apCross[0] != NULL && pFrom != NULL && pFrom != apCross[0] &&
                                           ( apCross[0]->GetStartVertex() == pAt || apCross[0]->GetEndVertex() == pAt ))))
                    { bDone = TRUE; break; }
                  if ( pNext == NULL )
                    {
                      for (ULONG e=0; e<sFaceEdges.GetSize(); e++)
                        {
                          SmEdge * pE = sFaceEdges[e];
                          if ( pE == pFrom || pE == apSide[0] || pE == apSide[1] || crBlended.count( pE ) || pE == apCross[1]
                               || ( pE->GetStartVertex() != pAt && pE->GetEndVertex() != pAt ))
                            { continue; }
                          if ( pNext != NULL ) { pNext = NULL; break; }  // a branch: no one way
                          pNext = pE;
                        }
                      if ( pNext == NULL )
                        { break; }
                    }
                  if ( pNext == apCross[0] && apCross[0] != NULL )
                    { pFrom = pNext; bDone = TRUE; break; }  // the start's crossing lies on this edge
                  sWay.m_vEdges.push_back( pNext );
                  pAt = pNext->GetOtherVertex( pAt );
                  sWay.m_vVerts.push_back( pAt );
                  pFrom = pNext;
                  pNext = NULL;
                }
              if ( bDone )
                { vWays.push_back( sWay ); }
            }
        }
      // whether a candidate loop's inside, beside the chain's middle, lies in the material, on the
      // side of the chain the blend takes off (toward the blended edge's nearest point, seen in
      // the plane: the material may lie on both sides, a lip over the end)
      auto fTakesOff = [&]( const sm_Extension & rExt, SmBoolean & rbIn ) -> SmStatus
        {
          std::vector<SmPoint2d> vUV;
          SER( sm_LoopPolygon( rExt.m_vLoop, sAt, sX, sY, vUV ));
          double dArea = 0.0;
          for (size_t p=0; p<vUV.size(); p++)
            {
              const SmPoint2d & rA = vUV[p];
              const SmPoint2d & rB = vUV[(p+1) % vUV.size()];
              dArea += rA.x * rB.y - rB.x * rA.y;
            }
          const sm_OffRun & rMidRun = vRuns[vChain.front().first];
          SmExtent1d sIvl = rMidRun.m_pCurve->GetNaturalInterval();
          SmVector3d asD[2];
          SER( rMidRun.m_pCurve->Evaluate( sIvl.Evaluate( 0.5 ), 1, TRUE, asD, FALSE ));
          SmVector3d sTan = vChain.front().second ? -asD[1] : asD[1];
          const double dLen = rMidRun.m_asEnd[0].DistanceBetween( rMidRun.m_asEnd[1] );
          SmVector3d sIn = sNormal * sTan;  // left of the run in the plane, as the loop runs
          if ( sIn.Unitize() != SM_SUCCESS )
            { SER( SM_ERR ); }
          const double dSide = ( dArea > 0.0 ) ? 1.0 : -1.0;
          SmPoint3d sProbe = SmPoint3d( asD[0] ) + sIn * ( dSide * 0.01 * dLen );
          SmBoolean bIn = FALSE, bOnBoundary = FALSE;
          SER( sm_InMaterial( crTarget, sProbe, bIn, bOnBoundary ));
          double dNearest = SM_BIG_DOUBLE;
          SmPoint3d sToward;
          for (SmEdge * pBlended : crBlended)
            {
              SmBoolean bOnCurve = FALSE;
              double    dParam = 0.0, dDist = SM_BIG_DOUBLE;
              SER( pBlended->GetCurve()->DropPoint( pBlended->GetInterval(), SmPoint3d( asD[0] ), NULL, SM_BIG_DOUBLE, NULL, bOnCurve, dParam, dDist ));
              if ( !bOnCurve || dDist >= dNearest )
                { continue; }
              dNearest = dDist;
              SER( pBlended->GetCurve()->EvaluatePoint( dParam, sToward ));
            }
          if ( dNearest == SM_BIG_DOUBLE )
            { SER( SM_ERR ); }
          SmVector3d sOff = sToward - SmPoint3d( asD[0] );
          sOff = sOff - sNormal * sOff.Dot( sNormal );
          if ( sOff.Dot( sIn ) * dSide <= 0.0 )
            { bIn = FALSE; }
          rbIn = bIn;
          return SM_SUCCESS;
        };
      // the chain's curves, as the loop runs, into a new candidate
      auto fChain = [&]( sm_Extension & rExt ) -> SmStatus
        {
          rExt.m_sAt     = sAt;
          rExt.m_sNormal = sNormal;
          for (const std::pair<size_t, bool> & rLink : vChain)
            {
              const sm_OffRun & rRun = vRuns[rLink.first];
              SmCurve * pCopy = NULL;
              SER( rRun.m_pCurve->Copy( crContext, pCopy ));
              const ULONG lFrom = rLink.second ? 1 : 0;
              SER( sm_AddCurve( rExt.m_vLoop, pCopy, rRun.m_asEnd[lFrom], rRun.m_asEnd[1 - lFrom] ));
            }
          return SM_SUCCESS;
        };
      std::vector<sm_Extension> vCandidates;
      // Both ends crossing one edge: back along that edge between them, either way round it
      // when it closes on itself (a pin's rim under a bar's bottom), or the one way along it.
      if ( apCross[0] != NULL && apCross[0] == apCross[1] )
        {
          std::vector<std::vector<SmCurve*>> vArcs;
          SER( sm_EdgeArcs( crContext, *apCross[1], asEnd[1], asEnd[0], dTol, vArcs ));
          for (std::vector<SmCurve*> & rArc : vArcs)
            {
              sm_Extension sExt;
              SER( fChain( sExt ));
              SmPoint3d sFrom = asEnd[1];
              for (SmCurve * pPiece : rArc)
                {
                  SmPoint3d sA, sB;
                  SER( sm_CurveEnds( *pPiece, sA, sB ));
                  const SmPoint3d sTo = ( sA.DistanceBetween( sFrom ) <= sB.DistanceBetween( sFrom )) ? sB : sA;
                  SER( sm_AddCurve( sExt.m_vLoop, pPiece, sFrom, sTo ));
                  sExt.m_vAlongPieces.push_back( std::make_pair( sFrom, sTo ));
                  sFrom = sTo;
                }
              rArc.clear();
              SmBoolean bIn = FALSE;
              SER( fTakesOff( sExt, bIn ));
              if ( bIn )
                { vCandidates.push_back( std::move( sExt )); }
            }
        }
      // each way's loop; the one whose inside beside the chain lies in the material
      for (const sm_Way & rWay : vWays)
        {
          sm_Extension sExt;
          SER( fChain( sExt ));
          // the end to the way's first vertex
          SmPoint3d sHere = asEnd[1];
          const SmPoint3d sFirst = rWay.m_vVerts.front()->GetPoint();
          if ( sHere.DistanceBetween( sFirst ) > dTol )
            {
              SmCurve * pPiece = NULL;
              if ( apCross[1] != NULL )
                { SER( sm_EdgePiece( crContext, *apCross[1], sHere, sFirst, dTol, pPiece )); }
              else
                {
                  SmBSplineCurve * pLine = NULL;
                  SER( SmBSplineCurve::CreateLineSegment( crContext, 3, sHere, sFirst, pLine ));
                  pPiece = pLine;
                }
              if ( pPiece != NULL )
                { SER( sm_AddCurve( sExt.m_vLoop, pPiece, sHere, sFirst )); }
              if ( apCross[1] != NULL )
                { sExt.m_vAlongPieces.push_back( std::make_pair( sHere, sFirst )); }
            }
          for (size_t e=0; e<rWay.m_vEdges.size(); e++)
            {
              SmCurve * pCurve = NULL;
              SER( sm_EdgeCurve( crContext, *rWay.m_vEdges[e], pCurve ));
              SER( sm_AddCurve( sExt.m_vLoop, pCurve, rWay.m_vVerts[e]->GetPoint(), rWay.m_vVerts[e+1]->GetPoint() ));
              sExt.m_vAlong.push_back( rWay.m_vEdges[e] );
            }
          // the way's last vertex to the start
          const SmPoint3d sLast = rWay.m_vVerts.back()->GetPoint();
          if ( sLast.DistanceBetween( asEnd[0] ) > dTol )
            {
              SmCurve * pPiece = NULL;
              if ( apCross[0] != NULL )
                { SER( sm_EdgePiece( crContext, *apCross[0], sLast, asEnd[0], dTol, pPiece )); }
              else
                {
                  SmBSplineCurve * pLine = NULL;
                  SER( SmBSplineCurve::CreateLineSegment( crContext, 3, sLast, asEnd[0], pLine ));
                  pPiece = pLine;
                }
              if ( pPiece != NULL )
                { SER( sm_AddCurve( sExt.m_vLoop, pPiece, sLast, asEnd[0] )); }
              if ( apCross[0] != NULL )
                { sExt.m_vAlongPieces.push_back( std::make_pair( sLast, asEnd[0] )); }
            }
          SmBoolean bIn = FALSE;
          SER( fTakesOff( sExt, bIn ));
          if ( bIn )
            { vCandidates.push_back( std::move( sExt )); }
        }
      if ( vCandidates.size() != 1 )
        {
          rcWhy = _T("Fillet Error: a blend's rail past a face's boundary does not close on the face in one way");
          return SM_ERR;
        }
      rExtensions.push_back( std::move( vCandidates[0] ));
    }
  return SM_SUCCESS;
}

// The face on the plane through crAt square to crNormal bounded by rExt's loop (whose curves
// it takes), in crBrep's infinite region: the plane over the loop's box, the loop anticlockwise.
SmStatus sm_MakeExtensionFace(const SmContext & crContext, SmBrep & rBrep, sm_Extension & rExt, double dTol, SmFace *& rpFace)
{
  rpFace = NULL;
  SmVector3d sX, sY;
  SER( sm_PlaneAxes( rExt.m_sNormal, sX, sY ));
  // the loop's points in the plane's coordinates, in its order: its box and signed area
  std::vector<SmPoint2d> vUV;
  SER( sm_LoopPolygon( rExt.m_vLoop, rExt.m_sAt, sX, sY, vUV ));
  if ( vUV.empty() )
    { SER( SM_ERR ); }
  SmExtent2d sBox( vUV[0], vUV[0] );
  double dArea = 0.0;
  for (size_t p=0; p<vUV.size(); p++)
    {
      sBox.AddPoint2d( vUV[p] );
      const SmPoint2d & rA = vUV[p];
      const SmPoint2d & rB = vUV[(p+1) % vUV.size()];
      dArea += rA.x * rB.y - rB.x * rA.y;
    }
  if ( dArea < 0.0 )
    {
      std::reverse( rExt.m_vLoop.begin(), rExt.m_vLoop.end() );
      for (sm_LoopCurve & rC : rExt.m_vLoop)
        { rC.m_eOrient = ( rC.m_eOrient == SM_OT_SAME ) ? SM_OT_OPPOSITE : SM_OT_SAME; }
    }
  const double dMargin = 0.01 * sBox.GetSize().Length() + 10.0 * dTol;
  SmExtent2d sDomain( sBox.GetMin().x - dMargin, sBox.GetMin().y - dMargin, sBox.GetMax().x + dMargin, sBox.GetMax().y + dMargin );
  SmPlane * pPlane = new (crContext) SmPlane( rExt.m_sAt, sX, sY, SmVector2d( 1.0, 1.0 ), sDomain, &crContext );
  NER(pPlane);
  SmTArray<ULONG>        sCounts;
  SmTArray<SmCurve*>     sCurves;
  SmTArray<SmOrientType> sOrients;
  SmTArray<SmPoint3d>    sNoPoints;
  sCounts.Add( (ULONG)rExt.m_vLoop.size() );
  for (sm_LoopCurve & rC : rExt.m_vLoop)
    {
      sCurves.Add( rC.m_pCurve );
      sOrients.Add( rC.m_eOrient );
      rC.m_pCurve = NULL;  // the face takes them
    }
  SmRegion * pNewRegion = NULL;
  SmShell  * pNewShell  = NULL;
  SER( rBrep.MakeFaceWithCurves( rBrep.GetInfiniteRegion(), sCounts, &sCurves, NULL, sOrients, sNoPoints, pPlane,
                                 pPlane->GetNaturalUVDomain(), SM_OT_SAME, pNewRegion, pNewShell, rpFace ));
  NER(rpFace);
  return SM_SUCCESS;
}

} // namespace

/*******************************************************************//**
PURPOSE: SolidCraft change. Whether a variable round on planar faces holds
    its radius law, between the filleted edge's ends.

NOTES: For TrimAtCreases. A rolling ball of radius r between two
    planes meeting at the interior angle a has its centre r / sin(a/2) off
    the edge on the wedge's bisector. The round touches it along the circle
    square to the centres' line, its centre moved along that line by
    -r dr/ds (s the centres' arc length) and its radius r sqrt(1 - (dr/ds)^2)
    (a constant radius: the great circle); that circle's point in the
    bisector plane, on the edge's side, lies on the round. The fillet
    surface is traced on the faces' continuations; a trace that went wrong
    past them (clamped, or run into the law's continuation where it shrinks
    to nothing) gives a surface that does not, though it is a valid face
    (SmFilletExecutive.cpp, sm_ExtendTraceToItsRails). That point must lie
    on the surface, at seven points along the edge, to twenty times the
    solver's tolerance and a ten-thousandth of the radius; of the two angles the planes make, and the bisector
    either way, the one that fits. Rounds of other solvers, on faces that
    are not both planes, or along an edge that is not straight, pass.
***********************************************************************/
static SmStatus sm_HoldsItsLaw
  (SmFilletSolver * pSolver,   // in : the solver
   SmBoolean      & rbHolds)   // out: FALSE = the round's surface misses its law
{
  rbHolds = TRUE;
  SmVariableRadiusFS * pVRFS = SM_CAST_PTR( SmVariableRadiusFS, pSolver );
  if ( pVRFS == NULL )
    { return SM_SUCCESS; }
  SmFilletSurfaceGenerator * pFSG = pSolver->GetFilletSurfaceGenerator();
  if ( pFSG == NULL || pFSG->GetFilletSurfaceGeneratorType() != SM_FSG_CIRCULAR )
    { return SM_SUCCESS; }
  SmVector3d asNormal[2];
  for (ULONG k=0; k<2; k++)
    {
      SmOffsetSurface * pOffset = pSolver->GetSurface(k);                 NER(pOffset);
      const SmSurface * pBase   = pOffset->GetBaseSurface();              NER(pBase);
      SmPoint3d sAt;
      if ( !pBase->IsPlanar( SM_EFF_ZERO, &sAt, &asNormal[k] ))
        { return SM_SUCCESS; }
    }
  double dPhi = acos( smos_Max( -1.0, smos_Min( 1.0, asNormal[0].Dot( asNormal[1] ))));
  SmTArray<SmFilletGeom*> sFGs;
  pSolver->GetFilletGeoms( sFGs );
  if ( sFGs.GetSize() != 1 || sFGs[0] == NULL || sFGs[0]->GetFilletSurface() == NULL )
    { return SM_SUCCESS; }
  const SmSurface * pFillet = sFGs[0]->GetFilletSurface();
  SmEdge * pEdge = pSolver->GetEdgeuse(0)->GetEdge();                     NER(pEdge);
  NER(pEdge->GetCurve());
  if ( !pEdge->GetCurve()->IsLinear( SM_EFF_ZERO ))
    { return SM_SUCCESS; }
  double dTol = 20.0 * pSolver->GetThisApproxTol3d();
  for (ULONG i=1; i<8 && rbHolds; i++)
    {
      double     dT = pEdge->GetInterval().Evaluate( i / 8.0 );
      SmVector3d asEdge[2];
      SER( pEdge->GetCurve()->Evaluate( dT, 1, TRUE, asEdge, FALSE ));
      SmPoint3d  sPnt( asEdge[0] );
      double     adLaw[3];
      SER( pVRFS->GetFilletLaw()->Evaluate( dT, pVRFS->GetFilletLawInterval(), pVRFS->GetFilletLawOrientation(), adLaw ));
      double dRadius = smos_Fabs( adLaw[0] );
      double dSlope  = ( adLaw[0] < 0.0 ) ? -adLaw[1] : adLaw[1]; // d|r|/dT
      // The ball's centre on the wedge's bisector, r / sin(a/2) off the edge, and the touch
      // circle's point on the bisector plane toward the edge. The planes' normals are both
      // outward (or inward), the bisector along their sum and the angle pi less theirs, or one of
      // each, the bisector along their difference.
      SmBoolean bFits = FALSE;
      for (ULONG c=0; c<4 && !bFits; c++)
        {
          SmVector3d sBisector = ( c < 2 ) ? asNormal[0] + asNormal[1] : asNormal[0] - asNormal[1];
          double     dAngle    = ( c < 2 ) ? SM_PI - dPhi : dPhi;
          double     dLen      = sBisector.Length();
          if ( dLen <= SM_EFF_ZERO_SQRT || dAngle <= SM_EFF_ZERO_SQRT )
            { continue; }
          SmVector3d sOut      = sBisector * ( (( c % 2 ) ? -1.0 : 1.0 ) / dLen ); // edge to centre
          double     dOff      = 1.0 / sin( dAngle / 2.0 );
          SmPoint3d  sCentre   = sPnt + sOut * ( dRadius * dOff );
          SmVector3d sAlong    = asEdge[1] + sOut * ( dSlope * dOff );           // dC/dT
          double     dSpeed    = sAlong.Length();
          if ( dSpeed <= SM_EFF_ZERO )
            { continue; }
          sAlong = sAlong / dSpeed;
          double     dRate     = dSlope / dSpeed;                                 // dr/ds
          if ( smos_Fabs( dRate ) >= 1.0 )
            { continue; }
          SmVector3d sToEdge   = -sOut - sAlong * ( -sOut ).Dot( sAlong );
          double     dToEdge   = sToEdge.Length();
          if ( dToEdge <= SM_EFF_ZERO )
            { continue; }
          SmPoint3d sMiddle = sCentre + sAlong * ( -dRadius * dRate ) + sToEdge * ( dRadius * sqrt( 1.0 - dRate * dRate ) / dToEdge );
          SmBoolean bDropped = FALSE, bMulti = FALSE;
          SmPoint2d sUV;
          double    dGap = SM_BIG_DOUBLE;
          SER( pFillet->DropPoint( sMiddle, pFillet->GetNaturalUVDomain(), NULL, bDropped, sUV, dGap, bMulti ));
          // the solver's tolerance, and its circular sections' relative accuracy (a non-rational
          // cubic across the round, within about 2e-5 of the radius)
          bFits = ( bDropped && dGap <= dTol + 1.0e-4 * dRadius ) ? TRUE : FALSE;
        }
      rbHolds = bFits;
    }
  return SM_SUCCESS;
} // end sm_HoldsItsLaw

/*******************************************************************//**
PURPOSE: SolidCraft change. Note an end of a fillet that its Nx1 corner left
    uncut in a cliff overflow (SmFilletNx1Corner::CalcCornerEdgeGeom).

NOTES: pEndFace's surface, continued, cuts the fillet off there
    (CutFilletAtOpenEnd).
***********************************************************************/
SmStatus SmFilletExecutive::NoteOpenEnd
  (SmFilletSolver * pSolver,   // in : the solver of the filleted edge
   SmVertex       * pVertex,   // in : the edge's vertex at that end
   SmFace         * pEndFace)  // in : the end face there
{
  NER(pSolver);
  NER(pVertex);
  NER(pEndFace);
  SmFilletOpenEnd sEnd;
  sEnd.m_pSolver  = pSolver;
  sEnd.m_pVertex  = pVertex;
  sEnd.m_pEndFace = pEndFace;
  m_vOpenEnds.push_back( sEnd );
  return SM_SUCCESS;
} // end SmFilletExecutive::NoteOpenEnd

/*******************************************************************//**
PURPOSE: SolidCraft change. Cut the fillet off at an end its Nx1 corner left
    uncut, by the end face's surface continued.

NOTES: The fillet face of the end's solver runs on past the end, its rails
    past their side edges. The end face, planar, continues as its plane: the
    plane sections the fillet face (SmSurface::CreatePlanarSectionCurves,
    exact for the fillet's analytic surfaces), its section nearest the
    edge's end (the end rib) splits it (SmBrep::MergeCurvesOnSurface), and of
    the two pieces the one beyond the
    end goes: the piece kept holds the fillet's point nearest the middle of
    the filleted edge. A rib that does not split the face in two is
    SM_FILERR_NO_ROOM, as is an end face that is not planar.
***********************************************************************/
SmStatus SmFilletExecutive::CutFilletAtOpenEnd
  (const SmFilletOpenEnd & crEnd)   // in : the end left uncut
{
  NER(crEnd.m_pSolver);
  NER(crEnd.m_pEndFace);
  SmFilletSolver * pSolver = crEnd.m_pSolver;
  double           dTol    = pSolver->GetThisApproxTol3d();

  const SmSurface * pEndSurface = crEnd.m_pEndFace->GetSurface();   NER(pEndSurface);
  SmPoint3d  sAt;
  SmVector3d sNormal;
  if ( !pEndSurface->IsPlanar( SM_EFF_ZERO, &sAt, &sNormal ))
    {
      GetFilletErrorInfo()->NoteError( SM_FILERR_NO_ROOM, _T("Fillet Error: a cliff overflow runs past a fillet end whose end face is not planar") );
      return SM_ERR;
    }
  SmPlane sCut( sAt, sNormal, &m_crContext );

  // the fillet's point nearest the middle of the filleted edge
  SmEdge * pFilEdge = pSolver->GetEdgeuse(0)->GetEdge();             NER(pFilEdge);
  SmPoint3d sEdgeMid;
  SER( pFilEdge->GetCurve()->EvaluatePoint( pFilEdge->GetInterval().Evaluate( 0.5 ), sEdgeMid ));

  SmTArray<SmFilletGeom*> sFGs;
  pSolver->GetFilletGeoms( sFGs );
  for (ULONG k=0; k<sFGs.GetSize(); k++)
    {
      SmFilletGeom * pFG = sFGs[k];
      if ( pFG == NULL || pFG->GetFilletSurface() == NULL )
        { continue; }
      SmFace * pFace = SM_CAST_PTR( SmFace, pFG->GetFilletSurface()->GetOwner() );
      if ( pFace == NULL || pFace->GetBrep() != m_pFilletBrep )
        { continue; }
      SmSurface * pSurface = pFace->GetSurface();                      NER(pSurface);

      SmBoolean bDropped = FALSE, bMulti = FALSE;
      SmPoint2d sMidUV;
      double    dGap = SM_BIG_DOUBLE;
      SER( pSurface->DropPoint( sEdgeMid, pFace->GetUVDomain(), NULL, bDropped, sMidUV, dGap, bMulti ));
      if ( !bDropped )
        { SER( SM_ERR ); }
      SmPoint3d sKeep;
      SER( pSurface->EvaluatePoint( sMidUV, sKeep ));

      // section the face by the end's plane, and split it there
      SmApproxTol3d       sApprox( dTol );
      SmTArray<SmCurve*>  s3DCurves, sUVCurves;
      SmObjsDelete<SmCurve*> sClean3D( &s3DCurves ), sCleanUV( &sUVCurves );
      SER( pSurface->CreatePlanarSectionCurves( m_crContext, pFace->GetUVDomain(), sAt, sNormal, &sApprox, NULL, &s3DCurves, &sUVCurves ));
      if ( s3DCurves.GetSize() == 0 )
        { continue; }  // this face does not reach the end
      // The end rib: the section nearest the edge's end (a fillet traced far past its end
      // can cross the plane again out there, where it is cut off anyway).
      ULONG  lRib = 0;
      double dNearest = SM_BIG_DOUBLE;
      for (ULONG c=0; c<s3DCurves.GetSize(); c++)
        {
          SmBoolean bOnCurve = FALSE;
          double    dParam = 0.0, dDist = SM_BIG_DOUBLE;
          SER( s3DCurves[c]->DropPoint( s3DCurves[c]->GetNaturalInterval(), crEnd.m_pVertex->GetPoint(), NULL, SM_BIG_DOUBLE, NULL, bOnCurve, dParam, dDist ));
          if ( bOnCurve && dDist < dNearest )
            {
              dNearest = dDist;
              lRib     = c;
            }
        }
      SmTArray<SmCurve*> sRib3D, sRibUV;
      sRib3D.Add( s3DCurves[lRib] );
      if ( sUVCurves.GetSize() == s3DCurves.GetSize() )
        { sRibUV.Add( sUVCurves[lRib] ); }
      SmTArray<SmFace*> sNewFaces;
      SmTArray<SmEdge*> sOnEdges;
      SER( m_pFilletBrep->MergeCurvesOnSurface( *pSurface, dTol, sRib3D,
                                                ( sRibUV.GetSize() == 1 ) ? &sRibUV : NULL,
                                                sNewFaces, sOnEdges, FALSE, &sCut ));
      if ( sNewFaces.GetSize() == 0 )
        { continue; }  // the section only touches it
      if ( sNewFaces.GetSize() != 1 )
        {
          GetFilletErrorInfo()->NoteError( SM_FILERR_NO_ROOM, _T("Fillet Error: the end face's plane does not cut a fillet that overflows its end in two") );
          return SM_ERR;
        }
      // keep the piece holding the fillet's middle, delete the one beyond the end
      SmFace * apPieces[2] = { pFace, sNewFaces[0] };
      SmBoolean abHolds[2] = { FALSE, FALSE };
      for (ULONG p=0; p<2; p++)
        {
          SmPointClassification sClass;
          SER( apPieces[p]->Point3DClassify( sKeep, SmTol::GetZoneTol3d( apPieces[p] ), TRUE, sClass ));
          abHolds[p] = ( sClass.GetPointClass() != SM_PC_UNKNOWN ) ? TRUE : FALSE;
        }
      if ( abHolds[0] == abHolds[1] )
        {
          GetFilletErrorInfo()->NoteError( SM_FILERR_NO_ROOM, _T("Fillet Error: the end face's plane does not part a fillet that overflows its end from its middle") );
          return SM_ERR;
        }
      SmFace * pBeyond = abHolds[0] ? apPieces[1] : apPieces[0];
      if ( abHolds[1] )
        {
          // the kept piece is the fillet surface's face from here on
          pFG->SetFilletSurface( apPieces[1]->GetSurface() );
        }
      SER( m_pFilletBrep->DeleteFace( pBeyond, TRUE ));
    }
  return SM_SUCCESS;
} // end SmFilletExecutive::CutFilletAtOpenEnd


/*******************************************************************//**
PURPOSE: SolidCraft change. How far along the blended edge pEdge, per unit of
    radius, the corner at its vertex pV takes it: a convex corner of three
    planar faces whose three edges are all blended at one radius ends its
    rounds on the plane through the corner's ball centre c = V + r u square to
    each edge, so the edge loses r <u, e> there (e its unit direction leaving
    V; u the point at unit distance inside each of the three planes, <n_f, u>
    = -1 for their outward unit normals). 0 at any other vertex: an end the
    round runs out to.
***********************************************************************/
static SmStatus sm_CornerTakesPerRadius
  (const SmVertex          & crV,       // in : the vertex
   const SmEdge            & crEdge,    // in : the blended edge
   const std::set<SmEdge*> & crBlended, // in : the edges blended (at one radius)
   double                  & rdPerR)    // out: the length taken per unit radius
{
  rdPerR = 0.0;
  SmTArray<SmEdge*> sEdges;
  crV.GetEdges(sEdges);
  if (sEdges.GetSize() != 3)
    { return SM_SUCCESS; }
  SmVector3d asDir[3];
  ULONG lThis = 3;
  for (ULONG i=0; i<3; i++)
    {
      if (!crBlended.count(sEdges[i]) || sEdges[i]->GetCurve() == NULL)
        { return SM_SUCCESS; }
      if (sEdges[i] == &crEdge)
        { lThis = i; }
      // the edge's direction leaving V
      const SmExtent1d sIvl = sEdges[i]->GetInterval();
      const SmBoolean bAtStart = (sEdges[i]->GetVertex() == &crV) ? TRUE : FALSE;
      SmVector3d sPV[2];
      SER(sEdges[i]->GetCurve()->Evaluate(bAtStart ? sIvl.GetMin() : sIvl.GetMax(), 1, bAtStart, sPV));
      asDir[i] = bAtStart ? sPV[1] : -sPV[1];
      SER(asDir[i].Unitize());
    }
  if (lThis == 3)
    { return SM_SUCCESS; }
  // each face between two of the edges: its plane's normal, outward (away from the third edge)
  SmVector3d asN[3];
  for (ULONG i=0; i<3; i++)
    {
      SmVector3d sN = asDir[(i+1)%3] * asDir[(i+2)%3];
      if (!(sN.Length() > SM_EFF_ZERO))
        { return SM_SUCCESS; }
      SER(sN.Unitize());
      if (sN.Dot(asDir[i]) > 0.0)
        { sN = -sN; }
      asN[i] = sN;
    }
  // the faces at V are planes: each face holding two of the edges
  SmTArray<SmFace*> sFaces;
  crV.GetFaces(sFaces);
  if (sFaces.GetSize() != 3)
    { return SM_SUCCESS; }
  for (ULONG f=0; f<3; f++)
    {
      if (SM_CAST_PTR(SmPlane, sFaces[f]->GetSurface()) == NULL)
        { return SM_SUCCESS; }
    }
  // <n_f, u> = -1 for the three planes (Cramer's rule)
  const SmVector3d & a = asN[0];
  const SmVector3d & b = asN[1];
  const SmVector3d & c = asN[2];
  const double dDet = a.Dot(b * c);
  if (!(smos_Fabs(dDet) > SM_EFF_ZERO))
    { return SM_SUCCESS; }
  const SmVector3d sU = ((b * c) + (c * a) + (a * b)) * (-1.0 / dDet);
  rdPerR = smos_Max(0.0, sU.Dot(asDir[lThis]));
  return SM_SUCCESS;
} // end sm_CornerTakesPerRadius

/*******************************************************************//**
PURPOSE: SolidCraft change. In a cliff overflow, trim the blends along
    their creases with the faces beyond the cliffs and sew them into the
    target: the blends built on their faces' continuations, the sharp edges
    moved to where they meet the body.

NOTES: See the file's head. Each sheet face is sectioned by every face of
    the target near it whose surface holds none of its boundary (the faces
    it rolls on and ends on are left out: it meets them along its rails and
    end ribs). A crease is kept where it lies inside the face it is on, and
    the sheet face is split along all of its creases at once. A piece
    whose inner point lies inside the target is kept; the others, beyond
    the cliffs, go. The kept pieces' lamina edges are then imprinted on
    the target faces they lie on (SmBrep::MergeCurvesOnSurface, which
    clips each to the face and cuts the edges it ends on); the target faces
    reached from the blended edges without crossing an imprinted edge are
    the material the blends take off, and are deleted; the sheet is merged
    in and stitched, and the body oriented. Refused (SM_FILERR_NO_ROOM): a
    variable round whose surface misses its law (sm_HoldsItsLaw). A sheet
    whose boundary does not close on the body (SM_FILERR_BAD_MERGE). A crease
    on a face that runs on tangentially from a face the sheet rolls on (the
    rail rolls over onto it: no cliff, SM_FILERR_NO_ROOM).
***********************************************************************/
SmStatus SmFilletExecutive::TrimAtCreases()
{
  NER(m_pFilletBrep);
  NER(m_pTargetBrep);

  for (ULONG i=0; i<m_vFilletSolvers.GetSize(); i++)
    {
      SmBoolean bHolds = TRUE;
      SER( sm_HoldsItsLaw( m_vFilletSolvers[i], bHolds ));
      if ( !bHolds )
        {
          GetFilletErrorInfo()->NoteError( SM_FILERR_NO_ROOM, _T("Fillet Error: a variable round traced on its faces' continuations does not hold its radius law") );
          return SM_ERR;
        }
    }

  // the ends the corners left uncut, cut off by their end faces' surfaces continued (an end's
  // event: the end face cuts the blend whatever lies past it, a lip over a wall's top)
  for (const SmFilletOpenEnd & rEnd : m_vOpenEnds)
    { SER( CutFilletAtOpenEnd( rEnd )); }

  double dTol = m_pFilletBrep->GetTolerance();
  for (ULONG i=0; i<m_vFilletSolvers.GetSize(); i++)
    { dTol = smos_Max( dTol, m_vFilletSolvers[i]->GetThisApproxTol3d() ); }
  const double dOnTol = 10.0 * dTol;
  double dTangency = 0.0;
  for (ULONG i=0; i<m_vFilletSolvers.GetSize(); i++)
    { dTangency = smos_Max( dTangency, m_vFilletSolvers[i]->GetTangencyTolerance() ); }

  // the blended edges' curves, to find their pieces after the imprint
  std::vector<SmCurve*> vBlended;
  SmTArray<SmCurve*> sBlendedOwned;
  SmObjsDelete<SmCurve*> sCleanOwned( &sBlendedOwned );
  for (ULONG i=0; i<m_vFilletSolvers.GetSize(); i++)
    {
      SmEdge * pEdge = m_vFilletSolvers[i]->GetEdgeuse(0)->GetEdge();   NER(pEdge);
      SmCurve * pCopy = NULL;
      SER( pEdge->GetCurve()->Copy( m_crContext, pCopy ));                NER(pCopy);
      sBlendedOwned.Add( pCopy );
      SmExtent1d sIvl = pEdge->GetInterval();
      if ( !( sIvl == pCopy->GetNaturalInterval() ))
        { SER( pCopy->Trim( sIvl )); }
      vBlended.push_back( pCopy );
    }

  SmTArray<SmFace*> sTargetFaces;
  m_pTargetBrep->GetFaces( sTargetFaces );
  // the faces the blends roll on: those of the blended edges (not the end faces they end on)
  std::set<SmFace*> sRolledOn;
  for (ULONG i=0; i<m_vFilletSolvers.GetSize(); i++)
    {
      SmTArray<SmFace*> sAtEdge;
      m_vFilletSolvers[i]->GetEdgeuse(0)->GetEdge()->GetFaces( sAtEdge );
      for (ULONG k=0; k<sAtEdge.GetSize(); k++)
        { sRolledOn.insert( sAtEdge[k] ); }
    }
  std::vector<SmExtent3d> vTargetBoxes( sTargetFaces.GetSize() );
  for (ULONG f=0; f<sTargetFaces.GetSize(); f++)
    { SER( sTargetFaces[f]->CalculateBoundingBox( vTargetBoxes[f] )); }

  // each blended edge's own fillet faces (not its corners')
  std::map<SmFace*, ULONG> sRoundOf;
  std::set<ULONG>          sKeptRounds;
  for (ULONG i=0; i<m_vFilletSolvers.GetSize(); i++)
    {
      SmTArray<SmFilletGeom*> sFGs;
      m_vFilletSolvers[i]->GetFilletGeoms( sFGs );
      for (ULONG k=0; k<sFGs.GetSize(); k++)
        {
          if ( sFGs[k] == NULL || sFGs[k]->GetFilletSurface() == NULL )
            { continue; }
          SmFace * pFace = SM_CAST_PTR( SmFace, sFGs[k]->GetFilletSurface()->GetOwner() );
          if ( pFace != NULL && pFace->GetBrep() == m_pFilletBrep )
            { sRoundOf[pFace] = i; }
        }
    }

  // 1, 2: each sheet face split along its creases, the pieces outside the body (beyond the
  // cliffs, or past an end) deleted, whether a crease split it or not
  SmTArray<SmFace*> sSheetFaces;
  m_pFilletBrep->GetFaces( sSheetFaces );
  for (ULONG s=0; s<sSheetFaces.GetSize(); s++)
    {
      SmFace * pSheet = sSheetFaces[s];
      SmExtent3d sBox;
      SER( pSheet->CalculateBoundingBox( sBox ));
      sBox.ExpandAbsolute( dOnTol );
      SmTArray<SmEdge*> sEdges;
      pSheet->GetEdges( sEdges );
      SmTArray<SmCurve*> sPieces;
      SmObjsDelete<SmCurve*> sCleanPieces( &sPieces );
      const SmSurface * pOnly = NULL;
      ULONG lFaces = 0;
      // the faces the sheet face rolls on or ends on: each holds one of its boundary edges
      std::set<SmFace*> sSupports;
      for (ULONG f=0; f<sTargetFaces.GetSize(); f++)
        {
          if ( vTargetBoxes[f].AreDisjoint( sBox, dOnTol ))
            { continue; }
          SmBoolean bHolds = FALSE;
          for (ULONG e=0; e<sEdges.GetSize() && !bHolds; e++)
            {
              if ( !sEdges[e]->IsLamina() )
                { continue; }
              SER( sm_EdgeOnSurface( *sEdges[e], *sTargetFaces[f], dOnTol, bHolds ));
            }
          if ( bHolds )
            { sSupports.insert( sTargetFaces[f] ); }
        }
      for (ULONG f=0; f<sTargetFaces.GetSize(); f++)
        {
          SmFace * pFace = sTargetFaces[f];
          if ( vTargetBoxes[f].AreDisjoint( sBox, dOnTol ) || sSupports.count( pFace ))
            { continue; }
          sm_Curves sCreases;
          // fitted to the body's own accuracy, so that creases meeting at a corner of the sheet
          // (on two faces past a cliff) meet there within it; exact sections stay exact
          SER( sm_Creases( m_crContext, *pSheet, *pFace, m_pTargetBrep->GetTolerance(), sCreases ));
          // each crease where it lies inside both faces
          SmTArray<SmCurve*> sInFace;
          SmObjsDelete<SmCurve*> sCleanInFace( &sInFace );
          for (ULONG c=0; c<sCreases.m_v3D.GetSize(); c++)
            { SER( sm_PiecesInFace( m_crContext, *sCreases.m_v3D[c], *pFace, dTol, sInFace )); }
          ULONG lBefore = sPieces.GetSize();
          for (ULONG c=0; c<sInFace.GetSize(); c++)
            { SER( sm_PiecesInFace( m_crContext, *sInFace[c], *pSheet, dTol, sPieces )); }
          if ( sPieces.GetSize() > lBefore )
            {
              // A face past a cliff meets the face the rail leaves (one the blend rolls on) at a
              // sharp edge. One that runs on from it tangentially (a round already on the face) is
              // no cliff: the rail would roll over onto it, and the ball does not fit between.
              SmTArray<SmEdge*> sFaceEdges;
              pFace->GetEdges( sFaceEdges );
              for (ULONG e=0; e<sFaceEdges.GetSize(); e++)
                {
                  SmTArray<SmFace*> sAt;
                  sFaceEdges[e]->GetFaces( sAt );
                  for (ULONG k=0; k<sAt.GetSize(); k++)
                    {
                      if ( sAt[k] != pFace && sSupports.count( sAt[k] ) && sRolledOn.count( sAt[k] )
                           && SmFilletSolver::IsTangentEdge( sFaceEdges[e], dTangency ))
                        {
                          GetFilletErrorInfo()->NoteError( SM_FILERR_NO_ROOM, _T("Fillet Error: a rail leaves the face it rolls on across a tangent edge, onto a face it would cut (the blend is larger than the face allows)") );
                          return SM_ERR;
                        }
                    }
                }
              lFaces++;
              pOnly = pFace->GetSurface();
            }
        }
      SmTArray<SmFace*> sNewFaces;
      SmTArray<SmEdge*> sOnEdges;
      if ( sPieces.GetSize() > 0 )
        {
          SER( m_pFilletBrep->MergeCurvesOnSurface( *pSheet->GetSurface(), dTol, sPieces, NULL, sNewFaces, sOnEdges, FALSE,
                                                    ( lFaces == 1 ) ? pOnly : NULL ));
        }
      // keep the pieces inside the target
      SmTArray<SmFace*> sAll;
      sAll.Add( pSheet );
      sAll.Append( sNewFaces );
      for (ULONG p=0; p<sAll.GetSize(); p++)
        {
          // a piece keeps the domain of the face it was split from: its own, from its trims
          SmExtent2d sDomain;
          SER( sAll[p]->CalculateUVDomainFromUVTrimCurves( sDomain, TRUE ));
          sAll[p]->SetUVDomain( sDomain );
          SmPoint3d sInner;
          SER( sAll[p]->CalculateAnInternalPoint( sInner ));
          SmBoolean bInside = FALSE, bOnBoundary = FALSE;
          SER( sm_InMaterial( *m_pTargetBrep, sInner, bInside, bOnBoundary ));
          if ( bOnBoundary )
            {
              GetFilletErrorInfo()->NoteError( SM_FILERR_BAD_MERGE, _T("Fillet Error: a blend past its faces lies on the body's boundary away from its creases") );
              return SM_ERR;
            }
          if ( !bInside )
            {
              SER( m_pFilletBrep->DeleteFace( sAll[p], TRUE ));
              sAll[p] = NULL;
            }
          else if ( sRoundOf.count( pSheet ))
            { sKeptRounds.insert( sRoundOf[pSheet] ); }
        }
    }

  // A blended edge whose blend keeps no face in the body has no length left: its corners take the
  // whole of it (a 20 cube's corner rounded R20: the corner's ball reaches the far faces, the
  // rounds lie past them). Larger than the faces allow, so a fit can reduce it.
  for (ULONG i=0; i<m_vFilletSolvers.GetSize(); i++)
    {
      SmBoolean bHas = FALSE;
      for (const std::pair<SmFace* const, ULONG> & rOf : sRoundOf)
        { bHas = ( bHas || rOf.second == i ) ? TRUE : FALSE; }
      if ( bHas && !sKeptRounds.count( i ))
        {
          GetFilletErrorInfo()->NoteError( SM_FILERR_NO_ROOM, _T("Fillet Error: a blend past its faces keeps no length of its edge: its corners take the whole of it (the blend is larger than the faces allow)") );
          // the largest radius that keeps a length of the edge: its corners take r <u, e> at
          // each end (sm_CornerTakesPerRadius), the edge's length L in all at r = L / (s0 + s1);
          // open (at it the round keeps no length), so less the fit margin (FitMargin)
          SmFilletSolver * pFS = m_vFilletSolvers[i];
          SmBoolean bOnEdge = FALSE;
          SmEdge * pEdge = ( pFS->GetEdgeuse(0) != NULL ) ? pFS->GetEdgeuse(0)->GetEdge() : NULL;
          if ( pEdge != NULL && pFS->RollsABall( bOnEdge ) && !bOnEdge )
            {
              std::set<SmEdge*> sBlended;
              for (ULONG k=0; k<m_vFilletSolvers.GetSize(); k++)
                {
                  SmBoolean bLaw = FALSE;
                  if ( m_vFilletSolvers[k]->RollsABall( bLaw ) && !bLaw
                       && smos_Fabs( m_vFilletSolvers[k]->GetSurface(0)->GetOffsetDistance() - pFS->GetSurface(0)->GetOffsetDistance() )
                            <= SM_EFF_ZERO * ( 1.0 + smos_Fabs( pFS->GetSurface(0)->GetOffsetDistance() )))
                    { sBlended.insert( m_vFilletSolvers[k]->GetEdgeuse(0)->GetEdge() ); }
                }
              double adPer[2] = { 0.0, 0.0 };
              SER( sm_CornerTakesPerRadius( *pEdge->GetVertex(), *pEdge, sBlended, adPer[0] ));
              SER( sm_CornerTakesPerRadius( *pEdge->GetOtherVertex( pEdge->GetVertex() ), *pEdge, sBlended, adPer[1] ));
              double dLength = 0.0;
              SER( pEdge->GetCurve()->Length( pEdge->GetInterval(), SM_EFF_ZERO, dLength ));
              if ( adPer[0] + adPer[1] > 0.0 )
                { NoteFitRadius( smos_Max( 0.0, dLength / ( adPer[0] + adPer[1] ) - FitMargin( *pEdge )), pEdge ); }
            }
          return SM_ERR;
        }
    }

  // the support faces continued to the kept rails that run past them into the material
  SmTArray<SmEdge*> sSheetEdges;
  m_pFilletBrep->GetEdges( sSheetEdges );

  std::set<SmEdge*> sBlendedEdges;
  for (ULONG i=0; i<m_vFilletSolvers.GetSize(); i++)
    { sBlendedEdges.insert( m_vFilletSolvers[i]->GetEdgeuse(0)->GetEdge() ); }
  std::vector<sm_Extension> vExtensions;
  SmTArray<SmCurve*> sAlong;
  SmObjsDelete<SmCurve*> sCleanAlong( &sAlong );
  for (ULONG f=0; f<sTargetFaces.GetSize(); f++)
    {
      std::vector<SmEdge*> vRails;
      for (ULONG e=0; e<sSheetEdges.GetSize(); e++)
        {
          if ( !sSheetEdges[e]->IsLamina() )
            { continue; }
          SmBoolean bOn = FALSE;
          SER( sm_EdgeOnSurface( *sSheetEdges[e], *sTargetFaces[f], dOnTol, bOn ));
          if ( bOn )
            { vRails.push_back( sSheetEdges[e] ); }
        }
      if ( vRails.empty() )
        { continue; }
      const size_t lBefore = vExtensions.size();
      const TCHAR * cWhy = NULL;
      SmStatus eStat = sm_Extend( m_crContext, *m_pTargetBrep, *sTargetFaces[f], vRails, sBlendedEdges, dTol, vExtensions, cWhy );
      if ( eStat != SM_SUCCESS && cWhy != NULL )
        {
          GetFilletErrorInfo()->NoteError( SM_FILERR_NO_ROOM, cWhy );
          return SM_ERR;
        }
      SER( eStat );
      for (size_t x=lBefore; x<vExtensions.size(); x++)
        {
          for (SmEdge * pEdge : vExtensions[x].m_vAlong)
            {
              SmCurve * pCurve = NULL;
              SER( sm_EdgeCurve( m_crContext, *pEdge, pCurve ));
              sAlong.Add( pCurve );
            }
        }
    }

  // 4: the kept sheet's boundary imprinted on the target faces it lies on
  std::set<SmEdge*> sImprinted;
  for (ULONG f=0; f<sTargetFaces.GetSize(); f++)
    {
      SmFace * pFace = sTargetFaces[f];
      SmTArray<SmCurve*> sOn;
      SmObjsDelete<SmCurve*> sCleanOn( &sOn );
      for (ULONG e=0; e<sSheetEdges.GetSize(); e++)
        {
          SmEdge * pEdge = sSheetEdges[e];
          if ( !pEdge->IsLamina() )
            { continue; }
          SmBoolean bOn = FALSE;
          SER( sm_EdgeOnSurface( *pEdge, *pFace, dOnTol, bOn ));
          if ( !bOn )
            { continue; }
          SmCurve * pCopy = NULL;
          SER( pEdge->GetCurve()->Copy( m_crContext, pCopy ));
          NER(pCopy);
          sOn.Add( pCopy );
          SmExtent1d sIvl = pEdge->GetInterval();
          if ( !( sIvl == pCopy->GetNaturalInterval() ))
            { SER( pCopy->Trim( sIvl )); }
        }
      if ( sOn.GetSize() == 0 )
        { continue; }
      SmTArray<SmFace*> sNewFaces;
      SmTArray<SmEdge*> sOnEdges;
      SER( m_pTargetBrep->MergeCurvesOnSurface( *pFace->GetSurface(), dTol, sOn, NULL, sNewFaces, sOnEdges, FALSE, NULL ));
      for (ULONG e=0; e<sOnEdges.GetSize(); e++)
        { sImprinted.insert( sOnEdges[e] ); }
    }

  // the faces round the blended edges, up to the imprint and the extensions' edges along the body:
  // the material the blends take off
  SmTArray<SmEdge*> sTargetEdges;
  m_pTargetBrep->GetEdges( sTargetEdges );
  for (ULONG e=0; e<sTargetEdges.GetSize(); e++)
    {
      SmEdge * pEdge = sTargetEdges[e];
      // a piece of an edge an extension runs along, cut at its crossing by the imprint
      for (const sm_Extension & rExt : vExtensions)
        {
          for (const std::pair<SmPoint3d, SmPoint3d> & rEnds : rExt.m_vAlongPieces)
            {
              const SmPoint3d sS = pEdge->GetStartVertex()->GetPoint(), sE = pEdge->GetEndVertex()->GetPoint();
              if (( sS.DistanceBetween( rEnds.first ) <= 10.0 * dTol && sE.DistanceBetween( rEnds.second ) <= 10.0 * dTol ) ||
                  ( sS.DistanceBetween( rEnds.second ) <= 10.0 * dTol && sE.DistanceBetween( rEnds.first ) <= 10.0 * dTol ))
                { sImprinted.insert( pEdge ); }
            }
        }
      if ( sAlong.GetSize() == 0 )
        { continue; }
      SmPoint3d sMid;
      SER( pEdge->GetCurve()->EvaluatePoint( pEdge->GetInterval().Evaluate( 0.5 ), sMid ));
      for (ULONG c=0; c<sAlong.GetSize(); c++)
        {
          SmBoolean bOnCurve = FALSE;
          double    dParam = 0.0, dDist = SM_BIG_DOUBLE;
          SER( sAlong[c]->DropPoint( sAlong[c]->GetNaturalInterval(), sMid, NULL, SM_BIG_DOUBLE, NULL, bOnCurve, dParam, dDist ));
          if ( bOnCurve && dDist <= dTol )
            { sImprinted.insert( pEdge ); break; }
        }
    }
  std::set<SmFace*> sGone;
  std::vector<SmFace*> vFront;
  for (ULONG e=0; e<sTargetEdges.GetSize(); e++)
    {
      SmEdge * pEdge = sTargetEdges[e];
      SmPoint3d sMid;
      SER( pEdge->GetCurve()->EvaluatePoint( pEdge->GetInterval().Evaluate( 0.5 ), sMid ));
      SmBoolean bBlended = FALSE;
      for (SmCurve * pCurve : vBlended)
        {
          SmBoolean bOnCurve = FALSE;
          double    dParam = 0.0, dDist = SM_BIG_DOUBLE;
          SER( pCurve->DropPoint( pCurve->GetNaturalInterval(), sMid, NULL, SM_BIG_DOUBLE, NULL, bOnCurve, dParam, dDist ));
          if ( bOnCurve && dDist <= dTol )
            { bBlended = TRUE; break; }
        }
      if ( !bBlended || sImprinted.count( pEdge ))
        { continue; }
      SmTArray<SmFace*> sFaces;
      pEdge->GetFaces( sFaces );
      for (ULONG k=0; k<sFaces.GetSize(); k++)
        {
          if ( sGone.insert( sFaces[k] ).second )
            { vFront.push_back( sFaces[k] ); }
        }
    }
  while ( !vFront.empty() )
    {
      SmFace * pFace = vFront.back();
      vFront.pop_back();
      SmTArray<SmEdge*> sEdges;
      pFace->GetEdges( sEdges );
      for (ULONG e=0; e<sEdges.GetSize(); e++)
        {
          if ( sImprinted.count( sEdges[e] ))
            { continue; }
          SmTArray<SmFace*> sFaces;
          sEdges[e]->GetFaces( sFaces );
          for (ULONG k=0; k<sFaces.GetSize(); k++)
            {
              if ( sGone.insert( sFaces[k] ).second )
                { vFront.push_back( sFaces[k] ); }
            }
        }
    }
  if ( sGone.empty() )
    {
      GetFilletErrorInfo()->NoteError( SM_FILERR_BAD_MERGE, _T("Fillet Error: a blend past its faces takes off nothing round its blended edges") );
      return SM_ERR;
    }
  if ( sGone.size() >= m_pTargetBrep->GetNumFaces() )
    {
      GetFilletErrorInfo()->NoteError( SM_FILERR_BAD_MERGE, _T("Fillet Error: a blend past its faces reaches every face of the body (its creases and rails enclose no part of it)") );
      return SM_ERR;
    }
  for (SmFace * pFace : sGone)
    { SER( m_pTargetBrep->DeleteFace( pFace, TRUE )); }

  // the fillet edges' attribute stays with the executive's own Breps
  for (ULONG i=0; i<sSheetEdges.GetSize(); i++)
    {
      if ( sSheetEdges[i]->FindAttribute( SM_AI_ALL_FILLET_EDGES ) == m_pAttributeE )
        { sSheetEdges[i]->RemoveAttribute( m_pAttributeE ); }
    }

  // The sheet's boundary is cut where an extension's loop meets it: a rib that passes under a lip
  // runs on the end face, under the lip (the extension's edge) and on the end face again, three
  // edges of the body for one of the sheet's.
  {
    for (const sm_Extension & rExt : vExtensions)
      {
        for (const sm_LoopCurve & rC : rExt.m_vLoop)
          {
            SmPoint3d sA, sB;
            SER( sm_CurveEnds( *rC.m_pCurve, sA, sB ));
            for (const SmPoint3d & rP : { sA, sB })
              {
                SmEdge * pFound = NULL;
                double   dAt = 0.0;
                SmTArray<SmEdge*> sEdges;
                m_pFilletBrep->GetEdges( sEdges );
                for (ULONG e=0; e<sEdges.GetSize() && pFound == NULL; e++)
                  {
                    SmEdge * pEdge = sEdges[e];
                    if ( !pEdge->IsLamina() || pEdge->GetStartVertex()->GetPoint().DistanceBetween( rP ) <= 10.0 * dTol
                         || pEdge->GetEndVertex()->GetPoint().DistanceBetween( rP ) <= 10.0 * dTol )
                      { continue; }
                    SmBoolean bOnCurve = FALSE;
                    double    dParam = 0.0, dDist = SM_BIG_DOUBLE;
                    SER( pEdge->GetCurve()->DropPoint( pEdge->GetInterval(), rP, NULL, SM_BIG_DOUBLE, NULL, bOnCurve, dParam, dDist ));
                    if ( bOnCurve && dDist <= 10.0 * dTol )
                      { pFound = pEdge; dAt = dParam; }
                  }
                if ( pFound != NULL )
                  {
                    SmEdge   * pOne = NULL, * pTwo = NULL;
                    SmVertex * pNew = NULL;
                    SER( m_pFilletBrep->MakeVertexSplitEdge( pFound, dAt, pOne, pTwo, pNew ));
                  }
              }
          }
      }
  }

  // the sheet sewn in, and the body oriented. The sheet is merged in as the fillet's normal path
  // merges it (SmBrep::MergeBrep, surfaces kept as they are): its edges keep their own
  // tolerances and the body its own (SmBrep::CopyFaces sets the body's tolerance to the sheet's,
  // which raises every edge of the body to it).
  SmTArray<SmFace*> sBefore, sAfter, sCopies;
  m_pTargetBrep->GetFaces( sBefore );
  std::set<SmFace*> sOld;
  for (ULONG f=0; f<sBefore.GetSize(); f++)
    { sOld.insert( sBefore[f] ); }
  SER( m_pTargetBrep->MergeBrep( *m_pFilletBrep, NULL, NULL, FALSE ));
  m_pTargetBrep->GetFaces( sAfter );
  for (ULONG f=0; f<sAfter.GetSize(); f++)
    {
      if ( !sOld.count( sAfter[f] ))
        { sCopies.Add( sAfter[f] ); }
    }
  SmBrep * pSheet = m_pFilletBrep;
  m_pFilletBrep = NULL;
  delete pSheet;
  for (sm_Extension & rExt : vExtensions)
    {
      SmFace * pFace = NULL;
      SER( sm_MakeExtensionFace( m_crContext, *m_pTargetBrep, rExt, dTol, pFace ));
      sCopies.Add( pFace );
    }
  ULONG  lStitched = 0, lLamina = 0;
  double dVertexGap = 0.0, dEdgeGap = 0.0;
  SER( m_pTargetBrep->StitchFaces( dTol, lStitched, lLamina, dVertexGap, dEdgeGap, &sCopies ));
  // the edges the sheet is sewn along hold it to its own accuracy
  for (ULONG f=0; f<sCopies.GetSize(); f++)
    {
      SmTArray<SmEdge*> sEdges;
      sCopies[f]->GetEdges( sEdges );
      for (ULONG e=0; e<sEdges.GetSize(); e++)
        { sEdges[e]->SetTolerance( smos_Max( dTol, smos_Max( dVertexGap, dEdgeGap ))); }
    }
  SmBoolean bMaybeOpen = FALSE;
  if ( lLamina == 0 )
    { SER( m_pTargetBrep->OrientTrimmedSurfaces( TRUE, bMaybeOpen, FALSE )); }
  if ( lLamina != 0 || bMaybeOpen || !m_pTargetBrep->IsManifoldSolid() )
    {
      GetFilletErrorInfo()->NoteError( SM_FILERR_BAD_MERGE, _T("Fillet Error: a blend past its faces does not sew into the body along its creases and rails") );
      return SM_ERR;
    }

  // Faces sewn on one surface are one face, as a boolean's post-process leaves them (only the
  // edges it made are looked at, SmMerge::BooleanPostProcess): a corner's end face on the end
  // wall's plane carried on over a rim is that wall, past the edges the blend sewed in. Every
  // other edge and vertex of the body stays.
  {
    std::set<SmEdge*>   sSewnEdges;
    std::set<SmVertex*> sSewnVertices;
    for (ULONG f=0; f<sCopies.GetSize(); f++)
      {
        SmTArray<SmEdge*> sEdges;
        sCopies[f]->GetEdges( sEdges );
        for (ULONG e=0; e<sEdges.GetSize(); e++)
          {
            sSewnEdges.insert( sEdges[e] );
            sSewnVertices.insert( sEdges[e]->GetStartVertex() );
            sSewnVertices.insert( sEdges[e]->GetEndVertex() );
          }
      }
    SmTArray<SmEdge*>   sAllEdges, sKeepEdges;
    SmTArray<SmVertex*> sAllVertices, sKeepVertices;
    m_pTargetBrep->GetEdges( sAllEdges );
    m_pTargetBrep->GetVertices( sAllVertices );
    for (ULONG e=0; e<sAllEdges.GetSize(); e++)
      {
        if ( !sSewnEdges.count( sAllEdges[e] ))
          { sKeepEdges.Add( sAllEdges[e] ); }
      }
    for (ULONG v=0; v<sAllVertices.GetSize(); v++)
      {
        if ( !sSewnVertices.count( sAllVertices[v] ))
          { sKeepVertices.Add( sAllVertices[v] ); }
      }
    SER( m_pTargetBrep->RemoveTopologicalEdgesAndVertices( &sKeepEdges, &sKeepVertices ));
  }
  return SM_SUCCESS;
} // end SmFilletExecutive::TrimAtCreases
