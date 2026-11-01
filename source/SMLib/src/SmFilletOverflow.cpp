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
*   1. finds each sheet face's creases: its surface sectioned by the surface of
*      each face of the body it may cross (exactly by a plane: a line on a
*      cylinder, a conic on a cone, a circle on a torus; a certified surface
*      intersection otherwise), each kept where it lies on that face;
*   2. splits the sheet face along them and keeps the pieces inside the body;
*   3. continues a support face on its plane to a kept rail that runs past the
*      face's foot into the material (a chamfer's leg past its wall's foot: the
*      wall carries on down to the leg), between the lines of its side edges;
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
  sm_Extension() = default;
  sm_Extension(const sm_Extension &) = delete;
  sm_Extension & operator=(const sm_Extension &) = delete;
  sm_Extension(sm_Extension && rOther) noexcept
    : m_sAt( rOther.m_sAt ), m_sNormal( rOther.m_sNormal ), m_vLoop( std::move( rOther.m_vLoop )), m_vAlong( std::move( rOther.m_vAlong ))
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

/* The extension of the planar support face crFace to the kept rails on its plane that run past
   its boundary into the material (a chamfer whose leg passes the foot of its wall: the wall's
   plane carries on down to the leg). Each chain of such rails ends at vertices on the lines of
   the face's side edges (the walls beside it meet it there); the loop runs along the chain, along
   each end's side line to the side edge's end, and along the face's boundary between those ends,
   away from the blended edges. rcWhy is set where the extension is not one this builds. */
SmStatus sm_Extend(const SmContext & crContext, const SmBrep & crTarget, const SmFace & crFace, const std::vector<SmEdge*> & crRails,
                   const std::set<SmEdge*> & crBlended, double dTol, std::vector<sm_Extension> & rExtensions, const TCHAR *& rcWhy)
{
  rcWhy = NULL;
  // the rails wholly off the face, in the material
  std::vector<SmEdge*> vOff;
  for (SmEdge * pRail : crRails)
    {
      SmCurve * pCurve = NULL;
      SER( sm_EdgeCurve( crContext, *pRail, pCurve ));
      SmObjDelete sClean( pCurve );
      SmCurveClassification sClass( pCurve, pCurve->GetNaturalInterval(), NULL, dTol );
      SER( crFace.CurveOnClassify( TRUE, sClass, TRUE ));
      SmBoolean bOff = FALSE, bOn = FALSE;
      for (ULONG j=0; j<sClass.GetSize(); j++)
        {
          SmExtent1d sRun = sClass.GetSubInterval( j );
          SmPoint3d sA, sB;
          SER( pCurve->EvaluatePoint( sRun.GetMin(), sA ));
          SER( pCurve->EvaluatePoint( sRun.GetMax(), sB ));
          if ( sA.DistanceBetween( sB ) <= dTol )
            { continue; }
          if ( sClass[j].m_vMid.GetPointClass() == SM_PC_UNKNOWN )
            {
              SmPoint3d sMid;
              SER( pCurve->EvaluatePoint( sRun.Evaluate( 0.5 ), sMid ));
              SmBoolean bIn = FALSE, bOnBoundary = FALSE;
              SER( sm_InMaterial( crTarget, sMid, bIn, bOnBoundary ));
              if ( bIn )
                { bOff = TRUE; }
              else
                { bOn = TRUE; }  // in the open: cut off along a crease elsewhere
            }
          else
            { bOn = TRUE; }
        }
      if ( bOff && bOn )
        {
          rcWhy = _T("Fillet Error: a blend's rail runs past a face's boundary into the material partway along it");
          return SM_ERR;
        }
      if ( bOff )
        { vOff.push_back( pRail ); }
    }
  if ( vOff.empty() )
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

  // the chains of off rails, by the sheet's vertices
  std::set<SmEdge*> sDone;
  for (SmEdge * pFirst : vOff)
    {
      if ( sDone.count( pFirst ))
        { continue; }
      std::vector<SmEdge*>   vChain( 1, pFirst );
      std::vector<SmVertex*> vVerts = { pFirst->GetStartVertex(), pFirst->GetEndVertex() };
      sDone.insert( pFirst );
      for (ULONG lSide=0; lSide<2; lSide++)
        {
          for (;;)
            {
              SmVertex * pAt = ( lSide == 0 ) ? vVerts.back() : vVerts.front();
              SmEdge   * pNext = NULL;
              for (SmEdge * pE : vOff)
                {
                  if ( !sDone.count( pE ) && ( pE->GetStartVertex() == pAt || pE->GetEndVertex() == pAt ))
                    { pNext = pE; break; }
                }
              if ( pNext == NULL )
                { break; }
              sDone.insert( pNext );
              if ( lSide == 0 )
                { vChain.push_back( pNext ); vVerts.push_back( pNext->GetOtherVertex( pAt )); }
              else
                { vChain.insert( vChain.begin(), pNext ); vVerts.insert( vVerts.begin(), pNext->GetOtherVertex( pAt )); }
            }
        }
      if ( vVerts.front() == vVerts.back() )
        {
          rcWhy = _T("Fillet Error: a blend's rail runs all round past a face's boundary into the material");
          return SM_ERR;
        }
      sm_Extension sExt;
      sExt.m_sAt     = sAt;
      sExt.m_sNormal = sNormal;
      for (size_t e=0; e<vChain.size(); e++)
        {
          SmCurve * pCurve = NULL;
          SER( sm_EdgeCurve( crContext, *vChain[e], pCurve ));
          SER( sm_AddCurve( sExt.m_vLoop, pCurve, vVerts[e]->GetPoint(), vVerts[e+1]->GetPoint() ));
        }
      // each end along its side line to the side edge's end on the face
      SmEdge   * apSide[2] = { NULL, NULL };
      SmVertex * apFoot[2] = { NULL, NULL };
      const SmPoint3d asEnd[2] = { vVerts.front()->GetPoint(), vVerts.back()->GetPoint() };
      for (ULONG k=0; k<2; k++)
        {
          SER( sm_SideLineThrough( crFace, asEnd[k], crBlended, 10.0 * dTol, apSide[k], apFoot[k] ));
          if ( apSide[k] == NULL )
            {
              rcWhy = _T("Fillet Error: a blend's rail past a face's boundary ends off the lines of the face's side edges");
              return SM_ERR;
            }
        }
      if ( asEnd[1].DistanceBetween( apFoot[1]->GetPoint() ) > dTol )
        {
          SmBSplineCurve * pLine = NULL;
          SER( SmBSplineCurve::CreateLineSegment( crContext, 3, asEnd[1], apFoot[1]->GetPoint(), pLine ));
          SER( sm_AddCurve( sExt.m_vLoop, pLine, asEnd[1], apFoot[1]->GetPoint() ));
        }
      // along the face's boundary from the end's foot to the start's, off the blended and side edges
      SmTArray<SmEdge*> sEdges;
      crFace.GetEdges( sEdges );
      SmVertex * pAt = apFoot[1];
      SmEdge   * pFrom = NULL;
      for (ULONG lSteps=0; pAt != apFoot[0]; lSteps++)
        {
          SmEdge * pNext = NULL;
          for (ULONG e=0; e<sEdges.GetSize() && lSteps<=sEdges.GetSize(); e++)
            {
              SmEdge * pE = sEdges[e];
              if ( pE == pFrom || pE == apSide[0] || pE == apSide[1] || crBlended.count( pE )
                   || ( pE->GetStartVertex() != pAt && pE->GetEndVertex() != pAt ))
                { continue; }
              if ( pNext != NULL )
                { pNext = NULL; break; }  // a branch: not one path
              pNext = pE;
            }
          if ( pNext == NULL )
            {
              rcWhy = _T("Fillet Error: a blend's rail past a face's boundary does not close on the face between its side edges");
              return SM_ERR;
            }
          SmVertex * pOther = pNext->GetOtherVertex( pAt );   NER(pOther);
          SmCurve  * pCurve = NULL;
          SER( sm_EdgeCurve( crContext, *pNext, pCurve ));
          SER( sm_AddCurve( sExt.m_vLoop, pCurve, pAt->GetPoint(), pOther->GetPoint() ));
          sExt.m_vAlong.push_back( pNext );
          pFrom = pNext;
          pAt   = pOther;
        }
      if ( asEnd[0].DistanceBetween( apFoot[0]->GetPoint() ) > dTol )
        {
          SmBSplineCurve * pLine = NULL;
          SER( SmBSplineCurve::CreateLineSegment( crContext, 3, apFoot[0]->GetPoint(), asEnd[0], pLine ));
          SER( sm_AddCurve( sExt.m_vLoop, pLine, apFoot[0]->GetPoint(), asEnd[0] ));
        }
      rExtensions.push_back( std::move( sExt ));
    }
  return SM_SUCCESS;
}

// The face on the plane through crAt square to crNormal bounded by rExt's loop (whose curves
// it takes), in crBrep's infinite region: the plane over the loop's box, the loop anticlockwise.
SmStatus sm_MakeExtensionFace(const SmContext & crContext, SmBrep & rBrep, sm_Extension & rExt, double dTol, SmFace *& rpFace)
{
  rpFace = NULL;
  SmVector3d sAxis( 1.0, 0.0, 0.0 );
  const SmVector3d & rN = rExt.m_sNormal;
  if ( smos_Fabs( rN.y ) < smos_Fabs( rN.x ) && smos_Fabs( rN.y ) <= smos_Fabs( rN.z ))
    { sAxis = SmVector3d( 0.0, 1.0, 0.0 ); }
  else if ( smos_Fabs( rN.z ) < smos_Fabs( rN.x ))
    { sAxis = SmVector3d( 0.0, 0.0, 1.0 ); }
  SmVector3d sX = rN * sAxis;
  if ( sX.Unitize() != SM_SUCCESS )
    { SER( SM_ERR ); }
  SmVector3d sY = rN * sX;
  // the loop's points in the plane's coordinates, in its order: its box and signed area
  std::vector<SmPoint2d> vUV;
  for (const sm_LoopCurve & rC : rExt.m_vLoop)
    {
      SmExtent1d sIvl = rC.m_pCurve->GetNaturalInterval();
      for (ULONG i=0; i<8; i++)
        {
          const double dT = ( rC.m_eOrient == SM_OT_SAME ) ? i / 8.0 : 1.0 - i / 8.0;
          SmPoint3d sPnt;
          SER( rC.m_pCurve->EvaluatePoint( sIvl.Evaluate( dT ), sPnt ));
          const SmVector3d sD = sPnt - rExt.m_sAt;
          vUV.push_back( SmPoint2d( sD.Dot( sX ), sD.Dot( sY )));
        }
    }
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
  std::vector<SmExtent3d> vTargetBoxes( sTargetFaces.GetSize() );
  for (ULONG f=0; f<sTargetFaces.GetSize(); f++)
    { SER( sTargetFaces[f]->CalculateBoundingBox( vTargetBoxes[f] )); }

  // 1, 2: each sheet face split along its creases, the pieces beyond the cliffs deleted
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
              // A face past a cliff meets the face the rail leaves at a sharp edge. One that
              // runs on from it tangentially (a round already on the face) is no cliff: the rail
              // would roll over onto it, and the ball does not fit between.
              SmTArray<SmEdge*> sFaceEdges;
              pFace->GetEdges( sFaceEdges );
              for (ULONG e=0; e<sFaceEdges.GetSize(); e++)
                {
                  SmTArray<SmFace*> sAt;
                  sFaceEdges[e]->GetFaces( sAt );
                  for (ULONG k=0; k<sAt.GetSize(); k++)
                    {
                      if ( sAt[k] != pFace && sSupports.count( sAt[k] )
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
      if ( sPieces.GetSize() == 0 )
        { continue; }
      SmTArray<SmFace*> sNewFaces;
      SmTArray<SmEdge*> sOnEdges;
      SER( m_pFilletBrep->MergeCurvesOnSurface( *pSheet->GetSurface(), dTol, sPieces, NULL, sNewFaces, sOnEdges, FALSE,
                                                ( lFaces == 1 ) ? pOnly : NULL ));
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
  for (ULONG e=0; e<sTargetEdges.GetSize() && sAlong.GetSize() > 0; e++)
    {
      SmEdge * pEdge = sTargetEdges[e];
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
  if ( sGone.empty() || sGone.size() >= m_pTargetBrep->GetNumFaces() )
    {
      GetFilletErrorInfo()->NoteError( SM_FILERR_BAD_MERGE, _T("Fillet Error: a blend past its faces does not close on the body (its creases and rails enclose no part of it)") );
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
  return SM_SUCCESS;
} // end SmFilletExecutive::TrimAtCreases
