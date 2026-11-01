// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

/**********************************************************************//**
* FILE NAME --- SmFilletKeepFeatures.cpp
* PURPOSE: SolidCraft change. Keep features: a feature standing inside a
*   round's band (the strip of face the round takes, from the edge to its
*   contact line) stays. A feature here is a boss or a cut whose footprint is
*   an inner loop of a face the round rolls on.
*
* The round is built along its whole edge as if the feature were not there
* (CreateFilletBrep: m_pFilletBrep is the sheet of fillet and corner faces).
* KeepFeaturesInBands then trims it in place:
*   1. the feature's walls are carried to the round: a cut's own wall faces
*      (they run on below the face into the body), a boss's side faces
*      continued down from its footprint into the body (its footprint swept
*      along the face's normal; the walls must be upright to a planar face);
*   2. each sheet face is sectioned by each wall (exactly by a planar wall: a
*      conic on the round; a certified surface intersection otherwise, a
*      quartic for a cylinder's wall on a cylindrical round), the crease kept
*      where it lies on both; a feature with a crease is one the band reaches,
*      noted in the trace (SM_ERR_MESSAGE);
*   3. the sheet face is split along the creases and the pieces bounded by one
*      feature's creases alone, over its footprint, go; a boss's continued
*      walls are split along them too and keep their pieces up to the
*      footprint;
*   4. the kept sheet's boundary (rails, end ribs, a cut's creases) is
*      imprinted on the body's faces (a cut's wall is cut back at the round),
*      the faces it encloses round the blended edges are deleted, stopping at
*      a kept boss's footprint, and the sheet and the boss's continued walls
*      are sewn in, each continued wall joined to the wall it continues.
* Nothing is subtracted. The approach is SmFilletOverflow.cpp's (the cliff
* overflow trimmed along its creases); its helpers are repeated here, kept
* apart, under kf_ names.
**********************************************************************/

#include "StdAfx.h"

#include <SmFilletExecutive.h>
#include <SmFilletGeom.h>
#include <SmFilletSolver.h>
#include <SmBrep.h>
#include <SmFace.h>
#include <SmEdge.h>
#include <SmLoop.h>
#include <SmVertex.h>
#include <SmRegion.h>
#include <SmSurface.h>
#include <SmPlane.h>
#include <SmCurve.h>
#include <SmCurveClass.h>
#include <SmPointClass.h>
#include <SmExtent3d.h>
#include <SmPrimitiveCreation.h>
#include <SmMessages.h>
#include <SmTol.h>

#include <algorithm>
#include <map>
#include <set>
#include <vector>

namespace
{

// The plane a planar surface lies on, by a point and its unit normal (rbPlanar FALSE when it is
// not planar).
SmStatus kf_PlaneOf(const SmSurface & crSurface, SmPoint3d & rAt, SmVector3d & rNormal, SmBoolean & rbPlanar)
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
SmStatus kf_DistanceToSurface(const SmFace & crFace, const SmPoint3d & crPnt, double & rdDist)
{
  rdDist = SM_BIG_DOUBLE;
  const SmSurface * pSurface = crFace.GetSurface();   NER(pSurface);
  SmPoint3d  sAt;
  SmVector3d sNormal;
  SmBoolean  bPlanar = FALSE;
  SER( kf_PlaneOf( *pSurface, sAt, sNormal, bPlanar ));
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
SmStatus kf_EdgeOnSurface(const SmEdge & crEdge, const SmFace & crFace, double dTol, SmBoolean & rbOn)
{
  rbOn = TRUE;
  const SmCurve * pCurve = crEdge.GetCurve();   NER(pCurve);
  const SmExtent1d sIvl = crEdge.GetInterval();
  for (ULONG i=0; i<=4 && rbOn; i++)
    {
      SmPoint3d sPnt;
      SER( pCurve->EvaluatePoint( sIvl.Evaluate( i / 4.0 ), sPnt ));
      double dDist = SM_BIG_DOUBLE;
      SER( kf_DistanceToSurface( crFace, sPnt, dDist ));
      rbOn = ( dDist <= dTol ) ? TRUE : FALSE;
    }
  return SM_SUCCESS;
}

// Whether the edge lies along one of crCurves: its ends and three points between within dTol.
SmStatus kf_EdgeAlong(const SmEdge & crEdge, const SmTArray<SmCurve*> & crCurves, double dTol, SmBoolean & rbAlong)
{
  rbAlong = FALSE;
  const SmCurve * pCurve = crEdge.GetCurve();   NER(pCurve);
  const SmExtent1d sIvl = crEdge.GetInterval();
  for (ULONG c=0; c<crCurves.GetSize() && !rbAlong; c++)
    {
      SmBoolean bOn = TRUE;
      for (ULONG i=0; i<=4 && bOn; i++)
        {
          SmPoint3d sPnt;
          SER( pCurve->EvaluatePoint( sIvl.Evaluate( i / 4.0 ), sPnt ));
          SmBoolean bDropped = FALSE;
          double    dParam = 0.0, dDist = SM_BIG_DOUBLE;
          SER( crCurves[c]->DropPoint( crCurves[c]->GetNaturalInterval(), sPnt, NULL, SM_BIG_DOUBLE, NULL, bDropped, dParam, dDist ));
          bOn = ( bDropped && dDist <= dTol ) ? TRUE : FALSE;
        }
      rbAlong = bOn;
    }
  return SM_SUCCESS;
}

// Owns the curves an intersection hands back.
struct kf_Curves
{
  SmTArray<SmCurve*> m_v3D, m_vUV1, m_vUV2;
  ~kf_Curves()
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

// The curves where crSheet's surface meets crWall's: a planar wall sections it exactly
// (SmSurface::CreatePlanarSectionCurves), any other wall meets it by SMLib's surface
// intersection over both faces' domains (within dTol), its touching points left out. Into
// rCurves (3D).
SmStatus kf_Creases(const SmContext & crContext, const SmFace & crSheet, const SmFace & crWall, double dTol, kf_Curves & rCurves)
{
  const SmSurface * pSheet = crSheet.GetSurface();   NER(pSheet);
  const SmSurface * pWall  = crWall.GetSurface();    NER(pWall);
  SmApproxTol3d sApprox( dTol );
  SmPoint3d  sAt;
  SmVector3d sNormal;
  SmBoolean  bPlanar = FALSE;
  SER( kf_PlaneOf( *pWall, sAt, sNormal, bPlanar ));
  if ( bPlanar )
    {
      SER( pSheet->CreatePlanarSectionCurves( crContext, crSheet.GetUVDomain(), sAt, sNormal, &sApprox, NULL, &rCurves.m_v3D, &rCurves.m_vUV1 ));
      return SM_SUCCESS;
    }
  SmBoolean abEdges[2] = { TRUE, TRUE };
  SmTArray<SmTsectCurveType> sTypes;
  kf_Curves sAll;
  SER( pSheet->GlobalSurfaceIntersect( crContext, crSheet.GetUVDomain(), *pWall, crWall.GetUVDomain(), abEdges, &sApprox, NULL,
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
SmStatus kf_PiecesInFace(const SmContext & crContext, const SmCurve & crCurve, const SmFace & crFace, double dTol,
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

// A copy of the edge's curve over the edge's interval.
SmStatus kf_EdgeCurve(const SmContext & crContext, const SmEdge & crEdge, SmCurve *& rpCurve)
{
  rpCurve = NULL;
  const SmCurve * pCurve = crEdge.GetCurve();   NER(pCurve);
  SER( pCurve->Copy( crContext, rpCurve ));     NER(rpCurve);
  SmExtent1d sIvl = crEdge.GetInterval();
  if ( !( sIvl == rpCurve->GetNaturalInterval() ))
    { SER( rpCurve->Trim( sIvl )); }
  return SM_SUCCESS;
}

// Whether every normal of crFace at a 5 x 5 grid over its domain is square to crUp: a wall
// made of lines along crUp (a prism's side, a cylinder about crUp).
SmStatus kf_Upright(const SmFace & crFace, const SmVector3d & crUp, double dAngTol, SmBoolean & rbUpright)
{
  rbUpright = TRUE;
  const SmSurface * pSurface = crFace.GetSurface();   NER(pSurface);
  SmExtent2d sBox  = crFace.GetUVDomain();
  SmPoint2d  sSize = sBox.GetSize();
  for (ULONG i=0; i<=4 && rbUpright; i++)
    {
      for (ULONG j=0; j<=4 && rbUpright; j++)
        {
          SmVector3d sNormal;
          SER( pSurface->EvaluateNormal( SmPoint2d( sBox.GetMin().x + sSize.x * i / 4.0, sBox.GetMin().y + sSize.y * j / 4.0 ),
                                         TRUE, TRUE, sNormal ));
          rbUpright = ( smos_Fabs( sNormal.Dot( crUp )) <= dAngTol ) ? TRUE : FALSE;
        }
    }
  return SM_SUCCESS;
}

// A feature on a face a round rolls on: its footprint (an inner loop of the face) and its walls.
struct kf_Feature
{
  SmFace             * m_pFace = NULL;      // the face it stands on (or is cut into)
  std::vector<SmEdge*> m_vRim;              // the footprint's edges
  SmTArray<SmCurve*>   m_vRimCurves;        // copies of them, owned
  SmBoolean            m_bBoss = FALSE;     // TRUE = a boss (its rim concave), FALSE = a cut
  SmVector3d           m_sUp;               // a boss: out of the face, along its walls
  std::vector<SmFace*> m_vWalls;            // the faces carried to the round: a cut's walls in the target,
                                            //   a boss's continued walls, built in the sheet (m_pFilletBrep)
  SmTArray<SmCurve*>   m_vCreases;          // where its walls meet the round, owned
  kf_Feature() = default;
  kf_Feature(const kf_Feature &) = delete;
  kf_Feature & operator=(const kf_Feature &) = delete;
  ~kf_Feature()
    {
      for (ULONG i=0; i<m_vRimCurves.GetSize(); i++) { delete m_vRimCurves[i]; }
      for (ULONG i=0; i<m_vCreases.GetSize(); i++)   { delete m_vCreases[i]; }
    }
};

// Whether crPnt lies in the material of crBrep (a solid region, not the infinite one or a void).
SmStatus kf_InMaterial(const SmBrep & crBrep, const SmPoint3d & crPnt, SmBoolean & rbIn)
{
  rbIn = FALSE;
  SmPointClassification sClass;
  SER( crBrep.Point3DClassify( crPnt, SmTol::GetZoneTol3d( &crBrep ), TRUE, sClass ));
  if ( sClass.GetPointClass() == SM_PC_REGION )
    {
      const SmRegion * pRegion = SM_CAST_PTR( SmRegion, sClass.GetObject() );
      rbIn = ( pRegion != NULL && !pRegion->IsInfiniteRegion() && !pRegion->IsVoid() ) ? TRUE : FALSE;
    }
  return SM_SUCCESS;
}

} // namespace

/*******************************************************************//**
PURPOSE: SolidCraft change. Keep features: trim the rounds in place round
    the features their bands reach, and sew them into the target.

NOTES: See the file's head. rbKept is FALSE, and nothing is changed, when
    no band reaches a feature: the fillet is merged as before. Refused: a
    boss whose walls are not upright to a planar face it stands on
    (SM_FILERR_NOT_EXACT: its walls continued are not built yet). A sheet
    whose boundary does not close on the body is SM_FILERR_BAD_MERGE.
***********************************************************************/
SmStatus SmFilletExecutive::KeepFeaturesInBands(SmBoolean & rbKept)
{
  rbKept = FALSE;
  NER(m_pFilletBrep);
  NER(m_pTargetBrep);

  double dTol = m_pFilletBrep->GetTolerance();
  double dAngTol = SM_EFF_ZERO_SQRT;
  for (ULONG i=0; i<m_vFilletSolvers.GetSize(); i++)
    {
      dTol    = smos_Max( dTol, m_vFilletSolvers[i]->GetThisApproxTol3d() );
      dAngTol = smos_Max( dAngTol, m_vFilletSolvers[i]->GetThisAngTolRad() );
    }
  const double dOnTol = 10.0 * dTol;
  // A crease off a planar wall is exact (a conic on the round); any other is fitted, certified
  // within a tenth of the fillet's tolerance, as the fillet's own fitted sections are
  // (SmGnSectionFSG): at the fillet's tolerance itself it would move the faces it bounds by as
  // much as the fillet may be off its own surface.
  const double dCreaseTol = 0.1 * dTol;

  // the faces the rounds roll on
  std::vector<SmFace*> vSupports;
  for (ULONG i=0; i<m_vFilletSolvers.GetSize(); i++)
    {
      SmFilletSolver * pSolver = m_vFilletSolvers[i];
      for (ULONG k=0; k<pSolver->m_vFilletGeoms.GetSize(); k++)
        {
          SmFilletGeom * pFG = pSolver->m_vFilletGeoms[k];
          if ( pFG == NULL || pFG->GetFilletGeomType() != SM_FG_DEFAULT )
            { continue; }
          for (ULONG ii=0; ii<2; ii++)
            {
              SmFilletEdge * pRail = pFG->GetRail( ii );
              SmFace * pFace = ( pRail != NULL ) ? pRail->GetOriginalFace() : NULL;
              if ( pFace != NULL && std::find( vSupports.begin(), vSupports.end(), pFace ) == vSupports.end() )
                { vSupports.push_back( pFace ); }
            }
        }
    }

  // the sheet's extent, and the depth below a face to which a boss's walls are continued:
  // past the whole sheet
  SmExtent3d sSheetBox;
  SER( m_pFilletBrep->CalculateBoundingBox( sSheetBox ));
  SmPoint3d asCorners[8];
  for (ULONG c=0; c<8; c++)
    {
      asCorners[c] = SmPoint3d( (c & 1) ? sSheetBox.GetMax().x : sSheetBox.GetMin().x,
                                (c & 2) ? sSheetBox.GetMax().y : sSheetBox.GetMin().y,
                                (c & 4) ? sSheetBox.GetMax().z : sSheetBox.GetMin().z );
    }
  const double dSheetSize = sSheetBox.GetMax().DistanceBetween( sSheetBox.GetMin() );

  // the sheet's own faces (a boss's continued walls are built in it below)
  SmTArray<SmFace*> sSheetFaces;
  m_pFilletBrep->GetFaces( sSheetFaces );

  // the blended edges: a loop holding one is the round's own edge (a branch's saddle on a pipe),
  // not a feature's footprint
  std::set<SmEdge*> sBlended;
  for (ULONG i=0; i<m_vFilletSolvers.GetSize(); i++)
    { sBlended.insert( m_vFilletSolvers[i]->GetEdgeuse(0)->GetEdge() ); }

  // 1: the features, each footprint an inner loop of a support face, and their walls
  std::vector<kf_Feature*> vFeatures;
  struct kf_Owner { std::vector<kf_Feature*> & m_r; ~kf_Owner() { for (kf_Feature * p : m_r) { delete p; } } } sOwner{ vFeatures };
  for (SmFace * pFace : vSupports)
    {
      SmTArray<SmLoop*> sLoops;
      pFace->GetLoops( sLoops );
      for (ULONG l=0; l<sLoops.GetSize(); l++)
        {
          if ( sLoops[l]->IsOuterLoop() )
            { continue; }
          SmTArray<SmEdge*> sRim;
          sLoops[l]->GetEdges( sRim );
          SmBoolean bBlendedRim = FALSE;
          for (ULONG e=0; e<sRim.GetSize() && !bBlendedRim; e++)
            { bBlendedRim = sBlended.count( sRim[e] ) ? TRUE : FALSE; }
          if ( sRim.GetSize() == 0 || bBlendedRim )
            { continue; }
          ULONG lConcave = 0;
          std::vector<SmFace*> vWalls;
          SmBoolean bManifold = TRUE;
          for (ULONG e=0; e<sRim.GetSize(); e++)
            {
              SmBoolean bConcave = FALSE;
              SER( SmFilletSolver::IsConcaveEdge( *sRim[e], bConcave ));
              if ( bConcave )
                { lConcave++; }
              SmTArray<SmFace*> sFaces;
              sRim[e]->GetFaces( sFaces );
              if ( sFaces.GetSize() != 2 )
                { bManifold = FALSE; break; }
              SmFace * pWall = ( sFaces[0] == pFace ) ? sFaces[1] : sFaces[0];
              if ( std::find( vWalls.begin(), vWalls.end(), pWall ) == vWalls.end() )
                { vWalls.push_back( pWall ); }
            }
          // a footprint whose rim is all of one kind: a boss's, or a cut's
          if ( !bManifold || ( lConcave != 0 && lConcave != sRim.GetSize() ))
            { continue; }
          kf_Feature * pFeature = new kf_Feature;
          vFeatures.push_back( pFeature );
          pFeature->m_pFace = pFace;
          pFeature->m_bBoss = ( lConcave != 0 ) ? TRUE : FALSE;
          for (ULONG e=0; e<sRim.GetSize(); e++)
            {
              pFeature->m_vRim.push_back( sRim[e] );
              SmCurve * pCurve = NULL;
              SER( kf_EdgeCurve( m_crContext, *sRim[e], pCurve ));
              pFeature->m_vRimCurves.Add( pCurve );
            }
          if ( !pFeature->m_bBoss )
            {
              pFeature->m_vWalls = vWalls;
              continue;
            }
          // A boss's walls continued down from its footprint: the footprint swept into the body,
          // along the planar face's normal, past the whole sheet.
          SmPoint3d  sAt;
          SmVector3d sNormal;
          SmBoolean  bPlanar = FALSE;
          SER( kf_PlaneOf( *pFace->GetSurface(), sAt, sNormal, bPlanar ));
          SmBoolean bUpright = bPlanar;
          for (SmFace * pWall : vWalls)
            {
              if ( bUpright )
                { SER( kf_Upright( *pWall, sNormal, dAngTol, bUpright )); }
            }
          SmPoint3d sInWall;
          SER( vWalls[0]->CalculateAnInternalPoint( sInWall ));
          pFeature->m_sUp = ( ( sInWall - sAt ).Dot( sNormal ) >= 0.0 ) ? sNormal : sNormal * -1.0;
          pFeature->m_vWalls.clear();
          if ( !bUpright )
            {
              // judged below, only if the band reaches it
              continue;
            }
          double dDepth = 0.0;
          for (ULONG c=0; c<8; c++)
            { dDepth = smos_Max( dDepth, ( sAt - asCorners[c] ).Dot( pFeature->m_sUp )); }
          dDepth += dSheetSize;
          // Each rim curve swept in two halves: every continued wall a patch whose crease with the
          // round runs from one of its sides to the other (a closed rim swept whole is a tube,
          // which its crease would go round).
          SmTArray<SmCurve*> sProfile;
          for (ULONG e=0; e<pFeature->m_vRimCurves.GetSize(); e++)
            {
              const SmCurve * pRim = pFeature->m_vRimCurves[e];
              const SmExtent1d sIvl = pRim->GetNaturalInterval();
              const double dMid = sIvl.Evaluate( 0.5 );
              for (ULONG h=0; h<2; h++)
                {
                  SmCurve * pHalf = NULL;
                  SER( pRim->Copy( m_crContext, pHalf ));
                  sProfile.Add( pHalf );
                  SmExtent1d sHalf;
                  SER( sHalf.SetMinMax( h == 0 ? sIvl.GetMin() : dMid, h == 0 ? dMid : sIvl.GetMax() ));
                  SER( pHalf->Trim( sHalf ));
                }
            }
          SmTArray<SmFace*> sBefore;
          m_pFilletBrep->GetFaces( sBefore );
          SmPrimitiveCreation sSweep( m_pFilletBrep->GetInfiniteRegion(), 0, TRUE );
          // the sweep leaves the Brep it builds in closed to editing: the sheet is edited on below
          const SmBoolean bEditing = m_pFilletBrep->m_bEditingEnabled;
          SmStatus eSwept = sSweep.CreateLinearSweep( sProfile, pFeature->m_sUp * -1.0, dDepth, 1, FALSE );
          m_pFilletBrep->m_bEditingEnabled = bEditing;
          SER( eSwept );
          SmTArray<SmFace*> sAfter;
          m_pFilletBrep->GetFaces( sAfter );
          for (ULONG w=0; w<sAfter.GetSize(); w++)
            {
              ULONG lAt = 0;
              if ( !sBefore.FindElement( sAfter[w], lAt ))
                { pFeature->m_vWalls.push_back( sAfter[w] ); }
            }
        }
    }
  if ( vFeatures.empty() )
    { return SM_SUCCESS; }

  // 2: the creases, each where a feature's wall meets a sheet face, kept on both
  std::vector<SmExtent3d> vSheetBoxes( sSheetFaces.GetSize() );
  for (ULONG s=0; s<sSheetFaces.GetSize(); s++)
    { SER( sSheetFaces[s]->CalculateBoundingBox( vSheetBoxes[s] )); }
  // per sheet face, its creases (owned by their features)
  std::map<SmFace*, SmTArray<SmCurve*> > sSheetCreases;
  // per continued wall of a boss, its creases
  std::map<SmFace*, SmTArray<SmCurve*> > sWallCreases;
  // per continued wall of a boss, the round's surface its creases lie on (when one)
  std::map<SmFace*, const SmSurface*> sWallRound;
  std::vector<kf_Feature*> vKept;
  for (kf_Feature * pFeature : vFeatures)
    {
      SmBoolean bUnbuilt = ( pFeature->m_bBoss && pFeature->m_vWalls.empty() ) ? TRUE : FALSE;
      if ( bUnbuilt )
        {
          // a boss whose walls are not upright: reached when its footprint's rim crosses a sheet face's box
          SmExtent3d sRimBox;
          for (ULONG e=0; e<pFeature->m_vRimCurves.GetSize(); e++)
            {
              const SmCurve * pRim = pFeature->m_vRimCurves[e];
              for (ULONG k=0; k<=8; k++)
                {
                  SmPoint3d sPnt;
                  SER( pRim->EvaluatePoint( pRim->GetNaturalInterval().Evaluate( k / 8.0 ), sPnt ));
                  sRimBox.AddPoint3d( sPnt );
                }
            }
          SmBoolean bNear = FALSE;
          for (ULONG s=0; s<sSheetFaces.GetSize() && !bNear; s++)
            { bNear = vSheetBoxes[s].AreDisjoint( sRimBox, dOnTol ) ? FALSE : TRUE; }
          if ( bNear )
            {
              GetFilletErrorInfo()->NoteError( SM_FILERR_NOT_EXACT, _T("Fillet Error: a boss inside a round's band, kept, whose walls are not upright to a planar face it stands on (its walls continued are not built yet)") );
              return SM_ERR;
            }
          continue;
        }
      for (SmFace * pWall : pFeature->m_vWalls)
        {
          SmExtent3d sWallBox;
          SER( pWall->CalculateBoundingBox( sWallBox ));
          for (ULONG s=0; s<sSheetFaces.GetSize(); s++)
            {
              if ( vSheetBoxes[s].AreDisjoint( sWallBox, dOnTol ))
                { continue; }
              kf_Curves sCreases;
              SER( kf_Creases( m_crContext, *sSheetFaces[s], *pWall, dCreaseTol, sCreases ));
              SmTArray<SmCurve*> sOnWall;
              SmObjsDelete<SmCurve*> sCleanOnWall( &sOnWall );
              for (ULONG c=0; c<sCreases.m_v3D.GetSize(); c++)
                { SER( kf_PiecesInFace( m_crContext, *sCreases.m_v3D[c], *pWall, dTol, sOnWall )); }
              SmTArray<SmCurve*> sOnBoth;
              for (ULONG c=0; c<sOnWall.GetSize(); c++)
                { SER( kf_PiecesInFace( m_crContext, *sOnWall[c], *sSheetFaces[s], dTol, sOnBoth )); }
              for (ULONG c=0; c<sOnBoth.GetSize(); c++)
                {
                  pFeature->m_vCreases.Add( sOnBoth[c] );
                  sSheetCreases[sSheetFaces[s]].Add( sOnBoth[c] );
                  if ( pFeature->m_bBoss )
                    {
                      sWallCreases[pWall].Add( sOnBoth[c] );
                      auto itRound = sWallRound.find( pWall );
                      if ( itRound == sWallRound.end() )
                        { sWallRound[pWall] = sSheetFaces[s]->GetSurface(); }
                      else if ( itRound->second != sSheetFaces[s]->GetSurface() )
                        { itRound->second = NULL; }
                    }
                }
            }
        }
      if ( pFeature->m_vCreases.GetSize() == 0 )
        { continue; }
      vKept.push_back( pFeature );
      // the event, in the trace
      SmPoint3d sAt;
      SER( pFeature->m_vRimCurves[0]->EvaluatePoint( pFeature->m_vRimCurves[0]->GetNaturalInterval().GetMin(), sAt ));
      TCHAR sMsg[SM_TBLOCK_SIZE];
      smos_sprintf( sMsg, _T("Fillet: keep features: a round's band reaches the footprint of a %s at (%.6g, %.6g, %.6g); its walls are carried to the round along %lu crease(s)"),
                    pFeature->m_bBoss ? _T("boss") : _T("cut"), sAt.x, sAt.y, sAt.z, static_cast<unsigned long>(pFeature->m_vCreases.GetSize()) );
      smos_ErrorMessage( SM_ERR_MESSAGE, _T(__FILE__), static_cast<ULONG>(__LINE__), sMsg );
    }
  // the continued walls of the bosses no band reaches leave the sheet
  for (kf_Feature * pFeature : vFeatures)
    {
      if ( !pFeature->m_bBoss || std::find( vKept.begin(), vKept.end(), pFeature ) != vKept.end() )
        { continue; }
      for (SmFace * pWall : pFeature->m_vWalls)
        { SER( m_pFilletBrep->DeleteFace( pWall, TRUE )); }
      pFeature->m_vWalls.clear();
    }
  if ( vKept.empty() )
    { return SM_SUCCESS; }

  // 3: each sheet face split along its creases; the pieces over a footprint, bounded by its
  // feature's creases alone, go
  for (ULONG s=0; s<sSheetFaces.GetSize(); s++)
    {
      SmFace * pSheet = sSheetFaces[s];
      auto it = sSheetCreases.find( pSheet );
      if ( it == sSheetCreases.end() )
        { continue; }
      SmTArray<SmFace*> sNewFaces;
      SmTArray<SmEdge*> sOnEdges;
      SER( m_pFilletBrep->MergeCurvesOnSurface( *pSheet->GetSurface(), dTol, it->second, NULL, sNewFaces, sOnEdges, FALSE, NULL ));
      SmTArray<SmFace*> sAll;
      sAll.Add( pSheet );
      sAll.Append( sNewFaces );
      for (ULONG p=0; p<sAll.GetSize(); p++)
        {
          SmExtent2d sDomain;
          SER( sAll[p]->CalculateUVDomainFromUVTrimCurves( sDomain, TRUE ));
          sAll[p]->SetUVDomain( sDomain );
          SmTArray<SmEdge*> sEdges;
          sAll[p]->GetEdges( sEdges );
          SmBoolean bOver = FALSE;
          for (kf_Feature * pFeature : vKept)
            {
              SmBoolean bAll = ( sEdges.GetSize() > 0 ) ? TRUE : FALSE;
              for (ULONG e=0; e<sEdges.GetSize() && bAll; e++)
                { SER( kf_EdgeAlong( *sEdges[e], pFeature->m_vCreases, dOnTol, bAll )); }
              if ( bAll )
                { bOver = TRUE; break; }
            }
          if ( bOver )
            { SER( m_pFilletBrep->DeleteFace( sAll[p], TRUE )); }
        }
    }
  // a boss's continued walls split along their creases, each keeping its piece up to the footprint
  for (kf_Feature * pFeature : vKept)
    {
      if ( !pFeature->m_bBoss )
        { continue; }
      for (SmFace * pWall : std::vector<SmFace*>( pFeature->m_vWalls ))
        {
          auto it = sWallCreases.find( pWall );
          if ( it == sWallCreases.end() )
            {
              GetFilletErrorInfo()->NoteError( SM_FILERR_BAD_MERGE, _T("Fillet Error: a kept boss's wall continued down meets no round") );
              return SM_ERR;
            }
          SmTArray<SmFace*> sNewFaces;
          SmTArray<SmEdge*> sOnEdges;
          SER( m_pFilletBrep->MergeCurvesOnSurface( *pWall->GetSurface(), dTol, it->second, NULL, sNewFaces, sOnEdges, FALSE,
                                                           sWallRound[pWall] ));
          SmTArray<SmFace*> sAll;
          sAll.Add( pWall );
          sAll.Append( sNewFaces );
          for (ULONG p=0; p<sAll.GetSize(); p++)
            {
              SmExtent2d sDomain;
              SER( sAll[p]->CalculateUVDomainFromUVTrimCurves( sDomain, TRUE ));
              sAll[p]->SetUVDomain( sDomain );
              SmTArray<SmEdge*> sEdges;
              sAll[p]->GetEdges( sEdges );
              SmBoolean bUpToFootprint = FALSE;
              for (ULONG e=0; e<sEdges.GetSize() && !bUpToFootprint; e++)
                { SER( kf_EdgeAlong( *sEdges[e], pFeature->m_vRimCurves, dOnTol, bUpToFootprint )); }
              if ( !bUpToFootprint )
                { SER( m_pFilletBrep->DeleteFace( sAll[p], TRUE )); }
            }
        }
    }

  // 4: the kept sheet's boundary imprinted on the target faces it lies on (a boss's creases and
  // its continued walls' tops at its footprint are left out: they meet no target face to cut)
  SmTArray<SmCurve*> sBossCreases;
  for (kf_Feature * pFeature : vKept)
    {
      if ( pFeature->m_bBoss )
        {
          sBossCreases.Append( pFeature->m_vCreases );
          sBossCreases.Append( pFeature->m_vRimCurves );
        }
    }
  SmTArray<SmEdge*> sSheetEdges;
  m_pFilletBrep->GetEdges( sSheetEdges );
  SmTArray<SmFace*> sTargetFaces;
  m_pTargetBrep->GetFaces( sTargetFaces );
  std::set<SmEdge*> sStops;
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
          SmBoolean bBossCrease = FALSE;
          SER( kf_EdgeAlong( *pEdge, sBossCreases, dOnTol, bBossCrease ));
          if ( bBossCrease )
            { continue; }
          SmBoolean bOn = FALSE;
          SER( kf_EdgeOnSurface( *pEdge, *pFace, dOnTol, bOn ));
          if ( !bOn )
            { continue; }
          SmCurve * pCopy = NULL;
          SER( kf_EdgeCurve( m_crContext, *pEdge, pCopy ));
          sOn.Add( pCopy );
        }
      if ( sOn.GetSize() == 0 )
        { continue; }
      SmTArray<SmFace*> sNewFaces;
      SmTArray<SmEdge*> sOnEdges;
      SER( m_pTargetBrep->MergeCurvesOnSurface( *pFace->GetSurface(), dTol, sOn, NULL, sNewFaces, sOnEdges, FALSE, NULL ));
      for (ULONG e=0; e<sOnEdges.GetSize(); e++)
        { sStops.insert( sOnEdges[e] ); }
    }
  // a kept boss stands on its footprint: the faces the blends take stop there
  SmTArray<SmEdge*> sTargetEdges;
  m_pTargetBrep->GetEdges( sTargetEdges );
  for (ULONG e=0; e<sTargetEdges.GetSize(); e++)
    {
      for (kf_Feature * pFeature : vKept)
        {
          SmBoolean bRim = FALSE;
          if ( pFeature->m_bBoss )
            { SER( kf_EdgeAlong( *sTargetEdges[e], pFeature->m_vRimCurves, dOnTol, bRim )); }
          if ( bRim )
            { sStops.insert( sTargetEdges[e] ); break; }
        }
    }

  // the faces round the blended edges, up to the stops: the material the blends take off
  std::vector<SmCurve*> vBlended;
  SmTArray<SmCurve*> sBlendedOwned;
  SmObjsDelete<SmCurve*> sCleanOwned( &sBlendedOwned );
  for (ULONG i=0; i<m_vFilletSolvers.GetSize(); i++)
    {
      SmEdge * pEdge = m_vFilletSolvers[i]->GetEdgeuse(0)->GetEdge();   NER(pEdge);
      SmCurve * pCopy = NULL;
      SER( kf_EdgeCurve( m_crContext, *pEdge, pCopy ));
      sBlendedOwned.Add( pCopy );
    }
  std::set<SmFace*> sGone;
  std::vector<SmFace*> vFront;
  for (ULONG e=0; e<sTargetEdges.GetSize(); e++)
    {
      SmEdge * pEdge = sTargetEdges[e];
      if ( sStops.count( pEdge ))
        { continue; }
      SmBoolean bBlended = FALSE;
      SER( kf_EdgeAlong( *pEdge, sBlendedOwned, dTol, bBlended ));
      if ( !bBlended )
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
          if ( sStops.count( sEdges[e] ))
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
      GetFilletErrorInfo()->NoteError( SM_FILERR_BAD_MERGE, _T("Fillet Error: a round kept round a feature does not close on the body (its rails and creases enclose no part of it)") );
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

  // the sheet and the bosses' continued walls sewn in, and the body oriented
  SmTArray<SmFace*> sKept;
  m_pFilletBrep->GetFaces( sKept );
  SER( m_pFilletBrep->CopyFaces( sKept, m_pTargetBrep, NULL, FALSE ));
  SmBrep * pSheet = m_pFilletBrep;
  m_pFilletBrep = NULL;
  delete pSheet;
  ULONG  lStitched = 0, lLamina = 0;
  double dVertexGap = 0.0, dEdgeGap = 0.0;
  SER( m_pTargetBrep->StitchFaces( dTol, lStitched, lLamina, dVertexGap, dEdgeGap ));
  SmBoolean bMaybeOpen = FALSE;
  if ( lLamina == 0 )
    { SER( m_pTargetBrep->OrientTrimmedSurfaces( TRUE, bMaybeOpen, FALSE )); }
  if ( lLamina != 0 || bMaybeOpen || !m_pTargetBrep->IsManifoldSolid() )
    {
      GetFilletErrorInfo()->NoteError( SM_FILERR_BAD_MERGE, _T("Fillet Error: a round kept round a feature does not sew into the body along its rails and creases") );
      return SM_ERR;
    }
  // each boss's wall and its continuation one face again: only the edges along its footprint go
  SmTArray<SmEdge*> sKeep;
  m_pTargetBrep->GetEdges( sTargetEdges );
  SmBoolean bAnyRim = FALSE;
  for (ULONG e=0; e<sTargetEdges.GetSize(); e++)
    {
      SmBoolean bRim = FALSE;
      for (kf_Feature * pFeature : vKept)
        {
          if ( pFeature->m_bBoss && !bRim )
            { SER( kf_EdgeAlong( *sTargetEdges[e], pFeature->m_vRimCurves, dOnTol, bRim )); }
        }
      if ( bRim )
        { bAnyRim = TRUE; }
      else
        { sKeep.Add( sTargetEdges[e] ); }
    }
  if ( bAnyRim )
    { SER( m_pTargetBrep->RemoveTopologicalEdgesAndVertices( &sKeep, NULL )); }
  rbKept = TRUE;
  return SM_SUCCESS;
} // end SmFilletExecutive::KeepFeaturesInBands
