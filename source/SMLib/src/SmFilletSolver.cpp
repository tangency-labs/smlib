// SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
// SPDX-License-Identifier: Apache-2.0

/**********************************************************************//**
* FILE NAME --- SmFilletSolver.cpp
* PURPOSE: Source code file for SmFilletSolver object.
**********************************************************************/

#include "StdAfx.h"

#ifndef __SMFILLETSOLVER_H__
#include <SmFilletSolver.h>
#endif

#ifndef __SMFILLETEXECUTIVE_H__
#include <SmFilletExecutive.h>
#endif

#ifndef __SMFILLETINTERSECTOR_H__
#include <SmFilletIntersector.h>
#include <SmFilletCorner.h>
#endif

#ifndef __SMFILLETSTANDARDSOLVER_H__
#include <SmFilletStandardSolver.h>
#endif

#include <SmLine.h>
#include <SmGraphicsExtern.h>
#include <SmGeomUtility.h>
#include <SmIsoCurve.h>
#include <SmCone.h>
#include <SmCylinder.h>
#include <SmSurfOfRevolution.h>
#include <SmPlane.h>
#include <SmFilletExact.h>
#include <SmRegion.h>
#include <SmCircle.h>
#include <SmSolutionArray.h>
#include <SmFace.h>
#include <SmTol.h>
#include <vector>

// Forward decl of local routine:
static SmBoolean sm_TestManifoldTangentRollover
  (SmEdge * pE, SmBoolean & rbOppositeNormals, double dTangencyTolRad, int iDebugLevel );


/*******************************************************************//**
PURPOSE: Constructor for the fillet solver object.  This is an
    abstract class.

NOTES:
***********************************************************************/
SmFilletSolver::SmFilletSolver
  (const SmContext & crContext,
   double            dThisApproxTol3d,
   double            dThisAngTolRad,
   double            dTangencyTolerance)
 : m_crContext              (crContext),
   m_dThisApproxTol3d       (dThisApproxTol3d),
   m_dThisAngTolRad         (dThisAngTolRad),       
   m_dTangencyTolerance     (dTangencyTolerance),
   m_dEdgeBlendFactor       (1.0),
   m_pFGForSurfaces         (NULL),
   m_dSurfaceExtensionFactor(0.5),
   m_dConversionTol         (SM_EFF_ZERO/10.0),
   m_bReverseTrim           (FALSE),
   m_pFSG                   (NULL),
   m_eStatus                (SM_FIL_UNPROCESSED)
{
    for (ULONG i=0; i<2; i++) 
      { m_pSurfaces[i] = NULL ; }

    SmFilletGeom * pFilletGeom = new SmFilletGeom(this);
    m_vFilletGeoms.Add(pFilletGeom);

} // end SmFilletSolver::SmFilletSolver


/*******************************************************************//**
PURPOSE: Destructor for the fillet solver object.  This is an
    abstract class.

NOTES:
***********************************************************************/
SmFilletSolver::~SmFilletSolver
  ()
{
  // delete filletGeoms
  for (ULONG i=0; i<m_vFilletGeoms.GetSize(); i++)
    {
      SM_ASSERT(m_vFilletGeoms[i] != NULL) ; 
      delete( m_vFilletGeoms[i] );
      m_vFilletGeoms[i] = NULL ;
    }

  // clear surfaces
  if (m_pSurfaces[0]) { delete m_pSurfaces[0]; m_pSurfaces[0] = NULL ; }
  if (m_pSurfaces[1]) { delete m_pSurfaces[1]; m_pSurfaces[1] = NULL ; }

  for (ULONG i=0; i<m_vLawCornerPieces.GetSize(); i++)
    { delete m_vLawCornerPieces[i]; m_vLawCornerPieces[i] = NULL; }

} // end SmFilletSolver::~SmFilletSolver destructor

/*******************************************************************//**
PURPOSE: SolidCraft change. The kept pieces, not matched yet, of a fillet
    surface that turns (a C0 knot across u) at each of their corners, NULL
    when there are none.
***********************************************************************/
SmStatus SmFilletSolver::FindLawCornerPieces
  (const SmSurface    & crSurface, // in : a fillet geom's surface
   SmLawCornerPieces *& rpPieces)  // out: its pieces, marked matched, or NULL
{
  rpPieces = NULL;
  SmTArray<double> sUKnots, sVKnots;
  SER(crSurface.FindKnotsByMaxContinuity(SM_CT_C0, FALSE, sUKnots, sVKnots));
  const SmExtent2d sDomain = crSurface.GetNaturalUVDomain();
  const double dTol = 1.0e-9 * (1.0 + smos_Fabs(sDomain.GetMin().x) + smos_Fabs(sDomain.GetMax().x));
  for (ULONG i=0; i<m_vLawCornerPieces.GetSize() && rpPieces == NULL; i++)
    {
      SmLawCornerPieces * pKept = m_vLawCornerPieces[i];
      if (pKept == NULL || pKept->m_bTaken || pKept->m_vJoinedEnds.GetSize() < 2)
        { continue; }
      // each corner (the pieces' inner ends) a turn of the surface; an extension at an end may
      // add a turn of its own where it joins
      SmBoolean bSame = TRUE;
      for (ULONG k=0; k+1<pKept->m_vJoinedEnds.GetSize() && bSame; k++)
        {
          SmBoolean bFound = FALSE;
          for (ULONG m=0; m<sUKnots.GetSize() && !bFound; m++)
            { bFound = (smos_Fabs(sUKnots[m] - pKept->m_vJoinedEnds[k]) <= dTol) ? TRUE : FALSE; }
          bSame = bFound;
        }
      if (bSame)
        {
          pKept->m_bTaken = TRUE;
          rpPieces = pKept;
        }
    }
  return SM_SUCCESS;

} // end SmFilletSolver::FindLawCornerPieces

/*******************************************************************//**
PURPOSE: SolidCraft change. Whether a round of this solver has a law corner
    that dips or peaks (built by SmFilletExecutive::BuildLawCorners).
***********************************************************************/
SmBoolean SmFilletSolver::HasBuiltLawCorners() const
{
  for (ULONG i=0; i<m_vLawCornerPieces.GetSize(); i++)
    {
      const SmLawCornerPieces * pKept = m_vLawCornerPieces[i];
      if (pKept != NULL && (pKept->m_vKinds.IsIn(SmLawCornerPieces::SM_LC_DIP) || pKept->m_vKinds.IsIn(SmLawCornerPieces::SM_LC_PEAK)))
        { return TRUE; }
    }
  return FALSE;

} // end SmFilletSolver::HasBuiltLawCorners

/*******************************************************************//**
PURPOSE: SolidCraft change. SmLawCornerPieces owns its pieces.
***********************************************************************/
SmLawCornerPieces::~SmLawCornerPieces()
{
  for (ULONG i=0; i<m_vPieces.GetSize(); i++)
    { delete m_vPieces[i]; m_vPieces[i] = NULL; }

} // end SmLawCornerPieces::~SmLawCornerPieces

/*******************************************************************//**
PURPOSE: SolidCraft change. A copy with copies of the pieces; the joined
    surface pointer is copied as it is.
***********************************************************************/
SmStatus SmLawCornerPieces::Copy
  (const SmContext     & crContext,  // in : context for the copies
   SmLawCornerPieces *& rpCopy) const // out: new copy
{
  rpCopy = NULL;
  SmLawCornerPieces * pCopy = new SmLawCornerPieces;
  pCopy->m_bTaken      = m_bTaken;
  pCopy->m_vJoinedEnds = m_vJoinedEnds;
  pCopy->m_vKinds      = m_vKinds;
  pCopy->m_vCentres    = m_vCentres;
  pCopy->m_vRadii      = m_vRadii;
  pCopy->m_dTol3d      = m_dTol3d;
  for (ULONG i=0; i<m_vPieces.GetSize(); i++)
    {
      SmSurface * pPiece = NULL;
      SmStatus eStatus = m_vPieces[i]->Copy(crContext, pPiece);
      SmBSplineSurface * pBSpline = SM_CAST_PTR(SmBSplineSurface, pPiece);
      if (eStatus != SM_SUCCESS || pBSpline == NULL)
        {
          delete pPiece;
          delete pCopy;
          return SM_ERR;
        }
      pCopy->m_vPieces.Add(pBSpline);
    }
  rpCopy = pCopy;
  return SM_SUCCESS;

} // end SmLawCornerPieces::Copy

/*******************************************************************//**
PURPOSE: Private method for base class, called by derived classes
    at end of their constructors.

NOTES: Installs offset surfaces into any FilletGeom's we have.
    Note, nobody here owns any base surfaces.
***********************************************************************/
SmStatus SmFilletSolver::FinishConstruction()
{
  ULONG lWhichFG, lWhichSrf;
  for ( lWhichFG = 0; lWhichFG < m_vFilletGeoms.GetSize(); lWhichFG++ )
  {
      for ( lWhichSrf = 0; lWhichSrf < 2; lWhichSrf++ )
      {
          SmOffsetSurface *pOurOffSrf = this->GetSurface( lWhichSrf );
          if ( pOurOffSrf == NULL )
              { continue; }

          // SmOffsetSurface::Copy() returns an SmSurface*, not SmOffsetSurface*
          SmSurface *pSrfCopy = NULL;
          SER( pOurOffSrf->Copy( *(pOurOffSrf->GetContext()), pSrfCopy ) );
          SmOffsetSurface *pOffsetCopy = SM_CAST_PTR(SmOffsetSurface, pSrfCopy);
          if ( pOffsetCopy != NULL )
          {
              SER( m_vFilletGeoms[lWhichFG]->SetOffsetSurface( lWhichSrf, pOffsetCopy ) );
          }
      }
  }
  return SM_SUCCESS;
} // end SmFilletSolver::FinishConstruction()

/*******************************************************************//**
PURPOSE: Retrieve the indicated offset surface.

NOTES:
***********************************************************************/
SmOffsetSurface * SmFilletSolver::GetSurface
 (ULONG lRailIndex) 
 const
{
  if ( m_pFGForSurfaces != NULL )
    {
      SmOffsetSurface *pRetVal = m_pFGForSurfaces->GetOffsetSurface(lRailIndex) ;
      if(pRetVal != NULL )
        { return pRetVal; }
    }

  return m_pSurfaces[lRailIndex];

} // end SmFilletSolver::GetSurface

/*******************************************************************//**
PURPOSE: Set the indicated offset surface.

NOTES:
***********************************************************************/
void SmFilletSolver::SetSurface(
        ULONG lRailIndex,
        SmOffsetSurface *pNewOffsetSurface
    )
{
    if ( m_pFGForSurfaces != NULL )
    {
        (void)m_pFGForSurfaces->SetOffsetSurface( lRailIndex, pNewOffsetSurface ); // SolidCraft: a void setter cannot report
    }
    else
    {
        m_pSurfaces[lRailIndex] = pNewOffsetSurface;
    }
}

/*******************************************************************//**
PURPOSE: Get the Self Intersection Handler of a fillet solver

NOTES: Will Append all fillet geoms to output array
***********************************************************************/
void SmFilletSolver::GetFilletGeoms
  (SmTArray<SmFilletGeom*> & rFilletGeoms)
{
    rFilletGeoms.ReSet();
    rFilletGeoms.Append(m_vFilletGeoms);

} // end SmFilletSolver::GetFilletGeoms

/*******************************************************************//**
PURPOSE: Get the Self Intersection Handler of a fillet solver

NOTES:
***********************************************************************/
SmSelfIntersectionHandler * SmFilletSolver::GetSelfIntersectionHandler
  ()
 const
{
    return m_pExecutive->m_pSelfIntersectionHandler;

} // end SmFilletSolver::GetSelfIntersectionHandler

/*******************************************************************//**
PURPOSE: Get the start (lIndex == 0) or the end (lIndex == 1) vertex
     of the fillet solver if they exist.  The start or end is defined
     by the orientation of the curve of the edge.  The start corresponds
     to the minimum parameter value of the curve.

NOTES:
***********************************************************************/
SmVertex * SmFilletSolver::GetVertex
  (ULONG lIndex)
 const
{
    SmEdgeuse *pEU = GetEdgeuse(0);
    if (!pEU) return NULL;
    SmEdge *pE = pEU->GetEdge();
    if (!pE) { SE(SM_ERR); return NULL; }
    SmVertex *pStartV = pE->GetStartVertex();
    if (lIndex == 0) return pStartV;
    return pE->GetOtherVertex(pStartV);

} // end SmFilletSolver::GetVertex

/*******************************************************************//**
PURPOSE: Get the start (lIndex == 0) or the end (lIndex == 1) FilletCorner
     of the fillet solver if they exist.  The start or end is defined
     by the orientation of the curve of the edge.  The start corresponds
     to the minimum parameter value of the curve.

NOTES:
***********************************************************************/
SmFilletCorner *SmFilletSolver::GetFilletCorner( ULONG lIndex ) const
{
  if ( lIndex > 1 ) { return NULL; }
  SmVertex * pVtx = this->GetVertex( lIndex );
  return this->GetFilletExecutive()->GetFilletCornerOfVertex( pVtx );

} // end SmFilletSolver::GetFilletCorner

/*******************************************************************//**
PURPOSE: Get the other FilletCorner than the one passed in.

NOTES: Will be the same as pThisCorner, iff this FilletSolver is closed.
   Will be Null if pThisCorner is not one of our FilletCorners.
***********************************************************************/
SmFilletCorner *SmFilletSolver::GetOtherFilletCorner( const SmFilletCorner *pThisCorner ) const
{
  SmFilletCorner *pCorner0 = this->GetFilletCorner( 0 );
  SmFilletCorner *pCorner1 = this->GetFilletCorner( 1 );
  if ( pThisCorner == pCorner0 )
    { return pCorner1; }
  if ( pThisCorner == pCorner1 )
    { return pCorner0; }

  // Error: pThisCorner is not one of ours.
  // SM_ASSERT( SM_ERR_INVALID_INPUT );
  return NULL;

} // end SmFilletSolver::GetOtherFilletCorner

/*******************************************************************//**
PURPOSE: Find the index (0-1) of a rail curve which lies on a given face.

NOTES:
***********************************************************************/
ULONG SmFilletSolver::FindIndexOfRailOnFace
  (SmFace * pOrigFace)
{
    SmEdgeuse * pEU = GetEdgeuse(0);
    if (!pEU) {
        SE(SM_ERR);
        return 0;
    }
    if (pEU->GetLoopuse()->GetFaceuse()->GetFace() == pOrigFace)
        return 0;

#ifdef SM_DEBUG_CODE
    if ( DebugLevel() > 20 ) {
        pOrigFace->GetSurface()->DrawUV(2,2);
        smgfx_SetColor(1,0,0);
        pEU->Draw();
        smgfx_SetColor(0,1,0);
        GetEdgeuse(1)->Draw();
        sm_GraphicsLoop();
    }
#endif

    return 1;
} // end ULONG SmFilletSolver::FindIndexOfRailOnFace


/*******************************************************************//**
PURPOSE: Find the index (0-1) of a rail curve which intersects with
    side edgeuse.

NOTES:
***********************************************************************/
ULONG SmFilletSolver::FindIndexOfRailXSideEdgeuse
    ( SmEdgeuse * pSideEU )  // Side edgeuse that intersects RAIL
{
  ULONG lRailIndex = 0;
  SmEdgeuse * pFilEdgeuse = GetEdgeuse(0);
  if (pSideEU->GetEdge()->IsClosed())
    { // Such as when in 2x1 cases
      SmVertex * pV     = pFilEdgeuse->GetVertexuse()->GetVertex();
      SmVertex * pSideV = pSideEU->GetVertexuse()->GetVertex();

      if (pSideV == pV)
        {
          if (pFilEdgeuse->GetCWEdgeuse() == pSideEU)
            {
              lRailIndex = 1;
            }
        }
      else
        { //pSideV != pV
          if (pFilEdgeuse->GetCCWEdgeuse()->GetMate() == pSideEU)
            {
              lRailIndex = 1;
            }
        }
    }
  else
    {
      // Work from the faces of our two Edgeuses.
      lRailIndex = 0;
      if (GetEdgeuse(0)->GetFace() != pSideEU->GetFace())
        {
          //if (GetEdgeuse(0)->GetLoopuse() != pSideEU->GetLoopuse()) {}
          lRailIndex = 1;
        }
    }

  return lRailIndex;
} // end ULONG SmFilletSolver::FindIndexOfRailXSideEdgeuse


/*******************************************************************//**
PURPOSE: This method sets the surface offsets back to zero making
    them equivalent to the base surfaces.

NOTES:
***********************************************************************/
void SmFilletSolver::SetSurfacesToZeroOffset
  ()
{
    if (GetSurface(0)) {
        m_dSavedOffsets[0] = GetSurface(0)->GetOffsetDistance();
        GetSurface(0)->SetOffsetDistance(0.0);
    }
    if (GetSurface(1)) {
        m_dSavedOffsets[1] = GetSurface(1)->GetOffsetDistance();
        GetSurface(1)->SetOffsetDistance(0.0);
    }
} // end void SmFilletSolver::SetSurfacesToZeroOffset

/*******************************************************************//**
PURPOSE: Reload the surface offset values from the saved offsets.

NOTES:
***********************************************************************/
void SmFilletSolver::ReloadSurfaceOffsets
  ()
{
    if (GetSurface(0)) {
        GetSurface(0)->SetOffsetDistance(m_dSavedOffsets[0]);
    }
    if (GetSurface(1)) {
        GetSurface(1)->SetOffsetDistance(m_dSavedOffsets[1]);
    }
} // end void SmFilletSolver::ReloadSurfaceOffsets

/****************************************************************
PURPOSE: Find a tsect point on fillet which is obtained from solving
            a 3D point on a cross-sectional plane

NOTES:
  Given a FilletSolver container for the edge and surfaces being filleted
        and an initial guess point,
  Find
****************************************************************/
static SmStatus sm_FindFilletTsectPnt
  (SmFilletSolver  * pFilSolver,      // in : Fillet Solver for target to-be-filleted edge
   const SmSurface * cpSurf1,         // in : Base Surface 1 connected to to-be-filleted edge
   const SmSurface * cpSurf2,         // in : Base Surface 2 connected to to-bo-filleted edge
   SmPoint3d       & rPntToDrop,      // in : guess point as close to desired point on fillet center curve as possible
   SmVector3d      & rTraceDir,       // in : fillet center curve tangent at guess point
   SmTsectPnt      & rNewTsectPnt)    // out: contains fillet intersect UVpnt values
                                      //        baseSurface(uvPnt)   = points on rail
                                      //        OffsetSurface(uvPnt) = point on fillet center curve and
                                      //                               plane = [rPntToDrop, rTraceDir]
{
  // check state
  SM_ASSERT(pFilSolver->GetEdgeuse(0)->GetFace()->GetSurface() == cpSurf1) ;
  SM_ASSERT(pFilSolver->GetSurface(0)->GetBaseSurface()        == cpSurf1) ;
  SM_ASSERT(pFilSolver->GetEdgeuse(1)->GetFace()->GetSurface() == cpSurf2) ;
  SM_ASSERT(pFilSolver->GetSurface(1)->GetBaseSurface()        == cpSurf2) ;

#ifdef SM_DEBUG_CODE
SmBoolean bDebug = FALSE;
  if (bDebug) {
      smgfx_SetColor(1,0,0);
      smgfx_SetPointSize(10);
      rPntToDrop.Draw();
      rTraceDir.Draw(&rPntToDrop);
      smgfx_SetColor(0,1,0);
      if (FALSE) {
          cpSurf1->DrawUV(5,5);
          cpSurf2->DrawUV(5,5);
      }
      sm_GraphicsLoop();
  }
#endif
  // locals
  double          dThisApproxTol3d = pFilSolver->GetThisApproxTol3d();
  SmExtent2d      sDomain1          = cpSurf1->GetNaturalUVDomain();
  SmExtent2d      sDomain2          = cpSurf2->GetNaturalUVDomain();
  SmSurfaceCache *pSC1              = smsurf_GetSurfaceCache(cpSurf1); NER(pSC1);
  SmSurfaceCache *pSC2              = smsurf_GetSurfaceCache(cpSurf2); NER(pSC2);

  SmCacheCheckOutIn sCheckIO1(pSC1);
  SmCacheCheckOutIn sCheckIO2(pSC2);

  // set scope for temporary state bit changes - state is restored when scope is exited
  {
    // Turn off SurfaceCache1 point testing so GlobalPointSolve() will keep all point solutions
    //   without classifying the solution point against the trim boundaries.
    // Turn on boundary curve processing to force LocalSolve to look for drop points
    //   on boundary curves.  This will find solutions where the curve comes close
    //   to the surface but does not actually intersect it.
    SmTemporaryChangeValue<SmBoolean> sStack10(pSC1->m_bPointTestEnabled,FALSE);
    SmTemporaryChangeValue<SmBoolean> sStack11(pSC1->m_bProcessBoundaryCurves,TRUE);

    // locals
    SmSolution      sSData[8];
    SmSolutionArray sSolutions(8,sSData);

    // Drop guess point to cpSurf1
    SER(cpSurf1->GlobalPointSolve(sDomain1, SM_SO_MINIMIZE,
                                  rPntToDrop, dThisApproxTol3d, NULL ,
                                  SM_SR_ALL, sSolutions));
    if (sSolutions.GetSize() < 1) { SER(SM_ERR); }

    // extract surface1 UVPoint
    SmPoint2d sUV1(sSolutions[0].m_vStart[0],sSolutions[0].m_vStart[1]);

    // Turn off SurfaceCache2 point testing so GlobalPointSolve() will keep all point solutions
    //   without classifying the solution point against the trim boundaries.
    // Turn on boundary curve processing to force LocalSolve to look for drop points
    //   on boundary curves. This will find solutions where the drop point comes close
    //   to the surface but does not have a surface normal vector which intersects it.
    SmTemporaryChangeValue<SmBoolean> sStack20(pSC2->m_bPointTestEnabled,FALSE);
    SmTemporaryChangeValue<SmBoolean> sStack21(pSC2->m_bProcessBoundaryCurves,TRUE);

    // Drop guess point to cpSurf2
    SER(cpSurf2->GlobalPointSolve(sDomain2,SM_SO_MINIMIZE,
        rPntToDrop,dThisApproxTol3d,NULL,SM_SR_ALL,sSolutions));
    if (sSolutions.GetSize() < 1) { SER(SM_ERR); }

    // extract surface2 UVPoint
    SmPoint2d sUV2(sSolutions[0].m_vStart[0],sSolutions[0].m_vStart[1]);

    SmBoolean bFoundSolution;

    // Find surface uvPnts on plane = [rPntToDrop, rTraceDir]
    //   such that BaseSurface(uvPnt)   = points on rail curves
    //     and     OffsetSurface(uvPnt) = point on fillet center curve and given plane
    SER(pFilSolver->PointOnPlaneSolve(rPntToDrop, rTraceDir,        // plane
                                      sDomain1,   sDomain2,         // domain search limits
                                      sUV1,       sUV2,             // guess points
                                      bFoundSolution,rNewTsectPnt));
    if (!bFoundSolution) { return SM_ERR; }

    // store trace direction in data structure
    rNewTsectPnt.CrvDeriv() = rTraceDir;

#ifdef SM_DEBUG_CODE
SmBoolean bDebugMe1 = FALSE;
    if (bDebugMe1) {
        sm_GraphicsLoop();
        smgfx_SetColor(1,0,0);
        sm_GraphicsLoop();
        smgfx_SetPointSize(6.0);
        SmPoint3d sPnt;
        smgfx_SetColor(1,0,0);
        cpSurf1->EvaluatePoint(rNewTsectPnt.UVPos(0),sPnt);
        sPnt.Draw();
        smgfx_SetColor(1,0,1);
        cpSurf2->EvaluatePoint(rNewTsectPnt.UVPos(1),sPnt);
        sPnt.Draw();
        sm_GraphicsLoop();
    }
#endif
  }

  // all done
  return SM_SUCCESS;

} // end sm_FindFilletTsectPnt

#ifdef SM_DEBUG_CODE
/****************************************************************                                 \
 PURPOSE: static helper to draw an SmTsectPnt, checking for unset values.

NOTES:
****************************************************************/
static void sm_DrawTsectPnt( const SmTsectPnt * cpTSP, const SmSurface *pSrf0, const SmSurface *pSrf1 )
{
  // The basic Draw method is ok if things are initialized.
  if ( cpTSP->m_v3DCurvePV[0].IsInitialized() )
  {
      cpTSP->Draw();
      return;
  }

  // If the surface uv points are defined, draw them.
  if ( cpTSP->m_vUVCurvePV[0][0].IsInitialized() )
  {
      SmPoint2d sPosUV = cpTSP->m_vUVCurvePV[0][0];
      SmPoint3d sPos3d;

      // Surfaces are probably offset surfaces.
      // If so, draw the base surface points.
      const SmOffsetSurface *pOffSrf = SM_CAST_PTR( SmOffsetSurface, pSrf0 );
      if ( pOffSrf != NULL )
      {
          const SmSurface *pBaseSrf = pOffSrf->GetBaseSurface();
          pBaseSrf->EvaluatePoint( sPosUV, sPos3d );
          sPos3d.Draw(); sm_GraphicsLoop();
          sm_GraphicsLoop();
      }
      pSrf0->EvaluatePoint( sPosUV, sPos3d );
      sPos3d.Draw(); sm_GraphicsLoop();
      sm_GraphicsLoop();


      sPosUV = cpTSP->m_vUVCurvePV[1][0];

      pOffSrf = SM_CAST_PTR( SmOffsetSurface, pSrf1 );
      if ( pOffSrf != NULL )
      {
          const SmSurface *pBaseSrf = pOffSrf->GetBaseSurface();
          pBaseSrf->EvaluatePoint( sPosUV, sPos3d );
          sPos3d.Draw(); sm_GraphicsLoop();
          sm_GraphicsLoop();
      }

      sPosUV = cpTSP->m_vUVCurvePV[1][0];
      pSrf1->EvaluatePoint( sPosUV, sPos3d );
      sPos3d.Draw(); sm_GraphicsLoop();
      sm_GraphicsLoop();
  }

} // end sm_DrawTsectPnt
#endif // SM_DEBUG_CODE

/*******************************************************************//**
PURPOSE: Calculate the fillet surface geometry of a
            given FilletSolver(edge)

NOTES:

METHOD ---
  1. build an ordered array of sample points stored as SmTsectPnts each built from
     a starter point on the to-be-filleted edge.  The starter point set includes:

     a. to-be-filleted edge endPoints (already computed and stored in SmFilletVertexuse objects)
     b. each g1 continuity point on the to-be-filleted edge
     c. when bSplitOriginalFillet == TRUE
         A. split point (snapped to internal G1 point when possible)

        where each SmTsectPnt contains
           1. fillet center-curve position and tangent
           2. associated surface uvPnts where
                baseSurface(uvPnt)   = points on rail curves
                offsetSurface(uvPnt) = point on fillet center-curve

  2. Clean up the array of SmTsectPnts so that endPoints are not too close
     to the 1st internal points

  3. Use ordered sample point sequence to build a fillet surface
     that interpolates the sampled rail curve positions and
     has the cross-section shape specified by the derived SmFilletSurfaceGenerator
     object stored in m_pFSG.

***********************************************************************/
SmStatus SmFilletSolver::CalcFilletGeom()
{
  // Set the extension flags for two ends of fillets
  // extend ends when vertices are not blended and not degenerate
  SetExtensionFlags();

  // local to-be-filleted edge data
  SmEdge         * pEdge      = GetEdgeuse(0)->GetEdge(); NER(pEdge);
  SmBSplineCurve * pEdgeCurve = SM_CAST_PTR(SmBSplineCurve, pEdge->GetCurve());
  NER(pEdgeCurve);
  SmExtent1d       sIvl       = pEdge->GetInterval();
  ULONG ii;

  // evaluate curve end points
  SmPoint3d sPVSt[2], sPVEn[2];
  SER(pEdgeCurve->Evaluate(sIvl.GetMin(),1,TRUE,sPVSt));
  SER(pEdgeCurve->Evaluate(sIvl.GetMax(),1,TRUE,sPVEn));

  // get fillet offset and base surfaces
  SmOffsetSurface * pOffsetSrf1 = GetSurface(0);
  SmOffsetSurface * pOffsetSrf2 = GetSurface(1);
  const SmSurface * pBaseSrf1  = pOffsetSrf1->GetBaseSurface(); NER(pBaseSrf1);
  const SmSurface * pBaseSrf2  = pOffsetSrf2->GetBaseSurface(); NER(pBaseSrf2);

#ifdef SM_DEBUG_CODE
int iDebugLevel = DebugLevel();
static ULONG lCount      = 1 ; lCount++ ;
static ULONG lDebugCount = 0 ;
TCHAR sBuff[SM_TBLOCK_SIZE];

  // draw Brep, to-be-filleted curve(red), start point
  if ( iDebugLevel > 0 || lCount == lDebugCount)
    {
      SmBrep *pBrep = pEdge->GetBrep() ;
      pEdge->Dump() ;

      smgfx_Erase() ;
      smgfx_SetLook(1,2, 0,0,0) ; if(pBrep) { pBrep->Draw(TRUE) ; } sm_GraphicsLoop() ;
      smgfx_SetLook(2,3, 0,1,1) ; pEdge->Draw() ; sm_GraphicsLoop() ;
      smgfx_SetLook(4,6, 1,0,0) ; sPVSt[0].Draw() ; sm_GraphicsLoop() ;
      smgfx_SetLook(4,6, 0,0,1) ; sPVEn[0].Draw() ; sm_GraphicsLoop() ;
      if ( iDebugLevel > 20 ) {
          smgfx_SetLook(1,2, 0,1,0) ; pOffsetSrf1->Draw() ; sm_GraphicsLoop() ;
          smgfx_SetLook(1,2, 1,1,0) ; pOffsetSrf2->Draw() ; sm_GraphicsLoop() ;
      }
      sm_GraphicsLoop();
    }
#endif // SM_DEBUG_CODE

  // get fillet geometry pointer and locals
  SmFilletGeom * pFilletGeom0         = m_vFilletGeoms[0];

#ifdef SM_DEBUG_CODE
static int sbDumpFVUs=0;
  if ( sbDumpFVUs>0 )
    { pFilletGeom0->DumpFilletVUs(); }
#endif // SM_DEBUG_CODE

  SmBoolean      bSplitOriginalFillet = pFilletGeom0->GetSplitFlag();
  SmTsectPnt   * pSplitTsectPnt = NULL;

  // locals for offsetSurface intersection points and temporary point array
  SmTsectPnt              * sP1Data[64];
  SmTArray<SmTsectPnt*>     sTsectPnts(64,sP1Data);
  SmTsectPnt              * bData[50];
  SmTArray<SmTsectPnt*>     sDeletePnts(50,bData);
  SmObjsDelete<SmTsectPnt*> sObjs(&sDeletePnts);

  // 1. build an ordered array of SmTsectPnts
  //    including one for each of the following:
  //   1. to-be-filleted edge Rail endpoints
  //      (already computed and stored in SmFilletVertexuse objects)
  //   2. internal G1 points on the to-be-filleted curve
  //   3. when bSplitOriginalFillet == TRUE
  //          A split point (snapped to internal G1 point when possible)
  //
  //   Points are ordered in uv space, by the parameter value along the
  //   edge->edgeuse UVTrimCurves.
  //
  // We can sort along either rail.  The only problem is if the
  // rail is closed *in uv space*.  (Don't care whether closed in 3d space.)
  // Iter both sector-bounding edgeuses,  break after 1st success.
  // If the first is closed, continue to the next.
  // If the first is not closed, sort along it and forget the other one.
  // If they're both closed, then don't sort at all,
  //   things get confused at the endPoints.

  SmBoolean bDoSorting = TRUE;
  ULONG lWhichSurf;
  for ( lWhichSurf=0; lWhichSurf<=1; lWhichSurf++ )
    {
      // init array being built
      sTsectPnts.ReSet();

      // get edgeuse->UVTrimCurve marking a rail curve and bracketing sector to be filleted
      SmEdgeuse      * pEUFilletedEdge = GetEdgeuse( lWhichSurf );
      NER(pEUFilletedEdge);
      SmBSplineCurve * pUVCurveFilletedEdge = NULL ;
      SER(pEUFilletedEdge->GetOrCreateUVTrimCurve( pUVCurveFilletedEdge ));

      // special case: don't sort SmTsectPnts on closed UVTrimCurves
      if ( pUVCurveFilletedEdge->IsClosed( sIvl ) )
        {
          if (lWhichSurf == 0) { continue; // try the other one.
                               }
          else                 { // MSG(_T("Possible sorting problem"));
                                  bDoSorting = FALSE; // both closed
                               }
        }

      // array of all TsectPnts (including potential duplicates):
      // from rail ends, and possibly G1 discontinuity points and split point.
      SmTArray<SmTsectPnt*> sAllTsectPnts;

      // when asked to split the edge - find and save a split parameter
      double    dSplitParam     = sIvl.Evaluate(0.5);
      SmBoolean bFoundSplitKnot = FALSE;
      if (bSplitOriginalFillet)
        {
          ULONG lTotalG1Knots = m_vG1Knots.GetSize();
          if (lTotalG1Knots > 0)
            {
              // Pick a middle G1 knot and split it right there
              dSplitParam = m_vG1Knots[lTotalG1Knots/2];
              bFoundSplitKnot = TRUE;
            }
          else
            {
              m_vG1Knots.Add(dSplitParam);
            }
        } // end need to split check

      // SolidCraft change: the corners of the radius law are through points too, so the trace
      // has a point where the fillet turns.
      SmTArray<double> sThroughKnots;
      sThroughKnots.Append(m_vG1Knots);
      SER(GetRadiusCornerParameters(sThroughKnots));

      // SolidCraft change: a rail that ends on an extended edge (SM_FV_RAIL_X_EXTENDED_EDGEUSE)
      // has its end found on the fillet surface once that is made, so it gives the trace no point
      // (see the rail ends below). The corner's cap edges meet at the filleted edge's end there
      // (SmFilletNx1Corner::MakeCornerTopology), and that rail ends where the cap edge on its face
      // leaves that end: the fillet runs at least to the edge's end, its cross-section there is a
      // through point, and the end's extension runs on past it. The other rail can end well
      // short of the edge's end (a round along a boss's base ending at a wall the boss overhangs).
      for ( ULONG lRail=0; lRail<2; lRail++ )
        {
          SmFilletEdge * pRail = pFilletGeom0->GetRail( lRail );
          NER(pRail);
          SmFilletVertex * pRailEnds[2] = { (SmFilletVertex*)pRail->GetVertex(), NULL };
          if ( pRailEnds[0] != NULL )
            { pRailEnds[1] = (SmFilletVertex*)pRail->GetOtherVertex( pRailEnds[0] ); }
          for ( ULONG lEnd=0; lEnd<2; lEnd++ )
            {
              SmFilletVertex * pFV = pRailEnds[ lEnd ];
              if (   pFV == NULL
                  || pFV->GetFilletVertexType() != SM_FV_RAIL_X_EXTENDED_EDGEUSE
                  || pFV->GetFilletCorner() == NULL )
                { continue; }
              const SmVertex * cpCornerVertex = pFV->GetFilletCorner()->GetFilletedVertex();
              if ( pEdge->GetStartVertex() == pEdge->GetEndVertex() )
                { continue; }  // a closed edge has no end the fillet stops at
              double dEndParam =   ( cpCornerVertex == pEdge->GetStartVertex() ) ? sIvl.GetMin()
                                 : ( cpCornerVertex == pEdge->GetEndVertex() )   ? sIvl.GetMax()
                                 :                                                 SM_BIG_DOUBLE;
              if ( dEndParam != SM_BIG_DOUBLE && ! sThroughKnots.IsIn( dEndParam ) )
                { sThroughKnots.Add( dEndParam ); }
            }
        }

      // Now, compute 'through points' from all internal G1 knots
      // and add them to sAllTsectPnts (and sDeletePnts arrays to stop mem leaks)
      for (ULONG k=0; k<sThroughKnots.GetSize(); k++)
        {
          // get position and tangent at G1 knot
          SmPoint3d sPV[2];
          double dKnot = sThroughKnots[k];
          SER(pEdgeCurve->Evaluate(dKnot,1,TRUE,sPV));

          // locals
          SmPoint3d  sPnt;
          SmVector3d sBinVec;
          SmVector3d sTangent;
          SmVector2d sUVParameter;
          SmEdgeuse *pEUOnSurf = GetEdgeuse(0);

          // get binormal pointing into the current edgeuse->surface
          //     tangent  pointing in edgeuse direction
          SER(pEUOnSurf->EvaluateBinormal(dKnot,FALSE,
                                          sPnt, sBinVec,&sTangent,
                                          NULL,&sUVParameter));

          // set binormal length = 1/2* 2nd surface offset distance
          double dOffset = smos_Fabs(GetSurface(1)->GetOffsetDistance());
          sBinVec = sBinVec*dOffset*0.5;
          SmVector2d sUVDir;

          // set guess point to Surf1 point nearest to BinVec endPoint
          SER(pBaseSrf1->DropVectors(sUVParameter,TRUE,TRUE,1,&sBinVec,&sUVDir));
          sUVParameter = sUVParameter+sUVDir;
          SER(pOffsetSrf1->EvaluatePoint(sUVParameter,sPnt));

          // Build a SmTsectPnt for each to-be-filleted curve g1 points.
          // Find fillet surface uvPnts for fillet center-curve sPnt and
          // sTangent values.  [sPnt, sTangent] = position and tangent of
          // best-guess for next fillet center-curve point.
          //   baseSurface(uvPnt)   = points on fillet rail curves
          //   offsetSurface(uvPnt) = point on fillet center curve and
          //                          on plane = [sPnt, sTangent]
          SmTsectPnt * pNewPnt = new SmTsectPnt;
          if (sm_FindFilletTsectPnt(this,
                                       pBaseSrf1,pBaseSrf2,
                                       sPnt,sTangent,
                                       *pNewPnt) == SM_SUCCESS)
            {
              // add new TsectPnt to sAllTsectPnts array
              sAllTsectPnts.Add(pNewPnt);
              sDeletePnts.Add(pNewPnt);
            }
          else
            { SER(SM_ERR);
            }

          // when splitting original fillet - save the appropriate TsectPnt
          if (   bSplitOriginalFillet
              && smos_Fabs(dKnot-dSplitParam) < SM_EFF_ZERO)
            {
              pSplitTsectPnt = pNewPnt;
            }
        } // end iter every G1 internal knot computing through points

      // restore G1Knots if a split param was added
      if (bSplitOriginalFillet && !bFoundSplitKnot)
        {
          m_vG1Knots.RemoveLast();
        }

      // Process rail_ends -
      // There are two rails, and two edgeuses on each rail
      //   (one for the base surface and one for the fillet surface).
      // Place rail[i]->edgeuse[jj]->FilletVertexuse->TsectPnt into sAllTsectPnts array
      //   Each rail end intersects with some stopping geometry.
      //   In general those 4 locations will happen at different filletEdge parameters
      //   note: in some special cases the railEnd intersection parameters at
      //         one end, the other end, or both ends may match up.
      //   note: The SmTsectPnt's fetched from the edgeuse->filletVertexuse
      //         corners are for the intersection of the offset surfaces:
      //         their uv values correspond to dropping from the spine curve
      //         to the base surfaces.
      //         They have only their surfUVPnt positions set.  All other
      //         SmTsectPnt data is uninitialized, including - 3dCurve position
      //         and tangent, surfUVPoint tang, and all 3dSurf value data.

      ULONG lWhichRail;
      for ( lWhichRail=0; lWhichRail<2; lWhichRail++ )
        {
          SmFilletEdge * pThisRail= pFilletGeom0->GetRail( lWhichRail );
          NER(pThisRail);
          SmTArray<SmEdgeuse*> sFilEdgeuses;
          SER( pThisRail->GetEdgeuses(sFilEdgeuses) );
          if (sFilEdgeuses.GetSize() != 2) { SER(SM_ERR); }

          // for both Rail->edgeuses
          for (ULONG thisEdgeuse=0; thisEdgeuse<2; thisEdgeuse++)
            {
              // from edgeuse get FilletVertexuse and FilletVertex
              //   the FilletVertexUse stores the SmTsectPoint for the
              //   to-be-filleted edge endPoints.
              SmEdgeuse         * pEU = sFilEdgeuses[ thisEdgeuse ];
              SmFilletVertexuse * pVU = (SmFilletVertexuse*)pEU->GetVertexuse();
              SmFilletVertex    * pFV = (SmFilletVertex*)pVU->GetVertex();

              // skip: mates just tag along with their mates, same param value.
              if (  pFV->GetFilletVertexType() == SM_FV_MATE
                 || pFV->GetFilletVertexType() == SM_FV_RAIL_X_EXTENDED_EDGEUSE)
                { continue; }

              // fetch TsectPnt from SmFilletVertexuse (already computed)
              SmTsectPnt & rTsectPnt = pVU->GetTsectPnt();
              if (pEU->GetOrientation() == SM_OT_SAME)
                {
                  sAllTsectPnts.InsertAt(0,&rTsectPnt);
                }
              else // SM_OT_OPPOSITE
                {
                  sAllTsectPnts.Add(&rTsectPnt);
                }

              // for both components convert numbers like 1.1200000000007 to 1.12
              for (ii=0; ii<2; ii++)
                {
                  rTsectPnt.UVPos(ii).x = smos_CleanUpNoise( rTsectPnt.UVPos(ii).x );
                  rTsectPnt.UVPos(ii).y = smos_CleanUpNoise( rTsectPnt.UVPos(ii).y );
                }
            } // end iter two edgeuses on this rail
        } // end iter both rails

      if ( bDoSorting )   //cbi: 090219: try sorting in 3d.
        {
          double dScale = smos_Fabs(sIvl.GetMin()) + smos_Fabs(sIvl.GetMax()) + 1.0;
          SmTArray<double> sProjParams;
          for (ULONG j=0; j<sAllTsectPnts.GetSize(); j++)
            {
              SmTsectPnt & rTsectPnt = *sAllTsectPnts[j];
              SmPoint3d sPnt;
              SmPoint2d sUV = rTsectPnt.UVPos(lWhichSurf);
              SmPoint3d sUVPnt(sUV);

              // skip points set to -SM_BIG_DOUBLE: we marked them below
              // in this loop, if we inserted them into the sTsectPnts list
              if (   sUV.x == -SM_BIG_DOUBLE
                  || sUV.y == -SM_BIG_DOUBLE)
                {
                  continue;
                }

              // Drop surfacePoint sUVPnt to pUVCurveFilletedEdge - set dT = drop point curve parameter
              SmSolution sSData[8];
              SmSolutionArray sSolutions(8,sSData);
              SER( pUVCurveFilletedEdge->GlobalPointSolve( sIvl, SM_SO_MINIMIZE,
                                                           sUVPnt, m_dThisApproxTol3d, 
                                                           NULL, NULL, SM_SR_SINGLE,
                                                           sSolutions ));
              if (sSolutions.GetSize() != 1) SER(SM_ERR);
              double dT = sSolutions[0].m_vStart[0];

              // when sUVPnt drops to filletEdge endPoint -
              //   it might be beyond edge interval end and edge needs to be extended
              if (sIvl.IsValueOnBoundary(dT))
                {
                  // set dT to sUVPnt drop to filletEdge-tangent-extension-line
                  SmPoint3d sPV[2];
                  SER(pUVCurveFilletedEdge->Evaluate(dT,1,1,sPV));
                  double dAdjustment;
                  SER(smgu_LineClosestPoint(sPV[0],sPV[1],sUVPnt,dAdjustment));
                  dT += dAdjustment;
                }

#ifdef SM_DEBUG_CODE
              if ( iDebugLevel > 0 ) 
                {
                  SmPoint3d sP1, sP2;
                  const SmSurface * pBaseSurf1 = GetSurface(lWhichSurf)->GetBaseSurface();
                  const SmSurface * pBaseSurf2 = GetSurface(1-lWhichSurf)->GetBaseSurface();
                  SmPoint2d         sUV2      = rTsectPnt.UVPos(1-lWhichSurf);
                  SmBrep          * pBrep     = pEdge->GetBrep() ;
                  pBaseSurf1->EvaluatePoint(sUV, sP1);
                  pBaseSurf2->EvaluatePoint(sUV2,sP2);
                  
                  smgfx_Erase();
                  smgfx_SetLook(1,2, 0,0,0) ; if(pBrep) { pBrep->Draw(TRUE) ; } sm_GraphicsLoop() ;
                  smgfx_SetLook(2,3, 0,1,1) ; pEdge->Draw() ; sm_GraphicsLoop() ;
                  smgfx_SetLook(4,6, 1,0,0) ; sPVSt[0].Draw() ; sm_GraphicsLoop() ;
                  smgfx_SetLook(4,6, 0,0,1) ; sPVEn[0].Draw() ; sm_GraphicsLoop() ;
                  smgfx_SetLook(1,6, 0,1,0) ; sP1.Draw(); sm_GraphicsLoop();
                  smgfx_SetLook(1,6, 0,0,0) ; sP2.Draw(); sm_GraphicsLoop();
                  if ( iDebugLevel >= 10 ) {
                      smgfx_SetLook( 4,6, 0,1,0) ; pEUFilletedEdge->Draw();  sm_GraphicsLoop();
                      smgfx_SetLook( 1,2, 1,1,0 ); pBaseSurf1->Draw(); sm_GraphicsLoop();

                      // Look at the curve that we dropped to, in 3d.
                      SmSurface * pSrf = pEUFilletedEdge->GetFaceuse()->GetFace()->GetSurface();
                      SmCrvOnSurf s3dCrv( *pUVCurveFilletedEdge, *pSrf );
                      smgfx_SetLook( 6,8, 1,0,0 ); s3dCrv.DrawParams(); sm_GraphicsLoop();
                  }
                  sm_GraphicsLoop();
              }
#endif // SM_DEBUG_CODE
              // sort by dT values - don't insert points already present
              ULONG lIndexOfInsert   = sProjParams.GetSize();
              SmBoolean bDoInsertion = TRUE;
              for (ULONG kk=0; kk<sProjParams.GetSize(); kk++)
                {
                  if (smos_Fabs(dT-sProjParams[kk]) < SM_EFF_ZERO_SQRT*dScale)
                    {
                      bDoInsertion = FALSE;
                      break;
                    }
                  else if (dT < sProjParams[kk])
                    {
                      lIndexOfInsert = kk;
                      break;
                    }
                } // end iter all sProjParams looking for duplicates and insert indices

              // when this point is to be inserted
              if (bDoInsertion)
                {
                  // Insert point - mark it by setting its curve parameter value to -SM_BIG_DOUBLE.
                  // Note, the value that was in there was not actually a curve parameter.
                  sProjParams.InsertAt(lIndexOfInsert, dT);
                  rTsectPnt.m_dCurveParameter = -SM_BIG_DOUBLE;
                  sTsectPnts.InsertAt(lIndexOfInsert, &rTsectPnt);
                }
            } // end iter all tsect points building ordered sTsectPnts array

          // quit when sTsectPnts is built
          if (sTsectPnts.GetSize() > 1)
            {
              break;
            }
        } // end sorting branch
      else // no sorting branch
        {
          // just add points to sTsectPnts
          sTsectPnts.Append(sAllTsectPnts);
        }
    } // end iter twice building the sTsectPnts array

#ifdef SM_DEBUG_CODE
  if ( sbDumpFVUs>0 )
    { pFilletGeom0->DumpFilletVUs(); }
#endif // SM_DEBUG_CODE


  // check state
  ULONG lTotalTsects = sTsectPnts.GetSize();
  if ( lTotalTsects == 0 )
    {
      SER(SM_ERR);
    }

  // Q: Should we keep the intersection point associated with the
  //    to-be-filleted edge endPoints even when they don't happen to
  //    be the first or last SmTsectPnt rather than just keeping the
  //    very last points?
  // A: It's not essential that we do, but it gives us a knot right at
  //    the edge end, which very often is where it gets trimmed to,
  //    which can be handy.

  // Now, determine if the last two Tsect points were very close
  // to each other. If so, remove the next to last one.
  if ( lTotalTsects > 2 )
    {
      SmPoint3d sPnt1, sPnt2;
      SER( pOffsetSrf1->EvaluatePoint( sTsectPnts[lTotalTsects-1]->UVPos(0), sPnt1 ) );
      SER( pOffsetSrf1->EvaluatePoint( sTsectPnts[lTotalTsects-2]->UVPos(0), sPnt2 ) );
      SmVector3d sV = sPnt1 - sPnt2;
      double dLen = sV.Length();

      if ( dLen < m_dThisApproxTol3d )
        {
          // remove the penultimate point.
          sTsectPnts[lTotalTsects-2] = sTsectPnts[lTotalTsects-1];
          sTsectPnts.RemoveLast();
          lTotalTsects--;
        }
    } // end enough points to remove duplicate endPoints check

  // Similarly, if the first two Tsect points are very close
  // to each other. remove the second one.
  if ( lTotalTsects > 2 )
    {
      SmPoint3d sPnt1, sPnt2;
      SER( pOffsetSrf1->EvaluatePoint( sTsectPnts[0]->UVPos(0), sPnt1 ) );
      SER( pOffsetSrf1->EvaluatePoint( sTsectPnts[1]->UVPos(0), sPnt2 ) );
      SmVector3d sV = sPnt1 - sPnt2;
      double dLen = sV.Length();

      if (dLen < m_dThisApproxTol3d)
        {
          // remove the 2nd point
          sTsectPnts.RemoveAt(1);
          lTotalTsects--;
        }
    } // end enough points to remove duplicate beginPoints check

#ifdef SM_DEBUG_CODE
  if ( iDebugLevel > 0 ) 
    {
      smos_sprintf( sBuff,
          _T("\n  In SmFilletSolver::CalcFilletGeom(): Finished sorting TsectPnts: %ld\n"), lTotalTsects );
      smos_WriteBuffer(sBuff);

      smgfx_SetLook( 3,5, 0,0,1 );
      for ( ii=0; ii<lTotalTsects; ii++ )
        {
          SmTsectPnt * pTsectPnt = sTsectPnts[ii];
          pTsectPnt->Dump();
          sm_DrawTsectPnt( pTsectPnt, m_pSurfaces[0], m_pSurfaces[1] ); sm_GraphicsLoop();
          sm_GraphicsLoop();
        }

      smos_sprintf( sBuff,_T("%s"),
          _T("\n  In SmFilletSolver::CalcFilletGeom(): Now creating the Fillet Geometry.\n") );
      smos_WriteBuffer(sBuff);
      smos_sprintf( sBuff,_T("%s"),
          _T("\n      First check whether Analytic is possible.\n") );
      smos_WriteBuffer(sBuff);
    }
#endif // SM_DEBUG_CODE


  // SolidCraft change: the round's exact class, from its faces' kinds and placement
  // (SmFilletExact.cpp): a round whose spine is rational is built exactly, every other traced.
  if ( lTotalTsects > 0 )
    {
      SmExactBlend sBlend;
      SER( ClassifyExactBlend( *sTsectPnts[0], sBlend ) );
      SmBoolean bBuilt = FALSE;
      SER( CalcExactFilletGeom( sBlend, sTsectPnts, NULL, FALSE, bBuilt ) );
      if ( bBuilt )
        {
#ifdef SM_DEBUG_CODE
          if ( iDebugLevel > 0 )
            {
              smos_sprintf( sBuff, _T("\n     SmFilletSolver::CalcFilletGeom(): Created exact fillet geometry, spine class %d.\n"),
                            static_cast<int>( sBlend.m_eSpine ) );
              smos_WriteBuffer( sBuff );
            }
#endif // SM_DEBUG_CODE
          return SM_SUCCESS;
        }
    }

#ifdef SM_DEBUG_CODE
  if ( iDebugLevel > 0 ) 
    {
      smos_sprintf( sBuff,_T("%s"),
          _T("\n  In SmFilletSolver::CalcFilletGeom(): Analytic not possible, doing Trace.\n") );
      smos_WriteBuffer(sBuff);
    }
#endif // SM_DEBUG_CODE

  // arrive here when we are unable to create an analytic fillet surface.
  // create a nurb fillet surface.

  // get fillet baseSurface domains
  SmExtent2d sDomain1 = pBaseSrf1->GetNaturalUVDomain();
  SmExtent2d sDomain2 = pBaseSrf2->GetNaturalUVDomain();

  // extend baseSurface domains by m_dSurfaceExtensionFactor * rDomain.GetSize()
  // don't extend in periodic directions nor through singular boundaries
  SetupOffsetExtension(pBaseSrf1, sDomain1);
  SetupOffsetExtension(pBaseSrf2, sDomain2);

  // Setup FilletIntersector : public SmAdvSurfaceIntersector : public SmSurfaceIntersector : public SmGlobalSolver
  SmFilletIntersector sFI(*pOffsetSrf1, sDomain1,
                          *pOffsetSrf2, sDomain2,
                          *this);

  // If we need to trace beyond the end and we are tracing along an arc,
  // Put a stop tracing criterion to avoid possible infinite looping
  double dRadius, dStartAng, dEndAng;
  SmAxis2Placement sRefFrame;

  // if extending the fillet, and pEdgeCurve is an arc:
  if (   (m_bExtendAfter || m_bExtendBefore)
      && pEdgeCurve->IsArc(5, m_dThisApproxTol3d, sRefFrame, dRadius, dStartAng, dEndAng))
    {
      // set max extension distance
      sFI.SetExtensionDistance(smos_Fabs(5.0*pOffsetSrf1->GetOffsetDistance()));
    }

  // locals
  SmBSplineCurve *p3DCurve = NULL ;
  SmTsectCurveType eCurveType;
  double dDeviation;
  double dThisApproxTol3d  = GetThisApproxTol3d();
  double dAngleTol         = GetThisAngTolRad();
  SmVector3d sDir          = sPVSt[1];                // intersection trace start direction

  // Check the last TsectPnt to see if it is degenerate (at a surf/surf tangency point),
  // If so, reverse the tracing to avoid its early termination
  SmBoolean bReverseTracing = FALSE;
  if (   lTotalTsects > 0
      && sTsectPnts[lTotalTsects-1]->m_ePointType == SM_IP_TANGENT_POINT)
    {
      sTsectPnts.ReverseArray( 0, lTotalTsects );
      sDir = -sPVEn[1];
      bReverseTracing = TRUE;
    }

  // SolidCraft change: a trace that starts at a tangent point starts toward the next through
  // point. The offset surfaces touch there: their normals' cross product gives no direction,
  // and the filleted edge's tangent at that end can be the line along which they touch (the
  // round/wall curve of a bar's top edge rounded into a cylinder tangent to the bar's side: it
  // ends along the bar's side, where a round of the bar's edge and the cylinder's offsets meet
  // to second order). A first step along that line lies on both offsets to within tolerance,
  // and the trace took it for the centre curve and ran off round the cylinder. The centre
  // curve leaves the point toward the next through point; the step from there finds it.
  if (   lTotalTsects > 1
      && sTsectPnts[0]->m_ePointType == SM_IP_TANGENT_POINT )
    {
      SmPoint3d sFrom, sTo;
      SER( pOffsetSrf1->EvaluatePoint( sTsectPnts[0]->UVPos(0), sFrom ));
      SER( pOffsetSrf1->EvaluatePoint( sTsectPnts[1]->UVPos(0), sTo ));
      SmVector3d sChord = sTo - sFrom;
      if ( sChord.Length() > dThisApproxTol3d )
        { sDir = sChord; }
    }

  // get the param value at start end of pEdgeCurve
  //double dParam = (bReverseTracing == TRUE)
  //                ? sIvl.GetMax()
  //                : sIvl.GetMin() ;

  // locals
  SmBoolean bExtendBefore   = m_bExtendBefore;
  SmBoolean bExtendAfter    = m_bExtendAfter;
  SmBoolean bFoundAllPoints = FALSE;
  ULONG lMaxLoopCount       = 10;
  ULONG lLoopCount          = 0;

  // Trace fillet centerlines until all FoundPoints are consumed.
  //   GWC: typically expect a single centerline to run from start to stop point in a single iteration
  //   but multiple iterations are possible when DoPointIntersection stepping is interrupted by
  //   running into a surface singularity or a surface cusp (made with a fillet radius
  //   larger than a surface curvature) where the step size or the surface tangent length
  //   gets so small that the stepping algorithm decides that stepping has stopped moving.
  //   When that happens this loop tries to recover by jumping the trace a bit and restarting
  //   to create multiple fillet curve sections which then get blended back together.
  while (!bFoundAllPoints)
    {
      lLoopCount++;
      if (lLoopCount > lMaxLoopCount)
        {
          // Either too many self-intersections or other
          // suspicious situations might have happened

          SetStatus( SM_FS_SURF_TRACING_FAILURE );
          RecordFilletError( SM_FILERR_BAD_INT, sTsectPnts,
              _T("Offset surface intersection failed, probably radius too large") );

          SER_MSG(SM_ERR, _T("CalcFilletGeom(): Ran into too many self-intersections or other problem"));
        }

#ifdef SM_DEBUG_CODE
      if ( iDebugLevel > 20 ) 
       {
         if(GetSurface(0))
           {
             smgfx_SetColor(.5,.8,0);
             GetSurface(0)->Draw() ; sm_GraphicsLoop();
             GetSurface(0)->DrawUV(4,4) ; sm_GraphicsLoop();
             if(GetSurface(0)->GetBaseSurface())
               {
                 GetSurface(0)->GetBaseSurface()->Draw() ; sm_GraphicsLoop();
                 GetSurface(0)->GetBaseSurface()->DrawUV(4,4) ; sm_GraphicsLoop();
               }
           }
         if(GetSurface(1))
           {
             smgfx_SetColor(.5,0,.8);
             GetSurface(1)->Draw() ; sm_GraphicsLoop();
             GetSurface(1)->DrawUV(4,4) ; sm_GraphicsLoop();
             if(GetSurface(1)->GetBaseSurface())
               {
                 GetSurface(1)->GetBaseSurface()->Draw() ; sm_GraphicsLoop();
                 GetSurface(1)->GetBaseSurface()->DrawUV(4,4) ; sm_GraphicsLoop();
               }
           }
         sm_GraphicsLoop();
       }
#endif // SM_DEBUG_CODE

      // Trace the rail curves and create fillet surface

      // Determine the tracing direction
      if (bReverseTracing == TRUE)
        {
          //sPV[1] = -sPV[1];
          bExtendBefore = m_bExtendAfter;
          bExtendAfter  = m_bExtendBefore;
        }

      // Do the offset-surface intersection.
      // Build the offset surface/surface intersection curve
      // and associated filletSurface and railCurves.
      // The base function SmSurfaceIntersector::DoPointIntersection
      // calls virtual SmFilletIntersector::FlushCurve() function
      // to cause the following side effects:
      //   1.  pFilletGeom = m_rFilletSolver.GetLastFilletGeom()
      //   2a. pFilletGeom->m_pFilletSurface    = newly created fillet surface
      //   2b. pFilletGeom->m_pCenterLineCurve  = newly created CenterLine
      //   2c. pFilletGeom->Rail[i]->Edgeuse->Mate->UVCurve =
      //         newly created rail UVCurves
      //   2d. pFilletGeom->Rail[i]->Curve      = newly created rail 3D Curves

      SmBSplineCurve *pSurfaceUVCurve1 = NULL;
      SmBSplineCurve *pSurfaceUVCurve2 = NULL;

#ifdef SM_DEBUG_CODE
  if ( sbDumpFVUs>0 )
    { pFilletGeom0->DumpFilletVUs(); }
#endif // SM_DEBUG_CODE

      SmStatus eStat = sFI.DoPointIntersection
          (m_crContext,        // in : context for new construction
           sTsectPnts,         // in : points on the intersection
           bExtendBefore,      // in : allow/prohibit extensions from
           bExtendAfter,       //      sTsectPnts 1st and last points
           &sDir,              // in : Opt start direction, 
           NULL,               // in : Opt end direction
           SM_CAST_APPROXTOL3D_PTR(&dThisApproxTol3d), // in : specify dist tol, NULL = use m_dThisApproxTol3d value
           &dAngleTol,         // in : specify ang tol,  NULL = use m_dThisAngTolRad    value
           p3DCurve,           // out: the fillet offsetSurface intersection curve
           pSurfaceUVCurve1,   // out: Resulting UVTrimCurve on Surface1           
           pSurfaceUVCurve2,   // out: Resulting UVTrimCurve on Surface2           
           eCurveType,         // out: specify the kind of intersection curve found
           dDeviation) ;       // out: max distance from through points to source surfaces


#ifdef SM_DEBUG_CODE
  // draw Breps(blue,green), Surfaces(cyan,yellow), faces(black), and newly created fillet geometry
  if (iDebugLevel > 0) 
    {
      SmFace *pFace1 = (SmFace *)sFI.GetSurface(0)->GetFace() ;
      SmFace *pFace2 = (SmFace *)sFI.GetSurface(1)->GetFace() ;
      SmBrep *pBrep1 = pFace1 ? pFace1->GetBrep() : NULL ;
      SmBrep *pBrep2 = pFace2 ? pFace2->GetBrep() : NULL ;

      smgfx_Erase();
      smgfx_SetLook(1,2, 0,0,1) ; if(pBrep1) { pBrep1->Draw(TRUE) ; } sm_GraphicsLoop() ;
      smgfx_SetLook(1,2, 0,1,0) ; if(pBrep2) { pBrep2->Draw(TRUE) ; } sm_GraphicsLoop() ;
      smgfx_SetLook(1,2, 0,1,1) ; if(sFI.GetSurface(0)) sFI.GetSurface(0)->DrawUV(4,4,FALSE,NULL,pBrep1==NULL); sm_GraphicsLoop();
      smgfx_SetLook(1,2, 1,1,0) ; if(sFI.GetSurface(1)) sFI.GetSurface(1)->DrawUV(4,4,FALSE,NULL,pBrep2==NULL); sm_GraphicsLoop();
      smgfx_SetLook(1,2, 0,0,0) ; if(pFace1) { pFace1->Draw(SM_DM_CROSSHATCH) ; } sm_GraphicsLoop() ;
      smgfx_SetLook(1,2, 0,0,0) ; if(pFace2) { pFace2->Draw(SM_DM_CROSSHATCH) ; } sm_GraphicsLoop() ;

      SmFilletGeom        *pCurrFilletGeom    =  sFI.GetCurrFilletGeom()
                                               ? sFI.GetCurrFilletGeom()
                                               : sFI.GetFilletSolver()->GetLastFilletGeom();
      // gwc: I think pSideSurface no longer has to be a BBSplineSurface - replace 1 line
      // rm : SmBSplineSurface * pFilletSurface = pCurrFilletGeom->GetFilletSurface() ;                                    
      SM_FILLETSURF_TYPE * pFilletSurface = pCurrFilletGeom->GetFilletSurface() ;                                    
      SmBSplineCurve   * pCenterLine    = pCurrFilletGeom->GetCenterLineCurve() ;
                                        
      SmFilletEdge     * pRailEdge0     = pCurrFilletGeom->GetRail(0) ;                                            
      SmBSplineCurve   * pRailUVCurve0  = pRailEdge0 ? pRailEdge0->GetPrimaryEdgeuse()->GetUVTrimCurve(): NULL ;
           
      SmFilletEdge     * pRailEdge1     = pCurrFilletGeom->GetRail(1) ;
      SmBSplineCurve   * pRailUVCurve1  = pRailEdge1 ? pRailEdge1->GetPrimaryEdgeuse()->GetUVTrimCurve() : NULL ;     
      
      smgfx_SetLook(2,3, 1,0,0) ; if(pCenterLine) pCenterLine->Draw(NULL, TRUE) ; sm_GraphicsLoop() ;
      smgfx_SetLook(2,3, 0,0,1) ; if(pRailEdge0) pRailEdge0->Draw() ; sm_GraphicsLoop() ;
      smgfx_SetLook(2,3, 0,1,0) ; if(pRailEdge1) pRailEdge1->Draw() ; sm_GraphicsLoop() ;
      smgfx_SetLook(1,2, 1,.5,0); if(pFilletSurface) pFilletSurface->DrawUV(8,8,FALSE,NULL,TRUE) ; sm_GraphicsLoop() ;
      smgfx_SetLook(3,4, 0,1,1) ; if(pRailUVCurve0) { SmCrvOnSurf sSurf1Curve( *pRailUVCurve0, (SmSurface &)*sFI.GetSurface(0)) ;
                                                      sSurf1Curve.Draw(NULL, TRUE) ; sm_GraphicsLoop() ;
                                                    }   
      smgfx_SetLook(3,4, 1,1,0) ; if(pRailUVCurve1) { SmCrvOnSurf sSurf2Curve( *pRailUVCurve1, (SmSurface &)*sFI.GetSurface(1)) ;      
                                                      sSurf2Curve.Draw(NULL, TRUE) ; sm_GraphicsLoop() ;
                                                    } 
      sm_GraphicsLoop();
    }
#endif // SM_DEBUG_CODE

      if ( eStat != SM_SUCCESS )
      {
          if ( m_eStatus == SM_FIL_UNPROCESSED )
          {
              SetStatus( SM_FS_SURF_TRACING_FAILURE );
          }
      }

      bFoundAllPoints = TRUE;

      // Check all points, after the first one, to make sure
      // that they were properly processed by DoPointIntersection().

      for (ii=1; ii<sTsectPnts.GetSize(); ii++)
        {
          SmTsectPnt * pTSP = sTsectPnts[ii];
          // any TsectPnt cruve param equal to -1 represents a through or end TsectPnt
          // that was missed in DoPointIntersection() - that's an error.  One expects
          // the xsect curve to run through all the TsectPnts given to DoPointIntersection().
          // Debugging note: check sm_IsPointOnCurve(), missing intersections too far away.
          if ( pTSP->m_dCurveParameter < 0.0 )
            {
#ifdef SM_DEBUG_CODE
              if ( iDebugLevel > 0 ) 
                {
                  smos_WriteBuffer(  _T("*** FILLET ERROR in SmFilletSolver::CalcFilletGeom(): ")) ;
                  smos_WriteBuffer(  _T("\n call to DoPointIntersection() stopped before stepping through all sTsectPnts.")) ; 
                  smos_WriteBuffer(  _T("\n   GWC: known to happen (sometimes but not always) when stepping into a ")) ;
                  smos_WriteBuffer(  _T("\n        centerline cusp made by a fillet radius larger than surface curvature.")) ;  
                }
#endif // SM_DEBUG_CODE
              // remember the failure to force another iteration and another piece of fillet curve to be stepped.
              bFoundAllPoints = FALSE;

              // on the next iteration trace from that startPoint in both directions.
              // this is part of the RanIntoACusp fix. We'll find a start point after
              // the current problem and trace that in both directions and we'll
              // try to hook the pieces back together.
              bExtendBefore   = TRUE;

              // remove the current 1st start point - we've tried walking from it.
              sTsectPnts.RemoveAt(0,ii);

              // done checking through points - on with the recovery
              break;
            }
        } // end iter sTsectPnts

      // quit with success when all FoundPoints have been processed
      if ( bFoundAllPoints )
        {
          break ; //Success
        }

      // when able to trace the the fillet to a surface natural boundary - marked in m_bCurrFilletTouchBoundary
      if ( sFI.IsCurrFilletTouchBoundary() )
        {
          // assume there were no self intersections: 
          // GWC: this is not absolutely true: my_testFillet_regression lCase = 2, lCount = 200 ;
          //      has two self intersections, one was navigated by the DoPointIntersections and
          //      one stopped the DoPointIntersection stepping sequence.  This bit of
          //      code assumes a self intersection always stops DoPointIntersection stepping.
          sFI.SetSelfIntersect( FALSE );

          // The newly created fillet touches the boundary edge of
          // a base surface.  Typically, it means that we have a G1-closed
          // base surface and the fillet needs to cross the boundary.
          // Trace again to create another fillet.
          if ( m_bExtendAfter ) 
            {
              // Something strange happening, return ERROR
              SER(SM_ERR);
            }

          // no need to trace next start point in both directions
          bExtendBefore = FALSE;
        }
      else // The newly created fillet did not trace to surface end - probably a surface cusp
        {
          // GWC: From my_testFillet_regression lCase = 2, lCount = 200
          //      we know that the DoPointIntersection stepping
          //      can be stopped by the cusp on a self-intersecting fillet caused
          //      by a fillet radius larger than the surface curvature.  The stepping
          //      speed decreases until a pair of points on the fillet are found
          //      so close together that the stepping stops commonly at the 
          //      first encounter of the cusp point where the tangent direction flips
          //      180 degrees. However, some cases are known to exist where the
          //      stepping manages to make it through a cusp.  So this is not
          //      a definitive test.
          // We have self-intersecting fillet. Will step forward from
          // where the solver stopped, insert a new TsectPnt and
          // re-trace the fillet. At the end, we will blend these two
          // fillets with a small patch.
          sFI.SetSelfIntersect( TRUE );

          SmFilletGeom   * pFG         = GetLastFilletGeom();
          SmBSplineCurve * pCenterLine = pFG->GetCenterLineCurve();

          // remember when none of fillet centerline curve has been stepped
          if ( pCenterLine == NULL )
            {
              SetStatus( SM_FS_SURF_TRACING_FAILURE );
              RecordFilletError( SM_FILERR_BAD_INT, sTsectPnts,
                                 _T("Offset surface intersection failed, probably radius too large") );
              NER( pCenterLine );
            }

          // get centerline's endPoint - it's located very near a surface cusp
          SmExtent1d sIvl1 = pCenterLine->GetNaturalInterval();
          SmPoint3d sPnt;
          SER(pCenterLine->EvaluatePoint(sIvl1.GetMax(),sPnt));

#ifdef SM_DEBUG_CODE
          if ( iDebugLevel > 0 ) 
            {
              sm_GraphicsLoop();
              smgfx_SetLook(2,1, 1,0,0); pCenterLine->Draw(); sm_GraphicsLoop();
              smgfx_SetLook(1,4, 0,0,1); sPnt.Draw(); sm_GraphicsLoop();
              sm_GraphicsLoop();
            }
#endif // SM_DEBUG_CODE

          // find closest filletEdge point to centerline end
          SmSolution sData[16];
          SmSolutionArray sSolutions(16,sData);
          SER(pEdgeCurve->GlobalPointSolve(sIvl, 
                                           SM_SO_MINIMIZE, 
                                           sPnt,
                                           dThisApproxTol3d, 
                                           NULL, NULL,
                                           SM_SR_SINGLE,sSolutions));
          
          // expect 1 solution. GWC: this might be before or after the cusp section!
          if (sSolutions.GetSize() != 1) SER(SM_ERR);

          // dT = param at which last centerline step ran into trouble (probable cusp)
          double dT = sSolutions[0].m_vStart[0];

          // pick initial skip size
          double dSkipSize = 0.1*sIvl.GetLength();
          if (bReverseTracing == TRUE)
            {
              dSkipSize *= -1.0;
            }

          // make an additional TSectPoint which may or may not be added to sTsectPnts
          // to handle self-intersection case - always add point to 
          // sDeletePnts array to make sure object is deleted when this function returns.
          SmTsectPnt *pNewTsectPnt = new SmTsectPnt() ;   // scope: while (!bFoundAllPoints) && SelfIntersecting Check
          sDeletePnts.Add(pNewTsectPnt);

          // for a sequence of bigger steps - try tracing the curve after the cusp region
          for (ii=0; ii<20; ii++)
            {
              // this iters startT = dt + scale * stepsize
              double dNextT = dT + dSkipSize*(1.0+0.1*ii);
              if (!sIvl.ContainsValue(dNextT)) 
                {
                  continue;
                }

              // Define a cross-sectional plane from surface offset.
              // This should produce better results for non-constant-radius filleting
              SmPoint3d  sPt;
              SmVector3d sBinVec;
              SmVector3d sTangent;
              SmVector2d sUVParameter;
              SmEdgeuse *pEUOnSurf = GetEdgeuse(0);

              // evaluate filletEdge->edgeuse to Surface0
              SER(pEUOnSurf->EvaluateBinormal(dNextT, 
                                              FALSE,           // in : TRUE  = when pUVTrimCurve==NULL, proj Edge pts to Surface else use pUVTrimCurve 
                                                               //      FALSE = get/create pUVTrimCurve to calc Surface points
                                              sPt,             // out: 3D pt on edge                                                 
                                              sBinVec,         // out: unit-vector pointing to SmFace interior from rBinormalPoint   
                                             &sTangent,        // out: non-unit Edgeuse tangent at rBinormalPoint                    
                                              NULL,            // out: unit Faceuse normal                                           
                                              &sUVParameter)); // out: Surface param value at edge point  
              
              // scale BinVec to dOffset length                                                          
              double dOffset = smos_Fabs(GetSurface(1)->GetOffsetDistance());
              sBinVec        = sBinVec*dOffset;

              // find Surface0 UVDirection moving to SmFace0 interior 
              SmVector2d sUVDir;
              SER(pBaseSrf1->DropVectors(sUVParameter, // in : Target UV Point                                                                   
                                         TRUE,         // in : if crUV is on U interval boundary, TRUE=Use upper interval, FALSE = lower interval
                                         TRUE,         // in : if crUV is on V interval boundary, TRUE=Use upper interval, FALSE = lower interval
                                         1,            // in : Number of vectors in following array                                              
                                        &sBinVec,      // in : Vectors to drop into parameter space.                                             
                                        &sUVDir));     // out: Corresponding parameter space vectors produced by drop  
              
              // get Surface0 guess for next trace start point                                                    
              sUVParameter   = sUVParameter + sUVDir;
              SER(GetSurface(0)->EvaluatePoint(sUVParameter,sPt));

              // Adjust tracing direction
              sDir = sTangent;

              // Solve the point on the cross-section plane
              // Will find and insert a new SmTsectPnt
              if(   sm_FindFilletTsectPnt(this,          // in : Fillet Solver for target to-be-filleted edge                             
                                          pBaseSrf1,     // in : Base Surface 1 connected to to-be-filleted edge                          
                                          pBaseSrf2,     // in : Base Surface 2 connected to to-bo-filleted edge                          
                                          sPt,           // in : guess point as close to desired point on fillet center curve as possible 
                                          sTangent,      // in : fillet center curve tangent at guess point                               
                                          *pNewTsectPnt) // out: contains fillet intersect UVpnt values                                   
                                                         //        baseSurface(uvPnt)   = points on rail                                  
                                                         //        OffsetSurface(uvPnt) = point on fillet center curve and                
                                                         //                               plane = [rPntToDrop, rTraceDir]                 
                 == SM_SUCCESS) 
                { // a FilletTSect point was found
                  
                  SmTArray<SmCurve*> sCurves;

                  // gather all fillet centerlines made so far
                  for (ULONG i=0; i<m_vFilletGeoms.GetSize(); i++)
                    {
                      SmFilletGeom *pFGM = m_vFilletGeoms[i];
                      SmBSplineCurve * pCenterLineCrv = pFGM->GetCenterLineCurve();
                      if (pCenterLineCrv != NULL)
                        {
                          sCurves.Add( pCenterLineCrv );
                        }
                    }

                  // When new point is not on an existing centerline
                  if (!sFI.IsPointOnCurve(sPt,sCurves))
                    {
                      // use it to start the next segment of the fillet trace
                      sTsectPnts.InsertAt(0,pNewTsectPnt);
                    }

                  // break from sequence of bigger steps until we are past the cusp
                  break;
                } // end sm_FindFilletTsectPnt for dNextT is a success check
            } // end iter 20 times trying larger dNextT steps
        } // end checking for a self-intersecting fillet - if so add a TSectPoint to bFoundAllPoints
    } // end while !bFoundAllPoints during fillet stepping algorithm

  // when just one filletGeom was created (no selfIntersections)
  if (m_vFilletGeoms.GetSize() == 1)
    {
      // See if we are required to split FG
      if (bSplitOriginalFillet == TRUE)
        {
          SmFilletGeom * pNewFilletGeom = NULL;
          SER(pFilletGeom0->TopologySplit(pNewFilletGeom,SM_FG_DEFAULT));
          this->m_vFilletGeoms.Add(pNewFilletGeom);
          double dSplitAtParam = 0.0;

          // Setup rail curves for new FilletGeom
          for (ULONG mm=0; mm<2; mm++)
            {
              SmFilletEdge    * pRail           = pFilletGeom0->GetRail(mm);
              SmCurve         * pRailCurve      = pRail->GetCurve(); NER(pRailCurve);
              SmExtent1d        sTrimIvl        = pRailCurve->GetNaturalInterval();
              SmExtent1d        sNewTrimIvl     = sTrimIvl;
              SmEdgeuse       * pMateEU         = pRail->GetPrimaryEdgeuse()->GetMate();
              SmCurve         * pUVCurve        = pMateEU->GetUVTrimCurvePointer(); NER(pUVCurve);
              SmFilletEdge    * pNewRail        = pNewFilletGeom->GetRail(mm);
              SmFilletEdgeuse * pNewMateEU      = (SmFilletEdgeuse*)pNewRail->GetPrimaryEdgeuse()->GetMate();
              SmCurve         * pNewRailCurve   = NULL;
              SmBSplineCurve  * pNewRailUVCurve = NULL;

              // Copy the original rail curve and UVTrimCurve
              SER(pRailCurve->Copy(m_crContext,pNewRailCurve));

              // Temp workaround to appease Linux gcc compiler
              // Original code with typecast has warning: dereferencing type-punned pointer will break strict-aliasing rules [-Wstrict-aliasing]
              // SER(pUVCurve->Copy(m_crContext,(SmCurve*&)pNewRailUVCurve));
              SmCurve* pCrv = NULL;
              SER(pUVCurve->Copy(m_crContext, pCrv));
              pNewRailUVCurve = (SmBSplineCurve*)pCrv;

              // Set NewRail->Curve = RailCurveCopy
              pNewRail->SetCurve(pNewRailCurve, TRUE ) ; // TRUE = delete preExisting Edge->Curve
                                                         // side effect: delete current pSurviveEdge->UVTrimCurve
              pNewRailCurve->SetOwner(pNewRail);

              // Set NewRail->UVTrimCurve = Rail->UVTrimCurveCopy
              pNewMateEU->SetUVCurve(pNewRailUVCurve, TRUE ) ; // TRUE = delete preExisting Edgeuse->UVTrimCurve - expect none here

#ifdef SM_USE_NEWTOL      
              SM_NEWTOL_LINE pNewRail->ClearLocalZoneTol3d( ) ; // what's going on with the MaxGaps here?
#else // SM_USE_OLDTOL 
              SM_OLDTOL_LINE pNewRail->SetTolerance(pRail->GetTolerance());
#endif // SM_USE_OLDTOL

              SmPoint2d sUV = pSplitTsectPnt->UVPos(mm);
              SmPoint3d sPnt;
              if (mm == 0)
                {
                  SER(pBaseSrf1->EvaluatePoint(sUV,sPnt));
                  // Find the split-at parameter(i.e. -u) of the fillet
                  double dDroppingTol = dThisApproxTol3d;
                  SmBoolean bSuccess = false;
                  double dDist;
                  for (double dScale=10.0; dScale<39.0; dScale+=10.0)
                    {
                      SER(pRailCurve->DropPoint(sTrimIvl,      // in : target curve allowed domain
                                                sPnt,          // in : Point to drop to curve
                                                NULL,          // in : For Snapping, NULL=NoSnap, Vec pointing to CrvInside from desired EndPt for snapping.
                                                               //      when DropPt is within dDistTol in cpVecToCrvInside dir of EndPt, Snap to Ivl EndPt.
                                                               //      Vec handles ambiguities and makes sure that curves are not just touching at the ends.
                                                dDroppingTol,  // in : Used for SnapDist, and SM_SO_INTERSECT sXSectTol3d distance.
                                                               //      sXSectTol3d = dDistTol ;      SM_SO_INTERSECT 
                                                               //      sXSectTol3d = SM_BIG_DOUBLE ; SM_SO_MINIMIZE, SM_SO_NORMALIZE
                                                               //      sSnapTol3d  = dDistTol ;      SM_SO_MINIMIZE, SM_SO_NORMALIZE, SM_SO_INTERSECT    
                                                NULL,          // in : If given, call cheap LocalPointSolve() else call expensive GlobalPointSolve()
                                                bSuccess,      // out: TRUE = found a drop point
                                                dSplitAtParam, // out: found drop curve param
                                                dDist)) ;      // out: found drop distance
                                                               // in : SM_SO_MINIMIZE = allow nonNormal drops near endPoints
                                                               //      SM_SO_INTERSECT= like Minimize but point must be within dDistTol
                                                               //      SM_SO_NORMALIZE= exclude nonNormal drops near endPoints
                                                               //      default:[SM_SO_MINIMIZE] to preserve original behavior
                      dDroppingTol = dThisApproxTol3d*dScale;
                      if (bSuccess) break;
                    }
                  if (!bSuccess) SER(SM_ERR);
                }
              else
                {
                  SER(pBaseSrf2->EvaluatePoint(sUV,sPnt));
                }
              SER( sTrimIvl.SetMinMax(sTrimIvl.GetMin(),dSplitAtParam) );
              SER(pRailCurve->Trim(sTrimIvl));  // may snap sIvl by tol to existing knots
              SER( sNewTrimIvl.SetMinMax(dSplitAtParam,sNewTrimIvl.GetMax()) );
              SER(pNewRailCurve->Trim(sNewTrimIvl));  // may snap sIvl by tol to existing knots
              SmFilletVertex * pNewVert = (SmFilletVertex*)pNewRail->GetStartVertex();
              pNewVert->SetPoint(sPnt);
              pNewVert->SetStatus(SM_FIL_PROCESSED);
              SmTArray<SmVertexuse*> sVUs;
              pNewVert->GetVertexuses(sVUs);
              for (ULONG jj=0; jj<sVUs.GetSize(); jj++)
                {
                  SmFilletVertexuse * pVU = (SmFilletVertexuse*)sVUs[jj];
                  SmFilletEdge * pE = (SmFilletEdge*)pVU->GetEdgeuse()->GetEdge();
                  if (pE->GetFilletEdgeType() != SM_FE_RAIL) continue;
                  pVU->SetTsectPnt(*pSplitTsectPnt);
                }
            }

          // Split fillet into two
          // gwc: I think pSideSurface no longer has to be a BBSplineSurface - replace 1 line
          // rm : SmBSplineSurface * pFilletSurface = pFilletGeom0->GetFilletSurface();
          SM_FILLETSURF_TYPE * pFilletSurface = pFilletGeom0->GetFilletSurface();
          NER(pFilletSurface);
          SmSurface * pSurL = NULL;
          SmSurface * pSurR = NULL;
          SER(pFilletSurface->SplitAt(m_crContext,dSplitAtParam,SM_SP_U,pSurL,pSurR));
          SM_ASSERT(pFilletSurface != NULL) ; delete pFilletSurface ; pFilletSurface = NULL ;
          
          pFilletGeom0->m_pFilletSurface   = SM_CAST_PTR(SmBSplineSurface,pSurL);
          pNewFilletGeom->m_pFilletSurface = SM_CAST_PTR(SmBSplineSurface,pSurR);

        } // end need to split fillet check
    } // end just one FilletGeom created check (no selfIntersections)

  // when more than one filletGeom was created (self intersections)
  if (m_vFilletGeoms.GetSize() > 1)
    {
      // Self-intersection cases
      for (ULONG j=0; j<m_vFilletGeoms.GetSize(); j++)
        {
          SmFilletGeom * pFilletGeom = m_vFilletGeoms[j];
          //if (pFilletGeom->GetFilletGeomType() == SM_FG_GAP_FILLER) {
          //    continue;
          //}
          if (pFilletGeom->GetFilletGeomType() == SM_FG_DEFAULT)
            {
              // Calculate fillet edges of type SM_FE_CROSS_SECTION
              SmTArray<SmFilletEdge*> sFGEdges;
              pFilletGeom->GetFilletEdges(sFGEdges);
              for (ULONG jj=0; jj<sFGEdges.GetSize(); jj++)
                {
                  SmFilletEdge * pE = sFGEdges[jj];
                  if (pE->GetFilletEdgeType() == SM_FE_CROSS_SECTION)
                    {
                      SER(pE->CalcCrossSection(m_crContext));
                    }
                }
              // Recompute FG if necessary - depends on the status flag
              //SER(pFilletGeom->ReCalcFilletGeom());
              //SER(pFilletGeom->TrimRailCurves()); // No trimming is neeedd here
            }
        }

      for (ULONG k=0; k<m_vFilletGeoms.GetSize(); k++)
        {
          SmFilletGeom * pFilletGeom = m_vFilletGeoms[k];
          if (pFilletGeom->GetFilletGeomType() == SM_FG_BLENDS)
            {
              // Calculate the surface for blending patch
              SER(pFilletGeom->CalcBlendingGeom());
            }
        }
    } // end more than one filletGeom created check

  // when trace was reversed (due to endPoint degeneracies) - restore output trace directions
  if (bReverseTracing)
    {
      // Reverse the fillet surface & rails
      for (ULONG j=0; j<m_vFilletGeoms.GetSize(); j++)
        {
          SmFilletGeom * pFilletGeom = m_vFilletGeoms[j];
          for (ULONG k=0; k<2; k++)
            {
              SmEdge         * pRail    = pFilletGeom->GetRail(k);
              SmBSplineCurve * p3DCrv = SM_CAST_PTR(SmBSplineCurve, pRail->GetCurve()); 
              NER(p3DCrv);
              SmEdgeuse      * pPrimEU  = pRail->GetPrimaryEdgeuse();
              SmEdgeuse      * pMateEU  = pPrimEU->GetMate();
              SmBSplineCurve * pUVCurve = NULL ;
              SER( pMateEU->GetOrCreateUVTrimCurve(pUVCurve) );
              NER(pUVCurve);
              SmExtent1d       sIvl1    = p3DCrv->GetNaturalInterval();
              SER(p3DCrv->ReverseParameterization(sIvl1,sIvl1));
              SER(pUVCurve->ReverseParameterization(sIvl1,sIvl1));
            }

          SmBSplineCurve   * pCenterline = pFilletGeom->GetCenterLineCurve();
          // gwc: I think pSideSurface no longer has to be a BBSplineSurface - replace 1 line
          // rm : SmBSplineSurface * pFilletSurf = pFilletGeom->GetFilletSurface();
          SM_FILLETSURF_TYPE * pFilletSurf = pFilletGeom->GetFilletSurface();
          NER(pCenterline);
          NER(pFilletSurf);

          SmExtent1d         sIvl2       = pCenterline->GetNaturalInterval();
          SER(pCenterline->ReverseParameterization(sIvl2,sIvl2));
          SER(pFilletSurf->Reverse(SM_SP_U));
        }
    } // end reverse trace check

#ifdef SM_DEBUG_CODE
  if ( sbDumpFVUs>0 )
    { pFilletGeom0->DumpFilletVUs(); }
#endif // SM_DEBUG_CODE


  // all done
  return SM_SUCCESS;

} // end SmFilletSolver::CalcFilletGeom

/***************************************************************
PURPOSE:  Evaluate Func                  = Func(X) and
                      Jacobian Matrix[i][j] = dFi/dxj (partial Func[i]/partial X[j])
            Given X values.

USAGE ---
  This function is used by a NewtonRaphson solver
  to find fillet edge points near initial guesses.

  Each of the different SmFilletSolver derived types evaluate a slightly
  different function set to find a point which is on the plane stored in
  this SmFilletSphereSolveENFO::m_vPlaneABC and ::m_dPlaneD and satisfies
  the geometric requirements of the derived SmFilletSolver type.

  Currently the derived SmFilletSolver types and their function sets are:


  SmFilletSolver derived type         crX      Functions: F ,G  = points on offset Surfaces
                                                          F2,G2 = points on baseSurfaces
                                                          nF,nG = baseSurface surface normals
                                                          R     = radius or offsetSurface offset distance
                                                          L     = given distance

    SmConstantDistanceFS     [u1 v1 u2 v2 R]   rF[0] = m_vPlaneABC.Dot(sF) + m_dPlaneD;
                                               rF[1] = (F2-G2+R*nF-R*nG).x
                                               rF[2] = (F2-G2+R*nF-R*nG).y
                                               rF[3] = (F2-G2+R*nF-R*nG).z
                                               rF[4] = (F2-G2)*(F-G) - L*L

    SmConstantRadiusFS       [u1 v1 u2 v2]     rF[0] = m_vPlaneABC.Dot(sF) + m_dPlaneD;
                                               rF[1] = (F-G).x
                                               rF[2] = (F-G).y
                                               rF[3] = (F-G).z

    SmEiffelTowerFS          [u1 v1 u2 v2]     rF[0] = m_vPlaneABC.Dot(sF) + m_dPlaneD;
                                               rF[1] = (F-G).(F-G) - diffR*diffR           - Distance between f & g = diffR
                                               rF[2] = (nFxnG).(F-G)                       - nF, nG, F, G are coplanar
                                                        nF     (F-G)
                                               rF[3] = ---- . ------- - cos(m_dTowerAngle) - Angle between nF & (F-G) = m_dTowerAngle
                                                       |nF|    |F-G|

    SmSecantConstantRadiusFS [u1 v1 u2 v2]     rF[0] = m_vPlaneABC.Dot(sF) + m_dPlaneD;            - V2 = F-G2
                                               rF[1] = V2.V2 - m_dFilletRadii[0]*m_dFilletRadii[1] - Distance between F and G2
                                               rF[2] = V2.N2 - m_dFilletRadii[1]*smos_Cosine(m_dSecantAngle)
                                               rF[3] = V2.N3                                       - V2 is perpendicular to N3(=N1xN2)

    SmVariableRadiusFS       [u1 v1 u2 v2 t]   rF[0] = m_vPlaneABC.Dot(sF) + m_dPlaneD;
                                               rF[1] = (F-G) . X     - X component of vector F-G
                                               rF[2] = (F-G) . Y     - Y component of vector F-G
                                               rF[3] = (F-G) . Z     - Z component of vector F-G
                                               rF[4] = C'(t) * (C(t) - F)

***************************************************************/

SmStatus SmFilletSphereSolveENFO::Evaluate
  (const SmTArray<double> & crX,           // in : x of Ax=F
   SmTArray<double>       & rF,            // out: F of Ax=F function values,
   SmMatrix               * pOptJacobian,  // out: Partial derivatives of the functions.
   SmBoolean              & rbFoundAnswer) // out: TRUE = all rF members are less than scaled tolerance
                                           //      FALSE= not converged
{
  // check state
  SM_ASSERT(   (pOptJacobian == NULL)
            || (rF.GetSize() == pOptJacobian->GetNumRows())) ;

  // init output: let rF = 0.0 and pOptJacobian = 0.0
  rbFoundAnswer = FALSE;
  if (pOptJacobian) { for (ULONG ll=0; ll<pOptJacobian->GetNumRows(); ll++)
                        { for (ULONG mm=0; mm<pOptJacobian->GetNumColumns(); mm++)
                            { (*pOptJacobian)[ll][mm] = 0.0;
                            }
                          rF[ll] = 0.0;
                        }
                    }
  else              { for (ULONG ll=0; ll<rF.GetSize(); ll++)
                        { rF[ll] = 0.0;
                        }
                    }

  // LoadJacobian() locals
  ULONG       lNumEquations  = 1;                               // number of equations to skip LoadJacobian() sets rF and pOptJacobian
  ULONG       lNumParameters = 2;                               // index of surf2 UVPnt in crX array
  SmSurface * pSurface       = m_crFilletSolver.GetSurface(0);  // 1st offset surface
  ULONG       lSurf1Offset   = 0;                               // index of surf1 UVPnt in crX array
  SmSurface * pSurface2      = NULL;                            // 2nd offset surface, NULL == index of surf2 UVPnt in crX = lSurf1Offset + 2
  ULONG       lSurf2Offset   = 0;                               // not used when pSurface2 == NULL
  SmBoolean   bFoundAnswer;
  ULONG       lRailIndex     = 0;                               // rail index for pSurface

  // set rows from rF and pOptJacobian
  //  [lNumEquations] to [lNumEquations+Derived SmFilletSolver type number of equations]
  //   SmConstantDistanceFS     = 4 rows
  //   SmConstantRadiusFS       = 3 rows
  //   SmEiffelTowerFS          = 3 rows
  //   SmSecantConstantRadiusFS = 3 rows
  //   SmVariableRadiusFS       = 4 rows
  if (m_crFilletSolver.LoadJacobian(lNumEquations, lNumParameters,
                                    lRailIndex,
                                    pSurface,  lSurf1Offset,
                                    pSurface2, lSurf2Offset,
                                    crX,
                                    rF, pOptJacobian, bFoundAnswer) != SM_SUCCESS)
    {
      return SM_ERR;
    }

  // previous call sets all but the 1st equation in rF and pOptJacobian

  // Set first row of rF and pOptJacobian

  // The first function measures the distance between the Surf1 point and the target plane
  //
  // f0 = PlaneNormal*F + m_dPlaneD;
  //    = N*F + D
  //    = A*x + B*y + C*z + D
  //
  //    where F           = [x y z]
  //                      = OffsetSurf1.Position(crX[0], crX[1])
  //          D           = m_dPlaneD
  //          PlaneNormal = [A B C]
  //
  // The Jacobian first row nonZero entries are the first two entries as:
  //
  // pOptJacobian[0][] = [ A*du.x + B*du.y + C*du.z    A*dv.x + B*dv.y + C*dv.z  0 0 . . . ]


  // calc OffsetSurface0 current solution surface values
  SmPoint2d sUV(crX[0],crX[1]);
  SmPoint3d sF;
  SmVector3d sDUF, sDVF;
  SER(m_crFilletSolver.GetSurface(0)->Evaluate1stDerivatives(sUV,TRUE,TRUE,sF,sDUF,sDVF));

#ifdef SM_DEBUG_CODE
int iDebugLevel = m_crFilletSolver.DebugLevel();
  if ( iDebugLevel > 60 ) {
      sm_GraphicsLoop();
      smgfx_SetColor(1,0,0);
      m_crFilletSolver.GetSurface(0)->DrawAt(sUV,1);
      smgfx_SetColor(0,0,1);
      SmPoint2d sUV2(crX[lSurf2Offset],crX[lSurf2Offset+1]);
      m_crFilletSolver.GetSurface(1)->DrawAt(sUV2,1);
      sm_GraphicsLoop();
  }
#endif

  // update the plane to fit the current solution
  if (m_bPlaneAdjusts)
    {
      // get OffsetSurface0 surface normal at current solution
      SmVector3d sNorm1 = sDUF * sDVF;
      if (sNorm1.LengthSquared() < SM_EFF_ZERO_SQ)
        {
          SER(m_crFilletSolver.GetSurface(0)->EvaluateNormal(sUV,TRUE,TRUE,sNorm1));
        }
      SER(sNorm1.Unitize());

      // calc OffsetSurface1 current solution surface values
      SmPoint2d sUV2(crX[lSurf2Offset],crX[lSurf2Offset+1]);
      SmPoint3d sG;
      SmVector3d sDUG, sDVG;
      SER(m_crFilletSolver.GetSurface(1)->Evaluate1stDerivatives(sUV2,TRUE,TRUE,sG,sDUG,sDVG));

      // get OffsetSurface1 surface normal at current solution
      SmVector3d sNorm2 = sDUG * sDVG;
      if (sNorm2.LengthSquared() < SM_EFF_ZERO_SQ)
        {
          SER(m_crFilletSolver.GetSurface(1)->EvaluateNormal(sUV2,TRUE,TRUE,sNorm2));
        }
      SER(sNorm2.Unitize());

      // define new plane
      //   planeOrigin = avg(OffsetSurface0.position, OffsetSurface1.position)
      //   planeNormal = crossProduct(OffsetSurface0.SurfaceNormal, OffsetSurface1.SurfaceNormal)
      SmPoint3d  sPlaneOrig = (sF + sG) / 2.0;
      SmVector3d sPlaneABC = sNorm1 * sNorm2;

      // preserve plane normal orientation
      if (sPlaneABC.Dot(m_vPlaneABC) < 0.0) { sPlaneABC = - sPlaneABC; }

      // save plane in Ax + By + Cz - D = 0 form
      m_vPlaneABC = sPlaneABC;
      m_dPlaneD   = - m_vPlaneABC.Dot(sPlaneOrig);

    } // end m_bPlaneAdjusts == TRUE check

  // set rF 1st row
  rF[0] = m_vPlaneABC.Dot(sF) + m_dPlaneD;

  // set pOptJacobian 1st row nonZero entries
  if (pOptJacobian)
    {
      (*pOptJacobian)[0][0] = m_vPlaneABC.Dot(sDUF);
      (*pOptJacobian)[0][1] = m_vPlaneABC.Dot(sDVF);
    }

  // bFoundAnswer is currently based on all but the 1st row of rF values
  // when bFoundAnswer is TRUE
  if (bFoundAnswer)
    {
      // check the 1swt row rF value as well
      double dScaledTol = m_crFilletSolver.GetConversionTol() * (1.0 + sF.GetMaxDimension());
      if (smos_Fabs(rF[0]) < dScaledTol)
        {
          rbFoundAnswer = TRUE;
        }
    }

#ifdef SM_DEBUG_CODE
  if ( iDebugLevel > 60 ) {
      SmPoint3d sPnt;
      SER(pSurface->EvaluatePoint(sUV,sPnt));
      smgfx_SetColor(1,0,0);
      sPnt.Draw();
      sm_GraphicsLoop();
      if (pOptJacobian) {
          pOptJacobian->Dump();
          SmMatrix sTmpJ(6,6);
          SmTArray<double> sTmpX(crX);
          for (ULONG ll=0; ll<sTmpJ.GetNumRows(); ll++) {
              for (ULONG mm=0; mm<sTmpJ.GetNumColumns(); mm++) {
                  sTmpJ[ll][mm] = 0.0;
              }
          }
          SER(ComputeJacobian(sTmpX,rF,sTmpJ));
          sTmpJ.Dump();
          smos_WriteBuffer(_T(" F --- \n"));
          rF.Dump();
      }
  }
#endif

  return SM_SUCCESS;

} // end SmFilletSphereSolveENFO::Evaluate

/*******************************************************************//**
PURPOSE: Find fillet point on a given plane.

NOTES: Finds and stores surface uvPnt values such that
       BaseSurface(uvPnt)   = points on rails
       OffsetSurface(uvPnt) = point on fillet center curve (offsetSurface intersection curve)
                              and on specified plane

METHOD ---
  Given a point on and a plane perpendicular to the
     fillet offset-surface intersection-curve,
  Find the two intersection points between that plane and the
     fillet rail curves.

  The intersection points are found satisfying the particular equation
  set implemented within the virtual SmFilletSolver::LoadJacobian() methods
  for each SmFilletSolver derived type.

***********************************************************************/
SmStatus SmFilletSolver::PointOnPlaneSolve
  (const SmPoint3d  & crPlaneOrig,   // in : pt on offset-surface xsect curve
   const SmVector3d & crPlaneNormal, // in : tang at that offset-surface xsect curve
   const SmExtent2d & crUVDomain1,   // in : fillet-surf1 UVdomain
   const SmExtent2d & crUVDomain2,   // in : fillet-surf2 UVdomain
   const SmVector2d & rUV1,          // in : fillet-surf1 UV guess value
   const SmVector2d & rUV2,          // in : fillet-surf2 UV guess value
   SmBoolean  & rbFoundSolution,     // out: TRUE = next rTsectPnt found with
                                     //             SmLocalSolveNd:SmFilletSphereSolveENFO::SolveIt().
                                     //      FALSE= no next point found.
   SmTsectPnt & rTsectPnt)           // out: gets UV values for both surfaces
                                     //        BaseSurface(uvPnt)   = points on rails
                                     //        OffsetSurface(uvPnt) = point on fillet center curve
{
  // Now Set up Solver
  ULONG lJacobianSize = GetJacobianSize();

  SmExtentNd sIntervals(lJacobianSize);

  sIntervals[0] = SmExtent1d(crUVDomain1.GetMin().x,crUVDomain1.GetMax().x);
  sIntervals[1] = SmExtent1d(crUVDomain1.GetMin().y,crUVDomain1.GetMax().y);
  sIntervals[2] = SmExtent1d(crUVDomain2.GetMin().x,crUVDomain2.GetMax().x);
  sIntervals[3] = SmExtent1d(crUVDomain2.GetMin().y,crUVDomain2.GetMax().y);

  SmBoolean alPerData[16];
  SmTArray<SmBoolean> sPeriodicities(16,alPerData);

  sPeriodicities.Add(FALSE);
  sPeriodicities.Add(FALSE);
  sPeriodicities.Add(FALSE);
  sPeriodicities.Add(FALSE);

  double adGuessData[16];
  SmTArray<double> sGuessT(16,adGuessData);
  sGuessT.Add(rUV1.x);
  sGuessT.Add(rUV1.y);
  sGuessT.Add(rUV2.x);
  sGuessT.Add(rUV2.y);

  SmSurface *pSurface1      = GetSurface(0);  // offset of to-be-filleted edge surf1
  SmSurface *pSurface2      = GetSurface(1);  // offset of to-be-filleted edge surf2
  ULONG      lDoSurf2Calcs  = 0;               // Don't do surface2 calculations
  ULONG      lSurf1Index    = 0;               // UV values of Surface1
  ULONG      lSurf2Index    = 2;               // UV values of Surface2
  ULONG      lRailIndex     = 0;               // rail index for pSurface1

#ifdef SM_DEBUG_CODE
  if ( DebugLevel() >= 20 ) {
      smgfx_Erase();
      smgfx_SetColor(1,0,0);
      crPlaneOrig.Draw();
      crPlaneNormal.Draw(&crPlaneOrig);

      SmOffsetSurface *pOffSrf1 = SM_CAST_PTR( SmOffsetSurface, pSurface1 );
      SmOffsetSurface *pOffSrf2 = SM_CAST_PTR( SmOffsetSurface, pSurface2 );
      smgfx_SetColor(0,0,1);
      SmBoolean bBaseOnly = TRUE;
      if ( bBaseOnly && pOffSrf1 && pOffSrf2 )
      {
          // Draw only the base surfaces.
          pOffSrf1->GetBaseSurface()->DrawUV(2,2);
          pOffSrf2->GetBaseSurface()->DrawUV(2,2);
      }
      else
      {
          pSurface1->DrawUV(2,2);
          pSurface2->DrawUV(2,2);
      }
      sm_GraphicsLoop();
  }
#endif

  // Note that in the future we could load the additional initial values
  // from the current point and also setup all of the intervals, periodicities,
  // and guess values outside of here and pass them in.  This may help
  // performance a little because we could eliminate LoadInitialValues
  // call each time we are here.

  // However for now we will just load the initial values again.
  // given current inputs - place a good guess for upcoming solve in sGuessT
  SER(LoadInitialValues(lDoSurf2Calcs, lRailIndex,
                        *pSurface1,    lSurf1Index,
                         pSurface2,    lSurf2Index,
                         sIntervals,   sPeriodicities, sGuessT));

  double           sdData[16];
  SmTArray<double> sSolutionVector(16,sdData);

  // get plane equation for given PlanePt and PlaneNormal
  double dPlaneD = - crPlaneOrig.Dot(crPlaneNormal);

  // allocate NewtonRaphson evaluator - contains equation definitions being solved
  SmFilletSphereSolveENFO sEvalFun(crPlaneNormal,dPlaneD,*this);

  // allocate generic NewtonRaphson solver for this evaluator
  SmLocalSolveNd sLS(sEvalFun,          // in : Define Eqns to set to Zero (defines the DOF COUNT)
                     &sIntervals,       // in : NULL or sized:[N], N = this rEvalFun DOF count
                                        //      bounds on problem parameters.
                                        //      Set m_cpIntervals[i].SetUnbounded() for ivl whose prob params don't get clamped
                     &sPeriodicities) ; // in : NULL or sized:[N], N = this rEvalFun DOF count
                                        //      ivl[i] == TRUE = Interval[i] is periodic and values get wrapped
                                        //      ivl[i] == FALSE= Interval[i] is NOT periodic and values get clamped

  // set solve parameters
  sLS.SetDesiredAccuracy(SM_EFF_ZERO_SQ);
  sLS.SetBoundaryHandler(SM_BH_TOTAL_BOUNDARY_HITS,8);
  sLS.SetMaximumIterations(30);

  // Solve equations with NewtonRaphson iteration
  SER(sLS.SolveIt(sGuessT,m_dThisApproxTol3d,rbFoundSolution,sSolutionVector));

  // Try it again and allow the plane to adjust as convergance happens
  if (!rbFoundSolution) {
      sEvalFun.SetPlaneAdjusts(TRUE);
      SER(sLS.SolveIt(sGuessT,m_dThisApproxTol3d,rbFoundSolution,sSolutionVector));
  }

  // when a solution point was found
  if (rbFoundSolution)
    {
      // save the surface uvPnts
      rTsectPnt.UVPos(0).x = sSolutionVector[lSurf1Index];
      rTsectPnt.UVPos(0).y = sSolutionVector[lSurf1Index+1];
      rTsectPnt.UVPos(1).x = sSolutionVector[lSurf2Index];
      rTsectPnt.UVPos(1).y = sSolutionVector[lSurf2Index+1];
      for (ULONG i=4; i<sSolutionVector.GetSize(); i++)
        {
          rTsectPnt.m_adUserDoubles[i-4] = sSolutionVector[i];
        }
  }

  return SM_SUCCESS;

} // end SmFilletSolver::PointOnPlaneSolve


/****************************************************************
PURPOSE: Given an X vector in the equation set Ax=f,
  compute f.  Also compute the jacobian matrix of
  1st derivatives of f with respect to the x parameters.

  This function is called repeated by the multivariable
  NewtonRaphson solver.

NOTES:

  // The unknowns are:
  // x[0] = parameter for sideEdge->UVTrimCurve
  // x[1] = u of UVPnt on FilletFace->UVPnt
  // x[2] = v of UVPnt on FilletFace->UVPnt
  //
  // The first two functions are the difference between the uv of the side edge
  // at the current point and the uv of the surface that it's in:
  //   f0 = C(t).x - x[1]
  //   f1 = C(t).y - x[2]
  // where C(t) is the uv curve of the side edge, and t is x[0].
  //
  // The first three columns of the Jacobian will contain
  // dt  du  dv as follows:
  //
  // C(t)'.x      - 1      0.0 ...
  // C(t)'.y       0.0     - 1 ...
  //
 // The remaining function values and Jacobian will be filled in
  // by m_crFilletSolver.LoadJacobian

****************************************************************/
SmStatus SmRailUVCurveIntersectENFO::Evaluate
  (const SmTArray<double> & crX,           // in : x of Ax=F
   SmTArray<double>       & rF,            // out: F of Ax=F function values,
   SmMatrix               * pOptJacobian,  // out: Partial derivatives of the functions.
   SmBoolean              & rbFoundAnswer) // out: TRUE = all F values < m_crFilletSolver.GetConversionTol()
                                           //      FALSE= some large F values remain
{
  // init output
  rbFoundAnswer = FALSE;

  // get UVPnt = UVTrimCurve(crX[0])
  SmVector3d sPV[2];
  SER(m_crUVCurve.Evaluate(crX[0],1,TRUE,sPV));

#ifdef SM_DEBUG_CODE
int iDebugLevel = m_crFilletSolver.DebugLevel();
  if ( iDebugLevel > 40 ) 
    {
      SmCrvOnSurf s3DUVCurve((SmCurve&)m_crUVCurve,(SmSurface&)m_crSurfaceOfCurve);
      smgfx_SetColor(1,0,0);
      s3DUVCurve.Draw();
      sm_GraphicsLoop();
    }
#endif

  // evaluate 1st 2 equations
  rF[0] = sPV[0].x - crX[1];
  rF[1] = sPV[0].y - crX[2];

  // load 1st two rows in Jacobian - init rest of matrix to 0.0
  if (pOptJacobian)
    {
      (*pOptJacobian)[0][0] = sPV[1].x;
      (*pOptJacobian)[0][1] = - 1.0;
      (*pOptJacobian)[0][2] =   0.0;

      (*pOptJacobian)[1][0] = sPV[1].y;
      (*pOptJacobian)[1][1] =   0.0;
      (*pOptJacobian)[1][2] = - 1.0;

      // init rest of matrix to zero
      for (ULONG kk=3; kk<pOptJacobian->GetNumColumns(); kk++)
        {
          (*pOptJacobian)[0][kk] = 0.0;
          (*pOptJacobian)[1][kk] = 0.0;
        }
      for (ULONG ll=2; ll<pOptJacobian->GetNumRows(); ll++)
        {
          for (ULONG mm=0; mm<pOptJacobian->GetNumColumns(); mm++)
            {
              (*pOptJacobian)[ll][mm] = 0.0;
            }
          rF[ll] = 0.0;
        }

    } // end pOptJacobian == TRUE check

  ULONG       lNumEquations  = 2;
  ULONG       lNumParameters = 3;
  SmSurface * pSurface       = m_crFilletSolver.GetSurface(m_lRailIndex);
  ULONG       lSurf1Offset   = 1;
  SmSurface * pSurface2      = NULL;
  ULONG       lSurf2Offset = 0;
  SmBoolean   bFoundAnswer;

  // pass the rest of the evaluation along
  //  compute last set of rF values
  //  compute 2nd half of jacopbian matrix
  if (m_crFilletSolver.LoadJacobian(lNumEquations, lNumParameters,
      m_lRailIndex, pSurface, lSurf1Offset, pSurface2, lSurf2Offset,
      crX, rF, pOptJacobian, bFoundAnswer) != SM_SUCCESS)
    {
      return SM_ERR;
    }

  // 2nd half of equations set answer looks good
  // now check 1st 2 eqns to see if we set rbFoundAnswer = TRUE or FALSE
  if (bFoundAnswer)
    {
      double dScaledTol = m_crFilletSolver.GetConversionTol() * (1.0 + sPV[0].GetMaxDimension());
      if (   smos_Fabs(rF[0]) < dScaledTol
          && smos_Fabs(rF[1]) < dScaledTol)
        {
          rbFoundAnswer = TRUE;
        }
    }

#ifdef SM_DEBUG_CODE
  if ( iDebugLevel > 60 ) 
    {
      SmOffsetSurface *pOffSrf1 = m_crFilletSolver.GetSurface(   m_lRailIndex );
      SmOffsetSurface *pOffSrf2 = m_crFilletSolver.GetSurface( 1-m_lRailIndex );
      const SmSurface *pBaseSrf1 = pOffSrf1->GetBaseSurface();
      const SmSurface *pBaseSrf2 = pOffSrf2->GetBaseSurface();
      SmEdge *pE = m_crFilletSolver.GetEdgeuse(0)->GetEdge();

      if ( FALSE ) {
          smgfx_Erase();

          SmObject *pOwnerObj = pBaseSrf1->GetFace();
          SmFace *pFace = SM_CAST_PTR( SmFace, pOwnerObj );
          SmBrep *pBrep = ( pFace == NULL ) ? NULL : pFace->GetBrep();

          if ( pBrep ) {
              smgfx_SetLook( 1,2, 0,0,0 ); pBrep->Draw(TRUE); sm_GraphicsLoop();
              sm_GraphicsLoop();
          }

          smgfx_SetLook( 3,5, 1,0,1 ); if ( pE ) { pE->Draw(); } sm_GraphicsLoop();
          sm_GraphicsLoop();
      }

      SmPoint2d sUV( crX[1], crX[2] );
      SmPoint3d sBasePnt, sOffPnt;
      SER(pBaseSrf1->EvaluatePoint( sUV, sBasePnt ));
      SER(pOffSrf1 ->EvaluatePoint( sUV, sOffPnt  ));

      sm_GraphicsLoop();
      smgfx_SetLook(1,2, 0,0,1); sBasePnt.DrawPointToPoint( sOffPnt ); sm_GraphicsLoop();
      sm_GraphicsLoop();

      sUV.Set( crX[3], crX[4] );
      SER(pBaseSrf2->EvaluatePoint( sUV, sBasePnt ));
      SER(pOffSrf2 ->EvaluatePoint( sUV, sOffPnt  ));

      smgfx_SetLook(1,2, 0,1,1); sBasePnt.DrawPointToPoint( sOffPnt ); sm_GraphicsLoop();
      sm_GraphicsLoop();

      if ( pE ) {
          SmCurve *pCrv = pE->GetCurve();
          smgfx_SetLook( 2,3, 1,0,0 ); pCrv->DrawAt( crX[5], 1 ); sm_GraphicsLoop();
          sm_GraphicsLoop();
      }

      if (pOptJacobian)
        {
          pOptJacobian->Dump();
          SmMatrix sTmpJ(6,6);
          SmTArray<double> sTmpX(crX);
          SER(ComputeJacobian(sTmpX,rF,sTmpJ));
          sTmpJ.Dump();
          smos_WriteBuffer(_T(" X --- \n"));
          sTmpX.Dump();
          smos_WriteBuffer(_T(" F --- \n"));
          rF.Dump();
        }
    }
#endif

  return SM_SUCCESS;

} // end SmRailUVCurveIntersectENFO::Evaluate


/*******************************************************************//**
PURPOSE: Compute the fillet point corresponding to the intersection
     of a rail and an edge given a starting point on the edge.

NOTES:  This method assumes that pEdgeuse connects to the surface
     corresponding to the lRailIndex.
     rTsectPnt contains the edge parameter and other information about
     the intersection of the edge and the rail.

     This method sets rTsectPnt.m_dCurveParameter to rdEdgeuserParameter,
     which is nonstandard usage for that field.
***********************************************************************/
SmStatus SmFilletSolver::RailEdgeuseIntersect
 (ULONG        lRailIndex,           // in : target rail index
  SmEdgeuse  * pSideEU,              // in : target edgeuse->edge->curve to intersect
  double       dGuessEdgeParameter,  // in : Edgeuse->curve guess parameter
  SmBoolean  & rbFoundIntersection,  // out: TRUE=found an intersection
  SmTsectPnt & rTsectPnt,            // out: Contains solution
  double     & rdEdgeuseParameter)   // out: edgeuse Parameter of intersection
{
  // init output
  rbFoundIntersection = FALSE;

  // check state - need both offset surfaces
  if (GetSurface(0) == NULL || GetSurface(1) == NULL) { return SM_ERR; }

  // get original Brep sideEdge->Curve to intersect
  SmEdge  *pSideEdge = pSideEU->GetEdge();
  

  // let edgeuse->Surface UV guess = edgeuse->UVTrimCurve(GuessEdgeParameter)
  SmBSplineCurve *pUVCurve = NULL ;
  SER(pSideEU->GetOrCreateUVTrimCurve(pUVCurve));
  SmPoint3d s2DPnt;
  SER(pUVCurve->EvaluatePoint(dGuessEdgeParameter,s2DPnt));
  SmPoint2d sUVGuess(s2DPnt.x,s2DPnt.y);

#ifdef SM_DEBUG_CODE
int iDebugLevel = DebugLevel();
static ULONG lCount      = 1 ; lCount++ ;
static ULONG lDebugCount = 0 ;

SmCurve* p3DCurve = pSideEdge->GetCurve();

  // draw edgeuse->Curve->GuessPoint(red), edgeuse->Curve(yellow), edgeuse->Face(blue)
  //      edgeuse->OffsetSurface->GuessPoint(magenta), offsetSurfaces(green)
  if ( iDebugLevel > 0 || lCount == lDebugCount )
    {
      if ( FALSE )
        { smgfx_Erase(); }

      sm_GraphicsLoop();
      smgfx_SetLook(1,2, 1,0,0); p3DCurve->DrawAt(dGuessEdgeParameter,1); sm_GraphicsLoop();
      smgfx_SetLook(1,3, 1,1,0); p3DCurve->DrawWDeriv(pSideEdge->GetInterval(),0); sm_GraphicsLoop();
      smgfx_SetLook(1,2, 0,0,1); pSideEU->GetFace()->Draw();                 sm_GraphicsLoop();
      smgfx_SetLook(1,2, 0,1,1); GetSurface(lRailIndex)->DrawAt(sUVGuess,2); sm_GraphicsLoop();
      sm_GraphicsLoop();

      if ( iDebugLevel > 10 ) {
          SmOffsetSurface *pOffSrf1 = GetSurface(   lRailIndex );
          SmOffsetSurface *pOffSrf2 = GetSurface( 1-lRailIndex );

          // draw base or offset surfaces
          SmBoolean bBaseOnly = TRUE;
          if ( bBaseOnly )
            {
              smgfx_SetLook( 1,2, 0,0,0 ); pOffSrf1->GetBaseSurface()->DrawUV(4,4); sm_GraphicsLoop();
              smgfx_SetLook( 1,2, 0,1,0 ); pOffSrf2->GetBaseSurface()->DrawUV(7,7); sm_GraphicsLoop();
            }
          else
            {
              smgfx_SetLook( 1,2, 0,0,0 ); pOffSrf1->DrawUV(4,4); sm_GraphicsLoop();
              smgfx_SetLook( 1,2, 0,1,0 ); pOffSrf2->DrawUV(7,7); sm_GraphicsLoop();
            }
        }
      sm_GraphicsLoop();
    }

static SmBoolean sbDumpFVUs=FALSE;
  if ( sbDumpFVUs ) {
      SmTArray< SmFilletGeom* > sFGs;
      this->GetFilletGeoms( sFGs );
      for ( ULONG ij=0; ij<sFGs.GetSize(); ij++ )
      {
          SmFilletGeom *pFG = sFGs[ij];
          pFG->DumpFilletVUs();
      }
  }
#endif // SM_DEBUG_CODE

  // Now Set up Solver
  ULONG lJacobianSize = GetJacobianSize() + 1;

  SmExtentNd sIntervals(lJacobianSize);
  sIntervals[0] = pUVCurve->GetNaturalInterval();

  const SmSurface * pBSP    = GetSurface(lRailIndex)->GetBaseSurface();
  SmExtent2d        sDomain = pBSP->GetNaturalUVDomain();
  SetupOffsetExtension(pBSP,sDomain);

  // SolidCraft change: a straight side edge on a plane continues exactly along its line, as far
  // as the base surface's domain is extended across the edge (its extension factor, or the
  // trace's reach, SetTraceReach). A rail that leaves its face across a cliff meets the side
  // edge's line past the edge's end; bounded to the edge's own interval, the solve found that
  // point only from a first guess already on it (a round on a 5 x 5 bar's edge found it at
  // r9, a guess r along the edge, and not at r10, where that guess is not made). So does a
  // straight side edge on a cylinder or cone, a ruling of it (its UV curve a line): the surface
  // continues as itself past its domain (SmOffsetSurface's sm_ExtendAsItself), along that line.
  // A closed rim's seam is one: a disc's rim rounded past its 2 mm wall meets the seam's line
  // below the wall, and bounded to the seam the closed corner found no point there.
  SmCurve * pSide3D = pSideEdge->GetCurve();
  SmBoolean bAlongLine = FALSE;
  if (   pSide3D != NULL && pSide3D->IsLinear( SM_EFF_ZERO ) && pUVCurve->IsLinear( SM_EFF_ZERO )
      && ( pBSP->IsPlanar( SM_EFF_ZERO ) || pBSP->IsKindOf( SmCone_TYPE )))
    {
      SmExtent1d sEdgeIvl = sIntervals[0];
      SmPoint3d  sFrom, sTo;
      SER(pSide3D->EvaluatePoint(pSideEdge->GetInterval().GetMin(),sFrom));
      SER(pSide3D->EvaluatePoint(pSideEdge->GetInterval().GetMax(),sTo));
      double dLength = sFrom.DistanceBetween(sTo);
      double dReach  = smos_Max( m_dTraceReach, m_dSurfaceExtensionFactor * dLength );
      if ( dLength > SM_EFF_ZERO_SQRT )
        {
          // Only past the end away from the corner: the side edge leaves the corner's vertex
          // into the sector the blend fills, and a point of its line behind that vertex lies on
          // the far side of the filleted edge, where no rail of the blend can be. There, with
          // the faces continued as themselves, lies the mirror image of the ball (a knife edge's
          // ball touching the curved wall's continuation past the tip), and the solve took it.
          const SmEdge * pFilleted = GetEdgeuse(0)->GetEdge();
          auto bOnFilleted = [pFilleted]( const SmVertex * cpV )
            { return cpV != NULL && pFilleted != NULL && ( cpV == pFilleted->GetStartVertex() || cpV == pFilleted->GetEndVertex() ); };
          // A tangent seam T (the side edge's faces meet tangentially at that vertex) is the
          // exception: past the vertex T runs on in the rail's face, which continues across it
          // into the end face, and a rail that leaves the filleted edge's end on that side (a
          // concave round on the wall above a flat side running tangent into it) crosses T there.
          const SmBoolean bCornerAtStart = bOnFilleted( pSideEdge->GetStartVertex() );
          const SmBoolean bCornerAtEnd   = bOnFilleted( pSideEdge->GetEndVertex() );
          const SmVertex * cpCornerV     = bCornerAtStart ? pSideEdge->GetStartVertex() : pSideEdge->GetEndVertex();
          const SmBoolean bTangentSeam   =    ( bCornerAtStart || bCornerAtEnd )
                                           && sm_CheckTangentialEdgeuseAtVertex( pSideEU, cpCornerV, SM_RAD2DEG( GetTangencyTolerance() ));
          double dPast = dReach * sEdgeIvl.GetLength() / dLength;
          double dBefore = ( bCornerAtStart && !bCornerAtEnd && !bTangentSeam ) ? 0.0 : dPast;
          double dAfter  = ( bCornerAtEnd && !bCornerAtStart && !bTangentSeam ) ? 0.0 : dPast;
          sIntervals[0] = SmExtent1d( sEdgeIvl.GetMin() - dBefore, sEdgeIvl.GetMax() + dAfter );
          bAlongLine = TRUE;
        }
    }
  // The UV curve, a line, is evaluated on past its ends along the line: clamped to its own
  // interval (SmBSplineCurve evaluates out of bounds only where that is enabled), the solve
  // could not leave the edge's end, and found the point past it only from a guess already
  // there (a curve read back from a file, its flag unset, evaluated on).
  SmTemporaryChangeValue<SmBoolean> sAlongLine( pUVCurve->GetOutOfBoundsEnabled(), bAlongLine ? TRUE : pUVCurve->GetOutOfBoundsEnabled() );
  sIntervals[1] = SmExtent1d(sDomain.GetMin().x,sDomain.GetMax().x);
  sIntervals[2] = SmExtent1d(sDomain.GetMin().y,sDomain.GetMax().y);

  SmBoolean alPerData[16];
  SmTArray<SmBoolean> sPeriodicities(16,alPerData);

  sPeriodicities.Add(FALSE);
  sPeriodicities.Add(FALSE);
  sPeriodicities.Add(FALSE);

  double adGuessData[16];
  SmTArray<double> sGuessT(16,adGuessData);
  sGuessT.Add(dGuessEdgeParameter);
  sGuessT.Add(sUVGuess.x);
  sGuessT.Add(sUVGuess.y);

  SmSurface *pSurface1 = GetSurface(lRailIndex);
  SmSurface *pSurface2 = NULL;
  ULONG lDoSurf2Calcs  = 1; // Do surface2 calculations.
  ULONG lSurf1Index    = 1; // UV values of Surface1
  ULONG lSurf2Index;

  // given current inputs - place a good guess for upcoming solve in sGuessT
  SER(LoadInitialValues(lDoSurf2Calcs,
                        lRailIndex,
                        *pSurface1, 
                        lSurf1Index, 
                        pSurface2, 
                        lSurf2Index,
                        sIntervals, 
                        sPeriodicities, 
                        sGuessT));

  double sdData[16];
  SmTArray<double> sSolutionVector(16,sdData);

  // specify the rail/curve intersection evaluate function for the NewtonRaphson Solver
  SmRailUVCurveIntersectENFO sEvalFun(*pUVCurve, 
                                       pUVCurve->GetNaturalInterval(),
                                      *pSurface1, 
                                       lRailIndex,
                                      *this);

  // construct the NewtonRaphson multi-variable solver
  SmLocalSolveNd sLS(sEvalFun,          // in : Define Eqns to set to Zero (defines the DOF COUNT)
                     &sIntervals,       // in : NULL or sized:[N], N = this rEvalFun DOF count
                                        //      bounds on problem parameters.
                                        //      Set m_cpIntervals[i].SetUnbounded() for ivl whose prob params don't get clamped
                     &sPeriodicities) ; // in : NULL or sized:[N], N = this rEvalFun DOF count
                                        //      ivl[i] == TRUE = Interval[i] is periodic and values get wrapped
                                        //      ivl[i] == FALSE= Interval[i] is NOT periodic and values get clamped

  sLS.SetDesiredAccuracy(SM_EFF_ZERO_SQ);
  sLS.SetBoundaryHandler(SM_BH_TOTAL_BOUNDARY_HITS,8);
  sLS.SetMaximumIterations(30);

#ifdef SM_DEBUG_CODE
SmBoolean bDumpSols=FALSE;
if ( bDumpSols )
  { SM_DUMP_TARRAY( sGuessT ); }
#endif // SM_DEBUG_CODE

  // Solve it
  SmBoolean bFoundSolution;
  SER(sLS.SolveIt(sGuessT,
                  m_dThisApproxTol3d + pSideEdge->GetTolerance(),
                  bFoundSolution,
                  sSolutionVector));

  // when a solution was found
  if (bFoundSolution)
    {
#ifdef SM_DEBUG_CODE
if ( bDumpSols )
  { SM_DUMP_TARRAY( sSolutionVector ); }
#endif // SM_DEBUG_CODE

      // set outputs
      rbFoundIntersection = TRUE;
      rdEdgeuseParameter  = sSolutionVector[0];

      // set intersection point
      rTsectPnt.UVPos(lRailIndex).x   = sSolutionVector[lSurf1Index];
      rTsectPnt.UVPos(lRailIndex).y   = sSolutionVector[lSurf1Index+1];
      rTsectPnt.UVPos(1-lRailIndex).x = sSolutionVector[lSurf2Index];
      rTsectPnt.UVPos(1-lRailIndex).y = sSolutionVector[lSurf2Index+1];
      rTsectPnt.m_dCurveParameter     = rdEdgeuseParameter; // Just throw this into here for now.
      // Normally, m_dCurveParameter is the param of the intersection curve.
      // It's passed as output arg rdEdgeuseParameter.  Not sure why we're putting this here.

      // Set this here, otherwise it will be ignored in LoadInitialValues().
      SmIntersectionPointType eIntType = rTsectPnt.PointType();
      if ( eIntType == SM_IP_UNDEFINED || eIntType == SM_IP_UNKNOWN )
        { rTsectPnt.m_ePointType = SM_IP_CROSSING; }

      // Load other fillet surface variables being solved for into the
      // fillet point.  This should be sufficient to reproduce what is
      // happening at that point and to use as seed values for stepping.
      ULONG lCnt = 0;
      for (ULONG i=lSurf2Index+2; i<sSolutionVector.GetSize(); i++)
        {
          if (lCnt == SM_MAX_USER_DOUBLES) SER(SM_ERR);
          rTsectPnt.m_adUserDoubles[lCnt++] = sSolutionVector[i];
        }

    } // end found a solution check

  // all done
  return SM_SUCCESS;

} // end SmFilletSolver::RailEdgeuseIntersect

/*******************************************************************//**
PURPOSE: Given Rail UVPoints on both surfaces, compute Surface Point,
         Curve position and derivative values.

NOTES:
  For curve: rTsectPnt.m_dTangentPlaneAngleRad  = angle between surface normal vectors
***********************************************************************/
SmStatus SmFilletSolver::ComputeSurfaceValues
 (SmPoint2d    aUVValues[2], // in : rail UVPoint values to evaluate
  SmTsectPnt & rTsectPnt)    // out: container for surface point, curve position and derivative values
{
  // Reset offset values to be compatible to those solved for when computing rTsectPnt.
  SER(SetupOffsetValues(rTsectPnt));

  // First compute values on each surface

  // for both surfaces
  for (ULONG lSrf=0; lSrf<=1; lSrf++)
    {
      // get surf position, uTangent, and vTangent for given UVPoint
      SER(GetSurface(lSrf)->Evaluate1stDerivatives( aUVValues[lSrf], 
                                                    TRUE, 
                                                    TRUE,
                                                    rTsectPnt.SrfPos(lSrf), 
                                                    rTsectPnt.SrfDu(lSrf), 
                                                    rTsectPnt.SrfDv(lSrf) ));
      
      // We don't need 2nd surface derivs, but don't leave them uninitialized.
      rTsectPnt.SrfDuu(lSrf).Set( 0,0,0 );
      rTsectPnt.SrfDuv(lSrf).Set( 0,0,0 );
      rTsectPnt.SrfDvv(lSrf).Set( 0,0,0 );

      //      // get surf position, uTangent, and vTangent for given UVPoint
      //      SER(GetSurface(lSrf)->Evaluate2ndDerivatives( aUVValues[lSrf], TRUE, TRUE,
      //          rTsectPnt.SrfPos(lSrf), rTsectPnt.SrfDu(lSrf),  rTsectPnt.SrfDv(lSrf),
      //          rTsectPnt.SrfDuu(lSrf), rTsectPnt.SrfDuv(lSrf), rTsectPnt.SrfDvv(lSrf)));

      // store unitized surface normal in rTsectPnt.SrfNorm(lSrf)
      SmVector3d sNorm = rTsectPnt.SrfDu(lSrf) * rTsectPnt.SrfDv(lSrf);
      if ( sNorm.LengthSquared() < SM_EFF_ZERO_SQ )
        {
          SER(GetSurface(lSrf)->EvaluateNormal(aUVValues[lSrf],
                                               TRUE,
                                               TRUE,
                                               sNorm));
          if (sNorm.LengthSquared() < SM_EFF_ZERO_SQ)
            { SER(SM_ERR); }
        }
      (void)sNorm.Unitize(); // SolidCraft: a zero normal stays zero
      rTsectPnt.SrfNorm(lSrf) = sNorm;
      rTsectPnt.UVPos(lSrf)   = aUVValues[lSrf];  // Assign UV point values

    } // end iter both surfaces

    // Now compute 3D curve point and vectors
    // use surf0 point as curve point (used to use surface point averages)
    rTsectPnt.CrvPos() = rTsectPnt.SrfPos(0);

    // Compute angle between normal vectors
    SER(rTsectPnt.SrfNorm(0).AngleBetween(rTsectPnt.SrfNorm(1), rTsectPnt.m_dTangentPlaneAngleRad ));

    if (rTsectPnt.m_dTangentPlaneAngleRad > SM_PI/2.0)
      {
        rTsectPnt.m_dTangentPlaneAngleRad = SM_PI - rTsectPnt.m_dTangentPlaneAngleRad;
      }

    // all done
    return SM_SUCCESS;

} // end SmFilletSolver::ComputeSurfaceValues

/*******************************************************************//**
PURPOSE: Compute the point and derivative values of the point.

NOTES:

METHOD ---
     set rTsectPnt.UVDeriv(jSrf) = Surface UVDomain tangent in 3D xSect curve tangent direction

  For curve: rTsectPnt.m_dTangentPlaneAngleRad  = angle between surface normal vectors


***********************************************************************/
SmStatus SmFilletSolver::ComputePointValues
  (SmPoint2d aUVValues[2],           // in : rail UVcurve point values
   const SmExtent2d & crUVDomain1,   // in : fillet surface1 UVDomain
   const SmExtent2d & crUVDomain2,   // in : fillet surface2 UVDomain
   double dCurveTraceDirection,      // in : offset-surface xsect curve trace direction [+/-1]
   SmTsectPnt & rTsectPnt,           // out: intersection point gets updated surface values
   SmTsectPnt * pOptPreviousPnt)     // in : last intersection point, used to set
                                     //      this point's xsectCurve tangent direction when
                                     //      this surfaces are tangent at this point.
                                     //      NULL to ignore
{
  // set rTsectPnt.m_vSurfacePV position, uTangent, vTangent, and surface Normal values
  //     rTsectPnt.m_dTangentPlaneAngleRad  = angle between surface normal vectors
  SER(ComputeSurfaceValues(aUVValues,rTsectPnt));

#ifdef SM_DEBUG_CODE
int iDebugLevel = DebugLevel();
#endif

  // Now set the derivative of the spine curve, rTsectPnt.CrvDeriv().
  // That's really all we do in the rest of this method, except for
  // calling ComputeUVVectors() after we have that derivative.

  // set intersection curve tangent
  rTsectPnt.CrvDeriv() =  ( rTsectPnt.SrfNorm(0) * rTsectPnt.SrfNorm(1) ) * dCurveTraceDirection;

  // when surfaces are tangent to within tolerance
  if ( rTsectPnt.CrvDeriv().LengthSquared() < 10000.0*SM_EFF_ZERO_SQ )
    {
      // label the tangency
      rTsectPnt.m_ePointType = SM_IP_TANGENT_POINT;

      // when given a previous point - use its xSectCurve tangent as this points curve's tangent
      if ( pOptPreviousPnt != NULL )
        {
          rTsectPnt.CrvDeriv() = pOptPreviousPnt->CrvDeriv();
        }
      else
        {
          return SM_SUCCESS;
        }
    } // end tangent surface check

  // normalize the curve tangent
  SER( rTsectPnt.CrvDeriv().Unitize() );

  // Set SmTsectPt tangent values by finite differences
  // when not building a constant radius fillet or at a surf/surf tangency point.
  //
  //    Handle special case of starting a curve trace that walks off a DomainBoundary
  //    by negating the trace direction by setting all the SmTsectPt tangent values
  //    from a negative dStep value
  if (   OffsetRadiiCanChange()
      || rTsectPnt.m_ePointType == SM_IP_TANGENT_POINT )
   {
       // set rTsectPnt.UVDeriv(jSrf) = UVTangent Vector point in direction of 3d xSect curve
       SER( ComputeUVVectors( rTsectPnt ));

      // Take a small forward step.
      // SolidCraft change: chosen for this point. It was a function-local static, set once by the
      // first point the process computed here: every later point took that point's step (a
      // tangent point 1e-5 after an earlier crossing point, a crossing point 1e-3 after an
      // earlier tangent point), and a blend's result depended on the blends before it.
      double dStep = (rTsectPnt.m_ePointType == SM_IP_TANGENT_POINT)
                     ? 1.0e-3
                     : 1.0e-5;  // (1e-8 is too small: noise)

      // Take a step.
      // Note: won't find a step value when stepping off a surface domain boundary.
      // Note: StepSolve() sets only the uv positions in the output SmTsectPnt.

      SmTsectPnt sStepForw;
      SmBoolean bFoundForw, bClipped, bBoundaryHit;
      SER(StepSolve( rTsectPnt, dStep, crUVDomain1, crUVDomain2,
                    bFoundForw, bClipped, bBoundaryHit, sStepForw));

      // If that fails, try a bigger step, in the small-step case.
      if ( ! bFoundForw && dStep < 1.0e-4 )
        {
          double dBigStep = 1.0e-3;
          SER(StepSolve( rTsectPnt, dStep, crUVDomain1, crUVDomain2,
                  bFoundForw, bClipped, bBoundaryHit, sStepForw ));
          // If the bigger step worked, use it for the backward step as well.
          if ( bFoundForw )
            { dStep = dBigStep; }
        }

      // Collect the surface points now.  We can't wait until after the next call
      // to StepSolve(), becuase that adjusts the offset distances in the offset surfaces.
      SmPoint3d sSrf0PtForw, sSrf1PtForw;
      if ( bFoundForw )
        {
          SER( GetSurface(0)->EvaluatePoint( sStepForw.UVPos(0), sSrf0PtForw ));
          SER( GetSurface(1)->EvaluatePoint( sStepForw.UVPos(1), sSrf1PtForw ));
        }

      // Take a small backward step.
      dStep = -dStep;

      SmTsectPnt sStepBack;
      SmBoolean bFoundBack;
      SER( StepSolve( rTsectPnt, dStep, crUVDomain1, crUVDomain2,
           bFoundBack, bClipped, bBoundaryHit, sStepBack ));

      SmPoint3d sSrf0PtBack, sSrf1PtBack;
      if ( bFoundBack )
        {
          SER( GetSurface(0)->EvaluatePoint( sStepBack.UVPos(0), sSrf0PtBack ));
          SER( GetSurface(1)->EvaluatePoint( sStepBack.UVPos(1), sSrf1PtBack ));
        }

      // State: we tried a forward step and a backward step.

      // Use finite differences to set center-curve tangent value.
      // Note: If one direction failed, we must be at the end of a domain.
      // In that case, just revert to single-step forward or backwards;
      // at the end of a domain is probably just an extension that will
      // not be used in the end anyway.
      // (Otherwise we could take a second step in the direction that works
      // and use 2nd order forward or backward difference.)

      // SolidCraft change: at a tangent point the intersection curve ends or branches; it does
      // not run on through it, so a central difference across it averages the curve with
      // whatever the surfaces do past it. Where the surfaces touch only to within the model's
      // accuracy (a sketch with solver residues, 1e-10 mm), a step back past the point finds
      // the near-tangent crossing beyond it, a branch at 18 degrees to the curve, and the
      // averaged tangent made every step after the point fail the angle test, the trace
      // creeping on in tiny steps (repeated control points, an internal pole). Where they
      // touch exactly that step finds nothing and the forward difference is used. So it is
      // here: the forward step runs along the trace (the previous point's tangent).
      SmBoolean bTangentPoint = ( rTsectPnt.m_ePointType == SM_IP_TANGENT_POINT );
      SmPoint3d sStepPtForw, sStepPtBack;
      if ( bFoundForw && bFoundBack && ! bTangentPoint )
        {
            sStepPtForw = ( sSrf0PtForw + sSrf1PtForw ) / 2.0;
            sStepPtBack = ( sSrf0PtBack + sSrf1PtBack ) / 2.0;
            dStep = 2.0 * dStep;
        }
      else if ( bFoundForw )
        {
            sStepPtForw = ( sSrf0PtForw + sSrf1PtForw ) / 2.0;
            sStepPtBack = rTsectPnt.CrvPos();
        }
      else if ( bFoundBack )
        {
            sStepPtForw = rTsectPnt.CrvPos();
            sStepPtBack = ( sSrf0PtBack + sSrf1PtBack ) / 2.0;
        }
      else
        { SER(SM_ERR); }

      dStep = smos_Fabs( dStep );

      SmVector3d sDiff     = sStepPtForw - sStepPtBack;
      rTsectPnt.CrvDeriv() = sDiff / dStep;

      // Check this: Why Unitize()?  Then it's no longer a 1st derivative.
      // And there would be no need to div by dStep.
      SER( rTsectPnt.CrvDeriv().Unitize() );

      // In the tangent case, use the direction of the cross product of the
      // stepped-off normals to see whether the curve direction should be reversed.
      if (rTsectPnt.m_ePointType == SM_IP_TANGENT_POINT)
        {
          SmVector3d sNorm1, sNorm2, sCurveDirection; 
          if ( bFoundForw )
          {
              SER( GetSurface(0)->EvaluateNormal( sStepForw.UVPos(0), TRUE, TRUE, sNorm1 ));
              SER( GetSurface(1)->EvaluateNormal( sStepForw.UVPos(1), TRUE, TRUE, sNorm2 ));
              sCurveDirection = sNorm1 * sNorm2;
          }
          else // (If both False, we would have already returned.)
          {
              SER( GetSurface(0)->EvaluateNormal( sStepBack.UVPos(0), TRUE, TRUE, sNorm1 ));
              SER( GetSurface(1)->EvaluateNormal( sStepBack.UVPos(1), TRUE, TRUE, sNorm2 ));
              sCurveDirection = sNorm1 * sNorm2;
          }
          if ( sCurveDirection.LengthSquared() < SM_EFF_ZERO_SQ )
            {
              //cbi This seems harsh.  If this is ever hit, try harder for divergent normals.
              return SM_ERR;
            }
          sCurveDirection = sCurveDirection * dCurveTraceDirection;
          if ( sCurveDirection.Dot(rTsectPnt.CrvDeriv()) < 0.0 )
            {
              rTsectPnt.CrvDeriv() = - rTsectPnt.CrvDeriv();
            }
        }

#ifdef SM_DEBUG_CODE
      if ( iDebugLevel > 0 ) {
          if ( FALSE )
            { smgfx_Erase(); }
          SmVector3d &rCurrPt = rTsectPnt.CrvPos();
          if ( rCurrPt.x!=SM_UNDEF_DOUBLE ) {
              smgfx_SetLook(1,4, 0,0,1);
              rCurrPt.Draw(); sm_GraphicsLoop();
              rTsectPnt.CrvDeriv().Draw( &rCurrPt ); sm_GraphicsLoop();
              sm_GraphicsLoop();
          }
          sm_GraphicsLoop();
          if ( pOptPreviousPnt != NULL ) {
              SmVector3d &rPrevPt = pOptPreviousPnt->CrvPos();
              if ( rPrevPt.x!=SM_UNDEF_DOUBLE ) {
                  smgfx_SetLook(1,4, 0,1,0);
                  rPrevPt.Draw(); sm_GraphicsLoop();
                  pOptPreviousPnt->CrvDeriv().Draw(& rPrevPt ); sm_GraphicsLoop();
                  sm_GraphicsLoop();
              }
          }
          smgfx_SetLook(1,2, 1,0,0);
          GetSurface(0)->DrawAt( rTsectPnt.UVPos(0), 1 ); sm_GraphicsLoop();
          GetSurface(1)->DrawAt( rTsectPnt.UVPos(1), 1 ); sm_GraphicsLoop();
          smgfx_SetLook(1,4, 0,1,1);
          sStepPtForw.Draw(); sm_GraphicsLoop();
          sStepPtBack.Draw(); sm_GraphicsLoop();
          sm_GraphicsLoop();
          if ( iDebugLevel > 50 ) {
              smgfx_SetLook(1,2, 0,0,1); GetSurface(0)->DrawUV(4,4); sm_GraphicsLoop();
              smgfx_SetLook(1,2, 0,0,1); GetSurface(1)->DrawUV(4,4); sm_GraphicsLoop();
              sm_GraphicsLoop();
          }

      }
#endif

    } // end if using finite differences: variable radius or tangent.


  rTsectPnt.m_ePointType = SM_IP_CROSSING;

  // use center-curve tangent to compute surface uvTangents
  SER( ComputeUVVectors( rTsectPnt ) );

  return SM_SUCCESS;

} // end SmFilletSolver::ComputePointValues

/*******************************************************************//**
PURPOSE: Compute uv tangent vector for each surface pointing in the
            direction of the 3d XSect curve.

NOTES:
  for both surfaces, set
     rTsectPnt.UVDeriv(jSrf) = UVTangent Vector point in direction of 3d xSect curve
***********************************************************************/
SmStatus SmFilletSolver::ComputeUVVectors
  (SmTsectPnt & rTsectPnt)                 // out: container for computed surface UVTangent vectors
{
  // for both surfaces
  for (ULONG jSrf=0; jSrf<=1; jSrf++)
    {
      SER(smsurf_DropVectors(rTsectPnt.SrfDu(jSrf),      // in : Surface U_dir tangent
                             rTsectPnt.SrfDv(jSrf),      // in : Surface V_dir tangent
                             1,                          // in : vector count
                             &rTsectPnt.CrvDeriv(),      // in : vector to drop
                             &rTsectPnt.UVDeriv(jSrf))); // out: UV parameter value of given vectors
    }

  return SM_SUCCESS;

} // end SmFilletSolver::ComputeUVVectors

/*******************************************************************//**
PURPOSE:

NOTES:
***********************************************************************/
class SmRailRailIntersectENFO : public SmEvalNFunctionsObject
{
protected:
    ULONG            m_lRailIndex;
    SmFilletSolver & m_crFilletSolver;
    ULONG            m_lOtherRailIndex;
    SmFilletSolver & m_crOtherFilletSolver;
public:
    SmRailRailIntersectENFO(ULONG            lRailIndex,
                            SmFilletSolver & crFilletSolver,
                            ULONG            lOtherRailIndex,
                            SmFilletSolver & crOtherFilletSolver)

      : m_lRailIndex(lRailIndex), m_crFilletSolver(crFilletSolver),
        m_lOtherRailIndex(lOtherRailIndex), m_crOtherFilletSolver(crOtherFilletSolver) {}
    virtual ~SmRailRailIntersectENFO() {}
    virtual SmStatus Evaluate(const SmTArray<double> & crX,             // in : x of Ax=F
                              SmTArray<double>       & rF,              // out: F of Ax=F function values,
                              SmMatrix               * pOptJacobian,    // out: Partial derivatives of the functions.
                              SmBoolean              & rbFoundAnswer);  // out: Not always used, when used
                                                                        //      TRUE = converged (F members are within tolerance of 0.0
                                                                        //      FALSE= Not Used or Not Converged
}; // end class SmRailRailIntersectENFO

/****************************************************************
PURPOSE:

NOTES:
****************************************************************/
SmStatus SmRailRailIntersectENFO::Evaluate
  (const SmTArray<double> & crX,           // in : x of Ax=F
   SmTArray<double>       & rF,            // out: F of Ax=F function values,
   SmMatrix               * pOptJacobian,  // out: Partial derivatives of the functions.
   SmBoolean              & rbFoundAnswer) // out: TRUE =
                                           //      FALSE=

{
    rbFoundAnswer = FALSE;

    if (pOptJacobian) {
        for (ULONG ll=0; ll<pOptJacobian->GetNumRows(); ll++) {
            for (ULONG mm=0; mm<pOptJacobian->GetNumColumns(); mm++) {
                (*pOptJacobian)[ll][mm] = 0.0;
            }
            rF[ll] = 0.0;
        }
    }

    ULONG lNumEquations = 0;
    ULONG lNumParameters = 2;
    SmSurface *pSurface = m_crFilletSolver.GetSurface(m_lRailIndex);
    ULONG lSurf1Offset = 0;
    SmSurface * pSurface2 = NULL;
    ULONG lSurf2Offset = 0;
    SmBoolean bFoundAnswer;

    if (m_crFilletSolver.LoadJacobian(lNumEquations, lNumParameters,
        m_lRailIndex, pSurface, lSurf1Offset, pSurface2, lSurf2Offset,
        crX, rF, pOptJacobian, bFoundAnswer) != SM_SUCCESS) {
        return SM_ERR;
    }

    pSurface2 = NULL;
    SmBoolean bFoundAnswer2;
    if (m_crOtherFilletSolver.LoadJacobian(lNumEquations, lNumParameters,
        m_lOtherRailIndex, pSurface, lSurf1Offset, pSurface2, lSurf2Offset,
        crX, rF, pOptJacobian, bFoundAnswer2) != SM_SUCCESS) {
        return SM_ERR;
    }


    // Jacobian answer looks good now check our answer to see if we return
    // TRUE for rbFoundAnswer
    if (bFoundAnswer && bFoundAnswer2) {
        rbFoundAnswer = TRUE;
    }

#ifdef SM_DEBUG_CODE
    if ( m_crFilletSolver.DebugLevel() > 60 ) {
        SmPoint2d sUV(crX[0],crX[1]);
        SmPoint3d sPnt;
        ((SmOffsetSurface*)pSurface)->SetOffsetDistance(0.0);
        SER(pSurface->EvaluatePoint(sUV,sPnt));
        smgfx_SetColor(1,0,0);
        sPnt.Draw();
        sm_GraphicsLoop();
        if (pOptJacobian) {
            pOptJacobian->Dump();
            SmMatrix sTmpJ(8,8);
            SmTArray<double> sTmpX(crX);
            SER(ComputeJacobian(sTmpX,rF,sTmpJ));
            sTmpJ.Dump();
            smos_WriteBuffer(_T(" X --- \n"));
            crX.Dump();
            smos_WriteBuffer(_T(" F --- \n"));
            rF.Dump();
        }
    }
#endif

    return SM_SUCCESS;

} // end SmRailRailIntersectENFO::Evaluate

/*******************************************************************//**
PURPOSE: Intersect two rails on the same surface from two different fillets.

NOTES:
   Sets up a 6x6 nonlinear system.  Variables are:
   [0],[1] uv of common surface
   [2],[3] uv of first side surface
   [4],[5] uv of second side surface
   Objective functions:
   [0],[1],[2] xyz of diff between point in common surface and point in first side surface
   [3],[4],[5] xyz of diff between point in common surface and point in second side surface
***********************************************************************/
SmStatus SmFilletSolver::RailRailIntersect
  (const SmPoint2d & crUVGuess,             // in : FiletSolver->RailSurface UVPoint guess
   ULONG             lRailIndex,            // in : 1st Target Rail index in this FilletSolver
   SmFilletSolver  * pOtherFillet,          // in : Other rail FilletSolver
   ULONG             lOtherRailIndex,       // in : 2nd Target Rail index in other FilletSolver
   SmBoolean       & rbFoundIntersection,   // out: TRUE = Found an intersection Point
   SmTsectPnt      & rTsectPnt,             // out: 
   SmTsectPnt      & rOtherTsectPnt,        // out: 
   SmBoolean         bHoldEndRadii)         // in : hold variable radii at their end values (a corner calculated as constant radius)
{
  // init output
  rbFoundIntersection = FALSE;

  // check inputs
  if(   GetSurface(0)               == NULL 
     || GetSurface(1)               == NULL  
     || pOtherFillet->GetSurface(0) == NULL
     || pOtherFillet->GetSurface(1) == NULL) { return SM_ERR; }

#ifdef SM_DEBUG_CODE
int iDebugLevel = DebugLevel();
static ULONG lCount      = 1 ; lCount++ ;
static ULONG lDebugCount = 0 ;
  if ( iDebugLevel > 0 || lDebugCount == lCount) 
    {
      SmSurface * pSurf = GetSurface(lRailIndex);
      SmFace    * pFace = (SmFace *)pSurf->GetFace() ;
      SmBrep    * pBrep = pFace ? pFace->GetBrep() : NULL ;
      SmPoint3d sCornerPt;
      pSurf->EvaluatePoint(crUVGuess, sCornerPt);
      
      smgfx_Erase();
      smgfx_SetLook(1,2, 0,0,1) ; if(pBrep) pBrep->Draw(TRUE) ; sm_GraphicsLoop() ;
      if ( iDebugLevel > 20 ) {
          smgfx_SetLook(1,2, 0,1,1) ; if(pSurf) pSurf->DrawUV(1,1) ; sm_GraphicsLoop() ;
          smgfx_SetLook(1,2, 0,0,0) ; if(pFace) pFace->Draw(SM_DM_CROSSHATCH); sm_GraphicsLoop() ;
      }
      smgfx_SetLook(3,4, 1,0,0) ; GetEdgeuse(lRailIndex)->Draw(); sm_GraphicsLoop() ;
      smgfx_SetLook(4,5, 1,0,1) ; pOtherFillet->GetEdgeuse(lOtherRailIndex)->Draw(); sm_GraphicsLoop() ;
      smgfx_SetLook(5,7, 0,1,0) ; sCornerPt.Draw(); sm_GraphicsLoop() ;
      sm_GraphicsLoop() ;
    }
#endif

  // Set up Solver Size
  ULONG lJacobianSize       = GetJacobianSize();
  ULONG lOtherJacobianSize  = pOtherFillet->GetJacobianSize();
  lJacobianSize             = lJacobianSize + lOtherJacobianSize - 2;

  // compute and save extended domain size for TargetRail->Surface
  const SmSurface * pBSP      = GetSurface(lRailIndex)->GetBaseSurface();
  SmExtent2d        sDomain   = pBSP->GetNaturalUVDomain();
  SetupOffsetExtension(pBSP,sDomain);
  SmExtentNd sIntervals(lJacobianSize);
  sIntervals[0] = SmExtent1d(sDomain.GetMin().x,sDomain.GetMax().x);
  sIntervals[1] = SmExtent1d(sDomain.GetMin().y,sDomain.GetMax().y);

  // mark periodicities for TargetRail->Surface 
  SmBoolean alPerData[16];
  SmTArray<SmBoolean> sPeriodicities(16,alPerData);
  sPeriodicities.Add(FALSE);
  sPeriodicities.Add(FALSE);

  // Set Solver GuessPoint
  double adGuessData[16];
  SmTArray<double> sGuessT(16,adGuessData);
  sGuessT.Add(crUVGuess.x);
  sGuessT.Add(crUVGuess.y);

  // Temporarily replace the variable-radius functions with a constant. [Reg Iter:306]
  // See notes in the class header.
  // SolidCraft change: only at a corner calculated as constant radius, as the rail/side-edge
  // search holds it (SmFilletVertex::CalcRailIntEdgeuse): there the radius laws are then mapped
  // to end at the rail ends found (SmFilletCorner::CalcCornerVertGeom). Held at any other corner,
  // the rails' meeting point is found at the laws' end radii, short of the edges' ends where the
  // laws are smaller: two variable rounds 0.5 to 2 meeting at a box's corner met at y = 8 with
  // radius 2 where the law gives 1.7, and the trace from that rail end, held to the law by the
  // fillet solved inside each span (4438859), found no span it could take (1001).
  SmTemporaryRadiusChange sTempRadThis ( bHoldEndRadii ? this : NULL, lRailIndex, pOtherFillet->GetEdgeuse( lOtherRailIndex ) );
  SmTemporaryRadiusChange sTempRadOther( bHoldEndRadii ? pOtherFillet : NULL, lOtherRailIndex, this->GetEdgeuse( lRailIndex ) );

  // Locals for LoadInitialValues() 
  SmSurface * pSurface1     = GetSurface(lRailIndex);
  SmSurface * pSurface2     = NULL;
  ULONG       lDoSurf2Calcs = 1;  // Do surface2 calculations.
  ULONG       lSurf1Index   = 0;  // index into sGuessT of srf 1's UV values
  ULONG       lSurf2Index;

  // Append a good lRailIndex->OtherSurface UVguessPoint into sGuessT
  SER(LoadInitialValues(lDoSurf2Calcs,    // not used
                        lRailIndex,       // in : target rail index
                       *pSurface1,        // in : offset Surface for target rail
                        lSurf1Index,      // in : index into sGuessT of srf 1's UV values
                        pSurface2,        // i/o: ptr to other offset surface, NULL on input.
                        lSurf2Index,      // out: index into rGuessT of srf2's UV values
                        sIntervals,       // i/o: u and v domains of srf2 appended       
                        sPeriodicities,   // i/o: u and v periodicities of srf2 appended.
                        sGuessT));        // i/o: u and v guesses of srf2 appended       

  // Append a good lOtherRailIndex->OtherSurface UVguessPoint into sGuessT
  pSurface1                 = pOtherFillet->GetSurface(lOtherRailIndex);
  SmSurface *pOtherSurface2 = NULL;
  ULONG lOtherSurf2Index;
  SER(pOtherFillet->LoadInitialValues(lDoSurf2Calcs,   lOtherRailIndex,
                                      *pSurface1,      lSurf1Index, 
                                       pOtherSurface2, lOtherSurf2Index,
                                       sIntervals, sPeriodicities, sGuessT));

  double sdData[16];
  SmTArray<double> sSolutionVector(16,sdData);

  // Set Solver Manager
  SmRailRailIntersectENFO sEvalFun(lRailIndex,*this,lOtherRailIndex,*pOtherFillet);
  SmLocalSolveNd sLS(sEvalFun,          // in : Define Eqns to set to Zero (defines the DOF COUNT)
                     &sIntervals,       // in : NULL or sized:[N], N = this rEvalFun DOF count
                                        //      bounds on problem parameters.
                                        //      Set m_cpIntervals[i].SetUnbounded() for ivl whose prob params don't get clamped
                     &sPeriodicities) ; // in : NULL or sized:[N], N = this rEvalFun DOF count
                                        //      ivl[i] == TRUE = Interval[i] is periodic and values get wrapped
                                        //      ivl[i] == FALSE= Interval[i] is NOT periodic and values get clamped

  sLS.SetDesiredAccuracy(SM_EFF_ZERO_SQ);
  sLS.SetBoundaryHandler(SM_BH_TOTAL_BOUNDARY_HITS, 8);
  sLS.SetMaximumIterations(30);

  // Solve it
  SmBoolean bFoundSolution;
  SER(sLS.SolveIt(sGuessT, 
                  m_dThisApproxTol3d + pOtherFillet->m_dThisApproxTol3d,
                  bFoundSolution,
                  sSolutionVector));

  // save found solution
  if (bFoundSolution) 
    {
      rbFoundIntersection = TRUE;

      rTsectPnt.UVPos(  lRailIndex).x = sSolutionVector[lSurf1Index];
      rTsectPnt.UVPos(  lRailIndex).y = sSolutionVector[lSurf1Index+1];
      rTsectPnt.UVPos(1-lRailIndex).x = sSolutionVector[lSurf2Index];
      rTsectPnt.UVPos(1-lRailIndex).y = sSolutionVector[lSurf2Index+1];

      // Load other fillet surface variables being solved for into the
      // fillet point.  This should be sufficient to reproduce what is
      // happening at that point and to use as seed values for stepping.
      // Note these values should be between the the end of Surface 2 of Fillet 1
      // and the beginning of Surface 2 of Fillet 2
      ULONG lCnt = 0;
      for (ULONG i=lSurf2Index+2; i<lOtherSurf2Index; i++) 
        {
          if (lCnt == SM_MAX_USER_DOUBLES) SER(SM_ERR);
          rTsectPnt.m_adUserDoubles[lCnt++] = sSolutionVector[i];
        }

      rOtherTsectPnt.UVPos(  lOtherRailIndex).x = sSolutionVector[lSurf1Index];
      rOtherTsectPnt.UVPos(  lOtherRailIndex).y = sSolutionVector[lSurf1Index+1];
      rOtherTsectPnt.UVPos(1-lOtherRailIndex).x = sSolutionVector[lOtherSurf2Index];
      rOtherTsectPnt.UVPos(1-lOtherRailIndex).y = sSolutionVector[lOtherSurf2Index+1];

      // Load values into rOtherTsectPnt.  These values should be the last
      // on the list after the surface2 values from fillet 2.
      lCnt = 0;
      for (ULONG ii=lOtherSurf2Index+2; ii<sSolutionVector.GetSize(); ii++) 
        {
          if (lCnt == SM_MAX_USER_DOUBLES) SER(SM_ERR);
          rOtherTsectPnt.m_adUserDoubles[lCnt++] = sSolutionVector[ii];
        }
    } // end found a solution check

  // all done
  return SM_SUCCESS;

} // end SmFilletSolver::RailRailIntersect

/*******************************************************************//**
PURPOSE: Set extension flags to TRUE for each end if there is not
    a blend being done at the corner.

NOTES:
***********************************************************************/
void SmFilletSolver::SetExtensionFlags
  ()
{
  // locals
  SmEdgeuse * pEU0   = GetEdgeuse(0);
  SmEdgeuse * pEU1   = GetEdgeuse(1);
  SmVertex  * pVert0 = pEU0->GetVertexuse()->GetVertex();
  SmVertex  * pVert1 = pEU1->GetVertexuse()->GetVertex();

  // swap vertices for opposite oriented edgeuses
  if (pEU0->GetOrientation() == SM_OT_OPPOSITE)
    {
      SmVertex * pTempVert = pVert0;
      pVert0               = pVert1;
      pVert1               = pTempVert;
    }

  // get fillet information for both vertices
  SmFilletExecutive * pFilExec = GetFilletExecutive();
  SmFilletCorner    * pCorner0 = pFilExec->GetFilletCornerOfVertex(pVert0);
  SmFilletCorner    * pCorner1 = pFilExec->GetFilletCornerOfVertex(pVert1);

  // extend ends connected to vertices which are not blended or degenerate
  m_bExtendBefore = (   !pCorner0->IsBlendingCorner()
                     && !pCorner0->IsDegenerateCorner()) ;
  m_bExtendAfter  = (   !pCorner1->IsBlendingCorner()
                     && !pCorner1->IsDegenerateCorner()) ;

} // end void SmFilletSolver::SetExtensionFlags


/*******************************************************************//**
PURPOSE: Restore the offset surfaces & directions

NOTES: Copy OffsetSurface and Orientation values from
   input SmFilletSolverOffsetsData object into this SmFilletSolver
***********************************************************************/
void SmFilletSolver::SetOffsetsData
  (SmFilletSolverOffsetsData & rOffsetsData) // in : source rOffsetsData object
{
  // for both sides of the fillet edge
  for (ULONG i=0; i<2; i++)
    {
      m_dOrientations[i]             = rOffsetsData.m_dOrientations[i] ;
      SmOffsetSurface * pOrigOffSurf = rOffsetsData.m_pSurfaces[i] ;
      if (GetSurface(i) != pOrigOffSurf)
        {
          SmObjDelete sDelete( GetSurface(i) );
          SetSurface( i, pOrigOffSurf );
        }
    } // end iter both filletEdge sides

  return;

} // end void SmFilletSolver::SetOffsetsData

/*******************************************************************//**
PURPOSE: Save the offset surfaces & directions

NOTES: This method (as well as SetOffsetsData) should
    only be used when processing RollOver. Original offset surfaces &
    dirctions (if existed) can be saved into the SmFilletSolverOffsetsData class
***********************************************************************/
void SmFilletSolver::GetOffsetsData
  (SmFilletSolverOffsetsData & rOffsetsData)
{
  // for both filletEdge sides
  for (ULONG i=0; i<2; i++)
    {
      // copy orientation and Surface data
      rOffsetsData.m_dOrientations[i] = m_dOrientations[i];
      rOffsetsData.m_pSurfaces[i]     = GetSurface(i);
    } // end iter both FilletEdge sides

  // all done
  return;

} // end void SmFilletSolver::GetOffsetsData

/*******************************************************************//**
PURPOSE: Determine the amount an offset surface can be extended for
         evaluation purposes during the fillet process.

NOTES:
    This function only computes the extended domain, it does not actually
    build the extended SmOffsetSurface::m_pExtendedSurface object.  That 
    happens as needed during the normal evaluation sequence.
    
    when pSurf is not periodic and has no singular boundaries,
    extends rDomain by m_dSurfaceExtensionFactor * rDomain.GetSize()
    at both its min and max boundaries.

    Does not extend the rDomain through a singular boundary.

    Formerly, this did no extension in a direction that is marked as periodic,
    but that caused errors.  Filleting often needs to converge slightly out of
    the domain.  In practice, extending by a smaller amount works well.
    [bd, 09-Nov-05; 050718]

***********************************************************************/
void SmFilletSolver::SetupOffsetExtension
  (const SmSurface * pSurf,         // in : target surface
   SmExtent2d      & rDomain)       // i/o: target surface domain of interest
{
  // get length of extensions
  SmVector2d sSize = GetDomainExtension( pSurf, rDomain ); // SolidCraft change: and the trace's reach

  // Extend by a smaller amount if periodic.
  // Note, this used to be 0; see also Fillet regressions 212, 213, 214.
  // It is possible that we don't want to do this if the fillet surface
  // itself is periodic, i.e., an Nx1 case where a single fillet surface
  // joins itself smoothly.
  double dPeriodicExtensionFactor =  0.2;
  double dLinearExtensionFactor   = 50.0;

static int cbiLinear=0;

  if (pSurf->IsPeriodic(rDomain,SM_SP_U))
    {
      //if (!m_bExtendBefore && !m_bExtendAfter)
      sSize.x *= dPeriodicExtensionFactor;
    }
  else if ( cbiLinear>0 &&    pSurf->GetDegree( SM_SP_U ) == 1 )
    {
      sSize.x *= dLinearExtensionFactor;
    }
  if (pSurf->IsPeriodic(rDomain,SM_SP_V))
    {
      //if (!m_bExtendBefore && !m_bExtendAfter)
      sSize.y *= dPeriodicExtensionFactor;
    }
  else if ( cbiLinear>0 &&    pSurf->GetDegree( SM_SP_V ) == 1 )
    {
      sSize.y *= dLinearExtensionFactor;
    }

  // process the min domain boundaries
  SmSurfParamType eSingDir;
  if (!pSurf->IsSingularity(rDomain.GetMin(),eSingDir))
    {
      // when surface is not singular - extend the min domain boundary
      rDomain.AddPoint2d(rDomain.GetMin() - sSize);
    }
  else if (eSingDir == SM_SP_U)
    {
      // singular in U direction only extend the u direction
      SmPoint2d sNewMin  = rDomain.GetMin();
      sNewMin.x         -= sSize.x;
      rDomain.AddPoint2d(sNewMin);
    }
  else
    { // eSingDir == SM_SP_V
      // singular in V direction only extend the v direction
      SmPoint2d sNewMin = rDomain.GetMin();
      sNewMin.y -= sSize.y;
      rDomain.AddPoint2d(sNewMin);
    }

  // process the max domain boundaries
  if (!pSurf->IsSingularity(rDomain.GetMax(),eSingDir))
    {
      // when surface is not singular - extend the max domain boundary
      rDomain.AddPoint2d(rDomain.GetMax() + sSize);
    }
  else if (eSingDir == SM_SP_U)
    {
      // singular in U direction only extend the u direction
      SmPoint2d sNewMax  = rDomain.GetMax();
      sNewMax.x         += sSize.x;
      rDomain.AddPoint2d(sNewMax);
    }
  else
    { // eSingDir == SM_SP_V
      // singular in V direction only extend the v direction
      SmPoint2d sNewMax = rDomain.GetMax();
      sNewMax.y += sSize.y;
      rDomain.AddPoint2d(sNewMax);
    }

} // end void SmFilletSolver::SetupOffsetExtension

/*******************************************************************//**
PURPOSE: Set offset surface data to agree with that in the given
  FilletGeom, if different.

NOTES:
  If this routine does change one or both surfaces, it does not delete
  the existing surface, it simply overwrites its pointer.
  It is assumed that the caller has backed up the offset data:
    // Backup original offsets data
    SmFilletSolverOffsetsData sOrigOffsets;
    GetOffsetsData( sOrigOffsets );
  and will retore them when done:
    SetOffsetsData( sOrigOffsets );

  SetOffsetsData() will delete any new offset surfaces.

LIMITATIONS: If the surfaces are different, they must be in adjacent
  faces, i.e., faces that share a common edge.  If that is not the
  case, the offset direction could be backwards.
***********************************************************************/

SmStatus SmFilletSolver::SetOffSurfsToFilletGeom( SmFilletGeom *pFG )
{
int iDebugLevel = 0;
#ifdef SM_DEBUG_CODE
iDebugLevel = DebugLevel();
#endif

  ULONG lSrfIdx;
  for ( lSrfIdx = 0; lSrfIdx < 2; lSrfIdx++ )
  {
      const SmSurface *pOurBaseSurf = GetSurface( lSrfIdx )->GetBaseSurface();
      SmSurface *pRailBaseSurf = pFG->GetRail(lSrfIdx)->GetOriginalFace()->GetSurface();

      if ( pRailBaseSurf != pOurBaseSurf )
      {
          // We have to find out whether FG's surface has the
          // same normal direction as FS's surface.
          SmBoolean bOppositeNormals = FALSE;

          SmFace *pFSFace = SM_CAST_PTR( SmFace, pOurBaseSurf->GetFace() );
          SmFace *pFGFace = pFG->GetRail(lSrfIdx)->GetOriginalFace();
          if ( pFSFace == NULL || pFGFace == NULL )
              continue;

          SmTArray<SmEdge*> aFSEdges, aFGEdges, aCommonEdges;
          pFSFace->GetEdges( aFSEdges );
          pFGFace->GetEdges( aFGEdges );
          aFSEdges.FindCommonElements( aFGEdges, aCommonEdges );
          SM_ASSERT( aCommonEdges.GetSize() == 1 );

          // If more than one, just use the first.
          if ( aCommonEdges.GetSize() > 0 )
          {
              sm_TestManifoldTangentRollover( aCommonEdges[0], bOppositeNormals, m_dTangencyTolerance, iDebugLevel );
          }

          // Create new offset surface
          double dOffsetDist = m_pSurfaces[ lSrfIdx ]->GetOffsetDistance();
          if ( bOppositeNormals )
          {
              dOffsetDist *= -1.0;
              m_dOrientations[ lSrfIdx ] *= -1.0;
          }
          SmOffsetSurface *pNewOffsetSurface = new (m_crContext) SmOffsetSurface(
                            dOffsetDist, *pRailBaseSurf );
          pNewOffsetSurface->SetAllowExtension( TRUE );
          SetSurface( lSrfIdx, pNewOffsetSurface );
      }
  }
  return SM_SUCCESS;

} // end SmFilletSolver::SetOffSurfsToFilletGeom

/*******************************************************************//**
PURPOSE: Split a Fillet Geom Object into two pieces
            at a side-edge/Rail intersection.

NOTES: 
   The input pSplitFilVU is one of the to-be-split FilletGeom's
   vertexuses, and lies on a side edge.

   The original FilletGeom pointer will be the low end of the original
   FilletGeom, as output from TopologySplit().
   Either the original or the new FilletGeom will lie on the face
   on the other side of the split-Vertexuse.  If lWhichEnd is 0,
   the low end (original FG) will be on the other face,
   otherwise the new FG will be on the other face.

***********************************************************************/
SmStatus SmFilletSolver::SplitFilletGeomAtSideEdge (
    SmFilletGeom * pFilletGeom,      // in: FG to be split
    SmFilletVertexuse * pSplitFilVU, // in: FVU on Rail where side edge hits
    ULONG lWhichEnd,                 // in: pSplitFilVU is at low (0)
                                     //     or high (1) end of pFilletGeom.
    SmFilletGeom * *ppRetFilletGeom  // out, opt: ptr to new FG
  )
{
int iDebugLevel = 0;
#ifdef SM_DEBUG_CODE
iDebugLevel = DebugLevel();

SmBoolean bDebug_VUs=FALSE;
    if ( bDebug_VUs )
     { pFilletGeom->DumpFilletVUs(); }
#endif  // SM_DEBUG_CODE

    // Get the rail being split.
    ULONG lRailIdx;
    SER( pFilletGeom->GetRailIndexOfVertexuse( pSplitFilVU, lRailIdx ));
    SmFilletEdge *pSplitRail = pFilletGeom->GetRail( lRailIdx );

    // The Fillet Vertex contains a PointClass object whose object is the side Edgeuse that we intersect.
    SmFilletVertex *pFilVtx = SM_CAST_PTR( SmFilletVertex, pSplitFilVU->GetVertex() );
    NER( pFilVtx );

    SmEdgeuse *pSideEU = SM_CAST_PTR( SmEdgeuse, pFilVtx->GetPointClassObject() );
    NER( pSideEU );
    SmEdge *pSideEdge = pSideEU->GetEdge();

    // Should be the same Face on both sides of the side Edge.
    SmFace *pOrigFace = pSplitRail->GetOriginalFace(); NER(pOrigFace);
    SmTArray<SmFace*> sFaces;
    pSideEdge->GetFaces( sFaces );

    SM_ASSERT( sFaces.GetSize() == 1 );
    SM_ASSERT( sFaces[0] == pOrigFace );

    // Find the adjacent face -- the one on the other side of
    // the side edge that we hit.

    // Can't do this: the Rail Edge is a Wire, and so is the cross-section Edge.
    //    SmFace *cbiFace = pSplitFilVU->GetEdgeuse()->GetRadial()->GetFace();
    //
    SmFace *pAdjFace = NULL;

    // If there is only one Face for the side edge, it must be closed.
    SmBoolean bClosedSideFace = FALSE;
    if ( sFaces.GetSize() == 1 )
    {
        // Closed Edge to start with.
        SM_ASSERT( sFaces[0] == pOrigFace );
        pAdjFace = sFaces[0];
        bClosedSideFace = TRUE;
    }
    else
    {
        // Find the other Face.
        // Can't use its Radial EU because it's a Wire Edge.
        SM_ASSERT( sFaces.GetSize() == 2 ); // Otherwise spine or wire.

        for ( ULONG jj=0; jj<sFaces.GetSize(); jj++ )
        {
            if ( sFaces[jj] != pOrigFace )
            {
                pAdjFace = sFaces[jj];
                break;
            }
        }
    }

    NER(pAdjFace);

    SmSurface * pAdjSurface = pAdjFace->GetSurface();


    // Create new offset surface for the adjacent face.

    // First get offset distance and relative direction.
    double dOffsetDist = GetSurface( lRailIdx )->GetOffsetDistance();

    // Same normals?  This routine will tell us that:
    SmBoolean bOppositeNormals = FALSE;
    if ( pAdjFace != pOrigFace )  // Mild assumption: closed surface is not Moebius...
    {
        if ( ! sm_TestManifoldTangentRollover( pSideEdge, bOppositeNormals, m_dTangencyTolerance, iDebugLevel ) )
          { SER( SM_ERR_INVALID_INPUT ); }
    }

    if ( bOppositeNormals )
    {
        dOffsetDist *= -1.0;
        //cbi: no: m_dOrientations[lRailIdx] *= -1.0;
    }

    SmOffsetSurface *pAdjOffSurf = new (m_crContext) SmOffsetSurface(
            dOffsetDist, *pAdjSurface );
    pAdjOffSurf->SetAllowExtension(TRUE);

    // If the base part is the rollover part, set its type before
    // call to split: makes a difference inside that routine.
    if ( lWhichEnd == 0 )
    {
        // Original FG becomes rollover, on adjacent face.
        pFilletGeom->SetFilletGeomType( SM_FG_TANGENT_ROLLOVER );
    }

#ifdef SM_DEBUG_CODE
    if ( bDebug_VUs )
      { pFilletGeom->DumpFilletVUs(); }
#endif

    // Do the split.
    SmFilletGeom *pNewFilletGeom;
    SER( pFilletGeom->TopologySplit( pNewFilletGeom ) );

#ifdef SM_DEBUG_CODE
    if ( bDebug_VUs )
      { pFilletGeom   ->DumpFilletVUs();
        pNewFilletGeom->DumpFilletVUs();
      }
#endif

    if ( ppRetFilletGeom != NULL )
        *ppRetFilletGeom = pNewFilletGeom;


    // The new fillet geom comes out as the 'high end',
    // i.e., the high end of the original meets the low end
    // of the new one.

    // Update offset surfaces, orig face, and type in the two FG's.
    // Copy offset surface from non-split side of old to new.
    // Note, FG::SetOffsetSurface() deletes what's there, if anything.

    SmSurface *pSurf = NULL;
    SmOffsetSurface *pOffSurf = NULL;
    SmOffsetSurface *pOrigOffSurf = pFilletGeom->GetOffsetSurface( 1-lRailIdx );
    if ( pOrigOffSurf != NULL ) {
        SER( pOrigOffSurf->Copy( m_crContext, pSurf ) );
        pOffSurf = SM_CAST_PTR( SmOffsetSurface, pSurf );
    }
    SER( pNewFilletGeom->SetOffsetSurface( 1-lRailIdx, pOffSurf ) );

    if ( lWhichEnd == 0 )
    {
        pFilletGeom->GetRail(lRailIdx)->SetOriginalFace( pAdjFace );

        pSurf = pOffSurf = NULL;
        pOrigOffSurf = pFilletGeom->GetOffsetSurface( lRailIdx );
        if ( pOrigOffSurf != NULL ) {
            SER( pOrigOffSurf->Copy( m_crContext, pSurf ) );
            pOffSurf = SM_CAST_PTR( SmOffsetSurface, pSurf );
        }
        SER( pNewFilletGeom->SetOffsetSurface( lRailIdx, pOffSurf ) );
        SER( pFilletGeom->   SetOffsetSurface( lRailIdx, pAdjOffSurf ) );
    }
    else
    {
        // New FG becomes rollover, on adjacent face.
        pNewFilletGeom->SetFilletGeomType( SM_FG_TANGENT_ROLLOVER );
        pNewFilletGeom->GetRail(lRailIdx)->SetOriginalFace( pAdjFace );

        SER( pNewFilletGeom->SetOffsetSurface( lRailIdx, pAdjOffSurf ) );
    }

    m_vFilletGeoms.Add( pNewFilletGeom );


    // Grab the two new Vertices.  TopologySplit() leaves the new ones
    // in the middle, which is the high end of the original FG.

    SmFilletVertex * pNewV[2] = { NULL, NULL };
    for (ULONG kk=0; kk<2; kk++)
    {
        SmFilletEdge *pThisRail = pFilletGeom->GetRail(kk);
        SmVertex * pV = pThisRail->GetVertex();
        pNewV[kk] = (SmFilletVertex*)( pThisRail->GetOtherVertex(pV) );
    }

    // Set them as each others' mates.
    pNewV[0]->SetMate(0,pNewV[1]);
    pNewV[1]->SetMate(0,pNewV[0]);

    // The new vertex on the split rail is now the vertex on the
    // side edge, where the input vertex was.  Move the geometry
    // from the original split vertex into the new one.
    pNewV[lRailIdx]->SetPoint(            pFilVtx->GetPoint() );
#ifdef SM_USE_NEWTOL      
    SM_NEWTOL_LINE pNewV[lRailIdx]->ClearLocalZoneTol3d( ) ;  // what's going on with the MaxGaps here?
#else // SM_USE_OLDTOL 
    SM_OLDTOL_LINE pNewV[lRailIdx]->SetTolerance(        pFilVtx->GetTolerance(), FALSE );
#endif // SM_USE_OLDTOL
    pNewV[lRailIdx]->SetFilletVertexType( pFilVtx->GetFilletVertexType() );
    pNewV[lRailIdx]->SetStatus( SM_FIL_UNPROCESSED /* cbi pFilVtx->GetStatus() */  );

    // We have to get the proper SideEdgeuse into the point class.
    // pNewV[lRailIdx]->CopyPointClass( pFilVtx->GetPointClassification() );

    // Find the Filleted Edgeuse on the same face as the rail curve.
    if ( bClosedSideFace )
    {
        pSideEU = pSideEU->GetRadial();
    }
    else
    {
        SmEdgeuse *pFilletedEU = this->GetEdgeuse( lRailIdx );
        SmFace *pRailFace = pFilletedEU->GetFace();
        if ( pSideEU->GetFace() != pRailFace )
          { pSideEU = pSideEU->GetRadial(); }
        SM_ASSERT( pSideEU->GetFace() == pRailFace );
    }
    pNewV[lRailIdx]->SetPointClass( pFilVtx->GetPointClass(), pSideEU );

    // The opposite new vertex is now this new vertex's mate.
    pNewV[1-lRailIdx]->SetFilletVertexType( SM_FV_MATE );

    // The original split-vertex, and its mate, are no longer valid.
    pFilVtx->SetStatus( SM_FIL_UNPROCESSED );
    SmFilletVertex *pOppFV = NULL;
    if ( lWhichEnd == 0 )
    {
        pOppFV = pFilletGeom->GetRailVertex( 0, 1-lRailIdx );
    }
    else
    {
        pOppFV = pNewFilletGeom->GetRailVertex( 1, 1-lRailIdx );
    }
    pOppFV->SetStatus( SM_FIL_UNPROCESSED );

    return SM_SUCCESS;

} // end SmFilletSolver::SplitFilletGeomAtSideEdge


/*******************************************************************//**
Static Functions:
***********************************************************************/

/****************************************************************
PURPOSE: Trim the rail curves to the specified interval.

NOTES:
****************************************************************/
static SmStatus sm_TrimRail
  (SmFilletEdge * pRail,
   const SmExtent1d & crTrimInterval)
{
    SmFilletEdgeuse * pPrimEU = (SmFilletEdgeuse*)pRail->GetPrimaryEdgeuse();
    SmFilletEdgeuse * pMateEU = (SmFilletEdgeuse*)pPrimEU->GetMate();
    SmBSplineCurve *pCurve = SM_CAST_PTR(SmBSplineCurve,pRail->GetCurve()); NER(pCurve);
    SmBSplineCurve *pUVCurve = pPrimEU->GetUVTrimCurve(); NER(pUVCurve);
    SmBSplineCurve *pOrigUVCurve = pMateEU->GetUVTrimCurve(); NER(pOrigUVCurve);

    SmExtent1d sTrim = crTrimInterval;
    SER(pCurve->Trim(sTrim));        // may snap sIvl by tol to existing knots
    SER(pUVCurve->Trim(sTrim));      // may snap sIvl by tol to existing knots
    SER(pOrigUVCurve->Trim(sTrim));  // may snap sIvl by tol to existing knots

    return SM_SUCCESS;
} // end sm_TrimRail

/****************************************************************
PURPOSE:  Find the next curve interval that is not outside of face

NOTES:
****************************************************************/
static SmStatus sm_FindNextGoodInterval
  (SmCurveClassification * pCC,
   ULONG                   lStartSearch,
   SmBoolean             & rbFound,
   ULONG                 & rlFoundIndex)
{
  rbFound = FALSE;
  for (ULONG i=lStartSearch; i<pCC->GetSize(); i++) 
    {
      SmCurveInterval & rIvl1 = (*pCC)[i];
      if(   rIvl1.m_vMid.GetPointClass() == SM_PC_FACE
         || rIvl1.m_vMid.GetPointClass() == SM_PC_EDGE) 
        {
          rbFound = TRUE;
          rlFoundIndex = i;
          return SM_SUCCESS;
        }
    }
  return SM_SUCCESS;
} // end sm_FindNextGoodInterval

/*******************************************************************//**
    END - Static functions
***********************************************************************/


/*******************************************************************//**
PURPOSE: Split the Fillet Geom Object into one or more pieces.
    This routine will fix problems which may occur in the initial
    creation of rail points.

NOTES: It is assumed that the interval lists have corresponding
   segments and that the segments are disjoint.  They also need to be
   ordered by increasing parameter value of the rails.

   When done, there will be one FilletGeom for each input interval,
   and rails, etc, are trimmed to the intervals.

   This method was written for surface-surface filleting.  More work might
   be required to make it work for topology-based filleting.
***********************************************************************/
SmStatus SmFilletSolver::SplitFilletGeom
 (ULONG                         lFilletGeomIndex,
  SmBoundaryTrimmingType        eTrimType,
  const SmTArray<SmExtent1d>  & crIvls1,          // intervals for splitting rail[0]
  const SmTArray<SmExtent1d>  & crIvls2,          // intervals for splitting rail[1]
  const SmCurveClassification & crCC1,            // describes rail[0]
  const SmCurveClassification & crCC2)            // describes rail[1]
{
  SM_ASSERT(crIvls1.GetSize() == crIvls2.GetSize());

  SmFilletGeom *pFilletGeom = m_vFilletGeoms[lFilletGeomIndex];

  SmTArray<SmFilletGeom*> sFilletGeoms;
  for (ULONG ii=0; ii<crIvls1.GetSize(); ii++) 
    {
      // Create a new Fillet Geometry object by copying old one.
      // For now just take one of the intervals.  In the future we may be
      // smarter and take the one which is closest to the start point
      // or perhaps take them all.
      SmFilletGeom *pNewFilletGeom;
      if (ii < crIvls1.GetSize() - 1) 
        {
          pNewFilletGeom = new SmFilletGeom(pFilletGeom);
          m_vFilletGeoms.Add(pNewFilletGeom);
        }
      else 
        {
          pNewFilletGeom = pFilletGeom;
        }
      sFilletGeoms.Add(pNewFilletGeom);
    }

  for (ULONG i=0; i<crIvls1.GetSize(); i++) 
    {
      SmFilletGeom *pNewFilletGeom = sFilletGeoms[i];

      // Found classifications of each interval
      SmPointClassification sPC00(SmTol::GetSrcZoneTol3d(&crCC1), &m_crContext), 
                            sPC01(SmTol::GetSrcZoneTol3d(&crCC1), &m_crContext), 
                            sPC10(SmTol::GetSrcZoneTol3d(&crCC2), &m_crContext), 
                            sPC11(SmTol::GetSrcZoneTol3d(&crCC2), &m_crContext);

      SmExtent1d sIvl1 = crIvls1[i];
      SmExtent1d sIvl2 = crIvls2[i];

      // If we split things then set the processed flag for interior
      // vertexuses to FALSE.
      SmFilletVertex *pV00 = pFilletGeom->GetRailVertex(0,0);
      SmFilletVertex *pV01 = pFilletGeom->GetRailVertex(0,1);
      SmFilletVertex *pV10 = pFilletGeom->GetRailVertex(1,0);
      SmFilletVertex *pV11 = pFilletGeom->GetRailVertex(1,1);

      if (i > 0) 
        {
          pV00->SetStatus(SM_FIL_UNPROCESSED);
          pV10->SetStatus(SM_FIL_UNPROCESSED);
        }
      if (i<crIvls1.GetSize()-1) 
        {
          pV01->SetStatus(SM_FIL_UNPROCESSED);
          pV11->SetStatus(SM_FIL_UNPROCESSED);
        }

      SER(crCC1.FindPointClass(sIvl1.GetMin(),0.0,sPC00));
      SER(crCC1.FindPointClass(sIvl1.GetMax(),0.0,sPC01));
      SER(crCC2.FindPointClass(sIvl2.GetMin(),0.0,sPC10));
      SER(crCC2.FindPointClass(sIvl2.GetMax(),0.0,sPC11));

      // Create the fillet points using the rail and parameter value

      // Here we will have to do something if we are doing topology filleting
      // to distinguish between data that we have and data that we do not.

      SER(pNewFilletGeom->SetUpRailPoint(0,0,sIvl1.GetMin(),&sPC00));
      SER(pNewFilletGeom->SetUpRailPoint(0,1,sIvl1.GetMax(),&sPC01));
      SER(pNewFilletGeom->SetUpRailPoint(1,0,sIvl2.GetMin(),&sPC10));
      SER(pNewFilletGeom->SetUpRailPoint(1,1,sIvl2.GetMax(),&sPC11));

      // Refine fillet points using numerical solvers to get exact answers
      // at the points to minimize gaps in downstream processing.
//            SER(pNewFilletGeom->RefineFilletPoint(0,0));
//            SER(pNewFilletGeom->RefineFilletPoint(0,1));
//            SER(pNewFilletGeom->RefineFilletPoint(1,0));
//            SER(pNewFilletGeom->RefineFilletPoint(1,1));

      // Recreate the fillet surface using the refined fillet points
      // to get the fillet surface to go exactly through the points at
      // edges and vertices.
//            SER(pNewFilletGeom->ReCalcFilletGeom());

      // Create the ending curve on the fillet
      SmFilletBrep * pBrep = m_pExecutive->m_pPseudoBrep;
      SER(pNewFilletGeom->CreateEndCurve(0,eTrimType,pBrep));
      SER(pNewFilletGeom->CreateEndCurve(1,eTrimType,pBrep));

      // Trim the rail and that should give us a nice set of boundary
      // for the creation of topology.
      SER(sm_TrimRail(pNewFilletGeom->GetRail(0),sIvl1));
      SER(sm_TrimRail(pNewFilletGeom->GetRail(1),sIvl2));
    }

  return SM_SUCCESS;

} // end SmFilletSolver::SplitFilletGeom


/*******************************************************************//**
PURPOSE: Split and trim the fillet surface using the information
    contained in the classifications attached to the fillet geometry.

NOTES:
***********************************************************************/
SmStatus SmFilletSolver::SplitAndTrimFilletSurface
  (ULONG lFilletGeomIndex,
   SmBoundaryTrimmingType eTrimType)
{
    SmFilletGeom * pFilletGeom = m_vFilletGeoms[lFilletGeomIndex];
    SmCurveClassification *pCC1 = pFilletGeom->GetRail(0)->GetCurveClass(); NER(pCC1);
    SmCurveClassification *pCC2 = pFilletGeom->GetRail(1)->GetCurveClass(); NER(pCC2);

    ULONG lStartSearch1 = 0;
    ULONG lStartSearch2 = 0;
    SmBoolean bFound1, bFound2;

#ifdef SM_DEBUG_CODE
    if ( DebugLevel() > 0 ) {
        pCC1->Dump();
        pCC2->Dump();
    }
#endif

    // Find first interval to be used
    SER(sm_FindNextGoodInterval(pCC1,lStartSearch1,bFound1,lStartSearch1));
    SER(sm_FindNextGoodInterval(pCC2,lStartSearch2,bFound2,lStartSearch2));

    SmTArray<ULONG> sFound1, sFound2;

    while (bFound1 && bFound2) {
        SmCurveInterval & rCIvl1 = (*pCC1)[lStartSearch1];
        SmCurveInterval & rCIvl2 = (*pCC2)[lStartSearch2];
        SmExtent1d sIvl1, sIvl2;
        // If we have disjoint intervals there is not much to be
        // done except to get the next good interval.
        if (!rCIvl1.m_vInterval.AreDisjoint(rCIvl2.m_vInterval)) {
            // multiple segments.
            if (eTrimType == SM_BT_MINIMAL) {
                SER( rCIvl1.m_vInterval.Intersect(rCIvl2.m_vInterval,sIvl1) );
                SER( rCIvl1.m_vInterval.Intersect(rCIvl2.m_vInterval,sIvl2) );
            }
            else if (eTrimType == SM_BT_MAXIMAL) {
                rCIvl1.m_vInterval.Union(rCIvl2.m_vInterval,sIvl1);
                rCIvl1.m_vInterval.Union(rCIvl2.m_vInterval,sIvl2);
            }
//            else if (eTrimType == SM_BT_ROLLOVER) {
//                // For Roll over we take the intersection of the interval
//                // except for on the ends where we take the maximal values.
//                rCIvl1.m_vInterval.Intersect(rCIvl2.m_vInterval,sIvl1);
//                rCIvl1.m_vInterval.Intersect(rCIvl2.m_vInterval,sIvl2);
//                if (lStartSearch1 == 0) {
//                    sIvl1.SetMinMax(rCIvl1.m_vInterval.GetMin(),sIvl1.GetMax());
//                }
//                if (lStartSearch2 == 0) {
//                    sIvl2.SetMinMax(rCIvl2.m_vInterval.GetMin(),sIvl2.GetMax());
//                }
//                if (lStartSearch1 == pCC1->GetSize()-1) {
//                    sIvl1.SetMinMax(sIvl1.GetMin(),rCIvl1.m_vInterval.GetMax());
//                }
//                if (lStartSearch2 == pCC2->GetSize()-1) {
//                    sIvl2.SetMinMax(sIvl2.GetMin(),rCIvl2.m_vInterval.GetMax());
//                }
//            }
            else {
                sIvl1 = rCIvl1.m_vInterval;
                sIvl2 = rCIvl2.m_vInterval;
            }
            // Make sure that length of the intervals is sufficient.
            const SmCurve *pCurve1 = pCC1->GetCurve();
            const SmCurve *pCurve2 = pCC2->GetCurve();
            double dTol = GetThisApproxTol3d();
            SmBoolean bTooSmall = FALSE;
            if (pCurve1->ApproximateLength(sIvl1,5) < dTol) {
                bTooSmall = TRUE;
            }
            if (pCurve2->ApproximateLength(sIvl2,5) < dTol) {
                bTooSmall = TRUE;
            }

            if (!bTooSmall) {
                sFound1.Add(lStartSearch1);
                sFound2.Add(lStartSearch2);
            }
        }
        bFound1 = FALSE;
        bFound2 = FALSE;
        // Move to next good interval(s)
        if (rCIvl1.m_vInterval.GetMax() < rCIvl2.m_vInterval.GetMax() + SM_EFF_ZERO_SQ) {
            lStartSearch1++;
            bFound2 = TRUE;
            SER(sm_FindNextGoodInterval(pCC1,lStartSearch1,bFound1,lStartSearch1));
        }
        if (rCIvl2.m_vInterval.GetMax() < rCIvl1.m_vInterval.GetMax() + SM_EFF_ZERO_SQ) {
            bFound1 = TRUE;
            lStartSearch2++;
            SER(sm_FindNextGoodInterval(pCC2,lStartSearch2,bFound2,lStartSearch2));
        }
    } // While


    SmTArray<SmExtent1d> sIvls1, sIvls2;

    SmBoolean bFoundIntersection = FALSE;
    for (ULONG i=0; i<sFound1.GetSize(); i++) {
        SmCurveInterval & rCIvl1 = (*pCC1)[sFound1[i]];
        SmCurveInterval & rCIvl2 = (*pCC2)[sFound2[i]];
        SmExtent1d sIvl1, sIvl2;
        // If we have disjoint intervals there is not much to be
        // done except to get the next good interval.
        if (!rCIvl1.m_vInterval.AreDisjoint(rCIvl2.m_vInterval)) {
            bFoundIntersection = TRUE;
            if (eTrimType == SM_BT_MINIMAL) {
                SER( rCIvl1.m_vInterval.Intersect(rCIvl2.m_vInterval,sIvl1) );
                SER( rCIvl1.m_vInterval.Intersect(rCIvl2.m_vInterval,sIvl2) );
            }
            else if (eTrimType == SM_BT_MAXIMAL) {
                rCIvl1.m_vInterval.Union(rCIvl2.m_vInterval,sIvl1);
                rCIvl1.m_vInterval.Union(rCIvl2.m_vInterval,sIvl2);
            }
            else if (eTrimType == SM_BT_NONE) {
                if (i>0) break;
                sIvl1 = pCC1->GetInterval();
                sIvl2 = pCC2->GetInterval();
            }
            else {
                sIvl1 = rCIvl1.m_vInterval;
                sIvl2 = rCIvl2.m_vInterval;
            }

            sIvls1.Add(sIvl1);
            sIvls2.Add(sIvl2);
        }  // Intervals are not disjoint
    }


    if (!bFoundIntersection) {
        // No Good classifications.  Delete the fillet surface.
        SM_ASSERT(pFilletGeom->GetFilletSurface() != NULL) ; delete pFilletGeom->GetFilletSurface() ;
        pFilletGeom->m_pFilletSurface = NULL;
    }
    else {
        SER(SplitFilletGeom(lFilletGeomIndex,eTrimType,sIvls1,sIvls2,*pCC1,*pCC2));
    }

    return SM_SUCCESS;

} // end SmFilletSolver::SplitAndTrimFilletSurface

/****************************************************************
PURPOSE:  The following class is for solving the offset radius of
             SmConstantDistanceFS when clipping by surface-domain occurred
NOTES:
****************************************************************/
class SmFindClippedRadiusEFO : public SmEvalFunctionObject
{
protected:
  SmOffsetSurface * m_pSurface;             // fillet offset surface whose domain boundary is being crosses
  SmOffsetSurface * m_pOtherSurface;        // fillet other offset surface
  SmSurfParamType   m_eSurfParam;           // offset surface domain boundary isoParameter direction
  double            m_dIsoParam;            // offset surface domain boundary isoParameter value
  double            m_dOrientation;         // orientation for offset surface
  double            m_dOtherOrientation;    // orientation for other offset surface
  SmExtent2d        m_OtherDomain;          // domain for other offset surface
  double            m_dGuessT;              // u or v value on isoParameter line of actual intersection location
  double            m_dThisApproxTol3d;     // max distance between
  SmPoint2d       & m_rUV;                  // starts as guess m_pSurface      uvPoint for m_pSurface/m_pOtherSurface intersection - then refined
  SmPoint2d       & m_rOtherUV;             // starts as guess m_pOtherSurface uvPoint for m_pSurface/m_pOtherSurface intersection - then refined
public:
  double            m_dTargetDistance;      // desired distance for ConstantDistanceFillets - the distance between railCurves
public:

  // constructor
  SmFindClippedRadiusEFO(SmFilletSolver * pFS,           // in : container for target edge filleting data
                         ULONG            lSurfIndex,    // in : index of offsetSurface whose boundary is being crossed by offsetSurf/offsetSurf intersection curve
                         SmExtent2d       sUVDomains[2], // in : uvDomains for fillet surfaces
                         SmSurfParamType  eSurfParam,    // in : isoParameter direction of offsetSurface boundary being crossed
                         double           dIsoParam,     // in : isoParameter value     of offsetSurface boundary being crossed
                         SmVector2d       sUV[2]);       // in : guess uvPoints near actual offsetSurf/offsetSurf intersection curve
  // destructor
  virtual ~SmFindClippedRadiusEFO() { }

  // evaluate
  virtual SmStatus Evaluate(double      dT,                   // in : target T value to query
                            double    & rdFOfT,               // out: F(T)     = Function value for given dT value
                            double    & rdFPrimeOfT,          // out: dF(T)/dT = Function derivative value for given dT value
                            SmBoolean & rbFoundAnswer,        // out: TRUE = Function value is within tolerance of zero, FALSE=Not
                            SmBoolean   bSignalErrors=TRUE) ; // NotUsed: in : TRUE = signal errors, FALSE=errors anticipated, don't signal, default:[TRUE]

} ; // end class SmFindClippedRadiusEFO

/****************************************************************
PURPOSE:

NOTES:
****************************************************************/
SmFindClippedRadiusEFO::SmFindClippedRadiusEFO
 (SmFilletSolver * pFS,                // in : container for target edge filleting data
  ULONG            lSurfIndex,         // in : index of offsetSurface whose boundary is being crossed by offsetSurf/offsetSurf intersection curve
  SmExtent2d       sUVDomains[2],      // in : uvDomains for fillet surfaces
  SmSurfParamType  eSurfParam,         // in : isoParameter direction of offsetSurface boundary being crossed
  double           dIsoParam,          // in : isoParameter value     of offsetSurface boundary being crossed
  SmVector2d       sUV[2])             // in : guess uvPoints near actual offsetSurf/offsetSurf intersection curve
 : m_eSurfParam(eSurfParam),
   m_dIsoParam(dIsoParam),
   m_rUV(sUV[lSurfIndex]),
   m_rOtherUV(sUV[1-lSurfIndex])
{
  SmConstantDistanceFS * pFilletSolver = (SmConstantDistanceFS*)pFS;
  m_pSurface                           = pFilletSolver->GetSurface(lSurfIndex);
  m_pOtherSurface                      = pFilletSolver->GetSurface(1-lSurfIndex);
  m_dOrientation                       = pFilletSolver->GetOrientation(lSurfIndex);
  m_dOtherOrientation                  = pFilletSolver->GetOrientation(1-lSurfIndex);
  m_OtherDomain                        = sUVDomains[1-lSurfIndex];
  m_dTargetDistance                    = pFilletSolver->GetDistance();
  m_dThisApproxTol3d                   = pFilletSolver->GetThisApproxTol3d();
  m_dGuessT                            = sUV[lSurfIndex].x;
} // end SmFindClippedRadiusEFO::SmFindClippedRadiusEFO

/****************************************************************
PURPOSE: find F(T) and dF(T)/dT for given T value where

     T  = current best guess for radius of a ConstantDistance Fillet
   F(T) = DistanceBetweenRailPointsAtIntersectionPoint - ConstantFilletDistance

   for a fillet center-curve point located at the intersection of one of
   the fillet offsetSurface boundaries with the other fillet offsetsurface.

NOTES:

  sets
    m_rUV      = m_pSurface      uvPoint for OffsetSurfaceBoundaryCurve/OtherOffsetSurface intersection point
    m_rOtherUV = m_pOtherSurface uvPoint for SurfaceBoundaryCurve/OtherOffsetSurface intersection point
  returns SM_ERR when isoParameter Curve does not intersect OtherSurface

METHOD ---
  1. set offsetSurface distance = Radius
  2. find OffsetSurfaceBoundaryCurve/OtherOffsetSurface intersection point
  3. find baseSurface RailPoints
  4. compute F(Radius)          = DistanceBetweenRailPoints - m_dTargetDistance
  5. compute dF(Radius)/dRadius with finite difference as
       dF(Radius)/dRadius = (F(radius+StepSize) - F(radius)) / StepSize
       get F(radius+StepSize) by repeating steps 1 to 4 after incrementing radius value
****************************************************************/
SmStatus SmFindClippedRadiusEFO::Evaluate
 (double      dRadius,        // in : guess radius value that sets rail curve positions and the rail curve distance value
  double    & rdFOfT,         // out: F(dRadius) = DistanceBetweenRailPoints(dRadius) - m_dTargetDistance
  double    & rdFPrimeOfT,    // out: dF(T)/dT   = Function derivative value for given dRadius value
  SmBoolean & rbFoundAnswer,  // out: TRUE = (DistanceBetweenRailPoints - m_dTargetDistance) is less than tolerance
                              //      FALSE= Not TRUE
  SmBoolean   bSignalErrors)  // NotUsed: in : TRUE = signal errors, FALSE=errors anticipated, don't signal, default:[TRUE]
{
  SM_REF1(bSignalErrors) ;
  // init output
  rdFOfT = 0.0;
  rdFPrimeOfT = 0.0;
  rbFoundAnswer = FALSE;

  // set offsetSurface offset values
  m_pSurface->SetOffsetDistance(dRadius*m_dOrientation);
  m_pOtherSurface->SetOffsetDistance(dRadius*m_dOtherOrientation);

  // build surface boundary isoParameter curve
  // Note, don't clamp the isoParam to the natural domain of the surface.
  SmIsoCurve sIsoCurve( *m_pSurface, m_eSurfParam, m_dIsoParam, TRUE, FALSE );
  sIsoCurve.SetContext(NULL);

  // Find isoParameterCurve/OtherSurface intersection point
  double    dDeviation, dTFound;
  SmPoint2d sOtherUVFound;
  SmBoolean bFoundGoodPoint = FALSE;
  SER(m_pOtherSurface->LocalCurveIntersect(m_OtherDomain,
                                           sIsoCurve, sIsoCurve.GetNaturalInterval(),
                                           m_dThisApproxTol3d, m_rOtherUV, m_dGuessT,
                                           bFoundGoodPoint, sOtherUVFound, dTFound, dDeviation));
  // failure: no isoParameter/OtherSurface intersection
  if (!bFoundGoodPoint)
    {
      return SM_ERR;
    }

  // set outputs - uvPoint values for intersection point
  m_rOtherUV = sOtherUVFound;
  m_dGuessT  = dTFound;
  if (m_eSurfParam == SM_SP_U) { m_rUV.x = m_dIsoParam;
                                 m_rUV.y = dTFound;
                               }
  else                         { m_rUV.x = dTFound;
                                 m_rUV.y = m_dIsoParam;
                               }

  // evaluate fillet rail points for this intersection
  SmPoint3d sPnt,sOtherPnt;
  const SmSurface * pBaseSurface      = m_pSurface->GetBaseSurface();
  const SmSurface * pOtherBaseSurface = m_pOtherSurface->GetBaseSurface();

  SmBoolean bTmp0 = FALSE, bTmp1 = FALSE ; 
  SmBSplineSurface * pBSplineSurface0 = SM_CAST_PTR(SmBSplineSurface, pBaseSurface     ) ; // [B236]
  SmBSplineSurface * pBSplineSurface1 = SM_CAST_PTR(SmBSplineSurface, pOtherBaseSurface) ;
  SmTemporaryChangeValue<SmBoolean> sClean0(pBSplineSurface0 ? pBSplineSurface0->GetOutOfBoundsEnabled() : bTmp0, TRUE) ;
  SmTemporaryChangeValue<SmBoolean> sClean1(pBSplineSurface1 ? pBSplineSurface1->GetOutOfBoundsEnabled() : bTmp1, TRUE) ;

  SER(pBaseSurface->EvaluatePoint(m_rUV,sPnt));
  SER(pOtherBaseSurface->EvaluatePoint(m_rOtherUV,sOtherPnt));

  // set function value = DistanceBetweenRailPointsAtIntersection - m_dTargetDistance
  double dDistance = sPnt.DistanceBetween(sOtherPnt);
  rdFOfT = dDistance-m_dTargetDistance;

  // check for convergence
  if (smos_Fabs(rdFOfT) < SM_EFF_ZERO)
    {
      rbFoundAnswer = TRUE;
    }

  // Compute function derivative value by finite difference
  // Compute new railPoint positions for a slightly increased radius value
  double dStepSize  = SM_ZONE_TOL_3D;
  dRadius += dStepSize;

  // set offset distances
  m_pSurface->SetOffsetDistance(dRadius*m_dOrientation);
  m_pOtherSurface->SetOffsetDistance(dRadius*m_dOtherOrientation);

  // find intersection point for slightly increased offset distances
  SER(m_pOtherSurface->LocalCurveIntersect(m_OtherDomain, sIsoCurve,
                                           sIsoCurve.GetNaturalInterval(),
                                           m_dThisApproxTol3d, m_rOtherUV, m_dGuessT,
                                           bFoundGoodPoint, sOtherUVFound, dTFound, dDeviation));
  // failure: quit when intersection fails
  if (!bFoundGoodPoint)
    {
      SER(SM_ERR);
    }

  // save the deviation uvPoint answers
  SmPoint2d sUVFound;
  if (m_eSurfParam == SM_SP_U) { sUVFound.x = m_dIsoParam;
                                 sUVFound.y = dTFound;
                               }
  else                         { sUVFound.x = dTFound;
                                 sUVFound.y = m_dIsoParam;
                               }

  // get the deviation rail point positions
  SER(pBaseSurface->EvaluatePoint(sUVFound,sPnt));
  SER(pOtherBaseSurface->EvaluatePoint(sOtherUVFound,sOtherPnt));

  // get distance between rail points
  double dDistance1 = sPnt.DistanceBetween(sOtherPnt);

  // set dF(Radius)/dRadius = F(radius+deviation) - F(radius) / dStepSize ;
  rdFPrimeOfT = (dDistance1-m_dTargetDistance - rdFOfT)/dStepSize;

  return SM_SUCCESS;

} // end SmFindClippedRadiusEFO::Evaluate

/*******************************************************************//**
PURPOSE: Step to the next point in the fillet given a uvDomain step size
     and compute the FilletPoint values.

NOTES:
  set
    3dCurve Point values
    Surface Point values where baseSurf(uvPnt)   = point on rail
                               offsetSurf(uvPnt) = point on fillet center-curve
     rNextPoint.UVPos(0) = surf0 UVPoint
     rNextPoint.UVPos(1) = surf1 UVPoint
     rNextPoint.m_adUserDoubles[0] = dRadiusFound;  // contstant distance fillets only
***********************************************************************/
SmStatus SmFilletSolver::StepSolve
  (const SmTsectPnt & crCurrentPoint,    // in : Starting Point
   double & rdStep,                      // i/o: Size of uvDomain step shortened if necessary to stay inside UVDomains
   const SmExtent2d & crUVDomain1,       // in : surf1 UVDomain - may be larger then current Srf domain because SmOffsetSurfaces can be extended 
   const SmExtent2d & crUVDomain2,       // in : surf2 UVDomain - may be larger then current Srf domain because SmOffsetSurfaces can be extended
   SmBoolean & rbFoundGoodPoint,         // out: TRUE = found a good point
                                         //      FALSE= didn't
   SmBoolean & rbClipped,                // out: TRUE = Step was clipped to stay in UVDomains
   SmBoolean & rbBoundaryHit,            // out: TRUE = CurrentPoint is on boundary; no Step possible.
   SmTsectPnt & rNextPoint)              // out: next intersection point, tangent, and associated surface values
{
  // init output
  rbFoundGoodPoint = FALSE;
  rbClipped        = FALSE;
  rbBoundaryHit    = FALSE;

  // locals
  SmTsectPnt *pCurrentPoint = SM_CONST_CAST( SmTsectPnt*, &crCurrentPoint );

  SmExtent2d sUVDomains[2];
  sUVDomains[0] = crUVDomain1;
  sUVDomains[1] = crUVDomain2;

  // First pass - just use first derivative of parameter space
  // curves to get start point.
  double dStepSize    = rdStep;
  double dOldStepSize = dStepSize;
  ULONG lSurfIndex    = 0;

  // see if uvLineSegment steps outside of UVDomain1
  // get lineSeg parameter for which comes first
  //    1. lineSeg endPoint
  //    2. lineSeg/extentBoundary intersection point nearest to startPoint
  dStepSize = sUVDomains[0].ClipLine2d(
          pCurrentPoint->UVPos(0),
          pCurrentPoint->UVDeriv(0), dOldStepSize );

  // remember when stepSize is truncated by UVDomain1
  if (smos_Fabs(dStepSize) < smos_Fabs(dOldStepSize) - SM_EFF_ZERO*smos_Fabs(dOldStepSize)
      || dStepSize * dOldStepSize < 0 )  // Also check direction change
    {
      dOldStepSize = dStepSize;
      lSurfIndex   = 0;
      rbClipped    = TRUE;
    }

  // see if (possibly shortened) lineSegment steps outside of UVDomain2
  dStepSize = sUVDomains[1].ClipLine2d(
          pCurrentPoint->UVPos(1),
          pCurrentPoint->UVDeriv(1), dOldStepSize );

  // remember when stepSize is truncated by UVDomain2
  if (smos_Fabs(dStepSize) < smos_Fabs(dOldStepSize) - SM_EFF_ZERO*smos_Fabs(dOldStepSize)
      || dStepSize * dOldStepSize < 0 )  // Also check direction change
    {
      dOldStepSize = dStepSize;
      lSurfIndex   = 1;
      rbClipped    = TRUE;
    }

  // avoid tolerance sized variations in dStepSize by setting it equal to dOldStepSize
  dStepSize = dOldStepSize;

  // quit - when clipped stepsize is zero, or direction change
  if (SM_IS_ZERO(dStepSize) || dStepSize * rdStep < 0 )
    {
      rbFoundGoodPoint = FALSE;
      rbBoundaryHit    = TRUE;
      return SM_SUCCESS;
    }

  // compute next sUV guess points for each surface
  SmPoint2d sUV[2];
  sUV[0] = pCurrentPoint->UVPos(0) + dStepSize * pCurrentPoint->UVDeriv(0);
  sUV[1] = pCurrentPoint->UVPos(1) + dStepSize * pCurrentPoint->UVDeriv(1);

  // clamp points to UVDomains (shouldn't need this now that dStepSize is clipped)
  sUV[0] = sUVDomains[0].ClampPoint2d(sUV[0]);
  sUV[1] = sUVDomains[1].ClampPoint2d(sUV[1]);

  // when next step LineSegment was clipped
  //  - find the next step point as an OffsetSurfaceDomainBoundaryCurve/Offsetsurface intersection
  //  - update rdStep Size value
  if (rbClipped)
    {
      // Do a local OffsetSurfaceIsoCurve/OffsetSurface intersection
      SmSurfParamType eSurfParam;
      double dIsoParam;
      double dGuess;

      // set isoParameter values for when uvPoint is on U extent boundary
      if (   smos_Fabs(sUV[lSurfIndex].x-sUVDomains[lSurfIndex].GetMin().x) < SM_EFF_ZERO
          || smos_Fabs(sUV[lSurfIndex].x-sUVDomains[lSurfIndex].GetMax().x) < SM_EFF_ZERO)
        {
          // set isoParameter
          dIsoParam  = sUV[lSurfIndex].x;
          dGuess     = sUV[lSurfIndex].y;
          eSurfParam = SM_SP_U;
        }
      else // set isoParameter values for when uvPoint is on V extent boundary
        {
          dIsoParam  = sUV[lSurfIndex].y;
          dGuess     = sUV[lSurfIndex].x;
          eSurfParam = SM_SP_V;
        }

      // locals
      double dTFound;
      double dRadius;
      SER( SetupOffsetValues( *pCurrentPoint, &dRadius ));

      // branch on SolverType
      if ( GetSolverType() == SM_FS_CONST_DIST )
        {
          // This is an implicit radius function: you have to do the solve
          // to find the appropriate radius value.
          //
          // set up the upcoming NewtonRaphson iteration evaluation function.
          //   F(radius) = distanceBetweenRailCurvePointsAtIntersection - DesiredDistance
          //     radius is varied until resulting railCurve Points are the proper distance
          //     apart.
          SmFindClippedRadiusEFO sEFO(this,lSurfIndex,sUVDomains,
                                      eSurfParam,dIsoParam,sUV);
          SmExtent1d     sInterval(0.0,sEFO.m_dTargetDistance*10.0);
          SmLocalSolve1d sLS(sEFO,sInterval,FALSE);
          double dRadiusFound;

          // use NewtonRaphson iteration to find radius value so that the
          //   offsetSurface intersection on the boundary curve has associated
          //   railCurvePoints the proper distance apart from one another.
          if (sLS.SolveIt(dRadius,1.0e-12,rbFoundGoodPoint,dRadiusFound) != SM_SUCCESS)
            {
              rbFoundGoodPoint = FALSE;
            }

          // save found point
          if (rbFoundGoodPoint)
            {
              rNextPoint.UVPos(0) = sUV[0];
              rNextPoint.UVPos(1) = sUV[1];
              rNextPoint.m_adUserDoubles[0] = dRadiusFound;
            }
        } // end implicit radius branch
      else
        {
          // These are explicit radius functions: the radius value can be
          // calculated up front.
          // For constant radius, there's nothing to do, but for variable,
          // we have to calculate it.
          if (GetSolverType() == SM_FS_VARIABLE_RADIUS)
            {
              // Find the current radius and set up surface offsets.
              // We really just need the parameter value for the
              // radius function, then SetupOffsetValues() takes care of it.
              // All we have here is a surface UV point.  We evaluate that and
              // drop it to the fillet's defining curve to get the radius parameter.

              SmVariableRadiusFS *pVarRadFS = (SmVariableRadiusFS*) this;
              SmPoint3d sSrfPt;
              SER( GetSurface(lSurfIndex)->EvaluatePoint( sUV[lSurfIndex], sSrfPt ));

              const SmCurve *p3DFilletEdgeCurve = pVarRadFS->GetOriginalCurve();
              SmSolution sSData[4];
              SmSolutionArray sSolutions(4,sSData);
              SmExtent1d sIvl = pVarRadFS->GetFilletLawInterval();
              // This will usually be outside of this domain, so:
              sIvl.ExpandRelative( 2.0 );
              SER( p3DFilletEdgeCurve->GlobalPointSolve( sIvl, SM_SO_MINIMIZE, sSrfPt,
                       m_dThisApproxTol3d, NULL, NULL, SM_SR_SINGLE,
                       sSolutions ));

              if (sSolutions.GetSize() != 1) { SER(SM_ERR); }

              // Parameter gets stored in here:
              rNextPoint.m_adUserDoubles[0] = sSolutions[0].m_vStart[0];

              SER( SetupOffsetValues( rNextPoint, &dRadius ));
            }

          // Create iso-curve
          // Note, don't clamp the isoParam to the natural domain of the surface
          SmIsoCurve sIsoCurve( *GetSurface(lSurfIndex), eSurfParam, dIsoParam,
                  TRUE, FALSE );
          sIsoCurve.SetContext(NULL);

          // get OtherOffsetSurface/IsoParameterCurve intersection
          double dDeviation;
          SmPoint2d sUVFound;
          SER(GetSurface(1-lSurfIndex)->LocalCurveIntersect(sUVDomains[1-lSurfIndex],
              sIsoCurve,sIsoCurve.GetNaturalInterval(),
              m_dThisApproxTol3d,sUV[1-lSurfIndex],dGuess,
              rbFoundGoodPoint,sUVFound,dTFound,dDeviation));

          // when an intersection point was found
          if (rbFoundGoodPoint)
            {
              // store solution surface uvPoints
              if (eSurfParam == SM_SP_U) { sUV[lSurfIndex].y = dTFound;
                                         }
              else                       { sUV[lSurfIndex].x = dTFound;
                                         }
              sUV[1-lSurfIndex] = sUVFound;

              // set rNextPoint outputs
              rNextPoint.UVPos(0) = sUV[0];
              rNextPoint.UVPos(1) = sUV[1];

              // for constant radius assisted fillets
              //   - refine railcurve points to lie on common plane containing
              //       center-curve point and railcurve surface normals
              if (GetSolverType() == SM_FS_CONST_RADIUS_ASSISTED)
                {
                  // Move offset surfaces to base surfaces.  This has the effect
                  // of allowing out-of-bounds evaluations.
                  // [Note, this will also cause sPlaneOrig not to be offset.
                  //  cbi Dec 06]
                  SetSurfacesToZeroOffset();

                  // define a plane through fillet center-curve point
                  // and perpendicular to tangent direction
                  //    = crossProduct(railCurve surface Normals)
                  SmPoint3d  sPlaneOrig;
                  SER(sIsoCurve.EvaluatePoint(dTFound,sPlaneOrig));

                  // get railcurve point surface normals
                  SmVector3d  sNormal1, sNormal2;
                  SER(GetSurface(0)->EvaluateNormal(sUV[0],TRUE,TRUE,sNormal1));
                  SER(GetSurface(1)->EvaluateNormal(sUV[1],TRUE,TRUE,sNormal2));
                  SmVector3d sPlaneNormal = sNormal1*sNormal2;
                  SER(sPlaneNormal.Unitize());

                  // find two railcurve points on plane at radius distance from plane's origin
                  SmBoolean bFoundSolution;
                  SER(PointOnPlaneSolve(sPlaneOrig, sPlaneNormal,
                                        sUVDomains[0], sUVDomains[1],
                                        sUV[0], sUV[1], bFoundSolution, rNextPoint));

                  // restore the surface offsets
                  // [Shouldn't this happen *before* PointOnPlaneSolve() ? cbi Dec 06]
                  ReloadSurfaceOffsets();

                  if (!bFoundSolution) { SER(SM_ERR); }

                }  // end (GetSolverType() == SM_FS_CONST_RADIUS_ASSISTED) check
            } // end (rbFoundGoodPoint == TRUE) check
        } // end (GetSolverType() != SM_FS_CONST_DIST) branch

      // when a next point was found - surf uvPoints are set
      if (rbFoundGoodPoint)
        {
          // set stepSize based on point found and relative distance from the old point.
          SmVector2d sDiff     = pCurrentPoint->UVPos(0) - rNextPoint.UVPos(0);
          double     dLengDiff = sDiff.Length();
          SmPoint2d  sUV0      = pCurrentPoint->UVPos(0) + rdStep * pCurrentPoint->UVDeriv(0);
          SmVector2d sDiffOld  = pCurrentPoint->UVPos(0) - sUV0;
          double     dLengTan  = sDiffOld.Length();
          if (dLengTan < SM_EFF_ZERO) { dLengTan = 1.0; }
          rdStep = rdStep * (dLengDiff / dLengTan);
          return SM_SUCCESS;
        }

    } // end stepSize was clipped check

  // arrive here when step size was not clipped making the next guess
  //   Find Next SmTsectPnt as solution to offsetSurface/offsetSurface intersection near current guess.
  //   Current guess stored in sUV[2] was estimated as
  //        sUV[0] = pCurrentPoint->UVPos(0) + dStepSize * pCurrentPoint->UVDeriv(0);
  //        sUV[1] = pCurrentPoint->UVPos(1) + dStepSize * pCurrentPoint->UVDeriv(1);

  // Set up Solver
  ULONG lJacobianSize = GetJacobianSize();

  SmExtentNd sIntervals(lJacobianSize);

  sIntervals[0] = SmExtent1d(crUVDomain1.GetMin().x,crUVDomain1.GetMax().x);
  sIntervals[1] = SmExtent1d(crUVDomain1.GetMin().y,crUVDomain1.GetMax().y);
  sIntervals[2] = SmExtent1d(crUVDomain2.GetMin().x,crUVDomain2.GetMax().x);
  sIntervals[3] = SmExtent1d(crUVDomain2.GetMin().y,crUVDomain2.GetMax().y);

  // set all bboundaries as nonPeriodic
  // no need for periodicity because steps spanning boundaries have been processed
  SmBoolean alPerData[16];
  SmTArray<SmBoolean> sPeriodicities(16,alPerData);

  sPeriodicities.Add(FALSE);
  sPeriodicities.Add(FALSE);
  sPeriodicities.Add(FALSE);
  sPeriodicities.Add(FALSE);

  // set 1st guess
  double adGuessData[16];
  SmTArray<double> sGuessT(16,adGuessData);
  sGuessT.Add(sUV[0].x);
  sGuessT.Add(sUV[0].y);
  sGuessT.Add(sUV[1].x);
  sGuessT.Add(sUV[1].y);

  // set solver initialization values
  SmSurface *pSurface1 = GetSurface(0);  // fillet offset surface1
  SmSurface *pSurface2 = GetSurface(1);  // fillet offset surface2
  ULONG lDoSurf2Calcs  = 0;               // Don't do surface2 calculations
  ULONG lSurf1Index    = 0;               // UV values of Surface1
  ULONG lSurf2Index    = 2;               // UV values of Surface2
  ULONG lRailIndex     = 0;

#ifdef SM_DEBUG_CODE
  if ( DebugLevel() > 0 ) {
      if ( FALSE )
          { smgfx_Erase(); }
      smgfx_SetLook(3,6, 0,0,0); pCurrentPoint->CrvPos().Draw(); sm_GraphicsLoop();
      pCurrentPoint->CrvDeriv().Draw(&pCurrentPoint->CrvPos()); sm_GraphicsLoop();
      smgfx_SetLook(4,6, 0,1,0); pSurface1->DrawAt(pCurrentPoint->UVPos(0),0); sm_GraphicsLoop();
      if ( FALSE )
          { pSurface1->DrawUV(2,2); }
      smgfx_SetLook(5,6, 1,0,0); pSurface2->DrawAt(pCurrentPoint->UVPos(1),0); sm_GraphicsLoop();
      if ( FALSE )
          { pSurface2->DrawUV(2,2); }
      sm_GraphicsLoop();
  }
#endif

  // Note that in the future we could load the additional initial values
  // from the current point and also setup all of the intervals, periodicities,
  // and guess values outside of here and pass them in.  This may help
  // performance a little because we could eliminate LoadInitialValues
  // call each time we are here.

  // However for now we will just load the initial values again.
  // given current inputs - place a good guess for upcoming solve in sGuessT
  SER(LoadInitialValues(lDoSurf2Calcs,lRailIndex,
                        *pSurface1, lSurf1Index, pSurface2, lSurf2Index,
                        sIntervals, sPeriodicities, sGuessT));

  double sdData[16];
  SmTArray<double> sSolutionVector(16,sdData);

  // set upcoming NewtonRaphson solver evaluate function
  //  solution found when functions = zero
  //    1.
  SmVector3d sPlaneABC  = pCurrentPoint->CrvDeriv();
  SmPoint3d  sPlaneOrig = pCurrentPoint->CrvPos() + rdStep * sPlaneABC;
  double     dPlaneD    = - sPlaneOrig.Dot(sPlaneABC);
  SmFilletSphereSolveENFO sEvalFun(sPlaneABC,dPlaneD,*this);

  // allocated NewtonRaphson solver
  SmLocalSolveNd sLS(sEvalFun,          // in : Define Eqns to set to Zero (defines the DOF COUNT)
                     &sIntervals,       // in : NULL or sized:[N], N = this rEvalFun DOF count
                                        //      bounds on problem parameters.
                                        //      Set m_cpIntervals[i].SetUnbounded() for ivl whose prob params don't get clamped
                     &sPeriodicities) ; // in : NULL or sized:[N], N = this rEvalFun DOF count
                                        //      ivl[i] == TRUE = Interval[i] is periodic and values get wrapped
                                        //      ivl[i] == FALSE= Interval[i] is NOT periodic and values get clamped

  sLS.SetDesiredAccuracy(SM_EFF_ZERO_SQ);
  sLS.SetBoundaryHandler(SM_BH_TOTAL_BOUNDARY_HITS,8);
  sLS.SetMaximumIterations(30);

  // User NewtonRaphson to iterate to solution
  SmBoolean bFoundSolution;
  SER(sLS.SolveIt(sGuessT,m_dThisApproxTol3d,bFoundSolution,sSolutionVector));

  // Try it again and allow the plane to adjust as convergance happens
  if (!bFoundSolution)
    {
      sEvalFun.SetPlaneAdjusts(TRUE);
      SER(sLS.SolveIt(sGuessT,m_dThisApproxTol3d,bFoundSolution,sSolutionVector));
    }

  // when a solution was found - set output
  if (bFoundSolution)
    {
      rbFoundGoodPoint = TRUE;

      rNextPoint.UVPos(0).x = sSolutionVector[lSurf1Index];
      rNextPoint.UVPos(0).y = sSolutionVector[lSurf1Index+1];
      rNextPoint.UVPos(1).x = sSolutionVector[lSurf2Index];
      rNextPoint.UVPos(1).y = sSolutionVector[lSurf2Index+1];

      // Load other fillet surface variables being solved for into the
      // fillet point.  This should be sufficient to reproduce what is
      // happening at that point and to use as seed values for stepping.
      ULONG lCount = 0;
      for (ULONG i=lSurf2Index+2; i<sSolutionVector.GetSize(); i++)
        {
          if (lCount == SM_MAX_USER_DOUBLES) { SER(SM_ERR); }
          rNextPoint.m_adUserDoubles[lCount++] = sSolutionVector[i];
        }

    }

  return SM_SUCCESS;

} // end SmFilletSolver::StepSolve

/*******************************************************************//**
PURPOSE: Trim the rails using the processed fillet vertices.

NOTES:
***********************************************************************/
SmStatus SmFilletSolver::TrimRails()
{
  for (ULONG k=0; k<m_vFilletGeoms.GetSize(); k++)
    {
      SmFilletGeom * pFilletGeom = m_vFilletGeoms[k];
      SER(pFilletGeom->TrimRailCurves());
    }

  // all done
  return SM_SUCCESS;

} // end SmFilletSolver::TrimRails

/*******************************************************************//**
PURPOSE: Helper function for SmFilletSolver::TestRollOver. Attempts to 
    refine the CurveClassification in the presence of small intervals.

NOTES: This method was inspired by [B662].  In that case, a variable radius
    fillet creates a rail curve that crosses an edge->curve, so a RollOver is 
    necessary. But the rail curve doesn't follow the variable radius function 
    exactly. Based on the variable radius function, the rail curve should not 
    have crossed the edge->curve, but touched it tangentially. 

    Here, we find the tangential intersection and set the CurveClassification
    to match. ProcessRollOver will modify the fillet surfaces to match
    the CurveClassification
***********************************************************************/
SmStatus sm_RefineRailCurveClassification (
    ULONG lRailIndex,
    const SmTArray<ULONG> & crSmallIntervals,
    SmFilletEdge          & rRail,
    SmFilletGeom          & rFilletGeom,
    SmCurveClassification & rRailCurveClassification // i/o:
)
{
    // Locals
    SmCurve        * pRailCurve           = rRail.GetCurve();
    SmFace         * pRailFace            = rRail.GetOriginalFace();
    SmEdgeuse      * pSideEU              = NULL;
    SmFilletSolver * pFilletSolver        = rFilletGeom.GetFilletSolver();
    ULONG            lNumSmallIvls        = crSmallIntervals.GetSize();
    ULONG            lNumIvls             = rRailCurveClassification.GetSize();
    SmSolutionArray  sSolutions;
    SmZoneTol3d      sZoneTol3d;
    SmTsectPnt       sTsectPntMin, sTsectPntMax;
    SmPoint3d        s3DPntMin, s3DPntMax;
    SmBoolean        bFoundIntersection;
    SmStatus         eRtn = SM_SUCCESS;
    
    // The algorithm that follows is lifted from the ProcessRollOver algorithm.
    // Fixing the fillet geometry modification in situ is not practical, so we 
    // test here to determine if ProcessRollOver would fail in a preventable way.
    // If so, we squeeze out the small interval to prevent the failure.
    // Check each small interval for this repairable case.
    for ( ULONG ii = lNumSmallIvls; ii > 0; --ii )
    {
        // The first and last intervals should have been squeezed by TestRollOver
        if ( ii == 0 || ii == lNumIvls )
        {
            SM_ASSERT( ii != 0 && ii != lNumIvls );
            eRtn = SM_ERR;
            continue;
        }

        // Locals
        SmCurveInterval& rCI        = rRailCurveClassification[ii];
        SmExtent1d       sIvl       = rCI.GetInterval();
        SmEdge         * pStartEdge = rCI.m_vStart.GetEdgeObject();
        SmEdge         * pEndEdge   = rCI.m_vEnd.GetEdgeObject();
        SmEdge         * pEdge      = pStartEdge;  // Renamed to be clear there is only one Edge
        double           dEUParamMin, dEUParamMax;
        double           dGuessParam = 0.0;
        
        // Must have the interval classify to the same edge on both ends.
        if ( !pEdge || pStartEdge != pEndEdge )
        {
            eRtn = SM_ERR;
            continue;
        }

        // Get the Edgeuse connecting the RailFace and the Edge that XSects the Rail
        pSideEU = pEdge->GetEdgeuseOfFace( pRailFace );
        if ( !pSideEU )
        {
            eRtn = SM_ERR;
            continue;
        }

        // Get the RailCurve local end points 
        SER( pRailCurve->EvaluatePoint( sIvl.GetMin(), s3DPntMin) );
        SER( pRailCurve->EvaluatePoint( sIvl.GetMax(), s3DPntMax) );

        // Locals
        SmExtent1d   sEdgeIvl   = pEdge->GetInterval();
        SmCurve    * pEdgeCurve = pEdge->GetCurve();
        sZoneTol3d = SmTol::GetZoneTol3d( pEdge );

        // Get a guess for the Edge parameter where it XSects the rail on the Ivl min
        SE( pEdgeCurve->GlobalPointSolve( sEdgeIvl, SM_SO_MINIMIZE, s3DPntMin,
                                          SM_EFF_ZERO_SQRT, NULL, NULL, SM_SR_SINGLE, sSolutions ) );
        if ( sSolutions.GetSize() < 1 ) 
        {
            eRtn = SM_ERR;
            continue;
        }
        dGuessParam = sSolutions[0].m_vStart[0];

        // Calculate the Edge parameter where it will XSect the rail-to-be-created in ProcessRollOver for Ivl min
        SE( pFilletSolver->RailEdgeuseIntersect( lRailIndex, pSideEU, dGuessParam, bFoundIntersection, sTsectPntMin, dEUParamMin ) );

        // Get a guess for the Edge parameter where it XSects the rail on the Ivl max
        SE( pEdgeCurve->GlobalPointSolve( sEdgeIvl, SM_SO_MINIMIZE, s3DPntMax,
                                          SM_EFF_ZERO_SQRT, NULL, NULL, SM_SR_SINGLE, sSolutions ) );
        if ( sSolutions.GetSize() < 1 ) 
        {
            eRtn = SM_ERR;
            continue;
        }
        dGuessParam = sSolutions[0].m_vStart[0];

        // Calculate the Edge parameter where it will XSect the rail-to-be-created in ProcessRollOver Ivl max
        SE( pFilletSolver->RailEdgeuseIntersect( lRailIndex, pSideEU, dGuessParam, bFoundIntersection, sTsectPntMax, dEUParamMax ) );

        // Evaluate the params so we can check if the points are within tolerance
        SER( pEdgeCurve->EvaluatePoint( dEUParamMin, s3DPntMin ) );
        SER( pEdgeCurve->EvaluatePoint( dEUParamMax, s3DPntMax ) );

        // Can't squeeze the interval if it's not smaller than tolerance
        if ( s3DPntMin.DistanceBetween( s3DPntMax ) > sZoneTol3d )
        {
            eRtn = SM_ERR;
            continue;
        }

        // Set the intervals of the CurveIntervals on each side of the short interval
        double dNewParam = (dEUParamMax + dEUParamMin) / 2.;
        SmExtent1d sMaxSideIvl = rRailCurveClassification[ii + 1].GetInterval();
        SmExtent1d sMinSideIvl = rRailCurveClassification[ii - 1].GetInterval();
        SER( sMaxSideIvl.SetMin( dNewParam ) );
        SER( sMinSideIvl.SetMax( dNewParam ) );
        rRailCurveClassification[ii + 1].SetInterval( sMaxSideIvl );
        rRailCurveClassification[ii - 1].SetInterval( sMinSideIvl );

        // Remove the small interval
        rRailCurveClassification.RemoveAt( ii );
    }

    return eRtn;

} // end sm_RefineRailCurveClassification

/*******************************************************************//**
PURPOSE: Test for and process rollover fillets.

NOTES:
  Check to see whether either of the fillet's rail curves intersects
  any edges of the base Brep.  If so, then the fillet rolls over
  onto, or into, another face, or off into space.
  If the edge hit by a rail is a tangent edge, then we roll the fillet
  onto the next face.  This requires creating another fillet geom,
  splitting where the rail hits the edge.
  If it's not a tangent edge, then it has either the same or opposite
  convexity to the edge being filleted.  If the convexities are the same,
  then the fillet surface will crash into the adjoining face, and we'll
  trim it.  Otherwise, the fillet surface will be hanging over a cliff,
  and we'll extend the side face (or add a new little face) to intersect
  the fillet.
  With a non-tangent other edge, the same fillet geometry will still
  work, we only have to split the rail edge.
***********************************************************************/
SmStatus SmFilletSolver::TestRollOver
  (SmTArray<SmFilletCorner*> & rAdjustedCorners) // out:
{
  SmTArray<SmFilletGeom*> sFilletGeoms;
  sFilletGeoms.Append(m_vFilletGeoms);

  // for every filletSolver->FilletGeom
  for ( ULONG k=0; k<sFilletGeoms.GetSize(); k++ )
  {
      SmFilletGeom * pFilletGeom = sFilletGeoms[k];

      // for both rails
      for ( ULONG ii=0; ii<2; ii++ )
      {
          SmFilletEdge * pRail  = pFilletGeom->GetRail(ii);
          SmCurve      * pRailCurve = pRail->GetCurve();
          SmExtent1d     sIvl   = pRailCurve->GetNaturalInterval();

          // Classify the rail curve against the face in the base Brep,
          // i.e., see if it hits any edges.

          // create and add SmCurveClassification to RailEdge
          SmEdgeuse      * pMateEU            = pRail->GetPrimaryEdgeuse()->GetMate();
          SmBSplineCurve * pOrigFaceUVCurve   = pMateEU->GetUVTrimCurvePointer();
          SmCurveClassification * pCurveClass = new(m_crContext) SmCurveClassification
                                                    (pRailCurve,
                                                     sIvl,
                                                     pOrigFaceUVCurve,
                                                     GetThisApproxTol3d());
          pRail->SetCurveClass(pCurveClass);

          // Do the classification.
          SmFace * pFace = pRail->GetOriginalFace();
          SER( pFace->CurveOnClassify( TRUE,         // in : TRUE =Do expensive PtClassification for Curves not XSecting any Face->Bndrys
                                      *pCurveClass,  // i/o: contains curve to classify, accumulates intervals as they are found.
                                      TRUE)) ;       // in : TRUE=increase Edge tols for new verts found to be slightly too far from their surfaces

#ifdef SM_DEBUG_CODE
          // draw rail->curve(red), FilletSurface(black)
          if ( DebugLevel() > 0 )
            {
              sm_GraphicsLoop();
              smgfx_SetLook(2,2, 1,0,0); pRailCurve->DrawWDeriv(sIvl); sm_GraphicsLoop();
              smgfx_SetLook(1,1, 0,0,0);
              if ( pFilletGeom->GetFilletSurface() != NULL )
                  { pFilletGeom->GetFilletSurface()->DrawUV(0,0); sm_GraphicsLoop(); }
              sm_GraphicsLoop();
              pCurveClass->Dump();
            }
#endif // SM_DEBUG_CODE

          // SolidCraft change: the rail leaves its face only where it runs outside it. Where it
          // touches the face's boundary at a vertex or runs along a piece of an edge of it (a
          // blend whose setback reaches exactly to the rim of a hole), it is still on the face:
          // neighbouring intervals on the face or its boundary are one interval, so the rollover
          // below sees only the intervals that really leave the face.
          for (ULONG jj = pCurveClass->GetSize(); jj-- > 1; )
          {
              SmCurveInterval & rPrev = (*pCurveClass)[jj-1];
              SmCurveInterval & rThis = (*pCurveClass)[jj];
              if (   rPrev.m_vMid.GetPointClass() == SM_PC_UNKNOWN
                  || rThis.m_vMid.GetPointClass() == SM_PC_UNKNOWN )
                { continue; }
              if ( rPrev.m_vMid.GetPointClass() != SM_PC_FACE )
                { rPrev.m_vMid = rThis.m_vMid; }
              SER( pCurveClass->ConnectIntervals( jj, SM_IP_START, FALSE ));
          }

          // If the rail did hit an edge (and so is in more than one piece),
          // check for tiny intervals.  Squeeze out short intervals, if
          // they are the first or last one.  For interior short intervals,
          // don't squeeze, but give a warning.

          // What do we call 'short'? This is currently fairly arbitrary,
          // but 50*tol works better in practice than 10*tol.
          // [bd 060323a]
          double dMinIvlSize = 50.0*m_dThisApproxTol3d;

          SmTArray<ULONG> sSmallIntervals(4);
          int lNumIvls = pCurveClass->GetSize();
          if ( lNumIvls > 1 )
          {
              for (int jj=lNumIvls-1; jj>=0; jj--)
              {
                  // Check the distance between this interval's end points.
                  SmExtent1d sInterval = pCurveClass->GetSubInterval( jj );
                  SmPoint3d sPnt1, sPnt2;
                  SER( pRailCurve->EvaluatePoint( sInterval.GetMin(), sPnt1 ));
                  SER( pRailCurve->EvaluatePoint( sInterval.GetMax(), sPnt2 ));

                  double dDist = sPnt1.DistanceBetween( sPnt2 );
                  if ( dDist > dMinIvlSize )
                      continue;

                  // Squeeze out tiny end intervals, if first or last.
                  if ( jj == 0 )
                  {
                      SER(pCurveClass->ConnectIntervals( 0, SM_IP_START, TRUE ));
                  }
                  else if ( jj == lNumIvls-1 )
                  {
                      SER(pCurveClass->ConnectIntervals( jj, SM_IP_END, TRUE ));
                  }
                  else // assume tiny intervals are on the end
                  {
                      sSmallIntervals.Add( jj );
                  }
              } // end iter every interval
          }

          // Try to refine the CurveClassification.
          // The rail curve is an approximation that may not exactly match the radius curve
          // which can conflict with ProcessRollOver [B662]
          if ( SM_SUCCESS != sm_RefineRailCurveClassification( ii, sSmallIntervals, *pRail, *pFilletGeom, *pCurveClass ) )
          { SM_DBG_WARN( _T( "Intersecting a FilletRail with its Face->Boundaries produced a short internal interval - case not handled yet" ) ); }
      } // end iter both rails

      // Use the classifications to rebuild the surface if
      // rolling over tangent surfaces or onto opposite-convexity edges.

      SER( ProcessRollOver( pFilletGeom, rAdjustedCorners ));

  } // end iter every FilletSolver->FilletGeom

  return SM_SUCCESS;
} // end SmFilletSolver::TestRollOver

/*******************************************************************//**
PURPOSE: SolidCraft change. How far the trace extends a base surface's domain
    past crDomain, in each parameter direction.

NOTES: The extension factor times the domain's size (upstream), and in the
    directions that run across the edge (SetTraceReach) at least as far as
    m_dTraceReach in 3D: the size in parameters scaled by the domain's length
    in 3D along that direction, through its middle. Along the edge it keeps
    the factor's, or the trace would run on past the ends, where a variable
    round's law grows on and asks for more reach still. A direction of no
    length in 3D (a degenerate side) keeps the factor's.
***********************************************************************/
SmVector2d SmFilletSolver::GetDomainExtension
  (const SmSurface  * pSurf,      // in : the base surface
   const SmExtent2d & crDomain)   // in : its domain of interest
  const
{
  SmVector2d sSize = m_dSurfaceExtensionFactor * crDomain.GetSize();
  if ( pSurf == NULL || !( m_dTraceReach > 0.0 ))
    { return sSize; }
  SmPoint2d sMid = crDomain.GetMid();
  SmPoint3d asEnds[4];
  if (   pSurf->EvaluatePoint( SmPoint2d( crDomain.GetMin().x, sMid.y ), asEnds[0] ) != SM_SUCCESS
      || pSurf->EvaluatePoint( SmPoint2d( crDomain.GetMax().x, sMid.y ), asEnds[1] ) != SM_SUCCESS
      || pSurf->EvaluatePoint( SmPoint2d( sMid.x, crDomain.GetMin().y ), asEnds[2] ) != SM_SUCCESS
      || pSurf->EvaluatePoint( SmPoint2d( sMid.x, crDomain.GetMax().y ), asEnds[3] ) != SM_SUCCESS )
    { return sSize; }
  // which of the solver's base surfaces this is: its directions across the edge
  ULONG lSurface = 2;
  for (ULONG k=0; k<2; k++)
    {
      if ( m_pSurfaces[k] != NULL && m_pSurfaces[k]->GetBaseSurface() == pSurf )
        { lSurface = k; }
    }
  SmBoolean bAcrossU = ( lSurface < 2 ) ? m_abTraceAcross[lSurface][0] : TRUE;
  SmBoolean bAcrossV = ( lSurface < 2 ) ? m_abTraceAcross[lSurface][1] : TRUE;
  SmVector2d sUVSize = crDomain.GetSize();
  double d3DU = asEnds[0].DistanceBetween( asEnds[1] );
  double d3DV = asEnds[2].DistanceBetween( asEnds[3] );
  if ( bAcrossU && d3DU > SM_EFF_ZERO_SQRT )
    { sSize.x = smos_Max( sSize.x, m_dTraceReach * sUVSize.x / d3DU ); }
  if ( bAcrossV && d3DV > SM_EFF_ZERO_SQRT )
    { sSize.y = smos_Max( sSize.y, m_dTraceReach * sUVSize.y / d3DV ); }
  return sSize;
} // end SmFilletSolver::GetDomainExtension

/*******************************************************************//**
PURPOSE: SolidCraft change. Whether an edge is concave: the body's material
    wedge between its two faces spans more than a half turn there (a wall's
    foot, an inner corner), rather than less (a cliff's rim, an outer edge).

NOTES: SmEdgeuse::IsConvexRadialSector tells whether the sector between an
    edgeuse's face and its radial partner's spans more than a half turn; of
    the edge's loop edgeuses, one on the material side (its shell's region
    not void) is asked. A tangent edge spans a half turn and is not concave.
    An edge with no loop edgeuse on the material side is an error.
***********************************************************************/
SmStatus SmFilletSolver::IsConcaveEdge
  (const SmEdge & crEdge,      // in : the edge
   SmBoolean    & rbConcave)   // out: TRUE = the material wedge spans more than a half turn
{
  rbConcave = FALSE;
  SmTArray<SmEdgeuse*> sUses;
  SER( crEdge.GetEdgeuses( sUses ));
  for (ULONG i=0; i<sUses.GetSize(); i++)
    {
      SmEdgeuse * pUse = sUses[i];
      if ( pUse == NULL || !pUse->IsLoopEdgeuse() )
        { continue; }
      SmShell * pShell = pUse->GetShell();
      if ( pShell == NULL || pShell->GetRegion() == NULL || pShell->GetRegion()->IsVoid() )
        { continue; }
      rbConcave = pUse->IsConvexRadialSector( 10 );
      return SM_SUCCESS;
    }
  SER( SM_ERR );
  return SM_ERR;
} // end SmFilletSolver::IsConcaveEdge

/*******************************************************************//**
PURPOSE: SolidCraft change. Whether a rail that leaves the face it rolls on
    near crPoint runs on into the body's material: beyond a wall, the face's
    continuation lies inside the body.

NOTES: For RailsFitTheirFaces. Judged by the face's boundary edge nearest
    crPoint (the one the rail crosses, or for a rail nowhere on its face the
    nearest it lies beyond): concave (a wall's foot) or not (a cliff's rim).
***********************************************************************/
static SmStatus sm_InsideMaterial
  (const SmFace    & crFace,      // in : the face the rail rolls on
   const SmPoint3d & crPoint,     // in : a point of the rail beyond that face, or on its boundary
   SmBoolean       & rbInside)    // out: TRUE = the rail runs on into the material
{
  rbInside = FALSE;
  SmTArray<SmEdge*> sEdges;
  crFace.GetEdges( sEdges );
  const SmEdge * pNearest = NULL;
  double         dNearest = SM_BIG_DOUBLE;
  for (ULONG i=0; i<sEdges.GetSize(); i++)
    {
      const SmCurve * pCurve = sEdges[i]->GetCurve();
      if ( pCurve == NULL )
        { continue; }
      SmBoolean bDropped = FALSE;
      double    dParam = 0.0, dDist = SM_BIG_DOUBLE;
      SER( pCurve->DropPoint( sEdges[i]->GetInterval(), crPoint, NULL, SM_BIG_DOUBLE, NULL, bDropped, dParam, dDist ));
      if ( bDropped && dDist < dNearest )
        {
          dNearest = dDist;
          pNearest = sEdges[i];
        }
    }
  NER(pNearest);
  SER( SmFilletSolver::IsConcaveEdge( *pNearest, rbInside ));
  return SM_SUCCESS;
} // end sm_InsideMaterial

/*******************************************************************//**
PURPOSE: SolidCraft change. Whether pEdge is a knife edge: its two faces
    tangent along it and folded back on each other, the body's material
    between them of zero angle (a hole whose wall touches the bottom face
    along a line).
***********************************************************************/
SmBoolean SmFilletSolver::IsKnifeEdge
  (SmEdge * pEdge,             // in : the edge, of two faces
   double   dTangencyTolRad)   // in : faces this close to tangent are tangent
{
  if ( pEdge == NULL || pEdge->GetPrimaryEdgeuse() == NULL || pEdge->GetPrimaryEdgeuse()->GetRadial() == NULL )
    { return FALSE; }
  // the sector between the faces, a whole turn or none (a smooth edge's is a half turn)
  double dAngle = SM_PI;
  if ( pEdge->GetPrimaryEdgeuse()->CalcSectorAngle( pEdge->GetInterval().Evaluate( 0.5 ), dAngle ) != SM_SUCCESS )
    { return FALSE; }
  return ( smos_Min( dAngle, SM_2PI - dAngle ) < dTangencyTolRad ) ? TRUE : FALSE;
} // end SmFilletSolver::IsKnifeEdge

/*******************************************************************//**
PURPOSE: SolidCraft change. Whether the fillet's end section at each Nx1
    end lies on the end face: rbFit is FALSE when it runs off it.

NOTES: At a vertex of three edges where only this edge is filleted (an
    Nx1 corner) the fillet ends on the third face E, the end face: its end
    section is where the fillet's surface meets E's, from rail 0's crossing
    of its side edge (the edge of its face and E at the vertex) to rail 1's
    (SmFilletNx1Corner::CalcCornerEdgeGeom builds it as the end rib). The
    end face must hold all of it. Where E is narrower than the fillet and
    something stands on it next to the vertex (a lip on a shelled box's rim,
    the section running under the lip's foot: off E across a concave edge of
    it), part of the section lies beyond E's boundary under material: the
    corner's end rib leaves its face and the fillet comes back unsewn. The
    fillet then overflows its end face, and is built as an overflow is (the
    executive's m_bCliffOverflow): cut off by E's surface continued, and the
    material it removes subtracted. (A section off E across a convex edge,
    a hole's rim or a cliff, is not judged here.) The section
    is the surfaces' exact intersection between the two crossings
    (sm_SrfSrfIntersection), classified against E (SmFace::CurveOnClassify).
    Judged where the corner takes the end face as it is: both side edges
    of the filleted edge's convexity, neither tangent. Where one is of the
    other convexity the corner continues the end face past its boundary
    (a chamfer ending on a boss's wall, a concave chamfer at a missing
    quadrant), and an end that is no Nx1 corner, or where a rail misses its
    side edge, has no end section of this kind.
***********************************************************************/
SmStatus SmFilletSolver::EndSectionsFitTheirFaces
  (SmBoolean & rbFit)   // out: FALSE = an end section runs off its end face
{
  rbFit = TRUE;
  SmEdgeuse * pFilEU   = GetEdgeuse(0);       NER(pFilEU);
  SmEdge    * pFilEdge = pFilEU->GetEdge();   NER(pFilEdge);
  SmBoolean   bConcave = FALSE;
  SER( IsConcaveEdge( *pFilEdge, bConcave ));
  for (ULONG k=0; k<m_vFilletGeoms.GetSize() && rbFit; k++)
    {
      SmFilletGeom * pFG = m_vFilletGeoms[k];
      if (pFG == NULL || pFG->GetFilletGeomType() != SM_FG_DEFAULT)
        { continue; }
      SmSurface * pFilletSurface = pFG->GetFilletSurface();   NER(pFilletSurface);
      for (SmVertex * pV : { pFilEdge->GetStartVertex(), pFilEdge->GetEndVertex() })
        {
          NER(pV);
          SmFilletCorner * pCorner = GetFilletExecutive()->GetFilletCornerOfVertex( pV );
          SmTArray<SmEdge*> sVertexEdges;
          pV->GetEdges( sVertexEdges );
          if ( pCorner == NULL || pCorner->GetCornerType() != SM_FCR_N_x_1 || sVertexEdges.GetSize() != 3 )
            { continue; }
          // each rail's side edge at pV, the end face they share, and where the rail crosses it
          SmFace    * pEndFace = NULL;
          SmPoint3d   asCross[2];
          SmBoolean   bBoth = TRUE;
          for (ULONG ii=0; ii<2 && bBoth; ii++)
            {
              SmFilletEdge * pRail  = pFG->GetRail(ii);              NER(pRail);
              SmCurve      * pCurve = pRail->GetCurve();             NER(pCurve);
              SmFace       * pFace  = pRail->GetOriginalFace();      NER(pFace);
              SmEdgeuse    * pMate  = pRail->GetPrimaryEdgeuse()->GetMate();   NER(pMate);
              SmEdge       * pSide  = NULL;
              for (ULONG e=0; e<sVertexEdges.GetSize(); e++)
                {
                  SmEdge * pE = sVertexEdges[e];
                  SmTArray<SmFace*> sFaces;
                  pE->GetFaces( sFaces );
                  if ( pE != pFilEdge && sFaces.GetSize() == 2 && ( sFaces[0] == pFace || sFaces[1] == pFace ))
                    {
                      pSide = pE;
                      SmFace * pOther = ( sFaces[0] == pFace ) ? sFaces[1] : sFaces[0];
                      if ( pEndFace != NULL && pEndFace != pOther )
                        { pSide = NULL; }
                      pEndFace = pOther;
                    }
                }
              if ( pSide == NULL || pEndFace == NULL )
                { bBoth = FALSE; continue; }
              // An end face the corner takes as it is: a side edge of the filleted edge's
              // convexity, and not tangent (where one is of the other convexity the corner
              // continues the end face past its boundary, SmFilletNx1Corner::MakeCornerTopology)
              SmBoolean bSideConcave = FALSE, bOpposite = FALSE;
              SER( IsConcaveEdge( *pSide, bSideConcave ));
              if ( bSideConcave != bConcave || sm_TestManifoldTangentRollover( pSide, bOpposite, m_dTangencyTolerance, 0 ))
                { bBoth = FALSE; continue; }
              // the rail's crossing of pSide nearest pV
              SmCurveClassification sClass(pCurve, pCurve->GetNaturalInterval(), pMate->GetUVTrimCurvePointer(), GetThisApproxTol3d());
              SER( pFace->CurveOnClassify( TRUE, sClass, TRUE ));
              double dNearest = SM_BIG_DOUBLE;
              for (ULONG j=0; j<sClass.GetSize(); j++)
                {
                  for (const SmPointClassification * pPC : { &sClass[j].m_vStart, &sClass[j].m_vEnd })
                    {
                      const double dT = ( pPC == &sClass[j].m_vStart ) ? sClass[j].GetMin() : sClass[j].GetMax();
                      if ( pPC->GetPointClass() != SM_PC_EDGE || pPC->GetObject() != pSide )
                        { continue; }
                      SmPoint3d sAt;
                      SER( pCurve->EvaluatePoint( dT, sAt ));
                      const double dDist = sAt.DistanceBetween( pV->GetPoint() );
                      if ( dDist < dNearest )
                        {
                          dNearest = dDist;
                          asCross[ii] = sAt;
                        }
                    }
                }
              if ( dNearest == SM_BIG_DOUBLE )
                { bBoth = FALSE; }
            }
          if ( !bBoth || pEndFace == NULL || pEndFace->GetSurface() == NULL )
            { continue; }
          // the end section, and whether the end face holds it
          const SmContext & crContext = GetCreationContext();
          SmFilletVertex sFrom, sTo;
          sFrom.SetPoint( asCross[0] );
          sTo.SetPoint( asCross[1] );
          SmBSplineCurve * p3D = NULL, * pUVFillet = NULL, * pUVEnd = NULL;
          SmStatus eStat = sm_SrfSrfIntersection( crContext, &sFrom, &sTo, NULL, NULL, pFilletSurface, pEndFace->GetSurface(),
                                                  GetThisApproxTol3d(), GetThisAngTolRad(), p3D, pUVFillet, pUVEnd );
          SmObjDelete sDelete3D(p3D), sDeleteUV1(pUVFillet), sDeleteUV2(pUVEnd);
          SER( eStat );
          NER(p3D);
          // the section from one crossing to the other (the intersection may run on past them)
          double adEnd[2] = { 0.0, 0.0 };
          for (ULONG ii=0; ii<2; ii++)
            {
              SmBoolean bDropped = FALSE;
              double    dDist = SM_BIG_DOUBLE;
              SER( p3D->DropPoint( p3D->GetNaturalInterval(), asCross[ii], NULL, SM_BIG_DOUBLE, NULL, bDropped, adEnd[ii], dDist ));
              if ( !bDropped || dDist > 10.0 * GetThisApproxTol3d() )
                { SER( SM_ERR ); }
            }
          SmExtent1d sBetween;
          SER( sBetween.SetMinMax( smos_Min( adEnd[0], adEnd[1] ), smos_Max( adEnd[0], adEnd[1] )));
          SmCurveClassification sSection(p3D, sBetween, pUVEnd, GetThisApproxTol3d());
          SER( pEndFace->CurveOnClassify( TRUE, sSection, TRUE ));
          for (ULONG j=0; j<sSection.GetSize() && rbFit; j++)
            {
              if ( sSection[j].m_vMid.GetPointClass() != SM_PC_UNKNOWN )
                { continue; }
              // The section's ends lie on the side edges: a first or last run off the face no
              // longer than the classification's tolerance is its end on that edge.
              if ( j == 0 || j + 1 == sSection.GetSize() )
                {
                  SmPoint3d sA, sB;
                  SER( p3D->EvaluatePoint( sSection[j].GetMin(), sA ));
                  SER( p3D->EvaluatePoint( sSection[j].GetMax(), sB ));
                  if ( sA.DistanceBetween( sB ) <= GetThisApproxTol3d() )
                    { continue; }
                }
              // Off the face across a wall's foot (a concave edge of the end face: something
              // standing on it), the section runs under material. Across a convex edge (a hole's
              // rim, a cliff) it runs into the open, which this does not judge.
              SmBoolean bUnder = FALSE;
              for (const SmPointClassification * pPC : { &sSection[j].m_vStart, &sSection[j].m_vEnd })
                {
                  if ( pPC->GetPointClass() != SM_PC_EDGE || pPC->GetObject() == NULL )
                    { continue; }
                  SmBoolean bWallFoot = FALSE;
                  SER( IsConcaveEdge( *(const SmEdge*)pPC->GetObject(), bWallFoot ));
                  if ( bWallFoot )
                    { bUnder = TRUE; }
                }
              if ( bUnder )
                { rbFit = FALSE; }
            }
        }
    }
  return SM_SUCCESS;
} // end SmFilletSolver::EndSectionsFitTheirFaces

/*******************************************************************//**
PURPOSE: SolidCraft change. Whether every rail of this solver's fillet stays on the face it
    rolls on, so that the blend fits.

NOTES: Each rail, untrimmed (CalcFilletGeom made it), is classified against its face with
    SmFace::CurveOnClassify. It may run off the face past the filleted edge's ends
    (LeavesPastAnEnd; the corner trims it there, also where the face goes on beyond what the
    rail passes there and the rail comes back onto it across the same end's side edge), and
    off across a
    tangent side edge (a
    rollover onto the next face, which TestRollOver processes; rbRollsOver is TRUE then). A
    rail that is nowhere on its face, or that leaves it between its ends across an edge that
    is not tangent (a cliff: an opening's rim, the body's outline; or a wall), means the
    blend is larger than the faces: rbFit is FALSE. rbIntoMaterial tells which: TRUE when such
    a rail runs on inside the body's material (beyond a wall, the face's continuation lies
    inside the body): the edge it leaves across is concave (IsConcaveEdge), or for a rail
    nowhere on its face the face's edge nearest its point nearest the edge's middle.
    Touching or running along the face's boundary is on the face.
***********************************************************************/
SmStatus SmFilletSolver::RailsFitTheirFaces
  (SmBoolean & rbFit,           // out: FALSE = a rail leaves its face
   SmBoolean & rbRollsOver,     // out: TRUE  = a rail rolls over a tangent edge onto the next face
   SmBoolean & rbIntoMaterial)  // out: TRUE  = a rail that leaves its face runs on inside the body
{
  rbFit = TRUE;
  rbRollsOver = FALSE;
  rbIntoMaterial = FALSE;
  for (ULONG k=0; k<m_vFilletGeoms.GetSize(); k++)
    {
      SmFilletGeom * pFG = m_vFilletGeoms[k];
      if (pFG == NULL || pFG->GetFilletGeomType() != SM_FG_DEFAULT)
        { continue; }
      for (ULONG ii=0; ii<2; ii++)
        {
          SmFilletEdge * pRail = pFG->GetRail(ii);                 NER(pRail);
          SmCurve      * pCurve = pRail->GetCurve();               NER(pCurve);
          SmFace       * pFace = pRail->GetOriginalFace();         NER(pFace);
          SmEdgeuse    * pEU = pRail->GetPrimaryEdgeuse();         NER(pEU);
          SmEdgeuse    * pMate = pEU->GetMate();                   NER(pMate);
          SmCurveClassification sClass(pCurve, pCurve->GetNaturalInterval(), pMate->GetUVTrimCurvePointer(), GetThisApproxTol3d());
          SER( pFace->CurveOnClassify( TRUE, sClass, TRUE ));
          ULONG lNumIvls = sClass.GetSize();
          SmBoolean bOnFace = FALSE;
          for (ULONG j=0; j<lNumIvls; j++)
            {
              SmCurveInterval & rIvl = sClass[j];
              if (rIvl.m_vMid.GetPointClass() != SM_PC_UNKNOWN)
                {
                  bOnFace = TRUE;
                  continue;
                }
              // A run of intervals off the face (split where the rail passes through a point of
              // the face's boundary from outside to outside) leaves the face once: judge it at
              // its first interval, from where it leaves to the rail's end or start.
              if (j > 0 && sClass[j-1].m_vMid.GetPointClass() == SM_PC_UNKNOWN)
                { continue; }
              ULONG jLast = j;
              while (jLast+1 < lNumIvls && sClass[jLast+1].m_vMid.GetPointClass() == SM_PC_UNKNOWN)
                { jLast++; }
              // Where the rail leaves the face: the first interval ends there, the others start.
              const SmPointClassification & rLeave = (j == 0) ? sClass[jLast].m_vEnd : rIvl.m_vStart;
              const double dLeave = (j == 0) ? sClass[jLast].GetMax() : rIvl.GetMin();
              if (j == 0 || jLast+1 == lNumIvls)
                {
                  // The rail beyond the face from there to its end: past an end of the
                  // filleted edge, where a corner trims it, only if it leaves the face there.
                  SmBoolean bPastAnEnd = FALSE;
                  SER( LeavesPastAnEnd( *pFace, *pCurve, dLeave, rLeave, bPastAnEnd ));
                  if (bPastAnEnd)
                    { continue; }
                }
              else
                {
                  // A run between two on the face: the face goes on beyond something the rail
                  // passes (a boss the filleted edge ends on, rising from the face) and the rail
                  // comes back onto it. Past an end when the rail leaves the face across that
                  // end's side edge and comes back across it, wholly beyond the end, where the
                  // corner trims it. (Past an end by its foot alone it is not: any point of a run
                  // beyond the edge's extent has its foot at an end, wherever it crosses the
                  // face's boundary.)
                  SmBoolean bLeavesPast = FALSE, bReturnsPast = FALSE, bLeavesBySide = FALSE, bReturnsBySide = FALSE;
                  const SmVertex * pLeaveEnd = NULL, * pReturnEnd = NULL;
                  SER( LeavesPastAnEnd( *pFace, *pCurve, dLeave, rLeave, bLeavesPast, &pLeaveEnd, &bLeavesBySide ));
                  SER( LeavesPastAnEnd( *pFace, *pCurve, sClass[jLast].GetMax(), sClass[jLast].m_vEnd, bReturnsPast, &pReturnEnd, &bReturnsBySide ));
                  if (bLeavesPast && bReturnsPast && bLeavesBySide && bReturnsBySide && pLeaveEnd == pReturnEnd)
                    { continue; }
                }
              SmEdge * pSide = (rLeave.GetPointClass() == SM_PC_EDGE) ? (SmEdge*)rLeave.GetObject() : NULL;
              SmBoolean bOppositeNormals = FALSE;
              if (pSide != NULL && sm_TestManifoldTangentRollover(pSide, bOppositeNormals, m_dTangencyTolerance, 0))
                {
                  rbRollsOver = TRUE;
                  continue;
                }
              rbFit = FALSE;
              // the edge it leaves across (or the face's edge nearest where it leaves)
              SmBoolean bInside = FALSE;
              if ( pSide != NULL )
                { SER( IsConcaveEdge( *pSide, bInside )); }
              else
                {
                  SmPoint3d sLeave;
                  SER( pCurve->EvaluatePoint( dLeave, sLeave ));
                  SER( sm_InsideMaterial( *pFace, sLeave, bInside ));
                }
              if (bInside)
                { rbIntoMaterial = TRUE; }
            }
          if (!bOnFace)
            {
              rbFit = FALSE;
              // the rail's point nearest the middle of the filleted edge
              SmEdge * pFilEdge = GetEdgeuse(0)->GetEdge();                   NER(pFilEdge);
              SmPoint3d sMid, sOnRail;
              SER( pFilEdge->GetCurve()->EvaluatePoint( pFilEdge->GetInterval().Evaluate( 0.5 ), sMid ));
              SmBoolean bDropped = FALSE;
              double dParam = 0.0, dDist = 0.0;
              SER( pCurve->DropPoint( pCurve->GetNaturalInterval(), sMid, NULL, SM_BIG_DOUBLE, NULL, bDropped, dParam, dDist ));
              if (!bDropped)
                { SER( SM_ERR ); }
              SER( pCurve->EvaluatePoint( dParam, sOnRail ));
              SmBoolean bInside = FALSE;
              SER( sm_InsideMaterial( *pFace, sOnRail, bInside ));
              if (bInside)
                { rbIntoMaterial = TRUE; }
            }
        }
    }
  return SM_SUCCESS;
} // end SmFilletSolver::RailsFitTheirFaces

/*******************************************************************//**
PURPOSE: SolidCraft change. A face's jet at a contact: its point, first to
    third derivatives, its unit normal toward the ball's centre and the
    normal's derivatives, and its fundamental forms there.

NOTES: The normal n = s * m / |m|, m = Su x Sv, s = +-1 toward the centre.
    Its derivatives follow from n |m| = m: n_a = (m_a - n l_a) / l and
    n_ab = (m_ab - n_a l_b - n_b l_a - n l_ab) / l, l = |m|, l_a = <n, m_a>,
    l_ab = <n_b, m_a> + <n, m_ab> (m's derivatives need S's third). E, F, G
    and L, M, N (L = <Suu, n>) are the first and second fundamental forms
    with n toward the ball. rbJudged is FALSE at a pole (Su x Sv = 0).
***********************************************************************/
namespace
{
  struct sm_FaceJet
  {
    SmVector3d s, su, sv, suu, suv, svv;
    SmVector3d n, nu, nv, nuu, nuv, nvv;
    double     E = 0.0, F = 0.0, G = 0.0, L = 0.0, M = 0.0, N = 0.0;
  };
}

static SmStatus sm_FaceJetAt
  (const SmSurface  & crSurface,  // in : the face's surface
   const SmPoint2d  & crUV,       // in : the contact
   const SmVector3d & crCentre,   // in : the ball's centre
   sm_FaceJet       & rJet,       // out: the jet
   SmBoolean        & rbJudged)   // out: FALSE at a pole
{
  rbJudged = FALSE;
  // Evaluate's 1d order with three derivatives each way: index u*4 + v
  SmVector3d d[16];
  SER(crSurface.Evaluate(crUV, 3, 3, TRUE, TRUE, TRUE, d, FALSE));
  rJet.s = d[0];  rJet.sv = d[1];  rJet.svv = d[2];
  rJet.su = d[4]; rJet.suv = d[5]; rJet.suu = d[8];
  const SmVector3d & sVvv = d[3];
  const SmVector3d & sUvv = d[6];
  const SmVector3d & sUuv = d[9];
  const SmVector3d & sUuu = d[12];

  const SmVector3d sM = rJet.su * rJet.sv;
  const double dL = sM.Length();
  if (!(dL > SM_EFF_ZERO * rJet.su.Length() * rJet.sv.Length()) || !(dL > 0.0))
    { return SM_SUCCESS; }
  const SmVector3d sMu  = rJet.suu * rJet.sv + rJet.su * rJet.suv;
  const SmVector3d sMv  = rJet.suv * rJet.sv + rJet.su * rJet.svv;
  const SmVector3d sMuu = sUuu * rJet.sv + 2.0 * (rJet.suu * rJet.suv) + rJet.su * sUuv;
  const SmVector3d sMuv = sUuv * rJet.sv + rJet.suu * rJet.svv + rJet.su * sUvv;
  const SmVector3d sMvv = sUvv * rJet.sv + 2.0 * (rJet.suv * rJet.svv) + rJet.su * sVvv;

  const SmVector3d sN = sM / dL;
  const double dLu = sN.Dot(sMu), dLv = sN.Dot(sMv);
  const SmVector3d sNu = (sMu - sN * dLu) / dL;
  const SmVector3d sNv = (sMv - sN * dLv) / dL;
  const double dLuu = sNu.Dot(sMu) + sN.Dot(sMuu);
  const double dLuv = sNv.Dot(sMu) + sN.Dot(sMuv);
  const double dLvv = sNv.Dot(sMv) + sN.Dot(sMvv);
  const SmVector3d sNuu = (sMuu - sNu * (2.0 * dLu) - sN * dLuu) / dL;
  const SmVector3d sNuv = (sMuv - sNu * dLv - sNv * dLu - sN * dLuv) / dL;
  const SmVector3d sNvv = (sMvv - sNv * (2.0 * dLv) - sN * dLvv) / dL;

  const double dSide = (sN.Dot(crCentre - rJet.s) < 0.0) ? -1.0 : 1.0;
  rJet.n   = sN   * dSide; rJet.nu  = sNu  * dSide; rJet.nv  = sNv  * dSide;
  rJet.nuu = sNuu * dSide; rJet.nuv = sNuv * dSide; rJet.nvv = sNvv * dSide;
  rJet.E = rJet.su.Dot(rJet.su);  rJet.F = rJet.su.Dot(rJet.sv);  rJet.G = rJet.sv.Dot(rJet.sv);
  rJet.L = rJet.suu.Dot(rJet.n);  rJet.M = rJet.suv.Dot(rJet.n);  rJet.N = rJet.svv.Dot(rJet.n);
  rbJudged = TRUE;
  return SM_SUCCESS;
} // end sm_FaceJetAt

/*******************************************************************//**
PURPOSE: SolidCraft change. The face's largest principal curvature toward
    the ball: the largest root k of (L - k E)(N - k G) - (M - k F)^2 = 0.

NOTES: The ball of radius r is locally enclosed by the face (every normal
    section curves less toward the ball than the ball does: (E - rL) x^2 +
    2 (F - rM) x y + (G - rN) y^2 positive definite) exactly when 1 - r k > 0 for this k: h(r) = (E - rL)(G - rN)
    - (F - rM)^2 first reaches 0, as r grows from 0, at r = 1 / k.
***********************************************************************/
static double sm_LargestCurvatureTowardBall
  (const sm_FaceJet & crJet)  // in : the face's jet at the contact
{
  const double dA = crJet.L * crJet.N - crJet.M * crJet.M;
  const double dB = crJet.E * crJet.N + crJet.G * crJet.L - 2.0 * crJet.F * crJet.M;
  const double dC = crJet.E * crJet.G - crJet.F * crJet.F;
  const double dDisc = dB * dB - 4.0 * dA * dC;
  return (dB + std::sqrt(dDisc > 0.0 ? dDisc : 0.0)) / (2.0 * dC);
} // end sm_LargestCurvatureTowardBall

/*******************************************************************//**
PURPOSE: SolidCraft change. Solve the 4 x 4 system A x = b by Gaussian
    elimination with partial pivoting; rbSingular where A is singular to
    working precision (its determinant over the product of its columns'
    lengths, Hadamard's bound, at most SM_EFF_ZERO).
***********************************************************************/
static void sm_Solve4
  (const double   cadA[4][4],  // in : A
   const double   cadB[4],     // in : b
   double         adX[4],      // out: x
   SmBoolean    & rbSingular)  // out: TRUE = A singular
{
  double adM[4][5];
  double dHadamard = 1.0;
  for (int j=0; j<4; j++)
    {
      double dCol = 0.0;
      for (int i=0; i<4; i++)
        { dCol += cadA[i][j] * cadA[i][j]; }
      dHadamard *= std::sqrt(dCol);
    }
  for (int i=0; i<4; i++)
    {
      for (int j=0; j<4; j++)
        { adM[i][j] = cadA[i][j]; }
      adM[i][4] = cadB[i];
    }
  double dDet = 1.0;
  for (int k=0; k<4; k++)
    {
      int iPivot = k;
      for (int i=k+1; i<4; i++)
        if (smos_Fabs(adM[i][k]) > smos_Fabs(adM[iPivot][k]))
          { iPivot = i; }
      if (iPivot != k)
        {
          for (int j=0; j<5; j++)
            { const double d = adM[k][j]; adM[k][j] = adM[iPivot][j]; adM[iPivot][j] = d; }
          dDet = -dDet;
        }
      dDet *= adM[k][k];
      if (adM[k][k] == 0.0)
        { break; }
      for (int i=k+1; i<4; i++)
        {
          const double dF = adM[i][k] / adM[k][k];
          for (int j=k; j<5; j++)
            { adM[i][j] -= dF * adM[k][j]; }
        }
    }
  rbSingular = (!(dHadamard > 0.0) || !(smos_Fabs(dDet) > SM_EFF_ZERO * dHadamard)) ? TRUE : FALSE;
  for (int i=0; i<4; i++)
    { adX[i] = 0.0; }
  if (rbSingular)
    { return; }
  for (int i=3; i>=0; i--)
    {
      double d = adM[i][4];
      for (int j=i+1; j<4; j++)
        { d -= adM[i][j] * adX[j]; }
      adX[i] = d / adM[i][i];
    }
} // end sm_Solve4

/*******************************************************************//**
PURPOSE: SolidCraft change. The rolling ball's local regularity at the
    station crPnt (its contacts' uv on the solver's two faces): the blend
    design's §6 checks.

NOTES: With the ball's centre c and radius r, its contacts s1, s2 and the
    faces' unit normals n1, n2 there toward c:
    1. Each face locally encloses the ball: 1 - r k_i > 0, k_i the face's
       largest principal curvature toward the ball (h_i(r) > 0 with the
       form positive definite). Equal (to a millionth, a ball as wide as
       the face's curve: the hemispherical end of a pin as wide as the
       round, two rounds' corner ball) is a focal station, where the ball
       is the blend and the rest is not judged.
    2. The spine is transversal to the faces' bisector: the centre moves,
       c' x (s2 - s1) != 0. The derivatives below solve the station's
       equations c = s_i + r n_i (i = 1, 2) with one more fixing the trace's
       parameter t, and their matrix is singular where this fails.
    3. Every point p of the section is progressive (the envelope does not
       fold back over itself there): for the balls' family f(t, p) = (|p - c|^2 - r^2) / 2,
       f_tt = |c'|^2 - r'^2 - r r'' - <p - c, c''> > 0. On the touch circle
       (the sphere's points with <p - s1, c'> = 0) from s1 to s2, the short
       way, f_tt is linear in p and its least value is in closed form. For a
       section other than the ball's circle only its ends, the rails, are on
       the ball: there f_tt is the contact conformity mu_i = <s_i', c'>, and
       those two are judged.
    The trace's parameter: for a law along the edge (rbOnEdge), the edge's
    t, the centre in the plane <C'(t), C(t) - c> = 0; for a constant ball a
    parameter along the fixed direction w = n1 x n2, <c', w> = 1 and
    <c'', w> = 0 (f_tt / |c'|^2 does not depend on it). c' and c'' come from
    differentiating the equations once and twice: the second needs the
    normals' second derivatives (sm_FaceJetAt). Where the contacts meet
    (the tangent end, an event of its own) only check 1 is judged; a pole is
    not judged.
***********************************************************************/
SmStatus SmFilletSolver::JudgeBallAt
  (const SmTsectPnt & crPnt,      // in : the station: UVPos(0), UVPos(1) and, for a law, its t
   SmBallStation    & rStation)   // out: the station's judgement
{
  rStation = SmBallStation();
  SmBoolean bBall = FALSE, bOnEdge = FALSE;
  double adR[3] = { 0.0, 0.0, 0.0 };
  SmVector3d asC[4];
  SER(GetBallLaw(crPnt, bBall, bOnEdge, adR, asC));
  SmOffsetSurface * apOffset[2] = { GetSurface(0), GetSurface(1) };
  NER(apOffset[0]); NER(apOffset[1]);
  const SmSurface * apBase[2] = { apOffset[0]->GetBaseSurface(), apOffset[1]->GetBaseSurface() };
  NER(apBase[0]); NER(apBase[1]);
  rStation.m_bLaw  = bOnEdge;
  rStation.m_bBall = bBall;

  // each face's offset point at the station's offsets: the ball's centre (both the same), or
  // where no ball rolls (unequal offsets) each face's own
  const SmBoolean bRadiusChanges = OffsetRadiiCanChange();
  if (bRadiusChanges)
    { SER(SetupOffsetValues(crPnt)); }
  SmPoint3d asOn[2];
  SmStatus eOn = SM_SUCCESS;
  for (ULONG k=0; k<2 && eOn == SM_SUCCESS; k++)
    { eOn = apOffset[k]->EvaluatePoint(crPnt.m_vUVCurvePV[k][0], asOn[k]); }
  if (bRadiusChanges)
    { ReloadSurfaceOffsets(); }
  SER(eOn);
  const SmVector3d sC(asOn[0]);

  // each face's surface as its offset evaluates it there: the base, or past the base's domain
  // its extension (made by the offset's own evaluation, which the trace solves on)
  sm_FaceJet asJet[2];
  double adOffset[2] = { 0.0, 0.0 };
  for (ULONG k=0; k<2; k++)
    {
      const SmPoint2d & crUV = crPnt.m_vUVCurvePV[k][0];
      const SmSurface * pActing = apBase[k];
      const SmSurface * pExtended = apOffset[k]->GetExtendedBaseSurface();
      if (   !apBase[k]->GetNaturalUVDomain().ContainsPoint2d(crUV)
          && pExtended != NULL && pExtended->GetNaturalUVDomain().ContainsPoint2d(crUV))
        { pActing = pExtended; }
      SmBoolean bJudged = FALSE;
      SER(sm_FaceJetAt(*pActing, crUV, asOn[k], asJet[k], bJudged));
      if (!bJudged)
        { return SM_SUCCESS; }
      adOffset[k] = (asOn[k] - asJet[k].s).Length();
    }
  const double dR = adOffset[0];
  const SmVector3d sChord = asJet[1].s - asJet[0].s;
  rStation.m_dRadius = bBall ? dR : smos_Max(adOffset[0], adOffset[1]);
  if (!(rStation.m_dRadius > SM_EFF_ZERO))
    { return SM_SUCCESS; }

  // 1. each face encloses the ball: 1 - d k > 0 for its offset d (the ball's radius), the
  // offset not folded (where no ball rolls, each face's own offset)
  const double dFocal = 1.0e-6;
  rStation.m_dEnclosed = SM_BIG_DOUBLE;
  for (ULONG k=0; k<2; k++)
    {
      const double dK = sm_LargestCurvatureTowardBall(asJet[k]);
      const double dMargin = 1.0 - adOffset[k] * dK;
      if (dMargin < rStation.m_dEnclosed)
        {
          rStation.m_dEnclosed = dMargin;
          if (!bBall)
            { rStation.m_dRadius = adOffset[k]; }
        }
      if (dK > 0.0 && 1.0 / dK < rStation.m_dFits)
        { rStation.m_dFits = 1.0 / dK; }
    }
  if (rStation.m_dEnclosed < -dFocal)
    {
      rStation.m_lCheck = SmBallStation::SM_BS_NOT_ENCLOSED;
      return SM_SUCCESS;
    }
  rStation.m_lCheck = SmBallStation::SM_BS_REGULAR;
  if (!bBall)
    { return SM_SUCCESS; } // no ball: its spine and section are not a rolling ball's
  // a focal station, or contacts that meet (the faces tangent there: a tangent end, an event of
  // its own), judge nothing further
  if (rStation.m_dEnclosed <= dFocal || !(sChord.Length() > GetThisApproxTol3d()))
    { return SM_SUCCESS; }

  // 2. the spine's tangent: O_i = s_i + r n_i, O_i' = O_iu u_i' + O_iv v_i' + r' n_i
  const double dRp = adR[1], dRpp = adR[2];
  SmVector3d asOu[2], asOv[2];
  for (ULONG k=0; k<2; k++)
    {
      asOu[k] = asJet[k].su + asJet[k].nu * dR;
      asOv[k] = asJet[k].sv + asJet[k].nv * dR;
    }
  SmVector3d sW;
  double dA1 = 0.0;
  if (bOnEdge)
    {
      sW  = asC[1];
      dA1 = asC[1].Dot(asC[1]) + asC[2].Dot(asC[0] - sC);
    }
  else
    {
      sW = asJet[0].n * asJet[1].n;
      const double dW = sW.Length();
      if (!(dW > SM_EFF_ZERO))
        {
          rStation.m_lCheck = SmBallStation::SM_BS_NOT_TRANSVERSAL;
          return SM_SUCCESS;
        }
      sW  = sW / dW;
      dA1 = 1.0;
    }
  double adA[4][4], adB[4], adX[4];
  for (ULONG i=0; i<3; i++)
    {
      adA[i][0] = asOu[0][i];  adA[i][1] = asOv[0][i];
      adA[i][2] = -asOu[1][i]; adA[i][3] = -asOv[1][i];
      adB[i] = -dRp * (asJet[0].n[i] - asJet[1].n[i]);
    }
  adA[3][0] = sW.Dot(asOu[0]); adA[3][1] = sW.Dot(asOv[0]); adA[3][2] = 0.0; adA[3][3] = 0.0;
  adB[3] = dA1 - dRp * sW.Dot(asJet[0].n);
  SmBoolean bSingular = FALSE;
  sm_Solve4(adA, adB, adX, bSingular);
  const SmVector3d sCp = asOu[0] * adX[0] + asOv[0] * adX[1] + asJet[0].n * dRp;
  const double dCp2 = sCp.LengthSquared();
  if (   bSingular || !(dCp2 > 0.0)
      || !((sCp * sChord).Length() > SM_EFF_ZERO * std::sqrt(dCp2) * sChord.Length()))
    {
      rStation.m_lCheck = SmBallStation::SM_BS_NOT_TRANSVERSAL;
      return SM_SUCCESS;
    }

  // 3. c'': the equations differentiated twice, K_i the terms without u_i'', v_i''
  const double adU[2][2] = { { adX[0], adX[1] }, { adX[2], adX[3] } };
  SmVector3d asK[2];
  for (ULONG k=0; k<2; k++)
    {
      const sm_FaceJet & crJ = asJet[k];
      const double dU = adU[k][0], dV = adU[k][1];
      const SmVector3d sOuu = crJ.suu + crJ.nuu * dR;
      const SmVector3d sOuv = crJ.suv + crJ.nuv * dR;
      const SmVector3d sOvv = crJ.svv + crJ.nvv * dR;
      asK[k] =   sOuu * (dU * dU) + sOuv * (2.0 * dU * dV) + sOvv * (dV * dV)
               + (crJ.nu * dU + crJ.nv * dV) * (2.0 * dRp) + crJ.n * dRpp;
    }
  double dA2 = 0.0;
  if (bOnEdge)
    { dA2 = asC[3].Dot(asC[0] - sC) + 2.0 * asC[2].Dot(asC[1] - sCp) + asC[1].Dot(asC[2]); }
  for (ULONG i=0; i<3; i++)
    { adB[i] = -(asK[0][i] - asK[1][i]); }
  adB[3] = dA2 - sW.Dot(asK[0]);
  sm_Solve4(adA, adB, adX, bSingular);
  if (bSingular)
    {
      rStation.m_lCheck = SmBallStation::SM_BS_NOT_TRANSVERSAL;
      return SM_SUCCESS;
    }
  const SmVector3d sCpp = asK[0] + asOu[0] * adX[0] + asOv[0] * adX[1];

  // f_tt on the section: A0 - <p - c, c''>
  const double dA0 = dCp2 - dRp * dRp - dR * dRpp;
  double dLeast = SM_BIG_DOUBLE;
  for (ULONG k=0; k<2; k++)
    {
      const double dEnd = dA0 - (asJet[k].s - sC).Dot(sCpp);
      if (dEnd < dLeast)
        { dLeast = dEnd; }
    }
  if (m_pFSG != NULL && m_pFSG->GetFilletSurfaceGeneratorType() == SM_FSG_CIRCULAR)
    {
      // the touch circle: centre m on the centre line, through s1 and s2
      const SmVector3d sM = sC + sCp * ((asJet[0].s - sC).Dot(sCp) / dCp2);
      SmVector3d sE1 = asJet[0].s - sM;
      const double dRho = sE1.Length();
      if (dRho > SM_EFF_ZERO * dR)
        {
          sE1 = sE1 / dRho;
          const SmVector3d sQ = asJet[1].s - sM;
          SmVector3d sE2 = sQ - sE1 * sQ.Dot(sE1);
          const double dE2 = sE2.Length();
          if (dE2 > SM_EFF_ZERO * dR)
            { sE2 = sE2 / dE2; }
          else
            {
              sE2 = sCp * sE1;
              SER(sE2.Unitize());
            }
          const double dPhi = std::atan2(sQ.Dot(sE2), sQ.Dot(sE1));
          // the most of a cos(phi) + b sin(phi) over [0, phi]
          const double dA = sE1.Dot(sCpp), dB = sE2.Dot(sCpp);
          double dMost = dA;
          const double dAtEnd = dA * std::cos(dPhi) + dB * std::sin(dPhi);
          if (dAtEnd > dMost)
            { dMost = dAtEnd; }
          const double dStar = std::atan2(dB, dA);
          if (dStar > 0.0 && dStar < dPhi)
            { dMost = std::sqrt(dA * dA + dB * dB); }
          const double dArc = dA0 - (sM - sC).Dot(sCpp) - dRho * dMost;
          if (dArc < dLeast)
            { dLeast = dArc; }
        }
    }
  rStation.m_dProgress = dLeast / dCp2;
  if (!(rStation.m_dProgress > 0.0))
    { rStation.m_lCheck = SmBallStation::SM_BS_REGRESSIVE; }
  return SM_SUCCESS;
} // end SmFilletSolver::JudgeBallAt

/*******************************************************************//**
PURPOSE: SolidCraft change. The solver's section in the plane square to its
    filleted edge at the edge's parameter dT (PointOnPlaneSolve on the
    faces' surfaces extended as the trace extends them), with the ball at
    the offsets the solver's surfaces have. Its guess is the edge point
    dropped on each face, or where that finds none, a section found near
    (continuing from it as the trace does: a wide ball's contacts lie far
    from the edge).
***********************************************************************/
static SmStatus sm_BallSectionAt
  (SmFilletSolver   & rSolver,     // in : the solver
   const SmCurve    & crEdgeCurve, // in : the filleted edge's curve
   double             dT,          // in : the edge's parameter
   const SmTsectPnt * pGuess,      // in : a section found near, whose contacts are the guess (NULL: none)
   SmBoolean        & rbFound,     // out: TRUE = the section was found
   SmTsectPnt      & rPnt)        // out: the section's contacts (and, for a law, its t)
{
  rbFound = FALSE;
  const SmSurface * apBase[2] = { rSolver.GetSurface(0)->GetBaseSurface(), rSolver.GetSurface(1)->GetBaseSurface() };
  NER(apBase[0]); NER(apBase[1]);
  // the domains the trace solves on: the natural ones, extended as for the trace
  SmExtent2d asDomain[2] = { apBase[0]->GetNaturalUVDomain(), apBase[1]->GetNaturalUVDomain() };
  rSolver.SetupOffsetExtension(apBase[0], asDomain[0]);
  rSolver.SetupOffsetExtension(apBase[1], asDomain[1]);
  SmVector3d sPV[2];
  SER(crEdgeCurve.Evaluate(dT, 1, TRUE, sPV));
  if (!(sPV[1].Length() > SM_EFF_ZERO))
    { return SM_SUCCESS; }
  SER(sPV[1].Unitize());
  SmVector2d asGuess[2];
  SmBoolean  bDropped = TRUE;
  for (ULONG k=0; k<2 && bDropped; k++)
    {
      SmBoolean bMulti = FALSE;
      SmPoint2d sUV;
      double    dGap = SM_BIG_DOUBLE;
      SER(apBase[k]->DropPoint(sPV[0], asDomain[k], NULL, bDropped, sUV, dGap, bMulti));
      asGuess[k] = SmVector2d(sUV.x, sUV.y);
    }
  if (bDropped)
    { SER(rSolver.PointOnPlaneSolve(sPV[0], sPV[1], asDomain[0], asDomain[1], asGuess[0], asGuess[1], rbFound, rPnt)); }
  if (!rbFound && pGuess != NULL)
    {
      SER(rSolver.PointOnPlaneSolve(sPV[0], sPV[1], asDomain[0], asDomain[1], pGuess->m_vUVCurvePV[0][0],
                                    pGuess->m_vUVCurvePV[1][0], rbFound, rPnt));
    }
  return SM_SUCCESS;
} // end sm_BallSectionAt

/*******************************************************************//**
PURPOSE: SolidCraft change. Whether the solver's ball is regular at the
    section at the edge's dT (JudgeBallAt), with the ball's radius dRadius
    for a constant ball (0: the solver's own), its surfaces' offsets back at
    the solver's own after.
***********************************************************************/
static SmStatus sm_RegularAt
  (SmFilletSolver  & rSolver,     // in : the solver
   const SmCurve   & crEdgeCurve, // in : the filleted edge's curve
   double            dT,          // in : the edge's parameter
   double            dRadius,     // in : a constant ball's trial radius, 0 for the solver's
   SmTsectPnt      & rNear,       // i/o: a section found near (the guess) and then this one, if found
   SmBoolean       & rbNear,      // i/o: whether rNear holds one
   SmBoolean       & rbRegular,   // out: TRUE = regular there (or not judged)
   SmBallStation   & rStation)    // out: the judgement
{
  rbRegular = TRUE;
  rStation  = SmBallStation();
  double dOwn = 0.0;
  if (dRadius > 0.0)
    {
      SmBoolean bOnEdge = FALSE;
      if (!rSolver.RollsABall(bOnEdge) || bOnEdge)
        { return SM_SUCCESS; }
      SmTsectPnt sAny;
      SmBoolean bBall = FALSE;
      double adR[3] = { 0.0, 0.0, 0.0 };
      SmVector3d asC[4];
      SER(rSolver.GetBallLaw(sAny, bBall, bOnEdge, adR, asC));
      dOwn = adR[0];
      if (!rSolver.SetBallRadius(dRadius))
        { return SM_SUCCESS; }
    }
  SmTsectPnt sPnt;
  SmBoolean  bFound = FALSE;
  SmStatus eStatus = sm_BallSectionAt(rSolver, crEdgeCurve, dT, rbNear ? &rNear : NULL, bFound, sPnt);
  if (eStatus == SM_SUCCESS && bFound)
    {
      eStatus = rSolver.JudgeBallAt(sPnt, rStation);
      rNear   = sPnt;
      rbNear  = TRUE;
    }
  if (dRadius > 0.0)
    { (void)rSolver.SetBallRadius(dOwn); } // the solver's own radius back, whatever the solve did
  SER(eStatus);
  if (dRadius > 0.0 && !bFound)
    { rbRegular = FALSE; } // no section at the trial radius: not one that works
  else if (   rStation.m_lCheck != SmBallStation::SM_BS_UNJUDGED
           && rStation.m_lCheck != SmBallStation::SM_BS_REGULAR)
    { rbRegular = FALSE; }
  return SM_SUCCESS;
} // end sm_RegularAt

/*******************************************************************//**
PURPOSE: SolidCraft change. Whether this solver's rolling ball is regular
    along its edge (JudgeBallAt at sections across the edge), and where it
    is not, the radius that works.

NOTES: For CreateFilletBrep, before the corners (which look for a rail's end
    where it cannot be when the ball does not fit). The sections are at the
    middles of 32 equal steps of the edge's parameter; where all of them are
    regular, each margin (m_dEnclosed, m_dProgress) is minimized about the
    section where it is least (golden section over the steps either side),
    so that a narrow dip between sections is found. rFirst is the station
    that fails (m_lCheck SM_BS_REGULAR where none does), with m_dWorks:
    - a constant ball: the largest radius regular at the failing station,
      by bisection on the radius with the section solved again at each (the
      contacts move with r), the least over the failing stations (0 where
      no radius down to a millionth of r is regular there);
    - a radius law: the radius where the law first stops being regular
      along the edge, by bisection on t from a regular section below.
    A section the solve does not find and a pole are not judged.
***********************************************************************/
SmStatus SmFilletSolver::BallIsRegular
  (SmBallStation & rFirst)  // out: the failing station and the radius that works
{
  rFirst = SmBallStation();
  rFirst.m_lCheck = SmBallStation::SM_BS_REGULAR;
  if (GetSurface(0) == NULL || GetSurface(1) == NULL || GetEdgeuse(0) == NULL)
    { return SM_SUCCESS; }
  const SmEdge  * pEdge  = GetEdgeuse(0)->GetEdge();  NER(pEdge);
  const SmCurve * pCurve = pEdge->GetCurve();         NER(pCurve);
  const SmExtent1d sIvl  = pEdge->GetInterval();
  // a solver that rolls no ball (unequal offsets) has only its offsets judged (check 1)
  SmBoolean bLawBall = FALSE;
  const SmBoolean bRollsABall = RollsABall(bLawBall);
  const ULONG lSteps = 32;
  const ULONG lBisections = 40;
  std::vector<SmTsectPnt> vSections(lSteps);  // each step's section, the guess about it (vFound: one was found)
  std::vector<SmBoolean>  vFound(lSteps, FALSE);

  // the radius that works at a failing section dT (judged crFail), from a regular one at
  // dPassT below it (bPass) for a law, the solves continuing from crNear (bHaveNear)
  auto Fail = [&](double dT, const SmBallStation & crFail, SmBoolean bPass, double dPassT,
                  const SmTsectPnt & crNear, SmBoolean bHaveNear) -> SmStatus
  {
    SmTsectPnt sNear = crNear;
    SmBoolean  bNear = bHaveNear;
    if (crFail.m_bLaw)
      {
        SmBallStation sOnset = crFail;
        if (bPass)
          {
            double dLo = dPassT, dHi = dT;
            for (ULONG i=0; i<lBisections; i++)
              {
                const double dMid = 0.5 * (dLo + dHi);
                SmBoolean bMid = TRUE;
                SmBallStation sMid;
                SER(sm_RegularAt(*this, *pCurve, dMid, 0.0, sNear, bNear, bMid, sMid));
                if (bMid)
                  { dLo = dMid; }
                else
                  {
                    dHi = dMid;
                    sOnset = sMid;
                  }
              }
          }
        if (rFirst.m_lCheck == SmBallStation::SM_BS_REGULAR)
          {
            rFirst = sOnset;
            rFirst.m_dWorks = sOnset.m_dRadius;
          }
        return SM_SUCCESS;
      }
    if (!bRollsABall)
      {
        // unequal offsets: the offset that folds must stay below the face's radius of curvature
        if (crFail.m_lCheck == SmBallStation::SM_BS_NOT_ENCLOSED
            && (rFirst.m_lCheck == SmBallStation::SM_BS_REGULAR || crFail.m_dFits < rFirst.m_dWorks))
          {
            rFirst = crFail;
            rFirst.m_dWorks = crFail.m_dFits;
          }
        return SM_SUCCESS;
      }
    const double dRadius = crFail.m_dRadius;
    double dLo = 0.0, dHi = dRadius;
    SmBallStation sLo;
    for (double dTry = 0.5 * dRadius; dTry > 1.0e-6 * dRadius; dTry *= 0.5)
      {
        SmBoolean bTry = TRUE;
        SmBallStation sTry;
        SER(sm_RegularAt(*this, *pCurve, dT, dTry, sNear, bNear, bTry, sTry));
        if (bTry)
          {
            dLo = dTry;
            sLo = sTry;
            break;
          }
        dHi = dTry;
      }
    if (dLo > 0.0)
      {
        for (ULONG i=0; i<lBisections; i++)
          {
            const double dMid = 0.5 * (dLo + dHi);
            SmBoolean bMid = TRUE;
            SmBallStation sMid;
            SER(sm_RegularAt(*this, *pCurve, dT, dMid, sNear, bNear, bMid, sMid));
            if (bMid)
              {
                dLo = dMid;
                sLo = sMid;
              }
            else
              { dHi = dMid; }
          }
      }
    // a ball that does not fit: the limit is where h(r) reaches 0 at the contacts of the largest
    // regular radius (the bisection itself stops at the focal tolerance, a millionth past it)
    if (dLo > 0.0 && crFail.m_lCheck == SmBallStation::SM_BS_NOT_ENCLOSED && sLo.m_dFits < SM_BIG_DOUBLE)
      { dLo = sLo.m_dFits; }
    if (rFirst.m_lCheck == SmBallStation::SM_BS_REGULAR || dLo < rFirst.m_dWorks)
      {
        rFirst = crFail;
        rFirst.m_dWorks = dLo;
      }
    return SM_SUCCESS;
  };

  SmTArray<double> sT;
  SmTArray<double> asMargin[2];
  SmBoolean bPassed = FALSE;
  double    dPassedT = 0.0;
  SmTsectPnt sNear;
  SmBoolean  bNear = FALSE;
  for (ULONG j=0; j<lSteps; j++)
    {
      const double dT = sIvl.Evaluate((j + 0.5) / lSteps);
      SmBoolean bRegular = TRUE;
      SmBallStation sStation;
      SER(sm_RegularAt(*this, *pCurve, dT, 0.0, sNear, bNear, bRegular, sStation));
      vSections[j] = sNear; // this step's section, or the last one found before it
      vFound[j]    = bNear;
      sT.Add(dT);
      asMargin[0].Add(sStation.m_dEnclosed);
      asMargin[1].Add(sStation.m_dProgress);
      if (bRegular)
        {
          bPassed  = TRUE;
          dPassedT = dT;
          continue;
        }
      SER(Fail(dT, sStation, bPassed, dPassedT, sNear, bNear));
      if (sStation.m_bLaw)
        { return SM_SUCCESS; } // a law: the first failing section along the edge
    }
  if (rFirst.m_lCheck != SmBallStation::SM_BS_REGULAR)
    { return SM_SUCCESS; }

  // every section regular: each margin's least between the sections about its least one
  for (ULONG m=0; m<2; m++)
    {
      ULONG lLeast = 0;
      for (ULONG j=1; j<sT.GetSize(); j++)
        if (asMargin[m][j] < asMargin[m][lLeast])
          { lLeast = j; }
      if (sT.GetSize() == 0 || !(asMargin[m][lLeast] < SM_BIG_DOUBLE))
        { continue; }
      double dA = (lLeast > 0) ? sT[lLeast-1] : sIvl.GetMin();
      double dB = (lLeast + 1 < sT.GetSize()) ? sT[lLeast+1] : sIvl.GetMax();
      const double dGolden = 0.5 * (std::sqrt(5.0) - 1.0);
      sNear = vSections[lLeast];
      bNear = vFound[lLeast];
      auto Margin = [&](double dT, double & rdMargin, SmBallStation & rStation) -> SmStatus
      {
        SmBoolean bRegular = TRUE;
        SER(sm_RegularAt(*this, *pCurve, dT, 0.0, sNear, bNear, bRegular, rStation));
        rdMargin = (m == 0) ? rStation.m_dEnclosed : rStation.m_dProgress;
        if (!bRegular)
          { rdMargin = -SM_BIG_DOUBLE; } // failing: as low as it gets
        return SM_SUCCESS;
      };
      double dX1 = dB - dGolden * (dB - dA), dX2 = dA + dGolden * (dB - dA);
      double dF1 = 0.0, dF2 = 0.0;
      SmBallStation s1, s2;
      SER(Margin(dX1, dF1, s1));
      SER(Margin(dX2, dF2, s2));
      for (ULONG i=0; i<lBisections && dF1 > -SM_BIG_DOUBLE && dF2 > -SM_BIG_DOUBLE; i++)
        {
          if (dF1 < dF2)
            {
              dB = dX2; dX2 = dX1; dF2 = dF1; s2 = s1;
              dX1 = dB - dGolden * (dB - dA);
              SER(Margin(dX1, dF1, s1));
            }
          else
            {
              dA = dX1; dX1 = dX2; dF1 = dF2; s1 = s2;
              dX2 = dA + dGolden * (dB - dA);
              SER(Margin(dX2, dF2, s2));
            }
        }
      const SmBoolean bFirstFails = (dF1 <= dF2) ? TRUE : FALSE;
      const double dFailT = bFirstFails ? dX1 : dX2;
      const SmBallStation & crFail = bFirstFails ? s1 : s2;
      if (   crFail.m_lCheck == SmBallStation::SM_BS_UNJUDGED
          || crFail.m_lCheck == SmBallStation::SM_BS_REGULAR)
        { continue; }
      // a regular section below the failing one: the least section's, or the one before it
      const double dPassT = (sT[lLeast] < dFailT) ? sT[lLeast] : ((lLeast > 0) ? sT[lLeast-1] : sIvl.GetMin());
      const SmBoolean bPass = (sT[lLeast] < dFailT || lLeast > 0) ? TRUE : FALSE;
      SER(Fail(dFailT, crFail, bPass, dPassT, vSections[lLeast], vFound[lLeast]));
    }
  return SM_SUCCESS;
} // end SmFilletSolver::BallIsRegular

/*******************************************************************//**
PURPOSE: SolidCraft change. The points where the body's edges cross the
    plane through crOrigin square to crNormal: each edge's curve intersected
    with the plane (SmSurface::GlobalCurveIntersect); an edge lying in the
    plane gives its ends there.
***********************************************************************/
static SmStatus sm_EdgesCrossingPlane
  (const SmContext     & crContext,  // in : for the plane
   const SmBrep        & crBody,     // in : the body
   const SmPoint3d     & crOrigin,   // in : a point of the plane
   const SmVector3d    & crNormal,   // in : its unit normal
   SmTArray<SmPoint3d> & rCrossings) // out: the crossings
{
  rCrossings.RemoveAll();
  SmVector3d sX = (smos_Fabs(crNormal.x) < 0.9) ? SmVector3d(1.0, 0.0, 0.0) * crNormal : SmVector3d(0.0, 1.0, 0.0) * crNormal;
  SER(sX.Unitize());
  const SmVector3d sY = crNormal * sX;
  SmPlane * pPlane = NULL;
  SER(SmPlane::CreateCanonical(crContext, SmAxis2Placement(crOrigin, sX, sY), pPlane));
  NER(pPlane);
  SmObjDelete sDeletePlane(pPlane);
  const double dTol = crBody.GetTolerance().val;
  SmTArray<SmEdge*> sEdges;
  crBody.GetEdges(sEdges);
  for (ULONG i=0; i<sEdges.GetSize(); i++)
    {
      const SmEdge  * pEdge  = sEdges[i];
      const SmCurve * pCurve = (pEdge != NULL) ? pEdge->GetCurve() : NULL;
      if (pCurve == NULL)
        { continue; }
      SmSolution aData[4];
      SmSolutionArray sSolutions(4, aData);
      SER(pPlane->GlobalCurveIntersect(pPlane->GetNaturalUVDomain(), *pCurve, pEdge->GetInterval(), dTol, sSolutions));
      for (ULONG j=0; j<sSolutions.GetSize(); j++)
        {
          SmPoint3d sAt;
          SER(pCurve->EvaluatePoint(sSolutions[j].m_vStart[0], sAt));
          rCrossings.Add(sAt);
          if (sSolutions[j].m_eSolutionType == SM_ST_RANGE_OF_VALUES)
            {
              SER(pCurve->EvaluatePoint(sSolutions[j].m_vEnd[0], sAt));
              rCrossings.Add(sAt);
            }
        }
    }
  return SM_SUCCESS;
} // end sm_EdgesCrossingPlane

/*******************************************************************//**
PURPOSE: SolidCraft change. Whether the part keeps material behind the
    section of a convex blend at a station: the region between the section
    and the ball's centre c (the ball's sector between the rails for a
    round, the triangle c, s1, s2 for a chamfer's chord).

NOTES: Kept where a body edge crosses the station's plane inside that region
    (the part's boundary reaches behind the section), or the section meets
    the part between its rails. The latter is decided exactly: the section
    (the ball's circular arc, or the chord, as an exact curve) is
    intersected with every face's surface (SmSurface::GlobalCurveIntersect)
    and each crossing classified on its face; a crossing inside a face is
    the section meeting the part's boundary (one on a face's edge, within
    the intersection's tolerance of it, decides nothing: where the section
    passes a body edge the edge's crossing of the plane decides, exactly,
    above). With none, the
    section between its rails lies wholly inside or wholly outside the part,
    and one point of it, classified, says which. Otherwise the blend takes
    all of the part between its rails there: the slice is removed whole.
***********************************************************************/
static SmStatus sm_KeepsMaterialBehind
  (const SmContext           & crContext,   // in : for the section's curve
   const SmBrep              & crBody,      // in : the body
   const SmPoint3d           & crCentre,    // in : the ball's centre
   const SmPoint3d           & crRail1,     // in : the section's ends
   const SmPoint3d           & crRail2,
   SmBoolean                   bChord,      // in : TRUE: the section is the chord, else the ball's arc
   const SmTArray<SmPoint3d> & crCrossings, // in : the body's edges crossing the station's plane
   SmBoolean                 & rbKeeps)     // out: TRUE = material behind the section
{
  rbKeeps = TRUE;
  const SmVector3d sA = crRail1 - crCentre, sB = crRail2 - crCentre;
  const double dR2 = smos_Max(sA.LengthSquared(), sB.LengthSquared());
  const double dAA = sA.Dot(sA), dAB = sA.Dot(sB), dBB = sB.Dot(sB);
  const double dDet = dAA * dBB - dAB * dAB;
  if (!(dDet > SM_EFF_ZERO * dAA * dBB))
    { return SM_SUCCESS; } // rails opposite or together: no region to judge
  const double dTol = 1.0e-9 * std::sqrt(dR2);
  for (ULONG i=0; i<crCrossings.GetSize(); i++)
    {
      const SmVector3d sX = crCrossings[i] - crCentre;
      const double dXA = sX.Dot(sA), dXB = sX.Dot(sB);
      const double dAlpha = (dXA * dBB - dXB * dAB) / dDet;
      const double dBeta  = (dXB * dAA - dXA * dAB) / dDet;
      if (dAlpha < -1.0e-12 || dBeta < -1.0e-12)
        { continue; }
      if (bChord ? (dAlpha + dBeta <= 1.0 + 1.0e-12) : (sX.Length() <= std::sqrt(dR2) + dTol))
        { return SM_SUCCESS; }
    }

  // the section as an exact curve: the chord, or the circle's arc from the first rail to the
  // second the short way
  SmCurve * pSection = NULL;
  SmPoint3d sMiddle;
  if (bChord)
    {
      SmLine * pLine = NULL;
      SER(SmLine::CreateLineSegment(crContext, 3, crRail1, crRail2, pLine));
      pSection = pLine;
      sMiddle = (crRail1 + crRail2) * 0.5;
    }
  else
    {
      const double dRho = std::sqrt(dAA);
      SmVector3d sE1 = sA / dRho;
      SmVector3d sE2 = sB - sA * (dAB / dAA);
      SER(sE2.Unitize());
      const double dPhi = std::atan2(sB.Dot(sE2), sB.Dot(sE1));
      pSection = new (crContext) SmCircle(crCentre, sE1, sE2, SmExtent1d(0.0, SM_RAD2DEG(dPhi)), dRho, 3, &crContext);
      sMiddle = crCentre + sE1 * (dRho * std::cos(0.5 * dPhi)) + sE2 * (dRho * std::sin(0.5 * dPhi));
    }
  NER(pSection);
  SmObjDelete sDeleteSection(pSection);
  const SmExtent1d sIvl = pSection->GetNaturalInterval();
  const double dEnds = crBody.GetTolerance().val;
  SmTArray<SmFace*> sFaces;
  crBody.GetFaces(sFaces);
  for (ULONG i=0; i<sFaces.GetSize(); i++)
    {
      const SmFace    * pFace = sFaces[i];
      const SmSurface * pSurf = (pFace != NULL) ? pFace->GetSurface() : NULL;
      if (pSurf == NULL)
        { continue; }
      SmSolution aData[4];
      SmSolutionArray sSolutions(4, aData);
      SER(pSurf->GlobalCurveIntersect(pSurf->GetNaturalUVDomain(), *pSection, sIvl, dEnds, sSolutions));
      for (ULONG j=0; j<sSolutions.GetSize(); j++)
        {
          const SmSolution & crSol = sSolutions[j];
          if (crSol.m_eSolutionType == SM_ST_RANGE_OF_VALUES)
            { return SM_SUCCESS; } // the section runs along the face's surface
          SmPoint3d sAt;
          SER(pSection->EvaluatePoint(crSol.m_vStart[0], sAt));
          if (sAt.DistanceBetween(crRail1) <= dEnds || sAt.DistanceBetween(crRail2) <= dEnds)
            { continue; } // the section's own ends, on the blended faces
          SmPointClassification sOnFace;
          SER(pFace->PointClassify(SmPoint2d(crSol.m_vStart[1], crSol.m_vStart[2]), SmTol::GetZoneTol3d(&crBody),
                                   TRUE, TRUE, sOnFace));
          if (sOnFace.GetPointClass() == SM_PC_FACE)
            { return SM_SUCCESS; } // the section crosses the part's boundary inside a face
        }
    }
  // no crossing: the section between its rails is wholly in the part or wholly outside it
  SmPointClassification sClass;
  SER(crBody.Point3DClassify(sMiddle, SmTol::GetZoneTol3d(&crBody), TRUE, sClass));
  if (sClass.GetPointClass() != SM_PC_REGION)
    { return SM_SUCCESS; }
  const SmRegion * pRegion = static_cast<const SmRegion *>(sClass.GetObject());
  if (pRegion != NULL && !pRegion->IsVoid())
    { return SM_SUCCESS; }
  rbKeeps = FALSE;
  return SM_SUCCESS;
} // end sm_KeepsMaterialBehind

/*******************************************************************//**
PURPOSE: SolidCraft change. Whether this convex blend (one that removes
    material) would take all of the part between its rails somewhere along
    its edge, cutting it into pieces (a round or chamfer through the bridge
    over a slot); and where it would, the largest size that keeps the part
    whole.

NOTES: For CreateFilletBrep, before the corners. At the same sections as
    BallIsRegular (the middles of 32 equal steps of the edge), with the
    body's edges crossing each section's plane found once. A section where
    no material lies behind the blend's section (sm_KeepsMaterialBehind) is
    a cut where sections either side of its run along the edge keep material
    (a run reaching an end of the edge, a tip thinner than the band, blends
    that end away). The limit is then found at each cutting section by
    bisection on the ball's radius, the section solved again at each, the
    least over the sections: rdRadius the ball's radius there, rdSetback the distance from
    the edge to the first rail (a chamfer's size). A constant ball only.
***********************************************************************/
SmStatus SmFilletSolver::BandCutsThePart
  (const SmBrep & crBody,     // in : the body blended
   SmBoolean    & rbCuts,     // out: TRUE = the blend would cut the part
   double       & rdRadius,   // out: the largest radius that keeps it whole
   double       & rdSetback)  // out: and that radius's setback from the edge
{
  rbCuts = FALSE;
  rdRadius = rdSetback = 0.0;
  if (GetSurface(0) == NULL || GetSurface(1) == NULL || GetEdgeuse(0) == NULL)
    { return SM_SUCCESS; }
  const SmEdge  * pEdge  = GetEdgeuse(0)->GetEdge();  NER(pEdge);
  const SmCurve * pCurve = pEdge->GetCurve();         NER(pCurve);
  const SmExtent1d sIvl  = pEdge->GetInterval();
  SmBoolean bOnEdge = FALSE;
  if (!RollsABall(bOnEdge) || bOnEdge)
    { return SM_SUCCESS; } // a constant ball only
  SmTsectPnt sProbe;
  SmBoolean bBall = FALSE;
  double adR[3] = { 0.0, 0.0, 0.0 };
  SmVector3d asC[4];
  SER(GetBallLaw(sProbe, bBall, bOnEdge, adR, asC));
  const SmBoolean bChord = (m_pFSG != NULL && m_pFSG->GetFilletSurfaceGeneratorType() == SM_FSG_LINEAR) ? TRUE : FALSE;

  // at the edge's dT with the ball's radius dTrial (0: the solver's): whether material is kept
  // behind the section, and the section's centre and rails
  auto Keeps = [&](double dT, double dTrial, const SmTArray<SmPoint3d> & crCrossings, SmTsectPnt & rNear,
                   SmBoolean & rbNear, SmBoolean & rbJudged, SmBoolean & rbKeeps, double & rdSetbackHere) -> SmStatus
  {
    rbJudged = FALSE;
    rbKeeps  = TRUE;
    if (dTrial > 0.0 && !SetBallRadius(dTrial))
      { return SM_SUCCESS; }
    SmTsectPnt sPnt;
    SmBoolean  bFound = FALSE;
    SmStatus eStatus = sm_BallSectionAt(*this, *pCurve, dT, rbNear ? &rNear : NULL, bFound, sPnt);
    SmPoint3d asAt[3];
    if (eStatus == SM_SUCCESS && bFound)
      {
        eStatus = GetSurface(0)->EvaluatePoint(sPnt.UVPos(0), asAt[0]);
        if (eStatus == SM_SUCCESS) { eStatus = GetSurface(0)->GetBaseSurface()->EvaluatePoint(sPnt.UVPos(0), asAt[1]); }
        if (eStatus == SM_SUCCESS) { eStatus = GetSurface(1)->GetBaseSurface()->EvaluatePoint(sPnt.UVPos(1), asAt[2]); }
      }
    // the rails past a base surface's domain lie on its extension, as the offset evaluates it
    for (ULONG k=0; k<2 && eStatus == SM_SUCCESS && bFound; k++)
      {
        const SmSurface * pExtended = GetSurface(k)->GetExtendedBaseSurface();
        if (   !GetSurface(k)->GetBaseSurface()->GetNaturalUVDomain().ContainsPoint2d(sPnt.UVPos(k))
            && pExtended != NULL && pExtended->GetNaturalUVDomain().ContainsPoint2d(sPnt.UVPos(k)))
          { eStatus = pExtended->EvaluatePoint(sPnt.UVPos(k), asAt[k+1]); }
      }
    if (dTrial > 0.0)
      { (void)SetBallRadius(adR[0]); }
    SER(eStatus);
    if (!bFound)
      { return SM_SUCCESS; }
    rNear  = sPnt;
    rbNear = TRUE;
    rbJudged = TRUE;
    SmPoint3d sP;
    SER(pCurve->EvaluatePoint(dT, sP));
    rdSetbackHere = sP.DistanceBetween(asAt[1]);
    SER(sm_KeepsMaterialBehind(m_crContext, crBody, asAt[0], asAt[1], asAt[2], bChord, crCrossings, rbKeeps));
    return SM_SUCCESS;
  };

  const ULONG lSteps = 32;
  const ULONG lBisections = 40;
  // each section: its parameter, plane crossings and section found near, and whether it keeps
  // material behind it (0), takes all of the part between the rails (1), or is not judged (2)
  std::vector<double>              vT(lSteps, 0.0);
  std::vector<SmTArray<SmPoint3d>> vCrossings(lSteps);
  std::vector<SmTsectPnt>          vNear(lSteps);
  std::vector<int>                 vState(lSteps, 2);
  SmTsectPnt sNear;
  SmBoolean  bNear = FALSE;
  for (ULONG j=0; j<lSteps; j++)
    {
      const double dT = sIvl.Evaluate((j + 0.5) / lSteps);
      vT[j] = dT;
      SmVector3d sPV[2];
      SER(pCurve->Evaluate(dT, 1, TRUE, sPV));
      if (!(sPV[1].Length() > SM_EFF_ZERO))
        { continue; }
      SER(sPV[1].Unitize());
      SER(sm_EdgesCrossingPlane(m_crContext, crBody, sPV[0], sPV[1], vCrossings[j]));
      SmBoolean bJudged = FALSE, bKeeps = TRUE;
      double dSetback = 0.0;
      SER(Keeps(dT, 0.0, vCrossings[j], sNear, bNear, bJudged, bKeeps, dSetback));
      vNear[j] = sNear;
      if (bJudged)
        { vState[j] = bKeeps ? 0 : 1; }
    }

  // a run of sections that take all of the part cuts it in two only where sections either side
  // of the run along the edge keep material: a run that reaches an end of the edge (a tip, where
  // the part runs out thinner than the band) blends that end away and leaves the part whole
  for (ULONG j=0; j<lSteps; j++)
    {
      if (vState[j] != 1)
        { continue; }
      ULONG lFirst = j, lLast = j;
      while (lLast + 1 < lSteps && vState[lLast+1] != 0)
        { lLast++; }
      j = lLast;
      SmBoolean bKeptBefore = FALSE, bKeptAfter = FALSE;
      for (ULONG k=0; k<lFirst; k++)
        if (vState[k] == 0) { bKeptBefore = TRUE; }
      for (ULONG k=lLast+1; k<lSteps; k++)
        if (vState[k] == 0) { bKeptAfter = TRUE; }
      if (!bKeptBefore || !bKeptAfter)
        { continue; }
      for (ULONG m=lFirst; m<=lLast; m++)
        {
          if (vState[m] != 1)
            { continue; }
          const double dT = vT[m];
          // the largest radius that keeps material behind its section here
          SmTsectPnt sAt = vNear[m];
          SmBoolean  bAt = TRUE;
          const double dRadius = adR[0];
          double dLo = 0.0, dHi = dRadius, dLoSetback = 0.0;
          for (double dTry = 0.5 * dRadius; dTry > 1.0e-6 * dRadius; dTry *= 0.5)
            {
              SmBoolean bTryJudged = FALSE, bTryKeeps = TRUE;
              double dTrySetback = 0.0;
              SER(Keeps(dT, dTry, vCrossings[m], sAt, bAt, bTryJudged, bTryKeeps, dTrySetback));
              if (bTryJudged && bTryKeeps)
                {
                  dLo = dTry;
                  dLoSetback = dTrySetback;
                  break;
                }
              dHi = dTry;
            }
          for (ULONG i=0; i<lBisections && dLo > 0.0; i++)
            {
              const double dMid = 0.5 * (dLo + dHi);
              SmBoolean bMidJudged = FALSE, bMidKeeps = TRUE;
              double dMidSetback = 0.0;
              SER(Keeps(dT, dMid, vCrossings[m], sAt, bAt, bMidJudged, bMidKeeps, dMidSetback));
              if (bMidJudged && bMidKeeps)
                {
                  dLo = dMid;
                  dLoSetback = dMidSetback;
                }
              else
                { dHi = dMid; }
            }
          if (!rbCuts || dLo < rdRadius)
            {
              rbCuts    = TRUE;
              rdRadius  = dLo;
              rdSetback = dLoSetback;
            }
        }
    }
  return SM_SUCCESS;
} // end SmFilletSolver::BandCutsThePart

/*******************************************************************//**
PURPOSE: SolidCraft change. Whether a rail leaves the face it rolls on, at crLeave (rail
    parameter dRailParam), past an end of the filleted edge, where the corner at that end
    trims it; otherwise it leaves the face beside the edge.

NOTES: At an end the fillet stops in one of two ways, and its rail leaves the face only after
    it:
    - The face ends there in its side edge (the edge before or after the filleted edge in the
      face's loop, at the end's vertex), where the end face meets it: the end face cuts the
      fillet off. A square end face does so in the cross section at the end (the plane
      normal to the edge there, SmFilletVertex::CalcVertOnCrossSection); one sloped either
      way, or curved, before or beyond it. The rail leaves across that side edge, or through
      one of its vertices (through the side edge's far corner).
    - The face goes on past the end (a concave corner, or a short side edge that the rail
      passes, onto the faces beyond it): the rail stays on the face up to the cross section
      at the end and beyond, and the corner trims it there. Where it leaves the face later
      does not matter: the edge point nearest to where it leaves is the end.
    A rail that leaves across any other edge, at a point beside the edge (whose nearest edge
    point is not an end), leaves the face before the fillet's end: rbPastAnEnd is FALSE.
    ppEnd (if not NULL) gets the filleted edge's vertex at the end it leaves past, pbBySide
    whether it leaves across that end's side edge (the first way).
***********************************************************************/
SmStatus SmFilletSolver::LeavesPastAnEnd
  (const SmFace                & crFace,      // in : the face the rail rolls on
   const SmCurve               & crRail,      // in : the rail
   double                        dRailParam,  // in : rail parameter where it leaves the face
   const SmPointClassification & crLeave,     // in : where it leaves the face (an edge or a vertex of its boundary)
   SmBoolean                   & rbPastAnEnd, // out: TRUE = past an end of the filleted edge
   const SmVertex             ** ppEnd,       // out: (optional) that end's vertex
   SmBoolean                   * pbBySide)    // out: (optional) TRUE = across or through that end's side edge
  const
{
  rbPastAnEnd = FALSE;
  if (ppEnd != NULL)
    { *ppEnd = NULL; }
  if (pbBySide != NULL)
    { *pbBySide = FALSE; }
  SmEdgeuse * pFilEU     = GetEdgeuse(0);        NER(pFilEU);
  SmEdge    * pFilEdge   = pFilEU->GetEdge();    NER(pFilEdge);
  SmCurve   * pEdgeCurve = pFilEdge->GetCurve(); NER(pEdgeCurve);

  // Across or through the face's side edge at an end.
  const SmEdge   * pLeaveEdge   = (crLeave.GetPointClass() == SM_PC_EDGE)   ? (const SmEdge*)crLeave.GetObject()   : NULL;
  const SmVertex * pLeaveVertex = (crLeave.GetPointClass() == SM_PC_VERTEX) ? (const SmVertex*)crLeave.GetObject() : NULL;
  SmTArray<SmEdgeuse*> sUses;
  SER( pFilEdge->GetEdgeuses( sUses ));
  for (ULONG i=0; i<sUses.GetSize() && !rbPastAnEnd; i++)
    {
      SmEdgeuse * pUse = sUses[i];
      if (pUse == NULL || pUse->GetFace() != &crFace || !pUse->IsLoopEdgeuse())
        { continue; }
      for (SmEdgeuse * pNext : { pUse->GetCWEdgeuse(), pUse->GetCCWEdgeuse() })
        {
          const SmEdge * pSide = (pNext != NULL) ? pNext->GetEdge() : NULL;
          if (pSide == NULL || pSide == pFilEdge || rbPastAnEnd)
            { continue; }
          rbPastAnEnd =    (pLeaveEdge != NULL && pLeaveEdge == pSide)
                        || (pLeaveVertex != NULL && (   pLeaveVertex == pSide->GetVertex()
                                                     || pLeaveVertex == pSide->GetOtherVertex( pSide->GetVertex() )));
          // the end: the filleted edge's vertex the side edge starts at
          if (rbPastAnEnd && ppEnd != NULL)
            {
              const SmVertex * pSideV = pSide->GetVertex();
              *ppEnd = (pSideV == pFilEdge->GetStartVertex() || pSideV == pFilEdge->GetEndVertex()) ? pSideV : pSide->GetOtherVertex( pSideV );
            }
        }
    }
  if (rbPastAnEnd)
    {
      if (pbBySide != NULL)
        { *pbBySide = TRUE; }
      return SM_SUCCESS;
    }

  // On or beyond the cross section at an end.
  SmExtent1d sEdgeIvl = pFilEdge->GetInterval();
  SmPoint3d  sLeave;
  SER( crRail.EvaluatePoint( dRailParam, sLeave ));
  SmBoolean bDropped   = FALSE;
  double    dEdgeParam = 0.0, dDist = 0.0;
  SER( pEdgeCurve->DropPoint( sEdgeIvl, sLeave, NULL, SM_BIG_DOUBLE, NULL, bDropped, dEdgeParam, dDist ));
  if (!bDropped)
    { SER( SM_ERR ); }
  SmPoint3d sFoot, sEnd;
  SER( pEdgeCurve->EvaluatePoint( dEdgeParam, sFoot ));
  for (ULONG i=0; i<2 && !rbPastAnEnd; i++)
    {
      SER( pEdgeCurve->EvaluatePoint( (i == 0) ? sEdgeIvl.GetMin() : sEdgeIvl.GetMax(), sEnd ));
      rbPastAnEnd = ( sFoot.DistanceBetween( sEnd ) <= GetThisApproxTol3d() );
    }
  // the end: the filleted edge's vertex at that end of its curve
  if (rbPastAnEnd && ppEnd != NULL)
    {
      const SmVertex * pStart = pFilEdge->GetStartVertex(), * pEnd = pFilEdge->GetEndVertex();
      NER(pStart); NER(pEnd);
      *ppEnd = ( pStart->GetPoint().DistanceBetween( sEnd ) <= pEnd->GetPoint().DistanceBetween( sEnd )) ? pStart : pEnd;
    }
  return SM_SUCCESS;
} // end SmFilletSolver::LeavesPastAnEnd

/*******************************************************************//**
Static Functions:
***********************************************************************/

/****************************************************************
PURPOSE: Calculate the intersection of a rail with a (side-)edge

NOTES:
   This puts the result into the TSectPnt of the FilletVertexuse of pFV.

   If the resulting uv point is on a seam of a closed surface, that
   uv value can be on the wrong side of the seam.  Unfortunately, we
   can't tell at this point which side of the seam is the correct side.
   We have to check for that in SmFilletGeom::ReCalcFilletGeom().  [B398]
****************************************************************/
static SmStatus sm_CalcRailXEdge
  (SmFilletVertex * pFV,
   SmFilletGeom * pFG,
   ULONG lRailIndex,
   SmFilletEdge * pRail,
   SmEdge * pSideEdge,
   SmPoint3d * pGuess3DPnt,
   int iDebugLevel)
{
    // Find the right edgeuse for solver to work with.
    //cbi398 Seams?
    SmFace * pRailFace = pRail->GetOriginalFace();
    SmEdgeuse * pSideEU = NULL;
    SmTArray<SmEdgeuse*> sEUs;
    SER( pSideEdge->GetEdgeuses(sEUs) );
    for (ULONG i=0; i<sEUs.GetSize(); i++) {
        if (sEUs[i]->GetFace() == pRailFace) {
            pSideEU = sEUs[i];  // Note, for a seam, all are on pRailFace. [B398]
            break;
        }
    }
    NER(pSideEU);
    SmExtent1d sIvl = pSideEdge->GetInterval();
    SmBoolean bHaveGuessT = FALSE;
    double dGuessParam = 0.0;
    if (pGuess3DPnt) {
        SmCurve * p3DCurve = pSideEdge->GetCurve();
        SmSolutionArray sSolutions;
        SER(p3DCurve->GlobalPointSolve(sIvl,SM_SO_MINIMIZE,*pGuess3DPnt,
            SM_EFF_ZERO_SQRT,NULL,NULL,SM_SR_SINGLE,sSolutions));
        if (sSolutions.GetSize() < 1) { SER(SM_ERR); }
        bHaveGuessT = TRUE;
        dGuessParam = sSolutions[0].m_vStart[0];
    }

    // Call solver
    SmFilletSolver    * pFS = pFG->GetFilletSolver();
    SmFilletVertexuse * pVU = pFV->GetVUAtRailEnd(pFG);
    SmTsectPnt & rTsectPnt = pVU->GetTsectPnt();
    SmBoolean bFoundIntersection = false;
    double dEdgeuseParameter;
    for (ULONG j=1; j<=9; j++) {
        if (!bHaveGuessT) {
            dGuessParam = sIvl.Evaluate(0.1*j);
        }
        SER(pFS->RailEdgeuseIntersect(lRailIndex,pSideEU,dGuessParam,
            bFoundIntersection,rTsectPnt,dEdgeuseParameter));
        if (bHaveGuessT || bFoundIntersection) break;
    }
    if (!bFoundIntersection) { SER(SM_ERR); }

    SmFace * pFace = pRail->GetOriginalFace(); NER(pFace);
    SmPoint2d sUV = rTsectPnt.UVPos(lRailIndex);
    SmPoint3d sPnt;
    SER(pFace->GetSurface()->EvaluatePoint(sUV,sPnt));
    pFV->SetOriginalUV(sUV);
    pFV->SetPoint(sPnt);
    pFV->SetStatus(SM_FIL_PROCESSED);
#ifdef SM_DEBUG_CODE
    if ( iDebugLevel > 0 ) {
        if (pGuess3DPnt) {
            smgfx_SetLook(3,5, 1,0,1); pGuess3DPnt->Draw(); sm_GraphicsLoop();
            smgfx_SetLook(3,4, 0,0,1); pSideEdge->Draw(); sm_GraphicsLoop();
        }
        smgfx_SetLook(4,6, 1,0,0); sPnt.Draw(); sm_GraphicsLoop();
        smgfx_SetLook(3,4, 1,1,0); if (pRail->GetCurve()) { pRail->Draw(); } sm_GraphicsLoop();
        sm_GraphicsLoop();
    }
#else
    SM_REF1(iDebugLevel);
#endif

    SmFilletVertex * pMate = pFV->GetMate(0);
    if (pMate) {
        SmFilletEdge * pRail1 = pFG->GetRail(1-lRailIndex);
        SmFace * pFace1 = pRail1->GetOriginalFace(); NER(pFace1);
        SmPoint2d sUV1 = rTsectPnt.UVPos(1-lRailIndex);
        SmPoint3d sPnt1;
        SER(pFace1->GetSurface()->EvaluatePoint(sUV1,sPnt1));
        pMate->SetOriginalUV(sUV1);
        pMate->SetPoint(sPnt1);
        pMate->SetStatus(SM_FIL_PROCESSED);
        SmFilletVertexuse * pMateVU = pMate->GetVUAtRailEnd(pFG);
        SmTsectPnt & rTsectPnt1 = pMateVU->GetTsectPnt();
        rTsectPnt1.UVPos(1-lRailIndex) = sUV1;

#ifdef SM_DEBUG_CODE
        if ( iDebugLevel > 0 ) {
            smgfx_SetLook(4,6, 0,1,1); sPnt1.Draw(); sm_GraphicsLoop();
            smgfx_SetLook(1,2, 1,1,0); pFace1->GetSurface()->DrawUV(1,1); sm_GraphicsLoop();
            smgfx_SetLook(3,4, 0,1,0); if (pRail1->GetCurve()) { pRail1->Draw(); } sm_GraphicsLoop();
            sm_GraphicsLoop();
        }
#endif
    }

    return SM_SUCCESS;

} // end sm_CalcRailXEdge

/****************************************************************
PURPOSE: Test if the input edge bounds two tangentially surfaces
            If so, are two surfaces in opposite orientations

NOTES:
****************************************************************/
static SmBoolean sm_TestManifoldTangentRollover
  (SmEdge    * pE,                 // in :
   SmBoolean & rbOppositeNormals,  // in :
   double      dTangencyTolRad,    // in :
   int         iDebugLevel )       // in :
{
    SmEdgeuse * pEU = pE->GetPrimaryEdgeuse();
    double dTestParam = pE->GetInterval().Evaluate(0.34);
    SmPoint3d sPnt, sPnt2;
    SmVector3d sBinVec2, sBinVec, sFaceuseNormal, sFaceuseNormal2;
    SE(pEU->EvaluateBinormal(dTestParam,
        FALSE,sPnt,sBinVec,NULL,&sFaceuseNormal));
    if (pEU->GetFaceuse()->GetOrientation() == SM_OT_OPPOSITE) {
        sFaceuseNormal = - sFaceuseNormal;
    }
    SmEdgeuse * pOtherEU = pEU->GetRadial()->GetMate();
    SE(pOtherEU->EvaluateBinormal(dTestParam,
        FALSE,sPnt2,sBinVec2,NULL,&sFaceuseNormal2));
    if (pOtherEU->GetFaceuse()->GetOrientation() == SM_OT_OPPOSITE) {
        sFaceuseNormal2 = - sFaceuseNormal2;
    }
#ifdef SM_DEBUG_CODE
    if ( iDebugLevel > 0 ) {
        smgfx_SetLook(2,2, 1,0,0); pEU->Draw(); sm_GraphicsLoop();
        smgfx_SetLook(2,2, 0,1,0); pOtherEU->Draw(); sm_GraphicsLoop();
        smgfx_SetLook(1,1, 0,0,0); pEU->GetFace()->DrawUV(2,2); sm_GraphicsLoop();
        smgfx_SetLook(1,1, 0,0,1); pEU->GetRadial()->GetFace()->DrawUV(2,2); sm_GraphicsLoop();
        smgfx_SetLook(2,2, 1,0,0); sFaceuseNormal .Draw(&sPnt); sm_GraphicsLoop();
        smgfx_SetLook(2,2, 0,1,0); sFaceuseNormal2.Draw(&sPnt2); sm_GraphicsLoop();
        sm_GraphicsLoop();
    }
#else
    SM_REF1(iDebugLevel);
#endif
    rbOppositeNormals = FALSE;
    if (sFaceuseNormal.Dot(sFaceuseNormal2) < 0.0)
        rbOppositeNormals = TRUE;

    // Test tangency.  Note, we're comparing vectors that would be perpendicular
    // when the Faces are tangent.
    double dDotLimit = smos_Cosine( SM_PI/2 - dTangencyTolRad );  // [B683]

    if ( smos_Fabs( sBinVec2.Dot( sFaceuseNormal )) > dDotLimit )
      { return FALSE; }

    return TRUE;
} // end sm_TestManifoldTangentRollover

/*******************************************************************//**
PURPOSE: SolidCraft change. Whether the two faces at the manifold edge
    pEdge meet tangentially there (sm_TestManifoldTangentRollover).
***********************************************************************/
SmBoolean SmFilletSolver::IsTangentEdge
  (SmEdge * pEdge,             // in : a manifold edge
   double   dTangencyTolRad)   // in : the angle tolerance
{
  SmBoolean bOpposite = FALSE;
  return ( pEdge != NULL && pEdge->IsManifold() && sm_TestManifoldTangentRollover( pEdge, bOpposite, dTangencyTolRad, 0 )) ? TRUE : FALSE;
}


/****************************************************************
PURPOSE:  Split a rail edge into two topologically

NOTES:
****************************************************************/
static SmStatus sm_SplitRail       
  (SmFilletGeom    * pFG,          // in :
   SmFilletEdge    * pRailToSplit, // in : Will remain as the first "half" of the rail
   SmFilletVertex *& rpNewVert,    // out:
   SmFilletEdge   *& rpNewRail)    // out: Newly created 'second half' of the rail
{
  // locals
  SmFilletBrep   * pPseudoBrep = pFG->GetFilletSolver()->GetFilletExecutive()->GetPseudoBrep();
  SmFilletVertex * pStartVert  = (SmFilletVertex*)pRailToSplit->GetVertex();
  SmFilletVertex * pEndVert    = (SmFilletVertex*)pRailToSplit->GetOtherVertex(pStartVert);
  SmEdgeuse      * pEU         = pRailToSplit->GetPrimaryEdgeuse();
  SmVertexuse    * pVU         = pEU->GetVertexuse();

  if (pVU->GetVertex() != pEndVert) 
    { 
      pEU = pEU->GetMate();
      pVU = pEU->GetVertexuse();
    }

  if(   pStartVert == pEndVert 
     && pEU->GetOrientation() != SM_OT_OPPOSITE) 
    {
      // Typically, this might happen when we had closed-rail,
      // Get the other one instead
      pEU = pEU->GetMate();
      pVU = pEU->GetVertexuse();
    }

  // Create a new vertex in the middle of the rail
  ((SmFilletEdgeuse*)pEU)->SetUVCurve(NULL, TRUE) ; // TRUE = delete preExisting m_pUVTrimCurve

  SER( pEndVert->Remove(pVU) );
  rpNewVert = new (pPseudoBrep) SmFilletVertex();
  SER(rpNewVert->PostInsert(pVU));

  // Make NewRail from NewVert to EndVert - has no Curve or UVTrimCurve data
  SER(SmFilletEdge::MakeFilletEdge(pPseudoBrep,  // in : target Brep to receive new topology objects
                                   rpNewVert,    // in : start of new FilletEdge
                                   pEndVert,     // in : end   of new FilletEdge
                                   rpNewRail));  // out: newly allocated FilletEdge - no Curve or UVTrimCurve data
                                                 // in : NewFilletEdge's filletCorner, gets stored in rpNewFilletEdge->m_cpCorner
                                                 //      NULL to ignore, default:[NULL]
  rpNewRail->SetFilletEdgeType(SM_FE_RAIL);
  rpNewRail->SetOriginalFace(pRailToSplit->GetOriginalFace());

  SmFilletEdgeuse * pNewEU = (SmFilletEdgeuse*)rpNewRail->GetPrimaryEdgeuse();
  if (pNewEU->GetVertexuse()->GetVertex() != pEndVert) 
    {
      pNewEU = (SmFilletEdgeuse*)pNewEU->GetMate();
    }
  pNewEU->SetFilletGeom(pFG);

  SmCurve * pOrigCurve = pRailToSplit->GetCurve();
  SmCurve * pNewCurve = NULL;
  SER(pOrigCurve->Copy(*pOrigCurve->GetContext(),pNewCurve));
  rpNewRail->SetCurve(pNewCurve, TRUE ) ; // TRUE = delete preExisting Edge->Curve - expect none here
                                          // side effect: delete current pSurviveEdge->UVTrimCurve
  pNewCurve->SetOwner(rpNewRail) ;
#ifdef SM_USE_NEWTOL      
  SM_NEWTOL_LINE rpNewRail->ClearLocalZoneTol3d( ) ;  // what's going on with the MaxGaps here?
#else // SM_USE_OLDTOL 
  SM_OLDTOL_LINE rpNewRail->SetTolerance(pRailToSplit->GetTolerance());
#endif // SM_USE_OLDTOL

  return SM_SUCCESS;
} // end sm_SplitRail

/****************************************************************
PURPOSE: Determine the convexity of the side-edge of the original brep.
            Tangent rollover should not happen when this routine is called.

NOTES:
****************************************************************/
SmBoolean sm_TestSideEdgeConvexity
  (SmEdge * pSideEdge,
   SmBoolean bIsTangentialSideEdge)
{
    SmTArray<SmEdgeuse*> sEUs;
    (void)pSideEdge->GetEdgeuses(sEUs); // SolidCraft: a list getter; it fails only out of memory, and the list stays empty
    for (ULONG ii=0; ii<sEUs.GetSize(); ii++) {
        SmEdgeuse * pEU = sEUs[ii];
        SmShell *pShell = pEU->GetShell();
        SmRegion * pInfiniteReg = pShell->GetBrep()->GetInfiniteRegion();
        if (pShell->GetRegion() == pInfiniteReg) {
            pEU = pEU->GetMate();
        }

        SmBoolean bStepOff = FALSE;
        if (bIsTangentialSideEdge) {
            bStepOff = TRUE;
        }
        if (pEU->IsConvexRadialSector(10,bStepOff)) {
            return TRUE;
        }
    } // end for each Edgeuse

    return FALSE;
} // end sm_TestSideEdgeConvexity

/****************************************************************
PURPOSE: Find the side-edge which intersected with rail

NOTES:
****************************************************************/
static SmStatus sm_FindRolloverSideEdge
  (SmCurveClassification * pCC,                      // in :
   ULONG                 & rlCurrIndx,               // in :
   SmFilletEdge          * pRail,                    // in :
   SmCurve               * pOrigRailCurve,           // NotUsed: in :
   SmCurve               * pOrigRailUVCurve,         // in :
   SmSurface             * pOrigFillet,              // in :
   SmBoolean               bConvexFillet,            // in :
   SmEdge               *& rpSideEdge,               // out:
   SmBoolean             & rbIsTangentialSideEdge,   // out:
   SmBoolean             & rbTrimSideFaces,          // out:
   SmBoolean             & rbOppositeNormals,        // out:
   SmEdgeuse            *& rpRolloverCornerEU,       // out:
   SmVertex             *& rpRolloverCornerV,        // out:
   double                  dTangencyTolRad,          // in :
   int                     iDebugLevel)              // in :
{
  SM_REF1(pOrigRailCurve) ; 
  SmEdge * pSideEdge = NULL;
  ULONG lTotalIntervals = pCC->GetSize();
  SmFilletVertex * pStartV = (SmFilletVertex*) pRail->GetVertex();
  SmFilletVertex * pEndV = (SmFilletVertex*) pRail->GetOtherVertex(pStartV);
  SmFace * pOrigFace = pRail->GetOriginalFace(); NER(pOrigFace);

  rbTrimSideFaces = FALSE;
  while (TRUE)
  {
      SmCurveInterval & rIvl = (*pCC)[rlCurrIndx];

      if (    rlCurrIndx == 0
          && !pStartV->IsProcessed()
          &&  rIvl.m_vStart.GetPointClass() == SM_PC_UNKNOWN )
      {
          rlCurrIndx++;  // Skip this interval, move on to next.
          continue;
      }

      SmPointClassification sStartPC = rIvl.m_vStart;
      SmPointClassification sEndPC = rIvl.m_vEnd;
      if ( rlCurrIndx < lTotalIntervals-1 )
      {
          if ( rIvl.m_vEnd.GetPointClass() != SM_PC_EDGE )
          {
              SER(SM_ERR); // Unknown case
          }
          pSideEdge = (SmEdge*)rIvl.m_vEnd.GetObject();

          // Check whether the side edge we hit is a tangent edge,
          // and whether it's convex.
          rbIsTangentialSideEdge = sm_TestManifoldTangentRollover(pSideEdge,
                                                                  rbOppositeNormals,
                                                                  dTangencyTolRad,
                                                                  iDebugLevel );
          SmBoolean bConvexSideEdge = sm_TestSideEdgeConvexity(pSideEdge,
                                                               rbIsTangentialSideEdge );

          // If the filleted edge and the side edge that the rail hit
          // have the same convexity, then the fillet surface will
          // crash into the adjoining face, and we'll trim it.
          // Otherwise, the fillet surface will be hanging over a cliff,
          // and we'll extend the side face (or add a new little face)
          // to intersect the fillet.
          if ( bConvexFillet == bConvexSideEdge )
          {
              // Trimming process
              rbTrimSideFaces = TRUE;
              rlCurrIndx++;
              continue;
          }

          if (sStartPC.GetPointClass() == SM_PC_EDGE &&
              sEndPC.GetPointClass() == SM_PC_EDGE &&
              rIvl.m_vMid.GetPointClass() == SM_PC_UNKNOWN) {
              SmEdge * pE1 = (SmEdge*)sStartPC.GetObject();
              SmEdge * pE2 = (SmEdge*)sEndPC.GetObject();
              SmEdgeuse * pEUOnOrigFace = pE1->GetPrimaryEdgeuse();
              if (pEUOnOrigFace->GetFace() != pOrigFace) {
                  pEUOnOrigFace = pEUOnOrigFace->GetRadial();
#ifdef SM_DEBUG_CODE
                  if ( iDebugLevel > 0 ) {
                      smgfx_SetColor(1,0,0);
                      pE1->Draw();
                      sm_GraphicsLoop();
                      pEUOnOrigFace->Draw();
                      sm_GraphicsLoop();
                      if (rpRolloverCornerEU) {
                          smgfx_SetColor(0,1,0);
                          rpRolloverCornerEU->Draw();
                          sm_GraphicsLoop();
                      }
                  }
#endif
              }

              double dStartParam = sStartPC.GetTParam();
              double dEndParam = sEndPC.GetTParam();
              SmCurve * pE1Curve = pE1->GetCurve();
              SmVector3d sPV[2];
              SER(pE1Curve->Evaluate(dStartParam,1,TRUE,sPV));
              SmVector3d sTangent = sPV[1];
              // Compute 'binormal' of rail curve
              double dRailParam = rIvl.m_vInterval.GetMin();
              SmVector3d sPnt,sPVV[2][2];
              SER(pOrigRailUVCurve->EvaluatePoint(dRailParam,sPnt));
              SmPoint2d sUV(sPnt.x,sPnt.y);
              SER(pOrigFillet->Evaluate(sUV,1,1,TRUE,TRUE,TRUE,sPVV[0]));
              SmVector3d sBinVec = sPVV[0][1];
              SmPoint2d sCenterUV = pOrigFillet->GetNaturalUVDomain().Evaluate(0.5,0.5);
              if (sUV.y > sCenterUV.y) {
                  sBinVec = - sBinVec;
              }

              double dAngle;
              SER(sBinVec.AngleBetween(sTangent,dAngle));
              SmOrientType eEUOrient = (dAngle<SM_PI/2.0) ? SM_OT_SAME : SM_OT_OPPOSITE;

              if (rpRolloverCornerEU) {
                  SmEdgeuse * pRadialEU = rpRolloverCornerEU->GetRadial();
                  SmEdgeuse * pEU = pRadialEU->GetCWEdgeuse();
                  if (rbIsTangentialSideEdge) {
                      sm_TestManifoldTangentRollover(pEU->GetEdge(),rbOppositeNormals, dTangencyTolRad, iDebugLevel);
                  }
                  if (pEU->GetEdge() != pE2) {
                      rpRolloverCornerV = pEU->GetVertexuse()->GetVertex();
                      rpRolloverCornerEU = pEU->GetCWEdgeuse();
                      rpSideEdge = pEU->GetEdge();
                      return SM_SUCCESS;
                  }
                  else {
                      rpRolloverCornerEU = NULL;
                      rpRolloverCornerV = NULL;
                  }
              }
              else
              {
                  if ( pEUOnOrigFace->GetOrientation() != eEUOrient )
                  {
                      pEUOnOrigFace = pEUOnOrigFace->GetMate();
                  }
                  SmEdgeuse * pEU = pEUOnOrigFace->GetRadial();
                  //SmVertex * pCornerV = NULL;
                  SmBoolean bRolloverCorner = (pE1 != pE2);
                  if ( pE1 == pE2 && pE1->IsClosed() )
                  {
                      if ((dStartParam < dEndParam && eEUOrient == SM_OT_OPPOSITE) ||
                          (dStartParam > dEndParam && eEUOrient == SM_OT_SAME)) {
                          bRolloverCorner = TRUE;
                      }
                  }
                  if ( bRolloverCorner )
                  {
                      rbIsTangentialSideEdge = sm_TestManifoldTangentRollover(
                          pEU->GetEdge(), rbOppositeNormals, dTangencyTolRad, iDebugLevel );
                      rpRolloverCornerV  = pEU->GetVertexuse()->GetVertex();
                      rpRolloverCornerEU = pEU->GetCWEdgeuse();
                      if ( rbIsTangentialSideEdge )
                      {
                          rpSideEdge = rpRolloverCornerEU->GetEdge();
                      }
                      else
                      {
                          rpSideEdge = pEU->GetEdge();
                      }
                      return SM_SUCCESS;
                  }
#ifdef SM_DEBUG_CODE
                  if ( iDebugLevel > 0 ) {
                      if ( rpSideEdge)
                        { smgfx_SetLook(2,4, 1,0,0); rpSideEdge->Draw(); sm_GraphicsLoop(); }
                      if (rpRolloverCornerV)
                        { smgfx_SetLook(2,4, 0,0,1); rpRolloverCornerV->Draw(); sm_GraphicsLoop(); }
                      if (rpRolloverCornerEU)
                        { smgfx_SetLook(2,4, 0,1,0); rpRolloverCornerEU->Draw(); sm_GraphicsLoop(); }
                      sm_GraphicsLoop();
                  }
#endif
              }
          }

          // If this is the one before the very last segment,
          // do the following check
          if ( rlCurrIndx+1 == lTotalIntervals-1 )
          {
              SmCurveInterval & rNextIvl = (*pCC)[rlCurrIndx+1];
              if (    rNextIvl.m_vEnd.GetPointClass() == SM_PC_UNKNOWN
                   && !pEndV->IsProcessed())
              {
                  rlCurrIndx++;
                  return SM_SUCCESS;
              }
          }
          rpRolloverCornerEU = NULL;
          rpSideEdge = pSideEdge;
          return SM_SUCCESS;
      }
      else // Process the very last segment
      {
          if ( rIvl.m_vStart.GetPointClass() != SM_PC_EDGE )
          {
              SER(SM_ERR); // Unknown case
          }
          if (   rIvl.m_vEnd.GetPointClass() == SM_PC_UNKNOWN
              && !pEndV->IsProcessed() )
          {
              return SM_SUCCESS;
          }
          pSideEdge = (SmEdge*)rIvl.m_vStart.GetObject();

          // Check whether the side edge we hit is tangent, or concave/convex.
          rbIsTangentialSideEdge = sm_TestManifoldTangentRollover(
              pSideEdge, rbOppositeNormals, dTangencyTolRad, iDebugLevel );
          SmBoolean bConvexSideEdge = sm_TestSideEdgeConvexity(
              pSideEdge, rbIsTangentialSideEdge );

          // Cliff or crash (comment above).
          if ( bConvexFillet == bConvexSideEdge )
          {
              // Trimming process
              rbTrimSideFaces = TRUE;
              return SM_SUCCESS;
          }
          rpSideEdge = pSideEdge;

          return SM_SUCCESS;
      }
  }

  return SM_SUCCESS;
} // end sm_FindRolloverSideEdge

/****************************************************************
PURPOSE: Computes the 3D Geometry (Curve and endVertexPoint) for pSupportEdge

NOTES: 1. Updates pSupportEdge->Curve = Trimmed XSectCurve(FilletSurf1, FilletSurf2)
                  pSupportEdge->EndVertex point and status to SM_FIL_PROCESSED

       2. preExisting pSupportEdge->Curve is deleted
       3. preExisting pSupportEdge->UVTrimCurves are deleted
****************************************************************/
static SmStatus sm_ComputeSupportEdge
  (const SmContext  & crContext,         // in : 
   SmFilletEdge     * pSupportEdge,      // in : Gets new 3DCurve (any old 3DCurve or UVTrimCurves are deleted)
   // gwc: I think pSideSurface no longer has to be a BBSplineSurface - replace 1 line
   // rm : SmBSplineSurface * pFilletSurface,    // in :
   SM_FILLETSURF_TYPE * pFilletSurface,    // in :
   double             dThisApproxTol3d,  // in :
   double             dAngleTolerance,   // in :
   int                iDebugLevel )      // in :
{
  // locals
  SmFilletVertex * pEdgeStartV = (SmFilletVertex*)pSupportEdge->GetVertex();

  // pEdgeStartV should have been calculated
  if (!pEdgeStartV->IsProcessed()) 
    { SER(SM_ERR); }

  SmPoint3d          sStartPnt = pEdgeStartV->GetPoint();
  SmFilletEdgeuse  * pPrimEU   = (SmFilletEdgeuse*)pSupportEdge->GetPrimaryEdgeuse();
  SmFilletEdgeuse  * pMateEU   = (SmFilletEdgeuse*)pPrimEU->GetMate();
  SmFilletGeom     * pFG1      = pPrimEU->GetFilletGeom(); NER(pFG1);
  SmFilletGeom     * pFG2      = pMateEU->GetFilletGeom(); NER(pFG2);
  
  // gwc: I think pSideSurface no longer has to be a BBSplineSurface - replace 2 line
  // rm : SmBSplineSurface * pSurf1    = pFG1->GetFilletSurface(); NER(pSurf1);
  // rm : SmBSplineSurface * pSurf2    = pFG2->GetFilletSurface(); NER(pSurf2);
  SM_FILLETSURF_TYPE * pSurf1    = pFG1->GetFilletSurface(); NER(pSurf1);
  SM_FILLETSURF_TYPE * pSurf2    = pFG2->GetFilletSurface(); NER(pSurf2);

  SmBSplineCurve   * p3DCurve  = NULL;
  SmBSplineCurve   * pUVCurve1 = NULL;
  SmBSplineCurve   * pUVCurve2 = NULL;
  if(   SM_SUCCESS != sm_SrfSrfIntersection(crContext,          // in : context for new object construction
                                            pEdgeStartV,        // in : 1st point known to be on intersection curve
                                            NULL,               // in : 2nd point known to be on intersection curve,
                                                                //      NULL to ignore
                                            NULL,               // in : expected general direction of intersection curve from 1st point
                                            NULL,               // in   expected intersection end direction
                                            pSurf1,             // in : 1st intersecting surface
                                            pSurf2,             // in : 2nd intersecting surface
                                            dThisApproxTol3d,  // in : max allowed distance between xSect Curve and surfaces
                                            dAngleTolerance,    // in : max allowed angle between consecutive xSect curve segment tangents
                                            p3DCurve,           // out: 3d intersection curve
                                            pUVCurve1,          // out: associated UVTrimCurve on 1st surface
                                            pUVCurve2,          // out: associated UVTrimCurve on 2nd surface
                                            FALSE,              // in : FALSE= try cheap sm_TwoPntsIntersection() before
                                                                //             resorting to expensive general surf/surf xSect solver
                                            NULL,               // in : reference point
                                            iDebugLevel)        // in : iDebugLevel, 0 = No Debug output
     || !p3DCurve)                                              
    {                                                           
      pSupportEdge->SetStatus(SM_FIL_SURF_INT_FAILURE);         
      SER(SM_ERR);
    }

  // not using the rail curves - delete them
  SM_ASSERT(pUVCurve1 != NULL) ; delete pUVCurve1 ; pUVCurve1 = NULL ;
  SM_ASSERT(pUVCurve2 != NULL) ; delete pUVCurve2 ; pUVCurve2 = NULL ;

  // Intersect the curve with the fillet for the pNewV's geomtry
  SmExtent1d      sTrimIvl = p3DCurve->GetNaturalInterval();
  SmSolution      aData[10];
  SmSolutionArray sSolutions(10,aData);
  SER(pFilletSurface->GlobalCurveIntersect(pFilletSurface->GetNaturalUVDomain(),
                                   *p3DCurve,
                                    sTrimIvl,
                                    dThisApproxTol3d,
                                    sSolutions));
  if (sSolutions.GetSize() < 1) 
    { SER(SM_ERR); }

  double dT = sSolutions[0].m_vStart[0];
  SmPoint3d sEndPnt;
  SER(p3DCurve->EvaluatePoint(dT,sEndPnt));

  // Update end vertex
  SmFilletVertex * pEdgeEndV = (SmFilletVertex*)pSupportEdge->GetOtherVertex(pEdgeStartV);
  pEdgeEndV->SetPoint(sEndPnt);
  pEdgeEndV->SetStatus(SM_FIL_PROCESSED);

#ifdef SM_DEBUG_CODE
  if ( iDebugLevel > 0 ) 
    {
      smgfx_SetLook(1,2, 1,0,0); sStartPnt.Draw(); sm_GraphicsLoop();
      smgfx_SetLook(1,2, 0,0,1); sEndPnt.Draw(); sm_GraphicsLoop();
      smgfx_SetLook(1,2, 0,1,0); if (p3DCurve) p3DCurve->DrawWDeriv(sTrimIvl,0); sm_GraphicsLoop();
      smgfx_SetLook(1,2, 0,1,1); pSurf1->DrawUV(4,4); sm_GraphicsLoop();
      smgfx_SetLook(1,2,.5,1,0); pSurf2->DrawUV(4,4); sm_GraphicsLoop();
      sm_GraphicsLoop();
    }
#endif // SM_DEBUG_CODE

  // Trim the 'support' edge by the end points
  SmBoolean bSuccess;
  double    dParam1 = 0.0, dParam2 = 0.0, dDist = 0.0;
  SER(p3DCurve->DropPoint(sTrimIvl,          // in : target curve allowed domain
                          sStartPnt,         // in : Point to drop to curve
                          NULL,              // in : For Snapping, NULL=NoSnap, Vec pointing to CrvInside from desired EndPt for snapping.
                                             //      when DropPt is within dDistTol in cpVecToCrvInside dir of EndPt, Snap to Ivl EndPt.
                                             //      Vec handles ambiguities and makes sure that curves are not just touching at the ends.
                          dThisApproxTol3d,  // in : Used for SnapDist, and SM_SO_INTERSECT sXSectTol3d distance.
                                             //      sXSectTol3d = dDistTol ;      SM_SO_INTERSECT 
                                             //      sXSectTol3d = SM_BIG_DOUBLE ; SM_SO_MINIMIZE, SM_SO_NORMALIZE
                                             //      sSnapTol3d  = dDistTol ;      SM_SO_MINIMIZE, SM_SO_NORMALIZE, SM_SO_INTERSECT    
                          NULL,              // in : If given, call cheap LocalPointSolve() else call expensive GlobalPointSolve()
                          bSuccess,          // out: TRUE = found a drop point
                          dParam1,           // out: found drop curve param
                          dDist)) ;          // out: found drop distance
                                             // in : SM_SO_MINIMIZE = allow nonNormal drops near endPoints
                                             //      SM_SO_INTERSECT= like Minimize but point must be within dDistTol
                                             //      SM_SO_NORMALIZE= exclude nonNormal drops near endPoints
                                             //      default:[SM_SO_MINIMIZE] to preserve original behavior
  if (bSuccess) 
    {
      SER(p3DCurve->DropPoint(sTrimIvl,          // in : target curve allowed domain
                              sEndPnt,           // in : Point to drop to curve
                              NULL,              // in : For Snapping, NULL=NoSnap, Vec pointing to CrvInside from desired EndPt for snapping.
                                                 //      when DropPt is within dDistTol in cpVecToCrvInside dir of EndPt, Snap to Ivl EndPt.
                                                 //      Vec handles ambiguities and makes sure that curves are not just touching at the ends.
                              dThisApproxTol3d,  // in : Used for SnapDist, and SM_SO_INTERSECT sXSectTol3d distance.
                                                 //      sXSectTol3d = dDistTol ;      SM_SO_INTERSECT 
                                                 //      sXSectTol3d = SM_BIG_DOUBLE ; SM_SO_MINIMIZE, SM_SO_NORMALIZE
                                                 //      sSnapTol3d  = dDistTol ;      SM_SO_MINIMIZE, SM_SO_NORMALIZE, SM_SO_INTERSECT    
                              NULL,              // in : If given, call cheap LocalPointSolve() else call expensive GlobalPointSolve()
                              bSuccess,          // out: TRUE = found a drop point
                              dParam2,           // out: found drop curve param
                              dDist)) ;          // out: found drop distance
                                                 // in : SM_SO_MINIMIZE = allow nonNormal drops near endPoints
                                                 //      SM_SO_INTERSECT= like Minimize but point must be within dDistTol
                                                 //      SM_SO_NORMALIZE= exclude nonNormal drops near endPoints
                                                 //      default:[SM_SO_MINIMIZE] to preserve original behavior
    }

  if (!bSuccess) 
    {
      SM_ASSERT(p3DCurve != NULL) ; delete p3DCurve ; p3DCurve = NULL ;
      SER(SM_ERR);
    }

  SmBoolean bReverseCurves = FALSE;
  if (dParam1 > dParam2) { bReverseCurves = TRUE;
                           (void)sTrimIvl.SetMinMax(dParam2,dParam1); // SolidCraft: ordered by the test, it cannot fail
                         }
  else                   { (void)sTrimIvl.SetMinMax(dParam1,dParam2); // SolidCraft: ordered by the test, it cannot fail
                         }

  SER(p3DCurve->Trim(sTrimIvl));  // may snap sIvl by tol to existing knots
  if (bReverseCurves) 
    {
      SER(p3DCurve->ReverseParameterization(sTrimIvl,sTrimIvl));
    }

  if (p3DCurve) 
    {
      pSupportEdge->SetCurve(p3DCurve, TRUE ) ; // TRUE = delete preExisting Edge->Curve - expect none here
                                                // side effect: delete current pSurviveEdge->UVTrimCurve
      p3DCurve->SetOwner(pSupportEdge);
    }

  // all done
  return SM_SUCCESS;

} // end sm_ComputeSupportEdge

/*******************************************************************//**
    END - Static functions
***********************************************************************/


/*******************************************************************//**
PURPOSE: Create side fillet geom for cliff rollover. The surface of
    this fillet geom is an extension of the side face where the ball rollover

NOTES: If bTrimSideFaces is TRUE, the rollover is a trimming process
    rather than 'adding' extended side-face(or side-FG)
***********************************************************************/
SmStatus SmFilletSolver::CreateCliffRollOverSideFG
  (SmFilletEdge  * pRailOnCliff,     // in : FilletRail on Cliff
   SmEdge        * pSideEdge,        // in : 
   SmFace        * pSideFace,        // in : 
   SmBoolean       bTrimSideFaces,   // in : TRUE = Side Faces set - set pRailOnCliff->OrigFace = pSideFace as is
                                     //                            - set CornerVert->pos        = CrossSectionCurve/SideFace XSectPt
                                     //                            - adjust CrossSection->Curve and UVTrimCurve intervals
                                     //      FALSE= Extend copy of pSideFace - set pRailOnCliff->OrigFace = pSideFace as is
                                     //                            - set CornerVert->pos = CrossSectionCurve/SideFace XSectPt
                                     //                            - adjust CrossSection->Curve and UVTrimCurve intervals
                                     //                            - build new SplitEdge if needed
                                     //                            - build new SideFG: Curve        = Trimmed Copy of pSideEdge
                                     //                                                UVTrimCurves = NULL
                                     //                                                OrigEdge     = pSideEdge
                                     //                                                OrigFace     = pSideFace
   double          dExtensionDist,   // in : Amount pSideFace->SurfaceCopy is extended, if extended
   SmFilletGeom *& rpSideFG)         // out: The new FilletGeom when bTrimSideFaces == FALSE
{
  // locals
  SmFilletVertex * pRailStartV    = (SmFilletVertex*)pRailOnCliff->GetVertex();
  SmFilletVertex * pRailEndV      = (SmFilletVertex*)pRailOnCliff->GetOtherVertex(pRailStartV);
  SmFilletBrep   * pPseudoBrep    = this->GetFilletExecutive()->GetPseudoBrep();
  SmFilletVertex * pAdjCornerVert =   (pRailStartV->GetFilletCorner() != NULL) ? pRailStartV // CornerVert geometry
                                    : (pRailEndV->GetFilletCorner()   != NULL) ? pRailEndV 
                                                                               : NULL;
int iDebugLevel = 0;
#ifdef SM_DEBUG_CODE
iDebugLevel = DebugLevel();
SmBoolean bDebugMe0=FALSE;
    if (bDebugMe0) 
      {
        SmEdge * pOrigEdge   = pRailOnCliff->GetOriginalEdge() ;
        SmFace * pOrigFace   = pRailOnCliff->GetOriginalFace() ;
        //SmBrep * pPseudoBrep = pRailOnCliff->GetBrep() ;
        SmBrep * pTargetBrep =   pOrigEdge ? pOrigEdge->GetBrep() 
                               : pOrigFace ? pOrigFace->GetBrep()
                               : NULL ;

        smgfx_Erase();
        smgfx_SetLook(1,2, 0,0,1) ; if(pTargetBrep) pTargetBrep->Draw(TRUE) ; sm_GraphicsLoop() ;
        smgfx_SetLook(5,6, 1,0,0); if (pRailStartV->IsProcessed()) pRailStartV->Draw(); sm_GraphicsLoop();
        smgfx_SetLook(5,6, 0,0,1); if (pRailEndV->IsProcessed()) pRailEndV->Draw(); sm_GraphicsLoop(); // It's OK if pRailEndV is not processed
        smgfx_SetLook(3,4, 0,1,0); pRailOnCliff->Draw(); sm_GraphicsLoop();
        smgfx_SetLook(3,4,.5,0,1); pSideEdge->Draw(); sm_GraphicsLoop();
        sm_GraphicsLoop();
      }
#endif // SM_DEBUG_CODE

  // Pick the SideSurface to extend - either pSideFace->GetSurface() or its copy
  
  // gwc: I think pSideSurface no longer has to be a BBSplineSurface - replace 1 line
  // rm :SmBSplineSurface * pSideSurface = NULL;
  SM_FILLETSURF_TYPE * pSideSurface = NULL;
  SmSurface * pAdjSurface  = pSideFace->GetSurface();
  SmObjDelete sCleanSurf;
  if ( bTrimSideFaces ) { // gwc: I think pSideSurface no longer has to be a BBSplineSurface - replace 1 line
                          // rm : pSideSurface = SM_CAST_PTR( SmBSplineSurface, pAdjSurface );
                          pSideSurface = SM_CAST_FILLETSURF_PTR(SM_FILLETSURF_TYPE, pAdjSurface) ; 
                          pRailOnCliff->SetOriginalFace(pSideFace);
                        }
  else                { NER(pAdjSurface);
                        SmSurface * pNewSurface = NULL ;
                        SmSurface * pCopySurface = NULL;
                        SER(pAdjSurface->CreateExtendedSurface(m_crContext,
                                                               dExtensionDist,
                                                               SM_CT_G1,
                                                               pNewSurface));
                        sCleanSurf.SetObj(pNewSurface);

                        // Copy pNewSurface, when possible as an analytic surface
                        SER(pNewSurface->CopyAndAddAnalytics(m_crContext, pCopySurface));
                          // gwc: I think pSideSurface no longer has to be a BBSplineSurface - replace 1 line
                          // rm : pSideSurface = SM_CAST_PTR(SmBSplineSurface, pCopySurface);
                          pSideSurface = SM_CAST_FILLETSURF_PTR(SM_FILLETSURF_TYPE, pCopySurface) ; 
                        pRailOnCliff->SetOriginalFace(NULL);
                      }

  // when given the CornerVert geometry - set CornerVert->Pos = CrossSection/ExtendedSideSurf XSect
  //                                    - update CrossSection->Curve and UVTrimCurve intervals to run 
  //                                      between current CrossSection->OtherVert and new XSectPoint
  if (pAdjCornerVert) 
    {
      SmTArray<SmEdge*> sCornerEdges;
      pAdjCornerVert->GetEdges(sCornerEdges);

      // for every Corner edge - find the SM_FE_CROSS_SECTION and SM_FE_SPLIT_FACE edges
      SmFilletEdge * pSplitEdge    = NULL;
      SmFilletEdge * pCrossSection = NULL;

      // search list for CROSS_SECTION and SPLIT_FACE FilletEdges - gwc?:does it matter that this only saves the last such edges found?
      for (ULONG ii=0; ii<sCornerEdges.GetSize(); ii++) 
        {
          SmFilletEdge * pFE = (SmFilletEdge*)sCornerEdges[ii];

          if      (pFE->GetFilletEdgeType() == SM_FE_CROSS_SECTION) { pCrossSection = pFE; }
          else if (pFE->GetFilletEdgeType() == SM_FE_SPLIT_FACE)    { pSplitEdge    = pFE; }
        }
      NER(pCrossSection);

      // CrossSection Curve locals
      SmCurve   * pCurve   = pCrossSection->GetCurve(); NER(pCurve);
      SmExtent1d  sTrimIvl = pCurve->GetNaturalInterval();

      // next: 1. Set CornerVert geom (i.e. pAdjCornerVert) = intersection of cross-section with the extended side surface
      //       2. Trim cross-section edge to new CornerVert position

      // Intersect CrossSection with extended SideSurface
      SmSolution      aData[10];
      SmSolutionArray sSolutions(10,aData);
      SER(pSideSurface->GlobalCurveIntersect(pSideSurface->GetNaturalUVDomain(),
                                            *pCurve,
                                             sTrimIvl,
                                             m_dThisApproxTol3d,
                                             sSolutions));
      if (sSolutions.GetSize() < 1) SER(SM_ERR);

      // get CrossSection/SideFace intersection point
      double dParam1 = sSolutions[0].m_vStart[0];
      SmPoint3d sCornerVertPnt;
      SER(pCurve->EvaluatePoint(dParam1,sCornerVertPnt));

      // Update corner vert geom (i.e. pAdjCornerVert) to CrossSection/SideFace intersection point
      pAdjCornerVert->SetPoint(sCornerVertPnt);

      // Drop CrossSection OtherVertex onto CrossSection->Curve - to define CrossSection->Ivl[dParam1 dParam2]
      SmVertex * pOtherV = pCrossSection->GetOtherVertex(pAdjCornerVert);
      double     dParam2, dDist;
      SmBoolean  bSuccess;
      SER(pCurve->DropPoint(sTrimIvl,             // in : target curve allowed domain
                            pOtherV->GetPoint(),  // in : Point to drop to curve
                            NULL,                 // in : For Snapping, NULL=NoSnap, Vec pointing to CrvInside from desired EndPt for snapping.
                                                  //      when DropPt is within dDistTol in cpVecToCrvInside dir of EndPt, Snap to Ivl EndPt.
                                                  //      Vec handles ambiguities and makes sure that curves are not just touching at the ends.
                            m_dThisApproxTol3d,   // in : Used for SnapDist, and SM_SO_INTERSECT sXSectTol3d distance.
                                                  //      sXSectTol3d = dDistTol ;      SM_SO_INTERSECT 
                                                  //      sXSectTol3d = SM_BIG_DOUBLE ; SM_SO_MINIMIZE, SM_SO_NORMALIZE
                                                  //      sSnapTol3d  = dDistTol ;      SM_SO_MINIMIZE, SM_SO_NORMALIZE, SM_SO_INTERSECT    
                            NULL,                 // in : If given, call cheap LocalPointSolve() else call expensive GlobalPointSolve()
                            bSuccess,             // out: TRUE = found a drop point
                            dParam2,              // out: found drop curve param
                            dDist)) ;             // out: found drop distance
                                                  // in : SM_SO_MINIMIZE = allow nonNormal drops near endPoints
                                                  //      SM_SO_INTERSECT= like Minimize but point must be within dDistTol
                                                  //      SM_SO_NORMALIZE= exclude nonNormal drops near endPoints
                                                  //      default:[SM_SO_MINIMIZE] to preserve original behavior
      // when DropPoint for Param2 succeeds - Trim existing CrosSection UVTrimCurve Ivls
      if (bSuccess) // GWC: should bSuccess != TRUE be a bug?
        {
          // order the desired CrossSection->Ivl params
          if (dParam1 > dParam2) { (void)sTrimIvl.SetMinMax(dParam2,dParam1); } // SolidCraft: ordered by the test, it cannot fail
          else                   { (void)sTrimIvl.SetMinMax(dParam1,dParam2); } // SolidCraft: ordered by the test, it cannot fail

          // Set CrossSection->Curve Ivl
          SER( pCurve->Trim(sTrimIvl) );  // may snap sIvl by tol to existing knots

          // Set CrossSection->UVTrimCurves Ivls
          SmEdgeuse      * pPrimEU   = pCrossSection->GetPrimaryEdgeuse();
          SmEdgeuse      * pMateEU   = pPrimEU->GetMate();
          SmBSplineCurve * pUVCurve1 = pPrimEU->GetUVTrimCurvePointer();
          SmBSplineCurve * pUVCurve2 = pMateEU->GetUVTrimCurvePointer();
          if (pUVCurve1) { SER( pUVCurve1->Trim(sTrimIvl) ); }   // may snap sIvl by tol to existing knots
          if (pUVCurve2) { SER( pUVCurve2->Trim(sTrimIvl) ); }   // may snap sIvl by tol to existing knots

        } // end DropPoint success check
      else
        {
          SM_ASSERT_MSG(bSuccess == TRUE, _T("SmFilletSolver::CreateCliffRollOverSideFG - CrossSection Vert failed to drop to CrossSection Curve and CrossSection Ivl was not updated - needs review - might be a bug")) ;
        }

      // when no SplitEdge was found and bTrimSideFaces = FALSE - make a new pSplitEdge
      if (!bTrimSideFaces && pSplitEdge == NULL) 
        {
          // Attach a support edge
          // Create a split FilletEdge between CrossSection/SideFace XSectPt and projection of that Pt to SideCurve
          SmFilletVertex * pEndVert = new (pPseudoBrep) SmFilletVertex();
          SER(SmFilletEdge::MakeFilletEdge(pPseudoBrep,    // in : target Brep to receive new topology objects
                                           pEndVert,       // in : start of new FilletEdge
                                           pAdjCornerVert, // in : end   of new FilletEdge
                                           pSplitEdge));   // out: newly allocated FilletEdge - no Curve or UVTrimCurve data
                                                           // in : The filletCorner that generated the need for this edge
          pSplitEdge->SetFilletEdgeType(SM_FE_SPLIT_FACE);

          // drop new CrossSection/SideFace XSectPT to pSideCurve
          SmCurve * pSideCurve = pSideEdge->GetCurve(); NER(pSideCurve);
          SER(pSideCurve->GlobalPointSolve(pSideEdge->GetInterval(),
                                           SM_SO_MINIMIZE,
                                           sCornerVertPnt,
                                           SM_EFF_ZERO_SQRT,
                                           NULL,
                                           NULL,
                                           SM_SR_SINGLE,
                                           sSolutions));
          if (sSolutions.GetSize() == 0) 
            { SER(SM_ERR); }
          
          // set new EndVert position = DropPoint and mark it SM_FIL_PROCESSED
          SmPoint3d sPV[2];
          SER(pSideCurve->Evaluate(sSolutions[0].m_vStart[0],1,TRUE,sPV));
          SmPoint3d sEndPnt = sPV[0];
          pEndVert->SetPoint(sEndPnt);
          pEndVert->SetStatus(SM_FIL_PROCESSED);

          // create a line segment[sCornerVertPnt, sEndPt] between two vertices
          // and project the line onto the extended surface through plane perp to the Line
          SmVector3d  sDir         = sCornerVertPnt - sEndPnt;
          SmVector3d  sPlaneNormal = (sPV[1]*sDir) * sDir;
          SER(sPlaneNormal.Unitize());
          SmPlane   * pPlane        = new (m_crContext) SmPlane(sEndPnt,sPlaneNormal); NER(pPlane);
          SmObjDelete sDelete(pPlane);
          SER(pPlane->MakeNurb());

          // get SideSurface/PlanePerpToLine[sCornerVertPtn sEndPnt] XSection Curve and UVTrimCurves 
          SmBSplineCurve * p3DCurve  = NULL;
          SmBSplineCurve * pUVCurve1 = NULL;
          SmBSplineCurve * pUVCurve2 = NULL;
          SER(sm_SrfSrfIntersection(m_crContext,         // in : context for new object construction
                                    pEndVert,            // in : 1st point known to be on intersection curve
                                    pAdjCornerVert,      // in : 2nd point known to be on intersection curve,
                                                         //      NULL to ignore
                                    &sDir,               // in : expected general direction of intersection curve from 1st point
                                    NULL,                // in   expected intersection end direction
                                    pSideSurface,        // in : 1st intersecting surface
                                    pPlane,              // in : 2nd intersecting surface
                                    m_dThisApproxTol3d,  // in : max allowed distance between xSect Curve and surfaces
                                    30.0*SM_PI/180.0,    // in : max allowed angle between consecutive xSect curve segment tangents
                                    p3DCurve,            // out: 3d intersection curve
                                    pUVCurve1,           // out: associated UVTrimCurve on 1st surface
                                    pUVCurve2,           // out: associated UVTrimCurve on 2nd surface
                                    FALSE,               // in : FALSE= try cheap sm_TwoPntsIntersection() before
                                                         //             resorting to expensive general surf/surf xSect solver
                                    NULL,                // in : reference point
                                    iDebugLevel ));      // in : iDebugLevel, 0 = No Debug output

          pSplitEdge->SetCurve(p3DCurve, TRUE ) ; // TRUE = delete preExisting Edge->Curve - expect none here
                                                  // side effect: delete current pSurviveEdge->UVTrimCurve
          p3DCurve->SetOwner(pSplitEdge);
          SM_ASSERT(pUVCurve1 != NULL) ; delete pUVCurve1 ; pUVCurve1 = NULL ; // not usint the SplitEdge->UVTrimCurves. gwc:Why not?
          SM_ASSERT(pUVCurve2 != NULL) ; delete pUVCurve2 ; pUVCurve2 = NULL ;

      } // end when no SplitEdge was found and bTrimSideFaces = FALSE so need to make a new pSplitEdge check

#ifdef SM_DEBUG_CODE
SmBoolean bDebugMe=FALSE;
      if (bDebugMe) 
        {
          SmEdge * pOrigEdge   = pRailOnCliff->GetOriginalEdge() ;
          SmFace * pOrigFace   = pRailOnCliff->GetOriginalFace() ;
          //SmBrep * pPseudoBrep = pRailOnCliff->GetBrep() ;
          SmBrep * pTargetBrep =   pOrigEdge ? pOrigEdge->GetBrep() 
                                 : pOrigFace ? pOrigFace->GetBrep()
                                 : NULL ;

          smgfx_Erase();
          smgfx_SetLook(1,2, 0,0,1) ; if(pTargetBrep) pTargetBrep->Draw(TRUE) ; sm_GraphicsLoop() ;
          smgfx_SetLook(5,6, 1,0,0); if (pRailStartV->IsProcessed()) pRailStartV->Draw(); sm_GraphicsLoop();
          smgfx_SetLook(5,6, 0,0,1); if (pRailEndV->IsProcessed()) pRailEndV->Draw(); sm_GraphicsLoop(); // It's OK if pRailEndV is not processed
          smgfx_SetLook(3,4, 0,1,0); pRailOnCliff->Draw(); sm_GraphicsLoop();
          smgfx_SetLook(3,4,.5,0,1); pSideEdge->Draw(); sm_GraphicsLoop();

          smgfx_SetLook(3,4, 1,0,0); pCrossSection->Draw(); sm_GraphicsLoop();
          smgfx_SetLook(1,2, 0,1,0); pSideSurface->DrawUV(2,2); sm_GraphicsLoop();
          smgfx_SetLook(7,8, 0,0,1); sCornerVertPnt.Draw(); sm_GraphicsLoop();
          smgfx_SetLook(3,4, 1,0,1); pSplitEdge->Draw(); sm_GraphicsLoop();
          sm_GraphicsLoop();
        }
#endif // SM_DEBUG_CODE
    } // end if (pAdjCornerVert)

  // when bTrimSideFaces == TRUE, done after Setting   pRailOnCliff->OriginalFace = pSideFace and 
  //                                         Adjusting VertCorner->Pos and 
  //                                         Adjusting CrossSectionEdge->Ivls 
  if (bTrimSideFaces) 
    { return SM_SUCCESS; }

  // arrive here when bTrimSideFaces == FALSE
  //    pSideSurface               = extended pSideFace->SurfaceCopy
  //    pRailOnCliff->OriginalFace = NULL ;
  //    CornerVert->pos            = CrossSectionEdge/SideSurface intersection Point
  //    CrossSectionEdge->Curve and UVTrimCurve Ivls set to run From CornerVert->Pos to CrossSection->OtherVert->Pos dropped to CrossSectionCurve Point
  //    New SplitEdge              = FilletEdge running from CornerVert->Pos to CornerVertDroppedToSideCurve->Pos

  // next: rpSideFG = New FilletGeom
  //           rpSideFG->Edgeuse       = pRailOnCliff->Edgeuse starting at pRailStartV
  //           rpSideFG->FilletSurface = pSideSurface
  //           rpSideFG->Rail0         = pRailOnCliff
  //           rpSideFG->Rail1         = New Rail: Curve          = trimmed copy of pSideEdge 
  //                                               UVTrimCurves   = NULL
  //                                               OrigEdge       = pSideEdge
  //                                               OrigFace       = pSideFace
  //                                               new Vertex1    = pRailStartV
  //                                               new Vertex2    = pRailEndV
  //                                               new Edgeuse1   = start at Vertexuse1
  //                                               new Edgeuse2   = start at Vertexuse2

  // Create new fillet geom added to FilletSolver->m_vFilletGeoms array
  rpSideFG = new SmFilletGeom(this,SM_FG_CLIFF_SIDE_PATCH);
  m_vFilletGeoms.Add(rpSideFG);

  SmFilletEdgeuse * pEU = (SmFilletEdgeuse*)pRailOnCliff->GetPrimaryEdgeuse();
  if (pEU->GetVertexuse()->GetVertex() != pRailStartV) 
    { pEU = (SmFilletEdgeuse*)pEU->GetMate(); }
  pEU->SetFilletGeom(rpSideFG);

  // Attach the surface
  rpSideFG->SetFilletSurface(pSideSurface);

  // Make some new topo objects (edge, edgeuse, vertexuse)
  rpSideFG->SetRail(0, pRailOnCliff); // Use existing rail as rail #1
  SmFilletEdge * pNewRail = new (pPseudoBrep) SmFilletEdge(pPseudoBrep);
  rpSideFG->SetRail(1, pNewRail);
  pNewRail->SetFilletEdgeType(SM_FE_RAIL);

  SmFilletEdgeuse   * pNewEU1 = new (pPseudoBrep) SmFilletEdgeuse(pPseudoBrep);
  SmFilletEdgeuse   * pNewEU2 = new (pPseudoBrep) SmFilletEdgeuse(pPseudoBrep);
  SER(pNewRail->PreInsert(pNewEU1));
  SER(pNewRail->PostInsert(pNewEU2));
  SmFilletVertexuse * pNewVU1 = new (pPseudoBrep) SmFilletVertexuse();
  SmFilletVertexuse * pNewVU2 = new (pPseudoBrep) SmFilletVertexuse();
  pNewVU1->SetProperty(pNewEU1);
  pNewVU2->SetProperty(pNewEU2);
  pNewEU1->SetFilletVertexuse(pNewVU1);
  pNewEU2->SetFilletVertexuse(pNewVU2);

  // Put pNewVU1 into rail start vertex
  SmTArray<SmEdge*> sEdges;
  SmVertex * pNewRailStartV = pRailStartV;
  pRailStartV->GetEdges(sEdges);

  // Determine if pRailStartV is connected with a 'support' edge
  // If so, will put pNewVU1 into the other end of that edge
  for (ULONG i=0; i<sEdges.GetSize(); i++) 
    {
      SmFilletEdge * pFE = (SmFilletEdge*)sEdges[i];
      if(   pFE->GetFilletEdgeType() == SM_FE_ON_EXTEND_EDGE
         || pFE->GetFilletEdgeType() == SM_FE_SPLIT_FACE) 
        {
          pNewRailStartV = pFE->GetOtherVertex(pRailStartV);
          // Attach new fillet geom to the support edge
          SmFilletEdgeuse * pFEU = (SmFilletEdgeuse*)pFE->GetPrimaryEdgeuse();
          if (pFEU->GetFilletGeom()) { pFEU = (SmFilletEdgeuse*)pFEU->GetMate(); }

          rpSideFG->AddSideFilletEdgeuse( pFEU );
          pFEU->SetFilletGeom(rpSideFG);
          break;
        }
    }
  SER(pNewRailStartV->PostInsert(pNewVU1));

  // Put pNewVU2 into rail end vertex
  SmVertex * pNewRailEndV = pRailEndV;
  pRailEndV->GetEdges(sEdges);

  // Determine if pRailEndV is connected with a 'support' edge
  // If so, will put pNewVU2 into the other end of that edge
  for (ULONG j=0; j<sEdges.GetSize(); j++) 
    {
      SmFilletEdge * pFE = (SmFilletEdge*)sEdges[j];
      if(   pFE->GetFilletEdgeType() == SM_FE_ON_EXTEND_EDGE
         || pFE->GetFilletEdgeType() == SM_FE_SPLIT_FACE) 
        {
          pNewRailEndV = pFE->GetOtherVertex(pRailEndV);
          SmFilletEdgeuse * pFEU = (SmFilletEdgeuse*)pFE->GetPrimaryEdgeuse();
          if (pFEU->GetFilletGeom()) {
            pFEU = (SmFilletEdgeuse*)pFEU->GetMate();
          }
          rpSideFG->AddSideFilletEdgeuse( pFEU );
          pFEU->SetFilletGeom(rpSideFG);
          break;
        }
    }
  SER(pNewRailEndV->PostInsert(pNewVU2));

  // Compute pRailOnEdge
  pNewRail->SetOriginalEdge(pSideEdge);
  pNewRail->SetOriginalFace(pSideFace);
  SER(pNewRail->CalcEdgeOnEdge(m_crContext));  // sets pNewRail->Curve        = trimmed copy of pNewRail->OrigEdge
                                               // sets pNewRail->UVTrimCurves = NULL
  // all done
  return SM_SUCCESS;

} // end SmFilletSolver::CreateCliffRollOverSideFG

/*******************************************************************//**
PURPOSE: Process fillets that rollover all the edges/cliffs of the face

NOTES: This routine will not handle tangent rollover. Currently,
    only closed fillet surfaces can be processed.
***********************************************************************/
SmStatus SmFilletSolver::ProcessFullRollOver
  (SmFilletGeom * pFilletGeom,  // in :
   ULONG lRollOverRailIndex)    // in :
{
int iDebugLevel = 0;
#ifdef SM_DEBUG_CODE
iDebugLevel = DebugLevel();
#endif // SM_DEBUG_CODE

  SmEdgeuse * pSolverEU = GetEdgeuse(0);

  // Determine if we need to 'add' or 'trim' side faces
  SmBoolean bTrimSideFaces = !(pSolverEU->IsConvexRadialSector(10));

  SmFilletBrep * pPseudoBrep = GetFilletExecutive()->GetPseudoBrep();
  pFilletGeom->SetFilletGeomType(SM_FG_CLIFF_ROLLOVER);

  // gwc: I think pSideSurface no longer has to be a BBSplineSurface - replace 1 line
  // rm : SmBSplineSurface * pFilletSurface  = pFilletGeom->GetFilletSurface(); NER(pFilletSurface);
  SM_FILLETSURF_TYPE * pFilletSurface  = pFilletGeom->GetFilletSurface(); NER(pFilletSurface);
  SmFilletEdge     * pRail           = pFilletGeom->GetRail(lRollOverRailIndex);
  SmCurve          * pOrigCurve      = pRail->GetCurve(); NER(pOrigCurve);
  SmFace           * pOrigFace       = pRail->GetOriginalFace(); NER(pOrigFace);
  SmFaceuse        * pFU             = pOrigFace->GetUpwardFaceuse();
  //SmSurface        * pOrigSurface = pOrigFace->GetSurface();

  // Collect the edgeuses of the outer loop
  SmLoopuse * aData[10];
  SmVertex  * aData2[20];
  SmTArray<SmLoopuse*> sLoopuses(10,aData);
  SmTArray<SmVertex*>  sVertices(20,aData2);
  SmTArray<SmEdgeuse*> sEdgeuses;
  pFU->GetLoopuses(sLoopuses);
  SmLoopuse *pLU = sLoopuses[0];
  pLU->GetEdgeuses(sEdgeuses);

  // Determine which edge the start-of-the-rail is 'faced' to
  // by finding the minimum distance from the rail start to each edge
  SmFilletVertex * pStartV       = (SmFilletVertex*)pRail->GetVertex();
  SmPoint3d        sRailStartPnt = pStartV->GetPoint();
  SmVector3d       sDir;
  SmPoint3d        sPnts[2];

  //ULONG lStartEdgeIndex = 0;
  double      dMinDist    = SM_BIG_DOUBLE;
  ULONG       lTotalEdges = sEdgeuses.GetSize();
  SmEdgeuse * pStartEU    = NULL;

  //
  for (ULONG i=0; i<lTotalEdges; i++) 
    {
      SmEdgeuse * pEU = sEdgeuses[i];
      SmEdge    * pSideEdge = pEU->GetEdge();
      SmBoolean bOppositeNormals;
      if (sm_TestManifoldTangentRollover(pSideEdge, bOppositeNormals, m_dTangencyTolerance, iDebugLevel )) 
        {
          MSG(_T("Unimplemneted Tangent Rollover Case"));
          SER(SM_ERR);
        }

      SmCurve  * pSideCurve = pSideEdge->GetCurve();
      SmExtent1d sIvl       = pSideEdge->GetInterval();
      SmSolution sData[16];
      SmSolutionArray sSolutions(16,sData);

      SER(pSideCurve->GlobalPointSolve(sIvl,
                                       SM_SO_MINIMIZE,
                                       sRailStartPnt,
                                       SM_EFF_ZERO_SQRT,
                                       NULL,
                                       NULL,
                                       SM_SR_SINGLE,
                                       sSolutions));
      if (sSolutions.GetSize() > 0) 
        {
          double dDist = sSolutions[0].m_vStart.m_dSolutionValue;
          if (dDist < dMinDist) 
            {
              dMinDist =  dDist;
              pStartEU = pEU;
              SER(pSideCurve->Evaluate(sSolutions[0].m_vStart[0],1,TRUE,sPnts));
              sDir = sPnts[1];
            }
        }
    }

  // Is the rail running the same direction as the loop
  SmExtent1d sOrigIvl = pRail->GetInterval();
  SER(pOrigCurve->Evaluate(sOrigIvl.GetMin(),1,TRUE,sPnts));
  //SmBoolean  bTracingForward = TRUE;

  SmOrientType eOrient = sPnts[1].Dot(sDir) > 0.0? SM_OT_SAME : SM_OT_OPPOSITE;

  // Get the edgeuse that goes along the right direction (i.e. rail's direction)
  if (eOrient != pStartEU->GetOrientation()) 
    {
      pStartEU = pStartEU->GetMate(); NER(pStartEU);
    }

  pRail->SetFilletEdgeType(SM_FE_CLIFF_RAIL);
  SmTArray<SmFilletVertex*> sAllNewVerts;
  sAllNewVerts.SetDataSize(lTotalEdges);

  // Traverse the loop
  SmEdgeuse * pCurrEU = pStartEU;
  for (ULONG j=0; j<=lTotalEdges; j++) 
    {
      SmEdgeuse * pNextEU   = pCurrEU->GetCCWEdgeuse();
      SmEdge    * pSideEdge = pCurrEU->GetEdge();
      SmEdgeuse * pRadialEU = pCurrEU->GetRadial();
      SmFace    * pSideFace = pRadialEU->GetFace();
      NER(pSideFace);
      SmVertex  * pCurrV    = pNextEU->GetVertexuse()->GetVertex();

#ifdef SM_DEBUG_CODE
SmBoolean bDebugMe = FALSE;
      if (bDebugMe) 
        {
          smgfx_SetColor(1,0,0); pCurrEU->Draw(); sm_GraphicsLoop();
          smgfx_SetColor(0,0,1); pRadialEU->GetCWEdgeuse()->Draw(); sm_GraphicsLoop();
          smgfx_SetColor(0,0,1); pSideFace->DrawUV(4,4); sm_GraphicsLoop();
          smgfx_SetLook(8,1, 0,0,1); pCurrV->Draw(); sm_GraphicsLoop();
          sm_GraphicsLoop();
      }
#endif // SM_DEBUG_CODE

      SmFilletEdge * pNewRail = NULL;
      if (j != lTotalEdges) // Not for the last segment
        { 
          // Split the rail now
          SmFilletVertex * pNewV    = NULL;
          SmFilletEdge   * pNewEdge = NULL;

          SER(sm_SplitRail(pFilletGeom,pRail,pNewV,pNewRail));

          pNewV->SetFilletVertexType(SM_FV_CLIFF_RAIL_X_EDGE);
          pNewV->SetPointClass(SM_PC_EDGEUSE,pRadialEU->GetCWEdgeuse());
          sAllNewVerts.Add(pNewV);

          if (lRollOverRailIndex == 0) { pFilletGeom->m_vOtherRails1.Add(pNewRail); }
          else                         { pFilletGeom->m_vOtherRails2.Add(pNewRail); }
          pNewRail->SetFilletEdgeType(SM_FE_CLIFF_RAIL);

          if (!bTrimSideFaces) 
            {
              // Create an edge connecting pNewV & pCurrV as the 'support' edge
              SmFilletVertex * pEndVert = new (pPseudoBrep) SmFilletVertex();
              pEndVert->SetPoint(pCurrV->GetPoint());
              pEndVert->SetStatus(SM_FIL_PROCESSED);
              SER(SmFilletEdge::MakeFilletEdge(pPseudoBrep,  // in : target Brep to receive new topology objects
                                               pEndVert,     // in : start of new FilletEdge
                                               pNewV,        // in : end   of new FilletEdge
                                               pNewEdge));   // out: newly allocated FilletEdge - no Curve or UVTrimCurve data
                                                             // in : NewFilletEdge's filletCorner, gets stored in rpNewFilletEdge->m_cpCorner
                                                             //      NULL to ignore, default:[NULL]
              pNewEdge->SetFilletEdgeType(SM_FE_ON_EXTEND_EDGE);
            }
        }

      // Now this rail become a rail along the cliff

      // Create a side fillet geom which will be coincided
      // with the extended side-face
      SmFilletGeom * pNewSideFG = NULL;
      pRail->SetFilletEdgeType(SM_FE_CLIFF_RAIL);

      double dExtensionDist = smos_Fabs(GetSurface(0)->GetOffsetDistance());
      SER(CreateCliffRollOverSideFG(pRail,
                                    pSideEdge,
                                    pSideFace,
                                    bTrimSideFaces,
                                    dExtensionDist,
                                    pNewSideFG));

      if (pNewRail) 
        { pRail = pNewRail; }
      pCurrEU = pNextEU;
    }

  // Calculate the newly created vertices
  for (ULONG k=0; k<lTotalEdges; k++) 
    {
      SmFilletVertex * pNewV = sAllNewVerts[k];
      if (bTrimSideFaces) 
        {
          // Intersect the curve with the fillet for the pNewV's geomtry
          SER(pNewV->CalcCliffRailIntEdge(pFilletSurface,m_dThisApproxTol3d));
        }
      else 
        {
          // Find the 'support' edge and compute its geometry
          SmTArray<SmEdge*> sEdges;
          pNewV->GetEdges(sEdges);
          SM_ASSERT(sEdges.GetSize() == 3);
          SmFilletEdge * pEdgeOnEdge = NULL;
          for (ULONG ii=0; ii<3; ii++) {
              SmFilletEdge * pFE = (SmFilletEdge*)sEdges[ii];
              if (pFE->GetFilletEdgeType() == SM_FE_ON_EXTEND_EDGE) 
                {
                  pEdgeOnEdge = pFE;
                  break;
                }
            }
          NER(pEdgeOnEdge);

          // Computes the 3D Geometry (Curve and endVertexPoint) for pSupportEdge - existing Curve and UVTrimCurves are deleted
          SER(sm_ComputeSupportEdge(m_crContext,
                                    pEdgeOnEdge,
                                    pFilletSurface,
                                    m_dThisApproxTol3d,
                                    m_dThisAngTolRad, 
                                    iDebugLevel ));
        }
    }

  // Finally, trim all rails of the fillet
  SER(pFilletGeom->TrimCliffRails());

  // all done
  return SM_SUCCESS;

} // end SmFilletSolver::ProcessFullRollOver

/*******************************************************************//**
PURPOSE: Process any rollovers and create new surfaces as needed.
     Rollover occurs when a fillet rail intersects some other edge
     of the target Brep along the way.

     There are three basic types of roll over, depending on the
     relative convexity of the other edge intersected by the rail.
     If the other edge has opposite convexity to the filleted edge,
     then the rail rolls over a cliff.  If the convexities are the
     same, then the fillet surface crashes into the adjoining face.
     If the other edge is a tangent edge, then the fillet will switch
     to rolling along the adjoining face.  See the Usage Notes for
     SmFilletSolver::TestRollOver() for more explanation.

     In the case of the cliff roll over there are two possible ways
     to handle the situation. 1) Roll ball along the edge.
     2) Trim fillet surface by extending the adjacent face.
     We currently use only the second method.

     In the case of a tangent roll over the best thing to do is to
     just roll onto that adjacent surface.  This requires creating
     a new FilletGeom, for filleting on the other face.

NOTES:
     Right now we will primarily worry about handling the situation
     where there is only one edge between the gaps of the roll over.
     Later we will extend to roll to multiple edges along a path.

     Note that the current logic only handles rolling over adjacent
     tangent surface.

***********************************************************************/
SmStatus SmFilletSolver::ProcessRollOver
  (SmFilletGeom * pFilletGeom,
   SmTArray<SmFilletCorner*> & rAdjustedCorners)
{
  // If we have a blending patch, skip processing
  if (pFilletGeom->GetFilletGeomType() == SM_FG_BLENDS)
    { return SM_SUCCESS; }

int iDebugLevel = 0;
#ifdef SM_DEBUG_CODE
iDebugLevel = DebugLevel();
#endif // SM_DEBUG_CODE

  // collect some info about the rails: their pointers, their
  // curve classifications, and the number of intervals in their CC's.
  SmFilletEdge          * pRail[2];
  SmCurveClassification * pCC[2];
  ULONG                   lNumIvls[2];
  for (ULONG k=0; k<2; k++)
    {
      // get rail->CurveClassification->IntervalCount
      pRail[k]    = pFilletGeom->GetRail(k);
      pCC[k]      = pRail[k]->GetCurveClass(); NER(pCC[k]);
      lNumIvls[k] = pCC[k]->GetSize();

#ifdef SM_DEBUG_CODE
      if ( iDebugLevel > 0 ) 
        {
          smgfx_SetLook(1,1, 1,0,0);
          pRail[k]->Draw();
          sm_GraphicsLoop();
        }
#endif // SM_DEBUG_CODE
    }

  // If both rails have only a single classification interval,
  // we can take care of everything here.
  if ( lNumIvls[0] == 1 && lNumIvls[1] == 1 )
    {
      SmCurveInterval & rIvl1 = (*pCC[0])[0];
      SmCurveInterval & rIvl2 = (*pCC[1])[0];

      // If both rails classify to something in their original Faces,
      // then there is not rollover problem
      if (   rIvl1.m_vMid.GetPointClass() != SM_PC_UNKNOWN
          && rIvl2.m_vMid.GetPointClass() != SM_PC_UNKNOWN )
        {
          return SM_SUCCESS;
        }

      // else Try FullRollover recovery for certain corner types only
      SmFilletVertex * pV            = (SmFilletVertex*)pRail[0]->GetVertex();
      SmFilletCorner * pCorner       = pV->GetFilletCorner(); NER(pCorner);
      SmFilletCornerType eCornerType = pCorner->GetCornerType();

      // These are the only types we can do:
      if (   eCornerType != SM_FCR_N_x_1_CLOSED
          && eCornerType != SM_FCR_1_x_1)
        {
          SER_MSG(SM_ERR, _T("ProcessRollOver(): Unimplemented Roll-Over Case"));
        }

      // gwc:change - changed things so ProcessFullRollOver
      //              is called on both rails when both
      //              are unclassified not just the 2nd rail.

      // Process rail1 rollover problems
      if ( rIvl1.m_vMid.GetPointClass() == SM_PC_UNKNOWN )
        {
          SER( ProcessFullRollOver( pFilletGeom, 0 ));
        }

      // Process rail2 rollover problems
      if ( rIvl2.m_vMid.GetPointClass() == SM_PC_UNKNOWN )
        {
          SER( ProcessFullRollOver( pFilletGeom, 1 ));
        }

      return SM_SUCCESS;

    } // end only 1 interval in each rail curveClassification check

  // Ok, at least one rail has more than one interval
  // (and so has hit an edge of the target Brep).

  // Determine if this is a convex fillet (e.g. boss-on-box).

  //  Aside, about concave/convex.  First, the fillet will always go into
  //  the concave sector -- it has to.  This check is for the convexity
  //  with respect to the Edgeuses that were selected to define the fillet.
  //  That selection is apparently done by picking the EUs that point "into"
  //  a solid body: on the solid side, not the infinite region.
  //  However, that concept is meaningless for filleting: the fillet will
  //  go into the concave sector, it doesn't matter whether it's inside or
  //  outside, solid or void.
  //  What does matter is relative convexity: whether two edges have the same
  //  or opposite convexities, looking from either side.  That's what we do here.

  SmEdgeuse * pSolverEU = GetEdgeuse(0);
  SmBoolean bFilletIsConvex = pSolverEU->IsConvexRadialSector(10);

  // Determine if this is an analytic fillet
  // gwc: I think pSideSurface no longer has to be a BBSplineSurface - replace 1 line with 2
  // rm : SmBSplineSurface * pFilletSurface = pFilletGeom->GetFilletSurface(); NER(pFilletSurface);
  SM_FILLETSURF_TYPE * pFilletSurface       = pFilletGeom->GetFilletSurface(); NER(pFilletSurface);
  SmBSplineSurface   * pBSplineFilletSurface = SM_CAST_PTR(SmBSplineSurface, pFilletSurface) ;

  SmSurfOfRevolution * pRev        = NULL;
  // rm : SmBoolean            bIsAnalytic = SmSurfOfRevolution::IsNurbSurfaceSurfOfRevolution(m_crContext, pFilletSurface, pRev);
  SmBoolean            bIsAnalytic =   pBSplineFilletSurface 
                                     ? SmSurfOfRevolution::IsNurbSurfaceSurfOfRevolution(m_crContext, pBSplineFilletSurface, pRev) 
                                     : FALSE ;
  if ( pRev ) 
    {
      delete pRev; pRev = NULL; // (Don't need this, just checking analytic.)
    }

  // Make copies of original rail curves and fillet surface.
  SmCurve   * pOrigRailCurve[2];
  SmCurve   * pOrigRailUVCurve[2];
  SmObjDelete sDelete3DCurve[2];
  SmObjDelete sDeleteUVCurve[2];

  for (ULONG m=0; m<2; m++) 
    {
      SER( pRail[m]->GetCurve()->Copy( m_crContext, pOrigRailCurve[m] ));
      sDelete3DCurve[m].SetObj( pOrigRailCurve[m] );
      SmCurve * pUVCurve = pRail[m]->GetPrimaryEdgeuse()->GetUVTrimCurvePointer();
      NER( pUVCurve );
      SER( pUVCurve->Copy( m_crContext, pOrigRailUVCurve[m] ));
      sDeleteUVCurve[m].SetObj( pOrigRailUVCurve[m] );
    }

  SmSurface * pOrigFillet = NULL;
  SER( pFilletSurface->Copy( m_crContext, pOrigFillet ));
  SmObjDelete sDeleteSurface( pOrigFillet );

  SmBoolean bSplitRailLastTime[2] = { FALSE, FALSE };
  SmBoolean bProcessingRollover   = FALSE;
  ULONG     lCurrIndx[2]          = { 0, 0 };
  //SmFace * pAdjFace[2]  = { NULL, NULL };
  //SmEdge * pSideEdge[2] = { NULL, NULL };
  SmEdgeuse * pRolloverCornerEU   = NULL;
  SmVertex  * pRolloverCornerV    = NULL;

  // Temporary array to store all support edges for cliffs
  SmTArray<SmFilletEdge*> sCliffSupportEdges;

  SmFilletGeom * pCurrFilletGeom = pFilletGeom;

  // Walk along both rails, segment by segment, until the end.
  while ( lCurrIndx[0] < lNumIvls[0]  &&  lCurrIndx[1] < lNumIvls[1] )
    {
      SmBoolean bProcessingFirstGeom = lCurrIndx[0] == 0              && lCurrIndx[1] == 0;
      SmBoolean bProcessingLastGeom  = lCurrIndx[0] == lNumIvls[0]-1  && lCurrIndx[1] == lNumIvls[1]-1 ;

      // Walk through each interval and determine where to split
      // this SmFilletGeom. Will adjust solver to handle rollover
      // to adjacent faces (even on both sides)

      // Backup original offsets data
      SmFilletSolverOffsetsData sOrigOffsets;
      GetOffsetsData( sOrigOffsets );

      SmFilletGeomType eFGType = SM_FG_DEFAULT;

      // Pointers to the side edge that rail[i] hit, and the adjoining
      // face on the other side of that edge;
      // null if rail[i] doesn't intersect any edges of the target Brep.
      SmEdge * pSideEdge[2] = { NULL, NULL };
      SmFace * pAdjFace[2]  = { NULL, NULL };

      // Whether pAdjFace's normal is opposite to that of the base face.
      SmBoolean bOppositeNormals[2] = { FALSE, FALSE };

      // whether pSideEdge[i] is a tangent edge.
      SmBoolean bTangentialSideEdge[2] = { FALSE, FALSE };

      // If the side edge has the same convexity as the filleted edge,
      // then the fillet surface will crash into the adjoining face,
      // and we'll trim it.
      // If the side edge has opposite convexity to the filleted edge,
      // the fillet surface will be hanging over a cliff,
      // and we'll extend the side face (or add a new little face)
      // to intersect the fillet.
      // Otherwise, the side edge is tangent, and we'll split the fillet into two.
      SmBoolean bTrimSideFaces[2] = { FALSE, FALSE };

      // Loop over both rails, and process the one that gets split.
      ULONG lWhichRail;
      for(lWhichRail=0; lWhichRail<2; lWhichRail++ )
        {
          SmEdge * pPrevSideEdge = pSideEdge[lWhichRail];
          pAdjFace [lWhichRail] = NULL;
          pSideEdge[lWhichRail] = NULL;

          // If this rail isn't split, skip it
          if ( lNumIvls[lWhichRail] == 1 )
            { continue; }

          // This next call does a lot more than find the side edge(s)
          // intersected by the rail(s): see comments for output args.

          SER( sm_FindRolloverSideEdge
                   (pCC[lWhichRail],
                    lCurrIndx[lWhichRail], // in/out: might get incremented.
                    pRail[lWhichRail],
                    pOrigRailCurve[lWhichRail], pOrigRailUVCurve[lWhichRail],
                    pOrigFillet,  bFilletIsConvex,
                    pSideEdge           [ lWhichRail ], // out: edge we hit
                    bTangentialSideEdge [ lWhichRail ], // out: is it a tangent edge?
                    bTrimSideFaces      [ lWhichRail ], // out: will fillet srf hit adjacent face?
                    bOppositeNormals    [ lWhichRail ], // out: of the two faces of side edge
                    pRolloverCornerEU,  // out
                    pRolloverCornerV,   // out: vtx at rail/edge intersection
                    m_dTangencyTolerance, // in
                    iDebugLevel));

          if ( pSideEdge[lWhichRail] == NULL )
            {
              // This rail didn't hit any edges in the target Brep.
              if ( lCurrIndx[lWhichRail] == lNumIvls[lWhichRail]-1 )
                {
                  // This is the final segment on this rail

                  // Note, we know that lNumIvls[lWhichRail] > 1,
                  // otherwise we wouldn't be here.

                  pSideEdge[lWhichRail] = pPrevSideEdge;
                  bProcessingLastGeom = TRUE;
                }
              continue;  // go on to next rail
            }

          // If this interval hangs outside of its face
          // (indicated by its mid being classified unknown),
          // determine its rollover type.
          SmCurveInterval & rIvl = (*pCC[lWhichRail])[lCurrIndx[lWhichRail]];

          if ( rIvl.m_vMid.GetPointClass() == SM_PC_UNKNOWN )
            {
              // Find the adjacent face -- the one on the other side of
              // the side edge that we hit.
              SmFace *pOrigFace = pRail[lWhichRail]->GetOriginalFace();
              NER(pOrigFace);
              SmTArray<SmFace*> sFaces;
              pSideEdge[lWhichRail]->GetFaces( sFaces );
              for ( ULONG jj=0; jj<sFaces.GetSize(); jj++ )
                {
                  if ( sFaces[jj] != pOrigFace )
                    {
                      pAdjFace[lWhichRail] = sFaces[jj];
                      break;
                    }
                }
              NER( pAdjFace[lWhichRail] ); // Adjacent face

              SmSurface * pAdjSurface = pAdjFace[lWhichRail]->GetSurface();

              // If tangent-rollover, alter the solver to solve for new fillet
              if ( !bTangentialSideEdge[lWhichRail] || bTrimSideFaces[lWhichRail] )
                {
                  // We'll not do tangent rollover when
                  // same convexity found on both sides
                  eFGType = SM_FG_CLIFF_ROLLOVER;
                }
              else
                {
                  eFGType = SM_FG_TANGENT_ROLLOVER;

                  // Create new offset surface
                  double dOffsetDist = GetSurface(lWhichRail)->GetOffsetDistance();
                  if ( bOppositeNormals[lWhichRail] )
                    {
                      dOffsetDist *= -1.0;
                      m_dOrientations[lWhichRail] *= -1.0;
                    }
                  SmOffsetSurface * pNewOffSrf = new (m_crContext) SmOffsetSurface(
                          dOffsetDist, *pAdjSurface );
                  pNewOffSrf->SetAllowExtension(TRUE);
                  SetSurface( lWhichRail, pNewOffSrf );

                  // Put a new offset surf into pCurrFilletGeom as well.
                  pNewOffSrf = new (m_crContext) SmOffsetSurface(
                          dOffsetDist, *pAdjSurface );
                  pNewOffSrf->SetAllowExtension(TRUE);
                  SER( pCurrFilletGeom->SetOffsetSurface( lWhichRail, pNewOffSrf ) );
                }
              pCurrFilletGeom->SetFilletGeomType( eFGType );
            } // end if mid UNKNOWN (hanging ouside face)
        } // end Loop over both rails, processing the one that gets split.

      // Check which rail's interval ends first (or both).
      // Set guessParam to the first hit,
      // and increment lCurrIndx of the rail that hit.
      double dParam1 = (*pCC[0])[lCurrIndx[0]].m_vInterval.GetMax();
      double dParam2 = (*pCC[1])[lCurrIndx[1]].m_vInterval.GetMax();
      double dParamTol = 1.0e-3;
      double dGuessParam = 0.0;
      if ( pRolloverCornerEU == NULL )
        {
          if ( smos_Fabs(dParam1-dParam2) < dParamTol )
            {
              // done with both intervals.
              dGuessParam = dParam1;
              lCurrIndx[0]++; lCurrIndx[1]++;
            }
          else if ( dParam1 < dParam2 )
            {
              // interval on rail[0] ended first.
              dGuessParam = dParam1;
              lCurrIndx[0]++;
            }
          else
            {
              dGuessParam = dParam2;
              lCurrIndx[1]++;
            }
        }

      // If we have no side edges, we might want to quit.
      if ( pSideEdge[0] == NULL  &&  pSideEdge[1] == NULL )
        {
          if ( !bProcessingLastGeom )
            {
              continue; // Proceed to the next segments.
            }
          else if ( !bProcessingRollover )
            {
              return SM_SUCCESS;
            }
        }

      bProcessingRollover = TRUE;

      // Now, we'll split the fillet geom or just split the rail
      // depending upon the rollover type.
      // If we do split the geom, then we'll set pNewFilletGeom;
      // if we split rails, then pNewRail1/pNewRail2 will be set.

      SmFilletGeom * pNewFilletGeom = NULL;
      SmFilletEdge * pNewRail1 = NULL;
      SmFilletEdge * pNewRail2 = NULL;
      SmFilletVertex * pNewV[2] = { NULL, NULL };

      if ( bProcessingLastGeom )
        {
          // If the LAST geom is tangent-rollover, set original face of rails,
          // but don't split anything.
          for ( ULONG ii=0; ii<2; ii++ )
            {
              if ( GetSurface(ii) != sOrigOffsets.m_pSurfaces[ii] )
                {
                  pRail[ii]->SetOriginalFace( pAdjFace[ii] );
                }
            }
        }
      else
        {
          // For cliff, split only the rail curve.
          if (  (eFGType == SM_FG_CLIFF_ROLLOVER || eFGType == SM_FG_DEFAULT)
                   && ( !bTangentialSideEdge[0] && !bTangentialSideEdge[1] ))
            {
              if ( pSideEdge[0] != NULL )
                {
                  // Split Rail1
                  SER( sm_SplitRail(pCurrFilletGeom, pRail[0], pNewV[0], pNewRail1) );
                  pCurrFilletGeom->m_vOtherRails1.Add( pNewRail1 );
                }
              if ( pSideEdge[1] != NULL )
                {
                  // Split Rail2
                  SER( sm_SplitRail(
                      pCurrFilletGeom, pRail[1], pNewV[1], pNewRail2 ));
                  pCurrFilletGeom->m_vOtherRails2.Add( pNewRail2 );
                }
            }
          else  // Tangent rollover: split the entire geom.
            {
              SER( pCurrFilletGeom->TopologySplit( pNewFilletGeom ));
              m_vFilletGeoms.Add( pNewFilletGeom );

              // Grab the two rail vertices at the split point into pNewV,
              // and set them as each others' mates.
              // Note, they're the ones at the top end of the low half.
              for (ULONG kk=0; kk<2; kk++)
                {
                  SmVertex * pV = pRail[kk]->GetVertex();
                  pNewV[kk] = (SmFilletVertex*)pRail[kk]->GetOtherVertex(pV);
                }
              pNewV[0]->SetMate(0,pNewV[1]);
              pNewV[1]->SetMate(0,pNewV[0]);

              // One of the FilletGeom's base surface needs to be changed
              // to that of the adjacent face that we rolled onto.
            }

          // For both rails, if appropriate, update the base face of the rail
          // and the classification of pNewV.  Also do the rail/edge
          // intersection if needed.
          for (ULONG ii=0; ii<2; ii++)
            {
              // If tangent-rollover, set original face of rails,
              // and offset surface in one of the FilletGeoms.
              // Note, in the tangent-rollover case, we have changed
              // our base surface to the adjacent face that we
              // rolled onto.  So we check whether that has changed.
              if ( GetSurface(ii) != sOrigOffsets.m_pSurfaces[ii] )
                {
                  pRail[ii]->SetOriginalFace( pAdjFace[ii] );

                  // Offset surface in FilletGeoms:
                  // If rolling from adjacent rollover face back onto
                  // original face, then pCurrFilletGeom's type is TANGENT_ROLLOVER.
                  // In that case, pCurrFilletGeom's offset surface has already been
                  // set accordingly, above, when we set eFGType to TANGENT_ROLLOVER.
                  // Then the new FG's offset srf needs to be set back to that of
                  // the original offset surface.
                  //
                  // On the other hand, if rolling from the original face onto
                  // the new adjacent face, we need to set pNewFilletGeom's offset
                  // surface to that of the adjacent face.  However, in that case,
                  // we haven't found the adjacent face.  We'll take care of that
                  // on the next iteration of the interval-segments loop.

                  // SmOffsetSurface::Copy() returns an SmSurface*, not SmOffsetSurface*
                  SmSurface *pSrfCopy = NULL;
                  SER( sOrigOffsets.m_pSurfaces[ii]->Copy( m_crContext, pSrfCopy ) );
                  SmOffsetSurface *pOffsetCopy = SM_CAST_PTR(SmOffsetSurface, pSrfCopy);
                  if ( pOffsetCopy != NULL )
                    {
                      SER( pNewFilletGeom->SetOffsetSurface( ii, pOffsetCopy ) );
                    }
                }

              // Compute where the rail should end
              if ( pSideEdge[ii] == NULL )
                { continue; }

//cbi: I don't know what pRolloverCornerEU is (or means).

              if ( bTangentialSideEdge[ii] || pRolloverCornerEU == NULL )
                {
                  NER( pNewV[ii] );
                  pNewV[ii]->SetFilletVertexType( SM_FV_RAIL_X_EDGEUSE );
                  if ( pNewV[1-ii] )
                    {
                      pNewV[1-ii]->SetFilletVertexType( SM_FV_MATE );
                    }
                  // Calculate SM_CD_RAIL_X_EDGE
                  SmPoint3d * pGuess3DPnt = NULL;
                  SmPoint3d s3DPnt;
                  if ( pRolloverCornerEU == NULL )
                    {
                      SER( pOrigRailCurve[ii]->EvaluatePoint(
                          dGuessParam, s3DPnt ));
                      pGuess3DPnt = &s3DPnt;
                    }
                  pNewV[ii]->SetPointClass( SM_PC_EDGE, pSideEdge[ii] );

                  // This puts the result into the TSectPnt of the FilletVertexuse of pNewV[ii].

                  // Note: this sets the uv value of the surface into the FilletVertex's
                  // SmPointClass object.  If it's on a seam of a closed surface, that
                  // uv value can be on the wrong side of the seam.  Unfortunately, we
                  // can't tell at this point which side of the seam is the correct side.
                  // We have to check for that in SmFilletGeom::ReCalcFilletGeom().  [B398]

                  SER( sm_CalcRailXEdge( pNewV[ii], 
                                         pCurrFilletGeom,
                                         ii, 
                                         pRail[ii], 
                                         pSideEdge[ii], 
                                         pGuess3DPnt, 
                                         iDebugLevel ));
                }
              else 
                {
                  // New vertex is at the intersection of fillet & pRolloverCornerEU
                  pNewV[ii]->SetFilletVertexType( SM_FV_CLIFF_RAIL_X_EDGE );
                  pNewV[ii]->SetPointClass( SM_PC_EDGEUSE, pRolloverCornerEU );
                }
            } // end loop over both rails
        }  // end if-else (not) processing last geom

      // If we split a rail on the previous segment,
      //   and the current segment of that rail is SM_FV_RAIL_X_EDGEUSE,
      //   and it's a tangent rollover,
      // then recalculate the rail/edge intersection.

      if ( !bProcessingFirstGeom )
        {
          for (ULONG ii=0; ii<2; ii++)
            {
              if ( !bSplitRailLastTime[ii] )
                { continue; }

              SmFilletVertex * pStartFV =
                  (SmFilletVertex*)pRail[ii]->GetVertex();

              if ( pStartFV->GetFilletVertexType() != SM_FV_RAIL_X_EDGEUSE )
                { continue; }

              // Recompute tsect point of this vertex of this fillet geom
              SmEdge * pE = (SmEdge*)pStartFV->GetPointClassObject();NER(pE);
              SmBoolean bFlag;
              SmPoint3d sGuessPnt = pStartFV->GetPoint();
              if ( sm_TestManifoldTangentRollover( pE, bFlag, m_dTangencyTolerance, iDebugLevel ))
                {
                  // This puts the result into the TSectPnt of the FilletVertexuse of pNewV[ii].

                  // Also, see note just above about seams.  [B398]

                  SER( sm_CalcRailXEdge( pStartFV, 
                                         pCurrFilletGeom,
                                         ii, pRail[ii], 
                                         pE, 
                                         &sGuessPnt, 
                                         iDebugLevel ));
                }
            }
        }

      if ( eFGType == SM_FG_CLIFF_ROLLOVER )
        {
          // Determine if we need to 'trim' or 'add' side faces.
          // In cliff rollover cases, one of the rails shall be the
          // intersections between fillet & (extended-)side face
          // In addition, the extension of the side face will fill in
          // the 'hole' as a side fillet geom

          // Create side fillet geom
          for ( ULONG ii=0; ii<2; ii++ )
            {
              if ( pSideEdge[ii] == NULL )
                { continue; } // this rail didn't hit anything.

              if ( pRolloverCornerEU && !bTrimSideFaces[ii] )
                {
                  SmFilletBrep * pPseudoBrep = GetFilletExecutive()->GetPseudoBrep();

                  // Rolling over a corner, add an edge which connects
                  // pNewV & pRolloverCornerV as the 'support' edge
                  SmFilletVertex * pEndVert = new (pPseudoBrep) SmFilletVertex();
                  pEndVert->SetPoint( pRolloverCornerV->GetPoint() );
                  pEndVert->SetFilletVertexType( SM_FV_ON_VERTEX );
                  pEndVert->SetStatus( SM_FIL_PROCESSED );

                  SmFilletEdge * pNewEdge = NULL;
                  SER( SmFilletEdge::MakeFilletEdge( pPseudoBrep, // in : target Brep to receive new topology objects
                                                     pEndVert,    // in : start of new FilletEdge
                                                     pNewV[ii],   // in : end   of new FilletEdge
                                                     pNewEdge )); // out: newly allocated FilletEdge - no Curve or UVTrimCurve data
                                                                  // in : NewFilletEdge's filletCorner, gets stored in rpNewFilletEdge->m_cpCorner
                                                                  //      NULL to ignore, default:[NULL]
                  pNewEdge->SetFilletEdgeType( SM_FE_ON_EXTEND_EDGE );
                  sCliffSupportEdges.Add( pNewEdge );
                }

              // Now this rail becomes a rail along the cliff
              SmFilletGeom * pNewSideFG = NULL;
              pRail[ii]->SetFilletEdgeType( SM_FE_CLIFF_RAIL );
              double dExtensionDist = smos_Fabs( sOrigOffsets.m_pSurfaces[ii]->GetOffsetDistance() );

              SER( CreateCliffRollOverSideFG( pRail[ii], pSideEdge[ii],
                                              pAdjFace[ii], bTrimSideFaces[ii],
                                              dExtensionDist, pNewSideFG ));
              break;
            }
        } // end if Cliff Rollover

      // 'else' removed: was missing 1st & last corners. [bd; 060216a]
      if ( bProcessingFirstGeom || bProcessingLastGeom )
        {
          // Non-cliff rollover cases. Will find
          // adjacent corner which need to be recomputed
          for ( ULONG lRailIdx=0; lRailIdx<2; lRailIdx++ )
            {
              SmFilletEdge *pFERail = pCurrFilletGeom->GetRail( lRailIdx );
              SmFilletVertex * pV = (SmFilletVertex*)pFERail->GetVertex();
              if ( bProcessingLastGeom )
                {
                  pV = (SmFilletVertex*)pFERail->GetOtherVertex(pV);
                }
              if ( pV->GetFilletVertexType() != SM_FV_RAIL_X_EDGEUSE )
                { continue; }

              SmFilletCorner * pAdjCorner = pV->GetFilletCorner();
              if ( pAdjCorner == NULL ) 
                { continue; }

              // Found adjacent fillet corner that needs adjustments.
              // Will recompute vertices of this corner here.
              // The edge will be recomputed later on by SmFilletExecutive.
              rAdjustedCorners.AddUnique( pAdjCorner );

#ifdef SM_DEBUG_CODE
SmBoolean bDebugMe1 = FALSE;
              if (bDebugMe1) 
                {
                  SmCurve * pCurve = pFERail->GetCurve();
                  SmExtent1d sIvl = pCurve->GetNaturalInterval();

                  smgfx_SetLook(1,6, 1,0,0); pV->Draw(); sm_GraphicsLoop();
                  smgfx_SetLook(3,6, 1,0,1); pCurve->DrawAt( sIvl.Evaluate(0.1), 1 ); sm_GraphicsLoop();
                  sm_GraphicsLoop();
                }
#endif // SM_DEBUG_CODE

              if ( eFGType == SM_FG_TANGENT_ROLLOVER )
                {
                  // Re-compute corner vertex
                  pV->SetStatus( SM_FIL_UNPROCESSED );
                  SmFilletVertex * pMate = pV->GetMate(0);
                  if ( pMate != NULL )
                    {
                      pMate->SetStatus( SM_FIL_UNPROCESSED );
                      if ( lRailIdx == 0 )
                          pMate->SetPointClass( SM_PC_FACE, pAdjFace[1] );
                      else
                          pMate->SetPointClass( SM_PC_FACE, pAdjFace[0] );
                    }
                  SER( pV->CalcRailIntEdgeuse( pCurrFilletGeom ));
                  if ( pV->GetStatus() == SM_FV_SOLVER_NOT_CONVERGE )
                    {
                      SER(SM_ERR);
                    }
                }
              break;
            } // end loop over both rails looking for corner vtx to recompute.
        } // end if processing first or last geom

      // If we have created a new fillet geom, then recalculate the
      // fillet surface to include new geometry we've created in all this.
      // We create new fillet geom for tangential side edges,
      // and also do it at the end no matter what.

      if (    bTangentialSideEdge[0]
           || bTangentialSideEdge[1]
           || bProcessingLastGeom)
        {
          // Blow away pCurrGeom's fillet surface
          // gwc: I think pSideSurface no longer has to be a BBSplineSurface - replace 1 line
          // rm : SmBSplineSurface * pSurf = pCurrFilletGeom->GetFilletSurface();
          SM_FILLETSURF_TYPE * pSurf = pCurrFilletGeom->GetFilletSurface();

          if ( pSurf )
            {
              SmObjDelete sDelete(pSurf);
              pCurrFilletGeom->SetFilletSurface( NULL );
            }

          // Compute a 3d direction for ReCalcFilletGeom()
          SmVector3d sPV[2];
          SmVector3d * p3DDir = NULL;

          SmPoint3d sPnt = pCurrFilletGeom->GetRail(0)->GetVertex()->GetPoint();
          SmSolutionArray sSolutions1;
          SER( pOrigRailCurve[0]->GlobalPointSolve(pOrigRailCurve[0]->GetNaturalInterval(),
                                                   SM_SO_MINIMIZE, 
                                                   sPnt, 
                                                   SM_EFF_ZERO_SQRT, 
                                                   NULL, 
                                                   NULL,
                                                   SM_SR_SINGLE, 
                                                   sSolutions1));
          if ( sSolutions1.GetSize() > 0 )
            {
              SmStatus eStat = pOrigRailCurve[0]->Evaluate(
                      sSolutions1[0].m_vStart[0], 1, TRUE, sPV );
              if ( eStat == SM_SUCCESS )
                {
                  p3DDir = &( sPV[1] );
                }
            }

          // Recreate the fillet surface using the refined fillet points
          // to get the fillet surface to go exactly through the points at
          // edges and vertices.
          SER( pCurrFilletGeom->ReCalcFilletGeom( bIsAnalytic, p3DDir ));

          // Compute support edges for cliffs first
          for (ULONG jj=0; jj<sCliffSupportEdges.GetSize(); jj++)
            {
              // Computes the 3D Geometry (Curve and endVertexPoint) for pSupportEdge - existing Curve and UVTrimCurves are deleted
              SER( sm_ComputeSupportEdge( m_crContext,
                                          sCliffSupportEdges[jj],
                                          pCurrFilletGeom->GetFilletSurface(),
                                          m_dThisApproxTol3d, 
                                          m_dThisAngTolRad, 
                                          iDebugLevel) ) ;
            }
          sCliffSupportEdges.ReSet();

          // Then trim all rails
          SER( pCurrFilletGeom->TrimRailCurves() );
          if ( eFGType == SM_FG_TANGENT_ROLLOVER )
            {
              // Restore original offsets data
              SetOffsetsData( sOrigOffsets );
            }
        } // end if tangent, or processing last geom

      // Set the rail and geom pointers for the next segment.
      // If we split for this segment, then they'll be updated,
      // else they stay the same as this time.
      // Also set a flag for the next segment, whether they were split.

      if ( bProcessingLastGeom )
        { break; }

      bSplitRailLastTime[0] = bSplitRailLastTime[1] = FALSE;
      if ( pNewFilletGeom )
        {
          pCurrFilletGeom = pNewFilletGeom;
          pRail[0] = pNewFilletGeom->GetRail(0);
          pRail[1] = pNewFilletGeom->GetRail(1);
          bSplitRailLastTime[0] = bSplitRailLastTime[1] = TRUE;
        }
      else
        {
          if ( pNewRail1 )
            {
              pRail[0] = pNewRail1;
              bSplitRailLastTime[0] = TRUE;
            }
          if ( pNewRail2 )
            {
              pRail[1] = pNewRail2;
              bSplitRailLastTime[1] = TRUE;
            }
        }

    } // end while walking rails segment by segment

  // all done
  return SM_SUCCESS;

} // end SmFilletSolver::ProcessRollOver

/*******************************************************************//**
PURPOSE: Tell the FilletExecutive about a fillet error.

USAGE ---
    This could arise from offset-surface intersection failure.
***********************************************************************/
void SmFilletSolver::RecordFilletError
 (SmFilletErrorType       eType,    // in :
  SmTArray<SmTsectPnt*> & rapTSPs,  // in :
  const TCHAR           * cComment) // in :
{
  SmEdge *pFilletedEdge = GetEdgeuse(0)->GetEdge();

  // Collect all existing fillet surfaces and rail curves
  SmTArray<SmFilletGeom*> apFilletGeoms;
  GetFilletGeoms( apFilletGeoms );
  ULONG lNumGeoms = apFilletGeoms.GetSize();
  
  // gwc: I think pSideSurface no longer has to be a BBSplineSurface - replace 1 line
  // rm : SmTArray<SmBSplineSurface*> apFilSurfs;
  SmTArray<SM_FILLETSURF_TYPE*> apFilSurfs;
  SmTArray<SmBSplineCurve*>   apFilCurvs;
  SmTArray<SmFilletEdge*>     apFilEdges;
  for ( ULONG i = 0; i < lNumGeoms; i++ )
    {
      SmFilletGeom * pThisGeom = apFilletGeoms[i];
      if ( pThisGeom->GetFilletSurface()     != NULL ) { apFilSurfs.AddUnique( pThisGeom->GetFilletSurface() ); }
      if ( pThisGeom->GetCenterLineCurve()   != NULL ) { apFilCurvs.AddUnique( pThisGeom->GetCenterLineCurve() ); }
      if ( pThisGeom->GetRail(0)->GetCurve() != NULL ) { apFilEdges.AddUnique( pThisGeom->GetRail(0) ); }
      if ( pThisGeom->GetRail(1)->GetCurve() != NULL ) { apFilEdges.AddUnique( pThisGeom->GetRail(0) ); }
    }

  // Just look at one TsectPnt position for now
  SmPoint3d sTSPPos;

  if ( rapTSPs.GetSize() > 0 )
    {
      if ( rapTSPs[0]->CrvPos().IsInitialized() )
          sTSPPos = rapTSPs[0]->CrvPos();
    }

  // pass the call along to FilletExecutive 
  this->GetFilletExecutive()->NoteFilletError(eType, 
                                              GetStatus(),
                                              pFilletedEdge, 
                                              apFilSurfs, 
                                              apFilCurvs, 
                                              apFilEdges, 
                                              sTSPPos, 
                                              (TCHAR *)cComment);
} // end SmFilletSolver::RecordFilletError

/*******************************************************************//**
PURPOSE:
NOTES:
***********************************************************************/
SmBoolean SmFilletSolver::IsKindOf( SM_TYPE t ) const
{
  return ((SmFilletSolver_TYPE == t) ? TRUE : SmObject::IsKindOf( (t) ));
}

/****************************************************************
PURPOSE: Debug dump.
****************************************************************/
void SmFilletSolver::Dump() const { DumpLevel( 0, NULL ); }

/****************************************************************
PURPOSE: Debug dump.
****************************************************************/
void SmFilletSolver::DumpLevel( int iDebugLevel, const TCHAR * cpMsg ) const
{
  TCHAR sBuff[SM_TBLOCK_SIZE];
  ULONG i, lNum;

  // message
  if(cpMsg != NULL) smos_WriteBuffer(cpMsg) ;

  // header
  smos_sprintf( sBuff, _T("\nDump of SmFilletSolver [%4ld]  Status: %d\n"),
                     this->FindIndexInExec(), m_eStatus) ;
  smos_WriteBuffer(sBuff);

  // tols
  smos_sprintf( sBuff, _T("  Approx tol %lf, Angle tol %lf, Tangency tol %lf\n"),
                     m_dThisApproxTol3d, m_dThisAngTolRad, m_dTangencyTolerance) ;
  smos_WriteBuffer(sBuff);

  // Blend Factor
  smos_sprintf( sBuff, _T("  Edge Blend Factor %lf\n"), m_dEdgeBlendFactor) ;
  smos_WriteBuffer(sBuff);

  // Offset Surface distances
  for(i=0;i<2;i++)
    {
      double dOffset = GetSurface(i)->GetOffsetDistance() ;

      if(smos_Fabs(dOffset) < SM_EFF_ZERO )
        { dOffset = m_dSavedOffsets[i] ; }
      smos_sprintf( sBuff, _T("  Base Surface %2ld, Offset Distance %lf\n"), i+1, dOffset );
      smos_WriteBuffer(sBuff);
    }

  // Extensions and Reverse Trim
  smos_sprintf(sBuff, _T("  Extend Before:[%s]  Extend After:[%s]  Reverse Trim:[%s]\n"),
             m_bExtendBefore ? _T("TRUE") : _T("FALSE"), 
             m_bExtendAfter  ? _T("TRUE") : _T("FALSE"), 
             m_bReverseTrim  ? _T("TRUE") : _T("FALSE") );
  smos_WriteBuffer(sBuff);

  // Num of internal G1 knots
  lNum = m_vG1Knots.GetSize(); 
  smos_sprintf( sBuff, _T("  %ld Internal G1 knots  "), lNum );
  smos_WriteBuffer(sBuff);

   // G1 Knots
  for(i=0;i<lNum;i++)
    {
      smos_sprintf( sBuff, _T("  %lf"), m_vG1Knots[i] );
      smos_WriteBuffer(sBuff);
    }
  smos_sprintf( sBuff,_T("%s"), _T("\n") );
  smos_WriteBuffer(sBuff);

  // Num of FilletGeoms
  lNum = m_vFilletGeoms.GetSize(); 
  smos_sprintf( sBuff, _T("  %ld FilletGeoms\n"), lNum );
  smos_WriteBuffer(sBuff);

  // Fillet Geoms
  for(i=0;i<lNum;i++)
    {
      m_vFilletGeoms[i]->DumpLevel( iDebugLevel, _T("    Fillet Geom:\n") );
    }

}  // end SmFilletSolver::DumpLevel()

/****************************************************************
PURPOSE: Debug draw.
****************************************************************/
void SmFilletSolver::Draw()
{
#ifdef SM_DEBUG_CODE
  smgfx_SetLook( 1,2, 1,0,0 );
  if ( m_pEdgeuses[0] != NULL ) { m_pEdgeuses[0]->Draw(); } sm_GraphicsLoop();
  smgfx_SetLook( 1,2, 0,0,1 );
  if ( m_pEdgeuses[1] != NULL ) { m_pEdgeuses[1]->Draw(); } sm_GraphicsLoop();

  ULONG ii, lNumFGs = m_vFilletGeoms.GetSize();
  for ( ii = 0; ii < lNumFGs; ii++ )
    { m_vFilletGeoms[ii]->Draw(); }
#endif
}  // end SmFilletSolver::Draw()

/*******************************************************************//**
PURPOSE: Return the index of a given FilletGeom in this FilletSolver.

NOTES: Used in debugging.
    Returns error 9999 if the given FilletGeom is not in our list.
***********************************************************************/
ULONG SmFilletSolver::FindFilGeomIndex( const SmFilletGeom *pFilGeom ) const
{
  ULONG lIdx;
  if ( m_vFilletGeoms.FindElement( SM_CONST_CAST(SmFilletGeom*, pFilGeom), lIdx ))
    { return lIdx; }
  return 9999;  // Error return
}
/*******************************************************************//**
PURPOSE: Return the index of this FilletSolver in our FilletExecutive.

NOTES: Used in debugging.
    Returns error 9999 if 'this' is not in our FilletExecutive.
***********************************************************************/
ULONG SmFilletSolver::FindIndexInExec() const
{
  if ( m_pExecutive == NULL ) { return 9999; }
  return m_pExecutive->FindSolverIndex( this );
}

#ifdef SM_DEBUG_CODE
/****************************************************************
PURPOSE: Debug level.
****************************************************************/
int SmFilletSolver::DebugLevel()
  { return m_pExecutive->DebugLevel(); }
#endif // SM_DEBUG_CODE


/****************************************************************
PURPOSE: Compute fillet radius given UV curve of rails and the curve
            parameter. In addition, we need a 3d point and a tangent vector
            along the filleted edge where the radius is calculated.

NOTES:
****************************************************************/
static SmStatus sm_ComputeFilletRadius
  (SmFilletSolver * pFS,         // in : FilletSolver containing fillet base surfaces
   SmBSplineCurve * pUVCrv1,     // in : rail UVCurve on surface0  
   SmBSplineCurve * pUVCrv2,     // in : rail UVCurve on surface1 
   double           dT,          // in : Curve param value of interest
   SmVector3d       sPV[],       // in : sPV[0] = pt on offset-surface centerline xsect curve        
                                 //      sPV[1] = tang at that offset-surface centerline xsect curve 
   double         * pdRadius)    // out: fillet radius needed between the pUVCrv1[dT] and pUVCrv2[dT]
{
  // fillet surface natural (not extended) UVDomains
  SmExtent2d sDomain1 = pFS->GetSurface(0)->GetBaseSurface()->GetNaturalUVDomain();
  SmExtent2d sDomain2 = pFS->GetSurface(1)->GetBaseSurface()->GetNaturalUVDomain();

  // locals
  SmPoint3d  sPnt;
  SmBoolean  bFoundSolution;
  SmTsectPnt sTsectPnt;
  
  // surface0 UVPoint 
  SER(pUVCrv1->EvaluatePoint(dT,sPnt));
  SmVector2d sUV1(sPnt.x,sPnt.y);

  // surface1 UVPoint
  SER(pUVCrv2->EvaluatePoint(dT,sPnt));
  SmVector2d sUV2(sPnt.x,sPnt.y);

  // find filletPoint on given plane satisfying geometry requirements 
  // implemented in derived SmFilletSolver class
  SER(pFS->PointOnPlaneSolve(sPV[0],
                             sPV[1],
                             sDomain1,
                             sDomain2,
                             sUV1,
                             sUV2,
                             bFoundSolution,
                             sTsectPnt));
  if (!bFoundSolution) SER(SM_ERR);

  // extract radius at this fillet point for the 
  SER(pFS->SetupOffsetValues(sTsectPnt, pdRadius));

  // all done
  return SM_SUCCESS;

} // end sm_ComputeFilletRadius

/*******************************************************************//**
PURPOSE: Handle a self intersection of fillets by doing a step
  back and blend of the curves then remake the surface.

NOTES: Joins together two pieces of a fillet.  
 The prev piece is stored in the pFilletGeom
 The next piece is stored in the input r3DCurves rUVCurves1, and rUVCurves2 arrays.

 Trims prev and next curves back from the join point by about 1.5 * the fillet radius, 
 and adds a blend curve accross the gap.  The intent is to cut out a a section of fillet
 curve that contains a cusp.

METHOD: choses a stepback distance from the end of the prev piece of the

SIDE EFFECTS: rebuilds the centerline, rail3D and railUV curves in pFilletGeom
***********************************************************************/
SmStatus SmMakeCurveBlendSIH::SmoothOutIntersection
  (SmFilletGeom              * pFilletGeom,   // i/o: gets rebuilt centerline, rail3D and railUV curves
   SmTArray<SmBSplineCurve*> & r3DCurves,     // i/o: in : r3DCurves[0]  = next piece of fillet centerline
                                              //           r3DCurves[1]  = next piece of Rail Curve surface0
                                              //           r3DCurves[2]  = next piece of Rail Curve surface1
                                              //      out: r3dCurves[0]  = joined prev, blend, next pieces of fillet centerline
                                              //           r3DCurves[1]  = joined prev, blend, next pieces of Rail Curve surface0
                                              //           r3DCurves[2]  = joined prev, blend, next pieces of Rail Curve surface1
   SmTArray<SmBSplineCurve*> & rUVCurves1,    // i/o: in : rUVCurves1[0] = next piece of rail UVCurve surface0
                                              //      out: rUVCurves1[0] = joined prev, blend, next pieces of UVCurve surface0
   SmTArray<SmBSplineCurve*> & rUVCurves2)    // i/o: in : rUVCurves2[0] = next piece of rail UVCurve surface1
                                              //      out: rUVCurves1[0] = joined prev, blend, next pieces of UVCurve surface0
{
  // check input
  NER(pFilletGeom); 

  // locals
  SmBSplineCurve * pNextCenterCurve ;
  SmBSplineCurve * apNextRailCurves[2] ;
  SmBSplineCurve * apNextUVCurves[2];
  SmFilletEdge   * pPrevRail[2];
  SmFace         * pFace1, * pFace2 ;
  SmSurface      * pSurf1, * pSurf2 ; 
  SmExtent2d       sDomain1, sDomain2 ;   // natural (not extended) surface domains

  // Get next piece of fillet geometry curves from the input arrays
  pNextCenterCurve     = r3DCurves[0] ; SmObjDelete sClean1(pNextCenterCurve) ;  // next centerline piece
  apNextRailCurves[0]  = r3DCurves[1] ; SmObjDelete sClean2(apNextRailCurves[0]) ;   // next rail0 piece
  apNextRailCurves[1]  = r3DCurves[2] ; SmObjDelete sClean3(apNextRailCurves[1]) ;   // next rail1 piece
  apNextUVCurves[0]    = rUVCurves1[0] ; SmObjDelete sClean4(apNextUVCurves[0]) ;    // next railUV0 piece
  apNextUVCurves[1]    = rUVCurves2[0] ; SmObjDelete sClean5(apNextUVCurves[1]) ;    // next railUV0 piece

  // previous piece of fillet geometry is in pFilletGeom object
  SmFilletSolver * pFS = pFilletGeom->GetFilletSolver() ;
  NER(pFS);

  pPrevRail[0] = pFilletGeom->GetRail(0);
  pPrevRail[1] = pFilletGeom->GetRail(1);
  pFace1       = pPrevRail[0]->GetOriginalFace(); NER(pFace1);
  pFace2       = pPrevRail[1]->GetOriginalFace(); NER(pFace2);
  pSurf1       = pFace1->GetSurface(); NER(pSurf1);
  pSurf2       = pFace2->GetSurface(); NER(pSurf2);
  sDomain1     = pSurf1->GetNaturalUVDomain();
  sDomain2     = pSurf2->GetNaturalUVDomain();

  // This method assumes that parameter values are increasing from the
  // existing curve to the new curves.

  SmExtent1d sPrevTrimIvl = pFilletGeom->GetCenterLineCurve()->GetNaturalInterval();
  SmExtent1d sNextTrimIvl = pNextCenterCurve->GetNaturalInterval();

  SmBoolean bFoundIntersection = FALSE;

  // Find intersection between corresponding UV rail curves
  for (ULONG ii=0; ii<2 && !bFoundIntersection; ii++) 
    {
      // prev and next UVTrimCurve pieces
      SmFilletEdge    * pPrevRailEdge = pFilletGeom->GetRail(ii);
      SmFilletEdgeuse * pPrevPrimEU   = (SmFilletEdgeuse*)pPrevRailEdge->GetPrimaryEdgeuse();
      SmFilletEdgeuse * pPrevMateEU   = (SmFilletEdgeuse*)pPrevPrimEU->GetMate();

      SmBSplineCurve  * pPrevUVCurve = pPrevMateEU->GetUVTrimCurve();
      SmBSplineCurve  * pNextUVCurve = apNextUVCurves[ii];

      // solution locals
      SmSolution sSData[16];
      SmSolutionArray sSolutions(16,sSData);
#ifdef SM_DEBUG_CODE
SmBoolean bDebugMe3 = FALSE;
      if (bDebugMe3) 
        {
          smgfx_SetColor(1,0,0); pPrevRailEdge->GetCurve()->DrawWDeriv(pPrevRailEdge->GetCurve()->GetNaturalInterval(),0); sm_GraphicsLoop();
          smgfx_SetColor(0,0,1);apNextRailCurves[ii]->DrawWDeriv(apNextRailCurves[ii]->GetNaturalInterval(),0); sm_GraphicsLoop();
          sm_GraphicsLoop();
        }
#endif // SM_DEBUG_CODE

      // Use this little half interval trick to eliminate the possibility of
      // finding closed intersections.  It may be a little dangerous.
      SmExtent1d sHalfInterval = pPrevUVCurve->GetNaturalInterval();
      (void)sHalfInterval.SetMinMax(sHalfInterval.Evaluate(0.5),sHalfInterval.GetMax()); // SolidCraft: from the middle to the end, it cannot fail

      // intersect prev UVTrimCurve piece with next piece
      SER(pPrevUVCurve->GlobalCurveIntersect(sHalfInterval,
                                            *pNextUVCurve,
                                             pNextUVCurve->GetNaturalInterval(),
                                             SM_EFF_ZERO_SQRT,
                                             sSolutions));

      // for every prev/next piece UVTrimCurve intersection 
      for (ULONG jj=0; jj<sSolutions.GetSize(); jj++) 
        {
          SmSolution & rSol    = sSolutions[jj];
          double       dPrevT  = rSol.m_vStart[0];
          double       dNextT  = rSol.m_vStart[1];

          // prev piece UVPoint
          SmVector3d sPrevPV[2], sNextPV[2] ;
          SER(pPrevRailEdge->GetCurve()->Evaluate(dPrevT,1,TRUE,sPrevPV));

          // Compute stepback distance using derivative of curve
          // This is just an approximation.
          double dRadius;
          SER(sm_ComputeFilletRadius(pFS,                // in : FilletSolver containing fillet base surfaces                
                                     apNextUVCurves[0],  // in : rail UVCurve on surface0                                    
                                     apNextUVCurves[1],  // in : rail UVCurve on surface1                                    
                                     dNextT,             // in : Curve param value of interest                               
                                     sPrevPV,            // in : sPV[0] = pt on offset-surface centerline xsect curve        
                                                         //      sPV[1] = tang at that offset-surface centerline xsect curve 
                                     &dRadius));         // out: fillet radius needed between the pUVCrv1[dT] and pUVCrv2[dT]
          
          // step distance                                          
          double dStepDistance = m_dStepBackDistFactor * dRadius;
          
          // step back param distance                                          
          double dTanLeng = sPrevPV[1].Length();
          if (SM_IS_ZERO(dTanLeng)) SER(SM_ERR);
          double dPrevStep = dPrevT - dStepDistance/dTanLeng;

          // If the step goes too far then step back 1/2 of distance to start of curve.
          if (dPrevStep <= sPrevTrimIvl.GetMin()) 
            {
              dPrevStep = (dPrevT + sPrevTrimIvl.GetMin()) / 2.0;
            }

          // step forward param distance
          SER(apNextRailCurves[ii]->Evaluate(dNextT,1,TRUE,sNextPV));
          dTanLeng = sNextPV[1].Length();
          if (SM_IS_ZERO(dTanLeng)) SER(SM_ERR);
          double dNextStep = dNextT + dStepDistance/dTanLeng;
          if (dNextStep >= sNextTrimIvl.GetMax()) 
            {
              dNextStep = (dNextT + sNextTrimIvl.GetMax()) / 2.0;
            }

          // adjust the Prev and Next curve intervals to the step back/forward points
          if (dPrevStep < sPrevTrimIvl.GetMax()) { (void)sPrevTrimIvl.SetMinMax(sPrevTrimIvl.GetMin(),dPrevStep); } // SolidCraft: ordered by the test, it cannot fail
          if (dNextStep > sNextTrimIvl.GetMin()) { (void)sNextTrimIvl.SetMinMax(dNextStep,sNextTrimIvl.GetMax()); } // SolidCraft: ordered by the test, it cannot fail
          bFoundIntersection = TRUE;
          
        } // end iter (jj) every next/prev UVtrimCurve intersection solution
    } // end iter (ii) both rail curves

  // when Prev and Next UVTrim curves did not intersect
  if (!bFoundIntersection) 
    {
      double dMin = SM_BIG_DOUBLE;
      ULONG  lMin = 0;
      double dPrevT = sPrevTrimIvl.GetMin();
      double dNextT = sNextTrimIvl.GetMax();

      // for both rail curves - 
      for (ULONG kkk=0; kkk<2; kkk++) 
        {
          SmCurve * pPrevRailCurve = pFilletGeom->GetRail(kkk)->GetCurve();

          SmSolutionArray sSolutions;
          // Use this little half interval trick to eliminate the possibility of
          // finding closed intersections.  It may be a little dangerous.
          SmExtent1d sHalfInterval = pPrevRailCurve->GetNaturalInterval();
          (void)sHalfInterval.SetMinMax(sHalfInterval.Evaluate(0.5),sHalfInterval.GetMax()); // SolidCraft: from the middle to the end, it cannot fail

          // find closest point between prev and Next rail curve pieces
          SER(pPrevRailCurve->GlobalCurveSolve(sHalfInterval,
                                              *apNextRailCurves[kkk],
                                               apNextRailCurves[kkk]->GetNaturalInterval(),
                                               SM_SO_MINIMIZE,
                                               SM_EFF_ZERO_SQRT,
                                               NULL,NULL,
                                               SM_SR_SINGLE,
                                               sSolutions));

          // save the closest points on prev/next rail curve pieces 
          if (sSolutions[0].m_vStart.m_dSolutionValue < dMin) 
            {
              dPrevT = sSolutions[0].m_vStart[0];
              dNextT = sSolutions[0].m_vStart[1];
#ifdef SM_DEBUG_CODE
SmBoolean bDebugMe4 = FALSE;
              if (bDebugMe4) 
                {
                  SmPoint3d sPnt1, sPnt2;
                  pPrevRailCurve->EvaluatePoint(dPrevT, sPnt1);
                  apNextRailCurves[kkk]->EvaluatePoint(dNextT, sPnt2);

                  sPnt1.Draw(); sm_GraphicsLoop();
                  sPnt2.Draw(); sm_GraphicsLoop();
                  sm_GraphicsLoop();
                }
#endif // SM_DEBUG_CODE

              dMin = sSolutions[0].m_vStart.m_dSolutionValue;
              lMin = kkk;
            }
        } // end iter both rail curves

      // prev railCurve nearest 3d point
      SmVector3d sPrevPV[2];
      SmFilletEdge *pPrevRailEdge = pFilletGeom->GetRail(lMin);
      SER(pPrevRailEdge->GetCurve()->Evaluate(dPrevT,1,TRUE,sPrevPV));

      // Compute stepback distance using derivative of curve
      // This is just an approximation.
      double dRadius;
      SER(sm_ComputeFilletRadius(pFS,                // in : FilletSolver containing fillet base surfaces                
                                 apNextUVCurves[0],  // in : rail UVCurve on surface0                                    
                                 apNextUVCurves[1],  // in : rail UVCurve on surface1                                    
                                 dNextT,             // in : Curve param value of interest                               
                                 sPrevPV,            // in : sPrevPV[0] = pt on offset-surface centerline xsect curve        
                                                     //      sPrevPV[1] = tang at that offset-surface centerline xsect curve 
                                 &dRadius));         // out: fillet radius needed between the pUVCrv1[dT] and pUVCrv2[dT]

      // step distance scaled by fillet radius
      double dStepDistance = m_dStepBackDistFactor * dRadius;

      if (dMin > 3.0*dStepDistance) 
        {
          return SM_ERR;
        }

      // stap back param distance
      double dTanLeng = sPrevPV[1].Length();
      if (SM_IS_ZERO(dTanLeng)) SER(SM_ERR);
      double dPrevStep = dPrevT - dStepDistance/dTanLeng;

      // If the step goes too far then step back 1/2 of distance to start of curve.
      if (dPrevStep <= sPrevTrimIvl.GetMin()) 
        {
          dPrevStep = (dPrevT + sPrevTrimIvl.GetMin()) / 2.0;
        }

      // step forward param distance
      SmVector3d pNextPV[2];
      SER(apNextRailCurves[lMin]->Evaluate(dNextT,1,TRUE,pNextPV));
      dTanLeng = pNextPV[1].Length();
      if (SM_IS_ZERO(dTanLeng)) SER(SM_ERR);
      double dNextStep = dNextT + dStepDistance/dTanLeng;

      // if step goes too face, set step to 1/2 distance to end of curve
      if (dNextStep >= sNextTrimIvl.GetMax()) 
        {
          dNextStep = (dNextT + sNextTrimIvl.GetMax()) / 2.0;
        }

      // Now see if we are actually close enough to do this
      if (dPrevStep < sPrevTrimIvl.GetMax()) { (void)sPrevTrimIvl.SetMinMax(sPrevTrimIvl.GetMin(),dPrevStep); } // SolidCraft: ordered by the test, it cannot fail
      if (dNextStep > sNextTrimIvl.GetMin()) { (void)sNextTrimIvl.SetMinMax(dNextStep,sNextTrimIvl.GetMax()); } // SolidCraft: ordered by the test, it cannot fail
    
    } // end Prev and Next UVTrim curves did not intersect check

  // Now Trim Back next and prev centerline curves with new intervals
  SER(pNextCenterCurve->Trim(sNextTrimIvl));                   // may snap sIvl by tol to existing knots
  SER(pFilletGeom->GetCenterLineCurve()->Trim(sPrevTrimIvl));  // may snap sIvl by tol to existing knots

  // locals
  SmBSplineCurve * pCenterBlend    = NULL ;
  SmBSplineCurve * pPrevCenterLine = pFilletGeom->GetCenterLineCurve();
  SmTArray<SmPoint3d>  sPoints;   // Centerline positions
  SmTArray<SmVector3d> sVectors;  // Centerline 1stDerivs
  SmVector3d sPrevPV[2], sNextPV[2], sPV[2];

  // add Prev trimmed Centerline end point and 1stderiv
  SER(pPrevCenterLine->Evaluate(sPrevTrimIvl.GetMax(),1,TRUE,sPrevPV));
  sPoints.Add(sPrevPV[0]);
  sVectors.Add(sPrevPV[1]);

  // add Next trimmed Centerline start point and 1stDeriv
  SER(pNextCenterCurve->Evaluate(sNextTrimIvl.GetMin(),1,TRUE,sNextPV));
  sPoints.Add(sNextPV[0]);
  sVectors.Add(sNextPV[1]);

  // adjust 1stDeriv mags to optimize curvature distribution of upcoming curve
  SER(SmBSplineCurve::FairBlendCurveDerivatives(sPoints,sVectors,NULL));

  // create blend Centerline piece as a BSpline Curve from PrevEndPoint to NextStartPoint 
  SER(SmBSplineCurve::CreateInterpolatingCurve(pFS->GetCreationContext(),
                                               SM_CP_UNIFORM,
                                               3,3,
                                               sPoints,
                                               sVectors,
                                               NULL,TRUE,
                                               pCenterBlend));
  
  // load sCurves array with Prev, Blend, and Next Centerline pieces
  SmTArray<SmBSplineCurve*> sCurves;
  sCurves.Add(pPrevCenterLine);
  sCurves.Add(pCenterBlend);
  sCurves.Add(pNextCenterCurve);

#ifdef SM_DEBUG_CODE
SmBoolean bDebugMe = FALSE;
  if (bDebugMe) 
    {
      smgfx_SetColor(1,0,0); pPrevCenterLine->DrawWDeriv(pPrevCenterLine->GetNaturalInterval(),0); sm_GraphicsLoop();
      smgfx_SetColor(0,1,0); pCenterBlend->DrawWDeriv(pCenterBlend->GetNaturalInterval(),0); sm_GraphicsLoop();
      smgfx_SetColor(0,1,0); pCenterBlend->DrawPolygon(); sm_GraphicsLoop();
      smgfx_SetColor(0,0,1); pNextCenterCurve->DrawWDeriv(pNextCenterCurve->GetNaturalInterval(),0); sm_GraphicsLoop();
      sm_GraphicsLoop();
    }
#endif // SM_DEBUG_CODE

  // Create these blend curves with 10 knots along the U direction
  SmTArray<double> sNewKnots;
  sNewKnots.Add(0.1);
  sNewKnots.Add(0.2);
  sNewKnots.Add(0.3);
  sNewKnots.Add(0.4);
  sNewKnots.Add(0.5);
  sNewKnots.Add(0.6);
  sNewKnots.Add(0.7);
  sNewKnots.Add(0.8);
  sNewKnots.Add(0.9);

  // add internal knots to CenterBlend curve
  SER( pCenterBlend->InsertKnots(sNewKnots) );

  double dBlendArcLeng = 0.0;
  SmObjDelete sClean21(pCenterBlend);

  // set output - trimmed prev piece of centerline
  r3DCurves.ReSet();
  r3DCurves.Add(pPrevCenterLine);

  // locals
  SmTArray<SmBSplineCurve*> sNewCurves;
  SmTArray<SmBSplineCurve*> sUVBlends;
  SmTArray<SmBSplineCurve*> s3DBlends;
  sNewCurves.Add(pPrevCenterLine);

  // sets the surface offsets back to zero making them equivalent to the base surfaces.
  pFS->SetSurfacesToZeroOffset();

  // for both rail curves
  for (ULONG ii=0; ii<2; ii++) 
    {
      SmFilletEdge    * pPrevRailEdge   = pFilletGeom->GetRail(ii);
      SmFilletEdgeuse * pPrevPrimEU     = (SmFilletEdgeuse*)pPrevRailEdge->GetPrimaryEdgeuse();
      SmFilletEdgeuse * pPrevMateEU     = (SmFilletEdgeuse*)pPrevPrimEU->GetMate();
      SmBSplineCurve  * pPrevUVCurve    = pPrevMateEU->GetUVTrimCurve(); NER(pPrevUVCurve);
      SmBSplineCurve  * pPrevRailCurve  = SM_CAST_PTR(SmBSplineCurve,pPrevRailEdge->GetCurve()); NER(pPrevRailCurve);

      // Trim PrevRailEdge(ii), UVCurve, and RailCurve intervals
      pPrevRailEdge->SetInterval(sPrevTrimIvl);
      SER(pPrevUVCurve->Trim(sPrevTrimIvl));      // may snap sIvl by tol to existing knots
      SER(pPrevRailCurve->Trim(sPrevTrimIvl));    // may snap sIvl by tol to existing knots

      // Trim NextRail(ii) UVCurve and RailCurve intervals
      SER(apNextUVCurves[ii]->Trim(sNextTrimIvl));     // may snap sIvl by tol to existing knots
      SER(apNextRailCurves[ii]->Trim(sNextTrimIvl));   // may snap sIvl by tol to existing knots

      // Create UV Blend curves and their 3D appproximations
      sPoints.ReSet();
      sVectors.ReSet();

      // Add PrevUVCurve endPoint and 1stDeriv
      SER(pPrevUVCurve->Evaluate(sPrevTrimIvl.GetMax(),1,TRUE,sPV));
      sPoints.Add(sPV[0]);
      sVectors.Add(sPV[1]);

      // Add NextUVCurve startPoint and 1stDeriv
      SER(apNextUVCurves[ii]->Evaluate(sNextTrimIvl.GetMin(),1,TRUE,sPV));
      sPoints.Add(sPV[0]);
      sVectors.Add(sPV[1]);

      SmBSplineCurve *pUVBlend = NULL ;
      
      // adjust 1stDeriv mags to optimize curvature distribution of upcoming curve
      SER(SmBSplineCurve::FairBlendCurveDerivatives(sPoints,sVectors,NULL));
      
      // create blend Centerline piece as a BSpline Curve from PrevEndPoint to NextStartPoint 
      SER(SmBSplineCurve::CreateInterpolatingCurve(pFS->GetCreationContext(),
                                                   SM_CP_UNIFORM,
                                                   2,3,
                                                   sPoints,
                                                   sVectors,
                                                   NULL,TRUE,
                                                   pUVBlend));
      SER( pUVBlend->InsertKnots(sNewKnots) );

      // add rail(ii) blend UVCurve to array
      sUVBlends.Add(pUVBlend);


      // lift UVCurve to 3D curve and add to Blends array
      SmBSplineCurve * p3DLift = NULL ;
      double dMaxDist;
      SER(pFS->GetSurface(ii)->LiftCurve(pFS->GetCreationContext(),
                                         pFS->GetSurface(ii)->GetNaturalUVDomain(),
                                         *pUVBlend,pUVBlend->GetNaturalInterval(),
                                         pFS->GetThisApproxTol3d()*10.0,
                                         dMaxDist,
                                         p3DLift));
      s3DBlends.Add(p3DLift);

      // accumulate arc length
      double dArcLeng = p3DLift->ApproximateLength(p3DLift->GetNaturalInterval(),20);
      dBlendArcLeng += dArcLeng;

    } // end iter (ii) both rail curves

  // parameterize blend on arclenght
  SmExtent1d sArcLengIvl(0.0,dBlendArcLeng/2.0);
  SER(pCenterBlend->EditParameterization(sArcLengIvl));

  // join 3 curves and place output in pCenterBlend
  SER(SmBSplineCurve::CreateByJoining(pFS->GetCreationContext(),sCurves,NULL,pCenterBlend));

  // replace Prev CenterLine piece Nurb with joint Nurb curve
  gw_CURVE *pNurb = pCenterBlend->GetOrCreateGwNurbPointer();
  SER( pPrevCenterLine->SetFromGwNurb(0, pNurb) );

  // for both rail curves
  for (ULONG iii=0; iii<2; iii++) 
    {
      SmFilletEdge    * pPrevRailEdge  = pFilletGeom->GetRail(iii);
      SmFilletEdgeuse * pPrevPrimEU    = (SmFilletEdgeuse*)pPrevRailEdge->GetPrimaryEdgeuse();
      SmFilletEdgeuse * pPrevMateEU    = (SmFilletEdgeuse*)pPrevPrimEU->GetMate();
      SmBSplineCurve  * pPrevUVCurve   = pPrevMateEU->GetUVTrimCurve(); NER(pPrevUVCurve);
      SmBSplineCurve  * pPrevRailCurve = SM_CAST_PTR(SmBSplineCurve,pPrevRailEdge->GetCurve()); NER(pPrevRailCurve);

      SmBSplineCurve * pUVBlend = sUVBlends[iii];
      SmObjDelete pClean7(pUVBlend);

      // adjust UVblend parameterization to arc length
      SER(pUVBlend->EditParameterization(sArcLengIvl));

      // adjust lifted curve parameritization to arc length
      SmBSplineCurve * p3DLift = s3DBlends[iii];
      SmObjDelete sClean8(p3DLift);
      SER(p3DLift->EditParameterization(sArcLengIvl));

      // load sCurves array with Prev, Blend, and Next UVCurves
      sCurves.ReSet();
      sCurves.Add(pPrevUVCurve);
      sCurves.Add(pUVBlend);
      sCurves.Add(apNextUVCurves[iii]);

      // join the 3 curves into pUVBlend
      SER(SmBSplineCurve::CreateByJoining(pFS->GetCreationContext(),sCurves,NULL,pUVBlend));
      SmObjDelete sClean9(pUVBlend);

      // replace PrevUVCurve with pUVBlend NURB definition
      gw_CURVE *pUVNurb = pUVBlend->GetOrCreateGwNurbPointer();
      SER( pPrevUVCurve->SetFromGwNurb(0, pUVNurb) );

      // add to sNewCurves array
      sNewCurves.Add(pPrevUVCurve);

      // set output
      if (iii==0) 
        {
          rUVCurves1.ReSet();
          rUVCurves1.Add(pPrevUVCurve);
        }
      else 
        {
          rUVCurves2.ReSet();
          rUVCurves2.Add(pPrevUVCurve);
        }

      // lift railUVCurve and add to to output
      double dMaxDist;
      SER(pFS->GetSurface(iii)->LiftCurve(pFS->GetCreationContext(),
                                          pFS->GetSurface(iii)->GetNaturalUVDomain(),
                                          *pPrevUVCurve,pPrevUVCurve->GetNaturalInterval(),
                                          pFS->GetThisApproxTol3d()*10.0,
                                          dMaxDist,
                                          p3DLift));
      SmObjDelete sClean10(p3DLift);

      gw_CURVE *p3DNurb = p3DLift->GetOrCreateGwNurbPointer();
      SER( pPrevRailCurve->SetFromGwNurb(0, p3DNurb) );

      r3DCurves.Add(pPrevRailCurve);
      sNewCurves.Add(pPrevRailCurve);

#ifdef SM_DEBUG_CODE
SmBoolean bDebugMe2 = FALSE;
      if (bDebugMe2) 
        {
          p3DLift->Dump();
          pUVBlend->Dump();
          pPrevRailCurve->Dump();

          smgfx_SetColor(0,0,1); p3DLift->DrawWDeriv(p3DLift->GetNaturalInterval(),0); sm_GraphicsLoop();
          smgfx_SetColor(0,0,1); p3DLift->DrawPolygon(); sm_GraphicsLoop();
          smgfx_SetColor(0,0,1); pUVBlend->DrawWDeriv(pUVBlend->GetNaturalInterval(),0); sm_GraphicsLoop();
          smgfx_SetColor(0,0,1); pUVBlend->DrawPolygon(); sm_GraphicsLoop();
          smgfx_SetColor(1,1,0); pPrevUVCurve->DrawWDeriv(pPrevUVCurve->GetNaturalInterval(),0); sm_GraphicsLoop();
          smgfx_SetColor(1,1,0); pPrevUVCurve->DrawPolygon(); sm_GraphicsLoop();
          smgfx_SetColor(0,1,1); pPrevRailCurve->DrawWDeriv(pPrevUVCurve->GetNaturalInterval(),0); sm_GraphicsLoop();
          smgfx_SetColor(0,1,1); pPrevRailCurve->DrawPolygon(); sm_GraphicsLoop();
          sm_GraphicsLoop();
        }
#endif // SM_DEBUG_CODE

    } // end iter both rail curves

  // Make sure all curves are syncronized
  SER(SmBSplineCurve::SyncronizeKnotsOfCurves(sNewCurves,SM_EFF_ZERO*100.0));
  for (ULONG kk=0; kk<sNewCurves.GetSize(); kk++) 
    {
      SmBSplineCurve *pBSC = sNewCurves[kk];
      SER( pBSC->DampenKnotSpacing(2.0,3) );
    }

  // restore FilletSolver offsets
  pFS->ReloadSurfaceOffsets();

  // all done
  return SM_SUCCESS;

} // end SmMakeCurveBlendSIH::SmoothOutIntersection

/*******************************************************************//**
PURPOSE: Detect when a fillet-surface is self-intersecting and find
            parameter values that enclose the self-intersecting portion
            of the surface

NOTES: NOT IMPLMENTED YET - this currently always returns bSelfIntersecting == FALSE
           Needs work depending on self intersections.

METHOD: Check a family of iso-parameter curves for self-intersection.
           For every self-intersection check the distances to the
           the centerLine to determine which side of the self-intersection
           is part of a good-fillet and which part is part of a
           self-intersecting fillet.  When cpFilletSurface is given, the
           rail curves and a family of cpFilletSurface IsoParameter lines
           are checked.  Without cpFilletSurface only the rail curves are checked.

ASSUMES: that cpFilletSurface is a fillet surface created from the
            rail curves so that the rail curves are the min and max
            constant_v isoparameter curves and that the cpFilletSurface
            u direction tracks the input CenterLine.

************************************************************************/
SmStatus SmSelfIntersectionHandler::IsSelfIntersecting
  (const SmContext & crContext,                     // in : context for temporary geometry
   const SmBSplineSurface *cpFilletSurface,         // in : fillet surface to examine, NULL to ignore]
   double dThisApproxTol3d,                         // in : tol used to build pFilletSurface, pFilletGeom->GetFilletSolver()->GetThisApproxTol3d()
   const SmBSplineCurve   *cpCenterLine,            // in : centerline used to define fillet surface
   const SmBSplineCurve   *cpRailCurve1,            // in : min constant_v value cpFilletSurface iso-parameter line
   const SmBSplineCurve   *cpRailCurve2,            // in : max constant_v value cpFilletSurface iso-parameter line
   double                  dFilletRadius,           // in : radius used to define fillet surface
   SmBoolean              &bSelfIntersecting,       // out: TRUE=is self-intersecting, FALSE = not
   SmPoint2d              &rStartPoint,             // NotUsed: out: start of surface self-intersecting region
   SmPoint2d              &rEndPoint)               // NotUsed: out: end   of surface self-intersecting
{
  SM_REF2(rStartPoint, rEndPoint) ; 
  // check state
  NER(cpFilletSurface) ;
  NER(cpCenterLine) ;
  SM_ASSERT(dFilletRadius > 0.0) ;

  // init output
  bSelfIntersecting = FALSE ;

  // locals
  int ii ;
  const int  lInsideCurveCount = 4 ;
  int        lCurveCount       = cpFilletSurface ? lInsideCurveCount + 2 : 2 ;
  SmExtent2d sUVDomain         = cpFilletSurface->GetNaturalUVDomain() ;
  double     dVMin             = sUVDomain.GetMin().y ;
  double     dVMax             = sUVDomain.GetMax().y ;
  double     dDV               = (dVMax - dVMin) / (lInsideCurveCount + 1) ;
  double     dV                = dVMin + dDV ;
  SmBSplineCurve *apIsoCurves[lInsideCurveCount+2] ;
  SmObjDelete     aIsoClean[lInsideCurveCount] ;

  // add the top and bottom bounding isoCurves (the rail curves)
  apIsoCurves[0] = (SmBSplineCurve *)cpRailCurve1 ;
  apIsoCurves[1] = (SmBSplineCurve *)cpRailCurve2 ;

  // when given a cpFilletSurface generate a family of constant_v internal IsoCurves
  for(ii=2;ii<lCurveCount;ii++,dV+=dDV)
    {
      // create/store the next isoparameter curve
      SER(cpFilletSurface->CreateIsoParametricCurve(crContext,
                                                    SM_SP_V,
                                                    dV,
                                                    0.0,
                                                    apIsoCurves[ii])) ;
      aIsoClean[ii-2].SetObj(apIsoCurves[ii]) ;

    } // end iter every internal curve

  // locals for curve self-intersection
  SmSolution sSData[16] ;
  SmSolutionArray sSolutions(16,sSData);

  // for every isoParameterCurve
  for(ii=0;ii<lCurveCount;ii++)
    {
      // get the target curve's interval
      SmExtent1d sInterval = apIsoCurves[ii]->GetNaturalInterval() ;

      // clear solutions array
      sSolutions.ReSet() ;

      // look for self-intersections
      SER(apIsoCurves[ii]->GlobalCurveSelfIntersect(sInterval,
          10.0*dThisApproxTol3d,sSolutions));


//      // GWC::FUTURE when Surface/Curve intersections are sophisticated enough to
//      // handle this situation:
//      if(cpFilletSurface)
//        {
//          SER(cpFilletSurface->GlobalCurveIntersect(sUVDomain, 
//                                                   *apIsoCurves[ii],
//                                                    sInterval,
//                                                    SM_EFF_ZERO,
//                                                    sSolutions)) ;
//        }
//      else
//        {
//          SER(apIsoCurves[ii]->GlobalCurveSelfIntersect(sInterval,
//                                                        10.0*dThisApproxTol3d,
//                                                        sSolutions));
//        }

#ifdef SM_DEBUG_CODE
SmBoolean bDebugMeE = FALSE;
   if (bDebugMeE) 
     {
       smgfx_SetLineWidth(1.0) ;
       smgfx_SetColor(0,1,0);
       apIsoCurves[ii]->Draw();
       sm_GraphicsLoop();
     }
#endif
      // for every solution
      for(ULONG jj=0; jj<sSolutions.GetSize(); jj++)
        {
#ifdef SM_DEBUG_CODE
SmBoolean bDebugMeF = FALSE;
   if (bDebugMeF) 
     {
       // store the parameter values
       SmSolution& rSol = sSolutions[jj];
       double dT1 = rSol.m_vStart[0];
       // double dT2 = rSol.m_vStart[1];

       SmPoint3d sXSectPoint1, sXSectPoint2 ;
       apIsoCurves[ii]->EvaluatePoint(dT1, sXSectPoint1) ;
       // apIsoCurves[ii]->EvaluatePoint(dT2, sXSectPoint2) ;
       smgfx_SetPointSize(4.0) ;
       smgfx_SetColor(1,0,0);
       sXSectPoint1.Draw();
       // sXSectPoint2.Draw();
       sm_GraphicsLoop();
     }
#endif

        } // end iter every solution
    } // end iter every isoParameterCurve

  // all done
  return(SM_SUCCESS) ;

} // end SmMakeSurfaceBlendsSIH::IsSelfIntersecting


/*******************************************************************//**
PURPOSE: Create the Geometry for a set:[Beginning Middle End] of 
         FilletGeoms that replace a single self-intersecting FilletGeom
         with two end segments and a blending segment in the middle.

NOTES: 1. The 3 FilletGeometries are input with topology data but no geometry data.
       2. Any preExisting Curve data in the input FilletGeoms is deleted and
           replaced by new objects.

METHOD ---
  1.

 GWC: This method needs review - 
      I was documenting Curve and UVTrimCurve management
      and saw things that looked suspicious but didn't have the time
      to review this method at that time.
***********************************************************************/
SmStatus SmMakeSurfaceBlendSIH::CreateBlends
  (SmFilletGeom     * pOrigFilletGeom,     // i/o: orig FilletGeom                     - beg of output FilletSequence
   SmFilletGeom     * pBlendFilletGeom,    // i/o: FilletGeom after 1 TopologySplit()  - mid of output FilletSequence
   SmFilletGeom     * pNewFilletGeom,      // i/o: FilletGeom after 2 topologySplits() - end of output FilletSequence
   
   SmBSplineSurface * pNewFilletSurface,   // in : Fillet surface with potential self-intersections
   SmBSplineCurve   * pCenterLine,         // in : Fillet surface center line

   SmBSplineCurve   * pRailCurve1,         // in : 1st 3d rail curve
   SmBSplineCurve   * pRailCurve2,         // in : 2nd 3d rail curve
   SmBSplineCurve   * pUVCurve1,           // in : 1st UV rail curve
   SmBSplineCurve   * pUVCurve2)           // in : 2nd UV rail curve
{
  // check state - must have an OrigFilletGeom
  NER(pOrigFilletGeom);
  SmFilletSolver * pFS = pOrigFilletGeom->GetFilletSolver();
  NER(pFS);

  // locals
  SmBSplineCurve * apRailCurves[2];
  SmBSplineCurve * apUVCurves[2];
  SmFilletEdge   * pOrigRail[2];

  apRailCurves[0] = pRailCurve1;
  apRailCurves[1] = pRailCurve2;
  apUVCurves[0]   = pUVCurve1;
  apUVCurves[1]   = pUVCurve2;
  pOrigRail[0]    = pOrigFilletGeom->GetRail(0);
  pOrigRail[1]    = pOrigFilletGeom->GetRail(1);

  // more locals
  SmFace    * pFace1   = pOrigRail[0]->GetOriginalFace(); NER(pFace1);
  SmFace    * pFace2   = pOrigRail[1]->GetOriginalFace(); NER(pFace2);
  SmSurface * pSurf1   = pFace1->GetSurface(); NER(pSurf1);
  SmSurface * pSurf2   = pFace2->GetSurface(); NER(pSurf2);
  SmExtent2d  sDomain1 = pSurf1->GetNaturalUVDomain();
  SmExtent2d  sDomain2 = pSurf2->GetNaturalUVDomain();

  // This method assumes that parameter values are increasing from the
  // existing curve to the new curves.

  SmBSplineCurve * pOrigCenterLine =  pOrigFilletGeom->GetCenterLineCurve() 
                                    ? pOrigFilletGeom->GetCenterLineCurve() 
                                    : pCenterLine ;
  SmExtent1d sTrimIvl1 = pOrigCenterLine->GetNaturalInterval();
  SmExtent1d sTrimIvl2 = pCenterLine->GetNaturalInterval();

  ULONG     lMin = 0;
  double    dT1  = 0.0 ;
  double    dT2  = 0.0;
  SmBoolean bFoundIntersection = FALSE;

  // Find intersection between corresponding UV rail curves
  for (ULONG i=0; i<2 && !bFoundIntersection; i++)
    {
      // get this rail edge's Edgeuses
      SmFilletEdgeuse * pPrimEU = (SmFilletEdgeuse*)pOrigRail[i]->GetPrimaryEdgeuse();
      SmFilletEdgeuse * pMateEU = (SmFilletEdgeuse*)pPrimEU->GetMate();

      // Get rail edge's UVTrimCurve
      SmBSplineCurve *pCurve1 = pMateEU->GetUVTrimCurve();
      if(!pCurve1)
        { continue; }  //  RCLxx

      // Get the corresponding input UVCurve
      SmBSplineCurve *pCurve2 = apUVCurves[i];
      SmSolution      sSData[16];
      SmSolutionArray sSolutions(16,sSData);

#ifdef SM_DEBUG_CODE
SmBoolean bDebugMe3 = FALSE;
      if (bDebugMe3) 
        {
          smgfx_SetColor(1,0,0); pOrigRail[i]->Draw(); sm_GraphicsLoop();
          smgfx_SetColor(0,0,1); apRailCurves[i]->DrawWDeriv(apRailCurves[i]->GetNaturalInterval(),0); sm_GraphicsLoop();
          sm_GraphicsLoop();
        }
#endif // SM_DEBUG_CODE

      // Use this little half interval trick to eliminate the possibility of
      // finding closed intersections.  It may be a little dangerous.
      SmExtent1d sHalfInterval = pCurve1->GetNaturalInterval();
      (void)sHalfInterval.SetMinMax(sHalfInterval.Evaluate(0.5),sHalfInterval.GetMax()); // SolidCraft: from the middle to the end, it cannot fail

      // look for self-intersections
      SER(pCurve1->GlobalCurveIntersect(sHalfInterval,
                                       *pCurve2,
                                        pCurve2->GetNaturalInterval(),
                                        SM_EFF_ZERO_SQRT,
                                        sSolutions));

      // for every solution
      for (ULONG j=0; j<sSolutions.GetSize(); j++)
        {
          // store the parameter values
          SmSolution & rSol  = sSolutions[j];

          dT1                = rSol.m_vStart[0];
          dT2                = rSol.m_vStart[1];
          lMin               = i;
          bFoundIntersection = TRUE;
        } // end iter both solutions
    } // end iter both rail curves

  // when no UV self intersections were found - check the 3D curve
  double dMin = SM_BIG_DOUBLE;
  if (!bFoundIntersection)
    {
      // for both railcurves
      for (ULONG kkk=0; kkk<2; kkk++)
        {
          SmCurve       * pRailCurve = pOrigRail[kkk]->GetCurve();
          SmSolutionArray sSolutions;

          // Use this little half interval trick to eliminate the possibility of
          // finding closed intersections.  It may be a little dangerous.
          NER( pRailCurve );  //cbi NO:  SER(!pRailCurve);  // RCLxx
          SmExtent1d sHalfInterval = pRailCurve->GetNaturalInterval();
          (void)sHalfInterval.SetMinMax(sHalfInterval.Evaluate(0.5),sHalfInterval.GetMax()); // SolidCraft: from the middle to the end, it cannot fail

          // search for self intersections
          SER(pRailCurve->GlobalCurveSolve(sHalfInterval,
                                           *apRailCurves[kkk],
                                           apRailCurves[kkk]->GetNaturalInterval(),
                                           SM_SO_MINIMIZE,
                                           SM_EFF_ZERO_SQRT,
                                           NULL,
                                           NULL,
                                           SM_SR_SINGLE,
                                           sSolutions));

          // for every solution
          if (   sSolutions.GetSize() > 0
              && sSolutions[0].m_vStart.m_dSolutionValue < dMin)
            {
              dT1 = sSolutions[0].m_vStart[0];
              dT2 = sSolutions[0].m_vStart[1];

#ifdef SM_DEBUG_CODE
SmBoolean bDebugMe4 = FALSE;
              if (bDebugMe4) 
               {
                  SmPoint3d sPnt1, sPnt2;
                  pRailCurve->EvaluatePoint(dT1, sPnt1);
                  apRailCurves[kkk]->EvaluatePoint(dT2, sPnt2);
                  sPnt1.Draw();
                  sPnt2.Draw();
                  sm_GraphicsLoop();
                }
#endif // SM_DEBUG_CODE

              dMin = sSolutions[0].m_vStart.m_dSolutionValue;
              lMin = kkk;
            } // end found a solution check
        } // end iter both rail curves

      // no self intersections were found
      if (dMin == SM_BIG_DOUBLE)
        {
          // One of the rail should intersect
          SER(SM_ERR);
        }
    } // end no intersection found check

  // Evaluate centerLine and MinRailCurve at self-intersection point
  SmVector3d sPV1[2];
  SER(pOrigCenterLine->Evaluate(dT1,1,TRUE,sPV1));

  SmPoint3d sPnt;
  SmCurve * pOrigRailCurve = pOrigRail[lMin]->GetCurve();
  SER(pOrigRailCurve->EvaluatePoint(dT1,sPnt));
  double dRadius = sPV1[0].DistanceBetween(sPnt);

  // pick a step-back value
  double dStepBackDistance = m_dStepBackDistFactor * dRadius;
  if(   dMin < SM_BIG_DOUBLE
     && dMin > 3.0*dStepBackDistance)
    {
      return SM_ERR;
    }  // end step-back dist check

  double dTanLeng = sPV1[1].Length();
  if (SM_IS_ZERO(dTanLeng)) SER(SM_ERR);
  double dT1Step  = dT1 - dStepBackDistance/dTanLeng;

  // If the step goes too far then step back 1/2 of distance to start
  // of curve.
  if(dT1Step <= sTrimIvl1.GetMin()) 
    { dT1Step = (dT1 + sTrimIvl1.GetMin()) / 2.0; }

  // Try to snap it to a valid knot of the fillet to remain accurate
  SmTArray<double> sOrigKnots;
  SER(pOrigCenterLine->GetKnots(sOrigKnots));
  for (ULONG m=sOrigKnots.GetSize()-1; m>0; m--)
    {
      double dParam = sOrigKnots[m];
      if (dT1Step > dParam) 
        {
          dT1Step = dParam;
          break;
        }
    }

  // evaluate rail curve at max self intersection point
  SmVector3d sPV2[2];
  SER(apRailCurves[lMin]->Evaluate(dT2,1,TRUE,sPV2));

  // pick a step back distance
  dTanLeng = sPV2[1].Length();
  if (SM_IS_ZERO(dTanLeng)) SER(SM_ERR);
  double dT2Step = dT2 + dStepBackDistance/dTanLeng;
  if(dT2Step >= sTrimIvl2.GetMax())
    { dT2Step = (dT2 + sTrimIvl2.GetMax()) / 2.0; }

  // Do snapping for the max self-intersection param value
  SmTArray<double> sKnots;
  SER(pCenterLine->GetKnots(sKnots));
  for(ULONG n=0; n<sKnots.GetSize(); n++) 
    {
      double dParam = sKnots[n];
      if (dT2Step < dParam) 
        {
          dT2Step = dParam;
          break;
        }
    }

  // Now Trim Back all curves with new intervals
  if (dT1Step < sTrimIvl1.GetMax()) { (void)sTrimIvl1.SetMinMax(sTrimIvl1.GetMin(),dT1Step); } // SolidCraft: ordered by the test, it cannot fail
  if (dT2Step > sTrimIvl2.GetMin()) { (void)sTrimIvl2.SetMinMax(dT2Step,sTrimIvl2.GetMax()); } // SolidCraft: ordered by the test, it cannot fail

  // Trim both CenterLine and pOrigCenterLIne to new Interval1 and Interval2
  // GWC:??? pCenterLine and pOrigCenterLine need to be different
  SER(pCenterLine->Trim(sTrimIvl2));      // may snap sIvl by tol to existing knots
  SER(pOrigCenterLine->Trim(sTrimIvl1));  // may snap sIvl by tol to existing knots

  SmFilletVertex * pBlendStartV[2];
  SmFilletVertex * pBlendEndV[2];

  // for all three segments of the FilletGeom:[Beg Mid End]
  for (ULONG ii=0; ii<2; ii++)
    {
      // get thisSegments railEdge's UVCurve and 3DCurve
      SmFilletEdgeuse * pPrimEU  = (SmFilletEdgeuse*)pOrigRail[ii]->GetPrimaryEdgeuse();
      SmFilletEdgeuse * pMateEU  = (SmFilletEdgeuse*)pPrimEU->GetMate();
      SmBSplineCurve  * pUVCurve = pMateEU->GetUVTrimCurve(); NER(pUVCurve);
      SmBSplineCurve  * p3DCurve = SM_CAST_PTR(SmBSplineCurve, pOrigRail[ii]->GetCurve()); NER(p3DCurve);

      // trim Segment's 3DCurve and UVCurve to sTrimIvl1 - // gwc: this does not look right - review
      pOrigRail[ii]->SetInterval(sTrimIvl1);               // we're iterating on Segments, but not on intervals
      SER(pUVCurve->Trim(sTrimIvl1));  // may snap sIvl by tol to existing knots
      SER(p3DCurve->Trim(sTrimIvl1));  // may snap sIvl by tol to existing knots

      // trim the edgeuse to sTrimIvl2
      SER(apUVCurves[ii]->Trim(sTrimIvl2));    // may snap sIvl by tol to existing knots   // gwc: this does not look right - review
      SER(apRailCurves[ii]->Trim(sTrimIvl2));  // may snap sIvl by tol to existing knots   // we're iterating on Segments, but not on intervals

      // get the blend start UV and XYZ positions
      SmPoint3d        sPnt1, sUVPnt1;
      SmFilletVertex * pV = (SmFilletVertex*)pOrigRail[ii]->GetVertex();
      pBlendStartV[ii]    = (SmFilletVertex*)pOrigRail[ii]->GetOtherVertex(pV);

      SER(p3DCurve->EvaluatePoint(sTrimIvl1.GetMax(),sPnt1));
      SER(pUVCurve->EvaluatePoint(sTrimIvl1.GetMax(),sUVPnt1));

      pBlendStartV[ii]->SetPoint(sPnt1);
      SmPoint2d sUV1(sUVPnt1.x,sUVPnt1.y);
      pBlendStartV[ii]->SetOriginalUV(sUV1);

      // Build pNewFilletGeom railEdge geometry
      SmFilletEdge    * pNewRail   = pNewFilletGeom->GetRail(ii);
      SmFilletEdgeuse * pNewPrimEU = (SmFilletEdgeuse*)pNewRail->GetPrimaryEdgeuse();
      SmFilletEdgeuse * pNewMateEU = (SmFilletEdgeuse*)pNewPrimEU->GetMate();

      // Set thisSegment's Curve
      pNewRail->SetCurve(apRailCurves[ii], TRUE ) ; // TRUE = delete preExisting Edge->Curve - expect none here
                                                    // side effect: delete current pSurviveEdge->UVTrimCurve
      apRailCurves[ii]->SetOwner(pNewRail);

      // Set thisSegment's UVtrimCurve
      pNewMateEU->SetUVCurve(apUVCurves[ii]);
      apUVCurves[ii]->SetOwner(pNewMateEU);

      pNewRail->SetInterval(sTrimIvl2);

      // Set thisSegment's Rail->Vertex->Positions
      SmPoint3d sPnt2, sUVPnt2;
      pBlendEndV[ii] = (SmFilletVertex*)pNewRail->GetVertex();
      SER(apRailCurves[ii]->EvaluatePoint(sTrimIvl2.GetMin(),sPnt2));
      SER(apUVCurves[ii]->EvaluatePoint(sTrimIvl2.GetMin(),sUVPnt2));

      pBlendEndV[ii]->SetPoint(sPnt2);
      SmPoint2d sUV2(sUVPnt2.x,sUVPnt2.y);
      pBlendEndV[ii]->SetOriginalUV(sUV2);

#ifdef SM_DEBUG_CODE
SmBoolean bDebugMe = FALSE;
      if (bDebugMe) 
        {
          smgfx_SetLook(1,2, 1,0,0); p3DCurve->DrawWDeriv(sTrimIvl1,0); sm_GraphicsLoop();
          smgfx_SetLook(1,2, 0,0,1); apRailCurves[ii]->DrawWDeriv(sTrimIvl2,0); sm_GraphicsLoop();
          smgfx_SetLook(1,2, 0,0,1); pBlendStartV[ii]->Draw(); sm_GraphicsLoop();
          smgfx_SetLook(1,2, 0,0,1); pBlendEndV[ii]->Draw(); sm_GraphicsLoop();
          sm_GraphicsLoop();
        }
#endif // SM_DEBUG_CODE

      // Set blending rail data
      SmFilletEdge *pBlendRail = pBlendFilletGeom->GetRail(ii);
      pBlendRail->SetFilletEdgeType(SM_FE_BLENDING_RAIL);
      pBlendRail->SetOriginalFace(pOrigRail[ii]->GetOriginalFace());

    } // end iter all three segments of the FilletGeom:[Beg Mid End]

  // store the centerLine and FilletSurface for the 1st part of the fillet
  pNewFilletGeom->SetFilletSurface(pNewFilletSurface);
  pNewFilletGeom->SetCenterLineCurve(pCenterLine);

  // Calculate all 4 fillet vertices for pBlendFilletGeom
  SmPoint3d sPnt1, sPnt2;
  SmBoolean bFoundSolution;

  // First, calculate pBlendStartV[0] & pBlendStartV[1]
  SmFilletVertexuse * pVU1       = pBlendStartV[0]->GetVUAtRailEnd(pOrigFilletGeom);
  SmTsectPnt        & rTsectPnt1 = pVU1->GetTsectPnt();
  SmVector2d          sUV1       = pBlendStartV[0]->GetOriginalUV();
  SmVector2d          sUV2       = pBlendStartV[1]->GetOriginalUV();

  // Solve the point on the cross-section plane
  SmVector3d sPV[2];
  SER(pOrigCenterLine->Evaluate(sTrimIvl1.GetMax(),1,TRUE,sPV));

  // find filletPoint on given plane
  //  satisfying geometry requirements implemented in derived SmFilletSolver class
  SER(pFS->PointOnPlaneSolve(sPV[0],sPV[1],sDomain1,sDomain2,
      sUV1,sUV2,bFoundSolution,rTsectPnt1));
  if (!bFoundSolution) SER(SM_ERR);

  // locals for the railPoints at the cutBack positions
  sUV1 = rTsectPnt1.UVPos(0);
  sUV2 = rTsectPnt1.UVPos(1);
  SER(pSurf1->EvaluatePoint(sUV1,sPnt1));

  pBlendStartV[0]->SetOriginalUV(sUV1);
  pBlendStartV[0]->SetPoint(sPnt1);
  pBlendStartV[0]->SetStatus(SM_FIL_PROCESSED);

  SER(pSurf2->EvaluatePoint(sUV2,sPnt2));

  pBlendStartV[1]->SetOriginalUV(sUV2);
  pBlendStartV[1]->SetPoint(sPnt2);
  pBlendStartV[1]->SetStatus(SM_FIL_PROCESSED);

  SmFilletVertexuse * pMateVU     = pBlendStartV[1]->GetVUAtRailEnd(pOrigFilletGeom);
  SmTsectPnt        & rTsectPnt11 = pMateVU->GetTsectPnt();
  rTsectPnt11 = rTsectPnt1;

  // Then, calculate pBlendEndV[0] & pBlendEndV[1]
  SmFilletVertexuse * pVU2       = pBlendEndV[0]->GetVUAtRailEnd(pNewFilletGeom);
  SmTsectPnt        & rTsectPnt2 = pVU2->GetTsectPnt();
  sUV1 = pBlendEndV[0]->GetOriginalUV();
  sUV2 = pBlendEndV[1]->GetOriginalUV();

  // Solve the point on the cross-section plane
  SER(pCenterLine->Evaluate(sTrimIvl2.GetMin(),1,TRUE,sPV));

  // find filletPoint on given plane
  //  satisfying geometry requirements implemented in derived SmFilletSolver class
  SER(pFS->PointOnPlaneSolve(sPV[0],
                             sPV[1],
                             sDomain1,
                             sDomain2,
                             sUV1,
                             sUV2,
                             bFoundSolution,
                             rTsectPnt2));
  if (!bFoundSolution) SER(SM_ERR);

  sUV1 = rTsectPnt2.UVPos(0);
  sUV2 = rTsectPnt2.UVPos(1);

  SER(pSurf1->EvaluatePoint(sUV1,sPnt1));

  pBlendEndV[0]->SetOriginalUV(sUV1);
  pBlendEndV[0]->SetPoint(sPnt1);
  pBlendEndV[0]->SetStatus(SM_FIL_PROCESSED);

  SER(pSurf2->EvaluatePoint(sUV2,sPnt2));

  pBlendEndV[1]->SetOriginalUV(sUV2);
  pBlendEndV[1]->SetPoint(sPnt2);
  pBlendEndV[1]->SetStatus(SM_FIL_PROCESSED);

  pMateVU = pBlendEndV[1]->GetVUAtRailEnd(pNewFilletGeom);

  SmTsectPnt & rTsectPnt22 = pMateVU->GetTsectPnt();
  rTsectPnt22 = rTsectPnt2;

#ifdef SM_DEBUG_CODE
SmBoolean bDebugMe1 = FALSE;
  if (bDebugMe1) 
    {
      smgfx_SetColor(1,0,0); pBlendStartV[0]->Draw(); sm_GraphicsLoop();
      smgfx_SetColor(1,0,0); pBlendStartV[1]->Draw(); sm_GraphicsLoop();
      smgfx_SetColor(0,0,1); pBlendEndV[0]->Draw(); sm_GraphicsLoop();
      smgfx_SetColor(0,0,1); pBlendEndV[1]->Draw(); sm_GraphicsLoop();
      sm_GraphicsLoop();
    }
#endif // SM_DEBUG_CODE

  // all done
  return SM_SUCCESS;

} // end SmMakeSurfaceBlendSIH::CreateBlends


