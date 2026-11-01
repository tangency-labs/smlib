// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

/**********************************************************************//**
* FILE NAME --- SmBlendCornerEnd.cpp
* PURPOSE: SolidCraft change. The End event of the blend corner: a run of
*   one blend ending at a vertex whose side faces there are not blended.
*   Here, the end of a rail at a corner of mixed convexity: where the
*   rail meets the end face's surface continued past the side edge it
*   would otherwise meet (the side edge of the blend's own convexity, which
*   the other side edge, of the opposite convexity, turns away from).
**********************************************************************/

#include "StdAfx.h"
#include <SmFilletGeom.h>
#include <SmFilletCorner.h>
#include <SmFilletSolver.h>
#include <SmSolutionArray.h>
#include <SmSurface.h>
#include <SmCurve.h>
#include <SmVertex.h>
#include <SmBlendCorner.h>
#include <SmFilletExecutive.h>
#include <SmEdgeuse.h>
#include <SmEdge.h>
#include <SmFace.h>
#include <SmPlane.h>
#include <SmGeomUtility.h>
#ifdef SM_DEBUG_CODE
  #include <SmGraphicsOutput.h>
#endif

/*******************************************************************//**
PURPOSE: The rail's end at a mixed-convexity End
    (SM_FV_RAIL_X_EXTENDED_EDGEUSE of an end rib, SM_FE_FILLET_X_SIDE_FACE):
    the rail curve's crossing with the end face's surface continued.

NOTES: The rail lies on the side face F of the blend whose side edge E at
   the vertex has the blend's convexity, and the other side edge has the
   opposite one: the blend runs on past the vertex, and E continued past it
   (on F and on the end face's surface, both continued) is where the rail
   ends. That point is the rail's crossing with the end face's surface: a
   curve-surface event solved on the rail as traced (run on past the
   vertex: the corner is no blending corner), at the crossing nearest the
   vertex. The end face's surface is the corner's continuation of it
   (SmFilletCorner::GetExtendedSurface: an analytic surface continues as
   itself), else the face's own. The end rib, from this end to the other
   rail's end on the other side edge, is then the blend surface's
   intersection with that surface between the two ends
   (SmFilletEdge::CalcFilletIntSideFace).
   The rails exist once the blend is traced: before that (the corner's
   first pass) the vertex is left unprocessed.
***********************************************************************/
SmStatus SmFilletVertex::CalcRailIntEndSurface()
{
  NER( m_cpCorner );

  // the end rib and this vertex's rail
  SmFilletEdge * pRib  = NULL;
  SmFilletEdge * pRail = NULL;
  SmTArray<SmEdge*> sEdges;
  GetEdges( sEdges );
  for (ULONG ii=0; ii<sEdges.GetSize(); ii++)
    {
      SmFilletEdge * pFE = (SmFilletEdge*)sEdges[ii];
      if ( pFE->GetFilletEdgeType() == SM_FE_FILLET_X_SIDE_FACE && pRib == NULL )
        { pRib = pFE; }
      else if ( pFE->GetFilletEdgeType() == SM_FE_RAIL && pRail == NULL )
        { pRail = pFE; }
    }
  NER( pRib );
  NER( pRail );
  SmFilletEdgeuse * pRibEU      = (SmFilletEdgeuse*)pRib->GetPrimaryEdgeuse();
  NER( pRibEU );
  SmFilletGeom    * pFilletGeom = pRibEU->GetFilletGeom();
  NER( pFilletGeom );
  SmCurve         * pRailCurve  = pRail->GetCurve();
  if ( pFilletGeom->GetFilletSurface() == NULL || pRailCurve == NULL )
    { return SM_SUCCESS; } // not traced yet: the corner's second pass finds it

  SmFace    * pEndFace = pRib->GetOriginalFace();
  NER( pEndFace );
  SmSurface * pEndSurface = m_cpCorner->GetExtendedSurface( pEndFace );
  if ( pEndSurface == NULL )
    { pEndSurface = pEndFace->GetSurface(); }
  NER( pEndSurface );

  SmFilletSolver * pFilSolver = pFilletGeom->GetFilletSolver();
  NER( pFilSolver );
  const double dTol = pFilSolver->GetThisApproxTol3d();

  SmSolutionArray sSolutions;
  SER( pEndSurface->GlobalCurveIntersect( pEndSurface->GetNaturalUVDomain(), *pRailCurve,
                                          pRailCurve->GetNaturalInterval(), dTol, sSolutions ));
  const SmPoint3d sVertPnt = m_cpCorner->GetFilletedVertex()->GetPoint();
  SmBoolean bFound = FALSE;
  SmPoint3d sBest;
  double    dBest = SM_BIG_DOUBLE;
  for (ULONG ii=0; ii<sSolutions.GetSize(); ii++)
    {
      SmPoint3d sPnt;
      SER( pRailCurve->EvaluatePoint( sSolutions[ii].m_vStart.m_adParameters[0], sPnt ));
      double dDist = sPnt.DistanceBetween( sVertPnt );
      if ( dDist < dBest )
        {
          dBest  = dDist;
          sBest  = sPnt;
          bFound = TRUE;
        }
    }
  if ( !bFound )
    {
      SetStatus( SM_FIL_SURF_INT_FAILURE );
      SER_MSG( SM_ERR, _T("CalcRailIntEndSurface(): the rail does not meet the end face's surface") );
    }

  SetPoint( sBest );
  SetStatus( SM_FIL_PROCESSED );
  return SM_SUCCESS;

} // end SmFilletVertex::CalcRailIntEndSurface

/*******************************************************************//**
PURPOSE: Whether a round ending at the vertex runs out there: its side
    faces continue each other through the vertex (their normals agree)

NOTES:
***********************************************************************/
SmBoolean SmBlendCorner::EndRunsOut
  (const SmVertex * pCornerVert,
   SmEdgeuse * pFilEdgeuse)
{
    // Get adjacent edgeuses of filleted edgeuse
    SmEdgeuse * pSideEU1 = pFilEdgeuse->GetCWEdgeuse();
    SmEdgeuse * pSideEU2 = pFilEdgeuse->GetRadial()->GetCCWEdgeuse();
    if (pSideEU1 == pSideEU2) {
        // likely 2x1 case
        pSideEU2 = pSideEU1->GetMate();
    }

    SmVertex * pVert1 = pSideEU1->GetVertexuse()->GetVertex();
    SmExtent1d sIvl1 = pSideEU1->GetEdge()->GetInterval();
    double dT1 = sIvl1.GetMax();
    if (   (pVert1 == pCornerVert && pSideEU1->GetOrientation() == SM_OT_SAME)
        || (pVert1 != pCornerVert && pSideEU1->GetOrientation() == SM_OT_OPPOSITE)) 
      {
        dT1 = sIvl1.GetMin();
      }
    SmVertex * pVert2 = pSideEU2->GetVertexuse()->GetVertex();
    SmExtent1d sIvl2 = pSideEU2->GetEdge()->GetInterval();
    double dT2 = sIvl2.GetMax();
    if (   (pVert2 == pCornerVert && pSideEU2->GetOrientation() == SM_OT_SAME)
        || (pVert2 != pCornerVert && pSideEU2->GetOrientation() == SM_OT_OPPOSITE)) 
      {
        dT2 = sIvl2.GetMin();
      }
    SmVector3d sBinPnt1, sBinPnt2, sBinVec1, sBinVec2, sFaceuseNormal1, sFaceuseNormal2;
    SER(pSideEU1->EvaluateBinormal(dT1,FALSE,sBinPnt1,sBinVec1,NULL,&sFaceuseNormal1));
    SER(pSideEU2->EvaluateBinormal(dT2,FALSE,sBinPnt2,sBinVec2,NULL,&sFaceuseNormal2));
    // Check if faceuse normals are parallel to each other
#ifdef SM_DEBUG_CODE
    int iDebugLevel = 0; // Note: can't call DebugLevel() from a static method.
    if ( iDebugLevel > 0 ) {
        smgfx_SetLook(4,10, 1,0,0); pSideEU1->Draw(); sm_GraphicsLoop();
        smgfx_SetLook(4,10, 0,0,1); sFaceuseNormal1.Draw(&sBinPnt1); sm_GraphicsLoop();
        smgfx_SetLook(4,10, 0,1,0); pSideEU2->Draw(); sm_GraphicsLoop();
        sFaceuseNormal2.Draw(&sBinPnt2); sm_GraphicsLoop();
        sm_GraphicsLoop();
    }
#endif // SM_DEBUG_CODE
    // SolidCraft change: the side faces are tangent at the corner and continue each other (their
    // normals parallel, the edge smooth there): the fillet runs out to the vertex. Normals
    // anti-parallel are a cusp, the faces folded back on each other (a knife edge): the ball
    // touches them back from the vertex and the fillet ends there as at any other corner.
    SmBoolean bAntiParallel = FALSE;
    if (   sFaceuseNormal1.IsParallelTo(sFaceuseNormal2,SM_EFF_ZERO_SQRT,&bAntiParallel)
        && !bAntiParallel) {
        return TRUE;
    }
    return FALSE;

} // end SmBlendCorner::EndRunsOut

/*******************************************************************//**
PURPOSE: The End on its end face(s) (SM_BEND_ON_END_FACE): the FilletEdge
    boundary of the incoming fillet. When needed define
    endFaces and extra FilletEdges as needed.

NOTES:
 For all cases:
 -  new SmFilletVertex objects are constructed and classified to each rail
        m_pTargetBrep->edgeuse coming into this corner.
 -  new SmFilletVertices are added to this SmFilletCorner::m_vVertices array.

 - 1 new SmFilletEdge and 2 new SmFilletEdgeuses are constructed to
          connect appropriate new SmFilletVertex pairs.
 - New SmFilletEdges are added to this SmFilletCorner::m_vEdges array.

 Every such End needs an end Face bounded by a FilletEdge (connecting
 the two side edges at the rail/sideEdge intersection points)
 and the side edges themselves from the original vertex to the
 to the rail/sideEdge intersections.  In various cases
 different bits of the endFace (surface, edges, and vertices) need
 to be defined or are already defined.  Every case
 needs 1 filletEdge defined as the intersection of the
 filletSurface and the endFace.  This edge is always terminated
 by a pair of filletVertices located at rail/sideEdge intersections.

 In addition to this vertex/edge/vertex set, sometimes a face will need
 to be extended (or a new end-cap face defined) because the fillet
 surface would not intersect anything even when extended.  Sometimes,
 in addition to the extended end face, one (but never both) of the
 side edges must be extended in order to intersect the corresponding
 rail edge of the fillet surface.  The cases depend on the relative
 convexities of the fillet edge and two topologically adjacent side
 edges.  (Note, if there are more than three edges, the convexities
 of the remaining 'remote' edges don't matter.)  If all three edges
 have the same convexity, then no extending needs to be done.  In
 mixed-convexity cases, it depends on whether one or both of the side
 edges have convexities opposite to the filleted edge:

 Convexity of side edges relative to the filleted edge:
 ----------
  both same: nothing -- uniform convexity case.
  both opp : extend end face(s).
  same/opp : extend end face(s), and extend the side edge w/ same cvxty.

 (Note, 'side faces' are the two faces that contain the filleted edge,
 and 'end face(s)' are all other faces around this vertex (just one, in
 the normal 3x1 case).  End face(s) are what the fillet surface will
 intersect.)

 In the mixed-convexity case, the end face of the fillet will be either
 an extension of the original end surface, if the new fillet edges from
 the original vertex on the filleted edge to the rail/sideEdge intersections
 lie in the original end surface, otherwise a new end patch will be
 constructed that is G1 with end faces(s) at the two side edges.

 Note, regardless of convexities, there are cases where in which the
 expected intersections of fillet rail curves with side edges, or of
 the fillet surface with end surface(s), will not occur, due to
 geometric problems.  These will happen if the existing geometry
 (side edges or end faces) are short (generally, shorter than the
 fillet radius), or bend away from the vertex too quickly.
 These cases are not (yet) handled by this routine.


***********************************************************************/
SmStatus SmBlendCorner::MakeEndOnEndFace()
{
  // locals
  SmEdgeuse      * pFilEdgeuse = (*m_pSolverEUs)[0];
  SmFilletSolver * pFilSolver  = GetFilletExecutive()->GetFilletSolverOfEdgeuse(pFilEdgeuse);
  NER(pFilSolver);
  SmFilletGeom   * pFilletGeom = pFilSolver->GetFirstFilletGeom();
  ULONG            lTotalEUs   = m_pAllEUs->GetSize();
  ULONG            lFilEUIndex = SM_BIG_ULONG;

  // get the FIlEdgeuse index in the m_pAllEUs array
  ULONG ii;
  for ( ii=0; ii<lTotalEUs; ii++ )
    {
      if ( (*m_pAllEUs)[ii] == pFilEdgeuse )
        {
          lFilEUIndex = ii;
          break;
        }
    }
  SM_ASSERT(lFilEUIndex != SM_BIG_ULONG);

  // Get Adjacent edge->edgeuses of filleted edge
  // (the edges sharing a face with this filleted edge).

  SmEdgeuse * pEU[2];
  ULONG lIndexA = ((lFilEUIndex+lTotalEUs)-1)%lTotalEUs ; // 'behind' fil index
  ULONG lIndex1 = (lFilEUIndex+1)%lTotalEUs ;           // 'ahead of' fil index
  SM_ASSERT((lIndexA + 1)%lTotalEUs         == lFilEUIndex) ;
  SM_ASSERT(((lIndex1+lTotalEUs)-1)%lTotalEUs == lFilEUIndex) ;
  pEU[0] = (*m_pAllEUs)[lIndexA];
  pEU[1] = (*m_pAllEUs)[lIndex1];

  // check for 2x1 on non-closed-edge case with closed but not periodic side edge.
  SmBoolean bClosedSideEdge = (pEU[1] == pEU[0]) ;

int iDebugLevel = 0;

#ifdef SM_DEBUG_CODE
  // draw FilEdgeuse (red) and neighbor edgeuses (green, yellow)
  // connected to edges sharing a common face
  iDebugLevel = DebugLevel();
  if ( iDebugLevel > 0 )
    {
      smgfx_SetLook(3,5, 1,0,0); pFilEdgeuse->Draw(); sm_GraphicsLoop();
      smgfx_SetLook(3,5, 0,1,0); pEU[0]->Draw();      sm_GraphicsLoop(); // adjacent edgeuses to filleted edge
      smgfx_SetLook(3,5, 1,1,0); pEU[1]->Draw();      sm_GraphicsLoop();
      sm_GraphicsLoop();
    }
#endif // SM_DEBUG_CODE

  // Check the relative convexities of the filleted edge
  // and the two topologically adjacent side edges.  (See notes above.)

  SmBoolean bSideEdgeSameCvxty[2];
  bSideEdgeSameCvxty[0] = same_convexity_edges(pFilEdgeuse, pEU[0], m_cpVertex, iDebugLevel );
  bSideEdgeSameCvxty[1] = same_convexity_edges(pFilEdgeuse, pEU[1], m_cpVertex, iDebugLevel );

  // SolidCraft: a side edge whose faces meet tangentially at the vertex (a bar's flat side
  // running onto the cylinder its round end lies on) has no convexity: same_convexity_edges()
  // on it reads only the sign of the rounding in the geometry (IsConvexRadialSector's notes: a
  // G1 edge needs the step off), so the two ends of one arc read differently with the sketch
  // solver's residues. It is never of the opposite convexity: the side face runs on across it.
  // (Edges reached past such an edge round the vertex are read the same way below.)
  // Such an edge is a tangent seam T, and the rail on its side ends at a tangent end there:
  // where it crosses T, on the edge or on T continued past the vertex (bTangentSide).
  SmBoolean bTangentSide[2] = { FALSE, FALSE };
  {
    double dTanTolDeg = SM_RAD2DEG( GetTangencyTolRadians() );
    for ( ULONG kk=0; kk<2; kk++ )
      {
        if ( sm_CheckTangentialEdgeuseAtVertex( pEU[kk], m_cpVertex, dTanTolDeg ) )
          {
            bSideEdgeSameCvxty[kk] = TRUE;
            bTangentSide[kk]       = sm_IsTangentSeamAtVertex( pFilEdgeuse, pEU[kk], m_cpVertex, dTanTolDeg );
          }
      }
  }

  // A side edge whose two faces meet tangentially at the vertex (a flat top
  // running smoothly onto a round) is no edge of the corner: the side face the
  // rail runs on continues across it. When the edge reached past it round the
  // vertex has the filleted edge's convexity and the other side edge the
  // opposite one, the corner is the mixed-convexity one of that edge
  // (pCornerEU): its end face lies between it and the other side edge, and it
  // is the edge extended into the side face. pEU stays the neighbour on the
  // fillet's side face, which the rail ends against. lCornerEUs counts the
  // corner's edges.
  // So is the edge that continues the filleted edge tangentially across the
  // vertex with the opposite convexity, on the side face (past any tangent
  // edges, such as the side face's seam): a boss's base circle running on past
  // a wall as the rim of the boss's overhang. The rail on the side face runs
  // alongside it and never meets it; the round ends on the end face past it,
  // whose extension cuts the side face.
  SmEdgeuse * pCornerEU[2] = { pEU[0], pEU[1] };
  ULONG       lCornerEUs   = lTotalEUs;
  if ( ! bClosedSideEdge && lTotalEUs > 3 )
    {
      double dTanTolDeg = SM_RAD2DEG( GetTangencyTolRadians() );
      for ( ULONG kk=0; kk<2; kk++ )
        {
          ULONG lOther = 1 - kk;
          if (   bTangentSide[ kk ]
              || bSideEdgeSameCvxty[ lOther ]
              || sm_CheckTangentialEdgeuseAtVertex( pEU[ lOther ], m_cpVertex, dTanTolDeg ) )
            { continue; }

          ULONG lOtherIdx = ( kk == 0 ) ? lIndex1 : lIndexA;
          ULONG lIdx      = ( kk == 0 ) ? lIndexA : lIndex1;
          ULONG lSkipped  = 0;
          SmBoolean bPassedContinuation = FALSE;
          SmEdgeuse * pNextEU = pEU[ kk ];
          while ( lTotalEUs - lSkipped > 3 )
            {
              if ( ! sm_CheckTangentialEdgeuseAtVertex( pNextEU, m_cpVertex, dTanTolDeg ) )
                {
                  // the filleted edge's continuation, of the opposite convexity, once
                  if (   bPassedContinuation
                      || pNextEU->GetEdge() == pFilEdgeuse->GetEdge()
                      || ! sm_TestG1Continuity( m_cpVertex, pFilEdgeuse->GetEdge(), pNextEU->GetEdge(),
                                                dTanTolDeg, NULL, iDebugLevel )
                      || same_convexity_edges_round_vertex( pFilEdgeuse, pNextEU->GetEdge() ) )
                    { break; }
                  bPassedContinuation = TRUE;
                }
              lIdx = ( kk == 0 ) ? ( (lIdx + lTotalEUs) - 1 ) % lTotalEUs : ( lIdx + 1 ) % lTotalEUs;
              if ( lIdx == lFilEUIndex || lIdx == lOtherIdx )
                { lSkipped = 0; break; }
              pNextEU = (*m_pAllEUs)[ lIdx ];
              lSkipped++;
            }

          // The edge reached may itself be one whose faces meet tangentially
          // (a bar's flat side running onto its round end, below a wall the
          // bar's top meets): it has no convexity, as its faces' normals agree
          // along it, and IsConvexRadialSector() on it reads only the sign of
          // the rounding in the geometry (its notes: a G1 edge needs the step
          // off). It is not of the opposite convexity, so it is the corner's
          // edge: the end face beyond it runs on tangentially into the face
          // the rail runs onto.
          if (   lSkipped > 0
              && (   sm_CheckTangentialEdgeuseAtVertex( pNextEU, m_cpVertex, dTanTolDeg )
                  || same_convexity_edges_round_vertex( pFilEdgeuse, pNextEU->GetEdge() ) ) )
            {
              pCornerEU[ kk ]           = pNextEU;
              bSideEdgeSameCvxty[ kk ]  = TRUE;
              bTangentSide[ kk ]        = FALSE;  // the end face lies past the tangent edge
              lCornerEUs               -= lSkipped;
              break;
            }
        }
    }

  SmBoolean bMixedCvxty = ! bSideEdgeSameCvxty[0] || ! bSideEdgeSameCvxty[1];


  // for both rails, create filletvertices for rail/filletEdge intersections
  SmFilletVertex * pFVArray[2];

  // when working with mixed convexity:
  SmEdge * pExtendEdge = NULL; // sideEdge to be extended 
  SmFace * pEndFace    = NULL; // Face to extend when endFace->surface can
                               //   be defined by a surface extension.
                               //   In the 3x1 case, this is a unique face,
                               //   the one between the two side EU's.

  for(ULONG lThisRailIdx=0; lThisRailIdx<2; lThisRailIdx++)
    {
      ULONG lOtherRailIdx = 1 - lThisRailIdx;

      // construct FilletVertex to set on rail, and add it to m_vVertices
      pFVArray[ lThisRailIdx ] = new (m_pPseudoBrep) SmFilletVertex( this );
      m_vVertices.Add( pFVArray[ lThisRailIdx ] );
      SmEdgeuse * pSideEU = NULL;

      // get m_pTargetBrep->edge->edgeUse that intersects this rail
      ULONG lSolverRailIndex;
      if ( bClosedSideEdge )
        {
          // Assign rail index first, then find side edgeuse
          lSolverRailIndex = lThisRailIdx;
          pSideEU = sm_FindClosedSideEdgeuse(
                          lSolverRailIndex, pFilSolver, pEU[ lThisRailIdx ], iDebugLevel );
        }
      else
        { // Non-closed side edge (i.e. regular cases)
          pSideEU = sm_FindSideEdgeuse( pFilEdgeuse, pEU[ lThisRailIdx ], iDebugLevel  );
          lSolverRailIndex = pFilSolver->FindIndexOfRailXSideEdgeuse( pSideEU );
        }

      // set filletVertex to classify to the m_pTargetBrep->sideEdgeuse
      pFVArray[ lThisRailIdx ]->SetPointClass( SM_PC_EDGEUSE, pSideEU );

      // Connect vertex to rail edge with connection sequence:
      //   RailFilletEdge->newFilletEdgeuse->newFilletVertexuse->FilletVertex
      SmFilletEdgeuse * pNewEdgeuse = NULL;

      SER( pFilletGeom->MakeRailEdgeuse( m_pPseudoBrep,
              lSolverRailIndex, pFVArray[ lThisRailIdx ], pNewEdgeuse ));

      // set FilletVertex type:
      //    FilletVertexPosition will be the intersection of the rail
      //    with the side edge.  Note, in the mixed-convexity case,
      //    if the convexity of the side edge is the same as that of
      //    the filleted edge (so the other side edge has opposite cvxty
      //    to both), then this side edge will have to be extended.

      SmFace * pFace = pEU[ lThisRailIdx ]->GetFace();

      // At a mixed-convexity corner, if this side-edge has the same convexity
      // as the filleted edge and the other side-edge has opposite convexity,
      // then we'll need to extend this side-edge ('around the corner' of the
      // other side-edge) so that it will intersect this rail.

      if ( bMixedCvxty  &&  bSideEdgeSameCvxty[ lThisRailIdx ] )
        {
          // At a tangent side edge the rail's end is the tangent end, where the rail crosses
          // the seam T: an event of the fillet's trace, found with the other rail ends by the
          // solver's rail/edge crossing (SmFilletSolver::RailEdgeuseIntersect, which follows T
          // on as itself past the vertex), and the end rib starts there. It is never found
          // afterwards on the fillet surface: there the fillet surface and the end face touch,
          // and the end rib leaves the rail tangentially or at the node of their intersection.
          pFVArray[ lThisRailIdx ]->SetFilletVertexType(  bTangentSide[ lThisRailIdx ]
                                                        ? SM_FV_RAIL_X_EDGEUSE
                                                        : SM_FV_RAIL_X_EXTENDED_EDGEUSE );
          pFVArray[ lThisRailIdx ]->SetTangentEnd( bTangentSide[ lThisRailIdx ] );

          // Use the next edgeuse around the vertex, past the current one
          //SmEdgeuse * pExtendEdgeuse =
          //    ( lThisRailIdx==1 )
          //        ? m_pAllEUs->GetAt( ((lFilEUIndex+lTotalEUs)-2)%lTotalEUs)
          //        : m_pAllEUs->GetAt( (lFilEUIndex+2)%lTotalEUs) ;
          //
          // No, extend the edge with the same convexity, because it shares
          // a face with the fillet rail-to-be which it will intersect.
          // [bd, 3 Feb 06]
          // (Past a tangent side edge, the corner's side edge: see pCornerEU.)
          SmEdgeuse * pExtendEdgeuse = pCornerEU[ lThisRailIdx ];


          // Find original EndFace to be extended:
          //   the EndFace between the extended edge and the other side edge.
          // At a tangent end, the face across the seam T from the rail's face, whose surface
          // the end rib lies on (faces on that surface split by a tangent edge between T and
          // the other side edge are one end face).
          if ( bTangentSide[ lThisRailIdx ] )
            {
              pEndFace = pExtendEdgeuse->GetFace();
              if (   pEndFace == pFilEdgeuse->GetFace()
                  || pEndFace == pFilEdgeuse->GetRadial()->GetFace() )
                { pEndFace = pExtendEdgeuse->GetRadial()->GetFace(); }
            }
          else
            {
              pEndFace = sm_GetFaceByTwoEdgeuses(
                              pExtendEdgeuse, pCornerEU[ lOtherRailIdx ], NULL );
            }
          pExtendEdge = pExtendEdgeuse->GetEdge();

          // Point Edgeuse to the Fillet (side) Face containing the cliff edge.
          // (i.e. the one that's *not* getting extended.)
          // (i.e. the one into which the edge will be extended.)
          if (   pFace != pFilEdgeuse->GetFace()
              && pFace != pFilEdgeuse->GetRadial()->GetFace())
            {
              pEU[ lThisRailIdx ] = pEU[ lThisRailIdx ]->GetRadial();
            }

#ifdef SM_DEBUG_CODE
          if ( DebugLevel() > 0 ) 
            {
              smgfx_SetLook(1,2, 0,1,0); pExtendEdgeuse->Draw(); sm_GraphicsLoop();
              smgfx_SetLook(1,2, 1,0,0); pEU[lThisRailIdx]->Draw(); sm_GraphicsLoop();
              sm_GraphicsLoop();
            }
#endif // SM_DEBUG_CODE

        } // end extend-side-edge branch
      else  // Regular case
        {
          pFVArray[ lThisRailIdx ]->SetFilletVertexType( SM_FV_RAIL_X_EDGEUSE );

          // Point Edgeuse to an EndFace (not a fillet side face).
          if (   pFace == pFilEdgeuse->GetFace()
              || pFace == pFilEdgeuse->GetRadial()->GetFace())
            {
              pEU[ lThisRailIdx ] = pEU[ lThisRailIdx ]->GetRadial();
            }
        }
    } // end iter both rail vertices


  // What's left to do: create the topology for the end rib curve
  // (intersection of fillet surface with the end face),
  // and, if at a mixed-convexity corner, create a capping surface
  // (or extend the end surface).
  // But, if the end rib is not a simple one-piece thing (i.e.,
  // more than one end face), don't do it here.
  // That would be the case if lTotalEUs > 3.

  if ( !bMixedCvxty && lCornerEUs > 3 )
    {
      // Simple convex corner, OR more than one end face,
      // just let global merge do the work.

      m_pExecutive->SetDoGlobalMerge( TRUE );
      return SM_SUCCESS;
    }


  // Ok, now we're either capping (mixed convexity),
  // or the simple case of <=3 total edges at this vertex.

  // compute how much to extend the endSurface to make sure it covers the 'hole'
  double dExtDist =   bMixedCvxty 
                    ? sm_FindSurfaceExtensionDistAtCorner(pFilEdgeuse, 
                                                          pFilSolver, 
                                                          m_cpVertex )
                    : 0.0;

  // Ok, get to work on the topology of the end of the fillet.
  // Allocate the FilletEdge and its edgeuses connecting the two
  // rail/filletEdge intersection vertices.  This edge will be
  // the intersection of the fillet face with the end face.
  //  - Add new FilletEdge to this SmFilletCorner::m_vEdges list.
  //  - The shape of the new FilletEdge is not yet set - just the
  //    topological connection.

  SmFilletEdge * pFilletEndRibEdge = NULL;
  SER( SmFilletEdge::MakeFilletEdge(m_pPseudoBrep,     // in : target Brep to receive new topology objects
                                    m_vVertices[0],    // in : start of new FilletEdge
                                    m_vVertices[1],    // in : end   of new FilletEdge
                                    pFilletEndRibEdge, // out: newly allocated FilletEdge - no Curve or UVTrimCurve data
                                    this ));           // in : NewFilletEdge's filletCorner, gets stored in rpNewFilletEdge->m_cpCorner
                                                       //      NULL to ignore, default:[NULL]
  pFilletEndRibEdge->SetFilletEdgeType( SM_FE_FILLET_X_SIDE_FACE );
  SmFilletEdgeuse * pPrimEU = (SmFilletEdgeuse*)pFilletEndRibEdge->GetPrimaryEdgeuse();
  pPrimEU->SetFilletGeom( pFilletGeom );
  m_vEdges.Add( pFilletEndRibEdge );

  // In the mixed-cvxty case, pEndFace will be NULL when we arrive here
  // if both side edges have opposite convexity to the filleted edge.
  // In that case, the face connecting the two side edges will need to be
  // used (or extended) to generate an endFace.

  SmBoolean   bMultiNeighborEndFace = FALSE ;
  SmEdgeuse * pOtherSideEdgeuse     = NULL ;

  if ( pEndFace == NULL )
    {
      pEndFace = sm_GetFaceByTwoEdgeuses(pEU[0],pEU[1],NULL);

      // when sm_GetFaceByTwoEdgeuses() fails to find a result
      // there are more than 3 edges in this corner and the neighbor
      // edges are not connected by a single face but rather a sequence
      // of faces.

      if ( pEndFace == NULL )
        {
          // Also note that, from the previous bailout (where we "just let
          // global merge do the work"), we must be in the mixed-convexity
          // case here.
          SM_ASSERT( bMixedCvxty );

          bMultiNeighborEndFace = TRUE ;

          // assume that the endCap can be created by extending the
          // face->Surface connected to the 1st cliff side edge.
          // gwc:TODO in the future, new Surfaces can be defined which are
          //   tangent continuous with both side edges.
          if ( lCornerEUs > 3 )
            {
              MSG(_T("Building a FilletEndFace between two different faces - may be a G0 endPatch"));
              // Use the side edge that is a cliff edge.
              // If the two side edges have opposite convexities to each other,
              // then the cliff edge will be the one with the same convexity as
              // the filleted edge; otherwise (side edges both have opposite
              // convexity to the filleted edge -- remember, we must be in the
              // mixed-convexity case here) they're both cliff edges.
              if ( bSideEdgeSameCvxty[0] )
                {
                  pEndFace          = pEU[0]->GetFace();
                  pExtendEdge       = pEU[0]->GetEdge();
                  pOtherSideEdgeuse = pEU[1] ;
                }
              else if ( bSideEdgeSameCvxty[1] )
                {
                  pEndFace          = pEU[1]->GetFace();
                  pExtendEdge       = pEU[1]->GetEdge();
                  pOtherSideEdgeuse = pEU[0] ;
                }
              else  // both are cliff edges, just pick one.
                {
                  pEndFace          = pEU[0]->GetFace();
                  pExtendEdge       = pEU[0]->GetEdge();
                  pOtherSideEdgeuse = pEU[1] ;
                }
            }
          else
            {
              SM_DBG_WARN(_T("Found an unexpected multiNeighborEndFace with less than 4 vertex->edges")) ;
              SER(SM_ERR); // Invalid brep topology
            }
        }
    } // end pEndFace == NULL check

  pFilletEndRibEdge->SetOriginalFace( pEndFace );

  // compute how much to extend the sideSurface to make sure it covers the 'hole'
  SmBoolean bDoExtend = ( pExtendEdge != NULL );
  if ( dExtDist > SM_EFF_ZERO )
    {
      if (   bDoExtend == FALSE
          && (   pFilSolver->GetSolverType() == SM_FS_CONST_RADIUS
              || pFilSolver->GetSolverType() == SM_FS_CONST_RADIUS_ASSISTED))
        {
          dExtDist *= 0.1;   // use only a small extension
          bDoExtend = TRUE;  // for constant radius only
        }
    }

  // When a surface extension is needed to define the endFace shape,
  // build and save it.

  SmSurface *pNewExtendSurface = NULL ;
  if ( bDoExtend )
    {
      SmSurface * pSurface = pEndFace->GetSurface();
      NER(pSurface);
      SmSurface * pNewSurface = NULL;

      // create an extended surface large enough to fill in the endCap requirements
      // At a tangent end the end face is continued as itself: the end rib's shape there
      // depends on how its curvature differs from the rail's face across the seam, which a
      // linear (G1) extension would not keep.
      SmContinuityType eExtension = (   bMixedCvxty
                                     && ( bTangentSide[0] || bTangentSide[1] ))
                                    ? SM_CT_CINFINITY
                                    : SM_CT_G1;
      SER( pSurface->CreateExtendedSurface(m_crContext, dExtDist, eExtension, pNewSurface ));
      SmObjDelete sCleanSurf( pNewSurface );

      // Copy pNewSurface, when possible as an analytic surface
      SER( pNewSurface->CopyAndAddAnalytics( m_crContext, pNewExtendSurface ));

      // register the OrigFace/ExtendedSurface pair
      // in the SmFilletExecutive::m_vExtendedSurfacesMap list
      AddExtendedOriginalSurface( pEndFace, pNewExtendSurface );

      // If AddExtendedOriginalSurface() deleted pNewExtendSurface, then
      // pEndFace already had an extended surface, and we don't have to
      // deal with it anymore.  But if pNewExtendSurface stayed in the
      // map, we will have to do further checks on it.  So check whether
      // it changed:
      SmSurface *pTmpSrfPtr = GetExtendedSurface( pEndFace );
      if ( pTmpSrfPtr != pNewExtendSurface )
        {
          pNewExtendSurface = NULL;  // deleted in AddExtendedOriginalSurface()
        }
    }

  // When adding an endSurface to fill the 'hole',
  // set up the topology for the end cap.
  if ( bMixedCvxty )
    {
      // add a filletVertex at the origVertex location
      SmFilletVertex * pFV = new (m_pPseudoBrep) SmFilletVertex(this);
      pFV->SetPointClass(SM_PC_VERTEX, (SmObject*)m_cpVertex);
      pFV->SetFilletVertexType( SM_FV_ON_VERTEX );
      m_vVertices.Add( pFV );

      // add two filletEdges from origVertex to rail/edge intersections

      // Add 1st fillet edge from expected rail2/sideEdge2 intersection to origVertex
      SmFilletEdge * pNewFilletEdge2 = NULL;
      SER(SmFilletEdge::MakeFilletEdge(m_pPseudoBrep,   // in : target Brep to receive new topology objects
                                       m_vVertices[1],  // in : start of new FilletEdge
                                       m_vVertices[2],  // in : end   of new FilletEdge
                                       pNewFilletEdge2, // out: newly allocated FilletEdge - no Curve or UVTrimCurve data
                                       this));          // in : NewFilletEdge's filletCorner, gets stored in rpNewFilletEdge->m_cpCorner
                                                        //      NULL to ignore, default:[NULL]
      // set FilletEdgeType based on analysis done and saved in filletVertexType
      if (   m_vVertices[1]->GetFilletVertexType() == SM_FV_RAIL_X_EXTENDED_EDGEUSE
          || ( bMixedCvxty && bTangentSide[1] ) )
        {
          pNewFilletEdge2->SetFilletEdgeType( SM_FE_ON_EXTEND_EDGE );
          pNewFilletEdge2->SetOriginalEdge( pExtendEdge );
        }
      else
        {
          pNewFilletEdge2->SetFilletEdgeType( SM_FE_ON_EDGE );

          if (lCornerEUs <= 3)
            {
              // Register this edge and process it later
              // when removing topological vertices
              m_pExecutive->AddTopoEdges( pNewFilletEdge2 );
            }
          pNewFilletEdge2->SetOriginalEdge( pEU[1]->GetEdge() );
        }

      pNewFilletEdge2->SetOriginalFace( pEU[1]->GetFace() );
      m_vEdges.Add( pNewFilletEdge2 );

      // Add 2nd fillet edge from origVertex to other expected rail1/sideEdge intersection
      SmFilletEdge *pNewFilletEdge1 = NULL ;
      SER(SmFilletEdge::MakeFilletEdge(m_pPseudoBrep,   // in : target Brep to receive new topology objects
                                       m_vVertices[2],  // in : start of new FilletEdge
                                       m_vVertices[0],  // in : end   of new FilletEdge
                                       pNewFilletEdge1, // out: newly allocated FilletEdge - no Curve or UVTrimCurve data
                                       this));          // in : NewFilletEdge's filletCorner, gets stored in rpNewFilletEdge->m_cpCorner
                                                        //      NULL to ignore, default:[NULL]
      // set FilletEdgeType based on analysis done and saved in filletVertexType
      if (   m_vVertices[0]->GetFilletVertexType() == SM_FV_RAIL_X_EXTENDED_EDGEUSE
          || ( bMixedCvxty && bTangentSide[0] ) )
        {
          pNewFilletEdge1->SetFilletEdgeType( SM_FE_ON_EXTEND_EDGE );
          pNewFilletEdge1->SetOriginalEdge( pExtendEdge );
        }
      else
        {
          pNewFilletEdge1->SetFilletEdgeType( SM_FE_ON_EDGE );
          if ( lCornerEUs <= 3 )
            {
              m_pExecutive->AddTopoEdges( pNewFilletEdge1 );
            }
          pNewFilletEdge1->SetOriginalEdge( pEU[0]->GetEdge() );
        }
      pNewFilletEdge1->SetOriginalFace( pEU[0]->GetFace() );
      m_vEdges.Add( pNewFilletEdge1 );

     // all the corner topology is in place

     // When endFace connects to more than the filletSurface and 1 other Face,
     // check for gaps
     if ( bMultiNeighborEndFace )
       {
          SM_ASSERT(pNewFilletEdge1->GetFilletEdgeType() == SM_FE_ON_EDGE) ;
          SM_ASSERT(pNewFilletEdge2->GetFilletEdgeType() == SM_FE_ON_EDGE) ;

          // in this case we can set FilletVertex and SideFilletEdge shapes
          // without having the filletSurface - do so

          // set filletVertex->Positions
          SER( CalcCornerVertGeom()  );

          // set FilletSideEdge->Curves
          SER( pNewFilletEdge1->CalcCornerEdgeGeom()  );
          SER( pNewFilletEdge2->CalcCornerEdgeGeom()  );

          // When proposing an extended surface
          if ( pNewExtendSurface )
            {
              // get extendedSurface/OtherSideEdge coincidence to within tolerance
              SmBoolean  bFoundCoincidence = FALSE ;
              SmSolution sSolution ;
              SER( pNewExtendSurface->SimpleCoincidenceChecker(
                  pNewExtendSurface->GetNaturalUVDomain(),
                  *pOtherSideEdgeuse->GetEdge()->GetCurve(),
                  pOtherSideEdgeuse->GetEdge()->GetInterval(),
                  m_dThisApproxTol3d,
                  bFoundCoincidence,
                  sSolution
              ) );

              // when extendedSurface/OtherSideEdge are not coincident
              if( bFoundCoincidence == FALSE )
                {
                  // can't define endFace->Surface as a sideFace->ExtendedSurface
                  m_bExtendedSurfacePatch = FALSE ;

                }
            } // end pNewExtendSurface existence check
        } // end endFace connects to more than the filletSurface and 1 other Face check

    } // end bMixedCvxty check

  // all done
  return SM_SUCCESS;

} // end SmBlendCorner::MakeEndOnEndFace

/*******************************************************************//**
PURPOSE: The End at a lamina edge (SM_BEND_OPEN): its end section (SM_FE_FILLET_END)
    between the rails' crossings of the side edges; topology of an array of
    SmFilletVertex and an array of SmFilletEdge. Since we are dealing with
    OPEN end, the corner actually consists of only one fillet edge.

NOTES: For now, it will only handle 3x1 case.
***********************************************************************/
SmStatus SmBlendCorner::MakeEndOpen
  ()
{
    SmTArray<SmEdge*> sEdges;
    m_cpVertex->GetEdges(sEdges);
    SmEdgeuse * pFilEdgeuse = (*m_pSolverEUs)[0];

    SmEdgeuse * pEU[2]; // Adjacent edgeuses of filleted edgeuse
    pEU[0] = pFilEdgeuse->GetCWEdgeuse();
    pEU[1] = pFilEdgeuse->GetRadial()->GetCCWEdgeuse();
#ifdef SM_DEBUG_CODE
    if ( DebugLevel() > 0 ) {
        smgfx_SetLook(1,2, 1,0,0); pFilEdgeuse->Draw(); sm_GraphicsLoop();
        smgfx_SetLook(1,2, 1,0,0); pFilEdgeuse->GetRadial()->Draw(); sm_GraphicsLoop();
        smgfx_SetLook(1,2, 0,0,1); pEU[0]->Draw(); sm_GraphicsLoop();
        smgfx_SetLook(1,2, 0,0,1); pEU[1]->Draw(); sm_GraphicsLoop();
        sm_GraphicsLoop();
    }
#endif // SM_DEBUG_CODE

    SmFilletSolver * pFilSolver = m_pExecutive->GetFilletSolverOfEdgeuse(pFilEdgeuse);
    NER(pFilSolver);
    SmFilletGeom * pFilletGeom = pFilSolver->GetFirstFilletGeom();
    ULONG ii;
    for (ii=0; ii<2; ii++) {
        SmFilletVertex * pFV = new (m_pPseudoBrep) SmFilletVertex(this);
        pFV->SetPointClass(SM_PC_EDGEUSE, pEU[ii]);
        pFV->SetFilletVertexType(SM_FV_RAIL_X_EDGEUSE);
        m_vVertices.Add(pFV);
        ULONG lRailIndex = pFilSolver->FindIndexOfRailXSideEdgeuse(pEU[ii]);
        // Create rail-edgeuse using this vertex
        SmFilletEdgeuse * pNewEdgeuse = NULL;
        SER(pFilletGeom->MakeRailEdgeuse(m_pPseudoBrep,
            lRailIndex,pFV,pNewEdgeuse));
    }
    SmFilletEdge * pFilletEdge = NULL;
    SER(SmFilletEdge::MakeFilletEdge(m_pPseudoBrep,  // in : target Brep to receive new topology objects
                                     m_vVertices[0], // in : start of new FilletEdge
                                     m_vVertices[1], // in : end   of new FilletEdge
                                     pFilletEdge,    // out: newly allocated FilletEdge - no Curve or UVTrimCurve data
                                     this));         // in : NewFilletEdge's filletCorner, gets stored in rpNewFilletEdge->m_cpCorner
                                                     //      NULL to ignore, default:[NULL]
    pFilletEdge->SetFilletEdgeType(SM_FE_FILLET_END);
    SmFilletEdgeuse * pPrimEU = (SmFilletEdgeuse*)pFilletEdge->GetPrimaryEdgeuse();
    pPrimEU->SetFilletGeom(pFilletGeom);
    m_vEdges.Add(pFilletEdge);

    return SM_SUCCESS;

} // end SmBlendCorner::MakeEndOpen

/*******************************************************************//**
PURPOSE: The End where the round runs out (SM_BEND_RUNS_OUT): its side faces continue
    each other through the vertex, so both rails end at the vertex; topology of an array of
    SmFilletVertex and an array of SmFilletEdge. Since we are dealing with
    OPEN end, the corner actually consists of only one fillet edge.

NOTES: For now, it will only handle 3x1 case.
***********************************************************************/
SmStatus SmBlendCorner::MakeEndRunningOut
  ()
{
    SmTArray<SmEdge*> sEdges;
    m_cpVertex->GetEdges(sEdges);
    SmEdgeuse * pFilEdgeuse = (*m_pSolverEUs)[0];

    SmEdgeuse * pEU[2]; // Adjacent edgeuses of filleted edgeuse
    pEU[0] = pFilEdgeuse->GetCWEdgeuse();
    pEU[1] = pFilEdgeuse->GetRadial()->GetCCWEdgeuse();

    SmFilletSolver * pFilSolver = m_pExecutive->GetFilletSolverOfEdgeuse(pFilEdgeuse);
    SmFilletVertex * pFV = new (m_pPseudoBrep) SmFilletVertex(this);
    pFV->SetPointClass(SM_PC_VERTEX, (SmObject*)m_cpVertex);
    pFV->SetFilletVertexType(SM_FV_ON_VERTEX);
    SmPoint3d sCornerPnt = m_cpVertex->GetPoint();
    pFV->SetPoint(sCornerPnt);
    pFV->SetStatus(SM_FIL_PROCESSED);
    m_vVertices.Add(pFV);

#ifdef SM_DEBUG_CODE
    if ( DebugLevel() > 0 ) {
        smgfx_SetLook(1,2, 1,0,0);
        if (FALSE) {
            pFilEdgeuse->Draw(); sm_GraphicsLoop();
            pFilEdgeuse->GetRadial()->Draw(); sm_GraphicsLoop();
            sm_GraphicsLoop();
        }
        smgfx_SetLook(1,2, 0,1,0);
        if (FALSE) {
            pEU[0]->Draw(); sm_GraphicsLoop();
            pEU[1]->Draw(); sm_GraphicsLoop();
            sm_GraphicsLoop();
        }
        smgfx_SetLook(4,10, 0,1,1); pFV->Draw(); sm_GraphicsLoop();
        sm_GraphicsLoop();
    }
#endif // SM_DEBUG_CODE
    // Drop sCornerPnt onto 2 surfaces
    SmPoint2d sUV[2];
    ULONG ii;
    for (ii=0; ii<2; ii++) {
        const SmSurface * pBase = pFilSolver->GetSurface(ii)->GetBaseSurface();
        SmSolution sSData[16];
        SmSolutionArray sSolutions(16,sSData);
#ifdef SM_DEBUG_CODE
        if ( DebugLevel() > 0 ) {
            smgfx_SetLook(1,2, 0,0,0); pBase->DrawUV(5,5); sm_GraphicsLoop();
            sm_GraphicsLoop();
        }
#endif // SM_DEBUG_CODE
        SmSurfaceCache *pSC = smsurf_GetSurfaceCache(pBase); NER(pSC);
        SmCacheCheckOutIn sCheckIO(pSC);
        {
        // Turn off point testing so GlobalPointSolve() will keep all point solutions
        //   without classifying the solution point against the trim boundaries.
        SmTemporaryChangeValue<SmBoolean> sChange(pSC->m_bPointTestEnabled,FALSE);

        SER(pBase->GlobalPointSolve(pBase->GetNaturalUVDomain(),SM_SO_INTERSECT,
            sCornerPnt,m_cpVertex->GetTolerance()*10.0,NULL,SM_SR_SINGLE,sSolutions));
        }
        if (sSolutions.GetSize() != 1) SER(SM_ERR);
        sUV[ii].x = sSolutions[0].m_vStart[0];
        sUV[ii].y = sSolutions[0].m_vStart[1];
    }

    SmFilletGeom * pFilletGeom = pFilSolver->GetFirstFilletGeom();
    for (ii=0; ii<2; ii++) {
        ULONG lRailIndex = pFilSolver->FindIndexOfRailXSideEdgeuse(pEU[ii]);
        // Create rail-edgeuse using this vertex
        SmFilletEdgeuse * pNewEdgeuse = NULL;
        SER(pFilletGeom->MakeRailEdgeuse(m_pPseudoBrep,
            lRailIndex,pFV,pNewEdgeuse));
        SmFilletVertexuse * pVU = (SmFilletVertexuse*)pNewEdgeuse->GetVertexuse();
        SmTsectPnt & rTsectPnt = pVU->GetTsectPnt();
        rTsectPnt.UVPos(   lRailIndex ) = sUV[   lRailIndex ];
        rTsectPnt.UVPos( 1-lRailIndex ) = sUV[ 1-lRailIndex ];
        rTsectPnt.m_ePointType = SM_IP_TANGENT_POINT;
    }

    return SM_SUCCESS;

} // end SmBlendCorner::MakeEndRunningOut

/*******************************************************************//**
PURPOSE: Calculate the corner edges' geometry, see
    SmFilletCorner::CalcCornerEdgeGeom().

NOTES: At a corner of uniform convexity and at most three edges, the one
  corner edge is the end rib: where the fillet meets the end face between
  the two side edges, from one rail's end on its side edge to the other's
  (MakeCornerTopology). Each rail stays on its face, or rolls over a
  tangent edge, up to the fillet's ends, unless it overflows a cliff
  beside the edge (SmFilletSolver::RailsFitTheirFaces, LeavesPastAnEnd:
  a concave fillet is refused then, a convex one a cliff overflow). A rail can
  still miss its side edge (SM_FV_NO_INT_RAIL_EU), two ways:
  - The side edge ends within the fillet's width (a top round's end arc
    shortened by another round): the rail runs past the end face onto the
    faces beyond the side edge's other end, so the fillet ends on more
    than one face. As at a corner of more than three edges, there is no
    end rib then: the global merge cuts the fillet off where it meets the
    body.
  - The side edge lies on the fillet's other side surface too (a hole's
    rim split at a vertex): it continues the filleted edge between the
    fillet's own two surfaces, the rail runs alongside it and the fillet
    runs on past the vertex, crossing no end face. No end exists there,
    and the end rib fails as before.
  In a cliff overflow (SmFilletExecutive::m_bCliffOverflow) the fillet's
  faces continue past their cliffs, and so does the end face: a rail that
  misses its side edge there (either rail, or both), but for the second
  way, runs past the side edge's far end on its face's continuation. The
  end rib is left out then and the fillet uncut there; the executive cuts
  it off by the end face's surface continued (TrimAtCreases, CutFilletAtOpenEnd).
***********************************************************************/
SmStatus SmBlendCorner::CalcEndEdgeGeom
  ()
{
  SmFilletEdge * pEndRib = ( m_vEdges.GetSize() == 1 ) ? m_vEdges[0] : NULL;
  if (   pEndRib != NULL
      && pEndRib->GetFilletEdgeType() == SM_FE_FILLET_X_SIDE_FACE )
    {
      SmFilletVertex * pStartFV = (SmFilletVertex*)pEndRib->GetVertex();
      NER(pStartFV);
      SmFilletVertex * pEndFV   = (SmFilletVertex*)pEndRib->GetOtherVertex( pStartFV );
      NER(pEndFV);
      SmFilletVertex * pMissFV  =   ( pStartFV->GetStatus() == SM_FV_NO_INT_RAIL_EU && pEndFV->IsProcessed() )
                                  ? pStartFV
                                  : (   ( pEndFV->GetStatus() == SM_FV_NO_INT_RAIL_EU && pStartFV->IsProcessed() )
                                      ? pEndFV
                                      : NULL );
      // SolidCraft change: in a cliff overflow a rail can miss its side edge with the
      // intersection clamped to the side edge's far vertex, off the fillet: that is a miss too;
      // and both rails can miss.
      SmBoolean abMiss[2] = { pStartFV->GetStatus() == SM_FV_NO_INT_RAIL_EU, pEndFV->GetStatus() == SM_FV_NO_INT_RAIL_EU };
      if ( m_pExecutive->GetCliffOverflow() )
        {
          SmFilletEdgeuse * pRibEU = (SmFilletEdgeuse*)pEndRib->GetPrimaryEdgeuse();
          SmFilletGeom    * pRibFG = ( pRibEU != NULL ) ? pRibEU->GetFilletGeom() : NULL;
          SmSurface       * pFilletSurface = ( pRibFG != NULL ) ? pRibFG->GetFilletSurface() : NULL;
          NER(pFilletSurface);
          SmFilletVertex * apFV[2] = { pStartFV, pEndFV };
          for (ULONG kk=0; kk<2; kk++)
            {
              if ( abMiss[kk] || !apFV[kk]->IsProcessed() )
                { continue; }
              SmBoolean bDropped = FALSE, bMulti = FALSE;
              SmPoint2d sUV;
              double    dGap = SM_BIG_DOUBLE;
              SER( pFilletSurface->DropPoint( apFV[kk]->GetPoint(), pFilletSurface->GetNaturalUVDomain(), NULL, bDropped, sUV, dGap, bMulti ));
              SmBoolean bMisses = ( !bDropped || dGap > 10.0 * m_dThisApproxTol3d ) ? TRUE : FALSE;
              // The rail's end found on the side edge's line continued past the edge's end
              // (SmFilletSolver::RailEdgeuseIntersect continues a straight side edge): the rail
              // misses the side edge itself, however close past the end.
              if ( !bMisses && apFV[kk]->GetPointClass() == SM_PC_EDGEUSE )
                {
                  SmEdgeuse * pSideEU = (SmEdgeuse*)apFV[kk]->GetPointClassObject();
                  SmEdge    * pSideEdge = ( pSideEU != NULL ) ? pSideEU->GetEdge() : NULL;
                  NER(pSideEdge);
                  NER(pSideEdge->GetCurve());
                  SmBoolean bOnSide = FALSE;
                  double    dParam = 0.0, dDist = SM_BIG_DOUBLE;
                  SER( pSideEdge->GetCurve()->DropPoint( pSideEdge->GetInterval(), apFV[kk]->GetPoint(), NULL, SM_BIG_DOUBLE, NULL, bOnSide, dParam, dDist ));
                  bMisses = ( !bOnSide || dDist > m_dThisApproxTol3d ) ? TRUE : FALSE;
                }
              if ( bMisses )
                {
                  abMiss[kk] = TRUE;
                  apFV[kk]->SetStatus( SM_FV_NO_INT_RAIL_EU );  // its rail is not trimmed there
                }
            }
          if ( abMiss[0] || abMiss[1] )
            { pMissFV = abMiss[0] ? pStartFV : pEndFV; }
        }
      SmBoolean bBothMiss = m_pExecutive->GetCliffOverflow() && abMiss[0] && abMiss[1];
      if ( pMissFV != NULL && pMissFV->GetPointClass() == SM_PC_EDGEUSE )
        {
          // The side edge the rail missed, and the fillet's side face it does not bound.
          SmEdgeuse * pSideEU     = (SmEdgeuse*)pMissFV->GetPointClassObject();
          NER(pSideEU);
          SmEdgeuse * pFilEdgeuse = (*m_pSolverEUs)[0];
          SmFace    * pSideFaces[2] = { pSideEU->GetFace(), pSideEU->GetRadial()->GetFace() };
          SmFace    * pOtherFace    = pFilEdgeuse->GetFace();
          if ( pOtherFace == pSideFaces[0] || pOtherFace == pSideFaces[1] )
            { pOtherFace = pFilEdgeuse->GetRadial()->GetFace(); }
          NER(pOtherFace);
          SmSurface * pOtherSurface = pOtherFace->GetSurface();
          NER(pOtherSurface);
          SmEdge    * pSideEdge     = pSideEU->GetEdge();
          NER(pSideEdge->GetCurve());

          SmBoolean  bRunsOn = FALSE;
          SmSolution sSolution;
          SER( pOtherSurface->SimpleCoincidenceChecker( pOtherSurface->GetNaturalUVDomain(),
                                                        *pSideEdge->GetCurve(),
                                                        pSideEdge->GetInterval(),
                                                        m_dThisApproxTol3d,
                                                        bRunsOn,
                                                        sSolution ) );
          if ( ! bRunsOn && m_pExecutive->GetCliffOverflow() )
            {
              // The end face: the side edge's face that is neither of the fillet's side faces.
              SmFace * pEndFace = ( pSideFaces[0] == pOtherFace || pSideFaces[0] == pFilEdgeuse->GetFace()
                                    || pSideFaces[0] == pFilEdgeuse->GetRadial()->GetFace() )
                                  ? pSideFaces[1] : pSideFaces[0];
              if (   pEndFace == pFilEdgeuse->GetFace()
                  || pEndFace == pFilEdgeuse->GetRadial()->GetFace() )
                { SER( SM_ERR ); }
              SmFilletSolver * pFilSolver = GetFilletExecutive()->GetFilletSolverOfEdgeuse( pFilEdgeuse );
              NER(pFilSolver);
              SER( m_pExecutive->NoteOpenEnd( pFilSolver, SM_CONST_CAST(SmVertex*, m_cpVertex), pEndFace ));
              SER( pEndRib->GetBrep()->DeleteEdge( pEndRib ) );
              m_vEdges.ReSet();
            }
          else if ( ! bRunsOn && ! bBothMiss )
            {
              m_pExecutive->SetDoGlobalMerge( TRUE );
              SER( pEndRib->GetBrep()->DeleteEdge( pEndRib ) );
              m_vEdges.ReSet();
            }
        }
    }

  return SmFilletCorner::CalcCornerEdgeGeom();

} // end SmBlendCorner::CalcEndEdgeGeom

/*******************************************************************//**
PURPOSE: Create 'end' patch from extended 'side-face' of original brep
    if needed and create Face (and Edges and Verts) in m_pFilletBrep 
    from m_vSurfaces with call 
      SmFilletCorner::CreateCornerPatch() or
      SmFilletCorner::CreateEndPatch()

NOTES:
***********************************************************************/
SmStatus SmBlendCorner::CreateEndPatchOfEnd
  ()
{
  if (m_vEdges.GetSize() == 3)
    {
      if(m_bExtendedSurfacePatch)
           { SER(CreateEndPatch()); }
      else
           { SER(SmFilletCorner::CreateCornerPatch()); }
    }
  else
    { SM_DBG_WARN(_T("Unexpected End endFace construction case")) ; }

  return SM_SUCCESS;

} // end SmBlendCorner::CreateEndPatchOfEnd
