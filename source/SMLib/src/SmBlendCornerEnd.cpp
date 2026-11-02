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
#include <SmOffsetSurface.h>
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

  // an End at a smooth vertex: its plane; else the end face continued, or as it is
  SmSurface * pEndSurface = m_cpCorner->GetEndCapSurface();
  if ( pEndSurface == NULL )
    {
      SmFace * pEndFace = pRib->GetOriginalFace();
      NER( pEndFace );
      pEndSurface = m_cpCorner->GetExtendedSurface( pEndFace );
      if ( pEndSurface == NULL )
        { pEndSurface = pEndFace->GetSurface(); }
    }
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
  if ( m_eEnd == SM_BEND_AT_SMOOTH )
    {
      // the cap's edges on the side faces (they need only the vertices), then the end rib
      SER( CalcEndCapEdges() );
      return SmFilletCorner::CalcCornerEdgeGeom();
    }
  // A concave round's End of mixed convexity (its end rib, then the cap's two edges): a rail that
  // misses the side edge of the other convexity passes over it, and the round would run on past
  // the end face, adding material over it: refused, naming the largest radius whose rail still
  // ends on that edge (RefuseEndOverflowingSideEdge).
  if (   m_eEnd == SM_BEND_ON_END_FACE
      && m_vEdges.GetSize() == 3
      && m_vEdges[0] != NULL
      && m_vEdges[0]->GetFilletEdgeType() == SM_FE_FILLET_X_SIDE_FACE
      && m_vVertices.GetSize() >= 2 )
    {
      SmEdgeuse * pFilEdgeuse = (*m_pSolverEUs)[0];
      SmBoolean   bConcave    = FALSE;
      SER( SmFilletSolver::IsConcaveEdge( *pFilEdgeuse->GetEdge(), bConcave ));
      for (ULONG kk=0; kk<2 && bConcave; kk++)
        {
          SmFilletVertex * pFV = m_vVertices[kk];
          if (   pFV == NULL
              || pFV->GetStatus() != SM_FV_NO_INT_RAIL_EU
              || pFV->GetPointClass() != SM_PC_EDGEUSE )
            { continue; }
          // the side edge from the vertex, continued across seams' vertices up to where the
          // rail's end was last looked for (SmFilletVertex::CalcRailIntEdgeuse)
          SmTArray<SmEdgeuse*> sChain;
          SER( SideEdgeChainTo( (SmEdgeuse*)pFV->GetPointClassObject(), sChain ));
          SmEdgeuse * pSideEU = sChain[0];
          NER( pSideEU );
          if ( same_convexity_edges( pFilEdgeuse, pSideEU, m_cpVertex, 0 ))
            { continue; }
          SmFace * pSideFaces[2] = { pSideEU->GetFace(), pSideEU->GetRadial()->GetFace() };
          SmFace * pOtherFace    = pFilEdgeuse->GetFace();
          if ( pOtherFace == pSideFaces[0] || pOtherFace == pSideFaces[1] )
            { pOtherFace = pFilEdgeuse->GetRadial()->GetFace(); }
          NER( pOtherFace );
          SER( RefuseEndOverflowingSideEdge( pFilEdgeuse, sChain, pOtherFace ));
        }
    }
  SmFilletEdge * pEndRib = ( m_vEdges.GetSize() == 1 ) ? m_vEdges[0] : NULL;
  // An end on more than one end face (MakeEndOnEndFace: a corner of more than three edges of one
  // convexity) has its two rail ends and no end rib: the global merge cuts the fillet off where
  // it meets those faces, which holds only where both rails reach their side edges.
  SmBoolean bOnEndFaces = (   m_eEnd == SM_BEND_ON_END_FACE
                           && m_vEdges.GetSize() == 0
                           && m_vVertices.GetSize() == 2 );
  if (   bOnEndFaces
      || (   pEndRib != NULL
          && pEndRib->GetFilletEdgeType() == SM_FE_FILLET_X_SIDE_FACE ))
    {
      SmFilletVertex * pStartFV = bOnEndFaces ? m_vVertices[0] : (SmFilletVertex*)pEndRib->GetVertex();
      NER(pStartFV);
      SmFilletVertex * pEndFV   = bOnEndFaces ? m_vVertices[1] : (SmFilletVertex*)pEndRib->GetOtherVertex( pStartFV );
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
      // A rail's end found on the side edge's line continued past the edge's end
      // (SmFilletSolver::RailEdgeuseIntersect continues a straight side edge) misses the side
      // edge itself, however close past the end. In a cliff overflow a rail can also miss its
      // side edge with the intersection clamped to the side edge's far vertex, off the fillet;
      // and both rails can miss there.
      SmBoolean bCliff = m_pExecutive->GetCliffOverflow();
      SmSurface * pFilletSurface = NULL;
      if ( bCliff )
        {
          NER(pEndRib);
          SmFilletEdgeuse * pRibEU = (SmFilletEdgeuse*)pEndRib->GetPrimaryEdgeuse();
          SmFilletGeom    * pRibFG = ( pRibEU != NULL ) ? pRibEU->GetFilletGeom() : NULL;
          pFilletSurface = ( pRibFG != NULL ) ? pRibFG->GetFilletSurface() : NULL;
          NER(pFilletSurface);
        }
      SmFilletVertex * apFV[2] = { pStartFV, pEndFV };
      for (ULONG kk=0; kk<2; kk++)
        {
          if ( abMiss[kk] || !apFV[kk]->IsProcessed() )
            { continue; }
          SmBoolean bMisses = FALSE;
          if ( bCliff )
            {
              SmBoolean bDropped = FALSE, bMulti = FALSE;
              SmPoint2d sUV;
              double    dGap = SM_BIG_DOUBLE;
              SER( pFilletSurface->DropPoint( apFV[kk]->GetPoint(), pFilletSurface->GetNaturalUVDomain(), NULL, bDropped, sUV, dGap, bMulti ));
              bMisses = ( !bDropped || dGap > 10.0 * m_dThisApproxTol3d ) ? TRUE : FALSE;
            }
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
              if ( bCliff )
                { apFV[kk]->SetStatus( SM_FV_NO_INT_RAIL_EU ); }  // its rail is not trimmed there
            }
        }
      if ( ( bCliff || pMissFV == NULL ) && ( abMiss[0] || abMiss[1] ))
        { pMissFV = abMiss[0] ? pStartFV : pEndFV; }
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
              // The rail runs past its side edge's far end P onto the faces beyond. Refused, naming
              // the radius at which the rail reaches P: the ball touching the rail's face F at P
              // (its centre P offset along F's normal as the solver offsets F) that also touches
              // the fillet's other face G, the root of dist(centre(rho), G) = rho on (0, r].
              SER( RefuseEndPastSideEdge( pFilEdgeuse, pSideEU, pOtherFace ));
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
  if ( m_eEnd == SM_BEND_AT_SMOOTH )
    {
      NER( m_pEndCap );
      SmSurface * pCap = NULL;
      SER( m_pEndCap->CopyAndAddAnalytics( m_crContext, pCap ));
      m_vSurfaces.Add( pCap );
      if ( MakeFaceBrep() != SM_SUCCESS )
        {
          TCHAR sBuff[SM_TBLOCK_SIZE];
          GetName( sBuff, SM_TBLOCK_SIZE );
          GetFilletExecutive()->NoteFilletError( SM_FILERR_VERTEX_PROBLEM, GetCornerType(), SM_FIL_FAILURE, sBuff, GetFilletedVertex(),
                                                 NULL, NULL, NULL, NULL, _T("Fillet Corner Error: the cap of a round ending at a smooth vertex") );
          SER( SM_ERR );
        }
      return SM_SUCCESS;
    }
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


/*******************************************************************//**
PURPOSE: Whether a blended edge runs on at the vertex, tangentially, into an
    edge of its convexity that is not blended.

NOTES: With tangent propagation off, only the selected edges are rounded;
   where the selected edge runs on, along one of its faces, into a tangent
   edge left out (leaving the vertex the opposite way: two edges leaving it
   the same way meet in a cusp, a knife edge's tip), the round ends at that
   smooth vertex. (An edge running on
   between other faces, the seam of two other blends, is the end of the
   side faces: the End on the end faces past it takes it.) A tangent edge of the other convexity (a
   boss's base running on as the rim of an overhang past a wall) is no such
   continuation: the End on the end face past it takes it.
***********************************************************************/
SmBoolean SmBlendCorner::EndAtSmoothVertex
  (const SmVertex             * cpCornerVert,
   SmEdgeuse                  * pFilEdgeuse,
   const SmTArray<SmEdgeuse*> & crAllEUs,
   double                       dTanTolDeg)
{
  for (ULONG ii=0; ii<crAllEUs.GetSize(); ii++)
    {
      SmEdgeuse * pEU = crAllEUs[ii];
      if ( pEU == NULL || pEU->GetEdge() == pFilEdgeuse->GetEdge() )
        { continue; }
      // the edge running on along one of the blended edge's faces
      SmFace * pFaces[2] = { pEU->GetFace(), pEU->GetRadial()->GetFace() };
      SmBoolean bAlongItsFace = FALSE;
      for (ULONG ff=0; ff<2; ff++)
        {
          if ( pFaces[ff] == pFilEdgeuse->GetFace() || pFaces[ff] == pFilEdgeuse->GetRadial()->GetFace() )
            { bAlongItsFace = TRUE; }
        }
      // the other edge leaves the vertex opposite to the blended edge (it runs on; leaving the
      // same way, the two edges meet in a cusp, a knife edge's tip, which nothing continues)
      SmVector3d asLeave[2];
      SmEdge   * apEdge[2] = { pFilEdgeuse->GetEdge(), pEU->GetEdge() };
      SmBoolean  bLeaves = TRUE;
      for (ULONG ee=0; ee<2 && bLeaves; ee++)
        {
          if ( apEdge[ee]->GetStartVertex() == apEdge[ee]->GetEndVertex() )
            { bLeaves = FALSE; continue; }
          const SmBoolean bStart = ( apEdge[ee]->GetStartVertex() == cpCornerVert ) ? TRUE : FALSE;
          if ( apEdge[ee]->GetEndTangent( bStart, asLeave[ee] ) != SM_SUCCESS || asLeave[ee].Unitize() != SM_SUCCESS )
            { bLeaves = FALSE; continue; }
          if ( !bStart )
            { asLeave[ee] = -asLeave[ee]; }
        }
      if (   bAlongItsFace
          && bLeaves
          && asLeave[0].Dot( asLeave[1] ) < -cos( SM_DEG2RAD( dTanTolDeg ))
          && same_convexity_edges_round_vertex( pFilEdgeuse, pEU->GetEdge() ) )
        { return TRUE; }
    }
  return FALSE;

} // end SmBlendCorner::EndAtSmoothVertex

/*******************************************************************//**
PURPOSE: The End at a smooth vertex (SM_BEND_AT_SMOOTH): the blend is cut
    by the plane P through the vertex square to its edge there.

NOTES: Each rail ends where it crosses P (SM_FV_RAIL_X_EXTENDED_EDGEUSE,
   the curve-surface event on the rail as traced, CalcRailIntEndSurface:
   the blend is traced on past the vertex, the corner is no blending one);
   the end rib is the blend's intersection with P between the rails' ends
   (SM_FE_FILLET_X_SIDE_FACE, CalcFilletIntSideFace on the corner's cap).
   The cap is P's piece between the end rib and the side faces' traces in P
   from the rails' ends to the vertex (SM_FV_ON_VERTEX): those traces are
   each side face's intersection with P (exact on a plane: a line), edges
   of the side faces (CalcEndCapEdges). The faces and edges past the vertex
   stay as they are.
***********************************************************************/
SmStatus SmBlendCorner::MakeEndAtSmoothVertex()
{
  SmEdgeuse      * pFilEdgeuse = (*m_pSolverEUs)[0];
  SmFilletSolver * pFilSolver  = GetFilletExecutive()->GetFilletSolverOfEdgeuse( pFilEdgeuse );
  NER( pFilSolver );
  SmFilletGeom   * pFilletGeom = pFilSolver->GetFirstFilletGeom();
  NER( pFilletGeom );
  SmFace * apFace[2] = { pFilEdgeuse->GetFace(), pFilEdgeuse->GetRadial()->GetFace() };
  NER( apFace[0] );
  NER( apFace[1] );

  // P: through the vertex, square to the edge
  SmEdge   * pEdge = pFilEdgeuse->GetEdge();
  SmExtent1d sIvl  = pEdge->GetInterval();
  SmVector3d sPV[2];
  SER( pEdge->GetCurve()->Evaluate( ( pEdge->GetStartVertex() == m_cpVertex ) ? sIvl.GetMin() : sIvl.GetMax(), 1, TRUE, sPV ));
  SmPoint3d  sOrigin = m_cpVertex->GetPoint();
  SmVector3d sNormal = sPV[1];
  SER( sNormal.Unitize() );
  SmVector3d sX = ( smos_Fabs( sNormal.x ) < 0.9 ) ? SmVector3d( 1, 0, 0 ) : SmVector3d( 0, 1, 0 );
  sX = sX - sNormal * sX.Dot( sNormal );
  SER( sX.Unitize() );
  SmVector3d sY = sNormal * sX;
  SER( sY.Unitize() );
  // P's extent: past every point of the side faces from the vertex (the blend lies between them)
  double dReach = 0.0;
  for (ULONG ii=0; ii<2; ii++)
    {
      SmExtent3d sBox;
      SER( apFace[ii]->CalculateTightBoundingBox( sBox ));
      dReach = smos_Max( dReach, sBox.GetMin().DistanceBetween( sOrigin ));
      dReach = smos_Max( dReach, sBox.GetMax().DistanceBetween( sOrigin ));
    }
  dReach = 2.0 * dReach;
  SmVector2d sUVScale( 1.0, 1.0 );
  SmExtent2d sUVDomain( -dReach, -dReach, 2.0 * dReach, 2.0 * dReach );
  SmPlane * pPlane = new( m_crContext ) SmPlane( sOrigin, sX, sY, sUVScale, sUVDomain );
  NER( pPlane );
  SmStatus eNurb = pPlane->MakeNurb();
  if ( eNurb != SM_SUCCESS )
    { delete pPlane; SER( eNurb ); }
  delete m_pEndCap;
  m_pEndCap = pPlane;

  // the rails' ends on P, and the vertex
  SmFilletVertex * apEnd[2];
  for (ULONG ii=0; ii<2; ii++)
    {
      apEnd[ii] = new (m_pPseudoBrep) SmFilletVertex( this );
      apEnd[ii]->SetPointClass( SM_PC_FACE, apFace[ii] );
      apEnd[ii]->SetFilletVertexType( SM_FV_RAIL_X_EXTENDED_EDGEUSE );
      m_vVertices.Add( apEnd[ii] );
      ULONG lRail = pFilSolver->FindIndexOfRailOnFace( apFace[ii] );
      if ( lRail > 1 )
        { SER( SM_ERR ); }
      SmFilletEdgeuse * pNewEdgeuse = NULL;
      SER( pFilletGeom->MakeRailEdgeuse( m_pPseudoBrep, lRail, apEnd[ii], pNewEdgeuse ));
    }
  SmFilletVertex * pOnVertex = new (m_pPseudoBrep) SmFilletVertex( this );
  pOnVertex->SetPointClass( SM_PC_VERTEX, (SmObject*)m_cpVertex );
  pOnVertex->SetFilletVertexType( SM_FV_ON_VERTEX );
  m_vVertices.Add( pOnVertex );

  // the end rib, then the cap's edges on the side faces: a loop of the cap
  SmFilletEdge * pRib = NULL;
  SER( SmFilletEdge::MakeFilletEdge( m_pPseudoBrep, apEnd[0], apEnd[1], pRib, this ));
  NER( pRib );
  pRib->SetFilletEdgeType( SM_FE_FILLET_X_SIDE_FACE );
  ((SmFilletEdgeuse*)pRib->GetPrimaryEdgeuse())->SetFilletGeom( pFilletGeom );
  m_vEdges.Add( pRib );
  SmFilletEdge * pOnSecond = NULL;
  SER( SmFilletEdge::MakeFilletEdge( m_pPseudoBrep, apEnd[1], pOnVertex, pOnSecond, this ));
  NER( pOnSecond );
  pOnSecond->SetFilletEdgeType( SM_FE_PRECOMPUTED );
  pOnSecond->SetOriginalFace( apFace[1] );
  m_vEdges.Add( pOnSecond );
  SmFilletEdge * pOnFirst = NULL;
  SER( SmFilletEdge::MakeFilletEdge( m_pPseudoBrep, pOnVertex, apEnd[0], pOnFirst, this ));
  NER( pOnFirst );
  pOnFirst->SetFilletEdgeType( SM_FE_PRECOMPUTED );
  pOnFirst->SetOriginalFace( apFace[0] );
  m_vEdges.Add( pOnFirst );

  m_bBlending = FALSE;
  return SM_SUCCESS;

} // end SmBlendCorner::MakeEndAtSmoothVertex

/*******************************************************************//**
PURPOSE: An End at a smooth vertex: its cap's edges on the side faces, each
    the face's intersection with the cap's plane between the rail's end
    and the vertex.
***********************************************************************/
SmStatus SmBlendCorner::CalcEndCapEdges()
{
  NER( m_pEndCap );
  for (ULONG ii=1; ii<m_vEdges.GetSize(); ii++)
    {
      SmFilletEdge * pFE = m_vEdges[ii];
      NER( pFE );
      if ( pFE->GetFilletEdgeType() != SM_FE_PRECOMPUTED || pFE->GetCurve() != NULL )
        { continue; }
      SmFace         * pFace  = pFE->GetOriginalFace();
      NER( pFace );
      SmFilletVertex * pStart = (SmFilletVertex*)pFE->GetVertex();
      SmFilletVertex * pEnd   = (SmFilletVertex*)pFE->GetOtherVertex( pStart );
      NER( pStart );
      NER( pEnd );
      if ( !pStart->IsProcessed() || !pEnd->IsProcessed() )
        { SER( SM_ERR ); }
      // The rail's end must lie on the side face (a side face ending tangentially at the vertex,
      // as a round at its contact line, can have the rail roll over onto the face beyond).
      SmFilletVertex * pRailEnd = ( pStart->GetFilletVertexType() == SM_FV_ON_VERTEX ) ? pEnd : pStart;
      SmPointClassification sOn;
      SER( pFace->Point3DClassify( pRailEnd->GetPoint(), SmTol::GetZoneTol3d( pFace ), TRUE, sOn ));
      if ( sOn.GetPointClass() == SM_PC_UNKNOWN )
        {
          // The rail ends off its side face, past the face's tangent edge at the vertex: it rolls
          // over into the End, which is not built. That edge runs through the vertex, on the cap's
          // plane, so the ball touching the side face on it at the plane is the vertex itself: no
          // radius keeps the rail on its face up to the plane.
          SmFilletEdgeuse * pRibEU = (SmFilletEdgeuse*)m_vEdges[0]->GetPrimaryEdgeuse();
          SmFilletSolver  * pFS    = ( pRibEU != NULL && pRibEU->GetFilletGeom() != NULL ) ? pRibEU->GetFilletGeom()->GetFilletSolver() : NULL;
          NER( pFS );
          m_pExecutive->GetFilletErrorInfo()->NoteError( SM_FILERR_NO_ROOM,
              _T("Fillet Error: the round's contact line rolls over its face's tangent edge where the round ends at a smooth vertex (a rail rolling over into an end is not built): no radius builds it") );
          m_pExecutive->NoteFitRadius( 0.0, pFS->GetEdgeuse(0)->GetEdge() );
          SER( SM_ERR );
        }
      SmVector3d       sDir   = pEnd->GetPoint() - pStart->GetPoint();
      SmBSplineCurve * p3D = NULL, * pUV1 = NULL, * pUV2 = NULL;
      SER( sm_SrfSrfIntersection( m_crContext, pStart, pEnd, &sDir, NULL, pFace->GetSurface(), m_pEndCap,
                                  m_dThisApproxTol3d, GetTangencyTolRadians(), p3D, pUV1, pUV2, FALSE, NULL, 0 ));
      NER( p3D );
      delete pUV1;
      delete pUV2;
      pFE->SetCurve( p3D, TRUE );
      p3D->SetOwner( pFE );
    }
  return SM_SUCCESS;

} // end SmBlendCorner::CalcEndCapEdges

/*******************************************************************//**
PURPOSE: Whether an End's rail at the fillet's radius runs past its side
    edge's far end P, and the radius at which the rail reaches P.

NOTES: P is the side edge's end away from the vertex, F the fillet's side
   face the side edge bounds (the rail's), G the fillet's other side face.
   A ball of radius rho touching F at P has its centre at P offset by rho
   along F's normal, the way the solver offsets F; it touches G where its
   centre lies rho from G. d(rho) = dist(centre(rho), G) - rho is positive
   at 0 (P lies off G); where it is negative at the fillet's radius the
   ball there touching F at P cuts into G, so the rail lies beyond P and
   misses the side edge. d's root, by bisection on the cause's own
   geometry, is the radius.
***********************************************************************/
SmStatus SmBlendCorner::RadiusReachingSideEdgeEnd
  (SmEdgeuse * pFilEdgeuse,  // in : the filleted edgeuse
   SmEdgeuse * pSideEU,      // in : a side edge
   SmFace    * pOtherFace,   // in : the fillet's other side face (G)
   SmBoolean & rbPast,       // out: the rail at the fillet's radius runs past P
   double    & rdFit)        // out: with rbPast, the radius at which the rail reaches P
{
  SmEdge * pSideEdge = pSideEU->GetEdge();
  NER( pSideEdge );
  const SmVertex * pFar = ( pSideEdge->GetStartVertex() == m_cpVertex ) ? pSideEdge->GetEndVertex() : pSideEdge->GetStartVertex();
  NER( pFar );
  SmBoolean bOn = FALSE, bWithin = FALSE;
  SER( BallTouchingAt( pFilEdgeuse, pSideEU, pOtherFace, pFar->GetPoint(), bOn, bWithin, rdFit ));
  rbPast = ( !bOn && bWithin ) ? TRUE : FALSE;
  if ( !rbPast )
    { rdFit = 0.0; }
  return SM_SUCCESS;

} // end SmBlendCorner::RadiusReachingSideEdgeEnd

/*******************************************************************//**
PURPOSE: The ball touching the fillet's side face F at a point P of F (F
    the solver's face that the side edge bounds) and the fillet's other
    face G: its radius rho, the root of d(rho) = dist(P + rho n_F, G) - rho
    (its centre P offset along F's normal the way the solver offsets F), by
    bisection on (0, r], r the fillet's radius.

NOTES: rbOnG: P lies on G (d(0) <= 0, rho 0). rbWithin: the root lies in
    (0, r) (d(r) < 0), else rho is at least r (rdRho r).
***********************************************************************/
SmStatus SmBlendCorner::BallTouchingAt
  (SmEdgeuse       * pFilEdgeuse,  // in : the filleted edgeuse
   SmEdgeuse       * pSideEU,      // in : a side edge, on F
   SmFace          * pOtherFace,   // in : the fillet's other side face (G)
   const SmPoint3d & crP,          // in : the point on F
   SmBoolean       & rbOnG,        // out: P lies on G
   SmBoolean       & rbWithin,     // out: the root lies below the fillet's radius
   double          & rdRho)        // out: the ball's radius
{
  SmFilletSolver * pFS = GetFilletExecutive()->GetFilletSolverOfEdgeuse( pFilEdgeuse );
  NER( pFS );
  ULONG lF = 2;
  for (ULONG ii=0; ii<2; ii++)
    {
      SmFace * pFace = pFS->GetEdgeuse( ii )->GetFace();
      if ( pFace != pOtherFace && ( pFace == pSideEU->GetFace() || pFace == pSideEU->GetRadial()->GetFace() ))
        { lF = ii; }
    }
  if ( lF > 1 )
    { SER( SM_ERR ); }
  SmOffsetSurface * pOffsetF = pFS->GetSurface( lF );
  NER( pOffsetF );
  const SmSurface * pF = pOffsetF->GetBaseSurface();
  NER( pF );
  const SmSurface * pG = pOtherFace->GetSurface();
  NER( pG );
  const double dRadius = smos_Fabs( pOffsetF->GetOffsetDistance() );
  const double dSign   = ( pOffsetF->GetOffsetDistance() < 0.0 ) ? -1.0 : 1.0;
  SmBoolean bDropped = FALSE, bMulti = FALSE;
  SmPoint2d sUV;
  double    dGap = 0.0;
  SER( pF->DropPoint( crP, pF->GetNaturalUVDomain(), NULL, bDropped, sUV, dGap, bMulti ));
  if ( !bDropped )
    { SER( SM_ERR ); }
  SmVector3d sNormal;
  SER( pF->EvaluateNormal( sUV, TRUE, TRUE, sNormal ));
  const auto distance = [&]( double dRho, double & rdD ) -> SmStatus
    {
      SmPoint3d sCentre = crP + sNormal * ( dSign * dRho );
      SmBoolean bOn = FALSE, bM = FALSE;
      SmPoint2d sAt;
      double    dDist = 0.0;
      SER( pG->DropPoint( sCentre, pG->GetNaturalUVDomain(), NULL, bOn, sAt, dDist, bM ));
      if ( !bOn )
        { SER( SM_ERR ); }
      rdD = dDist - dRho;
      return SM_SUCCESS;
    };
  double dLo = 0.0, dHi = dRadius, dAtLo = 0.0, dAtHi = 0.0;
  SER( distance( dLo, dAtLo ));
  SER( distance( dHi, dAtHi ));
  rbOnG    = ( dAtLo > 0.0 ) ? FALSE : TRUE;
  rbWithin = ( dAtHi < 0.0 ) ? TRUE : FALSE;
  rdRho    = rbOnG ? 0.0 : dRadius;
  if ( rbOnG || !rbWithin )
    { return SM_SUCCESS; }
  for (ULONG it=0; it<100 && dHi - dLo > SM_EFF_ZERO * ( 1.0 + dRadius ); it++)
    {
      const double dMid = 0.5 * ( dLo + dHi );
      double dAtMid = 0.0;
      SER( distance( dMid, dAtMid ));
      if ( dAtMid > 0.0 ) { dLo = dMid; } else { dHi = dMid; }
    }
  rdRho = 0.5 * ( dLo + dHi );
  return SM_SUCCESS;

} // end SmBlendCorner::BallTouchingAt

/*******************************************************************//**
PURPOSE: The side edge at the vertex continued across seams' vertices up
    to pLast (sm_SideEdgeuseAcrossSeamVertex): its edgeuses in order from
    the vertex. pLast itself when it meets the vertex.

NOTES: The rail's end is looked for along that chain from the vertex
    (SmFilletVertex::CalcRailIntEdgeuse), pLast the edge it was last
    looked for on: walked back to the vertex, each step across the end of
    the edge that continues one not yet taken.
***********************************************************************/
SmStatus SmBlendCorner::SideEdgeChainTo
  (SmEdgeuse            * pLast,   // in : the side edgeuse the rail's end was last looked for on
   SmTArray<SmEdgeuse*> & rChain)  // out: the chain from the vertex to pLast
{
  NER( pLast );
  const double dTolDeg = SM_RAD2DEG( GetTangencyTolRadians() );
  SmTArray<SmEdgeuse*> sBack;
  sBack.Add( pLast );
  SmEdgeuse * pCurr = pLast;
  for (ULONG lStep=0; lStep<64; lStep++)
    {
      SmEdge * pE = pCurr->GetEdge();
      NER( pE );
      if ( pE->GetStartVertex() == m_cpVertex || pE->GetEndVertex() == m_cpVertex )
        {
          rChain.ReSet();
          for (ULONG ii=sBack.GetSize(); ii>0; ii--)
            { rChain.Add( sBack[ii-1] ); }
          return SM_SUCCESS;
        }
      SmEdgeuse * pPrev = NULL;
      const SmVertex * apEnds[2] = { pE->GetStartVertex(), pE->GetEndVertex() };
      for (ULONG kk=0; kk<2 && pPrev == NULL; kk++)
        {
          SmEdgeuse * pX = sm_SideEdgeuseAcrossSeamVertex( pCurr, apEnds[kk], dTolDeg );
          SmBoolean bTaken = FALSE;
          for (ULONG ii=0; pX != NULL && ii<sBack.GetSize(); ii++)
            { bTaken = ( bTaken || sBack[ii]->GetEdge() == pX->GetEdge() ) ? TRUE : FALSE; }
          if ( pX != NULL && !bTaken )
            { pPrev = pX; }
        }
      if ( pPrev == NULL )
        { SER( SM_ERR ); }
      sBack.Add( pPrev );
      pCurr = pPrev;
    }
  return SM_ERR;

} // end SmBlendCorner::SideEdgeChainTo

/*******************************************************************//**
PURPOSE: A concave round's End of mixed convexity whose rail misses the
    side edge of the other convexity (the rail passes over it): refused,
    SM_FILERR_NO_ROOM, naming the largest radius whose rail still ends on
    that edge.

NOTES: The rail on F ends on the side edge E where the ball touching F at
    a point of E also touches the fillet's other face G (BallTouchingAt):
    its radius rho over E runs from 0 where E meets G up to its most, past
    which no ball touching F on E fits, and the rail passes over E: the
    round would run on past the end face, adding material over it (a
    concave overflow, refused; owner, 2026-11-02). The most, by golden
    section of rho over E's interval on the cause's own geometry, is where
    the two Ends' rails meet on E: open, so less the part's own tolerance
    (coordinator, 2026-11-02: not the in-round fit margin). E is the side
    edge from the vertex continued across seams' vertices, as the rail's
    end was looked for (SideEdgeChainTo): the most over all its edges.
***********************************************************************/
SmStatus SmBlendCorner::RefuseEndOverflowingSideEdge
  (SmEdgeuse                  * pFilEdgeuse,  // in : the filleted edgeuse
   const SmTArray<SmEdgeuse*> & crSideEUs,    // in : the side edge the rail missed, from the vertex
                                              //      on across seams' vertices (SideEdgeChainTo)
   SmFace                     * pOtherFace)   // in : the fillet's other side face (G)
{
  double dMost = 0.0;
  for (ULONG lPiece=0; lPiece<crSideEUs.GetSize(); lPiece++)
    {
      SmEdgeuse * pSideEU = crSideEUs[lPiece];
      NER( pSideEU );
      SmEdge * pSideEdge = pSideEU->GetEdge();
      NER( pSideEdge );
      SmCurve * pCurve = pSideEdge->GetCurve();
      NER( pCurve );
      const SmExtent1d sIvl = pSideEdge->GetInterval();
      const auto rho = [&]( double dT, double & rdRho ) -> SmStatus
        {
          SmPoint3d sP;
          SER( pCurve->EvaluatePoint( dT, sP ));
          SmBoolean bOn = FALSE, bWithin = FALSE;
          SER( BallTouchingAt( pFilEdgeuse, pSideEU, pOtherFace, sP, bOn, bWithin, rdRho ));
          return SM_SUCCESS;
        };
      const double dGold = 0.5 * ( smos_Sqrt( 5.0 ) - 1.0 );
      double dA = sIvl.GetMin(), dB = sIvl.GetMax();
      double dC = dB - dGold * ( dB - dA ), dD = dA + dGold * ( dB - dA );
      double dAtC = 0.0, dAtD = 0.0;
      SER( rho( dC, dAtC ));
      SER( rho( dD, dAtD ));
      for (ULONG it=0; it<200 && dB - dA > SM_EFF_ZERO * ( 1.0 + sIvl.GetLength() ); it++)
        {
          if ( dAtC >= dAtD )
            {
              dB = dD; dD = dC; dAtD = dAtC;
              dC = dB - dGold * ( dB - dA );
              SER( rho( dC, dAtC ));
            }
          else
            {
              dA = dC; dC = dD; dAtC = dAtD;
              dD = dA + dGold * ( dB - dA );
              SER( rho( dD, dAtD ));
            }
        }
      double dPieceMost = 0.0;
      SER( rho( 0.5 * ( dA + dB ), dPieceMost ));
      dMost = smos_Max( dMost, dPieceMost );
    } // end for each piece of the side edge
  // the margin under the open bound: the part's own tolerance (the input's)
  NER( m_pExecutive->GetTargetBrep() );
  const double dFit = smos_Max( 0.0, dMost - m_pExecutive->GetTargetBrep()->GetTolerance().val );
  TCHAR sWhy[SM_TBLOCK_SIZE];
  smos_sprintf( sWhy, _T("Fillet Error: the round's contact line passes over its side edge where the round ends on a face (the round would run on past that face, adding material over it): a radius below %.6g fits"), dFit );
  m_pExecutive->GetFilletErrorInfo()->NoteError( SM_FILERR_NO_ROOM, sWhy );
  m_pExecutive->NoteFitRadius( dFit, pFilEdgeuse->GetEdge() );
  return SM_ERR;

} // end SmBlendCorner::RefuseEndOverflowingSideEdge

/*******************************************************************//**
PURPOSE: Refuse an End whose rail runs past its side edge's far end P,
    naming the radius at which the rail reaches P
    (RadiusReachingSideEdgeEnd).
***********************************************************************/
SmStatus SmBlendCorner::RefuseEndPastSideEdge
  (SmEdgeuse * pFilEdgeuse,  // in : the filleted edgeuse
   SmEdgeuse * pSideEU,      // in : the side edge the rail missed
   SmFace    * pOtherFace)   // in : the fillet's other side face (G)
{
  SmBoolean bPast = FALSE;
  double    dFit  = 0.0;
  SER( RadiusReachingSideEdgeEnd( pFilEdgeuse, pSideEU, pOtherFace, bPast, dFit ));
  if ( !bPast )
    { SER( SM_ERR ); }
  TCHAR sWhy[SM_TBLOCK_SIZE];
  SER( DescribeEndPastSideEdge( dFit, sWhy ));
  m_pExecutive->GetFilletErrorInfo()->NoteError( SM_FILERR_NO_ROOM, sWhy );
  m_pExecutive->NoteFitRadius( dFit, pFilEdgeuse->GetEdge() );
  return SM_ERR;

} // end SmBlendCorner::RefuseEndPastSideEdge

/*******************************************************************//**
PURPOSE: The reason of a refusal by RefuseEndPastSideEdge.
***********************************************************************/
SmStatus SmBlendCorner::DescribeEndPastSideEdge
  (double  dFit,   // in : the radius at which the rail reaches the side edge's far end
   TCHAR * cWhy)   // out: SM_TBLOCK_SIZE characters
{
  smos_sprintf( cWhy, _T("Fillet Error: the round's contact line runs past the end of its side edge at the corner (onto the faces beyond): a radius below %.6g fits"), dFit );
  return SM_SUCCESS;

} // end SmBlendCorner::DescribeEndPastSideEdge

/*******************************************************************//**
PURPOSE: Before the corner is built (the executive's regularity check
    refuses a blend before its corners): at an End on its end face, of
    a constant round, the least radius at which a rail reaches its side
    edge's far end, over the side edges the rail at the round's radius
    runs past (RadiusReachingSideEdgeEnd). Side edges of the other
    convexity or tangent at the vertex (whose rails end on an extended
    edge or at a tangent end) and a closed side edge are not sized here.
***********************************************************************/
SmStatus SmBlendCorner::FitOfEndPastSideEdges
  (SmBoolean & rbFound,   // out: some rail runs past its side edge's far end
   double    & rdFit,     // out: with rbFound, the least radius reaching it
   SmEdge   *& rpEdge)    // out: with rbFound, the blended edge
{
  rbFound = FALSE;
  rdFit   = 0.0;
  rpEdge  = NULL;
  if ( !IsEndOnEndFace() || m_pSolverEUs == NULL || m_pSolverEUs->GetSize() != 1 || m_pAllEUs == NULL )
    { return SM_SUCCESS; }
  SmEdgeuse      * pFilEdgeuse = (*m_pSolverEUs)[0];
  SmFilletSolver * pFS         = GetFilletExecutive()->GetFilletSolverOfEdgeuse( pFilEdgeuse );
  NER( pFS );
  if (   pFS->GetSolverType() != SM_FS_CONST_RADIUS
      && pFS->GetSolverType() != SM_FS_CONST_RADIUS_ASSISTED )
    { return SM_SUCCESS; }
  const ULONG lTotalEUs = m_pAllEUs->GetSize();
  ULONG lFil = SM_BIG_ULONG;
  for (ULONG ii=0; ii<lTotalEUs; ii++)
    {
      if ( (*m_pAllEUs)[ii] == pFilEdgeuse )
        { lFil = ii; break; }
    }
  if ( lFil == SM_BIG_ULONG || lTotalEUs < 3 )
    { return SM_SUCCESS; }
  SmEdgeuse * pEU[2] = { (*m_pAllEUs)[ ( lFil + lTotalEUs - 1 ) % lTotalEUs ], (*m_pAllEUs)[ ( lFil + 1 ) % lTotalEUs ] };
  if ( pEU[0] == pEU[1] )
    { return SM_SUCCESS; }
  const double dTanTolDeg = SM_RAD2DEG( GetTangencyTolRadians() );
  for (ULONG kk=0; kk<2; kk++)
    {
      if (   !same_convexity_edges( pFilEdgeuse, pEU[kk], m_cpVertex, 0 )
          || sm_CheckTangentialEdgeuseAtVertex( pEU[kk], m_cpVertex, dTanTolDeg ) )
        { return SM_SUCCESS; }
    }
  for (ULONG kk=0; kk<2; kk++)
    {
      SmEdgeuse * pSideEU = sm_FindSideEdgeuse( pFilEdgeuse, pEU[kk], 0 );
      NER( pSideEU );
      SmFace * pSideFaces[2] = { pSideEU->GetFace(), pSideEU->GetRadial()->GetFace() };
      SmFace * pOtherFace    = pFilEdgeuse->GetFace();
      if ( pOtherFace == pSideFaces[0] || pOtherFace == pSideFaces[1] )
        { pOtherFace = pFilEdgeuse->GetRadial()->GetFace(); }
      NER( pOtherFace );
      SmBoolean bPast = FALSE;
      double    dFit  = 0.0;
      SER( RadiusReachingSideEdgeEnd( pFilEdgeuse, pSideEU, pOtherFace, bPast, dFit ));
      if ( bPast && ( !rbFound || dFit < rdFit ))
        {
          rbFound = TRUE;
          rdFit   = dFit;
          rpEdge  = pFilEdgeuse->GetEdge();
        }
    }
  return SM_SUCCESS;

} // end SmBlendCorner::FitOfEndPastSideEdges
