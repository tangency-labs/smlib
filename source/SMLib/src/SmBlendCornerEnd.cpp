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

#include <SmFilletGeom.h>
#include <SmFilletCorner.h>
#include <SmFilletSolver.h>
#include <SmSolutionArray.h>
#include <SmSurface.h>
#include <SmCurve.h>
#include <SmVertex.h>

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
