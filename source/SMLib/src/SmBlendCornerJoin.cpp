// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

/**********************************************************************//**
* FILE NAME --- SmBlendCornerJoin.cpp
* PURPOSE: SolidCraft change. The Join event of the blend corner: blended
*   edges that continue each other (G1) at a vertex, two of them or one
*   closed edge meeting itself, so that the blends continue each other too.
*   They end on one shared section, the section of the first blend in the
*   plane through the vertex square to its edge; each blend's rails end at
*   that section's two ends.
**********************************************************************/

#include <SmBlendCorner.h>
#include <SmFilletExecutive.h>
#include <SmFilletGeom.h>
#include <SmFilletSolver.h>
#include <SmOffsetSurface.h>
#include <SmEdgeuse.h>
#include <SmEdge.h>
#include <SmFace.h>
#include <SmCurve.h>

/*******************************************************************//**
PURPOSE: Whether the corner is a Join of one closed blended edge meeting
    itself at the vertex.
***********************************************************************/
SmBoolean SmBlendCorner::IsClosedJoin() const
{
  return ( m_eEvent == SM_BCE_JOIN && m_pSolverEUs != NULL && m_pSolverEUs->GetSize() == 1 ) ? TRUE : FALSE;
}

/*******************************************************************//**
PURPOSE: The faces on each side of each blended edgeuse at a Join, and the
    edgeuses between the blends on each side.

NOTES: Going round the vertex from the first blended edgeuse to the
   second, the first side's faces lie between them, the second side's from
   the second back to the first. Each blend's face on a side is the one it
   shares with its neighbour on that side: the same face for both blends
   when nothing stands between them there, else faces meeting tangentially
   along the edges between them. With two edges in all, or one closed
   edge, both blends roll on the first edgeuse's face (first side) and its
   radial's (second side).
***********************************************************************/
static SmStatus sm_JoinSides
  (const SmTArray<SmEdgeuse*> & crSolverEUs,    // in : the blended edgeuses (one, or two)
   const SmTArray<SmEdgeuse*> & crAllEUs,       // in : every edgeuse round the vertex, in order
   SmFace                     * apSide[2][2],   // out: [side][blend] the blend's face on that side
   SmTArray<SmEdgeuse*>         asBetween[2])   // out: [side] the edgeuses between the blends there
{
  SmEdgeuse * pEU1 = crSolverEUs[0];
  apSide[0][0] = apSide[0][1] = pEU1->GetFace();
  apSide[1][0] = apSide[1][1] = pEU1->GetRadial()->GetFace();
  if ( crSolverEUs.GetSize() == 1 )
    {
      // one closed edge meeting itself: the other edges lie on its faces' sides by their faces
      for (ULONG ii=0; ii<crAllEUs.GetSize(); ii++)
        {
          SmEdgeuse * pEU = crAllEUs[ii];
          if ( pEU->GetEdge() == pEU1->GetEdge() )
            { continue; }
          SmFace * pFaces[2] = { pEU->GetFace(), pEU->GetRadial()->GetFace() };
          if ( pFaces[0] == apSide[0][0] || pFaces[1] == apSide[0][0] )
            { asBetween[0].Add( pEU ); }
          else if ( pFaces[0] == apSide[1][0] || pFaces[1] == apSide[1][0] )
            { asBetween[1].Add( pEU ); }
          else
            { SER( SM_ERR ); }
        }
      return SM_SUCCESS;
    }
  if ( crAllEUs.GetSize() <= 2 )
    { return SM_SUCCESS; }
  SmEdgeuse * pEU2 = crSolverEUs[1];
  ULONG lTotal = crAllEUs.GetSize(), lIdx1 = SM_BIG_ULONG, lIdx2 = SM_BIG_ULONG;
  for (ULONG ii=0; ii<lTotal; ii++)
    {
      if ( crAllEUs[ii]->GetEdge() == pEU1->GetEdge() ) { lIdx1 = ii; }
      if ( crAllEUs[ii]->GetEdge() == pEU2->GetEdge() ) { lIdx2 = ii; }
    }
  if ( lIdx1 == SM_BIG_ULONG || lIdx2 == SM_BIG_ULONG )
    { SER( SM_ERR ); }
  for (ULONG kk=( lIdx1 + 1 ) % lTotal; kk!=lIdx2; kk=( kk + 1 ) % lTotal)
    { asBetween[0].Add( crAllEUs[kk] ); }
  for (ULONG kk=( lIdx2 + 1 ) % lTotal; kk!=lIdx1; kk=( kk + 1 ) % lTotal)
    { asBetween[1].Add( crAllEUs[kk] ); }
  // a side with edges between the blends: each blend's face is the one it shares with its
  // neighbour there; its face on the other side is its other face
  SmEdgeuse * apEU[2] = { pEU1, pEU2 };
  const ULONG lKnown = ( asBetween[0].GetSize() > 0 ) ? 0 : 1;
  if ( asBetween[lKnown].GetSize() == 0 )
    { SER( SM_ERR ); }
  SmEdgeuse * apNeighbour[2] = { ( lKnown == 0 ) ? asBetween[0][0] : asBetween[1].GetLast(),
                                 ( lKnown == 0 ) ? asBetween[0].GetLast() : asBetween[1][0] };
  for (ULONG bb=0; bb<2; bb++)
    {
      SmFace * pKnown = sm_GetFaceByTwoEdgeuses( apEU[bb], apNeighbour[bb], NULL );
      NER( pKnown );
      SmFace * pFaces[2] = { apEU[bb]->GetFace(), apEU[bb]->GetRadial()->GetFace() };
      if ( pKnown != pFaces[0] && pKnown != pFaces[1] )
        { SER( SM_ERR ); }
      apSide[lKnown][bb]     = pKnown;
      apSide[1 - lKnown][bb] = ( pKnown == pFaces[0] ) ? pFaces[1] : pFaces[0];
    }
  return SM_SUCCESS;

} // end sm_JoinSides

/*******************************************************************//**
PURPOSE: The Join's topology: the blends end on one shared section.

NOTES: One corner edge, SM_FE_CROSS_SECTION, the shared section: the first
   blend's on its primary edgeuse, the second's on the mate (the one
   blend's on both, when one closed edge meets itself). Each blend's rail
   on a side ends at the section's end on that side. The blends meet there
   and nowhere else: the corner is blending (no blend is traced on past the
   vertex), with no patch of its own. Where nothing stands between the
   blends, the section is the one through the vertex: its end on the first
   side SM_FV_JOIN_SECTION, on the other its mate (SM_FV_MATE). Where one
   edge stands between them on a side (the faces change there, meeting
   tangentially, or a face's seam runs from the vertex), the section is
   where the rails cross that edge (each blend's faces reach that far, and
   no farther: the section through the vertex can lie beyond them): its
   end there on the edge (SM_FV_RAIL_X_EDGEUSE), its mate on the other
   side (or, with such an edge on both sides, where the rails cross the
   other one too). More than one edge on a side: refused.
***********************************************************************/
SmStatus SmBlendCorner::MakeJoinTopology()
{
  const ULONG lBlends = m_pSolverEUs->GetSize();
  if ( lBlends != 1 && lBlends != 2 )
    { SER( SM_ERR ); }
  SmFace * apSide[2][2] = { { NULL, NULL }, { NULL, NULL } };
  SmTArray<SmEdgeuse*> asBetween[2];
  SER( sm_JoinSides( *m_pSolverEUs, *m_pAllEUs, apSide, asBetween ));
  // the section's ends: on a side with an edge between the blends (where they change face, or
  // cross a face's seam) where the rails cross it, else (the first side) through the vertex
  if ( asBetween[0].GetSize() > 1 || asBetween[1].GetSize() > 1 )
    {
      GetFilletExecutive()->GetFilletErrorInfo()->NoteError( m_cpVertex, SM_FIL_FAILURE, SM_FILERR_CORNER_UNSUPPORTED,
          _T("Fillet Error: blends continuing each other across more than one edge on a side") );
      SER( SM_ERR );
    }
  const SmBoolean abOnEdge[2] = { asBetween[0].GetSize() == 1, asBetween[1].GetSize() == 1 };
  const ULONG     lFirst      = ( abOnEdge[0] || !abOnEdge[1] ) ? 0 : 1;
  const ULONG     lOther      = 1 - lFirst;

  SmFilletSolver * apFS[2] = { NULL, NULL };
  SmFilletGeom   * apFG[2] = { NULL, NULL };
  for (ULONG ii=0; ii<lBlends; ii++)
    {
      apFS[ii] = GetFilletExecutive()->GetFilletSolverOfEdgeuse( (*m_pSolverEUs)[ii] );
      NER( apFS[ii] );
      apFG[ii] = apFS[ii]->GetFirstFilletGeom();
      NER( apFG[ii] );
    }

  // the section's ends
  SmFilletVertex * apEnd[2] = { NULL, NULL };
  const ULONG      alSide[2] = { lFirst, lOther };
  for (ULONG ee=0; ee<2; ee++)
    {
      const ULONG lSide = alSide[ee];
      apEnd[ee] = new (m_pPseudoBrep) SmFilletVertex( this );
      if ( abOnEdge[lSide] )
        {
          SmEdgeuse * pSideEU = sm_FindSideEdgeuse( (*m_pSolverEUs)[0], asBetween[lSide][0], 0 );
          NER( pSideEU );
          apEnd[ee]->SetPointClass( SM_PC_EDGEUSE, pSideEU );
          apEnd[ee]->SetFilletVertexType( SM_FV_RAIL_X_EDGEUSE );
        }
      else
        {
          apEnd[ee]->SetPointClass( SM_PC_FACE, apSide[lSide][0] );
          apEnd[ee]->SetFilletVertexType( ee == 0 ? SM_FV_JOIN_SECTION : SM_FV_MATE );
        }
      m_vVertices.Add( apEnd[ee] );
    }
  SmFilletVertex * pFirst = apEnd[0];
  SmFilletVertex * pOther = apEnd[1];
  if ( !abOnEdge[lOther] )
    {
      pOther->SetMate( 0, pFirst );
      pFirst->SetMate( 0, pOther );
    }

  // the shared section
  SmFilletEdge * pSection = NULL;
  SER( SmFilletEdge::MakeFilletEdge( m_pPseudoBrep, pFirst, pOther, pSection, this ));
  NER( pSection );
  pSection->SetFilletEdgeType( SM_FE_CROSS_SECTION );
  SmFilletEdgeuse * pPrimEU = (SmFilletEdgeuse*)pSection->GetPrimaryEdgeuse();
  NER( pPrimEU );
  pPrimEU->SetFilletGeom( apFG[0] );
  ((SmFilletEdgeuse*)pPrimEU->GetMate())->SetFilletGeom( apFG[ lBlends - 1 ] );
  m_vEdges.Add( pSection );

  // the rails' ends
  for (ULONG ii=0; ii<lBlends; ii++)
    {
      ULONG lRail = apFS[ii]->FindIndexOfRailOnFace( apSide[lFirst][ii] );
      if ( lRail > 1 )
        { SER( SM_ERR ); }
      SmFilletEdgeuse * pNewEdgeuse = NULL;
      if ( lBlends == 1 )
        {
          SER( apFG[ii]->MakeRailEdgeuse( m_pPseudoBrep, lRail,     pFirst, pNewEdgeuse, SM_OT_SAME ));
          SER( apFG[ii]->MakeRailEdgeuse( m_pPseudoBrep, lRail,     pFirst, pNewEdgeuse, SM_OT_OPPOSITE ));
          SER( apFG[ii]->MakeRailEdgeuse( m_pPseudoBrep, 1 - lRail, pOther, pNewEdgeuse, SM_OT_SAME ));
          SER( apFG[ii]->MakeRailEdgeuse( m_pPseudoBrep, 1 - lRail, pOther, pNewEdgeuse, SM_OT_OPPOSITE ));
        }
      else
        {
          SER( apFG[ii]->MakeRailEdgeuse( m_pPseudoBrep, lRail,     pFirst, pNewEdgeuse ));
          SER( apFG[ii]->MakeRailEdgeuse( m_pPseudoBrep, 1 - lRail, pOther, pNewEdgeuse ));
        }
    }

  m_bBlending = TRUE;
  m_bG1Case   = TRUE;
  return SM_SUCCESS;

} // end SmBlendCorner::MakeJoinTopology

/*******************************************************************//**
PURPOSE: A Join's shared section (SM_FV_JOIN_SECTION and its mate): each
    blend's section in the plane through the corner's vertex square to the
    first blend's edge there.

NOTES: Each blend's section there is its ball's: the point of the spine in
   the plane (SmFilletSolver::PointOnPlaneSolve, on the blend's own faces),
   started from the ball's feet at the vertex (the vertex stepped the
   ball's radius along each face's binormal of the edge). That section
   sets the blend's rail ends at this corner. The section's ends are the
   first blend's feet: on this vertex's face and on its mate's.
***********************************************************************/
SmStatus SmFilletVertex::CalcJoinSection()
{
  NER( m_cpCorner );
  SmFilletVertex * pMate = GetMate( 0 );
  NER( pMate );
  SmFace * pFirstFace = (SmFace*)GetPointClassObject();
  SmFace * pOtherFace = (SmFace*)pMate->GetPointClassObject();
  NER( pFirstFace );
  NER( pOtherFace );
  const SmVertex * cpVertex = m_cpCorner->GetFilletedVertex();
  NER( cpVertex );

  SmTArray<SmFilletGeom*> sGeoms;
  GetFilletGeoms( sGeoms );
  if ( sGeoms.GetSize() < 1 )
    { SER( SM_ERR ); }

  // the plane: through the vertex, square to the first blend's edge
  SmPoint3d  sPlaneOrig;
  SmVector3d sPlaneNormal;
  for (ULONG gg=0; gg<sGeoms.GetSize(); gg++)
    {
      SmFilletSolver * pFS = sGeoms[gg]->GetFilletSolver();
      NER( pFS );
      SmEdgeuse * pEU   = pFS->GetEdgeuse( 0 );
      NER( pEU );
      SmEdge    * pEdge = pEU->GetEdge();
      NER( pEdge );
      SmCurve   * pCurve = pEdge->GetCurve();
      NER( pCurve );
      SmExtent1d sIvl   = pEdge->GetInterval();
      double     dParam = ( pEdge->GetStartVertex() == cpVertex ) ? sIvl.GetMin() : sIvl.GetMax();
      if ( gg == 0 )
        {
          SmVector3d sPV[2];
          SER( pCurve->Evaluate( dParam, 1, TRUE, sPV ));
          sPlaneOrig   = sPV[0];
          sPlaneNormal = sPV[1];
        }

      // the ball's feet at the vertex, on each of the blend's faces
      SmVector2d sUVs[2];
      for (ULONG ii=0; ii<2; ii++)
        {
          SmEdgeuse       * pSideEU = pFS->GetEdgeuse( ii );
          NER( pSideEU );
          SmOffsetSurface * pOffset = pFS->GetSurface( ii );
          NER( pOffset );
          SmSurface       * pBase   = pOffset->GetBaseSurface();
          NER( pBase );
          SmPoint3d  sPnt;
          SmVector3d sBin;
          SmPoint2d  sUV;
          SER( pSideEU->EvaluateBinormal( dParam, FALSE, sPnt, sBin, NULL, NULL, &sUV ));
          sBin = sBin * smos_Fabs( pOffset->GetOffsetDistance() );
          SmVector2d sUVDir;
          SER( pBase->DropVectors( sUV, TRUE, TRUE, 1, &sBin, &sUVDir ));
          sUVs[ii] = SmVector2d( sUV.x + sUVDir.x, sUV.y + sUVDir.y );
        }

      SmFilletVertexuse * pVU = GetVUAtRailEnd( sGeoms[gg] );
      NER( pVU );
      SmTsectPnt & rTsectPnt = pVU->GetTsectPnt();
      SmBoolean bFound = FALSE;
      SER( pFS->PointOnPlaneSolve( sPlaneOrig, sPlaneNormal,
                                   pFS->GetSurface( 0 )->GetNaturalUVDomain(), pFS->GetSurface( 1 )->GetNaturalUVDomain(),
                                   sUVs[0], sUVs[1], bFound, rTsectPnt ));
      if ( !bFound )
        {
          SetStatus( SM_FV_SOLVER_NOT_CONVERGE );
          SER_MSG( SM_ERR, _T("CalcJoinSection(): no ball of the blend lies in the shared section's plane") );
        }
      // a closed blend's other rail end here
      SmFilletVertexuse * pOtherVU = GetVUAtRailEnd( sGeoms[gg], pVU );
      if ( pOtherVU != NULL )
        { pOtherVU->GetTsectPnt() = rTsectPnt; }

      if ( gg == 0 )
        {
          // the section's ends: the first blend's feet
          ULONG     lFirst   = ( pEU->GetFace() == pFirstFace ) ? 0 : 1;
          SmPoint2d sUVFirst = rTsectPnt.UVPos( lFirst );
          SmPoint2d sUVOther = rTsectPnt.UVPos( 1 - lFirst );
          SmPoint3d sFirst, sOther;
          SER( pFirstFace->GetSurface()->EvaluatePoint( sUVFirst, sFirst ));
          SER( pOtherFace->GetSurface()->EvaluatePoint( sUVOther, sOther ));
          SetOriginalUV( sUVFirst );
          SetPoint( sFirst );
          SetStatus( SM_FIL_PROCESSED );
          pMate->SetOriginalUV( sUVOther );
          pMate->SetPoint( sOther );
          pMate->SetStatus( SM_FIL_PROCESSED );
        }
    }
  return SM_SUCCESS;

} // end SmFilletVertex::CalcJoinSection
