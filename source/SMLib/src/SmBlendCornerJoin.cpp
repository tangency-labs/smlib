// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

/**********************************************************************//**
* FILE NAME --- SmBlendCornerJoin.cpp
* PURPOSE: SolidCraft change. The Join event of the blend corner: two
*   blends of one convexity whose edges continue each other (G1) at a
*   vertex between the same two faces, so that the blends continue each
*   other too. They end on one shared section.
**********************************************************************/

#include <SmBlendCorner.h>
#include <SmFilletExecutive.h>
#include <SmFilletGeom.h>
#include <SmFilletSolver.h>
#include <SmEdgeuse.h>
#include <SmEdge.h>
#include <SmFace.h>

/*******************************************************************//**
PURPOSE: Whether the corner is a Join built on one shared section: its two
    blended edges continue each other (the Join) between the same two
    faces, and every other edge at the vertex is a seam of one face (the
    same face on both its sides: a periodic face's seam, which the blends
    roll across as across the face).
***********************************************************************/
SmBoolean SmBlendCorner::IsJoinSection() const
{
  if ( m_eEvent != SM_BCE_JOIN || m_pSolverEUs == NULL || m_pSolverEUs->GetSize() != 2 || m_pAllEUs == NULL )
    { return FALSE; }
  SmEdgeuse * pEU1 = (*m_pSolverEUs)[0];
  SmEdgeuse * pEU2 = (*m_pSolverEUs)[1];
  if ( pEU1 == NULL || pEU2 == NULL )
    { return FALSE; }
  const SmFace * apFaces1[2] = { pEU1->GetFace(), pEU1->GetRadial()->GetFace() };
  const SmFace * apFaces2[2] = { pEU2->GetFace(), pEU2->GetRadial()->GetFace() };
  if ( !(   ( apFaces1[0] == apFaces2[0] && apFaces1[1] == apFaces2[1] )
         || ( apFaces1[0] == apFaces2[1] && apFaces1[1] == apFaces2[0] )))
    { return FALSE; }
  for (ULONG ii=0; ii<m_pAllEUs->GetSize(); ii++)
    {
      SmEdgeuse * pEU = (*m_pAllEUs)[ii];
      if ( pEU == NULL )
        { return FALSE; }
      if ( pEU->GetEdge() == pEU1->GetEdge() || pEU->GetEdge() == pEU2->GetEdge() )
        { continue; }
      if ( pEU->GetFace() != pEU->GetRadial()->GetFace() )
        { return FALSE; }
    }
  return TRUE;

} // end SmBlendCorner::IsJoinSection

/*******************************************************************//**
PURPOSE: The Join's topology: the two blends end on one shared section.

NOTES: The blends roll on the same two faces A and B (IsJoinSection). The
   section is the cross-section of the first blend in the plane through
   the vertex square to its edge (SM_FV_ON_CROSS_SECTION, on A); its end
   on B is that vertex's mate (SM_FV_MATE), at the same section of the
   second blend. One corner edge, SM_FE_CROSS_SECTION, between them: the
   first blend's on its primary edgeuse, the second's on the mate. Each
   blend's rail on A ends at the first vertex, its rail on B at the mate.
   The blends meet there and nowhere else: the corner is blending (no
   blend is traced on past the vertex), with no patch of its own.
***********************************************************************/
SmStatus SmBlendCorner::MakeJoinTopology()
{
  SmEdgeuse      * pEU1 = (*m_pSolverEUs)[0];
  SmEdgeuse      * pEU2 = (*m_pSolverEUs)[1];
  SmFilletSolver * pFS1 = GetFilletExecutive()->GetFilletSolverOfEdgeuse( pEU1 );
  SmFilletSolver * pFS2 = GetFilletExecutive()->GetFilletSolverOfEdgeuse( pEU2 );
  NER( pFS1 );
  NER( pFS2 );
  SmFilletGeom   * pFG1 = pFS1->GetFirstFilletGeom();
  SmFilletGeom   * pFG2 = pFS2->GetFirstFilletGeom();
  NER( pFG1 );
  NER( pFG2 );
  SmFace * pFaceA = pEU1->GetFace();
  SmFace * pFaceB = pEU1->GetRadial()->GetFace();
  NER( pFaceA );
  NER( pFaceB );

  // the section's ends: on A, and its mate on B
  SmFilletVertex * pOnA = new (m_pPseudoBrep) SmFilletVertex( this );
  pOnA->SetPointClass( SM_PC_FACE, pFaceA );
  pOnA->SetFilletVertexType( SM_FV_ON_CROSS_SECTION );
  m_vVertices.Add( pOnA );
  SmFilletVertex * pOnB = new (m_pPseudoBrep) SmFilletVertex( this );
  pOnB->SetPointClass( SM_PC_FACE, pFaceB );
  pOnB->SetFilletVertexType( SM_FV_MATE );
  pOnB->SetMate( 0, pOnA );
  pOnA->SetMate( 0, pOnB );
  m_vVertices.Add( pOnB );

  // the shared section
  SmFilletEdge * pSection = NULL;
  SER( SmFilletEdge::MakeFilletEdge( m_pPseudoBrep, pOnA, pOnB, pSection, this ));
  NER( pSection );
  pSection->SetFilletEdgeType( SM_FE_CROSS_SECTION );
  SmFilletEdgeuse * pPrimEU = (SmFilletEdgeuse*)pSection->GetPrimaryEdgeuse();
  NER( pPrimEU );
  pPrimEU->SetFilletGeom( pFG1 );
  ((SmFilletEdgeuse*)pPrimEU->GetMate())->SetFilletGeom( pFG2 );
  m_vEdges.Add( pSection );

  // each blend's rails end at the section's ends
  SmFilletSolver * apFS[2] = { pFS1, pFS2 };
  SmFilletGeom   * apFG[2] = { pFG1, pFG2 };
  for (ULONG ii=0; ii<2; ii++)
    {
      ULONG lRailA = apFS[ii]->FindIndexOfRailOnFace( pFaceA );
      if ( lRailA > 1 )
        { SER( SM_ERR ); }
      SmFilletEdgeuse * pNewEdgeuse = NULL;
      SER( apFG[ii]->MakeRailEdgeuse( m_pPseudoBrep, lRailA,     pOnA, pNewEdgeuse ));
      SER( apFG[ii]->MakeRailEdgeuse( m_pPseudoBrep, 1 - lRailA, pOnB, pNewEdgeuse ));
    }

  m_bBlending = TRUE;
  m_bG1Case   = TRUE;
  return SM_SUCCESS;

} // end SmBlendCorner::MakeJoinTopology
