// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

/**********************************************************************//**
* FILE NAME --- SmBlendCorner.cpp
* PURPOSE: SolidCraft change. The blend corner (SmBlendCorner.h): its
*   classification of a filleted vertex into a corner event, and the Swept
*   turn.
*
* The Swept turn at a straight convex side edge E (unit direction e, through
* P): each round i ends on its section through the side vertex S, square to
* its edge, the ball there centred at C_i at distance r_i from its faces.
* With O_i the foot of C_i on E and R_i = |C_i - O_i|, the ball turning about
* E from C_1 to C_2 keeps touching the common face (square to E, as both
* edges are) and E: its centre runs on the circle about O = O_1 = O_2 of
* radius R = R_1 = R_2, and the envelope is the torus
*   T(u, v) = O + (R + r cos v)(cos u X + sin u Y) + r sin v Z,
* X = (C_1 - O) / R, Z = X x (C_2 - O) / |X x (C_2 - O)| (= +-e), Y = Z x X,
* u in [0, theta] (theta the angle from C_1 - O to C_2 - O), v over the end
* section's arc from S (v = 180 degrees: S = O when R = r, the horn torus)
* to the plate vertex A. Along u = 0 and u = theta the torus holds the end
* sections, its normal there radial from C_i in their planes as the rounds'
* is: G1 to both. Along v = v_A it is the circle of the plate vertices about
* O's foot on the plate, its normal the plate's: G1 to the plate.
**********************************************************************/

#include "StdAfx.h"

#ifndef __SMBLENDCORNER_H__
#include <SmBlendCorner.h>
#endif

#ifndef __SMFILLETSOLVER_H__
#include <SmFilletSolver.h>
#endif

#include <SmAxis2Placement.h>
#include <SmOffsetSurface.h>
#include <SmTorus.h>

#ifdef SM_GFX_CODE
  #include <SmGraphicsOutput.h>
#endif // SM_GFX_CODE

/*******************************************************************//**
    The shared corner helpers (declared in SmBlendCorner.h)
***********************************************************************/

/********************************************************************
PURPOSE: Given a side edgeuse determine whether itself or its radial
            is the one which FilletSolver use in solving Rail_x_SideEdge
NOTES:
********************************************************************/
SmEdgeuse * sm_FindSideEdgeuse      // rtn: edgeuse
 (SmEdgeuse * pFilEdgeuse,           // in : target fillet edgeuse
  SmEdgeuse * pSideEU,               // in : target sideEdgeuse - connects to a filleted edge at an filletCorner
  int         iDebugLevel)
{
  // locals - filleted edge and CW neighbor to SideEU
  SmEdge    * pFilletedEdge = pFilEdgeuse->GetEdge();
  SmEdgeuse * pEU           = pSideEU;
  SmEdge    * pE            = pEU->GetCWEdgeuse()->GetEdge();

  // when edgeuses connect to different edges
  if (pFilletedEdge != pE)
    {
      // check the radial edgeuse
      pEU = pEU->GetRadial();
      pE  = pEU->GetCCWEdgeuse()->GetEdge();
      if (pFilletedEdge != pE)
          return NULL;
    }

#ifdef SM_DEBUG_CODE
  if ( iDebugLevel )
    {
      smgfx_SetLook(1,2, 1,0,0); pFilEdgeuse->Draw() ; sm_GraphicsLoop();
      smgfx_SetLook(1,2, 1,0,0); pEU->Draw(); sm_GraphicsLoop();
      sm_GraphicsLoop();
    }
#else
    SM_REF1(iDebugLevel);
#endif // SM_DEBUG_CODE

  return pEU;

} // end sm_FindSideEdgeuse

/********************************************************************
PURPOSE:  Determine if two adjacent edges are tangent to each other
             at a given vertex
NOTES:
********************************************************************/
SmBoolean sm_TestG1Continuity        // rtn: TRUE = Edge EndTangents at vertex are within tolerance of being tangent
  (const SmVertex * pVert,           // in : Target Vertex
   const SmEdge   * pEdge1,          // in : 1st edge connecting to vertex
   const SmEdge   * pEdge2,          // in : 2nd edge connectint to vertex
   double           dAngTolDeg,      // in : Angle tolerance in degrees
   double         * pdAngDeg,        // out: angle (degrees) between edges at vertex
   int              iDebugLevel)     // in :
{
  // no work - not given 3 valid pointers
  if (   pVert  == NULL
      || pEdge1 == NULL
      || pEdge2 == NULL)
   { return FALSE ; }

  // locals
  SmVector3d  sPV1[2] ;
  SmCurve   * pCrv1   = pEdge1->GetCurve() ;
  SmExtent1d  sIvl1   = pEdge1->GetInterval() ;
  double      dParam1 = (pEdge1->GetStartVertex() == pVert) ? sIvl1.GetMin() : sIvl1.GetMax() ;

  // get Edge1 tangent at vertex
  // SolidCraft change: a tangent that cannot be evaluated is no tangency
  if (pCrv1->Evaluate(dParam1,1,TRUE,sPV1) != SM_SUCCESS)
   { return FALSE ; }

  // locals
  SmVector3d  sPV2[2];
  SmCurve   * pCrv2   = pEdge2->GetCurve();
  SmExtent1d  sIvl2   = pEdge2->GetInterval();
  double      dParam2 = (pEdge2 == pEdge1 || pEdge2->GetStartVertex() != pVert) ? sIvl2.GetMax() : sIvl2.GetMin() ;

  // get Edge2 tangent at vertex
  if (pCrv2->Evaluate(dParam2,1,TRUE,sPV2) != SM_SUCCESS)
   { return FALSE ; }

#ifdef SM_DEBUG_CODE
    if ( iDebugLevel > 0 )
      {
        smgfx_SetLook(1,6, 1,0,0); pVert->Draw(); sm_GraphicsLoop();
        sPV1[1].Draw(sPV1); sm_GraphicsLoop();
        sPV2[1].Draw(sPV2); sm_GraphicsLoop();
        sm_GraphicsLoop();
      }
#else
  SM_REF1(iDebugLevel);
#endif // SM_DEBUG_CODE

  if (   sPV1[1].LengthSquared() > SM_EFF_ZERO_SQ
      && sPV2[1].LengthSquared() > SM_EFF_ZERO_SQ)
    {
      // get angle (degrees) between Edge1 and Edge2 endTangents
      double dAngDeg;
      if (sPV1[1].AngleBetween(sPV2[1],dAngDeg) != SM_SUCCESS) // SolidCraft: an angle that cannot be measured is no tangency
        { return FALSE; }
      if (dAngDeg > SM_PI/2.0) dAngDeg = SM_PI - dAngDeg;
      dAngDeg = SM_RAD2DEG(dAngDeg);

      // set optional output
      if (pdAngDeg) { *pdAngDeg = dAngDeg; }

      // check for in tolerance angle
      if (dAngDeg < dAngTolDeg)
        {
          return TRUE;
        }
    } // end nonDegenerate endTangent check

  // endTangents are not within tolerance
  return FALSE;

} // end sm_TestG1Continuity

/********************************************************************
PURPOSE:  return face connected to sectors of two given edgeuses,
             return NULL when no such face exists.

NOTES:
   use pEU3 to specify a face to skip, this forces the function
   to return another face that connects both edges or NULL.

   When pEU3 != NULL, the function returns NULL when the
     edges are not connected to a common face or when
     they are only connected to one common face which
     is equal to pEU3->Face.

METHOD ---
  4 faces are examined pEU1->SectorBeginFace, pEU1->SectorEndFace
                       pEU2->SectorBeginFace, pEU2->SectorEndFace
  and the common face is returned.
  Note in nonManifold modeling not all faces attaced to pEU->edges are examined.

********************************************************************/
SmFace * sm_GetFaceByTwoEdgeuses
 (SmEdgeuse * pEU1,       // in : bounding edge 1
  SmEdgeuse * pEU2,       // in : bounding edge 2
  SmEdgeuse * pEU3)       // in : opt edgeuse to specify a face to ignore.
                          //      function will not return pEU3->Face.
                          //      NULL to ignore.
{
  if (   pEU1 == NULL
      || pEU2 == NULL)
    { return NULL; }

  // locals
  SmTArray<SmFace*> sFaces1;
  SmTArray<SmFace*> sFaces2;

  // let sFaces1 = two faces attached to pEU1->Sector
  sFaces1.Add(pEU1->GetFace());
  sFaces1.Add(pEU1->GetRadial()->GetFace());

  // let sFaces2 = two faces attached to pEU2->Sector
  sFaces2.Add(pEU2->GetFace());
  sFaces2.Add(pEU2->GetRadial()->GetFace());

  // for every pEU1->Sector->Face
  ULONG ii, jj;
  for (ii=0; ii<sFaces1.GetSize(); ii++)
    {
      SmFace * pFace1 = sFaces1[ii];

      // when asked skip pEU3->Face
      if (pEU3 && pEU3->GetFace() == pFace1)
        { continue; }

      // for every pEU2->Sector->Face
      for (jj=0; jj<sFaces2.GetSize(); jj++)
        {
          SmFace * pFace2 = sFaces2[jj];

          // return 1st found common face
          if (pFace1 == pFace2)
            {
              return pFace2;
            }
        } // end iter every pEU2->Sector->Face
    } // end iter every pEU1->Sector->Face

  return NULL;

} // end sm_GetFaceByTwoEdgeuses

/*******************************************************************//**
    Local helpers
***********************************************************************/
namespace
{

// A turn the blend corner does not build: noted at the corner's vertex, then SM_ERR.
SmStatus sm_TurnNotBuilt(SmFilletExecutive * pExec, const SmVertex * cpVertex, SmFilletErrorType eError, const TCHAR * cpWhy)
{
  pExec->GetFilletErrorInfo()->NoteError( cpVertex, SM_FIL_FAILURE, eError, cpWhy );
  return SM_ERR;
}

// The unit tangent of pEdge's curve at its end on cpVertex.
SmStatus sm_TangentAtVertex(const SmEdge * pEdge, const SmVertex * cpVertex, SmVector3d & rTangent)
{
  NER( pEdge );
  NER( pEdge->GetCurve() );
  SmExtent1d sIvl = pEdge->GetInterval();
  SmVector3d sPV[2];
  SER( pEdge->GetCurve()->Evaluate( (pEdge->GetStartVertex() == cpVertex) ? sIvl.GetMin() : sIvl.GetMax(), 1, TRUE, sPV ));
  rTangent = sPV[1];
  SER( rTangent.Unitize() );
  return SM_SUCCESS;
}

// The angle (degrees) of w in the frame (X, Z) of a meridian plane: atan2(w.Z, w.X).
double sm_MeridianAngleDeg(const SmVector3d & w, const SmVector3d & crX, const SmVector3d & crZ)
{
  return SM_RAD2DEG( atan2( w.Dot( crZ ), w.Dot( crX )));
}


} // namespace

/*******************************************************************//**
PURPOSE: Constructor.
***********************************************************************/
SmBlendCorner::SmBlendCorner
  (const SmContext      & crContext,
   const SmVertex       * cpVertex,
   SmTArray<SmEdgeuse*> * pSolverEUs,
   SmTArray<SmEdgeuse*> * pAllEUs,
   SmFilletExecutive    * pExec,
   double                 dApproxTol3d,
   double                 dTangencyTolRadians,
   SmBlendCornerEvent     eEvent,
   SmBlendTurn            eTurn)
 : SmFilletCorner(crContext, cpVertex, pSolverEUs, pAllEUs, pExec, dApproxTol3d, dTangencyTolRadians),
   m_eEvent(eEvent),
   m_eTurn(eTurn),
   m_bG1Case(FALSE),
   m_bMitre(FALSE),
   m_pSideEU(NULL),
   m_pTurnSurface(NULL)
{
}

SmBlendCorner::~SmBlendCorner()
{
  DeleteSweptSurface();
}

void SmBlendCorner::DeleteSweptSurface()
{
  for (ULONG ii=0; ii<m_vTurn3d.GetSize(); ii++) { delete m_vTurn3d[ii]; }
  for (ULONG ii=0; ii<m_vTurnUV.GetSize(); ii++) { delete m_vTurnUV[ii]; }
  m_vTurn3d.ReSet();
  m_vTurnUV.ReSet();
  m_vTurnOrients.ReSet();
  if ( m_pTurnSurface != NULL ) { delete m_pTurnSurface; m_pTurnSurface = NULL; }
}

/*******************************************************************//**
PURPOSE: The corner event at a filleted vertex with two blended edgeuses of
    one convexity among three.

NOTES: A Turn when the blended edges do not continue each other there (not
    G1) and share a face. Mitre is the default; with GetSweptTurns() it is
    Swept, which the blend corner builds where:
    - the side edge is convex: at a sharp concave one the ball cannot pass
      (SM_FILERR_NO_ROOM, naming the radius that edge must be rounded at
      first, the owner's ruling of 2026-10-31);
    - the side edge is a straight line, sharp, square at the vertex to both
      blended edges (each round's end section then holds it), the common
      face is a plane, and the blends are constant rounds of one radius
      (elsewhere SM_FILERR_CORNER_UNSUPPORTED: the surface of revolution, a
      drafted wall's sweep and unequal radii' fourth side are not built yet).
***********************************************************************/
SmStatus SmBlendCorner::Classify
  (const SmVertex             * cpVertex,
   const SmTArray<SmEdgeuse*> & crSolverEUs,
   const SmTArray<SmEdgeuse*> & crAllEUs,
   SmFilletExecutive          * pExec,
   double                       dTangencyTolDeg,
   SmBlendCornerEvent         & reEvent,
   SmBlendTurn                & reTurn)
{
  reEvent = SM_BCE_NONE;
  reTurn  = SM_BT_MITRE;
  NER( cpVertex );
  NER( pExec );
  if ( crSolverEUs.GetSize() != 2 || ( crAllEUs.GetSize() != 3 && crAllEUs.GetSize() != 4 ))
    { return SM_SUCCESS; }
  SmEdgeuse * pEU1 = crSolverEUs[0];
  SmEdgeuse * pEU2 = crSolverEUs[1];
  NER( pEU1 );
  NER( pEU2 );
  // a Join: the blends continue each other, or, at a face corner of the other convexity, meet at
  // 5 degrees or less (they blend there: no mitre so thin)
  double dAngleDeg = 90.0;
  const SmBoolean bG1 = sm_TestG1Continuity( cpVertex, pEU1->GetEdge(), pEU2->GetEdge(), dTangencyTolDeg, &dAngleDeg, 0 );
  if ( bG1 || ( dAngleDeg <= 5.0 && sm_CheckFaceComplimentaryCornerConcavity( cpVertex, pEU1, pEU2, NULL, 0 )))
    {
      reEvent = SM_BCE_JOIN;
      return SM_SUCCESS;
    }
  reEvent = SM_BCE_TURN;
  if ( !pExec->GetSweptTurns() )
    { return SM_SUCCESS; }
  reTurn = SM_BT_SWEPT;

  // a Swept turn is built at one side edge between the blends' other faces
  if ( crAllEUs.GetSize() != 3 )
    {
      return sm_TurnNotBuilt( pExec, cpVertex, SM_FILERR_CORNER_UNSUPPORTED,
                              _T("Fillet Error: a swept turn is built where three edges meet (not yet at a corner of four)") );
    }
  SmEdgeuse * pRingEU = NULL;
  for (ULONG ii=0; ii<crAllEUs.GetSize(); ii++)
    {
      if ( crAllEUs[ii] != pEU1 && crAllEUs[ii] != pEU2 )
        { pRingEU = crAllEUs[ii]; }
    }
  NER( pRingEU );
  SmEdgeuse * pSideEU = sm_FindSideEdgeuse( pEU1, pRingEU, 0 );
  if ( pSideEU == NULL || sm_GetFaceByTwoEdgeuses( pEU1, pEU2, pSideEU ) == NULL )
    {
      return sm_TurnNotBuilt( pExec, cpVertex, SM_FILERR_CORNER_UNSUPPORTED,
                              _T("Fillet Error: a swept turn is built where the two rounds share a face") );
    }

  SmFilletSolver * apFS[2] = { pExec->GetFilletSolverOfEdgeuse( pEU1 ), pExec->GetFilletSolverOfEdgeuse( pEU2 ) };
  NER( apFS[0] );
  NER( apFS[1] );
  SmEdge * pSide = pSideEU->GetEdge();
  NER( pSide );
  SmBoolean bSideConcave = FALSE;
  SER( SmFilletSolver::IsConcaveEdge( *pSide, bSideConcave ));
  if ( bSideConcave )
    {
      const SmOffsetSurface * pOffset = apFS[0]->GetSurface(0);
      NER( pOffset );
      TCHAR cWhy[SM_TBLOCK_SIZE];
      smos_sprintf( cWhy, _T("Fillet Error: a swept turn cannot pass the sharp concave corner edge between the rounds; round that edge at a radius of at least %g first"),
                    smos_Fabs( pOffset->GetOffsetDistance() ));
      return sm_TurnNotBuilt( pExec, cpVertex, SM_FILERR_NO_ROOM, cWhy );
    }

  // what the exact torus needs
  const double dTol = smos_Min( apFS[0]->GetThisApproxTol3d(), apFS[1]->GetThisApproxTol3d() );
  double dRadius = -1.0;
  for (ULONG ii=0; ii<2; ii++)
    {
      SmFilletSurfaceGenerator * pFSG = apFS[ii]->GetFilletSurfaceGenerator();
      NER( pFSG );
      const SmOffsetSurface * pOff0 = apFS[ii]->GetSurface(0);
      const SmOffsetSurface * pOff1 = apFS[ii]->GetSurface(1);
      NER( pOff0 );
      NER( pOff1 );
      const double dR0 = smos_Fabs( pOff0->GetOffsetDistance() ), dR1 = smos_Fabs( pOff1->GetOffsetDistance() );
      if (   pFSG->GetFilletSurfaceGeneratorType() != SM_FSG_CIRCULAR || apFS[ii]->OffsetRadiiCanChange()
          || smos_Fabs( dR0 - dR1 ) > dTol || ( dRadius >= 0.0 && smos_Fabs( dR0 - dRadius ) > dTol ))
        {
          return sm_TurnNotBuilt( pExec, cpVertex, SM_FILERR_CORNER_UNSUPPORTED,
                                  _T("Fillet Error: a swept turn is built between constant rounds of one radius (not yet between unequal or variable rounds)") );
        }
      dRadius = dR0;
    }
  if ( sm_CheckTangentialEdgeuseAtVertex( pSideEU, cpVertex, dTangencyTolDeg ))
    {
      return sm_TurnNotBuilt( pExec, cpVertex, SM_FILERR_CORNER_UNSUPPORTED,
                              _T("Fillet Error: a swept turn is built at a sharp side edge (not at a tangent seam)") );
    }
  // each round's end section (square to its edge) holds the side edge where the edge is square
  // to it: off square by an angle whose sine is s, the section's arc of radius r leaves the edge's
  // plane by up to s r, held to the fillet's tolerance (not the tangency tolerance, which takes
  // a few degrees for tangent)
  const double dSinTol = dTol / dRadius;
  SmPoint3d  sLinePt;
  SmVector3d sAxis;
  SmFace   * pCommon = sm_GetFaceByTwoEdgeuses( pEU1, pEU2, pSideEU );
  NER( pCommon );
  NER( pCommon->GetSurface() );
  SmVector3d sT1, sT2;
  SER( sm_TangentAtVertex( pEU1->GetEdge(), cpVertex, sT1 ));
  SER( sm_TangentAtVertex( pEU2->GetEdge(), cpVertex, sT2 ));
  if (   pSide->GetCurve() == NULL || !pSide->GetCurve()->IsLine( 10, dTol, sLinePt, sAxis ) || sAxis.Unitize( FALSE ) != SM_SUCCESS
      || smos_Fabs( sAxis.Dot( sT1 )) > dSinTol || smos_Fabs( sAxis.Dot( sT2 )) > dSinTol
      || !pCommon->GetSurface()->IsPlanar( dTol ))
    {
      return sm_TurnNotBuilt( pExec, cpVertex, SM_FILERR_CORNER_UNSUPPORTED,
                              _T("Fillet Error: a swept turn is built where the side edge is straight and square to both rounded edges on a plane (not yet at a curved or drafted side edge)") );
    }
  return SM_SUCCESS;

} // end SmBlendCorner::Classify

/*******************************************************************//**
PURPOSE: The Swept turn's topology: the side vertex S (both rounds' rails on
    their other faces meet the side edge there), the plate vertices A and B
    (S's mates on each round's rail on the common face), and the corner
    edges S-A and B-S (the rounds' end sections) and A-B (where the
    corner's face meets the common face).
***********************************************************************/
SmStatus SmBlendCorner::MakeSweptTopology()
{
  NER( m_pSolverEUs );
  NER( m_pAllEUs );
  if ( !IsSwept() || m_pSolverEUs->GetSize() != 2 || m_pAllEUs->GetSize() != 3 )
    { SER( SM_ERR ); }

  SmEdgeuse * apEU[2] = { (*m_pSolverEUs)[0], (*m_pSolverEUs)[1] };
  SmEdgeuse * pRingEU = NULL;
  for (ULONG ii=0; ii<m_pAllEUs->GetSize(); ii++)
    {
      if ( (*m_pAllEUs)[ii] != apEU[0] && (*m_pAllEUs)[ii] != apEU[1] )
        { pRingEU = (*m_pAllEUs)[ii]; }
    }
  NER( pRingEU );
  m_pSideEU = sm_FindSideEdgeuse( apEU[0], pRingEU, 0 );
  NER( m_pSideEU );
  SmFace * pCommon = sm_GetFaceByTwoEdgeuses( apEU[0], apEU[1], m_pSideEU );
  NER( pCommon );

  SmFilletSolver * apFS[2];
  SmFilletGeom   * apFG[2];
  ULONG            alSide[2], alFace[2];
  SmEdgeuse      * apSideEU[2] = { m_pSideEU, m_pSideEU->GetRadial() };
  for (ULONG ii=0; ii<2; ii++)
    {
      apFS[ii] = GetFilletExecutive()->GetFilletSolverOfEdgeuse( apEU[ii] );
      NER( apFS[ii] );
      apFG[ii] = apFS[ii]->GetFirstFilletGeom();
      NER( apFG[ii] );
      alSide[ii] = apFS[ii]->FindIndexOfRailXSideEdgeuse( apSideEU[ii] );
      alFace[ii] = apFS[ii]->FindIndexOfRailOnFace( pCommon );
      if ( alSide[ii] == alFace[ii] )
        { SER( SM_ERR ); }
    }

  // S, on the side edge, at both rounds' rails on their other faces
  SmFilletEdgeuse * pNewEU = NULL;
  SmFilletVertex  * pS     = new (m_pPseudoBrep) SmFilletVertex( this );
  NER( pS );
  pS->SetPointClass( SM_PC_EDGEUSE, m_pSideEU );
  pS->SetFilletVertexType( SM_FV_RAIL_X_EDGEUSE );
  m_vVertices.Add( pS );
  for (ULONG ii=0; ii<2; ii++)
    { SER( apFG[ii]->MakeRailEdgeuse( m_pPseudoBrep, alSide[ii], pS, pNewEU )); }

  // A and B, on each round's rail on the common face, straight across its end section from S:
  // S's mates, which SmFilletVertex::CalcRailIntEdgeuse places with S, each on its own round
  SmFilletVertex * apPlate[2];
  for (ULONG ii=0; ii<2; ii++)
    {
      apPlate[ii] = new (m_pPseudoBrep) SmFilletVertex( this );
      NER( apPlate[ii] );
      apPlate[ii]->SetPointClass( SM_PC_FACE, pCommon );
      apPlate[ii]->SetFilletVertexType( SM_FV_MATE );
      apPlate[ii]->SetMate( 0, pS );
      pS->SetMate( ii, apPlate[ii] );
      m_vVertices.Add( apPlate[ii] );
      SER( apFG[ii]->MakeRailEdgeuse( m_pPseudoBrep, alFace[ii], apPlate[ii], pNewEU ));
    }

  // each round ends on its own section; the corner's face is a patch of its own
  m_bBlending = TRUE;
  m_bBevel    = FALSE;

  // the corner's edges in its face's loop order: S-A, A-B, B-S
  SmFilletVertex * apFrom[3] = { pS, apPlate[0], apPlate[1] };
  SmFilletVertex * apTo[3]   = { apPlate[0], apPlate[1], pS };
  SmFilletGeom   * apOwn[3]  = { apFG[0], NULL, apFG[1] };
  for (ULONG ii=0; ii<3; ii++)
    {
      SmFilletEdge * pEdge = NULL;
      SER( SmFilletEdge::MakeFilletEdge( m_pPseudoBrep, apFrom[ii], apTo[ii], pEdge, this ));
      NER( pEdge );
      SmFilletEdgeuse * pPrimEU = (SmFilletEdgeuse*)pEdge->GetPrimaryEdgeuse();
      NER( pPrimEU );
      pPrimEU->SetFilletGeom( apOwn[ii] );
      pEdge->SetFilletEdgeType( apOwn[ii] != NULL ? SM_FE_CROSS_SECTION : SM_FE_PRECOMPUTED );
      if ( apOwn[ii] == NULL )
        { pEdge->SetOriginalFace( pCommon ); }
      m_vEdges.Add( pEdge );
    }
  return SM_SUCCESS;

} // end SmBlendCorner::MakeSweptTopology

/*******************************************************************//**
PURPOSE: The torus of the Swept turn and its loop, from the processed
    vertices S, A, B and the balls of the two rounds at S (see the file's
    head). Its loop is the natural boundary of its patch, the pole at S
    left out: the two end sections and the arc A-B.
***********************************************************************/
SmStatus SmBlendCorner::BuildSweptSurface()
{
  DeleteSweptSurface();
  if ( m_vVertices.GetSize() != 3 )
    { SER( SM_ERR ); }
  SmFilletVertex * pS = m_vVertices[0];
  SmFilletVertex * apPlate[2] = { m_vVertices[1], m_vVertices[2] };
  NER( pS );
  NER( apPlate[0] );
  NER( apPlate[1] );
  if ( !pS->IsProcessed() || !apPlate[0]->IsProcessed() || !apPlate[1]->IsProcessed() )
    { SER( SM_ERR ); }
  const double dTol = GetThisApproxTol3d();

  // the balls at S: each round's offset surfaces at its section through S
  SmPoint3d sCentre[2];
  double    dMinor = -1.0;
  for (ULONG ii=0; ii<2; ii++)
    {
      SmTArray<SmFilletGeom*> sGeoms;
      apPlate[ii]->GetFilletGeoms( sGeoms );
      if ( sGeoms.GetSize() != 1 )
        { SER( SM_ERR ); }
      SmFilletVertexuse * pVU = pS->GetVUAtRailEnd( sGeoms[0] );
      NER( pVU );
      SmTsectPnt & rTsect = pVU->GetTsectPnt();
      SmPoint3d sOn[2];
      for (ULONG k=0; k<2; k++)
        {
          SmOffsetSurface * pOffset = sGeoms[0]->GetOffsetSurface( k );
          NER( pOffset );
          SER( pOffset->EvaluatePoint( rTsect.UVPos( k ), sOn[k] ));
        }
      if ( sOn[0].DistanceBetween( sOn[1] ) > dTol )
        { SER( SM_ERR ); } // the section's offset points are its ball's centre
      sCentre[ii] = sOn[0];
      const double dR = sCentre[ii].DistanceBetween( pS->GetPoint() );
      if ( dMinor >= 0.0 && smos_Fabs( dR - dMinor ) > dTol )
        { SER( SM_ERR ); } // Classify took the radii to be equal
      dMinor = dR;
    }

  // the axis: the side edge's line; the torus's centre the balls' common foot on it
  SmPoint3d  sLinePt;
  SmVector3d sAxis;
  NER( m_pSideEU );
  SmEdge * pSide = m_pSideEU->GetEdge();
  NER( pSide );
  NER( pSide->GetCurve() );
  if ( !pSide->GetCurve()->IsLine( 10, dTol, sLinePt, sAxis ))
    { SER( SM_ERR ); }
  SER( sAxis.Unitize() );
  SmPoint3d sFoot[2];
  for (ULONG ii=0; ii<2; ii++)
    { sFoot[ii] = sLinePt + sAxis * SmVector3d( sCentre[ii] - sLinePt ).Dot( sAxis ); }
  if ( sFoot[0].DistanceBetween( sFoot[1] ) > dTol )
    { SER( SM_ERR ); }
  const SmPoint3d sOrigin = sFoot[0];
  SmVector3d sX  = sCentre[0] - sOrigin;
  SmVector3d sX2 = sCentre[1] - sOrigin;
  double dMajor = sX.Length();
  if ( dMajor <= dTol || smos_Fabs( sX2.Length() - dMajor ) > dTol )
    { SER( SM_ERR ); }
  // S on the axis: the horn torus, its tube through the axis at S (as a revolved circle that
  // reaches the axis, SmRotationalSweepGeometry::AnalyticSurface). The two radii, computed apart,
  // differ in their last bits, and a spindle torus's tube would cross the axis there.
  if ( smos_Fabs( dMajor - dMinor ) <= dTol && SmVector3d( pS->GetPoint() - sOrigin ).Length() <= dTol )
    { dMajor = dMinor; }
  SER( sX.Unitize() );
  SER( sX2.Unitize() );
  SmVector3d sZ = sX * sX2;
  const double dSin = sZ.Length();
  if ( dSin <= SM_EFF_ZERO )
    { SER( SM_ERR ); } // no turn
  SER( sZ.Unitize() );
  const SmVector3d sY = sZ * sX;
  const double dTurnDeg = SM_RAD2DEG( atan2( dSin, sX.Dot( sX2 )));

  // the end section's arc: from S to A in the meridian (X, Z) through the first ball
  const double dVS = sm_MeridianAngleDeg( pS->GetPoint() - sCentre[0], sX, sZ );
  double       dVA = sm_MeridianAngleDeg( apPlate[0]->GetPoint() - sCentre[0], sX, sZ );
  while ( dVA - dVS >  180.0 ) { dVA -= 360.0; }
  while ( dVA - dVS < -180.0 ) { dVA += 360.0; }
  SmExtent2d sDomain;
  SER( sDomain.SetMinMax( 0.0, smos_Min( dVS, dVA ), dTurnDeg, smos_Max( dVS, dVA )));
  SmAxis2Placement sFrame;
  SER( sFrame.SetCanonical( sOrigin, sX, sY ));
  // The face's normal is its surface's (SmBrep::OrientTrimmedSurfaces turns a face the other way
  // by swapping its surface's parameters): out of the material, towards the ball where the rounds
  // add material (concave), away from it where they remove it. The torus's normal points out of
  // its tube unless its NURB parameters are swapped.
  SmBoolean bAddsMaterial = FALSE;
  NER( m_pSolverEUs );
  NER( (*m_pSolverEUs)[0] );
  SER( SmFilletSolver::IsConcaveEdge( *(*m_pSolverEUs)[0]->GetEdge(), bAddsMaterial ));
  SmTorus * pTorus = NULL;
  SER( SmTorus::CreateCanonical( m_crContext, sFrame, dMajor, dMinor, pTorus, bAddsMaterial, FALSE, NULL, &sDomain ));
  NER( pTorus );
  m_pTurnSurface = pTorus;

  // its corners are the corner's vertices: A and S at u = 0, B at u = theta
  SmPoint3d sAt;
  SER( pTorus->EvaluateSTEPPoint( SmPoint2d( 0.0, dVA ), sAt ));
  if ( sAt.DistanceBetween( apPlate[0]->GetPoint() ) > dTol )
    { SER( SM_ERR ); }
  SER( pTorus->EvaluateSTEPPoint( SmPoint2d( dTurnDeg, dVA ), sAt ));
  if ( sAt.DistanceBetween( apPlate[1]->GetPoint() ) > dTol )
    { SER( SM_ERR ); } // the second round's end section is the torus's at theta
  SER( pTorus->EvaluateSTEPPoint( SmPoint2d( 0.0, dVS ), sAt ));
  if ( sAt.DistanceBetween( pS->GetPoint() ) > dTol )
    { SER( SM_ERR ); }

  SER( pTorus->CreateNaturalUVTrimCurves( m_crContext, pTorus->GetNaturalUVDomain(), SM_SP_V, TRUE,
                                          m_vTurn3d, m_vTurnUV, m_vTurnOrients ));
  if ( m_vTurn3d.GetSize() != 3 || m_vTurnUV.GetSize() != 3 || m_vTurnOrients.GetSize() != 3 )
    { SER( SM_ERR ); } // the ends and the arc; the side at S a point
  return SM_SUCCESS;

} // end SmBlendCorner::BuildSweptSurface

/*******************************************************************//**
PURPOSE: The corner's edges: the arc A-B is the torus's (its iso-curve at
    A's tube angle, from A to B); the end sections are the rounds' own
    cross-sections (SmFilletEdge::CalcCrossSection).
***********************************************************************/
SmStatus SmBlendCorner::CalcCornerEdgeGeom()
{
  if ( !IsSwept() )
    { return SmFilletCorner::CalcCornerEdgeGeom(); }
  SER( BuildSweptSurface() );
  if ( m_vEdges.GetSize() != 3 )
    { SER( SM_ERR ); }
  SmFilletEdge * pArcEdge = m_vEdges[1];
  NER( pArcEdge );
  const SmPoint3d sA = pArcEdge->GetStartVertex()->GetPoint();
  const SmPoint3d sB = pArcEdge->GetOtherVertex( pArcEdge->GetStartVertex() )->GetPoint();
  const double dTol = GetThisApproxTol3d();
  SmCurve * pArc = NULL;
  for (ULONG ii=0; ii<m_vTurn3d.GetSize() && pArc == NULL; ii++)
    {
      SmCurve * pCurve = m_vTurn3d[ii];
      NER( pCurve );
      SmExtent1d sIvl = pCurve->GetNaturalInterval();
      SmPoint3d  sLo, sHi;
      SER( pCurve->EvaluatePoint( sIvl.GetMin(), sLo ));
      SER( pCurve->EvaluatePoint( sIvl.GetMax(), sHi ));
      const SmBoolean bForward  = ( sLo.DistanceBetween( sA ) <= dTol && sHi.DistanceBetween( sB ) <= dTol ) ? TRUE : FALSE;
      const SmBoolean bBackward = ( sLo.DistanceBetween( sB ) <= dTol && sHi.DistanceBetween( sA ) <= dTol ) ? TRUE : FALSE;
      if ( !bForward && !bBackward )
        { continue; }
      SER( pCurve->Copy( m_crContext, pArc ));
      NER( pArc );
      if ( bBackward )
        {
          SmExtent1d sNew;
          SmStatus eStatus = pArc->ReverseParameterization( sIvl, sNew );
          if ( eStatus != SM_SUCCESS )
            { delete pArc; SER( eStatus ); }
        }
    }
  NER( pArc );
  pArcEdge->SetCurve( pArc, TRUE );
  pArc->SetOwner( pArcEdge );
  return SmFilletCorner::CalcCornerEdgeGeom();

} // end SmBlendCorner::CalcCornerEdgeGeom

/*******************************************************************//**
PURPOSE: The corner's face: the torus bounded by its natural loop. Its
    edges and the rounds' and the common face's meet in the stitch and the
    merge.
***********************************************************************/
SmStatus SmBlendCorner::CreateSweptPatch()
{
  NER( m_pTurnSurface );
  SmBrep * pFilletBrep = m_pExecutive->GetFilletBrep();
  NER( pFilletBrep );
  SmTArray<ULONG> sCurveLoops;
  sCurveLoops.Add( m_vTurn3d.GetSize() );
  SmTArray<SmPoint3d> sLoopPoints;
  SmRegion * pNewRegion = NULL;
  SmShell  * pNewShell  = NULL;
  SmFace   * pNewFace   = NULL;
  SmTorus  * pSurface   = m_pTurnSurface;
  // the face takes the surface and the curves
  SmTArray<SmCurve*>        s3d( m_vTurn3d );
  SmTArray<SmBSplineCurve*> sUV( m_vTurnUV );
  m_pTurnSurface = NULL;
  m_vTurn3d.ReSet();
  m_vTurnUV.ReSet();
  SER( pFilletBrep->MakeFaceWithCurves( NULL, sCurveLoops, &s3d, &sUV, m_vTurnOrients, sLoopPoints, pSurface,
                                        pSurface->GetNaturalUVDomain(), SM_OT_SAME, pNewRegion, pNewShell, pNewFace ));
  NER( pNewFace );
  m_vSurfaces.Add( pSurface ); // the corner's patch face (SmFilletCorner::GetPatchFaces, the stitch)
  return SM_SUCCESS;

} // end SmBlendCorner::CreateSweptPatch

/*******************************************************************//**
PURPOSE: The Mitre turn's topology: the side vertex S, where both rounds'
    rails on their other faces meet the side edge, the vertex X where their
    rails on the common face cross, and the mitre X-S between the two rounds,
    each round traced on past the vertex (a bevel corner) to meet the other.
***********************************************************************/
SmStatus SmBlendCorner::MakeMitreTopology()
{
  NER( m_pSolverEUs );
  NER( m_pAllEUs );
  if ( !IsMitreTurn() || m_pSolverEUs->GetSize() != 2 )
    { SER( SM_ERR ); }
  SmEdgeuse * apEU[2] = { (*m_pSolverEUs)[0], (*m_pSolverEUs)[1] };
  NER( apEU[0] );
  NER( apEU[1] );
  SmEdgeuse * pRingEU = NULL;
  for (ULONG ii=0; ii<m_pAllEUs->GetSize(); ii++)
    {
      if ( (*m_pAllEUs)[ii] != apEU[0] && (*m_pAllEUs)[ii] != apEU[1] )
        { pRingEU = (*m_pAllEUs)[ii]; }
    }
  NER( pRingEU );
  SmEdgeuse * pSideEU = sm_FindSideEdgeuse( apEU[0], pRingEU, 0 );
  NER( pSideEU );
  SmFace * pCommon = sm_GetFaceByTwoEdgeuses( apEU[0], apEU[1], pSideEU );
  NER( pCommon );

  // the rounds' face corner on the common face: concave where the rounds meet beyond the side
  // edge's faces (a pocket's rim, a boss's base), the mitre at the face's concave corner
  m_bMitre = sm_CheckFaceComplimentaryCornerConcavity( m_cpVertex, apEU[0], apEU[1], NULL, 0 );
  m_bBevel    = TRUE;
  m_bBlending = FALSE;

  SmFilletGeom * apFG[2];
  ULONG          alSide[2], alFace[2];
  SmEdgeuse    * apSideEU[2] = { pSideEU, pSideEU->GetRadial() };
  for (ULONG ii=0; ii<2; ii++)
    {
      SmFilletSolver * pFS = GetFilletExecutive()->GetFilletSolverOfEdgeuse( apEU[ii] );
      NER( pFS );
      apFG[ii] = pFS->GetFirstFilletGeom();
      NER( apFG[ii] );
      alSide[ii] = pFS->FindIndexOfRailXSideEdgeuse( apSideEU[ii] );
      alFace[ii] = pFS->FindIndexOfRailOnFace( pCommon );
      if ( alSide[ii] == alFace[ii] )
        { SER( SM_ERR ); }
    }

  SmFilletEdgeuse * pNewEU = NULL;
  SmFilletVertex  * pS     = new (m_pPseudoBrep) SmFilletVertex( this );
  NER( pS );
  pS->SetPointClass( SM_PC_EDGEUSE, pSideEU );
  pS->SetFilletVertexType( SM_FV_RAIL_X_EDGEUSE );
  m_vVertices.Add( pS );
  SmFilletVertex * pX = new (m_pPseudoBrep) SmFilletVertex( this );
  NER( pX );
  pX->SetPointClass( SM_PC_FACE, pCommon );
  pX->SetFilletVertexType( SM_FV_RAIL_X_RAIL );
  m_vVertices.Add( pX );
  for (ULONG ii=0; ii<2; ii++)
    {
      SER( apFG[ii]->MakeRailEdgeuse( m_pPseudoBrep, alSide[ii], pS, pNewEU ));
      SER( apFG[ii]->MakeRailEdgeuse( m_pPseudoBrep, alFace[ii], pX, pNewEU ));
    }

  SmFilletEdge * pMitre = NULL;
  SER( SmFilletEdge::MakeFilletEdge( m_pPseudoBrep, pS, pX, pMitre, this ));
  NER( pMitre );
  SmFilletEdgeuse * pPrimEU = (SmFilletEdgeuse*)pMitre->GetPrimaryEdgeuse();
  NER( pPrimEU );
  SmFilletEdgeuse * pMateEU = (SmFilletEdgeuse*)pPrimEU->GetMate();
  NER( pMateEU );
  pPrimEU->SetFilletGeom( apFG[0] );
  pMateEU->SetFilletGeom( apFG[1] );
  pMitre->SetFilletEdgeType( SM_FE_FILLET_X_FILLET );
  m_vEdges.Add( pMitre );
  return SM_SUCCESS;

} // end SmBlendCorner::MakeMitreTopology

/*******************************************************************//**
PURPOSE: Where the two rounds' rails meet the side edge at two places, the
    side vertex splits (SplitMitreSideVertex): the larger round (the one
    meeting the edge farther from the vertex) meets the smaller one's side
    face between them.
***********************************************************************/
SmStatus SmBlendCorner::AdjustMitreTopology()
{
  if ( m_vVertices.GetSize() != 2 || m_vEdges.GetSize() != 1 )
    { SER( SM_ERR ); }
  SmFilletVertex * pS = m_vVertices[0];
  NER( pS );
  SmTArray<SmFilletGeom*> sGeoms;
  pS->GetFilletGeoms( sGeoms );
  if ( sGeoms.GetSize() != 2 )
    { SER( SM_ERR ); }
  SmEdgeuse         * pSideEU = (SmEdgeuse*)pS->GetPointClassObject();
  NER( pSideEU );
  SmFilletVertexuse * pVU1    = pS->GetVUAtRailEnd( sGeoms[0] );
  NER( pVU1 );
  SmFilletVertexuse * pVU2    = pS->GetVUAtRailEnd( sGeoms[1], pVU1 );
  NER( pVU2 );
  SmCurve * pSideCurve = pSideEU->GetEdge()->GetCurve();
  NER( pSideCurve );
  SmPoint3d sPnt1, sPnt2;
  SER( pSideCurve->EvaluatePoint( pVU1->GetTsectPnt().m_dCurveParameter, sPnt1 ));
  SER( pSideCurve->EvaluatePoint( pVU2->GetTsectPnt().m_dCurveParameter, sPnt2 ));
  if ( sPnt1.DistanceBetween( sPnt2 ) < m_dThisApproxTol3d )
    { return SM_SUCCESS; } // both rails meet the side edge at S
  // at a concave face corner the larger round runs on past the smaller one's side face, onto
  // its continuation
  return SplitMitreSideVertex( pS, m_bMitre );

} // end SmBlendCorner::AdjustMitreTopology

/*******************************************************************//**
PURPOSE: The side vertex S split where the two rounds' rails meet the side
    edge at two places: L, the larger round's (farther from the vertex), and
    the smaller round's, Sm.

NOTES: The larger round meets the smaller one's side face F along the curve
   from L down to where it meets the mitre (P):
   - at a convex face corner on F itself: P is where the smaller round's
     rail meets the edge, and the mitre runs P-X (edges L-P, P-X);
   - at a concave face corner on F's continuation past the edge: P is where
     the mitre meets that continuation, and the gap P-L-Sm between the side
     edge, the larger round and the smaller round's rail run on (P-Sm) is a
     patch of F's extended surface (edges X-P, P-L, L-Sm, Sm-P).
***********************************************************************/
SmStatus SmBlendCorner::SplitMitreSideVertex
  (SmFilletVertex * pS,                 // in : the side vertex
   SmBoolean        bGapOnContinuation) // in : the larger round meets the side face's continuation
{
  NER( pS );
  SmTArray<SmFilletGeom*> sGeoms;
  pS->GetFilletGeoms( sGeoms );
  if ( sGeoms.GetSize() != 2 )
    { SER( SM_ERR ); }
  SmEdgeuse         * pSideEU = (SmEdgeuse*)pS->GetPointClassObject();
  NER( pSideEU );
  SmFilletVertexuse * apVU[2];
  apVU[0] = pS->GetVUAtRailEnd( sGeoms[0] );
  NER( apVU[0] );
  apVU[1] = pS->GetVUAtRailEnd( sGeoms[1], apVU[0] );
  NER( apVU[1] );
  SmCurve * pSideCurve = pSideEU->GetEdge()->GetCurve();
  NER( pSideCurve );
  double    adParam[2];
  SmPoint3d asPnt[2];
  for (ULONG ii=0; ii<2; ii++)
    {
      adParam[ii] = apVU[ii]->GetTsectPnt().m_dCurveParameter;
      SER( pSideCurve->EvaluatePoint( adParam[ii], asPnt[ii] ));
    }
  // the larger round: the one meeting the edge farther from the vertex. The smaller one's side face
  // is that of its side edgeuse: pSideEU's for the first round, its radial's for the second.
  const SmPoint3d sVertPnt = m_cpVertex->GetPoint();
  const ULONG     lLarge   = ( asPnt[0].DistanceBetween( sVertPnt ) < asPnt[1].DistanceBetween( sVertPnt )) ? 1 : 0;
  const ULONG     lSmall   = 1 - lLarge;
  SmFace        * pSideFace = ( lSmall == 0 ) ? pSideEU->GetFace() : pSideEU->GetRadial()->GetFace();
  NER( pSideFace );
  SmFilletVertex * pX = m_vVertices[1];
  NER( pX );
  SmFilletEdge * pMitre = m_vEdges[0];
  NER( pMitre );
  m_bTopologyAdjusted = TRUE;

  if ( !bGapOnContinuation )
    {
      // S becomes L; the smaller round's rail end moves to P, on F where its rail meets the edge
      pS->SetPoint( asPnt[lLarge] );
      pS->SetOriginalTParam( adParam[lLarge] );
      pS->SetStatus( SM_FIL_PROCESSED );
      SmFilletVertex * pP = new (m_pPseudoBrep) SmFilletVertex( this );
      NER( pP );
      SmObjDelete sCleanP( pP );
      pP->SetPoint( asPnt[lSmall] );
      pP->SetPointClass( SM_PC_FACE, pSideFace );
      pP->SetFilletVertexType( SM_FV_FILLET_X_FILLET );
      SER( pS->Remove( apVU[lSmall] ));
      SER( pP->PostInsert( apVU[lSmall] ));
      m_vVertices.InsertAt( 1, pP );
      sCleanP.Clear();

      SmFilletEdge * pOnFace = NULL;
      SER( SmFilletEdge::MakeFilletEdge( m_pPseudoBrep, pS, pP, pOnFace, this ));
      NER( pOnFace );
      pOnFace->SetFilletEdgeType( SM_FE_FILLET_X_SIDE_FACE );
      pOnFace->SetOriginalFace( pSideFace );
      ((SmFilletEdgeuse*)pOnFace->GetPrimaryEdgeuse())->SetFilletGeom( sGeoms[lLarge] );

      SmFilletEdge * pNewMitre = NULL;
      SER( SmFilletEdge::MakeFilletEdge( m_pPseudoBrep, pP, pX, pNewMitre, this ));
      NER( pNewMitre );
      pNewMitre->SetFilletEdgeType( SM_FE_FILLET_X_FILLET );
      SmFilletEdgeuse * pPrimEU = (SmFilletEdgeuse*)pNewMitre->GetPrimaryEdgeuse();
      pPrimEU->SetFilletGeom( sGeoms[lLarge] );
      ((SmFilletEdgeuse*)pPrimEU->GetMate())->SetFilletGeom( sGeoms[lSmall] );

      m_vEdges.RemoveAt( 0 );
      delete pMitre;
      m_vEdges.Add( pOnFace );
      m_vEdges.Add( pNewMitre );
      return SM_SUCCESS;
    }

  // S becomes P, where the mitre meets F's continuation (found with the mitre); the rails' ends on
  // the edge are L and Sm
  SmFilletVertex * apEnd[2];
  for (ULONG ii=0; ii<2; ii++)
    {
      apEnd[ii] = new (m_pPseudoBrep) SmFilletVertex( this );
      NER( apEnd[ii] );
      apEnd[ii]->SetFilletVertexType( SM_FV_RAIL_X_EDGEUSE );
      m_vVertices.Add( apEnd[ii] );
    }
  SmFilletVertex * pL = apEnd[0], * pSm = apEnd[1];
  pL->SetPoint( asPnt[lLarge] );
  pL->SetOriginalTParam( adParam[lLarge] );
  pL->SetStatus( SM_FIL_PROCESSED );
  pSm->SetPoint( asPnt[lSmall] );
  pSm->SetOriginalTParam( adParam[lSmall] );
  pSm->SetStatus( SM_FIL_PROCESSED );
  pS->SetFilletVertexType( SM_FV_FILLET_X_FILLET );
  pS->SetStatus( SM_FIL_UNPROCESSED );
  SER( pS->Remove( apVU[lLarge] ));
  SER( pL->PostInsert( apVU[lLarge] ));
  SER( pS->Remove( apVU[lSmall] ));
  SER( pSm->PostInsert( apVU[lSmall] ));

  // F's continuation, where the larger round meets it and the gap's patch lies
  SmSurface * pExtended = NULL;
  NER( pSideFace->GetSurface() );
  SER( pSideFace->GetSurface()->CreateExtendedSurface( m_crContext, 10.0 * asPnt[0].DistanceBetween( asPnt[1] ), SM_CT_G1, pExtended ));
  NER( pExtended );
  SmObjDelete sCleanExtended( pExtended );
  SmSurface * pAnalytic = NULL;
  SER( pExtended->CopyAndAddAnalytics( m_crContext, pAnalytic ));
  AddExtendedOriginalSurface( pSideFace, pAnalytic );

  SmFilletEdge * pPL = NULL, * pLSm = NULL, * pSmP = NULL;
  SER( SmFilletEdge::MakeFilletEdge( m_pPseudoBrep, pS, pL, pPL, this ));
  NER( pPL );
  pPL->SetFilletEdgeType( SM_FE_FILLET_X_SIDE_FACE );
  pPL->SetOriginalFace( pSideFace );
  ((SmFilletEdgeuse*)pPL->GetPrimaryEdgeuse())->SetFilletGeom( sGeoms[lLarge] );
  SER( SmFilletEdge::MakeFilletEdge( m_pPseudoBrep, pL, pSm, pLSm, this ));
  NER( pLSm );
  pLSm->SetFilletEdgeType( SM_FE_ON_EDGE );
  pLSm->SetOriginalEdge( pSideEU->GetEdge() );
  pLSm->SetOriginalFace( pSideFace );
  SER( SmFilletEdge::MakeFilletEdge( m_pPseudoBrep, pSm, pS, pSmP, this ));
  NER( pSmP );
  pSmP->SetFilletEdgeType( SM_FE_RAIL_EXTENSION );
  pSmP->SetOriginalFace( pSideFace );
  ((SmFilletEdgeuse*)pSmP->GetPrimaryEdgeuse())->SetFilletGeom( sGeoms[lSmall] );
  m_vEdges.Add( pPL );
  m_vEdges.Add( pLSm );
  m_vEdges.Add( pSmP );
  m_bBevelGapPatch = TRUE;
  return SM_SUCCESS;

} // end SmBlendCorner::SplitMitreSideVertex

/*******************************************************************//**
PURPOSE: The Mitre turn's patch: none where the rounds meet on the side
    edge's faces; where they meet on a side face's continuation, the gap
    there (the edges after the mitre) as a patch of that continuation.
***********************************************************************/
SmStatus SmBlendCorner::CreateMitrePatch()
{
  if ( !m_bBevelGapPatch )
    { return SM_SUCCESS; }
  if ( m_vEdges.GetSize() != 4 )
    { SER( SM_ERR ); }
  SmTArray<SmFilletEdge*> sAll;
  sAll.Append( m_vEdges );
  m_vEdges.RemoveAt( 0 );
  SmStatus eGap = SmFilletCorner::CreateEndPatch();
  m_vEdges.ReSet();
  m_vEdges.Append( sAll );
  SER( eGap );
  return SM_SUCCESS;

} // end SmBlendCorner::CreateMitrePatch

/*******************************************************************//**
PURPOSE: The executive's corner phases, by the corner's event: the Swept
    turn's own construction, or the two blends' (Mitre turn, Join).
***********************************************************************/
SmStatus SmBlendCorner::MakeCornerTopology()
{
  return IsSwept() ? MakeSweptTopology() : IsMitreTurn() ? MakeMitreTopology() : MakeTwoBlendTopology();
}

SmStatus SmBlendCorner::AdjustTopology()
{
  return IsSwept() ? SM_SUCCESS : IsMitreTurn() ? AdjustMitreTopology() : AdjustTwoBlendTopology();
}

SmStatus SmBlendCorner::AdjustTangentCorner()
{
  return ( IsSwept() || IsMitreTurn() ) ? SM_SUCCESS : AdjustTwoBlendTangentCorner();
}

SmStatus SmBlendCorner::CreateCornerPatch()
{
  return IsSwept() ? CreateSweptPatch() : IsMitreTurn() ? CreateMitrePatch() : CreateTwoBlendPatch();
}

/********************************************************************
PURPOSE: Compute the angle between two given edgeuses

NOTES:
********************************************************************/
static double sm_FindAngleBetweenEdgeuses
 (const SmVertex * cpCornerVert,
  SmEdgeuse      * pEU1,
  SmEdgeuse      * pEU2,
  int              iDebugLevel )
{
#ifdef SM_DEBUG_CODE
  if ( iDebugLevel ) 
    {
      smgfx_SetLook(1,2, 1,0,0); pEU1->Draw(); sm_GraphicsLoop();
      pEU2->Draw(); sm_GraphicsLoop();
      sm_GraphicsLoop();
    }
#else
    SM_REF1(iDebugLevel);
#endif // SM_DEBUG_CODE

  SmBoolean bAtEdgeStart1 = TRUE;
  SmVertex * pVert1 = pEU1->GetVertexuse()->GetVertex();
  if(   (cpCornerVert == pVert1 && pEU1->GetOrientation() == SM_OT_OPPOSITE)
     || (cpCornerVert != pVert1 && pEU1->GetOrientation() == SM_OT_SAME)) 
    {
      bAtEdgeStart1 = FALSE;
    }

  SmVector3d sTangent1;
  SmEdge   * pE1 = pEU1->GetEdge();

  SER(pE1->GetEndTangent(bAtEdgeStart1,sTangent1));
  if (!bAtEdgeStart1) 
    {
      sTangent1 = -sTangent1;
    }

  SmBoolean  bAtEdgeStart2 = TRUE;
  SmVertex * pVert2        = pEU2->GetVertexuse()->GetVertex();

  if (   (cpCornerVert == pVert2 && pEU2->GetOrientation() == SM_OT_OPPOSITE)
      || (cpCornerVert != pVert2 && pEU2->GetOrientation() == SM_OT_SAME)) 
    {
      bAtEdgeStart2 = FALSE;
    }

  SmVector3d sTangent2;
  SmEdge * pE2 = pEU2->GetEdge();
  SER(pE2->GetEndTangent(bAtEdgeStart2,sTangent2));
  if (!bAtEdgeStart2) 
    {
      sTangent2 = -sTangent2;
    }

  double dAngle;
  SE(sTangent1.AngleBetween(sTangent2,dAngle));

  return dAngle;

} // end sm_FindAngleBetweenEdgeuses

/*******************************************************************//**
PURPOSE: SolidCraft change. The corner's two filleted edgeuses, outward from
    its vertex.
***********************************************************************/
SmStatus SmBlendCorner::GetFilletedEdgeuses
  (SmEdgeuse *& rpFirst,    // out: the first filleted edgeuse
   SmEdgeuse *& rpSecond)   // out: the second
  const
{
  rpFirst = rpSecond = NULL;
  NER(m_pSolverEUs);
  if ( m_pSolverEUs->GetSize() != 2 )
    { SER( SM_ERR ); }
  rpFirst  = (*m_pSolverEUs)[0];
  rpSecond = (*m_pSolverEUs)[1];
  NER(rpFirst);
  NER(rpSecond);
  return SM_SUCCESS;
} // end SmBlendCorner::GetFilletedEdgeuses

/*******************************************************************//**
PURPOSE: Make corner topology which consists of an array of
    SmFilletVertex and an array of SmFilletEdge.

NOTES:
***********************************************************************/
SmStatus SmBlendCorner::MakeTwoBlendTopology()
{
int iDebugLevel = 0;
#ifdef SM_DEBUG_CODE
iDebugLevel = DebugLevel();
#endif // SM_DEBUG_CODE

  // locals
  ULONG ii;
  ULONG lTotalEUs = m_pAllEUs->GetSize();

  // no work - unhandled case - Currently we handle only 3x2 and 4x2
  if ( lTotalEUs > 4)
    { SER(SM_ERR) ; } 
     
  // Nx2 corner, get the two filleted edgeuses from the base brep
  SmEdgeuse * pBaseEdgeuse1 = (*m_pSolverEUs)[0];
  SmEdgeuse * pBaseEdgeuse2 = (*m_pSolverEUs)[1];

  // for 4 EU case
  if (lTotalEUs == 4)
    {
#ifdef SM_DEBUG_CODE
      if ( iDebugLevel > 0 ) 
        {
          SmBrep *pBrep = GetFilletExecutive()->GetTargetBrep() ;
          smgfx_SetLook(1,2, 0,0,1); if(pBrep) pBrep->Draw(TRUE) ; sm_GraphicsLoop() ;
          smgfx_SetLook(3,4, 1,0,0); pBaseEdgeuse1->Draw(); sm_GraphicsLoop();
          smgfx_SetLook(3,4, 1,0,1); pBaseEdgeuse2->Draw(); sm_GraphicsLoop();
          sm_GraphicsLoop();
        }
#endif // SM_DEBUG_CODE

      // check state: to be filleted edgeuses are on list of outward edgeuses surrounding the corner
      ULONG lFoundIndex1 = 0, lFoundIndex2 = 0;
      if (   !m_pAllEUs->FindElement( pBaseEdgeuse1, lFoundIndex1 )
          || !m_pAllEUs->FindElement( pBaseEdgeuse2, lFoundIndex2 ) )
        { SER(SM_ERR) ; }

      // Do a little test to make sure two filleted edgeuses
      // are not adjacent to each other
      if ( (lFoundIndex1+lFoundIndex2)%2 != 0 )
        {
          // 2 filleted edgeuses appeared next to each other
          // in the all-edgeuses array
          SER(SM_ERR);
        }

      // Check the non-filleted (side) edges for mixed convexity
      for (ii=0; ii<lTotalEUs; ii++)
        {
          SmEdgeuse * pSideEU = (*m_pAllEUs)[ii];

          // skip edgeuses to be filleted - looking for side (not filleted edgeuses)
          if ( pSideEU == pBaseEdgeuse1 || pSideEU == pBaseEdgeuse2 )
            { continue ; } // (not a side edge)

          // Classify a face's corner at a given vertex as convex or concave
          //            FALSE = inside face angle runs from 0 to 180+tol   degrees
          //            TRUE  = inside face angle runs from 180+tol to 360 degrees
          //  
          //         |\\ - face    \\\\\|  CONCAVE        Special Case: tangential edgeuses
          //  CONVEX |\\\   inside -\\\\|   FACE CORNER     Classified as CONVEX CORNER
          //         +-----          \\\+-----               -----+-----
          /*     FACE CORNER          \\\\\\\\                \\\\\\\\\\                   */
          /*                            \\\\\\\                \\\\\\\\\                   */
          /*                                                                               */  
          SmBoolean bMixedCvxty = sm_CheckFaceComplimentaryCornerConcavity
                                    (m_cpVertex,    // in : target corner vertex
                                     pBaseEdgeuse1, // in : 1st target edgeuse - points outward from CornerVert
                                     pBaseEdgeuse2, // in : 2nd target edgeuse - points outward from CornerVert
                                     pSideEU,       // in : NULL    = EU1 and EU2 are immediate neighbors bounding same face
                                                    //      notNULL = OptRefEU lies 'between' EU1 and EU2 so that EU1->edge and
                                                    //        EU2->edge do not bound the same face. Needs to be a tangential edge.
                                     iDebugLevel ); // in :

          if ( bMixedCvxty )
            {
              return Make4x2ConcaveCornerTopology();
            }
        }
      m_bBevel = TRUE ;

    } // end 4 outward EUs from this corner vertex case check

  // Determine if two edges are tangent at corner vertex.
  // Set bIsTangentCase, dAngleDeg, and bMixedConvexity.

  double dAngleDeg;
  double dTangencyTolDeg = SM_RAD2DEG( m_dTangencyTolRadians );
  SmBoolean bIsTangentCase = sm_TestG1Continuity          // rtn: TRUE = Edge EndTangents at vertex are within tolerance of being tangent
                               (m_cpVertex,               // in : Target Vertex
                                pBaseEdgeuse1->GetEdge(), // in : 1st edge connecting to vertex
                                pBaseEdgeuse2->GetEdge(), // in : 2nd edge connectint to vertex
                                dTangencyTolDeg,          // in : Tangency tolerance in degrees
                                &dAngleDeg,               // out: angle (degrees) between edges at vertex
                                iDebugLevel) ;            // in : 

  // Store this value for later use.
  m_bG1Case = bIsTangentCase;

  SmBoolean bMixedConvexity = sm_CheckFaceComplimentaryCornerConcavity
                                (m_cpVertex,     // in : target corner vertex
                                 pBaseEdgeuse1,  // in : 1st target edgeuse - points outward from CornerVert
                                 pBaseEdgeuse2,  // in : 2nd target edgeuse - points outward from CornerVert
                                 NULL,           // in : NULL    = EU1 and EU2 are immediate neighbors bounding same face
                                                 //      notNULL = OptRefEU lies 'between' EU1 and EU2 so that EU1->edge and
                                                 //        EU2->edge do not bound the same face. Needs to be a tangential edge.
                                 iDebugLevel) ;  // in :

  // Test if we need to create a corner surface for concave cases
  // Note: Only when angle is more than 5 degrees, cause otherwise
  // it will be a very skinny surface & may cause problems in later process
//    if (!bIsTangentCase && dAngleDeg < 5.0)
//    if (!bIsTangentCase && dAngleDeg > 5.0)
//

  // Get Solvers and FilletGeoms for the two filleted edges at this vertex
  SmFilletSolver * pFilSolver1  = GetFilletExecutive()->GetFilletSolverOfEdgeuse(pBaseEdgeuse1);
  SmFilletSolver * pFilSolver2  = GetFilletExecutive()->GetFilletSolverOfEdgeuse(pBaseEdgeuse2);
  SmFilletGeom   * pFilletGeom1 = pFilSolver1->GetFirstFilletGeom();
  SmFilletGeom   * pFilletGeom2 = pFilSolver2->GetFirstFilletGeom();

  // For the non-filleted edge(s), we expect the rails to hit it (them).
  // Create a (fillet) vertex there, along with two rail edgeuses.
  for (ii=0; ii<lTotalEUs; ii++)
    {
      SmEdgeuse * pCurrEU = (*m_pAllEUs)[ii];

      if(   pCurrEU != pBaseEdgeuse1 
         && pCurrEU != pBaseEdgeuse2 )
        {
          // This edgeuse is in between two filletted edgeuses
          SmEdgeuse      * pSideEU    = sm_FindSideEdgeuse( pBaseEdgeuse1, pCurrEU, iDebugLevel );

          //
          SmFilletVertex * pNewSideFV = new (m_pPseudoBrep) SmFilletVertex( this );
          pNewSideFV->SetPointClass( SM_PC_EDGEUSE, pSideEU );
          pNewSideFV->SetFilletVertexType( SM_FV_RAIL_X_EDGEUSE );
          m_vVertices.Add( pNewSideFV );

          ULONG             lRailIndex      = pFilSolver1->FindIndexOfRailXSideEdgeuse( pSideEU );
          ULONG             lOtherRailIndex = pFilSolver2->FindIndexOfRailXSideEdgeuse(pSideEU->GetRadial() );
          SmFilletEdgeuse * pNewEdgeuse     = NULL;

          // Make a FilletEdgeuse and a FilletVertexuse to connect 
          // a given filletVertex to one of this FilletGeom's Fillet Rails.
          SER( pFilletGeom1->MakeRailEdgeuse( m_pPseudoBrep, lRailIndex,      pNewSideFV, pNewEdgeuse ));
          SER( pFilletGeom2->MakeRailEdgeuse( m_pPseudoBrep, lOtherRailIndex, pNewSideFV, pNewEdgeuse ));

        } // end pCurrEU is a side EU check
    } // end iter all outgoing edgeuses from this vertex 

  // Create the topology for the fillet corner.
  if ( lTotalEUs == 3 )
    {
      // In the normal uniform-convexity case, two rails will intersect
      // in the common face, and the other two rails should intersect the
      // un-filleted side edge at the same place.
      // We make vertices at those two places, and an edge between them.
      // For tangent and mixed-convexity cases ... see comments below.

      // gwc: assuming two rails should intersect the un-filleted side edge AT THE SAME PLACE isn't always true.  
      // This gets detected and processed with virtual call SmFilletCorner::AdjustTopology() made in SmFilletExecutive::CreateFilletBrep().

      // Find the side EU -- the un-filleted one.
      SmEdgeuse * pSideEU = NULL;
      for(ii=0;ii<lTotalEUs;ii++)
        {
          SmEdgeuse * pCurrEU = (*m_pAllEUs)[ii];
          if ( pCurrEU != pBaseEdgeuse1 && pCurrEU != pBaseEdgeuse2 )
            {
              pSideEU = pCurrEU;
              break;
            }
        } // end iter all EUs looking for the side (nonFilleted) edgeuse

      // locals
      SmFace          * pCommonFace     = sm_GetFaceByTwoEdgeuses(pBaseEdgeuse1, pBaseEdgeuse2, pSideEU );
      ULONG             lRailIndex      = pFilSolver1->FindIndexOfRailOnFace(pCommonFace);
      ULONG             lOtherRailIndex = pFilSolver2->FindIndexOfRailOnFace(pCommonFace);
      SmFilletVertex  * pNewSideFV      = m_vVertices[0];
      SmFilletVertex  * pNewFaceFV      = new (m_pPseudoBrep) SmFilletVertex(this);
      SmFilletEdgeuse * pNewEdgeuse     = NULL;

      pNewFaceFV->SetPointClass( SM_PC_FACE, pCommonFace );
      m_vVertices.Add( pNewFaceFV );

      // At a corner of three edges only the Join is built here: the blends continue each other,
      // or meet at 5 degrees or less where the face's corner is of the other convexity
      // (SmBlendCorner::Classify). A Turn is MakeMitreTopology's or MakeSweptTopology's.
      if ( !(   bIsTangentCase
             || (bMixedConvexity && dAngleDeg <= 5.0) ) )
        { SER( SM_ERR ); }
      m_bBlending = TRUE;
      // Two filleted edgeuses are tangent to each other.
      // The new face-vertex will just be straight across from the
      // rail-side edge intersection vertex.
      pNewFaceFV->SetFilletVertexType( SM_FV_MATE );
      pNewFaceFV->SetMate( 0, pNewSideFV );
      pNewSideFV->SetMate( 0, pNewFaceFV );

      // Setup rails of fillet
      //   Make a FilletEdgeuse and a FilletVertexuse to connect 
      //   a given filletVertex to one of this FilletGeom's Fillet Rails.
      SER( pFilletGeom1->MakeRailEdgeuse( m_pPseudoBrep, lRailIndex, pNewFaceFV, pNewEdgeuse));
      SER( pFilletGeom2->MakeRailEdgeuse( m_pPseudoBrep, lOtherRailIndex, pNewFaceFV, pNewEdgeuse));

    }  // end lTotalEUs == 3 branch

  else if ( lTotalEUs == 4 )
    {
      // Note, at this point we don't know what the topology will be.
      // There are side-edges between both pairs of rails.
      // (Note, we're assuming exactly one side-edge between each
      // potentially-intersecting rails: we exluded other case above.)
//cbi: haven't implemented this yet:
      // If the two side-edge-rail-rail intersections are straight across
      // from each other, then the two fillet surfaces will meet, at a
      // single rib.  Otherwise, we will create a third fillet face to join
      // the two.

      // Set m_bBlending and m_bBevel in 4x2 case as well. [bd 050718]
      if (   bIsTangentCase
                           //dAngleDeg <= 2.0
          || (bMixedConvexity && dAngleDeg <= 5.0) )
        {
          // Two filleted edgeuses are tangent to each other
          m_bBlending = TRUE;
          m_bBevel = FALSE;
        }
      else
        {
          // Regular BEVEL corner
          m_bBevel = TRUE;
          m_bBlending = FALSE;
        }
    }  // end lTotalEUs == 4 branch


  // At this point, we have only one fillet (cross) edge at this corner.
  // In the 3x2 case, it will run from the rail-rail intersection of the
  // two rails in the common face to the rail/rail/side-edge intersection.
  if(m_vVertices.GetSize() != 2)
    { SER(SM_ERR) ; }

  // Create the new cross-edge.
  SmFilletEdge * pNewFilletCrossEdge = NULL;
  SER( SmFilletEdge::MakeFilletEdge(m_pPseudoBrep,        // in : target Brep to receive new topology objects
                                    m_vVertices[0],       // in : start of new FilletEdge
                                    m_vVertices[1],       // in : end   of new FilletEdge
                                    pNewFilletCrossEdge,  // out: newly allocated FilletEdge - no Curve or UVTrimCurve data
                                    this ));              // in : NewFilletEdge's filletCorner, gets stored in rpNewFilletEdge->m_cpCorner
                                                          //      NULL to ignore, default:[NULL]
  SmFilletEdgeuse * pPrimEU = (SmFilletEdgeuse*)pNewFilletCrossEdge->GetPrimaryEdgeuse() ;
  SmFilletEdgeuse * pMateEU = (SmFilletEdgeuse*)pPrimEU->GetMate() ;
  pPrimEU->SetFilletGeom(pFilletGeom1) ;
  pMateEU->SetFilletGeom(pFilletGeom2) ;

  // Note: in version 5, this said 'if ( bIsTangentCase )', but was changed
  // to 'if ( m_bBlending )'.  In the 3x2 case, they're the same, but not
  // in the 4x2 case.  'Or' is correct here.
  // [bd 09 Nov 05, 050718]

  // if (m_bBlending)
  if ( bIsTangentCase || m_bBlending )
    {
      pNewFilletCrossEdge->SetFilletEdgeType( SM_FE_CROSS_SECTION );
    }
  else
    {
      pNewFilletCrossEdge->SetFilletEdgeType( SM_FE_FILLET_X_FILLET );
    }

  m_vEdges.Add( pNewFilletCrossEdge );

  // all done
  return SM_SUCCESS;

} // end SmBlendCorner::MakeTwoBlendTopology

/*******************************************************************//**
PURPOSE: If this FilletCorner has exactly two FilletSolvers, and their Edges
   meet with G1 continuity, return the other one than the passed-in argument.

NOTES:
   If the passed-in FilletSolver is not one of ours, return Null with an error.
   If there is only one FilletSolver, and it is the one passed in, return it.
***********************************************************************/
SmFilletSolver * SmBlendCorner::GetTangentFilletSolver( const SmFilletSolver *pFS )
{
  if ( m_pSolverEUs->GetSize() != 2 )  { SM_ASSERT_ERR; return NULL; }

  // We have to test its continuity.
  SmEdgeuse      * pEU0 = (*m_pSolverEUs)[0];
  SmEdge *pEdge0 = pEU0->GetEdge();
  SmEdgeuse      * pEU1 = (*m_pSolverEUs)[1];
  SmEdge *pEdge1 = pEU1->GetEdge();

  double dTangencyTolDeg = SM_RAD2DEG( m_dTangencyTolRadians );
  int iDebugLevel = 0;  // ... enhance as needed.
  if ( ! sm_TestG1Continuity( m_cpVertex, pEdge0, pEdge1, dTangencyTolDeg, NULL, iDebugLevel ))
    { return NULL; }   // not an error.


  SmFilletSolver * pFS0 = m_pExecutive->GetFilletSolverOfEdgeuse( pEU0 );
  SmFilletSolver * pFS1 = m_pExecutive->GetFilletSolverOfEdgeuse( pEU1 );

  if ( pFS0 == pFS )
    { return pFS1; }
  if ( pFS1 == pFS )
    { return pFS0; }

  SM_ASSERT_ERR; return NULL;

} // end SmBlendCorner::GetTangentFilletSolver

/*******************************************************************//**
PURPOSE: Make corner topology for 4x2-Concave Corner cases.
    These (rare-)cases are very similiar to 3x2-Concave cases except that
    an extra edge (very likely a topological edge) is connected to
    the corner vertex

NOTES:
***********************************************************************/
SmStatus SmBlendCorner::Make4x2ConcaveCornerTopology
  ()
{
int iDebugLevel = 0;
#ifdef SM_DEBUG_CODE
iDebugLevel = DebugLevel();
#endif // SM_DEBUG_CODE

    SmEdgeuse * pBaseEdgeuse1 = (*m_pSolverEUs)[0];
    SmEdgeuse * pBaseEdgeuse2 = (*m_pSolverEUs)[1];

    // Get Solvers
    SmFilletSolver * pFilSolver1 = GetFilletExecutive()->
        GetFilletSolverOfEdgeuse(pBaseEdgeuse1);
    SmFilletSolver * pFilSolver2 = GetFilletExecutive()->
        GetFilletSolverOfEdgeuse(pBaseEdgeuse2);
    SmFilletGeom * pFilletGeom1 = pFilSolver1->GetFirstFilletGeom();
    SmFilletGeom * pFilletGeom2 = pFilSolver2->GetFirstFilletGeom();

    SmFace * pCommonFace = NULL;
    SmFilletVertex * pFVOnCommonFace = NULL;
    SmFilletVertex * pFVOnConvexEdge = NULL;
    SmFilletGeom * pConvexFG = NULL;
    SmFilletGeom * pConcaveFG = NULL;

    // Make topology for all (both) of the non-filleted edges at this vertex.

    ULONG ii;
    for (ii=0; ii<4; ii++)
    {
        SmEdgeuse * pCurrEU = (*m_pAllEUs)[ii];
        if (pCurrEU == pBaseEdgeuse1 || pCurrEU == pBaseEdgeuse2)
        {
            continue;  // Do only the non-filleted ones.
        }

        // This edgeuse is in between two filletted edgeuses
        SmFilletVertex * pFV = new (m_pPseudoBrep) SmFilletVertex(this);
        m_vVertices.Add(pFV);

        SmEdgeuse * pSideEU = sm_FindSideEdgeuse(pBaseEdgeuse1, pCurrEU, iDebugLevel );
        pFV->SetPointClass(SM_PC_EDGEUSE,pSideEU);
        pFV->SetFilletVertexType(SM_FV_RAIL_X_EDGEUSE);

        SmFilletVertex * pFV2 = NULL;

        if ( sm_CheckFaceComplimentaryCornerConcavity(
                m_cpVertex, pBaseEdgeuse1, pBaseEdgeuse2, pCurrEU, iDebugLevel )
           )
        {
            double dAngle1 = sm_FindAngleBetweenEdgeuses(
                    m_cpVertex, pCurrEU, pBaseEdgeuse1, iDebugLevel );

            double dAngle2 = sm_FindAngleBetweenEdgeuses(
                    m_cpVertex, pCurrEU, pBaseEdgeuse2, iDebugLevel );

            //if (sm_CheckFaceComplimentaryCornerConcavity(
                    //m_cpVertex, pCurrEU, pBaseEdgeuse1, NULL, iDebugLevel ))

            if ( dAngle1 > dAngle2 )
            {
                pCommonFace = sm_GetFaceByTwoEdgeuses(
                    pCurrEU,pBaseEdgeuse1,NULL);
                pFV->SetPointClass(SM_PC_EDGEUSE,pSideEU->GetRadial());

                pFVOnCommonFace = new (m_pPseudoBrep) SmFilletVertex(this);
                pFVOnCommonFace->SetPointClass(SM_PC_FACE, pCommonFace);

                pFV2 = pFVOnCommonFace;
                pFV2->SetFilletVertexType(SM_FV_MATE);
                SM_SWAP_PTR(SmFilletVertex,pFV,pFV2);
                pConcaveFG = pFilletGeom1;
                pConvexFG = pFilletGeom2;
            }

            //else if (sm_CheckFaceComplimentaryCornerConcavity(
                    //m_cpVertex, pCurrEU, pBaseEdgeuse2, NULL, iDebugLevel ))

            else // if ( dAngle1 < dAngle2 )  [bd 060424]
            {
                pCommonFace = sm_GetFaceByTwoEdgeuses(
                    pCurrEU,pBaseEdgeuse2,NULL);

                pFVOnCommonFace = new (m_pPseudoBrep) SmFilletVertex(this);
                pFVOnCommonFace->SetPointClass(SM_PC_FACE, pCommonFace);

                pFV2 = pFVOnCommonFace;
                pFV2->SetFilletVertexType(SM_FV_MATE);
                pConcaveFG = pFilletGeom2;
                pConvexFG = pFilletGeom1;
            }
        }
        else
        {
            pFV2 = pFV;
            pFVOnConvexEdge = pFV;
        }

        ULONG lRailIndex = pFilSolver1->FindIndexOfRailXSideEdgeuse(pSideEU);
        SmFilletEdgeuse * pNewEdgeuse = NULL;
        SER(pFilletGeom1->MakeRailEdgeuse(m_pPseudoBrep,
            lRailIndex,pFV,pNewEdgeuse));
        ULONG lOtherRailIndex = pFilSolver2->FindIndexOfRailXSideEdgeuse(
            pSideEU->GetRadial());
        SER(pFilletGeom2->MakeRailEdgeuse(m_pPseudoBrep,
            lOtherRailIndex,pFV2,pNewEdgeuse));

    } // end loop on all (both) non-filleted edges.

    NER(pFVOnCommonFace);
    NER(pFVOnConvexEdge);
    pFVOnCommonFace->SetMate(0, pFVOnConvexEdge);
    pFVOnConvexEdge->SetMate(0, pFVOnCommonFace);
    m_vVertices.Add(pFVOnCommonFace);

    // Will create a 3-sided corner face
    SmFilletVertex * pPrevFV = pFVOnCommonFace; //i.e. m_vVertices[2]
    for (ii=0; ii<3; ii++)
    {
        SmFilletVertex * pCurrFV = m_vVertices[ii];
        SmFilletEdge * pFilletEdge = NULL;
        if (ii==1 || pCurrFV == pFVOnConvexEdge || pPrevFV == pFVOnConvexEdge)
        {
            // Make Fillet Edge
            SER(SmFilletEdge::MakeFilletEdge(m_pPseudoBrep,  // in : target Brep to receive new topology objects
                                             pPrevFV,        // in : start of new FilletEdge
                                             pCurrFV,        // in : end   of new FilletEdge
                                             pFilletEdge,    // out: newly allocated FilletEdge - no Curve or UVTrimCurve data
                                             this));         // in : NewFilletEdge's filletCorner, gets stored in rpNewFilletEdge->m_cpCorner
                                                             //      NULL to ignore, default:[NULL]
            pFilletEdge->SetFilletEdgeType(SM_FE_CROSS_SECTION);
            SmFilletEdgeuse * pPrimEU = (SmFilletEdgeuse*)pFilletEdge->
                GetPrimaryEdgeuse();
            if (ii == 1) {
                pPrimEU->SetFilletGeom(pConvexFG);
            }
            else {
                pPrimEU->SetFilletGeom(pConcaveFG);
            }
        }
        else
        {
            SER(SmFilletEdge::MakeFilletEdge(m_pPseudoBrep,  // in : target Brep to receive new topology objects
                                             pPrevFV,        // in : start of new FilletEdge
                                             pCurrFV,        // in : end   of new FilletEdge
                                             pFilletEdge,    // out: newly allocated FilletEdge - no Curve or UVTrimCurve data
                                             this));         // in : NewFilletEdge's filletCorner, gets stored in rpNewFilletEdge->m_cpCorner
                                                             //      NULL to ignore, default:[NULL]
            pFilletEdge->SetFilletEdgeType(SM_FE_RAIL_RAIL_INTERPOLATION);
            pFilletEdge->SetOriginalFace(pCommonFace);
        }
        m_vEdges.Add(pFilletEdge);
        pPrevFV = pCurrFV;
    }

    return SM_SUCCESS;

} // end SmBlendCorner::Make4x2ConcaveCornerTopology

/*******************************************************************//**
PURPOSE: Adjust the corner topology for the mixed-convexity cases where two fillets
   of the corner does not intersect each other on the side edge.
   One additional corner edge will be added to fill the 'gap'.

NOTES:
***********************************************************************/
SmStatus SmBlendCorner::AdjustConcaveCornerTopology()
{
  // Find the fillet vertex to be adjusted, pAdjustVert.
  // It's the only SM_FV_RAIL_X_EDGEUSE.

  SmFilletVertex * pAdjustVert = NULL;
  SmTArray<SmFilletGeom*> sGeoms;
  ULONG ii;
  for ( ii=0; ii<m_vVertices.GetSize(); ii++ )
  {
      SmFilletVertex * pV = m_vVertices[ii];
      if ( pV->GetFilletVertexType() != SM_FV_RAIL_X_EDGEUSE )
        { continue; }

      pV->GetFilletGeoms( sGeoms );
      if ( sGeoms.GetSize() != 2 )
      { continue; }  //SER(SM_ERR);
      pAdjustVert = pV;
      break;
  }
  NER(pAdjustVert);

  SmFilletGeom * pFilletGeom1 = sGeoms[0];
  SmFilletGeom * pFilletGeom2 = sGeoms[1];
  SmEdgeuse * pSideEU1 = (SmEdgeuse*)pAdjustVert->GetPointClassObject();
  SmFilletVertexuse * pVU1 = pAdjustVert->GetVUAtRailEnd( pFilletGeom1 );
  SmFilletVertexuse * pVU2 = pAdjustVert->GetVUAtRailEnd( pFilletGeom2, pVU1 );
  double dEdgeParam1 = pVU1->GetTsectPnt().m_dCurveParameter;
  double dEdgeParam2 = pVU2->GetTsectPnt().m_dCurveParameter;
  SmPoint3d sPnt1, sPnt2;
  SmCurve * pSideEdgeCurve = pSideEU1->GetEdge()->GetCurve();
  SER( pSideEdgeCurve->EvaluatePoint( dEdgeParam1,sPnt1 ) );
  SER( pSideEdgeCurve->EvaluatePoint( dEdgeParam2,sPnt2 ) );

  // See if we have coincident geometry, i.e. two rails intersect at this vertex.
  double dGapDist = sPnt1.DistanceBetween( sPnt2 );
  if ( dGapDist < m_dThisApproxTol3d )
  {
      // No adjustments are needed
      return SM_SUCCESS;
  }

#ifdef SM_DEBUG_CODE
  if ( DebugLevel() > 0 ) {
      smgfx_SetLook(1,2, 1,0,0); pAdjustVert->Draw(); sm_GraphicsLoop();
      sm_GraphicsLoop();
  }
#endif // SM_DEBUG_CODE

  // Create two new FilletVertices where each rail hits the side edge.
  SmFilletVertex * pNewVertL = new (m_pPseudoBrep) SmFilletVertex( this );  // "Large"
  pNewVertL->SetFilletVertexType( SM_FV_RAIL_X_EDGEUSE );
  SmFilletVertex * pNewVertS = new (m_pPseudoBrep) SmFilletVertex( this );  // "Small"
  pNewVertS->SetFilletVertexType( SM_FV_RAIL_X_EDGEUSE );
  m_vVertices.Add( pNewVertL );
  m_vVertices.Add( pNewVertS );

  pAdjustVert->SetFilletVertexType( SM_FV_FILLET_X_FILLET );
  pAdjustVert->SetStatus( SM_FIL_UNPROCESSED );

  // See which of the two rail-side edge intersections is closer
  // to the filleted vertex, to decide Large/Small.
  SmPoint3d sVertPnt = m_cpVertex->GetPoint();
  SmFilletGeom * pLargeFG = NULL;
  SmFilletGeom * pSmallFG = NULL;
  SmVertexuse * pVUOnLargeFG = NULL;
  SmVertexuse * pVUOnSmallFG = NULL;
  SmFace * pSideFace = NULL;
  if ( sPnt1.DistanceBetween( sVertPnt ) < sPnt2.DistanceBetween( sVertPnt ) )
  {
      pNewVertL->SetPoint( sPnt2 );
      pNewVertL->SetOriginalTParam( dEdgeParam2 );
      pNewVertS->SetPoint( sPnt1 );
      pNewVertS->SetOriginalTParam( dEdgeParam1 );
      pLargeFG = pFilletGeom2;
      pSmallFG = pFilletGeom1;
      pVUOnLargeFG = pVU2;
      pVUOnSmallFG = pVU1;
      pSideFace = pSideEU1->GetFace();
  }
  else
  {
      pNewVertL->SetPoint( sPnt1 );
      pNewVertL->SetOriginalTParam( dEdgeParam1 );
      pNewVertS->SetPoint( sPnt2 );
      pNewVertS->SetOriginalTParam( dEdgeParam2 );
      pLargeFG = pFilletGeom1;
      pSmallFG = pFilletGeom2;
      pVUOnLargeFG = pVU1;
      pVUOnSmallFG = pVU2;
      pSideFace = pSideEU1->GetRadial()->GetFace();
  }
  pNewVertL->SetStatus( SM_FIL_PROCESSED );
  pNewVertS->SetStatus( SM_FIL_PROCESSED );
  SER( pAdjustVert->Remove( pVUOnLargeFG ) );
  SER( pNewVertL->PostInsert( pVUOnLargeFG ));
  SER( pAdjustVert->Remove( pVUOnSmallFG ) );
  SER( pNewVertS->PostInsert( pVUOnSmallFG ));

  SmFilletEdge * pNewEdge = NULL;
  SER( SmFilletEdge::MakeFilletEdge(m_pPseudoBrep, // in : target Brep to receive new topology objects
                                    pAdjustVert,   // in : start of new FilletEdge
                                    pNewVertL,     // in : end   of new FilletEdge
                                    pNewEdge,      // out: newly allocated FilletEdge - no Curve or UVTrimCurve data
                                    this));        // in : NewFilletEdge's filletCorner, gets stored in rpNewFilletEdge->m_cpCorner
                                                   //      NULL to ignore, default:[NULL]
  pNewEdge->SetFilletEdgeType( SM_FE_FILLET_X_SIDE_FACE );
  pNewEdge->SetOriginalFace( pSideFace );
  SmFilletEdgeuse * pPrimEU = (SmFilletEdgeuse*)pNewEdge->GetPrimaryEdgeuse();
  pPrimEU->SetFilletGeom( pLargeFG );

  // Need to extend the side surface
  SmSurface * pSideSurface = pSideFace->GetSurface();
  SmSurface * pNewSurface = NULL;
  SER( pSideSurface->CreateExtendedSurface( m_crContext, 10.0*dGapDist, SM_CT_G1, pNewSurface) );
  SmObjDelete sCleanSurf( pNewSurface );
  //SmBSplineSurface * pSideSurface = SM_CAST_PTR( SmBSplineSurface, pSideFace->GetSurface() );
  //SmBSplineSurface * pNewSurface = NULL;
  //SER( pSideSurface->CreateExtendedSurface( m_crContext, 10.0*dGapDist, SM_CT_G1, pNewSurface) );
  //SmObjDelete sCleanSurf( pNewSurface );
  SmSurface *pNewAnalyticSurface;

  // Copy pNewSurface, when possible as an analytic surface.
  SER( pNewSurface->CopyAndAddAnalytics( m_crContext, pNewAnalyticSurface ));
  AddExtendedOriginalSurface(pSideFace,pNewAnalyticSurface);
  pNewAnalyticSurface = NULL; // to be safe: might have been deleted.

  SmFilletEdge * pNewEdge1 = NULL;
  SER( SmFilletEdge::MakeFilletEdge(m_pPseudoBrep,  // in : target Brep to receive new topology objects
                                    pNewVertL,      // in : start of new FilletEdge
                                    pNewVertS,      // in : end   of new FilletEdge
                                    pNewEdge1,      // out: newly allocated FilletEdge - no Curve or UVTrimCurve data
                                    this));         // in : NewFilletEdge's filletCorner, gets stored in rpNewFilletEdge->m_cpCorner
                                                    //      NULL to ignore, default:[NULL]
  pNewEdge1->SetFilletEdgeType( SM_FE_ON_EDGE );
  pNewEdge1->SetOriginalEdge( pSideEU1->GetEdge() );
  pNewEdge1->SetOriginalFace( pSideFace );
  // Set this edge as a topo edges so that it can be removed in the merge process
  //m_pExecutive->AddTopoEdges( pNewEdge1 );

  SmFilletEdge * pNewEdge2 = NULL;
  SER( SmFilletEdge::MakeFilletEdge(m_pPseudoBrep, // in : target Brep to receive new topology objects
                                    pNewVertS,     // in : start of new FilletEdge
                                    pAdjustVert,   // in : end   of new FilletEdge
                                    pNewEdge2,     // out: newly allocated FilletEdge - no Curve or UVTrimCurve data
                                    this));        // in : NewFilletEdge's filletCorner, gets stored in rpNewFilletEdge->m_cpCorner
                                                   //      NULL to ignore, default:[NULL]
  pNewEdge2->SetFilletEdgeType( SM_FE_RAIL_EXTENSION );
  pNewEdge2->SetOriginalFace( pSideFace );
  SmFilletEdgeuse * pPrimEU2 = (SmFilletEdgeuse*)pNewEdge2->GetPrimaryEdgeuse();
  pPrimEU2->SetFilletGeom( pSmallFG );

  m_vEdges.Add( pNewEdge  );
  m_vEdges.Add( pNewEdge1 );
  m_vEdges.Add( pNewEdge2 );

  // pAdjustVert comes in with two mates, in the common face.
  // Find them -- they are at the other end of the two FilletEdges
  // of type SM_FE_CROSS_SECTION -- and reset them as mates of
  // the two FilletVertices we just created.
  //SmVertexuse * pSmallFilletVU = NULL;
  //SmVertexuse * pLargeFilletVU = NULL;
  SmTArray<SmVertexuse*> sVUs;
  pAdjustVert->GetVertexuses( sVUs );
  for ( ii=0; ii<sVUs.GetSize(); ii++ )
  {
      SmVertexuse * pVU = sVUs[ii];
      SmFilletEdgeuse * pEU = (SmFilletEdgeuse*)pVU->GetEdgeuse();
      SmFilletEdge * pFilletE = (SmFilletEdge*)pEU->GetEdge();
      if ( pFilletE->GetFilletEdgeType() != SM_FE_CROSS_SECTION )
        { continue; }

      // FilletEdge of type SM_FE_CROSS_SECTION:
      SmFilletVertex * pOtherVert = (SmFilletVertex*)pFilletE->GetOtherVertex( pAdjustVert );
      SmFilletGeom * pFG = pEU->GetFilletGeom();
      if ( pFG == NULL )
      {
          pEU = (SmFilletEdgeuse*)pEU->GetMate();
          pFG = pEU->GetFilletGeom();
          NER( pFG );
      }
      if ( pFG == pSmallFG && pOtherVert->GetMate(0) != NULL )
      {
          pNewVertS->SetMate(0,pOtherVert);
          pOtherVert->SetMate(0,pNewVertS);
      }
      else if ( pFG == pLargeFG && pOtherVert->GetMate(0) != NULL )
      {
          pNewVertL->SetMate(0,pOtherVert);
          pOtherVert->SetMate(0,pNewVertL);
      }
  }
  pAdjustVert->SetMate( 0, NULL );
  pAdjustVert->SetMate( 1, NULL );

  // Allow fillets to be extended at the corner
  m_bBlending = FALSE;

  m_bTopologyAdjusted = TRUE;

  return SM_SUCCESS;

} // end SmBlendCorner::AdjustConcaveCornerTopology

/*******************************************************************//**
PURPOSE: Adjust a corner for the Nx2 tangent case with N > 3.

NOTES:
    In that case, both rails are probably rail-int-edgeuse.
    If everything matches up just right geometrically,
    then we can mate the two fillet surfaces together,
    but there's a good chance that things don't meet properly,
    because the fillets are not rolling on the same face.
    In that case, we'll have to put another fillet surface
    between the two.

cbi: Plan: make a new Corner type, SM_FCR_TANGENT, subclass SmFilletTangentCorner.
     (possibly subtypes 3x2Tangent and Nx2Tangent.)
     Then this method will be SmFilletTangentCorner::AdjustTopology().
     Then this method won't have to check whether it's a Nx2 tangent corner;
     that check will be in Create().
cbi: try a 5x2 corner?
***********************************************************************/
SmStatus SmBlendCorner::AdjustTwoBlendTangentCorner()
{
    // First check that this is the appropriate type.
    // Must be Nx2 tangent corner, and currently N must be 4.
    if ( m_pSolverEUs == NULL  ||  m_pSolverEUs->GetSize() != 2 ) { return SM_SUCCESS; }
    if ( m_pAllEUs    == NULL  ||  m_pAllEUs   ->GetSize() != 4 ) { return SM_SUCCESS; }

    // Check that the two filleted edges are tangent at this vertex.

    // Note: using NormalizedEvaluate() assumes that m_pSolverEUs go outward
    // from the filleted vertex.
    SmPoint3d sPt;
    SmVector3d s1stDeriv1, s1stDeriv2;
    SER( (*m_pSolverEUs)[0]->NormalizedEvaluate( 0.0, FALSE, sPt, &s1stDeriv1 ) );   // TRUE = UV Eval, FALSE = 3d Eval
    SER( (*m_pSolverEUs)[1]->NormalizedEvaluate( 0.0, FALSE, sPt, &s1stDeriv2 ) );   // TRUE = UV Eval, FALSE = 3d Eval

    double dAngTolDeg = 2.0; //cbi: some places use 5 degrees; sort that all out.
    if ( !s1stDeriv1.IsParallelTo( s1stDeriv2, dAngTolDeg)  )
        { return SM_SUCCESS; }


    // Ok, we are an Nx2-tangent corner.

    // loop on all of our FilletVertices
    //   looking for surface points being different.

    SmBoolean bCloseEnough = TRUE;

    // Locals
    SmTArray<SmFilletGeom*> sGeoms;

#ifdef SM_DEBUG_CODE
    int iDebugLevel = DebugLevel();
#endif // SM_DEBUG_CODE

    ULONG ii;
    for ( ii=0; ii<m_vVertices.GetSize(); ii++ )
    {
        SmFilletVertex * pAdjustVert = m_vVertices[ii];

        // Check whether this is the vertex that corresponds to the
        // rail hitting the side edge, and not the rails hitting
        // each other in a common face.
        if ( pAdjustVert->GetFilletVertexType() != SM_FV_RAIL_X_EDGEUSE )
          { continue; }

#ifdef SM_DEBUG_CODE
        if ( iDebugLevel > 0 ) {
            smgfx_SetLook(1,2, 1,0,0); pAdjustVert->Draw(); sm_GraphicsLoop();
            sm_GraphicsLoop();
        }
#endif // SM_DEBUG_CODE

        // Ok, this fillet vertex corresponds to where the
        // rails from the two fillets hit a side edge.
        // Check whether those rails intersect the side edge at the same place.
        // In addition, the other end of this rib will have to coincide
        // in order for the fillet surfaces to mate here.

        pAdjustVert->GetFilletGeoms(sGeoms);
        if ( sGeoms.GetSize() != 2 )
            continue;  // SER( SM_ERR );  maybe already adjusted...

        SmFilletGeom * pFilletGeom1 = sGeoms[0];
        SmFilletGeom * pFilletGeom2 = sGeoms[1];

        //cbi: fix this in the future
        // We need the rail curves to work from.
        // If not there yet, just wait for a later call.
        if (   pFilletGeom1->GetRail(0)->GetCurve() == NULL
           ||  pFilletGeom1->GetRail(1)->GetCurve() == NULL
           ||  pFilletGeom2->GetRail(0)->GetCurve() == NULL
           ||  pFilletGeom2->GetRail(1)->GetCurve() == NULL )
        {
            return SM_SUCCESS;
        }

        // Get the vertexuses corresponding to the two surfaces,
        // and their TSectPnts:
        SmFilletVertexuse * pVU1 = pAdjustVert->GetVUAtRailEnd(pFilletGeom1);
        SmFilletVertexuse * pVU2 = pAdjustVert->GetVUAtRailEnd(pFilletGeom2, pVU1);
        SmTsectPnt &rIntPtSrf1 = pVU1->GetTsectPnt();
        SmTsectPnt &rIntPtSrf2 = pVU2->GetTsectPnt();

        SmSurface *pSrf1A = (SmSurface*)( rIntPtSrf1.m_apUserPointer[0] );
        SmSurface *pSrf1B = (SmSurface*)( rIntPtSrf1.m_apUserPointer[1] );
        SmSurface *pSrf2A = (SmSurface*)( rIntPtSrf2.m_apUserPointer[0] );
        SmSurface *pSrf2B = (SmSurface*)( rIntPtSrf2.m_apUserPointer[1] );

#ifdef SM_DEBUG_CODE
        if ( iDebugLevel > 0 ) {
            smgfx_SetLook(2,1, 0,0,1); pSrf1A->DrawUV(4,4,TRUE); sm_GraphicsLoop();
            smgfx_SetLook(2,1, 0,1,1); pSrf1B->DrawUV(4,4,TRUE); sm_GraphicsLoop();
            smgfx_SetLook(2,1, 1,0,0); pSrf2A->DrawUV(4,4,TRUE); sm_GraphicsLoop();
            smgfx_SetLook(2,1, 1,0,1); pSrf2B->DrawUV(4,4,TRUE); sm_GraphicsLoop();
            sm_GraphicsLoop();
        }
#endif // SM_DEBUG_CODE
        SmBoolean bTmp0 = FALSE, bTmp1 = FALSE, bTmp2 = FALSE, bTmp3 = FALSE ; 
        SmBSplineSurface * pBSplineSurface0 = SM_CAST_PTR(SmBSplineSurface, pSrf1A) ;
        SmBSplineSurface * pBSplineSurface1 = SM_CAST_PTR(SmBSplineSurface, pSrf1B) ;
        SmBSplineSurface * pBSplineSurface2 = SM_CAST_PTR(SmBSplineSurface, pSrf2A) ;
        SmBSplineSurface * pBSplineSurface3 = SM_CAST_PTR(SmBSplineSurface, pSrf2B) ;
        SmTemporaryChangeValue<SmBoolean> sClean0(pBSplineSurface0 ? pBSplineSurface0->GetOutOfBoundsEnabled() : bTmp0, TRUE) ;
        SmTemporaryChangeValue<SmBoolean> sClean1(pBSplineSurface1 ? pBSplineSurface1->GetOutOfBoundsEnabled() : bTmp1, TRUE) ;
        SmTemporaryChangeValue<SmBoolean> sClean2(pBSplineSurface2 ? pBSplineSurface2->GetOutOfBoundsEnabled() : bTmp2, TRUE) ;
        SmTemporaryChangeValue<SmBoolean> sClean3(pBSplineSurface3 ? pBSplineSurface3->GetOutOfBoundsEnabled() : bTmp3, TRUE) ;

        SmPoint3d sSrf1Pt, sSrf2Pt;
        SER( pSrf1A->EvaluatePoint( rIntPtSrf1.UVPos( 0 ), sSrf1Pt ) );
        SER( pSrf2A->EvaluatePoint( rIntPtSrf2.UVPos( 0 ), sSrf2Pt ) );

#ifdef SM_DEBUG_CODE
        if ( iDebugLevel > 0 ) {
            smgfx_SetLook(1,6, 0,0,1);
            sSrf1Pt.Draw(); sm_GraphicsLoop();
            sSrf2Pt.Draw(); sm_GraphicsLoop();
            sm_GraphicsLoop();
        }
#endif // SM_DEBUG_CODE

        if ( sSrf1Pt.DistanceBetween( sSrf2Pt ) > m_dThisApproxTol3d )
        {
            bCloseEnough = FALSE;
            break;
        }

        // The mates must agree as well for 4x2
        SER( pSrf1B->EvaluatePoint( rIntPtSrf1.UVPos( 1 ), sSrf1Pt ) );
        SER( pSrf2B->EvaluatePoint( rIntPtSrf2.UVPos( 1 ), sSrf2Pt ) );

#ifdef SM_DEBUG_CODE
        if ( iDebugLevel > 0 ) {
            smgfx_SetLook(1,6, 0,1,0);
            sSrf1Pt.Draw(); sm_GraphicsLoop();
            sSrf2Pt.Draw(); sm_GraphicsLoop();
            sm_GraphicsLoop();
        }
#endif // SM_DEBUG_CODE

        if ( sSrf1Pt.DistanceBetween( sSrf2Pt ) > m_dThisApproxTol3d )
        {
            bCloseEnough = FALSE;
            break;
        }
    
    } // end of loop on all of our FilletVertices
      //   looking for surface points being different.

    if ( bCloseEnough )
    {
        return SM_SUCCESS;
    }


    // Ok, we have to make a new fillet face in between.

    // Loop over all of our side vertices again, splitting any
    // fillet geom that the vertex is in the middle of.
    ULONG lWhichFV, lWhichFilletGeom;

    // When we process a vertex, we don't want to process its mate.
    // Create a separate list so we can remove what we don't want to process.
    SmTArray<SmFilletVertex*> apFVsToProcess( m_vVertices );

    for ( lWhichFV = 0; lWhichFV < apFVsToProcess.GetSize(); lWhichFV++ )
    {
        SmFilletVertex * pThisFV = apFVsToProcess[lWhichFV];

        // Looking for a vertex where a rail hits a side edge:
        if ( pThisFV->GetFilletVertexType() != SM_FV_RAIL_X_EDGEUSE )
            continue;

        pThisFV->GetFilletGeoms( sGeoms );
        SM_ASSERT( sGeoms.GetSize() == 2 );
        if ( sGeoms.GetSize() != 2 )
            { continue; }

        // This Fillet Vertex has two Fillet Vertexuses, one on each FilletGeom.
        // If either is interior to its FilletGeom, split that fillet geom.
        SmFilletVertexuse * apVUs[2];
        apVUs[0] = pThisFV->GetVUAtRailEnd( sGeoms[0] );
        apVUs[1] = pThisFV->GetVUAtRailEnd( sGeoms[1], apVUs[0] );

        // Drop the vertex position to the corresponding rails in
        // each FilletGeom.  If it's internal to the rail,
        // split the FilletGeom.

        SmPoint3d sVertexPos = pThisFV->GetPoint();

        // This is something we'll need after the next loop:
        SmFilletVertex *pOppositeFV = NULL;

        // Look at both FilletGeoms
        for ( lWhichFilletGeom = 0; lWhichFilletGeom <2; lWhichFilletGeom++ )
        {
            SmFilletGeom *pOrigGeom = sGeoms[ lWhichFilletGeom ];

            SmFilletVertexuse *pThisFVUse = apVUs[ lWhichFilletGeom ];

            // Grap the correct rail
            ULONG lRailIdx;
            SER( pOrigGeom->GetRailIndexOfVertex( pThisFV, lRailIdx ) );
            SmFilletEdge *pThisRail = pOrigGeom->GetRail( lRailIdx );

            SmExtent1d sRailIvl = pThisRail->GetInterval();
            double *pdGuessParam = NULL; //cbi.
            SmBoolean bSuccess;
            double dSplitParam = 0.0;
            double dDist = 0.0;

            SER( pThisRail->GetCurve()->DropPoint(sRailIvl,         // in : target curve allowed domain
                                             sVertexPos,       // in : Point to drop to curve
                                             NULL,             // in : For Snapping, NULL=NoSnap, Vec pointing to CrvInside from desired EndPt for snapping.
                                                               //      when DropPt is within dDistTol in cpVecToCrvInside dir of EndPt, Snap to Ivl EndPt.
                                                               //      Vec handles ambiguities and makes sure that curves are not just touching at the ends.
                                             SM_BIG_DOUBLE,    // in : Used for SnapDist, and SM_SO_INTERSECT sXSectTol3d distance.
                                                               //      sXSectTol3d = dDistTol ;      SM_SO_INTERSECT 
                                                               //      sXSectTol3d = SM_BIG_DOUBLE ; SM_SO_MINIMIZE, SM_SO_NORMALIZE
                                                               //      sSnapTol3d  = dDistTol ;      SM_SO_MINIMIZE, SM_SO_NORMALIZE, SM_SO_INTERSECT    
                                             pdGuessParam,     // in : If given, call cheap LocalPointSolve() else call expensive GlobalPointSolve()
                                             bSuccess,         // out: TRUE = found a drop point
                                             dSplitParam,      // out: found drop curve param
                                             dDist,            // out: found drop distance
                                             SM_SO_MINIMIZE)  ); // in : SM_SO_MINIMIZE = allow nonNormal drops near endPoints
                                                               //      SM_SO_INTERSECT= like Minimize but point must be within dDistTol
                                                               //      SM_SO_NORMALIZE= exclude nonNormal drops near endPoints
                                                               //      default:[SM_SO_MINIMIZE] to preserve original behavior

            if ( !bSuccess )
                { continue; }

            if ( sRailIvl.IsValueOnBoundary( dSplitParam, m_dThisApproxTol3d ) )
            {
                continue;
            }

            // Ok, here we have a fillet geom that needs to be split
            // at this FilletVertex (pThisFVUse).

            // Which end of the FilletGeom are we at?  0 for start, 1 for end.
            ULONG lWhichEnd = 0;
            if ( pOrigGeom->GetRailVertexuse( lRailIdx, 1 ) == pThisFVUse )
              { lWhichEnd = 1;  }

            // It makes a difference to the split routine whether the original
            // FG, or the new one, is a tangent-rollover.  The assumptions are:
            // 1: 'this' FG is the low end, and the new one is the high end; so,
            // 2: if we're moving from our base surfaces onto the adjacent face,
            //    then the new FG will be TangentRollover, otherwise the
            //    original FG is TangentRollover.
            // So we set the type of pOrigGeom according to what the low-end
            // part will be.
            // So, (only) if the low end is being changed, we change
            // the original FG to TangentRollover.
            if ( lWhichEnd == 0 )
                { pOrigGeom->SetFilletGeomType( SM_FG_TANGENT_ROLLOVER ); }


            SmFilletGeom *pNewGeom = NULL;
            SmStatus eStat =
                pOrigGeom->GetFilletSolver()->SplitFilletGeomAtSideEdge(
                    pOrigGeom,
                    pThisFVUse,
                    lWhichEnd,
                    &pNewGeom );

            // If failure, the whole fillet could still succeed, so don't SER.
            if ( eStat != SM_SUCCESS )
              { continue; }


            // After splitting, pOrigGeom is the low end and pNewGeom is the
            // high end: the newly created fillet vertices are at the high
            // end of pOrigGeom and the low end of pNewGeom.
            //
            // The fillet vertex at which we just split, pThisFV, is left at
            // either the low end of pOrigGeom (if lWhichEnd == 0) or the
            // high end of pNewGeom (otherwise).
            //
            // One of the new vertices now takes over for pThisFV, being the
            // intersection of the rail and the (tangent) side edge.
            // The split routine we just called takes care of that.
            //
            // The split fillet vertex, pThisFV, is now just a mate to the
            // corresponding FV on the other rail.  Deal with that.

            pThisFV->SetFilletVertexType( SM_FV_MATE );
            pThisFV->SetStatus( SM_FIL_UNPROCESSED );

            // Set it as mate of its opposite fillet vertex.
            SmFilletVertexuse *pOppositeFVU = NULL;

            if ( lWhichEnd == 0 )
                pOppositeFVU = pOrigGeom->GetRailVertexuse( 1-lRailIdx, 0 );
            else
                pOppositeFVU = pNewGeom->GetRailVertexuse( 1-lRailIdx, 1 );

            pOppositeFV = SM_CAST_PTR( SmFilletVertex, pOppositeFVU->GetVertex() );

            if ( pOppositeFV != NULL )
            {
                pThisFV->SetMate( 0, pOppositeFV );
                pOppositeFV->SetMate( 0, pThisFV );

                // Also tell it that its geometry is kaput
                // (because pThisFV is on a different base surface).
                pOppositeFV->SetStatus( SM_FIL_UNPROCESSED );
            }

            // Also, pThisFV is now in a face.
            SmFace *pNewBaseFace = pNewGeom->GetRail( lRailIdx )->GetOriginalFace();
            pThisFV->SetPointClass( SM_PC_FACE, pNewBaseFace );

            // Edges are kaput too.
            for ( ULONG j = 0; j < m_vEdges.GetSize(); j++ )
              { m_vEdges[j]->SetStatus( SM_FIL_UNPROCESSED ); }

            // Also: there are now two new FilletVertex's, which belong
            // to this FilletCorner.  Add them to our list, and set
            // their back-pointer to us.
            // 'Same' is on the same side as pThisFV, 'Opp' as pOppositeFV.
            // Also recall that the new FilletVertex's are at the high end
            // of the original FG -- hence the 2nd arguments:

            SmFilletVertex * pNewFVSame = pOrigGeom->GetRailVertex(   lRailIdx, 1 );
            SmFilletVertex * pNewFVOpp  = pOrigGeom->GetRailVertex( 1-lRailIdx, 1 );

            // pNewFVSame should be all set, it's what pThisFV used to be.
            // But pNewFVOpp is now in a face.
            SmFace *pOppBaseFace = pOrigGeom->GetRail( 1-lRailIdx )->GetOriginalFace();
            pNewFVOpp->SetPointClass( SM_PC_FACE, pOppBaseFace );

            // We also have to take responsibility for the new edge.
            SmTArray<SmEdge*> aNewEdges;
            pNewFVSame->GetCommonEdges( pNewFVOpp, aNewEdges );
            SM_ASSERT( aNewEdges.GetSize() == 1 );
            if ( aNewEdges.GetSize() > 0 )
            {
                SmFilletEdge *pThisFilEdge = SM_CAST_PTR( SmFilletEdge, aNewEdges[0] );
                this->m_vEdges.Add( pThisFilEdge );
                pThisFilEdge->SetFilletCorner( this );
            }

            // We have to update the fillet surfaces in the fillet geoms here.
            // A prerequisite to that is recalculating the vertex geom.

            SER( this->CalcCornerVertGeom() );


            // We need a direction for ReCalcFilletGeom.
            // We already dropped this vertex to the rail being split,
            // so just evaluate the rail's tangent direction there.

            SmVector3d sPtDer[2]; // point and 1st derivative
            SER( pThisRail->GetCurve()->Evaluate( dSplitParam, 1, FALSE, sPtDer ) );


            // Possible enhancement: We could avoid one ReCalcFilletGeom call.
            // The FG piece that's on the original two base faces could
            // just be trimmed, instead of recalculated.
            // If lWhichEnd is 1, then we could just trim pOrigGeom.
            // Also, if lWhichEnd is 0, we could move the geometry from
            // pOrigGeom into pNewGeom, then just trim pNewGeom,
            // and recalc pOrigGeom.

//cbi: this SetFilletGeomForSurfaces(): do it right in ReCalc?
            // Tell FilletSolver to use surfaces in the FilletGeoms.
            // First back up what's there, if anything, to restore when done.
            SmFilletSolver *pFS = pOrigGeom->GetFilletSolver();
            SmFilletGeom *pBackupFG = pFS->GetFilletGeomForSurfaces();

            pFS->SetFilletGeomForSurfaces( pOrigGeom );
            SER( pOrigGeom->ReCalcFilletGeom( FALSE, &( sPtDer[1] ) ) );

            // Now the other FG
            pFS->SetFilletGeomForSurfaces( pNewGeom );
            SER( pNewGeom->ReCalcFilletGeom( FALSE, &( sPtDer[1] ) ) );

            pFS->SetFilletGeomForSurfaces( pBackupFG );


            // Now, the other FilletGeom, the one that didn't get split,
            // needs to be trimmed back to the far end of the split FG.
            SmFilletGeom *pOtherGeom = sGeoms[ 1-lWhichFilletGeom ];

            // Grap the correct rail
            ULONG lOtherRailIdx;
            SER( pOtherGeom->GetRailIndexOfVertex( pThisFV, lOtherRailIdx ) );
            SmFilletEdge *pOtherRail = pOtherGeom->GetRail( lOtherRailIdx );

            // The fillet vertex at which the other FG gets trimmed to
            // is still pointed to by pThisFV.  (See earlier comments.)
            SmPoint3d sTrimPt = pThisFV->GetPoint();

            SmExtent1d sOtherRailIvl = pOtherRail->GetInterval();

            double dOtherSplitParam = 0;
            SER( pOtherRail->GetCurve()->DropPoint(sOtherRailIvl,    // in : target curve allowed domain
                                              sTrimPt,          // in : Point to drop to curve
                                              NULL,             // in : For Snapping, NULL=NoSnap, Vec pointing to CrvInside from desired EndPt for snapping.
                                                                //      when DropPt is within dDistTol in cpVecToCrvInside dir of EndPt, Snap to Ivl EndPt.
                                                                //      Vec handles ambiguities and makes sure that curves are not just touching at the ends.
                                              SM_BIG_DOUBLE,    // in : Used for SnapDist, and SM_SO_INTERSECT sXSectTol3d distance.
                                                                //      sXSectTol3d = dDistTol ;      SM_SO_INTERSECT 
                                                                //      sXSectTol3d = SM_BIG_DOUBLE ; SM_SO_MINIMIZE, SM_SO_NORMALIZE
                                                                //      sSnapTol3d  = dDistTol ;      SM_SO_MINIMIZE, SM_SO_NORMALIZE, SM_SO_INTERSECT    
                                              pdGuessParam,     // in : If given, call cheap LocalPointSolve() else call expensive GlobalPointSolve()
                                              bSuccess,         // out: TRUE = found a drop point
                                              dOtherSplitParam, // out: found drop curve param
                                              dDist,            // out: found drop distance
                                              SM_SO_MINIMIZE)  ); // in : SM_SO_MINIMIZE = allow nonNormal drops near endPoints
                                                                //      SM_SO_INTERSECT= like Minimize but point must be within dDistTol
                                                                //      SM_SO_NORMALIZE= exclude nonNormal drops near endPoints
                                                                //      default:[SM_SO_MINIMIZE] to preserve original behavior

            if ( !bSuccess )
                { continue; }

            // Low or high end of other FG?
            SmFilletVertexuse *pOtherFVUse = apVUs[ 1-lWhichFilletGeom ];

            ULONG lOtherWhichEnd = 0;
            if ( pOtherGeom->GetRailVertexuse( lRailIdx, 1 ) == pOtherFVUse )
                lOtherWhichEnd = 1;

            // Trim sOtherRailIvl accordingly
            if ( lOtherWhichEnd == 0 )
            {
                SER( sOtherRailIvl.SetMinMax( dOtherSplitParam, sOtherRailIvl.GetMax() ) );
            }
            else
            {
                SER( sOtherRailIvl.SetMinMax( sOtherRailIvl.GetMin(), dOtherSplitParam ) );
            }

            // Finally, do the trim.
            SER( pOtherGeom->Trim( sOtherRailIvl ) ); // may snap sIvl by tol to existing knots

        }  // end loop over two FilletGeom's of this FilletVert.

        // Remove the opposite-rail FilletVertex from the list to be processed.
        // so that we don't try to split the other FilletGeom at it.
        ULONG lIdx;
        if (    pOppositeFV != NULL
             && apFVsToProcess.FindElement( pOppositeFV, lIdx ) )
        {
            apFVsToProcess.RemoveAt( lIdx );
            if ( lIdx <= lWhichFV )  // if we removed one before this one,
                lWhichFV--;          //   slide us back one.
        }
    }  // end loop over all FilletVert's.


    return SM_SUCCESS;
} // end SmBlendCorner::AdjustTwoBlendTangentCorner

/*******************************************************************//**
PURPOSE: Adjust the corner topology for the bevel cases where two fillets
    of the corner intersect the side edge at different locations.
    One additional edge will be added to corner to bound the fillet
    that hits farther away from the filleted vertex.

NOTES:
***********************************************************************/
SmStatus SmBlendCorner::AdjustTwoBlendTopology()
{
#ifdef SM_DEBUG_CODE
  if ( DebugLevel() > 0 ) 
    {
      smgfx_SetLook(2,4, 0,0,1); this->Draw(); sm_GraphicsLoop();
      sm_GraphicsLoop();
    }
#endif // SM_DEBUG_CODE

  // Can't do it if more than four total edges
  if(m_pAllEUs->GetSize() > 4)
    { SER( SM_ERR ) ; }

  // Check for mixed-convexity case: special routine for that.
  // For the 3x2 mixed-cvxty case, we create a 3-sided face.
  // We have created the edges and vertices for it,
  // so this is how we check for that case:

  ULONG lNumVertexEdges = m_vEdges.GetSize();
  if ( lNumVertexEdges == 3 )
    {   // Mixed-convexity corners
      SER( AdjustConcaveCornerTopology() );
      return SM_SUCCESS;
    }

  // Note: do this separately:
  // // Another special case is the Nx2 tangent case, with N > 3.
  // // In that case we might have to put another fillet surface
  // // in between the two that are there now.
  // if ( this->IsBlendingCorner()  &&  m_pAllEUs->GetSize() > 3 )
  //   {
  //     SER( AdjustTangentCorner() );
  //     return SM_SUCCESS;
  //   }

  // Examine the vertices that have been created for this corner.
  // Typically there will be two, one where the rails from the two
  // fillets intersect in the common face (between filleted edges),
  // and one where the other rail of each fillet intersects the
  // side edge (the un-filleted edge of this corner).
  // If there's only one side edge between their RAIL_X_EDGEUSE rails,
  // and the fillets are the same size, and things are generally square, 
  // then their rails will hit the side edge at the same point.
  // But, if these intersections are not the same point, then we have
  // to split the vertex into two, and insert a new edge to connect them.

  ULONG ii;
  for(ii=0;ii<m_vVertices.GetSize();ii++)
    {
      SmFilletVertex * pAdjustVert = m_vVertices[ii];

#ifdef SM_DEBUG_CODE
      if ( DebugLevel() > 0 ) 
        {
          smgfx_SetLook(1,8, 1,0,0); pAdjustVert->Draw(); sm_GraphicsLoop();
          sm_GraphicsLoop();
        }
#endif // SM_DEBUG_CODE

      // Check whether this is the vertex that corresponds to the
      // rail hitting the side edge, and not the rails hitting
      // each other in the common face.
      if(pAdjustVert->GetFilletVertexType() != SM_FV_RAIL_X_EDGEUSE)
        { continue ; }

      // Ok, this is the new fillet vertex corresponding to where the
      // rails from the two fillets are supposed to hit the side edge.
      // Check whether those rails intersect the side edge at the same place.

      SmTArray<SmFilletGeom*> sGeoms;
      pAdjustVert->GetFilletGeoms(sGeoms);
      if(sGeoms.GetSize() != 2)
        { continue ; }  // SER( SM_ERR );  maybe already adjusted...

      SmFilletGeom      * pFilletGeom1   = sGeoms[0] ;
      SmFilletGeom      * pFilletGeom2   = sGeoms[1] ;
      SmEdgeuse         * pSideEU1       = (SmEdgeuse*)pAdjustVert->GetPointClassObject() ;
      SmFilletVertexuse * pVU1           = pAdjustVert->GetVUAtRailEnd(pFilletGeom1) ;
      SmFilletVertexuse * pVU2           = pAdjustVert->GetVUAtRailEnd(pFilletGeom2, pVU1) ;
      double              dEdgeParam1    = pVU1->GetTsectPnt().m_dCurveParameter ;
      double              dEdgeParam2    = pVU2->GetTsectPnt().m_dCurveParameter ;
      SmPoint3d           sPnt1, sPnt2 ;  
      SmCurve           * pSideEdgeCurve = pSideEU1->GetEdge()->GetCurve() ;

      SER( pSideEdgeCurve->EvaluatePoint(dEdgeParam1,sPnt1)  );
      SER( pSideEdgeCurve->EvaluatePoint(dEdgeParam2,sPnt2)  );

      double dDist = sPnt1.DistanceBetween(sPnt2) ;

      // Check for the big-radius case, where we're rolling a large fillet
      // on the inside of a smaller-radius edge.
      // In that case, we'll create a 3x3 corner for the vertex.
      //
      // We'll also check for the Nx2 tangent case, with N > 3.
      if(this->IsBlendingCorner())     // Big-radius doesn't apply if Beveled.
        {
          SmBoolean           bConvertedToCorner;
          SmFilletExecutive * pExec = GetFilletExecutive();

          // This routine in the fillet exec will check for the
          // large-radius case, and convert the corner if so.
          SER( pExec->TestAndConvertBigRadFillet(pFilletGeom1, bConvertedToCorner ));
          if ( !bConvertedToCorner )
            {
              SER( pExec->TestAndConvertBigRadFillet( pFilletGeom2,
                      bConvertedToCorner ));
            }

          if ( bConvertedToCorner ) 
            { return SM_SUCCESS ; }

          if(dDist > 500.0 * m_dThisApproxTol3d)
            { MSG( _T("Big gap found. Will try to fix it later") ) ; }

          continue ;

        }  // End big-radius check.

      SM_ASSERT(m_bBevel == TRUE); // Bevel cases

      // If rails hit side edge at same point, no adjustment is needed.
      if ( dDist < m_dThisApproxTol3d )
        { continue ; }

      // Ok, if we're here, we're going to adjust the topology
      // of this corner.  We're going to make two separate vertices
      // where the rails from each of the two fillets hit the common
      // side edge, and a little edge between them.

      m_bTopologyAdjusted = TRUE;

      // Note on naming: these aren't necessarily smaller and larger radii.
      // For a given radius, the sharper the interior angle of the edge,
      // the farther away from the vertex the intersection will be.
      // But these names are descriptive, they make it easy to tell
      // which we're talking about: pSmallRadGeom means the one that
      // hits closer to the vertex.

      SmPoint3d      sVertPnt      = m_cpVertex->GetPoint() ;
      SmFilletGeom * pSmallRadGeom = NULL ;
      SmFilletGeom * pLargeRadGeom = NULL ;
      SmVertexuse  * pVUOnSideFace = NULL ;
      SmFace       * pSideFace     = NULL ;
      SmPoint3d      sPntOnSideFace ;

      // We'll split this vertex and keep the farther point (along the edge)
      // as the one SM_FV_RAIL_X_EDGEUSE.  The new vertex will be on
      // the side face and will be labeled as SM_FV_FILLET_X_FILLET
      if(sPnt1.DistanceBetween(sVertPnt) < sPnt2.DistanceBetween(sVertPnt))
        {
          pAdjustVert->SetPoint( sPnt2 ) ;
          pAdjustVert->SetOriginalTParam( dEdgeParam2 ) ;

          pSideFace      = pSideEU1->GetFace() ;
          pVUOnSideFace  = pVU1 ;
          pSmallRadGeom  = pFilletGeom1 ;
          pLargeRadGeom  = pFilletGeom2 ;
          sPntOnSideFace = sPnt1 ;
        } // end sPnt2 further away from VertPnt branch
      else // sPnt1 further away from VertPnt branch
        {
          pAdjustVert->SetPoint( sPnt1 ) ;
          pAdjustVert->SetOriginalTParam( dEdgeParam1 ) ;

          pSideFace      = pSideEU1->GetRadial()->GetFace() ;
          pVUOnSideFace  = pVU2 ;
          pSmallRadGeom  = pFilletGeom2 ;
          pLargeRadGeom  = pFilletGeom1 ;
          sPntOnSideFace = sPnt2 ;
        } // end sPnt1 further away from VertPnt branch

      // The vert that is farther from the filleted vertex, which is
      // the one where one rail hits the side edge, is now done.
      // Label it as 'processed'.
      pAdjustVert->SetStatus( SM_FIL_PROCESSED );

      SmFilletVertex * pVertOnSideEdge = pAdjustVert;

      // Create new fillet vertex on the 'side' face
      SmFilletVertex * pVertOnSideFace = new (m_pPseudoBrep) SmFilletVertex(this) ; NER( pVertOnSideFace );
      SmObjDelete sClean1( pVertOnSideFace );

      pVertOnSideFace->SetPoint( sPntOnSideFace );
      pVertOnSideFace->SetPointClass( SM_PC_FACE, pSideFace );
      pVertOnSideFace->SetFilletVertexType( SM_FV_FILLET_X_FILLET );

      SmBoolean        bInsertAfter = TRUE;
      ULONG            lIndex       = ii+1;
      SmFilletVertex * pOtherVert   = NULL;

      if(lIndex < m_vVertices.GetSize() )
        {
          pOtherVert = m_vVertices[lIndex];
        }
      else
        {
          // Insert new vertex before i
          lIndex--;
          bInsertAfter = FALSE;
          pOtherVert   = m_vVertices[ lIndex-1 ];
        }

      m_vVertices.InsertAt( lIndex, pVertOnSideFace );

      // Find the vertexuse of rail which is on the
      // 'side' face and reconnect it to the newly created vertex.
      SER( pVertOnSideEdge->Remove( pVUOnSideFace ) );
      SER(pVertOnSideFace->PostInsert( pVUOnSideFace ));

      // Recreate fillet edges
      SmFilletEdge * pFilletEdge1 = NULL;
      if ( bInsertAfter )
        {
          SER( SmFilletEdge::MakeFilletEdge( m_pPseudoBrep,   // in : target Brep to receive new topology objects
                                             pVertOnSideEdge, // in : start of new FilletEdge
                                             pVertOnSideFace, // in : end   of new FilletEdge
                                             pFilletEdge1,    // out: newly allocated FilletEdge - no Curve or UVTrimCurve data
                                             this ));         // in : NewFilletEdge's filletCorner, gets stored in rpNewFilletEdge->m_cpCorner
                                                              //      NULL to ignore, default:[NULL]
        }
      else
        {
          SER( SmFilletEdge::MakeFilletEdge( m_pPseudoBrep,   // in : target Brep to receive new topology objects
                                             pVertOnSideFace, // in : start of new FilletEdge
                                             pVertOnSideEdge, // in : end   of new FilletEdge
                                             pFilletEdge1,    // out: newly allocated FilletEdge - no Curve or UVTrimCurve data
                                             this ));         // in : NewFilletEdge's filletCorner, gets stored in rpNewFilletEdge->m_cpCorner
                                                              //      NULL to ignore, default:[NULL]
        }
      pFilletEdge1->SetFilletEdgeType( SM_FE_FILLET_X_SIDE_FACE ) ;

      SmFilletEdgeuse * pPrimEU = (SmFilletEdgeuse*)pFilletEdge1->GetPrimaryEdgeuse() ;
      pPrimEU->SetFilletGeom( pLargeRadGeom );
      pFilletEdge1->SetOriginalFace( pSideFace );

      SmFilletEdge * pFilletEdge2 = NULL;
      if ( bInsertAfter )
        {
          SER( SmFilletEdge::MakeFilletEdge( m_pPseudoBrep,    // in : target Brep to receive new topology objects
                                             pVertOnSideFace,  // in : start of new FilletEdge
                                             pOtherVert,       // in : end   of new FilletEdge
                                             pFilletEdge2,     // out: newly allocated FilletEdge - no Curve or UVTrimCurve data
                                             this ));          // in : NewFilletEdge's filletCorner, gets stored in rpNewFilletEdge->m_cpCorner
                                                               //      NULL to ignore, default:[NULL]
        }
      else 
        {
          SER( SmFilletEdge::MakeFilletEdge( m_pPseudoBrep,    // in : target Brep to receive new topology objects
                                             pOtherVert,       // in : start of new FilletEdge
                                             pVertOnSideFace,  // in : end   of new FilletEdge
                                             pFilletEdge2,     // out: newly allocated FilletEdge - no Curve or UVTrimCurve data
                                             this ));          // in : NewFilletEdge's filletCorner, gets stored in rpNewFilletEdge->m_cpCorner
                                                               //      NULL to ignore, default:[NULL]
        }
      pFilletEdge2->SetFilletEdgeType( SM_FE_FILLET_X_FILLET );

      pPrimEU = (SmFilletEdgeuse*)pFilletEdge2->GetPrimaryEdgeuse();
      SmFilletEdgeuse * pMateEU = (SmFilletEdgeuse*)pPrimEU->GetMate();

      pPrimEU->SetFilletGeom( pLargeRadGeom );
      pMateEU->SetFilletGeom( pSmallRadGeom );

      lIndex--;
      SmFilletEdge * pRemovedEdge = m_vEdges[lIndex];
      m_vEdges.RemoveAt( lIndex );
      SM_ASSERT( pRemovedEdge != NULL );
      if(pRemovedEdge) { delete pRemovedEdge ; pRemovedEdge = NULL; }

      if ( bInsertAfter )
        {
          m_vEdges.Add( pFilletEdge1 );
          m_vEdges.Add( pFilletEdge2 );
        }
      else
        {
          m_vEdges.Add( pFilletEdge2 );
          m_vEdges.Add( pFilletEdge1 );
        }
      sClean1.Clear();
      ii++;

    } // end iter ii, every m_vVertices vertex

  // all done
  return SM_SUCCESS;

} // end SmBlendCorner::AdjustTwoBlendTopology

/*******************************************************************//**
PURPOSE: Calculate geometry of a Nx2 concave corner patch where two
    adjacent rails may not meet at a single vertex (with type: SM_FV_RAIL_X_EDGEUSE).
    In that case, two patches shall be created:
    (1) a small patch to fill in the tiny gap &
    (2) a regular 3-sided COONS patch to blend two fillets.

NOTES: This corner should have either 3 or 6 boundary edges
***********************************************************************/
SmStatus SmBlendCorner::CreateTwoBlendPatch()
{
  ULONG ii, lTotalEdges = m_vEdges.GetSize();

  SmFilletEdge * apData1[3];
  SmFilletEdge * apData2[3];
  SmFilletVertex * apData3[2];
  SmTArray<SmFilletEdge*> sFirstThreeEdges(3,apData1);
  SmTArray<SmFilletEdge*> sLastThreeEdges(3,apData2);
  SmTArray<SmFilletVertex*> sLastTwoVerts(2,apData3);
  SmSurface * pEndPatch =  NULL;
  if (lTotalEdges > 3) 
    {
      for (ii=0; ii<3; ii++) 
        { sFirstThreeEdges.Add(m_vEdges[ii]) ; }

      for (ii=3; ii<m_vEdges.GetSize(); ii++) 
        { sLastThreeEdges.Add(m_vEdges[ii]) ; }

      if (lTotalEdges == 6) 
        {
          // Use the last 3 edges to create a small patch
          m_vEdges.ReSet();
          m_vEdges.Append(sLastThreeEdges);
          // Create a small corner face
          SER(SmFilletCorner::CreateEndPatch());

          // Temporarily remove the small patch
          pEndPatch = m_vSurfaces[0];
          m_vSurfaces.ReSet();
        }

      m_vEdges.ReSet();
      m_vEdges.Append(sFirstThreeEdges);

      // Temporarily remove the last two vertices
      for(ii=m_vVertices.GetSize()-1; ii>=3; ii--) 
        {
          sLastTwoVerts.Add(m_vVertices[ii]) ;
          m_vVertices.RemoveLast() ;
        }
    }

  // Create a 3-sided corner patch
  SER(SmFilletCorner::CreateCornerPatch());

  // Put the end patch back
  if (pEndPatch) 
    { m_vSurfaces.Add(pEndPatch) ; }

  if (sLastThreeEdges.GetSize() > 0) 
    {
      m_vEdges.Append(sLastThreeEdges) ;
      m_vVertices.Append(sLastTwoVerts) ;
    }

  return SM_SUCCESS;

} // end SmBlendCorner::CreateTwoBlendPatch
