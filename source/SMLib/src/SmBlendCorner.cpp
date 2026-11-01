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
   SmBlendTurn            eTurn)
 : SmFilletCorner(crContext, cpVertex, pSolverEUs, pAllEUs, pExec, dApproxTol3d, dTangencyTolRadians),
   m_eTurn(eTurn),
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
  if ( crSolverEUs.GetSize() != 2 || crAllEUs.GetSize() != 3 )
    { return SM_SUCCESS; }
  SmEdgeuse * pEU1 = crSolverEUs[0];
  SmEdgeuse * pEU2 = crSolverEUs[1];
  NER( pEU1 );
  NER( pEU2 );
  if ( sm_TestG1Continuity( cpVertex, pEU1->GetEdge(), pEU2->GetEdge(), dTangencyTolDeg, NULL, 0 ))
    { return SM_SUCCESS; } // a Join: the blends continue each other
  SmEdgeuse * pRingEU = NULL;
  for (ULONG ii=0; ii<crAllEUs.GetSize(); ii++)
    {
      if ( crAllEUs[ii] != pEU1 && crAllEUs[ii] != pEU2 )
        { pRingEU = crAllEUs[ii]; }
    }
  NER( pRingEU );
  SmEdgeuse * pSideEU = sm_FindSideEdgeuse( pEU1, pRingEU, 0 );
  if ( pSideEU == NULL || sm_GetFaceByTwoEdgeuses( pEU1, pEU2, pSideEU ) == NULL )
    { return SM_SUCCESS; }
  reEvent = SM_BCE_TURN;
  if ( !pExec->GetSweptTurns() )
    { return SM_SUCCESS; }
  reTurn = SM_BT_SWEPT;

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
  const double dSinTol = sin( SM_DEG2RAD( dTangencyTolDeg ));
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
SmStatus SmBlendCorner::MakeCornerTopology()
{
  NER( m_pSolverEUs );
  NER( m_pAllEUs );
  if ( m_eTurn != SM_BT_SWEPT || m_pSolverEUs->GetSize() != 2 || m_pAllEUs->GetSize() != 3 )
    { SER( SM_ERR ); } // a Mitre turn is still SmFilletNx2Corner's

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

} // end SmBlendCorner::MakeCornerTopology

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
SmStatus SmBlendCorner::CreateCornerPatch()
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

} // end SmBlendCorner::CreateCornerPatch
