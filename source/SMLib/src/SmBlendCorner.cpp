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
#include <SmFilletExact.h>

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

/*******************************************************************//**
PURPOSE: A side edge continued across its end pAt, where it is split only
    by a seam: the next edge's edgeuse on pEU's faceuse, or NULL.

NOTES: The rail on a side face F ends on the side edge E, between F and
    the end face; E ends at pAt. When pAt's only other edges are N, which
    continues E G1 at pAt between the same two faces, and a seam of F or
    of the end face (both its sides on that one face: a periodic
    surface's seam ending on E), pAt is no corner of the body: it splits
    one boundary of the two faces only because the seam ends there, and
    the rail that runs past it ends on N. A third edge that is no seam
    (a real edge between two faces) is a corner, and E is not continued.
***********************************************************************/
SmEdgeuse * sm_SideEdgeuseAcrossSeamVertex
  (SmEdgeuse      * pEU,         // in : the side edgeuse, on the rail's face
   const SmVertex * pAt,         // in : the side edge's end to continue it across
   double           dAngTolDeg)  // in : tangency tolerance (degrees)
{
  if ( pEU == NULL || pAt == NULL || pEU->GetRadial() == NULL )
    { return NULL; }
  SmEdge * pE = pEU->GetEdge();
  if (   pE == NULL || pE->IsClosed()
      || ( pE->GetStartVertex() != pAt && pE->GetEndVertex() != pAt ))
    { return NULL; }
  SmFace * pF = pEU->GetFace();
  SmFace * pG = pEU->GetRadial()->GetFace();
  if ( pF == NULL || pG == NULL || pF == pG )
    { return NULL; }

  SmTArray<SmEdge*> sEdges;
  pAt->GetEdges( sEdges );
  if ( sEdges.GetSize() != 3 )
    { return NULL; }
  SmEdge * pN = NULL;
  SmEdge * pSeam = NULL;
  for (ULONG ii=0; ii<3; ii++)
    {
      SmEdge * pX = sEdges[ii];
      if ( pX == pE )
        { continue; }
      SmTArray<SmFace*> sFaces;
      pX->GetFaces( sFaces );
      SmBoolean bOneFace = ( sFaces.GetSize() > 0 && ( sFaces[0] == pF || sFaces[0] == pG ) ) ? TRUE : FALSE;
      for (ULONG kk=1; kk<sFaces.GetSize(); kk++)
        { bOneFace = ( bOneFace && sFaces[kk] == sFaces[0] ) ? TRUE : FALSE; }
      if ( bOneFace && pX->IsSeam() )
        { pSeam = pX; }
      else
        { pN = pX; }
    }
  if ( pN == NULL || pSeam == NULL || pN->IsClosed() )
    { return NULL; }

  SmTArray<SmEdgeuse*> sEUs;
  if ( pN->GetEdgeuses( sEUs ) != SM_SUCCESS )
    { return NULL; }
  SmEdgeuse * pNext = NULL;
  for (ULONG ii=0; ii<sEUs.GetSize() && pNext == NULL; ii++)
    {
      SmEdgeuse * pX = sEUs[ii];
      if (   pX->GetFaceuse() == pEU->GetFaceuse()
          && pX->GetRadial() != NULL
          && pX->GetRadial()->GetFace() == pG )
        { pNext = pX; }
    }
  if ( pNext == NULL || ! sm_TestG1Continuity( pAt, pE, pN, dAngTolDeg, NULL, 0 ))
    { return NULL; }
  return pNext;

} // end sm_SideEdgeuseAcrossSeamVertex

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


// The point of crSurface nearest crPoint (on it), from the guess rUV: Newton on the foot's normal
// equations, to rounding. rbFound FALSE where it does not converge.
SmStatus sm_FootOnSurface(const SmSurface & crSurface, const SmPoint3d & crPoint, SmPoint2d & rUV, SmBoolean & rbFound)
{
  rbFound = FALSE;
  for ( ULONG lIter=0; lIter<40; lIter++ )
    {
      SmVector3d sD[2][2];
      SER( crSurface.Evaluate( rUV, 1, 1, TRUE, TRUE, TRUE, sD[0], FALSE ));
      const SmVector3d sR  = crPoint - sD[0][0];
      const SmVector3d sSu = sD[1][0], sSv = sD[0][1];
      const double dA = sSu.Dot( sSu ), dB = sSu.Dot( sSv ), dC = sSv.Dot( sSv );
      const double dDet = dA * dC - dB * dB;
      if ( dDet <= SM_EFF_ZERO * dA * dC )
        { return SM_SUCCESS; }
      const double dRu = sSu.Dot( sR ), dRv = sSv.Dot( sR );
      const double dDu = ( dC * dRu - dB * dRv ) / dDet, dDv = ( dA * dRv - dB * dRu ) / dDet;
      rUV.x += dDu;
      rUV.y += dDv;
      if ( smos_Fabs( dDu ) + smos_Fabs( dDv ) <= 4.0 * SM_EFF_ZERO * ( 1.0 + smos_Fabs( rUV.x ) + smos_Fabs( rUV.y )))
        {
          rbFound = TRUE;
          return SM_SUCCESS;
        }
    }
  return SM_SUCCESS;
}

// The product of two polynomials in Bernstein form on [0, 1], of degrees 2 (aQ) and 1 (aL): the
// cubic's Bernstein coefficients.
void sm_BernsteinProduct21(const double aQ[3], const double aL[2], double aP[4])
{
  aP[0] = aQ[0] * aL[0];
  aP[1] = ( aQ[0] * aL[1] + 2.0 * aQ[1] * aL[0] ) / 3.0;
  aP[2] = ( 2.0 * aQ[1] * aL[1] + aQ[2] * aL[0] ) / 3.0;
  aP[3] = aQ[2] * aL[1];
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
   m_eEnd(SM_BEND_ON_END_FACE),
   m_pEndCap(NULL),
   m_bG1Case(FALSE),
   m_bMitre(FALSE),
   m_pSideEU(NULL),
   m_pTurnSurface(NULL),
   m_lFillOrder(1),
   m_bConvexSetBack(TRUE),
   m_bSetBackOverflows(FALSE),
   m_pFill(NULL)
{
  m_apMixedSurface[0] = m_apMixedSurface[1] = NULL;
}

SmBlendCorner::~SmBlendCorner()
{
  delete m_pEndCap;
  DeleteSweptSurface();
  DeleteBallSurface();
  delete m_pFill;
  DeleteMixedSurfaces();
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
      // no radius builds a swept turn past a sharp concave side edge: its size is not the cause
      if ( pExec != NULL )
        { pExec->NoteFitRadius( 0.0, pEU1->GetEdge() ); }
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
SmStatus SmBlendCorner::CalcCornerVertGeom()
{
  SER( SmFilletCorner::CalcCornerVertGeom() );
  if ( IsMixed() )
    { return CalcMixedVertGeom(); }
  return IsBall() ? PlaceBallContacts() : SM_SUCCESS;
}

SmStatus SmBlendCorner::CalcCornerEdgeGeom()
{
  if ( m_eEvent == SM_BCE_END )
    { return CalcEndEdgeGeom(); }
  if ( IsMixed() )
    { return CalcMixedEdgeGeom(); }
  if ( IsMitreTurn() )
    {
      // the mitre between the rails' crossing and the side vertex, exact where the class allows;
      // its other edges (and a mitre of no exact class) as the edges' own kinds compute them
      for (ULONG ii=0; ii<m_vEdges.GetSize(); ii++)
        {
          SmFilletEdge * pEdge = m_vEdges[ii];
          if ( pEdge == NULL || pEdge->GetFilletEdgeType() != SM_FE_FILLET_X_FILLET || pEdge->GetCurve() != NULL )
            { continue; }
          SmBoolean bExact = FALSE;
          SER( ExactMitreCurve( pEdge, bExact ));
        }
      return SmFilletCorner::CalcCornerEdgeGeom();
    }
  if ( IsBall() )
    {
      SER( BuildBallSurface() );
      return SmFilletCorner::CalcCornerEdgeGeom();
    }
  if ( IsSetBack() )
    {
      SER( BuildSetBackFill() );
      return SmFilletCorner::CalcCornerEdgeGeom();
    }
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
PURPOSE: The mitre between the rails' crossing X and the side vertex S as
    the rounds' exact plane section, where they are two cylinders of one
    radius r whose axes meet (at C) and the plane they are symmetric about,
    through C square to the difference (or sum) of their axes, holds X and S.

NOTES: Cylinder 1, axis a through C, is C + r (cos t u + sin t w) + h a, u, w
   square to a. In the plane through C with normal m, h = -r (cos t u.m +
   sin t w.m) / a.m, so its section is C + r (cos t U + sin t W), U = u -
   (u.m / a.m) a, W = w - (w.m / a.m) a: the circle's affine image, a rational
   quadratic arc with the circle arc's weights (1, cos(t/2), 1) and its
   control points mapped. By the symmetry the plane's section of cylinder 2
   is the same curve. Its pcurves on the rounds are exact too (rational
   cubics, below). rbExact FALSE (the edge left as it was) where any of this does
   not hold: the rounds' certified intersection is the mitre then.
***********************************************************************/
SmStatus SmBlendCorner::ExactMitreCurve
  (SmFilletEdge * pMitre,   // in : the mitre edge, S-X or X-S
   SmBoolean    & rbExact)  // out: the exact curve was set
{
  rbExact = FALSE;
  NER( pMitre );
  SmFilletVertex * pStart = (SmFilletVertex*)pMitre->GetVertex();
  NER( pStart );
  SmFilletVertex * pEnd   = (SmFilletVertex*)pMitre->GetOtherVertex( pStart );
  NER( pEnd );
  const SmBoolean bFromS = ( pStart->GetFilletVertexType() == SM_FV_RAIL_X_EDGEUSE ) ? TRUE : FALSE;
  SmFilletVertex * pS = bFromS ? pStart : pEnd;
  SmFilletVertex * pX = bFromS ? pEnd : pStart;
  if (   pS->GetFilletVertexType() != SM_FV_RAIL_X_EDGEUSE || pX->GetFilletVertexType() != SM_FV_RAIL_X_RAIL
      || !pS->IsProcessed() || !pX->IsProcessed() )
    return SM_SUCCESS;
  SmFilletEdgeuse * apEU[2];
  apEU[0] = (SmFilletEdgeuse*)pMitre->GetPrimaryEdgeuse();
  NER( apEU[0] );
  apEU[1] = (SmFilletEdgeuse*)apEU[0]->GetMate();
  NER( apEU[1] );
  const double dTol = GetThisApproxTol3d();

  // both rounds constant circular ones on straight spines (cylinders), of one radius; the balls'
  // centres at X and the spines' directions
  SmPoint3d  asCentre[2];
  SmVector3d asAxis[2];
  double     adRadius[2];
  for (ULONG ii=0; ii<2; ii++)
    {
      SmFilletGeom * pFG = apEU[ii]->GetFilletGeom();
      NER( pFG );
      SmFilletSolver * pFS = pFG->GetFilletSolver();
      NER( pFS );
      SmFilletVertexuse * pVU = pX->GetVUAtRailEnd( pFG );
      NER( pVU );
      SmTsectPnt   sTsect = pVU->GetTsectPnt();
      SmExactBlend sBlend;
      SER( pFS->ClassifyExactBlend( sTsect, sBlend ));
      // Open gap: a round on a circular spine (a torus) has no exact mitre here yet; the mitre
      // is then the rounds' certified SSI curve, and since the fillet's own trim curves reach
      // the result (SmFilletGeom::MakeFaceBrep, SmFace::ShrinkGeometry) its fitted trim curves
      // are the result's, not projections.
      if ( sBlend.m_eSpine != SM_BS_LINE || pFS->GetFilletSurfaceGenerator() == NULL
          || pFS->GetFilletSurfaceGenerator()->GetFilletSurfaceGeneratorType() != SM_FSG_CIRCULAR )
        return SM_SUCCESS;
      SmPoint3d sOn[2];
      for (ULONG k=0; k<2; k++)
        {
          SmOffsetSurface * pOffset = pFG->GetOffsetSurface( k );
          NER( pOffset );
          SER( pOffset->EvaluatePoint( pVU->GetTsectPnt().UVPos( k ), sOn[k] ));
        }
      if ( sOn[0].DistanceBetween( sOn[1] ) > dTol )
        return SM_SUCCESS;
      asCentre[ii] = sOn[0];
      adRadius[ii] = sBlend.m_dRadius;
      SmBSplineCurve * pSpine = pFG->GetCenterLineCurve();
      NER( pSpine );
      SmPoint3d sLinePt;
      if ( !pSpine->IsLine( 10, dTol, sLinePt, asAxis[ii] ) || asAxis[ii].Unitize( FALSE ) != SM_SUCCESS )
        return SM_SUCCESS;
    }
  // Open gap: two cylinders of unequal radii meet in a quartic, not a plane conic, and have no
  // exact mitre here yet: the certified SSI curve, its fitted trim curves the result's (above).
  const double dR = adRadius[0];
  if ( smos_Fabs( adRadius[1] - dR ) > dTol || dR <= dTol )
    return SM_SUCCESS;

  // C: where the axes meet
  const SmVector3d sD  = asCentre[1] - asCentre[0];
  const double     dAB = asAxis[0].Dot( asAxis[1] );
  const double     dDen = 1.0 - dAB * dAB;
  if ( dDen <= SM_EFF_ZERO )
    return SM_SUCCESS; // parallel axes: no turn
  const double dS1 = ( sD.Dot( asAxis[0] ) - dAB * sD.Dot( asAxis[1] )) / dDen;
  const double dS2 = ( dAB * sD.Dot( asAxis[0] ) - sD.Dot( asAxis[1] )) / dDen;
  const SmPoint3d sC  = asCentre[0] + asAxis[0] * dS1;
  const SmPoint3d sC2 = asCentre[1] + asAxis[1] * dS2;
  if ( sC.DistanceBetween( sC2 ) > dTol )
    return SM_SUCCESS;

  // X and S exactly: on the planes the rounds roll on, X is the foot of C on the common face
  // (where the rails there cross: the corner ball touches it there) and S where each round's rail
  // on its other face (the line through C's foot on that face along the axis) meets the side
  // edge's line. The rails are trimmed at the corner's vertices (SmFilletGeom::TrimRailCurves),
  // so the rails, the mitre and the faces between them all end there.
  SmEdgeuse * pSideEU = (SmEdgeuse*)pS->GetPointClassObject();
  NER( pSideEU );
  SmFace * apWall[2] = { pSideEU->GetFace(), pSideEU->GetRadial()->GetFace() };
  SmFace * pCommon   = (SmFace*)pX->GetPointClassObject();
  NER( pCommon );
  SmBlendQuadric sCommonPlane, asWallPlane[2];
  SER( SmBlendQuadric::FromSurface( pCommon->GetSurface(), pCommon->GetUVDomain().GetMid(), dTol, sCommonPlane ));
  for (ULONG ii=0; ii<2; ii++)
    {
      NER( apWall[ii] );
      SER( SmBlendQuadric::FromSurface( apWall[ii]->GetSurface(), apWall[ii]->GetUVDomain().GetMid(), dTol, asWallPlane[ii] ));
      if ( asWallPlane[ii].m_eKind != SmBlendQuadric::SM_BQ_PLANE )
        return SM_SUCCESS;
    }
  SmPoint3d  sEdgePt;
  SmVector3d sEdgeDir;
  NER( pSideEU->GetEdge()->GetCurve() );
  if (   sCommonPlane.m_eKind != SmBlendQuadric::SM_BQ_PLANE
      || !pSideEU->GetEdge()->GetCurve()->IsLine( 10, dTol, sEdgePt, sEdgeDir ) || sEdgeDir.Unitize( FALSE ) != SM_SUCCESS )
    return SM_SUCCESS;
  auto foot = []( const SmPoint3d & crP, const SmBlendQuadric & crPlane )
    { return crP - crPlane.m_sDir * SmVector3d( crP - crPlane.m_sPoint ).Dot( crPlane.m_sDir ); };
  const SmPoint3d sXPnt = foot( sC, sCommonPlane );
  SmPoint3d asOnEdge[2];
  for (ULONG ii=0; ii<2; ii++)
    {
      // the rail on wall ii, Q + s a, meets the edge's line, P + t e
      const SmPoint3d  sQ  = foot( sC, asWallPlane[ii] );
      const SmVector3d sQP = sQ - sEdgePt;
      const double     dAE = asAxis[ii].Dot( sEdgeDir );
      const double     dD2 = 1.0 - dAE * dAE;
      if ( dD2 <= SM_EFF_ZERO )
        return SM_SUCCESS;
      const double dT = ( sQP.Dot( sEdgeDir ) - dAE * sQP.Dot( asAxis[ii] )) / dD2;
      const double dS = ( dAE * sQP.Dot( sEdgeDir ) - sQP.Dot( asAxis[ii] )) / dD2;
      asOnEdge[ii] = sEdgePt + sEdgeDir * dT;
      if ( asOnEdge[ii].DistanceBetween( sQ + asAxis[ii] * dS ) > dTol )
        return SM_SUCCESS; // the rail misses the edge
    }
  const SmPoint3d sSPnt = asOnEdge[0];
  if (   asOnEdge[1].DistanceBetween( sSPnt ) > dTol
      || sXPnt.DistanceBetween( pX->GetPoint() ) > dTol || sSPnt.DistanceBetween( pS->GetPoint() ) > dTol )
    return SM_SUCCESS;

  // the plane of symmetry holding X and S
  SmVector3d sM;
  SmBoolean  bPlane = FALSE;
  for (int iSign=-1; iSign<=1 && !bPlane; iSign+=2)
    {
      SmVector3d sTry = asAxis[0] + asAxis[1] * static_cast<double>( iSign );
      if ( sTry.Unitize( FALSE ) != SM_SUCCESS )
        { continue; }
      if ( smos_Fabs( SmVector3d( sXPnt - sC ).Dot( sTry )) <= dTol && smos_Fabs( SmVector3d( sSPnt - sC ).Dot( sTry )) <= dTol )
        { sM = sTry; bPlane = TRUE; }
    }
  if ( !bPlane )
    return SM_SUCCESS;
  const SmVector3d & sA = asAxis[0];
  const double dAM = sA.Dot( sM );
  if ( smos_Fabs( dAM ) <= SM_EFF_ZERO )
    return SM_SUCCESS;

  // the arc on cylinder 1 from X (t = 0) to S
  SmVector3d sXOff = sXPnt - sC;
  SmVector3d sU    = sXOff - sA * sXOff.Dot( sA );
  if ( smos_Fabs( sU.Length() - dR ) > dTol || sU.Unitize( FALSE ) != SM_SUCCESS )
    return SM_SUCCESS;
  const SmVector3d sW    = sA * sU;
  const SmVector3d sSOff = sSPnt - sC;
  const double     dT    = atan2( sSOff.Dot( sW ), sSOff.Dot( sU ));
  if ( smos_Fabs( dT ) <= SM_EFF_ZERO || smos_Fabs( dT ) >= 0.99 * SM_PI )
    return SM_SUCCESS;
  const SmVector3d sUm = sU - sA * ( sU.Dot( sM ) / dAM );
  const SmVector3d sWm = sW - sA * ( sW.Dot( sM ) / dAM );
  const double adX[3] = { 1.0, 1.0, cos( dT ) };
  const double adY[3] = { 0.0, tan( 0.5 * dT ), sin( dT ) };
  SmTArray<SmPoint3d> sPoints;
  SmTArray<double>    sWeights;
  for (ULONG jj=0; jj<3; jj++)
    {
      const ULONG k = bFromS ? 2 - jj : jj; // the curve runs as the edge does
      sPoints.Add( sC + ( sUm * adX[k] + sWm * adY[k] ) * dR );
      sWeights.Add( k == 1 ? cos( 0.5 * dT ) : 1.0 );
    }
  SmTArray<ULONG>  sMults;
  SmTArray<double> sKnots;
  sMults.Add( 3 ); sMults.Add( 3 );
  sKnots.Add( 0.0 ); sKnots.Add( 1.0 );
  SmBSplineCurve * pCurve = NULL;
  SER( SmBSplineCurve::CreateCanonical( m_crContext, 3, 2, sPoints, SM_CF_ELLIPTIC_ARC, sMults, sKnots, SM_KT_UNSPECIFIED,
                                        &sWeights, NULL, pCurve ));
  NER( pCurve );
  SmObjDelete sCleanCurve( pCurve );

  // Its pcurves on the two rounds, exact. A round on a straight spine is its section's rational
  // arc moved along the spine: along it the parameter is affine in the point's distance along the
  // axis, so on the mitre, where the point is the curve's rational quadratic C(s) = P(s) / W(s),
  // that parameter is N(s) / W(s), N quadratic; across it the arc's parameter p and the mitre's s
  // both run as the tangent of the half angle about the axis (each its arc's, linearly), so p is
  // linear fractional in s, (a s + b) / (c s + d). The pcurve is the rational cubic with
  // denominator W(s) (c s + d). Both are read off the round's own parameters at s = 0, 1/2, 1 (its
  // feet there, to rounding), and the pcurve holds where the round's surface at it is the mitre at
  // s = 1/4, 3/4 and 1/8, 7/8 (to 1e-12 of the size): else the round is not so parameterized and the
  // mitre is its certified intersection.
  SmBSplineCurve * apUV[2] = { NULL, NULL };
  SmObjDelete asCleanUV[2];
  const double adW[3] = { sWeights[0], sWeights[1], sWeights[2] };
  const double dSize  = 1.0 + sXPnt.GetMaxDimension() + dR;
  for (ULONG ii=0; ii<2; ii++)
    {
      SmSurface * pSurface = apEU[ii]->GetFilletGeom()->GetFilletSurface();
      NER( pSurface );
      // a start for the feet: the round's own (fitted) drop of the mitre
      double dMaxDrop = 0.0, dDev = 0.0;
      SmTArray<SmBSplineCurve*> sUVs;
      SER( pSurface->DropCurve( m_crContext, pSurface->GetNaturalUVDomain(), *pCurve, pCurve->GetNaturalInterval(), dTol,
                                dMaxDrop, dDev, sUVs, FALSE ));
      SmObjsDelete<SmBSplineCurve*> sCleanDrops( &sUVs );
      if ( sUVs.GetSize() != 1 || dMaxDrop > dTol )
        { return SM_SUCCESS; } // not on that round: its certified intersection
      const double adAt[7] = { 0.0, 0.5, 1.0, 0.25, 0.75, 0.125, 0.875 };
      SmPoint2d asUV[7];
      SmPoint3d asOn[7];
      for (ULONG k=0; k<7; k++)
        {
          SmPoint3d sGuess;
          SER( sUVs[0]->EvaluatePoint( sUVs[0]->GetNaturalInterval().Evaluate( adAt[k] ), sGuess ));
          SER( pCurve->EvaluatePoint( pCurve->GetNaturalInterval().Evaluate( adAt[k] ), asOn[k] ));
          asUV[k] = SmPoint2d( sGuess.x, sGuess.y );
          SmBoolean bFound = FALSE;
          SER( sm_FootOnSurface( *pSurface, asOn[k], asUV[k], bFound ));
          if ( !bFound )
            { return SM_SUCCESS; }
        }
      // which parameter runs along the spine (N / W) and which across it (linear fractional)
      SmBSplineCurve * pExact = NULL;
      for (ULONG lAcross=0; lAcross<2 && pExact == NULL; lAcross++)
        {
          auto coord = [&]( ULONG k, ULONG lDir ) { return lDir == 0 ? asUV[k].x : asUV[k].y; };
          const double dP0 = coord( 0, lAcross ), dPh = coord( 1, lAcross ), dP1 = coord( 2, lAcross );
          if ( smos_Fabs( dPh - dP1 ) <= SM_EFF_ZERO )
            { continue; }
          const double dG  = ( dP1 + dP0 - 2.0 * dPh ) / ( dPh - dP1 );   // p = (a s + b) / (g s + 1)
          const double adL[2]   = { 1.0, 1.0 + dG };
          const double adNum[2] = { dP0, dP1 * ( 1.0 + dG ) };
          if ( adL[1] <= 0.0 )
            { continue; }
          const ULONG lAlong = 1 - lAcross;
          const double dU0 = coord( 0, lAlong ), dUh = coord( 1, lAlong ), dU1 = coord( 2, lAlong );
          const double dWh = 0.25 * ( adW[0] + 2.0 * adW[1] + adW[2] );
          const double adN[3] = { dU0 * adW[0], 0.5 * ( 4.0 * dUh * dWh - dU0 * adW[0] - dU1 * adW[2] ), dU1 * adW[2] };
          double adDen[4], adAlong[4], adAcross[4];
          sm_BernsteinProduct21( adW, adL, adDen );
          sm_BernsteinProduct21( adN, adL, adAlong );
          sm_BernsteinProduct21( adW, adNum, adAcross );
          SmTArray<SmPoint3d> sUVPoints;
          SmTArray<double>    sUVWeights;
          for (ULONG j=0; j<4; j++)
            {
              const double dA = adAlong[j] / adDen[j], dC = adAcross[j] / adDen[j];
              sUVPoints.Add( lAcross == 1 ? SmPoint3d( dA, dC, 0.0 ) : SmPoint3d( dC, dA, 0.0 ));
              sUVWeights.Add( adDen[j] );
            }
          SmTArray<ULONG>  sUVMults;
          SmTArray<double> sUVKnots;
          sUVMults.Add( 4 ); sUVMults.Add( 4 );
          sUVKnots.Add( pCurve->GetNaturalInterval().GetMin() ); sUVKnots.Add( pCurve->GetNaturalInterval().GetMax() );
          SmBSplineCurve * pTry = NULL;
          SER( SmBSplineCurve::CreateCanonical( m_crContext, 2, 3, sUVPoints, SM_CF_UNSPECIFIED, sUVMults, sUVKnots, SM_KT_UNSPECIFIED,
                                                &sUVWeights, NULL, pTry ));
          NER( pTry );
          SmBoolean bHolds = TRUE;
          for (ULONG k=0; k<7 && bHolds; k++)
            {
              SmPoint3d sUV, sOnRound;
              SER( pTry->EvaluatePoint( pTry->GetNaturalInterval().Evaluate( adAt[k] ), sUV ));
              SER( pSurface->EvaluatePoint( SmPoint2d( sUV.x, sUV.y ), sOnRound ));
              bHolds = ( sOnRound.DistanceBetween( asOn[k] ) <= 1e-12 * dSize ) ? TRUE : FALSE;
            }
          if ( bHolds )
            { pExact = pTry; }
          else
            { delete pTry; }
        }
      if ( pExact == NULL )
        { return SM_SUCCESS; } // the round is not parameterized so: its certified intersection
      apUV[ii] = pExact;
      asCleanUV[ii].SetObj( pExact );
    }

  pMitre->SetCurve( pCurve, TRUE );
  pCurve->SetOwner( pMitre );
  sCleanCurve.Clear();
  pX->SetPoint( sXPnt );
  pS->SetPoint( sSPnt );
  {
    SmCurve * pSideCurve = pSideEU->GetEdge()->GetCurve();
    SmBoolean bFound = FALSE;
    double    dParam = 0.0, dDist = 0.0;
    SER( pSideCurve->DropPoint( pSideEU->GetEdge()->GetInterval(), sSPnt, NULL, dTol, NULL, bFound, dParam, dDist ));
    if ( bFound )
      { pS->SetOriginalTParam( dParam ); }
  }
  for (ULONG ii=0; ii<2; ii++)
    {
      apEU[ii]->SetUVCurve( apUV[ii] );
      asCleanUV[ii].Clear();
      apEU[ii]->GetFilletGeom()->AddSideFilletEdgeuse( apEU[ii] );
    }
  pMitre->SetFilletEdgeType( SM_FE_PRECOMPUTED );
  rbExact = TRUE;
  return SM_SUCCESS;

} // end SmBlendCorner::ExactMitreCurve

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
  if ( IsBall() )
    { return MakeBallTopology(); }
  if ( IsMixed() )
    { return MakeMixedTopology(); }
  if ( m_eEvent == SM_BCE_JOIN )
    { return MakeJoinTopology(); }
  if ( m_eEvent == SM_BCE_END )
    {
      return   ( m_eEnd == SM_BEND_OPEN )      ? MakeEndOpen()
             : ( m_eEnd == SM_BEND_RUNS_OUT )  ? MakeEndRunningOut()
             : ( m_eEnd == SM_BEND_AT_SMOOTH ) ? MakeEndAtSmoothVertex()
             :                                   MakeEndOnEndFace();
    }
  if ( IsSetBack() )
    { return MakeSetBackTopology(); }
  return IsSwept() ? MakeSweptTopology() : IsMitreTurn() ? MakeMitreTopology() : MakeTwoBlendTopology();
}

SmStatus SmBlendCorner::AdjustTopology()
{
  if ( IsBall() || IsSetBack() || IsMixed() || m_eEvent == SM_BCE_JOIN || m_eEvent == SM_BCE_END )
    { return SM_SUCCESS; }
  return IsSwept() ? SM_SUCCESS : IsMitreTurn() ? AdjustMitreTopology() : AdjustTwoBlendTopology();
}

SmStatus SmBlendCorner::CreateCornerPatch()
{
  if ( IsBall() )
    { return CreateBallPatch(); }
  if ( IsMixed() )
    { return CreateMixedPatches(); }
  if ( m_eEvent == SM_BCE_END )
    { return CreateEndPatchOfEnd(); }
  if ( m_eEvent == SM_BCE_JOIN )
    { return SM_SUCCESS; }
  if ( IsSetBack() )
    { return CreateSetBackPatches(); }
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

  // only a Turn at a corner of four edges comes here (SmBlendCorner::MakeCornerTopology): a Turn
  // of three is MakeMitreTopology's or MakeSweptTopology's, a G1 pair the Join's
  if ( lTotalEUs != 4)
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

  // The two blends are mitred: the pair is not G1, nor within 5 degrees at a face corner of the
  // other convexity (SmBlendCorner::Classify sends those to the Join).
  m_bBevel = TRUE;
  m_bBlending = FALSE;

  // One fillet (cross) edge at this corner, between the two side vertices: the mitre.
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

  pNewFilletCrossEdge->SetFilletEdgeType( SM_FE_FILLET_X_FILLET );

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
  // a closed Join's one blend continues itself
  if ( IsClosedJoin() )
    { return ( m_pExecutive->GetFilletSolverOfEdgeuse( (*m_pSolverEUs)[0] ) == pFS ) ? SM_CONST_CAST( SmFilletSolver*, pFS ) : NULL; }
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
