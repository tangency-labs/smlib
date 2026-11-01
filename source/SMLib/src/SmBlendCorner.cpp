// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

/**********************************************************************//**
* FILE NAME --- SmBlendCorner.cpp
* PURPOSE: SolidCraft change. The blend corner (SmBlendCorner.h): first,
*   the corner helpers it shares with the corner hierarchy.
**********************************************************************/

#include "StdAfx.h"

#ifndef __SMBLENDCORNER_H__
#include <SmBlendCorner.h>
#endif

#ifndef __SMFILLETSOLVER_H__
#include <SmFilletSolver.h>
#endif


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
