// SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
// SPDX-License-Identifier: Apache-2.0

/**********************************************************************//**
* FILE NAME --- SmFilletStandardSolver.cpp
* PURPOSE: Source code file for those standard SmFilletSolver objects.
*   e.g. SmConstantRadiusFS, SmConstantDistanceFS, SmVariableRadiusFS &
*        SmSurfaceSurfaceFS
**********************************************************************/

#include "StdAfx.h"

#ifndef __SMFILLETSTANDARDSOLVER_H__
#include <SmFilletStandardSolver.h>
#endif

#ifndef __SMFILLETEXECUTIVE_H__
#include <SmFilletExecutive.h>
#endif

#ifndef __SMFILLETINTERSECTOR_H__
#include <SmFilletIntersector.h>
#endif

#include <SmBSplineCurve.h>
#include <SmGraphicsExtern.h>
#include <SmGeomUtility.h>
#include <SmSurfOfExtrusion.h>

/*******************************************************************//**
PURPOSE: SolidCraft change: the side of each face of the sector between
    rpEdgeuses[0] and its radial on which the other face lies, at the
    middle of the edge: rdSides[i] = +1 when face (1-i) lies on the side
    face i's surface normal points to, -1 when it lies on the other side:
    the sign of face i's offset (SmOffsetSurface offsets along the surface
    normal).

NOTES: The faces leave the edge along their binormals b0 and b1. When b1
    is not tangent to face 0 the side is that of b1 (and the same for
    face 1). Where the faces are tangent at the edge and their binormals
    point the same way (a cusp: a knife edge, or its complement) they
    part to second order: along the binormal b, face j lies at height
    s^2/2 * kj * Nj (kj its normal curvature in direction b, Nj its unit
    surface normal), so face (1-i) lies on the side of face i given by
    k(1-i) N(1-i).ni - ki Ni.ni (ni face i's normal at the edge). Faces
    tangent at the edge whose binormals point apart continue each other
    (a smooth edge), and faces of a cusp
    whose curvatures match there do not part to second order: neither
    has a sector to fillet (SM_ERR). dLength is the size of the blend,
    the scale the parting of the faces is measured at.
***********************************************************************/
static SmStatus sm_SectorSides
  (SmEdgeuse * const rpEdgeuses[2],  // in : the edgeuse bounding the sector, and its radial
   double            dLength,        // in : the blend's size
   double            rdSides[2])     // out: the side of each face the other lies on
{
  NER(rpEdgeuses[0]); NER(rpEdgeuses[1]);
  SmEdge *pEdge = rpEdgeuses[0]->GetEdge(); NER(pEdge);
  double dParam = pEdge->GetInterval().Evaluate(0.5);

  SmPoint3d  sPnt[2];
  SmVector3d sBinVec[2], sNormal[2];
  SmPoint2d  sUV[2];
  SmVector2d sUVDir[2];
  for (ULONG i = 0; i < 2; i++)
    {
      SER(rpEdgeuses[i]->EvaluateBinormal(dParam, FALSE, sPnt[i], sBinVec[i], NULL, &sNormal[i], &sUV[i], &sUVDir[i]));
      if (rpEdgeuses[i]->GetLoopuse()->GetFaceuse()->GetOrientation() == SM_OT_OPPOSITE)
        { sNormal[i] = -sNormal[i]; } // the surface's normal
    }

  double dFirst[2] = { sBinVec[1].Dot(sNormal[0]), sBinVec[0].Dot(sNormal[1]) };
  if (smos_Fabs(dFirst[0]) >= SM_EFF_ZERO_SQRT && smos_Fabs(dFirst[1]) >= SM_EFF_ZERO_SQRT)
    {
      rdSides[0] = dFirst[0] < 0.0 ? -1.0 : 1.0;
      rdSides[1] = dFirst[1] < 0.0 ? -1.0 : 1.0;
      return SM_SUCCESS;
    }
  if (sBinVec[0].Dot(sBinVec[1]) <= 0.0)
    { return SM_ERR; } // smooth: the faces continue each other

  // a cusp: each face's normal curvature along its binormal, into the face
  double     dCurvature[2];
  SmVector3d sSurfaceNormal[2];
  for (ULONG i = 0; i < 2; i++)
    {
      const SmSurface *cpSurface = rpEdgeuses[i]->GetLoopuse()->GetFaceuse()->GetFace()->GetSurface();
      NER(cpSurface);
      SmVector3d sTangent;
      SER(cpSurface->EvaluateNormalSection(sUV[i], sUVDir[i].x >= 0.0, sUVDir[i].y >= 0.0, sBinVec[i],
                                           sTangent, dCurvature[i], sSurfaceNormal[i]));
    }
  for (ULONG i = 0; i < 2; i++)
    {
      ULONG j = 1 - i;
      double dPart =   dCurvature[j] * sSurfaceNormal[j].Dot(sNormal[i])
                     - dCurvature[i] * sSurfaceNormal[i].Dot(sNormal[i]);
      if (smos_Fabs(dPart) * dLength < SM_EFF_ZERO_SQRT)
        { return SM_ERR; } // the faces osculate: no sector to fillet
      rdSides[i] = dPart < 0.0 ? -1.0 : 1.0;
    }
  return SM_SUCCESS;

} // end sm_SectorSides

/*******************************************************************//**
PURPOSE: Constructor for Constant Radius Fillet object that allows
    different offsets for corresponding surfaces.

NOTES:
***********************************************************************/
SmConstantRadiusFS::SmConstantRadiusFS
  (const SmContext & crContext,
   double dThisApproxTol3d,
   double dAngleTolerance,
   double dTangencyTolerance,
   double dOffsetRadiusSurface1,
   double dOffsetRadiusSurface2,
   const SmSurface & crSurface1,
   const SmSurface & crSurface2,
   SmBoolean bOrientationSurface1,
   SmBoolean bOrientationSurface2)
 : SmFilletSolver(crContext,
                  dThisApproxTol3d,
                  dAngleTolerance,
                  dTangencyTolerance)
{
    m_dFilletRadii[0] = dOffsetRadiusSurface1;
    m_dFilletRadii[1] = dOffsetRadiusSurface2;

    m_dOrientations[0] = 1.0;
    if (bOrientationSurface1 == FALSE) m_dOrientations[0] = - 1.0;
    m_dOrientations[1] = 1.0;
    if (bOrientationSurface2 == FALSE) m_dOrientations[1] = - 1.0;

    m_pSurfaces[0] = new (crContext) SmOffsetSurface( m_dFilletRadii[0]*m_dOrientations[0],
                                                      SM_CONST_CAST(SmSurface &,crSurface1));
    m_pSurfaces[1] = new (crContext) SmOffsetSurface( m_dFilletRadii[1]*m_dOrientations[1],
                                                      SM_CONST_CAST(SmSurface &,crSurface2));
    m_pSurfaces[0]->SetAllowExtension(TRUE);
    m_pSurfaces[1]->SetAllowExtension(TRUE);

    (void)FinishConstruction();  // Take care of other things in the base class. // SolidCraft: a constructor cannot report

} // end SmConstantRadiusFS::SmConstantRadiusFS constructor

/*******************************************************************//**
PURPOSE: Constructor for Constant Radius Fillet object.

NOTES:
***********************************************************************/
SmConstantRadiusFS::SmConstantRadiusFS
  (const SmContext & crContext,
   double            dThisApproxTol3d,
   double            dAngleTolerance,
   double            dTangencyTolerance,
   double            dFilletRadius,
   SmEdgeuse       * pEdgeuse,
   double          * pdOptSecondRadius)
  : SmFilletSolver(crContext,
                   dThisApproxTol3d,
                   dAngleTolerance,
                   dTangencyTolerance)
{
    if (!pEdgeuse) SE(SM_ERR);
    m_dFilletRadii[0] = dFilletRadius;
    m_dFilletRadii[1] = dFilletRadius;
    if (pdOptSecondRadius) {
        m_dFilletRadii[1] = *pdOptSecondRadius;
    }
    m_pEdgeuses[0] = pEdgeuse;
    m_pEdgeuses[1] = pEdgeuse->GetRadial();

    // Create offset surfaces of appropriate radius in appropriate direction.

    // SolidCraft change: offset each surface to the side of it the other face lies on
    // (sm_SectorSides), a cusp's faces too.
    double dSides[2] = { 1.0, 1.0 };
    if (sm_SectorSides(m_pEdgeuses, smos_Max(smos_Fabs(m_dFilletRadii[0]), smos_Fabs(m_dFilletRadii[1])), dSides) != SM_SUCCESS)
    {
        // Smooth, or no sector: can't fillet this.
        m_pEdgeuses[0] = NULL;
        m_pEdgeuses[1] = NULL;
        SE(SM_ERR);
        return;
    }
    m_dOrientations[0] = dSides[0];
    m_dOrientations[1] = dSides[1];

    SmSurface *pSurface =
        m_pEdgeuses[0]->GetLoopuse()->GetFaceuse()->GetFace()->GetSurface();
    SmSurface *pSurface2 =
        m_pEdgeuses[1]->GetLoopuse()->GetFaceuse()->GetFace()->GetSurface();
    m_pSurfaces[0] = new (crContext) SmOffsetSurface(
        m_dFilletRadii[0]*m_dOrientations[0],*pSurface);
    m_pSurfaces[1] = new (crContext) SmOffsetSurface(
        m_dFilletRadii[1]*m_dOrientations[1],*pSurface2);
    m_pSurfaces[0]->SetAllowExtension( TRUE );
    m_pSurfaces[1]->SetAllowExtension( TRUE );

    (void)FinishConstruction();  // Take care of other things in the base class. // SolidCraft: a constructor cannot report

} // end SmConstantRadiusFS::SmConstantRadiusFS constructor


/*******************************************************************//**
PURPOSE: Destructor for Constant Radius Fillet Surface.

NOTES:
***********************************************************************/
SmConstantRadiusFS::~SmConstantRadiusFS
  ()
{
} // end SmConstantRadiusFS::~SmConstantRadiusFS destructor


/*******************************************************************//**
PURPOSE: Set the offset values to correspond to what they should be
     at this intersection point.

NOTES:
***********************************************************************/
SmStatus SmConstantRadiusFS::SetupOffsetValues
  (const SmTsectPnt &,      // in : not used in this function
   double * pOffsetDist)    // out: Optional fillet radius, NULL to ignore
{
  // get offset distance
  if (pOffsetDist) *pOffsetDist = m_dFilletRadii[0];

  // set offset surface offset distance - Orientation is +/- 1
  GetSurface(0)->SetOffsetDistance(m_dFilletRadii[0]*m_dOrientations[0]);
  GetSurface(1)->SetOffsetDistance(m_dFilletRadii[1]*m_dOrientations[1]);

  return SM_SUCCESS;

} // end SmConstantRadiusFS::SetupOffsetValues


/*******************************************************************//**
PURPOSE: SolidCraft change. The ball this solver rolls: one, of constant
    radius, when both faces are offset by the same radius.

NOTES: For SmFilletSolver::JudgeBallAt. The trace is free (its stations are
    planes square to the spine), so rbOnEdge is FALSE and asC unused.
***********************************************************************/
SmStatus SmConstantRadiusFS::GetBallLaw
  (const SmTsectPnt & /*crPnt*/,   // in : the station
   SmBoolean        & rbBall,      // out: TRUE = one ball rolls on both faces
   SmBoolean        & rbOnEdge,    // out: FALSE: a free trace
   double             adR[3],      // out: r, 0, 0
   SmVector3d         asC[4])      // out: unused
{
  SM_REF1(asC);
  const double dR = m_dFilletRadii[0];
  rbOnEdge = FALSE;
  rbBall   = ( dR > 0.0 && smos_Fabs( m_dFilletRadii[1] - dR ) <= SM_EFF_ZERO * ( 1.0 + dR )) ? TRUE : FALSE;
  adR[0] = dR;
  adR[1] = 0.0;
  adR[2] = 0.0;
  return SM_SUCCESS;
} // end SmConstantRadiusFS::GetBallLaw

/*******************************************************************//**
PURPOSE: Load the current Jacobian matrix in an incremental way
    and update the parameters.

NOTES:
  With: crX = [u1 v1 u2 v2]
         F  = OffsetSurf1.Position(u1, v1)
         G  = OffsetSurf2.Position(u2, v2)

  Sets:

  rF[]                 = done elsewhere
  rF[rlNumEquations+0] = (F-G).x
  rF[rlNumEquations+1] = (F-G).y
  rF[rlNumEquations+2] = (F-G).z

  pOptJacobian[][]                 = done elsewhere
  pOptJacobian[rlNumEquations+0][] = [ duf.X      dvf.X       -dug.X       -dvg.X ]
  pOptJacobian[rlNumEquations+1][] = [ duf.Y      duf.Y       -dug.Y       -dvg.Y ]
  pOptJacobian[rlNumEquations+2][] = [ duf.Z      duf.Z       -dug.Z       -dvg.Z ]

***********************************************************************/
SmStatus SmConstantRadiusFS::LoadJacobian
  (ULONG                  & rlNumEquations,       // i/o: number of equations to skip when setting
                                                  //      outputs rF and pOptJacobian
   ULONG                  & rlNumParameters,      // i/o: Gets incremented in an odd way to
                                                  //      to point to the input UVPnts.
   ULONG                    lRailIndex,           // in : Specifies rail index of first surface.
   SmSurface             *& rpSurface1,           // in : used only to map which crX[] uvPnt to surface1
   ULONG                  & rlSurf1Offset,        // i/o: index of rpSurface1 uvPnt in crX - ignored if rpSurface1 == NULL
   SmSurface             *& rpSurface2,           // in : used only to map which crX[] uvPnt to surface2
   ULONG                  & rlSurf2Offset,        // i/o: index of rpSurface2 uvPnt in crX - ignored if rpSurface2 == NULL
   const SmTArray<double> & crX,                  // in : current solution, [u1 v1 u2 v2 radius]
   SmTArray<double>       & rF,                   // out: sets rF[rlNumEquations] to rF[rlNumEquations+3] values
   SmMatrix               * pOptJacobian,         // out: sets pOptJacobian[rlNumEquations][] to pOptJacobian[rlNumEquations+3][] rows
   SmBoolean              & rbFoundAnswer)        // out: TRUE = all rF values less than scaled m_dConversionTol
{
  // init output
  rbFoundAnswer = FALSE;

  // Always need to set the offset surface offset distances prior to solving
  GetSurface(0)->SetOffsetDistance(m_dFilletRadii[0]*m_dOrientations[0]);
  GetSurface(1)->SetOffsetDistance(m_dFilletRadii[1]*m_dOrientations[1]);

  // get lSurf1Off index into crX array for lSurf1 parameters
  ULONG lSurf1Off = rlNumParameters;
  rlNumParameters = rlNumParameters + 2;
  if (   lRailIndex == 0
      && (   rpSurface1 != NULL
          && ((SmOffsetSurface*)rpSurface1)->GetBaseSurface() == GetSurface(0)->GetBaseSurface()))
    {
      lSurf1Off = rlSurf1Offset;
      rlNumParameters = rlNumParameters - 2;
    }
  else if (lRailIndex == 1 && (rpSurface2 != NULL &&
      ((SmOffsetSurface*)rpSurface2)->GetBaseSurface() ==
      GetSurface(0)->GetBaseSurface()))
    {
      lSurf1Off = rlSurf2Offset;
      rlNumParameters = rlNumParameters - 2;
    }

  // get offsetSurface1 F Point values for given crX values
  SmPoint3d sF;
  SmVector3d sDUF, sDVF;
  SmPoint2d sUV(crX[lSurf1Off],crX[lSurf1Off+1]);
  SER(GetSurface(0)->Evaluate1stDerivatives(sUV,TRUE,TRUE,sF,sDUF,sDVF));

#ifdef SM_DEBUG_CODE
SmBoolean bDebugMe = FALSE;
  if (bDebugMe) {
      sm_GraphicsLoop();
      smgfx_SetLook(1,2, 0,1,0); GetSurface(0)->DrawAt(sUV,1); sm_GraphicsLoop();
      sm_GraphicsLoop();
  }
#endif

  // get lSurf2Off index into crX array for lSurf2 parameters
  ULONG lSurf2Off = rlNumParameters;
  rlNumParameters = rlNumParameters + 2;
  if (lRailIndex == 1 && (rpSurface1 != NULL &&
      ((SmOffsetSurface*)rpSurface1)->GetBaseSurface() ==
      GetSurface(1)->GetBaseSurface())) {
      lSurf2Off = rlSurf1Offset;
      rlNumParameters = rlNumParameters - 2;
  }
  if (lRailIndex == 0 && (rpSurface2 != NULL &&
      ((SmOffsetSurface*)rpSurface2)->GetBaseSurface() ==
      GetSurface(1)->GetBaseSurface())) {
      lSurf2Off = rlSurf2Offset;
      rlNumParameters = rlNumParameters - 2;
  }

  // set rlSurfffset as needed
  if (lRailIndex == 0 && rpSurface1 == NULL) { rlSurf1Offset = lSurf1Off; }
  if (lRailIndex == 0 && rpSurface2 == NULL) { rlSurf2Offset = lSurf2Off; }
  if (lRailIndex == 1 && rpSurface1 == NULL) { rlSurf2Offset = lSurf1Off; }
  if (lRailIndex == 1 && rpSurface2 == NULL) { rlSurf1Offset = lSurf2Off; }

  // get offsetSurface2 G Point values for given crX values
  SmPoint3d sG;
  SmVector3d sDUG, sDVG;
  SmPoint2d sUV2(crX[lSurf2Off],crX[lSurf2Off+1]);
  SER(GetSurface(1)->Evaluate1stDerivatives(sUV2,TRUE,TRUE,sG,sDUG,sDVG));
#ifdef SM_DEBUG_CODE
  if (bDebugMe) {
      sm_GraphicsLoop();
      smgfx_SetColor(1,0,0);
      GetSurface(1)->DrawAt(sUV2,1);
      sm_GraphicsLoop();
  }
#endif

  // check surface tangent magnitudes
  double sDUFLenSq = sDUF.LengthSquared();
  double sDVFLenSq = sDVF.LengthSquared();

  double sDUGLenSq = sDUG.LengthSquared();
  double sDVGLenSq = sDVG.LengthSquared();

  // failure: zero length tangents
  if (   sDUFLenSq < SM_EFF_ZERO_SQ
      || sDVFLenSq < SM_EFF_ZERO_SQ
      || sDUGLenSq < SM_EFF_ZERO_SQ
      || sDVGLenSq < SM_EFF_ZERO_SQ)
    { return SM_ERR; }

  // set sDiff1 = distance between offset points
  SmVector3d sDiff1 = sF - sG;

#ifdef SM_DEBUG_CODE
  if (bDebugMe) {
      sm_GraphicsLoop();
      smgfx_SetLook(1,2, 1,0,0); sDiff1.Draw(&sG); sm_GraphicsLoop();
      sm_GraphicsLoop();
  }
#endif

  // Intersection tries to minimize the following six equations.
  // It is essentially intersecting the three surfaces.
  // (Three surfaces: each Fillet has two, one in common.)
  // fun[i+0] = (f-g) . X     - X component of vector f-g
  // fun[i+1] = (f-g) . Y     - Y component of vector f-g
  // fun[i+2] = (f-g) . Z     - Z component of vector f-g
  // Next 3 similar for (f-h).
  //
  // Where f is surface1, g is surface2 and h is surface3
  //
  // It produces the following Jacobian
  //  | duf.X      dvf.X       -dug.X       -dvg.X |
  //  | duf.Y      duf.Y       -dug.Y       -dvg.Y |
  //  | duf.Z      duf.Z       -dug.Z       -dvg.Z |

  if (pOptJacobian)
    {
      ULONG i = rlNumEquations;
      (*pOptJacobian)[i][lSurf1Off]   =   sDUF.x;
      (*pOptJacobian)[i][lSurf1Off+1] =   sDVF.x;
      (*pOptJacobian)[i][lSurf2Off]   = - sDUG.x;
      (*pOptJacobian)[i][lSurf2Off+1] = - sDVG.x;

      i++; // Move to next equation;
      (*pOptJacobian)[i][lSurf1Off]   =   sDUF.y;
      (*pOptJacobian)[i][lSurf1Off+1] =   sDVF.y;
      (*pOptJacobian)[i][lSurf2Off]   = - sDUG.y;
      (*pOptJacobian)[i][lSurf2Off+1] = - sDVG.y;

      i++;
      (*pOptJacobian)[i][lSurf1Off]   =   sDUF.z;
      (*pOptJacobian)[i][lSurf1Off+1] =   sDVF.z;
      (*pOptJacobian)[i][lSurf2Off]   = - sDUG.z;
      (*pOptJacobian)[i][lSurf2Off+1] = - sDVG.z;
    }

  // Compute function values
  rF[rlNumEquations++] = sDiff1.x;
  rF[rlNumEquations++] = sDiff1.y;
  rF[rlNumEquations++] = sDiff1.z;

  // See if we have converged
  double dScaledTol = GetConversionTol() * (1.0 + sF.GetMaxDimension());
  if (   smos_Fabs(sDiff1.x) < dScaledTol
      && smos_Fabs(sDiff1.y) < dScaledTol
      && smos_Fabs(sDiff1.z) < dScaledTol)
    {
      rbFoundAnswer = TRUE;
    }

#ifdef SM_DEBUG_CODE
  if (bDebugMe) {
      if (pOptJacobian) pOptJacobian->Dump();
      smos_WriteBuffer(_T(" F --- \n"));
      SM_DUMP_TARRAY( rF );
  }
#endif

  return SM_SUCCESS;

} // end SmConstantRadiusFS::LoadJacobian


/*******************************************************************//**
PURPOSE: Get the number of variables the Jacobian for the constant
     radius fillet contains.

NOTES:
***********************************************************************/
ULONG SmConstantRadiusFS::GetJacobianSize
  ()
 const
{
    return 4;

} // end SmConstantRadiusFS::GetJacobianSize

/*******************************************************************//**
PURPOSE: Load some initial values for the solver.

NOTES: Given the current input guesses on one surface,
  find a good guess on other surface for upcoming solve.
  Append the results to rIntervals (u- and v-domains of other surface),
  rPeriodicities (in u and v of other surface), and rGuessT (u and v
  guess values on the other surface, for the offset point determined
  by the input surface crSrf1Offset and the input uv values located
  at lSurface1Index in rGuessT).

***********************************************************************/
SmStatus SmConstantRadiusFS::LoadInitialValues
  (ULONG               & rbDoSurf2Calcs,  // in : if 0, don't have to calculate Surf values.
   ULONG                 lRailIndex,      // in : target rail index, 0 or 1
   const SmSurface     & crSrf1Offset,    // in : offset srf corresponding to lRailIndex
   ULONG                 lSurface1Index,  // in : index into rGuessT of srf 1's UV values
   SmSurface          *& rpSrf2Offset,    // out: pointer to other offset surface.
   ULONG               & rlSurface2Index, // out: index into rGuessT of srf2's UV values
   SmExtentNd          & rIntervals,      // i/o: u and v domains of srf2 appended
   SmTArray<SmBoolean> & rPeriodicities,  // i/o: u and v periodicities of srf2 appended.
   SmTArray<double>    & rGuessT)         // i/o: u and v guesses of srf2 appended
{
  // No work - second surface info already calculated.
  if( rbDoSurf2Calcs == 0 )
    { return SM_SUCCESS; }

  // get Other OffsetSurface
  ULONG lOtherRailIndex = 1 - lRailIndex;
  rpSrf2Offset          = this->GetSurface( lOtherRailIndex );

  // Locals
#ifdef SM_DEBUG_CODE
  SmEdgeuse       * pEdgeuse        = this->GetEdgeuse( lRailIndex );
#endif

  SmEdgeuse       * pOtherEdgeuse   = this->GetEdgeuse( lOtherRailIndex );
  const SmSurface * pSrf2Base       = ((SmOffsetSurface*)rpSrf2Offset)->GetBaseSurface();

  SmPoint3d         sSrf2GuessUV;   // This is what we'll be returning, in rGuessT.
  SmBoolean         bFoundSrf2UVGuess = FALSE;

  SmBoolean bTmp0 = FALSE ; 
  SmBSplineSurface * pBSplineSurface0 = SM_CAST_PTR(SmBSplineSurface, pSrf2Base) ; 
  SmTemporaryChangeValue<SmBoolean> sClean0(pBSplineSurface0 ? pBSplineSurface0->GetOutOfBoundsEnabled() : bTmp0, TRUE) ;

  SmSolution        sSData[16];
  SmSolutionArray   sSolutions(16,sSData);

  // Get given 3d guess point
  SmPoint3d sSrf1GuessPt3D;
  SmPoint2d sUV( rGuessT[lSurface1Index], rGuessT[lSurface1Index+1] );
  SER( crSrf1Offset.EvaluatePoint( sUV, sSrf1GuessPt3D ) );

#ifdef SM_DEBUG_CODE
SmBoolean bDebugMe = FALSE;
static ULONG lCount      = 1 ; lCount++ ;
static ULONG lDebugCount = 0 ;
  // Draw sSrf1GuessPt3D and see if it is on the offset Surface or not
  if (bDebugMe || lDebugCount == lCount)
    {
      SmFace    * pFace =  (SmFace *)crSrf1Offset.GetFace()  ? (SmFace *)crSrf1Offset.GetFace() 
                         : rpSrf2Offset                      ? (SmFace *)rpSrf2Offset->GetFace() 
                                                             : NULL ;
      SmEdge    * pEdge =  pEdgeuse      && pEdgeuse->GetEdge()      ? pEdgeuse->GetEdge() 
                         : pOtherEdgeuse && pOtherEdgeuse->GetEdge() ? pOtherEdgeuse->GetEdge() 
                                                                     : NULL ; 
      SmBrep    * pBrep =  pFace ? pFace->GetBrep() 
                         : pEdge ? pEdge->GetBrep() 
                                 : NULL ;

      smgfx_Erase();
      smgfx_SetLook(1,2, 0,0,1) ; if(pBrep) pBrep->Draw(TRUE) ; sm_GraphicsLoop() ;
      smgfx_SetLook(1,2, 0,1,1) ; crSrf1Offset.DrawUV() ; sm_GraphicsLoop() ;
      smgfx_SetLook(1,2, 1,1,0) ; rpSrf2Offset->DrawUV() ; sm_GraphicsLoop() ;
      smgfx_SetLook(1,2, 0,0,0) ; if(pFace) pFace->Draw(SM_DM_CROSSHATCH); sm_GraphicsLoop() ;
      smgfx_SetLook(3,4, 1,0,0) ; if(pEdgeuse) pEdgeuse->Draw(); sm_GraphicsLoop() ;
      smgfx_SetLook(3,4, 1,0,1) ; if(pOtherEdgeuse) pOtherEdgeuse->Draw(); sm_GraphicsLoop() ;
      smgfx_SetLook(5,7, 0,1,0) ; sSrf1GuessPt3D.Draw(); sm_GraphicsLoop() ;
      sm_GraphicsLoop() ;
    }
#endif // SM_DEBUG_CODE


  // In the tangent-rollover case, the proper uv guess points are set in the TsectPnts
  // of the FilletVertexuses.  They might be on the other side of a seam.

  SmFilletGeom *pFG = this->GetFilletGeomForSurfaces();
  if ( pFG != NULL  &&  pFG->GetFilletGeomType() == SM_FG_TANGENT_ROLLOVER )
  {

#ifdef SM_DEBUG_CODE
SmBoolean bDebugFVUs=FALSE;
      if ( bDebugFVUs )
        { pFG->DumpFilletVUs(); }
#endif

      SmFilletEdge *pRail = pFG->GetRail( lRailIndex );
      SmTArray< SmEdgeuse * > sAllEUs;
      SER( pRail->GetEdgeuses( sAllEUs ) );
      for ( ULONG ii=0; ii< sAllEUs.GetSize(); ii++ )
      {
          SmEdgeuse         *pEU  = sAllEUs[ ii ];
          SmVertexuse       *pVU  = pEU->GetVertexuse();
          SmFilletVertexuse *pFVU = SM_CAST_PTR( SmFilletVertexuse, pVU );
          if ( pFVU == NULL ) { continue; }
          SmVertex *pVtx = pFVU->GetVertex();
          SmFilletVertex *pFVtx = SM_CAST_PTR( SmFilletVertex, pVtx );
          if ( pFVtx == NULL ) { continue; }
          if ( pFVtx->GetFilletVertexType() != SM_FV_RAIL_X_EDGEUSE ) { continue; }

          // Here we have a FilletVertexuse from the stored FilletGeomForSurfaces.
          // Its uv should be set, and we should use it for srf2 guess.
          SmTsectPnt & rTsectPnt = pFVU->GetTsectPnt();
          SmIntersectionPointType eIntType = rTsectPnt.PointType();
          if ( eIntType == SM_IP_UNDEFINED || eIntType == SM_IP_UNKNOWN ) { continue; }
          sSrf2GuessUV = rTsectPnt.UVPos( 1-lRailIndex );

          bFoundSrf2UVGuess = TRUE;
          break;
      }

  } // end if our FilletGeomForSurfaces was set, and type TangentRollover.

  if ( ! bFoundSrf2UVGuess )
  {
      // get OtherSurface guess uv point

      // If we have edgeuses - use them.
      // If not - drop the GuessPoint to the other surface.

      // But first: if we are rolling on a different face (e.g., doing a rollover),
      // then this is the wrong Edgeuse, and we don't have an appropriate one to use.
      // [Fillet 408, 409]
      if ( pOtherEdgeuse != NULL )
      {
          SmFace *pF = pOtherEdgeuse->GetFace();
          if ( pF != NULL && pF->GetSurface() != pSrf2Base )
           { pOtherEdgeuse = NULL; }
      }

      if(   pOtherEdgeuse != NULL
         && pOtherEdgeuse->GetUVTrimCurve() )
      {
          // Drop pt to the filleted edge, in 3D
          SmCurve    *p3DFilletEdgeCurve = pOtherEdgeuse->GetEdge()->GetCurve();
          SmExtent1d  sIvl               = pOtherEdgeuse->GetEdge()->GetInterval();
          SER( p3DFilletEdgeCurve->GlobalPointSolve(sIvl, SM_SO_MINIMIZE,
                                                    sSrf1GuessPt3D,
                                                    m_dThisApproxTol3d, NULL, NULL,
                                                    SM_SR_SINGLE, sSolutions ));
          if ( sSolutions.GetSize() != 1 )
            { SER(SM_ERR); }

          // get OtherFilletEdge Parameter
          double dFilletEdgeT = sSolutions[0].m_vStart[0];

          // If the curve is closed, and the solution is within tol of an
          // end point, snap the solution to whichever end (start or finish)
          // such that the direction that moves towards the interior of the curve
          // also moves generally towards the interior of the surface.

          if ( p3DFilletEdgeCurve->IsClosed( sIvl ) )
          {
              SmPoint3d s3DPnt;
              SER( p3DFilletEdgeCurve->EvaluatePoint( dFilletEdgeT, s3DPnt ));
              SmVector3d sInwardVector;
              SER( crSrf1Offset.ComputeInwardVector( sUV, sInwardVector ));
              double dSnapTol = 0.001 * p3DFilletEdgeCurve->ApproximateLength( sIvl, 9 );
              double dDist, dParam;
              SmBoolean bSuccess;


              SER(p3DFilletEdgeCurve->DropPoint(sIvl,           // in : target curve allowed domain
                                                s3DPnt,         // in : Point to drop to curve
                                                &sInwardVector, // in : For Snapping, NULL=NoSnap, Vec pointing to CrvInside from desired EndPt for snapping.
                                                                //      when DropPt is within dDistTol in cpVecToCrvInside dir of EndPt, Snap to Ivl EndPt.
                                                                //      Vec handles ambiguities and makes sure that curves are not just touching at the ends.
                                                dSnapTol,       // in : Used for SnapDist in SM_SO_MINIZE case, which this is.
                                                NULL,           // in : If given, call cheap LocalPointSolve() else call expensive GlobalPointSolve()
                                                bSuccess,       // out: TRUE = found a drop point
                                                dParam,         // out: found drop curve param
                                                dDist)) ;       // out: found drop distance
                                                                // in : SM_SO_MINIMIZE = allow nonNormal drops near endPoints
                                                                //      SM_SO_INTERSECT= like Minimize but point must be within dDistTol
                                                                //      SM_SO_NORMALIZE= exclude nonNormal drops near endPoints
                                                                //      default:[SM_SO_MINIMIZE] to preserve original behavior
              if ( bSuccess )
                { dFilletEdgeT = dParam; }

          } // end if closed FilletEdgeCurve

          SmBSplineCurve *pOtherUVCurve = pOtherEdgeuse->GetUVTrimCurve();
          SER( pOtherUVCurve->EvaluatePoint( dFilletEdgeT, sSrf2GuessUV ));
          bFoundSrf2UVGuess = TRUE;

#ifdef SM_DEBUG_CODE
          // Draw sSrf1GuessPt3D and see if it is on the offset Surface or not
          if (bDebugMe)
          {
              SmPoint2d sUVDebug(sSrf2GuessUV.x, sSrf2GuessUV.y) ;
              SmFace    * pFace = (SmFace *)crSrf1Offset.GetFace() ;
              SmBrep    * pBrep = pFace ? pFace->GetBrep() : NULL ;
              SmPoint3d sEdgePoint, sBaseSurfPoint, sOffsetSurfPoint ;
              p3DFilletEdgeCurve->EvaluatePoint(dFilletEdgeT, sEdgePoint) ;
              pSrf2Base         ->EvaluatePoint( sUVDebug, sBaseSurfPoint) ;
              rpSrf2Offset      ->EvaluatePoint( sUVDebug, sOffsetSurfPoint) ;

              smgfx_Erase();
              smgfx_SetLook(1,2, 0,0,1) ; if(pBrep) pBrep->Draw(TRUE) ; sm_GraphicsLoop() ;
              smgfx_SetLook(1,2, 0,1,1) ; crSrf1Offset.DrawUV() ; sm_GraphicsLoop() ;
              smgfx_SetLook(1,2, 1,1,0) ; rpSrf2Offset->DrawUV() ; sm_GraphicsLoop() ;
              smgfx_SetLook(1,2, 0,0,0) ; if(pFace) pFace->Draw(SM_DM_CROSSHATCH); sm_GraphicsLoop() ;
              smgfx_SetLook(3,4, 1,0,0) ; if(pEdgeuse) pEdgeuse->Draw(); sm_GraphicsLoop() ;
              smgfx_SetLook(3,4, 1,0,1) ; if(pOtherEdgeuse) pOtherEdgeuse->Draw(); sm_GraphicsLoop() ;
              smgfx_SetLook(5,6, 0,1,0) ; sSrf1GuessPt3D.Draw(); sm_GraphicsLoop() ;
              smgfx_SetLook(6,7, 0,0,1) ; sEdgePoint.Draw(); sm_GraphicsLoop();
              smgfx_SetLook(7,8, 1,0,1) ; sBaseSurfPoint.Draw(); sm_GraphicsLoop();
              smgfx_SetLook(8,9, 1,0,0) ; sOffsetSurfPoint.Draw(); sm_GraphicsLoop();
              sm_GraphicsLoop() ;
      }
#endif // SM_DEBUG_CODE

      } // end given OtherFilletEdge check

      // If we found a uv, check that it's good enough.
      // Check how it looks on the other base surface.
      // If we are near a singularity then use global point solve.

      if ( bFoundSrf2UVGuess )
      {
          // Check parametric degeneracies:
          // (1) 1st derivs lining up,

          SmVector3d sDU, sDV;
          SmVector3d sMat[2][2];
          SmPoint2d sUVPnt( sSrf2GuessUV.x, sSrf2GuessUV.y );
          SER( pSrf2Base->Evaluate( sUVPnt, 1,1,TRUE,TRUE,TRUE, sMat[0] ));
          sDU = sMat[1][0];
          sDV = sMat[0][1];
          double dDeg;
          SER( sDU.AngleBetween( sDV, dDeg ));
          dDeg = SM_RAD2DEG( dDeg );
          if ( dDeg < 30 || dDeg > 150.0 )
          {
              bFoundSrf2UVGuess = FALSE;
          }

          // (2) one deriv going to zero.
          if (   sDU.GetMaxDimension() > 10000.0 * sDV.GetMaxDimension()
              || sDV.GetMaxDimension() > 10000.0 * sDU.GetMaxDimension() )
          {
              bFoundSrf2UVGuess = FALSE;
          }
      } // end found a UVGuessPoint from OtherEdgeuse check

      // If we didn't get a good enough uv on surface2, try a global solve
      if ( !bFoundSrf2UVGuess )
      {
          // Note: should we extend the domain (SetupOffsetExtension())?
          SER( pSrf2Base->GlobalPointSolve(pSrf2Base->GetNaturalUVDomain(),
                                           SM_SO_MINIMIZE, sSrf1GuessPt3D,
                                           SM_EFF_ZERO_SQ, NULL,
                                           SM_SR_SINGLE, sSolutions ));

          if ( sSolutions.GetSize() != 1 )
                { SER(SM_ERR); }

          bFoundSrf2UVGuess = TRUE;
          sSrf2GuessUV.x = sSolutions[0].m_vStart[0];
          sSrf2GuessUV.y = sSolutions[0].m_vStart[1];

      } // end need to use globalPointSolve to find OtherSurface UVGuessPoint check

      // Before we add the rGuessT values, drop point to other offset surface
      // to get a better guess
      SmBoolean  bFound;
      SmSolution sSol;
      SmPoint2d  sUVGuess( sSrf2GuessUV.x, sSrf2GuessUV.y );

      SmExtent2d sDomain = pSrf2Base->GetNaturalUVDomain();
      SetupOffsetExtension( pSrf2Base, sDomain ); // expands sDomain

      SER( rpSrf2Offset->LocalPointSolve(sDomain, SM_SO_MINIMIZE,
                                     sSrf1GuessPt3D, sUVGuess,
                                     bFound, sSol ));
      // when local solve worked
      if ( bFound )
      { // save solution
          double dDist = sSol.m_vStart.m_dSolutionValue;
          if ( dDist < (  smos_Fabs( m_dFilletRadii[0] )
                        + smos_Fabs( m_dFilletRadii[1] ) ) )
          {
              sSrf2GuessUV.x = sSol.m_vStart[0];
              sSrf2GuessUV.y = sSol.m_vStart[1];
          }
      }
      else // use global solve
      {
          // Note: should we extend the domain (use sDomain)?
          SER( pSrf2Base->GlobalPointSolve(pSrf2Base->GetNaturalUVDomain(),
                                           SM_SO_MINIMIZE, sSrf1GuessPt3D,
                                           SM_EFF_ZERO_SQ, NULL,
                                           SM_SR_SINGLE, sSolutions));

          if ( sSolutions.GetSize() == 1 )
          {
              sSrf2GuessUV.x = sSolutions[0].m_vStart[0];
              sSrf2GuessUV.y = sSolutions[0].m_vStart[1];
          }
      } // end drop GuessPoint to offsetSurface branches

  } // end if not bFoundSrf2UVGuess


  // Start setting the return values:
  // the next two entries in rIntervals and rPeriodicities.

  // Note: we're expanding rIntervals, but we also want to expand
  // sDomain, to allow the surface drop to converge slightly outside
  // the domain if that's where it is.
  // [BD; 25 Oct 05; 050718]

  SmExtent2d sDomain = pSrf2Base->GetNaturalUVDomain();
  SetupOffsetExtension( pSrf2Base, sDomain ); // expands sDomain

  rlSurface2Index = rGuessT.GetSize();
  rIntervals[ rlSurface2Index   ] = sDomain.GetUInterval();
  rIntervals[ rlSurface2Index+1 ] = sDomain.GetVInterval();

  rPeriodicities.Add(FALSE);
  rPeriodicities.Add(FALSE);

  // set output
  rGuessT.Add( sSrf2GuessUV.x );
  rGuessT.Add( sSrf2GuessUV.y );

#ifdef SM_DEBUG_CODE
  // Draw sSrf1GuessPt3D and see if it is on the offset Surface or not
  if (bDebugMe)
  {
      SmFace    * pFace = (SmFace *)crSrf1Offset.GetFace() ;
      SmBrep    * pBrep = pFace ? pFace->GetBrep() : NULL ;
      SmPoint2d sUVPnt(sSrf2GuessUV.x,sSrf2GuessUV.y);
      SmPoint3d s3DPnt;
      rpSrf2Offset->EvaluatePoint(sUVPnt,s3DPnt);

      smgfx_Erase();
      smgfx_SetLook(1,2, 0,0,1) ; if(pBrep) pBrep->Draw(TRUE) ; sm_GraphicsLoop() ;
      smgfx_SetLook(1,2, 0,1,1) ; crSrf1Offset.DrawUV() ; sm_GraphicsLoop() ;
      smgfx_SetLook(1,2, 1,1,0) ; rpSrf2Offset->DrawUV() ; sm_GraphicsLoop() ;
      smgfx_SetLook(1,2, 0,0,0) ; if(pFace) pFace->Draw(SM_DM_CROSSHATCH); sm_GraphicsLoop() ;
      smgfx_SetLook(3,4, 1,0,0) ; if(pEdgeuse) pEdgeuse->Draw(); sm_GraphicsLoop() ;
      smgfx_SetLook(3,4, 1,0,1) ; if(pOtherEdgeuse) pOtherEdgeuse->Draw(); sm_GraphicsLoop() ;
      smgfx_SetLook(5,7, 0,1,0) ; sSrf1GuessPt3D.Draw(); sm_GraphicsLoop() ;
      smgfx_SetLook(6,8, 1,0,0) ; s3DPnt.Draw(); sm_GraphicsLoop();
      sm_GraphicsLoop() ;
  }
#endif

  // all done
  return SM_SUCCESS;

} // end SmConstantRadiusFS::LoadInitialValues

/*******************************************************************//**
PURPOSE: Utility methods.
NOTES:
***********************************************************************/
SmBoolean SmConstantRadiusFS::IsKindOf( SM_TYPE t ) const
{
  return ((SmConstantRadiusFS_TYPE == t) ? TRUE : SmFilletSolver::IsKindOf( (t) ));
}
/*******************************************************************//**
PURPOSE: Dump a SmConstantRadiusFS.

NOTES: Unimplemented.
***********************************************************************/
void SmConstantRadiusFS::Dump() const {}


/*******************************************************************//**
PURPOSE: Constructor for SmConstantRadiusAssistedFS object where
    the fillet will require the help from SmConstantRadiusFS object

NOTES: The input pEdgeuse corresponds to the surface whose offset
    distance is dFilletRadius1.
***********************************************************************/
SmConstantRadiusAssistedFS::SmConstantRadiusAssistedFS
  (const SmContext & crContext,
   double dThisApproxTol3d,
   double dAngleTolerance,
   double dTangencyTolerance,
   double dFilletRadius1,
   double dFilletRadius2,
   SmEdgeuse *pEdgeuse)
  : SmConstantRadiusFS(crContext,
                       dThisApproxTol3d,
                       dAngleTolerance,
                       dTangencyTolerance,
                       dFilletRadius1,
                       pEdgeuse,
                       &dFilletRadius2)
{

} // end SmConstantRadiusAssistedFS::SmConstantRadiusAssistedFS


/*******************************************************************//**
PURPOSE: Compute the fillet values at given distance from existing point
***********************************************************************/
class SmGuessSphereSolveENFO : public SmFilletSphereSolveENFO
{
public:
    SmGuessSphereSolveENFO(const SmPoint3d & crSphereCenter,
                           double dSphereRadius,
                           SmFilletSolver & crFilletSolver)
      :SmFilletSphereSolveENFO(crSphereCenter,dSphereRadius,crFilletSolver) {}
    virtual SmStatus Evaluate(const SmTArray<double> & crX,             // in : x of Ax=F
                              SmTArray<double>       & rF,              // out: F of Ax=F function values,
                              SmMatrix               * pOptJacobian,    // out: Partial derivatives of the functions.
                              SmBoolean              & rbFoundAnswer);  // out: Not always used, when used
                                                                        //      TRUE = converged (F members are within tolerance of 0.0
                                                                        //      FALSE= Not Used or Not Converged
} ; // end class SmGuessSphereSolveENFO


/*******************************************************************//**
PURPOSE:

NOTES:
***********************************************************************/
SmStatus SmGuessSphereSolveENFO::Evaluate
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

    ULONG lNumEquations = 1;
    ULONG lNumParameters = 2;
    SmSurface *pSurface = m_crFilletSolver.GetSurface(0);
    ULONG lSurf1Offset = 0;
    SmSurface * pSurface2 = NULL;
    ULONG lSurf2Offset = 0;
    SmBoolean bFoundAnswer;
    ULONG lRailIndex = 0;

    SmConstantRadiusFS * cpFS = (SmConstantRadiusFS*)&m_crFilletSolver;
    if (cpFS->SmConstantRadiusFS::LoadJacobian(lNumEquations, lNumParameters,
        lRailIndex, pSurface, lSurf1Offset, pSurface2, lSurf2Offset,
        crX, rF, pOptJacobian, bFoundAnswer) != SM_SUCCESS) {
        return SM_ERR;
    }

    SmPoint2d sUV(crX[0],crX[1]);
    SmPoint3d sF;
    SmVector3d sDUF, sDVF;
    SER(m_crFilletSolver.GetSurface(0)->Evaluate1stDerivatives(sUV,TRUE,TRUE,sF,sDUF,sDVF));

#ifdef SM_DEBUG_CODE
SmBoolean bDebugMe2 = FALSE;
    if (bDebugMe2) {
        sm_GraphicsLoop();
        smgfx_SetColor(1,0,0);
        m_crFilletSolver.GetSurface(0)->DrawAt(sUV,1);
        smgfx_SetColor(0,0,1);
        SmPoint2d sUV2(crX[lSurf2Offset],crX[lSurf2Offset+1]);
        m_crFilletSolver.GetSurface(1)->DrawAt(sUV2,1);
        sm_GraphicsLoop();
    }
#endif
    // Let's try using the plane instead of the sphere equation
    // f0 = A*x + B*y + C*z + D
    // || A*du.x + B*du.y + C*du.z    A*dv.x + B*dv.y + C*dv.z

    if (m_bPlaneAdjusts) {
        SmVector3d sNorm1 = sDUF * sDVF;
        if (sNorm1.LengthSquared() < SM_EFF_ZERO_SQ) {
            SER(m_crFilletSolver.GetSurface(0)->EvaluateNormal(sUV,TRUE,TRUE,sNorm1));
        }
        SER(sNorm1.Unitize());

        SmPoint2d sUV2(crX[lSurf2Offset],crX[lSurf2Offset+1]);
        SmPoint3d sG;
        SmVector3d sDUG, sDVG;
        SER(m_crFilletSolver.GetSurface(1)->Evaluate1stDerivatives(sUV2,TRUE,TRUE,sG,sDUG,sDVG));

        SmVector3d sNorm2 = sDUG * sDVG;
        if (sNorm2.LengthSquared() < SM_EFF_ZERO_SQ) {
            SER(m_crFilletSolver.GetSurface(1)->EvaluateNormal(sUV2,TRUE,TRUE,sNorm2));
        }
        SER(sNorm2.Unitize());

        SmPoint3d sPlaneOrig = (sF + sG) / 2.0;
        SmVector3d sPlaneABC = sNorm1 * sNorm2;
        if (sPlaneABC.Dot(m_vPlaneABC) < 0.0) { sPlaneABC = - sPlaneABC; }
        m_vPlaneABC = sPlaneABC;
        m_dPlaneD = - m_vPlaneABC.Dot(sPlaneOrig);
     }

    rF[0] = m_vPlaneABC.Dot(sF) + m_dPlaneD;

    if (pOptJacobian) {
        (*pOptJacobian)[0][0] = m_vPlaneABC.Dot(sDUF);
        (*pOptJacobian)[0][1] = m_vPlaneABC.Dot(sDVF);
    }

    // Jacobian answer looks good now check our answer to see if we return
    // TRUE for rbFoundAnswer
    if (bFoundAnswer) {
        double dScaledTol = m_crFilletSolver.GetConversionTol() * (1.0 + sF.GetMaxDimension());
        if (smos_Fabs(rF[0]) < dScaledTol) {
            rbFoundAnswer = TRUE;
        }
    }

#ifdef SM_DEBUG_CODE
SmBoolean bDebugMe = FALSE;
    if (bDebugMe) {
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
            SM_DUMP_TARRAY( rF );
        }
    }
#endif

    return SM_SUCCESS;
} // end SmGuessSphereSolveENFO::Evaluate

/*******************************************************************//**
PURPOSE: Solve fillet point on a given plane

NOTES:
***********************************************************************/
SmStatus SmConstantRadiusAssistedFS::GuessPointOnPlaneSolve
  (const SmPoint3d & crPlaneOrig,
   const SmVector3d & crPlaneNormal,
   const SmExtent2d & crUVDomain1,
   const SmExtent2d & crUVDomain2,
   const SmVector2d & rUV1,
   const SmVector2d & rUV2,
   SmBoolean & rbFoundSolution,
   SmTsectPnt & rTsectPnt)
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

    //ULONG lNumVariables = 4; // Here is how many we have loaded so far
    ULONG lSurf1Index = 0; // UV values of Surface1
    ULONG lSurf2Index = 2; // UV values of Surface2
    //ULONG lRailIndex = 0;

#ifdef SM_DEBUG_CODE
SmBoolean bDebugMe = FALSE;
    if (bDebugMe) {
        smgfx_Erase();
        smgfx_SetColor(1,0,0);
        crPlaneNormal.Draw(&crPlaneOrig);
        smgfx_SetColor(0,0,1);
        GetSurface(0)->DrawUV(2,2);
        GetSurface(1)->DrawUV(2,2);
        sm_GraphicsLoop();
    }
#endif

    double sdData[16];
    SmTArray<double> sSolutionVector(16,sdData);

    double dPlaneD = - crPlaneOrig.Dot(crPlaneNormal);
    SmGuessSphereSolveENFO sEvalFun(crPlaneNormal,dPlaneD,*this);

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

    // Solve it
    SER(sLS.SolveIt(sGuessT,m_dThisApproxTol3d,rbFoundSolution,sSolutionVector));

    // Try it again and allow the plane to adjust as convergance happens
    if (!rbFoundSolution) {
        sEvalFun.SetPlaneAdjusts(TRUE);
        SER(sLS.SolveIt(sGuessT,m_dThisApproxTol3d,rbFoundSolution,sSolutionVector));
    }

    if (rbFoundSolution) {
        rTsectPnt.UVPos(0).x = sSolutionVector[lSurf1Index];
        rTsectPnt.UVPos(0).y = sSolutionVector[lSurf1Index+1];
        rTsectPnt.UVPos(1).x = sSolutionVector[lSurf2Index];
        rTsectPnt.UVPos(1).y = sSolutionVector[lSurf2Index+1];
        for (ULONG i=4; i<sSolutionVector.GetSize(); i++) {
            rTsectPnt.m_adUserDoubles[i-4] = sSolutionVector[i];
        }
    }

    return SM_SUCCESS;

} // end SmConstantRadiusAssistedFS::GuessPointOnPlaneSolve


/*******************************************************************//**
PURPOSE: Find fillet point on a given plane.

NOTES: Finds and stores surface uvPnt values such that
       BaseSurface(uvPnt)   = point on rails
       OffsetSurface(uvPnt) = point on fillet center curve (offsetSurface intersection curve)
                              and on specified plane

       Also finds and stores Surface position, tangent, and normal values
       at the solution as well as the center-curve position.

METHOD ---
  Given a point on and a plane perpendicular to the
     fillet offset-surface intersection-curve,
  Find the two intersection points between that plane and the
     fillet rail curves.

  The intersection points are found satisfying the particular equation
  set implemented within the virtual SmFilletSolver::LoadJacobian() methods
  for each SmFilletSolver derived type.
***********************************************************************/
SmStatus SmConstantRadiusAssistedFS::PointOnPlaneSolve
  (const SmPoint3d  & crPlaneOrig,         // in : pt on offset-surface xsect curve
   const SmVector3d & crPlaneNormal,       // in : tang at that offset-surface xsect curve
   const SmExtent2d & crUVDomain1,         // in : fillet-surf1 UVdomain
   const SmExtent2d & crUVDomain2,         // in : fillet-surf2 UVdomain
   const SmVector2d & rUV1,                // in : fillet-surf1 UV guess value
   const SmVector2d & rUV2,                // NotUsed: in : fillet-surf2 UV guess value
   SmBoolean        & rbFoundSolution,     // out: TRUE = next rTsectPnt found with
                                           //             SmLocalSolveNd:SmFilletSphereSolveENFO::SolveIt().
                                           //      FALSE= no next point found.
   SmTsectPnt       & rTsectPnt)           // out: gets UV values for both surfaces
                                           //        BaseSurface(uvPnt)   = points on rails
                                           //        OffsetSurface(uvPnt) = point on fillet center curve
{
  SM_REF1(rUV2) ;
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

  double adGuessData[16];
  SmTArray<double> sGuessT(16,adGuessData);
  sGuessT.Add(rUV1.x);
  sGuessT.Add(rUV1.y);

  SmSurface *pSurface1 = GetSurface(0);
  SmSurface *pSurface2 = NULL;
  ULONG lDoSurf2Calcs  = 1; // Do surface2 calculations.
  ULONG lSurf1Index    = 0;   // UV values of Surface1
  ULONG lSurf2Index;
  ULONG lRailIndex     = 0;

#ifdef SM_DEBUG_CODE
SmBoolean bDebugMe = FALSE;
  if (bDebugMe) {
      smgfx_Erase();
      smgfx_SetColor(1,0,0);
      crPlaneNormal.Draw(&crPlaneOrig);
      smgfx_SetColor(0,0,1);
      GetSurface(0)->DrawUV(2,2);
      GetSurface(1)->DrawUV(2,2);
      sm_GraphicsLoop();
  }
#endif

  // Load the initial values
  // given current inputs - place a good guess for upcoming solve in sGuessT
  SER(LoadInitialValues(lDoSurf2Calcs,lRailIndex,
      *pSurface1, lSurf1Index, pSurface2, lSurf2Index,
      sIntervals, sPeriodicities, sGuessT));

  double sdData[16];
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
      rTsectPnt.UVPos(0).x = sSolutionVector[lSurf1Index];
      rTsectPnt.UVPos(0).y = sSolutionVector[lSurf1Index+1];
      rTsectPnt.UVPos(1).x = sSolutionVector[lSurf2Index];
      rTsectPnt.UVPos(1).y = sSolutionVector[lSurf2Index+1];
      for (ULONG i=4; i<sSolutionVector.GetSize(); i++) {
          rTsectPnt.m_adUserDoubles[i-4] = sSolutionVector[i];
      }

      // for both surfaces - get and store solution point surface values
      for (ULONG lSrf=0; lSrf<=1; lSrf++)
        {
          // get and store the surface point values
          SmPoint2d sUV = rTsectPnt.UVPos(lSrf);
          SER(GetSurface(lSrf)->Evaluate1stDerivatives(sUV,TRUE,TRUE,
              rTsectPnt.SrfPos(lSrf),
              rTsectPnt.SrfDu(lSrf),
              rTsectPnt.SrfDv(lSrf)));

          // Store a unitized surface normal.
          SmVector3d sNorm =   rTsectPnt.SrfDu(lSrf)
                             * rTsectPnt.SrfDv(lSrf);
          if (sNorm.LengthSquared() < SM_EFF_ZERO_SQ) {
              SER(GetSurface(lSrf)->EvaluateNormal(sUV,TRUE,TRUE,
                  sNorm));
              if (sNorm.LengthSquared() < SM_EFF_ZERO_SQ) {
                  SER(SM_ERR);
              }
          }
          (void)sNorm.Unitize(); // SolidCraft: a zero vector stays zero
          rTsectPnt.SrfNorm(lSrf) = sNorm;

        } // end iter both surfaces

      // Now compute the point and vectors for the 3D curve
      // Point is average of surface points
      rTsectPnt.CrvPos() = rTsectPnt.SrfPos(0);
    }

  return SM_SUCCESS;

} // end SmConstantRadiusAssistedFS::PointOnPlaneSolve

/*******************************************************************//**
PURPOSE: Utility methods.

NOTES:
***********************************************************************/
SmBoolean SmConstantRadiusAssistedFS::IsKindOf( SM_TYPE t ) const
{
  return ((SmConstantRadiusAssistedFS_TYPE == t) ? TRUE : SmConstantRadiusFS::IsKindOf( (t) ));
}
/*******************************************************************//**
PURPOSE: Dump a SmConstantRadiusAssistedFS.

NOTES: Unimplemented.
***********************************************************************/
void SmConstantRadiusAssistedFS::Dump() const {}


/*******************************************************************//**
PURPOSE: Constructor for Constant Distance Fillet object.

NOTES:
***********************************************************************/
SmConstantDistanceFS::SmConstantDistanceFS
  (const SmContext & crContext,
   double            dThisApproxTol3d,
   double            dAngleTolerance,
   double            dTangencyTolerance,
   double            dDistance, // Euclidian distance between the two rails
   SmEdgeuse       * pEdgeuse)
 : SmFilletSolver(crContext,
                  dThisApproxTol3d,
                  dAngleTolerance,
                  dTangencyTolerance),
   m_dDistance(dDistance)
{
    SM_ASSERT(pEdgeuse != NULL);
    m_pEdgeuses[0] = pEdgeuse;
    m_pEdgeuses[1] = pEdgeuse->GetRadial();

    // Create offset surfaces of appropriate radius in appropriate
    // direction.  Note that we may need to adjust this if the two
    // surfaces go to tangency and then back down in the opposite direction.
    // For now assume that the orientations at the center of the edge
    // work for the rest of the edge.

    // SolidCraft change: offset each surface to the side of it the other face lies on
    // (sm_SectorSides), a cusp's faces too.
    double dSides[2] = { 1.0, 1.0 };
    if (sm_SectorSides(m_pEdgeuses, m_dDistance, dSides) != SM_SUCCESS)
    {
        // Smooth, or no sector: can't fillet this.
        m_pEdgeuses[0] = NULL;
        m_pEdgeuses[1] = NULL;
        SE(SM_ERR);
        return;
    }
    m_dOrientations[0] = dSides[0];
    m_dOrientations[1] = dSides[1];

    SmSurface *pSurface = m_pEdgeuses[0]->GetLoopuse()->GetFaceuse()->GetFace()->GetSurface();
    SmSurface *pSurface2 = m_pEdgeuses[1]->GetLoopuse()->GetFaceuse()->GetFace()->GetSurface();
    double dApproxRadius = m_dDistance * smos_Sqrt(2.0) / 2.0;
    m_pSurfaces[0] = new (crContext) SmOffsetSurface(dApproxRadius*m_dOrientations[0],*pSurface);
    m_pSurfaces[1] = new (crContext) SmOffsetSurface(dApproxRadius*m_dOrientations[1],*pSurface2);
    m_pSurfaces[0]->SetAllowExtension(TRUE);
    m_pSurfaces[1]->SetAllowExtension(TRUE);

    (void)FinishConstruction();  // Take care of other things in the base class. // SolidCraft: a constructor cannot report

#ifdef SM_DEBUG_CODE
SmBoolean bDebugMe = FALSE;
    if (bDebugMe) {
        sm_GraphicsLoop();
        smgfx_SetColor(0,0,0);
        pSurface->DrawUV(4,8);
        sm_GraphicsLoop();
        smgfx_SetColor(0,0,1);
        pSurface2->DrawUV(9,18);
        sm_GraphicsLoop();
    }
#endif


} // end SmConstantDistanceFS::SmConstantDistanceFS


/*******************************************************************//**
PURPOSE: Destructor for Constant Radius Distance fillet solver.

NOTES:
***********************************************************************/
SmConstantDistanceFS::~SmConstantDistanceFS
  ()
{

} // end SmConstantDistanceFS::~SmConstantDistanceFS destructor


/*******************************************************************//**
PURPOSE: Set the offset values to correspond to what they should be
     at this intersection point.

NOTES:
***********************************************************************/
SmStatus SmConstantDistanceFS::SetupOffsetValues
  (const SmTsectPnt & rTsectPnt,    // in : each TsectPnt stores its own radius to preserve constant distance Fillets
   double * pOffsetDist)            // out: optional fillet radius at this point
{
  // get fillet radius at this point
  double dRadius = rTsectPnt.m_adUserDoubles[0];

  // save fillet radius and set OffsetSurface Offset distances
  if (pOffsetDist) *pOffsetDist = dRadius;
  GetSurface(0)->SetOffsetDistance(dRadius*m_dOrientations[0]);
  GetSurface(1)->SetOffsetDistance(dRadius*m_dOrientations[1]);

  return SM_SUCCESS;

} // end SmConstantDistanceFS::SetupOffsetValues

/*******************************************************************//**
PURPOSE: Load the current Jacobian matrix in an incremental way
    and update the parameters.

NOTES: The Jacobian of this one uses a ball with a variable
    radius and the distance between the points equal to the length.

  With: crX = [u1 v1 u2 v2 Radius]
         F  = BaseSurf1.Position(u1, v1)
         G  = BaseSurf2.Position(u2, v2)
        nF  = BaseSurf1.UnitSurfaceNormal(u1, v1)
        nG  = BaseSurf2.UnitSurfaceNormal(u2, v2)
         R  = Radius
         L  = m_dDistance, specified distance between rail curves
  Sets:

  rF[]                 = done elsewhere
  rF[rlNumEquations+0] = (F-G+R*nF-R*nG).x
  rF[rlNumEquations+1] = (F-G+R*nF-R*nG).y
  rF[rlNumEquations+2] = (F-G+R*nF-R*nG).z
  rF[rlNumEquations+3] = (F-G)*(F-G) - L*L

  pOptJacobian[][]                 = done elsewhere
  pOptJacobian[rlNumEquations+0][] = [ duf.X+R(d|Nf|/du).X   dvf.X+R(d|Nf|/dv).X   -dug.X-R(d|Ng|/du).X   -dvg.X-R(d|Ng|/dv).X   Nf.x-Ng.x ]
  pOptJacobian[rlNumEquations+1][] = [ duf.Y+R(d|Nf|/du).Y   dvf.Y+R(d|Nf|/dv).Y   -dug.Y-R(d|Ng|/du).Y   -dvg.Y-R(d|Ng|/dv).Y   Nf.y-Ng.y ]
  pOptJacobian[rlNumEquations+2][] = [ duf.Z+R(d|Nf|/du).Z   dvf.Y+R(d|Nf|/dv).Z   -dug.Y-R(d|Ng|/du).Z   -dvg.Y-R(d|Ng|/dv).Z   Nf.z-Ng.z ]
  pOptJacobian[rlNumEquations+3][] = [ 2*(f-g)*duf           2*(f-g)*dvf           -2*(f-g)*dug           -2*(f-g)*dvg           0.0       ]

***********************************************************************/
SmStatus SmConstantDistanceFS::LoadJacobian
  (ULONG                  & rlNumEquations,       // i/o: number of equations to skip when setting
                                                  //      outputs rF and pOptJacobian
   ULONG                  & rlNumParameters,      // i/o: Gets incremented in an odd way to
                                                  //      to point to input UVPnts and Radius at appropriate times.
   ULONG                    lRailIndex,           // in : Specifies rail index of first surface.
   SmSurface             *& rpSurface1,           // in : used only to map which crX[] uvPnt to surface1
   ULONG                  & rlSurf1Offset,        // in : index of rpSurface1 uvPnt in crX - ignored if rpSurface1 == NULL
   SmSurface             *& rpSurface2,           // in : used only to map which crX[] uvPnt to surface2
   ULONG                  & rlSurf2Offset,        // in : index of rpSurface2 uvPnt in crX - ignored if rpSurface2 == NULL
   const SmTArray<double> & crX,                  // in : current solution, [u1 v1 u2 v2 radius]
   SmTArray<double>       & rF,                   // out: sets rF[rlNumEquations] to rF[rlNumEquations+3] values
   SmMatrix               * pOptJacobian,         // out: sets pOptJacobian[rlNumEquations][] to pOptJacobian[rlNumEquations+3][] rows
   SmBoolean              & rbFoundAnswer)        // out: TRUE = all rF values less than scaled m_dConversionTol
{
  // init output
  rbFoundAnswer = FALSE;

  // get lSurf1Off index into crX array for lSurf1 parameters
  ULONG lSurf1Off = rlNumParameters;
  if (   lRailIndex == 0
      && (   rpSurface1 != NULL
          && ((SmOffsetSurface*)rpSurface1)->GetBaseSurface() == GetSurface(0)->GetBaseSurface()))
    {
      lSurf1Off = rlSurf1Offset;
    }
  else if (   lRailIndex == 1
           && (   rpSurface2 != NULL
               && ((SmOffsetSurface*)rpSurface2)->GetBaseSurface() == GetSurface(0)->GetBaseSurface()))
    {
      lSurf1Off = rlSurf2Offset;
    }
  else
    {
      rlNumParameters = rlNumParameters + 2;
    }

  // get lSurf2Off index into crX array for lSurf2 parameters
  ULONG lSurf2Off = rlNumParameters;
  if (   lRailIndex == 1
      && (   rpSurface1 != NULL
          && ((SmOffsetSurface*)rpSurface1)->GetBaseSurface() == GetSurface(1)->GetBaseSurface()))
    {
      lSurf2Off = rlSurf1Offset;
    }
  else if (   lRailIndex == 0
           && (   rpSurface2 != NULL
               && ((SmOffsetSurface*)rpSurface2)->GetBaseSurface() == GetSurface(1)->GetBaseSurface()))
    {
      lSurf2Off = rlSurf2Offset;
    }
  else
    {
      rlNumParameters = rlNumParameters + 2;
    }

  // get baseSurface1 F point values for given crX values
  SmPoint3d  sF ;
  SmVector3d sDUF, sDVF, sDUVF, sDUUF, sDVVF ;
  SmPoint2d  sUV(crX[lSurf1Off],crX[lSurf1Off+1]) ;
  GetSurface(0)->SetOffsetDistance(0.0) ;
  SER(GetSurface(0)->Evaluate2ndDerivatives(sUV,TRUE,TRUE,sF,
                                             sDUF,sDVF,sDUVF,sDUUF,sDVVF));
#ifdef SM_DEBUG_CODE
SmBoolean bDebugMe = FALSE;
  if (bDebugMe) {
      sm_GraphicsLoop();
      smgfx_SetColor(0,1,0);
      GetSurface(0)->DrawAt(sUV,1);
      sm_GraphicsLoop();
  }
#endif

  // get baseSurface2 G point values for given crX values
  SmPoint3d sG;
  SmVector3d sDUG,sDVG,sDUVG,sDUUG,sDVVG;
  SmPoint2d sUV2(crX[lSurf2Off],crX[lSurf2Off+1]);
  GetSurface(1)->SetOffsetDistance(0.0);
  SER(GetSurface(1)->Evaluate2ndDerivatives(sUV2,TRUE,TRUE,sG,
                                              sDUG,sDVG,sDUVG,sDUUG,sDVVG));
#ifdef SM_DEBUG_CODE
  if (bDebugMe) {
      sm_GraphicsLoop();
      smgfx_SetColor(1,0,0);
      GetSurface(1)->DrawAt(sUV2,1);
      sm_GraphicsLoop();
  }
#endif

  // get current radius value
  double dRadius = crX[rlNumParameters];
  rlNumParameters ++;

  // check surface tangent magnitudes
  double sDUFLenSq = sDUF.LengthSquared();
  double sDVFLenSq = sDVF.LengthSquared();

  double sDUGLenSq = sDUG.LengthSquared();
  double sDVGLenSq = sDVG.LengthSquared();

  // failure: zero length tangents
  if (   sDUFLenSq < SM_EFF_ZERO_SQ
      || sDVFLenSq < SM_EFF_ZERO_SQ
      || sDUGLenSq < SM_EFF_ZERO_SQ
      || sDVGLenSq < SM_EFF_ZERO_SQ)
    { return SM_ERR; }

  // get baseSurface SurfaceNormals at current solution points
  SmVector3d sNF = sDUF * sDVF; // Normal on base surface #1
  SmVector3d sNG = sDUG * sDVG; // Normal on base surface #2

  // failure: zero length SurfaceNormal
  if (   sNF.LengthSquared() < SM_EFF_ZERO_SQ
      || sNG.LengthSquared() < SM_EFF_ZERO_SQ)
    { SER(SM_ERR); }

  // get derivatives of SurfaceNormals
  SmVector3d sNFdu = sDUF.DerivativeOfCrossProduct(sDVF,sDUUF,sDUVF);
  SmVector3d sNFdv = sDUF.DerivativeOfCrossProduct(sDVF,sDUVF,sDVVF);
  SmVector3d sNGdu = sDUG.DerivativeOfCrossProduct(sDVG,sDUUG,sDUVG);
  SmVector3d sNGdv = sDUG.DerivativeOfCrossProduct(sDVG,sDUVG,sDVVG);

  // normalize SurfaceNormal derivatives
  SmVector3d sUnitNFdu = m_dOrientations[0]*sNF.UnitizedDerivative(sNFdu);
  SmVector3d sUnitNFdv = m_dOrientations[0]*sNF.UnitizedDerivative(sNFdv);
  SmVector3d sUnitNGdu = m_dOrientations[1]*sNG.UnitizedDerivative(sNGdu);
  SmVector3d sUnitNGdv = m_dOrientations[1]*sNG.UnitizedDerivative(sNGdv);

  // orient and normalize SurfaceNormal vectors
  sNF = sNF*m_dOrientations[0]; SER(sNF.Unitize());
  sNG = sNG*m_dOrientations[1]; SER(sNG.Unitize());

  // set sDiff1 = vector between current rail points (desire sDiff1.Length == m_dDistance
  SmVector3d sDiff1 = sF - sG;

  // set sDiff2 = vector between points on the fillet center curve (desire sDiff2.Length == 0.0)
  SmVector3d sDiff2 = sDiff1 + dRadius*(sNF-sNG);

#ifdef SM_DEBUG_CODE
  if (bDebugMe) {
      sm_GraphicsLoop();
      smgfx_SetColor(0,0,1);
      SmPoint3d sPF = sF + dRadius*sNF;
      sPF.Draw();
      sm_GraphicsLoop();
      smgfx_SetColor(0,1,1);
      SmPoint3d sPG = sG + dRadius*sNG;
      sPG.Draw();
      sm_GraphicsLoop();
      smgfx_SetColor(1,0,0);
      sDiff1.Draw(&sPG);
      sm_GraphicsLoop();
  }
#endif

  // Intersection simply tries to zero the following four equations.
  // Let L be the distance between points.
  //
  // fun[i]   = (f-g) + R(Nf-Ng) . X     - X component
  // fun[i+1] = (f-g) + R(Nf-Ng) . Y     - Y component
  // fun[i+2] = (f-g) + R(Nf-Ng) . Z     - Z component
  // fun[i+3] = (f-g)(f-g) - L*L  - Distance between points is given distance
  //
  // Where f is surface1 original, g is surface2 original
  // and R is offset distance(or radius),
  // Nf is the oriented unitized normal of f,
  // Ng is the oriented unitized normal of g.
  //
  // It produces the following Jacobian with variables U, V, U2, V2, R
  //  | duf.X+R(d|Nf|/du).X   dvf.X+R(d|Nf|/dv).X   -dug.X-R(d|Ng|/du).X   -dvg.X-R(d|Ng|/dv).X   Nf.x-Ng.x |
  //  | duf.Y+R(d|Nf|/du).Y   dvf.Y+R(d|Nf|/dv).Y   -dug.Y-R(d|Ng|/du).Y   -dvg.Y-R(d|Ng|/dv).Y   Nf.y-Ng.y |
  //  | duf.Z+R(d|Nf|/du).Z   dvf.Y+R(d|Nf|/dv).Z   -dug.Y-R(d|Ng|/du).Z   -dvg.Y-R(d|Ng|/dv).Z   Nf.z-Ng.z |
  //  | 2*(f-g)*duf           2*(f-g)*dvf           -2*(f-g)*dug           -2*(f-g)*dvg           0.0       |
  //

  if (pOptJacobian) {
      ULONG i = rlNumEquations;
      ULONG lLastColumn = smos_Max(lSurf2Off,lSurf1Off) + 2;

      (*pOptJacobian)[i][lSurf1Off]   =   sDUF.x + dRadius*sUnitNFdu.x;
      (*pOptJacobian)[i][lSurf1Off+1] =   sDVF.x + dRadius*sUnitNFdv.x;
      (*pOptJacobian)[i][lSurf2Off]   = - sDUG.x - dRadius*sUnitNGdu.x;
      (*pOptJacobian)[i][lSurf2Off+1] = - sDVG.x - dRadius*sUnitNGdv.x;
      (*pOptJacobian)[i][lLastColumn] =   sNF.x  - sNG.x;


      i++; // Move to next equation;
      (*pOptJacobian)[i][lSurf1Off]   =   sDUF.y + dRadius*sUnitNFdu.y;
      (*pOptJacobian)[i][lSurf1Off+1] =   sDVF.y + dRadius*sUnitNFdv.y;
      (*pOptJacobian)[i][lSurf2Off]   = - sDUG.y - dRadius*sUnitNGdu.y;
      (*pOptJacobian)[i][lSurf2Off+1] = - sDVG.y - dRadius*sUnitNGdv.y;
      (*pOptJacobian)[i][lLastColumn] =   sNF.y  - sNG.y;

      i++;
      (*pOptJacobian)[i][lSurf1Off]   =   sDUF.z + dRadius*sUnitNFdu.z;
      (*pOptJacobian)[i][lSurf1Off+1] =   sDVF.z + dRadius*sUnitNFdv.z;
      (*pOptJacobian)[i][lSurf2Off]   = - sDUG.z - dRadius*sUnitNGdu.z;
      (*pOptJacobian)[i][lSurf2Off+1] = - sDVG.z - dRadius*sUnitNGdv.z;
      (*pOptJacobian)[i][lLastColumn] =   sNF.z  - sNG.z;

      i++;
      (*pOptJacobian)[i][lSurf1Off]   =  2.0 * sDiff1.Dot(sDUF);
      (*pOptJacobian)[i][lSurf1Off+1] =  2.0 * sDiff1.Dot(sDVF);
      (*pOptJacobian)[i][lSurf2Off]   = -2.0 * sDiff1.Dot(sDUG);
      (*pOptJacobian)[i][lSurf2Off+1] = -2.0 * sDiff1.Dot(sDVG);
      (*pOptJacobian)[i][lLastColumn] =  0.0;
  }

  // Compute function values - distance vector between offset points on fillet center curve
  //                           distance between points on rail curves
  rF[rlNumEquations++] = sDiff2.x;
  rF[rlNumEquations++] = sDiff2.y;
  rF[rlNumEquations++] = sDiff2.z;
  double dValue        = sDiff1.Dot(sDiff1) - m_dDistance*m_dDistance;
  rF[rlNumEquations++] = dValue;

  // failure: railPoints more than 10*m_Distance apart
  if (smos_Sqrt(smos_Fabs(dValue)) > 10.0 * m_dDistance)
    { return SM_ERR;  } // No convergance just skip out of here

  // set offset surface distances = current radius solution
  GetSurface(0)->SetOffsetDistance(dRadius*m_dOrientations[0]);
  GetSurface(1)->SetOffsetDistance(dRadius*m_dOrientations[1]);

  // See if we have converged
  double dScaledTol = GetConversionTol() * (1.0 + sF.GetMaxDimension());
  if (   smos_Fabs(sDiff2.x) < dScaledTol
      && smos_Fabs(sDiff2.y) < dScaledTol
      && smos_Fabs(sDiff2.z) < dScaledTol
      && smos_Fabs(dValue)   < dScaledTol)
    {
      //smos_Fabs(rF[rlNumEquations-1]) < dScaledTol)
      //smos_Fabs(sDiff1.Length()-m_dDistance) < dScaledTol)
      rbFoundAnswer = TRUE;
    }

#ifdef SM_DEBUG_CODE
  if (bDebugMe) {
      if (pOptJacobian) pOptJacobian->Dump();
      smos_WriteBuffer(_T(" F --- \n"));
      SM_DUMP_TARRAY( rF );
  }
#endif

  // all done
  return SM_SUCCESS;

} // end SmConstantDistanceFS::LoadJacobian

/*******************************************************************//**
PURPOSE: Get the number of variables the Jacobian for the constant
     radius fillet contains.

NOTES:
***********************************************************************/
ULONG SmConstantDistanceFS::GetJacobianSize
  ()
 const
{
    return 5;

} // end SmConstantDistanceFS::GetJacobianSize

/*******************************************************************//**
PURPOSE: Load the initial values for SmConstantDistanceFS.

NOTES:
***********************************************************************/
SmStatus SmConstantDistanceFS::LoadInitialValues
  (ULONG               & rbDoSurf2Calcs,     // in :
   ULONG                 lRailIndex,         // in :
   const SmSurface     & crSurface1,         // in :
   ULONG                 lSurface1Index,     // in :
   SmSurface          *& rpSurface2,         // out:
   ULONG               & rlSurface2Index,    // out:
   SmExtentNd          & rIntervals,         // out:
   SmTArray<SmBoolean> & rPeriodicities,     // out:
   SmTArray<double>    & rGuessT)            // i/o:
{
  if ( rbDoSurf2Calcs != 0 )
    {
      // Now find the UV guess on the other Rail's surface
      SmPoint3d sPnt;
      SmPoint2d sUV(rGuessT[lSurface1Index],rGuessT[lSurface1Index+1]);
      SER(crSurface1.EvaluatePoint(sUV,sPnt));
#ifdef SM_DEBUG_CODE
SmBoolean bDebugMe = FALSE;
      if (bDebugMe) {
          smgfx_SetColor(0,0,1);
          sPnt.Draw();
          sm_GraphicsLoop();
          smgfx_SetColor(1,0,0);
          SmPoint3d s3DPnt;
          SER(((SmOffsetSurface&)crSurface1).GetBaseSurface()->EvaluatePoint(sUV,s3DPnt));
          s3DPnt.Draw();
          sm_GraphicsLoop();
      }
#endif

      ULONG lOtherRailIndex = 1 - lRailIndex;
      rpSurface2            = GetSurface(lOtherRailIndex);
      SmPoint3d s2DPnt;

      SmSolution      sSData[16];
      SmSolutionArray sSolutions(16,sSData);

      const SmSurface *pSrf2Base = ((SmOffsetSurface*)rpSurface2)->GetBaseSurface();

      SmBoolean bTmp0 = FALSE ; 
      SmBSplineSurface * pBSplineSurface0 = SM_CAST_PTR(SmBSplineSurface, pSrf2Base) ; 
      SmTemporaryChangeValue<SmBoolean> sClean0(pBSplineSurface0 ? pBSplineSurface0->GetOutOfBoundsEnabled() : bTmp0, TRUE) ;

      // If we have edgeuses then use them to compute a point on the other surfaces
      // If not then drop the point to the other surface.
      SmBoolean bFoundUV = FALSE;
      SmEdgeuse *pOtherEU = m_pEdgeuses[lOtherRailIndex];
      if ( pOtherEU != NULL && pOtherEU->GetUVTrimCurve() == NULL )
        { pOtherEU = NULL; }

      // But: if we are rolling on a different face (e.g., doing a rollover),
      // then this is the wrong Edgeuse, and we don't have an appropriate one to use.
      if ( pOtherEU != NULL )
        {
          SmFace *pF = pOtherEU->GetFace();
          if ( pF != NULL && pF->GetSurface() != pSrf2Base )
            { pOtherEU = NULL; }
        }

      if ( pOtherEU != NULL ) 
        {
          SmCurve    *p3DFilletEdgeCurve = pOtherEU->GetEdge()->GetCurve();
          SmExtent1d  sIvl               = pOtherEU->GetEdge()->GetInterval();
          SER(p3DFilletEdgeCurve->GlobalPointSolve(sIvl,         SM_SO_MINIMIZE,
                                                   sPnt,         m_dThisApproxTol3d,
                                                   NULL,         NULL,
                                                   SM_SR_SINGLE, sSolutions));
          if (sSolutions.GetSize() != 1) SER(SM_ERR);
          double dFilletEdgeT = sSolutions[0].m_vStart[0];

          // If the curve is closed, and the solution is within tol of an
          // end point, snap the solution to whichever end (start or finish)
          // such that the direction that moves towards the interior of the
          // curve also moves generally towards the interior of the surface.

          if (p3DFilletEdgeCurve->IsClosed(sIvl))
            {
              SmPoint3d s3DPnt;
              SER( p3DFilletEdgeCurve->EvaluatePoint( dFilletEdgeT, s3DPnt ));
              SmVector3d sInwardVector;
              SER( crSurface1.ComputeInwardVector( sUV, sInwardVector ));
              double dDist, dParam;
              SmBoolean bSuccess;
              SER(p3DFilletEdgeCurve->DropPoint(sIvl,                 // in : target curve allowed domain
                                                s3DPnt,               // in : Point to drop to curve
                                                &sInwardVector,       // in : For Snapping, NULL=NoSnap, Vec pointing to CrvInside from desired EndPt for snapping.
                                                                      //      when DropPt is within dDistTol in cpVecToCrvInside dir of EndPt, Snap to Ivl EndPt.
                                                                      //      Vec handles ambiguities and makes sure that curves are not just touching at the ends.
                                                m_dThisApproxTol3d,   // in : Used for SnapDist, and SM_SO_INTERSECT sXSectTol3d distance.
                                                                      //      sXSectTol3d = dDistTol ;      SM_SO_INTERSECT 
                                                                      //      sXSectTol3d = SM_BIG_DOUBLE ; SM_SO_MINIMIZE, SM_SO_NORMALIZE
                                                                      //      sSnapTol3d  = dDistTol ;      SM_SO_MINIMIZE, SM_SO_NORMALIZE, SM_SO_INTERSECT    
                                                NULL,                 // in : If given, call cheap LocalPointSolve() else call expensive GlobalPointSolve()
                                                bSuccess,             // out: TRUE = found a drop point
                                                dParam,               // out: found drop curve param
                                                dDist)) ;             // out: found drop distance
                                                                      // in : SM_SO_MINIMIZE = allow nonNormal drops near endPoints
                                                                      //      SM_SO_INTERSECT= like Minimize but point must be within dDistTol
                                                                      //      SM_SO_NORMALIZE= exclude nonNormal drops near endPoints
                                                                      //      default:[SM_SO_MINIMIZE] to preserve original behavior
              if ( bSuccess )
                  dFilletEdgeT = dParam;
            }

          // Now determine edgeuses, surfaces and orientations for current
          // iteration of the solver.
          SmBSplineCurve *pOtherUVCurve = pOtherEU->GetUVTrimCurve();
          SER(pOtherUVCurve->EvaluatePoint(dFilletEdgeT,s2DPnt));
          bFoundUV = TRUE;
          SmSolution sSol;
          SmBoolean bFoundAnswer;
          SmPoint2d sGuess(s2DPnt.x,s2DPnt.y);
          SER(rpSurface2->LocalPointSolve(rpSurface2->GetNaturalUVDomain(),
              SM_SO_MINIMIZE,sPnt,sGuess,bFoundAnswer,sSol));

          if (bFoundAnswer)
            {
              double dDist = sSol.m_vStart.m_dSolutionValue;
              if (dDist < 2.0*m_dDistance)
                {
                  s2DPnt.x = sSol.m_vStart[0];
                  s2DPnt.y = sSol.m_vStart[1];
#ifdef SM_DEBUG_CODE
SmBoolean bDebugMe1 = FALSE;
                  if (bDebugMe1) {
                      SmPoint2d sUVTest(s2DPnt.x,s2DPnt.y);
                      SmPoint3d s3DPnt;
                      SER(pSrf2Base->EvaluatePoint(sUVTest,s3DPnt));
                      smgfx_SetColor(0,1,0);
                      s3DPnt.Draw();
                      sm_GraphicsLoop();
                      pSrf2Base->DrawUV(1,1);
                      sm_GraphicsLoop();
                  }
#endif
                }
            }
        }

      if (bFoundUV)
        {
          SmVector3d sDU, sDV;
          SmVector3d sMat[2][2];
          SmPoint2d sUVPnt(s2DPnt.x,s2DPnt.y);
          SER(pSrf2Base->Evaluate(sUVPnt,1,1,TRUE,TRUE,TRUE,sMat[0]));
          sDU = sMat[1][0];
          sDV = sMat[0][1];
          double dDeg;
          SER(sDU.AngleBetween(sDV,dDeg));
          dDeg = SM_RAD2DEG(dDeg);
          if (dDeg < 30 || dDeg > 150.0) {
              bFoundUV = FALSE;
          }
          // If we are near a singularity then use global point solve
          if (sDU.GetMaxDimension() > 10000.0 * sDV.GetMaxDimension() ||
              sDV.GetMaxDimension() > 10000.0 * sDU.GetMaxDimension()) {
              bFoundUV = FALSE;
          }
        }

      if (!bFoundUV)
        {
          SER(rpSurface2->GlobalPointSolve(rpSurface2->GetNaturalUVDomain(),
              SM_SO_MINIMIZE,sPnt,SM_EFF_ZERO_SQ,
              NULL,SM_SR_SINGLE,sSolutions));
          if (sSolutions.GetSize() != 1) SER(SM_ERR);

          s2DPnt.x = sSolutions[0].m_vStart[0];
          s2DPnt.y = sSolutions[0].m_vStart[1];
        }

      SmExtent2d sDomain = GetSurface(lOtherRailIndex)->GetBaseSurface()->GetNaturalUVDomain();
      SmVector2d sSize = GetDomainExtension( GetSurface(lOtherRailIndex)->GetBaseSurface(), sDomain ); // SolidCraft change: and the trace's reach

      rIntervals[rGuessT.GetSize()] = SmExtent1d(sDomain.GetMin().x-sSize.x,sDomain.GetMax().x+sSize.x);
      rIntervals[rGuessT.GetSize()+1] = SmExtent1d(sDomain.GetMin().y-sSize.y,sDomain.GetMax().y+sSize.y);

      rPeriodicities.Add(FALSE);
      rPeriodicities.Add(FALSE);

      rlSurface2Index = rGuessT.GetSize();

      rGuessT.Add(s2DPnt.x);
      rGuessT.Add(s2DPnt.y);
  }

  // Last of all load the distance
  rIntervals[rlSurface2Index+2] = SmExtent1d(0.0,m_dDistance*10.0);
  rPeriodicities.Add(FALSE);
  double dApproxRadius          = m_dDistance * smos_Sqrt(2.0) / 2.0;
  GetSurface(0)->SetOffsetDistance(dApproxRadius*m_dOrientations[0]);
  GetSurface(1)->SetOffsetDistance(dApproxRadius*m_dOrientations[1]);

  rGuessT.Add(dApproxRadius);

  return SM_SUCCESS;

} // end SmConstantDistanceFS::LoadInitialValues

/*******************************************************************//**
PURPOSE: Utility methods.

NOTES:
***********************************************************************/
SmBoolean SmConstantDistanceFS::IsKindOf( SM_TYPE t ) const
{
  return ((SmConstantDistanceFS_TYPE == t) ? TRUE : SmFilletSolver::IsKindOf( (t) ));
}
/*******************************************************************//**
PURPOSE: Dump a SmConstantDistanceFS.

NOTES: Unimplemented.
***********************************************************************/
void SmConstantDistanceFS::Dump() const {}


/*******************************************************************//**
PURPOSE: Constructor for Variable Radius Fillet object.

NOTES:
***********************************************************************/
SmVariableRadiusFS::SmVariableRadiusFS
  (const SmContext & crContext,
   double            dThisApproxTol3d,
   double            dAngleTolerance,
   double            dTangencyTolerance,
   double         /* dDistance */,       // Euclidian distance between the two rails
   SmEdgeuse       * pEdgeuse,
   SmFilletLaw     & rFilletLaw,
   SmBoolean         bLawOrientation)
 : SmCurveBasedFS(crContext, dThisApproxTol3d, dAngleTolerance, dTangencyTolerance),
   m_rLaw(&rFilletLaw),
   m_bLawOrientation(bLawOrientation)
{
    SM_ASSERT(pEdgeuse != NULL);
    m_pEdgeuses[0] = pEdgeuse;
    m_pEdgeuses[1] = pEdgeuse->GetRadial();

    SmEdge *pEdge = pEdgeuse->GetEdge();
    m_vLawInterval = pEdge->GetInterval();

    m_pOriginal = pEdge->GetCurve();

    if (m_pOriginal->GetDegree() > 1) {
        double dAccuracy = dThisApproxTol3d*10.0;
        double dEdgeArcLength = 0.0;
        SE(m_pOriginal->Length(m_vLawInterval,dAccuracy,dEdgeArcLength));
        rFilletLaw.SetEdgeArcLength(dEdgeArcLength);
    }

    // Create offset surfaces of appropriate radius in appropriate
    // direction.  Note that we may need to adjust this if the two
    // surfaces go to tangency and then back down in the opposite direction.
    // For now assume that the orientations at the center of the edge
    // work for the rest of the edge.

    double dVals[3];
    SE(rFilletLaw.Evaluate(m_vLawInterval.Evaluate(0.5),m_vLawInterval,m_bLawOrientation,dVals));
    double dApproxRadius = dVals[0];

    // SolidCraft change: offset each surface to the side of it the other face lies on
    // (sm_SectorSides), a cusp's faces too.
    double dSides[2] = { 1.0, 1.0 };
    if (sm_SectorSides(m_pEdgeuses, smos_Fabs(dApproxRadius), dSides) != SM_SUCCESS)
    {
        // Smooth, or no sector: can't fillet this.
        m_pEdgeuses[0] = NULL;
        m_pEdgeuses[1] = NULL;
        SE(SM_ERR);
        return;
    }
    m_dOrientations[0] = dSides[0];
    m_dOrientations[1] = dSides[1];

    SmSurface *pSurface = m_pEdgeuses[0]->GetLoopuse()->GetFaceuse()->GetFace()->GetSurface();
    SmSurface *pSurface2 = m_pEdgeuses[1]->GetLoopuse()->GetFaceuse()->GetFace()->GetSurface();
    m_pSurfaces[0] = new (crContext) SmOffsetSurface(dApproxRadius*m_dOrientations[0],*pSurface);
    m_pSurfaces[1] = new (crContext) SmOffsetSurface(dApproxRadius*m_dOrientations[1],*pSurface2);
    m_pSurfaces[0]->SetAllowExtension(TRUE);
    m_pSurfaces[1]->SetAllowExtension(TRUE);

    (void)FinishConstruction();  // Take care of other things in the base class. // SolidCraft: a constructor cannot report

#ifdef SM_DEBUG_CODE
SmBoolean bDebugMe = FALSE;
    if (bDebugMe) {
        sm_GraphicsLoop();
        smgfx_SetColor(0,0,0);
        pSurface->DrawUV(4,8);
        sm_GraphicsLoop();
        smgfx_SetColor(0,0,1);
        pSurface2->DrawUV(9,18);
        sm_GraphicsLoop();
    }
#endif


} // end SmVariableRadiusFS::SmVariableRadiusFS constructor

/*******************************************************************//**
PURPOSE: Constructor for Variable Radius Fillet object.

NOTES:
***********************************************************************/
SmVariableRadiusFS::SmVariableRadiusFS
  (const SmContext & crContext,
   double            dThisApproxTol3d,
   double            dAngleTolerance,
   double            dTangencyTolerance,
   SmSurface       * pSurface,
   SmSurface       * pSurface2,
   double            dOffsetOrientation,
   double            dOffsetOrientation2,
   SmCurve         * pCenterCurve,
   SmFilletLaw     & rFilletLaw,
   SmBoolean         bLawOrientation)
 : SmCurveBasedFS(crContext,
                  dThisApproxTol3d,
                  dAngleTolerance,
                  dTangencyTolerance),
   m_rLaw(&rFilletLaw),
   m_bLawOrientation(bLawOrientation)
{
  m_pOriginal    = pCenterCurve;
  m_pEdgeuses[0] = NULL;
  m_pEdgeuses[1] = NULL;

  m_vLawInterval = m_pOriginal->GetNaturalInterval();

  if (m_pOriginal->GetDegree() > 1) 
    {
      double dAccuracy = dThisApproxTol3d*10.0;
      double dEdgeArcLength = 0.0;
      SE(m_pOriginal->Length(m_vLawInterval,dAccuracy,dEdgeArcLength));
      rFilletLaw.SetEdgeArcLength(dEdgeArcLength);
    }

  // Create offset surfaces of appropriate radius in appropriate
  // direction.
  m_dOrientations[0] = dOffsetOrientation;
  m_dOrientations[1] = dOffsetOrientation2;

  double dVals[3];
  SE(rFilletLaw.Evaluate(m_vLawInterval.Evaluate(0.5),m_vLawInterval,m_bLawOrientation,dVals));
  double dApproxRadius = dVals[0];

  m_pSurfaces[0] = new (crContext) SmOffsetSurface(dApproxRadius*m_dOrientations[0],*pSurface);
  m_pSurfaces[1] = new (crContext) SmOffsetSurface(dApproxRadius*m_dOrientations[1],*pSurface2);
  m_pSurfaces[0]->SetAllowExtension(TRUE);
  m_pSurfaces[1]->SetAllowExtension(TRUE);
#ifdef SM_DEBUG_CODE
SmBoolean bDebugMe = FALSE;
  if (bDebugMe) 
    {
      sm_GraphicsLoop();
      smgfx_SetLook(1,2, 0,0,0); pSurface->DrawUV(4,8); sm_GraphicsLoop();
      smgfx_SetLook(1,2, 0,0,1); pSurface2->DrawUV(9,18); sm_GraphicsLoop();
      sm_GraphicsLoop();
    }
#endif // SM_DEBUG_CODE

} // end SmVariableRadiusFS::SmVariableRadiusFS constructor

/*******************************************************************//**
PURPOSE: Destructor for Variable Radius Distance fillet solver.

NOTES:
***********************************************************************/
SmVariableRadiusFS::~SmVariableRadiusFS
  ()
{

} // end SmVariableRadiusFS::~SmVariableRadiusFS destructor


/*******************************************************************//**
PURPOSE: Compute fillet radius at a given edge parameter.

NOTES:
***********************************************************************/
double SmVariableRadiusFS::GetFilletRadius
  (double dParam)
{
    double dValues[3];
    SER(m_rLaw->Evaluate(dParam,m_vLawInterval,m_bLawOrientation,dValues));
    return dValues[0];

} // end SmVariableRadiusFS::GetFilletRadius


/*******************************************************************//**
PURPOSE: Recreate the fillet surface (and fillet edges, if needed)
    Typically, this routine can be used to refine the fillet surfaces
    such that it will interpolate all the rail ends. (For example, when
    SmFilletGeom were split, it is often desirable to have the fillet
    surfaces pass through their rail ends) Another use of it is to
    compute the fillet surface patch for each SmFilletGeom when roll-over
    occurred.

NOTES: User may want to specify 'optional' pMarchDir if the
    geometry hasn't yet been computed
***********************************************************************/
SmStatus SmVariableRadiusFS::ReCalcFilletGeom
 (SmBoolean    bIsAnalyticFillet,  // in : 
  SmVector3d * pOptMarchDir,       // NotUsed: in : 
  SmCurve    * pOptRefCurve)       // NotUsed: in : 
{
  SM_REF2(bIsAnalyticFillet, pOptMarchDir) ;
  const SmContext & crContext = GetCreationContext();
    m_bExtendBefore = FALSE;
    m_bExtendAfter = FALSE;
    SmFilletGeom * pFilletGeom = m_vFilletGeoms[0];

    SmCurve * pCurve = NULL;
    SmBoolean b3DTest = TRUE;
    SmObjDelete sDelete;
//    if (m_eType == SM_FG_CORNER_FILLET)
    if (1) {
        pCurve = pOptRefCurve; NER(pCurve);
    }
    // else 
    // {
    //    // if (m_eType != SM_FG_TANGENT_ROLLOVER)
    //    //     b3DTest = FALSE;
    // 
    //     SmEdge * pEdge = GetEdgeuse(0)->GetEdge();
    //     pCurve = pEdge->GetCurve(); NER(pCurve);
    // }

    SmExtent1d sIvl = pCurve->GetNaturalInterval();
    SmBoolean bIsClosedSolverEdge = FALSE;
#if 0
    if (0 &&//m_eType != SM_FG_CORNER_FILLET && pCurve->IsClosed(sIvl)) 
    {
        bIsClosedSolverEdge = TRUE;
    }
#endif

#ifdef SM_DEBUG_CODE
SmBoolean bDebugMe0=FALSE;
    if (bDebugMe0) {
        smgfx_SetColor(0,1,0);
        pCurve->DrawWDeriv(sIvl);
        sm_GraphicsLoop();
    }
#endif
    SmOffsetSurface * pSurf1 = GetSurface(0);
    SmOffsetSurface * pSurf2 = GetSurface(1);
    double dTol = m_dThisApproxTol3d;
    //double dAngTol = GetThisAngTolRad();
    double dScale = sIvl.GetLength();


    // Sort the SmTsectPnt's of each rail-end, the rails should
    // go in the same 'direction' as the solver edge
    // First, determine if the sorting is in 2D or 3D.

    SmTsectPnt* sP1Data[64];
    SmTArray<SmTsectPnt*> sTsectPnts(64,sP1Data);
    SmTArray<double> sProjParams;
    for (ULONG i=0; i<2; i++) {
        if (!b3DTest) {
            SmEdgeuse * pEUOnSurf = GetEdgeuse(i); NER(pEUOnSurf);
            pCurve = pEUOnSurf->GetUVTrimCurve(); NER(pCurve);
        }

        SmTArray<SmFilletEdge *> sAllRails;
        sAllRails.Add(pFilletGeom->GetRail(i));
        if (i==0) {
//            sAllRails.Append(pFilletGeom->m_vOtherRails1);
        }
        else {
//            sAllRails.Append(pFilletGeom->m_vOtherRails2);
        }
        for (ULONG kk=0; kk<sAllRails.GetSize(); kk++) {
            SmFilletEdge * pRail= sAllRails[kk]; NER(pRail);
            SmTArray<SmEdgeuse*> sFilEUs;
            SER( pRail->GetEdgeuses(sFilEUs) );
            if (sFilEUs.GetSize() != 2) SER(SM_ERR);
            for (ULONG jjj=0; jjj<2; jjj++) {
                SmEdgeuse * pEU = sFilEUs[jjj];
                SmFilletVertexuse * pVU = (SmFilletVertexuse*)pEU->GetVertexuse();
                SmFilletVertex * pFV = (SmFilletVertex*)pVU->GetVertex();
                if (pFV->GetStatus() == SM_FIL_UNPROCESSED ||
                    pFV->GetFilletVertexType() == SM_FV_MATE ||
                    pFV->GetFilletVertexType() == SM_FV_RAIL_X_EXTENDED_EDGEUSE) {
                    continue;
                }
                SmTsectPnt & rTsectPnt = pVU->GetTsectPnt();
                SmPoint3d sPnt;
                if (!b3DTest) { // Do a 2D test
                    sPnt = SmPoint3d(rTsectPnt.UVPos(i));
                    if (sPnt.x == -SM_BIG_DOUBLE || sPnt.y == -SM_BIG_DOUBLE) {
                        continue;
                    }
                }
                else { // Do a 3D test
                    if (rTsectPnt.m_ePointType == SM_IP_UNKNOWN) {
//                        SmPoint2d sUV = rTsectPnt.UVPos(0);
//                        SER(SetupOffsetValues(rTsectPnt));
//                        SER(pSurf1->EvaluatePoint(sUV,rTsectPnt.CrvPos()));
                    }
                    sPnt = rTsectPnt.CrvPos();
                }

                // Drop 3D point on surface
                SmSolution sSData[8];
                SmSolutionArray sSolutions(8,sSData);
                SER(pCurve->GlobalPointSolve(sIvl,SM_SO_MINIMIZE,
                    sPnt,dTol,NULL,NULL,SM_SR_SINGLE,sSolutions));
                if (sSolutions.GetSize() != 1) SER(SM_ERR);
                double dT = sSolutions[0].m_vStart[0];
#ifdef SM_DEBUG_CODE
SmBoolean bDebugMe1=FALSE;
                if (bDebugMe1) {
                    smgfx_SetColor(0,1,0);
                    SmPoint3d sP1, sP2;
                    SmPoint2d sUV = rTsectPnt.UVPos(0);
                    GetSurface(0)->GetBaseSurface()->EvaluatePoint(sUV,sP1);
                    smgfx_SetPointSize(6);
                    sP1.Draw();
                    if (FALSE) {
                        GetSurface(0)->DrawUV(1,1);
                    }
                    sm_GraphicsLoop();
                    SmPoint2d sUV2 = rTsectPnt.UVPos(1);
                    GetSurface(1)->GetBaseSurface()->EvaluatePoint(sUV2,sP2);
                    smgfx_SetColor(1,0,0);
                    sP2.Draw();
                    if (FALSE) {
                        GetSurface(1)->DrawUV(1,1);
                    }
                    sm_GraphicsLoop();
                    rTsectPnt.CrvPos().Draw();
                    sm_GraphicsLoop();
                }
#endif
                //if (sIvl.IsValueOnBoundary(dT) && !b3DTest) {
                //    SmPoint3d sPV[2];
                //    SER(pUVCurveOnSurf->Evaluate(dT,1,1,sPV));
                //    double dAdjustment;
                //    SER(smgu_LineClosestPoint(sPV[0],sPV[1],sUVPnt,dAdjustment));
                //    dT += dAdjustment;
                //}
                if (b3DTest && sIvl.IsValueOnBoundary(dT) && bIsClosedSolverEdge) {
                    // Do this logical test if dT is on boundary of a
                    // closed edge. (Please noted that the rail is already
                    // in the same direction as the solver)
                    if (pEU->GetOrientation() == SM_OT_OPPOSITE) {
                        dT = sIvl.GetMax();
                    }
                }

                ULONG lIndexOfInsert = sProjParams.GetSize();
                SmBoolean bDoInsertion = TRUE;
                for (ULONG k=0; k<sProjParams.GetSize(); k++) {
                    if (smos_Fabs(dT-sProjParams[k]) <
                        SM_EFF_ZERO_SQRT*dScale) {
                        bDoInsertion = FALSE;
                        break;
                    }
                    else if (dT < sProjParams[k]) {
                        lIndexOfInsert = k;
                        break;
                    }
                }
                if (bDoInsertion) {
                    sProjParams.InsertAt(lIndexOfInsert, dT);
                    //rTsectPnt.m_dCurveParameter = -SM_BIG_DOUBLE;
                    sTsectPnts.InsertAt(lIndexOfInsert, &rTsectPnt);
                }
            } // for-jjj
        } // for-kk
    } // for-i

    SmExtent2d sDomain1 = pSurf1->GetNaturalUVDomain();
    SmExtent2d sDomain2 = pSurf2->GetNaturalUVDomain();
    //        if (m_eType != SM_FG_CORNER_FILLET) {
    SetupOffsetExtension(pSurf1,sDomain1);
    SetupOffsetExtension(pSurf2,sDomain2);
    //        }
    // Setup FilletIntersector
    SmFilletIntersector sFI(*pSurf1,sDomain1,
                            *pSurf2,sDomain2,
                            *this);


    SmBSplineCurve *p3DCurve = NULL ;
    SmTsectCurveType eCurveType;
    double dDeviation = 0.0;
    SmBSplineCurve *pSurfaceUVCurves[2];
//    SER(sFI.DoPointIntersection(crContext,sTsectPnts,m_bExtendBefore,
//        m_bExtendAfter,pOptMarchDir,NULL,&dTol,&dAngTol,p3DCurve,pSurfaceUVCurves[0],
//        pSurfaceUVCurves[1], eCurveType, dDeviation));
    SmBoolean bUseNormalAveraging = TRUE;
    SmBoolean bUniformSteps = TRUE;
    SER(sFI.DoLawIntersection(crContext,sTsectPnts,
                               *m_rLaw,
                               m_bLawOrientation,
                               *m_pOriginal,
                               m_vLawInterval,
                               m_dOrientations,
                               bUseNormalAveraging,
                               bUniformSteps,
                               SM_CAST_APPROXTOL3D_PTR(&m_dThisApproxTol3d),
                               &m_dThisAngTolRad,
                               p3DCurve,
                               pSurfaceUVCurves[0],
                               pSurfaceUVCurves[1],
                               eCurveType,
                               dDeviation));
    for (ULONG ii=1; ii<sTsectPnts.GetSize(); ii++) {
        SmTsectPnt * pTSP = sTsectPnts[ii];
        if (pTSP->m_dCurveParameter < 0.0) {
            return SM_ERR;
        }
    }


#ifdef SM_DEBUG_CODE
SmBoolean bDebugMe=FALSE;
    if (bDebugMe) {
        smgfx_SetColor(1,0,0);
        pFilletGeom->GetFilletSurface()->DrawUV(1,1);
        sm_GraphicsLoop();
    }
#endif

    return SM_SUCCESS;

} // end SmVariableRadiusFS::ReCalcFilletGeom

/*******************************************************************//**
PURPOSE: Check whether the radius can change.

NOTES:
  This must return False if it's actually constant.  [090511]
***********************************************************************/
SmBoolean SmVariableRadiusFS::OffsetRadiiCanChange
  ()
{
  SmBoolean bRet = TRUE;
  if ( m_pFSG->GetFilletSurfaceGeneratorType() == SM_FSG_CIRCULAR
      && m_rLaw->IsLinearLaw() == TRUE )
    {
      SmLinearFilletLaw *pLinearLaw = SM_CAST_PTR( SmLinearFilletLaw, m_rLaw );
      if ( pLinearLaw != NULL )
        {
          if ( smos_Fabs(  pLinearLaw->GetStartRadius()
                         - pLinearLaw->GetEndRadius() ) < SM_EFF_ZERO )
            { bRet = FALSE; }
        }
    }

  return bRet;

} // end SmVariableRadiusFS::OffsetRadiiCanChange

/*******************************************************************//**
PURPOSE: The edge parameters where the radius law has a corner.

NOTES: SolidCraft change. Appended to rParams, ascending.
***********************************************************************/
SmStatus SmVariableRadiusFS::GetRadiusCornerParameters
  (SmTArray<double> & rParams)       // out: corner parameters on the filleted edge
{
  NER(m_rLaw);
  SER(m_rLaw->GetCornerParameters(m_vLawInterval, m_bLawOrientation, rParams));
  return SM_SUCCESS;

} // end SmVariableRadiusFS::GetRadiusCornerParameters

/*******************************************************************//**
PURPOSE: Whether an intersection point lies where the radius law has a corner.

NOTES: SolidCraft change. The point's law parameter is m_adUserDoubles[0]
  (see SetupOffsetValues).
***********************************************************************/
SmBoolean SmVariableRadiusFS::IsAtRadiusCorner
  (const SmTsectPnt & crTsectPnt)    // in : intersection point
{
  SmTArray<double> sCorners;
  if (GetRadiusCornerParameters(sCorners) != SM_SUCCESS)
    { return FALSE; }
  double dTol = 1.0e-9 * (1.0 + smos_Fabs(m_vLawInterval.GetMin()) + smos_Fabs(m_vLawInterval.GetMax()));
  for (ULONG i=0; i<sCorners.GetSize(); i++)
    {
      if (smos_Fabs(crTsectPnt.m_adUserDoubles[0] - sCorners[i]) <= dTol)
        { return TRUE; }
    }
  return FALSE;

} // end SmVariableRadiusFS::IsAtRadiusCorner

/*******************************************************************//**
PURPOSE: The edge parameter of an intersection point, when the radius law
  has corners.

NOTES: SolidCraft change. The point's law parameter is m_adUserDoubles[0].
***********************************************************************/
SmBoolean SmVariableRadiusFS::GetRadiusParameter
  (const SmTsectPnt & crTsectPnt,    // in : intersection point
   double           & rdParam)       // out: its edge parameter
{
  SmTArray<double> sCorners;
  if (GetRadiusCornerParameters(sCorners) != SM_SUCCESS || sCorners.GetSize() == 0)
    { return FALSE; }
  rdParam = crTsectPnt.m_adUserDoubles[0];
  return TRUE;

} // end SmVariableRadiusFS::GetRadiusParameter

/*******************************************************************//**
PURPOSE: Set extension flags to TRUE for each end if there is not
    a blend being done at the corner.

NOTES:
***********************************************************************/
void SmVariableRadiusFS::SetExtensionFlags
  ()
{
    SmFilletSolver::SetExtensionFlags();

    // In order to help out computing fillet radius outside the curve boundary
    // when the filleted edge is not linear. Let's extend the edge.
    SmBSplineCurve *pBSC = SM_CAST_PTR(SmBSplineCurve,m_pOriginal);
    if (pBSC) {
        double dAnalyticTol = 1.0e-8;
        SmPoint3d sLinePoint;
        SmVector3d sLineVector;
        if (   m_rLaw->GetExtendedEdgeCurve() == NULL
            && m_pOriginal->GetDegree() > 1
            && !pBSC->IsLine(5,dAnalyticTol,sLinePoint,sLineVector))
        {
            // Should come in here only once.
            SmExtent1d sIvl = m_pOriginal->GetNaturalInterval();
            SmBSplineCurve *pExtended = NULL;
            SmBoolean bPreciseExtension = TRUE;  // Otherwise it adds another length of sIvl.
            double dNewT = sIvl.GetMin() - sIvl.GetLength();
            static constexpr SmContinuityType eExtCont = SM_CT_CINFINITY;
            SE(pBSC->CreateExtendedCurve( m_crContext,
                dNewT, eExtCont, pExtended, bPreciseExtension ));
            pBSC = pExtended;
            SmObjDelete sDelObj;
            if (pExtended) {
                sDelObj.SetObj(pExtended);
            }
            dNewT = sIvl.GetMax() + sIvl.GetLength();
            SE(pBSC->CreateExtendedCurve( m_crContext,
                dNewT, eExtCont, pExtended, bPreciseExtension ));
            m_rLaw->SetExtendedEdgeCurve(pExtended);

            // Set ours as well.
            //  No, this gets deleted in m_rLaw's destructor.
            //  If we want to set this we must make a copy.
            //m_pExtended = pExtended;

#ifdef SM_DEBUG_CODE
SmBoolean bDebugMe = FALSE;
            if (bDebugMe) {
                smgfx_SetLook( 3,5, 0,1,1 ); if( m_pOriginal ) { m_pOriginal->DrawParams(); sm_GraphicsLoop(); }
                smgfx_SetLook( 4,6, 1,0,0 ); if(   pExtended ) {   pExtended->DrawParams(); sm_GraphicsLoop(); }
                sm_GraphicsLoop();
            }
#endif // SM_DEBUG_CODE

        } // end if we want to extend the Edge curve.
    }

} // end SmVariableRadiusFS::SetExtensionFlags

/*******************************************************************//**
PURPOSE: Set the offset values to correspond to what they should be
     at this intersection point.

NOTES:
***********************************************************************/
SmStatus SmVariableRadiusFS::SetupOffsetValues
  (const SmTsectPnt & rTsectPnt,     // in : different TsectPnts get different radius fillets
   double * pdOffsetDist)            // out: optional fillet radius at this point, NULL to ignore
{
  // compute the radius value at this intersection point
  // (A failed evaluation leaves dValues unset: stop rather than offset by whatever is there.)
  double dValues[3];
  SER(m_rLaw->Evaluate(rTsectPnt.m_adUserDoubles[0], // in : target parameter - range:[crCurveInterval.Min,Max]
                  m_vLawInterval,                // in : crCurveInterval = interval defining range of fillet edge
                  m_bLawOrientation,             // in : bReverseOrientation: TRUE = parameters run from interval end to interval start
                  dValues));                     // out: dValues[0] = fillet-radius at dParameter value
                                                 //      dValues[1] = 1st derivative: d(fillet-radius)/d(Parameter)
                                                 //      dValues[2] = 2nd derivative: d2(fillet-radius)/d(Parameter)**2

  // save fillet radius and set OffsetSurface offset distances
  if (pdOffsetDist) *pdOffsetDist = dValues[0];
  GetSurface(0)->SetOffsetDistance(dValues[0]*m_dOrientations[0]);
  GetSurface(1)->SetOffsetDistance(dValues[0]*m_dOrientations[1]);

  return SM_SUCCESS;

} // end SmVariableRadiusFS::SetupOffsetValues

/*******************************************************************//**
PURPOSE: SolidCraft change. The ball this solver rolls at a station: its
    radius law r(t) and the filleted edge's curve C(t) there, t the
    station's law parameter (m_adUserDoubles[0]).

NOTES: For SmFilletSolver::JudgeBallAt: the ball's centre c lies in the
    plane square to the edge at C(t), <C'(t), C(t) - c> = 0 (LoadJacobian).
***********************************************************************/
SmStatus SmVariableRadiusFS::GetBallLaw
  (const SmTsectPnt & crPnt,       // in : the station
   SmBoolean        & rbBall,      // out: TRUE
   SmBoolean        & rbOnEdge,    // out: TRUE
   double             adR[3],      // out: r, r', r'' at t
   SmVector3d         asC[4])      // out: C, C', C'', C''' at t
{
  NER(m_rLaw);
  const double dT = crPnt.m_adUserDoubles[0];
  SER(m_rLaw->Evaluate(dT, m_vLawInterval, m_bLawOrientation, adR));
  SER(EvaluateBaseCurve(dT, 3, asC));
  rbBall   = TRUE;
  rbOnEdge = TRUE;
  return SM_SUCCESS;
} // end SmVariableRadiusFS::GetBallLaw

/*******************************************************************//**
PURPOSE: Load the current Jacobian matrix in an incremental way
    and update the parameters.

NOTES:
  // Intersection simply tries to zero the following 4 equations.

  // fun[i+0] = (f-g) . X     - X component of vector f-g
  // fun[i+1] = (f-g) . Y     - Y component of vector f-g
  // fun[i+2] = (f-g) . Z     - Z component of vector f-g
  // fun[i+3] = C'(t) * (C(t) - f)
  //
  // Where f is surface1 offset, g is surface2 offset and
  // where f2 is surface1 original and g2 is surface2 original
  // And C(t) is the filleted edge curve: offset point is in the plane.
  //
  // It produces the following Jacobian with variables U, V, U2, V2, T
  //  | duf.X      dvf.X       -dug.X       -dvg.X      L'N.x-L'N2.x |
  //  | duf.Y      dvf.Y       -dug.Y       -dvg.Y      L'N.y-L'N2.y |
  //  | duf.Z      dvf.Z       -dug.Z       -dvg.Z      L'N.z-L'N2.z |
  //  | C'.-duf    C'.-dvf      0.0          0.0   C''C + C'C' - C''f - C'L'N |
  //
  //  N is the unitized normal of f2 and N2 is unitized normal of g2


  With: crX   = [u1 v1 u2 v2 t]
         F    = OffsetSurf1.Position(u1, v1) for variableRadius = variableRadiusLaw(t)
         G    = OffsetSurf2.Position(u2, v2) for variableRadius = variableRadiusLaw(t)
         t    = BaseCurveParameter - typically a parameter point on the edge being filleted
        C(t)  = BaseCurve.Position
        C'(t) = BaseCurve.Tangent

  Sets:

  rF[]                 = done elsewhere
  rF[rlNumEquations+0] = (F-G) . X     - X component of vector F-G
  rF[rlNumEquations+1] = (F-G) . Y     - Y component of vector F-G
  rF[rlNumEquations+2] = (F-G) . Z     - Z component of vector F-G
  rF[rlNumEquations+3] = C'(t) * (C(t) - F)

  pOptJacobian[][]                 = done elsewhere
  pOptJacobian[rlNumEquations+0][] = [ duf.X      dvf.X       -dug.X       -dvg.X      L'N.x-L'N2.x          ]
  pOptJacobian[rlNumEquations+1][] = [ duf.Y      dvf.Y       -dug.Y       -dvg.Y      L'N.y-L'N2.y          ]
  pOptJacobian[rlNumEquations+2][] = [ duf.Z      dvf.Z       -dug.Z       -dvg.Z      L'N.z-L'N2.z          ]
  pOptJacobian[rlNumEquations+3][] = [ C'.-duf    C'.-dvf      0.0          0.0   C''C + C'C' - C''f - C'L'N ]

***********************************************************************/
SmStatus SmVariableRadiusFS::LoadJacobian
  (ULONG                  & rlNumEquations,       // i/o: number of equations to skip when setting
                                                  //      outputs rF and pOptJacobian
   ULONG                  & rlNumParameters,      // i/o: Gets incremented in an odd way to
                                                  //      to point to the input UVPnts.
   ULONG                    lRailIndex,           // in : Specifies rail index of first surface.
   SmSurface             *& rpSurface1,           // in : used only to map which crX[] uvPnt to surface1
   ULONG                  & rlSurf1Offset,        // i/o: index of rpSurface1 uvPnt in crX - ignored if rpSurface1 == NULL
   SmSurface             *& rpSurface2,           // in : used only to map which crX[] uvPnt to surface2
   ULONG                  & rlSurf2Offset,        // i/o: index of rpSurface2 uvPnt in crX - ignored if rpSurface2 == NULL
   const SmTArray<double> & crX,                  // in : current solution, [u1 v1 u2 v2 radius]
   SmTArray<double>       & rF,                   // out: sets rF[rlNumEquations] to rF[rlNumEquations+3] values
   SmMatrix               * pOptJacobian,         // out: sets pOptJacobian[rlNumEquations][] to pOptJacobian[rlNumEquations+3][] rows
   SmBoolean              & rbFoundAnswer)        // out: TRUE = all rF values less than scaled m_dConversionTol
{
  // init outputs
  rbFoundAnswer = FALSE;

  // get lSurf1Off index into crX array for lSurf1 parameters
  ULONG lSurf1Off = rlNumParameters;
  rlNumParameters = rlNumParameters + 2;
  if (lRailIndex == 0 && (rpSurface1 != NULL &&
      ((SmOffsetSurface*)rpSurface1)->GetBaseSurface() ==
      ((SmOffsetSurface*)GetSurface(0))->GetBaseSurface())) {
      lSurf1Off = rlSurf1Offset;
      rlNumParameters = rlNumParameters - 2;
  }
  if (lRailIndex == 1 && (rpSurface2 != NULL &&
      ((SmOffsetSurface*)rpSurface2)->GetBaseSurface() ==
      ((SmOffsetSurface*)GetSurface(0))->GetBaseSurface())) {
      lSurf1Off = rlSurf2Offset;
      rlNumParameters = rlNumParameters - 2;
  }

  // get lSurf2Off index into crX array for lSurf2 parameters
  ULONG lSurf2Off = rlNumParameters;
  rlNumParameters = rlNumParameters + 2;
  if (lRailIndex == 1 && (rpSurface1 != NULL &&
      ((SmOffsetSurface*)rpSurface1)->GetBaseSurface() ==
      ((SmOffsetSurface*)GetSurface(1))->GetBaseSurface())) {
      lSurf2Off = rlSurf1Offset;
      rlNumParameters = rlNumParameters - 2;
  }
  if (lRailIndex == 0 && (rpSurface2 != NULL &&
      ((SmOffsetSurface*)rpSurface2)->GetBaseSurface() ==
      ((SmOffsetSurface*)GetSurface(1))->GetBaseSurface())) {
      lSurf2Off = rlSurf2Offset;
      rlNumParameters = rlNumParameters - 2;
  }

  // get variableRadiusLaw Parameter
  double dT = crX[rlNumParameters];
  rlNumParameters ++;

  // evaluate radius its derivatives values
  double dValues[3];
  SER(m_rLaw->Evaluate(dT,m_vLawInterval,m_bLawOrientation,dValues));
  double dOffsetDistance = dValues[0];
  double dLp = dValues[1];

  // get OffsetSurface1 F point values for given crX values
  SmPoint3d sF;
  SmVector3d sDUF, sDVF;
  SmPoint2d sUV(crX[lSurf1Off],crX[lSurf1Off+1]);
  GetSurface(0)->SetOffsetDistance(dOffsetDistance*m_dOrientations[0]);
  SER(GetSurface(0)->Evaluate1stDerivatives(sUV,TRUE,TRUE,sF,sDUF,sDVF));
  SmVector3d sN = sDUF * sDVF;
  SER(sN.Unitize());
  sN = sN * m_dOrientations[0];

#ifdef SM_DEBUG_CODE
SmBoolean bDebugMe = FALSE;
  if (bDebugMe) {
      smgfx_Erase();
      sm_GraphicsLoop();
      smgfx_SetColor(0,1,0);
      GetSurface(0)->DrawAt(sUV,1);
      sm_GraphicsLoop();
  }
#endif

  // get OffsetSurface2 G point values for given crX values
  SmPoint3d sG;
  SmVector3d sDUG, sDVG;
  SmPoint2d sUV2(crX[lSurf2Off],crX[lSurf2Off+1]);
  GetSurface(1)->SetOffsetDistance(dOffsetDistance*m_dOrientations[1]);
  SER(GetSurface(1)->Evaluate1stDerivatives(sUV2,TRUE,TRUE,sG,sDUG,sDVG));
  SmVector3d sN2 = sDUG * sDVG;
  SER(sN2.Unitize());
  sN2 = sN2 * m_dOrientations[1];

#ifdef SM_DEBUG_CODE
  if (bDebugMe) {
      sm_GraphicsLoop();
      smgfx_SetColor(1,0,0);
      GetSurface(1)->DrawAt(sUV2,1);
      sm_GraphicsLoop();
  }
#endif

  // check surface tangent magnitudes
  double sDUFLenSq = sDUF.LengthSquared();
  double sDVFLenSq = sDVF.LengthSquared();

  double sDUGLenSq = sDUG.LengthSquared();
  double sDVGLenSq = sDVG.LengthSquared();

  // failure: zero length tangents
  if (   sDUFLenSq < SM_EFF_ZERO_SQ
      || sDVFLenSq < SM_EFF_ZERO_SQ
      || sDUGLenSq < SM_EFF_ZERO_SQ
      || sDVGLenSq < SM_EFF_ZERO_SQ)
    { return SM_ERR; }

  // set sDiff1 = vector between offsetSurface points
  SmVector3d sDiff1 = sF - sG;

#ifdef SM_DEBUG_CODE
  if (bDebugMe) {
      sm_GraphicsLoop();
      smgfx_SetColor(1,0,0);
      sDiff1.Draw(&sG);
      sm_GraphicsLoop();
  }
#endif

  ULONG lLastColumn = rlNumParameters - 1;

  // set sPV = BaseCurve[pos, tang, 2nd derivative] for dT value from crX
  SmVector3d sPV[3];
  SER(EvaluateBaseCurve(dT,2,sPV));

#ifdef SM_DEBUG_CODE
  if (bDebugMe) {
      sm_GraphicsLoop();
      smgfx_SetColor(0,0,1);
      sPV[0].Draw();
      smgfx_SetColor(0,1,1);
      sPV[1].Draw(&sPV[0]);
      sm_GraphicsLoop();
  }

SmBoolean bDebugMe2 = FALSE;
  if (bDebugMe2) {
//        smgfx_Erase();
      sm_GraphicsLoop();
      smgfx_SetColor(1,0,0);
      sPV[0].Draw();
      m_pOriginal->DrawWDeriv(m_pOriginal->GetNaturalInterval(),0);
      sm_GraphicsLoop();
      if (m_pExtended) m_pExtended->DrawWDeriv(m_pExtended->GetNaturalInterval(),0);
      sm_GraphicsLoop();
      smgfx_SetColor(0,0,1);
      GetSurface(0)->SetOffsetDistance(0.0);
      GetSurface(0)->DrawAt(sUV,0);
      sm_GraphicsLoop();
      smgfx_SetColor(1,0,1);
      GetSurface(1)->SetOffsetDistance(0.0);
      GetSurface(1)->DrawAt(sUV2,0);
      sm_GraphicsLoop();
  }
#endif


  // Intersection simply tries to zero the following 4 equations.
  // It is essentially intersecting the three surfaces by intersecting
  // tangent planes.
  // fun[i+0] = (f-g) . X     - X component of vector f-g
  // fun[i+1] = (f-g) . Y     - Y component of vector f-g
  // fun[i+2] = (f-g) . Z     - Z component of vector f-g
  // fun[i+3] = C'(t) * (C(t) - f)
  //
  // Where f is surface1 offset, g is surface2 offset and
  // where f2 is surface1 original and g2 is surface2 original
  //
  // It produces the following Jacobian with variables U, V, U2, V2, T
  //  | duf.X      dvf.X       -dug.X       -dvg.X      L'N.x-L'N2.x |
  //  | duf.Y      dvf.Y       -dug.Y       -dvg.Y      L'N.y-L'N2.y |
  //  | duf.Z      dvf.Z       -dug.Z       -dvg.Z      L'N.z-L'N2.z |
  //  | C'.-duf    C'.-dvf      0.0          0.0   C''C + C'C' - C''f - C'L'N |
  //
  //  N is the unitized normal of f2 and N2 is unitized normal of g2


  if (pOptJacobian)
    {
      ULONG i = rlNumEquations;

      (*pOptJacobian)[i][lSurf1Off]   =   sDUF.x;
      (*pOptJacobian)[i][lSurf1Off+1] =   sDVF.x;
      (*pOptJacobian)[i][lSurf2Off]   = - sDUG.x;
      (*pOptJacobian)[i][lSurf2Off+1] = - sDVG.x;
      (*pOptJacobian)[i][lLastColumn] =   dLp*(sN.x - sN2.x);


      i++; // Move to next equation;
      (*pOptJacobian)[i][lSurf1Off]   =   sDUF.y;
      (*pOptJacobian)[i][lSurf1Off+1] =   sDVF.y;
      (*pOptJacobian)[i][lSurf2Off]   = - sDUG.y;
      (*pOptJacobian)[i][lSurf2Off+1] = - sDVG.y;
      (*pOptJacobian)[i][lLastColumn] =   dLp*(sN.y - sN2.y);

      i++;
      (*pOptJacobian)[i][lSurf1Off]   =   sDUF.z;
      (*pOptJacobian)[i][lSurf1Off+1] =   sDVF.z;
      (*pOptJacobian)[i][lSurf2Off]   = - sDUG.z;
      (*pOptJacobian)[i][lSurf2Off+1] = - sDVG.z;
      (*pOptJacobian)[i][lLastColumn] =   dLp*(sN.z - sN2.z);

      i++;
      (*pOptJacobian)[i][lSurf1Off]   = - sPV[1].Dot(sDUF);
      (*pOptJacobian)[i][lSurf1Off+1] = - sPV[1].Dot(sDVF);
      (*pOptJacobian)[i][lSurf2Off]   =   0.0;
      (*pOptJacobian)[i][lSurf2Off+1] =   0.0;
      (*pOptJacobian)[i][lLastColumn] =   sPV[2].Dot(sPV[0]) + sPV[1].Dot(sPV[1])
                                        - sPV[2].Dot(sF)     - sPV[1].Dot(sN)*dLp;
    }

  // Compute function values
  rF[rlNumEquations++] = sDiff1.x;
  rF[rlNumEquations++] = sDiff1.y;
  rF[rlNumEquations++] = sDiff1.z;
  SmVector3d sDiff2    = sPV[0] - sF;
  rF[rlNumEquations++] = sPV[1].Dot(sDiff2);

  // See if we have converged
  double dScaledTol = GetConversionTol() * (1.0 + sF.GetMaxDimension());
  if (smos_Fabs(sDiff1.x) < dScaledTol &&
      smos_Fabs(sDiff1.y) < dScaledTol &&
      smos_Fabs(sDiff1.z) < dScaledTol &&
      smos_Fabs(rF[rlNumEquations-1]) < dScaledTol) {
      rbFoundAnswer = TRUE;
  }

#ifdef SM_DEBUG_CODE
  if (bDebugMe) {
      if (pOptJacobian) pOptJacobian->Dump();
      smos_WriteBuffer(_T(" F --- \n"));
      SM_DUMP_TARRAY( rF );
  }
#endif

#ifdef SM_DEBUG_CODE
SmBoolean bDebugMe3 = FALSE;
  if ( bDebugMe3 ) {
      SmOffsetSurface *pOffSrf1 = this->GetSurface(   lRailIndex );
      SmOffsetSurface *pOffSrf2 = this->GetSurface( 1-lRailIndex );
      const SmSurface *pBaseSrf1 = pOffSrf1->GetBaseSurface();
      const SmSurface *pBaseSrf2 = pOffSrf2->GetBaseSurface();
      SmEdge *pE = this->GetEdgeuse(0)->GetEdge();

      if ( FALSE ) {
          smgfx_Erase();

          SmObject *pOwnerObj = pBaseSrf1->GetOwner();
          SmFace *pFace = SM_CAST_PTR( SmFace, pOwnerObj );
          SmBrep *pBrep = ( pFace == NULL ) ? NULL : pFace->GetBrep();

          if ( pBrep ) {
              smgfx_SetLook( 1,2, 0,0,0 ); pBrep->Draw(TRUE); sm_GraphicsLoop();
              sm_GraphicsLoop();
          }

          smgfx_SetLook( 3,5, 1,0,1 ); if ( pE ) { pE->Draw(); } sm_GraphicsLoop();
          sm_GraphicsLoop();
      }

      SmPoint2d sUVDebug( crX[1], crX[2] );
      SmPoint3d sBasePnt1, sOffPnt1;
      SER(pBaseSrf1->EvaluatePoint( sUVDebug, sBasePnt1 ));
      SER(pOffSrf1 ->EvaluatePoint( sUVDebug, sOffPnt1  ));

      sm_GraphicsLoop();
      smgfx_SetLook(1,4, 0,0,1); sBasePnt1.DrawPointToPoint( sOffPnt1 ); sm_GraphicsLoop();
      sm_GraphicsLoop();

      sUVDebug.Set( crX[3], crX[4] );
      SmPoint3d sBasePnt2, sOffPnt2;
      SER(pBaseSrf2->EvaluatePoint( sUVDebug, sBasePnt2 ));
      SER(pOffSrf2 ->EvaluatePoint( sUVDebug, sOffPnt2  ));

      smgfx_SetLook(1,4, 0,1,1); sBasePnt2.DrawPointToPoint( sOffPnt2 ); sm_GraphicsLoop();
      sm_GraphicsLoop();

      if ( pE ) {
          SmCurve *pCrv = pE->GetCurve();
          smgfx_SetLook( 2,3, 1,0,0 ); pCrv->DrawAt( crX[5], 1 ); sm_GraphicsLoop();
          sm_GraphicsLoop();

          // SmPoint3d sEdgePt;
          // pCrv->EvaluatePoint( crX[5], sEdgePt );
          SmVector3d sPtVt[3];
          this->EvaluateBaseCurve( crX[5], 2, sPtVt );

          SmPoint3d sMid( ( sOffPnt1 + sOffPnt2 ) / 2 );

          smgfx_SetLook(1,4, 1,0,0); sPtVt[0].DrawPointToPoint( sMid ); sm_GraphicsLoop();
          sm_GraphicsLoop();
      }
  }
#endif

  return SM_SUCCESS;

} // end SmVariableRadiusFS::LoadJacobian

/*******************************************************************//**
PURPOSE: Get the number of variables the Jacobian for the constant
     radius fillet contains.

NOTES:
***********************************************************************/
ULONG SmVariableRadiusFS::GetJacobianSize
  ()
 const
{
    return 5;

} // end SmVariableRadiusFS::GetJacobianSize

/*******************************************************************//**
PURPOSE: Load the initial values for the solver.

NOTES:
***********************************************************************/
SmStatus SmVariableRadiusFS::LoadInitialValues
  (ULONG & rbDoSurf2Calcs,
   ULONG lRailIndex,
   const SmSurface & crSurface1,
   ULONG lSurface1Index,
   SmSurface *& pSurface2,
   ULONG & rlSurface2Index,
   SmExtentNd & rIntervals,
   SmTArray<SmBoolean> & rPeriodicities,
   SmTArray<double> & rGuessT)
{
    // Now find the UV guess on the other Rail's surface
    SmPoint3d sPnt;
    SmPoint2d sUV(rGuessT[lSurface1Index],rGuessT[lSurface1Index+1]);
    SER(crSurface1.EvaluatePoint(sUV,sPnt));
    // Use the extended curve if it's there.
    const SmCurve *p3DFilletEdgeCurve = ( m_pExtended ) ? m_pExtended : m_pOriginal;

    SmSolution sSData[16];
    SmSolutionArray sSolutions(16,sSData);

    // Use the curve's whole domain.
    SmExtent1d sIvl ( p3DFilletEdgeCurve->GetNaturalInterval() );

    SER( p3DFilletEdgeCurve->GlobalPointSolve( sIvl, SM_SO_MINIMIZE, sPnt,
        m_dThisApproxTol3d, NULL, NULL, SM_SR_ALL, sSolutions ));

    if ( sSolutions.GetSize() == 2 )  // [B391]
    {
        // Not handled yet: case of single closed fillet, where m_vLawInterval == sIvl.

        // Check whether one solution is closer to our interval.
        double dT0 = sSolutions[0].m_vStart[0];
        double dT1 = sSolutions[1].m_vStart[0];

        double dDist0 = m_vLawInterval.DistanceFrom( dT0 );
        double dDist1 = m_vLawInterval.DistanceFrom( dT1 );

        if ( dDist0 < dDist1 / 2.0 )
          { sSolutions.RemoveAt( 1 ); }
        else if ( dDist1 < dDist0 / 2.0 )
          { sSolutions.RemoveAt( 0 ); }
        else
          {
            // See whether one has a better solution value.
            double dVal0 = sSolutions[0].m_vStart.m_dSolutionValue;
            double dVal1 = sSolutions[1].m_vStart.m_dSolutionValue;
            if ( dVal0 < dVal1 * 0.9 )
              { sSolutions.RemoveAt( 1 ); }
            else if ( dVal1 < dVal0 * 0.9 )
              { sSolutions.RemoveAt( 0 ); }
            else
              {
                SE_MSG(SM_ERR, _T("Questionble Point Drop in LoadInitialValues(): Please investigate."));
                sSolutions.RemoveAt( 1 ); // This preserves previous behavior.
              }
          }
    } // end if two solutions


    if (sSolutions.GetSize() != 1)
      { SER(SM_ERR); }

    // Now determine edgeuses, surfaces and orientations for current iteration of the solver.
    double dFilletEdgeT = sSolutions[0].m_vStart[0];

    // If start point is beyond end of curve in one direction or the
    // other then use the parameter on the linearly extended curve.
    if ( ! sIvl.ContainsValue( dFilletEdgeT, -SM_EFF_ZERO ) ) // Neg tol: return False for on boundary
    {
        SmPoint3d sPV[2];
        SER(p3DFilletEdgeCurve->Evaluate( dFilletEdgeT, 1, TRUE, sPV ));
        double dParamExtension;
        SER(smgu_LineClosestPoint( sPV[0], sPV[1], sPnt, dParamExtension ));
        dFilletEdgeT = dFilletEdgeT + dParamExtension;
    }

    SmTsectPnt sTsectPnt;
    sTsectPnt.m_adUserDoubles[0] = dFilletEdgeT; // SetupOffsetValues() looks for this here.
    double dRadius;
    SER( SetupOffsetValues( sTsectPnt, &dRadius ));
    SER( crSurface1.EvaluatePoint( sUV, sPnt ));

    // Load surface 2 if neccesary
    if ( rbDoSurf2Calcs != 0 )
      {
        ULONG lOtherRailIndex = 1 - lRailIndex;
        SmExtent2d sDomain = GetSurface(1-lRailIndex)->GetBaseSurface()->GetNaturalUVDomain();
        SmVector2d sSize = GetDomainExtension( GetSurface(1-lRailIndex)->GetBaseSurface(), sDomain ); // SolidCraft change: and the trace's reach
        SmExtent2d sExtendedDomain(SmPoint2d(sDomain.GetMin().x-sSize.x,sDomain.GetMin().y-sSize.y),
            SmPoint2d(sDomain.GetMax().x+sSize.x, sDomain.GetMax().y+sSize.y) );

        pSurface2 = GetSurface(1-lRailIndex);

        SmPoint3d s2DPnt;

        // If we have edgeuses then use them to compute a point on the other surfaces
        // If not then drop the point to the other surface.
        SmEdgeuse *pOtherEU = m_pEdgeuses[lOtherRailIndex];
        if ( pOtherEU != NULL && pOtherEU->GetUVTrimCurve() == NULL )
          { pOtherEU = NULL; }

        // But: if we are rolling on a different face (e.g., doing a rollover),
        // then this is the wrong Edgeuse, and we don't have an appropriate one to use.
        if ( pOtherEU != NULL )
          {
            SmFace *pF = pOtherEU->GetFace();
            const SmSurface *pSrf2Base = ((SmOffsetSurface*)pSurface2)->GetBaseSurface();
            if ( pF != NULL && pF->GetSurface() != pSrf2Base )
              { pOtherEU = NULL; }
          }

        if ( pOtherEU != NULL)
          {
            SmBSplineCurve *pOtherUVCurve = m_pEdgeuses[lOtherRailIndex]->GetUVTrimCurve();
            SER(pOtherUVCurve->EvaluatePoint(sSolutions[0].m_vStart[0],s2DPnt));
            // Before we add it drop point to other surface to get a better guess
            SmBoolean bFound;
            SmSolution sSol;
            SmPoint2d sUVGuess(s2DPnt.x,s2DPnt.y);

            SER(pSurface2->LocalPointSolve(sExtendedDomain,SM_SO_MINIMIZE,sPnt,sUVGuess,bFound,sSol));
            if (bFound) {
                double dDist = sSol.m_vStart.m_dSolutionValue;
                if (dDist < 1.1*dRadius) {
                    s2DPnt.x = sSol.m_vStart[0];
                    s2DPnt.y = sSol.m_vStart[1];
                }
            }
        }
        else
        {
          // Do global solve
          SmSolution sSol[16];
          SmSolutionArray sSolutionArray( 16, sSol );
          SER( pSurface2->GlobalPointSolve( sExtendedDomain, SM_SO_MINIMIZE, sPnt, SM_EFF_ZERO_SQ, NULL, SM_SR_SINGLE, sSolutionArray ) );
          if(sSolutionArray.GetSize() == 1)
          {
            s2DPnt.x = sSolutionArray[0].m_vStart[0];
            s2DPnt.y = sSolutionArray[0].m_vStart[1];
          }
        }

        rIntervals[rGuessT.GetSize()] = SmExtent1d(sDomain.GetMin().x-sSize.x,sDomain.GetMax().x+sSize.x);
        rIntervals[rGuessT.GetSize()+1] = SmExtent1d(sDomain.GetMin().y-sSize.y,sDomain.GetMax().y+sSize.y);

        rPeriodicities.Add(FALSE);
        rPeriodicities.Add(FALSE);
        rlSurface2Index = rGuessT.GetSize();
        rGuessT.Add(s2DPnt.x);
        rGuessT.Add(s2DPnt.y);

#ifdef SM_DEBUG_CODE
SmBoolean bDebugMe = FALSE;
        if (bDebugMe) {
            sm_GraphicsLoop();
            sPnt.Draw();
            sm_GraphicsLoop();
            p3DFilletEdgeCurve->DrawAt(dFilletEdgeT,2);
            sm_GraphicsLoop();
            SmPoint2d sUVPt(s2DPnt.x,s2DPnt.y);
            pSurface2->DrawAt( sUVPt,2);
            sm_GraphicsLoop();
            GetSurface(lRailIndex)->DrawUV(4,8);
            sm_GraphicsLoop();
            smgfx_SetColor(0,0,1);
            GetSurface(1-lRailIndex)->DrawUV(9,18);
            sm_GraphicsLoop();
        }
#endif

    }

    // Now load the edge parameter data.
    // Make Edge extension factor go out further than surface extension (5x)
    double dEdgeExtension = 5.0 * m_dSurfaceExtensionFactor * (sIvl.GetMax() - sIvl.GetMin());
    rIntervals[rlSurface2Index+2] = SmExtent1d(sIvl.GetMin() - dEdgeExtension,
                                                 sIvl.GetMax() + dEdgeExtension);

    rPeriodicities.Add(FALSE);

    rGuessT.Add(dFilletEdgeT);


    return SM_SUCCESS;

} // end SmVariableRadiusFS::LoadInitialValues

/*******************************************************************//**
PURPOSE:
NOTES:
***********************************************************************/
SmBoolean SmVariableRadiusFS::IsKindOf( SM_TYPE t ) const
{
  return ((SmVariableRadiusFS_TYPE == t) ? TRUE : SmCurveBasedFS::IsKindOf( (t) ));
}

/*******************************************************************//**
PURPOSE: Dump a SmVariableRadiusFS.

NOTES: Unimplemented.
***********************************************************************/
void SmVariableRadiusFS::Dump() const {}


/*******************************************************************//**
PURPOSE: Evaluate either the original base curve or an extension of it.

NOTES: Typically the base curve is either an edge or an extension
    of the edge.  In the future it might be worth while to reapproximate
    the base curve to get smoother fillets.
***********************************************************************/
SmStatus SmCurveBasedFS::EvaluateBaseCurve
  (double dParameter,                  // in : target parameter
   ULONG lNumDerivatives,              // in : number of derivatives to output
   SmVector3d aPointAndDerivatives[])  // out: array of [pos, 1st deriv, 2nd deriv, ...]
                                       //      sized:[1+lNumDerivatives]
{
  // pick curve to use
  const SmCurve   *pCurve =  (m_pExtended)
                            ? m_pExtended
                            : m_pOriginal;

  // get curve interval
  SmExtent1d sIvl = pCurve->GetNaturalInterval();

  // when parameter is not within interval
  if (!sIvl.ContainsValue(dParameter))
    {
      // pick parameter space tolerance
      double dScaledTol = SM_EFF_ZERO * (1.0 + smos_Fabs(sIvl.GetMin()) + smos_Fabs(sIvl.GetMax()));

      // snap to begin point
      if (   dParameter < sIvl.GetMin()
          && dParameter + dScaledTol > sIvl.GetMin())
        {
          dParameter = sIvl.GetMin();
        }

      // snap to end point
      if (   dParameter > sIvl.GetMax()
          && dParameter - dScaledTol < sIvl.GetMax())
        {
          dParameter = sIvl.GetMax();
        }

      // when parameter is still outside interval
      if (!sIvl.ContainsValue(dParameter))
        {
          // Works better to extend the original: extending the extension
          // results in progressively smaller extensions.  [B551]
          pCurve = m_pOriginal;

          // when curve is a BSpline
          if (pCurve->IsKindOf(SmBSplineCurve_TYPE))
            {
              // create and store an extended curve
              SmBSplineCurve *pBSC = SM_CAST_PTR(SmBSplineCurve,pCurve); NER(pBSC);
              SmBSplineCurve *pExtended = NULL ;
              SmBoolean bPreciseExtension = TRUE;  // Otherwise it adds another length of sIvl.
              // But extend farther than requested, so that we won't have to do this too often.
              // This doubles the extension distance:
              double dExtendParam = ( dParameter<sIvl.GetMin() )
                          ? dParameter - ( sIvl.GetMin() - dParameter )
                          : dParameter + ( dParameter - sIvl.GetMax() );
              // But never do a really tiny extension: it might get snapped to the end,
              // and therefore not extended at all.  [Fillet 308]
              double dLen = sIvl.GetLength();
              if ( dExtendParam < sIvl.GetMin() )
              {
                  if ( dExtendParam > sIvl.GetMin() - dLen / 10 )
                    {  dExtendParam = sIvl.GetMin() - dLen / 10; }
              }
              else
              {
                  if ( dExtendParam < sIvl.GetMax() + dLen / 10 )
                    {  dExtendParam = sIvl.GetMax() + dLen / 10; }
              }
              SER(pBSC->CreateExtendedCurve( *pCurve->GetContext(), dExtendParam, SM_CT_G1,
                                              pExtended, bPreciseExtension ));
              if (m_pExtended) { delete m_pExtended; m_pExtended = NULL ; }
              m_pExtended = pExtended;
              pCurve      = pExtended;
            }
          else // Unable to extend non BSpline curves
            { SER(SM_ERR);
            }
        } // end parameter is still not within interval check after snapping endPoints
    } // end parameter is not within interval check

  // evaluate the curve
  SER(pCurve->Evaluate(dParameter,lNumDerivatives,TRUE,aPointAndDerivatives));

  // all done
  return SM_SUCCESS;

} // end SmCurveBasedFS::EvaluateBaseCurve


/*******************************************************************//**
PURPOSE: Destructor for Curve Based fillet solver.  It deletes any
     extended curves left.

NOTES: It assumes that it does not own the original curve.
***********************************************************************/
SmCurveBasedFS::~SmCurveBasedFS
  ()
{
    if (m_pExtended) { delete m_pExtended; m_pExtended = NULL ; }

} // end SmCurveBasedFS::~SmCurveBasedFS destructor

/*******************************************************************//**
PURPOSE: Utility methods.

NOTES:
***********************************************************************/
SmBoolean SmCurveBasedFS::IsKindOf( SM_TYPE t ) const
{
  return ((SmCurveBasedFS_TYPE == t) ? TRUE : SmFilletSolver::IsKindOf( (t) ));
}
/*******************************************************************//**
PURPOSE: Dump a SmCurveBasedFS.

NOTES: Unimplemented.
***********************************************************************/
void SmCurveBasedFS::Dump() const {}



/*******************************************************************//**
PURPOSE: Constructor for the surface surface fillet solver.

NOTES:
***********************************************************************/
SmSurfaceSurfaceFS::SmSurfaceSurfaceFS
  (const SmContext & crContext,
   double dThisApproxTol3d,
   double dAngleTolerance,
   double dTangencyTolerance,
   double dOffsetRadiusSurface1,
   double dOffsetRadiusSurface2,
   const SmSurface & crSurface1,
   const SmSurface & crSurface2,
   SmBoolean bOrientationSurface1,
   SmBoolean bOrientationSurface2,
   const SmPoint2d & crUVGuessSurface1,
   const SmPoint2d & crUVGuessSurface2)
  : SmConstantRadiusFS(crContext,
                       dThisApproxTol3d,
                       dAngleTolerance,
                       dTangencyTolerance,
                       dOffsetRadiusSurface1,
                       dOffsetRadiusSurface2,
                       crSurface1, 
                       crSurface2,
                       bOrientationSurface1, 
                       bOrientationSurface2),
    m_eTrimType(SM_BT_MINIMAL)

{
  m_vUVGuesses[0] = crUVGuessSurface1;
  m_vUVGuesses[1] = crUVGuessSurface2;


} // end SmSurfaceSurfaceFS::SmSurfaceSurfaceFS constructor

/*******************************************************************//**
PURPOSE: Destructor for the surface surface fillet solver.

NOTES: Offset Surfaces are deleted by SmConstantRadiusFS destructor.
***********************************************************************/
SmSurfaceSurfaceFS::~SmSurfaceSurfaceFS
  ()
{

} // end SmSurfaceSurfaceFS::~SmSurfaceSurfaceFS

// Not referenced anywhere
/*******************************************************************//**
PURPOSE: Make a rail edge for the surface surface fillet solver.

NOTES:
**********************************************************************
static SmStatus MakeRail(SmFilletBrep * pPseudoBrep,
                         SmFilletEdge *& rpNewEdge)
{
    SmFilletEdge * pRail = new (pPseudoBrep) SmFilletEdge(pPseudoBrep);
    NER(pRail); SmObjDelete sClean1(pRail);
    pRail->SetFilletEdgeType(SM_FE_RAIL);
    SmFilletEdgeuse *pNewEU1 = new (pPseudoBrep) SmFilletEdgeuse(pPseudoBrep);
    NER(pNewEU1); SmObjDelete sClean2(pNewEU1);
    SmFilletEdgeuse *pNewEU2 = new (pPseudoBrep) SmFilletEdgeuse(pPseudoBrep);
    NER(pNewEU2); SmObjDelete sClean3(pNewEU2);
    SmFilletVertexuse *pNewVU1 = new (pPseudoBrep) SmFilletVertexuse();
    NER(pNewVU1); SmObjDelete sClean4(pNewVU1);
    pNewVU1->SetProperty(pNewEU1);
    SmFilletVertexuse *pNewVU2 = new (pPseudoBrep) SmFilletVertexuse();
    NER(pNewVU2); SmObjDelete sClean5(pNewVU2);
    SmFilletVertex * pStartVert = new (pPseudoBrep) SmFilletVertex();
    NER(pStartVert); SmObjDelete sClean6(pStartVert);
    SmFilletVertex * pEndVert = new (pPseudoBrep) SmFilletVertex();
    NER(pEndVert); SmObjDelete sClean7(pEndVert);
     // Put VU1 into start vertex list and VU2 into end vertex
    SER(pStartVert->PostInsert(pNewVU1));
    SER(pEndVert->PostInsert(pNewVU2));
    pNewVU2->SetProperty(pNewEU2);
    SER(pRail->PreInsert(pNewEU1));
    SER(pRail->PostInsert(pNewEU2));
    pNewEU1->SetFilletVertexuse(pNewVU1);
    pNewEU2->SetFilletVertexuse(pNewVU2);
    // Clear cleanup objects
    sClean1.Clear(); sClean2.Clear(); sClean3.Clear();
    sClean4.Clear(); sClean5.Clear(); sClean6.Clear();
    sClean7.Clear();

    rpNewEdge = pRail;

    return SM_SUCCESS;
}
*/

/*******************************************************************//**
PURPOSE: Compute the geometry of a surface surface fillet solver by
    rolling a ball along the surfaces (intersecting the offset surfaces)
    and then trimming the results by the trim boundary of the face.

NOTES:
***********************************************************************/
SmStatus SmSurfaceSurfaceFS::CalcFilletGeom()
{
  SmSurface * pSurf1 = GetSurface(0);
  SmSurface * pSurf2 = GetSurface(1);
  SmExtent2d sDomain1 = pSurf1->GetNaturalUVDomain();
  SmExtent2d sDomain2 = pSurf2->GetNaturalUVDomain();

  // Need to solve for start point using guess point
  SmTsectPnt sCurrPnt;
  sCurrPnt.UVPos(0)   = m_vUVGuesses[0];
  sCurrPnt.UVDeriv(0) = SmPoint2d(0.0, 0.0);
  sCurrPnt.UVPos(1)   = m_vUVGuesses[1];
  sCurrPnt.UVDeriv(1) = SmPoint2d(0.0, 0.0);

  // Compute a plane that should intersect fillet edge and use it to solve
  SmPoint3d sPnt1, sPnt2;
  SER(pSurf1->EvaluatePoint(m_vUVGuesses[0],sPnt1));
  SER(pSurf2->EvaluatePoint(m_vUVGuesses[1],sPnt2));
  sCurrPnt.CrvPos() = (sPnt1+sPnt2)/2.0;

  SmVector3d sNorm1, sNorm2;
  SER(pSurf1->EvaluateNormal(m_vUVGuesses[0],TRUE,TRUE,sNorm1));
  SER(pSurf2->EvaluateNormal(m_vUVGuesses[1],TRUE,TRUE,sNorm2));
  sCurrPnt.CrvDeriv() = sNorm1 * sNorm2;

  // Bad guess points.  They must be points where the surfaces tangent planes are not locally parallel
  if (sCurrPnt.CrvDeriv().LengthSquared() < SM_EFF_ZERO_SQ) 
    { return SM_ERR; }

#ifdef SM_DEBUG_CODE
SmBoolean bDebugMe = FALSE;
  if (bDebugMe) 
    {
      smgfx_SetLook(2,4, 0,0,1); sNorm1.Draw(&sPnt1); sm_GraphicsLoop();
      smgfx_SetLook(2,4, 1,0,0); sNorm2.Draw(&sPnt2); sm_GraphicsLoop();
      smgfx_SetLook(1,2, 0,1,1); pSurf1->DrawUV(4,4); sm_GraphicsLoop();
      smgfx_SetLook(1,2, 1,0,1); pSurf2->DrawUV(6,6); sm_GraphicsLoop();
      sm_GraphicsLoop();
    }
#endif // SM_DEBUG_CODE

  double     dStep = 0.0001;
  SmBoolean  bFoundGoodPoint, bClipped, bHitBoundary;
  SmTsectPnt sNextPnt;
  SER(StepSolve(sCurrPnt,
                dStep,
                sDomain1,
                sDomain2,
                bFoundGoodPoint, 
                bClipped, 
                bHitBoundary, 
                sNextPnt));

  // If we did not find a good point, try again, then quit.
  if (!bFoundGoodPoint) 
    {
      SmBoolean  bFoundAnswer;
      SmSolution sSolution;
      SmPoint3d  sPnt3, sPnt4;

      for (ULONG ii=0; ii<3; ii++) 
        {
          SER(pSurf2->EvaluatePoint(sCurrPnt.UVPos(1),sPnt4));
          SER(pSurf1->LocalPointSolve(sDomain1,
                                      SM_SO_MINIMIZE,
                                      sPnt4,
                                      sCurrPnt.UVPos(0), 
                                      bFoundAnswer, 
                                      sSolution));
          if (bFoundAnswer) 
            {
              sCurrPnt.UVPos(0) = SmPoint2d(sSolution.m_vStart[0],sSolution.m_vStart[1]);
            }

          SER(pSurf1->EvaluatePoint(sCurrPnt.UVPos(0),sPnt3));
          SER(pSurf2->LocalPointSolve(sDomain2,
                                      SM_SO_MINIMIZE,
                                      sPnt3,
                                      sCurrPnt.UVPos(1), 
                                      bFoundAnswer, 
                                      sSolution));
          if (bFoundAnswer) 
            {
              sCurrPnt.UVPos(1) = SmPoint2d(sSolution.m_vStart[0],sSolution.m_vStart[1]);
            }
        }
      SER(StepSolve(sCurrPnt,
                    dStep,
                    sDomain1,
                    sDomain2,
                    bFoundGoodPoint, 
                    bClipped, 
                    bHitBoundary, 
                    sNextPnt));
      if (!bFoundGoodPoint) 
        {
          return SM_SUCCESS;
        }
    }

  SmTsectPnt * sP1Data[4];
  SmTArray<SmTsectPnt*> sTsectPnts(4,sP1Data);

  // Add start point
  sTsectPnts.Add(&sNextPnt);

  SetupOffsetExtension(pSurf1,sDomain1);
  SetupOffsetExtension(pSurf2,sDomain2);

  if (m_vFilletGeoms.GetSize() == 0) 
    { return SM_SUCCESS; }

  SM_ASSERT(m_vFilletGeoms.GetSize() == 1);
  SmFilletGeom * pFilletGeom = m_vFilletGeoms[0];
  SmFilletBrep * pBrep       = m_pExecutive->GetPseudoBrep();

  // Create Rails
  SmFace * pFace[2];
  pFace[0] = (SmFace*)((SmOffsetSurface*)pSurf1)->GetBaseSurface()->GetFace(); NER(pFace[0]);
  pFace[1] = (SmFace*)((SmOffsetSurface*)pSurf2)->GetBaseSurface()->GetFace(); NER(pFace[1]);

  for (ULONG ii=0; ii<2; ii++) 
    {
      SmFilletEdge * pRail = NULL;
      SER(pFilletGeom->MakeRailEdge(pBrep, ii, pRail));
      pRail->SetOriginalFace(pFace[ii]);
    }

  // Setup FilletIntersector
  SmFilletIntersector sFI(*pSurf1,
                           sDomain1,
                          *pSurf2,
                           sDomain2,
                          *this);

  SmBSplineCurve * p3DCurve = NULL ;
  SmBSplineCurve * pSurfaceUVCurve[2];
  SmTsectCurveType eCurveType;
  double           dDeviation;
  double           dAppTol = GetThisApproxTol3d();
  double           dAngTol = GetThisAngTolRad();

  // Trace the spine and rail curves and create fillet surafce.
  // This puts the results into the FilletGeom in 'this' FilletSolver
  // (and not into the output arguemnts).
  SER(sFI.DoPointIntersection(m_crContext,
                              sTsectPnts, 
                              TRUE, TRUE, NULL, NULL, 
                              SM_CAST_APPROXTOL3D_PTR(&dAppTol), 
                              &dAngTol, 
                              p3DCurve,
                              pSurfaceUVCurve[0], 
                              pSurfaceUVCurve[1], 
                              eCurveType, 
                              dDeviation));

  // Here we need to go backward through the process and derive
  // rail intersection points.
  SmSurface * pFilletSurf = pFilletGeom->GetFilletSurface(); NER(pFilletSurf);
  SmExtent2d  sUVDomain   = pFilletSurf->GetNaturalUVDomain();
  SmPoint2d   sUVMax = sUVDomain.GetMax();
  SmPoint2d   sUVMin = sUVDomain.GetMin();
  double      dV[2];// 'V' parameter of each rail(iso-curve)
  dV[0] = sUVDomain.GetMin().y;
  dV[1] = sUVDomain.GetMax().y;

  // Classify the rail curves relative to the faces.
  for (ULONG ii=0; ii<2; ii++) 
    {
      SmFilletEdge    * pRail   = pFilletGeom->GetRail(ii);
      SmFilletEdgeuse * pPrimEU = (SmFilletEdgeuse*)pRail->GetPrimaryEdgeuse();
      SmFilletEdgeuse * pMateEU = (SmFilletEdgeuse*)pPrimEU->GetMate();

      // Create a new CurveClass object and install it into pRail.
      SmCurve               * pCurve      = pRail->GetCurve();
      SmCurveClassification * pCurveClass = new(m_crContext) SmCurveClassification(pCurve, 
                                                                                   pCurve->GetNaturalInterval(),
                                                                                   pMateEU->GetUVTrimCurve(),
                                                                                   GetThisApproxTol3d());
      pRail->SetCurveClass(pCurveClass);

      // Classify rail curve relative to its owning face.
      SER(pRail->GetOriginalFace()->CurveOnClassify(TRUE,         // in : TRUE =Do expensive PtClassification for Curves not XSecting any Face->Bndrys
                                                    *pCurveClass, // i/o: contains curve to classify, accumulates intervals as they are found.
                                                    TRUE));       // in : TRUE=increase Edge tols for new verts found to be slightly too far from their surfaces

#ifdef SM_DEBUG_CODE
SmBoolean bDebugMe2 = FALSE;
      if (bDebugMe2) 
        {
          pCurveClass->Dump();
          pCurveClass->Draw();
          sm_GraphicsLoop();
        }
#endif // SM_DEBUG_CODE

      // Create UV-curve on fillet for this rail: isoparameter.
      SmPoint3d sStartUV( sUVMin.x, dV[ii], 0 );
      SmPoint3d sEndUV  ( sUVMax.x, dV[ii], 0 );

      SmBSplineCurve *pUVCurve = NULL ;
      SER(SmBSplineCurve::CreateLineSegment(m_crContext,
                                            2,
                                            sStartUV,
                                            sEndUV,
                                            pUVCurve));
      NER(pUVCurve);

      SmExtent1d sIvl(sUVMin.x, sUVMax.x);
      SER(pUVCurve->EditParameterization(sIvl));
      pPrimEU->SetUVCurve(pUVCurve);
    }

  // Use the classifications to rebuild the surface and
  // set up the trimming information.
  SER( SplitAndTrimFilletSurface( 0, m_eTrimType ));

  // all done
  return SM_SUCCESS;

} // end SmSurfaceSurfaceFS::CalcFilletGeom

/*******************************************************************//**
PURPOSE: Determine the amount an offset surface can be extended for
         evaluation purposes during the fillet process.

NOTES: Same as base class version except zero extension for periodic surfaces. [B566]
***********************************************************************/
void SmSurfaceSurfaceFS::SetupOffsetExtension
  (const SmSurface * pSurf,    // in : target surface
   SmExtent2d      & rDomain)  // i/o: target surface domain of interest
{
  // get length of extensions
  SmVector2d sSize = m_dSurfaceExtensionFactor * rDomain.GetSize();

  // Do not extend if periodic.
  double dPeriodicExtensionFactor = 0.0;

  if (pSurf->IsPeriodic(rDomain,SM_SP_U))
    {
      //if (!m_bExtendBefore && !m_bExtendAfter)
      sSize.x *= dPeriodicExtensionFactor;
    }
  if (pSurf->IsPeriodic(rDomain,SM_SP_V))
    {
      //if (!m_bExtendBefore && !m_bExtendAfter)
      sSize.y *= dPeriodicExtensionFactor;
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

} // end void SmSurfaceSurfaceFS::SetupOffsetExtension

/*******************************************************************//**
PURPOSE: Utility methods.

NOTES:
***********************************************************************/
SmBoolean SmSurfaceSurfaceFS::IsKindOf( SM_TYPE t ) const
{
  return ((SmSurfaceSurfaceFS_TYPE == t) ? TRUE : SmConstantRadiusFS::IsKindOf( (t) ));
}
/*******************************************************************//**
PURPOSE: Dump a SmSurfaceSurfaceFS.

NOTES: Unimplemented.
***********************************************************************/
void SmSurfaceSurfaceFS::Dump() const {}
