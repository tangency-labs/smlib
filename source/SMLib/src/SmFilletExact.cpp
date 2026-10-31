// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

/**********************************************************************//**
* FILE NAME --- SmFilletExact.cpp
* PURPOSE: SolidCraft change. Which constant rounds are exact, and their exact
*   surfaces.
*
* A constant round of radius r is the envelope of a ball of radius r whose
* centre runs on the spine, where the faces' offsets by r meet. The round is
* rational exactly where the spine is: a plane section of a natural quadric
* (a conic, or a line or a circle), one of the two plane conics two cones meet
* in when their offsets touch a common ball the same way round, or the quartic
* with a node two of them (or one and a sphere) meet in when they touch at one
* point the same way round. Each of these holds for every r when it holds for
* one: it is a matter of the faces' placement alone. Every other pair of
* natural quadrics, and every other surface, meets in a spine that is not
* rational, and its round is traced.
**********************************************************************/

#include "StdAfx.h"

#include <SmFilletExact.h>

#include <SmBrep.h>
#include <SmBSplineCurve.h>
#include <SmBSplineSurface.h>
#include <SmEdge.h>
#include <SmEdgeuse.h>
#include <SmSolutionArray.h>
#include <SmFilletGeom.h>
#include <SmFilletSolver.h>
#include <SmCone.h>
#include <SmCylinder.h>
#include <SmOffsetSurface.h>
#include <SmPlane.h>
#include <SmSphere.h>
#include <SmSurface.h>

#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

/*******************************************************************//**
PURPOSE: The component of crV square to the unit crDir.
***********************************************************************/
static SmVector3d sm_Across(const SmVector3d & crV, const SmVector3d & crDir)
{
  return crV - crDir * crV.Dot( crDir );
}

/*******************************************************************//**
PURPOSE: The natural quadric of an analytic surface, oriented by its normal
   at crUV, or SM_BQ_NONE. pSurface is no SmOffsetSurface.
***********************************************************************/
static SmStatus sm_QuadricOfSurface(const SmSurface * pSurface, const SmPoint2d & crUV, double dTol, SmBlendQuadric & rQuadric)
{
  rQuadric = SmBlendQuadric();

  SmPoint3d  sP;
  SmVector3d sN;
  SER( pSurface->EvaluatePoint( crUV, sP ) );
  SER( pSurface->EvaluateNormal( crUV, TRUE, TRUE, sN ) );

  const SmBSplineSurface * pBSS = SM_CAST_PTR( SmBSplineSurface, pSurface );
  if ( pBSS == NULL )
    { return SM_SUCCESS; }

  // the analytic kinds, or a B-spline that is one (a plane, a cone or cylinder, a sphere)
  const SmCone   * pCone   = SM_CAST_PTR( SmCone, pSurface );
  const SmSphere * pSphere = SM_CAST_PTR( SmSphere, pSurface );
  SmCone   * pNurbCone   = NULL;
  SmSphere * pNurbSphere = NULL;
  SmObjDelete sDelCone, sDelSphere;
  SmPoint3d  sPlanePoint;
  SmVector3d sPlaneNormal;
  SmBoolean  bPlane = pSurface->GetType() == SmPlane_TYPE;
  if ( !bPlane && pCone == NULL && pSphere == NULL )
    {
      if ( pBSS->IsPlanar( dTol, &sPlanePoint, &sPlaneNormal ) )
        { bPlane = TRUE; }
      else if ( SmCone::IsNurbSurfaceCone( *pBSS->GetContext(), pBSS, pNurbCone ) && pNurbCone != NULL )
        { sDelCone.SetObj( pNurbCone ); pCone = pNurbCone; }
      else if ( SmSphere::IsNurbSurfaceSphere( *pBSS->GetContext(), pBSS, pNurbSphere ) && pNurbSphere != NULL )
        { sDelSphere.SetObj( pNurbSphere ); pSphere = pNurbSphere; }
    }

  if ( bPlane )
    {
      rQuadric.m_eKind  = SmBlendQuadric::SM_BQ_PLANE;
      rQuadric.m_sPoint = sP;
      rQuadric.m_sDir   = sN;
    }
  else if ( pCone != NULL )
    {
      // Rad(v) = R0 + tan(semi-angle) v about the axis through the origin along z
      SmVector3d sZ( pCone->GetZAxis() );
      SER( sZ.Unitize() );
      const double dR0   = pCone->GetRadius( 0.0 );
      const double dTanB = pCone->GetRadius( 1.0 ) - dR0;
      if ( pCone->IsCylinder() )
        {
          rQuadric.m_eKind   = SmBlendQuadric::SM_BQ_CYLINDER;
          rQuadric.m_sPoint  = pCone->GetOrigin();
          rQuadric.m_sDir    = sZ;
          rQuadric.m_dRadius = smos_Fabs( dR0 );
        }
      else
        {
          rQuadric.m_eKind  = SmBlendQuadric::SM_BQ_CONE;
          rQuadric.m_sPoint = pCone->GetOrigin() + sZ * ( -dR0 / dTanB );
          rQuadric.m_sDir   = ( dTanB > 0.0 ) ? sZ : -sZ;
          rQuadric.m_dAngle = smos_ArcTangent( smos_Fabs( dTanB ) );
        }
      // the nappe the surface point is on
      if ( rQuadric.m_eKind == SmBlendQuadric::SM_BQ_CONE && ( sP - rQuadric.m_sPoint ).Dot( rQuadric.m_sDir ) < 0.0 )
        { rQuadric.m_sDir = -rQuadric.m_sDir; }
      SmVector3d sOut( rQuadric.Normal( sP ) );
      rQuadric.m_dSide = ( sOut.Dot( sN ) >= 0.0 ) ? 1.0 : -1.0;
    }
  else if ( pSphere != NULL )
    {
      rQuadric.m_eKind   = SmBlendQuadric::SM_BQ_SPHERE;
      rQuadric.m_sPoint  = pSphere->GetOrigin();
      rQuadric.m_dRadius = pSphere->GetRadius();
      rQuadric.m_dSide   = ( ( sP - rQuadric.m_sPoint ).Dot( sN ) >= 0.0 ) ? 1.0 : -1.0;
    }
  else
    { return SM_SUCCESS; }

  // what was read off the surface is the surface
  if ( rQuadric.Distance( sP ) > dTol || rQuadric.Normal( sP ).Dot( sN ) < 1.0 - SM_EFF_ZERO_SQRT )
    { SER_MSG( SM_ERR, _T("SmBlendQuadric::FromSurface(): the quadric read off a surface does not pass through it") ); }
  return SM_SUCCESS;
} // end static sm_QuadricOfSurface

/*******************************************************************//**
PURPOSE: The natural quadric pSurface lies on, oriented by its normal at crUV.
   An SmOffsetSurface is its base's quadric moved by its offset there.
***********************************************************************/
SmStatus SmBlendQuadric::FromSurface(const SmSurface * pSurface, const SmPoint2d & crUV, double dTol, SmBlendQuadric & rQuadric)
{
  rQuadric = SmBlendQuadric();
  NER( pSurface );
  const SmOffsetSurface * pOffset = SM_CAST_PTR( SmOffsetSurface, pSurface );
  if ( pOffset == NULL )
    { return sm_QuadricOfSurface( pSurface, crUV, dTol, rQuadric ); }

  const SmSurface * pBase = pOffset->GetBaseSurface(); NER( pBase );
  SmBlendQuadric sBase;
  SER( FromSurface( pBase, crUV, dTol, sBase ) );
  if ( sBase.m_eKind == SM_BQ_NONE )
    { return SM_SUCCESS; }

  // the offset there, along the base's normal
  SmPoint3d  sP, sQ;
  SmVector3d sN;
  SER( pBase->EvaluatePoint( crUV, sP ) );
  SER( pBase->EvaluateNormal( crUV, TRUE, TRUE, sN ) );
  SER( pOffset->EvaluatePoint( crUV, sQ ) );
  rQuadric = sBase.Offset( ( sQ - sP ).Dot( sN ), dTol );
  return SM_SUCCESS;
} // end SmBlendQuadric::FromSurface

/*******************************************************************//**
PURPOSE: This quadric moved dDist along its normal.
***********************************************************************/
SmBlendQuadric SmBlendQuadric::Offset(double dDist, double dTol) const
{
  SmBlendQuadric sO( *this );
  switch ( m_eKind )
    {
    case SM_BQ_PLANE:
      sO.m_sPoint = m_sPoint + m_sDir * dDist;
      break;
    case SM_BQ_SPHERE:
    case SM_BQ_CYLINDER:
      sO.m_dRadius = m_dRadius + m_dSide * dDist;
      if ( smos_Fabs( sO.m_dRadius ) <= dTol )
        {
          sO.m_eKind   = ( m_eKind == SM_BQ_SPHERE ) ? SM_BQ_POINT : SM_BQ_LINE;
          sO.m_dRadius = 0.0;
        }
      else if ( sO.m_dRadius < 0.0 )
        {
          // moved in past its centre or axis: the same surface the other way round
          sO.m_dRadius = -sO.m_dRadius;
          sO.m_dSide   = -m_dSide;
        }
      break;
    case SM_BQ_CONE:
      // the normal away from the axis is cos(a) radial - sin(a) d: the apex moves back along d
      sO.m_sPoint = m_sPoint - m_sDir * ( m_dSide * dDist / smos_Sine( m_dAngle ) );
      break;
    default:
      break;
    }
  return sO;
} // end SmBlendQuadric::Offset

/*******************************************************************//**
PURPOSE: The unit normal at crPoint, a point on the quadric, as the surface's.
***********************************************************************/
SmVector3d SmBlendQuadric::Normal(const SmPoint3d & crPoint) const
{
  SmVector3d sN( 0.0, 0.0, 0.0 );
  switch ( m_eKind )
    {
    case SM_BQ_PLANE:
      return m_sDir;
    case SM_BQ_SPHERE:
      sN = crPoint - m_sPoint;
      break;
    case SM_BQ_CYLINDER:
      sN = sm_Across( crPoint - m_sPoint, m_sDir );
      break;
    case SM_BQ_CONE:
      {
        SmVector3d sR( sm_Across( crPoint - m_sPoint, m_sDir ) );
        const double dLen = sR.Length();
        if ( dLen > 0.0 )
          { sN = sR * ( smos_Cosine( m_dAngle ) / dLen ) - m_sDir * smos_Sine( m_dAngle ); }
      }
      break;
    default:
      return sN;
    }
  const double dLen = sN.Length();
  return ( dLen > 0.0 ) ? sN * ( m_dSide / dLen ) : sN;
} // end SmBlendQuadric::Normal

/*******************************************************************//**
PURPOSE: The distance of crPoint from the quadric (from the cone's nappe, in
   the plane through the axis).
***********************************************************************/
double SmBlendQuadric::Distance(const SmPoint3d & crPoint) const
{
  const SmVector3d sW( crPoint - m_sPoint );
  switch ( m_eKind )
    {
    case SM_BQ_PLANE:
      return smos_Fabs( sW.Dot( m_sDir ) );
    case SM_BQ_SPHERE:
    case SM_BQ_POINT:
      return smos_Fabs( sW.Length() - m_dRadius );
    case SM_BQ_CYLINDER:
    case SM_BQ_LINE:
      return smos_Fabs( sm_Across( sW, m_sDir ).Length() - m_dRadius );
    case SM_BQ_CONE:
      return smos_Fabs( sm_Across( sW, m_sDir ).Length() * smos_Cosine( m_dAngle ) - sW.Dot( m_sDir ) * smos_Sine( m_dAngle ) );
    default:
      return SM_BIG_DOUBLE;
    }
} // end SmBlendQuadric::Distance

/*******************************************************************//**
PURPOSE: The implicit form x.A.x + 2 b.x + c = 0.
***********************************************************************/
void SmBlendQuadric::Implicit(double aA[3][3], double aB[3], double & rC) const
{
  const double aP[3] = { m_sPoint.x, m_sPoint.y, m_sPoint.z };
  const double aD[3] = { m_sDir.x, m_sDir.y, m_sDir.z };
  // |x - p|^2 kI - (d.(x - p))^2 kD - k0
  double dKI = 0.0, dKD = 0.0, dK0 = 0.0;
  switch ( m_eKind )
    {
    case SM_BQ_PLANE:
      for ( int i=0; i<3; i++ )
        {
          for ( int j=0; j<3; j++ ) { aA[i][j] = 0.0; }
          aB[i] = 0.5 * aD[i];
        }
      rC = -( aP[0]*aD[0] + aP[1]*aD[1] + aP[2]*aD[2] );
      return;
    case SM_BQ_SPHERE:
    case SM_BQ_POINT:
      dKI = 1.0; dK0 = m_dRadius * m_dRadius;
      break;
    case SM_BQ_CYLINDER:
    case SM_BQ_LINE:
      dKI = 1.0; dKD = 1.0; dK0 = m_dRadius * m_dRadius;
      break;
    case SM_BQ_CONE:
      dKI = smos_Cosine( m_dAngle ) * smos_Cosine( m_dAngle ); dKD = 1.0;
      break;
    default:
      break;
    }
  for ( int i=0; i<3; i++ )
    {
      for ( int j=0; j<3; j++ )
        { aA[i][j] = ( i == j ? dKI : 0.0 ) - dKD * aD[i] * aD[j]; }
    }
  // b = -A p
  for ( int i=0; i<3; i++ )
    { aB[i] = -( aA[i][0]*aP[0] + aA[i][1]*aP[1] + aA[i][2]*aP[2] ); }
  // c = p.A.p - k0
  double dPAP = 0.0;
  for ( int i=0; i<3; i++ )
    {
      for ( int j=0; j<3; j++ ) { dPAP += aP[i] * aA[i][j] * aP[j]; }
    }
  rC = dPAP - dK0;
} // end SmBlendQuadric::Implicit

/*******************************************************************//**
PURPOSE: The point where two lines meet within dTol (their nearest points'
   midpoint), FALSE when they are parallel or miss each other.
***********************************************************************/
static SmBoolean sm_LinesMeet(const SmPoint3d & crP0, const SmVector3d & crD0, const SmPoint3d & crP1, const SmVector3d & crD1,
                              double dTol, SmPoint3d & rMeet)
{
  const SmVector3d sW( crP0 - crP1 );
  const double dB   = crD0.Dot( crD1 );
  const double dDen = crD0.Dot( crD0 ) * crD1.Dot( crD1 ) - dB * dB;
  if ( dDen <= SM_EFF_ZERO_SQRT * crD0.Dot( crD0 ) * crD1.Dot( crD1 ) )
    { return FALSE; }
  const double dD = crD0.Dot( sW );
  const double dE = crD1.Dot( sW );
  const SmPoint3d sQ0 = crP0 + crD0 * ( ( dB * dE - crD1.Dot( crD1 ) * dD ) / dDen );
  const SmPoint3d sQ1 = crP1 + crD1 * ( ( crD0.Dot( crD0 ) * dE - dB * dD ) / dDen );
  if ( sQ0.DistanceBetween( sQ1 ) > dTol )
    { return FALSE; }
  rMeet = sQ0 + ( sQ1 - sQ0 ) * 0.5;
  return TRUE;
} // end static sm_LinesMeet

/*******************************************************************//**
PURPOSE: Whether two quadrics of revolution (or a plane, or a sphere's centre)
   share their axis within dTol at a distance dReach from it: then their
   offsets meet in circles about it. rPoint, rDir: the axis.
***********************************************************************/
static SmBoolean sm_Coaxial(const SmBlendQuadric & crA, const SmBlendQuadric & crB, double dTol, double dReach,
                            SmPoint3d & rPoint, SmVector3d & rDir)
{
  typedef SmBlendQuadric Q;
  const SmBoolean bPointA = crA.m_eKind == Q::SM_BQ_SPHERE || crA.m_eKind == Q::SM_BQ_POINT;
  const SmBoolean bPointB = crB.m_eKind == Q::SM_BQ_SPHERE || crB.m_eKind == Q::SM_BQ_POINT;
  const SmBoolean bAxisA  = crA.m_eKind == Q::SM_BQ_CYLINDER || crA.m_eKind == Q::SM_BQ_CONE || crA.m_eKind == Q::SM_BQ_LINE;
  const SmBoolean bAxisB  = crB.m_eKind == Q::SM_BQ_CYLINDER || crB.m_eKind == Q::SM_BQ_CONE || crB.m_eKind == Q::SM_BQ_LINE;
  const double    dAngTol = dTol / smos_Max( dReach, dTol );

  if ( bPointA && bPointB )
    {
      SmVector3d sD( crB.m_sPoint - crA.m_sPoint );
      if ( sD.Length() <= dTol ) { return FALSE; } // concentric: no circle
      if ( sD.Unitize() != SM_SUCCESS ) { return FALSE; }
      rPoint = crA.m_sPoint; rDir = sD;
      return TRUE;
    }
  if ( ( bPointA && bAxisB ) || ( bAxisA && bPointB ) )
    {
      const Q & crP = bPointA ? crA : crB;
      const Q & crX = bPointA ? crB : crA;
      if ( sm_Across( crP.m_sPoint - crX.m_sPoint, crX.m_sDir ).Length() > dTol ) { return FALSE; }
      rPoint = crX.m_sPoint; rDir = crX.m_sDir;
      return TRUE;
    }
  if ( bAxisA && bAxisB )
    {
      if ( ( crA.m_sDir * crB.m_sDir ).Length() > dAngTol ) { return FALSE; }
      if ( sm_Across( crB.m_sPoint - crA.m_sPoint, crA.m_sDir ).Length() > dTol ) { return FALSE; }
      rPoint = crA.m_sPoint; rDir = crA.m_sDir;
      return TRUE;
    }
  // a plane square to an axis, or through a sphere's centre's normal line
  const Q & crPl = ( crA.m_eKind == Q::SM_BQ_PLANE ) ? crA : crB;
  const Q & crX  = ( crA.m_eKind == Q::SM_BQ_PLANE ) ? crB : crA;
  if ( crPl.m_eKind != Q::SM_BQ_PLANE || crX.m_eKind == Q::SM_BQ_PLANE )
    { return FALSE; }
  if ( crX.m_eKind == Q::SM_BQ_SPHERE )
    {
      rPoint = crX.m_sPoint; rDir = crPl.m_sDir;
      return TRUE;
    }
  if ( bAxisA || bAxisB )
    {
      if ( ( crPl.m_sDir * crX.m_sDir ).Length() > dAngTol ) { return FALSE; }
      rPoint = crX.m_sPoint; rDir = crX.m_sDir;
      return TRUE;
    }
  return FALSE;
} // end static sm_Coaxial

/*******************************************************************//**
PURPOSE: The radius of the ball about crC, a point of the axis, that the
   cylinder or cone touches along a circle (negative: crC beyond the apex).
***********************************************************************/
static double sm_InscribedRadius(const SmBlendQuadric & crQ, const SmPoint3d & crC)
{
  if ( crQ.m_eKind == SmBlendQuadric::SM_BQ_CYLINDER )
    { return crQ.m_dRadius; }
  return ( crC - crQ.m_sPoint ).Dot( crQ.m_sDir ) * smos_Sine( crQ.m_dAngle );
}

/*******************************************************************//**
PURPOSE: Two cylinders or cones (offsets) touching a common ball, crC radius
   dRho, the same way round: the plane of the conic of their intersection
   through crCentre, as a point and unit normal.

NOTES: Each is the ball's sphere S less a multiple of the square of the plane
   P of the circle it touches the ball along: S - P^2 / cos^2(a) (a cylinder:
   a = 0, P through the centre square to the axis). Their difference is
   (P1/cos a1 - P0/cos a0)(P1/cos a1 + P0/cos a0): two planes, each meeting
   both quadrics in a conic.
***********************************************************************/
static SmBoolean sm_PairPlane(const SmBlendQuadric aQ[2], const SmPoint3d & crC, double dRho, const SmPoint3d & crCentre, double dTol,
                              SmPoint3d & rPoint, SmVector3d & rNormal)
{
  SmVector3d aU[2];
  double     adK[2], adS[2];
  for ( int i=0; i<2; i++ )
    {
      // P_i(x) = (x - c).u_i - k_i, u_i from c toward the apex (along -d), k_i = rho sin(a)
      const SmBoolean bCone = aQ[i].m_eKind == SmBlendQuadric::SM_BQ_CONE;
      aU[i]  = -aQ[i].m_sDir;
      adK[i] = bCone ? dRho * smos_Sine( aQ[i].m_dAngle ) : 0.0;
      adS[i] = bCone ? 1.0 / smos_Cosine( aQ[i].m_dAngle ) : 1.0;
    }
  double dBest = SM_BIG_DOUBLE;
  for ( double dSign : { 1.0, -1.0 } )
    {
      SmVector3d sN( aU[1] * adS[1] - aU[0] * ( dSign * adS[0] ) );
      const double dK = adK[1] * adS[1] - adK[0] * dSign * adS[0];
      const double dLen = sN.Length();
      if ( dLen <= SM_EFF_ZERO_SQRT )
        { continue; }
      const double dOff = smos_Fabs( sN.Dot( crCentre - crC ) - dK ) / dLen;
      if ( dOff < dBest )
        {
          dBest   = dOff;
          rNormal = sN * ( 1.0 / dLen );
          rPoint  = crC + rNormal * ( dK / dLen );
        }
    }
  return dBest <= dTol ? TRUE : FALSE;
} // end static sm_PairPlane

/*******************************************************************//**
PURPOSE: The points of quadric crQ where the ball touching it has contact
   direction crE (from its centre to the contact; the surface's normal there
   is crE times dContact): a line of a cylinder or cone, a point of a sphere.
   FALSE where there are none.
***********************************************************************/
static SmBoolean sm_PointsWithContact(const SmBlendQuadric & crQ, double dContact, const SmVector3d & crE,
                                      SmPoint3d & rPoint, SmVector3d & rDir, SmBoolean & rbLine)
{
  // the surface's normal there: crE / dContact
  const SmVector3d sN( crE * dContact );
  switch ( crQ.m_eKind )
    {
    case SmBlendQuadric::SM_BQ_SPHERE:
      rPoint = crQ.m_sPoint + sN * ( crQ.m_dSide * crQ.m_dRadius );
      rbLine = FALSE;
      return TRUE;
    case SmBlendQuadric::SM_BQ_CYLINDER:
      if ( smos_Fabs( sN.Dot( crQ.m_sDir ) ) > SM_EFF_ZERO_SQRT )
        { return FALSE; }
      rPoint = crQ.m_sPoint + sN * ( crQ.m_dSide * crQ.m_dRadius );
      rDir   = crQ.m_sDir;
      rbLine = TRUE;
      return TRUE;
    case SmBlendQuadric::SM_BQ_CONE:
      {
        // the normal is side (cos a radial - sin a d): along the generator in the plane of d and n
        if ( smos_Fabs( sN.Dot( crQ.m_sDir ) + crQ.m_dSide * smos_Sine( crQ.m_dAngle ) ) > SM_EFF_ZERO_SQRT )
          { return FALSE; }
        SmVector3d sG( sm_Across( crQ.m_sDir, sN ) );
        if ( sG.Unitize( FALSE ) != SM_SUCCESS )
          { return FALSE; }
        rPoint = crQ.m_sPoint;
        rDir   = sG;
        rbLine = TRUE;
        return TRUE;
      }
    default:
      return FALSE;
    }
} // end static sm_PointsWithContact

/*******************************************************************//**
PURPOSE: The unit directions e with e.d0 = s0 and e.d1 = s1 (none, or two).
***********************************************************************/
static int sm_DirectionsWithDots(const SmVector3d & crD0, double dS0, const SmVector3d & crD1, double dS1, SmVector3d aE[2])
{
  const double dC = crD0.Dot( crD1 );
  const double dDen = 1.0 - dC * dC;
  if ( dDen <= SM_EFF_ZERO_SQRT )
    { return 0; }
  const double dAl = ( dS0 - dC * dS1 ) / dDen;
  const double dBe = ( dS1 - dC * dS0 ) / dDen;
  SmVector3d sIn( crD0 * dAl + crD1 * dBe );
  SmVector3d sX( crD0 * crD1 );
  const double dG2 = ( 1.0 - sIn.Dot( sIn ) ) / sX.Dot( sX );
  if ( dG2 < 0.0 )
    { return 0; }
  const double dG = smos_Sqrt( dG2 );
  aE[0] = sIn + sX * dG;
  aE[1] = sIn - sX * dG;
  return 2;
} // end static sm_DirectionsWithDots

/*******************************************************************//**
PURPOSE: A point where two offsets (cylinders, cones, or one and a sphere)
   touch with the ball's contact directions agreeing, within dTol: the node of
   their quartic intersection. FALSE where there is none.
***********************************************************************/
static SmBoolean sm_TouchingPoint(const SmBlendQuadric aQ[2], const double adContact[2], double dTol, SmPoint3d & rNode)
{
  typedef SmBlendQuadric Q;
  // the candidate contact directions e: from the kinds' normal constraints
  SmVector3d aE[4];
  int nE = 0;
  const Q * pS = ( aQ[0].m_eKind == Q::SM_BQ_SPHERE ) ? &aQ[0] : ( aQ[1].m_eKind == Q::SM_BQ_SPHERE ? &aQ[1] : NULL );
  const int iS = ( pS == &aQ[0] ) ? 0 : 1;
  if ( pS != NULL )
    {
      const Q & crX = aQ[1 - iS];
      if ( crX.m_eKind != Q::SM_BQ_CYLINDER && crX.m_eKind != Q::SM_BQ_CONE )
        { return FALSE; }
      // e in the plane through the axis and the centre
      SmVector3d sW( sm_Across( pS->m_sPoint - crX.m_sPoint, crX.m_sDir ) );
      if ( sW.Unitize( FALSE ) != SM_SUCCESS )
        { return FALSE; } // the centre on the axis: coaxial
      const double dDot = ( crX.m_eKind == Q::SM_BQ_CYLINDER ) ? 0.0 : -adContact[1 - iS] * crX.m_dSide * smos_Sine( crX.m_dAngle );
      const double dAcross = smos_Sqrt( smos_Max( 0.0, 1.0 - dDot * dDot ) );
      aE[nE++] = crX.m_sDir * dDot + sW * dAcross;
      aE[nE++] = crX.m_sDir * dDot - sW * dAcross;
    }
  else
    {
      if (   ( aQ[0].m_eKind != Q::SM_BQ_CYLINDER && aQ[0].m_eKind != Q::SM_BQ_CONE )
          || ( aQ[1].m_eKind != Q::SM_BQ_CYLINDER && aQ[1].m_eKind != Q::SM_BQ_CONE ) )
        { return FALSE; }
      double adDot[2];
      for ( int i=0; i<2; i++ )
        { adDot[i] = ( aQ[i].m_eKind == Q::SM_BQ_CYLINDER ) ? 0.0 : -adContact[i] * aQ[i].m_dSide * smos_Sine( aQ[i].m_dAngle ); }
      nE = sm_DirectionsWithDots( aQ[0].m_sDir, adDot[0], aQ[1].m_sDir, adDot[1], aE );
    }

  for ( int k=0; k<nE; k++ )
    {
      SmVector3d sE( aE[k] );
      if ( sE.Unitize( FALSE ) != SM_SUCCESS )
        { continue; }
      SmPoint3d  aP[2];
      SmVector3d aD[2];
      SmBoolean  abLine[2];
      if (   !sm_PointsWithContact( aQ[0], adContact[0], sE, aP[0], aD[0], abLine[0] )
          || !sm_PointsWithContact( aQ[1], adContact[1], sE, aP[1], aD[1], abLine[1] ) )
        { continue; }
      SmPoint3d sQ;
      if ( abLine[0] && abLine[1] )
        {
          // both lines lie square to e: they meet when they lie in one plane
          if ( smos_Fabs( ( aP[1] - aP[0] ).Dot( sE ) ) > dTol )
            { continue; }
          if ( !sm_LinesMeet( aP[0], aD[0], aP[1], aD[1], dTol, sQ ) )
            { continue; }
        }
      else
        {
          const int iPt = abLine[0] ? 1 : 0;
          sQ = aP[iPt];
          if ( sm_Across( sQ - aP[1 - iPt], aD[1 - iPt] ).Length() > dTol )
            { continue; }
        }
      // on both, the ball touching each along e there
      SmBoolean bTouch = TRUE;
      for ( int i=0; i<2; i++ )
        {
          if (   aQ[i].Distance( sQ ) > dTol
              || ( aQ[i].m_eKind == Q::SM_BQ_CONE && ( sQ - aQ[i].m_sPoint ).Dot( aQ[i].m_sDir ) <= dTol )
              || ( aQ[i].Normal( sQ ) * adContact[i] ).Dot( sE ) < 1.0 - SM_EFF_ZERO_SQRT )
            { bTouch = FALSE; }
        }
      if ( bTouch )
        {
          rNode = sQ;
          return TRUE;
        }
    }
  return FALSE;
} // end static sm_TouchingPoint

/*******************************************************************//**
PURPOSE: The one classification of a constant round's spine by its faces'
   kinds and placement (see SmBlendSpine).
***********************************************************************/
SmStatus sm_ClassifyBlendSpine(const SmPoint3d & crCentre, double dTol, SmExactBlend & rBlend)
{
  typedef SmBlendQuadric Q;
  rBlend.m_eSpine = SM_BS_TRACED;
  const Q & crA = rBlend.m_aOffsets[0];
  const Q & crB = rBlend.m_aOffsets[1];
  if ( crA.m_eKind == Q::SM_BQ_NONE || crB.m_eKind == Q::SM_BQ_NONE )
    { return SM_SUCCESS; }

  // the reach: how far the spine lies from the quadrics' points (at least the ball's radius),
  // over which an angle of dTol / dReach moves a point by dTol
  double dReach = smos_Max( dTol, rBlend.m_dRadius );
  for ( const Q * pQ : { &crA, &crB } )
    {
      if ( pQ->m_eKind != Q::SM_BQ_PLANE )
        { dReach = smos_Max( dReach, crCentre.DistanceBetween( pQ->m_sPoint ) ); }
    }

  // a common axis: circles (a plane square to it, a sphere's centre on it, the same axis)
  SmPoint3d  sAxisPoint;
  SmVector3d sAxisDir;
  if ( sm_Coaxial( crA, crB, dTol, dReach, sAxisPoint, sAxisDir ) )
    {
      // two coaxial cylinders meet nowhere or all along: no spine
      if (   ( crA.m_eKind == Q::SM_BQ_CYLINDER || crA.m_eKind == Q::SM_BQ_LINE )
          && ( crB.m_eKind == Q::SM_BQ_CYLINDER || crB.m_eKind == Q::SM_BQ_LINE ) )
        { return SM_SUCCESS; }
      rBlend.m_eSpine = SM_BS_CIRCLE;
      rBlend.m_sPoint = sAxisPoint;
      rBlend.m_sDir   = sAxisDir;
      return SM_SUCCESS;
    }

  const double dAngTol = dTol / smos_Max( dReach, dTol );
  const Q::Kind eA = crA.m_eKind, eB = crB.m_eKind;

  // two planes: a line
  if ( eA == Q::SM_BQ_PLANE && eB == Q::SM_BQ_PLANE )
    {
      if ( ( crA.m_sDir * crB.m_sDir ).Length() > dAngTol )
        { rBlend.m_eSpine = SM_BS_LINE; }
      return SM_SUCCESS;
    }

  // a plane and a cylinder or cone: a line along a cylinder parallel to the plane, else a conic
  if ( eA == Q::SM_BQ_PLANE || eB == Q::SM_BQ_PLANE )
    {
      const Q & crPl = ( eA == Q::SM_BQ_PLANE ) ? crA : crB;
      const Q & crX  = ( eA == Q::SM_BQ_PLANE ) ? crB : crA;
      if ( crX.m_eKind == Q::SM_BQ_CYLINDER )
        {
          if ( smos_Fabs( crPl.m_sDir.Dot( crX.m_sDir ) ) <= dAngTol )
            { rBlend.m_eSpine = SM_BS_LINE; }
          else
            { rBlend.m_eSpine = SM_BS_CONIC; }
        }
      else if ( crX.m_eKind == Q::SM_BQ_CONE )
        {
          // a plane through the apex meets the cone in lines: not a blend's conic
          if ( crPl.Distance( crX.m_sPoint ) > dTol )
            { rBlend.m_eSpine = SM_BS_CONIC; }
        }
      if ( rBlend.m_eSpine == SM_BS_CONIC )
        {
          rBlend.m_sPoint = crPl.m_sPoint;
          rBlend.m_sDir   = crPl.m_sDir;
        }
      return SM_SUCCESS;
    }

  // two lines (cylinders as wide as the ball) meeting: a point
  if ( eA == Q::SM_BQ_LINE && eB == Q::SM_BQ_LINE )
    {
      SmPoint3d sMeet;
      if ( sm_LinesMeet( crA.m_sPoint, crA.m_sDir, crB.m_sPoint, crB.m_sDir, dTol, sMeet ) )
        {
          rBlend.m_eSpine = SM_BS_POINT;
          rBlend.m_sPoint = sMeet;
        }
      return SM_SUCCESS;
    }
  if ( eA == Q::SM_BQ_LINE || eB == Q::SM_BQ_LINE || eA == Q::SM_BQ_POINT || eB == Q::SM_BQ_POINT )
    { return SM_SUCCESS; }

  // two cylinders with parallel axes: lines
  if ( eA == Q::SM_BQ_CYLINDER && eB == Q::SM_BQ_CYLINDER && ( crA.m_sDir * crB.m_sDir ).Length() <= dAngTol )
    {
      rBlend.m_eSpine = SM_BS_LINE;
      return SM_SUCCESS;
    }

  // two cylinders or cones touching a common ball the same way round: two plane conics
  // (the inscribed radius moves by -side * contact * r as the ball's radius r does: alike for both)
  if ( eA != Q::SM_BQ_SPHERE && eB != Q::SM_BQ_SPHERE
      && crA.m_dSide * rBlend.m_adContact[0] == crB.m_dSide * rBlend.m_adContact[1] )
    {
      SmPoint3d sC;
      if ( sm_LinesMeet( crA.m_sPoint, crA.m_sDir, crB.m_sPoint, crB.m_sDir, dTol, sC ) )
        {
          const double dRhoA = sm_InscribedRadius( crA, sC );
          const double dRhoB = sm_InscribedRadius( crB, sC );
          if ( dRhoA > dTol && smos_Fabs( dRhoA - dRhoB ) <= dTol )
            {
              const Q aQ[2] = { crA, crB };
              if ( sm_PairPlane( aQ, sC, 0.5 * ( dRhoA + dRhoB ), crCentre, dTol, rBlend.m_sPoint, rBlend.m_sDir ) )
                {
                  rBlend.m_eSpine = SM_BS_CONIC_PAIR;
                  return SM_SUCCESS;
                }
            }
        }
    }

  // touching at one point the same way round: a quartic with a node there
  {
    const Q aQ[2] = { crA, crB };
    SmPoint3d sNode;
    if ( sm_TouchingPoint( aQ, rBlend.m_adContact, dTol, sNode ) )
      {
        rBlend.m_eSpine = SM_BS_NODAL_QUARTIC;
        rBlend.m_sPoint = sNode;
      }
  }
  return SM_SUCCESS;
} // end sm_ClassifyBlendSpine

/*******************************************************************//**
                          RATIONAL SPINES AND THEIR ROUNDS

A spine is a run of rational Bezier segments, homogeneous (M, w), continuous
in its homogeneous coordinates from segment to segment. Polynomials are kept
by their Bernstein coefficients on [0, 1].
***********************************************************************/

typedef std::vector<double> sm_Poly;
struct sm_VPoly { sm_Poly a[3]; };
struct sm_Segment { sm_VPoly M; sm_Poly W; };

static double sm_Binom(int n, int k)
{
  double d = 1.0;
  for ( int i=1; i<=k; i++ ) { d = d * ( n - k + i ) / i; }
  return d;
}

static sm_Poly sm_Elevate(const sm_Poly & crP, int iTo)
{
  sm_Poly sP( crP );
  while ( static_cast<int>( sP.size() ) - 1 < iTo )
    {
      const int n = static_cast<int>( sP.size() ) - 1;
      sm_Poly sQ( n + 2 );
      sQ[0] = sP[0]; sQ[n+1] = sP[n];
      for ( int i=1; i<=n; i++ )
        { sQ[i] = ( static_cast<double>( i ) / ( n + 1 ) ) * sP[i-1] + ( 1.0 - static_cast<double>( i ) / ( n + 1 ) ) * sP[i]; }
      sP = sQ;
    }
  return sP;
}

static sm_Poly sm_Mul(const sm_Poly & crA, const sm_Poly & crB)
{
  const int m = static_cast<int>( crA.size() ) - 1, n = static_cast<int>( crB.size() ) - 1;
  sm_Poly sC( m + n + 1, 0.0 );
  for ( int i=0; i<=m; i++ )
    {
      for ( int j=0; j<=n; j++ )
        { sC[i+j] += sm_Binom( m, i ) * sm_Binom( n, j ) / sm_Binom( m + n, i + j ) * crA[i] * crB[j]; }
    }
  return sC;
}

static sm_Poly sm_Add(const sm_Poly & crA, const sm_Poly & crB, double dB = 1.0)
{
  const int n = static_cast<int>( smos_Max( crA.size(), crB.size() ) ) - 1;
  sm_Poly sA( sm_Elevate( crA, n ) ), sB( sm_Elevate( crB, n ) );
  for ( int i=0; i<=n; i++ ) { sA[i] += dB * sB[i]; }
  return sA;
}

static sm_Poly sm_Scale(const sm_Poly & crA, double d)
{
  sm_Poly s( crA );
  for ( double & x : s ) { x *= d; }
  return s;
}

static double sm_Eval(const sm_Poly & crP, double t)
{
  sm_Poly s( crP );
  for ( size_t k=1; k<s.size(); k++ )
    {
      for ( size_t i=0; i+k<s.size(); i++ ) { s[i] = ( 1.0 - t ) * s[i] + t * s[i+1]; }
    }
  return s[0];
}

// halves of a polynomial (de Casteljau at 1/2)
static void sm_Split(const sm_Poly & crP, sm_Poly & rL, sm_Poly & rR)
{
  const size_t n = crP.size();
  sm_Poly s( crP );
  rL.assign( n, 0.0 ); rR.assign( n, 0.0 );
  rL[0] = s[0]; rR[n-1] = s[n-1];
  for ( size_t k=1; k<n; k++ )
    {
      for ( size_t i=0; i+k<n; i++ ) { s[i] = 0.5 * ( s[i] + s[i+1] ); }
      rL[k] = s[0]; rR[n-1-k] = s[n-1-k];
    }
}

static sm_Poly sm_Const(double d) { return sm_Poly( 1, d ); }

static sm_VPoly sm_VConst(const SmVector3d & crV)
{
  sm_VPoly s;
  s.a[0] = sm_Const( crV.x ); s.a[1] = sm_Const( crV.y ); s.a[2] = sm_Const( crV.z );
  return s;
}
static sm_VPoly sm_VAdd(const sm_VPoly & crA, const sm_VPoly & crB, double dB = 1.0)
{
  sm_VPoly s;
  for ( int k=0; k<3; k++ ) { s.a[k] = sm_Add( crA.a[k], crB.a[k], dB ); }
  return s;
}
static sm_VPoly sm_VMul(const sm_VPoly & crA, const sm_Poly & crB)
{
  sm_VPoly s;
  for ( int k=0; k<3; k++ ) { s.a[k] = sm_Mul( crA.a[k], crB ); }
  return s;
}
static sm_VPoly sm_VScale(const sm_VPoly & crA, double d)
{
  sm_VPoly s;
  for ( int k=0; k<3; k++ ) { s.a[k] = sm_Scale( crA.a[k], d ); }
  return s;
}
static sm_Poly sm_Dot(const sm_VPoly & crA, const sm_VPoly & crB)
{
  return sm_Add( sm_Add( sm_Mul( crA.a[0], crB.a[0] ), sm_Mul( crA.a[1], crB.a[1] ) ), sm_Mul( crA.a[2], crB.a[2] ) );
}
static sm_Poly sm_DotConst(const sm_VPoly & crA, const SmVector3d & crV)
{
  return sm_Add( sm_Add( sm_Scale( crA.a[0], crV.x ), sm_Scale( crA.a[1], crV.y ) ), sm_Scale( crA.a[2], crV.z ) );
}
static SmPoint3d sm_VEval(const sm_VPoly & crA, double t)
{
  return SmPoint3d( sm_Eval( crA.a[0], t ), sm_Eval( crA.a[1], t ), sm_Eval( crA.a[2], t ) );
}

/*******************************************************************//**
PURPOSE: The contact of the ball with one face along a spine segment, as
   polynomials: the unit direction from the centre to the contact is E / D.
   bW: D is the spine's own weight w (or 1), so the contact is (M + r E) / w
   (sphere, cylinder) or (M + r w E) / w (plane, D = 1); otherwise (a cone) it
   is (M D + r w E) / (w D).
***********************************************************************/
struct sm_Contact { sm_VPoly E; sm_Poly D; SmBoolean bW; SmBoolean bPlane; };

static sm_Contact sm_ContactOf(const SmBlendQuadric & crO, double dContact, const sm_Segment & crS)
{
  sm_Contact c;
  const double dS = dContact * crO.m_dSide;
  // V = M - p w
  sm_VPoly sV = sm_VAdd( crS.M, sm_VMul( sm_VConst( crO.m_sPoint ), crS.W ), -1.0 );
  switch ( crO.m_eKind )
    {
    case SmBlendQuadric::SM_BQ_PLANE:
      c.E = sm_VConst( crO.m_sDir * dContact ); c.D = sm_Const( 1.0 ); c.bW = TRUE; c.bPlane = TRUE;
      break;
    case SmBlendQuadric::SM_BQ_SPHERE:
      c.E = sm_VScale( sV, dS / crO.m_dRadius ); c.D = crS.W; c.bW = TRUE; c.bPlane = FALSE;
      break;
    case SmBlendQuadric::SM_BQ_CYLINDER:
      {
        sm_Poly sH = sm_DotConst( sV, crO.m_sDir );
        sm_VPoly sR = sm_VAdd( sV, sm_VMul( sm_VConst( crO.m_sDir ), sH ), -1.0 );
        c.E = sm_VScale( sR, dS / crO.m_dRadius ); c.D = crS.W; c.bW = TRUE; c.bPlane = FALSE;
      }
      break;
    default: // SM_BQ_CONE: n = ((V - H d) cos^2 a - H d sin^2 a) / (H sin a)
      {
        const double dC = smos_Cosine( crO.m_dAngle ), dSn = smos_Sine( crO.m_dAngle );
        sm_Poly sH = sm_DotConst( sV, crO.m_sDir );
        sm_VPoly sR = sm_VAdd( sV, sm_VMul( sm_VConst( crO.m_sDir ), sH ), -1.0 );
        sm_VPoly sN = sm_VAdd( sm_VScale( sR, dC * dC ), sm_VMul( sm_VConst( crO.m_sDir ), sH ), -dSn * dSn );
        c.E = sm_VScale( sN, dS ); c.D = sm_Scale( sH, dSn ); c.bW = FALSE; c.bPlane = FALSE;
      }
      break;
    }
  return c;
}

/*******************************************************************//**
PURPOSE: The rows of the round's surface over one spine segment, from face
   iA's contact to face iB's: homogeneous rows (X, W), three for the ball's
   arc, two for a chamfer's line. The arc is the rational quadratic from
   contact A through the meeting point of its end tangents to contact B, its
   weights 2 G L, G L and Da Db L with G = Da Db (1 + cos), the cosine of the
   arc's angle e_A.e_B = (Ea.Eb) / (Da Db): the middle weight squared is the
   product of the end ones times (1 + cos) / 2, the arc's. L clears the
   denominators. Face iA is a plane, sphere or cylinder when either is.
***********************************************************************/
static void sm_RoundRows(const sm_Contact & crA, const sm_Contact & crB, const sm_Segment & crS, double dR, SmBoolean bArc,
                         std::vector<sm_VPoly> & rX, std::vector<sm_Poly> & rW)
{
  rX.clear(); rW.clear();
  // contact numerators: A over w (bW) or w Da
  auto num = [&]( const sm_Contact & c ) {
    if ( c.bPlane )
      { return sm_VAdd( crS.M, sm_VMul( c.E, crS.W ), dR ); }
    if ( c.bW )
      { return sm_VAdd( crS.M, c.E, dR ); }
    return sm_VAdd( sm_VMul( crS.M, c.D ), sm_VMul( c.E, crS.W ), dR );
  };
  const sm_VPoly sNA = num( crA ), sNB = num( crB );
  if ( !bArc )
    {
      // a chamfer: the line from contact to contact, each row over its own denominator
      rX.push_back( sNA ); rW.push_back( crA.bW ? crS.W : sm_Mul( crS.W, crA.D ) );
      rX.push_back( sNB ); rW.push_back( crB.bW ? crS.W : sm_Mul( crS.W, crB.D ) );
      return;
    }
  const sm_Poly sDaDb = sm_Mul( crA.D, crB.D );
  const sm_Poly sG = sm_Add( sDaDb, sm_Dot( crA.E, crB.E ) );
  // r w (Ea Db + Eb Da), the middle row's offset from G M
  const sm_VPoly sMid = sm_VMul( sm_VAdd( sm_VMul( crA.E, crB.D ), sm_VMul( crB.E, crA.D ) ), sm_Scale( crS.W, dR ) );
  if ( crA.bW )
    {
      // L = w
      rX.push_back( sm_VScale( sm_VMul( sNA, sG ), 2.0 ) ); rW.push_back( sm_Scale( sm_Mul( sG, crS.W ), 2.0 ) );
      rX.push_back( sm_VAdd( sm_VMul( crS.M, sG ), sMid ) ); rW.push_back( sm_Mul( sG, crS.W ) );
      if ( crB.bW )
        { rX.push_back( sm_VMul( sNB, sDaDb ) ); }
      else
        { rX.push_back( sm_VMul( sNB, crA.D ) ); }
      rW.push_back( sm_Mul( sDaDb, crS.W ) );
    }
  else
    {
      // both cones: L = w Da
      rX.push_back( sm_VScale( sm_VMul( sNA, sG ), 2.0 ) ); rW.push_back( sm_Scale( sm_Mul( sm_Mul( sG, crS.W ), crA.D ), 2.0 ) );
      rX.push_back( sm_VMul( sm_VAdd( sm_VMul( crS.M, sG ), sMid ), crA.D ) ); rW.push_back( sm_Mul( sm_Mul( sG, crS.W ), crA.D ) );
      const sm_Poly sDa2 = sm_Mul( crA.D, crA.D );
      rX.push_back( sm_VMul( sNB, sDa2 ) ); rW.push_back( sm_Mul( sm_Mul( sDa2, crB.D ), crS.W ) );
    }
}

/*******************************************************************//**
PURPOSE: The conic in which a plane (crPoint, unit crNormal) meets the
   quadric crQ, as its centre and axes: ellipse C + cos t A1 + sin t A2,
   hyperbola C + s (cosh t A1 + sinh t A2) on the branch s = +-1, parabola
   C + t A2 + t^2 A1.
***********************************************************************/
enum sm_ConicKind { SM_CK_ELLIPSE, SM_CK_HYPERBOLA, SM_CK_PARABOLA };
struct sm_Conic { sm_ConicKind eKind; SmPoint3d sC; SmVector3d sA1, sA2; double dBranch; };

static SmStatus sm_PlaneSection(const SmPoint3d & crPoint, const SmVector3d & crNormal, const SmBlendQuadric & crQ, sm_Conic & rConic)
{
  // a frame of the plane
  SmVector3d sU( smos_Fabs( crNormal.x ) < 0.9 ? SmVector3d( 1, 0, 0 ) : SmVector3d( 0, 1, 0 ) );
  sU = sU - crNormal * sU.Dot( crNormal );
  SER( sU.Unitize() );
  SmVector3d sV( crNormal * sU );
  double aA[3][3], aB[3], dC;
  crQ.Implicit( aA, aB, dC );
  auto quad = [&]( const SmVector3d & a, const SmVector3d & b ) {
    const double av[3] = { a.x, a.y, a.z }, bv[3] = { b.x, b.y, b.z };
    double s = 0.0;
    for ( int i=0; i<3; i++ ) { for ( int j=0; j<3; j++ ) { s += av[i] * aA[i][j] * bv[j]; } }
    return s;
  };
  const SmVector3d sBv( aB[0], aB[1], aB[2] );
  // q(x, y) = [x y] A2 [x y] + 2 b2.[x y] + c2 about crPoint
  const double a11 = quad( sU, sU ), a12 = quad( sU, sV ), a22 = quad( sV, sV );
  const SmVector3d sG( SmVector3d( quad( SmVector3d(1,0,0), crPoint ), quad( SmVector3d(0,1,0), crPoint ), quad( SmVector3d(0,0,1), crPoint ) ) + sBv );
  const double b1 = sG.Dot( sU ), b2 = sG.Dot( sV );
  const double c2 = quad( crPoint, crPoint ) + 2.0 * sBv.Dot( crPoint ) + dC;
  // principal axes of A2
  const double dTr = a11 + a22, dDet = a11 * a22 - a12 * a12;
  const double dDisc = smos_Sqrt( smos_Max( 0.0, 0.25 * dTr * dTr - dDet ) );
  double adL[2] = { 0.5 * dTr + dDisc, 0.5 * dTr - dDisc };
  double adE[2][2];
  for ( int k=0; k<2; k++ )
    {
      double ex = a12, ey = adL[k] - a11;
      if ( smos_Fabs( ex ) + smos_Fabs( ey ) < SM_EFF_ZERO * ( smos_Fabs( a11 ) + smos_Fabs( a22 ) + 1.0 ) )
        { ex = ( k == 0 ) == ( a11 >= a22 ) ? 1.0 : 0.0; ey = 1.0 - ex; }
      const double l = smos_Sqrt( ex * ex + ey * ey );
      adE[k][0] = ex / l; adE[k][1] = ey / l;
    }
  const SmVector3d aE[2] = { sU * adE[0][0] + sV * adE[0][1], sU * adE[1][0] + sV * adE[1][1] };
  const double adBeta[2] = { b1 * adE[0][0] + b2 * adE[0][1], b1 * adE[1][0] + b2 * adE[1][1] };
  const double dScale = smos_Max( smos_Fabs( adL[0] ), smos_Fabs( adL[1] ) );
  if ( dScale <= 0.0 )
    { SER_MSG( SM_ERR, _T("sm_PlaneSection(): the plane meets the quadric in no conic") ); }
  // a parabola: one principal coefficient zero
  for ( int k=0; k<2; k++ )
    {
      if ( smos_Fabs( adL[k] ) <= SM_EFF_ZERO * dScale )
        {
          // L t^2 + 2 beta_k s + 2 beta_o t + c2 = 0 along e_k (s) and the other (t)
          const int o = 1 - k;
          if ( smos_Fabs( adBeta[k] ) <= SM_EFF_ZERO * dScale )
            { SER_MSG( SM_ERR, _T("sm_PlaneSection(): the plane meets the quadric in lines") ); }
          // s = -(L t^2 + 2 beta_o t + c2) / (2 beta_k): about t0 = -beta_o / L
          const double t0 = -adBeta[o] / adL[o];
          const double s0 = -( adL[o] * t0 * t0 + 2.0 * adBeta[o] * t0 + c2 ) / ( 2.0 * adBeta[k] );
          rConic.eKind   = SM_CK_PARABOLA;
          rConic.sC      = crPoint + aE[k] * s0 + aE[o] * t0;
          rConic.sA1     = aE[k] * ( -adL[o] / ( 2.0 * adBeta[k] ) );
          rConic.sA2     = aE[o];
          rConic.dBranch = 1.0;
          return SM_SUCCESS;
        }
    }
  // centred: (x - x0).A2.(x - x0) + q0 = 0
  const double adX0[2] = { -adBeta[0] / adL[0], -adBeta[1] / adL[1] };
  const double q0 = c2 + adBeta[0] * adX0[0] + adBeta[1] * adX0[1];
  rConic.sC = crPoint + aE[0] * adX0[0] + aE[1] * adX0[1];
  const double adQ[2] = { -q0 / adL[0], -q0 / adL[1] };
  if ( adQ[0] > 0.0 && adQ[1] > 0.0 )
    {
      rConic.eKind = SM_CK_ELLIPSE;
      rConic.sA1 = aE[0] * smos_Sqrt( adQ[0] );
      rConic.sA2 = aE[1] * smos_Sqrt( adQ[1] );
    }
  else if ( adQ[0] > 0.0 || adQ[1] > 0.0 )
    {
      const int k = adQ[0] > 0.0 ? 0 : 1;
      rConic.eKind = SM_CK_HYPERBOLA;
      rConic.sA1 = aE[k] * smos_Sqrt( adQ[k] );
      rConic.sA2 = aE[1-k] * smos_Sqrt( -adQ[1-k] );
    }
  else
    { SER_MSG( SM_ERR, _T("sm_PlaneSection(): the plane misses the quadric") ); }
  rConic.dBranch = 1.0;
  return SM_SUCCESS;
}

static SmPoint3d sm_ConicPoint(const sm_Conic & crK, double t)
{
  switch ( crK.eKind )
    {
    case SM_CK_ELLIPSE:   return crK.sC + crK.sA1 * smos_Cosine( t ) + crK.sA2 * smos_Sine( t );
    case SM_CK_HYPERBOLA: return crK.sC + ( crK.sA1 * std::cosh( t ) + crK.sA2 * std::sinh( t ) ) * crK.dBranch;
    default:              return crK.sC + crK.sA2 * t + crK.sA1 * ( t * t );
    }
}

// the conic's parameter of a point on it (or near it)
static double sm_ConicParam(const sm_Conic & crK, const SmPoint3d & crX)
{
  const SmVector3d sD( crX - crK.sC );
  const double x = sD.Dot( crK.sA1 ) / crK.sA1.Dot( crK.sA1 ), y = sD.Dot( crK.sA2 ) / crK.sA2.Dot( crK.sA2 );
  switch ( crK.eKind )
    {
    case SM_CK_ELLIPSE:   return smos_ArcTangent2( y, x );
    case SM_CK_HYPERBOLA: return std::asinh( y * crK.dBranch );
    default:              return y;
    }
}

// one rational quadratic segment of the conic from t0 to t1 (|t1 - t0| <= pi / 2 on an ellipse)
static sm_Segment sm_ConicSegment(const sm_Conic & crK, double t0, double t1)
{
  const SmPoint3d sP0 = sm_ConicPoint( crK, t0 ), sP2 = sm_ConicPoint( crK, t1 );
  const double tm = 0.5 * ( t0 + t1 ), h = 0.5 * ( t1 - t0 );
  SmPoint3d sP1;
  double w = 1.0;
  switch ( crK.eKind )
    {
    case SM_CK_ELLIPSE:
      w = smos_Cosine( h );
      sP1 = crK.sC + ( crK.sA1 * smos_Cosine( tm ) + crK.sA2 * smos_Sine( tm ) ) * ( 1.0 / w );
      break;
    case SM_CK_HYPERBOLA:
      w = std::cosh( h );
      sP1 = crK.sC + ( crK.sA1 * std::cosh( tm ) + crK.sA2 * std::sinh( tm ) ) * ( crK.dBranch / w );
      break;
    default:
      // P(t0) + (t1 - t0)/2 P'(t0)
      sP1 = sP0 + ( crK.sA2 + crK.sA1 * ( 2.0 * t0 ) ) * h;
      break;
    }
  sm_Segment s;
  const SmPoint3d aP[3] = { sP0, sP1 * w, sP2 };
  for ( int k=0; k<3; k++ ) { s.M.a[k].resize( 3 ); }
  for ( int i=0; i<3; i++ )
    {
      s.M.a[0][i] = aP[i].x; s.M.a[1][i] = aP[i].y; s.M.a[2][i] = aP[i].z;
    }
  s.W = { 1.0, w, 1.0 };
  return s;
}

/*******************************************************************//**
PURPOSE: The power-basis polynomial on [a, b] as Bernstein coefficients on
   [0, 1].
***********************************************************************/
static sm_Poly sm_PowerToBernstein(const std::vector<double> & crP, double a, double b)
{
  const int n = static_cast<int>( crP.size() ) - 1;
  // compose with s = a + (b - a) t: Taylor coefficients at a, scaled
  std::vector<double> sC( crP );
  std::vector<double> sT( n + 1, 0.0 );
  for ( int k=0; k<=n; k++ )
    {
      // k-th Taylor coefficient at a: sum_j C(j, k) p_j a^(j-k)
      double s = 0.0;
      for ( int j=k; j<=n; j++ ) { s += sm_Binom( j, k ) * sC[j] * std::pow( a, j - k ); }
      sT[k] = s * std::pow( b - a, k );
    }
  sm_Poly sB( n + 1, 0.0 );
  for ( int i=0; i<=n; i++ )
    {
      for ( int j=0; j<=i; j++ ) { sB[i] += sm_Binom( i, j ) / sm_Binom( n, j ) * sT[j]; }
    }
  return sB;
}

/*******************************************************************//**
PURPOSE: The nodal quartic two offsets meet in, through crNode where they
   touch, from the direction crFrom (from the node) to crTo, as power-basis
   coefficients in s in [0, 1] of X (x, y, z) and W: the second point of the
   line from the node along each direction of the cone of directions the two
   quadrics share there.

NOTES: Q1(q + s v) = 2 s g1.v + s^2 v.A1.v and Q2 likewise with g2 = k g1 at
   the node: they meet again where v.(A2 - k A1).v = 0, a cone. That cone is
   swept rationally by the second point on it of the lines through one of its
   rays v0 (a branch's tangent there: g1.v0 = 0) in the directions w(s) =
   (1 - s) from + s to (from and to on the cone, either way round: the two
   arcs between them); the point of the quartic along v is
   q - 2 (g1.v) / (v.A1.v) v. iBase picks the branch whose tangent is v0;
   aBranch: both branches' tangents.
***********************************************************************/
static SmStatus sm_NodalQuartic(const SmBlendQuadric & crO1, const SmBlendQuadric & crO2, const SmPoint3d & crNode,
                                const SmVector3d & crFrom, const SmVector3d & crTo, int iBase, SmVector3d aBranch[2],
                                std::vector<double> aX[3], std::vector<double> & rW, std::vector<double> & rG)
{
  double aA1[3][3], aB1[3], dC1, aA2[3][3], aB2[3], dC2;
  crO1.Implicit( aA1, aB1, dC1 );
  crO2.Implicit( aA2, aB2, dC2 );
  const double q[3] = { crNode.x, crNode.y, crNode.z };
  double g1[3], g2[3];
  for ( int i=0; i<3; i++ )
    {
      g1[i] = aB1[i]; g2[i] = aB2[i];
      for ( int j=0; j<3; j++ ) { g1[i] += aA1[i][j] * q[j]; g2[i] += aA2[i][j] * q[j]; }
    }
  const double dG11 = g1[0]*g1[0] + g1[1]*g1[1] + g1[2]*g1[2];
  if ( dG11 <= 0.0 )
    { SER_MSG( SM_ERR, _T("sm_NodalQuartic(): no normal at the node") ); }
  const double dK = ( g1[0]*g2[0] + g1[1]*g2[1] + g1[2]*g2[2] ) / dG11;
  double aC[3][3];
  for ( int i=0; i<3; i++ ) { for ( int j=0; j<3; j++ ) { aC[i][j] = aA2[i][j] - dK * aA1[i][j]; } }
  auto bil = [&]( const double (&M)[3][3], const double a[3], const double b[3] ) {
    double s = 0.0;
    for ( int i=0; i<3; i++ ) { for ( int j=0; j<3; j++ ) { s += a[i] * M[i][j] * b[j]; } }
    return s;
  };
  // a ray of the cone square to g1: a branch's tangent
  SmVector3d sN( g1[0], g1[1], g1[2] );
  SER( sN.Unitize() );
  SmVector3d sT1( smos_Fabs( sN.x ) < 0.9 ? SmVector3d( 1, 0, 0 ) : SmVector3d( 0, 1, 0 ) );
  sT1 = sT1 - sN * sT1.Dot( sN );
  SER( sT1.Unitize() );
  const SmVector3d sT2( sN * sT1 );
  const double t1[3] = { sT1.x, sT1.y, sT1.z }, t2[3] = { sT2.x, sT2.y, sT2.z };
  const double c11 = bil( aC, t1, t1 ), c12 = bil( aC, t1, t2 ), c22 = bil( aC, t2, t2 );
  // c11 cos^2 + 2 c12 cos sin + c22 sin^2 = 0: the two branches' tangents
  const double d = c12 * c12 - c11 * c22;
  if ( d < 0.0 || ( c11 == 0.0 && c22 == 0.0 ) )
    { SER_MSG( SM_ERR, _T("sm_NodalQuartic(): the node is isolated") ); }
  for ( int k=0; k<2; k++ )
    {
      const double dRoot = ( k == 0 ? 1.0 : -1.0 ) * smos_Sqrt( d );
      // (cos, sin) along (c22, -c12 + root) or (-c12 - root, c11), the better conditioned
      SmVector3d sB( smos_Fabs( c22 ) >= smos_Fabs( c11 ) ? sT1 * c22 + sT2 * ( -c12 + dRoot ) : sT1 * ( -c12 - dRoot ) + sT2 * c11 );
      SER( sB.Unitize() );
      aBranch[k] = sB;
    }
  const SmVector3d sV0( aBranch[iBase] );
  const double v0[3] = { sV0.x, sV0.y, sV0.z };
  const double w0[3] = { crFrom.x, crFrom.y, crFrom.z };
  const double w1[3] = { crTo.x - w0[0], crTo.y - w0[1], crTo.z - w0[2] };
  // v(s) = Qc(w) v0 - 2 B(v0, w) w, w = w0 + s w1
  const double q0c = bil( aC, w0, w0 ), q1c = 2.0 * bil( aC, w0, w1 ), q2c = bil( aC, w1, w1 );
  const double be0 = bil( aC, v0, w0 ), be1 = bil( aC, v0, w1 );
  double vc[3][3]; // [power][coord]
  for ( int i=0; i<3; i++ )
    {
      vc[0][i] = q0c * v0[i] - 2.0 * be0 * w0[i];
      vc[1][i] = q1c * v0[i] - 2.0 * ( be0 * w1[i] + be1 * w0[i] );
      vc[2][i] = q2c * v0[i] - 2.0 * be1 * w1[i];
    }
  // G = g1.v (quadratic), W = v.A1.v (quartic)
  rG.assign( 3, 0.0 );
  for ( int p=0; p<3; p++ ) { for ( int i=0; i<3; i++ ) { rG[p] += g1[i] * vc[p][i]; } }
  rW.assign( 5, 0.0 );
  for ( int p=0; p<3; p++ ) { for ( int r=0; r<3; r++ ) { rW[p+r] += bil( aA1, vc[p], vc[r] ); } }
  // X = q W - 2 G v
  for ( int i=0; i<3; i++ )
    {
      aX[i].assign( 5, 0.0 );
      for ( int p=0; p<5; p++ ) { aX[i][p] = q[i] * rW[p]; }
      for ( int p=0; p<3; p++ ) { for ( int r=0; r<3; r++ ) { aX[i][p+r] -= 2.0 * rG[p] * vc[r][i]; } }
    }
  return SM_SUCCESS;
}

static double sm_PowerEval(const std::vector<double> & crP, double s)
{
  double v = 0.0;
  for ( size_t i=crP.size(); i-- > 0; ) { v = v * s + crP[i]; }
  return v;
}

/*******************************************************************//**
PURPOSE: The spine of a rational class, as segments from the ball's centre at
   the trace's first point to its last (round the whole loop on a closed edge),
   extended at each end as the solver asks.
***********************************************************************/
SmStatus SmFilletSolver::CalcRationalSpineFilletGeom(const SmExactBlend & crBlend, SmTArray<SmTsectPnt*> & rTsectPnts,
                                                     SmFilletGeom * pFilletGeom)
{
  if ( pFilletGeom == NULL )
    { pFilletGeom = m_vFilletGeoms[0]; }
  const ULONG lLast = rTsectPnts.GetSize() - 1;
  SmEdge  * pEdge  = GetEdgeuse(0)->GetEdge(); NER( pEdge );
  SmCurve * pCurve = pEdge->GetCurve();         NER( pCurve );
  const SmBrep * pBrep = pEdge->GetBrep();
  const double dTol = ( pBrep != NULL ) ? static_cast<double>( pBrep->GetTolerance() ) : GetThisApproxTol3d();
  const double dR = crBlend.m_dRadius;
  const SmBoolean bArc = m_pFSG->GetFilletSurfaceGeneratorType() == SM_FSG_CIRCULAR;

  // the ball's centre at the trace's ends, and the edge's direction there
  SmPoint3d sStart, sEnd;
  SER( SetupOffsetValues( *rTsectPnts[0], NULL ) );
  SER( GetSurface(0)->EvaluatePoint( rTsectPnts[0]->UVPos(0), sStart ) );
  SER( SetupOffsetValues( *rTsectPnts[lLast], NULL ) );
  SER( GetSurface(0)->EvaluatePoint( rTsectPnts[lLast]->UVPos(0), sEnd ) );
  const SmBoolean bClosed = pCurve->IsClosed( pEdge->GetInterval() );
  // the edge's middle, and its direction at the trace's first point (for a closed edge, the way
  // round the trace runs)
  SmPoint3d  sEdgeMid;
  SmVector3d sEdgeDir;
  {
    SER( pCurve->EvaluatePoint( pEdge->GetInterval().Evaluate( 0.5 ), sEdgeMid ) );
    SmSolutionArray sSol;
    SER( pCurve->GlobalPointSolve( pEdge->GetInterval(), SM_SO_MINIMIZE, sStart, SM_EFF_ZERO_SQRT, NULL, NULL, SM_SR_SINGLE, sSol ) );
    if ( sSol.GetSize() < 1 ) { SER( SM_ERR ); }
    SmPoint3d sPV[2];
    SER( pCurve->Evaluate( sSol[0].m_vStart[0], 1, TRUE, sPV ) );
    sEdgeDir = sPV[1];
  }

  std::vector<sm_Segment> sSegs;
  const SmBlendQuadric & crO0 = crBlend.m_aOffsets[0];
  const SmBlendQuadric & crO1 = crBlend.m_aOffsets[1];
  if ( crBlend.m_eSpine == SM_BS_CONIC || crBlend.m_eSpine == SM_BS_CONIC_PAIR )
    {
      // the conic: the plane's section of the other quadric (either, for a pair)
      const SmBlendQuadric & crQ = ( crO0.m_eKind == SmBlendQuadric::SM_BQ_PLANE ) ? crO1 : crO0;
      sm_Conic sK;
      SER( sm_PlaneSection( crBlend.m_sPoint, crBlend.m_sDir, crQ, sK ) );
      if ( sK.eKind == SM_CK_HYPERBOLA && ( sStart - sK.sC ).Dot( sK.sA1 ) < 0.0 )
        { sK.dBranch = -1.0; }
      double t0 = sm_ConicParam( sK, sStart ), t1 = sm_ConicParam( sK, sEnd );
      // from the first point to the last: round the loop the way the edge runs on a closed edge,
      // else the way along an ellipse that passes the edge's middle
      double dSpan = t1 - t0;
      double dDir  = ( dSpan >= 0.0 ) ? 1.0 : -1.0;
      if ( sK.eKind == SM_CK_ELLIPSE )
        {
          if ( bClosed )
            {
              const double h = 1e-6;
              const SmVector3d sTan( sm_ConicPoint( sK, t0 + h ) - sm_ConicPoint( sK, t0 - h ) );
              dDir  = ( sTan.Dot( sEdgeDir ) >= 0.0 ) ? 1.0 : -1.0;
              dSpan = dDir * 2.0 * SM_PI;
            }
          else
            {
              // the arc t0 -> t1 one way or the other, by where the edge's middle lies
              const double dMid = sm_ConicParam( sK, sEdgeMid );
              auto within = [&]( double a, double span, double t ) {
                double d = ( t - a ) * ( span >= 0.0 ? 1.0 : -1.0 );
                d -= 2.0 * SM_PI * smos_Floor( d / ( 2.0 * SM_PI ) );
                return d <= smos_Fabs( span ); };
              while ( dSpan <= 0.0 ) { dSpan += 2.0 * SM_PI; }
              while ( dSpan > 2.0 * SM_PI ) { dSpan -= 2.0 * SM_PI; }
              if ( !within( t0, dSpan, dMid ) )
                { dSpan -= 2.0 * SM_PI; }
              dDir = ( dSpan >= 0.0 ) ? 1.0 : -1.0;
            }
        }
      // extensions
      if ( !bClosed )
        {
          double dExt = m_dSurfaceExtensionFactor * smos_Fabs( dSpan );
          if ( sK.eKind == SM_CK_ELLIPSE )
            { dExt = smos_Min( dExt, 0.25 * ( 2.0 * SM_PI - smos_Fabs( dSpan ) ) ); }
          if ( m_bExtendBefore ) { t0 -= dDir * dExt; dSpan += dDir * dExt; }
          if ( m_bExtendAfter )  { dSpan += dDir * dExt; }
        }
      const double dMax = ( sK.eKind == SM_CK_ELLIPSE ) ? 0.5 * SM_PI : 0.5;
      ULONG n = static_cast<ULONG>( smos_Ceil( smos_Fabs( dSpan ) / dMax ) );
      if ( sK.eKind == SM_CK_PARABOLA ) { n = 4; }
      if ( n < 2 ) { n = 2; }
      for ( ULONG i=0; i<n; i++ )
        { sSegs.push_back( sm_ConicSegment( sK, t0 + dSpan * i / n, t0 + dSpan * ( i + 1 ) / n ) ); }
    }
  else if ( crBlend.m_eSpine == SM_BS_NODAL_QUARTIC )
    {
      // from the node toward the trace's ends; an end at the node (the round fading out there)
      // leaves it along a branch's tangent
      const SmVector3d sFromPt( sStart - crBlend.m_sPoint ), sToPt( sEnd - crBlend.m_sPoint );
      const SmBoolean bFromNode = sFromPt.Length() <= dTol, bToNode = sToPt.Length() <= dTol;
      if ( bFromNode && bToNode )
        { SER_MSG( SM_ERR, _T("CalcRationalSpineFilletGeom(): the trace starts and ends at the node") ); }
      std::vector<double> aX[3], sW, sG;
      SmVector3d aBranch[2];
      SER( sm_NodalQuartic( crO0, crO1, crBlend.m_sPoint, sFromPt, sToPt, 0, aBranch, aX, sW, sG ) );
      // the arc from start to end that keeps off the node: G (g1.v) is not zero inside it, nor W
      SmBoolean bFound = FALSE;
      for ( int iCand=0; iCand<4 && !bFound; iCand++ )
        {
          const int iEndBranch = iCand / 2;
          const double dSign = ( iCand % 2 == 0 ) ? 1.0 : -1.0;
          const SmVector3d sFrom = ( bFromNode ? aBranch[iEndBranch] : sFromPt ) * dSign;
          const SmVector3d sTo   = bToNode ? aBranch[iEndBranch] : sToPt;
          if ( !bFromNode && !bToNode && iEndBranch == 1 )
            { continue; }
          // the base ray: the branch no end leaves along
          const int iBase = ( bFromNode || bToNode ) ? 1 - iEndBranch : 0;
          SER( sm_NodalQuartic( crO0, crO1, crBlend.m_sPoint, sFrom, sTo, iBase, aBranch, aX, sW, sG ) );
          SmBoolean bOk = TRUE;
          const double g = sm_PowerEval( sG, 0.5 ), w = sm_PowerEval( sW, 0.5 );
          for ( int k=1; k<64 && bOk; k++ )
            {
              const double t = k / 64.0;
              if ( !( sm_PowerEval( sG, t ) * g > 0.0 ) || !( sm_PowerEval( sW, t ) * w > 0.0 ) ) { bOk = FALSE; }
            }
          if ( !( sm_PowerEval( sW, 0.0 ) * w > 0.0 ) || !( sm_PowerEval( sW, 1.0 ) * w > 0.0 ) ) { bOk = FALSE; }
          bFound = bOk;
        }
      if ( !bFound )
        { SER_MSG( SM_ERR, _T("CalcRationalSpineFilletGeom(): the quartic's arc between the trace's ends runs through its node") ); }
      // s runs from the trace's first point to its last
      double s0 = 0.0, s1 = 1.0;
      {
        // extend, as far as the node and the point at infinity allow
        const double g0 = sm_PowerEval( sG, 0.5 ), w0 = sm_PowerEval( sW, 0.5 );
        auto clear = [&]( double s ) { return sm_PowerEval( sG, s ) * g0 > 0.0 && sm_PowerEval( sW, s ) * w0 > 0.0; };
        const double dExt = m_dSurfaceExtensionFactor;
        for ( int k=1; k<=32; k++ ) { const double s = -dExt * k / 32.0; if ( !m_bExtendBefore || !clear( s - dExt / 32.0 ) ) break; s0 = s; }
        for ( int k=1; k<=32; k++ ) { const double s = 1.0 + dExt * k / 32.0; if ( !m_bExtendAfter || !clear( s + dExt / 32.0 ) ) break; s1 = s; }
      }
      const double dWSign = ( sm_PowerEval( sW, 0.5 ) > 0.0 ) ? 1.0 : -1.0;
      const ULONG n = 16;
      for ( ULONG i=0; i<n; i++ )
        {
          const double a = s0 + ( s1 - s0 ) * i / n, b = s0 + ( s1 - s0 ) * ( i + 1 ) / n;
          sm_Segment sS;
          for ( int c=0; c<3; c++ ) { sS.M.a[c] = sm_Scale( sm_PowerToBernstein( aX[c], a, b ), dWSign ); }
          sS.W = sm_Scale( sm_PowerToBernstein( sW, a, b ), dWSign );
          sSegs.push_back( sS );
        }
    }
  else
    { SER( SM_ERR_INVALID_INPUT ); }

  // each segment reparametrised (its weights times c^i) so that its end weights agree, and scaled
  // to meet the one before: the same curve, its weights 1 at every joint
  {
    double dPrev = 1.0;
    for ( sm_Segment & rS : sSegs )
      {
        const int d = static_cast<int>( rS.W.size() ) - 1;
        if ( !( rS.W[0] > 0.0 ) || !( rS.W[d] > 0.0 ) )
          { SER_MSG( SM_ERR, _T("CalcRationalSpineFilletGeom(): a spine segment's end weight is not positive") ); }
        const double c = std::pow( rS.W[0] / rS.W[d], 1.0 / d );
        double dCi = dPrev / rS.W[0];
        for ( int i=0; i<=d; i++ )
          {
            rS.W[i] *= dCi;
            for ( int k=0; k<3; k++ ) { rS.M.a[k][i] *= dCi; }
            dCi *= c;
          }
        dPrev = rS.W[d];
      }
  }

  // the rows over each segment, segments halved until every weight is positive
  const int iA = ( crO0.m_eKind == SmBlendQuadric::SM_BQ_CONE && crO1.m_eKind != SmBlendQuadric::SM_BQ_CONE ) ? 1 : 0;
  const int iB = 1 - iA;
  std::vector<std::vector<sm_VPoly>> sRowsX;
  std::vector<std::vector<sm_Poly>>  sRowsW;
  std::vector<sm_Segment> sDone;
  for ( size_t k=0; k<sSegs.size(); k++ )
    {
      std::vector<std::pair<sm_Segment, int>> sTodo( 1, std::make_pair( sSegs[k], 0 ) );
      while ( !sTodo.empty() )
        {
          sm_Segment sS = sTodo.front().first;
          const int iDepth = sTodo.front().second;
          sTodo.erase( sTodo.begin() );
          std::vector<sm_VPoly> sX;
          std::vector<sm_Poly>  sWr;
          sm_RoundRows( sm_ContactOf( crBlend.m_aOffsets[iA], crBlend.m_adContact[iA], sS ),
                        sm_ContactOf( crBlend.m_aOffsets[iB], crBlend.m_adContact[iB], sS ), sS, dR, bArc, sX, sWr );
          SmBoolean bPositive = TRUE;
          for ( const sm_Poly & w : sWr ) { for ( double c : w ) { if ( !( c > 0.0 ) ) { bPositive = FALSE; } } }
          for ( double c : sS.W ) { if ( !( c > 0.0 ) ) { bPositive = FALSE; } }
          if ( !bPositive )
            {
              if ( iDepth >= 6 )
                { SER_MSG( SM_ERR, _T("CalcRationalSpineFilletGeom(): the round's weights do not keep positive") ); }
              sm_Segment sL, sR2;
              for ( int c=0; c<3; c++ ) { sm_Split( sS.M.a[c], sL.M.a[c], sR2.M.a[c] ); }
              sm_Split( sS.W, sL.W, sR2.W );
              // in order: the left half first
              sTodo.insert( sTodo.begin(), std::make_pair( sR2, iDepth + 1 ) );
              sTodo.insert( sTodo.begin(), std::make_pair( sL, iDepth + 1 ) );
              continue;
            }
          if ( iA == 1 )
            { std::reverse( sX.begin(), sX.end() ); std::reverse( sWr.begin(), sWr.end() ); }
          sRowsX.push_back( sX ); sRowsW.push_back( sWr ); sDone.push_back( sS );
        }
    }

  // the surface: u along the spine, its knots at the segments' ends spaced by their chords
  int iDeg = 0;
  for ( size_t k=0; k<sRowsW.size(); k++ )
    {
      for ( size_t j=0; j<sRowsW[k].size(); j++ )
        {
          iDeg = static_cast<int>( smos_Max( static_cast<size_t>( iDeg ), sRowsW[k][j].size() - 1 ) );
          for ( int c=0; c<3; c++ ) { iDeg = static_cast<int>( smos_Max( static_cast<size_t>( iDeg ), sRowsX[k][j].a[c].size() - 1 ) ); }
        }
    }
  const ULONG nV = static_cast<ULONG>( sRowsW[0].size() );
  const ULONG nSeg = static_cast<ULONG>( sRowsW.size() );
  const ULONG nU = nSeg * iDeg + 1;
  SmTArray<SmPoint3d> sCPs( nU * nV );
  SmTArray<double>    sWts( nU * nV );
  for ( ULONG i=0; i<nU*nV; i++ ) { sCPs.Add( SmPoint3d( 0, 0, 0 ) ); sWts.Add( 1.0 ); }
  SmTArray<double> sUKnots, sVKnots;
  SmTArray<ULONG>  sUMult, sVMult;
  double dU = 0.0;
  sUKnots.Add( dU ); sUMult.Add( iDeg + 1 );
  for ( ULONG k=0; k<nSeg; k++ )
    {
      for ( ULONG j=0; j<nV; j++ )
        {
          const sm_Poly sW = sm_Elevate( sRowsW[k][j], iDeg );
          sm_Poly sX[3];
          for ( int c=0; c<3; c++ ) { sX[c] = sm_Elevate( sRowsX[k][j].a[c], iDeg ); }
          for ( int i=0; i<=iDeg; i++ )
            {
              const ULONG lU = k * iDeg + i;
              sCPs[lU * nV + j] = SmPoint3d( sX[0][i] / sW[i], sX[1][i] / sW[i], sX[2][i] / sW[i] );
              sWts[lU * nV + j] = sW[i];
            }
        }
      const double w0 = sm_Eval( sDone[k].W, 0.0 ), w1 = sm_Eval( sDone[k].W, 1.0 );
      const double dChord = ( sm_VEval( sDone[k].M, 1.0 ) * ( 1.0 / w1 ) ).DistanceBetween( sm_VEval( sDone[k].M, 0.0 ) * ( 1.0 / w0 ) );
      dU += smos_Max( dChord, dTol );
      sUKnots.Add( dU ); sUMult.Add( k + 1 < nSeg ? iDeg : iDeg + 1 );
    }
  // v: from contact 0 to contact 1, about as long as the section
  const double dVLen = smos_Max( dR, dTol );
  sVKnots.Add( 0.0 ); sVKnots.Add( dVLen );
  sVMult.Add( nV ); sVMult.Add( nV );
  // the weights scaled about 1 (the surface is the same); they must stay within the range a
  // NURBS takes
  {
    double dMin = SM_BIG_DOUBLE, dMax = 0.0;
    for ( ULONG i=0; i<sWts.GetSize(); i++ ) { dMin = smos_Min( dMin, sWts[i] ); dMax = smos_Max( dMax, sWts[i] ); }
    if ( !( dMin > 0.0 ) || dMax > 1e4 * dMin )
      { SER_MSG( SM_ERR, _T("CalcRationalSpineFilletGeom(): the round's weights range too far for a NURBS") ); }
    const double dScale = 1.0 / smos_Sqrt( dMin * dMax );
    for ( ULONG i=0; i<sWts.GetSize(); i++ ) { sWts[i] *= dScale; }
  }
  SmBSplineSurface * pRound = NULL;
  SER( SmBSplineSurface::CreateCanonical( m_crContext, iDeg, nV - 1, sCPs, SM_SF_UNSPECIFIED, sUMult, sVMult, sUKnots, sVKnots,
                                          SM_KT_PIECEWISE_BEZIER_KNOTS, &sWts, NULL, pRound ) );
  NER( pRound );
  SmObjDelete sDelRound( pRound );

  // the round is the ball's: its rails on their faces, every point at r from the spine
  for ( ULONG k=0; k<nSeg; k++ )
    {
      for ( int a=0; a<=4; a++ )
        {
          const double t = a / 4.0;
          const double w = sm_Eval( sDone[k].W, t );
          const SmPoint3d sM = sm_VEval( sDone[k].M, t ) * ( 1.0 / w );
          const double u = sUKnots[k] + ( sUKnots[k+1] - sUKnots[k] ) * t;
          for ( int b=0; b<=4; b++ )
            {
              SmPoint3d sX;
              SER( pRound->EvaluatePoint( SmPoint2d( u, dVLen * b / 4.0 ), sX ) );
              const double dOff = bArc ? smos_Fabs( sX.DistanceBetween( sM ) - dR ) : 0.0;
              const double dOnFace = ( b == 0 ) ? crBlend.m_aFaces[0].Distance( sX ) : ( b == 4 ? crBlend.m_aFaces[1].Distance( sX ) : 0.0 );
              if ( dOff > dTol || dOnFace > dTol )
                { SER_MSG( SM_ERR, _T("CalcRationalSpineFilletGeom(): the round built is not the ball's") ); }
            }
        }
    }

  // store it, its rails the v = min and v = max isoparameter curves, and the spine
  sDelRound.Clear();
  pFilletGeom->m_pFilletSurface = pRound;
  SmExtent2d sUVDomain = pRound->GetNaturalUVDomain();
  const double adV[2] = { sUVDomain.GetMin().y, sUVDomain.GetMax().y };
  for ( ULONG kk=0; kk<2; kk++ )
    {
      SmFilletEdge   * pRail    = pFilletGeom->GetRail( kk );
      SmBSplineCurve * pRailCrv = NULL;
      SER( pRound->CreateIsoParametricCurve( m_crContext, SM_SP_V, adV[kk], 0.0, pRailCrv ) );
      pRail->SetCurve( pRailCrv, TRUE ); // TRUE = delete preExisting Edge->Curve
      pRailCrv->SetOwner( pRail );
      pRail->SetInterval( pRailCrv->GetNaturalInterval() );
    }
  {
    const int iSpineDeg = static_cast<int>( sDone[0].W.size() ) - 1;
    SmTArray<SmPoint3d> sSP;
    SmTArray<double>    sSW;
    SmTArray<ULONG>     sSM;
    for ( ULONG k=0; k<nSeg; k++ )
      {
        for ( int i=( k == 0 ? 0 : 1 ); i<=iSpineDeg; i++ )
          {
            const double w = sDone[k].W[i];
            sSP.Add( SmPoint3d( sDone[k].M.a[0][i] / w, sDone[k].M.a[1][i] / w, sDone[k].M.a[2][i] / w ) );
            sSW.Add( w );
          }
        sSM.Add( k == 0 ? iSpineDeg + 1 : iSpineDeg );
      }
    sSM.Add( iSpineDeg + 1 );
    SmBSplineCurve * pSpine = NULL;
    SER( SmBSplineCurve::CreateCanonical( m_crContext, 3, iSpineDeg, sSP, SM_CF_UNSPECIFIED, sSM, sUKnots,
                                          SM_KT_PIECEWISE_BEZIER_KNOTS, &sSW, NULL, pSpine ) );
    pFilletGeom->m_pCenterLineCurve = pSpine;
  }
  return SM_SUCCESS;
} // end SmFilletSolver::CalcRationalSpineFilletGeom

/*******************************************************************//**
PURPOSE: The exact class of this constant round, from its faces' kinds and
   placement where the ball touches them at rTsectPnt (a point of the trace,
   with each face's uv there).

NOTES: The faces are the offset surfaces' bases; each face's offset there is
   the ball's radius along its normal, either way. A section other than the
   circle (or the chamfer's line), a radius that varies, or a face of no
   natural quadric is SM_BS_TRACED; so is a rational spine for a solver whose
   offsets may vary along the edge (a chord round's radius follows the
   faces' angle, constant only where the angle is: on a line or a circle).
***********************************************************************/
SmStatus SmFilletSolver::ClassifyExactBlend(SmTsectPnt & rTsectPnt, SmExactBlend & rBlend)
{
  rBlend = SmExactBlend();
  if ( !CanMakeAnalytic() || m_pFSG == NULL )
    { return SM_SUCCESS; }
  const SmFilletSurfaceGeneratorType eSection = m_pFSG->GetFilletSurfaceGeneratorType();
  if ( eSection != SM_FSG_CIRCULAR && eSection != SM_FSG_LINEAR )
    { return SM_SUCCESS; }

  // placements are decided within the body's tolerance
  SmEdge * pEdge = GetEdgeuse(0)->GetEdge(); NER( pEdge );
  const SmBrep * pBrep = pEdge->GetBrep();
  const double dTol = ( pBrep != NULL ) ? static_cast<double>( pBrep->GetTolerance() ) : GetThisApproxTol3d();

  double dRadius = 0.0;
  SER( SetupOffsetValues( rTsectPnt, &dRadius ) );
  SmPoint3d sCentre;
  for ( ULONG i=0; i<2; i++ )
    {
      SmOffsetSurface * pOff  = GetSurface( i );           NER( pOff );
      const SmSurface * pBase = pOff->GetBaseSurface();    NER( pBase );
      const SmPoint2d   sUV   = rTsectPnt.UVPos( i );
      SER( SmBlendQuadric::FromSurface( pBase, sUV, dTol, rBlend.m_aFaces[i] ) );
      if ( rBlend.m_aFaces[i].m_eKind == SmBlendQuadric::SM_BQ_NONE )
        { return SM_SUCCESS; }
      SmPoint3d  sP, sQ;
      SmVector3d sN;
      SER( pBase->EvaluatePoint( sUV, sP ) );
      SER( pBase->EvaluateNormal( sUV, TRUE, TRUE, sN ) );
      SER( pOff->EvaluatePoint( sUV, sQ ) );
      const double dOff = ( sQ - sP ).Dot( sN );
      rBlend.m_aOffsets[i]  = rBlend.m_aFaces[i].Offset( dOff, dTol );
      rBlend.m_adContact[i] = ( dOff > 0.0 ) ? -1.0 : 1.0;
      rBlend.m_dRadius      = smos_Fabs( dOff );
      if ( i == 0 )
        { sCentre = sQ; }
    }

  SER( sm_ClassifyBlendSpine( sCentre, dTol, rBlend ) );

  if (   OffsetRadiiCanChange()
      && rBlend.m_eSpine != SM_BS_LINE && rBlend.m_eSpine != SM_BS_CIRCLE && rBlend.m_eSpine != SM_BS_POINT )
    { rBlend.m_eSpine = SM_BS_TRACED; }
  return SM_SUCCESS;
} // end SmFilletSolver::ClassifyExactBlend

/*******************************************************************//**
PURPOSE: The exact surface of a constant round by its class, stored in pFG
   (NULL: the solver's first fillet geom) with its rails. bRecalc: the
   surface rebuilt through refined points (SmFilletGeom::ReCalcFilletGeom).
   rbBuilt is FALSE when the round is traced.
***********************************************************************/
SmStatus SmFilletSolver::CalcExactFilletGeom(const SmExactBlend & crBlend, SmTArray<SmTsectPnt*> & rTsectPnts,
                                             SmFilletGeom * pFG, SmBoolean bRecalc, SmBoolean & rbBuilt)
{
  rbBuilt = FALSE;
  if ( rTsectPnts.GetSize() == 0 )
    { return SM_SUCCESS; }
  switch ( crBlend.m_eSpine )
    {
    case SM_BS_LINE:
      rbBuilt = ( bRecalc ? CalcConeFilletGeom( rTsectPnts, pFG ) : CalcCylinderFilletGeom( rTsectPnts, pFG ) ) == SM_SUCCESS;
      break;
    case SM_BS_CIRCLE:
      rbBuilt = CalcTorusFilletGeom( rTsectPnts, pFG ) == SM_SUCCESS;
      break;
    case SM_BS_POINT:
      rbBuilt = !bRecalc && CalcSphereFilletGeom( rTsectPnts, pFG ) == SM_SUCCESS;
      break;
    case SM_BS_NODAL_QUARTIC:
      {
        // A round that runs into the node ends there tangentially: its contacts merge, the end
        // event of the trace (the tangent end), which takes the traced round only. Such a round
        // is traced until that event takes the exact one.
        const SmBrep * pBrep = GetEdgeuse(0)->GetEdge()->GetBrep();
        const double dTol = ( pBrep != NULL ) ? static_cast<double>( pBrep->GetTolerance() ) : GetThisApproxTol3d();
        for ( ULONG k : { 0UL, static_cast<ULONG>( rTsectPnts.GetSize() - 1 ) } )
          {
            SmPoint3d sCentre;
            SER( SetupOffsetValues( *rTsectPnts[k], NULL ) );
            SER( GetSurface(0)->EvaluatePoint( rTsectPnts[k]->UVPos(0), sCentre ) );
            if ( sCentre.DistanceBetween( crBlend.m_sPoint ) <= dTol )
              { return SM_SUCCESS; }
          }
      }
      SER( CalcRationalSpineFilletGeom( crBlend, rTsectPnts, pFG ) );
      rbBuilt = TRUE;
      break;
    case SM_BS_CONIC:
    case SM_BS_CONIC_PAIR:
      SER( CalcRationalSpineFilletGeom( crBlend, rTsectPnts, pFG ) );
      rbBuilt = TRUE;
      break;
    default:
      break;
    }
  return SM_SUCCESS;
} // end SmFilletSolver::CalcExactFilletGeom
