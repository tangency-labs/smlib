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
#include <SmBSplineSurface.h>
#include <SmEdge.h>
#include <SmEdgeuse.h>
#include <SmFilletGeom.h>
#include <SmFilletSolver.h>
#include <SmCone.h>
#include <SmCylinder.h>
#include <SmOffsetSurface.h>
#include <SmPlane.h>
#include <SmSphere.h>
#include <SmSurface.h>

#include <cmath>

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
    default:
      break;
    }
  return SM_SUCCESS;
} // end SmFilletSolver::CalcExactFilletGeom
