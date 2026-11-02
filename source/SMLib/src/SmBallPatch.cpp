// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

/**********************************************************************//**
* FILE NAME --- SmBallPatch.cpp
* PURPOSE: SolidCraft change. A rational patch on a sphere by its
*   definition (SmBallPatch.h).
**********************************************************************/

#include "StdAfx.h"

#ifndef __SMBALLPATCH_H__
#include <SmBallPatch.h>
#endif

#include <SmAxis2Placement.h>
#include <SmBSplineCurve.h>
#include <SmDatabaseIO.h>

#include <cmath>
#include <istream>
#include <ostream>

SmBallPatch::SmBallPatch()
  : m_sCentre(0.0, 0.0, 0.0), m_dRadius(0.0)
{
}

SmBallPatch::SmBallPatch(const SmBSplineSurface & crNet, const SmPoint3d & crCentre, double dRadius)
  : SmBSplineSurface(crNet), m_sCentre(crCentre), m_dRadius(dRadius)
{
}

SmBallPatch::SmBallPatch(const SmBallPatch & crSource)
  : SmBSplineSurface(crSource), m_sCentre(crSource.m_sCentre), m_dRadius(crSource.m_dRadius)
{
}

SmStatus SmBallPatch::Create
  (const SmContext        & crContext,
   const SmBSplineSurface & crNet,
   const SmPoint3d        & crCentre,
   double                   dRadius,
   SmBallPatch           *& rpPatch)
{
  rpPatch = NULL;
  if ( !std::isfinite( dRadius ) || !( dRadius > 0.0 ) || !std::isfinite( crCentre.x ) || !std::isfinite( crCentre.y ) || !std::isfinite( crCentre.z ))
    { SER_MSG( SM_ERR_INVALID_INPUT, _T("SmBallPatch: a sphere's centre is finite and its radius finite and positive") ); }
  ULONG lU = 0, lV = 0;
  SmTArray<SmPoint3d> sPoints;
  SmTArray<double>    sWeights;
  SER( crNet.GetControlPointNet( lU, lV, sPoints, sWeights ));
  if ( !crNet.IsRational() || sWeights.GetSize() != sPoints.GetSize() )
    { SER_MSG( SM_ERR_INVALID_INPUT, _T("SmBallPatch: a patch on a sphere is rational") ); }
  for (ULONG ii=0; ii<sWeights.GetSize(); ii++)
    {
      if ( !( sWeights[ii] > 0.0 ))
        { SER_MSG( SM_ERR_INVALID_INPUT, _T("SmBallPatch: a patch on a sphere has positive weights") ); }
    }
  rpPatch = new (crContext) SmBallPatch( crNet, crCentre, dRadius );
  NER( rpPatch );
  return SM_SUCCESS;
}

namespace
{
double sm_Binomial(ULONG n, ULONG k)
{
  double d = 1.0;
  for (ULONG i=1; i<=k; i++) { d = d * double( n - k + i ) / double( i ); }
  return d;
}
} // namespace

SmStatus SmBallPatch::IsOnSphere
  (const SmBSplineSurface & crNet,
   const SmPoint3d        & crCentre,
   double                   dRadius,
   double                   dRelTol,
   SmBoolean              & rbOn)
{
  rbOn = FALSE;
  if ( !crNet.IsRational() || !( dRadius > 0.0 ))
    { return SM_SUCCESS; }
  // the net cut into Bezier pieces: every interior knot to the degree's multiplicity
  SmBSplineSurface sNet( crNet );
  const ULONG alDeg[2] = { sNet.GetDegree( SM_SP_U ), sNet.GetDegree( SM_SP_V ) };
  const SmSurfParamType aeDir[2] = { SM_SP_U, SM_SP_V };
  for (ULONG d=0; d<2; d++)
    {
      SmTArray<double> sKnots, sAdd;
      SmTArray<ULONG>  sMults;
      SER( sNet.GetKnots( aeDir[d], sKnots, &sMults ));
      for (ULONG k=1; k+1<sKnots.GetSize(); k++)
        {
          for (ULONG m=sMults[k]; m<alDeg[d]; m++)
            { sAdd.Add( sKnots[k] ); }
        }
      if ( sAdd.GetSize() > 0 )
        { SER( sNet.InsertKnots( aeDir[d], sAdd )); }
    }
  ULONG lU = 0, lV = 0;
  SmTArray<SmPoint3d> sPoints;
  SmTArray<double>    sWeights;
  SER( sNet.GetControlPointNet( lU, lV, sPoints, sWeights ));
  const ULONG p = alDeg[0], q = alDeg[1];
  if ( p == 0 || q == 0 || ( lU - 1 ) % p != 0 || ( lV - 1 ) % q != 0 || sWeights.GetSize() != lU * lV )
    { SER( SM_ERR ); }
  double dWMax = 0.0;
  for (ULONG k=0; k<sWeights.GetSize(); k++)
    { dWMax = smos_Max( dWMax, smos_Fabs( sWeights[k] )); }
  const double dScale = ( dRadius + crCentre.DistanceBetween( SmPoint3d( 0.0, 0.0, 0.0 ))) * dWMax;
  const double dTol = dRelTol * dScale * dScale;
  // per piece: the homogeneous coordinates' Bernstein coefficients H = w (P - C), W = w
  for (ULONG bu=0; bu+1<lU; bu+=p)
    for (ULONG bv=0; bv+1<lV; bv+=q)
      {
        // F = |H|^2 - R^2 W^2, of degree (2p, 2q)
        for (ULONG i=0; i<=2*p; i++)
          for (ULONG j=0; j<=2*q; j++)
            {
              double dF = 0.0;
              for (ULONG a=0; a<=p; a++)
                {
                  if ( a > i || i - a > p ) { continue; }
                  for (ULONG b=0; b<=q; b++)
                    {
                      if ( b > j || j - b > q ) { continue; }
                      const ULONG k1 = ( bu + a ) * lV + ( bv + b ), k2 = ( bu + i - a ) * lV + ( bv + j - b );
                      const double w1 = sWeights[k1], w2 = sWeights[k2];
                      const SmVector3d h1 = ( sPoints[k1] - crCentre ) * w1, h2 = ( sPoints[k2] - crCentre ) * w2;
                      const double dC = sm_Binomial( p, a ) * sm_Binomial( p, i - a ) * sm_Binomial( q, b ) * sm_Binomial( q, j - b );
                      dF += dC * ( h1.Dot( h2 ) - dRadius * dRadius * w1 * w2 );
                    }
                }
              dF /= sm_Binomial( 2 * p, i ) * sm_Binomial( 2 * q, j );
              if ( smos_Fabs( dF ) > dTol )
                { return SM_SUCCESS; }
            }
      }
  rbOn = TRUE;
  return SM_SUCCESS;
}

SmBoolean SmBallPatch::operator==(const SmSurface & crOther) const
{
  if ( this == &crOther )
    { return TRUE; }
  if ( !SmBSplineSurface::operator==( crOther ))
    { return FALSE; }
  const SmBallPatch * pOther = SM_CAST_NONNULL_PTR( SmBallPatch, &crOther );
  return ( pOther != NULL && m_dRadius == pOther->m_dRadius && m_sCentre.x == pOther->m_sCentre.x && m_sCentre.y == pOther->m_sCentre.y
           && m_sCentre.z == pOther->m_sCentre.z ) ? TRUE : FALSE;
}

SmStatus SmBallPatch::Copy(const SmContext & crContext, SmSurface *& rpNewSurface) const
{
  rpNewSurface = new (crContext) SmBallPatch( *this );
  NER( rpNewSurface );
  return SM_SUCCESS;
}

/*******************************************************************//**
PURPOSE: Which way the patch's normal points: +1 away from the centre, -1
    towards it, at the middle of its domain.
***********************************************************************/
SmStatus SmBallPatch::NormalSide(double & rdSide) const
{
  const SmPoint2d sMid = GetNaturalUVDomain().Evaluate( 0.5, 0.5 );
  SmPoint3d  sAt;
  SmVector3d sNormal;
  SER( SmBSplineSurface::EvaluatePoint( sMid, sAt ));
  SER( SmBSplineSurface::EvaluateNormal( sMid, TRUE, TRUE, sNormal ));
  rdSide = ( sNormal.Dot( sAt - m_sCentre ) >= 0.0 ) ? 1.0 : -1.0;
  return SM_SUCCESS;
}

SmBoolean SmBallPatch::OffsetCollapses(double dSignedOffsetDistance, double dTol) const
{
  double dSide = 1.0;
  if ( NormalSide( dSide ) != SM_SUCCESS )
    { return FALSE; }
  return ( m_dRadius + dSide * dSignedOffsetDistance <= dTol ) ? TRUE : FALSE;
}

/*******************************************************************//**
PURPOSE: The exact offset: the patch scaled about its centre onto the
    concentric sphere of radius R + s d, its weights unchanged.
***********************************************************************/
SmStatus SmBallPatch::CreateOffsetSurface
  (const SmContext & crContext,
   double            dSignedOffsetDistance,
   SmApproxTol3d     /*dThisApproxTol3d*/,
   SmSurface      *& rOffsetSurface) const
{
  rOffsetSurface = NULL;
  double dSide = 1.0;
  SER( NormalSide( dSide ));
  const double dNewRadius = m_dRadius + dSide * dSignedOffsetDistance;
  if ( !( dNewRadius > 0.0 ))
    { SER_MSG( SM_ERR_INVALID_INPUT, _T("SmBallPatch: the offset passes the sphere's centre (SmSurface::OffsetCollapses)") ); }
  const double dScale = dNewRadius / m_dRadius;
  // the net's B-spline part, its control points scaled about the centre
  SmBSplineSurface sNet( *this );
  ULONG lU = 0, lV = 0;
  SmTArray<SmPoint3d> sPoints;
  SmTArray<double>    sWeights;
  SER( GetControlPointNet( lU, lV, sPoints, sWeights ));
  if ( sPoints.GetSize() != lU * lV || sWeights.GetSize() != sPoints.GetSize() )
    { SER( SM_ERR ); }
  for (ULONG ii=0; ii<lU; ii++)
    {
      for (ULONG jj=0; jj<lV; jj++)
        {
          const ULONG k = ii * lV + jj;
          const SmPoint3d sScaled = m_sCentre + ( sPoints[k] - m_sCentre ) * dScale;
          SER( sNet.SetControlPoint( SM_CP_EUCLIDIAN_RATIONAL, ii, jj, sScaled, sWeights[k] ));
        }
    }
  SmBallPatch * pOffset = NULL;
  SER( Create( crContext, sNet, m_sCentre, dNewRadius, pOffset ));
  rOffsetSurface = pOffset;
  return SM_SUCCESS;
}

SmStatus SmBallPatch::Transform(const SmAxis2Placement & crRotateNMove, const SmVector3d * cpOptScale)
{
  double dScale = 1.0;
  if ( cpOptScale != NULL )
    {
      if ( !SM_ARE_SAME( cpOptScale->x, cpOptScale->y ) || !SM_ARE_SAME( cpOptScale->x, cpOptScale->z ) || !( cpOptScale->x > SM_EFF_ZERO ))
        { SER_MSG( SM_ERR, _T("SmBallPatch: a patch on a sphere is scaled only uniformly") ); }
      dScale = cpOptScale->x;
    }
  // the centre as the net's points move: a point curve at it, transformed alike
  SmBSplineCurve * pCentre = NULL;
  SER( SmBSplineCurve::CreatePointCurve( *GetContext(), m_sCentre, pCentre ));
  NER( pCentre );
  SmObjDelete sClean( pCentre );
  SER( pCentre->Transform( crRotateNMove, cpOptScale ));
  SER( pCentre->EvaluatePoint( pCentre->GetNaturalInterval().GetMin(), m_sCentre ));
  m_dRadius *= dScale;
  SER( SmBSplineSurface::Transform( crRotateNMove, cpOptScale ));
  return SM_SUCCESS;
}

/*******************************************************************//**
PURPOSE: Write the definition, then the net.
***********************************************************************/
SmStatus SmBallPatch::WriteToDB(SmDatabaseIO & rDB, ULONG lDBVersionNumber) const
{
  const double adValues[4] = { m_sCentre.x, m_sCentre.y, m_sCentre.z, m_dRadius };
  if ( rDB.GetFileType() == SM_ASCII )
    {
      std::ostream & rOut = *rDB.GetOutStreamPtr();
      for (int i=0; i<4; i++) { rOut << adValues[i] << " "; }
      rOut << " SmBallPatch centre, radius\n";
    }
  else
    { SER( rDB.WriteDoubles( adValues, 4 )); }
  SER( SmBSplineSurface::WriteToDB( rDB, lDBVersionNumber ));
  return SM_SUCCESS;
}

SmStatus SmBallPatch::ReadFromDB
  (SM_TYPE           /*lType*/,
   SmDatabaseIO    & rDB,
   const SmContext & crContext,
   SmSurface      *& rpNewSurface,
   ULONG             lDBVersionNumber)
{
  if ( rpNewSurface != NULL )
    { SER_MSG( SM_ERR, _T("SmBallPatch: read into a new surface only") ); }
  double adValues[4] = { 0.0, 0.0, 0.0, 0.0 };
  if ( rDB.GetFileType() == SM_ASCII )
    {
      std::istream & rIn = *rDB.GetInStreamPtr();
      for (int i=0; i<4; i++) { rIn >> adValues[i]; }
      rDB.GoToNextLine();
      if ( !rIn )
        { SER_MSG( SM_ERR, _T("SmBallPatch: its definition could not be read") ); }
    }
  else
    { SER( rDB.ReadDoubles( adValues, 4 )); }
  SmSurface * pNet = NULL;
  SER( SmBSplineSurface::ReadFromDB( SmBSplineSurface_TYPE, rDB, crContext, pNet, lDBVersionNumber ));
  SmObjDelete sDeleteNet( pNet );
  const SmBSplineSurface * pNetBS = SM_CAST_PTR( SmBSplineSurface, pNet );
  NER( pNetBS );
  SmBallPatch * pPatch = NULL;
  SER( Create( crContext, *pNetBS, SmPoint3d( adValues[0], adValues[1], adValues[2] ), adValues[3], pPatch ));
  rpNewSurface = pPatch;
  return SM_SUCCESS;
}

SmBoolean SmBallPatch::IsKindOf(SM_TYPE t) const
{
  return ( ( SmBallPatch_TYPE == t ) ? TRUE : SmBSplineSurface::IsKindOf( t ));
}

void SmBallPatch::Dump(void) const
{
  TCHAR sBuff[SM_TBLOCK_SIZE];
  smos_sprintf( sBuff, _T("\nSmBallPatch: centre (%g, %g, %g), radius %g"), m_sCentre.x, m_sCentre.y, m_sCentre.z, m_dRadius );
  smos_WriteBuffer( sBuff );
  SmBSplineSurface::Dump();
}
