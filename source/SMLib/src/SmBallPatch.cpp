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

SmBoolean SmBallPatch::operator==(const SmSurface & crOther) const
{
  if ( this == &crOther )
    { return TRUE; }
  if ( !SmBSplineSurface::operator==( crOther ))
    { return FALSE; }
  const SmBallPatch * pOther = SM_CAST_PTR( SmBallPatch, &crOther );
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
