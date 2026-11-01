// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

/**********************************************************************//**
* FILE NAME --- SmBallPatch.h
* PURPOSE: SolidCraft change. A rational patch on a sphere, stored with its
*   definition: the sphere's centre and radius. The blend corner's Ball
*   makes it (SmBlendCornerBall.cpp); its offset is exact.
**********************************************************************/

#ifndef __SMBALLPATCH_H__
#define __SMBALLPATCH_H__

#ifndef __SMBSPLINESURFACE_H__
#include <SmBSplineSurface.h>
#endif

/*******************************************************************//**
PURPOSE: A rational B-spline patch every point of which is on the sphere
    of centre C and radius R: its net is exact, not an approximation, and
    (C, R) is its definition, recorded where it is made.

NOTES: Its offset by d along its normal is the same patch scaled about C by
    (R + s d) / R (s = +1 where its normal points away from C, -1 where
    towards it), its weights unchanged: exactly on the concentric sphere of
    radius R + s d. The offset collapses when R + s d is not positive.
    A rigid or uniformly scaled transform moves (C, R) with the net.
***********************************************************************/
class SM_EXPORT SmBallPatch : public SmBSplineSurface
{
protected:
  SmPoint3d m_sCentre;
  double    m_dRadius;

  SmBallPatch();
  SmBallPatch(const SmBSplineSurface & crNet, const SmPoint3d & crCentre, double dRadius);

  // +1 where the normal points away from the centre, -1 where towards it
  SmStatus NormalSide(double & rdSide) const;

public:
  SmBallPatch(const SmBallPatch & crSource);
  virtual ~SmBallPatch() {}

  // The patch crNet (copied) on the sphere of centre crCentre and radius dRadius. Refused
  // (SM_ERR_INVALID_INPUT, named): a radius not finite and positive, a net not rational with
  // positive weights.
  static SmStatus Create(const SmContext        & crContext,
                         const SmBSplineSurface & crNet,
                         const SmPoint3d        & crCentre,
                         double                   dRadius,
                         SmBallPatch           *& rpPatch);

  const SmPoint3d & GetCentre() const { return m_sCentre; }
  double            GetRadius() const { return m_dRadius; }

  // SmSurface
  virtual SmBoolean operator==(const SmSurface & crOther) const;
  virtual SmStatus  Copy(const SmContext & crContext, SmSurface *& rpNewSurface) const;
  virtual SmStatus  CreateOffsetSurface(const SmContext & crContext, double dSignedOffsetDistance, SmApproxTol3d dThisApproxTol3d,
                                        SmSurface *& rOffsetSurface) const;
  virtual SmBoolean OffsetCollapses(double dSignedOffsetDistance, double dTol) const;
  virtual SmStatus  Transform(const SmAxis2Placement & crRotateNMove, const SmVector3d * cpOptScale = NULL);
  virtual SmStatus  WriteToDB(SmDatabaseIO & rDB, ULONG lDBVersionNumber) const;
  static  SmStatus  ReadFromDB(SM_TYPE lType, SmDatabaseIO & rDB, const SmContext & crContext, SmSurface *& rpNewSurface,
                               ULONG lDBVersionNumber);

  SM_COMMON(SmBallPatch, SmBSplineSurface, SmBallPatch_TYPE);
  virtual void Dump(SmBoolean bAbbrev) const { SmBSplineSurface::Dump(bAbbrev); }
  virtual void Dump(ULONG lLevel) const      { SmBSplineSurface::Dump(lLevel); }
  virtual void Dump(const TCHAR * message) const { SmBSplineSurface::Dump(message); }

} ; // end class SmBallPatch

#endif // __SMBALLPATCH_H__
