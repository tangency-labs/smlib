// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

/**********************************************************************//**
* FILE NAME --- SmFilletExact.h
* PURPOSE: SolidCraft change. Which constant rounds are exact: the faces of a
*   blend as the natural quadrics they lie on (plane, sphere, circular
*   cylinder, circular cone), and the one classification of the curve the
*   rolling ball's centre runs on (the spine, where the faces' offsets by the
*   radius meet) from their kinds and placement.
**********************************************************************/

#ifndef __SMFILLETEXACT_H__
#define __SMFILLETEXACT_H__

#ifndef __SMOS_TYPES_H__
#include <SmTypes.h>
#endif

#ifndef __SMVECTOR3D_H__
#include <SmVector3d.h>
#endif

class SmSurface;
class SmVector2d;

/*******************************************************************//**
PURPOSE: A surface as the natural quadric it lies on, with the surface's own
   normal: a plane, a sphere, a circular cylinder, a circular cone, or none.
   An offset of a cylinder or a sphere by its radius collapses onto its axis
   (a line) or its centre (a point).

NOTES:
  m_sPoint : plane, a point on it; sphere, point: the centre; cylinder, line:
             a point of the axis; cone: the apex
  m_sDir   : plane: the surface's unit normal; cylinder, line: the unit axis;
             cone: the unit axis from the apex into the nappe the surface is on
  m_dRadius: sphere, cylinder: the radius
  m_dAngle : cone: the half-angle, in radians, in (0, pi/2)
  m_dSide  : sphere, cylinder, cone: +1 when the surface's normal points away
             from the centre or the axis, -1 when toward it
***********************************************************************/
class SM_EXPORT SmBlendQuadric
{
public:
  enum Kind { SM_BQ_NONE, SM_BQ_PLANE, SM_BQ_SPHERE, SM_BQ_CYLINDER, SM_BQ_CONE, SM_BQ_LINE, SM_BQ_POINT };

  Kind       m_eKind   = SM_BQ_NONE;
  SmPoint3d  m_sPoint;
  SmVector3d m_sDir;
  double     m_dRadius = 0.0;
  double     m_dAngle  = 0.0;
  double     m_dSide   = 1.0;

  // The natural quadric pSurface lies on (SM_BQ_NONE if it lies on none), oriented by its normal
  // at crUV; an SmOffsetSurface is its base's quadric moved by the offset. A surface whose kind is
  // a natural quadric but which lies farther than dTol from the quadric read off it at crUV is an
  // error, not a surface of no kind.
  static SmStatus FromSurface(const SmSurface * pSurface, const SmPoint2d & crUV, double dTol, SmBlendQuadric & rQuadric);

  // This quadric moved dDist along its normal (a cylinder or sphere moved in by its radius, within
  // dTol, is its axis or its centre).
  SmBlendQuadric Offset(double dDist, double dTol) const;

  // The unit normal at crPoint (on the quadric), as the surface's.
  SmVector3d Normal(const SmPoint3d & crPoint) const;

  // The distance of crPoint from the quadric (from the cone's nappe, measured in the plane
  // through the axis).
  double Distance(const SmPoint3d & crPoint) const;

  // The implicit form x.A.x + 2 b.x + c = 0 (sphere, cylinder, cone: the double cone; plane:
  // A = 0, b the half normal).
  void Implicit(double aA[3][3], double aB[3], double & rC) const;

  SmBoolean IsRevolved() const { return m_eKind == SM_BQ_SPHERE || m_eKind == SM_BQ_CYLINDER || m_eKind == SM_BQ_CONE || m_eKind == SM_BQ_LINE; }
};

/*******************************************************************//**
PURPOSE: The curve a constant round's rolling ball's centre runs on, by the
   faces' kinds and placement. All but SM_BS_TRACED are rational, so the round
   is an exact rational surface.

NOTES:
  SM_BS_LINE         : a line (two planes; a plane, or a cylinder, parallel to
                       a cylinder's axis): the round is a cylinder
  SM_BS_CIRCLE       : a circle about a common axis (coaxial quadrics; a plane
                       square to the axis; a sphere with its centre on it):
                       the round is the ball's section turned about the axis
  SM_BS_POINT        : one point (two cylinders as wide as the ball, their
                       axes meeting): the round is a piece of the ball
  SM_BS_CONIC        : a conic, a plane section of the other quadric (a
                       plane and a cylinder or cone at any other angle)
  SM_BS_CONIC_PAIR   : one of two plane conics (two cylinders or cones whose
                       offsets touch a common ball the same way round)
  SM_BS_NODAL_QUARTIC: a quartic with a node where the offsets touch (two
                       cylinders or cones, or one and a sphere, touching at
                       one point the same way round)
  SM_BS_TRACED       : none of these: the spine is not rational, and the
                       round is traced and fitted
***********************************************************************/
enum SmBlendSpine
{
  SM_BS_TRACED,
  SM_BS_LINE,
  SM_BS_CIRCLE,
  SM_BS_POINT,
  SM_BS_CONIC,
  SM_BS_CONIC_PAIR,
  SM_BS_NODAL_QUARTIC
};

/*******************************************************************//**
PURPOSE: A blend's exact class: its faces' quadrics, their offsets (the ball's
   centre lies on both), and the spine's data.

NOTES:
  m_adContact[i]: the ball touches face i at its centre plus the radius times
                  m_adContact[i] times the face's normal there (+1 or -1)
  m_sPoint, m_sDir:
    SM_BS_CIRCLE       : a point of the common axis, its unit direction
    SM_BS_CONIC        : a point of the conic's plane, its unit normal
    SM_BS_CONIC_PAIR   : the same, for the conic the spine is
    SM_BS_NODAL_QUARTIC: the node (m_sDir unused)
***********************************************************************/
struct SM_EXPORT SmExactBlend
{
  SmBlendSpine   m_eSpine = SM_BS_TRACED;
  SmBlendQuadric m_aFaces[2];
  SmBlendQuadric m_aOffsets[2];
  double         m_adContact[2] = { 1.0, 1.0 };
  double         m_dRadius = 0.0;
  SmPoint3d      m_sPoint;
  SmVector3d     m_sDir;
};

// The one classification: the spine of a constant round of radius rBlend.m_dRadius between faces
// on rBlend.m_aFaces, given their offsets rBlend.m_aOffsets and crCentre, a point of the spine
// (the ball's centre at one end), from the faces' kinds and placement. Placements (axes meeting,
// a common ball, a touching point) are decided within dTol.
SM_EXPORT SmStatus sm_ClassifyBlendSpine(const SmPoint3d & crCentre, double dTol, SmExactBlend & rBlend);

#endif // __SMFILLETEXACT_H__
