// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

/**********************************************************************//**
* FILE NAME --- SmBlendSetback.h
* PURPOSE: SolidCraft change. The setback corner: where rounds meet at a
*   vertex and stop short of it, each at its setback, the hole they leave is
*   filled by patches that meet the rounds and the faces to the order asked
*   (tangent, G1, or curvature continuous, G2).
*
* The pieces, each usable on its own:
*   - Distances: each round's setback, by the owner's rule (the larger of
*     1.5 r and the least setback, where the rails of the rounds beside it
*     cross on their common face, plus 0.25 r);
*   - EaseSurface: a round's sections eased from its own into another
*     section (the curvature continuous one a G2 fill needs) over a stretch
*     of its rails;
*   - FillHole: the hole's patches, by the setback split: one setback patch
*     per round, between its end section, the two faces and an inner curve,
*     and an n-sided interior region within the inner curves, split into n
*     quadrilaterals at its centre.
*
* The interior region is one polynomial map F of a regular n-gon (the
* domain) into space, and each of its quadrilaterals F composed with the
* bilinear map of a quarter of the n-gon onto the unit square: the interior
* parts meet each other as one surface (no conditions between them, at their
* centre neither). A setback patch, of degree 2k + 1 across the hole (k the
* order), takes its first k + 1 rows from the round's own homogeneous
* derivatives at its end section and its last k + 1 from F's derivatives
* across the n-gon's edge, F composed with the map (u, v) -> e(v) + (1 - u)
* d(v) (e the edge, d a field of directions out of the n-gon along it, the
* two fields at an n-gon vertex opposite): it meets the round to order k in
* homogeneous coordinates and the interior parts to order k by that map, both
* exactly. It meets the faces to order k where its first k + 1 columns lie on
* the face's plane: on the round's side the round's do (it meets the face to
* that order), on the interior's side F is held to it, linear conditions on
* F's coefficients. F is otherwise the one of least bending energy (the sum
* over the patches of the integral of |R_uu|^2 + 2 |R_uv|^2 + |R_vv|^2), at
* the least degree whose conditions can be met: the conditions are met
* exactly, the energy only chooses among the maps that meet them.
*
* Faces are planes (a curved face is not filled here yet). Every patch lies
* in the convex hull of its control points (all weights positive), so
* Containment's bound is certified.
**********************************************************************/

#ifndef __SMBLENDSETBACK_H__
#define __SMBLENDSETBACK_H__

#ifndef __SMOS_TYPES_H__
#include <SmTypes.h>
#endif

#ifndef __SMVECTOR3D_H__
#include <SmVector3d.h>
#endif

#include <vector>

class SmBSplineSurface;
class SmContext;

// A stretch of a fillet surface's rail parameter over which its sections ease: w = 0 at
// m_dFrom (the round's own section), w = 1 at m_dTo (the one eased into), by the C2 law
// w(s) = 10 s^3 - 15 s^4 + 6 s^5 (its Bernstein coefficients 0, 0, 0, 1, 1, 1), which is flat
// to second order at both ends; w stays 1 beyond m_dTo.
struct SmEaseRamp
{
  double m_dFrom;
  double m_dTo;
};

class SM_EXPORT SmBlendSetback
{
public:
  // A face at the corner: a plane through m_sPoint, m_sNormal its unit normal out of the material.
  struct Face
  {
    SmPoint3d  m_sPoint;
    SmVector3d m_sNormal;
  };

  // A rounded edge at the corner, for the setbacks: m_sAway its unit direction away from the
  // vertex, m_dRadius the round's radius there, m_alFaces its two faces (indices into the faces).
  struct Edge
  {
    SmVector3d m_sAway;
    double     m_dRadius;
    ULONG      m_alFaces[2];
  };

  // Each edge's setback, the owner's rule: the larger of 1.5 r and its least setback plus
  // 0.25 r. Its least setback is the farthest along it, over its two faces, that its rail
  // crosses the rail of the other edge on that face: a round between planes whose normals are
  // psi apart has its rails r tan(psi / 2) from its edge, and on a face where edges i, j meet at
  // the angle alpha, rails d_i, d_j from them cross (d_j + d_i cos alpha) / sin alpha along
  // edge i. pvLeast (optional) gets the least setbacks. Refused (SM_ERR_INVALID_INPUT): an edge
  // whose faces are not two of the faces given, a face not shared by exactly two of the edges,
  // edges along one line on a face, a radius not positive.
  static SmStatus Distances(const std::vector<Face> & crFaces, const std::vector<Edge> & crEdges, std::vector<double> & rvSetbacks,
                            std::vector<double> * pvLeast = NULL);

  // The surface whose sections ease from crFrom's to crTo's over the ramps along eAlong, the
  // rails' direction, in which both run over the same interval with the same rail parameter:
  // P = (1 - w) From + w To in homogeneous coordinates, w the C2 law (SmEaseRamp), so P is From
  // to second order where w is 0 and To where it is 1. Exact: both surfaces are raised to
  // common degrees and knots, put in Bezier form along the rails, and each span's net multiplied
  // by the law's Bernstein form.
  static SmStatus EaseSurface(const SmContext & crContext, const SmBSplineSurface & crFrom, const SmBSplineSurface & crTo,
                              SmSurfParamType eAlong, const std::vector<SmEaseRamp> & crRamps, SmBSplineSurface *& rpEased);

  // A round at the corner, ending at its setback: m_pSurface its surface (one Bezier span
  // across its rails, any along them), m_eAlong its rails' direction, m_dEnd the rail
  // parameter of its end section, m_bCornerAbove whether the corner lies at rail parameters
  // above m_dEnd (the round below). The round must meet its faces to the fill's order at its
  // end section (a G1 round eased first for a G2 fill).
  struct Round
  {
    const SmBSplineSurface * m_pSurface;
    SmSurfParamType          m_eAlong;
    double                   m_dEnd;
    SmBoolean                m_bCornerAbove;
  };

  // The patches of a filled hole, in the cycle's order: rounds m = 0 .. n - 1 round the vertex,
  // round m between faces m and m + 1 (mod n).
  //
  // m_vPatches[m], m < n, round m's setback patch, (u, v) in [0, 1]^2:
  //   u = 0: round m's end section, v running from face m to face m + 1 (m_vReversed[m]: its
  //          across parameter runs the other way);
  //   v = 0: on face m (a spring curve's half), v = 1: on face m + 1;
  //   u = 1: the inner curve from face m's base point (v = 0) to face m + 1's (v = 1), its
  //          half v <= 1/2 interior part m's side s = 1 (t = 1 - 2 v), its half v >= 1/2
  //          interior part m + 1's side t = 1 (s = 2 v - 1).
  // m_vPatches[n + q], interior part q, (s, t) in [0, 1]^2: (0, 0) the centre, (1, 1) face q's
  //   base point; t = 0 the spoke to the middle of inner curve q, which interior part q + 1 has
  //   as its side s = 0 (its t the same as part q's s).
  // m_vRound[m]: the caller's index of the cycle's round m; m_vFace[m]: of face m.
  // m_lDegree: the interior map's degree. The patches are owned by the fill.
  struct Fill
  {
    std::vector<SmBSplineSurface *> m_vPatches;
    std::vector<ULONG>              m_vRound;
    std::vector<ULONG>              m_vFace;
    std::vector<SmBoolean>          m_vReversed;
    ULONG                           m_lDegree = 0;

    Fill() {}
    ~Fill();
    Fill(const Fill &) = delete;
    Fill & operator=(const Fill &) = delete;
  };

  // The hole the rounds leave at a corner of the faces, filled to lOrder (1: G1, 2: G2) by the
  // setback split (file header). dTol: the length the faces and the rounds' ends are matched to
  // (the rounds' rails on their faces, their end rows' first columns on the faces' planes). The
  // fill is not held to the faces' half-spaces here: Containment bounds that. Refused (SM_ERR_INVALID_INPUT, with its reason): rounds that do not close round
  // the vertex face to face, a round with more than one span across, a round that does not meet
  // its faces to lOrder at its end, an order other than 1 or 2. SM_ERR: conditions on the
  // interior map no degree up to 9 can meet, or a fill that folds.
  static SmStatus FillHole(const SmContext & crContext, const std::vector<Face> & crFaces, const std::vector<Round> & crRounds, ULONG lOrder,
                           double dTol, Fill & rFill);

  // The most any patch of crFill stands past a face out of the corner (positive: past it): out
  // of the material for a convex corner (bConvex), into it for a concave one. A certified bound:
  // every patch lies in the convex hull of its control points (weights positive).
  static SmStatus Containment(const Fill & crFill, const std::vector<Face> & crFaces, SmBoolean bConvex, double & rdWorst);
};

#endif // __SMBLENDSETBACK_H__
