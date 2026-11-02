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
* order), takes its first k + 1 rows from the round's homogeneous derivatives
* at its end section, taken of its geometry and not of its parameters (its
* coordinates over its rails' mean weight, by its rails' mean arc length along
* it: two parameters of one round give one fill) and its last k + 1 from F's derivatives
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
* A G1 fill (order 1) is the same layout as a graph over the plane square
* to the corner's axis: its interior map's projection is the base points'
* (each the middle of its face's spring curve, the parabola tangent to the
* two rails at their ends), the setback patches' projections are fixed by
* the rounds and the spring curves, and only heights along the axis are
* free (SmBlendSetback.cpp, sm_SbGraphG1).
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
  // A face at the corner: a plane through m_sPoint, m_sNormal its unit normal out of the material;
  // or a natural quadric (m_eKind): a cylinder or cone about the axis through m_sPoint along the unit
  // m_sAxis, of radius m_dRadius there and m_dRadius + m_dSlope h at h along the axis (a cylinder's
  // slope 0), or a sphere about m_sPoint of radius m_dRadius; m_dSide +1 where the material is inside
  // it, -1 outside. (m_sNormal is a curved face's normal out of the material at the corner's vertex.)
  enum Kind { SM_SB_PLANE, SM_SB_CYLINDER, SM_SB_CONE, SM_SB_SPHERE };
  struct Face
  {
    SmPoint3d  m_sPoint;
    SmVector3d m_sNormal;
    Kind       m_eKind   = SM_SB_PLANE;
    SmVector3d m_sAxis   = SmVector3d(0.0, 0.0, 1.0);
    double     m_dRadius = 0.0;
    double     m_dSlope  = 0.0;
    double     m_dSide   = 1.0;
  };

  // A face's signed distance at crX (negative in the material; a cone's in the plane through its
  // axis) and, in pNormal, its unit normal out of the material there.
  static double Distance(const Face & crFace, const SmPoint3d & crX, SmVector3d * pNormal = NULL);

  // The spring curve on a face from P to Q, tangent there to the unit directions crTP (leaving P) and
  // crTQ (arriving at Q): in the face's rational parameter plane (a plane's own coordinates; a
  // cylinder's or cone's tan(theta / 2) and its axial coordinate; a sphere's stereographic plane from
  // the point opposite the corner) the parabola through the images of P and Q tangent to the images
  // of their directions, mapped onto the face: rational, exactly on it, its degree 2 on a plane, 4 on
  // a cylinder (the axial coordinate a quartic with the parabola's end jets, its middle coefficient
  // the mean of its neighbours) or a sphere, 6 on a cone. Returned as its two halves meeting at its
  // middle (the base point), each in homogeneous Bernstein coefficients (x w, y w, z w, w), from its
  // rail end to the base point (the fill matches them to its rounds and to its interior map).
  // Refused (SM_ERR_INVALID_INPUT, ppReason) where the tangent lines' images cross behind P or Q (no
  // corner region between the rails) or a weight would not stay positive.
  // With each half, two homogeneous polynomials along it (Bernstein, (x, y, z, w)) that with the half
  // span the face's tangent plane there projectively: the face's parametric partials along the
  // parameter-plane parabola (a plane's two directions, degree 0; a cylinder's circumferential
  // direction, 4, and its axis, 0; a cone's circumferential direction and generator, 4 and 4; a
  // sphere's stereographic partials, 2 and 2). A column next to the half is G1 to the face exactly
  // where it differs from the half by alpha T1 + beta T2 + gamma H0 (Hahn's connecting functions).
  static SmStatus SpringCurve(const Face & crFace, const SmPoint3d & crP, const SmVector3d & crTP, const SmPoint3d & crQ, const SmVector3d & crTQ,
                              std::vector<double> & rFromP, std::vector<double> & rFromQ, const TCHAR ** ppReason = NULL,
                              std::vector<double> * pTanP = NULL, std::vector<double> * pTanQ = NULL);
  // (pTanP / pTanQ: T1's coefficients then T2's, 4 (k1 + 1) + 4 (k2 + 1) values, k1 and k2 as above.)

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
  // edge i. pvLeast (optional) gets the least setbacks, ppReason (optional) a refusal's reason.
  // Refused (SM_ERR_INVALID_INPUT): an edge
  // whose faces are not two of the faces given, a face not shared by exactly two of the edges,
  // edges along one line on a face, a radius not positive.
  static SmStatus Distances(const std::vector<Face> & crFaces, const std::vector<Edge> & crEdges, std::vector<double> & rvSetbacks,
                            std::vector<double> * pvLeast = NULL, const TCHAR ** ppReason = NULL);

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
  // m_lDegree: the interior map's degree. m_pReason: why FillHole refused, NULL when it did not.
  // m_lSolves: the quadratic programmes a G1 fill solved (its count, not its time).
  // m_bUncertified: FillHole refused a G1 fill as no fill up to degree 9 is certified convex.
  // The patches are owned by the fill.
  struct Fill
  {
    std::vector<SmBSplineSurface *> m_vPatches;
    std::vector<ULONG>              m_vRound;
    std::vector<ULONG>              m_vFace;
    std::vector<SmBoolean>          m_vReversed;
    ULONG                           m_lDegree = 0;
    const TCHAR *                   m_pReason = NULL;
    ULONG                           m_lSolves = 0;
    bool                            m_bUncertified = false;

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
  // interior map no degree up to 9 can meet, or a fill that folds. A G1 fill is taken only
  // certified convex (at a convex corner; the void's at a concave one): SM_ERR, its reason and
  // m_bUncertified, where none up to degree 9 is (its setbacks too short for it: CertifiedScale).
  static SmStatus FillHole(const SmContext & crContext, const std::vector<Face> & crFaces, const std::vector<Round> & crRounds, ULONG lOrder,
                           double dTol, Fill & rFill);

  // CertifiedScale's probe of a G1 fill at degree 9: whether its conditions can be met
  // (m_bFeasible), and the least uniform violation of its normalized sector conditions over the
  // heights (m_dMargin: <= 0 where met, m_bMargin asks it), one quadratic programme each; or
  // (m_bFill) the fill as FillHole's, its degrees tried from the least up only.
  struct Probe
  {
    bool   m_bFill = false;   // the fill itself, from the least degree up (no degree-9 proof first)
    bool   m_bMargin = false;
    bool   m_bFeasible = false;
    double m_dMargin = 0.0;
  };
  static SmStatus FillHoleProbe(const SmContext & crContext, const std::vector<Face> & crFaces, const std::vector<Round> & crRounds, double dTol,
                                Fill & rFill, Probe & rProbe);

  // The least scale k of the setbacks crSetBacks (one per edge, the edges and faces as for
  // Distances, every face a plane through the corner's vertex, every edge a straight line there),
  // each held at its most (pvMost, optional), at which the G1 fill of the corner's constant
  // circular rounds (radius m_dRadius between the edge's two faces, each round's end section
  // square to its edge at its setback) is certified convex (FillHole), to 1e-2: 1 where the
  // setbacks themselves are; else the certifying end of a bracket 1e-2 wide, found from the
  // margin by which the fill's conditions fail (FillHoleProbe: at 1, 1.02, their secant's root,
  // then the parabola's), up to 2; 0, with its reason (ppReason), where none up to 2 is. plSolves
  // (optional): the quadratic programmes solved. Refused (SM_ERR_INVALID_INPUT, its reason) where
  // the corner is not one of planes and straight edges.
  static SmStatus CertifiedScale(const SmContext & crContext, const std::vector<Face> & crFaces, const std::vector<Edge> & crEdges,
                                 const std::vector<double> & crSetBacks, double & rdScale, ULONG * plSolves = NULL, const TCHAR ** ppReason = NULL,
                                 const std::vector<double> * pvMost = NULL);

  // The most any patch of crFill stands past a face out of the corner (positive: past it): out
  // of the material for a convex corner (bConvex), into it for a concave one. A certified bound:
  // every patch lies in the convex hull of its control points (weights positive). Beside a
  // quadric face, certified by the face's implicit form over each patch (Bernstein, subdivision;
  // divided by v^2 or (1 - v)^2 along the patch's column on that face): rdWorst HUGE_VAL where
  // that is not certified.
  static SmStatus Containment(const Fill & crFill, const std::vector<Face> & crFaces, SmBoolean bConvex, double & rdWorst);

  // Owner's ruling (a) on one setback patch (a single Bezier span, as FillHole makes them) G1 to
  // crA along its v = 0 and to crB along v = 1: whether II >= dK I (dK < 0) over the whole patch,
  // II taken with the unit normal on crInto's side (into the fill's convex side), certified by
  // the faces' normal field (exact along the columns): riResult 1 certified, 0 refused (a value
  // below), -1 undecided. Refused (SM_ERR_INVALID_INPUT): not one span, dK not negative.
  // The work of the G1 fills beside curved faces since the last reset, process-wide (counts for
  // tests, never timings): raOut[0] the fills made (each a FillHole or a probe of one), [1] the
  // degrees tried (each one quadratic programme), [2] the sector rows built for them. bReset sets
  // them to 0 after reading.
  static void CurvedCounts(unsigned long long (&raOut)[3], bool bReset = false);

  // A constant round's least normal curvature along its end section (its m_dEnd), II with the unit
  // normal toward crInto, to rounding (the rounds' part of ruling (a)'s k: sm_SbRoundRowK); 0 where
  // the section is nowhere concave. Refused where the section is not a circular arc.
  static SmStatus RoundEndCurvature(const SmContext & crContext, const Round & crRound, const SmVector3d & crInto, double & rdK,
                                    const TCHAR ** ppReason = NULL);

  static SmStatus NoMoreConcave(const SmBSplineSurface & crPatch, const Face & crA, const Face & crB, const SmVector3d & crInto, double dK,
                                int & riResult);
};

#endif // __SMBLENDSETBACK_H__
