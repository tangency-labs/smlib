// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

/*******************************************************************//**
* FILE NAME --- SmProceduralBlend.h
* PURPOSE: SolidCraft change. A blend surface stored as its definition: the
*   two surfaces it rolls on, the ball's radius, the section rule and the
*   leader that fixes its parameter. It evaluates exactly from that
*   definition; its net is the certified B-spline approximation of it.
**********************************************************************/

#ifndef __SMPROCEDURALBLEND_H__
#define __SMPROCEDURALBLEND_H__

#ifndef __SMBSPLINESURFACE_H__
#include <SmBSplineSurface.h>
#endif

class SmBSplineCurve;

// the most definitions a blend's supports may nest (a round on a round on a round)
#define SM_BLEND_MAX_DEPTH 3

/*******************************************************************//**
PURPOSE: The counts of a blend's exact evaluations, the work its definition
    does: every evaluation of a point (with its partials) of any procedural
    blend, the station solves it takes, and the Newton steps of those
    solves. Machine independent, unlike timings.

NOTES: Shared by every thread and every blend; reset with
    SmProceduralBlend::ResetCounts.
***********************************************************************/
struct SmBlendCounts
{
  unsigned long long lEvaluations;     // points evaluated on a definition
  unsigned long long lSolves;          // station solves (the spine point and the feet); an evaluation at a
                                       // station just solved (the same s, another t) takes it from a memo
  unsigned long long lNewtonSteps;     // Newton steps of those solves
  unsigned long long lMostNewtonSteps; // the most one solve took
};

/*******************************************************************//**
PURPOSE: SolidCraft change. A blend of order n between two surfaces, by its
    definition.

    The definition: the supports S1, S2 (owned copies; a support may itself
    be a procedural blend, nested at most SM_BLEND_MAX_DEPTH deep), the ball's
    side of each (+-1 along its normal), the ball's radius (constant, or a
    law r(s) given as the x of a curve over the leader's parameter), the
    section's order n (1 the ball's circle, 2 to 4 the G^n section with
    thumbweight w), the leader L(s) (a curve near the spine whose planes
    square to it fix the parameter s) and seeds for the feet: the uv curves
    (u1, v1)(s), (u2, v2)(s) on the supports. Its net B is within dBound of
    it (the certified bound).

    Evaluation at (s, t): the feet x = (u1, v1, u2, v2) solved by Newton from
    the seeds on S1 + sigma1 r n1 = S2 + sigma2 r n2, (S1 + sigma1 r n1 - L(s))
    . L'(s) = 0; then, n = 1, the arc of the circle the balls' envelope
    touches the ball along (the great circle where r is constant), from the
    first rail (t = 0) to the second (t = 1), as the rational quadratic of
    weights 1, cos(theta / 2), 1; n >= 2, the G^n section of the faces' order-n
    jets along the rails (SmFilletSection.cpp). Partials in s from the
    implicit function theorem, order by order (truncated power series).

    The surface's (u, v) map to (s, t) affinely, and swapped if SwapUV was
    called; the map follows Reverse, Reparameterize and SwapUV, so the net
    and the definition keep the same parametrization. Evaluate,
    EvaluateSimple and EvaluatePoint are the definition's; the net (the
    SmBSplineSurface's own evaluators, called explicitly) is what bounds,
    seeds and tessellation use.

NOTES: The order-1 definition is exact for a constant radius, and for a
    radius law the touch circle square to the spine; n >= 2 needs a constant
    radius, as the G^n section does.
***********************************************************************/
class SM_EXPORT SmProceduralBlend : public SmBSplineSurface
{
protected:
  SmSurface      * m_apSupport[2];  // owned
  double           m_adSide[2];     // the ball's side of each support, +-1 along its normal
  double           m_adR[2];        // the ball's offsets from S1 and S2 (its radius, both), when m_pLaw is NULL
  SmBSplineCurve * m_pLaw;          // owned, or NULL: the radius r(s) as the x of a curve over the leader's parameter
  ULONG            m_lOrder;        // the section's order n, 1 to SM_GN_MAX_ORDER
  double           m_dThumbweight;  // n >= 2: the G^n section's thumbweight w < 1
  SmBSplineCurve * m_pLeader;       // owned: L(s)
  SmBSplineCurve * m_apSeed[2];     // owned: the feet's uv on each support, over s
  double           m_dBound;        // the net's certified distance from the definition, > 0
  // (s, t) = (m_adA[0] + m_adB[0] p, m_adA[1] + m_adB[1] q), (p, q) = (u, v), or (v, u) when m_bSwapped
  double           m_adA[2];
  double           m_adB[2];
  SmBoolean        m_bSwapped;
  // this definition's identity for the station memo (SmFilletSection.cpp): new whenever the
  // definition is made or moved, never reused
  unsigned long long m_lSerial;

  SmProceduralBlend();
  SmProceduralBlend(const SmBSplineSurface & crNet);

  void     Clear();
  static void CountEvaluation();
  static void CountSolve(ULONG lNewtonSteps);
  void        Renew();
  SmStatus CopyDefinition(const SmProceduralBlend & crSource);
  // u = dAlpha + dBeta u' (bU) or the same in v: the map's terms after the net's reparametrization
  void     ReparametrizeMap(SmBoolean bU, double dAlpha, double dBeta);

public:
  SmProceduralBlend(const SmProceduralBlend & crSource);
  virtual ~SmProceduralBlend();

  // A blend by its definition on the net crNet (both copied, nothing owned from the caller).
  // cpLaw NULL: the ball's offsets crOffset from S1 and S2 (both its radius; the G^n section's
  // system takes them apart, the circle's must agree). Refused (SM_ERR_INVALID_INPUT, named): an
  // order out of 1 to SM_GN_MAX_ORDER, n >= 2 with a law or a thumbweight not below 1, sides other
  // than +-1, offsets or a bound not finite and positive, a circle's offsets that differ, supports
  // nested deeper than SM_BLEND_MAX_DEPTH.
  static SmStatus Create(const SmContext        & crContext,
                         const SmBSplineSurface & crNet,
                         const SmSurface        * apSupport[2],
                         const double           (&crSide)[2],
                         const double           (&crOffset)[2],
                         const SmBSplineCurve   * cpLaw,
                         ULONG                    lOrder,
                         double                   dThumbweight,
                         const SmBSplineCurve   & crLeader,
                         const SmBSplineCurve   * apSeed[2],
                         double                   dBound,
                         SmProceduralBlend     *& rpBlend);

  // the definition
  const SmSurface      * GetSupport(ULONG i) const { return i < 2 ? m_apSupport[i] : NULL; }
  double                 GetSide(ULONG i) const    { return i < 2 ? m_adSide[i] : 0.0; }
  double                 GetOffset(ULONG i) const  { return i < 2 ? m_adR[i] : 0.0; }
  const SmBSplineCurve * GetLaw() const            { return m_pLaw; }
  ULONG                  GetOrder() const          { return m_lOrder; }
  double                 GetThumbweight() const    { return m_dThumbweight; }
  const SmBSplineCurve * GetLeader() const         { return m_pLeader; }
  const SmBSplineCurve * GetSeed(ULONG i) const    { return i < 2 ? m_apSeed[i] : NULL; }
  double                 GetBound() const          { return m_dBound; }
  // how deep definitions nest here: 1, or one more than the deepest procedural support
  ULONG                  GetDepth() const;
  // the definition's (s, t) at the surface's (u, v)
  SmPoint2d              DefinitionParameters(const SmPoint2d & crUV) const;

  // The definition's partials at (dS, dT): aDerivs[i * (lHighestT + 1) + j] = d^(i+j) S / ds^i dt^j,
  // i <= lHighestS, j <= lHighestT. The work is counted (SmBlendCounts).
  SmStatus EvaluateDefinition(double dS, double dT, ULONG lHighestS, ULONG lHighestT, SmVector3d * aDerivs) const;

  // the counts of every procedural blend's exact evaluations since the last reset
  static SmBlendCounts GetCounts();
  static void          ResetCounts();

  // SmSurface
  virtual SmBoolean operator==(const SmSurface & crOther) const;
  virtual SmStatus  Copy(const SmContext & crContext, SmSurface *& rpNewSurface) const;
  virtual SmStatus  Evaluate(const SmPoint2d & crUV, ULONG lHighestUDeriv, ULONG lHighestVDeriv, SmBoolean bUFromLeft,
                             SmBoolean bVFromLeft, SmBoolean bOnlyUpperHalf, SmVector3d * aDerivatives,
                             SmBoolean bNonZeroTangents = TRUE, SmBoolean bDoZeroSampling = TRUE) const;
  virtual SmStatus  EvaluateSimple(const SmPoint2d & crUV, ULONG lHighestUDeriv, ULONG lHighestVDeriv, SmBoolean bUFromLeft,
                                   SmBoolean bVFromLeft, SmBoolean bOnlyUpperHalf, SmVector3d * aDerivatives) const;
  virtual SmStatus  EvaluatePoint(const SmPoint2d & crUV, SmPoint3d & rPoint) const;
  virtual SmStatus  CalculateBoundingBox(const SmExtent2d & crUVDomain, SmExtent3d * pNormalBox = NULL, SmPseudoBox * pPseudoBox = NULL,
                                         SmPolarBox * pPolarBox = NULL, const SmPseudoBox * pOptPseudoBasisGuess = NULL,
                                         const SmPolarBox * pOptPolarBasisGuess = NULL, SmBoolean bExpandPosBoxesByZoneTol3d = TRUE) const;
  virtual SmStatus  Reparameterize(const SmExtent2d & crNewDomain);
  virtual SmStatus  Reverse(SmSurfParamType eSurfParam);
  virtual void      ToggleSwapUVBit();
  virtual SmStatus  SplitAt(const SmContext & crContext, double dParam, SmSurfParamType eSurfParam, SmSurface *& rpLeftSurface,
                            SmSurface *& rpRightSurface);
  virtual SmStatus  Transform(const SmAxis2Placement & crRotateNMove, const SmVector3d * cpOptScale = NULL);
  virtual SmStatus  WriteToDB(SmDatabaseIO & rDB, ULONG lDBVersionNumber) const;
  static  SmStatus  ReadFromDB(SM_TYPE lType, SmDatabaseIO & rDB, const SmContext & crContext, SmSurface *& rpNewSurface,
                               ULONG lDBVersionNumber);

  SM_COMMON(SmProceduralBlend, SmBSplineSurface, SmProceduralBlend_TYPE);
  virtual void Dump(SmBoolean bAbbrev) const { SmBSplineSurface::Dump(bAbbrev); }
  virtual void Dump(ULONG lLevel) const      { SmBSplineSurface::Dump(lLevel); }
  virtual void Dump(const TCHAR * message) const { SmBSplineSurface::Dump(message); }

} ; // end class SmProceduralBlend

#endif // __SMPROCEDURALBLEND_H__
