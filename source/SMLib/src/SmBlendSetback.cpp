// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

/**********************************************************************//**
* FILE NAME --- SmBlendSetback.cpp
* PURPOSE: SolidCraft change. The setback corner (SmBlendSetback.h): the
*   setback distances, the ease of a round into another section, and the
*   hole's fill by the setback split.
**********************************************************************/

#include "StdAfx.h"

#include <SmBlendSetback.h>

#include <SmBSplineSurface.h>
#include <SmMath.h>

#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

// A refusal: its reason kept where the caller reads it (SER_MSG's text is dropped in release
// builds), and its status returned.
#define SM_SB_FAIL(rpWhy, eStatus, pText) { (rpWhy) = (pText); return (eStatus); }

/*******************************************************************//**
    The ease of a round into another section before its setback.
***********************************************************************/

// Splits the Bezier coefficients crCoef (degree n = size - 1) at t by de Casteljau: rLeft on
// [0, t], rRight on [t, 1].
static void sm_SplitBezier(const SmTArray<double> & crCoef, double t, SmTArray<double> & rLeft, SmTArray<double> & rRight)
{
  ULONG n = crCoef.GetSize() - 1;
  SmTArray<double> sWork(crCoef);
  rLeft = crCoef;
  rRight = crCoef;
  rLeft[0] = sWork[0];
  rRight[n] = sWork[n];
  for (ULONG k = 1; k <= n; k++)
    {
      for (ULONG j = 0; j + k <= n; j++)
        { sWork[j] = (1.0 - t) * sWork[j] + t * sWork[j + 1]; }
      rLeft[k] = sWork[0];
      rRight[n - k] = sWork[n - k];
    }
}

// The Bernstein coefficients on [s0, s1] (0 <= s0 < s1 <= 1) of the polynomial whose
// coefficients on [0, 1] are rCoef: its part left of s1, then that part's right of s0 / s1.
static void sm_RestrictBezier(SmTArray<double> & rCoef, double s0, double s1)
{
  SmTArray<double> sLeft, sRight;
  sm_SplitBezier(rCoef, s1, sLeft, sRight);
  SmTArray<double> sLeft2, sRight2;
  sm_SplitBezier(sLeft, (s1 > 0.0) ? s0 / s1 : 0.0, sLeft2, sRight2);
  rCoef = sRight2;
}

// The law's Bernstein coefficients (degree 5) on the span [a, b] of the rail parameter, given
// the ramps: on a span inside one ramp its smoothstep restricted there; elsewhere the constant
// w is 1 on the side of a ramp's m_dTo, 0 on the side of its m_dFrom, the largest over the
// ramps. A span where two ramps disagree (one easing, another already 1) is refused.
static SmStatus sm_EaseLawOnSpan(const std::vector<SmEaseRamp> & crRamps, double a, double b, SmTArray<double> & rCoef)
{
  rCoef.ReSet();
  for (ULONG i = 0; i < 6; i++)
    { rCoef.Add(0.0); }
  double dConst = 0.0;
  SmBoolean bEasing = FALSE;
  double m = 0.5 * (a + b);
  for (ULONG r = 0; r < crRamps.size(); r++)
    {
      double lo = smos_Min(crRamps[r].m_dFrom, crRamps[r].m_dTo);
      double hi = smos_Max(crRamps[r].m_dFrom, crRamps[r].m_dTo);
      SmBoolean bRising = crRamps[r].m_dFrom < crRamps[r].m_dTo;
      if (m > lo && m < hi)
        {
          if (bEasing)
            { SER(SM_ERR); } // two ramps over one span
          bEasing = TRUE;
          SmTArray<double> sCoef;
          for (ULONG i = 0; i < 6; i++)
            { sCoef.Add((i >= 3) == bRising ? 1.0 : 0.0); }
          sm_RestrictBezier(sCoef, (a - lo) / (hi - lo), (b - lo) / (hi - lo));
          rCoef = sCoef;
        }
      else
        {
          SmBoolean bBeyondTo = bRising ? (m >= hi) : (m <= lo);
          if (bBeyondTo)
            { dConst = 1.0; }
        }
    }
  if (bEasing)
    {
      if (dConst != 0.0)
        { SER(SM_ERR); } // a ramp easing where another is already through
      return SM_SUCCESS;
    }
  for (ULONG i = 0; i < 6; i++)
    { rCoef[i] = dConst; }
  return SM_SUCCESS;
}

static double sm_Binomial(ULONG n, ULONG k)
{
  double d = 1.0;
  for (ULONG i = 1; i <= k; i++)
    { d = d * (double)(n - k + i) / (double)i; }
  return d;
}

// Inserts into rSurf, in direction eDir, the knot dKnot up to multiplicity lMult (counting
// knots within dTol of it as the same).
static SmStatus sm_RaiseKnot(SmBSplineSurface & rSurf, SmSurfParamType eDir, double dKnot, ULONG lMult, double dTol)
{
  ULONG lUDeg, lVDeg;
  SmTArray<SmPoint3d> sPts;
  SmBSplineSurfaceForm eForm;
  SmTArray<ULONG> sUMult, sVMult;
  SmTArray<double> sUKnots, sVKnots, sWeights;
  SmKnotType eKnotType;
  SER(rSurf.GetCanonical(lUDeg, lVDeg, sPts, eForm, sUMult, sVMult, sUKnots, sVKnots, eKnotType, sWeights));
  SmTArray<double> & rKnots = (eDir == SM_SP_U) ? sUKnots : sVKnots;
  SmTArray<ULONG> & rMult = (eDir == SM_SP_U) ? sUMult : sVMult;
  ULONG lHave = 0;
  for (ULONG i = 0; i < rKnots.GetSize(); i++)
    {
      if (smos_Fabs(rKnots[i] - dKnot) <= dTol)
        {
          lHave = rMult[i];
          dKnot = rKnots[i];
        }
    }
  if (lHave < lMult)
    { SER(rSurf.InsertOneKnot(dKnot, lMult - lHave, eDir)); }
  return SM_SUCCESS;
}

// The surface whose sections ease from crFrom's to crTo's over the ramps (SmEaseRamp) along
// eAlong, the rails' direction, in which both run over the same interval with the same rail
// parameter: P = (1 - w) From + w To in homogeneous coordinates, w the C2 law, so P is From
// to second order where w is 0 and To where it is 1. Exact: both surfaces are raised to common
// degrees and knots, put in Bezier form along the rails, and each span's net multiplied by the
// law's Bernstein form (the product of Bernstein polynomials). The span ends keep their knots
// at the new degree's multiplicity; the surface is as smooth across them as its parts are (C2
// where a ramp starts or ends: the law is flat to second order there).
SmStatus SmBlendSetback::EaseSurface(const SmContext & crContext, const SmBSplineSurface & crFrom, const SmBSplineSurface & crTo,
                                     SmSurfParamType eAlong, const std::vector<SmEaseRamp> & crRamps, SmBSplineSurface *& rpEased)
{
  rpEased = NULL;
  SmSurfParamType eAcross = (eAlong == SM_SP_U) ? SM_SP_V : SM_SP_U;
  SmSurface * pCopyA = NULL;
  SmSurface * pCopyB = NULL;
  SER(crFrom.Copy(crContext, pCopyA));
  SmObjDelete sDeleteA(pCopyA);
  SER(crTo.Copy(crContext, pCopyB));
  SmObjDelete sDeleteB(pCopyB);
  SmBSplineSurface * pA = SM_CAST_PTR(SmBSplineSurface, pCopyA); NER(pA);
  SmBSplineSurface * pB = SM_CAST_PTR(SmBSplineSurface, pCopyB); NER(pB);
  SmBSplineSurface * apBoth[2] = { pA, pB };

  // the same domain: along, the rail parameter in both; across, crFrom's
  SmExtent2d sDomA = pA->GetNaturalUVDomain();
  SmExtent2d sDomB = pB->GetNaturalUVDomain();
  SmExtent1d sAlongA = (eAlong == SM_SP_U) ? sDomA.GetUInterval() : sDomA.GetVInterval();
  SmExtent1d sAlongB = (eAlong == SM_SP_U) ? sDomB.GetUInterval() : sDomB.GetVInterval();
  double dParamTol = 1.0e-9 * (1.0 + sAlongA.GetLength());
  if (   smos_Fabs(sAlongA.GetMin() - sAlongB.GetMin()) > dParamTol
      || smos_Fabs(sAlongA.GetMax() - sAlongB.GetMax()) > dParamTol)
    { SER(SM_ERR); }
  SmExtent1d sAcrossA = (eAlong == SM_SP_U) ? sDomA.GetVInterval() : sDomA.GetUInterval();
  SmExtent2d sNewB = (eAlong == SM_SP_U) ? SmExtent2d(SmPoint2d(sAlongA.GetMin(), sAcrossA.GetMin()), SmPoint2d(sAlongA.GetMax(), sAcrossA.GetMax()))
                                         : SmExtent2d(SmPoint2d(sAcrossA.GetMin(), sAlongA.GetMin()), SmPoint2d(sAcrossA.GetMax(), sAlongA.GetMax()));
  SER(pB->Reparameterize(sNewB));
  // the same rail at the across direction's start in both: crTo's turned if it runs the other way
  {
    double t = sAlongA.Evaluate(0.5);
    SmPoint3d sA0, sB0, sB1;
    SmPoint2d sAt0 = (eAlong == SM_SP_U) ? SmPoint2d(t, sAcrossA.GetMin()) : SmPoint2d(sAcrossA.GetMin(), t);
    SmPoint2d sAt1 = (eAlong == SM_SP_U) ? SmPoint2d(t, sAcrossA.GetMax()) : SmPoint2d(sAcrossA.GetMax(), t);
    SER(pA->EvaluatePoint(sAt0, sA0));
    SER(pB->EvaluatePoint(sAt0, sB0));
    SER(pB->EvaluatePoint(sAt1, sB1));
    if (sA0.DistanceBetween(sB1) < sA0.DistanceBetween(sB0))
      {
        SER(pB->Reverse(eAcross));
        SER(pB->Reparameterize(sNewB));
        SmPoint3d sNow;
        SER(pB->EvaluatePoint(sAt0, sNow));
        if (sNow.DistanceBetween(sB1) > 1.0e-9 * (1.0 + sB0.DistanceBetween(sB1)))
          { SER(SM_ERR); }
      }
  }

  // common degrees, the law's 5 added along
  ULONG lAlongDeg = smos_Max(pA->GetDegree(eAlong), pB->GetDegree(eAlong));
  ULONG lAcrossDeg = smos_Max(pA->GetDegree(eAcross), pB->GetDegree(eAcross));
  for (SmBSplineSurface * pS : apBoth)
    {
      if (pS->GetDegree(eAlong) < lAlongDeg)
        { SER(pS->DegreeElevate(eAlong, lAlongDeg)); }
      if (pS->GetDegree(eAcross) < lAcrossDeg)
        { SER(pS->DegreeElevate(eAcross, lAcrossDeg)); }
    }

  // common knots: each surface gets the other's; along, every knot and the ramps' ends raised
  // to the degree (Bezier form)
  for (int iPass = 0; iPass < 2; iPass++)
    {
      SmBSplineSurface * pS = (iPass == 0) ? pA : pB;
      SmBSplineSurface * pO = (iPass == 0) ? pB : pA;
      ULONG lUDeg, lVDeg;
      SmTArray<SmPoint3d> sPts;
      SmBSplineSurfaceForm eForm;
      SmTArray<ULONG> sUMult, sVMult;
      SmTArray<double> sUKnots, sVKnots, sWeights;
      SmKnotType eKnotType;
      SER(pO->GetCanonical(lUDeg, lVDeg, sPts, eForm, sUMult, sVMult, sUKnots, sVKnots, eKnotType, sWeights));
      for (ULONG i = 1; i + 1 < sUKnots.GetSize(); i++)
        { SER(sm_RaiseKnot(*pS, SM_SP_U, sUKnots[i], sUMult[i], dParamTol)); }
      for (ULONG i = 1; i + 1 < sVKnots.GetSize(); i++)
        { SER(sm_RaiseKnot(*pS, SM_SP_V, sVKnots[i], sVMult[i], dParamTol)); }
    }
  for (SmBSplineSurface * pS : apBoth)
    {
      for (ULONG r = 0; r < crRamps.size(); r++)
        {
          const double adEnds[2] = { crRamps[r].m_dFrom, crRamps[r].m_dTo };
          for (double d : adEnds)
            {
              // a ramp ends inside the surface or on its end (the setback, where the trace stopped)
              if (d < sAlongA.GetMin() - dParamTol || d > sAlongA.GetMax() + dParamTol)
                { SER(SM_ERR); }
              if (d > sAlongA.GetMin() + dParamTol && d < sAlongA.GetMax() - dParamTol)
                { SER(sm_RaiseKnot(*pS, eAlong, d, lAlongDeg, dParamTol)); }
            }
        }
    }
  for (SmBSplineSurface * pS : apBoth)
    {
      ULONG lUDeg, lVDeg;
      SmTArray<SmPoint3d> sPts;
      SmBSplineSurfaceForm eForm;
      SmTArray<ULONG> sUMult, sVMult;
      SmTArray<double> sUKnots, sVKnots, sWeights;
      SmKnotType eKnotType;
      SER(pS->GetCanonical(lUDeg, lVDeg, sPts, eForm, sUMult, sVMult, sUKnots, sVKnots, eKnotType, sWeights));
      SmTArray<double> & rKnots = (eAlong == SM_SP_U) ? sUKnots : sVKnots;
      for (ULONG i = 1; i + 1 < rKnots.GetSize(); i++)
        { SER(sm_RaiseKnot(*pS, eAlong, rKnots[i], lAlongDeg, dParamTol)); }
    }

  // the two nets, which now have the same structure
  ULONG lUDegA, lVDegA, lUDegB, lVDegB;
  SmTArray<SmPoint3d> sPtsA, sPtsB;
  SmBSplineSurfaceForm eFormA, eFormB;
  SmTArray<ULONG> sUMultA, sVMultA, sUMultB, sVMultB;
  SmTArray<double> sUKnotsA, sVKnotsA, sWeightsA, sUKnotsB, sVKnotsB, sWeightsB;
  SmKnotType eKnotTypeA, eKnotTypeB;
  SER(pA->GetCanonical(lUDegA, lVDegA, sPtsA, eFormA, sUMultA, sVMultA, sUKnotsA, sVKnotsA, eKnotTypeA, sWeightsA));
  SER(pB->GetCanonical(lUDegB, lVDegB, sPtsB, eFormB, sUMultB, sVMultB, sUKnotsB, sVKnotsB, eKnotTypeB, sWeightsB));
  if (   lUDegA != lUDegB || lVDegA != lVDegB || sPtsA.GetSize() != sPtsB.GetSize()
      || sUKnotsA.GetSize() != sUKnotsB.GetSize() || sVKnotsA.GetSize() != sVKnotsB.GetSize())
    { SER(SM_ERR); }
  for (ULONG i = 0; i < sUKnotsA.GetSize(); i++)
    {
      if (sUMultA[i] != sUMultB[i] || smos_Fabs(sUKnotsA[i] - sUKnotsB[i]) > dParamTol)
        { SER(SM_ERR); }
    }
  for (ULONG i = 0; i < sVKnotsA.GetSize(); i++)
    {
      if (sVMultA[i] != sVMultB[i] || smos_Fabs(sVKnotsA[i] - sVKnotsB[i]) > dParamTol)
        { SER(SM_ERR); }
    }
  SmBoolean bRational = sWeightsA.GetSize() > 0 || sWeightsB.GetSize() > 0;
  SmTArray<double> & rAlongKnots = (eAlong == SM_SP_U) ? sUKnotsA : sVKnotsA;
  SmTArray<ULONG> & rAlongMult = (eAlong == SM_SP_U) ? sUMultA : sVMultA;
  SmTArray<double> & rAcrossKnots = (eAlong == SM_SP_U) ? sVKnotsA : sUKnotsA;
  SmTArray<ULONG> & rAcrossMult = (eAlong == SM_SP_U) ? sVMultA : sUMultA;
  ULONG lSpans = rAlongKnots.GetSize() - 1;
  ULONG lAcrossCount = 0;
  for (ULONG i = 0; i < rAcrossMult.GetSize(); i++)
    { lAcrossCount += rAcrossMult[i]; }
  lAcrossCount -= lAcrossDeg + 1;
  ULONG lAlongCount = lSpans * lAlongDeg + 1;
  if (lAlongCount * lAcrossCount != sPtsA.GetSize())
    { SER(SM_ERR); }
  for (ULONG i = 1; i < lSpans; i++)
    {
      if (rAlongMult[i] != lAlongDeg)
        { SER(SM_ERR); }
    }
  // the homogeneous control point [along i][across j] of a net (U-major lists)
  auto homog = [&](const SmTArray<SmPoint3d> & crPts, const SmTArray<double> & crW, ULONG i, ULONG j, double h[4])
    {
      ULONG k = (eAlong == SM_SP_U) ? i * lAcrossCount + j : j * lAlongCount + i;
      double w = (crW.GetSize() > 0) ? crW[k] : 1.0;
      h[0] = crPts[k].x * w; h[1] = crPts[k].y * w; h[2] = crPts[k].z * w; h[3] = w;
    };

  // each span's product with the law
  ULONG p = lAlongDeg;
  ULONG lNewDeg = p + 5;
  ULONG lNewAlongCount = lSpans * lNewDeg + 1;
  SmTArray<double> sNet(lNewAlongCount * lAcrossCount * 4);
  for (ULONG i = 0; i < lNewAlongCount * lAcrossCount * 4; i++)
    { sNet.Add(0.0); }
  for (ULONG k = 0; k < lSpans; k++)
    {
      SmTArray<double> sW;
      SER(sm_EaseLawOnSpan(crRamps, rAlongKnots[k], rAlongKnots[k + 1], sW));
      for (ULONG j = 0; j < lAcrossCount; j++)
        {
          for (ULONG m = 0; m <= lNewDeg; m++)
            {
              double h[4] = {0.0, 0.0, 0.0, 0.0};
              ULONG iMin = (m > 5) ? m - 5 : 0;
              ULONG iMax = smos_Min(p, m);
              for (ULONG i = iMin; i <= iMax; i++)
                {
                  ULONG l = m - i;
                  double dBeta = sm_Binomial(p, i) * sm_Binomial(5, l) / sm_Binomial(lNewDeg, m);
                  double ha[4], hb[4];
                  homog(sPtsA, sWeightsA, k * p + i, j, ha);
                  homog(sPtsB, sWeightsB, k * p + i, j, hb);
                  for (int c = 0; c < 4; c++)
                    { h[c] += dBeta * ((1.0 - sW[l]) * ha[c] + sW[l] * hb[c]); }
                }
              ULONG lRow = k * lNewDeg + m;
              for (int c = 0; c < 4; c++)
                { sNet[(lRow * lAcrossCount + j) * 4 + c] = h[c]; }
            }
        }
    }

  // the eased surface, back in Euclidean form
  SmTArray<SmPoint3d> sNewPts(lNewAlongCount * lAcrossCount);
  SmTArray<double> sNewWeights;
  for (ULONG i = 0; i < lNewAlongCount * lAcrossCount; i++)
    {
      sNewPts.Add(SmPoint3d(0.0, 0.0, 0.0));
      if (bRational)
        { sNewWeights.Add(1.0); }
    }
  for (ULONG i = 0; i < lNewAlongCount; i++)
    {
      for (ULONG j = 0; j < lAcrossCount; j++)
        {
          const double * h = &sNet[(i * lAcrossCount + j) * 4];
          if (!(h[3] > SM_EFF_ZERO))
            { SER(SM_ERR); }
          ULONG k = (eAlong == SM_SP_U) ? i * lAcrossCount + j : j * lNewAlongCount + i;
          sNewPts[k] = SmPoint3d(h[0] / h[3], h[1] / h[3], h[2] / h[3]);
          if (bRational)
            { sNewWeights[k] = h[3]; }
        }
    }
  SmTArray<ULONG> sNewAlongMult;
  for (ULONG i = 0; i <= lSpans; i++)
    { sNewAlongMult.Add((i == 0 || i == lSpans) ? lNewDeg + 1 : lNewDeg); }
  SmBSplineSurface * pEased = NULL;
  if (eAlong == SM_SP_U)
    {
      SER(SmBSplineSurface::CreateCanonical(crContext, lNewDeg, lAcrossDeg, sNewPts, SM_SF_UNSPECIFIED, sNewAlongMult, rAcrossMult,
                                            rAlongKnots, rAcrossKnots, SM_KT_UNSPECIFIED, bRational ? &sNewWeights : NULL, NULL, pEased));
    }
  else
    {
      SER(SmBSplineSurface::CreateCanonical(crContext, lAcrossDeg, lNewDeg, sNewPts, SM_SF_UNSPECIFIED, rAcrossMult, sNewAlongMult,
                                            rAcrossKnots, rAlongKnots, SM_KT_UNSPECIFIED, bRational ? &sNewWeights : NULL, NULL, pEased));
    }
  NER(pEased);
  rpEased = pEased;
  return SM_SUCCESS;
}

/*******************************************************************//**
    The setback distances.
***********************************************************************/

SmStatus SmBlendSetback::Distances(const std::vector<Face> & crFaces, const std::vector<Edge> & crEdges, std::vector<double> & rvSetbacks,
                                   std::vector<double> * pvLeast, const TCHAR ** ppReason)
{
  const TCHAR * pDummy = NULL;
  const TCHAR *& pWhy = (ppReason != NULL) ? *ppReason : pDummy;
  pWhy = NULL;
  rvSetbacks.clear();
  if (pvLeast != NULL)
    { pvLeast->clear(); }
  const size_t lF = crFaces.size(), lE = crEdges.size();
  // each edge's rail distance from it, the same on both its faces (a round between planes)
  std::vector<double> sRail(lE, 0.0);
  for (size_t i=0; i<lE; i++)
    {
      const Edge & crE = crEdges[i];
      if (!(crE.m_dRadius > 0.0) || !std::isfinite(crE.m_dRadius) || crE.m_alFaces[0] >= lF || crE.m_alFaces[1] >= lF
          || crE.m_alFaces[0] == crE.m_alFaces[1])
        { SM_SB_FAIL(pWhy, SM_ERR_INVALID_INPUT, _T("SmBlendSetback: an edge needs a positive radius and two of the faces given")); }
      const double dCos = smos_Max(-1.0, smos_Min(1.0, crFaces[crE.m_alFaces[0]].m_sNormal.Dot(crFaces[crE.m_alFaces[1]].m_sNormal)));
      sRail[i] = crE.m_dRadius * std::tan(0.5 * std::acos(dCos));
    }
  std::vector<double> sLeast(lE, 0.0);
  for (size_t f=0; f<lF; f++)
    {
      std::vector<size_t> sBeside;
      for (size_t i=0; i<lE; i++)
        if (crEdges[i].m_alFaces[0] == f || crEdges[i].m_alFaces[1] == f)
          { sBeside.push_back(i); }
      if (sBeside.size() != 2)
        { SM_SB_FAIL(pWhy, SM_ERR_INVALID_INPUT, _T("SmBlendSetback: a face at the corner must lie between exactly two of its rounded edges")); }
      const size_t i = sBeside[0], j = sBeside[1];
      const double dCos = smos_Max(-1.0, smos_Min(1.0, crEdges[i].m_sAway.Dot(crEdges[j].m_sAway)));
      const double dSin = smos_Sqrt(smos_Max(0.0, 1.0 - dCos * dCos));
      if (!(dSin > 1.0e-9))
        { SM_SB_FAIL(pWhy, SM_ERR_INVALID_INPUT, _T("SmBlendSetback: two edges on a face run along one line")); }
      sLeast[i] = smos_Max(sLeast[i], (sRail[j] + sRail[i] * dCos) / dSin);
      sLeast[j] = smos_Max(sLeast[j], (sRail[i] + sRail[j] * dCos) / dSin);
    }
  for (size_t i=0; i<lE; i++)
    {
      const double r = crEdges[i].m_dRadius;
      rvSetbacks.push_back(smos_Max(1.5 * r, sLeast[i] + 0.25 * r));
    }
  if (pvLeast != NULL)
    { *pvLeast = sLeast; }
  return SM_SUCCESS;
}

/*******************************************************************//**
    The hole's fill.
***********************************************************************/

SmBlendSetback::Fill::~Fill()
{
  for (size_t i=0; i<m_vPatches.size(); i++)
    { delete m_vPatches[i]; }
}

namespace
{

// A point in homogeneous coordinates (x w, y w, z w, w).
struct sm_SbH
{
  double a[4];
};

inline sm_SbH sm_SbZero() { sm_SbH h = {{0.0, 0.0, 0.0, 0.0}}; return h; }
inline sm_SbH sm_SbAxpy(const sm_SbH & crY, double s, const sm_SbH & crX)
{
  sm_SbH h;
  for (int c=0; c<4; c++) { h.a[c] = crY.a[c] + s * crX.a[c]; }
  return h;
}

// Bernstein coefficients of a polynomial on [0, 1] (degree size - 1).
typedef std::vector<double> sm_SbB;

sm_SbB sm_SbMul(const sm_SbB & crA, const sm_SbB & crB)
{
  const size_t p = crA.size() - 1, q = crB.size() - 1;
  sm_SbB sC(p + q + 1, 0.0);
  for (size_t i=0; i<=p; i++)
    for (size_t j=0; j<=q; j++)
      { sC[i + j] += sm_Binomial(static_cast<ULONG>(p), static_cast<ULONG>(i)) * sm_Binomial(static_cast<ULONG>(q), static_cast<ULONG>(j))
                     / sm_Binomial(static_cast<ULONG>(p + q), static_cast<ULONG>(i + j)) * crA[i] * crB[j]; }
  return sC;
}

// raised to degree lDeg (>= its own)
sm_SbB sm_SbElevate(const sm_SbB & crA, size_t lDeg)
{
  sm_SbB sA(crA);
  while (sA.size() - 1 < lDeg)
    {
      const size_t p = sA.size() - 1;
      sm_SbB sB(p + 2, 0.0);
      sB[0] = sA[0];
      sB[p + 1] = sA[p];
      for (size_t i=1; i<=p; i++)
        { sB[i] = (static_cast<double>(i) / static_cast<double>(p + 1)) * sA[i - 1] + (1.0 - static_cast<double>(i) / static_cast<double>(p + 1)) * sA[i]; }
      sA.swap(sB);
    }
  return sA;
}

// Tensor-product Bernstein coefficients, degree (p, q), c[i (q + 1) + j].
struct sm_SbT
{
  size_t p, q;
  std::vector<double> c;
};

sm_SbT sm_SbTMul(const sm_SbT & crA, const sm_SbT & crB)
{
  sm_SbT sC;
  sC.p = crA.p + crB.p;
  sC.q = crA.q + crB.q;
  sC.c.assign((sC.p + 1) * (sC.q + 1), 0.0);
  for (size_t i=0; i<=crA.p; i++)
    for (size_t j=0; j<=crA.q; j++)
      for (size_t k=0; k<=crB.p; k++)
        for (size_t l=0; l<=crB.q; l++)
          {
            const double f = sm_Binomial(static_cast<ULONG>(crA.p), static_cast<ULONG>(i)) * sm_Binomial(static_cast<ULONG>(crB.p), static_cast<ULONG>(k))
                             / sm_Binomial(static_cast<ULONG>(sC.p), static_cast<ULONG>(i + k))
                             * sm_Binomial(static_cast<ULONG>(crA.q), static_cast<ULONG>(j)) * sm_Binomial(static_cast<ULONG>(crB.q), static_cast<ULONG>(l))
                             / sm_Binomial(static_cast<ULONG>(sC.q), static_cast<ULONG>(j + l));
            sC.c[(i + k) * (sC.q + 1) + (j + l)] += f * crA.c[i * (crA.q + 1) + j] * crB.c[k * (crB.q + 1) + l];
          }
  return sC;
}

sm_SbT sm_SbTElevate(const sm_SbT & crA, size_t lP, size_t lQ)
{
  sm_SbT sA;
  sA.p = lP;
  sA.q = crA.q;
  sA.c.assign((lP + 1) * (crA.q + 1), 0.0);
  for (size_t j=0; j<=crA.q; j++)
    {
      sm_SbB sCol(crA.p + 1);
      for (size_t i=0; i<=crA.p; i++) { sCol[i] = crA.c[i * (crA.q + 1) + j]; }
      sCol = sm_SbElevate(sCol, lP);
      for (size_t i=0; i<=lP; i++) { sA.c[i * (crA.q + 1) + j] = sCol[i]; }
    }
  sm_SbT sB;
  sB.p = lP;
  sB.q = lQ;
  sB.c.assign((lP + 1) * (lQ + 1), 0.0);
  for (size_t i=0; i<=lP; i++)
    {
      sm_SbB sRow(sA.c.begin() + static_cast<long>(i * (crA.q + 1)), sA.c.begin() + static_cast<long>((i + 1) * (crA.q + 1)));
      sRow = sm_SbElevate(sRow, lQ);
      for (size_t j=0; j<=lQ; j++) { sB.c[i * (lQ + 1) + j] = sRow[j]; }
    }
  return sB;
}

// The Bernstein basis of degree p at t, and its first and second derivatives.
void sm_SbBasis(size_t p, double t, std::vector<double> & rB0, std::vector<double> & rB1, std::vector<double> & rB2)
{
  // the basis of every degree up to p, by the recurrence
  std::vector< std::vector<double> > sAll(p + 1);
  sAll[0].assign(1, 1.0);
  for (size_t d=1; d<=p; d++)
    {
      sAll[d].assign(d + 1, 0.0);
      for (size_t i=0; i<=d; i++)
        {
          const double dLeft = (i > 0) ? sAll[d - 1][i - 1] : 0.0;
          const double dHere = (i < d) ? sAll[d - 1][i] : 0.0;
          sAll[d][i] = t * dLeft + (1.0 - t) * dHere;
        }
    }
  rB0 = sAll[p];
  rB1.assign(p + 1, 0.0);
  rB2.assign(p + 1, 0.0);
  const double dp = static_cast<double>(p);
  for (size_t i=0; i<=p && p >= 1; i++)
    {
      const double dLeft = (i > 0) ? sAll[p - 1][i - 1] : 0.0;
      const double dHere = (i < p) ? sAll[p - 1][i] : 0.0;
      rB1[i] = dp * (dLeft - dHere);
    }
  for (size_t i=0; i<=p && p >= 2; i++)
    {
      const double d2 = (i >= 2) ? sAll[p - 2][i - 2] : 0.0;
      const double d1 = (i >= 1 && i - 1 <= p - 2) ? sAll[p - 2][i - 1] : 0.0;
      const double d0 = (i <= p - 2) ? sAll[p - 2][i] : 0.0;
      rB2[i] = dp * (dp - 1.0) * (d2 - 2.0 * d1 + d0);
    }
}

// A homogeneous tensor-product Bezier net, (p + 1) x (q + 1) points, P[i (q + 1) + j], i along u.
struct sm_SbNet
{
  size_t p, q;
  std::vector<sm_SbH> P;
};

// The point and its second derivatives R_uu, R_uv, R_vv (Euclidean) at (u, v).
void sm_SbSecond(const sm_SbNet & crN, double u, double v, double (&rR)[3][3])
{
  std::vector<double> sU0, sU1, sU2, sV0, sV1, sV2;
  sm_SbBasis(crN.p, u, sU0, sU1, sU2);
  sm_SbBasis(crN.q, v, sV0, sV1, sV2);
  // homogeneous derivatives: 0 value, 1 u, 2 v, 3 uu, 4 uv, 5 vv
  double h[6][4] = {{0}};
  for (size_t i=0; i<=crN.p; i++)
    for (size_t j=0; j<=crN.q; j++)
      {
        const double * a = crN.P[i * (crN.q + 1) + j].a;
        const double b[6] = { sU0[i] * sV0[j], sU1[i] * sV0[j], sU0[i] * sV1[j], sU2[i] * sV0[j], sU1[i] * sV1[j], sU0[i] * sV2[j] };
        for (int k=0; k<6; k++)
          for (int c=0; c<4; c++)
            { h[k][c] += b[k] * a[c]; }
      }
  const double w = h[0][3];
  for (int c=0; c<3; c++)
    {
      const double R = h[0][c] / w;
      const double Ru = (h[1][c] - R * h[1][3]) / w;
      const double Rv = (h[2][c] - R * h[2][3]) / w;
      rR[0][c] = (h[3][c] - 2.0 * Ru * h[1][3] - R * h[3][3]) / w;
      rR[1][c] = (h[4][c] - Ru * h[2][3] - Rv * h[1][3] - R * h[4][3]) / w;
      rR[2][c] = (h[5][c] - 2.0 * Rv * h[2][3] - R * h[5][3]) / w;
    }
}

// The point, its first derivatives (Euclidean) at (u, v).
void sm_SbFirst(const sm_SbNet & crN, double u, double v, SmVector3d & rP, SmVector3d & rU, SmVector3d & rV)
{
  std::vector<double> sU0, sU1, sU2, sV0, sV1, sV2;
  sm_SbBasis(crN.p, u, sU0, sU1, sU2);
  sm_SbBasis(crN.q, v, sV0, sV1, sV2);
  double h[3][4] = {{0}};
  for (size_t i=0; i<=crN.p; i++)
    for (size_t j=0; j<=crN.q; j++)
      {
        const double * a = crN.P[i * (crN.q + 1) + j].a;
        for (int c=0; c<4; c++)
          {
            h[0][c] += sU0[i] * sV0[j] * a[c];
            h[1][c] += sU1[i] * sV0[j] * a[c];
            h[2][c] += sU0[i] * sV1[j] * a[c];
          }
      }
  const double w = h[0][3];
  rP = SmVector3d(h[0][0] / w, h[0][1] / w, h[0][2] / w);
  rU = SmVector3d((h[1][0] - rP.x * h[1][3]) / w, (h[1][1] - rP.y * h[1][3]) / w, (h[1][2] - rP.z * h[1][3]) / w);
  rV = SmVector3d((h[2][0] - rP.x * h[2][3]) / w, (h[2][1] - rP.y * h[2][3]) / w, (h[2][2] - rP.z * h[2][3]) / w);
}

// A round's homogeneous derivatives along its rails at its end section, H[j] (j = 0 .. k) the
// across Bezier curve of the j-th derivative by the rail parameter, columns from its across
// start; and its across degree.
SmStatus sm_SbRoundJets(const SmContext & crContext, const SmBlendSetback::Round & crRound, ULONG k,
                        std::vector< std::vector<sm_SbH> > & rH, size_t & rlQ, const TCHAR *& rpReason)
{
  rH.clear();
  NER(crRound.m_pSurface);
  const SmSurfParamType eAlong = crRound.m_eAlong;
  if (eAlong != SM_SP_U && eAlong != SM_SP_V)
    { SM_SB_FAIL(rpReason, SM_ERR_INVALID_INPUT, _T("SmBlendSetback: a round's rails run along u or v")); }
  const SmSurfParamType eAcross = (eAlong == SM_SP_U) ? SM_SP_V : SM_SP_U;
  SmSurface * pCopy = NULL;
  SER(crRound.m_pSurface->Copy(crContext, pCopy));
  SmObjDelete sDelete(pCopy);
  SmBSplineSurface * pS = SM_CAST_PTR(SmBSplineSurface, pCopy); NER(pS);
  const ULONG p = pS->GetDegree(eAlong);
  const SmExtent2d sDom = pS->GetNaturalUVDomain();
  const SmExtent1d sAlong = (eAlong == SM_SP_U) ? sDom.GetUInterval() : sDom.GetVInterval();
  const double dEnd = crRound.m_dEnd;
  const double dParamTol = 1.0e-12 * (1.0 + sAlong.GetLength() + smos_Fabs(dEnd));
  // the round lies below its end when the corner lies above it
  if (crRound.m_bCornerAbove ? !(dEnd > sAlong.GetMin() + dParamTol && dEnd <= sAlong.GetMax() + dParamTol)
                             : !(dEnd >= sAlong.GetMin() - dParamTol && dEnd < sAlong.GetMax() - dParamTol))
    { SM_SB_FAIL(rpReason, SM_ERR_INVALID_INPUT, _T("SmBlendSetback: a round's end section lies outside its rails' interval")); }
  SER(sm_RaiseKnot(*pS, eAlong, dEnd, p, dParamTol));

  ULONG lUDeg, lVDeg;
  SmTArray<SmPoint3d> sPts;
  SmBSplineSurfaceForm eForm;
  SmTArray<ULONG> sUMult, sVMult;
  SmTArray<double> sUKnots, sVKnots, sWeights;
  SmKnotType eKnotType;
  SER(pS->GetCanonical(lUDeg, lVDeg, sPts, eForm, sUMult, sVMult, sUKnots, sVKnots, eKnotType, sWeights));
  {
    const SmTArray<double> & crK = (eAlong == SM_SP_U) ? sUKnots : sVKnots;
    ULONG ia = 0;
    while (ia < crK.GetSize() && smos_Fabs(crK[ia] - dEnd) > dParamTol) { ia++; }
    if (ia == crK.GetSize())
      { SER(SM_ERR); }
    // the span on the round's side, raised to Bezier form at its other end
    const ULONG ib = crRound.m_bCornerAbove ? ia - 1 : ia + 1;
    if (ib >= crK.GetSize())
      { SER(SM_ERR); }
    SER(sm_RaiseKnot(*pS, eAlong, crK[ib], p, dParamTol));
  }
  SER(pS->GetCanonical(lUDeg, lVDeg, sPts, eForm, sUMult, sVMult, sUKnots, sVKnots, eKnotType, sWeights));
  const SmTArray<double> & crAcrossK = (eAlong == SM_SP_U) ? sVKnots : sUKnots;
  if (crAcrossK.GetSize() != 2)
    { SM_SB_FAIL(rpReason, SM_ERR_INVALID_INPUT, _T("SmBlendSetback: a round must be one Bezier span across its rails")); }
  const SmTArray<double> & crK = (eAlong == SM_SP_U) ? sUKnots : sVKnots;
  const SmTArray<ULONG> & crM = (eAlong == SM_SP_U) ? sUMult : sVMult;
  ULONG ia = 0;
  while (ia < crK.GetSize() && smos_Fabs(crK[ia] - dEnd) > dParamTol) { ia++; }
  if (ia == crK.GetSize())
    { SER(SM_ERR); }
  auto rowOf = [&](ULONG lKnot) -> ULONG
    {
      if (lKnot == 0) { return 0; }
      ULONG s = 0;
      for (ULONG i=0; i<lKnot; i++) { s += crM[i]; }
      return s - 1;
    };
  const ULONG lQ = pS->GetDegree(eAcross);
  const ULONG lAlongCount = static_cast<ULONG>(sPts.GetSize()) / (lQ + 1);
  auto at = [&](ULONG i, ULONG j) -> sm_SbH
    {
      const ULONG idx = (eAlong == SM_SP_U) ? i * (lQ + 1) + j : j * lAlongCount + i;
      const double w = (sWeights.GetSize() > 0) ? sWeights[idx] : 1.0;
      sm_SbH h = {{ sPts[idx].x * w, sPts[idx].y * w, sPts[idx].z * w, w }};
      return h;
    };
  const ULONG lRow = rowOf(ia);
  const double dDelta = crRound.m_bCornerAbove ? crK[ia] - crK[ia - 1] : crK[ia + 1] - crK[ia];
  if (!(dDelta > 0.0))
    { SER(SM_ERR); }
  // the span's rows in increasing rail parameter, Q[0 .. p]
  std::vector< std::vector<sm_SbH> > sQ(p + 1, std::vector<sm_SbH>(lQ + 1));
  const ULONG lFirst = crRound.m_bCornerAbove ? lRow - p : lRow;
  if (crRound.m_bCornerAbove && lRow < p)
    { SER(SM_ERR); }
  for (ULONG i=0; i<=p; i++)
    for (ULONG j=0; j<=lQ; j++)
      {
        if (lFirst + i >= lAlongCount)
          { SER(SM_ERR); }
        sQ[i][j] = at(lFirst + i, j);
      }
  rH.assign(k + 1, std::vector<sm_SbH>(lQ + 1, sm_SbZero()));
  for (ULONG j=0; j<=lQ; j++)
    {
      // derivatives at the span's end the round meets its end section at: forward differences
      // from that end, signed by the parameter
      const bool bAtTop = crRound.m_bCornerAbove ? true : false;
      for (ULONG d=0; d<=k; d++)
        {
          if (d > p)
            { continue; }
          sm_SbH h = sm_SbZero();
          double dFall = 1.0;
          for (ULONG l=0; l<d; l++) { dFall *= static_cast<double>(p - l) / dDelta; }
          for (ULONG l=0; l<=d; l++)
            {
              const double dSign = ((d - l) % 2 == 0) ? 1.0 : -1.0;
              const double c = dFall * sm_Binomial(d, l) * dSign;
              // at the top: sum_l (-1)^(d-l) C(d,l) Q[p - d + l]; at the bottom: sum_l (-1)^(d-l) C(d,l) Q[l]
              const sm_SbH & crQ = bAtTop ? sQ[p - d + l][j] : sQ[l][j];
              h = sm_SbAxpy(h, c, crQ);
            }
          rH[d][j] = h;
        }
    }
  rlQ = lQ;
  return SM_SUCCESS;
}

// The monomials x^a y^b, a + b <= d, in a fixed order.
std::vector< std::pair<int, int> > sm_SbMonomials(int d)
{
  std::vector< std::pair<int, int> > sM;
  for (int t=0; t<=d; t++)
    for (int a=t; a>=0; a--)
      { sM.push_back(std::make_pair(a, t - a)); }
  return sM;
}

// x^a y^b along the line (x(v), y(v)) (Bernstein, degree 1 each), raised to degree lDeg.
sm_SbB sm_SbAlong(int a, int b, const sm_SbB & crX, const sm_SbB & crY, size_t lDeg)
{
  sm_SbB s(1, 1.0);
  for (int i=0; i<a; i++) { s = sm_SbMul(s, crX); }
  for (int i=0; i<b; i++) { s = sm_SbMul(s, crY); }
  return sm_SbElevate(s, lDeg);
}

// x^a y^b over the bilinear patch (x(s, t), y(s, t)), raised to degree (lDeg, lDeg).
sm_SbT sm_SbOver(int a, int b, const sm_SbT & crX, const sm_SbT & crY, size_t lDeg)
{
  sm_SbT s;
  s.p = 0; s.q = 0; s.c.assign(1, 1.0);
  for (int i=0; i<a; i++) { s = sm_SbTMul(s, crX); }
  for (int i=0; i<b; i++) { s = sm_SbTMul(s, crY); }
  return sm_SbTElevate(s, lDeg, lDeg);
}

// Gauss-Legendre, 8 points on [0, 1].
const double SM_SB_GX[8] = { 0.019855071751231856, 0.10166676129318664, 0.2372337950418355, 0.4082826787521751,
                             0.5917173212478249, 0.7627662049581645, 0.8983332387068134, 0.9801449282487681 };
const double SM_SB_GW[8] = { 0.05061426814518813, 0.11119051722668724, 0.15685332293894363, 0.18134189168918100,
                             0.18134189168918100, 0.15685332293894363, 0.11119051722668724, 0.05061426814518813 };

// The rows of A (m x n, row-major) reduced by Gaussian elimination with complete pivoting: the
// solutions of A x = b are x = rX + rZ y, rZ n x (n - rank) column-major by column. FALSE when
// the equations are inconsistent beyond dTol (rows scaled to unit largest entry).
SmBoolean sm_SbAffineSolutions(std::vector<double> A, std::vector<double> b, size_t m, size_t n, double dTol,
                               std::vector<double> & rX, std::vector< std::vector<double> > & rZ)
{
  for (size_t i=0; i<m; i++)
    {
      double dMax = 0.0;
      for (size_t j=0; j<n; j++) { dMax = smos_Max(dMax, smos_Fabs(A[i * n + j])); }
      if (dMax > 0.0)
        {
          for (size_t j=0; j<n; j++) { A[i * n + j] /= dMax; }
          b[i] /= dMax;
        }
    }
  std::vector<size_t> sCol(n);
  for (size_t j=0; j<n; j++) { sCol[j] = j; }
  size_t r = 0;
  for (; r<m && r<n; r++)
    {
      size_t pi = r, pj = r;
      double dBest = 0.0;
      for (size_t i=r; i<m; i++)
        for (size_t j=r; j<n; j++)
          if (smos_Fabs(A[i * n + j]) > dBest) { dBest = smos_Fabs(A[i * n + j]); pi = i; pj = j; }
      if (!(dBest > 1.0e-10))
        { break; }
      if (pi != r)
        {
          for (size_t j=0; j<n; j++) { std::swap(A[pi * n + j], A[r * n + j]); }
          std::swap(b[pi], b[r]);
        }
      if (pj != r)
        {
          for (size_t i=0; i<m; i++) { std::swap(A[i * n + pj], A[i * n + r]); }
          std::swap(sCol[pj], sCol[r]);
        }
      for (size_t i=0; i<m; i++)
        {
          if (i == r) { continue; }
          const double f = A[i * n + r] / A[r * n + r];
          if (f == 0.0) { continue; }
          for (size_t j=r; j<n; j++) { A[i * n + j] -= f * A[r * n + j]; }
          b[i] -= f * b[r];
        }
    }
  for (size_t i=r; i<m; i++)
    if (smos_Fabs(b[i]) > dTol)
      { return FALSE; }
  // pivot rows: A[i][i] x_i + sum_{free j} A[i][j] x_j = b_i
  rX.assign(n, 0.0);
  for (size_t i=0; i<r; i++) { rX[sCol[i]] = b[i] / A[i * n + i]; }
  rZ.clear();
  for (size_t f=r; f<n; f++)
    {
      std::vector<double> z(n, 0.0);
      z[sCol[f]] = 1.0;
      for (size_t i=0; i<r; i++) { z[sCol[i]] = -A[i * n + f] / A[i * n + i]; }
      rZ.push_back(z);
    }
  return TRUE;
}

// The least squares solution of M y = c (M rows x cols, row-major, rows >= cols) by
// Householder QR; FALSE when M's columns are dependent (relative to dTol).
SmBoolean sm_SbLeastSquares(std::vector<double> M, std::vector<double> c, size_t lRows, size_t lCols, double dTol, std::vector<double> & rY)
{
  std::vector<double> sDiag(lCols, 0.0);
  for (size_t k=0; k<lCols; k++)
    {
      double dNorm = 0.0;
      for (size_t i=k; i<lRows; i++) { dNorm += M[i * lCols + k] * M[i * lCols + k]; }
      dNorm = smos_Sqrt(dNorm);
      if (!(dNorm > 0.0))
        { return FALSE; }
      const double dAlpha = (M[k * lCols + k] > 0.0) ? -dNorm : dNorm;
      std::vector<double> v(lRows - k);
      for (size_t i=k; i<lRows; i++) { v[i - k] = M[i * lCols + k]; }
      v[0] -= dAlpha;
      double dV = 0.0;
      for (size_t i=0; i<v.size(); i++) { dV += v[i] * v[i]; }
      if (dV > 0.0)
        {
          for (size_t j=k; j<lCols; j++)
            {
              double s = 0.0;
              for (size_t i=k; i<lRows; i++) { s += v[i - k] * M[i * lCols + j]; }
              s = 2.0 * s / dV;
              for (size_t i=k; i<lRows; i++) { M[i * lCols + j] -= s * v[i - k]; }
            }
          double s = 0.0;
          for (size_t i=k; i<lRows; i++) { s += v[i - k] * c[i]; }
          s = 2.0 * s / dV;
          for (size_t i=k; i<lRows; i++) { c[i] -= s * v[i - k]; }
        }
      sDiag[k] = M[k * lCols + k];
    }
  double dMax = 0.0;
  for (size_t k=0; k<lCols; k++) { dMax = smos_Max(dMax, smos_Fabs(sDiag[k])); }
  for (size_t k=0; k<lCols; k++)
    if (!(smos_Fabs(sDiag[k]) > dTol * dMax))
      { return FALSE; }
  rY.assign(lCols, 0.0);
  for (size_t k=lCols; k-- > 0; )
    {
      double s = c[k];
      for (size_t j=k+1; j<lCols; j++) { s -= M[k * lCols + j] * rY[j]; }
      rY[k] = s / M[k * lCols + k];
    }
  return TRUE;
}

} // namespace

/*******************************************************************//**
PURPOSE: The hole's fill by the setback split (SmBlendSetback.h).

NOTES: The cycle: round m between faces m (its setback patch's v = 0) and
    m + 1. The domain: the regular n-gon with vertices B_m = (cos 2 pi m / n,
    sin 2 pi m / n), B_m the base point of face m; its edge m from B_m to
    B_m+1 is the inner curve of setback patch m; at B_m the direction field
    of edge m starts as sigma_m = t_m + t_m-1 (t the edges' unit directions)
    and that of edge m - 1 ends as -sigma_m, so the two setback patches leave
    face m's base point along one line (the spring curve is tangent there).
***********************************************************************/
SmStatus SmBlendSetback::FillHole(const SmContext & crContext, const std::vector<Face> & crFaces, const std::vector<Round> & crRounds, ULONG lOrder,
                                  double dTol, Fill & rFill)
{
  for (size_t i=0; i<rFill.m_vPatches.size(); i++) { delete rFill.m_vPatches[i]; }
  rFill.m_vPatches.clear();
  rFill.m_vRound.clear();
  rFill.m_vFace.clear();
  rFill.m_vReversed.clear();
  rFill.m_lDegree = 0;
  rFill.m_pReason = NULL;
  if (lOrder != 1 && lOrder != 2)
    { SM_SB_FAIL(rFill.m_pReason, SM_ERR_INVALID_INPUT, _T("SmBlendSetback: a setback fill is of order 1 (G1) or 2 (G2)")); }
  if (!(dTol > 0.0) || crRounds.size() < 3 || crFaces.size() != crRounds.size())
    { SM_SB_FAIL(rFill.m_pReason, SM_ERR_INVALID_INPUT, _T("SmBlendSetback: a setback fill needs three or more rounds, as many faces, and a positive tolerance")); }
  const size_t n = crRounds.size();
  const ULONG k = lOrder;
  const size_t N = 2 * k + 1; // the setback patches' degree across the hole

  // each round's jets at its end, and the faces its across ends lie on
  std::vector< std::vector< std::vector<sm_SbH> > > sJets(n);
  std::vector<size_t> sQ(n);
  std::vector<size_t> sFaceAt[2];
  sFaceAt[0].assign(n, n);
  sFaceAt[1].assign(n, n);
  for (size_t r=0; r<n; r++)
    {
      SER(sm_SbRoundJets(crContext, crRounds[r], k, sJets[r], sQ[r], rFill.m_pReason));
      for (int e=0; e<2; e++)
        {
          const sm_SbH & crH = sJets[r][0][e == 0 ? 0 : sQ[r]];
          const SmPoint3d sP(crH.a[0] / crH.a[3], crH.a[1] / crH.a[3], crH.a[2] / crH.a[3]);
          for (size_t f=0; f<n; f++)
            if (smos_Fabs((sP - crFaces[f].m_sPoint).Dot(crFaces[f].m_sNormal)) <= dTol)
              {
                if (sFaceAt[e][r] != n)
                  { SM_SB_FAIL(rFill.m_pReason, SM_ERR_INVALID_INPUT, _T("SmBlendSetback: a round's rail lies on two of the faces")); }
                sFaceAt[e][r] = f;
              }
          if (sFaceAt[e][r] == n)
            { SM_SB_FAIL(rFill.m_pReason, SM_ERR_INVALID_INPUT, _T("SmBlendSetback: a round's rail lies on none of the faces")); }
        }
    }
  // the cycle, face to face
  std::vector<ULONG> sRound, sFace;
  std::vector<SmBoolean> sRev;
  {
    std::vector<bool> sUsed(n, false);
    size_t r = 0;
    bool bRev = false;
    sFace.push_back(static_cast<ULONG>(sFaceAt[0][0]));
    for (size_t m=0; m<n; m++)
      {
        sUsed[r] = true;
        sRound.push_back(static_cast<ULONG>(r));
        sRev.push_back(bRev ? TRUE : FALSE);
        const size_t lNext = bRev ? sFaceAt[0][r] : sFaceAt[1][r];
        if (m + 1 == n)
          {
            if (lNext != sFace[0])
              { SM_SB_FAIL(rFill.m_pReason, SM_ERR_INVALID_INPUT, _T("SmBlendSetback: the rounds do not close round the vertex face to face")); }
            break;
          }
        sFace.push_back(static_cast<ULONG>(lNext));
        size_t lFound = n;
        for (size_t o=0; o<n; o++)
          if (!sUsed[o] && (sFaceAt[0][o] == lNext || sFaceAt[1][o] == lNext))
            {
              if (lFound != n)
                { SM_SB_FAIL(rFill.m_pReason, SM_ERR_INVALID_INPUT, _T("SmBlendSetback: a face lies between more than two of the rounds")); }
              lFound = o;
            }
        if (lFound == n)
          { SM_SB_FAIL(rFill.m_pReason, SM_ERR_INVALID_INPUT, _T("SmBlendSetback: the rounds do not close round the vertex face to face")); }
        bRev = (sFaceAt[1][lFound] == lNext);
        r = lFound;
      }
    std::vector<bool> sFaceUsed(n, false);
    for (size_t m=0; m<n; m++)
      {
        if (sFaceUsed[sFace[m]])
          { SM_SB_FAIL(rFill.m_pReason, SM_ERR_INVALID_INPUT, _T("SmBlendSetback: the rounds do not close round the vertex face to face")); }
        sFaceUsed[sFace[m]] = true;
      }
  }

  // the domain
  std::vector<SmPoint2d> sB(n);
  for (size_t m=0; m<n; m++)
    { sB[m] = SmPoint2d(std::cos(2.0 * SM_PI * static_cast<double>(m) / static_cast<double>(n)), std::sin(2.0 * SM_PI * static_cast<double>(m) / static_cast<double>(n))); }
  std::vector<SmVector2d> sT(n), sSigma(n);
  for (size_t m=0; m<n; m++)
    {
      sT[m] = sB[(m + 1) % n] - sB[m];
      SER(sT[m].Unitize());
    }
  for (size_t m=0; m<n; m++)
    { sSigma[m] = sT[m] + sT[(m + n - 1) % n]; }

  // the setback patches' first k + 1 rows, from the rounds: homogeneous Hermite data at u = 0,
  // the rail parameter run toward the corner at a speed lambda giving the patch a speed across
  // of half its spoke (from its end section's middle to the middle of the end sections' middles)
  std::vector< std::vector< std::vector<sm_SbH> > > sRoundRows(n); // [m][i][j], i = 0 .. k, degree sQ across
  {
    std::vector<SmPoint3d> sMid(n);
    SmVector3d sCentre(0.0, 0.0, 0.0);
    std::vector<double> sSpeed(n);
    for (size_t m=0; m<n; m++)
      {
        const Round & crR = crRounds[sRound[m]];
        SmVector3d sD[4];
        const SmExtent2d sDom = crR.m_pSurface->GetNaturalUVDomain();
        const SmExtent1d sAcross = (crR.m_eAlong == SM_SP_U) ? sDom.GetVInterval() : sDom.GetUInterval();
        const SmPoint2d sUV = (crR.m_eAlong == SM_SP_U) ? SmPoint2d(crR.m_dEnd, sAcross.Evaluate(0.5)) : SmPoint2d(sAcross.Evaluate(0.5), crR.m_dEnd);
        SER(crR.m_pSurface->Evaluate(sUV, 1, 1, TRUE, TRUE, FALSE, sD, FALSE));
        sMid[m] = SmPoint3d(sD[0]);
        sSpeed[m] = ((crR.m_eAlong == SM_SP_U) ? sD[2] : sD[1]).Length();
        sCentre = sCentre + sD[0] * (1.0 / static_cast<double>(n));
      }
    for (size_t m=0; m<n; m++)
      {
        const size_t r = sRound[m];
        const Round & crR = crRounds[r];
        if (!(sSpeed[m] > 0.0))
          { SM_SB_FAIL(rFill.m_pReason, SM_ERR_INVALID_INPUT, _T("SmBlendSetback: a round does not move along its rails at its end")); }
        double dLambda = 0.5 * SmVector3d(sMid[m] - SmPoint3d(sCentre)).Length() / sSpeed[m];
        if (!crR.m_bCornerAbove) { dLambda = -dLambda; }
        const size_t q = sQ[r];
        std::vector< std::vector<sm_SbH> > sH(k + 1, std::vector<sm_SbH>(q + 1));
        double dScale = 1.0;
        for (ULONG j=0; j<=k; j++)
          {
            for (size_t c=0; c<=q; c++)
              {
                const size_t cc = sRev[m] ? q - c : c;
                sH[j][c] = sJets[r][j][cc];
                for (int x=0; x<4; x++) { sH[j][c].a[x] *= dScale; }
              }
            dScale *= dLambda;
          }
        std::vector< std::vector<sm_SbH> > sRows(k + 1, std::vector<sm_SbH>(q + 1));
        for (size_t c=0; c<=q; c++)
          {
            sRows[0][c] = sH[0][c];
            if (k >= 1) { sRows[1][c] = sm_SbAxpy(sRows[0][c], 1.0 / static_cast<double>(N), sH[1][c]); }
            if (k >= 2)
              {
                sm_SbH h = sm_SbAxpy(sm_SbZero(), 2.0, sRows[1][c]);
                h = sm_SbAxpy(h, -1.0, sRows[0][c]);
                sRows[2][c] = sm_SbAxpy(h, 1.0 / static_cast<double>(N * (N - 1)), sH[2][c]);
              }
          }
        // they must meet the faces to order k: their first and last k + 1 columns on the faces' planes
        for (ULONG i=0; i<=k; i++)
          for (ULONG c=0; c<=k; c++)
            for (int e=0; e<2; e++)
              {
                const sm_SbH & crH = sRows[i][e == 0 ? c : q - c];
                const Face & crF = crFaces[sFace[(m + static_cast<size_t>(e)) % n]];
                if (!(crH.a[3] > 0.0))
                  { SM_SB_FAIL(rFill.m_pReason, SM_ERR_INVALID_INPUT, _T("SmBlendSetback: a round's weights must be positive")); }
                const SmPoint3d sP(crH.a[0] / crH.a[3], crH.a[1] / crH.a[3], crH.a[2] / crH.a[3]);
                if (smos_Fabs((sP - crF.m_sPoint).Dot(crF.m_sNormal)) > dTol)
                  { SM_SB_FAIL(rFill.m_pReason, SM_ERR_INVALID_INPUT, _T("SmBlendSetback: a round does not meet its faces to the fill's order at its end (a G1 round is eased first for a G2 fill)")); }
              }
        sRoundRows[m] = sRows;
      }
  }

  // the interior map's degree: the least from 2k + 1 at which its conditions can be met
  for (size_t d=N; d<=9; d++)
    {
      const std::vector< std::pair<int, int> > sMono = sm_SbMonomials(static_cast<int>(d));
      const size_t M = sMono.size();
      const size_t lUnknowns = 3 * M;

      // each edge's rows from F per monomial: D0, D1, D2 (degree d across)
      std::vector< std::vector< std::vector<sm_SbB> > > sEdgeBasis(n, std::vector< std::vector<sm_SbB> >(k + 1, std::vector<sm_SbB>(M)));
      for (size_t m=0; m<n; m++)
        {
          const SmPoint2d & crB0 = sB[m];
          const SmPoint2d & crB1 = sB[(m + 1) % n];
          const sm_SbB sX = { crB0.x, crB1.x }, sY = { crB0.y, crB1.y };
          const SmVector2d sD0 = sSigma[m], sD1 = -sSigma[(m + 1) % n];
          const sm_SbB sDX = { sD0.x, sD1.x }, sDY = { sD0.y, sD1.y };
          const sm_SbB sDXX = sm_SbMul(sDX, sDX), sDXY = sm_SbMul(sDX, sDY), sDYY = sm_SbMul(sDY, sDY);
          for (size_t j=0; j<M; j++)
            {
              const int a = sMono[j].first, b = sMono[j].second;
              sEdgeBasis[m][0][j] = sm_SbAlong(a, b, sX, sY, d);
              // D1: x^a y^b differentiated along the field: a x^(a-1) y^b dx + b x^a y^(b-1) dy
              sm_SbB s1(d + 1, 0.0);
              if (a >= 1)
                {
                  const sm_SbB t = sm_SbElevate(sm_SbMul(sm_SbAlong(a - 1, b, sX, sY, d - 1), sDX), d);
                  for (size_t i=0; i<=d; i++) { s1[i] += a * t[i]; }
                }
              if (b >= 1)
                {
                  const sm_SbB t = sm_SbElevate(sm_SbMul(sm_SbAlong(a, b - 1, sX, sY, d - 1), sDY), d);
                  for (size_t i=0; i<=d; i++) { s1[i] += b * t[i]; }
                }
              sEdgeBasis[m][1][j] = s1;
              if (k >= 2)
                {
                  sm_SbB s2(d + 1, 0.0);
                  if (a >= 2)
                    {
                      const sm_SbB t = sm_SbElevate(sm_SbMul(sm_SbAlong(a - 2, b, sX, sY, d - 2), sDXX), d);
                      for (size_t i=0; i<=d; i++) { s2[i] += a * (a - 1) * t[i]; }
                    }
                  if (a >= 1 && b >= 1)
                    {
                      const sm_SbB t = sm_SbElevate(sm_SbMul(sm_SbAlong(a - 1, b - 1, sX, sY, d - 2), sDXY), d);
                      for (size_t i=0; i<=d; i++) { s2[i] += 2.0 * a * b * t[i]; }
                    }
                  if (b >= 2)
                    {
                      const sm_SbB t = sm_SbElevate(sm_SbMul(sm_SbAlong(a, b - 2, sX, sY, d - 2), sDYY), d);
                      for (size_t i=0; i<=d; i++) { s2[i] += b * (b - 1) * t[i]; }
                    }
                  sEdgeBasis[m][2][j] = s2;
                }
            }
        }
      // each interior part's net per monomial, degree (d, d): corners (0, 0) the centre, (1, 0)
      // the middle of edge q, (1, 1) B_q, (0, 1) the middle of edge q - 1
      std::vector< std::vector<sm_SbT> > sPartBasis(n, std::vector<sm_SbT>(M));
      for (size_t q=0; q<n; q++)
        {
          const SmPoint2d sC(0.0, 0.0);
          const SmPoint2d sMq((sB[q].x + sB[(q + 1) % n].x) * 0.5, (sB[q].y + sB[(q + 1) % n].y) * 0.5);
          const SmPoint2d sMp((sB[q].x + sB[(q + n - 1) % n].x) * 0.5, (sB[q].y + sB[(q + n - 1) % n].y) * 0.5);
          sm_SbT sX, sY;
          sX.p = sX.q = sY.p = sY.q = 1;
          sX.c = { sC.x, sMp.x, sMq.x, sB[q].x };
          sY.c = { sC.y, sMp.y, sMq.y, sB[q].y };
          for (size_t j=0; j<M; j++) { sPartBasis[q][j] = sm_SbOver(sMono[j].first, sMono[j].second, sX, sY, d); }
        }

      // the nets for a given vector of F's coefficients x[c M + j]
      auto nets = [&](const std::vector<double> & crX, std::vector<sm_SbNet> & rNets)
        {
          rNets.assign(2 * n, sm_SbNet());
          for (size_t m=0; m<n; m++)
            {
              const size_t qv = smos_Max(sQ[sRound[m]], d);
              sm_SbNet & rN = rNets[m];
              rN.p = N;
              rN.q = qv;
              rN.P.assign((N + 1) * (qv + 1), sm_SbZero());
              // the round's rows, raised across to qv
              for (ULONG i=0; i<=k; i++)
                for (int c=0; c<4; c++)
                  {
                    sm_SbB s(sRoundRows[m][i].size());
                    for (size_t j=0; j<s.size(); j++) { s[j] = sRoundRows[m][i][j].a[c]; }
                    s = sm_SbElevate(s, qv);
                    for (size_t j=0; j<=qv; j++) { rN.P[i * (qv + 1) + j].a[c] = s[j]; }
                  }
              // F's: E0 = D0, E1 = -D1, E2 = D2 at u = 1
              std::vector< std::vector<sm_SbH> > sE(k + 1, std::vector<sm_SbH>(qv + 1, sm_SbZero()));
              for (ULONG l=0; l<=k; l++)
                for (int c=0; c<3; c++)
                  {
                    sm_SbB s(d + 1, 0.0);
                    for (size_t j=0; j<M; j++)
                      {
                        const double x = crX[static_cast<size_t>(c) * M + j];
                        if (x == 0.0) { continue; }
                        for (size_t i=0; i<=d; i++) { s[i] += x * sEdgeBasis[m][l][j][i]; }
                      }
                    s = sm_SbElevate(s, qv);
                    const double dSign = (l == 1) ? -1.0 : 1.0;
                    for (size_t j=0; j<=qv; j++) { sE[l][j].a[c] = dSign * s[j]; }
                  }
              for (size_t j=0; j<=qv; j++)
                {
                  sE[0][j].a[3] = 1.0;
                  sm_SbH sPN = sE[0][j];
                  rN.P[N * (qv + 1) + j] = sPN;
                  sm_SbH sPN1 = sm_SbAxpy(sPN, -1.0 / static_cast<double>(N), sE[1][j]);
                  rN.P[(N - 1) * (qv + 1) + j] = sPN1;
                  if (k >= 2)
                    {
                      sm_SbH h = sm_SbAxpy(sm_SbZero(), 2.0, sPN1);
                      h = sm_SbAxpy(h, -1.0, sPN);
                      rN.P[(N - 2) * (qv + 1) + j] = sm_SbAxpy(h, 1.0 / static_cast<double>(N * (N - 1)), sE[2][j]);
                    }
                }
            }
          for (size_t q=0; q<n; q++)
            {
              sm_SbNet & rN = rNets[n + q];
              rN.p = rN.q = d;
              rN.P.assign((d + 1) * (d + 1), sm_SbZero());
              for (size_t i=0; i<(d + 1) * (d + 1); i++) { rN.P[i].a[3] = 1.0; }
              for (int c=0; c<3; c++)
                for (size_t j=0; j<M; j++)
                  {
                    const double x = crX[static_cast<size_t>(c) * M + j];
                    if (x == 0.0) { continue; }
                    for (size_t i=0; i<(d + 1) * (d + 1); i++) { rN.P[i].a[c] += x * sPartBasis[q][j].c[i]; }
                  }
            }
        };
      // the conditions: F's rows of each setback patch, their first k + 1 columns on face m's
      // plane and their last k + 1 on face m + 1's
      auto conditions = [&](const std::vector<sm_SbNet> & crNets, std::vector<double> & rC)
        {
          rC.clear();
          for (size_t m=0; m<n; m++)
            {
              const sm_SbNet & crN = crNets[m];
              for (size_t i=N-k; i<=N; i++)
                for (size_t c=0; c<=k; c++)
                  for (int e=0; e<2; e++)
                    {
                      const sm_SbH & crH = crN.P[i * (crN.q + 1) + (e == 0 ? c : crN.q - c)];
                      const Face & crF = crFaces[sFace[(m + static_cast<size_t>(e)) % n]];
                      const SmVector3d & crNorm = crF.m_sNormal;
                      rC.push_back(crNorm.x * crH.a[0] + crNorm.y * crH.a[1] + crNorm.z * crH.a[2] - crNorm.Dot(SmVector3d(crF.m_sPoint)) * crH.a[3]);
                    }
            }
        };
      // the bending energy's terms, each squared and summed
      auto energy = [&](const std::vector<sm_SbNet> & crNets, std::vector<double> & rE)
        {
          rE.clear();
          for (size_t pIdx=0; pIdx<crNets.size(); pIdx++)
            for (int g1=0; g1<8; g1++)
              for (int g2=0; g2<8; g2++)
                {
                  double R[3][3];
                  sm_SbSecond(crNets[pIdx], SM_SB_GX[g1], SM_SB_GX[g2], R);
                  const double w = smos_Sqrt(SM_SB_GW[g1] * SM_SB_GW[g2]);
                  for (int c=0; c<3; c++)
                    {
                      rE.push_back(w * R[0][c]);
                      rE.push_back(w * smos_Sqrt(2.0) * R[1][c]);
                      rE.push_back(w * R[2][c]);
                    }
                }
        };

      // both affine in x: their matrices column by column
      std::vector<double> sX0(lUnknowns, 0.0);
      std::vector<sm_SbNet> sNets;
      std::vector<double> sC0, sE0;
      nets(sX0, sNets);
      conditions(sNets, sC0);
      energy(sNets, sE0);
      const size_t lC = sC0.size(), lE = sE0.size();
      std::vector<double> sAC(lC * lUnknowns), sAE(lE * lUnknowns);
      for (size_t u=0; u<lUnknowns; u++)
        {
          std::vector<double> sX(lUnknowns, 0.0);
          sX[u] = 1.0;
          std::vector<double> sC, sE;
          nets(sX, sNets);
          conditions(sNets, sC);
          energy(sNets, sE);
          for (size_t i=0; i<lC; i++) { sAC[i * lUnknowns + u] = sC[i] - sC0[i]; }
          for (size_t i=0; i<lE; i++) { sAE[i * lUnknowns + u] = sE[i] - sE0[i]; }
        }
      std::vector<double> sRhs(lC);
      double dSize = 0.0;
      for (size_t i=0; i<lC; i++) { sRhs[i] = -sC0[i]; dSize = smos_Max(dSize, smos_Fabs(sC0[i])); }
      std::vector<double> sXp;
      std::vector< std::vector<double> > sZ;
      if (!sm_SbAffineSolutions(sAC, sRhs, lC, lUnknowns, 1.0e-9 * (1.0 + dSize), sXp, sZ) || sZ.empty())
        { continue; }
      // the energy over x = xp + Z y
      const size_t lFree = sZ.size();
      std::vector<double> sM(lE * lFree, 0.0), sR(lE, 0.0);
      for (size_t i=0; i<lE; i++)
        {
          double s = sE0[i];
          for (size_t u=0; u<lUnknowns; u++) { s += sAE[i * lUnknowns + u] * sXp[u]; }
          sR[i] = -s;
          for (size_t f=0; f<lFree; f++)
            {
              double t = 0.0;
              for (size_t u=0; u<lUnknowns; u++) { t += sAE[i * lUnknowns + u] * sZ[f][u]; }
              sM[i * lFree + f] = t;
            }
        }
      std::vector<double> sY;
      if (!sm_SbLeastSquares(sM, sR, lE, lFree, 1.0e-12, sY))
        { continue; }
      std::vector<double> sX(sXp);
      for (size_t f=0; f<lFree; f++)
        for (size_t u=0; u<lUnknowns; u++)
          { sX[u] += sZ[f][u] * sY[f]; }
      nets(sX, sNets);
      std::vector<double> sC;
      conditions(sNets, sC);
      double dWorst = 0.0;
      for (size_t i=0; i<sC.size(); i++) { dWorst = smos_Max(dWorst, smos_Fabs(sC[i])); }
      if (!(dWorst <= 1.0e-3 * dTol))
        { continue; }

      // no fold: every patch's normal keeps to one side, and each setback patch's agrees with
      // its round's side's interior part where they meet
      std::vector<SmVector3d> sMidNormal(2 * n);
      for (size_t pIdx=0; pIdx<2*n; pIdx++)
        {
          SmVector3d sP, sU, sV;
          sm_SbFirst(sNets[pIdx], 0.5, 0.5, sP, sU, sV);
          sMidNormal[pIdx] = sU * sV;
          for (int a=0; a<=8; a++)
            for (int b=0; b<=8; b++)
              {
                sm_SbFirst(sNets[pIdx], a / 8.0, b / 8.0, sP, sU, sV);
                const SmVector3d sNrm = sU * sV;
                if (!(sNrm.Dot(sMidNormal[pIdx]) > 0.0))
                  { SM_SB_FAIL(rFill.m_pReason, SM_ERR, _T("SmBlendSetback: the fill folds")); }
              }
        }
      for (size_t m=0; m<n; m++)
        {
          SmVector3d sP, sU, sV, sQP, sQU, sQV;
          sm_SbFirst(sNets[m], 1.0, 0.25, sP, sU, sV);
          sm_SbFirst(sNets[n + m], 1.0, 0.5, sQP, sQU, sQV);
          if (!((sU * sV).Dot(sQU * sQV) > 0.0) || !(SmVector3d(sP - sQP).Length() <= dTol))
            { SM_SB_FAIL(rFill.m_pReason, SM_ERR, _T("SmBlendSetback: the fill folds where a setback patch meets the interior")); }
        }

      // the surfaces
      std::vector<SmBSplineSurface *> sOut;
      for (size_t pIdx=0; pIdx<2*n; pIdx++)
        {
          const sm_SbNet & crN = sNets[pIdx];
          SmTArray<SmPoint3d> sPts;
          SmTArray<double> sW;
          bool bRational = false;
          for (size_t i=0; i<crN.P.size(); i++)
            {
              const sm_SbH & crH = crN.P[i];
              if (!(crH.a[3] > 0.0))
                {
                  for (size_t o=0; o<sOut.size(); o++) { delete sOut[o]; }
                  SER(SM_ERR);
                }
              sPts.Add(SmPoint3d(crH.a[0] / crH.a[3], crH.a[1] / crH.a[3], crH.a[2] / crH.a[3]));
              sW.Add(crH.a[3]);
              if (crH.a[3] != 1.0) { bRational = true; }
            }
          SmTArray<ULONG> sUM, sVM;
          SmTArray<double> sUK, sVK;
          sUM.Add(static_cast<ULONG>(crN.p + 1)); sUM.Add(static_cast<ULONG>(crN.p + 1));
          sVM.Add(static_cast<ULONG>(crN.q + 1)); sVM.Add(static_cast<ULONG>(crN.q + 1));
          sUK.Add(0.0); sUK.Add(1.0);
          sVK.Add(0.0); sVK.Add(1.0);
          SmBSplineSurface * pSurface = NULL;
          const SmStatus eStatus = SmBSplineSurface::CreateCanonical(crContext, static_cast<ULONG>(crN.p), static_cast<ULONG>(crN.q), sPts,
                                                                     SM_SF_UNSPECIFIED, sUM, sVM, sUK, sVK, SM_KT_UNSPECIFIED,
                                                                     bRational ? &sW : NULL, NULL, pSurface);
          if (eStatus != SM_SUCCESS || pSurface == NULL)
            {
              for (size_t o=0; o<sOut.size(); o++) { delete sOut[o]; }
              SER(eStatus != SM_SUCCESS ? eStatus : SM_ERR);
            }
          sOut.push_back(pSurface);
        }
      rFill.m_vPatches = sOut;
      rFill.m_vRound = sRound;
      rFill.m_vFace = sFace;
      rFill.m_vReversed = sRev;
      rFill.m_lDegree = static_cast<ULONG>(d);
      return SM_SUCCESS;
    }
  SM_SB_FAIL(rFill.m_pReason, SM_ERR, _T("SmBlendSetback: no interior map up to degree 9 meets the setback patches' conditions"));
  return SM_ERR;
}

SmStatus SmBlendSetback::Containment(const Fill & crFill, const std::vector<Face> & crFaces, SmBoolean bConvex, double & rdWorst)
{
  rdWorst = -HUGE_VAL;
  if (crFill.m_vPatches.empty())
    { SER_MSG(SM_ERR_INVALID_INPUT, _T("SmBlendSetback: an empty fill")); }
  const double dSide = bConvex ? 1.0 : -1.0;
  for (size_t i=0; i<crFill.m_vPatches.size(); i++)
    {
      const SmBSplineSurface * pS = crFill.m_vPatches[i]; NER(pS);
      ULONG lUDeg, lVDeg;
      SmTArray<SmPoint3d> sPts;
      SmBSplineSurfaceForm eForm;
      SmTArray<ULONG> sUMult, sVMult;
      SmTArray<double> sUKnots, sVKnots, sWeights;
      SmKnotType eKnotType;
      SER(pS->GetCanonical(lUDeg, lVDeg, sPts, eForm, sUMult, sVMult, sUKnots, sVKnots, eKnotType, sWeights));
      for (ULONG j=0; j<sWeights.GetSize(); j++)
        if (!(sWeights[j] > 0.0))
          { SER_MSG(SM_ERR, _T("SmBlendSetback: a patch with a weight not positive")); }
      for (ULONG j=0; j<sPts.GetSize(); j++)
        for (size_t f=0; f<crFaces.size(); f++)
          { rdWorst = smos_Max(rdWorst, dSide * (sPts[j] - crFaces[f].m_sPoint).Dot(crFaces[f].m_sNormal)); }
    }
  return SM_SUCCESS;
}
