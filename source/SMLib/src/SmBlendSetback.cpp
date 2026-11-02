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
#include <SmFilletExact.h>

#include <algorithm>
#include <cstdlib>
#include <cstdio>
#include <array>
#include <atomic>
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

// the curved fill's work (SmBlendSetback::CurvedCounts)
static std::atomic<unsigned long long> s_lCurvedFills(0), s_lCurvedDegrees(0), s_lCurvedRows(0);

void SmBlendSetback::CurvedCounts(unsigned long long (&raOut)[3], bool bReset)
{
  raOut[0] = s_lCurvedFills.load();
  raOut[1] = s_lCurvedDegrees.load();
  raOut[2] = s_lCurvedRows.load();
  if (bReset)
    {
      s_lCurvedFills = 0;
      s_lCurvedDegrees = 0;
      s_lCurvedRows = 0;
    }
}

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
// across Bezier curve of the j-th derivative, columns from its across start; and its across
// degree. They are the round's geometry's, not its parameter's: its homogeneous coordinates
// divided by the mean weight of its two rails, W = (w_0 + w_q) / 2 (so that W is 1 along it),
// and its parameter along that of the mean arc length of its rails, s with
// ds/dt = (|r_0'| + |r_q'|) / 2 (r_0, r_q its rails): two parameters of one round, polynomial or
// rational, give the same jets. rdDtDs: dt/ds there.
SmStatus sm_SbRoundJets(const SmContext & crContext, const SmBlendSetback::Round & crRound, ULONG k,
                        std::vector< std::vector<sm_SbH> > & rH, size_t & rlQ, double & rdDtDs, const TCHAR *& rpReason)
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

  // the weights' share: G = H / W, W the rails' mean weight
  rdDtDs = 1.0;
  if (k > 2)
    { SER(SM_ERR); } // the jets are taken to second order
  std::vector<double> sW(3, 0.0);
  for (ULONG d=0; d<=k; d++)
    { sW[d] = 0.5 * (rH[d][0].a[3] + rH[d][lQ].a[3]); }
  if (!(sW[0] > 0.0))
    { SM_SB_FAIL(rpReason, SM_ERR_INVALID_INPUT, _T("SmBlendSetback: a round's weights must be positive")); }
  const double dRho0 = 1.0 / sW[0];
  const double dRho1 = -sW[1] / (sW[0] * sW[0]);
  const double dRho2 = (2.0 * sW[1] * sW[1] - sW[0] * sW[2]) / (sW[0] * sW[0] * sW[0]);
  std::vector< std::vector<sm_SbH> > sG(rH);
  for (ULONG j=0; j<=lQ; j++)
    {
      sm_SbH h0 = sm_SbZero(), h1 = sm_SbZero(), h2 = sm_SbZero();
      h0 = sm_SbAxpy(h0, dRho0, rH[0][j]);
      if (k >= 1)
        {
          h1 = sm_SbAxpy(h1, dRho1, rH[0][j]);
          h1 = sm_SbAxpy(h1, dRho0, rH[1][j]);
        }
      if (k >= 2)
        {
          h2 = sm_SbAxpy(h2, dRho2, rH[0][j]);
          h2 = sm_SbAxpy(h2, 2.0 * dRho1, rH[1][j]);
          h2 = sm_SbAxpy(h2, dRho0, rH[2][j]);
        }
      sG[0][j] = h0;
      if (k >= 1) { sG[1][j] = h1; }
      if (k >= 2) { sG[2][j] = h2; }
    }
  // the rails' points and their first two derivatives by t, from the homogeneous ones
  double dSpeed = 0.0, dAccel = 0.0;
  for (int e=0; e<2; e++)
    {
      const ULONG c = (e == 0) ? 0 : lQ;
      const double w0 = sG[0][c].a[3];
      const double w1 = (k >= 1) ? sG[1][c].a[3] : 0.0;
      const double w2 = (k >= 2) ? sG[2][c].a[3] : 0.0;
      SmVector3d P0(sG[0][c].a[0] / w0, sG[0][c].a[1] / w0, sG[0][c].a[2] / w0);
      SmVector3d P1(0.0, 0.0, 0.0), P2(0.0, 0.0, 0.0);
      if (k >= 1)
        { P1 = (SmVector3d(sG[1][c].a[0], sG[1][c].a[1], sG[1][c].a[2]) - P0 * w1) / w0; }
      if (k >= 2)
        { P2 = (SmVector3d(sG[2][c].a[0], sG[2][c].a[1], sG[2][c].a[2]) - P1 * (2.0 * w1) - P0 * w2) / w0; }
      const double dLen = P1.Length();
      if (!(dLen > 0.0))
        { SM_SB_FAIL(rpReason, SM_ERR_INVALID_INPUT, _T("SmBlendSetback: a round's rails do not move along it at its end")); }
      dSpeed += 0.5 * dLen;
      dAccel += 0.5 * P1.Dot(P2) / dLen;
    }
  // the parameter along: s, ds/dt = dSpeed; dt/ds = 1 / dSpeed, d2t/ds2 = -dAccel / dSpeed^3
  const double dTs = 1.0 / dSpeed, dTss = -dAccel / (dSpeed * dSpeed * dSpeed);
  for (ULONG j=0; j<=lQ; j++)
    {
      rH[0][j] = sG[0][j];
      if (k >= 1) { rH[1][j] = sm_SbAxpy(sm_SbZero(), dTs, sG[1][j]); }
      if (k >= 2) { rH[2][j] = sm_SbAxpy(sm_SbAxpy(sm_SbZero(), dTs * dTs, sG[2][j]), dTss, sG[1][j]); }
    }
  rdDtDs = dTs;
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

// The strictly convex quadratic programme min 1/2 y' G y + g' y subject to C y >= b, by the dual
// active-set method (from the unconstrained minimum, adding the most violated constraint each
// step and dropping those whose multipliers would turn negative): G n x n, C m x n, row-major.
// SM_SUCCESS: every constraint holds within dTol (the caller scales C's rows);
// SM_ERR_INVALID_INPUT: G is not positive definite; SM_ERR: no y meets the constraints (or no
// progress within the step limit).
SmStatus sm_SbDualQP(const std::vector<double> & crG, const std::vector<double> & crg, const std::vector<double> & crC,
                     const std::vector<double> & crb, size_t n, size_t m, double dTol, std::vector<double> & rY, bool * pbProven = NULL)
{
  if (pbProven != NULL) { *pbProven = false; }
  // G = L L', J = L^-T (J' G J = I), J[r n + c] its column c
  std::vector<double> L(n * n, 0.0);
  for (size_t j=0; j<n; j++)
    {
      double s = crG[j * n + j];
      for (size_t k=0; k<j; k++) { s -= L[j * n + k] * L[j * n + k]; }
      if (!(s > 0.0))
        { return SM_ERR_INVALID_INPUT; }
      L[j * n + j] = smos_Sqrt(s);
      for (size_t i=j+1; i<n; i++)
        {
          double t = crG[i * n + j];
          for (size_t k=0; k<j; k++) { t -= L[i * n + k] * L[j * n + k]; }
          L[i * n + j] = t / L[j * n + j];
        }
    }
  std::vector<double> J(n * n, 0.0);
  for (size_t c=0; c<n; c++)
    for (size_t i=n; i-- > 0; )
      {
        double s = (i == c) ? 1.0 : 0.0;
        for (size_t k=i+1; k<n; k++) { s -= L[k * n + i] * J[k * n + c]; }
        J[i * n + c] = s / L[i * n + i];
      }
  // the unconstrained minimum, y = -J J' g
  rY.assign(n, 0.0);
  {
    std::vector<double> t(n, 0.0);
    for (size_t c=0; c<n; c++)
      for (size_t r=0; r<n; r++) { t[c] += J[r * n + c] * crg[r]; }
    for (size_t r=0; r<n; r++)
      for (size_t c=0; c<n; c++) { rY[r] -= J[r * n + c] * t[c]; }
  }
  std::vector<size_t> sActive;
  std::vector<double> sU, R(n * n, 0.0), d(n), z(n), r(n);
  auto slack = [&](size_t i)
    {
      double s = -crb[i];
      for (size_t k=0; k<n; k++) { s += crC[i * n + k] * rY[k]; }
      return s;
    };
  auto givens = [](double & a, double & b, double & cs, double & sn)
    {
      const double h = std::hypot(a, b);
      if (h == 0.0) { cs = 1.0; sn = 0.0; return; }
      cs = a / h; sn = b / h; a = h; b = 0.0;
    };
  const size_t lMaxSteps = 50 * (m + n) + 1000;
  for (size_t lStep=0; lStep<lMaxSteps; lStep++)
    {
      size_t p = m;
      double dWorst = -dTol;
      for (size_t i=0; i<m; i++)
        {
          const double s = slack(i);
          if (s < dWorst) { dWorst = s; p = i; }
        }
      if (p == m)
        { return SM_SUCCESS; }
      double dUp = 0.0;
      for (;;)
        {
          const size_t q = sActive.size();
          for (size_t c=0; c<n; c++)
            {
              double s = 0.0;
              for (size_t k=0; k<n; k++) { s += J[k * n + c] * crC[p * n + k]; }
              d[c] = s;
            }
          // the step in y, z = J2 d2, and in the multipliers, r = R^-1 d1
          double dZ = 0.0;
          for (size_t k=0; k<n; k++)
            {
              double s = 0.0;
              for (size_t c=q; c<n; c++) { s += J[k * n + c] * d[c]; }
              z[k] = s;
              dZ = smos_Max(dZ, smos_Fabs(s));
            }
          for (size_t i=q; i-- > 0; )
            {
              double s = d[i];
              for (size_t j=i+1; j<q; j++) { s -= R[i * n + j] * r[j]; }
              r[i] = s / R[i * n + i];
            }
          double t1 = HUGE_VAL;
          size_t lDrop = q;
          for (size_t k=0; k<q; k++)
            if (r[k] > 0.0 && sU[k] / r[k] < t1) { t1 = sU[k] / r[k]; lDrop = k; }
          double zn = 0.0;
          for (size_t k=0; k<n; k++) { zn += z[k] * crC[p * n + k]; }
          const double t2 = (dZ > 1.0e-14 && zn > 0.0) ? -slack(p) / zn : HUGE_VAL;
          const double t = smos_Min(t1, t2);
          if (t == HUGE_VAL)
            {
              // no step: n_p = N_A r (z = 0) with every r_k <= 0. For any feasible y,
              // n_p . y = sum r_k n_k . y <= sum r_k b_k: proven infeasible when b_p exceeds that
              // (checked here on the rows themselves, not trusted from the factors)
              if (!(dZ <= 1.0e-14))
                { return SM_ERR; }
              double dResid = 0.0, dScale = 1.0, dBound = 0.0;
              for (size_t k=0; k<n; k++)
                {
                  double s = crC[p * n + k];
                  for (size_t a=0; a<q; a++) { s -= r[a] * crC[sActive[a] * n + k]; }
                  dResid = smos_Max(dResid, smos_Fabs(s));
                }
              for (size_t a=0; a<q; a++)
                {
                  dBound += r[a] * crb[sActive[a]];
                  dScale = smos_Max(dScale, smos_Fabs(r[a]));
                }
              if (pbProven != NULL) { *pbProven = dResid <= 1.0e-9 * dScale && crb[p] - dBound > 1.0e-9 * (dScale + smos_Fabs(crb[p])); }
              return SM_ERR;
            }
          if (t2 != HUGE_VAL)
            { for (size_t k=0; k<n; k++) { rY[k] += t * z[k]; } }
          for (size_t k=0; k<q; k++) { sU[k] -= t * r[k]; }
          dUp += t;
          if (t2 != HUGE_VAL && t == t2)
            {
              // p joins: d[q + 1 ..] rotated into d[q], J's columns alike
              for (size_t j=n-1; j>q; j--)
                {
                  double cs, sn;
                  givens(d[j - 1], d[j], cs, sn);
                  for (size_t k=0; k<n; k++)
                    {
                      const double a = J[k * n + j - 1], b = J[k * n + j];
                      J[k * n + j - 1] = cs * a + sn * b;
                      J[k * n + j] = -sn * a + cs * b;
                    }
                }
              if (!(smos_Fabs(d[q]) > 1.0e-14))
                { return SM_ERR; }
              for (size_t i=0; i<=q; i++) { R[i * n + q] = d[i]; }
              sActive.push_back(p);
              sU.push_back(dUp);
              break;
            }
          // lDrop leaves: its column out of R, R made triangular again by rotations of its rows
          // (J's columns alike)
          sActive.erase(sActive.begin() + static_cast<long>(lDrop));
          sU.erase(sU.begin() + static_cast<long>(lDrop));
          for (size_t i=0; i<q; i++)
            {
              for (size_t j=lDrop; j+1<q; j++) { R[i * n + j] = R[i * n + j + 1]; }
              R[i * n + q - 1] = 0.0;
            }
          for (size_t j=lDrop; j+1<q; j++)
            {
              double cs, sn;
              double a = R[j * n + j], b = R[(j + 1) * n + j];
              givens(a, b, cs, sn);
              R[j * n + j] = a;
              R[(j + 1) * n + j] = 0.0;
              for (size_t k=j+1; k+1<q; k++)
                {
                  const double x1 = R[j * n + k], x2 = R[(j + 1) * n + k];
                  R[j * n + k] = cs * x1 + sn * x2;
                  R[(j + 1) * n + k] = -sn * x1 + cs * x2;
                }
              for (size_t k=0; k<n; k++)
                {
                  const double x1 = J[k * n + j], x2 = J[k * n + j + 1];
                  J[k * n + j] = cs * x1 + sn * x2;
                  J[k * n + j + 1] = -sn * x1 + cs * x2;
                }
            }
        }
    }
  return SM_ERR;
}

// sm_SbDualQP over many rows, few of which bind: solved on a working set of rows, those violated
// at its solution added (the most violated first, up to a block at a time) until none is. Its
// answer is the full programme's: a solution meets every row; a set of rows with no solution is
// one the full set has none of either (*pbProven as sm_SbDualQP's).
SmStatus sm_SbDualQPRows(const std::vector<double> & crG, const std::vector<double> & crg, const std::vector<double> & crC,
                         const std::vector<double> & crb, size_t n, size_t m, double dTol, std::vector<double> & rY, bool * pbProven = NULL)
{
  if (pbProven != NULL) { *pbProven = false; }
  const size_t lBlock = 2000;
  std::vector<char> sIn(m, 0);
  std::vector<size_t> sSet;
  // y at the unconstrained minimum: the empty set's solution
  SER(sm_SbDualQP(crG, crg, std::vector<double>(), std::vector<double>(), n, 0, dTol, rY));
  for (int iRound=0; iRound<200; iRound++)
    {
      std::vector< std::pair<double, size_t> > sViolated;
      for (size_t i=0; i<m; i++)
        {
          if (sIn[i]) { continue; }
          double s = -crb[i];
          for (size_t k=0; k<n; k++) { s += crC[i * n + k] * rY[k]; }
          if (s < -dTol) { sViolated.push_back(std::make_pair(s, i)); }
        }
      if (sViolated.empty())
        { return SM_SUCCESS; }
      const size_t lTake = smos_Min(lBlock, sViolated.size());
      std::partial_sort(sViolated.begin(), sViolated.begin() + static_cast<long>(lTake), sViolated.end());
      for (size_t t=0; t<lTake; t++)
        {
          sIn[sViolated[t].second] = 1;
          sSet.push_back(sViolated[t].second);
        }
      std::vector<double> sC(sSet.size() * n), sb(sSet.size());
      for (size_t t=0; t<sSet.size(); t++)
        {
          for (size_t k=0; k<n; k++) { sC[t * n + k] = crC[sSet[t] * n + k]; }
          sb[t] = crb[sSet[t]];
        }
      const SmStatus eStatus = sm_SbDualQP(crG, crg, sC, sb, n, sSet.size(), dTol, rY, pbProven);
      if (eStatus != SM_SUCCESS)
        { return eStatus; }
    }
  return SM_ERR;
}

// Bivariate Bernstein polynomials (sm_SbT) on [0, 1]^2: sums, derivatives, halves.
sm_SbT sm_SbTMake(size_t p, size_t q)
{
  sm_SbT s;
  s.p = p;
  s.q = q;
  s.c.assign((p + 1) * (q + 1), 0.0);
  return s;
}

// crA + s crB, both raised to the larger degrees
sm_SbT sm_SbTAxpy(const sm_SbT & crA, double s, const sm_SbT & crB)
{
  const size_t p = smos_Max(crA.p, crB.p), q = smos_Max(crA.q, crB.q);
  sm_SbT sA = (crA.p == p && crA.q == q) ? crA : sm_SbTElevate(crA, p, q);
  const sm_SbT sB = (crB.p == p && crB.q == q) ? crB : sm_SbTElevate(crB, p, q);
  for (size_t i=0; i<sA.c.size(); i++) { sA.c[i] += s * sB.c[i]; }
  return sA;
}

sm_SbT sm_SbTDu(const sm_SbT & crA)
{
  if (crA.p == 0) { return sm_SbTMake(0, crA.q); }
  sm_SbT s = sm_SbTMake(crA.p - 1, crA.q);
  for (size_t i=0; i<crA.p; i++)
    for (size_t j=0; j<=crA.q; j++)
      { s.c[i * (crA.q + 1) + j] = static_cast<double>(crA.p) * (crA.c[(i + 1) * (crA.q + 1) + j] - crA.c[i * (crA.q + 1) + j]); }
  return s;
}

sm_SbT sm_SbTDv(const sm_SbT & crA)
{
  if (crA.q == 0) { return sm_SbTMake(crA.p, 0); }
  sm_SbT s = sm_SbTMake(crA.p, crA.q - 1);
  for (size_t i=0; i<=crA.p; i++)
    for (size_t j=0; j<crA.q; j++)
      { s.c[i * crA.q + j] = static_cast<double>(crA.q) * (crA.c[i * (crA.q + 1) + j + 1] - crA.c[i * (crA.q + 1) + j]); }
  return s;
}

// crA's halves at u = 1/2 (iDir 0) or v = 1/2 (iDir 1), by de Casteljau.
void sm_SbTSplit(const sm_SbT & crA, int iDir, sm_SbT & rLo, sm_SbT & rHi)
{
  rLo = crA;
  rHi = crA;
  const size_t n = (iDir == 0) ? crA.p : crA.q, lOther = (iDir == 0) ? crA.q : crA.p;
  std::vector<double> w(n + 1);
  for (size_t k=0; k<=lOther; k++)
    {
      auto idx = [&](size_t i) { return (iDir == 0) ? i * (crA.q + 1) + k : k * (crA.q + 1) + i; };
      for (size_t i=0; i<=n; i++) { w[i] = crA.c[idx(i)]; }
      for (size_t r=0; r<=n; r++)
        {
          rLo.c[idx(r)] = w[0];
          rHi.c[idx(n - r)] = w[n - r];
          for (size_t i=0; i+r<n; i++) { w[i] = 0.5 * (w[i] + w[i + 1]); }
        }
    }
}

// The pieces of crA after lLevels halvings each way, in a fixed order.
void sm_SbTPieces(const sm_SbT & crA, size_t lLevels, std::vector<sm_SbT> & rPieces)
{
  rPieces.assign(1, crA);
  for (size_t l=0; l<lLevels; l++)
    {
      std::vector<sm_SbT> sNext;
      for (size_t k=0; k<rPieces.size(); k++)
        {
          sm_SbT a, b;
          sm_SbTSplit(rPieces[k], 0, a, b);
          sm_SbT a0, a1, b0, b1;
          sm_SbTSplit(a, 1, a0, a1);
          sm_SbTSplit(b, 1, b0, b1);
          sNext.push_back(a0); sNext.push_back(a1); sNext.push_back(b0); sNext.push_back(b1);
        }
      rPieces.swap(sNext);
    }
}

// crA crB, as sm_SbTMul but the binomials taken once: the coefficients scaled by them,
// convolved, and unscaled.
sm_SbT sm_SbTMulFast(const sm_SbT & crA, const sm_SbT & crB)
{
  sm_SbT sC = sm_SbTMake(crA.p + crB.p, crA.q + crB.q);
  std::vector<double> sA(crA.c.size()), sB(crB.c.size());
  for (size_t i=0; i<=crA.p; i++)
    for (size_t j=0; j<=crA.q; j++)
      { sA[i * (crA.q + 1) + j] = crA.c[i * (crA.q + 1) + j] * sm_Binomial(static_cast<ULONG>(crA.p), static_cast<ULONG>(i)) * sm_Binomial(static_cast<ULONG>(crA.q), static_cast<ULONG>(j)); }
  for (size_t i=0; i<=crB.p; i++)
    for (size_t j=0; j<=crB.q; j++)
      { sB[i * (crB.q + 1) + j] = crB.c[i * (crB.q + 1) + j] * sm_Binomial(static_cast<ULONG>(crB.p), static_cast<ULONG>(i)) * sm_Binomial(static_cast<ULONG>(crB.q), static_cast<ULONG>(j)); }
  for (size_t i=0; i<=crA.p; i++)
    for (size_t j=0; j<=crA.q; j++)
      {
        const double a = sA[i * (crA.q + 1) + j];
        if (a == 0.0) { continue; }
        for (size_t k=0; k<=crB.p; k++)
          {
            double * pRow = &sC.c[(i + k) * (sC.q + 1) + j];
            const double * pB = &sB[k * (crB.q + 1)];
            for (size_t l=0; l<=crB.q; l++) { pRow[l] += a * pB[l]; }
          }
      }
  for (size_t i=0; i<=sC.p; i++)
    for (size_t j=0; j<=sC.q; j++)
      { sC.c[i * (sC.q + 1) + j] /= sm_Binomial(static_cast<ULONG>(sC.p), static_cast<ULONG>(i)) * sm_Binomial(static_cast<ULONG>(sC.q), static_cast<ULONG>(j)); }
  return sC;
}

// The determinant of a 3 x 3 matrix of polynomials.
sm_SbT sm_SbTDet3(const sm_SbT * const (&crM)[3][3])
{
  sm_SbT s = sm_SbTMake(0, 0);
  for (int j=0; j<3; j++)
    {
      const int j1 = (j + 1) % 3, j2 = (j + 2) % 3;
      const sm_SbT sMinor = sm_SbTAxpy(sm_SbTMulFast(*crM[1][j1], *crM[2][j2]), -1.0, sm_SbTMulFast(*crM[1][j2], *crM[2][j1]));
      s = sm_SbTAxpy(s, 1.0, sm_SbTMulFast(*crM[0][j], sMinor));
    }
  return s;
}

// Whether crA >= -dEps on the whole square, certified: 1 (every coefficient of each piece >= -dEps,
// halving each way up to lDepth times), 0 (a corner of a piece, a value of crA, below -dEps), -1
// (undecided at that depth).
int sm_SbTNonneg(const sm_SbT & crA, double dEps, size_t lDepth)
{
  const double dCorner[4] = { crA.c[0], crA.c[crA.q], crA.c[crA.p * (crA.q + 1)], crA.c.back() };
  for (int k=0; k<4; k++)
    if (dCorner[k] < -dEps) { return 0; }
  if (*std::min_element(crA.c.begin(), crA.c.end()) >= -dEps)
    { return 1; }
  if (lDepth == 0)
    { return -1; }
  sm_SbT a, b, a0, a1, b0, b1;
  sm_SbTSplit(crA, 0, a, b);
  sm_SbTSplit(a, 1, a0, a1);
  sm_SbTSplit(b, 1, b0, b1);
  int iResult = 1;
  const sm_SbT * pPieces[4] = { &a0, &a1, &b0, &b1 };
  for (int k=0; k<4; k++)
    {
      const int r = sm_SbTNonneg(*pPieces[k], dEps, lDepth - 1);
      if (r == 0) { return 0; }
      iResult = smos_Min(iResult, r);
    }
  return iResult;
}

double sm_SbTScale(const sm_SbT & crA)
{
  double d = 0.0;
  for (size_t i=0; i<crA.c.size(); i++) { d = smos_Max(d, smos_Fabs(crA.c[i])); }
  return d;
}

// A homogeneous net's (x w, y w, z w, w) curvature determinants: l = det(H, Hu, Hv, Huu),
// m = det(H, Hu, Hv, Huv), n = det(H, Hu, Hv, Hvv) (each -w^4 (Pu x Pv) . P.. for P = H / w) and
// o = det(H, Hu, Hv, (e, 0)) (-w^3 (Pu x Pv) . e).
void sm_SbFundamentals(const sm_SbNet & crN, const SmVector3d & crE, sm_SbT (&rLMNO)[4])
{
  sm_SbT h[4], hu[4], hv[4], huu[4], huv[4], hvv[4], one[4];
  for (int c=0; c<4; c++)
    {
      h[c] = sm_SbTMake(crN.p, crN.q);
      for (size_t i=0; i<crN.P.size(); i++) { h[c].c[i] = crN.P[i].a[c]; }
      hu[c] = sm_SbTDu(h[c]);
      hv[c] = sm_SbTDv(h[c]);
      huu[c] = sm_SbTDu(hu[c]);
      huv[c] = sm_SbTDv(hu[c]);
      hvv[c] = sm_SbTDv(hv[c]);
      one[c] = sm_SbTMake(0, 0);
      one[c].c[0] = (c == 0) ? crE.x : (c == 1) ? crE.y : (c == 2) ? crE.z : 0.0;
    }
  const sm_SbT * pLast[4] = { huu, huv, hvv, one };
  const sm_SbT * pRows[3] = { h, hu, hv };
  for (int k=0; k<4; k++)
    {
      // along the last row: det = sum_c (-1)^(3 + c) last[c] minor(3, c)
      sm_SbT s = sm_SbTMake(0, 0);
      for (int c=0; c<4; c++)
        {
          const sm_SbT * pM[3][3];
          for (int r=0; r<3; r++)
            {
              int cc = 0;
              for (int col=0; col<4; col++)
                if (col != c) { pM[r][cc++] = &pRows[r][col]; }
            }
          const sm_SbT sMinor = sm_SbTDet3(pM);
          s = sm_SbTAxpy(s, ((3 + c) % 2 == 0) ? 1.0 : -1.0, sm_SbTMulFast(pLast[k][c], sMinor));
        }
      rLMNO[k] = s;
    }
}

// A patch whose control points move only along one direction e, linearly in the unknowns: its
// net is m_sBase + sum_j y_j Z_j (e, 0), Z_j (m_vZ, of the net's degrees) for the unknowns
// m_vIndex. Its projection along e (and its weights) never moves.
struct sm_SbAffinePatch
{
  sm_SbNet              m_sBase;
  std::vector<size_t>   m_vIndex;
  std::vector<sm_SbT>   m_vZ;
};

sm_SbNet sm_SbAffineNet(const sm_SbAffinePatch & crP, const std::vector<double> & crY, const SmVector3d & crE)
{
  sm_SbNet sN = crP.m_sBase;
  for (size_t k=0; k<crP.m_vIndex.size(); k++)
    {
      const double y = crY[crP.m_vIndex[k]];
      if (y == 0.0) { continue; }
      for (size_t i=0; i<sN.P.size(); i++)
        {
          const double z = y * crP.m_vZ[k].c[i];
          sN.P[i].a[0] += z * crE.x;
          sN.P[i].a[1] += z * crE.y;
          sN.P[i].a[2] += z * crE.z;
        }
    }
  return sN;
}

// A net's patch as a rational (or polynomial, all weights 1) Bezier surface over [0, 1]^2.
SmStatus sm_SbSurface(const SmContext & crContext, const sm_SbNet & crN, SmBSplineSurface *& rpSurface)
{
  rpSurface = NULL;
  SmTArray<SmPoint3d> sPts;
  SmTArray<double> sW;
  bool bRational = false;
  for (size_t i=0; i<crN.P.size(); i++)
    {
      const sm_SbH & crH = crN.P[i];
      if (!(crH.a[3] > 0.0))
        { return SM_ERR; }
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
  return SmBSplineSurface::CreateCanonical(crContext, static_cast<ULONG>(crN.p), static_cast<ULONG>(crN.q), sPts, SM_SF_UNSPECIFIED, sUM, sVM, sUK, sVK,
                                           SM_KT_UNSPECIFIED, bRational ? &sW : NULL, NULL, rpSurface);
}

// The determinant of a 4 x 4 matrix of polynomials (rows crM[r]), along its last row.
sm_SbT sm_SbTDet4(const sm_SbT * const (&crM)[4][4])
{
  sm_SbT s = sm_SbTMake(0, 0);
  for (int c=0; c<4; c++)
    {
      const sm_SbT * pM[3][3];
      for (int r=0; r<3; r++)
        {
          int cc = 0;
          for (int col=0; col<4; col++)
            if (col != c) { pM[r][cc++] = crM[r][col]; }
        }
      s = sm_SbTAxpy(s, ((3 + c) % 2 == 0) ? 1.0 : -1.0, sm_SbTMulFast(*crM[3][c], sm_SbTDet3(pM)));
    }
  return s;
}

// An affine patch's curvature determinants l, m, n (k = 0, 1, 2) as affine functions of the
// unknowns: rQ0[k] at y = 0, rQj[k][t] the share of its unknown crP.m_vIndex[t]; and o (fixed).
// With H = H0 + Z (e, 0) and E = (e, 0), det(H, Hu, Hv, X) = det(H0, H0u, H0v, X0) +
// sum_r Z_r det(.., E in row r, ..) (two rows along E give 0): affine in Z.
void sm_SbAffineFundamentals(const sm_SbAffinePatch & crP, const SmVector3d & crE, sm_SbT (&rQ0)[3], std::vector<sm_SbT> (&rQj)[3], sm_SbT & rO)
{
  sm_SbT h[4], hu[4], hv[4], hx[3][4], one[4];
  for (int c=0; c<4; c++)
    {
      h[c] = sm_SbTMake(crP.m_sBase.p, crP.m_sBase.q);
      for (size_t i=0; i<crP.m_sBase.P.size(); i++) { h[c].c[i] = crP.m_sBase.P[i].a[c]; }
      hu[c] = sm_SbTDu(h[c]);
      hv[c] = sm_SbTDv(h[c]);
      hx[0][c] = sm_SbTDu(hu[c]);
      hx[1][c] = sm_SbTDv(hu[c]);
      hx[2][c] = sm_SbTDv(hv[c]);
      one[c] = sm_SbTMake(0, 0);
      one[c].c[0] = (c == 0) ? crE.x : (c == 1) ? crE.y : (c == 2) ? crE.z : 0.0;
    }
  const sm_SbT * pO[4][4];
  for (int c=0; c<4; c++) { pO[0][c] = &h[c]; pO[1][c] = &hu[c]; pO[2][c] = &hv[c]; pO[3][c] = &one[c]; }
  rO = sm_SbTDet4(pO);
  for (int k=0; k<3; k++)
    {
      const sm_SbT * pRows[4][4];
      for (int c=0; c<4; c++) { pRows[0][c] = &h[c]; pRows[1][c] = &hu[c]; pRows[2][c] = &hv[c]; pRows[3][c] = &hx[k][c]; }
      rQ0[k] = sm_SbTDet4(pRows);
      // the cofactors: row r replaced by E (row 3's is o)
      sm_SbT sCof[4];
      for (int r=0; r<3; r++)
        {
          const sm_SbT * pR[4][4];
          for (int rr=0; rr<4; rr++)
            for (int c=0; c<4; c++) { pR[rr][c] = (rr == r) ? &one[c] : pRows[rr][c]; }
          sCof[r] = sm_SbTDet4(pR);
        }
      sCof[3] = rO;
      rQj[k].clear();
      for (size_t t=0; t<crP.m_vZ.size(); t++)
        {
          const sm_SbT & crZ = crP.m_vZ[t];
          const sm_SbT sZu = sm_SbTDu(crZ), sZv = sm_SbTDv(crZ);
          const sm_SbT sZx = (k == 0) ? sm_SbTDu(sZu) : (k == 1) ? sm_SbTDv(sZu) : sm_SbTDv(sZv);
          sm_SbT s = sm_SbTMulFast(crZ, sCof[0]);
          s = sm_SbTAxpy(s, 1.0, sm_SbTMulFast(sZu, sCof[1]));
          s = sm_SbTAxpy(s, 1.0, sm_SbTMulFast(sZv, sCof[2]));
          s = sm_SbTAxpy(s, 1.0, sm_SbTMulFast(sZx, sCof[3]));
          rQj[k].push_back(s);
        }
    }
}

// The fill's shape certified from its nets: on every patch o = det(H, Hu, Hv, (e, 0)) of one sign
// s and s l, s n, l n - m^2 >= 0 (each to 1e-12 of its coefficients' scale), from Bernstein
// coefficients halving up to lDepth times. rlPatch: the first patch not certified.
bool sm_SbCertifiedConvex(const std::vector<sm_SbNet> & crNets, const SmVector3d & crE, size_t lDepth, size_t & rlPatch)
{
  for (size_t p=0; p<crNets.size(); p++)
    {
      rlPatch = p;
      sm_SbT sF[4];
      sm_SbFundamentals(crNets[p], crE, sF);
      const sm_SbT & crO = sF[3];
      const double dO = sm_SbTScale(crO);
      sm_SbT sNegO = crO;
      for (size_t i=0; i<sNegO.c.size(); i++) { sNegO.c[i] = -sNegO.c[i]; }
      const double s = (sm_SbTNonneg(crO, 1.0e-12 * dO, lDepth) == 1) ? 1.0 : (sm_SbTNonneg(sNegO, 1.0e-12 * dO, lDepth) == 1) ? -1.0 : 0.0;
      if (s == 0.0)
        { return false; }
      sm_SbT sL = sF[0], sN = sF[2];
      for (size_t i=0; i<sL.c.size(); i++) { sL.c[i] *= s; }
      for (size_t i=0; i<sN.c.size(); i++) { sN.c[i] *= s; }
      const sm_SbT sDisc = sm_SbTAxpy(sm_SbTMulFast(sF[0], sF[2]), -1.0, sm_SbTMulFast(sF[1], sF[1]));
      const sm_SbT * pChecks[3] = { &sL, &sN, &sDisc };
      for (int k=0; k<3; k++)
        if (sm_SbTNonneg(*pChecks[k], 1.0e-12 * sm_SbTScale(*pChecks[k]), lDepth) != 1)
          { return false; }
    }
  return true;
}

} // namespace

/*******************************************************************//**
PURPOSE: The G1 fill of a corner as a graph over the corner's axis.

NOTES: The fill keeps the setback split's layout (SmBlendSetback.h) and its
    G1 by construction (the rounds' rows, Hahn's connecting map, the faces'
    columns), but every patch is a graph over the plane square to the corner's
    axis e (into the material at a convex corner, into the void at a concave
    one): the interior map is F = xi_1 a_1 + xi_2 a_2 + h(xi) e over the
    projection of the base points (a_1, a_2, e orthonormal, h a polynomial in
    xi), and each setback patch's projection along e is fixed by its round's
    rows, its speed off the round and F's projection (the spring curves fix
    them, below); only heights along e
    are free. Then each patch's l = det(H, Hu, Hv, Huu), m and n (its second
    fundamental form times -w^4 |Pu x Pv|) are affine in the heights, and the
    side its normal faces, o = det(H, Hu, Hv, (e, 0)), depends on the
    projection only: certified of one sign, or the fill folds and is refused.
    The body stays convex where s [[l, m], [m, n]] (s o's sign) is positive
    semidefinite on every patch. A linear sufficient condition: every
    Bernstein coefficient matrix C of a patch's pieces (halved up to
    SM_SB_GRAPH_LEVELS times each way) has d_k' C d_k >= 0 and
    d_k' C d_k+1 >= 0 for directions d_0 .. d_K over a half turn (C is then
    nonnegative on each sector between them); an interior part's directions
    include sigma's at its base point, along which F's Hessian is forced
    singular by the faces' conditions. The heights are the least thin-plate
    energy meeting the faces' conditions (equalities) and these
    (inequalities), by the dual active-set method, at degree 9 first (proven
    infeasible there, the conditions are at every lower degree: the search
    stops), then up from the least degree, each result certified from the
    patches' own Bernstein coefficients (o one-signed, s l, s n and
    l n - m^2 >= 0, subdividing) before it is taken. Where none is the
    corner is refused (m_bUncertified): its setbacks are too short for a
    convex fill (CertifiedScale finds the scale that is not).

    The spring curves (geometry only): on face m, with P and Q its two rails'
    ends and X where the rails cross, the spring curve is the parabola
    (P, X, Q), tangent to both rails at their ends (the spring conic of
    weight 1: polynomial whatever the angle between the rails, and inside the
    face's flat corner region, which a convex body's face is). Its middle
    (P + 2 X + Q) / 4 is the base point; the two setback patches' sides on the
    face are its halves as cubics: the round's speed off its end at P is a
    third of |X - P| (the halves' second control points), linear across the
    round from one rail to the other, and the cross field at the base point is
    (P - Q) / 2, the parabola's tangent there.
***********************************************************************/
namespace
{

const size_t SM_SB_GRAPH_DEPTH = 6;       // the certificate's halvings
const size_t SM_SB_GRAPH_LEVELS = 1;      // the sector conditions' halvings
const size_t SM_SB_GRAPH_DIRECTIONS = 16; // their directions over a half turn

// The fill's surfaces from its nets.
SmStatus sm_SbFillFrom(const SmContext & crContext, const std::vector<sm_SbNet> & crNets, const std::vector<ULONG> & crRound, const std::vector<ULONG> & crFace,
                       const std::vector<SmBoolean> & crRev, size_t lDegree, SmBlendSetback::Fill & rFill)
{
  std::vector<SmBSplineSurface *> sOut;
  for (size_t p=0; p<crNets.size(); p++)
    {
      SmBSplineSurface * pSurface = NULL;
      const SmStatus eStatus = sm_SbSurface(crContext, crNets[p], pSurface);
      if (eStatus != SM_SUCCESS || pSurface == NULL)
        {
          for (size_t o=0; o<sOut.size(); o++) { delete sOut[o]; }
          SER(eStatus != SM_SUCCESS ? eStatus : SM_ERR);
        }
      sOut.push_back(pSurface);
    }
  rFill.m_vPatches = sOut;
  rFill.m_vRound = crRound;
  rFill.m_vFace = crFace;
  rFill.m_vReversed = crRev;
  rFill.m_lDegree = static_cast<ULONG>(lDegree);
  return SM_SUCCESS;
}

SmStatus sm_SbGraphG1(const SmContext & crContext, const std::vector<SmBlendSetback::Face> & crWorldFaces,
                       const std::vector<SmBlendSetback::Round> & crRounds, const std::vector< std::vector< std::vector<sm_SbH> > > & crWorldJets,
                       const std::vector<size_t> & crQ, const std::vector<ULONG> & crRound, const std::vector<ULONG> & crFace,
                       const std::vector<SmBoolean> & crRev, double dTol, SmBlendSetback::Fill & rFill,
                       SmBlendSetback::Probe * pProbe = NULL)
{
  const size_t n = crRound.size();
  const size_t Nc = 3;                    // the cubic across: rows 0, 1 the round's, 2, 3 F's
  // the fill is built and certified about the corner's vertex O, the faces' common point (least
  // squares): its rounding that of the corner's size, not of its place, and a face through O an
  // axis plane exactly where it is one in the world
  SmVector3d sO(0.0, 0.0, 0.0);
  {
    double A[3][3] = {{0.0}}, b[3] = {0.0, 0.0, 0.0};
    for (size_t f=0; f<crWorldFaces.size(); f++)
      {
        const SmVector3d & crN = crWorldFaces[f].m_sNormal;
        const double c[3] = { crN.x, crN.y, crN.z }, d = crN.Dot(SmVector3d(crWorldFaces[f].m_sPoint));
        for (int i=0; i<3; i++)
          {
            for (int j=0; j<3; j++) { A[i][j] += c[i] * c[j]; }
            b[i] += c[i] * d;
          }
      }
    const double dDet = A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) - A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0])
                        + A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);
    if (!(smos_Fabs(dDet) > 1.0e-12))
      { SM_SB_FAIL(rFill.m_pReason, SM_ERR_INVALID_INPUT, _T("SmBlendSetback: the corner's faces meet in no one point")); }
    double x[3];
    for (int k=0; k<3; k++)
      {
        double M[3][3];
        for (int i=0; i<3; i++)
          for (int j=0; j<3; j++) { M[i][j] = (j == k) ? b[i] : A[i][j]; }
        x[k] = (M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
                + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0])) / dDet;
      }
    sO = SmVector3d(x[0], x[1], x[2]);
  }
  std::vector< std::vector< std::vector<sm_SbH> > > crJets(crWorldJets);
  for (size_t r=0; r<crJets.size(); r++)
    for (size_t k=0; k<crJets[r].size(); k++)
      for (size_t c=0; c<crJets[r][k].size(); c++)
        {
          sm_SbH & rH = crJets[r][k][c];
          rH.a[0] -= rH.a[3] * sO.x;
          rH.a[1] -= rH.a[3] * sO.y;
          rH.a[2] -= rH.a[3] * sO.z;
        }
  std::vector<SmBlendSetback::Face> crFaces(crWorldFaces);
  for (size_t f=0; f<crFaces.size(); f++) { crFaces[f].m_sPoint = SmPoint3d(SmVector3d(crFaces[f].m_sPoint) - sO); }
  // the nets moved back, the faces' points snapped onto their planes there
  auto toWorld = [&](std::vector<sm_SbNet> & rNets)
    {
      for (size_t p=0; p<rNets.size(); p++)
        for (size_t i=0; i<rNets[p].P.size(); i++)
          {
            sm_SbH & rH = rNets[p].P[i];
            rH.a[0] += rH.a[3] * sO.x;
            rH.a[1] += rH.a[3] * sO.y;
            rH.a[2] += rH.a[3] * sO.z;
          }
      for (size_t m=0; m<n; m++)
        {
          sm_SbNet & rN = rNets[m];
          const size_t qv = rN.q;
          for (size_t i=2; i<=Nc; i++)
            for (size_t c=0; c<2; c++)
              for (int e=0; e<2; e++)
                {
                  sm_SbH & rH = rN.P[i * (qv + 1) + ((e == 0) ? c : qv - c)];
                  const SmBlendSetback::Face & crF = crWorldFaces[crFace[(m + static_cast<size_t>(e)) % n]];
                  const SmVector3d & crN = crF.m_sNormal;
                  const double dOff = crN.x * rH.a[0] + crN.y * rH.a[1] + crN.z * rH.a[2] - crN.Dot(SmVector3d(crF.m_sPoint)) * rH.a[3];
                  rH.a[0] -= dOff * crN.x;
                  rH.a[1] -= dOff * crN.y;
                  rH.a[2] -= dOff * crN.z;
                  for (int x=0; x<3; x++)
                    {
                      const double dNx = (x == 0) ? crN.x : (x == 1) ? crN.y : crN.z;
                      const double dPx = (x == 0) ? crF.m_sPoint.x : (x == 1) ? crF.m_sPoint.y : crF.m_sPoint.z;
                      if (smos_Fabs(dNx) == 1.0) { rH.a[x] = dPx * rH.a[3]; }
                    }
                }
        }
    };
  // each cycle patch's round rows (columns from face m to face m + 1), the rails' ends and
  // their directions toward the corner
  std::vector< std::vector<sm_SbH> > sH0(n), sH1(n);
  std::vector<SmPoint3d> sEnd0(n), sEnd1(n), sMid(n);
  std::vector<SmVector3d> sDir0(n), sDir1(n);
  std::vector<double> sSpeed0(n), sSpeed1(n);
  for (size_t m=0; m<n; m++)
    {
      const size_t r = crRound[m], q = crQ[r];
      const double dSign = crRounds[r].m_bCornerAbove ? 1.0 : -1.0;
      sH0[m].resize(q + 1);
      sH1[m].resize(q + 1);
      for (size_t c=0; c<=q; c++)
        {
          const size_t cc = crRev[m] ? q - c : c;
          sH0[m][c] = crJets[r][0][cc];
          sH1[m][c] = sm_SbAxpy(sm_SbZero(), dSign, crJets[r][1][cc]);
        }
      // the round's end columns exactly on their faces (within dTol): its rail end on the face
      // (row 0) and its derivative along it (row 1), as the fill's own face columns are
      for (int e=0; e<2; e++)
        {
          const SmBlendSetback::Face & crF = crFaces[crFace[(m + static_cast<size_t>(e)) % n]];
          const SmVector3d & crN = crF.m_sNormal;
          sm_SbH * pH[2] = { &sH0[m][e == 0 ? 0 : q], &sH1[m][e == 0 ? 0 : q] };
          for (int k=0; k<2; k++)
            {
              sm_SbH & rH = *pH[k];
              const double dOff = crN.x * rH.a[0] + crN.y * rH.a[1] + crN.z * rH.a[2] - crN.Dot(SmVector3d(crF.m_sPoint)) * rH.a[3];
              if (smos_Fabs(dOff) > dTol * smos_Max(1.0, smos_Fabs(sH0[m][e == 0 ? 0 : q].a[3])))
                { SM_SB_FAIL(rFill.m_pReason, SM_ERR_INVALID_INPUT, _T("SmBlendSetback: a round does not meet its faces at its end")); }
              rH.a[0] -= dOff * crN.x;
              rH.a[1] -= dOff * crN.y;
              rH.a[2] -= dOff * crN.z;
              for (int x=0; x<3; x++)
                {
                  const double dNx = (x == 0) ? crN.x : (x == 1) ? crN.y : crN.z;
                  const double dPx = (x == 0) ? crF.m_sPoint.x : (x == 1) ? crF.m_sPoint.y : crF.m_sPoint.z;
                  if (smos_Fabs(dNx) == 1.0) { rH.a[x] = dPx * rH.a[3]; }
                }
            }
        }
      for (int e=0; e<2; e++)
        {
          const sm_SbH & crA = sH0[m][e == 0 ? 0 : q];
          const sm_SbH & crB = sH1[m][e == 0 ? 0 : q];
          const SmPoint3d sP(crA.a[0] / crA.a[3], crA.a[1] / crA.a[3], crA.a[2] / crA.a[3]);
          SmVector3d sT((crB.a[0] - sP.x * crB.a[3]) / crA.a[3], (crB.a[1] - sP.y * crB.a[3]) / crA.a[3], (crB.a[2] - sP.z * crB.a[3]) / crA.a[3]);
          (e == 0 ? sSpeed0[m] : sSpeed1[m]) = sT.Length();
          SER(sT.Unitize());
          (e == 0 ? sEnd0[m] : sEnd1[m]) = sP;
          (e == 0 ? sDir0[m] : sDir1[m]) = sT;
        }
      // the section's middle
      sm_SbB sB(q + 1);
      sm_SbH sMh = sm_SbZero();
      for (size_t c=0; c<=q; c++)
        { sMh = sm_SbAxpy(sMh, sm_Binomial(static_cast<ULONG>(q), static_cast<ULONG>(c)) / std::pow(2.0, static_cast<double>(q)), sH0[m][c]); }
      sMid[m] = SmPoint3d(sMh.a[0] / sMh.a[3], sMh.a[1] / sMh.a[3], sMh.a[2] / sMh.a[3]);
    }
  // convex (the rounds take material away) or concave (they add it), and the axis e into the
  // convex side
  SmVector3d sSum(0.0, 0.0, 0.0);
  int iConvex = 0;
  for (size_t m=0; m<n; m++)
    {
      const SmBlendSetback::Face & crF = crFaces[crFace[m]];
      sSum = sSum + crF.m_sNormal;
      const int iSide = ((sMid[m] - crF.m_sPoint).Dot(crF.m_sNormal) < 0.0) ? 1 : -1;
      if (iConvex != 0 && iSide != iConvex)
        { SM_SB_FAIL(rFill.m_pReason, SM_ERR_INVALID_INPUT, _T("SmBlendSetback: a convex G1 fill needs rounds all convex or all concave")); }
      iConvex = iSide;
    }
  SmVector3d sE = sSum * (iConvex > 0 ? -1.0 : 1.0);
  SER(sE.Unitize());
  SmVector3d sA1(1.0, 0.0, 0.0);
  if (smos_Fabs(sE.y) < smos_Fabs(sE.x) && smos_Fabs(sE.y) <= smos_Fabs(sE.z)) { sA1 = SmVector3d(0.0, 1.0, 0.0); }
  else if (smos_Fabs(sE.z) < smos_Fabs(sE.x) && smos_Fabs(sE.z) < smos_Fabs(sE.y)) { sA1 = SmVector3d(0.0, 0.0, 1.0); }
  sA1 = sA1 - sE * sA1.Dot(sE);
  SER(sA1.Unitize());
  const SmVector3d sA2 = sE * sA1;
  auto inPlane = [&](const SmVector3d & crV) { return SmVector2d(crV.Dot(sA1), crV.Dot(sA2)); };

  // the base points, the domain (their projections), the cross fields and the speeds
  std::vector<SmPoint2d> sXi(n);
  std::vector<SmVector2d> sSigma(n);
  std::vector<double> sToX0(n), sToX1(n);
  for (size_t m=0; m<n; m++)
    {
      const size_t mp = (m + n - 1) % n;
      const SmPoint3d & crPa = sEnd0[m];   // patch m's rail on face m
      const SmPoint3d & crPb = sEnd1[mp];  // patch m - 1's rail on face m
      const SmVector3d & crTa = sDir0[m], & crTb = sDir1[mp];
      // the rails' crossing: Pa + s Ta = Pb + t Tb, least squares
      const double a11 = crTa.Dot(crTa), a12 = -crTa.Dot(crTb), a22 = crTb.Dot(crTb);
      const SmVector3d sD = crPb - crPa;
      const double b1 = crTa.Dot(sD), b2 = -crTb.Dot(sD);
      const double dDet = a11 * a22 - a12 * a12;
      if (!(smos_Fabs(dDet) > 1.0e-12))
        { SM_SB_FAIL(rFill.m_pReason, SM_ERR_INVALID_INPUT, _T("SmBlendSetback: two rails on a face run parallel (no corner region between them)")); }
      const double s = (b1 * a22 - a12 * b2) / dDet, t = (a11 * b2 - a12 * b1) / dDet;
      const SmVector3d sX = (SmVector3d(crPa) + crTa * s + SmVector3d(crPb) + crTb * t) * 0.5;
      const SmVector3d sMidC = (SmVector3d(crPa) + SmVector3d(crPb)) * 0.5;
      if (!(s > 0.0) || !(t > 0.0))
        { SM_SB_FAIL(rFill.m_pReason, SM_ERR_INVALID_INPUT, _T("SmBlendSetback: two rounds' rails on a face cross behind their ends (no corner region between them)")); }
      // the spring curve on face m: the parabola (Pa, X, Pb), as two cubic halves meeting at
      // its middle B = (Pa + 2 X + Pb) / 4, its tangent there (Pa - Pb) / 2 (sigma_m), each
      // half's second control point on its rail a third of the way to X (lambda)
      const SmVector3d sBase = sMidC + (sX - sMidC) * 0.5;
      const SmVector2d sXi2 = inPlane(sBase);
      sXi[m] = SmPoint2d(sXi2.x, sXi2.y);
      sSigma[m] = inPlane(crPa - crPb) * 0.5;
      sToX0[m] = s / sSpeed0[m];
      sToX1[mp] = t / sSpeed1[mp];
    }
  SmPoint2d sCentre(0.0, 0.0);
  double dRho = 0.0;
  for (size_t m=0; m<n; m++) { sCentre = SmPoint2d(sCentre.x + sXi[m].x / static_cast<double>(n), sCentre.y + sXi[m].y / static_cast<double>(n)); }
  for (size_t m=0; m<n; m++) { dRho = smos_Max(dRho, SmVector2d(sXi[m] - sCentre).Length()); }
  if (!(dRho > 0.0))
    { SER(SM_ERR); }

  // the domain scaled: x = (xi - centre) / rho
  std::vector<SmPoint2d> sB(n);
  for (size_t m=0; m<n; m++) { sB[m] = SmPoint2d((sXi[m].x - sCentre.x) / dRho, (sXi[m].y - sCentre.y) / dRho); }

  // The degrees: 9 first, then up from the least. A map of degree d is one of degree d + 1, and
  // its pieces' coefficient matrices at d + 1 are convex combinations of those at d (degree
  // raising), so the sector conditions met at d are met at d + 1: proven unmet at 9 (the QP's
  // certificate), they are unmet at every degree, and the search stops. Else the least degree
  // whose fill is certified, or 9's.
  std::vector<sm_SbNet> sTop;
  size_t lTopDegree = 0;
  const bool bProbe = (pProbe != NULL && !pProbe->m_bFill);
  std::vector<size_t> sDegrees;
  if (pProbe == NULL || bProbe) { sDegrees.push_back(9); }
  for (size_t d=Nc; d<=9 && !bProbe; d++)
    if (d < 9 || (pProbe != NULL && pProbe->m_bFill)) { sDegrees.push_back(d); }
  for (size_t lD=0; lD<sDegrees.size(); lD++)
    {
      const size_t d = sDegrees[lD];
      const std::vector< std::pair<int, int> > sMono = sm_SbMonomials(static_cast<int>(d));
      const size_t M = sMono.size();
      // the unknowns: h's coefficients, then each setback patch's free heights
      std::vector<size_t> sQv(n);
      const size_t lUnknowns = M;
      for (size_t m=0; m<n; m++) { sQv[m] = smos_Max(crQ[crRound[m]] + 1, d); }
      std::vector<sm_SbAffinePatch> sPatches(2 * n);
      for (size_t m=0; m<n; m++)
        {
          const size_t qv = sQv[m];
          // the cubic across: rows 0, 1 the round's, at the speed lambda(v) off it (linear across,
          // from each rail's distance to its rails' crossing: the spring curves' first legs), 2, 3 F's
          std::vector< std::vector<sm_SbH> > sRows(Nc + 1, std::vector<sm_SbH>(qv + 1, sm_SbZero()));
          const sm_SbB sLam = { sToX0[m] / static_cast<double>(Nc), sToX1[m] / static_cast<double>(Nc) };
          for (int c=0; c<4; c++)
            {
              sm_SbB s0(sH0[m].size()), s1(sH0[m].size());
              for (size_t j=0; j<s0.size(); j++) { s0[j] = sH0[m][j].a[c]; s1[j] = sH1[m][j].a[c]; }
              s1 = sm_SbMul(sLam, s1);
              s0 = sm_SbElevate(s0, qv);
              s1 = sm_SbElevate(s1, qv);
              for (size_t j=0; j<=qv; j++) { sRows[0][j].a[c] = s0[j]; sRows[1][j].a[c] = s0[j] + s1[j]; }
            }
          // F's projection along edge m and its derivative along the field d (linear in v)
          const SmPoint2d & crX0 = sXi[m];
          const SmPoint2d & crX1 = sXi[(m + 1) % n];
          const SmVector2d sD0 = sSigma[m], sD1 = sSigma[(m + 1) % n] * -1.0;
          for (int c=0; c<3; c++)
            {
              const double dA1 = (c == 0) ? sA1.x : (c == 1) ? sA1.y : sA1.z;
              const double dA2 = (c == 0) ? sA2.x : (c == 1) ? sA2.y : sA2.z;
              const sm_SbB sP = { crX0.x * dA1 + crX0.y * dA2, crX1.x * dA1 + crX1.y * dA2 };
              const sm_SbB sDP = { sD0.x * dA1 + sD0.y * dA2, sD1.x * dA1 + sD1.y * dA2 };
              const sm_SbB s3 = sm_SbElevate(sP, qv);
              const sm_SbB sDe = sm_SbElevate(sDP, qv);
              for (size_t j=0; j<=qv; j++)
                {
                  sRows[3][j].a[c] = s3[j];
                  sRows[2][j].a[c] = s3[j] + sDe[j] / static_cast<double>(Nc);
                }
            }
          for (size_t j=0; j<=qv; j++) { sRows[2][j].a[3] = 1.0; sRows[3][j].a[3] = 1.0; }
          // h's monomials along the edge (scaled coordinates) and along the field (d / rho)
          const sm_SbB sX = { sB[m].x, sB[(m + 1) % n].x }, sY = { sB[m].y, sB[(m + 1) % n].y };
          const sm_SbB sDX = { sD0.x / dRho, sD1.x / dRho }, sDY = { sD0.y / dRho, sD1.y / dRho };
          std::vector< std::vector<double> > sZ2(M), sZ3(M);
          for (size_t j=0; j<M; j++)
            {
              const int a = sMono[j].first, b = sMono[j].second;
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
              const sm_SbB v0 = sm_SbElevate(sm_SbAlong(a, b, sX, sY, d), qv), v1 = sm_SbElevate(s1, qv);
              sZ3[j] = v0;
              sZ2[j].resize(qv + 1);
              for (size_t i=0; i<=qv; i++) { sZ2[j][i] = v0[i] + v1[i] / static_cast<double>(Nc); }
            }
          // the net
          sm_SbAffinePatch & rP = sPatches[m];
          rP.m_sBase.p = Nc;
          rP.m_sBase.q = qv;
          rP.m_sBase.P.assign((Nc + 1) * (qv + 1), sm_SbZero());
          for (size_t j=0; j<=qv; j++)
            for (int c=0; c<4; c++)
              {
                for (size_t i=0; i<=Nc; i++) { rP.m_sBase.P[i * (qv + 1) + j].a[c] = sRows[i][j].a[c]; }
              }
          for (size_t jm=0; jm<M; jm++)
            {
              sm_SbT sZ = sm_SbTMake(Nc, qv);
              for (size_t j=0; j<=qv; j++)
                {
                  sZ.c[2 * (qv + 1) + j] = sZ2[jm][j];
                  sZ.c[3 * (qv + 1) + j] = sZ3[jm][j];
                }
              rP.m_vIndex.push_back(jm);
              rP.m_vZ.push_back(sZ);
            }
        }
      for (size_t q=0; q<n; q++)
        {
          // interior part q: corners (0, 0) the centre, (0, 1) the middle of edge q - 1, (1, 0)
          // the middle of edge q, (1, 1) B_q
          const SmPoint2d sC(0.0, 0.0);
          const SmPoint2d sMq((sB[q].x + sB[(q + 1) % n].x) * 0.5, (sB[q].y + sB[(q + 1) % n].y) * 0.5);
          const SmPoint2d sMp((sB[q].x + sB[(q + n - 1) % n].x) * 0.5, (sB[q].y + sB[(q + n - 1) % n].y) * 0.5);
          sm_SbT sX, sY;
          sX.p = sX.q = sY.p = sY.q = 1;
          sX.c = { sC.x, sMp.x, sMq.x, sB[q].x };
          sY.c = { sC.y, sMp.y, sMq.y, sB[q].y };
          sm_SbAffinePatch & rP = sPatches[n + q];
          rP.m_sBase.p = rP.m_sBase.q = d;
          rP.m_sBase.P.assign((d + 1) * (d + 1), sm_SbZero());
          const sm_SbT sXd = sm_SbTElevate(sX, d, d), sYd = sm_SbTElevate(sY, d, d);
          for (size_t i=0; i<(d + 1) * (d + 1); i++)
            {
              const double xi1 = sCentre.x + dRho * sXd.c[i], xi2 = sCentre.y + dRho * sYd.c[i];
              rP.m_sBase.P[i].a[0] = xi1 * sA1.x + xi2 * sA2.x;
              rP.m_sBase.P[i].a[1] = xi1 * sA1.y + xi2 * sA2.y;
              rP.m_sBase.P[i].a[2] = xi1 * sA1.z + xi2 * sA2.z;
              rP.m_sBase.P[i].a[3] = 1.0;
            }
          for (size_t j=0; j<M; j++)
            {
              rP.m_vIndex.push_back(j);
              rP.m_vZ.push_back(sm_SbOver(sMono[j].first, sMono[j].second, sX, sY, d));
            }
        }
      // the equalities: each setback patch's rows 2, 3, its first two columns on face m's
      // plane and its last two on face m + 1's
      std::vector<double> sAeq, sBeq;
      for (size_t m=0; m<n; m++)
        {
          const sm_SbAffinePatch & crP = sPatches[m];
          const size_t qv = crP.m_sBase.q;
          for (size_t i=2; i<=Nc; i++)
            for (size_t c=0; c<2; c++)
              for (int e=0; e<2; e++)
                {
                  const size_t col = (e == 0) ? c : qv - c, idx = i * (qv + 1) + col;
                  const SmBlendSetback::Face & crF = crFaces[crFace[(m + static_cast<size_t>(e)) % n]];
                  const sm_SbH & crH = crP.m_sBase.P[idx];
                  const SmVector3d & crN = crF.m_sNormal;
                  std::vector<double> sRow(lUnknowns, 0.0);
                  for (size_t k=0; k<crP.m_vIndex.size(); k++) { sRow[crP.m_vIndex[k]] += crP.m_vZ[k].c[idx] * crN.Dot(sE); }
                  sAeq.insert(sAeq.end(), sRow.begin(), sRow.end());
                  sBeq.push_back(-(crN.x * crH.a[0] + crN.y * crH.a[1] + crN.z * crH.a[2] - crN.Dot(SmVector3d(crF.m_sPoint)) * crH.a[3]));
                }
        }
      std::vector<double> sYp;
      std::vector< std::vector<double> > sZn;
      double dSize = 0.0;
      for (size_t i=0; i<sBeq.size(); i++) { dSize = smos_Max(dSize, smos_Fabs(sBeq[i])); }
      if (!sm_SbAffineSolutions(sAeq, sBeq, sBeq.size(), lUnknowns, 1.0e-9 * (1.0 + dSize), sYp, sZn) || sZn.empty())
        { continue; }
      const size_t lFree = sZn.size();


      // the energy, affine in y: its rows at y = 0 (sE0) and their change per unknown (sAE,
      // lE x lUnknowns)
      auto energyRows = [&](const sm_SbNet & crN, std::vector<double> & rRows)
        {
          rRows.clear();
          for (int g1=0; g1<8; g1++)
            for (int g2=0; g2<8; g2++)
              {
                double R[3][3];
                sm_SbSecond(crN, SM_SB_GX[g1], SM_SB_GX[g2], R);
                const double w = smos_Sqrt(SM_SB_GW[g1] * SM_SB_GW[g2]);
                for (int c=0; c<3; c++)
                  {
                    rRows.push_back(w * R[0][c]);
                    rRows.push_back(w * smos_Sqrt(2.0) * R[1][c]);
                    rRows.push_back(w * R[2][c]);
                  }
              }
        };
      std::vector<double> sE0;
      std::vector<size_t> sEStart;
      std::vector< std::vector<double> > sPatchE0(2 * n);
      for (size_t p=0; p<2*n; p++)
        {
          energyRows(sPatches[p].m_sBase, sPatchE0[p]);
          sEStart.push_back(sE0.size());
          sE0.insert(sE0.end(), sPatchE0[p].begin(), sPatchE0[p].end());
        }
      const size_t lE = sE0.size();
      std::vector<double> sAE(lE * lUnknowns, 0.0);
      {
        std::vector<double> sUnit(lUnknowns, 0.0), sRows;
        for (size_t p=0; p<2*n; p++)
          for (size_t t=0; t<sPatches[p].m_vIndex.size(); t++)
            {
              const size_t u = sPatches[p].m_vIndex[t];
              sUnit[u] = 1.0;
              energyRows(sm_SbAffineNet(sPatches[p], sUnit, sE), sRows);
              sUnit[u] = 0.0;
              for (size_t i=0; i<sRows.size(); i++) { sAE[(sEStart[p] + i) * lUnknowns + u] += sRows[i] - sPatchE0[p][i]; }
            }
      }
      // over y = yp + Zn t: B = AE Zn, c = E0 + AE yp; the energy |B t + c|^2
      std::vector<double> sBm(lE * lFree, 0.0), sCv(lE, 0.0);
      for (size_t i=0; i<lE; i++)
        {
          double s = sE0[i];
          for (size_t u=0; u<lUnknowns; u++) { s += sAE[i * lUnknowns + u] * sYp[u]; }
          sCv[i] = s;
          for (size_t f=0; f<lFree; f++)
            {
              double t = 0.0;
              for (size_t u=0; u<lUnknowns; u++) { t += sAE[i * lUnknowns + u] * sZn[f][u]; }
              sBm[i * lFree + f] = t;
            }
        }
      // the nets at t, their face points snapped exactly onto their planes
      auto netsAt = [&](const std::vector<double> & crT, std::vector<sm_SbNet> & rNets) -> SmStatus
        {
          std::vector<double> sY(sYp);
          for (size_t f=0; f<lFree; f++)
            for (size_t u=0; u<lUnknowns; u++) { sY[u] += sZn[f][u] * crT[f]; }
          rNets.assign(2 * n, sm_SbNet());
          for (size_t p=0; p<2*n; p++) { rNets[p] = sm_SbAffineNet(sPatches[p], sY, sE); }
          for (size_t m=0; m<n; m++)
            {
              sm_SbNet & rN = rNets[m];
              const size_t qv = rN.q;
              for (size_t i=2; i<=Nc; i++)
                for (size_t c=0; c<2; c++)
                  for (int e=0; e<2; e++)
                    {
                      sm_SbH & rH = rN.P[i * (qv + 1) + ((e == 0) ? c : qv - c)];
                      const SmBlendSetback::Face & crF = crFaces[crFace[(m + static_cast<size_t>(e)) % n]];
                      const SmVector3d & crN = crF.m_sNormal;
                      const double dOff = crN.x * rH.a[0] + crN.y * rH.a[1] + crN.z * rH.a[2] - crN.Dot(SmVector3d(crF.m_sPoint)) * rH.a[3];
                      if (smos_Fabs(dOff) > dTol * rH.a[3])
                        { return SM_ERR; }
                      rH.a[0] -= dOff * crN.x;
                      rH.a[1] -= dOff * crN.y;
                      rH.a[2] -= dOff * crN.z;
                      for (int x=0; x<3; x++)
                        {
                          const double dNx = (x == 0) ? crN.x : (x == 1) ? crN.y : crN.z;
                          const double dPx = (x == 0) ? crF.m_sPoint.x : (x == 1) ? crF.m_sPoint.y : crF.m_sPoint.z;
                          if (smos_Fabs(dNx) == 1.0) { rH.a[x] = dPx * rH.a[3]; }
                        }
                    }
            }
          return SM_SUCCESS;
        };

      // each patch's curvature determinants, affine in y, and the side its normal faces (o of
      // one sign, certified, or the fill folds)
      std::vector< std::array<sm_SbT, 3> > sQ0(2 * n);
      std::vector< std::array<std::vector<sm_SbT>, 3> > sQj(2 * n);
      std::vector<double> sS(2 * n, 0.0);
      for (size_t p=0; p<2*n; p++)
        {
          sm_SbT sQ0p[3], sO;
          std::vector<sm_SbT> sQjp[3];
          sm_SbAffineFundamentals(sPatches[p], sE, sQ0p, sQjp, sO);
          sm_SbT sNegO = sO;
          for (size_t i=0; i<sNegO.c.size(); i++) { sNegO.c[i] = -sNegO.c[i]; }
          const double dO = 1.0e-12 * sm_SbTScale(sO);
          sS[p] = (sm_SbTNonneg(sO, dO, SM_SB_GRAPH_DEPTH) == 1) ? 1.0 : (sm_SbTNonneg(sNegO, dO, SM_SB_GRAPH_DEPTH) == 1) ? -1.0 : 0.0;
          if (sS[p] == 0.0)
            { SM_SB_FAIL(rFill.m_pReason, SM_ERR, _T("SmBlendSetback: the fill folds (a patch's normal turns across the corner's axis)")); }
          for (int k=0; k<3; k++) { sQ0[p][k] = sQ0p[k]; sQj[p][k] = sQjp[k]; }
        }
      // the energy as a quadratic programme: G = 2 B'B (a least nudge where it leaves a direction
      // free: the energy only chooses), g = 2 B'c
      std::vector<double> sG(lFree * lFree, 0.0), sg(lFree, 0.0);
      double dTrace = 0.0;
      for (size_t a=0; a<lFree; a++)
        {
          for (size_t b=a; b<lFree; b++)
            {
              double s = 0.0;
              for (size_t i=0; i<lE; i++) { s += sBm[i * lFree + a] * sBm[i * lFree + b]; }
              sG[a * lFree + b] = sG[b * lFree + a] = 2.0 * s;
            }
          double s = 0.0;
          for (size_t i=0; i<lE; i++) { s += sBm[i * lFree + a] * sCv[i]; }
          sg[a] = 2.0 * s;
          dTrace += sG[a * lFree + a];
        }
      for (size_t a=0; a<lFree; a++) { sG[a * lFree + a] += 1.0e-12 * dTrace / static_cast<double>(lFree); }
      // each interior part's extra direction: at its base point (its (1, 1)) F's Hessian is forced
      // singular along sigma (the faces' conditions there); sigma in its (s, t)
      std::vector<double> sKernel(2 * n, -1.0);
      for (size_t q=0; q<n; q++)
        {
          const SmVector2d sDs((sB[q].x - sB[(q + n - 1) % n].x) * 0.5, (sB[q].y - sB[(q + n - 1) % n].y) * 0.5);
          const SmVector2d sDt((sB[q].x - sB[(q + 1) % n].x) * 0.5, (sB[q].y - sB[(q + 1) % n].y) * 0.5);
          const double dDet = sDs.x * sDt.y - sDs.y * sDt.x;
          if (!(smos_Fabs(dDet) > 0.0))
            { SER(SM_ERR); }
          const double a = (sSigma[q].x * sDt.y - sSigma[q].y * sDt.x) / dDet, b = (sDs.x * sSigma[q].y - sDs.y * sSigma[q].x) / dDet;
          double th = std::atan2(b, a);
          if (th < 0.0) { th += SM_PI; }
          if (th >= SM_PI) { th -= SM_PI; }
          sKernel[n + q] = th;
        }

      {
        const size_t lLevels = SM_SB_GRAPH_LEVELS;
        {
          // the sector conditions on every piece's coefficient matrices, in t: C t >= b
          std::vector<double> sC, sb;
          for (size_t p=0; p<2*n; p++)
            {
              const sm_SbAffinePatch & crP = sPatches[p];
              std::vector<double> sTheta;
              for (size_t k=0; k<SM_SB_GRAPH_DIRECTIONS; k++) { sTheta.push_back(SM_PI * static_cast<double>(k) / static_cast<double>(SM_SB_GRAPH_DIRECTIONS)); }
              if (sKernel[p] >= 0.0) { sTheta.push_back(sKernel[p]); }
              std::sort(sTheta.begin(), sTheta.end());
              sTheta.push_back(sTheta[0] + SM_PI);
              size_t P = 0, Q = 0;
              for (int k=0; k<3; k++) { P = smos_Max(P, sQ0[p][k].p); Q = smos_Max(Q, sQ0[p][k].q); }
              // l, m, n over t: at t = 0 (y = yp) and per free direction, at common degrees, in pieces
              std::vector< std::vector<sm_SbT> > sPc0(3), sPcf(3 * lFree);
              double dScale = 0.0;
              for (int k=0; k<3; k++)
                {
                  sm_SbT s0 = sm_SbTElevate(sQ0[p][k], P, Q);
                  std::vector<sm_SbT> sF(lFree, sm_SbTMake(P, Q));
                  for (size_t t=0; t<crP.m_vIndex.size(); t++)
                    {
                      const sm_SbT sJ = sm_SbTElevate(sQj[p][k][t], P, Q);
                      dScale = smos_Max(dScale, sm_SbTScale(sJ));
                      const size_t u = crP.m_vIndex[t];
                      for (size_t i=0; i<sJ.c.size(); i++) { s0.c[i] += sYp[u] * sJ.c[i]; }
                      for (size_t f=0; f<lFree; f++)
                        {
                          const double z = sZn[f][u];
                          if (z == 0.0) { continue; }
                          for (size_t i=0; i<sJ.c.size(); i++) { sF[f].c[i] += z * sJ.c[i]; }
                        }
                    }
                  dScale = smos_Max(dScale, sm_SbTScale(sQ0[p][k]));
                  sm_SbTPieces(s0, lLevels, sPc0[k]);
                  for (size_t f=0; f<lFree; f++) { sm_SbTPieces(sF[f], lLevels, sPcf[k * lFree + f]); }
                }
              std::vector<double> sRowT(lFree);
              for (size_t pc=0; pc<sPc0[0].size(); pc++)
                for (size_t i=0; i<sPc0[0][pc].c.size(); i++)
                  for (size_t a=0; a+1<sTheta.size(); a++)
                    for (int iPair=0; iPair<2; iPair++)
                      {
                        const double c0 = std::cos(sTheta[a]), s0 = std::sin(sTheta[a]);
                        const double c1 = (iPair == 0) ? c0 : std::cos(sTheta[a + 1]), s1 = (iPair == 0) ? s0 : std::sin(sTheta[a + 1]);
                        const double w[3] = { c0 * c1, c0 * s1 + s0 * c1, s0 * s1 };
                        double dConst = 0.0, dMax = 0.0;
                        for (int k=0; k<3; k++) { dConst += sS[p] * w[k] * sPc0[k][pc].c[i]; }
                        for (size_t f=0; f<lFree; f++)
                          {
                            double s = 0.0;
                            for (int k=0; k<3; k++) { s += sS[p] * w[k] * sPcf[k * lFree + f][pc].c[i]; }
                            sRowT[f] = s;
                            dMax = smos_Max(dMax, smos_Fabs(s));
                          }
                        if (!(dMax > 1.0e-11 * dScale))
                          { continue; } // no free height in it: the certificate decides it
                        for (size_t f=0; f<lFree; f++) { sC.push_back(sRowT[f] / dMax); }
                        sb.push_back((-dConst + 1.0e-12 * dScale) / dMax); // held a little inside
                      }
            }
          std::vector<double> sT;
          rFill.m_lSolves++;
          if (bProbe && pProbe->m_bMargin)
            {
              // the least tau with C t + tau >= b (tau >= -1): min tau + eps (|t|^2 + tau^2)
              const size_t lN = lFree + 1, lM = sb.size();
              std::vector<double> sGm(lN * lN, 0.0), sgm(lN, 0.0), sCm(lM * lN + lN, 0.0), sbm(sb);
              double dG = 0.0;
              for (size_t a=0; a<lFree; a++) { dG = smos_Max(dG, sG[a * lFree + a]); }
              for (size_t a=0; a<lN; a++) { sGm[a * lN + a] = 1.0e-9; }
              sgm[lFree] = 1.0;
              for (size_t i=0; i<lM; i++)
                {
                  for (size_t f=0; f<lFree; f++) { sCm[i * lN + f] = sC[i * lFree + f]; }
                  sCm[i * lN + lFree] = 1.0;
                }
              for (size_t f=0; f<lN; f++) { sCm[lM * lN + f] = (f == lFree) ? 1.0 : 0.0; }
              sbm.push_back(-1.0);
              std::vector<double> sTm;
              const SmStatus eM = sm_SbDualQPRows(sGm, sgm, sCm, sbm, lN, lM + 1, 1.0e-13, sTm);
              if (eM != SM_SUCCESS)
                { SER(eM); }
              pProbe->m_dMargin = sTm[lFree];
              (void)dG;
              return SM_SUCCESS;
            }
          bool bProven = false;
          const SmStatus eQP = sm_SbDualQPRows(sG, sg, sC, sb, lFree, sb.size(), 1.0e-13, sT, &bProven);
          if (bProbe)
            {
              pProbe->m_bFeasible = (eQP == SM_SUCCESS);
              return SM_SUCCESS;
            }
          if (eQP != SM_SUCCESS && bProven && d == 9)
            { break; }
          if (eQP != SM_SUCCESS)
            { continue; }
          std::vector<sm_SbNet> sNets;
          SER(netsAt(sT, sNets));
          size_t lBad = 0;
          if (!sm_SbCertifiedConvex(sNets, sE, SM_SB_GRAPH_DEPTH, lBad))
            { continue; }
          if (d == 9 && pProbe == NULL)
            {
              sTop = sNets;
              lTopDegree = d;
              continue;
            }
          toWorld(sNets);
          SER(sm_SbFillFrom(crContext, sNets, crRound, crFace, crRev, d, rFill));
          return SM_SUCCESS;
        }
      }
    }
  if (!sTop.empty())
    {
      toWorld(sTop);
      SER(sm_SbFillFrom(crContext, sTop, crRound, crFace, crRev, lTopDegree, rFill));
      return SM_SUCCESS;
    }
  rFill.m_bUncertified = true;
  SM_SB_FAIL(rFill.m_pReason, SM_ERR, _T("SmBlendSetback: no G1 fill up to degree 9 is certified convex at these setbacks"));
}

// det(H, Hu, Hv, (V, 0)) of a homogeneous net and a polynomial vector field V (3 components).
sm_SbT sm_SbTDetWith(const sm_SbNet & crN, const sm_SbT (&crV)[3])
{
  sm_SbT h[4], hu[4], hv[4];
  for (int c=0; c<4; c++)
    {
      h[c] = sm_SbTMake(crN.p, crN.q);
      for (size_t i=0; i<crN.P.size(); i++) { h[c].c[i] = crN.P[i].a[c]; }
      hu[c] = sm_SbTDu(h[c]);
      hv[c] = sm_SbTDv(h[c]);
    }
  const sm_SbT * pRows[3] = { h, hu, hv };
  sm_SbT s = sm_SbTMake(0, 0);
  for (int c=0; c<3; c++)
    {
      const sm_SbT * pM[3][3];
      for (int r=0; r<3; r++)
        {
          int cc = 0;
          for (int col=0; col<4; col++)
            if (col != c) { pM[r][cc++] = &pRows[r][col]; }
        }
      s = sm_SbTAxpy(s, ((3 + c) % 2 == 0) ? 1.0 : -1.0, sm_SbTMulFast(crV[c], sm_SbTDet3(pM)));
    }
  return s;
}

// A setback patch's normal field (ruling (a)'s certificate): V (3 components) and g > 0 with
// V / g = (1 - v) N_A + v N_B, N_F the unit normal field of face F (crA on the column v = 0, crB on
// v = 1) from F's implicit form at the patch's own points, turned to crE's side: N_F = G_F / q_F,
// G_F and q_F polynomial in the net, |G_F| = q_F wherever the point is on F. Cleared:
// V = (1 - v) qB GA + v qA GB, g = qA qB. Along a column where the patch is G1 to its face,
// V / g is the patch's unit normal there.
void sm_SbNormalField(const sm_SbNet & crN, const SmBlendSetback::Face & crA, const SmBlendSetback::Face & crB, const SmVector3d & crE,
                      sm_SbT (&rV)[3], sm_SbT & rG)
{
  sm_SbT sGs[2][3], sq[2];
  for (int e=0; e<2; e++)
    {
      const SmBlendSetback::Face & crF = (e == 0) ? crA : crB;
      const size_t j = (e == 0) ? 0 : crN.q;
      if (crF.m_eKind == SmBlendSetback::SM_SB_PLANE)
        {
          for (int c=0; c<3; c++)
            {
              sGs[e][c] = sm_SbTMake(0, 0);
              sGs[e][c].c[0] = (c == 0) ? crF.m_sNormal.x : (c == 1) ? crF.m_sNormal.y : crF.m_sNormal.z;
            }
          sq[e] = sm_SbTMake(0, 0);
          sq[e].c[0] = 1.0;
        }
      else
        {
          const double dP[3] = { crF.m_sPoint.x, crF.m_sPoint.y, crF.m_sPoint.z };
          const double dA[3] = { crF.m_sAxis.x, crF.m_sAxis.y, crF.m_sAxis.z };
          sm_SbT sU[3], sW = sm_SbTMake(crN.p, crN.q), sH = sm_SbTMake(crN.p, crN.q);
          for (size_t i=0; i<crN.P.size(); i++) { sW.c[i] = crN.P[i].a[3]; }
          for (int c=0; c<3; c++)
            {
              sU[c] = sm_SbTMake(crN.p, crN.q);
              for (size_t i=0; i<crN.P.size(); i++) { sU[c].c[i] = crN.P[i].a[c] - sW.c[i] * dP[c]; }
              for (size_t i=0; i<crN.P.size(); i++) { sH.c[i] += sU[c].c[i] * dA[c]; }
            }
          if (crF.m_eKind == SmBlendSetback::SM_SB_SPHERE)
            {
              for (int c=0; c<3; c++) { sGs[e][c] = sU[c]; }
              sq[e] = sm_SbTAxpy(sm_SbTMake(crN.p, crN.q), crF.m_dRadius, sW);
            }
          else
            {
              // cylinder (slope 0) or cone: rho w = R w + slope h; G = U - h a - slope rho w a,
              // |G| = sqrt(1 + slope^2) rho w on the face: the constant kept in q
              const sm_SbT sRho = sm_SbTAxpy(sm_SbTAxpy(sm_SbTMake(crN.p, crN.q), crF.m_dRadius, sW), crF.m_dSlope, sH);
              for (int c=0; c<3; c++)
                { sGs[e][c] = sm_SbTAxpy(sm_SbTAxpy(sU[c], -dA[c], sH), -crF.m_dSlope * dA[c], sRho); }
              sq[e] = sm_SbTAxpy(sm_SbTMake(crN.p, crN.q), std::sqrt(1.0 + crF.m_dSlope * crF.m_dSlope), sRho);
            }
        }
      // its side: G . e > 0 at the column's middle
      double dDot = 0.0;
      for (int c=0; c<3; c++)
        {
          const sm_SbT & crC = sGs[e][c];
          const size_t jj = (crC.q == 0) ? 0 : j;
          double dMid = 0.0;
          for (size_t i=0; i<=crC.p; i++)
            { dMid += crC.c[i * (crC.q + 1) + jj] * sm_Binomial(static_cast<ULONG>(crC.p), static_cast<ULONG>(i)); }
          dDot += dMid * ((c == 0) ? crE.x : (c == 1) ? crE.y : crE.z);
        }
      if (dDot < 0.0)
        for (int c=0; c<3; c++)
          for (size_t i=0; i<sGs[e][c].c.size(); i++) { sGs[e][c].c[i] = -sGs[e][c].c[i]; }
    }
  sm_SbT sLo = sm_SbTMake(0, 1), sHi = sm_SbTMake(0, 1);
  sLo.c = { 1.0, 0.0 };
  sHi.c = { 0.0, 1.0 };
  for (int c=0; c<3; c++)
    { rV[c] = sm_SbTAxpy(sm_SbTMulFast(sm_SbTMulFast(sLo, sq[1]), sGs[0][c]), 1.0, sm_SbTMulFast(sm_SbTMulFast(sHi, sq[0]), sGs[1][c])); }
  rG = sm_SbTMulFast(sq[0], sq[1]);
}

// Ruling (a)'s certificate on a patch G1 to its faces crA (along v = 0) and crB (v = 1): II >= k I
// (k < 0), II with the unit normal toward crE's side. With Phi = w^3 |Pu x Pv| and any Q > 0,
// Phi >= 2 Q Phi^2 / (Q^2 + Phi^2) ((Phi - Q)^2 >= 0), equal where Q = Phi. Q = s ot / g, ot =
// det(H, Hu, Hv, (V, 0)) (V, g: sm_SbNormalField), certified one-signed (s): Q = Phi along a column
// where the patch is G1 to its face. Cleared of (ot^2 + g^2 Phi^2) w^2 > 0, with S = w^2 Phi^2 =
// J0 J2 - J1^2 (J the net's first fundamental form times w^4):
// s [D w^3 (l, m, n) + 2 |k| ot g S J] PSD, D = w^2 ot^2 + g^2 S, certified as s C11, s C22,
// C11 C22 - C12^2 >= 0. Proof: C-design.md, "The normal-field bound". 1 certified, 0 a negative
// value found, -1 undecided (ot not certified one-signed counts as undecided).
int sm_SbNoMoreConcave(const sm_SbNet & crN, const SmBlendSetback::Face & crA, const SmBlendSetback::Face & crB, const SmVector3d & crE, double dK)
{
  sm_SbT sF[4];
  sm_SbFundamentals(crN, crE, sF);
  sm_SbT sV[3], sG;
  sm_SbNormalField(crN, crA, crB, crE, sV, sG);
  const sm_SbT sOt = sm_SbTDetWith(crN, sV);
  sm_SbT sNegOt = sOt;
  for (size_t i=0; i<sNegOt.c.size(); i++) { sNegOt.c[i] = -sNegOt.c[i]; }
  const double dOs = 1.0e-12 * sm_SbTScale(sOt);
  const double s = (sm_SbTNonneg(sOt, dOs, SM_SB_GRAPH_DEPTH) == 1) ? 1.0 : (sm_SbTNonneg(sNegOt, dOs, SM_SB_GRAPH_DEPTH) == 1) ? -1.0 : 0.0;
  if (s == 0.0)
    { return -1; }
  sm_SbT h[4], hu[4], hv[4];
  for (int c=0; c<4; c++)
    {
      h[c] = sm_SbTMake(crN.p, crN.q);
      for (size_t i=0; i<crN.P.size(); i++) { h[c].c[i] = crN.P[i].a[c]; }
      hu[c] = sm_SbTDu(h[c]);
      hv[c] = sm_SbTDv(h[c]);
    }
  sm_SbT sXu[3], sXv[3], sJ[3];
  for (int c=0; c<3; c++)
    {
      sXu[c] = sm_SbTAxpy(sm_SbTMulFast(h[3], hu[c]), -1.0, sm_SbTMulFast(hu[3], h[c]));
      sXv[c] = sm_SbTAxpy(sm_SbTMulFast(h[3], hv[c]), -1.0, sm_SbTMulFast(hv[3], h[c]));
    }
  for (int k=0; k<3; k++)
    {
      const sm_SbT * pA = (k == 2) ? sXv : sXu;
      const sm_SbT * pB = (k == 0) ? sXu : sXv;
      sJ[k] = sm_SbTMulFast(pA[0], pB[0]);
      for (int c=1; c<3; c++) { sJ[k] = sm_SbTAxpy(sJ[k], 1.0, sm_SbTMulFast(pA[c], pB[c])); }
    }
  const sm_SbT sW2 = sm_SbTMulFast(h[3], h[3]);
  const sm_SbT sS = sm_SbTAxpy(sm_SbTMulFast(sJ[0], sJ[2]), -1.0, sm_SbTMulFast(sJ[1], sJ[1]));
  const sm_SbT sD = sm_SbTAxpy(sm_SbTMulFast(sW2, sm_SbTMulFast(sOt, sOt)), 1.0, sm_SbTMulFast(sm_SbTMulFast(sG, sG), sS));
  const sm_SbT sDw3 = sm_SbTMulFast(sD, sm_SbTMulFast(sW2, h[3]));
  const sm_SbT sLift = sm_SbTMulFast(sm_SbTMulFast(sOt, sG), sS);
  // each entry's rounding is that of its two terms (where II + |k| I is singular they cancel)
  sm_SbT sCm[3];
  double dTerm[3];
  for (int k=0; k<3; k++)
    {
      const sm_SbT sA = sm_SbTMulFast(sDw3, sF[k]), sB = sm_SbTMulFast(sLift, sJ[k]);
      dTerm[k] = sm_SbTScale(sA) + 2.0 * smos_Fabs(dK) * sm_SbTScale(sB);
      sCm[k] = sm_SbTAxpy(sA, -2.0 * dK, sB);
    }
  sm_SbT sL = sCm[0], sNn = sCm[2];
  for (size_t i=0; i<sL.c.size(); i++) { sL.c[i] *= s; }
  for (size_t i=0; i<sNn.c.size(); i++) { sNn.c[i] *= s; }
  const sm_SbT sDisc = sm_SbTAxpy(sm_SbTMulFast(sCm[0], sCm[2]), -1.0, sm_SbTMulFast(sCm[1], sCm[1]));
  const sm_SbT * pChecks[3] = { &sL, &sNn, &sDisc };
  const double dEps[3] = { 1.0e-12 * dTerm[0], 1.0e-12 * dTerm[2], 1.0e-12 * (dTerm[0] * dTerm[2] + dTerm[1] * dTerm[1]) };
  int iResult = 1;
  for (int k=0; k<3; k++)
    {
      const int iNon = sm_SbTNonneg(*pChecks[k], dEps[k], SM_SB_GRAPH_DEPTH);
      if (iNon == 0) { return 0; }
      iResult = smos_Min(iResult, iNon);
    }
  return iResult;
}

// Whether a constant round between faces crA (its row's start crEnd[0] on it) and crB (crEnd[1])
// has a line for its spine, by the one classification (sm_ClassifyBlendSpine): the faces as their
// natural quadrics, offset through the ball's centre (the row's circumcentre).
SmStatus sm_SbLineSpine(const SmBlendSetback::Face & crA, const SmBlendSetback::Face & crB, const SmPoint3d (&crEnd)[2], const std::vector<sm_SbH> & crRow,
                        double dTol, bool & rbLine)
{
  rbLine = false;
  const size_t q = crRow.size() - 1;
  sm_SbH sMid = sm_SbZero();
  for (size_t j=0; j<=q; j++) { sMid = sm_SbAxpy(sMid, sm_Binomial(static_cast<ULONG>(q), static_cast<ULONG>(j)) / std::pow(2.0, static_cast<double>(q)), crRow[j]); }
  const SmVector3d A(crEnd[0]), B(sMid.a[0] / sMid.a[3], sMid.a[1] / sMid.a[3], sMid.a[2] / sMid.a[3]), Cc(crEnd[1]);
  const SmVector3d a = A - Cc, b = B - Cc, axb = a * b;
  const double dAxb2 = axb.Dot(axb);
  if (!(dAxb2 > 0.0))
    { return SM_SUCCESS; }
  const SmPoint3d sC(Cc + ((b * a.Dot(a) - a * b.Dot(b)) * axb) * (0.5 / dAxb2));
  SmExactBlend sBlend;
  sBlend.m_dRadius = (A - SmVector3d(sC)).Length();
  for (int e=0; e<2; e++)
    {
      const SmBlendSetback::Face & crF = (e == 0) ? crA : crB;
      SmBlendQuadric & rQ = sBlend.m_aFaces[e];
      switch (crF.m_eKind)
        {
        case SmBlendSetback::SM_SB_PLANE:    rQ.m_eKind = SmBlendQuadric::SM_BQ_PLANE; rQ.m_sPoint = crF.m_sPoint; rQ.m_sDir = crF.m_sNormal; break;
        case SmBlendSetback::SM_SB_SPHERE:   rQ.m_eKind = SmBlendQuadric::SM_BQ_SPHERE; rQ.m_sPoint = crF.m_sPoint; rQ.m_dRadius = crF.m_dRadius; rQ.m_dSide = crF.m_dSide; break;
        case SmBlendSetback::SM_SB_CYLINDER: rQ.m_eKind = SmBlendQuadric::SM_BQ_CYLINDER; rQ.m_sPoint = crF.m_sPoint; rQ.m_sDir = crF.m_sAxis; rQ.m_dRadius = crF.m_dRadius; rQ.m_dSide = crF.m_dSide; break;
        default:
          {
            // the apex where R + slope h = 0, the axis into the nappe the face is on
            rQ.m_eKind = SmBlendQuadric::SM_BQ_CONE;
            rQ.m_sPoint = SmPoint3d(SmVector3d(crF.m_sPoint) - crF.m_sAxis * (crF.m_dRadius / crF.m_dSlope));
            rQ.m_sDir = crF.m_sAxis * ((crF.m_dSlope > 0.0) ? 1.0 : -1.0);
            rQ.m_dAngle = std::atan(smos_Fabs(crF.m_dSlope));
            rQ.m_dSide = crF.m_dSide;
            break;
          }
        }
      const SmVector3d sN = rQ.Normal(crEnd[e]);
      const double dOff = (SmVector3d(sC) - SmVector3d(crEnd[e])).Dot(sN);
      sBlend.m_adContact[e] = (dOff < 0.0) ? 1.0 : -1.0;
      sBlend.m_aOffsets[e] = rQ.Offset(dOff, dTol);
    }
  SER(sm_ClassifyBlendSpine(sC, dTol, sBlend));
  rbLine = (sBlend.m_eSpine == SM_BS_LINE);
  return SM_SUCCESS;
}

// A round's least normal curvature along its end row (owner's ruling (a), the rounds' part of k),
// II with the unit normal toward crE's side, to rounding. crJ: its jets there (sm_SbRoundJets, 0 ..
// 2), the row a circular arc (a constant round's section about its ball's centre C, radius r: its
// unit normal (P - C) / r, exactly, so |Pu x Pv| is polynomial there). On the degree 2 net with
// those jets at u = 0 (P0 = H, P1 = H + Hs / 2, P2 = H + Hs + Hss / 2) restricted to u = 0:
// k is certified a lower bound where s w^4 r (l, m, n) - k sigma Dc J is PSD on the row (sigma Dc =
// |det(H, Hs, Ht, (H_xyz - w C, 0))| = w^4 r |Hs x Ht|-scaled, s the side of o); the largest such k
// (bisection to 1e-13 of it, each step certified) is taken, 0 where the row is certified not
// concave (an end on a plane face exactly 0 along its rail, crPlaneEnd; a line spine's round, a
// cylinder, 0 along its rails everywhere, bLineSpine). Refused
// where the row is not a circular arc or o or Dc is not one-signed on it.
SmStatus sm_SbRoundRowK(const std::vector< std::vector<sm_SbH> > & crJ, const SmVector3d & crE, const bool (&crPlaneEnd)[2], bool bLineSpine,
                        double & rdK, const TCHAR *& rpWhy)
{
  rdK = 0.0;
  if (crJ.size() < 3)
    { SER(SM_ERR); }
  const size_t q = crJ[0].size() - 1;
  sm_SbNet sN;
  sN.p = 2;
  sN.q = q;
  sN.P.assign(3 * (q + 1), sm_SbZero());
  for (size_t j=0; j<=q; j++)
    for (int c=0; c<4; c++)
      {
        const double h = crJ[0][j].a[c], hs = crJ[1][j].a[c], hss = crJ[2][j].a[c];
        sN.P[j].a[c] = h;
        sN.P[(q + 1) + j].a[c] = h + 0.5 * hs;
        sN.P[2 * (q + 1) + j].a[c] = h + hs + 0.5 * hss;
      }
  // the arc's centre: the circumcentre of its ends and its middle
  auto at = [&](double t) {
    sm_SbH sH = sm_SbZero();
    for (size_t j=0; j<=q; j++)
      sH = sm_SbAxpy(sH, sm_Binomial(static_cast<ULONG>(q), static_cast<ULONG>(j)) * std::pow(t, static_cast<double>(j)) * std::pow(1.0 - t, static_cast<double>(q - j)), crJ[0][j]);
    return SmVector3d(sH.a[0] / sH.a[3], sH.a[1] / sH.a[3], sH.a[2] / sH.a[3]);
  };
  const SmVector3d A = at(0.0), B = at(0.5), Cc = at(1.0);
  const SmVector3d a = A - Cc, b = B - Cc, axb = a * b;
  const double dAxb2 = axb.Dot(axb);
  if (!(dAxb2 > 0.0))
    { SM_SB_FAIL(rpWhy, SM_ERR_INVALID_INPUT, _T("SmBlendSetback: a round's end section is not a circular arc")); }
  const SmVector3d sC = Cc + ((b * a.Dot(a) - a * b.Dot(b)) * axb) * (0.5 / dAxb2);
  const double r = (A - sC).Length();
  sm_SbT sU[3], sW = sm_SbTMake(2, q);
  for (size_t i=0; i<sN.P.size(); i++) { sW.c[i] = sN.P[i].a[3]; }
  for (int c=0; c<3; c++)
    {
      const double dC = (c == 0) ? sC.x : (c == 1) ? sC.y : sC.z;
      sU[c] = sm_SbTMake(2, q);
      for (size_t i=0; i<sN.P.size(); i++) { sU[c].c[i] = sN.P[i].a[c] - sN.P[i].a[3] * dC; }
    }
  auto row = [&](const sm_SbT & crT) {
    sm_SbT sR = sm_SbTMake(0, crT.q);
    for (size_t j=0; j<=crT.q; j++) { sR.c[j] = crT.c[j]; }
    return sR;
  };
  {
    // on the sphere about C of radius r along the row (the section a circular arc)
    sm_SbT sQ = sm_SbTMake(0, 0);
    for (int c=0; c<3; c++) { sQ = sm_SbTAxpy(sQ, 1.0, sm_SbTMulFast(row(sU[c]), row(sU[c]))); }
    const sm_SbT sR2 = sm_SbTMulFast(row(sW), row(sW));
    const double dScale = sm_SbTScale(sQ) + r * r * sm_SbTScale(sR2);
    sQ = sm_SbTAxpy(sQ, -r * r, sR2);
    if (sm_SbTScale(sQ) > 1.0e-11 * dScale)
      { SM_SB_FAIL(rpWhy, SM_ERR_INVALID_INPUT, _T("SmBlendSetback: a round's end section is not a circular arc")); }
  }
  sm_SbT sF[4];
  sm_SbFundamentals(sN, crE, sF);
  const sm_SbT sDc = row(sm_SbTDetWith(sN, sU));
  const sm_SbT sO = row(sF[3]);
  auto sign = [&](const sm_SbT & crT) {
    sm_SbT sNeg = crT;
    for (size_t i=0; i<sNeg.c.size(); i++) { sNeg.c[i] = -sNeg.c[i]; }
    const double dEps = 1.0e-12 * sm_SbTScale(crT);
    return (sm_SbTNonneg(crT, dEps, SM_SB_GRAPH_DEPTH) == 1) ? 1.0 : (sm_SbTNonneg(sNeg, dEps, SM_SB_GRAPH_DEPTH) == 1) ? -1.0 : 0.0;
  };
  const double s = sign(sO), sigma = sign(sDc);
  if (s == 0.0 || sigma == 0.0)
    { SM_SB_FAIL(rpWhy, SM_ERR_INVALID_INPUT, _T("SmBlendSetback: a round's end section turns across the corner's axis")); }
  sm_SbT h[4], hu[4], hv[4];
  for (int c=0; c<4; c++)
    {
      h[c] = sm_SbTMake(2, q);
      for (size_t i=0; i<sN.P.size(); i++) { h[c].c[i] = sN.P[i].a[c]; }
      hu[c] = sm_SbTDu(h[c]);
      hv[c] = sm_SbTDv(h[c]);
    }
  sm_SbT sXu[3], sXv[3], sA[3], sB[3];
  for (int c=0; c<3; c++)
    {
      sXu[c] = sm_SbTAxpy(sm_SbTMulFast(h[3], hu[c]), -1.0, sm_SbTMulFast(hu[3], h[c]));
      sXv[c] = sm_SbTAxpy(sm_SbTMulFast(h[3], hv[c]), -1.0, sm_SbTMulFast(hv[3], h[c]));
    }
  const sm_SbT sW4 = row(sm_SbTMulFast(sm_SbTMulFast(h[3], h[3]), sm_SbTMulFast(h[3], h[3])));
  double dTerm[3];
  for (int k=0; k<3; k++)
    {
      const sm_SbT * pA = (k == 2) ? sXv : sXu;
      const sm_SbT * pB = (k == 0) ? sXu : sXv;
      sm_SbT sJk = sm_SbTMulFast(pA[0], pB[0]);
      for (int c=1; c<3; c++) { sJk = sm_SbTAxpy(sJk, 1.0, sm_SbTMulFast(pA[c], pB[c])); }
      sA[k] = sm_SbTMulFast(sW4, row(sF[k]));
      for (size_t i=0; i<sA[k].c.size(); i++) { sA[k].c[i] *= s * r; }
      sB[k] = sm_SbTMulFast(sDc, row(sJk));
      for (size_t i=0; i<sB[k].c.size(); i++) { sB[k].c[i] *= sigma; }
      dTerm[k] = sm_SbTScale(sA[k]) + sm_SbTScale(sB[k]) / (r > 0.0 ? r : 1.0);
    }
  // an end of the row on a PLANE face: there the round is G1 to the plane along its rail, so its
  // second fundamental form on the rail's direction is the plane's, zero (Meusnier): its (rail,
  // rail) and (rail, across) entries are 0 at that end exactly, set so (the end's Bernstein
  // coefficient is the value there)
  for (int e=0; e<2; e++)
    if (crPlaneEnd[e])
      for (int k=0; k<2; k++)
        { sA[k].c[(e == 0) ? 0 : sA[k].c.size() - 1] = 0.0; }
  // a round whose spine is a line (sm_ClassifyBlendSpine: two planes; a plane or cylinder parallel
  // to a cylinder's axis) is a cylinder, its rails its rulings: those entries are 0 all along
  if (bLineSpine)
    for (int k=0; k<2; k++)
      for (size_t i=0; i<sA[k].c.size(); i++) { sA[k].c[i] = 0.0; }
  // 1 where sA - k sB is certified PSD on the row
  auto certified = [&](double dK) {
    sm_SbT sCm[3];
    double dT[3];
    for (int k=0; k<3; k++)
      {
        sCm[k] = sm_SbTAxpy(sA[k], -dK, sB[k]);
        dT[k] = sm_SbTScale(sA[k]) + smos_Fabs(dK) * sm_SbTScale(sB[k]);
      }
    const sm_SbT sDisc = sm_SbTAxpy(sm_SbTMulFast(sCm[0], sCm[2]), -1.0, sm_SbTMulFast(sCm[1], sCm[1]));
    return sm_SbTNonneg(sCm[0], 1.0e-12 * dT[0], SM_SB_GRAPH_DEPTH) == 1 && sm_SbTNonneg(sCm[2], 1.0e-12 * dT[2], SM_SB_GRAPH_DEPTH) == 1
           && sm_SbTNonneg(sDisc, 1.0e-12 * (dT[0] * dT[2] + dT[1] * dT[1]), SM_SB_GRAPH_DEPTH) == 1;
  };
  (void)dTerm;
  if (certified(0.0))
    { return SM_SUCCESS; }
  double dLo = -1.0 / r, dHi = 0.0;
  for (int i=0; i<64 && !certified(dLo); i++) { dHi = dLo; dLo *= 2.0; }
  if (!certified(dLo))
    { SM_SB_FAIL(rpWhy, SM_ERR, _T("SmBlendSetback: a round's least normal curvature at its end cannot be bounded")); }
  while (dHi - dLo > 1.0e-13 * smos_Fabs(dLo))
    {
      const double dMid = 0.5 * (dLo + dHi);
      if (certified(dMid)) { dLo = dMid; } else { dHi = dMid; }
    }
  rdK = dLo;
  return SM_SUCCESS;
}

/*******************************************************************//**
PURPOSE: The G1 fill of a corner beside curved faces, as a graph over the
    corner's axis (SmBlendSetback.h; #85, the owner's rulings of
    2026-11-02).

NOTES: As sm_SbGraphG1, with three changes where a face is curved (every
    face of such a corner is taken this way, planes too):
    - Each face's spring curve is SmBlendSetback::SpringCurve's, exactly on
      the face; its half from round m's rail end is setback patch m's side
      on that face (a fixed column). It is matched to the round there by a
      scale c, the Mobius weight ratio rho and the round's speed lambda off
      its end (rho = w0 / (alpha - kappa mu), lambda = d c rho kappa, from
      the half's second control point alpha (P, 1) + kappa nu' and the
      round's along-derivative mu (P, 1) + nu': both along the rail), then
      multiplied at the base point by the cubic [1, 1, q2, q3] that makes
      its last two weights 1 and keeps its first two control points (so the
      polynomial interior map meets it there): degree + 3.
    - The column next to it is G1 to the face by construction (Hahn's
      connecting functions): column 1 = column 0 + (alpha T1 + beta T2 +
      gamma H0) / q, T1 and T2 the face's parametric partials along the
      spring curve (SpringCurve), alpha, beta polynomials, gamma linear;
      chosen first (its rows 0, 1 the round's, its projections and weights
      at rows N - 1, N the interior map's, the rest the least second
      differences), then fixed. The interior map's heights there must equal
      the fixed columns': equalities, as a plane's were.
    - The convexity conditions are owner's ruling (a): each patch no more
      concave than what it is G1 to, II + |k| I PSD, k the most negative
      normal curvature of its faces (0 where all are convex); sufficient,
      and affine in the heights: w^3 s (l, m, n) + |k| o J PSD, J the first
      fundamental form of the projection times w^4 (I >= I of the
      projection; w^4 |Pu x Pv| >= w |o|).
    The patches' degree across is the larger of their two sides' (a
    cylinder or sphere side 4 + 3, a cone side 6 + 3, a plane side 2 + 3),
    plus one for gamma.
***********************************************************************/
SmStatus sm_SbGraphG1Curved(const SmContext & crContext, const std::vector<SmBlendSetback::Face> & crWorldFaces,
                            const std::vector<SmBlendSetback::Round> & crRounds, const std::vector< std::vector< std::vector<sm_SbH> > > & crWorldJets,
                            const std::vector<size_t> & crQ, const std::vector<ULONG> & crRound, const std::vector<ULONG> & crFace,
                            const std::vector<SmBoolean> & crRev, double dTol, SmBlendSetback::Fill & rFill, bool bLeastFirst)
{
  typedef SmBlendSetback::Face Face;
  const size_t n = crRound.size();
  // the frame about the rails' ends' centre (rounding of the corner's size, not its place)
  SmVector3d sO(0.0, 0.0, 0.0);
  {
    size_t lCount = 0;
    for (size_t r=0; r<crWorldJets.size(); r++)
      for (int e=0; e<2; e++)
        {
          const sm_SbH & crH = crWorldJets[r][0][e == 0 ? 0 : crQ[r]];
          sO = sO + SmVector3d(crH.a[0] / crH.a[3], crH.a[1] / crH.a[3], crH.a[2] / crH.a[3]);
          lCount++;
        }
    sO = sO * (1.0 / static_cast<double>(lCount));
  }
  std::vector< std::vector< std::vector<sm_SbH> > > crJets(crWorldJets);
  for (size_t r=0; r<crJets.size(); r++)
    for (size_t k=0; k<crJets[r].size(); k++)
      for (size_t c=0; c<crJets[r][k].size(); c++)
        {
          sm_SbH & rH = crJets[r][k][c];
          rH.a[0] -= rH.a[3] * sO.x;
          rH.a[1] -= rH.a[3] * sO.y;
          rH.a[2] -= rH.a[3] * sO.z;
        }
  std::vector<Face> crFaces(crWorldFaces);
  for (size_t f=0; f<crFaces.size(); f++) { crFaces[f].m_sPoint = SmPoint3d(SmVector3d(crFaces[f].m_sPoint) - sO); }

  // each cycle patch's round rows, its rails' ends and their directions toward the corner
  std::vector< std::vector<sm_SbH> > sH0(n), sH1(n);
  std::vector<SmPoint3d> sEnd0(n), sEnd1(n), sMid(n);
  std::vector<SmVector3d> sDir0(n), sDir1(n);
  for (size_t m=0; m<n; m++)
    {
      const size_t r = crRound[m], q = crQ[r];
      const double dSign = crRounds[r].m_bCornerAbove ? 1.0 : -1.0;
      sH0[m].resize(q + 1);
      sH1[m].resize(q + 1);
      for (size_t c=0; c<=q; c++)
        {
          const size_t cc = crRev[m] ? q - c : c;
          sH0[m][c] = crJets[r][0][cc];
          sH1[m][c] = sm_SbAxpy(sm_SbZero(), dSign, crJets[r][1][cc]);
        }
      for (int e=0; e<2; e++)
        {
          const sm_SbH & crA = sH0[m][e == 0 ? 0 : q];
          const sm_SbH & crB = sH1[m][e == 0 ? 0 : q];
          const SmPoint3d sP(crA.a[0] / crA.a[3], crA.a[1] / crA.a[3], crA.a[2] / crA.a[3]);
          if (smos_Fabs(SmBlendSetback::Distance(crFaces[crFace[(m + static_cast<size_t>(e)) % n]], sP)) > dTol)
            { SM_SB_FAIL(rFill.m_pReason, SM_ERR_INVALID_INPUT, _T("SmBlendSetback: a round does not meet its faces at its end")); }
          SmVector3d sT((crB.a[0] - sP.x * crB.a[3]) / crA.a[3], (crB.a[1] - sP.y * crB.a[3]) / crA.a[3], (crB.a[2] - sP.z * crB.a[3]) / crA.a[3]);
          SER(sT.Unitize());
          (e == 0 ? sEnd0[m] : sEnd1[m]) = sP;
          (e == 0 ? sDir0[m] : sDir1[m]) = sT;
        }
      sm_SbH sMh = sm_SbZero();
      for (size_t c=0; c<=q; c++)
        { sMh = sm_SbAxpy(sMh, sm_Binomial(static_cast<ULONG>(q), static_cast<ULONG>(c)) / std::pow(2.0, static_cast<double>(q)), sH0[m][c]); }
      sMid[m] = SmPoint3d(sMh.a[0] / sMh.a[3], sMh.a[1] / sMh.a[3], sMh.a[2] / sMh.a[3]);
    }
  // convex or concave, and the axis e into the convex side
  SmVector3d sSum(0.0, 0.0, 0.0);
  int iConvex = 0;
  for (size_t m=0; m<n; m++)
    {
      const Face & crF = crFaces[crFace[m]];
      sSum = sSum + crF.m_sNormal;
      const int iSide = (SmBlendSetback::Distance(crF, sMid[m]) < 0.0) ? 1 : -1;
      if (iConvex != 0 && iSide != iConvex)
        { SM_SB_FAIL(rFill.m_pReason, SM_ERR_INVALID_INPUT, _T("SmBlendSetback: a convex G1 fill needs rounds all convex or all concave")); }
      iConvex = iSide;
    }
  SmVector3d sE = sSum * (iConvex > 0 ? -1.0 : 1.0);
  SER(sE.Unitize());
  SmVector3d sA1(1.0, 0.0, 0.0);
  if (smos_Fabs(sE.y) < smos_Fabs(sE.x) && smos_Fabs(sE.y) <= smos_Fabs(sE.z)) { sA1 = SmVector3d(0.0, 1.0, 0.0); }
  else if (smos_Fabs(sE.z) < smos_Fabs(sE.x) && smos_Fabs(sE.z) < smos_Fabs(sE.y)) { sA1 = SmVector3d(0.0, 0.0, 1.0); }
  sA1 = sA1 - sE * sA1.Dot(sE);
  SER(sA1.Unitize());
  const SmVector3d sA2 = sE * sA1;
  auto inPlane = [&](const SmVector3d & crV) { return SmVector2d(crV.Dot(sA1), crV.Dot(sA2)); };

  // the faces' least normal curvature on the fill's convex side (ruling (a)): a cylinder or sphere
  // whose material is on the other side than a convex body's, -1 / R; such a cone, -cos / rho over
  // the spring curves' least radius (bounded below by its control points' axial coordinates)
  std::vector<double> sFaceK(n, 0.0);
  // each round's least normal curvature along its end row (the rounds' part of k)
  std::vector<double> sRoundK(n, 0.0);
  for (size_t r=0; r<n; r++)
    {
      const TCHAR * pWhy = NULL;
      // the faces at the row's ends: which are planes, and the round's spine by the one
      // classification (sm_ClassifyBlendSpine) from the faces' kinds and placement
      bool bPlaneEnd[2] = { false, false };
      SmPoint3d sEnd[2];
      size_t lFaceAt[2] = { crFaces.size(), crFaces.size() };
      for (int e=0; e<2; e++)
        {
          const sm_SbH & crH = crJets[r][0][e == 0 ? 0 : crQ[r]];
          sEnd[e] = SmPoint3d(crH.a[0] / crH.a[3], crH.a[1] / crH.a[3], crH.a[2] / crH.a[3]);
          for (size_t f=0; f<crFaces.size(); f++)
            if (smos_Fabs(SmBlendSetback::Distance(crFaces[f], sEnd[e])) <= dTol)
              {
                lFaceAt[e] = f;
                bPlaneEnd[e] = (crFaces[f].m_eKind == SmBlendSetback::SM_SB_PLANE);
              }
        }
      bool bLine = false;
      if (lFaceAt[0] < crFaces.size() && lFaceAt[1] < crFaces.size())
        SER(sm_SbLineSpine(crFaces[lFaceAt[0]], crFaces[lFaceAt[1]], sEnd, crJets[r][0], dTol, bLine));
      if (sm_SbRoundRowK(crJets[r], sE, bPlaneEnd, bLine, sRoundK[r], pWhy) != SM_SUCCESS)
        { SM_SB_FAIL(rFill.m_pReason, SM_ERR_INVALID_INPUT, pWhy != NULL ? pWhy : _T("SmBlendSetback: a round's least normal curvature at its end cannot be found")); }
    }

  // the spring curves: per face m, patch m's half (its v = 0 column) and patch m - 1's (its v = 1
  // column), matched to their rounds and to the interior map; their tangent-plane directions
  struct Side
  {
    std::vector<sm_SbH> m_vCol;    // the column, degree m_lDeg
    size_t m_lDeg = 0;
    std::vector<sm_SbH> m_vT[2];  // the face's partials along it (Mobius applied)
    double m_dLambda = 0.0;        // the round's speed off its end there
    SmVector2d m_sD;               // the interior map's field there (the column's derivative at the base point, negated, in the plane)
  };
  std::vector<Side> sSide0(n), sSide1(n); // patch m's at v = 0 (face m), at v = 1 (face m + 1)
  std::vector<SmPoint2d> sXi(n);
  for (size_t m=0; m<n; m++)
    {
      const size_t mp = (m + n - 1) % n;
      const Face & crF = crFaces[crFace[m]];
      std::vector<double> sHP, sHQ, sTP, sTQ;
      const TCHAR * pWhy = NULL;
      const SmStatus eSpring = SmBlendSetback::SpringCurve(crF, sEnd0[m], sDir0[m], sEnd1[mp], -sDir1[mp], sHP, sHQ, &pWhy, &sTP, &sTQ);
      if (eSpring != SM_SUCCESS)
        { SM_SB_FAIL(rFill.m_pReason, eSpring, pWhy != NULL ? pWhy : _T("SmBlendSetback: a face's spring curve cannot be made")); }
      // the face's curvature bound
      if (crF.m_eKind != SmBlendSetback::SM_SB_PLANE && ((crF.m_dSide < 0.0) == (iConvex > 0)))
        {
          if (crF.m_eKind != SmBlendSetback::SM_SB_CONE)
            { sFaceK[crFace[m]] = -1.0 / crF.m_dRadius; }
          else
            {
              double dHMin = SM_BIG_DOUBLE;
              for (const std::vector<double> * pH : { &sHP, &sHQ })
                for (size_t i=0; i<pH->size() / 4; i++)
                  {
                    const SmVector3d sX((*pH)[4 * i] / (*pH)[4 * i + 3], (*pH)[4 * i + 1] / (*pH)[4 * i + 3], (*pH)[4 * i + 2] / (*pH)[4 * i + 3]);
                    dHMin = smos_Min(dHMin, (sX - SmVector3d(crF.m_sPoint)).Dot(crF.m_sAxis));
                  }
              const double dRho = crF.m_dRadius + crF.m_dSlope * dHMin;
              if (!(dRho > 0.0))
                { SM_SB_FAIL(rFill.m_pReason, SM_ERR_INVALID_INPUT, _T("SmBlendSetback: a spring curve reaches a cone's apex")); }
              sFaceK[crFace[m]] = -1.0 / (dRho * sqrt(1.0 + crF.m_dSlope * crF.m_dSlope));
            }
        }
      for (int s=0; s<2; s++)
        {
          const std::vector<double> & crH = (s == 0) ? sHP : sHQ;
          const std::vector<double> & crT = (s == 0) ? sTP : sTQ;
          const size_t lPatch = (s == 0) ? m : mp;
          Side & rSide = (s == 0) ? sSide0[m] : sSide1[mp];
          const size_t ds = crH.size() / 4 - 1;
          std::vector<sm_SbH> sC(ds + 1);
          for (size_t i=0; i<=ds; i++) for (int c=0; c<4; c++) { sC[i].a[c] = crH[4 * i + c]; }
          // the round there: its rail end H_r0 and its along-derivative H_r1
          const size_t q = crQ[crRound[lPatch]];
          const sm_SbH & crR0 = sH0[lPatch][s == 0 ? 0 : q];
          const sm_SbH & crR1 = sH1[lPatch][s == 0 ? 0 : q];
          const SmPoint3d sP(crR0.a[0] / crR0.a[3], crR0.a[1] / crR0.a[3], crR0.a[2] / crR0.a[3]);
          const double w0 = sC[0].a[3], dAlpha = sC[1].a[3];
          const SmVector3d sB(sC[1].a[0] - dAlpha * sP.x, sC[1].a[1] - dAlpha * sP.y, sC[1].a[2] - dAlpha * sP.z);
          const double dMu = crR1.a[3];
          const SmVector3d sNu(crR1.a[0] - dMu * sP.x, crR1.a[1] - dMu * sP.y, crR1.a[2] - dMu * sP.z);
          const double dKappa = sB.Dot(sNu) / sNu.Dot(sNu); // the half's tangent part over the round's (both along the rail)
          const double dRho = w0 / (dAlpha - dKappa * dMu), dC = crR0.a[3] / w0;
          const double dLambda = static_cast<double>(ds) * dC * dRho * dKappa;
          if (!(dKappa > 0.0) || !(dRho > 0.0) || !(dLambda > 0.0))
            { SM_SB_FAIL(rFill.m_pReason, SM_ERR_INVALID_INPUT, _T("SmBlendSetback: a spring curve cannot be matched to its round's end")); }
          double dPow = 1.0;
          for (size_t i=0; i<=ds; i++, dPow *= dRho)
            for (int c=0; c<4; c++) { sC[i].a[c] *= dC * dPow; }
          // the base point's weights 1, 1: times [1, 1, q2, q3]
          const double dW1 = sC[ds].a[3], dDW1 = static_cast<double>(ds) * (sC[ds].a[3] - sC[ds - 1].a[3]);
          const double q3 = 1.0 / dW1, q2 = q3 * (1.0 + dDW1 / (3.0 * dW1));
          if (!(q2 > 0.0))
            { SM_SB_FAIL(rFill.m_pReason, SM_ERR_INVALID_INPUT, _T("SmBlendSetback: a spring curve's weight would not stay positive")); }
          const sm_SbB sFactor = { 1.0, 1.0, q2, q3 };
          rSide.m_lDeg = ds + 3;
          rSide.m_vCol.assign(ds + 4, sm_SbZero());
          for (int c=0; c<4; c++)
            {
              sm_SbB sComp(ds + 1);
              for (size_t i=0; i<=ds; i++) { sComp[i] = sC[i].a[c]; }
              const sm_SbB sProd = sm_SbMul(sComp, sFactor);
              for (size_t i=0; i<=ds+3; i++) { rSide.m_vCol[i].a[c] = sProd[i]; }
            }
          rSide.m_vCol[ds + 3].a[3] = 1.0;
          rSide.m_vCol[ds + 2].a[3] = 1.0;
          // the partials, the same reparametrization (rho^i; a direction's scale is free)
          const size_t k1 = (crF.m_eKind == SmBlendSetback::SM_SB_PLANE) ? 0 : (crF.m_eKind == SmBlendSetback::SM_SB_SPHERE) ? 2 : 4;
          const size_t k2 = (crF.m_eKind == SmBlendSetback::SM_SB_PLANE || crF.m_eKind == SmBlendSetback::SM_SB_CYLINDER) ? 0 : k1;
          if (crT.size() != 4 * (k1 + 1) + 4 * (k2 + 1))
            { SER(SM_ERR); }
          for (int t=0; t<2; t++)
            {
              const size_t kk = (t == 0) ? k1 : k2, lOff = (t == 0) ? 0 : 4 * (k1 + 1);
              rSide.m_vT[t].assign(kk + 1, sm_SbZero());
              double dP = 1.0;
              for (size_t i=0; i<=kk; i++, dP *= dRho)
                for (int c=0; c<4; c++) { rSide.m_vT[t][i].a[c] = crT[lOff + 4 * i + c] * dP; }
            }
          rSide.m_dLambda = dLambda;
          // the field: minus the column's derivative at the base point, in the plane
          const sm_SbH & crL = rSide.m_vCol[ds + 3], & crL1 = rSide.m_vCol[ds + 2];
          const SmVector3d sDer = SmVector3d(crL.a[0] - crL1.a[0], crL.a[1] - crL1.a[1], crL.a[2] - crL1.a[2]) * static_cast<double>(ds + 3);
          rSide.m_sD = inPlane(sDer * -1.0);
          if (s == 0)
            {
              const SmVector2d sXi2 = inPlane(SmVector3d(crL.a[0], crL.a[1], crL.a[2]));
              sXi[m] = SmPoint2d(sXi2.x, sXi2.y);
            }
        }
    }
  SmPoint2d sCentre(0.0, 0.0);
  double dRho = 0.0;
  for (size_t m=0; m<n; m++) { sCentre = SmPoint2d(sCentre.x + sXi[m].x / static_cast<double>(n), sCentre.y + sXi[m].y / static_cast<double>(n)); }
  for (size_t m=0; m<n; m++) { dRho = smos_Max(dRho, SmVector2d(sXi[m] - sCentre).Length()); }
  if (!(dRho > 0.0))
    { SER(SM_ERR); }
  std::vector<SmPoint2d> sB(n);
  for (size_t m=0; m<n; m++) { sB[m] = SmPoint2d((sXi[m].x - sCentre.x) / dRho, (sXi[m].y - sCentre.y) / dRho); }
  auto heightOf = [&](const sm_SbH & crH) { return crH.a[0] * sE.x + crH.a[1] * sE.y + crH.a[2] * sE.z; };

  // the degrees: 9 first (proven unmet there, unmet at every degree), then up from the least
  std::vector<sm_SbNet> sTop;
  size_t lTopDegree = 0;
  std::vector<size_t> sDegrees;
  if (!bLeastFirst) { sDegrees.push_back(9); }
  for (size_t d=3; d<9; d++) { sDegrees.push_back(d); }
  if (bLeastFirst) { sDegrees.push_back(9); }
  for (size_t lD=0; lD<sDegrees.size(); lD++)
    {
      const size_t d = sDegrees[lD];
      const std::vector< std::pair<int, int> > sMono = sm_SbMonomials(static_cast<int>(d));
      const size_t M = sMono.size();
      std::vector<size_t> sNc(n), sQv(n), sFirstFree(n);
      size_t lUnknowns = M;
      for (size_t m=0; m<n; m++)
        {
          sNc[m] = smos_Max(sSide0[m].m_lDeg, sSide1[m].m_lDeg) + 3;
          sQv[m] = smos_Max(smos_Max(crQ[crRound[m]] + 1, d), static_cast<size_t>(5));
          sFirstFree[m] = lUnknowns;
          lUnknowns += (sNc[m] - 3) * (sQv[m] - 3); // rows 2 .. Nc - 2 of columns 2 .. qv - 2
        }
      std::vector<sm_SbAffinePatch> sPatches(2 * n);
      std::vector<double> sAeq, sBeq;
      bool bFailed = false;
      for (size_t m=0; m<n && !bFailed; m++)
        {
          const size_t Nc = sNc[m], qv = sQv[m];
          // rows 0, 1: the round's, at its speeds off its end (linear across)
          std::vector< std::vector<sm_SbH> > sRows(Nc + 1, std::vector<sm_SbH>(qv + 1, sm_SbZero()));
          const sm_SbB sLam = { sSide0[m].m_dLambda / static_cast<double>(Nc), sSide1[m].m_dLambda / static_cast<double>(Nc) };
          for (int c=0; c<4; c++)
            {
              sm_SbB s0(sH0[m].size()), s1(sH0[m].size());
              for (size_t j=0; j<s0.size(); j++) { s0[j] = sH0[m][j].a[c]; s1[j] = sH1[m][j].a[c]; }
              s1 = sm_SbMul(sLam, s1);
              s0 = sm_SbElevate(s0, qv);
              s1 = sm_SbElevate(s1, qv);
              for (size_t j=0; j<=qv; j++) { sRows[0][j].a[c] = s0[j]; sRows[1][j].a[c] = s0[j] + s1[j]; }
            }
          // rows Nc - 1, Nc: the interior map's projection along edge m and its field there
          const SmPoint2d & crX0 = sXi[m];
          const SmPoint2d & crX1 = sXi[(m + 1) % n];
          const SmVector2d sD0 = sSide0[m].m_sD, sD1 = sSide1[m].m_sD;
          for (int c=0; c<3; c++)
            {
              const double dA1 = (c == 0) ? sA1.x : (c == 1) ? sA1.y : sA1.z;
              const double dA2 = (c == 0) ? sA2.x : (c == 1) ? sA2.y : sA2.z;
              const sm_SbB sP = { crX0.x * dA1 + crX0.y * dA2, crX1.x * dA1 + crX1.y * dA2 };
              const sm_SbB sDP = { sD0.x * dA1 + sD0.y * dA2, sD1.x * dA1 + sD1.y * dA2 };
              const sm_SbB s3 = sm_SbElevate(sP, qv), sDe = sm_SbElevate(sDP, qv);
              for (size_t j=0; j<=qv; j++)
                {
                  sRows[Nc][j].a[c] = s3[j];
                  sRows[Nc - 1][j].a[c] = s3[j] + sDe[j] / static_cast<double>(Nc);
                }
            }
          for (size_t j=0; j<=qv; j++) { sRows[Nc - 1][j].a[3] = 1.0; sRows[Nc][j].a[3] = 1.0; }
          // rows 2 .. Nc - 2: the cubic with the same end rows and derivatives, raised to Nc (the
          // planar fill's own shape across; their heights free in columns 2 .. qv - 2)
          for (size_t j=0; j<=qv; j++)
            for (int c=0; c<4; c++)
              {
                const double dK = static_cast<double>(Nc) / 3.0;
                const sm_SbB sCubic = { sRows[0][j].a[c], sRows[0][j].a[c] + dK * (sRows[1][j].a[c] - sRows[0][j].a[c]),
                                        sRows[Nc][j].a[c] + dK * (sRows[Nc - 1][j].a[c] - sRows[Nc][j].a[c]), sRows[Nc][j].a[c] };
                const sm_SbB sUp = sm_SbElevate(sCubic, Nc);
                for (size_t i=2; i+2<=Nc; i++) { sRows[i][j].a[c] = sUp[i]; }
              }
          // the fixed columns: the spring halves (0 and qv; rows 0, 1 the round's) and the columns
          // G1 to the faces by construction (1 and qv - 1)
          std::vector<sm_SbH> sFixedHeight[2][2]; // [side][col 0 or 1]: rows Nc - 1, Nc targets
          for (int e=0; e<2 && !bFailed; e++)
            {
              const Side & crS = (e == 0) ? sSide0[m] : sSide1[m];
              const size_t jEdge = (e == 0) ? 0 : qv, jNext = (e == 0) ? 1 : qv - 1;
              // the side's column at Nc - 1 (for gamma's product) and at Nc
              std::vector<sm_SbH> sColLow(Nc - 2), sCol(Nc + 1);
              for (int c=0; c<4; c++)
                {
                  sm_SbB sComp(crS.m_vCol.size());
                  for (size_t i=0; i<sComp.size(); i++) { sComp[i] = crS.m_vCol[i].a[c]; }
                  const sm_SbB sLow = sm_SbElevate(sComp, Nc - 3), sHigh = sm_SbElevate(sComp, Nc);
                  for (size_t i=0; i<=Nc-3; i++) { sColLow[i].a[c] = sLow[i]; }
                  for (size_t i=0; i<=Nc; i++) { sCol[i].a[c] = sHigh[i]; }
                }
              sCol[0] = sRows[0][jEdge];
              sCol[1] = sRows[1][jEdge];
              for (size_t i=0; i+1<Nc; i++) { sRows[i][jEdge] = sCol[i]; }
              sFixedHeight[e][0] = { sCol[Nc - 1], sCol[Nc] };
              // the constructive column: unknowns alpha (Nc - k1 + 1), beta (Nc - k2 + 1), gamma (2)
              const size_t k1 = crS.m_vT[0].size() - 1, k2 = crS.m_vT[1].size() - 1;
              if (k1 > Nc || k2 > Nc)
                { bFailed = true; break; }
              const size_t lA = Nc - k1 + 1, lBt = Nc - k2 + 1, lU = lA + lBt + 4;
              // each unknown's contribution to the column (Bernstein of degree Nc, homogeneous)
              std::vector< std::vector<sm_SbH> > sBasis(lU, std::vector<sm_SbH>(Nc + 1, sm_SbZero()));
              for (size_t u=0; u<lU; u++)
                for (int c=0; c<4; c++)
                  {
                    sm_SbB sUnit, sOther;
                    if (u < lA)
                      {
                        sUnit.assign(lA, 0.0); sUnit[u] = 1.0;
                        sOther.resize(k1 + 1);
                        for (size_t i=0; i<=k1; i++) { sOther[i] = crS.m_vT[0][i].a[c]; }
                      }
                    else if (u < lA + lBt)
                      {
                        sUnit.assign(lBt, 0.0); sUnit[u - lA] = 1.0;
                        sOther.resize(k2 + 1);
                        for (size_t i=0; i<=k2; i++) { sOther[i] = crS.m_vT[1][i].a[c]; }
                      }
                    else
                      {
                        sUnit.assign(4, 0.0); sUnit[u - lA - lBt] = 1.0;
                        sOther.resize(Nc - 2);
                        for (size_t i=0; i<=Nc-3; i++) { sOther[i] = sColLow[i].a[c]; }
                      }
                    const sm_SbB sProd = sm_SbElevate(sm_SbMul(sUnit, sOther), Nc);
                    for (size_t i=0; i<=Nc; i++) { sBasis[u][i].a[c] = sProd[i] / static_cast<double>(qv); }
                  }
              // the conditions: rows 0, 1 the round's (4 each); rows Nc - 1, Nc the interior map's
              // projection and weight (3 each)
              std::vector<double> sA, sb;
              auto addRow = [&](size_t i, const double (&rW)[4], double dTarget)
                {
                  for (size_t u=0; u<lU; u++)
                    {
                      double s = 0.0;
                      for (int c=0; c<4; c++) { s += rW[c] * sBasis[u][i].a[c]; }
                      sA.push_back(s);
                    }
                  double dBase = 0.0;
                  for (int c=0; c<4; c++) { dBase += rW[c] * sCol[i].a[c]; }
                  sb.push_back(dTarget - dBase);
                };
              for (size_t i=0; i<2; i++)
                for (int c=0; c<4; c++)
                  {
                    double w[4] = { 0.0, 0.0, 0.0, 0.0 };
                    w[c] = 1.0;
                    addRow(i, w, sRows[i][jNext].a[c]);
                  }
              for (size_t i=Nc-1; i<=Nc; i++)
                {
                  const double w1[4] = { sA1.x, sA1.y, sA1.z, 0.0 }, w2[4] = { sA2.x, sA2.y, sA2.z, 0.0 }, ww[4] = { 0.0, 0.0, 0.0, 1.0 };
                  addRow(i, w1, sRows[i][jNext].a[0] * sA1.x + sRows[i][jNext].a[1] * sA1.y + sRows[i][jNext].a[2] * sA1.z);
                  addRow(i, w2, sRows[i][jNext].a[0] * sA2.x + sRows[i][jNext].a[1] * sA2.y + sRows[i][jNext].a[2] * sA2.z);
                  addRow(i, ww, 1.0);
                }
              std::vector<double> sY0;
              std::vector< std::vector<double> > sZ0;
              double dSz = 0.0;
              for (size_t i=0; i<sb.size(); i++) { dSz = smos_Max(dSz, smos_Fabs(sb[i])); }
              if (!sm_SbAffineSolutions(sA, sb, sb.size(), lU, 1.0e-9 * (1.0 + dSz), sY0, sZ0))
                {
                  bFailed = true; break;
                }
              // the rest: column 1's middle rows as near as they can be to the interpolated rows the
              // free columns start from (their projections and weights)
              if (!sZ0.empty())
                {
                  const size_t lF = sZ0.size();
                  std::vector<double> sMrow, sCrow;
                  for (size_t i=2; i+2<=Nc; i++)
                    {
                      const double wP1[4] = { sA1.x, sA1.y, sA1.z, 0.0 }, wP2[4] = { sA2.x, sA2.y, sA2.z, 0.0 }, wW[4] = { 0.0, 0.0, 0.0, 1.0 };
                      const double * apW[3] = { wP1, wP2, wW };
                      for (int g=0; g<3; g++)
                        {
                          double dTarget = 0.0, dBase = 0.0;
                          for (int c=0; c<4; c++) { dTarget += apW[g][c] * sRows[i][jNext].a[c]; dBase += apW[g][c] * sCol[i].a[c]; }
                          std::vector<double> sCoef(lU, 0.0);
                          for (size_t u=0; u<lU; u++)
                            for (int c=0; c<4; c++) { sCoef[u] += apW[g][c] * sBasis[u][i].a[c]; }
                          double dAtY0 = dBase;
                          for (size_t u=0; u<lU; u++) { dAtY0 += sCoef[u] * sY0[u]; }
                          for (size_t f=0; f<lF; f++)
                            {
                              double t = 0.0;
                              for (size_t u=0; u<lU; u++) { t += sCoef[u] * sZ0[f][u]; }
                              sMrow.push_back(t);
                            }
                          sCrow.push_back(dTarget - dAtY0);
                        }
                    }
                  // and a least nudge in every free direction (the rows may leave some free)
                  for (size_t g=0; g<lF; g++)
                    {
                      for (size_t f=0; f<lF; f++) { sMrow.push_back((f == g) ? 1.0e-9 : 0.0); }
                      sCrow.push_back(0.0);
                    }
                  std::vector<double> sT;
                  if (!sm_SbLeastSquares(sMrow, sCrow, sCrow.size(), lF, 1.0e-14, sT))
                    { bFailed = true; break; }
                  for (size_t f=0; f<lF; f++)
                    for (size_t u=0; u<lU; u++) { sY0[u] += sZ0[f][u] * sT[f]; }
                }
              std::vector<sm_SbH> sNext(Nc + 1);
              for (size_t i=0; i<=Nc; i++)
                {
                  sNext[i] = sCol[i];
                  for (size_t u=0; u<lU; u++) { sNext[i] = sm_SbAxpy(sNext[i], sY0[u], sBasis[u][i]); }
                }
              sNext[0] = sRows[0][jNext];
              sNext[1] = sRows[1][jNext];
              for (size_t i=0; i+1<Nc; i++) { sRows[i][jNext] = sNext[i]; }
              sFixedHeight[e][1] = { sNext[Nc - 1], sNext[Nc] };
            }
          if (bFailed)
            { break; }
          // h's monomials along the edge and along the field (rows Nc - 1, Nc)
          const sm_SbB sX = { sB[m].x, sB[(m + 1) % n].x }, sY = { sB[m].y, sB[(m + 1) % n].y };
          const sm_SbB sDX = { sD0.x / dRho, sD1.x / dRho }, sDY = { sD0.y / dRho, sD1.y / dRho };
          std::vector< std::vector<double> > sZ2(M), sZ3(M);
          for (size_t j=0; j<M; j++)
            {
              const int a = sMono[j].first, b = sMono[j].second;
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
              const sm_SbB v0 = sm_SbElevate(sm_SbAlong(a, b, sX, sY, d), qv), v1 = sm_SbElevate(s1, qv);
              sZ3[j] = v0;
              sZ2[j].resize(qv + 1);
              for (size_t i=0; i<=qv; i++) { sZ2[j][i] = v0[i] + v1[i] / static_cast<double>(Nc); }
            }
          sm_SbAffinePatch & rP = sPatches[m];
          rP.m_sBase.p = Nc;
          rP.m_sBase.q = qv;
          rP.m_sBase.P.assign((Nc + 1) * (qv + 1), sm_SbZero());
          for (size_t i=0; i<=Nc; i++)
            for (size_t j=0; j<=qv; j++) { rP.m_sBase.P[i * (qv + 1) + j] = sRows[i][j]; }
          for (size_t jm=0; jm<M; jm++)
            {
              sm_SbT sZ = sm_SbTMake(Nc, qv);
              for (size_t j=0; j<=qv; j++)
                {
                  sZ.c[(Nc - 1) * (qv + 1) + j] = sZ2[jm][j];
                  sZ.c[Nc * (qv + 1) + j] = sZ3[jm][j];
                }
              rP.m_vIndex.push_back(jm);
              rP.m_vZ.push_back(sZ);
            }
          // the free middle heights
          size_t lAt = sFirstFree[m];
          for (size_t i=2; i+2<=Nc; i++)
            for (size_t j=2; j+2<=qv; j++)
              {
                sm_SbT sZ = sm_SbTMake(Nc, qv);
                sZ.c[i * (qv + 1) + j] = 1.0;
                rP.m_vIndex.push_back(lAt++);
                rP.m_vZ.push_back(sZ);
              }
          // the equalities: the interior map's heights at rows Nc - 1, Nc of the fixed columns
          for (int e=0; e<2; e++)
            for (int jj=0; jj<2; jj++)
              {
                const size_t col = (e == 0) ? static_cast<size_t>(jj) : qv - static_cast<size_t>(jj);
                for (size_t k=0; k<2; k++)
                  {
                    const size_t row = Nc - 1 + k, idx = row * (qv + 1) + col;
                    std::vector<double> sRowEq(lUnknowns, 0.0);
                    for (size_t t=0; t<rP.m_vIndex.size(); t++) { sRowEq[rP.m_vIndex[t]] += rP.m_vZ[t].c[idx]; }
                    sAeq.insert(sAeq.end(), sRowEq.begin(), sRowEq.end());
                    sBeq.push_back(heightOf(sFixedHeight[e][jj][k]) - heightOf(rP.m_sBase.P[idx]));
                  }
              }
        }
      if (bFailed)
        { continue; }
      // the interior parts, as sm_SbGraphG1's
      for (size_t q=0; q<n; q++)
        {
          const SmPoint2d sC(0.0, 0.0);
          const SmPoint2d sMq((sB[q].x + sB[(q + 1) % n].x) * 0.5, (sB[q].y + sB[(q + 1) % n].y) * 0.5);
          const SmPoint2d sMp((sB[q].x + sB[(q + n - 1) % n].x) * 0.5, (sB[q].y + sB[(q + n - 1) % n].y) * 0.5);
          sm_SbT sXq, sYq;
          sXq.p = sXq.q = sYq.p = sYq.q = 1;
          sXq.c = { sC.x, sMp.x, sMq.x, sB[q].x };
          sYq.c = { sC.y, sMp.y, sMq.y, sB[q].y };
          sm_SbAffinePatch & rP = sPatches[n + q];
          rP.m_sBase.p = rP.m_sBase.q = d;
          rP.m_sBase.P.assign((d + 1) * (d + 1), sm_SbZero());
          const sm_SbT sXd = sm_SbTElevate(sXq, d, d), sYd = sm_SbTElevate(sYq, d, d);
          for (size_t i=0; i<(d + 1) * (d + 1); i++)
            {
              const double xi1 = sCentre.x + dRho * sXd.c[i], xi2 = sCentre.y + dRho * sYd.c[i];
              rP.m_sBase.P[i].a[0] = xi1 * sA1.x + xi2 * sA2.x;
              rP.m_sBase.P[i].a[1] = xi1 * sA1.y + xi2 * sA2.y;
              rP.m_sBase.P[i].a[2] = xi1 * sA1.z + xi2 * sA2.z;
              rP.m_sBase.P[i].a[3] = 1.0;
            }
          for (size_t j=0; j<M; j++)
            {
              rP.m_vIndex.push_back(j);
              rP.m_vZ.push_back(sm_SbOver(sMono[j].first, sMono[j].second, sXq, sYq, d));
            }
        }
      std::vector<double> sYp;
      std::vector< std::vector<double> > sZn;
      double dSize = 0.0;
      for (size_t i=0; i<sBeq.size(); i++) { dSize = smos_Max(dSize, smos_Fabs(sBeq[i])); }
      if (!sm_SbAffineSolutions(sAeq, sBeq, sBeq.size(), lUnknowns, 1.0e-9 * (1.0 + dSize), sYp, sZn) || sZn.empty())
        { continue; }
      const size_t lFree = sZn.size();

      // the energy (as sm_SbGraphG1's)
      auto energyRows = [&](const sm_SbNet & crN, std::vector<double> & rRows)
        {
          rRows.clear();
          for (int g1=0; g1<8; g1++)
            for (int g2=0; g2<8; g2++)
              {
                double R[3][3];
                sm_SbSecond(crN, SM_SB_GX[g1], SM_SB_GX[g2], R);
                const double w = smos_Sqrt(SM_SB_GW[g1] * SM_SB_GW[g2]);
                for (int c=0; c<3; c++)
                  {
                    rRows.push_back(w * R[0][c]);
                    rRows.push_back(w * smos_Sqrt(2.0) * R[1][c]);
                    rRows.push_back(w * R[2][c]);
                  }
              }
        };
      std::vector<double> sE0;
      std::vector<size_t> sEStart;
      std::vector< std::vector<double> > sPatchE0(2 * n);
      for (size_t p=0; p<2*n; p++)
        {
          energyRows(sPatches[p].m_sBase, sPatchE0[p]);
          sEStart.push_back(sE0.size());
          sE0.insert(sE0.end(), sPatchE0[p].begin(), sPatchE0[p].end());
        }
      const size_t lE = sE0.size();
      std::vector<double> sAE(lE * lUnknowns, 0.0);
      {
        std::vector<double> sUnit(lUnknowns, 0.0), sRowsE;
        for (size_t p=0; p<2*n; p++)
          for (size_t t=0; t<sPatches[p].m_vIndex.size(); t++)
            {
              const size_t u = sPatches[p].m_vIndex[t];
              sUnit[u] = 1.0;
              energyRows(sm_SbAffineNet(sPatches[p], sUnit, sE), sRowsE);
              sUnit[u] = 0.0;
              for (size_t i=0; i<sRowsE.size(); i++) { sAE[(sEStart[p] + i) * lUnknowns + u] += sRowsE[i] - sPatchE0[p][i]; }
            }
      }
      std::vector<double> sBm(lE * lFree, 0.0), sCv(lE, 0.0);
      for (size_t i=0; i<lE; i++)
        {
          double s = sE0[i];
          for (size_t u=0; u<lUnknowns; u++) { s += sAE[i * lUnknowns + u] * sYp[u]; }
          sCv[i] = s;
          for (size_t f=0; f<lFree; f++)
            {
              double t = 0.0;
              for (size_t u=0; u<lUnknowns; u++) { t += sAE[i * lUnknowns + u] * sZn[f][u]; }
              sBm[i * lFree + f] = t;
            }
        }
      auto netsAt = [&](const std::vector<double> & crT, std::vector<sm_SbNet> & rNets)
        {
          std::vector<double> sYv(sYp);
          for (size_t f=0; f<lFree; f++)
            for (size_t u=0; u<lUnknowns; u++) { sYv[u] += sZn[f][u] * crT[f]; }
          rNets.assign(2 * n, sm_SbNet());
          for (size_t p=0; p<2*n; p++) { rNets[p] = sm_SbAffineNet(sPatches[p], sYv, sE); }
        };

      // each patch's ruling-(a) matrix: w^3 (l, m, n) + |k| o J (J the projection's w^4 I), affine
      // in the heights; o's sign certified
      std::vector< std::array<sm_SbT, 3> > sQ0(2 * n);
      std::vector< std::array<std::vector<sm_SbT>, 3> > sQj(2 * n);
      std::vector<double> sS(2 * n, 0.0), sK(2 * n, 0.0);
      for (size_t m=0; m<n; m++)
        {
          sK[m] = smos_Min(smos_Min(sFaceK[crFace[m]], sFaceK[crFace[(m + 1) % n]]), sRoundK[crRound[m]]);
          sK[n + m] = smos_Min(sFaceK[crFace[m]], smos_Min(sK[m], sK[(m + n - 1) % n]));
        }
      for (size_t p=0; p<2*n; p++)
        {
          sm_SbT sQ0p[3], sOp;
          std::vector<sm_SbT> sQjp[3];
          sm_SbAffineFundamentals(sPatches[p], sE, sQ0p, sQjp, sOp);
          sm_SbT sNegO = sOp;
          for (size_t i=0; i<sNegO.c.size(); i++) { sNegO.c[i] = -sNegO.c[i]; }
          const double dO = 1.0e-12 * sm_SbTScale(sOp);
          sS[p] = (sm_SbTNonneg(sOp, dO, SM_SB_GRAPH_DEPTH) == 1) ? 1.0 : (sm_SbTNonneg(sNegO, dO, SM_SB_GRAPH_DEPTH) == 1) ? -1.0 : 0.0;
          if (sS[p] == 0.0)
            {
              SM_SB_FAIL(rFill.m_pReason, SM_ERR, _T("SmBlendSetback: the fill folds (a patch's normal turns across the corner's axis)"));
            }
          // w^3 and the projection's J
          const sm_SbNet & crN = sPatches[p].m_sBase;
          sm_SbT h[4], hu[4], hv[4];
          for (int c=0; c<4; c++)
            {
              h[c] = sm_SbTMake(crN.p, crN.q);
              for (size_t i=0; i<crN.P.size(); i++)
                {
                  double dV = crN.P[i].a[c];
                  if (c < 3)
                    {
                      const double dH = heightOf(crN.P[i]);
                      dV -= dH * ((c == 0) ? sE.x : (c == 1) ? sE.y : sE.z);
                    }
                  h[c].c[i] = dV;
                }
              hu[c] = sm_SbTDu(h[c]);
              hv[c] = sm_SbTDv(h[c]);
            }
          const sm_SbT sW2 = sm_SbTMulFast(h[3], h[3]), sW3 = sm_SbTMulFast(sW2, h[3]);
          sm_SbT sJ[3];
          {
            sm_SbT sXu[3], sXv[3];
            for (int c=0; c<3; c++)
              {
                sXu[c] = sm_SbTAxpy(sm_SbTMulFast(h[3], hu[c]), -1.0, sm_SbTMulFast(hu[3], h[c]));
                sXv[c] = sm_SbTAxpy(sm_SbTMulFast(h[3], hv[c]), -1.0, sm_SbTMulFast(hv[3], h[c]));
              }
            for (int k=0; k<3; k++)
              {
                const sm_SbT * pA = (k == 2) ? sXv : sXu;
                const sm_SbT * pB = (k == 0) ? sXu : sXv;
                sJ[k] = sm_SbTMulFast(pA[0], pB[0]);
                for (int c=1; c<3; c++) { sJ[k] = sm_SbTAxpy(sJ[k], 1.0, sm_SbTMulFast(pA[c], pB[c])); }
              }
          }
          for (int k=0; k<3; k++)
            {
              // (where k = 0, w^3 > 0 changes nothing but the sufficient conditions: left out)
              sQ0[p][k] = (sK[p] < 0.0) ? sm_SbTMulFast(sW3, sQ0p[k]) : sQ0p[k];
              if (sK[p] < 0.0)
                { sQ0[p][k] = sm_SbTAxpy(sQ0[p][k], -sK[p], sm_SbTMulFast(sOp, sJ[k])); }
              sQj[p][k].clear();
              for (size_t t=0; t<sQjp[k].size(); t++) { sQj[p][k].push_back((sK[p] < 0.0) ? sm_SbTMulFast(sW3, sQjp[k][t]) : sQjp[k][t]); }
            }
        }
      // the QP (as sm_SbGraphG1's, without the base points' kernel directions: the faces' columns
      // are fixed here)
      std::vector<double> sG(lFree * lFree, 0.0), sg(lFree, 0.0);
      double dTrace = 0.0;
      for (size_t a=0; a<lFree; a++)
        {
          for (size_t b=a; b<lFree; b++)
            {
              double s = 0.0;
              for (size_t i=0; i<lE; i++) { s += sBm[i * lFree + a] * sBm[i * lFree + b]; }
              sG[a * lFree + b] = sG[b * lFree + a] = 2.0 * s;
            }
          double s = 0.0;
          for (size_t i=0; i<lE; i++) { s += sBm[i * lFree + a] * sCv[i]; }
          sg[a] = 2.0 * s;
          dTrace += sG[a * lFree + a];
        }
      for (size_t a=0; a<lFree; a++) { sG[a * lFree + a] += 1.0e-12 * dTrace / static_cast<double>(lFree); }
      // each interior part's extra direction at its base point (as sm_SbGraphG1's): the field there
      std::vector<double> sKernel(2 * n, -1.0);
      for (size_t q=0; q<n; q++)
        {
          const SmVector2d sDs((sB[q].x - sB[(q + n - 1) % n].x) * 0.5, (sB[q].y - sB[(q + n - 1) % n].y) * 0.5);
          const SmVector2d sDt((sB[q].x - sB[(q + 1) % n].x) * 0.5, (sB[q].y - sB[(q + 1) % n].y) * 0.5);
          const double dDet = sDs.x * sDt.y - sDs.y * sDt.x;
          if (!(smos_Fabs(dDet) > 0.0))
            { SER(SM_ERR); }
          const SmVector2d & crSig = sSide0[q].m_sD;
          const double a = (crSig.x * sDt.y - crSig.y * sDt.x) / dDet, b = (sDs.x * crSig.y - sDs.y * crSig.x) / dDet;
          double th = std::atan2(b, a);
          if (th < 0.0) { th += SM_PI; }
          if (th >= SM_PI) { th -= SM_PI; }
          sKernel[n + q] = th;
        }
      std::vector<double> sC, sb;
      for (size_t p=0; p<2*n; p++)
        {
          const sm_SbAffinePatch & crP = sPatches[p];
          std::vector<double> sTheta;
          for (size_t k=0; k<SM_SB_GRAPH_DIRECTIONS; k++) { sTheta.push_back(SM_PI * static_cast<double>(k) / static_cast<double>(SM_SB_GRAPH_DIRECTIONS)); }
          if (sKernel[p] >= 0.0) { sTheta.push_back(sKernel[p]); }
          std::sort(sTheta.begin(), sTheta.end());
          sTheta.push_back(sTheta[0] + SM_PI);
          size_t P = 0, Q = 0;
          for (int k=0; k<3; k++) { P = smos_Max(P, sQ0[p][k].p); Q = smos_Max(Q, sQ0[p][k].q); }
          std::vector< std::vector<sm_SbT> > sPc0(3), sPcf(3 * lFree);
          double dScale = 0.0;
          for (int k=0; k<3; k++)
            {
              sm_SbT s0 = sm_SbTElevate(sQ0[p][k], P, Q);
              std::vector<sm_SbT> sF(lFree, sm_SbTMake(P, Q));
              for (size_t t=0; t<crP.m_vIndex.size(); t++)
                {
                  const sm_SbT sJt = sm_SbTElevate(sQj[p][k][t], P, Q);
                  dScale = smos_Max(dScale, sm_SbTScale(sJt));
                  const size_t u = crP.m_vIndex[t];
                  for (size_t i=0; i<sJt.c.size(); i++) { s0.c[i] += sYp[u] * sJt.c[i]; }
                  for (size_t f=0; f<lFree; f++)
                    {
                      const double z = sZn[f][u];
                      if (z == 0.0) { continue; }
                      for (size_t i=0; i<sJt.c.size(); i++) { sF[f].c[i] += z * sJt.c[i]; }
                    }
                }
              dScale = smos_Max(dScale, sm_SbTScale(sQ0[p][k]));
              sm_SbTPieces(s0, SM_SB_GRAPH_LEVELS, sPc0[k]);
              for (size_t f=0; f<lFree; f++) { sm_SbTPieces(sF[f], SM_SB_GRAPH_LEVELS, sPcf[k * lFree + f]); }
            }
          std::vector<double> sRowT(lFree);
          for (size_t pc=0; pc<sPc0[0].size(); pc++)
            for (size_t i=0; i<sPc0[0][pc].c.size(); i++)
              for (size_t a=0; a+1<sTheta.size(); a++)
                for (int iPair=0; iPair<2; iPair++)
                  {
                    const double c0 = std::cos(sTheta[a]), s0 = std::sin(sTheta[a]);
                    const double c1 = (iPair == 0) ? c0 : std::cos(sTheta[a + 1]), s1 = (iPair == 0) ? s0 : std::sin(sTheta[a + 1]);
                    const double w[3] = { c0 * c1, c0 * s1 + s0 * c1, s0 * s1 };
                    double dConst = 0.0, dMax = 0.0;
                    for (int k=0; k<3; k++) { dConst += sS[p] * w[k] * sPc0[k][pc].c[i]; }
                    for (size_t f=0; f<lFree; f++)
                      {
                        double s = 0.0;
                        for (int k=0; k<3; k++) { s += sS[p] * w[k] * sPcf[k * lFree + f][pc].c[i]; }
                        sRowT[f] = s;
                        dMax = smos_Max(dMax, smos_Fabs(s));
                      }
                    if (!(dMax > 1.0e-11 * dScale))
                      { continue; }
                    for (size_t f=0; f<lFree; f++) { sC.push_back(sRowT[f] / dMax); }
                    sb.push_back((-dConst + 1.0e-12 * dScale) / dMax);
                  }
        }
      std::vector<double> sT;
      rFill.m_lSolves++;
      s_lCurvedDegrees++;
      s_lCurvedRows += sb.size();
      bool bProven = false;
      const SmStatus eQP = sm_SbDualQPRows(sG, sg, sC, sb, lFree, sb.size(), 1.0e-13, sT, &bProven);
      if (eQP != SM_SUCCESS && bProven && d == 9)
        { break; }
      if (eQP != SM_SUCCESS)
        { continue; }
      std::vector<sm_SbNet> sNets;
      netsAt(sT, sNets);
      // the certificate, on the nets themselves: o one-signed, s C11, s C22, C11 C22 - C12^2 >= 0
      bool bCertified = true;
      for (size_t p=0; p<2*n && bCertified; p++)
        {
          sm_SbT sF[4];
          sm_SbFundamentals(sNets[p], sE, sF);
          const double dOs = sm_SbTScale(sF[3]);
          sm_SbT sNegO = sF[3];
          for (size_t i=0; i<sNegO.c.size(); i++) { sNegO.c[i] = -sNegO.c[i]; }
          const double s = (sm_SbTNonneg(sF[3], 1.0e-12 * dOs, SM_SB_GRAPH_DEPTH) == 1) ? 1.0 : (sm_SbTNonneg(sNegO, 1.0e-12 * dOs, SM_SB_GRAPH_DEPTH) == 1) ? -1.0 : 0.0;
          if (s == 0.0) { bCertified = false; break; }
          // beside a face where k < 0: exact along the faces' columns (sm_SbNoMoreConcave)
          if (p < n && sK[p] < 0.0)
            {
              const int iNo = sm_SbNoMoreConcave(sNets[p], crFaces[crFace[p]], crFaces[crFace[(p + 1) % n]], sE, sK[p]);
              if (iNo != 1)
                {
                  bCertified = false;
                }
              continue;
            }
          // the shift's J from the net's projection (as above) and w^3
          sm_SbT h[4], hu[4], hv[4];
          for (int c=0; c<4; c++)
            {
              h[c] = sm_SbTMake(sNets[p].p, sNets[p].q);
              for (size_t i=0; i<sNets[p].P.size(); i++)
                {
                  double dV = sNets[p].P[i].a[c];
                  if (c < 3) { dV -= heightOf(sNets[p].P[i]) * ((c == 0) ? sE.x : (c == 1) ? sE.y : sE.z); }
                  h[c].c[i] = dV;
                }
              hu[c] = sm_SbTDu(h[c]);
              hv[c] = sm_SbTDv(h[c]);
            }
          const sm_SbT sW3 = sm_SbTMulFast(sm_SbTMulFast(h[3], h[3]), h[3]);
          sm_SbT sXu[3], sXv[3], sCm[3];
          for (int c=0; c<3; c++)
            {
              sXu[c] = sm_SbTAxpy(sm_SbTMulFast(h[3], hu[c]), -1.0, sm_SbTMulFast(hu[3], h[c]));
              sXv[c] = sm_SbTAxpy(sm_SbTMulFast(h[3], hv[c]), -1.0, sm_SbTMulFast(hv[3], h[c]));
            }
          for (int k=0; k<3; k++)
            {
              sCm[k] = (sK[p] < 0.0) ? sm_SbTMulFast(sW3, sF[k]) : sF[k];
              if (sK[p] < 0.0)
                {
                  const sm_SbT * pA = (k == 2) ? sXv : sXu;
                  const sm_SbT * pB = (k == 0) ? sXu : sXv;
                  sm_SbT sJk = sm_SbTMulFast(pA[0], pB[0]);
                  for (int c=1; c<3; c++) { sJk = sm_SbTAxpy(sJk, 1.0, sm_SbTMulFast(pA[c], pB[c])); }
                  sCm[k] = sm_SbTAxpy(sCm[k], -sK[p], sm_SbTMulFast(sF[3], sJk));
                }
            }
          sm_SbT sL = sCm[0], sNn = sCm[2];
          for (size_t i=0; i<sL.c.size(); i++) { sL.c[i] *= s; }
          for (size_t i=0; i<sNn.c.size(); i++) { sNn.c[i] *= s; }
          const sm_SbT sDisc = sm_SbTAxpy(sm_SbTMulFast(sCm[0], sCm[2]), -1.0, sm_SbTMulFast(sCm[1], sCm[1]));
          const sm_SbT * pChecks[3] = { &sL, &sNn, &sDisc };
          for (int k=0; k<3 && bCertified; k++)
            {
              const int iNon = sm_SbTNonneg(*pChecks[k], 1.0e-12 * sm_SbTScale(*pChecks[k]), SM_SB_GRAPH_DEPTH);
              if (iNon != 1)
                {
                  bCertified = false;
                }
            }
        }
      if (!bCertified)
        { continue; }
      // back to the world
      for (size_t p=0; p<sNets.size(); p++)
        for (size_t i=0; i<sNets[p].P.size(); i++)
          {
            sm_SbH & rH = sNets[p].P[i];
            rH.a[0] += rH.a[3] * sO.x;
            rH.a[1] += rH.a[3] * sO.y;
            rH.a[2] += rH.a[3] * sO.z;
          }
      if (d == 9)
        {
          sTop = sNets;
          lTopDegree = d;
          continue;
        }
      SER(sm_SbFillFrom(crContext, sNets, crRound, crFace, crRev, d, rFill));
      return SM_SUCCESS;
    }
  if (!sTop.empty())
    {
      SER(sm_SbFillFrom(crContext, sTop, crRound, crFace, crRev, lTopDegree, rFill));
      return SM_SUCCESS;
    }
  rFill.m_bUncertified = true;
  SM_SB_FAIL(rFill.m_pReason, SM_ERR, _T("SmBlendSetback: no G1 fill beside a curved face up to degree 9 is certified (ruling (a)) at these setbacks"));
}

} // namespace

/*******************************************************************//**
PURPOSE: The hole's fill by the setback split (SmBlendSetback.h).

NOTES: The cycle: round m between faces m (its setback patch's v = 0) and
    m + 1. A G1 fill is a graph over the corner's axis (sm_SbGraphG1); the
    rest is the G2 fill's. Its domain: the regular n-gon with vertices B_m = (cos 2 pi m / n,
    sin 2 pi m / n), B_m the base point of face m; its edge m from B_m to
    B_m+1 is the inner curve of setback patch m; at B_m the direction field
    of edge m starts as sigma_m = t_m + t_m-1 (t the edges' unit directions)
    and that of edge m - 1 ends as -sigma_m, so the two setback patches leave
    face m's base point along one line (the spring curve is tangent there).
***********************************************************************/
SmStatus sm_SbFillHole(const SmContext & crContext, const std::vector<SmBlendSetback::Face> & crFaces, const std::vector<SmBlendSetback::Round> & crRounds,
                       ULONG lOrder, double dTol, SmBlendSetback::Fill & rFill, SmBlendSetback::Probe * pProbe);

SmStatus SmBlendSetback::FillHoleProbe(const SmContext & crContext, const std::vector<Face> & crFaces, const std::vector<Round> & crRounds, double dTol,
                                       Fill & rFill, Probe & rProbe)
{
  return sm_SbFillHole(crContext, crFaces, crRounds, 1, dTol, rFill, &rProbe);
}

SmStatus SmBlendSetback::FillHole(const SmContext & crContext, const std::vector<Face> & crFaces, const std::vector<Round> & crRounds, ULONG lOrder,
                                  double dTol, Fill & rFill)
{
  return sm_SbFillHole(crContext, crFaces, crRounds, lOrder, dTol, rFill, NULL);
}

SmStatus sm_SbFillHole(const SmContext & crContext, const std::vector<SmBlendSetback::Face> & crFaces, const std::vector<SmBlendSetback::Round> & crRounds,
                       ULONG lOrder, double dTol, SmBlendSetback::Fill & rFill, SmBlendSetback::Probe * pProbe)
{
  typedef SmBlendSetback::Round Round;
  typedef SmBlendSetback::Face Face;
  for (size_t i=0; i<rFill.m_vPatches.size(); i++) { delete rFill.m_vPatches[i]; }
  rFill.m_vPatches.clear();
  rFill.m_vRound.clear();
  rFill.m_vFace.clear();
  rFill.m_vReversed.clear();
  rFill.m_lDegree = 0;
  rFill.m_lSolves = 0;
  rFill.m_bUncertified = false;
  rFill.m_pReason = NULL;
  if (lOrder != 1 && lOrder != 2)
    { SM_SB_FAIL(rFill.m_pReason, SM_ERR_INVALID_INPUT, _T("SmBlendSetback: a setback fill is of order 1 (G1) or 2 (G2)")); }
  if (!(dTol > 0.0) || crRounds.size() < 3 || crFaces.size() != crRounds.size())
    { SM_SB_FAIL(rFill.m_pReason, SM_ERR_INVALID_INPUT, _T("SmBlendSetback: a setback fill needs three or more rounds, as many faces, and a positive tolerance")); }
  const size_t n = crRounds.size();
  const ULONG k = lOrder;
  const size_t N = 2 * k + 1; // the setback patches' degree across the hole
  bool bAnyCurved = false; // beside a curved face the rounds' second derivatives along too (their k)
  for (size_t f=0; f<crFaces.size(); f++) { bAnyCurved = bAnyCurved || crFaces[f].m_eKind != SmBlendSetback::SM_SB_PLANE; }

  // each round's jets at its end, and the faces its across ends lie on
  std::vector< std::vector< std::vector<sm_SbH> > > sJets(n);
  std::vector<size_t> sQ(n);
  std::vector<double> sDtDs(n, 1.0);
  std::vector<size_t> sFaceAt[2];
  sFaceAt[0].assign(n, n);
  sFaceAt[1].assign(n, n);
  for (size_t r=0; r<n; r++)
    {
      SER(sm_SbRoundJets(crContext, crRounds[r], bAnyCurved ? smos_Max(k, static_cast<ULONG>(2)) : k, sJets[r], sQ[r], sDtDs[r], rFill.m_pReason));
      for (int e=0; e<2; e++)
        {
          const sm_SbH & crH = sJets[r][0][e == 0 ? 0 : sQ[r]];
          const SmPoint3d sP(crH.a[0] / crH.a[3], crH.a[1] / crH.a[3], crH.a[2] / crH.a[3]);
          for (size_t f=0; f<n; f++)
            if (smos_Fabs(SmBlendSetback::Distance(crFaces[f], sP)) <= dTol)
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

  bool bCurved = false;
  for (size_t f=0; f<n; f++) { bCurved = bCurved || crFaces[f].m_eKind != SmBlendSetback::SM_SB_PLANE; }
  if (bCurved)
    {
      if (k != 1 || (pProbe != NULL && !pProbe->m_bFill))
        { SM_SB_FAIL(rFill.m_pReason, SM_ERR_INVALID_INPUT, _T("SmBlendSetback: beside a curved face only a G1 fill is built, and probed only as the fill itself")); }
      s_lCurvedFills++;
      return sm_SbGraphG1Curved(crContext, crFaces, crRounds, sJets, sQ, sRound, sFace, sRev, dTol, rFill, pProbe != NULL);
    }
  if (k == 1)
    { return sm_SbGraphG1(crContext, crFaces, crRounds, sJets, sQ, sRound, sFace, sRev, dTol, rFill, pProbe); }
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
  // of its spoke (from its end section's middle to the middle of the end sections' middles).
  // Half that left saddles in the box corner's setback patches (Gaussian curvature down to
  // -0.82 at r 2, against the corner ball's 0.25); the whole spoke, -0.25
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
        sSpeed[m] = ((crR.m_eAlong == SM_SP_U) ? sD[2] : sD[1]).Length() * sDtDs[sRound[m]]; // by s (sm_SbRoundJets)
        sCentre = sCentre + sD[0] * (1.0 / static_cast<double>(n));
      }
    for (size_t m=0; m<n; m++)
      {
        const size_t r = sRound[m];
        const Round & crR = crRounds[r];
        if (!(sSpeed[m] > 0.0))
          { SM_SB_FAIL(rFill.m_pReason, SM_ERR_INVALID_INPUT, _T("SmBlendSetback: a round does not move along its rails at its end")); }
        double dLambda = SmVector3d(sMid[m] - SmPoint3d(sCentre)).Length() / sSpeed[m];
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

/*******************************************************************//**
PURPOSE: The least scale of a corner's setbacks at which its G1 fill is
    certified convex (SmBlendSetback.h).

NOTES: The corner's rounds are made from its geometry: edge i's round, of
    radius r between its faces f and g (normals n out of the material), has
    its ball's centre C at V + t a + s r (n_f + n_g) / (1 + n_f . n_g) at t
    along the edge from the vertex V (s = -1 at a convex edge, the ball in
    the material; +1 at a concave one), its section the arc from C - s r n_f
    to C - s r n_g about C, the rational quadratic through the edge's point
    with the weight cos of half the arc, (1 + n_f . n_g) / 2 under the root.
    A round runs from its end section away from the vertex, linear along.
***********************************************************************/
namespace
{

SmStatus sm_SbModelRounds(const SmContext & crContext, const std::vector<SmBlendSetback::Face> & crFaces, const std::vector<SmBlendSetback::Edge> & crEdges,
                          const std::vector<double> & crSetBacks, std::vector<SmBSplineSurface *> & rvSurfaces, std::vector<SmBlendSetback::Round> & rvRounds,
                          const TCHAR *& rpWhy)
{
  rvSurfaces.clear();
  rvRounds.clear();
  const size_t n = crEdges.size();
  if (n < 3 || crFaces.size() != n || crSetBacks.size() != n)
    { SM_SB_FAIL(rpWhy, SM_ERR_INVALID_INPUT, _T("SmBlendSetback: a corner needs three or more edges, as many faces and setbacks")); }
  // the vertex: the faces' planes' common point (least squares)
  double A[3][3] = {{0.0}}, b[3] = {0.0, 0.0, 0.0};
  for (size_t f=0; f<n; f++)
    {
      const SmVector3d & crN = crFaces[f].m_sNormal;
      const double c[3] = { crN.x, crN.y, crN.z }, d = crN.Dot(SmVector3d(crFaces[f].m_sPoint));
      for (int i=0; i<3; i++)
        {
          for (int j=0; j<3; j++) { A[i][j] += c[i] * c[j]; }
          b[i] += c[i] * d;
        }
    }
  const double dDet = A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) - A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0])
                      + A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);
  if (!(smos_Fabs(dDet) > 1.0e-12))
    { SM_SB_FAIL(rpWhy, SM_ERR_INVALID_INPUT, _T("SmBlendSetback: the corner's faces meet in no one point")); }
  double x[3];
  for (int k=0; k<3; k++)
    {
      double M[3][3];
      for (int i=0; i<3; i++)
        for (int j=0; j<3; j++) { M[i][j] = (j == k) ? b[i] : A[i][j]; }
      x[k] = (M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
              + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0])) / dDet;
    }
  const SmVector3d sV(x[0], x[1], x[2]);
  for (size_t i=0; i<n; i++)
    {
      const SmBlendSetback::Edge & crE = crEdges[i];
      const ULONG f = crE.m_alFaces[0], g = crE.m_alFaces[1];
      if (f >= n || g >= n || !(crE.m_dRadius > 0.0) || !(crSetBacks[i] > 0.0))
        { SM_SB_FAIL(rpWhy, SM_ERR_INVALID_INPUT, _T("SmBlendSetback: an edge needs a positive radius, a positive setback and two of the faces given")); }
      const SmVector3d & crNf = crFaces[f].m_sNormal, & crNg = crFaces[g].m_sNormal;
      // convex where the other edge on face g runs into face f's material
      double dSide = 0.0;
      for (size_t j=0; j<n; j++)
        if (j != i && (crEdges[j].m_alFaces[0] == g || crEdges[j].m_alFaces[1] == g))
          { dSide = crEdges[j].m_sAway.Dot(crNf); }
      if (dSide == 0.0)
        { SM_SB_FAIL(rpWhy, SM_ERR_INVALID_INPUT, _T("SmBlendSetback: an edge's face holds no other edge of the corner")); }
      const double s = (dSide < 0.0) ? -1.0 : 1.0, r = crE.m_dRadius, c = crNf.Dot(crNg);
      if (!(1.0 + c > 1.0e-9))
        { SM_SB_FAIL(rpWhy, SM_ERR_INVALID_INPUT, _T("SmBlendSetback: an edge's faces fold onto each other")); }
      const SmVector3d sW = (crNf + crNg) * (s * r / (1.0 + c));
      const double w = smos_Sqrt(0.5 * (1.0 + c));
      const double dLen = 2.0 * r;
      SmTArray<SmPoint3d> sPts;
      SmTArray<double> sWts;
      for (int k=0; k<2; k++)
        {
          const SmVector3d sE = sV + crE.m_sAway * (crSetBacks[i] + k * dLen);
          const SmVector3d sC = sE + sW;
          sPts.Add(SmPoint3d(sC - crNf * (s * r)));
          sPts.Add(SmPoint3d(sE));
          sPts.Add(SmPoint3d(sC - crNg * (s * r)));
          sWts.Add(1.0);
          sWts.Add(w);
          sWts.Add(1.0);
        }
      SmTArray<ULONG> sUM, sVM;
      SmTArray<double> sUK, sVK;
      sUM.Add(2); sUM.Add(2);
      sVM.Add(3); sVM.Add(3);
      sUK.Add(crSetBacks[i]); sUK.Add(crSetBacks[i] + dLen);
      sVK.Add(0.0); sVK.Add(1.0);
      SmBSplineSurface * pSurface = NULL;
      const SmStatus eStatus = SmBSplineSurface::CreateCanonical(crContext, 1, 2, sPts, SM_SF_UNSPECIFIED, sUM, sVM, sUK, sVK, SM_KT_UNSPECIFIED, &sWts, NULL,
                                                                 pSurface);
      if (eStatus != SM_SUCCESS || pSurface == NULL)
        {
          for (size_t o=0; o<rvSurfaces.size(); o++) { delete rvSurfaces[o]; }
          rvSurfaces.clear();
          rvRounds.clear();
          SER(eStatus != SM_SUCCESS ? eStatus : SM_ERR);
        }
      rvSurfaces.push_back(pSurface);
      SmBlendSetback::Round sRound = { pSurface, SM_SP_U, crSetBacks[i], FALSE };
      rvRounds.push_back(sRound);
    }
  return SM_SUCCESS;
}

} // namespace

SmStatus SmBlendSetback::CertifiedScale(const SmContext & crContext, const std::vector<Face> & crFaces, const std::vector<Edge> & crEdges,
                                        const std::vector<double> & crSetBacks, double & rdScale, ULONG * plSolves, const TCHAR ** ppReason,
                                        const std::vector<double> * pvMost)
{
  rdScale = 0.0;
  if (plSolves != NULL) { *plSolves = 0; }
  const TCHAR * pWhy = NULL;
  if (ppReason != NULL) { *ppReason = NULL; }
  double dSize = 0.0;
  for (size_t i=0; i<crSetBacks.size(); i++) { dSize = smos_Max(dSize, smos_Fabs(crSetBacks[i])); }
  if (pvMost != NULL && pvMost->size() != crSetBacks.size())
    { SM_SB_FAIL(*(ppReason != NULL ? ppReason : &pWhy), SM_ERR_INVALID_INPUT, _T("SmBlendSetback: one most setback per edge")); }
  // the setbacks at scale k: each k times its own, held at its most
  auto at = [&](double k)
    {
      std::vector<double> sAt(crSetBacks);
      for (size_t i=0; i<sAt.size(); i++)
        {
          sAt[i] *= k;
          if (pvMost != NULL) { sAt[i] = smos_Min(sAt[i], (*pvMost)[i]); }
        }
      return sAt;
    };
  // the model's fill at scale k: its probe at degree 9 (pProbe), or the fill itself (taken certified)
  auto model = [&](double k, Probe * pProbe, bool & rbCertified) -> SmStatus
    {
      rbCertified = false;
      std::vector<SmBSplineSurface *> sSurfaces;
      std::vector<Round> sRounds;
      SmStatus eStatus = sm_SbModelRounds(crContext, crFaces, crEdges, at(k), sSurfaces, sRounds, pWhy);
      if (eStatus != SM_SUCCESS)
        { return eStatus; }
      Fill sFill;
      const double dTol = 1.0e-9 * (1.0 + k * dSize);
      Probe sAsFill;
      sAsFill.m_bFill = true;
      eStatus = FillHoleProbe(crContext, crFaces, sRounds, dTol, sFill, (pProbe != NULL) ? *pProbe : sAsFill);
      if (plSolves != NULL) { *plSolves += sFill.m_lSolves; }
      for (size_t o=0; o<sSurfaces.size(); o++) { delete sSurfaces[o]; }
      if (pProbe == NULL && eStatus != SM_SUCCESS && sFill.m_bUncertified)
        { return SM_SUCCESS; }
      if (eStatus != SM_SUCCESS)
        {
          pWhy = sFill.m_pReason;
          return eStatus;
        }
      rbCertified = (pProbe == NULL);
      return SM_SUCCESS;
    };
  auto margin = [&](double k, double & rdMargin) -> SmStatus
    {
      Probe sProbe;
      sProbe.m_bMargin = true;
      bool bUnused = false;
      SER(model(k, &sProbe, bUnused));
      rdMargin = sProbe.m_dMargin;
      return SM_SUCCESS;
    };
  auto fail = [&](SmStatus eStatus) -> SmStatus
    {
      if (ppReason != NULL) { *ppReason = pWhy; }
      return eStatus;
    };
  // The search, from the margin m(k) by which the fill's sector conditions fail at degree 9 (m <= 0:
  // met, and the conditions only grow easier with k): m at 1 and 1.02, then at their secant's root,
  // then the root of the parabola through the three; the certifying end checked a half percent
  // above it, the failing end the bracket's 1e-2 below.
  double m1 = 0.0;
  SmStatus eStatus = margin(1.0, m1);
  if (eStatus != SM_SUCCESS) { return fail(eStatus); }
  bool bCertified = false;
  if (m1 <= 0.0)
    {
      eStatus = model(1.0, NULL, bCertified);
      if (eStatus != SM_SUCCESS) { return fail(eStatus); }
      if (bCertified)
        {
          rdScale = 1.0;
          return SM_SUCCESS;
        }
    }
  std::vector<double> sK(1, 1.0), sM(1, m1);
  double dLo = 1.0; // the largest scale known to fail
  double dPredict = 0.0;
  for (int iStep=0; iStep<3; iStep++)
    {
      double k = 0.0;
      if (iStep == 0)
        { k = 1.02; }
      else if (iStep == 1)
        {
          const double dSlope = (sM[1] - sM[0]) / (sK[1] - sK[0]);
          k = (dSlope < 0.0) ? sK[1] - sM[1] / dSlope : 2.0;
        }
      else
        {
          // the parabola through the three points, its root past the last
          const double k0 = sK[0], k1 = sK[1], k2 = sK[2], m0 = sM[0], mm1 = sM[1], m2 = sM[2];
          const double d01 = (mm1 - m0) / (k1 - k0), d12 = (m2 - mm1) / (k2 - k1), a = (d12 - d01) / (k2 - k0);
          const double b = d12 - a * (k1 + k2), c = m2 - k2 * (a * k2 + b);
          double dRoot = 0.0;
          if (smos_Fabs(a) > 1.0e-12)
            {
              const double dDisc = b * b - 4.0 * a * c;
              if (dDisc >= 0.0)
                {
                  const double r1 = (-b - smos_Sqrt(dDisc)) / (2.0 * a), r2 = (-b + smos_Sqrt(dDisc)) / (2.0 * a);
                  // the root nearest past the last scale that fails
                  const double dA = smos_Min(r1, r2), dB = smos_Max(r1, r2);
                  dRoot = (dA > dLo) ? dA : dB;
                }
            }
          if (!(dRoot > dLo))
            { dRoot = (d12 < 0.0) ? k2 - m2 / d12 : 2.0; }
          dPredict = dRoot;
          break;
        }
      k = smos_Min(smos_Max(k, dLo * 1.001), 2.0);
      double m = 0.0;
      eStatus = margin(k, m);
      if (eStatus != SM_SUCCESS) { return fail(eStatus); }
      sK.push_back(k);
      sM.push_back(m);
      if (m > 0.0) { dLo = smos_Max(dLo, k); }
      if (iStep == 1) { dPredict = k; }
    }
  // the certifying end: a half percent past the prediction (and past the last failing scale),
  // stepping up by 1e-2 while it does not certify, up to 2
  double dHi = smos_Max(dPredict * 1.005, dLo * 1.01);
  for (;;)
    {
      if (dHi > 2.0)
        {
          dHi = 2.0;
        }
      eStatus = model(dHi, NULL, bCertified);
      if (eStatus != SM_SUCCESS) { return fail(eStatus); }
      if (bCertified)
        { break; }
      dLo = dHi;
      if (dHi >= 2.0)
        {
          if (ppReason != NULL) { *ppReason = _T("SmBlendSetback: no G1 fill is certified convex at the corner's setbacks up to twice them"); }
          return SM_SUCCESS;
        }
      dHi *= 1.01;
    }
  // the failing end within 1e-2 below it
  while (dHi > dLo * 1.01)
    {
      const double k = dHi / 1.01;
      double m = 0.0;
      eStatus = margin(k, m);
      if (eStatus != SM_SUCCESS) { return fail(eStatus); }
      if (m > 0.0)
        {
          dLo = k;
          break;
        }
      eStatus = model(k, NULL, bCertified);
      if (eStatus != SM_SUCCESS) { return fail(eStatus); }
      if (!bCertified)
        {
          dLo = k;
          break;
        }
      dHi = k;
    }
  // on a grid of 1e-3, up (certifying still: the conditions only grow easier with k), so that
  // one corner described two ways (its faces and edges in another order) gets one scale
  rdScale = std::ceil(dHi * 1000.0) / 1000.0;
  return SM_SUCCESS;
}

SmStatus SmBlendSetback::Containment(const Fill & crFill, const std::vector<Face> & crFaces, SmBoolean bConvex, double & rdWorst)
{
  rdWorst = -HUGE_VAL;
  if (crFill.m_vPatches.empty())
    { SER_MSG(SM_ERR_INVALID_INPUT, _T("SmBlendSetback: an empty fill")); }
  const double dSide = bConvex ? 1.0 : -1.0;
  bool bAnyCurved = false;
  for (size_t f=0; f<crFaces.size(); f++) { bAnyCurved = bAnyCurved || crFaces[f].m_eKind != SM_SB_PLANE; }
  if (bAnyCurved) { rdWorst = 0.0; }
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
          {
            if (bAnyCurved)
              { continue; }
            rdWorst = smos_Max(rdWorst, dSide * (sPts[j] - crFaces[f].m_sPoint).Dot(crFaces[f].m_sNormal));
          }
      // a quadric face: its implicit form at the patch, Qh = |U - h a|^2 - rho_w^2 (cylinder, cone) or
      // |U|^2 - R^2 w^2 (sphere), U = H_xyz - w X0, h = U . a, rho_w = R w + slope h, of the sign of
      // the face's distance times m_dSide: certified dSide m_dSide Qh <= 0 over the patch (Bernstein,
      // subdivision). Along a column on the face (the setback patch's v = 0 or 1) Qh vanishes to
      // second order (the column on the face, the patch G1 to it): there Qh = v^2 R (or (1 - v)^2 R),
      // the division exact in Bernstein form, and R is certified instead (its dropped coefficients
      // checked to be rounding).
      const size_t n = crFill.m_vFace.size();
      for (size_t f=0; f<crFaces.size(); f++)
        {
          const Face & crF = crFaces[f];
          if (!bAnyCurved)
            { continue; }
          const ULONG p = lUDeg, q = lVDeg;
          sm_SbT sW = sm_SbTMake(p, q), sH = sm_SbTMake(p, q), sU[3];
          const double dP[3] = { crF.m_sPoint.x, crF.m_sPoint.y, crF.m_sPoint.z }, dA[3] = { crF.m_sAxis.x, crF.m_sAxis.y, crF.m_sAxis.z };
          for (ULONG k=0; k<sPts.GetSize(); k++) { sW.c[k] = (sWeights.GetSize() == sPts.GetSize()) ? sWeights[k] : 1.0; }
          for (int c=0; c<3; c++)
            {
              sU[c] = sm_SbTMake(p, q);
              for (ULONG k=0; k<sPts.GetSize(); k++)
                {
                  const double x = (c == 0) ? sPts[k].x : (c == 1) ? sPts[k].y : sPts[k].z;
                  sU[c].c[k] = sW.c[k] * (x - dP[c]);
                }
              if (crF.m_eKind != SM_SB_SPHERE)
                for (ULONG k=0; k<sPts.GetSize(); k++) { sH.c[k] += sU[c].c[k] * dA[c]; }
            }
          sm_SbT sQ, sA = sm_SbTMake(0, 0), sB = sm_SbTMake(0, 0);
          double dScale = 0.0;
          if (crF.m_eKind == SM_SB_PLANE)
            {
              // a plane beside a curved face: w (X - X0) . N, Distance's sign already
              sQ = sm_SbTMake(p, q);
              for (int c=0; c<3; c++)
                {
                  const double dN = (c == 0) ? crF.m_sNormal.x : (c == 1) ? crF.m_sNormal.y : crF.m_sNormal.z;
                  sQ = sm_SbTAxpy(sQ, dN, sU[c]);
                  dScale += smos_Fabs(dN) * sm_SbTScale(sU[c]);
                }
              for (size_t k=0; k<sQ.c.size(); k++) { sQ.c[k] *= -dSide; }
            }
          else
          {
          for (int c=0; c<3; c++)
            {
              const sm_SbT sR = (crF.m_eKind == SM_SB_SPHERE) ? sU[c] : sm_SbTAxpy(sU[c], -dA[c], sH);
              sA = sm_SbTAxpy(sA, 1.0, sm_SbTMulFast(sR, sR));
            }
          const sm_SbT sRho = (crF.m_eKind == SM_SB_SPHERE) ? sm_SbTAxpy(sm_SbTMake(p, q), crF.m_dRadius, sW)
                                                             : sm_SbTAxpy(sm_SbTAxpy(sm_SbTMake(p, q), crF.m_dRadius, sW), crF.m_dSlope, sH);
          sB = sm_SbTMulFast(sRho, sRho);
          sQ = sm_SbTAxpy(sA, -1.0, sB);
          dScale = sm_SbTScale(sA) + sm_SbTScale(sB);
          for (size_t k=0; k<sQ.c.size(); k++) { sQ.c[k] *= -dSide * crF.m_dSide; } // >= 0 where contained
          }
          // a column of this patch on this face
          int iCol = -1;
          if (i < n)
            {
              if (crFill.m_vFace[i] == f) { iCol = 0; }
              else if (crFill.m_vFace[(i + 1) % n] == f) { iCol = 1; }
            }
          if (iCol >= 0)
            {
              const size_t Q = sQ.q;
              if (Q < 2)
                { SER(SM_ERR); }
              sm_SbT sR = sm_SbTMake(sQ.p, Q - 2);
              double dDropped = 0.0;
              for (size_t a=0; a<=sQ.p; a++)
                for (size_t j=0; j<=Q; j++)
                  {
                    const double c = sQ.c[a * (Q + 1) + j];
                    const size_t jj = (iCol == 0) ? j : Q - j; // the power of v (or 1 - v)
                    if (jj < 2) { dDropped = smos_Max(dDropped, smos_Fabs(c)); continue; }
                    const double f2 = static_cast<double>(Q) * static_cast<double>(Q - 1) / (static_cast<double>(jj) * static_cast<double>(jj - 1));
                    sR.c[a * (Q - 1) + ((iCol == 0) ? j - 2 : j)] = c * f2;
                  }
              if (dDropped > 1.0e-11 * dScale)
                { rdWorst = HUGE_VAL; return SM_SUCCESS; }
              sQ = sR;
            }
          if (sm_SbTNonneg(sQ, 1.0e-12 * dScale, SM_SB_GRAPH_DEPTH) != 1)
            { rdWorst = HUGE_VAL; return SM_SUCCESS; }
        }
    }
  return SM_SUCCESS;
}

/*******************************************************************//**
PURPOSE: The spring curve on a face (SmBlendSetback.h).

NOTES: The face's rational parameter plane: a plane's own coordinates
    along two unit directions in it; a cylinder's or cone's (t, h), t =
    tan(theta / 2) about the axis from the radial direction e1 halfway
    between P's and Q's, h along the axis; a sphere's stereographic plane
    (a, b) from the point opposite the middle of P's and Q's directions. A
    tangent T at a point maps by the parametrization's differential:
    dtheta = T . e_theta / rho, dt = (1 + t^2) dtheta / 2, dh = T . axis;
    on the sphere (a, b) = (u1, u2) / (R + u3), u the point's coordinates
    from the centre along (e1, e2, m), so da = (du1 (R + u3) - u1 du3) /
    (R + u3)^2. The parabola (P^, X^, Q^), X^ where P^ + a dP^ meets
    Q^ - b dQ^ (a, b > 0), mapped up:
      cylinder: W = 1 + t^2, (X, Y) = R (1 - t^2, 2 t), the axial h W of
        degree 6 replaced by the quartic with its end values and
        derivatives (its middle coefficient their neighbours' mean): the
        curve stays on the cylinder (X^2 + Y^2 = R^2 W^2 whatever the
        axial coordinate) and keeps its end tangents;
      cone: W = 1 + t^2, (X, Y) = rho(h) (1 - t^2, 2 t), the axial h W, all
        of degree 6;
      sphere: W = 1 + a^2 + b^2, R (2 a, 2 b, 1 - a^2 - b^2) about the
        centre, degree 4.
***********************************************************************/
SmStatus SmBlendSetback::SpringCurve(const Face & crFace, const SmPoint3d & crP, const SmVector3d & crTP, const SmPoint3d & crQ, const SmVector3d & crTQ,
                                     std::vector<double> & rFromP, std::vector<double> & rFromQ, const TCHAR ** ppReason,
                                     std::vector<double> * pTanP, std::vector<double> * pTanQ)
{
  rFromP.clear();
  rFromQ.clear();
  const TCHAR * pLocalWhy = NULL;
  const TCHAR *& rWhy = (ppReason != NULL) ? *ppReason : pLocalWhy;
  rWhy = NULL;
  const SmVector3d sO(crFace.m_sPoint);
  // the frame: e3 the axis (a plane's normal), e1 radial halfway between P and Q (a plane's: P to Q)
  SmVector3d sE3 = (crFace.m_eKind == SM_SB_PLANE) ? crFace.m_sNormal : crFace.m_sAxis;
  SER( sE3.Unitize() );
  SmVector3d sE1;
  if ( crFace.m_eKind == SM_SB_PLANE )
    { sE1 = crQ - crP; }
  else if ( crFace.m_eKind == SM_SB_SPHERE )
    {
      SmVector3d sUP = SmVector3d(crP) - sO, sUQ = SmVector3d(crQ) - sO;
      SER( sUP.Unitize() );
      SER( sUQ.Unitize() );
      sE3 = sUP + sUQ; // m: the middle of P's and Q's directions
      SER( sE3.Unitize() );
      sE1 = sUQ - sUP;
    }
  else
    {
      SmVector3d sRP = SmVector3d(crP) - sO, sRQ = SmVector3d(crQ) - sO;
      sRP = sRP - sE3 * sRP.Dot( sE3 );
      sRQ = sRQ - sE3 * sRQ.Dot( sE3 );
      SER( sRP.Unitize() );
      SER( sRQ.Unitize() );
      sE1 = sRP + sRQ;
    }
  sE1 = sE1 - sE3 * sE1.Dot( sE3 );
  SER( sE1.Unitize() );
  const SmVector3d sE2 = sE3 * sE1;
  const double R = crFace.m_dRadius;

  // a point and a tangent to the parameter plane
  auto toPlane = [&](const SmPoint3d & crX, const SmVector3d & crT, double (&rU)[2], double (&rD)[2]) -> SmStatus
    {
      const SmVector3d sU = SmVector3d(crX) - sO;
      const double u1 = sU.Dot( sE1 ), u2 = sU.Dot( sE2 ), u3 = sU.Dot( sE3 );
      const double d1 = crT.Dot( sE1 ), d2 = crT.Dot( sE2 ), d3 = crT.Dot( sE3 );
      switch ( crFace.m_eKind )
        {
        case SM_SB_PLANE:
          rU[0] = u1; rU[1] = u2; rD[0] = d1; rD[1] = d2;
          return SM_SUCCESS;
        case SM_SB_SPHERE:
          {
            const double dDen = R + u3;
            if ( !( dDen > SM_EFF_ZERO * R ))
              { SER( SM_ERR ); }
            rU[0] = u1 / dDen; rU[1] = u2 / dDen;
            rD[0] = ( d1 * dDen - u1 * d3 ) / ( dDen * dDen );
            rD[1] = ( d2 * dDen - u2 * d3 ) / ( dDen * dDen );
            return SM_SUCCESS;
          }
        default:
          {
            const double dRho = sqrt( u1 * u1 + u2 * u2 );
            if ( !( dRho > SM_EFF_ZERO ))
              { SER( SM_ERR ); }
            const double dTheta = atan2( u2, u1 );
            if ( smos_Fabs( dTheta ) >= 0.5 * SM_PI )
              { SER( SM_ERR ); } // more than a half turn between P and Q
            const double t = tan( 0.5 * dTheta );
            const double dDTheta = ( -sin( dTheta ) * d1 + cos( dTheta ) * d2 ) / dRho;
            rU[0] = t; rU[1] = u3;
            rD[0] = 0.5 * ( 1.0 + t * t ) * dDTheta; rD[1] = d3;
            return SM_SUCCESS;
          }
        }
    };
  double adP[2], adDP[2], adQ[2], adDQ[2];
  SER( toPlane( crP, crTP, adP, adDP ));
  SER( toPlane( crQ, crTQ, adQ, adDQ ));
  // X^: P^ + a dP^ = Q^ - b dQ^, i.e. a dP^ + b dQ^ = Q^ - P^
  const double dDet = adDP[0] * adDQ[1] - adDP[1] * adDQ[0];
  if ( !( smos_Fabs( dDet ) > SM_EFF_ZERO ))
    { rWhy = _T("SmBlendSetback: two rails on a face run parallel (no corner region between them)"); return SM_ERR_INVALID_INPUT; }
  const double r0 = adQ[0] - adP[0], r1 = adQ[1] - adP[1];
  const double a = ( r0 * adDQ[1] - r1 * adDQ[0] ) / dDet;
  const double b = ( adDP[0] * r1 - adDP[1] * r0 ) / dDet;
  if ( !( a > 0.0 ) || !( b > 0.0 ))
    { rWhy = _T("SmBlendSetback: two rounds' rails on a face cross behind their ends (no corner region between them)"); return SM_ERR_INVALID_INPUT; }
  const sm_SbB sU = { adP[0], adP[0] + a * adDP[0], adQ[0] };
  const sm_SbB sV = { adP[1], adP[1] + a * adDP[1], adQ[1] };

  // the homogeneous curve's components about the frame (c1, c2, c3 along e1, e2, e3, then w)
  sm_SbB c1, c2, c3, w;
  switch ( crFace.m_eKind )
    {
    case SM_SB_PLANE:
      c1 = sU; c2 = sV; c3 = sm_SbB( 3, 0.0 ); w = sm_SbB( 3, 1.0 );
      break;
    case SM_SB_SPHERE:
      {
        const sm_SbB sUU = sm_SbMul( sU, sU ), sVV = sm_SbMul( sV, sV );
        w.resize( 5 ); c1.resize( 5 ); c2.resize( 5 ); c3.resize( 5 );
        const sm_SbB sU4 = sm_SbElevate( sU, 4 ), sV4 = sm_SbElevate( sV, 4 );
        const sm_SbB sUU4 = sm_SbElevate( sUU, 4 ), sVV4 = sm_SbElevate( sVV, 4 );
        for (size_t i=0; i<5; i++)
          {
            w[i]  = 1.0 + sUU4[i] + sVV4[i];
            c1[i] = R * 2.0 * sU4[i];
            c2[i] = R * 2.0 * sV4[i];
            c3[i] = R * ( 1.0 - sUU4[i] - sVV4[i] );
          }
        break;
      }
    default:
      {
        // t quadratic: t^2, W = 1 + t^2 and 1 - t^2 quartic, 2 t raised to 4
        const sm_SbB sTT = sm_SbMul( sU, sU );
        sm_SbB sW( 5 ), sOneMinus( 5 );
        for (size_t i=0; i<5; i++) { sW[i] = 1.0 + sTT[i]; sOneMinus[i] = 1.0 - sTT[i]; }
        const sm_SbB sTwoT = sm_SbElevate( sm_SbB( { 2.0 * sU[0], 2.0 * sU[1], 2.0 * sU[2] } ), 4 );
        if ( crFace.m_eKind == SM_SB_CYLINDER )
          {
            w = sW;
            c1 = sOneMinus;
            c2 = sTwoT;
            for (size_t i=0; i<5; i++) { c1[i] *= R; c2[i] *= R; }
            // the axial coordinate: the quartic with h W's end values and derivatives (h W is of
            // degree 6), its middle coefficient its neighbours' mean
            const sm_SbB sHW = sm_SbMul( sV, sW );
            c3.assign( 5, 0.0 );
            c3[0] = sHW[0];
            c3[4] = sHW[6];
            c3[1] = sHW[0] + ( sHW[1] - sHW[0] ) * 6.0 / 4.0;
            c3[3] = sHW[6] - ( sHW[6] - sHW[5] ) * 6.0 / 4.0;
            c3[2] = 0.5 * ( c3[1] + c3[3] );
          }
        else
          {
            const sm_SbB sRho = { R + crFace.m_dSlope * sV[0], R + crFace.m_dSlope * sV[1], R + crFace.m_dSlope * sV[2] };
            w  = sm_SbElevate( sW, 6 );
            c1 = sm_SbMul( sRho, sOneMinus );
            c2 = sm_SbMul( sRho, sTwoT );
            c3 = sm_SbMul( sV, sW );
          }
        break;
      }
    }
  const size_t d = w.size() - 1;
  // homogeneous world coordinates
  std::vector<sm_SbH> sCurve( d + 1 );
  for (size_t i=0; i<=d; i++)
    {
      const SmVector3d sX = sO * w[i] + sE1 * c1[i] + sE2 * c2[i] + sE3 * c3[i];
      if ( !( w[i] > 0.0 ))
        { rWhy = _T("SmBlendSetback: a spring curve's weight would not stay positive"); return SM_ERR_INVALID_INPUT; }
      sCurve[i].a[0] = sX.x; sCurve[i].a[1] = sX.y; sCurve[i].a[2] = sX.z; sCurve[i].a[3] = w[i];
    }
  // the halves at s = 1/2 (de Casteljau), the second reversed (from Q), each times L so that its
  // weights at the base point are 1, 1: L(1) = 1 / W(1), L(0) = L(1) (1 + W'(1) / W(1))
  std::vector<sm_SbH> sLeft( d + 1 ), sRight( d + 1 ), sWork( sCurve );
  for (size_t k=0; k<=d; k++)
    {
      sLeft[k]  = sWork[0];
      sRight[k] = sWork[d - k]; // from Q (k = 0) to the middle (k = d)
      for (size_t i=0; i+k<d; i++)
        { sWork[i] = sm_SbAxpy( sm_SbAxpy( sm_SbZero(), 0.5, sWork[i] ), 0.5, sWork[i + 1] ); }
    }
  auto out = [&](const std::vector<sm_SbH> & crHalf, std::vector<double> & rOut)
    {
      rOut.assign( 4 * ( d + 1 ), 0.0 );
      for (size_t i=0; i<=d; i++)
        for (int c=0; c<4; c++) { rOut[4 * i + c] = crHalf[i].a[c]; }
    };
  out( sLeft, rFromP );
  out( sRight, rFromQ );

  // the tangent directions along each half: the parameter-plane parabola split at its middle
  // (from P, and from Q), the face's partials there
  auto tangents = [&](const sm_SbB & crU, const sm_SbB & crV, std::vector<double> & rOut)
    {
      std::vector< std::vector<SmVector3d> > sXYZ( 2 );
      std::vector< std::vector<double> > sW( 2 );
      switch ( crFace.m_eKind )
        {
        case SM_SB_PLANE:
          sXYZ[0] = { sE1 }; sW[0] = { 0.0 };
          sXYZ[1] = { sE2 }; sW[1] = { 0.0 };
          break;
        case SM_SB_SPHERE:
          // S(a, b) = O (1 + a^2 + b^2) + R (2 a e1 + 2 b e2 + (1 - a^2 - b^2) m), w = 1 + a^2 + b^2:
          // S_a = (2 a O + R (2 e1 - 2 a m), 2 a), S_b likewise; affine in (a, b), quadratic along it
          for (int k=0; k<2; k++)
            for (size_t i=0; i<3; i++)
              {
                const double dA = ( k == 0 ) ? crU[i] : crV[i];
                const SmVector3d & crE = ( k == 0 ) ? sE1 : sE2;
                sXYZ[k].push_back( sO * ( 2.0 * dA ) + ( crE * 2.0 - sE3 * ( 2.0 * dA )) * R );
                sW[k].push_back( 2.0 * dA );
              }
          break;
        default:
          {
            // the circumferential direction -2 t e1 + (1 - t^2) e2 (degree 4), and the cylinder's axis or
            // the cone's generator kappa (1 - t^2) e1 + 2 kappa t e2 + (1 + t^2) e3 (degree 4)
            const sm_SbB sTT = sm_SbMul( crU, crU );
            const sm_SbB sT4 = sm_SbElevate( crU, 4 );
            for (size_t i=0; i<5; i++)
              {
                sXYZ[0].push_back( sE1 * ( -2.0 * sT4[i] ) + sE2 * ( 1.0 - sTT[i] ));
                sW[0].push_back( 0.0 );
              }
            if ( crFace.m_eKind == SM_SB_CYLINDER )
              { sXYZ[1] = { sE3 }; sW[1] = { 0.0 }; }
            else
              for (size_t i=0; i<5; i++)
                {
                  const double k = crFace.m_dSlope;
                  sXYZ[1].push_back( sE1 * ( k * ( 1.0 - sTT[i] )) + sE2 * ( 2.0 * k * sT4[i] ) + sE3 * ( 1.0 + sTT[i] ));
                  sW[1].push_back( 0.0 );
                }
            break;
          }
        }
      rOut.clear();
      for (int k=0; k<2; k++)
        for (size_t i=0; i<sXYZ[k].size(); i++)
          { rOut.insert( rOut.end(), { sXYZ[k][i].x, sXYZ[k][i].y, sXYZ[k][i].z, sW[k][i] } ); }
    };
  const sm_SbB sUL = { sU[0], 0.5 * ( sU[0] + sU[1] ), 0.25 * ( sU[0] + 2.0 * sU[1] + sU[2] ) };
  const sm_SbB sVL = { sV[0], 0.5 * ( sV[0] + sV[1] ), 0.25 * ( sV[0] + 2.0 * sV[1] + sV[2] ) };
  const sm_SbB sUR = { sU[2], 0.5 * ( sU[1] + sU[2] ), sUL[2] };
  const sm_SbB sVR = { sV[2], 0.5 * ( sV[1] + sV[2] ), sVL[2] };
  if ( pTanP != NULL ) { tangents( sUL, sVL, *pTanP ); }
  if ( pTanQ != NULL ) { tangents( sUR, sVR, *pTanQ ); }
  return SM_SUCCESS;
}

/*******************************************************************//**
PURPOSE: A face's signed distance and its normal out of the material
    (SmBlendSetback.h).
***********************************************************************/
SmStatus SmBlendSetback::RoundEndCurvature(const SmContext & crContext, const Round & crRound, const SmVector3d & crInto, double & rdK,
                                           const TCHAR ** ppReason)
{
  rdK = 0.0;
  std::vector< std::vector<sm_SbH> > sJ;
  size_t lQ = 0;
  double dDtDs = 1.0;
  const TCHAR * pLocal = NULL;
  const TCHAR *& rWhy = (ppReason != NULL) ? *ppReason : pLocal;
  rWhy = NULL;
  SER(sm_SbRoundJets(crContext, crRound, 2, sJ, lQ, dDtDs, rWhy));
  const bool bPlaneEnd[2] = { false, false };
  return sm_SbRoundRowK(sJ, crInto, bPlaneEnd, false, rdK, rWhy);
}

SmStatus SmBlendSetback::NoMoreConcave(const SmBSplineSurface & crPatch, const Face & crA, const Face & crB, const SmVector3d & crInto, double dK,
                                       int & riResult)
{
  riResult = -1;
  ULONG lUDeg = 0, lVDeg = 0;
  SmTArray<SmPoint3d> sPts;
  SmTArray<ULONG> sUMult, sVMult;
  SmTArray<double> sUKnots, sVKnots, sWeights;
  SmBSplineSurfaceForm eForm;
  SmKnotType eKnotType;
  SER(crPatch.GetCanonical(lUDeg, lVDeg, sPts, eForm, sUMult, sVMult, sUKnots, sVKnots, eKnotType, sWeights));
  if (sUKnots.GetSize() != 2 || sVKnots.GetSize() != 2 || !(dK < 0.0) || sPts.GetSize() != (lUDeg + 1) * (lVDeg + 1))
    { SER(SM_ERR_INVALID_INPUT); }
  sm_SbNet sN;
  sN.p = lUDeg;
  sN.q = lVDeg;
  sN.P.resize(sPts.GetSize());
  for (ULONG i=0; i<sPts.GetSize(); i++)
    {
      const double w = (sWeights.GetSize() == sPts.GetSize()) ? sWeights[i] : 1.0;
      sN.P[i].a[0] = sPts[i].x * w;
      sN.P[i].a[1] = sPts[i].y * w;
      sN.P[i].a[2] = sPts[i].z * w;
      sN.P[i].a[3] = w;
    }
  riResult = sm_SbNoMoreConcave(sN, crA, crB, crInto, dK);
  return SM_SUCCESS;
}

double SmBlendSetback::Distance(const Face & crF, const SmPoint3d & crX, SmVector3d * pNormal)
{
  const SmVector3d sU = crX - crF.m_sPoint;
  SmVector3d sN( 0.0, 0.0, 0.0 );
  double dG = 0.0;
  switch ( crF.m_eKind )
    {
    case SM_SB_PLANE:
      sN = crF.m_sNormal;
      dG = sU.Dot( sN );
      break;
    case SM_SB_SPHERE:
      {
        const double dLen = sU.Length();
        sN = ( dLen > 0.0 ) ? sU * ( crF.m_dSide / dLen ) : sN;
        dG = crF.m_dSide * ( dLen - crF.m_dRadius );
        break;
      }
    default:
      {
        const double dH = sU.Dot( crF.m_sAxis );
        const SmVector3d sR = sU - crF.m_sAxis * dH;
        const double dRho = sR.Length();
        const double dCos = 1.0 / sqrt( 1.0 + crF.m_dSlope * crF.m_dSlope ), dSin = crF.m_dSlope * dCos;
        if ( dRho > 0.0 )
          { sN = ( sR * ( dCos / dRho ) - crF.m_sAxis * dSin ) * crF.m_dSide; }
        dG = crF.m_dSide * ( dRho - ( crF.m_dRadius + crF.m_dSlope * dH )) * dCos;
        break;
      }
    }
  if ( pNormal != NULL ) { *pNormal = sN; }
  return dG;
}
