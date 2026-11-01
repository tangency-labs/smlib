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

#include <algorithm>
#include <cmath>
#include <vector>

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
                                   std::vector<double> * pvLeast)
{
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
        { SER_MSG(SM_ERR_INVALID_INPUT, _T("SmBlendSetback: an edge needs a positive radius and two of the faces given")); }
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
        { SER_MSG(SM_ERR_INVALID_INPUT, _T("SmBlendSetback: a face at the corner must lie between exactly two of its rounded edges")); }
      const size_t i = sBeside[0], j = sBeside[1];
      const double dCos = smos_Max(-1.0, smos_Min(1.0, crEdges[i].m_sAway.Dot(crEdges[j].m_sAway)));
      const double dSin = smos_Sqrt(smos_Max(0.0, 1.0 - dCos * dCos));
      if (!(dSin > 1.0e-9))
        { SER_MSG(SM_ERR_INVALID_INPUT, _T("SmBlendSetback: two edges on a face run along one line")); }
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

SmStatus SmBlendSetback::FillHole(const SmContext & crContext, const std::vector<Face> & crFaces, const std::vector<Round> & crRounds, ULONG lOrder,
                                  SmBoolean bConvex, double dTol, Fill & rFill)
{
  (void)crContext; (void)crFaces; (void)crRounds; (void)lOrder; (void)bConvex; (void)dTol; (void)rFill;
  SER_MSG(SM_ERR, _T("SmBlendSetback: the setback fill is not built yet"));
  return SM_ERR;
}

SmStatus SmBlendSetback::Containment(const Fill & crFill, const std::vector<Face> & crFaces, SmBoolean bConvex, double & rdWorst)
{
  (void)crFill; (void)crFaces; (void)bConvex;
  rdWorst = HUGE_VAL;
  SER_MSG(SM_ERR, _T("SmBlendSetback: the setback fill is not built yet"));
  return SM_ERR;
}
