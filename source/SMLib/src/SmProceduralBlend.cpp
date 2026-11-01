// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

/**********************************************************************//**
* FILE NAME --- SmProceduralBlend.cpp
* PURPOSE: SolidCraft change. A blend surface stored as its definition
*   (SmProceduralBlend.h): its ownership, parametrization, the operations
*   that keep the net and the definition alike, and its persistence. Its
*   evaluation is in SmFilletSection.cpp, beside the G^n section's.
**********************************************************************/

#include "StdAfx.h"

#include <SmProceduralBlend.h>
#include <SmBlendNetScope.h>
#include <SmBSplineCurve.h>
#include <SmDatabaseIO.h>
#include <SmFilletSolver.h>
#include <SmPseudoBox.h>

#include <atomic>
#include <cmath>
#include <vector>

namespace
{
std::atomic<unsigned long long> s_lEvaluations(0), s_lSolves(0), s_lNewtonSteps(0), s_lMostNewtonSteps(0);
std::atomic<unsigned long long> s_lLastSerial(0);

SmBoolean sm_Finite(double d) { return std::isfinite(d) ? TRUE : FALSE; }

// a copy of a curve, as a B-spline that evaluates past its ends on its end spans
SmStatus sm_CopyCurve(const SmContext & crContext, const SmBSplineCurve & crCurve, SmBSplineCurve *& rpCopy)
{
  rpCopy = NULL;
  SmCurve * pCopy = NULL;
  SER(crCurve.Copy(crContext, pCopy));
  rpCopy = SM_CAST_PTR(SmBSplineCurve, pCopy);
  if (rpCopy == NULL)
    {
      delete pCopy;
      SER(SM_ERR);
    }
  rpCopy->SetOutOfBoundsEnabled(TRUE);
  return SM_SUCCESS;
}

SmStatus sm_WriteCurve(SmDatabaseIO & rDB, const SmBSplineCurve & crCurve, ULONG lDBVersionNumber)
{
  ULONG lDim = crCurve.GetDim();
  SER(rDB.WriteType(crCurve.GetType(), &lDim));
  SER(crCurve.WriteToDB(rDB, lDBVersionNumber));
  return SM_SUCCESS;
}

SmStatus sm_ReadCurve(SmDatabaseIO & rDB, const SmContext & crContext, ULONG lDBVersionNumber, SmBSplineCurve *& rpCurve)
{
  rpCurve = NULL;
  SM_TYPE lType;
  ULONG lDim = 0;
  SER(rDB.ReadType(lType, &lDim));
  SmCurve * pCurve = NULL;
  SER(SmCurve::ReadFromDB(lType, rDB, lDim, crContext, pCurve, lDBVersionNumber));
  rpCurve = SM_CAST_PTR(SmBSplineCurve, pCurve);
  if (rpCurve == NULL)
    {
      delete pCurve;
      SER_MSG(SM_ERR, _T("SmProceduralBlend: a curve of its definition is not a B-spline"));
    }
  rpCurve->SetOutOfBoundsEnabled(TRUE);
  return SM_SUCCESS;
}

SmBoolean sm_SameCurve(const SmBSplineCurve * a, const SmBSplineCurve * b)
{
  return (a == b || (a != NULL && b != NULL && *a == *b)) ? TRUE : FALSE;
}
} // namespace

namespace
{
thread_local int t_lNetScopes = 0;
}

SmBlendNetScope::SmBlendNetScope()
{
  t_lNetScopes++;
}

SmBlendNetScope::~SmBlendNetScope()
{
  t_lNetScopes--;
}

SmBoolean SmBlendNetScope::IsActive()
{
  return (t_lNetScopes > 0) ? TRUE : FALSE;
}

SmBlendCounts SmProceduralBlend::GetCounts()
{
  SmBlendCounts sCounts;
  sCounts.lEvaluations = s_lEvaluations.load();
  sCounts.lSolves = s_lSolves.load();
  sCounts.lNewtonSteps = s_lNewtonSteps.load();
  sCounts.lMostNewtonSteps = s_lMostNewtonSteps.load();
  return sCounts;
}

void SmProceduralBlend::ResetCounts()
{
  s_lEvaluations = 0;
  s_lSolves = 0;
  s_lNewtonSteps = 0;
  s_lMostNewtonSteps = 0;
}

void SmProceduralBlend::CountEvaluation()
{
  s_lEvaluations++;
}

void SmProceduralBlend::Renew()
{
  m_lSerial = ++s_lLastSerial;
}

void SmProceduralBlend::CountSolve(ULONG lNewtonSteps)
{
  s_lSolves++;
  s_lNewtonSteps += lNewtonSteps;
  unsigned long long lMost = s_lMostNewtonSteps.load();
  while (lNewtonSteps > lMost && !s_lMostNewtonSteps.compare_exchange_weak(lMost, lNewtonSteps)) {}
}

SmProceduralBlend::SmProceduralBlend()
  : m_pLaw(NULL), m_lOrder(0), m_dThumbweight(0.0), m_pLeader(NULL), m_dBound(0.0), m_bSwapped(FALSE)
{
  m_apSupport[0] = m_apSupport[1] = NULL;
  m_apSeed[0] = m_apSeed[1] = NULL;
  m_adSide[0] = m_adSide[1] = 1.0;
  m_adR[0] = m_adR[1] = 0.0;
  m_adA[0] = m_adA[1] = 0.0;
  m_adB[0] = m_adB[1] = 1.0;
  Renew();
}

SmProceduralBlend::SmProceduralBlend(const SmBSplineSurface & crNet)
  : SmBSplineSurface(crNet), m_pLaw(NULL), m_lOrder(0), m_dThumbweight(0.0), m_pLeader(NULL), m_dBound(0.0),
    m_bSwapped(FALSE)
{
  m_apSupport[0] = m_apSupport[1] = NULL;
  m_apSeed[0] = m_apSeed[1] = NULL;
  m_adSide[0] = m_adSide[1] = 1.0;
  m_adR[0] = m_adR[1] = 0.0;
  m_adA[0] = m_adA[1] = 0.0;
  m_adB[0] = m_adB[1] = 1.0;
  Renew();
}

SmProceduralBlend::SmProceduralBlend(const SmProceduralBlend & crSource)
  : SmBSplineSurface(crSource), m_pLaw(NULL), m_lOrder(0), m_dThumbweight(0.0), m_pLeader(NULL), m_dBound(0.0),
    m_bSwapped(FALSE)
{
  m_apSupport[0] = m_apSupport[1] = NULL;
  m_apSeed[0] = m_apSeed[1] = NULL;
  m_adSide[0] = m_adSide[1] = 1.0;
  m_adR[0] = m_adR[1] = 0.0;
  m_adA[0] = m_adA[1] = 0.0;
  m_adB[0] = m_adB[1] = 1.0;
  Renew();
  // a constructor cannot report: a failed copy leaves no definition, which Copy refuses
  (void)CopyDefinition(crSource);
}

SmProceduralBlend::~SmProceduralBlend()
{
  Clear();
}

void SmProceduralBlend::Clear()
{
  for (int i=0; i<2; i++)
    {
      delete m_apSupport[i];
      m_apSupport[i] = NULL;
      delete m_apSeed[i];
      m_apSeed[i] = NULL;
    }
  delete m_pLaw;
  m_pLaw = NULL;
  delete m_pLeader;
  m_pLeader = NULL;
}

SmStatus SmProceduralBlend::CopyDefinition(const SmProceduralBlend & crSource)
{
  Clear();
  const SmContext * pContext = GetContext();
  NER(pContext);
  for (int i=0; i<2; i++)
    {
      NER(crSource.m_apSupport[i]);
      NER(crSource.m_apSeed[i]);
      SER(crSource.m_apSupport[i]->Copy(*pContext, m_apSupport[i]));
      SER(sm_CopyCurve(*pContext, *crSource.m_apSeed[i], m_apSeed[i]));
      m_adSide[i] = crSource.m_adSide[i];
      m_adR[i] = crSource.m_adR[i];
      m_adA[i] = crSource.m_adA[i];
      m_adB[i] = crSource.m_adB[i];
    }
  NER(crSource.m_pLeader);
  SER(sm_CopyCurve(*pContext, *crSource.m_pLeader, m_pLeader));
  if (crSource.m_pLaw != NULL)
    { SER(sm_CopyCurve(*pContext, *crSource.m_pLaw, m_pLaw)); }
  m_lOrder = crSource.m_lOrder;
  m_dThumbweight = crSource.m_dThumbweight;
  m_dBound = crSource.m_dBound;
  m_bSwapped = crSource.m_bSwapped;
  Renew();
  return SM_SUCCESS;
}

/*******************************************************************//**
PURPOSE: A blend by its definition on the net crNet; every input copied.

NOTES: The refusals are the hostile records of a definition, made at
    construction (and so at load) and never at evaluation.
***********************************************************************/
SmStatus SmProceduralBlend::Create
  (const SmContext        & crContext,
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
   SmProceduralBlend     *& rpBlend)
{
  rpBlend = NULL;
  NER(apSupport[0]); NER(apSupport[1]); NER(apSeed[0]); NER(apSeed[1]);
  if (lOrder < 1 || lOrder > SM_GN_MAX_ORDER)
    { SER_MSG(SM_ERR_INVALID_INPUT, _T("SmProceduralBlend: a section's order is 1 to 4")); }
  if (lOrder >= 2 && (cpLaw != NULL || !SmGnSectionFSG::IsValid(lOrder, dThumbweight)))
    { SER_MSG(SM_ERR_INVALID_INPUT, _T("SmProceduralBlend: a G^n section needs a ball of one radius and a thumbweight below 1 or the equal area rule")); }
  for (int i=0; i<2; i++)
    if (crSide[i] != 1.0 && crSide[i] != -1.0)
      { SER_MSG(SM_ERR_INVALID_INPUT, _T("SmProceduralBlend: a support's side is +1 or -1")); }
  if (cpLaw == NULL)
    {
      for (int i=0; i<2; i++)
        if (!(sm_Finite(crOffset[i]) && crOffset[i] > 0.0))
          { SER_MSG(SM_ERR_INVALID_INPUT, _T("SmProceduralBlend: the ball's radius is finite and positive")); }
      if (lOrder == 1 && crOffset[0] != crOffset[1])
        { SER_MSG(SM_ERR_INVALID_INPUT, _T("SmProceduralBlend: a circle's ball has one radius")); }
    }
  if (!(sm_Finite(dBound) && dBound > 0.0))
    { SER_MSG(SM_ERR_INVALID_INPUT, _T("SmProceduralBlend: the net's bound is finite and positive")); }
  ULONG lDepth = 1;
  for (int i=0; i<2; i++)
    {
      const SmProceduralBlend * pInner = SM_CAST_PTR(SmProceduralBlend, apSupport[i]);
      if (pInner != NULL)
        { lDepth = smos_Max(lDepth, pInner->GetDepth() + 1); }
    }
  if (lDepth > SM_BLEND_MAX_DEPTH)
    { SER_MSG(SM_ERR_INVALID_INPUT, _T("SmProceduralBlend: blends nest at most 3 definitions deep")); }

  SmProceduralBlend * pBlend = new (crContext) SmProceduralBlend(crNet);
  NER(pBlend);
  SmObjDelete sDelete(pBlend);
  for (int i=0; i<2; i++)
    {
      SER(apSupport[i]->Copy(crContext, pBlend->m_apSupport[i]));
      SER(sm_CopyCurve(crContext, *apSeed[i], pBlend->m_apSeed[i]));
      pBlend->m_adSide[i] = crSide[i];
    }
  SER(sm_CopyCurve(crContext, crLeader, pBlend->m_pLeader));
  if (cpLaw != NULL)
    { SER(sm_CopyCurve(crContext, *cpLaw, pBlend->m_pLaw)); }
  for (int i=0; i<2; i++) { pBlend->m_adR[i] = (cpLaw == NULL) ? crOffset[i] : 0.0; }
  pBlend->m_lOrder = lOrder;
  pBlend->m_dThumbweight = (lOrder >= 2) ? dThumbweight : 0.0;
  pBlend->m_dBound = dBound;
  pBlend->Renew();
  sDelete.Clear();
  rpBlend = pBlend;
  return SM_SUCCESS;
} // end SmProceduralBlend::Create

ULONG SmProceduralBlend::GetDepth() const
{
  ULONG lDepth = 1;
  for (int i=0; i<2; i++)
    {
      const SmProceduralBlend * pInner = SM_CAST_PTR(SmProceduralBlend, m_apSupport[i]);
      if (pInner != NULL)
        { lDepth = smos_Max(lDepth, pInner->GetDepth() + 1); }
    }
  return lDepth;
}

SmPoint2d SmProceduralBlend::DefinitionParameters(const SmPoint2d & crUV) const
{
  const double p = m_bSwapped ? crUV.y : crUV.x, q = m_bSwapped ? crUV.x : crUV.y;
  return SmPoint2d(m_adA[0] + m_adB[0] * p, m_adA[1] + m_adB[1] * q);
}

void SmProceduralBlend::ReparametrizeMap(SmBoolean bU, double dAlpha, double dBeta)
{
  // u drives s unless swapped
  const int k = ((bU ? 1 : 0) != (m_bSwapped ? 1 : 0)) ? 0 : 1;
  m_adA[k] += m_adB[k] * dAlpha;
  m_adB[k] *= dBeta;
}

SmBoolean SmProceduralBlend::operator==(const SmSurface & crOther) const
{
  if (this == &crOther)
    { return TRUE; }
  if (!SmBSplineSurface::operator==(crOther))
    { return FALSE; }
  const SmProceduralBlend * pOther = SM_CAST_PTR(SmProceduralBlend, &crOther);
  if (pOther == NULL)
    { return FALSE; }
  for (int i=0; i<2; i++)
    {
      if (m_apSupport[i] == NULL || pOther->m_apSupport[i] == NULL || !(*m_apSupport[i] == *pOther->m_apSupport[i]))
        { return FALSE; }
      if (!sm_SameCurve(m_apSeed[i], pOther->m_apSeed[i]) || m_adSide[i] != pOther->m_adSide[i] || m_adR[i] != pOther->m_adR[i] || m_adA[i] != pOther->m_adA[i]
          || m_adB[i] != pOther->m_adB[i])
        { return FALSE; }
    }
  return (sm_SameCurve(m_pLeader, pOther->m_pLeader) && sm_SameCurve(m_pLaw, pOther->m_pLaw)
          && m_lOrder == pOther->m_lOrder && m_dThumbweight == pOther->m_dThumbweight && m_dBound == pOther->m_dBound
          && m_bSwapped == pOther->m_bSwapped) ? TRUE : FALSE;
}

SmStatus SmProceduralBlend::Copy(const SmContext & crContext, SmSurface *& rpNewSurface) const
{
  rpNewSurface = NULL;
  SmProceduralBlend * pCopy = new (crContext) SmProceduralBlend(*this);
  NER(pCopy);
  if (pCopy->m_pLeader == NULL)
    {
      delete pCopy;
      SER_MSG(SM_ERR, _T("SmProceduralBlend: its definition could not be copied"));
    }
  rpNewSurface = pCopy;
  return SM_SUCCESS;
}

/*******************************************************************//**
PURPOSE: The definition's partials at (u, v), in the surface's parameters:
    aDerivatives[i (lHighestVDeriv + 1) + j] = d^(i+j) S / du^i dv^j, as
    SmBSplineSurface::Evaluate lays them out. (u, v) is clamped to the
    domain as the net's evaluation clamps it.

NOTES: bNonZeroTangents: where a first partial vanishes (a degenerate
    section), its tolerance sized stand-in is the net's, whose direction
    is the sampled one SmBSplineSurface::Evaluate gives.
***********************************************************************/
SmStatus SmProceduralBlend::Evaluate
  (const SmPoint2d & crUV,
   ULONG             lHighestUDeriv,
   ULONG             lHighestVDeriv,
   SmBoolean         bUFromLeft,
   SmBoolean         bVFromLeft,
   SmBoolean         bOnlyUpperHalf,
   SmVector3d      * aDerivatives,
   SmBoolean         bNonZeroTangents,
   SmBoolean         bDoZeroSampling)
  const
{
  NER(aDerivatives);
  // meshing and display (SmBlendNetScope) take the net, the approximation made for them
  if (SmBlendNetScope::IsActive())
    {
      return SmBSplineSurface::Evaluate(crUV, lHighestUDeriv, lHighestVDeriv, bUFromLeft, bVFromLeft, bOnlyUpperHalf, aDerivatives,
                                        bNonZeroTangents, bDoZeroSampling);
    }
  if (bOnlyUpperHalf && lHighestUDeriv != lHighestVDeriv)
    { SER(SM_ERR_INVALID_INPUT); }
  SmPoint2d sUV(crUV);
  const SmExtent2d sDomain = GetNaturalUVDomain();
  if (!sDomain.ContainsPoint2d(sUV) && !IsOutOfBoundsEnabled())
    { sUV = sDomain.ClampPoint2d(sUV); }
  SER(EvaluateSimple(sUV, lHighestUDeriv, lHighestVDeriv, bUFromLeft, bVFromLeft, FALSE, aDerivatives));
  if (bNonZeroTangents && (lHighestUDeriv >= 1 || lHighestVDeriv >= 1))
    {
      const double dZero = SM_EFF_ZERO * (1.0 + aDerivatives[0].GetMaxDimension());
      const SmBoolean bU = (lHighestUDeriv >= 1 && aDerivatives[lHighestVDeriv + 1].Length() < dZero) ? TRUE : FALSE;
      const SmBoolean bV = (lHighestVDeriv >= 1 && aDerivatives[1].Length() < dZero) ? TRUE : FALSE;
      if (bU || bV)
        {
          SmVector3d asNet[4];
          SER(SmBSplineSurface::Evaluate(sUV, 1, 1, bUFromLeft, bVFromLeft, FALSE, asNet, TRUE, bDoZeroSampling));
          if (bU) { aDerivatives[lHighestVDeriv + 1] = asNet[2]; }
          if (bV) { aDerivatives[1] = asNet[1]; }
        }
    }
  return SM_SUCCESS;
} // end SmProceduralBlend::Evaluate

SmStatus SmProceduralBlend::EvaluateSimple
  (const SmPoint2d & crUV,
   ULONG             lHighestUDeriv,
   ULONG             lHighestVDeriv,
   SmBoolean         bUFromLeft,
   SmBoolean         bVFromLeft,
   SmBoolean         bOnlyUpperHalf,
   SmVector3d      * aDerivatives)
  const
{
  NER(aDerivatives);
  if (SmBlendNetScope::IsActive())
    { return SmBSplineSurface::EvaluateSimple(crUV, lHighestUDeriv, lHighestVDeriv, bUFromLeft, bVFromLeft, bOnlyUpperHalf, aDerivatives); }
  if (bOnlyUpperHalf && lHighestUDeriv != lHighestVDeriv)
    { SER(SM_ERR_INVALID_INPUT); }
  SmPoint2d sUV(crUV);
  const SmExtent2d sDomain = GetNaturalUVDomain();
  if (!sDomain.ContainsPoint2d(sUV))
    { sUV = sDomain.ClampPoint2d(sUV); }
  const SmPoint2d sST = DefinitionParameters(sUV);
  const ULONG lS = m_bSwapped ? lHighestVDeriv : lHighestUDeriv, lT = m_bSwapped ? lHighestUDeriv : lHighestVDeriv;
  std::vector<SmVector3d> sD(static_cast<size_t>((lS + 1) * (lT + 1)));
  SER(EvaluateDefinition(sST.x, sST.y, lS, lT, sD.data()));
  for (ULONG i=0; i<=lHighestUDeriv; i++)
    for (ULONG j=0; j<=lHighestVDeriv; j++)
      {
        // d/du = b_s d/ds (b_t d/dt when swapped), d/dv the other
        const ULONG a = m_bSwapped ? j : i, b = m_bSwapped ? i : j;
        const double dScale = smos_Pow(m_adB[0], static_cast<double>(a)) * smos_Pow(m_adB[1], static_cast<double>(b));
        aDerivatives[i * (lHighestVDeriv + 1) + j] = sD[static_cast<size_t>(a * (lT + 1) + b)] * dScale;
      }
  return SM_SUCCESS;
} // end SmProceduralBlend::EvaluateSimple

SmStatus SmProceduralBlend::EvaluatePoint(const SmPoint2d & crUV, SmPoint3d & rPoint) const
{
  if (SmBlendNetScope::IsActive())
    { return SmBSplineSurface::EvaluatePoint(crUV, rPoint); }
  const SmPoint2d sST = DefinitionParameters(crUV);
  SmVector3d sD;
  SER(EvaluateDefinition(sST.x, sST.y, 0, 0, &sD));
  rPoint = SmPoint3d(sD);
  return SM_SUCCESS;
}

/*******************************************************************//**
PURPOSE: The net's boxes grown by the bound: the definition lies within
    it of the net, so the boxes hold the definition's points.
***********************************************************************/
SmStatus SmProceduralBlend::CalculateBoundingBox
  (const SmExtent2d  & crUVDomain,
   SmExtent3d        * pNormalBox,
   SmPseudoBox       * pPseudoBox,
   SmPolarBox        * pPolarBox,
   const SmPseudoBox * pOptPseudoBasisGuess,
   const SmPolarBox  * pOptPolarBasisGuess,
   SmBoolean           bExpandPosBoxesByZoneTol3d)
  const
{
  SER(SmBSplineSurface::CalculateBoundingBox(crUVDomain, pNormalBox, pPseudoBox, pPolarBox, pOptPseudoBasisGuess, pOptPolarBasisGuess,
                                             bExpandPosBoxesByZoneTol3d));
  if (pNormalBox != NULL)
    { pNormalBox->ExpandAbsolute(m_dBound); }
  if (pPseudoBox != NULL)
    { pPseudoBox->ExpandAbsolute(m_dBound); }
  return SM_SUCCESS;
}

SmStatus SmProceduralBlend::Reparameterize(const SmExtent2d & crNewDomain)
{
  const SmExtent2d sOld = GetNaturalUVDomain();
  SER(SmBSplineSurface::Reparameterize(crNewDomain));
  for (int c=0; c<2; c++)
    {
      const double o0 = (c == 0) ? sOld.GetMin().x : sOld.GetMin().y, o1 = (c == 0) ? sOld.GetMax().x : sOld.GetMax().y;
      const double n0 = (c == 0) ? crNewDomain.GetMin().x : crNewDomain.GetMin().y, n1 = (c == 0) ? crNewDomain.GetMax().x : crNewDomain.GetMax().y;
      if (!(n1 != n0))
        { SER(SM_ERR_INVALID_INPUT); }
      // old = o0 + (new - n0) (o1 - o0) / (n1 - n0)
      const double dBeta = (o1 - o0) / (n1 - n0);
      ReparametrizeMap(c == 0 ? TRUE : FALSE, o0 - n0 * dBeta, dBeta);
    }
  return SM_SUCCESS;
}

SmStatus SmProceduralBlend::Reverse(SmSurfParamType eSurfParam)
{
  const SmExtent2d sDomain = GetNaturalUVDomain();
  SER(SmBSplineSurface::Reverse(eSurfParam));
  // old = min + max - new
  const SmBoolean bU = (eSurfParam == SM_SP_U) ? TRUE : FALSE;
  const double dSum = bU ? sDomain.GetMin().x + sDomain.GetMax().x : sDomain.GetMin().y + sDomain.GetMax().y;
  ReparametrizeMap(bU, dSum, -1.0);
  return SM_SUCCESS;
}

// SmBSplineSurface::SwapUV calls this once per swap of the net's parameters
void SmProceduralBlend::ToggleSwapUVBit()
{
  m_bSwapped = m_bSwapped ? FALSE : TRUE;
}

SmStatus SmProceduralBlend::SplitAt
  (const SmContext & crContext,
   double            dParam,
   SmSurfParamType   eSurfParam,
   SmSurface      *& rpLeftSurface,
   SmSurface      *& rpRightSurface)
{
  rpLeftSurface = rpRightSurface = NULL;
  SmSurface * apNets[2] = { NULL, NULL };
  SER(SmBSplineSurface::SplitAt(crContext, dParam, eSurfParam, apNets[0], apNets[1]));
  SmObjDelete asDelete[2];
  asDelete[0].SetObj(apNets[0]);
  asDelete[1].SetObj(apNets[1]);
  SmSurface * apOut[2] = { NULL, NULL };
  for (int i=0; i<2; i++)
    {
      const SmBSplineSurface * pNet = SM_CAST_PTR(SmBSplineSurface, apNets[i]);
      NER(pNet);
      // each piece keeps the parametrization of this one, so the same definition and map
      SmProceduralBlend * pPiece = new (crContext) SmProceduralBlend(*pNet);
      NER(pPiece);
      SmObjDelete sDelete(pPiece);
      SER(pPiece->CopyDefinition(*this));
      sDelete.Clear();
      apOut[i] = pPiece;
    }
  rpLeftSurface = apOut[0];
  rpRightSurface = apOut[1];
  return SM_SUCCESS;
}

/*******************************************************************//**
PURPOSE: A rigid motion, with a uniform scale: the supports, the leader
    and the net move alike; the radius and the bound scale.
***********************************************************************/
SmStatus SmProceduralBlend::Transform(const SmAxis2Placement & crRotateNMove, const SmVector3d * cpOptScale)
{
  double dScale = 1.0;
  if (cpOptScale != NULL)
    {
      if (!SM_ARE_SAME(cpOptScale->x, cpOptScale->y) || !SM_ARE_SAME(cpOptScale->x, cpOptScale->z) || !(cpOptScale->x > SM_EFF_ZERO))
        { SER_MSG(SM_ERR, _T("SmProceduralBlend: a blend by its definition is scaled only uniformly")); }
      dScale = cpOptScale->x;
    }
  NER(m_apSupport[0]); NER(m_apSupport[1]); NER(m_pLeader);
  for (int i=0; i<2; i++)
    { SER(m_apSupport[i]->Transform(crRotateNMove, cpOptScale)); }
  SER(m_pLeader->Transform(crRotateNMove, cpOptScale));
  if (m_pLaw != NULL && dScale != 1.0)
    {
      SmAxis2Placement sIdentity;
      SER(m_pLaw->Transform(sIdentity, cpOptScale));
    }
  m_adR[0] *= dScale;
  m_adR[1] *= dScale;
  m_dBound *= dScale;
  Renew();
  SER(SmBSplineSurface::Transform(crRotateNMove, cpOptScale));
  return SM_SUCCESS;
}

/*******************************************************************//**
PURPOSE: Write the definition, then the net.
***********************************************************************/
SmStatus SmProceduralBlend::WriteToDB(SmDatabaseIO & rDB, ULONG lDBVersionNumber) const
{
  NER(m_apSupport[0]); NER(m_apSupport[1]); NER(m_pLeader); NER(m_apSeed[0]); NER(m_apSeed[1]);
  const SmFileType eType = rDB.GetFileType();
  const SmBoolean bHasLaw = (m_pLaw != NULL) ? TRUE : FALSE;
  const double adValues[10] = { m_adSide[0], m_adSide[1], m_adR[0], m_adR[1], m_dThumbweight, m_dBound,
                                m_adA[0], m_adB[0], m_adA[1], m_adB[1] };
  if (eType == SM_ASCII)
    {
      std::ostream & rOut = *rDB.GetOutStreamPtr();
      rOut << m_lOrder << " " << m_bSwapped << " " << bHasLaw << " SmProceduralBlend order, swapped, has a law\n";
      for (int i=0; i<10; i++) { rOut << adValues[i] << " "; }
      rOut << " SmProceduralBlend sides, offsets, thumbweight, bound, map\n";
    }
  else
    {
      SER(rDB.WriteLong(static_cast<long>(m_lOrder)));
      SER(rDB.WriteBoolean(m_bSwapped));
      SER(rDB.WriteBoolean(bHasLaw));
      SER(rDB.WriteDoubles(adValues, 10));
    }
  for (int i=0; i<2; i++)
    {
      if (eType == SM_ASCII) { *rDB.GetOutStreamPtr() << " SmProceduralBlend->Support\n"; }
      SER(rDB.WriteType(m_apSupport[i]->GetType()));
      SER(m_apSupport[i]->WriteToDB(rDB, lDBVersionNumber));
    }
  if (eType == SM_ASCII) { *rDB.GetOutStreamPtr() << " SmProceduralBlend->Leader\n"; }
  SER(sm_WriteCurve(rDB, *m_pLeader, lDBVersionNumber));
  for (int i=0; i<2; i++)
    {
      if (eType == SM_ASCII) { *rDB.GetOutStreamPtr() << " SmProceduralBlend->Seed\n"; }
      SER(sm_WriteCurve(rDB, *m_apSeed[i], lDBVersionNumber));
    }
  if (bHasLaw)
    {
      if (eType == SM_ASCII) { *rDB.GetOutStreamPtr() << " SmProceduralBlend->Law\n"; }
      SER(sm_WriteCurve(rDB, *m_pLaw, lDBVersionNumber));
    }
  SER(SmBSplineSurface::WriteToDB(rDB, lDBVersionNumber));
  return SM_SUCCESS;
} // end SmProceduralBlend::WriteToDB

/*******************************************************************//**
PURPOSE: Read a blend's definition and net; a hostile record (an order out
    of range, a side other than +-1, a radius or bound not positive, a
    depth past the cap) is refused as Create refuses it.
***********************************************************************/
SmStatus SmProceduralBlend::ReadFromDB
  (SM_TYPE           /*lType*/,
   SmDatabaseIO    & rDB,
   const SmContext & crContext,
   SmSurface      *& rpNewSurface,
   ULONG             lDBVersionNumber)
{
  if (rpNewSurface != NULL)
    { SER_MSG(SM_ERR, _T("SmProceduralBlend: read into a new surface only")); }
  const SmFileType eType = rDB.GetFileType();
  ULONG lOrder = 0;
  SmBoolean bSwapped = FALSE, bHasLaw = FALSE;
  double adValues[10] = { 0.0 };
  if (eType == SM_ASCII)
    {
      std::istream & rIn = *rDB.GetInStreamPtr();
      rIn >> lOrder >> bSwapped >> bHasLaw;
      rDB.GoToNextLine();
      for (int i=0; i<10; i++) { rIn >> adValues[i]; }
      rDB.GoToNextLine();
      if (!rIn)
        { SER_MSG(SM_ERR, _T("SmProceduralBlend: its definition could not be read")); }
    }
  else
    {
      long lRead = 0;
      SER(rDB.ReadLong(lRead));
      if (lRead < 0)
        { SER_MSG(SM_ERR_INVALID_INPUT, _T("SmProceduralBlend: a section's order is 1 to 4")); }
      lOrder = static_cast<ULONG>(lRead);
      SER(rDB.ReadBoolean(bSwapped));
      SER(rDB.ReadBoolean(bHasLaw));
      SER(rDB.ReadDoubles(adValues, 10));
    }
  SmSurface * apSupport[2] = { NULL, NULL };
  SmBSplineCurve * pLeader = NULL, * apSeed[2] = { NULL, NULL }, * pLaw = NULL;
  SmObjDelete asDelete[6];
  for (int i=0; i<2; i++)
    {
      if (eType == SM_ASCII) { rDB.GoToNextLine(); }
      SM_TYPE lSupportType;
      SER(rDB.ReadType(lSupportType));
      SER(SmSurface::ReadFromDB(lSupportType, rDB, crContext, apSupport[i], lDBVersionNumber));
      asDelete[i].SetObj(apSupport[i]);
    }
  if (eType == SM_ASCII) { rDB.GoToNextLine(); }
  SER(sm_ReadCurve(rDB, crContext, lDBVersionNumber, pLeader));
  asDelete[2].SetObj(pLeader);
  for (int i=0; i<2; i++)
    {
      if (eType == SM_ASCII) { rDB.GoToNextLine(); }
      SER(sm_ReadCurve(rDB, crContext, lDBVersionNumber, apSeed[i]));
      asDelete[3 + i].SetObj(apSeed[i]);
    }
  if (bHasLaw)
    {
      if (eType == SM_ASCII) { rDB.GoToNextLine(); }
      SER(sm_ReadCurve(rDB, crContext, lDBVersionNumber, pLaw));
      asDelete[5].SetObj(pLaw);
    }
  SmSurface * pNet = NULL;
  SER(SmBSplineSurface::ReadFromDB(SmBSplineSurface_TYPE, rDB, crContext, pNet, lDBVersionNumber));
  SmObjDelete sDeleteNet(pNet);
  const SmBSplineSurface * pNetBS = SM_CAST_PTR(SmBSplineSurface, pNet);
  NER(pNetBS);
  for (int i=0; i<10; i++)
    if (!sm_Finite(adValues[i]))
      { SER_MSG(SM_ERR_INVALID_INPUT, _T("SmProceduralBlend: its definition holds a number that is not finite")); }
  if (!(adValues[7] != 0.0 && adValues[9] != 0.0))
    { SER_MSG(SM_ERR_INVALID_INPUT, _T("SmProceduralBlend: its parametrization is singular")); }
  const SmSurface * apConst[2] = { apSupport[0], apSupport[1] };
  const SmBSplineCurve * apConstSeed[2] = { apSeed[0], apSeed[1] };
  const double adSide[2] = { adValues[0], adValues[1] };
  const double adOffset[2] = { adValues[2], adValues[3] };
  SmProceduralBlend * pBlend = NULL;
  SER(Create(crContext, *pNetBS, apConst, adSide, adOffset, pLaw, lOrder, adValues[4], *pLeader, apConstSeed, adValues[5], pBlend));
  NER(pBlend);
  pBlend->m_adA[0] = adValues[6];
  pBlend->m_adB[0] = adValues[7];
  pBlend->m_adA[1] = adValues[8];
  pBlend->m_adB[1] = adValues[9];
  pBlend->m_bSwapped = bSwapped ? TRUE : FALSE;
  rpNewSurface = pBlend;
  return SM_SUCCESS;
} // end SmProceduralBlend::ReadFromDB

SmBoolean SmProceduralBlend::IsKindOf(SM_TYPE t) const
{
  return ((SmProceduralBlend_TYPE == t) ? TRUE : SmBSplineSurface::IsKindOf(t));
}

void SmProceduralBlend::Dump(void) const
{
  TCHAR sBuff[SM_TBLOCK_SIZE];
  smos_sprintf(sBuff, _T("\nSmProceduralBlend: order %lu, offsets %g %g%s, thumbweight %g, bound %g, depth %lu"),
               static_cast<unsigned long>(m_lOrder), m_adR[0], m_adR[1], m_pLaw != NULL ? _T(" (a law)") : _T(""), m_dThumbweight, m_dBound,
               static_cast<unsigned long>(GetDepth()));
  smos_WriteBuffer(sBuff);
  SmBSplineSurface::Dump();
}
