// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

/**********************************************************************//**
* FILE NAME --- SmFilletSection.cpp
* PURPOSE: SolidCraft change. The cross-section of a blend by the order n of
*   its contact with the faces along the rolling ball's contact lines (the
*   rails): SmGnSectionFSG.
*
* The G^n section (n >= 2), in the plane square to the leader (the centre
* line) at its parameter s: V(tau) is the point of that plane at distance
* tau r_i from face i (r_i the ball's offset from face i), F_i(tau) its foot
* on face i (tau = 1 at the rails, 0 at the edge), and
* T_i(u) = sum_{k<=n} c_ik (u - 1)^k the order-n jet at tau = 1 of a curve on
* face i: c_ik = F_i^(k)(1) / k! for k < n, and c_in the part along the
* face's normal of the n-th coefficient of S_i(x_i(tau)), x_i the feet's
* parameters to order n - 1 (sm_GnTopJets). The feet's n-th coefficient
* would take the face's partials of order n + 1; this one takes them to
* order n, so the section is continuous on a C^n face. The section is
*
*   S(t) = f(t) T1(1 - p1(t)) + (1 - f(t)) T2(p2(t)),  t in [0, 1],
*   f(t) = (1-t)^(n+1) / ((1-t)^(n+1) + t^(n+1)),
*   p1(t) = t / (1 - w + w t),  p2(t) = t (1 - w) / (1 - w t),
*
* w the thumbweight. f's first n derivatives vanish at both ends, so S has
* the faces' own derivatives across each rail up to order n. Over the common
* denominator W(t) = ((1-t)^(n+1) + t^(n+1)) (1 - w + w t)^n (1 - w t)^n, of
* degree 3n + 1 with positive Bernstein coefficients, S is a rational curve
* whose control points are linear in the jets c_ik and whose weights depend
* on t only. The blend is therefore the tensor product of the jets' curves
* along s (piecewise cubic Hermite interpolants of the jets' exact values and
* s-derivatives at the stations, all on the same knots) with that one
* rational section; a jet constant along s (a straight spine) is carried
* exactly. The equal area thumbweight varies along s, w(s) resolved at each
* station from its own wedge: the surface then interpolates the section's
* homogeneous rows W_j(w(s)) P_j(s) and W_j(w(s)), each a cubic Hermite
* curve of its values and s-derivatives at the stations (dw/ds from the
* implicit function I(w) = target), and the certificate bounds those rows'
* remainder too (sm_GnRowTermsOver).
*
* The jets and their s-derivatives come from the system G(x; s, tau) = 0 in
* the feet's parameters x = (u1, v1, u2, v2):
*   S1(x1) + tau r_1 s_1 n1(x1) - S2(x2) - tau r_2 s_2 n2(x2) = 0,
*   (S1(x1) + tau r_1 s_1 n1(x1) - L(s)) . L'(s) / |L'(s)| = 0,
* s_i the side of face i the ball is on. Its Taylor coefficients in
* (s, tau - 1) are solved order by order with the Jacobian at the station:
* truncated power series arithmetic, exact to rounding, no differences.
**********************************************************************/

#include <SmFilletSolver.h>
#include <SmBSplineCurve.h>
#include <SmBSplineSurface.h>
#include <SmOffsetSurface.h>
#include <SmAxis2Placement.h>
#include <SmFilletExecutive.h>
#include <SmFilletIntersector.h>
#include <SmGraphicsExtern.h>
#include <SmGeomUtility.h>
#include <SmBrep.h>
#include <SmProceduralBlend.h>
#include <SmFilletStandardSolver.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <vector>

namespace
{

inline SmBoolean sm_Positive(double x) { return x > 0.0 ? TRUE : FALSE; }
inline double    sm_Sqrt(double x)     { return smos_Sqrt(x); }

/*******************************************************************//**
   Rounded interval arithmetic: every result widened outward past its
   rounding error, which holds the exact result of the operation on any
   numbers in the operands. An empty or undefined result (a division by an
   interval holding 0, a NaN) is the whole line, which no test passes.
***********************************************************************/
struct sm_Ivl
{
  double lo, hi;
  sm_Ivl() : lo(0.0), hi(0.0) {}
  sm_Ivl(double d) : lo(d), hi(d) {}
  sm_Ivl(double a, double b) : lo(a), hi(b) {}
};

// An operation's result rounded to nearest is within half a unit in the last place of the exact one,
// |x| 2^-53, or 2^-1075 below the normal range; moving it out by twice that holds the exact result
// (the move's own rounding is inside the margin).
inline double sm_Dn(double x) { return x - (std::fabs(x) * 2.3e-16 + 1.0e-300); }
inline double sm_Up(double x) { return x + (std::fabs(x) * 2.3e-16 + 1.0e-300); }
inline sm_Ivl sm_Whole() { return sm_Ivl(-HUGE_VAL, HUGE_VAL); }
inline sm_Ivl sm_Checked(double a, double b)
{
  return (std::isnan(a) || std::isnan(b)) ? sm_Whole() : sm_Ivl(sm_Dn(a), sm_Up(b));
}
inline sm_Ivl operator+(const sm_Ivl & a, const sm_Ivl & b) { return sm_Checked(a.lo + b.lo, a.hi + b.hi); }
inline sm_Ivl operator-(const sm_Ivl & a, const sm_Ivl & b) { return sm_Checked(a.lo - b.hi, a.hi - b.lo); }
inline sm_Ivl operator*(const sm_Ivl & a, const sm_Ivl & b)
{
  const double p[4] = { a.lo * b.lo, a.lo * b.hi, a.hi * b.lo, a.hi * b.hi };
  double lo = p[0], hi = p[0];
  for (int i=1; i<4; i++)
    {
      if (std::isnan(p[i])) { return sm_Whole(); }
      lo = std::min(lo, p[i]);
      hi = std::max(hi, p[i]);
    }
  return sm_Checked(lo, hi);
}
inline sm_Ivl operator/(const sm_Ivl & a, const sm_Ivl & b)
{
  if (!(b.lo > 0.0 || b.hi < 0.0))
    { return sm_Whole(); }
  return a * sm_Checked(1.0 / b.hi, 1.0 / b.lo);
}
inline SmBoolean sm_Positive(const sm_Ivl & x) { return x.lo > 0.0 ? TRUE : FALSE; }
inline sm_Ivl    sm_Sqrt(const sm_Ivl & x)
{
  if (!(x.hi >= 0.0)) { return sm_Whole(); }
  return sm_Checked(smos_Sqrt(std::max(x.lo, 0.0)), smos_Sqrt(x.hi));
}
inline double sm_Mag(const sm_Ivl & x) { return std::max(std::fabs(x.lo), std::fabs(x.hi)); }
inline double sm_Mid(const sm_Ivl & x) { return 0.5 * (x.lo + x.hi); }
inline sm_Ivl sm_Hull(const sm_Ivl & a, const sm_Ivl & b) { return sm_Ivl(std::min(a.lo, b.lo), std::max(a.hi, b.hi)); }

/*******************************************************************//**
   Truncated power series in two variables, sigma (along the leader's
   parameter s) and eps (tau - 1): c[i][k] the coefficient of sigma^i eps^k,
   kept up to the orders of an sm_Trunc.
***********************************************************************/
const int SM_GN_MAXS = 7;
// the Taylor form's order of the certificate's fourth s-derivatives: their Taylor polynomial of
// degree SM_GN_TAYLOR - 1 at the span's middle, and the remainder bound over the span
const int SM_GN_TAYLOR = 2;
const int SM_GN_MAXE = 4;

struct sm_Trunc
{
  int s; // highest power of sigma kept
  int e; // highest power of eps kept
};

template <class T> struct sm_Ser
{
  T c[SM_GN_MAXS + 1][SM_GN_MAXE + 1];
};

template <class T> sm_Ser<T> sm_SerConst(const T & crV)
{
  sm_Ser<T> sR;
  for (int i=0; i<=SM_GN_MAXS; i++)
    for (int k=0; k<=SM_GN_MAXE; k++)
      sR.c[i][k] = T(0.0);
  sR.c[0][0] = crV;
  return sR;
}

template <class T> sm_Ser<T> sm_SerAdd(const sm_Ser<T> & a, const sm_Ser<T> & b, const sm_Trunc & tr)
{
  sm_Ser<T> sR = sm_SerConst<T>(T(0.0));
  for (int i=0; i<=tr.s; i++)
    for (int k=0; k<=tr.e; k++)
      sR.c[i][k] = a.c[i][k] + b.c[i][k];
  return sR;
}

template <class T> sm_Ser<T> sm_SerSub(const sm_Ser<T> & a, const sm_Ser<T> & b, const sm_Trunc & tr)
{
  sm_Ser<T> sR = sm_SerConst<T>(T(0.0));
  for (int i=0; i<=tr.s; i++)
    for (int k=0; k<=tr.e; k++)
      sR.c[i][k] = a.c[i][k] - b.c[i][k];
  return sR;
}

template <class T> sm_Ser<T> sm_SerScale(const sm_Ser<T> & a, const T & crF, const sm_Trunc & tr)
{
  sm_Ser<T> sR = sm_SerConst<T>(T(0.0));
  for (int i=0; i<=tr.s; i++)
    for (int k=0; k<=tr.e; k++)
      sR.c[i][k] = a.c[i][k] * crF;
  return sR;
}

template <class T> sm_Ser<T> sm_SerMul(const sm_Ser<T> & a, const sm_Ser<T> & b, const sm_Trunc & tr)
{
  sm_Ser<T> sR = sm_SerConst<T>(T(0.0));
  for (int i=0; i<=tr.s; i++)
    for (int k=0; k<=tr.e; k++)
      {
        T sSum(0.0);
        for (int i1=0; i1<=i; i1++)
          for (int k1=0; k1<=k; k1++)
            sSum = sSum + a.c[i1][k1] * b.c[i-i1][k-k1];
        sR.c[i][k] = sSum;
      }
  return sR;
}

// x^(-1/2) for a series whose constant term is positive: x0^(-1/2) (1 + z)^(-1/2), z = x / x0 - 1,
// summed as the binomial series in z, which has no constant term.
template <class T> SmStatus sm_SerInvSqrt(const sm_Ser<T> & x, const sm_Trunc & tr, sm_Ser<T> & rY)
{
  const T x0 = x.c[0][0];
  if (!sm_Positive(x0))
    { SER_MSG(SM_ERR, _T("SmGnSectionFSG: a face's normal vanishes at a foot")); }
  sm_Ser<T> z = sm_SerScale(x, T(1.0) / x0, tr);
  z.c[0][0] = T(0.0);
  sm_Ser<T> sPow = sm_SerConst<T>(T(1.0));
  sm_Ser<T> sSum = sm_SerConst<T>(T(1.0));
  double dBinom = 1.0;
  for (int m=1; m<=tr.s+tr.e; m++)
    {
      dBinom *= (-0.5 - (m - 1)) / m;
      sPow = sm_SerMul(sPow, z, tr);
      sSum = sm_SerAdd(sSum, sm_SerScale(sPow, T(dBinom), tr), tr);
    }
  rY = sm_SerScale(sSum, T(1.0) / sm_Sqrt(x0), tr);
  return SM_SUCCESS;
}

template <class T> struct sm_V3
{
  T x, y, z;
};

template <class T> sm_V3<T> sm_V3From(const SmVector3d & v)
{
  sm_V3<T> r = { T(v.x), T(v.y), T(v.z) };
  return r;
}

template <class T> struct sm_VSer
{
  sm_Ser<T> a[3];
};

template <class T> sm_VSer<T> sm_VSerConst(const sm_V3<T> & v)
{
  sm_VSer<T> r;
  r.a[0] = sm_SerConst<T>(v.x);
  r.a[1] = sm_SerConst<T>(v.y);
  r.a[2] = sm_SerConst<T>(v.z);
  return r;
}

template <class T> sm_VSer<T> sm_VSerAdd(const sm_VSer<T> & p, const sm_VSer<T> & q, const sm_Trunc & tr)
{
  sm_VSer<T> r;
  for (int j=0; j<3; j++) { r.a[j] = sm_SerAdd(p.a[j], q.a[j], tr); }
  return r;
}

template <class T> sm_VSer<T> sm_VSerSub(const sm_VSer<T> & p, const sm_VSer<T> & q, const sm_Trunc & tr)
{
  sm_VSer<T> r;
  for (int j=0; j<3; j++) { r.a[j] = sm_SerSub(p.a[j], q.a[j], tr); }
  return r;
}

template <class T> sm_VSer<T> sm_VSerTimes(const sm_VSer<T> & p, const sm_Ser<T> & s, const sm_Trunc & tr)
{
  sm_VSer<T> r;
  for (int j=0; j<3; j++) { r.a[j] = sm_SerMul(p.a[j], s, tr); }
  return r;
}

template <class T> sm_Ser<T> sm_VSerDot(const sm_VSer<T> & p, const sm_VSer<T> & q, const sm_Trunc & tr)
{
  sm_Ser<T> r = sm_SerMul(p.a[0], q.a[0], tr);
  r = sm_SerAdd(r, sm_SerMul(p.a[1], q.a[1], tr), tr);
  return sm_SerAdd(r, sm_SerMul(p.a[2], q.a[2], tr), tr);
}

template <class T> sm_VSer<T> sm_VSerCross(const sm_VSer<T> & p, const sm_VSer<T> & q, const sm_Trunc & tr)
{
  sm_VSer<T> r;
  r.a[0] = sm_SerSub(sm_SerMul(p.a[1], q.a[2], tr), sm_SerMul(p.a[2], q.a[1], tr), tr);
  r.a[1] = sm_SerSub(sm_SerMul(p.a[2], q.a[0], tr), sm_SerMul(p.a[0], q.a[2], tr), tr);
  r.a[2] = sm_SerSub(sm_SerMul(p.a[0], q.a[1], tr), sm_SerMul(p.a[1], q.a[0], tr), tr);
  return r;
}

/*******************************************************************//**
   A face's partial derivatives at a point: P[a][b] = d^(a+b) S / du^a dv^b,
   a + b <= the order asked for.
***********************************************************************/
const int SM_GN_MAXK = SM_GN_MAXS + SM_GN_MAXE + 2;

template <class T> struct sm_Partials
{
  sm_V3<T> P[SM_GN_MAXK + 1][SM_GN_MAXK + 1];
};

/*******************************************************************//**
   A NURBS surface's data, for evaluation of its partials at any (u, v): a
   point outside the knots' range is on the analytic continuation of the
   end span (a plane goes on as the plane, a cylinder's arc as its circle),
   which is the face's own surface beyond its trimmed domain.
***********************************************************************/
struct sm_Nurbs
{
  int                 p, q;   // degrees in u and v
  int                 nu, nv; // control points each way
  std::vector<double> U, V;   // full knot vectors
  std::vector<double> Pw;     // homogeneous control points (x w, y w, z w, w), u-major
};

SmStatus sm_NurbsOf(const SmBSplineSurface & crSurface, sm_Nurbs & rN)
{
  ULONG lKU = 0, lKV = 0, lNU = 0, lNV = 0;
  double * pU = NULL, * pV = NULL, * pCP = NULL;
  SER(crSurface.GetKnotsPointers(lKU, lKV, pU, pV));
  SER(crSurface.GetControlPointsPointer(lNU, lNV, pCP));
  NER(pU); NER(pV); NER(pCP);
  rN.p  = static_cast<int>(crSurface.GetDegree(SM_SP_U));
  rN.q  = static_cast<int>(crSurface.GetDegree(SM_SP_V));
  rN.nu = static_cast<int>(lNU);
  rN.nv = static_cast<int>(lNV);
  if (rN.p < 1 || rN.q < 1 || static_cast<int>(lKU) != rN.nu + rN.p + 1 || static_cast<int>(lKV) != rN.nv + rN.q + 1)
    { SER_MSG(SM_ERR, _T("SmGnSectionFSG: a face's NURBS data do not agree")); }
  rN.U.assign(pU, pU + lKU);
  rN.V.assign(pV, pV + lKV);
  rN.Pw.resize(static_cast<size_t>(4 * lNU * lNV));
  for (size_t i=0; i<static_cast<size_t>(lNU * lNV); i++)
    {
      const double w = pCP[4*i+3];
      const SmBoolean bRational = (w != NL_NOW) ? TRUE : FALSE;
      rN.Pw[4*i]   = pCP[4*i];
      rN.Pw[4*i+1] = pCP[4*i+1];
      rN.Pw[4*i+2] = pCP[4*i+2];
      rN.Pw[4*i+3] = bRational ? w : 1.0;
      if (!(rN.Pw[4*i+3] > 0.0))
        { SER_MSG(SM_ERR, _T("SmGnSectionFSG: a face's weights are not positive")); }
    }
  return SM_SUCCESS;
}

// The span k (U[k] <= u < U[k+1], p <= k < n) of a parameter, the end spans continued past the
// knots' range.
int sm_Span(const std::vector<double> & crU, int p, int n, double u)
{
  int k = p;
  while (k + 1 < n && u >= crU[static_cast<size_t>(k + 1)])
    { k++; }
  while (k > p && crU[static_cast<size_t>(k)] >= crU[static_cast<size_t>(k + 1)])
    { k--; }
  return k;
}

// The non-zero basis functions N_{k-p..k} and their derivatives up to lD at u in span k:
// rDers[d][j] (the NURBS book's A2.3, its knot differences taken from the knots directly).
template <class T> void sm_DersBasis(const std::vector<double> & crU, int k, int p, const T & u, int lD, std::vector< std::vector<T> > & rDers)
{
  std::vector< std::vector<T> > ndu(static_cast<size_t>(p + 1), std::vector<T>(static_cast<size_t>(p + 1), T(0.0)));
  std::vector<T> left(static_cast<size_t>(p + 1), T(0.0)), right(static_cast<size_t>(p + 1), T(0.0));
  ndu[0][0] = T(1.0);
  for (int j=1; j<=p; j++)
    {
      left[j]  = u - T(crU[static_cast<size_t>(k + 1 - j)]);
      right[j] = T(crU[static_cast<size_t>(k + j)]) - u;
      T saved(0.0);
      for (int r=0; r<j; r++)
        {
          ndu[j][r] = T(crU[static_cast<size_t>(k + r + 1)] - crU[static_cast<size_t>(k + 1 - j + r)]);
          const T temp = ndu[r][j-1] / ndu[j][r];
          ndu[r][j] = saved + right[r+1] * temp;
          saved = left[j-r] * temp;
        }
      ndu[j][j] = saved;
    }
  rDers.assign(static_cast<size_t>(lD + 1), std::vector<T>(static_cast<size_t>(p + 1), T(0.0)));
  for (int j=0; j<=p; j++) { rDers[0][j] = ndu[j][p]; }
  std::vector< std::vector<T> > a(2, std::vector<T>(static_cast<size_t>(p + 1), T(0.0)));
  const int lDu = (lD < p) ? lD : p;
  for (int r=0; r<=p; r++)
    {
      int s1 = 0, s2 = 1;
      a[0][0] = T(1.0);
      for (int kk=1; kk<=lDu; kk++)
        {
          T d(0.0);
          const int rk = r - kk, pk = p - kk;
          if (r >= kk)
            {
              a[s2][0] = a[s1][0] / ndu[pk+1][rk];
              d = a[s2][0] * ndu[rk][pk];
            }
          const int j1 = (rk >= -1) ? 1 : -rk;
          const int j2 = (r - 1 <= pk) ? kk - 1 : p - r;
          for (int j=j1; j<=j2; j++)
            {
              a[s2][j] = (a[s1][j] - a[s1][j-1]) / ndu[pk+1][rk+j];
              d = d + a[s2][j] * ndu[rk+j][pk];
            }
          if (r <= pk)
            {
              a[s2][kk] = (T(0.0) - a[s1][kk-1]) / ndu[pk+1][r];
              d = d + a[s2][kk] * ndu[r][pk];
            }
          rDers[kk][r] = d;
          const int t = s1; s1 = s2; s2 = t;
        }
    }
  double f = p;
  for (int kk=1; kk<=lDu; kk++)
    {
      for (int j=0; j<=p; j++) { rDers[kk][j] = rDers[kk][j] * T(f); }
      f *= (p - kk);
    }
}

// The partials P[a][b], a + b <= lK, of the NURBS at (u, v) (the NURBS book's A4.4).
template <class T> void sm_NurbsPartials(const sm_Nurbs & crN, int ku, int kv, const T & u, const T & v, int lK, sm_Partials<T> & rP)
{
  std::vector< std::vector<T> > Nu, Nv;
  sm_DersBasis(crN.U, ku, crN.p, u, lK, Nu);
  sm_DersBasis(crN.V, kv, crN.q, v, lK, Nv);
  // the homogeneous partials A (x, y, z) and w
  std::vector< std::vector< sm_V3<T> > > A(static_cast<size_t>(lK + 1), std::vector< sm_V3<T> >(static_cast<size_t>(lK + 1)));
  std::vector< std::vector<T> > W(static_cast<size_t>(lK + 1), std::vector<T>(static_cast<size_t>(lK + 1), T(0.0)));
  for (int a=0; a<=lK; a++)
    for (int b=0; a+b<=lK; b++)
      {
        T sx(0.0), sy(0.0), sz(0.0), sw(0.0);
        for (int i=0; i<=crN.p; i++)
          for (int j=0; j<=crN.q; j++)
            {
              const size_t lI = static_cast<size_t>(4 * ((ku - crN.p + i) * crN.nv + (kv - crN.q + j)));
              const T c = Nu[a][i] * Nv[b][j];
              sx = sx + c * T(crN.Pw[lI]);
              sy = sy + c * T(crN.Pw[lI+1]);
              sz = sz + c * T(crN.Pw[lI+2]);
              sw = sw + c * T(crN.Pw[lI+3]);
            }
        A[a][b].x = sx; A[a][b].y = sy; A[a][b].z = sz;
        W[a][b] = sw;
      }
  for (int a=0; a<=lK; a++)
    for (int b=0; a+b<=lK; b++)
      {
        sm_V3<T> v3 = A[a][b];
        double cka = 1.0;
        for (int i=0; i<=a; i++)
          {
            if (i > 0) { cka = cka * (a - i + 1) / i; }
            double cbj = 1.0;
            for (int j=0; j<=b; j++)
              {
                if (j > 0) { cbj = cbj * (b - j + 1) / j; }
                if (i == 0 && j == 0) { continue; }
                const T c = T(cka * cbj) * W[i][j];
                v3.x = v3.x - c * rP.P[a-i][b-j].x;
                v3.y = v3.y - c * rP.P[a-i][b-j].y;
                v3.z = v3.z - c * rP.P[a-i][b-j].z;
              }
          }
        rP.P[a][b].x = v3.x / W[0][0];
        rP.P[a][b].y = v3.y / W[0][0];
        rP.P[a][b].z = v3.z / W[0][0];
      }
}

// The partials of a face's surface at crUV up to total order lK: a NURBS (every analytic surface is
// one, evaluated in its net's parameters) by its net, continued past its knots; any other surface,
// and a blend by its definition (whose net is only its approximation), by its own evaluator.
// lSideU, lSideV: at a knot line (within 1e-9 of it), -1 the piece before it, +1 the piece after
// it, 0 the span sm_Span finds.
int sm_SidedSpan(const std::vector<double> & crU, int p, int n, double u, int lSide)
{
  int k = sm_Span(crU, p, n, u);
  const double dTol = 1.0e-9 * (1.0 + std::fabs(u));
  if (lSide < 0 && std::fabs(u - crU[static_cast<size_t>(k)]) <= dTol)
    {
      while (k > p && !(crU[static_cast<size_t>(k - 1)] < crU[static_cast<size_t>(k)])) { k--; }
      if (k > p) { k--; }
    }
  else if (lSide > 0 && k + 1 < n && std::fabs(u - crU[static_cast<size_t>(k + 1)]) <= dTol)
    {
      k++;
      while (k + 1 < n && !(crU[static_cast<size_t>(k)] < crU[static_cast<size_t>(k + 1)])) { k++; }
    }
  return k;
}

SmStatus sm_SurfacePartials(const SmSurface & crSurface, const SmPoint2d & crUV, int lK, sm_Partials<double> & rP, int lSideU = 0,
                            int lSideV = 0)
{
  if (lK < 0 || lK > SM_GN_MAXK)
    { SER(SM_ERR); }
  const SmBSplineSurface * pBS = (crSurface.IsKindOf(SmBSplineSurface_TYPE) && !crSurface.IsKindOf(SmProceduralBlend_TYPE))
                                 ? static_cast<const SmBSplineSurface *>(&crSurface) : NULL;
  if (pBS != NULL)
    {
      sm_Nurbs sN;
      SER(sm_NurbsOf(*pBS, sN));
      const int ku = sm_SidedSpan(sN.U, sN.p, sN.nu, crUV.x, lSideU), kv = sm_SidedSpan(sN.V, sN.q, sN.nv, crUV.y, lSideV);
      sm_NurbsPartials<double>(sN, ku, kv, crUV.x, crUV.y, lK, rP);
      return SM_SUCCESS;
    }
  std::vector<SmVector3d> sD(static_cast<size_t>((lK + 1) * (lK + 1)));
  SER(crSurface.Evaluate(crUV, static_cast<ULONG>(lK), static_cast<ULONG>(lK), TRUE, TRUE, TRUE, sD.data(), FALSE));
  for (int a=0; a<=lK; a++)
    for (int b=0; a+b<=lK; b++)
      rP.P[a][b] = sm_V3From<double>(sD[static_cast<size_t>(a * (lK + 1) + b)]);
  return SM_SUCCESS;
}

/*******************************************************************//**
   The blend's system at a station: the two faces' partials at the feet,
   the leader's point and unit tangent as series in sigma, the offsets r_i
   and the sides s_i. bHold: the offsets held at tau = 1 (no eps in them).
***********************************************************************/
template <class T> struct sm_GnSystem
{
  const sm_Partials<T> * apP[2];
  sm_VSer<T>             sLeader;
  sm_VSer<T>             sTangent;
  double                 adR[2];
  double                 adSide[2];
  SmBoolean              bHold;
  const sm_Ser<T>      * pR; // a ball whose radius varies: r(s) as a series in sigma (adR unused), or NULL
};

// G and the feet F_i as series, given the feet's parameters as series x[0..3] (u1, v1, u2, v2),
// whose constant terms are the point the partials were taken at; pV: the point V(tau) off face 1
// (the ball's centre at tau = 1); pN: the faces' unit normals at the feet.
template <class T> SmStatus sm_GnEvaluate(const sm_GnSystem<T> & crSys, const sm_Ser<T> (&x)[4], const sm_Trunc & tr,
                                          sm_Ser<T> (&rG)[4], sm_VSer<T> (&rF)[2], sm_VSer<T> * pV = NULL, sm_VSer<T> (*pN)[2] = NULL)
{
  sm_VSer<T> asV[2];
  const int lM = tr.s + tr.e;
  for (int f=0; f<2; f++)
    {
      sm_Ser<T> du = x[2*f], dv = x[2*f+1];
      du.c[0][0] = T(0.0);
      dv.c[0][0] = T(0.0);
      std::vector< sm_Ser<T> > sDu, sDv;
      sDu.push_back(sm_SerConst<T>(T(1.0)));
      sDv.push_back(sm_SerConst<T>(T(1.0)));
      for (int a=1; a<=lM; a++)
        {
          sDu.push_back(sm_SerMul(sDu.back(), du, tr));
          sDv.push_back(sm_SerMul(sDv.back(), dv, tr));
        }
      const sm_Partials<T> & crP = *crSys.apP[f];
      const sm_V3<T> sZero = { T(0.0), T(0.0), T(0.0) };
      sm_VSer<T> sS  = sm_VSerConst<T>(sZero);
      sm_VSer<T> sSu = sS, sSv = sS;
      double dFa = 1.0;
      for (int a=0; a<=lM; a++)
        {
          if (a > 0) { dFa *= a; }
          double dFb = 1.0;
          for (int b=0; a+b<=lM; b++)
            {
              if (b > 0) { dFb *= b; }
              const sm_Ser<T> sMono = sm_SerMul(sDu[static_cast<size_t>(a)], sDv[static_cast<size_t>(b)], tr);
              const T sF(1.0 / (dFa * dFb));
              const sm_V3<T> * apC[3] = { &crP.P[a][b], &crP.P[a+1][b], &crP.P[a][b+1] };
              sm_VSer<T> * apOut[3] = { &sS, &sSu, &sSv };
              for (int o=0; o<3; o++)
                {
                  const T acV[3] = { apC[o]->x, apC[o]->y, apC[o]->z };
                  for (int j=0; j<3; j++)
                    { apOut[o]->a[j] = sm_SerAdd(apOut[o]->a[j], sm_SerScale(sMono, acV[j] * sF, tr), tr); }
                }
            }
        }
      const sm_VSer<T> sN = sm_VSerCross(sSu, sSv, tr);
      sm_Ser<T> sInvLen;
      SER(sm_SerInvSqrt(sm_VSerDot(sN, sN, tr), tr, sInvLen));
      // the offset towards the ball: tau r_i s_i, tau = 1 + eps (r_i a series in sigma for a law)
      sm_Ser<T> sD = sm_SerConst<T>(T(crSys.adR[f] * crSys.adSide[f]));
      if (crSys.pR != NULL)
        for (int i=0; i<=tr.s; i++)
          sD.c[i][0] = crSys.pR->c[i][0] * T(crSys.adSide[f]);
      if (!crSys.bHold && tr.e >= 1)
        for (int i=0; i<=tr.s; i++)
          sD.c[i][1] = sD.c[i][0];
      asV[f] = sm_VSerAdd(sS, sm_VSerTimes(sN, sm_SerMul(sInvLen, sD, tr), tr), tr);
      if (pN != NULL)
        { (*pN)[f] = sm_VSerTimes(sN, sInvLen, tr); }
      rF[f] = sS;
    }
  const sm_VSer<T> sDiff = sm_VSerSub(asV[0], asV[1], tr);
  for (int j=0; j<3; j++) { rG[j] = sDiff.a[j]; }
  rG[3] = sm_VSerDot(sm_VSerSub(asV[0], crSys.sLeader, tr), crSys.sTangent, tr);
  if (pV != NULL)
    { *pV = asV[0]; }
  return SM_SUCCESS;
}

// The Jacobian dG/dx at x's constant terms, offsets held: G along x + e_j eps, first order.
template <class T> SmStatus sm_GnJacobian(const sm_GnSystem<T> & crSys, const sm_Ser<T> (&x)[4], T (&rJ)[4][4])
{
  const sm_Trunc tr = { 0, 1 };
  sm_GnSystem<T> sHeld = crSys;
  sHeld.bHold = TRUE;
  for (int j=0; j<4; j++)
    {
      sm_Ser<T> sX[4];
      for (int q=0; q<4; q++)
        {
          sX[q] = sm_SerConst<T>(x[q].c[0][0]);
          if (q == j) { sX[q].c[0][1] = T(1.0); }
        }
      sm_Ser<T> sG[4];
      sm_VSer<T> sF[2];
      SER(sm_GnEvaluate(sHeld, sX, tr, sG, sF));
      for (int i=0; i<4; i++) { rJ[i][j] = sG[i].c[0][1]; }
    }
  return SM_SUCCESS;
}

/*******************************************************************//**
   The section's jets of order lTop (its order n) and up, from the faces'
   order-n partials only. The feet's parameters x are solved to order n - 1
   in eps (they take the faces' partials to order n), and the face curve
   C_i(tau) = S_i(x_i(tau)), x_i so truncated, is a curve on face i whose jet
   agrees with the feet's to order n - 1. Its eps^n coefficient's part along
   the face's normal is the face's own (the second fundamental form and its
   order-n kin, along the foot's direction), and its part in the tangent
   plane is only the parametrization's: a curve on the face has any tangent
   n-th coefficient, by reparametrizing it (DeRose's G^n). The jet of order
   k >= n is that normal part, N (N . C_k), so the section meets the face
   in contact of order n + 1 (n + 2 between planes, where it is the feet's
   own, zero), and is as continuous as the face's n-th partials: on a C^n
   face (a cubic B-spline's C2 for n = 2) it is continuous across the face's
   knot lines. rF holds C (x truncated) on entry, the jets on return; crN the
   faces' unit normals at the feet.
***********************************************************************/
template <class T> void sm_GnTopJets(int lTop, const sm_Trunc & tr, const sm_VSer<T> (&crN)[2], sm_VSer<T> (&rF)[2])
{
  if (lTop < 1)
    { return; }
  const sm_Trunc sSigma = { tr.s, 0 };
  for (int f=0; f<2; f++)
    {
      sm_VSer<T> sN = crN[f];
      for (int c=0; c<3; c++)
        for (int i=0; i<=SM_GN_MAXS; i++)
          for (int k=1; k<=SM_GN_MAXE; k++)
            sN.a[c].c[i][k] = T(0.0);
      for (int k=lTop; k<=tr.e; k++)
        {
          sm_VSer<T> sC = sN;
          for (int c=0; c<3; c++)
            for (int i=0; i<=SM_GN_MAXS; i++)
              for (int k2=0; k2<=SM_GN_MAXE; k2++)
                sC.a[c].c[i][k2] = (k2 == 0 && i <= tr.s) ? rF[f].a[c].c[i][k] : T(0.0);
          const sm_VSer<T> sJ = sm_VSerTimes(sN, sm_VSerDot(sN, sC, sSigma), sSigma);
          for (int c=0; c<3; c++)
            for (int i=0; i<=tr.s; i++)
              rF[f].a[c].c[i][k] = sJ.a[c].c[i][0];
        }
    }
}

// Solves A y = b, 4 x 4, by Gaussian elimination with partial pivoting.
SmStatus sm_Solve4(const double (&crA)[4][4], const double (&crB)[4], double (&rY)[4])
{
  double a[4][5];
  for (int i=0; i<4; i++)
    {
      for (int j=0; j<4; j++) { a[i][j] = crA[i][j]; }
      a[i][4] = crB[i];
    }
  double dScale = 0.0;
  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++)
      dScale = smos_Max(dScale, smos_Fabs(a[i][j]));
  for (int c=0; c<4; c++)
    {
      int p = c;
      for (int r=c+1; r<4; r++)
        if (smos_Fabs(a[r][c]) > smos_Fabs(a[p][c])) { p = r; }
      if (!(smos_Fabs(a[p][c]) > 1.0e-13 * dScale))
        { SER_MSG(SM_ERR, _T("SmGnSectionFSG: the blend's system is singular at a station (its faces' normals agree there)")); }
      for (int j=0; j<5; j++) { double t = a[c][j]; a[c][j] = a[p][j]; a[p][j] = t; }
      for (int r=c+1; r<4; r++)
        {
          const double m = a[r][c] / a[c][c];
          for (int j=c; j<5; j++) { a[r][j] -= m * a[c][j]; }
        }
    }
  for (int c=3; c>=0; c--)
    {
      double v = a[c][4];
      for (int j=c+1; j<4; j++) { v -= a[c][j] * rY[j]; }
      rY[c] = v / a[c][c];
    }
  return SM_SUCCESS;
}

/*******************************************************************//**
   One station of the blend: the leader's parameter, the feet's series and
   the faces' series there.
***********************************************************************/
struct sm_GnStation
{
  double            dS;
  sm_Ser<double>    asX[4];  // u1, v1, u2, v2 in (sigma, eps)
  sm_VSer<double>   asF[2];  // the feet on faces 1 and 2
  sm_VSer<double>   sV;      // V, the point at tau r_i from the faces (the ball's centre at tau = 1)
  double            dW;      // the equal area thumbweight there, and its s-derivative (NaN: not resolved)
  double            dWds;
  int               lSteps;  // the Newton steps the feet took
  // the series above are the station's on the pieces after it along s (the leader's, a face's
  // where a foot is on one of its knot lines); these on the pieces before it, whose s-derivatives
  // the span ending here takes (the values agree; the derivatives differ where the leader or a
  // face is not smooth at the station); aiBefore: a face's side for them (sm_SidedSpan)
  int               aiBefore[4];
  SmBoolean         bSided;
  sm_Ser<double>    asXL[4];
  sm_VSer<double>   asFL[2];
  sm_VSer<double>   sVL;
  double            dWdsL;
};

// a station's series as the span ending there takes them
inline const sm_VSer<double> & sm_InF(const sm_GnStation & crSt, int f) { return crSt.bSided ? crSt.asFL[f] : crSt.asF[f]; }
inline const sm_Ser<double> & sm_InX(const sm_GnStation & crSt, int q) { return crSt.bSided ? crSt.asXL[q] : crSt.asX[q]; }
inline double sm_InWds(const sm_GnStation & crSt) { return crSt.bSided ? crSt.dWdsL : crSt.dWds; }

/*******************************************************************//**
   The leader (centre line) at s as series in sigma up to tr.s: its point,
   and its unit tangent (which takes one derivative more).
***********************************************************************/
SmStatus sm_LeaderSeries(const SmBSplineCurve & crLeader, double dS, const sm_Trunc & tr, sm_VSer<double> & rL, sm_VSer<double> & rT,
                         SmBoolean bBefore = FALSE)
{
  const sm_Trunc sSigma = { tr.s, 0 };
  SmVector3d asD[SM_GN_MAXS + 2];
  // a polynomial leader's derivatives past its degree are zero; a rational one's are evaluated
  ULONG lDerivs = static_cast<ULONG>(tr.s + 1);
  if (!crLeader.IsRational())
    { lDerivs = smos_Min(lDerivs, crLeader.GetDegree()); }
  for (ULONG j=lDerivs+1; j<=static_cast<ULONG>(tr.s + 1); j++) { asD[j] = SmVector3d(0.0, 0.0, 0.0); }
  // at a knot of the leader its piece after dS (upper interval), or the one before it with bBefore
  if (crLeader.Evaluate(dS, lDerivs, bBefore ? FALSE : TRUE, asD, FALSE) != SM_SUCCESS)
    { SER_MSG(SM_ERR, _T("SmGnSectionFSG: the leader's derivatives of this order are not evaluated")); }
  const sm_V3<double> sZero = { 0.0, 0.0, 0.0 };
  rL = sm_VSerConst<double>(sZero);
  sm_VSer<double> sD1 = rL;
  double dFact = 1.0;
  for (int j=0; j<=tr.s; j++)
    {
      if (j > 0) { dFact *= j; }
      const double acL[3] = { asD[j].x, asD[j].y, asD[j].z };
      const double acD[3] = { asD[j+1].x, asD[j+1].y, asD[j+1].z };
      for (int c=0; c<3; c++)
        {
          rL.a[c].c[j][0]  = acL[c] / dFact;
          sD1.a[c].c[j][0] = acD[c] / dFact;
        }
    }
  sm_Ser<double> sInv;
  SER(sm_SerInvSqrt(sm_VSerDot(sD1, sD1, sSigma), sSigma, sInv));
  rT = sm_VSerTimes(sD1, sInv, sSigma);
  return SM_SUCCESS;
}

/*******************************************************************//**
   The station at leader parameter dS: the feet solved by Newton from the
   seeds crSeed (u1, v1, u2, v2), then their series in (sigma, eps) up to tr,
   order by order with the station's Jacobian.
***********************************************************************/
SmStatus sm_GnSolveStation(const SmSurface * apFace[2], const SmBSplineCurve & crLeader, double dS, const double (&crSeed)[4],
                           const double (&crR)[2], const double (&crSide)[2], const sm_Trunc & tr, double dTol,
                           sm_GnStation & rStation, const sm_Ser<double> * pR = NULL, int lTop = 0, const int (*pSides)[4] = NULL,
                           sm_VSer<double> (*pC)[2] = NULL, sm_VSer<double> (*pN)[2] = NULL, SmBoolean bBefore = FALSE)
{
  rStation.dS = dS;
  rStation.dW = rStation.dWds = rStation.dWdsL = std::numeric_limits<double>::quiet_NaN();
  rStation.lSteps = 0;
  rStation.bSided = FALSE;
  for (int q=0; q<4; q++) { rStation.aiBefore[q] = 0; }
  const int lK = tr.s + tr.e + 1;
  sm_Partials<double> asP[2];
  sm_GnSystem<double> sSys;
  sSys.apP[0] = &asP[0];
  sSys.apP[1] = &asP[1];
  sSys.adR[0] = crR[0];
  sSys.adR[1] = crR[1];
  sSys.adSide[0] = crSide[0];
  sSys.adSide[1] = crSide[1];
  sSys.bHold = FALSE;
  sSys.pR = pR;
  SER(sm_LeaderSeries(crLeader, dS, tr, sSys.sLeader, sSys.sTangent, bBefore));

  // the feet at tau = 1, by Newton (the leader at sigma = 0 only)
  double x[4] = { crSeed[0], crSeed[1], crSeed[2], crSeed[3] };
  const sm_Trunc sPoint = { 0, 0 };
  SmBoolean bConverged = FALSE;
  for (int it=0; it<50 && !bConverged; it++)
    {
      rStation.lSteps++;
      for (int f=0; f<2; f++)
        { SER(sm_SurfacePartials(*apFace[f], SmPoint2d(x[2*f], x[2*f+1]), 2, asP[f], pSides != NULL ? (*pSides)[2*f] : 0, pSides != NULL ? (*pSides)[2*f+1] : 0)); }
      sm_Ser<double> sX[4];
      for (int q=0; q<4; q++) { sX[q] = sm_SerConst<double>(x[q]); }
      sm_Ser<double> sG[4];
      sm_VSer<double> sF[2];
      SER(sm_GnEvaluate(sSys, sX, sPoint, sG, sF));
      double J[4][4], b[4], dx[4];
      SER(sm_GnJacobian(sSys, sX, J));
      double dRes = 0.0;
      for (int i=0; i<4; i++) { b[i] = -sG[i].c[0][0]; dRes = smos_Max(dRes, smos_Fabs(b[i])); }
      SER(sm_Solve4(J, b, dx));
      double dStep = 0.0;
      for (int q=0; q<4; q++)
        {
          x[q] += dx[q];
          dStep = smos_Max(dStep, smos_Fabs(dx[q]) / (1.0 + smos_Fabs(x[q])));
        }
      // converged: the residual within the tolerance and the step at rounding
      if (dRes <= dTol && dStep <= 1.0e-12)
        { bConverged = TRUE; }
    }
  if (!bConverged)
    { SER_MSG(SM_ERR, _T("SmGnSectionFSG: the feet of a station did not converge")); }

  // the series, order by order: each coefficient of G is J x[i][k] plus the lower ones' (to order
  // lTop - 1 in eps with lTop, sm_GnTopJets)
  for (int f=0; f<2; f++)
    { SER(sm_SurfacePartials(*apFace[f], SmPoint2d(x[2*f], x[2*f+1]), lK, asP[f], pSides != NULL ? (*pSides)[2*f] : 0, pSides != NULL ? (*pSides)[2*f+1] : 0)); }
  for (int q=0; q<4; q++) { rStation.asX[q] = sm_SerConst<double>(x[q]); }
  double J[4][4];
  SER(sm_GnJacobian(sSys, rStation.asX, J));
  // (the coefficients of one total degree do not enter each other's: one evaluation solves them all)
  for (int lTot=1; lTot<=tr.s+tr.e; lTot++)
    {
    sm_Ser<double> sG[4];
    sm_VSer<double> sF[2];
    SER(sm_GnEvaluate(sSys, rStation.asX, tr, sG, sF));
    for (int i=0; i<=tr.s && i<=lTot; i++)
      {
        const int k = lTot - i;
        if (k > tr.e || (lTop >= 1 && k >= lTop))
          { continue; }
        double b[4], y[4];
        for (int j=0; j<4; j++) { b[j] = -sG[j].c[i][k]; }
        SER(sm_Solve4(J, b, y));
        for (int q=0; q<4; q++) { rStation.asX[q].c[i][k] = y[q]; }
      }
    }
  sm_Ser<double> sG[4];
  sm_VSer<double> asN[2];
  SER(sm_GnEvaluate(sSys, rStation.asX, tr, sG, rStation.asF, &rStation.sV, &asN));
  if (pC != NULL) { for (int f=0; f<2; f++) { (*pC)[f] = rStation.asF[f]; } }
  if (pN != NULL) { for (int f=0; f<2; f++) { (*pN)[f] = asN[f]; } }
  sm_GnTopJets(lTop, tr, asN, rStation.asF);
  return SM_SUCCESS;
}

/*******************************************************************//**
   Truncated Taylor series in one variable x, to order SM_GN_TAY: c[i] the
   coefficient of x^i, the i-th derivative over i!. Over intervals each
   coefficient holds that derivative at every point the constant term's
   interval stands for: the operations are the derivatives' own rules
   (Leibniz's, the chain rule's recurrences) applied at each point, so they
   hold pointwise (interval automatic differentiation). Every constant they
   take is an S, so that in intervals none is rounded unaccounted for.
***********************************************************************/
const int SM_GN_TAY = 4;

template <class S> struct sm_Tay
{
  S c[SM_GN_TAY + 1];
  sm_Tay() { for (int i=0; i<=SM_GN_TAY; i++) { c[i] = S(0.0); } }
  sm_Tay(double d) { for (int i=0; i<=SM_GN_TAY; i++) { c[i] = S(0.0); } c[0] = S(d); }
};

// x0 + x1 x
template <class S> sm_Tay<S> sm_TayOf(const S & x0, const S & x1)
{
  sm_Tay<S> r;
  r.c[0] = x0;
  r.c[1] = x1;
  return r;
}

template <class S> sm_Tay<S> operator+(const sm_Tay<S> & a, const sm_Tay<S> & b)
{
  sm_Tay<S> r;
  for (int i=0; i<=SM_GN_TAY; i++) { r.c[i] = a.c[i] + b.c[i]; }
  return r;
}

template <class S> sm_Tay<S> operator-(const sm_Tay<S> & a, const sm_Tay<S> & b)
{
  sm_Tay<S> r;
  for (int i=0; i<=SM_GN_TAY; i++) { r.c[i] = a.c[i] - b.c[i]; }
  return r;
}

template <class S> sm_Tay<S> operator*(const sm_Tay<S> & a, const sm_Tay<S> & b)
{
  sm_Tay<S> r;
  for (int i=0; i<=SM_GN_TAY; i++)
    {
      S sSum(0.0);
      for (int k=0; k<=i; k++) { sSum = sSum + a.c[k] * b.c[i-k]; }
      r.c[i] = sSum;
    }
  return r;
}

// a / b: q_i = (a_i - sum_{k=1..i} b_k q_(i-k)) / b_0 (an interval b_0 holding 0 gives the whole line)
template <class S> sm_Tay<S> operator/(const sm_Tay<S> & a, const sm_Tay<S> & b)
{
  sm_Tay<S> r;
  for (int i=0; i<=SM_GN_TAY; i++)
    {
      S sSum = a.c[i];
      for (int k=1; k<=i; k++) { sSum = sSum - b.c[k] * r.c[i-k]; }
      r.c[i] = sSum / b.c[0];
    }
  return r;
}

// The constant terms' functions, rounded outward in intervals: the C library's acos, sin and cos are
// within an ulp or two of the exact values, inside the margin of 1e-15 relative (four ulps) kept here.
inline sm_Ivl sm_Wide(double a, double b)
{
  if (std::isnan(a) || std::isnan(b))
    { return sm_Whole(); }
  return sm_Ivl(a - (std::fabs(a) * 1.0e-15 + 1.0e-300), b + (std::fabs(b) * 1.0e-15 + 1.0e-300));
}
inline double sm_Acos(double x) { return smos_ArcCosine(smos_Max(-1.0, smos_Min(1.0, x))); }
inline sm_Ivl sm_Acos(const sm_Ivl & x)
{
  if (!(x.lo >= -1.0 && x.hi <= 1.0))
    { return sm_Whole(); }
  return sm_Wide(smos_ArcCosine(x.hi), smos_ArcCosine(x.lo));
}
inline void sm_SinCos(double u, double & rS, double & rC) { rS = smos_Sine(u); rC = smos_Cosine(u); }
// on [0, pi/2] (the double below pi/2 its top), where sin rises and cos falls
inline void sm_SinCos(const sm_Ivl & u, sm_Ivl & rS, sm_Ivl & rC)
{
  if (!(u.lo >= 0.0 && u.hi <= 1.5707963267948966))
    {
      rS = rC = sm_Whole();
      return;
    }
  rS = sm_Wide(smos_Sine(u.lo), smos_Sine(u.hi));
  rC = sm_Wide(smos_Cosine(u.hi), smos_Cosine(u.lo));
}

// sqrt(y): r_0 = sqrt(y_0), r_i = (y_i - sum_{k=1..i-1} r_k r_(i-k)) / (2 r_0)
template <class S> sm_Tay<S> sm_TaySqrt(const sm_Tay<S> & y)
{
  sm_Tay<S> r;
  r.c[0] = sm_Sqrt(y.c[0]);
  for (int i=1; i<=SM_GN_TAY; i++)
    {
      S sSum = y.c[i];
      for (int k=1; k<i; k++) { sSum = sSum - r.c[k] * r.c[i-k]; }
      r.c[i] = sSum / (S(2.0) * r.c[0]);
    }
  return r;
}

// sin u and cos u: (sin u)' = cos u u', (cos u)' = -sin u u', coefficient by coefficient
template <class S> void sm_TaySinCos(const sm_Tay<S> & u, sm_Tay<S> & rS, sm_Tay<S> & rC)
{
  sm_SinCos(u.c[0], rS.c[0], rC.c[0]);
  for (int i=1; i<=SM_GN_TAY; i++)
    {
      S a(0.0), b(0.0);
      for (int k=1; k<=i; k++)
        {
          a = a + S(static_cast<double>(k)) * u.c[k] * rC.c[i-k];
          b = b + S(static_cast<double>(k)) * u.c[k] * rS.c[i-k];
        }
      rS.c[i] = a / S(static_cast<double>(i));
      rC.c[i] = (S(0.0) - b) / S(static_cast<double>(i));
    }
}

template <class Z> Z sm_PowN(const Z & a, ULONG n)
{
  Z r(1.0);
  for (ULONG m=0; m<n; m++) { r = r * a; }
  return r;
}

/*******************************************************************//**
   The area integral I(w) of the order-n section between two planes
   (EqualAreaThumbweight): the integral over t in [0, 1] of a b' - a' b,
   a = f (1 - p1), b = (1 - f) p2, which is a^2 (b/a)', that is
     (1 - w)^2 t^(n+1) (1 - t)^(n+1) ((n + 2) d1 d2 + t (1 - t) w (2 - w))
       / (D^2 d1^2 d2^2),
   D = (1 - t)^(n+1) + t^(n+1), d1 = 1 - w + w t, d2 = 1 - w t; by composite
   Simpson on SM_GN_AREA_N intervals. That sum is the I whose root the equal
   area thumbweight is: a rational function of w, smooth for w < 1, so its
   w-derivatives are its terms' (Z a series in w for them; an interval, or
   a series over one, for enclosures).
***********************************************************************/
const int SM_GN_AREA_N = 2048;

template <class Z> Z sm_GnAreaIntegrand(ULONG n, double dT, const Z & w)
{
  const Z one(1.0), t(dT);
  const Z u = one - t;
  const Z tn = sm_PowN(t, n + 1), un = sm_PowN(u, n + 1);
  const Z D = un + tn;
  const Z d1 = one - w + w * t, d2 = one - w * t, v = one - w;
  const Z num = v * v * tn * un * (Z(static_cast<double>(n + 2)) * d1 * d2 + t * u * w * (Z(2.0) - w));
  return num / (D * D * d1 * d1 * d2 * d2);
}

template <class Z> Z sm_GnAreaSum(ULONG n, const Z & w)
{
  Z sSum(0.0);
  for (int i=0; i<=SM_GN_AREA_N; i++)
    {
      const double c = (i == 0 || i == SM_GN_AREA_N) ? 1.0 : ((i % 2) ? 4.0 : 2.0);
      sSum = sSum + Z(c) * sm_GnAreaIntegrand(n, static_cast<double>(i) / SM_GN_AREA_N, w);
    }
  return sSum / Z(3.0 * SM_GN_AREA_N);
}

/*******************************************************************//**
   The equal area thumbweight along s. x is the cosine of the angle 2u
   between the ball's radii to its contacts, x = (F1 - V).(F2 - V) / (r1 r2)
   (V the ball's centre); the wedge's angle is alpha = pi - 2u. The ball's
   circle takes off what the section does where I(w) equals the target
   (sin u - u cos u) / sin^3 u (EqualAreaThumbweight), so w(s) solves
   I(w(s)) = target(u(s)), and its s-derivatives follow from the target's
   and I's w-derivatives order by order (the implicit function's).
***********************************************************************/
// x as a series in s from the feet's and V's (their s-coefficients to lOrder)
template <class S> sm_Tay<S> sm_GnWedgeCos(const sm_VSer<S> (&crF)[2], const sm_VSer<S> & crV, double dRR, int lOrder)
{
  sm_Tay<S> x;
  for (int c=0; c<3; c++)
    {
      sm_Tay<S> a[2];
      for (int f=0; f<2; f++)
        for (int i=0; i<=lOrder && i<=SM_GN_TAY; i++)
          a[f].c[i] = crF[f].a[c].c[i][0] - crV.a[c].c[i][0];
      x = x + a[0] * a[1];
    }
  const S sInv = S(1.0) / S(dRR);
  for (int i=0; i<=SM_GN_TAY; i++) { x.c[i] = x.c[i] * sInv; }
  return x;
}

// u = acos(x) / 2: u' = -x' / (2 sqrt(1 - x^2))
template <class S> sm_Tay<S> sm_GnHalfAngle(const sm_Tay<S> & x)
{
  const sm_Tay<S> q = sm_TaySqrt(sm_Tay<S>(1.0) - x * x);
  sm_Tay<S> dx;
  for (int i=0; i<SM_GN_TAY; i++) { dx.c[i] = x.c[i+1] * S(static_cast<double>(i + 1)); }
  const sm_Tay<S> v = dx / q;
  sm_Tay<S> u;
  u.c[0] = sm_Acos(x.c[0]) * S(0.5);
  for (int i=1; i<=SM_GN_TAY; i++) { u.c[i] = (S(0.0) - v.c[i-1]) / S(2.0 * i); }
  return u;
}

template <class S> sm_Tay<S> sm_GnTarget(const sm_Tay<S> & u)
{
  sm_Tay<S> s, c;
  sm_TaySinCos(u, s, c);
  return (s - u * c) / (s * s * s);
}

// w(s) from I(w(s)) = Tg(s), w(s) in w0: with dw = w - w(s), Tg_i is the coefficient of sigma^i in
// sum_m d_m dw^m, d_m = I^(m)(w(s)) / m! (held by crD over w0), in which w_i enters only through
// d_1 w_i.
template <class S> sm_Tay<S> sm_GnInvert(const sm_Tay<S> & crTg, const sm_Tay<S> & crD, const S & w0)
{
  sm_Tay<S> dw;
  for (int i=1; i<=SM_GN_TAY; i++)
    {
      S sP(0.0);
      sm_Tay<S> sPow = dw;
      for (int m=2; m<=i; m++)
        {
          sPow = sPow * dw;
          sP = sP + crD.c[m] * sPow.c[i];
        }
      dw.c[i] = (crTg.c[i] - sP) / crD.c[1];
    }
  dw.c[0] = w0;
  return dw;
}

/*******************************************************************//**
   Bernstein polynomials: coefficients of degree size() - 1, in any ring Z.
***********************************************************************/
typedef std::vector<double> sm_Bern;

// C(p, i), exact for the degrees here
inline double sm_Binom(size_t p, size_t i)
{
  double c = 1.0;
  for (size_t m=1; m<=i; m++) { c = c * static_cast<double>(p - i + m) / static_cast<double>(m); }
  return smos_Floor(c + 0.5);
}

template <class Z> std::vector<Z> sm_BernMul(const std::vector<Z> & a, const std::vector<Z> & b)
{
  const size_t p = a.size() - 1, q = b.size() - 1;
  std::vector<Z> r(p + q + 1, Z(0.0));
  for (size_t i=0; i<=p; i++)
    for (size_t j=0; j<=q; j++)
      {
        // C(p, i) C(q, j) / C(p + q, i + j)
        const Z c = Z(sm_Binom(p, i) * sm_Binom(q, j)) / Z(sm_Binom(p + q, i + j));
        r[i+j] = r[i+j] + c * a[i] * b[j];
      }
  return r;
}

template <class Z> std::vector<Z> sm_BernPow(const std::vector<Z> & a, ULONG n)
{
  std::vector<Z> r(1, Z(1.0));
  for (ULONG m=0; m<n; m++) { r = sm_BernMul(r, a); }
  return r;
}

/*******************************************************************//**
   The section of order n and thumbweight w over the common denominator:
   its weights W_j and, for each jet c_ik, the coefficient g_ik_j of c_ik in
   W_j P_j (j = 0..3n+1).
***********************************************************************/
template <class Z> struct sm_GnSectionFormT
{
  std::vector<Z>                W;
  std::vector< std::vector<Z> > aG[2]; // aG[i][k]
};
typedef sm_GnSectionFormT<double> sm_GnSectionForm;

template <class Z> sm_GnSectionFormT<Z> sm_GnFormT(ULONG n, const Z & w)
{
  sm_GnSectionFormT<Z> sForm;
  std::vector<Z> sA(n + 2, Z(0.0)), sB(n + 2, Z(0.0));
  sA[0] = Z(1.0);           // (1 - t)^(n+1)
  sB[n + 1] = Z(1.0);       // t^(n+1)
  const std::vector<Z> sOneMinusWt = { Z(1.0), Z(1.0) - w };  // 1 - w t
  const std::vector<Z> sOneMinusWpWt = { Z(1.0) - w, Z(1.0) }; // 1 - w + w t
  const std::vector<Z> sMinusT = { Z(0.0), Z(-1.0) };          // -t
  const std::vector<Z> sTMinus1 = { Z(-1.0), Z(0.0) };         // t - 1
  std::vector<Z> sAB(n + 2, Z(0.0));
  for (size_t j=0; j<sAB.size(); j++) { sAB[j] = sA[j] + sB[j]; }
  sForm.W = sm_BernMul(sm_BernMul(sAB, sm_BernPow(sOneMinusWpWt, n)), sm_BernPow(sOneMinusWt, n));
  for (ULONG k=0; k<=n; k++)
    {
      sForm.aG[0].push_back(sm_BernMul(sm_BernMul(sA, sm_BernPow(sOneMinusWt, n)),
                                       sm_BernMul(sm_BernPow(sMinusT, k), sm_BernPow(sOneMinusWpWt, n - k))));
      sForm.aG[1].push_back(sm_BernMul(sm_BernMul(sB, sm_BernPow(sOneMinusWpWt, n)),
                                       sm_BernMul(sm_BernPow(sTMinus1, k), sm_BernPow(sOneMinusWt, n - k))));
    }
  return sForm;
}

inline sm_GnSectionForm sm_GnForm(ULONG n, double w) { return sm_GnFormT<double>(n, w); }

// The form over W_0 (the weight at the first rail, (1 - w)^n, as at the second): the same section,
// its rows a common factor apart, with the weight 1 at the first rail whatever w is, so that rows
// interpolated along s keep that rail the jets' own curve.
template <class Z> sm_GnSectionFormT<Z> sm_GnFormNormalized(ULONG n, const Z & w)
{
  sm_GnSectionFormT<Z> sForm = sm_GnFormT<Z>(n, w);
  const Z sW0 = sForm.W[0];
  for (size_t j=0; j<sForm.W.size(); j++)
    {
      sForm.W[j] = sForm.W[j] / sW0;
      for (int i=0; i<2; i++)
        for (size_t k=0; k<sForm.aG[i].size(); k++)
          sForm.aG[i][k][j] = sForm.aG[i][k][j] / sW0;
    }
  return sForm;
}

/*******************************************************************//**
   A piecewise cubic Hermite curve on the stations: its control points
   (3 (N - 1) + 1), from values and derivatives in s at each.
***********************************************************************/
void sm_HermitePoints(const std::vector<double> & crS, const std::vector<SmVector3d> & crV, const std::vector<SmVector3d> & crD,
                      std::vector<SmVector3d> & rCPs)
{
  rCPs.clear();
  for (size_t j=0; j+1<crS.size(); j++)
    {
      const double h = crS[j+1] - crS[j];
      if (j == 0) { rCPs.push_back(crV[j]); }
      rCPs.push_back(crV[j] + crD[j] * (h / 3.0));
      rCPs.push_back(crV[j+1] - crD[j+1] * (h / 3.0));
      rCPs.push_back(crV[j+1]);
    }
}

// The same with each station's derivative as the span ending there takes it (crDIn) and as the
// span starting there does (crD): a station on a face's knot line has one of each.
void sm_HermitePoints(const std::vector<double> & crS, const std::vector<SmVector3d> & crV, const std::vector<SmVector3d> & crD,
                      const std::vector<SmVector3d> & crDIn, std::vector<SmVector3d> & rCPs)
{
  rCPs.clear();
  for (size_t j=0; j+1<crS.size(); j++)
    {
      const double h = crS[j+1] - crS[j];
      if (j == 0) { rCPs.push_back(crV[j]); }
      rCPs.push_back(crV[j] + crD[j] * (h / 3.0));
      rCPs.push_back(crV[j+1] - crDIn[j+1] * (h / 3.0));
      rCPs.push_back(crV[j+1]);
    }
}

// The knots of those curves: the stations, the inner ones three times (piecewise Bezier).
void sm_HermiteKnots(const std::vector<double> & crS, SmTArray<double> & rKnots, SmTArray<ULONG> & rMult)
{
  rKnots.ReSet();
  rMult.ReSet();
  for (size_t j=0; j<crS.size(); j++)
    {
      rKnots.Add(crS[j]);
      rMult.Add((j == 0 || j + 1 == crS.size()) ? 4 : 3);
    }
}


/*******************************************************************//**
   The certificate. For the rational section the blend's distance from its
   definition (the section of the true jets c_ik(s) at every s) is at most
   max_i sum_k |c_ik - C_ik|, C_ik the jets' Hermite curves: f is in [0, 1]
   and |u - 1| <= 1 in T_i. On a span of length h a cubic Hermite
   interpolant of a function's values and derivatives at both ends is off it
   by at most h^4 / 384 times the largest fourth derivative (the Hermite
   remainder). That bound comes from rounded interval arithmetic: the feet
   over the span are enclosed by Krawczyk's operator, which proves a unique
   root of G in its box for every s in the span when it maps the box into
   its interior and the preconditioned Jacobian's norm is below one; the
   feet's Taylor coefficients, to order 4 in s and n in tau, are then
   enclosed order by order through the same preconditioned system.
***********************************************************************/

// A NURBS curve's data (the leader), as sm_Nurbs.
struct sm_NurbsCurve
{
  int                 p, n;
  std::vector<double> U;
  std::vector<double> Pw; // (x w, y w, z w, w)
};

SmStatus sm_NurbsCurveOf(const SmBSplineCurve & crCurve, sm_NurbsCurve & rN)
{
  SmTArray<double> sKnots;
  SmTArray<ULONG> sMult;
  SER(crCurve.GetKnots(sKnots, &sMult));
  SmTArray<SmPoint3d> sCPs;
  SmTArray<double> sW;
  SER(crCurve.GetControlPolygon(sCPs, sW));
  rN.p = static_cast<int>(crCurve.GetDegree());
  rN.n = static_cast<int>(sCPs.GetSize());
  rN.U.clear();
  for (ULONG i=0; i<sKnots.GetSize(); i++)
    for (ULONG m=0; m<sMult[i]; m++)
      rN.U.push_back(sKnots[i]);
  if (rN.p < 1 || static_cast<int>(rN.U.size()) != rN.n + rN.p + 1)
    { SER_MSG(SM_ERR, _T("SmGnSectionFSG: the centre line's NURBS data do not agree")); }
  rN.Pw.resize(static_cast<size_t>(4 * rN.n));
  for (int i=0; i<rN.n; i++)
    {
      const ULONG lI = static_cast<ULONG>(i);
      const double w = (sW.GetSize() == sCPs.GetSize()) ? sW[lI] : 1.0;
      if (!(w > 0.0))
        { SER(SM_ERR); }
      rN.Pw[static_cast<size_t>(4*i)]   = sCPs[lI].x * w;
      rN.Pw[static_cast<size_t>(4*i+1)] = sCPs[lI].y * w;
      rN.Pw[static_cast<size_t>(4*i+2)] = sCPs[lI].z * w;
      rN.Pw[static_cast<size_t>(4*i+3)] = w;
    }
  return SM_SUCCESS;
}

// The derivatives D[0..lD] of the curve's polynomial piece lPiece (the span of the knots it starts
// at; -1: each piece S covers) for s in S, enclosed.
SmStatus sm_CurveDerivsOver(const sm_NurbsCurve & crN, const sm_Ivl & crS, int lD, std::vector< sm_V3<sm_Ivl> > & rD, int lPiece = -1)
{
  const int k0 = (lPiece >= 0) ? lPiece : sm_Span(crN.U, crN.p, crN.n, crS.lo);
  const int k1 = (lPiece >= 0) ? lPiece : sm_Span(crN.U, crN.p, crN.n, crS.hi);
  SmBoolean bFirst = TRUE;
  rD.assign(static_cast<size_t>(lD + 1), sm_V3<sm_Ivl>());
  for (int k=k0; k<=k1; k++)
    {
      if (!(crN.U[static_cast<size_t>(k)] < crN.U[static_cast<size_t>(k + 1)]))
        { continue; } // a knot repeated: no span
      const double a = (k == crN.p || lPiece >= 0) ? crS.lo : std::max(crS.lo, crN.U[static_cast<size_t>(k)]);
      const double b = (k == crN.n - 1 || lPiece >= 0) ? crS.hi : std::min(crS.hi, crN.U[static_cast<size_t>(k + 1)]);
      // an interval that ends on a knot is the span it lies in, not also the next one's end point
      if (a > b || (a == b && crS.lo < crS.hi))
        { continue; }
      // the homogeneous curve is a polynomial of degree p on the span: its derivatives over [a, b]
      // are its Taylor polynomials at the middle m, from exact point derivatives (rounding only)
      const double m = 0.5 * (a + b);
      std::vector< std::vector<sm_Ivl> > N;
      sm_DersBasis(crN.U, k, crN.p, sm_Ivl(m), crN.p, N);
      std::vector< sm_V3<sm_Ivl> > Am(static_cast<size_t>(crN.p + 1));
      std::vector<sm_Ivl> Wm(static_cast<size_t>(crN.p + 1));
      for (int d=0; d<=crN.p; d++)
        {
          sm_Ivl x(0.0), y(0.0), z(0.0), w(0.0);
          for (int j=0; j<=crN.p; j++)
            {
              const size_t lI = static_cast<size_t>(4 * (k - crN.p + j));
              x = x + N[d][j] * sm_Ivl(crN.Pw[lI]);
              y = y + N[d][j] * sm_Ivl(crN.Pw[lI+1]);
              z = z + N[d][j] * sm_Ivl(crN.Pw[lI+2]);
              w = w + N[d][j] * sm_Ivl(crN.Pw[lI+3]);
            }
          Am[d].x = x; Am[d].y = y; Am[d].z = z; Wm[d] = w;
        }
      const sm_Ivl dt = sm_Ivl(a, b) - sm_Ivl(m);
      std::vector< sm_V3<sm_Ivl> > A(static_cast<size_t>(lD + 1));
      std::vector<sm_Ivl> W(static_cast<size_t>(lD + 1));
      for (int d=0; d<=lD; d++)
        {
          sm_Ivl x(0.0), y(0.0), z(0.0), w(0.0), pw(1.0);
          double f = 1.0;
          for (int j=d; j<=crN.p; j++)
            {
              if (j > d) { pw = pw * dt; f *= (j - d); }
              const sm_Ivl c = pw * sm_Ivl(1.0 / f);
              x = x + Am[j].x * c; y = y + Am[j].y * c; z = z + Am[j].z * c; w = w + Wm[j] * c;
            }
          A[d].x = x; A[d].y = y; A[d].z = z; W[d] = w;
        }
      std::vector< sm_V3<sm_Ivl> > C(static_cast<size_t>(lD + 1));
      for (int d=0; d<=lD; d++)
        {
          sm_V3<sm_Ivl> v = A[d];
          double c = 1.0;
          for (int i=1; i<=d; i++)
            {
              c = c * (d - i + 1) / i;
              v.x = v.x - sm_Ivl(c) * W[i] * C[d-i].x;
              v.y = v.y - sm_Ivl(c) * W[i] * C[d-i].y;
              v.z = v.z - sm_Ivl(c) * W[i] * C[d-i].z;
            }
          C[d].x = v.x / W[0]; C[d].y = v.y / W[0]; C[d].z = v.z / W[0];
        }
      for (int d=0; d<=lD; d++)
        {
          if (bFirst) { rD[d] = C[d]; }
          else
            {
              rD[d].x = sm_Hull(rD[d].x, C[d].x);
              rD[d].y = sm_Hull(rD[d].y, C[d].y);
              rD[d].z = sm_Hull(rD[d].z, C[d].z);
            }
        }
      bFirst = FALSE;
    }
  if (bFirst)
    { SER(SM_ERR); }
  return SM_SUCCESS;
}

// The partials of one polynomial piece (span ku, kv) over [ua, ub] x [va, vb], enclosed: the
// homogeneous surface's Taylor polynomials at the box's middle, from exact point derivatives (it is
// of degree p, q there), then the rational quotient's (the NURBS book's A4.4) in intervals.
SmStatus sm_NurbsPartialsOverPiece(const sm_Nurbs & crN, int ku, int kv, double ua, double ub, double va, double vb, int lK,
                                   sm_Partials<sm_Ivl> & rP)
{
  const double um = 0.5 * (ua + ub), vm = 0.5 * (va + vb);
  std::vector< std::vector<sm_Ivl> > Nu, Nv;
  sm_DersBasis(crN.U, ku, crN.p, sm_Ivl(um), crN.p, Nu);
  sm_DersBasis(crN.V, kv, crN.q, sm_Ivl(vm), crN.q, Nv);
  const int lA = crN.p, lB = crN.q;
  // homogeneous point derivatives Hm[a][b] (x, y, z, w) at the middle
  std::vector< std::vector< std::array<sm_Ivl, 4> > > Hm(static_cast<size_t>(lA + 1), std::vector< std::array<sm_Ivl, 4> >(static_cast<size_t>(lB + 1)));
  for (int a=0; a<=lA; a++)
    for (int b=0; b<=lB; b++)
      {
        std::array<sm_Ivl, 4> h = { sm_Ivl(0.0), sm_Ivl(0.0), sm_Ivl(0.0), sm_Ivl(0.0) };
        for (int i=0; i<=crN.p; i++)
          for (int j=0; j<=crN.q; j++)
            {
              const size_t lI = static_cast<size_t>(4 * ((ku - crN.p + i) * crN.nv + (kv - crN.q + j)));
              const sm_Ivl c = Nu[a][i] * Nv[b][j];
              for (int e=0; e<4; e++) { h[e] = h[e] + c * sm_Ivl(crN.Pw[lI + e]); }
            }
        Hm[a][b] = h;
      }
  const sm_Ivl du = sm_Ivl(ua, ub) - sm_Ivl(um), dv = sm_Ivl(va, vb) - sm_Ivl(vm);
  std::vector<sm_Ivl> pu(static_cast<size_t>(lA + 1)), pv(static_cast<size_t>(lB + 1));
  pu[0] = sm_Ivl(1.0); pv[0] = sm_Ivl(1.0);
  for (int i=1; i<=lA; i++) { pu[i] = pu[i-1] * du * sm_Ivl(1.0 / i); }
  for (int j=1; j<=lB; j++) { pv[j] = pv[j-1] * dv * sm_Ivl(1.0 / j); }
  std::vector< std::vector< sm_V3<sm_Ivl> > > A(static_cast<size_t>(lK + 1), std::vector< sm_V3<sm_Ivl> >(static_cast<size_t>(lK + 1)));
  std::vector< std::vector<sm_Ivl> > W(static_cast<size_t>(lK + 1), std::vector<sm_Ivl>(static_cast<size_t>(lK + 1), sm_Ivl(0.0)));
  for (int a=0; a<=lK; a++)
    for (int b=0; a+b<=lK; b++)
      {
        std::array<sm_Ivl, 4> h = { sm_Ivl(0.0), sm_Ivl(0.0), sm_Ivl(0.0), sm_Ivl(0.0) };
        for (int i=a; i<=lA; i++)
          for (int j=b; j<=lB; j++)
            {
              const sm_Ivl c = pu[static_cast<size_t>(i - a)] * pv[static_cast<size_t>(j - b)];
              for (int e=0; e<4; e++) { h[e] = h[e] + Hm[i][j][e] * c; }
            }
        A[a][b].x = h[0]; A[a][b].y = h[1]; A[a][b].z = h[2]; W[a][b] = h[3];
      }
  if (!(W[0][0].lo > 0.0))
    { SER_MSG(SM_ERR, _T("SmGnSectionFSG: a face's weights are not positive over a box")); }
  for (int a=0; a<=lK; a++)
    for (int b=0; a+b<=lK; b++)
      {
        sm_V3<sm_Ivl> v3 = A[a][b];
        double cka = 1.0;
        for (int i=0; i<=a; i++)
          {
            if (i > 0) { cka = cka * (a - i + 1) / i; }
            double cbj = 1.0;
            for (int j=0; j<=b; j++)
              {
                if (j > 0) { cbj = cbj * (b - j + 1) / j; }
                if (i == 0 && j == 0) { continue; }
                const sm_Ivl c = sm_Ivl(cka * cbj) * W[i][j];
                v3.x = v3.x - c * rP.P[a-i][b-j].x;
                v3.y = v3.y - c * rP.P[a-i][b-j].y;
                v3.z = v3.z - c * rP.P[a-i][b-j].z;
              }
          }
        rP.P[a][b].x = v3.x / W[0][0];
        rP.P[a][b].y = v3.y / W[0][0];
        rP.P[a][b].z = v3.z / W[0][0];
      }
  return SM_SUCCESS;
}

// The leader's point and unit tangent over S as series in sigma to tr.s, enclosed.
SmStatus sm_LeaderSeriesOver(const sm_NurbsCurve & crLeader, const sm_Ivl & crS, const sm_Trunc & tr, sm_VSer<sm_Ivl> & rL, sm_VSer<sm_Ivl> & rT)
{
  const sm_Trunc sSigma = { tr.s, 0 };
  std::vector< sm_V3<sm_Ivl> > sD;
  // the span lies in one piece of the leader (the stations hold its knots)
  SER(sm_CurveDerivsOver(crLeader, crS, tr.s + 1, sD, sm_Span(crLeader.U, crLeader.p, crLeader.n, sm_Mid(crS))));
  const sm_V3<sm_Ivl> sZero = { sm_Ivl(0.0), sm_Ivl(0.0), sm_Ivl(0.0) };
  rL = sm_VSerConst<sm_Ivl>(sZero);
  sm_VSer<sm_Ivl> sD1 = rL;
  double dFact = 1.0;
  for (int j=0; j<=tr.s; j++)
    {
      if (j > 0) { dFact *= j; }
      const sm_Ivl f(1.0 / dFact);
      const sm_Ivl acL[3] = { sD[j].x, sD[j].y, sD[j].z }, acD[3] = { sD[j+1].x, sD[j+1].y, sD[j+1].z };
      for (int c=0; c<3; c++)
        {
          rL.a[c].c[j][0] = acL[c] * f;
          sD1.a[c].c[j][0] = acD[c] * f;
        }
    }
  sm_Ser<sm_Ivl> sInv;
  SER(sm_SerInvSqrt(sm_VSerDot(sD1, sD1, sSigma), sSigma, sInv));
  rT = sm_VSerTimes(sD1, sInv, sSigma);
  return SM_SUCCESS;
}

/*******************************************************************//**
   A ball whose radius varies: its law r(s) over the leader's parameter (the
   x of a B-spline curve), as a series in sigma at a point or over a span.
***********************************************************************/
struct sm_GnLaw
{
  const SmBSplineCurve * pCurve;
  sm_NurbsCurve          sN;
};

// r at dS to order lOrder in sigma: rR, and rpR pointing at it (NULL without a law)
SmStatus sm_LawSeries(const sm_GnLaw * pLaw, double dS, int lOrder, sm_Ser<double> & rR, const sm_Ser<double> *& rpR,
                      SmBoolean bBefore = FALSE)
{
  rpR = NULL;
  if (pLaw == NULL)
    { return SM_SUCCESS; }
  rR = sm_SerConst<double>(0.0);
  ULONG lDerivs = static_cast<ULONG>(lOrder);
  if (!pLaw->pCurve->IsRational())
    { lDerivs = smos_Min(lDerivs, pLaw->pCurve->GetDegree()); }
  SmVector3d asD[SM_GN_MAXS + 2];
  SER(pLaw->pCurve->Evaluate(dS, lDerivs, bBefore ? FALSE : TRUE, asD, FALSE));
  double dFact = 1.0;
  for (ULONG i=0; i<=lDerivs; i++)
    {
      if (i > 0) { dFact *= static_cast<double>(i); }
      rR.c[i][0] = asD[i].x / dFact;
    }
  if (!(rR.c[0][0] > 0.0))
    { SER_MSG(SM_ERR, _T("SmGnSectionFSG: the radius law is not positive")); }
  rpR = &rR;
  return SM_SUCCESS;
}

// r over the span crS (in one piece of the law: the stations hold its knots), to order lOrder
SmStatus sm_LawSeriesOver(const sm_GnLaw & crLaw, const sm_Ivl & crS, int lOrder, sm_Ser<sm_Ivl> & rR)
{
  rR = sm_SerConst<sm_Ivl>(sm_Ivl(0.0));
  std::vector< sm_V3<sm_Ivl> > sD;
  SER(sm_CurveDerivsOver(crLaw.sN, crS, lOrder, sD, sm_Span(crLaw.sN.U, crLaw.sN.p, crLaw.sN.n, sm_Mid(crS))));
  double dFact = 1.0;
  for (int i=0; i<=lOrder; i++)
    {
      if (i > 0) { dFact *= i; }
      rR.c[i][0] = sD[static_cast<size_t>(i)].x * sm_Ivl(1.0 / dFact);
    }
  if (!sm_Positive(rR.c[0][0]))
    { SER(SM_ERR); }
  return SM_SUCCESS;
}

// The inverse of a 4 x 4 matrix (the preconditioner; any approximation will do).
SmStatus sm_Inverse4(const double (&crA)[4][4], double (&rY)[4][4])
{
  for (int j=0; j<4; j++)
    {
      double e[4] = { 0.0, 0.0, 0.0, 0.0 }, y[4];
      e[j] = 1.0;
      SER(sm_Solve4(crA, e, y));
      for (int i=0; i<4; i++) { rY[i][j] = y[i]; }
    }
  return SM_SUCCESS;
}

// In the certificate an evaluation that fails (an enclosure too wide for a normal or a weight to
// be proved non-zero) proves nothing: the span's bound stays HUGE_VAL, and the span is halved.
#define SM_GN_TRY(a) { if ((a) != SM_SUCCESS) { return SM_SUCCESS; } }

/*******************************************************************//**
   The feet's series over S = [crA.dS, crB.dS] (a point when they are one
   station), to order lOrderS in s and n in tau, enclosed: rF (and pV, the
   point V's), with rbProved FALSE where Krawczyk's test fails on S.
***********************************************************************/
SmStatus sm_GnEnclose(const sm_Nurbs * apFace[2], const sm_NurbsCurve & crLeader, const sm_GnStation & crA, const sm_GnStation & crB,
                      const double (&crR)[2], const double (&crSide)[2], int n, int lOrderS, sm_VSer<sm_Ivl> (&rF)[2], SmBoolean & rbProved,
                      sm_Ser<sm_Ivl> (*pX)[4] = NULL, double (*pDS)[2] = NULL, sm_VSer<sm_Ivl> * pV = NULL, const sm_GnLaw * pLaw = NULL)
{
  rbProved = FALSE;
  const sm_Ivl sS(crA.dS, crB.dS);
  const sm_Trunc sTr = { lOrderS, n };
  sm_GnSystem<sm_Ivl> sSys;
  sm_Partials<sm_Ivl> asP[2];
  sSys.apP[0] = &asP[0];
  sSys.apP[1] = &asP[1];
  for (int f=0; f<2; f++) { sSys.adR[f] = crR[f]; sSys.adSide[f] = crSide[f]; }
  sSys.bHold = FALSE;
  sSys.pR = NULL;
  sm_Ser<sm_Ivl> sLawR;
  if (pLaw != NULL)
    {
      SM_GN_TRY(sm_LawSeriesOver(*pLaw, sS, lOrderS, sLawR));
      sSys.pR = &sLawR;
    }
  SM_GN_TRY(sm_LeaderSeriesOver(crLeader, sS, sTr, sSys.sLeader, sSys.sTangent));

  // each face's polynomial piece over the span: the one its middle lies in. A station stands where a
  // foot crosses a knot line, so every foot over the span lies in that piece, and the piece's
  // polynomial, continued over the rest of a box, is the face wherever a foot can be: the root
  // Krawczyk's test proves unique in a box is the foot
  int aku[2], akv[2];
  for (int f=0; f<2; f++)
    {
      aku[f] = sm_Span(apFace[f]->U, apFace[f]->p, apFace[f]->nu, 0.5 * (crA.asX[2*f].c[0][0] + crB.asX[2*f].c[0][0]));
      akv[f] = sm_Span(apFace[f]->V, apFace[f]->q, apFace[f]->nv, 0.5 * (crA.asX[2*f+1].c[0][0] + crB.asX[2*f+1].c[0][0]));
    }

  // the feet's box: Krawczyk's operator from the two stations' hull
  double xm[4];
  sm_Ivl X[4];
  for (int c=0; c<4; c++)
    {
      const double a = crA.asX[c].c[0][0], b = crB.asX[c].c[0][0];
      xm[c] = 0.5 * (a + b);
      const double w = 0.5 * std::fabs(b - a) + 1.0e-9 * (1.0 + std::fabs(xm[c]));
      X[c] = sm_Ivl(xm[c] - 2.0 * w, xm[c] + 2.0 * w);
    }
  double Y[4][4];
  sm_Ivl C[4][4];
  double q = HUGE_VAL;
  SmBoolean bProved = FALSE;
  for (int it=0; it<8 && !bProved; it++)
    {
      // G at the box's centre for every s in the span
      for (int f=0; f<2; f++)
        { SM_GN_TRY(sm_NurbsPartialsOverPiece(*apFace[f], aku[f], akv[f], xm[2*f], xm[2*f], xm[2*f+1], xm[2*f+1], 2, asP[f])); }
      sm_Ser<sm_Ivl> sXc[4];
      for (int c=0; c<4; c++) { sXc[c] = sm_SerConst<sm_Ivl>(sm_Ivl(xm[c])); }
      const sm_Trunc sPoint = { 0, 0 };
      sm_Ser<sm_Ivl> sG[4];
      sm_VSer<sm_Ivl> sF[2];
      SM_GN_TRY(sm_GnEvaluate(sSys, sXc, sPoint, sG, sF));
      // the Jacobian over the box
      for (int f=0; f<2; f++)
        { SM_GN_TRY(sm_NurbsPartialsOverPiece(*apFace[f], aku[f], akv[f], X[2*f].lo, X[2*f].hi, X[2*f+1].lo, X[2*f+1].hi, 2, asP[f])); }
      sm_Ser<sm_Ivl> sXb[4];
      for (int c=0; c<4; c++) { sXb[c] = sm_SerConst<sm_Ivl>(X[c]); }
      sm_Ivl J[4][4];
      SM_GN_TRY(sm_GnJacobian(sSys, sXb, J));
      double Jm[4][4];
      for (int i=0; i<4; i++)
        for (int j=0; j<4; j++)
          Jm[i][j] = sm_Mid(J[i][j]);
      if (sm_Inverse4(Jm, Y) != SM_SUCCESS)
        { return SM_SUCCESS; }
      q = 0.0;
      for (int i=0; i<4; i++)
        {
          double dRow = 0.0;
          for (int j=0; j<4; j++)
            {
              sm_Ivl c((i == j) ? 1.0 : 0.0);
              for (int k=0; k<4; k++) { c = c - sm_Ivl(Y[i][k]) * J[k][j]; }
              C[i][j] = c;
              dRow = sm_Up(dRow + sm_Mag(c));
            }
          q = std::max(q, dRow);
        }
      sm_Ivl K[4];
      SmBoolean bInside = (q < 1.0) ? TRUE : FALSE;
      for (int i=0; i<4; i++)
        {
          sm_Ivl k(xm[i]);
          for (int j=0; j<4; j++) { k = k - sm_Ivl(Y[i][j]) * sG[j].c[0][0]; }
          for (int j=0; j<4; j++) { k = k + C[i][j] * (X[j] - sm_Ivl(xm[j])); }
          K[i] = k;
          if (!(k.lo > X[i].lo && k.hi < X[i].hi)) { bInside = FALSE; }
        }
      if (bInside)
        {
          for (int i=0; i<4; i++) { X[i] = sm_Ivl(std::max(X[i].lo, K[i].lo), std::min(X[i].hi, K[i].hi)); }
          bProved = TRUE;
        }
      else
        {
          if (!(q < 1.0))
            { return SM_SUCCESS; } // the span is too long for its Jacobian to be one preconditioner's
          // the next box from the contraction: K's radii are at most d + |C| w for a box of radii w
          // about the centre, d those of Y G, so the box of radii w = (I - |C|)^-1 d holds its K;
          // a half more, to clear rounding and the change of C with the box
          double d[4], w[4];
          for (int i=0; i<4; i++)
            {
              sm_Ivl v(0.0);
              for (int j=0; j<4; j++) { v = v + sm_Ivl(Y[i][j]) * sG[j].c[0][0]; }
              d[i] = sm_Mag(v - sm_Ivl(sm_Mid(v))) + std::fabs(sm_Mid(v));
              w[i] = d[i];
            }
          for (int k=0; k<200; k++)
            {
              double wn[4];
              for (int i=0; i<4; i++)
                {
                  wn[i] = d[i];
                  for (int j=0; j<4; j++) { wn[i] += sm_Mag(C[i][j]) * w[j]; }
                }
              for (int i=0; i<4; i++) { w[i] = wn[i]; }
            }
          for (int i=0; i<4; i++)
            {
              const double r = 1.5 * w[i] + 1.0e-12 * (1.0 + std::fabs(xm[i]));
              if (!std::isfinite(r))
                { return SM_SUCCESS; }
              X[i] = sm_Ivl(xm[i] - r, xm[i] + r);
            }
        }
    }
  if (!bProved)
    { return SM_SUCCESS; }

  // the feet's Taylor coefficients over the span, order by order: x = -Y R + C x, |x| <= rho
  const int lK = sTr.s + sTr.e + 1;
  for (int f=0; f<2; f++)
    { SM_GN_TRY(sm_NurbsPartialsOverPiece(*apFace[f], aku[f], akv[f], X[2*f].lo, X[2*f].hi, X[2*f+1].lo, X[2*f+1].hi, lK, asP[f])); }
  sm_Ser<sm_Ivl> sX[4];
  for (int c=0; c<4; c++) { sX[c] = sm_SerConst<sm_Ivl>(X[c]); }
  for (int lTot=1; lTot<=sTr.s+sTr.e; lTot++)
    {
    sm_Ser<sm_Ivl> sG[4];
    sm_VSer<sm_Ivl> sF[2];
    SM_GN_TRY(sm_GnEvaluate(sSys, sX, sTr, sG, sF));
    for (int i=0; i<=sTr.s && i<=lTot; i++)
      {
        const int k = lTot - i;
        if (k > sTr.e || (n >= 1 && k >= n))
          { continue; } // the section's jets of order n are the face curve's (sm_GnTopJets)
        sm_Ivl yr[4];
        double dRho = 0.0;
        for (int r=0; r<4; r++)
          {
            sm_Ivl v(0.0);
            for (int c=0; c<4; c++) { v = v - sm_Ivl(Y[r][c]) * sG[c].c[i][k]; }
            yr[r] = v;
            dRho = std::max(dRho, sm_Mag(v));
          }
        dRho = sm_Up(dRho / (1.0 - q));
        if (!std::isfinite(dRho))
          { return SM_SUCCESS; }
        for (int r=0; r<4; r++)
          {
            sm_Ivl v = yr[r];
            for (int c=0; c<4; c++) { v = v + C[r][c] * sm_Ivl(-dRho, dRho); }
            sX[r].c[i][k] = v;
          }
      }
    }
  sm_Ser<sm_Ivl> sG[4];
  sm_VSer<sm_Ivl> asN[2];
  SM_GN_TRY(sm_GnEvaluate(sSys, sX, sTr, sG, rF, pV, &asN));
  sm_GnTopJets(n, sTr, asN, rF);
  if (pX != NULL)
    { for (int c=0; c<4; c++) { (*pX)[c] = sX[c]; } }
  if (pDS != NULL)
    {
      // the faces' first partials over the box: |S_u| + |S_v| bounds |dS| per unit uv step (max norm)
      for (int f=0; f<2; f++)
        {
          const sm_V3<sm_Ivl> & crU = asP[f].P[1][0], & crV = asP[f].P[0][1];
          const double dU = smos_Sqrt(sm_Mag(crU.x) * sm_Mag(crU.x) + sm_Mag(crU.y) * sm_Mag(crU.y) + sm_Mag(crU.z) * sm_Mag(crU.z));
          const double dV = smos_Sqrt(sm_Mag(crV.x) * sm_Mag(crV.x) + sm_Mag(crV.y) * sm_Mag(crV.y) + sm_Mag(crV.z) * sm_Mag(crV.z));
          (*pDS)[f] = dU + dV;
        }
    }
  rbProved = TRUE;
  return SM_SUCCESS;
}

/*******************************************************************//**
   Enclosures of I's Taylor coefficients in w (sm_GnAreaSum) over intervals
   of w, on a grid of cells over an edge's thumbweights, each computed once:
   cell k's series holds I^(m)(w) / m! for every w in it, and the hull of
   the cells an interval meets holds them over that interval. An interval
   off the grid is taken on its own.
***********************************************************************/
struct sm_GnAreaCells
{
  ULONG                         lOrder;
  std::vector<double>           sE;    // the cells' ends
  std::vector<int>              sHave;
  std::vector< sm_Tay<sm_Ivl> > sI;

  void Init(ULONG n, double dLo, double dHi, int lCells)
  {
    lOrder = n;
    sE.clear();
    for (int k=0; k<=lCells; k++) { sE.push_back(dLo + (dHi - dLo) * k / lCells); }
    sHave.assign(static_cast<size_t>(lCells), 0);
    sI.assign(static_cast<size_t>(lCells), sm_Tay<sm_Ivl>());
  }

  sm_Tay<sm_Ivl> Over(const sm_Ivl & W)
  {
    const int lCells = static_cast<int>(sE.size()) - 1;
    if (lCells < 1 || !(W.lo >= sE.front() && W.hi <= sE.back()))
      { return sm_GnAreaSum(lOrder, sm_TayOf(W, sm_Ivl(1.0))); }
    int k0 = static_cast<int>(std::upper_bound(sE.begin(), sE.end(), W.lo) - sE.begin()) - 1;
    int k1 = static_cast<int>(std::lower_bound(sE.begin(), sE.end(), W.hi) - sE.begin()) - 1;
    k0 = std::max(0, std::min(lCells - 1, k0));
    k1 = std::max(k0, std::min(lCells - 1, k1));
    sm_Tay<sm_Ivl> sR;
    for (int k=k0; k<=k1; k++)
      {
        const size_t lK = static_cast<size_t>(k);
        if (!sHave[lK])
          {
            sI[lK] = sm_GnAreaSum(lOrder, sm_TayOf(sm_Ivl(sE[lK], sE[lK + 1]), sm_Ivl(1.0)));
            sHave[lK] = 1;
          }
        for (int m=0; m<=SM_GN_TAY; m++) { sR.c[m] = (k == k0) ? sI[lK].c[m] : sm_Hull(sR.c[m], sI[lK].c[m]); }
      }
    return sR;
  }
};

/*******************************************************************//**
   The certificate of a section whose thumbweight varies along s (the equal
   area one, resolved per station). The surface interpolates the section's
   homogeneous rows H_j(s) = sum_ik g_ikj(w(s)) c_ik(s) and W_j(w(s)) (the
   form over W_0, sm_GnFormNormalized), each a cubic Hermite curve of its
   values and s-derivatives at the stations. With E the Hermite remainder
   (E[q] = Hermite(q) - q, |E[q]| <= K max |q''''| on a span of length h,
   K = h^4 / 384), a fixed thumbweight wb (the span's first station's) and
   any fixed point O, the rows relative to O split as
     H_j - W_j O = sum_ik g_ikj(wb) c'_ik + sum_ik dg_ikj c'_ik,
   c'_i0 = c_i0 - O (sum_i g_i0j = W_j: the section is affine), c'_ik = c_ik
   for k > 0, dg = g(w(s)) - g(wb). The surface's point less the section's
   is sum_j B_j (E[H_j - W_j O] - E[W_j] (S - O)) / sum_j B_j (W_j + E[W_j]),
   B_j the Bernstein basis in t, and the first part of E[H_j - W_j O],
   summed with B_j, is W(t; wb) sum_ik phi_ik(t; wb) E[c_ik], phi_ik the
   section's functions of the jets, so at most W(t; wb) A, A = max_i sum_k
   |E[c_ik]| (the fixed section's bound). A ratio of sums of positive terms
   is at most their terms' largest ratio, so the distance is at most
     max_j (W_j(wb) A + K X_j + K Y_j R) / (min W_j - K Y_j),
   X_j and Y_j the largest |(sum_ik dg_ikj c'_ik)''''| and |W_j(w(s))''''|
   over the span, R the largest |S - O| (at most max_i sum_k |c'_ik|: the
   section's functions are f (-p1)^k and (1 - f) (p2 - 1)^k), min W_j the
   least W_j(w(s)). A constant thumbweight makes X_j = Y_j = 0 and the
   bound A.
   The fourth derivatives come from interval series in s over each piece of
   the span: x from the enclosed feet and V, the target through acos, sin
   and cos, w through I(w) = target order by order with I's w-derivatives
   over w's enclosure, and the rows through the form's ring operations.
   w's enclosure is the interval Newton operator's: N(W) = g - (I(g) -
   T) / I'(W) inside W proves, for every T the target takes on the piece,
   one root of I(w) = T in W, and it lies in N(W).
***********************************************************************/
struct sm_GnVarying
{
  ULONG            n;
  double           dRR;    // r1 r2
  sm_GnAreaCells * pCells;
};

struct sm_GnRowTerms
{
  std::vector<double> X, Y, Wlo;
  double              R;
};

// The coefficients 0..SM_GN_TAY of a series over a piece in their Taylor form at its middle m:
// c_i(s) = sum_{j < K - i} C(i + j, i) c_(i+j)(m) d^j + C(K, i) c_K(xi) d^(K-i), d = s - m, from
// the middle's coefficients acMid[0..K-1] and the piece's crTop at order K (the remainder's). The
// form keeps what a quantity's enclosure over the piece loses to its parts' motion along it (the
// difference of two feet that move together).
sm_Tay<sm_Ivl> sm_Centred(const sm_Ivl * acMid, const sm_Ivl & crTop, int K, const sm_Ivl & crD)
{
  sm_Tay<sm_Ivl> r;
  for (int i=0; i<=SM_GN_TAY; i++)
    {
      sm_Ivl sSum(0.0), sPow(1.0);
      for (int j=0; i+j<K; j++)
        {
          if (j > 0) { sPow = sPow * crD; }
          sSum = sSum + sm_Ivl(sm_Binom(static_cast<size_t>(i + j), static_cast<size_t>(i))) * acMid[i+j] * sPow;
        }
      sPow = sPow * crD;
      r.c[i] = sSum + sm_Ivl(sm_Binom(static_cast<size_t>(K), static_cast<size_t>(i))) * crTop * sPow;
    }
  return r;
}

// The row terms over one piece, folded into rT (largest X, Y, R; least W): rbProved FALSE where
// w's enclosure or a term is not proved. crF, crCentre: the feet's and V's series over the piece
// (to order K = 4 + SM_GN_TAYLOR in s); crFm, crCm: at its middle (to K - 1); crD: s less the
// middle over the piece.
SmStatus sm_GnRowTermsOver(const sm_GnVarying & crVary, const sm_VSer<sm_Ivl> (&crF)[2], const sm_VSer<sm_Ivl> & crCentre,
                           const sm_VSer<sm_Ivl> (&crFm)[2], const sm_VSer<sm_Ivl> & crCm, const sm_Ivl & crD, double dGuess,
                           const sm_GnSectionFormT<sm_Ivl> & crBar, const double (&crO)[3], sm_GnRowTerms & rT, SmBoolean & rbProved)
{
  rbProved = FALSE;
  const ULONG n = crVary.n;
  const int K = 4 + SM_GN_TAYLOR;
  // x = (F1 - V).(F2 - V) / (r1 r2): at the middle to order K - 1, over the piece at order K
  sm_Ivl axMid[K], sxTop(0.0);
  for (int i=0; i<K; i++) { axMid[i] = sm_Ivl(0.0); }
  for (int q=0; q<3; q++)
    for (int i=0; i<=K; i++)
      for (int l=0; l<=i; l++)
        {
          if (i < K)
            {
              axMid[i] = axMid[i] + (crFm[0].a[q].c[l][0] - crCm.a[q].c[l][0]) * (crFm[1].a[q].c[i-l][0] - crCm.a[q].c[i-l][0]);
            }
          else
            {
              sxTop = sxTop + (crF[0].a[q].c[l][0] - crCentre.a[q].c[l][0]) * (crF[1].a[q].c[i-l][0] - crCentre.a[q].c[i-l][0]);
            }
        }
  const sm_Ivl sInvRR = sm_Ivl(1.0) / sm_Ivl(crVary.dRR);
  for (int i=0; i<K; i++) { axMid[i] = axMid[i] * sInvRR; }
  sxTop = sxTop * sInvRR;
  const sm_Tay<sm_Ivl> u = sm_GnHalfAngle(sm_Centred(axMid, sxTop, K, crD));
  const sm_Tay<sm_Ivl> sTg = sm_GnTarget(u);
  const sm_Ivl sPhi = sm_GnAreaSum<sm_Ivl>(n, sm_Ivl(dGuess)) - sTg.c[0];
  sm_Ivl W0;
  SmBoolean bRoot = FALSE;
  double dRho = 1.0e-12 * (1.0 + std::fabs(dGuess));
  for (int it=0; it<40 && !bRoot; it++, dRho *= 4.0)
    {
      const sm_Ivl W(dGuess - dRho, dGuess + dRho);
      if (!(W.hi < 1.0))
        { break; }
      const sm_Ivl d1 = crVary.pCells->Over(W).c[1];
      if (!(d1.hi < 0.0))
        { continue; }
      const sm_Ivl N = sm_Ivl(dGuess) - sPhi / d1;
      if (N.lo > W.lo && N.hi < W.hi)
        {
          W0 = N;
          bRoot = TRUE;
        }
    }
  if (!bRoot)
    { return SM_SUCCESS; }
  const sm_Tay<sm_Ivl> sW = sm_GnInvert(sTg, crVary.pCells->Over(W0), W0);
  const sm_GnSectionFormT< sm_Tay<sm_Ivl> > sForm = sm_GnFormNormalized(n, sW);
  // the jets c'_ik as series in s over the piece
  sm_Tay<sm_Ivl> c[2][SM_GN_MAXE + 1][3];
  double dR = 0.0;
  for (int f=0; f<2; f++)
    {
      double dSum = 0.0;
      for (ULONG k=0; k<=n; k++)
        {
          double dSq = 0.0;
          for (int q=0; q<3; q++)
            {
              sm_Ivl acMid[K];
              for (int i=0; i<K; i++) { acMid[i] = crFm[f].a[q].c[i][k]; }
              c[f][k][q] = sm_Centred(acMid, crF[f].a[q].c[K][k], K, crD);
              if (k == 0)
                { c[f][k][q].c[0] = c[f][k][q].c[0] - sm_Ivl(crO[q]); }
              dSq += sm_Mag(c[f][k][q].c[0]) * sm_Mag(c[f][k][q].c[0]);
            }
          dSum += smos_Sqrt(dSq);
        }
      dR = std::max(dR, dSum);
    }
  const size_t lD = sForm.W.size();
  for (size_t j=0; j<lD; j++)
    {
      double dSq = 0.0;
      for (int q=0; q<3; q++)
        {
          // the coefficient of sigma^4 of sum_ik dg_ikj c'_ik
          sm_Ivl s4(0.0);
          for (int f=0; f<2; f++)
            for (ULONG k=0; k<=n; k++)
              {
                sm_Tay<sm_Ivl> dg = sForm.aG[f][k][j];
                dg.c[0] = dg.c[0] - crBar.aG[f][k][j];
                for (int i=0; i<=SM_GN_TAY; i++) { s4 = s4 + dg.c[i] * c[f][k][q].c[SM_GN_TAY - i]; }
              }
          dSq += sm_Mag(s4) * sm_Mag(s4);
        }
      const double dX = 24.0 * smos_Sqrt(dSq), dY = 24.0 * sm_Mag(sForm.W[j].c[SM_GN_TAY]), dWlo = sForm.W[j].c[0].lo;
          if (!(std::isfinite(dX) && std::isfinite(dY) && dWlo > 0.0))
        { return SM_SUCCESS; }
      rT.X[j] = std::max(rT.X[j], dX);
      rT.Y[j] = std::max(rT.Y[j], dY);
      rT.Wlo[j] = std::min(rT.Wlo[j], dWlo);
    }
  if (!std::isfinite(dR))
    { return SM_SUCCESS; }
  rT.R = std::max(rT.R, dR);
  rbProved = TRUE;
  return SM_SUCCESS;
}

/*******************************************************************//**
   The ball's circle as rows: the rational quadratic arc from the rail rho1
   to rho2 about the ball's centre m, its control points rho1, P1, rho2 and
   weights 1, w, 1, w = cos(theta / 2) (theta the angle between the radii
   b_i = rho_i - m), P1 = m + (mid - m) / w^2 where the rails' tangents meet;
   as homogeneous rows H_j = W_j (P_j - O) about a fixed point O. The
   series are in s (to SM_GN_TAY), over a ring S (double, or intervals).
***********************************************************************/
template <class S> SmBoolean sm_ArcRows(const sm_Tay<S> (&crRho1)[3], const sm_Tay<S> (&crRho2)[3], const sm_Tay<S> (&crM)[3],
                                        const double (&crO)[3], sm_Tay<S> (&rH)[3][3], sm_Tay<S> (&rW)[3])
{
  sm_Tay<S> b11, b22, b12;
  for (int q=0; q<3; q++)
    {
      const sm_Tay<S> b1 = crRho1[q] - crM[q], b2 = crRho2[q] - crM[q];
      b11 = b11 + b1 * b1;
      b22 = b22 + b2 * b2;
      b12 = b12 + b1 * b2;
    }
  if (!sm_Positive(b11.c[0]) || !sm_Positive(b22.c[0]))
    { return FALSE; }
  // cos^2(theta / 2) = (1 + cos theta) / 2
  const sm_Tay<S> sHalf = (sm_Tay<S>(1.0) + b12 / sm_TaySqrt(b11 * b22)) * sm_Tay<S>(0.5);
  if (!sm_Positive(sHalf.c[0]))
    { return FALSE; } // a half turn or more: no single arc
  const sm_Tay<S> w = sm_TaySqrt(sHalf), sInv = sm_Tay<S>(1.0) / sHalf;
  rW[0] = sm_Tay<S>(1.0);
  rW[1] = w;
  rW[2] = sm_Tay<S>(1.0);
  for (int q=0; q<3; q++)
    {
      const sm_Tay<S> sO(crO[q]);
      const sm_Tay<S> sP1 = crM[q] + ((crRho1[q] + crRho2[q]) * sm_Tay<S>(0.5) - crM[q]) * sInv;
      rH[0][q] = crRho1[q] - sO;
      rH[1][q] = w * (sP1 - sO);
      rH[2][q] = crRho2[q] - sO;
    }
  return TRUE;
}

// The arc's row terms over one piece (sm_GnRowTermsOver's, every row varying): the rails' and the
// centre's series over the piece (crF, crCentre to order K) and at its middle (crFm, crCm, to K - 1).
// bTouch: the ball's radius varies, and the circle is the one it touches the envelope along, about
// m + ((rho1 - m) . m' / m' . m') m' (crCentre and crCm one order more, for m').
SmStatus sm_ArcRowTermsOver(const sm_VSer<sm_Ivl> (&crF)[2], const sm_VSer<sm_Ivl> & crCentre, const sm_VSer<sm_Ivl> (&crFm)[2],
                            const sm_VSer<sm_Ivl> & crCm, const sm_Ivl & crD, const double (&crO)[3], sm_GnRowTerms & rT, SmBoolean & rbProved,
                            SmBoolean bTouch = FALSE)
{
  rbProved = FALSE;
  const int K = 4 + SM_GN_TAYLOR;
  sm_Tay<sm_Ivl> asRho[2][3], asM[3];
  for (int q=0; q<3; q++)
    {
      sm_Ivl acMid[K];
      for (int f=0; f<2; f++)
        {
          for (int i=0; i<K; i++) { acMid[i] = crFm[f].a[q].c[i][0]; }
          asRho[f][q] = sm_Centred(acMid, crF[f].a[q].c[K][0], K, crD);
        }
      for (int i=0; i<K; i++) { acMid[i] = crCm.a[q].c[i][0]; }
      asM[q] = sm_Centred(acMid, crCentre.a[q].c[K][0], K, crD);
    }
  if (bTouch)
    {
      // m' from the centre's coefficients one order up: (i + 1) c_(i+1)
      sm_Tay<sm_Ivl> asDm[3];
      for (int q=0; q<3; q++)
        {
          sm_Ivl acMid[K];
          for (int i=0; i<K; i++) { acMid[i] = crCm.a[q].c[i+1][0] * sm_Ivl(static_cast<double>(i + 1)); }
          asDm[q] = sm_Centred(acMid, crCentre.a[q].c[K+1][0] * sm_Ivl(static_cast<double>(K + 1)), K, crD);
        }
      sm_Tay<sm_Ivl> sDD, sBD;
      for (int q=0; q<3; q++)
        {
          sDD = sDD + asDm[q] * asDm[q];
          sBD = sBD + (asRho[0][q] - asM[q]) * asDm[q];
        }
      if (!sm_Positive(sDD.c[0]))
        { return SM_SUCCESS; }
      const sm_Tay<sm_Ivl> sK = sBD / sDD;
      for (int q=0; q<3; q++) { asM[q] = asM[q] + sK * asDm[q]; }
    }
  sm_Tay<sm_Ivl> asH[3][3], asW[3];
  if (!sm_ArcRows(asRho[0], asRho[1], asM, crO, asH, asW))
    { return SM_SUCCESS; }
  // |S - O| <= |m - O| + the ball's radius over the piece
  double dMO = 0.0, dRad = 0.0;
  for (int q=0; q<3; q++)
    {
      const double d = sm_Mag(asM[q].c[0] - sm_Ivl(crO[q]));
      dMO += d * d;
      const double b = sm_Mag(asRho[0][q].c[0] - asM[q].c[0]);
      dRad += b * b;
    }
  const double dR = smos_Sqrt(dMO) + smos_Sqrt(dRad);
  for (size_t j=0; j<3; j++)
    {
      double dSq = 0.0;
      for (int q=0; q<3; q++) { dSq += sm_Mag(asH[j][q].c[SM_GN_TAY]) * sm_Mag(asH[j][q].c[SM_GN_TAY]); }
      const double dX = 24.0 * smos_Sqrt(dSq), dY = 24.0 * sm_Mag(asW[j].c[SM_GN_TAY]), dWlo = asW[j].c[0].lo;
      if (!(std::isfinite(dX) && std::isfinite(dY) && dWlo > 0.0))
        { return SM_SUCCESS; }
      rT.X[j] = std::max(rT.X[j], dX);
      rT.Y[j] = std::max(rT.Y[j], dY);
      rT.Wlo[j] = std::min(rT.Wlo[j], dWlo);
    }
  if (!std::isfinite(dR))
    { return SM_SUCCESS; }
  rT.R = std::max(rT.R, dR);
  rbProved = TRUE;
  return SM_SUCCESS;
}

/*******************************************************************//**
   The bound on the jets' fourth s-derivatives over [crA.dS, crB.dS], per
   face and jet (rM4[f][k], HUGE_VAL where not proved), in their Taylor form
   at the interval's middle m: c''''(s) is within
   sum_{j<M} c^(4+j)(m) (s - m)^j / j! plus the remainder
   c^(4+M)(interval) (h/2)^M / M!, M = SM_GN_TAYLOR, the derivatives at the
   middle enclosed on that one point (the root's box is the root to
   rounding) and the remainder's over the interval.
***********************************************************************/
SmStatus sm_GnFourthOver(const sm_Nurbs * apFace[2], const sm_NurbsCurve & crLeader, const SmSurface * apSurf[2], const SmBSplineCurve & crCenter,
                         const sm_GnStation & crA, const sm_GnStation & crB, const double (&crR)[2], const double (&crSide)[2], int n,
                         double (&rM4)[2][SM_GN_MAXE + 1], double (&rUV4)[2], sm_VSer<sm_Ivl> (*pF)[2] = NULL, sm_VSer<sm_Ivl> * pV = NULL,
                         sm_VSer<sm_Ivl> (*pFm)[2] = NULL, sm_VSer<sm_Ivl> * pVm = NULL, double * pdMid = NULL,
                         const sm_GnLaw * pLaw = NULL, int lExtra = 0)
{
  for (int f=0; f<2; f++)
    {
      for (int k=0; k<=SM_GN_MAXE; k++)
        rM4[f][k] = HUGE_VAL;
      rUV4[f] = HUGE_VAL;
    }
  const double h = crB.dS - crA.dS;
  sm_VSer<sm_Ivl> asF5[2], asF4[2];
  sm_Ser<sm_Ivl> asX5[4], asX4[4];
  double adDS[2] = { 0.0, 0.0 }, adDSm[2] = { 0.0, 0.0 };
  SmBoolean bProved = FALSE;
  SM_GN_TRY(sm_GnEnclose(apFace, crLeader, crA, crB, crR, crSide, n, 4 + SM_GN_TAYLOR + lExtra, asF5, bProved, &asX5, &adDS, pV, pLaw));
  if (!bProved)
    { return SM_SUCCESS; }
  if (pF != NULL)
    { for (int f=0; f<2; f++) { (*pF)[f] = asF5[f]; } }
  // the middle's station, and the series there
  const double dMid = 0.5 * (crA.dS + crB.dS);
  double adSeed[4];
  for (int c=0; c<4; c++) { adSeed[c] = 0.5 * (crA.asX[c].c[0][0] + crB.asX[c].c[0][0]); }
  SmPoint3d sM;
  SM_GN_TRY(crCenter.EvaluatePoint(dMid, sM));
  const sm_Trunc sPointTr = { 0, 0 };
  sm_GnStation sMid;
  sm_Ser<double> sMidR;
  const sm_Ser<double> * pMidR = NULL;
  SM_GN_TRY(sm_LawSeries(pLaw, dMid, 0, sMidR, pMidR));
  SM_GN_TRY(sm_GnSolveStation(apSurf, crCenter, dMid, adSeed, crR, crSide, sPointTr, 1.0e-9 * smos_Max(1.0, SmVector3d(sM).Length()), sMid, pMidR));
  SM_GN_TRY(sm_GnEnclose(apFace, crLeader, sMid, sMid, crR, crSide, n, 4 + SM_GN_TAYLOR - 1 + lExtra, asF4, bProved, &asX4, &adDSm, pVm, pLaw));
  if (!bProved)
    { return SM_SUCCESS; }
  if (pFm != NULL)
    { for (int f=0; f<2; f++) { (*pFm)[f] = asF4[f]; } }
  if (pdMid != NULL)
    { *pdMid = dMid; }
  // the uv curves' fourth derivatives, in the same Taylor form, times the faces' |S_u| + |S_v| over the span's box
  for (int f=0; f<2; f++)
    {
      double dMax = 0.0;
      for (int c=2*f; c<2*f+2; c++)
        {
          double m = 0.0, dFact = 24.0, dPow = 1.0, dJ = 1.0;
          for (int j=0; j<SM_GN_TAYLOR; j++)
            {
              if (j > 0) { dFact *= (4 + j); dPow *= 0.5 * h; dJ *= j; }
              m += dFact * sm_Mag(asX4[c].c[4+j][0]) * dPow / dJ;
            }
          dFact *= (4 + SM_GN_TAYLOR); dPow *= 0.5 * h; dJ *= SM_GN_TAYLOR;
          m += dFact * sm_Mag(asX5[c].c[4+SM_GN_TAYLOR][0]) * dPow / dJ;
          dMax = std::max(dMax, m);
        }
      rUV4[f] = dMax * adDS[f];
    }
  for (int f=0; f<2; f++)
    for (int k=0; k<=n; k++)
      {
        double dSq = 0.0;
        for (int c=0; c<3; c++)
          {
            // c^(4+j)(m) = (4+j)! F[4+j]: the Taylor polynomial's bound over |s - m| <= h/2, and the remainder
            double m = 0.0, dFact = 24.0, dPow = 1.0, dJ = 1.0;
            for (int j=0; j<SM_GN_TAYLOR; j++)
              {
                if (j > 0) { dFact *= (4 + j); dPow *= 0.5 * h; dJ *= j; }
                m += dFact * sm_Mag(asF4[f].a[c].c[4+j][k]) * dPow / dJ;
              }
            dFact *= (4 + SM_GN_TAYLOR); dPow *= 0.5 * h; dJ *= SM_GN_TAYLOR;
            m += dFact * sm_Mag(asF5[f].a[c].c[4+SM_GN_TAYLOR][k]) * dPow / dJ;
            dSq += m * m;
          }
        rM4[f][k] = smos_Sqrt(dSq);
      }
  return SM_SUCCESS;
}

/*******************************************************************//**
   The certified bound on the Hermite span between two stations: rdBound,
   or HUGE_VAL where the span is not proved. The Hermite remainder is
   h^4 / 384 times the jets' largest fourth s-derivative over the span,
   bounded on its pieces (sm_GnFourthOver on each): the span cut in 2^L equal
   pieces, L raised until every piece is proved and, over the target, while
   that halves the bound, so that a long span's bound follows its
   derivatives, not the width of one enclosure over it. The feet's uv curves
   are Hermite interpolants of the same stations, and the face at a uv
   curve's point is off the face's true foot by at most |S_u, S_v| times
   the uv error, h^4 / 384 |x''''|: that bound counts too (an edge's curve
   against its face). pVary: the thumbweight varies along s, and the rows'
   bound (sm_GnRowTermsOver) takes the place of the jets' sum.
***********************************************************************/
SmStatus sm_GnCertifySpan(const sm_Nurbs * apFace[2], const sm_NurbsCurve & crLeader, const SmSurface * apSurf[2], const SmBSplineCurve & crCenter,
                          const sm_GnStation & crA, const sm_GnStation & crB, const double (&crR)[2], const double (&crSide)[2], int n,
                          double dTarget, double & rdBound, const sm_GnVarying * pVary = NULL, const sm_GnLaw * pLaw = NULL)
{
  rdBound = HUGE_VAL;
  const double h = crB.dS - crA.dS;
  double dLast = HUGE_VAL;
  // the ball's circle (n = 1): its arc's rows all vary along s, about the first station's centre
  const SmBoolean bArc = (n == 1) ? TRUE : FALSE;
  const SmBoolean bRows = (pVary != NULL || bArc) ? TRUE : FALSE;
  const int lJets = bArc ? 0 : n; // the jets graded on their own: the arc's rails only
  // the varying thumbweight's fixed one wb and point O: the first station's
  sm_GnSectionFormT<sm_Ivl> sBar;
  double adO[3] = { 0.0, 0.0, 0.0 };
  if (pVary != NULL)
    {
      if (!(std::isfinite(crA.dW) && std::isfinite(crB.dW)))
        { SER(SM_ERR); }
      sBar = sm_GnFormNormalized(pVary->n, sm_Ivl(crA.dW));
      for (int q=0; q<3; q++) { adO[q] = crA.asF[0].a[q].c[0][0]; }
    }
  if (bArc)
    {
      sBar.W.assign(3, sm_Ivl(0.0));
      for (int q=0; q<3; q++) { adO[q] = crA.sV.a[q].c[0][0]; }
    }
  for (int lLevel=0; lLevel<=6; lLevel++)
    {
      const int lPieces = 1 << lLevel;
      // the pieces' end stations
      std::vector<sm_GnStation> sEnds(static_cast<size_t>(lPieces + 1));
      sEnds.front() = crA;
      sEnds.back() = crB;
      SmBoolean bSolved = TRUE;
      for (int i=1; i<lPieces && bSolved; i++)
        {
          const double t = static_cast<double>(i) / lPieces, dS = crA.dS + t * h;
          double adSeed[4];
          for (int c=0; c<4; c++) { adSeed[c] = (1.0 - t) * crA.asX[c].c[0][0] + t * crB.asX[c].c[0][0]; }
          SmPoint3d sM;
          SM_GN_TRY(crCenter.EvaluatePoint(dS, sM));
          const sm_Trunc sPointTr = { 0, 0 };
          sm_Ser<double> sEndR;
          const sm_Ser<double> * pEndR = NULL;
          SM_GN_TRY(sm_LawSeries(pLaw, dS, 0, sEndR, pEndR));
          if (sm_GnSolveStation(apSurf, crCenter, dS, adSeed, crR, crSide, sPointTr, 1.0e-9 * smos_Max(1.0, SmVector3d(sM).Length()),
                                sEnds[static_cast<size_t>(i)], pEndR) != SM_SUCCESS)
            { bSolved = FALSE; }
        }
      if (!bSolved)
        { continue; }
      double adM4[2][SM_GN_MAXE + 1], adUV4[2] = { 0.0, 0.0 };
      for (int f=0; f<2; f++)
        for (int k=0; k<=lJets; k++)
          adM4[f][k] = 0.0;
      SmBoolean bProved = TRUE;
      sm_GnRowTerms sRows;
      sRows.X.assign(sBar.W.size(), 0.0);
      sRows.Y.assign(sBar.W.size(), 0.0);
      sRows.Wlo.assign(sBar.W.size(), HUGE_VAL);
      sRows.R = 0.0;
      for (int i=0; i<lPieces && bProved; i++)
        {
          double adP[2][SM_GN_MAXE + 1], adPUV[2];
          sm_VSer<sm_Ivl> asF[2], sV, asFm[2], sVm;
          double dMid = 0.0;
          const sm_GnStation & crPA = sEnds[static_cast<size_t>(i)], & crPB = sEnds[static_cast<size_t>(i + 1)];
          // the touch circle of a varying ball takes the centre's s-derivative: one order more
          SM_GN_TRY(sm_GnFourthOver(apFace, crLeader, apSurf, crCenter, crPA, crPB, crR, crSide, lJets, adP, adPUV,
                                    bRows ? &asF : NULL, bRows ? &sV : NULL, bRows ? &asFm : NULL, bRows ? &sVm : NULL, &dMid,
                                    pLaw, (bArc && pLaw != NULL) ? 1 : 0));
          for (int f=0; f<2; f++)
            {
              for (int k=0; k<=lJets; k++)
                {
                  if (!(adP[f][k] < HUGE_VAL)) { bProved = FALSE; }
                  adM4[f][k] = std::max(adM4[f][k], adP[f][k]);
                }
              if (!(adPUV[f] < HUGE_VAL)) { bProved = FALSE; }
              adUV4[f] = std::max(adUV4[f], adPUV[f]);
            }
          if (bProved && pVary != NULL)
            {
              // the root's guess at the piece's middle: the stations' thumbweights, linearly
              const double t = (dMid - crA.dS) / h;
              const sm_Ivl sD = sm_Ivl(crPA.dS, crPB.dS) - sm_Ivl(dMid);
              SM_GN_TRY(sm_GnRowTermsOver(*pVary, asF, sV, asFm, sVm, sD, (1.0 - t) * crA.dW + t * crB.dW, sBar, adO, sRows, bProved));
            }
          if (bProved && bArc)
            {
              const sm_Ivl sD = sm_Ivl(crPA.dS, crPB.dS) - sm_Ivl(dMid);
              SM_GN_TRY(sm_ArcRowTermsOver(asF, sV, asFm, sVm, sD, adO, sRows, bProved, pLaw != NULL ? TRUE : FALSE));
            }
        }
      if (!bProved)
        { continue; }
      double dBound = 0.0, dA = 0.0;
      for (int f=0; f<2; f++)
        {
          double dSum = 0.0;
          for (int k=0; k<=lJets; k++) { dSum += h * h * h * h / 384.0 * adM4[f][k]; }
          dA = std::max(dA, dSum);
          dBound = std::max(dBound, std::max(bRows ? 0.0 : dSum, h * h * h * h / 384.0 * adUV4[f]));
        }
      // the arc: no fixed section, every row's remainder (sm_GnRowTermsOver's bound with A = 0)
      if (bArc)
        { dA = 0.0; }
      if (bRows)
        {
          const double K = h * h * h * h / 384.0;
          for (size_t j=0; j<sBar.W.size(); j++)
            {
              const double dDen = sRows.Wlo[j] - K * sRows.Y[j];
              const double dNum = sBar.W[j].hi * dA + K * sRows.X[j] + K * sRows.Y[j] * sRows.R;
              dBound = (dDen > 0.0) ? std::max(dBound, dNum / dDen) : HUGE_VAL;
              if (!(dDen > 0.0))
                { break; }
            }
        }
      if (!std::isfinite(dBound))
        { continue; }
      rdBound = std::min(rdBound, sm_Up(dBound * (1.0 + 1.0e-12)));
      // once proved, one finer level more only where the bound is over the target and the last
      // one halved it
      if (rdBound <= dTarget || dBound > 0.5 * dLast)
        { break; }
      dLast = dBound;
    }
  return SM_SUCCESS;
}


/*******************************************************************//**
   Truncated Taylor series in (sigma, delta t) of any orders, for the
   evaluation of a blend by its definition: c[i (E + 1) + j] the coefficient
   of sigma^i dt^j.
***********************************************************************/
struct sm_T2
{
  int S, E;
  std::vector<double> c;
  sm_T2(int lS = 0, int lE = 0, double d = 0.0) : S(lS), E(lE), c(static_cast<size_t>((lS + 1) * (lE + 1)), 0.0) { c[0] = d; }
  double & at(int i, int j) { return c[static_cast<size_t>(i * (E + 1) + j)]; }
  double   at(int i, int j) const { return c[static_cast<size_t>(i * (E + 1) + j)]; }
};
sm_T2 operator+(const sm_T2 & a, const sm_T2 & b) { sm_T2 r(a.S, a.E); for (size_t k=0; k<r.c.size(); k++) { r.c[k] = a.c[k] + b.c[k]; } return r; }
sm_T2 operator-(const sm_T2 & a, const sm_T2 & b) { sm_T2 r(a.S, a.E); for (size_t k=0; k<r.c.size(); k++) { r.c[k] = a.c[k] - b.c[k]; } return r; }
sm_T2 operator*(const sm_T2 & a, double d) { sm_T2 r(a.S, a.E); for (size_t k=0; k<r.c.size(); k++) { r.c[k] = a.c[k] * d; } return r; }
sm_T2 operator*(const sm_T2 & a, const sm_T2 & b)
{
  sm_T2 r(a.S, a.E);
  for (int i=0; i<=a.S; i++)
    for (int j=0; j<=a.E; j++)
      {
        double d = 0.0;
        for (int i1=0; i1<=i; i1++)
          for (int j1=0; j1<=j; j1++)
            d += a.at(i1, j1) * b.at(i - i1, j - j1);
        r.at(i, j) = d;
      }
  return r;
}
// x^p for a series whose constant term is positive: x0^p (1 + z)^p, z = x / x0 - 1, summed as the
// binomial series in z, which has no constant term (so S + E terms are all there are)
SmStatus sm_T2Pow(const sm_T2 & x, double p, sm_T2 & rY)
{
  const double x0 = x.c[0];
  if (!(x0 > 0.0) || !std::isfinite(x0))
    { SER_MSG(SM_ERR, _T("SmProceduralBlend: its section degenerates at this point (a weight or a length is zero)")); }
  sm_T2 z = x * (1.0 / x0);
  z.c[0] = 0.0;
  sm_T2 sPow(x.S, x.E, 1.0), sSum(x.S, x.E, 1.0);
  double dBinom = 1.0;
  for (int m=1; m<=x.S+x.E; m++)
    {
      dBinom *= (p - (m - 1)) / m;
      sPow = sPow * z;
      sSum = sSum + sPow * dBinom;
    }
  rY = sSum * smos_Pow(x0, p);
  return SM_SUCCESS;
}
struct sm_TV { sm_T2 a[3]; };
sm_TV sm_TVZero(int lS, int lE) { sm_TV r; for (int k=0; k<3; k++) { r.a[k] = sm_T2(lS, lE); } return r; }
sm_TV operator+(const sm_TV & p, const sm_TV & q) { sm_TV r; for (int k=0; k<3; k++) { r.a[k] = p.a[k] + q.a[k]; } return r; }
sm_TV operator-(const sm_TV & p, const sm_TV & q) { sm_TV r; for (int k=0; k<3; k++) { r.a[k] = p.a[k] - q.a[k]; } return r; }
sm_TV operator*(const sm_TV & p, const sm_T2 & s) { sm_TV r; for (int k=0; k<3; k++) { r.a[k] = p.a[k] * s; } return r; }
sm_TV operator*(const sm_TV & p, double d) { sm_TV r; for (int k=0; k<3; k++) { r.a[k] = p.a[k] * d; } return r; }
sm_T2 sm_TVDot(const sm_TV & p, const sm_TV & q) { return p.a[0] * q.a[0] + p.a[1] * q.a[1] + p.a[2] * q.a[2]; }
// a series in sigma of a station (its eps^k coefficients), up to lS, as an sm_T2 of orders (lS, lE)
sm_TV sm_TVOf(const sm_VSer<double> & crV, int k, int lS, int lE)
{
  sm_TV r = sm_TVZero(lS, lE);
  for (int c=0; c<3; c++)
    for (int i=0; i<=lS && i<=SM_GN_MAXS; i++)
      r.a[c].at(i, 0) = crV.a[c].c[i][k];
  return r;
}
// The stations a thread solved last for blends by their definitions: an evaluation at the same s
// (another t: a row of a grid, quadrature, an isocurve) takes its station from here. Keyed by the
// definition's serial (SmProceduralBlend::Renew), never reused, so a station is never another's.
struct sm_StationMemo
{
  unsigned long long lSerial; // 0: empty
  double             dS;
  int                ts, te;
  sm_GnStation       sStation;
  SmBoolean          bW;         // sW holds the station's equal area thumbweight series
  sm_Tay<double>     sW;
};
const int SM_STATION_MEMO = 4;
thread_local sm_StationMemo t_asMemo[SM_STATION_MEMO];
thread_local int t_lMemoNext = 0;

// p^(j)(t0) / j!, j <= lJ, of the polynomial of Bernstein coefficients crB
std::vector<double> sm_BernTaylor(const sm_Bern & crB, double t0, int lJ)
{
  std::vector<double> sOut(static_cast<size_t>(lJ + 1), 0.0);
  sm_Bern sD = crB;
  double dFact = 1.0, dJFact = 1.0;
  for (int j=0; j<=lJ && !sD.empty(); j++)
    {
      if (j > 0) { dJFact *= j; }
      // de Casteljau at t0
      std::vector<double> sW = sD;
      for (size_t r=1; r<sW.size(); r++)
        for (size_t q=0; q+r<sW.size(); q++)
          sW[q] = (1.0 - t0) * sW[q] + t0 * sW[q+1];
      sOut[static_cast<size_t>(j)] = dFact * sW[0] / dJFact;
      // the derivative's coefficients: deg (b_{q+1} - b_q)
      const double dDeg = static_cast<double>(sD.size() - 1);
      sm_Bern sNext;
      for (size_t q=0; q+1<sD.size(); q++) { sNext.push_back(sD[q+1] - sD[q]); }
      dFact *= dDeg;
      sD.swap(sNext);
    }
  return sOut;
}
} // namespace

/*******************************************************************//**
PURPOSE: A blend's section of order lOrder (0 to SM_GN_MAX_ORDER) and
    thumbweight dThumbweight (w < 1, or SM_GN_THUMBWEIGHT_EQUAL_AREA, used from order 2). Order 1's
    circles are built as SmArcSection's defaults say.
***********************************************************************/
SmGnSectionFSG::SmGnSectionFSG(ULONG lOrder, double dThumbweight)
  : m_lOrder(lOrder), m_dThumbweight(dThumbweight), m_dCertifiedBound(-1.0), m_lCertifiedSpans(0)
{
}

namespace
{
// The thumbweight at which the order-n section takes off what the ball's circle does, for half
// the angle u between the ball's radii (EqualAreaThumbweight): the root of I(w) = target(u).
SmStatus sm_GnEqualAreaW(ULONG lOrder, double dU, double & rdW)
{
  rdW = 0.0;
  if (lOrder < 2 || lOrder > SM_GN_MAX_ORDER || !(dU > 0.0 && dU < 0.5 * SM_PI))
    { SER(SM_ERR_INVALID_INPUT); }
  const double dTarget = sm_GnTarget(sm_Tay<double>(dU)).c[0];
  double dLo = -1.0, dHi = 1.0 - 1.0e-9;
  SmBoolean bHeld = (dTarget >= sm_GnAreaSum<double>(lOrder, dHi)) ? TRUE : FALSE;
  while (bHeld && sm_GnAreaSum<double>(lOrder, dLo) < dTarget)
    {
      dHi = dLo;
      dLo = 2.0 * dLo - 1.0;
      bHeld = (dLo >= -1.0e3) ? TRUE : FALSE;
    }
  if (!bHeld)
    { SER_MSG(SM_ERR_INVALID_INPUT, _T("SmGnSectionFSG: no thumbweight gives this wedge's section the ball's circle's area")); }
  for (int it=0; it<200 && dHi - dLo > 1.0e-15 * (1.0 - dLo); it++)
    {
      const double dMid = 0.5 * (dLo + dHi);
      if (sm_GnAreaSum<double>(lOrder, dMid) > dTarget) { dLo = dMid; } else { dHi = dMid; }
    }
  rdW = 0.5 * (dLo + dHi);
  return SM_SUCCESS;
}

// A station's equal area thumbweight and its s-derivative: x from the station's feet and V to
// first order in s, w the root at its u, dw/ds = target' / I'(w).
SmStatus sm_GnStationThumbweight(ULONG n, double dRR, sm_GnStation & rSt)
{
  const sm_Tay<double> u = sm_GnHalfAngle(sm_GnWedgeCos<double>(rSt.asF, rSt.sV, dRR, 1));
  double dW = 0.0;
  SER(sm_GnEqualAreaW(n, u.c[0], dW));
  const sm_Tay<double> sD = sm_GnAreaSum(n, sm_TayOf(dW, 1.0));
  if (!(sD.c[1] < 0.0))
    { SER_MSG(SM_ERR, _T("SmGnSectionFSG: the section's area does not fall with its thumbweight")); }
  const sm_Tay<double> sW = sm_GnInvert(sm_GnTarget(u), sD, dW);
  rSt.dW = dW;
  rSt.dWds = sW.c[1];
  return SM_SUCCESS;
}
}

/*******************************************************************//**
PURPOSE: The equal area thumbweight: the w at which the order-n section
    between two planes at the angle alpha takes off what the ball's circle
    does on the same rails.
NOTES: With e1, e2 the rails' points from the edge (|e_i| = d = r cot(alpha/2)),
    the section takes off (1/2) d^2 sin(alpha) I(w), I the integral of
    a b' - a' b (sm_GnAreaSum), and the circle r^2 (cot(alpha/2)
    - (pi - alpha)/2). So I(w) = (cot(alpha/2) - (pi - alpha)/2) /
    (cot(alpha/2) cos^2(alpha/2)), a target that depends on alpha only:
    with u = (pi - alpha)/2, half the angle between the ball's radii to its
    contacts, (sin u - u cos u) / sin^3 u. I falls as w rises (to 0 as w
    nears 1) and nears the chord's 1 as w falls, so the root is bisected on
    a bracket whose lower end steps down from -1 until it holds it: a sharp
    wedge's circle takes off nearly the chord's area, and its w is well
    below -1. A target the bracket does not reach by w = -1000 (a wedge so
    open or so shut that the circle takes off next to nothing, or all but
    the chord's area) is an error naming it.
***********************************************************************/
SmStatus SmGnSectionFSG::EqualAreaThumbweight(ULONG lOrder, double dAlpha, double & rdW)
{
  rdW = 0.0;
  if (!(dAlpha > 0.0 && dAlpha < SM_PI))
    { SER(SM_ERR_INVALID_INPUT); }
  SER(sm_GnEqualAreaW(lOrder, 0.5 * (SM_PI - dAlpha), rdW));
  return SM_SUCCESS;
}

/*******************************************************************//**
PURPOSE: The thumbweight at a station whose ball's radii to its contacts
    meet at the angle whose cosine is dCosRadii: the one given, or the
    equal area one of the wedge alpha = pi - (that angle). Exact for two
    planes; elsewhere the section of the station's own wedge.
***********************************************************************/
SmStatus SmGnSectionFSG::ThumbweightAt(double dCosRadii, double & rdW) const
{
  rdW = m_dThumbweight;
  if (m_dThumbweight != SM_GN_THUMBWEIGHT_EQUAL_AREA)
    { return SM_SUCCESS; }
  SER(sm_GnEqualAreaW(m_lOrder, 0.5 * sm_Acos(dCosRadii), rdW));
  return SM_SUCCESS;
}

/*******************************************************************//**
PURPOSE: The order-1 section, the ball's circle, built as crArcs says.
***********************************************************************/
SmGnSectionFSG::SmGnSectionFSG(const SmArcSection & crArcs)
  : m_sArcs(crArcs), m_lOrder(1), m_dThumbweight(0.0), m_dCertifiedBound(-1.0), m_lCertifiedSpans(0)
{
}


/*******************************************************************//**
PURPOSE: The interior knots where a fillet's center curve or one of its rails
    turns: its tangents before and after the knot differ by more than
    dAngTolRad.

NOTES: SolidCraft change. The three curves share their knots (FlushCurve
    builds them span by span); a turn comes from a corner in the radius.
***********************************************************************/
static SmStatus sm_FilletCurveCorners
  (const SmBSplineCurve & crCenterLine,  // in : fillet center curve
   const SmBSplineCurve & crRail1Curve,  // in : 1st rail
   const SmBSplineCurve & crRail2Curve,  // in : 2nd rail
   double                 dAngTolRad,    // in : angle above which a knot is a corner
   SmTArray<double>     & rCorners)      // out: corner knots, ascending
{
  rCorners.ReSet();
  SmTArray<double> sKnots;
  SER(crRail1Curve.GetKnots(sKnots));
  SmExtent1d sIvl = crRail1Curve.GetNaturalInterval();
  double dTol = SM_EFF_ZERO_SQRT * (1.0 + smos_Fabs(sIvl.GetMin()) + smos_Fabs(sIvl.GetMax()));
  const SmBSplineCurve * apCurves[3] = { &crCenterLine, &crRail1Curve, &crRail2Curve };
  for (ULONG i=0; i<sKnots.GetSize(); i++)
    {
      double dKnot = sKnots[i];
      if (   (i > 0 && sKnots[i-1] == dKnot)
          || dKnot < sIvl.GetMin() + dTol || dKnot > sIvl.GetMax() - dTol)
        { continue; }
      SmBoolean bCorner = FALSE;
      for (ULONG j=0; j<3 && !bCorner; j++)
        {
          if (!apCurves[j]->GetNaturalInterval().ContainsValue(dKnot))
            { continue; }
          SmVector3d sBefore[2], sAfter[2];
          SER(apCurves[j]->Evaluate(dKnot, 1, FALSE, sBefore, FALSE));
          SER(apCurves[j]->Evaluate(dKnot, 1, TRUE,  sAfter,  FALSE));
          if (   sBefore[1].LengthSquared() < SM_EFF_ZERO_SQ
              || sAfter[1].LengthSquared()  < SM_EFF_ZERO_SQ)
            { continue; }
          double dAngle;
          SER(sBefore[1].AngleBetween(sAfter[1], dAngle));
          if (dAngle > dAngTolRad)
            { bCorner = TRUE; }
        }
      if (bCorner)
        { rCorners.Add(dKnot); }
    }
  return SM_SUCCESS;

} // end sm_FilletCurveCorners


/*******************************************************************//**
PURPOSE: Whether the ball's circles along a centre line and its rails sweep a
    cylinder (all three parallel lines) or revolve a torus (coaxial arcs):
    the exact surfaces of order 1 the circle's surface is built as.
***********************************************************************/
static SmStatus sm_CircleSweptOrRevolved(const SmBSplineCurve & crCenterLine, const SmBSplineCurve & crRail1Curve,
                                         const SmBSplineCurve & crRail2Curve, double dThisApproxTol3d,
                                         SmBoolean & rbSweep, SmBoolean & rbRevolve, SmAxis2Placement & sRefFrame, double & dEndAng)
{
  rbSweep = FALSE;
  rbRevolve = FALSE;
  // see if newSurf == RuledSurface
  double    dAnalyticTol = 1.0e-8;
  SmPoint3d  sLinePoint;
  SmVector3d sLineVector;
  SmBoolean bIsDeg1Line =    crCenterLine.GetDegree() == 1
                          && crCenterLine.IsLine(11,dAnalyticTol,sLinePoint,sLineVector);
  if (bIsDeg1Line)
    {
      rbSweep = TRUE;
      SmPoint3d  sLinePoint2;
      SmVector3d sLineVector2;
      SmBoolean bIsLineR1 = crRail1Curve.IsLine(11,dAnalyticTol,sLinePoint2,sLineVector2);
      double dAngle;
      if ( bIsLineR1 )
        {
          SER(sLineVector2.AngleBetween(sLineVector,dAngle));
          if (dAngle > SM_EFF_ZERO_SQRT)
            {
              rbSweep = FALSE;
            }
        }
      SmBoolean bIsLineR2 = crRail2Curve.IsLine(11,dAnalyticTol,sLinePoint2,sLineVector2);
      if ( bIsLineR2 )
        {
          SER(sLineVector2.AngleBetween(sLineVector,dAngle));
          if (dAngle > SM_EFF_ZERO_SQRT)
            {
              rbSweep = FALSE;
            }
        }
    } // end NewSurf == SurfaceOfRevolution check

  // see if NewSurf == SurfaceOfRevolution
  double dArcRadius, dStartAng;
  SmBoolean bIsArc = crCenterLine.IsArc(11, dAnalyticTol,
                                       sRefFrame, dArcRadius,
                                       dStartAng, dEndAng);
  if (bIsArc)
    {
      //cbi this is never hit in prog_test.
      rbRevolve = TRUE;
      double dRadius2, dStartAng2, dEndAng2;
      double dRadius3, dStartAng3, dEndAng3;
      SmAxis2Placement sRefFrame2, sRefFrame3;
      SmBoolean bIsArc2 = crRail1Curve.IsArc(11, dAnalyticTol, sRefFrame2,
                                             dRadius2, dStartAng2, dEndAng2);
      if ( !bIsArc2 )
        { rbRevolve = FALSE; }
      else if ( smos_Fabs( dEndAng2 - dEndAng ) * dRadius2 > dThisApproxTol3d )
        { rbRevolve = FALSE; }
      else if ( ! sRefFrame2.AreCoaxial( sRefFrame, dThisApproxTol3d, dRadius2 ) )
        { rbRevolve = FALSE; }

      if ( rbRevolve )
      {
          SmBoolean bIsArc3 = crRail2Curve.IsArc(11, dAnalyticTol, sRefFrame3,
                                                 dRadius3, dStartAng3, dEndAng3);
          if ( !bIsArc3 )
            { rbRevolve = FALSE; }
          else if ( smos_Fabs( dEndAng3 - dEndAng ) * dRadius3 > dThisApproxTol3d )
            { rbRevolve = FALSE; }
          else if ( ! sRefFrame3.AreCoaxial( sRefFrame, dThisApproxTol3d, dRadius3 ) )
            { rbRevolve = FALSE; }
          // Also check against 1s Rail curve.
          else if ( smos_Fabs( dEndAng3 - dEndAng2 ) * dRadius3 > dThisApproxTol3d )
            { rbRevolve = FALSE; }
          else if ( ! sRefFrame3.AreCoaxial( sRefFrame2, dThisApproxTol3d, dRadius3 ) )
            { rbRevolve = FALSE; }
      }
    }

  return SM_SUCCESS;
} // end sm_CircleSweptOrRevolved

/*******************************************************************//**
PURPOSE: The variable radius solver's law as a curve over the edge's
    parameter t: the radius its x, exactly the solver's (its law curve
    mapped by the law's interval, edge map and orientation, all affine in t;
    a linear law as a line). NULL where the law is of another kind.
***********************************************************************/
static SmStatus sm_LawInParameter(const SmContext & crContext, SmVariableRadiusFS & rSolver, SmBSplineCurve *& rpLaw)
{
  rpLaw = NULL;
  SmFilletLaw * pLaw = rSolver.GetFilletLaw();
  NER(pLaw);
  const SmExtent1d sIvl = rSolver.GetFilletLawInterval();
  const SmBoolean bReverse = rSolver.GetFilletLawOrientation();
  if (!(sIvl.GetLength() > 0.0) || !(pLaw->m_vEdgeMap.GetLength() > 0.0))
    { return SM_SUCCESS; }
  // the law's parameter at t: m0 + (m1 - m0) T, T = (t - t0) / (t1 - t0), or (t1 - t) / (t1 - t0) reversed
  const double t0 = sIvl.GetMin(), t1 = sIvl.GetMax(), m0 = pLaw->m_vEdgeMap.GetMin(), m1 = pLaw->m_vEdgeMap.GetMax();
  const double dB = (m1 - m0) / (t1 - t0) * (bReverse ? -1.0 : 1.0);
  const double dA = bReverse ? m0 + (m1 - m0) * t1 / (t1 - t0) : m0 - (m1 - m0) * t0 / (t1 - t0);
  SmBSplineCurve * pCurve = NULL;
  if (pLaw->IsKindOf(SmBSplineFilletLaw_TYPE))
    {
      const SmBSplineCurve * pLawCurve = static_cast<SmBSplineFilletLaw *>(pLaw)->GetLawCurve();
      NER(pLawCurve);
      ULONG lDim = 0, lDeg = 0;
      SmTArray<SmPoint3d> sCPs;
      SmBSplineCurveForm eForm;
      SmTArray<ULONG> sMult;
      SmTArray<double> sKnots, sWts;
      SmKnotType eKnots;
      SER(pLawCurve->GetCanonical(lDim, lDeg, sCPs, eForm, sMult, sKnots, eKnots, sWts));
      // the knots taken to t (t = (k - a) / b), in their order along t
      SmTArray<SmPoint3d> sTCPs;
      SmTArray<ULONG> sTMult;
      SmTArray<double> sTKnots, sTWts;
      const ULONG lK = sKnots.GetSize(), lC = sCPs.GetSize();
      for (ULONG i=0; i<lK; i++)
        {
          const ULONG k = (dB > 0.0) ? i : lK - 1 - i;
          sTKnots.Add((sKnots[k] - dA) / dB);
          sTMult.Add(sMult[k]);
        }
      for (ULONG i=0; i<lC; i++)
        {
          const ULONG k = (dB > 0.0) ? i : lC - 1 - i;
          sTCPs.Add(sCPs[k]);
          if (sWts.GetSize() == lC) { sTWts.Add(sWts[k]); }
        }
      SER(SmBSplineCurve::CreateCanonical(crContext, lDim, lDeg, sTCPs, SM_CF_UNSPECIFIED, sTMult, sTKnots, SM_KT_UNSPECIFIED,
                                          sTWts.GetSize() == lC ? &sTWts : NULL, NULL, pCurve));
    }
  else if (pLaw->IsKindOf(SmLinearFilletLaw_TYPE))
    {
      double adR0[3], adR1[3];
      SER(pLaw->Evaluate(t0, sIvl, bReverse, adR0));
      SER(pLaw->Evaluate(t1, sIvl, bReverse, adR1));
      SmTArray<SmPoint3d> sCPs;
      sCPs.Add(SmPoint3d(adR0[0], 0.0, 0.0));
      sCPs.Add(SmPoint3d(adR1[0], 0.0, 0.0));
      SmTArray<ULONG> sMult;
      sMult.Add(2);
      sMult.Add(2);
      SmTArray<double> sKnots;
      sKnots.Add(t0);
      sKnots.Add(t1);
      SER(SmBSplineCurve::CreateCanonical(crContext, 3, 1, sCPs, SM_CF_UNSPECIFIED, sMult, sKnots, SM_KT_UNSPECIFIED, NULL, NULL, pCurve));
    }
  if (pCurve != NULL)
    { pCurve->SetOutOfBoundsEnabled(TRUE); }
  rpLaw = pCurve;
  return SM_SUCCESS;
}

/*******************************************************************//**
PURPOSE: Create a fillet surface with a circular cross section.  This
    is a lofted surface or in some special cases a cylinder or a torus.

NOTES:
***********************************************************************/
SmStatus SmGnSectionFSG::CircleSurface
  (SmFilletSolver & rFilletSolver,
   const SmTArray<SmTsectPnt*> & crFilletPoints,
   SmBSplineCurve * pCenterLine,
   SmBSplineCurve * pRail1Curve,
   SmBSplineCurve * pRail2Curve,
   SmSurface *pSurf1,
   SmBSplineCurve *pUV1,
   SmSurface *pSurf2,
   SmBSplineCurve *pUV2,
   SmBSplineSurface *& rpFilletSurface)
{
    double dThisApproxTol3d = rFilletSolver.GetThisApproxTol3d();

    SmBoolean bNeedG1Continuity = FALSE;
    if (rFilletSolver.GetSolverType() == SM_FS_CONST_RADIUS) {
        //bNeedG1Continuity = TRUE; // fix it later
    }
    // SolidCraft change: the solvers whose ball changes its radius along the round (a radius law,
    // or a constant chord between faces whose angle changes); the others' balls touch the round
    // along their great circles.
    const SmFilletSolverType eSolver = rFilletSolver.GetSolverType();
    const SmBoolean bBallRadiusVaries = (eSolver == SM_FS_VARIABLE_RADIUS || eSolver == SM_FS_CONST_DIST) ? TRUE : FALSE;

    // SolidCraft change: a skinned surface interpolates its cross sections smoothly, so where
    // the center curve or a rail turns (a corner in a variable radius) it would overshoot.
    // The surface is skinned piece by piece between such corners and the pieces joined there.
    // A turn is a corner past SM_CONTINUITY_ANGLE, the most an edge's curve may turn at a knot
    // inside it (SmBSplineCurve's check): a gentle corner (2 -> 3 -> 2 over 20 turns the rails
    // 11 degrees, under the solver's 20 degree tolerance) would otherwise be skinned over, the
    // round off its radius and its edges kinked inside.
    SmTArray<double> sCorners;
    SER(sm_FilletCurveCorners(*pCenterLine, *pRail1Curve, *pRail2Curve,
                              smos_Min(rFilletSolver.GetThisAngTolRad(), SM_DEG2RAD(SM_CONTINUITY_ANGLE)), sCorners));

    if (sCorners.GetSize() == 0)
    {
        // SolidCraft change: a ball of one radius between faces that are NURBS surfaces is the
        // round by its definition (Build: stations solved on the faces, the arc's rows certified
        // span by span), unless the circles sweep a cylinder or revolve a torus, which stay those
        // exact surfaces, or the rails meet at an end (the section fades to a point there, which
        // the definition does not take yet), or the solver's rails are not the ball's feet on the two
        // faces everywhere (a rail held on a face's boundary while the ball rolls on its edge).
        SmBoolean bDefinition = FALSE;
        double dRadius = 0.0;
        SmSurface * apBase[2] = { NULL, NULL };
        SmBSplineCurve * pLeader = NULL, * pLaw = NULL;
        SmObjDelete sDelLeader, sDelLaw;
        {
            const SmOffsetSurface * apOff[2] = { SM_CAST_PTR(SmOffsetSurface, pSurf1), SM_CAST_PTR(SmOffsetSurface, pSurf2) };
            bDefinition = (apOff[0] != NULL && apOff[1] != NULL && pUV1 != NULL && pUV2 != NULL
                           && !m_sArcs.bApproximateArcs && !m_sArcs.bMirror && !m_sArcs.bComplement) ? TRUE : FALSE;
            for (int f=0; f<2 && bDefinition; f++)
            {
                apBase[f] = apOff[f]->GetBaseSurface();
                if (   apBase[f] == NULL || !apBase[f]->IsKindOf(SmBSplineSurface_TYPE) || apBase[f]->IsKindOf(SmProceduralBlend_TYPE))
                    { bDefinition = FALSE; }
            }
            // the ball's radius: the constant radius solver's, the same for both faces; or the
            // variable radius solver's law over the edge's curve (sm_LawInParameter)
            SmConstantRadiusFS * pConstant = (rFilletSolver.GetSolverType() == SM_FS_CONST_RADIUS && !bBallRadiusVaries)
                                             ? static_cast<SmConstantRadiusFS *>(&rFilletSolver) : NULL;
            SmVariableRadiusFS * pVariable = (rFilletSolver.GetSolverType() == SM_FS_VARIABLE_RADIUS)
                                             ? static_cast<SmVariableRadiusFS *>(&rFilletSolver) : NULL;
            if (bDefinition && pConstant != NULL)
            {
                dRadius = pConstant->GetFilletRadius(0);
                if (!(dRadius > dThisApproxTol3d) || pConstant->GetFilletRadius(1) != dRadius)
                    { bDefinition = FALSE; }
            }
            else if (bDefinition && pVariable != NULL && pVariable->GetEdgeuse(0) != NULL)
            {
                // (an edge's round: a variable ball a corner rolls about its own curve stays the
                // corner's, step 4's)
                const SmBSplineCurve * pEdge = SM_CAST_PTR(SmBSplineCurve, pVariable->GetBaseCurve());
                if (pEdge == NULL || crFilletPoints.GetSize() < 2)
                    { bDefinition = FALSE; }
                else
                {
                    SmCurve * pCopy = NULL;
                    SER(pEdge->Copy(rFilletSolver.GetCreationContext(), pCopy));
                    sDelLeader.SetObj(pCopy);
                    pLeader = SM_CAST_PTR(SmBSplineCurve, pCopy);
                    NER(pLeader);
                    pLeader->SetOutOfBoundsEnabled(TRUE);
                    SER(sm_LawInParameter(rFilletSolver.GetCreationContext(), *pVariable, pLaw));
                    sDelLaw.SetObj(pLaw);
                    if (pLaw == NULL)
                        { bDefinition = FALSE; }
                    else
                    {
                        // the solver's stations run past the edge's ends, the law continued there: a ball
                        // of the law's radius exists over them only where it stays positive (an
                        // enclosure of the law over their range)
                        double t0 = SM_BIG_DOUBLE, t1 = -SM_BIG_DOUBLE;
                        for (ULONG i=0; i<crFilletPoints.GetSize(); i++)
                        {
                            NER(crFilletPoints[i]);
                            t0 = smos_Min(t0, crFilletPoints[i]->m_adUserDoubles[0]);
                            t1 = smos_Max(t1, crFilletPoints[i]->m_adUserDoubles[0]);
                        }
                        sm_NurbsCurve sLawN;
                        SER(sm_NurbsCurveOf(*pLaw, sLawN));
                        std::vector< sm_V3<sm_Ivl> > sR;
                        if (!(t1 > t0) || sm_CurveDerivsOver(sLawN, sm_Ivl(t0, t1), 0, sR) != SM_SUCCESS || !(sR[0].x.lo > dThisApproxTol3d))
                            { bDefinition = FALSE; }
                    }
                }
            }
            else
                { bDefinition = FALSE; }
            if (bDefinition)
            {
                SmBoolean bSweep = FALSE, bRevolve = FALSE;
                SmAxis2Placement sFrame;
                double dEnd = 0.0;
                SER(sm_CircleSweptOrRevolved(*pCenterLine, *pRail1Curve, *pRail2Curve, dThisApproxTol3d, bSweep, bRevolve, sFrame, dEnd));
                const SmExtent1d sIvl = pRail1Curve->GetNaturalInterval();
                for (int lEnd=0; lEnd<2 && bDefinition; lEnd++)
                {
                    SmPoint3d sEnd1, sEnd2;
                    const double dT = (lEnd == 0) ? sIvl.GetMin() : sIvl.GetMax();
                    SER(pRail1Curve->EvaluatePoint(dT, sEnd1));
                    SER(pRail2Curve->EvaluatePoint(dT, sEnd2));
                    if (sEnd1.DistanceBetween(sEnd2) < dThisApproxTol3d)
                        { bDefinition = FALSE; }
                }
                if (bSweep || bRevolve)
                    { bDefinition = FALSE; }
            }
        }
        if (bDefinition)
        {
            SmBoolean bNotTheFeet = FALSE;
            SmGnLawRound sLawRound = { pLeader, pLaw, &crFilletPoints };
            SER(Build(rFilletSolver.GetCreationContext(), *pCenterLine, *pRail1Curve, *pRail2Curve, dThisApproxTol3d,
                      apBase[0], pUV1, apBase[1], pUV2, TRUE, rpFilletSurface, dRadius, &bNotTheFeet, (pLaw != NULL) ? &sLawRound : NULL));
            if (!bNotTheFeet)
                { return SM_SUCCESS; }
        }
        SER(CircleSurfaceFromCurves(rFilletSolver.GetCreationContext(),
            *pCenterLine,*pRail1Curve,*pRail2Curve,dThisApproxTol3d,
            bNeedG1Continuity,bBallRadiusVaries,pSurf1,pUV1,pSurf2,pUV2,rpFilletSurface));
        return SM_SUCCESS;
    }

    const SmContext & crContext = rFilletSolver.GetCreationContext();
    SmExtent1d sIvl = pRail1Curve->GetNaturalInterval();
    SmTArray<double> sEnds;
    sEnds.Add(sIvl.GetMin());
    sEnds.Append(sCorners);
    sEnds.Add(sIvl.GetMax());

    SmBSplineSurface * pJoined = NULL;
    SmObjDelete sDeleteJoined;
    // a ball whose radius varies: the pieces are kept exact for the corners (SmLawCornerPieces)
    std::unique_ptr<SmLawCornerPieces> sKept(bBallRadiusVaries ? new SmLawCornerPieces : NULL);
    SmLawCornerPieces * pKept = sKept.get();
    for (ULONG i=0; i+1<sEnds.GetSize(); i++)
    {
        SmExtent1d sPiece(sEnds[i], sEnds[i+1]);
        SmBSplineCurve * apPieces[3] = { NULL, NULL, NULL };
        SmBSplineCurve * apWholes[3] = { pCenterLine, pRail1Curve, pRail2Curve };
        SmObjDelete asDelete[3];
        for (ULONG j=0; j<3; j++)
        {
            SmCurve * pCopy = NULL;
            SER(apWholes[j]->Copy(crContext, pCopy));
            asDelete[j].SetObj(pCopy);
            apPieces[j] = SM_CAST_PTR(SmBSplineCurve, pCopy); NER(apPieces[j]);
            SmExtent1d sTrim = sPiece;
            SER(apPieces[j]->Trim(sTrim));
        }

        SmBSplineSurface * pPiece = NULL;
        SER(CircleSurfaceFromCurves(crContext, *apPieces[0], *apPieces[1], *apPieces[2],
            dThisApproxTol3d, bNeedG1Continuity, bBallRadiusVaries, pSurf1, pUV1, pSurf2, pUV2, pPiece));
        NER(pPiece);
        SmObjDelete sDeletePiece(pPiece);
        if (pKept != NULL)
        {
            SmSurface * pCopy = NULL;
            SER(pPiece->Copy(crContext, pCopy));
            SmBSplineSurface * pExact = SM_CAST_PTR(SmBSplineSurface, pCopy);
            if (pExact == NULL) { delete pCopy; return SM_ERR; }
            pKept->m_vPieces.Add(pExact);
        }

        if (pJoined == NULL)
        {
            sDeletePiece.Clear();
            pJoined = pPiece;
            sDeleteJoined.SetObj(pJoined);
        }
        else
        {
            // join the piece on at its u-min side; no knot removal there, the surface turns
            SmBSplineSurface * pNext = NULL;
            SER(pJoined->JoinSurface(crContext, pPiece, SM_SP_U, 0.0, pNext));
            NER(pNext);
            delete pJoined;
            pJoined = pNext;
            sDeleteJoined.SetObj(pJoined);
        }
        if (pKept != NULL)
            { pKept->m_vJoinedEnds.Add(pJoined->GetNaturalUVDomain().GetMax().x); }
    }

    // Each corner: the corner ball (its centre on the centre line, through the rails) and how its
    // two sides' touch circles lie. The side before ends on its circle at the corner, the side
    // after starts on its own; both run through the rail points there. Along the rail, the circle
    // after lies ahead of the one before where the law dips (the ball's sphere between them), and
    // behind it where the law peaks (the sides cross).
    if (pKept != NULL)
    {
        for (ULONG i=0; i+1<pKept->m_vPieces.GetSize(); i++)
        {
            const double dT = sEnds[i+1];
            SmVector3d sCentPV[2], sRailPV[2];
            SER(pCenterLine->Evaluate(dT, 1, TRUE, sCentPV));
            SER(pRail1Curve->Evaluate(dT, 1, TRUE, sRailPV, FALSE));
            SmPoint3d sCentre(sCentPV[0]);
            SmVector3d sAlong = sRailPV[1];
            if (sAlong.Length() <= SM_EFF_ZERO) { return SM_ERR; }
            sAlong = sAlong / sAlong.Length();
            SmPoint3d asMid[2];
            for (ULONG j=0; j<2; j++)
            {
                const SmBSplineSurface * pSide = pKept->m_vPieces[i+j];
                SmExtent2d sDomain = pSide->GetNaturalUVDomain();
                SmPoint2d sUV(j == 0 ? sDomain.GetMax().x : sDomain.GetMin().x, 0.5 * (sDomain.GetMin().y + sDomain.GetMax().y));
                SER(pSide->EvaluatePoint(sUV, asMid[j]));
            }
            const double dAhead = (asMid[1] - asMid[0]).Dot(sAlong);
            pKept->m_vKinds.Add(dAhead > dThisApproxTol3d ? SmLawCornerPieces::SM_LC_DIP
                              : dAhead < -dThisApproxTol3d ? SmLawCornerPieces::SM_LC_PEAK : SmLawCornerPieces::SM_LC_MATCHING);
            pKept->m_vCentres.Add(sCentre);
            pKept->m_vRadii.Add((SmPoint3d(sRailPV[0]) - sCentre).Length());
        }
        // the pieces' sections are the balls' circles within the skin's tolerance, and within the
        // arcs' own where they are approximated (a share of the radius)
        pKept->m_dTol3d = dThisApproxTol3d;
        for (ULONG i=0; i<pKept->m_vRadii.GetSize() && m_sArcs.bApproximateArcs; i++)
            { pKept->m_dTol3d = smos_Max(pKept->m_dTol3d, 0.01 * m_sArcs.dTolPercentRadius * pKept->m_vRadii[i]); }
        rFilletSolver.AddLawCornerPieces(sKept.release());
    }

    sDeleteJoined.Clear();
    rpFilletSurface = pJoined;
    return SM_SUCCESS;

} // end SmGnSectionFSG::CircleSurface


/*******************************************************************//**
PURPOSE: Create a surface from a centerline curve and two rail curves.
    Users can specify whether G1 condition are required on both ends.
    Typically, in the cases of constant-radius filleting, it is desirable
    to have the cross-boundary derivatives(or tangent-fields) on both ends
    of the fillets be normal to the cross-section planes.

NOTES: All curves should be parameterized in the same direction
     with the same knot vector.

METHOD ---
   1. Generate CrossSections for every center/rail curve knot value.
      Every crossSection curve must have the same knot vector.
      a. When railPoints are on a circle centered on centerPoint,
         crossSection is a cricular arc.
      b. When railPoints are not on circle,
         crossSection is an interpolating curve made either by
         modifying the circle endPoints or by building an interpolating curve.
      c. When railPoints are coincident,
         cross section is a degenerate curve
   2. Make the new Surface
      a. When center, rail1, and rail2 curves are parallel lines,
         NewSurface = Ruled Surface from Start CrossSection to End CrossSection
      b. When center, rail1, and rail2 crves are arcs,
         NewSurface = Surface of Rotation
      c. Else
         NewSurface = skin surface that interpolates all cross sections
   3. When bNeedG1Continuity == TRUE and NewSurface is a skinSurface,
         The skin surface end CrossTangents are computed as an interpolation of the
         rail1 and rail2 endTangent values.

***********************************************************************/
SmStatus SmGnSectionFSG::CircleSurfaceFromCurves
  (const SmContext   & crContext,          // in : context for new object construction
   SmBSplineCurve    & rCenterLine,        // in : Loci of centerPoints
   SmBSplineCurve    & rRail1Curve,        // in : 1st NewSurface Boundary
   SmBSplineCurve    & rRail2Curve,        // in : 2nd NewSurface Boundary
   double              dThisApproxTol3d,   // in : max allowed dist between NewSurface and idealized surface
   SmBoolean           bNeedG1Continuity,  // in : TRUE =
                                           //      FALSE=
   SmBoolean           bBallRadiusVaries,  // in : TRUE = the ball's radius changes along the centre line (SolidCraft change)
   SmSurface         *,                    // in :
   SmBSplineCurve    *,                    // in :
   SmSurface         *,                    // in :
   SmBSplineCurve    *,                    // in :
   SmBSplineSurface *& rpFilletSurface)    // out: New Surface
{

  // Init outputs
  rpFilletSurface = NULL;

#ifdef SM_DEBUG_CODE
SmBoolean bDebugMe=FALSE;
  if (bDebugMe)
    {
      SmEdge *pCenterEdge = (SmEdge *)rCenterLine.GetEdge() ;
      SmEdge *pRailEdge1  = (SmEdge *)rRail1Curve.GetEdge() ;
      SmEdge *pRailEdge2  = (SmEdge *)rRail2Curve.GetEdge() ;

      SmBrep *pBrep1 = pCenterEdge ? pCenterEdge->GetBrep() : NULL ;
      SmBrep *pBrep2 = pRailEdge1  ? pRailEdge1->GetBrep()  : NULL ;
      SmBrep *pBrep3 = pRailEdge2  ? pRailEdge2->GetBrep()  : NULL ;

      smgfx_Erase() ;
      smgfx_SetLook(1,2, 0,0,1) ; if(pBrep1) pBrep1->Draw(TRUE) ; sm_GraphicsLoop() ;
      smgfx_SetLook(1,2, 0,1,0) ; if(pBrep2 && pBrep2 != pBrep1) pBrep2->Draw(TRUE) ; sm_GraphicsLoop() ;
      smgfx_SetLook(1,2, 1,1,0) ; if(pBrep3 && pBrep3 != pBrep2 && pBrep3 != pBrep1) pBrep3->Draw(TRUE) ; sm_GraphicsLoop() ;
      smgfx_SetLook(3,4, 1,0,0) ; rCenterLine.Draw() ; sm_GraphicsLoop() ;
      smgfx_SetLook(3,5, 0,1,0) ; rRail1Curve.Draw() ; sm_GraphicsLoop() ;
      smgfx_SetLook(3,6, 0,0,1) ; rRail2Curve.Draw() ; sm_GraphicsLoop() ;
      sm_GraphicsLoop() ;
    }
#endif

  // locals
  SmTArray<double> sBreaks;
  SmTArray<double> sKnots, sKnots2;
  SER(rRail1Curve.GetKnots(sKnots));
  SER(rRail2Curve.GetKnots(sKnots2));

  // make sure all 3 curves have compatible knot vectors
  if (sKnots2.GetSize() != sKnots.GetSize())
    {
      SmTArray<SmBSplineCurve*> sRails;
      sRails.Add(&rRail1Curve);
      sRails.Add(&rRail2Curve);
      sRails.Add(&rCenterLine);
      SER(SmBSplineCurve::SyncronizeKnotsOfCurves(sRails,SM_EFF_ZERO_SQRT));
    }

  // local array for cross section curves and parameters
  SmTArray<SmBSplineCurve*> sCrossSections(sKnots.GetSize()*2);
  SmObjsDelete<SmBSplineCurve*> sCleanCrossSections(&sCrossSections);
  SmTArray<double> sParametersOfCrossSections;

  // Surface creation requires 4 or more cross sections so add
  // some extra knots if needed.

  // insert a mid knot for 2 knot curves
  if (sKnots.GetSize() == 2)
    {
      double dStep = (sKnots[1]-sKnots[0])/3.0;
      double dNew1 = sKnots[0] + dStep;
      double dNew2 = dNew1 + dStep;
      SER( rRail1Curve.InsertOneKnot(dNew1,1) );
      SER( rRail2Curve.InsertOneKnot(dNew1,1) );
      SER( rRail1Curve.InsertOneKnot(dNew2,1) );
      SER( rRail2Curve.InsertOneKnot(dNew2,1) );
      SER( rRail1Curve.GetKnots(sKnots) );
    }

  ULONG nKnotsOld = sKnots.GetSize();

  // insert a mid interval knot into largest span
  while( sKnots.GetSize() < 4 ) {

      ULONG lMaxSpan = 0;
      double dMaxSpanValue = 0.0;
      for (ULONG k=1; k<sKnots.GetSize(); k++)
        {
          double dSpan = sKnots[k]-sKnots[k-1];
          if (dSpan > dMaxSpanValue)
            {
              lMaxSpan = k;
              dMaxSpanValue = dSpan;
            }
        }

      double dMid = (sKnots[lMaxSpan] + sKnots[lMaxSpan-1])/2.0;
      SER( rRail1Curve.InsertOneKnot(dMid,1) );
      SER( rRail2Curve.InsertOneKnot(dMid,1) );
      SER( rRail1Curve.GetKnots(sKnots) );

      ULONG nKnotsNew = sKnots.GetSize();

      // if no knots were added then bail to avoid endless loop
      if( nKnotsNew == nKnotsOld )
        { return( SM_ERR ); }

      nKnotsOld = nKnotsNew;

    } // end while inserting knots into curves to make at least 4

  // SolidCraft change: where the rails meet at an end of the fillet (a degenerate corner: the
  // two faces touch at its one fillet vertex and the cross section is a point, below), both
  // rails end at that one point. Each rail is lifted from its own surface
  // (SmFilletIntersector::FlushCurve), and at a vertex the two surfaces pass through it only to
  // within the model's accuracy, so their ends differ by that much (6.5e-12 mm on a sketch with
  // solver residues). The point section was the first rail's end while the second rail, a
  // boundary of the skinned surface, kept its own: the surface came out close to singular there
  // instead of exactly (SmSurface test 5). Both rail ends move to their midpoint, within the
  // distance the ends already agree to.
  for ( ULONG lEnd=0; lEnd<2; lEnd++ )
    {
      SmExtent1d sIvl = rRail1Curve.GetNaturalInterval();
      double     dT   = ( lEnd == 0 ) ? sIvl.GetMin() : sIvl.GetMax();
      SmPoint3d  sEnd1, sEnd2;
      SER( rRail1Curve.EvaluatePoint( dT, sEnd1 ) );
      SER( rRail2Curve.EvaluatePoint( dT, sEnd2 ) );
      if ( sEnd1.DistanceBetween( sEnd2 ) >= dThisApproxTol3d )
        { continue; }
      SmPoint3d sMeet = ( sEnd1 + sEnd2 ) / 2.0;
      SmBSplineCurve * apRails[2] = { &rRail1Curve, &rRail2Curve };
      for ( ULONG lRail=0; lRail<2; lRail++ )
        {
          ULONG     lIdx = ( lEnd == 0 ) ? 0 : apRails[lRail]->GetNumberControlPoints() - 1;
          SmPoint3d sCP;
          double    dW = 1.0;
          SER( apRails[lRail]->GetControlPoint( SM_CP_EUCLIDIAN_RATIONAL, lIdx, sCP, dW ) );
          SER( apRails[lRail]->SetControlPoint( SM_CP_EUCLIDIAN_RATIONAL, lIdx, sMeet, dW ) );
        }
    }

  // when center, rail1, and rail2 curves are parallel lines - NewSurf = RuledSurface
  // when center, rail1, and rail2 curves are arcs           - NewSurf = SurfaceOfRevolution
  // else                                                    - NewSurf = Skinned Surface
  SmBoolean bSweep   = FALSE;
  SmBoolean bRevolve = FALSE;

  SmAxis2Placement sRefFrame;
  double dEndAng = 0.0;
  SER(sm_CircleSweptOrRevolved(rCenterLine, rRail1Curve, rRail2Curve, dThisApproxTol3d, bSweep, bRevolve, sRefFrame, dEndAng));

  // track max/min number of knots in cross-section curves
  ULONG lMaxNumUniqueKnots = 0;
  ULONG lMinNumUniqueKnots = 1000000;

  ULONG ii, jj, kkk;

  // for two passes - 1st pass allocate curves
  //                  2nd pass make all curve knot vectors compatible
  for ( kkk=0; kkk<2; kkk++ )
    {
      // for every pRail1Curve->knot
      for ( ii=1; ii<sKnots.GetSize(); ii++ )
        {
          // skip to end cross section curve for sweeps and revolves
          if(   ii > 1
             && (   bSweep
                 || bRevolve))
            { ii = sKnots.GetSize() - 1 ; }

          // locals - center, rail1, and rail2 curve sample points
          SmVector3d sCentPV[2];
          SmPoint3d sP1, sP2;

          // for lNumBetween count locations between knot points
          // ok for lNumBetween == 0
          ULONG lNumBetween = 0;
          for (jj=0; jj<=lNumBetween+1; jj++)
            {
              double dT =   (jj==0) ? sKnots[ii-1]
                          : (jj==lNumBetween+1) ? sKnots[ii]
                          : (sKnots[ii-1] + jj * (sKnots[ii]-sKnots[ii-1]) / (lNumBetween+1.0)) ;

              if (jj==0 && ii != 1) { continue; } // Only do first point first time through
              // Compute the knot index

              ULONG lKnotIndex = (jj==0 && ii == 1)
                                 ? 0
                                 : ii;

              // get center, rail1, rail2 points (and their tangents) for dT
              SmVector3d sRail1PV[2], sRail2PV[2];
              SER(rCenterLine.Evaluate(dT,1,TRUE,sCentPV));
              SER(rRail1Curve.Evaluate(dT,1,TRUE,sRail1PV,FALSE));
              SER(rRail2Curve.Evaluate(dT,1,TRUE,sRail2PV,FALSE));
              sP1 = SmPoint3d(sRail1PV[0]);
              sP2 = SmPoint3d(sRail2PV[0]);

              // SolidCraft change: the section is the circle along which the balls' envelope
              // touches the ball, square to the centre line. Where the ball's radius varies, its
              // centre is the ball's centre moved along the centre line by -r r'/|c'| (the
              // envelope's condition <p - c, c'> = -r r' on the ball |p - c| = r), r and r' read
              // off the rails, which lie on the ball: at the ball's touch points that is their
              // projection on the centre line, and where a corner's rails only lie on the ball, its
              // radius law's own slope. At a constant radius it is zero, the great circle (whose
              // rails, where a corner moved them, would only add their own error).
              SmVector3d sV1 = sP1 - sCentPV[0];
              SmVector3d sV2 = sP2 - sCentPV[0];
              SmVector3d sAlong = sCentPV[1];
              double dSpeed = sAlong.Length();
              if (bBallRadiusVaries && dSpeed > SM_EFF_ZERO && sV1.Length() > SM_EFF_ZERO && sV2.Length() > SM_EFF_ZERO)
                {
                  double dR1 = sV1.Length(), dR2 = sV2.Length();
                  double dRate = 0.5 * (sV1.Dot(sRail1PV[1] - sCentPV[1]) / dR1 + sV2.Dot(sRail2PV[1] - sCentPV[1]) / dR2);
                  sAlong = sAlong / dSpeed;
                  sCentPV[0] += sAlong * (-0.5 * (dR1 + dR2) * dRate / dSpeed);
                  sV1 = sP1 - sCentPV[0];
                  sV2 = sP2 - sCentPV[0];
                }
              SmBSplineCurve * pArc = NULL;

#ifdef SM_DEBUG_CODE
              if (bDebugMe)
                {
                  sm_GraphicsLoop();
                  smgfx_SetLook(1,5, 1,0,1); sCentPV[0].Draw(); sm_GraphicsLoop();
                  smgfx_SetLook(1,5, 0,0,1); sP1.Draw(); sm_GraphicsLoop();
                  smgfx_SetLook(1,5, 0,1,1); sP2.Draw(); sm_GraphicsLoop();
                  smgfx_SetLook(1,5, 1,0,0); sCentPV[0].DrawPointToPoint( sP1 ); sm_GraphicsLoop();
                  smgfx_SetLook(1,5, 0,0,1); sCentPV[0].DrawPointToPoint( sP2 ); sm_GraphicsLoop();
                  sm_GraphicsLoop();
                }
#endif

              // when railCurve samplePoints are close together - make a degenerate cross section curve
              if (sP1.DistanceBetween(sP2) < dThisApproxTol3d)
                {
                  // On the 1st pass
                  if (kkk==0)
                    {
                      // Create Create Degenerate Point curve for each dT value on Rail1
                      SER(SmBSplineCurve::CreatePointCurve(crContext,sP1,pArc));  /* parameterized from 0 to 1 */

                      // track max/min cross section unique knot count
                      ULONG lNumCtrlPts  = pArc->GetNumberOfUniqueKnots();
                      lMaxNumUniqueKnots = smos_Max(lMaxNumUniqueKnots,lNumCtrlPts);
                      lMinNumUniqueKnots = smos_Min(lMinNumUniqueKnots,lNumCtrlPts);

                      // accumulate cross section curve and parameter value
                      sParametersOfCrossSections.Add(dT);
                      sCrossSections.Add(pArc);
                    }
                  else // on the 2nd pass
                    {
                      // Insert knots into degenerate Curve until UniqueKnotCount = MaxNumUniqueKnots
                      pArc               = sCrossSections[lKnotIndex];
                      ULONG lNumKnots    = lMaxNumUniqueKnots;
                      ULONG lNumNewKnots = lNumKnots - 2;
                      if (lNumNewKnots == 0) continue;
                      SmTArray<double> sNewKnots;
                      ULONG mm;
                      for (mm=0; mm<lNumNewKnots; mm++)
                        {
                          double dParam = (mm+1.0) / (lNumNewKnots+1.0);
                          sNewKnots.Add(pArc->GetNaturalInterval().Evaluate(dParam));
                        }
                      SER(pArc->InsertKnots(sNewKnots));
                    }
                  continue;
                } // end sameLength gapVector check

              // else the the railCurve points are far appart.
              // build and add a circular arc to the cross section array.
              double dRadius1 = sV1.Length();
              double dRadius2 = sV2.Length();

              // Check for case of zero radius on interior of surface
              // We do allow it on the ends of the surface - however it
              // may not be allowed in the creation - we'll see.
              if (   dRadius1  < SM_EFF_ZERO || dRadius2 < SM_EFF_ZERO)
                {
                  // When sample point with zero radius is an internal point
                  if (ii > 0 && ii < sKnots.GetSize()-1)
                       { // signal an error
                         SER_MSG(SM_ERR, _T("Zero radius cross-section curve")) ;
                       }
                  else { SER_MSG(SM_ERR, _T("Zero radius cross-section curve")) ; // Someday we will handle this at the ends
                       }
                } // end zero radius check

              // else if the radii differences are large
              else if (smos_Fabs(dRadius1 - dRadius2) > dRadius1 / 10.0)
                {
                  // skip adding a cross section for this parameter value
                  SM_DBG_WARN(_T("skipped adding a circular cross-section to a fillet surface because \n") 
                              _T("  railCurve points are not on a circle centered on the centerCurve point") ) ;

                  // But put something into the arrays (on 1st pass), to keep the indexing correct. [B663]
                  if ( kkk==0 )
                  {
                      sParametersOfCrossSections.Add( dT );
                      sCrossSections.Add( NULL );
                  }

                  continue;
                }

              // get the CCW angle from sV1 to sV2
              double dAngle;
              SER(sV1.Unitize());
              SER(sV2.Unitize());
              SmVector3d sZVec = sV1 * sV2;
              SER(sZVec.Unitize());
              SER(sZVec.CCWAngleBetween(sV1,sV2,dAngle));

              // quit for zero and 180 angles
              if (   smos_Fabs(dAngle) < SM_EFF_ZERO_SQRT          // Degenerate case here
                  || smos_Fabs(dAngle) > SM_PI - SM_EFF_ZERO_SQRT) // Nearly 180 degree case here
                { SER_MSG(SM_ERR, _T("Tangent or Degenerate Angle for cross-section curve")); }

              // negate normal for neg angles and recompute CCW angle
              if (dAngle < 0.0) { sZVec = - sZVec; }
              SER(sZVec.CCWAngleBetween(sV1,sV2,dAngle));
              if (dAngle < SM_EFF_ZERO) { SER(SM_ERR); }  // Can't happen ??

              // make unitNormal coordinate system
              SmVector3d sYVec = sZVec * sV1;
              SmVector3d sXVec = sYVec * sZVec;

              // when asked - make complementary circle
              if (m_sArcs.bComplement)
                {
                  sYVec = - sYVec;
                  dAngle = 2.0 * SM_PI - dAngle;
                }

              // make coordinate system object
              SER(sXVec.Unitize());
              SER(sYVec.Unitize());
              SmAxis2Placement sPlacement;
              SER(sPlacement.SetCanonical(sCentPV[0],sXVec,sYVec));

              // when approximating arcs
              if (m_sArcs.bApproximateArcs)
                {
                  // on the 1st pass - make an arc, track the min/max number of Unique knots
                  if (kkk==0)
                    {
                      SER(SmBSplineCurve::ApproximateArc(crContext,
                                                         3, NULL, sPlacement,
                                                         dRadius1, 0.0,
                                                         dAngle*180.0/SM_PI,
                                                         m_sArcs.dTolPercentRadius,
                                                         pArc));

                      // track max/min number of unique knots
                      ULONG lNumKnots = pArc->GetNumberOfUniqueKnots();
                      lMaxNumUniqueKnots = smos_Max( lMaxNumUniqueKnots, lNumKnots );
                      lMinNumUniqueKnots = smos_Min( lMinNumUniqueKnots, lNumKnots );
                    }
                  // on 2nd pass with equal gapVector radii
                  else if (smos_Fabs(dRadius1 - dRadius2) < dThisApproxTol3d)
                    {
                      // get current arc and uniqueKnot count
                      pArc = sCrossSections[lKnotIndex];
                      if ( pArc == NULL )
                        { continue; }

                      ULONG lNumUniqueKnots = pArc->GetNumberOfUniqueKnots();

                      // when needed - replace pArc with a new Arc with appropriate number of knots
                      if (lNumUniqueKnots != lMaxNumUniqueKnots)
                        {
                          ULONG lMaxCtrlPts = lMaxNumUniqueKnots + 2;
                          SER(SmBSplineCurve::ApproximateArc(crContext,
                                                             3, &lMaxCtrlPts, sPlacement,
                                                             dRadius1, 0.0,
                                                             dAngle*180.0/SM_PI,
                                                             m_sArcs.dTolPercentRadius,
                                                             pArc));
                          SM_ASSERT(sCrossSections[lKnotIndex] != NULL) ; delete sCrossSections[lKnotIndex] ; sCrossSections[lKnotIndex] = NULL ;
                          sCrossSections[lKnotIndex] = pArc;
                        }
                    }

                  else // on 2nd pass with unequal gapVector radii
                    { // Do need to fix knot vector of blend
                      pArc               = sCrossSections[lKnotIndex];
                      if ( pArc == NULL )
                        { continue; }
                      ULONG lNumKnots    = lMaxNumUniqueKnots;
                      ULONG lNumNewKnots = lNumKnots - 2;
                      if (lNumNewKnots == 0) continue;
                      SmTArray<double> sNewKnots;
                      ULONG mm;
                      for (mm=0; mm<lNumNewKnots; mm++)
                        {
                          double dParam = (mm+1.0) / (lNumNewKnots+1.0);
                          sNewKnots.Add(pArc->GetNaturalInterval().Evaluate(dParam));
                        }
                      SER(pArc->InsertKnots(sNewKnots));
                      continue;
                    }
                } // end m_sArcs.bApproximateArcs == TRUE check
              else // not approximating arcs
                {
                  // 1st pass build a circular segment
                  if (kkk == 0)
                    { SER(SmBSplineCurve::CreateCircleSegment(crContext, 3, sPlacement,
                                                              dRadius1, 0.0,
                                                              dAngle*180.0/SM_PI,
                                                              SM_CO_QUADRATIC,
                                                              pArc));
                    }
                }

#ifdef SM_DEBUG_CODE
              if (bDebugMe)
                {
                  sm_GraphicsLoop();
                  smgfx_SetLook(2,3, 0,1,1); pArc->DrawWDeriv(pArc->GetNaturalInterval(),0); sm_GraphicsLoop();
                  sm_GraphicsLoop();
                }
#endif
              // on 1st pass - See if we have a non-circle case
              if (kkk == 0 && smos_Fabs(dRadius1 - dRadius2) > 1.0e-8)
                {
                  // Tweak arc or create a blend curve using railPoints and current pArc tangents
                  // to interpolate the railCurve sample points
                  SmTArray<SmPoint3d> sPoints;
                  SmTArray<SmVector3d> sVectors;
                  SmPoint3d sPV[2];

                  // add railCurve1Pos/pArcStartTangent to point/vector arrays
                  SER(pArc->Evaluate(pArc->GetNaturalInterval().GetMin(),1,TRUE,sPV));
                  sPoints.Add(sP1);
                  sVectors.Add(sPV[1]);

                  // add railCurve2Pos/pArcEndTangent to point/vector arrays
                  SER(pArc->Evaluate(pArc->GetNaturalInterval().GetMax(),1,TRUE,sPV));
                  sPoints.Add(sP2);
                  sVectors.Add(sPV[1]);

                  // Here either tweak the arc a little by moving the ends or
                  // create an interpolating curve.
                  if (smos_Fabs(dRadius1 - dRadius2) < dThisApproxTol3d*4.0)
                    {
                      // move the arc endPoints to lie directly on the rail curves
                      // gwc: this step is the construction technique for this function
                      //      it is not an attempt to close an arbitrary gap - don't output debug checks
                      SER(pArc->EditEndPoint(sP1, TRUE,  NULL, SM_EE_CONSTRUCTED_END_POINT, FALSE));
                      SER(pArc->EditEndPoint(sP2, FALSE, NULL, SM_EE_CONSTRUCTED_END_POINT, FALSE));
                    }
                  else // replace the circular arc with an interpolating curve
                    {
                      SM_ASSERT(pArc != NULL) ; delete pArc ; pArc = NULL ;
                      SER(SmBSplineCurve::FairBlendCurveDerivatives(sPoints,sVectors,NULL));
                      SER(SmBSplineCurve::CreateInterpolatingCurve(crContext,SM_CP_UNIFORM,3,
                          3,sPoints,sVectors,NULL,TRUE,pArc));
                    }
#ifdef SM_DEBUG_CODE
SmBoolean bDebugMe4 = FALSE;
                  if (bDebugMe4)
                    {
                      //                    smgfx_Erase();
                      smgfx_SetColor(0,0,1);
                      pArc->DrawWDeriv(pArc->GetNaturalInterval(),0);
                      sm_GraphicsLoop();
                      smgfx_SetColor(1,0,0);
                      pArc->DrawPolygon();
                      sm_GraphicsLoop();
                    }
#endif
                } // end need to patch up cross-sections to interpolate rail curves check


              // When asked - mirror the cross section
              if (m_sArcs.bMirror)
                {
                  SmAxis2Placement sMirror;
                  SmVector3d sMirrorXAxis = sP2 - sP1;
                  SER(sMirrorXAxis.Unitize());
                  SER(sZVec.Unitize());
                  SER(sMirror.SetCanonical(sP1,sMirrorXAxis,sZVec));
                  SmCurve *pMirror = NULL;
                  SER(pArc->CreateMirrorCurve(crContext,sMirror,pMirror));
                  NER(pMirror);
#ifdef SM_DEBUG_CODE
SmBoolean bDebugMe2 = FALSE;
                  if (bDebugMe2)
                    {
                      smgfx_Erase();
                      smgfx_SetColor(0,0,1);
                      pArc->DrawWDeriv(pArc->GetNaturalInterval(),0);
                      sm_GraphicsLoop();
                      smgfx_SetColor(1,0,0);
                      smgfx_SetLineWidth(4.0);
                      pMirror->DrawWDeriv(pMirror->GetNaturalInterval(),0);
                      smgfx_SetLineWidth(2.0);
                      sm_GraphicsLoop();
                    }
#endif

                  SM_ASSERT(pArc != NULL) ; delete pArc ; pArc = NULL ;
                  pArc = SM_CAST_PTR(SmBSplineCurve,pMirror); NER(pArc);
                  if (kkk == 1)
                    {
                      sCrossSections[lKnotIndex] = pArc;
                    }

                } // end mirror cross section curve check

              // accumulate the cross section curves
              if (kkk == 0)
                {
                  sParametersOfCrossSections.Add(dT);
                  sCrossSections.Add(pArc);
                }
            } // end iter every between knot cross section curve (currently none)
        } // end iter every pRail1Curve->Knot

      // if (! (m_sArcs.bApproximateArcs && lMaxNumUniqueKnots != lMinNumUniqueKnots))
      if (   m_sArcs.bApproximateArcs == FALSE
          || lMaxNumUniqueKnots == lMinNumUniqueKnots)
        {
          // only need 1 pass
          break;
        }
    } // end iter 2 times

  // arrive here when all cross section curves and their params are stored in
  // sCrossSections and sParametersOfCrossSections

  // In case we put blanks in:  (Can't use CompressZeros() on <double> array.)
  for ( ii=0; ii<sCrossSections.GetSize(); ii++ )
  {
      if ( sCrossSections[ii] == NULL )
      {
          sCrossSections.RemoveAt( ii );
          sParametersOfCrossSections.RemoveAt( ii );
          ii--;
      }
  }



  SmBSplineSurface *pFilletSurface = NULL;

  // when surface can be constructed as a ruled surface (linear center and rail curves)
  if (bSweep)
    {
      if ( sCrossSections.GetSize() < 2 )  // [B448]
        { SER( SM_ERR ); }

      SmBSplineCurve *pStart = sCrossSections[0];
      SmBSplineCurve *pEnd   = sCrossSections.GetLast();
      SER(SmBSplineSurface::CreateRuledSurface(crContext, *pStart, *pEnd,
                                               SM_SP_U,pFilletSurface));
      SmExtent1d sIvl = rRail1Curve.GetNaturalInterval();
      SmExtent2d sDomain(SmPoint2d(sIvl.GetMin(),0.0),
                         SmPoint2d(sIvl.GetMax(),1.0));
      SER(pFilletSurface->Reparameterize(sDomain));

    } // end ruled surface case
  else if (bRevolve) // when surface can be constructed as a surfaceOfRevolution
    {
      if ( sCrossSections.GetSize() < 1 )  // [B448]
        { SER( SM_ERR ); }

      SmBSplineCurve *pGenCurve = sCrossSections[0];
      SER(SmBSplineSurface::CreateSurfOfRevolution(crContext, pGenCurve,
                                                   sRefFrame.GetOriginRef(),
                                                   sRefFrame.GetZAxis(),
                                                   dEndAng,
                                                   pFilletSurface));
      SmExtent1d sIvl = rRail1Curve.GetNaturalInterval();
      SmExtent2d sDomain(SmPoint2d(sIvl.GetMin(),0.0),
                         SmPoint2d(sIvl.GetMax(),1.0));
      SER(pFilletSurface->Reparameterize(sDomain));
    }
  else // build a skin surface to interpolate all crossSection curves and option end crossTangent values
    {
      if ( sCrossSections.GetSize() < 2 )  // [B448]
        { SER( SM_ERR ); }

      SmBSplineSurface * pDerivSurfs[2] = { NULL, NULL };

      // when asked for G1 continuity
      if (bNeedG1Continuity)
        {
          SmExtent1d sIvl = rCenterLine.GetNaturalInterval();

          // for NewSurface start and end boundaries
          for (kkk=0; kkk<2; kkk++)
            {
              // We need to provide cross-boundary derivative surface
              // for skinned-surface creation.

              // Create rulled surface using end cross-section
              // curve and the tangent vectors of rails
              SmBSplineSurface * pSurf     = NULL;
              SmBSplineCurve   * pNewCurve = NULL;
              SmBSplineCurve   * pGenCurve = sCrossSections[0];
              double dT = sIvl.GetMin();
              if (kkk == 1)
                {
                  pGenCurve = sCrossSections.GetLast();
                  dT = sIvl.GetMax();
                }

              // get railCurve endTangents
              SmVector3d sPVSt[2], sPVEnd[2];
              SER(rRail1Curve.Evaluate(dT,1,TRUE,sPVSt));
              // gwc:replaced one line to interpolate rail2Curve rather than CenterLine
              SER(rRail2Curve.Evaluate(dT,1,TRUE,sPVEnd));
              // SER(rCenterLine.Evaluate(dT,1,TRUE,sPVEnd));

              // make a curve that sweeps from rail1 to rail2 endTangent values
              SER(pGenCurve->CreateByScaleTransRot(crContext,sPVSt[1],sPVEnd[1],
                                                   0,pNewCurve));
              SmObjDelete sCleanNewCrv(pNewCurve);

              // make a Ruled Surface from the crossTangent Curve
              SER(SmBSplineSurface::CreateRuledSurface(crContext,*pGenCurve,*pNewCurve,SM_SP_U,pSurf));
#ifdef SM_DEBUG_CODE
//                    pSurf->Dump();
#endif
              SmExtent1d sCrvIvl = pGenCurve->GetNaturalInterval();
              SmPoint2d  sUVMin  = SmPoint2d(sCrvIvl.GetMin(),0.0);
              SmPoint2d  sUVMax  = SmPoint2d(sCrvIvl.GetMax(),1.0);
              SER( pSurf->AdjustSTEPUVDomain(SmExtent2d(sUVMin,sUVMax)) );
#ifdef SM_DEBUG_CODE
//                    pSurf->Dump();
#endif
              // save the crossTangent surface
              pDerivSurfs[kkk] = pSurf;

            } // end iter NewSurface start/ends to build crossTangent surfaces

          // add a new knot in the NewSurface start/end intervals
          // to allow the newSurface to interpolate the computed crossTangent values
          SmTArray<double> sKnts;
          SER( rRail1Curve.GetKnots(sKnts) );

          double dStartKnotMid = (sKnts[0] + sKnts[1]) / 2.0;
          SER( rRail1Curve.InsertOneKnot(dStartKnotMid,1) );
          SER( rRail2Curve.InsertOneKnot(dStartKnotMid,1) );

          double dEndKnotMid = (sKnts[sKnts.GetSize()-2] + sKnts.GetLast()) / 2.0;
          SER( rRail1Curve.InsertOneKnot(dEndKnotMid,1) );
          SER( rRail2Curve.InsertOneKnot(dEndKnotMid,1) );

        } // end need G1 continuity check

      // Build the Skinned Surface to interpolate crossSections and optional end Cross Tangents
      SER(SmBSplineSurface::CreateSkinnedSurface(crContext,
                                                sCrossSections,
                                                FALSE, SM_SP_U, dThisApproxTol3d,
                                                &rRail1Curve, &rRail2Curve,
                                                FALSE, &sParametersOfCrossSections,
                                                pDerivSurfs,
                                                pFilletSurface));
      // clean up the crossTangent surfaces
      for (kkk=0; kkk<2; kkk++)
        {
          if (pDerivSurfs[kkk])  { delete pDerivSurfs[kkk]; pDerivSurfs[kkk] = NULL ; }
        }


      // We have a case where the rail curves start off in the same direction,
      // because the surfaces are not only tangent at one end but also both locally
      // flat.  This means that the end control points of the first interior row of
      // control points will be identical.  However, the algorithm above creates
      // internal points on that first internal row that are not the same.  We fix this
      // by making them all the same.  This means that all of the u-derivatives at the
      // singularity will be parallel: that's legal and doesn't cause problems. [B597]
      ULONG lSing;  // Looking for UMIN or UMAX, or both.
      pFilletSurface->GetSingularities( NULL, NULL, &lSing );
      if ( lSing & SM_SS_UMIN )
        {
          // ULONG lNumU = pFilletSurface->GetNumberControlPoints( SM_SP_U );
          ULONG lNumV = pFilletSurface->GetNumberControlPoints( SM_SP_V );
          SmPoint3d sPt1, sPt2;
          double dW1, dW2;
          SmZoneTol3d dTol = SmTol::GetZoneTol3d( pFilletSurface );
          ULONG lIdxU = 1;  // First interior row.
          SER( pFilletSurface->GetControlPoint( SM_CP_EUCLIDIAN_RATIONAL, lIdxU,       0, sPt1, dW1 ) );
          SER( pFilletSurface->GetControlPoint( SM_CP_EUCLIDIAN_RATIONAL, lIdxU, lNumV-1, sPt2, dW2 ) );
          if ( sPt2.CloserThan( dTol, sPt1 ) )
            {
              sPt1 = ( sPt1 + sPt2 ) / 2;  // Make sure they're exactly the same.
              dW1  = ( dW1  + dW2  ) / 2;
              for ( ii=0; ii<lNumV; ii++ )
                { SER( pFilletSurface->SetControlPoint( SM_CP_EUCLIDIAN_RATIONAL, lIdxU, ii, sPt1, dW1 ) ); }
            }
        }
      if ( lSing & SM_SS_UMAX )
        {
          ULONG lNumU = pFilletSurface->GetNumberControlPoints( SM_SP_U );
          ULONG lNumV = pFilletSurface->GetNumberControlPoints( SM_SP_V );
          SmPoint3d sPt1, sPt2;
          double dW1, dW2;
          SmZoneTol3d dTol = SmTol::GetZoneTol3d( pFilletSurface );
          ULONG lIdxU = lNumU-2;  // First interior row.
          SER( pFilletSurface->GetControlPoint( SM_CP_EUCLIDIAN_RATIONAL, lIdxU,       0, sPt1, dW1 ) );
          SER( pFilletSurface->GetControlPoint( SM_CP_EUCLIDIAN_RATIONAL, lIdxU, lNumV-1, sPt2, dW2 ) );
          if ( sPt2.CloserThan( dTol, sPt1 ) )
            {
              sPt1 = ( sPt1 + sPt2 ) / 2;  // Make sure they're exactly the same.
              dW1  = ( dW1  + dW2  ) / 2;
              for ( ii=0; ii<lNumV; ii++ )
                { SER( pFilletSurface->SetControlPoint( SM_CP_EUCLIDIAN_RATIONAL, lIdxU, ii, sPt1, dW1 ) ); }
            }
        }

    } // end need to build a skin surface branch

  // set output
  rpFilletSurface = pFilletSurface;

#ifdef SM_DEBUG_CODE
SmBoolean bDebugMe2=FALSE;
  if (bDebugMe2)
    {
      smgfx_Erase();
      smgfx_SetLook(1,2, 1,0,0);
      rRail1Curve.DrawWDeriv(rRail1Curve.GetNaturalInterval()); sm_GraphicsLoop();
      rRail2Curve.DrawWDeriv(rRail2Curve.GetNaturalInterval()); sm_GraphicsLoop();
      smgfx_SetLook(1,2, 0,0,1);
      rCenterLine.DrawWDeriv(rCenterLine.GetNaturalInterval()); sm_GraphicsLoop();
      sm_GraphicsLoop();
      pFilletSurface->Dump();
      smgfx_SetLook(1,2, 0,1,1); pFilletSurface->DrawUV(4,4); sm_GraphicsLoop();
      sm_GraphicsLoop();
    }
#endif

  // all done
  return SM_SUCCESS;

} // end SmGnSectionFSG::CircleSurfaceFromCurves

/*******************************************************************//**
PURPOSE: The blend surface between the solver's rails. The solver's faces
    are its offset surfaces; the section is built on the faces under them.
    The rails and their uv curves are rebuilt on the section's stations from
    the same exact feet, so that the surface's boundaries are the rails.
***********************************************************************/
SmStatus SmGnSectionFSG::CreateSurface
  (SmFilletSolver              & rFilletSolver,
   SmBSplineCurve              * pCenterLine,
   SmBSplineCurve              * pRail1Curve,
   SmBSplineCurve              * pRail2Curve,
   const SmTArray<SmTsectPnt*> & crFilletPoints,
   SmSurface                   * pSurf1,
   SmBSplineCurve              * pUV1,
   SmSurface                   * pSurf2,
   SmBSplineCurve              * pUV2,
   SmBSplineSurface           *& rpFilletSurface)
{
  NER(pCenterLine); NER(pRail1Curve); NER(pRail2Curve); NER(pSurf1); NER(pSurf2);
  if (m_lOrder == 0)
    {
      SER(SmBSplineSurface::CreateRuledSurface(rFilletSolver.GetCreationContext(), *pRail1Curve, *pRail2Curve, SM_SP_V, rpFilletSurface));
      return SM_SUCCESS;
    }
  if (m_lOrder == 1)
    {
      SER(CircleSurface(rFilletSolver, crFilletPoints, pCenterLine, pRail1Curve, pRail2Curve, pSurf1, pUV1, pSurf2, pUV2, rpFilletSurface));
      return SM_SUCCESS;
    }
  const SmOffsetSurface * apOff[2] = { SM_CAST_PTR(SmOffsetSurface, pSurf1), SM_CAST_PTR(SmOffsetSurface, pSurf2) };
  SmSurface * apFace[2] = { apOff[0] != NULL ? apOff[0]->GetBaseSurface() : pSurf1,
                            apOff[1] != NULL ? apOff[1]->GetBaseSurface() : pSurf2 };
  NER(apFace[0]); NER(apFace[1]);
  SER(Build(rFilletSolver.GetCreationContext(), *pCenterLine, *pRail1Curve, *pRail2Curve, rFilletSolver.GetThisApproxTol3d(),
            apFace[0], pUV1, apFace[1], pUV2, TRUE, rpFilletSurface));
  return SM_SUCCESS;
} // end SmGnSectionFSG::CreateSurface

/*******************************************************************//**
PURPOSE: The blend surface between two rails on the faces' surfaces pSurf1
    and pSurf2 (their uv curves pUV1, pUV2, or NULL), about the centre line.
    The rails are not changed.
***********************************************************************/
SmStatus SmGnSectionFSG::CreateSurfaceFromCurves
  (const SmContext   & crContext,
   SmBSplineCurve    & rCenterLine,
   SmBSplineCurve    & rRail1Curve,
   SmBSplineCurve    & rRail2Curve,
   double              dThisApproxTol3d,
   SmBoolean           bNeedG1Continuity,
   SmSurface         * pSurf1,
   SmBSplineCurve    * pUV1,
   SmSurface         * pSurf2,
   SmBSplineCurve    * pUV2,
   SmBSplineSurface *& rpFilletSurface)
{
  if (m_lOrder == 0)
    {
      SER(SmBSplineSurface::CreateRuledSurface(crContext, rRail1Curve, rRail2Curve, SM_SP_V, rpFilletSurface));
      return SM_SUCCESS;
    }
  if (m_lOrder == 1)
    {
      // the curves of a ball of constant radius (an offset's edge, a corner's ease)
      SER(CircleSurfaceFromCurves(crContext, rCenterLine, rRail1Curve, rRail2Curve, dThisApproxTol3d, bNeedG1Continuity,
                                  FALSE, pSurf1, pUV1, pSurf2, pUV2, rpFilletSurface));
      return SM_SUCCESS;
    }
  SER(Build(crContext, rCenterLine, rRail1Curve, rRail2Curve, dThisApproxTol3d, pSurf1, pUV1, pSurf2, pUV2, FALSE, rpFilletSurface));
  return SM_SUCCESS;
} // end SmGnSectionFSG::CreateSurfaceFromCurves

/*******************************************************************//**
PURPOSE: The G^n blend surface (file header) between the rails, its stations
    the knots of the rails and the centre line. bReplaceRails: the rails and
    uv curves become the jets' own order-0 curves (the surface's boundaries).

NOTES: The surface's u is the rails' parameter, v the section's t, from the
    first rail (v = 0) to the second (v = 1).
***********************************************************************/
SmStatus SmGnSectionFSG::Build
  (const SmContext   & crContext,
   SmBSplineCurve    & rCenterLine,
   SmBSplineCurve    & rRail1Curve,
   SmBSplineCurve    & rRail2Curve,
   double              dTol,
   SmSurface         * pFace1,
   SmBSplineCurve    * pUV1,
   SmSurface         * pFace2,
   SmBSplineCurve    * pUV2,
   SmBoolean           bReplaceRails,
   SmBSplineSurface *& rpFilletSurface,
   double              dRadius,
   SmBoolean         * pbNotTheFeet,
   const SmGnLawRound * pLawRound)
{
  rpFilletSurface = NULL;
  if (pbNotTheFeet != NULL) { *pbNotTheFeet = FALSE; }
  NER(pFace1); NER(pFace2); NER(pUV1); NER(pUV2);
  // a ball whose radius varies (the ball's circle only): its law r(t) and the curve whose planes square
  // to it at t hold the ball's centre (the solver's), t the stations' parameter
  if (pLawRound != NULL && (m_lOrder != 1 || pLawRound->pLeader == NULL || pLawRound->pLaw == NULL || pLawRound->pPoints == NULL))
    { SER(SM_ERR_INVALID_INPUT); }
  sm_GnLaw sLaw;
  const sm_GnLaw * pGnLaw = NULL;
  if (pLawRound != NULL)
    {
      sLaw.pCurve = pLawRound->pLaw;
      SER(sm_NurbsCurveOf(*pLawRound->pLaw, sLaw.sN));
      pGnLaw = &sLaw;
    }
  const SmBSplineCurve & rLeader = (pLawRound != NULL) ? *pLawRound->pLeader : rCenterLine;
  const ULONG n = m_lOrder;
  if (!(n == 1 && !m_sArcs.bApproximateArcs && !m_sArcs.bMirror && !m_sArcs.bComplement) && !IsValid(n, m_dThumbweight))
    { SER_MSG(SM_ERR_INVALID_INPUT, _T("SmGnSectionFSG: the order must be 1 (the ball's circle) or 2 to 4 with a thumbweight below 1")); }
  const SmSurface * apFace[2] = { pFace1, pFace2 };
  SmBSplineCurve * apRail[2] = { &rRail1Curve, &rRail2Curve };
  SmBSplineCurve * apUV[2] = { pUV1, pUV2 };
  // the stations: the unique knots of the rails and the centre line; with a law, the solver's own
  // stations' t (their feet the seeds) and the leader's and the law's knots between
  SmTArray<double> sKnots;
  std::vector<double> sS;
  std::vector< std::array<double, 4> > sLawSeeds; // with a law: each station's seed
  std::vector<SmBoolean> sLawFoot;                // and whether it is one of the solver's stations
  std::vector<double> sLawCp;                     // and then its parameter on the solver's rails
  const SmBSplineCurve * apCurves[3] = { &rCenterLine, &rRail1Curve, &rRail2Curve };
  SmExtent1d sIvl = rCenterLine.GetNaturalInterval();
  if (pLawRound != NULL)
    {
      // each solver's station: its t, and its feet as the solver's uv curves hold them at its curve
      // parameter (the rails' and uv curves' own, which is not t)
      std::vector< std::pair<double, std::array<double, 5> > > sPts;
      for (ULONG i=0; i<pLawRound->pPoints->GetSize(); i++)
        {
          SmTsectPnt * pPnt = (*pLawRound->pPoints)[i];
          NER(pPnt);
          SmPoint3d asW[2];
          for (int f=0; f<2; f++) { SER(apUV[f]->EvaluatePoint(pPnt->m_dCurveParameter, asW[f])); }
          const std::array<double, 5> sX = { asW[0].x, asW[0].y, asW[1].x, asW[1].y, pPnt->m_dCurveParameter };
          sPts.push_back(std::make_pair(pPnt->m_adUserDoubles[0], sX));
        }
      std::sort(sPts.begin(), sPts.end(), [](const std::pair<double, std::array<double, 5> > & a, const std::pair<double, std::array<double, 5> > & b)
                { return a.first < b.first; });
      if (sPts.size() < 2 || !(sPts.back().first > sPts.front().first))
        { SER_MSG(SM_ERR, _T("SmGnSectionFSG: the round's stations have no length")); }
      sIvl = SmExtent1d(sPts.front().first, sPts.back().first);
      const double dLawScale = smos_Max(1.0, sIvl.GetLength());
      std::vector<double> sAt;
      for (size_t i=0; i<sPts.size(); i++) { sAt.push_back(sPts[i].first); }
      const SmBSplineCurve * apLawCurves[2] = { pLawRound->pLeader, pLawRound->pLaw };
      for (int c=0; c<2; c++)
        {
          SER(apLawCurves[c]->GetKnots(sKnots));
          for (ULONG j=0; j<sKnots.GetSize(); j++)
            if (sKnots[j] > sIvl.GetMin() + 1.0e-9 * dLawScale && sKnots[j] < sIvl.GetMax() - 1.0e-9 * dLawScale)
              { sAt.push_back(sKnots[j]); }
        }
      std::sort(sAt.begin(), sAt.end());
      for (size_t j=0; j<sAt.size(); j++)
        {
          if (!sS.empty() && sAt[j] - sS.back() <= 1.0e-9 * dLawScale)
            { continue; }
          // the seed: the solver's station at t, or between the two about it
          size_t k = 0;
          while (k + 1 < sPts.size() && sPts[k+1].first <= sAt[j]) { k++; }
          std::array<double, 5> sX = sPts[k].second;
          SmBoolean bFoot = (smos_Fabs(sPts[k].first - sAt[j]) <= 1.0e-9 * dLawScale) ? TRUE : FALSE;
          if (!bFoot && k + 1 < sPts.size())
            {
              const double w = (sAt[j] - sPts[k].first) / (sPts[k+1].first - sPts[k].first);
              for (int q=0; q<5; q++) { sX[static_cast<size_t>(q)] = (1.0 - w) * sPts[k].second[static_cast<size_t>(q)] + w * sPts[k+1].second[static_cast<size_t>(q)]; }
            }
          sS.push_back(sAt[j]);
          sLawSeeds.push_back({ sX[0], sX[1], sX[2], sX[3] });
          sLawFoot.push_back(bFoot);
          sLawCp.push_back(sX[4]);
        }
    }
  const double dScale = smos_Max(1.0, sIvl.GetLength());
  // the equal area thumbweight is resolved per station, from its own wedge
  const SmBoolean bEqualArea = (m_dThumbweight == SM_GN_THUMBWEIGHT_EQUAL_AREA) ? TRUE : FALSE;
  for (int c=0; c<3 && pLawRound == NULL; c++)
    {
      const SmExtent1d sCIvl = apCurves[c]->GetNaturalInterval();
      if (   smos_Fabs(sCIvl.GetMin() - sIvl.GetMin()) > 1.0e-9 * dScale
          || smos_Fabs(sCIvl.GetMax() - sIvl.GetMax()) > 1.0e-9 * dScale)
        { SER_MSG(SM_ERR, _T("SmGnSectionFSG: the rails and the centre line run over different intervals")); }
      SER(apCurves[c]->GetKnots(sKnots));
      for (ULONG j=0; j<sKnots.GetSize(); j++) { sS.push_back(sKnots[j]); }
    }
  if (pLawRound == NULL)
  {
    std::sort(sS.begin(), sS.end());
    std::vector<double> sUnique;
    for (size_t j=0; j<sS.size(); j++)
      if (sUnique.empty() || sS[j] - sUnique.back() > 1.0e-9 * dScale)
        { sUnique.push_back(sS[j]); }
    sUnique.front() = sIvl.GetMin();
    sUnique.back() = sIvl.GetMax();
    sS.swap(sUnique);
  }
  if (sS.size() < 2)
    { SER_MSG(SM_ERR, _T("SmGnSectionFSG: the rails have no length")); }

  // the offsets and the sides, at the first station: the centre line there is the ball's centre
  double adR[2], adSide[2];
  {
    // (at the curves' own start: with a law their parameter is not t)
    const double dStart = rCenterLine.GetNaturalInterval().GetMin();
    SmPoint3d sM;
    SER(rCenterLine.EvaluatePoint(dStart, sM));
    double dLawR = 0.0;
    if (pLawRound != NULL)
      {
        SmPoint3d sR;
        SER(pLawRound->pLaw->EvaluatePoint(sIvl.GetMin(), sR));
        dLawR = sR.x;
      }
    for (int f=0; f<2; f++)
      {
        SmPoint3d sRho, sW;
        SER(apRail[f]->EvaluatePoint(dStart, sRho));
        SER(apUV[f]->EvaluatePoint(dStart, sW));
        const SmPoint2d sUV(sW.x, sW.y);
        SmVector3d sN;
        SER(apFace[f]->EvaluateNormal(sUV, TRUE, TRUE, sN));
        SER(sN.Unitize());
        const SmVector3d sOff = sM - sRho;
        // the ball's radius where the solver gives it (its offsets'), else as the rails stand off
        adR[f] = (pLawRound != NULL) ? dLawR : ((dRadius > 0.0) ? dRadius : sOff.Length());
        adSide[f] = (sOff.Dot(sN) >= 0.0) ? 1.0 : -1.0;
        if (!(adR[f] > dTol))
          { SER_MSG(SM_ERR, _T("SmGnSectionFSG: the ball has no radius at the rails' start")); }
      }
  }

  // every station's feet and their series: values and s-derivatives (sigma order 1), jets to order n
  // (the ball's circle takes the rails and the centre only)
  // (and with a law the centre to one order more, for the touch circle's m')
  const sm_Trunc sTr = { (pLawRound != NULL) ? 2 : 1, (n >= 2) ? static_cast<int>(n) : 0 };
  std::vector<sm_GnStation> sStations(sS.size());
  for (size_t j=0; j<sS.size(); j++)
    {
      double adSeed[4];
      for (int f=0; f<2; f++)
        {
          SmPoint3d sW;
          if (pLawRound != NULL)
            { sW = SmPoint3d(sLawSeeds[j][static_cast<size_t>(2*f)], sLawSeeds[j][static_cast<size_t>(2*f+1)], 0.0); }
          else
            { SER(apUV[f]->EvaluatePoint(sS[j], sW)); }
          adSeed[2*f] = sW.x;
          adSeed[2*f+1] = sW.y;
        }
      // a constant ball: the station's centre on the centre line at the same offsets
      SmPoint3d sM;
      SER(rLeader.EvaluatePoint(sS[j], sM));
      const double dSolveTol = 1.0e-9 * smos_Max(1.0, SmVector3d(sM).Length());
      sm_Ser<double> sR;
      const sm_Ser<double> * pR = NULL;
      SER(sm_LawSeries(pGnLaw, sS[j], sTr.s, sR, pR));
      // with a ball of a given radius or law (pbNotTheFeet), a station that does not solve from the
      // solver's rails (Newton from within the tolerance of a regular root converges) says the rails
      // are not that ball's feet there
      if (pbNotTheFeet != NULL
          && sm_GnSolveStation(apFace, rLeader, sS[j], adSeed, adR, adSide, sTr, dSolveTol, sStations[j], pR, static_cast<int>(n)) != SM_SUCCESS)
        {
          *pbNotTheFeet = TRUE;
          return SM_SUCCESS;
        }
      if (pbNotTheFeet == NULL)
        { SER(sm_GnSolveStation(apFace, rLeader, sS[j], adSeed, adR, adSide, sTr, dSolveTol, sStations[j], pR, static_cast<int>(n))); }
      for (int f=0; f<2; f++)
        {
          const SmPoint3d sRho(sStations[j].asF[f].a[0].c[0][0], sStations[j].asF[f].a[1].c[0][0], sStations[j].asF[f].a[2].c[0][0]);
          // the ball's circle of a given radius (pbNotTheFeet): the solver's rails are this ball's feet
          // on the two faces, or the round is not their ball's (a rail held on a face's boundary, the
          // ball rolling on its edge), which the caller builds otherwise; with a law, the solver's
          // stations' feet are
          if (pbNotTheFeet != NULL)
            {
              if (pLawRound != NULL && !sLawFoot[j])
                { continue; }
              SmPoint3d sRail;
              SER(apRail[f]->EvaluatePoint((pLawRound != NULL) ? sLawCp[j] : sS[j], sRail));
              if (sRail.DistanceBetween(sRho) > dTol)
                {
                  *pbNotTheFeet = TRUE;
                  return SM_SUCCESS;
                }
              continue;
            }
          if (smos_Fabs(sM.DistanceBetween(sRho) - adR[f]) > smos_Max(dTol, 1.0e-6 * adR[f]))
            { SER_MSG(SM_ERR_INVALID_INPUT, _T("SmGnSectionFSG: a section of order 2 and up is built for a ball of one radius only")); }
        }
    }

  // the certificate: each span's bound (sm_GnCertifySpan) at most a tenth of the tolerance, a
  // span over it halved (a station solved at its middle) until none is
  {
    sm_Nurbs asN[2];
    const sm_Nurbs * apN[2] = { &asN[0], &asN[1] };
    for (int f=0; f<2; f++)
      {
        const SmBSplineSurface * pBS = SM_CAST_PTR(SmBSplineSurface, apFace[f]);
        if (pBS == NULL)
          { SER_MSG(SM_ERR, _T("SmGnSectionFSG: a face that is not a NURBS surface: its blend of order 2 and up cannot be certified")); }
        // a blend by its definition is enclosed only on its net, which is not the surface the feet are solved on
        if (apFace[f]->IsKindOf(SmProceduralBlend_TYPE))
          { SER_MSG(SM_ERR, _T("SmGnSectionFSG: a face that is a blend by its definition: a blend of order 2 and up on it cannot be certified until its definition is enclosed")); }
        SER(sm_NurbsOf(*pBS, asN[f]));
      }
    sm_NurbsCurve sLeader;
    SER(sm_NurbsCurveOf(rLeader, sLeader));
    const double dTarget = 0.1 * dTol;
    // a station wherever a foot crosses a knot line of its face: the face is only as smooth there
    // as the knot allows, and a span's enclosures are taken on one side of it
    for (int f=0; f<2; f++)
      for (int c=0; c<2; c++)
        {
          const std::vector<double> & crK = (c == 0) ? asN[f].U : asN[f].V;
          const int lDeg = (c == 0) ? asN[f].p : asN[f].q;
          std::vector<double> sInner;
          for (size_t m=static_cast<size_t>(lDeg + 1); m+static_cast<size_t>(lDeg + 1)<crK.size(); m++)
            if (sInner.empty() || crK[m] > sInner.back())
              { sInner.push_back(crK[m]); }
          std::vector<sm_GnStation> sNext;
          for (size_t j=0; j+1<sStations.size(); j++)
            {
              sNext.push_back(sStations[j]);
              const sm_GnStation & sA = sStations[j];
              const sm_GnStation & crB = sStations[j+1];
              for (size_t m=0; m<sInner.size(); m++)
                {
                  const double dK = sInner[m];
                  const double xa = sA.asX[2*f+c].c[0][0], xb = crB.asX[2*f+c].c[0][0];
                  const double dMargin = 1.0e-9 * (1.0 + std::fabs(dK));
                  if (!((xa < dK - dMargin && xb > dK + dMargin) || (xa > dK + dMargin && xb < dK - dMargin)))
                    { continue; }
                  // the parameter where that foot is on the knot line: the secant, kept in its bracket
                  double s0 = sA.dS, s1 = crB.dS, x0 = xa - dK, x1 = xb - dK;
                  sm_GnStation sAt;
                  for (int it=0; it<60; it++)
                    {
                      double dS = s0 - x0 * (s1 - s0) / (x1 - x0);
                      if (!(dS > s0 && dS < s1)) { dS = 0.5 * (s0 + s1); }
                      double adSeed[4];
                      const double t = (dS - sA.dS) / (crB.dS - sA.dS);
                      for (int q2=0; q2<4; q2++) { adSeed[q2] = (1.0 - t) * sA.asX[q2].c[0][0] + t * crB.asX[q2].c[0][0]; }
                      SmPoint3d sM;
                      SER(rLeader.EvaluatePoint(dS, sM));
                      sm_Ser<double> sR;
                      const sm_Ser<double> * pR = NULL;
                      SER(sm_LawSeries(pGnLaw, dS, sTr.s, sR, pR));
                      SER(sm_GnSolveStation(apFace, rLeader, dS, adSeed, adR, adSide, sTr, 1.0e-9 * smos_Max(1.0, SmVector3d(sM).Length()), sAt,
                                            pR, static_cast<int>(n)));
                      const double x = sAt.asX[2*f+c].c[0][0] - dK;
                      if (std::fabs(x) <= 1.0e-12 * (1.0 + std::fabs(dK)))
                        { break; }
                      if ((x < 0.0) == (x0 < 0.0)) { s0 = dS; x0 = x; } else { s1 = dS; x1 = x; }
                    }
                  // on the knot line the face's derivatives of order n + 1 differ on its two sides, and
                  // with them the jets' s-derivatives: the station is solved on the piece after it (the
                  // foot moves on to xb's side) and on the piece before, for the span ending there
                  {
                    const int lAfter = (xb > xa) ? 1 : -1;
                    int aiAfter[4] = { 0, 0, 0, 0 }, aiBefore[4] = { 0, 0, 0, 0 };
                    aiAfter[2*f+c] = lAfter;
                    aiBefore[2*f+c] = -lAfter;
                    double adOn[4];
                    for (int q2=0; q2<4; q2++) { adOn[q2] = sAt.asX[q2].c[0][0]; }
                    SmPoint3d sM;
                    SER(rLeader.EvaluatePoint(sAt.dS, sM));
                    const double dOnTol = 1.0e-9 * smos_Max(1.0, SmVector3d(sM).Length());
                    sm_Ser<double> sR;
                    const sm_Ser<double> * pR = NULL;
                    SER(sm_LawSeries(pGnLaw, sAt.dS, sTr.s, sR, pR));
                    SER(sm_GnSolveStation(apFace, rLeader, sAt.dS, adOn, adR, adSide, sTr, dOnTol, sAt, pR, static_cast<int>(n), &aiAfter));
                    for (int q2=0; q2<4; q2++) { sAt.aiBefore[q2] = aiBefore[q2]; }
                  }
                  sNext.push_back(sAt);
                }
            }
          sNext.push_back(sStations.back());
          std::sort(sNext.begin(), sNext.end(), [](const sm_GnStation & a, const sm_GnStation & b) { return a.dS < b.dS; });
          sStations.swap(sNext);
        }
    // the equal area thumbweight at every station, and I's enclosures on a grid over their range
    // (a little wider: the stations put in later lie between these)
    const double dRR = adR[0] * adR[1];
    sm_GnAreaCells sCells;
    sm_GnVarying sVary = { n, dRR, &sCells };
    if (bEqualArea)
      {
        double dWLo = HUGE_VAL, dWHi = -HUGE_VAL;
        for (size_t j=0; j<sStations.size(); j++)
          {
            SER(sm_GnStationThumbweight(n, dRR, sStations[j]));
            dWLo = smos_Min(dWLo, sStations[j].dW);
            dWHi = smos_Max(dWHi, sStations[j].dW);
          }
        const double dPad = 0.05 * (dWHi - dWLo) + 1.0e-6 * (1.0 + smos_Max(smos_Fabs(dWLo), smos_Fabs(dWHi)));
        sCells.Init(n, dWLo - dPad, smos_Min(dWHi + dPad, 0.5 * (1.0 + dWHi)), 64);
      }
    std::vector<double> sBounds(sStations.size() - 1, -1.0);
    ULONG lUnproved = 0, lOver = 0; // spans halved: Krawczyk's test failed on them, or their bound was over the target
    for (;;)
      {
        SmBoolean bAll = TRUE;
        for (size_t j=0; j+1<sStations.size(); j++)
          {
            if (sBounds[j] < 0.0)
              { SER(sm_GnCertifySpan(apN, sLeader, apFace, rLeader, sStations[j], sStations[j+1], adR, adSide, static_cast<int>(n), dTarget, sBounds[j],
                                     bEqualArea ? &sVary : NULL, pGnLaw)); }
            if (!(sBounds[j] <= dTarget))
              { bAll = FALSE; }
          }
        if (bAll)
          { break; }
        if (sStations.size() > 4096)
          { SER_MSG(SM_ERR, _T("SmGnSectionFSG: the blend could not be certified within its tolerance")); }
        std::vector<sm_GnStation> sNext;
        std::vector<double> sNextBounds;
        for (size_t j=0; j+1<sStations.size(); j++)
          {
            sNext.push_back(sStations[j]);
            if (sBounds[j] <= dTarget)
              {
                sNextBounds.push_back(sBounds[j]);
                continue;
              }
            if (sBounds[j] == HUGE_VAL) { lUnproved++; } else { lOver++; }
            const double dMid = 0.5 * (sStations[j].dS + sStations[j+1].dS);
            if (!(dMid > sStations[j].dS && dMid < sStations[j+1].dS))
              { SER_MSG(SM_ERR, _T("SmGnSectionFSG: the blend could not be certified within its tolerance")); }
            double adSeed[4];
            for (int c=0; c<4; c++) { adSeed[c] = 0.5 * (sStations[j].asX[c].c[0][0] + sStations[j+1].asX[c].c[0][0]); }
            SmPoint3d sM;
            SER(rLeader.EvaluatePoint(dMid, sM));
            sm_GnStation sMid;
            sm_Ser<double> sR;
            const sm_Ser<double> * pR = NULL;
            SER(sm_LawSeries(pGnLaw, dMid, sTr.s, sR, pR));
            SER(sm_GnSolveStation(apFace, rLeader, dMid, adSeed, adR, adSide, sTr, 1.0e-9 * smos_Max(1.0, SmVector3d(sM).Length()), sMid,
                                  pR, static_cast<int>(n)));
            if (bEqualArea)
              { SER(sm_GnStationThumbweight(n, dRR, sMid)); }
            sNext.push_back(sMid);
            sNextBounds.push_back(-1.0);
            sNextBounds.push_back(-1.0);
          }
        sNext.push_back(sStations.back());
        sStations.swap(sNext);
        sBounds.swap(sNextBounds);
      }
    m_dCertifiedBound = 0.0;
    for (size_t j=0; j<sBounds.size(); j++) { m_dCertifiedBound = smos_Max(m_dCertifiedBound, sBounds[j]); }
    m_lCertifiedSpans = static_cast<ULONG>(sBounds.size());
    double dWLo = m_dThumbweight, dWHi = m_dThumbweight;
    if (bEqualArea)
      {
        dWLo = HUGE_VAL;
        dWHi = -HUGE_VAL;
        for (size_t j=0; j<sStations.size(); j++)
          {
            dWLo = smos_Min(dWLo, sStations[j].dW);
            dWHi = smos_Max(dWHi, sStations[j].dW);
          }
      }
    TCHAR sMsg[SM_TBLOCK_SIZE];
    smos_sprintf(sMsg, _T("SmGnSectionFSG: order %lu, thumbweight %.6g to %.6g%s: %lu spans, certified within %.3g (target %.3g); halved %lu unproved, %lu over"),
                 static_cast<unsigned long>(n), dWLo, dWHi, bEqualArea ? _T(" (equal area)") : _T(""), static_cast<unsigned long>(m_lCertifiedSpans),
                 m_dCertifiedBound, dTarget, static_cast<unsigned long>(lUnproved), static_cast<unsigned long>(lOver));
    smos_ErrorMessage(SM_ERR_MESSAGE, _T(__FILE__), static_cast<ULONG>(__LINE__), sMsg);
  }
  sS.clear();
  for (size_t j=0; j<sStations.size(); j++) { sS.push_back(sStations[j].dS); }

  // each station's series on the pieces before it, for the span ending there: a Hermite span
  // interpolates the values and s-derivatives of its own pieces, which the certificate bounds
  for (size_t j=1; j<sStations.size(); j++)
    {
      sm_GnStation & rSt = sStations[j];
      double adOn[4];
      for (int q=0; q<4; q++) { adOn[q] = rSt.asX[q].c[0][0]; }
      SmPoint3d sM;
      SER(rLeader.EvaluatePoint(rSt.dS, sM));
      sm_GnStation sBefore;
      const int (&raiBefore)[4] = rSt.aiBefore;
      sm_Ser<double> sR;
      const sm_Ser<double> * pR = NULL;
      SER(sm_LawSeries(pGnLaw, rSt.dS, sTr.s, sR, pR, TRUE));
      SER(sm_GnSolveStation(apFace, rLeader, rSt.dS, adOn, adR, adSide, sTr, 1.0e-9 * smos_Max(1.0, SmVector3d(sM).Length()), sBefore,
                            pR, static_cast<int>(n), &raiBefore, NULL, NULL, TRUE));
      if (bEqualArea)
        {
          SER(sm_GnStationThumbweight(n, adR[0] * adR[1], sBefore));
          rSt.dWdsL = sBefore.dWds;
        }
      for (int q=0; q<4; q++) { rSt.asXL[q] = sBefore.asX[q]; }
      for (int f=0; f<2; f++) { rSt.asFL[f] = sBefore.asF[f]; }
      rSt.sVL = sBefore.sV;
      rSt.bSided = TRUE;
    }

  const size_t lD = (n >= 2) ? 3 * n + 2 : 3; // the section's control points across
  SmTArray<double> sUKnots, sVKnots;
  SmTArray<ULONG> sUMult, sVMult;
  sVKnots.Add(0.0); sVKnots.Add(1.0);
  sVMult.Add(static_cast<ULONG>(lD)); sVMult.Add(static_cast<ULONG>(lD));

  // the jets' Hermite curves along s, all on the stations' knots
  sm_HermiteKnots(sS, sUKnots, sUMult);
  std::vector< std::vector<SmVector3d> > asJets[2]; // [face][k] control points
  for (int f=0; f<2; f++)
    {
      asJets[f].resize(n + 1);
      for (ULONG k=0; k<=n; k++)
        {
          std::vector<SmVector3d> sV, sD, sDIn;
          for (size_t j=0; j<sS.size(); j++)
            {
              const sm_VSer<double> & crF = sStations[j].asF[f], & crFIn = sm_InF(sStations[j], f);
              sV.push_back(SmVector3d(crF.a[0].c[0][k], crF.a[1].c[0][k], crF.a[2].c[0][k]));
              sD.push_back(SmVector3d(crF.a[0].c[1][k], crF.a[1].c[1][k], crF.a[2].c[1][k]));
              sDIn.push_back(SmVector3d(crFIn.a[0].c[1][k], crFIn.a[1].c[1][k], crFIn.a[2].c[1][k]));
            }
          sm_HermitePoints(sS, sV, sD, sDIn, asJets[f][k]);
        }
    }
  const size_t lNU = asJets[0][0].size();

  // the homogeneous control points W P and weights W, [i along s][j across]
  std::vector<SmVector3d> sHomog(lNU * lD);
  std::vector<double> sW(lNU * lD);
  if (n == 1)
    {
      // the ball's circle: each station's arc rows (sm_ArcRows) and their s-derivatives, as the span
      // starting there and the one ending there take them, each a Hermite curve along s
      std::vector< std::vector<SmVector3d> > sHV(lD), sHD(lD), sHDIn(lD), sWV(lD), sWD(lD), sWDIn(lD);
      const double adO[3] = { 0.0, 0.0, 0.0 };
      for (size_t m=0; m<sStations.size(); m++)
        for (int lIn=0; lIn<2; lIn++)
          {
            const sm_GnStation & crSt = sStations[m];
            sm_Tay<double> asRho[2][3], asM[3];
            for (int q=0; q<3; q++)
              {
                for (int f=0; f<2; f++)
                  {
                    const sm_VSer<double> & crFD = lIn ? sm_InF(crSt, f) : crSt.asF[f];
                    asRho[f][q] = sm_TayOf(crSt.asF[f].a[q].c[0][0], crFD.a[q].c[1][0]);
                  }
                // the centre's s-derivatives, as the span the rows are for takes them
                const sm_VSer<double> & crV = (lIn && crSt.bSided) ? crSt.sVL : crSt.sV;
                asM[q] = sm_TayOf(crSt.sV.a[q].c[0][0], crV.a[q].c[1][0]);
                if (pLawRound != NULL)
                  { asM[q].c[2] = crV.a[q].c[2][0]; }
              }
            if (pLawRound != NULL)
              {
                // a varying ball: the circle it touches the envelope along, about m + k m',
                // k = (rho1 - m) . m' / m' . m'
                sm_Tay<double> asDm[3], sDD, sBD;
                for (int q=0; q<3; q++)
                  {
                    asDm[q] = sm_TayOf(asM[q].c[1], 2.0 * asM[q].c[2]);
                    asM[q].c[2] = 0.0;
                  }
                for (int q=0; q<3; q++)
                  {
                    sDD = sDD + asDm[q] * asDm[q];
                    sBD = sBD + (asRho[0][q] - asM[q]) * asDm[q];
                  }
                if (!(sDD.c[0] > 0.0))
                  { SER_MSG(SM_ERR, _T("SmGnSectionFSG: the balls' centres stand still")); }
                const sm_Tay<double> sK = sBD / sDD;
                for (int q=0; q<3; q++)
                  {
                    asM[q] = asM[q] + sK * asDm[q];
                    for (int i=2; i<=SM_GN_TAY; i++) { asM[q].c[i] = 0.0; }
                  }
              }
            sm_Tay<double> asH[3][3], asW[3];
            if (!sm_ArcRows(asRho[0], asRho[1], asM, adO, asH, asW))
              { SER_MSG(SM_ERR, _T("SmGnSectionFSG: a section of the ball's circle is a half turn or more")); }
            for (size_t j=0; j<lD; j++)
              {
                const SmVector3d sD(asH[j][0].c[1], asH[j][1].c[1], asH[j][2].c[1]);
                if (lIn)
                  {
                    sHDIn[j].push_back(sD);
                    sWDIn[j].push_back(SmVector3d(asW[j].c[1], 0.0, 0.0));
                    continue;
                  }
                sHV[j].push_back(SmVector3d(asH[j][0].c[0], asH[j][1].c[0], asH[j][2].c[0]));
                sHD[j].push_back(sD);
                sWV[j].push_back(SmVector3d(asW[j].c[0], 0.0, 0.0));
                sWD[j].push_back(SmVector3d(asW[j].c[1], 0.0, 0.0));
              }
          }
      for (size_t j=0; j<lD; j++)
        {
          std::vector<SmVector3d> sHP, sWP;
          sm_HermitePoints(sS, sHV[j], sHD[j], sHDIn[j], sHP);
          sm_HermitePoints(sS, sWV[j], sWD[j], sWDIn[j], sWP);
          for (size_t i=0; i<lNU; i++)
            {
              sHomog[i * lD + j] = sHP[i];
              sW[i * lD + j] = sWP[i].x;
            }
        }
    }
  else if (!bEqualArea)
    {
      // one section's form: its rows are linear in the jets, so the jets' curves carry them
      const sm_GnSectionForm sForm = sm_GnForm(n, m_dThumbweight);
      for (size_t i=0; i<lNU; i++)
        for (size_t j=0; j<lD; j++)
          {
            SmVector3d sP(0.0, 0.0, 0.0);
            for (int f=0; f<2; f++)
              for (ULONG k=0; k<=n; k++)
                sP = sP + asJets[f][k][i] * sForm.aG[f][k][j];
            sHomog[i * lD + j] = sP;
            sW[i * lD + j] = sForm.W[j];
          }
    }
  else
    {
      // each station's own form (over W_0): every row W_j P_j and weight W_j a Hermite curve of
      // its values and s-derivatives there (the form's in w(s) and the jets')
      std::vector< std::vector<SmVector3d> > sHV(lD), sHD(lD), sWV(lD), sWD(lD), sHDIn(lD), sWDIn(lD);
      for (size_t m=0; m<sStations.size(); m++)
        {
          const sm_GnStation & crSt = sStations[m];
          // the rows' s-derivatives as the span starting here takes them, then as the one ending here
          for (int lIn=0; lIn<2; lIn++)
            {
              const sm_GnSectionFormT< sm_Tay<double> > sF = sm_GnFormNormalized(n, sm_TayOf(crSt.dW, lIn ? sm_InWds(crSt) : crSt.dWds));
              for (size_t j=0; j<lD; j++)
                {
                  SmVector3d sV(0.0, 0.0, 0.0), sD(0.0, 0.0, 0.0);
                  for (int f=0; f<2; f++)
                    for (ULONG k=0; k<=n; k++)
                      {
                        const sm_VSer<double> & crF = crSt.asF[f], & crFD = lIn ? sm_InF(crSt, f) : crSt.asF[f];
                        const SmVector3d sC(crF.a[0].c[0][k], crF.a[1].c[0][k], crF.a[2].c[0][k]);
                        const SmVector3d sC1(crFD.a[0].c[1][k], crFD.a[1].c[1][k], crFD.a[2].c[1][k]);
                        const sm_Tay<double> & crG = sF.aG[f][k][j];
                        sV = sV + sC * crG.c[0];
                        sD = sD + sC * crG.c[1] + sC1 * crG.c[0];
                      }
                  if (lIn)
                    {
                      sHDIn[j].push_back(sD);
                      sWDIn[j].push_back(SmVector3d(sF.W[j].c[1], 0.0, 0.0));
                      continue;
                    }
                  sHV[j].push_back(sV);
                  sHD[j].push_back(sD);
                  sWV[j].push_back(SmVector3d(sF.W[j].c[0], 0.0, 0.0));
                  sWD[j].push_back(SmVector3d(sF.W[j].c[1], 0.0, 0.0));
                }
            }
        }
      for (size_t j=0; j<lD; j++)
        {
          std::vector<SmVector3d> sHP, sWP;
          sm_HermitePoints(sS, sHV[j], sHD[j], sHDIn[j], sHP);
          sm_HermitePoints(sS, sWV[j], sWD[j], sWDIn[j], sWP);
          for (size_t i=0; i<lNU; i++)
            {
              sHomog[i * lD + j] = sHP[i];
              sW[i * lD + j] = sWP[i].x;
            }
        }
    }
  double dWMin = SM_BIG_DOUBLE, dWMax = 0.0;
  for (size_t l=0; l<sW.size(); l++) { dWMin = smos_Min(dWMin, sW[l]); dWMax = smos_Max(dWMax, sW[l]); }
  if (!(dWMin > 0.0))
    { SER_MSG(SM_ERR, _T("SmGnSectionFSG: a weight of the section's rows is not positive")); }
  const double dWScale = 1.0 / smos_Sqrt(dWMin * dWMax);
  SmTArray<SmPoint3d> sCPs;
  SmTArray<double> sWts;
  for (size_t l=0; l<sW.size(); l++)
    {
      sCPs.Add(sHomog[l] * (1.0 / sW[l]));
      sWts.Add(sW[l] * dWScale);
    }
  SmBSplineSurface * pSurface = NULL;
  SER(SmBSplineSurface::CreateCanonical(crContext, 3, static_cast<ULONG>(lD - 1), sCPs, SM_SF_UNSPECIFIED, sUMult, sVMult, sUKnots,
                                        sVKnots, SM_KT_UNSPECIFIED, &sWts, NULL, pSurface));
  NER(pSurface);
  SmObjDelete sDelSurface(pSurface);

  // the feet's uv curves: the Hermite interpolants of the stations' feet, the definition's seeds
  SmBSplineCurve * apSeed[2] = { NULL, NULL };
  SmObjDelete asDelSeed[2];
  for (int f=0; f<2; f++)
    {
      std::vector<SmVector3d> sV, sD, sDIn, sUVCPs;
      for (size_t j=0; j<sS.size(); j++)
        {
          sV.push_back(SmVector3d(sStations[j].asX[2*f].c[0][0], sStations[j].asX[2*f+1].c[0][0], 0.0));
          sD.push_back(SmVector3d(sStations[j].asX[2*f].c[1][0], sStations[j].asX[2*f+1].c[1][0], 0.0));
          sDIn.push_back(SmVector3d(sm_InX(sStations[j], 2*f).c[1][0], sm_InX(sStations[j], 2*f+1).c[1][0], 0.0));
        }
      sm_HermitePoints(sS, sV, sD, sDIn, sUVCPs);
      SmTArray<SmPoint3d> sUVArr;
      for (size_t i=0; i<sUVCPs.size(); i++) { sUVArr.Add(sUVCPs[i]); }
      SmCurve * pCopy = NULL;
      SER(apUV[f]->Copy(crContext, pCopy));
      asDelSeed[f].SetObj(pCopy);
      apSeed[f] = SM_CAST_PTR(SmBSplineCurve, pCopy);
      NER(apSeed[f]);
      SER(apSeed[f]->SetCanonical(apUV[f]->GetDim(), 3, sUVArr, SM_CF_UNSPECIFIED, sUMult, sUKnots, SM_KT_UNSPECIFIED, NULL));
      if (bReplaceRails)
        {
          SmTArray<SmPoint3d> sRailCPs;
          for (size_t i=0; i<lNU; i++) { sRailCPs.Add(asJets[f][0][i]); }
          SER(apRail[f]->SetCanonical(3, 3, sRailCPs, SM_CF_UNSPECIFIED, sUMult, sUKnots, SM_KT_UNSPECIFIED, NULL));
          SER(apUV[f]->SetCanonical(apUV[f]->GetDim(), 3, sUVArr, SM_CF_UNSPECIFIED, sUMult, sUKnots, SM_KT_UNSPECIFIED, NULL));
        }
    }

  // with a law the rails now run over t: the centre line too, the balls' centres at the stations
  if (bReplaceRails && pLawRound != NULL)
    {
      std::vector<SmVector3d> sV, sD, sDIn, sCPs3;
      for (size_t j=0; j<sS.size(); j++)
        {
          const sm_VSer<double> & crV = sStations[j].sV, & crVIn = sStations[j].bSided ? sStations[j].sVL : sStations[j].sV;
          sV.push_back(SmVector3d(crV.a[0].c[0][0], crV.a[1].c[0][0], crV.a[2].c[0][0]));
          sD.push_back(SmVector3d(crV.a[0].c[1][0], crV.a[1].c[1][0], crV.a[2].c[1][0]));
          sDIn.push_back(SmVector3d(crVIn.a[0].c[1][0], crVIn.a[1].c[1][0], crVIn.a[2].c[1][0]));
        }
      sm_HermitePoints(sS, sV, sD, sDIn, sCPs3);
      SmTArray<SmPoint3d> sCentreCPs;
      for (size_t i=0; i<sCPs3.size(); i++) { sCentreCPs.Add(sCPs3[i]); }
      SER(rCenterLine.SetCanonical(3, 3, sCentreCPs, SM_CF_UNSPECIFIED, sUMult, sUKnots, SM_KT_UNSPECIFIED, NULL));
    }

  // the blend by its definition on this net: the faces, the ball's offsets and sides, the section,
  // the centre line as the leader the stations were solved in, the seeds, the certified bound
  const SmBSplineCurve * apConstSeed[2] = { apSeed[0], apSeed[1] };
  SmProceduralBlend * pBlend = NULL;
  SER(SmProceduralBlend::Create(crContext, *pSurface, apFace, adSide, adR, (pLawRound != NULL) ? pLawRound->pLaw : NULL, n, m_dThumbweight, rLeader, apConstSeed,
                                m_dCertifiedBound, pBlend));
  NER(pBlend);
  rpFilletSurface = pBlend;
  return SM_SUCCESS;
} // end SmGnSectionFSG::Build

/*******************************************************************//**
PURPOSE: The order of contact k of this section with face lRail (0 or 1)
    at the station of the centre line's parameter dS: the power at which
    the section stands off the face, c s^k, s along the face across the
    rail.

NOTES: S - C_i = (T_i - C_i)(1 - p1) + (1 - f)(T_o(p2) - T_i(1 - p1)) about
    the rail on face i (o the other face, C_i the curve on face i whose jet
    T_i is, sm_GnTopJets): 1 - f = t^(n+1) + ..., T_i - C_i = -C_i,n+1 (u -
    1)^(n+1) + ... (T_i's n-th coefficient differs from C_i's in the tangent
    plane only, which the standoff along N_i does not see) with u - 1 = -p1 =
    -t/(1 - w) + ..., so the coefficient of t^(n+1) in the standoff along the
    face's normal N_i is [(T_o(0) - rho_i) - C_i,n+1 (-1/(1 - w))^(n+1)] .
    N_i: k = n + 1 where it is not zero. On two planes the feet are linear in
    tau, C_i is the line of the feet, so C_i,n+1 = 0, and T_o(0) is the edge's
    point, which lies on face i: the coefficient vanishes and the next, (1 -
    w) (rho_o - e) . N_i (from p2's first term), does not: k = n + 2, the
    section G(n + 1).
***********************************************************************/
SmStatus SmGnSectionFSG::ContactOrderAt
  (const SmSurface      * apFace[2],   // in : the faces under the rails
   const SmBSplineCurve & crCenter,    // in : the centre line (the leader)
   double                 dS,          // in : its parameter at the station
   const double         (&crSeed)[4],  // in : the feet's (u1, v1, u2, v2), near
   const SmPoint3d      (&crRail)[2],  // in : the rails' points there
   ULONG                  lRail,       // in : the face
   ULONG                & rlOrder)     // out: k
  const
{
  rlOrder = m_lOrder + 1;
  if (m_lOrder < 2)
    { return SM_SUCCESS; } // the chord crosses the face; the circle touches it, and curves away
  NER(apFace[0]); NER(apFace[1]);
  if (lRail > 1 || !IsValid(m_lOrder, m_dThumbweight))
    { SER(SM_ERR_INVALID_INPUT); }
  const int n = static_cast<int>(m_lOrder);
  SmPoint3d sM;
  SER(crCenter.EvaluatePoint(dS, sM));
  double adR[2], adSide[2];
  for (int f=0; f<2; f++)
    {
      SmVector3d sN;
      SER(apFace[f]->EvaluateNormal(SmPoint2d(crSeed[2*f], crSeed[2*f+1]), TRUE, TRUE, sN));
      SER(sN.Unitize());
      const SmVector3d sOff = sM - crRail[f];
      adR[f] = sOff.Length();
      adSide[f] = (sOff.Dot(sN) >= 0.0) ? 1.0 : -1.0;
    }
  const sm_Trunc sTr = { 0, n + 1 };
  sm_GnStation sAt;
  sm_VSer<double> asC[2], asN[2];
  SER(sm_GnSolveStation(apFace, crCenter, dS, crSeed, adR, adSide, sTr, 1.0e-9 * smos_Max(1.0, SmVector3d(sM).Length()), sAt, NULL, n,
                        NULL, &asC, &asN));
  auto jet = [&](int f, int k) { return SmVector3d(sAt.asF[f].a[0].c[0][k], sAt.asF[f].a[1].c[0][k], sAt.asF[f].a[2].c[0][k]); };
  // the thumbweight of the station's own wedge: the ball's radii to its contacts from V
  const SmVector3d sV(sAt.sV.a[0].c[0][0], sAt.sV.a[1].c[0][0], sAt.sV.a[2].c[0][0]);
  double dW = 0.0;
  SER(ThumbweightAt((jet(0, 0) - sV).Dot(jet(1, 0) - sV) / (adR[0] * adR[1]), dW));
  const int i = static_cast<int>(lRail), o = 1 - i;
  SmVector3d sTo(0.0, 0.0, 0.0);
  for (int k=0; k<=n; k++) { sTo = sTo + jet(o, k) * ((k % 2 == 0) ? 1.0 : -1.0); }
  sm_Partials<double> sP;
  SER(sm_SurfacePartials(*apFace[i], SmPoint2d(sAt.asX[2*i].c[0][0], sAt.asX[2*i+1].c[0][0]), 1, sP));
  SmVector3d sNi = SmVector3d(sP.P[1][0].x, sP.P[1][0].y, sP.P[1][0].z) * SmVector3d(sP.P[0][1].x, sP.P[0][1].y, sP.P[0][1].z);
  SER(sNi.Unitize());
  const double dScale = smos_Pow(-1.0 / (1.0 - dW), static_cast<double>(n + 1));
  // the jet's n-th coefficient leaves out the face curve's tangent part tau_n: the section is off
  // the curve by -tau_n (u - 1)^n in the tangent plane at the rail, which the face's normal, turning
  // by N_1 (u - 1) along the curve, sees at order n + 1
  const SmVector3d sCn(asC[i].a[0].c[0][n], asC[i].a[1].c[0][n], asC[i].a[2].c[0][n]);
  const SmVector3d sTau = sCn - sNi * sCn.Dot(sNi);
  const SmVector3d sN1(asN[i].a[0].c[0][1], asN[i].a[1].c[0][1], asN[i].a[2].c[0][1]);
  const double dCoef = ((sTo - jet(i, 0)) - jet(i, n + 1) * dScale).Dot(sNi) - sTau.Dot(sN1) * dScale;
  const double dSize = jet(0, 0).DistanceBetween(jet(1, 0));
  rlOrder = (smos_Fabs(dCoef) > 1.0e-9 * smos_Max(1.0, dSize)) ? m_lOrder + 1 : m_lOrder + 2;
  return SM_SUCCESS;
} // end SmGnSectionFSG::ContactOrderAt

/*******************************************************************//**
PURPOSE: SolidCraft change. The blend's definition at (dS, dT) and its
    partials to lHighestS in s and lHighestT in t: aDerivs[i (lHighestT +
    1) + j] = d^(i+j) S / ds^i dt^j (SmProceduralBlend.h).

NOTES: One station solve: the feet by Newton from the seeds, then their
    series in s (and in tau for the G^n jets) order by order, as the G^n
    section's stations are built (sm_GnSolveStation). The section is then
    exact in t: n >= 2, sum_ik c_ik(s) phi_ik(t), phi_ik the section's
    rational functions (sm_GnForm) over its common weight; n = 1, the
    rational quadratic arc from the first rail to the second on the circle
    the balls' envelope touches the ball along: about m + ((rho - m) . m' /
    m' . m') m', in the plane square to the spine's tangent m' (the great
    circle where the radius is constant), weights 1, cos(theta / 2), 1. All
    of it as truncated series in (sigma, dt), so the partials are exact to
    rounding.
***********************************************************************/
SmStatus SmProceduralBlend::EvaluateDefinition(double dS, double dT, ULONG lHighestS, ULONG lHighestT, SmVector3d * aDerivs) const
{
  NER(aDerivs); NER(m_apSupport[0]); NER(m_apSupport[1]); NER(m_pLeader); NER(m_apSeed[0]); NER(m_apSeed[1]);
  CountEvaluation();
  const int lS = static_cast<int>(lHighestS), lT = static_cast<int>(lHighestT);
  const int n = static_cast<int>(m_lOrder);
  // the station: the jets to order n in tau (G^n), or the rails and the centre one order more in s (the circle's m')
  sm_Trunc sTr;
  sTr.s = (n >= 2) ? lS : lS + 1;
  sTr.e = (n >= 2) ? n : 0;
  if (lHighestS > 64 || lHighestT > 64 || sTr.s > SM_GN_MAXS || sTr.e > SM_GN_MAXE || sTr.s + sTr.e + 1 > SM_GN_MAXK)
    { SER_MSG(SM_ERR_INVALID_INPUT, _T("SmProceduralBlend: partials of this order are past what its definition's series hold")); }
  double adSeed[4];
  for (int f=0; f<2; f++)
    {
      SmPoint3d sW;
      SER(m_apSeed[f]->EvaluatePoint(dS, sW));
      adSeed[2*f] = sW.x;
      adSeed[2*f+1] = sW.y;
    }
  double adR[2] = { m_adR[0], m_adR[1] };
  sm_Ser<double> sR = sm_SerConst<double>(0.0);
  const sm_Ser<double> * pR = NULL;
  if (m_pLaw != NULL)
    {
      SmVector3d asD[SM_GN_MAXS + 1];
      SER(m_pLaw->Evaluate(dS, static_cast<ULONG>(sTr.s), TRUE, asD, FALSE));
      double dFact = 1.0;
      for (int i=0; i<=sTr.s; i++)
        {
          if (i > 0) { dFact *= i; }
          sR.c[i][0] = asD[i].x / dFact;
        }
      if (!(asD[0].x > 0.0))
        { SER_MSG(SM_ERR, _T("SmProceduralBlend: its radius law is not positive here")); }
      adR[0] = adR[1] = asD[0].x;
      pR = &sR;
    }
  sm_GnStation sAt;
  SmBoolean bMemo = FALSE;
  sm_StationMemo * pMemo = NULL;
  for (int m=0; m<SM_STATION_MEMO && !bMemo; m++)
    {
      sm_StationMemo & rMemo = t_asMemo[m];
      if (rMemo.lSerial == m_lSerial && rMemo.dS == dS && rMemo.ts == sTr.s && rMemo.te == sTr.e)
        {
          sAt = rMemo.sStation;
          pMemo = &rMemo;
          bMemo = TRUE;
        }
    }
  if (!bMemo)
    {
      SmPoint3d sL;
      SER(m_pLeader->EvaluatePoint(dS, sL));
      const SmSurface * apFace[2] = { m_apSupport[0], m_apSupport[1] };
      // (a nested support's own evaluations use the memo meanwhile: the station is kept only once solved)
      const SmStatus eSolved = sm_GnSolveStation(apFace, *m_pLeader, dS, adSeed, adR, m_adSide, sTr,
                                                 1.0e-9 * smos_Max(1.0, SmVector3d(sL).Length()), sAt, pR, n >= 2 ? n : 0);
      CountSolve(static_cast<ULONG>(sAt.lSteps));
      SER(eSolved);
      sm_StationMemo & rMemo = t_asMemo[t_lMemoNext];
      t_lMemoNext = (t_lMemoNext + 1) % SM_STATION_MEMO;
      rMemo.lSerial = m_lSerial;
      rMemo.dS = dS;
      rMemo.ts = sTr.s;
      rMemo.te = sTr.e;
      rMemo.sStation = sAt;
      rMemo.bW = FALSE;
      pMemo = &rMemo;
    }

  sm_TV sSum = sm_TVZero(lS, lT);
  if (n >= 2 && m_dThumbweight == SM_GN_THUMBWEIGHT_EQUAL_AREA)
    {
      // the equal area thumbweight of this station's own wedge, as a series in s (the builder's
      // rule, sm_GnStationThumbweight, to every order asked): w(s) solves I(w) = target(u(s)), u
      // half the angle between the ball's radii to its contacts
      if (lS > SM_GN_TAY)
        { SER_MSG(SM_ERR_INVALID_INPUT, _T("SmProceduralBlend: partials of this order are past what its thumbweight's series holds")); }
      // (kept with the station: the root is the costly part, and a row of points shares it)
      sm_Tay<double> sWs;
      if (pMemo != NULL && pMemo->bW)
        { sWs = pMemo->sW; }
      else
        {
          const sm_Tay<double> sU = sm_GnHalfAngle(sm_GnWedgeCos<double>(sAt.asF, sAt.sV, m_adR[0] * m_adR[1], lS));
          double dW0 = 0.0;
          SER(sm_GnEqualAreaW(m_lOrder, sU.c[0], dW0));
          const sm_Tay<double> sD = sm_GnAreaSum(m_lOrder, sm_TayOf(dW0, 1.0));
          if (!(sD.c[1] < 0.0))
            { SER_MSG(SM_ERR, _T("SmProceduralBlend: the section's area does not fall with its thumbweight")); }
          sWs = sm_GnInvert(sm_GnTarget(sU), sD, dW0);
          if (pMemo != NULL)
            {
              pMemo->sW = sWs;
              pMemo->bW = TRUE;
            }
        }
      const sm_GnSectionFormT< sm_Tay<double> > sForm = sm_GnFormT< sm_Tay<double> >(m_lOrder, sWs);
      // a row of the form in (sigma, dt): sum_j g_j(sigma) B_j(t0 + dt)
      const size_t lDeg = sForm.W.size();
      std::vector< std::vector<double> > sBasis(lDeg);
      for (size_t j=0; j<lDeg; j++)
        {
          sm_Bern sUnit(lDeg, 0.0);
          sUnit[j] = 1.0;
          sBasis[j] = sm_BernTaylor(sUnit, dT, lT);
        }
      auto row = [&](const std::vector< sm_Tay<double> > & crG)
      {
        sm_T2 sR(lS, lT);
        for (size_t j=0; j<lDeg; j++)
          for (int i=0; i<=lS; i++)
            for (int q=0; q<=lT; q++)
              sR.at(i, q) += crG[j].c[i] * sBasis[j][static_cast<size_t>(q)];
        return sR;
      };
      sm_T2 sInvW;
      SER(sm_T2Pow(row(sForm.W), -1.0, sInvW));
      for (int f=0; f<2; f++)
        for (int k=0; k<=n; k++)
          sSum = sSum + sm_TVOf(sAt.asF[f], k, lS, lT) * (row(sForm.aG[f][static_cast<size_t>(k)]) * sInvW);
    }
  else if (n >= 2)
    {
      const sm_GnSectionForm sForm = sm_GnForm(m_lOrder, m_dThumbweight);
      sm_T2 sW(0, lT);
      const std::vector<double> sWt = sm_BernTaylor(sForm.W, dT, lT);
      for (int j=0; j<=lT; j++) { sW.at(0, j) = sWt[static_cast<size_t>(j)]; }
      sm_T2 sInvW;
      SER(sm_T2Pow(sW, -1.0, sInvW));
      for (int f=0; f<2; f++)
        for (int k=0; k<=n; k++)
          {
            sm_T2 sN(0, lT);
            const std::vector<double> sNt = sm_BernTaylor(sForm.aG[f][static_cast<size_t>(k)], dT, lT);
            for (int j=0; j<=lT; j++) { sN.at(0, j) = sNt[static_cast<size_t>(j)]; }
            const sm_T2 sPhi1 = sN * sInvW;
            sm_T2 sPhi(lS, lT);
            for (int j=0; j<=lT; j++) { sPhi.at(0, j) = sPhi1.at(0, j); }
            sSum = sSum + sm_TVOf(sAt.asF[f], k, lS, lT) * sPhi;
          }
    }
  else
    {
      // the centre to one order more, for its tangent m'
      const sm_TV sM1 = sm_TVOf(sAt.sV, 0, lS + 1, lT);
      sm_TV sD = sm_TVZero(lS, lT), sM = sm_TVZero(lS, lT);
      for (int c=0; c<3; c++)
        for (int i=0; i<=lS; i++)
          {
            sD.a[c].at(i, 0) = (i + 1) * sM1.a[c].at(i + 1, 0);
            sM.a[c].at(i, 0) = sM1.a[c].at(i, 0);
          }
      const sm_TV sRho1 = sm_TVOf(sAt.asF[0], 0, lS, lT), sRho2 = sm_TVOf(sAt.asF[1], 0, lS, lT);
      sm_T2 sInvDD;
      SER(sm_T2Pow(sm_TVDot(sD, sD), -1.0, sInvDD));
      // the touch circle's centre: the ball's moved along m' by (rho - m) . m' / m' . m' (zero where r is constant)
      const sm_TV sMc = sM + sD * (sm_TVDot(sRho1 - sM, sD) * sInvDD);
      const sm_TV sB1 = sRho1 - sMc, sB2 = sRho2 - sMc;
      sm_T2 sI1, sI2;
      SER(sm_T2Pow(sm_TVDot(sB1, sB1), -0.5, sI1));
      SER(sm_T2Pow(sm_TVDot(sB2, sB2), -0.5, sI2));
      const sm_T2 sHalf = (sm_TVDot(sB1, sB2) * sI1 * sI2 + sm_T2(lS, lT, 1.0)) * 0.5; // cos^2(theta / 2)
      sm_T2 sW, sInvW2;
      SER(sm_T2Pow(sHalf, 0.5, sW));
      SER(sm_T2Pow(sHalf, -1.0, sInvW2));
      // the middle control point, where the tangents at the rails meet
      const sm_TV sP1 = sMc + ((sRho1 + sRho2) * 0.5 - sMc) * sInvW2;
      sm_T2 asB[3] = { sm_T2(lS, lT), sm_T2(lS, lT), sm_T2(lS, lT) };
      const double t0 = dT;
      const double adB[3][3] = { { (1.0 - t0) * (1.0 - t0), -2.0 * (1.0 - t0), 1.0 },
                                 { 2.0 * t0 * (1.0 - t0), 2.0 - 4.0 * t0, -2.0 },
                                 { t0 * t0, 2.0 * t0, 1.0 } };
      for (int b=0; b<3; b++)
        for (int j=0; j<=lT && j<3; j++)
          asB[b].at(0, j) = adB[b][j];
      const sm_T2 sWB1 = sW * asB[1];
      sm_T2 sInvDen;
      SER(sm_T2Pow(asB[0] + sWB1 + asB[2], -1.0, sInvDen));
      sSum = (sRho1 * asB[0] + sP1 * sWB1 + sRho2 * asB[2]) * sInvDen;
    }
  double dFi = 1.0;
  for (int i=0; i<=lS; i++)
    {
      if (i > 0) { dFi *= i; }
      double dFj = 1.0;
      for (int j=0; j<=lT; j++)
        {
          if (j > 0) { dFj *= j; }
          aDerivs[i * (lT + 1) + j] = SmVector3d(sSum.a[0].at(i, j), sSum.a[1].at(i, j), sSum.a[2].at(i, j)) * (dFi * dFj);
        }
    }
  return SM_SUCCESS;
} // end SmProceduralBlend::EvaluateDefinition
