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
* T_i(u) = sum_{k<=n} c_ik (u - 1)^k the order-n Taylor polynomial of F_i at
* tau = 1 (c_ik = F_i^(k)(1) / k!). The section is
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
* exactly.
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
#include <SmFilletIntersector.h>
#include <SmGraphicsExtern.h>
#include <SmGeomUtility.h>
#include <SmBrep.h>

#include <algorithm>
#include <array>
#include <cmath>
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
// one, evaluated in its net's parameters) by its net, continued past its knots; any other surface
// by its own evaluator.
SmStatus sm_SurfacePartials(const SmSurface & crSurface, const SmPoint2d & crUV, int lK, sm_Partials<double> & rP)
{
  if (lK < 0 || lK > SM_GN_MAXK)
    { SER(SM_ERR); }
  const SmBSplineSurface * pBS = crSurface.IsKindOf(SmBSplineSurface_TYPE) ? static_cast<const SmBSplineSurface *>(&crSurface) : NULL;
  if (pBS != NULL)
    {
      sm_Nurbs sN;
      SER(sm_NurbsOf(*pBS, sN));
      const int ku = sm_Span(sN.U, sN.p, sN.nu, crUV.x), kv = sm_Span(sN.V, sN.q, sN.nv, crUV.y);
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
};

// G and the feet F_i as series, given the feet's parameters as series x[0..3] (u1, v1, u2, v2),
// whose constant terms are the point the partials were taken at.
template <class T> SmStatus sm_GnEvaluate(const sm_GnSystem<T> & crSys, const sm_Ser<T> (&x)[4], const sm_Trunc & tr,
                                          sm_Ser<T> (&rG)[4], sm_VSer<T> (&rF)[2])
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
      // the offset towards the ball: tau r_i s_i, tau = 1 + eps
      sm_Ser<T> sD = sm_SerConst<T>(T(crSys.adR[f] * crSys.adSide[f]));
      if (!crSys.bHold && tr.e >= 1)
        { sD.c[0][1] = T(crSys.adR[f] * crSys.adSide[f]); }
      asV[f] = sm_VSerAdd(sS, sm_VSerTimes(sN, sm_SerMul(sInvLen, sD, tr), tr), tr);
      rF[f] = sS;
    }
  const sm_VSer<T> sDiff = sm_VSerSub(asV[0], asV[1], tr);
  for (int j=0; j<3; j++) { rG[j] = sDiff.a[j]; }
  rG[3] = sm_VSerDot(sm_VSerSub(asV[0], crSys.sLeader, tr), crSys.sTangent, tr);
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
};

/*******************************************************************//**
   The leader (centre line) at s as series in sigma up to tr.s: its point,
   and its unit tangent (which takes one derivative more).
***********************************************************************/
SmStatus sm_LeaderSeries(const SmBSplineCurve & crLeader, double dS, const sm_Trunc & tr, sm_VSer<double> & rL, sm_VSer<double> & rT)
{
  const sm_Trunc sSigma = { tr.s, 0 };
  SmVector3d asD[SM_GN_MAXS + 2];
  SER(crLeader.Evaluate(dS, static_cast<ULONG>(tr.s + 1), TRUE, asD, FALSE));
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
                           sm_GnStation & rStation)
{
  rStation.dS = dS;
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
  SER(sm_LeaderSeries(crLeader, dS, tr, sSys.sLeader, sSys.sTangent));

  // the feet at tau = 1, by Newton (the leader at sigma = 0 only)
  double x[4] = { crSeed[0], crSeed[1], crSeed[2], crSeed[3] };
  const sm_Trunc sPoint = { 0, 0 };
  SmBoolean bConverged = FALSE;
  for (int it=0; it<50 && !bConverged; it++)
    {
      for (int f=0; f<2; f++)
        { SER(sm_SurfacePartials(*apFace[f], SmPoint2d(x[2*f], x[2*f+1]), 2, asP[f])); }
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

  // the series, order by order: each coefficient of G is J x[i][k] plus the lower ones'
  for (int f=0; f<2; f++)
    { SER(sm_SurfacePartials(*apFace[f], SmPoint2d(x[2*f], x[2*f+1]), lK, asP[f])); }
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
        if (k > tr.e)
          { continue; }
        double b[4], y[4];
        for (int j=0; j<4; j++) { b[j] = -sG[j].c[i][k]; }
        SER(sm_Solve4(J, b, y));
        for (int q=0; q<4; q++) { rStation.asX[q].c[i][k] = y[q]; }
      }
    }
  sm_Ser<double> sG[4];
  SER(sm_GnEvaluate(sSys, rStation.asX, tr, sG, rStation.asF));
  return SM_SUCCESS;
}

/*******************************************************************//**
   Bernstein polynomials: coefficients of degree size() - 1.
***********************************************************************/
typedef std::vector<double> sm_Bern;

sm_Bern sm_BernMul(const sm_Bern & a, const sm_Bern & b)
{
  const size_t p = a.size() - 1, q = b.size() - 1;
  sm_Bern r(p + q + 1, 0.0);
  for (size_t i=0; i<=p; i++)
    for (size_t j=0; j<=q; j++)
      {
        double c = 1.0;
        // C(p, i) C(q, j) / C(p + q, i + j)
        double cp = 1.0, cq = 1.0, cpq = 1.0;
        for (size_t m=1; m<=i; m++) { cp *= static_cast<double>(p - i + m) / m; }
        for (size_t m=1; m<=j; m++) { cq *= static_cast<double>(q - j + m) / m; }
        for (size_t m=1; m<=i+j; m++) { cpq *= static_cast<double>(p + q - i - j + m) / m; }
        c = cp * cq / cpq;
        r[i+j] += c * a[i] * b[j];
      }
  return r;
}

sm_Bern sm_BernPow(const sm_Bern & a, ULONG n)
{
  sm_Bern r(1, 1.0);
  for (ULONG m=0; m<n; m++) { r = sm_BernMul(r, a); }
  return r;
}

/*******************************************************************//**
   The section of order n and thumbweight w over the common denominator:
   its weights W_j and, for each jet c_ik, the coefficient g_ik_j of c_ik in
   W_j P_j (j = 0..3n+1).
***********************************************************************/
struct sm_GnSectionForm
{
  sm_Bern                W;
  std::vector<sm_Bern>   aG[2]; // aG[i][k]
};

sm_GnSectionForm sm_GnForm(ULONG n, double w)
{
  sm_GnSectionForm sForm;
  sm_Bern sA(n + 2, 0.0), sB(n + 2, 0.0);
  sA[0] = 1.0;           // (1 - t)^(n+1)
  sB[n + 1] = 1.0;       // t^(n+1)
  const sm_Bern sOneMinusWt = { 1.0, 1.0 - w };  // 1 - w t
  const sm_Bern sOneMinusWpWt = { 1.0 - w, 1.0 }; // 1 - w + w t
  const sm_Bern sMinusT = { 0.0, -1.0 };          // -t
  const sm_Bern sTMinus1 = { -1.0, 0.0 };         // t - 1
  sm_Bern sAB(n + 2, 0.0);
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
   station), to order lOrderS in s and n in tau, enclosed: rF, with rbProved
   FALSE where Krawczyk's test fails on S.
***********************************************************************/
SmStatus sm_GnEnclose(const sm_Nurbs * apFace[2], const sm_NurbsCurve & crLeader, const sm_GnStation & crA, const sm_GnStation & crB,
                      const double (&crR)[2], const double (&crSide)[2], int n, int lOrderS, sm_VSer<sm_Ivl> (&rF)[2], SmBoolean & rbProved,
                      sm_Ser<sm_Ivl> (*pX)[4] = NULL, double (*pDS)[2] = NULL)
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
        if (k > sTr.e)
          { continue; }
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
  SM_GN_TRY(sm_GnEvaluate(sSys, sX, sTr, sG, rF));
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
                         double (&rM4)[2][SM_GN_MAXE + 1], double (&rUV4)[2])
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
  SM_GN_TRY(sm_GnEnclose(apFace, crLeader, crA, crB, crR, crSide, n, 4 + SM_GN_TAYLOR, asF5, bProved, &asX5, &adDS));
  if (!bProved)
    { return SM_SUCCESS; }
  // the middle's station, and the series there
  const double dMid = 0.5 * (crA.dS + crB.dS);
  double adSeed[4];
  for (int c=0; c<4; c++) { adSeed[c] = 0.5 * (crA.asX[c].c[0][0] + crB.asX[c].c[0][0]); }
  SmPoint3d sM;
  SM_GN_TRY(crCenter.EvaluatePoint(dMid, sM));
  const sm_Trunc sPointTr = { 0, 0 };
  sm_GnStation sMid;
  SM_GN_TRY(sm_GnSolveStation(apSurf, crCenter, dMid, adSeed, crR, crSide, sPointTr, 1.0e-9 * smos_Max(1.0, SmVector3d(sM).Length()), sMid));
  SM_GN_TRY(sm_GnEnclose(apFace, crLeader, sMid, sMid, crR, crSide, n, 4 + SM_GN_TAYLOR - 1, asF4, bProved, &asX4, &adDSm));
  if (!bProved)
    { return SM_SUCCESS; }
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
   against its face).
***********************************************************************/
SmStatus sm_GnCertifySpan(const sm_Nurbs * apFace[2], const sm_NurbsCurve & crLeader, const SmSurface * apSurf[2], const SmBSplineCurve & crCenter,
                          const sm_GnStation & crA, const sm_GnStation & crB, const double (&crR)[2], const double (&crSide)[2], int n,
                          double dTarget, double & rdBound)
{
  rdBound = HUGE_VAL;
  const double h = crB.dS - crA.dS;
  double dLast = HUGE_VAL;
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
          if (sm_GnSolveStation(apSurf, crCenter, dS, adSeed, crR, crSide, sPointTr, 1.0e-9 * smos_Max(1.0, SmVector3d(sM).Length()),
                                sEnds[static_cast<size_t>(i)]) != SM_SUCCESS)
            { bSolved = FALSE; }
        }
      if (!bSolved)
        { continue; }
      double adM4[2][SM_GN_MAXE + 1], adUV4[2] = { 0.0, 0.0 };
      for (int f=0; f<2; f++)
        for (int k=0; k<=n; k++)
          adM4[f][k] = 0.0;
      SmBoolean bProved = TRUE;
      for (int i=0; i<lPieces && bProved; i++)
        {
          double adP[2][SM_GN_MAXE + 1], adPUV[2];
          SM_GN_TRY(sm_GnFourthOver(apFace, crLeader, apSurf, crCenter, sEnds[static_cast<size_t>(i)], sEnds[static_cast<size_t>(i + 1)],
                                    crR, crSide, n, adP, adPUV));
          for (int f=0; f<2; f++)
            {
              for (int k=0; k<=n; k++)
                {
                  if (!(adP[f][k] < HUGE_VAL)) { bProved = FALSE; }
                  adM4[f][k] = std::max(adM4[f][k], adP[f][k]);
                }
              if (!(adPUV[f] < HUGE_VAL)) { bProved = FALSE; }
              adUV4[f] = std::max(adUV4[f], adPUV[f]);
            }
        }
      if (!bProved)
        { continue; }
      double dBound = 0.0;
      for (int f=0; f<2; f++)
        {
          double dSum = 0.0;
          for (int k=0; k<=n; k++) { dSum += h * h * h * h / 384.0 * adM4[f][k]; }
          dBound = std::max(dBound, std::max(dSum, h * h * h * h / 384.0 * adUV4[f]));
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

} // namespace

/*******************************************************************//**
PURPOSE: A blend's section of order lOrder (0 to SM_GN_MAX_ORDER) and
    thumbweight dThumbweight (-1 < w < 1, or SM_GN_THUMBWEIGHT_EQUAL_AREA, used from order 2). Order 1's
    circles are built as SmArcSection's defaults say.
***********************************************************************/
SmGnSectionFSG::SmGnSectionFSG(ULONG lOrder, double dThumbweight)
  : m_lOrder(lOrder), m_dThumbweight(dThumbweight), m_dCertifiedBound(-1.0), m_lCertifiedSpans(0)
{
}

namespace
{
// The integral of a b' - a' b over [0, 1] for the order-n section between two planes, where
// S(t) = a(t) e1 + b(t) e2 with e1, e2 the rails' points from the edge: a = f (1 - p1),
// b = (1 - f) p2 (SmFilletSection.cpp's header). Twice the area the section takes off, over
// |e1 x e2|; it does not depend on the planes' angle. Composite Simpson on 2048 intervals: the
// integrand is a rational function smooth on [0, 1].
double sm_GnSectionAreaIntegral(ULONG n, double w)
{
  const double dN1 = static_cast<double>(n + 1);
  auto ab = [&](double t, double & a, double & b, double & da, double & db)
  {
    const double dA = smos_Pow(1.0 - t, dN1), dB = smos_Pow(t, dN1);
    const double f = dA / (dA + dB);
    const double dfa = -dN1 * smos_Pow(1.0 - t, dN1 - 1.0), dfb = dN1 * smos_Pow(t, dN1 - 1.0);
    const double df = (dfa * dB - dA * dfb) / ((dA + dB) * (dA + dB));
    const double d1 = 1.0 - w + w * t, d2 = 1.0 - w * t;
    const double p1 = t / d1, dp1 = (1.0 - w) / (d1 * d1);
    const double p2 = t * (1.0 - w) / d2, dp2 = (1.0 - w) / (d2 * d2);
    a = f * (1.0 - p1);           da = df * (1.0 - p1) - f * dp1;
    b = (1.0 - f) * p2;           db = -df * p2 + (1.0 - f) * dp2;
  };
  const int N = 2048;
  double dSum = 0.0;
  for (int i=0; i<=N; i++)
    {
      const double t = static_cast<double>(i) / N;
      double a, b, da, db;
      ab(t, a, b, da, db);
      const double g = a * db - da * b;
      dSum += (i == 0 || i == N ? 1.0 : (i % 2 ? 4.0 : 2.0)) * g;
    }
  return dSum / (3.0 * N);
}
}

/*******************************************************************//**
PURPOSE: The equal area thumbweight: the w at which the order-n section
    between two planes at the angle alpha takes off what the ball's circle
    does on the same rails.
NOTES: With e1, e2 the rails' points from the edge (|e_i| = d = r cot(alpha/2)),
    the section takes off (1/2) d^2 sin(alpha) I(w), I the integral of
    a b' - a' b (sm_GnSectionAreaIntegral), and the circle r^2 (cot(alpha/2)
    - (pi - alpha)/2). So I(w) = (cot(alpha/2) - (pi - alpha)/2) /
    (cot(alpha/2) cos^2(alpha/2)), a target that depends on alpha only. I
    falls as w rises, so the root is bisected on (-1, 1). A target past
    what the range reaches (a wedge so open or so shut the circle takes off
    next to nothing, or nearly the chord's area) is an error naming it.
***********************************************************************/
SmStatus SmGnSectionFSG::EqualAreaThumbweight(ULONG lOrder, double dAlpha, double & rdW)
{
  rdW = 0.0;
  if (lOrder < 2 || lOrder > SM_GN_MAX_ORDER || !(dAlpha > 0.0 && dAlpha < SM_PI))
    { SER(SM_ERR_INVALID_INPUT); }
  const double dHalf = 0.5 * dAlpha;
  const double dCot = smos_Cosine(dHalf) / smos_Sine(dHalf);
  const double dTarget = (dCot - 0.5 * (SM_PI - dAlpha)) / (dCot * smos_Cosine(dHalf) * smos_Cosine(dHalf));
  double dLo = -1.0 + 1.0e-9, dHi = 1.0 - 1.0e-9;
  const double dILo = sm_GnSectionAreaIntegral(lOrder, dLo), dIHi = sm_GnSectionAreaIntegral(lOrder, dHi);
  if (!(dTarget <= dILo && dTarget >= dIHi))
    { SER_MSG(SM_ERR_INVALID_INPUT, _T("SmGnSectionFSG: no thumbweight gives this wedge's section the ball's circle's area")); }
  for (int it=0; it<200 && dHi - dLo > 1.0e-15; it++)
    {
      const double dMid = 0.5 * (dLo + dHi);
      if (sm_GnSectionAreaIntegral(lOrder, dMid) > dTarget) { dLo = dMid; } else { dHi = dMid; }
    }
  rdW = 0.5 * (dLo + dHi);
  return SM_SUCCESS;
}

/*******************************************************************//**
PURPOSE: The thumbweight used, the equal area one resolved once from the
    wedge at one station: the ball's centre crCentre and its contacts
    crRail1, crRail2 give the wedge's angle alpha = pi - (the angle between
    the radii). Exact for two planes; elsewhere the default of the station
    it is resolved at.
***********************************************************************/
SmStatus SmGnSectionFSG::ResolveThumbweight(const SmPoint3d & crCentre, const SmPoint3d & crRail1, const SmPoint3d & crRail2) const
{
  if (m_dThumbweight != SM_GN_THUMBWEIGHT_EQUAL_AREA)
    { return SM_SUCCESS; }
  SmVector3d sA = crRail1 - crCentre, sB = crRail2 - crCentre;
  const double dLa = sA.Length(), dLb = sB.Length();
  if (!(dLa > 0.0 && dLb > 0.0))
    { SER(SM_ERR_INVALID_INPUT); }
  double dCos = sA.Dot(sB) / (dLa * dLb);
  dCos = smos_Max(-1.0, smos_Min(1.0, dCos));
  const double dAlpha = SM_PI - smos_ArcCosine(dCos);
  double dW = 0.0;
  SER(EqualAreaThumbweight(m_lOrder, dAlpha, dW));
  m_dThumbweight = dW;
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
PURPOSE: Create a fillet surface with a circular cross section.  This
    is a lofted surface or in some special cases a cylinder or a torus.

NOTES:
***********************************************************************/
SmStatus SmGnSectionFSG::CircleSurface
  (SmFilletSolver & rFilletSolver,
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
        SER(CircleSurfaceFromCurves(rFilletSolver.GetCreationContext(),
            *pCenterLine,*pRail1Curve,*pRail2Curve,dThisApproxTol3d,
            bNeedG1Continuity,pSurf1,pUV1,pSurf2,pUV2,rpFilletSurface));
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
            dThisApproxTol3d, bNeedG1Continuity, pSurf1, pUV1, pSurf2, pUV2, pPiece));
        NER(pPiece);
        SmObjDelete sDeletePiece(pPiece);

        if (pJoined == NULL)
        {
            sDeletePiece.Clear();
            pJoined = pPiece;
            sDeleteJoined.SetObj(pJoined);
            continue;
        }

        // join the piece on at its u-min side; no knot removal there, the surface turns
        SmBSplineSurface * pNext = NULL;
        SER(pJoined->JoinSurface(crContext, pPiece, SM_SP_U, 0.0, pNext));
        NER(pNext);
        delete pJoined;
        pJoined = pNext;
        sDeleteJoined.SetObj(pJoined);
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

  // see if newSurf == RuledSurface
  double    dAnalyticTol = 1.0e-8;
  SmPoint3d  sLinePoint;
  SmVector3d sLineVector;
  SmBoolean bIsDeg1Line =    rCenterLine.GetDegree() == 1
                          && rCenterLine.IsLine(11,dAnalyticTol,sLinePoint,sLineVector);
  if (bIsDeg1Line)
    {
      bSweep = TRUE;
      SmPoint3d  sLinePoint2;
      SmVector3d sLineVector2;
      SmBoolean bIsLineR1 = rRail1Curve.IsLine(11,dAnalyticTol,sLinePoint2,sLineVector2);
      double dAngle;
      if ( bIsLineR1 )
        {
          SER(sLineVector2.AngleBetween(sLineVector,dAngle));
          if (dAngle > SM_EFF_ZERO_SQRT)
            {
              bSweep = FALSE;
            }
        }
      SmBoolean bIsLineR2 = rRail2Curve.IsLine(11,dAnalyticTol,sLinePoint2,sLineVector2);
      if ( bIsLineR2 )
        {
          SER(sLineVector2.AngleBetween(sLineVector,dAngle));
          if (dAngle > SM_EFF_ZERO_SQRT)
            {
              bSweep = FALSE;
            }
        }
    } // end NewSurf == SurfaceOfRevolution check

  // see if NewSurf == SurfaceOfRevolution
  double dArcRadius, dStartAng, dEndAng;
  SmAxis2Placement sRefFrame;
  SmBoolean bIsArc = rCenterLine.IsArc(11, dAnalyticTol,
                                       sRefFrame, dArcRadius,
                                       dStartAng, dEndAng);
  if (bIsArc)
    {
      //cbi this is never hit in prog_test.
      bRevolve = TRUE;
      double dRadius2, dStartAng2, dEndAng2;
      double dRadius3, dStartAng3, dEndAng3;
      SmAxis2Placement sRefFrame2, sRefFrame3;
      SmBoolean bIsArc2 = rRail1Curve.IsArc(11, dAnalyticTol, sRefFrame2,
                                             dRadius2, dStartAng2, dEndAng2);
      if ( !bIsArc2 )
        { bRevolve = FALSE; }
      else if ( smos_Fabs( dEndAng2 - dEndAng ) * dRadius2 > dThisApproxTol3d )
        { bRevolve = FALSE; }
      else if ( ! sRefFrame2.AreCoaxial( sRefFrame, dThisApproxTol3d, dRadius2 ) )
        { bRevolve = FALSE; }

      if ( bRevolve )
      {
          SmBoolean bIsArc3 = rRail2Curve.IsArc(11, dAnalyticTol, sRefFrame3,
                                                 dRadius3, dStartAng3, dEndAng3);
          if ( !bIsArc3 )
            { bRevolve = FALSE; }
          else if ( smos_Fabs( dEndAng3 - dEndAng ) * dRadius3 > dThisApproxTol3d )
            { bRevolve = FALSE; }
          else if ( ! sRefFrame3.AreCoaxial( sRefFrame, dThisApproxTol3d, dRadius3 ) )
            { bRevolve = FALSE; }
          // Also check against 1s Rail curve.
          else if ( smos_Fabs( dEndAng3 - dEndAng2 ) * dRadius3 > dThisApproxTol3d )
            { bRevolve = FALSE; }
          else if ( ! sRefFrame3.AreCoaxial( sRefFrame2, dThisApproxTol3d, dRadius3 ) )
            { bRevolve = FALSE; }
      }
    }

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

              // get center, rail1, rail2 points for dT
              SER(rCenterLine.Evaluate(dT,1,TRUE,sCentPV));
              SER(rRail1Curve.EvaluatePoint(dT,sP1));
              SER(rRail2Curve.EvaluatePoint(dT,sP2));

              // get sample Point gaps
              SmVector3d sV1 = sP1 - sCentPV[0];
              SmVector3d sV2 = sP2 - sCentPV[0];
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
   const SmTArray<SmTsectPnt*> & /*crFilletPoints*/,
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
      SER(CircleSurface(rFilletSolver, pCenterLine, pRail1Curve, pRail2Curve, pSurf1, pUV1, pSurf2, pUV2, rpFilletSurface));
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
      SER(CircleSurfaceFromCurves(crContext, rCenterLine, rRail1Curve, rRail2Curve, dThisApproxTol3d, bNeedG1Continuity,
                                  pSurf1, pUV1, pSurf2, pUV2, rpFilletSurface));
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
   SmBSplineSurface *& rpFilletSurface)
{
  rpFilletSurface = NULL;
  NER(pFace1); NER(pFace2); NER(pUV1); NER(pUV2);
  const ULONG n = m_lOrder;
  if (!IsValid(n, m_dThumbweight))
    { SER_MSG(SM_ERR_INVALID_INPUT, _T("SmGnSectionFSG: the order must be 2 to 4 and the thumbweight in (-1, 1)")); }
  const SmSurface * apFace[2] = { pFace1, pFace2 };
  SmBSplineCurve * apRail[2] = { &rRail1Curve, &rRail2Curve };
  SmBSplineCurve * apUV[2] = { pUV1, pUV2 };
  const double dScale = smos_Max(1.0, rCenterLine.GetNaturalInterval().GetLength());

  // the stations: the unique knots of the rails and the centre line
  SmTArray<double> sKnots;
  std::vector<double> sS;
  const SmBSplineCurve * apCurves[3] = { &rCenterLine, &rRail1Curve, &rRail2Curve };
  const SmExtent1d sIvl = rCenterLine.GetNaturalInterval();
  {
    // the default thumbweight resolved from the wedge at the middle station
    const double dMidS = 0.5 * (sIvl.GetMin() + sIvl.GetMax());
    SmPoint3d sC, sR1, sR2;
    SER(rCenterLine.EvaluatePoint(dMidS, sC));
    SER(rRail1Curve.EvaluatePoint(dMidS, sR1));
    SER(rRail2Curve.EvaluatePoint(dMidS, sR2));
    SER(ResolveThumbweight(sC, sR1, sR2));
  }
  for (int c=0; c<3; c++)
    {
      const SmExtent1d sCIvl = apCurves[c]->GetNaturalInterval();
      if (   smos_Fabs(sCIvl.GetMin() - sIvl.GetMin()) > 1.0e-9 * dScale
          || smos_Fabs(sCIvl.GetMax() - sIvl.GetMax()) > 1.0e-9 * dScale)
        { SER_MSG(SM_ERR, _T("SmGnSectionFSG: the rails and the centre line run over different intervals")); }
      SER(apCurves[c]->GetKnots(sKnots));
      for (ULONG j=0; j<sKnots.GetSize(); j++) { sS.push_back(sKnots[j]); }
    }
  std::sort(sS.begin(), sS.end());
  {
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
    SmPoint3d sM;
    SER(rCenterLine.EvaluatePoint(sS.front(), sM));
    for (int f=0; f<2; f++)
      {
        SmPoint3d sRho, sW;
        SER(apRail[f]->EvaluatePoint(sS.front(), sRho));
        SER(apUV[f]->EvaluatePoint(sS.front(), sW));
        const SmPoint2d sUV(sW.x, sW.y);
        SmVector3d sN;
        SER(apFace[f]->EvaluateNormal(sUV, TRUE, TRUE, sN));
        SER(sN.Unitize());
        const SmVector3d sOff = sM - sRho;
        adR[f] = sOff.Length();
        adSide[f] = (sOff.Dot(sN) >= 0.0) ? 1.0 : -1.0;
        if (!(adR[f] > dTol))
          { SER_MSG(SM_ERR, _T("SmGnSectionFSG: the ball has no radius at the rails' start")); }
      }
  }

  // every station's feet and their series: values and s-derivatives (sigma order 1), jets to order n
  const sm_Trunc sTr = { 1, static_cast<int>(n) };
  std::vector<sm_GnStation> sStations(sS.size());
  for (size_t j=0; j<sS.size(); j++)
    {
      double adSeed[4];
      for (int f=0; f<2; f++)
        {
          SmPoint3d sW;
          SER(apUV[f]->EvaluatePoint(sS[j], sW));
          adSeed[2*f] = sW.x;
          adSeed[2*f+1] = sW.y;
        }
      // a constant ball: the station's centre on the centre line at the same offsets
      SmPoint3d sM;
      SER(rCenterLine.EvaluatePoint(sS[j], sM));
      const double dSolveTol = 1.0e-9 * smos_Max(1.0, SmVector3d(sM).Length());
      SER(sm_GnSolveStation(apFace, rCenterLine, sS[j], adSeed, adR, adSide, sTr, dSolveTol, sStations[j]));
      for (int f=0; f<2; f++)
        {
          const SmPoint3d sRho(sStations[j].asF[f].a[0].c[0][0], sStations[j].asF[f].a[1].c[0][0], sStations[j].asF[f].a[2].c[0][0]);
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
        SER(sm_NurbsOf(*pBS, asN[f]));
      }
    sm_NurbsCurve sLeader;
    SER(sm_NurbsCurveOf(rCenterLine, sLeader));
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
                      SER(rCenterLine.EvaluatePoint(dS, sM));
                      SER(sm_GnSolveStation(apFace, rCenterLine, dS, adSeed, adR, adSide, sTr, 1.0e-9 * smos_Max(1.0, SmVector3d(sM).Length()), sAt));
                      const double x = sAt.asX[2*f+c].c[0][0] - dK;
                      if (std::fabs(x) <= 1.0e-12 * (1.0 + std::fabs(dK)))
                        { break; }
                      if ((x < 0.0) == (x0 < 0.0)) { s0 = dS; x0 = x; } else { s1 = dS; x1 = x; }
                    }
                  sNext.push_back(sAt);
                }
            }
          sNext.push_back(sStations.back());
          std::sort(sNext.begin(), sNext.end(), [](const sm_GnStation & a, const sm_GnStation & b) { return a.dS < b.dS; });
          sStations.swap(sNext);
        }
    std::vector<double> sBounds(sStations.size() - 1, -1.0);
    ULONG lUnproved = 0, lOver = 0; // spans halved: Krawczyk's test failed on them, or their bound was over the target
    for (;;)
      {
        SmBoolean bAll = TRUE;
        for (size_t j=0; j+1<sStations.size(); j++)
          {
            if (sBounds[j] < 0.0)
              { SER(sm_GnCertifySpan(apN, sLeader, apFace, rCenterLine, sStations[j], sStations[j+1], adR, adSide, static_cast<int>(n), dTarget, sBounds[j])); }
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
            SER(rCenterLine.EvaluatePoint(dMid, sM));
            sm_GnStation sMid;
            SER(sm_GnSolveStation(apFace, rCenterLine, dMid, adSeed, adR, adSide, sTr, 1.0e-9 * smos_Max(1.0, SmVector3d(sM).Length()), sMid));
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
    TCHAR sMsg[SM_TBLOCK_SIZE];
    smos_sprintf(sMsg, _T("SmGnSectionFSG: order %lu, thumbweight %.6g: %lu spans, certified within %.3g (target %.3g); halved %lu unproved, %lu over"),
                 static_cast<unsigned long>(n), m_dThumbweight, static_cast<unsigned long>(m_lCertifiedSpans), m_dCertifiedBound, dTarget,
                 static_cast<unsigned long>(lUnproved), static_cast<unsigned long>(lOver));
    smos_ErrorMessage(SM_ERR_MESSAGE, _T(__FILE__), static_cast<ULONG>(__LINE__), sMsg);
  }
  sS.clear();
  for (size_t j=0; j<sStations.size(); j++) { sS.push_back(sStations[j].dS); }

  // the section's form
  const sm_GnSectionForm sForm = sm_GnForm(n, m_dThumbweight);
  const size_t lD = sForm.W.size(); // 3n + 2 control points across
  double dWMin = SM_BIG_DOUBLE, dWMax = 0.0;
  for (size_t j=0; j<lD; j++) { dWMin = smos_Min(dWMin, sForm.W[j]); dWMax = smos_Max(dWMax, sForm.W[j]); }
  if (!(dWMin > 0.0))
    { SER(SM_ERR); }
  const double dWScale = 1.0 / smos_Sqrt(dWMin * dWMax);

  SmTArray<double> sUKnots, sVKnots;
  SmTArray<ULONG> sUMult, sVMult;
  sVKnots.Add(0.0); sVKnots.Add(1.0);
  sVMult.Add(static_cast<ULONG>(lD)); sVMult.Add(static_cast<ULONG>(lD));
  SmTArray<SmPoint3d> sCPs;
  SmTArray<double> sWts;

  // the jets' Hermite curves along s, all on the stations' knots
  sm_HermiteKnots(sS, sUKnots, sUMult);
  std::vector< std::vector<SmVector3d> > asJets[2]; // [face][k] control points
  for (int f=0; f<2; f++)
    {
      asJets[f].resize(n + 1);
      for (ULONG k=0; k<=n; k++)
        {
          std::vector<SmVector3d> sV, sD;
          for (size_t j=0; j<sS.size(); j++)
            {
              const sm_VSer<double> & crF = sStations[j].asF[f];
              sV.push_back(SmVector3d(crF.a[0].c[0][k], crF.a[1].c[0][k], crF.a[2].c[0][k]));
              sD.push_back(SmVector3d(crF.a[0].c[1][k], crF.a[1].c[1][k], crF.a[2].c[1][k]));
            }
          sm_HermitePoints(sS, sV, sD, asJets[f][k]);
        }
    }
  const size_t lNU = asJets[0][0].size();
  for (size_t i=0; i<lNU; i++)
    for (size_t j=0; j<lD; j++)
      {
        SmVector3d sP(0.0, 0.0, 0.0);
        for (int f=0; f<2; f++)
          for (ULONG k=0; k<=n; k++)
            sP = sP + asJets[f][k][i] * sForm.aG[f][k][j];
        sCPs.Add(sP * (1.0 / sForm.W[j]));
        sWts.Add(sForm.W[j] * dWScale);
      }
  SmBSplineSurface * pSurface = NULL;
  SER(SmBSplineSurface::CreateCanonical(crContext, 3, static_cast<ULONG>(lD - 1), sCPs, SM_SF_UNSPECIFIED, sUMult, sVMult, sUKnots,
                                        sVKnots, SM_KT_UNSPECIFIED, &sWts, NULL, pSurface));
  NER(pSurface);
  SmObjDelete sDelSurface(pSurface);

  if (bReplaceRails)
    {
      for (int f=0; f<2; f++)
        {
          SmTArray<SmPoint3d> sRailCPs;
          for (size_t i=0; i<lNU; i++) { sRailCPs.Add(asJets[f][0][i]); }
          SER(apRail[f]->SetCanonical(3, 3, sRailCPs, SM_CF_UNSPECIFIED, sUMult, sUKnots, SM_KT_UNSPECIFIED, NULL));
          {
              std::vector<SmVector3d> sV, sD, sUVCPs;
              for (size_t j=0; j<sS.size(); j++)
                {
                  sV.push_back(SmVector3d(sStations[j].asX[2*f].c[0][0], sStations[j].asX[2*f+1].c[0][0], 0.0));
                  sD.push_back(SmVector3d(sStations[j].asX[2*f].c[1][0], sStations[j].asX[2*f+1].c[1][0], 0.0));
                }
              sm_HermitePoints(sS, sV, sD, sUVCPs);
              SmTArray<SmPoint3d> sUVArr;
              for (size_t i=0; i<sUVCPs.size(); i++) { sUVArr.Add(sUVCPs[i]); }
              SER(apUV[f]->SetCanonical(apUV[f]->GetDim(), 3, sUVArr, SM_CF_UNSPECIFIED, sUMult, sUKnots, SM_KT_UNSPECIFIED, NULL));
            }
        }
    }

  sDelSurface.Clear();
  rpFilletSurface = pSurface;
  return SM_SUCCESS;
} // end SmGnSectionFSG::Build

/*******************************************************************//**
PURPOSE: The order of contact k of this section with face lRail (0 or 1)
    at the station of the centre line's parameter dS: the power at which
    the section stands off the face, c s^k, s along the face across the
    rail.

NOTES: S - F_i = (T_i - F_i)(1 - p1) + (1 - f)(T_o(p2) - T_i(1 - p1)) about
    the rail on face i (o the other face, F_i the foot curve): 1 - f = t^(n+1)
    + ..., T_i - F_i = -c_i,n+1 (u - 1)^(n+1) + ... with u - 1 = -p1 = -t/(1 - w)
    + ..., so the coefficient of t^(n+1) in the standoff along the face's
    normal N_i is [(T_o(0) - rho_i) - c_i,n+1 (-1/(1 - w))^(n+1)] . N_i:
    k = n + 1 where it is not zero. On two planes the feet are linear in tau,
    so c_i,n+1 = 0 and T_o(0) is the edge's point, which lies on face i: the
    coefficient vanishes and the next, (1 - w) (rho_o - e) . N_i (from p2's
    first term), does not: k = n + 2, the section G(n + 1).
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
  SER(ResolveThumbweight(sM, crRail[0], crRail[1]));
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
  SER(sm_GnSolveStation(apFace, crCenter, dS, crSeed, adR, adSide, sTr, 1.0e-9 * smos_Max(1.0, SmVector3d(sM).Length()), sAt));
  auto jet = [&](int f, int k) { return SmVector3d(sAt.asF[f].a[0].c[0][k], sAt.asF[f].a[1].c[0][k], sAt.asF[f].a[2].c[0][k]); };
  const int i = static_cast<int>(lRail), o = 1 - i;
  SmVector3d sTo(0.0, 0.0, 0.0);
  for (int k=0; k<=n; k++) { sTo = sTo + jet(o, k) * ((k % 2 == 0) ? 1.0 : -1.0); }
  sm_Partials<double> sP;
  SER(sm_SurfacePartials(*apFace[i], SmPoint2d(sAt.asX[2*i].c[0][0], sAt.asX[2*i+1].c[0][0]), 1, sP));
  SmVector3d sNi = SmVector3d(sP.P[1][0].x, sP.P[1][0].y, sP.P[1][0].z) * SmVector3d(sP.P[0][1].x, sP.P[0][1].y, sP.P[0][1].z);
  SER(sNi.Unitize());
  const double dScale = smos_Pow(-1.0 / (1.0 - m_dThumbweight), static_cast<double>(n + 1));
  const double dCoef = ((sTo - jet(i, 0)) - jet(i, n + 1) * dScale).Dot(sNi);
  const double dSize = jet(0, 0).DistanceBetween(jet(1, 0));
  rlOrder = (smos_Fabs(dCoef) > 1.0e-9 * smos_Max(1.0, dSize)) ? m_lOrder + 1 : m_lOrder + 2;
  return SM_SUCCESS;
} // end SmGnSectionFSG::ContactOrderAt
