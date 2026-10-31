// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

// The operations on bodies the C API cannot make yet (a cylinder through SM_API), and the
// SMLib defects found under ASan, each as a regression case. Run it under ASan/UBSan
// (README, "Fuzzing"): every case must stay clean there, not only return the right status.

#include "../src/smlib_shim.h"

#include <SmApiBrep.h>
#include <SmApiFillets.h>
#include <SmApiPrimitives.h>
#include <SmApiQueries.h>
#include <SmBSplineCurve.h>
#include <SmBSplineSurface.h>
#include <SmCone.h>
#include <SmCylinder.h>
#include <SmExtent3d.h>
#include <SmEdge.h>
#include <SmEdgeuse.h>
#include <SmFace.h>
#include <SmFaceuse.h>
#include <SmFilletExecutive.h>
#include <SmLoop.h>
#include <SmOffsetGeometryCreation.h>
#include <SmOffsetSurface.h>
#include <SmPlane.h>
#include <SmRegion.h>
#include <SmAxis2Placement.h>
#include <SmSolutionArray.h>
#include <SmSphere.h>
#include <SmSurface.h>
#include <SmSurfOfRevolution.h>
#include <SmTrimmingTools.h>
#include <SmLine.h>
#include <SmVector3d.h>
#include <SmVertex.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <iterator>
#include <string>
#include <utility>
#include <vector>

namespace
{

int failures = 0;

#define CHECK(cond)                                                                                \
    do                                                                                             \
    {                                                                                              \
        if (!(cond))                                                                               \
        {                                                                                          \
            std::fprintf(stderr, "%s:%d: CHECK(%s) failed\n", __FILE__, __LINE__, #cond);          \
            failures++;                                                                            \
        }                                                                                          \
    } while (0)

const double PI = 3.14159265358979323846;

smlib_brep *box(double x, double y, double z, double dx, double dy, double dz)
{
    const double o[3] = {x, y, z};
    const double s[3] = {dx, dy, dz};
    smlib_brep *b = nullptr;
    CHECK(smlib_make_box(o, s, &b) == SMLIB_OK);
    return b;
}

// A Z-axis cylinder through SM_API (the C API has no primitives besides the box yet).
smlib_brep *cylinder(double x, double y, double z, double r, double h)
{
    SmApiGetOrCreateContext();
    SmBrep *p = nullptr;
    CHECK(SmApiCreateCylinder(SmVector3d(x, y, z), r, h, p) == SM_SUCCESS && p != nullptr);
    return new smlib_brep{p};
}

double volume_of(const smlib_brep *b)
{
    double v = -1.0;
    CHECK(smlib_brep_volume(b, 1e-7, &v) == SMLIB_OK);
    return v;
}

// SMLib's own validity check, every test (level 2), must find nothing wrong with `b`.
bool sound(const smlib_brep *b)
{
    static char buf[1 << 16];
    int64_t len = 0;
    int64_t failures = -1;
    if (smlib_brep_check(b, 2, buf, sizeof buf, &len, &failures) != SMLIB_OK)
    {
        std::fprintf(stderr, "  smlib_brep_check did not run\n");
        return false;
    }
    if (failures != 0)
    {
        std::fprintf(stderr, "  %lld failed validity tests:\n%s\n", static_cast<long long>(failures), buf);
    }
    return failures == 0;
}

// Checks a result and frees it.
void done(smlib_brep *&r)
{
    CHECK(r != nullptr && sound(r));
    smlib_brep_free(r);
    r = nullptr;
}

bool near(double got, double want, double rel)
{
    bool ok = std::fabs(got - want) <= rel * std::fabs(want);
    if (!ok)
    {
        std::fprintf(stderr, "  got %.9g, want %.9g\n", got, want);
    }
    return ok;
}

// A sphere of radius `r` about the origin in six patches, each on its own sphere surface turned
// to it, as an importer hands over a cube projected onto the sphere: the topology of a box, its
// faces sphere patches (reference direction to the patch's middle, longitude and latitude in
// [-45, 45] degrees), its edges the great-circle arcs between the corners.
smlib_brep *cube_sphere(double r)
{
    using V = std::array<double, 3>;
    auto unit = [](V a) {
        double l = std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
        return V{a[0] / l, a[1] / l, a[2] / l};
    };
    smlib_brep *b = box(-1, -1, -1, 2, 2, 2);
    smlib_flat *f = nullptr;
    CHECK(smlib_brep_export(b, &f) == SMLIB_OK);
    smlib_brep_free(b);
    const int64_t *ev = nullptr;
    const int64_t *loops = nullptr;
    const int64_t *uses = nullptr;
    const int64_t *kinds = nullptr;
    const double *vp = nullptr;
    int64_t n_ev = 0, n_loops = 0, n_uses = 0, n_faces = 0, n_vp = 0;
    CHECK(smlib_flat_get_i64(f, SMLIB_FLAT_EDGE_VERTICES, &ev, &n_ev) == SMLIB_OK);
    CHECK(smlib_flat_get_i64(f, SMLIB_FLAT_LOOP_EDGEUSE_COUNT, &loops, &n_loops) == SMLIB_OK);
    CHECK(smlib_flat_get_i64(f, SMLIB_FLAT_EDGEUSE_EDGE, &uses, &n_uses) == SMLIB_OK);
    CHECK(smlib_flat_get_i64(f, SMLIB_FLAT_FACE_SURFACE_TYPE, &kinds, &n_faces) == SMLIB_OK);
    CHECK(smlib_flat_get_f64(f, SMLIB_FLAT_VERTEX_POSITION, &vp, &n_vp) == SMLIB_OK);
    CHECK(n_faces == 6 && n_loops == 6 && n_ev == 24 && n_vp == 24);
    std::vector<double> position, centre, axis, ref, radius, range;
    for (int64_t i = 0; i < n_vp / 3; ++i)
    {
        V u = unit({vp[3 * i], vp[3 * i + 1], vp[3 * i + 2]});
        for (double x : u)
        {
            position.push_back(r * x);
        }
    }
    for (int64_t e = 0; e < n_ev / 2; ++e)
    {
        V a = unit({vp[3 * ev[2 * e]], vp[3 * ev[2 * e] + 1], vp[3 * ev[2 * e] + 2]});
        V z = unit({vp[3 * ev[2 * e + 1]], vp[3 * ev[2 * e + 1] + 1], vp[3 * ev[2 * e + 1] + 2]});
        V n = unit({a[1] * z[2] - a[2] * z[1], a[2] * z[0] - a[0] * z[2], a[0] * z[1] - a[1] * z[0]});
        for (int k = 0; k < 3; ++k)
        {
            centre.push_back(0.0);
            axis.push_back(n[k]);
            ref.push_back(a[k]);
        }
        radius.push_back(r);
        range.push_back(0.0);
        range.push_back(std::acos(a[0] * z[0] + a[1] * z[1] + a[2] * z[2]));
    }
    std::vector<double> s_centre, s_axis, s_ref, s_radius, s_range;
    int64_t at = 0;
    for (int64_t i = 0; i < n_faces && at < n_uses; ++i)
    {
        // The patch's direction: the middle of its loop's corners (one loop per face, in order).
        V mid{0, 0, 0};
        for (int64_t k = 0; k < loops[i]; ++k)
        {
            int64_t e = uses[at + k];
            for (int j = 0; j < 2; ++j)
            {
                for (int q = 0; q < 3; ++q)
                {
                    mid[q] += vp[3 * ev[2 * e + j] + q];
                }
            }
        }
        at += loops[i];
        V d = unit(mid);
        V up = std::fabs(d[2]) > 0.5 ? V{1, 0, 0} : V{0, 0, 1};
        for (int q = 0; q < 3; ++q)
        {
            s_centre.push_back(0.0);
            s_axis.push_back(up[q]);
            s_ref.push_back(d[q]);
        }
        s_radius.push_back(r);
        for (double x : {-PI / 4.0, -PI / 4.0, PI / 4.0, PI / 4.0})
        {
            s_range.push_back(x);
        }
    }
    std::vector<int64_t> circles(static_cast<size_t>(n_ev / 2), SMLIB_FLAT_CURVE_CIRCLE);
    std::vector<int64_t> spheres(static_cast<size_t>(n_faces), SMLIB_FLAT_SURFACE_SPHERE);
    CHECK(smlib_flat_set_i64(f, SMLIB_FLAT_EDGE_CURVE_TYPE, circles.data(), static_cast<int64_t>(circles.size())) == SMLIB_OK);
    CHECK(smlib_flat_set_i64(f, SMLIB_FLAT_FACE_SURFACE_TYPE, spheres.data(), static_cast<int64_t>(spheres.size())) == SMLIB_OK);
    for (int id : {SMLIB_FLAT_EDGE_LINE_ORIGIN, SMLIB_FLAT_EDGE_LINE_DIRECTION, SMLIB_FLAT_FACE_PLANE_ORIGIN, SMLIB_FLAT_FACE_PLANE_AXIS,
                   SMLIB_FLAT_FACE_PLANE_REF_DIRECTION})
    {
        CHECK(smlib_flat_set_f64(f, id, nullptr, 0) == SMLIB_OK);
    }
    const std::pair<int, const std::vector<double> *> arrays[] = {
        {SMLIB_FLAT_VERTEX_POSITION, &position},       {SMLIB_FLAT_EDGE_CIRCLE_CENTER, &centre},
        {SMLIB_FLAT_EDGE_CIRCLE_AXIS, &axis},          {SMLIB_FLAT_EDGE_CIRCLE_REF_DIRECTION, &ref},
        {SMLIB_FLAT_EDGE_CIRCLE_RADIUS, &radius},      {SMLIB_FLAT_EDGE_RANGE, &range},
        {SMLIB_FLAT_FACE_SPHERE_CENTER, &s_centre},    {SMLIB_FLAT_FACE_SPHERE_AXIS, &s_axis},
        {SMLIB_FLAT_FACE_SPHERE_REF_DIRECTION, &s_ref}, {SMLIB_FLAT_FACE_SPHERE_RADIUS, &s_radius},
        {SMLIB_FLAT_FACE_RANGE, &s_range}};
    for (const auto &a : arrays)
    {
        CHECK(smlib_flat_set_f64(f, a.first, a.second->data(), static_cast<int64_t>(a.second->size())) == SMLIB_OK);
    }
    smlib_brep *out = nullptr;
    CHECK(smlib_brep_import(f, &out) == SMLIB_OK);
    smlib_flat_free(f);
    return out;
}

// A 120 x 80 x 40 box drafted 5 degrees inward with its four sloping edges rounded r 4, shelled 2
// open at the top (SolidCraft rig shell/3576). The offset's analytic copy of each round (a cone
// spanning a quarter turn) was grown round its axis to a full turn, which SmCone's NURBS rebuild
// refuses (SM_ERR_INVALID_INPUT); MateEdges failed at every round, the rounds' walls were
// left unmated and the piecewise merge hit a missing vertex (SM_ERR_NULL_POINTER, 1006). The
// copy is grown only along the rulings, where the offset's tangent extension is still the cone.
// The corner of corner_setback left open (smlib_blend_setback_open): SMLib's patch, six
// quadrilaterals meeting at its centre, goes with its six inner edges and its centre, so the
// sheet keeps the box's six faces and the three rounds (9), 36 - 6 edges and 23 - 1 vertices.
// Its one hole runs through the six section ends, each side (a section or a face curve) in the
// two halves SMLib split it into: twelve edges, each between the hole and one face, which are
// the sheet's only free edges, chained end to start in the direction the faces beside them run
// them. Too small a buffer makes nothing and names the length needed.
void corner_setback_open()
{
    smlib_brep *part = box(0, 0, 0, 60, 40, 20);
    const double edges[9] = {30, 40, 20, 60, 20, 20, 60, 40, 10};
    const double radii[3] = {5, 5, 5};
    const double corner[3] = {60, 40, 20};
    const double setback = 8.0;
    smlib_brep *r = nullptr;
    int64_t len = -1;
    CHECK(smlib_blend_setback_open(part, SMLIB_BLEND_ROUND, edges, radii, 3, corner, &setback, 1, 0.1, &r, nullptr, 0, &len, nullptr,
                                   nullptr) == SMLIB_ERR_BUFFER &&
          r == nullptr && len == 1 + 3 * 12);
    int64_t holes[64];
    CHECK(smlib_blend_setback_open(part, SMLIB_BLEND_ROUND, edges, radii, 3, corner, &setback, 1, 0.1, &r, holes, 64, &len, nullptr,
                                   nullptr) == SMLIB_OK &&
          r != nullptr && len == 1 + 3 * 12 && holes[0] == 12);
    if (r != nullptr)
    {
        int64_t counts[3] = {0, 0, 0};
        CHECK(smlib_brep_counts(r, counts) == SMLIB_OK);
        CHECK(counts[0] == 9 && counts[1] == 30 && counts[2] == 22);
        int32_t closed = -1;
        CHECK(smlib_brep_is_manifold_solid(r, &closed) == SMLIB_OK && closed == 0);
        smlib_flat *f = nullptr;
        CHECK(smlib_brep_export(r, &f) == SMLIB_OK);
        const int64_t *ev = nullptr;
        const int64_t *eu = nullptr;
        const double *vp = nullptr;
        int64_t n_ev = 0, n_eu = 0, n_vp = 0;
        CHECK(smlib_flat_get_i64(f, SMLIB_FLAT_EDGE_VERTICES, &ev, &n_ev) == SMLIB_OK && n_ev == 2 * 30);
        CHECK(smlib_flat_get_i64(f, SMLIB_FLAT_EDGEUSE_EDGE, &eu, &n_eu) == SMLIB_OK);
        CHECK(smlib_flat_get_f64(f, SMLIB_FLAT_VERTEX_POSITION, &vp, &n_vp) == SMLIB_OK && n_vp == 3 * 22);
        // The edges history says each edge bounds: the hole's own edges bound their one face.
        int64_t hlen = 0;
        CHECK(smlib_history_read(r, nullptr, 0, &hlen) == SMLIB_ERR_BUFFER);
        std::vector<int64_t> h(static_cast<size_t>(hlen));
        CHECK(smlib_history_read(r, h.data(), hlen, &hlen) == SMLIB_OK);
        std::vector<std::vector<int64_t>> bounds;
        size_t at = 3;
        for (int64_t i = 0; i < h[1]; ++i)
        {
            at += 2 + 2 * static_cast<size_t>(h[at + 1]);
        }
        for (int64_t i = 0; i < h[2]; ++i)
        {
            std::vector<int64_t> fs(h.begin() + static_cast<long>(at) + 1, h.begin() + static_cast<long>(at) + 1 + h[at]);
            bounds.push_back(fs);
            at += 1 + static_cast<size_t>(h[at]);
            at += 1 + 2 * static_cast<size_t>(h[at]);
        }
        int64_t free_edges = 0;
        for (const auto &b : bounds)
        {
            free_edges += b.size() == 1 ? 1 : 0;
        }
        CHECK(free_edges == 12);
        const double ends[6][3] = {{52, 35, 20}, {52, 40, 15}, {55, 32, 20}, {60, 32, 15}, {55, 40, 12}, {60, 35, 12}};
        int found = 0;
        for (int64_t k = 0; k < 12; ++k)
        {
            int64_t e = holes[1 + 3 * k], face = holes[2 + 3 * k], fwd = holes[3 + 3 * k];
            int64_t n = holes[1 + 3 * ((k + 1) % 12)], nfwd = holes[3 + 3 * ((k + 1) % 12)];
            CHECK(e >= 0 && e < 30 && face >= 0 && face < 9 && (fwd == 0 || fwd == 1));
            CHECK(bounds[static_cast<size_t>(e)].size() == 1 && bounds[static_cast<size_t>(e)][0] == face);
            int64_t end = ev[2 * e + (fwd != 0 ? 1 : 0)];
            int64_t next_start = ev[2 * n + (nfwd != 0 ? 0 : 1)];
            CHECK(end == next_start);
            for (const auto &p : ends)
            {
                double d = std::hypot(vp[3 * end] - p[0], vp[3 * end + 1] - p[1], vp[3 * end + 2] - p[2]);
                found += d <= 1e-6 ? 1 : 0;
            }
        }
        CHECK(found == 6);
        smlib_flat_free(f);
        done(r);
    }
    // Every edge rounded and all eight corners set back 7.5: eight holes, each round its own
    // corner (every vertex of it within the setback's reach of that corner, none of another's).
    const double all_edges[36] = {0, 0, 10, 60, 0, 10, 0, 40, 10, 60, 40, 10, 30, 0, 0, 30, 40, 0,
                                  30, 0, 20, 30, 40, 20, 0, 20, 0, 60, 20, 0, 0, 20, 20, 60, 20, 20};
    const double all_radii[12] = {5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5};
    double corners[24];
    double backs[8];
    for (int k = 0; k < 8; ++k)
    {
        corners[3 * k] = (k & 1) != 0 ? 60 : 0;
        corners[3 * k + 1] = (k & 2) != 0 ? 40 : 0;
        corners[3 * k + 2] = (k & 4) != 0 ? 20 : 0;
        backs[k] = 7.5;
    }
    int64_t all_holes[8 * 64];
    CHECK(smlib_blend_setback_open(part, SMLIB_BLEND_ROUND, all_edges, all_radii, 12, corners, backs, 8, 0.1, &r, all_holes, 8 * 64, &len,
                                   nullptr, nullptr) == SMLIB_OK &&
          r != nullptr && len == 8 * (1 + 3 * 12));
    if (r != nullptr)
    {
        smlib_flat *f = nullptr;
        CHECK(smlib_brep_export(r, &f) == SMLIB_OK);
        const int64_t *ev = nullptr;
        const double *vp = nullptr;
        int64_t n_ev = 0, n_vp = 0;
        CHECK(smlib_flat_get_i64(f, SMLIB_FLAT_EDGE_VERTICES, &ev, &n_ev) == SMLIB_OK);
        CHECK(smlib_flat_get_f64(f, SMLIB_FLAT_VERTEX_POSITION, &vp, &n_vp) == SMLIB_OK);
        std::vector<int64_t> seen;
        for (int k = 0; k < 8; ++k)
        {
            const int64_t *h = all_holes + k * (1 + 3 * 12);
            CHECK(h[0] == 12);
            double farthest = 0.0;
            for (int64_t i = 0; i < 12; ++i)
            {
                int64_t e = h[1 + 3 * i];
                seen.push_back(e);
                for (int end = 0; end < 2; ++end)
                {
                    const double *q = vp + 3 * ev[2 * e + end];
                    farthest = std::max(farthest, std::hypot(q[0] - corners[3 * k], q[1] - corners[3 * k + 1], q[2] - corners[3 * k + 2]));
                }
            }
            CHECK(farthest < 7.5 * std::sqrt(2.0));
        }
        std::sort(seen.begin(), seen.end());
        CHECK(std::unique(seen.begin(), seen.end()) == seen.end());
        smlib_flat_free(f);
        done(r);
    }
    smlib_brep_free(part);
}

// Every section half round an open set-back corner of `r` (`holes`, as smlib_blend_setback_open
// writes them, `n_corners` holes) is the exact arc of its round: a rational B-spline lying on the
// round within the body's tolerance, its tangent in the round's tangent plane. The sides on the
// planar faces are the faces' own curves, not looked at here. Returns how many halves it saw.
int sections_on_their_rounds(const smlib_brep *r, const int64_t *holes, int64_t n_corners, const char *what)
{
    SmBrep *b = r->p;
    SmTArray<SmEdge *> edges;
    SmTArray<SmFace *> faces;
    b->GetEdges(edges);
    b->GetFaces(faces);
    double tol = static_cast<double>(b->GetTolerance());
    int halves = 0;
    const int64_t *h = holes;
    for (int64_t k = 0; k < n_corners; ++k)
    {
        int64_t n = h[0];
        for (int64_t i = 0; i < n; ++i)
        {
            int64_t e = h[1 + 3 * i];
            int64_t f = h[2 + 3 * i];
            bool indexed = e >= 0 && static_cast<ULONG>(e) < edges.GetSize() && f >= 0 && static_cast<ULONG>(f) < faces.GetSize();
            CHECK(indexed);
            if (!indexed)
            {
                continue;
            }
            SmSurface *round = faces[static_cast<ULONG>(f)]->GetSurface();
            if (round == nullptr || round->IsPlanar(tol))
            {
                continue;
            }
            ++halves;
            SmEdge *edge = edges[static_cast<ULONG>(e)];
            SmCurve *c = edge->GetCurve();
            bool rational = c != nullptr && c->IsKindOf(SmBSplineCurve_TYPE) && c->IsRational();
            double off = 0.0;
            double tilt = 0.0;
            if (c != nullptr)
            {
                SmExtent1d ivl = edge->GetInterval();
                for (int s = 0; s <= 16; ++s)
                {
                    double t = ivl.GetMin() + (ivl.GetMax() - ivl.GetMin()) * s / 16.0;
                    SmVector3d pd[2];
                    CHECK(c->Evaluate(t, 1, TRUE, pd) == SM_SUCCESS);
                    SmSolutionArray sol;
                    CHECK(round->GlobalPointSolve(round->GetNaturalUVDomain(), SM_SO_MINIMIZE, SmPoint3d(pd[0]), SmZoneTol3d(1e-12), nullptr,
                                                  SM_SR_SINGLE, sol) == SM_SUCCESS &&
                          sol.GetSize() > 0);
                    if (sol.GetSize() == 0)
                    {
                        off = HUGE_VAL;
                        continue;
                    }
                    SmPoint2d uv(sol[0].m_vStart[0], sol[0].m_vStart[1]);
                    SmPoint3d q;
                    SmVector3d normal;
                    CHECK(round->EvaluatePoint(uv, q) == SM_SUCCESS && round->EvaluateNormal(uv, TRUE, TRUE, normal) == SM_SUCCESS);
                    off = std::max(off, (SmPoint3d(pd[0]) - q).Length());
                    double speed = pd[1].Length();
                    if (speed > 0.0)
                    {
                        tilt = std::max(tilt, std::asin(std::min(1.0, std::fabs(pd[1].Dot(normal)) / speed)));
                    }
                }
            }
            // The tangent's tilt out of the plane, no more than a point within `tol` of the
            // round can give it over the half's length.
            double length = (edge->GetStartVertex()->GetPoint() - edge->GetEndVertex()->GetPoint()).Length();
            bool ok = rational && off <= tol && tilt <= std::max(1e-9, 10.0 * tol / std::max(length, tol));
            if (!ok)
            {
                std::fprintf(stderr,
                             "  %s, hole %lld: section half (edge %lld, face %lld) rational %d degree %lu, %.3g off the round, tangent %.3g deg "
                             "out of its plane (tolerance %.3g)\n",
                             what, static_cast<long long>(k), static_cast<long long>(e), static_cast<long long>(f), rational ? 1 : 0,
                             c != nullptr && c->IsKindOf(SmBSplineCurve_TYPE)
                                 ? static_cast<unsigned long>(static_cast<SmBSplineCurve *>(c)->GetDegree())
                                 : 0ul,
                             off, tilt * 180.0 / PI, tol);
            }
            CHECK(ok);
        }
        h += 1 + 3 * n;
    }
    return halves;
}

// tangency-labs/solidcore#1377: a 40 x 40 x 40 box with every edge rounded r 4 and every corner
// set back 7, left open, and its corner (0, 0, 0) alone: each of the six section halves round
// every hole is its round's exact arc, never a polynomial spline fitted near it.
void setback_sections_are_exact_arcs()
{
    smlib_brep *part = box(0, 0, 0, 40, 40, 40);
    double edges[36];
    for (int k = 0; k < 4; ++k)
    {
        double a = (k & 1) != 0 ? 40 : 0;
        double c = (k & 2) != 0 ? 40 : 0;
        const double along_x[3] = {20, a, c};
        const double along_y[3] = {a, 20, c};
        const double along_z[3] = {a, c, 20};
        for (int j = 0; j < 3; ++j)
        {
            edges[9 * k + j] = along_x[j];
            edges[9 * k + 3 + j] = along_y[j];
            edges[9 * k + 6 + j] = along_z[j];
        }
    }
    double radii[12];
    std::fill(radii, radii + 12, 4.0);
    double corners[24];
    double backs[8];
    for (int k = 0; k < 8; ++k)
    {
        corners[3 * k] = (k & 1) != 0 ? 40 : 0;
        corners[3 * k + 1] = (k & 2) != 0 ? 40 : 0;
        corners[3 * k + 2] = (k & 4) != 0 ? 40 : 0;
        backs[k] = 7.0;
    }
    for (int64_t n_corners : {8, 1})
    {
        smlib_brep *r = nullptr;
        std::vector<int64_t> holes(static_cast<size_t>(8 * 64));
        int64_t len = -1;
        CHECK(smlib_blend_setback_open(part, SMLIB_BLEND_ROUND, edges, radii, 12, corners, backs, n_corners, 0.1, &r, holes.data(),
                                       static_cast<int64_t>(holes.size()), &len, nullptr, nullptr) == SMLIB_OK &&
              r != nullptr && len == n_corners * (1 + 3 * 12));
        if (r != nullptr)
        {
            CHECK(sections_on_their_rounds(r, holes.data(), n_corners, n_corners == 8 ? "eight corners" : "corner (0, 0, 0)") ==
                  6 * n_corners);
            done(r);
        }
    }
    smlib_brep_free(part);
}

// The inner corner of a 20 x 20 x 10 notch cut from the top corner (60, 40, 20) of the 60 x 40 x
// 20 box: its three concave edges, meeting at (40, 20, 10), rounded r5 and set back 7.5, left
// open (one hole of twelve edges, a sheet SMLib's check passes) and closed (SMLib's patch, a
// valid solid). At (60, 20, 10) two convex edges meet the concave one: a mixed corner, which
// SMLib sets back with the same-convexity topology and gets wrong, refused.
void concave_corner_setback()
{
    smlib_brep *part = box(0, 0, 0, 60, 40, 20);
    smlib_brep *block = box(40, 20, 10, 30, 30, 20);
    smlib_brep *notched = nullptr;
    CHECK(smlib_boolean(part, block, SMLIB_BOOLEAN_SUBTRACT, &notched) == SMLIB_OK && notched != nullptr);
    smlib_brep_free(part);
    smlib_brep_free(block);
    if (notched == nullptr)
    {
        return;
    }
    const double inner[9] = {50, 20, 10, 40, 30, 10, 40, 20, 15};
    const double radii[3] = {5, 5, 5};
    const double corner[3] = {40, 20, 10};
    const double setback = 7.5;
    smlib_brep *r = nullptr;
    int64_t holes[64];
    int64_t len = -1;
    CHECK(smlib_blend_setback_open(notched, SMLIB_BLEND_ROUND, inner, radii, 3, corner, &setback, 1, 0.1, &r, holes, 64, &len, nullptr,
                                   nullptr) == SMLIB_OK &&
          r != nullptr && len == 1 + 3 * 12 && holes[0] == 12);
    if (r != nullptr)
    {
        int32_t closed = -1;
        CHECK(smlib_brep_is_manifold_solid(r, &closed) == SMLIB_OK && closed == 0);
        done(r);
    }
    CHECK(smlib_blend_setback(notched, SMLIB_BLEND_ROUND, inner, radii, 3, corner, &setback, 1, 0.1, &r, nullptr, nullptr) == SMLIB_OK &&
          r != nullptr);
    if (r != nullptr)
    {
        done(r);
    }
    const double mixed[9] = {50, 20, 10, 60, 30, 10, 60, 20, 15};
    const double mixed_corner[3] = {60, 20, 10};
    CHECK(smlib_blend_setback(notched, SMLIB_BLEND_ROUND, mixed, radii, 3, mixed_corner, &setback, 1, 0.1, &r, nullptr, nullptr) ==
              SMLIB_ERR_UNSUPPORTED &&
          r == nullptr);
    CHECK(smlib_blend_setback_open(notched, SMLIB_BLEND_ROUND, mixed, radii, 3, mixed_corner, &setback, 1, 0.1, &r, holes, 64, &len, nullptr,
                                   nullptr) == SMLIB_ERR_UNSUPPORTED &&
          r == nullptr);
    // Both corners asked for, the mixed one second: the refusal names it by its point's index.
    const double both_edges[15] = {50, 20, 10, 40, 30, 10, 40, 20, 15, 60, 30, 10, 60, 20, 15};
    const double both_radii[5] = {5, 5, 5, 5, 5};
    const double both_corners[6] = {40, 20, 10, 60, 20, 10};
    const double both_setbacks[2] = {7.5, 7.5};
    for (int open = 0; open < 2; ++open)
    {
        int64_t failed_corner = -1;
        int32_t s = open == 1 ? smlib_blend_setback_open(notched, SMLIB_BLEND_ROUND, both_edges, both_radii, 5, both_corners, both_setbacks, 2, 0.1,
                                                         &r, holes, 64, &len, nullptr, &failed_corner)
                              : smlib_blend_setback(notched, SMLIB_BLEND_ROUND, both_edges, both_radii, 5, both_corners, both_setbacks, 2, 0.1,
                                                    &r, nullptr, &failed_corner);
        CHECK(s == SMLIB_ERR_UNSUPPORTED && r == nullptr && failed_corner == 1);
    }
    smlib_brep_free(notched);
}

// The vertex `p` lies on: its index in `vs`, or -1.
int64_t vertex_at(const SmTArray<SmVertex *> &vs, const SmPoint3d &p)
{
    for (ULONG k = 0; k < vs.GetSize(); ++k)
    {
        if ((vs[k]->GetPoint() - p).Length() <= 1e-6)
        {
            return static_cast<int64_t>(k);
        }
    }
    return -1;
}

// The curvature of `e`'s curve at its end on `v`.
double end_curvature(SmEdge *e, const SmVertex *v)
{
    SmExtent1d ivl = e->GetInterval();
    double t = e->GetStartVertex() == v ? ivl.GetMin() : ivl.GetMax();
    SmVector3d d[3];
    if (e->GetCurve() == nullptr || e->GetCurve()->Evaluate(t, 2, TRUE, d) != SM_SUCCESS)
    {
        return HUGE_VAL;
    }
    double speed = d[1].Length();
    return (d[1] * d[2]).Length() / (speed * speed * speed);
}

// The corner patch of a blended corner at `corner`: every B-spline face of `b` whose vertices
// all lie within `reach` of it. Returns the largest angle (degrees) between the faces either
// side of the patch's edges, at 31 stations inside each edge (its ends left out: a 3-sided
// patch's collapsed side is a point where it has no one normal), and sets `n_faces` to the
// patches found. A tangent seam reads 0.
double worst_corner_seam_angle(const smlib_brep *b, const SmPoint3d &corner, double reach, int &n_faces)
{
    n_faces = 0;
    double worst = 0.0;
    SmTArray<SmFace *> faces;
    b->p->GetFaces(faces);
    for (ULONG f = 0; f < faces.GetSize(); ++f)
    {
        if (faces[f]->GetSurface() == nullptr || !faces[f]->GetSurface()->IsKindOf(SmBSplineSurface_TYPE))
        {
            continue;
        }
        SmTArray<SmVertex *> vs;
        faces[f]->GetVertices(vs);
        bool near_corner = vs.GetSize() > 0;
        for (ULONG k = 0; k < vs.GetSize(); ++k)
        {
            near_corner = near_corner && (vs[k]->GetPoint() - corner).Length() <= reach;
        }
        if (!near_corner)
        {
            continue;
        }
        n_faces++;
        SmTArray<SmEdge *> es;
        faces[f]->GetEdges(es);
        for (ULONG k = 0; k < es.GetSize(); ++k)
        {
            SmEdgeuse *a = es[k]->GetPrimaryEdgeuse();
            SmEdgeuse *o = a != nullptr ? a->GetRadial() : nullptr;
            if (o == nullptr || o == a)
            {
                worst = HUGE_VAL;
                continue;
            }
            SmExtent1d ivl = es[k]->GetInterval();
            for (int s = 1; s < 32; ++s)
            {
                SmPoint3d pa, pb;
                SmVector3d ba, bb;
                if (a->EvaluateBinormal(ivl.Evaluate(s / 32.0), FALSE, pa, ba) != SM_SUCCESS ||
                    o->EvaluateBinormal(ivl.Evaluate(s / 32.0), FALSE, pb, bb) != SM_SUCCESS || ba.Unitize(FALSE) != SM_SUCCESS ||
                    bb.Unitize(FALSE) != SM_SUCCESS)
                {
                    worst = HUGE_VAL;
                    continue;
                }
                // tangent faces' binormals (into each face, square to the edge) are opposed
                worst = std::max(worst, std::acos(std::min(1.0, -ba.Dot(bb))) * 180.0 / PI);
            }
        }
    }
    return worst;
}

// The normal curvature, at the point of `b`'s non-planar faces nearest `p` (within 1e-6), of
// the section square to `along` (unit): the round's cross curvature there. HUGE_VAL if no
// such face.
double cross_curvature(const SmBrep *b, const SmPoint3d &p, const SmVector3d &along)
{
    SmTArray<SmFace *> faces;
    b->GetFaces(faces);
    for (ULONG i = 0; i < faces.GetSize(); ++i)
    {
        SmSurface *s = faces[i]->GetSurface();
        if (s == nullptr || SM_CAST_PTR(SmPlane, s) != nullptr)
        {
            continue;
        }
        SmSolutionArray sols;
        if (s->GlobalPointSolve(s->GetNaturalUVDomain(), SM_SO_MINIMIZE, p, 1e-7, nullptr, SM_SR_SINGLE, sols) != SM_SUCCESS ||
            sols.GetSize() < 1)
        {
            continue;
        }
        SmPoint2d uv(sols[0].m_vStart[0], sols[0].m_vStart[1]);
        SmPoint3d at;
        SmVector3d su, sv, suv, suu, svv;
        if (s->Evaluate2ndDerivatives(uv, TRUE, TRUE, at, su, sv, suv, suu, svv) != SM_SUCCESS || (at - p).Length() > 1e-6)
        {
            continue;
        }
        SmVector3d n = su * sv;
        n = n / n.Length();
        SmVector3d t = along * n;
        t = t / t.Length();
        double e = su.Dot(su), f = su.Dot(sv), g = sv.Dot(sv);
        double det = e * g - f * f;
        double a = (g * t.Dot(su) - f * t.Dot(sv)) / det;
        double c = (e * t.Dot(sv) - f * t.Dot(su)) / det;
        return std::fabs(a * a * suu.Dot(n) + 2.0 * a * c * suv.Dot(n) + c * c * svv.Dot(n));
    }
    return HUGE_VAL;
}

// The corner (40, 40, 40) of a 40 mm box, its three edges rounded r 8 and set back 12 with an
// ease of 4 (smlib_blend_setback_open_eased; SolidCraft ladder 07/015 with corner: g2). Over
// the 4 before the setback each round's section eases from the circle into the curvature-
// continuous quintic: on the top face's contact line of the round along X, at 16 and more from
// the vertex the cross curvature is the circle's 1/8, at 12 (the setback) 0, falling in between;
// across where the ease starts the round is curvature-continuous. The hole's six corners (the
// sections' ends) each meet a section half and a face-curve half, both with zero curvature
// there. The sheet passes SMLib's check and its only free edges are the hole's. An ease of 0 is
// today's body exactly; the ease is refused where it does not fit or is no length.
// The corner (60, 40, 20) of a 60 x 40 x 20 box, its three edges rounded r5 and set back 12 or
// 13 with an ease of 0.5, 1.25 or 5 (SolidCraft's G2 corner, its long setbacks). Each face curve
// between two sections is split at the side's middle, and both halves are merged into the face
// they lie on. One half's midpoint lay on the corner's face curve to scaled zero and the other's
// only to the edges' tolerance, so SmFilletEdge::FindFilletBrepEdge kept the first and dropped
// the second: that half was never merged, and the sheet came back with two free edges besides
// the hole's (SMLIB_ERR_NOT_SOLID). The sheet now holds one hole of twelve edges for every pair.
void corner_setback_eased_long_setbacks()
{
    smlib_brep *part = box(0, 0, 0, 60, 40, 20);
    const double edges[9] = {30, 40, 20, 60, 20, 20, 60, 40, 10};
    const double radii[3] = {5, 5, 5};
    const double corner[3] = {60, 40, 20};
    for (double setback : {12.0, 13.0})
    {
        for (double ease : {0.5, 1.25, 5.0})
        {
            int64_t holes[64];
            int64_t len = -1;
            smlib_brep *r = nullptr;
            const int32_t s = smlib_blend_setback_open_eased(part, SMLIB_BLEND_ROUND, edges, radii, 3, corner, &setback, &ease, 1, 0.1, &r,
                                                             holes, 64, &len, nullptr, nullptr);
            if (s != SMLIB_OK)
            {
                std::printf("set back %g, eased %g: status %d\n", setback, ease, static_cast<int>(s));
            }
            CHECK(s == SMLIB_OK && r != nullptr && len == 1 + 3 * 12 && holes[0] == 12);
            smlib_brep_free(r);
        }
    }
    smlib_brep_free(part);
}

// The unit tangent and the curvature vector of `e`'s curve at its end on `v`, as the curve
// leaves `v` (the tangent pointing away from it). False if it cannot be evaluated.
bool leaving(SmEdge *e, const SmVertex *v, SmVector3d &tangent, SmVector3d &curvature)
{
    SmExtent1d ivl = e->GetInterval();
    bool at_start = e->GetStartVertex() == v;
    SmVector3d d[3];
    if (e->GetCurve() == nullptr || e->GetCurve()->Evaluate(at_start ? ivl.GetMin() : ivl.GetMax(), 2, TRUE, d) != SM_SUCCESS)
    {
        return false;
    }
    double speed = d[1].Length();
    if (!(speed > 0.0))
    {
        return false;
    }
    tangent = d[1] / (at_start ? speed : -speed);
    // the curvature vector, the second derivative in arc length: square to the tangent, the
    // same whichever way the curve runs
    curvature = (d[2] - d[1] * (d[1].Dot(d[2]) / (speed * speed))) / (speed * speed);
    return true;
}

// The corner (60, 40, 20) of a 60 x 40 x 20 box, its three edges rounded r5, set back 8 to 15
// and eased 0.5 to 5 (smlib_blend_setback_open_eased; SolidCraft's explicit G2 corner patch, a
// biquintic net per patch, holds a hole edge exactly only as one quintic span). Each face curve
// between two sections is split at its middle, and each half is one Bezier span of degree 5 at
// most: the half of the side's own quintic, not the corner patch's iso-curve (degree 6 with
// interior knots, from the patch's cross-derivative fields). At the hole's corner a half meets
// the round's contact line on that face curvature-continuously (tangent along it, no curvature
// in the plane, as the straight contact line), and the two halves meet so at the middle. The
// sheet is valid and its only free edges are the hole's; closed by SMLib's patch the body is a
// valid closed solid.
void corner_setback_eased_face_curves_are_quintics()
{
    smlib_brep *part = box(0, 0, 0, 60, 40, 20);
    const double edges[9] = {30, 40, 20, 60, 20, 20, 60, 40, 10};
    const double radii[3] = {5, 5, 5};
    const double corner[3] = {60, 40, 20};
    for (double setback : {8.0, 10.0, 12.0, 13.0, 15.0})
    {
        for (double ease : {0.5, 1.25, 3.0, 5.0})
        {
            if (setback + ease >= 20.0)
            {
                continue; // the 20 edge has no room (SMLIB_ERR_NO_ROOM)
            }
            int64_t holes[64];
            int64_t len = -1;
            smlib_brep *r = nullptr;
            const int32_t s = smlib_blend_setback_open_eased(part, SMLIB_BLEND_ROUND, edges, radii, 3, corner, &setback, &ease, 1, 0.1, &r,
                                                             holes, 64, &len, nullptr, nullptr);
            CHECK(s == SMLIB_OK && r != nullptr && len == 1 + 3 * 12 && holes[0] == 12);
            if (s != SMLIB_OK || r == nullptr || len != 1 + 3 * 12)
            {
                std::printf("set back %g, eased %g: status %d\n", setback, ease, static_cast<int>(s));
                smlib_brep_free(r);
                continue;
            }
            SmTArray<SmEdge *> es;
            SmTArray<SmFace *> fs;
            r->p->GetEdges(es);
            r->p->GetFaces(fs);
            std::vector<SmEdge *> hole;
            for (int64_t k = 0; k < 12; ++k)
            {
                hole.push_back(es[static_cast<ULONG>(holes[1 + 3 * k])]);
            }
            int face_curves = 0, not_quintic = 0, met_rail = 0, met_half = 0;
            double worst_turn = 0.0, worst_bend = 0.0;
            for (int64_t k = 0; k < 12; ++k)
            {
                SmEdge *h = hole[static_cast<size_t>(k)];
                SmFace *f = fs[static_cast<ULONG>(holes[2 + 3 * k])];
                if (f->GetSurface() == nullptr || SM_CAST_PTR(SmPlane, f->GetSurface()) == nullptr)
                {
                    continue; // a section, on its round
                }
                ++face_curves;
                const SmBSplineCurve *c = SM_CAST_PTR(SmBSplineCurve, h->GetCurve());
                ULONG degree = c != nullptr ? c->GetDegree() : 0;
                ULONG spans = c != nullptr ? c->GetNumberOfUniqueKnots() - 1 : 0;
                if (c == nullptr || c->IsRational() || degree > 5 || spans != 1)
                {
                    ++not_quintic;
                    std::printf("set back %g, eased %g: a face curve of degree %lu, %lu spans, %lu points\n", setback, ease,
                                static_cast<unsigned long>(degree), static_cast<unsigned long>(spans),
                                static_cast<unsigned long>(c != nullptr ? c->GetNumberControlPoints() : 0));
                }
                // at each end, the other edge on this face there: the contact line (at the hole's
                // corner) or the face curve's other half (at the middle)
                for (SmVertex *v : {h->GetStartVertex(), h->GetEndVertex()})
                {
                    SmTArray<SmEdge *> at;
                    v->GetEdges(at);
                    for (ULONG j = 0; j < at.GetSize(); ++j)
                    {
                        SmEdge *o = at[j];
                        SmTArray<SmFace *> beside;
                        o->GetFaces(beside);
                        ULONG where = 0;
                        if (o == h || !beside.FindElement(f, where))
                        {
                            continue;
                        }
                        SmVector3d t1, k1, t2, k2;
                        CHECK(leaving(h, v, t1, k1) && leaving(o, v, t2, k2));
                        // continuing across v: opposed tangents, one curvature vector
                        worst_turn = std::max(worst_turn, (t1 + t2).Length());
                        worst_bend = std::max(worst_bend, (k1 - k2).Length());
                        bool half = std::find(hole.begin(), hole.end(), o) != hole.end();
                        ++(half ? met_half : met_rail);
                        if (!half)
                        {
                            // the contact line is straight: the face curve leaves it unturned
                            CHECK(k1.Length() <= 1e-9);
                        }
                    }
                }
            }
            CHECK(face_curves == 6 && not_quintic == 0);
            CHECK(met_rail == 6 && met_half == 6);
            if (worst_turn > 1e-9 || worst_bend > 1e-9)
            {
                std::printf("set back %g, eased %g: face curves meet their neighbours %.3g off in tangent, %.3g in curvature\n", setback, ease,
                            worst_turn, worst_bend);
            }
            CHECK(worst_turn <= 1e-9 && worst_bend <= 1e-9);
            // the sheet's free edges are exactly the hole's
            int free_edges = 0;
            for (ULONG k = 0; k < es.GetSize(); ++k)
            {
                SmTArray<SmFace *> beside;
                es[k]->GetFaces(beside);
                if (beside.GetSize() == 1)
                {
                    CHECK(std::find(hole.begin(), hole.end(), es[k]) != hole.end());
                    ++free_edges;
                }
            }
            CHECK(free_edges == 12);
            done(r);
            CHECK(smlib_blend_setback_eased(part, SMLIB_BLEND_ROUND, edges, radii, 3, corner, &setback, &ease, 1, 0.1, &r, nullptr, nullptr) ==
                      SMLIB_OK &&
                  r != nullptr);
            if (r != nullptr)
            {
                int32_t closed = -1;
                CHECK(smlib_brep_is_manifold_solid(r, &closed) == SMLIB_OK && closed == 1);
                done(r);
            }
        }
    }
    smlib_brep_free(part);
}

void corner_setback_eased()
{
    smlib_brep *part = box(0, 0, 0, 40, 40, 40);
    const double edges[9] = {20, 40, 40, 40, 20, 40, 40, 40, 20};
    const double radii[3] = {8, 8, 8};
    const double corner[3] = {40, 40, 40};
    const double setback = 12.0;
    const double ease = 4.0;
    const double no_ease = 0.0;
    // the setback sections' ends, the hole's six corners
    const double ends[6][3] = {{28, 32, 40}, {28, 40, 32}, {32, 28, 40}, {40, 28, 32}, {32, 40, 28}, {40, 32, 28}};
    int64_t holes[64], holes0[64];
    int64_t len = -1, len0 = -1;
    smlib_brep *r = nullptr;
    smlib_brep *r0 = nullptr;

    // ease 0: the same sheet and holes as without one, and the same closed body
    CHECK(smlib_blend_setback_open(part, SMLIB_BLEND_ROUND, edges, radii, 3, corner, &setback, 1, 0.1, &r0, holes0, 64, &len0, nullptr,
                                   nullptr) == SMLIB_OK);
    CHECK(smlib_blend_setback_open_eased(part, SMLIB_BLEND_ROUND, edges, radii, 3, corner, &setback, &no_ease, 1, 0.1, &r, holes, 64, &len,
                                         nullptr, nullptr) == SMLIB_OK);
    if (r != nullptr && r0 != nullptr)
    {
        CHECK(len == len0 && std::equal(holes, holes + len, holes0));
        smlib_flat *f = nullptr;
        smlib_flat *f0 = nullptr;
        CHECK(smlib_brep_export(r, &f) == SMLIB_OK && smlib_brep_export(r0, &f0) == SMLIB_OK);
        const double *vp = nullptr;
        const double *vp0 = nullptr;
        int64_t n_vp = 0, n_vp0 = 0;
        CHECK(smlib_flat_get_f64(f, SMLIB_FLAT_VERTEX_POSITION, &vp, &n_vp) == SMLIB_OK);
        CHECK(smlib_flat_get_f64(f0, SMLIB_FLAT_VERTEX_POSITION, &vp0, &n_vp0) == SMLIB_OK);
        CHECK(n_vp == n_vp0 && std::equal(vp, vp + n_vp, vp0));
        smlib_flat_free(f);
        smlib_flat_free(f0);
    }
    smlib_brep_free(r);
    smlib_brep_free(r0);
    r = r0 = nullptr;
    CHECK(smlib_blend_setback(part, SMLIB_BLEND_ROUND, edges, radii, 3, corner, &setback, 1, 0.1, &r0, nullptr, nullptr) == SMLIB_OK);
    CHECK(smlib_blend_setback_eased(part, SMLIB_BLEND_ROUND, edges, radii, 3, corner, &setback, &no_ease, 1, 0.1, &r, nullptr, nullptr) ==
          SMLIB_OK);
    if (r != nullptr && r0 != nullptr)
    {
        CHECK(volume_of(r) == volume_of(r0));
    }
    smlib_brep_free(r);
    smlib_brep_free(r0);
    r = r0 = nullptr;

    // ease 4, open
    CHECK(smlib_blend_setback_open_eased(part, SMLIB_BLEND_ROUND, edges, radii, 3, corner, &setback, &ease, 1, 0.1, &r, holes, 64, &len,
                                         nullptr, nullptr) == SMLIB_OK &&
          r != nullptr && len == 1 + 3 * 12 && holes[0] == 12);
    if (r != nullptr)
    {
        SmTArray<SmEdge *> es;
        SmTArray<SmVertex *> vs;
        r->p->GetEdges(es);
        r->p->GetVertices(vs);
        // the sheet's free edges are exactly the hole's
        std::vector<int64_t> free_edges, hole_edges;
        for (ULONG k = 0; k < es.GetSize(); ++k)
        {
            SmTArray<SmFace *> beside;
            es[k]->GetFaces(beside);
            if (beside.GetSize() == 1)
            {
                free_edges.push_back(static_cast<int64_t>(k));
            }
        }
        for (int64_t k = 0; k < 12; ++k)
        {
            hole_edges.push_back(holes[1 + 3 * k]);
        }
        std::sort(hole_edges.begin(), hole_edges.end());
        CHECK(free_edges == hole_edges);
        // at each of the hole's six corners a section and a face curve, both flat across there
        int met = 0;
        double worst = 0.0;
        for (const auto &e : ends)
        {
            int64_t v = vertex_at(vs, SmPoint3d(e[0], e[1], e[2]));
            CHECK(v >= 0);
            for (int64_t k = 0; k < 12 && v >= 0; ++k)
            {
                SmEdge *h = es[static_cast<ULONG>(holes[1 + 3 * k])];
                if (h->GetStartVertex() == vs[static_cast<ULONG>(v)] || h->GetEndVertex() == vs[static_cast<ULONG>(v)])
                {
                    worst = std::max(worst, end_curvature(h, vs[static_cast<ULONG>(v)]));
                    ++met;
                }
            }
        }
        std::fprintf(stderr, "  eased setback: worst curvature at the hole's corners %.3g\n", worst);
        CHECK(met == 12 && worst <= 1e-6);
        // the round along X, on its top contact line y 32, z 40: the circle's curvature up to
        // the ease, none at the setback, falling between, continuous where the ease starts
        const SmVector3d x(1, 0, 0);
        double k[9];
        const double at[9] = {12, 13, 14, 15, 16 - 1e-4, 16 + 1e-4, 18, 20, 24};
        for (int i = 0; i < 9; ++i)
        {
            k[i] = cross_curvature(r->p, SmPoint3d(40 - at[i], 32, 40), x);
        }
        std::fprintf(stderr, "  eased setback: cross curvature at 12..24: %.3g %.3g %.3g %.3g %.6g %.6g %.6g %.6g %.6g\n", k[0], k[1], k[2],
                     k[3], k[4], k[5], k[6], k[7], k[8]);
        CHECK(k[0] <= 1e-6);
        CHECK(k[0] < k[1] && k[1] < k[2] && k[2] < k[3] && k[3] < k[4]);
        CHECK(std::fabs(k[4] - 0.125) <= 1e-4 && std::fabs(k[5] - 0.125) <= 1e-4);
        for (int i = 6; i < 9; ++i)
        {
            CHECK(std::fabs(k[i] - 0.125) <= 1e-4);
        }
        done(r);
    }

    // ease 4, closed by SMLib's patch: a valid solid
    CHECK(smlib_blend_setback_eased(part, SMLIB_BLEND_ROUND, edges, radii, 3, corner, &setback, &ease, 1, 0.1, &r, nullptr, nullptr) ==
              SMLIB_OK &&
          r != nullptr);
    if (r != nullptr)
    {
        int32_t closed = -1;
        CHECK(smlib_brep_is_manifold_solid(r, &closed) == SMLIB_OK && closed == 1);
        done(r);
    }

    // the curvature-continuous blend: its sections are the quintic already, and eased its face
    // curves meet the contact lines curvature-continuously
    for (int32_t kind : {SMLIB_BLEND_CURVATURE})
    {
        CHECK(smlib_blend_setback_open_eased(part, kind, edges, radii, 3, corner, &setback, &ease, 1, 0.1, &r, holes, 64, &len, nullptr,
                                             nullptr) == SMLIB_OK &&
              r != nullptr && holes[0] == 12);
        if (r != nullptr)
        {
            SmTArray<SmEdge *> es;
            SmTArray<SmVertex *> vs;
            r->p->GetEdges(es);
            r->p->GetVertices(vs);
            double worst = 0.0;
            int met = 0;
            for (const auto &e : ends)
            {
                int64_t v = vertex_at(vs, SmPoint3d(e[0], e[1], e[2]));
                CHECK(v >= 0);
                for (int64_t k = 0; k < holes[0] && v >= 0; ++k)
                {
                    SmEdge *h = es[static_cast<ULONG>(holes[1 + 3 * k])];
                    if (h->GetStartVertex() == vs[static_cast<ULONG>(v)] || h->GetEndVertex() == vs[static_cast<ULONG>(v)])
                    {
                        worst = std::max(worst, end_curvature(h, vs[static_cast<ULONG>(v)]));
                        ++met;
                    }
                }
            }
            CHECK(met == 12);
            std::fprintf(stderr, "  eased setback, kind %d: worst curvature at the hole's corners %.3g\n", static_cast<int>(kind), worst);
            CHECK(worst <= 1e-6);
            done(r);
        }
    }

    // every edge rounded and every corner set back 12 with an ease of 3: each round eases at
    // both of its ends (12 + 3 twice is 30, inside the 40 edge), eight holes. (An ease of 4 on
    // this model is SMLIB_ERR_NOT_SOLID: the merge's stitch glues the two ends of the bottom
    // face curve at (40, 40, 0), 5.66 apart; with 3 or 6 it stitches. Open.)
    const double ease3 = 3.0;
    const double all_edges[36] = {20, 0, 0, 20, 40, 0, 20, 0, 40, 20, 40, 40, 0, 20, 0, 40, 20, 0,
                                  0, 20, 40, 40, 20, 40, 0, 0, 20, 40, 0, 20, 0, 40, 20, 40, 40, 20};
    const double all_radii[12] = {8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8};
    double all_corners[24];
    double all_backs[8];
    double all_eases[8];
    for (int k = 0; k < 8; ++k)
    {
        all_corners[3 * k] = (k & 1) != 0 ? 40 : 0;
        all_corners[3 * k + 1] = (k & 2) != 0 ? 40 : 0;
        all_corners[3 * k + 2] = (k & 4) != 0 ? 40 : 0;
        all_backs[k] = setback;
        all_eases[k] = ease3;
    }
    int64_t all_holes[8 * 64];
    CHECK(smlib_blend_setback_open_eased(part, SMLIB_BLEND_ROUND, all_edges, all_radii, 12, all_corners, all_backs, all_eases, 8, 0.1, &r,
                                         all_holes, 8 * 64, &len, nullptr, nullptr) == SMLIB_OK &&
          r != nullptr && len == 8 * (1 + 3 * 12));
    if (r != nullptr)
    {
        // on the round along X at y 40, z 40, the top contact line: flat across at both
        // setbacks, the circle's curvature in the middle
        const SmVector3d x(1, 0, 0);
        double k0 = cross_curvature(r->p, SmPoint3d(12, 32, 40), x);
        double k1 = cross_curvature(r->p, SmPoint3d(20, 32, 40), x);
        double k2 = cross_curvature(r->p, SmPoint3d(28, 32, 40), x);
        std::fprintf(stderr, "  eased setbacks at both ends: cross curvature %.3g %.6g %.3g\n", k0, k1, k2);
        CHECK(k0 <= 1e-6 && std::fabs(k1 - 0.125) <= 1e-4 && k2 <= 1e-6);
        done(r);
    }

    // refused: no length, too long for the edge (12 + 28 is the whole edge), no eases
    const double negative = -1.0;
    const double too_long = 28.0;
    CHECK(smlib_blend_setback_open_eased(part, SMLIB_BLEND_ROUND, edges, radii, 3, corner, &setback, &negative, 1, 0.1, &r, holes, 64, &len,
                                         nullptr, nullptr) == SMLIB_ERR_ARG &&
          r == nullptr);
    CHECK(smlib_blend_setback_eased(part, SMLIB_BLEND_ROUND, edges, radii, 3, corner, &setback, &too_long, 1, 0.1, &r, nullptr, nullptr) ==
              SMLIB_ERR_NO_ROOM &&
          r == nullptr);
    CHECK(smlib_blend_setback_eased(part, SMLIB_BLEND_ROUND, edges, radii, 3, corner, &setback, nullptr, 1, 0.1, &r, nullptr, nullptr) ==
              SMLIB_ERR_ARG &&
          r == nullptr);
    smlib_brep_free(part);
}

void shelled_tapered_box_rounds()
{
    const std::vector<int64_t> kinds = {SMLIB_CURVE_LINE, SMLIB_CURVE_LINE, SMLIB_CURVE_LINE, SMLIB_CURVE_LINE};
    const std::vector<double> reals = {-60, -40, 0, 60, -40, 0, 60, -40, 0, 60, 40, 0,
                                       60, 40, 0, -60, 40, 0, -60, 40, 0, -60, -40, 0};
    const double up[3] = {0, 0, 1};
    smlib_brep *box5 = nullptr, *round = nullptr, *shell = nullptr;
    CHECK(smlib_extrude(kinds.data(), 4, reals.data(), static_cast<int64_t>(reals.size()), up, 40.0, 5.0 * PI / 180.0,
                        SMLIB_DRAFT_CORNER_CURVATURE, 1, &box5, nullptr) == SMLIB_OK);
    if (box5 == nullptr)
    {
        return;
    }
    const double in = 20.0 * std::tan(5.0 * PI / 180.0); // the sloping edges' midpoints, at z 20
    const double edges[12] = {60 - in, 40 - in, 20, 60 - in, -(40 - in), 20, -(60 - in), 40 - in, 20, -(60 - in), -(40 - in), 20};
    CHECK(smlib_fillet(box5, edges, 4, 0.1, 4.0, &round, nullptr) == SMLIB_OK);
    if (round != nullptr)
    {
        const double top[3] = {0, 0, 40};
        CHECK(smlib_shell(round, top, 1, 0.1, -2.0, &shell, nullptr) == SMLIB_OK);
        CHECK(shell != nullptr && sound(shell));
    }
    smlib_brep_free(shell);
    smlib_brep_free(round);
    smlib_brep_free(box5);
}

// A dome (a quarter disc r 20 revolved about Z: a sphere and a plane) shelled 2 inward, open at
// its base (SolidCraft rig shell/4847). The offset meets the base plane at the inner sphere's
// equator, the circle r 18: its edge is that circle exactly, and the shell's volume is
// (2/3) pi (20^3 - 18^3). The offset used to intersect the two faces' tangent (G1R) B-spline
// extensions, marching a cubic within 6.5e-7 of the circle; SolidCraft read that back as a circle
// off the axis, whose seam then sat 1.1e-7 from the dome's in the union that joins the rim back in.
// Three R5 rounds at the corner (60, 40, 20) of a 60 x 40 x 20 box, the corner set back 8 along
// each edge (smlib_blend_setback). Each round runs whole to its cross-section 8 from the vertex,
// and the patch between those sections meets each face in a curve between the two rounds'
// contact lines: the six vertices where the patch meets the faces are those sections' ends
// (on the top face (52, 35, 20) and (55, 32, 20), and so on). Outside the corner cube C =
// [52, 60] x [32, 40] x [12, 20] the body is that of the corner without a setback, the sphere
// corner (ladder 07-fillets/009: V0 = 47377.120025868, of which 512 - 225 (1 - pi/4) -
// 125 (1 - pi/6) in C). The body is convex, so inside C it holds the hull of its points on C's
// faces (the three sections and the box's corners outside the patch: 299.0201, by a hull of
// them sampled, below the true hull) and lies within every plane tangent to it along the
// sections, the three rounds' prisms: 512 - 225 (1 - pi/4) - 125 (sqrt 2 - 1), the corner cube
// of r keeping the three cylinders' common octant, (2 - sqrt 2) r^3, in place of the sphere's.
void corner_setback()
{
    smlib_brep *part = box(0, 0, 0, 60, 40, 20);
    const double edges[9] = {30, 40, 20, 60, 20, 20, 60, 40, 10};
    const double radii[3] = {5, 5, 5};
    const double corner[3] = {60, 40, 20};
    const double setback = 8.0;
    smlib_brep *r = nullptr;
    int64_t failed = -1;
    CHECK(smlib_blend_setback(part, SMLIB_BLEND_ROUND, edges, radii, 3, corner, &setback, 1, 0.1, &r, nullptr, &failed) == SMLIB_OK);
    if (r != nullptr)
    {
        const double in_c_sphere = 512.0 - 225.0 * (1.0 - PI / 4.0) - 125.0 * (1.0 - PI / 6.0);
        const double outside = 47377.120025868 - in_c_sphere;
        const double hull = 299.0201;
        const double prisms = 512.0 - 225.0 * (1.0 - PI / 4.0) - 125.0 * (std::sqrt(2.0) - 1.0);
        double v = volume_of(r);
        std::fprintf(stderr, "  corner setback: volume %.9g in [%.9g, %.9g]\n", v, outside + hull, outside + prisms);
        CHECK(v > outside + hull && v < outside + prisms);
        // The box's six faces, the three rounds and the patch, which SMLib builds as two
        // quadrilaterals per side (2 x 3), each pair split at the middle of its side (its
        // section's or face curve's) and all six meeting at the patch's centre. Vertices: the
        // box's four corners untouched, two at each round's far end, the six section ends, the
        // six side middles and the centre (23); edges by Euler's formula, 23 + 15 - 2.
        int64_t counts[3] = {0, 0, 0};
        CHECK(smlib_brep_counts(r, counts) == SMLIB_OK);
        CHECK(counts[0] == 15 && counts[1] == 36 && counts[2] == 23);
        double lo[3];
        double hi[3];
        CHECK(smlib_brep_bounds(r, lo, hi) == SMLIB_OK);
        // The box's own: the patch stays inside the corner.
        const double want_hi[3] = {60, 40, 20};
        for (int k = 0; k < 3; ++k)
        {
            CHECK(std::fabs(lo[k]) <= 1e-9 && std::fabs(hi[k] - want_hi[k]) <= 1e-9);
        }
        // The patch's extent along each edge: where it meets the faces, the setback sections'
        // ends, each a vertex of the body.
        const double ends[6][3] = {{52, 35, 20}, {52, 40, 15}, {55, 32, 20}, {60, 32, 15}, {55, 40, 12}, {60, 35, 12}};
        SmTArray<SmVertex *> vs;
        r->p->GetVertices(vs);
        for (const auto &e : ends)
        {
            double best = HUGE_VAL;
            for (ULONG k = 0; k < vs.GetSize(); ++k)
            {
                best = std::min(best, (vs[k]->GetPoint() - SmPoint3d(e[0], e[1], e[2])).Length());
            }
            if (best > 1e-6)
            {
                std::fprintf(stderr, "  corner setback: no vertex at (%g, %g, %g), nearest %.3g off\n", e[0], e[1], e[2], best);
            }
            CHECK(best <= 1e-6);
        }
        done(r);
    }
    // One setback per corner, past where the rounds meet: a vertex no point picks, a point at
    // a vertex no blended edge meets, a setback as long as an edge there, a chamfer.
    const double far_corner[3] = {0, 0, 0};
    const double off[3] = {59, 39, 19};
    const double too_long = 20.0;
    CHECK(smlib_blend_setback(part, SMLIB_BLEND_ROUND, edges, radii, 3, off, &setback, 1, 0.1, &r, nullptr, &failed) == SMLIB_ERR_PICK_NONE &&
          failed == 0 && r == nullptr);
    CHECK(smlib_blend_setback(part, SMLIB_BLEND_ROUND, edges, radii, 3, far_corner, &setback, 1, 0.1, &r, nullptr, nullptr) == SMLIB_ERR_ARG &&
          r == nullptr);
    CHECK(smlib_blend_setback(part, SMLIB_BLEND_ROUND, edges, radii, 3, corner, &too_long, 1, 0.1, &r, nullptr, nullptr) == SMLIB_ERR_NO_ROOM &&
          r == nullptr);
    CHECK(smlib_blend_setback(part, SMLIB_BLEND_CHAMFER, edges, radii, 3, corner, &setback, 1, 0.1, &r, nullptr, nullptr) ==
              SMLIB_ERR_UNSUPPORTED &&
          r == nullptr);
    // Two of the corner's three edges rounded: SMLib builds a setback only where all are.
    CHECK(smlib_blend_setback(part, SMLIB_BLEND_ROUND, edges, radii, 2, corner, &setback, 1, 0.1, &r, nullptr, nullptr) ==
              SMLIB_ERR_UNSUPPORTED &&
          r == nullptr);
    smlib_brep_free(part);
}

void shelled_dome_rim()
{
    std::vector<int64_t> kinds = {SMLIB_CURVE_LINE, SMLIB_CURVE_ARC, SMLIB_CURVE_LINE};
    std::vector<double> reals = {0, 0, 0, 20, 0, 0,                                // (0,0,0) -> (20,0,0)
                                 0, 0, 0, 0, -1, 0, 1, 0, 0, 20, 0, PI / 2,      // the arc, about -Y from +X to +Z
                                 0, 0, 20, 0, 0, 0};                             // (0,0,20) -> (0,0,0)
    const double origin[3] = {0, 0, 0}, up[3] = {0, 0, 1};
    smlib_brep *dome = nullptr, *shell = nullptr;
    CHECK(smlib_revolve(kinds.data(), 3, reals.data(), static_cast<int64_t>(reals.size()), origin, up, 2 * PI, 1, &dome) == SMLIB_OK);
    if (dome == nullptr)
    {
        return;
    }
    const double base[3] = {10, 0, 0};
    CHECK(smlib_shell(dome, base, 1, 0.01, -2.0, &shell, nullptr) == SMLIB_OK);
    if (shell != nullptr)
    {
        char buf[4096];
        int64_t len = 0, bad = -1;
        CHECK(smlib_brep_check(shell, 2, buf, sizeof buf, &len, &bad) == SMLIB_OK && bad == 0);
        double volume = 0;
        CHECK(smlib_brep_volume(shell, 1e-8, &volume) == SMLIB_OK);
        CHECK(std::fabs(volume - 2.0 / 3.0 * PI * (8000.0 - 5832.0)) < 1e-5);

        // The rim's inner edge: on the base, at radius 18 to rounding.
        SmTArray<SmEdge *> edges;
        shell->p->GetEdges(edges);
        int rims = 0;
        for (ULONG i = 0; i < edges.GetSize(); i++)
        {
            SmExtent1d ivl = edges[i]->GetInterval();
            SmPoint3d mid;
            CHECK(edges[i]->GetCurve()->EvaluatePoint(ivl.Evaluate(0.5), mid) == SM_SUCCESS);
            if (std::fabs(mid.z) > 1e-6 || std::hypot(mid.x, mid.y) > 19.0)
            {
                continue;
            }
            rims++;
            double worst = 0;
            for (int k = 0; k <= 256; k++)
            {
                SmPoint3d p;
                CHECK(edges[i]->GetCurve()->EvaluatePoint(ivl.Evaluate(k / 256.0), p) == SM_SUCCESS);
                worst = std::max(worst, std::fabs(std::hypot(p.x, p.y) - 18.0) + std::fabs(p.z));
            }
            CHECK(worst < 1e-12);
        }
        CHECK(rims == 1);
    }
    smlib_brep_free(shell);
    smlib_brep_free(dome);
}

// A 20 x 20 x 10 box with a through hole of radius 3 at (10, 10), its top rim chamfered 0.5 (picked
// a quarter turn from the hole's seam): the
// chamfer is the cone between the circles r 3.5 on the top and r 3 half a millimetre down, an
// SmCone (it was a surface of revolution, which SolidCraft read back as a rational B-spline), and
// it takes the frustum less the hole's piece, pi (0.5 / 3) (3.5^2 + 3.5 * 3 + 3^2) - 9 pi 0.5.
void chamfered_hole_rim_is_a_cone()
{
    smlib_brep *b = box(0, 0, 0, 20, 20, 10), *c = cylinder(10, 10, -1, 3, 12), *holed = nullptr, *r = nullptr;
    CHECK(smlib_boolean(b, c, SMLIB_BOOLEAN_SUBTRACT, &holed) == SMLIB_OK);
    const double rim[3] = {10, 13, 10};
    CHECK(holed != nullptr && smlib_chamfer(holed, rim, 1, 0.1, 0.5, &r, nullptr) == SMLIB_OK);
    if (r != nullptr)
    {
        CHECK(sound(r));
        int cones = 0;
        SmTArray<SmFace *> faces;
        r->p->GetFaces(faces);
        for (ULONG i = 0; i < faces.GetSize(); i++)
        {
            const SmCone *cone = SM_CAST_PTR(SmCone, faces[i]->GetSurface());
            if (cone != nullptr && std::fabs(cone->GetBotRadius() - cone->GetTopRadius()) > 0.1)
            {
                cones++;
            }
        }
        CHECK(cones == 1);
        const double removed = PI * 0.5 / 3.0 * (3.5 * 3.5 + 3.5 * 3.0 + 9.0) - 9.0 * PI * 0.5;
        CHECK(near(volume_of(r), 4000.0 - 9.0 * PI * 10.0 - removed, 1e-7));
    }
    for (smlib_brep *x : {b, c, holed, r})
    {
        smlib_brep_free(x);
    }
}

// An 80 x 60 x 30 box, its four upright edges rounded R6 and shelled 2 inward, open at its top
// (SolidCraft's enclosure rim: robustness lip_on_a_rounded_shelled_rim). The inner walls are the
// planes 2 in and the rounds' offsets, cylinders R4 about the same axes: the rim's inner edges
// on the top are those lines and quarter circles, exact conics on the offsets. They used to be
// cubics marched over the offset surfaces, 2.4e-6 off the circles, so a sketch projected from the rim came apart at its corners and the
// lip refused the rim's edges.
void shelled_rounded_box_rim()
{
    smlib_brep *b = box(-40, -30, 0, 80, 60, 30), *round = nullptr, *shell = nullptr;
    const double edges[12] = {40, 30, 15, -40, 30, 15, -40, -30, 15, 40, -30, 15};
    CHECK(smlib_fillet(b, edges, 4, 0.1, 6.0, &round, nullptr) == SMLIB_OK);
    const double top[3] = {0, 0, 30};
    CHECK(round != nullptr && smlib_shell(round, top, 1, 0.1, -2.0, &shell, nullptr) == SMLIB_OK);
    // The worst distance from the top rim's arcs of `body` (on circles of radius r about
    // (+-34, +-24)) to their ideal circles, and how many arcs there are; each must be an SMLib
    // circle.
    auto arcs = [](const smlib_brep *body, double r, int &count)
    {
        double worst = 0;
        count = 0;
        SmTArray<SmEdge *> all;
        body->p->GetEdges(all);
        for (ULONG i = 0; i < all.GetSize(); i++)
        {
            SmExtent1d ivl = all[i]->GetInterval();
            SmPoint3d mid;
            CHECK(all[i]->GetCurve()->EvaluatePoint(ivl.Evaluate(0.5), mid) == SM_SUCCESS);
            const double off = std::fabs(std::hypot(std::fabs(mid.x) - 34.0, std::fabs(mid.y) - 24.0) - r);
            if (std::fabs(mid.z - 30.0) > 1e-6 || std::fabs(mid.x) < 34.0 || std::fabs(mid.y) < 24.0 || off > 1e-3)
            {
                continue;
            }
            count++;
            CHECK(all[i]->GetCurve()->IsKindOf(SmCircle_TYPE));
            for (int k = 0; k <= 256; k++)
            {
                SmPoint3d p;
                CHECK(all[i]->GetCurve()->EvaluatePoint(ivl.Evaluate(k / 256.0), p) == SM_SUCCESS);
                worst = std::max(worst, std::fabs(std::hypot(std::fabs(p.x) - 34.0, std::fabs(p.y) - 24.0) - r) + std::fabs(p.z - 30.0));
            }
        }
        return worst;
    };
    if (shell != nullptr)
    {
        CHECK(sound(shell));
        // The rounds' offsets: R4 to 1e-12. The shell extends them as cylinders; a tangent
        // extension, into which the rim's faces reached a hair, came back as cylinders 3.2e-10
        // off R4 by an amount that depended on the platform. The rim's arcs lie on them.
        double radius = 0.0;
        SmTArray<SmFace *> faces;
        shell->p->GetFaces(faces);
        for (ULONG i = 0; i < faces.GetSize(); i++)
        {
            const SmCone *c = SM_CAST_PTR(SmCone, faces[i]->GetSurface());
            if (c != nullptr && std::fabs(c->GetBotRadius() - 4.0) < 1e-3)
            {
                radius = c->GetBotRadius();
            }
        }
        int n_shell = 0;
        const double w_shell = arcs(shell, radius, n_shell);
        std::printf("shelled rounded box rim: arcs %.3g off radius 4 %+.3g\n", w_shell, radius - 4.0);
        CHECK(n_shell == 4 && std::fabs(radius - 4.0) < 1e-12);
        CHECK(w_shell < 1e-12);
    }
    smlib_brep_free(shell);
    smlib_brep_free(round);
    smlib_brep_free(b);
}

// A 40 x 40 x 20 box with one upright edge rounded R5, as SolidCraft hands it over
// (tests/data/rounded-corner-box.txt: the round an SmCylinder in one rational span), shelled 2
// inward, open at its top. The round's offset is the cylinder R3 about the same axis, which the
// shell extends round its axis as itself (SmOffsetExecutive::GetCreateExtendedSurface): the face
// on that extension comes back an SmCone, R3 to 1e-12. Both sides extended in one call left a
// single knot at the far end, which SmCone::IsNurbSurfaceCone does not take for a circle, and
// the face stayed a surface of revolution.
// That wall, in the shell as SolidCraft hands it back (rounded-corner-box-shelled.txt), moved 0.5
// into the cavity (R3 to R2.5) steps in from the two planes it was tangent to. The step walls
// are ruled between those planes' edges and the moved wall's: one of these lines came off the
// extended cylinder with its weight there, 0.99998828, and ruled against a line of weight 1 the
// wall was rational across, so the plane it became put its trims 5.9e-6 off, past the 5e-6
// tolerance. The result now passes SMLib's check, and the cavity loses the quarter ring between
// the radii, pi/4 (3^2 - 2.5^2) 18.
void shelled_single_round_stays_a_cylinder()
{
    auto read = [](const char *path)
    {
        smlib_flat *f = nullptr;
        smlib_brep *b = nullptr;
        CHECK(smlib_flat_read(path, &f) == SMLIB_OK && smlib_brep_import(f, &b) == SMLIB_OK);
        smlib_flat_free(f);
        return b;
    };
    smlib_brep *round = read(SMLIB_TEST_DATA "/rounded-corner-box.txt"), *shell = nullptr;
    const double top[3] = {20, 20, 20};
    CHECK(round != nullptr && smlib_shell(round, top, 1, 0.1, -2.0, &shell, nullptr) == SMLIB_OK);
    if (shell != nullptr)
    {
        CHECK(sound(shell));
        int inner = 0;
        double radius = 0.0;
        SmTArray<SmFace *> faces;
        shell->p->GetFaces(faces);
        for (ULONG i = 0; i < faces.GetSize(); i++)
        {
            const SmCone *c = SM_CAST_PTR(SmCone, faces[i]->GetSurface());
            if (c != nullptr && std::fabs(c->GetBotRadius() - 3.0) < 1e-3)
            {
                inner++;
                radius = c->GetBotRadius();
            }
        }
        std::printf("shelled single round: %d inner round(s), radius 3 %+.3g\n", inner, radius - 3.0);
        CHECK(inner == 1 && std::fabs(radius - 3.0) < 1e-12);
    }
    smlib_brep *handed = read(SMLIB_TEST_DATA "/rounded-corner-box-shelled.txt"), *moved = nullptr;
    const double wall[3] = {35.0 + 3.0 / std::sqrt(2.0), 35.0 + 3.0 / std::sqrt(2.0), 10.0};
    CHECK(handed != nullptr && sound(handed) && smlib_offset_faces(handed, wall, 1, 0.1, 0.5, &moved, nullptr) == SMLIB_OK);
    if (moved != nullptr)
    {
        CHECK(sound(moved));
        CHECK(near(volume_of(moved) - volume_of(handed), PI / 4.0 * (9.0 - 6.25) * 18.0, 1e-6));
    }
    for (smlib_brep *b : {round, shell, handed, moved})
    {
        smlib_brep_free(b);
    }
}

// A 60 x 40 x 15 box with a blind drilled hole (diameter 5, 8 deep from the top, its tip a 118
// degree cone) at (30, 20), shelled 2 inward, open at the side y = 0 and then at the top
// (SolidCraft robustness shell_after_a_drilled_hole). The cone's offset ends on the circle its
// apex's normals sweep out (r 2 cos 59 deg); the apex itself, a vertex with only the seam, had no
// offset and the cavity was open there, so the shell failed (1001). Its offset is the sphere of
// radius 2 about the tip.
void shelled_box_with_a_drilled_hole()
{
    const double tip = 2.5 / std::tan(59.0 * PI / 180.0);
    std::vector<int64_t> kinds = {SMLIB_CURVE_LINE, SMLIB_CURVE_LINE, SMLIB_CURVE_LINE, SMLIB_CURVE_LINE};
    std::vector<double> reals = {30, 20, 16, 32.5, 20, 16, 32.5, 20, 16, 32.5, 20, 7,
                                 32.5, 20, 7, 30, 20, 7 - tip, 30, 20, 7 - tip, 30, 20, 16};
    const double origin[3] = {30, 20, 0}, up[3] = {0, 0, 1};
    smlib_brep *b = box(0, 0, 0, 60, 40, 15), *drill = nullptr, *holed = nullptr;
    CHECK(smlib_revolve(kinds.data(), 4, reals.data(), static_cast<int64_t>(reals.size()), origin, up, 2 * PI, 1, &drill) == SMLIB_OK);
    CHECK(drill != nullptr && smlib_boolean(b, drill, SMLIB_BOOLEAN_SUBTRACT, &holed) == SMLIB_OK);
    if (holed != nullptr)
    {
        const double opens[2][3] = {{30, 0, 7.5}, {10, 10, 15}};
        for (const double *open : opens)
        {
            smlib_brep *r = nullptr;
            CHECK(smlib_shell(holed, open, 1, 0.1, -2.0, &r, nullptr) == SMLIB_OK);
            CHECK(r != nullptr && sound(r));
            // The tip's offset: the sphere of radius 2 about it.
            int caps = 0;
            if (r != nullptr)
            {
                SmTArray<SmFace *> faces;
                r->p->GetFaces(faces);
                for (ULONG i = 0; i < faces.GetSize(); i++)
                {
                    const SmSphere *sphere = SM_CAST_PTR(SmSphere, faces[i]->GetSurface());
                    caps += sphere != nullptr && std::fabs(sphere->GetRadius() - 2.0) < 1e-9 ? 1 : 0;
                }
            }
            CHECK(caps == 1);
            smlib_brep_free(r);
        }
    }
    smlib_brep_free(holed);
    smlib_brep_free(drill);
    smlib_brep_free(b);
}

// A 30 x 35 x 7.5 box notched 22.5 x 5 down to z = 2.5 at its far corner, shelled 1.5 inward
// open at the far side's strip under the notch (SolidCraft robustness
// shell_thicker_than_a_step_says_so): the wall would move the strip's top corner past the
// end of its edge, turning the strip inside out; that is SMLIB_ERR_UNSUPPORTED. Open at the
// notch's floor it shells.
void shell_thicker_than_a_step()
{
    smlib_brep *b = box(0, 0, 0, 30, 35, 7.5), *notch = box(0, 30, 2.5, 22.5, 20, 22.5), *stepped = nullptr;
    CHECK(smlib_boolean(b, notch, SMLIB_BOOLEAN_SUBTRACT, &stepped) == SMLIB_OK);
    if (stepped != nullptr)
    {
        const double strip[3] = {18.75, 35, 1.25}, floor[3] = {18.75, 32.5, 2.5};
        smlib_brep *r = nullptr;
        CHECK(smlib_shell(stepped, strip, 1, 0.1, -1.5, &r, nullptr) == SMLIB_ERR_UNSUPPORTED && r == nullptr);
        CHECK(smlib_shell(stepped, floor, 1, 0.1, -1.5, &r, nullptr) == SMLIB_OK && r != nullptr && sound(r));
        smlib_brep_free(r);
    }
    for (smlib_brep *x : {b, notch, stepped})
    {
        smlib_brep_free(x);
    }
}

// A stepped block, 40 x 30 x 10 with its half x < 20 raised to 20 (SolidCraft kernel case
// shell/cov-replace-face-step-onto-the-higher-face): the low top moved up onto the high top's
// plane closes the step wall between them to nothing (it was SMLIB_ERR_UNSUPPORTED: the wall's
// edges along the sides closed to points). The wall goes and the tops, meeting on one plane, are
// one face: the 40 x 30 x 20 box, 6 faces. The high top moved down onto the low one the same
// way: the 40 x 30 x 10 box.
void offset_faces_close_a_step_wall()
{
    smlib_brep *base = box(0, 0, 0, 40, 30, 10), *tower = box(0, 0, 0, 20, 30, 20), *stepped = nullptr;
    CHECK(smlib_boolean(base, tower, SMLIB_BOOLEAN_UNION, &stepped) == SMLIB_OK);
    CHECK(stepped != nullptr && stepped->p->GetNumFaces() == 8);
    if (stepped != nullptr)
    {
        struct Move
        {
            double at[3], distance, volume, area, top;
        };
        const Move moves[] = {{{30, 15, 10}, 10.0, 24000.0, 5200.0, 20.0}, {{10, 15, 20}, -10.0, 12000.0, 3800.0, 10.0}};
        for (const Move &m : moves)
        {
            smlib_brep *r = nullptr;
            CHECK(smlib_offset_faces(stepped, m.at, 1, 0.1, m.distance, &r, nullptr) == SMLIB_OK);
            if (r != nullptr)
            {
                double area = -1.0, lo[3] = {}, hi[3] = {};
                CHECK(near(volume_of(r), m.volume, 1e-9));
                CHECK(smlib_brep_area(r, 1e-9, &area) == SMLIB_OK && near(area, m.area, 1e-9));
                CHECK(smlib_brep_bounds(r, lo, hi) == SMLIB_OK && std::fabs(hi[2] - m.top) < 1e-9 && std::fabs(hi[0] - 40.0) < 1e-9);
                CHECK(r->p->GetNumFaces() == 6);
                done(r);
            }
        }
    }
    for (smlib_brep *x : {base, tower, stepped})
    {
        smlib_brep_free(x);
    }
}

// A pin r 10, 20 tall, its side drafted 5 degrees about its base circle (neutral z = 0, pull +z):
// the frustum from r 10 to 10 - 20 tan 5. A 30 x 30 x 10 block's through hole r 5 drafted the
// same way widens to 5 + 10 tan 5 at the top. (SolidCraft draft_on_a_cylinder_side_makes_a_frustum,
// draft_takes_a_cylinders_side_by_name.) A cone r 10 at its base to 8 at its top, 20 tall, drafted
// the same way turns its generators 5 degrees more: r 10 + 20 tan(atan(-0.1) - 5 deg) at the top.
void drafted_cylinders_are_cones()
{
    const double neutral[3] = {0, 0, 0}, up[3] = {0, 0, 1}, a = 5.0 * PI / 180.0;
    auto frustum = [](double h, double r0, double r1) { return PI * h / 3.0 * (r0 * r0 + r0 * r1 + r1 * r1); };
    smlib_brep *pin = cylinder(0, 0, 0, 10, 20), *r = nullptr;
    const double side[3] = {0, 10, 10};
    CHECK(smlib_draft(pin, side, 1, 0.1, neutral, up, up, a, &r, nullptr) == SMLIB_OK);
    if (r != nullptr)
    {
        CHECK(sound(r));
        CHECK(near(volume_of(r), frustum(20, 10, 10 - 20 * std::tan(a)), 1e-7));
    }
    smlib_brep_free(r);
    r = nullptr;
    smlib_brep *b = box(-15, -15, 0, 30, 30, 10), *c = cylinder(0, 0, -1, 5, 12), *holed = nullptr;
    CHECK(smlib_boolean(b, c, SMLIB_BOOLEAN_SUBTRACT, &holed) == SMLIB_OK);
    const double wall[3] = {0, 5, 5};
    CHECK(holed != nullptr && smlib_draft(holed, wall, 1, 0.1, neutral, up, up, a, &r, nullptr) == SMLIB_OK);
    if (r != nullptr)
    {
        CHECK(sound(r));
        CHECK(near(volume_of(r), 9000.0 - frustum(10, 5, 5 + 10 * std::tan(a)), 1e-7));
    }
    smlib_brep_free(r);
    r = nullptr;
    const double base[3] = {0, 0, 0}, flank[3] = {0, 9, 10};
    smlib_brep *cone = nullptr;
    CHECK(smlib_make_cone(base, 10, 8, 20, &cone) == SMLIB_OK);
    CHECK(cone != nullptr && smlib_draft(cone, flank, 1, 0.1, neutral, up, up, a, &r, nullptr) == SMLIB_OK);
    if (r != nullptr)
    {
        CHECK(sound(r));
        const double top = 10 + 20 * std::tan(std::atan(-0.1) - a);
        std::printf("drafted cone: volume %.12g vs %.12g\n", volume_of(r), frustum(20, 10, top));
        CHECK(near(volume_of(r), frustum(20, 10, top), 1e-7));
    }
    for (smlib_brep *x : {pin, b, c, holed, cone, r})
    {
        smlib_brep_free(x);
    }
}

// The analytic STEPInversion() answers: SM_SUCCESS on the surface, SM_ERR_OUTSIDE_OF_DOMAIN
// outside a cylinder's sweep, SM_ERR_NOT_WITHIN_TOLERANCE off the surface (a sphere's with the
// dropped point's parameters), each with its location; SMLib's callers pass any other status on.
void analytic_step_inversion()
{
    const SmContext &ctx = *SmApiGetOrCreateContext();
    const double tol = 1e-6;
    SmPoint2d uv;
    SmLocationType loc = SM_LT_INTERIOR;
    SmCylinder *cyl = new (ctx) SmCylinder(SmPoint3d(0, 0, 0), SmVector3d(1, 0, 0), SmVector3d(0, 1, 0), 1.0, 0.0, 90.0, 2.0, FALSE, FALSE,
                                           FALSE, &ctx);
    const SmExtent2d cd = cyl->GetSTEPUVDomain();
    const double h = std::sqrt(0.5);
    CHECK(cyl->STEPInversion(cd, SmPoint3d(h, h, 1), tol, uv, loc) == SM_SUCCESS && loc == SM_LT_INTERIOR);
    CHECK(std::fabs(uv.x - 45.0) < 1e-9 && std::fabs(uv.y - 1.0) < 1e-9);
    CHECK(cyl->STEPInversion(cd, SmPoint3d(-1, 0, 1), tol, uv, loc) == SM_ERR_OUTSIDE_OF_DOMAIN && loc == SM_LT_EXTERIOR);
    CHECK(std::fabs(uv.y - 1.0) < 1e-9);
    CHECK(cyl->STEPInversion(cd, SmPoint3d(2, 0.1, 1), tol, uv, loc) == SM_ERR_NOT_WITHIN_TOLERANCE && loc == SM_LT_EXTERIOR);
    delete cyl;

    SmSphere *sphere = nullptr;
    CHECK(SmSphere::CreateCanonical(ctx, SmAxis2Placement(SmPoint3d(0, 0, 0), SmVector3d(1, 0, 0), SmVector3d(0, 1, 0)), 2.0, sphere) ==
          SM_SUCCESS);
    if (sphere != nullptr)
    {
        const SmExtent2d sd = sphere->GetSTEPUVDomain();
        CHECK(sphere->STEPInversion(sd, SmPoint3d(0, 2, 0), tol, uv, loc) == SM_SUCCESS && loc == SM_LT_INTERIOR);
        CHECK(std::fabs(uv.x - 90.0) < 1e-9 && std::fabs(uv.y) < 1e-9);
        CHECK(sphere->STEPInversion(sd, SmPoint3d(0, 3, 0), tol, uv, loc) == SM_ERR_NOT_WITHIN_TOLERANCE && loc == SM_LT_INTERIOR);
        CHECK(std::fabs(uv.x - 90.0) < 1e-9 && std::fabs(uv.y) < 1e-9);
        CHECK(sphere->STEPInversion(sd, SmPoint3d(0, 0, 0), tol, uv, loc) == SM_ERR_NOT_WITHIN_TOLERANCE && loc == SM_LT_EXTERIOR);
        delete sphere;
    }
    CHECK(SmSurface::STEPInversionAnswered(SM_ERR_NOT_WITHIN_TOLERANCE) && SmSurface::STEPInversionAnswered(SM_ERR_OUTSIDE_OF_DOMAIN));
    CHECK(!SmSurface::STEPInversionAnswered(SM_ERR) && !SmSurface::STEPInversionAnswered(SM_ERR_NOT_CONVERGING));
}

// A cylinder whose offset construction gives no surface although the offset does not collapse
// (the cylinder's radius is larger than the distance): an offset that is missing for any other
// reason than a collapse.
class NoOffsetCylinder : public SmCylinder
{
  public:
    explicit NoOffsetCylinder(const SmCylinder &c) : SmCylinder(c)
    {
    }
    SmStatus CreateOffsetSurface(const SmContext &, double, SmApproxTol3d, SmSurface *&rOffset) const override
    {
        rOffset = nullptr;
        return SM_SUCCESS;
    }
};

// The largest gap, over every edgeuse's UV trim curve, between the curve lifted onto its face's
// surface and its edge's curve at the same parameter, as a multiple of the edge's tolerance
// (sampled at 33 points per edgeuse).
double worst_pcurve_gap(const smlib_brep *b)
{
    SmTArray<SmEdgeuse *> uses;
    b->p->GetEdgeuses(uses);
    double worst = 0.0;
    for (ULONG i = 0; i < uses.GetSize(); i++)
    {
        SmEdgeuse *use = uses[i];
        SmBSplineCurve *uv = use->GetUVTrimCurvePointer();
        if (uv == nullptr || use->GetFaceuse() == nullptr)
        {
            continue;
        }
        SmEdge *e = use->GetEdge();
        SmSurface *s = use->GetFaceuse()->GetFace()->GetSurface();
        SmExtent1d ivl = e->GetInterval();
        for (int k = 0; k <= 32; k++)
        {
            double t = ivl.GetMin() + ivl.GetLength() * k / 32.0;
            SmPoint3d c, q, w;
            CHECK(e->GetCurve()->EvaluatePoint(t, c) == SM_SUCCESS && uv->EvaluatePoint(t, w) == SM_SUCCESS &&
                  s->EvaluatePoint(SmPoint2d(w.x, w.y), q) == SM_SUCCESS);
            worst = std::max(worst, c.DistanceBetween(q) / e->GetTolerance());
        }
    }
    return worst;
}

// A planar profile (lines and arcs about +z, in the plane z0) extruded `h` up.
struct Profile
{
    std::vector<int64_t> kinds;
    std::vector<double> reals;
    Profile &line(double ax, double ay, double bx, double by, double z0)
    {
        kinds.push_back(SMLIB_CURVE_LINE);
        reals.insert(reals.end(), {ax, ay, z0, bx, by, z0});
        return *this;
    }
    // About +z (nz 1) counter-clockwise from a0 to a1, or about -z (nz -1).
    Profile &arc(double cx, double cy, double z0, double nz, double r, double a0, double a1)
    {
        kinds.push_back(SMLIB_CURVE_ARC);
        reals.insert(reals.end(), {cx, cy, z0, 0, 0, nz, 1, 0, 0, r, a0, a1});
        return *this;
    }
    smlib_brep *extrude(double h) const
    {
        const double up[3] = {0, 0, 1};
        smlib_brep *b = nullptr;
        CHECK(smlib_extrude(kinds.data(), static_cast<int64_t>(kinds.size()), reals.data(), static_cast<int64_t>(reals.size()), up, h, 0.0,
                            SMLIB_DRAFT_CORNER_CURVATURE, 1, &b, nullptr) == SMLIB_OK);
        return b;
    }
};

// A 50 x 40 x 10 block with its corners rounded r 4.5, and a lip on its top: a 1 x 1 square swept
// along half the block's rim, flush with it, round two corners (SolidCraft's lip swept flush round
// a corner). The sweep's walls are B-splines approximating the swept square; its outer bottom edge
// lies on the block's rim arc within tolerance and crosses it twice, 0.05 apart, at the arc's
// middle. The union split the rim arc there, and the edge where the two meet along the whole arc
// was re-pointed to the 0.05 piece in the middle, whose vertices are not its ends: the imprint
// found no vertex there (SM_ERR_NULL_POINTER). The lip sits on the block: the union is the sum,
// the block (50 x 40 - (4 - pi) 4.5^2) 10 and the lip 41 + 31 + 41 + 2 pi (4.5^2 - 3.5^2) / 4.
void lip_swept_flush_round_a_corner()
{
    Profile base;
    base.line(20.5, -20, -20.5, -20, 0).arc(-20.5, -15.5, 0, 1, 4.5, PI, 1.5 * PI).line(-25, -15.5, -25, 15.5, 0);
    base.arc(-20.5, 15.5, 0, 1, 4.5, 0.5 * PI, PI).line(-20.5, 20, 20.5, 20, 0).arc(20.5, 15.5, 0, 1, 4.5, 0, 0.5 * PI);
    base.line(25, 15.5, 25, -15.5, 0).arc(20.5, -15.5, 0, 1, 4.5, 1.5 * PI, 2 * PI);
    smlib_brep *block = base.extrude(10.0), *lip = nullptr, *r = nullptr;
    // The square in the plane x = 20.5, from (y, z) = (-20, 10) to (-19, 11).
    const std::vector<int64_t> sk = {SMLIB_CURVE_LINE, SMLIB_CURVE_LINE, SMLIB_CURVE_LINE, SMLIB_CURVE_LINE};
    const std::vector<double> sr = {20.5, -20, 10, 20.5, -19, 10, 20.5, -19, 10, 20.5, -19, 11,
                                    20.5, -19, 11, 20.5, -20, 11, 20.5, -20, 11, 20.5, -20, 10};
    Profile path;
    path.line(20.5, -20, -20.5, -20, 10).arc(-20.5, -15.5, 10, -1, 4.5, 0.5 * PI, PI).line(-25, -15.5, -25, 15.5, 10);
    path.arc(-20.5, 15.5, 10, -1, 4.5, PI, 1.5 * PI).line(-20.5, 20, 20.5, 20, 10);
    CHECK(smlib_sweep(sk.data(), static_cast<int64_t>(sk.size()), sr.data(), static_cast<int64_t>(sr.size()), path.kinds.data(),
                      static_cast<int64_t>(path.kinds.size()), path.reals.data(), static_cast<int64_t>(path.reals.size()), SMLIB_SWEEP_FOLLOW, 1,
                      &lip) == SMLIB_OK);
    CHECK(block != nullptr && lip != nullptr && smlib_boolean(block, lip, SMLIB_BOOLEAN_UNION, &r) == SMLIB_OK);
    if (r != nullptr)
    {
        CHECK(sound(r));
        const double blk = (2000.0 - (4.0 - PI) * 4.5 * 4.5) * 10.0, ring = 41.0 + 31.0 + 41.0 + 2.0 * PI * (4.5 * 4.5 - 3.5 * 3.5) / 4.0;
        std::printf("lip swept flush round a corner: %.9f, the block and the lip %.9f\n", volume_of(r), blk + ring);
        CHECK(near(volume_of(r), blk + ring, 1e-6));
    }
    for (smlib_brep *x : {block, lip, r})
    {
        smlib_brep_free(x);
    }
}

// A 41 1/6 x 21 plate and, beside it, the corner a circle of radius 21 about (TIP, -21) leaves
// against the plate's back line y = 0, extruded 10 (SolidCraft's knife-edge part, #851): the arc
// wall meets the back face tangentially at x = TIP, a knife edge (a cusp, the faces folded back
// on each other). SmEdge::GetBlendEdgeuse took the edge for a smooth one (faces tangent), the
// solver could not tell which way to offset its faces, and the corners at its ends were taken
// for corners where the faces run smoothly into each other (the fillet shrinking to a point).
// Rounded r, the ball 2 r wide touches the back face sqrt(4 R r) back from the tip: it takes
// off the quadrilateral tip - wall centre - ball centre - back contact, less the wall's sector
// and the ball's, all the way up. Mirrored, the edge and its faces' uses reverse: the same.
void knife_edge_rounds()
{
    const double R = 21.0, TIP = 62.0 + 1.0 / 6.0, H = 10.0;
    Profile p;
    p.line(0, 0, 0, -R, 0).line(0, -R, TIP - R, -R, 0).arc(TIP, -R, 0, -1, R, PI, 1.5 * PI).line(TIP, 0, 0, 0, 0);
    smlib_brep *part = p.extrude(H);
    const double before = (TIP - R) * R * H + R * R * (1.0 - PI / 4.0) * H;
    for (double r : {0.5, 2.0, 5.0})
    {
        const double beta = std::acos((R - r) / (R + r)), ox = TIP - std::sqrt(4.0 * R * r);
        const double quad = 0.5 * (TIP - ox) * (R + r);
        const double off = quad - 0.5 * R * R * beta - 0.5 * r * r * (PI - beta);
        for (int mirrored = 0; mirrored < 2; ++mirrored)
        {
            const double m[4][4] = {{mirrored ? -1.0 : 1.0, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};
            smlib_brep *b = nullptr, *round = nullptr;
            CHECK(part != nullptr && smlib_brep_transform(part, m, &b) == SMLIB_OK);
            const double tip[3] = {mirrored ? -TIP : TIP, 0, H / 2};
            CHECK(b != nullptr && smlib_fillet(b, tip, 1, 0.01, r, &round, nullptr) == SMLIB_OK);
            if (round != nullptr)
            {
                int64_t counts[3] = {0, 0, 0};
                CHECK(smlib_brep_counts(round, counts) == SMLIB_OK && counts[0] == 7);
                CHECK(sound(round));
                std::printf("knife edge rounded r%g%s: %.9f, expected %.9f\n", r, mirrored ? " (mirrored)" : "", volume_of(round), before - off * H);
                CHECK(near(volume_of(round), before - off * H, 1e-9));
            }
            smlib_brep_free(round);
            smlib_brep_free(b);
        }
    }
    smlib_brep_free(part);
}

// An overflow goes first. R21 on the 20 cube's edge along x at y = z = 0 overflows its faces (its
// circle about (y, z) = (21, 21) meets them only at 21, past the cube) and takes all of both:
// what is left is the square [0, 20]^2 inside that circle, with u = 21 - y,
// int_1^sqrt(440) (sqrt(441 - u^2) - 1) du, along 20: 6107.529370. Rounded with R5 on the
// corner's other two edges, which lie on those faces, the R5 edges are gone (dropped, flagged) and
// the body is R21's alone. With R5 on the back upright x = 0, y = 20 instead, the overflow cuts
// that edge off at z = 21 - sqrt(440) and R5 rounds the rest, ending on R21's cylinder: it takes
// the corner of the section outside the R5 circle about (5, 15), from the cylinder
// z_c(y) = 21 - sqrt(441 - (y - 21)^2) up to 20; with y = 15 + 5 sin t,
// int_0^pi/2 5 (1 - cos t) (20 - z_c(y)) 5 cos t dt = 106.612300588. Its end, where the R5 and
// R21 cylinders cross (a quartic, no exact B-spline), is held to the edge's tolerance tau: the
// measure integrates over the faces it trims, so the volume is good to 2 tau L rho / 3 there, L
// <= 5 pi / 2 + 1 the edge's length and rho <= 21 its distance from the origin.
void overflow_goes_first()
{
    auto exact_volume = [](const smlib_brep *b)
    {
        double m[11] = {};
        CHECK(smlib_brep_mass_properties(b, 1e-8, m) == SMLIB_OK);
        return m[1];
    };
    auto F = [](double u) { return 0.5 * (u * std::sqrt(441.0 - u * u) + 441.0 * std::asin(u / 21.0)); };
    const double w = std::sqrt(440.0), alone = 20.0 * (F(w) - F(1.0) - (w - 1.0));
    CHECK(near(alone, 6107.529370, 1e-10));
    // The R5 corner over the cylinder, by Simpson's rule in t (smooth there) on 4000 spans.
    auto strip = [](double t)
    {
        const double y = 15.0 + 5.0 * std::sin(t), zc = 21.0 - std::sqrt(441.0 - (y - 21.0) * (y - 21.0));
        return 5.0 * (1.0 - std::cos(t)) * (20.0 - zc) * 5.0 * std::cos(t);
    };
    const int n = 4000;
    const double h = 0.5 * PI / n;
    double sum = strip(0.0) + strip(0.5 * PI);
    for (int i = 1; i < n; ++i)
    {
        sum += (i % 2 ? 4.0 : 2.0) * strip(i * h);
    }
    const double corner = sum * h / 3.0;
    CHECK(near(corner, 106.612300588, 1e-9));
    smlib_brep *cube = box(0, 0, 0, 20, 20, 20);
    struct Case
    {
        const char *name;
        std::vector<double> points, radii;
        std::vector<int64_t> dropped;
        double volume;
        bool quartic_end;
    };
    const Case cases[] = {
        {"R21 with R5 on the corner's other edges", {10, 0, 0, 0, 10, 0, 0, 0, 10}, {21, 5, 5}, {0, 1, 1}, alone, false},
        {"R21 with R5 on the back upright", {10, 0, 0, 0, 20, 10}, {21, 5}, {0, 0}, alone - corner, true},
    };
    for (const Case &c : cases)
    {
        const int64_t k = static_cast<int64_t>(c.radii.size());
        std::vector<int64_t> dropped(c.radii.size(), -1);
        smlib_brep *r = nullptr;
        CHECK(smlib_blend_dropping(cube, SMLIB_BLEND_ROUND, c.points.data(), c.radii.data(), k, 1.0, &r, nullptr, dropped.data()) == SMLIB_OK);
        CHECK(dropped == c.dropped);
        if (r != nullptr)
        {
            CHECK(sound(r));
            double tau = 0.0;
            smlib_flat *f = nullptr;
            CHECK(smlib_brep_export(r, &f) == SMLIB_OK);
            const double *t = nullptr;
            int64_t nt = 0;
            CHECK(f != nullptr && smlib_flat_get_f64(f, SMLIB_FLAT_EDGE_TOLERANCE, &t, &nt) == SMLIB_OK);
            for (int64_t i = 0; i < nt; ++i)
            {
                tau = std::max(tau, t[i]);
            }
            smlib_flat_free(f);
            const double bound = c.quartic_end ? 2.0 * tau * (2.5 * PI + 1.0) * 21.0 / 3.0 : 1e-9 * c.volume;
            std::printf("overflow first, %s: %.9f, expected %.9f (within %.3g)\n", c.name, exact_volume(r), c.volume, bound);
            CHECK(std::fabs(exact_volume(r) - c.volume) <= bound);
        }
        smlib_brep_free(r);
    }
    smlib_brep_free(cube);
}

// The knife edge of knife_edge_rounds chamfered d: at a cusp no ball lies d back along both faces
// (their angle is 0), and the chamfer is the chord between the contacts of the ball of radius d
// (the R d round's): on the flat face sqrt(4 R d) back from the tip, on the wall where the line
// from its axis to the ball's centre meets it. It takes the triangle flat contact, tip, wall
// contact less the wall's circular segment R^2/2 (b - sin b), b = acos((R - d)/(R + d)): at
// d = 2, 3.822 wide, 159.23645 removed (SolidCraft rig fillets/851-chamfer-knife-edge-tip-*).
// At d 0.5 the chord meets the flat face at 81 degrees, which the trim took for no side at all.
void knife_edge_chamfers()
{
    const double R = 21.0, TIP = 62.0 + 1.0 / 6.0, H = 10.0;
    Profile p;
    p.line(0, 0, 0, -R, 0).line(0, -R, TIP - R, -R, 0).arc(TIP, -R, 0, -1, R, PI, 1.5 * PI).line(TIP, 0, 0, 0, 0);
    smlib_brep *part = p.extrude(H);
    const double before = (TIP - R) * R * H + R * R * (1.0 - PI / 4.0) * H;
    for (double d : {0.5, 0.8, 1.0, 2.0, 5.0})
    {
        const double flat = TIP - std::sqrt(4.0 * R * d);
        const double wx = TIP - std::sqrt(4.0 * R * d) * R / (R + d), wy = -R + (R - d) * R / (R + d);
        const double triangle = 0.5 * std::fabs((TIP - flat) * wy);
        const double b = std::acos((R - d) / (R + d));
        const double removed = (triangle - 0.5 * R * R * (b - std::sin(b))) * H;
        if (d == 2.0)
        {
            CHECK(near(removed, 159.23645074290118, 1e-12) && near(std::hypot(wx - flat, wy), 3.8221323436236347, 1e-12));
        }
        for (int mirrored = 0; mirrored < 2; ++mirrored)
        {
            const double m[4][4] = {{mirrored ? -1.0 : 1.0, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};
            smlib_brep *body = nullptr, *cut = nullptr;
            CHECK(part != nullptr && smlib_brep_transform(part, m, &body) == SMLIB_OK);
            const double tip[3] = {mirrored ? -TIP : TIP, 0, H / 2};
            CHECK(body != nullptr && smlib_chamfer(body, tip, 1, 0.01, d, &cut, nullptr) == SMLIB_OK);
            if (cut != nullptr)
            {
                CHECK(sound(cut));
                std::printf("knife edge chamfered %g%s: %.9f, expected %.9f\n", d, mirrored ? " (mirrored)" : "", volume_of(cut), before - removed);
                CHECK(near(volume_of(cut), before - removed, 1e-9));
            }
            smlib_brep_free(cut);
            smlib_brep_free(body);
        }
    }
    smlib_brep_free(part);
}

// A block 40 x 20, 20 deep, whose top rises to a crease at x = 20, each half sloping `slope`
// degrees, so the faces turn by twice that across it. Faces within half a degree of each other are
// tangent, a smooth edge with nothing to round (SolidCraft's angle): the 0.3 degree crease is
// refused. The 1 degree crease was refused too, as within 2 degrees; it is an edge the model means
// and rounds: R2 and R20 each lie r (sec(d/2) - 1) inside the crease at its top, d its turn, and
// remove 20 r^2 (tan(d/2) - d/2).
void crease_rounds_from_half_a_degree()
{
    for (double turn : {0.3, 1.0, 3.0})
    {
        const double d = turn * PI / 180.0, h = 20.0 * std::tan(0.5 * d);
        Profile p;
        p.line(0, 0, 40, 0, 0).line(40, 0, 40, 20, 0).line(40, 20, 20, 20 + h, 0).line(20, 20 + h, 0, 20, 0).line(0, 20, 0, 0, 0);
        smlib_brep *block = p.extrude(20.0);
        const double crease[3] = {20, 20 + h, 10};
        auto exact_volume = [](const smlib_brep *b)
        {
            double m[11] = {};
            CHECK(smlib_brep_mass_properties(b, 1e-8, m) == SMLIB_OK);
            return m[1];
        };
        auto top = [](const smlib_brep *b)
        {
            double lo[3] = {}, hi[3] = {};
            CHECK(smlib_brep_bounds(b, lo, hi) == SMLIB_OK);
            return hi[1];
        };
        const double before = block != nullptr ? exact_volume(block) : 0.0;
        for (double r : {2.0, 20.0})
        {
            smlib_brep *round = nullptr;
            const int32_t s = block != nullptr ? smlib_fillet(block, crease, 1, 0.01, r, &round, nullptr) : SMLIB_ERR_EMPTY;
            if (turn < 0.5)
            {
                CHECK(s != SMLIB_OK && round == nullptr);
                continue;
            }
            CHECK(s == SMLIB_OK && round != nullptr && sound(round));
            if (round != nullptr)
            {
                const double removed = 20.0 * r * r * (std::tan(0.5 * d) - 0.5 * d), got = before - exact_volume(round);
                const double sag = r * (1.0 / std::cos(0.5 * d) - 1.0), lowered = top(block) - top(round);
                std::printf("crease %g deg rounded r%g: removed %.9g, expected %.9g; top lowered %.9g, expected %.9g\n", turn, r, got, removed,
                            lowered, sag);
                CHECK(std::fabs(lowered - sag) < 1e-9);
                CHECK(std::fabs(got - removed) < 1e-6);
            }
            smlib_brep_free(round);
        }
        smlib_brep_free(block);
    }
}

// The 60 x 40 x 20 box, its top edge at y = 40 rounded r5, its top edge at x = 60 rounded r3,
// then the vertical edge at x = 60, y = 40 rounded `radius` with the edges tangent to it up into
// the corner (main has the case), or null.
smlib_brep *rounded_corner(double radius)
{
    smlib_brep *slab = box(0, 0, 0, 60, 40, 20);
    smlib_brep *top = nullptr;
    smlib_brep *both = nullptr;
    smlib_brep *corner = nullptr;
    const double top_edge[3] = {30, 40, 20};
    const double end_edge[3] = {60, 20, 20};
    const double up[9] = {60, 40, 7.5, 60, 35.0 + std::sqrt(24.0), 16, 58.931891825700852, 37.55984948554854, 19.294016994391871};
    CHECK(smlib_fillet(slab, top_edge, 1, 0.1, 5.0, &top, nullptr) == SMLIB_OK);
    CHECK(top != nullptr && smlib_fillet(top, end_edge, 1, 0.1, 3.0, &both, nullptr) == SMLIB_OK);
    CHECK(both != nullptr && smlib_fillet(both, up, 3, 0.1, radius, &corner, nullptr) == SMLIB_OK);
    for (smlib_brep *x : {both, top, slab})
    {
        smlib_brep_free(x);
    }
    return corner;
}

// The constant-distance blend (SMLIB_BLEND_CHORD) of tangency-labs/solidcore#453's back rim, on
// the 20 cube bored r5 along y with the bore's wall on its bottom face: refused (main has the
// case), after its solver has computed points across the rim.
void failed_chord_blend()
{
    smlib_brep *cube = box(0, 0, 0, 20, 20, 20);
    SmBrep *bore = nullptr;
    CHECK(SmApiCreateCylinder(SmVector3d(0, 0, 0), 5.0, 20.0, bore) == SM_SUCCESS && bore != nullptr);
    CHECK(SmApiRotate(bore, SmVector3d(0, 0, 0), SmVector3d(1, 0, 0), 90.0) == SM_SUCCESS);
    CHECK(SmApiTranslate(bore, SmVector3d(10, 20, 5)) == SM_SUCCESS);
    smlib_brep tool{bore};
    smlib_brep *holed = nullptr;
    CHECK(smlib_boolean(cube, &tool, SMLIB_BOOLEAN_SUBTRACT, &holed) == SMLIB_OK);
    const double back_rim[3] = {10, 20, 10};
    const double size = 3.0;
    smlib_brep *r = nullptr;
    CHECK(holed != nullptr && smlib_blend(holed, SMLIB_BLEND_CHORD, back_rim, &size, 1, 0.1, &r, nullptr) != SMLIB_OK && r == nullptr);
    smlib_brep_free(holed);
    smlib_brep_free(cube);
    delete bore;
}

// A run of this program as a child (tangency-labs/solidcore#1173): the corner of rounded_corner
// rounded r2.5, either first thing in the process (`order` "first") or after the failed blend
// above ("after"), written to `path`. 0 if the body is sound.
int corner_child(const char *order, const char *path)
{
    if (std::strcmp(order, "after") == 0)
    {
        failed_chord_blend();
    }
    smlib_brep *corner = rounded_corner(2.5);
    const bool ok = corner != nullptr && sound(corner);
    CHECK(corner != nullptr && smlib_brep_write(corner, path, 1) == SMLIB_OK);
    smlib_brep_free(corner);
    return ok && failures == 0 ? 0 : 1;
}

// Runs `self` as a corner_child in `order` and returns what it wrote ("" if it failed); notes
// whether the body it made was sound in `sound_out`.
std::string corner_child_run(const char *self, const char *order, bool &sound_out)
{
    const std::string path = std::string("solidcore-1173-") + order + ".smb";
    std::string command = std::string("\"") + self + "\" --solidcore-1173 " + order + " " + path;
#ifdef _WIN32
    command = "\"" + command + "\"";
#endif
    std::remove(path.c_str());
    sound_out = std::system(command.c_str()) == 0;
    std::ifstream in(path, std::ios::binary);
    std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    in.close();
    std::remove(path.c_str());
    return text;
}

// A 40 x 20 x 20 box less a cylinder r 8 along y, its axis 4 below the top (SolidCraft
// a_curved_edge_blend_ends_in_a_face_at_an_angle): the groove's arc on the front face y = 0
// spans 240 degrees and ends in the top at 60 degrees to it. Rounded r 0.5, the round is a torus
// about the groove's axis. The cylinder is turned so its seam points up, out of the part: the
// arc is one edge, an SmCircle framed on the seam with its analytic domain from 60 to 300
// degrees, and the torus's span (the turn from the start section to the end one) had the 60
// taken off: the torus stopped 60 degrees short of the arc's end and its end in the top found no
// fillet surface there (1001). The round takes the section's corner (0 <= y <= r and
// 8 <= rho <= 8 + r, outside the ball about (r, 8 + r)) at every angle where it lies in the box,
// z = 16 + rho sin(phi) <= 20: pi + 2 asin(4 / rho) of the turn, so with
// F(rho) = pi rho^2 / 2 + rho^2 asin(4 / rho) + 4 sqrt(rho^2 - 16) its volume is
// the integral over y of F(8 + r - sqrt(r^2 - (y - r)^2)) - F(8), 1.8160502. The round's ends
// in the top are intersection curves fitted within the solver's 3D tolerance (4e-5 here) of both
// surfaces, not exact: the volume is within 1e-4 of its size of that (SMLib's is 5e-5 off).
void groove_arc_round_ends_in_the_top()
{
    smlib_brep *slab = box(0, 0, 0, 40, 20, 20), *bore = cylinder(0, 0, 0, 8.0, 22.0), *along_y = nullptr, *groove = nullptr;
    const double m[4][4] = {{0, 0, 1, 0}, {1, 0, 0, 0}, {0, 1, 0, 0}, {20, -1, 16, 1}};
    CHECK(bore != nullptr && smlib_brep_transform(bore, m, &along_y) == SMLIB_OK);
    CHECK(along_y != nullptr && smlib_boolean(slab, along_y, SMLIB_BOOLEAN_SUBTRACT, &groove) == SMLIB_OK);
    const double r = 0.5, arc[3] = {20, 0, 8};
    auto F = [](double rho) { return PI * rho * rho / 2.0 + rho * rho * std::asin(4.0 / rho) + 4.0 * std::sqrt(rho * rho - 16.0); };
    // y = r - r sin t, t from 0 to pi/2: the corner's outer edge 8 + r - r cos t, dy = r cos t dt.
    const int n = 20000;
    double taken = 0.0;
    for (int i = 0; i < n; ++i)
    {
        const double t = (i + 0.5) / n * PI / 2.0;
        taken += (F(8.0 + r - r * std::cos(t)) - F(8.0)) * r * std::cos(t) * (PI / 2.0) / n;
    }
    smlib_brep *round = nullptr;
    CHECK(groove != nullptr && smlib_fillet(groove, arc, 1, 0.1, r, &round, nullptr) == SMLIB_OK);
    if (round != nullptr)
    {
        int64_t before[3] = {0, 0, 0}, after[3] = {0, 0, 0};
        CHECK(smlib_brep_counts(groove, before) == SMLIB_OK && smlib_brep_counts(round, after) == SMLIB_OK);
        CHECK(after[0] == before[0] + 1);
        CHECK(sound(round));
        double v0 = 0.0, v1 = 0.0;
        CHECK(smlib_brep_volume(groove, 1e-12, &v0) == SMLIB_OK && smlib_brep_volume(round, 1e-12, &v1) == SMLIB_OK);
        const double removed = v0 - v1;
        std::printf("groove arc rounded r%g: removed %.9f, expected %.9f\n", r, removed, taken);
        CHECK(near(removed, taken, 1e-4));
    }
    for (smlib_brep *x : {slab, bore, along_y, groove, round})
    {
        smlib_brep_free(x);
    }
}

// SolidCraft rig case fillets/report-fillet-corner-after-round-running-into-cylinder-solver-residues
// (tests/data/corner-after-round.txt, SMLib's own r 2 round of the front top edge of a bar whose
// sketch left residues of 1e-10 mm): the bar's vertical front left edge, its left top edge and
// the first round's end arc on the left face rounded r 2 together. Along the end arc, of the
// ball's own radius, the earlier round's rail collapses onto the arc's axis, and the round there
// is a sphere patch built as a torus whose rail is a point. The residues left that rail point
// 1e-11 off the axis and the surface nearly singular there, not exactly (smlib_brep_check: Bad
// Surf Pole Approximation); the rail point is the axis's point. Removed: the spandrel
// A = r^2 (1 - pi/4) along the vertical edge's 8 and the top edge's 14, and in the corner's cube
// of side r what three rounds take there, r^3 (1 - pi/6), less the first round's 2 A:
// 20 A + 8 (1 - pi/6).
void corner_ball_after_a_round_is_exactly_singular()
{
    smlib_flat *flat = nullptr;
    smlib_brep *bar = nullptr, *round = nullptr;
    CHECK(smlib_flat_read(SMLIB_TEST_DATA "/corner-after-round.txt", &flat) == SMLIB_OK && smlib_brep_import(flat, &bar) == SMLIB_OK);
    smlib_flat_free(flat);
    const double points[9] = {0.4893, 0.0579, 4.0, 0.4893, 9.058, 10.0, 0.4893, 0.6437, 9.4142};
    const double sizes[3] = {2.0, 2.0, 2.0};
    CHECK(bar != nullptr && smlib_blend(bar, SMLIB_BLEND_ROUND, points, sizes, 3, 2.77, &round, nullptr) == SMLIB_OK);
    if (round != nullptr)
    {
        double v0 = 0.0, v1 = 0.0;
        CHECK(smlib_brep_volume(bar, 1e-12, &v0) == SMLIB_OK && smlib_brep_volume(round, 1e-12, &v1) == SMLIB_OK);
        const double want = 20.0 * 4.0 * (1.0 - PI / 4.0) + 8.0 * (1.0 - PI / 6.0);
        std::printf("corner ball after a round: removed %.9f, expected %.9f\n", v0 - v1, want);
        CHECK(sound(round));
        CHECK(near(v0 - v1, want, 1e-8));
    }
    smlib_brep_free(bar);
    smlib_brep_free(round);
}

// Rounds meeting at a face's concave corner beside an unpicked sharp edge are mitred (the owner's
// ruling, 2026-10-30): each runs on past the corner until it meets the other, where SMLib swept a
// 3-sided face round the side edge (four more faces). Each edge loses (convex) or gains (concave)
// the spandrel r^2 (1 - pi/4) along its length; in the quadrant beyond both walls the mitred
// rounds take or add what both of their spandrel prisms hold, the cube r^3 less two quarter
// cylinders plus their Steinmetz eighth, r^3 (5/3 - pi/2). The mitre's curve is SMLib's fitted
// fillet-fillet intersection, as at a box's convex corner (3.5e-5 to 6e-5 mm^3 per corner here):
// within 2e-8 of the body.
//  - A 60 x 40 x 10 plate with a 30 x 20 x 5 pocket (sharp vertical corners), its four rim edges
//    rounded (SolidCraft rig fillets/every-pocket-rim-edges-r2, -r5-as-deep-as-pocket): convex
//    rounds beside the concave corner edge. 15 faces; at r 5, as deep as the pocket, 11.
//  - A 60 x 60 x 5 plate with a 20 x 20 x 15 boss, its four root edges rounded r 4: concave
//    rounds beside the boss's convex corner edge. 15 faces.
void rounds_at_a_concave_face_corner_are_mitred()
{
    smlib_brep *plate = box(0, 0, 0, 60, 40, 10), *cut = box(15, 10, 5, 30, 20, 6), *pocket = nullptr;
    smlib_brep *base = box(0, 0, 0, 60, 60, 5), *block = box(20, 20, 5, 20, 20, 15), *boss = nullptr;
    CHECK(smlib_boolean(plate, cut, SMLIB_BOOLEAN_SUBTRACT, &pocket) == SMLIB_OK);
    CHECK(smlib_boolean(base, block, SMLIB_BOOLEAN_UNION, &boss) == SMLIB_OK);
    const double rim[12] = {30, 10, 10, 30, 30, 10, 15, 20, 10, 45, 20, 10};
    const double root[12] = {30, 20, 5, 40, 30, 5, 30, 40, 5, 20, 30, 5};
    struct Case
    {
        const char *what;
        const smlib_brep *body;
        const double *edges;
        double r, v0, length, sign;
        int64_t faces;
    };
    const Case cases[] = {{"pocket rim", pocket, rim, 2.0, 21000.0, 100.0, -1.0, 15},
                          {"pocket rim", pocket, rim, 5.0, 21000.0, 100.0, -1.0, 11},
                          {"boss root", boss, root, 4.0, 24000.0, 80.0, 1.0, 15}};
    for (const Case &c : cases)
    {
        smlib_brep *round = nullptr;
        CHECK(c.body != nullptr && smlib_fillet(c.body, c.edges, 4, 0.1, c.r, &round, nullptr) == SMLIB_OK);
        if (round == nullptr)
        {
            continue;
        }
        const double want = c.v0 + c.sign * (c.length * c.r * c.r * (1.0 - PI / 4.0) + 4.0 * c.r * c.r * c.r * (5.0 / 3.0 - PI / 2.0));
        int64_t counts[3] = {0, 0, 0};
        double v = 0.0;
        CHECK(smlib_brep_counts(round, counts) == SMLIB_OK && smlib_brep_volume(round, 1e-12, &v) == SMLIB_OK);
        std::printf("%s rounded r%g: %lld faces, volume %.9f, expected %.9f\n", c.what, c.r, static_cast<long long>(counts[0]), v, want);
        CHECK(counts[0] == c.faces);
        CHECK(sound(round));
        CHECK(near(v, want, 2e-8));
        smlib_brep_free(round);
    }
    for (smlib_brep *x : {plate, cut, pocket, base, block, boss})
    {
        smlib_brep_free(x);
    }
}

// SolidCraft rig fillets/1080-fillet-sloped-basin-1-mm-rim-r0-5-4-edges: a 100 x 60 x 30 slab less
// a basin lofted from 20..80 x 1..50 at the top to 30..70 x 20..40 at z = 10, its front and right
// rim edges, then all four, rounded r 0.5. The walls lean differently (the front 19 over 20, the
// others 10 over 20), so at each corner of two rounds the rails meet the sloped corner edge at
// two places: the mitre leaves a gap there on the side face's continuation past the corner edge
// (beyond the wall the smaller round rolls on, inside the material behind the other wall),
// filled by a patch on that wall's plane that joins the wall. SMLib cut the larger round by the
// side face's own surface, which does not reach past the corner edge, and failed (SM_ERR). Each
// round's section at height z spans from its wall to its arc, its width w(z), and runs between its
// neighbours' walls extended, L(z); in the quadrant beyond both walls of a mitred corner the two
// take what both take, w1(z) w2(z): removed Int (sum w L + sum over corners w1 w2) dz, taken over
// the arc's angle from the top (z = 30 - r + r cos phi, smooth there) and split where each round
// meets its wall (Simpson, converged to 1e-10 mm^3). V0 = 144533.333333 (a prismatoid:
// 20/6 (A_top + A_bottom + 4 A_mid)).
// 13 faces, then 15: the gaps' patches join the walls.
void basin_rim_rounds_mitred_with_a_gap_at_the_corner()
{
    smlib_brep *slab = box(0, 0, 0, 100, 60, 30), *pocket = nullptr, *basin = nullptr;
    const double rims[2][5] = {{20, 1, 80, 50, 30}, {30, 20, 70, 40, 10}};
    int64_t kinds[8];
    double reals[48];
    for (int s = 0; s < 2; ++s)
    {
        const double *q = rims[s];
        const double xy[4][2] = {{q[0], q[1]}, {q[2], q[1]}, {q[2], q[3]}, {q[0], q[3]}};
        for (int i = 0; i < 4; ++i)
        {
            const int j = (i + 1) % 4;
            const double seg[6] = {xy[i][0], xy[i][1], q[4], xy[j][0], xy[j][1], q[4]};
            kinds[4 * s + i] = SMLIB_CURVE_LINE;
            std::copy(seg, seg + 6, reals + 6 * (4 * s + i));
        }
    }
    const int64_t sizes[2] = {4, 4};
    CHECK(smlib_loft(kinds, 8, reals, 48, sizes, 2, 1, 1, &pocket) == SMLIB_OK);
    CHECK(pocket != nullptr && smlib_boolean(slab, pocket, SMLIB_BOOLEAN_SUBTRACT, &basin) == SMLIB_OK);
    const double r = 0.5, front = 19.0 / 20.0, side = 10.0 / 20.0;
    // A round of a rim edge whose wall leans k (horizontally per unit down): the wall's line at
    // depth s = 30 - z is k s from the rim, the arc's centre r below the top and r from the wall.
    auto width = [r](double k, double z) {
        const double q = std::sqrt(1.0 + k * k), uc = k * r - r * q, zc = 30.0 - r;
        if (z < zc + r * k / q || z > 30.0)
        {
            return 0.0;
        }
        return k * (30.0 - z) - (uc + std::sqrt(std::max(0.0, r * r - (z - zc) * (z - zc))));
    };
    auto simpson = [](const std::function<double(double)> &f, double a, double b) {
        const int n = 2000;
        const double h = (b - a) / n;
        double s = f(a) + f(b);
        for (int i = 1; i < n; ++i)
        {
            s += (i % 2 ? 4.0 : 2.0) * f(a + i * h);
        }
        return s * h / 3.0;
    };
    // the edges' leans and lengths between their neighbours' walls: front, right, back, left
    const double k[4] = {front, side, side, side};
    auto length = [&](int e, double z) {
        const double s = 30.0 - z;
        return e % 2 == 0 ? 60.0 - 2.0 * side * s : 49.0 - (front + side) * s;
    };
    const double edges[12] = {50, 1, 30, 80, 27, 30, 50, 50, 30, 20, 27, 30};
    const double v0 = 3000.0 * 60.0 - 20.0 / 6.0 * (60.0 * 49.0 + 40.0 * 20.0 + 4.0 * 50.0 * 34.5);
    for (int n : {2, 4})
    {
        auto removed = [&](double z) {
            double sum = 0.0;
            for (int e = 0; e < n; ++e)
            {
                sum += width(k[e], z) * length(e, z);
            }
            for (int e = 0; e < 4; ++e)
            {
                const int f = (e + 1) % 4;
                if (e < n && f < n)
                {
                    sum += width(k[e], z) * width(k[f], z);
                }
            }
            return sum;
        };
        // over the arc's angle phi from the top, split where each round meets its wall
        // (cos phi = k / sqrt(1 + k^2)): dz = r sin phi dphi
        auto by_angle = [&](double phi) { return removed(30.0 - r + r * std::cos(phi)) * r * std::sin(phi); };
        const double at_front = std::acos(front / std::sqrt(1.0 + front * front)), at_side = std::acos(side / std::sqrt(1.0 + side * side));
        const double want = v0 - simpson(by_angle, 0.0, at_front) - simpson(by_angle, at_front, at_side);
        smlib_brep *round = nullptr;
        CHECK(basin != nullptr && smlib_fillet(basin, edges, n, 0.1, r, &round, nullptr) == SMLIB_OK);
        if (round == nullptr)
        {
            continue;
        }
        int64_t counts[3] = {0, 0, 0};
        double v = 0.0;
        CHECK(smlib_brep_counts(round, counts) == SMLIB_OK && smlib_brep_volume(round, 1e-12, &v) == SMLIB_OK);
        std::printf("sloped basin, %d rim edges rounded r%g: %lld faces, volume %.9f, expected %.9f\n", n, r, static_cast<long long>(counts[0]), v,
                    want);
        CHECK(counts[0] == 11 + n);
        CHECK(sound(round));
        CHECK(near(v, want, 1e-9));
        smlib_brep_free(round);
    }
    for (smlib_brep *x : {slab, pocket, basin})
    {
        smlib_brep_free(x);
    }
}

// SolidCraft issue_871_tests::rounds_that_fit_are_closed: a 14 x 21.5 x 10 plate with a narrower
// block on top (y from a to b, 10 high) and a 5.2 x 7.1 window through both, every edge of both
// end faces (x = 9.5 and 23.5) rounded r 0.5. The rounds meeting at the end faces' concave
// corners (the two shoulders, the window's four corners) are mitred: their rails run on past
// those corners until they meet. The check that a blend's rail stays out of the band another
// blend takes off a face (BandsClearOtherRails) read the rails before they were trimmed, and
// the window's top rail, run on past its corner, crossed the band of the block's side, 1.26 away
// (SMLIB_ERR_NO_ROOM). The band is the rail's part between its end vertices. Each end face's
// round section at depth u into the round is the end face's outline moved in by
// d = r - sqrt(r^2 - (r - u)^2), every corner mitred: the area less L d (L = 107.6, the outline's
// and the window's length together), plus d^2 at each of the six convex corners and less d^2 at
// each of the six concave ones, so both ends remove 2 L r^2 (1 - pi/4), to the fit of the 24
// corners' bevel curves (1e-7 of the body). At r 1 the rounds of the window and of the block's
// sides overlap across the 1.26 wall: refused.
void end_face_rounds_mitred_at_concave_corners()
{
    const double a = 12.74437365538526, b = 22.69275963366936;
    auto block = [](double x0, double y0, double z0, double x1, double y1, double z1) { return box(x0, y0, z0, x1 - x0, y1 - y0, z1 - z0); };
    smlib_brep *base = block(9.5, 8, 0, 23.5, 29.5, 10), *top = block(9.5, a, 10, 23.5, b, 20), *joined = nullptr;
    smlib_brep *window = block(9.5, 14, 7, 23.5, 19.2, 14.1), *part = nullptr;
    CHECK(smlib_boolean(base, top, SMLIB_BOOLEAN_UNION, &joined) == SMLIB_OK);
    CHECK(joined != nullptr && smlib_boolean(joined, window, SMLIB_BOOLEAN_SUBTRACT, &part) == SMLIB_OK);
    const double yz[12][2] = {{a, 15}, {10.37218682769263, 10}, {8, 5}, {18.75, 0}, {29.5, 5}, {26.09637981683468, 10},
                              {b, 15}, {17.71856664452731, 20}, {16.6, 7}, {14, 10.55}, {16.6, 14.1}, {19.2, 10.55}};
    double edges[72];
    for (int i = 0; i < 24; ++i)
    {
        edges[3 * i] = i < 12 ? 23.5 : 9.5;
        edges[3 * i + 1] = yz[i % 12][0];
        edges[3 * i + 2] = yz[i % 12][1];
    }
    const double v0 = 14.0 * (21.5 * 10.0 + (b - a) * 10.0 - 5.2 * 7.1), r = 0.5;
    const double want = v0 - 2.0 * (2.0 * (21.5 + 20.0) + 2.0 * (5.2 + 7.1)) * r * r * (1.0 - PI / 4.0);
    smlib_brep *round = nullptr, *overlap = nullptr;
    CHECK(part != nullptr && smlib_fillet(part, edges, 24, 0.1, r, &round, nullptr) == SMLIB_OK);
    if (round != nullptr)
    {
        double v = 0.0;
        CHECK(smlib_brep_volume(round, 1e-12, &v) == SMLIB_OK);
        std::printf("end faces' edges rounded r%g: volume %.9f, expected %.9f\n", r, v, want);
        CHECK(sound(round));
        CHECK(near(v, want, 1e-7));
    }
    CHECK(part != nullptr && smlib_fillet(part, edges, 24, 0.1, 1.0, &overlap, nullptr) == SMLIB_ERR_NO_ROOM && overlap == nullptr);
    for (smlib_brep *x : {base, top, joined, window, part, round})
    {
        smlib_brep_free(x);
    }
}

// SolidCraft tests::fillet_tabbed_plate_outer_loop_and_partial_selection_after_hole: a plate 2
// thick, its outline an R25 arc to the sides x = +-4 of a tab whose top is 4 wide with R2
// corners, every top edge of the outline rounded r 0.5. At the concave corners (+-4, y) the arc's
// round and the tab side's are mitred: the trace of the arc's round runs on past the arc's end
// to where its rails meet the tab side's, on the R25 wall's continuation. The wall is an
// SmCylinder, and its continuation was a tangent (G1R) extension of its NURB, 0.0036 off the
// circle 0.44 past the arc's end: the round's section came from there, every section of the
// torus was that far off, and the mitre's curve found no fillet surface through its ends
// (1001). Analytic surfaces now continue as themselves. Removed: the outline's area less the
// area of the outline moved in by d = r (1 - cos phi) at the height r sin phi into the round
// (each element parallel to itself, the arc's radius 25 - d, the tab's run moved in by d and
// the moved elements meeting where they cross at the mitres), over that height (midpoint rule
// in phi, 8192 steps, which 4096 agree with to 1e-7); the area by Green's theorem. The fit of the
// rounds' and mitres' curves leaves 1.6e-7 of it.
void tabbed_plate_rounds_mitred_on_the_arcs_continuation()
{
    const double theta = std::acos(4.0 / 25.0), y = std::sqrt(625.0 - 16.0), cy = 280.0 / 9.0;
    Profile p;
    p.arc(0, 0, 0, 1, 25, PI - theta, 2.0 * PI + theta).line(4, y, 4, cy, 0).arc(2, cy, 0, 1, 2, 0, PI / 2.0);
    p.line(2, cy + 2.0, -2, cy + 2.0, 0).arc(-2, cy, 0, 1, 2, PI / 2.0, PI).line(-4, cy, -4, y, 0);
    smlib_brep *plate = p.extrude(2.0), *round = nullptr;
    // the outline's area with the arc moved in by da and the tab's run by dt
    auto area = [cy](double da, double dt) {
        const double big = 25.0 - da, rho = 2.0 - dt, tab = rho * rho * PI + 2.0 * rho * (2.0 + cy) + 4.0 * (cy + rho);
        const double yc = std::sqrt(big * big - (4.0 - dt) * (4.0 - dt)), at = std::atan2(yc, 4.0 - dt);
        return 0.5 * (big * big * (PI + 2.0 * at) + tab) + (4.0 - dt) * (cy - yc);
    };
    const double r = 0.5, dphi = (PI / 2.0) / 8192.0;
    double want = 0.0;
    for (int i = 0; i < 8192; ++i)
    {
        const double phi = (i + 0.5) * dphi, d = r * (1.0 - std::cos(phi));
        want += (area(0.0, 0.0) - area(d, d)) * r * std::cos(phi) * dphi;
    }
    const double top[18] = {0, -25, 2, 4, 27.8945, 2, 3.41421356, 32.5253247, 2, 0, cy + 2.0, 2, -3.41421356, 32.5253247, 2, -4, 27.8945, 2};
    CHECK(plate != nullptr && near(volume_of(plate), 2.0 * area(0.0, 0.0), 1e-10));
    CHECK(plate != nullptr && smlib_fillet(plate, top, 6, 0.01, r, &round, nullptr) == SMLIB_OK);
    if (round != nullptr)
    {
        const double removed = volume_of(plate) - volume_of(round);
        std::printf("tabbed plate's outline rounded r%g: removed %.9f, expected %.9f\n", r, removed, want);
        CHECK(sound(round));
        CHECK(near(removed, want, 1e-6));
    }
    smlib_brep_free(plate);
    smlib_brep_free(round);
}

// SolidCraft blend_corner::tests::issue1557_connected_rounds_keep_their_radii (its reference): a
// 25 x 25 x 10 box, its four bottom edges rounded r 2 and the vertical edge (25, 0) r 0.5 in one
// call. SMLib's own check found a 6.9e-5 gap between the corner patch's edge (23, 2, 0)-(24.5, 0, 2)
// and the patch, on which it lies to 4e-15: the check drops the edge's points onto the patch
// with the last foot point as a guess, and the local solve stopped 0.0016 short in u, where the
// patch's u runs at 0.043 mm per unit near its collapsed side; the drop took any point less than
// 1e-4 mm off along the surface for the foot point. A foot point is now off along the surface by
// no more than the face's tolerance. The body is unchanged.
void corner_patch_edge_lies_on_its_patch()
{
    smlib_brep *b = box(0, 0, 0, 25, 25, 10), *r = nullptr;
    const double points[15] = {12.5, 0, 0, 25, 12.5, 0, 12.5, 25, 0, 0, 12.5, 0, 25, 0, 6};
    const double sizes[5] = {2, 2, 2, 2, 0.5};
    CHECK(b != nullptr && smlib_blend(b, SMLIB_BLEND_ROUND, points, sizes, 5, 0.1, &r, nullptr) == SMLIB_OK);
    if (r != nullptr)
    {
        int64_t counts[3] = {0, 0, 0};
        CHECK(smlib_brep_counts(r, counts) == SMLIB_OK && counts[0] == 12);
        CHECK(sound(r));
    }
    smlib_brep_free(b);
    smlib_brep_free(r);
}

// A dropped point's foot from a guess is one the local solve converged to. The drop took a local
// solve's result for the foot when it was off along the surface by less than the face's
// tolerance, so on a face of a looser tolerance a stalled solve was taken (on a raked wedge one
// 6.05e-4 mm off, read as a gap of 7.07e-4 where it is about 3.7e-4). The corner patch of the
// #1557 box above, its tolerance set to 1e-3: the edge (23, 2, 0)-(24.5, 0, 2) lies on it, and its
// points dropped onto it with the last foot point as the guess (as SMLib's check drops them),
// first from the guess (0.01, 0) where the local solve stalls 6.9e-5 off, are on it.
void corner_patch_drop_from_a_guess_is_its_foot()
{
    smlib_brep *b = box(0, 0, 0, 25, 25, 10), *r = nullptr;
    const double points[15] = {12.5, 0, 0, 25, 12.5, 0, 12.5, 25, 0, 0, 12.5, 0, 25, 0, 6};
    const double sizes[5] = {2, 2, 2, 2, 0.5};
    CHECK(b != nullptr && smlib_blend(b, SMLIB_BLEND_ROUND, points, sizes, 5, 0.1, &r, nullptr) == SMLIB_OK && r != nullptr);
    int patches = 0;
    double worst = HUGE_VAL;
    if (r != nullptr)
    {
        worst = 0.0;
        SmTArray<SmEdge *> es;
        r->p->GetEdges(es);
        const SmPoint3d a(23, 2, 0), z(24.5, 0, 2);
        for (ULONG i = 0; i < es.GetSize(); ++i)
        {
            const SmPoint3d p0 = es[i]->GetStartVertex()->GetPoint(), p1 = es[i]->GetEndVertex()->GetPoint();
            if (!((p0.DistanceBetween(a) < 1e-6 && p1.DistanceBetween(z) < 1e-6) || (p0.DistanceBetween(z) < 1e-6 && p1.DistanceBetween(a) < 1e-6)))
            {
                continue;
            }
            SmTArray<SmFace *> fs;
            es[i]->GetFaces(fs);
            for (ULONG j = 0; j < fs.GetSize(); ++j)
            {
                SmSurface *srf = fs[j]->GetSurface();
                if (srf == nullptr || srf->IsKindOf(SmCone_TYPE) || srf->IsKindOf(SmTorus_TYPE) || srf->IsKindOf(SmPlane_TYPE))
                {
                    continue; // the patch is a plain B-spline
                }
                ++patches;
                fs[j]->SetTolerance(SM_CAST_ZONETOL3D(1.0e-3), FALSE, FALSE);
                const SmExtent2d dom = srf->GetNaturalUVDomain();
                const SmExtent1d ivl = es[i]->GetInterval();
                SmPoint2d guess(0.01, 0.0);
                for (int k = 0; k <= 200; ++k)
                {
                    SmPoint3d q;
                    CHECK(es[i]->GetCurve()->EvaluatePoint(ivl.Evaluate(k / 200.0), q) == SM_SUCCESS);
                    SmBoolean ok = FALSE, multi = FALSE;
                    SmPoint2d uv;
                    double gap = HUGE_VAL;
                    CHECK(srf->DropPoint(q, dom, &guess, ok, uv, gap, multi) == SM_SUCCESS && ok);
                    SmPoint3d on;
                    CHECK(srf->EvaluatePoint(uv, on) == SM_SUCCESS);
                    worst = std::max(worst, on.DistanceBetween(q));
                    guess = uv;
                }
            }
        }
    }
    std::printf("#1557 corner patch, tolerance 1e-3: its edge's points dropped from the last foot are on it to %.3g\n", worst);
    CHECK(patches == 1);
    CHECK(worst < 1e-9);
    smlib_brep_free(b);
    smlib_brep_free(r);
}

} // namespace

// A disc r 15, 5 thick, one circular edge rounded r 1.5 and shelled 1.6 or 2 inward, open at one
// end (SolidCraft shell_of_a_disc_rounded_smaller_than_the_wall, #319). The round is a fillet's
// surface of revolution of an arc; moved in past its radius it has no offset, as a torus's (a
// round thinner than the wall), and the floor's and side's offsets meet across it in a circle.
// The arc's offset used to be built as the arc turned inside out, r 0.1 or 0.5 round its centre,
// and the offset was no closed solid (1001). The cavity is the cylinder r 15 - t, 5 - t deep; the
// round takes r^2 (1 - pi/4) from the disc's section, its centroid r (10 - 3 pi) / (3 (4 - pi)) in
// from the rim, swept round the axis (Pappus).
void shelled_disc_round_thinner_than_the_wall()
{
    const double rd = 15.0, h = 5.0, r = 1.5;
    const double corner = r * r * (1.0 - PI / 4.0) * 2.0 * PI * (rd - r * (10.0 - 3.0 * PI) / (3.0 * (4.0 - PI)));
    const double top[3] = {0, 0, 0}, bottom[3] = {0, 0, -h};
    struct Case
    {
        double edge_z;
        const double *open;
    };
    for (Case c : {Case{-h, top}, Case{-h, bottom}, Case{0.0, bottom}})
    {
        smlib_brep *disc = cylinder(0, 0, -h, rd, h), *rounded = nullptr;
        const double edge[3] = {0, rd, c.edge_z};
        CHECK(smlib_fillet(disc, edge, 1, 0.1, r, &rounded, nullptr) == SMLIB_OK);
        for (double t : {1.6, 2.0})
        {
            smlib_brep *s = nullptr;
            CHECK(rounded != nullptr && smlib_shell(rounded, c.open, 1, 0.1, -t, &s, nullptr) == SMLIB_OK);
            if (s != nullptr)
            {
                const double want = PI * rd * rd * h - corner - PI * (rd - t) * (rd - t) * (h - t);
                std::printf("disc rounded at z %g, open at z %g, shelled %g: volume %.9f, expected %.9f\n", c.edge_z, c.open[2], t, volume_of(s), want);
                CHECK(sound(s));
                CHECK(near(volume_of(s), want, 1e-7));
            }
            smlib_brep_free(s);
        }
        smlib_brep_free(disc);
        smlib_brep_free(rounded);
    }
}

// The lines `lines` (each from its first point to its second) ordered into loops by
// SmTrimmingTools::OrderCurvesIntoLoops: its status, and each loop's area as it runs (shoelace
// over its lines), smallest first, each loop checked closed.
static SmStatus ordered_loop_areas(const std::vector<std::pair<SmPoint3d, SmPoint3d>> &lines, std::vector<double> &areas)
{
    const SmContext &ctx = *SmApiGetOrCreateContext();
    areas.clear();
    SmTArray<SmCurve *> curves;
    for (const auto &l : lines)
    {
        SmLine *c = nullptr;
        CHECK(SmLine::CreateLineSegment(ctx, 3, l.first, l.second, c) == SM_SUCCESS && c != nullptr);
        curves.Add(c);
    }
    SmTArray<SmCurve *> ordered;
    SmTArray<ULONG> counts;
    SmTArray<SmOrientType> orients;
    SmPoint3d at;
    SmVector3d normal;
    SmExtent3d box;
    SmStatus s = SmTrimmingTools::OrderCurvesIntoLoops(curves, 1e-6, ordered, counts, orients, at, normal, box);
    if (s == SM_SUCCESS && ordered.GetSize() == curves.GetSize())
    {
        ULONG first = 0;
        for (ULONG c = 0; c < counts.GetSize(); ++c)
        {
            double area = 0.0;
            SmPoint3d end;
            for (ULONG k = first; k < first + counts[c]; ++k)
            {
                SmPoint3d a, b;
                ordered[k]->GetEnds(a, b);
                if (orients[k] != SM_OT_SAME)
                {
                    std::swap(a, b);
                }
                CHECK(k == first || a.DistanceBetween(end) < 1e-12);
                area += 0.5 * (a.x * b.y - b.x * a.y);
                end = b;
            }
            SmPoint3d start, other;
            ordered[first]->GetEnds(start, other);
            CHECK(end.DistanceBetween(orients[first] == SM_OT_SAME ? start : other) < 1e-12);
            areas.push_back(area);
            first += counts[c];
        }
        std::sort(areas.begin(), areas.end());
    }
    else if (s == SM_SUCCESS)
    {
        s = SM_ERR;
    }
    for (ULONG k = 0; k < curves.GetSize(); ++k)
    {
        delete curves[k];
    }
    return s;
}

// The lines of the closed polygon `p` (n corners), reversed or not, added to `lines`.
static void polygon_lines(std::vector<std::pair<SmPoint3d, SmPoint3d>> &lines, const double (*p)[2], int n, bool reversed)
{
    for (int i = 0; i < n; ++i)
    {
        SmPoint3d a(p[i][0], p[i][1], 0.0), b(p[(i + 1) % n][0], p[(i + 1) % n][1], 0.0);
        lines.push_back(reversed ? std::make_pair(b, a) : std::make_pair(a, b));
    }
}

// The list `lines` rotated by `turn` (its curve `turn` first).
static std::vector<std::pair<SmPoint3d, SmPoint3d>> rotated(const std::vector<std::pair<SmPoint3d, SmPoint3d>> &lines, size_t turn)
{
    std::vector<std::pair<SmPoint3d, SmPoint3d>> out;
    for (size_t k = 0; k < lines.size(); ++k)
    {
        out.push_back(lines[(k + turn) % lines.size()]);
    }
    return out;
}

// SolidCraft rig shell/616-shell-open-a-rib-side-face-of-a-lid (tests/data/lid-rib-open-at-its-
// side.txt, as smlib_shell was handed it): #616's lid shelled 0.5 open at a rib's outer side
// (x = 102.8), which stands 0.15 in from the skirt's outer wall (102.65). The skirt top's offset
// (z = 4.9) has the rib's foot as a hole; the hole's side along the open face stays at 102.8, the
// outer loop's side along the skirt's wall moves to 103.15, and between them, under the opening,
// the loops wind round the other way: the cavity's floor, the cavity above it, a face of the skirt
// top's offset turned over. It was dropped and the cavity left open (1001). The owner's ruling of
// 2026-10-31: it builds, with the 0.15 lip beside the opening.
static void shell_open_beside_a_step_keeps_the_floor_under_the_opening()
{
    smlib_flat *f = nullptr;
    smlib_brep *part = nullptr, *r = nullptr;
    CHECK(smlib_flat_read(SMLIB_TEST_DATA "/lid-rib-open-at-its-side.txt", &f) == SMLIB_OK && smlib_brep_import(f, &part) == SMLIB_OK);
    smlib_flat_free(f);
    const double pick[3] = {102.8, 17, 6.85};
    CHECK(part != nullptr && smlib_shell(part, pick, 1, 0.11862452211935036, -0.5, &r, nullptr) == SMLIB_OK);
    if (r != nullptr)
    {
        CHECK(sound(r));
        SmTArray<SmFace *> faces;
        r->p->GetFaces(faces);
        int floors = 0;
        for (ULONG i = 0; i < faces.GetSize(); i++)
        {
            SmExtent3d e;
            CHECK(faces[i]->CalculateBoundingBox(e, TRUE) == SM_SUCCESS);
            if (std::fabs(e.GetMin(0) - 102.8) < 1e-6 && std::fabs(e.GetMax(0) - 103.15) < 1e-6 && std::fabs(e.GetMin(2) - 4.9) < 1e-6 &&
                std::fabs(e.GetMax(2) - 4.9) < 1e-6)
            {
                floors++;
            }
        }
        std::printf("lid shelled 0.5 open at a rib's side: volume %.9f, %d floor face(s) under the opening\n", volume_of(r), floors);
        CHECK(floors == 1);
        done(r);
    }
    smlib_brep_free(part);
}

// SolidCraft rig shell/cov-shell-trough-whose-floor-falls-to-a-drain and
// cov-shell-pocket-whose-floor-falls-to-a-drain (tests/data/trough-falling-to-a-drain.txt,
// pocket-falling-to-a-drain.txt, as smlib_shell was handed them): a 100 x 60 x 40 trough whose
// floor is four facets falling 1 degree to a drain at its centre, a tab joined on one end,
// shelled 2 open at the top; and a 140 x 90 x 50 slab less a pocket whose floor falls so, shelled
// 3 open at the top. Four faces meet at each corner of the floor (two walls, two facets) and at
// the drain (four facets), and their offsets do not meet in one point: each offset vertex is two,
// joined by a new edge. The volumes are the solutions' own closed forms (the offset faces
// extended to meet: the cavity integrated exactly over the cells where each moved facet is
// highest); both were 1001.
static void four_faces_meeting_at_a_vertex_shell()
{
    struct
    {
        const char *data;
        double pick[3], tolerance, wall, volume;
    } cases[2] = {{SMLIB_TEST_DATA "/trough-falling-to-a-drain.txt", {75.458529763090851, 18.269689917857015, 40}, 0.14015151326053074, 2.0,
                   37962.27221527102},
                  {SMLIB_TEST_DATA "/pocket-falling-to-a-drain.txt", {13.333333333333332, 60, 50}, 0.17378247196982766, 3.0, 151815.7900294403}};
    for (const auto &c : cases)
    {
        smlib_flat *f = nullptr;
        smlib_brep *part = nullptr, *r = nullptr;
        CHECK(smlib_flat_read(c.data, &f) == SMLIB_OK && smlib_brep_import(f, &part) == SMLIB_OK);
        smlib_flat_free(f);
        CHECK(part != nullptr && smlib_shell(part, c.pick, 1, c.tolerance, -c.wall, &r, nullptr) == SMLIB_OK);
        if (r != nullptr)
        {
            double v = volume_of(r);
            std::printf("%s shelled %g: %.9f, the closed form %.9f\n", c.data, c.wall, v, c.volume);
            CHECK(near(v, c.volume, 1e-8));
            done(r);
        }
        smlib_brep_free(part);
    }
}

// SolidCraft rig shell/cov-shell-plate-with-a-hole-through-its-arc-side (tangency-labs/solidcore#1795;
// tests/data/plate-hole-through-arc-side.txt, as smlib_shell was handed it): a 42.5 x 20 x 27.5
// plate closed by a half-round side (r10 about z through (0, 10)), a blind hole r5 drilled in
// along x through that side (axis y = z = 12.5, bottom at x = 7.5), shelled 2 open at the far
// end (x = 42.5). The hole's offset (r7) pokes out through the half-round's offset (r8) where
// y > 18: their intersection is no longer one closed curve, so the two offsets cannot be mated
// along the edge they came from (MateEdges); the merge then has to intersect them, and it was
// told they already shared that edge (1001). The cavity is the plate inset by 2 (the half disc
// r8 and the rectangle 0..42.5 x 2..18, z 2..25.5) less the hole's offset where it lies in it:
// over the disc r7, y <= 18, x from the r8 side to 9.5.
static void hole_whose_offset_leaves_the_round_side_shell()
{
    auto simpson = [](auto f, double a, double b)
    {
        const int n = 20000;
        double h = (b - a) / n, t = f(a) + f(b);
        for (int i = 1; i < n; i++)
        {
            t += (i % 2 ? 4.0 : 2.0) * f(a + i * h);
        }
        return t * h / 3.0;
    };
    // The hole's offset in the cavity, its length over the disc integrated in y (substituted at
    // either end, where the square roots fall to zero, so the integrand is smooth).
    auto across = [](double y)
    {
        return 2.0 * std::sqrt(std::max(0.0, 49.0 - (y - 12.5) * (y - 12.5))) * (9.5 + std::sqrt(std::max(0.0, 64.0 - (y - 10.0) * (y - 10.0))));
    };
    double tube = simpson([&](double v) { return across(5.5 + v * v) * 2.0 * v; }, 0.0, std::sqrt(6.5)) +
                  simpson([&](double u) { return across(18.0 - u * u) * 2.0 * u; }, 0.0, std::sqrt(6.0));
    double cavity = (PI * 32.0 + 42.5 * 16.0) * 23.5 - tube;
    smlib_flat *f = nullptr;
    smlib_brep *plate = nullptr, *r = nullptr;
    CHECK(smlib_flat_read(SMLIB_TEST_DATA "/plate-hole-through-arc-side.txt", &f) == SMLIB_OK && smlib_brep_import(f, &plate) == SMLIB_OK);
    smlib_flat_free(f);
    const double end[3] = {42.5, 3.333333333333333, 21.388888888888886};
    CHECK(plate != nullptr && smlib_shell(plate, end, 1, 0.062550980015984026, -2.0, &r, nullptr) == SMLIB_OK);
    if (r != nullptr && plate != nullptr)
    {
        double v = volume_of(r), v0 = volume_of(plate);
        std::printf("plate with a hole through its arc side shelled 2: %.9f, the plate less its cavity %.9f\n", v, v0 - cavity);
        CHECK(near(v, v0 - cavity, 1e-6));
        int64_t counts[3] = {0, 0, 0};
        CHECK(smlib_brep_counts(r, counts) == SMLIB_OK && counts[0] == 15);
        done(r);
    }
    smlib_brep_free(plate);
}

// SolidCraft rig shell/545-shell-ring-cut-frame-open-at-an-outer-wall
// (tests/data/ring-cut-frame-thin-floor.txt, as smlib_shell was handed it): a 100 x 100 x 30
// frame with a recess x, y 6..94 down to z 18 from the top and a ring cut between the squares
// 2.5..97.5 and 16..84 up to z 17.5 from the bottom, a 0.5 floor between them, shelled 0.5 open
// at the ring cut's outer wall x = 97.5. The floor's two offsets cross; the region between them
// is excess, not cavity, and the cavity is two pieces: the core's [16.5, 83.5]^2 x [0.5, 17.5]
// and the frame's (the ring between the squares 0.5..99.5 and 2..98 for z 0.5..18, the ring
// between 0.5..99.5 and 5.5..94.5 for z 18..29.5, and x 97.5..98, y 2..98, z 0.5..18 behind the
// open wall). It was 1001 (the region between the floor's offsets kept, touching both pieces
// along edges).
static void floor_thinner_than_two_walls_shell()
{
    smlib_flat *f = nullptr;
    smlib_brep *frame = nullptr, *r = nullptr;
    CHECK(smlib_flat_read(SMLIB_TEST_DATA "/ring-cut-frame-thin-floor.txt", &f) == SMLIB_OK && smlib_brep_import(f, &frame) == SMLIB_OK);
    smlib_flat_free(f);
    const double open[3] = {97.5, 89.583333333333329, 11.666666666666666};
    CHECK(frame != nullptr && smlib_shell(frame, open, 1, 0.14456932294800962, -0.5, &r, nullptr) == SMLIB_OK);
    if (r != nullptr)
    {
        const double cavity = 67.0 * 67.0 * 17.0 + (99.0 * 99.0 - 96.0 * 96.0) * 17.5 + (99.0 * 99.0 - 89.0 * 89.0) * 11.5 + 0.5 * 96.0 * 17.5;
        double v = volume_of(r);
        std::printf("ring-cut frame shelled 0.5: %.9f, the closed form %.9f\n", v, 130054.5 - cavity);
        CHECK(near(v, 130054.5 - cavity, 1e-8));
        done(r);
    }
    smlib_brep_free(frame);
}

// SolidCraft rig shell/cov-shell-up-down-joined-squares-open-at-a-side-wall
// (tests/data/up-down-joined-squares.txt): four 10 x 10 squares about the origin, three extruded
// 10 up and three 10 down and joined, so near the origin the body is space less two opposite
// octants, shelled 0.346 open at a side wall (y = 0, x 0..10, z 0..10). The offsets round the two
// empty octants overlap in a box (-t, t) x (0, t) x (-t, t) that the merge closes, touching the
// cavity along six edges; it is no cavity (1001 while it was kept). With every face axis-aligned
// the sharp offset is the body less every point whose cube of half-side t leaves the body other
// than through the open face; integrated exactly over the cells the planes 0, +-t, 2t, +-10 and
// +-(10 - t) cut, the cavity is 5232.766153056.
static void up_down_joined_squares_shell()
{
    smlib_flat *f = nullptr;
    smlib_brep *part = nullptr, *r = nullptr;
    CHECK(smlib_flat_read(SMLIB_TEST_DATA "/up-down-joined-squares.txt", &f) == SMLIB_OK && smlib_brep_import(f, &part) == SMLIB_OK);
    smlib_flat_free(f);
    const double side[3] = {8.3333333333333321, 0, 6.6666666666666661};
    CHECK(part != nullptr && smlib_shell(part, side, 1, 0.034642016151377554, -0.346, &r, nullptr) == SMLIB_OK);
    if (r != nullptr)
    {
        double v = volume_of(r);
        std::printf("up-down joined squares shelled 0.346: %.9f, the closed form %.9f\n", v, 6000.0 - 5232.766153056);
        CHECK(near(v, 6000.0 - 5232.766153056, 1e-9));
        done(r);
    }
    smlib_brep_free(part);
}

// A 10 x 6 plate's outline and a triangular hole whose corner touches it at (0, 3), as eight
// lines, ordered into loops (SmTrimmingTools::OrderCurvesIntoLoops) from every rotation of the
// list, with the outline or the hole reversed. Four curve ends meet at (0, 3). Each loop must
// close and never cross itself there: either the outline and the hole apart (areas 60 and 1) or
// one loop through the corner twice round the material (area 59), never the hole run the wrong
// way round inside the outline (61). Each line's own direction is kept, so reversing an outline
// or hole reverses its lines.
static void loops_touching_at_a_vertex_in_any_order()
{
    const double outline[5][2] = {{0, 0}, {10, 0}, {10, 6}, {0, 6}, {0, 3}}, hole[3][2] = {{0, 3}, {2, 3.5}, {2, 2.5}};
    for (int variant = 0; variant < 3; ++variant)
    {
        std::vector<std::pair<SmPoint3d, SmPoint3d>> lines;
        polygon_lines(lines, outline, 5, variant == 1);
        polygon_lines(lines, hole, 3, variant == 2);
        for (size_t turn = 0; turn < lines.size(); ++turn)
        {
            std::vector<double> areas;
            SmStatus s = ordered_loop_areas(rotated(lines, turn), areas);
            for (double &a : areas)
            {
                a = std::fabs(a);
            }
            std::sort(areas.begin(), areas.end());
            bool apart = areas.size() == 2 && std::fabs(areas[0] - 1.0) < 1e-12 && std::fabs(areas[1] - 60.0) < 1e-12;
            bool joined = areas.size() == 1 && std::fabs(areas[0] - 59.0) < 1e-12;
            if (s != SM_SUCCESS || (!apart && !joined))
            {
                std::fprintf(stderr, "loops touching at a vertex, variant %d, turn %zu: status %ld, %zu loops, largest area %g\n", variant, turn,
                             static_cast<long>(s), areas.size(), areas.empty() ? 0.0 : areas.back());
                failures++;
            }
        }
    }
}

// Two unit squares touching at a corner, (1, 1), from outside each other, both counter-clockwise
// or both clockwise, ordered into loops from the first square's lines then the second's and the
// other way round, each list from every rotation. Four curve ends meet at (1, 1), and the loop
// arriving there turns about one normal for the whole list: each square closes as itself, the
// way its lines run (areas 1 and 1, or -1 and -1), never one loop through the corner twice
// round both (area 2). Before, the normal there was the first pair of rays found by the ends'
// order in the list, so its sign, and which of the two the list gave, changed with the order.
static void loops_touching_from_outside_close_apart()
{
    const double a[4][2] = {{0, 0}, {1, 0}, {1, 1}, {0, 1}}, b[4][2] = {{1, 1}, {2, 1}, {2, 2}, {1, 2}};
    for (bool reversed : {false, true})
    {
        for (bool b_first : {false, true})
        {
            std::vector<std::pair<SmPoint3d, SmPoint3d>> lines;
            polygon_lines(lines, b_first ? b : a, 4, reversed);
            polygon_lines(lines, b_first ? a : b, 4, reversed);
            for (size_t turn = 0; turn < lines.size(); ++turn)
            {
                std::vector<double> areas;
                SmStatus s = ordered_loop_areas(rotated(lines, turn), areas);
                double sense = reversed ? -1.0 : 1.0;
                bool apart = areas.size() == 2 && std::fabs(areas[0] - sense) < 1e-12 && std::fabs(areas[1] - sense) < 1e-12;
                if (s != SM_SUCCESS || !apart)
                {
                    std::fprintf(stderr, "loops touching from outside, %s, %s first, turn %zu: status %ld, %zu loops, areas %g .. %g\n",
                                 reversed ? "clockwise" : "counter-clockwise", b_first ? "second" : "first", turn, static_cast<long>(s), areas.size(),
                                 areas.empty() ? 0.0 : areas.front(), areas.empty() ? 0.0 : areas.back());
                    failures++;
                }
            }
        }
    }
}

// A body from a SolidCraft flat file of capi/tests/data.
smlib_brep *read_flat(const char *path)
{
    smlib_flat *flat = nullptr;
    smlib_brep *b = nullptr;
    CHECK(smlib_flat_read(path, &flat) == SMLIB_OK && smlib_brep_import(flat, &b) == SMLIB_OK);
    smlib_flat_free(flat);
    return b;
}

// SolidCraft's sloped basin (solidcraft#1080): a 100 x 60 x 30 slab less the loft from the
// rectangle (20, m)-(80, 50) at the top z = 30 to (30, 20)-(70, 40) at z = 10, its walls planar
// B-spline faces, the front rim m wide (data: sloped-basin-1-mm-rim, -0-5-mm-rim, as SolidCraft
// hands them over). Its front rim edge chamfered d > m runs its leg on the top off the rim, over
// the slab's front face (a cliff): built as if the top went on and cut by the front face. At each
// end of the edge the pocket's corner is an Nx1 corner of mixed convexity (the end wall's top
// edge convex, the pocket's corner edge concave): it ends the chamfer on the end wall's surface
// extended past the vertex, along the end wall's top edge continued back over the rim, and that
// end face, in the sheet, has its edge along the pocket's corner edge on the surfaces of both
// walls. It was refused (SMLIB_ERR_NO_ROOM, "a cliff overflow's fillet boundary lies on the
// surfaces of two faces"). Removed: the section, the triangle of the edge (m, 30), the top's
// contact (m - d, 30) and the wall's d down it, clipped to y >= 0, along 30 + z between the
// sloped end walls (shoelace for its area and first moment).
void chamfer_past_a_sloped_basin_rim()
{
    auto removed = [](double m, double d)
    {
        const double len = std::sqrt((20.0 - m) * (20.0 - m) + 400.0);
        const double e[2] = {m, 30.0}, t[2] = {m - d, 30.0}, w[2] = {m + d * (20.0 - m) / len, 30.0 - d * 20.0 / len};
        std::vector<std::array<double, 2>> poly;
        if (t[0] >= 0.0)
        {
            poly = {{e[0], e[1]}, {t[0], t[1]}, {w[0], w[1]}};
        }
        else
        {
            const double f = -t[0] / (w[0] - t[0]);
            poly = {{e[0], e[1]}, {0.0, 30.0}, {0.0, t[1] + f * (w[1] - t[1])}, {w[0], w[1]}};
        }
        double area = 0.0, moment = 0.0;
        for (size_t i = 0; i < poly.size(); ++i)
        {
            const auto &a = poly[i];
            const auto &b = poly[(i + 1) % poly.size()];
            const double c = a[0] * b[1] - b[0] * a[1];
            area += c / 2.0;
            moment += (a[1] + b[1]) * c / 6.0;
        }
        return std::fabs(area) * 30.0 + std::fabs(moment);
    };
    const struct
    {
        const char *file;
        double m, d;
    } cases[3] = {{SMLIB_TEST_DATA "/sloped-basin-1-mm-rim.txt", 1.0, 2.0},
                  {SMLIB_TEST_DATA "/sloped-basin-1-mm-rim.txt", 1.0, 3.0},
                  {SMLIB_TEST_DATA "/sloped-basin-0-5-mm-rim.txt", 0.5, 4.0}};
    for (const auto &c : cases)
    {
        smlib_brep *basin = read_flat(c.file), *r = nullptr;
        const double rim[3] = {50.0, c.m, 30.0};
        const int32_t s = basin != nullptr ? smlib_blend(basin, SMLIB_BLEND_CHAMFER, rim, &c.d, 1, 0.1, &r, nullptr) : -1;
        if (s != SMLIB_OK || r == nullptr || !near(volume_of(r), volume_of(basin) - removed(c.m, c.d), 1e-8) || !sound(r))
        {
            std::fprintf(stderr, "  chamfer %g past the sloped basin's %g mm rim: status %d\n", c.d, c.m, static_cast<int>(s));
            failures++;
        }
        smlib_brep_free(r);
        smlib_brep_free(basin);
    }
}

// A loop of edges every one of which is blended, past its rim: the chamfers' boundary on the
// face they surround closes on that face, round the loop, and the region's face there lies
// between the two, one loop inside the other. It was refused (SMLIB_ERR_NO_ROOM, "a cliff
// overflow's fillet boundary closes on one face").
// - SolidCraft's equal-slope basin (solidcraft#1080, data: equal-slope-basin-1-mm-rim, the
//   loft from (20, 1)-(80, 50) at z = 30 to (30, 11)-(70, 40) at z = 10 out of a 100 x 60 x 30
//   slab, every wall sloping 10 over 20), its rim loop chamfered 2: the front leg runs off the
//   1 mm rim over the slab's front face. Sliced at depth t each chamfer is a strip
//   w(t) = d (1 - t / delta) wide beside its wall, delta = d / sqrt(1 + k^2), k = 1/2; mitred
//   at the rim's corners, the section removed is the pocket's rectangle at t grown by w(t) on
//   each side and clipped to y >= 0, less the rectangle: quadratic in t between delta and where
//   the front strip crosses y = 0, so Simpson's rule is exact on each piece.
// - A 60 x 40 x 10 plate less a 30 x 20 pocket 5 deep, its rim loop chamfered 5.5: the legs on
//   the walls run 0.5 past the floor. Removed: 50 d^2 + (4/3) d^3 (the pocket's perimeter 100
//   along the triangle d^2 / 2, and the four mitred corners), the groove round the floor
//   included.
void chamfer_of_a_rim_loop_past_its_rim()
{
    {
        smlib_brep *basin = read_flat(SMLIB_TEST_DATA "/equal-slope-basin-1-mm-rim.txt"), *r = nullptr;
        const double d = 2.0, k = 0.5, delta = d / std::sqrt(1.0 + k * k);
        auto area = [&](double t)
        {
            const double w = d * (1.0 - t / delta);
            const double x0 = 20.0 + k * t, x1 = 80.0 - k * t, y0 = 1.0 + k * t, y1 = 50.0 - k * t;
            return (x1 - x0 + 2.0 * w) * (y1 + w - std::max(y0 - w, 0.0)) - (x1 - x0) * (y1 - y0);
        };
        const double crossing = (d - 1.0) / (k + d / delta);
        double removed = 0.0;
        for (const auto &piece : {std::array<double, 2>{0.0, crossing}, std::array<double, 2>{crossing, delta}})
        {
            const double a = piece[0], b = piece[1];
            removed += (b - a) / 6.0 * (area(a) + 4.0 * area(0.5 * (a + b)) + area(b));
        }
        const double rim[12] = {50, 1, 30, 80, 25.5, 30, 50, 50, 30, 20, 25.5, 30};
        const double sizes[4] = {d, d, d, d};
        CHECK(basin != nullptr && near(volume_of(basin), 140333.333333333, 1e-10));
        CHECK(basin != nullptr && smlib_blend(basin, SMLIB_BLEND_CHAMFER, rim, sizes, 4, 0.1, &r, nullptr) == SMLIB_OK);
        if (r != nullptr)
        {
            CHECK(near(volume_of(r), volume_of(basin) - removed, 1e-8));
            done(r);
        }
        smlib_brep_free(basin);
    }
    {
        smlib_brep *plate = box(0, 0, 0, 60, 40, 10), *cut = box(15, 10, 5, 30, 20, 6), *pocket = nullptr, *r = nullptr;
        CHECK(smlib_boolean(plate, cut, SMLIB_BOOLEAN_SUBTRACT, &pocket) == SMLIB_OK);
        const double rim[12] = {30, 10, 10, 30, 30, 10, 15, 20, 10, 45, 20, 10};
        const double d = 5.5, sizes[4] = {d, d, d, d};
        CHECK(pocket != nullptr && smlib_blend(pocket, SMLIB_BLEND_CHAMFER, rim, sizes, 4, 0.1, &r, nullptr) == SMLIB_OK);
        if (r != nullptr)
        {
            CHECK(near(volume_of(r), 21000.0 - (50.0 * d * d + 4.0 / 3.0 * d * d * d), 1e-8));
            done(r);
        }
        for (smlib_brep *b : {plate, cut, pocket})
        {
            smlib_brep_free(b);
        }
    }
}

// Two variable rounds, 0.5 at their far ends to 2 at the corner (10, 10, 10) of a box of 10 they
// meet at: the top edges along y at x = 10 and along x at y = 10 (the third edge there is not
// rounded). The removed volume, or -1 when the rounds fail or the body is unsound.
static double variable_rounds_meeting_volume()
{
    SmApiGetOrCreateContext();
    SmBrep *p = nullptr;
    CHECK(SmApiCreateBox(SmVector3d(0, 0, 0), 10, 10, 10, p) == SM_SUCCESS && p != nullptr);
    smlib_brep *b = new smlib_brep{p};
    SmTArray<SmEdge *> all, edges;
    p->GetEdges(all);
    const SmPoint3d corner(10, 10, 10), far0(10, 0, 10), far1(0, 10, 10);
    for (ULONG i = 0; i < all.GetSize(); ++i)
    {
        const SmPoint3d a = all[i]->GetStartVertex()->GetPoint(), e = all[i]->GetEndVertex()->GetPoint();
        if (e.DistanceBetween(corner) < 1e-9 && (a.DistanceBetween(far0) < 1e-9 || a.DistanceBetween(far1) < 1e-9))
        {
            edges.Add(all[i]);
        }
    }
    CHECK(edges.GetSize() == 2);
    double v = -1.0;
    if (edges.GetSize() == 2 && SmApiVariableRadiusFillet(p, edges, 0.5, 2.0, 1, 1) == SM_SUCCESS && sound(b))
    {
        v = 1000.0 - volume_of(b);
    }
    smlib_brep_free(b);
    return v;
}

// The three edges of an s x s x s box at the origin rounded r through SmApiVariableRadiusFillet
// (constant law, circular sections, G1), its corner a Coons patch; the volume removed, or -1
// when the round fails, the body is unsound or it has not the 10 faces (6, 3 rounds, the corner).
static double variable_corner_removed(double s, double r)
{
    SmApiGetOrCreateContext();
    SmBrep *p = nullptr;
    CHECK(SmApiCreateBox(SmVector3d(0, 0, 0), s, s, s, p) == SM_SUCCESS && p != nullptr);
    smlib_brep *b = new smlib_brep{p};
    SmTArray<SmVertex *> vertices;
    p->GetVertices(vertices);
    SmTArray<SmEdge *> edges;
    for (ULONG i = 0; i < vertices.GetSize(); ++i)
    {
        if (vertices[i]->GetPoint().Length() < 1e-9 * s)
        {
            vertices[i]->GetEdges(edges);
        }
    }
    CHECK(edges.GetSize() == 3);
    double v = -1.0;
    if (edges.GetSize() == 3 && SmApiVariableRadiusFillet(p, edges, r, r, 1, 1) == SM_SUCCESS && sound(b) && p->GetNumFaces() == 10)
    {
        v = s * s * s - volume_of(b);
    }
    smlib_brep_free(b);
    return v;
}

// Rounds of a sloped basin's rim (solidcraft#1080, the rig's fillets/1080-fillet-* cases, data:
// sloped-basin-1-mm-rim, -0-5-mm-rim and equal-slope-basin-1-mm-rim) rounded R4: the front round
// runs its rail on the top off the thin front rim, over the slab's front face (a cliff), and goes
// first; the other rounds fit, and each meets the front one, or its other neighbours, in a mitre
// at the pocket's concave corner. Mitred rounds are one corner, built together whichever
// overflows. They were split (the front round alone, then the others on what it left, where the
// right round's rail met the front round's cut, a concave edge, beyond its wall: SMLIB_ERR_NO_ROOM,
// "a round's rail leaves the face it rolls on into the body's material"). The volumes are the
// rig's analytic answers: each round's section at height z, from its wall to its arc and not
// past y = 0, along its edge between its neighbours' walls (the sloped end walls extended),
// Int width(z) length(z) dz, and at each mitred corner the quadrant beyond both walls,
// Int width_1(z) width_2(z) dz. Within 1e-8 of the body: the mitre's curve is SMLib's fitted
// fillet-fillet intersection.
void rounds_mitred_to_one_past_a_rim()
{
    const struct
    {
        const char *name, *file;
        double m;
        int edges;
        double volume;
    } cases[4] = {
        {"1 mm rim, front and right", SMLIB_TEST_DATA "/sloped-basin-1-mm-rim.txt", 1.0, 2, 144461.19634789118},
        {"1 mm rim, all four", SMLIB_TEST_DATA "/sloped-basin-1-mm-rim.txt", 1.0, 4, 144347.32581768231},
        {"0.5 mm rim, all four", SMLIB_TEST_DATA "/sloped-basin-0-5-mm-rim.txt", 0.5, 4, 144084.31810540144},
        {"equal slope, all four", SMLIB_TEST_DATA "/equal-slope-basin-1-mm-rim.txt", 1.0, 4, 140114.73514508773},
    };
    for (const auto &c : cases)
    {
        smlib_brep *basin = read_flat(c.file), *r = nullptr;
        const double mid = 0.5 * (c.m + 50.0);
        const double rim[12] = {50, c.m, 30, 80, mid, 30, 50, 50, 30, 20, mid, 30};
        const double sizes[4] = {4, 4, 4, 4};
        int64_t dropped[4] = {-1, -1, -1, -1};
        const int32_t s = basin != nullptr ? smlib_blend_dropping(basin, SMLIB_BLEND_ROUND, rim, sizes, c.edges, 0.5, &r, nullptr, dropped) : -1;
        if (s != SMLIB_OK || r == nullptr || !near(volume_of(r), c.volume, 1e-8) || !sound(r))
        {
            std::fprintf(stderr, "  R4 round of the basin's rim, %s: status %d\n", c.name, static_cast<int>(s));
            failures++;
        }
        for (int k = 0; k < c.edges; ++k)
        {
            CHECK(dropped[k] == 0);
        }
        smlib_brep_free(r);
        smlib_brep_free(basin);
    }
}

// SolidCraft's rounded rectangle (solidcraft#2498: a 26.4 x 8.6 outline, its corners R4.29,
// extruded 8.4), its top rim rounded R9: too large for the 8.6 top, refused cleanly. The big-radius
// conversion (SmFilletExecutive::TestAndConvertBigRadFillet) deleted a fillet vertex its corner
// still listed; where the conversion then failed the corner deleted it again (SIGBUS).
void big_round_of_a_rounded_rectangle_is_refused()
{
    const double w = 26.4, h = 8.6, r = 4.29;
    std::vector<int64_t> kinds;
    std::vector<double> reals;
    auto line = [&](double x0, double y0, double x1, double y1)
    {
        kinds.push_back(SMLIB_CURVE_LINE);
        reals.insert(reals.end(), {10 + x0, 10 + y0, 0, 10 + x1, 10 + y1, 0});
    };
    auto arc = [&](double cx, double cy, double from, double to)
    {
        kinds.push_back(SMLIB_CURVE_ARC);
        reals.insert(reals.end(), {10 + cx, 10 + cy, 0, 0, 0, 1, 1, 0, 0, r, from, to});
    };
    line(r, 0, w - r, 0);
    arc(w - r, r, -0.5 * PI, 0);
    line(w, r, w, h - r);
    arc(w - r, h - r, 0, 0.5 * PI);
    line(w - r, h, r, h);
    arc(r, h - r, 0.5 * PI, PI);
    line(0, h - r, 0, r);
    arc(r, r, PI, 1.5 * PI);
    const double up[3] = {0, 0, 1};
    double corner[3];
    smlib_brep *bar = nullptr, *out = nullptr;
    CHECK(smlib_extrude(kinds.data(), static_cast<int64_t>(kinds.size()), reals.data(), static_cast<int64_t>(reals.size()), up, 8.4, 0, 0, 1,
                        &bar, corner) == SMLIB_OK);
    const double rim[24] = {23.2, 10.0, 8.4, 35.14348809129028, 11.25651190870971, 8.4, 11.256511908709715, 11.256511908709712, 8.4,
                            10.0, 14.3, 8.4, 11.256511908709713, 17.343488091290293, 8.4, 23.2, 18.6, 8.4,
                            35.14348809129028, 17.34348809129029, 8.4, 36.4, 14.3, 8.4};
    const double sizes[8] = {9, 9, 9, 9, 9, 9, 9, 9};
    CHECK(bar != nullptr && smlib_blend(bar, SMLIB_BLEND_ROUND, rim, sizes, 8, 0.5, &out, nullptr) != SMLIB_OK && out == nullptr);
    smlib_brep_free(out);
    smlib_brep_free(bar);
}

// Simpson's rule for f on [a, b] in n (even) spans.
template <class F> double simpson(F f, double a, double b, int n)
{
    const double h = (b - a) / n;
    double sum = f(a) + f(b);
    for (int i = 1; i < n; ++i)
    {
        sum += (i % 2 ? 4.0 : 2.0) * f(a + i * h);
    }
    return sum * h / 3.0;
}

// A closed rim blended past the far end of its wall (the rig's fillets/190-fillet-disc-rim-past-
// its-2mm-wall, every-hole-top-rim-r12-past-plate-bottom, 2893-chamfer-two-distances-past-the-
// bottom; the owner's ruling: built as if the wall went on, then cut by the face beyond). The
// wall's rail runs round the wall's cylinder continued past its end; the closed corner meets the
// seam's line there, and the region the blend removes closes on the wall's cylinder continued
// as itself, the band from the rim to the rail. Before: the closed corner found no point on the
// seam (1001, "Unable to Compute Geometry"), and past it the region found no exact continuation
// of the wall (SMLIB_ERR_NOT_EXACT) nor a face closing round a closed edge alone.
// - A disc R13 x 2, its top rim rounded R2.5 (centre at rho 10.5, z -0.5): pi 338 less
//   pi Int (169 - (10.5 + 2.5 cos t)^2) 2.5 cos t dt, t from asin(0.2) to pi/2. The wall is used
//   up: top, round and bottom.
// - A 40 x 40 x 10 plate with a through hole R5, its top rim rounded R12 (the circle's centre at
//   rho 17, z -2): the spandrel clipped to the plate, at t = rho - 5 down min(12 - sqrt(144 -
//   (12 - t)^2), 10), 2 pi Int (5 + t) depth dt, the bottom face trimmed round the hole to
//   rho 12 - sqrt 140 + 5. Seven faces.
// - A cylinder R10 x 20, its top rim chamfered 2 across the top and 25 down the wall, 5 past the
//   bottom: the frustum from R8 at the top to R9.6 at the bottom. Three faces.
void closed_rims_blended_past_their_walls()
{
    auto faces = [](const smlib_brep *b)
    {
        int64_t c[3] = {0, 0, 0};
        CHECK(smlib_brep_counts(b, c) == SMLIB_OK);
        return c[0];
    };
    {
        const double o[3] = {0, 0, 0};
        smlib_brep *disc = nullptr, *r = nullptr;
        CHECK(smlib_make_cylinder(o, 13.0, 2.0, &disc) == SMLIB_OK);
        const double taken = PI * simpson([](double t) { const double rho = 10.5 + 2.5 * std::cos(t); return (169.0 - rho * rho) * 2.5 * std::cos(t); },
                                          std::asin(0.2), 0.5 * PI, 2000);
        const double rim[3] = {13.0 * std::cos(1.0), 13.0 * std::sin(1.0), 2.0}, radius = 2.5;
        CHECK(disc != nullptr && smlib_blend(disc, SMLIB_BLEND_ROUND, rim, &radius, 1, 0.5, &r, nullptr) == SMLIB_OK);
        if (r != nullptr)
        {
            CHECK(near(volume_of(r), PI * 338.0 - taken, 1e-9));
            CHECK(faces(r) == 3);
            done(r);
        }
        smlib_brep_free(disc);
    }
    {
        smlib_brep *plate = box(0, 0, 0, 40, 40, 10), *pin = nullptr, *holed = nullptr, *r = nullptr;
        const double o[3] = {20, 20, -1};
        CHECK(smlib_make_cylinder(o, 5.0, 12.0, &pin) == SMLIB_OK);
        CHECK(smlib_boolean(plate, pin, SMLIB_BOOLEAN_SUBTRACT, &holed) == SMLIB_OK);
        auto ring = [](double t) { return 2.0 * PI * (5.0 + t) * std::min(12.0 - std::sqrt(144.0 - (12.0 - t) * (12.0 - t)), 10.0); };
        const double kink = 12.0 - std::sqrt(140.0);
        const double taken = simpson(ring, 0.0, kink, 200) + simpson(ring, kink, 12.0, 20000);
        const double rim[3] = {20.0 + 5.0 * std::cos(1.0), 20.0 + 5.0 * std::sin(1.0), 10.0}, radius = 12.0;
        CHECK(holed != nullptr && smlib_blend(holed, SMLIB_BLEND_ROUND, rim, &radius, 1, 0.5, &r, nullptr) == SMLIB_OK);
        if (r != nullptr)
        {
            CHECK(near(volume_of(r), 16000.0 - 250.0 * PI - taken, 1e-9));
            CHECK(faces(r) == 7);
            done(r);
        }
        for (smlib_brep *b : {plate, pin, holed})
        {
            smlib_brep_free(b);
        }
    }
    {
        const double o[3] = {0, 0, 0}, up[3] = {0, 0, 1};
        smlib_brep *can = nullptr, *r = nullptr;
        CHECK(smlib_make_cylinder(o, 10.0, 20.0, &can) == SMLIB_OK);
        const double rim[3] = {10.0 * std::cos(1.0), 10.0 * std::sin(1.0), 20.0};
        CHECK(can != nullptr && smlib_chamfer_sides(can, rim, 1, 0.5, 2.0, 25.0, 0, up, &r, nullptr) == SMLIB_OK);
        if (r != nullptr)
        {
            CHECK(near(volume_of(r), PI * 20.0 * (64.0 + 76.8 + 92.16) / 3.0, 1e-9));
            CHECK(faces(r) == 3);
            done(r);
        }
        smlib_brep_free(can);
    }
}

// A surface of revolution whose generator is a circle carries the circle on past its domain
// (sm_ExtendSurface, as a torus does): SMLib's own concave r2 round of the boss's base arc in
// bar-boss-over-sides (the bar 100 x 30 x 10 and the boss of diameter 40 at the origin) is the
// round's quarter circle revolved about the boss's axis, a tube of radius 2 about the circle of
// radius 22 at z = 12. Its extension along the generator was its tangent continuation, and a
// round of the bar's edge chained onto the round's end (SolidCraft rig fillets/1207 r0.5) traced
// its rail past the round's rim on it: the rail stood still there, and started its span on the
// round with a tangent 1.2 degrees off. Offset 0.5, every point of the extension past either end
// of the generator lies on the tube of radius 2 +- 0.5 (the side the surface's normal says; past
// the rim the tube runs on inside the boss), and moves with the parameter.
void revolved_circle_extends_as_its_circle()
{
    smlib_flat *flat = nullptr;
    smlib_brep *joined = nullptr, *r = nullptr;
    CHECK(smlib_flat_read(SMLIB_TEST_DATA "/bar-boss-over-sides.txt", &flat) == SMLIB_OK && smlib_brep_import(flat, &joined) == SMLIB_OK);
    smlib_flat_free(flat);
    const double arc[6] = {20 * std::cos(0.73), 20 * std::sin(0.73), 10, 20 * std::cos(-0.12), 20 * std::sin(-0.12), 10};
    const double sizes[2] = {2, 2};
    CHECK(joined != nullptr && smlib_blend(joined, SMLIB_BLEND_ROUND, arc, sizes, 2, 0.1, &r, nullptr) == SMLIB_OK && r != nullptr);
    if (r == nullptr)
    {
        smlib_brep_free(joined);
        return;
    }
    SmTArray<SmFace *> faces;
    r->p->GetFaces(faces);
    int tested = 0;
    double worst = 0.0, least_step = 1e300;
    for (ULONG i = 0; i < faces.GetSize(); ++i)
    {
        SmSurface *s = faces[i]->GetSurface();
        if (s == nullptr || s->GetType() != SmSurfOfRevolution_TYPE)
        {
            continue;
        }
        ++tested;
        SmOffsetSurface *off = new (*s->GetContext()) SmOffsetSurface(0.5, *s);
        off->SetAllowExtension(TRUE);
        const SmExtent2d d = s->GetNaturalUVDomain();
        // the generator: the direction other than the rotation round the axis
        SmSurfParamType rot = SM_SP_U;
        (void)static_cast<SmSurfOfRevolution *>(s)->GetRotateParamExtent(&rot);
        const int g = rot == SM_SP_U ? 1 : 0;
        const double lo[2] = {d.GetMin().x, d.GetMin().y}, size[2] = {d.GetSize().x, d.GetSize().y};
        for (double past : {0.1, 0.3, 0.6})
        {
            for (int side : {-1, 1})
            {
                SmPoint3d p, q;
                double at[2] = {lo[0] + 0.5 * size[0], lo[1] + 0.5 * size[1]};
                at[g] = side < 0 ? lo[g] - past * size[g] : lo[g] + (1.0 + past) * size[g];
                double at2[2] = {at[0], at[1]};
                at2[g] += side * 0.1 * size[g];
                CHECK(off->EvaluatePoint(SmPoint2d(at[0], at[1]), p) == SM_SUCCESS && off->EvaluatePoint(SmPoint2d(at2[0], at2[1]), q) == SM_SUCCESS);
                const double rho = std::hypot(p.x, p.y);
                const double tube = std::hypot(rho - 22.0, p.z - 12.0);
                worst = std::max(worst, std::min(std::fabs(tube - 2.5), std::fabs(tube - 1.5)));
                least_step = std::min(least_step, p.DistanceBetween(q));
            }
        }
        delete off;
    }
    std::printf("r2 round of the boss's arc offset 0.5, past its generator's ends: off the tube by at most %.3g, steps at least %.3g\n",
                worst, least_step);
    CHECK(tested == 2 && worst < 1e-9 && least_step > 1e-3);
    smlib_brep_free(r);
    smlib_brep_free(joined);
}

// SolidCraft rig case fillets/926: the top edge x = 40, z = 20 of a 40 x 40 x 20 box joined with
// an R12 boss 40 high centred on its side (40, 20), chamfered 0.8. The edge runs from the front
// face to the boss's wall, which rises above the box. The chamfer's rail on the top face crosses
// the boss's base arc just past the edge's end (y = 8.0267 against 8) and comes back onto the
// top face beyond the boss: a run off the face between two on it, which was taken for a rail
// leaving its face beside the edge, so the chamfer was built as a cliff overflow and refused
// (SMLIB_ERR_NO_ROOM, "a cliff overflow's fillet boundary lies on the surfaces of two faces").
// V0 = 32000 + 144 pi 40 - 144 pi 20 / 2; the chamfer takes the section h = x - 39.2 swept
// along y to the wall at y = 20 - sqrt(144 - (x - 40)^2).
void chamfer_ending_on_a_boss_wall_above_it()
{
    smlib_brep *b = box(0, 0, 0, 40, 40, 20), *boss = cylinder(40, 20, 0, 12, 40), *joined = nullptr, *r = nullptr;
    CHECK(smlib_boolean(b, boss, SMLIB_BOOLEAN_UNION, &joined) == SMLIB_OK && joined != nullptr);
    const int n = 20000;
    double removed = 0.0;
    for (int i = 0; i < n; ++i)
    {
        const double x = 39.2 + (i + 0.5) * 0.8 / n;
        removed += (x - 39.2) * (20.0 - std::sqrt(144.0 - (x - 40.0) * (x - 40.0))) * 0.8 / n;
    }
    const double want = 32000.0 + PI * 144.0 * 40.0 - PI * 72.0 * 20.0 - removed;
    const double edge[3] = {40, 4, 20};
    const double d = 0.8;
    const int32_t s = joined != nullptr ? smlib_blend(joined, SMLIB_BLEND_CHAMFER, edge, &d, 1, 0.1, &r, nullptr) : -1;
    if (s != SMLIB_OK || r == nullptr || !near(volume_of(r), want, 1e-9) || !sound(r))
    {
        std::fprintf(stderr, "  chamfer 0.8 of the box's top edge ending on the boss's wall: status %d\n", static_cast<int>(s));
        failures++;
    }
    for (smlib_brep *x : {b, boss, joined, r})
    {
        smlib_brep_free(x);
    }
}

// Rails trimmed at an end vertex within the zone tolerance of one of their own knots (marching
// points): the trim left an end span shorter than the tolerance, its control points coincident,
// which SMLib's check refuses as an internal point singularity (Test #15) of the rail and of the
// face curves trimmed with it. As SolidCraft hands them over:
// - rig case fillets/1207 r0.5 (bar-boss-arc-rounded.smb: the bar and boss of bar-boss-over-sides
//   with the boss's base arc over the bar rounded r2): the bar's front top edge and the r2
//   round's end curve on the front face, its tangent chain, rounded r0.5. The rail on the r2
//   round stands still at the round's boundary past its end (where it meets the boss), and its
//   vertex dropped onto it 1e-8 before the knot where it starts to move.
// - rig case fillets/cov-fillet-arc-edge-on-an-oblique-cut-of-a-bar (oblique-cut-bar: a bar with
//   a concave R3 arc in its profile cut by a plane tilted 21 degrees): the elliptical edge where
//   the arc's cylinder meets the cut rounded r0.1. One rail crosses its side edge 6e-7 mm before
//   the other, where the march put its point.
// Each is one closed valid solid.
void rails_trimmed_next_to_a_knot()
{
    const struct
    {
        const char *file;
        std::vector<double> points;
        double r, pick;
        // Still refused by SMLib's check for another cause: the end rib on the boss, where the
        // round and the boss are tangent, is marched into its tangential start in spans too short
        // for the tolerance. Listed so that this test fails once it builds, and the mark comes off.
        bool known_failing;
    } cases[2] = {{SMLIB_TEST_DATA "/bar-boss-arc-rounded.smb", {58.046738469715542, -15, 10, 14.289893382452973, -14.999999999999996, 10.475338853397696}, 0.5, 0.1, true},
                  {SMLIB_TEST_DATA "/oblique-cut-bar.txt", {149.60309403593678, -458.87269130746967, 848.84285884468318}, 0.1, 1.0, false}};
    for (const auto &c : cases)
    {
        // The 1207 body is SMLib's own (a native file, as SolidCraft's dump wrote it): through flat
        // arrays its rail drops onto the stationary stretch elsewhere, clear of the knot.
        smlib_brep *b = nullptr;
        if (std::strstr(c.file, ".smb") != nullptr)
        {
            CHECK(smlib_brep_read(c.file, 1, &b) == SMLIB_OK);
        }
        else
        {
            b = read_flat(c.file);
        }
        smlib_brep *r = nullptr;
        const int64_t k = static_cast<int64_t>(c.points.size() / 3);
        const std::vector<double> sizes(static_cast<size_t>(k), c.r);
        const int32_t s = b != nullptr ? smlib_blend(b, SMLIB_BLEND_ROUND, c.points.data(), sizes.data(), k, c.pick, &r, nullptr) : -1;
        const bool built = s == SMLIB_OK && r != nullptr && sound(r);
        if (built == c.known_failing)
        {
            std::fprintf(stderr, built ? "  round %g of %s builds now: take it off the known failing\n" : "  round %g of %s: status %d\n", c.r, c.file, static_cast<int>(s));
            failures++;
        }
        smlib_brep_free(r);
        smlib_brep_free(b);
    }
}

// A vertex where the body touches itself split into one vertex per fan of faces round it, and
// joined again (SmBrep::SplitVertexByFans, JoinVertexFans): the up-down squares (data:
// up-down-squares, as in concave_chamfers_at_a_vertex_the_body_touches_itself_at) touch
// themselves at the origin, where the faces round the two missing octants make two fans of
// three edges each. Split, the origin is two vertices of three edges each and the body one
// closed valid solid with one vertex more; joined, it is the body it was. A box's corner has
// one fan and stays as it is.
void vertex_fans_split_and_joined()
{
    smlib_brep *b = read_flat(SMLIB_TEST_DATA "/up-down-squares.txt");
    if (b == nullptr)
    {
        CHECK(!"the up-down squares read");
        return;
    }
    auto at = [](SmBrep *brep, const SmPoint3d &p) {
        std::vector<SmVertex *> out;
        SmTArray<SmVertex *> vs;
        brep->GetVertices(vs);
        for (ULONG i = 0; i < vs.GetSize(); ++i)
        {
            if (vs[i]->GetPoint().DistanceBetween(p) < 1e-9)
            {
                out.push_back(vs[i]);
            }
        }
        return out;
    };
    auto edges = [](SmVertex *v) {
        SmTArray<SmEdge *> es;
        v->GetEdges(es);
        return es.GetSize();
    };
    const ULONG before = b->p->GetNumVertices();
    const double volume = volume_of(b);
    std::vector<SmVertex *> origin = at(b->p, SmPoint3d(0, 0, 0));
    CHECK(origin.size() == 1 && edges(origin[0]) == 6);
    SmTArray<SmVertex *> made;
    CHECK(!origin.empty() && b->p->SplitVertexByFans(origin[0], made) == SM_SUCCESS);
    CHECK(made.GetSize() == 1 && b->p->GetNumVertices() == before + 1);
    origin = at(b->p, SmPoint3d(0, 0, 0));
    CHECK(origin.size() == 2 && edges(origin[0]) == 3 && edges(origin[1]) == 3);
    CHECK(b->p->IsManifoldSolid() && sound(b) && near(volume_of(b), volume, 1e-12));
    CHECK(origin.size() == 2 && b->p->JoinVertexFans(origin[0], origin[1]) == SM_SUCCESS);
    origin = at(b->p, SmPoint3d(0, 0, 0));
    CHECK(origin.size() == 1 && edges(origin[0]) == 6 && b->p->GetNumVertices() == before);
    CHECK(sound(b) && near(volume_of(b), volume, 1e-12));
    smlib_brep_free(b);

    smlib_brep *cube = box(0, 0, 0, 10, 10, 10);
    std::vector<SmVertex *> corner = at(cube->p, SmPoint3d(10, 10, 10));
    made.RemoveAll();
    CHECK(corner.size() == 1 && cube->p->SplitVertexByFans(corner[0], made) == SM_SUCCESS && made.GetSize() == 0 && cube->p->GetNumVertices() == 8);
    smlib_brep_free(cube);
}

// SolidCraft rig case fillets/cov-chamfer-concave-edge-through-up-down-joined-squares
// (up-down-squares): four 10 x 10 squares drawn edge to edge, three extruded 10 up and three 10
// down and joined. The down part lacks the square x, y < 0 and the up part x, y > 0, so the line
// x = y = 0 is a concave edge of each, and the origin a vertex where the body touches itself:
// two cones of faces, one round each missing octant. The two concave edges chamfered 1 fill
// their corners with the triangle 1 x 1 / 2 along 10, on opposite sides of the line:
// V = 6000 + 10. Each chamfer ends on its own cone of faces at the origin; one corner of both
// was 1001 (SmFilletNx2Corner: unable to make topology).
void concave_chamfers_at_a_vertex_the_body_touches_itself_at()
{
    smlib_brep *b = read_flat(SMLIB_TEST_DATA "/up-down-squares.txt"), *r = nullptr;
    const double edges[6] = {0, 0, -5, 0, 0, 5};
    const double d[2] = {1, 1};
    const int32_t s = b != nullptr ? smlib_blend(b, SMLIB_BLEND_CHAMFER, edges, d, 2, 0.1, &r, nullptr) : -1;
    if (s != SMLIB_OK || r == nullptr || !near(volume_of(r), 6010.0, 1e-9) || !sound(r))
    {
        std::fprintf(stderr, "  the up and down parts' concave edges chamfered 1: status %d\n", static_cast<int>(s));
        failures++;
    }
    smlib_brep_free(r);
    smlib_brep_free(b);
}

// SolidCraft rig case fillets/cov-chamfer-outer-corner-edge-of-a-shelled-box-with-a-lip
// (shelled-box-lip): an 80 x 50 x 30 box shelled 2 open at the top, a 1 x 2 lip on the inner half
// of its rim. Its outer corner edge x = 0, y = 50 chamfered 3: the chamfer's plane
// x + (50 - y) = 3 clears the inside corner (2, 48), and its end section at the top crosses the
// lip's foot. The chamfer face came back unsewn: not a closed solid (SMLIB_ERR_NOT_SOLID).
// The section's triangle (0, 50), (3, 50), (0, 47) lies in the walls all the way up, 4.5 x 30 =
// 135 taken off V0 = 22608. Where it ends at the top is a ruling (one line, below): the chamfer
// cut off by its end face's plane z = 30, the lip left over it (a 0.5 triangle of its underside
// open there); or run on up through the lip to its top, taking that triangle's 0.5 x 2 too.
void chamfer_of_a_shelled_box_corner_ending_under_its_lip()
{
    const bool through_the_lip = false; // the owner's ruling: false = cut off at the end face's plane
    const double want = 22608.0 - 135.0 - (through_the_lip ? 1.0 : 0.0);
    smlib_brep *b = read_flat(SMLIB_TEST_DATA "/shelled-box-lip.txt"), *r = nullptr;
    const double edge[3] = {0, 50, 15};
    const double d = 3;
    const int32_t s = b != nullptr ? smlib_blend(b, SMLIB_BLEND_CHAMFER, edge, &d, 1, 0.1, &r, nullptr) : -1;
    if (s != SMLIB_OK || r == nullptr || !near(volume_of(r), want, 1e-9) || !sound(r))
    {
        std::fprintf(stderr, "  the shelled box's corner edge chamfered 3 up to its lip: status %d\n", static_cast<int>(s));
        failures++;
    }
    smlib_brep_free(r);
    smlib_brep_free(b);
}

int main(int argc, char **argv)
{
    if (argc == 4 && std::strcmp(argv[1], "--solidcore-1173") == 0)
    {
        return corner_child(argv[2], argv[3]);
    }

    loops_touching_at_a_vertex_in_any_order();
    four_faces_meeting_at_a_vertex_shell();
    shell_open_beside_a_step_keeps_the_floor_under_the_opening();
    hole_whose_offset_leaves_the_round_side_shell();
    floor_thinner_than_two_walls_shell();
    up_down_joined_squares_shell();
    loops_touching_from_outside_close_apart();

    smlib_brep *cube = box(0, 0, 0, 20, 20, 20);
    smlib_brep *r = nullptr;

    // SMLib defect 1 (SmProtoEdge's destructor left a pointer to itself in a ProtoVertex that
    // bounds both its ends; ~SmProtoTopologyManager then wrote through it): a box minus a
    // through cylinder.
    // A stadium ring base (arcs r 3 and 5, centres 20 apart, 2 thick) joined with two half-annulus
    // lips (r 4..5, 18 tall) standing on it with their walls 0.0001 off the base's (SolidCraft rig
    // booleans/2928): the base's top arc is split 3.5e-6 of its parameter before the knot of its
    // UV curve, on a plane whose UV runs 38 times shorter than its 3D size. The split left that
    // short span of the UV curve; trimming used to drop it as degenerate (a 3D tolerance read in
    // UV units) and the pcurve lifted 2.5 times the edge tolerance off its edge.
    {
        const double e = 1e-4;
        Profile base, lip1, lip2;
        base.line(5, -10, 5, 10, 0).arc(0, 10, 0, 1, 5, 0, PI).line(-5, 10, -5, -10, 0).arc(0, -10, 0, 1, 5, PI, 2 * PI);
        base.line(3, -10, 3, 10, 0).arc(0, 10, 0, 1, 3, 0, PI).line(-3, 10, -3, -10, 0).arc(0, -10, 0, 1, 3, PI, 2 * PI);
        lip1.arc(e, 10, 2, 1, 5, 0, PI).line(-5 + e, 10, -4 + e, 10, 2).arc(e, 10, 2, -1, 4, -PI, 0).line(4 + e, 10, 5 + e, 10, 2);
        lip2.arc(-e, -10, 2, 1, 5, PI, 2 * PI).line(5 - e, -10, 4 - e, -10, 2).arc(-e, -10, 2, -1, 4, -2 * PI, -PI).line(-4 - e, -10, -5 - e, -10, 2);
        smlib_brep *b = base.extrude(2.0), *l1 = lip1.extrude(18.0), *l2 = lip2.extrude(18.0), *one = nullptr, *two = nullptr;
        CHECK(smlib_boolean(b, l1, SMLIB_BOOLEAN_UNION, &one) == SMLIB_OK);
        CHECK(one != nullptr && smlib_boolean(one, l2, SMLIB_BOOLEAN_UNION, &two) == SMLIB_OK);
        if (two != nullptr)
        {
            // Base 2 * (80 + 16 pi), lips 2 * 18 * 4.5 pi, not overlapping.
            CHECK(near(volume_of(two), 160.0 + 194.0 * PI, 1e-6));
            CHECK(sound(two));
            double worst = worst_pcurve_gap(two);
            std::printf("lips on a stadium ring: worst pcurve gap %.3g edge tolerances\n", worst);
            CHECK(worst <= 1.0);
        }
        for (smlib_brep *x : {b, l1, l2, one, two})
        {
            smlib_brep_free(x);
        }
    }

    // A box 80 x 79 x 32.5 and a flange plate 100 x 79 x 2.5 under it, both with their upright
    // corners rounded r 3, joined (SolidCraft boolean_regression step_export_of_a_flanged_enclosure).
    // The plate's long sides lie along the box's planes y = 0 and y = 79 and touch its corner
    // cylinders along their ends' lines. SolidCraft's plate came with its planes turned 1e-13 rad;
    // SmCone::IntersectWithPlane took a plane as tangent to a cylinder only when its two crossing
    // lines were within tolerance of each other, so the plane cut each corner cylinder in two lines
    // a few microns apart: a sliver face (turned 1e-11) or an unsound body (1e-9). Turned about z
    // through its centre, the plate stays within 1e-7 of its planes. Volume: the box's section
    // (80 x 79 less (4 - pi) 9) 30 high on the plate's (100 x 79 less the same) 2.5 high.
    {
        auto rounded_rect = [](double x0, double x1, double y1, double h)
        {
            const double r = 3.0;
            Profile p;
            p.line(x0 + r, 0, x1 - r, 0, 0).arc(x1 - r, r, 0, 1, r, -PI / 2, 0).line(x1, r, x1, y1 - r, 0).arc(x1 - r, y1 - r, 0, 1, r, 0, PI / 2);
            p.line(x1 - r, y1, x0 + r, y1, 0).arc(x0 + r, y1 - r, 0, 1, r, PI / 2, PI).line(x0, y1 - r, x0, r, 0).arc(x0 + r, r, 0, 1, r, PI, 1.5 * PI);
            return p.extrude(h);
        };
        smlib_brep *box_ = rounded_rect(0, 80, 79, 32.5), *plate = rounded_rect(-10, 90, 79, 2.5);
        const double corners = (4.0 - PI) * 9.0;
        const double want = (80.0 * 79.0 - corners) * 30.0 + (100.0 * 79.0 - corners) * 2.5;
        for (double turn : {0.0, 1e-13, 1e-11, 1e-9})
        {
            const double c = std::cos(turn), sn = std::sin(turn);
            const double m[4][4] = {{c, -sn, 0, 40 - 40 * c + 39.5 * sn}, {sn, c, 0, 39.5 - 40 * sn - 39.5 * c}, {0, 0, 1, 0}, {0, 0, 0, 1}};
            smlib_brep *turned = nullptr, *joined = nullptr;
            CHECK(smlib_brep_transform(plate, m, &turned) == SMLIB_OK);
            CHECK(turned != nullptr && smlib_boolean(box_, turned, SMLIB_BOOLEAN_UNION, &joined) == SMLIB_OK);
            if (joined != nullptr)
            {
                CHECK(near(volume_of(joined), want, 1e-9));
                // no sliver: every face spans more than 1e-3 in two directions
                SmTArray<SmFace *> faces;
                joined->p->GetFaces(faces);
                int slivers = 0;
                for (ULONG i = 0; i < faces.GetSize(); i++)
                {
                    SmExtent3d fb;
                    faces[i]->CalculateBoundingBox(fb);
                    SmVector3d d = fb.GetMax() - fb.GetMin();
                    slivers += (d.x > 1e-3) + (d.y > 1e-3) + (d.z > 1e-3) < 2 ? 1 : 0;
                }
                std::printf("flange plate turned %g: %lu faces, %d slivers\n", turn, static_cast<unsigned long>(faces.GetSize()), slivers);
                CHECK(slivers == 0);
                done(joined);
            }
            smlib_brep_free(turned);
        }
        smlib_brep_free(box_);
        smlib_brep_free(plate);
    }

    // A revolved hub (r 8..30, flange 6 tall, boss r 8..18 to 22) and its flange's top annulus
    // extruded 5 up from its edges, each circle two half arcs (SolidCraft rig sweeps/494), joined.
    // The ring's outer wall halves merge into the flange's cylinder; the line between them at 180
    // degrees then has that face on both sides, and used to stay as a slit (it is no intersection
    // edge, which were all BooleanPostProcess checked). Left: the three seams. Pappus:
    // pi (30^2 - 8^2) 6 + pi (18^2 - 8^2) 16 + pi (30^2 - 18^2) 5 = 12056 pi.
    {
        const double outline[6][2] = {{8, 0}, {30, 0}, {30, 6}, {18, 6}, {18, 22}, {8, 22}};
        std::vector<int64_t> kinds;
        std::vector<double> reals;
        for (int i = 0; i < 6; i++)
        {
            const double *p = outline[i], *q = outline[(i + 1) % 6];
            kinds.push_back(SMLIB_CURVE_LINE);
            reals.insert(reals.end(), {p[0], 0, p[1], q[0], 0, q[1]});
        }
        const double origin[3] = {0, 0, 0}, up[3] = {0, 0, 1};
        smlib_brep *hub = nullptr, *ring = nullptr, *joined = nullptr;
        CHECK(smlib_revolve(kinds.data(), 6, reals.data(), static_cast<int64_t>(reals.size()), origin, up, 2 * PI, 1, &hub) == SMLIB_OK);
        Profile annulus;
        annulus.arc(0, 0, 6, 1, 30, 0, PI).arc(0, 0, 6, 1, 30, PI, 2 * PI).arc(0, 0, 6, -1, 18, -PI, 0).arc(0, 0, 6, -1, 18, -2 * PI, -PI);
        ring = annulus.extrude(5.0);
        CHECK(hub != nullptr && ring != nullptr && smlib_boolean(hub, ring, SMLIB_BOOLEAN_UNION, &joined) == SMLIB_OK);
        if (joined != nullptr)
        {
            CHECK(near(volume_of(joined), 12056.0 * PI, 1e-6));
            CHECK(sound(joined));
            SmTArray<SmEdge *> edges;
            joined->p->GetEdges(edges);
            int inside = 0;
            for (ULONG i = 0; i < edges.GetSize(); i++)
            {
                SmTArray<SmFace *> sides;
                edges[i]->GetFaces(sides);
                inside += sides.GetSize() == 1 ? 1 : 0;
            }
            std::printf("hub joined with a half-arc ring: %lu edges, %d with one face\n", static_cast<unsigned long>(edges.GetSize()), inside);
            CHECK(inside == 3);
        }
        for (smlib_brep *x : {hub, ring, joined})
        {
            smlib_brep_free(x);
        }
    }

    // Four r 2 disks in one plane at the corners of a 3 m sheet's hole pattern, joined in pairs and
    // then all four (SolidCraft display_mesh, a 3 m sheet's 2D boolean). Boolean2D rebuilds the
    // second operand's faces on copies of the first's plane after growing it to cover them; a
    // plane only grew 100 face widths (B-spline reach), so faces 1.5 m apart kept a plane 6 mm
    // wide, their edges up to 1470 off it: tolerances 5875 and 11756, and the four-disk join 1006.
    {
        const double centres[4][2] = {{30, 20}, {1500, 20}, {2970, 980}, {30, 980}};
        smlib_brep *disks[4] = {nullptr, nullptr, nullptr, nullptr};
        for (int i = 0; i < 4; i++)
        {
            const int64_t kind = SMLIB_CURVE_ARC;
            const double arc[12] = {centres[i][0], centres[i][1], 0, 0, 0, 1, 1, 0, 0, 2, 0, 2 * PI};
            CHECK(smlib_planar_faces(&kind, 1, arc, 12, &disks[i]) == SMLIB_OK);
        }
        smlib_brep *low = nullptr, *high = nullptr, *all = nullptr;
        CHECK(smlib_boolean2d(disks[0], disks[1], SMLIB_BOOLEAN_UNION, &low) == SMLIB_OK);
        CHECK(smlib_boolean2d(disks[2], disks[3], SMLIB_BOOLEAN_UNION, &high) == SMLIB_OK);
        CHECK(low != nullptr && high != nullptr && smlib_boolean2d(low, high, SMLIB_BOOLEAN_UNION, &all) == SMLIB_OK);
        for (smlib_brep *x : {low, high, all})
        {
            if (x == nullptr)
            {
                continue;
            }
            smlib_flat *f = nullptr;
            CHECK(smlib_brep_export(x, &f) == SMLIB_OK);
            const double *tol = nullptr;
            int64_t n = 0;
            CHECK(f != nullptr && smlib_flat_get_f64(f, SMLIB_FLAT_EDGE_TOLERANCE, &tol, &n) == SMLIB_OK);
            double worst = 0.0;
            for (int64_t i = 0; i < n; i++)
            {
                worst = std::max(worst, tol[i]);
            }
            std::printf("disks 1.5 m apart joined in 2D: %lld edges, largest tolerance %.3g\n", static_cast<long long>(n), worst);
            CHECK(worst <= 1e-4);
            CHECK(worst_pcurve_gap(x) <= 1.0);
            smlib_flat_free(f);
        }
        CHECK(all != nullptr && all->p->GetNumFaces() == 4);
        for (smlib_brep *x : {disks[0], disks[1], disks[2], disks[3], low, high, all})
        {
            smlib_brep_free(x);
        }
    }

    // A plate perforated by k x k r 2 holes on a 20 mm pitch, the holes united in balanced rounds
    // and cut from the plate in 2D (SolidCraft's perforated sheet, a thousand holes). Boolean2D
    // grows every face's plane to cover the whole sheet, so the coincident-surface check passed
    // for every pair of faces and each pair classified its edges on the other face, and each
    // hole imprinted in the plate classified every hole imprinted before it: time in n^2 (a
    // thousand holes over 20 minutes in SolidCraft, 256 here 14 times as long as 64). Four times
    // the holes must take well under sixteen times as long, and the plate keeps its exact area.
    {
        const double radius = 2.0, pitch = 20.0;
        auto plate_seconds = [&](int k) -> double {
            const auto t0 = std::chrono::steady_clock::now();
            std::vector<smlib_brep *> holes;
            for (int i = 0; i < k; i++)
            {
                for (int j = 0; j < k; j++)
                {
                    const int64_t kind = SMLIB_CURVE_ARC;
                    const double arc[12] = {pitch * (i + 1), pitch * (j + 1), 0, 0, 0, 1, 1, 0, 0, radius, 0, 2 * PI};
                    smlib_brep *d = nullptr;
                    CHECK(smlib_planar_faces(&kind, 1, arc, 12, &d) == SMLIB_OK);
                    holes.push_back(d);
                }
            }
            while (holes.size() > 1)
            {
                std::vector<smlib_brep *> next;
                for (size_t h = 0; h + 1 < holes.size(); h += 2)
                {
                    smlib_brep *u = nullptr;
                    CHECK(smlib_boolean2d(holes[h], holes[h + 1], SMLIB_BOOLEAN_UNION, &u) == SMLIB_OK);
                    smlib_brep_free(holes[h]);
                    smlib_brep_free(holes[h + 1]);
                    next.push_back(u);
                }
                if (holes.size() % 2 == 1)
                {
                    next.push_back(holes.back());
                }
                holes = next;
            }
            const double w = pitch * (k + 1);
            const int64_t kinds[4] = {SMLIB_CURVE_LINE, SMLIB_CURVE_LINE, SMLIB_CURVE_LINE, SMLIB_CURVE_LINE};
            const double sides[24] = {0, 0, 0, w, 0, 0, w, 0, 0, w, w, 0, w, w, 0, 0, w, 0, 0, w, 0, 0, 0, 0};
            smlib_brep *plate = nullptr, *cut = nullptr;
            CHECK(smlib_planar_faces(kinds, 4, sides, 24, &plate) == SMLIB_OK);
            CHECK(holes.size() == 1 && holes[0] != nullptr && holes[0]->p->GetNumFaces() == static_cast<ULONG>(k * k));
            CHECK(plate != nullptr && smlib_boolean2d(plate, holes[0], SMLIB_BOOLEAN_SUBTRACT, &cut) == SMLIB_OK);
            const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
            double area = 0.0;
            CHECK(cut != nullptr && smlib_brep_area(cut, 1e-9, &area) == SMLIB_OK);
            const double exact = w * w - k * k * PI * radius * radius;
            std::printf("plate with %d holes in 2D: %.3f s, area %.6f, exact %.6f\n", k * k, seconds, area, exact);
            CHECK(near(area, exact, 1e-7));
            for (smlib_brep *x : {holes[0], plate, cut})
            {
                smlib_brep_free(x);
            }
            return seconds;
        };
        const double t64 = plate_seconds(8), t256 = plate_seconds(16);
        CHECK(t256 <= 8.0 * t64 + 0.1);
    }

    smlib_brep *rod = cylinder(10, 10, -1, 4, 22);
    CHECK(smlib_boolean(cube, rod, SMLIB_BOOLEAN_SUBTRACT, &r) == SMLIB_OK);
    CHECK(near(volume_of(r), 8000.0 - PI * 16.0 * 20.0, 1e-6));
    CHECK(sound(r));
    smlib_brep *holed = r;
    r = nullptr;
    CHECK(smlib_boolean(cube, rod, SMLIB_BOOLEAN_UNION, &r) == SMLIB_OK);
    CHECK(near(volume_of(r), 8000.0 + PI * 16.0 * 2.0, 1e-6));
    done(r);
    // A blind hole, ending inside.
    smlib_brep *pin = cylinder(10, 10, 10, 4, 15);
    CHECK(smlib_boolean(cube, pin, SMLIB_BOOLEAN_SUBTRACT, &r) == SMLIB_OK);
    CHECK(near(volume_of(r), 8000.0 - PI * 16.0 * 10.0, 1e-6));
    done(r);

    // The rim of the hole is one closed circular edge.
    // At 45 degrees: the cylinder's seam meets the rim at 0 degrees, and a point there is
    // as near to the seam as to the rim.
    const double rim[3] = {10 + 4 * std::sqrt(0.5), 10 + 4 * std::sqrt(0.5), 20};
    int32_t s = smlib_fillet(holed, rim, 1, 0.1, 1.0, &r, nullptr);
    if (s == SMLIB_OK)
    {
        // A round of r = 1 on a hole of radius R = 4: the cross-section (1 − π/4)·r² goes, its
        // centroid at R + r·(10 − 3π)/(12 − 3π) from the axis (Pappus).
        double area = (1.0 - PI / 4.0);
        double c = 4.0 + (10.0 - 3.0 * PI) / (12.0 - 3.0 * PI);
        CHECK(near(volume_of(r), volume_of(holed) - 2.0 * PI * c * area, 1e-5));
        done(r);
    }
    else
    {
        std::fprintf(stderr, "fillet of the hole rim failed: %d\n", static_cast<int>(s));
        failures++;
    }
    s = smlib_chamfer(holed, rim, 1, 0.1, 1.0, &r, nullptr);
    if (s == SMLIB_OK)
    {
        // A 1 × 1 triangle, centroid at R + 1/3.
        CHECK(near(volume_of(r), volume_of(holed) - 2.0 * PI * (4.0 + 1.0 / 3.0) * 0.5, 1e-5));
        done(r);
    }
    else
    {
        std::fprintf(stderr, "chamfer of the hole rim failed: %d\n", static_cast<int>(s));
        failures++;
    }
    // tangency-labs/solidcore#379: SmApiChamferFillet (SMLib's constant-distance solver, a
    // chord of 1 across each chamfer) on every edge of the holed cube, the hole's rims too,
    // gives a closed solid, the one the shim's chamfer (constant radius, linear section) gives
    // at the setback a chord of 1 has on these square edges, 1 / sqrt(2).
    {
        SmBrep *all = nullptr;
        CHECK(SmApiBrepCopy(holed->p, all) == SM_SUCCESS && all != nullptr);
        if (all != nullptr)
        {
            CHECK(SmApiChamferFillet(all, 1.0) == SM_SUCCESS);
            smlib_brep *chamfered = new smlib_brep{all};
            CHECK(all->IsManifoldSolid() && sound(chamfered));
            double edges[42];
            int n = 0;
            for (double a : {0.0, 20.0})
            {
                for (double b : {0.0, 20.0})
                {
                    const double rows[3][3] = {{10, a, b}, {a, 10, b}, {a, b, 10}};
                    for (const auto &row : rows)
                    {
                        for (double x : row)
                        {
                            edges[n++] = x;
                        }
                    }
                }
            }
            for (double z : {0.0, 20.0})
            {
                edges[n++] = rim[0];
                edges[n++] = rim[1];
                edges[n++] = z;
            }
            CHECK(smlib_chamfer(holed, edges, 14, 0.1, 1.0 / std::sqrt(2.0), &r, nullptr) == SMLIB_OK);
            if (r != nullptr)
            {
                CHECK(near(volume_of(chamfered), volume_of(r), 1e-8));
                done(r);
            }
            smlib_brep_free(chamfered);
        }
    }
    const double top[3] = {2, 2, 20};
    if (smlib_shell(holed, top, 1, 0.1, -1.0, &r, nullptr) == SMLIB_OK)
    {
        done(r);
    }

    // Every edge of the cube, then a tangent chain.
    double mids[36];
    int k = 0;
    const double ends[4][2] = {{0, 0}, {20, 0}, {0, 20}, {20, 20}};
    for (auto &e : ends)
    {
        double rows[3][3] = {{10, e[0], e[1]}, {e[0], 10, e[1]}, {e[0], e[1], 10}};
        for (auto &row : rows)
        {
            for (double x : row)
            {
                mids[k++] = x;
            }
        }
    }
    CHECK(smlib_fillet(cube, mids, 12, 0.1, 2.0, &r, nullptr) == SMLIB_OK);
    done(r);
    // tangency-labs/solidcore#1176: the clean-up after a fillet kept the target's edges and faces
    // by address. Filleting deletes some of them, and an edge or face it then made at a freed
    // address passed for the target's own: BooleanPostProcess passed over it without its
    // topological test, and so without the UV trims that test builds. The cube's 12 edges
    // chamfered 1/sqrt(2) came out with 96 UV trims or 94 (two edgeuses without, not always the
    // same two) as the heap lay. Run again with blocks of other sizes allocated first, some of them
    // freed, the body written out must be the same byte for byte.
    {
        auto chamfered = [&]() -> std::string
        {
            const char *path = "solidcore-1176.smb";
            smlib_brep *c = nullptr;
            std::string text;
            CHECK(smlib_chamfer(cube, mids, 12, 0.1, 1.0 / std::sqrt(2.0), &c, nullptr) == SMLIB_OK);
            if (c != nullptr && smlib_brep_write(c, path, 1) == SMLIB_OK)
            {
                std::ifstream in(path, std::ios::binary);
                text.assign(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
            }
            std::remove(path);
            CHECK(!text.empty());
            done(c);
            return text;
        };
        const std::string first = chamfered();
        for (size_t layout = 1; layout <= 5; ++layout)
        {
            std::vector<void *> blocks;
            for (size_t round = 0; round < layout; ++round)
            {
                for (size_t size = 16; size <= 1024; size += 8)
                {
                    blocks.push_back(std::malloc(size));
                }
            }
            for (size_t i = 0; i < blocks.size(); i += layout + 1)
            {
                std::free(blocks[i]);
                blocks[i] = nullptr;
            }
            const bool same = chamfered() == first;
            if (!same)
            {
                std::fprintf(stderr, "  the chamfered cube differs after heap layout %zu\n", layout);
            }
            CHECK(same);
            for (void *b : blocks)
            {
                std::free(b);
            }
        }
    }
    const double vertical[3] = {20, 20, 10};
    smlib_brep *round = nullptr;
    CHECK(smlib_fillet(cube, vertical, 1, 0.1, 5.0, &round, nullptr) == SMLIB_OK && sound(round));
    const double straight[3] = {7.5, 20, 20};
    CHECK(smlib_fillet(round, straight, 1, 0.1, 1.0, &r, nullptr) == SMLIB_OK);
    done(r);
    CHECK(smlib_chamfer(round, straight, 1, 0.1, 1.0, &r, nullptr) == SMLIB_OK);
    done(r);

    // tangency-labs/solidcore#420: a chamfer as wide as the corner rounds it runs round. A box
    // 40 x 40 x 20 with its vertical edges rounded R1, its bottom loop chamfered 1: on each
    // corner the round's arc gets a cone closing to a point on the round's axis (SMLib traced
    // it as an offset intersection, which collapses to that point, and failed). Each straight
    // side loses a 1 x 1 triangle along 38, each corner the triangle turned a quarter about the
    // axis: (pi / 2) * (1 / 3).
    {
        smlib_brep *plate = box(-20, -20, 0, 40, 40, 20);
        const double corners[12] = {20, 20, 10, -20, 20, 10, 20, -20, 10, -20, -20, 10};
        smlib_brep *rounded = nullptr;
        CHECK(smlib_fillet(plate, corners, 4, 0.1, 1.0, &rounded, nullptr) == SMLIB_OK);
        // The whole bottom loop: its four sides and the four arcs of the corner rounds.
        const double h = 19.0 + std::sqrt(0.5);
        const double bottom[24] = {0, -20, 0, 20, 0, 0, 0, 20, 0, -20, 0, 0, h, h, 0, -h, h, 0, h, -h, 0, -h, -h, 0};
        if (rounded != nullptr && smlib_chamfer(rounded, bottom, 8, 0.1, 1.0, &r, nullptr) == SMLIB_OK)
        {
            CHECK(near(volume_of(r), volume_of(rounded) - 4.0 * 38.0 * 0.5 - 4.0 * (PI / 2.0) / 3.0, 1e-7));
            done(r);
        }
        else
        {
            CHECK(!"the bottom loop of the rounded plate chamfered 1");
        }
        smlib_brep_free(rounded);
        smlib_brep_free(plate);
    }

    // tangency-labs/solidcore#417: a box 297 x 100 x 60 with its vertical edges rounded R5, then
    // its top rim rounded R5: on each corner the round of the corner round's arc is a ball
    // centered on that round's axis, which SMLib's torus turned the long way round (the
    // section's rail points, taken where the ball's offset faces collapse onto the axis, carry
    // no angle; on a smaller box it happened to turn the right way). Each straight edge loses
    // (1 - pi/4) 25 along 287 or 90, each corner the square 5 x 5 less the quarter disc, turned
    // a quarter about the axis: (pi/2) (62.5 - 125/3).
    {
        smlib_brep *block = box(0, 0, 0, 297, 100, 60);
        const double corners[12] = {0, 0, 30, 297, 0, 30, 0, 100, 30, 297, 100, 30};
        smlib_brep *rounded = nullptr;
        CHECK(smlib_fillet(block, corners, 4, 0.1, 5.0, &rounded, nullptr) == SMLIB_OK);
        // The whole top rim: its four sides and the four arcs of the corner rounds.
        const double a = 5.0 - 5.0 * std::sqrt(0.5);
        const double rim[24] = {148.5, 0, 60, 297, 50, 60, 148.5, 100, 60, 0, 50, 60, a, a, 60, 297 - a, a, 60, a, 100 - a, 60, 297 - a, 100 - a, 60};
        if (rounded != nullptr && smlib_fillet(rounded, rim, 8, 0.1, 5.0, &r, nullptr) == SMLIB_OK)
        {
            double removed = 754.0 * 25.0 * (1.0 - PI / 4.0) + 4.0 * (PI / 2.0) * (62.5 - 125.0 / 3.0);
            CHECK(near(volume_of(r), volume_of(rounded) - removed, 1e-7));
            done(r);
        }
        else
        {
            CHECK(!"the top rim of the block with R5 corners rounded R5");
        }
        smlib_brep_free(rounded);
        smlib_brep_free(block);
    }

    // tangency-labs/solidcore#416: a box 25 x 25 x 10 with its bottom edges rounded R2, then its
    // vertical edges R2. Each vertical edge's tangent chain runs on down the seam where two
    // bottom rounds meet, between two cylinders as wide as the round: the ball touching both
    // sits where their axes meet, the corner's ball, and SMLib traced the offset cylinders'
    // intersection (their axes, crossing at that one point) and failed. Each vertical edge
    // loses (1 - pi/4) 4 along 8, each corner the octant of the two bottom rounds' crossing
    // (16/3, a Steinmetz solid's eighth) less the corner ball's (4 pi / 3).
    {
        smlib_brep *tray = box(0, 0, 0, 25, 25, 10);
        const double bottom[12] = {12.5, 0, 0, 25, 12.5, 0, 12.5, 25, 0, 0, 12.5, 0};
        smlib_brep *rounded = nullptr;
        CHECK(smlib_fillet(tray, bottom, 4, 0.1, 2.0, &rounded, nullptr) == SMLIB_OK);
        // Each vertical edge and the seam where the two bottom rounds below it meet, which
        // continues it tangentially (a point of it: x = y = z = 2 - sqrt 2 at the origin corner).
        const double s2 = 2.0 - std::sqrt(2.0);
        const double vertical[24] = {0, 0, 6, 25, 0, 6, 25, 25, 6, 0, 25, 6, s2, s2, s2, 25 - s2, s2, s2, 25 - s2, 25 - s2, s2, s2, 25 - s2, s2};
        if (rounded != nullptr && smlib_fillet(rounded, vertical, 8, 0.1, 2.0, &r, nullptr) == SMLIB_OK)
        {
            double removed = 4.0 * (32.0 * (1.0 - PI / 4.0) + 16.0 / 3.0 - 4.0 * PI / 3.0);
            CHECK(near(volume_of(r), volume_of(rounded) - removed, 1e-7));
            done(r);
        }
        else
        {
            CHECK(!"the vertical edges of the tray with R2 bottom rounds rounded R2");
        }
        smlib_brep_free(rounded);
        smlib_brep_free(tray);
    }

    // The same tray as SolidCraft hands it over (rig fillets/1557): extruded from its outline and
    // read back through the flat arrays, its bottom rounds read back as cylinders (SMLib's own
    // cone of semi-angle 0, as the record keeps it), which puts the rounds' axes a rounding apart
    // from SMLib's own. At the bottom face's corner, where the contact circles of the corner's
    // ball meet, the trace's point on the seam had its centre 2.8e-4 off where the axes meet,
    // so the ball was not made, and the traced seam round failed (1001). The ball's centre is
    // where the axes meet: the same tray, the same volume.
    {
        // `b` written to the flat arrays and read back, its cones of semi-angle 0 as cylinders
        auto read_back = [](const smlib_brep *b, bool cylinders) {
            smlib_flat *f = nullptr;
            CHECK(smlib_brep_export(b, &f) == SMLIB_OK);
            const int64_t *k = nullptr;
            int64_t n = 0;
            CHECK(smlib_flat_get_i64(f, SMLIB_FLAT_FACE_SURFACE_TYPE, &k, &n) == SMLIB_OK);
            std::vector<int64_t> kinds(k, k + n);
            if (cylinders && std::find(kinds.begin(), kinds.end(), int64_t{SMLIB_FLAT_SURFACE_CONE}) != kinds.end())
            {
                std::replace(kinds.begin(), kinds.end(), int64_t{SMLIB_FLAT_SURFACE_CONE}, int64_t{SMLIB_FLAT_SURFACE_CYLINDER});
                const int cone[4] = {SMLIB_FLAT_FACE_CONE_ORIGIN, SMLIB_FLAT_FACE_CONE_AXIS, SMLIB_FLAT_FACE_CONE_REF_DIRECTION, SMLIB_FLAT_FACE_CONE_RADIUS};
                const int cylinder[4] = {SMLIB_FLAT_FACE_CYLINDER_ORIGIN, SMLIB_FLAT_FACE_CYLINDER_AXIS, SMLIB_FLAT_FACE_CYLINDER_REF_DIRECTION,
                                         SMLIB_FLAT_FACE_CYLINDER_RADIUS};
                for (int i = 0; i < 4; ++i)
                {
                    const double *d = nullptr;
                    int64_t len = 0;
                    CHECK(smlib_flat_get_f64(f, cone[i], &d, &len) == SMLIB_OK);
                    std::vector<double> v(d, d + len);
                    CHECK(smlib_flat_set_f64(f, cylinder[i], v.data(), len) == SMLIB_OK);
                    CHECK(smlib_flat_set_f64(f, cone[i], nullptr, 0) == SMLIB_OK);
                }
                CHECK(smlib_flat_set_f64(f, SMLIB_FLAT_FACE_CONE_SEMI_ANGLE, nullptr, 0) == SMLIB_OK);
                CHECK(smlib_flat_set_i64(f, SMLIB_FLAT_FACE_SURFACE_TYPE, kinds.data(), n) == SMLIB_OK);
            }
            smlib_brep *out = nullptr;
            CHECK(smlib_brep_import(f, &out) == SMLIB_OK);
            smlib_flat_free(f);
            return out;
        };
        const int64_t kinds[4] = {SMLIB_CURVE_LINE, SMLIB_CURVE_LINE, SMLIB_CURVE_LINE, SMLIB_CURVE_LINE};
        const double outline[24] = {0, 0, 0, 25, 0, 0, 25, 0, 0, 25, 25, 0, 25, 25, 0, 0, 25, 0, 0, 25, 0, 0, 0, 0};
        const double up[3] = {0, 0, 1};
        smlib_brep *extruded = nullptr;
        CHECK(smlib_extrude(kinds, 4, outline, 24, up, 10.0, 0.0, 0, 1, &extruded, nullptr) == SMLIB_OK);
        smlib_brep *tray = extruded != nullptr ? read_back(extruded, false) : nullptr;
        // SolidCraft's pick tolerance for the tray, and the points it picks the seams by
        const double tol = 1.8381173070873837;
        const double bottom[12] = {12.5, 0, 0, 25, 12.5, 0, 12.5, 25, 0, 0, 12.5, 0};
        const double sizes[8] = {2, 2, 2, 2, 2, 2, 2, 2};
        smlib_brep *rounded = nullptr;
        CHECK(tray != nullptr && smlib_blend(tray, SMLIB_BLEND_ROUND, bottom, sizes, 4, tol, &rounded, nullptr) == SMLIB_OK);
        smlib_brep *handed = rounded != nullptr ? read_back(rounded, true) : nullptr;
        const double vertical[24] = {0, 0, 6, 25, 0, 6, 25, 25, 6, 0, 25, 6,
                                     0.71320933207705417, 0.71320933207757631, 0.47032672374040252,
                                     24.286791299790224, 0.71320870020977245, 0.47032422391299583,
                                     24.28679129979017, 24.286791299790359, 0.47032422391295658,
                                     0.7132087002097981, 24.286791299790536, 0.47032422391285106};
        if (handed != nullptr && smlib_blend(handed, SMLIB_BLEND_ROUND, vertical, sizes, 8, tol, &r, nullptr) == SMLIB_OK)
        {
            // the inner 21 x 21 x 12 box grown by a ball of R2, cut at z = 10: the slab up to
            // z = 2 (the core 21 x 21 x 2, quarter discs along its 84 of edges, corners a hemisphere) and the
            // rounded square above it, 8 high
            const double want = 21.0 * 21.0 * 2.0 + 42.0 * PI * 2.0 + 16.0 * PI / 3.0 + 8.0 * (625.0 - (4.0 - PI) * 4.0);
            CHECK(near(volume_of(r), want, 1e-7));
            done(r);
        }
        else
        {
            CHECK(!"the vertical edges of the extruded tray with R2 bottom rounds as cylinders rounded R2");
        }
        for (smlib_brep *b : {handed, rounded, tray, extruded})
        {
            smlib_brep_free(b);
        }
    }

    // SolidCraft issue 851's part: a 41 1/6 x 21 plate and the corner a circle of radius 21
    // leaves, 10 high, the arc wall running out into the back face in a knife edge. Its top front
    // edge rounded R2, then the arc wall's top edge R2: that round ends on the front round where
    // both are tangent to the top face, a singular point of their intersection with four
    // directions out of it. The start direction kept the first of them, not the one pointing
    // into both surfaces nearest the given direction, unless it was found later in the list: on
    // the part read back through the flat arrays the first one ran along the top face square to
    // the round's end, the end came out a point, and the round failed (1001). Natively the
    // first direction happened to be the right one reversed. Both bodies now give the round, the
    // same.
    {
        const double tip = 62.0 + 1.0 / 6.0, wall = 21.0, x0 = tip - wall;
        auto read_back = [](const smlib_brep *b) {
            smlib_flat *f = nullptr;
            CHECK(smlib_brep_export(b, &f) == SMLIB_OK);
            smlib_brep *out = nullptr;
            CHECK(smlib_brep_import(f, &out) == SMLIB_OK);
            smlib_flat_free(f);
            return out;
        };
        smlib_brep *plate = box(0, -wall, 0, x0, wall, 10);
        smlib_brep *corner = box(x0, -wall, 0, wall, wall, 10);
        const double axis[3] = {tip, -wall, 0};
        smlib_brep *disc = nullptr;
        CHECK(smlib_make_cylinder(axis, wall, 10, &disc) == SMLIB_OK);
        smlib_brep *cut = nullptr;
        CHECK(smlib_boolean(corner, disc, SMLIB_BOOLEAN_SUBTRACT, &cut) == SMLIB_OK);
        smlib_brep *part = nullptr;
        CHECK(cut != nullptr && smlib_boolean(plate, cut, SMLIB_BOOLEAN_UNION, &part) == SMLIB_OK);
        // SolidCraft's pick tolerance for the part, and the points it picks the edges by
        const double tol = 3.3197702709152845;
        const double front[3] = {20.583333333333332, -21, 10};
        const double arc[3] = {48.045114044109852, -5.4635486396119797, 10};
        const double size = 2.0;
        double volumes[2] = {0.0, 0.0};
        for (int read = 0; read < 2 && part != nullptr; ++read)
        {
            smlib_brep *rounded = nullptr;
            CHECK(smlib_blend(part, SMLIB_BLEND_ROUND, front, &size, 1, tol, &rounded, nullptr) == SMLIB_OK);
            smlib_brep *start = rounded != nullptr && read == 1 ? read_back(rounded) : rounded;
            if (start != nullptr && smlib_blend(start, SMLIB_BLEND_ROUND, arc, &size, 1, tol, &r, nullptr) == SMLIB_OK)
            {
                volumes[read] = volume_of(r);
                CHECK(volumes[read] < volume_of(start));
                done(r);
            }
            else
            {
                CHECK(!"the arc wall's top edge rounded R2 after the top front edge");
            }
            if (start != rounded)
            {
                smlib_brep_free(start);
            }
            smlib_brep_free(rounded);
        }
        CHECK(near(volumes[1], volumes[0], 1e-7));
        for (smlib_brep *b : {part, cut, disc, corner, plate})
        {
            smlib_brep_free(b);
        }
    }

    // tangency-labs/solidcore#453: a 20 cube with a hole r5 along y whose wall touches the
    // bottom face along one line, its back rim chamfered 3. With the rim's tangent chain the
    // bottom edge came too, four filleted edges met at the touching point among five, a corner
    // no corner type handles, and the smooth-sequence walk called it as a corner (a crash). No
    // chamfer of that size exists there (the rim has no face below it): refused, every kind.
    {
        SmBrep *bore = nullptr;
        CHECK(SmApiCreateCylinder(SmVector3d(0, 0, 0), 5.0, 20.0, bore) == SM_SUCCESS && bore != nullptr);
        CHECK(SmApiRotate(bore, SmVector3d(0, 0, 0), SmVector3d(1, 0, 0), 90.0) == SM_SUCCESS);
        CHECK(SmApiTranslate(bore, SmVector3d(10, 20, 5)) == SM_SUCCESS);
        smlib_brep tool{bore};
        smlib_brep *holed_low = nullptr;
        CHECK(smlib_boolean(cube, &tool, SMLIB_BOOLEAN_SUBTRACT, &holed_low) == SMLIB_OK);
        CHECK(near(volume_of(holed_low), 8000.0 - PI * 25.0 * 20.0, 1e-6));
        const double back_rim[3] = {10, 20, 10};
        CHECK(smlib_chamfer(holed_low, back_rim, 1, 0.1, 3.0, &r, nullptr) != SMLIB_OK && r == nullptr);
        for (int32_t kind : {SMLIB_BLEND_ROUND, SMLIB_BLEND_CURVATURE, SMLIB_BLEND_CHORD, SMLIB_BLEND_CHAMFER})
        {
            const double size = 3.0;
            int64_t failed = -1;
            CHECK(smlib_blend(holed_low, kind, back_rim, &size, 1, 0.1, &r, &failed) != SMLIB_OK && r == nullptr);
        }
        smlib_brep_free(holed_low);
        delete bore;
    }

    // tangency-labs/solidcore#418: a bar of stadium section (25 x 19 and two half discs r9.5)
    // 14 long along x, a slot x 3.5..9.5 cutting away the y < 0 half of one round end, and the
    // edge between the flat top and the slot floor (6 long) rounded r. At each of its ends the
    // top runs on tangentially onto the round end; the corner past that tangent edge is the
    // slot wall, whose extension the round ends on. SMLib took the tangent edge for the side
    // edge and the round for the end face (1001). Each round takes (1 - pi/4) r^2 along 6, a
    // chamfer of setback 2 the triangle of 2 x 2.
    {
        const int64_t kinds[4] = {SMLIB_CURVE_LINE, SMLIB_CURVE_ARC, SMLIB_CURVE_LINE, SMLIB_CURVE_ARC};
        const double reals[36] = {0, 0,  -9.5, 0, 25, -9.5, 0, 25, 0,   1, 0, 0, 0, 0, -1, 9.5, 0, PI,
                                  0, 25, 9.5,  0, 0,  9.5,  0, 0,  0,   1, 0, 0, 0, 0, 1,  9.5, 0, PI};
        const double along[3] = {1, 0, 0};
        smlib_brep *bar = nullptr;
        CHECK(smlib_extrude(kinds, 4, reals, 36, along, 14.0, 0.0, SMLIB_DRAFT_CORNER_TANGENT, 1, &bar, nullptr) == SMLIB_OK);
        smlib_brep *slot = box(3.5, -12, -10, 6, 12, 20);
        smlib_brep *slotted = nullptr;
        CHECK(smlib_boolean(bar, slot, SMLIB_BOOLEAN_SUBTRACT, &slotted) == SMLIB_OK && sound(slotted));
        const double floor_edge[3] = {6.5, 0, 9.5};
        for (double radius : {0.5, 2.0, 4.0})
        {
            if (slotted != nullptr && smlib_fillet(slotted, floor_edge, 1, 0.1, radius, &r, nullptr) == SMLIB_OK)
            {
                CHECK(near(volume_of(r), volume_of(slotted) - 6.0 * radius * radius * (1.0 - PI / 4.0), 1e-8));
                done(r);
            }
            else
            {
                CHECK(!"the slot floor's top edge of the stadium bar rounded");
            }
        }
        if (slotted != nullptr && smlib_chamfer(slotted, floor_edge, 1, 0.1, 2.0, &r, nullptr) == SMLIB_OK)
        {
            CHECK(near(volume_of(r), volume_of(slotted) - 6.0 * 2.0, 1e-8));
            done(r);
        }
        else
        {
            CHECK(!"the slot floor's top edge of the stadium bar chamfered");
        }
        for (smlib_brep *b : {bar, slot, slotted})
        {
            smlib_brep_free(b);
        }
    }

    // A constant radius law (r -> r) on the top back edge of a 60 x 40 x 20 box whose vertical
    // edge at (60, 40) is rounded R5. The round runs on past the edge's end until it meets the R5
    // cylinder, so the law is evaluated past its end; built with no slope (the shim's law curve
    // had radius in x and 0 in y), it could not be extended there, and the round failed (1001).
    // Removed: (1 - pi/4) r^2 along 55, and past x = 55 the square r x r less the quarter disc,
    // as far as the cylinder reaches: the integral of (r - sqrt(r^2 - u^2)) sqrt(25 - (5 - r +
    // u)^2) over 0 < u < r.
    {
        smlib_brep *block = box(0, 0, 0, 60, 40, 20);
        const double corner_edge[3] = {60, 40, 10};
        smlib_brep *rounded = nullptr;
        CHECK(smlib_fillet(block, corner_edge, 1, 0.1, 5.0, &rounded, nullptr) == SMLIB_OK);
        const double back_top[3] = {27.5, 40, 20};
        const double from[3] = {0, 40, 20};
        const double knots[4] = {0, 0, 1, 1};
        for (double radius : {1.0, 2.0, 3.0})
        {
            const double law[2] = {radius, radius};
            if (rounded != nullptr && smlib_fillet_variable(rounded, back_top, 0.1, from, 1, knots, 4, law, 2, &r) == SMLIB_OK)
            {
                double past = 0.0;
                const int steps = 200000;
                for (int i = 0; i < steps; ++i)
                {
                    double u = (i + 0.5) * radius / steps;
                    double w = 5.0 - radius + u;
                    past += (radius - std::sqrt(radius * radius - u * u)) * std::sqrt(std::max(0.0, 25.0 - w * w)) * radius / steps;
                }
                CHECK(near(volume_of(r), volume_of(rounded) - 55.0 * radius * radius * (1.0 - PI / 4.0) - past, 1e-6));
                done(r);
            }
            else
            {
                CHECK(!"a constant law into the R5 corner");
            }
        }
        smlib_brep_free(rounded);
        smlib_brep_free(block);
    }

    // A variable round 2 -> w -> 2 (degree 1, a corner at the middle) on the top edge along y of a
    // 10 x 20 x 30 box. Corners in the rails were found by the solver's 20 degree angle, so a
    // gentle one (w 3: the rails turn 11 degrees) was skinned over: the round off its law, its
    // edges kinked inside and refused by the check. Each half is a cone-like round with r linear
    // along 10: removed about 2 (1 - pi/4) 10 (4 + 2 w + w^2) / 3 (the ball's contact leans off
    // the cross-section plane as r changes, to 1e-4 of it here).
    {
        smlib_brep *block = box(0, 0, 0, 10, 20, 30);
        const double at[3] = {10, 10, 30}, from[3] = {10, 0, 30};
        const double knots[5] = {0, 0, 0.5, 1, 1};
        for (double w : {2.5, 3.0, 3.2, 4.0})
        {
            const double law[3] = {2, w, 2};
            if (smlib_fillet_variable(block, at, 0.1, from, 1, knots, 5, law, 3, &r) == SMLIB_OK)
            {
                CHECK(near(6000.0 - volume_of(r), 20.0 * (1.0 - PI / 4.0) * (4.0 + 2.0 * w + w * w) / 3.0, 1e-4));
                CHECK(r->p->GetNumFaces() == 8);
                done(r);
            }
            else
            {
                CHECK(!"a variable round with a corner in its law");
            }
        }
        smlib_brep_free(block);
    }

    // tangency-labs/solidcore#458: a variable round whose law SMLib cannot evaluate past the edge's
    // ends (a constant law curve laid out with no slope, radius 2 as x and 0 as y, as the shim
    // once built it; extending a curve of zero speed fails). SmVariableRadiusFS::SetupOffsetValues
    // and SmTemporaryRadiusChange::SetUp ignored the failed evaluation and offset the surfaces by
    // an uninitialised radius (it went on to fail, so only a memory checker sees the read). Now
    // the evaluation's failure stops it: the round fails, or is the round of radius 2.
    {
        SmBrep *b = nullptr;
        CHECK(SmApiBrepCopy(cube->p, b) == SM_SUCCESS && b != nullptr);
        SmTArray<SmEdge *> all;
        b->GetEdges(all);
        SmTArray<SmEdge *> top_edge;
        for (ULONG i = 0; i < all.GetSize(); ++i)
        {
            SmVertex *v0 = nullptr;
            SmVertex *v1 = nullptr;
            all[i]->GetVertices(v0, v1);
            SmPoint3d m = (v0->GetPoint() + v1->GetPoint()) * 0.5;
            if (std::fabs(m.x - 10) < 1e-9 && std::fabs(m.y - 20) < 1e-9 && std::fabs(m.z - 20) < 1e-9)
            {
                top_edge.Add(all[i]);
            }
        }
        CHECK(top_edge.GetSize() == 1);
        const SmContext &ctx = *SmApiGetOrCreateContext();
        SmTArray<SmPoint3d> cps;
        cps.Add(SmPoint3d(2, 0, 0));
        cps.Add(SmPoint3d(2, 0, 0));
        SmTArray<ULONG> mults;
        mults.Add(2);
        mults.Add(2);
        SmTArray<double> unique;
        unique.Add(0.0);
        unique.Add(1.0);
        SmBSplineCurve *flat = nullptr;
        CHECK(SmBSplineCurve::CreateCanonical(ctx, 3, 1, cps, SM_CF_UNSPECIFIED, mults, unique, SM_KT_UNSPECIFIED, nullptr, nullptr,
                                              flat) == SM_SUCCESS && flat != nullptr);
        if (top_edge.GetSize() == 1 && flat != nullptr)
        {
            SmBSplineFilletLaw law(flat, SmExtent1d(0.0, 1.0));
            SmTArray<SmFilletLaw *> laws;
            laws.Add(&law);
            SmTArray<double> values;
            values.Add(2.0);
            SmTArray<SmFilletSurfaceGeneratorType> sections;
            sections.Add(SM_FSG_CIRCULAR);
            SmTArray<SmFilletSolverType> solvers;
            solvers.Add(SM_FS_VARIABLE_RADIUS);
            SmFilletExecutive exec(ctx, b);
            SmStatus st = exec.SetFilletParameters(b, top_edge, sections, solvers, values, &laws, -1.0, -1.0, -1.0, -1.0, 1, 1.0);
            if (st == SM_SUCCESS)
            {
                st = exec.DoFillet();
            }
            if (st == SM_SUCCESS && b->StitchAndOrient() == SM_SUCCESS)
            {
                smlib_brep got{b};
                CHECK(near(volume_of(&got), 8000.0 - 4.0 * 20.0 * (1.0 - PI / 4.0), 1e-6));
            }
        }
        delete b;
    }

    // tangency-labs/solidcore#430: a variable round whose radius law has a corner. SMLib traced
    // the round smoothly across the corner, averaging the tangents of its two sides, and skinned
    // one smooth surface over it: it missed the law past the corner and failed its own check.
    // The trace now ends a step at the corner with the tangents of each side, the surface is
    // skinned piece by piece, and the face is split there. On the cube's edge y = z = 20, the
    // round takes (1 - pi/4) r^2 per unit length; for 4 -> 2 -> 4 that is the same as for a
    // straight 4 -> 2 along the whole edge.
    {
        const double PI4 = 1.0 - PI / 4.0;
        const double edge[3] = {10, 20, 20};
        const double from0[3] = {0, 20, 20};
        const double from1[3] = {20, 20, 20};
        const double v_knots[5] = {0, 0, 0.5, 1, 1};
        const double v_radii[3] = {4, 2, 4};
        const double removed = PI4 * 20.0 * (16.0 + 8.0 + 4.0) / 3.0;
        for (const double *from : {from0, from1})
        {
            CHECK(smlib_fillet_variable(cube, edge, 0.1, from, 1, v_knots, 5, v_radii, 3, &r) == SMLIB_OK);
            if (r != nullptr)
            {
                CHECK(near(8000.0 - volume_of(r), removed, 1e-4));
                done(r);
            }
        }
        // Two corners (4, 2, 4, 2 at 0, 0.25, 0.75, 1): the straight 4 -> 2 in three pieces.
        const double w_knots[6] = {0, 0, 0.25, 0.75, 1, 1};
        const double w_radii[4] = {4, 2, 4, 2};
        CHECK(smlib_fillet_variable(cube, edge, 0.1, from0, 1, w_knots, 6, w_radii, 4, &r) == SMLIB_OK);
        if (r != nullptr)
        {
            CHECK(near(8000.0 - volume_of(r), removed, 1e-4));
            done(r);
        }
        // The same law as a quadratic with a corner at 0.5 (the knot repeated twice): 4, 3, 2
        // then 2, 3, 4, evenly spaced, are the straight halves again.
        const double q_knots[8] = {0, 0, 0, 0.5, 0.5, 1, 1, 1};
        const double q_radii[5] = {4, 3, 2, 3, 4};
        CHECK(smlib_fillet_variable(cube, edge, 0.1, from0, 2, q_knots, 8, q_radii, 5, &r) == SMLIB_OK);
        if (r != nullptr)
        {
            CHECK(near(8000.0 - volume_of(r), removed, 1e-4));
            done(r);
        }
        // A radius of zero inside the edge would pinch the round to a point: refused.
        const double p_radii[3] = {4, 0, 4};
        CHECK(smlib_fillet_variable(cube, edge, 0.1, from0, 1, v_knots, 5, p_radii, 3, &r) == SMLIB_ERR_UNSUPPORTED && r == nullptr);
    }

    // A variable round whose radius falls fast between its ends: the cubic Bezier law 6, 6, 1, 1
    // along the cube's edge y = z = 20 (its slope reaches 3/8). On the cube's planes any uv curve
    // lies on the face, so the trace's span test held for any step up to its angle tolerance,
    // and the rails between the traced points cut across the law: the round came out up to
    // 0.05 mm off its radius and took 0.2% too much. Each span's rails are now held to the
    // fillet solved inside the span. The round takes (1 - pi/4) r^2 per unit length, integrated
    // here (Simpson, 2000 panels, exact to 1e-12 on this polynomial); the sections are cubics
    // within 0.05% of the radius, so to 1e-4.
    {
        const double edge[3] = {10, 20, 20};
        const double from0[3] = {0, 20, 20};
        const double knots[8] = {0, 0, 0, 0, 1, 1, 1, 1};
        const double radii[4] = {6, 6, 1, 1};
        auto law = [&](double t) {
            const double u = 1.0 - t;
            return radii[0] * u * u * u + 3.0 * radii[1] * u * u * t + 3.0 * radii[2] * u * t * t + radii[3] * t * t * t;
        };
        const int n = 2000;
        double sum = 0.0;
        for (int i = 0; i <= n; ++i)
        {
            const double t = static_cast<double>(i) / n, rr = law(t);
            sum += (i == 0 || i == n ? 1.0 : (i % 2 == 1 ? 4.0 : 2.0)) * rr * rr;
        }
        const double removed = (1.0 - PI / 4.0) * 20.0 * sum / (3.0 * n);
        CHECK(smlib_fillet_variable(cube, edge, 0.1, from0, 3, knots, 8, radii, 4, &r) == SMLIB_OK);
        if (r != nullptr)
        {
            CHECK(near(8000.0 - volume_of(r), removed, 1e-4));
            done(r);
        }
    }

    // A variable round 2 -> r1 on the top edge x = 10, z = 30 of a 10 x 20 x 30 box. The top is
    // 10 wide: up to r1 = 10 the round fits, and takes (1 - pi/4) L (a^2 + ab + b^2) / 3. Past
    // it the contact line on the top leaves it across its far edge x = 0 before the round's end
    // (and past r1 = 30 the one on the side leaves it across its bottom edge): a convex round
    // overflowing cliffs, built on the faces' continuations and subtracted from the box, from
    // either end, as is a constant R12. The analytic volume integrates, along y, the section of
    // the round's removal region (the r x r square at the edge less the quarter disc about the
    // ball's centre (10 - r, 30 - r), SMLib's circular section in the planes y = const) inside
    // the box's section [0, 10] x [0, 30] (closed form in local coordinates a, b from the disc's
    // centre: a, b in [0, r], a^2 + b^2 >= r^2, a >= r - 10, b >= r - 30), by Gauss-Legendre in
    // pieces split where r passes 10 and 30. SMLib's variable round sections are cubics within
    // about 2e-5 of the radius: the variable rounds are checked to 1e-4, the constant to 1e-7.
    // (The corner without an end rib, for a short side edge, had built r1 = 20 and 30 across the
    // top's far edge as valid bodies; then they were refused, SMLIB_ERR_NO_ROOM.)
    {
        smlib_brep *slab = box(0, 0, 0, 10, 20, 30);
        const double edge[3] = {10, 10, 30};
        const double from0[3] = {10, 0, 30};
        const double from1[3] = {10, 20, 30};
        const double knots[4] = {0, 0, 1, 1};
        const double fits[2] = {2, 10};
        CHECK(smlib_fillet_variable(slab, edge, 0.1, from0, 1, knots, 4, fits, 2, &r) == SMLIB_OK);
        if (r != nullptr)
        {
            CHECK(near(6000.0 - volume_of(r), (1.0 - PI / 4.0) * 20.0 * (4.0 + 20.0 + 100.0) / 3.0, 1e-4));
            done(r);
        }
        // the integral of sqrt(r^2 - a^2) da
        auto disc = [](double rr, double a)
        {
            a = std::max(-rr, std::min(rr, a));
            return 0.5 * (a * std::sqrt(std::max(rr * rr - a * a, 0.0)) + rr * rr * std::asin(a / rr));
        };
        // the removal region's section inside the box's
        auto section = [&](double rr)
        {
            const double al = std::max(0.0, rr - 10.0), be = std::max(0.0, rr - 30.0);
            if (al >= rr || be >= rr)
            {
                return 0.0;
            }
            const double top = std::max(al, std::sqrt(std::max(rr * rr - be * be, 0.0)));
            return (rr - al) * rr - (disc(rr, top) - disc(rr, al) + be * (rr - top));
        };
        // 20-point Gauss-Legendre on [lo, hi]
        auto gauss = [](double lo, double hi, const auto &f)
        {
            static const double x[10] = {0.0765265211334973, 0.2277858511416451, 0.3737060887154195, 0.5108670019508271,
                                         0.6360536807265150, 0.7463319064601508, 0.8391169718222188, 0.9122344282513259,
                                         0.9639719272779138, 0.9931285991850949};
            static const double w[10] = {0.1527533871307258, 0.1491729864726037, 0.1420961093183820, 0.1316886384491766,
                                         0.1181945319615184, 0.1019301198172404, 0.0832767415767048, 0.0626720483341091,
                                         0.0406014298003869, 0.0176140071391521};
            const double m = 0.5 * (lo + hi), h = 0.5 * (hi - lo);
            double sum = 0.0;
            for (int i = 0; i < 10; ++i)
            {
                sum += w[i] * (f(m - h * x[i]) + f(m + h * x[i]));
            }
            return h * sum;
        };
        // the volume removed by the round going linearly r0 -> r1 along y in [0, 20]
        auto removed = [&](double r0, double r1)
        {
            auto at = [&](double y) { return section(r0 + (r1 - r0) * y / 20.0); };
            std::vector<double> cuts = {0.0, 20.0};
            for (double c : {10.0, 30.0})
            {
                const double y = (r1 != r0) ? 20.0 * (c - r0) / (r1 - r0) : -1.0;
                if (y > 0.0 && y < 20.0)
                {
                    cuts.push_back(y);
                }
            }
            std::sort(cuts.begin(), cuts.end());
            double sum = 0.0;
            for (size_t i = 0; i + 1 < cuts.size(); ++i)
            {
                // ten pieces per stretch between kinks
                for (int k = 0; k < 10; ++k)
                {
                    const double lo = cuts[i] + (cuts[i + 1] - cuts[i]) * k / 10.0, hi = cuts[i] + (cuts[i + 1] - cuts[i]) * (k + 1) / 10.0;
                    sum += gauss(lo, hi, at);
                }
            }
            return sum;
        };
        for (double r1 : {11.0, 20.0, 30.0, 50.0})
        {
            const double law[2] = {2, r1};
            const double back[2] = {r1, 2};
            const double want = 6000.0 - removed(2.0, r1);
            for (int from = 0; from < 2; ++from)
            {
                const int32_t s = from == 0 ? smlib_fillet_variable(slab, edge, 0.1, from0, 1, knots, 4, law, 2, &r)
                                            : smlib_fillet_variable(slab, edge, 0.1, from1, 1, knots, 4, back, 2, &r);
                if (s != SMLIB_OK || r == nullptr || !near(volume_of(r), want, 1e-4))
                {
                    std::fprintf(stderr, "  variable round 2 -> %g from y = %d: status %d\n", r1, from * 20, static_cast<int>(s));
                    failures++;
                }
                if (r != nullptr)
                {
                    done(r);
                }
            }
        }
        // the constant R12: 5384.17828518
        CHECK(smlib_fillet(slab, edge, 1, 0.1, 12.0, &r, nullptr) == SMLIB_OK);
        if (r != nullptr)
        {
            CHECK(near(volume_of(r), 6000.0 - removed(12.0, 12.0), 1e-7));
            done(r);
        }
        smlib_brep_free(slab);
    }

    // A round wider than its faces, on the long top edge of a bar whose section is W along the
    // top by H down the side: the round's rail on each face lies r from the edge, past the
    // face's far edge where r > W or r > H, a cliff overflow built on the faces continued and
    // subtracted. Its section takes off, at distance a from the edge along the top (a <= r), the
    // depth min(H, r - sqrt(r^2 - (a - r)^2)) down the side: in closed form with
    // D(x) = (x sqrt(r^2 - x^2) + r^2 asin(x / r)) / 2. An 80 x 5 x 5 bar rounded r9, r10 and
    // r12 (r10 and r12 failed in the fillet surface, with no rail ends: the corners found where
    // a rail meets the end face's side edge only on the edge itself, or from a guess already on
    // its line past the end); a 40 x 20 x 10 box's edge on its 10 high side rounded r10.05,
    // r10.07 and r12 (r10.05 and r10.07 built an end rib off the fillet: the rail's end, just
    // past the side edge's end, was placed on that end).
    {
        auto removed = [](double rr, double w, double h)
        {
            const double l = std::min(w, rr);
            const double a0 = rr > h ? rr - std::sqrt(rr * rr - (rr - h) * (rr - h)) : 0.0;
            auto d = [rr](double x) { return 0.5 * (x * std::sqrt(std::max(rr * rr - x * x, 0.0)) + rr * rr * std::asin(std::max(-1.0, std::min(1.0, x / rr)))); };
            auto f = [&](double a) { return rr * a - d(a - rr); };
            return h * a0 + f(l) - f(a0);
        };
        smlib_brep *bar = box(0, 0, 0, 80, 5, 5);
        const double bar_edge[3] = {40, 0, 5};
        for (double radius : {9.0, 10.0, 12.0})
        {
            CHECK(smlib_fillet(bar, bar_edge, 1, 0.1, radius, &r, nullptr) == SMLIB_OK && sound(r));
            if (r != nullptr)
            {
                CHECK(near(volume_of(r), 80.0 * (25.0 - removed(radius, 5.0, 5.0)), 1e-7));
                done(r);
            }
        }
        smlib_brep_free(bar);
        smlib_brep *slab = box(0, 0, 0, 40, 20, 10);
        const double slab_edge[3] = {20, 0, 10};
        for (double radius : {10.05, 10.07, 12.0})
        {
            CHECK(smlib_fillet(slab, slab_edge, 1, 0.1, radius, &r, nullptr) == SMLIB_OK && sound(r));
            if (r != nullptr)
            {
                CHECK(near(volume_of(r), 40.0 * (200.0 - removed(radius, 20.0, 10.0)), 1e-7));
                done(r);
            }
        }
        smlib_brep_free(slab);
    }

    // A 20 cube less a hole of radius 5 down from its top, its axis 5 from the side y = 0, so
    // its wall touches that side along a line: the top's rim and the side's top edge meet at the
    // touching point, where they run on tangentially into each other, a vertex no blend corner
    // handles. Chamfered 3 along the rim and that edge, the blend is refused at that corner,
    // SMLIB_ERR_CORNER_UNSUPPORTED (it failed with SMLib's generic error, its reason not kept).
    {
        smlib_brep *cube = box(0, 0, 0, 20, 20, 20);
        smlib_brep *rod = cylinder(10, 5, -1, 5, 22);
        smlib_brep *holed = nullptr;
        CHECK(smlib_boolean(cube, rod, SMLIB_BOOLEAN_SUBTRACT, &holed) == SMLIB_OK);
        const double rim[12] = {14.426091400901679, 2.7004358211515118, 20, 5.5730421003354138, 7.2997714423711244, 20, 5, 0, 20, 15, 0, 20};
        const double sizes[4] = {3, 3, 3, 3};
        CHECK(holed != nullptr && smlib_blend(holed, SMLIB_BLEND_CHAMFER, rim, sizes, 4, 0.5, &r, nullptr) == SMLIB_ERR_CORNER_UNSUPPORTED && r == nullptr);
        smlib_brep_free(holed);
        smlib_brep_free(rod);
        smlib_brep_free(cube);
    }

    // A chamfer of setback 2 on the top edge y = 0 of a prism 20 high on (0,0) (40,0)
    // (40 - 20s, 20) (20s, 20): its end faces slope by s, so the chamfer's contact line on the
    // top leaves it across an end face's top edge (the top's side edge at that end) before the
    // cross section at the edge's end. That is past the end, where the end face cuts the
    // chamfer off, not a blend too large for the top: the triangle 2 x 2 along 40 - 4s/3 comes
    // off.
    for (double s : {-0.5, 0.5})
    {
        const double xy[8] = {0, 0, 40, 0, 40 - 20 * s, 20, 20 * s, 20};
        int64_t kinds[4];
        double reals[24];
        for (int i = 0; i < 4; ++i)
        {
            const int j = (i + 1) % 4;
            const double seg[6] = {xy[2 * i], xy[2 * i + 1], 0, xy[2 * j], xy[2 * j + 1], 0};
            kinds[i] = SMLIB_CURVE_LINE;
            std::copy(seg, seg + 6, reals + 6 * i);
        }
        const double up[3] = {0, 0, 1};
        smlib_brep *prism = nullptr;
        CHECK(smlib_extrude(kinds, 4, reals, 24, up, 20.0, 0.0, SMLIB_DRAFT_CORNER_CURVATURE, 1, &prism, nullptr) == SMLIB_OK);
        const double edge[3] = {20, 0, 20};
        if (prism != nullptr && smlib_chamfer(prism, edge, 1, 0.1, 2.0, &r, nullptr) == SMLIB_OK)
        {
            CHECK(near(volume_of(prism) - volume_of(r), 2.0 * (40.0 - 4.0 * s / 3.0), 1e-8));
            done(r);
        }
        else
        {
            CHECK(!"a chamfer between sloped end faces");
        }
        smlib_brep_free(prism);
    }

    // A 100 x 60 x 30 slab less a basin lofted from 20..80 x 4..50 at the top to 30..70 x 20..40
    // at z = 10, its rim edge y = 4 chamfered 2. Its contact line on the top runs on over the top
    // past the rim's concave corners and leaves it at the slab's sides x = 0 and 100, beyond the
    // cross sections at the edge's ends; the one on the sloped wall leaves it across the wall's
    // sloped corner edges, before them. Both are past the ends. The section's drop below the top
    // is 2 / sqrt(1 + 0.8^2) and it comes off along 60 - drop / 3.
    {
        smlib_brep *slab = box(0, 0, 0, 100, 60, 30);
        const double rims[2][5] = {{20, 4, 80, 50, 30}, {30, 20, 70, 40, 10}};
        int64_t kinds[8];
        double reals[48];
        for (int s = 0; s < 2; ++s)
        {
            const double *q = rims[s];
            const double xy[4][2] = {{q[0], q[1]}, {q[2], q[1]}, {q[2], q[3]}, {q[0], q[3]}};
            for (int i = 0; i < 4; ++i)
            {
                const int j = (i + 1) % 4;
                const double seg[6] = {xy[i][0], xy[i][1], q[4], xy[j][0], xy[j][1], q[4]};
                kinds[4 * s + i] = SMLIB_CURVE_LINE;
                std::copy(seg, seg + 6, reals + 6 * (4 * s + i));
            }
        }
        const int64_t sizes[2] = {4, 4};
        smlib_brep *pocket = nullptr;
        smlib_brep *basin = nullptr;
        CHECK(smlib_loft(kinds, 8, reals, 48, sizes, 2, 1, 1, &pocket) == SMLIB_OK);
        CHECK(pocket != nullptr && smlib_boolean(slab, pocket, SMLIB_BOOLEAN_SUBTRACT, &basin) == SMLIB_OK);
        const double rim[3] = {50, 4, 30};
        const double drop = 2.0 / std::sqrt(1.0 + 0.64);
        if (basin != nullptr && smlib_chamfer(basin, rim, 1, 0.1, 2.0, &r, nullptr) == SMLIB_OK)
        {
            CHECK(near(volume_of(r), volume_of(basin) - drop * (60.0 - drop / 3.0), 1e-8));
            done(r);
        }
        else
        {
            CHECK(!"a chamfer on the sloped basin's rim");
        }
        for (smlib_brep *b : {slab, pocket, basin})
        {
            smlib_brep_free(b);
        }
    }

    // SolidCraft ladder 07-fillets/016: a 60 x 40 x 20 box, its top back edge rounded R5, its top
    // right edge R3, then the vertical edge x = 60, y = 40 (z 0..15). The R3 round leaves the R5
    // round's end arc on the x = 60 face 0.46 long; the vertical round's contact line there rolls
    // over the R3 round's tangent edge past the round's end and never meets the arc. The corner
    // has no end rib and the round is cut off where it meets the body: removed, outside the
    // vertical round, the integral of min(R5, R3 top) over the corner column.
    {
        smlib_brep *block = box(0, 0, 0, 60, 40, 20);
        const double back_top[3] = {30, 40, 20};
        const double right_top[3] = {60, 20, 20};
        const double vertical[3] = {60, 40, 7.5};
        smlib_brep *r5 = nullptr;
        smlib_brep *r3 = nullptr;
        CHECK(smlib_fillet(block, back_top, 1, 0.1, 5.0, &r5, nullptr) == SMLIB_OK);
        CHECK(r5 != nullptr && smlib_fillet(r5, right_top, 1, 0.1, 3.0, &r3, nullptr) == SMLIB_OK);
        const double radii[3] = {7.0, 4.0, 2.0};
        const double volumes[3] = {47418.926881, 47546.786538, 47591.811067};
        for (int i = 0; i < 3 && r3 != nullptr; ++i)
        {
            if (smlib_fillet(r3, vertical, 1, 0.1, radii[i], &r, nullptr) == SMLIB_OK)
            {
                CHECK(near(volume_of(r), volumes[i], 1e-8));
                done(r);
            }
            else
            {
                CHECK(!"ladder 016's vertical round under the short R5 end arc");
            }
        }
        for (smlib_brep *b : {block, r5, r3})
        {
            smlib_brep_free(b);
        }
    }

    // tangency-labs/solidcore#423: a sphere r2 in six patches (each on its own sphere surface),
    // cut flat at z = 1, its rim rounded R1. The clean-up after filleting removed every
    // topological edge in the body, the seams between the sphere's own patches too, and merged
    // the five patches into one face that no patch's parameterization could carry: the result
    // failed the check (loop orientation, rails half a millimetre off). It now leaves the
    // target's own faces apart. The round's section is the square 1 x 1 less the quarter disc,
    // centred on the circle of radius 1 at z = 0: 2 pi - pi^2 / 2 comes off the 9 pi.
    {
        smlib_brep *ball = cube_sphere(2.0);
        CHECK(ball != nullptr && sound(ball) && near(volume_of(ball), 32.0 * PI / 3.0, 1e-6));
        smlib_brep *lid = box(-3, -3, 1, 6, 6, 3);
        smlib_brep *cut = nullptr;
        CHECK(smlib_boolean(ball, lid, SMLIB_BOOLEAN_SUBTRACT, &cut) == SMLIB_OK && sound(cut));
        // The rim's four arcs, one on each side patch.
        const double q = std::sqrt(3.0);
        const double rim[12] = {q, 0, 1, 0, q, 1, -q, 0, 1, 0, -q, 1};
        if (cut != nullptr && smlib_fillet(cut, rim, 4, 0.5, 1.0, &r, nullptr) == SMLIB_OK)
        {
            CHECK(near(volume_of(r), 7.0 * PI + PI * PI / 2.0, 5e-5));
            done(r);
        }
        else
        {
            CHECK(!"the rim of the cut sphere in patches rounded R1");
        }
        for (smlib_brep *b : {ball, lid, cut})
        {
            smlib_brep_free(b);
        }
    }

    // tangency-labs/solidcore#422: a plate 40 x 20 x 10 with a through hole r3 at (20, 14.5),
    // its top edge y = 20 rounded R4. The round's contact line on the top (y = 16) runs across
    // the hole: the round is built as if the top went on over the hole, and what it removes is
    // subtracted, so the hole breaks into the round and nothing bridges it (rig case
    // 324-fillet-over-hole-overlap): 8000 - 90 pi - 40 (16 - 4 pi) + the overlap counted twice,
    // the integral over y in [16, 17.5] of (4 - sqrt(16 - (y - 16)^2)) 2 sqrt(9 - (y - 14.5)^2),
    // 0.3758. (Then it was refused, SMLIB_ERR_NO_ROOM; before that built across the hole.)
    {
        smlib_brep *plate = box(0, 0, 0, 40, 20, 10);
        smlib_brep *bore = cylinder(20, 14.5, -1, 3, 12);
        smlib_brep *holed_plate = nullptr;
        CHECK(smlib_boolean(plate, bore, SMLIB_BOOLEAN_SUBTRACT, &holed_plate) == SMLIB_OK);
        const double edge[3] = {5, 20, 10};
        CHECK(holed_plate != nullptr && smlib_fillet(holed_plate, edge, 1, 0.1, 4.0, &r, nullptr) == SMLIB_OK);
        if (r != nullptr)
        {
            CHECK(near(volume_of(r), 7580.287315360674, 1e-7));
            done(r);
        }
        for (smlib_brep *b : {plate, bore, holed_plate})
        {
            smlib_brep_free(b);
        }
    }

    // tangency-labs/solidcraft#1080: a cube 20 with three square openings 10 x 10 through it,
    // 5 in from the faces; its top front edge (y = 0, z = 20) chamfered or rounded 6 or 9 runs
    // past the 5 mm between the edge and the openings. Each is built as if the faces went on
    // over the openings, then cut by them: the chamfer's triangle y + 20 - z <= d (the round's
    // corner square less the quarter disc about (r, 20 - r)) along x, less its parts inside the
    // openings (the x opening y, z in [5, 15] everywhere along x; the y and z openings z in
    // [5, 15], y in [5, 15] for x in [5, 15]), integrated exactly (the chamfer) or as the round's
    // section in closed form by pieces (Gauss-Legendre), from 4000 left by the openings: chamfer
    // 6 3650, chamfer 9 3350, round 6 3846.04504440, round 9 3676.80702897. (Then they were
    // refused, SMLIB_ERR_NO_ROOM; before that built across the openings.)
    {
        smlib_brep *cube = box(0, 0, 0, 20, 20, 20);
        for (const std::array<double, 6> &o : {std::array<double, 6>{-1, 5, 5, 22, 10, 10}, std::array<double, 6>{5, -1, 5, 10, 22, 10},
                                                std::array<double, 6>{5, 5, -1, 10, 10, 22}})
        {
            smlib_brep *bar = box(o[0], o[1], o[2], o[3], o[4], o[5]);
            smlib_brep *cut = nullptr;
            CHECK(smlib_boolean(cube, bar, SMLIB_BOOLEAN_SUBTRACT, &cut) == SMLIB_OK);
            smlib_brep_free(cube);
            smlib_brep_free(bar);
            cube = cut;
        }
        const double edge[3] = {10, 0, 20};
        const std::array<std::pair<double, double>, 2> chamfers = {{{6.0, 3650.0}, {9.0, 3350.0}}};
        for (const auto &c : chamfers)
        {
            CHECK(smlib_chamfer(cube, edge, 1, 0.1, c.first, &r, nullptr) == SMLIB_OK);
            if (r != nullptr)
            {
                CHECK(near(volume_of(r), c.second, 1e-7));
                done(r);
            }
        }
        const std::array<std::pair<double, double>, 2> rounds = {{{6.0, 3846.0450444039}, {9.0, 3676.8070289689}}};
        for (const auto &c : rounds)
        {
            CHECK(smlib_fillet(cube, edge, 1, 0.1, c.first, &r, nullptr) == SMLIB_OK);
            if (r != nullptr)
            {
                CHECK(near(volume_of(r), c.second, 1e-6));
                done(r);
            }
        }
        smlib_brep_free(cube);
    }

    // tangency-labs/solidcraft#3435: a cylinder r7 h48.2 with a boss r2.7 on its top, its top rim
    // rounded 6.9: the round's contact line on the top (r = 0.1) lies under the boss, inside the
    // body's material. Built and cut, the round would meet the boss's underside tangentially,
    // in a knife edge: refused (SMLIB_ERR_NO_ROOM).
    {
        smlib_brep *can = cylinder(0, 0, 0, 7, 48.2);
        smlib_brep *boss = cylinder(0, 0, 48.2, 2.7, 2);
        smlib_brep *bossed = nullptr;
        CHECK(smlib_boolean(can, boss, SMLIB_BOOLEAN_UNION, &bossed) == SMLIB_OK);
        const double at = 100.0 * PI / 180.0;
        const double rim[3] = {7.0 * std::cos(at), 7.0 * std::sin(at), 48.2};
        CHECK(bossed != nullptr && smlib_fillet(bossed, rim, 1, 0.1, 6.9, &r, nullptr) == SMLIB_ERR_NO_ROOM && r == nullptr);
        for (smlib_brep *b : {can, boss, bossed})
        {
            smlib_brep_free(b);
        }
    }

    // A concave round that overflows stays refused: an L of a floor 20 x 20 x 2 and a wall 2
    // thick and 10 high along x = 0, its inner edge (x = 2, z = 2) rounded 9: the contact line on
    // the wall (z = 11) leaves it across its top. A round that adds material, built past the
    // wall, would fill across what lies beyond it (SMLIB_ERR_NO_ROOM). R6 fits: the floor's
    // 800 and the wall's 2 x 8 x 20 above it, and the round's fill (1 - pi/4) 6^2 along 20.
    {
        smlib_brep *floor = box(0, 0, 0, 20, 20, 2);
        smlib_brep *wall = box(0, 0, 0, 2, 20, 10);
        smlib_brep *l = nullptr;
        CHECK(smlib_boolean(floor, wall, SMLIB_BOOLEAN_UNION, &l) == SMLIB_OK);
        const double inner[3] = {2, 10, 2};
        CHECK(l != nullptr && smlib_fillet(l, inner, 1, 0.1, 6.0, &r, nullptr) == SMLIB_OK);
        if (r != nullptr)
        {
            CHECK(near(volume_of(r), 800.0 + 320.0 + (1.0 - PI / 4.0) * 36.0 * 20.0, 1e-6));
            done(r);
        }
        CHECK(l != nullptr && smlib_fillet(l, inner, 1, 0.1, 9.0, &r, nullptr) == SMLIB_ERR_NO_ROOM && r == nullptr);
        for (smlib_brep *b : {floor, wall, l})
        {
            smlib_brep_free(b);
        }
    }

    // A cylinder r10 h20, its top rim chamfered or rounded as wide as the top: each blend's rail
    // on the top is the top's centre, so the blend takes the whole top. A chamfer's line meets
    // the top at an angle, so its cone would end in a point there, and is refused
    // (SMLIB_ERR_NO_ROOM) in every form: the equal chamfer, by smlib_blend too, two distances
    // with either face first, a distance and an angle (SolidCraft's
    // chamfer_sides_curved_tests::too_large_or_unsupported_is_refused). The round's arc is
    // tangent to the top at its rail, so the R10 round closes the end smoothly as a hemisphere
    // and is built: the side, the bottom and the spherical end, the cylinder's first 10 and the
    // half ball r10 (SolidCraft's rig case fillets/5236). 9.99 fits every form and takes the
    // triangle (the round: the square less the quarter disc) of its section revolved about the
    // axis (Pappus).
    {
        smlib_brep *can = cylinder(0, 0, 0, 10, 20);
        const double whole = PI * 100.0 * 20.0;
        const double at = 100.0 * PI / 180.0;
        const double top_rim[3] = {10.0 * std::cos(at), 10.0 * std::sin(at), 20.0};
        const double up[3] = {0, 0, 1};
        const double down[3] = {0, 0, -1};
        // The ring the triangle with corners (r, z) sweeps about the axis.
        auto ring = [](double r0, double z0, double r1, double z1, double r2, double z2)
        {
            double area = 0.5 * std::fabs((r1 - r0) * (z2 - z0) - (r2 - r0) * (z1 - z0));
            return 2.0 * PI * (r0 + r1 + r2) / 3.0 * area;
        };
        for (double d : {10.0, 9.99})
        {
            auto builds = [&](const char *what, int32_t s, double removed)
            {
                if (s != SMLIB_OK || !near(volume_of(r), whole - removed, 1e-7))
                {
                    std::fprintf(stderr, "  %s %g: status %d\n", what, d, static_cast<int>(s));
                    failures++;
                }
                done(r);
            };
            auto chamfer = [&](const char *what, int32_t s, double removed)
            {
                if (d < 10.0)
                {
                    builds(what, s, removed);
                }
                else if (s != SMLIB_ERR_NO_ROOM || r != nullptr)
                {
                    std::fprintf(stderr, "  %s %g: status %d, not refused for want of room\n", what, d, static_cast<int>(s));
                    failures++;
                    smlib_brep_free(r);
                    r = nullptr;
                }
            };
            const double equal = ring(10, 20, 10 - d, 20, 10, 20 - d);
            const double sides = ring(10, 20, 10 - d, 20, 10, 17);
            const double slope = std::tan(PI / 6.0);
            chamfer("chamfer", smlib_chamfer(can, top_rim, 1, 0.1, d, &r, nullptr), equal);
            chamfer("blend chamfer", smlib_blend(can, SMLIB_BLEND_CHAMFER, top_rim, &d, 1, 0.1, &r, nullptr), equal);
            chamfer("chamfer across the top, 3 down", smlib_chamfer_sides(can, top_rim, 1, 0.1, d, 3.0, 0, up, &r, nullptr), sides);
            chamfer("chamfer 3 down, across the top", smlib_chamfer_sides(can, top_rim, 1, 0.1, 3.0, d, 0, down, &r, nullptr), sides);
            chamfer("chamfer across the top at 30 degrees", smlib_chamfer_sides(can, top_rim, 1, 0.1, d, PI / 6.0, 1, up, &r, nullptr),
                    ring(10, 20, 10 - d, 20, 10, 20 - d * slope));
            chamfer("chamfer 3 down at 60 degrees", smlib_chamfer_sides(can, top_rim, 1, 0.1, d / std::tan(PI / 3.0), PI / 3.0, 1, down, &r, nullptr),
                    ring(10, 20, 10 - d, 20, 10, 20 - d / std::tan(PI / 3.0)));
            // The round's section: the square d x d at the rim less the quarter disc, its
            // centroid d (10 - 3 pi) / (12 - 3 pi) in from the rim. At 10 the body left is the
            // hemisphere-ended cylinder, pi 100 10 + (2/3) pi 1000.
            int32_t s = smlib_fillet(can, top_rim, 1, 0.1, d, &r, nullptr);
            if (d == 10.0 && r != nullptr)
            {
                CHECK(near(volume_of(r), PI * 100.0 * 10.0 + 2.0 / 3.0 * PI * 1000.0, 1e-7));
                CHECK(r->p->GetNumFaces() == 3);
            }
            builds("fillet", s, 2.0 * PI * (10.0 - d * (10.0 - 3.0 * PI) / (12.0 - 3.0 * PI)) * (1.0 - PI / 4.0) * d * d);
        }
        smlib_brep_free(can);
    }

    // A round wider than the rim it rounds: a ball of r2 rolling inside a pin of radius R < 2,
    // tangent to its side and its top, has its centre on the far side of the axis, so no ball of
    // that size fits; the round's rail on the top (the rim offset in by 2, a circle of radius
    // R - 2 < 0) turns back against the rim and the round would be a spindle torus through
    // itself. Refused, SMLIB_ERR_NO_ROOM (SMLib built the folded round and failed in the merge
    // after it with its generic error). At R = 2 the ball's centre is on the axis: the
    // hemispherical end, pi 4 8 + (2/3) pi 8.
    {
        for (double R : {1.99, 1.5, 1.0})
        {
            smlib_brep *pin = cylinder(0, 0, 0, R, 10);
            const double rim[3] = {0, R, 10};
            smlib_brep *rounded = nullptr;
            CHECK(smlib_fillet(pin, rim, 1, 0.1, 2.0, &rounded, nullptr) == SMLIB_ERR_NO_ROOM && rounded == nullptr);
            smlib_brep_free(pin);
        }
        smlib_brep *pin = cylinder(0, 0, 0, 2, 10);
        const double rim[3] = {0, 2, 10};
        smlib_brep *rounded = nullptr;
        CHECK(smlib_fillet(pin, rim, 1, 0.1, 2.0, &rounded, nullptr) == SMLIB_OK);
        if (rounded != nullptr)
        {
            CHECK(near(volume_of(rounded), PI * 4.0 * 8.0 + 2.0 / 3.0 * PI * 8.0, 1e-7));
            done(rounded);
        }
        smlib_brep_free(pin);
    }

    // Both vertical edges of the 15 wide end of an 80 x 15 x 20 box rounded r: up to r = 7.5
    // the two rounds' bands on the end face stay apart (at 7.5 they meet, a full round, removing
    // 2 (1 - pi/4) r^2 20); past it each rail lies inside the other's band and the rounds are
    // larger than the face allows, SMLIB_ERR_NO_ROOM (r7.6 and r10 were built as bodies that are
    // not closed, SMLIB_ERR_NOT_SOLID).
    {
        smlib_brep *bar = box(0, 0, 0, 80, 15, 20);
        const double ends[6] = {80, 0, 10, 80, 15, 10};
        for (double radius : {7.4, 7.5})
        {
            CHECK(smlib_fillet(bar, ends, 2, 0.1, radius, &r, nullptr) == SMLIB_OK && sound(r));
            if (r != nullptr)
            {
                CHECK(near(volume_of(r), 24000.0 - 2.0 * (1.0 - PI / 4.0) * radius * radius * 20.0, 1e-9));
                done(r);
            }
        }
        for (double radius : {7.6, 10.0})
        {
            CHECK(smlib_fillet(bar, ends, 2, 0.1, radius, &r, nullptr) == SMLIB_ERR_NO_ROOM && r == nullptr);
        }
        smlib_brep_free(bar);
    }

    // A round run up a box's vertical edge into the corner where two smaller rounds of the top
    // meet (SolidCraft ladder 07-fillets/016): the 60 x 40 x 20 box, its top edge at y = 40
    // rounded r5, 60 25 (1 - pi/4) off, then its top edge at x = 60 rounded r3 (the r5 round
    // ends it at y = 35: between 35 and 40 lengths of 9 (1 - pi/4) off), then the vertical edge
    // at x = 60, y = 40 with the edges tangent to it up into the corner: the side's edge on the
    // r5 round (a circular arc of radius 5 about the round's axis, y = 35, z = 15) and the edge
    // between the two rounds. A round of r there rolls on the side and on the r5 round, its rail
    // on the side the arc offset in by r, a circle of radius 5 - r: past r = 5 it lies beyond
    // the axis and runs against the rail on the round (the round a spindle torus through
    // itself, its ball curving less than the r5 round around it). On the edge between the
    // rounds the ball rolls around the r3 round, which past r = 3 curves tighter than the ball
    // (the two rounds' offsets, of radii 5 - r and r - 3 about axes 2 apart, then meet only in
    // a point). No ball over r3 fits: r3.1, r4, r5 and r7 are SMLIB_ERR_NO_ROOM (r7 failed in
    // the corner's vertex, looking on the side for a rail end that is not there, with SMLib's
    // generic error, and so did r4 and r5); r1, r2, r2.5 and r2.9 build and take off at least
    // the 15 straight lengths of r^2 (1 - pi/4) below the top rounds and less than twice that.
    // (tangency-labs/solidcore#1173: r2.5 and r2.9 came out unsound here, after the earlier
    // blends of this run, and sound in a fresh process: the finite-difference step at a tangent
    // point was a static, set by the first point any blend of the process computed. Run as a
    // child first thing in a process and after a failed constant-distance blend, r2.5 must make
    // the same body byte for byte.)
    {
        const double quarter = 1.0 - PI / 4.0;
        smlib_brep *slab = box(0, 0, 0, 60, 40, 20);
        const double top_edge[3] = {30, 40, 20};
        const double end_edge[3] = {60, 20, 20};
        smlib_brep *top = nullptr;
        smlib_brep *both = nullptr;
        CHECK(smlib_fillet(slab, top_edge, 1, 0.1, 5.0, &top, nullptr) == SMLIB_OK && sound(top));
        CHECK(near(volume_of(top), 48000.0 - 60.0 * 25.0 * quarter, 1e-9));
        CHECK(smlib_fillet(top, end_edge, 1, 0.1, 3.0, &both, nullptr) == SMLIB_OK && sound(both));
        const double rounded = volume_of(both);
        CHECK(rounded < volume_of(top) - 35.0 * 9.0 * quarter && rounded > volume_of(top) - 40.0 * 9.0 * quarter);
        const double up[9] = {60, 40, 7.5, 60, 35.0 + std::sqrt(24.0), 16, 58.931891825700852, 37.55984948554854, 19.294016994391871};
        for (double radius : {3.1, 4.0, 5.0, 7.0})
        {
            smlib_brep *corner = nullptr;
            CHECK(smlib_fillet(both, up, 3, 0.1, radius, &corner, nullptr) == SMLIB_ERR_NO_ROOM && corner == nullptr);
        }
        for (double radius : {1.0, 2.0, 2.5, 2.9})
        {
            smlib_brep *corner = nullptr;
            CHECK(smlib_fillet(both, up, 3, 0.1, radius, &corner, nullptr) == SMLIB_OK);
            if (corner != nullptr)
            {
                const double off = rounded - volume_of(corner);
                CHECK(off > 15.0 * radius * radius * quarter && off < 30.0 * radius * radius * quarter);
                done(corner);
            }
        }
        smlib_brep_free(both);
        smlib_brep_free(top);
        smlib_brep_free(slab);
        bool first_sound = false;
        bool after_sound = false;
        const std::string first = corner_child_run(argv[0], "first", first_sound);
        const std::string after = corner_child_run(argv[0], "after", after_sound);
        std::printf("corner r2.5 in a fresh process: sound %d, after a failed chord blend: sound %d, same body %d\n", first_sound, after_sound,
                    !first.empty() && first == after);
        CHECK(first_sound && after_sound);
        CHECK(!first.empty() && first == after);
    }

    // tangency-labs/solidcore#419: an offset surface evaluated beyond its
    // base surface's domain extends the base, and the extension's shape depends on its length.
    // It was made to reach just past the first point evaluated there, and extended again from
    // itself for the next, so two copies of one offset surface evaluated in another order (the
    // fillet solver's, which finds a round's end points, and the trace's) gave points 0.03 mm
    // apart at the same parameters, and the trace refused the end points. Now the extension's
    // reach depends only on how far out the point lies. A half cylinder r5 offset by 4, beyond
    // both ends of its arc: once straight at two points, once after two nearer points first.
    {
        const SmContext &ctx = *SmApiGetOrCreateContext();
        SmCylinder *half = new (ctx) SmCylinder(SmPoint3d(0, 0, 0), SmVector3d(1, 0, 0), SmVector3d(0, 1, 0), 5.0, 0.0, 180.0, 10.0,
                                                FALSE, FALSE, FALSE, &ctx);
        SmOffsetSurface *first = new (ctx) SmOffsetSurface(4.0, *half);
        SmOffsetSurface *second = new (ctx) SmOffsetSurface(4.0, *half);
        first->SetAllowExtension(TRUE);
        second->SetAllowExtension(TRUE);
        SmExtent2d d = half->GetNaturalUVDomain();
        double length = d.GetMax().x - d.GetMin().x;
        double v = 0.5 * (d.GetMin().y + d.GetMax().y);
        SmPoint3d p;
        // inside the arc, before and after the extension: the same points (its pcurves hold)
        const double inside[3] = {d.GetMin().x + 0.1 * length, d.GetMin().x + 0.5 * length, d.GetMin().x + 0.93 * length};
        SmPoint3d before[3];
        for (int k = 0; k < 3; ++k)
        {
            CHECK(first->EvaluatePoint(SmPoint2d(inside[k], v), before[k]) == SM_SUCCESS);
        }
        const double nearer[2] = {d.GetMax().x + 0.067 * length, d.GetMin().x - 0.067 * length};
        for (double u : nearer)
        {
            CHECK(first->EvaluatePoint(SmPoint2d(u, v), p) == SM_SUCCESS);
        }
        for (int k = 0; k < 3; ++k)
        {
            CHECK(first->EvaluatePoint(SmPoint2d(inside[k], v), p) == SM_SUCCESS && p.DistanceBetween(before[k]) < 1e-12);
        }
        const double farther[2] = {d.GetMax().x + 0.116 * length, d.GetMin().x - 0.116 * length};
        for (double u : farther)
        {
            SmPoint3d a;
            SmPoint3d b;
            CHECK(second->EvaluatePoint(SmPoint2d(u, v), b) == SM_SUCCESS);
            CHECK(first->EvaluatePoint(SmPoint2d(u, v), a) == SM_SUCCESS);
            CHECK(a.DistanceBetween(b) < 1e-9);
        }
        // The offset continues as the circle of radius 9 about the cylinder's axis (z), out to
        // three domains (a turn and a half) past either end of the arc.
        double worst = 0.0;
        for (double past : {0.07, 0.5, 1.0, 2.0, 3.0})
        {
            for (double u : {d.GetMax().x + past * length, d.GetMin().x - past * length})
            {
                SmOffsetSurface *fresh = new (ctx) SmOffsetSurface(4.0, *half);
                fresh->SetAllowExtension(TRUE);
                CHECK(fresh->EvaluatePoint(SmPoint2d(u, v), p) == SM_SUCCESS);
                worst = std::max(worst, std::fabs(std::hypot(p.x, p.y) - 9.0));
                delete fresh;
            }
        }
        std::printf("half cylinder r5 offset 4, 0.07 to 3 domains past its ends: off radius 9 by at most %.3g\n", worst);
        CHECK(worst < 1e-9);
        // A quarter turn's span per half domain: half a domain past the end at 180 degrees is
        // at 270, a domain past it at 360.
        const double quarter[2][3] = {{0.5, 0.0, -9.0}, {1.0, 9.0, 0.0}};
        for (const auto &q : quarter)
        {
            SmOffsetSurface *fresh = new (ctx) SmOffsetSurface(4.0, *half);
            fresh->SetAllowExtension(TRUE);
            CHECK(fresh->EvaluatePoint(SmPoint2d(d.GetMax().x + q[0] * length, v), p) == SM_SUCCESS);
            std::printf("  %g domain past the end: (%.9f, %.9f, %.9f)\n", q[0], p.x, p.y, p.z);
            CHECK(std::fabs(p.x - q[1]) < 1e-9 && std::fabs(p.y - q[2]) < 1e-9);
            delete fresh;
        }
        delete first;
        delete second;
        delete half;
    }

    // A 294 deg partial cylinder r5 offset 4, evaluated past its ends (the round solver's first
    // requests are 18 and 4.5 deg past them). The extension reached a domain past either end, 294
    // deg each way, and covered the turn two and a half times over. Now it closes the turn, 33 deg
    // past either end (half the gap), and wraps parameters past that like a periodic surface's.
    {
        const SmContext &ctx = *SmApiGetOrCreateContext();
        SmCylinder *arc = new (ctx) SmCylinder(SmPoint3d(0, 0, 0), SmVector3d(1, 0, 0), SmVector3d(0, 1, 0), 5.0, 0.0, 294.0, 10.0,
                                               FALSE, FALSE, FALSE, &ctx);
        const SmExtent2d d = arc->GetNaturalUVDomain();
        const double length = d.GetMax().x - d.GetMin().x;
        const double per_deg = length / 294.0;
        const double turn = 360.0 * per_deg;
        const double v = 0.5 * (d.GetMin().y + d.GetMax().y);
        SmOffsetSurface *off = new (ctx) SmOffsetSurface(4.0, *arc);
        off->SetAllowExtension(TRUE);
        SmPoint3d p;
        SmPoint3d inside;
        CHECK(off->EvaluatePoint(SmPoint2d(d.GetMin().x + 0.4 * length, v), inside) == SM_SUCCESS);
        double worst = 0.0;
        for (double past : {18.0, 4.5, 33.0, 40.0, 147.0, 294.0, 588.0, 900.0})
        {
            for (double u : {d.GetMax().x + past * per_deg, d.GetMin().x - past * per_deg})
            {
                SmPoint3d q;
                CHECK(off->EvaluatePoint(SmPoint2d(u, v), p) == SM_SUCCESS);
                CHECK(off->EvaluatePoint(SmPoint2d(u - turn, v), q) == SM_SUCCESS);
                worst = std::max({worst, std::fabs(std::hypot(p.x, p.y) - 9.0), p.DistanceBetween(q)});
            }
        }
        const SmExtent2d e = off->GetNaturalUVDomain();
        std::printf("294 deg cylinder r5 offset 4, 4.5 to 900 deg past its ends: off radius 9 or a turn back by at most %.3g, "
                    "extended over %.12g deg\n", worst, (e.GetMax().x - e.GetMin().x) / per_deg);
        CHECK(worst < 1e-9);
        // one turn, the gap's half past either end
        CHECK(std::fabs(e.GetMax().x - e.GetMin().x - turn) < 1e-9 * turn);
        CHECK(std::fabs(e.GetMax().x - d.GetMax().x - 33.0 * per_deg) < 1e-9 * turn);
        // the closure, at -33 deg from either side; inside the arc the same point as before
        for (double u : {d.GetMax().x + 33.0 * per_deg, d.GetMin().x - 33.0 * per_deg})
        {
            CHECK(off->EvaluatePoint(SmPoint2d(u, v), p) == SM_SUCCESS);
            const double a = -33.0 * SM_PI / 180.0;
            CHECK(std::fabs(p.x - 9.0 * std::cos(a)) < 1e-9 && std::fabs(p.y - 9.0 * std::sin(a)) < 1e-9);
        }
        CHECK(off->EvaluatePoint(SmPoint2d(d.GetMin().x + 0.4 * length, v), p) == SM_SUCCESS && p.DistanceBetween(inside) < 1e-12);
        delete off;
        delete arc;
    }

    // SolidCraft rig case fillets/4174: a block 40 x 30 x 5 and a cylinder r10.5 standing on its
    // top, its axis on the block's front edge, so half of it overhangs. The concave arc where the
    // cylinder meets the top, rounded r, ends at the front face, where the arc runs on tangentially
    // as the convex rim of the overhang's underside. The rail on the cylinder runs alongside that
    // rim and never meets it; the corner took the rim for the side edge the rail ends on (1001).
    // The round ends on the front face's extension: Pappus over the half turn, the section
    // (1 - pi/4) r^2 with its centroid r (5/6 - pi/4) / (1 - pi/4) from the cylinder.
    {
        smlib_brep *block = box(-20, 0, 0, 40, 30, 5);
        smlib_brep *post = cylinder(0, 0, 5, 10.5, 30);
        smlib_brep *joined = nullptr;
        CHECK(smlib_boolean(block, post, SMLIB_BOOLEAN_UNION, &joined) == SMLIB_OK && sound(joined));
        const double base_arc[6] = {10.5 * std::cos(PI / 4), 10.5 * std::sin(PI / 4), 5, -10.5 * std::cos(PI / 4), 10.5 * std::sin(PI / 4), 5};
        for (double radius : {1.0, 2.0, 3.0})
        {
            const double section = (1.0 - PI / 4.0) * radius * radius;
            const double centroid = radius * (5.0 / 6.0 - PI / 4.0) / (1.0 - PI / 4.0);
            if (joined != nullptr && smlib_fillet(joined, base_arc, 2, 0.1, radius, &r, nullptr) == SMLIB_OK)
            {
                CHECK(near(volume_of(r), volume_of(joined) + PI * (10.5 + centroid) * section, 1e-8));
                done(r);
            }
            else
            {
                CHECK(!"the cylinder's base arc rounded up to the front face");
            }
        }
        for (smlib_brep *b : {block, post, joined})
        {
            smlib_brep_free(b);
        }
    }

    // SolidCraft rig case fillets/4369: a 60 x 40 x 10 block rounded r10 in plan and a boss r10 x 10
    // on its top at (12, 32), overhanging the block's corner. The boss's base arc over the block,
    // rounded R2, ends at the block's corner wall and back wall, where the arc runs on as the rim
    // of the boss's underside. On the corner wall the top's rail meets the wall 42 degrees round
    // the boss before the arc's end, the boss's rail only at it: the round's surface, made to
    // span its rail ends found then, stopped short of the arc's end, and the cap's rib up the
    // wall ran off it (1001). The round adds the material within 2 of the boss's foot, under the
    // quarter disc, over the block: the integral of 2 - sqrt(4 - (12 - rho)^2) over the plan
    // annulus 10 < rho < 12 about the boss's axis inside the block's outline.
    {
        const int64_t kinds[8] = {SMLIB_CURVE_LINE, SMLIB_CURVE_ARC, SMLIB_CURVE_LINE, SMLIB_CURVE_ARC,
                                  SMLIB_CURVE_LINE, SMLIB_CURVE_ARC, SMLIB_CURVE_LINE, SMLIB_CURVE_ARC};
        const double reals[72] = {10, 0,  0, 50, 0,  0, 50, 10, 0, 0, 0, 1, 1, 0, 0, 10, -PI / 2, 0,
                                  60, 10, 0, 60, 30, 0, 50, 30, 0, 0, 0, 1, 1, 0, 0, 10, 0,       PI / 2,
                                  50, 40, 0, 10, 40, 0, 10, 30, 0, 0, 0, 1, 1, 0, 0, 10, PI / 2,  PI,
                                  0,  30, 0, 0,  10, 0, 10, 10, 0, 0, 0, 1, 1, 0, 0, 10, PI,      3 * PI / 2};
        const double up[3] = {0, 0, 1};
        smlib_brep *plate = nullptr;
        CHECK(smlib_extrude(kinds, 8, reals, 72, up, 10.0, 0.0, SMLIB_DRAFT_CORNER_TANGENT, 1, &plate, nullptr) == SMLIB_OK);
        smlib_brep *boss = cylinder(12, 32, 10, 10, 10);
        smlib_brep *joined = nullptr;
        CHECK(plate != nullptr && smlib_boolean(plate, boss, SMLIB_BOOLEAN_UNION, &joined) == SMLIB_OK && sound(joined));
        auto in_outline = [](double x, double y) {
            double qx = std::max({10.0 - x, 0.0, x - 50.0});
            double qy = std::max({10.0 - y, 0.0, y - 30.0});
            return x >= 0 && x <= 60 && y >= 0 && y <= 40 && qx * qx + qy * qy <= 100.0;
        };
        double added = 0.0;
        const int n_rho = 1000;
        const int n_turn = 8000;
        for (int i = 0; i < n_rho; ++i)
        {
            double rho = 10.0 + (i + 0.5) * 2.0 / n_rho;
            double height = 2.0 - std::sqrt(4.0 - (12.0 - rho) * (12.0 - rho));
            for (int j = 0; j < n_turn; ++j)
            {
                double t = (j + 0.5) * 2.0 * PI / n_turn;
                if (in_outline(12.0 + rho * std::cos(t), 32.0 + rho * std::sin(t)))
                {
                    added += height * rho * (2.0 / n_rho) * (2.0 * PI / n_turn);
                }
            }
        }
        const double foot[6] = {12, 22, 10, 12 + 10 * std::cos(PI / 6), 32 + 10 * std::sin(PI / 6), 10};
        if (joined != nullptr && smlib_fillet(joined, foot, 2, 0.1, 2.0, &r, nullptr) == SMLIB_OK)
        {
            CHECK(near(volume_of(r), volume_of(joined) + added, 1e-7));
            done(r);
        }
        else
        {
            CHECK(!"the boss's base arc rounded up to the block's walls");
        }
        for (smlib_brep *b : {plate, boss, joined})
        {
            smlib_brep_free(b);
        }
    }

    // SolidCraft rig case fillets/1207 (tests/data/bar-boss-over-sides.txt, as SolidCraft hands it
    // to smlib_blend): a bar 100 x 30 x 10 (y -15..15) joined with a boss of diameter 40, 20 tall,
    // at the origin (seamed at 0.6131 rad), standing out over both of the bar's sides. The bar's
    // front top edge and the boss's base arc over the bar (two edges, split at the seam) rounded
    // together meet at the vertex where the top, the front face and the boss meet, whose third
    // edge, the boss's line on the front face, is concave and runs down from it. The 3x2 mixed
    // corner there took the concave round's rail to cross that line, as in an L-shaped box, where
    // the third edge is convex and runs up; it crosses only the line's continuation, and the corner
    // patch's boundary ran off to the far side of the boss (1001). Away from the corner the body
    // is the bar and boss, plus the concave round's spandrel over the bar (strips ending on the
    // side planes), less the convex round's spandrel out to the boss's wall; the corner patch lies
    // within the cube of side 2 r about the vertex.
    {
        smlib_flat *flat = nullptr;
        smlib_brep *joined = nullptr;
        CHECK(smlib_flat_read(SMLIB_TEST_DATA "/bar-boss-over-sides.txt", &flat) == SMLIB_OK && smlib_brep_import(flat, &joined) == SMLIB_OK);
        smlib_flat_free(flat);
        for (const double rad : {1.0, 2.0})
        {
            const int n = 2000;
            double added = 0.0;
            for (int i = 0; i < n; ++i)
            {
                double rho = 20.0 + (i + 0.5) * rad / n;
                added += (rad - std::sqrt(rad * rad - (20.0 + rad - rho) * (20.0 + rad - rho))) * rho * 2.0 * std::asin(15.0 / rho) * rad / n;
            }
            double removed = 0.0;
            const double h = rad / n;
            for (int i = 0; i < n; ++i)
            {
                double y = -15.0 + (i + 0.5) * h;
                removed += (rad - std::sqrt(rad * rad - (y + 15.0 - rad) * (y + 15.0 - rad))) * (100.0 - std::sqrt(400.0 - y * y)) * h;
            }
            const double edges[9] = {60, -15, 10, 20 * std::cos(0.73), 20 * std::sin(0.73), 10, 20 * std::cos(-0.12), 20 * std::sin(-0.12), 10};
            const double sizes[3] = {rad, rad, rad};
            if (joined != nullptr && smlib_blend(joined, SMLIB_BLEND_ROUND, edges, sizes, 3, 0.1, &r, nullptr) == SMLIB_OK)
            {
                double want = volume_of(joined) + added - removed;
                CHECK(near(volume_of(r), want, 8.0 * rad * rad * rad / want));
                // The corner patch is tangent to both rounds and to the boss's side along all its
                // edges. Its Coons construction averaged the corners' twists, which each round's
                // curvature fixes differently, and left it 15.6 degrees off the bar's round at r 2.
                int patches = 0;
                const double angle = worst_corner_seam_angle(r, SmPoint3d(std::sqrt(175.0), -15, 10), 3.0 * rad, patches);
                CHECK(patches == 1 && angle < 1.0);
                if (!(angle < 1.0))
                {
                    std::fprintf(stderr, "  r %g: corner patch seams %.3g degrees\n", rad, angle);
                }
                done(r);
            }
            else
            {
                CHECK(!"the bar's top edge and the boss's base arc rounded together");
            }
        }
        smlib_brep_free(joined);
    }

    // The same 3x2 mixed corner in an L-shaped box, its third edge convex (the riser's front
    // edge): the 60 x 40 x 20 box less the notch (40, 20, 10)-(70, 50, 30), the notch floor's
    // concave edge and its convex side edge rounded r 2 and 5. Its patch was 14.7 degrees off the
    // convex round, and 38 degrees off the side face near the end of its rail.
    {
        smlib_brep *part = box(0, 0, 0, 60, 40, 20);
        smlib_brep *block = box(40, 20, 10, 30, 30, 20);
        smlib_brep *notched = nullptr;
        CHECK(smlib_boolean(part, block, SMLIB_BOOLEAN_SUBTRACT, &notched) == SMLIB_OK && notched != nullptr);
        smlib_brep_free(part);
        smlib_brep_free(block);
        for (const double rad : {2.0, 5.0})
        {
            const double edges[6] = {50, 20, 10, 60, 30, 10};
            const double sizes[2] = {rad, rad};
            if (notched != nullptr && smlib_blend(notched, SMLIB_BLEND_ROUND, edges, sizes, 2, 0.1, &r, nullptr) == SMLIB_OK)
            {
                int patches = 0;
                const double angle = worst_corner_seam_angle(r, SmPoint3d(60, 20, 10), 2.0 * rad, patches);
                CHECK(patches == 1 && angle < 1.0);
                if (!(angle < 1.0))
                {
                    std::fprintf(stderr, "  L-shaped box r %g: corner patch seams %.3g degrees\n", rad, angle);
                }
                done(r);
            }
            else
            {
                CHECK(!"an L-shaped box's notch floor edges rounded together");
            }
        }
        smlib_brep_free(notched);
    }

    // SolidCraft ladder 07-fillets/012: a variable round, r0 at x = 0 to r1 at x = 55, on the top
    // back edge of a 60 x 40 x 20 box whose vertical edge at (60, 40) is rounded R5. With r1 = 5
    // the rail on the top meets the R5 corner's arc only at its end (60, 35), where the law,
    // extended past the edge's end, is larger: there it leaves the top across the right face's
    // edge. The second search for that rail's end held the radius at r1 even though the corner
    // is not calculated as constant radius (the first search holds it only then), found the
    // arc's end, and the variable round's trace could not pass through it (1001). Removed: the
    // section (1 - pi/4) r^2 along 55, and past x = 55 (s along it, the law going on as
    // r = r1 + (r1 - r0) s / 55) the square r x r less the quarter disc as far as the R5
    // cylinder reaches, y - 35 < sqrt(25 - s^2).
    {
        smlib_brep *block = box(0, 0, 0, 60, 40, 20);
        const double corner_edge[3] = {60, 40, 10};
        smlib_brep *rounded = nullptr;
        CHECK(smlib_fillet(block, corner_edge, 1, 0.1, 5.0, &rounded, nullptr) == SMLIB_OK);
        const double back_top[3] = {27.5, 40, 20};
        const double from[3] = {0, 40, 20};
        const double knots[4] = {0, 0, 1, 1};
        auto under_arc = [](double rr, double a) { // the integral of rr - sqrt(rr^2 - w^2) over 0 < w < a
            a = std::max(0.0, std::min(a, rr));
            return rr * a - (0.5 * a * std::sqrt(rr * rr - a * a) + 0.5 * rr * rr * std::asin(a / rr));
        };
        for (double r0 : {1.0, 2.0})
        {
            const double r1 = 5.0;
            const double law[2] = {r0, r1};
            double removed = (1.0 - PI / 4.0) * 55.0 * (r1 * r1 * r1 - r0 * r0 * r0) / (3.0 * (r1 - r0));
            const int steps = 200000;
            for (int i = 0; i < steps; ++i)
            {
                double s = (i + 0.5) * 5.0 / steps;
                double rr = r1 + (r1 - r0) * s / 55.0;
                removed += under_arc(rr, std::sqrt(25.0 - s * s) - 5.0 + rr) * 5.0 / steps;
            }
            if (rounded != nullptr && smlib_fillet_variable(rounded, back_top, 0.1, from, 1, knots, 4, law, 2, &r) == SMLIB_OK)
            {
                CHECK(near(volume_of(r), volume_of(rounded) - removed, 1e-6));
                done(r);
            }
            else
            {
                CHECK(!"a variable round reaching R5 at the R5 corner");
            }
        }
        smlib_brep_free(rounded);
        smlib_brep_free(block);
    }

    // Booleans that split, vanish or touch.
    smlib_brep *cut = box(8, -1, -1, 4, 22, 22);
    CHECK(smlib_boolean(cube, cut, SMLIB_BOOLEAN_SUBTRACT, &r) == SMLIB_OK);
    done(r);
    smlib_brep *big = box(-5, -5, -5, 30, 30, 30);
    CHECK(smlib_boolean(cube, big, SMLIB_BOOLEAN_SUBTRACT, &r) == SMLIB_ERR_EMPTY);
    smlib_brep *far = box(50, 0, 0, 10, 10, 10);
    CHECK(smlib_boolean(cube, far, SMLIB_BOOLEAN_UNION, &r) == SMLIB_OK);
    done(r);
    smlib_brep *touching = box(20, 0, 0, 10, 20, 20);
    if (smlib_boolean(cube, touching, SMLIB_BOOLEAN_UNION, &r) == SMLIB_OK)
    {
        CHECK(near(volume_of(r), 12000.0, 1e-9));
        done(r);
    }

    // History through the boolean shim: tagged operands in, the result's history out. A slot
    // across the top of the cube: every result face comes from a face of one operand, the
    // slot's walls and floor from the tool (operand 1), the rest from the cube (operand 0).
    {
        smlib_brep *slot = box(5, -1, 10, 10, 22, 15);
        smlib_brep *ta = nullptr;
        smlib_brep *tb = nullptr;
        CHECK(smlib_history_tag(cube, 0, &ta) == SMLIB_OK);
        CHECK(smlib_history_tag(slot, 1, &tb) == SMLIB_OK);
        CHECK(smlib_boolean(ta, tb, SMLIB_BOOLEAN_SUBTRACT, &r) == SMLIB_OK);
        int64_t counts[3] = {0, 0, 0};
        CHECK(smlib_brep_counts(r, counts) == SMLIB_OK && counts[0] == 10);
        int64_t len = 0;
        CHECK(smlib_history_read(r, nullptr, 0, &len) == SMLIB_ERR_BUFFER && len > 3);
        std::vector<int64_t> h(static_cast<size_t>(len));
        CHECK(smlib_history_read(r, h.data(), len, &len) == SMLIB_OK);
        CHECK(h[0] == SMLIB_HISTORY_VERSION && h[1] == counts[0] && h[2] == counts[1]);
        int from[2] = {0, 0};
        size_t at = 3;
        for (int64_t f = 0; f < h[1] && at + 1 < h.size(); ++f)
        {
            int64_t n = h[at + 1];
            CHECK(n >= 1 && h[at] != SMLIB_CHANGE_GENERATED);
            for (int64_t i = 0; i < n && at + 3 + 2 * i < h.size(); ++i)
            {
                int64_t operand = h[at + 2 + 2 * i];
                CHECK(operand == 0 || operand == 1);
                from[operand == 1]++;
            }
            at += 2 + 2 * static_cast<size_t>(n);
        }
        // The cube's 6 faces give 7 (the top is cut in two); the tool gives the 3 slot faces.
        CHECK(from[0] == 7 && from[1] == 3);
        done(r);
        for (smlib_brep *b : {slot, ta, tb})
        {
            smlib_brep_free(b);
        }
    }

    // SMLib defects 2-5 (curves of a failed imprint freed under their Edge, a failed imprint
    // carried on, a squeezed vertex reused, a combined-away ProtoEdge left in its ProtoFaces):
    // inward walls as thick as the cube's inradius or thicker, refused before SMLib sees them
    // (no cavity can remain).
    const double lid[3] = {10, 10, 20};
    CHECK(smlib_shell(cube, lid, 1, 0.1, -10.0, &r, nullptr) == SMLIB_ERR_NO_ROOM && r == nullptr);
    CHECK(smlib_shell(cube, lid, 1, 0.1, -15.0, &r, nullptr) == SMLIB_ERR_NO_ROOM);
    CHECK(smlib_shell(cube, nullptr, 0, 0.1, -10.0, &r, nullptr) == SMLIB_ERR_NO_ROOM);
    if (smlib_shell(cube, lid, 1, 0.1, -9.5, &r, nullptr) == SMLIB_OK)
    {
        CHECK(near(volume_of(r), 8000.0 - 1.0 * 1.0 * 10.5, 1e-6));
        done(r);
    }
    // tangency-labs/solidcore#380: SMLib itself, past the shim's guard. A closed inward shell
    // whose offset leaves nothing handed the cube back unhollowed as a success; it fails now,
    // open or closed, and leaves the input as it was.
    {
        SmTArray<SmFace *> faces;
        cube->p->GetFaces(faces);
        SmTArray<const SmFace *> lidface;
        for (ULONG i = 0; i < faces.GetSize(); ++i)
        {
            SmPoint3d lo;
            SmPoint3d hi;
            if (SmApiFaceBoundingBox(faces[i], TRUE, lo, hi) == SM_SUCCESS && lo[2] > 19.9)
            {
                lidface.Add(faces[i]);
            }
        }
        CHECK(lidface.GetSize() == 1);
        SmTArray<const SmFace *> none;
        for (double wall : {10.0, 12.0, 15.0})
        {
            for (const SmTArray<const SmFace *> *open : {&none, &lidface})
            {
                SmBrep *res = nullptr;
                CHECK(SmApiShellBrepFull(cube->p, -wall, TRUE, TRUE, TRUE, *open, res) != SM_SUCCESS && res == nullptr);
                delete res;
            }
        }
        CHECK(near(volume_of(cube), 8000.0, 1e-9) && sound(cube));
    }

    // SMLib defect 6 (SmOffsetExecutive::RemoveExcessFaces deleted faces its face maps still
    // held, which PiecewiseMerge then read), found by fuzz_ops: a closed inward shell of the
    // cube with two of its edges chamfered 0.1 mm. Any status, never a crash.
    {
        const double two[6] = {0, 0, 10, 0, 20, 10};
        smlib_brep *nicked = nullptr;
        CHECK(smlib_chamfer(cube, two, 2, 0.01, 0.1, &nicked, nullptr) == SMLIB_OK && sound(nicked));
        if (smlib_shell(nicked, nullptr, 0, 0.01, -2.624267578125, &r, nullptr) == SMLIB_OK)
        {
            done(r);
        }
        smlib_brep_free(nicked);
    }

    // tangency-labs/solidcore#380: the guard against walls with no room inside is the body's
    // own inradius, not its box. An L of two 6 mm legs, 40 long and 20 tall: its box allows 10,
    // its inradius is 6 * sqrt(2) / (1 + sqrt(2)) = 3.51 (at the inner corner).
    {
        smlib_brep *legs[2] = {box(0, 0, 0, 40, 6, 20), box(0, 0, 0, 6, 40, 20)};
        smlib_brep *ell = nullptr;
        CHECK(smlib_boolean(legs[0], legs[1], SMLIB_BOOLEAN_UNION, &ell) == SMLIB_OK && sound(ell));
        for (double wall : {5.0, 4.0, 3.6})
        {
            CHECK(smlib_shell(ell, nullptr, 0, 0.1, -wall, &r, nullptr) == SMLIB_ERR_NO_ROOM && r == nullptr);
        }
        // Open at the top, the wall only needs room from the sides and the floor: the same.
        const double top[3] = {20, 3, 20};
        CHECK(smlib_shell(ell, top, 1, 0.1, -3.6, &r, nullptr) == SMLIB_ERR_NO_ROOM);
        if (smlib_shell(ell, top, 1, 0.1, -1.0, &r, nullptr) == SMLIB_OK)
        {
            // Floor 1 thick under a cavity 19 tall, its plan the L inset by 1: 2 * 38 * 4 - 4 * 4.
            CHECK(std::fabs(volume_of(r) - (volume_of(ell) - 19.0 * (2.0 * 38.0 * 4.0 - 16.0))) < 1e-4);
            done(r);
        }
        else
        {
            CHECK(!"the L opened at the top, 1 mm walls");
        }
        smlib_brep_free(ell);
        smlib_brep_free(legs[0]);
        smlib_brep_free(legs[1]);
    }

    // tangency-labs/solidcore#403: a hemisphere r10 (a sphere cut by a box, its rim two arcs of
    // 270 and 90 degrees) shelled 2 open at the flat face. The rim's offset is a closed circle,
    // which both arcs were trimmed to as the same 270 degrees: the flat face got the long arc
    // twice and the cavity was no solid. Each arc now takes the part of the circle it runs along.
    {
        const double centre[3] = {0, 0, 0};
        smlib_brep *ball = nullptr;
        CHECK(smlib_make_sphere(centre, 10.0, &ball) == SMLIB_OK);
        smlib_brep *upper = box(-20, -20, 0, 40, 40, 20);
        smlib_brep *half = nullptr;
        CHECK(smlib_boolean(ball, upper, SMLIB_BOOLEAN_INTERSECT, &half) == SMLIB_OK && sound(half));
        const double flat[3] = {3, 2, 0};
        CHECK(smlib_shell(half, flat, 1, 0.1, -2.0, &r, nullptr) == SMLIB_OK);
        if (r != nullptr)
        {
            CHECK(near(volume_of(r), 2.0 / 3.0 * PI * (1000.0 - 512.0), 1e-5));
            done(r);
        }
        for (smlib_brep *b : {ball, upper, half})
        {
            smlib_brep_free(b);
        }
    }

    // tangency-labs/solidcore#403: the same shell gives the same answer every time, wherever the
    // heap puts SMLib's objects. SMLib's pointer maps (SmMapPtrToPtr, SmRelation) handed out
    // their keys in hash order, which is address order, so the offset processed edges and faces
    // in an order that changed from call to call: a 47.5 x 45 x 30 box joined with a cylinder
    // r7.5 along x whose seam lies on the box's front face, shelled 1 open at the bottom,
    // failed with 1001 or 1006 by turns. Heap blocks of varying size between the calls move
    // the objects about.
    {
        smlib_brep *block = box(0, 0, 0, 47.5, 45, 30);
        SmBrep *bar = nullptr;
        CHECK(SmApiCreateCylinder(SmVector3d(0, 0, 0), 7.5, 30.0, bar) == SM_SUCCESS && bar != nullptr);
        CHECK(SmApiRotate(bar, SmVector3d(0, 0, 0), SmVector3d(0, 1, 0), 90.0) == SM_SUCCESS);
        CHECK(SmApiTranslate(bar, SmVector3d(47.5, 0, 10)) == SM_SUCCESS);
        smlib_brep side{bar};
        smlib_brep *joined = nullptr;
        CHECK(smlib_boolean(block, &side, SMLIB_BOOLEAN_UNION, &joined) == SMLIB_OK && sound(joined));
        const double bottom[3] = {23.75, 22.5, 0};
        std::vector<std::vector<char>> spacers;
        int32_t first = 0;
        double first_volume = 0.0;
        for (int i = 0; i < 6; ++i)
        {
            spacers.emplace_back(static_cast<size_t>(16 + 97 * i * i));
            int32_t got = smlib_shell(joined, bottom, 1, 0.1, -1.0, &r, nullptr);
            double v = got == SMLIB_OK ? volume_of(r) : 0.0;
            if (i == 0)
            {
                first = got;
                first_volume = v;
            }
            CHECK(got == first && v == first_volume);
            smlib_brep_free(r);
            r = nullptr;
        }
        delete bar;
        for (smlib_brep *b : {block, joined})
        {
            smlib_brep_free(b);
        }
    }

    // Rig case shell/1494, the same box and cylinder, its seam turned off the front face: the
    // shell failed (1001). Where the box's side x = 47.5 and the base of the cylinder (the half
    // disc in front of the box, facing -x) meet the front face at the cylinder's rim (z = 2.5
    // and 17.5), the front face's loop runs from its edge with the side to its edge with the
    // base, whose offsets (x = 46.5 and 48.5) are parallel: they never meet, and the vertex kept
    // by the one lay 2 from the other's curve, as in the cylinder's loop. The merge grew that
    // vertex's tolerance to 4, it captured another vertex, and SmMerge::MakeFaceuse met an edge
    // whose ends were not its vertices' mates. The gap is now closed by the straight edge between
    // the two ends (on the front face's offset and on the cylinder's, a ruling). The cavity is
    // the box's inset, the cylinder's inset beyond the base's offset (x 48.5 to 76.5), and
    // between the two offsets of the side the part of the cylinder's inset disc behind the front
    // face's offset (y >= 1): the circular segment of r 6.5 beyond the chord 1 from its centre,
    // 2 long. (With the seam through one of those two vertices, as above, the cylinder's loop
    // passes the vertex twice and its rim's offsets stop at the seam, not at the front face's
    // offset: still refused.)
    {
        const double segment = 6.5 * 6.5 * std::acos(1.0 / 6.5) - std::sqrt(6.5 * 6.5 - 1.0);
        const double saddle = 47.5 * 45.0 * 30.0 + PI * 7.5 * 7.5 * 30.0 - (45.5 * 43.0 * 29.0 + 2.0 * segment + PI * 6.5 * 6.5 * 28.0);
        smlib_brep *block = box(0, 0, 0, 47.5, 45, 30);
        for (double seam : {30.0, 55.0, 125.0, 200.0})
        {
            SmBrep *bar = nullptr;
            CHECK(SmApiCreateCylinder(SmVector3d(0, 0, 0), 7.5, 30.0, bar) == SM_SUCCESS && bar != nullptr);
            CHECK(SmApiRotate(bar, SmVector3d(0, 0, 0), SmVector3d(0, 0, 1), seam) == SM_SUCCESS);
            CHECK(SmApiRotate(bar, SmVector3d(0, 0, 0), SmVector3d(0, 1, 0), 90.0) == SM_SUCCESS);
            CHECK(SmApiTranslate(bar, SmVector3d(47.5, 0, 10)) == SM_SUCCESS);
            smlib_brep side{bar};
            smlib_brep *joined = nullptr;
            CHECK(smlib_boolean(block, &side, SMLIB_BOOLEAN_UNION, &joined) == SMLIB_OK);
            const double bottom[3] = {23.75, 22.5, 0};
            CHECK(joined != nullptr && smlib_shell(joined, bottom, 1, 0.1, -1.0, &r, nullptr) == SMLIB_OK);
            if (r != nullptr)
            {
                CHECK(near(volume_of(r), saddle, 1e-7));
                done(r);
            }
            delete bar;
            smlib_brep_free(joined);
        }
        smlib_brep_free(block);
    }

    // tangency-labs/solidcore#397: an inward shell of a box 30 x 20 x 10 whose edges are rounded
    // thinner than the 2 mm wall. Each round's inward offset collapses (SMLib's offset surface is
    // NULL), which failed the shell (1001); the faces beside it are now mated across it, so the
    // inner walls meet at sharp edges: the cavity is the box shrunk by the wall on every side
    // (but the open top), whatever the radius. All twelve edges rounded (the corner balls
    // collapse too) and only the four vertical ones (their ends meet the top and bottom faces,
    // where the offset edges along them shrink to a corner). A round thicker than the wall (r 3)
    // still shells as before, to a cavity rounded r - t.
    {
        const double a = 30, b = 20, c = 10, t = 2;
        // a x b x c with all edges rounded (corners balls), as the core box grown by the radius
        auto rounded = [](double x, double y, double z, double rad)
        {
            double cx = x - 2 * rad, cy = y - 2 * rad, cz = z - 2 * rad;
            return cx * cy * cz + 2 * rad * (cx * cy + cy * cz + cx * cz) + PI * rad * rad * (cx + cy + cz) + 4.0 / 3.0 * PI * rad * rad * rad;
        };
        std::vector<double> all, vertical;
        for (double y : {0.0, b})
        {
            for (double z : {0.0, c})
            {
                all.insert(all.end(), {a / 2, y, z});
            }
        }
        for (double x : {0.0, a})
        {
            for (double z : {0.0, c})
            {
                all.insert(all.end(), {x, b / 2, z});
            }
        }
        for (double x : {0.0, a})
        {
            for (double y : {0.0, b})
            {
                all.insert(all.end(), {x, y, c / 2});
                vertical.insert(vertical.end(), {x, y, c / 2});
            }
        }
        const double top[3] = {a / 2, b / 2, c};
        smlib_brep *block = box(0, 0, 0, a, b, c);
        struct Case
        {
            const std::vector<double> *edges;
            double radius;
        };
        for (Case k : {Case{&all, 1.0}, Case{&all, 0.5}, Case{&vertical, 1.0}, Case{&vertical, 0.5}, Case{&all, 3.0}})
        {
            const double rad = k.radius;
            const bool every = k.edges == &all;
            smlib_brep *round = nullptr;
            CHECK(smlib_fillet(block, k.edges->data(), static_cast<int64_t>(k.edges->size() / 3), 0.1, rad, &round, nullptr) == SMLIB_OK);
            if (round == nullptr)
            {
                continue;
            }
            const double outer = every ? rounded(a, b, c, rad) : a * b * c - 4.0 * (1.0 - PI / 4.0) * rad * rad * c;
            CHECK(near(volume_of(round), outer, 1e-7));
            const double inner = rad > t ? rounded(a - 2 * t, b - 2 * t, c - 2 * t, rad - t) : (a - 2 * t) * (b - 2 * t) * (c - 2 * t);
            std::printf("box %s rounded r%g, inward shell %g: closed", every ? "all edges" : "vertical edges", rad, t);
            int32_t got = smlib_shell(round, nullptr, 0, 0.1, -t, &r, nullptr);
            std::printf(" %d", static_cast<int>(got));
            CHECK(got == SMLIB_OK);
            if (r != nullptr)
            {
                // (a cavity rounded r - t is SMLib's approximation of the offset, to its tolerance)
                CHECK(near(volume_of(r), outer - inner, rad > t ? 1e-6 : 1e-7));
                done(r);
            }
            if (rad <= t)
            {
                got = smlib_shell(round, top, 1, 0.1, -t, &r, nullptr);
                std::printf(", open at the top %d", static_cast<int>(got));
                CHECK(got == SMLIB_OK);
                if (r != nullptr)
                {
                    CHECK(near(volume_of(r), outer - (a - 2 * t) * (b - 2 * t) * (c - t), 1e-7));
                    done(r);
                }
            }
            std::printf("\n");
            smlib_brep_free(round);
        }
        smlib_brep_free(block);
    }

    // An offset is no surface exactly where it shrinks its surface onto its spine (a cylinder or
    // cone onto its axis, a torus onto its centre circle, a sphere onto its centre): radius 2,
    // offset 2.1 inward collapses, 1.9 does not, and outward never; the cone (radii 2 and 1 over
    // 5) collapses only when both its radii do, at 2.1 (cos of its semi-angle 0.98), not at 1.5,
    // where only its top passes the axis.
    {
        SmContext *ctx = SmApiGetOrCreateContext();
        SmSurface *cyl = nullptr, *cone = nullptr, *tor = nullptr;
        CHECK(SmApiCreateCylinder(SmVector3d(0, 0, 0), 2.0, 5.0, cyl) == SM_SUCCESS && cyl != nullptr);
        CHECK(SmApiCreateCone(SmVector3d(0, 0, 0), 2.0, 1.0, 5.0, cone) == SM_SUCCESS && cone != nullptr);
        CHECK(SmApiCreateTorus(SmVector3d(0, 0, 0), 10.0, 2.0, tor) == SM_SUCCESS && tor != nullptr);
        SmBrep *ball = nullptr;
        CHECK(SmApiCreateSphere(SmVector3d(0, 0, 0), 2.0, ball) == SM_SUCCESS && ball != nullptr);
        SmTArray<SmFace *> ball_faces;
        ball->GetFaces(ball_faces);
        CHECK(ball_faces.GetSize() > 0);
        const SmSurface *sph = ball_faces.GetSize() > 0 ? ball_faces[0]->GetSurface() : nullptr;
        const double tol = 1e-4;
        for (const SmSurface *s : {static_cast<const SmSurface *>(cyl), static_cast<const SmSurface *>(cone), static_cast<const SmSurface *>(tor), sph})
        {
            if (s == nullptr)
            {
                continue;
            }
            int collapsed = 0;
            for (double d : {-2.1, 2.1, -1.9, 1.9, -1.5, 1.5})
            {
                SmSurface *off = nullptr;
                CHECK(s->CreateOffsetSurface(*ctx, d, tol, off) == SM_SUCCESS);
                const bool collapses = s->OffsetCollapses(d, tol) != FALSE;
                CHECK((off == nullptr) == collapses);
                if (collapses)
                {
                    collapsed++;
                    CHECK(std::fabs(d) > 2.0);
                }
                delete off;
            }
            CHECK(collapsed == 1);
        }
        delete cyl;
        delete cone;
        delete tor;
        delete ball;
    }

    // The offset of a face: no surface only for a collapse, flagged as one; an offset construction
    // that gives none otherwise fails. A cylinder r5: inset 6 it collapses (one of its two
    // faceuses offsets toward the axis), inset 2 it is a cylinder r3; with a construction that
    // gives nothing, inset 2 fails on both sides and is no collapse.
    {
        smlib_brep *rod = cylinder(0, 0, 0, 5, 10);
        SmTArray<SmFace *> faces;
        rod->p->GetFaces(faces);
        SmFace *side = nullptr;
        for (ULONG i = 0; i < faces.GetSize(); i++)
        {
            if (faces[i]->GetSurface()->IsKindOf(SmCylinder_TYPE))
            {
                side = faces[i];
            }
        }
        CHECK(side != nullptr);
        if (side != nullptr)
        {
            const SmContext &ctx = *SmApiGetOrCreateContext();
            const SmFaceuse *uses[2] = {side->GetUpwardFaceuse(), side->GetUpwardFaceuse()->GetMate()};
            auto offset = [&](double d, const SmFaceuse *fu, SmBoolean &collapses, SmSurface *&surface)
            {
                SmOffsetGeometryCreation make(d, rod->p->GetTolerance());
                SmSSIData crossings;
                surface = nullptr;
                collapses = FALSE;
                return make.FaceuseOffset(ctx, fu, TRUE, surface, crossings, collapses);
            };
            int collapsed = 0;
            for (const SmFaceuse *fu : uses)
            {
                SmBoolean collapses = FALSE;
                SmSurface *surface = nullptr;
                CHECK(offset(6.0, fu, collapses, surface) == SM_SUCCESS);
                CHECK((surface == nullptr) == (collapses != FALSE));
                collapsed += collapses ? 1 : 0;
                delete surface;
                CHECK(offset(2.0, fu, collapses, surface) == SM_SUCCESS && !collapses && surface != nullptr);
                delete surface;
            }
            CHECK(collapsed == 1);
            SmSurface *own = side->GetSurface();
            const SmExtent2d domain = side->GetUVDomain();
            NoOffsetCylinder *none = new (ctx) NoOffsetCylinder(*static_cast<const SmCylinder *>(own));
            side->SetSurface(none);
            for (const SmFaceuse *fu : uses)
            {
                SmBoolean collapses = FALSE;
                SmSurface *surface = nullptr;
                CHECK(offset(2.0, fu, collapses, surface) != SM_SUCCESS && !collapses && surface == nullptr);
            }
            side->SetSurface(own);
            side->SetUVDomain(domain);
            delete none;
            CHECK(sound(rod));
        }
        smlib_brep_free(rod);
    }

    // A box 20 x 10 x 30 open at the top with walls 4 thick: the sides need room 4 from each
    // other (10 > 8), the floor's is far, so the shell goes through even though the box's
    // smallest half extent is 5; walls of 5 leave none.
    {
        smlib_brep *slab = box(0, 0, 0, 20, 10, 30);
        const double top[3] = {10, 5, 30};
        CHECK(smlib_shell(slab, top, 1, 0.1, -4.0, &r, nullptr) == SMLIB_OK);
        if (r != nullptr)
        {
            CHECK(std::fabs(volume_of(r) - (6000.0 - 12.0 * 2.0 * 26.0)) < 1e-6);
            done(r);
        }
        CHECK(smlib_shell(slab, top, 1, 0.1, -5.0, &r, nullptr) == SMLIB_ERR_NO_ROOM);
        smlib_brep_free(slab);
    }

    // Offset faces: the cube's top up by 5, then down by 5; its side and top together; a
    // cylinder's wall out by 2. Down by the cube's height leaves nothing and is refused.
    {
        const double top[3] = {10, 10, 20};
        CHECK(smlib_offset_faces(cube, top, 1, 0.1, 5.0, &r, nullptr) == SMLIB_OK);
        if (r != nullptr)
        {
            CHECK(std::fabs(volume_of(r) - 10000.0) < 1e-6);
            done(r);
        }
        CHECK(smlib_offset_faces(cube, top, 1, 0.1, -5.0, &r, nullptr) == SMLIB_OK);
        if (r != nullptr)
        {
            CHECK(std::fabs(volume_of(r) - 6000.0) < 1e-6);
            done(r);
        }
        const double two[6] = {10, 10, 20, 20, 10, 10};
        CHECK(smlib_offset_faces(cube, two, 2, 0.1, 2.0, &r, nullptr) == SMLIB_OK);
        if (r != nullptr)
        {
            CHECK(std::fabs(volume_of(r) - 22.0 * 20.0 * 22.0) < 1e-6);
            done(r);
        }
        // Up by three times the cube's height: the sides are extended far enough to meet
        // the moved top (they were sized from the cube alone and stopped short of it).
        CHECK(smlib_offset_faces(cube, top, 1, 0.1, 60.0, &r, nullptr) == SMLIB_OK);
        if (r != nullptr)
        {
            CHECK(std::fabs(volume_of(r) - 20.0 * 20.0 * 80.0) < 1e-6);
            done(r);
        }
        CHECK(smlib_offset_faces(cube, top, 1, 0.1, -20.0, &r, nullptr) == SMLIB_ERR_NO_ROOM);
        CHECK(smlib_offset_faces(cube, top, 1, 0.1, 0.0, &r, nullptr) == SMLIB_ERR_ARG);
        const double wall[3] = {14, 10, 5};
        if (smlib_offset_faces(rod, wall, 1, 0.1, 2.0, &r, nullptr) == SMLIB_OK)
        {
            done(r);
        }
        // tangency-labs/solidcore#401: a sphere r10 out by 2 is a sphere r12 whose edges lie
        // within their tolerances (the gap was 5.015e-6 over 5e-6 before the results' tolerances
        // were settled).
        const double centre[3] = {0, 0, 0};
        smlib_brep *ball = nullptr;
        CHECK(smlib_make_sphere(centre, 10.0, &ball) == SMLIB_OK);
        const double pole[3] = {10, 0, 0};
        CHECK(smlib_offset_faces(ball, pole, 1, 0.1, 2.0, &r, nullptr) == SMLIB_OK);
        if (r != nullptr)
        {
            CHECK(near(volume_of(r), 4.0 / 3.0 * PI * 1728.0, 1e-6));
            done(r);
        }
        smlib_brep_free(ball);
    }

    // tangency-labs/solidcore#402: a face moved away from a tangent neighbour that keeps its
    // surface no longer meets it. SMLib's local operation dropped the edge across the gap and
    // returned a body whose edge tolerance had grown to the gap. A face moved along its normal
    // stays parallel to a tangent neighbour that is no blend: a step, closed by a wall between
    // the two (one piece of a split top moved along its normal: a 60 x 40 x 20 block, the top
    // split at x = 30, the piece x > 30 up by 1 adds 30 x 40 x 1, down by 1 takes it off).
    //
    // A blend tangent to a moved face follows it (SmOffsetFacesExecutive): taken off, the faces
    // moved, put back at its own radius, tangent to both faces. A boss of radius R, 10 high,
    // on a 60 x 40 x 10 plate: a round r at its foot adds, and one at its rim takes off, the
    // corner region of the section, A = r^2 (1 - pi/4), revolved about the axis at R + d or
    // R - d, d = r^3 (5/6 - pi/4) / A its centroid's distance from the wall (Pappus).
    // tangency-labs/solidcraft#3468: the r10 boss's wall out by 1 with its r2 foot round,
    // 27863.06545 and 9 faces (the round ends at r11, not at r10 under a step); in by 1 and 3,
    // the same with R = 9 and 7. A rim round r2, the wall in to R = 3: the ball rolls inside a
    // boss of radius 3, the top keeps a disc of radius 1; to R = 2: the ball's centre is on the
    // axis, the top goes and the round is the hemispherical end (8 faces), the same formula's
    // limit; to R = 1.5: no ball of r2 fits inside a boss of radius 1.5 (its centre would lie
    // beyond the axis, the round a spindle torus through itself), SMLIB_ERR_NO_ROOM. Both rounds
    // together follow the wall too. The boss's top moved instead takes its rim round along: as
    // a boss built that high and rounded (the local operation continues the wall's cylinder
    // past its domain by a G1 extension, within 3e-7 of the cylinder here, and so is not exact:
    // 2e-6 relative).
    {
        smlib_brep *plate = box(0, 0, 0, 60, 40, 10);
        const double at[3] = {30, 20, 10};
        smlib_brep *boss = nullptr;
        smlib_brep *joined = nullptr;
        smlib_brep *foot = nullptr;
        smlib_brep *rim = nullptr;
        smlib_brep *both = nullptr;
        CHECK(smlib_make_cylinder(at, 10.0, 10.0, &boss) == SMLIB_OK);
        CHECK(smlib_boolean(plate, boss, SMLIB_BOOLEAN_UNION, &joined) == SMLIB_OK);
        const double base[3] = {30, 30, 10};
        const double top_rim[3] = {30, 30, 20};
        const double both_rims[6] = {30, 30, 10, 30, 30, 20};
        CHECK(smlib_fillet(joined, base, 1, 0.1, 2.0, &foot, nullptr) == SMLIB_OK && sound(foot));
        CHECK(smlib_fillet(joined, top_rim, 1, 0.1, 2.0, &rim, nullptr) == SMLIB_OK && sound(rim));
        CHECK(smlib_fillet(joined, both_rims, 2, 0.1, 2.0, &both, nullptr) == SMLIB_OK && sound(both));
        const double corner = 4.0 * (1.0 - PI / 4.0);
        const double centroid = 8.0 * (5.0 / 6.0 - PI / 4.0) / corner;
        auto boss_of = [&](double R, bool at_foot, bool at_rim)
        {
            return 24000.0 + PI * R * R * 10.0 + (at_foot ? 2.0 * PI * (R + centroid) * corner : 0.0) -
                   (at_rim ? 2.0 * PI * (R - centroid) * corner : 0.0);
        };
        CHECK(near(boss_of(11, true, false), 27863.06545197909, 1e-12));
        const double wall[3] = {30, 30, 15};
        for (double R : {11.0, 9.0, 7.0})
        {
            CHECK(smlib_offset_faces(foot, wall, 1, 0.1, R - 10.0, &r, nullptr) == SMLIB_OK);
            if (r != nullptr)
            {
                CHECK(near(volume_of(r), boss_of(R, true, false), 1e-9));
                CHECK(r->p->GetNumFaces() == 9);
                done(r);
            }
        }
        for (double R : {3.0, 2.0})
        {
            CHECK(smlib_offset_faces(rim, wall, 1, 0.1, R - 10.0, &r, nullptr) == SMLIB_OK);
            if (r != nullptr)
            {
                CHECK(near(volume_of(r), boss_of(R, false, true), 1e-9));
                CHECK(r->p->GetNumFaces() == (R == 2.0 ? 8u : 9u));
                done(r);
            }
        }
        CHECK(smlib_offset_faces(rim, wall, 1, 0.1, 1.5 - 10.0, &r, nullptr) == SMLIB_ERR_NO_ROOM && r == nullptr);
        for (double R : {11.0, 7.0})
        {
            CHECK(smlib_offset_faces(both, wall, 1, 0.1, R - 10.0, &r, nullptr) == SMLIB_OK);
            if (r != nullptr)
            {
                CHECK(near(volume_of(r), boss_of(R, true, true), 1e-9));
                CHECK(r->p->GetNumFaces() == 10);
                done(r);
            }
        }
        const double cap[3] = {30, 20, 20};
        for (double up : {1.0, -1.0})
        {
            CHECK(smlib_offset_faces(rim, cap, 1, 0.1, up, &r, nullptr) == SMLIB_OK);
            if (r != nullptr)
            {
                CHECK(near(volume_of(r), boss_of(10, false, true) + PI * 100.0 * up, 2e-6));
                CHECK(r->p->GetNumFaces() == 9);
                done(r);
            }
        }
        smlib_brep *block = box(0, 0, 0, 60, 40, 20);
        smlib_brep *cutter = box(30, -1, 10, 40, 42, 20);
        smlib_brep *split = nullptr;
        CHECK(smlib_imprint(block, cutter, &split) == SMLIB_OK);
        const double piece[3] = {45, 20, 20};
        CHECK(smlib_offset_faces(split, piece, 1, 0.1, 1.0, &r, nullptr) == SMLIB_OK && sound(r));
        CHECK(near(volume_of(r), 48000.0 + 1200.0, 1e-6));
        smlib_brep_free(r);
        r = nullptr;
        CHECK(smlib_offset_faces(split, piece, 1, 0.1, -1.0, &r, nullptr) == SMLIB_OK && sound(r));
        CHECK(near(volume_of(r), 48000.0 - 1200.0, 1e-6));
        smlib_brep_free(r);
        r = nullptr;
        for (smlib_brep *b : {plate, boss, joined, foot, rim, both, block, cutter, split})
        {
            smlib_brep_free(b);
        }
    }

    // Blends with ends follow too. A 60 x 40 x 20 block, its top front edge rounded r2: the top
    // up or down 1 adds or takes 60 x 40 x 1, the round (its section's corner region, A, along
    // the 60) as it was. A 40 x 30 x 20 block, its four upright edges rounded r3, the side
    // x = 40 moved in: the rounds on it follow, a W x 30 rounded rectangle 20 high,
    // (30 W - 4 A) 20 with A = 9 (1 - pi/4). At W = 10 the front and back keep 4 of their flat;
    // at W = 6 they keep none, the rounds meet tangentially (8 faces); at W = 5 the rounds of
    // r3 no longer fit beside each other (a ball of r3 in a slab 5 thick): SMLIB_ERR_NO_ROOM.
    // At W = 1.5 the side would cut into the rounds that stay at x = 0: the local operation
    // refuses it (as it does for a block with no rounds on the moved side).
    {
        smlib_brep *block = box(0, 0, 0, 60, 40, 20);
        smlib_brep *rounded = nullptr;
        const double front_top[3] = {30, 0, 20};
        CHECK(smlib_fillet(block, front_top, 1, 0.1, 2.0, &rounded, nullptr) == SMLIB_OK);
        const double top[3] = {30, 20, 20};
        for (double up : {1.0, -1.0})
        {
            CHECK(smlib_offset_faces(rounded, top, 1, 0.1, up, &r, nullptr) == SMLIB_OK);
            if (r != nullptr)
            {
                CHECK(near(volume_of(r), 2400.0 * (20.0 + up) - 60.0 * 4.0 * (1.0 - PI / 4.0), 1e-9));
                CHECK(r->p->GetNumFaces() == 7);
                done(r);
            }
        }
        smlib_brep *slab = box(0, 0, 0, 40, 30, 20);
        smlib_brep *corners = nullptr;
        const double uprights[12] = {0, 0, 10, 40, 0, 10, 40, 30, 10, 0, 30, 10};
        CHECK(smlib_fillet(slab, uprights, 4, 0.1, 3.0, &corners, nullptr) == SMLIB_OK);
        const double side[3] = {40, 15, 10};
        for (double W : {10.0, 6.0})
        {
            CHECK(smlib_offset_faces(corners, side, 1, 0.1, W - 40.0, &r, nullptr) == SMLIB_OK);
            if (r != nullptr)
            {
                CHECK(near(volume_of(r), (30.0 * W - 36.0 * (1.0 - PI / 4.0)) * 20.0, 1e-9));
                CHECK(r->p->GetNumFaces() == (W == 6.0 ? 8u : 10u));
                done(r);
            }
        }
        CHECK(smlib_offset_faces(corners, side, 1, 0.1, 5.0 - 40.0, &r, nullptr) == SMLIB_ERR_NO_ROOM && r == nullptr);
        CHECK(smlib_offset_faces(corners, side, 1, 0.1, 1.5 - 40.0, &r, nullptr) != SMLIB_OK && r == nullptr);
        // tangency-labs/solidcraft#4492: a 40 x 30 x 10 block, its four top edges rounded r3,
        // mitred where they meet, the top moved: the rounds, which meet at their ends, come off
        // together (the planes round them extended) and go back on the moved top. Each round
        // takes its section's corner region A = 9 (1 - pi/4) along its edge; at a corner the
        // two overlap in the corner cube: of the 2 A r both count there, the mitred corner
        // takes r^3 / 3 (the half of the cube nearer each face, the region outside the round's
        // quarter disc weighted by its distance from the other face: r^3 / 6 each). Up 19 (29
        // high) and down 6 (4 high); the same as the block built that high and rounded.
        smlib_brep *top_rounded = nullptr;
        smlib_brep *low = box(0, 0, 0, 40, 30, 10);
        const double top_edges[12] = {20, 0, 10, 20, 30, 10, 0, 15, 10, 40, 15, 10};
        CHECK(smlib_fillet(low, top_edges, 4, 0.1, 3.0, &top_rounded, nullptr) == SMLIB_OK);
        const double low_top[3] = {20, 15, 10};
        const double A3 = 9.0 * (1.0 - PI / 4.0);
        for (double up : {19.0, -6.0})
        {
            CHECK(smlib_offset_faces(top_rounded, low_top, 1, 0.1, up, &r, nullptr) == SMLIB_OK);
            if (r != nullptr)
            {
                // to the volume's own accuracy (1e-7: the mitre faces' integration)
                CHECK(near(volume_of(r), 1200.0 * (10.0 + up) - A3 * 140.0 + 4.0 * (2.0 * A3 * 3.0 - 9.0), 1e-7));
                CHECK(r->p->GetNumFaces() == 10);
                done(r);
            }
        }
        for (smlib_brep *b : {low, top_rounded})
        {
            smlib_brep_free(b);
        }
        for (smlib_brep *b : {block, rounded, slab, corners})
        {
            smlib_brep_free(b);
        }
    }

    // The blend that follows keeps the history of the blend it replaces: the r2 foot round of
    // the boss above, tagged, its wall out by 1: one result face comes from the round, kept
    // (its boundary the round's own edges), and its rails from the round's.
    {
        smlib_brep *plate = box(0, 0, 0, 60, 40, 10);
        const double at[3] = {30, 20, 10};
        smlib_brep *boss = nullptr;
        smlib_brep *joined = nullptr;
        smlib_brep *foot = nullptr;
        smlib_brep *tagged = nullptr;
        CHECK(smlib_make_cylinder(at, 10.0, 10.0, &boss) == SMLIB_OK);
        CHECK(smlib_boolean(plate, boss, SMLIB_BOOLEAN_UNION, &joined) == SMLIB_OK);
        const double base[3] = {30, 30, 10};
        CHECK(smlib_fillet(joined, base, 1, 0.1, 2.0, &foot, nullptr) == SMLIB_OK);
        CHECK(smlib_history_tag(foot, 0, &tagged) == SMLIB_OK);
        const double wall[3] = {30, 30, 15};
        CHECK(smlib_offset_faces(tagged, wall, 1, 0.1, 1.0, &r, nullptr) == SMLIB_OK);
        if (r != nullptr)
        {
            std::vector<int64_t> h(100000);
            int64_t len = 0;
            CHECK(smlib_history_read(r, h.data(), static_cast<int64_t>(h.size()), &len) == SMLIB_OK);
            // [version, faces, edges, per face: change, n, (operand, index) x n, ...]
            int64_t nf = h[1];
            int64_t generated = 0;
            size_t at_rec = 3;
            for (int64_t f = 0; f < nf; ++f)
            {
                int64_t change = h[at_rec];
                int64_t n = h[at_rec + 1];
                if (n == 0 || change == SMLIB_CHANGE_GENERATED)
                {
                    generated++;
                }
                at_rec += 2 + 2 * static_cast<size_t>(n);
            }
            // every result face comes from an input face: none is new
            CHECK(generated == 0);
            CHECK(r->p->GetNumFaces() == 9);
            done(r);
        }
        for (smlib_brep *b : {plate, boss, joined, foot, tagged})
        {
            smlib_brep_free(b);
        }
    }

    // Draft: the cube's side x = 20 tilted 5 degrees about its line on z = 0, leaning in going
    // up (+z): a wedge of 20 x 20 tan(5 deg) / 2 x 20 comes off; -5 degrees adds it. The top
    // isn't parallel to the pull: not supported. The rod's wall (r 4, z -1 to 21) becomes the
    // frustum from r 4 + tan 5 at its foot to 4 - 21 tan 5 at its top.
    {
        const double side[3] = {20, 10, 10};
        const double origin[3] = {0, 0, 0};
        const double up[3] = {0, 0, 1};
        const double deg = 3.14159265358979323846 / 180.0;
        const double wedge = 0.5 * 20.0 * 20.0 * std::tan(5.0 * deg) * 20.0;
        CHECK(smlib_draft(cube, side, 1, 0.1, origin, up, up, 5.0 * deg, &r, nullptr) == SMLIB_OK);
        if (r != nullptr)
        {
            CHECK(std::fabs(volume_of(r) - (8000.0 - wedge)) < 1e-6);
            done(r);
        }
        CHECK(smlib_draft(cube, side, 1, 0.1, origin, up, up, -5.0 * deg, &r, nullptr) == SMLIB_OK);
        if (r != nullptr)
        {
            CHECK(std::fabs(volume_of(r) - (8000.0 + wedge)) < 1e-6);
            done(r);
        }
        const double top[3] = {10, 10, 20};
        CHECK(smlib_draft(cube, top, 1, 0.1, origin, up, up, 5.0 * deg, &r, nullptr) == SMLIB_ERR_UNSUPPORTED);
        const double wall[3] = {14, 10, 5};
        CHECK(smlib_draft(rod, wall, 1, 0.1, origin, up, up, 5.0 * deg, &r, nullptr) == SMLIB_OK);
        if (r != nullptr)
        {
            const double r0 = 4.0 + std::tan(5.0 * deg), r1 = 4.0 - 21.0 * std::tan(5.0 * deg);
            CHECK(std::fabs(volume_of(r) - 3.14159265358979323846 * 22.0 / 3.0 * (r0 * r0 + r0 * r1 + r1 * r1)) < 1e-6);
            done(r);
        }
        const double flat[3] = {1, 0, 0};
        CHECK(smlib_draft(cube, side, 1, 0.1, origin, flat, up, 5.0 * deg, &r, nullptr) == SMLIB_ERR_ARG);
        CHECK(smlib_draft(cube, side, 1, 0.1, origin, up, up, 0.0, &r, nullptr) == SMLIB_ERR_ARG);
        CHECK(smlib_draft(cube, side, 1, 0.1, origin, up, up, 90.0 * deg, &r, nullptr) == SMLIB_ERR_ARG);
        CHECK(smlib_draft(cube, side, 1, 0.1, nullptr, up, up, 5.0 * deg, &r, nullptr) == SMLIB_ERR_ARG);
    }

    // A cone frustum whose side is two half-cone faces (a profile of two half circles r 20,
    // drafted inward), shelled 2 mm open at the top (rig case shell/3103). The inner rim meets
    // the top plane where the offsets of the side and the top were already mated along it, and
    // the offset's self-intersection merge intersected them again: the plane/offset-cone curve
    // ended 0.0036 mm short of the inner seam, which split each half of the rim there (three
    // extra vertices and edges, then the shell's difference left partially coincident edges).
    // The two halves of each wall and their rims share edges; the cavity is the frustum moved
    // in by the wall normal to the cone: r(z) = 20 - z tan(a) - t / cos(a), floor at z = t.
    for (double deg : {9.4623222080, 15.0})
    {
        const double a = deg * PI / 180.0, t = 2.0, h = 30.0, r0 = 20.0;
        const int64_t kinds[2] = {SMLIB_CURVE_ARC, SMLIB_CURVE_ARC};
        const double reals[24] = {0, 0, 0, 0, 0, 1, 1, 0, 0, r0, 0, PI, 0, 0, 0, 0, 0, 1, 1, 0, 0, r0, PI, 2 * PI};
        const double up[3] = {0, 0, 1};
        smlib_brep *frustum = nullptr;
        CHECK(smlib_extrude(kinds, 2, reals, 24, up, h, a, SMLIB_DRAFT_CORNER_TANGENT, 1, &frustum, nullptr) == SMLIB_OK);
        auto cone = [](double ra, double rb, double len) { return PI * len / 3.0 * (ra * ra + ra * rb + rb * rb); };
        auto inner = [&](double z) { return r0 - z * std::tan(a) - t / std::cos(a); };
        const double want = cone(r0, r0 - h * std::tan(a), h) - cone(inner(t), inner(h), h - t);
        const double lid[3] = {0, 0, h};
        CHECK(frustum != nullptr && smlib_shell(frustum, lid, 1, 0.1, -t, &r, nullptr) == SMLIB_OK);
        if (r != nullptr)
        {
            int64_t counts[3] = {0, 0, 0};
            CHECK(smlib_brep_counts(r, counts) == SMLIB_OK);
            CHECK(counts[0] == 7 && counts[1] == 12 && counts[2] == 8);
            CHECK(near(volume_of(r), want, 1e-5));
            done(r);
        }
        smlib_brep_free(frustum);
    }

    // A 40 x D x 10 plate with a 4 mm through hole at (x, y), shelled t inward open at the top:
    // where the hole lies nearer than t to a side, the offset of its wall (radius 2 + t) crosses
    // the side's offset, and the offset floor's hole loop crossed its outer loop (SmFace test
    // #38). The merge failed on that face (SM_ERR_NULL_POINTER, SmEdgeDefinition::Imprint); the
    // face is now rebuilt from the stretches of its loops that bound it, and comes apart in two
    // when the hole crosses both sides (D = 8). The cavity is the inner box less the part of the
    // hole's offset disk inside the floor, over the height 10 - t. A hole whose offset only
    // touches a side's (x = 5: the hole 3 = 2 t from the side's offset, its offset circle tangent
    // to the side's at (1.5, 20)) leaves the offset walls meeting along a line, which bounds no
    // solid: refused, SMLIB_ERR_WALLS_TOUCH (it was SM_ERR_BAD_INTERSECTIONS, "loops leave a
    // crossing more than one way", which does not say why). At x = 4.99 the hole's
    // offset crosses the side's 0.01 mm deep, at 4.3 degrees: the cut-off rims of that crossing
    // were welded by the offset's stitch, whose tolerance escalated to the gap between them, and
    // the shell failed (1001).
    {
        struct Plate
        {
            double x, y, d, t;
        };
        // Area of a disk of radius r between the lines at signed distances a < b from its centre.
        auto band = [](double rr, double a, double b)
        {
            auto F = [rr](double u)
            {
                u = std::max(-rr, std::min(rr, u));
                return u * std::sqrt(rr * rr - u * u) + rr * rr * std::asin(u / rr);
            };
            return F(b) - F(a);
        };
        for (Plate p : {Plate{4.5, 20, 40, 1.5}, Plate{4.5, 20, 40, 1.0}, Plate{20, 4, 8, 1.5}, Plate{5, 20, 40, 1.5},
                        Plate{4.99, 20, 40, 1.5}})
        {
            const double o[3] = {0, 0, 0}, size[3] = {40, p.d, 10}, axis[3] = {p.x, p.y, -1}, lid[3] = {35, p.d / 2, 10};
            smlib_brep *box = nullptr, *pin = nullptr, *plate = nullptr;
            CHECK(smlib_make_box(o, size, &box) == SMLIB_OK && smlib_make_cylinder(axis, 2.0, 12.0, &pin) == SMLIB_OK);
            CHECK(smlib_boolean(box, pin, SMLIB_BOOLEAN_SUBTRACT, &plate) == SMLIB_OK);
            const double rr = 2.0 + p.t;
            int32_t s = plate != nullptr ? smlib_shell(plate, lid, 1, 0.1, -p.t, &r, nullptr) : SMLIB_ERR_ARG;
            if (p.x - rr == p.t)
            {
                CHECK(s == SMLIB_ERR_WALLS_TOUCH && r == nullptr);
            }
            else
            {
                double inside = band(rr, std::max(p.t - p.y, -rr), std::min(p.d - p.t - p.y, rr));
                if (p.x - rr < p.t)
                {
                    inside -= band(rr, -rr, p.t - p.x);
                }
                const double cavity = ((40.0 - 2.0 * p.t) * (p.d - 2.0 * p.t) - inside) * (10.0 - p.t);
                CHECK(s == SMLIB_OK);
                if (r != nullptr)
                {
                    CHECK(near(volume_of(r), 40.0 * p.d * 10.0 - PI * 4.0 * 10.0 - cavity, 1e-7));
                    done(r);
                }
            }
            for (smlib_brep *b : {box, pin, plate})
            {
                smlib_brep_free(b);
            }
        }
    }

    // Rig case shell/4711: a 60 x 40 x 10 block with an r3 through hole at (22.5, 7) and a 10 x 10
    // notch at x 25..35 on its front edge, shelled t inward open at the bottom. The hole lies
    // 4 = 2 x 2 from the front face, so at t = 2 the hole's offset circle (r5) is tangent to the
    // front's offset y = 2 at (22.5, 2): the walls meet along the line x = 22.5, y = 2 and the
    // offset bounds no solid there, SMLIB_ERR_WALLS_TOUCH (it was SM_ERR_BAD_INTERSECTIONS, as
    // "loops meet at a loop vertex" here, where the circle's seam lies on the notch's offset at
    // (27.5, 7) and its crossing there was looked at first, and as "loops leave a crossing more
    // than one way" on SolidCraft's copy of the body). Thinner, the walls stay apart; thicker,
    // they cross and merge; both are solids.
    {
        const double o[3] = {0, 0, 0}, size[3] = {60, 40, 10}, axis[3] = {22.5, 7, -5}, corner[3] = {25, -5, -5}, notch_size[3] = {10, 10, 20};
        const double bottom[3] = {50, 30, 0};
        smlib_brep *block = nullptr, *pin = nullptr, *notch = nullptr, *holed = nullptr, *part = nullptr;
        CHECK(smlib_make_box(o, size, &block) == SMLIB_OK && smlib_make_cylinder(axis, 3.0, 20.0, &pin) == SMLIB_OK &&
              smlib_make_box(corner, notch_size, &notch) == SMLIB_OK);
        CHECK(smlib_boolean(block, pin, SMLIB_BOOLEAN_SUBTRACT, &holed) == SMLIB_OK &&
              smlib_boolean(holed, notch, SMLIB_BOOLEAN_SUBTRACT, &part) == SMLIB_OK);
        CHECK(part != nullptr && smlib_shell(part, bottom, 1, 0.1, -2.0, &r, nullptr) == SMLIB_ERR_WALLS_TOUCH && r == nullptr);
        for (double t : {1.5, 2.5})
        {
            CHECK(part != nullptr && smlib_shell(part, bottom, 1, 0.1, -t, &r, nullptr) == SMLIB_OK);
            if (r != nullptr)
            {
                done(r);
            }
        }
        for (smlib_brep *b : {block, pin, notch, holed, part})
        {
            smlib_brep_free(b);
        }
    }

    // Rig case shell/4567: a vase revolved about y from a fit-point spline (its six cubic spans
    // here, from the axis at the origin to the rim at (6, 46)), shelled 2 inward open at the top.
    // The spline meets the axis at a slant, so the bottom's inward offset crosses the axis, 1.6
    // above the frame's origin. SmSurfOfRevolution::CreateOffsetSurface trims the offset
    // generatrix at the axis, but looked for it along a segment from the frame's origin one way
    // only, along the frame's z axis: revolved about -y (as the rig's vase is), the crossing lies
    // the other way, it was missed, the offset was revolved 0.12 past the axis, and the bottom's
    // offset face is no cone point but a hole round it. Revolved about +y it was found. The offset surface now ends on the axis, at the exact offset's crossing.
    // The profile and its exact inward offset p + t n (n the unit normal into the material)
    // are evaluated here directly.
    {
        const double t = 2.0, top = 46.0;
        static const double S[6][4][2] = {{{0, 0}, {5.0046051747861462, 0.28969795384163027}, {10.009210349572294, 0.57939590768326088}, {14, 3}},
                                          {{14, 3}, {17.791122091620693, 5.2994962033305715}, {20.667337490527675, 9.5220052633338028}, {22, 14}},
                                          {{22, 14}, {23.175751555936337, 17.950744649667318}, {23.150034916853265, 22.100353826034457}, {22, 26}},
                                          {{22, 26}, {20.772697258333796, 30.161653155655173}, {18.264920922540426, 34.038626225241856}, {14, 36}},
                                          {{14, 36}, {11.021250018636229, 37.369882865619019}, {7.1853563364202362, 37.805327050833874}, {6, 40}},
                                          {{6, 40}, {5.2048387962612725, 41.472231370814619}, {5.6024193981306363, 43.736115685407313}, {6, 46}}};
        // span k of the profile at s: (r, y), its first and second derivatives
        auto bez = [](int k, double s, double d[3][2])
        {
            const double(*c)[2] = S[k];
            const double u = 1 - s;
            for (int i = 0; i < 2; ++i)
            {
                d[0][i] = u * u * u * c[0][i] + 3 * u * u * s * c[1][i] + 3 * u * s * s * c[2][i] + s * s * s * c[3][i];
                d[1][i] = 3 * (u * u * (c[1][i] - c[0][i]) + 2 * u * s * (c[2][i] - c[1][i]) + s * s * (c[3][i] - c[2][i]));
                d[2][i] = 6 * (u * (c[2][i] - 2 * c[1][i] + c[0][i]) + s * (c[3][i] - 2 * c[2][i] + c[1][i]));
            }
        };
        // the exact offset q = p + t (-y', r') / |p'| and its derivative
        auto off = [&](int k, double s, double q[2], double dq[2])
        {
            double d[3][2];
            bez(k, s, d);
            const double len = std::hypot(d[1][0], d[1][1]), dlen = (d[1][0] * d[2][0] + d[1][1] * d[2][1]) / len;
            q[0] = d[0][0] - t * d[1][1] / len;
            q[1] = d[0][1] + t * d[1][0] / len;
            dq[0] = d[1][0] - t * (d[2][1] / len - d[1][1] * dlen / (len * len));
            dq[1] = d[1][1] + t * (d[2][0] / len - d[1][0] * dlen / (len * len));
        };
        // where g changes sign in [a, b] (bisection)
        auto root = [](const std::function<double(double)> &g, double a, double b)
        {
            for (int i = 0; i < 200; ++i)
            {
                const double m = 0.5 * (a + b);
                ((g(a) < 0) == (g(m) < 0) ? a : b) = m;
            }
            return 0.5 * (a + b);
        };
        // the offset of the first span crosses the axis at s0, at height apex[1]
        const double s0 = root([&](double s) { double q[2], dq[2]; off(0, s, q, dq); return q[0]; }, 0.0, 0.5);
        double apex[2], dapex[2];
        off(0, s0, apex, dapex);
        // pi times the integral of r^2 dy along f from s = a to b (Simpson's rule): the vase, and
        // the cavity inside the offset from the axis to the open top
        auto turned = [](const std::function<void(double, double &, double &)> &f, double a, double b)
        {
            const int n = 4000;
            double sum = 0.0;
            for (int i = 0; i <= n; ++i)
            {
                double rr = 0.0, dy = 0.0;
                f(a + (b - a) * i / n, rr, dy);
                sum += (i == 0 || i == n ? 1 : (i % 2 ? 4 : 2)) * rr * rr * dy;
            }
            return PI * sum * (b - a) / n / 3.0;
        };
        const double s1 = root([&](double s) { double q[2], dq[2]; off(5, s, q, dq); return q[1] - top; }, 0.5, 1.0);
        double vase_volume = 0.0, cavity_volume = 0.0;
        for (int k = 0; k < 6; ++k)
        {
            vase_volume += turned([&](double s, double &rr, double &dy) { double d[3][2]; bez(k, s, d); rr = d[0][0]; dy = d[1][1]; }, 0.0, 1.0);
            cavity_volume += turned([&](double s, double &rr, double &dy) { double q[2], dq[2]; off(k, s, q, dq); rr = q[0]; dy = dq[1]; },
                                    k == 0 ? s0 : 0.0, k == 5 ? s1 : 1.0);
        }

        // the profile as curves: the six spans, the rim, the axis
        std::vector<int64_t> kinds;
        std::vector<double> reals;
        for (int k = 0; k < 6; ++k)
        {
            kinds.insert(kinds.end(), {SMLIB_CURVE_NURBS, 3, 4});
            for (int j = 0; j < 4; ++j)
            {
                reals.insert(reals.end(), {S[k][j][0], S[k][j][1], 0, 1});
            }
            reals.insert(reals.end(), {0, 0, 0, 0, 1, 1, 1, 1});
        }
        kinds.insert(kinds.end(), {SMLIB_CURVE_LINE, SMLIB_CURVE_LINE});
        reals.insert(reals.end(), {6, top, 0, 0, top, 0, 0, top, 0, 0, 0, 0});
        for (double turn : {1.0, -1.0})
        {
            const double origin[3] = {0, 0, 0}, axis[3] = {0, turn, 0};
            smlib_brep *vase = nullptr;
            CHECK(smlib_revolve(kinds.data(), static_cast<int64_t>(kinds.size()), reals.data(), static_cast<int64_t>(reals.size()), origin, axis, 2 * PI, 1,
                                &vase) == SMLIB_OK);
            if (vase == nullptr)
            {
                continue;
            }
            // The bottom's inward offset surface (the side that moves its middle up) ends on the
            // axis: one of its boundary lines is the point apex.
            SmTArray<SmFace *> faces;
            vase->p->GetFaces(faces);
            int bottoms = 0;
            for (ULONG i = 0; i < faces.GetSize(); ++i)
            {
                SmSurface *surf = faces[i]->GetSurface();
                SmExtent3d e;
                CHECK(faces[i]->CalculateBoundingBox(e, TRUE) == SM_SUCCESS);
                if (!surf->IsKindOf(SmSurfOfRevolution_TYPE) || e.GetMin(1) > 1e-9)
                {
                    continue;
                }
                bottoms++;
                for (double d : {t, -t})
                {
                    SmSurface *o = nullptr;
                    CHECK(surf->CreateOffsetSurface(*SmApiGetOrCreateContext(), d, SM_CAST_APPROXTOL3D(1e-6), o) == SM_SUCCESS && o != nullptr);
                    if (o == nullptr)
                    {
                        continue;
                    }
                    const SmExtent2d dom = o->GetNaturalUVDomain();
                    SmPoint3d mid_o, mid_s;
                    CHECK(o->EvaluatePoint(dom.GetMid(), mid_o) == SM_SUCCESS && surf->EvaluatePoint(dom.GetMid(), mid_s) == SM_SUCCESS);
                    if (mid_o.y > mid_s.y)
                    {
                        // the boundary line nearest the axis: its farthest point from it, its height
                        double nearest = HUGE_VAL, height = 0.0;
                        for (int side = 0; side < 4; ++side)
                        {
                            double far = 0.0, y = 0.0;
                            for (int j = 0; j <= 8; ++j)
                            {
                                const double w = j / 8.0;
                                const SmPoint2d uv = side < 2 ? dom.Evaluate(side, w) : dom.Evaluate(w, side - 2);
                                SmPoint3d p;
                                CHECK(o->EvaluatePoint(uv, p) == SM_SUCCESS);
                                far = std::max(far, std::hypot(p.x, p.z));
                                y += p.y / 9.0;
                            }
                            if (far < nearest)
                            {
                                nearest = far;
                                height = y;
                            }
                        }
                        std::printf("vase about %+gy: the bottom's inward offset ends %.3g from the axis at y %.9f (the exact offset crosses at %.9f)\n", turn,
                                    nearest, height, apex[1]);
                        CHECK(nearest < 1e-6 && std::fabs(height - apex[1]) < 1e-5);
                    }
                    delete o;
                }
            }
            CHECK(bottoms == 1);

            // The shell. The bottom's offset face was clipped to its offset surface's domain,
            // trimmed at the axis, by a 2D boolean of planar faces made from the face's UV loops,
            // taken from the inset's faceuse: clockwise, so the face was inside out, the
            // intersection empty, and the bottom got no offset at all (1001, no closed cavity).
            // The shell is the vase less the cavity inside the offset, up to the open top.
            const double rim[3] = {0, top, 3};
            smlib_brep *cup = nullptr;
            CHECK(smlib_shell(vase, rim, 1, 0.05, -t, &cup, nullptr) == SMLIB_OK);
            if (cup != nullptr)
            {
                double v = 0.0;
                CHECK(smlib_brep_volume(cup, 1e-10, &v) == SMLIB_OK);
                std::printf("vase about %+gy shelled 2 open at the top: %.9f, the profile's quadrature %.9f\n", turn, v, vase_volume - cavity_volume);
                CHECK(near(v, vase_volume - cavity_volume, 1e-8));
                done(cup);
            }
            smlib_brep_free(vase);
        }
    }

    // SolidCraft's ladder part 09/011 (tests/data/part-chamfers-shell.txt, as smlib_shell was
    // handed it): a block with a bezier-arched window 1 mm below a round, an M5 hole 0.35 mm from
    // that round, a chamfer, holes and a slot, shelled 1.5 inward open at the bottom. The round's
    // offset is tangent to the wall's offset along z = 32 and to the top's along y = 33, so the
    // window's and the hole's offsets meet them tangentially where they cross: those contacts were
    // found inside the curves instead of at their ends, and their sliver edges failed the merge.
    // The result is checked point by point: an inward offset is sharp at every edge, so each face
    // lies on the part's boundary or 1.5 or more from its faces but the open one, exactly 1.5
    // except near the mitres of the part's concave edges; sampled at one interior point per face
    // (where on a face that point falls decides whether it reads 1.5 or more).
    {
        smlib_flat *f = nullptr;
        smlib_brep *part = nullptr;
        CHECK(smlib_flat_read(SMLIB_TEST_DATA "/part-chamfers-shell.txt", &f) == SMLIB_OK && smlib_brep_import(f, &part) == SMLIB_OK);
        smlib_flat_free(f);
        const double bottom[3] = {40.745625596812403, 30.55585248496563, 0};
        CHECK(part != nullptr && smlib_shell(part, bottom, 1, 0.080157097709407002, -1.5, &r, nullptr) == SMLIB_OK);
        if (r != nullptr)
        {
            SmBrep *walls = nullptr;
            CHECK(SmApiBrepCopy(part->p, walls) == SM_SUCCESS && walls != nullptr);
            walls->m_bEditingEnabled = TRUE;
            SmTArray<SmFace *> wf;
            walls->GetFaces(wf);
            for (ULONG i = 0; i < wf.GetSize(); i++)
            {
                SmExtent3d e;
                CHECK(wf[i]->CalculateBoundingBox(e, TRUE) == SM_SUCCESS);
                if (std::fabs(e.GetMin(2)) < 1e-9 && std::fabs(e.GetMax(2)) < 1e-9)
                {
                    CHECK(walls->DeleteFace(wf[i], TRUE, FALSE) == SM_SUCCESS);
                }
            }
            SmTArray<SmFace *> faces;
            r->p->GetFaces(faces);
            int on_part = 0, at_wall = 0, past_wall = 0;
            for (ULONG i = 0; i < faces.GetSize(); i++)
            {
                SmPoint3d p, q;
                double d_part = -1.0, d_walls = -1.0;
                CHECK(faces[i]->CalculateAnInternalPoint(p) == SM_SUCCESS);
                CHECK(SmApiGetClosestPoint(part->p, p, q, d_part) == SM_SUCCESS);
                CHECK(SmApiGetClosestPoint(walls, p, q, d_walls) == SM_SUCCESS);
                if (d_part <= 1e-9)
                {
                    on_part++;
                }
                else
                {
                    CHECK(d_walls >= 1.5 - 1e-5);
                    (std::fabs(d_walls - 1.5) <= 1e-5 ? at_wall : past_wall)++;
                }
            }
            std::printf("part 09/011 shell: %d faces on the part, %d at the wall, %d past it (at a mitre)\n", on_part, at_wall, past_wall);
            CHECK(on_part == 29 && at_wall + past_wall == 26 && at_wall > past_wall);
            delete walls;
            done(r);
        }
        // Shelled 1.0, the window's offset arch rises exactly to the wall's offset top edge,
        // where the wall meets the round's offset (the window's top is 1 below the round): the
        // arch touches that edge tangentially at (30, 31, 32) and the walls meet along a line,
        // SMLIB_ERR_WALLS_TOUCH (it was SM_ERR_BAD_INTERSECTIONS).
        CHECK(part != nullptr && smlib_shell(part, bottom, 1, 0.080157097709407002, -1.0, &r, nullptr) == SMLIB_ERR_WALLS_TOUCH && r == nullptr);
        smlib_brep_free(part);
    }

    // tangency-labs/solidcore#865: shelled closed (no face open), the plates above whose hole's
    // offset crosses a side failed (1001), as did part 09/011; its whole offset came back with
    // no volume. Before the merge, SmOffsetExecutive::RemoveExcessFaces took every offset face
    // with an edge nearer the body than the offset distance, but the offset faces are not
    // trimmed yet then: the side's offset runs on to where the floor's offset meets it, past
    // where the hole's offset cuts it, so it has such an edge and went, with the hole's and the
    // floor's (part 09/011 lost its five faces round the window), and the merged offset closed
    // nothing. The test is made on the merged offset now, whose faces lie wholly on the offset
    // or wholly in the excess. A closed plate's cavity is the inset box less the hole's offset
    // disc, over 10 - 2t; the part's whole offset is the part less its closed shell.
    {
        auto band = [](double rr, double a, double b)
        {
            auto F = [rr](double u)
            {
                u = std::max(-rr, std::min(rr, u));
                return u * std::sqrt(rr * rr - u * u) + rr * rr * std::asin(u / rr);
            };
            return F(b) - F(a);
        };
        struct Plate
        {
            double x, y, d, t;
        };
        for (Plate p : {Plate{4.5, 20, 40, 1.5}, Plate{4.5, 20, 40, 1.0}, Plate{20, 4, 8, 1.5}, Plate{12, 20, 40, 1.5}})
        {
            const double o[3] = {0, 0, 0}, size[3] = {40, p.d, 10}, axis[3] = {p.x, p.y, -1};
            smlib_brep *block = nullptr, *pin = nullptr, *plate = nullptr;
            CHECK(smlib_make_box(o, size, &block) == SMLIB_OK && smlib_make_cylinder(axis, 2.0, 12.0, &pin) == SMLIB_OK);
            CHECK(smlib_boolean(block, pin, SMLIB_BOOLEAN_SUBTRACT, &plate) == SMLIB_OK);
            const double rr = 2.0 + p.t;
            double inside = band(rr, std::max(p.t - p.y, -rr), std::min(p.d - p.t - p.y, rr));
            if (p.x - rr < p.t)
            {
                inside -= band(rr, -rr, p.t - p.x);
            }
            const double cavity = ((40.0 - 2.0 * p.t) * (p.d - 2.0 * p.t) - inside) * (10.0 - 2.0 * p.t);
            CHECK(plate != nullptr && smlib_shell(plate, nullptr, 0, 0.1, -p.t, &r, nullptr) == SMLIB_OK);
            if (r != nullptr)
            {
                CHECK(near(volume_of(r), 40.0 * p.d * 10.0 - PI * 4.0 * 10.0 - cavity, 1e-7));
                done(r);
            }
            for (smlib_brep *b : {block, pin, plate})
            {
                smlib_brep_free(b);
            }
        }

        smlib_flat *f = nullptr;
        smlib_brep *part = nullptr;
        CHECK(smlib_flat_read(SMLIB_TEST_DATA "/part-chamfers-shell.txt", &f) == SMLIB_OK && smlib_brep_import(f, &part) == SMLIB_OK);
        smlib_flat_free(f);
        CHECK(part != nullptr && smlib_shell(part, nullptr, 0, 0.080157097709407002, -1.5, &r, nullptr) == SMLIB_OK);
        if (r != nullptr)
        {
            double v_part = 0.0, v_shell = 0.0, v_offset = 0.0;
            CHECK(smlib_brep_volume(part, 1e-10, &v_part) == SMLIB_OK && smlib_brep_volume(r, 1e-10, &v_shell) == SMLIB_OK);
            SmBrep *w = nullptr, *o = nullptr;
            CHECK(SmApiBrepCopy(part->p, w) == SM_SUCCESS && SmApiOffsetBrepFull(w, -1.5, TRUE, TRUE, o) == SM_SUCCESS);
            CHECK(o != nullptr && o->IsManifoldSolid());
            // (the offset leaves its input allocated beside its result, as smlib_thicken notes)
            if (w != o)
            {
                delete w;
            }
            if (o != nullptr)
            {
                smlib_brep inset{o};
                CHECK(smlib_brep_volume(&inset, 1e-10, &v_offset) == SMLIB_OK);
                std::printf("part 09/011 shelled closed: %.9f; its offset %.9f, the part less the shell %.9f\n", v_shell, v_offset, v_part - v_shell);
                CHECK(near(v_offset, v_part - v_shell, 1e-6));
                delete o;
            }
            done(r);
        }
        smlib_brep_free(part);
    }

    // SolidCraft's cross-hole cut (tests/data/cross-pin.txt less cross-hole.txt, as SolidCraft
    // handed them to smlib_boolean): a pin of radius 2 and height 6 about z, a hole of radius 1
    // through it along y, both seamed at 0.6131 rad. Each intersection loop is a closed curve that
    // touches the hole's seam line at its start, which is its end: moving that contact to the
    // loop's other end passed a one-way measurement of the whole turn against the line's one
    // point, and the cut failed (1001). The cavity is 24 pi less the integral of
    // 4 sqrt(1 - x^2) sqrt(4 - x^2) over [-1, 1] (x = sin a: 4 cos^2 a sqrt(4 - sin^2 a), smooth,
    // by Simpson's rule).
    {
        smlib_flat *fa = nullptr, *fb = nullptr;
        smlib_brep *pin = nullptr, *hole = nullptr;
        CHECK(smlib_flat_read(SMLIB_TEST_DATA "/cross-pin.txt", &fa) == SMLIB_OK && smlib_brep_import(fa, &pin) == SMLIB_OK);
        CHECK(smlib_flat_read(SMLIB_TEST_DATA "/cross-hole.txt", &fb) == SMLIB_OK && smlib_brep_import(fb, &hole) == SMLIB_OK);
        smlib_flat_free(fa);
        smlib_flat_free(fb);
        const int n = 2000;
        double lens = 0.0;
        for (int i = 0; i <= n; i++)
        {
            const double a = -PI / 2 + PI * i / n, c = std::cos(a), s = std::sin(a);
            lens += (i == 0 || i == n ? 1.0 : (i % 2 ? 4.0 : 2.0)) * 4.0 * c * c * std::sqrt(4.0 - s * s);
        }
        lens *= PI / n / 3.0;
        CHECK(pin != nullptr && hole != nullptr && smlib_boolean(pin, hole, SMLIB_BOOLEAN_SUBTRACT, &r) == SMLIB_OK);
        if (r != nullptr)
        {
            CHECK(near(volume_of(r), 24.0 * PI - lens, 1e-6));
            done(r);
        }
        smlib_brep_free(pin);
        smlib_brep_free(hole);
    }

    // Thicken: a sheet (the cube's top face alone, SMLib's own way of making one is not in the
    // C API; a solid is refused).
    CHECK(smlib_thicken(cube, 1.0, 0, &r) == SMLIB_ERR_ARG);
    CHECK(smlib_thicken(nullptr, 1.0, 0, &r) == SMLIB_ERR_ARG);
    CHECK(smlib_thicken(cube, NAN, 0, &r) == SMLIB_ERR_ARG);

    // Rig case shell/5326: an L-shaped sheet, two walls 30 and 20 long at a right angle and 20
    // deep, thickened 2 into its inside corner and 2 out of it: an L-shaped prism, 8 faces. The
    // thicken builds one side wall per edge of the sheet, and the two along each end of the L lie
    // in one plane: they came back as two faces, split along the edge from the corner to its
    // offset (10 faces). The walls now come out as one face where consecutive walls lie on one
    // surface, each put on its plane first (a wall beside the corner is a trapezoid, a bilinear
    // patch), so the volumes are exact: (30 + 20 -+ 2) * 2 * 20. The sheet's own splits stay: a
    // flat sheet of two faces side by side (a profile of two lines on one line, extruded)
    // thickened 1 keeps both faces and their offsets apart, its walls along the split edges one
    // face each: 8 faces, 20 * 20 * 1.
    {
        const int64_t kinds[2] = {SMLIB_CURVE_LINE, SMLIB_CURVE_LINE};
        const double ell[12] = {0, 0, 0, 30, 0, 0, 30, 0, 0, 30, 20, 0}, split[12] = {0, 0, 0, 10, 0, 0, 10, 0, 0, 20, 0, 0};
        const double up[3] = {0, 0, 1};
        for (const double *profile : {ell, split})
        {
            smlib_brep *sheet = nullptr;
            CHECK(smlib_extrude(kinds, 2, profile, 12, up, 20.0, 0.0, 0, 0, &sheet, nullptr) == SMLIB_OK);
            for (double t : {-2.0, 2.0})
            {
                if (profile == split)
                {
                    t /= 2.0;
                }
                CHECK(sheet != nullptr && smlib_thicken(sheet, t, 0, &r) == SMLIB_OK);
                if (r == nullptr)
                {
                    continue;
                }
                int64_t counts[3] = {0, 0, 0};
                double v = 0.0;
                CHECK(smlib_brep_counts(r, counts) == SMLIB_OK && smlib_brep_volume(r, 1e-10, &v) == SMLIB_OK);
                const double want = profile == ell ? (50.0 + (t < 0.0 ? -2.0 : 2.0)) * 2.0 * 20.0 : 400.0;
                std::printf("%s sheet thickened %+g: %lld faces, volume %.9f\n", profile == ell ? "L" : "split", t, static_cast<long long>(counts[0]), v);
                CHECK(counts[0] == 8 && near(v, want, 1e-12));
                done(r);
            }
            smlib_brep_free(sheet);
        }
    }

    // Rig case surfaces/cov-thicken-cylinder-sheet-symmetric-tube: SolidCraft's cylinder wall,
    // radius 10 and 20 high (tests/data/cylinder-sheet-r10.txt, one face on a cone of semi-angle
    // 0, as SolidCraft hands it over), thickened 2: a tube, 4 faces, its ends flat rings. Each end
    // is the wall between the sheet's circle and its offset, which lie in one plane; it came back
    // a surface of revolution of the line between them (SmSurfOfRevolution), which SolidCraft
    // reads as a B-spline. Such a wall is now on its plane (SmPlane), bounded by its two circles
    // (its seam, which the closed surface needed, gone). Symmetric, the tube is r 9 to
    // 11: 800 pi; by +2, r 10 to 12: 880 pi; by -2, r 8 to 10: 720 pi.
    {
        smlib_flat *f = nullptr;
        smlib_brep *sheet = nullptr;
        CHECK(smlib_flat_read(SMLIB_TEST_DATA "/cylinder-sheet-r10.txt", &f) == SMLIB_OK && smlib_brep_import(f, &sheet) == SMLIB_OK);
        smlib_flat_free(f);
        struct Thicken
        {
            double t;
            int symmetric;
            double volume;
        };
        const Thicken cases[] = {{2.0, 1, 800.0 * PI}, {2.0, 0, 880.0 * PI}, {-2.0, 0, 720.0 * PI}};
        for (const Thicken &c : cases)
        {
            CHECK(sheet != nullptr && smlib_thicken(sheet, c.t, c.symmetric, &r) == SMLIB_OK);
            if (r == nullptr)
            {
                continue;
            }
            SmTArray<SmFace *> faces;
            r->p->GetFaces(faces);
            int planes = 0, round = 0;
            for (ULONG i = 0; i < faces.GetSize(); i++)
            {
                const SmSurface *s = faces[i]->GetSurface();
                // a flat ring on its plane is bounded by its two circles, no seam between them
                SmTArray<SmLoop *> loops;
                faces[i]->GetLoops(loops);
                planes += s->IsKindOf(SmPlane_TYPE) && loops.GetSize() == 2 ? 1 : 0;
                round += s->IsKindOf(SmCone_TYPE) || s->IsKindOf(SmCylinder_TYPE) ? 1 : 0;
            }
            double v = volume_of(r);
            std::printf("cylinder sheet thickened %+g%s: %lu faces, %d planes, %d cylinders, volume %.9f\n", c.t,
                        c.symmetric ? " symmetric" : "", static_cast<unsigned long>(faces.GetSize()), planes, round, v);
            // (to the volume's own accuracy for faces of revolution, 1e-9)
            CHECK(faces.GetSize() == 4 && planes == 2 && round == 2 && near(v, c.volume, 1e-9));
            done(r);
        }
        smlib_brep_free(sheet);
    }

    // smlib_offset_sheet of faces that close on themselves: a sphere's one face (radius 10) moved
    // out by 2 is the sphere of radius 12, a closed solid (it was refused as not a sheet); moved in
    // to its centre, refused.
    {
        const double centre[3] = {0, 0, 0};
        smlib_brep *ball = nullptr;
        CHECK(smlib_make_sphere(centre, 10.0, &ball) == SMLIB_OK);
        CHECK(ball != nullptr && smlib_offset_sheet(ball, 2.0, &r) == SMLIB_OK);
        if (r != nullptr)
        {
            int32_t closed = 0;
            CHECK(smlib_brep_is_manifold_solid(r, &closed) == SMLIB_OK && closed == 1);
            std::printf("sphere offset out by 2: volume %.9f\n", volume_of(r));
            CHECK(near(volume_of(r), 4.0 / 3.0 * PI * 1728.0, 1e-6));
            done(r);
        }
        CHECK(smlib_offset_sheet(ball, -10.0, &r) != SMLIB_OK && r == nullptr);
        smlib_brep_free(ball);
    }

    // A sheet on a B-spline with simple inner knots (a gentle cubic wave of 6 control points,
    // knots 0 0 0 0 1/3 2/3 1 1 1 1, extruded 10 along y) offset by 0.5: its faces are on SMLib's offset
    // surface, which export approximates as cubic Bezier patches. The approximation took the
    // root B-spline's simple knots for full ones, ran its patches past the control points and
    // answered success with no surface: the export failed (1001).
    {
        const int64_t kinds[3] = {SMLIB_CURVE_NURBS, 3, 6};
        const double reals[34] = {0, 0, 0, 1, 2, 0, 0.4, 1, 4, 0, -0.4, 1, 6, 0, 0.4, 1, 8, 0, -0.4, 1, 10, 0, 0, 1,
                                       0, 0, 0, 0, 1.0 / 3.0, 2.0 / 3.0, 1, 1, 1, 1};
        const double along[3] = {0, 1, 0};
        smlib_brep *wave = nullptr;
        CHECK(smlib_extrude(kinds, 3, reals, 34, along, 10.0, 0.0, 0, 0, &wave, nullptr) == SMLIB_OK);
        CHECK(wave != nullptr && smlib_offset_sheet(wave, 0.5, &r) == SMLIB_OK);
        if (r != nullptr)
        {
            smlib_flat *flat = nullptr;
            CHECK(smlib_brep_export(r, &flat) == SMLIB_OK && flat != nullptr);
            smlib_flat_free(flat);
            done(r);
        }
        smlib_brep_free(wave);
    }

    // smlib_untrim_face in place on a sheet: a triangle's planar face back to the box of its span,
    // the face replaced in the body. Deleting the old face is an edit of the brep, which SMLib
    // refused outside edit mode (1001); the face is the one the copy out gives.
    {
        const int64_t kinds[3] = {SMLIB_CURVE_LINE, SMLIB_CURVE_LINE, SMLIB_CURVE_LINE};
        const double tri[18] = {0, 0, 0, 4, 0, 0, 4, 0, 0, 0, 3, 0, 0, 3, 0, 0, 0, 0};
        smlib_brep *sheet = nullptr;
        CHECK(smlib_planar_faces(kinds, 3, tri, 18, &sheet) == SMLIB_OK);
        const double inside[3] = {1.0, 1.0, 0.0};
        smlib_brep *apart = nullptr;
        smlib_brep *placed = nullptr;
        CHECK(sheet != nullptr && smlib_untrim_face(sheet, inside, 0.1, nullptr, -1, 1, 0, &apart) == SMLIB_OK);
        CHECK(sheet != nullptr && smlib_untrim_face(sheet, inside, 0.1, nullptr, -1, 1, 1, &placed) == SMLIB_OK);
        double a_apart = -1.0, a_placed = -2.0;
        int64_t counts[3] = {0, 0, 0};
        CHECK(apart != nullptr && smlib_brep_area(apart, 1e-9, &a_apart) == SMLIB_OK);
        CHECK(placed != nullptr && smlib_brep_area(placed, 1e-9, &a_placed) == SMLIB_OK && smlib_brep_counts(placed, counts) == SMLIB_OK);
        std::printf("triangle untrimmed: copied out %.9f, in place %.9f\n", a_apart, a_placed);
        CHECK(counts[0] == 1 && a_apart > 6.0 && near(a_placed, a_apart, 1e-12));
        smlib_brep_free(apart);
        smlib_brep_free(placed);
        smlib_brep_free(sheet);
    }

    // A sphere's cap thickened: the sphere of radius 10 split by the plane z = 5, the face on the
    // sphere above it copied out alone (its rim two half circles and a seam), thickened out by 1:
    // the sector of the shell between radii 10 and 11 over solid angle pi, pi/3 * (11^3 - 10^3).
    // Each half of the rim's offset lies on the whole offset circle, and its second trim point
    // was taken at that circle's seam, the first point again: each half took the whole circle,
    // the loop had a gap and the thicken failed (1001).
    {
        const double centre[3] = {0, 0, 0};
        smlib_brep *ball = nullptr;
        CHECK(smlib_make_sphere(centre, 10.0, &ball) == SMLIB_OK);
        const double o[3] = {0, 0, 5}, x[3] = {1, 0, 0}, y[3] = {0, 1, 0}, u[2] = {-50, 50};
        smlib_surface *cut = nullptr;
        CHECK(smlib_surface_plane(o, x, y, u, u, &cut) == SMLIB_OK);
        smlib_brep *below = nullptr, *above = nullptr;
        CHECK(ball != nullptr && cut != nullptr && smlib_split_body(ball, cut, &below, &above) == SMLIB_OK);
        smlib_surface_free(cut);
        SmBrep *cap = nullptr;
        if (above != nullptr)
        {
            SmTArray<SmFace *> faces;
            above->p->GetFaces(faces);
            SmTArray<SmFace *> one;
            for (ULONG i = 0; i < faces.GetSize(); ++i)
            {
                if (faces[i]->GetSurface()->IsKindOf(SmSphere_TYPE))
                {
                    one.Add(faces[i]);
                }
            }
            CHECK(one.GetSize() == 1);
            cap = new (*SmApiGetOrCreateContext()) SmBrep();
            cap->SetTolerance(above->p->GetTolerance());
            SmTArray<SmFace *> copied;
            CHECK(above->p->CopyFaces(one, cap, &copied, TRUE) == SM_SUCCESS && copied.GetSize() == 1);
        }
        smlib_brep sheet{cap};
        CHECK(cap != nullptr && smlib_thicken(&sheet, 1.0, 0, &r) == SMLIB_OK);
        if (r != nullptr)
        {
            const double want = PI / 3.0 * (1331.0 - 1000.0);
            std::printf("sphere cap thickened 1: volume %.9f, want %.9f\n", volume_of(r), want);
            CHECK(near(volume_of(r), want, 1e-6));
            done(r);
        }
        delete cap;
        for (smlib_brep *b : {ball, below, above})
        {
            smlib_brep_free(b);
        }
    }

    // A cube 0..10 cut by a sheet of two faces meeting at a right angle (the top z = 5 and the
    // side x = 5 of the box -20..5 x -20..20 x -20..5, copied off it): the corner x < 5, z < 5
    // (250) behind the sheet, the rest (750) in front, whichever face the sides are read from;
    // a body that is not a sheet is refused.
    {
        smlib_brep *cube = box(0, 0, 0, 10, 10, 10);
        smlib_brep *slab = box(-20, -20, -20, 25, 40, 25);
        SmTArray<SmFace *> faces, two;
        slab->p->GetFaces(faces);
        for (ULONG i = 0; i < faces.GetSize(); ++i)
        {
            SmPoint3d p;
            if (SmApiFaceInternalPoint(faces[i], p) == SM_SUCCESS && (std::fabs(p.z - 5.0) < 1e-9 || std::fabs(p.x - 5.0) < 1e-9))
            {
                two.Add(faces[i]);
            }
        }
        CHECK(two.GetSize() == 2);
        SmBrep *l = new (*SmApiGetOrCreateContext()) SmBrep();
        l->SetTolerance(slab->p->GetTolerance());
        SmTArray<SmFace *> copied;
        CHECK(slab->p->CopyFaces(two, l, &copied, TRUE) == SM_SUCCESS && copied.GetSize() == 2 && l->IsSheet());
        smlib_brep sheet{l};
        for (int64_t face = 0; face < 2; ++face)
        {
            smlib_brep *back = nullptr, *front = nullptr;
            CHECK(smlib_split_body_by_sheet(cube, &sheet, face, &back, &front) == SMLIB_OK && back != nullptr && front != nullptr);
            if (back != nullptr && front != nullptr)
            {
                const double vb = volume_of(back), vf = volume_of(front);
                std::printf("cube split by an L sheet (face %lld): %.9f + %.9f\n", static_cast<long long>(face), vb, vf);
                CHECK(near(std::fmin(vb, vf), 250.0, 1e-9) && near(std::fmax(vb, vf), 750.0, 1e-9));
                CHECK(sound(back) && sound(front));
            }
            smlib_brep_free(back);
            smlib_brep_free(front);
        }
        smlib_brep *back = nullptr, *front = nullptr;
        CHECK(smlib_split_body_by_sheet(cube, slab, 0, &back, &front) == SMLIB_ERR_ARG && back == nullptr);
        CHECK(smlib_split_body_by_sheet(cube, &sheet, 2, &back, &front) == SMLIB_ERR_ARG && back == nullptr);
        // The top face alone ends at x = 5, inside the cube: it nicks it and cuts nothing.
        SmTArray<SmFace *> top;
        for (ULONG i = 0; i < two.GetSize(); ++i)
        {
            SmPoint3d p;
            if (SmApiFaceInternalPoint(two[i], p) == SM_SUCCESS && std::fabs(p.z - 5.0) < 1e-9)
            {
                top.Add(two[i]);
            }
        }
        SmBrep *lid = new (*SmApiGetOrCreateContext()) SmBrep();
        lid->SetTolerance(slab->p->GetTolerance());
        SmTArray<SmFace *> one;
        CHECK(top.GetSize() == 1 && slab->p->CopyFaces(top, lid, &one, TRUE) == SM_SUCCESS && one.GetSize() == 1);
        smlib_brep nick{lid};
        CHECK(smlib_split_body_by_sheet(cube, &nick, 0, &back, &front) == SMLIB_ERR_NO_CUT && back == nullptr);
        delete lid;
        delete l;
        smlib_brep_free(cube);
        smlib_brep_free(slab);
    }

    // Two holes 6 apart cut one after the other through a disc: the second hole's wall meets the
    // first's (trimmed to the disc) in two lines. Parallel cylinders are intersected by
    // SmCone::IntersectWithCone, which once took the lines' span from the whole cylinders, so
    // the lines started below the trimmed wall and were dropped (the boolean failed, 1001).
    {
        smlib_brep *disc = cylinder(0, 0, 0, 20, 10);
        smlib_brep *first = cylinder(6, 0, -1, 5, 12);
        smlib_brep *second = cylinder(3, 3 * std::sqrt(3.0), -1, 5, 12);
        smlib_brep *one = nullptr;
        smlib_brep *two = nullptr;
        CHECK(smlib_boolean(disc, first, SMLIB_BOOLEAN_SUBTRACT, &one) == SMLIB_OK && sound(one));
        CHECK(smlib_boolean(one, second, SMLIB_BOOLEAN_SUBTRACT, &two) == SMLIB_OK && sound(two));
        // The union of two r5 discs 6 apart: twice a disc less their lens.
        const double lens = 2.0 * 25.0 * std::acos(0.6) - 3.0 * 8.0;
        CHECK(near(volume_of(two), (PI * 400.0 - (2.0 * PI * 25.0 - lens)) * 10.0, 1e-6));
        for (smlib_brep *b : {disc, first, second, one, two})
        {
            smlib_brep_free(b);
        }
    }

    for (smlib_brep *b : {cube, rod, holed, pin, round, cut, big, far, touching})
    {
        smlib_brep_free(b);
    }

    analytic_step_inversion();
    shelled_dome_rim();
    chamfered_hole_rim_is_a_cone();
    shelled_rounded_box_rim();
    shelled_single_round_stays_a_cylinder();
    shelled_box_with_a_drilled_hole();
    shell_thicker_than_a_step();
    offset_faces_close_a_step_wall();
    drafted_cylinders_are_cones();
    lip_swept_flush_round_a_corner();
    knife_edge_rounds();
    knife_edge_chamfers();
    overflow_goes_first();
    crease_rounds_from_half_a_degree();
    corner_setback();
    corner_setback_open();
    setback_sections_are_exact_arcs();
    concave_corner_setback();
    corner_setback_eased();
    corner_setback_eased_long_setbacks();
    corner_setback_eased_face_curves_are_quintics();
    shelled_tapered_box_rounds();
    shelled_disc_round_thinner_than_the_wall();
    groove_arc_round_ends_in_the_top();
    corner_ball_after_a_round_is_exactly_singular();
    rounds_at_a_concave_face_corner_are_mitred();
    basin_rim_rounds_mitred_with_a_gap_at_the_corner();
    end_face_rounds_mitred_at_concave_corners();
    chamfer_past_a_sloped_basin_rim();
    chamfer_of_a_rim_loop_past_its_rim();
    rounds_mitred_to_one_past_a_rim();
    closed_rims_blended_past_their_walls();
    big_round_of_a_rounded_rectangle_is_refused();
    tabbed_plate_rounds_mitred_on_the_arcs_continuation();
    corner_patch_edge_lies_on_its_patch();
    corner_patch_drop_from_a_guess_is_its_foot();
    revolved_circle_extends_as_its_circle();
    chamfer_ending_on_a_boss_wall_above_it();
    vertex_fans_split_and_joined();
    rails_trimmed_next_to_a_knot();
    concave_chamfers_at_a_vertex_the_body_touches_itself_at();
    chamfer_of_a_shelled_box_corner_ending_under_its_lip();

    // SolidCraft surface.offset of a ball's one face: the face taken off as a sheet is a closed
    // shell whose regions are all voids, which smlib_offset_sheet refused as a solid
    // (IsManifoldSolid). Offset 2 out it is the sphere of radius 22, area 4 pi 22^2. (The ball
    // itself offsets too: see the sphere's one face above.)
    {
        const double c[3] = {0, 0, 0};
        smlib_brep *ball = nullptr, *sheet = nullptr, *moved = nullptr;
        CHECK(smlib_make_sphere(c, 20.0, &ball) == SMLIB_OK && ball != nullptr);
        CHECK(ball != nullptr && smlib_brep_copy(ball, &sheet) == SMLIB_OK && sheet != nullptr);
        if (sheet != nullptr)
        {
            SmTArray<SmRegion *> regions;
            sheet->p->GetRegions(regions);
            for (ULONG i = 0; i < regions.GetSize(); ++i)
            {
                regions[i]->SetIsVoid(TRUE);
            }
        }
        CHECK(sheet != nullptr && smlib_offset_sheet(sheet, 2.0, &moved) == SMLIB_OK && moved != nullptr);
        if (moved != nullptr)
        {
            double area = 0.0;
            CHECK(smlib_brep_area(moved, 1e-8, &area) == SMLIB_OK);
            std::printf("closed sheet of a ball offset 2: area %.9f, expected %.9f\n", area, 4.0 * PI * 22.0 * 22.0);
            CHECK(near(area, 4.0 * PI * 22.0 * 22.0, 1e-7));
        }
        for (smlib_brep *x : {ball, sheet, moved})
        {
            smlib_brep_free(x);
        }
    }

    // Two variable rounds 0.5 to 2 meeting at a box's corner (SM_API test TestSmVariableRadiusFillet
    // rounds four top edges so). The rails' meeting point was found with both radii held at their
    // end value 2, at y = 8 where the law gives 1.7, and the trace from that rail end, its spans
    // held to the fillet solved inside them, found none it could take (1001). Each round alone
    // takes (1 - pi/4) times the integral of r^2, (1 - pi/4) 17.5; together they take more than
    // one and less than both, which share the corner.
    {
        const double one = (1.0 - PI / 4.0) * 17.5;
        const double removed = variable_rounds_meeting_volume();
        std::printf("variable rounds meeting at a corner: removed %.9g (one round %.9g)\n", removed, one);
        CHECK(removed > one && removed < 2.0 * one);
    }

    // A three-edge corner of variable rounds (constant law) builds at any size: the box 10 with
    // r 1 and the box 0.01 with r 0.001, the same corner 1000 times smaller. sm_CoonsSideAngle
    // measured the corner patch's normal along its sides with Unitize, whose absolute floor (a
    // length of 1e-12) refused a normal near a collapsing side of the small patch that its own
    // test, relative to the derivatives, had found: the small corner failed
    // (SM_ERR_INVALID_INPUT) and the round with it. SMLib's fitting tolerances are lengths, not
    // shares of the part, so the small corner's removed volume, scaled by 1e9, is the big one's
    // only to them (0.5% at r 0.001, 6.3224 against 6.2914).
    {
        const double big = variable_corner_removed(10.0, 1.0);
        const double small = variable_corner_removed(0.01, 0.001);
        std::printf("variable-round corner: removed %.9g at size 10, %.9g x 1e9 at size 0.01\n", big, small * 1e9);
        CHECK(big > 0.0 && small > 0.0);
        CHECK(near(small * 1e9, big, 1e-2));
    }

    if (failures != 0)
    {
        std::fprintf(stderr, "%d check(s) failed\n", failures);
        return 1;
    }
    std::printf("ops_kernel: all checks passed\n");
    return 0;
}
