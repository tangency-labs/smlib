// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

// Flat B-rep arrays through the C API: for a set of solids, export -> import -> compare volume,
// counts and manifoldness, then export again and compare the arrays with the first export.
// Then hostile arrays: truncated, out of range, NaN, huge counts, random mutations; each must
// return an error status, never crash.
//
// The solids are built with SM_API directly and handed to the C API as native .smb files.
//
//   smlib_flat_roundtrip <scratch dir>

#include "smlib.h"

#include <SmApiBrep.h>
#include <SmApiGeneral.h>
#include <SmApiImportExport.h>
#include <SmApiPrimitives.h>
#include <SmApiQueries.h>
#include <SmBrep.h>
#include <SmEdge.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <limits>
#include <map>
#include <string>
#include <vector>

namespace
{

int failures = 0;
std::string dir = ".";

#define CHECK(cond)                                                                                \
    do                                                                                             \
    {                                                                                              \
        if (!(cond))                                                                               \
        {                                                                                          \
            std::fprintf(stderr, "%s:%d: CHECK(%s) failed\n", __FILE__, __LINE__, #cond);          \
            failures++;                                                                            \
        }                                                                                          \
    } while (0)

// Every id: int64 ids 1..35, double ids 101..165.
std::vector<int32_t> IntIds()
{
    std::vector<int32_t> v;
    for (int32_t id = SMLIB_FLAT_REGION_SHELL_COUNT; id <= SMLIB_FLAT_WIRE_EDGE_NURBS_FORM; ++id)
    {
        v.push_back(id);
    }
    return v;
}

std::vector<int32_t> RealIds()
{
    std::vector<int32_t> v;
    for (int32_t id = SMLIB_FLAT_TOLERANCE; id <= SMLIB_FLAT_VERTEX_TOLERANCE; ++id)
    {
        v.push_back(id);
    }
    return v;
}

// A copy of every array of a flat set.
struct Arrays
{
    std::map<int32_t, std::vector<int64_t>> i;
    std::map<int32_t, std::vector<double>> f;
};

Arrays Snapshot(const smlib_flat *flat)
{
    Arrays a;
    for (int32_t id : IntIds())
    {
        const int64_t *d = nullptr;
        int64_t n = 0;
        CHECK(smlib_flat_get_i64(flat, id, &d, &n) == SMLIB_OK);
        a.i[id] = std::vector<int64_t>(d, d + n);
    }
    for (int32_t id : RealIds())
    {
        const double *d = nullptr;
        int64_t n = 0;
        CHECK(smlib_flat_get_f64(flat, id, &d, &n) == SMLIB_OK);
        a.f[id] = std::vector<double>(d, d + n);
    }
    return a;
}

smlib_flat *FromArrays(const Arrays &a)
{
    smlib_flat *flat = smlib_flat_new();
    for (auto &kv : a.i)
    {
        CHECK(smlib_flat_set_i64(flat, kv.first, kv.second.data(), static_cast<int64_t>(kv.second.size())) == SMLIB_OK);
    }
    for (auto &kv : a.f)
    {
        CHECK(smlib_flat_set_f64(flat, kv.first, kv.second.data(), static_cast<int64_t>(kv.second.size())) == SMLIB_OK);
    }
    return flat;
}

// Ints exactly, reals within 1e-9 relative (1e-9 absolute below magnitude 1). Returns the
// largest relative real difference, or infinity on a mismatch of sizes or ints.
double Compare(const Arrays &a, const Arrays &b, const char *what)
{
    double worst = 0.0;
    for (auto &kv : a.i)
    {
        if (kv.second != b.i.at(kv.first))
        {
            std::fprintf(stderr, "  %s: int array %d differs\n", what, kv.first);
            return INFINITY;
        }
    }
    for (auto &kv : a.f)
    {
        const std::vector<double> &x = kv.second, &y = b.f.at(kv.first);
        // EXTENT is SMLib's loose (non-tight) bounding box, informational and ignored by import;
        // it can differ between a B-rep and its re-import (blend1: the blend face's box comes
        // out 4.6e-5 wider), with an all-NURBS export too. Its length still has to match.
        if (kv.first == SMLIB_FLAT_EXTENT && x.size() == y.size())
        {
            continue;
        }
        // An edge's or vertex's own tolerance is what SMLib finds it needs: building the edge
        // again from the exported curves may need more (blend1: 2e-5 then 1.1e-4), never less.
        if ((kv.first == SMLIB_FLAT_EDGE_TOLERANCE || kv.first == SMLIB_FLAT_VERTEX_TOLERANCE) && x.size() == y.size())
        {
            for (size_t k = 0; k < x.size(); ++k)
            {
                if (y[k] < x[k] * (1.0 - 1e-9))
                {
                    std::fprintf(stderr, "  %s: real array %d [%zu]: %.17g then less, %.17g\n", what, kv.first, k, x[k], y[k]);
                    return INFINITY;
                }
            }
            continue;
        }
        if (x.size() != y.size())
        {
            std::fprintf(stderr, "  %s: real array %d has %zu then %zu values\n", what, kv.first, x.size(), y.size());
            return INFINITY;
        }
        for (size_t k = 0; k < x.size(); ++k)
        {
            double d = std::fabs(x[k] - y[k]) / std::max({1.0, std::fabs(x[k]), std::fabs(y[k])});
            if (d > worst)
            {
                worst = d;
            }
            if (d > 1e-9)
            {
                std::fprintf(stderr, "  %s: real array %d [%zu]: %.17g vs %.17g\n", what, kv.first, k, x[k], y[k]);
                return INFINITY;
            }
        }
    }
    return worst;
}

// Hands an SM_API B-rep to the C API through a native file.
smlib_brep *ToCApi(SmBrep *b, const char *name)
{
    if (b == nullptr)
    {
        return nullptr;
    }
    std::string path = dir + "/" + name + ".smb";
    CHECK(SmApiWriteBrepToFile(b, path.c_str(), FALSE) == SM_SUCCESS);
    delete b;
    smlib_brep *r = nullptr;
    CHECK(smlib_brep_read(path.c_str(), 0, &r) == SMLIB_OK);
    return r;
}

struct Facts
{
    int64_t counts[3] = {0, 0, 0};
    double volume = 0.0;
    int32_t manifold = -1;
};

Facts Measure(const smlib_brep *b)
{
    Facts f;
    CHECK(smlib_brep_counts(b, f.counts) == SMLIB_OK);
    CHECK(smlib_brep_volume(b, 1e-7, &f.volume) == SMLIB_OK);
    CHECK(smlib_brep_is_manifold_solid(b, &f.manifold) == SMLIB_OK);
    return f;
}

// Expected curve kinds per solid (nurbs, line, circle, ellipse). Export keeps SMLib's own curve
// classes: SmApiCreateCylinder, SmApiCreateCone and SmApiCreateTorus make their rims and seam
// circles SmCircles, and so does a boolean that cuts a cylinder with a plane square to its axis;
// a boolean's other intersection curves are B-splines.
const std::map<std::string, std::vector<int64_t>> g_kinds = {
    {"box", {0, 12, 0, 0}},      // twelve SmLines
    {"cylinder", {0, 1, 2, 0}},  // SmCircle rims, a line seam
    {"cone", {0, 1, 2, 0}},      // SmCircle rims, a line seam
    {"sphere", {0, 0, 1, 0}},    // the seam: an SmCircle half circle
    {"torus", {0, 0, 2, 0}},     // SmCircle seam circles
    {"fillet", {0, 24, 24, 0}},  // the box's trimmed lines, the fillets' SmCircle arcs
    {"blend", {48, 0, 0, 0}},    // blend-curve fillets: all B-splines
    {"blend1", {4, 11, 0, 0}},   // the blend face's four B-spline edges, the box's eleven lines
    {"boolean", {0, 13, 2, 0}},  // the hole's rims: SmCircles
    {"incline", {0, 1, 1, 2}},   // SmCircle bottom rim, line seam, two SmEllipse arcs
};

size_t g_analyticTrims = 0;

// The round trip; returns the first export (caller frees) for the hostile cases.
//
// `normalizes`: the import rewrites this solid into an equivalent canonical form, so the first
// re-export differs from the first export. SmuConvert's import does this to cones with a negative
// semi-angle (SmApiCreateCone's narrowing cone): it reverses the axis and both STEP parameters
// and rebuilds the cone from Vmin. Then the volume (integrated at 1e-7 over a different
// parameterisation) need only agree to 1e-7, and the arrays must reach a fixed point: the second
// and third exports must match.
//
// `volumeTolerance`: 1e-9 unless an edge curve is rebuilt with another parameterisation, which
// SMLib then integrates (and trims) differently, by up to its 1e-7 accuracy: the cone (above),
// the sphere, whose seam SmCircle is rebuilt over [-90, 90] degrees instead of [0, 1] (its
// volume moves by 1.9e-8 at accuracy 1e-7 and 8.5e-8 at 1e-8, where an all-NURBS export keeps it
// exactly), and the inclined cut, whose SmEllipse arcs are rebuilt likewise and whose cylinder
// trims SMLib rebuilds. The geometry is the same: export checks the rebuilt curve against the
// original pointwise.
smlib_flat *RoundTrip(const char *name, smlib_brep *b, bool normalizes = false, double volumeTolerance = 1e-9)
{
    if (b == nullptr)
    {
        std::fprintf(stderr, "%s: no solid\n", name);
        failures++;
        return nullptr;
    }
    Facts f0 = Measure(b);

    smlib_flat *flat1 = nullptr;
    int32_t s = smlib_brep_export(b, &flat1);
    CHECK(s == SMLIB_OK);
    if (s != SMLIB_OK)
    {
        std::fprintf(stderr, "%s: export -> %d\n", name, s);
        return nullptr;
    }
    smlib_brep *b2 = nullptr;
    s = smlib_brep_import(flat1, &b2);
    CHECK(s == SMLIB_OK);
    if (s != SMLIB_OK)
    {
        std::fprintf(stderr, "%s: import -> %d\n", name, s);
        return flat1;
    }
    Facts f1 = Measure(b2);
    double v8a = 0.0, v8b = 0.0;
    smlib_brep_volume(b, 1e-8, &v8a);
    smlib_brep_volume(b2, 1e-8, &v8b);
    std::printf("%-9s volume at 1e-8: %.12g -> %.12g (rel %.2e)\n", name, v8a, v8b, std::fabs(v8b - v8a) / std::max(1e-300, std::fabs(v8a)));
    double dv = std::fabs(f1.volume - f0.volume) / std::fabs(f0.volume);
    CHECK(dv <= volumeTolerance);
    CHECK(f1.counts[0] == f0.counts[0] && f1.counts[1] == f0.counts[1] && f1.counts[2] == f0.counts[2]);
    CHECK(f0.manifold == 1 && f1.manifold == 1);

    smlib_flat *flat2 = nullptr;
    CHECK(smlib_brep_export(b2, &flat2) == SMLIB_OK);
    Arrays a1 = Snapshot(flat1), a2 = Snapshot(flat2);
    double worst = 0.0;
    if (normalizes)
    {
        smlib_brep *b3 = nullptr;
        smlib_flat *flat4 = nullptr;
        CHECK(smlib_brep_import(flat2, &b3) == SMLIB_OK);
        CHECK(b3 != nullptr && smlib_brep_export(b3, &flat4) == SMLIB_OK);
        if (flat4 != nullptr)
        {
            Facts f2 = Measure(b3);
            CHECK(std::fabs(f2.volume - f1.volume) <= 1e-9 * std::fabs(f1.volume));
            worst = Compare(a2, Snapshot(flat4), name);
        }
        smlib_flat_free(flat4);
        smlib_brep_free(b3);
    }
    else
    {
        worst = Compare(a1, a2, name);
    }
    CHECK(worst <= 1e-9);

    // The text form reads back to the same arrays, exactly.
    std::string path = dir + "/" + name + ".flat.txt";
    CHECK(smlib_flat_write(flat1, path.c_str()) == SMLIB_OK);
    smlib_flat *flat3 = nullptr;
    CHECK(smlib_flat_read(path.c_str(), &flat3) == SMLIB_OK);
    if (flat3 != nullptr)
    {
        Arrays a3 = Snapshot(flat3);
        CHECK(a3.i == a1.i && a3.f == a1.f);
    }

    // Curve kinds of the first export: NURBS, line, circle, ellipse (edges and wire edges).
    int64_t kinds[4] = {0, 0, 0, 0};
    for (int32_t id : {SMLIB_FLAT_EDGE_CURVE_TYPE, SMLIB_FLAT_WIRE_EDGE_CURVE_TYPE})
    {
        for (int64_t t : a1.i[id])
        {
            CHECK(t >= 0 && t < 4);
            if (t >= 0 && t < 4)
            {
                kinds[t]++;
            }
        }
    }
    auto expected = g_kinds.find(name);
    if (expected != g_kinds.end())
    {
        const std::vector<int64_t> &k = expected->second;
        if (!(kinds[0] == k[0] && kinds[1] == k[1] && kinds[2] == k[2] && kinds[3] == k[3]))
        {
            std::fprintf(stderr, "%s: curve kinds nurbs/line/circle/ellipse %lld/%lld/%lld/%lld, expected %lld/%lld/%lld/%lld\n", name,
                         (long long)kinds[0], (long long)kinds[1], (long long)kinds[2], (long long)kinds[3], (long long)k[0], (long long)k[1],
                         (long long)k[2], (long long)k[3]);
            failures++;
        }
    }

    // UV trims carried on edges exported as analytic curves (their knots were re-mapped).
    size_t nAnalyticTrims = 0;
    if (!a1.i[SMLIB_FLAT_EDGEUSE_NURBS_VERTEX_COUNT].empty())
    {
        for (size_t k = 0; k < a1.i[SMLIB_FLAT_EDGEUSE_EDGE].size(); ++k)
        {
            nAnalyticTrims += a1.i[SMLIB_FLAT_EDGEUSE_NURBS_VERTEX_COUNT][k] != 0 &&
                              a1.i[SMLIB_FLAT_EDGE_CURVE_TYPE][static_cast<size_t>(a1.i[SMLIB_FLAT_EDGEUSE_EDGE][k])] != SMLIB_FLAT_CURVE_NURBS;
        }
    }
    g_analyticTrims += nAnalyticTrims;

    size_t nNurbsFaces = a1.i[SMLIB_FLAT_FACE_NURBS_U_ORDER].size();
    size_t nUV = 0;
    for (int64_t n : a1.i[SMLIB_FLAT_EDGEUSE_NURBS_VERTEX_COUNT])
    {
        nUV += n != 0;
    }
    std::printf(
        "%-9s faces/edges/vertices %lld/%lld/%lld -> %lld/%lld/%lld, manifold %d -> %d, volume %.12g -> %.12g (rel %.2e), "
        "NURBS faces %zu, UV trims %zu (%zu on analytic edges), curves nurbs/line/circle/ellipse %lld/%lld/%lld/%lld, re-export max rel diff %.2e\n",
        name, (long long)f0.counts[0], (long long)f0.counts[1], (long long)f0.counts[2], (long long)f1.counts[0], (long long)f1.counts[1],
        (long long)f1.counts[2], f0.manifold, f1.manifold, f0.volume, f1.volume, dv, nNurbsFaces, nUV, nAnalyticTrims, (long long)kinds[0], (long long)kinds[1],
        (long long)kinds[2], (long long)kinds[3], worst
    );

    smlib_flat_free(flat3);
    smlib_flat_free(flat2);
    smlib_brep_free(b2);
    smlib_brep_free(b);
    return flat1;
}

// A region's shells may come in any order: import puts the one that encloses the others first,
// as SMLib requires (solidcore#388). An island inside a void: a cube of 50 hollowed by a cube of
// 40, joined with a cube of 20 inside the cavity, which a cube of 10 hollows in turn. Each order
// of the regions (the infinite one first) and of the shells within each region must import as a
// valid body of the same volume.
void NestedShellOrder()
{
    auto box = [](double lo, double size) {
        const double o[3] = {lo, lo, lo}, s[3] = {size, size, size};
        smlib_brep *b = nullptr;
        CHECK(smlib_make_box(o, s, &b) == SMLIB_OK);
        return b;
    };
    smlib_brep *b50 = box(-25, 50), *b40 = box(-20, 40), *b20 = box(-10, 20), *b10 = box(-5, 10);
    smlib_brep *hollow = nullptr, *joined = nullptr, *nested = nullptr;
    CHECK(smlib_boolean(b50, b40, SMLIB_BOOLEAN_SUBTRACT, &hollow) == SMLIB_OK);
    CHECK(hollow && smlib_boolean(hollow, b20, SMLIB_BOOLEAN_UNION, &joined) == SMLIB_OK);
    CHECK(joined && smlib_boolean(joined, b10, SMLIB_BOOLEAN_SUBTRACT, &nested) == SMLIB_OK);
    smlib_flat *flat = nullptr;
    CHECK(nested && smlib_brep_export(nested, &flat) == SMLIB_OK);
    for (smlib_brep *b : {b50, b40, b20, b10, hollow, joined, nested})
    {
        smlib_brep_free(b);
    }
    if (flat == nullptr)
    {
        return;
    }
    const double want = 50.0 * 50 * 50 - 40.0 * 40 * 40 + 20.0 * 20 * 20 - 10.0 * 10 * 10;
    const Arrays a = Snapshot(flat);
    smlib_flat_free(flat);

    // Each shell's block: its counts and its faceuses. Faceuse shells only here.
    struct Shell
    {
        int64_t faceuses, wires, point;
        std::vector<int64_t> face, orientation;
    };
    const std::vector<int64_t> &regionShells = a.i.at(SMLIB_FLAT_REGION_SHELL_COUNT);
    std::vector<std::vector<Shell>> regions;
    size_t shell = 0, fu = 0;
    for (int64_t n : regionShells)
    {
        regions.emplace_back();
        for (int64_t k = 0; k < n; ++k, ++shell)
        {
            Shell s;
            s.faceuses = a.i.at(SMLIB_FLAT_SHELL_FACEUSE_COUNT)[shell];
            s.wires = a.i.at(SMLIB_FLAT_SHELL_WIRE_EDGE_COUNT)[shell];
            s.point = a.i.at(SMLIB_FLAT_SHELL_POINT_TYPE)[shell];
            CHECK(s.faceuses > 0 && s.wires == 0 && s.point == 0);
            const std::vector<int64_t> &face = a.i.at(SMLIB_FLAT_FACEUSE_FACE), &orient = a.i.at(SMLIB_FLAT_FACEUSE_ORIENTATION);
            s.face.assign(face.begin() + static_cast<ptrdiff_t>(fu), face.begin() + static_cast<ptrdiff_t>(fu + s.faceuses));
            s.orientation.assign(orient.begin() + static_cast<ptrdiff_t>(fu), orient.begin() + static_cast<ptrdiff_t>(fu + s.faceuses));
            fu += static_cast<size_t>(s.faceuses);
            regions.back().push_back(s);
        }
    }
    // The infinite region with the 50 cube's outside; the 50 hollow; the void between 40 and 20;
    // the 20 hollow; the void inside 10.
    CHECK(regions.size() == 5);
    size_t multi = 0;
    for (auto &r : regions)
    {
        multi += r.size() > 1;
    }
    CHECK(multi == 3);

    // Every order: regions 1.. permuted, then the shells of every region.
    std::vector<size_t> regionOrder;
    for (size_t r = 1; r < regions.size(); ++r)
    {
        regionOrder.push_back(r);
    }
    int orders = 0, good = 0;
    do
    {
        std::vector<std::vector<size_t>> shellOrder(regions.size());
        for (size_t r = 0; r < regions.size(); ++r)
        {
            for (size_t k = 0; k < regions[r].size(); ++k)
            {
                shellOrder[r].push_back(k);
            }
        }
        for (;;)
        {
            Arrays p = a;
            std::vector<int64_t> &rsc = p.i[SMLIB_FLAT_REGION_SHELL_COUNT], &rt = p.i[SMLIB_FLAT_REGION_TYPE];
            std::vector<int64_t> &sfc = p.i[SMLIB_FLAT_SHELL_FACEUSE_COUNT], &swc = p.i[SMLIB_FLAT_SHELL_WIRE_EDGE_COUNT],
                                 &spt = p.i[SMLIB_FLAT_SHELL_POINT_TYPE], &ff = p.i[SMLIB_FLAT_FACEUSE_FACE],
                                 &fo = p.i[SMLIB_FLAT_FACEUSE_ORIENTATION];
            rsc.clear();
            rt.clear();
            sfc.clear();
            swc.clear();
            spt.clear();
            ff.clear();
            fo.clear();
            std::vector<size_t> order = {0};
            order.insert(order.end(), regionOrder.begin(), regionOrder.end());
            for (size_t r : order)
            {
                rsc.push_back(static_cast<int64_t>(regions[r].size()));
                rt.push_back(a.i.at(SMLIB_FLAT_REGION_TYPE)[r]);
                for (size_t k : shellOrder[r])
                {
                    const Shell &s = regions[r][k];
                    sfc.push_back(s.faceuses);
                    swc.push_back(s.wires);
                    spt.push_back(s.point);
                    ff.insert(ff.end(), s.face.begin(), s.face.end());
                    fo.insert(fo.end(), s.orientation.begin(), s.orientation.end());
                }
            }
            smlib_flat *f = FromArrays(p);
            smlib_brep *b = nullptr;
            int32_t s = smlib_brep_import(f, &b);
            double v = 0.0;
            int64_t len = 0, bad = -1;
            char buf[4096] = {0};
            if (s == SMLIB_OK)
            {
                CHECK(smlib_brep_volume(b, 1e-9, &v) == SMLIB_OK);
                CHECK(smlib_brep_check(b, 2, buf, sizeof buf, &len, &bad) == SMLIB_OK);
            }
            orders++;
            if (s == SMLIB_OK && bad == 0 && std::fabs(v - want) <= 1e-9 * want)
            {
                good++;
            }
            else
            {
                std::fprintf(stderr, "nested shells, order %d: import %d, volume %.12g (want %.12g), failed checks %lld %s\n", orders,
                             static_cast<int>(s), v, want, static_cast<long long>(bad), buf);
                failures++;
            }
            smlib_brep_free(b);
            smlib_flat_free(f);

            // the next shell order, odometer-wise over the regions
            size_t r = 0;
            while (r < shellOrder.size() && !std::next_permutation(shellOrder[r].begin(), shellOrder[r].end()))
            {
                r++;
            }
            if (r == shellOrder.size())
            {
                break;
            }
        }
    } while (std::next_permutation(regionOrder.begin(), regionOrder.end()));
    std::printf("nested shells: %d of %d region and shell orders import as a valid body of volume %.12g\n", good, orders, want);
}

// ---- One-face sheets given without UV trims: the trims SMLib makes on export ----------------

// A NURBS curve in flat form: order, full knots, cartesian control points (`dim` coordinates
// each) and weights.
struct Spline
{
    int order = 0;
    std::vector<double> knots, cvs, weights;
};

// The point at `t` (clamped to the knot range) of a rational B-spline with `dim` coordinates per
// control point: de Boor's algorithm on the homogeneous points.
std::vector<double> SplinePoint(const Spline &c, int dim, double t)
{
    const size_t p = static_cast<size_t>(c.order - 1), n = c.weights.size(), d1 = static_cast<size_t>(dim) + 1;
    t = std::min(std::max(t, c.knots[p]), c.knots[n]);
    size_t k = p;
    while (k + 1 < n && t >= c.knots[k + 1])
    {
        k++;
    }
    std::vector<std::vector<double>> d(p + 1, std::vector<double>(d1));
    for (size_t j = 0; j <= p; ++j)
    {
        const size_t i = j + k - p;
        for (size_t a = 0; a + 1 < d1; ++a)
        {
            d[j][a] = c.cvs[i * (d1 - 1) + a] * c.weights[i];
        }
        d[j][d1 - 1] = c.weights[i];
    }
    for (size_t r = 1; r <= p; ++r)
    {
        for (size_t j = p; j >= r; --j)
        {
            const double lo = c.knots[j + k - p], hi = c.knots[j + 1 + k - r];
            const double alpha = hi > lo ? (t - lo) / (hi - lo) : 0.0;
            for (size_t a = 0; a < d1; ++a)
            {
                d[j][a] = (1.0 - alpha) * d[j - 1][a] + alpha * d[j][a];
            }
        }
    }
    std::vector<double> out(d1 - 1);
    for (size_t a = 0; a + 1 < d1; ++a)
    {
        out[a] = d[p][a] / d[p][d1 - 1];
    }
    return out;
}

// A NURBS surface in flat form: control points U-major (point (i, j) at i * nv + j).
struct Patch
{
    int nu = 0, nv = 0;
    Spline u, v; // orders and knots only
    std::vector<double> cvs, weights;
};

// S(u, v): each U row's homogeneous points evaluated in v, then that column in u.
std::vector<double> PatchPoint(const Patch &s, double u, double v)
{
    Spline column;
    column.order = s.u.order;
    column.knots = s.u.knots;
    for (int i = 0; i < s.nu; ++i)
    {
        Spline row;
        row.order = s.v.order;
        row.knots = s.v.knots;
        for (int j = 0; j < s.nv; ++j)
        {
            const size_t k = static_cast<size_t>(i * s.nv + j);
            const double w = s.weights[k];
            row.cvs.insert(row.cvs.end(), {s.cvs[3 * k] * w, s.cvs[3 * k + 1] * w, s.cvs[3 * k + 2] * w, w});
            row.weights.push_back(1.0);
        }
        const std::vector<double> h = SplinePoint(row, 4, v);
        column.cvs.insert(column.cvs.end(), {h[0] / h[3], h[1] / h[3], h[2] / h[3]});
        column.weights.push_back(h[3]);
    }
    return SplinePoint(column, 3, u);
}

double Dist(const std::vector<double> &a, const std::vector<double> &b)
{
    return std::sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2]));
}

// One edge of a sheet: its curve, range, end vertices and tolerance (in flat terms: twice its
// zone tolerance).
struct SheetEdge
{
    Spline curve;
    double t0 = 0.0, t1 = 1.0;
    int64_t from = 0, to = 0;
    double tolerance = 1e-5;
};

// The flat arrays of a one-face open sheet: `patch` over its whole domain, bounded by one loop of
// `edges`, each used along its curve and given without a UV trim (SMLib makes them), written as
// SMLib's own (SOURCE 1: no healing): how SolidCraft hands SMLib a face whose pcurves it needs.
Arrays SheetArrays(const Patch &patch, const std::vector<SheetEdge> &edges, const std::vector<double> &vertices,
                   const std::vector<double> &vertexTolerances)
{
    Arrays a;
    const int64_t n = static_cast<int64_t>(edges.size());
    a.i[SMLIB_FLAT_REGION_SHELL_COUNT] = {1};
    a.i[SMLIB_FLAT_REGION_TYPE] = {1};
    a.i[SMLIB_FLAT_SHELL_FACEUSE_COUNT] = {2};
    a.i[SMLIB_FLAT_SHELL_WIRE_EDGE_COUNT] = {0};
    a.i[SMLIB_FLAT_SHELL_POINT_TYPE] = {0};
    a.i[SMLIB_FLAT_FACEUSE_FACE] = {0, 0};
    a.i[SMLIB_FLAT_FACEUSE_ORIENTATION] = {1, 0};
    a.i[SMLIB_FLAT_FACE_LOOP_COUNT] = {1};
    a.i[SMLIB_FLAT_FACE_SURFACE_TYPE] = {SMLIB_FLAT_SURFACE_NURBS};
    a.i[SMLIB_FLAT_FACE_TRIM_TYPE] = {0};
    a.i[SMLIB_FLAT_LOOP_EDGEUSE_COUNT] = {n};
    a.i[SMLIB_FLAT_LOOP_VERTEX] = {-1};
    a.i[SMLIB_FLAT_SOURCE] = {1};
    a.i[SMLIB_FLAT_FACE_NURBS_U_VERTEX_COUNT] = {patch.nu};
    a.i[SMLIB_FLAT_FACE_NURBS_V_VERTEX_COUNT] = {patch.nv};
    a.i[SMLIB_FLAT_FACE_NURBS_U_ORDER] = {patch.u.order};
    a.i[SMLIB_FLAT_FACE_NURBS_V_ORDER] = {patch.v.order};
    a.f[SMLIB_FLAT_FACE_NURBS_CVS] = patch.cvs;
    a.f[SMLIB_FLAT_FACE_NURBS_U_KNOTS] = patch.u.knots;
    a.f[SMLIB_FLAT_FACE_NURBS_V_KNOTS] = patch.v.knots;
    a.f[SMLIB_FLAT_FACE_NURBS_WEIGHTS] = patch.weights;
    a.f[SMLIB_FLAT_FACE_RANGE] = {patch.u.knots.front(), patch.v.knots.front(), patch.u.knots.back(), patch.v.knots.back()};
    double tol = 0.0;
    for (int64_t k = 0; k < n; ++k)
    {
        const SheetEdge &e = edges[static_cast<size_t>(k)];
        a.i[SMLIB_FLAT_EDGEUSE_EDGE].push_back(k);
        a.i[SMLIB_FLAT_EDGEUSE_ORIENTATION].push_back(1);
        a.i[SMLIB_FLAT_EDGEUSE_NEXT_RADIAL].push_back(k);
        a.i[SMLIB_FLAT_EDGEUSE_RADIAL_ENTRY].push_back(0);
        a.i[SMLIB_FLAT_EDGEUSE_NURBS_VERTEX_COUNT].push_back(0);
        a.i[SMLIB_FLAT_EDGEUSE_NURBS_ORDER].push_back(0);
        a.i[SMLIB_FLAT_EDGE_CURVE_TYPE].push_back(SMLIB_FLAT_CURVE_NURBS);
        a.i[SMLIB_FLAT_EDGE_VERTICES].insert(a.i[SMLIB_FLAT_EDGE_VERTICES].end(), {e.from, e.to});
        a.i[SMLIB_FLAT_EDGE_NURBS_VERTEX_COUNT].push_back(static_cast<int64_t>(e.curve.weights.size()));
        a.i[SMLIB_FLAT_EDGE_NURBS_ORDER].push_back(e.curve.order);
        std::vector<double> &cvs = a.f[SMLIB_FLAT_EDGE_NURBS_CVS], &knots = a.f[SMLIB_FLAT_EDGE_NURBS_KNOTS],
                            &weights = a.f[SMLIB_FLAT_EDGE_NURBS_WEIGHTS];
        cvs.insert(cvs.end(), e.curve.cvs.begin(), e.curve.cvs.end());
        knots.insert(knots.end(), e.curve.knots.begin(), e.curve.knots.end());
        weights.insert(weights.end(), e.curve.weights.begin(), e.curve.weights.end());
        a.f[SMLIB_FLAT_EDGE_RANGE].insert(a.f[SMLIB_FLAT_EDGE_RANGE].end(), {e.t0, e.t1});
        a.f[SMLIB_FLAT_EDGE_TOLERANCE].push_back(e.tolerance);
        tol = std::max(tol, e.tolerance);
    }
    for (size_t k = 0; k < vertices.size() / 3; ++k)
    {
        a.i[SMLIB_FLAT_VERTEX_POINT_TYPE].push_back(1);
        tol = std::max(tol, vertexTolerances[k]);
    }
    a.f[SMLIB_FLAT_VERTEX_POSITION] = vertices;
    a.f[SMLIB_FLAT_VERTEX_TOLERANCE] = vertexTolerances;
    a.f[SMLIB_FLAT_TOLERANCE] = {tol};
    return a;
}

// How an exported use's UV trim lifts onto its face: the largest distance (over 129 samples)
// between S(trim) and the edge's point at the same place along the edge; the edge's range and its
// zone tolerance (half its EDGE_TOLERANCE), as exported.
struct TrimFit
{
    double gap = INFINITY, t0 = 0.0, t1 = 0.0, tolerance = 0.0;
};

// A block of `n` elements of `per` values from `at` (in elements) of `v`.
std::vector<double> Block(const std::vector<double> &v, size_t at, size_t n, size_t per = 1)
{
    return std::vector<double>(v.begin() + static_cast<ptrdiff_t>(at * per), v.begin() + static_cast<ptrdiff_t>((at + n) * per));
}

std::vector<TrimFit> TrimFits(const Arrays &x)
{
    Patch s;
    s.nu = static_cast<int>(x.i.at(SMLIB_FLAT_FACE_NURBS_U_VERTEX_COUNT).at(0));
    s.nv = static_cast<int>(x.i.at(SMLIB_FLAT_FACE_NURBS_V_VERTEX_COUNT).at(0));
    s.u.order = static_cast<int>(x.i.at(SMLIB_FLAT_FACE_NURBS_U_ORDER).at(0));
    s.v.order = static_cast<int>(x.i.at(SMLIB_FLAT_FACE_NURBS_V_ORDER).at(0));
    s.u.knots = x.f.at(SMLIB_FLAT_FACE_NURBS_U_KNOTS);
    s.v.knots = x.f.at(SMLIB_FLAT_FACE_NURBS_V_KNOTS);
    s.cvs = x.f.at(SMLIB_FLAT_FACE_NURBS_CVS);
    s.weights = x.f.at(SMLIB_FLAT_FACE_NURBS_WEIGHTS);
    // every edge's curve (these sheets have NURBS edges only)
    std::vector<Spline> curves;
    {
        size_t cv = 0, kn = 0;
        const std::vector<int64_t> &count = x.i.at(SMLIB_FLAT_EDGE_NURBS_VERTEX_COUNT), &order = x.i.at(SMLIB_FLAT_EDGE_NURBS_ORDER);
        for (size_t e = 0; e < count.size(); ++e)
        {
            Spline c;
            c.order = static_cast<int>(order[e]);
            const size_t m = static_cast<size_t>(count[e]);
            c.cvs = Block(x.f.at(SMLIB_FLAT_EDGE_NURBS_CVS), cv, m, 3);
            c.weights = Block(x.f.at(SMLIB_FLAT_EDGE_NURBS_WEIGHTS), cv, m);
            c.knots = Block(x.f.at(SMLIB_FLAT_EDGE_NURBS_KNOTS), kn, m + static_cast<size_t>(c.order));
            cv += m;
            kn += m + static_cast<size_t>(c.order);
            curves.push_back(c);
        }
    }
    std::vector<TrimFit> out;
    const std::vector<int64_t> &useEdge = x.i.at(SMLIB_FLAT_EDGEUSE_EDGE), &same = x.i.at(SMLIB_FLAT_EDGEUSE_ORIENTATION);
    const std::vector<int64_t> &tc = x.i.at(SMLIB_FLAT_EDGEUSE_NURBS_VERTEX_COUNT), &to = x.i.at(SMLIB_FLAT_EDGEUSE_NURBS_ORDER);
    size_t cv = 0, kn = 0;
    for (size_t k = 0; k < useEdge.size(); ++k)
    {
        TrimFit fit;
        const size_t e = static_cast<size_t>(useEdge[k]);
        fit.t0 = x.f.at(SMLIB_FLAT_EDGE_RANGE).at(2 * e);
        fit.t1 = x.f.at(SMLIB_FLAT_EDGE_RANGE).at(2 * e + 1);
        fit.tolerance = 0.5 * x.f.at(SMLIB_FLAT_EDGE_TOLERANCE).at(e);
        const size_t m = static_cast<size_t>(tc.at(k));
        if (m == 0 || e >= curves.size())
        {
            out.push_back(fit);
            continue;
        }
        Spline trim;
        trim.order = static_cast<int>(to[k]);
        trim.cvs = Block(x.f.at(SMLIB_FLAT_EDGEUSE_NURBS_CVS), cv, m, 2);
        trim.weights = Block(x.f.at(SMLIB_FLAT_EDGEUSE_NURBS_WEIGHTS), cv, m);
        trim.knots = Block(x.f.at(SMLIB_FLAT_EDGEUSE_NURBS_KNOTS), kn, m + static_cast<size_t>(trim.order));
        cv += m;
        kn += m + static_cast<size_t>(trim.order);
        const double s0 = trim.knots[static_cast<size_t>(trim.order - 1)], s1 = trim.knots[m];
        fit.gap = 0.0;
        for (int q = 0; q <= 128; ++q)
        {
            const double f = q / 128.0;
            const std::vector<double> uv = SplinePoint(trim, 2, s0 + f * (s1 - s0));
            const double t = same[k] != 0 ? fit.t0 + f * (fit.t1 - fit.t0) : fit.t1 - f * (fit.t1 - fit.t0);
            fit.gap = std::max(fit.gap, Dist(PatchPoint(s, uv[0], uv[1]), SplinePoint(curves[e], 3, t)));
        }
        out.push_back(fit);
    }
    return out;
}

// Imports `a` and exports it (SMLib makes every trim it was not given) into `exported`.
bool ImportExport(const Arrays &a, const char *what, Arrays &exported)
{
    smlib_flat *f = FromArrays(a), *x = nullptr;
    smlib_brep *b = nullptr;
    const int32_t si = smlib_brep_import(f, &b);
    const int32_t se = si == SMLIB_OK ? smlib_brep_export(b, &x) : si;
    if (si != SMLIB_OK || se != SMLIB_OK)
    {
        std::fprintf(stderr, "%s: import %d, export %d\n", what, static_cast<int>(si), static_cast<int>(se));
        failures++;
    }
    else
    {
        exported = Snapshot(x);
    }
    smlib_flat_free(x);
    smlib_brep_free(b);
    smlib_flat_free(f);
    return si == SMLIB_OK && se == SMLIB_OK;
}

// The patch of the trim tests: 10 x 10 mm, rational bicubic, with a crease (a triple u knot) at
// u = 0.08. On its v = 0 row (weights 1 there) S(u, 0) = (10u, 0, max(0, u - 0.08)): two lines.
Patch TrimTestPatch()
{
    Patch s;
    s.nu = 7;
    s.nv = 4;
    s.u.order = s.v.order = 4;
    s.u.knots = {0, 0, 0, 0, 0.08, 0.08, 0.08, 1, 1, 1, 1};
    s.v.knots = {0, 0, 0, 0, 1, 1, 1, 1};
    // Greville abscissae, so that x = 10u, y = 10v and the row's z are linear in u on each side
    const double gv[4] = {0.0, 1.0 / 3.0, 2.0 / 3.0, 1.0};
    const double z[4] = {0.0, 0.0, 2.0 / 3.0, 2.0}, w[4] = {1.0, 0.8, 0.8, 1.0};
    for (size_t i = 0; i < 7; ++i)
    {
        const double gu = (s.u.knots[i + 1] + s.u.knots[i + 2] + s.u.knots[i + 3]) / 3.0;
        for (int j = 0; j < 4; ++j)
        {
            s.cvs.insert(s.cvs.end(), {10.0 * gu, 10.0 * gv[j], std::max(0.0, gu - 0.08) + z[j]});
            s.weights.push_back(w[j]);
        }
    }
    return s;
}

// The patch's row `index` (v fixed) or column `index` (u fixed) as a curve, run backwards when
// `reversed`.
Spline PatchIso(const Patch &s, bool row, int index, bool reversed)
{
    Spline c;
    c.order = row ? s.u.order : s.v.order;
    const int n = row ? s.nu : s.nv;
    const std::vector<double> &k = row ? s.u.knots : s.v.knots;
    for (int q = 0; q < n; ++q)
    {
        const int p = reversed ? n - 1 - q : q;
        const size_t at = static_cast<size_t>(row ? p * s.nv + index : index * s.nv + p);
        c.cvs.insert(c.cvs.end(), {s.cvs[3 * at], s.cvs[3 * at + 1], s.cvs[3 * at + 2]});
        c.weights.push_back(s.weights[at]);
    }
    for (size_t q = 0; q < k.size(); ++q)
    {
        c.knots.push_back(reversed ? k.front() + k.back() - k[k.size() - 1 - q] : k[q]);
    }
    return c;
}

// The patch's u = 1 column, v = 1 row and u = 0 column, as the loop's edges 1 to 3 (vertex 1 to
// vertex 2, 2 to 3, 3 to 0).
void PatchSides(const Patch &s, std::vector<SheetEdge> &edges)
{
    SheetEdge e1, e2, e3;
    e1.curve = PatchIso(s, false, s.nu - 1, false);
    e1.from = 1;
    e1.to = 2;
    e2.curve = PatchIso(s, true, s.nv - 1, true);
    e2.from = 2;
    e2.to = 3;
    e3.curve = PatchIso(s, false, 0, true);
    e3.from = 3;
    e3.to = 0;
    edges.insert(edges.end(), {e1, e2, e3});
}

// solidcore#1383: an edge along the patch's v = 0 row, creased with it at u = 0.08 (a triple knot),
// parameterised unlike the surface's u (its control points at s^0.8 of the way along each line
// for their Greville values s), that starts 0.05 mm past the patch's u = 0 side, on the row's
// extension, its vertex taking no more tolerance than the file gave it. At the edge's tolerance
// SMLib's drop does not take that start, and comes back as the pieces on the patch, split at the
// crease; the longer one covers 96% of the edge and was stretched onto the whole edge: the trim
// ran at the wrong pace, 0.8 mm off the edge at its start. The trim lifts onto the edge everywhere
// within the edge's tolerance (as SMLib sizes it); the edge stays whole.
void BoundaryRowTrimPace()
{
    const Patch s = TrimTestPatch();
    SheetEdge e0;
    e0.curve.order = 4;
    e0.curve.knots = {0, 0, 0, 0, 0.04, 0.04, 0.04, 0.3, 0.6, 1, 1, 1, 1};
    const double start[3] = {-0.05, 0, 0}, crease[3] = {0.8, 0, 0}, end[3] = {10, 0, 0.92};
    for (size_t i = 0; i + 4 < e0.curve.knots.size(); ++i)
    {
        const double g = (e0.curve.knots[i + 1] + e0.curve.knots[i + 2] + e0.curve.knots[i + 3]) / 3.0;
        const bool first = g <= 0.04;
        const double *a = first ? start : crease, *b = first ? crease : end;
        const double f = std::pow(first ? g / 0.04 : (g - 0.04) / 0.96, 0.8);
        e0.curve.cvs.insert(e0.curve.cvs.end(), {a[0] + f * (b[0] - a[0]), a[1] + f * (b[1] - a[1]), a[2] + f * (b[2] - a[2])});
        e0.curve.weights.push_back(1.0);
    }
    e0.from = 0;
    e0.to = 1;
    e0.tolerance = 1e-4;
    std::vector<SheetEdge> edges = {e0};
    PatchSides(s, edges);
    const std::vector<double> vertices = {-0.05, 0, 0, 10, 0, 0.92, 10, 10, 2.92, 0, 10, 2};
    Arrays x;
    if (!ImportExport(SheetArrays(s, edges, vertices, {1e-5, 1e-5, 1e-5, 1e-5}), "boundary row trim pace", x))
    {
        return;
    }
    const std::vector<TrimFit> fits = TrimFits(x);
    CHECK(fits.size() == 4);
    for (size_t k = 0; k < fits.size(); ++k)
    {
        const TrimFit &f = fits[k];
        std::printf("boundary row trim pace: edge %zu range [%.9g, %.9g], its trim %.3g mm off it, tolerance %.3g mm\n", k, f.t0, f.t1, f.gap,
                    f.tolerance);
        CHECK(f.t0 == 0.0 && f.t1 == 1.0);
        CHECK(f.gap <= f.tolerance);
    }
}

// An edge along the patch's v = 0 row whose curve runs 2.5 um past its end vertex (the patch's
// corner) and the patch's u = 1 side, on the row's extension, with a tolerance of 0.5 um (as
// Poppy's Congé19 has them). Making its trim, SMLib trims the edge back to the patch, to its
// vertex: the edge does not lie within its tolerance of the face. smlib_brep_edge_face_distance
// measures the 2.5 um; given that as its tolerance, the edge keeps its interval, its trim within it.
void EdgePastPatch()
{
    const Patch s = TrimTestPatch();
    const double past = 2.5e-3; // mm
    SheetEdge e0;
    e0.curve.order = 4;
    e0.curve.knots = {0, 0, 0, 0, 0.04, 0.04, 0.04, 0.3, 0.6, 1, 1, 1, 1};
    // the row's two lines, the second one carried on past u = 1
    const double start[3] = {0, 0, 0}, crease[3] = {0.8, 0, 0}, end[3] = {10 + past, 0, 0.92 + 0.1 * past};
    for (size_t i = 0; i + 4 < e0.curve.knots.size(); ++i)
    {
        const double g = (e0.curve.knots[i + 1] + e0.curve.knots[i + 2] + e0.curve.knots[i + 3]) / 3.0;
        const bool first = g <= 0.04;
        const double *a = first ? start : crease, *b = first ? crease : end;
        const double f = first ? g / 0.04 : (g - 0.04) / 0.96;
        e0.curve.cvs.insert(e0.curve.cvs.end(), {a[0] + f * (b[0] - a[0]), a[1] + f * (b[1] - a[1]), a[2] + f * (b[2] - a[2])});
        e0.curve.weights.push_back(1.0);
    }
    e0.from = 0;
    e0.to = 1;
    e0.tolerance = 1e-3;
    std::vector<SheetEdge> edges = {e0};
    PatchSides(s, edges);
    const std::vector<double> vertices = {0, 0, 0, 10, 0, 0.92, 10, 10, 2.92, 0, 10, 2};
    const Arrays given = SheetArrays(s, edges, vertices, {1e-5, 1e-5, 1e-5, 1e-5});

    // how far each edge lies from the face's surface
    double d[4] = {-1.0, -1.0, -1.0, -1.0};
    {
        smlib_flat *f = FromArrays(given);
        smlib_brep *b = nullptr;
        CHECK(smlib_brep_import(f, &b) == SMLIB_OK);
        for (int64_t k = 0; b != nullptr && k < 4; ++k)
        {
            CHECK(smlib_brep_edge_face_distance(b, k, 0, &d[k]) == SMLIB_OK);
        }
        double x = 0.0;
        CHECK(smlib_brep_edge_face_distance(b, 4, 0, &x) == SMLIB_ERR_ARG);
        CHECK(smlib_brep_edge_face_distance(b, 0, 1, &x) == SMLIB_ERR_ARG);
        CHECK(smlib_brep_edge_face_distance(b, -1, 0, &x) == SMLIB_ERR_ARG);
        CHECK(smlib_brep_edge_face_distance(b, 0, 0, nullptr) == SMLIB_ERR_ARG);
        CHECK(smlib_brep_edge_face_distance(nullptr, 0, 0, &x) == SMLIB_ERR_ARG);
        smlib_brep_free(b);
        smlib_flat_free(f);
    }
    std::printf("edge past its patch: edges %.3g, %.3g, %.3g, %.3g mm from the face's surface\n", d[0], d[1], d[2], d[3]);
    // the edge's end, 2.5 um along the row's extension: sqrt(1 + 0.1^2) * 2.5 um from the patch's corner
    CHECK(std::fabs(d[0] - past * std::sqrt(1.01)) <= 1e-9);
    CHECK(d[1] >= 0.0 && d[1] <= 1e-9 && d[2] >= 0.0 && d[2] <= 1e-9 && d[3] >= 0.0 && d[3] <= 1e-9);

    Arrays x;
    if (ImportExport(given, "edge past its patch", x))
    {
        const TrimFit f = TrimFits(x).at(0);
        std::printf("edge past its patch, at 0.5 um: range [%.9g, %.9g], its trim %.3g mm off it, tolerance %.3g mm\n", f.t0, f.t1, f.gap, f.tolerance);
        // trimmed back to the patch, its end to its vertex
        CHECK(f.t0 == 0.0 && f.t1 < 1.0 && f.t1 > 0.999);
    }
    Arrays measured = given;
    measured.f[SMLIB_FLAT_EDGE_TOLERANCE][0] = 2.0 * d[0] * (1.0 + 1e-6);
    measured.f[SMLIB_FLAT_TOLERANCE][0] = std::max(measured.f[SMLIB_FLAT_TOLERANCE][0], measured.f[SMLIB_FLAT_EDGE_TOLERANCE][0]);
    if (ImportExport(measured, "edge past its patch, measured", x))
    {
        const TrimFit f = TrimFits(x).at(0);
        std::printf("edge past its patch, at its distance: range [%.9g, %.9g], its trim %.3g mm off it, tolerance %.3g mm\n", f.t0, f.t1, f.gap,
                    f.tolerance);
        CHECK(f.t0 == 0.0 && f.t1 == 1.0);
        CHECK(f.gap <= f.tolerance);
    }
}

// ---- A torus face given as its label, then as SMLib's export of it --------------------------

// A face given as its label, imported and exported (SMLib gives it its NURBS form), then that
// export imported and exported twice more: the form comes back as it went in, to the bit; the
// label's reference direction (`refId`) is SMLib's, from the form (the direction to the NURBS's
// start about the axis `axisId`: IsNurbSurfaceTorus and IsNurbSurfaceSphere find the frame in the
// form), a unit vector square to the same axis; everything else is the same, and so is the next
// export.
void FormRoundTrip(const Arrays &given, const char *kind, const char *name, int32_t axisId, int32_t refId)
{
    Arrays first, second, third;
    if (!(ImportExport(given, name, first) && ImportExport(first, name, second) && ImportExport(second, name, third)))
    {
        std::printf("%s %s: its export does not import again\n", kind, name);
        return;
    }
    bool sameForm = true;
    for (int32_t id : {SMLIB_FLAT_FACE_NURBS_CVS, SMLIB_FLAT_FACE_NURBS_U_KNOTS, SMLIB_FLAT_FACE_NURBS_V_KNOTS, SMLIB_FLAT_FACE_NURBS_WEIGHTS})
    {
        sameForm = sameForm && first.f.at(id) == second.f.at(id);
    }
    const std::vector<double> &axis = second.f.at(axisId), &xr = second.f.at(refId);
    const double square = std::fabs(axis[0] * xr[0] + axis[1] * xr[1] + axis[2] * xr[2]);
    const double unit = std::fabs(std::sqrt(xr[0] * xr[0] + xr[1] * xr[1] + xr[2] * xr[2]) - 1.0);
    Arrays firstInFrame = first;
    firstInFrame.f[refId] = xr;
    const double d = Compare(firstInFrame, second, name), fixed = Compare(second, third, name);
    std::printf("%s %s: its export imports again, form %s, then the same but the reference direction (worst %.3g), then a fixed point "
                "(worst %.3g)\n",
                kind, name, sameForm ? "kept" : "changed", d, fixed);
    CHECK(sameForm && square <= 1e-12 && unit <= 1e-12 && d <= 1e-9 && fixed <= 1e-12);
}

// A one-face torus sheet as a STEP reader hands it to SMLib: the torus label (frame o, z, x,
// radii R, r) over the STEP angles [u0, u1] x [v0, v1] (radians), no NURBS form, bounded by its
// four iso circles (circle edges, no UV trims). The torus is
// S(u, v) = o + (R + r cos v)(cos u x + sin u y) + r sin v z, with y = z cross x.
Arrays TorusSheet(const double o[3], const double z[3], const double x[3], double R, double r, double u0, double u1, double v0, double v1)
{
    const double y[3] = {z[1] * x[2] - z[2] * x[1], z[2] * x[0] - z[0] * x[2], z[0] * x[1] - z[1] * x[0]};
    const auto radial = [&](double u, int k) { return std::cos(u) * x[k] + std::sin(u) * y[k]; };
    const auto point = [&](double u, double v, int k) { return o[k] + (R + r * std::cos(v)) * radial(u, k) + r * std::sin(v) * z[k]; };
    Arrays a;
    a.i[SMLIB_FLAT_REGION_SHELL_COUNT] = {1};
    a.i[SMLIB_FLAT_REGION_TYPE] = {1};
    a.i[SMLIB_FLAT_SHELL_FACEUSE_COUNT] = {2};
    a.i[SMLIB_FLAT_SHELL_WIRE_EDGE_COUNT] = {0};
    a.i[SMLIB_FLAT_SHELL_POINT_TYPE] = {0};
    a.i[SMLIB_FLAT_FACEUSE_FACE] = {0, 0};
    a.i[SMLIB_FLAT_FACEUSE_ORIENTATION] = {1, 0};
    a.i[SMLIB_FLAT_FACE_LOOP_COUNT] = {1};
    a.i[SMLIB_FLAT_FACE_SURFACE_TYPE] = {SMLIB_FLAT_SURFACE_TORUS};
    a.i[SMLIB_FLAT_FACE_TRIM_TYPE] = {0};
    a.i[SMLIB_FLAT_LOOP_EDGEUSE_COUNT] = {4};
    a.i[SMLIB_FLAT_LOOP_VERTEX] = {-1};
    a.i[SMLIB_FLAT_SOURCE] = {1};
    a.f[SMLIB_FLAT_FACE_TORUS_ORIGIN] = {o[0], o[1], o[2]};
    a.f[SMLIB_FLAT_FACE_TORUS_AXIS] = {z[0], z[1], z[2]};
    a.f[SMLIB_FLAT_FACE_TORUS_REF_DIRECTION] = {x[0], x[1], x[2]};
    a.f[SMLIB_FLAT_FACE_TORUS_MAJOR_RADIUS] = {R};
    a.f[SMLIB_FLAT_FACE_TORUS_MINOR_RADIUS] = {r};
    a.f[SMLIB_FLAT_FACE_RANGE] = {u0, v0, u1, v1};
    // Vertices 0 to 3 at (u0, v0), (u1, v0), (u1, v1), (u0, v1). Edges: the major circles at v0
    // and v1 (u0 to u1), the minor circles at u1 and u0 (v0 to v1); the loop runs counter-
    // clockwise in (u, v), along edges 0 and 1, against 2 and 3.
    const double corner[4][2] = {{u0, v0}, {u1, v0}, {u1, v1}, {u0, v1}};
    for (const auto &c : corner)
    {
        for (int k = 0; k < 3; ++k)
        {
            a.f[SMLIB_FLAT_VERTEX_POSITION].push_back(point(c[0], c[1], k));
        }
        a.i[SMLIB_FLAT_VERTEX_POINT_TYPE].push_back(1);
        a.f[SMLIB_FLAT_VERTEX_TOLERANCE].push_back(1e-5);
    }
    const auto majorCircle = [&](double v)
    {
        for (int k = 0; k < 3; ++k)
        {
            a.f[SMLIB_FLAT_EDGE_CIRCLE_CENTER].push_back(o[k] + r * std::sin(v) * z[k]);
            a.f[SMLIB_FLAT_EDGE_CIRCLE_AXIS].push_back(z[k]);
            a.f[SMLIB_FLAT_EDGE_CIRCLE_REF_DIRECTION].push_back(x[k]);
        }
        a.f[SMLIB_FLAT_EDGE_CIRCLE_RADIUS].push_back(R + r * std::cos(v));
        a.f[SMLIB_FLAT_EDGE_RANGE].insert(a.f[SMLIB_FLAT_EDGE_RANGE].end(), {u0, u1});
    };
    // The minor circle at u: centre o + R radial(u), reference direction radial(u), axis
    // radial(u) cross z, so that its y axis is z and its angle is v.
    const auto minorCircle = [&](double u)
    {
        const double t[3] = {radial(u, 0), radial(u, 1), radial(u, 2)};
        const double n[3] = {t[1] * z[2] - t[2] * z[1], t[2] * z[0] - t[0] * z[2], t[0] * z[1] - t[1] * z[0]};
        for (int k = 0; k < 3; ++k)
        {
            a.f[SMLIB_FLAT_EDGE_CIRCLE_CENTER].push_back(o[k] + R * t[k]);
            a.f[SMLIB_FLAT_EDGE_CIRCLE_AXIS].push_back(n[k]);
            a.f[SMLIB_FLAT_EDGE_CIRCLE_REF_DIRECTION].push_back(t[k]);
        }
        a.f[SMLIB_FLAT_EDGE_CIRCLE_RADIUS].push_back(r);
        a.f[SMLIB_FLAT_EDGE_RANGE].insert(a.f[SMLIB_FLAT_EDGE_RANGE].end(), {v0, v1});
    };
    majorCircle(v0);
    minorCircle(u1);
    majorCircle(v1);
    minorCircle(u0);
    a.i[SMLIB_FLAT_EDGE_VERTICES] = {0, 1, 1, 2, 3, 2, 0, 3};
    for (int64_t k = 0; k < 4; ++k)
    {
        a.i[SMLIB_FLAT_EDGE_CURVE_TYPE].push_back(SMLIB_FLAT_CURVE_CIRCLE);
        a.f[SMLIB_FLAT_EDGE_TOLERANCE].push_back(2e-5);
        a.i[SMLIB_FLAT_EDGEUSE_EDGE].push_back(k);
        a.i[SMLIB_FLAT_EDGEUSE_ORIENTATION].push_back(k < 2 ? 1 : 0);
        a.i[SMLIB_FLAT_EDGEUSE_NEXT_RADIAL].push_back(k);
        a.i[SMLIB_FLAT_EDGEUSE_RADIAL_ENTRY].push_back(0);
    }
    a.f[SMLIB_FLAT_TOLERANCE] = {4e-5};
    return a;
}

// A torus face SMLib made from its label (SmTorus::CreateCanonical over the face's angles) comes
// back from export as that label and its NURBS form; the form imported again must be kept, to the
// bit (the face's trims are in its parameters), and export -> import must be the identity but for
// the label's reference direction. The torus SMLib rebuilt from the circles it fits to the form
// was ~1e-11 off it on a small torus far out, and the import refused it. These are the faces a
// STEP file has, parts of a torus: a quarter, one off the seam and across the minor angle 0, the
// inner side, a half about the seam, the off-seam one about an axis turned over (its minor
// circles running the other way in space), and a small one far from the origin.
void TorusFormRoundTrip()
{
    const double pi = 3.14159265358979323846, deg = pi / 180.0;
    const double o[3] = {1, 2, 3}, up[3] = {0, 0, 1}, down[3] = {0, 0, -1}, ref[3] = {1, 0, 0};
    // Poppy's Stock-head has a fillet torus R 2.85, r 0.2 about the y axis, 90 mm out: SMLib's
    // torus rebuilt from the circles it fits to that NURBS had weights 7e-12 off it.
    const double far[3] = {0, 88.18, -14.75}, side[3] = {0, 1, 0}, across[3] = {-1, 0, 0};
    struct Case
    {
        const char *name;
        const double *o, *z, *x;
        double R, r, u0, u1, v0, v1;
    };
    const Case cases[] = {
        {"quarter", o, up, ref, 10.0, 3.0, 0.0, 90.0, 0.0, 90.0},
        {"off the seam", o, up, ref, 10.0, 3.0, 30.0, 100.0, -50.0, 70.0},
        {"inner side", o, up, ref, 10.0, 3.0, 17.0, 115.0, 100.0, 260.0},
        {"half about the seam", o, up, ref, 10.0, 3.0, -90.0, 90.0, -40.0, 40.0},
        {"axis turned over", o, down, ref, 10.0, 3.0, 30.0, 100.0, -50.0, 70.0},
        {"small, far out", far, side, across, 2.85, 0.2, 0.0, 90.0, 90.0, 180.0},
    };
    // The form SMLib exported for that Stock-head face (its control points are themselves 1e-12
    // off a plane square to the axis), given as SMLib writes it: the torus label (passed over on
    // import) and the form, bounded by the form's iso curves, no UV trims. Its first import was
    // refused; now it is kept to the bit.
    {
        Patch f;
        f.nu = f.nv = 3;
        f.u.order = f.v.order = 3;
        f.u.knots = {3, 3, 3, 4, 4, 4};
        f.v.knots = {4, 4, 4, 5, 5, 5};
        f.weights = {1, 0.70710678118654757, 1, 0.70710678118654813, 0.50000000000000044, 0.70710678118654813, 1, 0.70710678118654757, 1};
        f.cvs = {-2.8500000000000139,     91.036626133921814, -14.772670798218494, -2.8500000000000139,     91.024416426014824, -14.573043838534121,
                 -2.6500000000000137,     91.024416426014838, -14.573043838534121, -2.8500000000011951,     88.191941958423001, -14.946659135892627,
                 -2.85000000000118,       88.179732250516011, -14.747032176208254, -2.650000000001195,      88.379359210200406, -14.734822468301282,
                 -4.7201503979908935e-12, 88.191941958421822, -14.946659135892693, -4.7064028649897626e-12, 88.179732250514846, -14.74703217620832,
                 -4.7202196957604515e-12, 88.379359210199226, -14.734822468301347};
        std::vector<SheetEdge> edges(1);
        edges[0].curve = PatchIso(f, true, 0, false);
        edges[0].from = 0;
        edges[0].to = 1;
        PatchSides(f, edges);
        for (SheetEdge &e : edges)
        {
            e.t0 = e.curve.knots.front();
            e.t1 = e.curve.knots.back();
        }
        const std::vector<double> vertices = {f.cvs[0], f.cvs[1], f.cvs[2], f.cvs[18], f.cvs[19], f.cvs[20],
                                              f.cvs[24], f.cvs[25], f.cvs[26], f.cvs[6], f.cvs[7], f.cvs[8]};
        Arrays given = SheetArrays(f, edges, vertices, {1e-5, 1e-5, 1e-5, 1e-5});
        // The label: the sweep's centre and axis from its first row's control triangle.
        const double c[3] = {f.cvs[0] + f.cvs[18] - f.cvs[9], f.cvs[1] + f.cvs[19] - f.cvs[10], f.cvs[2] + f.cvs[20] - f.cvs[11]};
        const double p[3] = {f.cvs[0] - c[0], f.cvs[1] - c[1], f.cvs[2] - c[2]}, q[3] = {f.cvs[18] - c[0], f.cvs[19] - c[1], f.cvs[20] - c[2]};
        double n[3] = {p[1] * q[2] - p[2] * q[1], p[2] * q[0] - p[0] * q[2], p[0] * q[1] - p[1] * q[0]}, x[3] = {p[0], p[1], p[2]};
        const double nl = std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]), R = std::sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
        for (int k = 0; k < 3; ++k)
        {
            n[k] /= nl;
            x[k] /= R;
        }
        given.i[SMLIB_FLAT_FACE_SURFACE_TYPE] = {SMLIB_FLAT_SURFACE_TORUS};
        given.i[SMLIB_FLAT_FACE_NURBS_FORM] = {1};
        given.f[SMLIB_FLAT_FACE_TORUS_ORIGIN] = {c[0], c[1], c[2]};
        given.f[SMLIB_FLAT_FACE_TORUS_AXIS] = {n[0], n[1], n[2]};
        given.f[SMLIB_FLAT_FACE_TORUS_REF_DIRECTION] = {x[0], x[1], x[2]};
        given.f[SMLIB_FLAT_FACE_TORUS_MAJOR_RADIUS] = {R};
        given.f[SMLIB_FLAT_FACE_TORUS_MINOR_RADIUS] = {0.2};
        Arrays first, second;
        const bool ok = ImportExport(given, "Stock-head torus form", first) && ImportExport(first, "Stock-head torus form", second);
        bool kept = ok;
        for (int32_t id : {SMLIB_FLAT_FACE_NURBS_CVS, SMLIB_FLAT_FACE_NURBS_U_KNOTS, SMLIB_FLAT_FACE_NURBS_V_KNOTS, SMLIB_FLAT_FACE_NURBS_WEIGHTS})
        {
            kept = kept && given.f.at(id) == first.f.at(id) && first.f.at(id) == second.f.at(id);
        }
        const double fixed = ok ? Compare(first, second, "Stock-head torus form") : INFINITY;
        std::printf("torus face given as Stock-head's form: imports %s, form %s, then a fixed point (worst %.3g)\n", ok ? "ok" : "refused",
                    kept ? "kept" : "changed", fixed);
        CHECK(ok && kept && fixed <= 1e-12);
    }
    // A torus face given with a form that runs once round the tube from its bottom (tube angle
    // -90 degrees: SolidCraft forms a face covering most of its tube so), the face half of it, 80 mm
    // from the origin (an alternator's fillet, R 13.1934, r 1.25). Recognising the torus in the
    // form reads the tube's ends from the form's first and last points, each to the rounding of
    // its coordinates (8e-13 degrees here): the whole turn was taken for a sliver 2e-12 degrees
    // long, and the face refused.
    {
        const double R = 13.1934, r = 1.25, h = std::sqrt(0.5);
        const double c[3] = {-28.640080455029999, 39.419688939149999, -63.4255};
        const double z[3] = {-0.58778525229251344, 0.80901699437491859, 0}, x[3] = {0, 0, -1};
        // The sweep: a quarter turn from x to y = z cross x. The tube: a whole turn from -90
        // degrees, in four quarters. The control points are the alternator's, to the bit (the
        // ends' rounding decides): c + rho (x or y, or x + y in the middle) + height z, with rho
        // and height R, -r; R + r, -r; R + r, 0; ... about the tube's centre circle.
        Patch f;
        f.nu = 3;
        f.nv = 9;
        f.u.order = f.v.order = 3;
        f.u.knots = {3, 3, 3, 4, 4, 4};
        f.v.knots = {-1, -1, -1, 0, 0, 1, 1, 2, 2, 3, 3, 3};
        f.cvs = {-27.905348889664353, 38.408417696181353, -76.618899999999996, -27.905348889664353, 38.408417696181353, -77.868899999999996,
                 -28.640080455029995, 39.419688939149999, -77.868899999999996, -29.374812020395638, 40.430960182118646, -77.868899999999996,
                 -29.374812020395638, 40.430960182118646, -76.618899999999996, -29.374812020395641, 40.430960182118646, -75.368899999999996,
                 -28.640080455029999, 39.419688939149999, -75.368899999999996, -27.905348889664356, 38.408417696181353, -75.368899999999996,
                 -27.905348889664353, 38.408417696181353, -76.618899999999996, -38.579033703250403, 30.653531748585308, -76.618899999999996,
                 -39.590304946219049, 29.918800183219666, -77.868899999999996, -40.325036511584692, 30.930071426188313, -77.868899999999996,
                 -41.059768076950334, 31.94134266915696, -77.868899999999996, -40.048496833981687, 32.676074234522602, -76.618899999999996,
                 -39.037225591013041, 33.410805799888244, -75.368899999999996, -38.302494025647398, 32.399534556919598, -75.368899999999996,
                 -37.567762460281756, 31.388263313950951, -75.368899999999996, -38.579033703250403, 30.653531748585308, -76.618899999999996,
                 -38.579033703250403, 30.653531748585308, -63.4255, -39.590304946219057, 29.918800183219666, -63.4255,
                 -40.325036511584699, 30.930071426188313, -63.4255, -41.059768076950341, 31.941342669156956, -63.4255,
                 -40.048496833981687, 32.676074234522602, -63.4255, -39.037225591013041, 33.410805799888244, -63.4255,
                 -38.302494025647398, 32.399534556919598, -63.4255, -37.567762460281756, 31.388263313950951, -63.4255,
                 -38.579033703250403, 30.653531748585308, -63.4255};
        for (double ws : {1.0, h, 1.0})
        {
            for (int j = 0; j < 9; ++j)
            {
                f.weights.push_back(ws * (j % 2 == 0 ? 1.0 : h));
            }
        }
        // Bounded by the tube circles at the sweep's ends (v 0 to 2: the tube from its outer
        // equator over the top to its inner one) and the sweep arcs on those equators (control
        // rows 2 and 6).
        std::vector<SheetEdge> edges(4);
        edges[0].curve = PatchIso(f, true, 2, false);
        edges[0].from = 0;
        edges[0].to = 1;
        edges[1].curve = PatchIso(f, false, f.nu - 1, false);
        edges[1].from = 1;
        edges[1].to = 2;
        edges[2].curve = PatchIso(f, true, 6, true);
        edges[2].from = 2;
        edges[2].to = 3;
        edges[3].curve = PatchIso(f, false, 0, true);
        edges[3].from = 3;
        edges[3].to = 0;
        for (SheetEdge &e : edges)
        {
            const bool tubeCircle = e.curve.knots.size() == f.v.knots.size();
            e.t0 = tubeCircle ? 0.0 : e.curve.knots.front();
            e.t1 = tubeCircle ? 2.0 : e.curve.knots.back();
        }
        std::vector<double> vertices;
        for (const auto &uv : {std::array<double, 2>{3, 0}, {4, 0}, {4, 2}, {3, 2}})
        {
            const std::vector<double> p = PatchPoint(f, uv[0], uv[1]);
            vertices.insert(vertices.end(), p.begin(), p.end());
        }
        Arrays given = SheetArrays(f, edges, vertices, {1e-5, 1e-5, 1e-5, 1e-5});
        given.f[SMLIB_FLAT_FACE_RANGE] = {3, 0, 4, 2};
        given.i[SMLIB_FLAT_FACE_SURFACE_TYPE] = {SMLIB_FLAT_SURFACE_TORUS};
        given.i[SMLIB_FLAT_FACE_NURBS_FORM] = {1};
        given.f[SMLIB_FLAT_FACE_TORUS_ORIGIN] = {c[0], c[1], c[2]};
        given.f[SMLIB_FLAT_FACE_TORUS_AXIS] = {z[0], z[1], z[2]};
        given.f[SMLIB_FLAT_FACE_TORUS_REF_DIRECTION] = {x[0], x[1], x[2]};
        given.f[SMLIB_FLAT_FACE_TORUS_MAJOR_RADIUS] = {R};
        given.f[SMLIB_FLAT_FACE_TORUS_MINOR_RADIUS] = {r};
        Arrays first;
        const bool ok = ImportExport(given, "whole-tube torus form", first);
        bool kept = ok;
        for (int32_t id : {SMLIB_FLAT_FACE_NURBS_CVS, SMLIB_FLAT_FACE_NURBS_U_KNOTS, SMLIB_FLAT_FACE_NURBS_V_KNOTS, SMLIB_FLAT_FACE_NURBS_WEIGHTS})
        {
            kept = kept && given.f.at(id) == first.f.at(id);
        }
        const double dR = ok ? std::fabs(first.f.at(SMLIB_FLAT_FACE_TORUS_MAJOR_RADIUS).at(0) - R) : INFINITY;
        const double dr = ok ? std::fabs(first.f.at(SMLIB_FLAT_FACE_TORUS_MINOR_RADIUS).at(0) - r) : INFINITY;
        std::printf("torus face on a form once round its tube from the bottom, 80 mm out: imports %s, form %s, radii off by %.3g and %.3g\n",
                    ok ? "ok" : "refused", kept ? "kept" : "changed", dR, dr);
        CHECK(ok && kept && dR <= 1e-9 && dr <= 1e-9);
    }
    for (const Case &c : cases)
    {
        FormRoundTrip(TorusSheet(c.o, c.z, c.x, c.R, c.r, c.u0 * deg, c.u1 * deg, c.v0 * deg, c.v1 * deg), "torus face", c.name,
                      SMLIB_FLAT_FACE_TORUS_AXIS, SMLIB_FLAT_FACE_TORUS_REF_DIRECTION);
    }
}

// A one-face sphere sheet as a STEP reader hands it to SMLib: the sphere label (centre o, frame
// z, x, radius r) over the longitudes [u0, u1] (radians) from pole to pole, no NURBS form,
// bounded by its two meridians (circle edges from the south pole to the north one, no UV trims),
// its other two sides collapsed at the poles. S(u, v) = o + r cos v (cos u x + sin u y) + r sin v z.
Arrays SphereLune(const double o[3], const double z[3], const double x[3], double r, double u0, double u1)
{
    const double pi = 3.14159265358979323846;
    const double y[3] = {z[1] * x[2] - z[2] * x[1], z[2] * x[0] - z[0] * x[2], z[0] * x[1] - z[1] * x[0]};
    Arrays a;
    a.i[SMLIB_FLAT_REGION_SHELL_COUNT] = {1};
    a.i[SMLIB_FLAT_REGION_TYPE] = {1};
    a.i[SMLIB_FLAT_SHELL_FACEUSE_COUNT] = {2};
    a.i[SMLIB_FLAT_SHELL_WIRE_EDGE_COUNT] = {0};
    a.i[SMLIB_FLAT_SHELL_POINT_TYPE] = {0};
    a.i[SMLIB_FLAT_FACEUSE_FACE] = {0, 0};
    a.i[SMLIB_FLAT_FACEUSE_ORIENTATION] = {1, 0};
    a.i[SMLIB_FLAT_FACE_LOOP_COUNT] = {1};
    a.i[SMLIB_FLAT_FACE_SURFACE_TYPE] = {SMLIB_FLAT_SURFACE_SPHERE};
    a.i[SMLIB_FLAT_FACE_TRIM_TYPE] = {0};
    a.i[SMLIB_FLAT_LOOP_EDGEUSE_COUNT] = {2};
    a.i[SMLIB_FLAT_LOOP_VERTEX] = {-1};
    a.i[SMLIB_FLAT_SOURCE] = {1};
    a.f[SMLIB_FLAT_FACE_SPHERE_CENTER] = {o[0], o[1], o[2]};
    a.f[SMLIB_FLAT_FACE_SPHERE_AXIS] = {z[0], z[1], z[2]};
    a.f[SMLIB_FLAT_FACE_SPHERE_REF_DIRECTION] = {x[0], x[1], x[2]};
    a.f[SMLIB_FLAT_FACE_SPHERE_RADIUS] = {r};
    a.f[SMLIB_FLAT_FACE_RANGE] = {u0, -0.5 * pi, u1, 0.5 * pi};
    // Vertices: the south pole, the north pole.
    for (double side : {-1.0, 1.0})
    {
        for (int k = 0; k < 3; ++k)
        {
            a.f[SMLIB_FLAT_VERTEX_POSITION].push_back(o[k] + side * r * z[k]);
        }
        a.i[SMLIB_FLAT_VERTEX_POINT_TYPE].push_back(1);
        a.f[SMLIB_FLAT_VERTEX_TOLERANCE].push_back(1e-5);
    }
    // The meridian at u: reference direction radial(u), axis radial(u) cross z (its y axis is z),
    // angle v. The loop runs counter-clockwise in (u, v): north along the meridian at u1, south
    // along (against) the one at u0.
    for (double u : {u1, u0})
    {
        const double t[3] = {std::cos(u) * x[0] + std::sin(u) * y[0], std::cos(u) * x[1] + std::sin(u) * y[1], std::cos(u) * x[2] + std::sin(u) * y[2]};
        const double n[3] = {t[1] * z[2] - t[2] * z[1], t[2] * z[0] - t[0] * z[2], t[0] * z[1] - t[1] * z[0]};
        for (int k = 0; k < 3; ++k)
        {
            a.f[SMLIB_FLAT_EDGE_CIRCLE_CENTER].push_back(o[k]);
            a.f[SMLIB_FLAT_EDGE_CIRCLE_AXIS].push_back(n[k]);
            a.f[SMLIB_FLAT_EDGE_CIRCLE_REF_DIRECTION].push_back(t[k]);
        }
        a.f[SMLIB_FLAT_EDGE_CIRCLE_RADIUS].push_back(r);
        a.f[SMLIB_FLAT_EDGE_RANGE].insert(a.f[SMLIB_FLAT_EDGE_RANGE].end(), {-0.5 * pi, 0.5 * pi});
    }
    a.i[SMLIB_FLAT_EDGE_VERTICES] = {0, 1, 0, 1};
    for (int64_t k = 0; k < 2; ++k)
    {
        a.i[SMLIB_FLAT_EDGE_CURVE_TYPE].push_back(SMLIB_FLAT_CURVE_CIRCLE);
        a.f[SMLIB_FLAT_EDGE_TOLERANCE].push_back(2e-5);
        a.i[SMLIB_FLAT_EDGEUSE_EDGE].push_back(k);
        a.i[SMLIB_FLAT_EDGEUSE_ORIENTATION].push_back(k == 0 ? 1 : 0);
        a.i[SMLIB_FLAT_EDGEUSE_NEXT_RADIAL].push_back(k);
        a.i[SMLIB_FLAT_EDGEUSE_RADIAL_ENTRY].push_back(0);
    }
    a.f[SMLIB_FLAT_TOLERANCE] = {4e-5};
    return a;
}

// The same for a sphere face from pole to pole: SmSphere::IsNurbSurfaceSphere took a half meridian
// it measured in the form at 180.0000000000024 degrees (a lune of a radius 2 sphere 500 mm out, on
// an engine's valve cover) for more than a half circle and said no, and the import refused it.
void SphereFormRoundTrip()
{
    const double pi = 3.14159265358979323846, deg = pi / 180.0;
    const double o[3] = {1, 2, 3}, far[3] = {107.214217174, 350.149503622, -367.5}, up[3] = {0, 0, 1}, ref[3] = {1, 0, 0};
    FormRoundTrip(SphereLune(o, up, ref, 10.0, 0.0, 90.0 * deg), "sphere face", "quarter lune", SMLIB_FLAT_FACE_SPHERE_AXIS,
                  SMLIB_FLAT_FACE_SPHERE_REF_DIRECTION);
    FormRoundTrip(SphereLune(o, up, ref, 10.0, 30.0 * deg, 150.0 * deg), "sphere face", "lune off the seam", SMLIB_FLAT_FACE_SPHERE_AXIS,
                  SMLIB_FLAT_FACE_SPHERE_REF_DIRECTION);
    FormRoundTrip(SphereLune(far, up, ref, 2.0, 0.0, 90.0 * deg), "sphere face", "small lune, far out", SMLIB_FLAT_FACE_SPHERE_AXIS,
                  SMLIB_FLAT_FACE_SPHERE_REF_DIRECTION);
}

// A cylinder mirrored in a plane (r 5 x 30 on (20, 20, -5) up z, mirrored in x = 30): its NURBS
// sweeps clockwise about its linear direction, so the cylinder is inside out, its STEP origin on
// the NURBS's top circle (40, 20, 25) and its STEP z axis -z. Its export imports again as the same
// cylinder: SmCylinder::IsNurbSurfaceCylinder put the origin on the bottom circle with z along the
// linear direction, a cylinder sweeping clockwise about its own z whose generator ran 30 mm off
// the NURBS, and SMLib's check refused it.
void MirroredCylinderRoundTrip()
{
    const double pi = 3.14159265358979323846;
    const double o[3] = {20, 20, -5};
    const double mirror[4][4] = {{-1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {60, 0, 0, 1}};
    smlib_brep *cyl = nullptr, *mirrored = nullptr, *again = nullptr;
    smlib_flat *x = nullptr, *y = nullptr;
    CHECK(smlib_make_cylinder(o, 5, 30, &cyl) == SMLIB_OK);
    CHECK(smlib_brep_transform(cyl, mirror, &mirrored) == SMLIB_OK);
    const int32_t sx = smlib_brep_export(mirrored, &x);
    // SMLib's primitive side is a cone of semi-angle 0; it goes back labelled a cylinder, as
    // SolidCraft's records and STEP name it.
    Arrays first = sx == SMLIB_OK ? Snapshot(x) : Arrays{};
    for (int64_t &t : first.i[SMLIB_FLAT_FACE_SURFACE_TYPE])
    {
        t = t == SMLIB_FLAT_SURFACE_CONE ? SMLIB_FLAT_SURFACE_CYLINDER : t;
    }
    for (const auto &[cone, cylinder] : {std::pair{SMLIB_FLAT_FACE_CONE_ORIGIN, SMLIB_FLAT_FACE_CYLINDER_ORIGIN},
                                         std::pair{SMLIB_FLAT_FACE_CONE_AXIS, SMLIB_FLAT_FACE_CYLINDER_AXIS},
                                         std::pair{SMLIB_FLAT_FACE_CONE_REF_DIRECTION, SMLIB_FLAT_FACE_CYLINDER_REF_DIRECTION},
                                         std::pair{SMLIB_FLAT_FACE_CONE_RADIUS, SMLIB_FLAT_FACE_CYLINDER_RADIUS}})
    {
        first.f[cylinder] = first.f[cone];
        first.f[cone].clear();
    }
    first.f[SMLIB_FLAT_FACE_CONE_SEMI_ANGLE].clear();
    smlib_flat *given = FromArrays(first);
    const int32_t sa = sx == SMLIB_OK ? smlib_brep_import(given, &again) : sx;
    const int32_t sy = sa == SMLIB_OK ? smlib_brep_export(again, &y) : sa;
    char buf[4096] = {0};
    int64_t len = 0, bad = -1;
    double v = 0.0;
    if (sy == SMLIB_OK)
    {
        CHECK(smlib_brep_check(again, 2, buf, sizeof buf, &len, &bad) == SMLIB_OK);
        CHECK(smlib_brep_volume(again, 1e-8, &v) == SMLIB_OK);
    }
    const Arrays second = sy == SMLIB_OK ? Snapshot(y) : Arrays{};
    const std::vector<double> want_origin = {40, 20, 25}, want_axis = {0, 0, -1};
    double frame = sy == SMLIB_OK ? 0.0 : INFINITY;
    for (const Arrays *a : {static_cast<const Arrays *>(&first), &second})
    {
        for (const auto &[id, want] : {std::pair{SMLIB_FLAT_FACE_CYLINDER_ORIGIN, want_origin}, std::pair{SMLIB_FLAT_FACE_CYLINDER_AXIS, want_axis}})
        {
            // Every face of the cylinder's side, each its own entry.
            const auto at = a->f.find(id);
            const bool some = at != a->f.end() && !at->second.empty() && at->second.size() % 3 == 0;
            frame = some ? frame : INFINITY;
            for (size_t k = 0; some && k < at->second.size(); ++k)
            {
                frame = std::max(frame, std::fabs(at->second[k] - want[k % 3]));
            }
        }
    }
    const double want = pi * 25.0 * 30.0;
    std::printf("mirrored cylinder: export %d, imported again %d, exported again %d, failed checks %lld %s, volume %.12g (want %.12g), "
                "STEP frame off by %.3g\n",
                static_cast<int>(sx), static_cast<int>(sa), static_cast<int>(sy), static_cast<long long>(bad), buf, v, want, frame);
    CHECK(sy == SMLIB_OK && bad == 0 && std::fabs(v - want) <= 1e-8 * want && frame <= 1e-12);
    smlib_flat_free(y);
    smlib_flat_free(given);
    smlib_flat_free(x);
    for (smlib_brep *b : {cyl, mirrored, again})
    {
        smlib_brep_free(b);
    }
}

// Imports `a`; expects a failure status (and no crash). Returns the status.
int32_t ExpectReject(const Arrays &a, const char *what)
{
    smlib_flat *flat = FromArrays(a);
    smlib_brep *b = nullptr;
    int32_t s = smlib_brep_import(flat, &b);
    if (s == SMLIB_OK)
    {
        std::fprintf(stderr, "  hostile case accepted: %s\n", what);
        failures++;
    }
    CHECK(s == SMLIB_OK || b == nullptr);
    smlib_brep_free(b);
    smlib_flat_free(flat);
    return s;
}

// Truncations, bad indices, bad values and huge counts, all of which must be rejected.
int Hostile(const char *name, const smlib_flat *flat)
{
    const Arrays base = Snapshot(flat);
    int cases = 0;
    char what[128];

    // Each non-empty array one element short (SOURCE may be empty: it is optional).
    for (auto &kv : base.i)
    {
        if (kv.second.empty() || kv.first == SMLIB_FLAT_SOURCE)
        {
            continue;
        }
        Arrays a = base;
        a.i[kv.first].pop_back();
        std::snprintf(what, sizeof what, "%s: int array %d truncated", name, kv.first);
        ExpectReject(a, what);
        cases++;
    }
    // (EDGE_ and VERTEX_TOLERANCE may be empty: they are optional.)
    for (auto &kv : base.f)
    {
        bool optional = kv.first == SMLIB_FLAT_EDGE_TOLERANCE || kv.first == SMLIB_FLAT_VERTEX_TOLERANCE;
        if (kv.second.empty() || (optional && kv.second.size() == 1))
        {
            continue;
        }
        Arrays a = base;
        a.f[kv.first].pop_back();
        std::snprintf(what, sizeof what, "%s: real array %d truncated", name, kv.first);
        ExpectReject(a, what);
        cases++;
    }

    // Index, kind and count arrays: values out of range and huge.
    const int32_t indexIds[] = {SMLIB_FLAT_REGION_SHELL_COUNT, SMLIB_FLAT_REGION_TYPE,       SMLIB_FLAT_SHELL_FACEUSE_COUNT,
                                SMLIB_FLAT_FACEUSE_FACE,       SMLIB_FLAT_FACEUSE_ORIENTATION, SMLIB_FLAT_FACE_LOOP_COUNT,
                                SMLIB_FLAT_FACE_SURFACE_TYPE,  SMLIB_FLAT_LOOP_EDGEUSE_COUNT, SMLIB_FLAT_LOOP_VERTEX,
                                SMLIB_FLAT_EDGEUSE_EDGE,       SMLIB_FLAT_EDGEUSE_NEXT_RADIAL, SMLIB_FLAT_EDGE_CURVE_TYPE,
                                SMLIB_FLAT_EDGE_VERTICES,      SMLIB_FLAT_EDGE_NURBS_VERTEX_COUNT, SMLIB_FLAT_EDGE_NURBS_ORDER,
                                SMLIB_FLAT_FACE_NURBS_U_VERTEX_COUNT, SMLIB_FLAT_FACE_NURBS_V_ORDER, SMLIB_FLAT_EDGEUSE_NURBS_ORDER,
                                SMLIB_FLAT_VERTEX_POINT_TYPE};
    const int64_t badValues[] = {-1, -2, 7, 99, int64_t(1) << 40, std::numeric_limits<int64_t>::max(), std::numeric_limits<int64_t>::min()};
    for (int32_t id : indexIds)
    {
        if (base.i.at(id).empty())
        {
            continue;
        }
        for (int64_t v : badValues)
        {
            Arrays a = base;
            int64_t &x = a.i[id][a.i[id].size() / 2];
            if (x == v || (id == SMLIB_FLAT_LOOP_VERTEX && v == -1))
            {
                continue;
            }
            x = v;
            std::snprintf(what, sizeof what, "%s: int array %d set to %lld", name, id, (long long)v);
            ExpectReject(a, what);
            cases++;
        }
    }

    // Non-finite and absurd reals everywhere.
    const double badReals[] = {NAN, INFINITY, -INFINITY, 1e300};
    for (auto &kv : base.f)
    {
        if (kv.second.empty() || kv.first == SMLIB_FLAT_EXTENT)
        {
            continue;
        }
        for (double v : badReals)
        {
            Arrays a = base;
            a.f[kv.first][0] = v;
            std::snprintf(what, sizeof what, "%s: real array %d [0] set to %g", name, kv.first, v);
            ExpectReject(a, what);
            cases++;
        }
    }
    {
        // Tolerance zero, an empty set, and a mis-sized extent.
        Arrays a = base;
        a.f[SMLIB_FLAT_TOLERANCE][0] = 0.0;
        ExpectReject(a, "zero tolerance");
        Arrays e;
        ExpectReject(e, "empty");
        Arrays x = base;
        x.f[SMLIB_FLAT_EXTENT].push_back(1.0);
        ExpectReject(x, "extent of 7");
        // Export writes a tolerance per edge and per vertex; one too many is refused.
        CHECK(base.f.at(SMLIB_FLAT_EDGE_TOLERANCE).size() == base.i.at(SMLIB_FLAT_EDGE_CURVE_TYPE).size());
        CHECK(base.f.at(SMLIB_FLAT_VERTEX_TOLERANCE).size() == base.i.at(SMLIB_FLAT_VERTEX_POINT_TYPE).size());
        Arrays t = base;
        t.f[SMLIB_FLAT_EDGE_TOLERANCE].push_back(1e-5);
        ExpectReject(t, "an edge tolerance too many");
        Arrays z = base;
        if (!z.f[SMLIB_FLAT_VERTEX_TOLERANCE].empty())
        {
            z.f[SMLIB_FLAT_VERTEX_TOLERANCE][0] = 0.0;
            ExpectReject(z, "a zero vertex tolerance");
        }
        cases += 5;
    }
    return cases;
}

// Random single-value mutations: any status is fine, a crash is not.
void Mutations(const smlib_flat *flat, uint64_t seed, int rounds, int &accepted)
{
    const Arrays base = Snapshot(flat);
    std::vector<int32_t> ints, reals;
    for (auto &kv : base.i)
    {
        if (!kv.second.empty())
        {
            ints.push_back(kv.first);
        }
    }
    for (auto &kv : base.f)
    {
        if (!kv.second.empty())
        {
            reals.push_back(kv.first);
        }
    }
    uint64_t r = seed;
    auto next = [&r]
    {
        r ^= r << 13;
        r ^= r >> 7;
        r ^= r << 17;
        return r;
    };
    for (int k = 0; k < rounds; ++k)
    {
        Arrays a = base;
        if (next() % 2 == 0)
        {
            std::vector<int64_t> &v = a.i[ints[next() % ints.size()]];
            int64_t &x = v[next() % v.size()];
            switch (next() % 4)
            {
            case 0: x += 1; break;
            case 1: x -= 1; break;
            case 2: x = static_cast<int64_t>(next() % 64); break;
            default: x = static_cast<int64_t>(next()); break;
            }
        }
        else
        {
            std::vector<double> &v = a.f[reals[next() % reals.size()]];
            double &x = v[next() % v.size()];
            switch (next() % 4)
            {
            case 0: x = -x; break;
            case 1: x *= 1.0 + 1e-3; break;
            case 2: x += static_cast<double>(next() % 1000) / 100.0 - 5.0; break;
            default: x = 0.0; break;
            }
        }
        smlib_flat *m = FromArrays(a);
        smlib_brep *b = nullptr;
        if (smlib_brep_import(m, &b) == SMLIB_OK)
        {
            accepted++;
            double vol = 0.0;
            smlib_brep_volume(b, 1e-3, &vol);
        }
        smlib_brep_free(b);
        smlib_flat_free(m);
    }
}

} // namespace

int main(int argc, char **argv)
{
    if (argc > 1)
    {
        dir = argv[1];
    }
    CHECK(smlib_abi_version() == SMLIB_ABI_VERSION && SMLIB_ABI_VERSION == 4);
    CHECK(smlib_init() == SMLIB_OK);

    // The solids. With SMLIB_FLAT_SOLIDS set to a directory of <name>.smb files from an earlier
    // run, they are read from there instead of built (for sanitizer runs: SMLib's boolean and
    // fillet code has memory errors of its own that would stop the run before the converter).
    std::vector<std::pair<const char *, smlib_brep *>> solids;
    if (const char *from = std::getenv("SMLIB_FLAT_SOLIDS"))
    {
        for (const char *name : {"box", "cylinder", "cone", "sphere", "torus", "fillet", "blend", "blend1", "boolean", "incline"})
        {
            smlib_brep *b = nullptr;
            CHECK(smlib_brep_read((std::string(from) + "/" + name + ".smb").c_str(), 0, &b) == SMLIB_OK);
            solids.emplace_back(name, b);
        }
    }
    else
    {
    {
        SmBrep *b = nullptr;
        SmApiCreateBox(SmVector3d(1, 2, 3), 10, 20, 30, b);
        solids.emplace_back("box", ToCApi(b, "box"));
    }
    {
        SmBrep *b = nullptr;
        SmApiCreateCylinder(SmVector3d(0, 0, 0), 5, 12, b);
        solids.emplace_back("cylinder", ToCApi(b, "cylinder"));
    }
    {
        SmBrep *b = nullptr;
        SmApiCreateCone(SmVector3d(1, 1, 1), 6, 2, 9, b);
        solids.emplace_back("cone", ToCApi(b, "cone"));
    }
    {
        SmBrep *b = nullptr;
        SmApiCreateSphere(SmVector3d(0, 0, 0), 7, b);
        solids.emplace_back("sphere", ToCApi(b, "sphere"));
    }
    {
        SmBrep *b = nullptr;
        SmApiCreateTorus(SmVector3d(0, 0, 0), 10, 3, b);
        solids.emplace_back("torus", ToCApi(b, "torus"));
    }
    {
        SmBrep *b = nullptr;
        SmApiCreateBox(SmVector3d(0, 0, 0), 10, 10, 10, b);
        SmTArray<SmEdge *> edges;
        SmApiGetEdges(b, edges);
        SmApiStatus s = SmApiFilletEdges(b, edges, 1.0, 1, 1, 1.0);
        if (s != SM_SUCCESS)
        {
            std::fprintf(stderr, "fillet: SmApiFilletEdges -> %d\n", (int)s);
            failures++;
        }
        solids.emplace_back("fillet", ToCApi(b, "fillet"));
    }
    {
        // Blend-curve cross sections give NURBS fillet faces, so UV trims are carried too.
        SmBrep *b = nullptr;
        SmApiCreateBox(SmVector3d(0, 0, 0), 10, 10, 10, b);
        SmTArray<SmEdge *> edges;
        SmApiGetEdges(b, edges);
        SmApiStatus s = SmApiFilletEdges(b, edges, 1.0, 2, 2, 1.0);
        if (s != SM_SUCCESS)
        {
            std::fprintf(stderr, "blend: SmApiFilletEdges -> %d\n", (int)s);
            failures++;
        }
        solids.emplace_back("blend", ToCApi(b, "blend"));
    }
    {
        // One blended edge: one NURBS face with UV trims. (Its four boundary edges come out as
        // B-splines, so this still carries no trim on an analytic edge.)
        SmBrep *b = nullptr;
        SmApiCreateBox(SmVector3d(0, 0, 0), 10, 10, 10, b);
        SmTArray<SmEdge *> edges, one;
        SmApiGetEdges(b, edges);
        one.Add(edges[0]);
        SmApiStatus s = SmApiFilletEdges(b, one, 2.0, 2, 2, 1.0);
        if (s != SM_SUCCESS)
        {
            std::fprintf(stderr, "blend1: SmApiFilletEdges -> %d\n", (int)s);
            failures++;
        }
        solids.emplace_back("blend1", ToCApi(b, "blend1"));
    }
    {
        SmBrep *box = nullptr, *cyl = nullptr, *r = nullptr;
        SmApiCreateBox(SmVector3d(0, 0, 0), 20, 20, 20, box);
        SmApiCreateCylinder(SmVector3d(10, 10, -5), 4, 30, cyl);
        SmApiStatus s = SmApiBooleanDifference(box, cyl, r);
        if (s != SM_SUCCESS)
        {
            std::fprintf(stderr, "boolean: SmApiBooleanDifference -> %d\n", (int)s);
            failures++;
            delete cyl;
            delete box;
            r = nullptr;
        }
        solids.emplace_back("boolean", ToCApi(r, "boolean"));
    }
    {
        // A cylinder minus a tilted box: the inclined rim is made of ellipse arcs. (SmApiCut with
        // an inclined plane also makes an ellipse, but leaves the body open: it does not cap.)
        SmBrep *cyl = nullptr, *box = nullptr, *r = nullptr;
        SmApiCreateCylinder(SmVector3d(0, 0, 0), 5, 20, cyl);
        SmApiCreateBox(SmVector3d(-20, -20, 10), 40, 40, 40, box);
        SmApiRotate(box, SmVector3d(0, 0, 10), SmVector3d(0, 1, 0), 25.0);
        SmApiStatus s = SmApiBooleanDifference(cyl, box, r);
        if (s != SM_SUCCESS)
        {
            std::fprintf(stderr, "incline: SmApiBooleanDifference -> %d\n", (int)s);
            failures++;
            delete cyl;
            delete box;
            r = nullptr;
        }
        solids.emplace_back("incline", ToCApi(r, "incline"));
    }
    }

    std::vector<std::pair<const char *, smlib_flat *>> exports;
    for (auto &s : solids)
    {
        const bool cone = std::strcmp(s.first, "cone") == 0;
        const bool reparameterised = cone || std::strcmp(s.first, "sphere") == 0 || std::strcmp(s.first, "incline") == 0;
        smlib_flat *f = RoundTrip(s.first, s.second, cone, reparameterised ? 1e-7 : 1e-9);
        if (f != nullptr)
        {
            exports.emplace_back(s.first, f);
        }
    }

    // An analytic face's angular range may start a turn past the principal one (solidcore#428):
    // the torus with its major, minor or both ranges moved on by 2 pi is the same torus. Its
    // tube circle was once rebuilt over a clamped, empty parameter interval, which gave NaN
    // knots, an overflow in N_IgesClampSrf (#427) and then a refused import.
    //
    // The angles are the range of a face given without its NURBS form (as SolidCraft designs and
    // other producers write a torus); export gives every torus face its form, whose range is in the
    // NURBS's parameters, where a turn on is off the surface. So the face goes in as such a
    // producer writes it: canonical only, its range the STEP angles that start at its seam vertex,
    // no UV trims (SMLib makes them from the edges). Shift 0 is that face as it is.
    for (auto &e : exports)
    {
        if (std::strcmp(e.first, "torus") != 0)
        {
            continue;
        }
        const double pi = 3.14159265358979323846;
        Arrays canonical = Snapshot(e.second);
        CHECK(canonical.i[SMLIB_FLAT_FACE_SURFACE_TYPE] == std::vector<int64_t>{SMLIB_FLAT_SURFACE_TORUS});
        CHECK(canonical.f[SMLIB_FLAT_VERTEX_POSITION].size() == 3);
        if (canonical.i[SMLIB_FLAT_FACE_SURFACE_TYPE].size() != 1 || canonical.f[SMLIB_FLAT_VERTEX_POSITION].size() != 3)
        {
            continue;
        }
        {
            // The seam vertex's angles in the torus's frame.
            const std::vector<double> &o = canonical.f[SMLIB_FLAT_FACE_TORUS_ORIGIN], &z = canonical.f[SMLIB_FLAT_FACE_TORUS_AXIS],
                                      &x = canonical.f[SMLIB_FLAT_FACE_TORUS_REF_DIRECTION], &p = canonical.f[SMLIB_FLAT_VERTEX_POSITION];
            const double R = canonical.f[SMLIB_FLAT_FACE_TORUS_MAJOR_RADIUS][0];
            const double y[3] = {z[1] * x[2] - z[2] * x[1], z[2] * x[0] - z[0] * x[2], z[0] * x[1] - z[1] * x[0]};
            const double d[3] = {p[0] - o[0], p[1] - o[1], p[2] - o[2]};
            auto dot = [](const double *a, const double *b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; };
            const double u0 = std::atan2(dot(d, y), dot(d, x.data()));
            const double radial[3] = {std::cos(u0) * x[0] + std::sin(u0) * y[0], std::cos(u0) * x[1] + std::sin(u0) * y[1],
                                      std::cos(u0) * x[2] + std::sin(u0) * y[2]};
            const double q[3] = {d[0] - R * radial[0], d[1] - R * radial[1], d[2] - R * radial[2]};
            const double v0 = std::atan2(dot(q, z.data()), dot(q, radial));
            canonical.f[SMLIB_FLAT_FACE_RANGE] = {u0, v0, u0 + 2.0 * pi, v0 + 2.0 * pi};
        }
        for (int32_t id : {SMLIB_FLAT_FACE_NURBS_FORM, SMLIB_FLAT_FACE_NURBS_U_VERTEX_COUNT, SMLIB_FLAT_FACE_NURBS_V_VERTEX_COUNT, SMLIB_FLAT_FACE_NURBS_U_ORDER,
                           SMLIB_FLAT_FACE_NURBS_V_ORDER, SMLIB_FLAT_EDGEUSE_NURBS_VERTEX_COUNT, SMLIB_FLAT_EDGEUSE_NURBS_ORDER})
        {
            canonical.i[id].clear();
        }
        for (int32_t id : {SMLIB_FLAT_FACE_NURBS_CVS, SMLIB_FLAT_FACE_NURBS_U_KNOTS, SMLIB_FLAT_FACE_NURBS_V_KNOTS, SMLIB_FLAT_FACE_NURBS_WEIGHTS,
                           SMLIB_FLAT_EDGEUSE_NURBS_CVS, SMLIB_FLAT_EDGEUSE_NURBS_KNOTS, SMLIB_FLAT_EDGEUSE_NURBS_WEIGHTS})
        {
            canonical.f[id].clear();
        }
        for (int shift = 0; shift <= 3; ++shift)
        {
            Arrays a = canonical;
            std::vector<double> &range = a.f[SMLIB_FLAT_FACE_RANGE];
            for (int c = 0; c < 2; ++c)
            {
                if (shift & (1 << c))
                {
                    range[c] += 2.0 * pi;
                    range[2 + c] += 2.0 * pi;
                }
            }
            smlib_flat *f = FromArrays(a);
            smlib_brep *b = nullptr;
            int32_t s = smlib_brep_import(f, &b);
            const double want = 2.0 * pi * pi * 10.0 * 9.0;
            double v = 0.0;
            int64_t len = 0, bad = -1;
            char buf[4096] = {0};
            if (s == SMLIB_OK)
            {
                v = Measure(b).volume;
                CHECK(smlib_brep_check(b, 2, buf, sizeof buf, &len, &bad) == SMLIB_OK);
            }
            std::printf("torus, canonical, %s range a turn on: import %d, volume %.12g (want %.12g), failed checks %lld %s\n",
                        shift == 0 ? "no" : shift == 1 ? "major" : shift == 2 ? "minor" : "both", static_cast<int>(s), v, want, static_cast<long long>(bad), buf);
            CHECK(s == SMLIB_OK && std::fabs(v - want) <= 1e-7 * want && bad == 0);
            // Its export (the label and its NURBS form) imports again.
            smlib_flat *x = nullptr;
            smlib_brep *again = nullptr;
            const int32_t sx = s == SMLIB_OK ? smlib_brep_export(b, &x) : s;
            const int32_t sa = sx == SMLIB_OK ? smlib_brep_import(x, &again) : sx;
            std::printf("torus, canonical, %s range a turn on: export %d, imported again %d\n",
                        shift == 0 ? "no" : shift == 1 ? "major" : shift == 2 ? "minor" : "both", static_cast<int>(sx), static_cast<int>(sa));
            CHECK(sa == SMLIB_OK);
            smlib_brep_free(again);
            smlib_flat_free(x);
            smlib_brep_free(b);
            smlib_flat_free(f);
        }
    }

    NestedShellOrder();
    BoundaryRowTrimPace();
    EdgePastPatch();
    TorusFormRoundTrip();
    SphereFormRoundTrip();
    MirroredCylinderRoundTrip();

    // The accessors' argument checks.
    {
        smlib_flat *f = smlib_flat_new();
        const int64_t *pi = nullptr;
        const double *pf = nullptr;
        int64_t n = 0;
        const int64_t one = 1;
        const double half = 0.5;
        for (int32_t id : {0, 36, 100, 166, -1, 1000})
        {
            CHECK(smlib_flat_get_i64(f, id, &pi, &n) == SMLIB_ERR_ARG);
            CHECK(smlib_flat_get_f64(f, id, &pf, &n) == SMLIB_ERR_ARG);
            CHECK(smlib_flat_set_i64(f, id, &one, 1) == SMLIB_ERR_ARG);
            CHECK(smlib_flat_set_f64(f, id, &half, 1) == SMLIB_ERR_ARG);
        }
        CHECK(smlib_flat_get_f64(f, SMLIB_FLAT_REGION_TYPE, &pf, &n) == SMLIB_ERR_ARG);
        CHECK(smlib_flat_get_i64(f, SMLIB_FLAT_TOLERANCE, &pi, &n) == SMLIB_ERR_ARG);
        CHECK(smlib_flat_set_f64(f, SMLIB_FLAT_REGION_TYPE, &half, 1) == SMLIB_ERR_ARG);
        CHECK(smlib_flat_set_i64(f, SMLIB_FLAT_TOLERANCE, &one, 1) == SMLIB_ERR_ARG);
        CHECK(smlib_flat_set_i64(f, SMLIB_FLAT_REGION_TYPE, &one, -1) == SMLIB_ERR_ARG);
        CHECK(smlib_flat_set_i64(f, SMLIB_FLAT_REGION_TYPE, &one, SMLIB_FLAT_MAX_LEN + 1) == SMLIB_ERR_ARG);
        CHECK(smlib_flat_set_i64(f, SMLIB_FLAT_REGION_TYPE, nullptr, 3) == SMLIB_ERR_ARG);
        CHECK(smlib_flat_set_i64(f, SMLIB_FLAT_REGION_TYPE, nullptr, 0) == SMLIB_OK);
        CHECK(smlib_flat_set_i64(f, SMLIB_FLAT_REGION_TYPE, &one, 1) == SMLIB_OK);
        CHECK(smlib_flat_get_i64(f, SMLIB_FLAT_REGION_TYPE, &pi, &n) == SMLIB_OK && n == 1 && pi[0] == 1);
        CHECK(smlib_flat_get_i64(nullptr, SMLIB_FLAT_REGION_TYPE, &pi, &n) == SMLIB_ERR_ARG);
        smlib_brep *b = nullptr;
        CHECK(smlib_brep_import(f, &b) == SM_ERR_INVALID_INPUT && b == nullptr);
        CHECK(smlib_brep_import(nullptr, &b) == SMLIB_ERR_ARG);
        CHECK(smlib_brep_export(nullptr, &f) == SMLIB_ERR_ARG);
        CHECK(smlib_flat_read((dir + "/no-such-file").c_str(), &f) == SM_ERR_INVALID_INPUT);
        smlib_flat_free(f);
        smlib_flat_free(nullptr);
    }

    // A failing import dumps its arrays; the dump reads back as the same arrays.
    if (!exports.empty())
    {
        std::string dumps = dir + "/dumps";
        std::error_code ec;
        std::filesystem::remove_all(dumps, ec);
        CHECK(std::filesystem::create_directories(dumps, ec) && !ec);
        CHECK(smlib_set_dump_dir(dumps.c_str()) == SMLIB_OK);
        Arrays bad = Snapshot(exports[0].second);
        bad.i[SMLIB_FLAT_EDGEUSE_EDGE][0] = 1000;
        smlib_flat *f = FromArrays(bad);
        smlib_brep *b = nullptr;
        CHECK(smlib_brep_import(f, &b) == SM_ERR_INVALID_INPUT);
        CHECK(smlib_set_dump_dir(nullptr) == SMLIB_OK);
        smlib_flat *replay = nullptr;
        bool found = false;
        for (int n = 1; n < 100000 && !found; ++n)
        {
            std::string path = dumps + "/" + std::to_string(n) + "-smlib_brep_import-flat.txt";
            if (smlib_flat_read(path.c_str(), &replay) == SMLIB_OK)
            {
                found = true;
            }
        }
        CHECK(found);
        if (found)
        {
            Arrays r = Snapshot(replay);
            CHECK(r.i == bad.i && r.f == bad.f);
            CHECK(smlib_brep_import(replay, &b) == SM_ERR_INVALID_INPUT);
        }
        smlib_flat_free(replay);
        smlib_flat_free(f);
    }

    int hostile = 0, accepted = 0, mutated = 0;
    uint64_t seed = 0x9E3779B97F4A7C15ull;
    for (auto &e : exports)
    {
        hostile += Hostile(e.first, e.second);
        Mutations(e.second, seed++, 300, accepted);
        mutated += 300;
    }
    std::printf("hostile cases: %d rejected as expected; random mutations: %d, of which %d still imported\n", hostile, mutated, accepted);

    for (auto &e : exports)
    {
        smlib_flat_free(e.second);
    }
    if (failures == 0)
    {
        std::printf("smlib flat round trip: ok\n");
    }
    return failures == 0 ? 0 : 1;
}
