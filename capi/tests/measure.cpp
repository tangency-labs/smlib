// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

// smlib_brep_mass_properties and smlib_brep_area against closed-form values: boxes, cylinders,
// spheres, cones, frustums and tori, at the origin and far from it. Prints the worst relative
// error of each quantity.

#include "../src/smlib_shim.h"

#include <SmApiBrep.h>
#include <SmApiPrimitives.h>
#include <SmApiTrimmedSurfaces.h>
#include <SmAxis2Placement.h>
#include <SmBSplineSurface.h>
#include <SmCircle.h>
#include <SmLine.h>

#include <chrono>
#include <cmath>
#include <cstdio>
#include <functional>
#include <vector>

static int failures = 0;

#define CHECK(cond)                                                                                \
    do                                                                                             \
    {                                                                                              \
        if (!(cond))                                                                               \
        {                                                                                          \
            std::fprintf(stderr, "%s:%d: CHECK(%s) failed\n", __FILE__, __LINE__, #cond);          \
            failures++;                                                                            \
        }                                                                                          \
    } while (0)

static const double PI = 3.14159265358979323846;

// Closed-form values: area, volume, centroid, Ixx Iyy Izz and Iyz Izx Ixy about the centroid.
struct Expected
{
    const char *name;
    std::function<SmApiStatus(SmBrep *&)> make;
    double v[11];
};

static Expected box(const char *name, SmVector3d o, double a, double b, double c)
{
    double V = a * b * c;
    return {name,
            [=](SmBrep *&p) { return SmApiCreateBox(o, a, b, c, p); },
            {2 * (a * b + b * c + c * a), V, o.x + a / 2, o.y + b / 2, o.z + c / 2, V * (b * b + c * c) / 12,
             V * (a * a + c * c) / 12, V * (a * a + b * b) / 12, 0, 0, 0}};
}

static Expected cylinder(const char *name, SmVector3d o, double r, double h)
{
    double V = PI * r * r * h;
    double Ixx = V * (3 * r * r + h * h) / 12;
    return {name,
            [=](SmBrep *&p) { return SmApiCreateCylinder(o, r, h, p); },
            {2 * PI * r * r + 2 * PI * r * h, V, o.x, o.y, o.z + h / 2, Ixx, Ixx, V * r * r / 2, 0, 0, 0}};
}

static Expected sphere(const char *name, SmVector3d o, double r)
{
    double V = 4.0 / 3.0 * PI * r * r * r;
    double I = 0.4 * V * r * r;
    return {name, [=](SmBrep *&p) { return SmApiCreateSphere(o, r, p); }, {4 * PI * r * r, V, o.x, o.y, o.z, I, I, I, 0, 0, 0}};
}

// A frustum (or a cone, rt = 0) on its base centre o, radius rb at z = 0 to rt at z = h. Its
// volume integrals are polynomials in z of degree at most 4, so three-point Gauss-Legendre gives
// them exactly.
static Expected cone(const char *name, SmVector3d o, double rb, double rt, double h)
{
    const double gx[3] = {-std::sqrt(0.6), 0.0, std::sqrt(0.6)}, gw[3] = {5.0 / 9, 8.0 / 9, 5.0 / 9};
    double V = 0, Mz = 0, Izz = 0, Zz = 0;
    for (int i = 0; i < 3; ++i)
    {
        double z = h / 2 * (1 + gx[i]), w = gw[i] * h / 2;
        double rho = rb + (rt - rb) * z / h, disc = PI * rho * rho;
        V += w * disc;
        Mz += w * disc * z;
        Zz += w * disc * z * z;
        Izz += w * PI * rho * rho * rho * rho / 2;
    }
    double zc = Mz / V;
    double Ixx = Izz / 2 + Zz - V * zc * zc;
    double s = std::sqrt((rb - rt) * (rb - rt) + h * h);
    return {name,
            [=](SmBrep *&p) { return SmApiCreateCone(o, rb, rt, h, p); },
            {PI * rb * rb + PI * rt * rt + PI * (rb + rt) * s, V, o.x, o.y, o.z + zc, Ixx, Ixx, Izz, 0, 0, 0}};
}

// A torus about the z axis through o (Pappus for area and volume).
static Expected torus(const char *name, SmVector3d o, double R, double r)
{
    double V = 2 * PI * PI * R * r * r;
    double Ixx = V * (R * R / 2 + 5 * r * r / 8);
    return {name,
            [=](SmBrep *&p) { return SmApiCreateTorus(o, R, r, p); },
            {4 * PI * PI * R * r, V, o.x, o.y, o.z, Ixx, Ixx, V * (R * R + 3 * r * r / 4), 0, 0, 0}};
}

int main()
{
    CHECK(smlib_init() == SMLIB_OK);
    const SmVector3d at0(0, 0, 0), far(1000.0, -2000.0, 500.0);
    Expected cases[] = {
        // tests/oracle (SolidCraft): 01-box, 02-cylinder, 38-coincident-cylinders.
        box("oracle 01-box 40x30x20", at0, 40, 30, 20),
        cylinder("oracle 02-cylinder r15 h50", at0, 15, 50),
        cylinder("oracle 38 cylinder r10 h50", at0, 10, 50),
        // The primitives of oracle 40-box-sphere-torus, whole.
        sphere("sphere r12", at0, 12),
        torus("torus R22 r4", at0, 22, 4),
        cone("frustum r10..4 h20", at0, 10, 4, 20),
        cone("cone r10 h20", at0, 10, 0, 20),
        cone("inverted frustum r3..9 h7", at0, 3, 9, 7),
        // The same far from the origin, and thin or flat ones.
        box("box far", far, 3, 50, 0.5),
        cylinder("cylinder far", far, 0.25, 80),
        sphere("sphere far", far, 7),
        torus("torus far", far, 30, 0.5),
        cone("frustum far", far, 6, 2, 9),
    };
    // Relative errors: area, volume and inertia against themselves; the centroid against the
    // body's box diagonal (as tests/oracle measures it), the products against its largest moment.
    const char *labels[] = {"area", "volume", "centroid", "inertia", "products"};
    double worst[5] = {0, 0, 0, 0, 0};
    std::printf("%-30s %10s %10s %10s %10s %10s\n", "body", labels[0], labels[1], labels[2], labels[3], labels[4]);
    for (const Expected &e : cases)
    {
        SmBrep *p = nullptr;
        CHECK(e.make(p) == SM_SUCCESS && p != nullptr);
        smlib_brep b{p};
        double m[11] = {0};
        int32_t s = smlib_brep_mass_properties(&b, 1e-8, m);
        CHECK(s == SMLIB_OK);
        double area = 0;
        CHECK(smlib_brep_area(&b, 1e-8, &area) == SMLIB_OK);
        SmExtent3d box;
        CHECK(p->CalculateBoundingBox(box) == SM_SUCCESS);
        double size = box.GetSize().Length();
        double imax = std::fmax(e.v[5], std::fmax(e.v[6], e.v[7]));
        double err[5] = {std::fabs(m[0] / e.v[0] - 1), std::fabs(m[1] / e.v[1] - 1), 0, 0, 0};
        err[0] = std::fmax(err[0], std::fabs(area / e.v[0] - 1));
        for (int i = 2; i < 5; ++i)
        {
            err[2] = std::fmax(err[2], std::fabs(m[i] - e.v[i]) / size);
        }
        for (int i = 5; i < 8; ++i)
        {
            err[3] = std::fmax(err[3], std::fabs(m[i] / e.v[i] - 1));
        }
        for (int i = 8; i < 11; ++i)
        {
            err[4] = std::fmax(err[4], std::fabs(m[i] - e.v[i]) / imax);
        }
        std::printf("%-30s", e.name);
        for (int i = 0; i < 5; ++i)
        {
            std::printf(" %10.1e", err[i]);
            worst[i] = std::fmax(worst[i], err[i]);
        }
        std::printf("\n");
        delete p;
    }
    // Asked for 1e-8; tests/oracle accepts 1e-3.
    for (int i = 0; i < 5; ++i)
    {
        if (!(worst[i] < 1e-6))
        {
            std::fprintf(stderr, "%s: worst relative error %.2e\n", labels[i], worst[i]);
            failures++;
        }
    }

    // A sheet has no material: its area, no volume or inertia, and the centroid of its area
    // (a rectangle's: the middle of its box).
    SmBrep *sheet = nullptr;
    CHECK(SmApiCreatePlane(at0, 10.0, 20.0, sheet) == SM_SUCCESS && sheet != nullptr);
    smlib_brep sb{sheet};
    double m[11], a = 0, lo[3], hi[3];
    CHECK(smlib_brep_mass_properties(&sb, 1e-4, m) == SMLIB_OK && std::fabs(m[0] - 200.0) < 1e-9 && m[1] == 0.0);
    CHECK(smlib_brep_bounds(&sb, lo, hi) == SMLIB_OK);
    for (int i = 0; i < 3; ++i)
    {
        CHECK(std::fabs(m[2 + i] - 0.5 * (lo[i] + hi[i])) < 1e-9);
    }
    for (int i = 5; i < 11; ++i)
    {
        CHECK(m[i] == 0.0);
    }
    CHECK(smlib_brep_area(&sb, 1e-4, &a) == SMLIB_OK && std::fabs(a - 200.0) < 1e-9);
    delete sheet;

    // A sheet hung on a solid (SolidCraft's blend sheet stitched between two boxes): a 10 mm cube
    // and a 10 x 10 sheet on from its top front edge, merged into one body. Its area is every
    // face's, the cube's 600 and the sheet's 100, not the cube's boundary alone; its volume the
    // cube's.
    {
        SmBrep *cube = nullptr, *flap = nullptr, *both = nullptr;
        CHECK(SmApiCreateBox(SmVector3d(0, 0, 0), 10.0, 10.0, 10.0, cube) == SM_SUCCESS && cube != nullptr);
        CHECK(SmApiCreatePlane(SmVector3d(10, 0, 10), 10.0, 10.0, flap) == SM_SUCCESS && flap != nullptr);
        CHECK(cube != nullptr && flap != nullptr && SmApiBooleanMerge(cube, flap, both) == SM_SUCCESS && both != nullptr);
        if (both != nullptr)
        {
            smlib_brep merged{both};
            double hm[11];
            CHECK(smlib_brep_mass_properties(&merged, 1e-8, hm) == SMLIB_OK);
            std::printf("cube with a sheet hung on it: area %.9f volume %.9f\n", hm[0], hm[1]);
            CHECK(std::fabs(hm[0] - 700.0) < 1e-6 && std::fabs(hm[1] - 1000.0) < 1e-6);
            merged.p = nullptr;
            delete both;
        }
    }

    // The tight box of part of a NURBS arc is the part's: a quarter circle of radius 50 over the
    // first third of its parameter reaches y = 25 or so, not the arc's last knot at (0, 50). (An
    // imported sphere cap's seam, 30 degrees of a quarter-circle NURBS, boxed the cap down to the
    // equator's side.)
    {
        const SmContext &ctx = *SmApiGetOrCreateContext();
        SmAxis2Placement frame(SmPoint3d(0, 0, 0), SmVector3d(1, 0, 0), SmVector3d(0, 1, 0));
        const SmExtent1d quarter(0.0, 90.0);
        SmCircle *arc = nullptr;
        CHECK(SmCircle::CreateCanonical(ctx, frame, 50.0, arc, &quarter) == SM_SUCCESS && arc != nullptr);
        if (arc != nullptr)
        {
            CHECK(arc->EditParameterization(quarter) == SM_SUCCESS);
            // The arc's own parameter is not the angle: its box is from (50, 0) to its point at 30.
            SmExtent3d box;
            SmPoint3d end;
            CHECK(arc->CalculateTightBoundingBox(SmExtent1d(0.0, 30.0), &box) == SM_SUCCESS);
            CHECK(arc->EvaluatePoint(30.0, end) == SM_SUCCESS && end.y > 20.0 && end.y < 30.0);
            CHECK(std::fabs(box.GetUMin() - end.x) < 1e-9 && std::fabs(box.GetUMax() - 50.0) < 1e-9);
            CHECK(std::fabs(box.GetVMin()) < 1e-9 && std::fabs(box.GetVMax() - end.y) < 1e-9);
            delete arc;
        }
    }

    // The tight box of a trimmed face is the face's, not its surface's: a crowned patch whose
    // surface peaks outside its trim. A bicubic height field z(x, y) over [-20, 20]^2 (8 uniform
    // spans each way, control points at their Greville abscissae, so its parameters are x and y)
    // is zero but for two control points: the middle one, height 9, whose basis function is
    // nonzero on (-10, 10) each way and 2/3 at its middle, and the corner one, height 50, nonzero
    // on [-20, -15) each way only. Trimmed to the square [-10, 10]^2, z = 0 on the trim and the
    // face rises to 9 * (2/3)^2 = 4 at (0, 0); the surface's highest point, z = 50 at the corner
    // (-20, -20), is outside the trim. The box is (-10, -10, 0) to (10, 10, 4), not its edges'
    // (z up to 0) as when the surface's peak alone was classified.
    {
        const SmContext &ctx = *SmApiGetOrCreateContext();
        const ULONG n = 11;
        SmTArray<double> unique, greville;
        SmTArray<ULONG> mult;
        std::vector<double> knots;
        for (int k = 0; k <= 8; ++k)
        {
            unique.Add(-20.0 + 5.0 * k);
            mult.Add(k == 0 || k == 8 ? 4 : 1);
            for (ULONG m = 0; m < (k == 0 || k == 8 ? 4u : 1u); ++m)
            {
                knots.push_back(-20.0 + 5.0 * k);
            }
        }
        for (ULONG i = 0; i < n; ++i)
        {
            greville.Add((knots[i + 1] + knots[i + 2] + knots[i + 3]) / 3.0);
        }
        SmTArray<SmPoint3d> net;
        net.SetSize(n * n);
        for (ULONG i = 0; i < n; ++i)
        {
            for (ULONG j = 0; j < n; ++j)
            {
                const double z = (i == 5 && j == 5) ? 9.0 : (i == 0 && j == 0) ? 50.0 : 0.0;
                net[i * n + j] = SmPoint3d(greville[i], greville[j], z);
            }
        }
        SmBSplineSurface *field = nullptr;
        CHECK(SmBSplineSurface::CreateCanonical(ctx, 3, 3, net, SM_SF_UNSPECIFIED, mult, mult, unique, unique, SM_KT_UNSPECIFIED,
                                                nullptr, nullptr, field) == SM_SUCCESS && field != nullptr);
        const SmPoint3d corner[4] = {SmPoint3d(-10, -10, 0), SmPoint3d(10, -10, 0), SmPoint3d(10, 10, 0), SmPoint3d(-10, 10, 0)};
        SmTArray<SmCurve *> sides;
        for (int k = 0; k < 4; ++k)
        {
            SmLine *side = nullptr;
            CHECK(SmLine::CreateLineSegment(ctx, 3, corner[k], corner[(k + 1) % 4], side) == SM_SUCCESS && side != nullptr);
            sides.Add(side);
        }
        SmTArray<ULONG> loops;
        loops.Add(4);
        SmBrep *patch = nullptr;
        CHECK(field != nullptr && SmApiTrimSurfaceWith3dCurves(field, loops, sides, patch) == SM_SUCCESS && patch != nullptr);
        if (patch != nullptr)
        {
            smlib_brep pb{patch};
            double plo[3], phi[3];
            CHECK(smlib_brep_bounds(&pb, plo, phi) == SMLIB_OK);
            std::printf("crowned patch peaking outside its trim: box (%.9f %.9f %.9f) to (%.9f %.9f %.9f)\n", plo[0], plo[1], plo[2], phi[0],
                        phi[1], phi[2]);
            const double want_lo[3] = {-10, -10, 0}, want_hi[3] = {10, 10, 4};
            for (int i = 0; i < 3; ++i)
            {
                CHECK(std::fabs(plo[i] - want_lo[i]) < 1e-6 && std::fabs(phi[i] - want_hi[i]) < 1e-6);
            }
            pb.p = nullptr;
            delete patch;
        }
    }

    // The tight box of a threaded part is the part's: a thread's walls (rational screw surfaces of
    // many spans, smlib_screw_sweep) are cut off by the part's ends, so each wall's highest and
    // lowest points are on its edges there, the wall rising past them along the helix. External:
    // the ISO 68-1 groove of an M10x1.5 thread (as tests/sweeps.cpp cuts it) running from below
    // z = 0 to above 12, cut from a shaft of radius 5 and 12 long: (-5, -5, 0) to (5, 5, 12).
    // Internal: the groove of the M10x1.5 nut thread (P/8 wide at the major radius, widening by
    // 2 tan 30 per unit of radius inwards) cut into the bore of radius 4.3 of a 20 x 20 x 8 block:
    // the block's box.
    {
        auto polygon = [](const std::vector<double> &rz, std::vector<int64_t> &kinds, std::vector<double> &reals)
        {
            const size_t n = rz.size() / 2;
            for (size_t i = 0; i < n; ++i)
            {
                const size_t j = (i + 1) % n;
                kinds.push_back(SMLIB_CURVE_LINE);
                reals.insert(reals.end(), {rz[2 * i], 0.0, rz[2 * i + 1], rz[2 * j], 0.0, rz[2 * j + 1]});
            }
        };
        const double p = 1.5, h = p * std::sqrt(3.0) / 2.0, slope = std::tan(PI / 6.0), tol = p * 1e-4;
        const double O[3] = {0, 0, 0}, Z[3] = {0, 0, 1};
        struct Threaded
        {
            const char *name;
            bool external;
            std::vector<double> rz;
            double turns;
            double lo[3], hi[3];
        };
        const double r0 = 5.0 - 5.0 * h / 8.0, r1 = 5.0 + 0.05 * p, w0 = p / 8.0, w1 = (p / 4.0 + 2.0 * (r1 - r0) * slope) / 2.0;
        const double ri = 4.2, ro = 5.0, wi0 = p / 16.0, wi1 = wi0 + (ro - ri) * slope;
        const Threaded parts[] = {
            {"M10x1.5 shaft 12", true, {r0, -p - w0, r1, -p - w1, r1, -p + w1, r0, -p + w0}, 15.0 / p, {-5, -5, 0}, {5, 5, 12}},
            {"M10x1.5 nut 8", false, {ri, -p - wi1, ro, -p - wi0, ro, -p + wi0, ri, -p + wi1}, 11.0 / p, {-10, -10, 0}, {10, 10, 8}},
        };
        for (const Threaded &t : parts)
        {
            std::vector<int64_t> kinds;
            std::vector<double> reals;
            polygon(t.rz, kinds, reals);
            smlib_brep *groove = nullptr, *part = nullptr, *cut = nullptr;
            CHECK(smlib_screw_sweep(kinds.data(), static_cast<int64_t>(kinds.size()), reals.data(), static_cast<int64_t>(reals.size()), O, Z, p,
                                    t.turns, tol, &groove) == SMLIB_OK);
            if (t.external)
            {
                CHECK(smlib_make_cylinder(O, 5.0, 12.0, &part) == SMLIB_OK);
            }
            else
            {
                smlib_brep *block = nullptr, *bore = nullptr;
                const double corner[3] = {-10, -10, 0}, size[3] = {20, 20, 8}, below[3] = {0, 0, -1};
                CHECK(smlib_make_box(corner, size, &block) == SMLIB_OK);
                CHECK(smlib_make_cylinder(below, 4.3, 10.0, &bore) == SMLIB_OK);
                CHECK(block != nullptr && bore != nullptr && smlib_boolean(block, bore, SMLIB_BOOLEAN_SUBTRACT, &part) == SMLIB_OK);
                smlib_brep_free(block);
                smlib_brep_free(bore);
            }
            CHECK(groove != nullptr && part != nullptr && smlib_boolean(part, groove, SMLIB_BOOLEAN_SUBTRACT, &cut) == SMLIB_OK && cut != nullptr);
            if (cut != nullptr)
            {
                double lo[3], hi[3];
                const auto start = std::chrono::steady_clock::now();
                CHECK(smlib_brep_bounds(cut, lo, hi) == SMLIB_OK);
                const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
                std::printf("%s: box (%.9f %.9f %.9f) to (%.9f %.9f %.9f) in %.2f s\n", t.name, lo[0], lo[1], lo[2], hi[0], hi[1], hi[2], secs);
                for (int i = 0; i < 3; ++i)
                {
                    CHECK(std::fabs(lo[i] - t.lo[i]) < 1e-6 && std::fabs(hi[i] - t.hi[i]) < 1e-6);
                }
            }
            smlib_brep_free(groove);
            smlib_brep_free(part);
            smlib_brep_free(cut);
        }
    }

    // Point classification: exact, at the body's own zone tolerance. A 10 x 20 x 30 box at the
    // origin: inside, outside, on a face, on an edge, on a vertex. A cylinder of radius 5 about +Z:
    // off its curved wall by half its tolerance is on it, by ten tolerances inside or outside.
    {
        auto where = [](const smlib_brep *b, double x, double y, double z)
        {
            int32_t c = -1;
            const int32_t s = smlib_brep_classify_point(b, x, y, z, &c);
            return s == SMLIB_OK ? c : -100 - s;
        };
        smlib_brep *bx = nullptr, *cyl = nullptr;
        const double o[3] = {0, 0, 0}, size[3] = {10, 20, 30};
        CHECK(smlib_make_box(o, size, &bx) == SMLIB_OK && bx != nullptr);
        CHECK(smlib_make_cylinder(o, 5.0, 10.0, &cyl) == SMLIB_OK && cyl != nullptr);
        if (bx != nullptr && cyl != nullptr)
        {
            CHECK(where(bx, 5, 10, 15) == SMLIB_POINT_INSIDE);
            CHECK(where(bx, 5, 10, 31) == SMLIB_POINT_OUTSIDE);
            CHECK(where(bx, -1, -1, -1) == SMLIB_POINT_OUTSIDE);
            CHECK(where(bx, 5, 10, 30) == SMLIB_POINT_ON_BOUNDARY);
            CHECK(where(bx, 10, 10, 30) == SMLIB_POINT_ON_BOUNDARY);
            CHECK(where(bx, 10, 20, 30) == SMLIB_POINT_ON_BOUNDARY);
            CHECK(where(bx, 0, 0, 0) == SMLIB_POINT_ON_BOUNDARY);

            const double tol = static_cast<double>(SmTol::GetZoneTol3d(cyl->p));
            std::printf("cylinder zone tolerance %.3g mm\n", tol);
            CHECK(tol > 0.0 && tol < 1e-2);
            const double c45 = std::sqrt(0.5);
            for (double r : {5.0 - 0.5 * tol, 5.0 + 0.5 * tol})
            {
                CHECK(where(cyl, r * c45, r * c45, 5) == SMLIB_POINT_ON_BOUNDARY);
            }
            CHECK(where(cyl, (5.0 - 10 * tol) * c45, (5.0 - 10 * tol) * c45, 5) == SMLIB_POINT_INSIDE);
            CHECK(where(cyl, (5.0 + 10 * tol) * c45, (5.0 + 10 * tol) * c45, 5) == SMLIB_POINT_OUTSIDE);
            CHECK(where(cyl, 0, 0, 5) == SMLIB_POINT_INSIDE);
            CHECK(where(cyl, 0, 0, 10) == SMLIB_POINT_ON_BOUNDARY);

            int32_t c = 7;
            CHECK(smlib_brep_classify_point(nullptr, 0, 0, 0, &c) == SMLIB_ERR_ARG);
            CHECK(smlib_brep_classify_point(bx, NAN, 0, 0, &c) == SMLIB_ERR_ARG);
            CHECK(smlib_brep_classify_point(bx, 0, INFINITY, 0, &c) == SMLIB_ERR_ARG);
            CHECK(smlib_brep_classify_point(bx, 0, 0, 0, nullptr) == SMLIB_ERR_ARG);
            CHECK(c == 7);
        }
        smlib_brep_free(bx);
        smlib_brep_free(cyl);
    }

    // Rejected arguments.
    CHECK(smlib_brep_mass_properties(nullptr, 1e-4, m) == SMLIB_ERR_ARG);
    CHECK(smlib_brep_area(nullptr, 1e-4, &a) == SMLIB_ERR_ARG);
    smlib_brep *bx = nullptr;
    const double o3[3] = {0, 0, 0}, s3[3] = {1, 2, 3};
    CHECK(smlib_make_box(o3, s3, &bx) == SMLIB_OK);
    CHECK(smlib_brep_mass_properties(bx, NAN, m) == SMLIB_ERR_ARG);
    CHECK(smlib_brep_mass_properties(bx, 1e-4, nullptr) == SMLIB_ERR_ARG);
    CHECK(smlib_brep_area(bx, INFINITY, &a) == SMLIB_ERR_ARG);
    CHECK(smlib_brep_area(bx, 1e-4, nullptr) == SMLIB_ERR_ARG);
    smlib_brep_free(bx);

    if (failures == 0)
    {
        std::printf("smlib_brep_mass_properties test: ok\n");
    }
    return failures == 0 ? 0 : 1;
}
