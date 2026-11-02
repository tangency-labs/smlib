// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

// smlib_stitch: a sheet sewn onto a solid along edges two of the solid's faces already share.

#include "../src/smlib_shim.h"

#include <SmApiBrep.h>
#include <SmApiPrimitives.h>
#include <SmApiTrimmedSurfaces.h>
#include <SmBSplineSurface.h>
#include <SmLine.h>

#include <cmath>
#include <cstdio>
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

// A 10 mm cube and a cover over its top in one body, not yet sewn: the cover a height field
// z(x, y) over [0, 10]^2 with `degree` (2 or 3) and uniform spans, its control points at their
// Greville abscissae (so its parameters are x and y), z = 10 but for the middle one, 10 + `rise`,
// trimmed to the square, so its rim is the cube's top edges. NULL when it cannot be made.
static SmBrep *cube_and_cover(ULONG degree, ULONG spans, double rise)
{
    const SmContext &ctx = *SmApiGetOrCreateContext();
    const ULONG n = spans + degree;
    const double step = 10.0 / spans;
    SmTArray<double> unique, greville;
    SmTArray<ULONG> mult;
    std::vector<double> knots;
    for (ULONG k = 0; k <= spans; ++k)
    {
        const ULONG m = k == 0 || k == spans ? degree + 1 : 1;
        unique.Add(step * k);
        mult.Add(m);
        knots.insert(knots.end(), m, step * k);
    }
    for (ULONG i = 0; i < n; ++i)
    {
        double sum = 0.0;
        for (ULONG d = 1; d <= degree; ++d)
        {
            sum += knots[i + d];
        }
        greville.Add(sum / degree);
    }
    SmTArray<SmPoint3d> net;
    net.SetSize(n * n);
    for (ULONG i = 0; i < n; ++i)
    {
        for (ULONG j = 0; j < n; ++j)
        {
            net[i * n + j] = SmPoint3d(greville[i], greville[j], (2 * i + 1 == n && 2 * j + 1 == n) ? 10.0 + rise : 10.0);
        }
    }
    SmBSplineSurface *field = nullptr;
    CHECK(SmBSplineSurface::CreateCanonical(ctx, degree, degree, net, SM_SF_UNSPECIFIED, mult, mult, unique, unique, SM_KT_UNSPECIFIED,
                                            nullptr, nullptr, field) == SM_SUCCESS &&
          field != nullptr);
    const SmPoint3d corner[4] = {SmPoint3d(0, 0, 10), SmPoint3d(10, 0, 10), SmPoint3d(10, 10, 10), SmPoint3d(0, 10, 10)};
    SmTArray<SmCurve *> sides;
    for (int k = 0; k < 4; ++k)
    {
        SmLine *side = nullptr;
        CHECK(SmLine::CreateLineSegment(ctx, 3, corner[k], corner[(k + 1) % 4], side) == SM_SUCCESS && side != nullptr);
        sides.Add(side);
    }
    SmTArray<ULONG> loops;
    loops.Add(4);
    SmBrep *cover = nullptr, *cube = nullptr;
    CHECK(field != nullptr && SmApiTrimSurfaceWith3dCurves(field, loops, sides, cover) == SM_SUCCESS && cover != nullptr);
    CHECK(SmApiCreateBox(SmVector3d(0, 0, 0), 10.0, 10.0, 10.0, cube) == SM_SUCCESS && cube != nullptr);
    if (cover == nullptr || cube == nullptr)
    {
        delete cover;
        delete cube;
        return nullptr;
    }
    CHECK(cube->MergeBrep(*cover) == SM_SUCCESS);
    delete cover;
    return cube;
}

// A biquadratic Bezier cover rising 9 at its middle control point, z = 10 + 9 * 2s(1-s) * 2t(1-t)
// (s = x/10, t = y/10), leaves the cube's top edges at an angle (but at the corners); the space
// between it and the cube's top encloses 9 * (10/3)^2 = 100 mm^3. Sewn with the cube, that space
// is part of the solid and the cube's top, now inside it, goes: a closed manifold solid of the
// cube's five other faces and the cover, 1100 mm^3.
static void cover_sewn_onto_a_cube()
{
    SmBrep *both = cube_and_cover(2, 1, 9.0);
    if (both == nullptr)
    {
        return;
    }
    smlib_brep b{both};
    smlib_brep *sewn = nullptr;
    int32_t closed = -1;
    CHECK(smlib_stitch(&b, 1e-6, &closed, &sewn) == SMLIB_OK && sewn != nullptr);
    if (sewn != nullptr)
    {
        int64_t counts[3] = {0, 0, 0};
        double m[11] = {0};
        CHECK(smlib_brep_counts(sewn, counts) == SMLIB_OK);
        CHECK(smlib_brep_mass_properties(sewn, 1e-8, m) == SMLIB_OK);
        std::printf("cover sewn onto a cube: closed %d, %lld faces, volume %.9f\n", closed, static_cast<long long>(counts[0]), m[1]);
        CHECK(closed == 1);
        CHECK(counts[0] == 6);
        CHECK(std::fabs(m[1] - 1100.0) < 1e-6);
        smlib_brep_free(sewn);
    }
    b.p = nullptr;
    delete both;
}

// A bicubic cover over 8 spans whose only raised control point is the middle one lies on the
// cube's top for 2.5 mm in from every edge (the basis functions of the points that are 10 are
// all that is nonzero there): it overlaps the top, and there is no sector between the faces at
// the rim to sew it into. Refused, not sewn so that the cube's inside joins the outside.
static void cover_lying_on_the_cube_is_refused()
{
    SmBrep *both = cube_and_cover(3, 8, 9.0);
    if (both == nullptr)
    {
        return;
    }
    smlib_brep b{both};
    smlib_brep *sewn = nullptr;
    int32_t closed = -1;
    const int32_t s = smlib_stitch(&b, 1e-6, &closed, &sewn);
    std::printf("cover lying on the cube's top: status %d\n", s);
    CHECK(s == SMLIB_ERR_UNSUPPORTED && sewn == nullptr);
    smlib_brep_free(sewn);
    b.p = nullptr;
    delete both;
}

// A 10 x 10 sheet on from the cube's top edge at x = 10 (the plane z = 10 on past it) closes
// nothing with the cube: it stays the sheet beside it that the manifold stitch leaves, the cube's
// 1000 mm^3 and six faces as they were, not glued into an edge three faces would share.
static void flap_on_a_cube_edge_stays_beside_it()
{
    SmBrep *cube = nullptr, *flap = nullptr;
    CHECK(SmApiCreateBox(SmVector3d(0, 0, 0), 10.0, 10.0, 10.0, cube) == SM_SUCCESS && cube != nullptr);
    CHECK(SmApiCreatePlane(SmVector3d(10, 0, 10), 10.0, 10.0, flap) == SM_SUCCESS && flap != nullptr);
    if (cube == nullptr || flap == nullptr)
    {
        delete cube;
        delete flap;
        return;
    }
    CHECK(cube->MergeBrep(*flap) == SM_SUCCESS);
    delete flap;
    smlib_brep b{cube};
    smlib_brep *sewn = nullptr;
    int32_t closed = -1;
    CHECK(smlib_stitch(&b, 1e-6, &closed, &sewn) == SMLIB_OK && sewn != nullptr);
    if (sewn != nullptr)
    {
        int64_t counts[3] = {0, 0, 0};
        double m[11] = {0};
        CHECK(smlib_brep_counts(sewn, counts) == SMLIB_OK);
        CHECK(smlib_brep_mass_properties(sewn, 1e-8, m) == SMLIB_OK);
        std::printf("flap on a cube edge: closed %d, %lld faces, volume %.9f, area %.9f\n", closed, static_cast<long long>(counts[0]), m[1], m[0]);
        CHECK(closed == 0 && counts[0] == 7);
        CHECK(std::fabs(m[1] - 1000.0) < 1e-6 && std::fabs(m[0] - 700.0) < 1e-6);
        smlib_brep_free(sewn);
    }
    b.p = nullptr;
    delete cube;
}

int main()
{
    cover_sewn_onto_a_cube();
    cover_lying_on_the_cube_is_refused();
    flap_on_a_cube_edge_stays_beside_it();
    if (failures == 0)
    {
        std::printf("smlib_stitch test: ok\n");
    }
    return failures == 0 ? 0 : 1;
}
